"""Stable pipeline entrypoint for epistasis Wald/BH testing.

Purpose:
    Run the reusable M2 and M3 Wald/BH procedures against one exported task0
    result directory and write generic output tables for downstream figure and
    manuscript scripts.

Inputs:
    - task0-style exported CSV directories under output/<subdir>
    - mutation-rate method name
    - BH alpha thresholds

Outputs:
    - generic CSV summaries under output/<subdir>/epistasis_testing

Assumptions:
    - Statistical calculations remain implemented in Python.
    - analysis_for_paper_2.R consumes these tables; it does not invoke this
      script automatically.
"""

from __future__ import annotations

import argparse
import os
import platform
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
TASK7_DIR = PROJECT_ROOT / "code" / "result_analysis" / "multiple_testing_helpers"
if str(TASK7_DIR) not in sys.path:
    sys.path.insert(0, str(TASK7_DIR))

from epistasis_testing_helpers import DEFAULT_RUN_ROOT, resolve_result_dir  # noqa: E402
from run_m3_expected_vs_observed_wt_one_sided import (  # noqa: E402
    run_expected_vs_observed_analysis,
)
from run_m3_regime_multiple_testing import run_m3_regime_analysis  # noqa: E402
from run_pairwise_multiple_testing import run_model_order  # noqa: E402


PRIMARY_ALPHA = 0.05
RELAXED_ALPHA = 0.10


def default_run_root() -> Path:
    """Purpose: choose the default exported run root. Inputs: optional environment variable. Outputs: output/<subdir> path. Assumptions: task0_indel_exclusion remains the fallback exported result set."""

    output_subdir = os.getenv("LUNG_SMOKING_OUTPUT_SUBDIR", "").strip()
    if output_subdir:
        return PROJECT_ROOT / "output" / output_subdir
    return DEFAULT_RUN_ROOT


def parse_args() -> argparse.Namespace:
    """Purpose: parse CLI inputs. Inputs: command-line arguments. Outputs: normalized pipeline options. Assumptions: run-root points to an exported output directory or method-specific result directory."""

    parser = argparse.ArgumentParser(
        description="Run the stable M2/M3 epistasis Wald/BH testing pipeline."
    )
    parser.add_argument("--run-root", default=str(default_run_root()))
    parser.add_argument("--method", default="variant", choices=["variant", "cesR"])
    parser.add_argument("--alpha", type=float, default=PRIMARY_ALPHA)
    parser.add_argument("--relaxed-alpha", type=float, default=RELAXED_ALPHA)
    return parser.parse_args()


def copy_if_exists(source: Path, target: Path) -> None:
    """Purpose: copy a generated file to a stable name if present. Inputs: source and target paths. Outputs: copied file on disk. Assumptions: source was created earlier in the same pipeline run."""

    if source.exists():
        shutil.copy2(source, target)


def write_manifest(
    output_root: Path,
    run_root: Path,
    result_dir: Path,
    method: str,
    alpha: float,
    relaxed_alpha: float,
    m2_summary: dict[str, object],
    m3_results: dict[str, pd.DataFrame],
    evo_results: dict[str, pd.DataFrame],
) -> None:
    """Purpose: record stable pipeline provenance. Inputs: paths, thresholds, and in-memory summaries. Outputs: markdown manifest. Assumptions: concise provenance supports reviewer audit and reproducibility."""

    lines = [
        "# Epistasis Testing Pipeline Manifest",
        "",
        "## Inputs",
        f"- Source run root: `{run_root}`",
        f"- Source result dir: `{result_dir}`",
        f"- Method: `{method}`",
        f"- Primary alpha: `{alpha}`",
        f"- Relaxed alpha: `{relaxed_alpha}`",
        "",
        "## Generated Families",
        "- M2 pairwise WT-reference Wald/BH contrasts.",
        "- M3 regime-specific WT-reference and descriptive double-vs-single families.",
        "- M3 ES-LUAD expected-vs-observed one-sided WT family.",
        "",
        "## Software",
        f"- Python executable: `{sys.executable}`",
        f"- Python version: `{platform.python_version()}`",
        f"- pandas version: `{pd.__version__}`",
        f"- numpy version: `{np.__version__}`",
        "",
        "## Row Counts",
        f"- M2 testable contrasts: `{m2_summary.get('n_testable', 'NA')}`",
        f"- M2 BH-significant contrasts: `{m2_summary.get('n_bh_signif', 'NA')}`",
        f"- M3 regime triad-target tests: `{len(m3_results['primary_tests'])}`",
        f"- M3 descriptive double-vs-single rows: `{len(m3_results['descriptive_tests'])}`",
        f"- M3 expected-vs-observed candidate rows: `{len(evo_results['candidates'])}`",
        "",
        "## Stable Outputs",
        "- `M2_pairwise_primary_tests.csv`",
        "- `M2_pairwise_summary.csv`",
        "- `M3_regime_primary_wt_tests.csv`",
        "- `M3_regime_descriptive_single_reference_tests.csv`",
        "- `M3_regime_summary.csv`",
        "- `M3_expected_vs_observed_wt_one_sided_tests.csv`",
    ]
    (output_root / "epistasis_testing_manifest.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    """Purpose: execute the stable epistasis testing pipeline end to end. Inputs: CLI options. Outputs: stable generic CSVs under output/<subdir>/epistasis_testing. Assumptions: task0-derived M2/M3 CSV exports already exist."""

    args = parse_args()
    output_root = run_epistasis_testing(
        run_root=args.run_root,
        method=args.method,
        alpha=args.alpha,
        relaxed_alpha=args.relaxed_alpha,
    )
    print(f"Wrote stable epistasis testing outputs to {output_root}")


def run_epistasis_testing(
    run_root: str | Path | None = None,
    method: str = "variant",
    alpha: float = PRIMARY_ALPHA,
    relaxed_alpha: float = RELAXED_ALPHA,
) -> Path:
    """Purpose: execute the stable epistasis testing pipeline from Python code. Inputs: run reference, method, and BH thresholds. Outputs: output directory path. Assumptions: task0-derived CSV exports already exist for the requested run."""

    run_ref = str(run_root) if run_root is not None else str(default_run_root())
    run_root, result_dir = resolve_result_dir(run_ref, method)
    output_root = run_root / "epistasis_testing"
    plots_dir = output_root / "plots"
    output_root.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    m2_summary = run_model_order(
        result_dir=result_dir,
        comparison_dir=run_root / "comparison",
        output_root=output_root,
        plots_dir=plots_dir,
        method=method,
        model_order=2,
        alpha=alpha,
    )
    if m2_summary.get("status") != "completed":
        raise RuntimeError(f"M2 pairwise analysis did not complete: {m2_summary}")
    copy_if_exists(output_root / "M2_primary_contrasts.csv", output_root / "M2_pairwise_primary_tests.csv")
    copy_if_exists(output_root / "M2_summary_joint.csv", output_root / "M2_pairwise_summary.csv")

    m3_results = run_m3_regime_analysis(
        result_dir=result_dir,
        output_root=output_root,
        run_root=run_root,
        method=method,
        alpha=alpha,
    )

    evo_results = run_expected_vs_observed_analysis(
        result_dir=result_dir,
        output_root=output_root,
        run_root=run_root,
        method=method,
        alpha=alpha,
        relaxed_alpha=relaxed_alpha,
    )

    write_manifest(
        output_root=output_root,
        run_root=run_root,
        result_dir=result_dir,
        method=method,
        alpha=alpha,
        relaxed_alpha=relaxed_alpha,
        m2_summary=m2_summary,
        m3_results=m3_results,
        evo_results=evo_results,
    )
    return output_root


if __name__ == "__main__":
    main()
