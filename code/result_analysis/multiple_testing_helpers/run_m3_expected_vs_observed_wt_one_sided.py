"""Alternative ES-LUAD M=3 WT-reference tests for expected-vs-observed figure.

This script is downstream of task0 model CSV exports. It constructs the
figure-specific testing family requested for the reviewer response: ES-LUAD
triad-targets with either two synergistic component pairwise effects (S+S) or
no significant component pairwise effects (N+N), tested only for
double-mutant selection greater than WT selection.
"""

from __future__ import annotations

import argparse
import math
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from run_m3_regime_multiple_testing import (
    MANUSCRIPT_INTERACTIONS,
    MODEL_ORDER,
    build_m3_regime_rows,
)
from epistasis_testing_helpers import (
    DEFAULT_RUN_ROOT,
    PROJECT_ROOT,
    add_wald_columns,
    bh_adjust,
    load_transitions,
    resolve_result_dir,
)


DEFAULT_OUTPUT_ROOT = (
    PROJECT_ROOT
    / "output"
    / "task7_multiple_testing"
    / "m3_expected_vs_observed_wt_one_sided"
)
ES_LUAD_KEY = "smoking_plus"
PRIMARY_ALPHA = 0.05
RELAXED_ALPHA = 0.10
MIN_PLOT_RATIO = 1e-1
INCLUDED_PAIRWISE_CONTEXTS = ("S+S", "N+N")


def parse_args() -> argparse.Namespace:
    """Purpose: parse CLI inputs. Inputs: command-line args. Outputs: run options. Assumptions: standard project layout and task0-compatible exports."""

    parser = argparse.ArgumentParser(
        description=(
            "Run ES-LUAD M=3 expected-vs-observed one-sided WT-reference "
            "BH correction for S+S and N+N candidate interactions."
        )
    )
    parser.add_argument("--run-root", default=str(DEFAULT_RUN_ROOT))
    parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT))
    parser.add_argument("--method", default="variant", choices=["variant", "cesR"])
    parser.add_argument("--alpha", type=float, default=PRIMARY_ALPHA)
    parser.add_argument("--relaxed-alpha", type=float, default=RELAXED_ALPHA)
    return parser.parse_args()


def one_sided_upper_p_value(z_value: float) -> float:
    """Purpose: compute a one-sided upper-tail normal p-value. Inputs: Wald z statistic. Outputs: P(Z >= z). Assumptions: z follows an approximately standard normal null distribution."""

    if pd.isna(z_value):
        return np.nan
    return 0.5 * math.erfc(float(z_value) / math.sqrt(2.0))


def add_figure_annotations(candidates: pd.DataFrame) -> pd.DataFrame:
    """Purpose: add expected-vs-observed figure quantities. Inputs: candidate WT-reference tests. Outputs: candidates with descriptive plotting annotations. Assumptions: expected effect is the product of the two single-mutant-vs-WT ratios and plotting ratios are floored at 0.1 as in the manuscript validation script."""

    candidates = candidates.copy()
    candidates["component_pairwise_effect_1"] = (
        candidates["single_context_1"].astype(str) + "->" + candidates["mutated_gene"].astype(str)
    )
    candidates["component_pairwise_effect_2"] = (
        candidates["single_context_2"].astype(str) + "->" + candidates["mutated_gene"].astype(str)
    )
    candidates["component_pairwise_class_1"] = candidates["single_class_1"]
    candidates["component_pairwise_class_2"] = candidates["single_class_2"]
    candidates["component_pairwise_ratio_1"] = candidates["single_ratio_1"]
    candidates["component_pairwise_ratio_2"] = candidates["single_ratio_2"]
    candidates["candidate_pairwise_context"] = candidates["regime"].astype(str)
    candidates["n_component_pairwise_ci_nonoverlap"] = (
        candidates[["single_legacy_signif_1", "single_legacy_signif_2"]]
        .fillna(False)
        .astype(bool)
        .sum(axis=1)
    )
    candidates["expected_pairwise_ratio"] = (
        candidates["single_ratio_1"] * candidates["single_ratio_2"]
    )
    candidates["expected_for_plot"] = candidates["expected_pairwise_ratio"].clip(lower=MIN_PLOT_RATIO)
    candidates["observed_for_plot"] = candidates["ratio"].clip(lower=MIN_PLOT_RATIO)
    candidates["obs_over_exp_ratio"] = (
        candidates["observed_for_plot"] / candidates["expected_for_plot"]
    )
    candidates["higher_order_ci_nonoverlap_any_direction"] = candidates["legacy_signif"].fillna(False)
    candidates["higher_order_ci_above_wt"] = (
        candidates["epi_gamma_ci_low"] > candidates["wt_gamma_ci_high"]
    )
    candidates["higher_order_ci_below_wt"] = (
        candidates["wt_gamma_ci_low"] > candidates["epi_gamma_ci_high"]
    )
    candidates["which_signif_figure_logic"] = np.where(
        (candidates["candidate_pairwise_context"] == "S+S")
        & candidates["higher_order_ci_nonoverlap_any_direction"],
        "All",
        np.where(
            candidates["candidate_pairwise_context"] == "S+S",
            "Both pairwise only",
            np.where(
                candidates["higher_order_ci_nonoverlap_any_direction"],
                "Higher-order only",
                "None",
            ),
        ),
    )
    candidates["in_expected_vs_observed_selected_logic"] = (
        (candidates["candidate_pairwise_context"] == "S+S")
        | (
            (candidates["candidate_pairwise_context"] == "N+N")
            & candidates["higher_order_ci_nonoverlap_any_direction"]
        )
    )
    candidates["visible_in_current_expected_vs_observed_axes_approx"] = (
        candidates["in_expected_vs_observed_selected_logic"]
        & (candidates["expected_for_plot"] >= 1.0)
        & (candidates["expected_for_plot"] <= 1300.0)
        & (candidates["observed_for_plot"] >= 1.0)
    )
    return candidates


def build_candidate_tests(transitions: pd.DataFrame) -> pd.DataFrame:
    """Purpose: construct ES-LUAD S+S/N+N double-mutant-vs-WT tests. Inputs: M=3 transition table. Outputs: one row per eligible triad-target. Assumptions: candidate classification uses the same CI non-overlap single-effect rules as the M=3 regime analysis."""

    primary_tests, _, _ = build_m3_regime_rows(transitions)
    if primary_tests.empty:
        return primary_tests

    candidates = primary_tests.loc[
        (primary_tests["key"] == ES_LUAD_KEY)
        & (primary_tests["regime"].isin(INCLUDED_PAIRWISE_CONTEXTS))
    ].copy()
    candidates = add_figure_annotations(candidates)
    candidates["in_combined_expected_vs_observed_one_sided_family"] = True
    return candidates.reset_index(drop=True)


def add_one_sided_bh_columns(
    candidates: pd.DataFrame,
    alpha: float,
    relaxed_alpha: float,
) -> pd.DataFrame:
    """Purpose: add one-sided Wald p-values and combined BH calls. Inputs: candidate tests and alpha thresholds. Outputs: candidate table with p/q/significance columns. Assumptions: BH correction pools all testable ES-LUAD S+S and N+N candidates in one family."""

    if candidates.empty:
        return candidates

    candidates = add_wald_columns(candidates, "null_facing", "one_sided")
    candidates["one_sided_p_raw"] = candidates["one_sided_z_wald"].map(one_sided_upper_p_value)
    candidates["one_sided_raw_signif_alpha_0_05"] = (
        candidates["one_sided_testable"].fillna(False)
        & (candidates["one_sided_p_raw"] <= alpha)
    )
    candidates["one_sided_q_bh_combined"] = np.nan
    testable_idx = candidates.index[candidates["one_sided_testable"].fillna(False)]
    candidates.loc[testable_idx, "one_sided_q_bh_combined"] = bh_adjust(
        candidates.loc[testable_idx, "one_sided_p_raw"]
    )
    candidates["one_sided_bh_signif_alpha_0_05"] = (
        candidates["one_sided_testable"].fillna(False)
        & (candidates["one_sided_q_bh_combined"] <= alpha)
    )
    candidates["one_sided_bh_signif_alpha_0_1"] = (
        candidates["one_sided_testable"].fillna(False)
        & (candidates["one_sided_q_bh_combined"] <= relaxed_alpha)
    )
    candidates["one_sided_q_bh_within_pairwise_context"] = np.nan
    for context, idx in candidates.groupby("candidate_pairwise_context", sort=True).groups.items():
        context_testable_idx = candidates.loc[idx].index[
            candidates.loc[idx, "one_sided_testable"].fillna(False)
        ]
        candidates.loc[context_testable_idx, "one_sided_q_bh_within_pairwise_context"] = bh_adjust(
            candidates.loc[context_testable_idx, "one_sided_p_raw"]
        )
    candidates["one_sided_bh_within_pairwise_context_alpha_0_05"] = (
        candidates["one_sided_testable"].fillna(False)
        & (candidates["one_sided_q_bh_within_pairwise_context"] <= alpha)
    )
    candidates["one_sided_bh_within_pairwise_context_alpha_0_1"] = (
        candidates["one_sided_testable"].fillna(False)
        & (candidates["one_sided_q_bh_within_pairwise_context"] <= relaxed_alpha)
    )
    candidates["one_sided_bh_and_ci_above_wt_alpha_0_05"] = (
        candidates["one_sided_bh_signif_alpha_0_05"]
        & candidates["higher_order_ci_above_wt"].fillna(False)
    )
    candidates["one_sided_bh_and_ci_above_wt_alpha_0_1"] = (
        candidates["one_sided_bh_signif_alpha_0_1"]
        & candidates["higher_order_ci_above_wt"].fillna(False)
    )
    candidates["one_sided_bh_within_context_and_ci_above_wt_alpha_0_05"] = (
        candidates["one_sided_bh_within_pairwise_context_alpha_0_05"]
        & candidates["higher_order_ci_above_wt"].fillna(False)
    )
    candidates["one_sided_bh_within_context_and_ci_above_wt_alpha_0_1"] = (
        candidates["one_sided_bh_within_pairwise_context_alpha_0_1"]
        & candidates["higher_order_ci_above_wt"].fillna(False)
    )
    candidates["legacy_vs_one_sided_bh_alpha_0_05"] = np.where(
        ~candidates["one_sided_testable"].fillna(False),
        "untestable",
        np.where(
            candidates["higher_order_ci_above_wt"] & candidates["one_sided_bh_signif_alpha_0_05"],
            "ci_above_wt_and_bh",
            np.where(
                candidates["higher_order_ci_above_wt"]
                & ~candidates["one_sided_bh_signif_alpha_0_05"],
                "ci_above_wt_only",
                np.where(
                    ~candidates["higher_order_ci_above_wt"]
                    & candidates["one_sided_bh_signif_alpha_0_05"],
                    "bh_only",
                    "neither",
                ),
            ),
        ),
    )
    candidates["manuscript_reference"] = candidates.apply(
        lambda row: MANUSCRIPT_INTERACTIONS.get(
            (row["key"], row["gene_set"], row["mutated_gene"]),
            row.get("manuscript_reference", ""),
        ),
        axis=1,
    )
    candidates["manuscript_subset"] = candidates["manuscript_reference"].fillna("").ne("")
    return candidates.sort_values(
        [
            "one_sided_bh_signif_alpha_0_05",
            "one_sided_q_bh_combined",
            "candidate_pairwise_context",
            "gene_set",
            "mutated_gene",
        ],
        ascending=[False, True, True, True, True],
        na_position="last",
    ).reset_index(drop=True)


def summarize_tests(
    candidates: pd.DataFrame,
    alpha: float,
    relaxed_alpha: float,
) -> pd.DataFrame:
    """Purpose: summarize the alternative expected-vs-observed test family. Inputs: annotated candidate tests and alpha thresholds. Outputs: overall and context-level count table. Assumptions: all testable candidates are corrected together even when summarized by context."""

    rows = []
    groups: list[tuple[str, pd.DataFrame]] = [("combined_S+S_and_N+N", candidates)]
    groups.extend(
        (str(context), df)
        for context, df in candidates.groupby("candidate_pairwise_context", sort=True)
    )
    for label, df in groups:
        rows.append(
            {
                "summary_group": label,
                "n_total_tests": int(len(df)),
                "n_testable": int(df["one_sided_testable"].fillna(False).sum()),
                "n_untestable": int((~df["one_sided_testable"].fillna(False)).sum()),
                "n_ci_nonoverlap_any_direction": int(
                    df["higher_order_ci_nonoverlap_any_direction"].fillna(False).sum()
                ),
                "n_ci_above_wt": int(df["higher_order_ci_above_wt"].fillna(False).sum()),
                "n_ci_below_wt": int(df["higher_order_ci_below_wt"].fillna(False).sum()),
                "n_raw_one_sided_p_lt_alpha_0_05": int(
                    (df["one_sided_p_raw"] <= alpha).fillna(False).sum()
                ),
                "n_bh_signif_alpha_0_05": int(
                    df["one_sided_bh_signif_alpha_0_05"].fillna(False).sum()
                ),
                "n_bh_signif_alpha_0_1": int(
                    df["one_sided_bh_signif_alpha_0_1"].fillna(False).sum()
                ),
                "n_bh_within_context_signif_alpha_0_05": int(
                    df["one_sided_bh_within_pairwise_context_alpha_0_05"]
                    .fillna(False)
                    .sum()
                ),
                "n_bh_within_context_signif_alpha_0_1": int(
                    df["one_sided_bh_within_pairwise_context_alpha_0_1"]
                    .fillna(False)
                    .sum()
                ),
                "n_bh_and_ci_above_wt_alpha_0_05": int(
                    df["one_sided_bh_and_ci_above_wt_alpha_0_05"].fillna(False).sum()
                ),
                "n_bh_and_ci_above_wt_alpha_0_1": int(
                    df["one_sided_bh_and_ci_above_wt_alpha_0_1"].fillna(False).sum()
                ),
                "n_bh_within_context_and_ci_above_wt_alpha_0_05": int(
                    df["one_sided_bh_within_context_and_ci_above_wt_alpha_0_05"]
                    .fillna(False)
                    .sum()
                ),
                "n_bh_within_context_and_ci_above_wt_alpha_0_1": int(
                    df["one_sided_bh_within_context_and_ci_above_wt_alpha_0_1"]
                    .fillna(False)
                    .sum()
                ),
                "n_in_expected_vs_observed_selected_logic": int(
                    df["in_expected_vs_observed_selected_logic"].fillna(False).sum()
                ),
                "n_visible_in_current_expected_vs_observed_axes_approx": int(
                    df["visible_in_current_expected_vs_observed_axes_approx"]
                    .fillna(False)
                    .sum()
                ),
                "n_manuscript_subset": int(df["manuscript_subset"].fillna(False).sum()),
                "alpha": alpha,
                "relaxed_alpha": relaxed_alpha,
            }
        )
    return pd.DataFrame(rows)


def select_significant_interactions(candidates: pd.DataFrame, mode: str) -> pd.DataFrame:
    """Purpose: create a compact significant-interaction table. Inputs: annotated candidate tests and BH mode. Outputs: BH-significant rows at alpha 0.1 with alpha 0.05 flag retained. Assumptions: alpha 0.1 is the broad reviewer-facing sensitivity threshold."""

    if mode == "combined":
        alpha_0_1_col = "one_sided_bh_signif_alpha_0_1"
    elif mode == "within_pairwise_context":
        alpha_0_1_col = "one_sided_bh_within_pairwise_context_alpha_0_1"
    else:
        raise ValueError(f"Unknown significant-interaction mode: {mode}")

    columns = [
        "key",
        "gene_set",
        "mutated_gene",
        "candidate_pairwise_context",
        "component_pairwise_effect_1",
        "component_pairwise_class_1",
        "component_pairwise_ratio_1",
        "component_pairwise_effect_2",
        "component_pairwise_class_2",
        "component_pairwise_ratio_2",
        "epistatic_gt",
        "ratio",
        "delta_hat",
        "one_sided_se_delta",
        "one_sided_z_wald",
        "one_sided_p_raw",
        "one_sided_q_bh_combined",
        "one_sided_q_bh_within_pairwise_context",
        "one_sided_bh_signif_alpha_0_05",
        "one_sided_bh_signif_alpha_0_1",
        "one_sided_bh_within_pairwise_context_alpha_0_05",
        "one_sided_bh_within_pairwise_context_alpha_0_1",
        "higher_order_ci_above_wt",
        "higher_order_ci_nonoverlap_any_direction",
        "expected_for_plot",
        "observed_for_plot",
        "obs_over_exp_ratio",
        "which_signif_figure_logic",
        "in_expected_vs_observed_selected_logic",
        "visible_in_current_expected_vs_observed_axes_approx",
        "manuscript_subset",
        "manuscript_reference",
    ]
    return candidates.loc[
        candidates[alpha_0_1_col].fillna(False),
        columns,
    ].copy()


def write_manifest(
    output_root: Path,
    run_root: Path,
    result_dir: Path,
    method: str,
    candidates: pd.DataFrame,
    alpha: float,
    relaxed_alpha: float,
) -> None:
    """Purpose: record analysis provenance. Inputs: paths, method, candidate table, thresholds. Outputs: markdown manifest. Assumptions: concise provenance supports reviewer audit."""

    lines = [
        "# M=3 Expected-vs-Observed One-Sided WT Testing Manifest",
        "",
        "## Inputs",
        f"- Source run root: `{run_root}`",
        f"- Source result dir: `{result_dir}`",
        f"- Method: `{method}`",
        f"- Cohort: `{ES_LUAD_KEY}` only",
        "",
        "## Testing Family",
        "- Inferential unit: one M=3 triad-target double-mutant transition.",
        "- Included component pairwise contexts: `S+S` and `N+N`.",
        "- Excluded contexts: `N+E`, `S+A`, `A+A`, nonsmoking, and incomplete transition sets.",
        "- Contrast: target selection in the double-mutant background versus WT selection in the same M=3 fit.",
        "- Directional hypothesis: `H0: log(gamma_double / gamma_WT) <= 0` vs `H1: log(gamma_double / gamma_WT) > 0`.",
        "- Standard error rule: same null-facing CI-derived log-scale Wald SE used in Task 7.",
        "- Multiple testing: Benjamini-Hochberg applied once across all testable ES-LUAD `S+S` and `N+N` candidates.",
        f"- Primary alpha: `{alpha}`; relaxed alpha: `{relaxed_alpha}`.",
        "",
        "## Output Counts",
        f"- Candidate rows: `{len(candidates)}`",
        f"- Testable rows: `{int(candidates['one_sided_testable'].fillna(False).sum())}`",
        f"- BH significant at alpha 0.05: `{int(candidates['one_sided_bh_signif_alpha_0_05'].sum())}`",
        f"- BH significant at alpha 0.1: `{int(candidates['one_sided_bh_signif_alpha_0_1'].sum())}`",
        f"- Within-context BH significant at alpha 0.05: `{int(candidates['one_sided_bh_within_pairwise_context_alpha_0_05'].sum())}`",
        f"- Within-context BH significant at alpha 0.1: `{int(candidates['one_sided_bh_within_pairwise_context_alpha_0_1'].sum())}`",
        "",
        "## Software",
        f"- Python executable: `{sys.executable}`",
        f"- Python version: `{platform.python_version()}`",
        f"- pandas version: `{pd.__version__}`",
        f"- numpy version: `{np.__version__}`",
        "",
        "## Main Tables",
        "- `M3_expected_vs_observed_wt_one_sided_tests.csv`: row-level candidate tests and figure annotations.",
        "- `M3_expected_vs_observed_wt_one_sided_summary.csv`: combined and context-level counts.",
        "- `M3_expected_vs_observed_wt_one_sided_significant_interactions.csv`: BH-significant rows at alpha 0.1.",
        "- `M3_expected_vs_observed_wt_one_sided_within_context_significant_interactions.csv`: within-context BH-significant rows at alpha 0.1.",
    ]
    (output_root / "M3_expected_vs_observed_wt_one_sided_manifest.md").write_text(
        "\n".join(lines) + "\n"
    )


def run_expected_vs_observed_analysis(
    result_dir: Path,
    output_root: Path,
    run_root: Path,
    method: str,
    alpha: float,
    relaxed_alpha: float,
) -> dict[str, pd.DataFrame]:
    """Purpose: run the ES-LUAD expected-vs-observed one-sided WT workflow as a reusable pipeline step. Inputs: resolved result/output paths, method, alpha thresholds. Outputs: named result tables written to disk and returned in memory. Assumptions: M=3 task0 CSV exports are present in result_dir."""

    output_root.mkdir(parents=True, exist_ok=True)

    transitions = load_transitions(result_dir, method, MODEL_ORDER)
    candidates = build_candidate_tests(transitions)
    if candidates.empty:
        raise RuntimeError("No ES-LUAD S+S/N+N M=3 candidate tests were constructed.")

    candidates = add_one_sided_bh_columns(candidates, alpha, relaxed_alpha)
    summary = summarize_tests(candidates, alpha, relaxed_alpha)
    significant = select_significant_interactions(candidates, mode="combined")
    significant_within_context = select_significant_interactions(
        candidates,
        mode="within_pairwise_context",
    )

    candidates.to_csv(
        output_root / "M3_expected_vs_observed_wt_one_sided_tests.csv",
        index=False,
    )
    summary.to_csv(
        output_root / "M3_expected_vs_observed_wt_one_sided_summary.csv",
        index=False,
    )
    significant.to_csv(
        output_root / "M3_expected_vs_observed_wt_one_sided_significant_interactions.csv",
        index=False,
    )
    significant_within_context.to_csv(
        output_root
        / "M3_expected_vs_observed_wt_one_sided_within_context_significant_interactions.csv",
        index=False,
    )
    write_manifest(
        output_root,
        run_root,
        result_dir,
        method,
        candidates,
        alpha,
        relaxed_alpha,
    )

    return {
        "candidates": candidates,
        "summary": summary,
        "significant": significant,
        "significant_within_context": significant_within_context,
    }


def main() -> None:
    """Purpose: run the alternative expected-vs-observed testing workflow end to end. Inputs: CLI options. Outputs: task-local CSVs and manifest. Assumptions: M=3 CSV exports were generated from task0 model outputs."""

    args = parse_args()
    run_root, result_dir = resolve_result_dir(args.run_root, args.method)
    output_root = Path(args.output_root).resolve()
    results = run_expected_vs_observed_analysis(
        result_dir=result_dir,
        output_root=output_root,
        run_root=run_root,
        method=args.method,
        alpha=args.alpha,
        relaxed_alpha=args.relaxed_alpha,
    )

    print(f"Wrote M=3 expected-vs-observed one-sided outputs to {output_root}")
    print(results["summary"].to_string(index=False))


if __name__ == "__main__":
    main()
