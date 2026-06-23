"""Pairwise downstream multiple-testing analysis for exported epistasis results."""

from __future__ import annotations

import argparse
import ast
import math
import os
import platform
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from epistasis_testing_helpers import (
    DEFAULT_RUN_ROOT,
    DELTA_TOL,
    PROJECT_ROOT,
    add_wald_columns,
    bh_adjust,
    combo_name,
    load_transitions,
    resolve_result_dir,
)


DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "output" / "task7_multiple_testing"
MPLCONFIGDIR = Path(tempfile.gettempdir()) / "matplotlib-task7"
MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPLCONFIGDIR))
os.environ.setdefault("MPLBACKEND", "Agg")
CLASS_ORDER = ["both", "legacy_only", "bh_only", "neither", "untestable"]
CLASS_COLORS = {
    "both": "#1b9e77",
    "legacy_only": "#d95f02",
    "bh_only": "#7570b3",
    "neither": "#d9d9d9",
    "untestable": "#636363",
}
MATRIX_CLASSES = [
    ("ci_and_bh_005", "non-CI overlap and\nadj. p < 0.05", "non-CI overlap and\nadj. $p$ < 0.05", "#1b9e77"),
    ("ci_and_bh_010", "non-CI overlap and\nadj. p < 0.1", "non-CI overlap and\nadj. $p$ < 0.1", "#8fd17f"),
    ("ci_only", "non-CI overlap only", "non-CI overlap only", "#d95f02"),
    ("bh_005_only", "adj. p < 0.05 only", "adj. $p$ < 0.05 only", "#4c78a8"),
    ("neither", "neither", "neither", "#d9d9d9"),
    ("untestable", "untestable", "untestable", "#2f2f2f"),
]
MATRIX_ORDER = [item[0] for item in MATRIX_CLASSES]
MATRIX_LABELS = {item[0]: item[1] for item in MATRIX_CLASSES}
MATRIX_LEGEND_LABELS = {item[0]: item[2] for item in MATRIX_CLASSES}
MATRIX_COLORS = {item[0]: item[3] for item in MATRIX_CLASSES}


def parse_args() -> argparse.Namespace:
    """Purpose: parse CLI inputs. Inputs: argv. Outputs: analysis options. Assumption: run path is resolvable."""

    parser = argparse.ArgumentParser(
        description="Run pairwise BH correction on exported epistasis results."
    )
    parser.add_argument("--run-root", default=str(DEFAULT_RUN_ROOT))
    parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT))
    parser.add_argument("--method", default="variant", choices=["variant", "cesR"])
    parser.add_argument("--model-orders", nargs="+", type=int, default=[2], choices=[2, 3])
    parser.add_argument("--alpha", type=float, default=0.05)
    return parser.parse_args()


def build_contrasts(transitions: pd.DataFrame, model_order: int) -> pd.DataFrame:
    """Purpose: compare each mutant context to WT for the same acquired mutation. Inputs: transition table, model order. Outputs: one row per directional contrast. Assumption: each group has one WT-origin transition."""

    records = []
    for _, group in transitions.groupby(["key", "gene_set", "mutated_gene"], sort=False):
        wt_rows = group.loc[group["from_gt"] == "WT"]
        if len(wt_rows) != 1:
            continue
        wt = wt_rows.iloc[0]
        for _, ctx in group.loc[group["from_gt"] != "WT"].iterrows():
            records.append(
                {
                    "method": ctx["method"],
                    "model_order": model_order,
                    "key": ctx["key"],
                    "gene_set": ctx["gene_set"],
                    "mutated_gene": ctx["mutated_gene"],
                    "epistatic_gt": ctx["from_gt"],
                    "tested_combo": f"{ctx['from_gt']}_{ctx['mutated_gene']}",
                    "combo_name": combo_name(ctx["mutated_gene"], ctx["from_gt"]),
                    "context_size": int(ctx["context_size"]),
                    "contrast_type": "wt_reference",
                    "wt_gamma_mle": wt["gamma_mle"],
                    "wt_gamma_ci_low": wt["gamma_ci_low"],
                    "wt_gamma_ci_high": wt["gamma_ci_high"],
                    "wt_from_count": wt["from_count"],
                    "wt_to_count": wt["to_count"],
                    "epi_gamma_mle": ctx["gamma_mle"],
                    "epi_gamma_ci_low": ctx["gamma_ci_low"],
                    "epi_gamma_ci_high": ctx["gamma_ci_high"],
                    "epi_from_count": ctx["from_count"],
                    "epi_to_count": ctx["to_count"],
                    "ratio": ctx["gamma_mle"] / wt["gamma_mle"] if wt["gamma_mle"] > 0 and ctx["gamma_mle"] > 0 else np.nan,
                    "delta_hat": math.log(ctx["gamma_mle"] / wt["gamma_mle"]) if wt["gamma_mle"] > 0 and ctx["gamma_mle"] > 0 else np.nan,
                    "legacy_signif": bool(
                        (ctx["gamma_ci_low"] > wt["gamma_ci_high"])
                        or (wt["gamma_ci_low"] > ctx["gamma_ci_high"])
                    ),
                }
            )
    return pd.DataFrame(records)


def finalize_contrasts(contrasts: pd.DataFrame, alpha: float) -> pd.DataFrame:
    """Purpose: add BH calls and comparison labels. Inputs: contrast table with Wald columns, FDR threshold. Outputs: finalized contrast table. Assumption: BH family is all primary-testable contrasts."""

    contrasts = contrasts.copy()
    primary = contrasts["primary_testable"].fillna(False)
    conservative = contrasts["maxhalf_testable"].fillna(False)
    contrasts["p_raw_signif"] = False
    contrasts.loc[primary, "p_raw_signif"] = contrasts.loc[primary, "primary_p_raw"] <= alpha
    contrasts["q_bh"] = np.nan
    contrasts.loc[primary, "q_bh"] = bh_adjust(contrasts.loc[primary, "primary_p_raw"])
    contrasts["bh_signif"] = False
    contrasts.loc[primary, "bh_signif"] = contrasts.loc[primary, "q_bh"] <= alpha
    contrasts["q_bh_max_half_width"] = np.nan
    contrasts.loc[conservative, "q_bh_max_half_width"] = bh_adjust(contrasts.loc[conservative, "maxhalf_p_raw"])
    contrasts["bh_signif_max_half_width"] = False
    contrasts.loc[conservative, "bh_signif_max_half_width"] = contrasts.loc[conservative, "q_bh_max_half_width"] <= alpha
    contrasts["direction"] = np.where(
        contrasts["delta_hat"] > DELTA_TOL,
        "synergistic",
        np.where(contrasts["delta_hat"] < -DELTA_TOL, "antagonistic", "neutral"),
    )
    contrasts["legacy_vs_bh"] = np.where(
        ~primary,
        "untestable",
        np.where(
            contrasts["legacy_signif"] & contrasts["bh_signif"],
            "both",
            np.where(
                contrasts["legacy_signif"] & ~contrasts["bh_signif"],
                "legacy_only",
                np.where(~contrasts["legacy_signif"] & contrasts["bh_signif"], "bh_only", "neither"),
            ),
        ),
    )
    contrasts["abs_log_ratio"] = contrasts["delta_hat"].abs()
    return contrasts.sort_values(
        ["key", "bh_signif", "primary_p_raw", "abs_log_ratio"],
        ascending=[True, False, True, False],
        na_position="last",
    ).reset_index(drop=True)


def matrix_status(row: pd.Series) -> str:
    """Purpose: assign matrix display class. Inputs: one finalized contrast row. Outputs: matrix class key. Assumption: matrix tiers are based on CI non-overlap and BH-adjusted q-value thresholds."""

    if not row["primary_testable"] or pd.isna(row["q_bh"]):
        return "untestable"
    if row["legacy_signif"] and row["q_bh"] <= 0.05:
        return "ci_and_bh_005"
    if row["legacy_signif"] and row["q_bh"] <= 0.1:
        return "ci_and_bh_010"
    if row["legacy_signif"]:
        return "ci_only"
    if row["q_bh"] <= 0.05:
        return "bh_005_only"
    return "neither"


def build_summaries(contrasts: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Purpose: summarize Task 7 results. Inputs: finalized contrasts. Outputs: per-cohort, joint, and sensitivity summaries. Assumption: table already contains primary and conservative calls."""

    cohort_rows = []
    for key, key_df in contrasts.groupby("key", sort=True):
        bh_df = key_df.loc[key_df["bh_signif"]]
        cohort_rows.append(
            {
                "model_order": int(key_df["model_order"].iloc[0]),
                "key": key,
                "n_total_contrasts": int(len(key_df)),
                "n_testable": int(key_df["primary_testable"].sum()),
                "n_untestable": int((~key_df["primary_testable"]).sum()),
                "n_legacy_signif": int(key_df["legacy_signif"].sum()),
                "n_raw_signif": int(key_df["p_raw_signif"].sum()),
                "n_bh_signif": int(key_df["bh_signif"].sum()),
                "n_both": int((key_df["legacy_vs_bh"] == "both").sum()),
                "n_legacy_only": int((key_df["legacy_vs_bh"] == "legacy_only").sum()),
                "n_bh_only": int((key_df["legacy_vs_bh"] == "bh_only").sum()),
                "n_neither": int((key_df["legacy_vs_bh"] == "neither").sum()),
                "n_synergistic_bh": int((bh_df["direction"] == "synergistic").sum()),
                "n_antagonistic_bh": int((bh_df["direction"] == "antagonistic").sum()),
                "n_bh_signif_max_half_width": int(key_df["bh_signif_max_half_width"].sum()),
            }
        )

    bh_df = contrasts.loc[contrasts["bh_signif"]]
    joint = pd.DataFrame(
        [
            {
                "model_order": int(contrasts["model_order"].iloc[0]),
                "n_total_contrasts": int(len(contrasts)),
                "n_testable": int(contrasts["primary_testable"].sum()),
                "n_untestable": int((~contrasts["primary_testable"]).sum()),
                "n_legacy_signif": int(contrasts["legacy_signif"].sum()),
                "n_raw_signif": int(contrasts["p_raw_signif"].sum()),
                "n_bh_signif": int(contrasts["bh_signif"].sum()),
                "n_both": int((contrasts["legacy_vs_bh"] == "both").sum()),
                "n_legacy_only": int((contrasts["legacy_vs_bh"] == "legacy_only").sum()),
                "n_bh_only": int((contrasts["legacy_vs_bh"] == "bh_only").sum()),
                "n_neither": int((contrasts["legacy_vs_bh"] == "neither").sum()),
                "n_synergistic_bh": int((bh_df["direction"] == "synergistic").sum()),
                "n_antagonistic_bh": int((bh_df["direction"] == "antagonistic").sum()),
                "n_bh_signif_max_half_width": int(contrasts["bh_signif_max_half_width"].sum()),
            }
        ]
    )
    sensitivity = pd.DataFrame(
        [
            {
                "model_order": int(contrasts["model_order"].iloc[0]),
                "n_primary_bh_signif": int(contrasts["bh_signif"].sum()),
                "n_conservative_bh_signif": int(contrasts["bh_signif_max_half_width"].sum()),
                "n_preserved_under_conservative": int((contrasts["bh_signif"] & contrasts["bh_signif_max_half_width"]).sum()),
                "n_primary_only": int((contrasts["bh_signif"] & ~contrasts["bh_signif_max_half_width"]).sum()),
                "n_conservative_only": int((~contrasts["bh_signif"] & contrasts["bh_signif_max_half_width"]).sum()),
            }
        ]
    )
    return pd.DataFrame(cohort_rows), joint, sensitivity


def validate_legacy(contrasts: pd.DataFrame, comparison_dir: Path, method: str, model_order: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Purpose: check reconstructed legacy calls. Inputs: finalized contrasts, comparison dir, method, model order. Outputs: summary and mismatch table. Assumption: source comparison CSV uses the same directional identifiers."""

    path = comparison_dir / f"M{model_order}_epistasis_comparison.csv"
    empty = pd.DataFrame([{"model_order": model_order, "n_shared_rows": 0, "n_exact_matches": 0, "n_mismatches": 0}])
    if not path.exists():
        return empty, pd.DataFrame()

    existing = pd.read_csv(path)
    keys = ["method", "model_order", "key", "gene_set", "mutated_gene", "epistatic_gt"]
    if not set(keys).issubset(existing.columns):
        return empty, pd.DataFrame()

    existing = existing.loc[(existing["method"] == method) & (existing["model_order"] == model_order), keys + ["signif_revision"]]
    merged = contrasts.merge(
        existing.rename(columns={"signif_revision": "legacy_signif_existing"}),
        on=keys,
        how="inner",
    )
    merged["legacy_match"] = merged["legacy_signif"] == merged["legacy_signif_existing"]
    summary = pd.DataFrame(
        [
            {
                "model_order": model_order,
                "n_shared_rows": int(len(merged)),
                "n_exact_matches": int(merged["legacy_match"].sum()),
                "n_mismatches": int((~merged["legacy_match"]).sum()),
            }
        ]
    )
    return summary, merged.loc[~merged["legacy_match"]].copy()


def plot_joint_diagnostics(contrasts: pd.DataFrame, plots_dir: Path, prefix: str, alpha: float) -> None:
    """Purpose: write histogram, QQ, and volcano plots. Inputs: finalized contrasts, plot dir, prefix, alpha. Outputs: PNG diagnostics. Assumption: only primary-testable contrasts are plotted."""

    import matplotlib.pyplot as plt

    data = contrasts.loc[contrasts["primary_testable"], "primary_p_raw"].dropna().sort_values()
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(data, bins=25, color="#4c78a8", edgecolor="white")
    ax.set(xlabel="Raw Wald p-value", ylabel="Number of contrasts", title="Task 7 raw p-value histogram")
    fig.tight_layout()
    fig.savefig(plots_dir / f"{prefix}_joint_pvalue_histogram.png", dpi=300)
    plt.close(fig)

    if not data.empty:
        observed = -np.log10(np.clip(data.to_numpy(dtype=float), 1e-300, 1.0))
        expected = -np.log10((np.arange(1, len(data) + 1) - 0.5) / len(data))
        fig, ax = plt.subplots(figsize=(5.5, 5.5))
        ax.scatter(expected, observed, s=14, alpha=0.75, color="#1b9e77")
        limit = max(expected.max(), observed.max())
        ax.plot([0, limit], [0, limit], linestyle="--", color="black", linewidth=1)
        ax.set(xlabel="Expected -log10(p)", ylabel="Observed -log10(p)", title="Task 7 QQ plot")
        fig.tight_layout()
        fig.savefig(plots_dir / f"{prefix}_joint_qq.png", dpi=300)
        plt.close(fig)

    volcano = contrasts.loc[contrasts["primary_testable"]].copy()
    if not volcano.empty:
        volcano["neg_log10_p"] = -np.log10(np.clip(volcano["primary_p_raw"], 1e-300, 1.0))
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.scatter(
            volcano["delta_hat"],
            volcano["neg_log10_p"],
            s=18,
            alpha=0.75,
            c=np.where(volcano["bh_signif"], "#d95f02", "#4c78a8"),
        )
        ax.axvline(0.0, color="black", linestyle="--", linewidth=1)
        ax.axhline(-math.log10(alpha), color="gray", linestyle=":", linewidth=1)
        ax.set(xlabel="Log epistatic ratio", ylabel="-log10(raw p-value)", title="Task 7 volcano plot")
        fig.tight_layout()
        fig.savefig(plots_dir / f"{prefix}_joint_volcano.png", dpi=300)
        plt.close(fig)


def plot_pairwise_matrices(contrasts: pd.DataFrame, output_root: Path, plots_dir: Path, prefix: str) -> None:
    """Purpose: write M2 legacy-vs-BH status matrices. Inputs: finalized contrasts and output paths. Outputs: matrix CSVs and PNGs. Assumption: only single-gene backgrounds are shown."""

    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.patches import Patch

    genes = sorted(set(contrasts["mutated_gene"]).union(set(contrasts["epistatic_gt"])))
    genes = [gene for gene in genes if gene != "WT" and "_" not in gene]
    cmap = ListedColormap([MATRIX_COLORS[label] for label in MATRIX_ORDER])
    norm = BoundaryNorm(np.arange(-0.5, len(MATRIX_ORDER) + 0.5, 1), cmap.N)
    code = {label: i for i, label in enumerate(MATRIX_ORDER)}

    for key in sorted(contrasts["key"].unique()):
        matrix = pd.DataFrame(index=genes, columns=genes, dtype=object).fillna("untestable")
        for _, row in contrasts.loc[contrasts["key"] == key].iterrows():
            if row["epistatic_gt"] == "WT" or "_" in row["epistatic_gt"]:
                continue
            if row["mutated_gene"] in matrix.index and row["epistatic_gt"] in matrix.columns:
                matrix.loc[row["mutated_gene"], row["epistatic_gt"]] = matrix_status(row)
        matrix.rename(index=MATRIX_LABELS, columns=MATRIX_LABELS)
        matrix_display = matrix.replace(MATRIX_LABELS)
        matrix_display.to_csv(output_root / f"{prefix}_{key}_legacy_vs_bh_matrix.csv", index=True)

        fig, ax = plt.subplots(figsize=(9, 8))
        ax.pcolormesh(
            np.arange(len(matrix.columns)),
            np.arange(len(matrix.index)),
            matrix.apply(lambda col: col.map(code)).to_numpy(dtype=float)[::-1],
            cmap=cmap,
            norm=norm,
            shading="auto",
        )
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks(np.arange(len(matrix.columns)))
        ax.set_xticklabels(matrix.columns, rotation=90, fontsize=10, style="italic")
        ax.set_yticks(np.arange(len(matrix.index)))
        ax.set_yticklabels(matrix.index[::-1], fontsize=10, style="italic")
        ax.set(xlabel="Context gene", ylabel="Mutation under selection", title=f"Task 7 CI vs BH status: {key}")
        ax.set_xticks(np.arange(0.5, len(matrix.columns), 1), minor=True)
        ax.set_yticks(np.arange(0.5, len(matrix.index), 1), minor=True)
        ax.grid(which="minor", color="white", linestyle="-", linewidth=1)
        ax.legend(
            handles=[
                Patch(facecolor=MATRIX_COLORS[label], edgecolor="none", label=MATRIX_LEGEND_LABELS[label])
                for label in MATRIX_ORDER
            ],
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
        )
        fig.tight_layout()
        fig.savefig(plots_dir / f"{prefix}_{key}_legacy_vs_bh_matrix.png", dpi=300)
        plt.close(fig)


def write_manifest(output_root: Path, run_root: Path, result_dir: Path, method: str, model_orders: list[int], summaries: list[dict[str, object]]) -> None:
    """Purpose: record run metadata. Inputs: output path, source paths, method, model orders, summaries. Outputs: markdown manifest. Assumption: concise provenance is sufficient."""

    lines = [
        "# Task 7 Run Manifest",
        "",
        "## Inputs",
        f"- Source run root: `{run_root}`",
        f"- Source result dir: `{result_dir}`",
        f"- Method: `{method}`",
        f"- Requested model orders: `{', '.join(str(x) for x in model_orders)}`",
        "",
        "## Fixed analysis choices",
        "- Primary family: all testable contrasts within model order, corrected jointly across `smoking_plus` and `nonsmoking_plus`.",
        "- Primary inferential scale: log epistatic ratio.",
        "- Standard error rule: null-facing 95% CI half-width on the log scale.",
        "- Multiple-testing procedure: Benjamini-Hochberg.",
        "- Primary exclusion rule: any contrast with nonpositive gamma bound or invalid interval.",
        "- Sensitivity analysis: repeat BH using the conservative max-half-width CI-to-SE rule.",
        "",
        "## Software",
        f"- Python executable: `{sys.executable}`",
        f"- Python version: `{platform.python_version()}`",
        f"- pandas version: `{pd.__version__}`",
        f"- numpy version: `{np.__version__}`",
        "",
        "## Model-order outcomes",
    ]
    for summary in summaries:
        if summary["status"] == "completed":
            lines.extend(
                [
                    f"- M{summary['model_order']}: completed.",
                    f"  Testable contrasts: `{summary['n_testable']}` / `{summary['n_total']}`.",
                    f"  BH-significant contrasts: `{summary['n_bh_signif']}`.",
                    f"  Legacy/BH overlap: `{summary['n_both']}` both, `{summary['n_legacy_only']}` legacy-only, `{summary['n_bh_only']}` BH-only.",
                ]
            )
        else:
            lines.append(f"- M{summary['model_order']}: skipped (`{summary['status']}`: {summary['message']}).")
    (output_root / "task7_run_manifest.md").write_text("\n".join(lines) + "\n")


def run_model_order(result_dir: Path, comparison_dir: Path, output_root: Path, plots_dir: Path, method: str, model_order: int, alpha: float) -> dict[str, object]:
    """Purpose: execute one model order end to end. Inputs: source dirs, output dirs, method, order, alpha. Outputs: files on disk and manifest summary row. Assumption: requested CSV exports exist."""

    gamma_path = result_dir / f"M{model_order}_gene_gammas.csv"
    sample_path = result_dir / f"M{model_order}_samples_per_combination.csv"
    if not gamma_path.exists() or not sample_path.exists():
        return {"model_order": model_order, "status": "missing_input", "message": f"Missing {gamma_path.name} or {sample_path.name}"}

    contrasts = build_contrasts(load_transitions(result_dir, method, model_order), model_order)
    if contrasts.empty:
        return {"model_order": model_order, "status": "no_contrasts", "message": "No WT-reference contrasts were reconstructed."}

    contrasts = finalize_contrasts(
        add_wald_columns(add_wald_columns(contrasts, "null_facing", "primary"), "max_half_width", "maxhalf"),
        alpha,
    )
    excluded = contrasts.loc[~contrasts["primary_testable"]].copy()
    comparison = contrasts[
        [
            "method",
            "model_order",
            "key",
            "gene_set",
            "mutated_gene",
            "epistatic_gt",
            "tested_combo",
            "combo_name",
            "ratio",
            "delta_hat",
            "direction",
            "legacy_signif",
            "primary_testable",
            "primary_exclusion_reason",
            "primary_se_delta",
            "primary_z_wald",
            "primary_p_raw",
            "q_bh",
            "p_raw_signif",
            "bh_signif",
            "legacy_vs_bh",
            "maxhalf_se_delta",
            "maxhalf_p_raw",
            "q_bh_max_half_width",
            "bh_signif_max_half_width",
            "wt_gamma_mle",
            "wt_gamma_ci_low",
            "wt_gamma_ci_high",
            "epi_gamma_mle",
            "epi_gamma_ci_low",
            "epi_gamma_ci_high",
            "wt_from_count",
            "epi_from_count",
            "epi_to_count",
        ]
    ].copy()
    cohort_summary, joint_summary, sensitivity = build_summaries(contrasts)
    legacy_summary, legacy_mismatches = validate_legacy(contrasts, comparison_dir, method, model_order)

    prefix = f"M{model_order}"
    contrasts.to_csv(output_root / f"{prefix}_primary_contrasts.csv", index=False)
    excluded.to_csv(output_root / f"{prefix}_boundary_untestable_contrasts.csv", index=False)
    comparison.to_csv(output_root / f"{prefix}_comparison_to_legacy.csv", index=False)
    cohort_summary.to_csv(output_root / f"{prefix}_summary_by_cohort.csv", index=False)
    joint_summary.to_csv(output_root / f"{prefix}_summary_joint.csv", index=False)
    sensitivity.to_csv(output_root / f"{prefix}_sensitivity_max_half_width_summary.csv", index=False)
    legacy_summary.to_csv(output_root / f"{prefix}_legacy_validation_summary.csv", index=False)
    if not legacy_mismatches.empty:
        legacy_mismatches.to_csv(output_root / f"{prefix}_legacy_validation_mismatches.csv", index=False)

    plot_joint_diagnostics(contrasts, plots_dir, prefix, alpha)
    if model_order == 2:
        plot_pairwise_matrices(contrasts, output_root, plots_dir, prefix)

    joint = joint_summary.iloc[0]
    return {
        "model_order": model_order,
        "status": "completed",
        "n_total": int(joint["n_total_contrasts"]),
        "n_testable": int(joint["n_testable"]),
        "n_bh_signif": int(joint["n_bh_signif"]),
        "n_both": int(joint["n_both"]),
        "n_legacy_only": int(joint["n_legacy_only"]),
        "n_bh_only": int(joint["n_bh_only"]),
    }


def main() -> None:
    """Purpose: run Task 7. Inputs: CLI options. Outputs: task-local CSVs, plots, and manifest. Assumption: analysis is strictly downstream of the source run."""

    args = parse_args()
    run_root, result_dir = resolve_result_dir(args.run_root, args.method)
    output_root = Path(args.output_root).resolve()
    plots_dir = output_root / "plots"
    output_root.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    summaries = []
    for model_order in args.model_orders:
        summary = run_model_order(
            result_dir=result_dir,
            comparison_dir=run_root / "comparison",
            output_root=output_root,
            plots_dir=plots_dir,
            method=args.method,
            model_order=model_order,
            alpha=args.alpha,
        )
        summaries.append(summary)
        if summary["status"] == "completed":
            print(f"M{model_order}: {summary['n_bh_signif']} BH-significant contrasts out of {summary['n_testable']} testable.")
        else:
            print(f"M{model_order}: skipped ({summary['status']} - {summary['message']}).")

    write_manifest(output_root, run_root, result_dir, args.method, args.model_orders, summaries)


if __name__ == "__main__":
    main()
