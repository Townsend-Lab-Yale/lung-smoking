"""Compare run results against each other."""

from __future__ import annotations

import argparse
import ast
import math
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch


CODE_DIR = Path(__file__).resolve().parents[1]
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

import matrix_plotting as matrix_plotting_module


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BASELINE_ROOT = PROJECT_ROOT / "code" / "result_analysis"
COMPARISON_DIRNAME = "comparison"
PLOTS_DIRNAME = "plots"

M1_CLASSIFICATION_ORDER = [
    "preserved_significant_same_direction",
    "attenuated_same_direction",
    "reversed_significant",
    "reversed_non_significant",
    "newly_significant_same_direction",
    "newly_significant_reversal",
    "stable_non_significant",
]

EPI_CLASSIFICATION_ORDER = [
    "preserved_significant_same_direction",
    "attenuated_same_direction",
    "reversed_significant",
    "data_sparse_reversed_significant",
    "reversed_non_significant",
    "newly_significant_same_direction",
    "newly_significant_reversal",
    "data_sparse_newly_significant_reversal",
    "stable_non_significant",
]

CLASSIFICATION_COLORS = {
    "preserved_significant_same_direction": "#1b9e77",
    "attenuated_same_direction": "#7570b3",
    "reversed_significant": "#d95f02",
    "data_sparse_reversed_significant": "#f4a261",
    "reversed_non_significant": "#e7298a",
    "newly_significant_same_direction": "#66a61e",
    "newly_significant_reversal": "#e6ab02",
    "data_sparse_newly_significant_reversal": "#f6d55c",
    "stable_non_significant": "#999999",
}

PAIRWISE_MATRIX_KEYS = ("nonsmoking_plus", "smoking_plus")
SIGNIFICANT_REVERSAL_CLASSES = {
    "reversed_significant",
    "newly_significant_reversal",
    "data_sparse_reversed_significant",
    "data_sparse_newly_significant_reversal",
}
DATA_SPARSE_SIGNIFICANT_REVERSAL_CLASSES = {
    "data_sparse_reversed_significant",
    "data_sparse_newly_significant_reversal",
}


@dataclass(frozen=True)
class RunLayout:
    """Resolved locations for a result-comparison run.

    Inputs
    ------
    path:
        User-provided run reference.

    Outputs
    -------
    run_root:
        The root directory for the run.
    result_dir:
        The method-specific result directory under the run.

    Assumptions
    -----------
    The run has already been exported to CSV via ``numpy_to_csv.py``.
    """

    run_root: Path
    result_dir: Path


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments for a revision comparison.

    Returns
    -------
    argparse.Namespace
        Parsed arguments for baseline/revision paths and comparison options.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Compare a revision run against a baseline run and write "
            "task-scoped summary tables, plots, and comparisons.md."
        )
    )
    parser.add_argument(
        "--baseline-root",
        default=str(DEFAULT_BASELINE_ROOT),
        help=(
            "Baseline run reference. Can be a run root, a result directory, "
            "or the legacy code/result_analysis directory."
        ),
    )
    parser.add_argument(
        "--revision-root",
        required=True,
        help=(
            "Revision run reference. Can be a run root or a result directory. "
            "Comparison outputs are written beside this run."
        ),
    )
    parser.add_argument(
        "--method",
        choices=["variant", "cesR"],
        default="variant",
        help="Mutation-rate method to compare.",
    )
    parser.add_argument(
        "--key-a",
        default="smoking_plus",
        help="First key for M1 differential-comparison summaries.",
    )
    parser.add_argument(
        "--key-b",
        default="nonsmoking_plus",
        help="Second key for M1 differential-comparison summaries.",
    )
    parser.add_argument(
        "--model-orders",
        nargs="+",
        type=int,
        default=[2, 3],
        choices=[2, 3],
        help="Interaction model orders to compare.",
    )
    parser.add_argument(
        "--top-labels",
        type=int,
        default=10,
        help="Maximum number of flagged points to label per scatter plot.",
    )
    return parser.parse_args()


def resolve_run_layout(run_ref: str, method: str) -> RunLayout:
    """Resolve a run reference to its root and method-specific result dir.

    Parameters
    ----------
    run_ref:
        User-provided run path or run name.
    method:
        Mutation-rate method. Determines which result directory is required.

    Returns
    -------
    RunLayout
        Resolved run root and method-specific result directory.

    Assumptions
    -----------
    A valid run reference either exists directly or can be found relative to
    the project root or the ``output/`` directory.
    """

    method_dirname = f"{method}_results"
    raw_path = Path(run_ref)
    candidate_paths = []

    if raw_path.is_absolute():
        candidate_paths.append(raw_path)
    else:
        candidate_paths.extend(
            [
                raw_path,
                PROJECT_ROOT / raw_path,
                PROJECT_ROOT / "output" / raw_path,
                PROJECT_ROOT / "code" / "result_analysis" / raw_path,
            ]
        )

    for candidate in candidate_paths:
        resolved = _coerce_run_layout(candidate, method_dirname)
        if resolved is not None:
            return resolved

    raise FileNotFoundError(
        f"Could not resolve run reference '{run_ref}' for method '{method}'."
    )


def _coerce_run_layout(path: Path, method_dirname: str) -> RunLayout | None:
    """Try to interpret a path as a run root or result directory."""

    if not path.exists():
        return None

    path = path.resolve()

    if path.is_dir() and (path / method_dirname).is_dir():
        return RunLayout(run_root=path, result_dir=path / method_dirname)

    if path.is_dir() and path.name == method_dirname:
        return RunLayout(run_root=path.parent, result_dir=path)

    expected_csv = path / "M1_gene_gammas.csv"
    if path.is_dir() and expected_csv.exists():
        return RunLayout(run_root=path.parent, result_dir=path)

    legacy_result_dir = path / "code" / "result_analysis" / method_dirname
    if legacy_result_dir.is_dir():
        return RunLayout(
            run_root=(path / "code" / "result_analysis").resolve(),
            result_dir=legacy_result_dir.resolve(),
        )

    return None


def read_csv(result_dir: Path, file_name: str) -> pd.DataFrame:
    """Read a required result CSV from a method-specific result directory."""

    csv_path = result_dir / file_name
    if not csv_path.exists():
        raise FileNotFoundError(f"Missing required result file: {csv_path}")
    return pd.read_csv(csv_path)


def log_inner_join_coverage(
    label: str,
    left_df: pd.DataFrame,
    right_df: pd.DataFrame,
    entity_cols: list[str],
) -> None:
    """Log retained and dropped entities for an inner join."""

    def entity_tuples(df: pd.DataFrame) -> set[tuple[object, ...]]:
        return {
            tuple(row)
            for row in df.loc[:, entity_cols].drop_duplicates().itertuples(index=False, name=None)
        }

    def entity_label(entity: tuple[object, ...]) -> str:
        if len(entity) == 1:
            return str(entity[0])
        return " | ".join(str(item) for item in entity)

    left_entities = entity_tuples(left_df)
    right_entities = entity_tuples(right_df)
    shared_entities = left_entities & right_entities
    left_only = sorted(entity_label(entity) for entity in left_entities - right_entities)
    right_only = sorted(entity_label(entity) for entity in right_entities - left_entities)

    print(
        f"{label}: {len(left_entities)} left, {len(right_entities)} right, "
        f"{len(shared_entities)} shared ({len(left_only)} left-only, "
        f"{len(right_only)} right-only dropped)"
    )
    if left_only:
        print(f"{label} left-only entities dropped: {', '.join(left_only)}")
    if right_only:
        print(f"{label} right-only entities dropped: {', '.join(right_only)}")


def load_m1_results(result_dir: Path, method: str) -> pd.DataFrame:
    """Load M1 selection, flux, mutation-rate, and frequency summaries.

    Outputs
    -------
    pandas.DataFrame
        One row per ``key`` and ``gene`` with M1 estimates and sample counts.
    """

    gammas = read_csv(result_dir, "M1_gene_gammas.csv")
    fluxes = read_csv(result_dir, "M1_gene_fluxes.csv")
    mutation_rates = read_csv(result_dir, "mutation_rates.csv")
    samples = read_csv(result_dir, "M1_samples_per_combination.csv")

    mutation_rates = mutation_rates.loc[mutation_rates["method"] == method].copy()
    samples = samples.rename(columns={"(0,)": "wt_count", "(1,)": "mut_count"})
    samples["freq"] = samples["mut_count"] / (samples["wt_count"] + samples["mut_count"])

    merged = (
        gammas.merge(fluxes, on=["key", "gene"], how="left")
        .merge(mutation_rates[["method", "key", "gene", "mu"]],
               on=["method", "key", "gene"], how="left")
        .merge(samples[["key", "gene", "wt_count", "mut_count", "freq"]],
               on=["key", "gene"], how="left")
    )

    merged["method"] = method
    return merged


def compare_m1_runs(
    baseline_m1: pd.DataFrame,
    revision_m1: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare M1 estimates gene-by-gene across baseline and revision runs.

    Returns
    -------
    tuple[pandas.DataFrame, pandas.DataFrame]
        Full gene-level comparison table and key-level summary table.
    """

    merge_cols = ["method", "key", "gene"]
    log_inner_join_coverage(
        "M1 within-key comparison", baseline_m1, revision_m1, merge_cols
    )
    comparison = baseline_m1.merge(
        revision_m1,
        on=merge_cols,
        how="inner",
        suffixes=("_baseline", "_revision"),
    )

    for metric in ["gamma_mle", "flux_mle", "mu", "freq"]:
        baseline_col = f"{metric}_baseline"
        revision_col = f"{metric}_revision"
        ratio_col = f"{metric}_ratio"
        log_col = f"log10_{metric}_ratio"
        if baseline_col in comparison.columns and revision_col in comparison.columns:
            comparison[ratio_col] = safe_ratio(
                comparison[revision_col], comparison[baseline_col]
            )
            comparison[log_col] = safe_log10_ratio(
                comparison[revision_col], comparison[baseline_col]
            )

    comparison["abs_log10_gamma_ratio"] = comparison["log10_gamma_mle_ratio"].abs()
    comparison["twofold_gamma_change"] = (
        comparison["abs_log10_gamma_ratio"] >= math.log10(2)
    )

    summary_rows = []
    for key, key_df in comparison.groupby("key", sort=True):
        summary_rows.append(
            {
                "key": key,
                "n_genes": int(len(key_df)),
                "spearman_log10_gamma": safe_spearman(
                    np.log10(key_df["gamma_mle_baseline"]),
                    np.log10(key_df["gamma_mle_revision"]),
                ),
                "median_abs_log10_gamma_ratio": key_df[
                    "abs_log10_gamma_ratio"
                ].median(),
                "max_abs_log10_gamma_ratio": key_df["abs_log10_gamma_ratio"].max(),
                "genes_with_ge_2fold_gamma_change": int(
                    key_df["twofold_gamma_change"].sum()
                ),
                "median_abs_frequency_change": (
                    key_df["freq_revision"] - key_df["freq_baseline"]
                ).abs().median(),
            }
        )

    summary = pd.DataFrame(summary_rows)
    return comparison, summary


def build_m1_directional_table(
    m1_results: pd.DataFrame,
    key_a: str,
    key_b: str,
) -> pd.DataFrame:
    """Build within-run directional comparison between two M1 cohorts.

    Purpose
    -------
    Quantify whether a gene shows higher inferred selection in one cohort than
    another, and whether the difference is significant under the CI non-overlap
    rule used in the manuscript workflow.
    """

    left = m1_results.loc[m1_results["key"] == key_a].copy()
    right = m1_results.loc[m1_results["key"] == key_b].copy()

    shared_cols = [
        "method",
        "gene",
        "gamma_mle",
        "gamma_ci_low",
        "gamma_ci_high",
        "flux_mle",
        "freq",
        "wt_count",
        "mut_count",
    ]
    left = left[shared_cols].rename(columns=lambda col: f"{col}_{key_a}" if col != "gene" else col)
    right = right[shared_cols].rename(columns=lambda col: f"{col}_{key_b}" if col != "gene" else col)

    log_inner_join_coverage(
        f"M1 directional comparison ({key_a} vs {key_b})",
        left,
        right,
        ["gene"],
    )
    merged = left.merge(right, on="gene", how="inner")
    merged["comparison"] = f"{key_a}_vs_{key_b}"
    merged["gamma_ratio"] = safe_ratio(
        merged[f"gamma_mle_{key_a}"], merged[f"gamma_mle_{key_b}"]
    )
    merged["log2_gamma_ratio"] = safe_log2_ratio(
        merged[f"gamma_mle_{key_a}"], merged[f"gamma_mle_{key_b}"]
    )
    merged["direction"] = np.where(
        merged["gamma_ratio"] > 1,
        key_a,
        np.where(merged["gamma_ratio"] < 1, key_b, "neutral"),
    )
    merged["signif"] = (
        merged[f"gamma_ci_low_{key_a}"] > merged[f"gamma_ci_high_{key_b}"]
    ) | (
        merged[f"gamma_ci_low_{key_b}"] > merged[f"gamma_ci_high_{key_a}"]
    )
    return merged


def compare_directional_tables(
    baseline_df: pd.DataFrame,
    revision_df: pd.DataFrame,
    entity_cols: list[str],
    ratio_col: str,
    baseline_signif_col: str = "signif",
    comparison_label: str = "Directional comparison",
) -> pd.DataFrame:
    """Compare directional summaries between baseline and revision.

    Parameters
    ----------
    baseline_df, revision_df:
        Directional summary tables for the same type of analysis.
    entity_cols:
        Columns that uniquely identify a result entity.
    ratio_col:
        Column that determines direction via whether the ratio is above or
        below 1.

    Returns
    -------
    pandas.DataFrame
        Comparison table with directional-classification labels.
    """

    log_inner_join_coverage(
        comparison_label, baseline_df, revision_df, entity_cols
    )
    merged = baseline_df.merge(
        revision_df,
        on=entity_cols,
        how="inner",
        suffixes=("_baseline", "_revision"),
    )

    baseline_sign = np.sign(np.log(merged[f"{ratio_col}_baseline"]))
    revision_sign = np.sign(np.log(merged[f"{ratio_col}_revision"]))
    same_direction = (
        (baseline_sign == revision_sign)
        | (baseline_sign == 0)
        | (revision_sign == 0)
    )

    baseline_signif = merged[f"{baseline_signif_col}_baseline"].fillna(False)
    revision_signif = merged[f"{baseline_signif_col}_revision"].fillna(False)

    classification = np.full(len(merged), "stable_non_significant", dtype=object)
    classification[baseline_signif & revision_signif & same_direction] = (
        "preserved_significant_same_direction"
    )
    classification[baseline_signif & (~revision_signif) & same_direction] = (
        "attenuated_same_direction"
    )
    classification[baseline_signif & revision_signif & (~same_direction)] = (
        "reversed_significant"
    )
    classification[baseline_signif & (~revision_signif) & (~same_direction)] = (
        "reversed_non_significant"
    )
    classification[(~baseline_signif) & revision_signif & same_direction] = (
        "newly_significant_same_direction"
    )
    classification[(~baseline_signif) & revision_signif & (~same_direction)] = (
        "newly_significant_reversal"
    )

    merged["baseline_direction_sign"] = baseline_sign
    merged["revision_direction_sign"] = revision_sign
    merged[f"log2_{ratio_col}_baseline"] = np.log2(
        merged[f"{ratio_col}_baseline"].where(merged[f"{ratio_col}_baseline"] > 0)
    )
    merged[f"log2_{ratio_col}_revision"] = np.log2(
        merged[f"{ratio_col}_revision"].where(merged[f"{ratio_col}_revision"] > 0)
    )
    merged["same_direction"] = same_direction
    merged["classification"] = pd.Categorical(
        classification,
        categories=M1_CLASSIFICATION_ORDER,
        ordered=True,
    )
    return merged


def load_samples_long(result_dir: Path, model_order: int) -> pd.DataFrame:
    """Load sample counts per genotype state for M2 or M3."""

    samples = read_csv(result_dir, f"M{model_order}_samples_per_combination.csv")
    gene_cols = ["first_gene", "second_gene"]
    if model_order == 3:
        gene_cols.append("third_gene")

    state_cols = [col for col in samples.columns if col.startswith("(")]
    samples["gene_set"] = samples[gene_cols].astype(str).agg("_".join, axis=1)
    long_df = samples.melt(
        id_vars=["key", "gene_set"] + gene_cols,
        value_vars=state_cols,
        var_name="state",
        value_name="count",
    )
    long_df["state_key"] = long_df["state"].map(normalize_state_label)
    return long_df[["key", "gene_set", "state_key", "count"]]


def build_epistasis_table(
    result_dir: Path,
    method: str,
    model_order: int,
) -> pd.DataFrame:
    """Build raw epistasis comparisons from exported M2/M3 gamma tables.

    Purpose
    -------
    Reproduce the core comparison logic of the existing R workflow while
    avoiding the plotting-layer lower bound. Each returned row corresponds to
    one mutant-context comparison against the WT-origin transition for the same
    gene set and mutated gene.
    """

    gamma_df = read_csv(result_dir, f"M{model_order}_gene_gammas.csv")
    gamma_df = gamma_df.loc[gamma_df["method"] == method].copy()

    gene_cols = ["first_gene", "second_gene"]
    if model_order == 3:
        gene_cols.append("third_gene")

    gamma_df["gene_set"] = gamma_df[gene_cols].astype(str).agg("_".join, axis=1)
    parsed = gamma_df["mutation"].map(parse_transition)
    gamma_df["from_state"] = parsed.map(lambda item: state_key(item[0]))
    gamma_df["to_state"] = parsed.map(lambda item: state_key(item[1]))
    gamma_df["mutated_index"] = parsed.map(lambda item: item[2])
    gamma_df["mutated_gene"] = gamma_df.apply(
        lambda row: row[gene_cols[int(row["mutated_index"])]], axis=1
    )
    gamma_df["from_gt"] = gamma_df.apply(
        lambda row: state_to_genotype_label(
            [row[col] for col in gene_cols],
            parse_transition(row["mutation"])[0],
        ),
        axis=1,
    )

    sample_long = load_samples_long(result_dir, model_order)
    gamma_df = gamma_df.merge(
        sample_long.rename(columns={"count": "from_count"}),
        left_on=["key", "gene_set", "from_state"],
        right_on=["key", "gene_set", "state_key"],
        how="left",
    ).drop(columns=["state_key"])
    gamma_df = gamma_df.merge(
        sample_long.rename(columns={"count": "to_count"}),
        left_on=["key", "gene_set", "to_state"],
        right_on=["key", "gene_set", "state_key"],
        how="left",
    ).drop(columns=["state_key"])

    records = []
    group_cols = ["key", "gene_set", "mutated_gene"]
    for _, group in gamma_df.groupby(group_cols, sort=False):
        wt_rows = group.loc[group["from_gt"] == "WT"]
        if len(wt_rows) != 1:
            continue
        wt_row = wt_rows.iloc[0]
        mutant_rows = group.loc[group["from_gt"] != "WT"]

        for _, mutant_row in mutant_rows.iterrows():
            ratio = (
                mutant_row["gamma_mle"] / wt_row["gamma_mle"]
                if wt_row["gamma_mle"] > 0
                else np.nan
            )
            signif = (
                mutant_row["gamma_ci_low"] > wt_row["gamma_ci_high"]
                or wt_row["gamma_ci_low"] > mutant_row["gamma_ci_high"]
            )
            epistatic_gt = mutant_row["from_gt"]
            records.append(
                {
                    "method": method,
                    "model_order": model_order,
                    "key": mutant_row["key"],
                    "gene_set": mutant_row["gene_set"],
                    "mutated_gene": mutant_row["mutated_gene"],
                    "epistatic_gt": epistatic_gt,
                    "tested_combo": f"{epistatic_gt}_{mutant_row['mutated_gene']}",
                    "combo_name": combo_name(mutant_row["mutated_gene"], epistatic_gt),
                    "ratio": ratio,
                    "signif": bool(signif),
                    "wt_gamma_mle": wt_row["gamma_mle"],
                    "wt_gamma_ci_low": wt_row["gamma_ci_low"],
                    "wt_gamma_ci_high": wt_row["gamma_ci_high"],
                    "epi_gamma_mle": mutant_row["gamma_mle"],
                    "epi_gamma_ci_low": mutant_row["gamma_ci_low"],
                    "epi_gamma_ci_high": mutant_row["gamma_ci_high"],
                    "wt_from_count": wt_row["from_count"],
                    "epi_from_count": mutant_row["from_count"],
                    "epi_to_count": mutant_row["to_count"],
                }
            )

    return pd.DataFrame(records)


def build_selection_dicts(
    result_dir: Path,
    method: str,
    key: str,
    model_order: int = 2,
) -> tuple[dict[tuple[str, ...], dict[tuple[tuple[int, ...], tuple[int, ...]], float]],
           dict[tuple[str, ...], dict[tuple[tuple[int, ...], tuple[int, ...]], tuple[float, float]]]]:
    """Reconstruct nested selection dictionaries from exported gamma CSVs.

    Purpose
    -------
    Feed the existing matrix plotting code with the same dictionary structure
    used by the original analysis scripts, without depending on the original
    numpy result files.
    """

    gamma_df = read_csv(result_dir, f"M{model_order}_gene_gammas.csv")
    gamma_df = gamma_df.loc[
        (gamma_df["method"] == method) & (gamma_df["key"] == key)
    ].copy()

    gene_cols = ["first_gene", "second_gene"]
    if model_order == 3:
        gene_cols.append("third_gene")

    results = {}
    results_cis = {}
    for _, row in gamma_df.iterrows():
        genes = tuple(row[col] for col in gene_cols)
        from_state, to_state, _ = parse_transition(row["mutation"])
        transition = (from_state, to_state)
        results.setdefault(genes, {})[transition] = row["gamma_mle"]
        results_cis.setdefault(genes, {})[transition] = (
            row["gamma_ci_low"],
            row["gamma_ci_high"],
        )

    return results, results_cis


def build_pairwise_gene_order(m1_results: pd.DataFrame, fallback_keys: tuple[str, ...]) -> list[str]:
    """Choose a stable gene order for pairwise matrix plots.

    Assumptions
    -----------
    The baseline ``pan_data`` ordering should be used when available so that
    original and revision matrices are visually aligned.
    """

    if "pan_data" in set(m1_results["key"]):
        source = m1_results.loc[m1_results["key"] == "pan_data"].copy()
    else:
        source = m1_results.loc[m1_results["key"].isin(fallback_keys)].copy()
        source = source.sort_values(["key", "gamma_mle"], ascending=[True, False])

    source = source.sort_values("gamma_mle", ascending=False)
    return source["gene"].drop_duplicates().tolist()


def write_pairwise_ratio_matrices(
    baseline_layout: RunLayout,
    revision_layout: RunLayout,
    method: str,
    gene_order: list[str],
    plots_dir: Path,
) -> list[str]:
    """Write baseline and revision pairwise epistatic-ratio matrix plots."""

    created_files = []
    original_location = matrix_plotting_module.location_figures
    matrix_plotting_module.location_figures = str(plots_dir)

    try:
        baseline_results, baseline_cis = {}, {}
        revision_results, revision_cis = {}, {}
        for key in PAIRWISE_MATRIX_KEYS:
            baseline_results[key], baseline_cis[key] = build_selection_dicts(
                baseline_layout.result_dir, method, key, model_order=2
            )
            revision_results[key], revision_cis[key] = build_selection_dicts(
                revision_layout.result_dir, method, key, model_order=2
            )

        if all(baseline_results[key] for key in PAIRWISE_MATRIX_KEYS):
            baseline_name = "baseline_pairwise_epistatic_ratio_matrices.png"
            matrix_plotting_module.plot_epistatic_ratios_2_matrices(
                baseline_results["nonsmoking_plus"],
                baseline_cis["nonsmoking_plus"],
                baseline_results["smoking_plus"],
                baseline_cis["smoking_plus"],
                gene_order,
                plot_name=baseline_name,
                axis_label_size=12,
                axis_title_size=14,
            )
            created_files.append(baseline_name)

        if all(revision_results[key] for key in PAIRWISE_MATRIX_KEYS):
            revision_name = "revision_pairwise_epistatic_ratio_matrices.png"
            matrix_plotting_module.plot_epistatic_ratios_2_matrices(
                revision_results["nonsmoking_plus"],
                revision_cis["nonsmoking_plus"],
                revision_results["smoking_plus"],
                revision_cis["smoking_plus"],
                gene_order,
                plot_name=revision_name,
                axis_label_size=12,
                axis_title_size=14,
            )
            created_files.append(revision_name)
    finally:
        matrix_plotting_module.location_figures = original_location

    return created_files


def build_pairwise_status_matrix(
    key_df: pd.DataFrame,
    gene_order: list[str],
) -> pd.DataFrame:
    """Build an asymmetric pairwise status matrix from M2 comparison rows."""

    matrix = pd.DataFrame(
        data=np.nan,
        index=gene_order,
        columns=gene_order,
        dtype=object,
    )

    for _, row in key_df.iterrows():
        context_gene = row["epistatic_gt"]
        mutated_gene = row["mutated_gene"]
        if context_gene not in matrix.columns or mutated_gene not in matrix.index:
            continue
        matrix.loc[mutated_gene, context_gene] = row["classification"]

    return matrix


def plot_pairwise_status_matrix(
    status_matrix: pd.DataFrame,
    key: str,
    output_path: Path,
) -> None:
    """Plot a categorical pairwise status matrix with manuscript-relevant colors."""

    categories = EPI_CLASSIFICATION_ORDER
    category_to_code = {category: index for index, category in enumerate(categories)}
    numeric_matrix = status_matrix.apply(lambda column: column.map(category_to_code))
    numeric_matrix = numeric_matrix.astype(float)

    cmap = ListedColormap([CLASSIFICATION_COLORS[category] for category in categories])
    cmap.set_bad(np.array([0.7, 0.7, 0.7, 1.0]), 1.0)
    norm = BoundaryNorm(np.arange(-0.5, len(categories) + 0.5, 1), cmap.N)

    fig, ax = plt.subplots(figsize=(9, 8))
    mesh = ax.pcolormesh(
        np.arange(len(status_matrix.columns)),
        np.arange(len(status_matrix.index)),
        numeric_matrix.to_numpy()[::-1],
        cmap=cmap,
        norm=norm,
        shading="auto",
    )

    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks(np.arange(len(status_matrix.columns)))
    ax.set_xticklabels(
        status_matrix.columns, rotation=90, fontsize=12, style="italic"
    )
    ax.set_xlabel("Context (mutated gene in somatic genotype)", fontsize=14)

    ax.set_yticks(np.arange(len(status_matrix.index)))
    ax.set_yticklabels(status_matrix.index[::-1], fontsize=12, style="italic")
    ax.set_ylabel("Mutation under selection", fontsize=14)

    ax.set_xticks(np.arange(0.5, len(status_matrix.columns), 1), minor=True)
    ax.set_yticks(np.arange(0.5, len(status_matrix.index), 1), minor=True)
    ax.grid(which="minor", color="gray", linestyle="-", linewidth=1)
    ax.set_title(f"Pairwise directional/significance changes: {key}", fontsize=14)

    legend_handles = [
        Patch(facecolor=CLASSIFICATION_COLORS[category], edgecolor="none", label=category)
        for category in categories
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        frameon=False,
        fontsize=9,
    )

    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def compare_epistasis_runs(
    baseline_epi: pd.DataFrame,
    revision_epi: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compare epistatic interaction summaries across baseline and revision."""

    entity_cols = [
        "method",
        "model_order",
        "key",
        "gene_set",
        "mutated_gene",
        "epistatic_gt",
        "tested_combo",
        "combo_name",
    ]
    comparison = compare_directional_tables(
        baseline_epi,
        revision_epi,
        entity_cols=entity_cols,
        ratio_col="ratio",
        comparison_label=(
            f"M{int(baseline_epi['model_order'].iloc[0])} epistasis comparison"
        ),
    )
    data_sparse_mask = (
        comparison["classification"].isin(
            ["reversed_significant", "newly_significant_reversal"]
        )
        & comparison["epi_to_count_revision"].fillna(-1).eq(0)
        & comparison["epi_to_count_baseline"].fillna(0).gt(0)
    )
    comparison["data_sparse_reversal"] = data_sparse_mask
    classification = comparison["classification"].astype(str)
    classification = classification.mask(
        data_sparse_mask & classification.eq("reversed_significant"),
        "data_sparse_reversed_significant",
    )
    classification = classification.mask(
        data_sparse_mask & classification.eq("newly_significant_reversal"),
        "data_sparse_newly_significant_reversal",
    )
    comparison["classification"] = pd.Categorical(
        classification,
        categories=EPI_CLASSIFICATION_ORDER,
        ordered=True,
    )

    summary_rows = []
    for (model_order, key), key_df in comparison.groupby(
        ["model_order", "key"], sort=True
    ):
        counts = key_df["classification"].value_counts().reindex(
            EPI_CLASSIFICATION_ORDER, fill_value=0
        )
        summary_row = {
            "model_order": int(model_order),
            "key": key,
            "n_interactions_compared": int(len(key_df)),
            "n_baseline_significant": int(key_df["signif_baseline"].sum()),
            "n_revision_significant": int(key_df["signif_revision"].sum()),
        }
        summary_row.update({str(label): int(count) for label, count in counts.items()})
        summary_rows.append(summary_row)

    summary = pd.DataFrame(summary_rows)
    return comparison, summary


def summarize_epistasis_prevalence(
    comparison_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Summarize synergy/antagonism prevalence before and after a rerun.

    Purpose
    -------
    Separate the biological interpretation question (how much synergy versus
    antagonism is present in each cohort) from the rerun-stability question
    (whether specific interactions preserve direction/significance).

    Inputs
    ------
    comparison_df:
        Output of :func:`compare_epistasis_runs`.

    Outputs
    -------
    tuple[pandas.DataFrame, pandas.DataFrame]
        A long-form prevalence table for baseline and revision runs, plus a
        wide change table with revision-minus-baseline deltas.

    Assumptions
    -----------
    - Synergy is defined as epistatic ratio > 1.
    - Antagonism is defined as epistatic ratio < 1.
    - Prevalence is reported both over all tested interactions and within the
      subset that meets the existing CI non-overlap significance rule.
    """

    prevalence_rows = []
    for (model_order, key), key_df in comparison_df.groupby(
        ["model_order", "key"], sort=True
    ):
        for run_label, ratio_col, signif_col in [
            ("baseline", "ratio_baseline", "signif_baseline"),
            ("revision", "ratio_revision", "signif_revision"),
        ]:
            ratio = key_df[ratio_col]
            signif = key_df[signif_col].fillna(False)

            n_total = int(len(key_df))
            n_significant = int(signif.sum())
            n_synergistic_all = int((ratio > 1).sum())
            n_antagonistic_all = int((ratio < 1).sum())
            n_neutral_or_missing_all = (
                n_total - n_synergistic_all - n_antagonistic_all
            )
            n_synergistic_significant = int(((ratio > 1) & signif).sum())
            n_antagonistic_significant = int(((ratio < 1) & signif).sum())

            prevalence_rows.append(
                {
                    "model_order": int(model_order),
                    "key": key,
                    "run": run_label,
                    "n_interactions_total": n_total,
                    "n_significant": n_significant,
                    "n_synergistic_all": n_synergistic_all,
                    "n_antagonistic_all": n_antagonistic_all,
                    "n_neutral_or_missing_all": n_neutral_or_missing_all,
                    "synergistic_prevalence_all": (
                        n_synergistic_all / n_total if n_total > 0 else np.nan
                    ),
                    "antagonistic_prevalence_all": (
                        n_antagonistic_all / n_total if n_total > 0 else np.nan
                    ),
                    "n_synergistic_significant": n_synergistic_significant,
                    "n_antagonistic_significant": n_antagonistic_significant,
                    "synergistic_prevalence_significant": (
                        n_synergistic_significant / n_significant
                        if n_significant > 0 else np.nan
                    ),
                    "antagonistic_prevalence_significant": (
                        n_antagonistic_significant / n_significant
                        if n_significant > 0 else np.nan
                    ),
                }
            )

    prevalence_df = pd.DataFrame(prevalence_rows)
    if prevalence_df.empty:
        return prevalence_df, prevalence_df

    change_df = (
        prevalence_df.pivot(
            index=["model_order", "key"],
            columns="run",
            values=[
                "n_interactions_total",
                "n_significant",
                "n_synergistic_all",
                "n_antagonistic_all",
                "synergistic_prevalence_all",
                "antagonistic_prevalence_all",
                "n_synergistic_significant",
                "n_antagonistic_significant",
                "synergistic_prevalence_significant",
                "antagonistic_prevalence_significant",
            ],
        )
        .sort_index(axis=1)
        .reset_index()
    )
    change_df.columns = [
        "_".join(str(part) for part in col if part).strip("_")
        if isinstance(col, tuple) else col
        for col in change_df.columns
    ]

    for metric in [
        "n_significant",
        "n_synergistic_all",
        "n_antagonistic_all",
        "synergistic_prevalence_all",
        "antagonistic_prevalence_all",
        "n_synergistic_significant",
        "n_antagonistic_significant",
        "synergistic_prevalence_significant",
        "antagonistic_prevalence_significant",
    ]:
        baseline_col = f"{metric}_baseline"
        revision_col = f"{metric}_revision"
        if baseline_col in change_df.columns and revision_col in change_df.columns:
            change_df[f"{metric}_delta_revision_minus_baseline"] = (
                change_df[revision_col] - change_df[baseline_col]
            )

    return prevalence_df, change_df


def save_dataframe(df: pd.DataFrame, output_path: Path) -> None:
    """Write a DataFrame to CSV with reproducible row order."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)


def plot_m1_gamma_scatter(m1_comparison: pd.DataFrame, output_path: Path) -> None:
    """Plot baseline vs revision M1 gamma estimates by key."""

    keys = list(m1_comparison["key"].dropna().unique())
    n_cols = 2
    n_rows = math.ceil(len(keys) / n_cols) or 1
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(12, 4 * n_rows), squeeze=False
    )

    for index, key in enumerate(keys):
        ax = axes.flat[index]
        key_df = m1_comparison.loc[m1_comparison["key"] == key].copy()
        ax.scatter(
            key_df["gamma_mle_baseline"],
            key_df["gamma_mle_revision"],
            alpha=0.7,
            s=25,
            color="#1f77b4",
        )
        lower = min(key_df["gamma_mle_baseline"].min(), key_df["gamma_mle_revision"].min())
        upper = max(key_df["gamma_mle_baseline"].max(), key_df["gamma_mle_revision"].max())
        lower = max(lower, 1e-12)
        ax.plot([lower, upper], [lower, upper], linestyle="--", color="black", linewidth=1)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(key)
        ax.set_xlabel("Baseline M1 gamma")
        ax.set_ylabel("Revision M1 gamma")

    for ax in list(axes.flat)[len(keys):]:
        ax.axis("off")

    fig.suptitle("M1 gamma stability across reruns", y=0.98)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_directional_scatter(
    comparison_df: pd.DataFrame,
    x_col: str,
    y_col: str,
    label_col: str,
    title: str,
    x_label: str,
    y_label: str,
    output_path: Path,
    top_labels: int,
) -> None:
    """Plot directional stability for a comparison table."""

    if comparison_df.empty:
        return

    fig, ax = plt.subplots(figsize=(7, 6))
    for classification, subset in comparison_df.groupby(
        "classification", sort=False, observed=False
    ):
        ax.scatter(
            subset[x_col],
            subset[y_col],
            alpha=0.75,
            s=28,
            label=str(classification),
            color=CLASSIFICATION_COLORS.get(str(classification), "#555555"),
        )

    limits = [
        np.nanmin([comparison_df[x_col].min(), comparison_df[y_col].min()]),
        np.nanmax([comparison_df[x_col].max(), comparison_df[y_col].max()]),
    ]
    if np.isfinite(limits).all():
        ax.plot(limits, limits, linestyle="--", color="black", linewidth=1)
        ax.axhline(0, linestyle=":", color="gray", linewidth=1)
        ax.axvline(0, linestyle=":", color="gray", linewidth=1)

    flagged = comparison_df.loc[
        comparison_df["classification"].isin(SIGNIFICANT_REVERSAL_CLASSES)
    ].copy()
    flagged = flagged.sort_values(
        by=[x_col, y_col], key=np.abs, ascending=False
    ).head(top_labels)
    for _, row in flagged.iterrows():
        ax.annotate(
            str(row[label_col]),
            (row[x_col], row[y_col]),
            fontsize=8,
            alpha=0.85,
        )

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def summarize_classifications(
    comparison_df: pd.DataFrame,
    group_cols: list[str],
    order: list[str],
) -> pd.DataFrame:
    """Count directional-classification labels within grouping variables."""

    summary = (
        comparison_df.groupby(
            group_cols + ["classification"], dropna=False, observed=False
        )
        .size()
        .rename("n_results")
        .reset_index()
    )
    summary["classification"] = pd.Categorical(
        summary["classification"], categories=order, ordered=True
    )
    return summary.sort_values(group_cols + ["classification"])


def collect_flagged_findings(
    m1_directional: pd.DataFrame,
    epi_comparisons: Iterable[pd.DataFrame],
    key_a: str,
    key_b: str,
) -> pd.DataFrame:
    """Collect high-attention findings that may alter manuscript statements."""

    def build_flagged_frames(exclude_data_sparse: bool) -> list[pd.DataFrame]:
        flagged_frames = []

        if not m1_directional.empty:
            flagged_m1 = m1_directional.loc[
                m1_directional["classification"].isin(
                    ["reversed_significant", "newly_significant_reversal"]
                )
            ].copy()
            if not flagged_m1.empty:
                flagged_m1["analysis_type"] = "M1_directional"
                flagged_m1["entity_label"] = flagged_m1["gene"]
                flagged_m1["manuscript_attention"] = (
                    f"Smoking-differential selection between {key_a} and {key_b} changes direction."
                )
                flagged_frames.append(
                    flagged_m1[
                        [
                            "analysis_type",
                            "entity_label",
                            "classification",
                            "signif_baseline",
                            "signif_revision",
                            "gamma_ratio_baseline",
                            "gamma_ratio_revision",
                            f"gamma_mle_{key_a}_revision",
                            f"gamma_mle_{key_b}_revision",
                            "manuscript_attention",
                        ]
                    ].rename(
                        columns={
                            "gamma_ratio_baseline": "diff_sel_ratio_baseline",
                            "gamma_ratio_revision": "diff_sel_ratio_revision",
                        }
                    )
                )

        for epi_df in epi_comparisons:
            if epi_df.empty:
                continue
            flagged_epi = epi_df.loc[
                epi_df["classification"].isin(SIGNIFICANT_REVERSAL_CLASSES)
            ].copy()
            if exclude_data_sparse:
                flagged_epi = flagged_epi.loc[
                    ~flagged_epi["classification"].isin(
                        DATA_SPARSE_SIGNIFICANT_REVERSAL_CLASSES
                    )
                ]
            if flagged_epi.empty:
                continue
            flagged_epi["analysis_type"] = flagged_epi["model_order"].map(
                lambda order: f"M{int(order)}_epistasis"
            )
            flagged_epi["entity_label"] = flagged_epi["combo_name"]
            flagged_epi["manuscript_attention"] = np.select(
                [
                    flagged_epi["classification"].isin(
                        [
                            "reversed_significant",
                            "data_sparse_reversed_significant",
                        ]
                    ),
                    flagged_epi["classification"].isin(
                        [
                            "newly_significant_reversal",
                            "data_sparse_newly_significant_reversal",
                        ]
                    ),
                ],
                [
                    "Original significant epistatic direction reverses in the rerun.",
                    "A new significant opposite-direction interaction appears in the rerun.",
                ],
                default="Opposite-direction significant interaction detected in the rerun.",
            )
            flagged_frames.append(
                flagged_epi[
                    [
                        "analysis_type",
                        "entity_label",
                        "classification",
                        "signif_baseline",
                        "signif_revision",
                        "ratio_baseline",
                        "ratio_revision",
                        "manuscript_attention",
                    ]
                ].rename(
                    columns={
                        "ratio_baseline": "epistatic_ratio_baseline",
                        "ratio_revision": "epistatic_ratio_revision",
                    }
                )
            )

        return flagged_frames

    flagged_frames = []
    for flagging_rule, exclude_data_sparse in [
        ("all_significant_reversals", False),
        ("excluding_data_sparse_significant_reversals", True),
    ]:
        frames = build_flagged_frames(exclude_data_sparse=exclude_data_sparse)
        if not frames:
            continue
        combined = pd.concat(frames, ignore_index=True)
        combined["flagging_rule"] = flagging_rule
        flagged_frames.append(combined)

    if not flagged_frames:
        return pd.DataFrame(
            columns=[
                "analysis_type",
                "entity_label",
                "classification",
                "signif_baseline",
                "signif_revision",
                "diff_sel_ratio_baseline",
                "diff_sel_ratio_revision",
                "epistatic_ratio_baseline",
                "epistatic_ratio_revision",
                "manuscript_attention",
                "flagging_rule",
            ]
        )

    return pd.concat(flagged_frames, ignore_index=True)


def write_markdown_report(
    output_path: Path,
    baseline_layout: RunLayout,
    revision_layout: RunLayout,
    method: str,
    key_a: str,
    key_b: str,
    m1_summary: pd.DataFrame,
    m1_direction_summary: pd.DataFrame,
    epi_summaries: list[pd.DataFrame],
    epi_prevalence_changes: list[pd.DataFrame],
    flagged_findings: pd.DataFrame,
    matrix_output_files: list[str],
) -> None:
    """Write the task-scoped Markdown comparison report."""

    lines = [
        "# Revision Result Comparison",
        "",
        "## Run Metadata",
        f"- Generated: {datetime.now().isoformat(timespec='seconds')}",
        f"- Method: `{method}`",
        f"- Baseline root: `{baseline_layout.run_root}`",
        f"- Revision root: `{revision_layout.run_root}`",
        f"- M1 directional comparison: `{key_a}` vs `{key_b}`",
        "",
        "## Interpretation Rule",
        "- Original significant findings are evaluated first for directional consistency.",
        "- Same-direction loss of significance is treated as attenuation, not reversal.",
        "- Opposite-direction rerun estimates, especially if significant, are flagged as high attention because they can change manuscript statements.",
        "",
        "## M1 Within-Key Stability",
        dataframe_to_markdown(m1_summary.round(4)),
        "",
        "## M1 Directional Stability",
        dataframe_to_markdown(m1_direction_summary),
        "",
    ]

    for summary_df in epi_summaries:
        if summary_df.empty:
            continue
        model_order = int(summary_df["model_order"].iloc[0])
        lines.extend(
            [
                f"## M{model_order} Epistasis Directional Stability",
                dataframe_to_markdown(summary_df),
                "",
            ]
        )

    for prevalence_df in epi_prevalence_changes:
        if prevalence_df.empty:
            continue
        model_order = int(prevalence_df["model_order"].iloc[0])
        all_interactions_df = prevalence_df[
            [
                "key",
                "n_synergistic_all_baseline",
                "n_synergistic_all_revision",
                "n_synergistic_all_delta_revision_minus_baseline",
                "synergistic_prevalence_all_baseline",
                "synergistic_prevalence_all_revision",
                "n_antagonistic_all_baseline",
                "n_antagonistic_all_revision",
                "n_antagonistic_all_delta_revision_minus_baseline",
                "antagonistic_prevalence_all_baseline",
                "antagonistic_prevalence_all_revision",
            ]
        ].copy()
        all_interactions_df = all_interactions_df.rename(
            columns={
                "key": "cohort",
                "n_synergistic_all_baseline": "syn_n_base",
                "n_synergistic_all_revision": "syn_n_rev",
                "n_synergistic_all_delta_revision_minus_baseline": "syn_n_delta",
                "synergistic_prevalence_all_baseline": "syn_prev_base",
                "synergistic_prevalence_all_revision": "syn_prev_rev",
                "n_antagonistic_all_baseline": "ant_n_base",
                "n_antagonistic_all_revision": "ant_n_rev",
                "n_antagonistic_all_delta_revision_minus_baseline": "ant_n_delta",
                "antagonistic_prevalence_all_baseline": "ant_prev_base",
                "antagonistic_prevalence_all_revision": "ant_prev_rev",
            }
        ).round(4)
        significant_df = prevalence_df[
            [
                "key",
                "n_significant_baseline",
                "n_significant_revision",
                "synergistic_prevalence_significant_baseline",
                "synergistic_prevalence_significant_revision",
                "synergistic_prevalence_significant_delta_revision_minus_baseline",
                "antagonistic_prevalence_significant_baseline",
                "antagonistic_prevalence_significant_revision",
                "antagonistic_prevalence_significant_delta_revision_minus_baseline",
            ]
        ].copy()
        significant_df = significant_df.rename(
            columns={
                "key": "cohort",
                "n_significant_baseline": "n_sig_base",
                "n_significant_revision": "n_sig_rev",
                "synergistic_prevalence_significant_baseline": "syn_prev_sig_base",
                "synergistic_prevalence_significant_revision": "syn_prev_sig_rev",
                "synergistic_prevalence_significant_delta_revision_minus_baseline": "syn_prev_sig_delta",
                "antagonistic_prevalence_significant_baseline": "ant_prev_sig_base",
                "antagonistic_prevalence_significant_revision": "ant_prev_sig_rev",
                "antagonistic_prevalence_significant_delta_revision_minus_baseline": "ant_prev_sig_delta",
            }
        ).round(4)
        lines.extend(
            [
                f"## M{model_order} Synergy/Antagonism Prevalence",
                "Counts and prevalences are shown for each cohort before and after the rerun.",
                "Direction is defined by epistatic ratio > 1 (synergy) versus < 1 (antagonism).",
                "",
                "All tested pairwise interactions:",
                dataframe_to_markdown(all_interactions_df),
                "",
                "Significant pairwise interactions only:",
                dataframe_to_markdown(significant_df),
                "",
            ]
        )

    lines.append("## High-Attention Findings")
    if flagged_findings.empty:
        lines.append("- No opposite-direction significant reversals were detected by this comparison run.")
    else:
        for flagging_rule, header in [
            ("all_significant_reversals", "All significant reversals:"),
            (
                "excluding_data_sparse_significant_reversals",
                "Significant reversals after filtering data-sparse entries:",
            ),
        ]:
            lines.append(header)
            subset = flagged_findings.loc[
                flagged_findings["flagging_rule"] == flagging_rule
            ]
            if subset.empty:
                lines.append("- None.")
                lines.append("")
                continue
            for _, row in subset.iterrows():
                ratio_baseline = row.get("diff_sel_ratio_baseline", np.nan)
                if pd.isna(ratio_baseline):
                    ratio_baseline = row.get("epistatic_ratio_baseline", np.nan)
                ratio_revision = row.get("diff_sel_ratio_revision", np.nan)
                if pd.isna(ratio_revision):
                    ratio_revision = row.get("epistatic_ratio_revision", np.nan)
                extra_context = ""
                if row["analysis_type"] == "M1_directional" and pd.isna(ratio_revision):
                    gamma_key_a_revision = row.get(f"gamma_mle_{key_a}_revision", np.nan)
                    gamma_key_b_revision = row.get(f"gamma_mle_{key_b}_revision", np.nan)
                    extra_context = (
                        f" Revision gamma values: {key_a}="
                        f"{format_float(gamma_key_a_revision)}, {key_b}="
                        f"{format_float(gamma_key_b_revision)}."
                    )
                lines.append(
                    "- "
                    f"{row['analysis_type']}: `{row['entity_label']}` is classified as "
                    f"`{row['classification']}` "
                    f"(baseline ratio={format_float(ratio_baseline)}, "
                    f"revision ratio={format_float(ratio_revision)}). "
                    f"{row['manuscript_attention']}{extra_context}"
                )
            lines.append("")

    lines.extend(
        [
            "",
            "## Pairwise Matrix Outputs",
        ]
    )

    if matrix_output_files:
        for matrix_file in matrix_output_files:
            lines.append(f"- `{COMPARISON_DIRNAME}/{PLOTS_DIRNAME}/{matrix_file}`")
    else:
        lines.append("- Pairwise matrix outputs were not generated for this run.")

    lines.extend(
        [
            "",
            "## Output Files",
            f"- Detailed tables: `{COMPARISON_DIRNAME}/`",
            f"- Plots: `{COMPARISON_DIRNAME}/{PLOTS_DIRNAME}/`",
            "- This report is task-scoped and should be kept with the revision run folder.",
        ]
    )

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def dataframe_to_markdown(df: pd.DataFrame, max_rows: int = 20) -> str:
    """Convert a DataFrame to a compact Markdown table without extra deps."""

    if df.empty:
        return "_No rows available._"

    display_df = df.head(max_rows).copy().astype(object)
    display_df = display_df.where(pd.notna(display_df), "")
    headers = [str(col) for col in display_df.columns]
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for _, row in display_df.iterrows():
        values = [str(value) for value in row.tolist()]
        lines.append("| " + " | ".join(values) + " |")

    if len(df) > max_rows:
        lines.append("")
        lines.append(f"_Showing first {max_rows} of {len(df)} rows._")

    return "\n".join(lines)


def safe_ratio(numerator: pd.Series, denominator: pd.Series) -> pd.Series:
    """Compute a ratio while preserving missing values."""

    denominator = denominator.replace(0, np.nan)
    return numerator / denominator


def safe_log10_ratio(numerator: pd.Series, denominator: pd.Series) -> pd.Series:
    """Compute log10 revision/baseline ratios for positive inputs."""

    ratio = safe_ratio(numerator, denominator)
    ratio = ratio.where(ratio > 0)
    return np.log10(ratio)


def safe_log2_ratio(numerator: pd.Series, denominator: pd.Series) -> pd.Series:
    """Compute log2 revision/baseline ratios for positive inputs."""

    ratio = safe_ratio(numerator, denominator)
    ratio = ratio.where(ratio > 0)
    return np.log2(ratio)


def safe_spearman(x: pd.Series, y: pd.Series) -> float:
    """Compute Spearman correlation with graceful failure on degenerate data."""

    valid = pd.DataFrame({"x": x, "y": y}).dropna()
    if len(valid) < 2:
        return np.nan
    return valid["x"].rank().corr(valid["y"].rank())


def normalize_state_label(raw_label: str) -> str:
    """Normalize a state label such as ``(0, 1)`` to ``0, 1``."""

    state = ast.literal_eval(raw_label)
    if not isinstance(state, tuple):
        state = (state,)
    return state_key(state)


def parse_transition(mutation: str) -> tuple[tuple[int, ...], tuple[int, ...], int]:
    """Parse a mutation transition string.

    Returns
    -------
    tuple
        ``(from_state, to_state, mutated_index)``
    """

    from_state, to_state = ast.literal_eval(mutation)
    from_state = tuple(int(value) for value in from_state)
    to_state = tuple(int(value) for value in to_state)
    changed_indices = [
        idx for idx, (from_val, to_val) in enumerate(zip(from_state, to_state))
        if from_val != to_val
    ]
    if len(changed_indices) != 1:
        raise ValueError(f"Expected exactly one mutated index in transition {mutation}")
    return from_state, to_state, changed_indices[0]


def state_key(state: tuple[int, ...]) -> str:
    """Render a genotype state tuple as the CSV-compatible key."""

    return ", ".join(str(value) for value in state)


def state_to_genotype_label(genes: list[str], state: tuple[int, ...]) -> str:
    """Convert a binary state vector to a genotype label.

    Outputs
    -------
    str
        ``WT`` if no genes are mutated, otherwise genes joined by underscores.
    """

    selected = sorted(gene for gene, bit in zip(genes, state) if bit)
    return "WT" if not selected else "_".join(selected)


def combo_name(mutated_gene: str, epistatic_gt: str) -> str:
    """Create the manuscript-style label used in pairwise plots."""

    if epistatic_gt == "WT":
        return mutated_gene
    return f"{mutated_gene} [{epistatic_gt.replace('_', '+')}]"


def format_float(value: float) -> str:
    """Format a float for concise Markdown reporting."""

    if pd.isna(value):
        return "NA"
    if abs(value) >= 100 or (0 < abs(value) < 0.01):
        return f"{value:.2e}"
    return f"{value:.3f}"


def main() -> None:
    """Run the revision comparison workflow."""

    args = parse_args()
    baseline_layout = resolve_run_layout(args.baseline_root, args.method)
    revision_layout = resolve_run_layout(args.revision_root, args.method)

    comparison_dir = revision_layout.run_root / COMPARISON_DIRNAME
    plots_dir = comparison_dir / PLOTS_DIRNAME
    comparison_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    baseline_m1 = load_m1_results(baseline_layout.result_dir, args.method)
    revision_m1 = load_m1_results(revision_layout.result_dir, args.method)
    pairwise_gene_order = build_pairwise_gene_order(
        baseline_m1, fallback_keys=PAIRWISE_MATRIX_KEYS
    )

    m1_comparison, m1_summary = compare_m1_runs(baseline_m1, revision_m1)
    save_dataframe(m1_comparison, comparison_dir / "m1_gene_level_comparison.csv")
    save_dataframe(m1_summary, comparison_dir / "m1_within_key_summary.csv")
    plot_m1_gamma_scatter(m1_comparison, plots_dir / "m1_gamma_scatter.png")

    baseline_directional = build_m1_directional_table(
        baseline_m1, args.key_a, args.key_b
    )
    revision_directional = build_m1_directional_table(
        revision_m1, args.key_a, args.key_b
    )
    m1_directional = compare_directional_tables(
        baseline_directional,
        revision_directional,
        entity_cols=["gene", "comparison"],
        ratio_col="gamma_ratio",
        comparison_label="M1 directional stability across runs",
    )
    m1_direction_summary = summarize_classifications(
        m1_directional, group_cols=["comparison"], order=M1_CLASSIFICATION_ORDER
    )
    save_dataframe(
        m1_directional, comparison_dir / "m1_directional_comparison.csv"
    )
    save_dataframe(
        m1_direction_summary, comparison_dir / "m1_directional_summary.csv"
    )
    plot_directional_scatter(
        m1_directional,
        x_col="log2_gamma_ratio_baseline",
        y_col="log2_gamma_ratio_revision",
        label_col="gene",
        title=f"M1 directional stability: {args.key_a} vs {args.key_b}",
        x_label="Baseline log2 selection ratio",
        y_label="Revision log2 selection ratio",
        output_path=plots_dir / "m1_directional_scatter.png",
        top_labels=args.top_labels,
    )

    epi_comparisons = []
    epi_summaries = []
    epi_prevalence_changes = []
    matrix_output_files = write_pairwise_ratio_matrices(
        baseline_layout=baseline_layout,
        revision_layout=revision_layout,
        method=args.method,
        gene_order=pairwise_gene_order,
        plots_dir=plots_dir,
    )
    for model_order in args.model_orders:
        baseline_epi = build_epistasis_table(
            baseline_layout.result_dir, args.method, model_order
        )
        revision_epi = build_epistasis_table(
            revision_layout.result_dir, args.method, model_order
        )
        if baseline_epi.empty or revision_epi.empty:
            continue

        epi_comparison, epi_summary = compare_epistasis_runs(
            baseline_epi, revision_epi
        )
        epi_comparisons.append(epi_comparison)
        epi_summaries.append(epi_summary)

        save_dataframe(
            epi_comparison,
            comparison_dir / f"M{model_order}_epistasis_comparison.csv",
        )
        save_dataframe(
            epi_summary,
            comparison_dir / f"M{model_order}_epistasis_summary.csv",
        )
        epi_prevalence, epi_prevalence_change = summarize_epistasis_prevalence(
            epi_comparison
        )
        save_dataframe(
            epi_prevalence,
            comparison_dir / f"M{model_order}_epistasis_prevalence.csv",
        )
        save_dataframe(
            epi_prevalence_change,
            comparison_dir / f"M{model_order}_epistasis_prevalence_change.csv",
        )
        epi_prevalence_changes.append(epi_prevalence_change)

        for key, key_df in epi_comparison.groupby("key", sort=True):
            if model_order == 2:
                status_matrix = build_pairwise_status_matrix(
                    key_df=key_df,
                    gene_order=pairwise_gene_order,
                )
                status_matrix_path = (
                    comparison_dir / f"M2_{key}_status_matrix.csv"
                )
                status_matrix.to_csv(status_matrix_path, index=True)
                status_plot_name = f"M2_{key}_status_matrix.png"
                plot_pairwise_status_matrix(
                    status_matrix=status_matrix,
                    key=key,
                    output_path=plots_dir / status_plot_name,
                )
                matrix_output_files.append(status_plot_name)

            plot_directional_scatter(
                key_df,
                x_col="log2_ratio_baseline",
                y_col="log2_ratio_revision",
                label_col="combo_name",
                title=f"M{model_order} epistasis stability: {key}",
                x_label="Baseline log2 epistasis ratio",
                y_label="Revision log2 epistasis ratio",
                output_path=plots_dir / f"M{model_order}_{key}_epistasis_scatter.png",
                top_labels=args.top_labels,
            )

    flagged_findings = collect_flagged_findings(
        m1_directional, epi_comparisons, args.key_a, args.key_b
    )
    save_dataframe(flagged_findings, comparison_dir / "flagged_findings.csv")

    write_markdown_report(
        revision_layout.run_root / "comparisons.md",
        baseline_layout=baseline_layout,
        revision_layout=revision_layout,
        method=args.method,
        key_a=args.key_a,
        key_b=args.key_b,
        m1_summary=m1_summary,
        m1_direction_summary=m1_direction_summary,
        epi_summaries=epi_summaries,
        epi_prevalence_changes=epi_prevalence_changes,
        flagged_findings=flagged_findings,
        matrix_output_files=matrix_output_files,
    )


if __name__ == "__main__":
    main()
