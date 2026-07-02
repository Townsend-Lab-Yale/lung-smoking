"""Shared helpers for epistasis Wald/BH testing.

Purpose:
    Centralize path resolution, transition loading, CI-derived Wald testing,
    and BH adjustment used by pairwise and higher-order epistasis analyses.

Inputs/outputs:
    Functions in this module operate on exported task0 CSVs and pandas
    DataFrames, returning normalized transition or contrast tables.

Assumptions:
    - Exported gamma confidence intervals are 95% intervals on the original
      gamma scale.
    - The log-ratio Wald approximation treats the two compared log-gamma
      estimates as independent.
"""

from __future__ import annotations

import ast
import math
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RUN_ROOT = PROJECT_ROOT / "output"
Z_975 = 1.959963984540054
DELTA_TOL = 1e-12


def resolve_result_dir(run_ref: str, method: str) -> tuple[Path, Path]:
    """Purpose: resolve run root/result dir. Inputs: run reference, method. Outputs: run root and method dir. Assumption: standard project layout."""

    method_dir = f"{method}_results"
    raw = Path(run_ref)
    candidates = [raw] if raw.is_absolute() else [
        raw,
        PROJECT_ROOT / raw,
        PROJECT_ROOT / "output" / raw,
        PROJECT_ROOT / "code" / "result_analysis" / raw,
    ]
    for candidate in candidates:
        if not candidate.exists():
            continue
        path = candidate.resolve()
        if (path / method_dir).is_dir():
            return path, path / method_dir
        if path.name == method_dir and path.is_dir():
            return path.parent, path
        if path.is_dir() and (path / "M2_gene_gammas.csv").exists():
            return path.parent, path
    raise FileNotFoundError(
        f"Could not resolve run reference '{run_ref}' for method '{method}'."
    )


def parse_transition(mutation: str) -> tuple[tuple[int, ...], tuple[int, ...], int]:
    """Purpose: parse one transition. Inputs: exported mutation string. Outputs: from-state, to-state, changed index. Assumption: exactly one coordinate changes."""

    from_state, to_state = ast.literal_eval(mutation)
    from_state = tuple(int(x) for x in from_state)
    to_state = tuple(int(x) for x in to_state)
    changed = [i for i, (a, b) in enumerate(zip(from_state, to_state)) if a != b]
    if len(changed) != 1:
        raise ValueError(f"Expected exactly one mutated index in transition {mutation}")
    return from_state, to_state, changed[0]


def state_key(state: tuple[int, ...]) -> str:
    """Purpose: normalize genotype state. Inputs: binary tuple. Outputs: CSV join key. Assumption: state order matches gene columns."""

    return ", ".join(str(x) for x in state)


def genotype_label(genes: list[str], state: tuple[int, ...]) -> str:
    """Purpose: build genotype label. Inputs: ordered genes, binary state. Outputs: WT or underscore-joined genes. Assumption: labels should be order-invariant."""

    hits = sorted(gene for gene, bit in zip(genes, state) if bit)
    return "WT" if not hits else "_".join(hits)


def combo_name(mutated_gene: str, epistatic_gt: str) -> str:
    """Purpose: create compact contrast label. Inputs: mutated gene, background genotype. Outputs: reviewer-facing label. Assumption: underscores in backgrounds denote combined contexts."""

    return mutated_gene if epistatic_gt == "WT" else f"{mutated_gene} [{epistatic_gt.replace('_', '+')}]"


def load_transitions(result_dir: Path, method: str, model_order: int) -> pd.DataFrame:
    """Purpose: read exported transitions with counts. Inputs: result dir, method, model order. Outputs: one row per transition. Assumption: gamma and sample CSVs come from the same run."""

    gene_cols = ["first_gene", "second_gene"] + (["third_gene"] if model_order == 3 else [])
    gammas = pd.read_csv(result_dir / f"M{model_order}_gene_gammas.csv")
    gammas = gammas.loc[gammas["method"] == method].copy()
    if gammas.empty:
        raise ValueError(f"No rows found for method '{method}' in M{model_order}_gene_gammas.csv.")

    samples = pd.read_csv(result_dir / f"M{model_order}_samples_per_combination.csv")
    state_cols = [col for col in samples.columns if col.startswith("(")]
    sample_long = (
        samples.assign(gene_set=samples[gene_cols].astype(str).agg("_".join, axis=1))
        .melt(
            id_vars=["key", "gene_set"] + gene_cols,
            value_vars=state_cols,
            var_name="state_label",
            value_name="count",
        )
        .assign(state_key=lambda df: df["state_label"].map(lambda x: state_key(ast.literal_eval(x))))
        [["key", "gene_set", "state_key", "count"]]
    )

    parsed = gammas["mutation"].map(parse_transition)
    gammas["gene_set"] = gammas[gene_cols].astype(str).agg("_".join, axis=1)
    gammas["from_state"] = parsed.map(lambda x: x[0])
    gammas["to_state"] = parsed.map(lambda x: x[1])
    gammas["from_state_key"] = gammas["from_state"].map(state_key)
    gammas["to_state_key"] = gammas["to_state"].map(state_key)
    gammas["mutated_index"] = parsed.map(lambda x: x[2])
    gammas["mutated_gene"] = gammas.apply(
        lambda row: row[gene_cols[int(row["mutated_index"])]],
        axis=1,
    )
    gammas["from_gt"] = gammas.apply(
        lambda row: genotype_label([row[col] for col in gene_cols], row["from_state"]),
        axis=1,
    )
    gammas["context_size"] = gammas["from_state"].map(sum).astype(int)
    gammas = gammas.merge(
        sample_long.rename(columns={"state_key": "from_state_key", "count": "from_count"}),
        on=["key", "gene_set", "from_state_key"],
        how="left",
    ).merge(
        sample_long.rename(columns={"state_key": "to_state_key", "count": "to_count"}),
        on=["key", "gene_set", "to_state_key"],
        how="left",
    )

    return gammas[
        [
            "method",
            "key",
            "gene_set",
            "mutation",
            "mutated_gene",
            "from_gt",
            "context_size",
            "gamma_mle",
            "gamma_ci_low",
            "gamma_ci_high",
            "from_count",
            "to_count",
        ]
    ].copy()


def add_wald_columns(contrasts: pd.DataFrame, mode: str, prefix: str) -> pd.DataFrame:
    """Purpose: reconstruct Wald SEs/p-values from CIs. Inputs: contrast table, CI rule, prefix. Outputs: table with test columns appended. Assumption: log-gamma estimates are approximately independent."""

    def interval_reason(mle: float, lo: float, hi: float) -> str | None:
        if any(pd.isna(x) for x in [mle, lo, hi]):
            return "missing_value"
        if mle <= 0:
            return "nonpositive_mle"
        if lo <= 0 or hi <= 0:
            return "nonpositive_bound"
        if lo > mle or hi < mle or lo > hi:
            return "invalid_interval"
        return None

    def half_width(mle: float, lo: float, hi: float, side: str) -> float:
        lower = math.log(mle) - math.log(lo)
        upper = math.log(hi) - math.log(mle)
        return {"lower": lower, "upper": upper, "max": max(lower, upper)}[side]

    def one_row(row: pd.Series) -> tuple[bool, str | None, float, float, float]:
        wt_reason = interval_reason(row["wt_gamma_mle"], row["wt_gamma_ci_low"], row["wt_gamma_ci_high"])
        epi_reason = interval_reason(row["epi_gamma_mle"], row["epi_gamma_ci_low"], row["epi_gamma_ci_high"])
        if wt_reason is not None:
            return False, f"wt_{wt_reason}", np.nan, np.nan, np.nan
        if epi_reason is not None:
            return False, f"epi_{epi_reason}", np.nan, np.nan, np.nan
        if pd.isna(row["delta_hat"]):
            return False, "undefined_delta", np.nan, np.nan, np.nan

        if mode == "null_facing":
            if abs(row["delta_hat"]) <= DELTA_TOL:
                wt_side = epi_side = "max"
            elif row["delta_hat"] > 0:
                wt_side, epi_side = "upper", "lower"
            else:
                wt_side, epi_side = "lower", "upper"
        elif mode == "max_half_width":
            wt_side = epi_side = "max"
        else:
            raise ValueError(f"Unknown CI mode '{mode}'.")

        se = math.sqrt(
            (half_width(row["wt_gamma_mle"], row["wt_gamma_ci_low"], row["wt_gamma_ci_high"], wt_side) / Z_975) ** 2
            + (half_width(row["epi_gamma_mle"], row["epi_gamma_ci_low"], row["epi_gamma_ci_high"], epi_side) / Z_975) ** 2
        )
        if not np.isfinite(se) or se <= 0:
            return False, "nonpositive_se", np.nan, np.nan, np.nan
        z = row["delta_hat"] / se
        p = math.erfc(abs(z) / math.sqrt(2.0))
        return True, None, se, z, p

    stats = contrasts.apply(one_row, axis=1, result_type="expand")
    stats.columns = [
        f"{prefix}_testable",
        f"{prefix}_exclusion_reason",
        f"{prefix}_se_delta",
        f"{prefix}_z_wald",
        f"{prefix}_p_raw",
    ]
    return pd.concat([contrasts, stats], axis=1)


def bh_adjust(p_values: pd.Series) -> pd.Series:
    """Purpose: apply Benjamini-Hochberg. Inputs: raw p-values. Outputs: q-values on the same index. Assumption: missing values were removed before calling."""

    if p_values.empty:
        return pd.Series(dtype=float, index=p_values.index)
    ordered = p_values.sort_values(kind="mergesort")
    ranks = np.arange(1, len(ordered) + 1, dtype=float)
    adjusted = np.minimum.accumulate((ordered.to_numpy(dtype=float) * len(ordered) / ranks)[::-1])[::-1]
    return pd.Series(np.clip(adjusted, 0.0, 1.0), index=ordered.index).reindex(p_values.index)
