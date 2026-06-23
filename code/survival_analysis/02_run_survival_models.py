"""
02_run_survival_models.py — restricted modeling workflow for the GENIE BPC
LUAD survival analysis.

This entrypoint combines:
    1. the primary left-truncated Cox interaction analysis, and
    2. the primary-tumor-only metastatic-exclusion sensitivity analysis.

Outputs (under code/survival_analysis/output/):
    cox_results_interaction.csv
    cox_results_primary_only_sensitivity.csv
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from survival_core import (
    PRIMARY_MIN_TOTAL,
    build_primary_targets,
    fit_interaction,
    prepare_pair_data,
)

OUTDIR = Path(__file__).resolve().parent / "output"

MIN_N_AB = 5
MIN_N_TOTAL = 30
MIN_EVENTS = 10


def run_primary_analysis() -> pd.DataFrame:
    """Fit the primary interaction Cox models for all cohort-matched targets.

    Inputs
    ------
    None. Reads `survival_table.csv`, `panel_coverage.csv`,
    `genotype_matrix.csv`, and `feasibility_table.csv` from `output/`.

    Outputs
    -------
    DataFrame of per-pair primary-model results, with BH-adjusted
    `q_interaction_fdr`.

    Assumptions
    -----------
    Only pairs that are both M2-significant and feasible in the matching
    smoking stratum are tested.
    """
    surv = pd.read_csv(OUTDIR / "survival_table.csv")
    coverage = pd.read_csv(OUTDIR / "panel_coverage.csv", index_col=0).astype(bool)
    geno = pd.read_csv(OUTDIR / "genotype_matrix.csv", index_col=0)
    feas = pd.read_csv(OUTDIR / "feasibility_table.csv")

    primary_targets = build_primary_targets(feas)
    print(f"  primary (pair, cohort) tests queued: {len(primary_targets)}")

    rows = []
    for _, row in primary_targets.iterrows():
        gene_a = row["gene_A"]
        gene_b = row["gene_B"]
        cohort = row["cohort"]

        df = prepare_pair_data(surv, coverage, geno, gene_a, gene_b, cohort=cohort)
        if df is None or len(df) < PRIMARY_MIN_TOTAL:
            continue

        fit = fit_interaction(df, drop_smoking=True)
        if fit is None:
            continue

        rows.append({"gene_A": gene_a, "gene_B": gene_b, "cohort": cohort, **fit})

    primary = pd.DataFrame(rows)
    if len(primary):
        pvals = primary["p_interaction"].values
        mask = np.isfinite(pvals)
        qvals = np.full_like(pvals, np.nan, dtype=float)
        if mask.sum() > 0:
            qvals[mask] = multipletests(pvals[mask], method="fdr_bh")[1]
        primary["q_interaction_fdr"] = qvals

    primary.to_csv(OUTDIR / "cox_results_interaction.csv", index=False)
    return primary


def fit_subset(
    surv_subset: pd.DataFrame,
    coverage: pd.DataFrame,
    geno: pd.DataFrame,
    gene_a: str,
    gene_b: str,
    cohort: str,
) -> dict:
    """Fit one subset-specific primary interaction model when feasible.

    Inputs
    ------
    surv_subset : restricted survival subset
    coverage : sample x gene boolean panel-coverage matrix
    geno : sample x gene mutation matrix
    gene_a, gene_b : focal genes
    cohort : `ES` or `NS`

    Outputs
    -------
    Dictionary of hazard-ratio estimates, subset sample counts, and a
    `skip_reason` when the subset is too sparse to fit.

    Assumptions
    -----------
    The subset analysis is descriptive; no additional multiplicity correction
    is applied within these refits.
    """
    out = {
        "HR": np.nan,
        "HR_lo": np.nan,
        "HR_hi": np.nan,
        "p": np.nan,
        "logHR": np.nan,
        "n_total": 0,
        "n_events": 0,
        "n_AB": 0,
        "skip_reason": "",
    }

    df = prepare_pair_data(surv_subset, coverage, geno, gene_a, gene_b, cohort)
    if df is None or len(df) == 0:
        out["skip_reason"] = "no samples after subset"
        return out

    n_ab = int(df["geno_AB"].sum())
    n_events = int(df["OS_event"].sum())
    n_total = len(df)
    out.update(n_total=n_total, n_events=n_events, n_AB=n_ab)

    if n_ab < MIN_N_AB:
        out["skip_reason"] = f"n_AB={n_ab} < {MIN_N_AB}"
        return out
    if n_total < MIN_N_TOTAL:
        out["skip_reason"] = f"n_total={n_total} < {MIN_N_TOTAL}"
        return out
    if n_events < MIN_EVENTS:
        out["skip_reason"] = f"n_events={n_events} < {MIN_EVENTS}"
        return out

    fit = fit_interaction(df, drop_smoking=True)
    if fit is None:
        out["skip_reason"] = "Cox fit failed"
        return out

    out.update(
        HR=fit.get("HR_interaction", np.nan),
        HR_lo=fit.get("HR_int_lo", np.nan),
        HR_hi=fit.get("HR_int_hi", np.nan),
        p=fit.get("p_interaction", np.nan),
        logHR=fit.get("logHR_interaction", np.nan),
    )
    return out


def sign_concord(log_full: float, log_sub: float) -> float:
    """Return 1 when the subset and full-cohort log-HR signs agree."""
    if not (np.isfinite(log_full) and np.isfinite(log_sub)):
        return np.nan
    if log_full == 0 or log_sub == 0:
        return np.nan
    return float(np.sign(log_full) == np.sign(log_sub))


def run_primary_only_sensitivity(primary: pd.DataFrame) -> pd.DataFrame:
    """Refit the primary model in primary-tumor-only sensitivity subsets.

    Inputs
    ------
    primary : full-cohort primary-analysis results

    Outputs
    -------
    DataFrame comparing full-cohort and subset-specific interaction estimates.

    Assumptions
    -----------
    `primary` must contain the pair list and full-cohort estimates written by
    `run_primary_analysis()`.
    """
    surv = pd.read_csv(OUTDIR / "survival_table.csv")
    coverage = pd.read_csv(OUTDIR / "panel_coverage.csv", index_col=0).astype(bool)
    geno = pd.read_csv(OUTDIR / "genotype_matrix.csv", index_col=0)

    subsets = {
        "primAll": surv[surv["Sample_metastatic"] == 0].copy(),
        "primLow": surv[
            (surv["Sample_metastatic"] == 0) & (surv["Stage_advanced"] == 0)
        ].copy(),
        "primHigh": surv[
            (surv["Sample_metastatic"] == 0) & (surv["Stage_advanced"] == 1)
        ].copy(),
    }
    print("\nSubset sizes (events / n):")
    for tag, subset in subsets.items():
        print(f"  {tag:8s}: {int(subset['OS_event'].sum())} / {len(subset)}")

    primary = primary.copy()
    primary["logHR_full"] = np.log(primary["HR_interaction"])

    rows = []
    for _, r in primary.iterrows():
        gene_a = r["gene_A"]
        gene_b = r["gene_B"]
        cohort = r["cohort"]

        row = {
            "gene_A": gene_a,
            "gene_B": gene_b,
            "cohort": cohort,
            "HR_full": r["HR_interaction"],
            "p_full": r["p_interaction"],
            "q_full": r["q_interaction_fdr"],
            "n_AB_full": int(r["n_AB"]),
            "n_total_full": int(r["n_total"]),
        }

        for tag, subset in subsets.items():
            res = fit_subset(subset, coverage, geno, gene_a, gene_b, cohort)
            row.update(
                {
                    f"HR_{tag}": res["HR"],
                    f"HR_{tag}_lo": res["HR_lo"],
                    f"HR_{tag}_hi": res["HR_hi"],
                    f"p_{tag}": res["p"],
                    f"n_total_{tag}": res["n_total"],
                    f"n_events_{tag}": res["n_events"],
                    f"n_AB_{tag}": res["n_AB"],
                    f"skip_{tag}": res["skip_reason"],
                    f"sign_concord_{tag}": sign_concord(
                        r["logHR_full"], res["logHR"]
                    ),
                }
            )
        rows.append(row)

    out = pd.DataFrame(rows)
    out.to_csv(OUTDIR / "cox_results_primary_only_sensitivity.csv", index=False)
    return out


def main() -> None:
    """Run the restricted survival-model workflow end to end."""
    print("=" * 72)
    print("GENIE BPC LUAD — primary and sensitivity survival models")
    print("=" * 72)

    print("\n[1/2] Fitting primary interaction models…")
    primary = run_primary_analysis()
    print(f"\nWrote {OUTDIR / 'cox_results_interaction.csv'}")

    print("\n[2/2] Fitting metastatic-exclusion sensitivity models…")
    sensitivity = run_primary_only_sensitivity(primary)
    print(f"\nWrote {OUTDIR / 'cox_results_primary_only_sensitivity.csv'}")

    if len(sensitivity):
        cols = [
            "gene_A",
            "gene_B",
            "cohort",
            "HR_full",
            "HR_primAll",
            "HR_primHigh",
            "HR_primLow",
            "n_AB_primAll",
            "n_AB_primHigh",
            "n_AB_primLow",
        ]
        print("\nTop rows by full-cohort q-value:")
        print(
            sensitivity.sort_values("q_full")[cols].to_string(
                index=False, float_format="%.3f"
            )
        )

        for tag in ["primAll", "primHigh", "primLow"]:
            sc = sensitivity[f"sign_concord_{tag}"]
            n_eval = sc.notna().sum()
            n_concord = int((sc == 1).sum())
            n_skip = int((sensitivity[f"skip_{tag}"] != "").sum())
            print(
                f"  direction concordance in {tag}: "
                f"{n_concord}/{n_eval} (skipped {n_skip})"
            )

    print("\nDone.")


if __name__ == "__main__":
    main()
