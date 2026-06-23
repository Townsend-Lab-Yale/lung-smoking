"""
Shared modeling utilities for the restricted GENIE BPC LUAD survival workflow.

Supported workflow
------------------
This module supports only the retained analyses:
  1. the primary left-truncated Cox interaction analysis, and
  2. the metastatic-sample exclusion sensitivity analysis.

Functions here are intentionally limited to the data preparation and model
fitting needed by those analyses.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

BASE_COVARIATES = [
    "Age",
    "Sex_male",
    "Stage_advanced",
    "Smoking_ever",
    "Sample_metastatic",
]

PRIMARY_MIN_TOTAL = 20

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


def institution_dummies(df: pd.DataFrame) -> list[str]:
    """Return institution-indicator columns present in a model frame.

    Inputs
    ------
    df : analysis DataFrame

    Outputs
    -------
    List of column names beginning with `Inst_`.

    Assumptions
    -----------
    Institution indicators have already been created upstream in the survival
    table and use the `Inst_` prefix.
    """
    return [c for c in df.columns if c.startswith("Inst_")]


def build_primary_targets(feasibility: pd.DataFrame) -> pd.DataFrame:
    """Return the cohort-matched list of pair-level primary tests.

    Inputs
    ------
    feasibility : per-pair feasibility table from `01_prepare_genie_cohort.py`

    Outputs
    -------
    DataFrame with one row per (gene_A, gene_B, cohort) primary test, where the
    pair was epistasis-significant in that cohort and meets the cohort-specific
    `n_AB` threshold for survival modeling.

    Assumptions
    -----------
    `M2_cohort == "ES"` rows correspond to epistasis significance in
    ever-smokers and must satisfy `testable_ES`; likewise for the NS rows.
    """
    es = feasibility[
        (feasibility["M2_cohort"] == "ES") & feasibility["testable_ES"]
    ].copy()
    es["cohort"] = "ES"

    ns = feasibility[
        (feasibility["M2_cohort"] == "NS") & feasibility["testable_NS"]
    ].copy()
    ns["cohort"] = "NS"

    targets = pd.concat([es, ns], ignore_index=True)
    if len(targets) == 0:
        return pd.DataFrame(columns=["gene_A", "gene_B", "cohort"])

    return (
        targets[["gene_A", "gene_B", "cohort"]]
        .drop_duplicates()
        .sort_values(["cohort", "gene_A", "gene_B"])
        .reset_index(drop=True)
    )


def prepare_pair_data(
    surv: pd.DataFrame,
    coverage: pd.DataFrame,
    geno: pd.DataFrame,
    gene_a: str,
    gene_b: str,
    cohort: str = "pooled",
) -> pd.DataFrame | None:
    """Build the left-truncation-ready model frame for one gene pair.

    Inputs
    ------
    surv : patient-level survival table. Expected units are months for
        `OS_time_months` and `Entry_time_months`.
    coverage : sample x gene boolean matrix; True means the assay covers the
        gene in that sample.
    geno : sample x gene mutation matrix with 1 = mutated, 0 = wild-type,
        NaN = gene not covered by the assay.
    gene_a, gene_b : focal genes for the pairwise model.
    cohort : `pooled`, `ES`, or `NS`.

    Outputs
    -------
    DataFrame containing only dual-panel-covered samples in the requested
    cohort, with binary genotype columns `A`, `B`, `geno_AB`, `geno_A_only`,
    and `geno_B_only`, plus the covariates needed for the primary Cox model.

    Assumptions
    -----------
    Samples missing OS time, event indicator, left-truncation entry time, or
    retained covariates are excluded at the end of frame construction.
    """
    if gene_a not in geno.columns or gene_b not in geno.columns:
        return None

    cov_ab = coverage[gene_a] & coverage[gene_b]
    panel_samples = cov_ab[cov_ab].index
    surv_samples = pd.Index(surv["Sample_ID"])
    samples = panel_samples.intersection(surv_samples)
    if len(samples) == 0:
        return None

    a = geno.reindex(samples)[gene_a].astype(int)
    b = geno.reindex(samples)[gene_b].astype(int)

    df = surv.set_index("Sample_ID").reindex(samples).copy()
    df["A"] = a
    df["B"] = b
    df["geno_AB"] = ((a == 1) & (b == 1)).astype(int)
    df["geno_A_only"] = ((a == 1) & (b == 0)).astype(int)
    df["geno_B_only"] = ((a == 0) & (b == 1)).astype(int)

    if cohort == "ES":
        df = df[df["Smoking_ever"] == 1]
    elif cohort == "NS":
        df = df[df["Smoking_ever"] == 0]

    keep = [
        "OS_time_months",
        "OS_event",
        "Entry_time_months",
        "A",
        "B",
        "geno_AB",
        "geno_A_only",
        "geno_B_only",
    ] + BASE_COVARIATES + institution_dummies(df)
    keep = [c for c in keep if c in df.columns]
    df = df[keep].dropna()

    if len(df) == 0:
        return None
    return df


def _fit_cox(df: pd.DataFrame, covariates: list[str]) -> CoxPHFitter | None:
    """Fit a left-truncated Cox model on the provided covariate set.

    Inputs
    ------
    df : model frame with OS time, event indicator, and entry time in months
    covariates : covariate columns to include in the linear predictor

    Outputs
    -------
    Fitted `CoxPHFitter` or `None` if the fit fails.

    Assumptions
    -----------
    Zero-variance covariates are dropped before fitting because they make the
    Cox solver singular.
    """
    cols = ["OS_time_months", "OS_event", "Entry_time_months"] + covariates
    df_use = df[cols].copy()

    drop = [c for c in covariates if df_use[c].nunique() < 2]
    covariates = [c for c in covariates if c not in drop]
    df_use = df_use.drop(columns=drop)

    cph = CoxPHFitter(penalizer=0.0)
    try:
        cph.fit(
            df_use,
            duration_col="OS_time_months",
            event_col="OS_event",
            entry_col="Entry_time_months",
            show_progress=False,
        )
    except Exception:
        return None
    return cph


def _hazard_row(cph: CoxPHFitter, name: str) -> dict[str, float]:
    """Return hazard-ratio summary fields for one fitted coefficient."""
    if name not in cph.params_.index:
        return {
            f"HR_{name}": np.nan,
            f"HR_{name}_lo": np.nan,
            f"HR_{name}_hi": np.nan,
            f"p_{name}": np.nan,
        }

    s = cph.summary.loc[name]
    return {
        f"HR_{name}": float(np.exp(s["coef"])),
        f"HR_{name}_lo": float(np.exp(s["coef lower 95%"])),
        f"HR_{name}_hi": float(np.exp(s["coef upper 95%"])),
        f"p_{name}": float(s["p"]),
    }


def fit_interaction(df: pd.DataFrame, drop_smoking: bool = False) -> dict | None:
    """Fit the retained primary Cox interaction model.

    Inputs
    ------
    df : pair-specific survival frame from `prepare_pair_data`
    drop_smoking : whether to omit `Smoking_ever` from the covariates; this
        should be True when fitting within ES or NS where smoking is constant

    Outputs
    -------
    Dictionary containing sample counts, main-effect HRs for `A` and `B`, and
    the primary interaction-term HR for `A:B`.

    Assumptions
    -----------
    The covariates enter linearly on the log-hazard scale. The returned hazard
    ratios are exponentiated Cox coefficients.
    """
    model_df = df.copy()
    model_df["AB_int"] = model_df["A"] * model_df["B"]

    covs = ["A", "B", "AB_int"] + [
        c for c in BASE_COVARIATES if not (drop_smoking and c == "Smoking_ever")
    ] + institution_dummies(model_df)

    cph = _fit_cox(model_df, covs)
    if cph is None:
        return None

    out = {
        "n_total": len(model_df),
        "n_events": int(model_df["OS_event"].sum()),
        "n_AB": int(model_df["geno_AB"].sum()),
        "n_A_only": int(model_df["geno_A_only"].sum()),
        "n_B_only": int(model_df["geno_B_only"].sum()),
        "n_neither": int(
            (
                (
                    model_df["geno_AB"]
                    + model_df["geno_A_only"]
                    + model_df["geno_B_only"]
                )
                == 0
            ).sum()
        ),
        "concordance": float(cph.concordance_index_),
    }
    out.update(_hazard_row(cph, "A"))
    out.update(_hazard_row(cph, "B"))

    if "AB_int" not in cph.params_.index:
        out.update(
            {
                "HR_interaction": np.nan,
                "HR_int_lo": np.nan,
                "HR_int_hi": np.nan,
                "p_interaction": np.nan,
                "logHR_interaction": np.nan,
                "se_interaction": np.nan,
            }
        )
        return out

    s = cph.summary.loc["AB_int"]
    out.update(
        {
            "HR_interaction": float(np.exp(s["coef"])),
            "HR_int_lo": float(np.exp(s["coef lower 95%"])),
            "HR_int_hi": float(np.exp(s["coef upper 95%"])),
            "p_interaction": float(s["p"]),
            "logHR_interaction": float(s["coef"]),
            "se_interaction": float(s["se(coef)"]),
        }
    )
    return out
