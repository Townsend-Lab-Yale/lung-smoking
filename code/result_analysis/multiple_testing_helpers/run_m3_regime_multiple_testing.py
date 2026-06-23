"""Task 7 M=3 regime-specific multiple-testing analysis.

This script is downstream of the task0 model outputs. It constructs the
higher-order epistasis test families requested for the reviewer response and
reuses the Task 7 CI-derived log-Wald machinery from the M=2 analysis.
"""

from __future__ import annotations

import argparse
import math
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from epistasis_testing_helpers import (
    DEFAULT_RUN_ROOT,
    PROJECT_ROOT,
    add_wald_columns,
    bh_adjust,
    combo_name,
    load_transitions,
    resolve_result_dir,
)


DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "output" / "task7_multiple_testing" / "m3_regime_correction"
MODEL_ORDER = 3
ALPHA = 0.05
RELAXED_ALPHA = 0.10
DELTA_TOL = 1e-12
REGIME_ORDER = ["N+E", "S+S", "S+A", "N+N"]
AA_EXCLUDED_REGIME = "A+A"
DVSS_WT_BH_GATED_REGIMES = {"N+E", "S+S"}
DVSS_UNGATED_REGIMES = {"S+A"}
DVSS_REGIME_ORDER = ["N+E", "S+S", "S+A"]


MANUSCRIPT_INTERACTIONS = {
    ("smoking_plus", "TP53_STK11_SMARCA4", "SMARCA4"): "text: SMARCA4 in STK11-TP53, ES-LUAD",
    ("nonsmoking_plus", "TP53_EGFR_RBM10", "RBM10"): "text/figure: RBM10 in TP53-EGFR, NS-LUAD",
    ("smoking_plus", "TP53_KRAS_ARID1A", "ARID1A"): "text/figure: ARID1A in KRAS-TP53, ES-LUAD",
    ("smoking_plus", "KRAS_STK11_ATM", "ATM"): "text: ATM in KRAS-STK11, ES-LUAD",
    ("smoking_plus", "KEAP1_STK11_SMARCA4", "SMARCA4"): "text: SMARCA4 in KEAP1-STK11, ES-LUAD",
    ("smoking_plus", "KRAS_KEAP1_SMARCA4", "SMARCA4"): "text: SMARCA4 in KRAS-KEAP1, ES-LUAD",
    ("smoking_plus", "TP53_KRAS_EGFR", "EGFR"): "text/figure: EGFR in TP53-KRAS, ES-LUAD",
    ("smoking_plus", "KEAP1_STK11_APC", "APC"): "text/figure: APC in KEAP1-STK11, ES-LUAD",
    ("nonsmoking_plus", "EGFR_PIK3CA_ARID1A", "ARID1A"): "text/figure: ARID1A in EGFR-PIK3CA, NS-LUAD",
    ("smoking_plus", "TP53_EGFR_STK11", "STK11"): "text: STK11 in TP53-EGFR, ES-LUAD",
    ("smoking_plus", "KRAS_SMARCA4_APC", "APC"): "figure: APC in KRAS-SMARCA4, ES-LUAD",
}


def parse_args() -> argparse.Namespace:
    """Purpose: parse CLI inputs. Inputs: command-line args. Outputs: run options. Assumptions: standard project layout."""

    parser = argparse.ArgumentParser(
        description="Run M=3 regime-specific BH correction for higher-order epistasis."
    )
    parser.add_argument("--run-root", default=str(DEFAULT_RUN_ROOT))
    parser.add_argument("--output-root", default=str(DEFAULT_OUTPUT_ROOT))
    parser.add_argument("--method", default="variant", choices=["variant", "cesR"])
    parser.add_argument("--alpha", type=float, default=ALPHA)
    return parser.parse_args()


def classify_single_effect(row: pd.Series) -> str:
    """Purpose: classify one single-mutant background. Inputs: single-vs-WT row. Outputs: N/S/A class. Assumptions: ratio > 1 is synergistic and ratio < 1 is antagonistic when CI non-overlap is present."""

    if not bool(row["single_legacy_signif"]):
        return "N"
    if row["single_ratio"] > 1.0:
        return "S"
    if row["single_ratio"] < 1.0:
        return "A"
    return "N"


def classify_regime(single_classes: list[str]) -> tuple[str, str]:
    """Purpose: assign reviewer-requested M=3 regime. Inputs: two single-mutant classes. Outputs: pooled regime and detailed regime. Assumptions: A+A is excluded from the primary manuscript-referenced analysis."""

    counts = {label: single_classes.count(label) for label in ["N", "S", "A"]}
    if counts["A"] == 2:
        return AA_EXCLUDED_REGIME, AA_EXCLUDED_REGIME
    if counts["N"] == 2:
        return "N+N", "N+N"
    if counts["S"] == 2:
        return "S+S", "S+S"
    if counts["S"] == 1 and counts["A"] == 1:
        return "S+A", "S+A"
    if counts["N"] == 1 and counts["S"] == 1:
        return "N+E", "N+S"
    if counts["N"] == 1 and counts["A"] == 1:
        return "N+E", "N+A"
    return "unclassified", "+".join(sorted(single_classes))


def ci_nonoverlap(row_a: pd.Series, row_b: pd.Series) -> bool:
    """Purpose: reproduce legacy CI non-overlap for two gamma estimates. Inputs: two transition rows. Outputs: boolean. Assumptions: intervals are on the original gamma scale."""

    return bool(
        (row_a["gamma_ci_low"] > row_b["gamma_ci_high"])
        or (row_b["gamma_ci_low"] > row_a["gamma_ci_high"])
    )


def ratio_and_delta(numerator: float, denominator: float) -> tuple[float, float]:
    """Purpose: compute interpretable ratio and log effect. Inputs: numerator and denominator estimates. Outputs: ratio and log ratio. Assumptions: nonpositive values are marked missing."""

    if numerator > 0 and denominator > 0:
        ratio = numerator / denominator
        return ratio, math.log(ratio)
    return np.nan, np.nan


def manuscript_key(key: str, gene_set: str, mutated_gene: str) -> tuple[str, str, str]:
    """Purpose: normalize manuscript subset lookup. Inputs: cohort, gene_set, target. Outputs: tuple key. Assumptions: gene_set is underscore-joined in model order."""

    return key, gene_set, mutated_gene


def build_m3_regime_rows(transitions: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Purpose: construct primary WT-reference tests and descriptive single-reference rows. Inputs: M=3 transitions. Outputs: primary, descriptive, excluded tables. Assumptions: each fitted triad-target has one WT, two single-background, and one double-background transition."""

    primary_records: list[dict[str, object]] = []
    descriptive_records: list[dict[str, object]] = []
    excluded_records: list[dict[str, object]] = []

    group_cols = ["key", "gene_set", "mutated_gene"]
    for (key, gene_set, mutated_gene), group in transitions.groupby(group_cols, sort=False):
        wt_rows = group.loc[group["context_size"] == 0]
        single_rows = group.loc[group["context_size"] == 1].copy()
        double_rows = group.loc[group["context_size"] == 2]

        base_record = {
            "method": group["method"].iloc[0],
            "model_order": MODEL_ORDER,
            "key": key,
            "gene_set": gene_set,
            "mutated_gene": mutated_gene,
        }
        if len(wt_rows) != 1 or len(single_rows) != 2 or len(double_rows) != 1:
            excluded_records.append(
                {
                    **base_record,
                    "regime": "incomplete_transition_set",
                    "regime_detail": "incomplete_transition_set",
                    "exclusion_reason": (
                        f"expected 1 WT, 2 single, 1 double; observed "
                        f"{len(wt_rows)} WT, {len(single_rows)} single, {len(double_rows)} double"
                    ),
                }
            )
            continue

        wt = wt_rows.iloc[0]
        double = double_rows.iloc[0]
        single_rows = single_rows.sort_values("from_gt").reset_index(drop=True)

        single_summaries = []
        for _, single in single_rows.iterrows():
            single_ratio, single_delta = ratio_and_delta(single["gamma_mle"], wt["gamma_mle"])
            single_legacy_signif = ci_nonoverlap(single, wt)
            summary = {
                "single_context": single["from_gt"],
                "single_gamma_mle": single["gamma_mle"],
                "single_gamma_ci_low": single["gamma_ci_low"],
                "single_gamma_ci_high": single["gamma_ci_high"],
                "single_from_count": single["from_count"],
                "single_to_count": single["to_count"],
                "single_ratio": single_ratio,
                "single_delta_hat": single_delta,
                "single_legacy_signif": single_legacy_signif,
            }
            summary["single_class"] = classify_single_effect(pd.Series(summary))
            single_summaries.append(summary)

        single_classes = [item["single_class"] for item in single_summaries]
        regime, regime_detail = classify_regime(single_classes)
        double_ratio, double_delta = ratio_and_delta(double["gamma_mle"], wt["gamma_mle"])
        manuscript_lookup = manuscript_key(key, gene_set, mutated_gene)
        manuscript_reference = MANUSCRIPT_INTERACTIONS.get(manuscript_lookup, "")

        common = {
            **base_record,
            "regime": regime,
            "regime_detail": regime_detail,
            "single_context_1": single_summaries[0]["single_context"],
            "single_class_1": single_summaries[0]["single_class"],
            "single_ratio_1": single_summaries[0]["single_ratio"],
            "single_legacy_signif_1": single_summaries[0]["single_legacy_signif"],
            "single_context_2": single_summaries[1]["single_context"],
            "single_class_2": single_summaries[1]["single_class"],
            "single_ratio_2": single_summaries[1]["single_ratio"],
            "single_legacy_signif_2": single_summaries[1]["single_legacy_signif"],
            "single_classes": "+".join(single_classes),
            "epistatic_gt": double["from_gt"],
            "tested_combo": f"{double['from_gt']}_{mutated_gene}",
            "combo_name": combo_name(mutated_gene, double["from_gt"]),
            "context_size": int(double["context_size"]),
            "manuscript_subset": bool(manuscript_reference),
            "manuscript_reference": manuscript_reference,
        }

        primary_record = {
            **common,
            "contrast_type": "double_mutant_vs_wt",
            "included_in_primary_bh": regime in REGIME_ORDER,
            "wt_gamma_mle": wt["gamma_mle"],
            "wt_gamma_ci_low": wt["gamma_ci_low"],
            "wt_gamma_ci_high": wt["gamma_ci_high"],
            "wt_from_count": wt["from_count"],
            "wt_to_count": wt["to_count"],
            "epi_gamma_mle": double["gamma_mle"],
            "epi_gamma_ci_low": double["gamma_ci_low"],
            "epi_gamma_ci_high": double["gamma_ci_high"],
            "epi_from_count": double["from_count"],
            "epi_to_count": double["to_count"],
            "ratio": double_ratio,
            "delta_hat": double_delta,
            "legacy_signif": ci_nonoverlap(double, wt),
        }

        if regime in REGIME_ORDER:
            primary_records.append(primary_record)
        else:
            excluded_records.append(
                {
                    **primary_record,
                    "exclusion_reason": "A+A excluded from reviewer-requested manuscript-referenced regimes"
                    if regime == AA_EXCLUDED_REGIME
                    else "unclassified regime",
                }
            )

        for single_summary in single_summaries:
            single = single_rows.loc[single_rows["from_gt"] == single_summary["single_context"]].iloc[0]
            descriptive_ratio, descriptive_delta = ratio_and_delta(double["gamma_mle"], single["gamma_mle"])
            descriptive_records.append(
                {
                    **common,
                    "contrast_type": "double_mutant_vs_single_descriptive",
                    "included_in_primary_bh": False,
                    "reference_single_context": single_summary["single_context"],
                    "reference_single_class": single_summary["single_class"],
                    "wt_gamma_mle": single["gamma_mle"],
                    "wt_gamma_ci_low": single["gamma_ci_low"],
                    "wt_gamma_ci_high": single["gamma_ci_high"],
                    "wt_from_count": single["from_count"],
                    "wt_to_count": single["to_count"],
                    "epi_gamma_mle": double["gamma_mle"],
                    "epi_gamma_ci_low": double["gamma_ci_low"],
                    "epi_gamma_ci_high": double["gamma_ci_high"],
                    "epi_from_count": double["from_count"],
                    "epi_to_count": double["to_count"],
                    "ratio": descriptive_ratio,
                    "delta_hat": descriptive_delta,
                    "legacy_signif": ci_nonoverlap(double, single),
                }
            )

    return (
        pd.DataFrame(primary_records),
        pd.DataFrame(descriptive_records),
        pd.DataFrame(excluded_records),
    )


def add_regime_bh_columns(tests: pd.DataFrame, alpha: float) -> pd.DataFrame:
    """Purpose: add regime-specific primary BH calls. Inputs: WT-reference tests and alpha. Outputs: finalized tests. Assumptions: primary family pools cohorts within each regime, with by-cohort q-values reported as sensitivity."""

    tests = tests.copy()
    tests = add_wald_columns(tests, "null_facing", "primary")
    tests = add_wald_columns(tests, "max_half_width", "maxhalf")
    tests["p_raw_signif"] = tests["primary_testable"].fillna(False) & (tests["primary_p_raw"] <= alpha)
    tests["q_bh_within_regime"] = np.nan
    tests["bh_signif_within_regime"] = False
    tests["q_bh_within_regime_key"] = np.nan
    tests["bh_signif_within_regime_key"] = False
    tests["q_bh_maxhalf_within_regime"] = np.nan
    tests["bh_signif_maxhalf_within_regime"] = False

    for regime, idx in tests.groupby("regime", sort=False).groups.items():
        testable_idx = tests.loc[idx].index[tests.loc[idx, "primary_testable"].fillna(False)]
        tests.loc[testable_idx, "q_bh_within_regime"] = bh_adjust(tests.loc[testable_idx, "primary_p_raw"])
        maxhalf_idx = tests.loc[idx].index[tests.loc[idx, "maxhalf_testable"].fillna(False)]
        tests.loc[maxhalf_idx, "q_bh_maxhalf_within_regime"] = bh_adjust(tests.loc[maxhalf_idx, "maxhalf_p_raw"])

    for (regime, key), idx in tests.groupby(["regime", "key"], sort=False).groups.items():
        testable_idx = tests.loc[idx].index[tests.loc[idx, "primary_testable"].fillna(False)]
        tests.loc[testable_idx, "q_bh_within_regime_key"] = bh_adjust(tests.loc[testable_idx, "primary_p_raw"])

    tests["bh_signif_within_regime"] = (
        tests["primary_testable"].fillna(False) & (tests["q_bh_within_regime"] <= alpha)
    )
    tests["bh_signif_within_regime_key"] = (
        tests["primary_testable"].fillna(False) & (tests["q_bh_within_regime_key"] <= alpha)
    )
    tests["bh_signif_maxhalf_within_regime"] = (
        tests["maxhalf_testable"].fillna(False) & (tests["q_bh_maxhalf_within_regime"] <= alpha)
    )
    tests["direction"] = np.where(
        tests["delta_hat"] > DELTA_TOL,
        "synergistic",
        np.where(tests["delta_hat"] < -DELTA_TOL, "antagonistic", "neutral"),
    )
    tests["legacy_vs_bh_within_regime"] = np.where(
        ~tests["primary_testable"].fillna(False),
        "untestable",
        np.where(
            tests["legacy_signif"] & tests["bh_signif_within_regime"],
            "both",
            np.where(
                tests["legacy_signif"] & ~tests["bh_signif_within_regime"],
                "legacy_only",
                np.where(
                    ~tests["legacy_signif"] & tests["bh_signif_within_regime"],
                    "bh_only",
                    "neither",
                ),
            ),
        ),
    )
    tests["abs_log_ratio"] = tests["delta_hat"].abs()
    tests["regime"] = pd.Categorical(tests["regime"], categories=REGIME_ORDER, ordered=True)
    return tests.sort_values(
        ["regime", "key", "manuscript_subset", "bh_signif_within_regime", "primary_p_raw"],
        ascending=[True, True, False, False, True],
        na_position="last",
    ).reset_index(drop=True)


def add_descriptive_wald_columns(descriptive: pd.DataFrame) -> pd.DataFrame:
    """Purpose: compute unadjusted descriptive double-vs-single Wald columns. Inputs: descriptive contrasts. Outputs: annotated table. Assumptions: these rows are not part of the primary BH family."""

    if descriptive.empty:
        return descriptive
    descriptive = add_wald_columns(descriptive, "null_facing", "descriptive")
    descriptive = add_wald_columns(descriptive, "max_half_width", "descriptive_maxhalf")
    descriptive["direction"] = np.where(
        descriptive["delta_hat"] > DELTA_TOL,
        "synergistic",
        np.where(descriptive["delta_hat"] < -DELTA_TOL, "antagonistic", "neutral"),
    )
    return descriptive


def add_dvss_bh_followup_columns(
    tests: pd.DataFrame,
    descriptive: pd.DataFrame,
    alpha: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Purpose: add conditional double-vs-single BH follow-up families. Inputs: WT-reference tests and double-vs-single contrasts. Outputs: annotated primary and descriptive tables. Assumptions: N+E/S+S follow-up families are gated on WT-reference BH significance; S+A includes all double-vs-single contrasts; N+N remains WT-reference-only."""

    tests = tests.copy()
    descriptive = descriptive.copy()
    key_cols = ["key", "gene_set", "mutated_gene"]

    if descriptive.empty:
        descriptive["included_in_dvss_bh_family"] = False
        descriptive["dvss_bh_family"] = ""
        descriptive["dvss_p_raw"] = np.nan
        descriptive["dvss_q_bh_within_regime_family"] = np.nan
        descriptive["dvss_bh_signif"] = False
        descriptive["dvss_bh_signif_alpha_0_1"] = False
        descriptive["dvss_bh_and_ci_same_comparison"] = False
        descriptive["dvss_bh_and_ci_same_comparison_alpha_0_1"] = False
    else:
        descriptive["included_in_dvss_bh_family"] = False
        descriptive["dvss_bh_family"] = ""
        descriptive["dvss_p_raw"] = descriptive["descriptive_p_raw"]
        descriptive["dvss_q_bh_within_regime_family"] = np.nan
        descriptive["dvss_bh_signif"] = False
        descriptive["dvss_bh_signif_alpha_0_1"] = False

        gated_keys = set(
            tests.loc[
                tests["regime"].astype(str).isin(DVSS_WT_BH_GATED_REGIMES)
                & tests["bh_signif_within_regime"].fillna(False),
                key_cols,
            ].itertuples(index=False, name=None)
        )
        descriptive_keys = list(descriptive[key_cols].itertuples(index=False, name=None))
        gated_mask = (
            descriptive["regime"].astype(str).isin(DVSS_WT_BH_GATED_REGIMES)
            & pd.Series([item in gated_keys for item in descriptive_keys], index=descriptive.index)
        )
        ungated_mask = descriptive["regime"].astype(str).isin(DVSS_UNGATED_REGIMES)
        descriptive.loc[gated_mask | ungated_mask, "included_in_dvss_bh_family"] = True
        descriptive.loc[gated_mask, "dvss_bh_family"] = (
            descriptive.loc[gated_mask, "regime"].astype(str) + "_wt_bh_gated"
        )
        descriptive.loc[ungated_mask, "dvss_bh_family"] = (
            descriptive.loc[ungated_mask, "regime"].astype(str) + "_all"
        )

        for regime in DVSS_REGIME_ORDER:
            idx = descriptive.index[
                descriptive["included_in_dvss_bh_family"]
                & (descriptive["regime"].astype(str) == regime)
                & descriptive["descriptive_testable"].fillna(False)
            ]
            if len(idx) > 0:
                descriptive.loc[idx, "dvss_q_bh_within_regime_family"] = bh_adjust(
                    descriptive.loc[idx, "dvss_p_raw"]
                )

        descriptive["dvss_bh_signif"] = (
            descriptive["included_in_dvss_bh_family"]
            & descriptive["descriptive_testable"].fillna(False)
            & (descriptive["dvss_q_bh_within_regime_family"] <= alpha)
        )
        descriptive["dvss_bh_signif_alpha_0_1"] = (
            descriptive["included_in_dvss_bh_family"]
            & descriptive["descriptive_testable"].fillna(False)
            & (descriptive["dvss_q_bh_within_regime_family"] <= RELAXED_ALPHA)
        )
        descriptive["dvss_bh_and_ci_same_comparison"] = (
            descriptive["dvss_bh_signif"] & descriptive["legacy_signif"].fillna(False)
        )
        descriptive["dvss_bh_and_ci_same_comparison_alpha_0_1"] = (
            descriptive["dvss_bh_signif_alpha_0_1"] & descriptive["legacy_signif"].fillna(False)
        )

    all_single_summary = pd.DataFrame()
    if descriptive.empty:
        all_single_summary = pd.DataFrame(columns=key_cols)
    else:
        all_single_summary = (
            descriptive.groupby(key_cols, dropna=False)
            .agg(
                n_single_reference_tests=("legacy_signif", "size"),
                n_single_ci_nonoverlap=("legacy_signif", "sum"),
                any_single_ci_nonoverlap=("legacy_signif", "any"),
                min_single_reference_p_raw=("descriptive_p_raw", "min"),
            )
            .reset_index()
        )

    tests = tests.merge(all_single_summary, on=key_cols, how="left")
    tests["n_single_reference_tests"] = tests["n_single_reference_tests"].fillna(0).astype(int)
    tests["n_single_ci_nonoverlap"] = tests["n_single_ci_nonoverlap"].fillna(0).astype(int)
    tests["any_single_ci_nonoverlap"] = tests["any_single_ci_nonoverlap"].fillna(False).astype(bool)

    if descriptive.empty:
        family_summary = pd.DataFrame(
            columns=key_cols
            + [
                "dvss_bh_test_family",
                "n_dvss_tests_in_family",
                "n_dvss_bh_signif",
                "any_dvss_bh_signif",
                "n_dvss_bh_signif_alpha_0_1",
                "any_dvss_bh_signif_alpha_0_1",
                "n_dvss_ci_nonoverlap",
                "any_dvss_ci_nonoverlap",
                "n_dvss_bh_and_ci_same_comparison",
                "any_dvss_bh_and_ci_same_comparison",
                "n_dvss_bh_and_ci_same_comparison_alpha_0_1",
                "any_dvss_bh_and_ci_same_comparison_alpha_0_1",
                "min_dvss_p_raw",
                "min_dvss_q_bh",
            ]
        )
    else:
        family_summary = (
            descriptive.loc[descriptive["included_in_dvss_bh_family"]]
            .groupby(key_cols, dropna=False)
            .agg(
                dvss_bh_test_family=("dvss_bh_family", "first"),
                n_dvss_tests_in_family=("dvss_p_raw", "size"),
                n_dvss_bh_signif=("dvss_bh_signif", "sum"),
                any_dvss_bh_signif=("dvss_bh_signif", "any"),
                n_dvss_bh_signif_alpha_0_1=("dvss_bh_signif_alpha_0_1", "sum"),
                any_dvss_bh_signif_alpha_0_1=("dvss_bh_signif_alpha_0_1", "any"),
                n_dvss_ci_nonoverlap=("legacy_signif", "sum"),
                any_dvss_ci_nonoverlap=("legacy_signif", "any"),
                n_dvss_bh_and_ci_same_comparison=("dvss_bh_and_ci_same_comparison", "sum"),
                any_dvss_bh_and_ci_same_comparison=("dvss_bh_and_ci_same_comparison", "any"),
                n_dvss_bh_and_ci_same_comparison_alpha_0_1=(
                    "dvss_bh_and_ci_same_comparison_alpha_0_1",
                    "sum",
                ),
                any_dvss_bh_and_ci_same_comparison_alpha_0_1=(
                    "dvss_bh_and_ci_same_comparison_alpha_0_1",
                    "any",
                ),
                min_dvss_p_raw=("dvss_p_raw", "min"),
                min_dvss_q_bh=("dvss_q_bh_within_regime_family", "min"),
            )
            .reset_index()
        )

    tests = tests.merge(family_summary, on=key_cols, how="left")
    fill_int_cols = [
        "n_dvss_tests_in_family",
        "n_dvss_bh_signif",
        "n_dvss_bh_signif_alpha_0_1",
        "n_dvss_ci_nonoverlap",
        "n_dvss_bh_and_ci_same_comparison",
        "n_dvss_bh_and_ci_same_comparison_alpha_0_1",
    ]
    for col in fill_int_cols:
        tests[col] = tests[col].fillna(0).astype(int)
    fill_bool_cols = [
        "any_dvss_bh_signif",
        "any_dvss_bh_signif_alpha_0_1",
        "any_dvss_ci_nonoverlap",
        "any_dvss_bh_and_ci_same_comparison",
        "any_dvss_bh_and_ci_same_comparison_alpha_0_1",
    ]
    for col in fill_bool_cols:
        tests[col] = tests[col].map(lambda value: bool(value) if pd.notna(value) else False)
    tests["dvss_bh_test_family"] = tests["dvss_bh_test_family"].fillna("not_in_dvss_family")

    tests["wt_and_dvss_bh_signif"] = (
        tests["regime"].astype(str).isin(DVSS_WT_BH_GATED_REGIMES)
        & tests["bh_signif_within_regime"].fillna(False)
        & tests["any_dvss_bh_signif"]
    )
    tests["wt_and_dvss_bh_signif_alpha_0_1"] = (
        tests["regime"].astype(str).isin(DVSS_WT_BH_GATED_REGIMES)
        & tests["primary_testable"].fillna(False)
        & (tests["q_bh_within_regime"] <= RELAXED_ALPHA)
        & tests["any_dvss_bh_signif_alpha_0_1"]
    )
    tests["wt_and_dvss_ci_nonoverlap"] = (
        tests["regime"].astype(str).isin(DVSS_WT_BH_GATED_REGIMES)
        & tests["legacy_signif"].fillna(False)
        & tests["any_single_ci_nonoverlap"]
    )
    tests["reported_signif_revised"] = np.where(
        tests["regime"].astype(str).isin(DVSS_REGIME_ORDER),
        tests["any_dvss_bh_signif"],
        np.where(
            tests["regime"].astype(str) == "N+N",
            tests["bh_signif_within_regime"],
            False,
        ),
    ).astype(bool)
    tests["reported_signif_revised_alpha_0_1"] = np.where(
        tests["regime"].astype(str).isin(DVSS_REGIME_ORDER),
        tests["any_dvss_bh_signif_alpha_0_1"],
        np.where(
            tests["regime"].astype(str) == "N+N",
            tests["primary_testable"].fillna(False) & (tests["q_bh_within_regime"] <= RELAXED_ALPHA),
            False,
        ),
    ).astype(bool)
    tests["reported_signif_rule"] = np.where(
        tests["regime"].astype(str).isin(DVSS_REGIME_ORDER),
        "double_vs_single_bh",
        np.where(tests["regime"].astype(str) == "N+N", "wt_bh", "excluded"),
    )
    return tests, descriptive


def summarize_primary_tests(tests: pd.DataFrame, excluded: pd.DataFrame, alpha: float) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Purpose: summarize M=3 regime tests. Inputs: finalized tests, exclusions, alpha. Outputs: regime, regime/cohort, and exclusion summaries. Assumptions: A+A exclusions are counted outside the primary BH families."""

    regime_rows = []
    for regime in REGIME_ORDER:
        df = tests.loc[tests["regime"].astype(str) == regime]
        regime_rows.append(
            {
                "regime": regime,
                "n_total_tests": int(len(df)),
                "n_testable": int(df["primary_testable"].sum()),
                "n_untestable": int((~df["primary_testable"]).sum()),
                "n_legacy_signif": int(df["legacy_signif"].sum()),
                "n_ci_nonoverlap_signif": int(df["legacy_signif"].sum()),
                "n_raw_p_lt_alpha": int(df["p_raw_signif"].sum()),
                "n_bh_signif_within_regime": int(df["bh_signif_within_regime"].sum()),
                "n_wt_bh_signif_alpha_0_1": int(
                    (
                        df["primary_testable"].fillna(False)
                        & (df["q_bh_within_regime"] <= RELAXED_ALPHA)
                    ).sum()
                ),
                "n_dvss_bh_family_comparisons": int(df["n_dvss_tests_in_family"].sum()),
                "n_dvss_bh_family_triad_targets": int((df["n_dvss_tests_in_family"] > 0).sum()),
                "n_dvss_ci_nonoverlap_comparisons": int(df["n_dvss_ci_nonoverlap"].sum()),
                "n_dvss_ci_nonoverlap_triad_targets": int(df["any_dvss_ci_nonoverlap"].sum()),
                "n_dvss_bh_signif_comparisons": int(df["n_dvss_bh_signif"].sum()),
                "n_dvss_bh_signif_triad_targets": int(df["any_dvss_bh_signif"].sum()),
                "n_dvss_bh_signif_comparisons_alpha_0_1": int(
                    df["n_dvss_bh_signif_alpha_0_1"].sum()
                ),
                "n_dvss_bh_signif_triad_targets_alpha_0_1": int(
                    df["any_dvss_bh_signif_alpha_0_1"].sum()
                ),
                "n_dvss_bh_ci_overlap_comparisons": int(df["n_dvss_bh_and_ci_same_comparison"].sum()),
                "n_dvss_bh_ci_overlap_triad_targets": int(
                    (df["any_dvss_bh_signif"] & df["any_dvss_ci_nonoverlap"]).sum()
                ),
                "n_dvss_bh_ci_overlap_comparisons_alpha_0_1": int(
                    df["n_dvss_bh_and_ci_same_comparison_alpha_0_1"].sum()
                ),
                "n_dvss_bh_ci_overlap_triad_targets_alpha_0_1": int(
                    (df["any_dvss_bh_signif_alpha_0_1"] & df["any_dvss_ci_nonoverlap"]).sum()
                ),
                "n_dvss_bh_ci_same_comparison_triad_targets": int(
                    df["any_dvss_bh_and_ci_same_comparison"].sum()
                ),
                "n_wt_and_dvss_bh_signif": int(df["wt_and_dvss_bh_signif"].sum()),
                "n_wt_and_dvss_bh_signif_alpha_0_1": int(
                    df["wt_and_dvss_bh_signif_alpha_0_1"].sum()
                ),
                "n_wt_and_dvss_ci_nonoverlap": int(df["wt_and_dvss_ci_nonoverlap"].sum()),
                "n_reported_signif_revised": int(df["reported_signif_revised"].sum()),
                "n_reported_signif_revised_alpha_0_1": int(
                    df["reported_signif_revised_alpha_0_1"].sum()
                ),
                "n_bh_signif_maxhalf_within_regime": int(df["bh_signif_maxhalf_within_regime"].sum()),
                "n_both": int((df["legacy_vs_bh_within_regime"] == "both").sum()),
                "n_legacy_only": int((df["legacy_vs_bh_within_regime"] == "legacy_only").sum()),
                "n_bh_only": int((df["legacy_vs_bh_within_regime"] == "bh_only").sum()),
                "n_neither": int((df["legacy_vs_bh_within_regime"] == "neither").sum()),
                "n_synergistic_bh": int((df["bh_signif_within_regime"] & (df["direction"] == "synergistic")).sum()),
                "n_antagonistic_bh": int((df["bh_signif_within_regime"] & (df["direction"] == "antagonistic")).sum()),
                "n_manuscript_subset": int(df["manuscript_subset"].sum()),
                "bh_alpha": alpha,
            }
        )

    cohort_rows = []
    for (regime, key), df in tests.groupby(["regime", "key"], observed=True, sort=True):
        cohort_rows.append(
            {
                "regime": str(regime),
                "key": key,
                "n_total_tests": int(len(df)),
                "n_testable": int(df["primary_testable"].sum()),
                "n_legacy_signif": int(df["legacy_signif"].sum()),
                "n_ci_nonoverlap_signif": int(df["legacy_signif"].sum()),
                "n_bh_signif_within_regime_primary": int(df["bh_signif_within_regime"].sum()),
                "n_wt_bh_signif_alpha_0_1": int(
                    (
                        df["primary_testable"].fillna(False)
                        & (df["q_bh_within_regime"] <= RELAXED_ALPHA)
                    ).sum()
                ),
                "n_dvss_bh_family_comparisons": int(df["n_dvss_tests_in_family"].sum()),
                "n_dvss_bh_signif_triad_targets": int(df["any_dvss_bh_signif"].sum()),
                "n_dvss_bh_signif_triad_targets_alpha_0_1": int(
                    df["any_dvss_bh_signif_alpha_0_1"].sum()
                ),
                "n_dvss_ci_nonoverlap_triad_targets": int(df["any_dvss_ci_nonoverlap"].sum()),
                "n_reported_signif_revised": int(df["reported_signif_revised"].sum()),
                "n_reported_signif_revised_alpha_0_1": int(
                    df["reported_signif_revised_alpha_0_1"].sum()
                ),
                "n_bh_signif_within_regime_key_sensitivity": int(df["bh_signif_within_regime_key"].sum()),
                "n_manuscript_subset": int(df["manuscript_subset"].sum()),
            }
        )

    if excluded.empty:
        excluded_summary = pd.DataFrame()
    else:
        excluded_summary = (
            excluded.groupby(["regime", "regime_detail"], dropna=False)
            .size()
            .reset_index(name="n_excluded_triad_targets")
        )
    return pd.DataFrame(regime_rows), pd.DataFrame(cohort_rows), excluded_summary


def make_bh_results_table(regime_summary: pd.DataFrame) -> pd.DataFrame:
    """Purpose: create a compact reviewer-facing BH table. Inputs: regime summary. Outputs: selected counts by regime. Assumptions: N+E/S+S use WT-BH-gated double-vs-single BH follow-up; S+A uses ungated double-vs-single BH; N+N uses WT-reference BH."""

    columns = [
        "regime",
        "n_total_tests",
        "n_testable",
        "n_ci_nonoverlap_signif",
        "n_bh_signif_within_regime",
        "n_wt_bh_signif_alpha_0_1",
        "n_dvss_bh_family_comparisons",
        "n_dvss_bh_family_triad_targets",
        "n_dvss_ci_nonoverlap_comparisons",
        "n_dvss_ci_nonoverlap_triad_targets",
        "n_dvss_bh_signif_comparisons",
        "n_dvss_bh_signif_triad_targets",
        "n_dvss_bh_signif_comparisons_alpha_0_1",
        "n_dvss_bh_signif_triad_targets_alpha_0_1",
        "n_dvss_bh_ci_overlap_comparisons",
        "n_dvss_bh_ci_overlap_triad_targets",
        "n_dvss_bh_ci_overlap_comparisons_alpha_0_1",
        "n_dvss_bh_ci_overlap_triad_targets_alpha_0_1",
        "n_wt_and_dvss_bh_signif",
        "n_wt_and_dvss_bh_signif_alpha_0_1",
        "n_wt_and_dvss_ci_nonoverlap",
        "n_reported_signif_revised",
        "n_reported_signif_revised_alpha_0_1",
        "n_both",
        "n_legacy_only",
        "n_bh_only",
        "n_neither",
        "bh_alpha",
    ]
    table = regime_summary.loc[:, columns].copy()
    table = table.rename(
        columns={
            "n_legacy_only": "n_ci_nonoverlap_only",
            "n_bh_signif_within_regime": "n_wt_bh_signif",
        }
    )
    return table


def summarize_descriptive_tests(descriptive: pd.DataFrame, alpha: float) -> pd.DataFrame:
    """Purpose: summarize double-vs-single comparisons. Inputs: descriptive contrast table and alpha. Outputs: regime-level counts. Assumptions: each triad-target contributes two single-reference rows; BH families are flagged row-wise before summarization."""

    if descriptive.empty:
        return pd.DataFrame()

    rows = []
    for regime, df in descriptive.groupby("regime", sort=True):
        family_df = df.loc[df["included_in_dvss_bh_family"].fillna(False)]
        bh_targets = family_df.loc[family_df["dvss_bh_signif"].fillna(False), ["key", "gene_set", "mutated_gene"]]
        bh_targets_alpha_0_1 = family_df.loc[
            family_df["dvss_bh_signif_alpha_0_1"].fillna(False),
            ["key", "gene_set", "mutated_gene"],
        ]
        ci_targets = family_df.loc[family_df["legacy_signif"].fillna(False), ["key", "gene_set", "mutated_gene"]]
        bh_ci_targets = family_df.loc[
            family_df["dvss_bh_signif"].fillna(False) & family_df["legacy_signif"].fillna(False),
            ["key", "gene_set", "mutated_gene"],
        ]
        bh_ci_targets_alpha_0_1 = family_df.loc[
            family_df["dvss_bh_signif_alpha_0_1"].fillna(False)
            & family_df["legacy_signif"].fillna(False),
            ["key", "gene_set", "mutated_gene"],
        ]
        rows.append(
            {
                "regime": regime,
                "n_descriptive_single_reference_comparisons": int(len(df)),
                "n_distinct_triad_targets": int(
                    df[["key", "gene_set", "mutated_gene"]].drop_duplicates().shape[0]
                ),
                "n_legacy_nonoverlap_vs_single": int(df["legacy_signif"].sum()),
                "n_raw_p_lt_alpha_vs_single": int(
                    (df["descriptive_testable"].fillna(False) & (df["descriptive_p_raw"] <= alpha)).sum()
                ),
                "n_dvss_bh_family_comparisons": int(len(family_df)),
                "n_dvss_bh_family_triad_targets": int(
                    family_df[["key", "gene_set", "mutated_gene"]].drop_duplicates().shape[0]
                ),
                "n_dvss_bh_signif_comparisons": int(family_df["dvss_bh_signif"].sum()),
                "n_dvss_bh_signif_triad_targets": int(bh_targets.drop_duplicates().shape[0]),
                "n_dvss_bh_signif_comparisons_alpha_0_1": int(
                    family_df["dvss_bh_signif_alpha_0_1"].sum()
                ),
                "n_dvss_bh_signif_triad_targets_alpha_0_1": int(
                    bh_targets_alpha_0_1.drop_duplicates().shape[0]
                ),
                "n_dvss_ci_nonoverlap_comparisons_in_family": int(family_df["legacy_signif"].sum()),
                "n_dvss_ci_nonoverlap_triad_targets_in_family": int(ci_targets.drop_duplicates().shape[0]),
                "n_dvss_bh_ci_overlap_comparisons": int(
                    (family_df["dvss_bh_signif"].fillna(False) & family_df["legacy_signif"].fillna(False)).sum()
                ),
                "n_dvss_bh_ci_overlap_triad_targets_same_comparison": int(
                    bh_ci_targets.drop_duplicates().shape[0]
                ),
                "n_dvss_bh_ci_overlap_comparisons_alpha_0_1": int(
                    (
                        family_df["dvss_bh_signif_alpha_0_1"].fillna(False)
                        & family_df["legacy_signif"].fillna(False)
                    ).sum()
                ),
                "n_dvss_bh_ci_overlap_triad_targets_same_comparison_alpha_0_1": int(
                    bh_ci_targets_alpha_0_1.drop_duplicates().shape[0]
                ),
                "n_dvss_bh_ci_overlap_triad_targets_any_comparison": int(
                    len(
                        set(bh_targets.itertuples(index=False, name=None))
                        & set(ci_targets.itertuples(index=False, name=None))
                    )
                ),
                "n_dvss_bh_ci_overlap_triad_targets_any_comparison_alpha_0_1": int(
                    len(
                        set(bh_targets_alpha_0_1.itertuples(index=False, name=None))
                        & set(ci_targets.itertuples(index=False, name=None))
                    )
                ),
                "n_manuscript_subset_comparisons": int(df["manuscript_subset"].sum()),
                "alpha": alpha,
            }
        )
    return pd.DataFrame(rows)


def make_interaction_subset_tables(primary_tests: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Purpose: select reviewer-facing interaction subsets. Inputs: annotated primary tests. Outputs: revised reported calls, WT+double-vs-single BH calls, and equivalent CI-overlap calls. Assumptions: WT+double-vs-single subsets are only defined for N+E and S+S."""

    subset_cols = [
        "key",
        "gene_set",
        "mutated_gene",
        "regime",
        "regime_detail",
        "epistatic_gt",
        "ratio",
        "primary_p_raw",
        "q_bh_within_regime",
        "bh_signif_within_regime",
        "reported_signif_revised_alpha_0_1",
        "legacy_signif",
        "n_single_reference_tests",
        "n_single_ci_nonoverlap",
        "any_single_ci_nonoverlap",
        "dvss_bh_test_family",
        "n_dvss_tests_in_family",
        "n_dvss_bh_signif",
        "any_dvss_bh_signif",
        "n_dvss_bh_signif_alpha_0_1",
        "any_dvss_bh_signif_alpha_0_1",
        "n_dvss_ci_nonoverlap",
        "any_dvss_ci_nonoverlap",
        "n_dvss_bh_and_ci_same_comparison",
        "any_dvss_bh_and_ci_same_comparison",
        "min_dvss_p_raw",
        "min_dvss_q_bh",
        "reported_signif_rule",
        "reported_signif_revised",
        "wt_and_dvss_bh_signif",
        "wt_and_dvss_bh_signif_alpha_0_1",
        "wt_and_dvss_ci_nonoverlap",
        "manuscript_subset",
        "manuscript_reference",
    ]
    revised_reported = primary_tests.loc[primary_tests["reported_signif_revised"], subset_cols].copy()
    wt_dvss_bh = primary_tests.loc[primary_tests["wt_and_dvss_bh_signif"], subset_cols].copy()
    wt_dvss_ci = primary_tests.loc[primary_tests["wt_and_dvss_ci_nonoverlap"], subset_cols].copy()
    return revised_reported, wt_dvss_bh, wt_dvss_ci


def write_manifest(
    output_root: Path,
    run_root: Path,
    result_dir: Path,
    method: str,
    primary_tests: pd.DataFrame,
    descriptive_tests: pd.DataFrame,
    excluded: pd.DataFrame,
) -> None:
    """Purpose: record M=3 analysis provenance. Inputs: paths and result tables. Outputs: markdown manifest. Assumptions: concise methods summary supports reviewer audit."""

    lines = [
        "# M=3 Regime Multiple-Testing Manifest",
        "",
        "## Inputs",
        f"- Source run root: `{run_root}`",
        f"- Source result dir: `{result_dir}`",
        f"- Method: `{method}`",
        "",
        "## Primary Test Definition",
        "- Inferential unit: one fitted triad-target double-mutant transition.",
        "- Primary contrast: selection for the target mutation in the double-mutant background versus WT, within the same M=3 model.",
        "- Regime classification: based on the two single-mutant-vs-WT effects for the same triad-target.",
        "- Included regimes: `N+E` pooled from `N+S` and `N+A`, `S+S`, `S+A`, and `N+N`.",
        "- Excluded regime: `A+A`, because it was not referenced in the manuscript subset requested for this response.",
        "- Double-mutant-vs-single-mutant contrasts are tested as a conditional follow-up family for selected regimes.",
        "",
        "## Multiple Testing",
        "- Wald p-values use the same CI-derived log-scale null-facing standard error as the M=2 Task 7 analysis.",
        "- WT-reference BH correction is applied within each regime, pooling `smoking_plus` and `nonsmoking_plus` as in the M=2 Task 7 correction family.",
        "- A by-cohort within-regime q-value is also reported as a sensitivity column, not as the primary call.",
        "- Double-mutant-vs-single-mutant BH correction is comparison-level within each regime-specific follow-up family.",
        "- `N+E` and `S+S` double-mutant-vs-single-mutant families are gated on WT-reference BH significance.",
        "- `S+A` double-mutant-vs-single-mutant tests include all comparisons in that regime.",
        "- `N+N` remains WT-reference-only.",
        "- CI non-overlap overlap with double-mutant-vs-single-mutant BH calls is reported separately, not used as the revised BH call criterion.",
        "- Counts using the same q-values at relaxed `alpha = 0.1` are reported separately from the primary `alpha = 0.05` calls.",
        "",
        "## Software",
        f"- Python executable: `{sys.executable}`",
        f"- Python version: `{platform.python_version()}`",
        f"- pandas version: `{pd.__version__}`",
        f"- numpy version: `{np.__version__}`",
        "",
        "## Output Counts",
        f"- Primary tests: `{len(primary_tests)}`",
        f"- Descriptive double-vs-single rows: `{len(descriptive_tests)}`",
        f"- Excluded triad-target rows: `{len(excluded)}`",
        "",
        "## Main Tables",
        "- `M3_regime_bh_results_table.csv`: compact regime-level WT and double-mutant-vs-single BH counts with CI overlap counts.",
        "- `M3_regime_primary_wt_tests.csv`: row-level WT-reference tests, revised reported calls, and triad-target follow-up summaries.",
        "- `M3_regime_descriptive_single_reference_tests.csv`: comparison-level double-mutant-vs-single tests and follow-up q-values.",
        "- `M3_regime_revised_reported_significant_interactions.csv`: final revised reported-significant triad-targets.",
    ]
    (output_root / "M3_regime_run_manifest.md").write_text("\n".join(lines) + "\n")


def run_m3_regime_analysis(
    result_dir: Path,
    output_root: Path,
    run_root: Path,
    method: str,
    alpha: float,
) -> dict[str, pd.DataFrame]:
    """Purpose: run the M=3 regime workflow as a reusable pipeline step. Inputs: resolved result/output paths, method, alpha. Outputs: named result tables written to disk and returned in memory. Assumptions: M=3 task0 CSV exports are present in result_dir."""

    output_root.mkdir(parents=True, exist_ok=True)

    transitions = load_transitions(result_dir, method, MODEL_ORDER)
    primary_tests, descriptive_tests, excluded = build_m3_regime_rows(transitions)
    if primary_tests.empty:
        raise RuntimeError("No M=3 primary tests were constructed.")

    primary_tests = add_regime_bh_columns(primary_tests, alpha)
    descriptive_tests = add_descriptive_wald_columns(descriptive_tests)
    primary_tests, descriptive_tests = add_dvss_bh_followup_columns(
        primary_tests, descriptive_tests, alpha
    )
    regime_summary, cohort_summary, excluded_summary = summarize_primary_tests(
        primary_tests, excluded, alpha
    )
    bh_results_table = make_bh_results_table(regime_summary)
    descriptive_summary = summarize_descriptive_tests(descriptive_tests, alpha)
    revised_reported, wt_dvss_bh, wt_dvss_ci = make_interaction_subset_tables(primary_tests)
    manuscript_subset = primary_tests.loc[primary_tests["manuscript_subset"]].copy()
    manuscript_descriptive_subset = descriptive_tests.loc[
        descriptive_tests["manuscript_subset"]
    ].copy()

    primary_tests.to_csv(output_root / "M3_regime_primary_wt_tests.csv", index=False)
    descriptive_tests.to_csv(
        output_root / "M3_regime_descriptive_single_reference_tests.csv",
        index=False,
    )
    excluded.to_csv(output_root / "M3_regime_excluded_triad_targets.csv", index=False)
    manuscript_subset.to_csv(
        output_root / "M3_regime_manuscript_subset_wt_tests.csv",
        index=False,
    )
    manuscript_descriptive_subset.to_csv(
        output_root / "M3_regime_manuscript_subset_single_reference_descriptive.csv",
        index=False,
    )
    revised_reported.to_csv(
        output_root / "M3_regime_revised_reported_significant_interactions.csv",
        index=False,
    )
    wt_dvss_bh.to_csv(
        output_root / "M3_regime_wt_and_dvss_bh_significant_interactions.csv",
        index=False,
    )
    wt_dvss_ci.to_csv(
        output_root / "M3_regime_wt_and_dvss_ci_nonoverlap_interactions.csv",
        index=False,
    )
    regime_summary.to_csv(output_root / "M3_regime_summary.csv", index=False)
    bh_results_table.to_csv(output_root / "M3_regime_bh_results_table.csv", index=False)
    cohort_summary.to_csv(output_root / "M3_regime_summary_by_cohort.csv", index=False)
    excluded_summary.to_csv(output_root / "M3_regime_excluded_summary.csv", index=False)
    descriptive_summary.to_csv(
        output_root / "M3_regime_descriptive_single_reference_summary.csv",
        index=False,
    )

    write_manifest(
        output_root,
        run_root,
        result_dir,
        method,
        primary_tests,
        descriptive_tests,
        excluded,
    )

    return {
        "primary_tests": primary_tests,
        "descriptive_tests": descriptive_tests,
        "excluded": excluded,
        "regime_summary": regime_summary,
        "bh_results_table": bh_results_table,
        "cohort_summary": cohort_summary,
        "excluded_summary": excluded_summary,
        "descriptive_summary": descriptive_summary,
        "revised_reported": revised_reported,
        "wt_dvss_bh": wt_dvss_bh,
        "wt_dvss_ci": wt_dvss_ci,
    }


def main() -> None:
    """Purpose: run M=3 regime correction end to end. Inputs: CLI options. Outputs: task-local CSVs and manifest. Assumptions: M=3 CSV exports were generated from task0 model outputs."""

    args = parse_args()
    run_root, result_dir = resolve_result_dir(args.run_root, args.method)
    output_root = Path(args.output_root).resolve()
    results = run_m3_regime_analysis(
        result_dir=result_dir,
        output_root=output_root,
        run_root=run_root,
        method=args.method,
        alpha=args.alpha,
    )

    print(f"Wrote M=3 regime correction outputs to {output_root}")
    print(results["regime_summary"].to_string(index=False))


if __name__ == "__main__":
    main()
