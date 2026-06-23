"""
01_prepare_genie_cohort.py — restricted input-preparation workflow for the
GENIE BPC LUAD survival analysis.

This entrypoint combines:
    1. Synapse download of the GENIE BPC NSCLC v2.0-public release files
       needed for the survival analysis, and
    2. construction of the LUAD survival cohort, genotype matrix, panel
       coverage matrix, directional epistasis table, and pair-feasibility
       table used by downstream Cox models.

Outputs (under code/survival_analysis/output/):
    survival_table.csv
    genotype_matrix.csv
    panel_coverage.csv
    epistasis_direction.csv
    feasibility_table.csv
"""

from __future__ import annotations

import getpass
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import helpers_genie as h  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data" / "genie_bpc_nsclc_v2.0"
OUTDIR = Path(__file__).resolve().parent / "output"

MIN_DOUBLE_MUTANT = 10
BOUNDARY_EPS = 1e-10  # gamma MLE below this is treated as boundary (gamma≈0)

# Synapse IDs from the BPC NSCLC v2.0-public release file index.
FILES = {
    "syn30358096": "patient_level_dataset.csv",
    "syn30358090": "cancer_level_dataset_index.csv",
    "syn30358091": "cancer_level_dataset_non_index.csv",
    "syn30358092": "cancer_panel_test_level_dataset.csv",
    "syn30358097": "regimen_cancer_level_dataset.csv",
    "syn30358095": "pathology_report_level_dataset.csv",
    "syn30358100": "data_clinical_patient.txt",
    "syn30358101": "data_clinical_sample.txt",
    "syn30358120": "data_mutations_extended.txt",
    "syn30358099": "data_CNA.txt",
    "syn30358105": "data_fusions.txt",
    "syn30358106": "data_gene_matrix.txt",
    "syn30358107": "data_gene_panel_DFCI-ONCOPANEL-1.txt",
    "syn30358108": "data_gene_panel_DFCI-ONCOPANEL-2.txt",
    "syn30358109": "data_gene_panel_DFCI-ONCOPANEL-3.txt",
    "syn30358110": "data_gene_panel_MSK-IMPACT341.txt",
    "syn30358112": "data_gene_panel_MSK-IMPACT410.txt",
    "syn30358113": "data_gene_panel_MSK-IMPACT468.txt",
    "syn30358114": "data_gene_panel_UHN-48-V1.txt",
    "syn30358116": "data_gene_panel_UHN-50-V2.txt",
    "syn39802620": "data_gene_panel_VICC-01-SOLIDTUMOR.txt",
    "syn30358118": "data_gene_panel_VICC-01-T5A.txt",
    "syn30358119": "data_gene_panel_VICC-01-T7.txt",
}


def resolve_m2_path(root: Path) -> Path:
    """Return the local path to `M2_gene_gammas.csv`.

    Inputs
    ------
    root : repository root

    Outputs
    -------
    Existing file path for the directional M2 epistasis estimates.

    Assumptions
    -----------
    The repository may store `variant_results/` either at repo root or under
    `code/`.
    """
    candidates = [
        root / "variant_results" / "M2_gene_gammas.csv",
        root / "code" / "variant_results" / "M2_gene_gammas.csv",
    ]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError("Could not locate M2_gene_gammas.csv in expected paths.")


def ensure_genie_download() -> None:
    """Download required GENIE BPC files when they are not already present.

    Inputs
    ------
    None. Uses the fixed Synapse file map defined above.

    Outputs
    -------
    Populates `data/genie_bpc_nsclc_v2.0/` with the required release files.

    Assumptions
    -----------
    The user has a Synapse account with the appropriate data-use agreement and
    can provide a personal access token interactively when files are missing.
    """
    missing = [fname for fname in FILES.values() if not (DATA_DIR / fname).exists()]
    if not missing:
        print("\n[1/6] GENIE BPC raw files already present; skipping download.")
        return

    try:
        import synapseclient
    except ImportError as exc:
        raise ImportError(
            "synapseclient is required to download missing GENIE BPC files."
        ) from exc

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    print("\n[1/6] Downloading missing GENIE BPC raw files from Synapse…")
    print(f"  missing files: {len(missing)}")

    token = getpass.getpass("Synapse personal access token: ").strip()
    if not token:
        raise SystemExit("No Synapse token provided.")

    syn = synapseclient.Synapse()
    syn.login(authToken=token)

    for syn_id, fname in FILES.items():
        out_path = DATA_DIR / fname
        if out_path.exists():
            continue
        print(f"  {syn_id}  {fname}")
        syn.get(syn_id, downloadLocation=str(DATA_DIR), ifcollision="overwrite.local")


def compute_epistasis_direction(m2_path: Path) -> pd.DataFrame:
    """Summarize directional M2 epistasis by CI non-overlap.

    Inputs
    ------
    m2_path : path to `M2_gene_gammas.csv`

    Outputs
    -------
    One row per (cohort, gene_A, gene_B) with gamma estimates in the wild-type
    and altered backgrounds, plus directional-significance flags.

    Assumptions
    -----------
    Significance is defined by non-overlap of the corresponding gamma
    confidence intervals, matching the directional epistasis rule used
    elsewhere in this project.
    """
    gammas = pd.read_csv(m2_path)
    cohort_map = {"smoking_plus": "ES", "nonsmoking_plus": "NS"}

    keys_needed = [
        "((0, 0), (1, 0))",
        "((0, 1), (1, 1))",
        "((0, 0), (0, 1))",
        "((1, 0), (1, 1))",
    ]

    records = []
    for key, kg in gammas.groupby("key"):
        if key not in cohort_map:
            continue
        m2_cohort = cohort_map[key]
        for (fg, sg), pdf in kg.groupby(["first_gene", "second_gene"]):
            d = {row["mutation"]: row for _, row in pdf.iterrows()}
            if not all(k in d for k in keys_needed):
                continue

            a_in_wt = d["((0, 0), (1, 0))"]
            a_in_b = d["((0, 1), (1, 1))"]
            b_in_wt = d["((0, 0), (0, 1))"]
            b_in_a = d["((1, 0), (1, 1))"]

            b_inf_a = (
                a_in_wt["gamma_ci_low"] > a_in_b["gamma_ci_high"]
                or a_in_b["gamma_ci_low"] > a_in_wt["gamma_ci_high"]
            )
            a_inf_b = (
                b_in_wt["gamma_ci_low"] > b_in_a["gamma_ci_high"]
                or b_in_a["gamma_ci_low"] > b_in_wt["gamma_ci_high"]
            )

            if a_inf_b and b_inf_a:
                focal = "both"
            elif a_inf_b:
                focal = "A->B"
            elif b_inf_a:
                focal = "B->A"
            else:
                focal = "none"

            records.append(
                {
                    "M2_cohort": m2_cohort,
                    "gene_A": fg,
                    "gene_B": sg,
                    "gamma_A_in_WT": a_in_wt["gamma_mle"],
                    "gamma_A_in_B": a_in_b["gamma_mle"],
                    "gamma_B_in_WT": b_in_wt["gamma_mle"],
                    "gamma_B_in_A": b_in_a["gamma_mle"],
                    "A_influences_B_selection": bool(a_inf_b),
                    "B_influences_A_selection": bool(b_inf_a),
                    "focal_direction": focal,
                    "has_boundary_mle": bool(
                        a_in_b["gamma_mle"] < BOUNDARY_EPS
                        or b_in_a["gamma_mle"] < BOUNDARY_EPS
                    ),
                    "any_significant": bool(a_inf_b or b_inf_a),
                }
            )

    return pd.DataFrame(records)


def assess_pair(
    geno_pair: pd.DataFrame,
    smoking: pd.Series,
    gene_a: str,
    gene_b: str,
    m2_cohort: str,
) -> dict:
    """Count genotype-group sample sizes for one gene pair.

    Inputs
    ------
    geno_pair : sample-indexed DataFrame with columns `[gene_a, gene_b]`
    smoking : sample-indexed ever-smoker indicator
    gene_a, gene_b : focal genes
    m2_cohort : cohort in which the pair was M2-significant (`ES` or `NS`)

    Outputs
    -------
    Dictionary of pooled, ever-smoker, and never-smoker genotype counts and
    testability flags.

    Assumptions
    -----------
    Samples with `NaN` in either gene are not dual-covered and are excluded
    from feasibility counts.
    """
    a = geno_pair[gene_a]
    b = geno_pair[gene_b]
    both_covered = a.notna() & b.notna()

    sub = pd.DataFrame(
        {
            "A": a[both_covered].astype(int),
            "B": b[both_covered].astype(int),
            "Smoking_ever": smoking.reindex(a.index)[both_covered].astype(int),
        }
    )

    def _counts(df: pd.DataFrame) -> dict:
        return {
            "n_total": len(df),
            "n_AB": int(((df["A"] == 1) & (df["B"] == 1)).sum()),
            "n_A_only": int(((df["A"] == 1) & (df["B"] == 0)).sum()),
            "n_B_only": int(((df["A"] == 0) & (df["B"] == 1)).sum()),
            "n_neither": int(((df["A"] == 0) & (df["B"] == 0)).sum()),
        }

    pooled = _counts(sub)
    es = _counts(sub[sub["Smoking_ever"] == 1])
    ns = _counts(sub[sub["Smoking_ever"] == 0])

    return {
        "gene_A": gene_a,
        "gene_B": gene_b,
        "M2_cohort": m2_cohort,
        "n_panel": pooled["n_total"],
        "n_AB": pooled["n_AB"],
        "n_A_only": pooled["n_A_only"],
        "n_B_only": pooled["n_B_only"],
        "n_neither": pooled["n_neither"],
        "testable_pooled": pooled["n_AB"] >= MIN_DOUBLE_MUTANT,
        "n_AB_ES": es["n_AB"],
        "n_total_ES": es["n_total"],
        "testable_ES": es["n_AB"] >= MIN_DOUBLE_MUTANT,
        "n_AB_NS": ns["n_AB"],
        "n_total_NS": ns["n_total"],
        "testable_NS": ns["n_AB"] >= MIN_DOUBLE_MUTANT,
    }


def build_survival_inputs() -> None:
    """Construct the restricted survival-analysis inputs from local raw data.

    Inputs
    ------
    None. Reads the previously downloaded GENIE BPC release files and the
    local M2 output file.

    Outputs
    -------
    Writes the cohort table, genotype matrix, panel coverage matrix,
    epistasis-direction table, and feasibility table under `output/`.

    Assumptions
    -----------
    Clinical QC, mutation calling, and panel-aware missingness are handled by
    `helpers_genie.py`, which defines the analysis cohort used throughout the
    restricted survival workflow.
    """
    OUTDIR.mkdir(parents=True, exist_ok=True)
    m2_path = resolve_m2_path(ROOT)

    print("\n[2/6] Loading harmonized clinical table…")
    clin = h.load_clinical()
    print(f"  retained {len(clin)} samples (one per LUAD-adeno patient)")

    print("\n[3/6] Loading M2 epistasis directional flags…")
    direction = compute_epistasis_direction(m2_path)
    direction.to_csv(OUTDIR / "epistasis_direction.csv", index=False)
    print(f"  {len(direction)} (gene_A, gene_B, cohort) M2 entries")
    print(f"  any_significant: {direction['any_significant'].sum()}")
    print(
        f"  focal_direction counts: "
        f"{direction['focal_direction'].value_counts().to_dict()}"
    )

    sig_pairs = direction[direction["any_significant"]].copy()
    gene_universe = sorted(set(sig_pairs["gene_A"]) | set(sig_pairs["gene_B"]))
    print(f"  gene universe for significant pairs: {len(gene_universe)} genes")

    print("\n[4/6] Building panel coverage and genotype matrices…")
    coverage = h.load_panel_coverage(clin["Sample_ID"], gene_universe)
    coverage.to_csv(OUTDIR / "panel_coverage.csv")
    print(f"  panel_coverage: {coverage.shape}")
    print(f"  mean coverage across (sample, gene) cells: {coverage.values.mean():.3f}")

    geno = h.load_genotype_matrix(clin["Sample_ID"], gene_universe, coverage=coverage)
    geno.to_csv(OUTDIR / "genotype_matrix.csv")
    print(
        f"  genotype_matrix: {geno.shape}, "
        f"NaN frac (non-covered cells): {geno.isna().values.mean():.3f}"
    )

    print("\n[5/6] Writing survival_table.csv (clinical + genotype merged)…")
    surv = clin.merge(geno, left_on="Sample_ID", right_index=True, how="left")
    surv.to_csv(OUTDIR / "survival_table.csv", index=False)
    print(f"  survival_table: {surv.shape}")

    print("\n[6/6] Per-pair feasibility…")
    smoking = clin.set_index("Sample_ID")["Smoking_ever"]
    feas_records = []
    for _, row in sig_pairs.iterrows():
        gene_a = row["gene_A"]
        gene_b = row["gene_B"]
        if gene_a not in geno.columns or gene_b not in geno.columns:
            continue
        rec = assess_pair(geno[[gene_a, gene_b]], smoking, gene_a, gene_b, row["M2_cohort"])
        rec["has_boundary_mle"] = bool(row["has_boundary_mle"])
        rec["focal_direction"] = row["focal_direction"]
        feas_records.append(rec)

    feas = pd.DataFrame(feas_records)
    feas.to_csv(OUTDIR / "feasibility_table.csv", index=False)
    print(f"  feasibility_table: {len(feas)} significant pairs evaluated")
    print(f"  testable (pooled, n_AB >= {MIN_DOUBLE_MUTANT}): {int(feas['testable_pooled'].sum())}")
    print(f"  testable (ES):  {int(feas['testable_ES'].sum())}")
    print(f"  testable (NS):  {int(feas['testable_NS'].sum())}")


def main() -> None:
    """Run the restricted input-preparation workflow end to end."""
    print("=" * 72)
    print("GENIE BPC LUAD — download and cohort construction")
    print("=" * 72)
    ensure_genie_download()
    build_survival_inputs()
    print("\nDone. Outputs in", OUTDIR)


if __name__ == "__main__":
    main()
