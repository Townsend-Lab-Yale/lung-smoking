"""
Shared loaders for the GENIE BPC NSCLC v2.0-public LUAD survival analysis.

All loaders return data restricted to the LUAD cohort defined as:
    cpt_oncotree_code == "LUAD"  (sample-level OncoTree)
intersected at the patient level with:
    ca_hist_adeno_squamous == "Adenocarcinoma"  (curator-assigned histology).

Sample-per-patient rule: earliest primary tumor by `dob_cpt_report_days`;
fall back to earliest non-primary if no primary exists for that patient
(metastatic inclusion is documented in the plan as a noted caveat — captured
by the Sample_metastatic indicator).

Time scales are in months throughout. OS clock is from cancer diagnosis
(`tt_os_dx_mos`); left-truncation entry time is `dx_cpt_rep_days/30.44`.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

DATA_DIR = Path(__file__).resolve().parents[2] / "data" / "genie_bpc_nsclc_v2.0"

DAYS_PER_MONTH = 30.4375  # average month length used by BPC for *_mos conversions

# Variant_Classification values that count as a "mutation" for genotype calls.
# Standard non-silent coding filter; matches plan §2.
MUTATION_CLASSES = {
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Splice_Site",
    "Translation_Start_Site",
    "Nonstop_Mutation",
}

# Drug name → class (case-insensitive substring match against `regimen_drugs`).
PLATINUM_AGENTS = {"Carboplatin", "Cisplatin"}
PLATINUM_PARTNERS = {
    "Pemetrexed", "Paclitaxel", "Docetaxel", "Gemcitabine",
    "Vinorelbine", "Etoposide",
}
IO_AGENTS = {
    "Pembrolizumab", "Nivolumab", "Atezolizumab", "Durvalumab",
    "Avelumab", "Ipilimumab", "Cemiplimab",
}
TKI_AGENTS = {
    "Erlotinib", "Gefitinib", "Afatinib", "Osimertinib", "Dacomitinib",
    "Crizotinib", "Ceritinib", "Alectinib", "Brigatinib", "Lorlatinib",
    "Entrectinib", "Dabrafenib", "Trametinib", "Capmatinib", "Tepotinib",
    "Selpercatinib", "Pralsetinib", "Sotorasib", "Adagrasib",
}

ADVANCED_STAGES = {"Stage III", "Stage IV"}
EVER_SMOKER_LEVELS = {
    "Current user", "Former user (quit <1 year)",
    "Former user (quit >1 year)", "Former user (unknown time)",
}
NEVER_SMOKER_LEVELS = {"Never used"}

# Institution dummies; reference = MSK (largest).
INST_LEVELS = ["DFCI", "MSK", "VICC", "UHN"]


# -------------------- Treatment flags --------------------

def _drug_str_contains(drug_str: str, agents: set[str]) -> bool:
    if not isinstance(drug_str, str):
        return False
    s = drug_str.lower()
    return any(a.lower() in s for a in agents)


def build_treatment_flags(regimen_df: pd.DataFrame) -> pd.DataFrame:
    """Per (record_id, ca_seq), aggregate ever-received treatment indicators
    across all regimens for that index cancer.

    Returns columns:
        record_id, ca_seq, Tx_platinum_doublet, Tx_IO, Tx_targeted_TKI,
        N_lines_therapy
    """
    df = regimen_df.copy()

    has_platinum = df["regimen_drugs"].apply(
        lambda s: _drug_str_contains(s, PLATINUM_AGENTS))
    has_partner = df["regimen_drugs"].apply(
        lambda s: _drug_str_contains(s, PLATINUM_PARTNERS))
    df["_is_platinum_doublet"] = (has_platinum & has_partner).astype(int)
    df["_is_IO"] = df["regimen_drugs"].apply(
        lambda s: int(_drug_str_contains(s, IO_AGENTS)))
    df["_is_TKI"] = df["regimen_drugs"].apply(
        lambda s: int(_drug_str_contains(s, TKI_AGENTS)))

    agg = df.groupby(["record_id", "ca_seq"]).agg(
        Tx_platinum_doublet=("_is_platinum_doublet", "max"),
        Tx_IO=("_is_IO", "max"),
        Tx_targeted_TKI=("_is_TKI", "max"),
        N_lines_therapy=("regimen_number_within_cancer", "max"),
    ).reset_index()
    return agg


# -------------------- Clinical loader --------------------

def load_clinical() -> pd.DataFrame:
    """Build one row per LUAD patient with harmonized clinical/treatment data.

    Selection rules:
        1. Index cancer with ca_hist_adeno_squamous == "Adenocarcinoma".
        2. CPT row with cpt_oncotree_code == "LUAD" for the same (record_id, ca_seq).
        3. One sample per record_id: prefer Primary tumor, then earliest
           dob_cpt_report_days. Tag Sample_metastatic accordingly.
        4. Drop rows missing OS, age, sex, or stage; drop rows where
           dx_cpt_rep_days >= tt_os_dx_days (left-truncation impossible).

    Returns DataFrame with columns:
        Sample_ID, Patient_ID, OS_time_months, OS_event, Entry_time_months,
        Age, Sex_male, Stage_advanced, Smoking_ever, Smoking_detailed,
        Sample_metastatic, Institution, SEQ_ASSAY_ID,
        Tx_platinum_doublet, Tx_IO, Tx_targeted_TKI, N_lines_therapy
    """
    patient = pd.read_csv(DATA_DIR / "patient_level_dataset.csv", engine="python")
    cancer = pd.read_csv(
        DATA_DIR / "cancer_level_dataset_index.csv",
        engine="python",
    )
    cpt = pd.read_csv(
        DATA_DIR / "cancer_panel_test_level_dataset.csv",
        engine="python",
    )
    # The GENIE BPC regimen file can contain rows that the default C parser
    # handles inconsistently across local exports; the Python engine is slower
    # but materially more robust for this treatment-aggregation input.
    regimen = pd.read_csv(
        DATA_DIR / "regimen_cancer_level_dataset.csv",
        engine="python",
    )

    # 1. LUAD-adeno index cancers.
    luad_cancer = cancer[
        cancer["ca_hist_adeno_squamous"] == "Adenocarcinoma"
    ].copy()

    # 2. CPT rows must be LUAD on OncoTree; join with LUAD-adeno cancers.
    luad_cpt = cpt[cpt["cpt_oncotree_code"] == "LUAD"].copy()
    cpt_cancer = luad_cpt.merge(
        luad_cancer[[
            "record_id", "ca_seq", "age_dx", "stage_dx", "ca_lung_cigarette",
            "os_dx_status", "tt_os_dx_mos",
        ]],
        on=["record_id", "ca_seq"],
        how="inner",
    )

    # 3. Sample selection: one per record_id.
    cpt_cancer["is_primary"] = (
        cpt_cancer["sample_type"] == "Primary tumor").astype(int)
    # Sort: primary first, then earliest sequencing report.
    cpt_cancer = cpt_cancer.sort_values(
        ["record_id", "is_primary", "dob_cpt_report_days"],
        ascending=[True, False, True],
    )
    one_per_pt = cpt_cancer.drop_duplicates("record_id", keep="first").copy()

    # 4. Build covariates.
    one_per_pt["Sample_metastatic"] = (one_per_pt["is_primary"] == 0).astype(int)
    one_per_pt["Stage_advanced"] = one_per_pt["stage_dx"].isin(
        ADVANCED_STAGES).astype(int)

    # Smoking: ever / never; "Unknown" → NaN (drop downstream).
    smoke = one_per_pt["ca_lung_cigarette"]
    one_per_pt["Smoking_ever"] = np.where(
        smoke.isin(EVER_SMOKER_LEVELS), 1,
        np.where(smoke.isin(NEVER_SMOKER_LEVELS), 0, np.nan),
    )
    one_per_pt["Smoking_detailed"] = smoke

    # Sex from patient table.
    one_per_pt = one_per_pt.merge(
        patient[["record_id", "naaccr_sex_code"]],
        on="record_id", how="left",
    )
    one_per_pt["Sex_male"] = (
        one_per_pt["naaccr_sex_code"] == "Male").astype(int)

    # OS clock and left-truncation entry.
    one_per_pt["OS_time_months"] = one_per_pt["tt_os_dx_mos"]
    one_per_pt["OS_event"] = one_per_pt["os_dx_status"].astype(int)
    one_per_pt["Entry_time_months"] = (
        one_per_pt["dx_cpt_rep_days"] / DAYS_PER_MONTH)

    # Drop QC-failing rows.
    before = len(one_per_pt)
    one_per_pt = one_per_pt.dropna(subset=[
        "OS_time_months", "OS_event", "age_dx", "Sex_male",
        "Stage_advanced", "Smoking_ever", "Entry_time_months",
    ])
    one_per_pt = one_per_pt[
        one_per_pt["Entry_time_months"] < one_per_pt["OS_time_months"]
    ]
    print(f"  load_clinical: {before} → {len(one_per_pt)} samples after QC")

    # Treatment.
    tx = build_treatment_flags(regimen)
    one_per_pt = one_per_pt.merge(tx, on=["record_id", "ca_seq"], how="left")
    for c in ["Tx_platinum_doublet", "Tx_IO", "Tx_targeted_TKI"]:
        one_per_pt[c] = one_per_pt[c].fillna(0).astype(int)
    one_per_pt["N_lines_therapy"] = one_per_pt["N_lines_therapy"].fillna(0)

    # Final select / rename.
    out = pd.DataFrame({
        "Sample_ID": one_per_pt["cpt_genie_sample_id"],
        "Patient_ID": one_per_pt["record_id"],
        "OS_time_months": one_per_pt["OS_time_months"],
        "OS_event": one_per_pt["OS_event"],
        "Entry_time_months": one_per_pt["Entry_time_months"],
        "Age": one_per_pt["age_dx"],
        "Sex_male": one_per_pt["Sex_male"],
        "Stage_advanced": one_per_pt["Stage_advanced"],
        "Smoking_ever": one_per_pt["Smoking_ever"].astype(int),
        "Smoking_detailed": one_per_pt["Smoking_detailed"],
        "Sample_metastatic": one_per_pt["Sample_metastatic"],
        "Institution": one_per_pt["institution"],
        "SEQ_ASSAY_ID": one_per_pt["cpt_seq_assay_id"],
        "Tx_platinum_doublet": one_per_pt["Tx_platinum_doublet"],
        "Tx_IO": one_per_pt["Tx_IO"],
        "Tx_targeted_TKI": one_per_pt["Tx_targeted_TKI"],
        "N_lines_therapy": one_per_pt["N_lines_therapy"],
    })

    # Institution dummies (ref = largest present institution).
    inst_counts = out["Institution"].value_counts()
    ref_inst = inst_counts.idxmax()
    for inst in inst_counts.index:
        if inst != ref_inst:
            out[f"Inst_{inst}"] = (out["Institution"] == inst).astype(int)
    print(f"  load_clinical: institution dummies (ref={ref_inst}): "
          f"{[c for c in out.columns if c.startswith('Inst_')]}")

    return out.reset_index(drop=True)


# -------------------- Panel coverage --------------------

def load_panel_gene_lists() -> dict[str, set[str]]:
    """Parse each data_gene_panel_*.txt into {SEQ_ASSAY_ID -> set(genes)}."""
    panels = {}
    for f in sorted(DATA_DIR.glob("data_gene_panel_*.txt")):
        with open(f) as fh:
            stable_id = None
            genes: set[str] = set()
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("stable_id:"):
                    stable_id = line.split(":", 1)[1].strip()
                elif line.startswith("gene_list:"):
                    rest = line.split(":", 1)[1]
                    genes = {g.strip() for g in rest.split("\t") if g.strip()}
        if stable_id is None:
            raise ValueError(f"No stable_id in {f}")
        panels[stable_id] = genes
    return panels


def load_panel_coverage(samples: pd.Series, genes: list[str]) -> pd.DataFrame:
    """Return sample × gene boolean DataFrame indicating which (sample, gene)
    pairs are covered by the sample's mutation panel.

    Args:
        samples: Sample IDs to include (rows).
        genes: Gene symbols to include (columns).
    """
    panels = load_panel_gene_lists()
    matrix = pd.read_csv(DATA_DIR / "data_gene_matrix.txt", sep="\t")
    matrix = matrix.set_index("SAMPLE_ID")
    # Only need mutation-panel coverage (CNA panel handled separately).
    sample_to_panel = matrix["mutations"].to_dict()

    rows = []
    idx = []
    for s in samples:
        panel = sample_to_panel.get(s)
        covered = panels.get(panel, set()) if panel else set()
        rows.append([g in covered for g in genes])
        idx.append(s)
    return pd.DataFrame(rows, index=idx, columns=genes, dtype=bool)


# -------------------- Genotype matrix --------------------

def load_genotype_matrix(
    samples: pd.Series,
    genes: list[str],
    coverage: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build sample × gene mutation indicator matrix from the MAF.

    Cells are 1 if the sample carries any qualifying variant in the gene,
    0 if covered by the panel but no qualifying variant, NaN if the gene
    is not on the sample's mutation panel (so it cannot be called).

    Args:
        samples: Sample IDs (rows).
        genes: Gene symbols (cols).
        coverage: Optional precomputed sample×gene boolean coverage.
            If None, computed via load_panel_coverage().
    """
    maf = pd.read_csv(
        DATA_DIR / "data_mutations_extended.txt", sep="\t",
        usecols=["Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode"],
        low_memory=False,
    )
    # Normalize the case-mismatch noted in the plan.
    maf["Variant_Classification"] = maf["Variant_Classification"].replace(
        {"Frame_Shift_DEL": "Frame_Shift_Del"})
    maf = maf[maf["Variant_Classification"].isin(MUTATION_CLASSES)]
    maf = maf[maf["Hugo_Symbol"].isin(genes)]
    maf = maf[maf["Tumor_Sample_Barcode"].isin(set(samples))]

    # Pivot: sample × gene boolean (any qualifying mutation).
    hits = maf.assign(_hit=1).pivot_table(
        index="Tumor_Sample_Barcode",
        columns="Hugo_Symbol",
        values="_hit",
        aggfunc="max",
        fill_value=0,
    )
    # Reindex to full sample × gene.
    hits = hits.reindex(index=list(samples), columns=genes, fill_value=0)

    # Mask non-covered genes to NaN.
    if coverage is None:
        coverage = load_panel_coverage(samples, genes)
    coverage = coverage.reindex(index=list(samples), columns=genes, fill_value=False)
    geno = hits.astype(float).where(coverage, other=np.nan)
    geno.index.name = "Sample_ID"
    return geno


# -------------------- CNA matrix --------------------

def load_cna_matrix(samples: pd.Series, genes: list[str]) -> pd.DataFrame:
    """Sample × gene matrix encoding |GISTIC|==2 events as 1 (deep deletion or
    high-level amplification), 0 if covered with no event, NaN if not covered.

    GENIE BPC `data_CNA.txt` is gene × sample with NaN already used for
    non-coverage; we transpose and binarize.
    """
    cna = pd.read_csv(DATA_DIR / "data_CNA.txt", sep="\t", index_col=0)
    cna = cna.T  # sample × gene
    keep_genes = [g for g in genes if g in cna.columns]
    cna = cna.reindex(index=list(samples), columns=keep_genes)
    binary = (cna.abs() == 2).astype(float)
    binary = binary.where(cna.notna(), other=np.nan)
    # Reintroduce missing genes as all-NaN columns.
    for g in genes:
        if g not in binary.columns:
            binary[g] = np.nan
    binary = binary[genes]
    binary.index.name = "Sample_ID"
    return binary


# -------------------- Convenience --------------------

def m2_pair_genes(m2_path: Path | str) -> list[str]:
    """Return sorted unique gene symbols appearing in either column of
    M2_gene_gammas.csv."""
    df = pd.read_csv(m2_path)
    return sorted(set(df["first_gene"]) | set(df["second_gene"]))
