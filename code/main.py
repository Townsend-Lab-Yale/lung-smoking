import pandas as pd
import numpy as np
from count_combinations import compute_samples
from cancer_epistasis import numbers_positive_lambdas
from cancer_epistasis import estimate_lambdas
from cancer_epistasis import asymp_CI_lambdas
from cancer_epistasis import convert_lambdas_to_dict

from locations import merged_maf_clinical_file_name
from locations import mutation_rates_file
from locations import fluxes_mles_file_name
from locations import fluxes_cis_file_name

## We fix M=3 so
numbers_positive_lambdas = numbers_positive_lambdas[3]

db = pd.read_csv(merged_maf_clinical_file_name)
db = db[db['Variant_Classification'] != 'Silent']
db = db[~pd.isnull(db['Mutation'])]

#ranked_gene_list = ['TP53', 'KRAS', 'EGFR', 'STK11', 'ATM', 'NF1', 'BRAF', 'PIK3CA', 'CDKN2A', 'ERBB4', 'APC', 'NTRK3', 'EPHA3', 'MET', 'EPHA5', 'RB1', 'ERBB2', 'KDR', 'ALK', 'SMAD4', 'GNAS', 'PDGFRA', 'PTEN', 'FLT1', 'CTNNB1', 'BRCA2', 'NOTCH1', 'NTRK1', 'ATR', 'NOTCH2', 'FLT4', 'RET', 'PIK3CG', 'NOTCH4', 'BRCA1', 'KIT', 'TSC2', 'SMO', 'NOTCH3', 'PIK3C2G', 'FLT3', 'JAK3', 'DDR2', 'EPHB1', 'PDGFRB', 'JAK2', 'NTRK2', 'BARD1', 'MAP3K1', 'FBXW7', 'WT1', 'TSHR', 'CSF1R', 'ABL1', 'FGFR4', 'ERBB3', 'CBL', 'IGF1R', 'IRS2', 'MSH6', 'FGFR1', 'RICTOR', 'FGFR2', 'NRAS', 'TSC1', 'PIK3R1', 'MSH2', 'CHEK2', 'BTK', 'PIK3C3', 'RUNX1', 'PTPN11', 'CDH1', 'BAP1', 'AKT3', 'ARAF', 'FGFR3', 'BCL6', 'AXL', 'ERG', 'RAF1', 'JAK1', 'NF2', 'SMAD2', 'TGFBR2', 'MPL', 'MLH1', 'MAP2K1', 'MAP3K13', 'MYCN', 'MDM4', 'MAP2K4', 'MEN1', 'CCNE1', 'IKBKE', 'AKT1', 'FGF3', 'CCND2', 'PIK3R2', 'ETV6', 'SYK', 'CHEK1', 'SRC', 'AKT2', 'MYC', 'GATA1', 'MDM2', 'REL', 'VHL', 'HRAS', 'MAP2K2', 'CDKN1B', 'CDK8', 'CDK6', 'PRKAR1A', 'GSK3B', 'SUFU', 'CCND1', 'FGF4', 'AURKB', 'JUN', 'CCND3', 'CDKN2C', 'CDK4', 'BCL2', 'CRKL', 'CDKN2B']
ranked_gene_list = ['TP53', 'KRAS', 'EGFR', 'STK11', 'ATM', 'NF1', 'BRAF', 'PIK3CA', 'ERBB4', 'APC', 'NTRK3', 'EPHA3', 'MET', 'EPHA5', 'RB1', 'ERBB2', 'KDR', 'ALK', 'SMAD4', 'GNAS', 'PDGFRA', 'PTEN', 'FLT1', 'CTNNB1', 'BRCA2', 'NOTCH1', 'NTRK1', 'ATR', 'NOTCH2', 'FLT4', 'RET', 'PIK3CG', 'NOTCH4', 'BRCA1', 'KIT', 'TSC2', 'SMO', 'NOTCH3', 'PIK3C2G', 'FLT3', 'JAK3', 'DDR2', 'EPHB1', 'PDGFRB', 'JAK2', 'NTRK2', 'BARD1', 'MAP3K1', 'FBXW7', 'WT1', 'TSHR', 'CSF1R', 'ABL1', 'FGFR4', 'ERBB3', 'CBL', 'IGF1R', 'IRS2', 'MSH6', 'FGFR1', 'RICTOR', 'FGFR2', 'NRAS', 'TSC1', 'PIK3R1', 'MSH2', 'CHEK2', 'BTK', 'PIK3C3', 'RUNX1', 'PTPN11', 'CDH1', 'BAP1', 'AKT3', 'ARAF', 'FGFR3', 'BCL6', 'AXL', 'ERG', 'RAF1', 'JAK1', 'NF2', 'SMAD2', 'TGFBR2', 'MPL', 'MLH1', 'MAP2K1', 'MAP3K13', 'MYCN', 'MDM4', 'MAP2K4', 'MEN1', 'CCNE1', 'IKBKE', 'AKT1', 'FGF3', 'CCND2', 'PIK3R2', 'ETV6', 'SYK', 'CHEK1', 'SRC', 'AKT2', 'MYC', 'GATA1', 'MDM2', 'REL', 'VHL', 'HRAS', 'MAP2K2', 'CDKN1B', 'CDK8', 'CDK6', 'PRKAR1A', 'GSK3B', 'SUFU', 'CCND1', 'FGF4', 'AURKB', 'JUN', 'CCND3', 'CDKN2C', 'CDK4', 'BCL2', 'CRKL', 'CDKN2B']


#mutation_rates_file = pd.read_csv(mutation_rates_file)
#mutation_rates = dict(zip(mutation_rates_file['gene'], mutation_rates_file['rate']))

def lambdas_from_samples(samples):
    bounds = 1
    print(f"Bounds for fluxes: {bounds}")
    MLE = estimate_lambdas(samples, draws=bounds,
                           upper_bound_prior=1,
                           kwargs={'return_raw':True})
    print(f"MLE: {MLE[0]['lambdas']}")

    bound_changes = 0
    while MLE[1].fun < -1e+20:
        print(f"Algorithm did not converge changing bounds...")
        bound_changes += 1
        bounds = np.array(
            [1/(2**bound_changes),  # (0, 0, 0) -> (0, 0, 1)
             1,                     # (0, 0, 0) -> (0, 1, 0)
             1,                     # (0, 0, 1) -> (0, 1, 1)
             1/(2**bound_changes),  # (0, 1, 0) -> (0, 1, 1)
             1,                     # (0, 0, 0) -> (1, 0, 0)
             1,                     # (0, 0, 1) -> (1, 0, 1)
             1/(2**bound_changes),  # (1, 0, 0) -> (1, 0, 1)
             1,                     # (0, 1, 0) -> (1, 1, 0)
             1,                     # (1, 0, 0) -> (1, 1, 0)
             1,                     # (0, 1, 1) -> (1, 1, 1)
             1,                     # (1, 0, 1) -> (1, 1, 1)
             1/(2**bound_changes)]) # (1, 1, 0) -> (1, 1, 1)
        print(f"Proposed bounds: {bounds}")
        MLE = estimate_lambdas(samples, draws=1,
                                    upper_bound_prior=bounds,
                                    kwargs={'return_raw':True})
        print(f"MLE: {MLE[0]['lambdas']}")

    return MLE[0]


def main():
    lambdas_mles = {}
    lambdas_cis = {}

    for i, gene in enumerate(ranked_gene_list[2:]):
        print(f"Running model with third gene {gene} "
              f"(gene number {i}/{len(ranked_gene_list[2:])})")
        genes = ['TP53', 'KRAS', gene]
        samples = compute_samples(db, mutations=genes)

        print("Estimating fluxes...")
        mle = lambdas_from_samples(samples)
        lambdas_mles[gene] = convert_lambdas_to_dict(mle)
        np.save(fluxes_mles_file_name, lambdas_mles)

        print("Estimating asymptotic confidence intervals...")
        cis = asymp_CI_lambdas(mle['lambdas'], samples)
        lambdas_cis[gene] = convert_lambdas_to_dict(cis)
        np.save(fluxes_cis_file_name, lambdas_cis)

        print("")

    return lambdas_mles, lambdas_cis


if __name__ == "__main__":
    main()
