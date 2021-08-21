import pandas as pd
import numpy as np
from count_combinations import compute_samples
from cancer_epistasis import numbers_positive_lambdas
from cancer_epistasis import estimate_lambdas
from cancer_epistasis import asymp_CI_lambdas
from locations import merged_maf_clinical_file_name
from locations import mutation_rates_file

## We fix M=3 so
numbers_positive_lambdas = numbers_positive_lambdas[3]

db = pd.read_csv(merged_maf_clinical_file_name)
db = db[db['Variant_Classification'] != 'Silent']
db = db[~pd.isnull(db['Mutation'])]

#ranked_gene_list = ['TP53', 'KRAS', 'EGFR', 'STK11', 'ATM', 'NF1', 'BRAF', 'PIK3CA', 'CDKN2A', 'ERBB4', 'APC', 'NTRK3', 'EPHA3', 'MET', 'EPHA5', 'RB1', 'ERBB2', 'KDR', 'ALK', 'SMAD4', 'GNAS', 'PDGFRA', 'PTEN', 'FLT1', 'CTNNB1', 'BRCA2', 'NOTCH1', 'NTRK1', 'ATR', 'NOTCH2', 'FLT4', 'RET', 'PIK3CG', 'NOTCH4', 'BRCA1', 'KIT', 'TSC2', 'SMO', 'NOTCH3', 'PIK3C2G', 'FLT3', 'JAK3', 'DDR2', 'EPHB1', 'PDGFRB', 'JAK2', 'NTRK2', 'BARD1', 'MAP3K1', 'FBXW7', 'WT1', 'TSHR', 'CSF1R', 'ABL1', 'FGFR4', 'ERBB3', 'CBL', 'IGF1R', 'IRS2', 'MSH6', 'FGFR1', 'RICTOR', 'FGFR2', 'NRAS', 'TSC1', 'PIK3R1', 'MSH2', 'CHEK2', 'BTK', 'PIK3C3', 'RUNX1', 'PTPN11', 'CDH1', 'BAP1', 'AKT3', 'ARAF', 'FGFR3', 'BCL6', 'AXL', 'ERG', 'RAF1', 'JAK1', 'NF2', 'SMAD2', 'TGFBR2', 'MPL', 'MLH1', 'MAP2K1', 'MAP3K13', 'MYCN', 'MDM4', 'MAP2K4', 'MEN1', 'CCNE1', 'IKBKE', 'AKT1', 'FGF3', 'CCND2', 'PIK3R2', 'ETV6', 'SYK', 'CHEK1', 'SRC', 'AKT2', 'MYC', 'GATA1', 'MDM2', 'REL', 'VHL', 'HRAS', 'MAP2K2', 'CDKN1B', 'CDK8', 'CDK6', 'PRKAR1A', 'GSK3B', 'SUFU', 'CCND1', 'FGF4', 'AURKB', 'JUN', 'CCND3', 'CDKN2C', 'CDK4', 'BCL2', 'CRKL', 'CDKN2B']
ranked_gene_list = ['TP53', 'KRAS', 'EGFR', 'STK11', 'ATM', 'NF1', 'BRAF', 'PIK3CA', 'ERBB4', 'APC', 'NTRK3', 'EPHA3', 'MET', 'EPHA5', 'RB1', 'ERBB2', 'KDR', 'ALK', 'SMAD4', 'GNAS', 'PDGFRA', 'PTEN', 'FLT1', 'CTNNB1', 'BRCA2', 'NOTCH1', 'NTRK1', 'ATR', 'NOTCH2', 'FLT4', 'RET', 'PIK3CG', 'NOTCH4', 'BRCA1', 'KIT', 'TSC2', 'SMO', 'NOTCH3', 'PIK3C2G', 'FLT3', 'JAK3', 'DDR2', 'EPHB1', 'PDGFRB', 'JAK2', 'NTRK2', 'BARD1', 'MAP3K1', 'FBXW7', 'WT1', 'TSHR', 'CSF1R', 'ABL1', 'FGFR4', 'ERBB3', 'CBL', 'IGF1R', 'IRS2', 'MSH6', 'FGFR1', 'RICTOR', 'FGFR2', 'NRAS', 'TSC1', 'PIK3R1', 'MSH2', 'CHEK2', 'BTK', 'PIK3C3', 'RUNX1', 'PTPN11', 'CDH1', 'BAP1', 'AKT3', 'ARAF', 'FGFR3', 'BCL6', 'AXL', 'ERG', 'RAF1', 'JAK1', 'NF2', 'SMAD2', 'TGFBR2', 'MPL', 'MLH1', 'MAP2K1', 'MAP3K13', 'MYCN', 'MDM4', 'MAP2K4', 'MEN1', 'CCNE1', 'IKBKE', 'AKT1', 'FGF3', 'CCND2', 'PIK3R2', 'ETV6', 'SYK', 'CHEK1', 'SRC', 'AKT2', 'MYC', 'GATA1', 'MDM2', 'REL', 'VHL', 'HRAS', 'MAP2K2', 'CDKN1B', 'CDK8', 'CDK6', 'PRKAR1A', 'GSK3B', 'SUFU', 'CCND1', 'FGF4', 'AURKB', 'JUN', 'CCND3', 'CDKN2C', 'CDK4', 'BCL2', 'CRKL', 'CDKN2B']


#mutation_rates_file = pd.read_csv(mutation_rates_file)
#mutation_rates = dict(zip(mutation_rates_file['gene'], mutation_rates_file['rate']))

def main(samples, bounds, draws = 1):
    print(f"Bounds: {bounds}")
    MAP = estimate_lambdas(samples, draws=1,
                        upper_bound_prior=bounds,
                        kwargs={'return_raw':True})
    print(f"MAP: {MAP[0]['lambdas']}")
    #hit_bound = np.isclose(MAP[0]['lambdas'], bounds)
    bound_changes = 0
    while MAP[1].fun == -1e+100:
        bound_changes += 1
        prop_bounds = np.array([1 / (2 ** bound_changes), 1, 1, 1 / (2 ** bound_changes), 1, 1, 1 / (2 ** bound_changes), 1, 1, 1, 1, 1 / (2 ** bound_changes)])
        #prop_bounds = np.array([1 * (bound_changes + 1), 1, 1, 1 * (bound_changes + 1), 1, 1, 1 * (bound_changes + 1), 1, 1, 1, 1, 1 * (bound_changes + 1)])

        print(f"Proposed bounds: {prop_bounds}")
        prop_MAP = estimate_lambdas(samples, draws=1,
                                    upper_bound_prior=prop_bounds,
                                    kwargs={'return_raw':True})
        print(f"Proposed MAP: {prop_MAP[0]['lambdas']}")
        if prop_MAP[1].fun > MAP[1].fun:
            MAP = prop_MAP
            bounds = prop_bounds
            #hit_bound = (MAP[0]['lambdas'] == bounds)
        else:
            print("Proposed bounds have lower likelihood")
    return MAP[0]['lambdas']

bounds = 1

with open('3_genes_1_bound_lambdas.txt','a') as lambda_output:
    lambda_output.write('3rd_gene,lambda0,lambda1,lambda2,lambda3,lambda4,lambda5,lambda6,lambda7,lambda8,lambda9,lambda10,lambda11,lambda0_ci_low,lambda0_ci_high,lambda1_ci_low,lambda1_ci_high,lambda2_ci_low,lambda2_ci_high,lambda3_ci_low,lambda3_ci_high,lambda4_ci_low,lambda4_ci_high,lambda5_ci_low,lambda5_ci_high,lambda6_ci_low,lambda6_ci_high,lambda7_ci_low,lambda7_ci_high,lambda8_ci_low,lambda8_ci_high,lambda9_ci_low,lambda9_ci_high,lambda10_ci_low,lambda10_ci_high,lambda11_ci_low,lambda11_ci_high,' + '\n')
    for i in ranked_gene_list[2:]:
        genes = ['TP53', 'KRAS', i]
        samples = compute_samples(db, mutations=genes)

        lambda_output.write(i + ',')
        lambda_list = main(samples, bounds)

        for j in lambda_list:
            lambda_output.write(str(j) + ',')
        for k in asymp_CI_lambdas(lambda_list, samples):
            lambda_output.write(str(k[0]) + ',' + str(k[1]) + ',')
        lambda_output.write('\n')