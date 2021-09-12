import pandas as pd
import numpy as np
from pandas._libs.missing import NA
from count_combinations import compute_samples
from count_combinations import are_all_fluxes_computable
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

genie_panels_used = pd.read_csv('../gene_panels/genie_panels_used.txt').set_index('Sample Identifier').to_dict()
msk2017_panels_used = pd.read_csv('../gene_panels/msk2017_panels_used.txt').set_index('Sample Identifier').to_dict()
msk2018_panels_used = pd.read_csv('../gene_panels/msk2018_panels_used.txt').set_index('Sample Identifier').to_dict()

db = pd.read_csv('../output/merged_luad_maf.txt')
db = db[db['Variant_Classification'] != 'Silent']
db = db[~pd.isnull(db['Mutation'])]
#for now in here, but can be cut and pasted into the importing_maf_data file
db.loc[db['Source'] == 'Genie','Panel'] = db['Sample ID'].map(genie_panels_used['Sequence Assay ID'])
db.loc[db['Source'] == 'MSK2017','Panel'] = db['Sample ID'].map(msk2017_panels_used['Gene Panel'])
db.loc[db['Source'] == 'MSK2018','Panel'] = db['Sample ID'].map(msk2018_panels_used['Gene Panel'])
db.loc[db['Source'] == 'TSP', 'Panel'] = 'TSP'
db.loc[db['Source'] == 'FM-AD', 'Panel'] = 'FoundationOne'

all_panel_genes = pd.read_csv('../gene_panels/all_panel_genes.txt')
included_panels = pd.unique(all_panel_genes['SEQ_ASSAY_ID'])

panels_to_remove_for_tp53_kras = []
for panel in included_panels:
    if 'TP53' not in all_panel_genes[all_panel_genes['SEQ_ASSAY_ID'] == panel]['Hugo_Symbol'].tolist() or 'KRAS' not in all_panel_genes[all_panel_genes['SEQ_ASSAY_ID'] == panel]['Hugo_Symbol'].tolist():
        panels_to_remove_for_tp53_kras.append(panel)
db = db[~db['Panel'].isin(panels_to_remove_for_tp53_kras)]

smoking_sample_ids = pd.read_csv('../data/smoking_sample_ids.txt', header=None).iloc[:,0].tolist()
nonsmoking_sample_ids = pd.read_csv('../data/nonsmoking_sample_ids.txt', header=None).iloc[:,0].tolist()

#DATA ONLY CONTAINING EXOME/GENOME SEQUENCED AND DETERMINED TO HAVE NO SMOKING SIGNATURE SAMPLES
#db = db[db['Sample ID'].isin(nonsmoking_sample_ids)]
#print(db.shape)

#ranked_gene_list = ['TP53', 'KRAS', 'EGFR', 'STK11', 'ATM', 'NF1', 'BRAF', 'PIK3CA', 'CDKN2A', 'ERBB4', 'APC', 'NTRK3', 'EPHA3', 'MET', 'EPHA5', 'RB1', 'ERBB2', 'KDR', 'ALK', 'SMAD4', 'GNAS', 'PDGFRA', 'PTEN', 'FLT1', 'CTNNB1', 'BRCA2', 'NOTCH1', 'NTRK1', 'ATR', 'NOTCH2', 'FLT4', 'RET', 'PIK3CG', 'NOTCH4', 'BRCA1', 'KIT', 'TSC2', 'SMO', 'NOTCH3', 'PIK3C2G', 'FLT3', 'JAK3', 'DDR2', 'EPHB1', 'PDGFRB', 'JAK2', 'NTRK2', 'BARD1', 'MAP3K1', 'FBXW7', 'WT1', 'TSHR', 'CSF1R', 'ABL1', 'FGFR4', 'ERBB3', 'CBL', 'IGF1R', 'IRS2', 'MSH6', 'FGFR1', 'RICTOR', 'FGFR2', 'NRAS', 'TSC1', 'PIK3R1', 'MSH2', 'CHEK2', 'BTK', 'PIK3C3', 'RUNX1', 'PTPN11', 'CDH1', 'BAP1', 'AKT3', 'ARAF', 'FGFR3', 'BCL6', 'AXL', 'ERG', 'RAF1', 'JAK1', 'NF2', 'SMAD2', 'TGFBR2', 'MPL', 'MLH1', 'MAP2K1', 'MAP3K13', 'MYCN', 'MDM4', 'MAP2K4', 'MEN1', 'CCNE1', 'IKBKE', 'AKT1', 'FGF3', 'CCND2', 'PIK3R2', 'ETV6', 'SYK', 'CHEK1', 'SRC', 'AKT2', 'MYC', 'GATA1', 'MDM2', 'REL', 'VHL', 'HRAS', 'MAP2K2', 'CDKN1B', 'CDK8', 'CDK6', 'PRKAR1A', 'GSK3B', 'SUFU', 'CCND1', 'FGF4', 'AURKB', 'JUN', 'CCND3', 'CDKN2C', 'CDK4', 'BCL2', 'CRKL', 'CDKN2B']
#ranked_gene_list = ['TP53', 'KRAS','EGFR', 'STK11', 'ATM', 'NF1', 'BRAF', 'PIK3CA', 'ERBB4', 'APC', 'NTRK3', 'EPHA3', 'MET', 'EPHA5', 'RB1', 'ERBB2', 'KDR', 'ALK', 'SMAD4', 'GNAS', 'PDGFRA', 'PTEN', 'FLT1', 'CTNNB1', 'BRCA2', 'NOTCH1', 'NTRK1', 'ATR', 'NOTCH2', 'FLT4', 'RET', 'PIK3CG', 'NOTCH4', 'BRCA1', 'KIT', 'TSC2', 'SMO', 'NOTCH3', 'PIK3C2G', 'FLT3', 'JAK3', 'DDR2', 'EPHB1', 'PDGFRB', 'JAK2', 'NTRK2', 'BARD1', 'MAP3K1', 'FBXW7', 'WT1', 'TSHR', 'CSF1R', 'ABL1', 'FGFR4', 'ERBB3', 'CBL', 'IGF1R', 'IRS2', 'MSH6', 'FGFR1', 'RICTOR', 'FGFR2', 'NRAS', 'TSC1', 'PIK3R1', 'MSH2', 'CHEK2', 'BTK', 'PIK3C3', 'RUNX1', 'PTPN11', 'CDH1', 'BAP1', 'AKT3', 'ARAF', 'FGFR3', 'BCL6', 'AXL', 'ERG', 'RAF1', 'JAK1', 'NF2', 'SMAD2', 'TGFBR2', 'MPL', 'MLH1', 'MAP2K1', 'MAP3K13', 'MYCN', 'MDM4', 'MAP2K4', 'MEN1', 'CCNE1', 'IKBKE', 'AKT1', 'FGF3', 'CCND2', 'PIK3R2', 'ETV6', 'SYK', 'CHEK1', 'SRC', 'AKT2', 'MYC', 'GATA1', 'MDM2', 'REL', 'VHL', 'HRAS', 'MAP2K2', 'CDKN1B', 'CDK8', 'CDK6', 'PRKAR1A', 'GSK3B', 'SUFU', 'CCND1', 'FGF4', 'AURKB', 'JUN', 'CCND3', 'CDKN2C', 'CDK4', 'BCL2', 'CRKL', 'CDKN2B']
with open('../data/genes_list.txt') as gene_list_input:
    gene_list = gene_list_input.read().split('\n')

genes_with_uncomputable_fluxes = []
for i, gene in enumerate(gene_list[2:]):
    print(f"(gene number {i}/{len(gene_list[2:])})")
    panels_to_remove = []
    for panel in included_panels:
        if gene not in all_panel_genes[all_panel_genes['SEQ_ASSAY_ID'] == panel]['Hugo_Symbol'].tolist():
            panels_to_remove.append(panel)
    subsetted_db = db[~db['Panel'].isin(panels_to_remove)]

    if not are_all_fluxes_computable(subsetted_db, mutations = ['TP53', 'KRAS', gene]):
        genes_with_uncomputable_fluxes.append(gene)

for gene in genes_with_uncomputable_fluxes:
    gene_list.remove(gene)

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

    for i, gene in enumerate(gene_list[2:]):
        print(f"Running model with third gene {gene} "
              f"(gene number {i}/{len(gene_list[2:])})")
        genes = ['TP53', 'KRAS', gene]

        panels_to_remove = []
        for panel in included_panels:
            if gene not in all_panel_genes[all_panel_genes['SEQ_ASSAY_ID'] == panel]['Hugo_Symbol'].tolist():
                panels_to_remove.append(panel)
        subsetted_db = db[~db['Panel'].isin(panels_to_remove)]
        print('Panels excluded because they did not sequence ' + gene + ':' + str(panels_to_remove)[1:-1])
        samples = compute_samples(subsetted_db, mutations=genes)

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
