import pandas as pd
import numpy as np
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
from locations import merged_maf_file_name
from locations import location_output

import os

#CHOOSE A DATASET FOR LAMBDA CALCULATIONS, BY DEFAULT IT IS PAN-DATA
possible_sample_options = ['pan_data','smoking','nonsmoking']
samples_used = possible_sample_options[0]

## We fix M=3 so
numbers_positive_lambdas = numbers_positive_lambdas[3]

genie_panels_used = pd.read_csv('../gene_panels/genie_panels_used.txt').set_index('Sample Identifier').to_dict()
msk2017_panels_used = pd.read_csv('../gene_panels/msk2017_panels_used.txt').set_index('Sample Identifier').to_dict()
msk2018_panels_used = pd.read_csv('../gene_panels/msk2018_panels_used.txt').set_index('Sample Identifier').to_dict()

db = pd.read_csv(merged_maf_file_name)
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

if samples_used == 'smoking':
    db = db[db['Sample ID'].isin(smoking_sample_ids)]
elif samples_used == 'nonsmoking':
    db = db[db['Sample ID'].isin(nonsmoking_sample_ids)]

with open('../data/genes_list.txt','r') as gene_list_input:
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
        np.save(os.path.join(location_output, str(samples_used + '_fluxes_mles.npy')), lambdas_mles)

        print("Estimating asymptotic confidence intervals...")
        cis = asymp_CI_lambdas(mle['lambdas'], samples)
        lambdas_cis[gene] = convert_lambdas_to_dict(cis)
        np.save(os.path.join(location_output, str(samples_used + '_fluxes_cis.npy')), lambdas_cis)

        print("")

    return lambdas_mles, lambdas_cis


if __name__ == "__main__":
    main()
