import os
import pandas as pd
import numpy as np
from count_combinations import compute_samples
from count_combinations import are_all_fluxes_computable

from cancer_epistasis import estimate_lambdas
from cancer_epistasis import asymp_CI_lambdas
from cancer_epistasis import convert_lambdas_to_dict
# from figures import plot_lambdas_gammas

from locations import location_data
from locations import location_gene_panels
from locations import all_panel_genes_file_name
from locations import gene_list_file
from locations import merged_maf_file_name
from locations import location_output
from locations import results_keys

from filter_data import prefiltered_dbs
from filter_data import filter_db_for_gene

#CHOOSE A DATASET FOR LAMBDA CALCULATIONS, BY DEFAULT IT IS PAN-DATA
samples_used = results_keys[0]

db = prefiltered_dbs[samples_used]

gene_list = list(pd.read_csv(gene_list_file, header=None)[0])

genes_with_uncomputable_fluxes = []
for i, gene in enumerate(gene_list[2:]):
    print(f"(gene number {i}/{len(gene_list[2:])})")
    panels_to_remove = []
    for panel in included_panels:
        if gene not in all_panel_genes[
                all_panel_genes['SEQ_ASSAY_ID'] == panel]['Hugo_Symbol'].tolist():
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


counts = {key:{} for key in results_keys}

def main(save_results=True):
    lambdas_mles = {}
    lambdas_cis = {}

    for i, gene in enumerate(gene_list[2:]):
        print(f"Running model with third gene {gene} "
              f"(gene number {i}/{len(gene_list[2:])})")

        subsetted_db = filter_db_for_gene(gene, db, print_info=True)

        samples = compute_samples(subsetted_db, mutations=['TP53', 'KRAS', gene])

        print("Estimating fluxes...")
        mle = lambdas_from_samples(samples)
        lambdas_mles[gene] = convert_lambdas_to_dict(mle)
        if save_results:
            np.save(os.path.join(location_output, str(samples_used + '_fluxes_mles.npy')),
                    lambdas_mles)

        print("Estimating asymptotic confidence intervals...")
        cis = asymp_CI_lambdas(mle['lambdas'], samples)
        lambdas_cis[gene] = convert_lambdas_to_dict(cis)
        if save_results:
            np.save(os.path.join(location_output, str(samples_used + '_fluxes_cis.npy')),
                    lambdas_cis)

        print("")

    return lambdas_mles, lambdas_cis


if __name__ == "__main__":
    main()
