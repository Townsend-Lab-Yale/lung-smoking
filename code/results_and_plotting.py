import numpy as np
import pandas as pd
import os

from itertools import chain

from filter_data import key_filtered_dbs
from filter_data import filter_samples_for_genes

from count_combinations import updated_compute_samples
from main import at_least_000_to_001_and_110_to_111
from cancer_epistasis import asymp_CI_lambdas
from cancer_epistasis import compute_CI_gamma
from cancer_epistasis import convert_lambdas_to_dict
from cancer_epistasis import estimate_lambdas
from cancer_epistasis import compute_gammas
from theory import build_S_as_array
from theory import order_pos_lambdas

from locations import full_mutation_rate_file_names

location_results = '../to_delete'

mutation_rate_dict = {
    key: pd.read_csv(full_mutation_rate_file_names[key]).set_index('gene').to_dict()['rate']
    for key in full_mutation_rate_file_names.keys()
}

## copied from cancer_epistasis
def convert_samples_to_dict(samples):
    """Convert a samples array to a dictionary.

    :type samples: list
    :param samples: Number of patients in each mutation combination
        (as returned by :func:`count_pts_per_combination`).

    :rtype: dict
    :return: Dictionary with the samples, indexed by tuples of 1's and
        0's representing whether the mutation occur fot the gene or
        not.

    """

    M = int(np.log2(len(samples))) # 2^M is the number mutation combinations

    S = build_S_as_array(M)

    results_as_dict = {tuple(x):value
                       for x, value in zip(S, samples)}

    return results_as_dict


## modified from cancer_epistasis so that it can do pathways too
def convert_mus_to_dict(genes, dataset, pathway=False):
    """Convert the dictionary with :const:`mutation_rates` to a dictionary
    indexed by the order of `genes`.

    :type genes: list
    :param genes: List with the names of the mutations.

    """
    M = len(genes)

    if not pathway:
        mus = {(i*(0,) + (1,) + (M-i-1)*(0,)):mutation_rate_dict[dataset][genes[i]]
            for i in range(M)}
    else:
        gene_indices = [i for i, included in enumerate(genes.values()) if len(included) == 1]
        mus = {(i*(0,) + (1,) + (M-i-1)*(0,)):mutation_rate_dict[dataset][list(genes.keys())[i]] if i in gene_indices
                    else sum([mutation_rate_dict[dataset][gene] for gene in list(genes.values())[i]])
                for i in range(M)}


    return mus


## new
def at_least_000_to_010_and_001_to_011(samples):
    return np.all(samples[[0,1]] > 0)

'''
TODO:
- check if necessary sample counts are fulfilled based on user input
- make bound changes for 4+ genes
'''

## new (should probably even go to cancer_epistasis)
def produce_results_for_one_gene_set(genes, dataset):
    '''
    Calculates lambdas and gammas for a single set of genes, NOT for multiple maps
    '''
    db = filter_samples_for_genes(genes, key_filtered_dbs[dataset])

    results = {}

    print("Counting samples on each mutation combination for " + ', '.join(genes) + "...")
    samples = updated_compute_samples(db, mutations=genes, print_info=True)
    results["samples"] = convert_samples_to_dict(samples)
    print("")

    print("Estimating fluxes MLE...")
    bounds = 1
    print(f"Bounds for fluxes: {bounds}")
    mle = estimate_lambdas(samples, draws=1,
                            upper_bound_prior=bounds,
                            kwargs={'return_raw':True})

    bound_changes = 0
    while mle[1].fun < -1e+20:
        ## We figure out that when this happens for TP53+KRAS+third
        ## gene models it is because our initial estimates for fluxes
        ## to the third gene are way too large (the initial bounds are
        ## [0, 1] so the initial estimate is 0.5)
        if len(genes) != 3:
            raise Exception("Lambdas incalculable at bounds=1, suitable bound changes for more or less than 3 genes are yet to be implemented")
        if bound_changes == 4:
            return "incomputable"
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
        mle = estimate_lambdas(samples, draws=1,
                                upper_bound_prior=bounds,
                            kwargs={'return_raw':True})

    if len(genes) == 3:
        S_3 = build_S_as_array(3)
        states_with_zero = [tuple(x) for x in S_3[samples == 0]]
        indices_with_zero = [order_pos_lambdas(S_3).index((x, tuple(y))) for x, y in order_pos_lambdas(S_3) if x in states_with_zero]
        for x in indices_with_zero: mle[0]['lambdas'][x] = np.nan
    else:
        print("Warning: Combinations with 0 samples cannot have the fluxes from them discarded at this point for n != 3")

    results["lambdas"] = convert_lambdas_to_dict(mle[0])
    print("...done")
    print("")

    print("Computing fluxes confidence intervals (CI)...")
    results["lambdas_cis"] = convert_lambdas_to_dict(
        asymp_CI_lambdas(mle[0]['lambdas'],
                         samples,
                         print_progress=False))
    print("...done")

    print("Importing mutation rates...")
    if 'plus' in dataset:
        results["mus"] = convert_mus_to_dict(genes, dataset[:dataset.index('_plus')])
    else:
        results["mus"] = convert_mus_to_dict(genes, dataset)
    print("...done")
    print("")

    print("Computing selection coefficients...")
    results["gammas"] = compute_gammas(results["lambdas"],
                                       results["mus"])
    results["gammas_cis"] = compute_CI_gamma(results["lambdas_cis"],
                                           results["mus"])
    print("...done")

    return results


def produce_results(gene_combinations, dataset, save=True):
    '''Takes NESTED LISTS for gene_combinations where within the list
    is each set of genes to be tested

    '''
    results = {
        '_'.join(gene_set):produce_results_for_one_gene_set(
            gene_set, dataset)
        for gene_set in gene_combinations}
    if save:
        print("")
        print("Saving results...")
        location_dataset = os.path.join(location_results, dataset)
        for gene_set, results_dict in results.items():
            if results_dict != "incomputable":
                if not os.path.isdir(location_dataset):
                    os.mkdir(location_dataset)
                location_gene_set = os.path.join(location_dataset, gene_set)
                if not os.path.isdir(location_gene_set):
                    os.mkdir(location_gene_set)

                for key, value in results_dict.items():
                    np.save(os.path.join(location_gene_set,
                                        f'{key}.npy'), value)
        print("...done")
    return results


def produce_results_for_pathway_analysis(samples, genes_and_pathways, dataset, save=True):
    '''
    genes_and_pathways is a dictionary with the keys as the genes and or pathway names
    and the values as a LIST of the gene(s) included. Singular genes should still be entered
    as a list value.
    '''
    results = {}

    print("Estimating fluxes MLE...")
    bounds = 1
    print(f"Bounds for fluxes: {bounds}")
    mle = estimate_lambdas(samples, draws=1,
                            upper_bound_prior=bounds,
                            kwargs={'return_raw':True})

    bound_changes = 0
    while mle[1].fun < -1e+20:
        if len(genes_and_pathways) != 3:
            raise Exception("Lambdas incalculable at bounds=1, suitable bound changes for more or less than are yet to be implemented")
        if bound_changes == 4:
            return "incomputable"
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
        mle = estimate_lambdas(samples, draws=1,
                                upper_bound_prior=bounds,
                            kwargs={'return_raw':True})

    if len(genes_and_pathways) == 3:
        S_3 = build_S_as_array(3)
        states_with_zero = [tuple(x) for x in S_3[samples == 0]]
        indices_with_zero = [order_pos_lambdas(S_3).index((x, tuple(y))) for x, y in order_pos_lambdas(S_3) if x in states_with_zero]
        for x in indices_with_zero: mle[0]['lambdas'][x] = np.nan
    else:
        print("Combinations with 0 samples cannot have the fluxes from them discarded at this point for n != 3")

    results["lambdas"] = convert_lambdas_to_dict(mle[0])
    print("...done")
    print("")

    print("Computing fluxes confidence intervals (CI)...")
    results["lambdas_cis"] = convert_lambdas_to_dict(
        asymp_CI_lambdas(mle[0]['lambdas'],
                         samples,
                         print_progress=False))
    print("...done")


    print("Importing mutation rates...")
    if 'plus' in dataset:
        results["mus"] = convert_mus_to_dict(genes_and_pathways, dataset[:dataset.index('_plus')], pathway=True)
    else:
        results["mus"] = convert_mus_to_dict(genes_and_pathways, dataset, pathway=True)
    print("...done")
    print("")

    print("Computing selection coefficients...")
    results["gammas"] = compute_gammas(results["lambdas"],
                                    results["mus"])
    results["gammas_cis"] = compute_CI_gamma(results["lambdas_cis"],
                                        results["mus"])
    print("...done")

    if save:
        print("")
        print("Saving results...")
        location_pathway_output = os.path.join(location_results, dataset, '_'.join(genes_and_pathways))
        if not os.path.isdir(location_pathway_output):
            os.makedirs(location_pathway_output)
        for key, value in results.items():
            np.save(os.path.join(location_pathway_output,
                                f'{key}.npy'), value)
        print("...done")

    return results

results_to_save = ["samples", "lambdas", "lambdas_cis", "mus",
                   "gammas", "gammas_cis"]



def load_results(dataset):
    results = {analysis:
                {key: np.load(os.path.join(location_results,
                                        dataset, analysis,
                                        f'{key}.npy'),
                           allow_pickle=True).item()
                for key in results_to_save}
               for analysis in [f for f in os.listdir(
                   os.path.join(location_results, dataset))
                   if not f.startswith('.')]}

    return results

'''
EGFR_pathway_genes = ['BRAF','MEK','ERK','MAPK','MAP2K1','PIK3CA','PTEN','AKT1','NFKB','RAS','RAF','RAC','TSC1','TSC2','NF1','JAK1','JAK2','SHC','SOS','GRB2','STAT','MYC','FOXO3A','IRS1','IRS2','PDK','AMPK1','STK11']
RAS_pathway_genes = pd.read_csv('~/Downloads/PID_RAS_PATHWAY.v7.5.1.tsv', sep = '\t').iloc[18,1].split(',')
PI3K_AKT_pathway_genes = pd.read_csv('~/Downloads/PID_PI3KCI_AKT_PATHWAY.v7.5.1.tsv', sep = '\t').iloc[18,1].split(',')
MTOR_pathway_genes = pd.read_csv('~/Downloads/PID_MTOR_4PATHWAY.v7.5.1.tsv', sep = '\t').iloc[18,1].split(',')

EGFR_pathway_genes = list(set(EGFR_pathway_genes + RAS_pathway_genes + PI3K_AKT_pathway_genes + MTOR_pathway_genes))
#set(EGFR_pathway_genes) - set(RAS_pathway_genes + PI3K_AKT_pathway_genes + MTOR_pathway_genes)
EGFR_pathway_genes.pop(EGFR_pathway_genes.index('KRAS'))

from locations import gene_list_file
gene_list = list(pd.read_csv(gene_list_file, header=None)[0])
gene_list = [gene.upper() for gene in gene_list]
gene_list = gene_list[:103]

EGFR_pathway_genes = [gene for gene in EGFR_pathway_genes if gene in gene_list]

mutations = {'KRAS':['KRAS'],'EGFR':['EGFR'],'EGFR_pathway':EGFR_pathway_genes}
key = 'smoking_plus'

db = filter_samples_for_genes(list(chain(*mutations.values())), key_filtered_dbs[key])
#create new column that represents if any gene in the EGFR pathway is mutated
pathway_grouped_db = db.assign(EGFR_pathway=db[EGFR_pathway_genes].sum(axis='columns').apply(lambda x: 1 if x > 1 else x))

samples = updated_compute_samples(pathway_grouped_db, mutations = list(mutations.keys()), print_info=True)

np.save(os.path.join(location_results, key,
                        '_'.join(mutations.keys()),
                        'samples.npy'),
                    convert_samples_to_dict(samples))

produce_results_for_pathway_analysis(samples, mutations, key)

key = 'nonsmoking_plus'

db = filter_samples_for_genes(list(chain(*mutations.values())), key_filtered_dbs[key])
#create new column that represents if any gene in the EGFR pathway is mutated
pathway_grouped_db = db.assign(EGFR_pathway=db[EGFR_pathway_genes].sum(axis='columns').apply(lambda x: 1 if x > 1 else x))

samples = updated_compute_samples(pathway_grouped_db, mutations = list(mutations.keys()), print_info=True)

np.save(os.path.join(location_results, key,
                        '_'.join(mutations.keys()),
                        'samples.npy'),
                    convert_samples_to_dict(samples))

produce_results_for_pathway_analysis(samples, mutations, key)
'''
'''
EGFR_pathway_genes = ['BRAF','MEK','ERK','MAPK','MAP2K1','PIK3CA','PTEN','AKT1','NFKB','RAS','RAF','RAC','TSC1','TSC2','NF1','JAK1','JAK2','SHC','SOS','GRB2','STAT','MYC','FOXO3A','IRS1','IRS2','PDK','AMPK1','STK11']
RAS_pathway_genes = pd.read_csv('~/Downloads/PID_RAS_PATHWAY.v7.5.1.tsv', sep = '\t').iloc[18,1].split(',')
PI3K_AKT_pathway_genes = pd.read_csv('~/Downloads/PID_PI3KCI_AKT_PATHWAY.v7.5.1.tsv', sep = '\t').iloc[18,1].split(',')
MTOR_pathway_genes = pd.read_csv('~/Downloads/PID_MTOR_4PATHWAY.v7.5.1.tsv', sep = '\t').iloc[18,1].split(',')

EGFR_pathway_genes = list(set(EGFR_pathway_genes + RAS_pathway_genes + PI3K_AKT_pathway_genes + MTOR_pathway_genes))
#set(EGFR_pathway_genes) - set(RAS_pathway_genes + PI3K_AKT_pathway_genes + MTOR_pathway_genes)
EGFR_pathway_genes.pop(EGFR_pathway_genes.index('KRAS'))

from locations import gene_list_file
gene_list = list(pd.read_csv(gene_list_file, header=None)[0])
gene_list = [gene.upper() for gene in gene_list]
gene_list = gene_list[:103]

EGFR_pathway_genes = [gene for gene in EGFR_pathway_genes if gene in gene_list]

for key in ['pan_data','smoking_plus','nonsmoking_plus']:
    produce_results([['KRAS','EGFR',gene] for gene in EGFR_pathway_genes if gene != 'CCNE1'], key)
'''
for key in ['pan_data','smoking_plus','nonsmoking_plus']:
    produce_results([['TP53','KRAS','STK11','KEAP1']], key)
