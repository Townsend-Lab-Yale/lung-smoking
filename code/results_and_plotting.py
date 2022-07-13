import numpy as np
import pandas as pd
import os

from filter_data import key_filtered_dbs
from filter_data import filter_samples_for_genes

from count_combinations import updated_compute_samples
from main import are_all_fluxes_computable
from cancer_epistasis import asymp_CI_lambdas, compute_CI_gamma, convert_lambdas_to_dict, estimate_lambdas, compute_gammas
from plotting import plot_all
from theory import build_S_as_array
from locations import full_mutation_rate_file_names, results_keys

location_results = '../to_delete'

mutation_rate_dict = {
    key: pd.read_csv(full_mutation_rate_file_names[key]).set_index('gene').to_dict()['rate']
    for key in full_mutation_rate_file_names.keys()
}

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

def convert_mus_to_dict(genes, dataset):
    """Convert the dictionary with :const:`mutation_rates` to a dictionary
    indexed by the order of `genes`.

    :type genes: list
    :param genes: List with the names of the mutations.

    """
    M = len(genes)
    mus = {(i*(0,) + (1,) + (M-i-1)*(0,)):mutation_rate_dict[dataset][genes[i]]
           for i in range(M)}
    return mus

def at_least_000_to_010_and_001_to_011(samples):
    return np.all(samples[[0,1]] > 0)

def produce_results_for_one_gene_set(genes, dataset):
    '''
    Calculates lambdas and gammas for a single set of genes, NOT for multiple maps
    '''

    db = filter_samples_for_genes(genes, key_filtered_dbs[dataset])

    results = {}

    print("Counting samples on each mutation combination for " + genes + "...")
    samples = updated_compute_samples(db, mutations=genes)
    results["samples"] = convert_samples_to_dict(samples)
    print("...done")
    print("")

    print("Estimating fluxes MLE...")
    mle = estimate_lambdas(samples, draws=1)
    results["lambdas"] = convert_lambdas_to_dict(mle)
    print("...done")
    print("")

    print("Computing fluxes confidence intervals (CI)...")
    results["lambdas_cis"] = convert_lambdas_to_dict(
        asymp_CI_lambdas(mle['lambdas'],
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
    results = {gene_set: produce_results_for_one_gene_set(gene_set, dataset) for gene_set in gene_combinations}
    if save:
        print("")
        print("Saving results...")
        location_dataset = os.path.join(location_results, dataset)
        if not os.path.isdir(location_dataset):
            os.mkdir(location_dataset)
        for gene_set, results_dict in results.items():
            location_gene_set = os.path.join(location_results, '_'.join(gene_set))
            if not os.path.isdir(location_gene_set):
                os.mkdir(location_gene_set)

            for key, value in results_dict.items():
                np.save(os.path.join(location_gene_set,
                                    f'{key}.npy'), value)
        print("...done")


results_to_save = ["samples", "lambdas", "lambdas_cis", "mus",
                   "gammas", "gammas_cis"]



def load_results(dataset):
    results = {key:np.load(os.path.join(location_results,
                                        f'{dataset}_{key}.npy'),
                           allow_pickle=True).item()
               for key in results_to_save}

    return results

'''
genes = ['TP53','KRAS','EGFR']

#produce_results(genes, 'pan_data')
#produce_results(genes, 'smoking_plus')
produce_results(genes, 'nonsmoking_plus')
'''
