import os
import numpy as np
import pandas as pd

from locations import location_output
from locations import pts_by_mutation_file
from locations import results_keys
from locations import samples_per_combination_files

## * Load mutation rates

mutation_rates = {
    key:pd.read_csv(os.path.join(location_output,
                                 f'{key}_mutation_rates.txt'),
                    index_col=0)['rate'].to_dict()
    for key in results_keys}


## * Load selection rates, without epistasis

selection_ne_mles = {
    key:pd.read_csv(
        os.path.join(location_output,
                     f'{key}_selections_no_epistasis.txt'),
        index_col=0)['selection_intensity'].to_dict()
    for key in results_keys}


selection_ne_cis = {
    key:pd.read_csv(
        os.path.join(location_output,
                     f'{key}_selections_no_epistasis.txt'),
        index_col=0)[['ci_low_95', "ci_high_95"]].apply(
            lambda x: [x[0], x[1]], axis=1).to_dict()
    for key in results_keys}


## * Compute fluxes, without epistasis

def compute_lambdas(gammas, mu):
    if isinstance(gammas, float):
        return gammas*mu
    elif isinstance(gammas, list):
        return [gammas[0]*mu, gammas[1]*mu]
    elif isinstance(gammas, dict):
        return {x_y:compute_gammas(the_gamma, mu)
                for x_y, the_gamma in gammas.items()}

fluxes_ne_mles = {
    key:{gene:compute_lambdas(selection_ne_mles[key][gene],
                              mutation_rates[key][gene])
         for gene in selection_ne_mles[key].keys()}
    for key in results_keys}


fluxes_ne_cis = {
    key:{gene:compute_lambdas(selection_ne_cis[key][gene],
                              mutation_rates[key][gene])
         for gene in selection_ne_cis[key].keys()}
    for key in results_keys}


## * Load fluxes with epistasis

fluxes_mles = {
    key:np.load(os.path.join(location_output,
                             f'{key}_fluxes_mles.npy'),
                allow_pickle=True).item()
    for key in results_keys}


fluxes_cis = {
    key:np.load(os.path.join(location_output,
                             f'{key}_fluxes_cis.npy'),
                allow_pickle=True).item()
    for key in results_keys}


## * Compute selection coefficients with epistasis

def compute_gammas(lambdas, mu):
    if isinstance(lambdas, float):
        return lambdas/mu
    elif isinstance(lambdas, list):
        return [lambdas[0]/mu, lambdas[1]/mu]
    elif isinstance(lambdas, dict):
        return {x_y:compute_gammas(the_lambda, mu)
                for x_y, the_lambda in lambdas.items()}


selection_mles = {
    key:{gene:compute_gammas(fluxes_mles[key][gene],
                             mutation_rates[key][gene])
         for gene in set.intersection(
                 set(fluxes_mles[key].keys()),
                 set([gene.upper() for gene in mutation_rates[key].keys()]))}
    for key in results_keys}


selection_cis = {
    key:{gene:compute_gammas(fluxes_cis[key][gene],
                             mutation_rates[key][gene])
         for gene in set.intersection(
                 set(fluxes_cis[key].keys()),
                 set([gene.upper() for gene in mutation_rates[key].keys()]))}
    for key in results_keys}


## * Helper function to filter the results with epistasis

def filter_estimates(all_estimates, from_x_to_y, genes=None):
    if genes is None:
        genes = list(all_estimates.keys())
    return {gene:estimates[from_x_to_y]
            for gene, estimates in all_estimates.items()
            if gene in genes}


def filter_110_to_111(all_estimates, genes=None):
    return filter_estimates(all_estimates,
                            ((1, 1, 0), (1, 1, 1)),
                            genes)

def filter_000_to_001(all_estimates, genes=None):
    return filter_estimates(all_estimates,
                            ((0, 0, 0), (0, 0, 1)),
                            genes)



def provide_all_relevant_lambdas_and_gammas():
    """Construct a dictionary with all relevant results for fluxes and
    selection coefficients.

    Keys of the dictionary are tuples of the form:

         (result_key, es, what)

    result_key can be any of:
        - pan_data
        - smoking
        - nonsmoking

    es is the epistasis status:
        - 'no_epi': for no epistasis considered
        - 'epi': if epistasis is considered (in this case we consider
          fluxes and selections from KRAS+TP53 to KRAS+TP53+ the third
          gene in the model)

    what refers to the estimation:
        - 'mles': for the maximum likehood estimator
        - 'cis': for the 95% confidence interval (given as a two item list)

    This function returns a tuple with the lambdas and the
    gammas. Each value of lambdas and gammas is another dictionary
    with the third gene as key and respective estimate as value.

    """

    lambdas = {(key, 'no_epi', 'mles'):fluxes_ne_mles[key]
               for key in results_keys}
    lambdas.update({(key, 'no_epi', 'cis'):fluxes_ne_cis[key]
               for key in results_keys})

    lambdas.update({(key, 'from_110', 'mles'):filter_110_to_111(fluxes_mles[key])
                    for key in results_keys})
    lambdas.update({(key, 'from_110', 'cis'):filter_110_to_111(fluxes_cis[key])
                    for key in results_keys})
    lambdas.update({(key, 'from_normal', 'mles'):filter_000_to_001(fluxes_mles[key])
                    for key in results_keys})
    lambdas.update({(key, 'from_normal', 'cis'):filter_000_to_001(fluxes_cis[key])
                    for key in results_keys})

    gammas = {(key, 'no_epi', 'mles'):selection_ne_mles[key]
              for key in results_keys}
    gammas.update({(key, 'no_epi', 'cis'):selection_ne_cis[key]
               for key in results_keys})

    gammas.update({(key, 'from_110', 'mles'):filter_110_to_111(selection_mles[key])
                    for key in results_keys})
    gammas.update({(key, 'from_110', 'cis'):filter_110_to_111(selection_cis[key])
                    for key in results_keys})
    gammas.update({(key, 'from_normal', 'mles'):filter_000_to_001(selection_mles[key])
                    for key in results_keys})
    gammas.update({(key, 'from_normal', 'cis'):filter_000_to_001(selection_cis[key])
                    for key in results_keys})

    return lambdas, gammas



## * Number of patients with mutation per gene

pts_per_mutation = pd.read_csv(pts_by_mutation_file, index_col=0)


## * Patients per mutation combination for all TP53, KRAS, and third gene models

samples_per_combination = {
    key:pd.read_csv(samples_per_combination_files[key],
                    index_col='third gene')
    for key in results_keys}
