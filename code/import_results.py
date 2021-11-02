import os
import numpy as np
import pandas as pd

from locations import location_output
from locations import pts_by_mutation_file

results_keys = ["pan_data", "smoking", "nonsmoking"]


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


## * Number of patients with mutation per gene

pts_per_mutation = pd.read_csv(pts_by_mutation_file, index_col=0)
