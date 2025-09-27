import numpy as np
import sys

sys.path.insert(0, '../')
from theory import order_pos_lambdas
from theory import build_S_as_array

from landscape_plotting import plot_landscape

from main import ave_tmb_and_samples_per_genotype
from cancer_epistasis import convert_mus_to_dict

def convert_samples_to_dict(samples):
    M = int(np.log2(len(samples))) # 2^M is the number mutation combinations
    S = build_S_as_array(M)
    results_as_dict = {tuple(x):value
                       for x, value in zip(S, samples)}
    return results_as_dict

def plot_trajectory(gene_list, dataset_key, param="selection",
                    fluxes_mles=None,selection_mles=None,mu_dict=None,all_samples=None,all_tmbs=None,
                    scale_circle_areas=0.02,multiplier_font_size=3,scale_arrows=None,mutation_colors=None):
    if not isinstance(gene_list,list):
        raise TypeError("`gene_list` must be a list of genes")
    if param not in ["fixation","mutation","selection"]: 
        raise ValueError("`param` must be one of 'fixation','mutation','selection'")
    if scale_arrows is None:
        if param == "fixation": scale_arrows = 1
        elif param == "mutation": scale_arrows = 0.5*10**6
        elif param == "selection": scale_arrows = 0.25*10**(-6)
    
    gene_tuple = tuple(gene_list)

    if param == "fixation": values = fluxes_mles[dataset_key][gene_tuple]
    elif param == "mutation": 
        tmbs_for_mus, samples_for_mus = ave_tmb_and_samples_per_genotype(
                all_tmbs,
                combo=gene_list,
                key=dataset_key,
                source="TCGA")
        values = convert_mus_to_dict(mu_dict[dataset_key], gene_list, tmbs_for_mus, samples_for_mus)
    elif param == "selection": values = selection_mles[dataset_key][gene_tuple]

    p = plot_landscape(
        arrows = values,
        circle_areas = convert_samples_to_dict(all_samples[dataset_key][gene_tuple]),
        mutation_names = gene_tuple,
        scale_arrows=scale_arrows,
        positions="left_to_right",
        include_n_circles=True,
        scale_circle_areas=scale_circle_areas,
        multiplier_font_size=multiplier_font_size,
        mutation_colors=mutation_colors,
        plot_name='trajectory' + '_' + dataset_key + '_' + '_'.join(gene_tuple) + '_' + param
    )
    return p