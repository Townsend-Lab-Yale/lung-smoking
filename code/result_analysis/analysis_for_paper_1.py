#!/usr/bin/env python
# coding: utf-8

model_results_extension = "model_results"
mu_method = "variant"

# Write data in numpy files to CSV for analysis in R

from numpy_to_csv import convert_all_numpy_to_csv
from run_epistasis_testing_pipeline import run_epistasis_testing
convert_all_numpy_to_csv(extension=model_results_extension, mu_method=mu_method)
run_epistasis_testing(method=mu_method)

# Plot evolutionary trajectories

import os
import sys

sys.path.insert(0, '../')
from load_results import load_results
from plot_trajectories import plot_trajectory
from landscape_plotting import color_blue, color_orange, color_green, color_purple, color_gray, color_red, color_brown, color_light_purple, color_chartreuse, color_light_blue, color_yellow
from main import compute_tmb

subset_extension = os.path.join(model_results_extension, "subset")

fluxes_mles = load_results('fluxes', 'mles',subset_extension)
fluxes_cis = load_results('fluxes','cis',subset_extension)

selection_mles = load_results('selections', 'mles',subset_extension)
selection_cis = load_results('selections', 'cis',subset_extension)

if mu_method == "cesR":
    mu_dict = load_results('mutations','cesR',extension=model_results_extension)
elif mu_method == "variant":
    mu_dict = load_results('mutations','variant',extension=model_results_extension)

all_samples = load_results('samples', extension=subset_extension)
all_tmbs = compute_tmb(source="TCGA")

for dataset_key in ["smoking_plus","nonsmoking_plus"]:
    for param in ["fixation","mutation","selection"]:
        plot_trajectory(['TP53','RB1'],dataset_key=dataset_key,param=param,
                        fluxes_mles=fluxes_mles, selection_mles=selection_mles, mu_dict=mu_dict, all_samples=all_samples, all_tmbs=all_tmbs,
                        mutation_colors=[color_blue,color_brown])

for dataset_key in ["smoking_plus","nonsmoking_plus"]:
    for param in ["fixation","mutation","selection"]:
        plot_trajectory(['TP53','RB1'],dataset_key=dataset_key,param=param,
                        fluxes_mles=fluxes_mles, selection_mles=selection_mles, mu_dict=mu_dict, all_samples=all_samples, all_tmbs=all_tmbs,
                        mutation_colors=[color_blue,color_brown])

for param in ["fixation","mutation","selection"]:
    plot_trajectory(['KRAS','KEAP1','STK11'],dataset_key="smoking_plus",param=param,
                    fluxes_mles=fluxes_mles, selection_mles=selection_mles, mu_dict=mu_dict, all_samples=all_samples, all_tmbs=all_tmbs,
                    mutation_colors=[color_orange,color_chartreuse,color_light_blue])

plot_trajectory(['STK11','ATM','ALK'],"smoking_plus",param='selection',
                selection_mles=selection_mles,all_samples=all_samples,
                mutation_colors=[color_light_blue,color_gray,color_red])

plot_trajectory(['TP53','KRAS','ARID1A'],"smoking_plus",param='selection',
                selection_mles=selection_mles,all_samples=all_samples,
                mutation_colors=[color_blue,color_orange,color_purple])

plot_trajectory(['TP53','KRAS','EGFR'],"smoking_plus",param='selection',
                selection_mles=selection_mles,all_samples=all_samples,
                mutation_colors=[color_blue,color_orange,color_green])

plot_trajectory(['KEAP1','STK11','APC'],"smoking_plus",param='selection',
                selection_mles=selection_mles,all_samples=all_samples,
                mutation_colors=[color_chartreuse,color_light_blue,color_light_purple])

# Plot matrix of pairwise epistatic effects

import sys
sys.path.insert(0, '../')
from matrix_plotting import plot_epistatic_ratios_separate
from order_results import order_genes_by_result_values

gene_list_by_selection = order_genes_by_result_values(selection_mles['pan_data'])

plot_epistatic_ratios_separate(selection_mles['nonsmoking_plus'],
                                 selection_cis['nonsmoking_plus'],
                                 selection_mles['smoking_plus'],
                                 selection_cis['smoking_plus'],
                                 gene_list_by_selection,
                                 plot_name="epistatic_ratios_matrix_separate")
