import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from locations import location_data
from locations import fluxes_mles_file_name
from locations import fluxes_cis_file_name
from locations import location_figures


lambdas_mles = np.load(fluxes_mles_file_name, allow_pickle=True).item()
lambdas_cis = np.load(fluxes_cis_file_name, allow_pickle=True).item()


mut_rates_keys = ["pan-data", "exome", "smoking", "nonsmoking"]
mut_rates = {key:pd.read_csv(os.path.join(location_data,
                                          f'{key}_mutation_rates.txt'),
                        index_col=0)
             for key in mut_rates_keys}


def compute_gammas_mles(key):
    gammas_mles = {gene:{x_y:mle/float(mut_rates[key].loc[gene])
                         for x_y, mle in lambdas.items()}
                   for gene, lambdas in lambdas_mles.items()}
    return gammas_mles


def compute_gammas_cis(key):
    gammas_cis = {gene:{x_y:[ci[0]/float(mut_rates[key].loc[gene]),
                             ci[1]/float(mut_rates[key].loc[gene])]
                        for x_y, ci in lambdas.items()}
                  for gene, lambdas in lambdas_cis.items()}
    return gammas_cis


def plot_estimates_ordered(all_mles, all_cis, xaxis_args=None, plot_name=None):
    mles_110_to_111 = {gene:mles[((1, 1, 0), (1, 1, 1))]
                       for gene, mles in all_mles.items()}
    cis_110_to_111 = {gene:cis[((1, 1, 0), (1, 1, 1))]
                      for gene, cis in all_cis.items()}

    ## Dictionaries preserve insertion order in Python 3.7+
    mles_110_to_111 = dict(sorted(mles_110_to_111.items(),
                                    key=lambda item: item[1]))
    genes_ordered = list(mles_110_to_111.keys())

    fig = plt.figure(figsize=[5.5, 20])
    ax = fig.add_subplot(111)

    for y_pos, gene in enumerate(genes_ordered):
        ax.scatter(mles_110_to_111[gene], y_pos, color='Blue')
        ax.scatter(cis_110_to_111[gene][0], y_pos, marker='|', color='Blue')
        ax.scatter(cis_110_to_111[gene][1], y_pos, marker='|', color='Blue')

    ax.set_yticks(range(len(genes_ordered)))
    ax.set_yticklabels(genes_ordered)
    ax.set_ylim(-0.5, len(genes_ordered)-0.5)
    ax.set_ylabel("gene")

    if xaxis_args is None:
       xaxis_args = {}

    if 'tick_each' in xaxis_args.keys():
        ax.set_xticks(
            np.arange(0,
                      xaxis_args['lim'],
                      xaxis_args['tick_each']))

    if 'tick_minor_each' in xaxis_args.keys():
        ax.set_xticks(
            np.arange(0,
                      xaxis_args['lim']+xaxis_args['tick_minor_each'],
                      xaxis_args['tick_minor_each']),
            minor=True)

    if 'tick_labels' in xaxis_args.keys():
        ax.set_xticklabels(xaxis_args['tick_labels'])

    if 'lim' in xaxis_args.keys():
        ax.set_xlim(0, xaxis_args['lim'])


    if 'label' in xaxis_args.keys():
        ax.set_xlabel(xaxis_args['label'])


    plt.tight_layout()
    if plot_name is None:
        plot_name = "estimates_ordered.png"

    fig.savefig(os.path.join(location_figures,
                             f"{plot_name}"))
    plt.close('all')
    return cis_110_to_111


def plot_all():
    plot_estimates_ordered(lambdas_mles, lambdas_cis,
                           plot_name='fluxes_ordered',
                           xaxis_args={'lim':0.42,
                                       'label':"Flux from TP53+KRAS to TP53+KRAS+'gene'",
                                       'tick_each':0.1,
                                       'tick_minor_each':0.01})
    for key in mut_rates_keys:
        plot_estimates_ordered(
            compute_gammas_mles(key),
            compute_gammas_cis(key),
            plot_name=f"selection_coeffs_ordered_{key}.png",
            xaxis_args={'lim':270*10**3,
                        'label':("Selection from "
                                 "TP53+KRAS to TP53+KRAS+'gene' "
                                 "(thousands)"),
                        'tick_each':50*10**3,
                        'tick_minor_each':10*10**3,
                        'tick_labels':["0", "50K", "100K", "150K", "200K", "250K"]})


if __name__ == "__figures__":
    plot_all()
