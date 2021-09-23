import os
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

from locations import location_data
from locations import fluxes_mles_file_name
from locations import fluxes_cis_file_name
from locations import location_figures


selection_mles = np.load(fluxes_mles_file_name, allow_pickle=True).item()
selection_cis = np.load(fluxes_cis_file_name, allow_pickle=True).item()


mut_rates_keys = ["pan-data", "exome", "smoking", "nonsmoking"]
mut_rates = {key:pd.read_csv(os.path.join(location_data,
                                          f'{key}_mutation_rates.txt'),
                        index_col=0)
             for key in mut_rates_keys}


def compute_gammas_mles(key):
    """We assume that selection coefficients do not change with exposure
    to smoking. So if `key` is `smoking` or `nonsmoking` we reset to
    `pan-data`.

    """
    if key in ["smoking", "nonsmoking"]:
        key = "pan-data"
    gammas_mles = {gene:{x_y:mle/float(mut_rates[key].loc[gene])
                         for x_y, mle in lambdas.items()}
                   for gene, lambdas in selection_mles.items()}
    return gammas_mles


def compute_gammas_cis(key):
    """We assume that selection coefficients do not change with exposure
    to smoking. So if `key` is `smoking` or `nonsmoking` we reset to
    `pan-data`.

    """
    if key in ["smoking", "nonsmoking"]:
        key = "pan-data"
    gammas_cis = {gene:{x_y:[ci[0]/float(mut_rates[key].loc[gene]),
                             ci[1]/float(mut_rates[key].loc[gene])]
                        for x_y, ci in lambdas.items()}
                  for gene, lambdas in selection_cis.items()}
    return gammas_cis


def compute_lambdas_mles(key):
    if key in ["pan-data", "exome"]:
        return selection_mles
    else:
        return {gene:{x_y:mle*float(mut_rates[key].loc[gene])
                      for x_y, mle in gammas.items()}
                for gene, gammas in compute_gammas_mles(key).items()}


def compute_lambdas_cis(key):
    if key in ["pan-data", "exome"]:
        return selection_cis
    else:
        return {gene:{x_y:[ci[0]*float(mut_rates[key].loc[gene]),
                           ci[1]*float(mut_rates[key].loc[gene])]
                      for x_y, ci in gammas.items()}
                for gene, gammas in compute_gammas_cis(key).items()}


def filter_110_to_111(all_estimates):
    return {gene:estimates[((1, 1, 0), (1, 1, 1))]
            for gene, estimates in all_estimates.items()}



axis_args_lambdas= {
    'lim':0.38,
    'label':"Flux from TP53+KRAS to TP53+KRAS+'gene'",
    'tick_each':0.05,
    'tick_minor_each':0.01,
    'tick_labels':["0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35"]}

axis_args_gammas= {
    'lim':260*10**3,
    'label':("Selection from "
             "TP53+KRAS to TP53+KRAS+'gene' "
             "(thousands)"),
    'tick_each':50*10**3,
    'tick_minor_each':10*10**3,
    'tick_labels':["0", "50K", "100K", "150K", "200K", "250K"]}


def prettify_axis(ax, which, axis_args):
    """Set ticks and labels for the xaxis or yaxis of ax.

    :type which: str
    :param which: Either 'x' or 'y' to modify the xaxis o yaxis,
        respectively.

    """

    if which == 'x':
        set_ticks = ax.set_xticks
        set_ticklabels = ax.set_xticklabels
        set_lim = ax.set_xlim
        set_label = ax.set_xlabel
    elif which == 'y':
        set_ticks = ax.set_yticks
        set_ticklabels = ax.set_yticklabels
        set_lim = ax.set_ylim
        set_label = ax.set_ylabel
    else:
        raise Exception("`which` has to be either 'x' or 'y'")


    if 'tick_each' in axis_args.keys():
        set_ticks(
            np.arange(0,
                      axis_args['lim'],
                      axis_args['tick_each']))

    if 'tick_minor_each' in axis_args.keys():
        set_ticks(
            np.arange(0,
                      axis_args['lim']+axis_args['tick_minor_each'],
                      axis_args['tick_minor_each']),
            minor=True)


    if 'tick_labels' in axis_args.keys():
        set_ticklabels(axis_args['tick_labels'])


    if 'lim' in axis_args.keys():
        set_lim(0, axis_args['lim'])


    if 'label' in axis_args.keys():
        set_label(axis_args['label'])

    return ax



def plot_estimates_ordered(all_mles, all_cis, xaxis_args=None,
                           plot_name=None):
    mles_110_to_111 = filter_110_to_111(all_mles)
    cis_110_to_111 = filter_110_to_111(all_cis)

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

    if xaxis_args is not None:
       ax = prettify_axis(ax, 'x', xaxis_args)

    plt.tight_layout()

    if plot_name is None:
        plot_name = "estimates_ordered.png"
    fig.savefig(os.path.join(location_figures,
                             plot_name))
    plt.close('all')
    return fig


def plot_estimates_comparison(plot_name=None):
    """

    """

    genes_ordered = list(
        dict(sorted(filter_110_to_111(
            compute_lambdas_mles('pan-data')).items(),
        key=lambda item: item[1])).keys())

    lambdas_mles = {
        key:filter_110_to_111(compute_lambdas_mles(key))
        for key in ['smoking', 'nonsmoking']}
    lambdas_cis = {
        key:filter_110_to_111(compute_lambdas_cis(key))
        for key in ['smoking', 'nonsmoking']}

    fig = plt.figure(figsize=[5.5, 20])
    ax = fig.add_subplot(111)

    for y_pos, gene in enumerate(genes_ordered):
        ax.scatter(lambdas_mles['smoking'][gene], y_pos, color='Red')
        ax.plot(lambdas_cis['smoking'][gene],
                [y_pos, y_pos],
                color='Red',
                lw=0.5,
                alpha=0.5)
        ax.scatter(lambdas_mles['nonsmoking'][gene], y_pos, color='Blue')
        ax.plot(lambdas_cis['nonsmoking'][gene],
                [y_pos, y_pos],
                color='Blue',
                lw=0.5,
                alpha=0.5)

    ax.set_yticks(range(len(genes_ordered)))
    ax.set_yticklabels(genes_ordered)
    ax.set_ylim(-0.5, len(genes_ordered)-0.5)
    ax.set_ylabel("gene")


    axis_args = axis_args_lambdas.copy()

    axis_args['lim'] = 0.43
    axis_args['tick_labels'] = axis_args['tick_labels'] + ["0.4"]

    ax = prettify_axis(ax, 'x', axis_args)

    plt.tight_layout()

    if plot_name is None:
        plot_name = "fluxes_comparison.png"
    fig.savefig(os.path.join(location_figures,
                             plot_name))
    plt.close('all')
    return fig



def interesting_genes(lambdas,
                      gammas,
                      high_lambda=0.15,
                      low_lambda=0.1,
                      high_gamma=100*10**3,
                      low_gamma=50*10**3,
                      for_sure_include=None):

    high_lambda_genes = {gene for gene in lambdas.keys()
                         if lambdas[gene] >= high_lambda}
    high_gamma_genes = {gene for gene in gammas.keys()
                         if gammas[gene] >= high_gamma}

    if for_sure_include is None:
        for_sure_include = set()

    return set.union(high_lambda_genes, high_gamma_genes, for_sure_include)



def plot_lambdas_gammas(key, with_cis=True, plot_name=None):

    lambdas_mles = filter_110_to_111(compute_lambdas_mles(key))
    gammas_mles = filter_110_to_111(compute_gammas_mles(key))

    ## sanity check to make sure keys are in the same order
    if list(lambdas_mles.keys()) != list(gammas_mles.keys()):
        raise Exception("lambdas and gammas keys "
                        "are not in the same order")
    slope = curve_fit(lambda x, s: s*x,
                      list(lambdas_mles.values()),
                      list(gammas_mles.values()))[0][0]

    if with_cis:
        lambdas_cis = filter_110_to_111(compute_lambdas_cis(key))
        gammas_cis = filter_110_to_111(compute_gammas_cis(key))

    bottom_genes = {'EPHA5', 'NTRK3', 'EPHB1', 'NTRK2', 'PIK3CG', 'NOTCH2', 'SMAD4', 'FLT1', 'PDGFRA'}
    left_top_genes = {'RAF1', 'CTNNB1', 'CHEK2', 'NOTCH3', 'FGFR2', 'PTEN', 'PIK3C2G', 'NTRK1'}
    left_bottom_genes = {'GNAS'}

    genes_to_annotate = interesting_genes(
        lambdas_mles, gammas_mles,
        high_lambda=0.125,
        high_gamma=70*10**3,
        for_sure_include=set.union(bottom_genes, left_top_genes, left_bottom_genes))


    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot([0, axis_args_lambdas['lim']],
            [0, slope*axis_args_lambdas['lim']],
            "--",
            color="Red",
            alpha=0.1)

    for gene in lambdas_mles.keys():
        if gene in genes_to_annotate:
            ax.text((lambdas_mles[gene]
                     + (-0.0005 if gene in set.union(left_top_genes, left_bottom_genes)
                        else 0.002)),
                    (gammas_mles[gene]
                     + (-2.5*10**3 if gene in set.union(bottom_genes, left_bottom_genes)
                        else 900)),
                    gene,
                    ha=('right' if gene in set.union(left_top_genes, left_bottom_genes)
                        else 'left'),
                    va=('top' if gene in set.union(bottom_genes, left_bottom_genes)
                        else 'bottom'),
                    fontsize=(7 if len(gene) > 5 else 8),
                    zorder=3)
        if with_cis:
            ax.plot([lambdas_mles[gene], lambdas_mles[gene]],
                    gammas_cis[gene],
                    color='Grey',
                    lw=0.5,
                    alpha=0.2,
                    zorder=4)
            ax.plot(lambdas_cis[gene],
                    [gammas_mles[gene], gammas_mles[gene]],
                    color='Grey',
                    lw=0.5,
                    alpha=0.2,
                    zorder=4)
        ax.scatter(lambdas_mles[gene],
                   gammas_mles[gene],
                   marker='.',
                   color='Blue',
                   alpha=0.5,
                   zorder=5)

    ax = prettify_axis(ax, 'x',
                       axis_args_lambdas)
    ax = prettify_axis(ax, 'y',
                       axis_args_gammas)

    if plot_name is None:
        plot_name = f"fluxes_and_selections_{key}.png"
    fig.savefig(os.path.join(location_figures,
                             plot_name))
    plt.close('all')

    return fig


def plot_all():
    plot_estimates_ordered(selection_mles, selection_cis,
                           plot_name='fluxes_ordered',
                           xaxis_args=axis_args_lambdas)
    for key in mut_rates_keys:
        plot_estimates_ordered(
            compute_gammas_mles(key),
            compute_gammas_cis(key),
            plot_name=f"selection_coeffs_ordered_{key}.png",
            xaxis_args=axis_args_gammas)

    for key in mut_rates_keys:
        plot_lambdas_gammas(key)



if __name__ == "__figures__":
    plot_all()
