import os
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

from locations import location_figures

from import_results import provide_all_relevant_lambdas_and_gammas
from import_results import results_keys
from import_results import pts_per_mutation

all_lambdas, all_gammas = provide_all_relevant_lambdas_and_gammas()


default_axis_args_lambdas= {
    'lim':1,
    'label':"Flux to 'gene'",
    'tick_each':0.2,
    'tick_minor_each':0.05,
    'scale':1}


default_axis_args_gammas= {
    'lim':400,
    'label':("Selection of 'gene' "
             "(thousands)"),
    'tick_each':100,
    'tick_minor_each':50,
    'scale':10**3}


def make_tick_labels(axis_args):
    scale = axis_args['scale']

    round_to = None
    if scale == 1:
        scale_label = ''
        round_to = 1
    elif scale == 1000:
        scale_label = 'K'
    elif scale == 10*10**6:
        scale_label = 'M'
    else:
        raise Exception(f"Scale {scale} not recognized.")


    tick_labels = (
        ["0"] + [f"{round(x, round_to)}{scale_label}"
                 for x in np.arange(axis_args['tick_each'],
                                    (axis_args['lim']
                                     + axis_args['tick_minor_each']),
                                    axis_args['tick_each'])])

    return tick_labels


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
                      axis_args['lim']+axis_args['tick_minor_each'],
                      axis_args['tick_each'])*axis_args['scale'])

    if 'tick_minor_each' in axis_args.keys():
        set_ticks(
            np.arange(0,
                      axis_args['lim']+axis_args['tick_minor_each'],
                      axis_args['tick_minor_each'])*axis_args['scale'],
            minor=True)


    set_ticklabels(make_tick_labels(axis_args))


    # if 'lim' in axis_args.keys():
    #     set_lim(0, axis_args['lim'])


    if 'label' in axis_args.keys():
        set_label(axis_args['label'])

    return ax



# def plot_estimates_ordered(all_mles, all_cis, xaxis_args=None,
#                            plot_name=None):
#     mles_110_to_111 = filter_110_to_111(all_mles)
#     cis_110_to_111 = filter_110_to_111(all_cis)

#     ## Dictionaries preserve insertion order in Python 3.7+
#     mles_110_to_111 = dict(sorted(mles_110_to_111.items(),
#                                     key=lambda item: item[1]))
#     genes_ordered = list(mles_110_to_111.keys())

#     fig = plt.figure(figsize=[5.5, 20])
#     ax = fig.add_subplot(111)

#     for y_pos, gene in enumerate(genes_ordered):
#         ax.scatter(mles_110_to_111[gene], y_pos, color='Blue')
#         ax.scatter(cis_110_to_111[gene][0], y_pos, marker='|', color='Blue')
#         ax.scatter(cis_110_to_111[gene][1], y_pos, marker='|', color='Blue')

#     ax.set_yticks(range(len(genes_ordered)))
#     ax.set_yticklabels(genes_ordered)
#     ax.set_ylim(-0.5, len(genes_ordered)-0.5)
#     ax.set_ylabel("gene")

#     if xaxis_args is not None:
#        ax = prettify_axis(ax, 'x', xaxis_args)

#     plt.tight_layout()

#     if plot_name is None:
#         plot_name = "estimates_ordered.png"
#     fig.savefig(os.path.join(location_figures,
#                              plot_name))
#     plt.close('all')
#     return fig


# def plot_estimates_comparison(plot_name=None):
#     """

#     """

#     genes_ordered = list(
#         dict(sorted(filter_110_to_111(
#             fluxes_mles['pan_data']).items(),
#         key=lambda item: item[1])).keys())

#     lambdas_mles = {
#         key:filter_110_to_111(fluxes_mles[key])
#         for key in ['smoking', 'nonsmoking']}
#     lambdas_cis = {
#         key:filter_110_to_111(fluxes_cis[key])
#         for key in ['smoking', 'nonsmoking']}

#     fig = plt.figure(figsize=[5.5, 20])
#     ax = fig.add_subplot(111)

#     for y_pos, gene in enumerate(genes_ordered):
#         ax.scatter(lambdas_mles['smoking'][gene], y_pos, color='Red')
#         ax.plot(lambdas_cis['smoking'][gene],
#                 [y_pos, y_pos],
#                 color='Red',
#                 lw=0.5,
#                 alpha=0.5)
#         ax.scatter(lambdas_mles['nonsmoking'][gene], y_pos, color='Blue')
#         ax.plot(lambdas_cis['nonsmoking'][gene],
#                 [y_pos, y_pos],
#                 color='Blue',
#                 lw=0.5,
#                 alpha=0.5)

#     ax.set_yticks(range(len(genes_ordered)))
#     ax.set_yticklabels(genes_ordered)
#     ax.set_ylim(-0.5, len(genes_ordered)-0.5)
#     ax.set_ylabel("gene")


#     axis_args = axis_args_lambdas.copy()

#     axis_args['lim'] = 0.43
#     axis_args['tick_labels'] = axis_args['tick_labels'] + ["0.4"]

#     ax = prettify_axis(ax, 'x', axis_args)

#     plt.tight_layout()

#     if plot_name is None:
#         plot_name = "fluxes_comparison.png"
#     fig.savefig(os.path.join(location_figures,
#                              plot_name))
#     plt.close('all')
#     return fig



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







def top_genes(rates, top=3):
     return dict(sorted(rates.items(),
                        key=lambda item: item[1],
                        reverse=True)[:top])

top_genes_no_epi = set.union(*[set(top_genes(lambdas, top=5).keys())
                               for lambdas in [all_lambdas[key, 'no_epi', 'mles']
                                                         for key in results_keys]])

# for multi_key, lambdas in all_lambdas.items():
#     if multi_key[2] == 'mles':
#         print(multi_key[:2], set(top_genes(lambdas, 5).keys()))


def plot_lambdas_gammas(key=None, include_epistasis=False, genes=None,
                        with_cis=True, axis_args_lambdas=None,
                        axis_args_gammas=None, plot_name=None, genes_to_annotate=None):
    """Create a scatter plot with the fluxes and the scaled selection
    coefficients for each gene.

    :type key: str or NoneType
    :param key: What estimates to use. Can be one of:
        - pan_data (default)
        - smoking
        - nonsmoking

    :type include_epistasis: bool or str
    :param include_epistasis: If False (or the string 'no_epi'), do
        not consider epistastic effects. Otherwise consider epistatic
        effects in a model with 3 genes that includes TP53, KRAS and a
        third gene and plot the estimates of fluxes and selection
        coefficients from the TP53+KRAS genotype to the mutation of
        the third gene.

    :type genes: list or NoneType
    :param genes: Restrict points to only these genes.

    :type with_cis: bool
    :param with_cis: If True, include confidence intervals as lines in
        the plot. Only available when `no_epistasis` is False.

    :type plot_name: str or NoneType
    :param plot_name: Name used to save the plot. Default is
        fluxes_and_selections_key.png.

    """
    if key is None:
        key = "pan_data"

    if include_epistasis == 'no_epi' or not include_epistasis:
        epi_key = 'no_epi'
    else:
        epi_key = 'epi'

    lambdas_mles = all_lambdas[(key, epi_key, 'mles')]
    gammas_mles = all_gammas[(key, epi_key, 'mles')]

    if with_cis:
        lambdas_cis = all_lambdas[(key, epi_key, 'cis')]
        gammas_cis = all_gammas[(key, epi_key, 'cis')]

    if genes is None:
        genes = set.intersection(set(lambdas_mles.keys()),
                                 set(gammas_mles.keys()))
        print(f'For {key} and {epi_key}:')
        print(f'  genes in lambdas MLEs: {len(lambdas_mles.keys())}')
        print(f'  genes in gammas MLEs: {len(gammas_mles.keys())}')

        if with_cis:
            genes = set.intersection(genes,
                                     set(lambdas_cis.keys()),
                                     set(gammas_cis.keys()))
            print(f'  genes in lambdas CIs: {len(lambdas_cis.keys())}')
            print(f'  genes in gammas CIs: {len(gammas_cis.keys())}')

        print(f'  genes to plot (intersection of above): {len(genes)}')


    if axis_args_lambdas is None:
        axis_args_lambdas = default_axis_args_lambdas
    if axis_args_gammas is None:
        axis_args_gammas = default_axis_args_gammas


    # slope = curve_fit(lambda x, s: s*x,
    #                   list(lambdas_mles.values()),
    #                   list(gammas_mles.values()))[0][0]

    # bottom_genes = {'EPHA5', 'NTRK3', 'EPHB1', 'NTRK2', 'PIK3CG', 'NOTCH2', 'SMAD4', 'FLT1', 'PDGFRA'}
    # left_top_genes = {'RAF1', 'CTNNB1', 'CHEK2', 'NOTCH3', 'FGFR2', 'PTEN', 'PIK3C2G', 'NTRK1'}
    # left_bottom_genes = {'GNAS'}

    # genes_to_annotate = interesting_genes(
    #     lambdas_mles, gammas_mles,
    #     high_lambda=0.125,
    #     high_gamma=70*10**3,
    #     for_sure_include=set.union(bottom_genes, left_top_genes, left_bottom_genes))
    # if genes is not None:
    #     genes_to_annotate = set.union(genes_to_annotate, set(genes))

    bottom_genes = set()
    left_top_genes = set()
    left_bottom_genes = set()
    if genes_to_annotate is None:
        genes_to_annotate = set()
    elif genes_to_annotate == 'all':
        genes_to_annotate = genes

    fig = plt.figure()
    ax = fig.add_subplot(111)

    # ax.plot([0, axis_args_lambdas['lim']],
    #         [0, slope*axis_args_lambdas['lim']],
    #         "--",
    #         color="Red",
    #         alpha=0.1)

    genes_to_remove = set()
    for gene in genes:
        if gene not in lambdas_mles.keys():
            print(f"Gene {gene} not in lambdas, removing from plot")
            genes_to_remove.add(gene)
        if gene not in gammas_mles.keys():
            print(f"Gene {gene} not in gammas, removing from plot")
            genes_to_remove.add(gene)

    genes = set.difference(genes, genes_to_remove)

    for gene in genes:
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
        plot_name = f"fluxes_and_selections_{key}_{epi_key}.png"
    fig.savefig(os.path.join(location_figures,
                             plot_name))
    plt.close('all')

    return fig


def plot_all():
    # plot_estimates_ordered(selection_mles, selection_cis,
    #                        plot_name='fluxes_ordered',
    #                        xaxis_args=axis_args_lambdas)
    # for key in results_keys:
    #     plot_estimates_ordered(
    #         compute_gammas_mles(key),
    #         compute_gammas_cis(key),
    #         plot_name=f"selection_coeffs_ordered_{key}.png",
    #         xaxis_args=axis_args_gammas)


    # default_axis_args_lambdas= {
    #     'lim':1,
    #     'label':"Flux from TP53+KRAS to TP53+KRAS+'gene'",
    #     'tick_each':0.2,
    #     'tick_minor_each':0.05,
    #     'tick_labels':["0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35"]
    # }

    scatter_plots = {}

    for key in results_keys:
        for epi_key in ['no_epi'# , 'epi'
                        ]:
            print(f"Plotting {key}, {epi_key}")
            scatter_plots[(key, epi_key)] = plot_lambdas_gammas(
                key,
                epi_key,
                genes=top_genes_no_epi,
                genes_to_annotate='all')

    return scatter_plots


if __name__ == "__figures__":
    plot_all()
