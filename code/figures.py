import os
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from statsmodels.graphics.mosaicplot import mosaic

from locations import location_figures

from import_results import provide_all_relevant_lambdas_and_gammas
from import_results import results_keys
from import_results import pts_per_mutation
from import_results import samples_per_combination

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
    elif scale == 10**6:
        scale_label = 'M'
        round_to = 1
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


    if 'lim' in axis_args.keys():
        set_lim(0, axis_args['lim']*axis_args['scale'])


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

top_genes_from_110 = set.union(
    set(top_genes(all_lambdas['pan_data', 'from_110', 'mles'],
                  top=4).keys()).difference(
                      {'CSMD3', 'TTN'}),
    set(top_genes(all_lambdas['smoking', 'from_110', 'mles'],
                  top=5).keys()).difference(
                   {'CSMD3', 'TTN'}),
    set(top_genes(all_lambdas['nonsmoking', 'from_110', 'mles'],
                  top=3).keys()))

top_genes_from_normal = set.union(
    set(top_genes(all_lambdas['pan_data', 'from_normal', 'mles'],
                  top=4).keys()).difference(
                   {'CSMD3', 'TTN'}),
    set(top_genes(all_lambdas['smoking', 'from_normal', 'mles'],
                  top=5).keys()).difference(
                   {'CSMD3', 'TTN'}),
    set(top_genes(all_lambdas['nonsmoking', 'from_normal', 'mles'],
                  top=5).keys()).difference(
                   {'CSMD3', 'TTN'}))


top_genes_no_epi = set()
top_genes_from_110 = set()
top_genes_from_normal = set(pts_per_mutation[:10].index)

top_genes_epi = set.union(top_genes_from_110, top_genes_from_normal)


# for multi_key, lambdas in all_lambdas.items():
#     if multi_key[2] == 'mles':
#         print(multi_key[:2], set(top_genes(lambdas, 5).keys()))


def plot_lambdas_gammas(key=None,
                        epi_key='no_epi',
                        genes=None,
                        with_cis=True,
                        axis_args_lambdas=None,
                        axis_args_gammas=None,
                        genes_to_annotate=None,
                        bottom_genes=None,
                        left_top_genes=None,
                        left_bottom_genes=None,
                        title=None,
                        plot_name=None):
    """Create a scatter plot with the fluxes and the scaled selection
    coefficients for each gene.

    :type key: str or NoneType
    :param key: What estimates to use. Can be one of:
        - pan_data (default)
        - smoking
        - nonsmoking

    :type epi_key: str or NoneType
    :param epi_key: If 'no_epi' (default), do not
        consider epistastic effects. Otherwise consider epistatic
        effects in a model with 3 genes that includes TP53, KRAS and a
        third gene and plot the estimates of fluxes and selection
        coefficients from either the normal genotype (if `epi_key` is
        'from_normal') or the TP53+KRAS genotype (if `epi_key` is
        'from_110') to the mutation of the third gene.

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



    # bottom_genes = {'EPHA5', 'NTRK3', 'EPHB1', 'NTRK2', 'PIK3CG',
    #                 'NOTCH2', 'SMAD4', 'FLT1', 'PDGFRA'}
    # left_top_genes = {'RAF1', 'CTNNB1', 'CHEK2', 'NOTCH3', 'FGFR2',
    #                   'PTEN', 'PIK3C2G', 'NTRK1'}
    # left_bottom_genes = {'GNAS'}

    # genes_to_annotate = interesting_genes(
    #     lambdas_mles, gammas_mles,
    #     high_lambda=0.125,
    #     high_gamma=70*10**3,
    #     for_sure_include=set.union(bottom_genes, left_top_genes, left_bottom_genes))
    # if genes is not None:
    #     genes_to_annotate = set.union(genes_to_annotate, set(genes))

    if bottom_genes is None:
        bottom_genes = set()
    if left_top_genes is None:
        left_top_genes = set()
    if left_bottom_genes is None:
        left_bottom_genes = set()

    if genes_to_annotate is None:
        genes_to_annotate = set()
    elif genes_to_annotate == 'all':
        genes_to_annotate = genes

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ## This needs to make sure that lambdas and gammas are in the same order
    # genes_in_both_sets = set.union(lambdas_mles.keys(), gammas_mles.keys())
    # slope = curve_fit(lambda x, s: s*x,
    #                   list(lambdas_mles.values()),
    #                   list(gammas_mles.values()))[0][0]
    # ax.plot([0, axis_args_lambdas['lim']],
    #         [0, slope*axis_args_lambdas['lim']],
    #         "--",
    #         color="Purple",
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
        if gene in set.intersection(top_genes_epi, top_genes_no_epi):
            color = 'Purple'
        elif gene in top_genes_no_epi:
            color = 'Blue'
        elif gene in top_genes_from_normal:
            color = 'Red'
        else:
            color = 'Pink'
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
                    fontsize=(7 if len(gene) > 4 else 8),
                    zorder=3,
                    color=color)
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
                   color=color,
                   alpha=0.5,
                   zorder=5)

    if title is not None:
        ax.set_title(title)

    ax = prettify_axis(ax, 'x',
                       axis_args_lambdas)
    ax = prettify_axis(ax, 'y',
                       axis_args_gammas)

    if epi_key != 'no_epi':
        ## Add axes with epi model and emphasis

        def isemphasis(origin, destiny):
            if epi_key == "from_normal":
                test = (origin == (0, 0, 0) and
                        destiny == (0, 0, 1))

            elif epi_key == "from_110":
                test = (origin == (1, 1, 0) and
                        destiny == (1, 1, 1))

            return test

        axins = fig.add_axes([0.61, 0.59, 0.28, 0.28])
        h_dist = 3
        v_dist = 3
        circle_radius = 1.3
        nodes_locs = {(0, 0, 0):[0, 0],
                      (1, 0, 0):[h_dist, v_dist],
                      (0, 1, 0):[h_dist, 0],
                      (0, 0, 1):[h_dist, -v_dist],
                      (1, 1, 0):[2*h_dist, v_dist],
                      (1, 0, 1):[2*h_dist, 0],
                      (0, 1, 1):[2*h_dist, -v_dist],
                      (1, 1, 1):[3*h_dist, 0]}
        nodes_names = {(0, 0, 0):"normal",
                       (1, 0, 0):"TP53",
                       (0, 1, 0):"KRAS",
                       (0, 0, 1):"gene",
                       (1, 1, 0):"TP53\nKRAS",
                       (1, 0, 1):"TP53\ngene",
                       (0, 1, 1):"KRAS\ngene",
                       (1, 1, 1):"all"}

        for node_origin, loc_origin in nodes_locs.items():
            ## Draw arrows
            for node_destin, loc_destin in nodes_locs.items():
                if (np.sum(np.array(node_destin) - np.array(node_origin)) == 1
                    and np.all(np.array(node_destin) >= np.array(node_origin))):
                    dxy = np.array(loc_destin) - np.array(loc_origin)
                    dxy_norm = dxy/np.sqrt(np.sum(dxy**2))
                    origin = np.array(loc_origin) + dxy_norm*circle_radius
                    d_destin = dxy - 2*dxy_norm*circle_radius
                    emphasis = isemphasis(node_origin, node_destin)
                    axins.arrow(*list(origin),
                                *list(d_destin),
                                width=0.002 if emphasis else 0.001,
                                head_length=0.3,
                                head_width=0.2,
                                length_includes_head=True,
                                alpha=0.8,
                                color="Red" if emphasis else "Gray")

            ## Draw circles
            axins.add_patch(plt.Circle(loc_origin, circle_radius,
                                       color="LightGray", alpha=0.3))
            ## Add names
            axins.text(loc_origin[0], loc_origin[1],
                       nodes_names[node_origin], va='center', ha='center',
                       fontsize=8)

        axins.axes.xaxis.set_visible(False)
        axins.axes.yaxis.set_visible(False)


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

    ## * No epistasis include TP53 and KRAS
    ## ** Pan data

    key, epi_key = ('pan_data', 'no_epi')
    print(f"Plotting {key}, {epi_key}...")

    axis_args_lambdas= {
        'lim':1.7,
        'label':"Flux to gene",
        'tick_each':0.5,
        'tick_minor_each':0.1,
        'scale':1}

    axis_args_gammas= {
        'lim':775,
        'label':("Selection of gene "
                 "(thousands)"),
        'tick_each':100,
        'tick_minor_each':50,
        'scale':10**3}

    scatter_plots[(key, epi_key)] = plot_lambdas_gammas(
        key,
        epi_key,
        genes=set.union(top_genes_no_epi, top_genes_epi),
        genes_to_annotate=top_genes_no_epi,
        axis_args_lambdas=axis_args_lambdas,
        axis_args_gammas=axis_args_gammas,
        left_top_genes={'KEAP1', 'SPATA3'},
        bottom_genes={'RYR2'},
        plot_name = f"fluxes_and_selections_{key}_{epi_key}_with_KRAS_TP53.png")

    print("...done.")
    print("")


    ## * No epistasis include TP53 and KRAS
    ## ** Smoking

    key, epi_key = ('smoking', 'no_epi')
    print(f"Plotting {key}, {epi_key}...")

    scatter_plots[(key, epi_key)] = plot_lambdas_gammas(
        key,
        epi_key,
        genes=set.union(top_genes_no_epi, top_genes_epi),
        genes_to_annotate=top_genes_no_epi,
        axis_args_lambdas=axis_args_lambdas,
        axis_args_gammas=axis_args_gammas,
        left_top_genes={'KEAP1', 'BRAF'},
        left_bottom_genes={'SPATA3'},
        bottom_genes={'STK11', 'RYR2'},
        title="Smoking",
        plot_name = f"fluxes_and_selections_{key}_{epi_key}_with_KRAS_TP53.png")

    print("...done.")
    print("")


    ## * No epistasis include TP53 and KRAS but not EGFR
    ## ** Non-smoking

    key, epi_key = ('nonsmoking', 'no_epi')
    print(f"Plotting {key}, {epi_key}...")

    scatter_plots[(key, epi_key)] = plot_lambdas_gammas(
        key,
        epi_key,
        genes=set.union(top_genes_no_epi, top_genes_epi),
        genes_to_annotate=top_genes_no_epi,
        axis_args_lambdas=axis_args_lambdas,
        axis_args_gammas=axis_args_gammas,
        bottom_genes={'RYR2', 'KRAS', 'BRAF'},
        title="Non-smoking",
        plot_name = f"fluxes_and_selections_{key}_{epi_key}_with_KRAS_TP53_no_EGRF.png")

    print("...done.")
    print("")


    ## * No epistasis include TP53 and KRAS with EGFR
    ## ** Non-smoking

    axis_args_lambdas= {
        'lim':1.7,
        'label':"Flux to gene",
        'tick_each':0.5,
        'tick_minor_each':0.1,
        'scale':1}

    axis_args_gammas= {
        'lim':1.55,
        'label':("Selection of gene "
                 "(millions)"),
        'tick_each':0.5,
        'tick_minor_each':0.1,
        'scale':10**6}

    scatter_plots[(key, epi_key)] = plot_lambdas_gammas(
        key,
        epi_key,
        genes=set.union(top_genes_no_epi, top_genes_epi),
        genes_to_annotate=top_genes_no_epi,
        axis_args_lambdas=axis_args_lambdas,
        axis_args_gammas=axis_args_gammas,
        bottom_genes={'RYR2', 'BRAF'},
        left_top_genes={'KRAS', 'SPATA3'},
        title="Non-smoking",
        plot_name = f"fluxes_and_selections_{key}_{epi_key}_with_KRAS_TP53.png")

    print("...done.")
    print("")


    ## * No epistasis zoom in
    ## ** Pan data

    key, epi_key = ('pan_data', 'no_epi')
    print(f"Plotting {key}, {epi_key}...")

    axis_args_lambdas= {
        'lim':0.44,
        'label':"Flux to gene",
        'tick_each':0.1,
        'tick_minor_each':0.01,
        'scale':1}

    axis_args_gammas= {
        'lim':160,
        'label':("Selection of gene "
                 "(thousands)"),
        'tick_each':50,
        'tick_minor_each':10,
        'scale':10**3}

    scatter_plots[(key, epi_key)] = plot_lambdas_gammas(
        key,
        epi_key,
        genes=set.union(top_genes_no_epi, top_genes_epi).difference(
            {'TP53', 'KRAS'}),
        genes_to_annotate=top_genes_no_epi,
        axis_args_lambdas=axis_args_lambdas,
        axis_args_gammas=axis_args_gammas,
        left_top_genes={'KEAP1', 'SPATA3'})

    print("...done.")
    print("")


    ## * No epistasis zoom in
    ## ** Smoking

    key, epi_key = ('smoking', 'no_epi')
    print(f"Plotting {key}, {epi_key}...")

    scatter_plots[(key, epi_key)] = plot_lambdas_gammas(
        key,
        epi_key,
        genes=set.union(top_genes_no_epi, top_genes_epi).difference(
            {'TP53', 'KRAS'}),
        genes_to_annotate=top_genes_no_epi,
        axis_args_lambdas=axis_args_lambdas,
        axis_args_gammas=axis_args_gammas,
        left_top_genes={'BRAF', 'SPATA3'},
        title="Smoking",
        bottom_genes={'STK11'})

    print("...done.")
    print("")


    ## * No epistasis zoom in
    ## ** Non-smoking

    key, epi_key = ('nonsmoking', 'no_epi')
    print(f"Plotting {key}, {epi_key}...")

    scatter_plots[(key, epi_key)] = plot_lambdas_gammas(
        key,
        epi_key,
        genes=set.union(top_genes_no_epi, top_genes_epi).difference(
            {'TP53', 'KRAS'}),
        genes_to_annotate=top_genes_no_epi,
        axis_args_lambdas=axis_args_lambdas,
        axis_args_gammas=axis_args_gammas,
        title="Non-smoking")

    print("...done.")
    print("")


    ## * Epistasis from normal
    ## ** Pan data

    key, epi_key = ('pan_data', 'from_normal')
    print(f"Plotting {key}, {epi_key}...")

    axis_args_lambdas= {
        'lim':0.6,
        'label':"Flux from normal genotype to gene",
        'tick_each':0.1,
        'tick_minor_each':0.01,
        'scale':1}

    axis_args_gammas= {
        'lim':135,
        'label':("Selection of gene from normal genotype "
                 "(thousands)"),
        'tick_each':25,
        'tick_minor_each':5,
        'scale':10**3}

    scatter_plots[(key, epi_key)] = plot_lambdas_gammas(
        key,
        epi_key,
        genes=set.union(top_genes_no_epi, top_genes_epi),
        genes_to_annotate=set.union(top_genes_from_normal,
                                    top_genes_no_epi),
        axis_args_lambdas=axis_args_lambdas,
        axis_args_gammas=axis_args_gammas,
        bottom_genes={'USH2A'},
        left_top_genes={'KEAP1'},
        left_bottom_genes={'LRP1B'})

    print("...done.")
    print("")


    ## * Epistasis from normal
    ## ** Smoking

    key, epi_key = ('smoking', 'from_normal')
    print(f"Plotting {key}, {epi_key}...")

    scatter_plots[(key, epi_key)] = plot_lambdas_gammas(
        key,
        epi_key,
        genes=set.union(top_genes_no_epi, top_genes_epi),
        genes_to_annotate=set.union(top_genes_from_normal,
                                    top_genes_no_epi),
        axis_args_lambdas=axis_args_lambdas,
        axis_args_gammas=axis_args_gammas,
        bottom_genes={'ZFHX4'},
        left_top_genes=None,
        left_bottom_genes={'RYR2'},
        title="Smoking",)

    print("...done.")
    print("")


    ## * Epistasis from normal
    ## ** Non-smoking

    key, epi_key = ('nonsmoking', 'from_normal')
    print(f"Plotting {key}, {epi_key}...")

    scatter_plots[(key, epi_key)] = plot_lambdas_gammas(
        key,
        epi_key,
        genes=set.union(top_genes_no_epi, top_genes_epi),
        genes_to_annotate=set.union(top_genes_from_normal,
                                    top_genes_no_epi),
        axis_args_lambdas=axis_args_lambdas,
        axis_args_gammas=axis_args_gammas,
        bottom_genes={'STK11', 'USH2A'},
        left_top_genes={'KEAP1'},
        left_bottom_genes=None,
        title="Non-smoking")

    print("...done.")
    print("")


    ## * Epistasis from TP53+KRAS not annotated
    ## ** Pan data

    key, epi_key = ('pan_data', 'from_110')
    print(f"Plotting {key}, {epi_key}...")

    axis_args_lambdas= {
        'lim':2.7,
        'label':"Flux from TP53+KRAS genotype to gene",
        'tick_each':0.5,
        'tick_minor_each':0.1,
        'scale':1}

    axis_args_gammas= {
        'lim':2.4,
        'label':("Selection of gene from TP53+KRAS genotype "
                 "(millions)"),
        'tick_each':0.5,
        'tick_minor_each':0.1,
        'scale':10**6}

    scatter_plots[(key, epi_key)] = plot_lambdas_gammas(
        key,
        epi_key,
        genes=set.union(top_genes_no_epi, top_genes_epi),
        genes_to_annotate=None,
        axis_args_lambdas=axis_args_lambdas,
        axis_args_gammas=axis_args_gammas,
        bottom_genes={'USH2A'},
        left_top_genes={'KEAP1'},
        left_bottom_genes={'LRP1B'},
        plot_name = f"fluxes_and_selections_{key}_{epi_key}_not_annotated.png")

    print("...done.")
    print("")


    ## * Epistasis from TP53+KRAS not annotated
    ## ** Smoking

    key, epi_key = ('smoking', 'from_110')
    print(f"Plotting {key}, {epi_key}...")

    scatter_plots[(key, epi_key)] = plot_lambdas_gammas(
        key,
        epi_key,
        genes=set.union(top_genes_no_epi, top_genes_epi),
        genes_to_annotate=None,
        axis_args_lambdas=axis_args_lambdas,
        axis_args_gammas=axis_args_gammas,
        bottom_genes={'ZFHX4'},
        left_top_genes={'RYR2'},
        left_bottom_genes=None,
        title="Smoking",
        plot_name = f"fluxes_and_selections_{key}_{epi_key}_not_annotated.png")

    print("...done.")
    print("")


    ## * Epistasis from TP53+KRAS not annotated
    ## ** Non-smoking

    key, epi_key = ('nonsmoking', 'from_110')
    print(f"Plotting {key}, {epi_key}...")

    scatter_plots[(key, epi_key)] = plot_lambdas_gammas(
        key,
        epi_key,
        genes=set.union(top_genes_no_epi, top_genes_epi),
        genes_to_annotate=None,
        axis_args_lambdas=axis_args_lambdas,
        axis_args_gammas=axis_args_gammas,
        bottom_genes={'STK11'},
        left_top_genes={'KEAP1'},
        left_bottom_genes=None,
        title="Non-smoking",
        plot_name = f"fluxes_and_selections_{key}_{epi_key}_not_annotated.png")

    print("...done.")
    print("")

    ## * Epistasis from TP53+KRAS not annotated
    ## ** Smoking Plus Panel Data

    key, epi_key = ('smoking_plus', 'from_110')
    print(f"Plotting {key}, {epi_key}...")

    scatter_plots[(key, epi_key)] = plot_lambdas_gammas(
        key,
        epi_key,
        genes=set.union(top_genes_no_epi, top_genes_epi),
        genes_to_annotate=None,
        axis_args_lambdas=axis_args_lambdas,
        axis_args_gammas=axis_args_gammas,
        bottom_genes={'ZFHX4'},
        left_top_genes={'RYR2'},
        left_bottom_genes=None,
        title="Smoking with Panel Data",
        plot_name = f"fluxes_and_selections_{key}_{epi_key}_not_annotated.png")

    print("...done.")
    print("")


    ## * Epistasis from TP53+KRAS not annotated
    ## ** Non-smoking Plus Panel Data

    key, epi_key = ('nonsmoking_plus', 'from_110')
    print(f"Plotting {key}, {epi_key}...")

    scatter_plots[(key, epi_key)] = plot_lambdas_gammas(
        key,
        epi_key,
        genes=set.union(top_genes_no_epi, top_genes_epi),
        genes_to_annotate=None,
        axis_args_lambdas=axis_args_lambdas,
        axis_args_gammas=axis_args_gammas,
        bottom_genes={'STK11'},
        left_top_genes={'KEAP1'},
        left_bottom_genes=None,
        title="Non-smoking with Panel Data",
        plot_name = f"fluxes_and_selections_{key}_{epi_key}_not_annotated.png")

    print("...done.")
    print("")


    ## * Epistasis from TP53+KRAS
    ## ** Pan data

    axis_args_lambdas= {
        'lim':0.85,
        'label':"Flux from TP53+KRAS genotype to gene",
        'tick_each':0.1,
        'tick_minor_each':0.05,
        'scale':1}

    axis_args_gammas= {
        'lim':550,
        'label':("Selection of gene from TP53+KRAS genotype "
                 "(thousands)"),
        'tick_each':100,
        'tick_minor_each':50,
        'scale':10**3}

    key, epi_key = ('pan_data', 'from_110')
    print(f"Plotting {key}, {epi_key}...")

    scatter_plots[(key, epi_key)] = plot_lambdas_gammas(
        key,
        epi_key,
        genes=set.union(top_genes_no_epi, top_genes_epi),
        genes_to_annotate=set.union(top_genes_epi,
                                    top_genes_no_epi),
        axis_args_lambdas=axis_args_lambdas,
        axis_args_gammas=axis_args_gammas,
        bottom_genes={'RYR2', 'ZFHX4'},
        left_top_genes={'RYR1', 'USH2A'},
        left_bottom_genes={'PAPPA2', 'CSMD1'})

    print("...done.")
    print("")


    ## * Epistasis from TP53+KRAS
    ## ** Smoking

    key, epi_key = ('smoking', 'from_110')
    print(f"Plotting {key}, {epi_key}...")

    scatter_plots[(key, epi_key)] = plot_lambdas_gammas(
        key,
        epi_key,
        genes=set.union(top_genes_no_epi, top_genes_epi),
        genes_to_annotate=set.union(top_genes_epi,
                                    top_genes_no_epi),
        axis_args_lambdas=axis_args_lambdas,
        axis_args_gammas=axis_args_gammas,
        bottom_genes={'ZFHX4', 'LRP1B'},
        left_top_genes={'RYR1', 'PAPPA2', 'NAV3'},
        left_bottom_genes={'USH2A', 'RYR2'},
        title="Smoking")

    print("...done.")
    print("")


    ## * Epistasis from TP53+KRAS
    ## ** Non-smoking

    key, epi_key = ('nonsmoking', 'from_110')
    print(f"Plotting {key}, {epi_key}...")

    scatter_plots[(key, epi_key)] = plot_lambdas_gammas(
        key,
        epi_key,
        genes=set.union(top_genes_no_epi, top_genes_epi),
        genes_to_annotate=set.union(top_genes_epi,
                                    top_genes_no_epi),
        axis_args_lambdas=axis_args_lambdas,
        axis_args_gammas=axis_args_gammas,
        bottom_genes={'STK11'},
        left_top_genes=None,
        left_bottom_genes={'RYR2', 'LRP1B'},
        title="Non-smoking")

    print("...done.")
    print("")


    return scatter_plots


## * Mosaic plot of proportions of lung cancer cases

## Based on Table 2 from Kenfield SA, Wei EK, Stampfer MJ, Rosner BA,
## Colditz GA. Comparison of aspects of smoking among the four
## histological types of lung cancer. Tob Control. 2008

kenfield_table_2 = {"L-C":[11, 7, 19, 23],
                    "SqC":[7, 1, 4, 34, 76],
                    "S-C":[2, 1, 5, 60, 110],
	            "Adenocarcinoma":[85, 4, 13, 83, 185]}

kenfield_data = {(cancer_type, 'Non-smoker'):table_values[0]
                 for cancer_type, table_values in kenfield_table_2.items()}
kenfield_data.update({(cancer_type, 'Smoker'):np.sum(table_values[1:])
                      for cancer_type, table_values in kenfield_table_2.items()})

cancer_type_colors = [cm.get_cmap("tab10")(x)
                      for x in np.linspace(0, 1, 10)]

# Create a mosaic plot using statsmodels

def create_mosaic_plot(plot_name=None):

    def tile_props(key):
        cancer_type, smoking_status = key
        colors = {"L-C":cancer_type_colors[7],
                  "SqC":cancer_type_colors[4],
                  "S-C":cancer_type_colors[3],
	          "Adenocarcinoma":cancer_type_colors[0]}
        if smoking_status == "Non-smoker":
            return {'color':colors[cancer_type], 'hatch':"/", 'alpha':0.4}
        else:
            return {'color':colors[cancer_type], 'hatch':"-", 'alpha':0.95}

    fig = fig = plt.figure(figsize=[3.61*1.1, 2.35*1.1])
    ax = fig.add_subplot(111)

    mosaic(kenfield_data,
           ax=ax,
           axes_label=True,
           gap=0.015,
           properties=tile_props,
           labelizer=lambda k: '')

    ax.yaxis.set_tick_params(rotation=90,
                             labelright=True,
                             right=True,
                             labelleft=False,
                             left=False)
    for tick_label in ax.yaxis.get_ticklabels():
        tick_label.set_verticalalignment('center')

    if plot_name is None:
        plot_name = "mosaic_plot.png"
    fig.savefig(os.path.join(location_figures,
                             plot_name),
                dpi=1000)
    plt.close('all')
    return fig



# import plotly.express as px

# # Create DataFrame from the given data
# data = {
#     'SmokingStatus': ['Non-smokers', 'Smokers', 'Non-smokers', 'Smokers', 'Non-smokers', 'Smokers', 'Non-smokers', 'Smokers'],
#     'CancerType': ['Squamous Cell Carcinoma', 'Squamous Cell Carcinoma', 'Small Cell Carcinoma', 'Small Cell Carcinoma', 'Adenocarcinoma', 'Adenocarcinoma', 'Large Cell Carcinoma', 'Large Cell Carcinoma'],
#     'Counts': [7, 115, 2, 176, 85, 285, 11, 49]
# }

# df = pd.DataFrame(data)

# # Create a sunburst plot
# fig = px.sunburst(df, path=['SmokingStatus', 'CancerType'], values='Counts', title='Cancer Types by Smoking Status')
# fig.show()

if __name__ == "__figures__":
    plot_all()
