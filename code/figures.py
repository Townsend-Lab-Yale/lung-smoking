import os
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from statsmodels.graphics.mosaicplot import mosaic

from locations import location_figures
from locations import location_output

# from import_results import provide_all_relevant_lambdas_and_gammas
# from import_results import results_keys
# from import_results import pts_per_mutation
# from import_results import samples_per_combination
# from import_results import mutation_rates


from load_results import load_results

from landscape_plotting import plot_landscape

all_samples = load_results('samples')
all_lambdas = load_results('fluxes')
all_mus = load_results('mutations')
all_gammas = load_results('selections')


# all_lambdas, all_gammas = provide_all_relevant_lambdas_and_gammas(["smoking_plus", "nonsmoking_plus"])

# results_TP53_KRAS_smoking = np.load(
#     os.path.join(location_output,
#                  "results_TP53_KRAS_model_smoking_plus.npy"),
#     allow_pickle=True).item()

# results_TP53_KRAS_nonsmoking = np.load(
#     os.path.join(location_output,
#                  "results_TP53_KRAS_model_nonsmoking_plus.npy"),
#     allow_pickle=True).item()


# results_TP53_KRAS_EGFR_smoking = np.load(
#     os.path.join(location_output,
#                  "results_TP53_KRAS_EGFR_model_smoking_plus.npy"),
#     allow_pickle=True).item()

# results_TP53_KRAS_EGFR_nonsmoking = np.load(
#     os.path.join(location_output,
#                  "results_TP53_KRAS_EGFR_model_nonsmoking_plus.npy"),
#     allow_pickle=True).item()



from cancer_epistasis import convert_samples_to_dict
from cancer_epistasis import convert_mus_to_dict


def plot_figures_presentation_montreal():
    all_plots = {}

    ## * Epistasis TP53, KRAS, EGFR in smokers
    print(f"Plotting TP53, KRAS, EGFR model")

    all_plots['tp53_kras_egfr_landscape_smoking', 'lambdas'] = plot_landscape(
        all_lambdas['smoking_plus'][('TP53', 'KRAS', 'EGFR')],
        all_samples['smoking_plus'][('TP53', 'KRAS', 'EGFR')],
        mutation_names=['TP53', 'KRAS', 'EGFR'],
        positions='left_to_right',
        scale_arrows=1,
        scale_circle_areas=0.027,
        include_n_circles=False,
        subplot_label="Substitution /",
        subplot_label_ha='right',
        subplot_label_va='top',
        multiplier_figsize=0.4,
        multiplier_font_size=0.4,
        plot_name='landscape_lambdas_tp53_kras_egfr_smoking')


    all_plots['tp53_kras_egfr_landscape_smoking', 'mus'] = plot_landscape(
        convert_mus_to_dict(all_mus['smoking_plus'], ['TP53', 'KRAS', 'EGFR']),
        all_samples['smoking_plus'][('TP53', 'KRAS', 'EGFR')],
        mutation_names=['TP53', 'KRAS', 'EGFR'],
        positions='left_to_right',
        scale_arrows=0.25*10**(6),
        scale_circle_areas=0.027,
        include_n_circles=False,
        subplot_label="Mutation rate",
        subplot_label_ha='right',
        subplot_label_va='top',
        multiplier_figsize=0.4,
        multiplier_font_size=0.4,
        plot_name='landscape_mus_tp53_kras_egfr_smoking')


    all_plots['tp53_kras_egfr_landscape_smoking', 'gammas'] = plot_landscape(
        all_gammas['smoking_plus'][('TP53', 'KRAS', 'EGFR')],
        all_samples['smoking_plus'][('TP53', 'KRAS', 'EGFR')],
        mutation_names=['TP53', 'KRAS', 'EGFR'],
        positions='left_to_right',
        scale_arrows=0.25*10**(-6),
        scale_circle_areas=0.027,
        include_n_circles=False,
        subplot_label="Selection =",
        subplot_label_ha='left',
        subplot_label_va='top',
        multiplier_figsize=0.4,
        multiplier_font_size=0.4,
        plot_name='landscape_gammas_tp53_kras_egfr_smoking')


    all_plots['braf_ctnnb1_setd2_landscape_smoking', 'gammas'] = plot_landscape(
        all_gammas['smoking_plus'][('BRAF', 'CTNNB1', 'SETD2')],
        all_samples['smoking_plus'][('BRAF', 'CTNNB1', 'SETD2')],
        mutation_names=['BRAF', 'CTNNB1', 'SETD2'],
        positions='left_to_right',
        scale_arrows=0.25*10**(-6),
        scale_circle_areas=0.0165,
        include_n_circles=False,
        # subplot_label="Selection =",
        # subplot_label_ha='left',
        # subplot_label_va='top',
        multiplier_figsize=0.75,
        multiplier_font_size=0.9,
        plot_name='landscape_gammas_braf_ctnnb1_setd2_smoking')

    return all_plots














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
        round_to = 2
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


## TODO: Needs to fix bug on import_figures, where keys are not matching for no_epi
# top_genes_no_epi = set.union(*[set(top_genes(lambdas, top=5).keys())
#                                for lambdas in [all_lambdas[key, 'no_epi', 'mles']
#                                                          for key in results_keys]])

# top_genes_from_110 = set.union(
#     set(top_genes(all_lambdas['pan_data', 'from_110', 'mles'],
#                   top=4).keys()).difference(
#                       {'CSMD3', 'TTN'}),
#     set(top_genes(all_lambdas['smoking', 'from_110', 'mles'],
#                   top=5).keys()).difference(
#                    {'CSMD3', 'TTN'}),
#     set(top_genes(all_lambdas['nonsmoking', 'from_110', 'mles'],
#                   top=3).keys()))

# top_genes_from_normal = set.union(
#     set(top_genes(all_lambdas['pan_data', 'from_normal', 'mles'],
#                   top=4).keys()).difference(
#                    {'CSMD3', 'TTN'}),
#     set(top_genes(all_lambdas['smoking', 'from_normal', 'mles'],
#                   top=5).keys()).difference(
#                    {'CSMD3', 'TTN'}),
#     set(top_genes(all_lambdas['nonsmoking', 'from_normal', 'mles'],
#                   top=5).keys()).difference(
#                    {'CSMD3', 'TTN'}))


# top_genes_no_epi = set()
# top_genes_from_110 = set()
# top_genes_from_normal = set(pts_per_mutation[:10].index)

# top_genes_epi = set.union(top_genes_from_110, top_genes_from_normal)


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

    all_plots = {}

    ## TODO: Needs to fix bug on import_figures, where keys are not matching for no_epi
    ## * No epistasis include TP53 and KRAS
    ## ** Pan data

    # key, epi_key = ('pan_data', 'no_epi')
    # print(f"Plotting {key}, {epi_key}...")

    # axis_args_lambdas= {
    #     'lim':1.7,
    #     'label':"Flux to gene",
    #     'tick_each':0.5,
    #     'tick_minor_each':0.1,
    #     'scale':1}

    # axis_args_gammas= {
    #     'lim':775,
    #     'label':("Selection of gene "
    #              "(thousands)"),
    #     'tick_each':100,
    #     'tick_minor_each':50,
    #     'scale':10**3}

    # all_plots[('scatter', key, epi_key)] = plot_lambdas_gammas(
    #     key,
    #     epi_key,
    #     genes=set.union(top_genes_no_epi, top_genes_epi),
    #     genes_to_annotate=top_genes_no_epi,
    #     axis_args_lambdas=axis_args_lambdas,
    #     axis_args_gammas=axis_args_gammas,
    #     left_top_genes={'KEAP1', 'SPATA3'},
    #     bottom_genes={'RYR2'},
    #     plot_name = f"fluxes_and_selections_{key}_{epi_key}_with_KRAS_TP53.png")

    # print("...done.")
    # print("")


    ## TODO: Needs to fix bug on import_figures, where keys are not matching for no_epi
    ## * No epistasis include TP53 and KRAS
    ## ** Smoking

    # key, epi_key = ('smoking', 'no_epi')
    # print(f"Plotting {key}, {epi_key}...")

    # all_plots[('scatter', key, epi_key)] = plot_lambdas_gammas(
    #     key,
    #     epi_key,
    #     genes=set.union(top_genes_no_epi, top_genes_epi),
    #     genes_to_annotate=top_genes_no_epi,
    #     axis_args_lambdas=axis_args_lambdas,
    #     axis_args_gammas=axis_args_gammas,
    #     left_top_genes={'KEAP1', 'BRAF'},
    #     left_bottom_genes={'SPATA3'},
    #     bottom_genes={'STK11', 'RYR2'},
    #     title="Smoking",
    #     plot_name = f"fluxes_and_selections_{key}_{epi_key}_with_KRAS_TP53.png")

    # print("...done.")
    # print("")


    ## TODO: Needs to fix bug on import_figures, where keys are not matching for no_epi
    ## * No epistasis include TP53 and KRAS but not EGFR
    ## ** Non-smoking

    # key, epi_key = ('nonsmoking', 'no_epi')
    # print(f"Plotting {key}, {epi_key}...")

    # all_plots[('scatter', key, epi_key)] = plot_lambdas_gammas(
    #     key,
    #     epi_key,
    #     genes=set.union(top_genes_no_epi, top_genes_epi),
    #     genes_to_annotate=top_genes_no_epi,
    #     axis_args_lambdas=axis_args_lambdas,
    #     axis_args_gammas=axis_args_gammas,
    #     bottom_genes={'RYR2', 'KRAS', 'BRAF'},
    #     title="Non-smoking",
    #     plot_name = f"fluxes_and_selections_{key}_{epi_key}_with_KRAS_TP53_no_EGRF.png")

    # print("...done.")
    # print("")


    ## TODO: Needs to fix bug on import_figures, where keys are not matching for no_epi
    ## * No epistasis include TP53 and KRAS with EGFR
    ## ** Non-smoking

    # axis_args_lambdas= {
    #     'lim':1.7,
    #     'label':"Flux to gene",
    #     'tick_each':0.5,
    #     'tick_minor_each':0.1,
    #     'scale':1}

    # axis_args_gammas= {
    #     'lim':1.55,
    #     'label':("Selection of gene "
    #              "(millions)"),
    #     'tick_each':0.5,
    #     'tick_minor_each':0.1,
    #     'scale':10**6}

    # all_plots[('scatter', key, epi_key)] = plot_lambdas_gammas(
    #     key,
    #     epi_key,
    #     genes=set.union(top_genes_no_epi, top_genes_epi),
    #     genes_to_annotate=top_genes_no_epi,
    #     axis_args_lambdas=axis_args_lambdas,
    #     axis_args_gammas=axis_args_gammas,
    #     bottom_genes={'RYR2', 'BRAF'},
    #     left_top_genes={'KRAS', 'SPATA3'},
    #     title="Non-smoking",
    #     plot_name = f"fluxes_and_selections_{key}_{epi_key}_with_KRAS_TP53.png")

    # print("...done.")
    # print("")


    ## TODO: Needs to fix bug on import_figures, where keys are not matching for no_epi
    ## * No epistasis zoom in
    ## ** Pan data

    # key, epi_key = ('pan_data', 'no_epi')
    # print(f"Plotting {key}, {epi_key}...")

    # axis_args_lambdas= {
    #     'lim':0.44,
    #     'label':"Flux to gene",
    #     'tick_each':0.1,
    #     'tick_minor_each':0.01,
    #     'scale':1}

    # axis_args_gammas= {
    #     'lim':160,
    #     'label':("Selection of gene "
    #              "(thousands)"),
    #     'tick_each':50,
    #     'tick_minor_each':10,
    #     'scale':10**3}

    # all_plots[('scatter', key, epi_key)] = plot_lambdas_gammas(
    #     key,
    #     epi_key,
    #     genes=set.union(top_genes_no_epi, top_genes_epi).difference(
    #         {'TP53', 'KRAS'}),
    #     genes_to_annotate=top_genes_no_epi,
    #     axis_args_lambdas=axis_args_lambdas,
    #     axis_args_gammas=axis_args_gammas,
    #     left_top_genes={'KEAP1', 'SPATA3'})

    # print("...done.")
    # print("")


    ## TODO: Needs to fix bug on import_figures, where keys are not matching for no_epi
    ## * No epistasis zoom in
    ## ** Smoking

    # key, epi_key = ('smoking', 'no_epi')
    # print(f"Plotting {key}, {epi_key}...")

    # all_plots[('scatter', key, epi_key)] = plot_lambdas_gammas(
    #     key,
    #     epi_key,
    #     genes=set.union(top_genes_no_epi, top_genes_epi).difference(
    #         {'TP53', 'KRAS'}),
    #     genes_to_annotate=top_genes_no_epi,
    #     axis_args_lambdas=axis_args_lambdas,
    #     axis_args_gammas=axis_args_gammas,
    #     left_top_genes={'BRAF', 'SPATA3'},
    #     title="Smoking",
    #     bottom_genes={'STK11'})

    # print("...done.")
    # print("")


    ## TODO: Needs to fix bug on import_figures, where keys are not matching for no_epi
    ## * No epistasis zoom in
    ## ** Non-smoking

    # key, epi_key = ('nonsmoking', 'no_epi')
    # print(f"Plotting {key}, {epi_key}...")

    # all_plots[('scatter', key, epi_key)] = plot_lambdas_gammas(
    #     key,
    #     epi_key,
    #     genes=set.union(top_genes_no_epi, top_genes_epi).difference(
    #         {'TP53', 'KRAS'}),
    #     genes_to_annotate=top_genes_no_epi,
    #     axis_args_lambdas=axis_args_lambdas,
    #     axis_args_gammas=axis_args_gammas,
    #     title="Non-smoking")

    # print("...done.")
    # print("")


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

    all_plots[('scatter', key, epi_key)] = plot_lambdas_gammas(
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

    all_plots[('scatter', key, epi_key)] = plot_lambdas_gammas(
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

    all_plots[('scatter', key, epi_key)] = plot_lambdas_gammas(
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

    all_plots[('scatter', key, epi_key)] = plot_lambdas_gammas(
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

    all_plots[('scatter', key, epi_key)] = plot_lambdas_gammas(
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

    all_plots[('scatter', key, epi_key)] = plot_lambdas_gammas(
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

    all_plots[('scatter', key, epi_key)] = plot_lambdas_gammas(
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

    all_plots[('scatter', key, epi_key)] = plot_lambdas_gammas(
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

    all_plots[('scatter', key, epi_key)] = plot_lambdas_gammas(
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

    all_plots[('scatter', key, epi_key)] = plot_lambdas_gammas(
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

    all_plots[('scatter', key, epi_key)] = plot_lambdas_gammas(
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


    return all_plots



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

    plt.gcf().set_facecolor('none')

    if plot_name is None:
        plot_name = "lc_types_by_smoking_mosaic_plot.png"
    fig.savefig(os.path.join(location_figures,
                             plot_name),
                dpi=1000,
                transparent=True)
    plt.close('all')
    return fig



def plot_xs_vs_ys_scatter(xs_key,
                          ys_key=None,
                          xs="lambdas",
                          ys="lambdas",
                          xs_epi_key='from_normal',
                          ys_epi_key='from_110',
                          genes=None,
                          with_cis=True,
                          axis_args_xs=None,
                          axis_args_ys=None,
                          genes_to_annotate=None,
                          always_annotate_genes=None,
                          xs_cut_for_annotate=np.inf,
                          ys_cut_for_annotate=np.inf,
                          annotate_adjust_x=0,
                          annotate_adjust_y=0,
                          bottom_genes=None,
                          left_top_genes=None,
                          left_bottom_genes=None,
                          title=None,
                          include_x_equal_y_line=True,
                          plot_name=None,
                          dpi=300):
    """Create a scatter plot with the xs and the ys values for each
    gene.

    :type xs_key: str or NoneType
    :param xs_key: What estimates to use. Can be one of:
        - 'pan_data'
        - 'smoking'
        - 'nonsmoking'
        - 'smoking_plus'
        - 'nonsmoking_plus'

    :type ys_key: str or NoneType
    :param ys_key: What estimates to use. Can be one of:
        - 'pan_data'
        - 'smoking'
        - 'nonsmoking'
        - 'smoking_plus'
        - 'nonsmoking_plus'
        If None is provided, then use `xs_key`

    :type xs: str
    :param xs: What to use for xs values. Can be one of:
        - 'gammas' (selection coefficients, default)
        - 'lambdas' (fluxes)
        - 'mus' (mutation rates)

    :type ys: str
    :param ys: What to use for xs values. Can be one of:
        - 'gammas' (selection coefficients, default)
        - 'lambdas' (fluxes)
        - 'mus' (mutation rates)

    :type xs_epi_key: str
    :param xs_epi_key: If 'no_epi' (default), do not
        consider epistastic effects. Otherwise consider epistatic
        effects in a model with 3 genes that includes TP53, KRAS and a
        third gene and plot the estimates of fluxes and selection
        coefficients from either the normal genotype (if `epi_key` is
        'from_normal') or the TP53+KRAS genotype (if `epi_key` is
        'from_110') to the mutation of the third gene.

    :type ys_epi_key: str
    :param ys_epi_key: If 'no_epi' (default), do not
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

    :type include_x_equal_y_line: bool
    :param include_x_equal_y_line: If True, include an x=y line.

    :type plot_name: str or NoneType
    :param plot_name: Name used to save the plot. If None, use
        f'{xs}_{xs_epi_key}_vs_{ys}_{ys_epi_key}_{key}.png'.

    """
    if ys_key is None:
        ys_key = xs_key


    def translate_epi_key(epi_key):
        if epi_key == 'from_normal':
            return "from normal"
        elif epi_key == 'from_110':
            return "from TP53+KRAS"

    if xs == "gammas":
        xs_mles = all_gammas[(xs_key, xs_epi_key, 'mles')]
        if with_cis:
            xs_cis = all_gammas[(xs_key, xs_epi_key, 'cis')]
        if axis_args_xs is None:
            axis_args_xs = default_axis_args_gammas


    elif xs == "lambdas":
        xs_mles = all_lambdas[(xs_key, xs_epi_key, 'mles')]
        if with_cis:
            xs_cis = all_lambdas[(xs_key, xs_epi_key, 'cis')]
        if axis_args_xs is None:
            axis_args_xs = default_axis_args_lambdas.copy()
            axis_args_xs['label'] = f"Flux {translate_epi_key(xs_epi_key)} to 'gene'"

    elif xs == "mus":
        xs_mles = all_mus[(key, xs_epi_key, 'mles')]


    if ys == "gammas":
        ys_mles = all_gammas[(key, ys_epi_key, 'mles')]
        if with_cis:
            ys_cis = all_gammas[(key, ys_epi_key, 'cis')]
        if axis_args_ys is None:
            axis_args_ys = default_axis_args_gammas

    elif ys == "lambdas":
        ys_mles = all_lambdas[(key, ys_epi_key, 'mles')]
        if with_cis:
            ys_cis = all_lambdas[(key, ys_epi_key, 'cis')]
        if axis_args_ys is None:
            axis_args_ys = default_axis_args_lambdas.copy()
            axis_args_ys['label'] = f"Flux {translate_epi_key(ys_epi_key)} to 'gene'"

    elif ys == "mus":
        ys_mles = all_mus[(key, ys_epi_key, 'mles')]
        if with_cis:
            ys_cis = all_mus[(key, ys_epi_key, 'cis')]


    if genes is None:
        genes = set.intersection(set(xs_mles.keys()),
                                 set(ys_mles.keys()))
        print(f'For {key} and {xs_epi_key}:')
        print(f'  genes in xs MLEs: {len(xs_mles.keys())}')
        print(f'  genes in ys MLEs: {len(ys_mles.keys())}')

        if with_cis:
            genes = set.intersection(genes,
                                     set(xs_cis.keys()),
                                     set(ys_cis.keys()))
            print(f'  genes in xs CIs: {len(xs_cis.keys())}')
            print(f'  genes in ys CIs: {len(ys_cis.keys())}')

        print(f'  genes to plot (intersection of above): {len(genes)}')


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

    if always_annotate_genes is None:
        always_annotate_genes = set()

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ## This needs to make sure that xs and ys are in the same order
    # genes_in_both_sets = set.union(xs_mles.keys(), ys_mles.keys())
    # slope = curve_fit(lambda x, s: s*x,
    #                   list(xs_mles.values()),
    #                   list(ys_mles.values()))[0][0]
    # ax.plot([0, axis_args_xs['lim']],
    #         [0, slope*axis_args_xs['lim']],
    #         "--",
    #         color="Purple",
    #         alpha=0.1)

    genes_to_remove = set()
    for gene in genes:
        if gene not in xs_mles.keys():
            print(f"Gene {gene} not in xs, removing from plot")
            genes_to_remove.add(gene)
        if gene not in ys_mles.keys():
            print(f"Gene {gene} not in ys, removing from plot")
            genes_to_remove.add(gene)

    genes = set.difference(genes, genes_to_remove)

    for gene in genes:
        if xs_mles[gene] < ys_mles[gene]:
            color = 'Red'
            if ys_mles[gene] < ys_cut_for_annotate:
                genes_to_annotate = set.difference(genes_to_annotate, {gene})
        elif xs_mles[gene] > ys_mles[gene]:
            color = 'Blue'
            if xs_mles[gene] < xs_cut_for_annotate:
                genes_to_annotate = set.difference(genes_to_annotate, {gene})
        else:
            color = 'Purple'

        if gene in always_annotate_genes:
            genes_to_annotate = set.union(genes_to_annotate, {gene})


        if gene in genes_to_annotate:
            print(f"Annotating {gene} in {color}, value: ({xs_mles[gene]}, {ys_mles[gene]})")
            ax.text((xs_mles[gene] - annotate_adjust_x
                     if gene in set.union(left_top_genes, left_bottom_genes)
                     else xs_mles[gene] + annotate_adjust_x),
                    (ys_mles[gene] - 3*annotate_adjust_y
                     if gene in set.union(bottom_genes, left_bottom_genes)
                     else ys_mles[gene] + annotate_adjust_y),
                    gene,
                    ha=('right' if gene in set.union(left_top_genes, left_bottom_genes)
                        else 'left'),
                    va=('top' if gene in set.union(bottom_genes, left_bottom_genes)
                        else 'bottom'),
                    fontsize=6 if gene in {'ERBB4', 'GRM1', 'SLIT3', 'KMT2D'} else 7,
                    zorder=3,
                    color=color)

            if with_cis:
                ax.plot([xs_mles[gene], xs_mles[gene]],
                        ys_cis[gene],
                        color=color,
                        lw=0.5,
                        alpha=0.1,
                        zorder=4)
                ax.plot(xs_cis[gene],
                        [ys_mles[gene], ys_mles[gene]],
                        color=color,
                        lw=0.5,
                        alpha=0.2,
                        zorder=4)
            ax.scatter(xs_mles[gene],
                       ys_mles[gene],
                       marker='.',
                       color=color,
                       alpha=0.5,
                       zorder=5)
        else:
            ax.scatter(xs_mles[gene],
                       ys_mles[gene],
                       marker='.',
                       color=color,
                       alpha=0.25,
                       zorder=5)


    print("")
    if include_x_equal_y_line:
        x_equal_y_xs = np.linspace(0,
                                   max(axis_args_xs['lim'], axis_args_ys['lim']),
                                   1000)

        ax.plot(x_equal_y_xs, x_equal_y_xs, "--", color="Purple", alpha=0.7)

    ax.spines['bottom'].set_position('zero')  # Set x-axis position to zero
    ax.spines['bottom'].set_zorder(1)  # Set a lower zorder for the x-axis

    # Hide the top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


    if title is not None:
        ax.set_title(title)

    ax = prettify_axis(ax, 'x',
                       axis_args_xs)
    ax = prettify_axis(ax, 'y',
                       axis_args_ys)


    if plot_name is None:
        plot_name = f'{xs}_{xs_epi_key}_{xs_key}_vs_{ys}_{ys_epi_key}_{ys_key}.png'
    fig.savefig(os.path.join(location_figures,
                             plot_name),
                dpi=dpi)
    plt.close('all')

    return fig





if __name__ == "__figures__":
    plot_all()


mutations_per_gene = {}

for key in ['smoking', 'nonsmoking']:

    mutations_per_gene[key] = {gene: (samples_per_combination[key].loc[gene].loc['(0, 0, 1)']
                                      + samples_per_combination[key].loc[gene].loc['(1, 0, 1)']
                                      + samples_per_combination[key].loc[gene].loc['(0, 1, 1)']
                                      + samples_per_combination[key].loc[gene].loc['(1, 1, 1)'])
                                   for gene in samples_per_combination[key].index}


    mutations_per_gene[key].update({'TP53':(samples_per_combination[key].loc['EGFR'].loc['(1, 0, 0)']
                                            + samples_per_combination[key].loc['EGFR'].loc['(1, 0, 1)']
                                            + samples_per_combination[key].loc['EGFR'].loc['(1, 1, 0)']
                                            + samples_per_combination[key].loc['EGFR'].loc['(1, 1, 1)']),
                                    'KRAS':(samples_per_combination[key].loc['EGFR'].loc['(0, 1, 0)']
                                            + samples_per_combination[key].loc['EGFR'].loc['(0, 1, 1)']
                                            + samples_per_combination[key].loc['EGFR'].loc['(1, 1, 0)']
                                            + samples_per_combination[key].loc['EGFR'].loc['(1, 1, 1)'])})


    mutations_per_gene[key] = {k: v for k, v in sorted(mutations_per_gene[key].items(), key=lambda item: item[1])}


def plot_figures_poster():
    all_plots = {}

    ## * Epistasis only TP53 and KRAS
    print(f"Plotting TP53 and KRAS only model")

    all_plots['tp53_kras_landscape_nonsmoking', 'gammas'] = plot_landscape(
        results_TP53_KRAS_nonsmoking[('gammas', 'mle')],
        results_TP53_KRAS_nonsmoking['samples'],
        ['TP53', 'KRAS'],
        positions='left_to_right',
        scale_arrows=0.5*10**(-6),
        scale_circles=0.025,
        # name_y_offsets=name_y_offsets,
        # name_x_offsets=name_x_offsets,
        plot_name='landscape_gammas_tp53_kras_nonsmoking')


    all_plots['tp53_kras_landscape_smoking', 'gammas'] = plot_landscape(
        results_TP53_KRAS_smoking[('gammas', 'mle')],
        results_TP53_KRAS_smoking['samples'],
        ['TP53', 'KRAS'],
        positions='left_to_right',
        scale_arrows=0.5*10**(-6),
        scale_circles=0.025,
        # name_y_offsets=name_y_offsets,
        # name_x_offsets=name_x_offsets,
        plot_name='landscape_gammas_tp53_kras_smoking')



    ## * Epistasis TP53, KRAS, EGFR
    print(f"Plotting TP53, KRAS, EGFR model")

    ## We overide this value that was actually non-computable
    results_TP53_KRAS_EGFR_nonsmoking['gammas', 'mle'][((0, 1, 1), (1, 1, 1))] = 0

    all_plots['tp53_kras_egfr_landscape_nonsmoking', 'gammas'] = plot_landscape(
        results_TP53_KRAS_EGFR_nonsmoking[('gammas', 'mle')],
        results_TP53_KRAS_EGFR_nonsmoking['samples'],
        ['TP53', 'KRAS', 'EGFR'],
        positions='left_to_right',
        scale_arrows=0.5*10**(-6),
        scale_circles=0.025,
        # name_y_offsets=name_y_offsets,
        # name_x_offsets=name_x_offsets,
        plot_name='landscape_gammas_tp53_kras_egfr_nonsmoking')


    all_plots['tp53_kras_egfr_landscape_smoking', 'gammas'] = plot_landscape(
        results_TP53_KRAS_EGFR_smoking[('gammas', 'mle')],
        results_TP53_KRAS_EGFR_smoking['samples'],
        ['TP53', 'KRAS', 'EGFR'],
        positions='left_to_right',
        scale_arrows=0.5*10**(-6),
        scale_circles=0.025,
        # name_y_offsets=name_y_offsets,
        # name_x_offsets=name_x_offsets,
        plot_name='landscape_gammas_tp53_kras_egfr_smoking')



    ## * Flux from normal vs TP53+KRAS smoking

    print(f"Plotting flux from normal vs from TP53+KRAS for smoking")

    axis_args_xs= {
        'lim':0.45,
        'label':"Flux from normal to 'gene'",
        'tick_each':0.1,
        'tick_minor_each':0.01,
        'scale':1}
    axis_args_ys= {
        'lim':0.7,
        'label':"Flux to from TP53+KRAS to 'gene'",
        'tick_each':0.1,
        'tick_minor_each':0.01,
        'scale':1}

    plot_xs_vs_ys_scatter("smoking_plus",
                          genes_to_annotate='all',
                          axis_args_xs=axis_args_xs,
                          axis_args_ys=axis_args_ys,
                          xs_cut_for_annotate=0.2,
                          ys_cut_for_annotate=0.4,
                          annotate_adjust_x=0.002,
                          annotate_adjust_y=0.002,
                          always_annotate_genes={'NAV3', 'RYR1', 'PXDNL', 'STK11'},
                          # bottom_genes={'KEAP1'},
                          left_top_genes={'STK11'},
                          dpi=150
                          # dpi=1500
)


    print(f"Plotting flux from normal vs from TP53+KRAS for nonsmoking")

    axis_args_xs= {
        'lim':0.07,
        'label':"Flux from normal to 'gene'",
        'tick_each':0.01,
        'tick_minor_each':0.005,
        'scale':1}
    axis_args_ys= {
        'lim':0.7,
        'label':"Flux to from TP53+KRAS to 'gene'",
        'tick_each':0.1,
        'tick_minor_each':0.01,
        'scale':1}

    plot_xs_vs_ys_scatter("nonsmoking_plus",
                          genes_to_annotate='all',
                          axis_args_xs=axis_args_xs,
                          axis_args_ys=axis_args_ys,
                          annotate_adjust_x=0.0001,
                          annotate_adjust_y=0.002,
                          xs_cut_for_annotate=0.03,
                          ys_cut_for_annotate=0.4,
                          always_annotate_genes={'RYR1', 'SETBP1', 'FAT1', 'APC', 'NAV3', 'LRP1B'},
                          bottom_genes={'RBM10', 'RB1'},
                          left_top_genes={'CTNNB1'},
                          left_bottom_genes={'SMAD4', 'FAT1'},
                          dpi=150
                          # dpi=1500
                          )






    return all_plots
