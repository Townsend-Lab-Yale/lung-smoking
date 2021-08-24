import numpy as np
import matplotlib.pyplot as plt

from locations import fluxes_mles_file_name
from locations import fluxes_cis_file_name
from locations import location_figures

lambdas_mles = np.load(fluxes_mles_file_name, allow_pickle=True).item()
lambdas_cis = np.load(fluxes_cis_file_name, allow_pickle=True).item()

# mus = {gene:0.1 for gene in lambdas_mles.keys()} ## need to load mut_rates here
# gammas_mles = {gene:{x_y:mle/mus[gene]
#                      for x_y, mle in lambdas.items()}
#                for gene, lambdas in lambdas_mles.items()}
# gammas_cis = {gene:{x_y:ci/mus[gene]
#                     for x_y, ci in lambdas.items()}
#               for gene, lambdas in lambdas_cis.items()}


def plot_estimates_ordered(all_mles, all_cis, plot_name=None):
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

    x_lim = 0.4
    x_tick_each = 0.1
    x_minor_tick_each = 0.05
    ax.set_xticks(np.arange(0, x_lim+x_tick_each, x_tick_each))
    ax.set_xticks(np.arange(0, x_lim+x_minor_tick_each, x_minor_tick_each),
                  minor=True)
    ax.set_xlim(0, x_lim)
    ax.set_xlabel("Flux from TP53+KRAS to TP53+KRAS+'gene'")


    plt.tight_layout()
    if plot_name is None:
        fig.savefig(os.path.join(location_figures,
                                 "estimates_ordered.png"))
    else:
        fig.savefig(os.path.join(location_figures,
                                 f"{plot_name}.png"))
    plt.close('all')
    return mles_110_to_111


def plot_all():
    plot_estimates_ordered(lambdas_mles, lambdas_cis,
                           plot_name='fluxes_ordered.png')
    # plot_estimates_ordered(gammas_mles, gammas_cis,
    #                        plot_name='selection_coeffs_ordered.png')


if __name__ == "__figures__":
    plot_all()
