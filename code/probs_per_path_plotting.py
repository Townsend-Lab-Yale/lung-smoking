import os
import numpy as np

import matplotlib.pyplot as plt

from cancer_epistasis import compute_probability_paths

from locations import location_figures


plt.rcParams['text.usetex'] = True # If getting an error about ! LaTeX
                                   # Error: File `type1ec.sty' not
                                   # found, then install cm-super
                                   # LaTeX package, on debian systems:
                                   # sudo apt install cm-super


def path_to_string(path, genes, joiner, before, after):
    jumps = [xy[1] for xy in path]
    jumps = ([list(jumps[0])] +
             [list(np.array(y)-np.array(x))
              for x, y in zip(jumps[:-1], jumps[1:])])
    gene_jumps = [before + genes[jump.index(1)] + after
                  for jump in jumps]
    return joiner.join(gene_jumps)



def plot_probs_per_path(lambdas,
                        genes,
                        tolerance_to_not_show_path=None,
                        save_fig_name=None):
    """Plot the probabilities of every path to every somatic genotype
    with more than two mutations.

    :type lambdas: dict
    :param lambdas: A dictionary with the fluxes estimates with keys
        being tuples representing the somatic genotypes from and to of
        the respective flux (as for example the output of
        :func:`estimate_lambdas` with draws=1).

    :type genes: NoneType or list or tuple
    :param genes: List or tuple with the genes as strings. If None
        (default) then the paths will be represented as they are in
        :func:`generate_paths`.

    :type tolerance_to_not_show_path: float or NoneType
    :param tolerance_to_not_show_path: Paths with probability below
        this tolerance will not be shown in the results. If None
        (default) plot all paths even if they have a very small
        probability (that will look like a value of zero).

    :type save_fig_name: NoneType or str
    :param save_fig_name: File name to use to save the figure, if None
        is provided then do not save the figure just show it.

    :rtype: matplotlib.figure.Figure
    :return: Figure with the plot.

    """

    ## Compute probabilities of each path
    all_probs = compute_probability_paths(lambdas)


    ## Remove the ones that are too small
    if tolerance_to_not_show_path is not None:
        all_probs = {destination:{path:value for
                                  path, value in probs.items()
                                  if value > tolerance_to_not_show_path}
                     for destination, probs in all_probs.items()}

    ## Convert paths into strings with arrows
    all_probs = {destination:{path_to_string(path,
                                             genes,
                                             joiner=r"\; $\rightarrow$ ",
                                             before=r"$\textit{",
                                             after=r"}$"):value
                              for path, value in probs.items()}
                 for destination, probs in all_probs.items()}


    ## Actual plotting
    fig = plt.figure(figsize=(1.57*345/72.27, 8))
    ax = fig.add_subplot(111)

    yticks = [0]
    ytick_labels = []

    for genotype, probs_per_path in all_probs.items():

        genotype_y_label = yticks[-1] + (len(probs_per_path) - 1)/2

        genotype_label = " ".join([r"$\textit{" + genes[index] + r"}$"
                                   for index, i in enumerate(genotype)
                                   if i == 1])

        ax.text(1.04,
                genotype_y_label,
                genotype_label,
                va="center")

        for i , (path, prob) in enumerate(probs_per_path.items()):
            ax.barh(yticks[-1],
                    prob,
                    height=0.8,
                    color="gray")

            ytick_labels.append(path)

            if i != len(probs_per_path) - 1:
                yticks.append(yticks[-1]+1)
            else:
                ax.axhline(yticks[-1]+1.5, linestyle="--", color="gray")

        yticks.append(yticks[-1]+3)

    ax.text(1.04, yticks[-1]-1, "Resulting somatic genotype", va="center")
    ax.text(0, yticks[-1]-1, "Evolutionary trajectory", va="center", ha="right")

    ax.set_yticks(yticks[:-1])
    ax.set_yticklabels(ytick_labels)

    ax.set_ylim(-1.5, yticks[-1]-1.5)

    ax.set_xlim(0, 1.02)
    ax.set_xticks(np.arange(0, 1.1, 0.1), minor=True)
    ax.set_xticks(np.arange(0, 1.1, 0.5))
    ax.set_xticklabels([f"{int(round(x*100))}\%"
                        for x in np.arange(0, 1.1, 0.5)])
    ax.set_xlabel("Probability of trajectory")

    fig.subplots_adjust(left=0.3, right=0.72)

    if save_fig_name is not None:
        fig.savefig(os.path.join(location_figures,
                                 f"{save_fig_name}.png"),
                    transparent=True,
                    dpi=400)
        plt.close('all')
    else:
        plt.show()
        plt.close('all')

    return fig





## * Example of how to use the code

# from load_results import load_results
# all_lambdas_ = load_results('fluxes', 'mles')

# genes_= ('TP53', 'EGFR', 'KEAP1')
# key_ = 'smoking'

# fig_ = plot_probs_per_path(all_lambdas_[key_][genes_],
#                            genes=genes_,
#                            tolerance_to_not_show_path=10**(-20),
#                            save_fig_name=f"probs_per_path_{key_}_" + "_".join(genes_))


# ## If you want the actual values you can do

# probs_per_path_ = compute_probability_paths(all_lambdas_[key_][genes_])

# ## and for the keys to be easier to read:
# probs_per_path_ = {genotype:{path_to_string_(path, genes_, "->", "", ""):probs
#                              for path,probs in value.items()}
#                    for genotype, value in probs_per_path_.items()}
