"""Matrix plotting -- Jorge A. Alfaro-Murillo

This module provide functions to produce figures where a matrix of
results is plotted as a heatmap.

We have used it so far to provide a visual representation of pair-wise
epistatic ratios for a set of genes, by representing what we define as
the epistatic ratio: the ratio between

"""

import os
import numpy as np

from epistatic_ratios import epistatic_ratios_2_matrix

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec

from locations import location_figures


class MidpointNormalize(mcolors. Normalize):
    def __init__(self, vmin=None, vmax=None, vcenter=None, clip=False):
        self.vcenter = vcenter
        super().__init__(vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.vcenter, self.vmax], [0, 0.5, 1.]
        return np.ma.masked_array(np.interp(value, x, y,
                                            left=-np.inf, right=np.inf))

    def inverse(self, value):
        y, x = [self.vmin, self.vcenter, self.vmax], [0, 0.5, 1]
        return np.interp(value, x, y, left=-np.inf, right=np.inf)



def plot_epistatic_ratios_2_matrices(results_nonsmoking,
                                     results_nonsmoking_cis,
                                     results_smoking,
                                     results_smoking_cis,
                                     gene_list_by_selection,
                                     third_gene_effects=None,
                                     plot_name=None,
                                     axis_title_size=18,
                                     axis_label_size=16):
    """Creat a plot with the epistatic ratios comparing non-smoking
    with smoking.

    This plot has a the matrix for smoking and non-smoking side to
    side with the colorbar in the middle.

    third_gene_effects if provided should be a tuple with third gene
    effects for nonsmoking and smoking, each element as return by
    epistatic_ratios_3rd_gene_effects. If included then the plot puts
    a small triangle whenever there is a third gene interaction that
    makes the epistatic ratio no longer be (statistically
    significantly) synergistic (above 1) or antagonistic (bellow 1).

    """

    matrix_smoking = epistatic_ratios_2_matrix(results_smoking,
                                               gene_list_by_selection,
                                               results_smoking_cis)
    matrix_nonsmoking = epistatic_ratios_2_matrix(results_nonsmoking,
                                                  gene_list_by_selection,
                                                  results_nonsmoking_cis)



    if plot_name is None:
        plot_name = "epistatic_ratios_matrices"
        if third_gene_effects is not None:
            plot_name = plot_name + "_with_third_gene_effects"

    fig = plt.figure(figsize=(15, 8))

    # Axes for the heatmaps, with center space for a colorbar
    gs = gridspec.GridSpec(1, 3,
                           width_ratios=[1, 0.05, 1])
    ax_smoking = fig.add_subplot(gs[0])
    ax_nonsmoking = fig.add_subplot(gs[2])


    # Color map
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "truncated_seismic_r",
        plt.cm.seismic_r(np.linspace(0.15, 0.85, 256)))
    cmap.set_bad(np.array([0.7, 0.7, 0.7, 1]), 1.0) # set to gray nan
                                                    # values; i.e. the
                                                    # diagonal

    midnorm = MidpointNormalize(vmin=0, vcenter=1, vmax=30)

    # Pseudocolor plots are the heatmaps (matrices represented by colors)
    pcm_smoking = ax_smoking.pcolormesh(
        np.arange(len(gene_list_by_selection)),
        np.arange(len(gene_list_by_selection)),
        matrix_smoking[::-1],
        rasterized=True,
        norm=midnorm,
        cmap=cmap,
        shading='auto')

    pcm_nonsmoking = ax_nonsmoking.pcolormesh(
        np.arange(len(gene_list_by_selection)),
        np.arange(len(gene_list_by_selection)),
        matrix_nonsmoking[::-1],
        rasterized=True,
        norm=midnorm,
        cmap=cmap,
        shading='auto')


    # ax_smoking.set_title("smokers", fontsize=18)
    ax_smoking.set_aspect('equal', adjustable='box')

    ax_smoking.set_xticks(np.arange(len(gene_list_by_selection)))
    ax_smoking.set_xticklabels(gene_list_by_selection, rotation=90, fontsize=axis_label_size, style='italic')
    ax_smoking.set_xlabel("Context (mutated gene in somatic genotype)", fontsize=axis_title_size)

    ax_smoking.set_yticks(np.arange(len(gene_list_by_selection)))
    ax_smoking.set_yticklabels(gene_list_by_selection[::-1], fontsize=axis_label_size, style='italic')
    ax_smoking.set_ylabel("Mutation under selection", fontsize=axis_title_size)

    ax_smoking.set_xticks(np.arange(0.5, len(gene_list_by_selection), 1), minor=True)
    ax_smoking.set_yticks(np.arange(0.5, len(gene_list_by_selection), 1), minor=True)
    ax_smoking.grid(which='minor', color='gray', linestyle='-', linewidth=1)

    if third_gene_effects is not None:
        for genes, value in third_gene_effects[0].items():
            x = gene_list_by_selection.index(genes[0])
            y = gene_list_by_selection[::-1].index(genes[1])
            xs = [x-0.5, x, x+0.5]
            ys = [y+0.5, y, y+0.5]
            ax_nonsmoking.fill(xs, ys, color=pcm_nonsmoking.to_rgba(value), edgecolor='none')

        for genes, value in third_gene_effects[1].items():
            x = gene_list_by_selection.index(genes[0])
            y = gene_list_by_selection[::-1].index(genes[1])
            xs = [x-0.5, x, x+0.5]
            ys = [y+0.5, y, y+0.5]
            ax_smoking.fill(xs, ys, color=pcm_smoking.to_rgba(value), edgecolor='none')


    # ax_nonsmoking.set_title("non-smokers", fontsize=18)
    ax_nonsmoking.set_aspect('equal', adjustable='box')

    ax_nonsmoking.set_xticks(np.arange(len(gene_list_by_selection)))
    ax_nonsmoking.set_xticklabels(gene_list_by_selection, rotation=90, fontsize=axis_label_size, style='italic')
    ax_nonsmoking.set_xlabel("Context (mutated gene in somatic genotype)", fontsize=axis_title_size)

    ax_nonsmoking.set_yticks(np.arange(len(gene_list_by_selection)))
    ax_nonsmoking.set_yticklabels(gene_list_by_selection[::-1], fontsize=axis_label_size, style='italic')
    ax_nonsmoking.set_ylabel("Mutation under selection", fontsize=axis_title_size)

    ax_nonsmoking.set_xticks(np.arange(0.5, len(gene_list_by_selection), 1), minor=True)
    ax_nonsmoking.set_yticks(np.arange(0.5, len(gene_list_by_selection), 1), minor=True)
    ax_nonsmoking.grid(which='minor', color='gray', linestyle='-', linewidth=1)

    ax_nonsmoking.yaxis.set_ticks_position("right")
    ax_nonsmoking.yaxis.set_label_position("right")


    ## Color bar in the center of the two matrices

    cb_ax = fig.add_axes([0.505, 0.22, 0.01, 0.56])  # Adjust these values as needed

    cb = fig.colorbar(pcm_smoking, cax=cb_ax,
                      orientation='vertical',
                      fraction=0.01,
                      pad=0.1)
    cb.set_ticks([0, 1, 10, 20, 30])
    cb.set_label('Epistatic ratio')
    cb.ax.yaxis.set_label_position('left')



    plt.subplots_adjust(wspace=0.13)


    dpi = 200
    fig.savefig(os.path.join(location_figures,
                             plot_name),
                dpi=dpi)
    plt.close('all')

    return fig




def plot_epistatic_ratios_2_matrices_poster(results_nonsmoking,
                                            results_nonsmoking_cis,
                                            results_smoking,
                                            results_smoking_cis,
                                            gene_list_by_selection,
                                            third_gene_effects=None,
                                            plot_name=None):
    """Creat a plot with the epistatic ratios.

    This plot has a vertical orientation.

    third_gene_effects if provided should be a tuple with third gene
    effects for nonsmoking and smoking, each element as return by
    epistatic_ratios_3rd_gene_effects. If included then the plot puts
    a small triangle whenever there is a third gene interaction that
    makes the epistatic ratio no longer be (statistically
    significantly) synergistic (above 1) or antagonistic (bellow 1).

    """

    matrix_smoking = epistatic_ratios_2_matrix(results_smoking,
                                               gene_list_by_selection,
                                               results_smoking_cis)
    matrix_nonsmoking = epistatic_ratios_2_matrix(results_nonsmoking,
                                                  gene_list_by_selection,
                                                  results_nonsmoking_cis)



    if plot_name is None:
        plot_name = "epistatic_ratios_matrices_poster"
        if third_gene_effects is not None:
            plot_name = plot_name + "_with_third_gene_effects"

    fig = plt.figure(figsize=(8, 15))

    # Axes for the heatmaps, with center space for a colorbar
    gs = gridspec.GridSpec(3, 1,
                           height_ratios=[1, 0.05, 1])
    ax_smoking = fig.add_subplot(gs[0])
    ax_nonsmoking = fig.add_subplot(gs[2])


    # Color map
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "truncated_seismic_r",
        plt.cm.seismic_r(np.linspace(0.15, 0.85, 256)))
    cmap.set_bad(np.array([0.7, 0.7, 0.7, 1]), 1.0) # set to gray nan
                                                    # values; i.e. the
                                                    # diagonal

    midnorm = MidpointNormalize(vmin=0, vcenter=1, vmax=30)

    # Pseudocolor plots are the heatmaps (matrices represented by colors)
    pcm_smoking = ax_smoking.pcolormesh(
        np.arange(len(gene_list_by_selection)),
        np.arange(len(gene_list_by_selection)),
        matrix_smoking[::-1],
        rasterized=True,
        norm=midnorm,
        cmap=cmap,
        shading='auto')

    pcm_nonsmoking = ax_nonsmoking.pcolormesh(
        np.arange(len(gene_list_by_selection)),
        np.arange(len(gene_list_by_selection)),
        matrix_nonsmoking[::-1],
        rasterized=True,
        norm=midnorm,
        cmap=cmap,
        shading='auto')


    # ax_smoking.set_title("Ever smokers")
    ax_smoking.set_aspect('equal', adjustable='box')

    ax_smoking.set_xticks(np.arange(len(gene_list_by_selection)))
    ax_smoking.set_xticklabels(gene_list_by_selection, rotation=90)
    ax_smoking.set_xlabel("Context (mutated gene in somatic genotype)")

    ax_smoking.set_yticks(np.arange(len(gene_list_by_selection)))
    ax_smoking.set_yticklabels(gene_list_by_selection[::-1])
    ax_smoking.set_ylabel("Gene mutation under selection")

    ax_smoking.set_xticks(np.arange(0.5, len(gene_list_by_selection), 1), minor=True)
    ax_smoking.set_yticks(np.arange(0.5, len(gene_list_by_selection), 1), minor=True)
    ax_smoking.grid(which='minor', color='gray', linestyle='-', linewidth=1)

    ax_smoking.xaxis.set_ticks_position("top")
    ax_smoking.xaxis.set_label_position("top")

    if third_gene_effects is not None:

        for genes, value in third_gene_effects[0].items():
            x = gene_list_by_selection.index(genes[0])
            y = gene_list_by_selection[::-1].index(genes[1])
            xs = [x-0.5, x, x+0.5]
            ys = [y+0.5, y, y+0.5]
            ax_nonsmoking.fill(xs, ys, color=pcm_nonsmoking.to_rgba(value), edgecolor='none')

        for genes, value in third_gene_effects[1].items():
            x = gene_list_by_selection.index(genes[0])
            y = gene_list_by_selection[::-1].index(genes[1])
            xs = [x-0.5, x, x+0.5]
            ys = [y+0.5, y, y+0.5]
            ax_smoking.fill(xs, ys, color=pcm_smoking.to_rgba(value), edgecolor='none')


    # ax_nonsmoking.set_title("Never-smokers")
    ax_nonsmoking.set_aspect('equal', adjustable='box')

    ax_nonsmoking.set_xticks(np.arange(len(gene_list_by_selection)))
    ax_nonsmoking.set_xticklabels(gene_list_by_selection, rotation=90)
    ax_nonsmoking.set_xlabel("Context (mutated gene in somatic genotype)")

    ax_nonsmoking.set_yticks(np.arange(len(gene_list_by_selection)))
    ax_nonsmoking.set_yticklabels(gene_list_by_selection[::-1])
    ax_nonsmoking.set_ylabel("Gene mutation under selection")

    ax_nonsmoking.set_xticks(np.arange(0.5, len(gene_list_by_selection), 1), minor=True)
    ax_nonsmoking.set_yticks(np.arange(0.5, len(gene_list_by_selection), 1), minor=True)
    ax_nonsmoking.grid(which='minor', color='gray', linestyle='-', linewidth=1)

    # ax_nonsmoking.yaxis.set_ticks_position("right")
    # ax_nonsmoking.yaxis.set_label_position("right")


    ## Color bar in the center of the two matrices

    cb_ax = fig.add_axes([0.1, 0.50, 0.8, 0.01])  # Adjust these values as needed

    cb = fig.colorbar(pcm_smoking, cax=cb_ax,
                      orientation='horizontal',
                      fraction=0.01,
                      pad=0.1)
    cb.set_ticks([0, 1, 10, 20, 30])
    cb.set_label('Epistatic ratio')
    cb.ax.yaxis.set_label_position('left')



    plt.subplots_adjust(wspace=0.13)


    dpi = 200
    fig.savefig(os.path.join(location_figures,
                             plot_name),
                dpi=dpi)
    plt.close('all')

    return fig




def plot_epistatic_ratios_separate(results_nonsmoking,
                                   results_nonsmoking_cis,
                                   results_smoking,
                                   results_smoking_cis,
                                   gene_list_by_selection,
                                   third_gene_effects=None,
                                   plot_name=None):
    """Create a plot with the epistatic ratios.

    This function plots the matrices in different figures, so that it
    can easily placed when building a composed figure.

    third_gene_effects if provided should be a tuple with third gene
    effects for nonsmoking and smoking, each element as return by
    epistatic_ratios_3rd_gene_effects. If included then the plot puts
    a small triangle whenever there is a third gene interaction that
    makes the epistatic ratio no longer be (statistically
    significantly) synergistic (above 1) or antagonistic (bellow 1).

    """
    matrix_smoking = epistatic_ratios_2_matrix(results_smoking,
                                               gene_list_by_selection,
                                               results_smoking_cis)
    matrix_nonsmoking = epistatic_ratios_2_matrix(results_nonsmoking,
                                                  gene_list_by_selection,
                                                  results_nonsmoking_cis)

    if plot_name is None:
        plot_name = "epistatic_ratios_matrix"
        if third_gene_effects is not None:
            plot_name = plot_name + "_with_third_gene_effects"

    # Color map
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "truncated_seismic_r",
        plt.cm.seismic_r(np.linspace(0.15, 0.85, 256)))
    cmap.set_bad(np.array([0.7, 0.7, 0.7, 1]), 1.0)  # set to gray for nan values; i.e. the diagonal

    midnorm = MidpointNormalize(vmin=0, vcenter=1, vmax=30)

    # Create figure for smokers
    fig_smoking, ax_smoking = plt.subplots(figsize=(10, 8))
    pcm_smoking = ax_smoking.pcolormesh(
        np.arange(len(gene_list_by_selection)),
        np.arange(len(gene_list_by_selection)),
        matrix_smoking[::-1],
        rasterized=True,
        norm=midnorm,
        cmap=cmap,
        shading='auto')

    # ax_smoking.set_title("Smokers", fontsize=18)
    ax_smoking.set_aspect('equal', adjustable='box')

    ax_smoking.set_xticks(np.arange(len(gene_list_by_selection)))
    ax_smoking.set_xticklabels(gene_list_by_selection, rotation=90)
    ax_smoking.set_xlabel("Context (mutated gene in somatic genotype)")

    ax_smoking.set_yticks(np.arange(len(gene_list_by_selection)))
    ax_smoking.set_yticklabels(gene_list_by_selection[::-1])
    ax_smoking.set_ylabel("Gene mutation under selection")

    ax_smoking.set_xticks(np.arange(0.5, len(gene_list_by_selection), 1), minor=True)
    ax_smoking.set_yticks(np.arange(0.5, len(gene_list_by_selection), 1), minor=True)
    ax_smoking.grid(which='minor', color='gray', linestyle='-', linewidth=1)

    if third_gene_effects is not None:
        for genes, value in third_gene_effects[1].items():
            x = gene_list_by_selection.index(genes[0])
            y = gene_list_by_selection[::-1].index(genes[1])
            xs = [x - 0.5, x, x + 0.5]
            ys = [y + 0.5, y, y + 0.5]
            ax_smoking.fill(xs, ys, color=pcm_smoking.to_rgba(value), edgecolor='none')

    fig_smoking.subplots_adjust(bottom=0.2)
    fig_smoking.savefig(os.path.join(location_figures, f"{plot_name}_smoking.png"), dpi=200)
    plt.close(fig_smoking)

    # Create figure for non-smokers
    fig_nonsmoking, ax_nonsmoking = plt.subplots(figsize=(10, 8))
    pcm_nonsmoking = ax_nonsmoking.pcolormesh(
        np.arange(len(gene_list_by_selection)),
        np.arange(len(gene_list_by_selection)),
        matrix_nonsmoking[::-1],
        rasterized=True,
        norm=midnorm,
        cmap=cmap,
        shading='auto')

    # ax_nonsmoking.set_title("Non-smokers", fontsize=18)
    ax_nonsmoking.set_aspect('equal', adjustable='box')

    ax_nonsmoking.set_xticks(np.arange(len(gene_list_by_selection)))
    ax_nonsmoking.set_xticklabels(gene_list_by_selection, rotation=90)
    ax_nonsmoking.set_xlabel("Context (mutated gene in somatic genotype)")

    ax_nonsmoking.set_yticks(np.arange(len(gene_list_by_selection)))
    ax_nonsmoking.set_yticklabels(gene_list_by_selection[::-1])
    ax_nonsmoking.set_ylabel("Gene mutation under selection")

    ax_nonsmoking.set_xticks(np.arange(0.5, len(gene_list_by_selection), 1), minor=True)
    ax_nonsmoking.set_yticks(np.arange(0.5, len(gene_list_by_selection), 1), minor=True)
    ax_nonsmoking.grid(which='minor', color='gray', linestyle='-', linewidth=1)

    if third_gene_effects is not None:
        for genes, value in third_gene_effects[0].items():
            x = gene_list_by_selection.index(genes[0])
            y = gene_list_by_selection[::-1].index(genes[1])
            xs = [x - 0.5, x, x + 0.5]
            ys = [y + 0.5, y, y + 0.5]
            ax_nonsmoking.fill(xs, ys, color=pcm_nonsmoking.to_rgba(value), edgecolor='none')

    fig_nonsmoking.subplots_adjust(bottom=0.2)
    fig_nonsmoking.savefig(os.path.join(location_figures, f"{plot_name}_nonsmoking.png"), dpi=200)
    plt.close(fig_nonsmoking)

    # Create figure for colorbar
    fig_colorbar, ax_colorbar = plt.subplots(figsize=(2, 8))
    cb = fig_colorbar.colorbar(pcm_smoking, cax=ax_colorbar, orientation='vertical')
    cb.set_ticks([0, 1, 10, 20, 30])
    cb.set_label('Epistatic ratio')
    fig_colorbar.subplots_adjust(right=0.78)
    fig_colorbar.savefig(os.path.join(location_figures, f"{plot_name}_colorbar.png"), dpi=200)
    plt.close(fig_colorbar)


    return fig_smoking, fig_nonsmoking, fig_colorbar
