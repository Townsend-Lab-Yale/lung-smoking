import os
import numpy as np
import pymc3 as pm

from scipy.special import binom

from itertools import permutations

import matplotlib.pyplot as plt
import matplotlib.cm as cm

from theory import numbers_positive_lambdas
from theory import build_S_as_array
from theory import build_S_with_tuples

from locations import location_figures
from locations import default_mutation_names


mutation_colors = [cm.get_cmap("tab20c")(x)
                   for x in np.linspace(0, 1, 20)]
"""Colors to use for each mutation in the figures. Can handle up to 20
mutations."""
mutation_colors = (mutation_colors[0::4]
                   + mutation_colors[3::4]
                   + mutation_colors[1::4]
                   + mutation_colors[2::4])



def positions_landscape_left_to_right(M):

    positions = {M*(0,):np.array([0, 0])}

    S = build_S_with_tuples(M)

    max_edges_per_col = np.max([np.sum(binom(M, i))
                                for i in range(M+1)])
    step_unit = 0.5*M/(+ max_edges_per_col//2
                       - (max_edges_per_col+1)%2 * 1/2)

    for m in range(1, M+1):
        x_with_m = [x for x in S if np.sum(x) == m]

        start_pos = (+ len(x_with_m)//2
                     - (len(x_with_m)+1)%2 * 1/2)*step_unit

        for i, x in enumerate(x_with_m):
            positions[x] = np.array([m, start_pos - i*step_unit])

            if m == 1:
                positions[tuple(-np.array(x))] = np.array([0, 0])

    return positions



def positions_landscape_out_to_in(M):

    positions = {M*(1,):np.array([0, 0])}

    if M == 1:
        positions[(0,)] = np.array([-1.1, 0])
        positions[(-1,)] = np.array([-1, 0])

    elif M == 2:
        positions[(0, 0)] = np.array([2.2, 0])
        positions[(-1, 0)] = np.array([-2, 0])
        positions[(0, -1)] = np.array([0, 2])
        positions[(1, 0)] = np.array([-1, 0])
        positions[(0, 1)] = np.array([0, 1])

    elif M == 3:
        positions[(0, 0, 0)] = np.array([0.8*3, 0.8*3])
        for i in range(M):
            angle = (np.pi/2 + 2*np.pi/3) + i*2*np.pi/3
            positions[(i*(0,) + (-1,) + (M-i-1)*(0,))] = (
                np.array([3*np.cos(angle),
                          3*np.sin(angle)]))
            positions[(i*(0,) + (1,) + (M-i-1)*(0,))] = (
                np.array([2*np.cos(angle),
                          2*np.sin(angle)]))
        positions[(1, 1, 0)] = (positions[(1, 0, 0)]
                                + positions[(0, 1, 0)])/2
        positions[(1, 0, 1)] = (positions[(1, 0, 0)]
                                + positions[(0, 0, 1)])/2
        positions[(0, 1, 1)] = (positions[(0, 1, 0)]
                                + positions[(0, 0, 1)])/2

    elif M == 4:
        positions[(0, 0, 0, 0)] = 1.1*np.array(
            [4*np.cos(np.pi/4), 4*np.sin(np.pi/4)])

        positions[(-1,  0,  0,  0)] = [-4,  0]
        positions[( 0, -1,  0,  0)] = [ 4,  0]
        positions[( 0,  0, -1,  0)] = [ 0,  4]
        positions[( 0,  0,  0, -1)] = [ 0, -4]

        positions[( 1,  0,  0,  0)] = [-3,  0]
        positions[( 0,  1,  0,  0)] = [ 3,  0]
        positions[( 0,  0,  1,  0)] = [ 0,  3]
        positions[( 0,  0,  0,  1)] = [ 0, -3]

        positions[( 1,  1,  0,  0)] = [-2,  0]
        positions[( 1,  0,  1,  0)] = np.array([-2,  2])/np.sqrt(2)
        positions[( 1,  0,  0,  1)] = np.array([-2, -2])/np.sqrt(2)
        positions[( 0,  1,  1,  0)] = np.array([ 2,  2])/np.sqrt(2)
        positions[( 0,  1,  0,  1)] = np.array([ 2, -2])/np.sqrt(2)
        positions[( 0,  0,  1,  1)] = [ 0, -2]

        positions[( 1,  1,  1,  0)] = [ 0,  1]
        positions[( 1,  1,  0,  1)] = [ 0, -1]
        positions[( 1,  0,  1,  1)] = [-1,  0]
        positions[( 0,  1,  1,  1)] = [ 1,  0]

        positions[(( 0,  1,  0,  0), (1, 1, 0, 0))] = [0,  1/2]
        positions[(( 0,  0,  1,  0), (0, 0, 1, 1))] = [1/2,  0]

    positions = {key:np.array(value)
                 for key, value in positions.items()}

    return positions



def positions_landscape(M, positions):
    """Establish the positions for the mutation landscape.

    :type positions: str
    :param positions: Positions of the mutation combinations for the
        landscape. Can either be 'out_to_in' or 'left_to_right' so
        that mutations (and fitness) are gain either from the outside
        to the inside or from left to right, of the landscape.

    :rtype: dict
    :return: Dictionary with the positions in the landscape for each
        mutation combination (dictionary keys, as tuples).

    """


    if positions == "out_to_in":
        return positions_landscape_out_to_in(M)
    elif positions == "left_to_right":
        return positions_landscape_left_to_right(M)
    else:
        raise Exception("Unknown 'positions'. Available options "
                        "are 'out_ot_in' or 'left_to_right'.")



def plot_landscape(arrows, circles, mutation_names=None,
                   scale_arrows=0.3, scale_circles=0.5,
                   include_n_circles=True,
                   positions="out_to_in",
                   subplot_label=None,
                   subplot_label_ha='left',
                   subplot_label_va='top',
                   name_x_offsets=None,
                   name_y_offsets=None,
                   multiplier_font_size=1,
                   plot_name=None):
    """Plot a mutation landscape with the fluxes, selection
    coefficients or mutation rates (arrow widths) and the counts of
    each of the mutation combinations (circle areas).

    :type arrows: dict
    :param arrows:

    :type circles: dict
    :param circles: Areas of the circles to plot. Generally, we are
        using this variable for the number of patients in each
        mutation combination, but it could be used to plot something
        else, for example the mus or gammas if the arrows are the
        lambdas.

    :type mutation_names: list or NoneType
    :param mutation_names: List with the names of mutations to be
        plotted. If None is provided use :const:`default_mutation_names`.

    :type positions: str or dict
    :param positions: Dictionary with the positions of the mutation
        combinations for the landscape with tuples as keys. Can also
        be a string equal to 'out_to_in' or 'left_to_right' so that
        the function :fun:`positions_landscape` is called to produce
        the positions.

    :type plot_name: str or NoneType
    :param plot_name: Name of the file to save the plot. If None
        (default), just call it 'landscape.png'.

    :rtype:
    :return:

    """

    M = numbers_positive_lambdas.index(len(arrows))
    S = build_S_as_array(M)

    if mutation_names is None:
        mutation_names = default_mutation_names[:M]
        mutation_names_sep = ''
    else:
        mutation_names_sep = '\n'

    colors = {i*(0,)+(1,)+(M-i-1)*(0,):mutation_colors[i]
              for i in range(M)}

    fig_width_pt = 345.0  # Get this from LaTeX using \showthe\columnwidth
    inches_per_pt = 1.0/72.27               # Convert pt to inches
    fig_width = fig_width_pt*inches_per_pt  # width in inches

    fig = plt.figure(figsize=[1.57*fig_width,
                              1.57*fig_width*4.3/4.62],
                     # the weird ratio above is to keep
                     # circles as circles it relates to
                     # the x and y lims of the ax
                     dpi=400)

    ax = fig.add_subplot(111)

    if isinstance(positions, str):
        pts = positions_landscape(M, positions)
    else:
        pts = positions


    for wz, arrow_wz in arrows.items():
        ## w is where the arrow comes from, z is where it goes
        w = wz[0]
        z = wz[1]

        to_xy = pts[z] ## xy are coordinates in the two-dimensional
                       ## space of the actual plot

        if np.sum(w) == 0:
            from_xy = pts[tuple(-np.array(z))]
        else:
            from_xy = pts[w]

        if wz in pts.keys():
            mid_xy = pts[wz]
            ax.arrow(from_xy[0],
                     from_xy[1],
                     (mid_xy[0]-from_xy[0]),
                     (mid_xy[1]-from_xy[1]),
                     width=scale_arrows*arrow_wz,
                     length_includes_head=True,
                     head_width=0,
                     head_length=0,
                     facecolor=colors[tuple(np.array(z)-np.array(w))],
                     linewidth=0,
                     alpha=0.7,
                     zorder=3)
            # ax.arrow(from_xy[0],
            #          from_xy[1],
            #          (mid_xy[0]-from_xy[0]),
            #          (mid_xy[1]-from_xy[1]),
            #          width=scale_arrows*mus[tuple(np.array(z)-np.array(w))],
            #          length_includes_head=True,
            #          head_width=0,
            #          head_length=0,
            #          facecolor='grey',
            #          linewidth=0,
            #          alpha=0.6)
            # ax.add_patch(
            #     plt.Circle(mid_xy,
            #                scale_arrows*arrow_wz/2,
            #                facecolor=colors[tuple(np.array(z)-np.array(w))],
            #                alpha=0.7))

            from_xy = mid_xy


        ax.arrow(from_xy[0],
                 from_xy[1],
                 (to_xy[0]-from_xy[0]),
                 (to_xy[1]-from_xy[1]),
                 width=scale_arrows*arrow_wz,
                 length_includes_head=True,
                 head_width=(1.5*scale_arrows*arrow_wz
                             if arrow_wz*scale_arrows > 0.01 else 0),
                 head_length=0.25,
                 facecolor=colors[tuple(np.array(z)-np.array(w))],
                 linewidth=0,
                 alpha=0.7,
                 zorder=3)

    ## Draw circles
    for m, circle in circles.items():
        ax.add_patch(
            plt.Circle(pts[m],
                       scale_circles*np.sqrt(circle),
                       color="LightGray",
                       zorder=2))

    ## Set labels for pts
    names = {tuple(Sj):mutation_names_sep.join(
        [f"{mutation_names[i]}"
         for i in range(M)
         if Sj[i] == 1])
             for Sj in S}


    names[M*(0,)] = 'normal'
    for x, name in names.items():
        ax.text((pts[x][0]+name_x_offsets[x]
                 if name_x_offsets is not None else
                 pts[x][0]),
                (pts[x][1]+name_y_offsets[x]
                 if name_y_offsets is not None else
                 pts[x][1]),
                (name + f"\n($n={circles[x]}$)" if include_n_circles
                 else name),
                # color=x if x != (1, 1, 1) else "Black",
                color="Black",
                ha='center',
                va='center',
                fontsize=(8 if multiplier_font_size == 1
                          else 6*multiplier_font_size))

    ## Draw a circle outside
    if M >= 3 and positions == "out_to_in":
        thetas = np.linspace(0, M*np.pi, 2000)
        ax.plot(M*np.cos(thetas),
                M*np.sin(thetas),
                color="Gray")

    # lims = None
    # if lims is None:
    #     if positions == "out_to_in":
    #         xlims = 1.01*np.array([-M, M])
    #         ylims = xlims
    #     else:
    #         xlims = np.array([-scale_circles*circles[M*(0,)]*3/2,
    #                           1*M + scale_circles*circles[M*(0,)]*3/2])
    #         ylims = 1.22*np.array([-M/2, M/2])
    # else:
    xlims = [-0.52, 4.1]
    ylims = [-2.2, 2.1]

    ax.set_xlim(xlims)
    ax.set_ylim(ylims)

    ## Clear ticks and spines
    ax.set_xticks([])
    ax.set_yticks([])
    ax.tick_params(labelcolor='w',
                   top='off',
                   bottom='off',
                   left='off',
                   right='off')
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')

    if subplot_label is not None:
        if subplot_label_ha == 'left':
            subplot_x_pos = 0
        elif subplot_label_ha == 'right':
            subplot_x_pos = 1
        else:
            raise Exception("Unknown 'subplot_label_ha'. "
                            "Available options "
                            "are 'left' or 'right'.")

        ax.text(0. if subplot_label_ha == 'left'
                else 1.05,
                1 if subplot_label_va == 'top'
                else 0,
                subplot_label,
                transform=ax.transAxes,
                va=subplot_label_va,
                ha=subplot_label_ha,
                fontsize=16*multiplier_font_size,
                weight='bold')


    fig.tight_layout()
    if plot_name is None:
        fig.savefig(os.path.join(location_figures,
                                 "landscape.png"),
                    transparent=True)
    else:
        fig.savefig(os.path.join(location_figures,
                                 f"{plot_name}.png"),
                    transparent=True)
    plt.close('all')

    return fig