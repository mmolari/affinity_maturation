import numpy as np
import copy

from am_sim.model_parameters import high_en_exp_cutoff

detpop_color = 'C0'


def plot_exp_normalized_pop(ax, det_pf):
    '''
    This function plots the binding energy distribution of the population
    function on the specified axis. The distribution is normalized only for on
    the part below the high-energy experimental detection limit, so that it can
    be compared with experimental measurements.

    Args:
    - ax (matplotlib ax object): ax on which the binding energy distribution
        is plotted.
    - det_pf (det_pop object): det_pop whose binding energy distribution must
        be plotted.
    '''
    # domain and discretization step
    x, dx = det_pf.x, det_pf.dx
    # distribution
    y = det_pf.varphi
    # mask to select the range of energies lower than the experimental cutoff
    norm_mask = x <= high_en_exp_cutoff
    # normalization factor for the part below-cutoff
    norm = np.sum(y[norm_mask]) * dx
    # normalized distribution
    y_norm = y / norm
    # plotting function
    ax.plot(x, y_norm, color=detpop_color, label='theory')
