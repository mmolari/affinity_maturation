import numpy as np

from am_sim.model_parameters import high_en_exp_cutoff

stochpop_color = 'C2'


def plot_error_area(ax, stoch_pop_avgs_list, N_GCs, x):
    '''
    This function plots on the specified axis the confidence interval for
    stochastic simulations. The interval is plotted as a green shaded area,
    spanning one standard deviation above and below the mean. The standard
    deviation is corrected with the square root of the number of GCs that we
    expect to measure.

    Args:
    - ax: ax object on which to plot.
    - stoch_pop_avgs_list (list of lists of float): list of list of
        measurements of the quantity, which can be either average binding
        energy or high-affinity fraction.
    - N_GCs (int): number of GCs we expect to measure on each mouse.
    - x (list of floats): x-domain of the plot. Same dimension as the number of
        lists in stoch_pop_avgs_list.
    '''
    # containers for means and standard deviations
    means, stds = [], []
    # cycle through all the stochastic simulation sets
    for avgs in stoch_pop_avgs_list:
        # for every set evaluate the mean and standard deviation
        means.append(np.mean(avgs))
        # the standard deviation is normalized by the square root of the
        # number of GCs considered
        stds.append(np.std(avgs) / np.sqrt(N_GCs))
    # transform the lists in numpy arrays.
    means, stds = np.array(means), np.array(stds)
    # color the area one standard deviation above and below the mean
    ax.fill_between(x, means - stds, means + stds,
                    color=stochpop_color, alpha=0.3, label='stoch. sim')


def plot_stoch_avg_energy_distr_exp_norm(ax, stoch_pops, bins):
    '''
    Plots the average binding energy distribution for a series of stochastic
    simulations of the same immunization scheme. For good comparison with
    experimental measurments the distribution is normalzied taking into account
    the experimental cutoff limit.

    Args:
    - ax: ax object on which to plot.
    - stoch_pops: list of stochastic populations resulting from the simulations
        of the same immunization scheme.
    - bins: binning for the histogram.
    '''
    bins = np.array(bins)
    # number of bins and of simulations
    B, N = bins.size - 1, len(stoch_pops)
    # center of bins
    x_bins = 0.5 * (bins[:-1] + bins[1:])
    # discretization interval
    dx = bins[1] - bins[0]
    # initialize histogram weights
    weights = np.zeros(B)
    for sp in stoch_pops:
        hist, _ = np.histogram(sp.en, bins=bins, density=True)
        weights += hist
    # normalize by the part below the detection threshold
    norm = np.sum(weights[x_bins <= high_en_exp_cutoff]) * dx
    weights /= norm
    ax.hist(x_bins, weights=weights, bins=bins, histtype='step',
            label='stoch. sim.', color=stochpop_color)
