import numpy as np

stochpop_color = 'C2'


def plot_error_area(ax, stoch_pop_avgs_list, N_GCs, x):
    means, stds = [], []
    for avgs in stoch_pop_avgs_list:
        means.append(np.mean(avgs))
        stds.append(np.std(avgs) / np.sqrt(N_GCs))
    means, stds = np.array(means), np.array(stds)
    ax.fill_between(x, means - stds, means + stds,
                    color=stochpop_color, alpha=0.3, label='stoch. sim')
