import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import am_sim as ams


def stoch_GC_clonotype_evo(par, D, T):
    '''
    This function is used to simulate the evolution of a stochastic germinal
    center, keeping track of the clonal families. It takes as input the model
    parameters, the injected Ag dosage and the total simulation time. It
    returns a database containing the number of clone for each clonal family
    and each evolution time, and the binding energy of the founder clones for
    each family.

    Args:
    - par: dictionary of model parameters, used for the simulation
    - D (float): injected Ag dosage in micrograms
    - T (int): duration of the simulation, in days from injection.

    Returns:
    - ctype_count_df (pandas dataframe): dataframe whose index is days, columns
        are clonal families and entries are number of clones per family at each
        evolution round.
    - en_founders (array of floats): binding energy of the founder clone for
        each clonal family.
    '''
    # convert test dosage into concentration
    D_test = 1.
    C_test = D_test / par['alpha_C']
    # initialize GC
    GC = ams.germinal_center(GC_type='stochastic', C_inj=C_test, par=par)
    # add clonotype tracing to GC stochastic population
    GC.pop = ams.stoch_pop(par=par, trace_clonotypes=True)

    # initial energy of founder clones
    en_founders = GC.pop.en_founders
    N_clones = en_founders.size

    # evolve GC and save time and clonotype array
    t_days, ctype_list = [], []
    for t in range(int(T // par['days_per_turn'] + 1)):

        if GC.t_days >= par['T_GC_formation_days']:
            # at formation time start saving the results at every round
            t_days.append(GC.t_days)
            ctype_list.append(np.copy(GC.pop.clonotype))

        _ = GC.evolve_one_round()

        if GC.state == 'extinct' or GC.t_days > T:
            break

    # create matrix with clonotype count, and put it in a pandas dataframe
    # line index is the time, and columns are the relative abundance of clonotypes
    clt_count = np.zeros((len(t_days), N_clones), dtype=np.uint)
    for nt, clt in enumerate(ctype_list):
        clt_id, clt_n = np.unique(clt, return_counts=True)
        clt_count[nt, clt_id] += clt_n.astype(np.uint)

    # build dataframe containing clonotype numbers per round
    ctype_count_df = pd.DataFrame(data=clt_count, index=t_days,
                                  columns=np.arange(N_clones))

    return ctype_count_df, en_founders


#  --- plotting functions

def ctype_cmap(v_min, v_max, en_founders):
    '''
    This function assigns a color to each clonal family based on the energy of
    its founder clones 'en_founders'. The color scale varies between v_min and
    v_max, and a ScalarMappable object containing details on the colorbar is
    returned.

    Args:
    - v_min, v_max (floats): minimum and maximum binding energy of the
        colorscale.
    - en_founders (list): list of founders binding energies.

    Returns:
    - colors (list): list of colors assigned to each founder clone.
    - mapp (ScalarMappable object): ScalarMappable object, can be used to build
        a colorbar.
    '''
    # pick the colormap
    cmap = plt.get_cmap('jet')
    # create a normalization object and use it to assign colors to energies
    norm = mpl.colors.Normalize(vmin=v_min, vmax=v_max)
    colors = [cmap(norm(ei)) for ei in en_founders]
    # create and return a ScalarMappable object
    mapp = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    mapp.set_array(en_founders)
    return colors, mapp


def add_text_max_frequency(ctype_count_df, ax):
    '''
    This function adds a text on the axis specifying the final fraction of the
    population represented by the most expanded clonal family. It takes as
    argument the axis on which to plot and the dataframe of clone counts
    over time.

    Args:
    - ctype_count_df (pandas DataFrame object): pandas dataframe containing the
        clone count over time.
    - ax (ax object): ax on which the text is added.
    '''
    # evaluate the final population fraction
    ctype_n_final = ctype_count_df.to_numpy()[-1]
    f_max = ctype_n_final.max() / ctype_n_final.sum()
    # print the results on the axis
    ax.text(1., 1.07,
            f'final fraction of most expanded clonal family = {f_max:.2}',
            horizontalalignment='right', verticalalignment='top',
            transform=ax.transAxes)


def plot_GC_ctype_evo(ax, ctype_count_df, en_founders, cmap_range):
    '''
    This function plots the evolution of the clonal families over time for a
    stochastic GC simulation. The families are represented in different colors,
    depending on the binding energy of their founder clone. A ScalarMappable
    object is returned to later plot the colorbar used to represent the
    energies.

    Args:
    - ax (ax object): ax on which to plot the stochastic evolution of the
        clonal families.
    - ctype_count_df (pandas dataframe): dataframe whose index is days, columns
        are clonal families and entries are number of clones per family at each
        evolution round.
    - en_founders (list of floats): binding energy of the founder clones for
        each clonal family.
    - cmap_range (list of two floats): this list, in the format [vmin, vmax],
        is used to assign the lower and upper limit for the colorbar

    Returns:
    - mapp (ScalarMappable object): this object can be used to plot the
        colorbar used to represent the binding energy of each clonal family.
    '''
    v_min, v_max = cmap_range
    colors, mapp = ctype_cmap(v_min, v_max, en_founders)
    ctype_count_df.plot(kind='area', legend=[],
                        color=colors, linewidth=0, ax=ax)
    add_text_max_frequency(ctype_count_df, ax)
    ax.set_xlabel('time (days)')
    ax.set_ylabel('clonal families size')
    t_min, t_max = ctype_count_df.index.min(), ctype_count_df.index.max()
    ax.set_xlim(t_min, t_max)
    N_t = ctype_count_df.sum(axis=1)
    ax.set_xticks(np.arange(t_min, t_max + 1, 2).astype(int))
    return mapp
