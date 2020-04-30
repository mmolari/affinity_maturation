import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import am_sim as ams


def stoch_GC_clonotype_evo(par, D, T):
    '''
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
    # frequency of most abundant clone over time
    f_max = clt_count.max(axis=1) / clt_count.sum(axis=1)

    # build dataframe containing clonotype numbers per round
    ctype_evo_df = pd.DataFrame(data=clt_count, index=t_days,
                                columns=np.arange(N_clones))

    return ctype_evo_df, en_founders, f_max, t_days


def plot_GC_ctype_evo(df, en_founders, savename=None, cmap_range=None):
    '''
    '''
    cmap = plt.get_cmap('jet')
    if cmap_range is None:
        vmin, vmax = en_founders.min(), en_founders.max()
    else:
        vmin, vmax = cmap_range
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    mapp = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    mapp.set_array(en_founders)
    c = [cmap(norm(ei)) for ei in en_founders]
    df.plot(kind='area', legend=[], figsize=(10, 5), color=c, linewidth=0)
    plt.colorbar(mapp, label='clonotype initial binding energy')
    f_max = df.to_numpy()[-1].max() / df.to_numpy()[-1].sum()
    plt.text(
        0.6, 0.95, f'final max frequency = {f_max:.2}', transform=plt.gca().transAxes)
    plt.xlabel('time (days)')
    plt.ylabel('clonotype abundance')
    if savename is not None:
        plt.savefig(savename)
    plt.show()
