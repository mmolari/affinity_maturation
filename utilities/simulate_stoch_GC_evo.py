import numpy as np
import copy

import am_sim as ams


def initialize_GC_with_stoch_pop(D_inj, par, st_pop):
    '''
    Given the injected dosage, the parameters set and a stochastic population,
    this function initializes a stochastic GC having as initial population the
    specified stochastic population.

    Args:
    - D_inj (float): injected dosage in micrograms of Ag
    - par (dict): model parameters dictionary
    - st_pop (stoch_pop object): initial population of the GC

    Returns:
    - GC (germinal_center object): GC initialized with the given parameters and
        injected concentration, and having the stochastic population as initial
        population.
    '''
    # conversion Dose -> Concentration
    C_inj = D_inj / par['alpha_C']

    # create stochastic GC
    GC = ams.germinal_center(GC_type='stochastic',
                             par=par, C_inj=C_inj, t_inj=0)

    # overwrite stochastic population with copy
    GC.pop = copy.deepcopy(st_pop)

    return GC


def evolve_GC_till_extinction(par, GC, T_max=1000):
    '''
    This function performs the evolution of the specified GC untill its
    extinction, and returns information on the Ag evolution, the evolution of
    the B-cell population and the final generated MC/PC population. Results are
    returned in the form of a dictionary.

    Args:
    - par (dict): model parameters dictionary
    - GC (germinal_center object): germinal_center object to evolve until
        extinction. It must be a newly initialized object.
    - T_max (int, optional): sets a maximum number of evolution rounds, after
        which the evolution is stopped even if the GC is not extinct. Set to
        1000 by default.

    Returns:
    - results (dictionary): dictionary containing info on the evolution. Its
        entries are:
        - 'ag_t' : time of Ag evolution in days
        - 'ag_C' : evolution of Ag concentration, relative to 'ag_t' time
        - 'GC_t': time of GC evolution
        - 'GC_pop_en': population binding energies for each evolution step.
            Relative to 'GC_t' time.
        - 'GC_avg_eps': average binding energy of the population. Relative to
            'GC_t' time.
        - 'GC_N': Population size. Relative to 'GC_t' time.
        - 'MC_en': binding energies of the final MC population
        - 'PC_en': binding energies of the final PC population
    '''

    # MC and PC populations containers
    MC, PC = ams.stoch_pop.create_empty(par), ams.stoch_pop.create_empty(par)

    # results containers
    save_t = []
    save_pop_en = []
    save_avg_eps = []
    save_N = []

    # maximum number of simulated rounds
    for t in range(T_max):
        if GC.state == 'mature':
            # if the GC is mature save results
            save_t.append(GC.t_days)
            save_pop_en.append(np.copy(GC.pop.energies()))
            save_avg_eps.append(GC.pop.mean_en())
            save_N.append(GC.pop.N_cells())
        elif GC.state == 'extinct':
            # end simulation when GC is extinct
            break

        # evolve GC
        MC_gc, PC_gc = GC.evolve_one_round()
        # collect MCs and PCs
        MC.merge_with(MC_gc)
        PC.merge_with(PC_gc)

    # assemble results in a dictionary
    results = {
        # antigen evolution
        'ag_t': np.array(GC.ag.history['t']),
        'ag_C': np.array(GC.ag.history['C']),
        # GC evolution
        'GC_t': np.array(save_t),
        'GC_pop_en': np.array(save_pop_en),
        'GC_avg_eps': np.array(save_avg_eps),
        'GC_N': np.array(save_N),
        # MC / PC populations
        'MC_en': MC.en,
        'PC_en': PC.en,

    }

    return results
