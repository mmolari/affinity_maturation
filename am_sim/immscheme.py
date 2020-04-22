import numpy as np

from .utils import evaluate_responders
from .pop_stochastic import stoch_pop
from .pop_deterministic import det_pop
from .germinal_center import germinal_center


def simulate_immscheme(sim_type, D_inj, T_delay, meas_prot, par):
    '''
    This function simulates an immunization scheme and returns the elicited
    population of responder cells.

    Args:
    - sim_type (string): either 'stochastic' or 'deterministic'
    - D_inj (float array): list of injected dosages
    - T_delay (int array): list of time delays between injections and
        before measurement (in days).
    - meas_prot (string): measurement protocol. Either '1d' or '4d' according
        to weather the measurement is performed 1 day after boost or 4 days
        after injection. This controls the fraction of MCs in the responder
        population, according to model parameters.
    - par (dict): model parameters dictionary

    Returns:
    - resp_pop (stoch_pop/det_pop object): resulting population of responder
        cells. It is either a stoch_pop or det_pop object, depending on the
        simulation type.
    '''

    # --- immunization scheme setup ---

    # check that the correct number of dosage and delays has been provided
    assert len(D_inj) == len(T_delay), 'Error: the number of injections and\
    delays is not the same'

    # conversion Dose -> Concentration
    C_inj = np.array(D_inj) / par['alpha_C']
    T_cum = np.cumsum(T_delay)
    T_inj = np.concatenate([[0], T_cum[:-1]])  # injection times
    t_end = T_cum[-1]  # end of simulation, prior to measurement
    N_inj = len(D_inj)  # number of injections

    # define g (responder memory fraction) based on the measurement
    # protocol
    if meas_prot == '1d':
        g_mem = par['g_1d']
    elif meas_prot == '4d':
        g_mem = par['g_4d']
    else:
        raise Exception('meas_prot must be either 1d or 4d')

    # either stochastic or deterministic simulation
    if sim_type == 'stochastic':
        pop_class = stoch_pop
    elif sim_type == 'deterministic':
        pop_class = det_pop
    else:
        raise Exception(
            'simulation sim_type must be either stochastic or deterministic')

    # if all dosages are zero then return naive population
    if np.all(C_inj == 0):
        resp_pop = pop_class(par)
        return resp_pop
    elif np.any(C_inj == 0):
        raise Exception('warning: one of the dosages is zero. This case is not\
        defined in the model')

    # instantiate empty MC / PC populations
    MC_tot, PC_tot = pop_class.create_empty(par), pop_class.create_empty(par)

    # initialize first GC performing first injection
    GC_1 = germinal_center(GC_type=sim_type, par=par, C_inj=C_inj[0], t_inj=0)
    GC_list = [GC_1]
    # parameters for next injection
    next_inj_id = 1
    next_inj_T = T_inj[next_inj_id]
    next_inj_C = C_inj[next_inj_id]

    # --- simulate immunization scheme ---

    while GC_1.t_days <= t_end:

        # at the given time perform new injection and add new GC:
        if next_inj_T is not None and GC_1.t_days == next_inj_T:
            # perform new injection and create new GC
            new_GC = germinal_center(
                GC_type=sim_type, par=par, C_inj=next_inj_C,
                t_inj=next_inj_T, mc_seed=MC_tot)
            # add it to the list of existing GCs
            GC_list.append(new_GC)
            # setup next injection if there is one left to do
            next_inj_id += 1
            if next_inj_id < N_inj:
                next_inj_T = C_inj[next_inj_id]
                next_inj_C = T_inj[next_inj_id]
            else:
                # last injection performed, no new injections scheduled
                next_inj_T = None
                next_inj_C = None

        # evolve all GCs and add MC/PC outcome to the total population
        for GC in GC_list:
            MC_gc, PC_gc = GC.evolve_one_round()
            MC_tot.merge_with(MC_gc)
            PC_tot.merge_with(PC_gc)

    # merge MCs and PCs into responder populations and return it
    resp_pop = evaluate_responders(MC_tot, PC_tot,
                                   g_mem, sim_type,
                                   N_res=par['N_i'])
    return resp_pop


def simulate_immscheme_from_dset(sim_type, par, dset):
    # TODO: implement this wrapper
    pass
