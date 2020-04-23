import am_sim as ams


def simulate_immscheme_from_dset(sim_type, par, dset):
    '''
    Wrapper function to execute simulate a given immunization scheme according
    to a dataset.

    Args:
    - sim_type (string): either 'stochastic' or 'deterministic'
    - dset (dataset object): dataset object containing the immunization scheme
         to simulate.
    - par (dict): model parameters dictionary

    Returns:
    - resp_pop (stoch_pop/det_pop object): resulting population of responder
        cells. It is either a stoch_pop or det_pop object, depending on the
        simulation type.
    '''
    D_inj, T_delay, meas_prot = dset.D_inj, dset.T_delay, dset.meas_prot
    resp_pop = ams.simulate_immscheme(sim_type, D_inj, T_delay, meas_prot, par)
    return resp_pop
