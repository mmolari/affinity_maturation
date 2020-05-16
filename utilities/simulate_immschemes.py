import am_sim as ams


def simulate_scheme_1(sim_type, par, D):
    '''
    Function to simulate the first immunization scheme at a specified Ag dosage.
    It returns the resulting population function.

    Args:
    - sim_type (str): either 'stochastic' or 'deterministic', depending on the
        desired simulation type.
    - par (dict): dictionary of model parameters.
    - D (float): injected Ag dosage.
    '''
    # function that simulates the first immunization scheme for a given dosage
    rp = ams.simulate_immscheme(sim_type=sim_type, D_inj=[D],
                                T_delay=[28], meas_prot='4d', par=par)
    return rp


def simulate_scheme_2(sim_type, par, D):
    '''
    Function to simulate the second immunization scheme at a specified Ag dosage.
    It returns the resulting population function.

    Args:
    - sim_type (str): either 'stochastic' or 'deterministic', depending on the
        desired simulation type.
    - par (dict): dictionary of model parameters.
    - D (float): injected Ag dosage.
    '''
    rp = ams.simulate_immscheme(sim_type=sim_type, D_inj=[D, D],
                                T_delay=[28, 28], meas_prot='1d', par=par)
    return rp


def simulate_scheme_3(sim_type, par, T):
    '''
    Function to simulate the second immunization scheme at a specified Ag dosage.
    It returns the resulting population function.

    Args:
    - sim_type (str): either 'stochastic' or 'deterministic', depending on the
        desired simulation type.
    - par (dict): dictionary of model parameters.
    - T (int): time delay between injections in days.
    '''
    rp = ams.simulate_immscheme(sim_type=sim_type, D_inj=[10, 10],
                                T_delay=[T, 28], meas_prot='1d', par=par)
    return rp


# dictionary of simulation functions, relates the number of immunization scheme
# to the corresponding simulation function
simulate_scheme_n = {
    1: simulate_scheme_1,
    2: simulate_scheme_2,
    3: simulate_scheme_3,
}


def simulate_stoch_scheme_N_times(scheme_n, N_sims, par, variable):
    '''
    Function that performs the stochastic simulation of immunization scheme
    'scheme_n' a number 'N_sims' of times. It takes as arguments the simulation
    parameters dictionary 'par' and the variable of the immunization scheme
    considered (Ag dosage in micrograms for scheme 1 and 2, and injection time
    delay in days for scheme 3). It returns a list of responders population.
    One per simulation performed

    Args:
    - scheme_n (int): number of the scheme to simulate (1,2 or 3).
    - N_sims (int): number of independent stochastic simulations.
    - par (dict): dictionary of model parameters
    - variable (int or float): Variable specifying the condition under which
        the scheme is tested. This corresponds to Ag dosage for scheme 1 or 2,
        or the time delay between injection for scheme 3.

    Returns:
    - rp_list (list of stoch_pop objects): list of responder populations, one
        per simulation performed.
    '''
    # simulation function for the scheme selected
    sim_f = simulate_scheme_n[scheme_n]
    # container of responder populations
    rp_list = []
    for n_sim in range(N_sims):
        # simulate the scheme and append the results
        rp = sim_f('stochastic', par, variable)
        rp_list.append(rp)
    # return the results list
    return rp_list
