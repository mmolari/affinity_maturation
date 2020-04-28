import am_sim as ams


def simulate_scheme_1(sim_type, par, D):
    # function that simulates the first immunization scheme for a given dosage
    rp = ams.simulate_immscheme(sim_type=sim_type, D_inj=[D],
                                T_delay=[28], meas_prot='4d', par=par)
    return rp


def simulate_scheme_2(sim_type, par, D):
    # function that simulates the second immunization scheme for a given dosage
    rp = ams.simulate_immscheme(sim_type=sim_type, D_inj=[D, D],
                                T_delay=[28, 28], meas_prot='1d', par=par)
    return rp


def simulate_scheme_3(sim_type, par, T):
    # function that simulates the third immunization scheme for a given time
    # delay
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
