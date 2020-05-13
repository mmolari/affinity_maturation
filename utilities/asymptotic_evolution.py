import numpy as np
import copy
import matplotlib.pyplot as plt

import am_sim as ams
from am_sim.utils import Tsel_psurv


def asymptotic_pop_evo(det_pop, par, C_const):
    '''
    This function performs one evolution round on the population passed as
    argument, at constant Ag concentration and without enforcing B-selection or
    carrying capacity.

    Args:
    - det_pop (det_pop object): the population to evolve, as a det_pop object.
    - par (dict): the parameters dictionary
    - C_const (float): the Ag concentration
    '''
    # --- evolve population
    # differentiation (simply remove a fraction par['diff_prob'] of the population)
    det_pop.N *= (1. - par['diff_prob'])

    # duplication and mutation
    det_pop.expand(par)

    # T selection
    en = det_pop.energies()
    bareps = det_pop.bareps()
    psurv_T = Tsel_psurv(en=en, bareps=bareps, C=C_const, par=par)
    det_pop.select_with_psurv(psurv_T)


def evolve_pop_at_constant_C(C_const, par, T_skip, T_save):
    '''
    This function simulates the evolution of a population under constant Ag
    concentration for a specified number of evolution rounds.

    Args:
    - C_const (float): Ag concentration.
    - par (dict): model parameters dictionary
    - T_skip (int): number of pre-evoluiont rounds to be simulated but not
        included in the results dictionary.
    - T_save (int): number of evolution rounds to be simulated and included
        in the results dictionary.

    Returns:
    - results (dict): a dictionary containing the results of the simulation.
        Its keys are:
        - 't' : list of times at which the quantities are referred,
            in rounds of evolution.
        - 'N' : list of population sizes for each time.
        - 'avg_eps' : list of average population energy for each time.
        - 'distr_y' : list of binding energy distributions of the population
            for each time.
        - 'distr_x' : energy domain of the distributions.
    '''
    # initialize population
    pop = ams.det_pop(par)

    # evolve some rounds without saving the results to get closer to asymptotic state
    for t in range(T_skip):
        asymptotic_pop_evo(pop, par, C_const)

    # containers for results
    N, avg_eps, distr_y = [], [], []
    t_arr = np.arange(T_save)

    for t in t_arr:
        # append results
        N.append(pop.N_cells())
        avg_eps.append(pop.mean_en())
        distr_y.append(np.copy(pop.varphi))

        # evolve one round
        asymptotic_pop_evo(pop, par, C_const)

    # pack results in a dictionary
    results = {
        't': t_arr,
        'N': np.array(N),
        'avg_eps': np.array(avg_eps),
        'distr_y': np.array(distr_y),
        'distr_x': pop.energies(),
    }

    return results


def L1_recentered_distr_distance(vp1, vp2, x):
    '''
    This utility function evaluates the L1 distance between two distributions,
    recentering them first around their mean.

    Args:
    - vp1, vp2 (lists of float): y-values of the two distributions.
    - x (list of float): common domain of the two distributions.

    Returns:
    - L1_norm (float): L1 distance between the two recentered distributions.
    '''
    # evaluate the means
    dx = x[1] - x[0]
    mean_1 = np.dot(vp1, x) * dx
    mean_2 = np.dot(vp2, x) * dx
    # evaluate interpolation with mean re-centered in zero
    vp1_interp = np.interp(x, x - mean_1, vp1, left=0, right=0)
    vp1_interp /= np.sum(vp1_interp) * dx  #  normalize
    vp2_interp = np.interp(x, x - mean_2, vp2, left=0, right=0)
    vp2_interp /= np.sum(vp2_interp) * dx  #  normalize
    # L1 distance between the re-centered distributions
    L1_norm = np.sum(np.abs(vp2_interp - vp1_interp)) * dx
    return L1_norm


def check_domain_distance(x, vp, par, sigma_min_dist):
    '''
    This utility function checks if the distribution mean is at least a number
    'sigma_min_dist' of standard deviations away from the simulation domain
    boundaries.

    Args:
    - x (list of float): distribution domain
    - vp (list of float): distribution values
    - par (dict): the model parameters dictionary
    - sigma_min_dist (float): minimum distance from the simulation boundaries
        in units of the distribution's standard deviation

    Returns:
    - check (bool): the result of the check.
    '''
    dx = x[2] - x[1]
    # mean of the distribution
    avg_eps = np.dot(vp, x) * dx
    # standard deviation of the distribution
    sigma = np.sqrt(np.dot(np.square(x - avg_eps), vp) * dx)
    # distance of the distribution from right and left domain boundaries, in
    # units of the standard deviation
    sigma_dist_left = (avg_eps - par['xlim_minus']) / sigma
    sigma_dist_right = (par['xlim_plus'] - avg_eps) / sigma
    check = (sigma_dist_left > sigma_min_dist) & (
        sigma_dist_right > sigma_min_dist)
    return check


def check_phi_u_relative_update(phi_list, u_list, threshold):
    '''
    This utility function checks if the relative update of u and phi is below
    a specified threshold.

    Args:
    - phi_list (list of float): list of values of phi, the last two are used
        for the check
    - u_list (list of float): list of values of u, the last two are used for
        the check
    - threshold (float): the specified threshold

    Returns:
    - check (bool): the result of the check.
    '''
    # relative update of u and phi
    phi_relative_update = (phi_list[-1] - phi_list[-2]) / phi_list[-1]
    u_relative_update = (u_list[-1] - u_list[-2]) / u_list[-1]
    check = (u_relative_update < threshold) & (phi_relative_update < threshold)
    return check


def check_shifted_distribution_L1_convergence(x, vp_list, threshold):
    '''
    This utility function checks if the L1 distance of the distributions,
    recentered around their mean, is below a given threshold.

    Args:
    - x (list of float): distributions domain
    - vp_list (list): list of distributions, the last two are used for the check
    - threshold (float): the specified threshold

    Returns:
    - check (bool): the result of the check.
    '''
    L1_dist = L1_recentered_distr_distance(vp_list[-1], vp_list[-2], x)
    check = L1_dist < threshold
    return check


def asymptotic_phi_and_u(C_const, par_to_copy, T_max,
                         sigma_min_dist=5, u_phi_relative_update=5e-5,
                         L1_norm_threshold=5e-5, xlim_m=-100, xlim_p=50,
                         dx=0.01):
    '''
    This function evaluates the asymptotic population growth rate 'phi' and
    asymptotic wave speed 'u' for a given value of the Ag concentration. The
    values are found by simulating the population evolution under constant Ag
    concentration and removing the carrying capacity constraint and Ag-binding
    selection. Convergence is reached when three conditions are met: the L1
    distance between two successive updates of the recentered binding energy
    distribution is under a specified theshold, the relative update of u and
    phi is below a specified threshold, and the distribution has a minimum
    distance from the simulation domain boundaries. If these conditions are not
    met the function returns the value of phi and u obtained after the maximum
    allowed number of iterations, together with a warning message. If the
    distribution does not respect the minimum distance from the simulation
    boundaries then the function returns None values and a warning message.

    Args:
    - C_const (float): the value of the Ag concentration.
    - par_to_copy (dict): dictionary of model parameters.
    - T_max (int): maximum allowed number of iterations of the algorithm.
    - sigma_min_dist (float, optional): minimum required distance of the
        distribution average value from the simulation domain boundaries, in
        units of the distribution's standard deviation. If this minimum
        distance is not respected the function stops and returns None values.
        If not specified is set to 5.
    - u_phi_relative_update (float, optional): maximum allowed value of the
        relative update of u and phi between two rounds of evolution that is
        compatible with convergence. Default value is 5e-5.
    - L1_norm_threshold (float, optional): maximum allowed distance between
        binding energy distributions before and after application of the
        evolution operator, that is compatible with convergence. The distance
        is evaluated as L1 norm between the two distributions recentered around
        their mean. Default value is 5e-5.
    - xlim_m (float, optional): minimum energy simulation boundary domain.
        Defalut value is -100.
    - xlim_p (float, optional): maximum energy simulation boundary domain.
        Defalut value is 50.
    - dx (float, optional): simulation domain discretization step. Defalut
        value is 0.01.

    Returns:
    - phi, u (floats): asymptotic values of the growth rate and wave speed.
    '''

    par = copy.deepcopy(par_to_copy)

    # extend simulation domain and change discretization
    par['xlim_minus'] = xlim_m
    par['xlim_plus'] = xlim_p
    par['dx'] = dx
    # initialize the distribution to a standard gaussian
    par['mu_i'] = 0
    par['sigma_i'] = 1

    # initialize population
    pop = ams.det_pop(par)
    pop.N = 1.  #  set the population size to one
    # add an attribute to the population: the logarithm of N
    # this helps keeping track of the size of the population under exponential
    # expansion
    pop.logN = np.log(pop.N)

    # containers for results
    logN, avg_eps, vp = [], [], []
    u, phi = [], []

    # binding energy distribution domain
    x = pop.x

    # append results
    logN.append(pop.logN)
    avg_eps.append(pop.mean_en())
    vp.append(np.copy(pop.varphi))

    # iteratively simulate evolution rounds until convergence is reached
    for t in range(T_max):

        # evolve one round
        asymptotic_pop_evo(pop, par, C_const)
        # update logN and reset N
        pop.logN += np.log(pop.N)
        pop.N = 1.

        # append results
        logN.append(pop.logN)
        avg_eps.append(pop.mean_en())
        vp.append(np.copy(pop.varphi))
        phi.append(logN[-1] - logN[-2])
        u.append(avg_eps[-1] - avg_eps[-2])

        # check for convergence only after at least 100 rounds of evolution
        if t < 100:
            continue

        # check the three conditions
        domain_check = check_domain_distance(x, vp[-1], par, sigma_min_dist)
        phi_u_check = check_phi_u_relative_update(
            phi, u, u_phi_relative_update)
        L1_norm_check = check_shifted_distribution_L1_convergence(
            x, vp, L1_norm_threshold)

        # if all are satisfied then break and return the results
        if domain_check and phi_u_check and L1_norm_check:
            print('Successful convergence to the desired precision ' +
                  f'after {t} iterations')

            return phi[-1], u[-1]

        elif not domain_check:
            # if domain borders are reached then print warning of not convergence
            print('Warning: the distribution reached the borders of' +
                  f' the simulation domain at iteration t = {t}.')

            # and return None
            return None, None

    # if the simulation reached the maximum number of iterations print a warining
    # and return the results
    print('Warning: maximum number of iterations reached.')
    return phi[-1], u[-1]
