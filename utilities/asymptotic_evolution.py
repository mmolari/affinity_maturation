import numpy as np
import copy
import matplotlib.pyplot as plt

import am_sim as ams
from am_sim.utils import Tsel_psurv


def asymptotic_pop_evo(det_pop, par, C_const):
    '''
    Evolve the population passed as argument, at constant Ag concentration and
    without enforcing B-selection or carrying capacity.
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


def L1_distr_distance(vp1, vp2, x):
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
    dx = x[2] - x[1]
    avg_eps = np.dot(vp, x) * dx
    sigma = np.sqrt(np.dot(np.square(x - avg_eps), vp) * dx)
    sigma_dist_left = (avg_eps - par['xlim_minus']) / sigma
    sigma_dist_right = (par['xlim_plus'] - avg_eps) / sigma
    return (sigma_dist_left > sigma_min_dist) & (sigma_dist_right > sigma_min_dist)


def check_phi_u_relative_update(phi_list, u_list, threshold):
    phi_relative_update = (phi_list[-1] - phi_list[-2]) / phi_list[-1]
    u_relative_update = (u_list[-1] - u_list[-2]) / u_list[-1]
    return (u_relative_update < threshold) & (phi_relative_update < threshold)


def check_shifted_distribution_L1_convergence(x, vp_list, threshold):
    L1_dist = L1_distr_distance(vp_list[-1], vp_list[-2], x)
    return L1_dist < threshold


def asymptotic_phi_and_u(C_const, par_to_copy, T_max,
                         sigma_min_dist=5, u_phi_relative_update=5e-5,
                         L1_norm_threshold=5e-5):

    par = copy.deepcopy(par_to_copy)

    # extend simulation domain and change discretization
    par['xlim_minus'] = -100
    par['xlim_plus'] = 50
#    par['dx'] = 0.01
    # initialize the distribution to a standard gaussian
    par['mu_i'] = 0
    par['sigma_i'] = 1

    # initialize population
    pop = ams.det_pop(par)

    # containers for results
    N, avg_eps, vp = [], [], []
    u, phi = [], []

    # binding energy distribution domain
    x = pop.x

    # append results
    N.append(pop.N_cells())
    avg_eps.append(pop.mean_en())
    vp.append(np.copy(pop.varphi))

    for t in range(T_max):

        # evolve one round
        asymptotic_pop_evo(pop, par, C_const)

        # append results
        N.append(pop.N_cells())
        avg_eps.append(pop.mean_en())
        vp.append(np.copy(pop.varphi))
        phi.append(np.log(N[-1] / N[-2]))
        u.append(avg_eps[-1] - avg_eps[-2])

        # check for convergence only after at least 10 rounds of evolution
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
            print(
                f'Successful convergence to the desired precision after {t} iterations')

            # plot to remove
            fig, ax = plt.subplots(2, 2, sharex=True, figsize=(5, 4))
            ax[0, 0].plot(phi)
            ax[0, 0].set_ylabel(r'$\phi$')
            ax[0, 0].axhline(0, c='k')
            ax[0, 1].plot(u)
            ax[0, 1].set_ylabel(r'$u$')
            ax[0, 1].axhline(0, c='k')
            ax[1, 0].plot(np.log(N))
            ax[1, 0].set_ylabel(r'$\log N$')
            ax[1, 1].plot(avg_eps)
            ax[1, 1].set_ylabel(r'$\langle \epsilon \rangle$')
            plt.tight_layout()
            plt.show()

            plt.plot(x - avg_eps[-1], vp[-1])
            plt.plot(x - avg_eps[-3], vp[-3])
            plt.plot(x - avg_eps[-5], vp[-5])
            plt.xlim(-10, 30)
            plt.show()

            return phi[-1], u[-1]
        elif not domain_check:
            # if borders are broken then print warning of not convergence
            print(
                'Warning: the distribution reached the borders of the simulation domain.')
            # and return None
            return None, None

    print('Warning: maximum number of iterations reached.')
    return phi[-1], u[-1]
