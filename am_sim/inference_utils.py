import numpy as np
import os
import pickle as pkl
import copy

from .utils import resize_to_exp_limits_det
from .immscheme import simulate_immscheme


def dset_logl(dset, det_pf):
    '''
    Gives the log-likelihood of the data for a given population function. The
    experimental sensitivity limits are taken into account inside the function.

    Args:
    - dset (dataset object): the experimental dataset
    - det_pf (det_pop object): deterministic population function, corresponding
        to the model predicted responders pop for this immunization scheme.

    Returns:
    - logl (float): total log-likelihood of the data given the model.
    '''
    # resize population function to experimental sensitivity limits
    x, dx, vp = resize_to_exp_limits_det(det_pf)
    # interpolate restricted varphi on experimental measurements
    exp_en = dset.all_en()
    likl = np.interp(exp_en, x, vp, left=0, right=0)
    logl = np.log(likl)
    # return total log likelihood
    return np.sum(logl)


def responders_from_dset(dset, par, sim_type):
    '''
    Given a dataset and a set of parameters it returns the responders
    population obtained from simulation of the specified scheme.

    Args:
    - dset (dataset object): the dataset containing the details of the
        immunization scheme to simulate
    - par: model parameters dictionary
    - sim_type (str): either 'stochastic' or 'deterministic'

    Returns:
    - resp_pop (stoch_pop/det_pop object): responders population. Its type
        depends on the value of sim_type argument
    '''
    # extract immunization scheme parameters from the dataset
    D_inj, T_delay, meas_prot = dset.D_inj, dset.T_delay, dset.meas_prot
    # simulate the protocol and return the responders population
    resp_pop = simulate_immscheme(sim_type, D_inj, T_delay, meas_prot, par)
    return resp_pop


def dset_list_logl(par, dset_list):
    '''
    For a given list of datasets this function evaluates the total
    log-likelihood of the parameter set. This is defined by evaluating for each
    dataset the model prediction for the responder population, and summing up
    the log-likelihood of the data given the prediction.

    Args:
    - par: model parameters dictionary
    - dset_list (list of dataset objects): list of all the datasets considered

    Returns:
    - tot_logl (float): total log-likelihood
    '''
    tot_logl = 0
    # for each dataset in the list
    for dset in dset_list:
        # evalutate the model prediction for the responders population
        resp_pop = responders_from_dset(dset, par, sim_type='deterministic')
        # evaluate the log-likelihood of the data w.r.t. the prediction
        tot_logl += dset_logl(dset, resp_pop)
    # return the total log-likelihood
    return tot_logl


def init_search_directory(dir_path):
    '''
    This function is used to initialize the directiory in which to save the
    results of the likelihood-maximization procedure. The function check that
    the directory is empty or non-existent (in which case it creates it)

    Args:
    - dir_path (str): path of the directory
    '''
    # check if the folder already exists
    if os.path.exists(dir_path):
        # check if it containts every file
        if not len(os.listdir(dir_path)) == 0:
            # if it contains files raise exceptions to avoid overwriting
            raise Exception(
                'To avoid overwriting the user should provide an empty or' +
                ' non-existent folder. Insted the specified folder exists' +
                f' and contains the following files:\n {os.listdir(dir_path)}')
    else:
        # if the folder does not exists create it
        os.makedirs(dir_path, exist_ok=False)


def save_search_initial_setup(pt):
    '''
    '''
    # save initial parameters
    with open(os.path.join(pt.save_folder, 'par_i.pkl'), 'wb') as f:
        pkl.dump(pt.pars[0], f)
        f.close()

    # save dataset
    with open(os.path.join(pt.save_folder, 'dataset_list.pkl'), 'wb') as f:
        pkl.dump(pt.dsets, f)
        f.close()

    # save search setup
    with open(os.path.join(pt.save_folder, 'search_setup.txt'), 'w') as f:
        f.write('params to vary:\n' + str(pt.pars_to_mutate) + '\n')
        f.write(f'n rounds = {pt.T_max}\n')
        f.write(f'n layers = {pt.n_layers}\n')
        f.write(f'inverse temperatures = {pt.betas}\n')
        f.write(f'mutate one param. at a time = {pt.mut_sing}\n')
        f.write(f'mut strength = {pt.mut_str}\n')
        f.write(f'initial parameters total logl = {pt.logls[0]}\n\n')
        f.write('parameters initial value:\n----------------------\n')
        for key in pt.pars[0]:
            f.write(f'{key:<30} - {pt.pars[0][key]}\n')
        f.close()


def generate_variated_par(par, keys_to_mut, mut_str, mut_sing):
    # create copy of the original parameters
    par_mut = copy.deepcopy(par)

    # decide which parameters to mutate at each turn
    if mut_sing:
        mutate_set = [np.random.choice(keys_to_mut)]
    else:
        mutate_set = list(keys_to_mut)

    # implement mutations
    for k in mutate_set:
        if k in ('f_mem_reinit', 'g_1d', 'g_4d', 'a_selection', 'b_selection'):
            # mutate fractions: random walk, keep between 0 and 1
            par_mut[k] += (2. * np.random.rand() - 1.) * mut_str
            par_mut[k] = np.min([1., np.max([0., par_mut[k]])])
        elif k in ('sigma_i', 'alpha_C', 'k_consumption_per_day'):
            # mutate positive values: keep > 0
            par_mut[k] *= 1. + (2. * np.random.rand() - 1.) * mut_str
            par_mut[k] = np.max([0., par_mut[k]])
        elif k in ('mu_i', 'eps_B'):
            # mutate real values in the space of binding energy: random walk
            par_mut[k] += (2. * np.random.rand() - 1.) * 10. * mut_str
        else:
            raise Exception(
                f'error while trying to mutate parameter {k}: mutation rule not specified')

    # in case of contemporaneous mutation a + b one could incurr in the case
    # a + b > 1. Here we correct for this:
    delta = par_mut['a_selection'] + par_mut['b_selection'] - 1
    if delta > 0:
        par_mut['a_selection'] -= delta / 2.
        par_mut['b_selection'] -= delta / 2.

    return par_mut


def mc_accept(logls, new_logls, betas):
    rand = np.random.rand(logls.size)
    boltz_fact = betas * (new_logls - logls)
    is_accepted = rand < np.exp(boltz_fact)
    return is_accepted


def mc_switch(logls, betas):
    N = logls.size
    logls_cpy = np.copy(logls)
    order = np.arange(N)

    # attempt exchange layer N - n with layer N-n-1
    for n in range(N - 1):
        # high and low temeperature indices
        idx_hT = N - n - 1
        idx_lT = N - n - 2
        # evaluate delta log-likelihood and temperature
        delta_logl = logls_cpy[idx_hT] - logls_cpy[idx_lT]
        delta_beta = betas[idx_hT] - betas[idx_lT]
        boltz_fact = - delta_beta * delta_logl
        # randomly accept switch
        if np.random.rand() < np.exp(boltz_fact):
            # if accepted invert likelihood and update order
            logls_cpy[idx_hT], logls_cpy[idx_lT] = logls_cpy[idx_lT], logls_cpy[idx_hT]
            order[idx_hT], order[idx_lT] = order[idx_lT], order[idx_hT]

    # which layer's order has changed
    is_switched = order != np.arange(N)

    return is_switched, order
