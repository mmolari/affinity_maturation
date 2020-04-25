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
        tot_logl += dset_logl(dset, det_pf)
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


def initialize_layers(par_i, n_layers, T_max, dset_list, pars_to_mutate,
                      beta_list, mut_strength_list, mut_single_list):
    '''
    # TODO: specify order of parameters
    # TODO: implement checks on arguments
    '''

    # create a list of parameter sets, all equal to the initial set
    par_list = np.array([copy.deepcopy(par_i) for _ in range(n_layers)])
    # evaluate log-likelihood of the initial parameters set
    logl_0 = dset_list_logl(par_i, dset_list)
    # initialize array of log-likelihoods
    logl_list = np.ones(n_layers) * logl_0
    # initialize id of parameters
    par_ids = np.arange(n_layers)

    # initialize search parameters dictionary
    search_par = {}
    # save container for best parameter so far and best logl
    search_par['par_best'] = copy.deepcopy(par_i)
    search_par['best_logl'] = logl_0
    # save initial parameters and log-likelihood
    search_par['par_i'] = copy.deepcopy(par_i)
    search_par['logl_i'] = logl_0
    # save search time
    search_par['T_max'] = T_max
    # save parameter keys to mutate
    search_par['pars_to_mut'] = pars_to_mutate
    # number of layers
    search_par['n_layers'] = n_layers

    # initialize inverse temperature
    if beta_list is None:
        # (log-distributed between 10^3 and 10^-3)
        search_par['beta'] = np.logspace(-3, 3, n_layers)[::-1]
    else:
        search_par['beta'] = beta_list
    # initialize mutation strengths per layer
    if mut_strength_list is None:
        search_par['mut_strength'] = np.logspace(-2, -1, n_layers)
    else:
        search_par['mut_strength'] = mut_strength_list
    # initialize which layers feature a single mutation (half)
    if mut_single_list is None:
        search_par['mut_sing'] = np.zeros(n_layers, dtype=np.bool)
        search_par['mut_sing'][:(n_layers // 2) + 1] = True
    else:
        search_par['mut_sing'] = mut_single_list

    # returns the "search dictionary"
    return par_list, logl_list, par_ids, search_par


def save_initial_setup(dset_list, search_par, save_folder):
    '''
    '''
    # save initial parameters
    with open(os.path.join(save_folder, 'par_i.pkl'), 'wb') as f:
        pkl.dump(dset_list, f)
        f.close()

    # save dataset
    with open(os.path.join(save_folder, 'dataset_list.pkl'), 'wb') as f:
        pkl.dump(dset_list, f)
        f.close()

    # save search setup
    keys = ['betas', 'pars_to_mut', 'T_max', 'n_layers']
    betas, ptm, T_max, n_layers = [search_par[k] for k in keys]
    keys = ['mut_sing', 'mut_strength', 'logl_i', 'par_i']
    mut_single, mut_strength, logl_i, par_i = [search_par[k] for k in keys]

    with open(os.path.join(save_folder, 'search_setup.txt'), 'w') as f:
        f.write('search params:\n' + str(ptm) + '\n')
        f.write(f'n rounds = {T_max}\n')
        f.write(f'n layers = {n_layers}\n')
        f.write(f'inverse temperatures = {betas}\n')
        f.write(f'mutate one param at a time = {mut_single}\n')
        f.write(f'mut strength = {mut_strength}\n')
        f.write(f'avg logl init = {logl_i}\n\n')
        f.write('parameters initial value:\n----------------------\n')
        for key in par_i:
            f.write(f'{key:<30} - {par_i[key]}\n')
        f.close()


# def save_parameter(save_folder, par, logl, layer_n, beta, switch, round_n, traj_id):
#     # TODO: implement
#     par_cpy = copy.deepcopy(par)
#     par_cpy['logl'] = logl
#     par_cpy['layer'] = layer_n
#     par_cpy['temp'] = 1. / beta
#     par_cpy['switch'] = switch
#     par_cpy['round'] = round_n
#     par_cpy['par_traj_id'] = traj_id
#     path = os.path.join(
#         save_folder, f'par_layer_{layer_n}_round_{round_n:05}.pkl')
#     with open(path, 'wb') as f:
#         pkl.dump(par_cpy, f)
#         f.close()
#
#
# def save_layers(par_array, logl_arr, par_id, is_changed, is_switched,
#                 savefile_path, round_n, sw_order, betas, par_id):
#     # TODO: implement
#     N = par_array.size
#     is_modified = np.logical_or(is_changed, is_switched)
#     for n in range(N):
#         if is_modified[n]:
#             save_test_results(par_array[n], logl_arr[n],
#                               savefile_path, round_n, sw_order[n], n, betas[n], par_id[n])


def generate_variated_par(par, keys_to_mut, mut_str, mut_sing):
    # create copy of the original parameters
    par_mut = copy.deepcopy(par_0)

    # decide which parameters to mutate at each turn
    if mutate_single:
        mutate_set = [np.random.choice(keys_to_mutate)]
    else:
        mutate_set = list(keys_to_mutate)

    # implement mutations
    for k in mutate_set:
        if k in ('f_mem_reinit', 'g_1d', 'g_4d', 'a_selection', 'b_selection'):
            # mutate fractions: random walk, keep between 0 and 1
            par_mut[k] += (2. * np.random.rand() - 1.) * mut_strength
            par_mut[k] = np.min([1., np.max([0., par_mut[k]])])
        elif k in ('sigma_i', 'alpha_C', 'k_consumption_per_day'):
            # mutate positive values: keep > 0
            par_mut[k] *= 1. + (2. * np.random.rand() - 1.) * mut_strength
            par_mut[k] = np.max([0., par_mut[k]])
        elif k in ('mu_i', 'eps_B'):
            # mutate real values in the space of binding energy: random walk
            par_mut[k] += (2. * np.random.rand() - 1.) * 10. * mut_strength
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


def generate_variated_parlist(par_list, search_par):
    new_par_list = []
    for n_par, par in par_list:
        new_par = generate_variated_par(par,
                                        keys_to_mut=search_par['pars_to_mut'],
                                        mut_str=search_par['mut_strength'],
                                        mut_sing=search_par['mut_sing'])
        new_par_list.append(new_par)
    return np.array(new_par_list)


def accept_variated_parameters(logl_list, new_logl_list, betas):
    pass


def tempering_layer_switch(logl_list, betas):
    pass


# #  old library
#
#
#
# def MC_accept_pars(logl_array, mut_logl_array, beta_arr):
#     N = logl_array.size
#     rand = np.random.rand(N)
#     boltz_fact = beta_arr * (mut_logl_array - logl_array)
#     is_accepted = rand < np.exp(boltz_fact)
#     return is_accepted
#
#
# def MC_switch_pars(logl_arr, beta_arr):
#     N = logl_arr.size
#     logl_arr_cpy = np.copy(logl_arr)
#     order = np.arange(N)
#     is_switched = np.zeros(N, dtype=np.bool)
#
#     for n in range(N - 1):
#         # attempt exchange layer N - n with layer N-n-1
#         idx_hT = N - n - 1
#         idx_lT = N - n - 2
#         delta_logl = logl_arr_cpy[idx_hT] - logl_arr_cpy[idx_lT]
#         delta_beta = beta_arr[idx_hT] - beta_arr[idx_lT]
#         boltz_fact = - delta_beta * delta_logl
#         if np.random.rand() < np.exp(boltz_fact):
#             # accept the change
#             logl_arr_cpy[idx_hT], logl_arr_cpy[idx_lT] = logl_arr_cpy[idx_lT], logl_arr_cpy[idx_hT]
#             is_switched[idx_hT], is_switched[idx_lT] = True, True
#             order[idx_hT], order[idx_lT] = order[idx_lT], order[idx_hT]
#
#     return logl_arr_cpy, is_switched, order
#
#

#
#
