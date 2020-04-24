import numpy as np


def dset_logl(dset, det_pf):
    '''
    Gives the log-likelihood of the data for a given population function. The
    experimental sensitivity limits are taken into account.

    Args:
    - dset (dataset object): the experimental dataset
    - det_pf (det_pop object): deterministic population function, corresponding
        to the model prediction for this immunization scheme.

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


# #  old library
#
#
# def generate_mutated_parameters(par_0, keys_to_mutate, mutate_single,
#                                 mut_strength):
#     '''
#     generates a mutated set of parameters to test
#     TODO: fix the a,b part. Maybe leave more freedom?
#     '''
#     # create copy of the original parameters
#     par_mut = copy.deepcopy(par_0)
#
#     # decide which parameters to mutate at each turn
#     if mutate_single:
#         mutate_set = [np.random.choice(keys_to_mutate)]
#     else:
#         mutate_set = list(keys_to_mutate)
#
#     # implement mutations
#     for k in mutate_set:
#         if k in ('f_mem_reinit', 'g_1d', 'g_4d', 'a_selection', 'b_selection'):
#             # mutate fractions: random walk, keep between 0 and 1
#             par_mut[k] += (2. * np.random.rand() - 1.) * mut_strength
#             par_mut[k] = np.min([1., np.max([0., par_mut[k]])])
#         elif k in ('sigma_i', 'alpha_C', 'k_consumption_per_day'):
#             # mutate positive values: keep > 0
#             par_mut[k] *= 1. + (2. * np.random.rand() - 1.) * mut_strength
#             par_mut[k] = np.max([0., par_mut[k]])
#         elif k in ('mu_i', 'eps_B'):
#             # mutate real values in the space of binding energy: random walk
#             par_mut[k] += (2. * np.random.rand() - 1.) * 10. * mut_strength
#         # elif k == 'a_selection':
#         #     # fractions that must sum to 1: keep bound between the 0 and 1 - b:
#         #     par_mut[k] += (2. * np.random.rand() - 1.) * mut_strength
#         #     par_mut[k] = np.min(
#         #         [1. - original_par['b_selection'], np.max([0., par_mut[k]])])
#         # elif k == 'b_selection':
#         #     # keep bound between the 0 and 1 - a:
#         #     par_mut[k] += (2. * np.random.rand() - 1.) * mut_strength
#         #     par_mut[k] = np.min(
#         #         [1. - original_par['a_selection'], np.max([0., par_mut[k]])])
#         else:
#             raise Exception(
#                 f'error while trying to mutate parameter {k}: mutation rule not specified')
#
#     # in case of contemporaneous mutation a + b one could incurr in the case
#     # a + b > 1. Here we correct for this:
#     delta = par_mut['a_selection'] + par_mut['b_selection'] - 1
#     if delta > 0:
#         par_mut['a_selection'] -= delta / 2.
#         par_mut['b_selection'] -= delta / 2.
#
#     return par_mut
#
#
# def generate_mutated_par_array(par_array, pars_to_mutate, mut_strength_arr, mut_single_arr):
#     # TODO: comment
#     mut_par_array = []
#     for n in range(par_array.size):
#         mut_par_array.append(generate_mutated_parameters(par_array[n], pars_to_mutate, mut_single_arr[n],
#                                                          mut_strength_arr[n]))
#     return np.array(mut_par_array)
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
# def save_test_results(par, logl, savedir, round_n, switch, layer_n, beta, par_id):
#     '''
#     Save the single round parameter results along with their likelihood.
#     '''
#     par_cpy = copy.deepcopy(par)
#     par_cpy['logl'] = logl
#     par_cpy['layer'] = layer_n
#     par_cpy['temp'] = 1. / beta
#     par_cpy['switch'] = switch
#     par_cpy['round'] = round_n
#     par_cpy['par_traj_id'] = par_id
#     path = os.path.join(savedir, f'par_layer_{layer_n}_round_{round_n:05}.pkl')
#     with open(path, 'wb') as f:
#         pkl.dump(par_cpy, f)
#         f.close()
#
#
# def save_changed_res(par_array, logl_arr, is_changed, is_switched, savefile_path, round_n, sw_order, betas, par_id):
#     N = par_array.size
#     is_modified = np.logical_or(is_changed, is_switched)
#     for n in range(N):
#         if is_modified[n]:
#             save_test_results(par_array[n], logl_arr[n],
#                               savefile_path, round_n, sw_order[n], n, betas[n], par_id[n])
