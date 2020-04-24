import numpy as np
import multiprocessing as mpc
from copy import deepcopy

from .utils import dset_logl


def parallel_tempering(dset, par_i, n_layers, beta_arr, T_max, pars_to_mutate,
                       save_folder=None):
    # if save_folder create the folder
    # initialize layers
    # save all parameters of the search, plus the initial parameter choice
    # start maximization cycle:
        # produce random parameters
        # in parallel make simulation and evaluate logl (deterministic)
        # randomly decide wether to accept and whether to switch
    # at the end of the search save the best parameter set found (+ layer and temp)
    pass

#
#
#
# def parallel_tempering(dataset, par_0, n_layers, beta_arr, mut_strength_arr, T_max, mut_single_arr,
#                        pars_to_mutate, savefile_path):
#
#     # initialize directory
#     if os.path.exists(savefile_path):
#         if not len(os.listdir(savefile_path)) == 0:
#             raise Exception(f'For safety reasons the user should provide an empty \
#             or non-existent folder. Insted the specified folder exists and \
#             contains the following files:\n {os.listdir(savefile_path)}')
#     else:
#         os.makedirs(savefile_path, exist_ok=False)
#
#     print('Evaluating initial logl...')
#     # initialize layers
#     par_array = np.array([copy.deepcopy(par_0) for _ in range(n_layers)])
#     par_0_logl = parameter_set_loglikelihood(dataset=dataset, par=par_0)
#     logl_array = np.ones(n_layers) * par_0_logl
#     par_id = np.arange(n_layers)
#
#     # layer zero has lower temperature (highest beta)
#     beta_order = np.argsort(beta_arr)[::-1]
#     beta_arr = beta_arr[beta_order]
#     mut_strength_arr = mut_strength_arr[beta_order]
#     mut_single_arr = mut_single_arr[beta_order]
#
#     # save search parameters
#     with open(os.path.join(savefile_path, 'search_params.txt'), 'w') as f:
#         f.write('search params:\n' + str(pars_to_mutate) + '\n')
#         f.write(f'n rounds = {T_max}\n')
#         f.write(f'n layers = {n_layers}\n')
#         f.write(f'temperatures = {beta_arr}\n')
#         f.write(f'mutate one param at a time = {mut_single_arr}\n')
#         f.write(f'mut strength = {mut_strength_arr}\n')
#         f.write(f'avg logl init = {par_0_logl}\n\n')
#         f.write('parameters initial value:\n----------------------\n')
#         for key in par_0:
#             f.write(f'{key:<30} - {par_0[key]}\n')
#         f.close()
#
#     # save the dataset, if needed for future checks
#     with open(os.path.join(savefile_path, 'dataset.pkl'), 'wb') as f:
#         pkl.dump(dataset, f)
#         f.close()
#
#     # save pars
#     print('Saving initial pars...')
#     save_changed_res(par_array, logl_array, is_changed=np.ones(n_layers, dtype=np.bool),
#                      is_switched=np.zeros(n_layers, dtype=np.bool),
#                      savefile_path=savefile_path, round_n=0, sw_order=np.arange(n_layers),
#                      betas=beta_arr, par_id=par_id)
#
#     # spawn pool and start the cycle
#     with mpc.Pool(processes=n_layers) as pool:
#
#         for t in range(T_max):
#             print(f'round {t}')
#
#             # generate mutated params
#             mut_par_array = generate_mutpar_array(
#                 par_array, pars_to_mutate, mut_strength_arr, mut_single_arr)
#
#             # parallel evaluate logl
#             mut_logl_array = pool.map(fct.partial(parameter_set_loglikelihood,
#                                                   dataset=dataset), mut_par_array)
#             mut_logl_array = np.array(mut_logl_array)
#
#             # decide and apply acceptance
#             is_changed = MC_accept_pars(logl_array, mut_logl_array, beta_arr)
#             # print(mut_logl_array - logl_array)
#             # print(mut_logl_array > logl_array)
#             # print(is_changed)
#             par_array[is_changed] = mut_par_array[is_changed]
#             # for nl in range(n_layers):
#             #     print(np.all([par_array[nl][k] == mut_par_array[nl][k]
#             #                   for k in par_array[nl].keys()]))
#             logl_array[is_changed] = mut_logl_array[is_changed]
#
#             # decide tempering layer switch
#             logl_array, is_switched, order = MC_switch_pars(
#                 logl_array, beta_arr)
#             par_array = par_array[order]
#             par_id = par_id[order]
#
#             # save params
#             save_changed_res(par_array, logl_array, is_changed,
#                              is_switched, savefile_path, t + 1, order,
#                              beta_arr, par_id)
#
#         pool.close()
