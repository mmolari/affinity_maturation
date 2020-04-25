import numpy as np
import multiprocessing as mpc
from copy import deepcopy
import functools as fct

from .inference_utils import init_search_directory, initialize_layers
from .inference_utils import save_initial_setup, dset_list_logl
from .inference_utils import generate_variated_parlist


def parallel_tempering(dset_list, par_i, n_layers, T_max,
                       pars_to_mutate, save_folder, verbose=True,
                       beta_list=None, mut_strength_list=None,
                       mut_single_list=None):
    '''
    # TODO: add parameter to save only the final result?
    Describe requirements on parameters passed
    - order of temperatures and other parameters

    Describe all files produced
    - save search setup
    - save initial parameters
    - save list of datasets
    - save all parameters change
    - save best parameter

    Args:
    - dset_list
    - par_i
    - n_layers
    - beta_arr
    - T_max
    - pars_to_mutate
    - save_folder (string): folder in which to save all the search parameters.
        For precaution it must be an empty or non-existent folder. In the
        latter case it will be created.
    '''

    # if save folder does not exists create it
    print('Initializing directory')
    init_search_directory(save_folder)

    # initialize layers and search parameters.
    print('Initializing layer and evaluating initial logl')
    par_list, logl_list, par_ids, search_par = initialize_layers(
        par_i, n_layers, T_max, dset_list, pars_to_mutate,
        beta_list, mut_strength_list, mut_single_list
    )

    # save search parameters, dataset and initial parameter choice
    print('Saving initial parameters, dataset list and search setup')
    save_initial_setup(dset_list, search_par, save_folder)

    # spawn pool of workers for parallel evaluation
    n_procs = np.min([n_layers, mpc.cpu_count()])
    print(f'Generating a pool of {n_procs} workers')

    # define function to evaluate posterior log-likelihood of parameters
    par_logl_fun = fct.partial(dset_list_logl, dset_list=dset_list)

    with mpc.Pool(processes=n_procs) as pool:

        # start maximization cycle:
        for t in range(T_max):
            print(f'round {t + 1} / {T_max}')
            # produce random parameters
            new_par_list = generate_variated_parlist(par_list, search_par)

            # in parallel make simulation and evaluate logl (deterministic)
            new_logl_list = pool.map(par_logl_fun, new_par_list)
            new_logl_list = np.array(new_logl_list)

            # monte-carlo step to accept variated parameters
            is_accepted = accept_variated_parameters(
                logl_list, new_logl_list, betas=search_par['beta'])
            par_list[is_changed] = new_par_list[is_changed]
            logl_array[is_changed] = mut_logl_array[is_changed]

            # parallel tempering step to switch layers
            logl_list, is_switched, order = tempering_layer_switch(
                logl_list, betas=search_par['beta'])
            par_array = par_array[order]
            par_id = par_id[order]

            # save all parameter changes

            # check if found a better likelihood, update best parameters

        pool.close()
    # at the end of the search save the best parameter set found (+ layer and temp)


#
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
