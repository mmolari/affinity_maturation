from .inference_utils import init_search_directory


class parallel_tempering:

    def __init__(self, dset_list, par_i, n_layers, T_max, pars_to_mutate,
                 save_folder, beta_list=None, mut_strength_list=None,
                 mut_single_list=None):

        # if save folder does not exists create it
        print('Initializing directory')
        init_search_directory(save_folder)

        # initialize search parameters
        self.dsets = dset_list
        self.n_layers = n_layers
        self.T_max = T_max
        self.pars_to_mutate = pars_to_mutate
        self.save_folder = save_folder

        # initialize layers

        # save initial search setup

    def search(self):

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
