import numpy as np
import copy
import multiprocessing as mpc
import pandas as pd

from .inference_utils import init_search_directory, dset_list_logl
from .inference_utils import save_search_initial_setup
from .inference_utils import generate_variated_par
from .inference_utils import mc_accept, mc_switch


class parallel_tempering:

    def __init__(self, dset_list, par_i, n_layers, T_max, pars_to_mutate,
                 save_folder, beta_list=None, mut_strength_list=None,
                 mut_single_list=None, save_every=100):

        # if save folder does not exists create it
        print('Initializing directory')
        init_search_directory(save_folder)

        # initialize search parameters
        self.dsets = dset_list  # list of datasets
        self.n_layers = n_layers  # number of parallel tempering layers
        self.T_max = T_max  #  maximum number of search steps
        self.pars_to_mutate = pars_to_mutate  # list of parameters to vary
        self.save_folder = save_folder  #  save directory
        self.save_every = save_every  #  number of iteration between two saves

        # list of layer temperatures.
        if beta_list is None:
            self.betas = np.logspace(-3, 3, n_layers)[::-1]
        else:
            self.betas = np.array(beta_list)

        # mutation strength
        if mut_strength_list is None:
            self.mut_str = np.logspace(-2, -1, n_layers)
        else:
            self.mut_str = np.array(mut_strength_list)

        # in which layer mutate a single parameter
        if mut_single_list is None:
            self.mut_sing = np.zeros(n_layers, dtype=np.bool)
            self.mut_sing[:(n_layers // 2) + 1] = True
        else:
            self.mut_sing = np.array(mut_single_list)

        # initialize layers and save initial setup
        self.init_layers_and_save(par_i)

    def search(self):

        # initialize empty search history archive
        self.df = pd.DataFrame()

        # save initial state for all layers
        print('Saving initial state of all layers')
        self.history_append_state(t=0, is_accepted=np.ones(self.n_layers),
                                  is_switched=np.zeros(self.n_layers))

        # define function to evaluate posterior log-likelihood of parameters
        logl_funct = fct.partial(dset_list_logl, dset_list=self.dsets)

        # spawn pool of workers for parallel evaluation
        n_procs = np.min([n_layers, mpc.cpu_count()])
        print(f'Generating a pool of {n_procs} workers')

        with mpc.Pool(processes=n_procs) as pool:

            # start maximization cycle:
            for t in range(1, T_max + 1):
                print(f'round {t} / {T_max}')

                # every one hundred iterations save search history
                if t % self.save_every == 0:
                    print(f'Save search history at search round t = {t}')
                    self.save_search_history(t=t)

                # produce random parameters
                new_pars = self.vary_pars()

                # in parallel make simulation and evaluate logl (deterministic)
                new_logls = pool.map(logl_funct, new_pars)
                new_logls = np.array(new_logls)

                # monte-carlo step to accept variated parameters
                is_accepted = mc_accept(self.logls, new_logls, self.betas)
                self.pars[is_accepted] = new_pars[is_accepted]
                self.logls[is_accepted] = new_logls[is_accepted]

                # parallel tempering step to switch layers
                is_switched, order = mc_switch(self.logls, betas=self.betas)
                # update the new order
                self.logls = self.logls[order]
                self.pars = self.pars[order]
                self.traj_id = self.traj_id[order]

                # save all parameter changes
                self.history_append_state(t, is_accepted, is_switched)

            pool.close()

        # save final version of the search history
        self.save_search_history(t='final')

    def init_layers_and_save(self, par_i):

        # create a list of parameter sets, all equal to the initial one
        self.pars = [copy.deepcopy(par_i) for _ in range(self.n_layers)]
        self.pars = np.array(self.pars)

        # evaluate log-likelihood of the initial parameters set
        logl_0 = dset_list_logl(par_i, self.dsets)

        # initialize array of log-likelihoods
        self.logls = np.ones(n_layers) * logl_0

        # initialize id of parameter sets
        self.traj_id = np.arange(n_layers)

        # save search parameters, dataset and initial parameter choice
        save_search_initial_setup(self)

    def history_append_state(self, t, is_accepted, is_switched):
        '''
        # TODO: save on the same pandas database?
        '''
        # which parameters have changed since last round
        is_changed = np.logical_or(is_accepted, is_switched)

        # create a list of parameters that have been updated since last round
        updated_par_list = []
        # indices of parameters that have changed
        idx_ch = np.argwhere(is_changed).flatten()
        for idx in idx_ch:
            # create a copy of the parameter and add additional entries
            # related to the search state
            par = copy.deepcopy(self.pars[idx])
            par['logl'] = self.logls[idx]
            par['layer'] = idx
            par['temp'] = 1. / self.betas[idx]
            par['switch'] = is_switched[idx]
            par['round'] = t
            par['traj_id'] = self.traj_id[idx]
            updated_par_list.append(par)

        # append parameters to dataframe
        self.df.append(updated_par_list, ignore_index=True)

    def save_search_history(self, t):
        # find the old save if present in the folder
        files = os.listdir(self.save_folder)
        files = [f for f in files if f.endswith('search_history.csv')]
        old_filename = None
        if len(files) > 0:
            if len(files) == 1:
                old_filename = files[0]
            else:
                print('WARNING: multiple search_history files in folder.\n' +
                      'as a precaution not erasing previous history when' +
                      ' saving the new.')
        # save current history
        current_filename = f't_{t}_search_history.csv'
        self.df.to_csv(os.path.join(self.save_folder, current_filename))

        # remove previous save file if present
        if old_filename is not None:
            print(f'Removing old save file: {old_filename}')
            os.remove(os.path.join(self.save_folder, old_filename))

    def vary_pars(self):
        new_pars = []
        for n_par, par in self.pars:
            new_par = generate_variated_par(par,
                                            keys_to_mut=self.pars_to_mutate,
                                            mut_str=self.mut_str,
                                            mut_sing=self.mut_sing)
            new_pars.append(new_par)
        return np.array(new_pars)
