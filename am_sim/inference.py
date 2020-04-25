import numpy as np
import copy
import multiprocessing as mpc
import pandas as pd
import functools as fct
import os

from .inference_utils import init_search_directory, dset_list_logl
from .inference_utils import save_search_initial_setup
from .inference_utils import generate_variated_par
from .inference_utils import mc_accept, mc_switch


class parallel_tempering:
    '''
    Class that implements the likelihood maximization inference procedure.
    It displays two main methods: the class initializer and the search method.

    The initializer methods creates an empty folder in the specified location,
    in which it saves the search parameters (search_setup.txt), the
    experimental datasets (dataset_list.pkl) and the initial value of the
    model parameters (par_i.pkl).

    The search method is parallelized and implements the parallel termpering
    technique to find the maximum-likelihood value of the model parameters.
    During the search results are progressively saved in a csv file, whose name
    ends with 'search_history.csv'. This file stores all the updates of all
    the parameters for every layer. Next to each parameter set also other
    values are saved, namely:
    - logl: the log-likelihood of the parameter set
    - layer: the layer to whom the parameter set belongs
    - temp: the temperature associated to the layer
    - switch: wether the parameter set was switched with another layer in this
        round
    - round: the round at which the parameter set was saved.
    - traj_id: the trajectory to whom the parameter set belongs. Trajectory ids
        are switched together with the parameter sets.
    '''

    def __init__(self, dset_list, par_i, n_layers, T_max, pars_to_mutate,
                 save_folder, beta_list=None, mut_strength_list=None,
                 mut_single_list=None, save_every=100):
        '''
        Class initializer for parallel_tempering class. This function
        initializes the search parameters, creates the folder in which to save
        data, and saves the initial state of the parameters, the experimental
        dataset and the initial value of the parameters.

        Args:
        - dset_list (list of dataset objects): list of dataset objects
            containing all the experimental measurements and schemes of which
            the likelihood must be maximized.
        - par_i (model parameters dictionary): initial value of the model
            parameters in the maximization algorithm.
        - n_layers (int): number of parallel tempering layers.
        - T_max (int): total number of iterations of the likelihood
            maximization procedure.
        - pars_to_mutate (list of str): list containing the parameters whose
            maximum-likelihood value must be found.
        - save_folder (str): name of the folder in which details and results of
            the search will be saved. As a precaution the folder must be either
            non-existent (in which case it will be created) or empty.
        - beta_list (optional, list of float): list of inverse temperatures per
            layer. If specified it must have the same dimensions of the number
            of layers. If not specified is initialized as log-spaced between
            10^3 and 10^-3.
        - mut_strength_list (optional, list of float): values of the parameters
            mutation strength, quantifying the variation magnitude of the
            parameters for each layer. It should be a small number (0.1~0.01).
            If specified it must have the same dimension of the number of
            layers, otherwise it is initialized as log-spaced between 0.1 and
            0.01.
        - mut_single_list (otpional, list of bool): wether parameters variation
            concerns all search parameters, or one at a time chosen randomly.
            If specified it must have the same dimension as the
        - save_every (optional, int): number of rounds between two successive
            update of the save-file. Default value is 100.
        '''

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

        # mutation strength for each layer
        if mut_strength_list is None:
            self.mut_str = np.logspace(-2, -1, n_layers)
        else:
            self.mut_str = np.array(mut_strength_list)

        # in which layer parameter variation concerns a single parameter at a time
        if mut_single_list is None:
            self.mut_sing = np.zeros(n_layers, dtype=np.bool)
            self.mut_sing[:(n_layers // 2) + 1] = True
        else:
            self.mut_sing = np.array(mut_single_list)

        # initialize layers and save initial setup
        print('Initializing layers')
        self.init_layers_and_save(par_i)

    def search(self):
        '''
        Search function. This function launches the likelihood-maximization
        algorithm. The algorithm is parallel and runs on a number of cores
        equal (or minor if not enough cores are avalilable) to the number of
        parallel tempering layers. Every 'save_every' number of rounds all the
        parameters variation for all the layers are saved in a csv file
        whose name ends in 'search_history.csv'.
        '''

        # initialize empty search history archive
        self.hist_df = pd.DataFrame()
        self.temp_hist = []

        # save initial state for all layers
        print('Saving initial state of all layers')
        self.history_append_state(t=0, is_accepted=np.ones(self.n_layers),
                                  is_switched=np.zeros(self.n_layers))

        # define function to evaluate posterior log-likelihood of parameters
        logl_funct = fct.partial(dset_list_logl, dset_list=self.dsets)

        # spawn pool of workers for parallel evaluation
        n_procs = np.min([self.n_layers, mpc.cpu_count()])
        print(f'Generating a pool of {n_procs} workers')

        with mpc.Pool(processes=n_procs) as pool:

            # start maximization cycle:
            for t in range(1, self.T_max + 1):
                print(f'round {t} / {self.T_max}')

                # produce random parameters
                new_pars = self.vary_pars()

                # in parallel evaluate logl of parameter sets for all layers
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

                # every one hundred iterations save search history
                if t % self.save_every == 0:
                    print(f'Save search history at search round t = {t}')
                    self.save_search_history(t=t)

            pool.close()

        # save final version of the search history
        self.save_search_history(t='final')

    def init_layers_and_save(self, par_i):
        '''
        Utility function used to initialize the layers and save the initial
        state of the search. It takes as argument the initial value of the
        model parameters.
        '''

        # create a list of parameter sets, all equal to the initial one
        self.pars = [copy.deepcopy(par_i) for _ in range(self.n_layers)]
        self.pars = np.array(self.pars)

        # evaluate log-likelihood of the initial parameters set
        logl_0 = dset_list_logl(par_i, self.dsets)

        # initialize array of log-likelihoods
        self.logls = np.ones(self.n_layers) * logl_0

        # initialize id of parameter sets
        self.traj_id = np.arange(self.n_layers)

        # save search parameters, dataset and initial parameter choice
        save_search_initial_setup(self)

    def history_append_state(self, t, is_accepted, is_switched):
        '''
        Utility function to save the current round of the search. It takes as
        argument the round number, and two boolean arrays. These array specify
        wether the parameter set corresponding to every layer have been changed
        or switched during the round. In any of these two cases the parameter
        set is saved, together with some extra information related to the
        search, in the temporary history list. This list will be regularly
        emptied when parameters are saved into the '.csv' file
        '''
        # which parameters have changed since last round
        is_changed = np.logical_or(is_accepted, is_switched)

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

            # add the changes to the temporary history
            self.temp_hist.append(par)

    def save_search_history(self, t):
        '''
        Utility function to empty the temporary history list and save the
        search results into the '.csv' file. It takes as argument the iteration
        round t, which is used to add a signature to the save-file. If a prior
        version of the save-file is present then it is removed.
        '''
        # append temporary history to full history database
        self.hist_df = self.hist_df.append(self.temp_hist, ignore_index=True)
        # empty temporary history
        self.temp_hist = []

        # find the old save file if present in the folder
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
        self.hist_df.to_csv(os.path.join(self.save_folder, current_filename))

        # remove previous save file if present
        if old_filename is not None:
            print(f'Removing old save file: {old_filename}')
            os.remove(os.path.join(self.save_folder, old_filename))

    def vary_pars(self):
        '''
        Utility function that from the values of the parameter set in every
        layer generates and returns a mutatated version of the parameters.
        Mutations is introduced in the form of a small parameter variation,
        whose intensity depends on the mutation strength of the layers, and
        wether a single or multiple mutation is allowed.
        '''
        # returns a list of parameters generated from the previous ones with
        # a small variation
        new_pars = []
        for n, par in enumerate(self.pars):
            new_par = generate_variated_par(par,
                                            keys_to_mut=self.pars_to_mutate,
                                            mut_str=self.mut_str[n],
                                            mut_sing=self.mut_sing[n])
            new_pars.append(new_par)
        return np.array(new_pars)
