import numpy as np

import am_sim as ams

from utilities.dataset_interface import load_all_datasets

# setup numpy random seed for reproducibility
np.random.seed(1)

# initial value of model parameters
par_i = ams.st_par()
par_i['k_consumption_per_day'] = 2.e-05
par_i['mu_i'] = -14.6
par_i['sigma_i'] = 1.6
par_i['g_1d'] = 0.5
par_i['g_4d'] = 0.5
par_i['a_selection'] = 0.2
par_i['b_selection'] = 0.2
par_i['alpha_C'] = 0.025

# load all datasets
dsets = load_all_datasets()

# setup search parameters
save_folder = 'inference_search_results'
pars_to_mutate = ['k_consumption_per_day', 'mu_i', 'sigma_i',
                  'g_1d', 'g_4d', 'a_selection', 'b_selection', 'alpha_C']
T_max = 10000
n_layers = 10

# launch search
inference = ams.parallel_tempering(dset_list=dsets,
                                   par_i=par_i,
                                   n_layers=n_layers,
                                   T_max=T_max,
                                   pars_to_mutate=pars_to_mutate,
                                   save_folder=save_folder,
                                   save_every=100)
inference.search()
