'''
This file contains a dictionary with the standard value of the model parameters.
Unelss otherwise specified the parameters will take this value.
'''

import numpy as np
import copy

# ---- Evaluate mutation probability
L = 50                      # lenght of AA chain involved in binding
pmut_per_base_per_division = 1e-3
p_sil = 0.5                 # probability of silent mutation
p_aa = 0.2                  # probability of affinity-affecting mutation
p_let = 0.3                 # probability of lethal mutation
assert p_sil + p_aa + p_let == 1

# probability of at least one mutation in one nucleotide
p_mut = 1 - np.power(1 - pmut_per_base_per_division, 3 * L)
# effective probabilities of mutation per duplication event
p_sil_eff = p_sil * p_mut + (1 - p_mut)  # no mutation + silent mutation
p_aa_eff = p_aa * p_mut
p_let_eff = p_let * p_mut

# High energy threshold, low and high energy cutoffs
low_en_threshold = np.log(50e-9)
high_en_exp_cutoff = np.log(500e-9)
low_en_exp_cutoff = np.log(0.1e-9)

# standard parameters value
std_par_val = {
    # --- Ag dynamics
    'k_decay_per_day': 1.22248e-2,
    'k_consumption_per_day': 2.0681836905240632e-05,
    'k_release_per_day': 0.97856,
    # conversion factor Ag. C to dosage, in unit of micrograms of Ag
    'alpha_C': 0.024705238382797094,
    # --- GC time specifications
    'days_per_turn': 0.5,  #  days per turn
    'T_GC_formation_days': 6,
    'GC_carrying_capacity': 2500,  # carrying capacity
    # mc seeding fraction. If set to 'pop' then the fraction depends on the mc population size
    'f_mem_reinit': 'pop',
    # --- mutations
    'ker_xlim': 20,  # numerical mutation kernel energy half-extension
    'ker_ln_mu': 1.9,
    'ker_ln_sigma': 0.5,
    'ker_ln_offset': -3.,
    'p_sil_eff': p_sil_eff,
    'p_aa_eff': p_aa_eff,
    'p_let_eff': p_let_eff,
    'n_duplications': 2,  # number of dupl. (and therefore mut.) per round
    # --- B-selection
    'B_sel': False,
    'eps_B': -13.6,
    # --- T-selection
    'T_sel': True,
    'a_selection': 0.1332447195425596,  # selection permissivity
    'b_selection': 0.6609819950474214,  # selection additional rejection rate
    # --- differentiation
    # differentiation probability (MC + PC)
    'diff_prob': 0.1,
    # fraction of residual MC/PC differentiation after/before the switch
    'diff_residual_fraction': 0.,
    # MC/PC switch time in days. If None the switch is not performed and
    # differentiation occurs with 1/2 probability for each fate.
    'diff_switch_time': 11,
    # width of the sigmoid in days
    'diff_switch_sigma': 2,
    # --- initial naive cell distribution
    'mu_i': -14.625558761797766,
    'sigma_i': 1.6552173695143662,
    'N_i': 2500,
    'N_founders': 100,  #  n. of founder clones (only for stochastic sim.)
    # --- simulation energy discretization and limits
    'dx': 0.01,
    'xlim_minus': -50,
    'xlim_plus': 20,
    # --- measurement PC/MC mixture
    'g_1d': 0.5723468195195263,  #  MC/PC ratio, measurement 1 day after boost
    'g_4d': 0.0,  # MC/PC ratio, measurement 4 days after injection
}


def st_par():
    '''
    Returns a copy of the standard parameters dictionary. Otherwise it risks
    getting modified while running the code.
    '''
    return copy.deepcopy(std_par_val)
