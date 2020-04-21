import numpy as np

from .utils import evaluate_responders
from .pop_stochastic import stoch_pop
from .pop_deterministic import det_pop
from .germinal_center import germinal_center

'''
This function implements immunization schemes. Be general or specific? Maybe
specific...

An immunization scheme is defined by 3 things:
- C_inj (list)
- T_delay (list)
- measure (how the measurement is performed, 1d vs 4d)

We are interested in knowing:
- sometimes the evolution
- during the inference only the outcome
'''


def simulate_immscheme(self, kind, D_inj, T_delay, measure, par):

    # conversion Dose -> Concentration
    C_inj = np.array(D_inj) / par['alpha_C']
    T_cum = np.cumsum(T_delay)

    # define g (responder memory fraction) based on the measurement
    # protocol
    if measure == '1d':
        g_mem = par['g_1d']
    elif measure == '4d':
        g_mem = par['g_4d']
    else:
        raise Exception('measure must be either 1d or 4d')

    # either stochastic or deterministic simulation
    if kind == 'stochastic':
        pop_class = stoch_pop
    elif kind == 'deterministi':
        pop_class = det_pop
    else:
        raise Exception(
            'simulation kind must be either stochastic or deterministic')

    # instantiate empty MC / PC populations
    MC, PC = pop_class.create_empty(par), pop_class.create_empty(par)

    # simulate immunization scheme
    # first injection at t = 0
    # cicle through subsequent injections
    # add memory cells and plasma cells

    # merge MCs and PCs into responder populations
    resp_pop = evaluate_responders(MC, PC, g_mem, kind)
    return resp_pop
