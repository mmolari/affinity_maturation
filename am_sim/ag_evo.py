import numpy as np

from .utils import next_ag_concentration, metadict_append


class antigen:
    '''
    Class implementing Ag concentration evolution.
    '''

    def __init__(self, C_inj, t_inj, par_dict):
        '''
        Initializer of the concentration evolver. It must be initialized with
        the amount of injected concentration, the time of injection (to keep
        track of external time) and the model parameters.
        '''
        self.C = 0.  #  available Ag concentration
        self.C_res = C_inj  # reservoir Ag concentration
        self.t = t_inj  # evolution time in days
        self.days_per_turn = par_dict['days_per_turn']
        self.k_dec_turn = par_dict['k_decay_per_day'] * self.days_per_turn
        self.k_cons_turn = par_dict['k_consumption_per_day'] * \
            self.days_per_turn
        self.k_rel_turn = par_dict['k_release_per_day'] * self.days_per_turn
        self.history = {
            'C': [0],  # concentration history
            't': [t_inj],  # time in days history
            'C_inj': C_inj,  # injected dosage
            't_inj': t_inj,  # injection time in days
        }

    def evolve(self, NB):
        '''
        Evolves the Ag concentration for a turn of evolution. The number of
        B-cells 'NB' controls the Ag consumption rate.
        '''
        # evaluate total removal rate
        km = self.k_dec_turn + NB * self.k_cons_turn
        # evolve ag concentration one turn
        self.C, self.C_res = next_ag_concentration(
            C_av=self.C, C_res=self.C_res, k_minus=km, k_plus=self.k_rel_turn)
        # evolve time
        self.t += self.days_per_turn
        # save results
        metadict_append(self.history, {'C': self.C, 't': self.t})
