import numpy as np
import copy

from .utils import prob_mc_pc_differentiation, Bsel_psurv, Tsel_psurv
from .ag_evo import antigen
from .pop_stochastic import stoch_pop
from .pop_deterministic import det_pop


class germinal_center:
    '''
    This class implements a single germinal centers. It displays two methods:
    the class initializer and the method to perform a single evolution round.
    The evolver method returns the population of MCs and PCs generated during
    the round. These could also be empty.
    '''

    def __init__(self, GC_type, par, C_inj, t_inj=0, mc_seed=None):
        '''
        The germinal center (GC) object is created by specifying the kind of
        evolution ('deterministic' or 'stochastic'), the model parameters
        dictionary and the concentration of Ag injected. Optionally the time of
        injection can be specified, and also the population of MCs from which
        seeder clones can be extracted.

        Parameters:
        - GC_type (string): either 'stochastic' or 'deterministic' according to
            the desired evolution
        - par: dictionary containing the model parameters
        - C_inj (float): injected Ag concentration
        - t_inj (float, optional): time of injection in days. If not specified
            then it is set to zero.
        - mc_seed (pop_class object, optional): population of MC from which
            the initial population can be partially seeded.
        '''
        # save parameters (makes a copy)
        self.par = copy.deepcopy(par)
        # initialize Ag concentration
        self.ag = antigen(C_inj=C_inj, t_inj=t_inj, par_dict=par)
        # current time in days and rounds
        self.t_days = t_inj
        self.t_rounds = 0
        # day of GC formation
        self.t_formation = par['T_GC_formation_days'] + t_inj
        # initialize population, either stochastic or deterministic evolution
        if GC_type == 'stochastic':
            self.pop_class = stoch_pop
        elif GC_type == 'deterministic':
            self.pop_class = det_pop
        else:
            raise Exception(
                'The GC type must be either stochastic or deterministic')
        self.pop = self.pop_class(par, mc_seed)
        # instantiate state (maturing - mature - extinct)
        self.state = 'maturing'

    def __evolve_maturing(self):
        '''
        Utility method to evolve a maturing GC. It simply evolves the Ag
        concentration considering an exponentially increasing population of
        B-cells.
        '''
        # consider an exponentially increasing number of cells
        exponent = self.t_rounds * \
            self.par['days_per_turn'] / self.par['T_GC_formation_days']
        NB = np.power(self.par['N_i'], exponent)
        # evolve concentration accordingly
        self.ag.evolve(NB=NB)

    def __evolve_mature(self):
        '''
        Utility function to perform an evolution round. It implements cell
        differentiation, duplication, mutation, selection, and GC carrying
        capacity. It returns the populations of differentiated MCs and PCs.
        '''
        # --- evolve population
        # differentiation
        p_mc, p_pc = prob_mc_pc_differentiation(self.par, self.t_rounds)
        MC, PC = self.pop.differentiate(p_mc, p_pc)

        # duplication and mutation
        self.pop.expand(self.par)

        # B selection
        en = self.pop.energies()
        psurv_B = Bsel_psurv(en=en, C=self.ag.C, par=self.par)
        self.pop.select_with_psurv(psurv_B)

        # T selection
        if self.pop.N_cells() > 0:
            en = self.pop.energies()
            bareps = self.pop.bareps()
            psurv_T = Tsel_psurv(en=en, bareps=bareps,
                                 C=self.ag.C, par=self.par)
            self.pop.select_with_psurv(psurv_T)

        # carrying capacity
        self.pop.carrying_cap(self.par)

        # --- evolve Ag concentration
        self.ag.evolve(NB=self.pop.N_cells())

        return MC, PC

    def evolve_one_round(self):
        '''
        Performs a single evolution round for the GC. If the state of the GC is
        'maturing' then only Ag concentration is updated, until the formation
        time is elapsed. At this point a single differentiation operation is
        performed and the state is set to 'mature'. From then on evolution
        consists in the usual round of duplication with mutation, selection,
        differentiation and enforcement of carrying capacity. The function
        returns the MCs and PCs produced in this round. These could also be
        empty.

        Returns:
        - MCs, PCs (pop_class objects): either stochastic or deterministic
            populations containing MCs/PCs generated during the round. These
            could be empty populations.
        '''
        # create empty MC and PC populations
        MC = self.pop_class.create_empty(self.par)
        PC = self.pop_class.create_empty(self.par)
        # if maturing then evolve Ag concentration
        if self.state == 'maturing':
            self.__evolve_maturing()
            # check if maturation has been reached
            if self.t_days >= self.t_formation:
                # if so then update state
                self.state = 'mature'
        # if mature perform standard evolution round
        elif self.state == 'mature':
            MC, PC = self.__evolve_mature()
            # if extinct change state
            if self.pop.N_cells() < 1:
                self.state = 'extinct'
                # erase remaining population (deterministic case)
                self.pop = self.pop_class.create_empty(self.par)
        # evolves time
        self.t_days += self.par['days_per_turn']
        self.t_rounds += 1
        # returns PCs and MCs
        return MC, PC
