import numpy as np

from .utils import pick_founders_en, generate_stoch_mutations


class stoch_pop:
    '''
    Stochastic population function class. It consists of a population of
    cells, each one characterized by its binding energy for the Ag.
    The class features methods for duplication, mutation and selection of the
    population.

    The class displays the following methods:
    - __init__: default class constructor
    - create_from_en: initializes a population from a list of energies
    - create_empty: initializes an empty population
    - merge_with: modifies the population by merging it with the one passed as
        argument
    - duplicate: implements cell duplication without mutation
    - keep_only: keeps only the specified cells in the mask and eliminates the
        rest
    - select_with_psurv: given a survival probability it stochastically removes
        cells according to that probability
    - differentiate: implements cell differentiation and return the
        differentiated MC/PC populations
    - carrying_cap: implements a finite carrying capacity
    - add_mutations: implements random cell mutations (silent + lethal + a.a.)
    - expand: implements the combination of duplications and mutations that
        occur during a single evolution round.
    - bareps: returns the current value of bar-epsilon for the population
    - N: returns the current population size
    - energies: return the population's binding energies.
    '''

    def __init__(self, par, mc_seed=None, trace_clonotypes=False):
        '''
        Initializes the population using the parameter set specified.

        Parameters:
        - par: the parameters dictionary
        - mc_seed (stoch_pop object, optional): MC population. If specified
            then the founders are extracted randomly either from the naive
            distribution or from the MC population.
        - trace_clonotypes (bool, optional): if set to true then the clonotype
            of the cells is traced during evolution
        '''
        Nf = par['N_founders']
        Ni = par['N_i']
        # energies of the founder clones
        if mc_seed is None:
            # either randomly extracted from the naive distribution
            self.en_founders = np.random.normal(
                loc=par['mu_i'], scale=par['sigma_i'], size=Nf)
        else:
            # or extracted randomly from both recalled MCs and naive clones.
            self.en_founders = pick_founders_en(par, mc_seed.en)
        # energies of the cells in the population, created repeating the
        # energies of the founders
        self.en = np.tile(self.en_founders, Ni // Nf + 1)[:Ni]
        self.tct = trace_clonotypes  # whether to trace the clonotype families
        if self.tct:
            # clonotype of cells in the population, created in the same way
            # as energies
            self.clonotype = np.tile(np.arange(Nf), Ni // Nf + 1)[:Ni]

    @classmethod
    def create_from_en(cls, en):
        '''
        Initialize a population object from simply the energy of the cells.

        Parameters:
        - en (array): list of energies of the newly-instantiated population.
        '''
        pop = cls.__new__(cls)
        pop.en = np.array(en)
        self.tct = False
        return pop

    @classmethod
    def create_empty(cls, *args):
        '''
        Initialize an empty population object. The extra argument is needed for compatibility with pop_deterministic but has not effect.
        '''
        pop = cls.__new__(cls)
        pop.en = np.array([])
        self.tct = False
        return pop

    def merge_with(self, pop_add):
        '''
        Merges the energy of two populations. NB: this operation is not
        compatible with clonotype tracing.

        Parameters:
        - pop_add (stoch_pop object): population of cells whose energy is to
            be merged with the current population.
        '''
        assert not self.tct, \
            'clonotype tracing is not well-defined if new cells are added'
        self.en = np.concatenate(self.en, pop_add.en)

    def duplicate(self):
        '''
        Function that encodes population duplication. Each cells divides into
        two daughter cells with the same affinity.
        '''
        self.en = np.repeat(self.en, 2)
        if self.tct:
            self.clonotype = np.repeat(self.clonotype, 2)

    def keep_only(self, mask):
        '''
        Utility function to keep only the cells specified in the mask and
        eliminate the rest.

        Parameters:
        - mask (array of bool): boolean array of cells to keep.
        '''
        self.en = self.en[mask]
        if self.tct:
            self.clonotype = self.clonotype[mask]

    def select_with_psurv(self, psurv):
        '''
        given a list 'psurv' containing the probability of survival for each
        cell in the population, it implements random selection with the
        specified probability.

        Parameters:
        - psurv (array of float): array containing the survival probability for
            each cell in the population. NB: probabilities must be in the same
            order as the cells.
        '''
        N = self.en.size
        surv_mask = np.random.rand(N) < psurv
        self.keep_only(surv_mask)

    def differentiate(self, prob_mc, prob_pc):
        '''
        This function takes as input the probability of MC and PC
        differentiation, stochastically performs differentiation by removing
        the differentiated cells from the population and returns the energies
        of the differentiated PCs and MCs.

        Parameters:
        - prob_mc, prob_pc (float): probability of MC / PC fate.
        '''
        # 0 for no differentiation, 1 for MC fate and 2 for PC fate
        outcome = np.random.choice(
            [0, 1, 2],
            size=self.en.size,
            p=[1 - prob_mc - prob_pc, prob_mc, prob_pc]
        )
        # create MC and PC populations
        MC_pop = stoch_pop.create_from_en(self.en[outcome == 1])
        PC_pop = stoch_pop.create_from_en(self.en[outcome == 2])
        # remove differentiated cells from population
        self.keep_only(outcome == 0)
        # return MCs and PCs pop
        return MC_pop, PC_pop

    def carrying_cap(self, par):
        '''
        This function implements a finite carrying capacity. If the population
        size exceed the 'GC_carry_capacity' threshold then the excess cells are
        discarded randomly.

        Parameters:
        - par: model parameters dictionary
        '''
        N_en = self.en.size  # population current size
        if N_en > par['GC_carry_capacity']:
            # randomly select the id of cells that will be kept
            surv_id = np.random.choice(
                np.arange(self.en.size), size=par['GC_carry_capacity'], replace=False)
            mask_in = np.zeros_like(self.en, dtype=np.bool)
            mask_in[surv_id] = True
            self.keep_only(mask_in)

    def add_mutations(self, par):
        '''
        This function adds mutations to the population. Mutations are applied
        according to the model definition and can be either silent, lethal or
        affinity-affecting.

        Parameters:
        - par: model parameters dictionary
        '''
        # also lethal mutations
        sil, aa, let = par['p_sil_eff'], par['p_aa_eff'], par['p_let_eff']
        outcome = np.random.choice(
            [0, 1, 2], size=self.en.size, p=[sil, aa, let])
        # add aa mutation energy (before removing!!!)
        mask_mut = outcome == 1
        delta_en = generate_stoch_mutations(self.par, mask_mut.sum())
        self.en[mask_mut] += delta_en
        # clean out lethal mutations
        mask_out = outcome == 2
        self.keep_only(~mask_out)

    def expand(self, par):
        '''
        This function implements the population expansion part of the evolution
        round, featuring multiple cell duplications with mutation.

        Parameters:
        - par: model parameters dictionary
        '''
        for _ in range(par['n_duplications']):
            self.duplicate()
            self.add_mutations(par)

    def bareps(self):
        '''
        This function return the current value of bar-epsilon for the
        population.
        '''
        beps = -np.log(np.exp(-self.en).mean())
        return beps

    def N_cells(self):
        '''
        This function returns the current population size.
        '''
        return self.en.size

    def energies(self):
        '''
        returns the population energies. NB: they are returned by reference.
        Therefore one must be careful not to modify them!
        '''
        return self.en
