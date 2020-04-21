import numpy as np

from .utils import gaussian_pdf, mutation_kernel


class det_pop:
    '''
    Deterministic population function class. This class implements the
    deterministic evolution of a population of cells, as defined in our model.
    A population is represented as a continuous distribution in the binding
    energy space. The class features methods to perform cell duplication and
    mutation, selection and differentiation.

    The class displays the following methods:
    - __init__: default class constructor
    - create_empty: initializes an empty population
    - create_with_explicit_attributes: creates an object having specified
        population size and distribution
    - create_copy_without_kernel: creates a copy of an object, copying every
        attribute but the mutation kernel.
    - merge_with: modifies the population by merging it with the one passed as
        argument
    - select_with_psurv: given a survival probability it models the effect of
        selection on the population according to this survival probability.
    - differentiate: implements cell differentiation and returns the
        differentiated MC/PC populations
    - carrying_cap: implements a finite carrying capacity
    - expand: implements the combination of duplications and mutations that
        occur during a single evolution round.
    - bareps: returns the current value of bar-epsilon for the population
    - N: returns the current population size
    - energies: returns the energy domain of the distribution (by reference!)
    - mean_en: returns the mean energy of the population
    '''

    def __init__(self, par, mc_seed=None):
        '''
        Initializes the population using the parameter set specified.

        Args:
        - par: the model parameters dictionary
        - mc_seed (stoch_pop object, optional): MC seed population. If
            specified then the population is seeded by a weighted mixture of
            reactivated memory and naive cells, according to the weight
            specified in the parameters.
        '''
        xlims, dx, Ni = par['xlims'], par['dx'], par['N_i']
        mu_i, sigma_i = par['mu_i'], par['sigma_i']
        # distribution domain and discretization step
        self.x = np.arange(xlims[0], xlims[1], dx)
        self.dx = dx
        # number of cells in the population
        self.N = Ni
        # naive cells normalized distribution
        self.varphi = gaussian_pdf(x=self.x, mu=mu_i, sigma=sigma_i)
        # if mc_seed specified then initialize distribution with a mixture of
        # naive and memory cells
        if mc_seed is not None:
            # the weight is specified in the parameters dictionary
            if par['f_mem_reinit'] == 'pop':
                # it either depends on the amount of MCs collected so far
                w = mc_seed.N / (self.N + mc_seed.N)
            else:
                # or it is a constant fraction
                w = par['f_mem_reinit']
            self.varphi = self.varphi * (1. - w) + mc_seed.varphi * w
        # build mutation kernel
        _, self.ker = mutation_kernel(par)

    @classmethod
    def create_empty(cls, par):
        '''
        Initialize an empty population. Both the distribution and the
        population size are set to zero.

        Args:
        - par: model parameters dictionary.
        '''
        pop = cls.__new__(cls)
        pop.N = 0  #  zero population size
        # create distribution domain according to model parameters.
        pop.x = np.arange(par['xlims'][0], par['xlims'][1], par['dx'])
        pop.dx = par['dx']
        pop.varphi = np.zeros_like(pop.x)  #  null distribution
        return pop

    @classmethod
    def create_with_explicit_attributes(cls, N, x, dx, varphi):
        '''
        Creates a new object having the attributes passed as argugment. Lists
        are copied in the process.

        Args:
        - N (float): population size
        - x (float array): distribution energy domain
        - dx (float): discretization interval of the energy domain
        - varphi (float array): values of the normalized distribution
        '''
        pop = cls.__new__(cls)
        # initialize parameters with the arguments specified
        pop.N = N
        pop.x = np.copy(x)  # creates a copy
        pop.dx = dx
        pop.varphi = np.copy(varphi)  # creates a copy
        return pop

    def create_copy_without_kernel(self):
        '''
        Creates a copy of the caller. It copies everything attribute except the
        mutation kernel, which is usually not needed in the copy.
        '''
        pop = det_pop.create_with_explicit_attributes(
            self.N, self.x, self.dx, self.varphi)
        return pop

    def merge_with(self, pop_add):
        '''
        Function that merges the current population with the population passed
        as argument.

        Args:
        - pop_add (det_pop object): population to be merged with the caller.
        '''
        # weight of the normalized distribution sum
        w = self.N / (self.N + pop_add.N)
        # merge distributions
        self.varphi = self.varphi * w + pop_add.varphi * (1. - w)
        # add up sizes
        self.N += pop_add.N

    def __renormalize_varphi(self):
        '''
        Renormalize the distribution after an operation that changes its size,
        and report the modification to the population size. This method should
        remain private.
        '''
        # evaluate the current normalization of the distribution
        N_factor = np.sum(self.varphi) * self.dx
        # update population size with the resulting factor
        self.N *= N_factor
        # renormalize the distribution
        self.varphi /= N_factor

    def select_with_psurv(self, psurv_x):
        '''
        Given a probability of survival, this method applies it to the
        population.

        Args:
        - psurv_x (float array): this array should contain the survival
            probability as a function of the energy domain of the distribution.
        '''
        # multiply the distribution by the probability of survival
        self.varphi *= psurv_x
        # renormalize the distribution and update population size
        self.__renormalize_varphi()

    def differentiate(self, prob_mc, prob_pc):
        '''
        This function implements differentiation. It returns the resulting
        populations of MCs and PCs.

        Args:
        - prob_mc, prob_pc (float): probabilities of respectivelt MC and PC
            differentiation.

        Returns:
        - MC_pop. PC_pop (det_pop objects): populations of differentiated
            MCs and PCs
        '''
        # create differentiated MC population from a copy of the current pop
        MC_pop = self.create_copy_without_kernel()
        # multiplied by the probability of differentiation
        MC_pop.N *= prob_mc
        # same for the plasma cell population
        PC_pop = self.create_copy_without_kernel()
        PC_pop.N *= prob_pc
        # remove the differentiated cells from the population size
        self.N *= (1 - prob_mc - prob_pc)

        return MC_pop, PC_pop

    def carrying_cap(self, par):
        '''
        This function implements a finite carry capacity.

        Args:
        - par: model parameters dictionary
        '''
        # if population size exceeds the carrying capacity remove the excess
        self.N = np.min([self.N, par['GC_carrying_capacity']])

    def expand(self, *args):
        '''
        This function it implements population expansion and mutation according
        to the model parameters.
        '''
        # perform convolution (amplification + mutation multiple times)
        self.varphi = np.convolve(self.ker, self.varphi,
                                  'same') * self.dx
        # renormalize the distribution and update population size
        self.__renormalize_varphi()

    def bareps(self):
        '''
        This function evaluate and returns the current value of bar-epsilon
        for the population.
        '''
        beps = -np.log(np.dot(self.varphi, np.exp(-self.x)) * self.dx)
        return beps

    def N_cells(self):
        '''
        This function returns the current population size.
        '''
        return self.N

    def energies(self):
        '''
        returns the distribution domain. NB: it is returned by reference.
        Therefore one must be careful not to modify them!
        '''
        return self.x

    def mean_en(self):
        '''
        returns the mean binding energy of the population. It returns None if
        the population is empty.
        '''
        norm = np.sum(self.varphi) * self.dx
        if norm > 0:
            return np.dot(self.x, self.varphi) * self.dx / norm
        else:
            return None
