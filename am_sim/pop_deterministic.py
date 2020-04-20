import numpy as np

from .utils import gaussian_pdf, mutation_kernel


class det_pop:
    '''
    '''

    def __init__(self, par, mc_seed=None):
        '''
        '''
        xlims, dx, Ni = par['xlims'], par['dx'], par['N_i']
        mu_i, sigma_i = par['mu_i'], par['sigma_i']
        # distribution domain and discretization step
        self.x = np.arange(xlims[0], xlims[1], dx)
        self.dx = dx
        # distribution size
        self.N = Ni
        # naive distribution
        self.varphi = gaussian_pdf(x=self.x, mu=mu_i, sigma=sigma_i)
        # if mc_seed specified then add MCs to the initial distribution
        if mc_seed is not None:
            if par['f_mem_reinit'] == 'pop':
                w = mc_seed.N / (self.N + mc_seed.N)
            else:
                w = par['f_mem_reinit']
            self.varphi = self.varphi * (1. - w) + mc_seed.varphi * w
        # build mutation kernel
        _, self.ker = mutation_kernel(par)

    @classmethod
    def create_empty(cls, par):
        '''
        Initialize an empty population
        '''
        pop = cls.__new__(cls)
        pop.N = 0  #  zero population size
        pop.x = np.arange(par['xlims'][0], par['xlims'][1], par['dx'])
        pop.dx = par['dx']
        pop.varphi = np.zeros_like(pop.x)  #  null distribution
        return pop

    def merge_with(self, pop_add):
        '''
        Function that merges the current population with the population passed
        as argument.
        '''
        # weight of the normalized distribution sum
        w = self.N / (self.N + pop_add.N)
        # merge distributions
        self.varphi = self.varphi * w + pop_add.varphi * (1. - w)
        # add up sizes
        self.N += pop_add.N

    def __renormalize_varphi(self):
        '''
        Renormalize varphi after an operation that changes it and modify N
        accordingly.
        '''
        N_factor = np.sum(self.varphi) * self.dx
        self.N *= N_factor
        self.varphi /= N_factor

    def select_with_psurv(self, psurv_x):
        '''
        Given a probability of survival applies it to the population
        '''
        self.varphi *= psurv_x
        self.__renormalize_varphi()

    def differentiate(self, prob_mc, prob_pc):
        '''
        '''
        MC_pop = copy.deepcopy(self)
        MC_pop.N *= prob_mc
        PC_pop = copy.deepcopy(self)
        PC_pop.N *= prob_pc
        self.N *= (1 - prob_mc - prob_pc)

    def carrying_cap(self, par):
        '''
        '''
        self.N = np.min([self.N, par['GC_carr_cap']])

    def expand(self, *args):
        '''
        '''
        # perform convolution (amplification + mutation multiple times)
        self.varphi = np.convolve(self.ker, self.varphi,
                                  'same') * self.dx
        self.__renormalize_varphi()

    def bareps(self):
        '''
        This function returns the current value of bar-epsilon for the
        population.
        '''
        beps = -np.log(np.dot(varphi, np.exp(-self.x)) * self.dx)
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
