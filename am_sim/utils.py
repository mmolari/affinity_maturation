'''
This file containst utilities to be used in the other libraries.
'''
import numpy as np
import scipy.stats as sps
import copy

from .model_parameters import high_en_exp_cutoff, low_en_exp_cutoff, low_en_threshold


# --- meta-dictionary

def metadict_append(meta_dict, el):
    '''
    Appends the elements of a results dictionary to corresponding lists in the
    meta-dictionary. If the meta_dict is empty then it also creates the lists.
    '''
    if len(meta_dict) == 0:
        for k in el.keys():
            meta_dict[k] = [el[k]]
    else:
        for k in el.keys():
            meta_dict[k].append(el[k])


# --- distributions


def gaussian_pdf(x, mu, sigma):
    '''
    return a gaussian distribution pdf given mean and standard deviation
    '''
    dx = x[1] - x[0]
    rho_i = sps.norm.pdf(x=x, loc=mu, scale=sigma)
    return rho_i / (np.sum(rho_i) * dx)


def lognorm_pdf(x, dx, mu, sigma, offset):
    '''
    return a lognormal distribution pdf given mu, sigma and offset
    '''
    lognorm = sps.lognorm.pdf(x, s=sigma, loc=offset, scale=np.exp(mu))
    lognorm /= np.sum(lognorm) * dx
    return lognorm


# --- selection functions

def sigmoid_psurv(en, en_thr, a, b, C):
    '''
    Utility function implementing the sigmoid survival probability function.

    Parameters:
    - en (array): energies for which the survival probability must be evaluated
    - en_thr (float): threshold selection energy
    - a,b (float): stochasticity selection parameters
    - C (float): Ag concentration
    '''
    return a + (1. - a - b) / (1. + np.exp(en - en_thr) / C)


def Bsel_psurv(en, C, par):
    '''
    Ag-binding selection survival probability.

    Parameters:
    - en (array): energies for which the survival probability must be evaluated
    - C (float): Ag concentration
    - par: model parameters array
    '''
    if par['B_sel']:
        return sigmoid_psurv(en, en_thr=par['eps_B'], a=0, b=0, C=C)
    else:
        return np.ones_like(en)


def Tsel_psurv(en, bareps, C, par):
    '''
    T-cell help selection survival probability.

    Parameters:
    - en (array): energies for which the survival probability must be evaluated
    - bareps (float): population's bar-epsilon.
    - C (float): Ag concentration
    - par: model parameters array
    '''
    if par['T_sel']:
        a, b = par['a_selection'], par['b_selection']
        return sigmoid_psurv(en, en_thr=bareps, a=a, b=b, C=C)
    else:
        return np.ones_like(en)


# --- concentration evolution


def next_ag_concentration(C_av, C_res, k_minus, k_plus):
    '''
    Perform one evolution step for the available and reservoir concentrations.
    Notice that k_minus and k_plus must be in units of turns, not of days.
    '''
    nrc = C_res * np.exp(-k_plus)
    nac = C_av * np.exp(-k_minus)
    nac += C_res * (k_plus / (k_plus - k_minus)) * \
        (np.exp(-k_minus) - np.exp(-k_plus))
    return nac, nrc


# --- differentiation


def prob_mc_pc_differentiation(par, t_rounds):
    '''
    Given the set of parameters and the evolution time in rounds returns the
    probability of mc and pc differentiation.

    Parameters:
    - par: parameters dictionary
    - t_rounds (int): evolution time in rounds

    Returns:
    - p_mc, p_pc (float): probabilities of MC and PC differentiation.
    '''
    p_diff = par['diff_prob']
    days_per_round = par['days_per_turn']
    sw_t = par['diff_switch_time'] / days_per_round
    sigma_t = par['diff_switch_sigma'] / days_per_round
    residual_f = par['diff_residual_fraction']
    # if no switch time then same probability of MC/PC fate
    if sw_t is None:
        return p_diff / 2., p_diff / 2.
    # if switch time but no sigma then hard switch
    elif sigma_t is None:
        p_main, p_res = p_diff * (1. - residual_f), p_diff * residual_f
        return (p_main, p_res) if t_rounds <= sw_t else (p_res, p_main)
    # else sigmoid switch
    else:
        fr_mc = residual_f + (1. - 2. * residual_f) / \
            (1. + np.exp((t_rounds - sw_t) / sigma_t))
        return p_diff * fr_mc, p_diff * (1. - fr_mc)


# --- GC seeding (stochastic GC)

def pick_founders_en(par, mc_seed_energies):
    '''
    Utility function for determining the founder clones population of a GC.
    It takes as argument the parameter dictionary and the list of MCs
    previously collected during evolution.

    It returns the list of founder clones, randomly picked between memory and
    naive cells according to the model specifications.

    Parameters:
    - par: model parameters dictionary
    - mc_seed_energies (array): list of energies for the MCs collected so far
        in evolution.
    '''
    par_mc_reinit = par['f_mem_reinit']
    Ni = par['N_i']
    Nf = par['N_founders']
    Nmc = mc_seed_energies.size
    # evaluate probability that a clone comes from the memory pool
    if par_mc_reinit == 'pop':
        # proportional to the size of the MC population
        pr_mc = Nmc / (Ni + Nmc)
    else:
        # constant
        pr_mc = par_mc_reinit

    # pick founders among MC + Naive cells
    N_mem_founders = np.random.binomial(n=Nf, p=pr_mc)
    en_founders = np.zeros(Nf)
    # add memory founders
    en_founders[:N_mem_founders] = np.random.choice(
        mc_seed_energies, N_mem_founders, replace=Nmc < N_mem_founders)
    # add naive founders
    en_founders[N_mem_founders:] = np.random.normal(
        loc=par['mu_i'], scale=par['sigma_i'], size=Nf - N_mem_founders)

    return en_founders


# --- mutations

def generate_stoch_mutations(par, N_mut):
    '''
    Generates log-normal distributed random mutations.

    Parameters:
    - par: model parameters dictionary
    - N_mut (int): number of mutations to be generated

    Returns:
    - delta_en (array): list of energy differences caused by mutation.
    '''
    delta_en = np.random.lognormal(
        mean=par['ker_ln_mu'], sigma=par['ker_ln_sigma'],
        size=N_mut) + par['ker_ln_offset']
    return delta_en


def mutation_kernel(par):
    '''
    Builds the total mutation kernel...
    '''
    # build x-bins
    dx = par['dx']
    ker_x = np.arange(0., par['ker_xlim'], dx)
    ker_x = np.concatenate((-ker_x[:0:-1], ker_x))

    # build affinity-affecting mutations kernel (lognormal distribution)
    ker_aa = lognorm_pdf(x=ker_x, dx=dx,
                         mu=par['ker_ln_mu'],
                         sigma=par['ker_ln_sigma'],
                         offset=par['ker_ln_offset'])

    # build kernel for silent mutations (delta on zero)
    nxk = len(ker_x)
    delta = np.zeros(nxk)
    delta[nxk // 2] = 1. / dx

    # building total kernel for a single mutation
    ker_one = par['p_aa_eff'] * ker_aa + par['p_sil_eff'] * delta

    # include the effect of duplication
    ker_one *= 2

    # build total kernel for n mutations and duplication (kernel self-convolution)
    ker_tot = np.copy(ker_one)
    for m in range(par['n_duplications'] - 1):
        ker_tot = np.convolve(ker_tot, ker_one, 'same') * dx

    return ker_x, ker_tot


# --- evaluate responders population

def evaluate_responders(MC, PC, g_mem, sim_type, N_res):
    '''
    This function evaluates the population of responder cells elicited by the
    immunization scheme. It is defined as a weighted mixture of MCs and PCs,
    containing a fraction 'g_mem' of memory cells.

    Args:
    - MC, PC (stoch_pop/det_pop objects): memory and plasma cell populations
        collected during the immunization scheme
    - g_mem (float): memory cell fraction of the responder population.
    - sim_type (string): either 'stochastic' or 'deterministic', depending on
        the class of the MC/PC populations
    - N_res (int): responder population desired size.

    Returns:
    - resp_pop (stoch_pop/det_pop object): population of responder cells.
    '''
    if sim_type == 'stochastic':
        resp_pop = stoch_responders(MC, PC, g_mem, N_res=N_res)
    elif sim_type == 'deterministic':
        resp_pop = det_responders(MC, PC, g_mem, N_res=N_res)
    else:
        raise Exception('sim_type must be either stochastic or deterministic')
    return resp_pop


def stoch_responders(MC, PC, g_mem, N_res):
    '''
    Generates a mixture of MC and PC populations, with fraction g_mem of memory
    cells. A total of N_res cells is randomly picked from the two populations,
    with replacement if necessary.

    Args:
    - MC, PC: (stoch_pop objects): memory and plasma cell populations collected
        during the immunization scheme
    - g_mem (float): memory cell fraction of the responder population.
    - N_res (int): total number of cells in the responder population.

    Returns:
    - resp_pop (stoch_pop object): responder population
    '''
    # set number of MCs in resp_pop (only MCs if no PCs are present)
    if PC.N_cells() == 0:
        N_mc = N_res
    else:
        N_mc = np.round(N_res * g_mem).astype(np.int)
    # number of PCs in responding population
    N_pc = N_res - N_mc

    # extract MC and PC energies
    MC_en = np.random.choice(MC.en, size=N_mc, replace=MC.N_cells() < N_mc)
    PC_en = np.random.choice(PC.en, size=N_pc, replace=PC.N_cells() < N_pc)

    # construct responders pop
    resp_pop = MC.create_empty()
    resp_pop.en = np.concatenate([MC_en, PC_en])

    return resp_pop


def det_responders(MC, PC, g_mem, N_res):
    '''
    Generates a mixture of MC and PC populations, with fraction g_mem of memory
    cells.

    Args:
    - MC, PC: (det_pop objects): memory and plasma cell populations collected
        during the immunization scheme
    - g_mem (float): memory cell fraction of the responder population.
    - N_res (int): size of the responder population.

    Returns:
    - resp_pop (det_pop object): responder population
    '''
    # create a new population from a copy of the MC pop
    resp_pop = MC.create_copy_without_kernel()
    # set population size
    resp_pop.N = N_res
    # perform weighted average
    resp_pop.varphi = MC.varphi * g_mem + PC.varphi * (1. - g_mem)

    return resp_pop


# --- experimental limits

def prob_low_en_measurement(det_pf):
    '''
    Given a deterministic population object this function returns the
    probability of a below-experimental-threshold measurement, which is the
    probability for a measurement to have energy < low_en_exp_cutoff, given
    that its energy is detectable and therefore < high_en_exp_cutoff

    Args:
    - det_pf (det_pop object): determinisitc population function whose left
        tails must be quantified.

    Returns:
    - p_low_en (float): probability of a measurement < low_en_exp_cutoff, given
        that it is < high_en_exp_cutoff.
    '''
    # evaluate probability that a measurement is detectable
    mask_detectable = det_pf.x < high_en_exp_cutoff
    p_detectable = np.sum(det_pf.varphi[mask_detectable])
    # evaluate probability that a measurement is below low threshold
    mask_low = det_pf.x < low_en_exp_cutoff
    p_low = np.sum(det_pf.varphi[mask_low])
    return p_low / p_detectable


def resize_to_exp_limits_det(det_pf):
    '''
    Given a deterministic population function it restricts its domain between
    the experimental detection limits and renormalizes it. The restricted
    domain and distribution are returned.

    Args:
    - det_pf (det_pop object): determinisitc population function to resize and
        renormalize.

    Returns:
    - res_x (list of float): resized domain of the binding energy distribution
    - res_varphi (list of float): renormalized binding energy distribution
    '''
    # distribution domain and discretization step
    x, dx = det_pf.x, det_pf.dx
    # select the subset of the domain between the experimental sensitivity
    # limits.
    mask = (x < (high_en_exp_cutoff + dx)) & (x > (low_en_exp_cutoff - dx))
    # restrict the domain and the distribution
    res_x = x[mask]
    res_varphi = det_pf.varphi[mask]
    # renormalize the distribution on the restricted domain
    res_varphi /= (np.sum(res_varphi) * dx)
    # return the results
    return res_x, res_varphi


def apply_exp_limits_to_en_list(en_list):
    '''
    Given a list of binding energies this function applies the experimental
    measurement limits on it. It returns a copy of the list in which all
    energies higher than the experimental detection threshold
    'high_en_exp_cutoff' are removed, and all energies lower than the
    experimental threshold 'low_en_exp_cutoff' are set equal to the threshold.
    NB: the returned array could be empty.

    Args:
    - en_list (array of float): list of binding energies.

    Returns:
    - exp_en (array of float): list of binding energies with experimental
        detection limits applied.
    '''
    # set the energies of cells with en > low exp cutoff equal to the cutoff
    exp_en = np.maximum(en_list, low_en_exp_cutoff)

    # remove all the cells with energies higher than the high cutoff
    exp_en = exp_en[exp_en <= high_en_exp_cutoff]

    return exp_en

# --- high affinity fraction


def r_haff_from_en_list(en_list):
    '''
    Given a list of energy measurements, this function evaluates the
    high-affiniy fraction, defined as the fraction of cells with binding energy
    smaller than 'low_en_threshold'.

    Args:
    - en_list (array of floats): list of binding energies.

    Returns:
    - r_haff (float): high affinity fraction.
    '''
    # boolean mask: whether the binding energy is lower than the high-affinity
    # threshold
    h_aff_mask = np.array(en_list) <= low_en_threshold
    return h_aff_mask.mean()
