"""Cancer epistasis -- Jorge A. Alfaro-Murillo

This module estimates fluxes from and to combinations of mutations.

The total number of mutations can be set with the parameter n.

A combination of mutations is represented by a vector of 0's and 1's,
that has a 1 or a 0 in the i-th position if the individual has or does
not have the i-th mutation.

The likelihood is a Multinomial distribution where the parameters are
given by the relationship between the fluxes and the probability of
being at each state at time t=T, where at t=0 everyone is starting
without mutations. The observed variables are the samples, and the
priors for the fluxes are uninformative uniform distributions.

"""

import os
import numpy as np


import pymc as pm
import pytensor.tensor as tt
from pytensor.compile.ops import as_op


from theory import numbers_positive_lambdas
from theory import build_S_as_array
from theory import build_S_with_tuples
from theory import obtain_pos_lambdas_indices
from theory import order_pos_lambdas
from theory import generate_paths
from theory import epistatic_comparisons

from scipy.optimize import fsolve

from scipy.stats import chi2
from scipy.stats import multinomial

random_seed = 777
"""Random seed to feed the random generators, to be able to replicate
results."""
np.random.RandomState(np.random.MT19937(np.random.SeedSequence(random_seed)))


T = 1
"""Time at which the probabilities are finally evaluated"""


## We will compute numerically integrals according to this resolution
resolution = 10000
times = np.linspace(0, T, resolution)


## * Methods
## ** Computing methods

def establish_non_definable_lambdas(samples):
    """Establish which lambdas are actually definable. If the samples of the
    combination mutations x and y are both 0, then the flux from x
    to y is not definable.

    :type samples: numpy.ndarray
    :param samples: A one dimensional array of size 2^M containing the
        number of individuals at each state of S.

    :rtype: numpy.ndarray
    :return: One dimensional array of the same size as :const`S` with
        the computed probilities.

    """
    M = int(np.log2(len(samples))) # 2^M is the number mutation combinations

    S = build_S_as_array(M)

    non_definable = [[True
                      if np.sum(S[j] >= S[i]) == M
                      and np.sum(S[j] - S[i]) == 1
                      and i < len(S)-1 # until here positive lambdas
                      and samples[i] == 0
                      and samples[j] == 0
                      else False
                      for i in range(len(S))]
                     for j in range(len(S))]
    return np.array(non_definable)


def compute_Ps_at_T(positive_lambdas):
    """Compute for each state in S the probability of being at that state
    at time t=T if starting without mutations.

    :type positive_lambdas: numpy.ndarray
    :param positive_lambdas: One dimensional array with the lambda's
        that are positive. It should be of size equal to the number of
        True indeces in `const:positive_lambdas_indices`.

    :rtype: numpy.ndarray
    :return: One dimensional array of the same size as :const`S` with
        the computed probilities.

    """

    M = numbers_positive_lambdas.index(len(positive_lambdas))

    S = build_S_as_array(M)

    lambdas = np.zeros([len(S), len(S)])

    positive_lambdas_indices = obtain_pos_lambdas_indices(S)

    lambdas[positive_lambdas_indices] = positive_lambdas

    lambdas[np.eye(len(S), len(S), dtype=bool)] = -np.sum(lambdas, axis=0)

    Ps = np.zeros((len(S), resolution))

    Ps[0] = np.exp(lambdas[0, 0]*times)

    mutations_per_x = np.sum(S, axis=1)

    lambda_xx = lambdas[0, 0]
    lambda_xys = lambdas[mutations_per_x == 1, 0]
    lambda_yys = lambdas[mutations_per_x == 1, mutations_per_x == 1]
    Ps[mutations_per_x == 1] = (
        lambda_xys[:, np.newaxis] / (lambda_xx - lambda_yys[:, np.newaxis])
        * (np.exp(lambda_xx*times) - np.exp(lambda_yys[:, np.newaxis]*times)))

    for m in range(2, M):
        for k, y in enumerate(S):
            if np.sum(y) == m:
                y_minus_indices = np.array(
                    [True
                     if np.sum(y >= S[i]) == M
                     and np.sum(y - S[i]) == 1
                     else False
                     for i in range(len(S))])
                lambda_yy = lambdas[k, k]
                lambda_y_minus_eis_y = lambdas[k, y_minus_indices]
                ## Notice that Ps[y_minus_indices] were already computed
                to_int = np.exp(-lambda_yy*times) * Ps[y_minus_indices]
                ## trapezoidal rule for integration
                integral = ((np.cumsum(to_int, axis=-1) - to_int/2)
                            * T/resolution)
                Ps[k] = np.sum((lambda_y_minus_eis_y[:, np.newaxis]
                                * np.exp(lambda_yy*times)
                                * integral),
                               axis=0)

    ## This saves time on the last calculation, as we know that the sum
    ## of P_x(t) on x should be 1
    Ps[-1] = 1 - np.sum(Ps[:-1], axis=0)

    return Ps[:, -1]


## as_op is a decorator that transforms functions so that they can be
## used in pymc3 models. It could use the syntatic sugar:
##
## @as_op
## compute_Ps_at_T
##
## and this would have redefined compute_Ps_at_T, but I'll name a new
## function to use the normal compute_Ps_at_T as well
compute_Ps_at_T_tens = as_op(itypes=[tt.dvector],
                             otypes=[tt.dvector])(compute_Ps_at_T)


## ** Main method

def estimate_lambdas(samples, upper_bound_prior=3, draws=10000,
                     burn=500, chains=8, save_name=None, kwargs=None):
    """Run the main simulation to find the estimates of the fluxes.

    If draws is equal to 1 find the maximum a posteriori (MAP)
    estimate, else perform a Metropolis-Hastings Markov chain Monte
    Carlo to evaluate the posterior of the fluxes.

    The fluxes are assumed to have a non-informative [0,1] uniform as
    prior.

    :type samples: numpy.ndarray
    :param samples: One dimensional array with the samples. It should
        be of the same size as S, and have in each entry the number of
        individuals that have the respective mutation combination.

    :type upper_bound_prior: float
    :param upper_bound_prior: Upper bound for the uniform prior used
        for the each lambda.

    :type draws: int
    :param draws: Number of samples to draw.

    :type burn: int
    :param burn: Number of samples to be discarded as the procedure is
        tuning in. Only takes effect if `draws` > 1.

    :type chains: int
    :param chains: Number of chains to run in parallel. Only takes
        effect if `draws` > 1.

    :type kwargs: dict
    :param kwargs: Dictionary of keyword arguments to pass to the pymc3
        find_MAP function if `draws` = 1, or to the sample function if
        `draws` > 1.

    :type save_name: str or NoneType
    :param save_name: Save results under this name. If None (default),
        do not save the results. If `draws` = 1, that is, if the
        results are the MAP estimate, then the dictionary with the
        results is saved as a npz file. If `draws` > 1, then each
        chain of the trace goes inside a directory with this name.

    :rtype: pymc3.backends.base.MultiTrace or dict or tuple
    :return: Estimates of the fluxes as a MultiTrace pymc3 object or,
        if `draws`=1, a dictionary with the MAP (also MLE, because the
        priors are uniform) estimates, or if `draws`=1 and
        'return_raw' is passed to kwargs as True, a tuple with the
        full output of scipy.optimize.minimize.

    """

    M = int(np.log2(len(samples))) # 2^M is the number mutation combinations

    if M == 1 and draws == 1:
        ## then we have a close formula for the result
        ## by maximizing the log of the likelihood
        results = {'lambdas':np.array([np.log(1+samples[1]/samples[0])])}

    else:

        S = build_S_as_array(M)

        positive_lambdas_indices = obtain_pos_lambdas_indices(S)

        number_positive_lambdas = np.sum(positive_lambdas_indices)

        number_samples = np.sum(samples)

        with pm.Model():

            ## We set an uninformative prior for the lambdas:
            positive_lambdas = pm.Uniform(
                name="lambdas",
                lower=0,
                upper=upper_bound_prior,
                shape=int(number_positive_lambdas))

            Ps = compute_Ps_at_T_tens(positive_lambdas)

            likelihood = pm.Multinomial(name="samples",
                                        p=Ps,
                                        n=number_samples,
                                        observed=samples)
            if kwargs is None:
                kwargs = {}

            if draws == 1:
                results = pm.find_MAP(**kwargs)

            else:
                results = pm.sample(int(draws/chains),
                                    cores=chains,
                                    tune=burn,
                                    step=pm.Metropolis(),
                                    random_seed=random_seed,
                                    **kwargs)

    if save_name is not None:
        if draws == 1:
            np.savez(save_name, **results)
        else:
            pm.save_trace(results, save_name, overwrite=True)


    return results


## Comparison of differences between two groups

def p_value_same_lambda_xy(samples1,
                           samples2,
                           xy,
                           lambdas1_h1=None,
                           lambdas2_h1=None,
                           upper_bound_priors1=1,
                           upper_bound_priors2=1,
                           upper_bound_prior_shared=1,
                           return_lambdas_estimates=False,
                           verbose=False,
                           kwargs=None):
    """Estimate the probability of a null hypothesis of lambdas_xy
    being the same when estimating the parameters of two models
    informed by samples1 and samples2.

    :type samples1: numpy.ndarray
    :param samples1: One dimensional array with the first set of
        samples to compare. It should be of the same size as S, and
        have in each entry the number of individuals that have the
        respective mutation combination.

    :type samples2: numpy.ndarray
    :param samples2: One dimensional array with the second set of
        samples to compare. It should be of the same size as
        `samples1'.

    :type xy: tuple
    :param xy: A tuple with the somatic genotype from and to that will
        make the comparison. Each element of the tuple should be tuple
        of 0s and 1s of size M, where the samples are divided into
        combinations of mutations in M genes (or categories).

    :type lambdas1_h1: numpy.ndarray or NoneType
    :param lambdas1_h1: Estimates of the fluxes informed by `samples1'
        under H_1. It is optional, otherwise they are computed.

    :type lambdas2_h1: numpy.ndarray or NoneType
    :param lambdas2_h1: Estimates of the fluxes informed by `samples2'
        under H_1. It is optional, otherwise they are computed.

    :type upper_bound_priors1: float
    :param upper_bound_priors1: Upper bound for the uniform prior used
        for each lambda for the model informed by samples1.

    :type upper_bound_priors2: float
    :param upper_bound_priors2: Upper bound for the uniform prior used
        for each lambda for the model informed by samples2.

    :type upper_bound_shared: float
    :param upper_bound_shared: Upper bound for the uniform prior used
        for the shared lambda.

    :type kwargs: dict
    :param kwargs: Dictionary of keyword arguments to pass to the pymc3
        find_MAP function.

    :type return_lambdas_estimates: bool
    :param return_lambdas_estimates: If True, return the lambda
        estimates H_1

    :type verbose: bool
    :param verbose: If True, print information on the lambda and
        loglikehood estimates under H_0 and H_1

    :rtype: numpy.float64 or tuple
    :return: p-value of the null hypothesis, or if
        return_lambdas_estimates is True, a tuple with the p-value,
        and the lambdas estimates for samples1 and samples2 under the
        H0 hypothesis.

    """

    if len(samples1) != len(samples2):
        raise ValueError("The lengths of samples1 and samples2 must be equal.")

    M = int(np.log2(len(samples1))) # 2^M is the number mutation combinations
    S = build_S_as_array(M)

    positive_lambdas_indices =  obtain_pos_lambdas_indices(S)
    ordered_pos_lambdas = order_pos_lambdas(S)
    xy_index = ordered_pos_lambdas.index(xy)

    number_samples1 = int(np.sum(samples1))
    number_samples2 = int(np.sum(samples2))


    ## Our H_0 is that the specific lambda_xy is equal for both
    ## models. So in the likelihood with H0 is the product of the
    ## likelihoods of the models with samples1 aand samples2 but
    ## fixing a shared lambda. First, we compute the lambdas that
    ## would maximize such a likelihood


    if M == 1:
        ## then we have a close formula for the result, and only one
        ## lambda, the shared one
        results_map_h0 = {'lambda_shared':np.array(
            [np.log(1+(samples1[1]+samples2[1])/(samples1[0]+samples2[0]))])}

        lambdas1_h0 = results_map_h0['lambda_shared']
        lambdas2_h0 = results_map_h0['lambda_shared']


    else:
        ## when M > 1 then there is one shared lambda and other
        ## lambdas that are model specific

        number_positive_lambdas = int(np.sum(positive_lambdas_indices)-1)

        with pm.Model():

            ## We set an uninformative prior for the lambdas:
            positive_lambdas1 = pm.Uniform(
                name="lambdas1",
                lower=0,
                upper=upper_bound_priors1,
                shape=number_positive_lambdas)

            positive_lambdas2 = pm.Uniform(
                name="lambdas2",
                lower=0,
                upper=upper_bound_priors2,
                shape=number_positive_lambdas)

            lambda_shared = pm.Uniform(
                name="lambda_shared",
                lower=0,
                upper=upper_bound_prior_shared,
                shape=1)

            concatenated_lambdas1 = tt.concatenate([
                positive_lambdas1[:xy_index],
                lambda_shared,
                positive_lambdas1[xy_index:]])

            concatenated_lambdas2 = tt.concatenate([
                positive_lambdas2[:xy_index],
                lambda_shared,
                positive_lambdas2[xy_index:]])

            Ps1 = compute_Ps_at_T_tens(concatenated_lambdas1)
            Ps2 = compute_Ps_at_T_tens(concatenated_lambdas2)

            likelihood1 = pm.Multinomial(name="samples1",
                                         p=Ps1,
                                         n=number_samples1,
                                         observed=samples1)
            likelihood2 = pm.Multinomial(name="samples2",
                                         p=Ps2,
                                         n=number_samples2,
                                         observed=samples2)
            if kwargs is None:
                kwargs = {}

            results_map_h0 = pm.find_MAP(**kwargs)

        lambdas1_h0 = np.concatenate([
            results_map_h0['lambdas1'][:xy_index],
            results_map_h0['lambda_shared'],
            results_map_h0['lambdas1'][xy_index:]])


        lambdas2_h0 = np.concatenate([
            results_map_h0['lambdas2'][:xy_index],
            results_map_h0['lambda_shared'],
            results_map_h0['lambdas2'][xy_index:]])

    if verbose:
        print(f"lambdas1 under H_0: {lambdas1_h0}")
        print(f"lambdas2 under H_0: {lambdas2_h0}")


    ## Now that we have the lambdas that maximize the likelihood under
    ## H0, we compute the actual likelihood

    logp1_h0 = multinomial.logpmf(samples1, number_samples1, compute_Ps_at_T(lambdas1_h0))
    logp2_h0 = multinomial.logpmf(samples2, number_samples2, compute_Ps_at_T(lambdas2_h0))
    logp_h0 = logp1_h0 + logp2_h0

    if verbose:
        print(f"loglikelihood for samples1 under h0: {logp1_h0}")
        print(f"loglikelihood for samples2 under h0: {logp2_h0}")


    ## Our H_1 is that the lambda_xy is different, which would come
    ## from separate MLE calculation with samples1 and samples2

    ## If lambdas1 or lambdas2 are not provided we need to compute them
    if lambdas1_h1 is None:
        lambdas1_h1 = estimate_lambdas(samples1, draws=1,
                                       ## this will fail if bounds are
                                       ## set specifically for each lambda
                                       upper_bound_prior=upper_bound_priors1)['lambdas']
        if verbose:
            print(f"lambdas1 under H_1: {lambdas1_h1}")

    if lambdas2_h1 is None:
        lambdas2_h1 = estimate_lambdas(samples2, draws=1,
                                       ## this will fail if bounds are
                                       ## set specifically for each lambda
                                       upper_bound_prior=upper_bound_priors2)['lambdas']
        if verbose:
            print(f"lambdas2 under H_1: {lambdas2_h1}")

    logp1_h1 = multinomial.logpmf(samples1, number_samples1, compute_Ps_at_T(lambdas1_h1))
    logp2_h1 = multinomial.logpmf(samples2, number_samples2, compute_Ps_at_T(lambdas2_h1))
    logp_h1 = logp1_h1 + logp2_h1
    if verbose:
        print(f"loglikehood for samples1 under h1: {logp1_h1}")
        print(f"loglikehood for samples2 under h1: {logp2_h1}")

    ## Now we compute the probability of H0 using Wilk's theorem
    D = 2*(logp_h1-logp_h0)
    if return_lambdas_estimates:
        return (chi2.sf(D, 1), lambdas1_h1, lambdas2_h1)
    else:
        return chi2.sf(D, 1)



def p_value_same_gamma_xy(samples1,
                          samples2,
                          xy,
                          mu1_y_minux_x,
                          mu2_y_minux_x,
                          lambdas1_h1=None,
                          lambdas2_h1=None,
                          upper_bound_priors1=1,
                          upper_bound_priors2=1,
                          upper_bound_prior_shared=1,
                          return_lambdas_estimates=False,
                          verbose=False,
                          kwargs=None):
    """Estimate the probability of a null hypothesis of gammas_xy
    being the same when estimating the parameters of two models
    informed by samples1 and samples2.

    :type samples1: numpy.ndarray
    :param samples1: One dimensional array with the first set of
        samples to compare. It should be of the same size as S, and
        have in each entry the number of individuals that have the
        respective mutation combination.

    :type samples2: numpy.ndarray
    :param samples2: One dimensional array with the second set of
        samples to compare. It should be of the same size as
        `samples1'.

    :type xy: tuple
    :param xy: A tuple with the somatic genotype from and to that will
        make the comparison. Each element of the tuple should be tuple
        of 0s and 1s of size M, where the samples are divided into
        combinations of mutations in M genes (or categories).

    :type mu1_y_minux_x: float
    :param mu1_y_minux_x: Mutation rate for `samples1' of the mutation
        gained by moving to the somatic genotype y when coming from
        the genotype x.

    :type mu2_y_minux_x: float
    :param mu2_y_minux_x: Mutation rate for `samples2' of the mutation
        gained by moving to the somatic genotype y when coming from
        the genotype x.

    :type lambdas1_h1: numpy.ndarray or NoneType
    :param lambdas1_h1: Estimates of the fluxes informed by `samples1'
        under H_1. It is optional, otherwise they are computed.

    :type lambdas2_h1: numpy.ndarray or NoneType
    :param lambdas2_h1: Estimates of the fluxes informed by `samples2'
        under H_1. It is optional, otherwise they are computed.

    :type upper_bound_priors1: float
    :param upper_bound_priors1: Upper bound for the uniform prior used
        for each lambda for the model informed by samples1.

    :type upper_bound_priors2: float
    :param upper_bound_priors2: Upper bound for the uniform prior used
        for each lambda for the model informed by samples2.

    :type upper_bound_shared: float
    :param upper_bound_shared: Upper bound for the uniform prior used
        for the shared lambda.

    :type kwargs: dict
    :param kwargs: Dictionary of keyword arguments to pass to the pymc3
        find_MAP function.

    :type return_lambdas_estimates: bool
    :param return_lambdas_estimates: If True, return the lambda
        estimates H_1

    :type verbose: bool
    :param verbose: If True, print information on the lambda and
        loglikehood estimates under H_0 and H_1

    :rtype: numpy.float64 or tuple
    :return: p-value of the null hypothesis, or if
        return_lambdas_estimates is True, a tuple with the p-value,
        and the lambdas estimates for samples1 and samples2 under the
        H0 hypothesis.

    """

    if len(samples1) != len(samples2):
        raise ValueError("The lengths of samples1 and samples2 must be equal.")

    M = int(np.log2(len(samples1))) # 2^M is the number mutation combinations
    S = build_S_as_array(M)

    positive_lambdas_indices =  obtain_pos_lambdas_indices(S)
    ordered_pos_lambdas = order_pos_lambdas(S)
    xy_index = ordered_pos_lambdas.index(xy)

    number_samples1 = int(np.sum(samples1))
    number_samples2 = int(np.sum(samples2))


    ## Our H_0 is that the specific gamma_xy is equal for both
    ## models. So in the likelihood with H0 is the product of the
    ## likelihoods of the models with samples1 and samples2 but
    ## fixing a shared gamma. First, we compute the lambdas that
    ## would maximize such a likelihood


    if M == 1:
        ## then the two lambdas are determined by one shared gamma and
        ## we have an equation in one variable to solve for when
        ## maximizing the likelihood (but no analytic solution)
        def equation_to_solve(lambda1, mu2_over_mu1, n1, N1_minus_n1, n2, N2_minus_n2):
            return (+ n1 / (np.exp(lambda1) - 1)
                    + mu2_over_mu1 * n2 / (np.exp(mu2_over_mu1 * lambda1) - 1)
                    - N1_minus_n1 - mu2_over_mu1 * N2_minus_n2)

        lambda1_initial_guess = 0.2

        lambda1_from_gamma_shared = fsolve(equation_to_solve,
                                           lambda1_initial_guess,
                                           args=(mu2_y_minux_x/mu1_y_minux_x,
                                                 samples1[1],  # n1
                                                 samples1[0],  # N1-n1
                                                 samples2[1],  # n2
                                                 samples2[0])) # N2-n2

        lambdas1_h0 = lambda1_from_gamma_shared

        lambdas2_h0 = lambda1_from_gamma_shared * mu2_y_minux_x/mu1_y_minux_x

        gamma_xy_h0 = lambda1_from_gamma_shared/mu1_y_minux_x

    else:
        ## when M > 1 then there is one shared gamma that determines
        ## two lambdas, and other lambdas that are model specific

        number_positive_lambdas = int(np.sum(positive_lambdas_indices)-1)

        with pm.Model():

            ## We set an uninformative prior for the lambdas:
            positive_lambdas1 = pm.Uniform(
                name="lambdas1",
                lower=0,
                upper=upper_bound_priors1,
                shape=number_positive_lambdas)

            positive_lambdas2 = pm.Uniform(
                name="lambdas2",
                lower=0,
                upper=upper_bound_priors2,
                shape=number_positive_lambdas)

            lambda1_from_gamma_shared = pm.Uniform(
                name="lambda1_from_gamma_shared",
                lower=0,
                upper=upper_bound_prior_shared,
                shape=1)

            lambda2_from_gamma_shared = lambda1_from_gamma_shared/mu1_y_minux_x * mu2_y_minux_x

            concatenated_lambdas1 = tt.concatenate([
                positive_lambdas1[:xy_index],
                lambda1_from_gamma_shared,
                positive_lambdas1[xy_index:]])

            concatenated_lambdas2 = tt.concatenate([
                positive_lambdas2[:xy_index],
                lambda2_from_gamma_shared,
                positive_lambdas2[xy_index:]])

            Ps1 = compute_Ps_at_T_tens(concatenated_lambdas1)
            Ps2 = compute_Ps_at_T_tens(concatenated_lambdas2)

            likelihood1 = pm.Multinomial(name="samples1",
                                         p=Ps1,
                                         n=number_samples1,
                                         observed=samples1)
            likelihood2 = pm.Multinomial(name="samples2",
                                         p=Ps2,
                                         n=number_samples2,
                                         observed=samples2)
            if kwargs is None:
                kwargs = {}

            results_map_h0 = pm.find_MAP(**kwargs)

        lambdas1_h0 = np.concatenate([
            results_map_h0['lambdas1'][:xy_index],
            results_map_h0['lambda1_from_gamma_shared'],
            results_map_h0['lambdas1'][xy_index:]])


        lambdas2_h0 = np.concatenate([
            results_map_h0['lambdas2'][:xy_index],
            results_map_h0['lambda1_from_gamma_shared']/mu1_y_minux_x * mu2_y_minux_x,
            results_map_h0['lambdas2'][xy_index:]])

        gamma_xy_h0 = results_map_h0['lambda1_from_gamma_shared']/mu1_y_minux_x

    if verbose:
        print(f"lambdas1 under H_0: {lambdas1_h0}")
        print(f"lambdas2 under H_0: {lambdas2_h0}")
        print(f"gamma_xy under H_0: {gamma_xy_h0}")


    ## Now that we have the lambdas that maximize the likelihood under
    ## H0, we compute the actual likelihood

    logp1_h0 = multinomial.logpmf(samples1, number_samples1, compute_Ps_at_T(lambdas1_h0))
    logp2_h0 = multinomial.logpmf(samples2, number_samples2, compute_Ps_at_T(lambdas2_h0))
    logp_h0 = logp1_h0 + logp2_h0

    if verbose:
        print(f"loglikelihood for samples1 under h0: {logp1_h0}")
        print(f"loglikelihood for samples2 under h0: {logp2_h0}")


    ## Our H_1 is that the gamma_xy is different, which would come
    ## from separate MLE calculation with samples1 and samples2

    ## If lambdas1 or lambdas2 are not provided we need to compute them
    if lambdas1_h1 is None:
        lambdas1_h1 = estimate_lambdas(samples1, draws=1,
                                       ## this will fail if bounds are
                                       ## set specifically for each lambda
                                       upper_bound_prior=upper_bound_priors1)['lambdas']
        if verbose:
            print(f"lambdas1 under H_1: {lambdas1_h1}")

    if lambdas2_h1 is None:
        lambdas2_h1 = estimate_lambdas(samples2, draws=1,
                                       ## this will fail if bounds are
                                       ## set specifically for each lambda
                                       upper_bound_prior=upper_bound_priors2)['lambdas']
        if verbose:
            print(f"lambdas2 under H_1: {lambdas2_h1}")

    logp1_h1 = multinomial.logpmf(samples1, number_samples1, compute_Ps_at_T(lambdas1_h1))
    logp2_h1 = multinomial.logpmf(samples2, number_samples2, compute_Ps_at_T(lambdas2_h1))
    logp_h1 = logp1_h1 + logp2_h1
    if verbose:
        print(f"loglikehood for samples1 under h1: {logp1_h1}")
        print(f"loglikehood for samples2 under h1: {logp2_h1}")

    ## Now we compute the probability of H0 using Wilk's theorem
    D = 2*(logp_h1-logp_h0)
    if return_lambdas_estimates:
        return (chi2.sf(D, 1), lambdas1_h1, lambdas2_h1)
    else:
        return chi2.sf(D, 1)


def p_value_gamma_xy_equal_1(samples,
                             xy,
                             mu_y_minus_x,
                             lambdas_h1=None,
                             upper_bound_priors=1,
                             return_lambdas_estimates=False,
                             verbose=False,
                             kwargs=None):
    """Estimate the probability of a null hypothesis of a particular
    gamma_xy being 1 (that is, there is no selection positive or
    negative).

    :type samples: numpy.ndarray
    :param samples: One dimensional array with the samples. It should
        be of the same size as S, and have in each entry the number of
        individuals that have the respective mutation combination.

    :type xy: tuple
    :param xy: A tuple with the somatic genotype from and to that will
        be compared. Each element of the tuple should be tuple of 0s
        and 1s of size M, where the samples are divided into
        combinations of mutations in M genes (or categories).

    :type mu_y_minux_x: float
    :param mu_y_minux_x: Mutation rate of the mutation gained by
        moving to the somatic genotype y when coming from the genotype
        x.

    :type lambdas_h1: numpy.ndarray or NoneType
    :param lambdas_h1: Estimates of the fluxes informed by `samples'
        under H_1. It is optional, otherwise they are computed.

    :type upper_bound_priors: float
    :param upper_bound_priors: Upper bound for the uniform prior used
        for each lambda for the model.

    :type kwargs: dict
    :param kwargs: Dictionary of keyword arguments to pass to the pymc3
        find_MAP function.

    :type return_lambdas_estimates: bool
    :param return_lambdas_estimates: If True, return the lambda
        estimates H_1

    :type verbose: bool
    :param verbose: If True, print information on the lambda and
        loglikehood estimates under H_0 and H_1

    :rtype: numpy.float64 or tuple
    :return: p-value of the null hypothesis, or if
        return_lambdas_estimates is True, a tuple with the p-value,
        and the lambdas estimates under the H0 hypothesis.

    """

    M = int(np.log2(len(samples))) # 2^M is the number mutation combinations
    S = build_S_as_array(M)

    positive_lambdas_indices =  obtain_pos_lambdas_indices(S)
    ordered_pos_lambdas = order_pos_lambdas(S)
    xy_index = ordered_pos_lambdas.index(xy)

    number_samples = int(np.sum(samples))

    ## Our H_0 is that the specific gamma_xy is equal to 1. First, we
    ## compute the lambdas that would maximize such the likelihood
    ## under H0


    if M == 1:
        ## then the there is no extra lambda to compute
        lambdas_h0 = [mu_y_minux_x]

    else:
        ## when M > 1 then one lambda is equal to mu, and other
        ## lambdas that are model specific

        number_positive_lambdas = int(np.sum(positive_lambdas_indices)-1)

        with pm.Model():

            ## We set an uninformative prior for the lambdas:
            positive_lambdas = pm.Uniform(
                name="lambdas",
                lower=0,
                upper=upper_bound_priors,
                shape=number_positive_lambdas)

            concatenated_lambdas = tt.concatenate([
                positive_lambdas[:xy_index],
                tt.as_tensor_variable([mu_y_minus_x]),
                positive_lambdas[xy_index:]])

            Ps = compute_Ps_at_T_tens(concatenated_lambdas)

            likelihood = pm.Multinomial(name="samples",
                                        p=Ps,
                                        n=number_samples,
                                        observed=samples)
            if kwargs is None:
                kwargs = {}

            results_map_h0 = pm.find_MAP(**kwargs)

        lambdas_h0 = np.concatenate([
            results_map_h0['lambdas'][:xy_index],
            [mu_y_minus_x],
            results_map_h0['lambdas'][xy_index:]])

    if verbose:
        print(f"lambdas under H_0: {lambdas_h0}")


    ## Now that we have the lambdas that maximize the likelihood under
    ## H0, we compute the actual likelihood

    logp_h0 = multinomial.logpmf(samples, number_samples, compute_Ps_at_T(lambdas_h0))

    if verbose:
        print(f"loglikelihood under h0: {logp_h0}")


    ## Our H_1 is that the gamma_xy is different, which would come
    ## from separate MLE calculation with samples1 and samples2

    ## If lambdas are not provided we need to compute them
    if lambdas_h1 is None:
        lambdas_h1 = estimate_lambdas(samples, draws=1,
                                      ## this will fail if bounds are
                                      ## set specifically for each lambda
                                      upper_bound_prior=upper_bound_priors)['lambdas']
        if verbose:
            print(f"lambdas under H_1: {lambdas_h1}")

    logp_h1 = multinomial.logpmf(samples, number_samples, compute_Ps_at_T(lambdas_h1))
    if verbose:
        print(f"loglikehood under h1: {logp_h1}")

    ## Now we compute the probability of H0 using Wilk's theorem
    D = 2*(logp_h1-logp_h0)
    if return_lambdas_estimates:
        return (chi2.sf(D, 1), lambdas_h1)
    else:
        return chi2.sf(D, 1)



## ** Other methods to handle the results

def convert_samples_to_dict(samples):
    """Convert a samples array to a dictionary.

    :type samples: list
    :param samples: Number of patients in each mutation combination
        (as returned by :func:`count_pts_per_combination`).

    :rtype: dict
    :return: Dictionary with the samples, indexed by tuples of 1's and
        0's representing whether the mutation occur fot the gene or
        not.

    """

    M = int(np.log2(len(samples))) # 2^M is the number mutation combinations

    S = build_S_as_array(M)

    results_as_dict = {tuple(x):value
                       for x, value in zip(S, samples)}

    return results_as_dict


def convert_mus_to_dict(mus, genes):
    """Convert a dictionary of mus with all mus per gene, to another
    dictionary with keys equal to each possible jumps under the
    assumption that the mus are not genotype dependent.

    :type mus: dict
    :param mus: Dictionary with mutation rates per gene.

    :type genes: list
    :param genes: List with the genes to extract.

    :rtype: dict
    :return: Dictionary with the mutation rates, indexed by a pair of
        tuples representing the mutation combination where the flux is
        coming from and going to.

    """

    M = len(genes)

    S = build_S_as_array(M)

    xys = order_pos_lambdas(S)

    results_as_dict = {xy:mus[genes[list(np.array(xy[1])-np.array(xy[0])).index(1)]]
                       for xy in xys}

    return results_as_dict


def convert_lambdas_to_dict(results):
    """Convert a fluxes array to a dictionary.

    The indexes of the dictionary are the subscripts of the
    lambdas. So if lambda[x, y] represents the flux from x to y, it
    can be obtain from the result of this function easily by using
    output[x, y], with x and y written as tuples.

    :type results: pymc3.backends.base.MultiTrace or dict or list
    :param results: Estimates of the fluxes as a MultiTrace pymc3
        object or a dictionary with the ordered MAP estimates (both as
        returned by :func:`estimate_lambdas`), or as a list of
        confidence intervals for the fluxes as returned by
        asymp_CI_lambdas.

    :rtype: dict
    :return: Dictionary with the fluxes, indexed by a pair of tuples
        representing the mutation combination where the flux is coming
        from and going to.

    """
    if isinstance(results, dict):
        # if results is a dictionary, it must be the return from
        # estimate_lambdas with draws=1, that is, the MAP (and MLE)
        results = results['lambdas']
    elif isinstance(results, pm.backends.base.MultiTrace):
        # if results is a pymc3.backends.base.MultiTrace, it must be
        # the return from estimate_lambdas with draws>1, that is, the
        # posterior estimates
        results = results.get_values('lambdas').T
    elif not isinstance(results, list):
        # otherwise results should be the confidence intervales and
        # thus should be a list because that is the return of
        # asymp_CI_lambdas
        raise Exception("`results` type not compatible")

    M = numbers_positive_lambdas.index(len(results))

    S = build_S_as_array(M)

    subscripts = order_pos_lambdas(S)

    results_as_dict = {subscript_pair:value
                       for subscript_pair, value in zip(subscripts,
                                                        results)}

    return results_as_dict


def compute_gammas(lambdas, mus):
    """Compute the selection coefficients.

    :type lambdas: dict
    :param lambdas: Dictionary with the fluxes, indexed by a pair of
        tuples representing the mutation combination where the flux is
        coming from and going to.

    :type mus: dict
    :param mus: Dictionary with the mutation rates indexed by a tuple
        of 0's and one 1, representing which mutation the respective
        rate represents. It is assumed that the mutation rates do not
        change via epistatic effects

    :rtype: dict
    :return: Dictionary with the selection coefficients, indexed by a
        pair of tuples representing the mutation combination where the
        flux is coming from and going to.

    """
    gammas = {
        xy:flux/mus[tuple(np.array(xy[1])-np.array(xy[0]))]
        for xy, flux in lambdas.items()}
    return gammas


def compute_CI_gamma(lambda_cis, mus):
    gamma_CIs = {
        xy:[bound/mus[tuple(np.array(xy[1])-np.array(xy[0]))] for bound in ci]
        for xy, ci in lambda_cis.items()}
    return gamma_CIs


def compute_log_lh(positive_lambdas, samples):
    """Not really the log likelihood but that plus a constant (that
    is, the log of number that is proportional to the likelihood).

    :type positive_lambdas: numpy.ndarray
    :param positive_lambdas: One dimensional array with the lambda's
        that are positive. It should be of size equal to the number of
        True indeces in `const:positive_lambdas_indices`.

    :type samples: numpy.ndarray
    :param samples: A one dimensional array of size 2^M containing the
        number of individuals at each state of S.

    :rtype:
    :return: Dictionary with the fluxes, indexed by a pair of tuples
        representing the mutation combination where the flux is coming
        from and going to.

    """
    Ps = compute_Ps_at_T(positive_lambdas)
    return np.sum(samples*np.log(Ps))


def asymp_CI_lambda(lambda_index, lambdas_mle, samples, ci=0.95, tol=10**(-6)):
    """Compute an asymptomatic confidence interval for the flux indexed by
    `lambda_index`.

    An asymptotic confidence interval for a particular lambda is given by:

    {lambda | log_lh(lambda, rest of lambdas_mle) >= log_lh(lambdas_mle) - Chi^2(1,ci)/2}

    where Chi^2(1,ci) is the ci-quantile of the Chi-squared distribution with 1 dof.

    (proof at https://statproofbook.github.io/P/ci-wilks.html)

    """
    rhs = (compute_log_lh(lambdas_mle, samples) - chi2.ppf(ci, 1)/2)

    upper_est = lambdas_mle[lambda_index]
    to_compare = 2*upper_est

    lambdas = lambdas_mle.copy()
    lambdas[lambda_index] = to_compare

    diff = compute_log_lh(lambdas, samples) - rhs

    while np.abs(diff) > tol:
        if diff > 0:
            ## lambdas in CI
            upper_est = to_compare
            to_compare = 2*to_compare
        else:
            ## lambdas not in CI
            to_compare = (to_compare + upper_est)/2

        lambdas[lambda_index] = to_compare
        diff = compute_log_lh(lambdas, samples) - rhs

    upper_est = to_compare


    lower_est = lambdas_mle[lambda_index]
    if lower_est > tol:
        ## otherwise it's practically 0
        to_compare = lower_est/2

        lambdas[lambda_index] = to_compare

        diff = compute_log_lh(lambdas, samples) - rhs

        while np.abs(diff) > tol and to_compare > tol:
            if diff > 0:
                ## lambdas in CI
                lower_est = to_compare
                to_compare = to_compare/2
            else:
                ## lambdas not in CI
                to_compare = (to_compare + lower_est)/2

            lambdas[lambda_index] = to_compare
            diff = compute_log_lh(lambdas, samples) - rhs

        lower_est = to_compare

    return [lower_est, upper_est]


def asymp_CI_lambdas(lambdas_mle, samples, ci=0.95, tol=10**(-6),
                     print_progress=False):

    CIs = []

    for lambda_index in range(len(lambdas_mle)):
        if print_progress:
            print(f"Computing CI for lambda {lambda_index}/{len(lambdas_mle)}")
        CIs.append(asymp_CI_lambda(lambda_index,
                                   lambdas_mle,
                                   samples,
                                   ci=ci,
                                   tol=tol))

    return CIs


## * Probability of each trajectory

def compute_probability_paths(lambdas,
                              destiny='all'):

    """Compute the probability of each possible sequence of mutations
    to lead to a specific somatic genotype.

    :type lambdas: dict
    :param lambdas: A dictionary with the fluxes estimates with keys
        being tuples representing the somatic genotypes from and to of
        the respective flux (as for example the output of
        :func:`estimate_lambdas` with draws=1).

    :type destiny: tuple or str
    :param destiny: The destination somatic genotype as a tuple, or
        'all' (default) which would give results for each of the
        possible destination somatic genotypes.

    :rtype: dict
    :return: Dictionary indexed by the paths with values equal to
        their respective probability of ending up a the somatic
        genotype specified by `destiny`. If `destiny` is 'all', then a
        meta dictionary indexed by every somatic genotype with more
        than two mutated genes.

    """

    if destiny == 'all':

        M = len(next(iter(lambdas))[0]) # any tuple of any tuple in
                                        # the keys of lambdas has
                                        # exactly M entries

        S_except_trivial = [x for x in build_S_with_tuples(M)
                            if np.sum(x) > 1]

        S_except_trivial = sorted(S_except_trivial, key=sum)

        return {x:compute_probability_paths(lambdas=lambdas,
                                            destiny=x)
                for x in S_except_trivial}

    else:

        paths = generate_paths(destiny)

        prod_lambdas = {tuple(path):np.prod([lambdas[xy] for xy in path])
                        for path in paths}

        sum_prod_lambdas = np.sum(list(prod_lambdas.values()))

        probs = {path:value/sum_prod_lambdas
                 for path, value in prod_lambdas.items()}

        return probs


## Ordering of results

def find_model(genes, main_gene_list):
    """If results where produced from `main_gene_list' then for any
    non-ordered tuple or list `genes', return the rightly ordered
    tuple that will index the results.

    For example, if main_gene_list = ['TP53', 'KRAS', 'EGFR', 'BRAF', ...]

    and ('TP53', 'BRAF', 'EGFR') is provided as genes

    it will return

    ('TP53', 'EGFR', 'BRAF')

    """

    gene_positions = {gene: index for index, gene in enumerate(main_gene_list)}

    # Sort the genes based on their positions in the main_gene_list
    sorted_genes = sorted(genes, key=lambda gene: gene_positions[gene])

    # Convert the sorted list back to a tuple and return
    return tuple(sorted_genes)



def order_genes_by_result_values(results):
    """Give a list of the genes ordered by values.

    :type results: dict
    :param results: Dictionary with the results for which will be
        ordered by result values. It should be indexed by tuples of
        genes that were included in the model, and where here we will
        only take the M=1 (no epistasis) results for the ordering. Use
        for example as input: load_results('selections')['pan_data']
        to order by selection, where load_results is the function from
        load_results.py

    :rtype: list
    :return: A list with the genes, ordered by values in `results'.

    """

    gene_list = {gene[0]:max(value.values()) # max doesn't
                                             # matter there
                                             # is only one
                                             # value because
                                             # M = 1
                 for gene, value in results.items()
                 if len(gene) == 1}
    gene_list = sorted(gene_list,
                       key=lambda k: gene_list[k],
                       reverse=True)

    return gene_list


def epistatic_ratios(results, M, results_cis=None):

    """Compute ratios between values of the `results' in a somatic
    genotype vs another one with one less mutated gene.


    :type results: dict
    :param results: Dictionary with the results for which we will
        compute epistatic ratios. It should be indexed by tuples
        of genes that were included in the model. Here we can put
        either selections (gammas), fluxes (lambdas, but Jeff doesn't
        like this approach), or later even mutation rates (mus).

    :type M: int
    :param M: Number of genes to consider in the models

    :type results_cis: dict or None
    :param results_cis: Dictionary with the results confidence
       intervals. Results of ratios are set to 1 if they are not
       statistically significant. Statistical significance is
       conservatively determined by non overlapping confidence
       intervals. If None is provided (default), then do consider
       statistical significance for the ratios.



    """

    comparisons = epistatic_comparisons(M)

    if results_cis is None:
        ratios = {genes:{comparison:(
            results[genes][comparison[1]]/results[genes][comparison[0]])
                         for comparison in comparisons}
                  for genes in results.keys() if len(genes) == M}
    else:
        ratios = {genes:{comparison:(
            results[genes][
                comparison[1]]/results[genes][comparison[0]])
                         if (max(results_cis[genes][comparison[0]][0],
                                 results_cis[genes][comparison[1]][0]) >
                             min(results_cis[genes][comparison[0]][1],
                                 results_cis[genes][comparison[1]][1]))
                         else 1
                         for comparison in comparisons}
                  for genes in results.keys() if len(genes) == M}

    return ratios


def epistatic_ratios_2_matrix(results, results_cis, genes_ordered_list):

    ratios_dict = epistatic_ratios(results, 2, results_cis)

    n = len(genes_ordered_list)

    matrix = np.array(n*[n*[np.nan]])

    for genes, values in ratios_dict.items():
        for comparison, value in values.items():
            gene_selected = genes[comparison[0][1].index(1)]
            context_gene = genes[comparison[1][0].index(1)]
            matrix[genes_ordered_list.index(gene_selected)][
                genes_ordered_list.index(context_gene)] = value

    return matrix



def epistatic_ratios_3rd_gene_effects(results, results_cis):
    """When does the presence of a third gene in a somatic genotype,
    alter the epistatic effect of the selection of a gene in the
    presence of another.

    Results are a dictionary where the keys are of the form

        ('gene_always_previously_mutated_in_comparison', 'gene_selected', 'third_gene')

    and the values are the gene ratios of

        gene_always_there+third_gene -> gene_always_there+third_gene+gene_selected

    and

       gene_always_there -> gene_always_there+gene_selected

    """

    ratios_dict = epistatic_ratios(results, 3, results_cis)

    effects = {}

    for genes, values in ratios_dict.items():
        for comparison, value in values.items():
            if value != 1:
                from_gene = genes[(comparison[0][0]).index(1)]
                to_gene_tuple = tuple(np.array(comparison[0][1]) -
                                      np.array(comparison[0][0]))
                to_gene = genes[to_gene_tuple.index(1)]
                gene_3rd_tuple = tuple(np.array(comparison[1][0]) -
                                       np.array(comparison[0][0]))
                gene_3rd = genes[gene_3rd_tuple.index(1)]
                effects[(from_gene, to_gene, gene_3rd)] = value


    return effects


# from load_results import load_results
# gammas_ = load_results('selections')
# gammas_cis_ = load_results('selections', 'cis')
