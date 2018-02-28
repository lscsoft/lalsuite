from __future__ import division
import numpy
import bisect
import sys

from glue.ligolw import lsctables
from glue.ligolw import dbtables
from pylal import rate


def margLikelihoodMonteCarlo(VTs, lambs, mu, mcerrs=None):
    '''
    This function marginalizes the loudest event likelihood over unknown
    Monte Carlo errors, assumed to be independent between each experiment.
    '''
    if mcerrs is None:
        mcerrs = [0]*len(VTs)

    # combine experiments, propagating statistical uncertainties
    # in the measured efficiency
    likely = 1
    for vA,dvA,mc in zip(VTs,lambs,mcerrs):
        if mc == 0:
            # we have perfectly measured our efficiency in this mass bin
            # so the posterior is given by eqn (11) in Biswas et al. [arXiv:0710.0465]
            likely *= (1+mu*vA*dvA)*numpy.exp(-mu*vA)
        else:
            # we have uncertainty in our efficiency in this mass bin and
            # want to marginalize it out using eqn (24) of Biswas et al.
            k = (vA/mc)**2 # k is 1./fractional_error**2
            likely *= (1+mu*vA*(1/k+dvA))*(1+mu*vA/k)**(-(k+1))

    return likely


def margLikelihood(VTs, lambs, mu, calerr=0, mcerrs=None):
    '''
    This function marginalizes the loudest event likelihood over unknown
    Monte Carlo and calibration errors. The vector VTs is the sensitive
    volumes for independent searches and lambs is the vector of loudest
    event likelihood. The statistical errors are assumed to be independent
    between each experiment while the calibration errors are applied
    the same in each experiment.
    '''
    if calerr == 0:
        return margLikelihoodMonteCarlo(VTs,lambs,mu,mcerrs)

    std = calerr
    mean = 0 # median volume = 1

    fracerrs = numpy.linspace(0.33,3,1e2) # assume we got the volume to a factor of three or better
    errdist = numpy.exp(-(numpy.log(fracerrs)-mean)**2/(2*std**2))/(fracerrs*std) # log-normal pdf
    errdist /= errdist.sum() #normalize

    likely = sum([ pd*margLikelihoodMonteCarlo(delta*VTs,lambs,mu,mcerrs) for delta, pd in zip(fracerrs,errdist)]) #marginalize over errors

    return likely


def integral_element(mu, pdf):
    '''
    Returns an array of elements of the integrand dP = p(mu) dmu
    for a density p(mu) defined at sample values mu ; samples need
    not be equally spaced.  Uses a simple trapezium rule.
    Number of dP elements is 1 - (number of mu samples).
    '''
    dmu = mu[1:] - mu[:-1]
    bin_mean = (pdf[1:] + pdf[:-1]) /2
    return dmu * bin_mean


def normalize_pdf(mu, pofmu):
    """
    Takes a function pofmu defined at rate sample values mu and
    normalizes it to be a suitable pdf. Both mu and pofmu must be
    arrays or lists of the same length.
    """
    if min(pofmu) < 0:
        raise ValueError, "Probabilities cannot be negative, don't \
        ask me to normalize a function with negative values!"
    if min(mu) < 0:  # maybe not necessary ?
        raise ValueError, "Rates cannot be negative, don't \
        ask me to normalize a function over a negative domain!"

    dp = integral_element(mu, pofmu)
    return mu, pofmu/sum(dp)


def compute_upper_limit(mu_in, post, alpha = 0.9):
    """
    Returns the upper limit mu_high of confidence level alpha for a
    posterior distribution post on the given parameter mu.
    The posterior need not be normalized.
    """
    if 0 < alpha < 1:
        dp = integral_element(mu_in, post)
        high_idx = bisect.bisect_left( dp.cumsum()/dp.sum(), alpha )
        # if alpha is in (0,1] and post is non-negative, bisect_left
        # will always return an index in the range of mu since
        # post.cumsum()/post.sum() will always begin at 0 and end at 1
        mu_high = mu_in[high_idx]
    elif alpha == 1:
        mu_high = numpy.max(mu_in[post > 0])
    else:
        raise ValueError, "Confidence level must be in (0,1]."

    return mu_high


def compute_lower_limit(mu_in, post, alpha = 0.9):
    """
    Returns the lower limit mu_low of confidence level alpha for a
    posterior distribution post on the given parameter mu.
    The posterior need not be normalized.
    """
    if 0 < alpha < 1:
        dp = integral_element(mu_in, post)
        low_idx = bisect.bisect_right( dp.cumsum()/dp.sum(), 1-alpha )
        # if alpha is in [0,1) and post is non-negative, bisect_right
        # will always return an index in the range of mu since
        # post.cumsum()/post.sum() will always begin at 0 and end at 1
        mu_low = mu_in[low_idx]
    elif alpha == 1:
        mu_low = numpy.min(mu_in[post > 0])
    else:
        raise ValueError, "Confidence level must be in (0,1]."

    return mu_low


def confidence_interval_min_width( mu, post, alpha = 0.9 ):
    '''
    Returns the minimal-width confidence interval [mu_low, mu_high] of
    confidence level alpha for a posterior distribution post on the parameter mu.
    '''
    if not 0 < alpha < 1:
        raise ValueError, "Confidence level must be in (0,1)."

    # choose a step size for the sliding confidence window
    alpha_step = 0.01

    # initialize the lower and upper limits
    mu_low = numpy.min(mu)
    mu_high = numpy.max(mu)

    # find the smallest window (by delta-mu) stepping by dalpha
    for ai in numpy.arange( 0, 1-alpha, alpha_step ):
        ml = compute_lower_limit( mu, post, 1 - ai )
        mh = compute_upper_limit( mu, post, alpha + ai)
        if mh - ml < mu_high - mu_low:
            mu_low = ml
            mu_high = mh

    return mu_low, mu_high


def hpd_coverage(mu, pdf, thresh):
    '''
    Integrates a pdf over mu taking only bins where
    the mean over the bin is above a given threshold
    This gives the coverage of the HPD interval for
    the given threshold.
    '''
    dp = integral_element(mu, pdf)
    bin_mean = (pdf[1:] + pdf[:-1]) /2

    return dp[bin_mean > thresh].sum()


def hpd_threshold(mu_in, post, alpha, tol):
    '''
    For a PDF post over samples mu_in, find a density
    threshold such that the region having higher density
    has coverage of at least alpha, and less than alpha
    plus a given tolerance.
    '''
    norm_post = normalize_pdf(mu_in, post)
    # initialize bisection search
    p_minus = 0.0
    p_plus = max(post)
    while abs(hpd_coverage(mu_in, post, p_minus)-hpd_coverage(mu_in, post, p_plus)) >= tol:
        test = (p_minus + p_plus) /2
        if hpd_coverage(mu_in, post, test) >= alpha: # test value was too low or just right
            p_minus = p_test
        else:                                        # test value was too high
            p_plus = p_test
    # p_minus never goes above the required threshold and p_plus never goes below
    # thus on exiting p_minus is at or below the required threshold and the
    # difference in coverage is within tolerance

    return p_minus


def hpd_credible_interval(mu_in, post, alpha = 0.9, tolerance = 1e-3):
    '''
    Returns the minimum and maximum rate values of the HPD
    (Highest Posterior Density) credible interval for a posterior
    post defined at the sample values mu_in.  Samples need not be
    uniformly spaced and posterior need not be normalized.

    Will not return a correct credible interval if the posterior
    is multimodal and the correct interval is not contiguous;
    in this case will over-cover by including the whole range from
    minimum to maximum mu.
    '''
    if alpha == 1:
        nonzero_samples = mu_in[post > 0]
        mu_low = numpy.min(nonzero_samples)
        mu_high = numpy.max(nonzero_samples)
    elif 0 < alpha < 1:
        # determine the highest PDF for which the region with
        # higher density has sufficient coverage
        pthresh = hpd_threshold(mu_in, post, alpha, tol = tolerance)
        samples_over_threshold = mu_in[post > pthresh]
        mu_low = numpy.min(samples_over_threshold)
        mu_high = numpy.max(samples_over_threshold)

    return mu_low, mu_high


def integrate_efficiency(dbins, eff, err=0, logbins=False):
    # NB logbins is only called in the unit tests
    if logbins:
        logd = numpy.log(dbins)
        dlogd = logd[1:] - logd[:-1]
        dreps = numpy.exp( (numpy.log(dbins[1:]) + numpy.log(dbins[:-1])) / 2. ) # log midpoint
        vol = numpy.sum(4. * numpy.pi * (dreps ** 3.) * eff * dlogd)
        verr = numpy.sqrt( numpy.sum((4. * numpy.pi * (dreps ** 3.) * err * dlogd) ** 2.) ) #propagate errors in eff to errors in v
    else:
        dd = dbins[1:] - dbins[:-1]
        dreps = (dbins[1:] + dbins[:-1]) / 2. # midpoint
        vol = numpy.sum(4. * numpy.pi * (dreps ** 2.) * eff * dd )
        verr = numpy.sqrt( numpy.sum((4. * numpy.pi * (dreps ** 2.) * err * dd) ** 2.) ) #propagate errors in eff to errors in v

    return vol, verr


def compute_efficiency(f_dist, m_dist, dbins):
    '''
    Compute the efficiency as a function of distance for the given sets of found
    and missed injection distances.
    Note that injections that do not fit into any dbin get lost :(.
    '''
    efficiency = numpy.zeros(len(dbins) - 1)
    error = numpy.zeros(len(dbins) - 1)
    for j, dlow in enumerate(dbins[:-1]):
        dhigh = dbins[j + 1]
        found = numpy.sum( (dlow <= f_dist) * (f_dist < dhigh) )
        missed = numpy.sum( (dlow <= m_dist) * (m_dist < dhigh) )
        if found + missed == 0: missed = 1.0     # avoid divide by 0 in empty bins
        efficiency[j] = found / (found + missed) # NB division is imported from __future__ !
        error[j] = numpy.sqrt( efficiency[j] * (1 - efficiency[j]) / (found + missed) )

    return efficiency, error


def mean_efficiency_volume(found, missed, dbins):

    if len(found) == 0: # no efficiency here
        return numpy.zeros(len(dbins)-1), numpy.zeros(len(dbins)-1), 0, 0

    # only need distances
    f_dist = numpy.array([l.distance for l in found])
    m_dist = numpy.array([l.distance for l in missed])

    # compute the efficiency and its variance
    eff, err = compute_efficiency(f_dist, m_dist, dbins)
    vol, verr = integrate_efficiency(dbins, eff, err)

    return eff, err, vol, verr


def volume_montecarlo(found, missed, distribution_param, distribution, limits_param, max_param=None, min_param=None):
    '''
    Compute the sensitive volume and standard error using a direct Monte Carlo integral

    * distribution_param, D: parameter of the injections used to generate a distribution over distance
      - may be 'distance', 'chirp_distance"
    * distribution: form of the distribution over the parameter
      - 'log' (uniform in log D), 'uniform' (uniform in D), 'distancesquared' (uniform in D**2),
        'volume' (uniform in D***3)
      - It is assumed that injections were carried out over a range of D such that sensitive
        volume due to signals at distances < D_min is negligible and efficiency at distances
        > D_max is negligibly small
    * limits_param, Dlim: parameter specifying limits in which injections were made
      - may be 'distance', 'chirp_distance'
    * max_param: maximum value of Dlim out to which injections were made, if None the maximum 
      value among all found and missed injections will be used
    * min_param: minimum value of Dlim at which injections were made - needed to normalize
      the log distance integral correctly.  If None, for the log distribution the minimum
      value among all found and missed injections will be used
    '''
    d_power = {
        'log'             : 3.,
        'uniform'         : 2.,
        'distancesquared' : 1.,
        'volume'          : 0.
    }[distribution]
    mchirp_power = {
        'log'             : 0.,
        'uniform'         : 5. / 6.,
        'distancesquared' : 5. / 3.,
        'volume'          : 5. / 2.
    }[distribution]

    found_d = numpy.array([inj.distance for inj in found])
    missed_d = numpy.array([inj.distance for inj in missed])

    # establish maximum physical distance: first in case of chirp distance distribution
    if limits_param == 'chirp_distance':
        mchirp_standard_bns = 1.4 * (2. ** (-1. / 5.))
        found_mchirp = numpy.array([inj.mchirp for inj in found])
        missed_mchirp = numpy.array([inj.mchirp for inj in missed])
        all_mchirp = numpy.concatenate((found_mchirp, missed_mchirp))
        max_mchirp = numpy.max(all_mchirp)
        if max_param is not None:
            # use largest actually injected mchirp for conversion
            max_distance = max_param * (max_mchirp / mchirp_standard_bns) ** (5. / 6.)
    elif limits_param == 'distance':
        max_distance = max_param
    else: raise NotImplementedError("%s is not a recognized parameter" % limits_param)

    # if no max distance given, use maximum distance actually injected
    if max_param == None:
        max_distance = max(numpy.max(found_d), numpy.max(missed_d))

    # volume of sphere
    montecarlo_vtot = (4. / 3.) * numpy.pi * max_distance ** 3.

    # arrays of weights for the MC integral
    if distribution_param == 'distance':
        found_weights = found_d ** d_power
        missed_weights = missed_d ** d_power
    elif distribution_param == 'chirp_distance':
        # weight by a power of mchirp to rescale injection density to the
        # target mass distribution
        found_weights = found_d ** d_power * \
                        found_mchirp ** mchirp_power
        missed_weights = missed_d ** d_power * \
                         missed_mchirp ** mchirp_power
    else: raise NotImplementedError("%s is not a recognized distance parameter" % distance_param)

    all_weights = numpy.concatenate((found_weights, missed_weights))

    # MC integral is volume of sphere * (sum of found weights)/(sum of all weights)
    # over injections covering the sphere
    mc_weight_samples = numpy.concatenate((found_weights, 0 * missed_weights))
    mc_sum = sum(mc_weight_samples)

    if limits_param == 'distance':
        mc_norm = sum(all_weights)
    elif limits_param == 'chirp_distance':
        # if injections are made up to a maximum chirp distance, account for
        # extra missed injections that would occur when injecting up to
        # maximum physical distance : this works out to a 'chirp volume' factor
        mc_norm = sum(all_weights * (max_mchirp / all_mchirp) ** (5. / 2.))

    # take out a constant factor
    mc_prefactor = montecarlo_vtot / mc_norm

    # count the samples
    if limits_param == 'distance':
        Ninj = len(mc_weight_samples)
    elif limits_param == 'chirp_distance':
        # find the total expected number after extending from maximum chirp
        # dist up to maximum physical distance
        if distribution == 'log':
            # need minimum distance only in this case
            if min_param is not None:
                min_distance = min_param * (numpy.min(all_mchirp) / mchirp_standard_bns) ** (5. / 6.)
            else:
                min_distance = min(numpy.min(found_d), numpy.min(missed_d))
            logrange = numpy.log(max_distance / min_distance)
            Ninj = len(mc_weight_samples) + (5. / 6.) * sum(numpy.log(max_mchirp / all_mchirp) / logrange)
        else:
            Ninj = sum((max_mchirp / all_mchirp) ** mchirp_power)

    # sample variance of injection efficiency: mean of the square - square of the mean
    mc_sample_variance = sum(mc_weight_samples ** 2.) / Ninj - (mc_sum / Ninj) ** 2.

    # return MC integral and its standard deviation; variance of mc_sum scales
    # relative to sample variance by Ninj (Bienayme' rule)
    return mc_prefactor * mc_sum, mc_prefactor * (Ninj * mc_sample_variance) ** 0.5


def filter_injections_by_mass(injs, mbins, bin_num , bin_type, bin_num2=None):
    '''
    For a given set of injections (sim_inspiral rows), return the subset
    of injections that fall within the given mass range.
    '''
    if bin_type == "Mass1_Mass2":
        m1bins = numpy.concatenate((mbins.lower()[0],numpy.array([mbins.upper()[0][-1]])))
        m1lo = m1bins[bin_num]
        m1hi = m1bins[bin_num+1]
        m2bins = numpy.concatenate((mbins.lower()[1],numpy.array([mbins.upper()[1][-1]])))
        m2lo = m2bins[bin_num2]
        m2hi = m2bins[bin_num2+1]
        newinjs = [l for l in injs if ( (m1lo<= l.mass1 <m1hi and m2lo<= l.mass2 <m2hi) or (m1lo<= l.mass2 <m1hi and m2lo<= l.mass1 <m2hi))]
        return newinjs

    mbins = numpy.concatenate((mbins.lower()[0],numpy.array([mbins.upper()[0][-1]])))
    mlow = mbins[bin_num]
    mhigh = mbins[bin_num+1]
    if bin_type == "Chirp_Mass":
        newinjs = [l for l in injs if (mlow <= l.mchirp < mhigh)]
    elif bin_type == "Total_Mass":
        newinjs = [l for l in injs if (mlow <= l.mass1+l.mass2 < mhigh)]
    elif bin_type == "Component_Mass": #it is assumed that m2 is fixed
        newinjs = [l for l in injs if (mlow <= l.mass1 < mhigh)]
    elif bin_type == "BNS_BBH":
        if bin_num == 0 or bin_num == 2: #BNS/BBH case
            newinjs = [l for l in injs if (mlow <= l.mass1 < mhigh and mlow <= l.mass2 < mhigh)]
        else:
            newinjs = [l for l in injs if (mbins[0] <= l.mass1 < mbins[1] and mbins[2] <= l.mass2 < mbins[3])] #NSBH
            newinjs += [l for l in injs if (mbins[0] <= l.mass2 < mbins[1] and mbins[2] <= l.mass1 < mbins[3])] #BHNS

    return newinjs


def compute_volume_vs_mass(found, missed, mass_bins, bin_type, dbins=None,
                           distribution_param=None, distribution=None, limits_param=None,
                           max_param=None, min_param=None):
    """
    Compute the average luminosity an experiment was sensitive to given the sets
    of found and missed injections and assuming luminosity is uniformly distributed
    in space.

    If distance_param and distance_distribution are not None, use an unbinned
    Monte Carlo integral (which optionally takes max_distance and min_distance
    parameters) for the volume and error
    Otherwise use a simple efficiency * distance**2 binned integral
    In either case use distance bins to return the efficiency in each bin
    """
    # initialize MC volume integral method: binned or unbinned
    if distribution_param is not None and distribution is not None and limits_param is not None:
        mc_unbinned = True
    else:
        mc_unbinned = False

    # mean and std estimate for sensitive volume averaged over search time
    volArray = rate.BinnedArray(mass_bins)
    vol2Array = rate.BinnedArray(mass_bins)

    # found/missed stats
    foundArray = rate.BinnedArray(mass_bins)
    missedArray = rate.BinnedArray(mass_bins)

    # efficiency over distance and its standard (binomial) error
    effvmass = []
    errvmass = []

    if bin_type == "Mass1_Mass2": # two-d case first
        for j,mc1 in enumerate(mass_bins.centres()[0]):
            for k,mc2 in enumerate(mass_bins.centres()[1]):

                # filter out injections not in this mass bin
                newfound = filter_injections_by_mass(found, mass_bins, j, bin_type, k)
                newmissed = filter_injections_by_mass(missed, mass_bins, j, bin_type, k)

                foundArray[(mc1,mc2)] = len(newfound)
                missedArray[(mc1,mc2)] = len(newmissed)

                # compute the volume using this injection set
                meaneff, efferr, meanvol, volerr = mean_efficiency_volume(newfound, newmissed, dbins)
                effvmass.append(meaneff)
                errvmass.append(efferr)
                if mc_unbinned:
                    meanvol, volerr = volume_montecarlo(newfound, newmissed, distribution_param,
                                            distribution, limits_param, max_param, min_param)
                volArray[(mc1,mc2)] = meanvol
                vol2Array[(mc1,mc2)] = volerr

        return volArray, vol2Array, foundArray, missedArray, effvmass, errvmass


    for j,mc in enumerate(mass_bins.centres()[0]):

        # filter out injections not in this mass bin
        newfound = filter_injections_by_mass(found, mass_bins, j, bin_type)
        newmissed = filter_injections_by_mass(missed, mass_bins, j, bin_type)

        foundArray[(mc,)] = len(newfound)
        missedArray[(mc,)] = len(newmissed)

        # compute the volume using this injection set
        meaneff, efferr, meanvol, volerr = mean_efficiency_volume(newfound, newmissed, dbins)
        effvmass.append(meaneff)
        errvmass.append(efferr)
        # if the unbinned MC calculation is available, do it
        if mc_unbinned:
            meanvol, volerr = volume_montecarlo(newfound, newmissed, distribution_param,
                                    distribution, limits_param, max_param, min_param)
        volArray[(mc,)] = meanvol
        vol2Array[(mc,)] = volerr

    return volArray, vol2Array, foundArray, missedArray, effvmass, errvmass


def log_volume_derivative_fit(x, vols):
    '''
    Performs a linear least squares to log(vols) as a function of x.
    '''
    if numpy.min(vols) == 0:
        print >> sys.stderr, "Warning: cannot fit log volume derivative, one or more volumes are zero!"
        print >> sys.stderr, vols
        return (0,0)

    coeffs, resids, rank, svs, rcond = numpy.polyfit(x, numpy.log(vols), 1, full=True)
    if coeffs[0] < 0: #negative derivatives may arise from rounding error
        print >> sys.stderr, "Warning: Derivative fit resulted in Lambda =", coeffs[0]
        coeffs[0] = 0
        print >> sys.stderr, "The value Lambda = 0 was substituted"

    return coeffs


def get_loudest_event(connection, coinc_table="coinc_inspiral", datatype="exclude_play"):

    far_threshold_query = """
SELECT coinc_event.instruments, MIN(combined_far) FROM %s JOIN coinc_event ON (%s.coinc_event_id == coinc_event.coinc_event_id) JOIN experiment_map ON (coinc_event.coinc_event_id == experiment_map.coinc_event_id) JOIN experiment_summary ON ( experiment_summary.experiment_summ_id == experiment_map.experiment_summ_id) WHERE experiment_summary.datatype == "%s" GROUP BY coinc_event.instruments;
""" % (coinc_table, coinc_table, datatype)

    for inst, far in connection.cursor().execute(far_threshold_query):
        inst = frozenset(lsctables.instrument_set_from_ifos(inst))
        yield inst, far

    return
