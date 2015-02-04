import sys
import math
import bisect
import itertools
from collections import defaultdict

import numpy
from scipy import integrate, interpolate
import itertools
import functools

import healpy

from statutils import cumvar

class NanOrInf(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class MCSampler(object):
    """
    Class to define a set of parameter names, limits, and probability densities.
    """

    @staticmethod
    def match_params_from_args(args, params):
        """
        Given two unordered sets of parameters, one a set of all "basic" elements (strings) possible, and one a set of elements both "basic" strings and "combined" (basic strings in tuples), determine whether the sets are equivalent if no basic element is repeated.

        e.g. set A ?= set B

        ("a", "b", "c") ?= ("a", "b", "c") ==> True
        (("a", "b", "c")) ?= ("a", "b", "c") ==> True
        (("a", "b"), "d")) ?= ("a", "b", "c") ==> False  # basic element 'd' not in set B
        (("a", "b"), "d")) ?= ("a", "b", "d", "c") ==> False  # not all elements in set B represented in set A
        """
        not_common = set(args) ^ set(params)
        if len(not_common) == 0:
            # All params match
            return True
        if all([not isinstance(i, tuple) for i in not_common]):
            # The only way this is possible is if there are
            # no extraneous params in args
            return False

        to_match, against = filter(lambda i: not isinstance(i, tuple), not_common), filter(lambda i: isinstance(i, tuple), not_common)

        matched = []
        for i in range(2, max(map(len, against))+1):
            matched.extend([t for t in itertools.permutations(to_match, i) if t in against])
        return (set(matched) ^ set(against)) == set()


    def __init__(self):
        # Total number of samples drawn
        self.ntotal = 0
        # Parameter names
        self.params = set()
        # parameter -> pdf function object
        self.pdf = {}
        # If the pdfs aren't normalized, this will hold the normalization 
        # constant
        self._pdf_norm = defaultdict(lambda: 1)
        # Cache for the sampling points
        self._rvs = {}
        # parameter -> cdf^{-1} function object
        self.cdf = {}
        self.cdf_inv = {}
        # params for left and right limits
        self.llim, self.rlim = {}, {}
        # Keep track of the adaptive parameters
        self.adaptive = []

        # Keep track of the adaptive parameter 1-D marginalizations
        self._hist = {}

        # MEASURES (=priors): ROS needs these at the sampler level, to clearly separate their effects
        # ASSUMES the user insures they are normalized
        self.prior_pdf = {}

    def clear(self):
        """
        Clear out the parameters and their settings, as well as clear the sample cache.
        """
        self.params = set()
        self.pdf = {}
        self._pdf_norm = defaultdict(lambda: 1.0)
        self._rvs = {}
        self._hist = {}
        self.cdf = {}
        self.cdf_inv = {}
        self.llim = {}
        self.rlim = {}
        self.adaptive = []

    def add_parameter(self, params, pdf,  cdf_inv=None, left_limit=None, right_limit=None, prior_pdf=None, adaptive_sampling=False):
        """
        Add one (or more) parameters to sample dimensions. params is either a string describing the parameter, or a tuple of strings. The tuple will indicate to the sampler that these parameters must be sampled together. left_limit and right_limit are on the infinite interval by default, but can and probably should be specified. If several params are given, left_limit, and right_limit must be a set of tuples with corresponding length. Sampling PDF is required, and if not provided, the cdf inverse function will be determined numerically from the sampling PDF.
        """
        self.params.add(params)
        if isinstance(params, tuple):
            assert all(map(lambda lim: lim[0] < lim[1], zip(left_limit, right_limit)))
            if left_limit is None:
                self.llim[params] = list(float("-inf"))*len(params)
            else:
                self.llim[params] = left_limit
            if right_limit is None:
                self.rlim[params] = list(float("+inf"))*len(params)
            else:
                self.rlim[params] = right_limit
        else:
            assert left_limit < right_limit
            if left_limit is None:
                self.llim[params] = float("-inf")
            else:
                self.llim[params] = left_limit
            if right_limit is None:
                self.rlim[params] = float("+inf")
            else:
                self.rlim[params] = right_limit
        self.pdf[params] = pdf
        # FIXME: This only works automagically for the 1d case currently
        self.cdf_inv[params] = cdf_inv or self.cdf_inverse(params)
        if not isinstance(params, tuple):
            self.cdf[params] =  self.cdf_function(params)
        if prior_pdf is None:
            self.prior_pdf[params] = lambda x: 1.0
        else:
            self.prior_pdf[params] = prior_pdf

        if adaptive_sampling:
            self.adaptive.append(params)


    def cdf_function(self, param):
        """
        Numerically determine the  CDF from a given sampling PDF. If the PDF itself is not normalized, the class will keep an internal record of the normalization and adjust the PDF values as necessary. Returns a function object which is the interpolated CDF.
        """
        # Solve P'(x) == p(x), with P[lower_boun] == 0
        def dP_cdf(p, x):
            if x > self.rlim[param] or x < self.llim[param]:
                return 0
            return self.pdf[param](x)
        x_i = numpy.linspace(self.llim[param], self.rlim[param], 1000)
        # Integrator needs to have a step size which doesn't step over the
        # probability mass
        # TODO: Determine h_max.
        cdf = integrate.odeint(dP_cdf, [0], x_i, hmax=0.01*(self.rlim[param]-self.llim[param])).T[0]
        if cdf[-1] != 1.0: # original pdf wasn't normalized
            self._pdf_norm[param] = cdf[-1]
            cdf /= cdf[-1]
        # Interpolate the inverse
        return interpolate.interp1d( x_i,cdf)

    def cdf_inverse(self, param):
        """
        Numerically determine the inverse CDF from a given sampling PDF. If the PDF itself is not normalized, the class will keep an internal record of the normalization and adjust the PDF values as necessary. Returns a function object which is the interpolated CDF inverse.
        """
        # Solve P'(x) == p(x), with P[lower_boun] == 0
        def dP_cdf(p, x):
            if x > self.rlim[param] or x < self.llim[param]:
                return 0
            return self.pdf[param](x)
        x_i = numpy.linspace(self.llim[param], self.rlim[param], 1000)
        # Integrator needs to have a step size which doesn't step over the
        # probability mass
        # TODO: Determine h_max.
        cdf = integrate.odeint(dP_cdf, [0], x_i, hmax=0.01*(self.rlim[param]-self.llim[param])).T[0]
        if cdf[-1] != 1.0: # original pdf wasn't normalized
            self._pdf_norm[param] = cdf[-1]
            cdf /= cdf[-1]
        # Interpolate the inverse
        return interpolate.interp1d(cdf, x_i)

    def draw(self, rvs, *args, **kwargs):
        """
        Draw a set of random variates for parameter(s) args. Left and right limits are handed to the function. If args is None, then draw *all* parameters. 'rdict' parameter is a boolean. If true, returns a dict matched to param name rather than list. rvs must be either a list of uniform random variates to transform for sampling, or an integer number of samples to draw.
        """
        if len(args) == 0:
            args = self.params

        if isinstance(rvs, int) or isinstance(rvs, float):
            #
            # Convert all arguments to tuples
            #
            # FIXME: UGH! Really? This was the most elegant thing you could come
            # up with?
            rvs_tmp = [numpy.random.uniform(0,1,(len(p), rvs)) for p in map(lambda i: (i,) if not isinstance(i, tuple) else i, args)]
            rvs_tmp = numpy.array([self.cdf_inv[param](*rv) for (rv, param) in zip(rvs_tmp, args)], dtype=numpy.object)
        else:
            rvs_tmp = numpy.array(rvs)

        # FIXME: ELegance; get some of that...
        # This is mainly to ensure that the array can be "splatted", e.g.
        # separated out into its components for matching with args. The case of
        # one argument has to be handled specially.
        res = []
        for (cdf_rv, param) in zip(rvs_tmp, args):
            if len(cdf_rv.shape) == 1:
                res.append((self.pdf[param](numpy.float64(cdf_rv)).astype(numpy.float64)/self._pdf_norm[param], self.prior_pdf[param](cdf_rv), cdf_rv.astype(numpy.float64)))
            else:
                # NOTE: the "astype" is employed here because the arrays can be
                # irregular and thus assigned the 'object' type. Since object
                # arrays can't be splatted, we have to force the conversion
                res.append((self.pdf[param](*cdf_rv.astype(numpy.float64))/self._pdf_norm[param], self.prior_pdf[param](*cdf_rv.astype(numpy.float64)), cdf_rv.astype(numpy.float64)))

        #
        # Cache the samples we chose
        #
        if len(self._rvs) == 0:
            self._rvs = dict(zip(args, rvs_tmp.astype(numpy.float64)))
        else:
            rvs_tmp = dict(zip(args, rvs_tmp.astype(numpy.float64)))
            #for p, ar in self._rvs.iteritems():
            for p in self.params:
                self._rvs[p] = numpy.hstack( (self._rvs[p], rvs_tmp[p]) ).astype(numpy.float64)

        #
        # Pack up the result if the user wants a dictonary instead
        #
        if kwargs.has_key("rdict"):
            return dict(zip(args, res))
        return zip(*res)

    #
    # FIXME: The priors are not strictly part of the MC integral, and so any
    # internal reference to them needs to be moved to a subclass which handles
    # the incovnenient part os doing the \int p/p_s L d\theta integral.
    #
    def integrate(self, func, *args, **kwargs):
        """
        Integrate func, by using n sample points. Right now, all params defined must be passed to args must be provided, but this will change soon.

        Limitations:
            func's signature must contain all parameters currently defined by the sampler, and with the same names. This is required so that the sample values can be passed consistently.

        kwargs:
        nmax -- total allowed number of sample points, will throw a warning if this number is reached before neff.
        neff -- Effective samples to collect before terminating. If not given, assume infinity
        n -- Number of samples to integrate in a 'chunk' -- default is 1000
        save_integrand -- Save the evaluated value of the integrand at the sample points with the sample point
        history_mult -- Number of chunks (of size n) to use in the adaptive histogramming: only useful if there are parameters with adaptation enabled
        tempering_exp -- Exponent to raise the weights of the 1-D marginalized histograms for adaptive sampling prior generation, by default it is 0 which will turn off adaptive sampling regardless of other settings
        floor_level -- *total probability* of a uniform distribution, averaged with the weighted sampled distribution, to generate a new sampled distribution
        n_adapt -- number of chunks over which to allow the pdf to adapt. Default is zero, which will turn off adaptive sampling regardless of other settings
        convergence_tests - dictionary of function pointers, each accepting self._rvs and self.params as arguments. CURRENTLY ONLY USED FOR REPORTING

        Pinning a value: By specifying a kwarg with the same of an existing parameter, it is possible to "pin" it. The sample draws will always be that value, and the sampling prior will use a delta function at that value.
        """

        #
        # Pin values
        #
        tempcdfdict, temppdfdict, temppriordict, temppdfnormdict = {}, {}, {}, {}
        temppdfnormdict = defaultdict(lambda: 1.0)
        for p, val in kwargs.iteritems():
            if p in self.params:
                # Store the previous pdf/cdf in case it's already defined
                tempcdfdict[p] = self.cdf_inv[p]
                temppdfdict[p] = self.pdf[p]
                temppdfnormdict[p] = self._pdf_norm[p]
                temppriordict[p] = self.prior_pdf[p]
                # Set a new one to always return the same value
                self.pdf[p] = functools.partial(delta_func_pdf_vector, val)
                self._pdf_norm[p] = 1.0
                self.prior_pdf[p] = functools.partial(delta_func_pdf_vector, val)
                self.cdf_inv[p] = functools.partial(delta_func_samp_vector, val)

        # This is a semi-hack to ensure that the integrand is called with
        # the arguments in the right order
        # FIXME: How dangerous is this?
        args = func.func_code.co_varnames[:func.func_code.co_argcount]
        if not MCSampler.match_params_from_args(args, self.params):
            raise ValueError("All integrand variables must be represented by integral parameters.")
        
        #
        # Determine stopping conditions
        #
        nmax = kwargs["nmax"] if kwargs.has_key("nmax") else float("inf")
        neff = kwargs["neff"] if kwargs.has_key("neff") else numpy.float128("inf")
        n = kwargs["n"] if kwargs.has_key("n") else min(1000, nmax)
        convergence_tests = kwargs["convergence_tests"] if kwargs.has_key("convergence_tests") else None

        #
        # Adaptive sampling parameters
        #
        n_history = int(kwargs["history_mult"]*n) if kwargs.has_key("history_mult") else None
        tempering_exp = kwargs["tempering_exp"] if kwargs.has_key("tempering_exp") else 0.0
        n_adapt = int(kwargs["n_adapt"]*n) if kwargs.has_key("n_adapt") else 0
        floor_integrated_probability = kwargs["floor_level"] if kwargs.has_key("floor_level") else 0

        save_intg = kwargs["save_intg"] if kwargs.has_key("save_intg") else False
        # FIXME: The adaptive step relies on the _rvs cache, so this has to be
        # on in order to work
        if n_adapt > 0 and tempering_exp > 0.0:
            save_intg = True

        deltalnL = kwargs['igrand_threshold_deltalnL'] if kwargs.has_key('igrand_threshold_deltalnL') else float("Inf") # default is to return all
        deltaP = kwargs["igrand_threshold_p"] if kwargs.has_key('igrand_threshold_p') else 0 # default is to omit 1e-7 of probability

        show_evaluation_log = kwargs['verbose'] if kwargs.has_key('verbose') else False
        if show_evaluation_log:
            print " .... mcsampler : providing verbose output ..... "

        int_val1 = numpy.float128(0)
        self.ntotal = 0
        maxval = -float("Inf")
        maxlnL = -float("Inf")
        eff_samp = 0
        mean, var = None, None
        last_convergence_test = defaultdict(lambda: False)   # initialize record of tests

        if show_evaluation_log:
            print "iteration Neff  sqrt(2*lnLmax) sqrt(2*lnLmarg) ln(Z/Lmax) int_var"

        while (eff_samp < neff and self.ntotal < nmax): #  and (not bConvergenceTests):
            # Draw our sample points
            p_s, p_prior, rv = self.draw(n, *self.params)
                        
            # Calculate the overall p_s assuming each pdf is independent
            joint_p_s = numpy.prod(p_s, axis=0)
            joint_p_prior = numpy.prod(p_prior, axis=0)

            #
            # Prevent zeroes in the sampling prior
            #
            # FIXME: If we get too many of these, we should bail
            if (isinstance(joint_p_s, numpy.ndarray) and any(joint_p_s <= 0)) \
              or (not isinstance(joint_p_s, numpy.ndarray) and joint_p_s <= 0):
                for p in self.params:
                    self._rvs[p] = numpy.resize(self._rvs[p], len(self._rvs[p])-n)
                print >>sys.stderr, "Zero prior value detected, skipping."
                continue

            #
            # Unpack rvs and evaluate integrand
            #
            if len(rv[0].shape) != 1:
                rv = rv[0]

            params = []
            for item in self.params:
                if isinstance(item, tuple):
                    params.extend(item)
                else:
                    params.append(item)
            unpacked = numpy.hstack([r.flatten() for r in rv]).reshape(len(args), -1)
            unpacked = dict(zip(params, unpacked))
            fval = func(**unpacked)

            #
            # Check if there is any practical contribution to the integral
            #
            # FIXME: While not technically a fatal error, this will kill the 
            # adaptive sampling
            if fval.sum() == 0:
                for p in self.params:
                    self._rvs[p] = numpy.resize(self._rvs[p], len(self._rvs[p])-n)
                print >>sys.stderr, "No contribution to integral, skipping."
                continue

            if save_intg:
                # FIXME: The joint_prior, if not specified is set to one and
                # will come out as a scalar here, hence the hack
                if not isinstance(joint_p_prior, numpy.ndarray):
                    joint_p_prior = numpy.ones(fval.shape)*joint_p_prior

                # FIXME: See warning at beginning of function. The prior values
                # need to be moved out of this, as they are not part of MC
                # integration
                if self._rvs.has_key("integrand"):
                    self._rvs["integrand"] = numpy.hstack( (self._rvs["integrand"], fval) )
                    self._rvs["joint_prior"] = numpy.hstack( (self._rvs["joint_prior"], joint_p_prior) )
                    self._rvs["joint_s_prior"] = numpy.hstack( (self._rvs["joint_s_prior"], joint_p_s) )
                    self._rvs["weights"] = numpy.hstack( (self._rvs["joint_s_prior"], fval*joint_p_prior/joint_p_s) )
                else:
                    self._rvs["integrand"] = fval
                    self._rvs["joint_prior"] = joint_p_prior
                    self._rvs["joint_s_prior"] = joint_p_s
                    self._rvs["weights"] = fval*joint_p_prior/joint_p_s

            # Calculate the integral over this chunk
            int_val = fval * joint_p_prior / joint_p_s

            # Calculate max L (a useful convergence feature) for debug 
            # reporting.  Not used for integration
            # Try to avoid nan's
            maxlnL = numpy.log(numpy.max([numpy.exp(maxlnL), numpy.max(fval),numpy.exp(-100)]))   # note if f<0, this will return nearly 0

            # Calculate the effective samples via max over the current 
            # evaluations
            maxval = [max(maxval, int_val[0]) if int_val[0] != 0 else maxval]
            for v in int_val[1:]:
                maxval.append( v if v > maxval[-1] and v != 0 else maxval[-1] )

            # running variance
            var = cumvar(int_val, mean, var, self.ntotal)[-1]
            # running integral
            int_val1 += int_val.sum()
            # running number of evaluations
            self.ntotal += n
            # FIXME: Likely redundant with int_val1
            mean = int_val1/self.ntotal
            maxval = maxval[-1]

            eff_samp = int_val1/maxval

            # Throw exception if we get infinity or nan
            if math.isnan(eff_samp):
                raise NanOrInf("Effective samples = nan")
            if maxlnL is float("Inf"):
                raise NanOrInf("maxlnL = inf")

            if show_evaluation_log:
                print " :",  self.ntotal, eff_samp, numpy.sqrt(2*maxlnL), numpy.sqrt(2*numpy.log(int_val1/self.ntotal)), numpy.log(int_val1/self.ntotal)-maxlnL, numpy.sqrt(var*self.ntotal)/int_val1

            if (not convergence_tests) and self.ntotal >= nmax and neff != float("inf"):
                print >>sys.stderr, "WARNING: User requested maximum number of samples reached... bailing."

            #
            # Convergence tests
            #
            if convergence_tests:
                converged = True
                for key in convergence_tests.keys():
                    last_convergence_test[key] = convergence_tests[key](self._rvs, self.params)
                    converged &= las_convergence_test[key]

            if convergence_tests and show_evaluation_log:  # Print status of each test
                for key in convergence_tests:
                    print "   -- Convergence test status : ", key, last_convergence_test[key]

            #
            # The total number of adaptive steps is reached
            #
            # FIXME: We need a better stopping condition here
            if self.ntotal > n_adapt:
                continue

            #
            # Iterate through each of the parameters, updating the sampling
            # prior PDF according to the 1-D marginalization
            #
            for itr, p in enumerate(self.params):
                # FIXME: The second part of this condition should be made more
                # specific to pinned parameters
                if p not in self.adaptive or p in kwargs.keys():
                    continue
                points = self._rvs[p][-n_history:]
                weights = (self._rvs["integrand"][-n_history:]/self._rvs["joint_s_prior"][-n_history:]*self._rvs["joint_prior"][-n_history:])**tempering_exp

                self._hist[p], edges = numpy.histogram(points, bins = 100, range = (self.llim[p], self.rlim[p]), weights = weights, normed = True)
                # FIXME: numpy.hist can't normalize worth a damn
                self._hist[p] /= self._hist[p].sum()

                # Mix with uniform distribution
                self._hist[p] = (1-floor_integrated_probability)*self._hist[p] + numpy.ones(len(self._hist[p]))*floor_integrated_probability/len(self._hist[p])

                # FIXME: These should be called the bin centers
                edges = [(e0+e1)/2.0 for e0, e1 in zip(edges[:-1], edges[1:])]
                edges.append(edges[-1] + (edges[-1] - edges[-2]))
                edges.insert(0, edges[0] - (edges[-1] - edges[-2]))

                # FIXME: KS test probably has a place here. Not sure yet where.
                #from scipy.stats import kstest
                #d, pval = kstest(self._rvs[p][-n_history:], self.cdf[p])
                #print p, d, pval

                self.pdf[p] = interpolate.interp1d(edges, [0] + list(self._hist[p]) + [0])
                self.cdf[p] = self.cdf_function(p)
                self.cdf_inv[p] = self.cdf_inverse(p)

        # If we were pinning any values, undo the changes we did before
        self.cdf_inv.update(tempcdfdict)
        self.pdf.update(temppdfdict)
        self._pdf_norm.update(temppdfnormdict)
        self.prior_pdf.update(temppriordict)

        # Clean out the _rvs arrays for 'irrelevant' points
        #   - find and remove samples with  lnL less than maxlnL - deltalnL (latter user-specified)
        #   - create the cumulative weights
        #   - find and remove samples which contribute too little to the cumulative weights
        # FIXME: This needs *a lot* of work
        if "integrand" in self._rvs:
            self._rvs["sample_n"] = numpy.arange(len(self._rvs["integrand"]))  # create 'iteration number'        
            # Step 1: Cut out any sample with lnL belw threshold
            indx_list = [k for k, value in enumerate( (self._rvs["integrand"] > maxlnL - deltalnL)) if value] # threshold number 1
            # FIXME: This is an unncessary initial copy, the second step (cum i
            # prob) can be accomplished with indexing first then only pare at
            # the end
            for key in self._rvs.keys():
                if isinstance(key, tuple):
                    self._rvs[key] = self._rvs[key][:,indx_list]
                else:
                    self._rvs[key] = self._rvs[key][indx_list]
            # Step 2: Create and sort the cumulative weights, among the remaining points, then use that as a threshold
            wt = self._rvs["integrand"]*self._rvs["joint_prior"]/self._rvs["joint_s_prior"]
            idx_sorted_index = numpy.lexsort((numpy.arange(len(wt)), wt))  # Sort the array of weights, recovering index values
            indx_list = numpy.array( [[k, wt[k]] for k in idx_sorted_index])     # pair up with the weights again
            cum_sum = numpy.cumsum(indx_list[:,1])  # find the cumulative sum
            cum_sum = cum_sum/cum_sum[-1]          # normalize the cumulative sum
            indx_list = [indx_list[k, 0] for k, value in enumerate(cum_sum > deltaP) if value]  # find the indices that preserve > 1e-7 of total probability
            # FIXME: See previous FIXME
            for key in self._rvs.keys():
                if isinstance(key, tuple):
                    self._rvs[key] = self._rvs[key][:,indx_list]
                else:
                    self._rvs[key] = self._rvs[key][indx_list]

        # Create extra dictionary to return things
        dict_return ={}
        if convergence_tests is not None:
            dict_return["convergence_test_results"] = last_convergence_test

        return int_val1/self.ntotal, var/self.ntotal, eff_samp, dict_return

class HealPixSampler(object):
    """
    Class to sample the sky using a FITS healpix map. Equivalent to a joint 2-D pdf in RA and dec.
    """

    @staticmethod
    def thph2decra(th, ph):
        """
        theta/phi to RA/dec
        theta (north to south) (0, pi)
        phi (east to west) (0, 2*pi)
        declination: north pole = pi/2, south pole = -pi/2
        right ascension: (0, 2*pi)
        
        dec = pi/2 - theta
        ra = phi
        """
        return numpy.pi/2-th, ph

    @staticmethod
    def decra2thph(dec, ra):
        """
        theta/phi to RA/dec
        theta (north to south) (0, pi)
        phi (east to west) (0, 2*pi)
        declination: north pole = pi/2, south pole = -pi/2
        right ascension: (0, 2*pi)
        
        theta = pi/2 - dec
        ra = phi
        """
        return numpy.pi/2-dec, ra

    def __init__(self, skymap, massp=1.0):
        self.skymap = skymap
        self._massp = massp
        self.renormalize()

    @property
    def massp(self):
        return self._massp

    @massp.setter
    def massp(self, value):
        assert 0 <= value <= 1
        self._massp = value
        norm = self.renormalize()

    def renormalize(self):
        """
        Identify the points contributing to the overall cumulative probability distribution, and set the proper normalization.
        """
        res = healpy.npix2nside(len(self.skymap))
        self.pdf_sorted = sorted([(p, i) for i, p in enumerate(self.skymap)], reverse=True)
        self.valid_points_decra = []
        cdf, np = 0, 0
        for p, i in self.pdf_sorted:
            if p == 0:
                continue # Can't have a zero prior
            self.valid_points_decra.append(HealPixSampler.thph2decra(*healpy.pix2ang(res, i)))
            cdf += p
            if cdf > self._massp:
                break
        self._renorm = cdf
        # reset to indicate we'd need to recalculate this
        self.valid_points_hist = None
        return self._renorm

    def __expand_valid(self, min_p=1e-7):
        #
        # Determine what the 'quanta' of probabilty is
        #
        if self._massp == 1.0:
            # This is to ensure we don't blow away everything because the map
            # is very spread out
            min_p = min(min_p, max(self.skymap))
        else:
            # NOTE: Only valid if CDF descending order is kept
            min_p = self.pseudo_pdf(*self.valid_points_decra[-1])

        self.valid_points_hist = []
        ns = healpy.npix2nside(len(self.skymap))

        # Renormalize first so that the vector histogram is properly normalized
        self._renorm = 0
        # Account for probability lost due to cut off
        for i, v in enumerate(self.skymap >= min_p):
            self._renorm += self.skymap[i] if v else 0

        for pt in self.valid_points_decra:
            th, ph = HealPixSampler.decra2thph(pt[0], pt[1])
            pix = healpy.ang2pix(ns, th, ph)
            if self.skymap[pix] < min_p:
                continue
            self.valid_points_hist.extend([pt]*int(round(self.pseudo_pdf(*pt)/min_p)))
        self.valid_points_hist = numpy.array(self.valid_points_hist).T

    def pseudo_pdf(self, dec_in, ra_in):
        """
        Return pixel probability for a given dec_in and ra_in. Note, uses healpy functions to identify correct pixel.
        """
        th, ph = HealPixSampler.decra2thph(dec_in, ra_in)
        res = healpy.npix2nside(len(self.skymap))
        return self.skymap[healpy.ang2pix(res, th, ph)]/self._renorm

    def pseudo_cdf_inverse(self, dec_in=None, ra_in=None, ndraws=1, stype='vecthist'):
        """
        Select points from the skymap with a distribution following its corresponding pixel probability. If dec_in, ra_in are suupplied, they are ignored except that their shape is reproduced. If ndraws is supplied, that will set the shape. Will return a 2xN numpy array of the (dec, ra) values.
        stype controls the type of sampling done to retrieve points. Valid choices are
        'rejsamp': Rejection sampling: accurate but slow
        'vecthist': Expands a set of points into a larger vector with the multiplicity of the points in the vector corresponding roughly to the probability of drawing that point. Because this is not an exact representation of the proability, some points may not be represented at all (less than quantum of minimum probability) or inaccurately (a significant fraction of the fundamental quantum).
        """
        # FIXME: Add the multinomial selection here

        if ra_in is not None:
            ndraws = len(ra_in)
        if ra_in is None:
            ra_in, dec_in = numpy.zeros((2, ndraws))

        if stype == 'rejsamp':
            # FIXME: This is only valid under descending ordered CDF summation
            ceiling = max(self.skymap)
            i, np = 0, len(self.valid_points_decra)
            while i < len(ra_in):
                rnd_n = numpy.random.randint(0, np)
                trial = numpy.random.uniform(0, ceiling)
                if trial <= self.pseudo_pdf(*self.valid_points_decra[rnd_n]):
                    dec_in[i], ra_in[i] = self.valid_points_decra[rnd_n]
                    i += 1
            return numpy.array([dec_in, ra_in])
        elif stype == 'vecthist':
            if self.valid_points_hist is None:
                self.__expand_valid()
            np = self.valid_points_hist.shape[1]
            rnd_n = numpy.random.randint(0, np, len(ra_in))
            dec_in, ra_in = self.valid_points_hist[:,rnd_n]
            return numpy.array([dec_in, ra_in])
        else:
            raise ValueError("%s is not a recgonized sampling type" % stype)

### UTILITIES: Predefined distributions

def uniform_samp(a, b, x):
    if x > a and x < b:
        return 1.0/(b-a)
    else:
        return 0
uniform_samp_vector = numpy.vectorize(uniform_samp,otypes=[numpy.float])

# syntatic sugar : predefine the most common distributions
uniform_samp_phase = numpy.vectorize(lambda x: 1/(2*numpy.pi))
uniform_samp_psi = numpy.vectorize(lambda x: 1/(numpy.pi))
uniform_samp_theta = numpy.vectorize(lambda x: numpy.sin(x)/(2))
uniform_samp_dec = numpy.vectorize(lambda x: numpy.cos(x)/(2))

def quadratic_samp(rmax, x):
    if x < rmax:
        return x**2/(3*rmax**3)
    else:
        return 0
quadratic_samp_vector = numpy.vectorize(quadratic_samp, otypes=[numpy.float])

def inv_uniform_cdf(a, b, x):
    return (b-a)*x+a

def gauss_samp(mu, std, x):
    return 1.0/numpy.sqrt(2*numpy.pi*std**2)*numpy.exp(-(x-mu)**2/2/std**2)

def cos_samp(x):
    return numpy.sin(x)/2   # x from 0, pi
cos_samp_vector = numpy.vectorize(cos_samp,otypes=[numpy.float])

def dec_samp(x):
    return numpy.sin(x+numpy.pi/2)/2   # x from 0, pi
dec_samp_vector = numpy.vectorize(dec_samp,otypes=[numpy.float])

def pseudo_dist_samp(r0,r):
    return r*r*numpy.exp( - (r0/r)*(r0/r)/2. + r0/r)+0.01  # put a floor on probability, so we converge. Note this floor only cuts out NEARBY distances
pseudo_dist_samp_vector = numpy.vectorize(pseudo_dist_samp,otypes=[numpy.float])

def delta_func_pdf(x_0, x):
    return 1.0 if x == x_0 else 0.0
delta_func_pdf_vector = numpy.vectorize(delta_func_pdf, otypes=[numpy.float])

def delta_func_samp(x_0, x):
    return x_0
delta_func_samp_vector = numpy.vectorize(delta_func_samp, otypes=[numpy.float])
