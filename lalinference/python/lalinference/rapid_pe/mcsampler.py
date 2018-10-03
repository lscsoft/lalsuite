# Copyright (C) 2013, 2015  Chris Pankow, Richard O'Shaughnessy
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

"""
Collection of classes and numpy functions to facilitate a random sampler Monte-Carlo integrator.
"""

from __future__ import print_function

import sys
import math
import bisect
import itertools
import time
from collections import defaultdict

import numpy
from scipy import integrate, interpolate
import itertools
import functools

import healpy

from . import statutils, synchlib
from lal import rate

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

        # Used if networking is useful
        self._address, self._port = None, None
        self.socket = None

    def init_socket(self, address, port):
        """
        Initialize a network connection to 'address' at 'port'.
        """
        self._address, self._port = address, port
        try:
            self.socket = synchlib.init_socket(self._address, self._port)
        except pysocket.error as sockerr:
            pass

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

    def add_parameter(self, params, pdf, cdf_inv=None, left_limit=None, right_limit=None, prior_pdf=None, adaptive_sampling=False):
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
            self._rvs = dict(zip(args, rvs_tmp))
        else:
            rvs_tmp = dict(zip(args, rvs_tmp))
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
        nmin -- minimum number of samples to allow, by default will be set to the end of the 'burnin' at n_adapt * n
        neff -- Effective samples to collect before terminating. If not given, assume infinity
        n -- Number of samples to integrate in a 'chunk' -- default is 1000
        save_integrand -- Save the evaluated value of the integrand at the sample points with the sample point
        history_mult -- Number of chunks (of size n) to use in the adaptive histogramming: only useful if there are parameters with adaptation enabled
        tempering_exp -- Exponent to raise the weights of the 1-D marginalized histograms for adaptive sampling prior generation, by default it is 0 which will turn off adaptive sampling regardless of other settings
        n_adapt -- number of chunks over which to allow the pdf to adapt. Default is zero, which will turn off adaptive sampling regardless of other settings
        convergence_tests - dictionary of function pointers, each accepting self._rvs and self.params as arguments. CURRENTLY ONLY USED FOR REPORTING
        maxval - Guess at the maximum value of the integrand -- used as a seed for the maxval counter

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
        tempering_exp = kwargs["tempering_exp"] if kwargs.has_key("tempering_exp") else 1.0
        n_adapt = int(kwargs["n_adapt"]*n) if kwargs.has_key("n_adapt") else 0
        nmax = kwargs["nmax"] if kwargs.has_key("nmax") else n_adapt
        nmin = kwargs["nmin"] if kwargs.has_key("nmin") else n_adapt

        save_intg = kwargs["save_intg"] if kwargs.has_key("save_intg") else False
        nkeep = kwargs["save_intg"] if kwargs.has_key("save_intg") else None
        # Corner case: we want all the samples, and we don't want them messed
        # with, so everything is saved, but no sort is done
        if nkeep is True:
            nkeep = None

        # FIXME: The adaptive step relies on the _rvs cache, so this has to be
        # on in order to work
        if n_adapt > 0 and tempering_exp > 0.0:
            save_intg = True

        deltalnL = kwargs['igrand_threshold_deltalnL'] if kwargs.has_key('igrand_threshold_deltalnL') else None # default is to return all
        deltaP = kwargs["igrand_threshold_p"] if kwargs.has_key('igrand_threshold_p') else 0 # default is to omit 1e-7 of probability

        show_evaluation_log = kwargs['verbose'] if kwargs.has_key('verbose') else False
        if show_evaluation_log:
            print(" .... mcsampler : providing verbose output ..... ")

        int_val1 = numpy.float128(0)
        self.ntotal = 0
        maxval = kwargs["maxval"] if "maxval" in kwargs else -float("Inf")
        old_maxval = maxval
        maxlnL = -float("Inf")
        eff_samp = 0
        mean, var = None, None
        last_convergence_test = defaultdict(lambda: False)   # initialize record of tests

        if show_evaluation_log:
            print("walltime : iteration Neff  ln(maxweight) lnLmarg ln(Z/Lmax) int_var")

        socket = None
        while self.ntotal < nmin or (eff_samp < neff and self.ntotal < nmax):
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
                print("Zero prior value detected, skipping.", file=sys.stderr)
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
                print("No contribution to integral, skipping.", file=sys.stderr)
                continue

            sample_n = numpy.arange(self.ntotal, self.ntotal + len(fval))

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
                    self._rvs["weights"] = numpy.hstack( (self._rvs["weights"], fval*joint_p_prior/joint_p_s) )
                    self._rvs["sample_n"] = numpy.hstack( (self._rvs["sample_n"], sample_n) )
                else:
                    self._rvs["integrand"] = fval
                    self._rvs["joint_prior"] = joint_p_prior
                    self._rvs["joint_s_prior"] = joint_p_s
                    self._rvs["weights"] = fval*joint_p_prior/joint_p_s
                    self._rvs["sample_n"] = sample_n

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
            var = statutils.cumvar(int_val, mean, var, self.ntotal)[-1]
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
                print("{0:.3f} : {1:d} {2:.5f} {3:.2f} {4:.2f} {5:.2f} {6:.3f}".format(time.time(), self.ntotal, eff_samp, math.log(maxval), numpy.log(int_val1 / self.ntotal), numpy.log(int_val1 / self.ntotal) - maxlnL, numpy.sqrt(var * self.ntotal) / int_val1))

            if (not convergence_tests) and self.ntotal >= nmin and self.ntotal >= nmax and neff != float("inf"):
                print("WARNING: User requested maximum number of samples reached... bailing.", file=sys.stderr)

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
                    print("   -- Convergence test status : ", key, last_convergence_test[key])

            self._address, self._port = "pcdev2.nemo.phys.uwm.edu", 1890
            #if self._address is not None:
            if False:
                dims = ("distance", "inclination", "right_ascension",
                        "declination", "integrand", "joint_prior", "joint_s_prior")
                send_data = synchlib.prepare_data(self._rvs, dims, self.ntotal - n)
                self.socket = synchlib.send_samples(send_data, self._address, self._port, verbose=True, socket=self.socket)

            #
            # The total number of adaptive steps is reached
            #
            # FIXME: We need a better stopping condition here
            if self.ntotal >= n_adapt and maxval == old_maxval:
                # Downsample points
                if save_intg and nkeep is not None:
                    pt_sort = self._rvs["weights"].argsort()[-nkeep:]
                    for key in self._rvs:
                        if len(self._rvs[key].shape) > 1:
                            self._rvs[key] = self._rvs[key][:,pt_sort]
                        else:
                            self._rvs[key] = self._rvs[key][pt_sort]
                continue
            old_maxval = maxval

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
                weights = (self._rvs["integrand"][-n_history:]*self._rvs["joint_prior"][-n_history:])**tempering_exp

                self._hist[p] = statutils.get_adaptive_binning(points, (self.llim[p], self.rlim[p]))
                for pt, w in zip(points, weights):
                    self._hist[p][pt,] += w
                self._hist[p].array += self._hist[p].array.mean()
                rate.filter_array(self._hist[p].array, rate.tophat_window(3))
                norm = numpy.sum(self._hist[p].array * self._hist[p].bins.volumes())
                self._hist[p].array /= norm
                # FIXME: Stupid pet trick while numpy version is lacking
                self.pdf[p] = numpy.frompyfunc(rate.InterpBinnedArray(self._hist[p]), 1, 1)
                #with open("%s_%d_hist.txt" % (p, self.ntotal), "w") as fout:
                    #for c, pdf in zip(self._hist[p].centres()[0], self._hist[p].array):
                        #print >>fout, "%f %g" % (c, pdf)

                self.cdf[p] = self.cdf_function(p)
                self.cdf_inv[p] = self.cdf_inverse(p)

        # If we were pinning any values, undo the changes we did before
        self.cdf_inv.update(tempcdfdict)
        self.pdf.update(temppdfdict)
        self._pdf_norm.update(temppdfnormdict)
        self.prior_pdf.update(temppriordict)

        # Clean out the _rvs arrays for 'irrelevant' points
        #   - create the cumulative weights
        if "weights" in self._rvs and deltaP > 0:
            # Sort the weights with the deltaL masked applied
            sorted_weights = self._rvs["weights"].argsort()
            # Make the (unnormalized) CDF
            total_weight = self._rvs["weights"][sorted_weights].cumsum()
            # Find the deltaP cutoff index
            idx = numpy.searchsorted(total_weight, deltaP*total_weight[-1], 'left')
            sorted_weights = sorted_weights[idx:]
            # Remove all samples which contribute to smallest 1e-3 of cumulative
            # probability
            for key in self._rvs.keys():
                if isinstance(key, tuple):
                    self._rvs[key] = self._rvs[key][:,sorted_weights]
                else:
                    self._rvs[key] = self._rvs[key][sorted_weights]

        if "integrand" in self._rvs and deltalnL is not None:
            deltal_mask = numpy.log(self._rvs["integrand"]) > (maxlnL - deltalnL)
            # Remove all samples which do not have L > maxlnL - deltalnL
            for key in self._rvs.keys():
                if isinstance(key, tuple):
                    self._rvs[key] = self._rvs[key][:,deltal_mask]
                else:
                    self._rvs[key] = self._rvs[key][deltal_mask]

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
        return numpy.asarray((numpy.pi/2-th, ph))

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
        return numpy.asarray((numpy.pi/2-dec, ra))

    def __init__(self, skymap):
        self.skymap = skymap
        self._argsort = skymap.argsort()
        self._probmap = skymap[self._argsort].cumsum()

    def pseudo_pdf(self, dec_in, ra_in):
        """
        Return pixel probability for a given dec_in and ra_in. Note, uses healpy functions to identify correct pixel.
        """
        th, ph = HealPixSampler.decra2thph(dec_in, ra_in)
        res = healpy.npix2nside(len(self.skymap))
        return self.skymap[healpy.ang2pix(res, th, ph)]

    def pseudo_cdf_inverse(self, dec_in=None, ra_in=None, ndraws=1):
        """
        Select points from the skymap with a distribution following its corresponding pixel probability. If dec_in, ra_in are suupplied, they are ignored except that their shape is reproduced. If ndraws is supplied, that will set the shape. Will return a 2xN numpy array of the (dec, ra) values.
        """
        res = healpy.npix2nside(len(self.skymap))

        if ra_in is not None:
            ndraws = len(ra_in)
        elif ndraws is not None:
            ra_in, dec_in = numpy.zeros((2, ndraws))
        else:
            raise ValueError("Either ra / dec or ndraws must be specified.")

        idx = numpy.searchsorted(self._probmap, \
                    numpy.random.uniform(0, 1, ndraws))
        return HealPixSampler.thph2decra(*healpy.pix2ang(res, self._argsort[idx]))

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

if __name__ == "__main__":

    ns = healpy.nside2npix(1024)
    smap = numpy.zeros(ns)
    rnd_idx = numpy.random.randint(len(smap))
    smap[rnd_idx] = 1.

    hps = HealPixSampler(smap)

    pts = hps.pseudo_cdf_inverse(None, ndraws=int(1e6))
    assert numpy.all(healpy.ang2pix(1024, *HealPixSampler.decra2thph(*pts)) == rnd_idx)
