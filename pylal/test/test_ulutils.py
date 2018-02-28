#!/usr/bin/env python

import unittest
import numpy
from numpy import random

from pylal import upper_limit_utils

class test_ulutils(unittest.TestCase):

    def test_normalize_pdf(self):
        '''
        Check the normalization of pdfs
        '''
        mu = numpy.logspace(-10,-4,1e5)
        likely = numpy.exp(-mu*3.14*10**6)  # just some number far from unity
        junk, prob_norm = upper_limit_utils.normalize_pdf(mu, likely)
        dmu = mu[1:] - mu[:-1]
        prob_in_bin = (prob_norm[1:] + prob_norm[:-1]) /2
        prob_integral = sum(dmu*prob_in_bin)
        self.assertTrue( abs(prob_integral - 1.0) < 0.00001 )

    def test_gaussian_upper_limit(self):
        '''
        Give the upper_limit function some known distributions with known 95% upper limits.
        '''
        # Gaussian over *positive* mu: 95% limit is equal to sqrt(2)*erfc^-1(0.05)
        mu = numpy.linspace(0,5,1e6)
        post = numpy.exp(-(mu**2)/2)
        muhi = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.95)
        self.assertTrue( abs(muhi - 1.95996) < 0.0001 ) # get upper limit to 4 decimal places

    def test_exponential_upper_limit(self):
        # Exponential
        mu = numpy.linspace(0,15,1e6)
        post = numpy.exp(-mu)
        muhi = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.95)
        self.assertTrue( abs( muhi - numpy.log(20) ) < 0.0001 ) # get upper limit to 4 decimal places

    def test_uniform_upper_limit(self):
        # Uniform posterior
        mumax = 15
        mu = numpy.linspace(0,mumax,1e6)
        post = numpy.ones(len(mu))
        alphas = numpy.arange(0.1,1,0.1)
        for a in alphas:
            muhi = upper_limit_utils.compute_upper_limit(mu, post, alpha = a)
            self.assertTrue( a*mumax - 0.0001 < muhi < a*mumax + 0.0001) # get upper limit to 4 decimal places

    def test_logspacing_rate_upper_limit(self):
        '''
        Give the upper_limit function some known distributions with known 95% upper limits.
        '''
        # Exponential
        mu = numpy.logspace(-5,2,1e6)
        post = numpy.exp(-mu)
        muhi = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.95)
        self.assertTrue( abs( muhi - numpy.log(20) ) < 0.0001 ) # get upper limit to 4 decimal places

        # Gaussian over positive mu
        post = numpy.exp(-(mu**2)/2)
        muhi = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.95)
        self.assertTrue( abs(muhi - 1.95996) < 0.0001 ) # get upper limit to 4 decimal places

    def test_volume_lambda(self):
        '''
        Check the dependence of upper limits on volume and lambda.
        '''
        # volumes to test: these range over 5 orders of magnitude
        volumes = numpy.logspace(-3,2,50)

        for vol in volumes:
            # take a large, fixed set of mu samples to bracket the 1/vol values
            mu = numpy.logspace(-5,5,1e5)

            # lambda = 0
            likely = upper_limit_utils.margLikelihood([vol], [0], mu)
            post = likely  # uniform prior; NB compute_upper_limit works for unnormalized posteriors
            muhi = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.90)
            self.assertTrue( 2.30/vol < muhi < 2.31/vol )

            # lambda = 1
            likely = upper_limit_utils.margLikelihood([vol], [1], mu)
            post = likely  # uniform prior
            muhi = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.90)
            self.assertTrue( 3.27/vol < muhi < 3.28/vol )

            # lambda = infinity
            likely = upper_limit_utils.margLikelihood([vol], [1e6], mu)
            post = likely  # uniform prior
            muhi = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.90)
            self.assertTrue( 3.885/vol < muhi < 3.895/vol )

    def test_zero_volume_search(self):
        '''
        Check that running a search with zero volume has
        no effect on upper limits.
        '''
        # no prior specified, unit volume
        mu = numpy.linspace(0,100,1e4)
        likely = upper_limit_utils.margLikelihood([1], [0], mu)
        post = likely/likely.sum() #uniform prior
        muhi_np = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.90)

        # posterior for prior, zero volume
        likely = upper_limit_utils.margLikelihood([0],[0], mu)
        post *= likely  #uniform post for prior
        post /= post.sum()
        muhi_zv = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.90)
        self.assertTrue( abs(muhi_zv - muhi_np) < 0.01 )

    def test_volume_additivity(self):
        '''
        For a series of independent searches with 0 lambda,
        the search volumes add. Check that upper limits computed
        by summing volumes is same as upper limit obtained by
        iteratively applying the individual searches.
        '''
        # volume additivity check
        mu = numpy.linspace(0,100,1e4)
        likely = upper_limit_utils.margLikelihood([1], [0], mu)
        post = likely/likely.sum() #uniform prior
        post *= upper_limit_utils.margLikelihood([9], [0], mu) # use post for prior
        post /= post.sum()
        muhi_1p9 = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.90)

        likely = upper_limit_utils.margLikelihood([10], [0], mu)
        post = likely/likely.sum() #uniform prior
        muhi_10 = upper_limit_utils.compute_upper_limit(mu, post, alpha = 0.90)
        self.assertTrue( abs( muhi_1p9 - muhi_10 ) < 0.01 )

    def test_integrate_efficiency_lind(self):
        '''
        Check that the numerical accuracy of the integration
        methods are sufficient.
        '''
        xbins = numpy.linspace(0,1,50)
        centres = (xbins[:-1]+xbins[1:])/2
        mockeff = numpy.exp(-centres)/(4*numpy.pi*centres**2)
        v, verr = upper_limit_utils.integrate_efficiency(xbins, mockeff, logbins=False)
        vexpect = 1 - numpy.exp(-1)
        self.assertTrue(abs(v -vexpect ) < 0.01)

    def test_integrate_efficiency_logd(self):
        '''
        Check that the numerical accuracy of the integration
        methods are sufficient.
        '''
        xbins = numpy.logspace(-2,0,50)
        centres = numpy.exp((numpy.log(xbins[1:])+numpy.log(xbins[:-1]))/2) # log midpoint
        mockeff = numpy.exp(-centres)/(4*numpy.pi*centres**2)
        v, verr = upper_limit_utils.integrate_efficiency(xbins, mockeff, logbins=True)
        vexpect = 1 - numpy.exp(-1)
        self.assertTrue(abs(v -vexpect ) < 0.01)

    def test_integrate_realistic_efficiency(self):
        '''
        Check that the volume calculation gives the expected results in the
        ideal case of an analytically known efficiency curve.
        '''
        Rmax = 150
        def mockeff(r):
            # This efficiency curve has the general features
            # of a real efficiency curve seen empirically.
            # It also integrates nicely and so the results of
            # integrating this efficiency curve numerically can be
            # compared against an analytic reference.
            return (1-r/Rmax)**4

        rbins = numpy.linspace(0,Rmax,20)
        centres = (rbins[1:]+rbins[:-1])/2
        v, verr = upper_limit_utils.integrate_efficiency(rbins, mockeff(centres))
        vexpect = (1./35)*(4*numpy.pi/3)*(Rmax**3)
        self.assertTrue(abs(1-v/vexpect) < 0.01)

    def test_mean_efficiency(self):
        '''
        Check the mean efficiency calculation in a controlled way.
        '''
        # Assume a sigmoidal mock efficiency curve.
        def eff_model(inj, rchar = 25.0, order = 6.0):
            return 1./(1+(inj/rchar)**order)

        # think of these as log-d bins between 1 and 100 Mpc
        bins = numpy.logspace(0, 2, num=50)
        centres = 10**((numpy.log10(bins[1:])+numpy.log10(bins[:-1]))/2) # log midpoint

        # generate a bunch of injections uniform in log-d
        injections = random.uniform(1, 100, size=10000)
        class MiniInj(object):
            def __init__(self, distance):
                self.distance = distance

        found = []
        missed = []
        for inj in injections:
            if random.binomial(1, eff_model(inj) ):
                found.append( MiniInj(inj) )
            else:
                missed.append( MiniInj(inj) )

        eff, err, meanvol, volerr = upper_limit_utils.mean_efficiency_volume(found, missed, bins)
        # the computed mean efficiency should agree with the efficiency
        # model to at least ~5% (though this can fluctuate)
        self.assertTrue( (eff - eff_model(centres)).sum()/len(eff) < 0.05 )

    def test_compute_many_posterior(self):
        # for 0 lambda's, volumes add
        mu = numpy.linspace(0,100,1e4)
        post1 = upper_limit_utils.margLikelihood([5,10,4,6],[0,0,0,0],mu)
        post2 = upper_limit_utils.margLikelihood([15,10],[0,0],mu)
        post3 = upper_limit_utils.margLikelihood([25],[0],mu)

        mu_90_1 = upper_limit_utils.compute_upper_limit(mu,post1/post1.sum(),0.90)
        mu_90_2 = upper_limit_utils.compute_upper_limit(mu,post2/post2.sum(),0.90)
        mu_90_3 = upper_limit_utils.compute_upper_limit(mu,post3/post3.sum(),0.90)

        self.assertTrue( (mu_90_2-mu_90_1) < 0.05 )
        self.assertTrue( (mu_90_3-mu_90_1) < 0.05 )
        self.assertTrue( (mu_90_3-mu_90_2) < 0.05 )

    def test_zero_padded_uniform(self):
        # Uniform posterior
        mumax = 15
        mu = numpy.linspace(0,2*mumax,2*1e6)
        post = numpy.concatenate((numpy.ones(len(mu)/2),numpy.zeros(len(mu)/2)))
        a = 1
        muhi = upper_limit_utils.compute_upper_limit(mu, post, alpha = 1)
        self.assertTrue( a*mumax - 0.0001 < muhi < a*mumax + 0.0001) # get upper limit to 4 sig figs

    def test_confidence_interval(self):
        '''
        The confidence interval function computes a minimal width interval.
        For a symmetric distribution, this interval is the same as computing
        upper and lower limits at half-confidences. For a monotonically
        decreasing distribution, this interval is the same as the upper limit
        at the same confidence level.
        '''
        exes = numpy.arange(-10,10,0.01)
        post = numpy.exp(-exes**2/2)
        lo, hi = upper_limit_utils.confidence_interval_min_width(exes, post, 0.9)
        lo2 = upper_limit_utils.compute_lower_limit(exes, post, 0.95)
        hi2 = upper_limit_utils.compute_upper_limit(exes, post, 0.95)

        self.assertTrue( abs(lo - lo2) < 0.001 )
        self.assertTrue( abs(hi - hi2) < 0.001 )


# construct and run the test suite.
suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_ulutils))
unittest.TextTestRunner(verbosity=2).run(suite)
