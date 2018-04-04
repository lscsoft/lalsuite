#!/usr/bin/env python

from itertools import izip
import unittest

import numpy as np
np.seterr(all="raise")
from scipy import integrate

from pylal import moving_histogram as mh
from pylal import rate

num_tests = 100

class test_MovingHistogramFixedN(unittest.TestCase):
    bins = rate.LinearBins(0, 1, 10)
    max_hist_size = 100

    def setUp(self):
        self.hist = mh.MovingHistogramFixedN(self.bins, self.max_hist_size)

    def test_monotonicity_check(self):
        for i in range(num_tests):
            self.setUp()
            a, b = np.random.randint(100, size=2)
            self.hist.update(a, 0.5)
            if b < a:
                self.assertRaises(ValueError, lambda: self.hist.update(b, 0.5))
            else:
                self.hist.update(b, 0.5)

    def test_pdf_normalization(self):
        for i in range(num_tests):
            self.setUp()
            for t, s in enumerate(np.random.random(size=100)):
                self.hist.update(t, s)
            x = np.linspace(0, 1, 100)
            y = np.array([self.hist.get_pdf(a) for a in x])
            integral = integrate.simps(y, x)
            self.assertTrue(abs(integral - 1.0) < 0.01)

    def test_cdf_normalization(self):
        for i in range(num_tests):
            self.setUp()
            for t, s in enumerate(np.random.random(size=100)):
                self.hist.update(t, s)
            self.assertAlmostEqual(self.hist.get_cdf(self.bins.max), 1)

    def test_sf_normalization(self):
        for i in range(num_tests):
            self.setUp()
            for t, s in enumerate(np.random.random(size=100)):
                self.hist.update(t, s)
            self.assertAlmostEqual(self.hist.get_sf(self.bins.min), 1)

    def test_hist_size_limit(self):
        for t, s in enumerate(np.random.random(size=self.max_hist_size + num_tests)):
            self.hist.update(t, s)
            self.assertTrue(len(self.hist) <= self.max_hist_size)

    def test_hist_discards_oldest(self):
        for t, s in enumerate(np.random.random(size=self.max_hist_size + num_tests)):
            self.hist.update(t, s)
            self.assertEqual(self.hist.get_oldest_timestamp(), max(0, t - self.max_hist_size + 1))

    def test_matches_naive_hist(self):
        rand_nums = np.random.random(size=self.max_hist_size + num_tests)
        for t, s in enumerate(rand_nums):
            self.hist.update(t, s)
            naive_hist = np.zeros(len(self.bins), dtype=int)
            for n in rand_nums[max(0, t - self.max_hist_size + 1):t + 1]:
                naive_hist[self.bins[n]] += 1
            self.assertTrue((naive_hist == self.hist.counts).all())


# construct and run the test suite
suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_MovingHistogramFixedN))
unittest.TextTestRunner(verbosity=2).run(suite)
