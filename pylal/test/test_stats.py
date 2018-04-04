#!/usr/bin/env python

import unittest

import numpy as np
np.seterr(all="raise")

from pylal.stats import rankdata
from scipy.stats import rankdata as scipy_rankdata

class test_stats(unittest.TestCase):
    def test_identical_continuous(self):
        for _ in xrange(1000):
            data = np.random.uniform(0, 1000, size=1000)
            self.assertTrue((rankdata(data) == scipy_rankdata(data)).all())

    def test_identical_discrete(self):
        for _ in xrange(1000):
            data = np.random.randint(0, 1000, size=1000).astype(float)
            self.assertTrue((rankdata(data) == scipy_rankdata(data)).all())

suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_stats))
unittest.TextTestRunner(verbosity=2).run(suite)
