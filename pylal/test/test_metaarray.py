#!/usr/bin/env python

import operator
import unittest

import numpy as np
np.seterr(all="raise")

from glue.segments import segment, segmentlist
from pylal import metaarray

TimeSeries_metadata = {
    "name": "test",
    "dt": 0.1,
    "segments": segmentlist([segment(1, 5)]),
    "comments": []
}
class test_TimeSeries(unittest.TestCase):
    def test_identity(self):
        """
        See that the TimeSeries wrapping doesn't touch array data
        """
        arr = np.arange(100, dtype=np.float32)
        spec = metaarray.TimeSeries(arr, TimeSeries_metadata)
        self.assertTrue((arr == spec.A).all())
        arr = np.arange(100, dtype=np.float64)
        spec = metaarray.TimeSeries(arr, TimeSeries_metadata)
        self.assertTrue((arr == spec.A).all())
        arr = np.ones(100, dtype=np.bool8)
        spec = metaarray.TimeSeries(arr, TimeSeries_metadata)
        self.assertTrue((arr == spec.A).all())
        arr = np.arange(100, dtype=np.int32)
        spec = metaarray.TimeSeries(arr, TimeSeries_metadata)
        self.assertTrue((arr == spec.A).all())
        arr = np.arange(100, dtype=np.int64)
        spec = metaarray.TimeSeries(arr, TimeSeries_metadata)
        self.assertTrue((arr == spec.A).all())
        arr = np.arange(100, dtype=np.int64)
        spec = metaarray.TimeSeries(arr, TimeSeries_metadata)
        self.assertTrue((arr == spec.A).all())

    def test_multiply_type(self):
        """
        See that TimeSeries metadata are retained by a multiplication.
        """
        a = np.random.random(1000)
        b = np.random.random(1000)

        a_TimeSeries = metaarray.TimeSeries(a, TimeSeries_metadata)
        b_TimeSeries = metaarray.TimeSeries(b, TimeSeries_metadata)

        # multiply
        self.assertTrue(isinstance(a_TimeSeries * b_TimeSeries, metaarray.TimeSeries))
        self.assertTrue(isinstance(a * b_TimeSeries, metaarray.TimeSeries))
        self.assertTrue(isinstance(a_TimeSeries * b, metaarray.TimeSeries))

        # in-place multiply
        c = np.array(a, copy=True)
        c *= b_TimeSeries
        self.assertTrue(isinstance(c, np.ndarray))
        d = metaarray.TimeSeries(a_TimeSeries, copy=True)
        d *= b
        self.assertTrue(isinstance(d, metaarray.TimeSeries))

    def test_multiply_correctness(self):
        a = np.random.random(1000)
        b = np.random.random(1000)

        a_TimeSeries = metaarray.TimeSeries(a, TimeSeries_metadata)
        b_TimeSeries = metaarray.TimeSeries(b, TimeSeries_metadata)

        # multiply
        self.assertTrue((a * b == (a_TimeSeries * b_TimeSeries).A).all())
        self.assertTrue((a * b == (a * b_TimeSeries).A).all())
        self.assertTrue((a * b == (a_TimeSeries * b).A).all())

        # in-place multiply
        c = np.array(a, copy=True)
        c *= b_TimeSeries
        self.assertTrue((a * b == c).all())
        d = metaarray.TimeSeries(a_TimeSeries, copy=True)
        d *= b
        self.assertTrue((a * b == d.A).all())

    def test_compatible_metadata_merge(self):
        a = metaarray.TimeSeries(np.arange(1000), TimeSeries_metadata)
        b = metaarray.TimeSeries(np.arange(1000), TimeSeries_metadata)
        c = a * b
        self.assertTrue(c.metadata.dt == TimeSeries_metadata["dt"])
        self.assertTrue(c.metadata.comments == TimeSeries_metadata["comments"])
        self.assertTrue(c.metadata.segments == TimeSeries_metadata["segments"])

    def test_incompatible_metadata_merge(self):
        a = metaarray.TimeSeries(np.arange(1000), TimeSeries_metadata)
        b = metaarray.TimeSeries(np.arange(1000), TimeSeries_metadata)

        b.metadata.dt = 2 * a.metadata.dt
        self.assertRaises(AssertionError, operator.mul, a, b)

Spectrum_metadata = {
    "name": "test",
     "df": 0.1,
     "f_low": 40,
     "segments": segmentlist([segment(1, 5)]),
     "comments": []
}
class test_Spectrum(unittest.TestCase):
    def test_identity(self):
        """
        See that the Spectrum wrapping doesn't touch array data
        """
        arr = np.arange(100, dtype=np.float32)
        spec = metaarray.Spectrum(arr, Spectrum_metadata)
        self.assertTrue((arr == spec.A).all())
        arr = np.arange(100, dtype=np.float64)
        spec = metaarray.Spectrum(arr, Spectrum_metadata)
        self.assertTrue((arr == spec.A).all())
        arr = np.ones(100, dtype=np.bool8)
        spec = metaarray.Spectrum(arr, Spectrum_metadata)
        self.assertTrue((arr == spec.A).all())
        arr = np.arange(100, dtype=np.int32)
        spec = metaarray.Spectrum(arr, Spectrum_metadata)
        self.assertTrue((arr == spec.A).all())
        arr = np.arange(100, dtype=np.int64)
        spec = metaarray.Spectrum(arr, Spectrum_metadata)
        self.assertTrue((arr == spec.A).all())
        arr = np.arange(100, dtype=np.int64)
        spec = metaarray.Spectrum(arr, Spectrum_metadata)
        self.assertTrue((arr == spec.A).all())

    def test_multiply_type(self):
        """
        See that Spectrum metadata are retained by a multiplication.
        """
        a = np.random.random(1000)
        b = np.random.random(1000)

        a_Spectrum = metaarray.Spectrum(a, Spectrum_metadata)
        b_Spectrum = metaarray.Spectrum(b, Spectrum_metadata)

        # multiply
        self.assertTrue(isinstance(a_Spectrum * b_Spectrum, metaarray.Spectrum))
        self.assertTrue(isinstance(a * b_Spectrum, metaarray.Spectrum))
        self.assertTrue(isinstance(a_Spectrum * b, metaarray.Spectrum))

        # in-place multiply
        c = np.array(a, copy=True)
        c *= b_Spectrum
        self.assertTrue(isinstance(c, np.ndarray))
        d = metaarray.Spectrum(a_Spectrum, copy=True)
        d *= b
        self.assertTrue(isinstance(d, metaarray.Spectrum))

    def test_multiply_correctness(self):
        a = np.random.random(1000)
        b = np.random.random(1000)

        a_Spectrum = metaarray.Spectrum(a, Spectrum_metadata)
        b_Spectrum = metaarray.Spectrum(b, Spectrum_metadata)

        # multiply
        self.assertTrue((a * b == (a_Spectrum * b_Spectrum).A).all())
        self.assertTrue((a * b == (a * b_Spectrum).A).all())
        self.assertTrue((a * b == (a_Spectrum * b).A).all())

        # in-place multiply
        c = np.array(a, copy=True)
        c *= b_Spectrum
        self.assertTrue((a * b == c).all())
        d = metaarray.Spectrum(a_Spectrum, copy=True)
        d *= b
        self.assertTrue((a * b == d.A).all())

    def test_compatible_metadata_merge(self):
        a = metaarray.Spectrum(np.arange(1000), Spectrum_metadata)
        b = metaarray.Spectrum(np.arange(1000), Spectrum_metadata)
        c = a * b
        self.assertTrue(c.metadata.df == Spectrum_metadata["df"])
        self.assertTrue(c.metadata.f_low == Spectrum_metadata["f_low"])
        self.assertTrue(c.metadata.comments == Spectrum_metadata["comments"])
        self.assertTrue(c.metadata.segments == Spectrum_metadata["segments"])

    def test_incompatible_metadata_merge(self):
        a = metaarray.Spectrum(np.arange(1000), Spectrum_metadata)
        b = metaarray.Spectrum(np.arange(1000), Spectrum_metadata)

        b.metadata.df = 2 * a.metadata.df
        self.assertRaises(AssertionError, operator.mul, a, b)
        b.metadata.df = a.metadata.df

        b.metadata.f_low = a.metadata.f_low + 10
        self.assertRaises(AssertionError, operator.mul, a, b)
        b.metadata.f_low = a.metadata.f_low


# construct and run the test suite
suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_Spectrum))
suite.addTest(unittest.makeSuite(test_TimeSeries))
unittest.TextTestRunner(verbosity=2).run(suite)
