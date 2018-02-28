#!/usr/bin/env python

import unittest
import numpy
import os

from pylal import Fr

class test_Fr(unittest.TestCase):
    def test_1d_default_roundtrip(self):
        """ roundtrip test call with default values """
        a = Fr.frgetvect1d("./test.dat","Adc1")
        Fr.frputvect('writetest.gwf', [{'name':'Adc1', 'data':a[0],
            'start':a[1], 'dx':a[3], 'kind':'ADC', 'x_unit':a[4],
            'y_unit':a[5]}])
        b = Fr.frgetvect1d("writetest.gwf", "Adc1")
        self.assert_(numpy.alltrue(a[0] == b[0]))
        self.assert_(numpy.alltrue(a[1:] == b[1:]))
        os.remove("writetest.gwf")

    def test_1d_keywords_roundtrip(self):
        """ roundtrip test call with keyword arguments """
        a = Fr.frgetvect1d("./test.dat", "Adc1", span=1)
        Fr.frputvect('writetest.gwf', [{'name':'Adc1', 'data':a[0],
        'start':a[1], 'dx':a[3], 'kind':'ADC', 'x_unit':a[4],
        'y_unit':a[5]}])
        b = Fr.frgetvect1d("writetest.gwf", "Adc1")
        self.assert_(numpy.alltrue(a[0] == b[0]))
        self.assert_(numpy.alltrue(a[1:] == b[1:]))
        os.remove("writetest.gwf")

    def test_1d_two_channels_roundtrip(self):
        """ roundtrip test call with two channels in a frame """
        a = Fr.frgetvect1d("./test.dat","Adc1")
        Fr.frputvect('writetest.gwf', [{'name':'Adc1', 'data':a[0],
        'start':a[1], 'dx':a[3], 'kind':'ADC', 'x_unit':a[4],
        'y_unit':a[5]},{'name':'reverse', 'data':a[0][::-1], 'start':a[1],
        'dx':a[3], 'kind':'ADC', 'x_unit':a[4], 'y_unit': a[5]}])
        b = Fr.frgetvect1d("writetest.gwf", "Adc1")
        self.assert_(numpy.alltrue(a[0] == b[0]))
        self.assert_(numpy.alltrue(a[1:] == b[1:]))

        c = Fr.frgetvect1d("writetest.gwf", "reverse")
        self.assert_(numpy.alltrue(a[0][::-1] == c[0]))
        self.assert_(numpy.alltrue(a[1:] == c[1:]))
        os.remove("writetest.gwf")

    def test_frgetevent_known(self):
        """ test that we can pull the known contents from a test MBTA file """
        a = Fr.frgetevent("MbtaFake-930909680-16.gwf")
        self.assertEqual(len(a), 1)
        self.assertEqual(a[0]["name"], "MbtaHLV_Chi2OK")
        self.assertAlmostEqual(a[0]["H1:end_time"], 930493014.980663, 6)
        self.assertAlmostEqual(a[0]["H1:SNR"], 9.0)
        self.assertAlmostEqual(a[0]["H1:mass1"], 3.14696, 5)
        self.assertAlmostEqual(a[0]["H1:mass2"], 10.2256, 4)
        self.assertAlmostEqual(a[0]["H1:eff_distance"], 230.718, 3)
        self.assertAlmostEqual(a[0]["L1:end_time"], 930493014.984246, 6)
        self.assertAlmostEqual(a[0]["L1:SNR"], 8.0)
        self.assertAlmostEqual(a[0]["L1:mass1"], 3.44696, 5)
        self.assertAlmostEqual(a[0]["L1:mass2"], 10.5256, 4)
        self.assertAlmostEqual(a[0]["L1:eff_distance"], 582.025, 3)
        self.assertAlmostEqual(a[0]["V1:end_time"], 930493014.989135, 6)
        self.assertAlmostEqual(a[0]["V1:SNR"], 7.0)
        self.assertAlmostEqual(a[0]["V1:mass1"], 2.84696, 5)
        self.assertAlmostEqual(a[0]["V1:mass2"], 9.92556, 5)
        self.assertAlmostEqual(a[0]["V1:eff_distance"], 200.113, 3)

# construct and run the test suite.
suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_Fr))
unittest.TextTestRunner(verbosity=2).run(suite)
