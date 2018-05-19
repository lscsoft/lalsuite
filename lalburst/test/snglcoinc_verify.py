import unittest
import sys

from lalburst import snglcoinc


class TestCoincRates(unittest.TestCase):
    def test_all_instrument_combos(self):
        rates = snglcoinc.CoincRates(("H1", "L1", "V1"), 0.005, 1)
        self.assertSetEqual(
            set(rates.all_instrument_combos),
            {frozenset(['H1']),
             frozenset(['L1']),
             frozenset(['V1']),
             frozenset(['H1', 'L1']),
             frozenset(['H1', 'V1']),
             frozenset(['L1', 'V1']),
             frozenset(['H1', 'V1', 'L1'])},
        )

        rates = snglcoinc.CoincRates(("H1", "L1", "V1"), 0.005, 2)
        self.assertSetEqual(
            set(rates.all_instrument_combos),
            {frozenset(['H1', 'L1']),
             frozenset(['H1', 'V1']),
             frozenset(['L1', 'V1']),
             frozenset(['H1', 'L1', 'V1'])},
        )

    def test_coinc_rates(self):
        rates = snglcoinc.CoincRates(("H1", "L1", "V1"), 0.005, 2)
        self.assertDictEqual(
            rates.coinc_rates(H1=0.001, L1=0.002, V1=0.003),
            {frozenset(['H1', 'L1']): 6.00513846088957e-08,
             frozenset(['H1', 'V1']): 1.9372787960306537e-07,
             frozenset(['L1', 'V1']): 3.77380092200718e-07,
             frozenset(['H1', 'L1', 'V1']): 1.0125819710267318e-11},
        )
        self.assertDictEqual(
            rates.coinc_rates(H1=0.001, L1=0.002, V1=0.002),
            {frozenset(['H1', 'L1']): 6.00513846088957e-08,
             frozenset(['H1', 'V1']): 1.291519197353769e-07,
             frozenset(['L1', 'V1']): 2.5158672813381197e-07,
             frozenset(['H1', 'L1', 'V1']): 6.750546473511545e-12},
        )
        self.assertDictEqual(
            rates.coinc_rates(H1=0.001, L1=0.002, V1=0.001),
            {frozenset(['H1', 'L1']): 6.00513846088957e-08,
             frozenset(['H1', 'V1']): 6.457595986768845e-08,
             frozenset(['L1', 'V1']): 1.2579336406690598e-07,
             frozenset(['H1', 'L1', 'V1']): 3.3752732367557724e-12},
        )

    def test_strict_coinc_rates(self):
        rates = snglcoinc.CoincRates(("H1", "L1", "V1"), 0.005, 2)
        self.assertDictEqual(
            rates.strict_coinc_rates(H1=0.001, L1=0.002, V1=0.003),
            {frozenset(['H1', 'L1']): 6.004125878918543e-08,
             frozenset(['H1', 'V1']): 1.937177537833551e-07,
             frozenset(['L1', 'V1']): 3.7736996638100773e-07,
             frozenset(['H1', 'L1', 'V1']): 1.0125819710267318e-11},
        )
        self.assertDictEqual(
            rates.strict_coinc_rates(H1=0.001, L1=0.002, V1=0.002),
            {frozenset(['H1', 'L1']): 6.004463406242219e-08,
             frozenset(['H1', 'V1']): 1.2914516918890337e-07,
             frozenset(['L1', 'V1']): 2.5157997758733847e-07,
             frozenset(['H1', 'L1', 'V1']): 6.750546473511545e-12},
        )
        self.assertDictEqual(
            rates.strict_coinc_rates(H1=0.001, L1=0.002, V1=0.001),
            {frozenset(['H1', 'L1']): 6.004800933565894e-08,
             frozenset(['H1', 'V1']): 6.457258459445168e-08,
             frozenset(['L1', 'V1']): 1.2578998879366924e-07,
             frozenset(['H1', 'L1', 'V1']): 3.3752732367557724e-12},
        )

    def test_lnP_instruments(self):
        rates = snglcoinc.CoincRates(("H1", "L1", "V1"), 0.005, 2)
        self.assertDictEqual(
            rates.lnP_instruments(H1=0.001, L1=0.002, V1=0.003),
            {frozenset(['H1', 'L1']): -2.352494317162074,
             frozenset(['H1', 'V1']): -1.181124067253893,
             frozenset(['L1', 'V1']): -0.5143002401188091,
             frozenset(['H1', 'L1', 'V1']): -11.040192999777876},
        )


suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(TestCoincRates))
sys.exit(not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful())
