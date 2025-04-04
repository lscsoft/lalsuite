from numpy.testing import assert_allclose
import sys
import unittest
import warnings

from lalburst import snglcoinc


# FIXME:  make this more strict when we can rely on a new-enough scipy
rtol = 1e-3
def assertDictEqual(x, y, rtol = rtol):
	# keys are the same
	assert set(x) == set(y)
	# values are the same
	for key in x:
		assert abs((x[key] - y[key]) / x[key]) <= rtol


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
			 frozenset(['H1', 'L1', 'V1'])},
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
		crates = rates.coinc_rates(H1=0.001, L1=0.002, V1=0.003)
		expected = {
			frozenset(['H1', 'L1']): 6.00513846088957e-08,
			frozenset(['H1', 'V1']): 1.9372787960306537e-07,
			frozenset(['L1', 'V1']): 3.77380092200718e-07,
			frozenset(['H1', 'L1', 'V1']): 1.0125819710267318e-11,
		}
		for key in crates:
			assert_allclose(crates[key], expected[key], rtol = rtol)

		crates = rates.coinc_rates(H1=0.001, L1=0.002, V1=0.002)
		expected = {
			frozenset(['H1', 'L1']): 6.00513846088957e-08,
			frozenset(['H1', 'V1']): 1.291519197353769e-07,
			frozenset(['L1', 'V1']): 2.5158672813381197e-07,
			frozenset(['H1', 'L1', 'V1']): 6.750546473511545e-12,
		}
		for key in crates:
			assert_allclose(crates[key], expected[key], rtol = rtol)

		crates = rates.coinc_rates(H1=0.001, L1=0.002, V1=0.001)
		expected = {
			frozenset(['H1', 'L1']): 6.00513846088957e-08,
			frozenset(['H1', 'V1']): 6.457595986768845e-08,
			frozenset(['L1', 'V1']): 1.2579336406690598e-07,
			frozenset(['H1', 'L1', 'V1']): 3.3752732367557724e-12,
		}
		for key in crates:
			assert_allclose(crates[key], expected[key], rtol = rtol)

	def test_strict_coinc_rates(self):
		rates = snglcoinc.CoincRates(("H1", "L1", "V1"), 0.005, 2)
		scrates = rates.strict_coinc_rates(H1=0.001, L1=0.002, V1=0.003)
		expected = {
			frozenset(['H1', 'L1']): 6.004125878918543e-08,
			frozenset(['H1', 'V1']): 1.937177537833551e-07,
			frozenset(['L1', 'V1']): 3.7736996638100773e-07,
			frozenset(['H1', 'L1', 'V1']): 1.0125819710267318e-11,
		}
		for key in scrates:
			assert_allclose(scrates[key], expected[key], rtol = rtol)

		scrates = rates.strict_coinc_rates(H1=0.001, L1=0.002, V1=0.002)
		expected = {
			frozenset(['H1', 'L1']): 6.004463406242219e-08,
			frozenset(['H1', 'V1']): 1.2914516918890337e-07,
			frozenset(['L1', 'V1']): 2.5157997758733847e-07,
			frozenset(['H1', 'L1', 'V1']): 6.750546473511545e-12,
		}
		for key in scrates:
			assert_allclose(scrates[key], expected[key], rtol = rtol)

		scrates = rates.strict_coinc_rates(H1=0.001, L1=0.002, V1=0.001)
		expected = {
			frozenset(['H1', 'L1']): 6.004800933565894e-08,
			frozenset(['H1', 'V1']): 6.457258459445168e-08,
			frozenset(['L1', 'V1']): 1.2578998879366924e-07,
			frozenset(['H1', 'L1', 'V1']): 3.3752732367557724e-12,
		}
		for key in scrates:
			assert_allclose(scrates[key], expected[key], rtol = rtol)

	def test_lnP_instruments(self):
		rates = snglcoinc.CoincRates(("H1", "L1", "V1"), 0.005, 2)
		# FIXME:  use the built-in method when a new enough scipy
		# is available
		#self.assertDictEqual(
		assertDictEqual(
			rates.lnP_instruments(H1=0.001, L1=0.002, V1=0.003),
			{frozenset(['H1', 'L1']): -2.352494317162074,
			 frozenset(['H1', 'V1']): -1.181124067253893,
			 frozenset(['L1', 'V1']): -0.5143002401188091,
			 frozenset(['H1', 'L1', 'V1']): -11.040192999777876},
		)


suite = unittest.TestSuite()
suite.addTest(unittest.defaultTestLoader.loadTestsFromTestCase(TestCoincRates))
sys.exit(not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful())
