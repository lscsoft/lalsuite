import doctest
from pylal import date
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
import unittest


#
# Define the components of the test suite.
#

class test_LIGOTimeGPS(unittest.TestCase):
	def test__init__(self):
		correct = LIGOTimeGPS(100, 500000000)
		tests = [
			(100.5,),
			(100.500000000,),
			(100.50000000000000000000000,),
			("100.5",),
			("100.500000000",),
			("100.50000000000000000000000",),
			(0, 100500000000),
			(99, 1500000000),
			(-10, 110500000000),
		]
		for num, test in enumerate(tests):
			try:
				self.assertEqual(correct, LIGOTimeGPS(*test))
			except AssertionError, e:
				raise AssertionError, "Test %d failed: " % (num) + str(e)

	def test__float__(self):
		self.assertEqual(100.5, float(LIGOTimeGPS(100.5)))

	def test__int__(self):
		self.assertEqual(100, int(LIGOTimeGPS(100.1)))
		self.assertEqual(100, int(LIGOTimeGPS(100.9)))

	def test__nonzero__(self):
		self.assertEqual(True, bool(LIGOTimeGPS(100.5)))
		self.assertEqual(False, bool(LIGOTimeGPS(0)))

	def test__add__(self):
		self.assertEqual(LIGOTimeGPS(110.5), LIGOTimeGPS(100.5) + LIGOTimeGPS(10))

	def test__mul__(self):
		self.assertEqual(LIGOTimeGPS(10), LIGOTimeGPS(5) * 2)
		self.assertEqual(LIGOTimeGPS(10), LIGOTimeGPS(20) * 0.5)

	def test__div__(self):
		self.assertEqual(LIGOTimeGPS(10), LIGOTimeGPS(20) / 2)
		self.assertEqual(LIGOTimeGPS(10), LIGOTimeGPS(5) / 0.5)
		self.assertEqual(LIGOTimeGPS(333333333,333333333), LIGOTimeGPS(1000000000) / 3)
		self.assertEqual(LIGOTimeGPS(666666666,666666667), LIGOTimeGPS(2000000000) / 3)

	def test__mod__(self):
		self.assertEqual(LIGOTimeGPS(3), LIGOTimeGPS(13) % 5.0)


#
# Construct and run the test suite.
#

suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_LIGOTimeGPS))

unittest.TextTestRunner(verbosity=2).run(suite)

doctest.testmod(date)
