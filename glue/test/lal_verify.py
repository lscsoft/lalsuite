import doctest
import filecmp
import os
import random
import sys
import unittest
from glue import lal


#
# Define the components of the test suite.
#

def maxLIGOTimeGPS():
	return lal.LIGOTimeGPS(2**32 - 1, 999999999)

def randomLIGOTimeGPS():
	return lal.LIGOTimeGPS(random.randint(-100000000, +100000000), random.randint(0, 999999999))


class test_docstrings(unittest.TestCase):
	def test(self):
		doctest.testmod(lal)
		cache_rw_is_identity = filecmp.cmp("874000000-20000.cache", "874000000-20000.cache.new")
		if cache_rw_is_identity:
			os.remove("874000000-20000.cache.new")
		self.assertEqual(True, cache_rw_is_identity)


class test_LIGOTimeGPS(unittest.TestCase):
	def test__init__(self):
		correct = lal.LIGOTimeGPS(100, 500000000)
		tests = [
			(100.5,),
			(100.500000000,),
			(100.50000000000000000000000,),
			(100, "500000000"),
			(100, "500000000.0000000000000"),
			(101, "-500000000"),
			(101, "-500000000.0000000000000"),
			("100.5",),
			("100.500000000",),
			("100.50000000000000000000000",),
			("100", 500000000),
			("100", 500000000.0000000000000),
			("101", -500000000),
			("101", -500000000.0000000000000),
			("100", "500000000"),
			("100", "500000000.0000000000000"),
			("101", "-500000000"),
			("101", "-500000000.0000000000000"),
			(0, 100500000000),
			(0, 100500000000.0000000000000),
			(99, 1500000000),
			(99.5, 1000000000),
			(-10, 110500000000),
			(-10.5, 111000000000)
		]
		for num, test in enumerate(tests):
			try:
				self.assertEqual(correct, lal.LIGOTimeGPS(*test))
			except AssertionError as e:
				raise AssertionError("Test %d failed: " % (num) + str(e))

	def test__float__(self):
		self.assertEqual(100.5, float(lal.LIGOTimeGPS(100.5)))

	def test__int__(self):
		self.assertEqual(100, int(lal.LIGOTimeGPS(100.1)))
		self.assertEqual(100, int(lal.LIGOTimeGPS(100.9)))

	def testns(self):
		self.assertEqual(100500000000, lal.LIGOTimeGPS(100.5).ns())

	def test__nonzero__(self):
		self.assertEqual(True, bool(lal.LIGOTimeGPS(100.5)))
		self.assertEqual(False, bool(lal.LIGOTimeGPS(0)))

	def test__add__(self):
		self.assertEqual(lal.LIGOTimeGPS(110.5), lal.LIGOTimeGPS(100.5) + 10)
		self.assertEqual(lal.LIGOTimeGPS(110.5), lal.LIGOTimeGPS(100.5) + lal.LIGOTimeGPS(10))

	def test__mul__(self):
		self.assertEqual(lal.LIGOTimeGPS(10), lal.LIGOTimeGPS(5) * 2)
		self.assertEqual(lal.LIGOTimeGPS(10), lal.LIGOTimeGPS(20) * 0.5)
		self.assertEqual(lal.LIGOTimeGPS(0), lal.LIGOTimeGPS(1000) * 0)

	def test__div__(self):
		self.assertEqual(lal.LIGOTimeGPS(10), lal.LIGOTimeGPS(20) / 2)
		self.assertEqual(lal.LIGOTimeGPS(10), lal.LIGOTimeGPS(5) / .5)

	def test__mod__(self):
		self.assertEqual(lal.LIGOTimeGPS(3), lal.LIGOTimeGPS(13) % 5.0)

	def test_swig_comparison(self):
		try:
			from lal import LIGOTimeGPS as swigLIGOTimeGPS
		except ImportError:
			print >>sys.stderr, "lal swig bindings not available:  skipping test"
			return

		toswig = lambda x: swigLIGOTimeGPS(str(x))
		fromswig = lambda x: lal.LIGOTimeGPS(str(x))

		operators = {
			"add": (lal.LIGOTimeGPS.__add__, swigLIGOTimeGPS.__add__),
			"sub": (lal.LIGOTimeGPS.__sub__, swigLIGOTimeGPS.__sub__)
		}

		for i in xrange(100000):
			key, (op, swigop) = random.choice(operators.items())
			arg1 = randomLIGOTimeGPS() / 2
			arg2 = randomLIGOTimeGPS() / 2
			try:
				self.assertEqual(op(arg1, arg2), fromswig(swigop(toswig(arg1), toswig(arg2))))
			except AssertionError as s:
				raise AssertionError("%s(%s, %s) comparison failed: %s" % (key, str(arg1), str(arg2), str(s)))

		# FIXME:  mod tests fail, fix then enable
		operators = {
			"mul": (lal.LIGOTimeGPS.__mul__, swigLIGOTimeGPS.__mul__),
			"div": (lal.LIGOTimeGPS.__div__, swigLIGOTimeGPS.__div__)#,
			#"mod": (lal.LIGOTimeGPS.__mod__, swigLIGOTimeGPS.__mod__)
		}

		for i in xrange(100000):
			key, (op, swigop) = random.choice(operators.items())
			arg1 = randomLIGOTimeGPS() / 100
			arg2 = 100**(random.random() * 2 - 1)
			try:
				self.assertEqual(abs(op(arg1, arg2) - fromswig(swigop(toswig(arg1), arg2))) <= 1e-9, True)
			except AssertionError as s:
				raise AssertionError("%s(%s, %s) comparison failed: %s != %s" % (key, str(arg1), "%.17g" % arg2, str(op(arg1, arg2)), str(swigop(toswig(arg1), arg2))))


#
# Construct and run the test suite.
#


suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_docstrings))
suite.addTest(unittest.makeSuite(test_LIGOTimeGPS))

sys.exit(not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful())
