import doctest
import filecmp
import os
import sys
import unittest

from lal.utils import cache


#
# Define the components of the test suite.
#


class test_docstrings(unittest.TestCase):
	def test(self):
		failures = doctest.testmod(cache)[0]
		if failures:
			sys.exit(bool(failures))
		cache_rw_is_identity = filecmp.cmp(os.path.join(os.environ.get("LAL_TEST_SRCDIR", "."), "874000000-20000.cache"), "874000000-20000.cache.new")
		if cache_rw_is_identity:
			os.remove("874000000-20000.cache.new")
		self.assertEqual(True, cache_rw_is_identity)


#
# Construct and run the test suite.
#


suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_docstrings))

sys.exit(not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful())
