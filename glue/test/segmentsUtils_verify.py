import random
import StringIO
import sys
import unittest

from glue import iterutils
from glue import segments
from glue import segmentsUtils
import verifyutils
verifyutils.segments = segments


#
#  How many times to repeat the algebraic tests
#


algebra_repeats = 8000
algebra_listlength = 200


#
# Define the components of the test suite.
#


class test_segwizard(unittest.TestCase):
	def test_fromsegwizard(self):
		"""
		Test segwizard parsing.
		"""
		data = StringIO.StringIO("""# This is a comment
 # This is another comment
	# Again a comment
1  10 100 90
2 110 120 10# Here's a comment
3 125 130 5 # Another one

4   0 200 200""")
		correct = segments.segmentlist([segments.segment(10, 100), segments.segment(110, 120), segments.segment(125, 130), segments.segment(0, 200)])
		self.assertEqual(correct, segmentsUtils.fromsegwizard(data, strict=True))

	def test_tofromseqwizard(self):
		"""
		Check that the segwizard writing routine's output is parsed
		correctly.
		"""
		data = StringIO.StringIO()
		correct = segments.segmentlist([segments.segment(10, 100), segments.segment(110, 120), segments.segment(125, 130), segments.segment(0, 200)])
		segmentsUtils.tosegwizard(data, correct)
		data.seek(0)
		self.assertEqual(correct, segmentsUtils.fromsegwizard(data, strict=True))


class test_vote(unittest.TestCase):
	def test_vote(self):
		"""
		Test vote().
		"""
		for i in range(algebra_repeats):
			seglists = []
			for j in range(random.randint(0, 10)):
				seglists.append(verifyutils.random_coalesced_list(algebra_listlength))
			n = random.randint(0, len(seglists))
			correct = reduce(lambda x, y: x | y, (votes and reduce(lambda a, b: a & b, votes) or segments.segmentlist() for votes in iterutils.choices(seglists, n)), segments.segmentlist())
			self.assertEqual(correct, segmentsUtils.vote(seglists, n))


#
# Construct and run the test suite.
#


suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(test_segwizard))
suite.addTest(unittest.makeSuite(test_vote))

sys.exit(not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful())
