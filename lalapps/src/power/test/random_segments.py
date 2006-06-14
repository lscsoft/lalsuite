"""
Write a random segment list in segwizard format to stdout.
"""

import random
import sys
from glue import segments
from glue import segmentsUtils

def randomlist(n):
	def r():
		return random.randint(1,5000)
	if n < 1:
		raise ValueError, "randomlist(n): n must be >= 1"
	x = r()
	l = segments.segmentlist([segments.segment(x, x + r())])
	for i in range(n - 1):
		x = l[-1][1] + r()
		l.append(segments.segment(x, x + r()))
	return l

segmentsUtils.tosegwizard(sys.stdout, randomlist(1000))
