import random


#
# segmentlist generation utilities
#


def random_coalesced_list(n):
	"""
	Return a coalesced segmentlist of n elements with random boundaries.
	"""
	def r():
		return random.randint(1, 127) / 128.0
	if n < 1:
		raise ValueError, n
	l = segments.segmentlist([None] * n)
	x = r()
	l[0] = segments.segment(x, x + r())
	x = l[0][1] + r()
	for i in xrange(1, n):
		l[i] = segments.segment(x, x + r())
		x = l[i][1] + r()
	return l


def random_uncoalesced_list(n):
	"""
	Return an uncoalesced segmentlist of n elements with random
	boundaries.
	"""
	def r():
		return float(random.randint(1, 999)) / 1000
	if n < 1:
		raise ValueError, n
	x = r()
	l = segments.segmentlist([segments.segment(x, x + r() / 100.0)])
	for i in xrange(n - 1):
		x = r()
		l.append(segments.segment(x, x + r() / 100.0))
	return l


def iscoalesced(l):
	"""
	Return True if the segmentlist l is coalesced.
	"""
	for a, b in zip(l, l[1:]):
		if a[1] >= b[0]:
			return False
	return True
