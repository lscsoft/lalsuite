# Copyright (C) 2006  Kipp Cannon
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


#
# NOTE:  the logic in this code is unintuitively complicated.  Small,
# apparently irrelevant, changes to conditionals can have subtly unexpected
# consequences to the behaviour of the class methods.  ALWAYS make sure that
# the test suite returns OK on ALL tests after any changes you make.
#


"""
This module defines the segment and segmentlist objects, as well as the
infinity object used to define semi-infinite and infinite segments.

See also:

glue.segmentsUtils
"""


from bisect import bisect_left as _bisect_left
from bisect import bisect_right as _bisect_right
from copy import copy as _shallowcopy


from glue import git_version


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                   infinity
#
# =============================================================================
#


class infinity(object):
	"""
	The infinity object possess the algebraic properties necessary for
	use as a bound on semi-infinite and infinite segments.

	This class uses comparison-by-identity rather than
	comparison-by-value.  What this means, is there are only ever two
	instances of this class, representing positive and negative
	infinity respectively.  All other "instances" of this class are
	infact simply references to one of these two, and comparisons are
	done by checking which one you've got.  This improves speed and
	reduces memory use, and is similar in implementation to Python's
	boolean True and False objects.

	The normal way to obtain references to positive or negative
	infinity is to do infinity() or -infinity() respectively.  It is
	also possible to select the sign by passing a single numeric
	argument to the constructor.  The sign of the argument causes a
	reference to either positive or negative infinity to be returned,
	respectively.  For example infinity(-1) is equivalent to
	-infinity().  However, this feature is a little slower and not
	recommended for normal use;  it is provided only to simplify the
	pickling and unpickling of instances of the class.

	Example:

	>>> x = infinity()
	>>> x > 0
	True
	>>> x + 10 == x
	True
	>>> segment(-10, 10) - segment(-x, 0)
	segment(0, 10)
	"""
	__slots__ = []

	def __new__(cls, *args):
		if args:
			# pickle support
			(sign,) = args
			if sign > 0:
				return PosInfinity
			if sign < 0:
				return NegInfinity
			raise ValueError(sign)
		return PosInfinity

	def __repr__(self):
		"""
		Returns a string.
		"""
		if self is PosInfinity:
			return "infinity"
		return "-infinity"

	__str__ = __repr__

	# pickle support

	def __reduce__(self):
		"""
		Pickle support.
		"""
		if self is PosInfinity:
			return (infinity, (1,))
		# self is NegInfinity
		return (infinity, (-1,))

	# tests

	def __cmp__(self, other):
		"""
		Positive infinity compares as greater than everything
		except itself, negative infinity compares as less than
		everything except itself.
		"""
		if self is other:
			return 0
		if self is PosInfinity:
			return 1
		# self is NegInfinity
		return -1

	def __nonzero__(self):
		"""
		Returns True.
		"""
		return True

	# arithmetic

	def __neg__(self):
		"""
		Returns -self.
		"""
		if self is PosInfinity:
			return NegInfinity
		# self is NegInfinity
		return PosInfinity

	def __pos__(self):
		"""
		Returns self.
		"""
		return self

	def __add__(self, other):
		"""
		Returns self.
		"""
		return self

	def __radd__(self, other):
		"""
		Returns self.
		"""
		return self

	def __sub__(self, other):
		"""
		Returns self.
		"""
		return self

	def __rsub__(self, other):
		"""
		Returns -self.
		"""
		if self is PosInfinity:
			return NegInfinity
		# self is NegInfinity
		return PosInfinity


PosInfinity = object.__new__(infinity)
NegInfinity = object.__new__(infinity)


#
# =============================================================================
#
#                                   segment
#
# =============================================================================
#


class segment(tuple):
	"""
	The segment class defines objects that represent a range of values.
	A segment has a start and an end, and is taken to represent the
	range of values in the semi-open interval [start, end).  Some
	limited arithmetic operations are possible with segments, but
	because the set of (single) segments is not closed under the
	sensible definitions of the standard arithmetic operations, the
	behaviour of the arithmetic operators on segments may not be as you
	would expect.  For general arithmetic on segments, use segmentlist
	objects.  The methods for this class exist mostly for purpose of
	simplifying the implementation of the segmentlist class.

	The segment class is a subclass of the tuple built-in class
	provided by Python.  This means segments are immutable --- you
	cannot modify a segment object after creating it, to change the
	boundaries of a segment you must create a new segment object with
	the desired boundaries.  Like tuples, segments can be used as
	dictionary keys, and like tuples the comparison used to find a
	segment in the dictionary is done by value not by ID.  And, like
	tuples, a segment can be created from any sequence-like object by
	passing it to the constructor (the sequence must have exactly two
	elements in it).

	Example:

	>>> segment(0, 10) & segment(5, 15)
	segment(5, 10)
	>>> segment(0, 10) | segment(5, 15)
	segment(0, 15)
	>>> segment(0, 10) - segment(5, 15)
	segment(0, 5)
	>>> segment(0, 10) < segment(5, 15)
	True
	>>> segment(1, 2) in segment(0, 10)
	True
	>>> bool(segment(0, 0))
	False
	>>> segment("AAA Towing", "York University") & segment("Pool", "Zoo")
	segment('Pool', 'York University')
	>>> x = [0, 1]
	>>> segment(x)
	segment(0, 1)
	>>> y = segment(0, 1)
	>>> y == x
	True
	>>> y is x
	False
	>>> x in y
	True
	>>> z = {x: ["/path/to/file1", "/path/to/file2"]}
	>>> y in z
	True
	>>> z[y]
	['/path/to/file1', '/path/to/file2']
	"""

	# basic class methods

	def __new__(cls, *args):
		if len(args) == 1:
			args = args[0]
		if len(args) != 2:
			raise TypeError("__new__() takes 2 arguments, or 1 argument when it is a sequence of length 2")
		if args[0] <= args[1]:
			return tuple.__new__(cls, args)
		else:
			return tuple.__new__(cls, (args[1], args[0]))

	def __repr__(self):
		return "segment(" + repr(self[0]) + ", " + repr(self[1]) + ")"

	def __str__(self):
		return "[" + str(self[0]) + " ... " + str(self[1]) + ")"

	# accessors

	def __abs__(self):
		"""
		Returns the length of the interval represented by the
		segment.  Requires the bounds to support the subtract
		operation.
		"""
		return self[1] - self[0]

	# comparisons

	def __nonzero__(self):
		"""
		Return True if the segment's boudaries are not equal, False
		if they are equal.
		"""
		return self[0] != self[1]

	def disjoint(self, other):
		"""
		Returns >0 if self covers an interval above other's
		interval, <0 if self covers an interval below other's, or 0
		if the two intervals are not disjoint (intersect or touch).
		A return value of 0 indicates the two segments would
		coalesce.
		"""
		if self[0] > other[1]:
			return 1
		if self[1] < other[0]:
			return -1
		return 0

	def __lt__(self, other):
		if isinstance(other, tuple):
			return tuple.__lt__(self, other)
		return self[0] < other

	def __le__(self, other):
		if isinstance(other, tuple):
			return tuple.__le__(self, other)
		return self[0] <= other

	def __eq__(self, other):
		if isinstance(other, tuple):
			return tuple.__eq__(self, other)
		return self[0] == other

	def __ne__(self, other):
		if isinstance(other, tuple):
			return tuple.__ne__(self, other)
		return self[0] != other

	def __gt__(self, other):
		if isinstance(other, tuple):
			return tuple.__gt__(self, other)
		return self[0] > other

	def __ge__(self, other):
		if isinstance(other, tuple):
			return tuple.__ge__(self, other)
		return self[0] >= other

	# some arithmetic operations that (mostly) make sense for segments

	def __and__(self, other):
		"""
		Return the segment that is the intersection of the given
		segments.  Raises ValueError if the result cannot be
		presented as a single segment.
		"""
		if (self[1] <= other[0]) or (self[0] >= other[1]):
			# self and other don't intersect
			raise ValueError(other)
		return tuple.__new__(self.__class__, (max(self[0], other[0]), min(self[1], other[1])))

	def __or__(self, other):
		"""
		Return the segment that is the union of the given segments.
		Raises ValueError if the result cannot be represented as a
		single segment.
		"""
		if (self[1] < other[0]) or (self[0] > other[1]):
			# self and other are disjoint
			raise ValueError(other)
		return tuple.__new__(self.__class__, (min(self[0], other[0]), max(self[1], other[1])))

	# addition is union
	__add__ = __or__

	def __sub__(self, other):
		"""
		Return the segment that is that part of self which is not
		contained in other.  Raises ValueError if the result cannot
		be represented as a single segment.
		"""
		if (self[1] <= other[0]) or (self[0] >= other[1]):
			# self and other do not intersect
			return self
		if (self in other) or ((self[0] < other[0]) and (self[1] > other[1])):
			# result is not exactly 1 segment
			raise ValueError(other)
		if self[0] < other[0]:
			return tuple.__new__(self.__class__, (self[0], other[0]))
		return tuple.__new__(self.__class__, (other[1], self[1]))

	# check for proper intersection and subsetness

	def intersects(self, other):
		"""
		Return True if the intersection of self and other is not a
		null segment.
		"""
		return (self[1] > other[0]) and (self[0] < other[1])

	def __contains__(self, other):
		"""
		Return True if other is wholly contained in self.  If other
		is an instance of the segment class or an instance of a
		subclass of segment then it is treated as an interval whose
		upper and lower bounds must not be outside of self,
		otherwise other is compared to the bounds of self as a
		scalar.
		"""
		try:
			a, b = other
		except ValueError:
			return self[0] <= other < self[1]
		else:
			return (self[0] <= a) and (self[1] >= b)

	# protraction and contraction and shifting

	def protract(self, x):
		"""
		Return a new segment whose bounds are given by subtracting
		x from the segment's lower bound and adding x to the
		segment's upper bound.
		"""
		return self.__class__(self[0] - x, self[1] + x)

	def contract(self, x):
		"""
		Return a new segment whose bounds are given by adding x to
		the segment's lower bound and subtracting x from the
		segment's upper bound.
		"""
		return self.__class__(self[0] + x, self[1] - x)

	def shift(self, x):
		"""
		Return a new segment whose bounds are given by adding x to
		the segment's upper and lower bounds.
		"""
		return tuple.__new__(self.__class__, (self[0] + x, self[1] + x))


#
# =============================================================================
#
#                                 segmentlist
#
# =============================================================================
#


class segmentlist(list):
	"""
	The segmentlist class defines a list of segments, and is an
	extension of the built-in list class.  This class provides
	addtional methods that assist in the manipulation of lists of
	segments.  In particular, arithmetic operations such as union and
	intersection are provided.  Unlike the segment class, the
	segmentlist class is closed under all supported arithmetic
	operations.

	All standard Python sequence-like operations are supported, like
	slicing, iteration and so on, including arithmetic operations.
	However, the arithmetic and other methods for this class generally
	require the segmentlist to be in what is refered to as a
	"coalesced" state --- consisting solely of disjoint segments listed
	in ascending order.  Using the standard Python sequence-like
	operations, a segmentlist can be easily constructed that is not in
	this state;  for example by simply appending a segment to the end
	of the list that overlaps some other segment already in the list.
	The use of methods that require coalesced lists with lists that are
	not coalesced has undefined results.  The class provides the
	.coalesce() method that can be called to put a segmentlist in the
	coalesced state.  All arithmetic methods return coalesced results,
	so typically the .coalesce() method will be executed once after
	importing a segmentlist from an untrusted source, then there is
	never a need to call the .coalesce() method again as long as the
	segmentlists are manipulated exclusively via the arithmetic
	operators.

	Example:

	>>> x = segmentlist([segment(-10, 10)])
	>>> x |= segmentlist([segment(20, 30)])
	>>> x -= segmentlist([segment(-5, 5)])
	>>> print x
	[segment(-10, -5), segment(5, 10), segment(20, 30)]
	>>> print ~x
	[segment(-infinity, -10), segment(-5, 5), segment(10, 20), segment(30, infinity)]
	"""

	# container method over-rides.

	def __contains__(self, item):
		"""
		Returns True if the given object is wholly contained within
		the segments in self.  If self has length n, then if item
		is a scalar or a segment this operation is O(log n), if
		item is a segmentlist of m segments this operation is O(m
		log n).

		Note the difference between this operator and the standard
		Python "in" operator for sequence-like objects:  in the
		case of standard sequence-like objects the in operator
		checks for an exact match between the given item and one of
		the contents of the list; for segmentlists, the in operator
		checks if the given item is contained within any of the
		segments in the segmentlist.
		"""
		if isinstance(item, self.__class__):
			return all(seg in self for seg in item)
		i = _bisect_left(self, item)
		return ((i != 0) and (item in self[i-1])) or ((i != len(self)) and (item in self[i]))

	# supplementary accessors

	def __abs__(self):
		"""
		Return the sum of the durations of all segments in self.
		Does not require the segmentlist to be coalesced.
		"""
		return sum(abs(seg) for seg in self)

	def extent(self):
		"""
		Return the segment whose end-points denote the maximum and
		minimum extent of the segmentlist.  Does not require the
		segmentlist to be coalesced.
		"""
		if not len(self):
			raise ValueError("empty list")
		min, max = self[0]
		for lo, hi in self:
			if min > lo:
				min = lo
			if max < hi:
				max = hi
		return segment(min, max)

	def find(self, item):
		"""
		Return the smallest i such that i is the index of an
		element that wholly contains item.  Raises ValueError if no
		such element exists.  Does not require the segmentlist to
		be coalesced.
		"""
		for i, seg in enumerate(self):
			if item in seg:
				return i
		raise ValueError(item)

	# arithmetic operations that are sensible with segment lists

	def __iand__(self, other):
		"""
		Replace the segmentlist with the intersection of itself and
		another.  If the two lists have lengths n and m
		respectively, this operation is O(n + m).
		"""
		return self.__isub__(~other)

	def __and__(self, other):
		"""
		Return the intersection of the segmentlist and another.  If
		the two lists have lengths n and m respectively, this
		operation is O(n + m).
		"""
		if len(self) >= len(other):
			return self.__class__(self).__iand__(other)
		return self.__class__(other).__iand__(self)

	def __ior__(self, other):
		"""
		Replace the segmentlist with the union of itself and
		another.  If the two lists have numbers of elements n and m
		respectively, then for m << n the algorithm is O(m log n),
		otherwise it is O((n + m) log (n + m)).
		"""
		if len(other) > len(self) / 2:
			self.extend(other)
			return self.coalesce()
		if other is self:
			return self
		i = 0
		for seg in other:
			i = j = _bisect_right(self, seg, i)
			lo, hi = seg
			if i and self[i - 1][1] >= lo:
				i -= 1
				lo = self[i][0]
			n = len(self)
			while j < n and self[j][0] <= hi:
				j += 1
			if j > i:
				self[i] = segment(lo, max(hi, self[j - 1][1]))
				del self[i + 1 : j]
			else:
				self.insert(i, seg)
			i += 1
		return self

	def __or__(self, other):
		"""
		Return the union of the segment list and another.  The
		algorithm has the same scaling as in the
		segmentlist.__ior__() method, except that the lists are
		reordered to attempt to use the O(m log n) case.
		"""
		if len(self) >= len(other):
			return self.__class__(self).__ior__(other)
		return self.__class__(other).__ior__(self)

	def __xor__(self, other):
		"""
		Return the segmentlist that is the list of all intervals
		contained in exactly one of this and another list.  This
		operation is O(n log n).
		"""
		l = self - other
		l.extend(other - self)
		l.sort()
		return l

	# addition is union
	__iadd__ = __ior__
	__add__ = __or__

	def __isub__(self, other):
		"""
		Replace the segmentlist with the difference between itself
		and another.  For lists of length m and n respectively,
		this operation is O(n + m).
		"""
		if not other:
			return self
		if other is self:
			del self[:]
			return self
		i = j = 0
		other_lo, other_hi = other[j]
		while i < len(self):
			self_lo, self_hi = self[i]
			while other_hi <= self_lo:
				j += 1
				if j >= len(other):
					return self
				other_lo, other_hi = other[j]
			if self_hi <= other_lo:
				i += 1
			elif other_lo <= self_lo:
				if other_hi >= self_hi:
					del self[i]
				else:
					self[i] = segment(other_hi, self_hi)
			else:
				self[i] = segment(self_lo, other_lo)
				i += 1
				if other_hi < self_hi:
					self.insert(i, segment(other_hi, self_hi))
		return self

	def __sub__(self, other):
		"""
		Return the difference between the segmentlist and another.
		This operation is O(n).
		"""
		return self.__class__(self).__isub__(other)

	def __invert__(self):
		"""
		Return the segmentlist that is the inversion of the given
		list.  This operation is O(n).
		"""
		if not len(self):
			return self.__class__([segment(NegInfinity, PosInfinity)])
		l = self.__class__()
		if self[0][0] > NegInfinity:
			l.append(segment(NegInfinity, self[0][0]))
		last = self[0][1]
		for i in xrange(1, len(self)):
			l.append(segment(last, self[i][0]))
			last = self[i][1]
		if last < PosInfinity:
			l.append(segment(last, PosInfinity))
		return l

	# other operations

	def intersects_segment(self, other):
		"""
		Returns True if the intersection of self and the segment
		other is not the null set, otherwise returns False.  The
		algorithm is O(log n).  Requires the list to be coalesced.
		"""
		i = _bisect_left(self, other)
		return ((i != 0) and (other[0] < self[i-1][1])) or ((i != len(self)) and (other[1] > self[i][0]))

	def intersects(self, other):
		"""
		Returns True if the intersection of self and the
		segmentlist other is not the null set, otherwise returns
		False.  The algorithm is O(n), but faster than explicit
		calculation of the intersection, i.e. by testing bool(self
		& other).  Requires both lists to be coalesced.
		"""
		# if either has zero length, the answer is False
		if not (self and other):
			return False
		# walk through both lists in order, searching for a match
		i = j = 0
		seg = self[0]
		otherseg = other[0]
		while True:
			if seg[1] <= otherseg[0]:
				i += 1
				if i >= len(self):
					return False
				seg = self[i]
			elif otherseg[1] <= seg[0]:
				j += 1
				if j >= len(other):
					return False
				otherseg = other[j]
			else:
				return True

	def coalesce(self):
		"""
		Sort the elements of the list into ascending order, and merge
		continuous segments into single segments.  Segmentlist is
		modified in place.  This operation is O(n log n).
		"""
		self.sort()
		i = j = 0
		n = len(self)
		while j < n:
			lo, hi = self[j]
			j += 1
			while j < n and hi >= self[j][0]:
				hi = max(hi, self[j][1])
				j += 1
			if lo != hi:
				self[i] = segment(lo, hi)
				i += 1
		del self[i : ]
		return self

	def protract(self, x):
		"""
		Execute the .protract() method on each segment in the list
		and coalesce the result.  Segmentlist is modified in place.
		"""
		for i in xrange(len(self)):
			self[i] = self[i].protract(x)
		return self.coalesce()

	def contract(self, x):
		"""
		Execute the .contract() method on each segment in the list
		and coalesce the result.  Segmentlist is modified in place.
		"""
		for i in xrange(len(self)):
			self[i] = self[i].contract(x)
		return self.coalesce()

	def shift(self, x):
		"""
		Execute the .shift() method on each segment in the list.
		The algorithm is O(n) and does not require the list to be
		coalesced nor does it coalesce the list.  Segmentlist is
		modified in place.
		"""
		for i in xrange(len(self)):
			self[i] = self[i].shift(x)
		return self


#
# =============================================================================
#
#                               segmentlistdict
#
# =============================================================================
#


class _offsets(dict):
	"""
	Implements the segmentlist offset book-keeping in the
	segmentlistdict class.  Not intended for use outside of the
	segmentlistdict class.
	"""
	def __new__(cls, parent):
		return dict.__new__(cls)

	def __init__(self, parent):
		dict.__init__(self)
		self.__parent = parent

	def __reduce__(self):
		return _offsets, (self.__parent,), None, None, iter(self.items())

	def __setitem__(self, key, value):
		"""
		Set an offset.  If the new offset is identical to the
		current offset this is a no-op, otherwise the corresponding
		segmentlist object is shifted.
		"""
		try:
			delta = value - self[key]
		except KeyError:
			dict.__setitem__(self, key, value)
			return
		if delta:
			self.__parent[key].shift(delta)
			dict.__setitem__(self, key, self[key] + delta)

	def update(self, d):
		"""
		From a dictionary of offsets, apply each offset to the
		corresponding segmentlist.  NOTE:  it is acceptable for the
		offset dictionary to contain entries for which there is no
		matching segmentlist; no error will be raised, but the
		offset will be ignored.  This simplifies the case of
		updating several segmentlistdict objects from a common
		offset dictionary, when one or more of the segmentlistdicts
		contains only a subset of the keys.
		"""
		for key, value in d.iteritems():
			if key in self:
				self[key] = value

	def clear(self):
		"""
		Remove the offsets from all segmentlists.
		"""
		for key in self:
			self[key] = 0.0

	# stubs to prevent bugs
	def __delitem__(*args):
		raise NotImplementedError
	def fromkeys(*args):
		raise NotImplementedError
	def pop(*args):
		raise NotImplementedError
	def popitem(*args):
		raise NotImplementedError


class segmentlistdict(dict):
	"""
	A dictionary associating a unique label and numeric offset with
	each of a set of segmentlist objects.

	This class implements a standard mapping interface, with additional
	features added to assist with the manipulation of a collection of
	segmentlist objects.  In particular, methods for taking unions and
	intersections of the lists in the dictionary are available, as well
	as the ability to record and apply numeric offsets to the
	boundaries of the segments in each list.

	The numeric offsets are stored in the "offsets" attribute, which
	itself is a dictionary, associating a number with each key in the
	main dictionary.  Assigning to one of the entries of the offsets
	attribute has the effect of shifting the corresponding segmentlist
	from its original position (not its current position) by the given
	amount.

	Example:

	>>> x = segmentlistdict()
	>>> x["H1"] = segmentlist([segment(0, 10)])
	>>> print x
	{'H1': [segment(0, 10)]}
	>>> x.offsets["H1"] = 6
	>>> print x
	{'H1': [segment(6.0, 16.0)]}
	>>> x.offsets.clear()
	>>> print x
	{'H1': [segment(0.0, 10.0)]}
	>>> x["H2"] = segmentlist([segment(5, 15)])
	>>> x.intersection(["H1", "H2"])
	[segment(5, 10.0)]
	>>> x.offsets["H1"] = 6
	>>> x.intersection(["H1", "H2"])
	[segment(6.0, 15)]
	>>> c = x.extract_common(["H1", "H2"])
	>>> c.offsets.clear()
	>>> c
	{'H2': [segment(6.0, 15)], 'H1': [segment(0.0, 9.0)]}
	"""
	def __new__(cls, *args):
		self = dict.__new__(cls, *args)
		self.offsets = _offsets(self)
		return self

	def __init__(self, *args):
		dict.__init__(self, *args)
		dict.clear(self.offsets)
		for key in self:
			dict.__setitem__(self.offsets, key, 0.0)
		if args and isinstance(args[0], self.__class__):
			dict.update(self.offsets, args[0].offsets)

	def copy(self, keys = None):
		"""
		Return a copy of the segmentlistdict object.  The return
		value is a new object with a new offsets attribute, with
		references to the original keys, and shallow copies of the
		segment lists.  Modifications made to the offset dictionary
		or segmentlists in the object returned by this method will
		not affect the original, but without using much memory
		until such modifications are made.  If the optional keys
		argument is not None, then should be an iterable of keys
		and only those segmentlists will be copied (KeyError is
		raised if any of those keys are not in the
		segmentlistdict).

		More details.  There are two "built-in" ways to create a
		copy of a segmentlist object.  The first is to initialize a
		new object from an existing one with

		>>> old = segmentlistdict()
		>>> new = segmentlistdict(old)

		This creates a copy of the dictionary, but not of its
		contents.  That is, this creates new with references to the
		segmentlists in old, therefore changes to the segmentlists
		in either new or old are reflected in both.  The second
		method is

		>>> new = old.copy()

		This creates a copy of the dictionary and of the
		segmentlists, but with references to the segment objects in
		the original segmentlists.  Since segments are immutable,
		this effectively creates a completely independent working
		copy but without the memory cost of a full duplication of
		the data.
		"""
		if keys is None:
			keys = self
		new = self.__class__()
		for key in keys:
			new[key] = _shallowcopy(self[key])
			dict.__setitem__(new.offsets, key, self.offsets[key])
		return new

	def __setitem__(self, key, value):
		"""
		Set the segmentlist associated with a key.  If key is not
		already in the dictionary, the corresponding offset is
		initialized to 0.0, otherwise it is left unchanged.
		"""
		dict.__setitem__(self, key, value)
		if key not in self.offsets:
			dict.__setitem__(self.offsets, key, 0.0)

	def __delitem__(self, key):
		dict.__delitem__(self, key)
		dict.__delitem__(self.offsets, key)

	# supplementary accessors

	def map(self, func):
		"""
		Return a dictionary of the results of func applied to each
		of the segmentlist objects in self.

		Example:

		>>> x = segmentlistdict()
		>>> x["H1"] = segmentlist([segment(0, 10)])
		>>> x["H2"] = segmentlist([segment(5, 15)])
		>>> x.map(lambda l: 12 in l)
		{'H2': True, 'H1': False}
		"""
		return dict((key, func(value)) for key, value in self.iteritems())

	def __abs__(self):
		"""
		Return a dictionary of the results of running .abs() on
		each of the segmentlists.
		"""
		return self.map(abs)

	def extent(self):
		"""
		Return a dictionary of the results of running .extent() on
		each of the segmentlists.
		"""
		return self.map(segmentlist.extent)

	def extent_all(self):
		"""
		Return the result of running .extent() on the union of all
		lists in the dictionary.
		"""
		segs = tuple(seglist.extent() for seglist in self.values() if seglist)
		if not segs:
			raise ValueError("empty list")
		return segment(min(seg[0] for seg in segs), max(seg[1] for seg in segs))

	def find(self, item):
		"""
		Return a dictionary of the results of running .find() on
		each of the segmentlists.

		Example:

		>>> x = segmentlistdict()
		>>> x["H1"] = segmentlist([segment(0, 10)])
		>>> x["H2"] = segmentlist([segment(5, 15)])
		>>> x.find(7)
		{'H2': 0, 'H1': 0}

		NOTE:  all segmentlists must contain the item or KeyError
		is raised.
		"""
		return self.map(lambda x: x.find(item))

	def keys_at(self, x):
		"""
		Return a list of the keys for the segment lists that
		contain x.

		Example:

		>>> x = segmentlistdict()
		>>> x["H1"] = segmentlist([segment(0, 10)])
		>>> x["H2"] = segmentlist([segment(5, 15)])
		>>> x.keys_at(12)
		['H2']
		"""
		return [key for key, segs in self.items() if x in segs]

	# list-by-list arithmetic

	def __iand__(self, other):
		for key, value in other.iteritems():
			if key in self:
				self[key] &= value
			else:
				self[key] = segmentlist()
		return self

	def __and__(self, other):
		if sum(len(s) for s in self.values()) <= sum(len(s) for s in other.values()):
			return self.copy().__iand__(other)
		return other.copy().__iand__(self)

	def __ior__(self, other):
		for key, value in other.iteritems():
			if key in self:
				self[key] |= value
			else:
				self[key] = _shallowcopy(value)
		return self

	def __or__(self, other):
		if sum(len(s) for s in self.values()) >= sum(len(s) for s in other.values()):
			return self.copy().__ior__(other)
		return other.copy().__ior__(self)

	__iadd__ = __ior__
	__add__ = __or__

	def __isub__(self, other):
		for key, value in other.iteritems():
			if key in self:
				self[key] -= value
		return self

	def __sub__(self, other):
		return self.copy().__isub__(other)

	def __ixor__(self, other):
		for key, value in other.iteritems():
			if key in self:
				self[key] ^= value
			else:
				self[key] = _shallowcopy(value)
		return self

	def __xor__(self, other):
		if sum(len(s) for s in self.values()) <= sum(len(s) for s in other.values()):
			return self.copy().__ixor__(other)
		return other.copy().__ixor__(self)

	def __invert__(self):
		new = self.copy()
		for key, value in new.items():
			dict.__setitem__(new, key, ~value)
		return new

	# other list-by-list operations

	def intersects_segment(self, seg):
		"""
		Returns True if any segmentlist in self intersects the
		segment, otherwise returns False.
		"""
		return any(value.intersects_segment(seg) for value in self.itervalues())

	def intersects(self, other):
		"""
		Returns True if there exists a segmentlist in self that
		intersects the corresponding segmentlist in other;  returns
		False otherwise.

		See also:

		.intersects_all(), .all_intersects(), .all_intersects_all()
		"""
		return any(key in self and self[key].intersects(value) for key, value in other.iteritems())

	def intersects_all(self, other):
		"""
		Returns True if each segmentlist in other intersects the
		corresponding segmentlist in self;  returns False
		if this is not the case, or if other is empty.

		See also:

		.intersects(), .all_intersects(), .all_intersects_all()
		"""
		return all(key in self and self[key].intersects(value) for key, value in other.iteritems()) and bool(other)

	def all_intersects(self, other):
		"""
		Returns True if each segmentlist in self intersects the
		corresponding segmentlist in other;  returns False
		if this is not the case or if self is empty.

		See also:

		.intersects, .intersects_all(), .all_intersects_all()
		"""
		return all(key in other and other[key].intersects(value) for key, value in self.iteritems()) and bool(self)

	def all_intersects_all(self, other):
		"""
		Returns True if self and other have the same keys, and each
		segmentlist intersects the corresponding segmentlist in the
		other;  returns False if this is not the case or if either
		dictionary is empty.

		See also:

		.intersects(), .all_intersects(), .intersects_all()
		"""
		return set(self) == set(other) and all(other[key].intersects(value) for key, value in self.iteritems()) and bool(self)

	def extend(self, other):
		"""
		Appends the segmentlists from other to the corresponding
		segmentlists in self, adding new segmentslists to self as
		needed.
		"""
		for key, value in other.iteritems():
			if key not in self:
				self[key] = _shallowcopy(value)
			else:
				self[key].extend(value)

	def coalesce(self):
		"""
		Run .coalesce() on all segmentlists.
		"""
		for value in self.itervalues():
			value.coalesce()
		return self

	def contract(self, x):
		"""
		Run .contract(x) on all segmentlists.
		"""
		for value in self.itervalues():
			value.contract(x)
		return self

	def protract(self, x):
		"""
		Run .protract(x) on all segmentlists.
		"""
		for value in self.itervalues():
			value.protract(x)
		return self

	def extract_common(self, keys):
		"""
		Return a new segmentlistdict containing only those
		segmentlists associated with the keys in keys, with each
		set to their mutual intersection.  The offsets are
		preserved.
		"""
		keys = set(keys)
		new = self.__class__()
		intersection = self.intersection(keys)
		for key in keys:
			dict.__setitem__(new, key, _shallowcopy(intersection))
			dict.__setitem__(new.offsets, key, self.offsets[key])
		return new

	# multi-list operations

	def is_coincident(self, other, keys = None):
		"""
		Return True if any segment in any list in self intersects
		any segment in any list in other.  If the optional keys
		argument is not None, then it should be an iterable of keys
		and only segment lists for those keys will be considered in
		the test (instead of raising KeyError, keys not present in
		both segment list dictionaries will be ignored).  If keys
		is None (the default) then all segment lists are
		considered.

		This method is equivalent to the intersects() method, but
		without requiring the keys of the intersecting segment
		lists to match.
		"""
		if keys is not None:
			keys = set(keys)
			self = tuple(self[key] for key in set(self) & keys)
			other = tuple(other[key] for key in set(other) & keys)
		else:
			self = tuple(self.values())
			other = tuple(other.values())
		# make sure inner loop is smallest
		if len(self) < len(other):
			self, other = other, self
		return any(a.intersects(b) for a in self for b in other)

	def intersection(self, keys):
		"""
		Return the intersection of the segmentlists associated with
		the keys in keys.
		"""
		keys = set(keys)
		if not keys:
			return segmentlist()
		seglist = _shallowcopy(self[keys.pop()])
		for key in keys:
			seglist &= self[key]
		return seglist

	def union(self, keys):
		"""
		Return the union of the segmentlists associated with the
		keys in keys.
		"""
		keys = set(keys)
		if not keys:
			return segmentlist()
		seglist = _shallowcopy(self[keys.pop()])
		for key in keys:
			seglist |= self[key]
		return seglist


#
# =============================================================================
#
#                          Use C Version If Possible
#
# =============================================================================
#


try:
	from __segments import *
except ImportError:
	pass


#
# =============================================================================
#
#                                Pickle Support
#
# =============================================================================
#


import copy_reg

copy_reg.pickle(segment, lambda x: (segment, tuple(x)))
copy_reg.pickle(segmentlist, lambda x: (segmentlist, (), None, iter(x)))
