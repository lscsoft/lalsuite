# Copyright (C) 2010--2013,2015,2016,2018  Kipp Cannon
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
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


import itertools
import warnings


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
from .git_version import date as __date__
from .git_version import version as __version__


#
# =============================================================================
#
#                                Offset Vector
#
# =============================================================================
#


class offsetvector(dict):
	"""
	Subclass of the dict built-in type for storing mappings of
	instrument to time offset.

	Examples:

	>>> x = offsetvector({"H1": 0, "L1": 10, "V1": 20})
	>>> x["H1"]
	0
	>>> not any(x.values())	# check for "zero-lag"
	False

	The motivation for introducing this class, instead of using
	dictionaries, is that it provides a number of tools for comparing
	offset vectors besides strict value-for-value equality.  For
	example the Python cmp() operation compares two offset vectors by
	the relative offsets between instruments rather than their absolute
	offsets, whereas the == operation compares two offset vectors by
	demanding strict equality.  There is also the ability to check if
	one offset vector is a subset of another one.
	"""
	@property
	def refkey(self):
		"""
		= min(self)

		Raises ValueError if the offsetvector is empty.
		"""
		# min() emits ValueError when the list is empty, but it
		# might also emit a ValueError if the comparison operations
		# inside it fail, so we can't simply wrap it in a
		# try/except pair or we might mask genuine failures
		if not self:
			raise ValueError("offsetvector is empty")
		return min(self)

	@property
	def deltas(self):
		"""
		Dictionary of relative offsets.  The keys in the result are
		pairs of keys from the offset vector, (a, b), and the
		values are the relative offsets, (offset[b] - offset[a]).
		Raises ValueError if the offsetvector is empty (WARNING:
		this behaviour might change in the future).

		Example:

		>>> x = offsetvector({"H1": 0, "L1": 10, "V1": 20})
		>>> assert x.deltas == {('H1', 'L1'): 10, ('H1', 'V1'): 20, ('H1', 'H1'): 0}
		>>> y = offsetvector({'H1': 100, 'L1': 110, 'V1': 120})
		>>> y.deltas == x.deltas
		True

		Note that the result always includes a "dummy" entry,
		giving the relative offset of self.refkey with respect to
		itself, which is always 0.

		See also .fromdeltas().

		BUGS:  I think the keys in each tuple should be reversed.
		I can't remember why I put them in the way they are.
		Expect them to change in the future.
		"""
		# FIXME:  instead of raising ValueError when the
		# offsetvector is empty this should return an empty
		# dictionary.  the inverse, .fromdeltas() accepts
		# empty dictionaries
		# NOTE:  the arithmetic used to construct the offsets
		# *must* match the arithmetic used by
		# time_slide_component_vectors() so that the results of the
		# two functions can be compared to each other without worry
		# of floating-point round off confusing things.
		refkey = self.refkey
		refoffset = self[refkey]
		return dict(((refkey, key), self[key] - refoffset) for key in self)

	def __str__(self, compact = False):
		"""
		Return a human-readable string representation of an offset
		vector.

		Example:

		>>> a = offsetvector({"H1": -10.1234567, "L1": 0.125})
		>>> str(a)
		'H1 = -10.1234567 s, L1 = +0.125 s'
		>>> a.__str__(compact = True)
		'H1=-10.123,L1=0.125'
		"""
		if compact:
			return ",".join(("%s=%.5g" % x) for x in sorted(self.items()))
		return ", ".join(("%s = %+.16g s" % x) for x in sorted(self.items()))

	def __repr__(self):
		"""
		Return a string representation of the offset vector.
		Running eval() on the result reconstructs the offsetvector.

		Example:

		>>> a = offsetvector({"H1": -10.1234567, "L1": 0.1})
		>>> b = eval(repr(a))
		>>> b == a
		True
		>>> b is a
		False
		>>> c = offsetvector({"H1": -10.1234567})
		>>> repr(c)
		"offsetvector({'H1': -10.1234567})"
		"""
		return "%s(%s)" % (self.__class__.__name__, dict.__repr__(self))

	def __abs__(self):
		"""
		Returns max(offset) - min(offset).

		Example:

		>>> abs(offsetvector({"H1": 0.0, "H2": 0.0, "L1": 0.0}))
		0.0
		>>> abs(offsetvector({"H1": 10.0, "H2": 10.0, "L1": 10.0}))
		0.0
		>>> abs(offsetvector({'H1': 10.0, 'L1': 0.0, 'V1': -10.0}))
		20.0
		"""
		return max(self.values()) - min(self.values())

	def __cmp__(self, other):
		"""
		Compare two offset vectors by their relative offsets.  The
		return value is 0 if the relative offsets are all equal,
		nonzero otherwise.

		This method is deprecated and will be removed in a future
		release because the cmp() method has been removed from
		Python 3.

		Instead of cmp(a, b), use (a.deltas != b.deltas).

		Example:

		>>> a = offsetvector({"H1": 0.0, "H2": 0.0, "L1": 0.0})
		>>> b = offsetvector({"H1": 10.0, "H2": 10.0, "L1": 10.0})
		>>> # the cmp() method is only present on Python 2:
		>>> import sys
		>>> if sys.version_info.major <= 2:
		...     assert cmp(a, b) == 0
		>>> a == b
		False
		>>> # cmp() is deprecated, use this instead
		>>> a.deltas == b.deltas
		True

		Note that cmp() and testing for equality are different
		tests!  The equality test returns False because the offset
		vectors are not identical, however the cmp() function
		returns 0 because the relative offsets are all equal.
		"""
		warnings.warn("Support for calling cmp(a, b) where a and b are instances of offesetvector is deprecated. Use a.deltas != b.deltas instead.")
		return cmp(self.deltas, other.deltas)

	def contains(self, other):
		"""
		Returns True if offset vector other can be found in self,
		False otherwise.  An offset vector is "found in" another
		offset vector if the latter contains all of the former's
		instruments and the relative offsets among those
		instruments are equal (the absolute offsets need not be).

		Example:

		>>> a = offsetvector({"H1": 10, "L1": 20, "V1": 30})
		>>> b = offsetvector({"H1": 20, "V1": 40})
		>>> a.contains(b)
		True

		Note the distinction between this and the "in" operator:

		>>> "H1" in a
		True
		"""
		return offsetvector((key, offset) for key, offset in self.items() if key in other).deltas == other.deltas

	def normalize(self, **kwargs):
		"""
		Adjust the offsetvector so that a particular instrument has
		the desired offset.  All other instruments have their
		offsets adjusted so that the relative offsets are
		preserved.  The instrument to noramlize, and the offset one
		wishes it to have, are provided as a key-word argument.
		The return value is the time slide dictionary, which is
		modified in place.

		If more than one key-word argument is provided the keys are
		sorted and considered in order until a key is found that is
		in the offset vector.  The offset vector is normalized to
		that value.  This function is a no-op if no key-word
		argument is found that applies.

		Example:

		>>> a = offsetvector({"H1": -10, "H2": -10, "L1": -10})
		>>> assert a.normalize(L1 = 0) == offsetvector({'H2': 0, 'H1': 0, 'L1': 0})
		>>> a = offsetvector({"H1": -10, "H2": -10})
		>>> assert a.normalize(L1 = 0, H2 = 5) == offsetvector({'H2': 5, 'H1': 5})
		"""
		# FIXME:  should it be performed in place?  if it should
		# be, the should there be no return value?
		for key, offset in sorted(kwargs.items()):
			if key in self:
				delta = offset - self[key]
				for key in self.keys():
					self[key] += delta
				break
		return self

	@classmethod
	def fromdeltas(cls, deltas):
		"""
		Construct an offsetvector from a dictionary of offset
		deltas as returned by the .deltas attribute.

		Example:

		>>> x = offsetvector({"H1": 0, "L1": 10, "V1": 20})
		>>> y = offsetvector.fromdeltas(x.deltas)
		>>> y == x
		True

		See also .deltas, .fromkeys()
		"""
		return cls((key, value) for (refkey, key), value in deltas.items())


#
# =============================================================================
#
#                                  Utilities
#
# =============================================================================
#


def component_offsetvectors(offsetvectors, n):
	"""
	Given an iterable of offset vectors, return the shortest list of
	the unique n-instrument offset vectors from which all the vectors
	in the input iterable can be constructed.  This can be used to
	determine the minimal set of n-instrument coincs required to
	construct all of the coincs for all of the requested instrument and
	offset combinations in a set of offset vectors.

	It is assumed that the coincs for the vector {"H1": 0, "H2": 10,
	"L1": 20} can be constructed from the coincs for the vectors {"H1":
	0, "H2": 10} and {"H2": 0, "L1": 10}, that is only the relative
	offsets are significant in determining if two events are
	coincident, not the absolute offsets.
	"""
	#
	# collect unique instrument set / deltas combinations
	#

	delta_sets = {}
	for vect in offsetvectors:
		for instruments in itertools.combinations(sorted(vect), n):
			# NOTE:  the arithmetic used to construct the
			# offsets *must* match the arithmetic used by
			# offsetvector.deltas so that the results of the
			# two can be compared to each other without worry
			# of floating-point round off confusing things.
			delta_sets.setdefault(instruments, set()).add(tuple(vect[instrument] - vect[instruments[0]] for instrument in instruments))

	#
	# translate into a list of normalized n-instrument offset vectors
	#

	return [offsetvector(zip(instruments, deltas)) for instruments, delta_set in delta_sets.items() for deltas in delta_set]
