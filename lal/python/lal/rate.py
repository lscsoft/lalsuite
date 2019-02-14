# Copyright (C) 2006--2017  Kipp Cannon
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


"""
This module provides facilities for studying impulsive events.  A number of
multi-dimensional binning functions are provided, as well as code to
convolve binned data with integral- and phase-preserving window functions
to produce smoothed representations of data sets.  This is particularly
well suited for use in computing moving-average rate data from collections
of impulsive events, elliminating binning artifacts from histograms, and
smoothing contour plots.
"""


from functools import reduce
try:
	from fpconst import PosInf, NegInf
except ImportError:
	# fpconst is not part of the standard library and might not
	# be available
	PosInf = float("+inf")
	NegInf = float("-inf")
import itertools
import math
import numpy
import random
import scipy
__numpy__version__ = tuple(map(int, numpy.__version__.strip().split(".")[:2]))
__scipy__version__ = tuple(map(int, scipy.__version__.strip().split(".")[:2]))
# FIXME Uncomment these lines when the interpolator problem is fixed or when we
# figure out the correct version numbers to check for
'''
if __scipy__version__ >= (0, 9) and __numpy__version__ >= (1, 7):
	from scipy.interpolate import interp1d, interp2d, LinearNDInterpolator
else:
	# pre scipy/numpy 0.9/1.7 had busted/missing interpolation code.
	# replacements are provided below
	pass
'''
from scipy.signal import signaltools

from ligo import segments

from glue import iterutils
from glue.ligolw import ligolw
from glue.ligolw import array as ligolw_array
from glue.ligolw import param as ligolw_param
from glue.ligolw import types as ligolw_types
import lal
from . import git_version


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                     Bins
#
# =============================================================================
#


class Bins(object):
	"""
	Parent class for 1-dimensional binnings.  This class is not
	intended to be used directly, but to be subclassed for use in real
	bins classes.
	"""
	def __init__(self):
		"""
		Initialize a Bins instance.  Subclasses must override this
		method.
		"""
		raise NotImplementedError

	def __len__(self):
		"""
		The number of bins in the binning.  Subclasses must
		override this method.
		"""
		raise NotImplementedError

	def __eq__(self, other):
		"""
		Two binnings are the same if they are instances of the same
		class, and describe the same binnings.  Subclasses should
		override this method but need not.
		"""
		raise NotImplementedError

	def __ne__(self, other):
		# this can be removed when we require Python 3, since in
		# Python 3 this relationship between the operators is the
		# default.
		return not self.__eq__(other)

	def __getitem__(self, x):
		"""
		Convert a co-ordinate to a bin index.  The co-ordinate can
		be a single value, or a Python slice instance describing a
		range of values.  If a single value is given, it is mapped
		to the bin index corresponding to that value.  If a slice
		is given, it is converted to a slice whose lower bound is
		the index of the bin in which the slice's lower bound
		falls, and whose upper bound is 1 greater than the index of
		the bin in which the slice's upper bound falls.  Steps are
		not supported in slices.

		Subclasses must override this method, but may chain to this
		to handle slices:

		def __getitem__(self, x):
			if type(x) is slice:
				return super(type(self), self).__getitem__(x)
			# now handle non-slices ...
		"""
		# assumes x is a slice.  works with anything that defines
		# .start, .stop and .step (but .step must be None).
		if x.step is not None and x.step != 1:
			raise NotImplementedError("step not supported: %s" % repr(x))
		return slice(self[x.start] if x.start is not None else 0, self[x.stop] + 1 if x.stop is not None else len(self))

	def __iter__(self):
		"""
		If __iter__ does not exist, Python uses __getitem__ with
		range(0) as input to define iteration. This is nonsensical
		for bin objects, so explicitly unsupport iteration.
		Subclasses do not need to override this method.
		"""
		raise NotImplementedError

	def lower(self):
		"""
		Return an array containing the locations of the lower
		boundaries of the bins.  Subclasses should override this
		method.
		"""
		raise NotImplementedError

	def centres(self):
		"""
		Return an array containing the locations of the bin
		centres.  Subclasses should override this method.
		"""
		raise NotImplementedError

	def upper(self):
		"""
		Return an array containing the locations of the upper
		boundaries of the bins.  Subclasses should override this
		method.
		"""
		raise NotImplementedError

	#
	# Sample from the binning
	#

	def randcoord(self, n = 1., domain = slice(None, None)):
		"""
		Generator yielding a sequence of x, ln(P(x)) tuples where x
		is a randomly-chosen co-ordinate and P(x) is the PDF from
		which x has been drawn evaluated at x.  Each co-ordinate is
		drawn uniformly from within a bin, which has been drawn
		from a distribution whose CDF goes as [bin index]^{n}.  For
		more information on how bins are drawn, see
		glue.iterutils.randindex.

		If a domain is given, the values returned fall within
		[start, stop].  If start or stop is None, the corresponding
		end of the binning is used.  If start or stop does not
		correspond exactly to a bin boundary, the probability of
		drawing a value from that bin is unchanged, but the values
		drawn from that bin will be restricted to the allowed part
		of the bin (the PDF is adjusted to reflect this).  No
		values are returned from bins with infinite size (after
		clipping them to the requested domain), and the PDF is
		adjusted to reflect this.

		Example:

		>>> import math
		>>> # natural log of 1/10
		>>> print("%.15g" % math.log(1./10))
		-2.30258509299405
		>>> # linear bins spanning [0, 10]
		>>> x = LinearBins(0, 10, 5).randcoord().next
		>>> # draw a random value, ln P(value) = ln 1/10
		>>> x()	# doctest: +ELLIPSIS
		(..., -2.3025850929940455)
		>>> # binning with infinite boundaries
		>>> x = ATanBins(-1, +1, 4)
		>>> # will ask for values in [0.5, +inf], i.e. the last two
		>>> # bins, but values from final bin will be disallowed, so
		>>> # return values will be uniform in part of the second
		>>> # last bin, [0.5, 0.6366]
		>>> print("%.15g" % math.log(1. / (x.upper()[-2] - 0.5)))
		1.99055359585182
		>>> x = x.randcoord(domain = slice(0.5, None)).next
		>>> x() # doctest: +ELLIPSIS
		(..., 1.9905535958518226)
		>>> # things that aren't supported:
		>>> # domain slice with a step
		>>> LinearBins(0, 10, 1).randcoord(domain = slice(None, None, 2)).next()
		Traceback (most recent call last):
			...
		NotImplementedError: step not supported: slice(None, None, 2)

		Subclasses should not override this method.
		"""
		if len(self) < 1:
			raise ValueError("empty binning")
		if domain.step is not None:
			raise NotImplementedError("step not supported: %s" % repr(domain))
		# avoid symbol look-ups in the sampling loop
		isinf = math.isinf
		uniform = random.uniform
		# determine boundaries and index range
		l = self.lower()
		u = self.upper()
		lo, hi, _ = self[domain].indices(len(l))
		if domain.start is not None:
			assert l[lo] <= domain.start
			l[lo] = domain.start
		if domain.stop is not None:
			assert u[hi - 1] >= domain.stop
			u[hi - 1] = domain.stop
		if isinf(u[lo] - l[lo]):
			lo += 1
		if isinf(u[hi - 1] - l[hi - 1]):
			hi -= 1
		if not lo < hi:
			raise ValueError("slice too small")
		# log() implicitly checks that the boundary adjustments
		# above haven't made any bins <= 0 in size.  converting
		# everything to tuples makes the sampling loop faster
		ln_dx = tuple(numpy.log(u - l))
		l = tuple(l)
		u = tuple(u)
		# one last safety check
		if any(map(isinf, ln_dx[lo:hi])):
			raise ValueError("unavoidable infinite bin detected")
		# generate samples
		for i, ln_Pi in iterutils.randindex(lo, hi, n = n):
			yield uniform(l[i], u[i]), ln_Pi - ln_dx[i]

	#
	# XML I/O related methods and data
	#

	@staticmethod
	def xml_bins_name_enc(name, suffix = u"pylal_rate_bins"):
		"""
		For internal use by XML I/O code.  Subclasses should not
		override this method.
		"""
		return u"%s:%s" % (name, suffix)

	@staticmethod
	def xml_bins_name_dec(name, suffix = u"pylal_rate_bins"):
		"""
		For internal use by XML I/O code.  Subclasses should not
		override this method.
		"""
		name = name.rsplit(u":", 1)
		if name[-1] != suffix:
			raise ValueError(name)
		return name[0]

	@classmethod
	def xml_bins_check(cls, elem, name):
		"""
		For internal use by XML I/O code.  Subclasses should not
		override this method.
		"""
		return elem.tagName == ligolw.Param.tagName and elem.hasAttribute(u"Name") and name == cls.xml_bins_name_dec(elem.Name)

	def to_xml(self):
		"""
		Construct a LIGO Light Weight XML representation of the
		Bins instance.  Subclasses must override this method to be
		serializable to LIGO Light Weight XML, otherwise they need
		not override it.
		"""
		raise NotImplementedError

	@classmethod
	def from_xml(cls, xml):
		"""
		From the XML Param element at xml, return the Bins object
		it describes.  Subclasses must override this method to be
		de-serializable from LIGO Light Weight XML, otherwise they
		need not override it.
		"""
		raise NotImplementedError


class LoHiCountBins(Bins):
	"""
	Base class to help implement binnings that can be defined by a
	lower bound, an upper bound, and a count of bins.  This is not a
	binning.
	"""
	def __init__(self, min, max, n):
		"""
		The three arguments are the minimum and maximum of the
		values spanned by the bins, and the number of bins to place
		between them.
		"""
		if not isinstance(n, int):
			raise TypeError(n)
		if n < 1:
			raise ValueError(n)
		if max <= min:
			raise ValueError((min, max))
		self.min = min
		self.max = max
		self.n = n

	def __len__(self):
		return self.n

	def __eq__(self, other):
		"""
		Two binnings are the same if they are instances of the same
		class, have the same lower and upper bounds, and the same
		count of bins.
		"""
		return isinstance(other, type(self)) and (self.min, self.max, self.n) == (other.min, other.max, other.n)

	#
	# XML I/O related methods and data
	#

	def to_xml(self):
		"""
		Construct a LIGO Light Weight XML representation of the
		Bins instance.  Subclasses must define the .xml_bins_name
		class attribute.
		"""
		return ligolw_param.Param.from_pyvalue(self.xml_bins_name_enc(self.xml_bins_name), u"%s,%s,%s" % (ligolw_types.FormatFunc[u"real_8"](self.min), ligolw_types.FormatFunc[u"real_8"](self.max), ligolw_types.FormatFunc[u"int_8s"](self.n)))

	@classmethod
	def from_xml(cls, xml):
		"""
		From the XML Param element at xml, return the Bins object
		it describes.  Subclasses must define the .xml_bins_name
		class attribute.
		"""
		if not cls.xml_bins_check(xml, cls.xml_bins_name):
			raise ValueError("not a %s" % repr(cls))
		lo, hi, n = xml.pcdata.split(u",")
		lo = ligolw_types.ToPyType[u"real_8"](lo)
		hi = ligolw_types.ToPyType[u"real_8"](hi)
		n = ligolw_types.ToPyType[u"int_8s"](n)
		return cls(lo, hi, n)


class IrregularBins(Bins):
	"""
	Bins with arbitrary, irregular spacing.  We only require strict
	monotonicity of the bin boundaries.  N boundaries define N-1 bins.

	Example:

	>>> x = IrregularBins([0.0, 11.0, 15.0, numpy.inf])
	>>> len(x)
	3
	>>> x[1]
	0
	>>> x[1.5]
	0
	>>> x[13]
	1
	>>> x[25]
	2
	>>> x[4:17]
	slice(0, 3, None)
	>>> IrregularBins([0.0, 15.0, 11.0])
	Traceback (most recent call last):
		...
	ValueError: non-monotonic boundaries provided
	>>> y = IrregularBins([0.0, 11.0, 15.0, numpy.inf])
	>>> x == y
	True
	>>> import sys
	>>> x.to_xml().write(sys.stdout)	# doctest: +NORMALIZE_WHITESPACE
	<Param Type="lstring" Name="irregularbins:pylal_rate_bins:param">0,11,15,inf</Param>
	>>> IrregularBins.from_xml(x.to_xml()) == x
	True
	"""
	def __init__(self, boundaries):
		"""
		Initialize a set of custom bins with the bin boundaries.
		This includes all left edges plus the right edge.  The
		boundaries must be monotonic and there must be at least two
		elements.
		"""
		# check pre-conditions
		if len(boundaries) < 2:
			raise ValueError("less than two boundaries provided")
		self.boundaries = numpy.array(boundaries)
		if (self.boundaries[:-1] > self.boundaries[1:]).any():
			raise ValueError("non-monotonic boundaries provided")
		self.lo, self.hi = float(self.boundaries[0]), float(self.boundaries[-1])

	def __eq__(self, other):
		"""
		Two binnings are the same if they are instances of the same
		class, and have the same boundaries.
		"""
		return isinstance(other, type(self)) and (self.boundaries == other.boundaries).all()

	def __len__(self):
		return len(self.boundaries) - 1

	def __getitem__(self, x):
		# slice cannot be sub-classed so no need to use
		# isinstance()
		if type(x) is slice:
			return super(IrregularBins, self).__getitem__(x)
		if self.lo <= x < self.hi:
			return self.boundaries.searchsorted(x, side = "right") - 1
		# special measure-zero edge case
		if x == self.hi:
			return len(self.boundaries) - 2
		raise IndexError(x)

	def lower(self):
		return self.boundaries[:-1]

	def upper(self):
		return self.boundaries[1:]

	def centres(self):
		return (self.lower() + self.upper()) / 2.0

	#
	# XML I/O related methods and data
	#

	xml_bins_name = u"irregularbins"

	def to_xml(self):
		"""
		Construct a LIGO Light Weight XML representation of the
		Bins instance.
		"""
		return ligolw_param.Param.from_pyvalue(self.xml_bins_name_enc(self.xml_bins_name), u",".join(map(ligolw_types.FormatFunc[u"real_8"], self.boundaries)))

	@classmethod
	def from_xml(cls, xml):
		"""
		From the XML Param element at xml, return the Bins object
		it describes.
		"""
		if not cls.xml_bins_check(xml, cls.xml_bins_name):
			raise ValueError("not a %s" % repr(cls))
		return cls(map(ligolw_types.ToPyType[u"real_8"], xml.pcdata.split(u",")))


class LinearBins(LoHiCountBins):
	"""
	Linearly-spaced bins.  There are n bins of equal size, the first
	bin starts on the lower bound and the last bin ends on the upper
	bound inclusively.

	Example:

	>>> x = LinearBins(1.0, 25.0, 3)
	>>> x.lower()
	array([  1.,   9.,  17.])
	>>> x.upper()
	array([  9.,  17.,  25.])
	>>> x.centres()
	array([  5.,  13.,  21.])
	>>> x[1]
	0
	>>> x[1.5]
	0
	>>> x[10]
	1
	>>> x[25]
	2
	>>> x[0:27]
	Traceback (most recent call last):
		...
	IndexError: 0
	>>> x[1:25]
	slice(0, 3, None)
	>>> x[:25]
	slice(0, 3, None)
	>>> x[10:16.9]
	slice(1, 2, None)
	>>> x[10:17]
	slice(1, 3, None)
	>>> x[10:]
	slice(1, 3, None)
	>>> import sys
	>>> x.to_xml().write(sys.stdout)	# doctest: +NORMALIZE_WHITESPACE
	<Param Type="lstring" Name="linbins:pylal_rate_bins:param">1,25,3</Param>
	>>> LinearBins.from_xml(x.to_xml()) == x
	True
	"""
	def __init__(self, min, max, n):
		super(LinearBins, self).__init__(min, max, n)
		self.delta = float(max - min) / n

	def __getitem__(self, x):
		# slice cannot be sub-classed so no need to use
		# isinstance()
		if type(x) is slice:
			return super(LinearBins, self).__getitem__(x)
		if self.min <= x < self.max:
			return int(math.floor((x - self.min) / self.delta))
		if x == self.max:
			# special "measure zero" corner case
			return len(self) - 1
		raise IndexError(x)

	def lower(self):
		return numpy.linspace(self.min, self.max - self.delta, len(self))

	def centres(self):
		return numpy.linspace(self.min + self.delta / 2., self.max - self.delta / 2., len(self))

	def upper(self):
		return numpy.linspace(self.min + self.delta, self.max, len(self))

	#
	# XML I/O related methods and data
	#

	xml_bins_name = u"linbins"


class LinearPlusOverflowBins(LoHiCountBins):
	"""
	Linearly-spaced bins with overflow at the edges.  There are n-2
	bins of equal size.  The bin 1 starts on the lower bound and bin
	n-2 ends on the upper bound.  Bins 0 and n-1 are overflow going
	from -infinity to the lower bound and from the upper bound to
	+infinity respectively.  Must have n >= 3.

	Example:

	>>> x = LinearPlusOverflowBins(1.0, 25.0, 5)
	>>> x.centres()
	array([-inf,   5.,  13.,  21.,  inf])
	>>> x.lower()
	array([-inf,   1.,   9.,  17.,  25.])
	>>> x.upper()
	array([  1.,   9.,  17.,  25.,  inf])
	>>> x[float("-inf")]
	0
	>>> x[0]
	0
	>>> x[1]
	1
	>>> x[10]
	2
	>>> x[24.99999999]
	3
	>>> x[25]
	4
	>>> x[100]
	4
	>>> x[float("+inf")]
	4
	>>> x[float("-inf"):9]
	slice(0, 3, None)
	>>> x[9:float("+inf")]
	slice(2, 5, None)
	>>> import sys
	>>> x.to_xml().write(sys.stdout)	# doctest: +NORMALIZE_WHITESPACE
	<Param Type="lstring" Name="linplusoverflowbins:pylal_rate_bins:param">1,25,5</Param>
	>>> LinearPlusOverflowBins.from_xml(x.to_xml()) == x
	True
	"""
	def __init__(self, min, max, n):
		if n < 3:
			raise ValueError("n must be >= 3")
		super(LinearPlusOverflowBins, self).__init__(min, max, n)
		self.delta = float(max - min) / (n - 2)

	def __getitem__(self, x):
		# slice cannot be sub-classed so no need to use
		# isinstance()
		if type(x) is slice:
			return super(LinearPlusOverflowBins, self).__getitem__(x)
		if self.min <= x < self.max:
			return int(math.floor((x - self.min) / self.delta)) + 1
		if x >= self.max:
			# +infinity overflow bin
			return len(self) - 1
		if x < self.min:
			# -infinity overflow bin
			return 0
		raise IndexError(x)

	def lower(self):
		return numpy.concatenate((numpy.array([NegInf]), numpy.linspace(self.min, self.max - self.delta, len(self) - 1)))

	def centres(self):
		return numpy.concatenate((numpy.array([NegInf]), numpy.linspace(self.min + self.delta / 2., self.max - self.delta / 2., len(self) - 2), numpy.array([PosInf])))

	def upper(self):
		return numpy.concatenate((numpy.linspace(self.min + self.delta, self.max, len(self) - 1), numpy.array([PosInf])))

	#
	# XML I/O related methods and data
	#

	xml_bins_name = u"linplusoverflowbins"


class LogarithmicBins(LoHiCountBins):
	"""
	Logarithmically-spaced bins.  There are n bins, each of whose upper
	and lower bounds differ by the same factor.  The first bin starts
	on the lower bound, and the last bin ends on the upper bound
	inclusively.

	Example:

	>>> x = LogarithmicBins(1.0, 25.0, 3)
	>>> x[1]
	0
	>>> x[5]
	1
	>>> x[25]
	2
	>>> import sys
	>>> x.to_xml().write(sys.stdout)	# doctest: +NORMALIZE_WHITESPACE
	<Param Type="lstring" Name="logbins:pylal_rate_bins:param">1,25,3</Param>
	>>> LogarithmicBins.from_xml(x.to_xml()) == x
	True
	"""
	def __init__(self, min, max, n):
		super(LogarithmicBins, self).__init__(min, max, n)
		self.logmin = math.log(min)
		self.logmax = math.log(max)
		self.delta = (self.logmax - self.logmin) / n

	def __getitem__(self, x):
		# slice cannot be sub-classed so no need to use
		# isinstance()
		if type(x) is slice:
			return super(LogarithmicBins, self).__getitem__(x)
		if self.min <= x < self.max:
			return int(math.floor((math.log(x) - self.logmin) / self.delta))
		if x == self.max:
			# special "measure zero" corner case
			return len(self) - 1
		raise IndexError(x)

	def lower(self):
		return numpy.exp(numpy.linspace(self.logmin, self.logmax - self.delta, len(self)))

	def centres(self):
		return numpy.exp(numpy.linspace(self.logmin, self.logmax - self.delta, len(self)) + self.delta / 2.)

	def upper(self):
		return numpy.exp(numpy.linspace(self.logmin + self.delta, self.logmax, len(self)))

	#
	# XML I/O related methods and data
	#

	xml_bins_name = u"logbins"


class LogarithmicPlusOverflowBins(LoHiCountBins):
	"""
	Logarithmically-spaced bins plus one bin at each end that goes to
	zero and positive infinity respectively.  There are n-2 bins each
	of whose upper and lower bounds differ by the same factor.  Bin 1
	starts on the lower bound, and bin n-2 ends on the upper bound
	inclusively.  Bins 0 and n-1 are overflow bins extending from 0 to
	the lower bound and from the upper bound to +infinity respectively.
	Must have n >= 3.

	Example:

	>>> x = LogarithmicPlusOverflowBins(1.0, 25.0, 5)
	>>> x[0]
	0
	>>> x[1]
	1
	>>> x[5]
	2
	>>> x[24.999]
	3
	>>> x[25]
	4
	>>> x[100]
	4
	>>> x.lower()
	array([  0.        ,   1.        ,   2.92401774,   8.54987973,  25.        ])
	>>> x.upper()
	array([  1.        ,   2.92401774,   8.54987973,  25.        ,          inf])
	>>> x.centres()
	array([  0.        ,   1.70997595,   5.        ,  14.62008869,          inf])
	>>> import sys
	>>> x.to_xml().write(sys.stdout)	# doctest: +NORMALIZE_WHITESPACE
	<Param Type="lstring" Name="logplusoverflowbins:pylal_rate_bins:param">1,25,5</Param>
	>>> LogarithmicPlusOverflowBins.from_xml(x.to_xml()) == x
	True
	"""
	def __init__(self, min, max, n):
		if n < 3:
			raise ValueError("n must be >= 3")
		super(LogarithmicPlusOverflowBins, self).__init__(min, max, n)
		self.logmin = math.log(min)
		self.logmax = math.log(max)
		self.delta = (self.logmax - self.logmin) / (n - 2)

	def __getitem__(self, x):
		# slice cannot be sub-classed so no need to use
		# isinstance()
		if type(x) is slice:
			return super(LogarithmicPlusOverflowBins, self).__getitem__(x)
		if self.min <= x < self.max:
			return 1 + int(math.floor((math.log(x) - self.logmin) / self.delta))
		if x >= self.max:
			# infinity overflow bin
			return len(self) - 1
		if x < self.min:
			# zero overflow bin
			return 0
		raise IndexError(x)

	def lower(self):
		return numpy.concatenate((numpy.array([0.]), numpy.exp(numpy.linspace(self.logmin, self.logmax, len(self) - 1))))

	def centres(self):
		return numpy.concatenate((numpy.array([0.]), numpy.exp(numpy.linspace(self.logmin, self.logmax - self.delta, len(self) - 2) + self.delta / 2.), numpy.array([PosInf])))

	def upper(self):
		return numpy.concatenate((numpy.exp(numpy.linspace(self.logmin, self.logmax, len(self) - 1)), numpy.array([PosInf])))

	#
	# XML I/O related methods and data
	#

	xml_bins_name = u"logplusoverflowbins"


class ATanBins(LoHiCountBins):
	"""
	Bins spaced uniformly in tan^-1 x.  Provides approximately linear
	binning in the middle portion, with the bin density dropping
	asymptotically to 0 as x goes to +/- \\infty.  The min and max
	parameters set the bounds of the region of approximately
	uniformly-spaced bins.  In a sense, these are where the roll-over
	from uniformly-spaced bins to asymptotically diminishing bin
	density occurs.  There is a total of n bins.

	Example:

	>>> x = ATanBins(-1.0, +1.0, 11)
	>>> x[float("-inf")]
	0
	>>> x[0]
	5
	>>> x[float("+inf")]
	10
	>>> x.centres()
	array([-4.42778777, -1.39400285, -0.73469838, -0.40913068, -0.18692843,
	        0.        ,  0.18692843,  0.40913068,  0.73469838,  1.39400285,
                4.42778777])
	>>> import sys
	>>> x.to_xml().write(sys.stdout)	# doctest: +NORMALIZE_WHITESPACE
	<Param Type="lstring" Name="atanbins:pylal_rate_bins:param">-1,1,11</Param>
	>>> ATanBins.from_xml(x.to_xml()) == x
	True
	"""
	def __init__(self, min, max, n):
		super(ATanBins, self).__init__(min, max, n)
		self.mid = (min + max) / 2.0
		self.scale = math.pi / float(max - min)
		self.delta = 1.0 / n

	def __getitem__(self, x):
		# slice cannot be sub-classed so no need to use
		# isinstance()
		if type(x) is slice:
			return super(ATanBins, self).__getitem__(x)
		# map to the domain [0, 1]
		x = math.atan(float(x - self.mid) * self.scale) / math.pi + 0.5
		if x < 1.:
			return int(math.floor(x / self.delta))
		# x == 1, special "measure zero" corner case
		return len(self) - 1

	def lower(self):
		x = numpy.tan(numpy.linspace(-math.pi / 2., +math.pi / 2., len(self), endpoint = False)) / self.scale + self.mid
		x[0] = NegInf
		return x

	def centres(self):
		offset = 0.5 * math.pi * self.delta
		return numpy.tan(numpy.linspace(-math.pi / 2. + offset, +math.pi / 2. + offset, len(self), endpoint = False)) / self.scale + self.mid

	def upper(self):
		offset = math.pi * self.delta
		x = numpy.tan(numpy.linspace(-math.pi / 2. + offset, +math.pi / 2. + offset, len(self), endpoint = False)) / self.scale + self.mid
		x[-1] = PosInf
		return x

	#
	# XML I/O related methods and data
	#

	xml_bins_name = u"atanbins"


class ATanLogarithmicBins(LoHiCountBins, IrregularBins):
	"""
	Provides the same binning as the ATanBins class but in the
	logarithm of the variable.  The min and max parameters set the
	bounds of the interval of approximately logarithmically-spaced
	bins.  In a sense, these are where the roll-over from
	logarithmically-spaced bins to asymptotically diminishing bin
	density occurs.

	Example:

	>>> x = ATanLogarithmicBins(+1.0, +1000.0, 11)
	>>> x[0]
	0
	>>> x[30]
	5
	>>> x[float("+inf")]
	10
	>>> x.centres()
	array([  7.21636246e-06,   2.56445876e-01,   2.50007148e+00,
		 7.69668960e+00,   1.65808715e+01,   3.16227766e+01,
		 6.03104608e+01,   1.29925988e+02,   3.99988563e+02,
		 3.89945831e+03,   1.38573971e+08])
	>>> import sys
	>>> x.to_xml().write(sys.stdout)	# doctest: +NORMALIZE_WHITESPACE
	<Param Type="lstring" Name="atanlogbins:pylal_rate_bins:param">1,1000,11</Param>
	>>> ATanLogarithmicBins.from_xml(x.to_xml()) == x
	True

	It is relatively easy to choose limits and a count of bins that
	result in numerical overflows and underflows when computing bin
	boundaries.  When this happens, one or more bins at the ends of the
	binning ends up with identical upper and lower boundaries (either
	0, or +inf), and this class behaves as though those bins simply
	don't exist.  That is, the actual number of bins can be less than
	the number requested.  len() returns the actual number of bins ---
	how large an array the binning corresponds to.
	"""
	def __init__(self, min, max, n):
		if not isinstance(n, int):
			raise TypeError(n)
		if n < 1:
			raise ValueError(n)
		if max <= min:
			raise ValueError((min, max))
		self.mid = (math.log(min) + math.log(max)) / 2.0
		self.scale = math.pi / (math.log(max) - math.log(min))
		self.delta = 1.0 / n
		boundaries = numpy.tan(-math.pi / 2 + math.pi * self.delta * numpy.arange(n)) / self.scale + self.mid
		with numpy.errstate(over = "ignore"):
			boundaries = numpy.exp(boundaries)
		boundaries = numpy.hstack((boundaries, [PosInf, 0.]))
		keepers = boundaries[:-1] != boundaries[1:]
		IrregularBins.__init__(self, boundaries[:-1][keepers])
		self.keepers = keepers[:-1]
		self.min = min
		self.max = max
		self.n = n

	__len__ = IrregularBins.__len__

	def centres(self):
		offset = 0.5 * math.pi * self.delta
		centres = numpy.tan(numpy.linspace(-math.pi / 2. + offset, +math.pi / 2. + offset, self.n, endpoint = False)) / self.scale + self.mid
		with numpy.errstate(over = "ignore"):
			return numpy.exp(centres)[self.keepers]

	#
	# XML I/O related methods and data
	#

	xml_bins_name = u"atanlogbins"


class Categories(Bins):
	"""
	Categories is a many-to-one mapping from a value to an integer
	category index.  A value belongs to a category if it is contained
	in the category's defining collection.  If a value is contained in
	more than one category's defining collection, it belongs to the
	category with the smallest index.  IndexError is raised if a value
	is not contained in any category's defining collection.

	Example with discrete values:

	>>> categories = Categories([
	...	set((frozenset(("H1", "L1")), frozenset(("H1", "V1")))),
	...	set((frozenset(("H1", "L1", "V1")),))
	... ])
	>>> categories[set(("H1", "L1"))]
	0
	>>> categories[set(("H1", "V1"))]
	0
	>>> categories[set(("H1", "L1", "V1"))]
	1

	Example with continuous values:

	>>> from ligo.segments import *
	>>> categories = Categories([
	...	segmentlist([segment(1, 3), segment(5, 7)]),
	...	segmentlist([segment(0, PosInfinity)])
	... ])
	>>> categories[2]
	0
	>>> categories[4]
	1
	>>> categories[-1]
	Traceback (most recent call last):
		...
	IndexError: -1

	This last example demonstrates the behaviour when the intersection
	of the categories is not the empty set.
	"""
	def __init__(self, categories):
		"""
		categories is an iterable of containers (objects that
		support the "in" operator) defining the categories.
		Objects will be mapped to the integer index of the
		container that contains them.
		"""
		# make immutable copy
		self.containers = tuple(categories)

	def __len__(self):
		return len(self.containers)

	def __getitem__(self, value):
		"""
		Return i if value is contained in i-th container.  If value
		is not contained in any of the containers, raise an
		IndexError.  This is O(n).
		"""
		for i, s in enumerate(self.containers):
			if value in s:
				return i
		raise IndexError(value)

	def __eq__(self, other):
		return isinstance(other, type(self)) and self.containers == other.containers

	def centres(self):
		return self.containers

	#
	# XML I/O related methods and data
	#

	xml_bins_name = u"categorybins"

	def to_xml(self):
		"""
		Construct a LIGO Light Weight XML representation of the
		Bins instance.
		"""
		return ligolw_param.Param.build(self.xml_bins_name_enc(self.xml_bins_name), u"yaml", self.containers)

	@classmethod
	def from_xml(cls, xml):
		"""
		From the XML Param element at xml, return the Bins object
		it describes.
		"""
		if not cls.xml_bins_check(xml, cls.xml_bins_name):
			raise ValueError("not a %s" % repr(cls))
		return cls(xml.pcdata)


class HashableBins(Categories):
	"""
	Maps hashable objects (things that can be used as dictionary keys)
	to integers.

	Example:
	>>> x = HashableBins([
	...    frozenset(("H1", "L1")),
	...    frozenset(("H1", "V1")),
	...    frozenset(("L1", "V1")),
	...    frozenset(("H1", "L1", "V1"))
	... ])
	>>> x[frozenset(("H1", "L1"))]
	0
	>>> x[set(("H1", "L1"))]	# equal, but not hashable
	Traceback (most recent call last):
		...
	IndexError: set(['H1', 'L1'])
	>>> x.centres()[2]
	frozenset(['V1', 'L1'])
	"""
	def __init__(self, hashables):
		super(HashableBins, self).__init__(hashables)
		self.mapping = dict(zip(self.containers, range(len(self.containers))))

	def __getitem__(self, value):
		try:
			return self.mapping[value]
		except (KeyError, TypeError):
			raise IndexError(value)

	xml_bins_name = u"hashablebins"


class NDBins(tuple):
	"""
	Multi-dimensional co-ordinate binning.  An instance of this object
	is used to convert a tuple of co-ordinates into a tuple of bin
	indices.  This can be used to allow the contents of an array object
	to be accessed with real-valued coordinates.

	NDBins is a subclass of the tuple builtin, and is initialized with
	an iterable of instances of subclasses of Bins.  Each Bins subclass
	instance describes the binning to apply in the corresponding
	co-ordinate direction, and the number of them sets the dimensions
	of the binning.

	Example:

	>>> x = NDBins((LinearBins(1, 25, 3), LogarithmicBins(1, 25, 3)))
	>>> x[1, 1]
	(0, 0)
	>>> x[1.5, 1]
	(0, 0)
	>>> x[10, 1]
	(1, 0)
	>>> x[1, 5]
	(0, 1)
	>>> x[1, 1:5]
	(0, slice(0, 2, None))
	>>> x.centres()
	(array([  5.,  13.,  21.]), array([  1.70997595,   5.        ,  14.62008869]))
	>>> y = NDBins((LinearBins(1, 25, 3), LogarithmicBins(1, 25, 3)))
	>>> x == y
	True
	>>> y = NDBins((LinearBins(1, 25, 4), LogarithmicBins(1, 25, 3)))
	>>> x == y
	False
	>>> y = NDBins((LogarithmicBins(1, 25, 3), LogarithmicBins(1, 25, 3)))
	>>> x == y
	False
	>>> from glue.ligolw.ligolw import LIGO_LW
	>>> import sys
	>>> elem = x.to_xml(LIGO_LW())
	>>> elem.write(sys.stdout)	# doctest: +NORMALIZE_WHITESPACE
	<LIGO_LW>
		<Param Type="lstring" Name="linbins:pylal_rate_bins:param">1,25,3</Param>
		<Param Type="lstring" Name="logbins:pylal_rate_bins:param">1,25,3</Param>
	</LIGO_LW>
	>>> NDBins.from_xml(elem) == x
	True

	Note that the co-ordinates to be converted must be a tuple, even if
	it is only a 1-dimensional co-ordinate.
	"""
	def __init__(self, binnings):
		self._getitems = tuple(binning.__getitem__ for binning in binnings)
		# instances cannot define a .__call__() attribute to make
		# themselves callable, Python always looks up .__call__()
		# on the class.  so we define .__realcall__() here and then
		# have .__call__() chain to it.  Python3 does not transfer
		# the current variable scope into exec() so we have to do
		# it for it ... whatever.
		define__realcall__ = """def __realcall__(self, %s):
	_getitems = self._getitems
	return %s""" % (", ".join("x%d" % i for i in range(len(binnings))), ", ".join("_getitems[%d](x%d)" % (i, i) for i in range(len(binnings))))
		l = {}
		exec(define__realcall__, globals(), l)
		__realcall__ = l["__realcall__"]
		self.__realcall__ = __realcall__.__get__(self)

	def __getitem__(self, coords):
		"""
		When coords is a tuple this is a synonym for self(*coords),
		otherwise coords is interpreted as an index into ourselves
		and the Bins object at that index is returned

		Example:

		>>> x = NDBins((LinearBins(1, 25, 3), LogarithmicBins(1, 25, 3)))
		>>> x[1, 1]
		(0, 0)
		>>> # slices can be given syntactically
		>>> x[10:12, 1]
		(slice(1, 2, None), 0)
		>>> type(x[1])
		<class 'lal.rate.LogarithmicBins'>

		Note that if the argument is to be interpreted as a
		co-ordinate it must be a tuple even if it is only a
		1-dimensional co-ordinate.

		Example:

		>>> x = NDBins((LinearBins(1, 25, 3),))
		>>> x[1,]
		(0,)
		"""
		return self(*coords) if isinstance(coords, tuple) else tuple.__getitem__(self, coords)

	def __call__(self, *coords):
		"""
		Convert an N-dimensional co-ordinate to an N-tuple of bin
		indices using the Bins instances in this object.  Calling
		the NDBins instance instead of using indexing (the "[]"
		operator) provides a more direct, faster, interface to the
		Bins instances contained herein, but slices cannot be given
		syntactically in the argument list

		Example:

		>>> x = NDBins((LinearBins(1, 25, 3), LogarithmicBins(1, 25, 3)))
		>>> x[1, 1]	# access using index operator
		(0, 0)
		>>> x(1, 1)	# access by calling
		(0, 0)
		>>> x = NDBins((LinearBins(1, 25, 3),))
		>>> x(1)	# explicit tuples not required for 1D
		(0,)
		>>> x = NDBins((LinearBins(1, 25, 1000),))
		>>> x(slice(10, 12))	# slices (constructed manually)
		(slice(375, 459, None),)
		>>> x = NDBins((Categories([set(("Cow", "Chicken", "Goat")), set(("Tractor", "Plough")), set(("Barn", "House"))]),))
		>>> x("Cow")
		(0,)

		Each co-ordinate can be anything the corresponding Bins
		instance will accept.
		"""
		return self.__realcall__(*coords)

	@property
	def shape(self):
		"""
		Tuple of the number of bins along each dimension.

		Example:

		>>> NDBins((LinearBins(0, 6, 3), LinearBins(0, 10, 5))).shape
		(3, 5)
		"""
		return tuple(len(b) for b in self)

	def lower(self):
		"""
		Return a tuple of arrays, where each array contains the
		locations of the lower boundaries of the bins in the
		corresponding dimension.

		Example:

		>>> NDBins((LinearBins(0, 6, 3), LinearBins(0, 10, 5))).lower()
		(array([ 0.,  2.,  4.]), array([ 0.,  2.,  4.,  6.,  8.]))
		"""
		return tuple(b.lower() for b in self)

	def centres(self):
		"""
		Return a tuple of arrays, where each array contains the
		locations of the bin centres for the corresponding
		dimension.

		Example:

		>>> NDBins((LinearBins(0, 6, 3), LinearBins(0, 10, 5))).centres()
		(array([ 1.,  3.,  5.]), array([ 1.,  3.,  5.,  7.,  9.]))
		"""
		return tuple(b.centres() for b in self)

	def upper(self):
		"""
		Return a tuple of arrays, where each array contains the
		locations of the upper boundaries of the bins in the
		corresponding dimension.

		Example:

		>>> NDBins((LinearBins(0, 6, 3), LinearBins(0, 10, 5))).upper()
		(array([ 2.,  4.,  6.]), array([  2.,   4.,   6.,   8.,  10.]))
		"""
		return tuple(b.upper() for b in self)

	def volumes(self):
		"""
		Return an n-dimensional array of the bin volumes.

		Example:

		>>> # 3x5 grid of bins, each 2 units by 2 units
		>>> x = NDBins((LinearBins(0, 6, 3), LinearBins(0, 10, 5)))
		>>> x.volumes()
		array([[ 4.,  4.,  4.,  4.,  4.],
		       [ 4.,  4.,  4.,  4.,  4.],
		       [ 4.,  4.,  4.,  4.,  4.]])
		"""
		volumes = tuple(u - l for u, l in zip(self.upper(), self.lower()))
		if len(volumes) == 1:
			# 1D short-cut
			return volumes[0]
		try:
			return numpy.einsum(",".join("abcdefghijklmnopqrstuvwxyz"[:len(volumes)]), *volumes)
		except AttributeError:
			# numpy < 1.6
			result = reduce(numpy.outer, volumes)
			result.shape = tuple(len(v) for v in volumes)
			return result

	def randcoord(self, ns = None, domain = None):
		"""
		Generator yielding a sequence of (x0, x1, ...), ln(P(x0,
		x1, ...)) tuples where (x0, x1, ...) is a randomly-chosen
		co-ordinate in the N-dimensional binning and P(x0, x1, ...)
		is the PDF from which the co-ordinate tuple has been drawn
		evaluated at those co-ordinates.  If ns is not None it must
		be a sequence of floats whose length matches the dimension
		of the binning.  The floats will set the exponents, in
		order, of the CDFs for the generators used for each
		co-ordinate.  If domain is not None it must be a sequence
		of slice objects whose length matches the dimension of the
		binning.  The slice objects will be passed, in order, as
		the domain keyword argument to the .randcoord() method
		corresponding to each dimension.  For more information on
		how each of the co-ordinates is drawn, see
		Bins.randcoord().

		Example:

		>>> binning = NDBins((LinearBins(0, 10, 5), LinearBins(0, 10, 5)))
		>>> coord = binning.randcoord().next
		>>> coord()	# doctest: +ELLIPSIS
		((..., ...), -4.6051701859880909)
		"""
		if ns is None:
			ns = (1.,) * len(self)
		if domain is None:
			domain = (slice(None, None),) * len(self)
		coordgens = tuple(iter(binning.randcoord(n, domain = d)).next for binning, n, d in zip(self, ns, domain))
		while 1:
			seq = sum((coordgen() for coordgen in coordgens), ())
			yield seq[0::2], sum(seq[1::2])

	#
	# XML I/O methods and data
	#

	xml_bins_name_mapping = dict((cls.xml_bins_name, cls) for cls in (LinearBins, LinearPlusOverflowBins, LogarithmicBins, LogarithmicPlusOverflowBins, ATanBins, ATanLogarithmicBins, Categories, HashableBins))
	xml_bins_name_mapping.update(list(zip(xml_bins_name_mapping.values(), xml_bins_name_mapping.keys())))

	def to_xml(self, elem):
		"""
		Construct a LIGO Light Weight XML representation of the
		NDBins instance.  The representation is an in-order list of
		Param elements, typically inserted inside a LIGO_LW
		element.  The elem argument provides the XML element to
		which the Param elements should be appended as children.

		NOTE:  The decoding process will require a specific parent
		element to be provided, and all Param elements that are
		immediate children of that element and contain the correct
		suffix will be used to reconstruct the NDBins.  At this
		time the suffix is undocumented, so to guarantee
		compatibility the Param elements should not be inserted
		into an element that might contain other, unrelated, Params
		among its immediate children.
		"""
		for binning in self:
			elem.appendChild(binning.to_xml())
		return elem

	@classmethod
	def from_xml(cls, xml):
		"""
		From the XML document tree rooted at xml construct an
		return an NDBins object described by the Param elements
		therein.  Note, the XML element must be the immediate
		parent of the Param elements describing the NDBins.
		"""
		params = []
		for elem in xml.childNodes:
			if elem.tagName != ligolw.Param.tagName:
				continue
			try:
				Bins.xml_bins_name_dec(elem.Name)
			except ValueError:
				continue
			params.append(elem)
		if not params:
			raise ValueError("no Param elements found at '%s'" % repr(xml))
		return cls([cls.xml_bins_name_mapping[Bins.xml_bins_name_dec(elem.Name)].from_xml(elem) for elem in params])


#
# =============================================================================
#
#                              Segments and Bins
#
# =============================================================================
#


def bins_spanned(bins, seglist):
	"""
	Input is a Bins subclass instance and a ligo.segments.segmentlist
	instance.  The output is an array object the length of the binning,
	which each element in the array set to the interval in the
	corresponding bin spanned by the segment list.

	Example:

	>>> from ligo.segments import *
	>>> s = segmentlist([segment(1.5, 10.333), segment(15.8, 24)])
	>>> b = LinearBins(0, 30, 100)
	>>> bins_spanned(b, s)
	array([ 0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.3  ,  0.3  ,  0.3  ,
	        0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,
	        0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,
	        0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,
	        0.3  ,  0.3  ,  0.133,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
	        0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
	        0.   ,  0.   ,  0.   ,  0.   ,  0.1  ,  0.3  ,  0.3  ,  0.3  ,
	        0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,
	        0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,
	        0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,  0.3  ,
	        0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
	        0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,  0.   ,
	        0.   ,  0.   ,  0.   ,  0.   ])
	"""
	lower = bins.lower()
	upper = bins.upper()
	# performance improvement:  pre-clip segments to the domain of the
	# binning
	seglist = seglist & segments.segmentlist([segments.segment(lower[0], upper[-1])])
	return numpy.fromiter((abs(seglist & segments.segmentlist([seg])) for seg in zip(lower, upper)), dtype = lower.dtype, count = len(bins))


#
# =============================================================================
#
#                                 Binned Array
#
# =============================================================================
#


class BinnedArray(object):
	"""
	A convenience wrapper, using the NDBins class to provide access to
	the elements of an array object.  Technical reasons preclude
	providing a subclass of the array object, so the array data is made
	available as the "array" attribute of this class.

	Examples:

	Note that even for 1 dimensional arrays the index must be a tuple.

	>>> x = BinnedArray(NDBins((LinearBins(0, 10, 5),)))
	>>> x.at_centres()
	array([ 0.,  0.,  0.,  0.,  0.])
	>>> x[0,] += 1
	>>> x[0.5,] += 1
	>>> x.at_centres()
	array([ 2.,  0.,  0.,  0.,  0.])
	>>> x.argmax()
	(1.0,)

	Note the relationship between the binning limits, the bin centres,
	and the co-ordinates of the BinnedArray

	>>> x = BinnedArray(NDBins((LinearBins(-0.5, 1.5, 2), LinearBins(-0.5, 1.5, 2))))
	>>> x.bins.centres()
	(array([ 0.,  1.]), array([ 0.,  1.]))
	>>> x[0, 0] = 0
	>>> x[0, 1] = 1
	>>> x[1, 0] = 2
	>>> x[1, 1] = 4
	>>> x.at_centres()
	array([[ 0.,  1.],
	       [ 2.,  4.]])
	>>> x[0, 0]
	0.0
	>>> x[0, 1]
	1.0
	>>> x[1, 0]
	2.0
	>>> x[1, 1]
	4.0
	>>> x.argmin()
	(0.0, 0.0)
	>>> x.argmax()
	(1.0, 1.0)

	A BinnedArray can be initialized from an existing array if the
	array's shape is the same as the binning's.

	>>> import numpy
	>>> x = BinnedArray(NDBins((LinearBins(0, 10, 5),)), array = numpy.zeros((5,)))
	>>> x = BinnedArray(NDBins((LinearBins(0, 10, 5),)), array = numpy.zeros((5,1)))
	Traceback (most recent call last):
		...
	ValueError: bins (shape = (5,)) and array (shape = (5, 1)), if supplied, must have the same shape

	A BinnedArray can be serialized to LIGO Light Weight XML.

	>>> import sys
	>>> x = BinnedArray(NDBins((LinearBins(-0.5, 1.5, 2), LinearBins(-0.5, 1.5, 2))))
	>>> elem = x.to_xml(u"test")
	>>> elem.write(sys.stdout)	# doctest: +NORMALIZE_WHITESPACE
	<LIGO_LW Name="test:pylal_rate_binnedarray">
		<Param Type="lstring" Name="linbins:pylal_rate_bins:param">-0.5,1.5,2</Param>
		<Param Type="lstring" Name="linbins:pylal_rate_bins:param">-0.5,1.5,2</Param>
		<Array Type="real_8" Name="array:array">
			<Dim>2</Dim>
			<Dim>2</Dim>
			<Stream Delimiter=" " Type="Local">
				0 0
				0 0
			</Stream>
		</Array>
	</LIGO_LW>
	>>> y = BinnedArray.from_xml(elem, u"test")
	>>> y.bins == x.bins
	True
	>>> (y.array == x.array).all()
	True
	"""
	def __init__(self, bins, array = None, dtype = "double"):
		self.bins = bins
		if array is None:
			self.array = numpy.zeros(bins.shape, dtype = dtype)
		elif array.shape != bins.shape:
			raise ValueError("bins (shape = %s) and array (shape = %s), if supplied, must have the same shape" % (str(bins.shape), str(array.shape)))
		else:
			self.array = array

	def __getitem__(self, coords):
		return self.array[self.bins(*coords)]

	def __setitem__(self, coords, val):
		self.array[self.bins(*coords)] = val

	def __len__(self):
		return len(self.array)

	def __iadd__(self, other):
		"""
		Add the contents of another BinnedArray object to this one.
		Both must have identical binnings.

		Example:

		>>> x = BinnedArray(NDBins((LinearBins(-0.5, 1.5, 2), LinearBins(-0.5, 1.5, 2))))
		>>> x[0, 0] = 0
		>>> x[0, 1] = 1
		>>> x[1, 0] = 2
		>>> x[1, 1] = 4
		>>> x.at_centres()
		array([[ 0.,  1.],
		       [ 2.,  4.]])
		>>> x += x
		>>> x.at_centres()
		array([[ 0.,  2.],
		       [ 4.,  8.]])
		"""
		if self.bins != other.bins:
			raise TypeError("incompatible binning: %s" % repr(other))
		self.array += other.array
		return self

	def __add__(self, other):
		"""
		Add two BinnedArray objects together.

		Example:

		>>> x = BinnedArray(NDBins((LinearBins(-0.5, 1.5, 2), LinearBins(-0.5, 1.5, 2))))
		>>> x[0, 0] = 0
		>>> x[0, 1] = 1
		>>> x[1, 0] = 2
		>>> x[1, 1] = 4
		>>> x.at_centres()
		array([[ 0.,  1.],
		       [ 2.,  4.]])
		>>> (x + x).at_centres()
		array([[ 0.,  2.],
		       [ 4.,  8.]])
		"""
		self = self.copy()
		self += other
		return self

	def copy(self):
		"""
		Return a copy of the BinnedArray.  The .bins attribute is
		shared with the original.
		"""
		return type(self)(self.bins, self.array.copy())

	def centres(self):
		"""
		Return a tuple of arrays containing the bin centres for
		each dimension.
		"""
		return self.bins.centres()

	def at_centres(self):
		"""
		Return an array of the BinnedArray's value evaluated at the
		bin centres.  In many cases this is simply a reference to
		the internal array object, but for subclasses that override
		item retrieval and assignment some additional work might be
		required to obtain this array.  In those cases, this method
		is a convenience wrapper to avoid coding the evaluation
		logic in the calling code.

		Because subclasses expect to be able to override this, in
		almost all cases calling code that wishes to access the
		values stored in the internal array directly should
		probably use this method to do so.

		NOTE:

		- The return value might be a newly-constructed object or a
		  reference to an internal object.
		"""
		return self.array

	def argmin(self):
		"""
		Return the co-ordinates of the bin centre containing the
		minimum value.  Same as numpy.argmin(), converting the
		indexes to bin co-ordinates.
		"""
		array = self.at_centres()
		return tuple(centres[index] for centres, index in zip(self.centres(), numpy.unravel_index(array.argmin(), array.shape)))

	def argmax(self):
		"""
		Return the co-ordinates of the bin centre containing the
		maximum value.  Same as numpy.argmax(), converting the
		indexes to bin co-ordinates.
		"""
		array = self.at_centres()
		return tuple(centres[index] for centres, index in zip(self.centres(), numpy.unravel_index(array.argmax(), array.shape)))

	def to_pdf(self):
		"""
		Normalize the internal array's contents so that when
		multiplied by the corresponding bin volumes the result sums
		to 1 (neglecting bins with infinite volume).

		NOTE:

		- This is a legacy method that has been superceded by the
		  BinnedDensity and BinnedLnPDF classes.  You almost
		  certainly want to be using those instead of whatever
		  you're doing that needs this method.
		"""
		# zero bins whose volumes are infinite so the rest will
		# appear to be normalized
		self.array[numpy.isinf(self.bins.volumes())] = 0.
		# make sum = 1
		self.array /= self.array.sum()
		# make integral = 1
		self.array /= self.bins.volumes()

	def logregularize(self, epsilon = 2**-1074):
		"""
		Find bins <= 0, and set them to epsilon, This has the
		effect of allowing the logarithm of the array to be
		evaluated without error.

		NOTE:

		- This is a legacy method that has been superceded by the
		  BinnedDensity and BinnedLnPDF classes.  You almost
		  certainly want to be using those instead of whatever
		  you're doing that needs this method.
		"""
		self.array[self.array <= 0] = epsilon
		return self

	def to_xml(self, name):
		"""
		Retrun an XML document tree describing a rate.BinnedArray
		object.
		"""
		elem = ligolw.LIGO_LW()
		elem.Name = u"%s:pylal_rate_binnedarray" % name
		self.bins.to_xml(elem)
		elem.appendChild(ligolw_array.Array.build(u"array", self.array))
		return elem

	@classmethod
	def get_xml_root(cls, xml, name):
		name = u"%s:pylal_rate_binnedarray" % name
		elem = [elem for elem in xml.getElementsByTagName(ligolw.LIGO_LW.tagName) if elem.hasAttribute(u"Name") and elem.Name == name]
		try:
			elem, = elem
		except ValueError:
			raise ValueError("XML tree at '%s' must contain exactly one '%s' LIGO_LW element" % (repr(xml), name))
		return elem

	@classmethod
	def from_xml(cls, xml, name):
		"""
		Search for the description of a rate.BinnedArray object
		named "name" in the XML document tree rooted at xml, and
		construct and return a new rate.BinnedArray object from the
		data contained therein.

		NOTE:  the .array attribute is a reference to the .array
		attribute of the XML element.  Changes to the contents of
		the BinnedArray object affect the XML document tree.
		"""
		elem = cls.get_xml_root(xml, name)
		self = cls(NDBins.from_xml(elem), array = ligolw_array.get_array(elem, u"array").array)
		# sanity check
		if self.bins.shape != self.array.shape:
			raise ValueError("'%s' binning shape does not match array shape:  %s != %s" % (name, self.bins.shape, self.array.shape))
		# done
		return self


#
# =============================================================================
#
#                          Binned Array Interpolator
#
# =============================================================================
#


def InterpBinnedArray(binnedarray, fill_value = 0.0):
	"""
	Wrapper constructing a scipy.interpolate interpolator from the
	contents of a BinnedArray.  Only piecewise linear interpolators are
	supported.  In 1 and 2 dimensions, scipy.interpolate.interp1d and
	.interp2d is used, respectively.  In more than 2 dimensions
	scipy.interpolate.LinearNDInterpolator is used.

	Example:

	One dimension

	>>> x = BinnedArray(NDBins((LinearBins(-0.5, 2.5, 3),)))
	>>> x[0,] = 0
	>>> x[1,] = 1
	>>> x[2,] = 3
	>>> y = InterpBinnedArray(x)
	>>> y(0)
	0.0
	>>> y(1)
	1.0
	>>> y(2)
	3.0
	>>> y(0.5)
	0.5
	>>> y(1.5)
	2.0

	Two dimensions

	>>> x = BinnedArray(NDBins((LinearBins(-0.5, 2.5, 3), LinearBins(-0.5, 1.5, 2))))
	>>> x[0, 0] = 0
	>>> x[0, 1] = 1
	>>> x[1, 0] = 2
	>>> x[1, 1] = 4
	>>> x[2, 0] = 2
	>>> x[2, 1] = 4
	>>> y = InterpBinnedArray(x)
	>>> y(0, 0)
	0.0
	>>> y(0, 1)
	1.0
	>>> y(1, 0)
	2.0
	>>> y(1, 1)
	4.0
	>>> y(2, 0)
	2.0
	>>> y(2, 1)
	4.0
	>>> y(0, 0.25)
	0.25
	>>> y(0, 0.75)
	0.75
	>>> y(0.25, 0)
	0.5
	>>> y(0.75, 0)
	1.5
	>>> y(0.25, 1)
	1.75
	>>> y(0.75, 1)
	3.25
	>>> y(1, 0.25)
	2.5
	>>> y(1, 0.75)
	3.5

	BUGS:  Due to bugs in some versions of scipy and numpy, if an old
	version of scipy and/or numpy is detected this code falls back to
	home-grown piece-wise linear interpolator code for 1- and 2
	dimensions that is slow, and in 3- and higher dimensions the
	fall-back is to nearest-neighbour "interpolation".
	"""
	# the upper and lower boundaries of the binnings are added as
	# additional co-ordinates with the array being assumed to equal
	# fill_value at those points.  this solves the problem of providing
	# a valid function in the outer halves of the first and last bins.

	# coords[0] = co-ordinates along 1st dimension,
	# coords[1] = co-ordinates along 2nd dimension,
	# ...
	coords = tuple(numpy.hstack((l[0], c, u[-1])) for l, c, u in zip(binnedarray.bins.lower(), binnedarray.bins.centres(), binnedarray.bins.upper()))

	# pad the contents of the binned array with 1 element of fill_value
	# on each side in each dimension
	z = binnedarray.at_centres()
	z = numpy.pad(z, [(1, 1)] * z.ndim, mode = "constant", constant_values = [(fill_value, fill_value)] * z.ndim)

	# if any co-ordinates are infinite, remove them.  also remove
	# degenerate co-ordinates from ends
	slices = []
	for c in coords:
		finite_indexes, = numpy.isfinite(c).nonzero()
		assert len(finite_indexes) != 0

		lo, hi = finite_indexes.min(), finite_indexes.max()

		while lo < hi and c[lo + 1] == c[lo]:
			lo += 1
		while lo < hi and c[hi - 1] == c[hi]:
			hi -= 1
		assert lo < hi

		slices.append(slice(lo, hi + 1))
	coords = tuple(c[s] for c, s in zip(coords, slices))
	z = z[tuple(slices)]

	# build the interpolator from the co-ordinates and array data.
	# scipy/numpy interpolators return an array-like thing so we have
	# to wrap them, in turn, in a float cast
	if len(coords) == 1:
		try:
			interp1d
		except NameError:
			# FIXME:  remove when we can rely on a new-enough scipy
			coords0 = coords[0]
			lo, hi = coords[0][0], coords[0][-1]
			with numpy.errstate(invalid = "ignore"):
				dz_over_dcoords0 = (z[1:] - z[:-1]) / (coords0[1:] - coords0[:-1])
			isinf = math.isinf
			def interp(x):
				if not lo < x < hi:
					return fill_value
				i = coords0.searchsorted(x) - 1
				if isinf(z[i]):
					return z[i]
				return z[i] + (x - coords0[i]) * dz_over_dcoords0[i]
			return interp
		interp = interp1d(coords[0], z, kind = "linear", copy = False, bounds_error = False, fill_value = fill_value)
		return lambda *coords: float(interp(*coords))
	elif len(coords) == 2:
		try:
			interp2d
		except NameError:
			# FIXME:  remove when we can rely on a new-enough scipy
			coords0 = coords[0]
			coords1 = coords[1]
			lox, hix = coords0[0], coords0[-1]
			loy, hiy = coords1[0], coords1[-1]
			dcoords0 = coords0[1:] - coords0[:-1]
			dcoords1 = coords1[1:] - coords1[:-1]
			with numpy.errstate(invalid = "ignore"):
				dz0 = z[1:,:] - z[:-1,:]
				dz1 = z[:,1:] - z[:,:-1]
			isinf = math.isinf
			def interp(x, y):
				if not (lox < x < hix and loy < y < hiy):
					return fill_value
				i = coords0.searchsorted(x) - 1
				j = coords1.searchsorted(y) - 1
				dx = (x - coords0[i]) / dcoords0[i]
				dy = (y - coords1[j]) / dcoords1[j]
				if dx + dy <= 1.:
					if isinf(z[i, j]):
						return z[i, j]
					return z[i, j] + dx * dz0[i, j] + dy * dz1[i, j]
				if isinf(z[i + 1, j + 1]):
					return z[i + 1, j + 1]
				return z[i + 1, j + 1] + (1. - dx) * -dz0[i, j + 1] + (1. - dy) * -dz1[i + 1, j]
			return interp
		interp = interp2d(coords[0], coords[1], z.T, kind = "linear", copy = False, bounds_error = False, fill_value = fill_value)
		return lambda *coords: float(interp(*coords))
	else:
		try:
			LinearNDInterpolator
		except NameError:
			# FIXME:  remove when we can rely on a new-enough scipy
			def interp(*coords):
				try:
					return binnedarray[coords]
				except IndexError:
					return fill_value
			return interp
		interp = LinearNDInterpolator(list(itertools.product(*coords)), z.flat, fill_value = fill_value)
		return lambda *coords: float(interp(*coords))


#
# =============================================================================
#
#                                BinnedDensity
#
# =============================================================================
#


class BinnedDensity(BinnedArray):
	"""
	Variant of the BinnedArray type that interprets its contents as a
	density.  When initialized, the volumes of the NDBins used to
	create the BinnedDensity are computed and stored in the .volume
	attribute.  When a value is retrieved from a bin, the value
	reported is the value stored in the bin divided by the bin's
	volume.  When a value is assigned to a bin, the value recorded is
	the assigned value multiplied by the bin's volume.  The true values
	recorded in the bins can be accessed via the .count attribute,
	which is a BinnedArray object wrapping the same array and binning.

	Because the internal array stores counts, not the densities, when
	the data in an instance of this class is processed with the
	filter_array() function the result is the density of smoothed
	counts, not the smoothed density.  This is best illustrated with an
	example.

	Example:  linear bins

	>>> # bins are 2 units, each
	>>> x = BinnedDensity(NDBins((LinearBins(0, 10, 5),)))
	>>> x.volume
	array([ 2.,  2.,  2.,  2.,  2.])
	>>> # set count at 5 to 1
	>>> x.count[5.0,] = 1
	>>> # internal array set to 1 in that bin
	>>> x.array
	array([ 0.,  0.,  1.,  0.,  0.])
	>>> # density reported is 0.5
	>>> print x[5.0,]
	0.5
	>>> # convolve counts with 3-bin top hat window
	>>> filter_array(x.array, tophat_window(3))
	array([ 0.        ,  0.33333333,  0.33333333,  0.33333333,  0.        ])
	>>> # density at 5 is now 1/6 counts / unit interval
	>>> print x[5.0,]
	0.166666666667
	>>> # total count has been preserved
	>>> print x.array.sum()
	1.0

	Example:  logarithmic bins

	>>> # bins increase in size by a factor of 2, each
	>>> x = BinnedDensity(NDBins((LogarithmicBins(1, 32, 5),)))
	>>> x.volume
	array([  1.,   2.,   4.,   8.,  16.])
	>>> # set count at 5 to 1
	>>> x.count[5.0,] = 1
	>>> # internal array set to 1 in that bin
	>>> x.array
	array([ 0.,  0.,  1.,  0.,  0.])
	>>> # density reported is 0.25
	>>> print x[5.0,]
	0.25
	>>> # convolve counts with 3-bin top hat window
	>>> filter_array(x.array, tophat_window(3))
	array([ 0.        ,  0.33333333,  0.33333333,  0.33333333,  0.        ])
	>>> # density at 5 is now 1/12 counts / unit interval
	>>> print x[5.0,]
	0.0833333333333
	>>> # density is 1/6 in bin below
	>>> print x[3,]
	0.166666666667
	>>> # and 1/24 in bin above
	>>> print x[10,]
	0.0416666666667
	>>> # total count has been preserved
	>>> print x.array.sum()
	1.0

	Explanation.  In the linear case there are five bins spanning the
	interval [0, 10], making each bin 2 "units" in size.  A single
	count is placed at 5.0, which is bin number 2.  The averaging
	filter is a top-hat window 3 bins wide.  The single count in bin
	#2, when averaged over the three bins around it, becomes an average
	of 1/3 count per bin, each of which is 2 units in size, so the
	average density is 1/6 events / unit.  The count has been spread
	out over several bins but the integral of the density, the total
	count, is unchanged.  The logarithmic case is identical but because
	the bin sizes are non-uniform, when the single count is spread
	across the three bins the density is non-uniform.  The integral is
	still preserved.

	Some binnings have infinite-sized bins.  For such bins, the counts
	may be manipulated directly, as usual, but for all finite counts a
	density of 0 will be reported for those bins (and NaN for infinite
	counts), and ValueError will be raised if an attempt is made to
	assign a density to those bins.

	NOTES:

	- While it is technically possible to modify the binning parameters
	  after creating an instance of this class, the steps required to
	  bring all internal data back into consistency following such a
	  change are undocumented.  One should consider the metadata
	  carried by these objects to be immutable.
	"""
	def __init__(self, *args, **kwargs):
		super(BinnedDensity, self).__init__(*args, **kwargs)
		self.count = BinnedArray(self.bins, array = self.array)
		self.volume = self.bins.volumes()

	def __getitem__(self, coords):
		coords = self.bins(*coords)
		return self.array[coords] / self.volume[coords]

	def __setitem__(self, coords, val):
		coords = self.bins(*coords)
		vol = self.volume[coords]
		if numpy.isinf(vol).any():
			raise ValueError("cannot assign density values to infinite-volume bins, try assigning a count instead")
		self.array[coords] = val * vol

	def at_centres(self):
		return self.array / self.volume

	def marginalize(self, dim):
		"""
		Return a new BinnedDensity object containing the density
		integrated over dimension dim.

		Example:

		>>> # 5x5 mesh of bins, each with volume = 4
		>>> x = BinnedDensity(NDBins((LinearBins(0, 10, 5), LinearBins(0, 10, 5))))
		>>> # set count at 5,5 to 1
		>>> x.count[5.0, 5.0] = 1
		>>> # density in central bin is 1/4
		>>> x.at_centres()
		array([[ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ],
		       [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ],
		       [ 0.  ,  0.  ,  0.25,  0.  ,  0.  ],
		       [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ],
		       [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ]])
		>>> # convolve counts with 3-bin top-hat window
		>>> filter_array(x.array, tophat_window(3, 3))
		array([[ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ],
		       [ 0.        ,  0.11111111,  0.11111111,  0.11111111,  0.        ],
		       [ 0.        ,  0.11111111,  0.11111111,  0.11111111,  0.        ],
		       [ 0.        ,  0.11111111,  0.11111111,  0.11111111,  0.        ],
		       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ]])
		>>> # density is now 1/(4 * 9) = 1/36 in 9 central bins
		>>> x.at_centres()
		array([[ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ],
		       [ 0.        ,  0.02777778,  0.02777778,  0.02777778,  0.        ],
		       [ 0.        ,  0.02777778,  0.02777778,  0.02777778,  0.        ],
		       [ 0.        ,  0.02777778,  0.02777778,  0.02777778,  0.        ],
		       [ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ]])
		>>> # densities still sum (not integrate) to 1/4
		>>> x.at_centres().sum()
		0.25
		>>> # integrate over dimension 1
		>>> x = x.marginalize(1)
		>>> # bin volumes are now 2
		>>> # densities in 3 central bins are = 1/(2 * 3) = 1/6
		>>> x.at_centres()
		array([ 0.        ,  0.16666667,  0.16666667,  0.16666667,  0.        ])
		>>> # densities sum (not integrate) to 1/2
		>>> x.at_centres().sum()
		0.5
		"""
		return type(self)(NDBins(self.bins[:dim] + self.bins[dim+1:]), self.array.sum(axis = dim))


#
# A discretely sampled PDF.
#


class BinnedLnPDF(BinnedDensity):
	"""
	Variant of the BinnedDensity class that (i) tracks a normalization
	which it uses to rescale the reported density so that its integral
	is one, and (ii) reports the natural logarithm of that normalized
	density.

	The .normalize() method needs to be invoked after any manipulation
	of the contents of the array before the results of .__getitem__()
	are meaningful.

	How the normalization is tracked should be assumed to be
	undocumented.  In this class the .norm attribute contains the
	natural log of the sum of the counts in all bins and this value is
	subtracted from the natural log of the (unnormalized) density in
	the .__getitem__() method, but subclasses are free to implement
	whatever unique mechanism is appropriate for themselves.

	As with the BinnedDensity class, the internal array contains counts
	(not densities, nor natural logarithms of densities), and the
	.count attribute continues to be a BinnedArray interface to those
	counts.  The intention is for the counts themselves to provide an
	additional degree of freedom apart from the normalized density.
	For example, see the .__iadd__() method where it is assumed that
	the total count encodes a relative weight to be used when
	marginalizing over two PDFs.

	Example:

	>>> # 5x5 mesh of bins each with volume = 4
	>>> x = BinnedLnPDF(NDBins((LinearBins(0, 10, 5), LinearBins(0, 10, 5))))
	>>> # set count at 5,5 to 36 and normalize
	>>> x.count[5.0, 5.0] = 36
	>>> x.normalize()
	>>> # log probability density = ln 1/4 = -1.3862943611198906
	>>> x.at_centres()
	array([[       -inf,        -inf,        -inf,        -inf,        -inf],
	       [       -inf,        -inf,        -inf,        -inf,        -inf],
	       [       -inf,        -inf, -1.38629436,        -inf,        -inf],
	       [       -inf,        -inf,        -inf,        -inf,        -inf],
	       [       -inf,        -inf,        -inf,        -inf,        -inf]])
	>>> # convolve with 3x3 top-hat window.  in general one must renormalize after this, but we'll skip that here because the demo is constructed so as to not need it
	>>> filter_array(x.array, tophat_window(3, 3))
	array([[ 0.,  0.,  0.,  0.,  0.],
	       [ 0.,  4.,  4.,  4.,  0.],
	       [ 0.,  4.,  4.,  4.,  0.],
	       [ 0.,  4.,  4.,  4.,  0.],
	       [ 0.,  0.,  0.,  0.,  0.]])
	>>> # ln probability density = ln 1/(4 * 9) = -3.58351893845611
	>>> x.at_centres()
	array([[       -inf,        -inf,        -inf,        -inf,        -inf],
	       [       -inf, -3.58351894, -3.58351894, -3.58351894,        -inf],
	       [       -inf, -3.58351894, -3.58351894, -3.58351894,        -inf],
	       [       -inf, -3.58351894, -3.58351894, -3.58351894,        -inf],
	       [       -inf,        -inf,        -inf,        -inf,        -inf]])
	>>> # .marginzlize() preserves normalization
	>>> y = x.marginalize(1)
	>>> y.count.at_centres()
	array([  0.,  12.,  12.,  12.,   0.])
	>>> # ln probability density = ln 1/(2 * 3) = -1.791759469228055
	>>> y.at_centres()
	array([       -inf, -1.79175947, -1.79175947, -1.79175947,        -inf])
	>>> # assuming \\sqrt{N} counting fluctuations, compute the fractional uncertainty
	>>> import numpy
	>>> d = BinnedArray(x.bins, 1. / numpy.sqrt(x.count.at_centres()))
	>>> d.at_centres()
	array([[ inf,  inf,  inf,  inf,  inf],
	       [ inf,  0.5,  0.5,  0.5,  inf],
	       [ inf,  0.5,  0.5,  0.5,  inf],
	       [ inf,  0.5,  0.5,  0.5,  inf],
	       [ inf,  inf,  inf,  inf,  inf]])
	"""
	def __init__(self, *args, **kwargs):
		super(BinnedLnPDF, self).__init__(*args, **kwargs)
		self.normalize()

	def __getitem__(self, coords):
		return numpy.log(super(BinnedLnPDF, self).__getitem__(coords)) - self.norm

	def __setitem__(self, coords, val):
		#
		# the relationship between the density, p, the
		# volume, Delta, the count, x, and the overall
		# normalization, N, for a bin i is
		#
		# p_i = x_i / (Delta_i * N)
		#
		# where N = sum_j x_j, and 0 <= p_i * Delta_i <= 1.
		#
		# to maintain the requirement that the density integrate to
		# 1, it is not possible to change the density in one bin
		# without changing the densities in all other bins.
		# because there is no unique way to assign new densities to
		# all other bins that preserves the normalization invariant
		# there is an ambiguity in how to implement the assignment
		# operation.  instead of choosing one, we elect to forbid
		# this operation.
		#
		# one interpretation could be to require the counts in all
		# other bins to be preserved, and solve for the new count
		# for the bin in question (and new total count).  because
		# the counts in other bins are held constant the likelihood
		# ratios among them are preserved.
		#
		# solving for x_i,
		#
		#         p_i * Delta_i
		# x_i = ----------------- sum_(j != i) x_j.
		#       1 - p_i * Delta_i
		#
		# the normalization would then need to be recomputed given
		# the new count for bin i.  forbidden cases include:
		# infinite bin size, total count is initially 0,
		# p_i*Delta_i = 1 unless bin i is the only non-zero bin to
		# begin with.
		#
		# another interpretation could be to require the total
		# count to be preserved and also the likelihood ratios
		# among all other bins.  because the total count is being
		# preserved, the new x_i can be obtained immediately as
		#
		# x_i' = N * p_i * Delta_i
		#
		# the counts in all other bins are obtained by rescaling by
		# the after-to-before ratio of the non-bin-i count.
		#
		# x_(j != i) = x_(j != i) * (N - x_i') / (N - x_i)
		#
		# because they are scaled by a common factor, the ratios
		# between them are preserved.  forbidden cases include:
		# infinite bin size, total count is initially 0.
		#
		raise NotImplementedError("item assignment operation not defined.  assign to .count then invoke .normalize()")

	def mkinterp(self):
		"""
		Return an interpolator to evaluate the density as a smooth
		function of co-ordinates.  If the density has not been
		normalized the interpolator's behaviour is undefined.

		NOTE:  the interpolator is the InterpBinnedArray object
		which might have limitations in 3 or more dimensions.  See
		its documentation for more information.

		NOTE:  in the future this is likely to be replaced with
		some sort of internal mechanism.

		Example:

		>>> # 5-bin linear mesh of bins each with volume = 2
		>>> x = BinnedLnPDF(NDBins((LinearBins(0, 10, 5), )))
		>>> # set some counts and normalize
		>>> x.count[3,] = x.count[7,] = 1
		>>> x.count[5,] = 2
		>>> x.normalize()
		>>> x.count.at_centres()
		array([ 0.,  1.,  2.,  1.,  0.])
		>>> x.at_centres()
		array([       -inf, -2.07944154, -1.38629436, -2.07944154,        -inf])
		>>> # construct interpolator object
		>>> x = x.mkinterp()
		>>> # log(P) is interpolated linearly, so it can be counter-intuitive where the density drops to 0 in neighbouring bins
		>>> x(3)
		-inf
		>>> # otherwise behaviour is as expected
		>>> x(4)
		-1.732867951399863
		>>> x(4.5)
		-1.5595811562598769
		>>> x(5)
		-1.3862943611198906
		"""
		return InterpBinnedArray(self)

	def at_centres(self):
		with numpy.errstate(divide = "ignore", invalid = "ignore"):
			return numpy.log(super(BinnedLnPDF, self).at_centres()) - self.norm

	def marginalize(self, dim):
		new = super(BinnedLnPDF, self).marginalize(dim)
		new.norm = self.norm
		return new

	def __iadd__(self, other):
		"""
		Adds the counts array and normalization of other to self.
		If the original two PDFs were normalized then the result is
		also normalized, otherwise .normalize() needs to be invoked
		to normalize the result.  Once normalized, the result is
		the PDF marginalized over the two original PDFs (the sum of
		the two PDFs weighted by the relative total frequency of
		events in each).

		Example:

		>>> # 5-bin linear mesh of bins each with volume = 2
		>>> x = BinnedLnPDF(NDBins((LinearBins(0, 10, 5), )))
		>>> y = BinnedLnPDF(NDBins((LinearBins(0, 10, 5), )))
		>>> x.count[5,] = 2
		>>> y.count[3,] = 2
		>>> x.normalize()
		>>> y.normalize()
		>>> x.at_centres()
		array([       -inf,        -inf, -0.69314718,        -inf,        -inf])
		>>> x += y
		>>> x.at_centres()
		array([       -inf, -1.38629436, -1.38629436,        -inf,        -inf])
		"""
		super(BinnedLnPDF, self).__iadd__(other)
		# c = a + b
		# if b is smaller than a:
		# c = a (1 + b/a)
		# log c = log a + log1p(b / a)
		# log c = log a + log1p(exp(log b - log a))
		# otherwise, swap a and b.
		if self.norm >= other.norm:
			self.norm += math.log1p(math.exp(other.norm - self.norm))
		else:
			self.norm = other.norm + math.log1p(math.exp(self.norm - other.norm))
		return self

	def __add__(self, other):
		self = super(BinnedLnPDF, self).__add__(other)
		self.normalize()
		return self

	def copy(self):
		new = super(BinnedLnPDF, self).copy()
		new.norm = self.norm
		return new

	def normalize(self):
		"""
		Updates the internal normalization.  Subclasses override
		this with custom normalization mechanisms if required.

		Note that counts contained in infinite-sized bins are
		included in the normalization.  In this way the relative
		frequency with which counts occur in those bins is
		accounted for in the normalization although the density
		reported for those bins will be 0.
		"""
		self.norm = self.array.sum()
		assert self.norm >= 0. or math.isnan(self.norm)
		self.norm = math.log(self.norm) if self.norm != 0. else NegInf

	def to_xml(self, *args, **kwargs):
		elem = super(BinnedLnPDF, self).to_xml(*args, **kwargs)
		elem.appendChild(ligolw_param.Param.from_pyvalue("norm", self.norm))
		return elem

	@classmethod
	def from_xml(cls, xml, name):
		elem = cls.get_xml_root(xml, name)
		self = super(BinnedLnPDF, cls).from_xml(elem, name)
		self.norm = ligolw_param.get_pyvalue(elem, "norm")
		return self


#
# =============================================================================
#
#                                   Windows
#
# =============================================================================
#


def gaussian_window(*bins, **kwargs):
	"""
	Generate a normalized (integral = 1) Gaussian window in N
	dimensions.  The bins parameters set the width of the window in bin
	counts in each dimension.  The optional keyword argument sigma,
	which defaults to 10, sets the size of the array in all dimensions
	in units of the width in each dimension.  The sizes are adjusted so
	that the array has an odd number of samples in each dimension, and
	the Gaussian is peaked on the middle sample.

	Example:

	>>> # 2D window with width of 1.5 bins in first dimension,
	>>> # 1 bin in second dimension, 3 widths long (rounded to odd
	>>> # integer = 5 x 3 bins) in each dimension
	>>> gaussian_window(1.5, 1, sigma = 3)
	array([[ 0.00161887,  0.01196189,  0.00161887],
	       [ 0.02329859,  0.17215456,  0.02329859],
	       [ 0.05667207,  0.41875314,  0.05667207],
	       [ 0.02329859,  0.17215456,  0.02329859],
	       [ 0.00161887,  0.01196189,  0.00161887]])
	"""
	if not bins:
		raise ValueError("function requires at least 1 width")
	if min(bins) < 0.:
		raise ValueError("widths must be non-negative, got %s" % str(bins))
	sigma = kwargs.pop("sigma", 10)
	if kwargs:
		raise ValueError("unrecognized keyword argument(s): %s" % ",".join(kwargs))
	windows = []
	for b in bins:
		l = int(math.floor(sigma * b / 2.0)) * 2
		w = lal.CreateGaussREAL8Window(l + 1, l / float(b) if b else PosInf)
		windows.append(w.data.data / w.sum)
	if len(windows) == 1:
		# 1D short-cut
		return windows[0]
	try:
		# only works upto 26 dimensions, but that's 2 trillion bins
		# if there are just 3 bins along each side, so it's
		# unlikely to be a practical limitation;  for a while at
		# least
		assert len(windows) <= 26
		return numpy.einsum(",".join("abcdefghijklmnopqrstuvwxyz"[:len(windows)]), *windows)
	except AttributeError:
		# numpy < 1.6
		window = reduce(numpy.outer, windows)
		window.shape = tuple(len(w) for w in windows)
		return window


def tophat_window(*bins):
	"""
	Generate a normalized (integral = 1) rectangular window in N
	dimensions.  The bins parameters set the width of the window in bin
	counts in each dimension, each of which must be positive and will
	be rounded up to the nearest odd integer.

	Example:

	>>> tophat_window(4)
	array([ 0.2,  0.2,  0.2,  0.2,  0.2])
	>>> tophat_window(4, 4)
	array([[ 0.04,  0.04,  0.04,  0.04,  0.04],
	       [ 0.04,  0.04,  0.04,  0.04,  0.04],
	       [ 0.04,  0.04,  0.04,  0.04,  0.04],
	       [ 0.04,  0.04,  0.04,  0.04,  0.04],
	       [ 0.04,  0.04,  0.04,  0.04,  0.04]])
	"""
	if not bins:
		raise ValueError("function requires at least 1 width")
	if min(bins) <= 0:
		raise ValueError("widths must be positive, got %s" % str(bins))
	windows = []
	for b in bins:
		w = lal.CreateRectangularREAL8Window(int(math.floor(b / 2.0)) * 2 + 1)
		windows.append(w.data.data / w.sum)
	if len(windows) == 1:
		# 1D short-cut
		return windows[0]
	try:
		# only works upto 26 dimensions, but that's 2 trillion bins
		# if there are just 3 bins along each side, so it's
		# unlikely to be a practical limitation;  for a while at
		# least
		assert len(windows) <= 26
		return numpy.einsum(",".join("abcdefghijklmnopqrstuvwxyz"[:len(windows)]), *windows)
	except AttributeError:
		# numpy < 1.6
		window = reduce(numpy.outer, windows)
		window.shape = tuple(len(w) for w in windows)
		return window


#
# =============================================================================
#
#                                  Filtering
#
# =============================================================================
#


def filter_array(a, window, cyclic = False, use_fft = True):
	"""
	Filter an array using the window function.  The transformation is
	done in place.  The data are assumed to be 0 outside of their
	domain of definition.  The window function must have an odd number
	of samples in each dimension;  this is done so that it is always
	clear which sample is at the window's centre, which helps prevent
	phase errors.  If the window function's size exceeds that of the
	data in one or more dimensions, the largest allowed central portion
	of the window function in the affected dimensions will be used.
	This is done silently;  to determine if window function truncation
	will occur, check for yourself that your window function is smaller
	than your data in all dimensions.

	If use_fft is True, the window is convolved with the array using
	FFT convolution.  This is done by processing the array in bands of
	relatively small dynamic range each to work around numerical
	dynamic range limitation in the FFTs, but the window function is
	not treated similarly.  As a result the window is effectively
	limited to a few orders of magnitude of dynamic range.  If use_fft
	is False there are no dynamic range issues but the convolution can
	be quite slow.
	"""
	assert not cyclic	# no longer supported, maybe in future
	# check that the window and the data have the same number of
	# dimensions
	dims = len(a.shape)
	if dims != len(window.shape):
		raise ValueError("array and window dimensions mismatch")
	# check that all of the window's dimensions have an odd size
	if 0 in map((1).__and__, window.shape):
		raise ValueError("window size is not an odd integer in at least 1 dimension")
	# determine how much of the window function can be used
	window_slices = []
	for d in range(dims):
		if window.shape[d] > a.shape[d]:
			# largest odd integer <= size of a
			n = ((a.shape[d] + 1) // 2) * 2 - 1
			first = (window.shape[d] - n) // 2
			window_slices.append(slice(first, first + n))
		else:
			window_slices.append(slice(0, window.shape[d]))
	window = window[tuple(window_slices)]

	if use_fft:
		# this loop works around dynamic range limits in the FFT
		# convolution code.  we move data 4 orders of magnitude at
		# a time from the original array into a work space,
		# convolve the work space with the filter, zero the
		# workspace in any elements that are more than 14 orders of
		# magnitude below the maximum value in the result, and add
		# the result to the total.
		result = numpy.zeros_like(a)
		while a.any():
			# copy contents of input array to work space
			workspace = numpy.copy(a)

			# mask = indexes of elements of work space not more
			# than 4 orders of magnitude larger than the
			# smallest non-zero element.  these are the
			# elements to be processed in this iteration
			abs_workspace = abs(workspace)
			mask = abs_workspace <= abs_workspace[abs_workspace > 0].min() * 1e4
			del abs_workspace

			# zero the masked elements in the input array, zero
			# everything except the masked elements in the work
			# space
			a[mask] = 0.
			workspace[~mask] = 0.
			del mask

			# convolve the work space with the kernel
			workspace = signaltools.fftconvolve(workspace, window, mode = "same")

			# determine the largest value in the work space,
			# and set to zero anything more than 14 orders of
			# magnitude smaller
			abs_workspace = abs(workspace)
			workspace[abs_workspace < abs_workspace.max() * 1e-14] = 0.
			del abs_workspace

			# sum what remains into the result
			result += workspace
			del workspace
	else:
		result = signaltools.convolve(a, window, mode = "same")
	# overwrite the input with the result
	numpy.copyto(a, result, casting = "no")

	return a
