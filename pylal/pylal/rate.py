# Copyright (C) 2006--2014  Kipp Cannon
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


from bisect import bisect_right
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


from glue import iterutils
from glue import segments
from glue.ligolw import ligolw
from glue.ligolw import array as ligolw_array
from glue.ligolw import param as ligolw_param
from glue.ligolw import types as ligolw_types
import lal
from pylal import git_version


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
	def __init__(self, min, max, n):
		"""
		Initialize a Bins instance.  The three arguments are the
		minimum and maximum of the values spanned by the bins, and
		the number of bins to place between them.  Subclasses may
		require additional arguments, or different arguments
		altogether.
		"""
		# convenience code to do some common initialization and
		# input checking
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

	def __cmp__(self, other):
		"""
		Two binnings are the same if they are instances of the same
		class, have the same lower and upper bounds, and the same
		count of bins.
		"""
		if not isinstance(other, type(self)):
			return -1
		return cmp((type(self), self.min, self.max, len(self)), (type(other), other.min, other.max, len(other)))

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
		"""
		if isinstance(x, slice):
			if x.step is not None:
				raise NotImplementedError("step not supported: %s" % repr(x))
			return slice(self[x.start] if x.start is not None else 0, self[x.stop] + 1 if x.stop is not None else len(self))
		raise NotImplementedError

	def __iter__(self):
		"""
		If __iter__ does not exist, Python uses __getitem__ with
		range(0) as input to define iteration. This is nonsensical
		for bin objects, so explicitly unsupport iteration.
		"""
		raise NotImplementedError

	def lower(self):
		"""
		Return an array containing the locations of the lower
		boundaries of the bins.
		"""
		raise NotImplementedError

	def centres(self):
		"""
		Return an array containing the locations of the bin
		centres.
		"""
		raise NotImplementedError

	def upper(self):
		"""
		Return an array containing the locations of the upper
		boundaries of the bins.
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
		>>> print "%.15g" % math.log(1./10)
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
		>>> print "%.15g" % math.log(1. / (x.upper()[-2] - 0.5))
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
		For internal use by XML I/O code.
		"""
		return u"%s:%s" % (name, suffix)

	@staticmethod
	def xml_bins_name_dec(name, suffix = u"pylal_rate_bins"):
		"""
		For internal use by XML I/O code.
		"""
		name = name.rsplit(u":", 1)
		if name[-1] != suffix:
			raise ValueError(name)
		return name[0]

	@classmethod
	def xml_bins_check(cls, elem, name):
		"""
		For internal use by XML I/O code.
		"""
		return elem.tagName == ligolw.Param.tagName and elem.hasAttribute(u"Name") and name == cls.xml_bins_name_dec(elem.Name)

	def to_xml(self):
		"""
		Construct a LIGO Light Weight XML representation of the
		Bins instance.
		"""
		raise NotImplementedError

	@classmethod
	def from_xml(cls, xml):
		"""
		From the XML Param element at xml, return the Bins object
		it describes.
		"""
		raise NotImplementedError


class LoHiCountToFromXMLMixin(object):
	"""
	For internal use by XML I/O code.
	"""
	def to_xml(self):
		"""
		Construct a LIGO Light Weight XML representation of the
		Bins instance.
		"""
		return ligolw_param.from_pyvalue(self.xml_bins_name_enc(self.xml_bins_name), u"%s,%s,%s" % (ligolw_types.FormatFunc[u"real_8"](self.min), ligolw_types.FormatFunc[u"real_8"](self.max), ligolw_types.FormatFunc[u"int_8s"](self.n)))

	@classmethod
	def from_xml(cls, xml):
		"""
		From the XML Param element at xml, return the Bins object
		it describes.
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
		boundaries = tuple(boundaries)
		if any(a > b for a, b in zip(boundaries[:-1], boundaries[1:])):
			raise ValueError("non-monotonic boundaries provided")

		self.boundaries = boundaries

	def __cmp__(self, other):
		"""
		Two binnings are the same if they are instances of the same
		class, and have the same boundaries.
		"""
		if not isinstance(other, type(self)):
			return -1
		return cmp(self.boundaries, other.boundaries)

	def __len__(self):
		return len(self.boundaries) - 1

	def __getitem__(self, x):
		if isinstance(x, slice):
			return super(IrregularBins, self).__getitem__(x)
		if self.boundaries[0] <= x < self.boundaries[-1]:
			return bisect_right(self.boundaries, x) - 1
		# special measure-zero edge case
		if x == self.boundaries[-1]:
			return len(self.boundaries) - 2
		raise IndexError(x)

	def lower(self):
		return numpy.array(self.boundaries[:-1])

	def upper(self):
		return numpy.array(self.boundaries[1:])

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
		return ligolw_param.from_pyvalue(self.xml_bins_name_enc(self.xml_bins_name), u",".join(map(ligolw_types.FormatFunc[u"real_8"], self.boundaries)))

	@classmethod
	def from_xml(cls, xml):
		"""
		From the XML Param element at xml, return the Bins object
		it describes.
		"""
		if not cls.xml_bins_check(xml, cls.xml_bins_name):
			raise ValueError("not a %s" % repr(cls))
		return cls(map(ligolw_types.ToPyType[u"real_8"], xml.pcdata.split(u",")))


class LinearBins(LoHiCountToFromXMLMixin, Bins):
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
		if isinstance(x, slice):
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


class LinearPlusOverflowBins(LoHiCountToFromXMLMixin, Bins):
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
		if isinstance(x, slice):
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
		return numpy.concatenate((numpy.array([NegInf]), self.min + self.delta * numpy.arange(len(self) - 2), numpy.array([self.max])))

	def centres(self):
		return numpy.concatenate((numpy.array([NegInf]), self.min + self.delta * (numpy.arange(len(self) - 2) + 0.5), numpy.array([PosInf])))

	def upper(self):
		return numpy.concatenate((numpy.array([self.min]), self.min + self.delta * (numpy.arange(len(self) - 2) + 1), numpy.array([PosInf])))

	#
	# XML I/O related methods and data
	#

	xml_bins_name = u"linplusoverflowbins"


class LogarithmicBins(LoHiCountToFromXMLMixin, Bins):
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
		self.delta = (math.log(max) - math.log(min)) / n

	def __getitem__(self, x):
		if isinstance(x, slice):
			return super(LogarithmicBins, self).__getitem__(x)
		if self.min <= x < self.max:
			return int(math.floor((math.log(x) - math.log(self.min)) / self.delta))
		if x == self.max:
			# special "measure zero" corner case
			return len(self) - 1
		raise IndexError(x)

	def lower(self):
		return numpy.exp(numpy.linspace(math.log(self.min), math.log(self.max) - self.delta, len(self)))

	def centres(self):
		return numpy.exp(numpy.linspace(math.log(self.min), math.log(self.max) - self.delta, len(self)) + self.delta / 2.)

	def upper(self):
		return numpy.exp(numpy.linspace(math.log(self.min) + self.delta, math.log(self.max), len(self)))

	#
	# XML I/O related methods and data
	#

	xml_bins_name = u"logbins"


class LogarithmicPlusOverflowBins(LoHiCountToFromXMLMixin, Bins):
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
		self.delta = (math.log(max) - math.log(min)) / (n - 2)

	def __getitem__(self, x):
		if isinstance(x, slice):
			return super(LogarithmicPlusOverflowBins, self).__getitem__(x)
		if self.min <= x < self.max:
			return 1 + int(math.floor((math.log(x) - math.log(self.min)) / self.delta))
		if x >= self.max:
			# infinity overflow bin
			return len(self) - 1
		if x < self.min:
			# zero overflow bin
			return 0
		raise IndexError(x)

	def lower(self):
		return numpy.concatenate((numpy.array([0.]), numpy.exp(numpy.linspace(math.log(self.min), math.log(self.max), len(self) - 1))))

	def centres(self):
		return numpy.concatenate((numpy.array([0.]), numpy.exp(numpy.linspace(math.log(self.min), math.log(self.max) - self.delta, len(self) - 2) + self.delta / 2.), numpy.array([PosInf])))

	def upper(self):
		return numpy.concatenate((numpy.exp(numpy.linspace(math.log(self.min), math.log(self.max), len(self) - 1)), numpy.array([PosInf])))

	#
	# XML I/O related methods and data
	#

	xml_bins_name = u"logplusoverflowbins"


class ATanBins(LoHiCountToFromXMLMixin, Bins):
	"""
	Bins spaced uniformly in tan^-1 x.  Provides approximately linear
	binning in the middle portion, with the bin density dropping
	asymptotically to 0 as x goes to +/- \infty.  The min and max
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
		if isinstance(x, slice):
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


class ATanLogarithmicBins(LoHiCountToFromXMLMixin, IrregularBins):
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
		super(ATanLogarithmicBins, self).__init__(boundaries[:-1][keepers])
		self.keepers = keepers[:-1]
		self.min = min
		self.max = max
		self.n = n

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

	>>> from glue.segments import *
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
		categories is an iterable of containers defining the
		categories.  (Recall that containers are collections that
		support the "in" operator.) Objects will be mapped to the
		integer index of the container that contains them.
		"""
		self.containers = tuple(categories)  # need to set an order and len

	def __len__(self):
		return len(self.containers)

	def __getitem__(self, value):
		"""
		Return i if value is contained in i-th container. If value
		is not contained in any of the containers, raise an
		IndexError.
		"""
		for i, s in enumerate(self.containers):
			if value in s:
				return i
		raise IndexError(value)

	def __cmp__(self, other):
		if not isinstance(other, type(self)):
			return -1
		return cmp(self.containers, other.containers)

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
		# FIXME:  make use of new "pickle" type for params when we
		# can rely on a new-enough glue
		#return ligolw_param.Param.build(self.xml_bins_name_enc(self.xml_bins_name), u"pickle", self.containers)
		import pickle
		return ligolw_param.from_pyvalue(self.xml_bins_name_enc(self.xml_bins_name), pickle.dumps(self.containers))

	@classmethod
	def from_xml(cls, xml):
		"""
		From the XML Param element at xml, return the Bins object
		it describes.
		"""
		if not cls.xml_bins_check(xml, cls.xml_bins_name):
			raise ValueError("not a %s" % repr(cls))
		# FIXME:  replace with commented-out code when we can rely
		# on new "pickle" type for params
		#return cls(xml.pcdata)
		import pickle
		return cls(pickle.loads(xml.pcdata))


class HashableBins(Categories):
	"""
	Maps hashable objects (things that can be used as dictionary keys) to integers.

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
		<class 'pylal.rate.LogarithmicBins'>

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
		indices using the Bins instances in this object.

		Example:

		>>> x = NDBins((LinearBins(1, 25, 3), LogarithmicBins(1, 25, 3)))
		>>> x(1, 1)
		(0, 0)
		>>> x = NDBins((LinearBins(1, 25, 3),))
		>>> x(1)
		(0,)
		>>> x = NDBins((LinearBins(1, 25, 1000),))
		>>> # slices require manual construction
		>>> x(slice(10, 12))
		(slice(375, 459, None),)
		>>> x = NDBins((Categories([set(("Cow", "Chicken", "Goat")), set(("Tractor", "Plough")), set(("Barn", "House"))]),))
		>>> x("Cow")
		(0,)

		Each co-ordinate can be anything the corresponding Bins
		instance will accept.
		"""
		if len(coords) != len(self):
			raise ValueError("dimension mismatch")
		return tuple(b[c] for b, c in zip(self, coords))

	@property
	def shape(self):
		return tuple(len(b) for b in self)

	def lower(self):
		"""
		Return a tuple of arrays, where each array contains the
		locations of the lower boundaries of the bins in the
		corresponding dimension.
		"""
		return tuple(b.lower() for b in self)

	def centres(self):
		"""
		Return a tuple of arrays, where each array contains the
		locations of the bin centres for the corresponding
		dimension.
		"""
		return tuple(b.centres() for b in self)

	def upper(self):
		"""
		Return a tuple of arrays, where each array contains the
		locations of the upper boundaries of the bins in the
		corresponding dimension.
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
	xml_bins_name_mapping.update(zip(xml_bins_name_mapping.values(), xml_bins_name_mapping.keys()))

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


def bins_spanned(bins, seglist, dtype = "double"):
	"""
	Input is a Bins subclass instance and a glue.segments.segmentlist
	instance.  The output is an array object the length of the binning,
	which each element in the array set to the interval in the
	corresponding bin spanned by the segment list.

	Example:

	>>> from glue.segments import *
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
	# make an intersection of the segment list with the extend of the bins
	# need to use lower/upper instead of min/max because the latter sometimes
	# merely correspond to low and high parameters used to construct the binning
	# (see, for example, the atan binning)
	seglist = seglist & segments.segmentlist([segments.segment(lower[0], upper[-1])])
	array = numpy.zeros((len(bins),), dtype = dtype)
	for i, (a, b) in enumerate(zip(lower, upper)):
		array[i] = abs(seglist & segments.segmentlist([segments.segment(a, b)]))
	return array


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
	>>> x.array
	array([ 0.,  0.,  0.,  0.,  0.])
	>>> x[0,] += 1
	>>> x[0.5,] += 1
	>>> x.array
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
	>>> x.array
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
	ValueError: input array and input bins must have the same shape:  (5, 1) != (5,)

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
		else:
			if array.shape != bins.shape:
				raise ValueError("input array and input bins must have the same shape:  %s != %s" % (str(array.shape), str(bins.shape)))
			self.array = array

	def __getitem__(self, coords):
		return self.array[self.bins[coords]]

	def __setitem__(self, coords, val):
		self.array[self.bins[coords]] = val

	def __len__(self):
		return len(self.array)

	def __iadd__(self, other):
		"""
		Add the contents of another BinnedArray object to this one.
		Both must have identical binnings.
		"""
		if self.bins != other.bins:
			raise TypeError("incompatible binning: %s" % repr(other))
		self.array += other.array
		return self
		# here's an implementation that allows the binnings to
		# differ.  each bin in other is added to whichever bin in
		# self its centre is found in.  this behaviour probably
		# leads to undesirable results if other's binning is less
		# dense than self's.  would need to spread other's bins'
		# contents out somehow.  probably there's no behaviour that
		# is correct for all use cases.
		#for coords in iterutils.MultiIter(*other.bins.centres()):
		#	self[coords] += other[coords]

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

	def argmin(self):
		"""
		Return the co-ordinates of the bin centre containing the
		minimum value.  Same as numpy.argmin(), converting the
		indexes to bin co-ordinates.
		"""
		return tuple(centres[index] for centres, index in zip(self.centres(), numpy.unravel_index(self.array.argmin(), self.array.shape)))

	def argmax(self):
		"""
		Return the co-ordinates of the bin centre containing the
		maximum value.  Same as numpy.argmax(), converting the
		indexes to bin co-ordinates.
		"""
		return tuple(centres[index] for centres, index in zip(self.centres(), numpy.unravel_index(self.array.argmax(), self.array.shape)))

	def to_density(self):
		"""
		Divide each bin's value by the volume of the bin.
		"""
		self.array /= self.bins.volumes()

	def to_pdf(self):
		"""
		Convert into a probability density.
		"""
		# zero bins whose volumes are infinite so the rest will
		# appear to be normalized
		self.array[numpy.isinf(self.bins.volumes())] = 0.
		# make sum = 1
		self.array /= self.array.sum()
		# make integral = 1
		self.to_density()

	def logregularize(self, epsilon = 2**-1074):
		"""
		Find bins <= 0, and set them to epsilon, This has the
		effect of allowing the logarithm of the array to be
		evaluated without error.
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
		elem.appendChild(ligolw_array.from_array(u"array", self.array))
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
		name = u"%s:pylal_rate_binnedarray" % name
		elem = [elem for elem in xml.getElementsByTagName(ligolw.LIGO_LW.tagName) if elem.hasAttribute(u"Name") and elem.Name == name]
		try:
			elem, = elem
		except ValueError:
			raise ValueError("XML tree at '%s' must contain exactly one '%s' LIGO_LW element" % (repr(xml), name))
		self = cls(NDBins.from_xml(elem), array = ligolw_array.get_array(elem, u"array").array)
		# sanity check
		if self.bins.shape != self.array.shape:
			raise ValueError("'%s' binning shape does not match array shape:  %s != %s" % (name, self.bins.shape, self.array.shape))
		# done
		return self


class BinnedRatios(object):
	"""
	Like BinnedArray, but provides a numerator array and a denominator
	array.  The incnumerator() method increments a bin in the numerator
	by the given weight, and the incdenominator() method increments a
	bin in the denominator by the given weight.  There are no methods
	provided for setting or decrementing either, but the they are
	accessible as the numerator and denominator attributes, which are
	both BinnedArray objects.
	"""
	def __init__(self, bins, dtype = "double"):
		self.numerator = BinnedArray(bins, dtype = dtype)
		self.denominator = BinnedArray(bins, dtype = dtype)

	def __getitem__(self, coords):
		return self.numerator[coords] / self.denominator[coords]

	def bins(self):
		return self.numerator.bins

	def __iadd__(self, other):
		"""
		Add the weights from another BinnedRatios object's
		numerator and denominator to the numerator and denominator
		of this one.  Note that this is not the same as adding the
		ratios.  It is not necessary for the binnings to be
		identical, but an integer number of the bins in other must
		fit into each bin in self.
		"""
		try:
			self.numerator += other.numerator
			self.denominator += other.denominator
		except TypeError:
			raise TypeError("incompatible binning: %s" % repr(other))
		return self

	def incnumerator(self, coords, weight = 1):
		"""
		Add weight to the numerator bin at coords.
		"""
		self.numerator[coords] += weight

	def incdenominator(self, coords, weight = 1):
		"""
		Add weight to the denominator bin at coords.
		"""
		self.denominator[coords] += weight

	def ratio(self):
		"""
		Compute and return the array of ratios.
		"""
		return self.numerator.array / self.denominator.array

	def regularize(self):
		"""
		Find bins in the denominator that are 0, and set them to 1.
		Presumably the corresponding bin in the numerator is also
		0, so this has the effect of allowing the ratio array to be
		evaluated without error, returning zeros in those bins that
		have had no weight added to them.
		"""
		self.denominator.array[self.denominator.array == 0] = 1
		return self

	def logregularize(self, epsilon = 2**-1074):
		"""
		Find bins in the denominator that are 0, and set them to 1,
		while setting the corresponding bin in the numerator to
		float epsilon.  This has the effect of allowing the
		logarithm of the ratio array to be evaluated without error.
		"""
		self.numerator.array[self.denominator.array == 0] = epsilon
		self.denominator.array[self.denominator.array == 0] = 1
		return self

	def centres(self):
		"""
		Return a tuple of arrays containing the bin centres for
		each dimension.
		"""
		return self.numerator.bins.centres()

	def used(self):
		"""
		Return the number of bins with non-zero denominator.
		"""
		return numpy.sum(self.denominator.array != 0)

	def to_pdf(self):
		"""
		Convert the numerator and denominator into a pdf.
		"""
		self.numerator.to_pdf()
		self.denominator.to_pdf()

	def to_xml(self, name):
		"""
		Return an XML document tree describing a rate.BinnedRatios object.
		"""
		xml = ligolw.LIGO_LW({u"Name": u"%s:pylal_rate_binnedratios" % name})
		xml.appendChild(self.numerator.to_xml(u"numerator"))
		xml.appendChild(self.denominator.to_xml(u"denominator"))
		return xml

	@classmethod
	def from_xml(cls, xml, name):
		"""
		Search for the description of a rate.BinnedRatios object named
		"name" in the XML document tree rooted at xml, and construct and
		return a new rate.BinnedRatios object from the data contained
		therein.
		"""
		xml, = [elem for elem in xml.getElementsByTagName(ligolw.LIGO_LW.tagName) if elem.hasAttribute(u"Name") and elem.Name == u"%s:pylal_rate_binnedratios" % name]
		self = cls(NDBins())
		self.numerator = BinnedArray.from_xml(xml, u"numerator")
		self.denominator = BinnedArray.from_xml(xml, u"denominator")
		# normally they share a single NDBins instance
		self.denominator.bins = self.numerator.bins
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
	try:
		z = numpy.pad(binnedarray.array, [(1, 1)] * len(binnedarray.array.shape), mode = "constant", constant_values = [(fill_value, fill_value)] * len(binnedarray.array.shape))
	except AttributeError:
		# numpy < 1.7 didn't have pad().  FIXME:  remove when we
		# can rely on a newer numpy
		z = numpy.empty(tuple(l + 2 for l in binnedarray.array.shape))
		z.fill(fill_value)
		z[(slice(1, -1),) * len(binnedarray.array.shape)] = binnedarray.array

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
	z = z[slices]

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
	sigma = kwargs.pop("sigma", 10)
	if kwargs:
		raise ValueError("unrecognized keyword argument(s): %s" % ",".join(kwargs))
	windows = []
	for b in bins:
		if b <= 0:
			raise ValueError("negative width: %s" % repr(b))
		l = int(math.floor(sigma * b / 2.0)) * 2
		w = lal.CreateGaussREAL8Window(l + 1, l / float(b))
		windows.append(w.data.data / w.sum)
	if len(windows) == 1:
		# 1D short-cut
		return windows[0]
	try:
		return numpy.einsum(",".join("abcdefghijklmnopqrstuvwxyz"[:len(windows)]), *windows)
	except AttributeError:
		# numpy < 1.6
		window = reduce(numpy.outer, windows)
		window.shape = tuple(len(w) for w in windows)
		return window


def tophat_window(bins):
	"""
	Generate a normalized (integral = 1) top-hat window in 1 dimension.
	bins sets the width of the window in bin counts, which is rounded
	up to the nearest odd integer.

	Example:

	>>> tophat_window(4)
	array([ 0.2,  0.2,  0.2,  0.2,  0.2])
	"""
	if bins <= 0:
		raise ValueError(bins)
	w = lal.CreateRectangularREAL8Window(int(math.floor(bins / 2.0)) * 2 + 1)
	return w.data.data / w.sum


def tophat_window2d(bins_x, bins_y):
	"""
	Generate a normalized (integral = 1) top-hat window in 2
	dimensions.  bins_x and bins_y set the widths of the window in bin
	counts, which are both rounded up to the nearest odd integer.  The
	result is a rectangular array, with an elliptical pattern of
	elements set to a constant value centred on the array's mid-point,
	and all other elements set to 0.
	"""
	if bins_x <= 0:
		raise ValueError(bins_x)
	if bins_y <= 0:
		raise ValueError(bins_y)

	# This might appear to be using a screwy, slow, algorithm but it's
	# the only way I have found to get a window with the correct bins
	# set and cleared as appropriate.  I'd love this to be replaced by
	# something that's easier to know is correct.

	# fill rectangle with ones, making the number of bins odd in each
	# direction
	window = numpy.ones((int(bins_x / 2.0) * 2 + 1, int(bins_y / 2.0) * 2 + 1), "Float64")

	# zero the bins outside the window
	for x, y in iterutils.MultiIter(*map(range, window.shape)):
		if ((x - window.shape[0] // 2) / float(bins_x) * 2.0)**2 + ((y - window.shape[1] // 2) / float(bins_y) * 2.0)**2 > 1.0:
			window[x, y] = 0.0

	# normalize
	window /= window.sum()

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
	for d in xrange(dims):
		if window.shape[d] > a.shape[d]:
			# largest odd integer <= size of a
			n = ((a.shape[d] + 1) // 2) * 2 - 1
			first = (window.shape[d] - n) // 2
			window_slices.append(slice(first, first + n))
		else:
			window_slices.append(slice(0, window.shape[d]))
	window = window[window_slices]

	if use_fft:
		# this loop works around dynamic range limits in the FFT
		# convolution code.  we move data 4 orders of magnitude at a time
		# from the original array into a work space, convolve the work
		# space with the filter, zero the workspace in any elements that
		# are more than 14 orders of magnitude below the maximum value in
		# the result, and add the result to the total.
		result = numpy.zeros_like(a)
		while a.any():
			# copy contents of input array to work space
			workspace = numpy.copy(a)

			# mask = indexes of elements of work space not more than 4
			# orders of magnitude larger than the smallest non-zero
			# element.  these are the elements to be processed in this
			# iteration
			abs_workspace = abs(workspace)
			mask = abs_workspace <= abs_workspace[abs_workspace > 0].min() * 1e4
			del abs_workspace

			# zero the masked elements in the input array, zero
			# everything except the masked elements in the work space
			a[mask] = 0.
			workspace[~mask] = 0.
			del mask

			# convolve the work space with the kernel
			workspace = signaltools.fftconvolve(workspace, window, mode = "same")

			# determine the largest value in the work space, and set to
			# zero anything more than 14 orders of magnitude smaller
			abs_workspace = abs(workspace)
			workspace[abs_workspace < abs_workspace.max() * 1e-14] = 0.
			del abs_workspace

			# sum what remains into the result
			result += workspace
			del workspace
	else:
		result = signaltools.convolve(a, window, mode = "same")
	# overwrite the input with the result
	# FIXME:  in numpy >= 1.7.0 there is a copyto() function
	a.flat = result.flat

	return a


def filter_binned_ratios(ratios, window, cyclic = False):
	"""
	Convolve the numerator and denominator of a BinnedRatios instance
	each with the same window function.  This has the effect of
	interpolating the ratio of the two between bins where it has been
	measured, weighting bins by the number of measurements made in
	each.  For example, consider a 1-dimensional binning, with zeros in
	the denominator and numerator bins everywhere except in one bin
	where both are set to 1.0.  The ratio is 1.0 in that bin, and
	undefined everywhere else, where it has not been measured.
	Convolving both numerator and denominator with a Gaussian window
	will replace the "delta function" in each with a smooth hill
	spanning some number of bins.  Since the same smooth hill will be
	seen in both the numerator and the denominator bins, the ratio of
	the two is now 1.0 --- the ratio from the bin where a measurement
	was made --- everywhere the window function had support.  Contrast
	this to the result of convolving the ratio with a window function.

	Convolving the numerator and denominator bins separately preserves
	the integral of each.  In other words the total number of events in
	each of the denominator and numerator is conserved, only their
	locations are shuffled about.  Convolving, instead, the ratios with
	a window function would preserve the integral of the ratio, which
	is probably meaningless.

	Note that you should be using the window functions defined in this
	module, which are carefully designed to be norm preserving (the
	integrals of the numerator and denominator bins are preserved), and
	phase preserving.

	Note, also, that you should apply this function *before* using
	either of the regularize() methods of the BinnedRatios object.
	"""
	filter_array(ratios.numerator.array, window, cyclic = cyclic)
	filter_array(ratios.denominator.array, window, cyclic = cyclic)


#
# =============================================================================
#
#                                    Rates
#
# =============================================================================
#


def to_moving_mean_density(binned_array, filterdata, cyclic = False):
	"""
	Convolve a BinnedArray with a filter function, then divide all bins
	by their volumes.  The result is the density function smoothed by
	the filter.  The default is to assume 0 values beyond the ends of
	the array when convolving with the filter function.  Set the
	optional cyclic parameter to True for periodic boundaries.

	Example:

	>>> x = BinnedArray(NDBins((LinearBins(0, 10, 5),)))
	>>> x[5.0,] = 1
	>>> x.array
	array([ 0.,  0.,  1.,  0.,  0.])
	>>> to_moving_mean_density(x, tophat_window(3))
	>>> x.array
	array([ 0.        ,  0.16666667,  0.16666667,  0.16666667,  0.        ])

	Explanation.  There are five bins spanning the interval [0, 10],
	making each bin 2 "units" in size.  A single count is placed at
	5.0, which is bin number 2.  The averaging filter is a top-hat
	window 3 bins wide.  The single count in bin #2, when averaged over
	the three bins around it, is equivalent to a mean density of 1/6
	events / unit.

	Example:

	>>> x = BinnedArray(NDBins((LinearBins(0, 10, 5),)))
	>>> x[1,] = 1
	>>> x[3,] = 1
	>>> x[5,] = 1
	>>> x[7,] = 1
	>>> x[9,] = 1
	>>> x.array
	array([ 1.,  1.,  1.,  1.,  1.])
	>>> to_moving_mean_density(x, tophat_window(3))
	>>> x.array
	array([ 0.33333333,  0.5       ,  0.5       ,  0.5       ,  0.33333333])

	We have uniformly distributed events at 2 unit intervals (the first
	is at 1, the second at 3, etc.).  The event density is 0.5 events /
	unit, except at the edges where the smoothing window has picked up
	zero values from beyond the ends of the array.
	"""
	filter_array(binned_array.array, filterdata, cyclic = cyclic)
	binned_array.to_density()


def marginalize(pdf, dim):
	"""
	From a BinnedArray object containing probability density data (bins
	whose volume integral is 1), return a new BinnedArray object
	containing the probability density marginalized over dimension
	dim.
	"""
	dx = pdf.bins[dim].upper() - pdf.bins[dim].lower()
	dx_shape = [1] * len(pdf.bins)
	dx_shape[dim] = len(dx)
	dx.shape = dx_shape

	result = BinnedArray(NDBins(pdf.bins[:dim] + pdf.bins[dim+1:]))
	result.array = (pdf.array * dx).sum(axis = dim)

	return result


def marginalize_ratios(likelihood, dim):
	"""
	Marginalize the numerator and denominator of a BinnedRatios object
	containing likelihood-ratio data (i.e., the numerator and
	denominator both contain probability density data) over dimension
	dim.
	"""
	result = BinnedRatios(NDBins())
	result.numerator = marginalize(likelihood.numerator, dim)
	result.denominator = marginalize(likelihood.denominator, dim)
	# normally they share an NDBins instance
	result.denominator.bins = result.numerator.bins
	return result
