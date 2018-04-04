# Copyright (C) 2006-2013  Kipp Cannon
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


"""
This module contains bits and pieces of use when interacting with LAL and
LAL-derived code (eg. LALApps programs)
"""


import fnmatch
import math
import os
import re
import sys
import urlparse
import warnings

try:  # python < 3
    long
except NameError:  # python >= 3
    long = int

from glue import git_version
from glue import segments


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                          High precision time object
#
# =============================================================================
#


#
# Python version in case LAL isn't available
#


class LIGOTimeGPS(object):
	"""
	An object for storing times with nanosecond resolution.  LAL defines an
	equivalent object which is used through-out the search algorithms to
	represent times.  Many LALApps routines input and output times in a
	manner that meshes well with this object.

	Internally the time is represented as a signed integer "seconds" part
	and an unsigned integer "nanoseconds" part.  The actual time is always
	constructed by adding the nanoseconds to the seconds.  So -0.5 s is
	represented by setting seconds = -1, and nanoseconds to 500000000.
	That's the way LAL does it.
	"""

	# basic class methods

	def __init__(self, seconds, nanoseconds = 0):
		"""
		Create a LIGOTimeGPS instance.  The first parameter is the
		count of seconds, and the second (optional) parameter is the
		count of nanoseconds.  If the nanoseconds parameter is not
		supplied, it is assumed to be 0.  Either parameter can be
		a numeric type or an ASCII string.

		Example:

		>>> LIGOTimeGPS(100.5)
		LIGOTimeGPS(100, 500000000)
		>>> LIGOTimeGPS("100.5")
		LIGOTimeGPS(100, 500000000)
		>>> LIGOTimeGPS(100, 500000000)
		LIGOTimeGPS(100, 500000000)
		>>> LIGOTimeGPS(0, 100500000000)
		LIGOTimeGPS(100, 500000000)
		>>> LIGOTimeGPS(100.2, 300000000)
		LIGOTimeGPS(100, 500000000)
		>>> LIGOTimeGPS("0.000000001")
		LIGOTimeGPS(0, 1)
		>>> LIGOTimeGPS("0.0000000012")
		LIGOTimeGPS(0, 1)
		>>> LIGOTimeGPS("0.0000000018")
		LIGOTimeGPS(0, 2)
		>>> LIGOTimeGPS("-0.8")
		LIGOTimeGPS(-1, 200000000)
		>>> LIGOTimeGPS("-1.2")
		LIGOTimeGPS(-2, 800000000)
		"""
		if not isinstance(nanoseconds, (float, int, long)):
			try:
				nanoseconds = float(nanoseconds)
			except:
				raise TypeError(nanoseconds)
		if isinstance(seconds, float):
			ns, seconds = math.modf(seconds)
			seconds = int(seconds)
			nanoseconds += ns * 1e9
		elif not isinstance(seconds, (int, long)):
			if isinstance(seconds, (str, unicode)):
				sign = -1 if seconds.lstrip().startswith("-") else +1
				try:
					if "." in seconds:
						seconds, ns = seconds.split(".")
						ns = round(sign * float("." + ns) * 1e9)
					else:
						ns = 0
					seconds = int(seconds)
				except:
					raise TypeError("invalid literal for LIGOTimeGPS(): %s" % seconds)
				nanoseconds += ns
			elif hasattr(seconds, "gpsSeconds") and hasattr(seconds, "gpsNanoSeconds"):
				# handle LIGOTimeGPS(x) where x is an
				# object with gpsSeconds and gpsNanoSeconds
				# fields.
				nanoseconds += seconds.gpsNanoSeconds
				seconds = seconds.gpsSeconds
			elif hasattr(seconds, "seconds") and hasattr(seconds, "nanoseconds"):
				# handle LIGOTimeGPS(x) where x is an
				# object with seconds and nanoseconds
				# fields.
				nanoseconds += seconds.nanoseconds
				seconds = seconds.seconds
			else:
				raise TypeError(seconds)
		self.__seconds = seconds + int(nanoseconds // 1000000000)
		self.__nanoseconds = int(nanoseconds % 1000000000)

	seconds = gpsSeconds = property(lambda self: self.__seconds)
	nanoseconds = gpsNanoSeconds = property(lambda self: self.__nanoseconds)

	def __repr__(self):
		return "LIGOTimeGPS(%d, %u)" % (self.__seconds, self.__nanoseconds)
	
	def __str__(self):
		"""
		Return an ASCII string representation of a LIGOTimeGPS.
		"""
		if (self.__seconds >= 0) or (self.__nanoseconds == 0):
			s = "%d.%09u" % (self.__seconds, self.__nanoseconds)
		elif self.__seconds < -1:
			s = "%d.%09u" % (self.__seconds + 1, 1000000000 - self.__nanoseconds)
		else:
			s = "-0.%09u" % (1000000000 - self.__nanoseconds)
		return s.rstrip("0").rstrip(".")

	# type conversion

	def __float__(self):
		"""
		Convert a LIGOTimeGPS to seconds as a float.

		Example:

		>>> float(LIGOTimeGPS(100.5))
		100.5
		"""
		return self.__seconds + self.__nanoseconds * 1e-9

	def __int__(self):
		"""
		Return the integer part (seconds) of a LIGOTimeGPS as an int.

		Example:

		>>> int(LIGOTimeGPS(100.5))
		100
		"""
		return self.__seconds

	def __long__(self):
		"""
		Return the integer part (seconds) of a LIGOTimeGPS as a long.

		Example:

		>>> long(LIGOTimeGPS(100.5))
		100L
		"""
		return long(self.__seconds)

	def ns(self):
		"""
		Convert a LIGOTimeGPS to a count of nanoseconds as a long.

		Example:

		>>> LIGOTimeGPS(100.5).ns()
		100500000000
		"""
		return self.__seconds * 1000000000 + self.__nanoseconds

	# comparison

	def __hash__(self):
		return self.__seconds ^ self.__nanoseconds

	def __cmp__(self, other):
		"""
		Compare a value to a LIGOTimeGPS.  If the value being compared
		to the LIGOTimeGPS is not also a LIGOTimeGPS, then an attempt
		is made to convert it to a LIGOTimeGPS.

		Example:

		>>> LIGOTimeGPS(100.5) < LIGOTimeGPS(200)
		True
		>>> LIGOTimeGPS(100.5) < 200
		True
		>>> LIGOTimeGPS(100.5) < "200"
		True
		"""
		if not isinstance(other, LIGOTimeGPS):
			try:
				other = LIGOTimeGPS(other)
			except TypeError:
				return NotImplemented
		return cmp(self.__seconds, other.seconds) or cmp(self.__nanoseconds, other.nanoseconds)

	def __nonzero__(self):
		"""
		Return True if the LIGOTimeGPS is nonzero.

		Example:

		>>> bool(LIGOTimeGPS(100.5))
		True
		"""
		return self.__seconds or self.__nanoseconds

	# arithmetic

	def __add__(self, other):
		"""
		Add a value to a LIGOTimeGPS.  If the value being added to the
		LIGOTimeGPS is not also a LIGOTimeGPS, then an attempt is made
		to convert it to a LIGOTimeGPS.

		Example:

		>>> LIGOTimeGPS(100.5) + LIGOTimeGPS(3)
		LIGOTimeGPS(103, 500000000)
		>>> LIGOTimeGPS(100.5) + 3
		LIGOTimeGPS(103, 500000000)
		>>> LIGOTimeGPS(100.5) + "3"
		LIGOTimeGPS(103, 500000000)
		"""
		if not isinstance(other, LIGOTimeGPS):
			other = LIGOTimeGPS(other)
		return LIGOTimeGPS(self.__seconds + other.seconds, self.__nanoseconds + other.nanoseconds)

	# addition is commutative.
	__radd__ = __add__

	def __sub__(self, other):
		"""
		Subtract a value from a LIGOTimeGPS.  If the value being
		subtracted from the LIGOTimeGPS is not also a LIGOTimeGPS, then
		an attempt is made to convert it to a LIGOTimeGPS.

		Example:

		>>> LIGOTimeGPS(100.5) - LIGOTimeGPS(3)
		LIGOTimeGPS(97, 500000000)
		>>> LIGOTimeGPS(100.5) - 3
		LIGOTimeGPS(97, 500000000)
		>>> LIGOTimeGPS(100.5) - "3"
		LIGOTimeGPS(97, 500000000)
		"""
		if not isinstance(other, LIGOTimeGPS):
			other = LIGOTimeGPS(other)
		return LIGOTimeGPS(self.__seconds - other.seconds, self.__nanoseconds - other.nanoseconds)

	def __rsub__(self, other):
		"""
		Subtract a LIGOTimeGPS from a value.
		"""
		if not isinstance(other, LIGOTimeGPS):
			other = LIGOTimeGPS(other)
		return LIGOTimeGPS(other.seconds - self.__seconds, other.nanoseconds - self.__nanoseconds)

	def __mul__(self, other):
		"""
		Multiply a LIGOTimeGPS by a number.

		Example:

		>>> LIGOTimeGPS(100.5) * 2
		LIGOTimeGPS(201, 0)
		"""
		seconds = self.__seconds
		nanoseconds = self.__nanoseconds

		if seconds < 0 and nanoseconds > 0:
			seconds += 1
			nanoseconds -= 1000000000
		elif seconds > 0 and nanoseconds < 0:
			seconds -=1
			nanoseconds += 1000000000

		slo = seconds % 131072
		shi = seconds - slo
		olo = other % 2**(int(math.log(other, 2)) - 26) if other else 0
		ohi = other - olo

		nanoseconds *= float(other)
		seconds = 0.
		for addend in (slo * olo, shi * olo, slo * ohi, shi * ohi):
			n, s = math.modf(addend)
			seconds += s
			nanoseconds += n * 1e9

		return LIGOTimeGPS(seconds, round(nanoseconds))

	# multiplication is commutative
	__rmul__ = __mul__

	def __div__(self, other):
		"""
		Divide a LIGOTimeGPS by a number.

		Example:

		>>> LIGOTimeGPS(100.5) / 2
		LIGOTimeGPS(50, 250000000)
		"""
		quotient = LIGOTimeGPS(float(self) / float(other))
		for n in xrange(100):
			residual = float(self - quotient * other) / float(other)
			quotient += residual
			if abs(residual) <= 0.5e-9:
				break
		return quotient

	def __mod__(self, other):
		"""
		Compute the remainder when a LIGOTimeGPS is divided by a number.

		Example:

		>>> LIGOTimeGPS(100.5) % 3
		LIGOTimeGPS(1, 500000000)
		"""
		quotient = int(self / other)
		return self - quotient * other

	# unary arithmetic

	def __pos__(self):
		return self

	def __neg__(self):
		return LIGOTimeGPS(0, -self.ns())

	def __abs__(self):
		if self.__seconds >= 0:
			return self
		return -self


#
# =============================================================================
#
#                         LAL Cache File Manipulation
#
# =============================================================================
#


#
# Representation of a line in a LAL cache file
#


class CacheEntry(object):
	"""
	An object representing one line in a LAL cache file.

	The LAL cache format is defined elsewhere, and what follows is
	meant only to be informative, not an official specification.  Each
	line in a LAL cache identifies a single file, and the line consists
	of five columns of white-space delimited text.

	The first column, "observatory", generally stores the name of an
	observatory site or one or more instruments (preferably delimited
	by ",", but often there is no delimiter between instrument names in
	which case they should be 2 characters each).

	The second column, "description", stores a short string tag that is
	usually all capitals with "_" separating components, in the style
	of the description part of the LIGO-Virgo frame filename format.

	The third and fourth columns store the start time and duration in
	GPS seconds of the interval spanned by the file identified by the
	cache line.  When the file does not start on an integer second or
	its duration is not an integer number of seconds, the conventions
	of the LIGO-Virgo frame filename format apply.

	The fifth (last) column stores the file's URL.

	The values for these columns are stored in the .observatory,
	.description, .segment and .url attributes of instances of this
	class, respectively.  The .segment attribute stores a
	glue.segments.segment object describing the interval spanned by the
	file.  Any of these attributes except the URL is allowed to be
	None.

	This module also provides a Cache class, which some codes use to
	represent an entire LAL cache file.  However, for most use cases a
	simple Python list or set of CacheEntry objects is sufficient for
	representing a LAL cache file.

	Example (parse a string):

	>>> c = CacheEntry("H1 S5 815901601 576.5 file://localhost/home/kipp/tmp/1/H1-815901601-576.xml")
	>>> c.scheme
	'file'
	>>> c.host
	'localhost'

	Example (one-liners to read and write a cache file):

	>>> filename = "874000000-20000.cache"
	>>> cache = map(CacheEntry, open(filename))
	>>> f = open(filename + ".new", "w")
	>>> for cacheentry in cache: print >>f, str(cacheentry)

	Example (extract segmentlist dictionary from LAL cache):

	>>> from glue import segments
	>>> seglists = segments.segmentlistdict()
	>>> for cacheentry in cache:
	...	seglists |= cacheentry.segmentlistdict
	...

	See also:

	glue.segmentsUtils.fromlalcache()
	"""
	# How to parse a line in a LAL cache file.  Five white-space
	# delimited columns.
	_regex = re.compile(r"\A\s*(?P<obs>\S+)\s+(?P<dsc>\S+)\s+(?P<strt>\S+)\s+(?P<dur>\S+)\s+(?P<url>\S+)\s*\Z")
	_url_regex = re.compile(r"\A((.*/)*(?P<obs>[^/]+)-(?P<dsc>[^/]+)-(?P<strt>[^/]+)-(?P<dur>[^/\.]+)\.[^/]+)\Z")

	def __init__(self, *args, **kwargs):
		"""
		Intialize a CacheEntry object.  The arguments can take two
		forms:  a single string argument, which is interpreted and
		parsed as a line from a LAL cache file, or four arguments
		used to explicitly initialize the observatory, description,
		segment and URL in that order.  When parsing a single line
		of text from a LAL cache, an optional key-word argument
		"coltype" can be provided to set the type the start and
		durations are parsed as.  The default is
		glue.lal.LIGOTimeGPS.

		Example:

		>>> c = CacheEntry("H1", "S5", segments.segment(815901601, 815902177.5), "file://localhost/home/kipp/tmp/1/H1-815901601-576.xml")
		>>> print c.segment
		[815901601 ... 815902177.5)
		>>> print str(c)
		H1 S5 815901601 576.5 file://localhost/home/kipp/tmp/1/H1-815901601-576.xml
		>>> c = CacheEntry("H1 S5 815901601 576.5 file://localhost/home/kipp/tmp/1/H1-815901601-576.xml")
		>>> print c.segment
		[815901601 ... 815902177.5)
		>>> print CacheEntry("H1 S5 815901601 576.5 file://localhost/home/kipp/tmp/1/H1-815901601-576.xml", coltype = float).segment
		[815901601.0 ... 815902177.5)

		See also the .from_T050017() class method for an
		alternative initialization mechanism.
		"""
		if len(args) == 1:
			# parse line of text as an entry in a cache file
			match = self._regex.search(args[0])
			try:
				match = match.groupdict()
			except AttributeError:
				raise ValueError("could not convert %s to CacheEntry" % repr(args[0]))
			self.observatory = match["obs"]
			self.description = match["dsc"]
			start = match["strt"]
			duration = match["dur"]
			coltype = kwargs.pop("coltype", LIGOTimeGPS)
			if start == "-" and duration == "-":
				# no segment information
				self.segment = None
			else:
				start = coltype(start)
				self.segment = segments.segment(start, start + coltype(duration))
			self.url = match["url"]
			if kwargs:
				raise TypeError("unrecognized keyword arguments: %s" % ", ".join(kwargs))
		elif len(args) == 4:
			# parse arguments as observatory, description,
			# segment, url
			if kwargs:
				raise TypeError("invalid arguments: %s" % ", ".join(kwargs))
			self.observatory, self.description, self.segment, self.url = args
		else:
			raise TypeError("invalid arguments: %s" % args)

		# "-" indicates an empty column
		if self.observatory == "-":
			self.observatory = None
		if self.description == "-":
			self.description = None


	def __str__(self):
		"""
		Convert the CacheEntry to a string in the format of a line
		in a LAL cache.  Used to write the CacheEntry to a file.

		Example:

		>>> c = CacheEntry("H1 S5 815901601 576.5 file://localhost/home/kipp/tmp/1/H1-815901601-576.xml")
		>>> str(c)
		'H1 S5 815901601 576.5 file://localhost/home/kipp/tmp/1/H1-815901601-576.xml'
		"""
		if self.segment is not None:
			start = str(self.segment[0])
			duration = str(abs(self.segment))
		else:
			start = "-"
			duration = "-"
		return "%s %s %s %s %s" % (self.observatory or "-", self.description or "-", start, duration, self.url)

	def __cmp__(self, other):
		"""
		Compare two CacheEntry objects by observatory, then
		description, then segment, then URL.
		"""
		if not isinstance(other, CacheEntry):
			raise TypeError("can only compare CacheEntry to CacheEntry")
		return cmp((self.observatory, self.description, self.segment, self.url), (other.observatory, other.description, other.segment, other.url))

	def __hash__(self):
		return hash((self.observatory, self.description, self.segment, self.url))

	@property
	def url(self):
		"""
		The cache entry's URL.  The URL is constructed from the
		values of the scheme, host, and path attributes.  Assigning
		a value to the URL attribute causes the value to be parsed
		and the scheme, host and path attributes updated.
		"""
		return urlparse.urlunparse((self.scheme, self.host, self.path, None, None, None))

	@url.setter
	def url(self, url):
		self.scheme, self.host, self.path = urlparse.urlparse(url)[:3]

	@property
	def segmentlistdict(self):
		"""
		A segmentlistdict object describing the instruments and
		time spanned by this CacheEntry.  A new object is
		constructed each time this attribute is accessed (segments
		are immutable so there is no reason to try to share a
		reference to the CacheEntry's internal segment;
		modifications of one would not be reflected in the other
		anyway).

		Example:

		>>> c = CacheEntry(u"H1 S5 815901601 576.5 file://localhost/home/kipp/tmp/1/H1-815901601-576.xml")
		>>> c.segmentlistdict
		{u'H1': [segment(LIGOTimeGPS(815901601, 0), LIGOTimeGPS(815902177, 500000000))]}

		The \"observatory\" column of the cache entry, which is
		frequently used to store instrument names, is parsed into
		instrument names for the dictionary keys using the same
		rules as glue.ligolw.lsctables.instrument_set_from_ifos().

		Example:

		>>> c = CacheEntry(u"H1H2, S5 815901601 576.5 file://localhost/home/kipp/tmp/1/H1H2-815901601-576.xml")
		>>> c.segmentlistdict
		{u'H1H2': [segment(LIGOTimeGPS(815901601, 0), LIGOTimeGPS(815902177, 500000000))]}
		"""
		# the import has to be done here to break the cyclic
		# dependancy
		from glue.ligolw.lsctables import instrument_set_from_ifos
		instruments = instrument_set_from_ifos(self.observatory) or (None,)
		return segments.segmentlistdict((instrument, segments.segmentlist(self.segment is not None and [self.segment] or [])) for instrument in instruments)

	@classmethod
	def from_T050017(cls, url, coltype = LIGOTimeGPS):
		"""
		Parse a URL in the style of T050017-00 into a CacheEntry.
		The T050017-00 file name format is, essentially,

		observatory-description-start-duration.extension

		Example:

		>>> c = CacheEntry.from_T050017("file://localhost/data/node144/frames/S5/strain-L2/LLO/L-L1_RDS_C03_L2-8365/L-L1_RDS_C03_L2-836562330-83.gwf")
		>>> c.observatory
		'L'
		>>> c.host
		'localhost'
		>>> os.path.basename(c.path)
		'L-L1_RDS_C03_L2-836562330-83.gwf'
		"""
		match = cls._url_regex.search(url)
		if not match:
			raise ValueError("could not convert %s to CacheEntry" % repr(url))
		observatory = match.group("obs")
		description = match.group("dsc")
		start = match.group("strt")
		duration = match.group("dur")
		if start == "-" and duration == "-":
			# no segment information
			segment = None
		else:
			segment = segments.segment(coltype(start), coltype(start) + coltype(duration))
		return cls(observatory, description, segment, url)


#
# An object representing a LAL cache file
#


class Cache(list):
	"""
	An object representing a LAL cache file. Currently it is possible to
	add anything to a Cache. This method should check that the thing you
	are adding is a CacheEntry and throw and error if it is not.
	"""
	entry_class = CacheEntry
	
	# methods to create new Cache objects
	def fromfile(cls, fileobj, coltype=LIGOTimeGPS):
		"""
		Return a Cache object whose entries are read from an open file.
		"""
		c = [cls.entry_class(line, coltype=coltype) for line in fileobj]
		return cls(c)
	fromfile = classmethod(fromfile)

	def fromfilenames(cls, filenames, coltype=LIGOTimeGPS):
		"""
		Read Cache objects from the files named and concatenate the results into a
		single Cache.
		"""
		cache = cls()
		for filename in filenames:
			cache.extend(cls.fromfile(open(filename), coltype=coltype))
		return cache
	fromfilenames = classmethod(fromfilenames)

	def from_urls(cls, urllist, coltype=LIGOTimeGPS):
		"""
		Return a Cache whose entries are inferred from the URLs
		in urllist, if possible.  PFN lists will also work; for PFNs, the path
		will be absolutized and "file://" and "localhost" will be assumed
		for the schemes and hosts.
		
		The filenames must be in the format set forth by DASWG in T050017-00.
		"""
		def pfn_to_url(url):
			scheme, host, path, dummy, dummy = urlparse.urlsplit(url)
			if scheme == "": path = os.path.abspath(path)
			return urlparse.urlunsplit((scheme or "file", host or "localhost",
			                            path, "", ""))
		return cls([cls.entry_class.from_T050017(pfn_to_url(f), coltype=coltype) \
		            for f in urllist])
	from_urls = classmethod(from_urls)

	# some set arithmetic to make life better
	def __isub__(self, other):
		"""
		Remove elements from self that are in other.
		"""
		end = len(self) - 1
		for i, elem in enumerate(self[::-1]):
			if elem in other:
				del self[end - i]
		return self

	def __sub__(self, other):
		"""
		Return a Cache containing the entries of self that are not in other.
		"""
		return self.__class__([elem for elem in self if elem not in other])

	def __ior__(self, other):
		"""
		Append entries from other onto self without introducing (new) duplicates.
		"""
		self.extend(other - self)
		return self
	
	def __or__(self, other):
		"""
		Return a Cache containing all entries of self and other.
		"""
		return self.__class__(self[:]).__ior__(other)

	def __iand__(self, other):
		"""
		Remove elements in self that are not in other.
		"""
		end = len(self) - 1
		for i, elem in enumerate(self[::-1]):
			if elem not in other:
				del self[end - i]
		return self
	
	def __and__(self, other):
		"""
		Return a Cache containing the entries of self that are also in other.
		"""
		return self.__class__([elem for elem in self if elem in other])

	def unique(self):
		"""
		Return a Cache which has every element of self, but without
		duplication.  Preserve order.  Does not hash, so a bit slow.
		"""
		new = self.__class__([])
		for elem in self:
			if elem not in new:
				new.append(elem)
		return new

	# other useful manipulations
	def tofile(self, fileobj):
		"""
		write a cache object to the fileobj as a lal cache file
		"""
		for entry in self:
			fileobj.write('%s\n' % str(entry))
		fileobj.close()

	def topfnfile(self, fileobj):
		"""
		write a cache object to filename as a plain text pfn file
		"""
		for entry in self:
			fileobj.write('%s\n' % entry.path)
		fileobj.close()

	def to_segmentlistdict(self):
		"""
		Return a segmentlistdict object describing the instruments
		and times spanned by the entries in this Cache.  The return
		value is coalesced.
		"""
		d = segments.segmentlistdict()
		for entry in self:
			d |= entry.segmentlistdict
		return d

	def sieve(self, ifos=None, description=None, segment=None,
		segmentlist=None, exact_match=False):
		"""
		Return a Cache object with those CacheEntries that
		contain the given patterns (or overlap, in the case of
		segment or segmentlist).  If exact_match is True, then
		non-None ifos, description, and segment patterns must
		match exactly, and a non-None segmentlist must contain
		a segment which matches exactly).

		It makes little sense to specify both segment and
		segmentlist arguments, but it is not prohibited.

		Bash-style wildcards (*?) are allowed for ifos and description.
		"""
		if exact_match:
			segment_func = lambda e: e.segment == segment
			segmentlist_func = lambda e: e.segment in segmentlist
		else:
			if ifos is not None: ifos = "*" + ifos + "*"
			if description is not None: description = "*" + description + "*"
			segment_func = lambda e: segment.intersects(e.segment)
			segmentlist_func = lambda e: segmentlist.intersects_segment(e.segment)
		
		c = self
		
		if ifos is not None:
			ifos_regexp = re.compile(fnmatch.translate(ifos))
			c = [entry for entry in c if ifos_regexp.match(entry.observatory) is not None]
		
		if description is not None:
			descr_regexp = re.compile(fnmatch.translate(description))
			c = [entry for entry in c if descr_regexp.match(entry.description) is not None]
		
		if segment is not None:
			c = [entry for entry in c if segment_func(entry)]
		
		if segmentlist is not None:
			# must coalesce for intersects_segment() to work
			segmentlist.coalesce()
			c = [entry for entry in c if segmentlist_func(entry)]

		return self.__class__(c)

	def pfnlist(self):
		"""
		Return a list of physical file names.
		"""
		return [entry.path for entry in self]

	def checkfilesexist(self, on_missing="warn"):
		'''
		Runs through the entries of the Cache() object and checks each entry
		if the file which it points to exists or not. If the file does exist then 
		it adds the entry to the Cache() object containing found files, otherwise it
		adds the entry to the Cache() object containing all entries that are missing. 
		It returns both in the follwing order: Cache_Found, Cache_Missed.
		
		Pass on_missing to control how missing files are handled:
		  "warn": print a warning message saying how many files
		          are missing out of the total checked.
		  "error": raise an exception if any are missing
		  "ignore": do nothing
		'''  
		if on_missing not in ("warn", "error", "ignore"):
			raise ValueError("on_missing must be \"warn\", \"error\", or \"ignore\".")
		
		c_found = []
		c_missed = []
		for entry in self:
			if os.path.isfile(entry.path):
				c_found.append(entry)
			else:
				c_missed.append(entry)
		
		if len(c_missed) > 0:
			msg = "%d of %d files in the cache were not found "\
			    "on disk" % (len(c_missed), len(self))
			if on_missing == "warn":
				sys.stderr.write("warning: %s\n" % msg)
			elif on_missing == "error":
				raise ValueError(msg)
			elif on_missing == "ignore":
				pass
			else:
				raise ValueError("Why am I here? "\
				      "Please file a bug report!")
		return self.__class__(c_found), self.__class__(c_missed)

	def __getslice__(self, i, j):
		return self.__class__(super(Cache, self).__getslice__(i, j))
