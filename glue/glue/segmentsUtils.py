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


"""
This module provides additional utilities for use with segments.segmentlist
objects.
"""


import re


from glue import git_version
from glue.lal import CacheEntry
from lal import LIGOTimeGPS
from glue import segments


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                     I/O
#
# =============================================================================
#


#
# A list of file names
#


def fromfilenames(filenames, coltype = int):
	"""
	Return a segmentlist describing the intervals spanned by the files
	whose names are given in the list filenames.  The segmentlist is
	constructed by parsing the file names, and the boundaries of each
	segment are coerced to type coltype.

	The file names are parsed using a generalization of the format
	described in Technical Note LIGO-T010150-00-E, which allows the
	start time and duration appearing in the file name to be
	non-integers.

	NOTE:  the output is a segmentlist as described by the file names;
	if the file names are not in time order, or describe overlaping
	segments, then thusly shall be the output of this function.  It is
	recommended that this function's output be coalesced before use.
	"""
	pattern = re.compile(r"-([\d.]+)-([\d.]+)\.[\w_+#]+\Z")
	l = segments.segmentlist()
	for name in filenames:
		[(s, d)] = pattern.findall(name.strip().rstrip(".gz"))
		s = coltype(s)
		d = coltype(d)
		l.append(segments.segment(s, s + d))
	return l


#
# LAL cache files
#


def fromlalcache(cachefile, coltype = int):
	"""
	Construct a segmentlist representing the times spanned by the files
	identified in the LAL cache contained in the file object file.  The
	segmentlist will be created with segments whose boundaries are of
	type coltype, which should raise ValueError if it cannot convert
	its string argument.

	Example:

	>>> from lal import LIGOTimeGPS
	>>> cache_seglists = fromlalcache(open(filename), coltype = LIGOTimeGPS).coalesce()

	See also:

	glue.lal.CacheEntry
	"""
	return segments.segmentlist(CacheEntry(l, coltype = coltype).segment for l in cachefile)


#
# Segwizard-formated segment list text files
#


def fromsegwizard(file, coltype = int, strict = True):
	"""
	Read a segmentlist from the file object file containing a segwizard
	compatible segment list.  Parsing stops on the first line that
	cannot be parsed (which is consumed).  The segmentlist will be
	created with segment whose boundaries are of type coltype, which
	should raise ValueError if it cannot convert its string argument.
	Two-column, three-column, and four-column segwizard files are
	recognized, but the entire file must be in the same format, which
	is decided by the first parsed line.  If strict is True and the
	file is in three- or four-column format, then each segment's
	duration is checked against that column in the input file.

	NOTE:  the output is a segmentlist as described by the file;  if
	the segments in the input file are not coalesced or out of order,
	then thusly shall be the output of this function.  It is
	recommended that this function's output be coalesced before use.
	"""
	commentpat = re.compile(r"\s*([#;].*)?\Z", re.DOTALL)
	twocolsegpat = re.compile(r"\A\s*([\d.+-eE]+)\s+([\d.+-eE]+)\s*\Z")
	threecolsegpat = re.compile(r"\A\s*([\d.+-eE]+)\s+([\d.+-eE]+)\s+([\d.+-eE]+)\s*\Z")
	fourcolsegpat = re.compile(r"\A\s*([\d]+)\s+([\d.+-eE]+)\s+([\d.+-eE]+)\s+([\d.+-eE]+)\s*\Z")
	format = None
	l = segments.segmentlist()
	for line in file:
		line = commentpat.split(line)[0]
		if not line:
			continue
		try:
			[tokens] = fourcolsegpat.findall(line)
			num = int(tokens[0])
			seg = segments.segment(map(coltype, tokens[1:3]))
			duration = coltype(tokens[3])
			this_line_format = 4
		except ValueError:
			try:
				[tokens] = threecolsegpat.findall(line)
				seg = segments.segment(map(coltype, tokens[0:2]))
				duration = coltype(tokens[2])
				this_line_format = 3
			except ValueError:
				try:
					[tokens] = twocolsegpat.findall(line)
					seg = segments.segment(map(coltype, tokens[0:2]))
					duration = abs(seg)
					this_line_format = 2
				except ValueError:
					break
		if strict:
			if abs(seg) != duration:
				raise ValueError("segment '%s' has incorrect duration" % line)
			if format is None:
				format = this_line_format
			elif format != this_line_format:
				raise ValueError("segment '%s' format mismatch" % line)
		l.append(seg)
	return l


def tosegwizard(file, seglist, header = True, coltype = int):
	"""
	Write the segmentlist seglist to the file object file in a
	segwizard compatible format.  If header is True, then the output
	will begin with a comment line containing column names.  The
	segment boundaries will be coerced to type coltype and then passed
	to str() before output.
	"""
	if header:
		file.write("# seg\tstart    \tstop     \tduration\n")
	for n, seg in enumerate(seglist):
		file.write("%d\t%s\t%s\t%s\n" % (n, str(coltype(seg[0])), str(coltype(seg[1])), str(coltype(abs(seg)))))


#
# TAMA-formated segment list text files
#


def fromtama(file, coltype = LIGOTimeGPS):
	"""
	Read a segmentlist from the file object file containing TAMA
	locked-segments data.  Parsing stops on the first line that cannot
	be parsed (which is consumed).  The segmentlist will be created
	with segments whose boundaries are of type coltype, which should
	raise ValueError if it cannot convert its string argument.

	NOTE:  TAMA locked-segments files contain non-integer start and end
	times, so the default column type is set to LIGOTimeGPS.

	NOTE:  the output is a segmentlist as described by the file;  if
	the segments in the input file are not coalesced or out of order,
	then thusly shall be the output of this function.  It is
	recommended that this function's output be coalesced before use.
	"""
	segmentpat = re.compile(r"\A\s*\S+\s+\S+\s+\S+\s+([\d.+-eE]+)\s+([\d.+-eE]+)")
	l = segments.segmentlist()
	for line in file:
		try:
			[tokens] = segmentpat.findall(line)
			l.append(segments.segment(map(coltype, tokens[0:2])))
		except ValueError:
			break
	return l


#
# Command line or config file strings
#


def from_range_strings(ranges, boundtype = int):
	"""
	Parse a list of ranges expressed as strings in the form "value" or
	"first:last" into an equivalent glue.segments.segmentlist.  In the
	latter case, an empty string for "first" and(or) "last" indicates a
	(semi)infinite range.  A typical use for this function is in
	parsing command line options or entries in configuration files.

	NOTE:  the output is a segmentlist as described by the strings;  if
	the segments in the input file are not coalesced or out of order,
	then thusly shall be the output of this function.  It is
	recommended that this function's output be coalesced before use.

	Example:

	>>> text = "0:10,35,100:"
	>>> from_range_strings(text.split(","))
	[segment(0, 10), segment(35, 35), segment(100, infinity)]
	"""
	# preallocate segmentlist
	segs = segments.segmentlist([None] * len(ranges))

	# iterate over strings
	for i, range in enumerate(ranges):
		parts = range.split(":")
		if len(parts) == 1:
			parts = boundtype(parts[0])
			segs[i] = segments.segment(parts, parts)
			continue
		if len(parts) != 2:
			raise ValueError(range)
		if parts[0] == "":
			parts[0] = segments.NegInfinity
		else:
			parts[0] = boundtype(parts[0])
		if parts[1] == "":
			parts[1] = segments.PosInfinity
		else:
			parts[1] = boundtype(parts[1])
		segs[i] = segments.segment(parts[0], parts[1])

	# success
	return segs


def to_range_strings(seglist):
	"""
	Turn a segment list into a list of range strings as could be parsed
	by from_range_strings().  A typical use for this function is in
	machine-generating configuration files or command lines for other
	programs.

	Example:

	>>> from glue.segments import *
	>>> segs = segmentlist([segment(0, 10), segment(35, 35), segment(100, infinity())])
	>>> ",".join(to_range_strings(segs))
	'0:10,35,100:'
	"""
	# preallocate the string list
	ranges = [None] * len(seglist)

	# iterate over segments
	for i, seg in enumerate(seglist):
		if not seg:
			ranges[i] = str(seg[0])
		elif (seg[0] is segments.NegInfinity) and (seg[1] is segments.PosInfinity):
			ranges[i] = ":"
		elif (seg[0] is segments.NegInfinity) and (seg[1] is not segments.PosInfinity):
			ranges[i] = ":%s" % str(seg[1])
		elif (seg[0] is not segments.NegInfinity) and (seg[1] is segments.PosInfinity):
			ranges[i] = "%s:" % str(seg[0])
		elif (seg[0] is not segments.NegInfinity) and (seg[1] is not segments.PosInfinity):
			ranges[i] = "%s:%s" % (str(seg[0]), str(seg[1]))
		else:
			raise ValueError(seg)

	# success
	return ranges


def segmentlistdict_to_short_string(seglists):
	"""
	Return a string representation of a segmentlistdict object.  Each
	segmentlist in the dictionary is encoded using to_range_strings()
	with "," used to delimit segments.  The keys are converted to
	strings and paired with the string representations of their
	segmentlists using "=" as a delimiter.  Finally the key=value
	strings are combined using "/" to delimit them.

	Example:

	>>> from glue.segments import *
	>>> segs = segmentlistdict({"H1": segmentlist([segment(0, 10), segment(35, 35), segment(100, infinity())]), "L1": segmentlist([segment(5, 15), segment(45, 60)])})
	>>> segmentlistdict_to_short_string(segs)
	'H1=0:10,35,100:/L1=5:15,45:60'

	This function, and its inverse segmentlistdict_from_short_string(),
	are intended to be used to allow small segmentlistdict objects to
	be encoded in command line options and config files.  For large
	segmentlistdict objects or when multiple sets of segmentlists are
	required, the LIGO Light Weight XML encoding available through the
	glue.ligolw library should be used.
	"""
	return "/".join(["%s=%s" % (str(key), ",".join(to_range_strings(value))) for key, value in seglists.items()])


def segmentlistdict_from_short_string(s, boundtype = int):
	"""
	Parse a string representation of a set of named segmentlists into a
	segmentlistdict object.  The string encoding is that generated by
	segmentlistdict_to_short_string().  The optional boundtype argument
	will be passed to from_range_strings() when parsing the segmentlist
	objects from the string.

	Example:

	>>> segmentlistdict_from_short_string("H1=0:10,35,100:/L1=5:15,45:60")
	{'H1': [segment(0, 10), segment(35, 35), segment(100, infinity)], 'L1': [segment(5, 15), segment(45, 60)]}

	This function, and its inverse segmentlistdict_to_short_string(),
	are intended to be used to allow small segmentlistdict objects to
	be encoded in command line options and config files.  For large
	segmentlistdict objects or when multiple sets of segmentlists are
	required, the LIGO Light Weight XML encoding available through the
	glue.ligolw library should be used.
	"""
	d = segments.segmentlistdict()
	for token in s.strip().split("/"):
		key, ranges = token.strip().split("=")
		d[key.strip()] = from_range_strings(ranges.strip().split(","), boundtype = boundtype)
	return d


def from_bitstream(bitstream, start, dt, minlen = 1):
	"""
	Convert consecutive True values in a bit stream (boolean-castable
	iterable) to a stream of segments. Require minlen consecutive True
	samples to comprise a segment.

	Example:

	>>> list(from_bitstream((True, True, False, True, False), 0, 1))
	[segment(0, 2), segment(3, 4)]
	>>> list(from_bitstream([[], [[]], [[]], [], []], 1013968613, 0.125))
	[segment(1013968613.125, 1013968613.375)]
	"""
	bitstream = iter(bitstream)
	i = 0
	while 1:
		if bitstream.next():
			# found start of True block; find the end
			j = i + 1
			try:
				while bitstream.next():
					j += 1
			finally:  # make sure StopIteration doesn't kill final segment
				if j - i >= minlen:
					yield segments.segment(start + i * dt, start + j * dt)
			i = j  # advance to end of block
		i += 1


#
# =============================================================================
#
#                    Pre-defined Segments and Segment Lists
#
# =============================================================================
#


def S2playground(extent):
	"""
	Return a segmentlist identifying the S2 playground times within the
	interval defined by the segment extent.

	Example:

	>>> from glue import segments
	>>> S2playground(segments.segment(874000000, 874010000))
	[segment(874000013, 874000613), segment(874006383, 874006983)]
	"""
	lo = int(extent[0])
	lo -= (lo - 729273613) % 6370
	hi = int(extent[1]) + 1
	return segments.segmentlist(segments.segment(t, t + 600) for t in range(lo, hi, 6370)) & segments.segmentlist([extent])


def segmentlist_range(start, stop, period):
	"""
	Analogous to Python's range() builtin, this generator yields a
	sequence of continuous adjacent segments each of length "period"
	with the first starting at "start" and the last ending not after
	"stop".  Note that the segments generated do not form a coalesced
	list (they are not disjoint).  start, stop, and period can be any
	objects which support basic arithmetic operations.

	Example:

	>>> from glue.segments import *
	>>> segmentlist(segmentlist_range(0, 15, 5))
	[segment(0, 5), segment(5, 10), segment(10, 15)]
	>>> segmentlist(segmentlist_range('', 'xxx', 'x'))
	[segment('', 'x'), segment('x', 'xx'), segment('xx', 'xxx')]
	"""
	n = 1
	b = start
	while True:
		a, b = b, start + n * period
		if b > stop:
			break
		yield segments.segment(a, b)
		n += 1


#
# =============================================================================
#
#                         Extra Manipulation Routines
#
# =============================================================================
#


def Fold(seglist1, seglist2):
	"""
	An iterator that generates the results of taking the intersection
	of seglist1 with each segment in seglist2 in turn.  In each result,
	the segment start and stop values are adjusted to be with respect
	to the start of the corresponding segment in seglist2.  See also
	the segmentlist_range() function.

	This has use in applications that wish to convert ranges of values
	to ranges relative to epoch boundaries.  Below, a list of time
	intervals in hours is converted to a sequence of daily interval
	lists with times relative to midnight.

	Example:

	>>> from glue.segments import *
	>>> x = segmentlist([segment(0, 13), segment(14, 20), segment(22, 36)])
	>>> for y in Fold(x, segmentlist_range(0, 48, 24)): print y
	...
	[segment(0, 13), segment(14, 20), segment(22, 24)]
	[segment(0, 12)]
	"""
	for seg in seglist2:
		yield (seglist1 & segments.segmentlist([seg])).shift(-seg[0])


def vote(seglists, n):
	"""
	Given a sequence of segmentlists, returns the intervals during
	which at least n of them intersect.  The input segmentlists must be
	coalesced, the output is coalesced.

	Example:

	>>> from glue.segments import *
	>>> w = segmentlist([segment(0, 15)])
	>>> x = segmentlist([segment(5, 20)])
	>>> y = segmentlist([segment(10, 25)])
	>>> z = segmentlist([segment(15, 30)])
	>>> vote((w, x, y, z), 3)
	[segment(10, 20)]

	The sequence of segmentlists is only iterated over once, and the
	segmentlists within it are only iterated over once;  they can all
	be generators.  If there are a total of N segments in M segment
	lists and the final result has L segments the algorithm is O(N M) +
	O(L).
	"""
	# check for no-op

	if n < 1:
		return segments.segmentlist()

	# digest the segmentlists into an ordered sequence of off-on and
	# on-off transitions with the vote count for each transition
	# FIXME:  this generator is declared locally for now, is it useful
	# as a stand-alone generator?

	def pop_min(l):
		# remove and return the smallest value from a list
		val = min(l)
		for i in xrange(len(l) - 1, -1, -1):
			if l[i] is val:
				return l.pop(i)
		assert False	# cannot get here

	def vote_generator(seglists):
		queue = []
		for seglist in seglists:
			segiter = iter(seglist)
			try:
				seg = segiter.next()
			except StopIteration:
				continue
			# put them in so that the smallest boundary is
			# closest to the end of the list
			queue.append((seg[1], -1, segiter))
			queue.append((seg[0], +1, None))
		if not queue:
			return
		queue.sort(reverse = True)
		bound = queue[-1][0]
		votes = 0
		while queue:
			this_bound, delta, segiter = pop_min(queue)
			if this_bound == bound:
				votes += delta
			else:
				yield bound, votes
				bound = this_bound
				votes = delta
			if segiter is not None:
				try:
					seg = segiter.next()
				except StopIteration:
					continue
				queue.append((seg[1], -1, segiter))
				queue.append((seg[0], +1, None))
		yield bound, votes

	# compute the cumulative sum of votes, and assemble a segmentlist
	# from the intervals when the vote count is equal to or greater
	# than n

	result = segments.segmentlist()
	votes = 0
	for bound, delta in vote_generator(seglists):
		if delta > 0 and n - delta <= votes < n:
			start = bound
		elif delta < 0 and n <= votes < n - delta:
			result.append(segments.segment(start, bound))
			del start	# detect stops that aren't preceded by starts
		votes += delta
	assert votes == 0	# detect failed cumulative sum

	return result
