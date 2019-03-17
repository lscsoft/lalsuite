# Copyright (C) 2006-2010,2012--2013  Kipp Cannon
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
LIGO Light-Weight XML coincidence analysis front end.
"""


from __future__ import print_function


import math
import sys


from lal import LIGOTimeGPS
from lal.utils import CacheEntry
from ligo import segments


from . import offsetvector
from . import packing


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
from .git_version import date as __date__
from .git_version import version as __version__


#
# =============================================================================
#
#                                    Input
#
# =============================================================================
#


def load_cache(filename, verbose = False):
	"""
	Parse a LAL cache file named filename into a list of
	lal.utils.CacheEntry objects.  If filename is None then input is
	taken from stdin.
	"""
	if verbose:
		print("reading %s ..." % (filename or "stdin"), file=sys.stderr)
	if filename is not None:
		f = open(filename)
	else:
		f = sys.stdin
	return [CacheEntry(line) for line in f]


def cache_to_seglistdict(cache):
	"""
	Construct a coalesced segmentlistdict object from a list of
	lal.utils.CacheEntry objects.
	"""
	s = segments.segmentlistdict()
	for c in cache:
		s |= c.segmentlistdict
	return s


#
# =============================================================================
#
#                             Performance Helpers
#
# =============================================================================
#


def segmentlistdict_normalize(seglistdict, origin):
	"""
	Convert the times in a segmentlist dictionary to floats relative to
	origin.  The purpose is to allow segment lists stored as
	LIGOTimeGPS times to be manipulated more quickly without loss of
	precision.  The modification is done in place.
	"""
	for seglist in seglistdict.itervalues():
		seglist[:] = (segments.segment(float(seg[0] - origin), float(seg[1] - origin)) for seg in seglist)


def get_coincident_segmentlistdict(seglistdict, offset_vectors):
	"""
	Compute the segments for which data is required in order to perform
	a complete coincidence analysis given the segments for which data
	is available and the list of offset vectors to be applied to the
	data during the coincidence analysis.

	seglistdict is a segmentlistdict object defining the instruments
	and times for which data is available.  offset_vectors is a list of
	offset vectors to be applied to the data --- dictionaries of
	instrument/offset pairs.

	The offset vectors in offset_vectors are applied to the input
	segments one by one and the interesection of the shifted segments
	is computed.  The segments surviving the intersection are unshifted
	to their original positions and stored.  The return value is the
	union of the results of this operation.

	In all cases all pair-wise intersections are computed, that is if
	an offset vector lists three instruments then this function returns
	the times when any two of those isntruments are on, including times
	when all three are on.

	For example, let us say that "input" is a segmentlistdict object
	containing segment lists for three instruments, "H1", "H2" and
	"L1".  And let us say that "slides" is a list of dictionaries, and
	is equal to [{"H1":0, "H2":0, "L1":0}, {"H1":0, "H2":10}].  Then if

	output = get_coincident_segmentlistdict(input, slides)

	output will contain, for each of the three instruments, the
	segments (or parts thereof) from the original lists that are
	required in order to perform a triple- and double-coincident
	analyses at zero lag with the three instruments, *and* a
	double-coincident analysis between H1 and H2 with H2 offset by 10
	seconds.

	The segmentlistdict object returned by this function has its
	offsets set to those of the input segmentlistdict.
	"""
	# don't modify original
	seglistdict = seglistdict.copy()
	all_instruments = set(seglistdict)

	# save original offsets
	origoffsets = dict(seglistdict.offsets)

	# compute result
	coincseglists = segments.segmentlistdict()
	for offset_vector in offsetvector.component_offsetvectors(offset_vectors, 2):
		if set(offset_vector).issubset(all_instruments):
			seglistdict.offsets.update(offset_vector)
			intersection = seglistdict.extract_common(offset_vector.keys())
			intersection.offsets.clear()
			coincseglists |= intersection

	# restore original offsets
	coincseglists.offsets.update(origoffsets)

	# done
	return coincseglists


def segmentlistdict_unnormalize(seglistdict, origin):
	"""
	The opposite of segmentlistdict_normalize(), restores the times in
	a segmentlist dictionary to absolute times.  The modification is
	done in place.
	"""
	for seglist in seglistdict.itervalues():
		seglist[:] = (segments.segment(origin + seg[0], origin + seg[1]) for seg in seglist)


#
# =============================================================================
#
#                             Output Cache Packing
#
# =============================================================================
#


class LALCacheBin(packing.Bin):
	"""
	Subclass of the packing.Bin class representing a LAL file cache.
	The files contained in the bin are available in the .objects
	attribute, which is a list of lal.utils.CacheEntry objects.  The
	.size attribute holds a ligo.segments.segmentlistdict object giving
	the times spanned by the files in the bin.  The .extent attribute
	holds the result of running .extent_all() on the .size attribute.
	"""
	def __init__(self):
		packing.Bin.__init__(self)
		self.size = segments.segmentlistdict()
		self.extent = None

	def add(self, cache_entry):
		packing.Bin.add(self, cache_entry, cache_entry.segmentlistdict)
		self.extent = self.size.extent_all()
		return self

	def __iadd__(self, *args):
		packing.Bin.__iadd__(self, *args)
		self.extent = self.size.extent_all()
		return self

	def __cmp__(self, other):
		return cmp(self.extent, other.extent)

	def __str__(self):
		return "\n".join(map(str, self.objects))


class CafePacker(packing.Packer):
	"""
	Packing algorithm implementing the ligolw_cafe file list packing
	algorithm.
	"""
	def set_offset_vectors(self, offset_vectors):
		"""
		Set the list of offset vectors to be considered when
		deciding the bins in which each file belongs.  Must be
		called before packing any files.  The input is a list of
		dictionaries, each mapping instruments to offsets.
		"""
		#
		# sort the offset vectors to reduce the number of
		# arithmetic operations performed while applying them
		#

		self.offset_vectors = list(offset_vectors)
		self.offset_vectors.sort(key = lambda offset_vector: sorted(offset_vector.items()))

		#
		# determine the largest gap that can conceivably be closed
		# by the time slides
		#

		min_offset = min(min(offset_vector.values()) for offset_vector in offset_vectors)
		max_offset = max(max(offset_vector.values()) for offset_vector in offset_vectors)
		self.max_gap = max_offset - min_offset
		assert self.max_gap >= 0

	def pack(self, cache_entry):
		"""
		Find all bins in which this lal.utils.CacheEntry instance
		belongs, merge them, and add this cache entry to the
		result.  Create a new bin for this cache entry if it does
		not belong in any of the existing bins.

		The cache entry "belongs" in a bin if after each of the
		preset offset vectors (see the .set_offset_vectors()
		method) is applied to both the contents of a bin and the
		cache entry, any of the segment lists of the bin and cache
		entry are found to intersect.  When checking for
		intersection, only the segment lists whose instrument names
		are listed in the offset vector are compared.
		"""
		#
		# add the cache entry to a new bin by itself
		#

		new = LALCacheBin()
		new.add(cache_entry)

		#
		# assemble a list of bins in which the cache entry belongs.
		# iterate over existing bins backwards so that we record
		# the indeces of matching bins in descending order.  bail
		# out when we find a bin that precedes the new one
		#

		matching_bins = []
		for n in range(len(self.bins) - 1, -1, -1):
			bin = self.bins[n]
			if bin.extent[1] < new.extent[0] - self.max_gap:
				break
			for offset_vector in self.offset_vectors:
				new.size.offsets.update(offset_vector)
				bin.size.offsets.update(offset_vector)
				if bin.size.is_coincident(new.size, keys = offset_vector.keys()):
					matching_bins.append(n)
					break
			bin.size.offsets.clear()
		new.size.offsets.clear()

		#
		# add new cache entry to bins
		#

		if not matching_bins:
			#
			# no existing bins match, add a new one
			#

			self.bins.append(new)
		else:
			#
			# put cache entry into earliest bin that was found
			# to match.  if cache entry belongs in more than
			# one bin, merge them.  note that the matching bin
			# indexes are given in descending order so the last
			# is the earliest bin, and after that popping them
			# in order does not affet the indexes of the
			# remaining, matching, bins.
			#

			dest = self.bins[matching_bins.pop(-1)]
			dest += new
			for n in matching_bins:
				dest += self.bins.pop(n)

		#
		# time-order the bins so the bail-out above works next time
		# this method is called
		#

		self.bins.sort()


def split_bins(cafepacker, extentlimit, verbose = False):
	"""
	Split bins in CafePacker so that each bin has an extent no longer
	than extentlimit.
	"""

	#
	# loop over all bins in cafepacker.bins.  loop is backwards because
	# list grows in size as bins are split
	#

	for idx in range(len(cafepacker.bins) - 1, -1, -1):
		#
		# retrieve bin
		#

		origbin = cafepacker.bins[idx]

		#
		# how many pieces?  if bin doesn't need splitting move to
		# next
		#

		n = int(math.ceil(float(abs(origbin.extent)) / extentlimit))
		if n <= 1:
			continue

		#
		# calculate the times of the splits, and then build
		# segmentlistdicts for clipping.
		#

		extents = [origbin.extent[0]] + [LIGOTimeGPS(origbin.extent[0] + i * float(abs(origbin.extent)) / n) for i in range(1, n)] + [origbin.extent[1]]
		if verbose:
			print("\tsplitting cache spanning %s at %s" % (str(origbin.extent), ", ".join(str(extent) for extent in extents[1:-1])), file=sys.stderr)
		extents = [segments.segment(*bounds) for bounds in zip(extents[:-1], extents[1:])]

		#
		# build new bins, pack objects from origbin into new bins
		#

		newbins = []
		for extent in extents:
			#
			# append new bin
			#

			newbins.append(LALCacheBin())

			#
			# test each cache entry in original bin
			#

			extent_plus_max_gap = extent.protract(cafepacker.max_gap)
			for cache_entry in origbin.objects:
				#
				# quick check of gap
				#

				if cache_entry.segment.disjoint(extent_plus_max_gap):
					continue

				#
				# apply each offset vector
				#

				cache_entry_segs = cache_entry.segmentlistdict
				for offset_vector in cafepacker.offset_vectors:
					cache_entry_segs.offsets.update(offset_vector)

					#
					# test against bin
					#

					if cache_entry_segs.intersects_segment(extent):
						#
						# object is coicident with
						# bin
						#

						newbins[-1].add(cache_entry)
						break

			#
			# override the bin's extent
			#

			newbins[-1].extent = extent

		#
		# replace original bin with split bins.
		#

		cafepacker.bins[idx:idx+1] = newbins

	#
	# done
	#


#
# =============================================================================
#
#                                    Output
#
# =============================================================================
#


def write_caches(base, bins, instruments = None, verbose = False):
	filenames = []
	if len(bins):
		pattern = "%%s%%0%dd.cache" % int(math.log10(len(bins)) + 1)
	for n, bin in enumerate(bins):
		filename = pattern % (base, n)
		filenames.append(filename)
		if verbose:
			print("writing %s ..." % filename, file=sys.stderr)
		f = open(filename, "w")
		for cacheentry in bin.objects:
			if instruments is None or (instruments & set(cacheentry.segmentlistdict)):
				print(str(cacheentry), file=f)
	return filenames


def write_single_instrument_caches(base, bins, instruments, verbose = False):
	for instrument in instruments:
		write_caches("%s%s_" % (base, instrument), bins, set([instrument]), verbose)


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#


def ligolw_cafe(cache, offset_vectors, verbose = False, extentlimit = None):
	"""
	Transform a LAL cache into a list of caches each of whose contents
	can be subjected to a coincidence analysis independently of the
	contents of the other caches, assuming the coincidence analyses
	will involve the application of the given offset vectors.

	cache is a sequence (e.g., list, tuple, etc.) of
	lal.utils.CacheEntry objects.  offset_vectors is a sequence of
	instrument/offset dictionaries describing the offset vectors to
	consider.  Set verbose to True for verbosity.

	The output is a two-element tuple.  The first element is a
	ligo.segments.segmentlistdict object describing the times for which
	coincident data is available (derived from the segment metadata of
	the input cache).  The second element is a list of LALCacheBin
	objects, providing the file groups.
	"""
	#
	# Construct a segment list dictionary from the cache
	#

	if verbose:
		print("computing segment list ...", file=sys.stderr)
	seglists = cache_to_seglistdict(cache)

	#
	# For each instrument compute the times for which it will (could)
	# contribute to a coincidence analysis.
	#

	epoch = min([min(seg[0] for seg in seglist) for seglist in seglists.values() if seglist] or [None])
	segmentlistdict_normalize(seglists, epoch)
	seglists = get_coincident_segmentlistdict(seglists, offset_vectors)
	segmentlistdict_unnormalize(seglists, epoch)

	#
	# Remove files that will not participate in a coincidence.  Take
	# care not to modify the calling code's data.  Note that because we
	# have established that this segment list describes exactly the
	# times spanned by the input files that are coincident under at
	# least one time slide, a file participates in a multi-instrument
	# coincidence if and only if it intersects these times.
	#

	if verbose:
		print("filtering input cache ...", file=sys.stderr)
	cache = [c for c in cache if seglists.intersects_all(c.segmentlistdict)]

	#
	# Optimization: adding files to bins in time order keeps the number
	# of bins from growing larger than needed.
	#

	if verbose:
		print("sorting input cache ...", file=sys.stderr)
	cache.sort(key = lambda x: x.segment)

	#
	# Pack cache entries into output caches.  Having reduced the file
	# list to just those that participate in coincidences, it only
	# remains to determine which other files each must be grouped with.
	#

	outputcaches = []
	packer = CafePacker(outputcaches)
	packer.set_offset_vectors(offset_vectors)
	if verbose:
		print("packing files (considering %s offset vectors) ..." % len(offset_vectors), file=sys.stderr)
	for n, cacheentry in enumerate(cache):
		if verbose and not n % 13:
			print("\t%.1f%%\t(%d files, %d caches)\r" % (100.0 * n / len(cache), n + 1, len(outputcaches)), end=' ', file=sys.stderr)
		packer.pack(cacheentry)
	if verbose:
		print("\t100.0%%\t(%d files, %d caches)" % (len(cache), len(outputcaches)), file=sys.stderr)

	#
	# Split caches with extent more than extentlimit
	#

	if extentlimit is not None:
		if verbose:
			print("splitting caches with extent greater than %g s ..." % extentlimit, file=sys.stderr)
		split_bins(packer, extentlimit, verbose = verbose)
		if verbose:
			print("\t\t(%d files, %d caches)" % (len(cache), len(outputcaches)), file=sys.stderr)

	#
	# Sort output caches
	#

	if verbose:
		print("sorting output caches ...", file=sys.stderr)
	for cache in outputcaches:
		cache.objects.sort()

	#
	# Done.
	#

	return seglists, outputcaches
