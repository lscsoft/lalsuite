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


import itertools


from glue import offsetvector


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
from git_version import date as __date__
from git_version import version as __version__


#
# =============================================================================
#
#                             Command Line Parsing
#
# =============================================================================
#


def parse_slidespec(slidespec):
	"""
	Accepts a string in the format
	instrument=first:last:step[,first:last:step]...  and returns the
	tuple (instrument, [offset1, offset2, ....]) where the offsets are
	the sorted list of unique numbers described by the ranges.  A range
	with a positive step describes the offsets (first + n * step) where
	n is an integer such that first <= offset <= last.  A range with a
	negative step describes the offsets (first + n * step) where n is
	an integer such that last <= offset <= first.

	Example:

	>>> parse_slidespec("H1=-5:+5:0.5")
	('H1', [-5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])
	>>> parse_slidespec("H1=+5:-5:-0.5")
	('H1', [-5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])
	>>> parse_slidespec("H1=-5:+5:8")
	('H1', [-5.0, 3.0])
	"""
	try:
		instrument, rangespec = slidespec.split("=")
	except ValueError:
		raise ValueError("cannot parse time slide '%s'" % slidespec)
	offsets = set()
	for rng in [s.strip() for s in rangespec.split(",")]:
		try:
			first, last, step = map(float, rng.split(":"))
		except ValueError:
			raise ValueError("malformed range '%s' in '%s'" % (rng, rangespec))
		if step == 0:
			if first != last:
				raise ValueError("divide by zero in range '%s'" % rng)
			offsets.add(first)
			continue
		elif last > first if step < 0 else last < first:
			raise ValueError("step has wrong sign in range '%s'" % rng)

		for i in itertools.count():
			x = first + i * step
			if not min(first, last) <= x <= max(first, last):
				break
			offsets.add(x)
	return instrument.strip(), sorted(offsets)


def parse_slides(slides):
	"""
	Accepts a list of strings of the format understood by
	parse_slidespec(), and returns a dictionary mapping instrument
	names to sorted lists of unique offsets.

	Example:

	>>> parse_slides(["H1=-1:+1:+1", "H2=-1:+1:+1", "L1=0:0:0"])
	{'H2': [-1.0, 0.0, 1.0], 'H1': [-1.0, 0.0, 1.0], 'L1': [0.0]}
	"""
	d = {}
	# store the offsets for each instrument as sets to uniquify the
	# numbers
	for slidespec in slides:
		instrument, offsets = parse_slidespec(slidespec)
		try:
			d[instrument] |= set(offsets)
		except KeyError:
			d[instrument] = set(offsets)
	# convert offsets back to sorted lists
	return dict((instrument, sorted(offsets)) for instrument, offsets in d.items())


def parse_inspiral_num_slides_slidespec(slidespec):
	"""
	Accepts a string in the format
	count:instrument=offset[,instrument=offset...] and returns the
	tuple (count, {instrument: offset, ...})

	Example:

	>>> parse_inspiral_num_slides_slidespec("3:H1=0,H2=5,L1=10")
	(3, offsetvector({'H2': 5.0, 'H1': 0.0, 'L1': 10.0}))
	"""
	count, offsets = slidespec.strip().split(":")
	offsetvect = offsetvector.offsetvector((instrument.strip(), float(offset)) for instrument, offset in (token.strip().split("=") for token in offsets.strip().split(",")))
	return int(count), offsetvect


#
# =============================================================================
#
#                              Build Time Slides
#
# =============================================================================
#


def SlidesIter(slides):
	"""
	Accepts a dictionary mapping instrument --> list-of-offsets (for
	example, as returned by parse_slides()), and iterates over the
	cartesian (outer) product of the offset lists, yielding all
	possible N-way instrument --> offset mappings.

	Example:

	>>> slides = {"H1": [-1, 0, 1], "H2": [-1, 0, 1], "L1": [0]}
	>>> list(SlidesIter(slides))
	[offsetvector({'H2': -1, 'H1': -1, 'L1': 0}), offsetvector({'H2': -1, 'H1': 0, 'L1': 0}), offsetvector({'H2': -1, 'H1': 1, 'L1': 0}), offsetvector({'H2': 0, 'H1': -1, 'L1': 0}), offsetvector({'H2': 0, 'H1': 0, 'L1': 0}), offsetvector({'H2': 0, 'H1': 1, 'L1': 0}), offsetvector({'H2': 1, 'H1': -1, 'L1': 0}), offsetvector({'H2': 1, 'H1': 0, 'L1': 0}), offsetvector({'H2': 1, 'H1': 1, 'L1': 0})]
	"""
	if not slides:
		# things get a little odd in the even that no
		# instrument/offset-list pairs are given. instead of
		# yielding an empty sequence, itertools.product(*()) yields
		# a sequence containing a single empty tuple, so instead of
		# yielding no offsetvectors this function yields one empty
		# one.  that's not what calling codes generally expect the
		# response to be so we trap the case and return an empty
		# sequence
		return
	instruments = slides.keys()
	for slide in itertools.product(*slides.values()):
		yield offsetvector.offsetvector(zip(instruments, slide))


def Inspiral_Num_Slides_Iter(count, offsets):
	"""
	This generator yields a sequence of time slide dictionaries in the
	style of lalapps_thinca's time slides.  Each resulting dictionary
	maps instrument to offset.  The input is a count of time slides (an
	integer), and a dictionary mapping instrument to offset.  The
	output dictionaries describe time slides that are integer multiples
	of the input time shifts.

	Example:

	>>> list(Inspiral_Num_Slides_Iter(3, {"H1": 0.0, "H2": 5.0,"L1": 10.0}))
	[offsetvector({'H2': -15.0, 'H1': -0.0, 'L1': -30.0}), offsetvector({'H2': -10.0, 'H1': -0.0, 'L1': -20.0}), offsetvector({'H2': -5.0, 'H1': -0.0, 'L1': -10.0}), offsetvector({'H2': 0.0, 'H1': 0.0, 'L1': 0.0}), offsetvector({'H2': 5.0, 'H1': 0.0, 'L1': 10.0}), offsetvector({'H2': 10.0, 'H1': 0.0, 'L1': 20.0}), offsetvector({'H2': 15.0, 'H1': 0.0, 'L1': 30.0})]

	The output time slides are all integer multiples of the input time
	shift vector in the range [-count, +count], and are returned in
	increasing order of mupltiplier.
	"""
	offsets = offsets.items()
	for n in range(-count, +count + 1):
		yield offsetvector.offsetvector((instrument, offset * n) for instrument, offset in offsets)
