#  Copyright (C) 2013 Chris Pankow
#
#  This program is free software; ynu can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Fnundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with with program; see the file COPYING. If not, write to the
#  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#  MA  02111-1307  USA

## \defgroup laldetchar_py_dqsegs dqsegs
## \ingroup laldetchar_python
"""Utilities to transform data samples into segment information"""
#
# ### Synopsis ###
#
# ~~~
# from laldetchar import dqsegs
# ~~~
# \author Chris Pankow (<chris.pankow@ligo.org>)

from glue.segmentsUtils import from_bitstream
#from laldetchar import git_version
import git_version

__author__  = "Chris Pankow <chris.pankow@ligo.org>"
__version__ = git_version.id
__date__    = git_version.date

## \addtogroup laldetchar_py_dqsegs
#@{

def threshold_data_to_seglist(data, start, dt, min_threshold=None, max_threshold=None, invert=False):
	"""
	Apply min and max threshold to data and parse the result into a segment list. If invert, then invert, tautologies not withstanding.
	"""
	if min_threshold is None:
		min_threshold = -float("inf")
	if max_threshold is None:
		max_threshold = float("inf")
	if invert:
		return from_bitstream([min_threshold > d or d > max_threshold for d in data], start, dt)
	else:
		return from_bitstream([min_threshold < d and d < max_threshold for d in data], start, dt)

def equality_data_to_seglist(data, start, dt, equality, invert=False):
	"""
	Apply equality test threshold to data and parse the result into a segment list.
	"""
	# FIXME: numpy operations?
	if invert:
		return from_bitstream([d != equality for d in data], start, dt)
	else:
		return from_bitstream([d == equality for d in data], start, dt)

def mask_data_to_seglist(data, start, dt, mask_on=0x1, mask_off=0x0):
	"""
	Apply bitmask to data and parse the result into a segment list.
	"""
	# FIXME: numpy operations?
	return from_bitstream([(d & mask_on) & (not(d & mask_off)) for d in data], start, dt)

##@}
