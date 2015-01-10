# Copyright (C) 2013 Chris Pankow
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

## \addtogroup laldetchar_py_hveto
"""Convenience utitilies for accessing the SWIG wrapped glib functions.

Potential developers, please note: this package should eventually be supplanted by the gobject bindings for glib if they ever become available. """
#
# Synopsis
# ~~~
# from laldetchar.hveto import glib_utils
# ~~~
# \author Chris Pankow (<chris.pankow@ligo.org>)

import sys

try:
	from glue.segments import segment, segmentlist, segmentlistdict
	from glue.ligolw import ligolw, utils, lsctables, ilwd, table
except ImportError:
	print >>sys.stderr, "Glue installation is required to import this module."
	sys.exit(-1)

import laldetchar
import lal
# FIXME: Will this be needed with Karl's bindings?
# Yes... lalburst doesn't import the wrapped version of SnglBurst, so trying to
# get its attributes will fail until they are defined with this
# But that's been fixed, so we should try it without
import lalmetaio

from laldetchar.hveto import lalburst_sb_to_glue_sb
from laldetchar import git_version

__author__ = "Chris Pankow <chris.pankow@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date

## \addtogroup laldetchar_py_hveto_glib_utils
#@{

# Functions to convert various SWIG wrapped LAL functions to and from their
# (more convenient) equivalents in glue.

def seq_to_sbtable(sblist, get_col=None):
	"""
	Convert a GSequence of SnglBurst to a SnglBurstTable.
	"""
	desired_columns = ["peak_time_ns", "start_time_ns", "channel",
	  "process_id", "ifo", "peak_time", "start_time",
	  "duration", "central_freq", "search",
	  "bandwidth", "event_id", "snr"]
	sbtable = lsctables.New(lsctables.SnglBurstTable, desired_columns)

	sblistitr = laldetchar.GSequenceBegin(sblist)

	col_ar = []
	while sblistitr is not None:
		next_sb = laldetchar.GetGSeqSnglBurst(sblistitr)
		if next_sb is None:
			break
		if get_col is not None:
			col_ar.append(getattr(next_sb, get_col))
		sbtable.append(lalburst_sb_to_glue_sb(next_sb, desired_columns))
		if not laldetchar.GSequenceNext(sblistitr):
                        sblistitr = None

	if get_col is not None:
		return col_ar
	return sbtable

def lalseg_from_seglist(segl):
	"""
	Convert a glue.segments segmentlist to a LALSegList structure.
	"""
	lalsegl = lal.SegList()
	lal.SegListInit(lalsegl)
	for i, seg in enumerate(segl):
		lalseg = lal.Seg()
		lal.SegSet(lalseg,
			lal.LIGOTimeGPS(float(seg[0])),
			lal.LIGOTimeGPS(float(seg[1])),
			i
		)
		lal.SegListAppend( lalsegl, lalseg )
	return lalsegl

def seglist_from_lalseg(lalsegl):
	"""
	Convert a LALSegList structure to a glue.segmentlist.
	"""
	segl = segmentlist()
	for i in range(lalsegl.length):
		lalseg = lal.SegListGet(lalsegl, i)
		segl.append(segment(lalseg.start,lalseg.end))
	return segl

# Functions to convert various SWIG wrapped glib and LAL structures to their
# (more convenient) equivalents in python.

def ghash_key_in_dict(ghash, dkey):
	"""
	Check if dkey is in ghash.
	"""
	i = 0
	while True:
		key = laldetchar.GetGHashTableKey(ghash, i)
		if key is None:
			break
		elif key == dkey:
			return True
		i += 1
	return False


def ghash_to_dict(ghash):
	"""
	Convert a GHashTable into a python dictionary. Note that all the GHashTables in use by laldetchar map a string to a fixed type (e.g. double or string), so the dtype argument should match that.
	"""
	i = 0
	keys = []
	while True:
		key = laldetchar.GetGHashTableKey(ghash, i)
		if key is None:
			break
		i += 1
		keys.append(key)

	dtype = laldetchar.GetGHashTableType(ghash)
	if dtype == laldetchar.LALGTYPE_INT:
		rfunc = laldetchar.GetGHashTblInt
	elif dtype == laldetchar.LALGTYPE_DBL:
		rfunc = laldetchar.GetGHashTblDbl
	elif dtype == laldetchar.LALGTYPE_STR:
		rfunc = laldetchar.GetGHashTblStr
	else:
		raise TypeError("Don't have right function to retrieve key for LALGType=%i" % dtype)
	return dict([ (key, rfunc(ghash, key)) for key in keys ])

def get_chan_list( gseq ):
	"""
	Get a listing of the channels within a trigger sequence.
	"""
	# This is necessarily complicated because there is no "set" type in glib
	# thus the closest thing is a hash. This retrieves the hash, then
	# converts its keys into a proper set.
	return set(ghash_to_dict( laldetchar.GetChannelList(gseq) ).keys())

##@}
