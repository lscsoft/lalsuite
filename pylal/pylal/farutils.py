# Copyright (C) 2010 Chad Hanna
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

from glue import iterutils
from glue import segments
from glue.ligolw import lsctables
from glue.ligolw import dbtables
from glue.ligolw.utils import search_summary as ligolw_search_summary
from glue.ligolw.utils import segments as ligolw_segments
from pylal import SnglBurstUtils
from pylal import db_thinca_rings

# get choices from a set (useful for on/off ifos)
def detector_combos( instruments ):
	out = []
	instruments = tuple(instruments)
	for i in range(1, len(instruments)):
		X = list(iterutils.choices(instruments, i))
		Y = list(iterutils.choices(instruments, len(instruments) - i))
		Y.reverse()
		out.extend(zip(X,Y)) #out.extend(zip(X, Y))
	instruments = list(instruments)
	instruments.sort()
	instruments = tuple(instruments)
	out.append((instruments, ()))
	return out

def background_livetime_nonring_by_slide(connection, seglists, veto_segments=None, coinc_segments=None, verbose = False):
	# get the segment lists and live time
	# FIXME veto segments not handled yet
	seglists = seglists.copy()
	if veto_segments is not None:
		seglists -= veto_segments

	zero_lag_time_slides, background_time_slides = SnglBurstUtils.get_time_slides(connection)
	instruments = frozenset(seglists.keys())
	background_livetime = {}
	for on_inst, off_inst in detector_combos(list(instruments)):
		on_inst = frozenset(on_inst)
		off_inst = frozenset(off_inst)
		key = on_inst
		old_offsets = seglists.offsets.copy()
		background_livetime.setdefault(key, {})
		for id, time_slide in background_time_slides.items():
			seglists.offsets.update(time_slide)
			segs=seglists.intersection(list(on_inst))-seglists.union(list(off_inst))
			if coinc_segments is not None:
				segs &= coinc_segments
			tskey = frozenset(time_slide.items())
			background_livetime[key].setdefault(tskey,0)
			background_livetime[key][tskey] += float(abs(segs))
		seglists.offsets.update(old_offsets)
	return background_livetime

def background_livetime_ring_by_slide(connection, live_time_program, seglists, veto_segments, verbose = False):
	background_livetime = {}
	instruments = frozenset(seglists.keys())
	offset_vectors = db_thinca_rings.get_background_offset_vectors(connection)
	# first work out time slide live time
	for on_instruments, livetimes in db_thinca_rings.get_thinca_livetimes(db_thinca_rings.get_thinca_rings_by_available_instruments(connection, program_name = live_time_program), veto_segments, offset_vectors, verbose = verbose).items():
		on_instruments = frozenset(on_instruments)#lsctables.ifos_from_instrument_set(on_instruments)
		for offset, lt in zip(offset_vectors,livetimes):
			background_livetime.setdefault(on_instruments,{})
			key = frozenset(offset.items())
			background_livetime[on_instruments].setdefault(key, 0)
			background_livetime[on_instruments][key] += lt

	return background_livetime

def add_background_livetime(connection, live_time_program, seglists, veto_segments, coinc_segments=None, verbose=False):
	if live_time_program == "thinca": lt = background_livetime_ring_by_slide(connection, live_time_program, seglists, veto_segments, verbose)
	if live_time_program == "gstlal_inspiral": lt = background_livetime_nonring_by_slide(connection, seglists, veto_segments, coinc_segments, verbose)
	if live_time_program == "lalapps_ring": lt = background_livetime_nonring_by_slide(connection, seglists, veto_segments, coinc_segments, verbose)
	out = {}
	for k, v in lt.items():
		out.setdefault(k,0)
		out[k] += sum(v.values())
	return out

def playground_nonplayground_livetime(seglists, playground_segs=None, verbose=False):
	playground_livetime = {}
	nonplayground_livetime = {}
	instruments = frozenset(seglists.keys())
	for on_inst, off_inst in detector_combos(list(instruments)):
		on_inst = frozenset(on_inst)
		off_inst = frozenset(off_inst)
		key = lsctables.ifos_from_instrument_set(on_inst)
		selected_segs = seglists.intersection(list(on_inst))-seglists.union(list(off_inst))
		if playground_segs:
			playground_livetime[on_inst] = float(abs(selected_segs & playground_segs))
			nonplayground_livetime[on_inst] = float(abs(selected_segs - playground_segs))
		else:
			playground_livetime[on_inst] = 0
			nonplayground_livetime[on_inst] = float(abs(selected_segs))

	return playground_livetime, nonplayground_livetime

def get_veto_segments(connection, program_name, xmldoc=None, veto_segments_name=None):
	veto_segments = segments.segmentlistdict()
	#FIXME only handles thinca case
	if not veto_segments_name: return veto_segments
	if program_name == "thinca": veto_segments = db_thinca_rings.get_veto_segments(connection, veto_segments_name)
	if program_name == "rinca": veto_segments = ligolw_segments.segmenttable_get_by_name(xmldoc, veto_segments_name).coalesce()
	return veto_segments


def get_segments(connection, xmldoc, program_name):
	seglists = segments.segmentlistdict()
	if program_name == "thinca":
		seglists = db_thinca_rings.get_thinca_zero_lag_segments(connection, program_name)
	if program_name == "gstlal_inspiral" or program_name == "lalapps_ring":
		seglists = ligolw_search_summary.segmentlistdict_fromsearchsummary(xmldoc, program_name).coalesce()
	return seglists

def get_background_livetime_by_slide(connection, program_name, seglists, veto_segments=None, verbose = False):

	if program_name == "thinca":
		return background_livetime_ring_by_slide(connection, program_name, seglists, veto_segments, verbose)

	if program_name == "gstlal_inspiral":
		return background_livetime_nonring_by_slide(connection, seglists, veto_segments, verbose)

def get_livetime_by_offset(background_by_slide, zero_live_time=None, det='L1'):
	lt = {}
	if zero_live_time:
		for inst in zero_live_time.keys():
			lt.setdefault(inst,{})
			lt[inst][0] = zero_live_time[inst]

	for inst, slide in background_by_slide.items():
		lt.setdefault(inst,{})
		for offsets, value in slide.items():
			for ifo, offset in offsets:
				if ifo == det:
					lt[inst][offset] = value
	return lt
