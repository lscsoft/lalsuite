# Copyright (C) 2013 Chris Pankow
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation,.

## \defgroup laldetchar_py_hveto hveto utility modules
## \ingroup laldetchar_python
"""A collection of tools for using, plotting, and interpreting HVeto results
as produced by the laldetchar-hveto program.

glib_utils: A collection of functions to better handle the GLib SWIG wrapped parts of hveto.
plot_utils: Some functions to plot triggers and the like from HVeto rounds.
"""
# \author Chris Pankow (<chris.pankow@ligo.org>)
# ### Synopsis ###
# ~~~
# from laldetchar import hveto
# ~~~

import sys
import re
import json
import itertools
from collections import defaultdict

from glue.lal import Cache
from glue.segments import segment, segmentlist, segmentlistdict
from glue.segmentsUtils import fromsegwizard
from glue.ligolw import ligolw, utils, lsctables, ilwd, table
from glue.ligolw.utils import process
from glue.ligolw.utils import segments as ligolw_segments
from glue.segmentdb.query_engine import LdbdQueryEngine
from glue.segmentdb import segmentdb_utils

# SWIG bindings
import laldetchar
import lal
# FIXME: Will this be needed with Karl's bindings?
# Yes... lalburst doesn't import the wrapped version of SnglBurst, so trying to
# get its attributes will fail until they are defined with this
# But that's been fixed, so we should try it without
import lalmetaio

from laldetchar import git_version as version
__author__ = "Chris Pankow <chris.pankow@ligo.org>"
__version__ = version.id
__date__ = version.date

## \addtogroup laldetchar_py_hveto
#@{

#
# NOTE: This is mostly a re-implementation of some of the framework within
# ligolw_burca because its very explicit in its dislike of intra-instrument
# coincidences. Since we construct coincidences explicitly through the algorihtm
# this is mostly just to record them
#
try:
	from pylal.ligolw_burca import make_multi_burst
	from pylal import snglcoinc
except ImportError:
	sys.exit( "Couldn't import from pylal, please install first." )

#
# CoincDefs
#
HVetoBBCoincDef = lsctables.CoincDef(search = u"hveto", search_coinc_type = 0, description = u"sngl_burst<-->sngl_burst coincidences")

class HVetoCoincTables(snglcoinc.CoincTables):
	def __init__(self, xmldoc):
		snglcoinc.CoincTables.__init__(self, xmldoc)

		# find the multi_burst table or create one if not found
		try:
			self.multibursttable = lsctables.table.get_table(xmldoc, lsctables.MultiBurstTable.tableName)
		except ValueError:
			self.multibursttable = lsctables.New(lsctables.MultiBurstTable, ("process_id", "duration", "central_freq", "bandwidth", "snr", "confidence", "amplitude", "coinc_event_id"))
			xmldoc.childNodes[0].appendChild(self.multibursttable)

	def append_coinc(self, process_id, time_slide_id, coinc_def_id, events):
	        coinc = snglcoinc.CoincTables.append_coinc(self, process_id, time_slide_id, coinc_def_id, events)
		# FIXME: Reenable if we add ms colums to sb
	        #self.multibursttable.append(make_multi_burst(process_id, coinc.coinc_event_id, events, self.time_slide_index[time_slide_id]))
	        return coinc

class HVetoResult(object):
	ROUND_KEY = re.compile("\((\d+), (\d*\.\d+), (\d*\.\d+)\)")
	@staticmethod
	def make_key(n, snr_t, wind):
		return str((n, snr_t, wind))

	@staticmethod
	def decode_key(kstr):
		m = re.match(HVetoResult.ROUND_KEY, kstr)
		try:
			return (int(m.group(1)), float(m.group(2)), float(m.group(3)))
		except IndexError, AttributeError:
			raise "Could not parse round values from key string."

	def __init__(self):
		# Yep. That is in fact a dictionary of dictionaries that return 0 by default. Confused yet?
		self.counts = defaultdict(dict)
		self.coincs = defaultdict(dict)

	def to_json(self, fobj):
		json.dump({"count": self.counts, "coinc": self.coincs}, fobj)

	def from_json(self, fobj):
		tmpdata = json.load(fobj)
		self.coincs.update(tmpdata["coinc"])
		self.counts.update(tmpdata["count"])

	def add_round(self, n, snr_t, wind, chan, chancount, chancoinc):
		"""
		Add a coinc and count matrix (a dictionary, keyed by channel to integer values, with a key corresponding to the round number (n), SNR threshold (snr_t) and time window (wind).
		"""
		self.counts[HVetoResult.make_key(n, snr_t, wind)][chan] = chancount
		self.coincs[HVetoResult.make_key(n, snr_t, wind)][chan] = chancoinc
	def get_keys(self):
		return set(self.counts.keys() + self.coincs.keys())

	def get_chan_list(self, n, snr_t, wind):
		"""
		Returns the union of the channels present in the internal coincidence and count dictionaries. The result is sorted so to preserve order.
		FIXME: This will only get the channels that are defined in the rows, not columns.
		"""
		return sorted(set(self.counts[HVetoResult.make_key(n, snr_t, wind)].keys() + self.coincs[HVetoResult.make_key(n, snr_t, wind)].keys()))

	def get_round_info(self, n, snr_t, wind, chan=None, crosschan=None, transpose=False):
		"""
		Retrieve information from a specific correlation matrix set, keyed by round number (n), SNR threshold (snr_t), and time window (wind). If neither chan nor crosschan is specified, then the entire matrix will be returned. If only chan is specified, a row will be returned. If both chan and crosschan are given, then the entry for that matrix will be returned. Optionally, setting transpose to True will couple the count and coinc values, e.g. they will be returned as a set, rather than independently.
		NOTE: Uses self.get_chan_list internally to retrieve values, so deterministic order is guaranteed.
		"""
		subkey = HVetoResult.make_key(n, snr_t, wind)
		count, coinc = self.counts[subkey], self.coincs[subkey]
		if count is None or coinc is None:
			raise ValueError("Could not retrieve count or coinc for round specifier %d %f %f" % (n, snr_t, wind))
		if chan is not None:
			if count.has_key(chan):
				if crosschan is not None and count[chan].has_key(crosschan):
					count = count[chan][crosschan]
				elif crosschan is not None:
					# We don't have a value for this
					count = 0
				else:
					count = [count[chan][cchan] for cchan in sorted(count[chan].keys())]
			if coinc.has_key(chan):
				if crosschan is not None and coinc[chan].has_key(crosschan):
					coinc = coinc[chan][crosschan]
				elif crosschan is not None:
					# We don't have a value for this
					coinc = 0
				else:
					coinc = [coinc[chan][cchan] for cchan in sorted(coinc[chan].keys())]
		else:
			chanlist = self.get_chan_list(n, snr_t, wind)
			count = [count[c1][c2] for c2 in sorted(count[c1].keys()) for c1 in chanlist]
			coinc = [coinc[c1][c2] for c2 in sorted(coinc[c1].keys()) for c1 in chanlist]

		if transpose:
			# Matrix?
			result = []
			if type(coinc) is list and type(coinc[0]) is list and type(count) is list and type(count[0]) is list:
				for c1row, c2row in zip(count, coinc):
					result.append(zip(c1row, c2ro2))
			# Row?
			elif type(coinc) is list and type(coinc[0]) is not list and type(count) is list and type(count[0]) is not list:
				result = zip(count, coinc)
			# Values
			else:
				result = (count, coinc)
			return result

		return count, coinc


def write_coinc_tables( vetotrigs, xmldoc, refchannel, twind, time_slide_id=None):
	"""
	Write a set of coinc tables for this round. We only write coincidences for coincs with refchannel. Note: This is probably gonna be slow... aaaand that's why we implemented the real algorithm in C.
	"""
	# Retrieve process information
	process = [ p for p in table.get_table( xmldoc, lsctables.ProcessTable.tableName ) if p.program == "laldetchar-hveto" ][0]
	process_id = process.process_id

	# Insert a time slide ID. It's not yet really necessary
	if time_slide_id is None:
		timeslidetable = lsctables.New(lsctables.TimeSlideTable)
		time_slide = timeslidetable.RowType
		time_slide.process_id = process_id
		time_slide.time_slide_id = time_slide_id = ilwd.ilwdchar( "time_slide:time_slide_id:0" )
		time_slide.instrument = opt.instrument
		time_slide.offset = 0.0
		timeslidetable.append(time_slide)
		xmldoc.childNodes[0].appendChild( timeslidetable )

	# Set up coinc tables
	coinc_def = HVetoBBCoincDef
	coincdeftable = lsctables.New(lsctables.CoincDefTable)
	coinc_def.coinc_def_id = coinc_def_id = coincdeftable.get_next_id()
	coincdeftable.append( coinc_def )
	xmldoc.childNodes[0].appendChild( coincdeftable )

	coinc_def = HVetoCoincTables( xmldoc )
	reftrigs = [ (segment( sb.get_peak()-twind/2.0, sb.get_peak()+twind/2.0 ), sb) for sb in vetotrigs if sb.channel == refchannel ]
	for vt in vetotrigs:
		if vt.channel == refchannel:
			continue
		for (s, t) in reftrigs:
			if vt.get_peak() in s:
				coinc_def.append_coinc( process_id, time_slide_id, coinc_def_id, (t, vt))
	return xmldoc

#
# Segment querying utilities
#

def query_segments_xml( xml_location, gps_start, gps_end, spec ):
	"""
	Retrieve the segment table from a location, and clip segments to (gps_start, gps_end). If spec is given, retrieve only segments with this definer, otherwise, get all of them.
	"""
	if spec is None:
		spec = True
	else:
		#ifo, definer, version = spec.split(":")
		definer = spec.split(":")
		ifo, definer, version = definer[0], ":".join(definer[1:-1]), definer[-1]
	xmldoc = utils.load_filename( xml_location )
	segment_definer = table.get_table( xmldoc, lsctables.SegmentDefTable.tableName )
	# FIXME: ifo in ifos? What does a segment for a set of ifos even mean?
	seg_def_id = [ sd.segment_def_id for sd in segment_definer if spec and ifo in sd.get_ifos() and definer == sd.name ]
	if len(seg_def_id) != 1:
		raise ValueError( "Need exactly one definer row for %s:%s:%s, got %d" % (ifo, definer, version, len(seg_def_id)) )
	seg_def_id = seg_def_id[0]

	segment = table.get_table( xmldoc, lsctables.SegmentTable.tableName )
	return segmentlist([s.get() for s in segment if s.segment_def_id == seg_def_id])

def query_segments_db( db_location, gps_start, gps_end, spec ):
	"""
	Query db_location to get segments between (gps_start, gps_end), with definer spec.
	"""
	engine = LdbdQueryEngine( segmentdb_utils.setup_database(db_location) )
	ifo, definer, version = spec.split(":")
	definer_args = [[ifo, definer, int(version), gps_start, gps_end, 0, 0]]
	result = segmentdb_utils.query_segments( engine, "segment", definer_args )
	return segmentlist(result[0])

def append_summ_vars(xmldoc, procid, **summvars):
  """
  Append round information in the form of SummVars.
  """

  try:
    summtable = table.get_table(xmldoc, lsctables.SearchSummVarsTable.tableName)
    procid = summtable[0].process_id
  except ValueError:
    summtable = lsctables.New(lsctables.SearchSummVarsTable, lsctables.SearchSummVarsTable.validcolumns.keys())

  for name, value in summvars.iteritems():
    summvar = summtable.RowType()
    summvar.name = name
    if isinstance(value, str):
      summvar.string = str(value)
      summvar.value = -1.0
    else:
      summvar.string = str(value)
      summvar.value = float(value)
    summvar.process_id = procid
    summvar.search_summvar_id = summtable.get_next_id()
    summtable.append(summvar)

  xmldoc.childNodes[0].appendChild(summtable)

def write_round_xml( vetosegs, vetotrigs, winner, ifo, opts ):
	"""
	Write out the products from this round of hveto: veto segments and the vetoed triggers.
	"""
	# Create a new document
	xmldoc = ligolw.Document()
	xmldoc.appendChild( ligolw.LIGO_LW() )

	# Append the process information
	procrow = utils.process.append_process( xmldoc, program="laldetchar-hveto" )
	utils.process.append_process_params( xmldoc, procrow, utils.process.process_params_from_dict(opts) )

	# Add the vetoed triggers
	xmldoc.childNodes[0].childNodes.append( vetotrigs )

	# Append the veto segments
	segl = segmentlistdict()
	segl[ifo] = vetosegs
	lwsegs = ligolw_segments.LigolwSegments( xmldoc )
	lwsegs.insert_from_segmentlistdict( segl, "hveto" )
	lwsegs.finalize(procrow)

	return xmldoc

#
# SWIG utilities
#

def lalburst_sb_to_glue_sb( sb_in, desired_columns ):
	"""
	Convert a lalburst SnglBurst structure to a SnglBurst row from glue.
	"""
	sb = lsctables.SnglBurstTable.RowType()
	for att in desired_columns:
		if att == "start_time":
			start_time = float( sb_in.start_time )
			sb.start_time = int(start_time)
			sb.start_time_ns = int(1e9*(start_time - int(start_time) ))
			continue
		elif att == "peak_time":
			peak_time = float( sb_in.peak_time )
			sb.peak_time = int(peak_time)
			sb.peak_time_ns = int(1e9*(peak_time - int(peak_time) ))
			continue
		elif att == "process_id":
			sb.process_id = ilwd.ilwdchar( "process:process_id:%d" % getattr( sb_in, att ) )
			continue
		elif att == "event_id":
			sb.event_id = ilwd.ilwdchar( "sngl_burst:sngl_burst_id:%d" % getattr( sb_in, att ) )
			continue

		try:
			setattr( sb, att, getattr( sb_in, att ) )
		except AttributeError:
			pass
	return sb

def write_wiki_page():
	"""
	TODO: Write this.
	"""
	raise NotImplementedError

# FIXME: These are here because glib_utils requires a function from this module
from . import glib_utils
from . import plot_utils


# close doxygen
##
#	\defgroup laldetchar_py_hveto_glib_utils	Glib Utils
#	\defgroup laldetchar_py_hveto_plot_utils	Plotting Utils
#@}
