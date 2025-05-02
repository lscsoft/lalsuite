##python
#
# Copyright (C) 2006,2013  Kipp Cannon
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


from optparse import OptionParser
import sys


from igwn_ligolw import ligolw
from igwn_ligolw import lsctables
from igwn_ligolw import utils as ligolw_utils
from igwn_ligolw.utils import segments as ligolw_segments
from igwn_ligolw.utils import process as ligolw_process
from igwn_ligolw.utils import search_summary as ligolw_search_summary
from lalburst import git_version
import igwn_segments as segments


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


process_program_name = "lalburst_cut"


#
# =============================================================================
#
#                                 Command Line
#
# =============================================================================
#


def parse_command_line():
	parser = OptionParser(
		version = "Name: %%prog\n%s" % git_version.verbose_msg,
		usage = "%prog [options] [file ...]",
		description = "Removes sngl_burst events from XML files according to a variety of criteria.  Files named on the command line are read one-by-one, and over-written with the new files.  If no files are named on the command line, input is read from stdin and written to stdout.  Note that, for the most part, this program does not understand coincidence information, and so if an injection or burst event is removed that participates in a coincidence, this program simply deletes the entire coincidence as well (before applying the --coinc-only cut)."
	)
	parser.add_option("--coinc-only", action = "store_true", help = "Discard events that are not participating in a coincident event.")
	parser.add_option("--comment", metavar = "text", help = "Set the comment string to be recorded in the process table for this job (default = None).")
	parser.add_option("--inj-made-only", action = "store_true", help = "Discard injections outside the search summary out segments.")
	parser.add_option("--min-amplitude", metavar = "value", type = "float", help = "Discard events below the given amplitude.")
	parser.add_option("--max-amplitude", metavar = "value", type = "float", help = "Discard events above the given amplitude.")
	parser.add_option("--min-bandwidth", metavar = "Hz", type = "float", help = "Discard events narrower than the given bandwidth.")
	parser.add_option("--max-bandwidth", metavar = "Hz", type = "float", help = "Discard events wider than the given bandwidth.")
	parser.add_option("--min-central-freq", metavar = "Hz", type = "float", help = "Discard events with central frequency lower than that given.")
	parser.add_option("--max-central-freq", metavar = "Hz", type = "float", help = "Discard events with central frequency higher than that given.")
	parser.add_option("--min-confidence", metavar = "value", type = "float", help = "Discard events below the given confidence.")
	parser.add_option("--max-confidence", metavar = "value", type = "float", help = "Discard events above the given confidence.")
	parser.add_option("--min-duration", metavar = "seconds", type = "float", help = "Discard events shorter than the given duration.")
	parser.add_option("--max-duration", metavar = "seconds", type = "float", help = "Discard events longer than the given duration.")
	parser.add_option("--min-fhigh", metavar = "Hz", type = "float", help = "Discard events with highest frequency below the given frequency.")
	parser.add_option("--max-fhigh", metavar = "Hz", type = "float", help = "Discard events with highest frequency above the given frequency.")
	parser.add_option("--min-flow", metavar = "Hz", type = "float", help = "Discard events with lowest frequency below the given frequency.")
	parser.add_option("--max-flow", metavar = "Hz", type = "float", help = "Discard events with loest frequency above the given frequency.")
	parser.add_option("--min-hrss", metavar = "value", type = "float", help = "Discard events with h_rss below the given value.")
	parser.add_option("--max-hrss", metavar = "value", type = "float", help = "Discard events with h_rss above the given value.")
	parser.add_option("--cut-instrument", metavar = "name", action = "append", default = [], help = "Discard events from given instrument.")
	parser.add_option("--min-peak-time", metavar = "seconds", help = "Discard events with peak time before the given GPS time.")
	parser.add_option("--max-peak-time", metavar = "seconds", help = "Discard events with peak time after the given GPS time.")
	parser.add_option("--min-start-time", metavar = "seconds", help = "Discard events starting before the given GPS time.")
	parser.add_option("--max-start-time", metavar = "seconds", help = "Discard events starting after the given GPS time.")
	parser.add_option("--min-stop-time", metavar = "seconds", help = "Discard events ending before the given GPS time.")
	parser.add_option("--max-stop-time", metavar = "seconds", help = "Discard events ending after the given GPS time.")
	parser.add_option("--min-snr", metavar = "value", type = "float", help = "Discard events below the given SNR.")
	parser.add_option("--max-snr", metavar = "value", type = "float", help = "Discard events above the given SNR.")
	parser.add_option("--program", metavar = "name", help = "Process events generated by the given program.")
	parser.add_option("--veto-file", metavar = "filename", help = "Veto events using the veto segment list extracted from this XML file.  The file must contain segment tables, and the veto list will be constructed from the segments named \"sngl_burst_veto\".")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()

	paramdict = options.__dict__.copy()

	if options.inj_made_only and not options.program:
		raise ValueError("must set --program when --inj-made-only is set")
	options.cut_instrument = set(options.cut_instrument)

	if options.min_peak_time is not None:
		options.min_peak_time = lsctables.LIGOTimeGPS(options.min_peak_time)
	if options.max_peak_time is not None:
		options.max_peak_time = lsctables.LIGOTimeGPS(options.max_peak_time)
	if options.min_start_time is not None:
		options.min_start_time = lsctables.LIGOTimeGPS(options.min_start_time)
	if options.max_start_time is not None:
		options.max_start_time = lsctables.LIGOTimeGPS(options.max_start_time)
	if options.min_stop_time is not None:
		options.min_stop_time = lsctables.LIGOTimeGPS(options.min_stop_time)
	if options.max_stop_time is not None:
		options.max_stop_time = lsctables.LIGOTimeGPS(options.max_stop_time)

	return options, paramdict, (filenames or [None])


#
# =============================================================================
#
#                                    Input
#
# =============================================================================
#


#
# Content handler
#


#
# =============================================================================
#
#                                 Preparation
#
# =============================================================================
#


def load_veto_segments(filename, verbose = False, contenthandler = None):
	return ligolw_segments.segmenttable_get_by_name(ligolw_utils.load_filename(filename, verbose = verbose, contenthandler = contenthandler), "sngl_burst_veto").coalesce()


class DocContents(object):
	def __init__(self, xmldoc, program = None):
		#
		# Find the out segments
		#

		self.outsegs = ligolw_search_summary.segmentlistdict_fromsearchsummary(xmldoc, program).coalesce()

		#
		# Find the sngl_burst table
		#

		self.snglbursttable = lsctables.SnglBurstTable.get_table(xmldoc)

		#
		# Get the list of process IDs we care about
		#

		self.process_ids = set(self.snglbursttable.getColumnByName("process_id"))
		if program is not None:
			self.process_ids &= lsctables.ProcessTable.get_table(xmldoc).get_ids_by_program(program)

		#
		# Find the sim_burst table, or make a fake one
		#

		try:
			self.simbursttable = lsctables.SimBurstTable.get_table(xmldoc)
		except:
			self.simbursttable = []

		#
		# Find the coinc tables, or make fake ones
		#

		try:
			self.coinctable = lsctables.CoincTable.get_table(xmldoc)
			self.coincmaptable = lsctables.CoincMapTable.get_table(xmldoc)
			self.multibursttable = lsctables.MultiBurstTable.get_table(xmldoc)
		except:
			self.coinctable = []
			self.coincmaptable = []
			self.multibursttable = []


#
# =============================================================================
#
#                                     Cuts
#
# =============================================================================
#


def remove_events_by_segment(contents, veto_segments):
	ids = set()
	for i in xrange(len(contents.snglbursttable) - 1, -1, -1):
		burst = contents.snglbursttable[i]
		if burst.process_id in contents.process_ids and burst.ifo in veto_segments and veto_segments[burst.ifo].intersects_segment(burst.period):
			ids.add(burst.event_id)
			del contents.snglbursttable[i]
	return ids


def remove_events_by_parameters(contents, testfunc):
	ids = set()
	for i in xrange(len(contents.snglbursttable) - 1, -1, -1):
		burst = contents.snglbursttable[i]
		if burst.process_id in contents.process_ids and not testfunc(burst):
			ids.add(burst.event_id)
			del contents.snglbursttable[i]
	return ids


def remove_non_coincidences(contents):
	coinc_burst_ids = set(row.event_id for row in contents.coincmaptable.getColumnByName("event_id") if row.table_name == "sngl_burst")
	for i in xrange(len(contents.snglbursttable) - 1, -1, -1):
		burst = contents.snglbursttable[i]
		if burst.process_id in contents.process_ids and burst.event_id not in coinc_burst_ids:
			del contents.snglbursttable[i]


def remove_skipped_injections(contents):
	ids = set()
	for i in xrange(len(contents.simbursttable) - 1, -1, -1):
		sim = contents.simbursttable[i]
		if True not in (sim.time_at_instrument(instrument) in seglist for instrument, seglist in contents.outsegs.iteritems()):
			ids.add(sim.simulation_id)
			del contents.simbursttable[i]
	return ids


def clean_coinc_tables(contents, removed_ids):
	# FIXME FIXME FIXME:  this is broken since the conversion to
	# integer IDs.  need to include table name constraints in these
	# tests.

	# remove dangling coinc_event_map rows
	removed_coinc_ids = set()
	for i in xrange(len(contents.coincmaptable) - 1, -1, -1):
		if contents.coincmaptable[i].event_id in removed_ids:
			removed_coinc_ids.add(contents.coincmaptable[i].coinc_event_id)
			del contents.coincmaptable[i]

	# remove broken coinc_event rows
	for i in xrange(len(contents.coinctable) - 1, -1, -1):
		if contents.coinctable[i].coinc_event_id in removed_coinc_ids:
			del contents.coinctable[i]

	# remove dangling coinc_event_map rows
	for i in xrange(len(contents.coincmaptable) - 1, -1, -1):
		if contents.coincmaptable[i].coinc_event_id in removed_coinc_ids:
			del contents.coincmaptable[i]

	# remove dangling multi_burst rows
	for i in xrange(len(contents.multibursttable) - 1, -1, -1):
		if contents.multibursttable[i].coinc_event_id in removed_coinc_ids:
			del contents.multibursttable[i]

	# recurse (e.g., removes injection coincs that point to the
	# coinc_events that got deleted)
	if removed_coinc_ids:
		clean_coinc_tables(contents, removed_coinc_ids)


#
# =============================================================================
#
#                                 Library API
#
# =============================================================================
#


def apply_filters(contents, burst_test_func, veto_segments, del_non_coincs = False, del_skipped_injections = False, verbose = False):
	removed_ids = set()
	if veto_segments:
		if verbose:
			print("applying veto segment list ...", file=sys.stderr)
		removed_ids |= remove_events_by_segment(contents, veto_segments)
	if verbose:
		print("filtering sngl_burst rows by parameters ...", file=sys.stderr)
	removed_ids |= remove_events_by_parameters(contents, burst_test_func)
	if del_skipped_injections:
		if verbose:
			print("removing injections that weren't performed ...", file=sys.stderr)
		remove_skipped_injections(contents)
	if verbose:
		print("removing broken coincidences ...", file=sys.stderr)
	clean_coinc_tables(contents, removed_ids)
	if del_non_coincs:
		if verbose:
			print("removing non-coincident events ...", file=sys.stderr)
		remove_non_coincidences(contents)


def ligolw_bucut(xmldoc, burst_test_func, veto_segments = segments.segmentlistdict(), del_non_coincs = False, del_skipped_injections = False, program = None, comment = None, verbose = False):
	process = ligolw_process.register_to_xmldoc(xmldoc, process_program_name, paramdict, version = __version__, cvs_repository = "lscsoft", cvs_entry_time = __date__, comment = comment)

	contents = DocContents(xmldoc, program)

	apply_filters(contents, burst_test_func, veto_segments, del_non_coincs = del_non_coincs, del_skipped_injections = del_skipped_injections, verbose = verbose)

	seg = contents.outsegs.extent_all()
	ligolw_search_summary.append_search_summary(xmldoc, process, inseg = seg, outseg = seg, nevents = len(contents.snglbursttable))

	process.set_end_time_now()

	return xmldoc


#
# =============================================================================
#
#                                     Main
#
# =============================================================================
#


#
# Parse command line.
#


options, paramdict, filenames = parse_command_line()


#
# Define sngl_burst test function.
#


def make_keep_this_sngl_burst(options):
	def add_test(func, t):
		return lambda arg: t(arg) and func(arg)

	keep_this_sngl_burst = lambda burst: True

	if options.min_amplitude is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.amplitude >= options.min_amplitude)
	if options.max_amplitude is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.amplitude <= options.max_amplitude)
	if options.min_bandwidth is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.bandwidth >= options.min_bandwidth)
	if options.max_bandwidth is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.bandwidth <= options.max_bandwidth)
	if options.min_central_freq is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.central_freq >= options.min_central_freq)
	if options.max_central_freq is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.central_freq <= options.max_central_freq)
	if options.min_confidence is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.confidence >= options.min_confidence)
	if options.max_confidence is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.confidence <= options.max_confidence)
	if options.min_duration is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.duration >= options.min_duration)
	if options.max_duration is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.duration <= options.max_duration)
	if options.min_fhigh is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.fhigh >= options.min_fhigh)
	if options.max_fhigh is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.fhigh <= options.max_fhigh)
	if options.min_flow is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.flow >= options.min_flow)
	if options.max_flow is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.flow <= options.max_flow)
	if options.min_hrss is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.hrss >= options.min_hrss)
	if options.max_hrss is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.hrss <= options.max_hrss)
	if options.min_snr is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.snr >= options.min_snr)
	if options.max_snr is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.snr <= options.max_snr)
	if options.cut_instrument:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.ifo not in options.cut_instrument)
	if options.min_peak_time is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.peak >= options.min_peak_time)
	if options.max_peak_time is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.peak <= options.max_peak_time)
	if options.min_start_time is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.start >= options.min_start_time)
	if options.max_start_time is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.start <= options.max_start_time)
	if options.min_stop_time is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.stop >= options.min_stop_time)
	if options.max_stop_time is not None:
		keep_this_sngl_burst = add_test(keep_this_sngl_burst, lambda burst: burst.stop <= options.max_stop_time)

	return keep_this_sngl_burst


keep_this_sngl_burst = make_keep_this_sngl_burst(options)



#
# Get veto segment information.
#


if options.veto_file:
	veto_segments = load_veto_segments(options.veto_file, verbose = options.verbose)
else:
	veto_segments = segments.segmentlistdict()


#
# Do the work.
#


for filename in filenames:
	xmldoc = ligolw_utils.load_filename(filename, verbose = options.verbose)
	xmldoc = ligolw_bucut(xmldoc, keep_this_sngl_burst, veto_segments = veto_segments, del_non_coincs = options.coinc_only, del_skipped_injections = options.inj_made_only, program = options.program, comment = options.comment, verbose = options.verbose)
	ligolw_utils.write_filename(xmldoc, filename, verbose = options.verbose)
	xmldoc.unlink()
