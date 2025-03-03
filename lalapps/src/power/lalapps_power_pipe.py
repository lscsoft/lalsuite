##python
#
# Copyright (C) 2005--2010,2013,2015,2016,2019  Kipp Cannon
# Copyright (C) 2004--2006  Saikat Ray-Majumder
# Copyright (C) 2003--2005  Duncan Brown
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
Excess power offline pipeline construction script.
"""


import math
from optparse import OptionParser
import os
import sys
import tempfile
from configparser import (ConfigParser, NoOptionError)


import igwn_segments as segments
from igwn_segments import utils as segmentsUtils
# pipeline.py is required, but the version available in lalsuite is
# currently broken and we can't add a dependency on gstlal because it would
# create a cycle, so we make the import optional.  the program cannot run
# to completion without it, but it can print a --help message and exit with
# a success code, which is enough to make lalsuite's build scripts happy.
try:
	from gstlal import pipeline
except ImportError:
	pass
from lal import LIGOTimeGPS
from lal.utils import CacheEntry
from lalburst import cafe
from lalburst import timeslides
from lalburst import power


__author__ = "Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"
__version__ = "$Revision$"


#
# =============================================================================
#
#                                 Command Line
#
# =============================================================================
#


def parse_command_line():
	parser = OptionParser(
		version = "%prog CVS $Id$",
		description = "%prog builds an excess power pipeline DAG suitable for running at the various LSC Data Grid sites.  The script requires a configuration file.  An example file can be found in the LALApps CVS."
	)
	parser.add_option("--condor-log-dir", metavar = "path", default = ".", help = "Set the directory for Condor log files (default = \".\").")
	parser.add_option("--config-file", metavar = "filename", default = "power.ini", help = "Set .ini configuration file name (default = \"power.ini\").")
	parser.add_option("--full-segments", action = "store_true", help = "Analyze all data from segment lists, not just coincident times.")
	parser.add_option("--minimum-gap", metavar = "seconds", type = "float", default = 60.0, help = "Merge jobs analyzing data from the same instrument if the gap between them is less than this many seconds (default = 60).")
	parser.add_option("--variant", metavar = "[injections|noninjections|both]", default = "both", help = "Select the variant of the pipeline to construct.  \"injections\" produces a simulations-only version of the pipeline, \"noninjections\" produces a version with no simulation jobs, and \"both\" produces a full pipeline with both simulation and non-simulation jobs.")
	parser.add_option("--background-time-slides", metavar = "filename", default = [], action = "append", help = "Set file from which to obtain the time slide table for use in the background branch of the pipeline (default = \"background_time_slides.xml.gz\").  Provide this argument multiple times to provide multiple time slide files, each will result in a separate set of lalburst_coinc jobs.")
	parser.add_option("--injection-time-slides", metavar = "filename", help = "Set file from which to obtain the time slide table for use in the injection branch of the pipeline (default = \"injection_time_slides.xml.gz\").")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()

	if options.variant not in ("injections", "noninjections", "both"):
		raise ValueError("unrecognized --variant %s" % options.variant)
	options.do_injections = options.variant in ("injections", "both")
	options.do_noninjections = options.variant in ("noninjections", "both")

	if options.do_injections and not options.injection_time_slides:
		raise ValueError("missing required --injection-time-slides argument")
	if options.do_noninjections and not options.background_time_slides:
		raise ValueError("missing required --background-time-slides argument")

	# simplifies life later by allowing the background and injection
	# branches of the dag to be constructed with nearly identical code
	options.injection_time_slides = [options.injection_time_slides]

	return options, (filenames or ["power.dag"])


#
# =============================================================================
#
#                                    Config
#
# =============================================================================
#


def parse_config_file(options):
	if options.verbose:
		print("reading %s ..." % options.config_file, file=sys.stderr)
	config = ConfigParser()
	config.read(options.config_file)

	options.tag = config.get("pipeline", "user_tag")
	options.enable_clustering = config.getboolean("pipeline", "enable_clustering")

	seglistdict = segments.segmentlistdict()
	tiling_phase = {}
	for ifo in config.get("pipeline", "ifos").split():
		seglistdict[ifo] = segmentsUtils.fromsegwizard(open(config.get("pipeline", "seglist_%s" % ifo)), coltype = LIGOTimeGPS).coalesce()
		try:
			offset = config.getfloat("pipeline", "tiling_phase_%s" % ifo)
		except NoOptionError:
			offset = 0.0
		if offset:
			tiling_phase[ifo] = offset

	options.psds_per_power = config.getint("pipeline", "psds_per_power")
	options.psds_per_injection = config.getint("pipeline", "psds_per_injection")
	options.timing_params = power.TimingParameters(config)

	return seglistdict, tiling_phase, config


#
# =============================================================================
#
#                            Determine Segment List
#
# =============================================================================
#


def compute_segment_lists(seglistdict, time_slides, minimum_gap, timing_params, full_segments = True, verbose = False):
	if verbose:
		print("constructing segment list ...", file=sys.stderr)

	seglistdict = seglistdict.copy()

	if not full_segments:
		# cull too-short single-instrument segments from the input
		# segmentlist dictionary;  this can significantly increase
		# the speed of the get_coincident_segmentlistdict()
		# function when the input segmentlists have had many data
		# quality holes poked out of them
		power.remove_too_short_segments(seglistdict, timing_params)

		# extract the segments that are coincident under the time
		# slides
		new = cafe.get_coincident_segmentlistdict(seglistdict, time_slides)

		# adjust surviving segment lengths up to the next integer
		# number of PSDs
		for seglist in new.values():
			# Try Adjusting Upper Bounds:

			# count the number of PSDs in each segment
			psds = [power.psds_from_job_length(timing_params, float(abs(seg))) for seg in seglist]

			# round up to the nearest integer.
			psds = [int(math.ceil(max(n, 1.0))) for n in psds]

			# compute the duration of each job
			durations = [power.job_length_from_psds(timing_params, n) for n in psds]

			# update segment list
			for i, seg in enumerate(seglist):
				seglist[i] = segments.segment(seg[0], seg[0] + durations[i])

			# and take intersection with original segments to
			# not exceed original bounds
			new &= seglistdict

			# Try Adjusting Lower Bounds:

			# count the number of PSDs in each segment
			psds = [power.psds_from_job_length(timing_params, float(abs(seg))) for seg in seglist]

			# round up to the nearest integer.
			psds = [int(math.ceil(max(n, 1.0))) for n in psds]

			# compute the duration of each job
			durations = [power.job_length_from_psds(timing_params, n) for n in psds]

			# update segment list
			for i, seg in enumerate(seglist):
				seglist[i] = segments.segment(seg[1] - durations[i], seg[1])

			# and take intersection with original segments to
			# not exceed original bounds
			new &= seglistdict


		# try to fill gaps between jobs
		new.protract(minimum_gap / 2).contract(minimum_gap / 2)

		# and take intersection with original segments to not
		# exceed original bounds
		seglistdict &= new

	# remove segments that are too short
	power.remove_too_short_segments(seglistdict, timing_params)

	# done
	return seglistdict


#
# =============================================================================
#
#                               DAG Construction
#
# =============================================================================
#


#
# Command line
#


options, filenames = parse_command_line()


#
# Parse .ini file, loading the single-instrument segment lists while at it.
#


seglistdict, tiling_phase, config_parser = parse_config_file(options)


#
# Define .sub files
#


power.init_job_types(config_parser)


#
# Using time slide information, construct segment lists describing times
# requiring trigger construction.
#


if options.verbose:
	print("Computing segments for which lalapps_power jobs are required ...", file=sys.stderr)

background_time_slides = {}
background_seglistdict = segments.segmentlistdict()
if options.do_noninjections:
	for filename in options.background_time_slides:
		cache_entry = CacheEntry(None, None, None, "file://localhost" + os.path.abspath(filename))
		background_time_slides[cache_entry] = timeslides.load_time_slides(filename, verbose = options.verbose).values()
		background_seglistdict |= compute_segment_lists(seglistdict, background_time_slides[cache_entry], options.minimum_gap, options.timing_params, full_segments = options.full_segments, verbose = options.verbose)


injection_time_slides = {}
injection_seglistdict = segments.segmentlistdict()
if options.do_injections:
	for filename in options.injection_time_slides:
		cache_entry = CacheEntry(None, None, None, "file://localhost" + os.path.abspath(filename))
		injection_time_slides[cache_entry] = timeslides.load_time_slides(filename, verbose = options.verbose).values()
		injection_seglistdict |= compute_segment_lists(seglistdict, injection_time_slides[cache_entry], options.minimum_gap, options.timing_params, full_segments = options.full_segments, verbose = options.verbose)


# apply time shifts to segment lists to shift tiling phases, but take
# intersection with original segments to stay within allowed times.  Note:
# can't use segmentlistdict's offset mechanism to do this because we need
# the offsets to still be 0 for coincidence testing later.


for key, offset in tiling_phase.items():
	if key in background_seglistdict:
		background_seglistdict[key].shift(offset)
	if key in injection_seglistdict:
		injection_seglistdict[key].shift(offset)
background_seglistdict &= seglistdict
injection_seglistdict &= seglistdict


#
# Start DAG
#


power.make_dag_directories(config_parser)
dag = pipeline.CondorDAG(tempfile.mkstemp(".log", "power_", options.condor_log_dir)[1])
dag.set_dag_file(os.path.splitext(filenames[0])[0])


#
# Build datafind jobs.
#


datafinds = power.make_datafind_stage(dag, injection_seglistdict | background_seglistdict, verbose = options.verbose)


#
# Main analysis
#


def make_coinc_branch(dag, datafinds, seglistdict, time_slides, timing_params, psds_per_power, enable_clustering, tag, do_injections = False, verbose = False):
	# injection list


	if do_injections:
		assert len(time_slides) == 1
		if verbose:
			print("Building lalapps_binj jobs ...", file=sys.stderr)
		binjnodes = power.make_binj_fragment(dag, seglistdict.extent_all(), time_slides.keys()[0], tag, 0.0, float(power.powerjob.get_opts()["low-freq-cutoff"]), float(power.powerjob.get_opts()["low-freq-cutoff"]) + float(power.powerjob.get_opts()["bandwidth"]))
		# add binj nodes as parents of the datafinds to force the binj's to
		# be run first.  this ensures that once a datafind has run the
		# power jobs that follow it will immediately be able to run, which
		# helps depth-first dagman do smarter things.
		for node in datafinds:
			for binjnode in binjnodes:
				node.add_parent(binjnode)
	else:
		binjnodes = set()


	# single-instrument trigger generation


	trigger_nodes = power.make_single_instrument_stage(dag, datafinds, seglistdict, tag, timing_params, psds_per_power, binjnodes = binjnodes, verbose = verbose)
	if enable_clustering:
		if verbose:
			print("building pre-lladd bucluster jobs ...", file=sys.stderr)
		trigger_nodes = power.make_bucluster_fragment(dag, trigger_nodes, "PRELLADD_%s" % tag, verbose = verbose)


	# coincidence analysis


	coinc_nodes = set()
	binj_cache = set([cache_entry for node in binjnodes for cache_entry in node.get_output_cache()])
	# otherwise too many copies of the offset vector will be fed into
	# burca
	assert len(binj_cache) < 2
	for n, (time_slides_cache_entry, these_time_slides) in enumerate(time_slides.items()):
		if verbose:
			print("%s %d/%d (%s):" % (tag, n + 1, len(time_slides), time_slides_cache_entry.path), file=sys.stderr)
		tisi_cache = set([time_slides_cache_entry])
		if do_injections:
			# lalapps_binj has already copied the time slide
			# document into its own output
			extra_input_cache = set()
		else:
			# ligolw_add needs to copy the time slide document
			# into is output
			extra_input_cache = tisi_cache
		nodes = set()
		for seg, parents, cache, clipseg in power.group_coinc_parents(trigger_nodes, these_time_slides, verbose = verbose):
			nodes |= power.make_lladd_fragment(dag, parents | binjnodes, "%s_%d" % (tag, n), segment = seg, input_cache = cache | binj_cache, extra_input_cache = extra_input_cache, remove_input = do_injections, preserve_cache = binj_cache | tisi_cache)
		if enable_clustering:
			if verbose:
				print("building post-lladd bucluster jobs ...", file=sys.stderr)
			nodes = power.make_bucluster_fragment(dag, nodes, "POSTLLADD_%s_%d" % (tag, n), verbose = verbose)
		if verbose:
			print("building burca jobs ...", file=sys.stderr)
		coinc_nodes |= power.make_burca_fragment(dag, nodes, "%s_%d" % (tag, n), verbose = verbose)
		if verbose:
			print("done %s %d/%d" % (tag, n + 1, len(time_slides)), file=sys.stderr)


	# injection identification


	if do_injections:
		if verbose:
			print("building binjfind jobs ...", file=sys.stderr)
		coinc_nodes = power.make_binjfind_fragment(dag, coinc_nodes, tag, verbose = verbose)


	# conversion to SQLite database files


	if verbose:
		print("building sqlite jobs ...", file=sys.stderr)
	coinc_nodes = power.make_sqlite_fragment(dag, coinc_nodes, tag, verbose = verbose)


	# done


	power.write_output_cache(coinc_nodes, "%s_%s_output.cache" % (os.path.splitext(dag.get_dag_file())[0], tag))
	return coinc_nodes


coinc_nodes = make_coinc_branch(dag, datafinds, background_seglistdict, background_time_slides, options.timing_params, options.psds_per_power, options.enable_clustering, options.tag, do_injections = False, verbose = options.verbose)
inj_coinc_nodes = make_coinc_branch(dag, datafinds, injection_seglistdict, injection_time_slides, options.timing_params, options.psds_per_injection, options.enable_clustering, "INJECTIONS_RUN_0_%s" % options.tag, do_injections = True, verbose = options.verbose)


#
# Output
#


if options.verbose:
	print("writing dag ...", file=sys.stderr)
dag.write_sub_files()
dag.write_dag()
