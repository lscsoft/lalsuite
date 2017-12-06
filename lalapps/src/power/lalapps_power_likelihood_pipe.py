#
# Copyright (C) 2007  Kipp C. Cannon
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
Excess power offline pipeline's likelihood stage construction script.
"""


import ConfigParser
import glob
from optparse import OptionParser
import sys
import tempfile


from glue import iterutils
from glue import pipeline
from glue import segmentsUtils
from glue.lal import CacheEntry
from pylal.date import LIGOTimeGPS
from lalapps import power


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
		usage = "%prog [options]",
		description = "Constructs the likelihood-ratio based coincidence stage for an excess power analysis.  The input consists of one or more LAL caches listing the sqlite database trigger files, and a list of segments giving the time intervals that should be considered to be independent.  The LAL caches list all trigger files together, that is injections, time slides, and zero-lag.  The individual trigger files are self-describing, so the analysis codes can autodetect their type.  Each segment will be analyzed using the files that intersect it:  the likelihood ratios will be constructed from the injections and time-lag triggers contained in files that intersect the segment, and that data used to assign likelihoods to the injections, time-lag, and zero-lag coincs in all files that intersect the same segment."
	)
	parser.add_option("--input-cache", metavar = "filename", action = "append", default = [], help = "Add the contents of this cache file to the list of files from which to draw statistics.")
	parser.add_option("--round-robin-cache", metavar = "filename", action = "append", default = [], help = "Add the contents of this cache file to the list of files from which to draw injection statistics in a round-robin way.")
	parser.add_option("--condor-log-dir", metavar = "path", default = ".", help = "Set the directory for Condor log files (default = \".\").")
	parser.add_option("--config-file", metavar = "filename", default = "power.ini", help = "Set .ini configuration file name (default = \"power.ini\").")
	parser.add_option("--distribution-segments", metavar = "filename", help = "Read boundaries for distribution data intervals from this segwizard format segments file (required).")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()

	if options.distribution_segments is None:
		raise ValueError, "missing required argument --distribution-segments"
	options.distribution_segments = segmentsUtils.fromsegwizard(file(options.distribution_segments), coltype = LIGOTimeGPS)

	options.input_cache = set([CacheEntry(line) for filename in options.input_cache for line in file(filename)])
	options.round_robin_cache = [set(map(CacheEntry, file(filename))) for filename in options.round_robin_cache]

	return options, (filenames or [])


#
# =============================================================================
#
#                                    Config
#
# =============================================================================
#


def parse_config_file(options):
	if options.verbose:
		print >>sys.stderr, "reading %s ..." % options.config_file
	config = ConfigParser.SafeConfigParser()
	config.read(options.config_file)

	options.tag = config.get("pipeline", "user_tag")
	options.likelihood_data_cache_base = config.get("pipeline", "likelihood_data_cache_base")

	return config


#
# =============================================================================
#
#                                 Place Holder
#
# =============================================================================
#


class PlaceHolder(object):
	def __init__(self):
		self.cache = set()

	def add_input_cache(self, cache):
		self.cache |= cache

	def get_output_cache(self):
		return self.cache


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
# Parse .ini file, input cache(s), and segment list.
#


config_parser = parse_config_file(options)


#
# Define .sub files
#


power.init_job_types(config_parser)


#
# Start DAG
#


power.make_dag_directories(config_parser)
dag = pipeline.CondorDAG(tempfile.mkstemp(".log", "power_likelihood_", options.condor_log_dir)[1])
dag.set_dag_file("power_likelihood")


#
# Generate likelihood data
#


input_cache_nodes = set()
round_robin_cache_nodes = [set() for cache in options.round_robin_cache]
for seg in options.distribution_segments:
	if options.verbose:
		print >>sys.stderr, "generating distribution measurement jobs for %s ..." % str(seg)
	input_cache_nodes |= power.make_burca_tailor_fragment(dag, set([entry for entry in options.input_cache if entry.segmentlistdict.intersects_segment(seg)]), seg, "LIKELIHOOD_MAIN")
	for i, (nodes, cache) in enumerate(zip(round_robin_cache_nodes, options.round_robin_cache)):
		nodes |= power.make_burca_tailor_fragment(dag, set([entry for entry in cache if entry.segmentlistdict.intersects_segment(seg)]), seg, "LIKELIHOOD_RR%02d" % i)


#
# Compute likelihood ratios for coincs
#


if options.verbose:
	print >>sys.stderr, "generating likelihood assignment jobs for main group ..."
parents = reduce(lambda a, b: a | b, round_robin_cache_nodes, input_cache_nodes)
nodes = power.make_burca2_fragment(dag, options.input_cache, parents, "LIKELIHOOD_MAIN")


def round_robin(round_robin_cache_nodes, round_robin_cache):
	parents = list(iterutils.choices(round_robin_cache_nodes, len(round_robin_cache_nodes) - 1))
	parents.reverse()
	parents = [reduce(lambda a, b: a | b, seq) for seq in parents]
	return zip(parents, [cache for (cache,) in iterutils.choices(round_robin_cache, 1)])

for i, (parents, apply_to_cache) in enumerate(round_robin(round_robin_cache_nodes, options.round_robin_cache)):
	if options.verbose:
		print >>sys.stderr, "generating likelihood assignment jobs for round-robin group %d ..." % i
	nodes |= power.make_burca2_fragment(dag, apply_to_cache, parents | input_cache_nodes, "LIKELIHOOD_RR%02d" % i)


#
# Output
#


if options.verbose:
	print >>sys.stderr, "writing dag ..."
dag.write_sub_files()
dag.write_dag()
