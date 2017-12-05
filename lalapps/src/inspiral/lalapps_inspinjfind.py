#
# Copyright (C) 2006--2009,2013,2014,2016,2017  Kipp Cannon
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
Command-line interface to CBC injection identification code.
"""


from optparse import OptionParser
import sys


from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils
from glue.ligolw.utils import process as ligolw_process
from lalinspiral import inspinjfind


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
from lalapps.git_version import date as __date__
from lalapps.git_version import version as __version__


process_program_name = "lalapps_inspinjfind"


#
# Must use injection finder's custom sngl_inspiral table row type.
#


lsctables.SnglInspiralTable.RowType = lsctables.SnglInspiral = inspinjfind.SnglInspiral


#
# =============================================================================
#
#                                 Command Line
#
# =============================================================================
#


def parse_command_line():
	parser = OptionParser(
		version = "Name: %%prog\n%s" % __version__,
		usage = "%prog [options] [file ...]",
		description = "Accepts as input one or more LIGO Light Weight XML files, each containing CBC candidates and a list of injections, and adds entries to the coincidence tables indicating which events match which injections."
	)
	parser.add_option("-f", "--force", action = "store_true", help = "Process even if file has already been processed.")
	parser.add_option("--comment", metavar = "text", help = "Set the comment string to be written to the process table (default = None).")
	parser.add_option("-c", "--match-algorithm", metavar = "[inspiral]", default = "inspiral", help = "Set the algorithm used to match event candidates with injections (default = \"inspiral\").")
	parser.add_option("-r", "--revert", action = "store_true", help = "Revert a previously-processed XML file to its pre-lalapps_inspinjfind state, allowing it to be re-processed.  This is achieved by removing lalapps_inspinjfind's entries from the process metadata tables, removing injection-related entries from the coinc_definer table, and removing all coincs created by lalapps_inspinjfind.  NOTE:  generally this is only possible if lalapps_inspinjfind was the last transformation applied to the file, otherwise data that has been added subsequently might be broken due to dangling cross-references to deleted data.  It is left as an exercise to the user to ensure that reverting the file is sensible.")
	parser.add_option("-w", "--time-window", metavar = "seconds", type = 'float', default = 9., help = "Set the time window used by the \"inspiral\" match algorithm (default = 9.0).")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()

	paramdict = options.__dict__.copy()

	if options.match_algorithm is None:
		raise ValueError("missing required --match-algorithm option")
	if options.match_algorithm not in ("inspiral",):
		raise ValueError("unrecognized --match-algorithm \"%s\"" % options.match_algorithm)

	if options.time_window < 0.:
		raise ValueError("--time-window must be non-negative")

	return options, paramdict, (filenames or [None])


#
# =============================================================================
#
#                                     Main
#
# =============================================================================
#


#
# command line
#


options, paramdict, filenames = parse_command_line()

# must match columns in sngl_inspiral table
search = {
	"inspiral": "inspiral"
}[options.match_algorithm]


#
# set compare functions
#


def InspiralSnglCompare(sim, sngl, twindow = lsctables.LIGOTimeGPS(options.time_window)):
	"""
	Returns 0 if sim and sngl are less than twindow apart.
	"""
	dt = sim.time_geocent - sngl.end
	return 0 if abs(dt) < twindow else dt


NearCoincCompare = InspiralSnglCompare


snglcomparefunc, nearcoinccomparefunc = {
	"inspiral": (InspiralSnglCompare, NearCoincCompare)
}[options.match_algorithm]


#
# loop over files
#


for n, filename in enumerate(filenames, start = 1):
	#
	# load the document
	#

	if options.verbose:
		print >>sys.stderr, "%d/%d:" % (n, len(filenames)),
	xmldoc = ligolw_utils.load_filename(filename, contenthandler = inspinjfind.LIGOLWContentHandler, verbose = options.verbose)

	#
	# process
	#

	if options.revert:
		inspinjfind.revert(xmldoc, process_program_name, verbose = options.verbose)
	else:
		#
		# have we already processed it?
		#

		if ligolw_process.doc_includes_process(xmldoc, process_program_name):
			if options.verbose:
				print >>sys.stderr, "warning: %s already processed," % (filename or "stdin"),
			if not options.force:
				if options.verbose:
					print >>sys.stderr, "skipping (use --force to force)"
				continue
			if options.verbose:
				print >>sys.stderr, "continuing by --force"

		#
		# add process metadata to document
		#

		process = ligolw_process.register_to_xmldoc(xmldoc, process_program_name, paramdict, version = __version__, cvs_repository = u"lscsoft", cvs_entry_time = __date__, comment = options.comment)

		#
		# run inspinjfind algorithm
		#

		inspinjfind.ligolw_inspinjfind(
			xmldoc,
			process,
			search,
			snglcomparefunc,
			nearcoinccomparefunc,
			end_time_bisect_window = 1.5 * options.time_window,
			verbose = options.verbose
		)

		#
		# close out the process metadata
		#

		ligolw_process.set_process_end_time(process)

	#
	# done
	#

	ligolw_utils.write_filename(xmldoc, filename, verbose = options.verbose, gz = (filename or "stdout").endswith(".gz"))
	xmldoc.unlink()
	lsctables.reset_next_ids(lsctables.TableByName.values())
