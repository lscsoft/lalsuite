##python
#
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


from optparse import OptionParser
import sys


from igwn_ligolw import ligolw
from igwn_ligolw import lsctables
from igwn_ligolw import utils as ligolw_utils
from igwn_ligolw.utils import process as ligolw_process
from lalburst import git_version
from lalburst import burca
from igwn_segments import utils as segmentsUtils


@lsctables.use_in
class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
	pass

process_program_name = "lalburst_coinc"


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


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
		description = "%prog implements the excess power and string cusp coincidence algorithms for use in performing trigger-based multi-instrument searches for gravitational wave events.  The LIGO Light Weight XML files listed on the command line are processed one by one in order, and over-written with the results.  If no files are named, then input is read from stdin and output written to stdout.  Gzipped files will be autodetected on input, if a file's name ends in \".gz\" it will be gzip-compressed on output."
	)
	parser.add_option("-c", "--comment", metavar = "text", help = "Set comment string in process table (default = None).")
	parser.add_option("-f", "--force", action = "store_true", help = "Process document even if it has already been processed.")
	parser.add_option("-a", "--coincidence-algorithm", metavar = "[excesspower|stringcusp]", default = None, help = "Select the coincidence test algorithm to use (required).")
	parser.add_option("-s", "--coincidence-segments", metavar = "start:end[,start:end,...]", help = "Set the GPS segments in which to retain coincidences.  Multiple segments can be specified by separating them with commas.  If either start or end is absent from a segment then the interval is unbounded on that side, for example \"874000000:\" causes all coincidences starting at 874000000 to be retained.  The \"time\" of a coincidence is ambiguous, and is different for different search algorithms, but a deterministic algorithm is used in all cases so the same coincidence of events will always be assigned the same time.  This feature is intended to allow large input files to be analyzed;  the procedure is to make n copies of the file and run n instances of burca specifying disjoint --coincidence-segments for each.")
	parser.add_option("-m", "--min-instruments", metavar = "N", type = "int", default = 2, help = "Set the minimum number of instruments required to form a coincidence (default = 2).")
	parser.add_option("-t", "--threshold", metavar = "threshold", type = "float", help = "Set the coincidence algorithm's threshold.  For excesspower this parameter is not used and must not be set.  For stringcusp, this parameter sets the maximum peak time difference in seconds not including light travel time (which will be added internally).")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()

	paramdict = options.__dict__.copy()

	#
	# check and convert and bunch of arguments
	#

	if options.coincidence_algorithm is None:
		raise ValueError("missing required argument --coincidence-algorithm")
	if options.coincidence_algorithm not in ("excesspower", "stringcusp"):
		raise ValueError("unrecognized --coincidence-algorithm %s" % options.coincidence_algorithm)
	if options.coincidence_segments is not None:
		options.coincidence_segments = segmentsUtils.from_range_strings(options.coincidence_segments.split(","), boundtype = lsctables.LIGOTimeGPS).coalesce()
		if filenames is not None and len(filenames) > 1:
			raise ValueError("refusing to allow use of --coincidence-segments with more than one input file")
	if options.min_instruments < 1:
		raise ValueError("--min-instruments must be >= 1")

	if options.coincidence_algorithm == "stringcusp" and options.threshold is None:
		raise ValueError("--threshold is required for --coincidence-algorithm stringcusp")
	elif options.coincidence_algorithm == "excesspower" and options.threshold is not None:
		raise ValueError("--threshold is meaningless for --coincidence-algorithm excesspower")

	#
	# done
	#

	return options, (filenames or [None]), paramdict


#
# =============================================================================
#
#                                     Main
#
# =============================================================================
#


#
# Command line
#


options, filenames, paramdict = parse_command_line()


#
# For excesspower and stringcusp methods, select the appropriate
# event comparison and book-keeping functions.
#


if options.coincidence_segments is not None:
	def coinc_segs_ntuple_comparefunc(events, offset_vector, min_instruments = options.min_instruments, coinc_segs = options.coincidence_segments):
		# for the purposes of the coinc segs feature, the coinc
		# time is minimum of event peak times.  this is a fast, and
		# guaranteed reproducible definition
		return min(event.peak + offset_vector[event.ifo] for event in events) not in coinc_segs
else:
	def coinc_segs_ntuple_comparefunc(*args):
		return False


if options.coincidence_algorithm == "excesspower":
	coincgen_doubles = burca.ep_coincgen_doubles
	ntuple_comparefunc = coinc_segs_ntuple_comparefunc
	CoincTables = burca.ExcessPowerCoincTables
	CoincDef = burca.ExcessPowerBBCoincDef
elif options.coincidence_algorithm == "stringcusp":
	coincgen_doubles = burca.string_coincgen_doubles
	ntuple_comparefunc = lambda *args: coinc_segs_ntuple_comparefunc(*args) or burca.StringCuspCoincTables.ntuple_comparefunc(*args)
	CoincTables = burca.StringCuspCoincTables
	CoincDef = burca.StringCuspBBCoincDef
else:
	raise Exception("should never get here")


#
# Iterate over files.
#


for n, filename in enumerate(filenames):
	#
	# Load the file.
	#

	if options.verbose:
		print("%d/%d:" % (n + 1, len(filenames)), end=' ', file=sys.stderr)
	xmldoc = ligolw_utils.load_filename(filename, verbose = options.verbose, contenthandler = LIGOLWContentHandler)

	#
	# Have we already processed it?
	#

	if ligolw_process.doc_includes_process(xmldoc, process_program_name):
		if options.verbose:
			print("warning: %s already processed," % (filename or "stdin"), end=' ', file=sys.stderr)
		if not options.force:
			if options.verbose:
				print("skipping", file=sys.stderr)
			continue
		if options.verbose:
			print("continuing by --force", file=sys.stderr)

	#
	# Add an entry to the process table.
	#

	process = ligolw_process.register_to_xmldoc(xmldoc, process_program_name, paramdict,
		comment = options.comment,
		version = __version__,
		cvs_repository = "lscsoft",
		cvs_entry_time = __date__
	)

	#
	# Run coincidence algorithm.
	#


	if options.coincidence_algorithm == "excesspower":
		delta_t = 2. * bruca.ep_coincgen_doubles.get_coincs.max_edge_peak_delta(lsctables.SnglBurstTable.get_table(xmldoc))
	elif options.coincidence_algorithm == "stringcusp":
		delta_t = options.threshold
	else:
		raise Exception("should never get here")
	burca.burca(
		xmldoc = xmldoc,
		process_id = process.process_id,
		coincgen_doubles = coincgen_doubles,
		CoincTables = CoincTables,
		coinc_definer_row = CoincDef,
		delta_t = delta_t,
		ntuple_comparefunc = ntuple_comparefunc,
		min_instruments = options.min_instruments,
		verbose = options.verbose
	)

	#
	# Close out the process table.
	#

	process.set_end_time_now()

	#
	# Write back to disk, and clean up.
	#

	ligolw_utils.write_filename(xmldoc, filename, verbose = options.verbose)
	xmldoc.unlink()
	lsctables.reset_next_ids(lsctables.TableByName.values())
