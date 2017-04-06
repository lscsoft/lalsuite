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


import glob
from optparse import OptionParser
import sys


from glue import segmentsUtils
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils
from glue.ligolw.utils import process as ligolw_process
from lalburst import git_version
from lalburst import burca
from pylal import snglcoinc


@lsctables.use_in
class LIGOLWContentHandler(ligolw.LIGOLWContentHandler):
	pass

lsctables.SnglBurstTable.RowType = burca.SnglBurst

process_program_name = "lalapps_burca"


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


class Thresholds(object):
	"""
	Fake dictionary that returns the same thing for all keys.
	"""
	def __init__(self, val):
		self.val = val

	def __setitem__(self, key, val):
		self.val = val

	def __getitem__(self, key):
		return self.val


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
	parser.add_option("-t", "--threshold", metavar = "threshold", default = None, help = "Set the coincidence algorithm's threshold.  For excesspower this parameter is not used.  For stringcusp, this parameter sets the maximum peak time difference not including light travel time (which will be added internally).")
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
		options.coinc_segs = segmentsUtils.from_range_strings(options.coincidence_segments.split(","), boundtype = lsctables.LIGOTimeGPS).coalesce()
		if filenames is not None and len(filenames) > 1:
			raise ValueError("refusing to allow use of --coincidence-segments with more than one input file")
	else:
		options.coinc_segs = None
	if options.min_instruments < 1:
		raise ValueError("--min-instruments must be >= 1")

	#
	# parse the --thresholds arguments
	#

	if options.coincidence_algorithm == "stringcusp":
		if options.threshold is None:
			raise ValueError("--threshold is required for --coincidence-algorithm stringcusp")
		options.threshold = Thresholds(float(options.threshold))

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


if options.coinc_segs is not None:
	def coinc_segs_ntuple_comparefunc(events, offset_vector, min_instruments = options.min_instruments, coinc_segs = options.coinc_segs):
		# for the purposes of the coinc segs feature, the coinc
		# time is minimum of event peak times.  this is a fast, and
		# guaranteed reproducible definition
		return len(events) < min_instruments or min(event.peak for event in events) not in coinc_segs
else:
	def coinc_segs_ntuple_comparefunc(events, offset_vector, min_instruments = options.min_instruments):
		return len(events) < min_instruments


if options.coincidence_algorithm == "excesspower":
	EventListType = burca.ExcessPowerEventList
	comparefunc = burca.ExcessPowerCoincCompare
	ntuple_comparefunc = coinc_segs_ntuple_comparefunc
	CoincTables = burca.ExcessPowerCoincTables
	CoincDef = burca.ExcessPowerBBCoincDef
elif options.coincidence_algorithm == "stringcusp":
	EventListType = burca.StringEventList
	comparefunc = burca.StringCoincCompare
	ntuple_comparefunc = lambda *args: coinc_segs_ntuple_comparefunc(*args) or burca.StringNTupleCoincCompare(*args)
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
		print >>sys.stderr, "%d/%d:" % (n + 1, len(filenames)),
	xmldoc = ligolw_utils.load_filename(filename, verbose = options.verbose, contenthandler = LIGOLWContentHandler)

	#
	# Have we already processed it?
	#

	if ligolw_process.doc_includes_process(xmldoc, process_program_name):
		if options.verbose:
			print >>sys.stderr, "warning: %s already processed," % (filename or "stdin"),
		if not options.force:
			if options.verbose:
				print >>sys.stderr, "skipping"
			continue
		if options.verbose:
			print >>sys.stderr, "continuing by --force"

	#
	# Add an entry to the process table.
	#

	process = ligolw_process.register_to_xmldoc(xmldoc, process_program_name, paramdict,
		comment = options.comment,
		version = __version__,
		cvs_repository = u"lscsoft",
		cvs_entry_time = __date__
	)

	#
	# Run coincidence algorithm.
	#

	burca.burca(
		xmldoc = xmldoc,
		process_id = process.process_id,
		EventListType = EventListType,
		CoincTables = CoincTables,
		coinc_definer_row = CoincDef,
		event_comparefunc = comparefunc,
		thresholds = options.threshold,
		ntuple_comparefunc = ntuple_comparefunc,
		verbose = options.verbose
	)

	#
	# Close out the process table.
	#

	ligolw_process.set_process_end_time(process)

	#
	# Write back to disk, and clean up.
	#

	ligolw_utils.write_filename(xmldoc, filename, verbose = options.verbose, gz = (filename or "stdout").endswith(".gz"))
	xmldoc.unlink()
	lsctables.reset_next_ids(lsctables.TableByName.values())
