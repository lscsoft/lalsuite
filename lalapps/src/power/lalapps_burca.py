#
# Copyright (C) 2006--2010  Kipp Cannon
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
from glue.ligolw import utils
from glue.ligolw.utils import process as ligolw_process
from lalburst import git_version
from lalburst import ligolw_burca
from pylal import snglcoinc


lsctables.use_in(ligolw.LIGOLWContentHandler)


#
# Use interning row builder to save memory.
#


lsctables.table.RowBuilder = lsctables.table.InterningRowBuilder


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


def parse_thresholdstrings(thresholdstrings):
	"""
	Turn a list of strings of the form
	inst1,inst2=threshold1[,threshold2,...] into a dictionary with
	(inst1, inst2) 2-tuples as keys and the values being the thresholds
	parsed into lists of strings split on the "," character.

	For each pair of instruments present among the input strings, the
	two possible orders are considered independent:  the input strings
	are allowed to contain one set of thresholds for (inst1, inst2),
	and a different set of thresholds for (inst2, inst1).  Be aware
	that no input checking is done to ensure the user has not provided
	duplicate, incompatible, thresholds.  This is considered the
	responsibility of the application program to verify.

	The output dictionary contains threshold sets for both instrument
	orders.  If, for some pair of instruments, the input strings
	specified thresholds for only one of the two possible orders, the
	thresholds for the other order are copied from the one that was
	provided.

	Whitespace is removed from the start and end of all strings.

	A typical use for this function is in parsing command line
	arguments or configuration file entries.

	Example:

	>>> from pylal.snglcoinc import parse_thresholds
	>>> parse_thresholds(["H1,H2=X=0.1,Y=100", "H1,L1=X=.2,Y=100"])
	{('H1', 'H2'): ['X=0.1', 'Y=100'], ('H1', 'L1'): ['X=.2', 'Y=100'], ('H2', 'H1'): ['X=0.1', 'Y=100'], ('L1', 'H1'): ['X=.2', 'Y=100']}
	"""
	thresholds = {}
	for pair, delta in [s.split("=", 1) for s in thresholdstrings]:
		try:
			A, B = [s.strip() for s in pair.split(",")]
		except Exception:
			raise ValueError("cannot parse instruments '%s'" % pair)
		thresholds[(A, B)] = [s.strip() for s in delta.split(",")]
	for (A, B), value in thresholds.items():
		if (B, A) not in thresholds:
			thresholds[(B, A)] = value
	return thresholds


def parse_thresholds(options):
	#
	# parse --thresholds options into a dictionary of instrument pairs
	# and components
	#

	try:
		thresholds = parse_thresholdstrings(options.thresholds)
	except Exception as e:
		raise ValueError("error parsing --thresholds: %s" % str(e))

	#
	# parse the components from --thresholds options
	#

	if options.coincidence_algorithm == "excesspower":
		#
		# excess power does not use adjustable thresholds, but for
		# speed it helps to pre-compute the light travel time
		# between the instruments involved in the analysis
		#

		for pair in thresholds.keys():
			thresholds[pair] = snglcoinc.light_travel_time(*pair)

	elif options.coincidence_algorithm == "stringcusp":
		#
		# parse thresholds into dt values
		#

		try:
			thresholds = dict((instrumentpair, (float(dt),)) for instrumentpair, (dt,) in thresholds.iteritems())
		except Exception as e:
			raise ValueError("error parsing --thresholds: %s" % str(e))

	else:
		#
		# unrecognized coincidence algorithm
		#

		raise ValueError(options.coincidence_algorithm)

	#
	# Done
	#

	return thresholds


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
	parser.add_option("-t", "--thresholds", metavar = "inst1,inst2=[threshold,...]", action = "append", default = [], help = "Set the coincidence algorithm's thresholds for an instrument pair.  For excesspower there are no thresholds.  For stringcusp, each instrument pair has a single threshold setting dt.  One set of thresholds must be provided for each instrument combination that will be compared, even if there are no thresholds to set.")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()

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

	#
	# parse the --thresholds arguments
	#

	options.thresholds = parse_thresholds(options)

	#
	# done
	#

	return options, (filenames or [None])


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


options, filenames = parse_command_line()


#
# For excesspower and stringcusp methods, select the appropriate
# event comparison and book-keeping functions.
#


if options.coinc_segs is not None:
	def coinc_segs_ntuple_comparefunc(events, offset_vector, coinc_segs = options.coinc_segs):
		# sort so we always do arithmetic in the same order
		events = sorted(events, key = lambda event: event.peak)
		# coinc time is SNR-weighted mean of event peak times
		epoch = events[0].peak + offset_vector[events[0].ifo]
		coinc_time = epoch + sum(float(event.peak + offset_vector[event.ifo] - epoch) * event.snr for event in events) / sum(event.snr for event in events)
		return coinc_time not in coinc_segs


if options.coincidence_algorithm == "excesspower":
	EventListType = ligolw_burca.ExcessPowerEventList
	comparefunc = ligolw_burca.ExcessPowerCoincCompare
	if options.coinc_segs is not None:
		ntuple_comparefunc = coinc_segs_ntuple_comparefunc
	else:
		ntuple_comparefunc = lambda *args: False	# keep everything
	CoincTables = ligolw_burca.ExcessPowerCoincTables
	CoincDef = ligolw_burca.ExcessPowerBBCoincDef
elif options.coincidence_algorithm == "stringcusp":
	EventListType = ligolw_burca.StringEventList
	comparefunc = ligolw_burca.StringCoincCompare
	if options.coinc_segs is not None:
		ntuple_comparefunc = lambda *args: coinc_segs_ntuple_comparefunc(*args) or ligolw_burca.StringNTupleCoincCompare(*args)
	else:
		ntuple_comparefunc = ligolw_burca.StringNTupleCoincCompare
	CoincTables = ligolw_burca.StringCuspCoincTables
	CoincDef = ligolw_burca.StringCuspBBCoincDef
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
	xmldoc = utils.load_filename(filename, verbose = options.verbose, contenthandler = ligolw.LIGOLWContentHandler)
	lsctables.table.InterningRowBuilder.strings.clear()

	#
	# Have we already processed it?
	#

	if ligolw_process.doc_includes_process(xmldoc, ligolw_burca.process_program_name):
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

	process = ligolw_burca.append_process(xmldoc, **options.__dict__)

	#
	# Run coincidence algorithm.
	#

	ligolw_burca.ligolw_burca(
		xmldoc = xmldoc,
		process_id = process.process_id,
		EventListType = EventListType,
		CoincTables = CoincTables,
		coinc_definer_row = CoincDef,
		event_comparefunc = comparefunc,
		thresholds = options.thresholds,
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

	utils.write_filename(xmldoc, filename, verbose = options.verbose, gz = (filename or "stdout").endswith(".gz"))
	xmldoc.unlink()
	lsctables.table.reset_next_ids(lsctables.TableByName.values())
