#
# Copyright (C) 2008--2017  Kipp Cannon
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


from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils
from glue.ligolw.utils import process as ligolw_process
from glue.ligolw.utils import segments as ligolw_segments
from glue import segmentsUtils
import lal
from lalinspiral import thinca


lsctables.use_in(ligolw.LIGOLWContentHandler)


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
from lalapps.git_version import date as __date__
from lalapps.git_version import version as __version__


#
# Must use coincidence code's custom sngl_inspiral table row type.
#


lsctables.SnglInspiralTable.RowType = lsctables.SnglInspiral = thinca.SnglInspiral


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
		description = "%prog implements the inspiral coincidence algorithm for use in performing trigger-based multi-instrument searches for gravitational wave events.  The LIGO Light Weight XML files listed on the command line are processed one by one in order, and over-written with the results.  If no files are named, then input is read from stdin and output written to stdout.  Gzipped files will be autodetected on input, if a file's name ends in \".gz\" it will be gzip-compressed on output."
	)
	parser.add_option("-c", "--comment", metavar = "text", help = "Set comment string in process table (default = None).")
	parser.add_option("-f", "--force", action = "store_true", help = "Process document even if it has already been processed.")
	parser.add_option("-m", "--match", metavar = "algorithm", default = "exact", help = "Select the coincidence test.  Allowed values are:  \"exact\" (require exact template match, the default).")
	parser.add_option("-t", "--threshold", metavar = "float", type = "float", help = "Set the coincidence threshold (required).  The meaning is defined by the match algorithm.  For --match=exact this sets the Delta t window in seconds.")
	parser.add_option("--min-instruments", metavar = "number", default = "2", type = "int", help = "Set the minimum number of instruments that must participate in a coincidence (default = 2).  The value must be greater than 0.")
	parser.add_option("--vetoes-name", metavar = "string", default = "vetoes", help = "From the input document, exatract the segment list having this name to use as the veto segments (default = \"vetoes\").  Warning:  if no segments by this name are found in the document then vetoes will not be applied, this is not an error condition.")
	parser.add_option("--coinc-end-time-segment", metavar = "seg", help = "The segment of time to retain coincident triggers from. Uses segmentUtils.from_range_strings() format \"START:END\" for an interval of the form [START,END), \"START:\" for an interval of the form [START,INF), and \":END\" for an interval of the form (-INF,END).")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()

	process_params = options.__dict__.copy()

	#
	# check arguments
	#

	required_options = ["threshold"]
	missing_options = [option for option in required_options if getattr(options, option) is None]
	if missing_options:
		raise ValueError("missing required option(s) %s" % ", ".join("--%s" % option.replace("_", "-") for option in missing_options))
	if options.match not in ("exact",):
		raise ValueError("unrecognized value for --match: \"%s\"" % options.match)
	if options.min_instruments < 1:
		raise ValueError("invalid --min-instruments: \"%s\"" % options.min_instruments)

	if options.coinc_end_time_segment is not None:
		if ',' in options.coinc_end_time_segment:
			raise ValueError("--coinc-end-time-segment may only contain a single segment")
		options.coinc_end_time_segs = segmentsUtils.from_range_strings([options.coinc_end_time_segment], boundtype = lal.LIGOTimeGPS).coalesce()
	else:
		options.coinc_end_time_segs = None

	#
	# done
	#

	return options, process_params, (filenames or [None])


#
# =============================================================================
#
#                             Process Information
#
# =============================================================================
#


process_program_name = u"lalapps_thinca"


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


options, process_params, filenames = parse_command_line()


#
# Select ntuple_comparefunc form
#


if options.coinc_end_time_segs is not None:
	# Custom Ntuple Comparison Function
	def ntuple_comparefunc(events, offset_vector, seg = options.coinc_end_time_segs):
		"""
		Return False (ntuple should be retained) if the end time of
		the coinc is in the segmentlist segs.
		"""
		return thinca.coinc_inspiral_end_time(events, offset_vector) not in seg
else:
	ntuple_comparefunc = thinca.InspiralCoincTables.ntuple_comparefunc


#
# Iterate over files.
#


for n, filename in enumerate(filenames, start = 1):
	#
	# Load the file.
	#

	if options.verbose:
		print >>sys.stderr, "%d/%d:" % (n, len(filenames)),
	xmldoc = ligolw_utils.load_filename(filename, verbose = options.verbose, contenthandler = ligolw.LIGOLWContentHandler)

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

	process = ligolw_process.register_to_xmldoc(
		xmldoc,
		process_program_name,
		process_params,
		comment = options.comment,
		version = __version__,
		cvs_repository = u"lscsoft",
		cvs_entry_time = __date__
	)

	#
	# LAL writes all triggers with event_id = 0.  Detect this and barf
	# on it.
	#

	tbl = lsctables.SnglInspiralTable.get_table(xmldoc)
	assert len(set(tbl.getColumnByName("event_id"))) == len(tbl), "degenerate sngl_inspiral event_id detected"

	#
	# Extract veto segments if present.
	#
	# FIXME:  using the tools in the glue.ligolw.utils.segments module
	# it's not hard to modify the veto segments in the .xml to be just
	# those that intersect the search summary segments.  That way, if
	# multiple documents are inserted into the same database, or merged
	# with ligolw_add, the veto lists will not get duplicated.
	#

	if not ligolw_segments.has_segment_tables(xmldoc):
		if options.verbose:
			print >>sys.stderr, "warning: no segment definitions found, vetoes will not be applied"
		vetoes = None
	elif not ligolw_segments.has_segment_tables(xmldoc, name = options.vetoes_name):
		if options.verbose:
			print >>sys.stderr, "warning: document contains segment definitions but none named \"%s\", vetoes will not be applied" % options.vetoes_name
		vetoes = None
	else:
		vetoes = ligolw_segments.segmenttable_get_by_name(xmldoc, options.vetoes_name).coalesce()

	#
	# Run coincidence algorithm.
	#

	thinca.ligolw_thinca(
		xmldoc,
		process_id = process.process_id,
		coinc_definer_row = thinca.InspiralCoincDef,
		thresholds = options.threshold,
		ntuple_comparefunc = ntuple_comparefunc,
		seglists = None,	# FIXME
		veto_segments = vetoes,
		min_instruments = options.min_instruments,
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
