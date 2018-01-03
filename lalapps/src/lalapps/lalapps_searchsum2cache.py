#
# Copyright (C) 2006  Kipp Cannon
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
Build a LAL cache from a list of LIGO LW XML files containing search
summary tables.
"""


import glob
from optparse import OptionParser
import os
import sys

from glue.lal import CacheEntry
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
#from lalapps import git_version


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
#__version__ = "git id %s" % git_version.id
#__date__ = git_version.date


#
# =============================================================================
#
#                                 Command Line
#
# =============================================================================
#


def parse_command_line():
	parser = OptionParser(
		#version = "Name: %%prog\n%s" % git_version.verbose_msg,
		usage = "%prog [options] filenames ...",
		description = "Generates a LAL format cache file describing a collection of LIGO light-weight XML files.  The cache is constructed by parsing the search_summary table in each file to extract the instruments and time each file spans.  To allow long file lists to be processed, the filenames are interpreted as shell patterns (wildcard expansion is performed)."
	)
	parser.add_option("--description", metavar = "string", help = "Set all descriptions to this string.  Use \"-\" for no description.  If not given then the description will be extracted from the search summary rows, and if the search summary rows do not provide a unique description an error is raised.")
	parser.add_option("--observatory", metavar = "string", help = "Set all observatories to this string.  Use \"-\" for no observatory.  If not given then the union of the instruments from the search summary rows will be used to construct an \"observatory\" string.")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	parser.add_option("-o", "--output", metavar = "filename", help = "Write output to this file (default = stdout).")
	parser.add_option("-p", "--program", metavar = "name", action = "append", help = "Obtain instruments, starts, durations, and descriptions from the search summary rows for this program (default = use all search summary rows).  Can be given multiple times to select rows from more than one program.")
	options, filenames = parser.parse_args()

	if options.output:
		options.output = file(options.output, "w")
	else:
		options.output = sys.stdout

	if options.program is not None:
		options.program = set(options.program)

	if filenames is None:
		filenames = []
	else:
		filenames = [filename for g in filenames for filename in glob.glob(g)]

	return options, filenames


#
# =============================================================================
#
#                                    Input
#
# =============================================================================
#


def element_filter(name, attrs):
	"""
	Return True if name & attrs describe a search summary table or a
	process table.
	"""
	return lsctables.IsTableProperties(lsctables.SearchSummaryTable, name, attrs) or lsctables.IsTableProperties(lsctables.ProcessTable, name, attrs)


class ContentHandler(ligolw.PartialLIGOLWContentHandler):
	def __init__(self, xmldoc):
		ligolw.PartialLIGOLWContentHandler.__init__(self, xmldoc, element_filter)


try:
	lsctables.use_in(ContentHandler)
except AttributeError:
	# old glue did not allow .use_in().
	# FIXME:  remove when we can require the latest version of glue
	pass


#
# =============================================================================
#
#                                     Main
#
# =============================================================================
#


options, filenames = parse_command_line()


for n, filename in enumerate(filenames):
	# load document and extract search summary table
	if options.verbose:
		print >>sys.stderr, "%d/%d:" % (n + 1, len(filenames)),
	xmldoc = utils.load_filename(filename, verbose = options.verbose, contenthandler = ContentHandler)
	searchsumm = table.get_table(xmldoc, lsctables.SearchSummaryTable.tableName)

	# extract process_ids for the requested program
	if options.program is not None:
		process_table = table.get_table(xmldoc, lsctables.ProcessTable.tableName)
		process_ids = reduce(lambda a, b: a | b, map(process_table.get_ids_by_program, options.program))
	else:
		process_ids = None

	# extract segment lists
	seglists = searchsumm.get_out_segmentlistdict(process_ids).coalesce()
	if not seglists:
		raise ValueError, "%s: no matching rows found in search summary table" % filename
	if None in seglists:
		if options.program is not None:
			raise ValueError, "%s: null value in ifos column in search_summary table" % filename
		raise ValueError, "%s: null value in ifos column in search_summary table, try using --program" % filename

	# extract observatory
	observatory = (options.observatory and options.observatory.strip()) or "+".join(sorted(seglists))

	# extract description
	if options.description:
		description = options.description
	else:
		if process_ids is None:
			description = set(searchsumm.getColumnByName("comment"))
		else:
			description = set(row.comment for row in searchsumm if row.process_id in process_ids)
		if len(description) < 1:
			raise ValueError, "%s: no matching rows found in search summary table" % filename
		if len(description) > 1:
			raise ValueError, "%s: comments in matching rows of search summary table are not identical" % filename
		description = description.pop().strip() or None

	# set URL
	url = "file://localhost" + os.path.abspath(filename)

	# write cache entry
	print >>options.output, str(CacheEntry(observatory, description, seglists.extent_all(), url))

	# allow garbage collection
	xmldoc.unlink()
