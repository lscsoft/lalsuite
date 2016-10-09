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
LIGO Light-Weight XML Coincidence Analysis Front End.
"""


from optparse import OptionParser
import sys


from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils
from lalburst import git_version
from lalburst import cafe


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
		usage = "%prog --time-slides time_slides_filename [options] [filename ...]",
		description = "%prog is a generic coincidence analysis front-end for trigger-based searches for gravitational wave events.  This program answers the question \"given that I have a collection of files containing event lists, and given that I wish to perform a coincidence analysis by time-shifting the events to simulate a background, what are the smallest groups of trigger files such that if each group is analyzed separately no coincidences will be missed?\"  The inputs consist of one or more LAL cache files describing the collection of trigger files to be analyzed, and a LIGO Light Weight XML file containing a time_slide table describing the time shifts to be applied to the triggers.  If no cache files are named on the command line, then the cache is read from stdin.  The output is a collection of LAL cache files, one each for the file groups identified from the input.  See also lalapps_ll2cache (a program that constructs a LAL cache file from a list of LIGO Light Weight XML trigger files by parsing each file's search_summary table), and ligolw_add (a program that combines LIGO Light Weight XML trigger files, and is capable of reading its input file list from a LAL cache)."
	)
	parser.add_option("-s", "--single-instrument", action = "store_true", help = "Select single instrument mode.  In this mode, after grouping the files as usual, for each group a separate LAL cache file is written for each instrument, rather than a single cache listing all input files.")
	parser.add_option("-t", "--time-slides", metavar = "filename", help = "Read the time slide table from this file.")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	parser.add_option("-b", "--base", metavar = "base", default = "cafe_", help = "Set base for output caches (default = \"cafe_\").")
	options, cachenames = parser.parse_args()

	if not options.time_slides:
		raise ValueError("--time-slides required")

	return options, (cachenames or [None])


#
# =============================================================================
#
#                                     Main
#
# =============================================================================
#


options, cachenames = parse_command_line()


cache = []
for filename in cachenames:
	cache.extend(cafe.load_cache(filename, options.verbose))


@lsctables.use_in
class LIGOLWContentHandler(lsctables.ligolw.LIGOLWContentHandler):
	pass

seglists, outputcaches = cafe.cafe(cache, lsctables.TimeSlideTable.get_table(ligolw_utils.load_filename(options.time_slides, verbose = options.verbose, contenthandler = LIGOLWContentHandler)).as_dict().values(), options.verbose)
instruments = set(seglists.keys())


if options.single_instrument:
	cafe.write_single_instrument_caches(options.base, outputcaches, instruments, options.verbose)
else:
	cafe.write_caches(options.base, outputcaches, instruments, options.verbose)
