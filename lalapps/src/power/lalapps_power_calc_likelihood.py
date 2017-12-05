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


from optparse import OptionParser
import sys
import sqlite3


from glue.ligolw import dbtables
from lal.utils import CacheEntry
from lalburst import git_version
from lalburst import burca_tailor
from lalburst import calc_likelihood
from lalburst import SnglBurstUtils
from lalburst.SimBurstUtils import MW_CENTER_J2000_RA_RAD, MW_CENTER_J2000_DEC_RAD


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
		description = "%prog uses likelihood ratio data stored in LIGO Light-Weight XML files to compute likelihood ratio values for excess power coincs in SQLite databases."
	)
	parser.add_option("-c", "--comment", metavar = "text", help = "Set comment string in process table (default = None).")
	parser.add_option("-p", "--program", metavar = "name", help = "Set the name of the program that generated the events as it appears in the process table (required).  The program name is used to extract live time information from the search summary tables in the input files.")
	parser.add_option("--likelihood-data", metavar = "filename", default = [], action = "append", help = "Read likelihood data from this XML file.  (use lalapps_burca_tailor to generate these files)")
	parser.add_option("--likelihood-data-cache", metavar = "filename", help = "Read likelihood data from the XML files described by this LAL cache.  For each trigger file, the live time of the trigger file is established and all likelihood data files whose segments intersect the trigger file's live time are loaded and merged into a single distribution data set.  (use lalapps_burca_tailor to generate these files)")
	parser.add_option("--tmp-space", metavar = "path", help = "Path to a directory suitable for use as a work area while manipulating the database file.  The database file will be worked on in this directory, and then moved to the final location when complete.  This option is intended to improve performance when running in a networked environment, where there might be a local disk with higher bandwidth than is available to the filesystem on which the final output will reside.")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()

	#
	# check and convert and bunch of arguments
	#

	options.likelihood_data = set(options.likelihood_data)
	if (not options.likelihood_data) and (options.likelihood_data_cache is None):
		raise ValueError("must set one of --likelihood-data or --likelihood-data-cache")
	if options.likelihood_data and (options.likelihood_data_cache is not None):
		raise ValueError("cannot set both --likelihood-data and --likelihood-data-cache")
	if options.likelihood_data_cache:
		options.likelihood_data_cache = set([CacheEntry(line) for line in file(options.likelihood_data_cache)])
	else:
		options.likelihood_data_cache = set()
	if options.program is None:
		raise ValueError("missing required argument --program")

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
# How to load likelihood data
#


def load_likelihood_data(filenames, verbose = False):
	distributions, ignored = burca_tailor.EPGalacticCoreCoincParamsDistributions.from_filenames(filenames, u"lalapps_burca_tailor", verbose = verbose)
	distributions.finish()
	return distributions


#
# Iterate over files.
#


cached_likelihood_files = set()


for n, filename in enumerate(filenames):
	#
	# Open the file.
	#


	if options.verbose:
		print >>sys.stderr, "%d/%d: %s" % (n + 1, len(filenames), filename)
	working_filename = dbtables.get_connection_filename(filename, tmp_path = options.tmp_space, verbose = options.verbose)
	connection = sqlite3.connect(working_filename)
	connection.execute("PRAGMA temp_store_directory = '%s';" % dbtables.tempfile.gettempdir())
	database = SnglBurstUtils.CoincDatabase(connection, options.program)
	if options.verbose:
		SnglBurstUtils.summarize_coinc_database(database)


	#
	# Retrieve appropriate likelihood data.
	#


	if options.likelihood_data_cache:
		likelihood_files = set(c.path for c in options.likelihood_data_cache if c.segmentlistdict.intersects(database.seglists))
	else:
		likelihood_files = options.likelihood_data
	if likelihood_files != cached_likelihood_files:
		distributions = load_likelihood_data(likelihood_files, verbose = options.verbose)
		cached_likelihood_files = likelihood_files


	#
	# Run likelihood ratio calculation.
	#


	calc_likelihood.ligolw_burca2(database, distributions, distributions.coinc_params, verbose = options.verbose, params_func_extra_args = (MW_CENTER_J2000_RA_RAD, MW_CENTER_J2000_DEC_RAD))


	#
	# Done with this file.
	#


	connection.close()
	dbtables.put_connection_filename(filename, working_filename, verbose = options.verbose)
