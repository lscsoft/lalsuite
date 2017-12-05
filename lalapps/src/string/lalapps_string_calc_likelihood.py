#
# Copyright (C) 2010  Kipp Cannon
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


import math
import numpy
from optparse import OptionParser
import sqlite3
import sys


from glue.ligolw import dbtables
from lal.utils import CacheEntry
from lalburst import git_version
from lalburst import calc_likelihood
from lalburst import SnglBurstUtils
from lalburst import stringutils


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
		version = "Name: %%prog\n%s" % git_version.verbose_msg
	)
	parser.add_option("-c", "--input-cache", metavar = "filename", help = "Also process the files named in this LAL cache.  See lalapps_path2cache for information on how to produce a LAL cache file.")
	parser.add_option("-l", "--likelihood-file", metavar = "filename", action = "append", help = "Set the name of the likelihood ratio data file to use.  Can be given more than once.")
	parser.add_option("--likelihood-cache", metavar = "filename", help = "Also load the likelihood ratio data files listsed in this LAL cache.  See lalapps_path2cache for information on how to produce a LAL cache file.")
	parser.add_option("-t", "--tmp-space", metavar = "path", help = "Path to a directory suitable for use as a work area while manipulating the database file.  The database file will be worked on in this directory, and then moved to the final location when complete.  This option is intended to improve performance when running in a networked environment, where there might be a local disk with higher bandwidth than is available to the filesystem on which the final output will reside.")
	parser.add_option("--vetoes-name", metavar = "name", help = "Set the name of the segment lists to use as vetoes (default = do not apply vetoes).")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()

	options.likelihood_filenames = []
	if options.likelihood_file is not None:
		options.likelihood_filenames += options.likelihood_file
	if options.likelihood_cache is not None:
		options.likelihood_filenames += [CacheEntry(line).path for line in open(options.likelihood_cache)]
	if not options.likelihood_filenames:
		raise ValueError("no ranking statistic likelihood data files specified")

	if options.input_cache:
		filenames += [CacheEntry(line).path for line in open(options.input_cache)]
	if not filenames:
		raise ValueError("no candidate databases specified")

	return options, filenames


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


options, filenames = parse_command_line()


#
# load likelihood data
#


coincparamsdistributions, likelihood_seglists = stringutils.load_likelihood_data(options.likelihood_filenames, verbose = options.verbose)
if options.verbose:
	print >>sys.stderr, "computing event densities ..."
coincparamsdistributions.finish()


#
# iterate over files
#


for n, filename in enumerate(filenames):
	#
	# Open the database file.
	#

	if options.verbose:
		print >>sys.stderr, "%d/%d: %s" % (n + 1, len(filenames), filename)

	working_filename = dbtables.get_connection_filename(filename, tmp_path = options.tmp_space, verbose = options.verbose)
	connection = sqlite3.connect(working_filename)
	if options.tmp_space is not None:
		dbtables.set_temp_store_directory(connection, options.tmp_space, verbose = options.verbose)

	#
	# Summarize the database.
	#

	contents = SnglBurstUtils.CoincDatabase(connection, live_time_program = "StringSearch", search = "StringCusp", veto_segments_name = options.vetoes_name)
	if options.verbose:
		SnglBurstUtils.summarize_coinc_database(contents)
	if not contents.seglists and options.verbose:
		print >>sys.stderr, "\twarning:  no segments found"

	#
	# Build triangulators.  The timing uncertainty of +/- 8e-5 s was
	# measured with lalapps_string_plot_binj and is essentially
	# identical for H1, H2, L1, and V1.
	#

	triangulators = stringutils.triangulators(dict.fromkeys(contents.instruments, 8e-5))

	#
	# Run likelihood ratio calculation.
	#

	calc_likelihood.ligolw_burca2(contents, coincparamsdistributions, coincparamsdistributions.coinc_params, verbose = options.verbose, params_func_extra_args = (triangulators,))

	#
	# Clean up.
	#

	contents.xmldoc.unlink()
	connection.close()
	dbtables.put_connection_filename(filename, working_filename, verbose = options.verbose)
