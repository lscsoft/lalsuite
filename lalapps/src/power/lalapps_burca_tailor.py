#
# Copyright (C) 2007-2010  Kipp Cannon
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
from optparse import OptionParser
import sqlite3
import string
import sys

from lal.utils import CacheEntry

from glue import segments
from glue.ligolw import dbtables
from glue.ligolw import utils
from lalburst import ligolw_burca_tailor
from lalburst import SnglBurstUtils
from lalburst.SimBurstUtils import MW_CENTER_J2000_RA_RAD, MW_CENTER_J2000_DEC_RAD
from lalburst import git_version


# characters allowed to appear in the description string
T010150_letters = set(string.ascii_lowercase + string.ascii_uppercase + string.digits + "_+#")


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
		usage = "%prog [options] [filename ...]",
		description = "%prog analyzes a collection of SQLite3 database files containing ligolw_burca outputs, and measures probability distributions for a variety of parameters computed from the coincidences therein.  The distributions are written to a likelihood data file in XML format, which can be used by ligolw_burca for the excesspower2 algorithm in which a second pass assigns likelihoods to each coincidence.  The command line arguments are used to provide shell patterns for the files from which to obtain injection and backgroun coincidences.  If file names are given on the command line following the arguments, then likelihood data is loaded from those files and added to the output."
	)
	parser.add_option("--add-from", metavar = "filename", default = [], action = "append", help = "Also add likelihood data from this XML file.")
	parser.add_option("--add-from-cache", metavar = "filename", help = "Also add likelihood data from all XML files listed in this LAL cache.")
	parser.add_option("-o", "--output", metavar = "filename", default = None, help = "Set the name of the likelihood control file to write (default = stdout).")
	parser.add_option("-t", "--tmp-space", metavar = "path", help = "Path to a directory suitable for use as a work area while manipulating the database file.  The database file will be worked on in this directory, and then moved to the final location when complete.  This option is intended to improve performance when running in a networked environment, where there might be a local disk with higher bandwidth than is available to the filesystem on which the final output will reside.")
	parser.add_option("--T010150", metavar = "description", default = None, help = "Write the output to a file whose name is compatible with the file name format described in LIGO-T010150-00-E, \"Naming Convention for Frame Files which are to be Processed by LDAS\".  The description string will be used to form the second field in the file name.")
	parser.add_option("-p", "--live-time-program", metavar = "program", default = "lalapps_power", help = "Program from which to draw the livetime segments. (Necessary in case of giving --T010150.")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()

	if options.T010150 is not None:
		if options.output is not None:
			raise ValueError("cannot set both --T010150 and --output")
		if options.T010150 == "":
			options.T010150 = "EXCESSPOWER_LIKELIHOOD"
		elif set(options.T010150) - T010150_letters:
			raise ValueError("invalid characters in description \"%s\"" % options.T010150)

	if options.add_from_cache:
		options.add_from += [CacheEntry(line).path for line in file(options.add_from_cache)]

	return options, filenames


#
# =============================================================================
#
#                                     Main
#
# =============================================================================
#


#
# Command line.
#


options, filenames = parse_command_line()


#
# Coinc params
#


distributions = ligolw_burca_tailor.EPGalacticCoreCoincParamsDistributions()
segs = segments.segmentlistdict()


#
# Load pre-computed likelihood data.
#


if options.add_from:
	c, s = distributions.from_filenames(options.add_from, u"ligolw_burca_tailor", verbose = options.verbose)
	distributions += c
	segs |= s
	del c
	del s


#
# Iterate over files
#


for n, filename in enumerate(filenames):
	#
	# Open the database file.
	#

	if options.verbose:
		print >>sys.stderr, "%d/%d: %s" % (n + 1, len(filenames), filename)

	working_filename = dbtables.get_connection_filename(filename, tmp_path = options.tmp_space, verbose = options.verbose)
	connection = sqlite3.connect(working_filename)
	connection.execute("PRAGMA synchronous = OFF;")
	connection.execute("PRAGMA temp_store_directory = '%s';" % dbtables.tempfile.gettempdir())

	#
	# Summarize the database.
	#

	database = SnglBurstUtils.CoincDatabase(connection, options.live_time_program)
	if options.verbose:
		SnglBurstUtils.summarize_coinc_database(database)
	segs |= database.seglists

	#
	# Record statistics.  Assume all files with sim_burst tables are
	# the outputs of injection runs, and others aren't.
	#

	if database.sim_burst_table is None:
		# iterate over burst<-->burst coincs
		for is_background, events, offsetvector in ligolw_burca_tailor.get_noninjections(database):
			params = distributions.coinc_params(events, offsetvector, MW_CENTER_J2000_RA_RAD, MW_CENTER_J2000_DEC_RAD)
			if params is not None:
				if is_background:
					distributions.add_background(params)
				else:
					distributions.add_zero_lag(params)
	else:
		# iterate over burst<-->burst coincs matching injections
		# "exactly"
		for sim, events, offsetvector in ligolw_burca_tailor.get_injections(database):
			params = distributions.coinc_params(events, offsetvector, MW_CENTER_J2000_RA_RAD, MW_CENTER_J2000_DEC_RAD)
			if params is not None:
				distributions.add_injection(params)

	#
	# Clean up.
	#

	del database, connection
	dbtables.discard_connection_filename(filename, working_filename, verbose = options.verbose)


#
# Output.
#


def T010150_basename(description, seglists):
	seg = seglists.extent_all()
	return "%s-%s-%s-%s" % ("+".join(sorted(seglists.keys())), description, str(int(seg[0])), str(int(math.ceil(abs(seg)))))


if options.T010150:
	filename = T010150_basename(options.T010150, segs) + ".xml.gz"
else:
	filename = options.output


xmldoc = ligolw_burca_tailor.gen_likelihood_control(distributions, segs)
utils.write_filename(xmldoc, filename, verbose = options.verbose, gz = (filename or "stdout").endswith(".gz"))
