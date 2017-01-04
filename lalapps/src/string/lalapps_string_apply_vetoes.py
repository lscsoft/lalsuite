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


from optparse import OptionParser
import sqlite3
import sys


from glue.ligolw import dbtables
from glue.ligolw.utils import segments as ligolwsegments
from lalburst import git_version


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
		description = "%prog constructs or clears and then populates the string_vetoed_sngl table providing a list of the sngl_burst event_ids of vetoed triggers.  This table is used by other post-processing tools like lalapps_string_meas_likelihood to determine which triggers have been vetoed and which haven't."
	)
	parser.add_option("-t", "--tmp-space", metavar = "path", help = "Path to a directory suitable for use as a work area while manipulating the database file.  The database file will be worked on in this directory, and then moved to the final location when complete.  This option is intended to improve performance when running in a networked environment, where there might be a local disk with higher bandwidth than is available to the filesystem on which the final output will reside.")
	parser.add_option("-n", "--vetoes-name", metavar = "name", default = "vetoes", help = "Set the name of the segment lists to use as vetoes (default = \"vetoes\").")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()

	if not filenames:
		raise ValueError("no input files!")

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
# Iterate over files
#


def create_string_sngl_is_vetoed_function(connection, veto_segments_name = None):
	"""
	Creates a function named string_sngl_is_vetoed in the database at
	connection.  The function accepts three parameters --- the
	instrument name, and the integer and integer nanoseconds components
	of a time --- and returns true if the instrument is vetoed at that
	time or false otherwise.  veto_segments_name sets the name of the
	segment lists used to define the vetoes.

	If veto_segments_name is None then a no-op function is created that
	always returns False.

	Note:  this funtion requires glue.ligolw.dbtables and
	glue.ligolw.utils.segments to be imported as dbtables and
	ligolwsegments respectively.
	"""
	if veto_segments_name is None:
		connection.create_function("string_sngl_is_vetoed", 3, lambda instrument, peak_time, peak_time_ns: False)
		return
	xmldoc = dbtables.get_xml(connection)
	seglists = ligolwsegments.segmenttable_get_by_name(xmldoc, options.vetoes_name).coalesce()
	xmldoc.unlink()
	def is_vetoed(instrument, peak_time, peak_time_ns, seglists = seglists):
		return instrument in seglists and dbtables.lsctables.LIGOTimeGPS(peak_time, peak_time_ns) in seglists[instrument]
	connection.create_function("string_sngl_is_vetoed", 3, is_vetoed)


for n, filename in enumerate(filenames):
	#
	# Open the database file.
	#

	if options.verbose:
		print >>sys.stderr, "%d/%d: %s" % (n + 1, len(filenames), filename)

	working_filename = dbtables.get_connection_filename(filename, tmp_path = options.tmp_space, verbose = options.verbose)
	connection = sqlite3.connect(working_filename)

	#
	# Define the is_vetoed() SQL function.
	#

	create_string_sngl_is_vetoed_function(connection, options.vetoes_name)

	#
	# (Re)build table.
	#

	cursor = connection.cursor()
	cursor.execute("""
DROP TABLE IF EXISTS string_vetoed_sngl;
	""")
	cursor.execute("""
CREATE TABLE string_vetoed_sngl (event_id TEXT PRIMARY KEY)
	""")
	cursor.execute("""
INSERT OR REPLACE INTO
	string_vetoed_sngl
SELECT
	sngl_burst.event_id AS event_id
FROM
	sngl_burst
WHERE
	string_sngl_is_vetoed(sngl_burst.ifo, sngl_burst.peak_time, sngl_burst.peak_time_ns)
	""")
	cursor.close()

	#
	# Clean up.

	connection.close()
	dbtables.discard_connection_filename(filename, working_filename, verbose = options.verbose)


#
# Done.
#
