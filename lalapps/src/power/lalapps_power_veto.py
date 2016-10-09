#
# Copyright (C) 2007  Kipp Cannon
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


from glue import segments
from glue.ligolw import dbtables
from glue.ligolw import utils
from glue.ligolw.utils import segments as ligolw_segments
from lalburst import git_version
from lalburst import SnglBurstUtils


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
		usage = "%prog [options] --veto-segments-db=filename files ...",
		description = "%prog probably does something..."
	)
	parser.add_option("--no-vacuum", action = "store_true", default = False, help = "Don't vacuum the database after vetos.  Vacuuming reduces the database size and improves the efficiency of indices, but can take a while so you might want to skip it.")
	parser.add_option("-t", "--tmp-space", metavar = "path", help = "Path to a directory suitable for use as a work area while manipulating the database file.  The database file will be worked on in this directory, and then moved to the final location when complete.  This option is intended to improve performance when running in a networked environment, where there might be a local disk with higher bandwidth than is available to the filesystem on which the final output will reside.")
	parser.add_option("--veto-segments-db", metavar = "filename", help = "Load veto segments from this SQLite DB file (required).")
	parser.add_option("--veto-segments-name", metavar = "name", default = "sngl_burst_veto", help = "Use the segment lists named this to define the veto times (default = \"sngl_burst_veto\").")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()

	if options.veto_segments_db is None:
		raise ValueError("missing required --veto-segments-db option")

	return options, (filenames or [None])


#
# =============================================================================
#
#                                Load Segments
#
# =============================================================================
#


def load_segments(filename, name, verbose = False):
	if verbose:
		print >>sys.stderr, "loading \"%s\" segments ... " % name,
	connection = sqlite3.connect(filename)
	segs = ligolw_segments.segmenttable_get_by_name(dbtables.get_xml(connection), name).coalesce()
	connection.close()
	if verbose:
		print >>sys.stderr, "done."
		for ifo in segs:
			print >>sys.stderr, "loaded %d veto segment(s) for %s totalling %g s" % (len(segs[ifo]), ifo, float(abs(segs[ifo])))
	return segs


#
# =============================================================================
#
#                              Excess Power Veto
#
# =============================================================================
#


def apply_excess_power_veto(contents, veto_segs, verbose = False):
	if verbose:
		SnglBurstUtils.summarize_coinc_database(contents)

	def sngl_burst_is_vetoed(ifo, start, start_ns, duration, veto_segs = veto_segs):
		start = dbtables.lsctables.LIGOTimeGPS(start, start_ns)
		return ifo in veto_segs and veto_segs[ifo].intersects_segment(segments.segment(start, start + duration))
	contents.connection.create_function("sngl_burst_is_vetoed", 4, sngl_burst_is_vetoed)
	cursor = contents.connection.cursor()

	#
	# Delete burst <--> burst coincs containing vetoed bursts
	#

	if verbose:
		print >>sys.stderr, "applying excess power event veto strategy:"
		print >>sys.stderr, "\tremoving vetoed burst <--> burst coincs ..."
	cursor.execute("""
DELETE FROM
	coinc_event
WHERE
	coinc_def_id == ?
	AND coinc_event_id IN (
		SELECT
			coinc_event_id
		FROM
			coinc_event_map
			JOIN sngl_burst ON (
				coinc_event_map.table_name == 'sngl_burst'
				AND coinc_event_map.event_id == sngl_burst.event_id
			)
		WHERE
			sngl_burst_is_vetoed(sngl_burst.ifo, sngl_burst.start_time, sngl_burst.start_time_ns, sngl_burst.duration)
	)
	""", (contents.bb_definer_id,))

	#
	# Delete sim <--> coinc coincs pointing to deleted coincs
	#

	if contents.sc_definer_id is not None:
		if verbose:
			print >>sys.stderr, "\tremoving vetoed sim <--> coinc coincs ..."
		cursor.execute("""
DELETE FROM
	coinc_event
WHERE
	coinc_def_id == ?
	AND coinc_event_id IN (
		SELECT
			coinc_event_id
		FROM
			coinc_event_map
		WHERE
			table_name == 'coinc_event'
			AND event_id NOT IN (
				SELECT
					coinc_event_id
				FROM 
					coinc_event
			)
	)
		""", (contents.sc_definer_id,))

	#
	# Now that we no longer need to form links through the
	# coinc_def_map table, delete vetoed burst rows from it
	#

	if verbose:
		print >>sys.stderr, "\tremoving vetoed bursts from coinc_def_map table ..."
	cursor.execute("""
DELETE FROM
	coinc_event_map
WHERE
	coinc_event_map.table_name == 'sngl_burst'
	AND coinc_event_map.event_id IN (
		SELECT
			event_id
		FROM
			sngl_burst
		WHERE
			sngl_burst_is_vetoed(ifo, start_time, start_time_ns, duration)
	)
	""")

	#
	# The coinc_def_map table no longer contains any rows pointing to
	# vetoed bursts, update the event counts in sim <--> burst coincs
	# and delete any for which this is now 0
	#

	if contents.sb_definer_id is not None:
		if verbose:
			print >>sys.stderr, "\tupdating sim <--> burst event counts ..."
		cursor.execute("""
UPDATE
	coinc_event
SET
	nevents = (SELECT COUNT(*) FROM coinc_event_map WHERE coinc_event_map.coinc_event_id == coinc_event.coinc_event_id AND coinc_event_map.table_name == 'sngl_burst')
WHERE
	coinc_event.coinc_def_id == ?
		""", (contents.sb_definer_id,))
		if verbose:
			print >>sys.stderr, "\tremoving empty sim <--> burst coincs ..."
		cursor.execute("""
DELETE FROM
	coinc_event
WHERE
	coinc_def_id == ?
	AND nevents == 0
		""", (contents.sb_definer_id,))

	#
	# We have removed all the rows from the coinc_event table that will
	# be removed, so now remove multi_burst and coinc_event_map rows
	# that point to non-existent coinc_events
	#

	if verbose:
		print >>sys.stderr, "\ttrimming coinc_event_map table ..."
	cursor.execute("""
DELETE FROM
	coinc_event_map
WHERE
	coinc_event_map.coinc_event_id NOT IN (
		SELECT
			coinc_event_id
		FROM 
			coinc_event
	)
	""")

	if verbose:
		print >>sys.stderr, "\ttrimming multi_burst table ..."
	cursor.execute("""
DELETE FROM
	multi_burst
WHERE
	multi_burst.coinc_event_id NOT IN (
		SELECT
			coinc_event_id
		FROM 
			coinc_event
	)
	""")

	#
	# Tables are now consistent, and contain no coincs involving a
	# vetoed burst
	#

	if verbose:
		print >>sys.stderr, "\tdone (excess power event vetos)."


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
# Load veto segments.
#


veto_segs = load_segments(options.veto_segments_db, options.veto_segments_name, verbose = options.verbose)


#
# Iterate over files
#


for n, filename in enumerate(filenames):
	#
	# Open the database file.
	#

	if options.verbose:
		print >>sys.stderr, "%d/%d: %s" % (n + 1, len(filenames), filename),

	working_filename = dbtables.get_connection_filename(filename, tmp_path = options.tmp_space, verbose = options.verbose)
	connection = sqlite3.connect(working_filename)
	connection.execute("PRAGMA synchronous = OFF;")
	connection.execute("PRAGMA temp_store_directory = '%s';" % dbtables.tempfile.gettempdir())

	#
	# Apply vetoes
	#

	apply_excess_power_veto(SnglBurstUtils.CoincDatabase(connection, "lalapps_power"), veto_segs, verbose = options.verbose)

	#
	# Clean up
	#

	if options.verbose:
		print >>sys.stderr, "committing ..."
	connection.commit()
	if not options.no_vacuum:
		if options.verbose:
			print >>sys.stderr, "vacuuming ..."
		connection.cursor().execute("VACUUM;")
	connection.close()
	del connection
	dbtables.put_connection_filename(filename, working_filename, verbose = options.verbose)
	if options.verbose:
		print >>sys.stderr, "done (%s)." % filename


#
# Done
#


if options.verbose:
	print >>sys.stderr, "done."
