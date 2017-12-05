# Copyright (C) 2005-2013,2016,2017  Kipp Cannon
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


import itertools
import math
import matplotlib
matplotlib.rcParams.update({
	"font.size": 8.0,
	"axes.titlesize": 10.0,
	"axes.labelsize": 10.0,
	"xtick.labelsize": 8.0,
	"ytick.labelsize": 8.0,
	"legend.fontsize": 8.0,
	"figure.dpi": 600,
	"savefig.dpi": 600,
	"text.usetex": True	# render all text with TeX
})
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import re
import sys


from glue.ligolw import ilwd
from glue.ligolw import lsctables
from glue.ligolw import dbtables
from glue.ligolw.utils import search_summary as ligolw_search_summary
from glue.ligolw.utils import segments as ligolw_segments
from glue.offsetvector import offsetvector


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
from git_version import date as __date__
from git_version import version as __version__


#
# =============================================================================
#
#                                   Database
#
# =============================================================================
#


class CoincDatabase(object):
	def __init__(self, connection, live_time_program, search = "excesspower", veto_segments_name = None):
		"""
		Compute and record some summary information about the
		database.  Call this after all the data has been inserted,
		and before you want any of this information.
		"""

		self.connection = connection
		self.xmldoc = dbtables.get_xml(connection)

		# find the tables
		try:
			self.sngl_burst_table = lsctables.SnglBurstTable.get_table(self.xmldoc)
		except ValueError:
			self.sngl_burst_table = None
		try:
			self.sim_burst_table = lsctables.SimBurstTable.get_table(self.xmldoc)
		except ValueError:
			self.sim_burst_table = None
		try:
			self.coinc_def_table = lsctables.CoincDefTable.get_table(self.xmldoc)
			self.coinc_table = lsctables.CoincTable.get_table(self.xmldoc)
			self.time_slide_table = lsctables.TimeSlideTable.get_table(self.xmldoc)
		except ValueError:
			self.coinc_def_table = None
			self.coinc_table = None
			self.time_slide_table = None
		try:
			self.multi_burst_table = lsctables.MultiBurstTable.get_table(self.xmldoc)
		except ValueError:
			self.multi_burst_table = None

		# get the segment lists
		self.seglists = ligolw_search_summary.segmentlistdict_fromsearchsummary(self.xmldoc, live_time_program).coalesce()
		self.instruments = set(self.seglists.keys())
		if veto_segments_name is not None:
			self.vetoseglists = ligolw_segments.segmenttable_get_by_name(self.xmldoc, veto_segments_name).coalesce()
		else:
			self.vetoseglists = ligolw_segments.segments.segmentlistdict()

		# determine a few coinc_definer IDs
		# FIXME:  don't hard-code the numbers
		if self.coinc_def_table is not None:
			try:
				self.bb_definer_id = self.coinc_def_table.get_coinc_def_id(search, 0, create_new = False)
			except KeyError:
				self.bb_definer_id = None
			try:
				self.sb_definer_id = self.coinc_def_table.get_coinc_def_id(search, 1, create_new = False)
			except KeyError:
				self.sb_definer_id = None
			try:
				self.sce_definer_id = self.coinc_def_table.get_coinc_def_id(search, 2, create_new = False)
			except KeyError:
				self.sce_definer_id = None
			try:
				self.scn_definer_id = self.coinc_def_table.get_coinc_def_id(search, 3, create_new = False)
			except KeyError:
				self.scn_definer_id = None
		else:
			self.bb_definer_id = None
			self.sb_definer_id = None
			self.sce_definer_id = None
			self.scn_definer_id = None


	def get_noninjections(self):
		"""
		Generator function to return

			is_background, event_list, offsetvector

		tuples by querying the coinc_event and sngl_burst tables in
		the database.  Only coincs corresponding to
		sngl_burst<-->sngl_burst coincs will be retrieved.
		"""
		cursor = self.connection.cursor()
		for coinc_event_id, time_slide_id in self.connection.cursor().execute("""
	SELECT
		coinc_event_id,
		time_slide_id
	FROM
		coinc_event
	WHERE
		coinc_def_id == ?
		""", (self.bb_definer_id,)):
			rows = [(self.sngl_burst_table.row_from_cols(row), row[-1]) for row in cursor.execute("""
	SELECT
		sngl_burst.*,
		time_slide.offset
	FROM
		coinc_event_map
		JOIN sngl_burst ON (
			coinc_event_map.table_name == 'sngl_burst'
			AND sngl_burst.event_id == coinc_event_map.event_id
		)
		JOIN time_slide ON (
			time_slide.instrument == sngl_burst.ifo
		)
	WHERE
		coinc_event_map.coinc_event_id == ?
		AND time_slide.time_slide_id == ?
			""", (coinc_event_id, time_slide_id))]
			offsets = offsetvector((event.ifo, offset) for event, offset in rows)
			yield any(offsets.values()), [event for event, offset in rows], offsets
		cursor.close()


	def get_injections(self):
		"""
		Generator function to return

			sim, event_list, offsetvector

		tuples by querying the sim_burst, coinc_event and
		sngl_burst tables in the database.  Only coincs
		corresponding to "exact" sim_burst<-->coinc_event coincs
		will be retrieved.
		"""
		cursor = self.connection.cursor()
		for values in self.connection.cursor().execute("""
	SELECT
		sim_burst.*,
		burst_coinc_event_map.event_id
	FROM
		sim_burst
		JOIN coinc_event_map AS sim_coinc_event_map ON (
			sim_coinc_event_map.table_name == 'sim_burst'
			AND sim_coinc_event_map.event_id == sim_burst.simulation_id
		)
		JOIN coinc_event AS sim_coinc_event ON (
			sim_coinc_event.coinc_event_id == sim_coinc_event_map.coinc_event_id
		)
		JOIN coinc_event_map AS burst_coinc_event_map ON (
			burst_coinc_event_map.coinc_event_id == sim_coinc_event_map.coinc_event_id
			AND burst_coinc_event_map.table_name == 'coinc_event'
		)
	WHERE
		sim_coinc_event.coinc_def_id == ?
		""", (self.sce_definer_id,)):
			# retrieve the injection and the coinc_event_id
			sim = self.sim_burst_table.row_from_cols(values)
			coinc_event_id = values[-1]

			# retrieve the list of the sngl_bursts in this
			# coinc, and their time slide dictionary
			rows = [(self.sngl_burst_table.row_from_cols(row), row[-1]) for row in cursor.execute("""
	SELECT
		sngl_burst.*,
		time_slide.offset
	FROM
		sngl_burst
		JOIN coinc_event_map ON (
			coinc_event_map.table_name == 'sngl_burst'
			AND coinc_event_map.event_id == sngl_burst.event_id
		)
		JOIN coinc_event ON (
			coinc_event.coinc_event_id == coinc_event_map.coinc_event_id
		)
		JOIN time_slide ON (
			coinc_event.time_slide_id == time_slide.time_slide_id
			AND time_slide.instrument == sngl_burst.ifo
		)
	WHERE
		coinc_event.coinc_event_id == ?
			""", (coinc_event_id,))]
			# pass the events to whatever wants them
			yield sim, [event for event, offset in rows], offsetvector((event.ifo, offset) for event, offset in rows)
		cursor.close()


def summarize_coinc_database(contents, filename = None):
	if filename is None:
		filename = ""
	else:
		filename = "%s: " % filename
	cursor = contents.connection.cursor()
	print >>sys.stderr, "%sdatabase stats:" % filename
	for instrument, seglist in sorted(contents.seglists.items()):
		print >>sys.stderr, "\t%s%s livetime: %g s (%g%% vetoed)" % (filename, instrument, abs(seglist), 100.0 * float(abs(instrument in contents.vetoseglists and (seglist & contents.vetoseglists[instrument]) or 0.0)) / float(abs(seglist)))
	if contents.sngl_burst_table is not None:
		print >>sys.stderr, "\t%sburst events: %d" % (filename, len(contents.sngl_burst_table))
	if contents.sim_burst_table is not None:
		print >>sys.stderr, "\t%sburst injections: %d" % (filename, len(contents.sim_burst_table))
	if contents.time_slide_table is not None:
		print >>sys.stderr, "\t%stime slides: %d" % (filename, cursor.execute("SELECT COUNT(DISTINCT(time_slide_id)) FROM time_slide").fetchone()[0])
	if contents.coinc_def_table is not None:
		for description, n in cursor.execute("SELECT description, COUNT(*) FROM coinc_definer NATURAL JOIN coinc_event GROUP BY coinc_def_id ORDER BY description"):
			print >>sys.stderr, "\t%s%s: %d" % (filename, description, n)
	cursor.close()


def coinc_sngl_bursts(contents, coinc_event_id):
	for values in contents.connection.cursor().execute("""
SELECT sngl_burst.* FROM
	sngl_burst
	JOIN coinc_event_map ON (
		sngl_burst.event_id == coinc_event_map.event_id
		AND coinc_event_map.table_name == 'sngl_burst'
	)
WHERE
	coinc_event_map.coinc_event_id == ?
	""", (coinc_event_id,)):
		yield contents.sngl_burst_table.row_from_cols(values)


#
# =============================================================================
#
#                               Live Time Tools
#
# =============================================================================
#


def get_time_slides(connection):
	"""
	Query the database for the IDs and offsets of all time slides, and
	return two dictionaries one containing the all-zero time slides and
	the other containing the not-all-zero time slides.
	"""
	time_slides = dbtables.TimeSlideTable(connection = connection).as_dict()

	zero_lag_time_slides = dict((time_slide_id, offsetvector) for time_slide_id, offsetvector in time_slides.items() if not any(offsetvector.values()))
	background_time_slides = dict((time_slide_id, offsetvector) for time_slide_id, offsetvector in time_slides.items() if any(offsetvector.values()))

	return zero_lag_time_slides, background_time_slides


def time_slides_livetime(seglists, time_slides, verbose = False):
	"""
	Given a sequence of time slides (each of which is an instrument -->
	offset dictionary), use the segmentlistdict dictionary of segment
	lists to compute the live time in each time slide.  Return the sum
	of the live times from all slides.
	"""
	livetime = 0.0
	old_offsets = seglists.offsets.copy()
	N = len(time_slides)
	if verbose:
		print >>sys.stderr, "computing the live time for %d time slides:" % N
	for n, time_slide in enumerate(time_slides):
		if verbose:
			print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / N),
		seglists.offsets.update(time_slide)
		livetime += float(abs(seglists.intersection(time_slide.keys())))
	seglists.offsets.update(old_offsets)
	if verbose:
		print >>sys.stderr, "\t100.0%"
	return livetime


#
# =============================================================================
#
#                                 TeX Helpers
#
# =============================================================================
#


floatpattern = re.compile("([+-]?[.0-9]+)[Ee]([+-]?[0-9]+)")

def latexnumber(s):
	"""
	Convert a string of the form "d.dddde-dd" to "d.dddd \times
	10^{-dd}"

	Example:

	>>> import math
	>>> print latexnumber("%.16e" % math.pi)
	3.1415926535897931 \\times 10^{0}
	"""
	m, e = floatpattern.match(s).groups()
	return r"%s \times 10^{%d}" % (m, int(e))


#
# =============================================================================
#
#                                    Plots
#
# =============================================================================
#


def make_burst_plot(x_label, y_label, width = 165.0):
	"""
	width is in mm
	"""
	fig = figure.Figure()
	FigureCanvas(fig)
	# width mm wide, golden ratio high
	fig.set_size_inches(width / 25.4, width / 25.4 / ((1 + math.sqrt(5)) / 2))
	axes = fig.gca()
	axes.grid(True)
	axes.set_xlabel(x_label)
	axes.set_ylabel(y_label)
	return fig, axes
