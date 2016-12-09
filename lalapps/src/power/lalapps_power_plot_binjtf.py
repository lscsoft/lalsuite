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


import math
from matplotlib import cm, colors, collections
import numpy
from optparse import OptionParser
import os
import sqlite3
import sys


from glue import segments
from glue.ligolw import dbtables
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
		version = "Name: %%prog\n%s" % git_version.verbose_msg
	)
	parser.add_option("-b", "--base", metavar = "base", default = "plotbinjtf_", help = "set the prefix for output filenames (default = plotbinj_)")
	parser.add_option("-f", "--format", metavar = "format", default = "png", help = "set the output image format (default = png)")
	parser.add_option("-v", "--verbose", action = "store_true", help = "be verbose")
	options, filenames = parser.parse_args()

	return options, (filenames or [None])

options, filenames = parse_command_line()


#
# =============================================================================
#
#                             Time-Frequency Plane
#
# =============================================================================
#

def time_freq_plot(database, instrument, sim):
	fig = SnglBurstUtils.figure.Figure()
	SnglBurstUtils.FigureCanvas(fig)
	# 6.5" wide, golden ratio high
	fig.set_size_inches(6.5, 6.5 / ((1 + math.sqrt(5)) / 2))

	#
	# top plot --- triggers that matched injection
	#

	t_sim = sim.time_at_instrument(instrument)

	axes = fig.add_subplot(211)
	axes.grid(True)
	#axes.set_xlabel("$t - t_{\mathrm{injection}}$ (s)")
	axes.set_ylabel("$f - f_{\mathrm{injection}}$ (Hz)")
	axes.set_title("%s Triggers Matching %g Hz Injection at GPS %s" % (instrument, sim.frequency, t_sim))

	xmin = xmax = 0.0
	ymin = ymax = 0.0
	verts = []
	colours = []
	peakx = []
	peaky = []
	match_ids = []
	# find triggers from the desired instrument that are coincident
	# with the injection, and iterate over them in order from least to
	# most confident
	for burst in map(database.sngl_burst_table.row_from_cols, database.connection.cursor().execute("""
SELECT sngl_burst.* FROM
	sngl_burst
	JOIN coinc_event_map AS a ON (
		sngl_burst.event_id == a.event_id
		AND a.table_name == 'sngl_burst'
	)
	JOIN coinc_event_map AS b ON (
		a.coinc_event_id == b.coinc_event_id
		AND b.table_name == 'sim_burst'
	)
WHERE
	sngl_burst.ifo == ?
	AND b.event_id == ?
ORDER BY
	sngl_burst.confidence ASC
	""", (instrument, sim.simulation_id))):
		match_ids.append(burst.event_id)

		# Add time-frequency tile to collection
		tmin = float(burst.start - t_sim)
		tmax = float(burst.start + burst.duration - t_sim)
		fmin = burst.central_freq - burst.bandwidth / 2 - sim.frequency
		fmax = burst.central_freq + burst.bandwidth / 2 - sim.frequency
		verts.append(((tmin, fmin), (tmax, fmin), (tmax, fmax), (tmin, fmax)))
		colours.append(burst.confidence)

		try:
			# draw most significant tile if there is one
			tmin = float(burst.ms_start - t_sim)
			tmax = float(burst.ms_start + burst.ms_duration - t_sim)
			fmin = burst.ms_flow - sim.frequency
			fmax = burst.ms_flow + burst.ms_bandwidth - sim.frequency
			verts.append(((tmin, fmin), (tmax, fmin), (tmax, fmax), (tmin, fmax)))
			colours.append(burst.ms_confidence)
		except AttributeError:
			pass

		peakx.append(float(burst.peak - t_sim))
		try:
			# use peak_frequency col if it exists
			peaky.append(burst.peak_frequency - sim.frequency)
		except AttributeError:
			peaky.append(burst.central_freq - sim.frequency)

		# update bounding box
		tmin = float(burst.start - t_sim)
		tmax = float(burst.start + burst.duration - t_sim)
		fmin = burst.central_freq - burst.bandwidth / 2 - sim.frequency
		fmax = burst.central_freq + burst.bandwidth / 2 - sim.frequency
		xmin = min(xmin, tmin)
		xmax = max(xmax, tmax)
		ymin = min(ymin, fmin)
		ymax = max(ymax, fmax)

	polys = collections.PolyCollection(verts)
	polys.set_array(numpy.array(colours))
	polys.set_alpha(0.3)
	polys.set_cmap(cm.get_cmap())
	polys.set_norm(colors.normalize())
	axes.add_collection(polys)

	axes.plot(peakx, peaky, "k+")

	axes.axvline(0, color = "k")
	axes.axhline(0, color = "k")

	# set the bounding box
	axes.set_xlim([1.4 * xmin, 1.4 * xmax])
	axes.set_ylim([1.4 * ymin, 1.4 * ymax])

	#
	# bottom plot --- triggers near injection
	#

	axes = fig.add_subplot(212)
	axes.grid(True)
	axes.set_xlabel("$t - t_{\mathrm{injection}}$ (s)")
	axes.set_ylabel("$f - f_{\mathrm{injection}}$ (Hz)")

	xmin = xmax = 0.0
	ymin = ymax = 0.0
	verts = []
	colours = []
	edgecolours = []
	peakx = []
	peaky = []
	# find triggers from the desired instrument that are near the
	# injection, and iterate over them in order from least to most
	# confident
	for burst in map(database.sngl_burst_table.row_from_cols, database.connection.cursor().execute("""
SELECT * FROM
	sngl_burst
WHERE
	ifo == ?
	AND start_time BETWEEN ? AND ?
	AND central_freq BETWEEN ? AND ?
ORDER BY
	sngl_burst.confidence ASC
	""", (instrument, int(t_sim - 2), int(t_sim + 2), sim.frequency - 300, sim.frequency + 300))):
		# Add time-frequency tile to collection
		tmin = float(burst.start - t_sim)
		tmax = float(burst.start + burst.duration - t_sim)
		fmin = burst.central_freq - burst.bandwidth / 2 - sim.frequency
		fmax = burst.central_freq + burst.bandwidth / 2 - sim.frequency
		verts.append(((tmin, fmin), (tmax, fmin), (tmax, fmax), (tmin, fmax)))
		colours.append(burst.confidence)
		if burst.event_id in match_ids:
			edgecolours.append("g")
		else:
			edgecolours.append("k")

		peakx.append(float(burst.peak - t_sim))
		try:
			# use peak_frequency col if it exists
			peaky.append(burst.peak_frequency - sim.frequency)
		except:
			peaky.append(burst.central_freq - sim.frequency)

		# update bounding box
		xmin = min(xmin, tmin)
		xmax = max(xmax, tmax)
		ymin = min(ymin, fmin)
		ymax = max(ymax, fmax)

	polys = collections.PolyCollection(verts, edgecolors = edgecolours)
	polys.set_array(numpy.array(colours))
	polys.set_alpha(0.3)
	polys.set_cmap(cm.get_cmap())
	polys.set_norm(colors.normalize())
	axes.add_collection(polys)

	axes.plot(peakx, peaky, "k+")

	axes.axvline(0, color = "k")
	axes.axhline(0, color = "k")

	# set the bounding box
	axes.set_xlim([1.4 * xmin, 1.4 * xmax])
	axes.set_ylim([1.4 * ymin, 1.4 * ymax])


	return fig


#
# =============================================================================
#
#                                     Plot
#
# =============================================================================
#


def found_injections(contents, instrument):
	for values in contents.connection.cursor().execute("""
SELECT
	*
FROM
	sim_burst
WHERE
	EXISTS (
		-- Find a link through the coinc_event_map table to a row
		-- in the sngl_burst table with the correct ifo value.
		SELECT
			*
		FROM
			coinc_event_map AS a
			JOIN coinc_event_map AS b ON (
				a.coinc_event_id == b.coinc_event_id
			)
			JOIN sngl_burst ON (
				b.table_name == 'sngl_burst'
				AND b.event_id == sngl_burst.event_id
			)
		WHERE
			a.table_name == 'sim_burst'
			AND a.event_id == sim_burst.simulation_id
			AND sngl_burst.ifo == ?
	)
	""", (instrument,)):
		yield contents.sim_burst_table.row_from_cols(values)


for n, filename in enumerate(filenames):
	if options.verbose:
		print >>sys.stderr, "%d/%d: %s" % (n + 1, len(filenames), filename)
	database = SnglBurstUtils.CoincDatabase(sqlite3.connect(filename), "lalapps_power")
	if options.verbose:
		SnglBurstUtils.summarize_coinc_database(database)
	for instrument in database.instruments:
		for sim in found_injections(database, instrument):
			plotname = "%s%d_%s.%s" % (options.base, sim.time_at_instrument(instrument).seconds, instrument, options.format)
			if options.verbose:
				print >>sys.stderr, "--> %s" % plotname
			time_freq_plot(database, instrument, sim).savefig(plotname)
	database.connection.close()
