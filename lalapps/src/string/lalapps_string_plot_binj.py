#
# Copyright (C) 2009  Kipp Cannon
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


from lal.utils import CacheEntry


from glue import segments
from glue.ligolw import dbtables
from glue.ligolw import utils as ligolw_utils
from lalburst import git_version
from pylal import rate
from lalburst import SimBurstUtils
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
	parser.add_option("-b", "--base", metavar = "base", default = "string_plotbinj_", help = "Set the prefix for output filenames (default = \"string_plotbinj_\").")
	parser.add_option("-f", "--format", metavar = "format", action = "append", default = [], help = "Set the output image format (default = \"png\").  Option can be given multiple times to generate plots in multiple formats.")
	parser.add_option("--amplitude-func", metavar = "det|wave", default = "det", help = "Select the amplitude to show on the plots.  \"det\" = the amplitude expected in the detector, \"wave\" = the amplitude of the wave (default = \"det\").")
	parser.add_option("--input-cache", metavar = "filename", action = "append", default = [], help = "Get list of files from this LAL cache.")
	parser.add_option("-l", "--live-time-program", metavar = "program", default = "StringSearch", help = "Set the name of the program, as it appears in the process table, whose search summary entries define the search live time (default = \"StringSearch\").")
	parser.add_option("--plot", metavar = "number", action = "append", default = None, help = "Generate the given plot number (default = make all plots).")
	parser.add_option("-t", "--tmp-space", metavar = "path", help = "Path to a directory suitable for use as a work area while manipulating the database file.  The database file will be worked on in this directory, and then moved to the final location when complete.  This option is intended to improve performance when running in a networked environment, where there might be a local disk with higher bandwidth than is available to the filesystem on which the final output will reside.")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()

	if options.amplitude_func == "wave":
		options.amplitude_func = lambda sim, instrument, offsetvector: sim.amplitude
		options.amplitude_lbl = r"$A^{\mathrm{wave}}$"
	elif options.amplitude_func == "det":
		options.amplitude_func = SimBurstUtils.string_amplitude_in_instrument
		options.amplitude_lbl = r"$A^{\mathrm{det}}$"
	else:
		raise ValueError("unrecognized --amplitude-func %s" % options.amplitude_func)

	if not options.format:
		options.format = ["png"]

	if options.plot is not None:
		# use sorted(set(...)) to ensure the indexes are ordered
		# and unique
		options.plot = sorted(set(map(int, options.plot)))

	filenames = filenames or []
	for cache in options.input_cache:
		if options.verbose:
			print >>sys.stderr, "reading '%s' ..." % cache
		filenames += [CacheEntry(line).path for line in file(cache)]
	if not filenames:
		raise ValueError("Nothing to do!")

	return options, filenames


#
# =============================================================================
#
#                           Trigger Count Histogram
#
# =============================================================================
#


class TriggerCountHistogram(object):
	def __init__(self):
		self.found = 0
		self.bins = []

	def add_contents(self, contents):
		for nevents, in contents.connection.cursor().execute("""
SELECT
	nevents
FROM
	coinc_event
WHERE
	coinc_def_id == ?
		""", (contents.sb_definer_id,)):
			self.found += 1
			while nevents + 1 >= len(self.bins):
				self.bins.append(0)
			self.bins[nevents] += 1

	def finish(self):
		fig, axes = SnglBurstUtils.make_burst_plot("Number of Triggers Coincident with Injection", "Count")
		axes.semilogy()
		axes.plot(range(len(self.bins)), self.bins, "ko-")
		axes.set_title("Triggers per Found Injection\n(%d Found Injections)" % self.found)

		return fig,


#
# =============================================================================
#
#                         Recovered vs. Injected h_rss
#
# =============================================================================
#


class RecoveredVsInjectedAmplitude(object):
	class Data(object):
		def __init__(self):
			self.found = 0
			self.x = []
			self.y = []
			self.c = []

	def __init__(self, amplitude_func, amplitude_lbl):
		self.amplitude_func = amplitude_func
		self.amplitude_lbl = amplitude_lbl
		self.data = {}

	def add_contents(self, contents):
		offsetvectors = contents.time_slide_table.as_dict()
		for values in contents.connection.cursor().execute("""
SELECT
	sim_burst.*,
	sngl_burst.ifo,
	sngl_burst.amplitude
FROM
	sim_burst
	JOIN sim_burst_map AS map ON (
		map.simulation_id == sim_burst.simulation_id
	)
	JOIN sngl_burst ON (
		sngl_burst.event_id == map.event_id
	)
		"""):
			sim = contents.sim_burst_table.row_from_cols(values)
			instrument, amplitude_rec = values[-2:]
			try:
				data = self.data[instrument]
			except KeyError:
				data = self.data[instrument] = self.Data()
			data.found += 1
			data.x.append(abs(self.amplitude_func(sim, instrument, offsetvectors[sim.time_slide_id])))
			data.y.append(abs(amplitude_rec))
			data.c.append(math.log(sim.frequency))

	def finish(self):
		for instrument, data in sorted(self.data.items()):
			fig, axes = SnglBurstUtils.make_burst_plot("Injected %s" % self.amplitude_lbl, r"Recovered $A$")
			axes.loglog()
			fig.set_size_inches(8, 8)
			axes.scatter(data.x, data.y, c = data.c, s = 16)
			#xmin, xmax = axes.get_xlim()
			#ymin, ymax = axes.get_ylim()
			xmin, xmax = min(data.x), max(data.x)
			ymin, ymax = min(data.y), max(data.y)
			xmin = ymin = min(xmin, ymin)
			xmax = ymax = max(xmax, ymax)
			axes.plot([xmin, xmax], [ymin, ymax], "k-")
			axes.set_xlim([xmin, xmax])
			axes.set_ylim([ymin, ymax])
			axes.set_title(r"Recovered $A$ vs.\ Injected %s in %s (%d Recovered Injections)" % (self.amplitude_lbl, instrument, data.found))
			yield fig


#
# =============================================================================
#
#                            Recovered Time Offset
#
# =============================================================================
#


class RecoveredTimeOffset(object):
	class Data(object):
		def __init__(self, binning):
			self.found = 0
			self.offsets = rate.BinnedDensity(binning)

	def __init__(self, interval, width):
		# 21 bins per filter width
		bins = int(float(abs(interval)) / width) * 21
		self.binning = rate.NDBins((rate.LinearBins(interval[0], interval[1], bins),))
		self.data = {}

	def add_contents(self, contents):
		offsetvectors = contents.time_slide_table.as_dict()
		for values in contents.connection.cursor().execute("""
SELECT
	sim_burst.*,
	sngl_burst.ifo,
	sngl_burst.snr,
	sngl_burst.peak_time, sngl_burst.peak_time_ns
FROM
	sim_burst
	JOIN sim_burst_map AS map ON (
		map.simulation_id == sim_burst.simulation_id
	)
	JOIN sngl_burst ON (
		sngl_burst.event_id == map.event_id
	)
		"""):
			sim = contents.sim_burst_table.row_from_cols(values)
			instrument = values[-4]
			burst_snr = values[-3]
			burst_peak = dbtables.lsctables.LIGOTimeGPS(*values[-2:])
			sim_peak = sim.time_at_instrument(instrument, offsetvectors[sim.time_slide_id])

			try:
				data = self.data[instrument]
			except KeyError:
				data = self.data[instrument] = self.Data(self.binning)

			data.found += 1
			try:
				data.offsets.count[float(burst_peak - sim_peak),] += 1.0
			except IndexError:
				# outside plot range
				pass

	def finish(self):
		for instrument, data in sorted(self.data.items()):
			fig, axes = SnglBurstUtils.make_burst_plot(r"$t_{\mathrm{recovered}} - t_{\mathrm{injected}}$ (s)", "Triggers per Unit Offset")
			axes.semilogy()
			axes.set_title("Trigger Peak Time - Injection Peak Time in %s\n(%d Found Injections)" % (instrument, data.found))
			# 21 bins per filter width
			rate.filter_array(data.offsets.array, rate.gaussian_window(21))
			axes.plot(data.offsets.centres()[0], data.offsets.at_centres(), "k")
			#axes.legend(["%s residuals" % instrument, "SNR-weighted mean of residuals in all instruments"], loc = "lower right")
			yield fig


#
# =============================================================================
#
#                    Recovered vs. Injected Frequency Plots
#
# =============================================================================
#


class RecoveredFreq(object):
	class DataA(object):
		def __init__(self):
			self.found = 0
			self.x = []
			self.y = []

	class DataB(object):
		def __init__(self):
			self.found = 0
			self.x = []
			self.y = []
			self.c = []

	def __init__(self):
		self.dataA = {}
		self.dataB = {}

	def add_contents(self, contents):
		simulation_ids = {}
		for simulation_id, sim_frequency, instrument, burst_frequency, burst_snr in contents.connection.cursor().execute("""
SELECT
	sim_burst.simulation_id,
	sim_burst.frequency,
	sngl_burst.ifo,
	-- central_freq = (f_cut + f_bankstart) / 2
	-- bandwidth = f_cut - f_bankstart
	-- therefore f_cut = ...
	sngl_burst.central_freq + sngl_burst.bandwidth / 2,
	sngl_burst.snr
FROM
	sim_burst
	JOIN sim_burst_map AS map ON (
		map.simulation_id == sim_burst.simulation_id
	)
	JOIN sngl_burst ON (
		sngl_burst.event_id == map.event_id
	)
		"""):
			simulation_ids.setdefault(instrument, set()).add(simulation_id)
			try:
				dataA = self.dataA[instrument]
				dataB = self.dataB[instrument]
			except KeyError:
				dataA = self.dataA[instrument] = self.DataA()
				dataB = self.dataB[instrument] = self.DataB()

			dataA.x.append(burst_snr)
			dataA.y.append((burst_frequency - sim_frequency) / sim_frequency)

			dataB.x.append(sim_frequency)
			dataB.y.append(burst_frequency)
			dataB.c.append(math.log(burst_snr))

		# count number of distinct injections that had matching
		# burst triggers
		for instrument, simulation_ids in simulation_ids.items():
			self.dataA[instrument].found += len(simulation_ids)
			self.dataB[instrument].found += len(simulation_ids)

	def finish(self):
		for instrument, data in sorted(self.dataA.items()):
			fig, axes = SnglBurstUtils.make_burst_plot(r"SNR in %s" % instrument, r"$(f_{\mathrm{recovered}} - f_{\mathrm{injected}})/ f_{\mathrm{injected}}$")

			axes.set_title("Cut-off Frequency Residual vs.\ SNR in %s\n(%d Found Injections)" % (instrument, data.found))
			axes.loglog()
			axes.plot(data.x, data.y, "k+")

			yield fig

		for instrument, data in sorted(self.dataB.items()):
			fig, axes = SnglBurstUtils.make_burst_plot(r"$f_{\mathrm{injected}}$ (Hz)", r"$f_{\mathrm{recovered}}$ (Hz)")
			axes.loglog()
			fig.set_size_inches(8, 8)
			axes.scatter(data.x, data.y, c = data.c, s = 16)
			xmin, xmax = min(data.x), max(data.x)
			ymin, ymax = min(data.y), max(data.y)
			xmin = ymin = min(xmin, ymin)
			xmax = ymax = max(xmax, ymax)
			axes.plot([xmin, xmax], [ymin, ymax], "k-")
			axes.set_xlim([xmin, xmax])
			axes.set_ylim([ymin, ymax])
			axes.set_title(r"Recovered Cut-off Frequency vs.\ Injected Cut-off Frequency in %s\\(%d Injections; red = high recovered SNR, blue = low recovered SNR)" % (instrument, data.found))

			yield fig


#
# =============================================================================
#
#                        \chi^{2} vs. \rho Scatter Plot
#
# =============================================================================
#


class Chi2VsRho(object):
	class Data(object):
		def __init__(self):
			self.injections_x = []
			self.injections_y = []
			self.noninjections_x = []
			self.noninjections_y = []

	def __init__(self):
		self.data = {}

	def add_contents(self, contents):
		for instrument, chi2, chi2dof, snr, is_injection in contents.connection.cursor().execute("""
SELECT
	sngl_burst.ifo,
	sngl_burst.chisq,
	sngl_burst.chisq_dof,
	sngl_burst.snr,
	sngl_burst.event_id IN (
		SELECT
			coinc_event_map.event_id
		FROM
			coinc_event_map
			JOIN coinc_event ON (
				coinc_event.coinc_event_id == coinc_event_map.coinc_event_id
			)
		WHERE
			coinc_event_map.table_name == 'sngl_burst'
			AND coinc_event.coinc_def_id == ?
	)
FROM
	sngl_burst
		""", (contents.sb_definer_id,)):
			if is_injection:
				try:
					self.data[instrument].injections_x.append(snr)
					self.data[instrument].injections_y.append(chi2 / chi2dof)
				except KeyError:
					self.data[instrument] = self.Data()
					self.data[instrument].injections_x.append(snr)
					self.data[instrument].injections_y.append(chi2 / chi2dof)
			else:
				try:
					self.data[instrument].noninjections_x.append(snr)
					self.data[instrument].noninjections_y.append(chi2 / chi2dof)
				except KeyError:
					self.data[instrument] = self.Data()
					self.data[instrument].noninjections_x.append(snr)
					self.data[instrument].noninjections_y.append(chi2 / chi2dof)

	def finish(self):
		for instrument, data in sorted(self.data.items()):
			fig, axes = SnglBurstUtils.make_burst_plot(r"SNR", r"$\chi^{2} / \mathrm{DOF}$")

			axes.set_title("$\chi^{2} / \mathrm{DOF}$ vs.\ SNR in %s" % instrument)
			axes.loglog()
			axes.plot(data.injections_x, data.injections_y, "r+")
			axes.plot(data.noninjections_x, data.noninjections_y, "k+")

			yield fig


#
# =============================================================================
#
#                                     Main
#
# =============================================================================
#


#
# Parse command line
#


options, filenames = parse_command_line()


#
# Initialize plots
#


plots = [
	(0, True, TriggerCountHistogram()),
	(1, True, RecoveredVsInjectedAmplitude(options.amplitude_func, options.amplitude_lbl)),
	(2, True, RecoveredTimeOffset(segments.segment(-0.01, +0.01), 0.00005)),
	(3, True, RecoveredFreq()),
	(4, False, Chi2VsRho())
]

if options.plot is not None:
	plots = [plots[i] for i in options.plot]


#
# Process files
#


for n, filename in enumerate(ligolw_utils.sort_files_by_size(filenames, options.verbose, reverse = True)):
	if options.verbose:
		print >>sys.stderr, "%d/%d: %s" % (n + 1, len(filenames), filename)
	working_filename = dbtables.get_connection_filename(filename, tmp_path = options.tmp_space, verbose = options.verbose)
	database = SnglBurstUtils.CoincDatabase(sqlite3.connect(working_filename), options.live_time_program, search = "StringCusp")
	if options.verbose:
		SnglBurstUtils.summarize_coinc_database(database)
	is_injection_db = "sim_burst" in dbtables.get_table_names(database.connection)
	if is_injection_db:
		database.connection.cursor().execute("""
CREATE TEMPORARY TABLE
	sim_burst_map
AS
	SELECT
		a.event_id AS simulation_id,
		a.coinc_event_id AS coinc_event_id,
		b.event_id AS event_id
	FROM
		coinc_event_map AS a
		JOIN coinc_event ON (
			coinc_event.coinc_event_id == a.coinc_event_id
		)
		JOIN coinc_event_map AS b ON (
			b.table_name == 'sngl_burst'
			AND b.coinc_event_id == a.coinc_event_id
		)
	WHERE
		a.table_name == 'sim_burst'
		AND coinc_event.coinc_def_id == ?
		""", (database.sb_definer_id,))
	for n, requires_injection_db, plot in plots:
		if requires_injection_db and not is_injection_db:
			continue
		if options.verbose:
			print >>sys.stderr, "adding to plot group %d ..." % n
		plot.add_contents(database)
	database.connection.close()
	dbtables.discard_connection_filename(filename, working_filename, verbose = options.verbose)


#
# finish and write plots, deleting them as we go to save memory
#


format = "%%s%%0%dd%%s.%%s" % (int(math.log10(len(plots) or 1)) + 1)
while len(plots):
	n, requires_injection_db, plot = plots.pop(0)
	if options.verbose:
		print >>sys.stderr, "finishing plot group %d ..." % n
	try:
		figs = plot.finish()
	except ValueError as e:
		print >>sys.stderr, "can't finish plot group %d: %s" % (n, str(e))
	else:
		for f, fig in enumerate(figs):
			for extension in options.format:
				filename = format % (options.base, n, chr(97 + f), extension)
				if options.verbose:
					print >>sys.stderr, "writing %s ..." % filename
				fig.savefig(filename)
