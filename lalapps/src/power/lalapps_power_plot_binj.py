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


import math
import numpy
from optparse import OptionParser
import sqlite3
import sys


from glue import segments
from glue.ligolw import dbtables
from glue.ligolw import utils
from lal.utils import CacheEntry
from lalburst import git_version
from lalburst import SimBurstUtils
from lalburst import SnglBurstUtils
from pylal import rate


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
	parser.add_option("--made-only", action = "store_true", default = False, help = "Plot only injections that were made.")
	parser.add_option("-b", "--base", metavar = "base", default = "plotbinj_", help = "Set the prefix for output filenames (default = \"plotbinj_\")")
	parser.add_option("-f", "--format", metavar = "format", default = "png", help = "Set the output image format (default = \"png\")")
	parser.add_option("--amplitude-func", metavar = "hrsswave|hrssdet|E", default = "hrsswave", help = "Select the amplitude to show on the plots.  \"hrsswave\" = the h_rss of the wave, \"hrssdet\" = the h_rss in the detector, \"E\" = the radiated energy over r^2.")
	parser.add_option("--input-cache", metavar = "filename", action = "append", default = [], help = "Get list of trigger files from this LAL cache file.")
	parser.add_option("-l", "--live-time-program", metavar = "program", default = "lalapps_power", help = "Set the name, as it appears in the process table, of the program whose search summary entries define the search live time (default = \"lalapps_power\").")
	parser.add_option("--plot", metavar = "number", action = "append", default = None, help = "Generate the given plot number (default = make all plots).  Use \"none\" to disable plots.")
	parser.add_option("--coinc-plot", metavar = "number", action = "append", default = None, help = "Generate the given coinc plot number (default = make all coinc plots).  Use \"none\" to disable coinc plots.")
	parser.add_option("-t", "--tmp-space", metavar = "path", help = "Path to a directory suitable for use as a work area while manipulating the database file.  The database file will be worked on in this directory, and then moved to the final location when complete.  This option is intended to improve performance when running in a networked environment, where there might be a local disk with higher bandwidth than is available to the filesystem on which the final output will reside.")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()

	if options.amplitude_func == "hrsswave":
		options.amplitude_func = lambda sim, instrument, offsetvector: sim.hrss
		options.amplitude_lbl = r"$h_{\mathrm{rss}}$"
	elif options.amplitude_func == "hrssdet":
		options.amplitude_func = SimBurstUtils.hrss_in_instrument
		options.amplitude_lbl = r"$h_{\mathrm{rss}}^{\mathrm{det}}$"
	elif options.amplitude_func == "E":
		options.amplitude_func = lambda sim, instrument, offsetvector: sim.egw_over_rsquared
		options.amplitude_lbl = r"$M_{\odot} / \mathrm{pc}^{2}$"
	else:
		raise ValueError("unrecognized --amplitude-func %s" % options.amplitude_func)

	if options.plot is None:
		options.plot = range(10)
	elif "none" in options.plot:
		options.plot = []
	else:
		options.plot = map(int, options.plot)
	if options.coinc_plot is None:
		options.coinc_plot = range(2)
	elif "none" in options.coinc_plot:
		options.coinc_plot = []
	else:
		options.coinc_plot = map(int, options.coinc_plot)

	filenames = filenames or []
	for cache in options.input_cache:
		if options.verbose:
			print >>sys.stderr, "reading '%s' ..." % cache
		filenames += [CacheEntry(line).path for line in file(cache)]

	return options, filenames


#
# =============================================================================
#
#                              Frequency vs. Time
#
# =============================================================================
#


class FreqVsTime(object):
	def __init__(self, instrument):
		self.fig, self.axes = SnglBurstUtils.make_burst_plot("GPS Time (s)", "Frequency (Hz)")
		self.axes.semilogy()
		self.instrument = instrument
		self.num_injections = 0
		self.injected_x = []
		self.injected_y = []
		self.missed_x = []
		self.missed_y = []
		self.seglist = segments.segmentlist()

	def add_contents(self, contents):
		self.seglist |= contents.seglists[self.instrument]
		for values in contents.connection.cursor().execute("""
SELECT
	*,
	EXISTS (
		SELECT
			*
		FROM
			sim_burst_map
			JOIN sngl_burst ON (
				sngl_burst.event_id == sim_burst_map.event_id
			)
		WHERE
			sim_burst_map.simulation_id == sim_burst.simulation_id
			AND sngl_burst.ifo == ?
	)
FROM
	sim_burst
		""", (self.instrument,)):
			sim = contents.sim_burst_table.row_from_cols(values)
			found = values[-1]
			if SimBurstUtils.injection_was_made(sim, contents.seglists, (self.instrument,)):
				self.num_injections += 1
				peak = float(sim.time_at_instrument(self.instrument))
				self.injected_x.append(peak)
				self.injected_y.append(sim.frequency)
				if not found:
					self.missed_x.append(peak)
					self.missed_y.append(sim.frequency)

	def finish(self):
		self.axes.plot(self.injected_x, self.injected_y, "k+")
		if not options.made_only:
			self.axes.plot(self.missed_x, self.missed_y, "rx")
		for seg in ~self.seglist & segments.segmentlist([segments.segment(self.axes.get_xlim())]):
			self.axes.axvspan(float(seg[0]), float(seg[1]), facecolor = "k", alpha = 0.2)
		self.axes.set_ylim([min(self.injected_y), max(self.injected_y)])
		self.axes.set_title("Injection Locations\n(%d Injections)" % self.num_injections)


#
# =============================================================================
#
#                           Amplitude vs. Frequency
#
# =============================================================================
#


class HrssVsFreqScatter(object):
	def __init__(self, instrument, amplitude_func, amplitude_lbl):
		self.fig, self.axes = SnglBurstUtils.make_burst_plot("Frequency (Hz)", amplitude_lbl)
		self.axes.loglog()
		self.instrument = instrument
		self.amplitude_func = amplitude_func
		self.amplitude_lbl = amplitude_lbl
		self.num_injections = 0
		self.injected_x = []
		self.injected_y = []
		self.missed_x = []
		self.missed_y = []

	def add_contents(self, contents):
		offsetvectors = contents.time_slide_table.as_dict()
		for values in contents.connection.cursor().execute("""
SELECT
	*,
	EXISTS (
		SELECT
			*
		FROM
			sim_burst_map
			JOIN sngl_burst ON (
				sngl_burst.event_id == sim_burst_map.event_id
			)
		WHERE
			sim_burst_map.simulation_id == sim_burst.simulation_id
			AND sngl_burst.ifo == ?
	)
FROM
	sim_burst
		""", (self.instrument,)):
			sim = contents.sim_burst_table.row_from_cols(values)
			found = values[-1]
			if SimBurstUtils.injection_was_made(sim, contents.seglists, (self.instrument,)):
				self.num_injections += 1
				amplitude = self.amplitude_func(sim, self.instrument, offsetvectors[sim.time_slide_id])
				self.injected_x.append(sim.frequency)
				self.injected_y.append(amplitude)
				if not found:
					self.missed_x.append(sim.frequency)
					self.missed_y.append(amplitude)

	def finish(self):
		self.axes.plot(self.injected_x, self.injected_y, "k+")
		if not options.made_only:
			self.axes.plot(self.missed_x, self.missed_y, "rx")
		self.axes.set_xlim([min(self.injected_x), max(self.injected_x)])
		self.axes.set_ylim([min(self.injected_y), max(self.injected_y)])
		self.axes.set_title("Injection " + self.amplitude_lbl + " vs. Frequency\n(%d Injections)" % self.num_injections)


#
# =============================================================================
#
#                           Trigger Count Histogram
#
# =============================================================================
#


class TriggerCountHistogram(object):
	def __init__(self, instrument):
		self.fig, self.axes = SnglBurstUtils.make_burst_plot("Number of Triggers Coincident with Injection", "Count")
		self.axes.semilogy()
		self.instrument = instrument
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
		self.axes.plot(range(len(self.bins)), self.bins, "ko-")
		self.axes.set_title("Triggers per Found Injection\n(%d Found Injections)" % self.found)


#
# =============================================================================
#
#                         Recovered vs. Injected h_rss
#
# =============================================================================
#


class RecoveredVsInjectedhrss(object):
	def __init__(self, instrument, amplitude_func, amplitude_lbl):
		self.fig, self.axes = SnglBurstUtils.make_burst_plot(r"Injected " + amplitude_lbl, r"Recovered $h_{\mathrm{rss}}$")
		self.axes.loglog()
		self.fig.set_size_inches(8, 8)
		self.instrument = instrument
		self.amplitude_func = amplitude_func
		self.amplitude_lbl = amplitude_lbl
		self.matches = 0
		self.x = []
		self.y = []
		self.c = []

	def add_contents(self, contents):
		offsetvectors = contents.time_slide_table.as_dict()
		for values in contents.connection.cursor().execute("""
SELECT
	sim_burst.*,
	sngl_burst.peak_frequency,
	sngl_burst.ms_hrss
FROM
	sim_burst
	JOIN sim_burst_map ON (
		sim_burst_map.simulation_id == sim_burst.simulation_id
	)
	JOIN sngl_burst ON (
		sngl_burst.event_id == sim_burst_map.event_id
	)
WHERE
	sngl_burst.ifo == ?
		""", (self.instrument,)):
			sim = contents.sim_burst_table.row_from_cols(values)
			freq_rec, hrss_rec = values[-2:]
			self.matches += 1
			self.x.append(self.amplitude_func(sim, self.instrument, offsetvectors[sim.time_slide_id]))
			self.y.append(hrss_rec)
			self.c.append(math.log(freq_rec))

	def finish(self):
		self.axes.scatter(self.x, self.y, c = self.c, s = 16)
		#xmin, xmax = self.axes.get_xlim()
		#ymin, ymax = self.axes.get_ylim()
		xmin, xmax = min(self.x), max(self.x)
		ymin, ymax = min(self.y), max(self.y)
		xmin = ymin = min(xmin, ymin)
		xmax = ymax = max(xmax, ymax)
		self.axes.plot([xmin, xmax], [ymin, ymax], "k-")
		self.axes.set_xlim([xmin, xmax])
		self.axes.set_ylim([ymin, ymax])
		self.axes.set_title(r"Recovered $h_{\mathrm{rss}}$ vs.\ Injected " + self.amplitude_lbl + " (%d Events Matching Injections)" % self.matches)


class RecoveredPerInjectedhrssVsFreq(object):
	def __init__(self, instrument, amplitude_func, amplitude_lbl):
		self.fig, self.axes = SnglBurstUtils.make_burst_plot(r"$f_{\mathrm{injected}}$ (Hz)", r"Recovered $h_{\mathrm{rss}}$ / Injected "+ amplitude_lbl)
		self.axes.loglog()
		self.fig.set_size_inches(8, 8)
		self.instrument = instrument
		self.amplitude_func = amplitude_func
		self.amplitude_lbl = amplitude_lbl
		self.matches = 0
		self.x = []
		self.y = []
		self.c = []

	def add_contents(self, contents):
		offsetvectors = contents.time_slide_table.as_dict()
		for values in contents.connection.cursor().execute("""
SELECT
	sim_burst.*,
	sngl_burst.peak_frequency,
	sngl_burst.ms_hrss
FROM
	sim_burst
	JOIN sim_burst_map ON (
		sim_burst_map.simulation_id == sim_burst.simulation_id
	)
	JOIN sngl_burst ON (
		sngl_burst.event_id == sim_burst_map.event_id
	)
WHERE
	sngl_burst.ifo == ?
		""", (self.instrument,)):
			sim = contents.sim_burst_table.row_from_cols(values)
			freq_rec, hrss_rec = values[-2:]
			self.matches += 1
			self.x.append(sim.frequency)
			self.y.append(hrss_rec / self.amplitude_func(sim, self.instrument, offsetvectors[sim.time_slide_id]))
			self.c.append(math.log(freq_rec))

	def finish(self):
		self.axes.scatter(self.x, self.y, c = self.c, s = 16)
		self.axes.set_xlim([min(self.x), max(self.x)])
		self.axes.set_ylim([min(self.y), max(self.y)])
		self.axes.set_title(r"Ratio of Recovered $h_{\mathrm{rss}}$ to Injected " + self.amplitude_lbl + r" vs.\ Frequency (%d Events Matching Injections)" % self.matches)


class RecoveredPerInjectedhrssVsBandwidth(object):
	def __init__(self, instrument, amplitude_func, amplitude_lbl):
		self.fig, self.axes = SnglBurstUtils.make_burst_plot(r"$\Delta f_{\mathrm{recovered}}$ (Hz)", r"Recovered $h_{\mathrm{rss}}$ / Injected " + amplitude_lbl)
		self.axes.loglog()
		self.fig.set_size_inches(8, 8)
		self.instrument = instrument
		self.amplitude_func = amplitude_func
		self.amplitude_lbl = amplitude_lbl
		self.matches = 0
		self.x = []
		self.y = []
		self.c = []

	def add_contents(self, contents):
		offsetvectors = contents.time_slide_table.as_dict()
		for values in contents.connection.cursor().execute("""
SELECT
	sim_burst.*,
	sngl_burst.bandwidth,
	sngl_burst.peak_frequency,
	sngl_burst.ms_hrss
FROM
	sim_burst
	JOIN sim_burst_map ON (
		sim_burst_map.simulation_id == sim_burst.simulation_id
	)
	JOIN sngl_burst ON (
		sngl_burst.event_id == sim_burst_map.event_id
	)
WHERE
	sngl_burst.ifo == ?
		""", (self.instrument,)):
			sim = contents.sim_burst_table.row_from_cols(values)
			bandwidth, freq_rec, hrss_rec = values[-3:]
			self.matches += 1
			self.x.append(bandwidth)
			self.y.append(hrss_rec / self.amplitude_func(sim, self.instrument, offsetvectors[sim.time_slide_id]))
			self.c.append(math.log(freq_rec))

	def finish(self):
		self.axes.scatter(self.x, self.y, c = self.c, s = 16)
		self.axes.set_xlim([min(self.x), max(self.x)])
		self.axes.set_ylim([min(self.y), max(self.y)])
		self.axes.set_title(r"Ratio of Recovered $h_{\mathrm{rss}}$ to Injected " + self.amplitude_lbl + r" vs.\ Recovered Bandwidth (%d Events Matching Injections)" % self.matches)


#
# =============================================================================
#
#                            Recovered Time Offset
#
# =============================================================================
#


class RecoveredTimeOffset(object):
	def __init__(self, instrument, interval, width):
		self.fig, self.axes = SnglBurstUtils.make_burst_plot(r"$t_{\mathrm{recovered}} - t_{\mathrm{injected}}$ (s)", "Triggers per Unit Offset")
		self.axes.semilogy()
		self.instrument = instrument
		self.found = 0
		# 21 bins per filter width
		bins = int(float(abs(interval)) / width) * 21
		binning = rate.NDBins((rate.LinearBins(interval[0], interval[1], bins),))
		self.offsets = rate.BinnedArray(binning)
		self.coinc_offsets = rate.BinnedArray(binning)

	def add_contents(self, contents):
		# this outer loop assumes each injection can only be found
		# in at most one coinc, otherwise the "found" count is
		# wrong.
		for values in contents.connection.cursor().execute("""
SELECT
	sim_burst.*,
	coinc_event.coinc_event_id
FROM
	coinc_event
	JOIN coinc_event_map ON (
		coinc_event_map.coinc_event_id == coinc_event.coinc_event_id
	)
	JOIN sim_burst ON (
		coinc_event_map.table_name == 'sim_burst'
		AND coinc_event_map.event_id == sim_burst.simulation_id
	)
WHERE
	coinc_def_id == ?
		""", (contents.sb_definer_id,)):
			sim = contents.sim_burst_table.row_from_cols(values)
			coinc_event_id = values[-1]
			sim_peak = sim.time_at_instrument(self.instrument)
			self.found += 1
			bursts = tuple(SnglBurstUtils.coinc_sngl_bursts(contents, coinc_event_id))
			coinc_dt = 0
			for burst in bursts:
				dt = float(burst.peak - sim_peak)
				if burst.ifo == self.instrument:
					try:
						self.offsets[dt,] += 1.0
					except IndexError:
						# outside plot range
						pass
				coinc_dt += dt * burst.ms_snr
			coinc_dt /= sum(burst.ms_snr for burst in bursts)
			try:
				self.coinc_offsets[coinc_dt,] += 1.0
			except IndexError:
				# outside plot range
				pass

	def finish(self):
		self.axes.set_title("Trigger Peak Time - Injection Peak Time\n(%d Found Injections)" % self.found)
		# 21 bins per filter width
		filter = rate.gaussian_window(21)
		rate.to_moving_mean_density(self.offsets, filter)
		rate.to_moving_mean_density(self.coinc_offsets, filter)
		self.axes.plot(self.offsets.centres()[0], self.offsets.array, "k")
		self.axes.plot(self.coinc_offsets.centres()[0], self.coinc_offsets.array, "r")
		self.axes.legend(["%s residuals" % self.instrument, "SNR-weighted mean of residuals in all instruments"], loc = "lower right")


#
# =============================================================================
#
#                          Recovered Frequency Offset
#
# =============================================================================
#


class RecoveredFrequencyOffset(object):
	def __init__(self, instrument, interval, width):
		self.fig, self.axes = SnglBurstUtils.make_burst_plot(r"$f_{\mathrm{recovered}} / f_{\mathrm{injected}}$", "Event Number Density")
		self.axes.loglog()
		self.instrument = instrument
		self.found = 0
		# 21 bins per filter width
		bins = int(float(abs(interval)) / width) * 21
		binning = rate.NDBins((rate.LinearBins(interval[0], interval[1], bins),))
		self.offsets = rate.BinnedArray(binning)
		self.coinc_offsets = rate.BinnedArray(binning)

	def add_contents(self, contents):
		# this outer loop assumes each injection can only be found
		# in at most one coinc, otherwise the "found" count is
		# wrong.
		for coinc_event_id, sim_frequency in contents.connection.cursor().execute("""
SELECT
	coinc_event.coinc_event_id,
	sim_burst.frequency
FROM
	coinc_event
	JOIN coinc_event_map ON (
		coinc_event_map.coinc_event_id == coinc_event.coinc_event_id
	)
	JOIN sim_burst ON (
		coinc_event_map.table_name == 'sim_burst'
		AND sim_burst.simulation_id == coinc_event_map.event_id
	)
WHERE
	coinc_event.coinc_def_id == ?
		""", (contents.sb_definer_id,)):
			self.found += 1
			bursts = tuple(SnglBurstUtils.coinc_sngl_bursts(contents, coinc_event_id))
			for burst in bursts:
				if burst.ifo == self.instrument:
					df = math.log(burst.peak_frequency / sim_frequency, 10)
					try:
						self.offsets[df,] += 1.0
					except IndexError:
						# outside plot range
						pass
			# snr-weighted mean of peak frequencies
			coinc_freq = sum(burst.peak_frequency * burst.ms_snr for burst in bursts) / sum(burst.ms_snr for burst in bursts)
			df = math.log(coinc_freq / sim_frequency, 10)
			try:
				self.coinc_offsets[df,] += 1.0
			except IndexError:
				# outside plot range
				pass

	def finish(self):
		self.axes.set_title("Trigger Peak Frequency / Injection Centre Frequency\n(%d Found Injections)" % self.found)
		# 21 bins per filter width
		filter = rate.gaussian_window(21)
		rate.to_moving_mean_density(self.offsets, filter)
		rate.to_moving_mean_density(self.coinc_offsets, filter)
		self.axes.plot(10**self.offsets.centres()[0], self.offsets.array, "k")
		self.axes.plot(10**self.coinc_offsets.centres()[0], self.coinc_offsets.array, "r")
		self.axes.legend(["%s triggers" % self.instrument, "SNR-weighted mean of all matching triggers"], loc = "lower right")
		ymin, ymax = self.axes.get_ylim()
		if ymax / ymin > 1e6:
			ymin = ymax / 1e6
			self.axes.set_ylim((ymin, ymax))


#
# =============================================================================
#
#                       Recovered vs. Injected Frequency
#
# =============================================================================
#


class RecoveredVsInjectedFreq(object):
	def __init__(self, instrument, amplitude_func):
		self.fig, self.axes = SnglBurstUtils.make_burst_plot(r"$f_{\mathrm{injected}}$ (Hz)", r"$f_{\mathrm{recovered}}$ (Hz)")
		self.axes.loglog()
		self.fig.set_size_inches(8, 8)
		self.instrument = instrument
		self.amplitude_func = amplitude_func
		self.matches = 0
		self.x = []
		self.y = []
		self.c = []

	def add_contents(self, contents):
		offsetvectors = contents.time_slide_table.as_dict()
		for values in contents.connection.cursor().execute("""
SELECT
	sim_burst.*,
	sngl_burst.peak_frequency,
	sngl_burst.ms_hrss
FROM
	sim_burst
	JOIN sim_burst_map ON (
		sim_burst_map.simulation_id == sim_burst.simulation_id
	)
	JOIN sngl_burst ON (
		sngl_burst.event_id == sim_burst_map.event_id
	)
WHERE
	sngl_burst.ifo == ?
		""", (self.instrument,)):
			sim = contents.sim_burst_table.row_from_cols(values)
			freq_rec, hrss_rec = values[-2:]
			self.matches += 1
			self.x.append(sim.frequency)
			self.y.append(freq_rec)
			self.c.append(math.log(hrss_rec / self.amplitude_func(sim, self.instrument, offsetvectors[sim.time_slide_id])))

	def finish(self):
		self.axes.scatter(self.x, self.y, c = self.c, s = 16)
		xmin, xmax = min(self.x), max(self.x)
		ymin, ymax = min(self.y), max(self.y)
		xmin = ymin = min(xmin, ymin)
		xmax = ymax = max(xmax, ymax)
		self.axes.plot([xmin, xmax], [ymin, ymax], "k-")
		self.axes.set_xlim([xmin, xmax])
		self.axes.set_ylim([ymin, ymax])
		self.axes.set_title(r"Recovered Frequency vs.\ Injected Frequency (%d Events Matching Injections)" % self.matches)


class CoincRecoveredVsInjectedFreq(object):
	def __init__(self, instruments):
		self.fig, self.axes = SnglBurstUtils.make_burst_plot(r"$f_{\mathrm{injected}}$ (Hz)", r"$f_{\mathrm{recovered}}$ (Hz)")
		self.axes.loglog()
		self.fig.set_size_inches(8, 8)
		self.instruments = set(instruments)
		self.matches = 0
		self.x = []
		self.y = []
		self.c = []

	def add_contents(self, contents):
		# FIXME: this query doesn't check for the correct
		# instruments (should be done via coinc_def_id)
		for values in contents.connection.cursor().execute("""
SELECT
	sim_burst.*,
	multi_burst.central_freq,
	multi_burst.amplitude
FROM
	sim_burst
	JOIN sim_coinc_map ON (
		sim_coinc_map.simulation_id == sim_burst.simulation_id
	)
	JOIN coinc_event AS sim_coinc_event ON (
		sim_coinc_event.coinc_event_id == sim_coinc_map.coinc_event_id
	)
	JOIN multi_burst ON (
		multi_burst.coinc_event_id == sim_coinc_map.event_id
	)
WHERE
	sim_coinc_event.coinc_def_id == ?
		""", (contents.sce_definer_id,)):
			sim = contents.sim_burst_table.row_from_cols(values)
			freq_rec, hrss_rec = values[-2:]
			self.matches += 1
			self.x.append(sim.frequency)
			self.y.append(freq_rec)
			self.c.append(math.log(hrss_rec / sim.hrss))

	def finish(self):
		self.axes.scatter(self.x, self.y, c = self.c, s = 16)
		xmin, xmax = min(self.x), max(self.x)
		ymin, ymax = min(self.y), max(self.y)
		xmin = ymin = min(xmin, ymin)
		xmax = ymax = max(xmax, ymax)
		self.axes.plot([xmin, xmax], [ymin, ymax], "k-")
		self.axes.set_xlim([xmin, xmax])
		self.axes.set_ylim([ymin, ymax])
		self.axes.set_title(r"Recovered Frequency vs.\ Injected Frequency (%d Events Matching Injections)" % self.matches)


#
# =============================================================================
#
#                                     Plot
#
# =============================================================================
#


#
# How to create new plots
#


def new_plots(instrument, amplitude_func, amplitude_lbl, plots):
	l = (
		FreqVsTime(instrument),
		HrssVsFreqScatter(instrument, amplitude_func, amplitude_lbl),
		SimBurstUtils.Efficiency_hrss_vs_freq((instrument,), amplitude_func, amplitude_lbl, 0.1),
		TriggerCountHistogram(instrument),
		RecoveredVsInjectedhrss(instrument, amplitude_func, amplitude_lbl),
		RecoveredPerInjectedhrssVsFreq(instrument, amplitude_func, amplitude_lbl),
		RecoveredPerInjectedhrssVsBandwidth(instrument, amplitude_func, amplitude_lbl),
		RecoveredTimeOffset(instrument, segments.segment(-0.03, +0.03), 0.00015),
		RecoveredFrequencyOffset(instrument, segments.segment(-1.0, +1.0), .002),
		RecoveredVsInjectedFreq(instrument, amplitude_func)
	)
	return [l[i] for i in plots]


def new_coinc_plots(instruments, amplitude_func, amplitude_lbl, plots):
	l = (
		SimBurstUtils.Efficiency_hrss_vs_freq(instruments, amplitude_func, amplitude_lbl, 0.1),
		CoincRecoveredVsInjectedFreq(instruments)
	)
	return [l[i] for i in plots]


#
# Parse command line
#


options, filenames = parse_command_line()
if not options.plot and not options.coinc_plot or not filenames:
	print >>sys.stderr, "Nothing to do!"
	sys.exit(0)


#
# Process files
#


plots = {}
coincplots = new_coinc_plots(("H1", "H2", "L1"), options.amplitude_func, options.amplitude_lbl, options.coinc_plot)

for n, filename in enumerate(utils.sort_files_by_size(filenames, options.verbose, reverse = True)):
	if options.verbose:
		print >>sys.stderr, "%d/%d: %s" % (n + 1, len(filenames), filename)
	working_filename = dbtables.get_connection_filename(filename, tmp_path = options.tmp_space, verbose = options.verbose)
	database = SnglBurstUtils.CoincDatabase(sqlite3.connect(working_filename), options.live_time_program)
	if options.verbose:
		SnglBurstUtils.summarize_coinc_database(database)
	if options.plot:
		if database.coinc_table is not None:
			database.connection.cursor().execute("""
CREATE TEMPORARY VIEW
	sim_burst_map
AS
	SELECT
		a.event_id AS simulation_id,
		a.coinc_event_id AS coinc_event_id,
		b.event_id AS event_id
	FROM
		coinc_event_map AS a
		JOIN coinc_event_map AS b ON (
			a.table_name == 'sim_burst'
			AND b.table_name == 'sngl_burst'
			AND b.coinc_event_id == a.coinc_event_id
		)
			""")
		for instrument in database.instruments:
			if instrument not in plots:
				plots[instrument] = new_plots(instrument, options.amplitude_func, options.amplitude_lbl, options.plot)
			for n, plot in zip(options.plot, plots[instrument]):
				if options.verbose:
					print >>sys.stderr, "adding to %s plot %d ..." % (instrument, n)
				plot.add_contents(database)
	if options.coinc_plot:
		database.connection.cursor().execute("""
CREATE TEMPORARY VIEW
	sim_coinc_map
AS
	SELECT
		a.event_id AS simulation_id,
		a.coinc_event_id AS coinc_event_id,
		b.event_id AS event_id
	FROM
		coinc_event_map AS a
		JOIN coinc_event_map AS b ON (
			a.table_name == 'sim_burst'
			AND b.table_name == 'coinc_event'
			AND b.coinc_event_id == a.coinc_event_id
		)
		""")
		for n, plot in enumerate(coincplots):
			if options.verbose:
				print >>sys.stderr, "adding to coinc plot %d ..." % options.coinc_plot[n]
			plot.add_contents(database)
	database.connection.close()
	dbtables.discard_connection_filename(filename, working_filename, verbose = options.verbose)


#
# compute the binning for the efficiency contour plots
#


def make_binning(plots):
	plots = [plot for instrument in plots.keys() for plot in plots[instrument] if isinstance(plot, SimBurstUtils.Efficiency_hrss_vs_freq)]
	if not plots:
		return None
	minx = min([min(plot.injected_x) for plot in plots])
	maxx = max([max(plot.injected_x) for plot in plots])
	miny = min([min(plot.injected_y) for plot in plots])
	maxy = max([max(plot.injected_y) for plot in plots])
	return rate.NDBins((rate.LogarithmicBins(minx, maxx, 512), rate.LogarithmicBins(miny, maxy, 512)))


binning = make_binning(plots)


#
# finish and write regular plots, deleting them as we go to save memory
#


efficiencies = []
for instrument in plots:
	n = 0
	format = "%%s%s_%%0%dd.%%s" % (instrument, int(math.log10(max(options.plot) or 1)) + 1)
	while len(plots[instrument]):
		plot = plots[instrument].pop(0)
		filename = format % (options.base, options.plot[n], options.format)
		if options.verbose:
			print >>sys.stderr, "finishing %s plot %d ..." % (instrument, options.plot[n])
		try:
			if isinstance(plot, SimBurstUtils.Efficiency_hrss_vs_freq):
				plot.finish(binning = binning)
				efficiencies.append(plot)
				fig = SimBurstUtils.plot_Efficiency_hrss_vs_freq(plot)
			else:
				plot.finish()
				fig = plot.fig
		except ValueError as e:
			print >>sys.stderr, "can't finish %s plot %d: %s" % (instrument, options.plot[n], str(e))
		else:
			if options.verbose:
				print >>sys.stderr, "writing %s ..." % filename
			fig.savefig(filename)
		n += 1


#
# finish and write coinc plots, deleting them as we go to save memory
#


for n, plot in enumerate(coincplots):
	format = "%%s%s_%%0%dd.%%s" % ("coinc", int(math.log10(max(options.coinc_plot) or 1)) + 1)
	filename = format % (options.base, options.coinc_plot[n], options.format)
	if options.verbose:
		print >>sys.stderr, "finishing coinc plot %d ..." % options.coinc_plot[n]
	try:
		if isinstance(plot, SimBurstUtils.Efficiency_hrss_vs_freq):
			plot.finish(binning = binning)
			fig = SimBurstUtils.plot_Efficiency_hrss_vs_freq(plot)
		else:
			plot.finish()
			fig = plot.fig
	except ValueError as e:
		print >>sys.stderr, "can't finish coinc plot %d: %s" % (options.coinc_plot[n], str(e))
	else:
		if options.verbose:
			print >>sys.stderr, "writing %s ..." % filename
		fig.savefig(filename)


#
# generate and write theoretical coincident detection efficiency
#


def plot_multi_Efficiency_hrss_vs_freq(efficiencies):
	fig, axes = SnglBurstUtils.make_burst_plot("Frequency (Hz)", r"$h_{\mathrm{rss}}$")
	axes.loglog()

	e = efficiencies[0]
	xcoords, ycoords = e.efficiency.centres()
	zvals = e.efficiency.ratio()
	error = e.error
	for n, e in enumerate(efficiencies[1:]):
		error += e.error
		other_xcoords, other_ycoords = e.efficiency.centres()
		if (xcoords != other_xcoords).any() or (ycoords != other_ycoords).any():
			# binnings don't match, can't compute product of
			# efficiencies
			raise ValueError("binning mismatch")
		zvals *= e.efficiency.ratio()
	error /= len(efficiencies)

	nfound = numpy.array([len(e.found_x) for e in efficiencies], dtype = "double")
	ninjected = numpy.array([len(e.injected_x) for e in efficiencies], dtype = "double")

	# the model for guessing the ninjected in the coincidence case is
	# to assume that the injections done in the instrument with the
	# most injections were done into all three with a probability given
	# by the ratio of the number actually injected into each instrument
	# to the number injected into the instrument with the most
	# injected, and then to assume that these are independent random
	# occurances and that to be done in coincidence an injection must
	# be done in all three instruments.

	ninjected_guess = (ninjected / ninjected.max()).prod() * ninjected.min()

	# the model for guessing the nfound in the coincidence case is to
	# assume that each injection is found or missed in each instrument
	# at random, and to be found in the coincidence case it must be
	# found in all three.

	nfound_guess = (nfound / ninjected).prod() * ninjected_guess

	ninjected_guess = int(round(ninjected_guess))
	nfound_guess = int(round(nfound_guess))

	instruments = r" \& ".join(sorted("+".join(sorted(e.instruments)) for e in efficiencies))

	cset = axes.contour(xcoords, ycoords, numpy.transpose(zvals), (.1, .2, .3, .4, .5, .6, .7, .8, .9))
	cset.clabel(inline = True, fontsize = 5, fmt = r"$%%g \pm %g$" % error, colors = "k")

	axes.set_title(r"%s Estimated Injection Detection Efficiency ($\sim$%d of $\sim$%d Found)" % (instruments, nfound_guess, ninjected_guess))

	return fig


if efficiencies:
	if options.verbose:
		print >>sys.stderr, "computing theoretical coincident detection efficiency ..."
	fig = plot_multi_Efficiency_hrss_vs_freq(efficiencies)
	filename = "%scoincidence.png" % options.base
	if options.verbose:
		print >>sys.stderr, "writing %s ..." % filename
	fig.savefig(filename)
