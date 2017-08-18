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


"""
String cusp search final output rendering tool.
"""


import bisect
import copy_reg
import heapq
import itertools
from optparse import OptionParser
import math
from matplotlib import patches
import numpy
import os
import cPickle as pickle
import random
import select
from scipy import interpolate
from scipy import optimize
import Queue as queue
import sqlite3
import sys
import traceback


from glue import segments
from glue import segmentsUtils
from glue.ligolw import dbtables
from glue.ligolw.utils import process as ligolwprocess
import lal
from lal import rate
from lal.utils import CacheEntry
from lalburst import git_version
from lalburst import packing
from lalburst import SimBurstUtils
from lalburst import SnglBurstUtils
from lalburst import stringutils


SnglBurstUtils.matplotlib.rcParams.update({
	"font.size": 10.0,
	"text.usetex": True,
	"axes.titlesize": 10.0,
	"axes.labelsize": 10.0,
	"xtick.labelsize": 8.0,
	"ytick.labelsize": 8.0,
	"legend.fontsize": 8.0,
	"grid.linestyle": "-"
})


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


copy_reg.pickle(lal.LIGOTimeGPS, lambda gps: (lal.LIGOTimeGPS, (gps.gpsSeconds, gps.gpsNanoSeconds)))


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
		description = "%prog performs the final, summary, stages of the upper-limit string cusp search.  Input consists of a list of all sqlite format database files produced by all injection and non-injection runs of the analysis pipeline.  The file names can be given on the command line and/or provided in a LAL cache file."
	)
	parser.add_option("--cal-uncertainty", metavar = "fraction", type = "float", help = "Set the fractional uncertainty in amplitude due to calibration uncertainty (eg. 0.08).  This option is required, use 0 to disable calibration uncertainty.")
	parser.add_option("--injections-bin-size", metavar = "bins", type = "float", default = 16.7, help = "Set bin width for injection efficiency curves (default = 16.7).")
	parser.add_option("-c", "--input-cache", metavar = "filename", action = "append", help = "Also process the files named in this LAL cache.  See lalapps_path2cache for information on how to produce a LAL cache file.  This option can be given multiple times.")
	parser.add_option("--import-dump", metavar = "filename", action = "append", help = "Import additional rate vs. threshold or efficiency data from this dump file.  Dump files are one of the data products produced by this program.  Whether the file provides rate vs. threshold data or efficiency data will be determined automatically.  This option can be given multiple times")
	parser.add_option("--image-formats", metavar = "ext[,ext,...]", default = "png,pdf", help = "Set list of graphics formats to produce by providing a comma-delimited list of the filename extensions (default = \"png,pdf\").")
	parser.add_option("-p", "--live-time-program", metavar = "program", default = "StringSearch", help = "Set the name, as it appears in the process table, of the program whose search summary entries define the search live time (default = StringSearch).")
	parser.add_option("-o", "--open-box", action = "store_true", help = "Perform open-box analysis.  In a closed-box analysis (the default), information about the events seen at zero-lag is concealed:  the rate vs. threshold plot only shows the rate of events seen in the background, the detection threshold used to measure the efficiency curves is obtained from n-th loudest background event where n is (the integer closest to) the ratio of background livetime to zero-lag livetime, and messages to stdout and stderr that contain information about event counts at zero-lag are silenced.")
	parser.add_option("-t", "--tmp-space", metavar = "path", help = "Path to a directory suitable for use as a work area while manipulating the database file.  The database file will be worked on in this directory, and then moved to the final location when complete.  This option is intended to improve performance when running in a networked environment, where there might be a local disk with higher bandwidth than is available to the filesystem on which the final output will reside.")
	parser.add_option("--vetoes-name", metavar = "name", help = "Set the name of the segment lists to use as vetoes (default = do not apply vetoes).")
	parser.add_option("--detection-threshold", metavar = "log likelihood ratio", type = "float", help = "Override the detection threshold.  Only injection files will be processed, and the efficiency curve measured.")
	parser.add_option("--record-background", metavar = "N", type = "int", default = 10000000, help = "Set the number of background log likelihood ratios to hold in memory for producing the rate vs. threshold plot (default = 10000000).")
	parser.add_option("--record-candidates", metavar = "N", type = "int", default = 100, help = "Set the number of highest-ranked zero-lag candidates to dump to the candidate file (default = 100).")
	parser.add_option("--threads", metavar = "N", type = "int", default = 1, help = "Set the maximum number of parallel threads to use for processing files (default = 1).  Contention for the global Python interpreter lock will throttle the true number that can run.  The number of threads will be automatically adjusted downwards if the number requested exceeds the number of input files.")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()

	if options.cal_uncertainty is None:
		raise ValueError("must set --cal-uncertainty (use 0 to ignore calibration uncertainty)")

	options.image_formats = options.image_formats.split(",")

	if options.input_cache:
		filenames += [CacheEntry(line).path for input_cache in options.input_cache for line in file(input_cache)]

	if options.threads < 1:
		raise ValueError("--threads must be >= 1")

	return options, filenames


#
# =============================================================================
#
#                              Rate vs. Threshold
#
# =============================================================================
#


def ratevsthresh_bounds(background_x, background_y, zero_lag_x, zero_lag_y):
	if len(zero_lag_x):
		if len(background_x):
			minX, maxX = min(min(zero_lag_x), min(background_x)), max(max(zero_lag_x), max(background_x))
			minY, maxY = min(min(zero_lag_y), min(background_y)), max(max(zero_lag_y), max(background_y))
		else:
			minX, maxX = min(zero_lag_x), max(zero_lag_x)
			minY, maxY = min(zero_lag_y), max(zero_lag_y)
	else:
		# don't check for background, if there's no zero-lag and no
		# background we're screwed anyway
		minX, maxX = min(background_x), max(background_x)
		minY, maxY = min(background_y), max(background_y)

	# override the plot's X axis lower bound
	minX = 1e-2	# FIXME:  don't hard-code
	if len(zero_lag_x):
		if len(background_x):
			maxY = max(zero_lag_y[bisect.bisect_left(zero_lag_x, minX)], background_y[bisect.bisect_left(background_x, minX)])
		else:
			maxY = zero_lag_y[bisect.bisect_left(zero_lag_x, minX)]
	else:
		maxY = background_y[bisect.bisect_left(background_x, minX)]

	return minX, maxX, minY, maxY


def expected_count(x, background_x, background_y):
	if x > background_x[-1]:
		return 0
	return background_y[bisect.bisect_left(background_x, x)]


def compress_ratevsthresh_curve(x, y, yerr):
	#
	# construct a background mask to retain the highest-ranked 10,000
	# elements, then every 10th until the 100,000th highest-ranked
	# element, then every 100th after that.  this is for reducing the
	# dataset size so matplotlib can handle it and vector graphics
	# output isn't ridiculous in size.
	#

	mask = numpy.arange(len(x))[::-1]
	mask = (mask < 10000) | ((mask < 100000) & (mask % 10 == 0)) | (mask % 100 == 0)

	return x.compress(mask), y.compress(mask), yerr.compress(mask)


class RateVsThreshold(object):
	def __init__(self, open_box, record_background, record_candidates):
		self.zero_lag = []
		self.background = []
		self.n_background = 0
		self.record_background = record_background
		self.most_significant_background = []
		self.candidates = []
		self.record_candidates = record_candidates
		self.zero_lag_time = 0.0
		self.background_time = 0.0
		self.open_box = open_box

	def __iadd__(self, other):
		self.zero_lag += other.zero_lag
		self.n_background += other.n_background
		self.background[:] = heapq.nlargest(itertools.chain(self.background, other.background), self.record_background)
		self.most_significant_background[:] = heapq.nlargest(itertools.chain(self.most_significant_background, other.most_significant_background), self.record_candidates)
		self.candidates[:] = heapq.nlargest(itertools.chain(self.candidates, other.candidates), self.record_candidates)
		self.zero_lag_time += other.zero_lag_time
		self.background_time += other.background_time
		return self

	def add_contents(self, contents, verbose = False):
		#
		# retrieve the offset vectors, retain only instruments that
		# are available
		#

		zero_lag_time_slides, background_time_slides = SnglBurstUtils.get_time_slides(contents.connection)
		assert len(zero_lag_time_slides) == 1

		#
		# compute the live time
		#

		seglists = contents.seglists - contents.vetoseglists
		self.zero_lag_time += stringutils.time_slides_livetime(seglists, zero_lag_time_slides.values(), 2, clip = contents.coincidence_segments)
		self.background_time += stringutils.time_slides_livetime(seglists, background_time_slides.values(), 2, clip = contents.coincidence_segments)
		if set(("H1", "H2")).issubset(set(contents.seglists)):
			# we have segments for both H1 and H2, remove time
			# when exactly that pair are on
			self.zero_lag_time -= stringutils.time_slides_livetime_for_instrument_combo(seglists, zero_lag_time_slides.values(), ("H1", "H2"), clip = contents.coincidence_segments)
			self.background_time -= stringutils.time_slides_livetime_for_instrument_combo(seglists, background_time_slides.values(), ("H1", "H2"), clip = contents.coincidence_segments)

		#
		# count events
		#

		for ln_likelihood_ratio, instruments, coinc_event_id, peak_time, is_background in contents.connection.cursor().execute("""
SELECT
	coinc_event.likelihood,
	coinc_event.instruments,
	coinc_event.coinc_event_id,
	(
		SELECT
			AVG(sngl_burst.peak_time) + 1e-9 * AVG(sngl_burst.peak_time_ns)
		FROM
			sngl_burst
			JOIN coinc_event_map ON (
				coinc_event_map.coinc_event_id == coinc_event.coinc_event_id
				AND coinc_event_map.table_name == 'sngl_burst'
				AND coinc_event_map.event_id == sngl_burst.event_id
			)
	),
	EXISTS (
		SELECT
			*
		FROM
			time_slide
		WHERE
			time_slide.time_slide_id == coinc_event.time_slide_id
			AND time_slide.offset != 0
	)
FROM
	coinc_event
WHERE
	coinc_event.coinc_def_id == ?
		""", (contents.bb_definer_id,)):
			# likelihood ratio must be listed first to
			# act as the sort key
			record = (ln_likelihood_ratio, contents.filename, coinc_event_id, dbtables.lsctables.instrumentsproperty.get(instruments), peak_time)
			if ln_likelihood_ratio is None:
				# coinc got vetoed (unable to compute
				# likelihood)
				pass
			elif is_background:
				# non-vetoed background
				self.n_background += 1
				if len(self.background) < self.record_background:
					heapq.heappush(self.background, ln_likelihood_ratio)
				else:
					heapq.heappushpop(self.background, ln_likelihood_ratio)
				if len(self.most_significant_background) < self.record_candidates:
					heapq.heappush(self.most_significant_background, record)
				else:
					heapq.heappushpop(self.most_significant_background, record)
			else:
				# non-vetoed zero lag
				self.zero_lag.append(ln_likelihood_ratio)
				if len(self.candidates) < self.record_candidates:
					heapq.heappush(self.candidates, record)
				else:
					heapq.heappushpop(self.candidates, record)

	def finish(self):
		#
		# dump info about highest-ranked zero-lag and background
		# events
		#

		self.most_significant_background.sort()
		self.candidates.sort()

		f = file("string_most_significant_background.txt", "w")
		print >>f, "Highest-Ranked Background Events"
		print >>f, "================================"
		for ln_likelihood_ratio, filename, coinc_event_id, instruments, peak_time in self.most_significant_background:
			print >>f
			print >>f, "%s in %s:" % (str(coinc_event_id), filename)
			print >>f, "Recovered in: %s" % ", ".join(sorted(instruments or []))
			print >>f, "Mean peak time:  %.16g s GPS" % peak_time
			print >>f, "\\log \\Lambda:  %.16g" % ln_likelihood_ratio

		f = file("string_candidates.txt", "w")
		print >>f, "Highest-Ranked Zero-Lag Events"
		print >>f, "=============================="
		if self.open_box:
			for ln_likelihood_ratio, filename, coinc_event_id, instruments, peak_time in self.candidates:
				print >>f
				print >>f, "%s in %s:" % (str(coinc_event_id), filename)
				print >>f, "Recovered in: %s" % ", ".join(sorted(instruments or []))
				print >>f, "Mean peak time:  %.16g s GPS" % peak_time
				print >>f, "\\log \\Lambda:  %.16g" % ln_likelihood_ratio
		else:
			print >>f
			print >>f, "List suppressed:  box is closed"

		#
		# sort the ranking statistics and convert to arrays.
		#

		self.background.sort()
		self.zero_lag.sort()
		self.zero_lag = numpy.array(self.zero_lag, dtype = "double")
		background_x = numpy.array(self.background, dtype = "double")

		#
		# convert counts to rates and their uncertainties
		#

		# background count expected in zero-lag and \sqrt{N}
		# standard deviation
		background_y = numpy.arange(len(background_x), 0.0, -1.0, dtype = "double") / self.background_time * self.zero_lag_time
		background_yerr = numpy.sqrt(background_y)

		# count observed in zero-lag
		zero_lag_y = numpy.arange(len(self.zero_lag), 0.0, -1.0, dtype = "double")

		# compute zero-lag - expected residual in units of \sqrt{N}
		zero_lag_residual = numpy.zeros((len(self.zero_lag),), dtype = "double")
		for i in xrange(len(self.zero_lag)):
			n = expected_count(self.zero_lag[i], background_x, background_y)
			zero_lag_residual[i] = (zero_lag_y[i] - n) / math.sqrt(n)

		# convert counts to rates
		background_y /= self.zero_lag_time
		background_yerr /= self.zero_lag_time
		zero_lag_y /= self.zero_lag_time

		#
		# determine the horizontal and vertical extent of the plot
		#

		if self.open_box:
			minX, maxX, minY, maxY = ratevsthresh_bounds(background_x, background_y, self.zero_lag, zero_lag_y)
		else:
			minX, maxX, minY, maxY = ratevsthresh_bounds(background_x, background_y, [], [])

		#
		# compress the background data
		#

		background_x, background_y, background_yerr = compress_ratevsthresh_curve(background_x, background_y, background_yerr)

		#
		# start the rate vs. threshold plot
		#

		ratefig, axes = SnglBurstUtils.make_burst_plot(r"Ranking Statistic Threshold, $\log \Lambda$", "Event Rate (Hz)", width = 108.0)
		axes.semilogy()
		axes.set_position([0.125, 0.15, 0.83, 0.75])
		axes.xaxis.grid(True, which = "major,minor")
		axes.yaxis.grid(True, which = "major,minor")
		if self.open_box:
			axes.set_title(r"Event Rate vs.\ Ranking Statistic Threshold")
		else:
			axes.set_title(r"Event Rate vs.\ Ranking Statistic Threshold (Closed Box)")

		# warning:  the error bar polygon is not *really* clipped
		# to the axes' bounding box, the result will be incorrect
		# if the density of data points is small where the polygon
		# encounters the axes' bounding box.

		poly_x = numpy.concatenate((background_x, background_x[::-1]))
		poly_y = numpy.concatenate((background_y + 1 * background_yerr, (background_y - 1 * background_yerr)[::-1]))
		axes.add_patch(patches.Polygon(zip(poly_x, numpy.clip(poly_y, minY, maxY)), edgecolor = "k", facecolor = "k", alpha = 0.3))
		line1, = axes.semilogy(background_x.repeat(2)[:-1], background_y.repeat(2)[1:], color = "k", linestyle = "--")
		if self.open_box:
			line2, = axes.semilogy(self.zero_lag.repeat(2)[:-1], zero_lag_y.repeat(2)[1:], color = "k", linestyle = "-", linewidth = 2)
			axes.legend((line1, line2), (r"Expected background", r"Zero-lag"), loc = "lower left")
		else:
			axes.legend((line1,), (r"Expected background",), loc = "lower left")

		axes.set_xlim((minX, maxX))
		axes.set_ylim((minY, maxY))

		#
		# start the count residual vs. threshold plot
		#

		residualfig, axes = SnglBurstUtils.make_burst_plot(r"Ranking Statistic Threshold, $\log \Lambda$", r"Event Count Residual / $\sqrt{N}$", width = 108.0)
		axes.set_position([0.125, 0.15, 0.83, 0.75])
		axes.xaxis.grid(True, which = "major,minor")
		axes.yaxis.grid(True, which = "major,minor")
		if self.open_box:
			axes.set_title(r"Event Count Residual vs.\ Ranking Statistic Threshold")
		else:
			axes.set_title(r"Event Count Residual vs.\ Ranking Statistic Threshold (Closed Box)")

		axes.add_patch(patches.Polygon(((minX, -1), (maxX, -1), (maxX, +1), (minX, +1)), edgecolor = "k", facecolor = "k", alpha = 0.3))
		if self.open_box:
			line1, = axes.plot(self.zero_lag, zero_lag_residual, color = "k", linestyle = "-", linewidth = 2)

		axes.set_xlim((minX, maxX))

		#
		# done
		#

		return ratefig, residualfig


#
# =============================================================================
#
#                                  Efficiency
#
# =============================================================================
#


def slope(x, y):
	"""
	From the x and y arrays, compute the slope at the x co-ordinates
	using a first-order finite difference approximation.
	"""
	slope = numpy.zeros((len(x),), dtype = "double")
	slope[0] = (y[1] - y[0]) / (x[1] - x[0])
	for i in xrange(1, len(x) - 1):
		slope[i] = (y[i + 1] - y[i - 1]) / (x[i + 1] - x[i - 1])
	slope[-1] = (y[-1] - y[-2]) / (x[-1] - x[-2])
	return slope


def upper_err(y, yerr, deltax):
	z = y + yerr
	deltax = int(deltax)
	upper = numpy.zeros((len(yerr),), dtype = "double")
	for i in xrange(len(yerr)):
		upper[i] = max(z[max(i - deltax, 0) : min(i + deltax, len(z))])
	return upper - y


def lower_err(y, yerr, deltax):
	z = y - yerr
	deltax = int(deltax)
	lower = numpy.zeros((len(yerr),), dtype = "double")
	for i in xrange(len(yerr)):
		lower[i] = min(z[max(i - deltax, 0) : min(i + deltax, len(z))])
	return y - lower


def write_efficiency(fileobj, bins, eff, yerr):
	print >>fileobj, "# A	e	D[e]"
	for A, e, De in zip(bins.centres()[0], eff, yerr):
		print >>fileobj, "%.16g	%.16g	%.16g" % (A, e, De)


def render_data_from_bins(dump_file, axes, efficiency_num, efficiency_den, cal_uncertainty, filter_width, colour = "k", erroralpha = 0.3, linestyle = "-"):
	# extract array of x co-ordinates, and the factor by which x
	# increases from one sample to the next.
	(x,) = efficiency_den.centres()
	x_factor_per_sample = efficiency_den.bins[0].delta

	# compute the efficiency, the slope (units = efficiency per
	# sample), the y uncertainty (units = efficiency) due to binomial
	# counting fluctuations, and the x uncertainty (units = samples)
	# due to the width of the smoothing filter.
	eff = efficiency_num.array / efficiency_den.array
	dydx = slope(numpy.arange(len(x), dtype = "double"), eff)
	yerr = numpy.sqrt(eff * (1. - eff) / efficiency_den.array)
	xerr = numpy.array([filter_width / 2.] * len(yerr))

	# compute the net y err (units = efficiency) by (i) multiplying the
	# x err by the slope, (ii) dividing the calibration uncertainty
	# (units = percent) by the fractional change in x per sample and
	# multiplying by the slope, (iii) adding the two in quadradure with
	# the y err.
	net_yerr = numpy.sqrt((xerr * dydx)**2. + yerr**2. + (cal_uncertainty / x_factor_per_sample * dydx)**2.)

	# compute net xerr (units = percent) by dividing yerr by slope and
	# then multiplying by the fractional change in x per sample.
	net_xerr = net_yerr / dydx * x_factor_per_sample

	# write the efficiency data to file
	write_efficiency(dump_file, efficiency_den.bins, eff, net_yerr)

	# plot the efficiency curve and uncertainty region
	patch = patches.Polygon(zip(numpy.concatenate((x, x[::-1])), numpy.concatenate((eff + upper_err(eff, yerr, filter_width / 2.), (eff - lower_err(eff, yerr, filter_width / 2.))[::-1]))), edgecolor = colour, facecolor = colour, alpha = erroralpha)
	axes.add_patch(patch)
	line, = axes.plot(x, eff, colour + linestyle)

	# compute 50% point and its uncertainty
	A50 = optimize.bisect(interpolate.interp1d(x, eff - 0.5), x[0], x[-1], xtol = 1e-40)
	A50_err = interpolate.interp1d(x, net_xerr)(A50)

	# mark 50% point on graph
	#axes.axvline(A50, color = colour, linestyle = linestyle)

	# print some analysis FIXME:  this calculation needs attention
	num_injections = efficiency_den.array.sum()
	num_samples = len(efficiency_den.array)
	print >>sys.stderr, "Bins were %g samples wide, ideal would have been %g" % (filter_width, (num_samples / num_injections / interpolate.interp1d(x, dydx)(A50)**2.0)**(1.0/3.0))
	print >>sys.stderr, "Average number of injections in each bin = %g" % efficiency_den.array.mean()

	return line, A50, A50_err


class Efficiency(object):
	def __init__(self, detection_threshold, cal_uncertainty, filter_width, open_box):
		self.detection_threshold = detection_threshold
		self.cal_uncertainty = cal_uncertainty
		self.filter_width = filter_width
		self.seglists = segments.segmentlistdict()
		self.vetoseglists = segments.segmentlistdict()
		self.found = []
		self.n_diagnostics = 100	# keep 100 loudest missed and quietest found injections
		self.loudest_missed = []
		self.quietest_found = []
		self.all = []
		self.open_box = open_box

	def __iadd__(self, other):
		assert other.detection_threshold == self.detection_threshold
		assert other.open_box == self.open_box
		self.seglists |= other.seglists
		self.vetoseglists |= other.vetoseglists
		self.found += other.found
		self.loudest_missed[:] = heapq.nlargest(itertools.chain(self.loudest_missed, other.loudest_missed), self.n_diagnostics)
		self.quietest_found[:] = heapq.nlargest(itertools.chain(self.quietest_found, other.quietest_found), self.n_diagnostics)
		self.all += other.all
		return self

	def add_contents(self, contents, verbose = False):
		#
		# update segment information
		#

		self.seglists |= contents.seglists
		self.vetoseglists |= contents.vetoseglists
		seglists = contents.seglists - contents.vetoseglists
		if set(("H1", "H2")).issubset(set(seglists)):
			# we have segments for both H1 and H2, remove time
			# when exactly that pair are on
			seglists -= segments.segmentlistdict.fromkeys(seglists, seglists.intersection(("H1", "H2")) - seglists.union(set(seglists) - set(("H1", "H2"))))

		# for each injection, retrieve the highest likelihood ratio
		# of the burst coincs that were found to match the
		# injection or null if no burst coincs matched the
		# injection
		offsetvectors = contents.time_slide_table.as_dict()
		stringutils.create_recovered_ln_likelihood_ratio_table(contents.connection, contents.bb_definer_id)
		for values in contents.connection.cursor().execute("""
SELECT
	sim_burst.*,
	recovered_ln_likelihood_ratio.ln_likelihood_ratio
FROM
	sim_burst
	LEFT JOIN recovered_ln_likelihood_ratio ON (
		recovered_ln_likelihood_ratio.simulation_id == sim_burst.simulation_id
	)
		"""):
			sim = contents.sim_burst_table.row_from_cols(values[:-1])
			ln_likelihood_ratio = values[-1]
			found = ln_likelihood_ratio is not None
			# were at least 2 instruments on when the injection
			# was made?
			if len(SimBurstUtils.on_instruments(sim, seglists, offsetvectors[sim.time_slide_id])) >= 2:
				# yes
				self.all.append(sim)
				if found and ln_likelihood_ratio > self.detection_threshold:
					self.found.append(sim)
					# 1/amplitude needs to be first so
					# that it acts as the sort key
					record = (1.0 / sim.amplitude, sim, offsetvectors[sim.time_slide_id], contents.filename, ln_likelihood_ratio)
					if len(self.quietest_found) < self.n_diagnostics:
						heapq.heappush(self.quietest_found, record)
					else:
						heapq.heappushpop(self.quietest_found, record)
				else:
					# amplitude needs to be first so
					# that it acts as the sort key
					record = (sim.amplitude, sim, offsetvectors[sim.time_slide_id], contents.filename, ln_likelihood_ratio)
					if len(self.loudest_missed) < self.n_diagnostics:
						heapq.heappush(self.loudest_missed, record)
					else:
						heapq.heappushpop(self.loudest_missed, record)
			elif found:
				# no, but it was found anyway
				print >>sys.stderr, "%s: odd, injection %s was found but not injected ..." % (contents.filename, sim.simulation_id)

	def finish(self):
		fig, axes = SnglBurstUtils.make_burst_plot(r"Injection Amplitude (\(\mathrm{s}^{-\frac{1}{3}}\))", "Detection Efficiency", width = 108.0)
		axes.set_title(r"Detection Efficiency vs.\ Amplitude")
		axes.semilogx()
		axes.set_position([0.10, 0.150, 0.86, 0.77])

		# set desired yticks
		axes.set_yticks((0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0))
		axes.set_yticklabels((r"\(0\)", r"\(0.1\)", r"\(0.2\)", r"\(0.3\)", r"\(0.4\)", r"\(0.5\)", r"\(0.6\)", r"\(0.7\)", r"\(0.8\)", r"\(0.9\)", r"\(1.0\)"))
		axes.xaxis.grid(True, which = "major,minor")
		axes.yaxis.grid(True, which = "major,minor")

		# put made and found injections in the denominators and
		# numerators of the efficiency bins
		bins = rate.NDBins((rate.LogarithmicBins(min(sim.amplitude for sim in self.all), max(sim.amplitude for sim in self.all), 400),))
		efficiency_num = rate.BinnedArray(bins)
		efficiency_den = rate.BinnedArray(bins)
		for sim in self.found:
			efficiency_num[sim.amplitude,] += 1
		for sim in self.all:
			efficiency_den[sim.amplitude,] += 1

		# generate and plot trend curves.  adjust window function
		# normalization so that denominator array correctly
		# represents the number of injections contributing to each
		# bin:  make w(0) = 1.0.  note that this factor has no
		# effect on the efficiency because it is common to the
		# numerator and denominator arrays.  we do this for the
		# purpose of computing the Poisson error bars, which
		# requires us to know the counts for the bins
		windowfunc = rate.gaussian_window(self.filter_width)
		windowfunc /= windowfunc[len(windowfunc) / 2 + 1]
		rate.filter_array(efficiency_num.array, windowfunc)
		rate.filter_array(efficiency_den.array, windowfunc)

		# regularize:  adjust unused bins so that the efficiency is
		# 0, not NaN
		assert (efficiency_num.array <= efficiency_den.array).all()
		efficiency_den.array[(efficiency_num.array == 0) & (efficiency_den.array == 0)] = 1

		line1, A50, A50_err = render_data_from_bins(file("string_efficiency.dat", "w"), axes, efficiency_num, efficiency_den, self.cal_uncertainty, self.filter_width, colour = "k", linestyle = "-", erroralpha = 0.2)
		print >>sys.stderr, "Pipeline's 50%% efficiency point for all detections = %g +/- %g%%\n" % (A50, A50_err * 100)

		# add a legend to the axes
		axes.legend((line1,), (r"\noindent Injections recovered with $\log \Lambda > %.2f$" % self.detection_threshold,), loc = "lower right")

		# adjust limits
		axes.set_xlim([3e-22, 3e-19])
		axes.set_ylim([0.0, 1.0])

		#
		# dump some information about the highest-amplitude missed
		# and quietest-amplitude found injections
		#

		self.loudest_missed.sort(reverse = True)
		self.quietest_found.sort(reverse = True)

		f = file("string_loud_missed_injections.txt", "w")
		print >>f, "Highest Amplitude Missed Injections"
		print >>f, "==================================="
		for amplitude, sim, offsetvector, filename, ln_likelihood_ratio in self.loudest_missed:
			print >>f
			print >>f, "%s in %s:" % (str(sim.simulation_id), filename)
			if ln_likelihood_ratio is None:
				print >>f, "Not recovered"
			else:
				print >>f, "Recovered with \\log \\Lambda = %.16g, detection threshold was %.16g" % (ln_likelihood_ratio, self.detection_threshold)
			for instrument in self.seglists:
				print >>f, "In %s:" % instrument
				print >>f, "\tInjected amplitude:\t%.16g" % SimBurstUtils.string_amplitude_in_instrument(sim, instrument, offsetvector)
				print >>f, "\tTime of injection:\t%s s" % sim.time_at_instrument(instrument, offsetvector)
			print >>f, "Amplitude in waveframe:\t%.16g" % sim.amplitude
			t = sim.get_time_geocent()
			print >>f, "Time at geocentre:\t%s s" % t
			print >>f, "Segments within 60 seconds:\t%s" % segmentsUtils.segmentlistdict_to_short_string(self.seglists & segments.segmentlistdict((instrument, segments.segmentlist([segments.segment(t-offsetvector[instrument]-60, t-offsetvector[instrument]+60)])) for instrument in self.seglists))
			print >>f, "Vetoes within 60 seconds:\t%s" % segmentsUtils.segmentlistdict_to_short_string(self.vetoseglists & segments.segmentlistdict((instrument, segments.segmentlist([segments.segment(t-offsetvector[instrument]-60, t-offsetvector[instrument]+60)])) for instrument in self.vetoseglists))

		f = file("string_quiet_found_injections.txt", "w")
		print >>f, "Lowest Amplitude Found Injections"
		print >>f, "================================="
		for inv_amplitude, sim, offsetvector, filename, ln_likelihood_ratio in self.quietest_found:
			print >>f
			print >>f, "%s in %s:" % (str(sim.simulation_id), filename)
			if ln_likelihood_ratio is None:
				print >>f, "Not recovered"
			else:
				print >>f, "Recovered with \\log \\Lambda = %.16g, detection threshold was %.16g" % (ln_likelihood_ratio, self.detection_threshold)
			for instrument in self.seglists:
				print >>f, "In %s:" % instrument
				print >>f, "\tInjected amplitude:\t%.16g" % SimBurstUtils.string_amplitude_in_instrument(sim, instrument, offsetvector)
				print >>f, "\tTime of injection:\t%s s" % sim.time_at_instrument(instrument, offsetvector)
			print >>f, "Amplitude in waveframe:\t%.16g" % sim.amplitude
			t = sim.get_time_geocent()
			print >>f, "Time at geocentre:\t%s s" % t
			print >>f, "Segments within 60 seconds:\t%s" % segmentsUtils.segmentlistdict_to_short_string(self.seglists & segments.segmentlistdict((instrument, segments.segmentlist([segments.segment(t-offsetvector[instrument]-60, t-offsetvector[instrument]+60)])) for instrument in self.seglists))
			print >>f, "Vetoes within 60 seconds:\t%s" % segmentsUtils.segmentlistdict_to_short_string(self.vetoseglists & segments.segmentlistdict((instrument, segments.segmentlist([segments.segment(t-offsetvector[instrument]-60, t-offsetvector[instrument]+60)])) for instrument in self.vetoseglists))

		#
		# done
		#

		return fig,


#
# =============================================================================
#
#                                     Main
#
# =============================================================================
#


def group_files(filenames, verbose = False):
	# figure out which files are injection runs and which aren't
	injection_filenames = []
	noninjection_filenames = []
	for n, filename in enumerate(filenames):
		if verbose:
			print >>sys.stderr, "%d/%d: %s" % (n + 1, len(filenames), filename)
		connection = sqlite3.connect(filename)
		if "sim_burst" in dbtables.get_table_names(connection):
			if verbose:
				print >>sys.stderr, "\t--> injections"
			injection_filenames.append(filename)
		else:
			if verbose:
				print >>sys.stderr, "\t--> non-injections"
			noninjection_filenames.append(filename)
		connection.close()
	return injection_filenames, noninjection_filenames


def pack_files(filenames, n_threads, verbose = False):
	bins = packing.BiggestIntoEmptiest([packing.Bin() for n in range(n_threads)])
	bins.packlist((os.stat(filename).st_size, filename) for filename in filenames)
	return [bin.objects for bin in bins.bins]


def import_dumps(filenames, verbose = False):
	rate_vs_threshold = None
	efficiency = None
	for filename in filenames:
		if verbose:
			print >>sys.stderr, "\t%s ..." % filename,
		dump = pickle.load(open(filename))
		if type(dump) is RateVsThreshold:
			if verbose:
				print >>sys.stderr, "found rate vs. threshold data"
			if rate_vs_threshold is None:
				rate_vs_threshold = dump
			else:
				rate_vs_threshold += dump
		elif type(dump) is Efficiency:
			if verbose:
				print >>sys.stderr, "found efficiency data"
			if efficiency is None:
				efficiency = dump
			else:
				efficiency += dump
		else:
			raise ValueError("cannot determine contents of %s" % filename)
	return rate_vs_threshold, efficiency


def process_file(filename, products, live_time_program, tmp_path = None, veto_segments_name = None, verbose = False):
	#
	# connect to database and summarize contents
	#

	working_filename = dbtables.get_connection_filename(filename, tmp_path = tmp_path, verbose = verbose)
	contents = SnglBurstUtils.CoincDatabase(sqlite3.connect(working_filename), live_time_program, search = "StringCusp", veto_segments_name = veto_segments_name)
	if verbose:
		SnglBurstUtils.summarize_coinc_database(contents, filename = working_filename)

	#
	# augment summary with extra stuff we need.  the filename
	# is recorded for dumping debuggin information related to
	# missed injections.  if burca was run with the
	# --coincidence-segments option then the value is copied
	# into a segmentlistdict to facilitate the computation of
	# livetime
	#

	contents.filename = filename

	contents.coincidence_segments = ligolwprocess.get_process_params(contents.xmldoc, "lalapps_burca", "--coincidence-segments")
	if contents.coincidence_segments:
		# as a side-effect, this enforces the rule that
		# burca has been run on the input file exactly once
		contents.coincidence_segments, = contents.coincidence_segments
		contents.coincidence_segments = segments.segmentlistdict.fromkeys(contents.seglists, segmentsUtils.from_range_strings(contents.coincidence_segments.split(","), boundtype = dbtables.lsctables.LIGOTimeGPS).coalesce())
	else:
		contents.coincidence_segments = None

	#
	# process contents
	#

	for n, product in enumerate(products):
		if verbose:
			print >>sys.stderr, "%s: adding to product %d ..." % (working_filename, n)
		product.add_contents(contents, verbose = verbose)

	#
	# close
	#

	contents.connection.close()
	dbtables.discard_connection_filename(filename, working_filename, verbose = verbose)


def process_files(filenames, products, live_time_program, tmp_path = None, veto_segments_name = None, verbose = False):
	for n, filename in enumerate(filenames):
		if verbose:
			print >>sys.stderr, "%d/%d: %s" % (n + 1, len(filenames), filename)
		process_file(filename, products, live_time_program, tmp_path = tmp_path, veto_segments_name = veto_segments_name, verbose = verbose)


def write_products(products, prefix, image_formats, verbose = False):
	format = "%%s%%0%dd%%s.%%s" % (int(math.log10(max(len(products) - 1, 1))) + 1)
	n = 1
	while products:
		if verbose:
			print >>sys.stderr, "finishing product %d ..." % (n - 1)
		product = products.pop(0)
		# write dump of raw data
		filename = "%sdump_%d.pickle" % (prefix, n)
		if verbose:
			print >>sys.stderr, "\twriting %s ..." % filename,
		pickle.dump(product, open(filename, "w"))
		if verbose:
			print >>sys.stderr, "done"
		# write plots
		for m, fig in enumerate(product.finish()):
			for ext in image_formats:
				filename = format % (prefix, n, chr(m + 97), ext)
				if verbose:
					print >>sys.stderr, "\twriting %s ..." % filename,
				fig.savefig(filename)
				if verbose:
					print >>sys.stderr, "done"
		n += 1


options, filenames = parse_command_line()


print >>sys.stderr, """Command line:

$ %s
""" % " ".join(sys.argv)
if options.open_box:
	print >>sys.stderr, """

---=== !! BOX IS OPEN !! ===---

PRESS CTRL-C SOON IF YOU DIDN'T MEAN TO OPEN THE BOX

"""


#
# figure out which files are from injection runs and which aren't
#


if options.verbose:
	print >>sys.stderr, "Identifying injection and non-injection databases ..."
injection_filenames, noninjection_filenames = group_files(filenames, verbose = options.verbose)


#
# intialize the data collection objects
#


if options.import_dump:
	if options.verbose:
		print >>sys.stderr, "Loading dump files ..."
	rate_vs_threshold, efficiency = import_dumps(options.import_dump, verbose = options.verbose)
	# override box openness in rate-vs-threshold data
	if rate_vs_threshold is not None:
		rate_vs_threshold.open_box = options.open_box
	if efficiency is not None and options.open_box and not efficiency.open_box:
		raise ValueError("Box is open but one or more imjported efficiency dump file was generated in closed-box mode.  Efficiency must be re-measured in open-box mode to use correct detection threshold.")
else:
	rate_vs_threshold, efficiency = None, None


#
# collect zero-lag and background statistics
#


if options.detection_threshold is None:
	if options.verbose:
		print >>sys.stderr, "Collecting background and zero-lag statistics ..."

	children = {}
	rate_vs_thresholds = []
	# group files into bins of about the same total number of bytes.
	# don't try to create more groups than there are files
	for filenames in pack_files(noninjection_filenames, min(options.threads, len(noninjection_filenames)), verbose = options.verbose):
		# shuffle file names to avoid copying all the big files
		# into the scratch space at the same time
		random.shuffle(filenames)
		# launch thread
		r, w = os.pipe()
		r, w = os.fdopen(r, "r", 0), os.fdopen(w, "w", 0)
		pid = os.fork()
		if pid == 0:
			# we're the child process
			r.close()
			try:
				# create a new book-keeping object
				rate_vs_threshold = RateVsThreshold(options.open_box, record_background = options.record_background, record_candidates = options.record_candidates)
				process_files(filenames, [rate_vs_threshold], options.live_time_program, tmp_path = options.tmp_space, veto_segments_name = options.vetoes_name, verbose = options.verbose)
				pickle.dump(rate_vs_threshold, w)
			except:
				pickle.dump(traceback.format_exc(), w)
				w.close()
				sys.exit(1)
			w.close()
			sys.exit(0)
		# we're not the child process
		w.close()
		children[r] = pid
	# collect all children, report any exceptions that occured, combine
	# results
	while children:
		for r in select.select(children.keys(), [], [])[0]:
			process_output = pickle.load(r)
			r.close()
			pid, exit_status = os.waitpid(children.pop(r), 0)
			if isinstance(process_output, RateVsThreshold):
				if rate_vs_threshold is None:
					rate_vs_threshold = process_output
				else:
					rate_vs_threshold += process_output
			else:
				print >>sys.stderr, process_output
			del process_output
			if exit_status:
				sys.exit(exit_status)
	if rate_vs_threshold is None:
		raise ValueError("no non-injection data available:  cannot determine detection threshold.  provide non-injection data or set an explicit detection treshold.  Consult --help for more information.")

	# write data products
	write_products([rate_vs_threshold], "string_rate_", options.image_formats, verbose = options.verbose)

	# print summary information, and set detection threshold
	if options.open_box:
		try:
			print >>sys.stderr, "Zero-lag events: %d" % len(rate_vs_threshold.zero_lag)
		except OverflowError:
			# python < 2.5 can't handle list-like things with more than 2**31 elements
			# FIXME:  remove when we can rely on python >= 2.5
			print >>sys.stderr, "Zero-lag events: %d" % rate_vs_threshold.zero_lag.n
	print >>sys.stderr, "Total time in zero-lag segments: %s s" % str(rate_vs_threshold.zero_lag_time)
	print >>sys.stderr, "Time-slide events: %d" % rate_vs_threshold.n_background
	print >>sys.stderr, "Total time in time-slide segments: %s s" % str(rate_vs_threshold.background_time)
	if options.open_box:
		detection_threshold = rate_vs_threshold.zero_lag[-1]
		print >>sys.stderr, "Likelihood ratio for highest-ranked zero-lag survivor: %.9g" % detection_threshold
	else:
		# if the background and zero-lag live times are identical,
		# then the loudest zero-lag event is simulated by the
		# loudest background event so we want entry -1 in the
		# sorted list; if the background livetime is twice the
		# zero-lag live time then the loudest zero-lag event is
		# simulated by the next-to-loudest background event so we
		# want entry -2 in the sorted list.
		detection_threshold = sorted(rate_vs_threshold.background)[-int(round(rate_vs_threshold.background_time / rate_vs_threshold.zero_lag_time))]
		print >>sys.stderr, "Simulated \\log \\Lambda for highest-ranked zero-lag survivor: %.9g" % detection_threshold
else:
	detection_threshold = options.detection_threshold
	print >>sys.stderr, "Likelihood ratio for highest-ranked zero-lag survivor from command line: %.9g" % detection_threshold


#
# measure detection efficiency based on loudest (simulated, if box is
# closed) zero-lag survivor
#


if options.verbose:
	print >>sys.stderr, "Collecting efficiency statistics ..."
children = {}
# group files into bins of about the same total number of bytes.  don't try
# to create more groups than there are files
for filenames in pack_files(injection_filenames, min(options.threads, len(injection_filenames)), verbose = options.verbose):
	# shuffle file names to avoid copying all the big files
	# into the scratch space at the same time
	random.shuffle(filenames)
	# launch thread
	r, w = os.pipe()
	r, w = os.fdopen(r, "r", 0), os.fdopen(w, "w", 0)
	pid = os.fork()
	if pid == 0:
		# we're the child process
		r.close()
		try:
			# create a new book-keeping object
			efficiency = Efficiency(detection_threshold, options.cal_uncertainty, options.injections_bin_size, options.open_box)
			process_files(filenames, [efficiency], options.live_time_program, tmp_path = options.tmp_space, veto_segments_name = options.vetoes_name, verbose = options.verbose)
			pickle.dump(efficiency, w)
		except:
			pickle.dump(traceback.format_exc(), w)
			w.close()
			sys.exit(1)
		w.close()
		sys.exit(0)
	# we're not the child process
	w.close()
	children[r] = pid
# collect all children, report any exceptions that occured, combine
# results
while children:
	for r in select.select(children.keys(), [], [])[0]:
		process_output = pickle.load(r)
		r.close()
		pid, exit_status = os.waitpid(children.pop(r), 0)
		if isinstance(process_output, Efficiency):
			if efficiency is None:
				efficiency = process_output
			else:
				efficiency += process_output
		else:
			print >>sys.stderr, process_output
		del process_output
		if exit_status:
			sys.exit(exit_status)
# write data products
if efficiency is not None:
	write_products([efficiency], "string_efficiency_", options.image_formats, verbose = options.verbose)
else:
	if options.verbose:
		print >>sys.stderr, "no injection data available, not writing efficiency data products."

if options.verbose:
	print >>sys.stderr, "done."
