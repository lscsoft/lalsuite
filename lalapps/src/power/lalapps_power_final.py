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


# Information extracted from data set:
#
# - from zero lag obtain zero lag event rate vs. amplitude threshold curve,
# and amplitude of loudest event to use as final threshold
#
# - from time slides, obtain background event rate vs. amplitude threshold
# curve.  this is \mu_{0}(x) in Brady et al. upper limit paper.
#
# - numerically differentiate that curve to obtain \mu_{0}'(x).  this is
# done by first computing a spline approximation of the \mu_{0}(x) curve,
# and then differentiating that curve.
#
# - from histogram of the times between events, confirm that background is
# a Poisson process.  but should this be of all background events, or only
# background events with amplitudes above the final threshold?  again, the
# latter could be a problem if there aren't enough events in the time
# slides above threshold
#
# - in each M_{sun} vs. centre frequency bin, measure the injection
# recovery probability at the final threshold.  this is \epsilon(x) in
# Brady et al., and required two 2D M_{sun} vs. frequency binnings to
# compute.  all made injections go into the first 2D binning, and all
# injections recovered with an amplitude above threshold go into the second
# 2D binning.  the ratio of the second binning's counts to those in the
# first yields the detection efficiency (found / made).
#
# - in each M_{sun} vs. centre frequency bin, measure the injection
# recovery probability density at the final threshold.  this is
# \epsilon'(x) in Brady et al., and requires one additional 2D M_{sun} vs.
# frequency binning to compute.  all recovered injections regardless of
# amplitude go into this binning weighted by an amplitude window function.
# the ratio of this third binnings' counts to those in the first above
# yields the efficiency density if the amplitude window has been normalized
# properly (so as to yield count per amplitude interval).


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


import bisect
try:
	from fpconst import PosInf, NegInf
except ImportError:
	# fpconst is not part of the standard library and might not be
	# available
	PosInf = float("+inf")
	NegInf = float("-inf")
import glob
import math
from matplotlib import patches
import numpy
from scipy import interpolate
from scipy import optimize
from optparse import OptionParser
import sqlite3
import sys


from glue import iterutils
from glue import segments
from glue.ligolw import dbtables
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
		version = "Name: %%prog\n%s" % git_version.verbose_msg,
		usage = "%prog [options] -i|--injections-glob pattern -b|--background-glob pattern",
		description = "%prog performs the final, summary, stages of an upper-limit excess power search for burst-like gravitational waves."
	)
	parser.add_option("-b", "--background-glob", metavar = "pattern", default = [], action = "append", help = "Shell filename pattern for non-injection files (required).")
	parser.add_option("-i", "--injections-glob", metavar = "pattern", default = [], action = "append", help = "Shell filename pattern for injection files (required).")
	parser.add_option("-l", "--live-time-program", metavar = "program", default = "lalapps_power", help = "Set the name, as it appears in the process table, of the program whose search summary \"out\" segments define the search live time (default = \"lalapps_power\")")
	parser.add_option("--open-box", action = "store_true", help = "Open the box.  If this is not set, then instead of setting a threshold based on the loudest zero-lag coincidence, the threshold is set at the predicted loudest-event amplitude as derived from the background.")
	parser.add_option("--confidence-contour-slope", metavar = "slope", type = "float", default = -12.0, help = "Set the slope of the confidence-likelihood joint threshold contours (default = -12.0).")
	parser.add_option("--dump-scatter-data", action = "store_true", help = "Instead of computing upper limit, dump data files for confidence-likelihood scatter diagnostic scatter plot.")
	parser.add_option("--plot-scatter-data", action = "store_true", help = "Instead of computing upper limit, load data files for confidence-likelihood scatter diagnostic scatter plot and generate scatter plot.")
	parser.add_option("-t", "--tmp-space", metavar = "path", help = "Path to a directory suitable for use as a work area while manipulating the database file.  The database file will be worked on in this directory, and then moved to the final location when complete.  This option is intended to improve performance when running in a networked environment, where there might be a local disk with higher bandwidth than is available to the filesystem on which the final output will reside.")
	parser.add_option("--upper-limit-confidence", metavar = "percent", type = "float", default = 90, help = "Set the confidence in percent of the upper limit to derive (default = 90.0).")
	parser.add_option("--upper-limit-scale", metavar = "E|hrss", default = "E", help = "Select the scale on which to measure the upper limit.  Allowed values are 'E' for an upper limit on rate of bursts by energy, and 'hrss' for an upper limit on rate of bursts by h_{rss}.  (default = 'E')")
	parser.add_option("--zero-lag-survivors", metavar = "number", type = "int", default = 0, help = "Set the final confidence-likelihood joint threshold so that this many events survive at zero lag.  If n survivors are requested, then the threshold is set at exactly the amplitude of the n+1 loudest event, and events are discarded if their amplitudes are <= the threshold.  A \"loudest-event\" upper limit is obtained by setting this to 0.  The default is 0.")
	parser.add_option("-v", "--verbose", action = "store_true", help = "be verbose")
	options, filenames = parser.parse_args()

	if not options.plot_scatter_data:
		if not options.background_glob:
			raise ValueError("missing required --background-glob argument")
		if not options.injections_glob:
			raise ValueError("missing required --injections-glob argument")

	options.upper_limit_confidence /= 100.0

	if options.upper_limit_scale not in ("E", "hrss"):
		raise ValueError("invalid value '%s' for --upper-limit-scale" % options.upper_limit_scale)

	return options, (filenames or [None])


#
# =============================================================================
#
#                              Re-usable Queries
#
# =============================================================================
#


#
# An iterator that returns (id, likelihood, confidence, is_background)
# tuples for all sngl_burst <--> sngl_burst coincs.
#


def bb_id_likelihood_confidence_background(database):
	return database.connection.cursor().execute("""
SELECT
	coinc_event.coinc_event_id,
	coinc_event.likelihood,
	multi_burst.confidence,
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
	JOIN multi_burst ON (
		multi_burst.coinc_event_id == coinc_event.coinc_event_id
	)
WHERE
	coinc_event.coinc_def_id == ?
	""", (database.bb_definer_id,))


#
# Construct the temporary sim_coinc_map view that allows for quick mapping
# from sim_burst to sngl_burst<-->sngl_burst coincs
#


def create_sim_coinc_map_view(connection):
	connection.cursor().execute("""
CREATE TEMPORARY VIEW
	sim_coinc_map
AS
	SELECT
		sim_burst.simulation_id AS simulation_id,
		sim_coinc_event.coinc_def_id AS sim_coinc_def_id,
		burst_coinc_event.coinc_event_id AS burst_coinc_event_id,
		burst_coinc_event.likelihood AS burst_coinc_likelihood,
		multi_burst.confidence AS burst_coinc_confidence,
		coinc_detection_statistic(burst_coinc_event.likelihood, multi_burst.confidence) AS burst_coinc_amplitude
	FROM
		sim_burst
		JOIN coinc_event_map AS a ON (
			a.table_name == 'sim_burst'
			AND a.event_id == sim_burst.simulation_id
		)
		JOIN coinc_event AS sim_coinc_event ON (
			sim_coinc_event.coinc_event_id == a.coinc_event_id
		)
		JOIN coinc_event_map AS b ON (
			b.coinc_event_id == a.coinc_event_id
		)
		JOIN coinc_event AS burst_coinc_event ON (
			b.table_name == 'coinc_event'
			AND b.event_id == burst_coinc_event.coinc_event_id
		)
		JOIN multi_burst ON (
			multi_burst.coinc_event_id == burst_coinc_event.coinc_event_id
		)
	""")


#
# =============================================================================
#
#            Measuring Confidence-Likelihood Density Contour Slope
#
# =============================================================================
#


#
# Phase 1:  dump data files.
#


def dump_confidence_likelihood_scatter_data(globs, live_time_program = "lalapps_power", tmp_path = None, verbose = False):
	#
	# Collect file names.
	#

	if verbose:
		print >>sys.stderr, "building file list ..."
	filenames = sorted(filename for g in globs for filename in glob.glob(g))

	#
	# Initialize storage.
	#

	injections = iterutils.Highest(max = 1e6)
	background = iterutils.Highest(max = 1e6)
	zero_lag = iterutils.Highest(max = 1e6)

	#
	# Iterate over files.
	#

	for n, filename in enumerate(filenames):
		#
		# Open the database file.
		#

		if verbose:
			print >>sys.stderr, "%d/%d: %s" % (n + 1, len(filenames), filename)
		working_filename = dbtables.get_connection_filename(filename, tmp_path = tmp_path, verbose = verbose)
		connection = sqlite3.connect(working_filename)
		connection.create_function("coinc_detection_statistic", 2, coinc_detection_statistic)
		database = SnglBurstUtils.CoincDatabase(connection, live_time_program)
		if verbose:
			SnglBurstUtils.summarize_coinc_database(database)

		#
		# Process database contents.  Assume all files with
		# sim_burst tables are the outputs of injection runs, and
		# others aren't.
		#

		if database.sim_burst_table is None:
			# non-injections
			background_buffer = []
			zero_lag_buffer = []
			for id, l, c, is_background in bb_id_likelihood_confidence_background(database):
				if is_background:
					background_buffer.append((coinc_detection_statistic(l, c), l, c))
				else:
					zero_lag_buffer.append((coinc_detection_statistic(l, c), l, c))
			background.extend(background_buffer)
			zero_lag.extend(zero_lag_buffer)
		else:
			# injections
			injections_buffer = []
			create_sim_coinc_map_view(database.connection)
			for a, l, c in database.connection.cursor().execute("""
SELECT
	burst_coinc_amplitude,
	burst_coinc_likelihood,
	burst_coinc_confidence
FROM
	sim_coinc_map
WHERE
	sim_coinc_def_id == ?
			""", (database.sce_definer_id,)):
				injections_buffer.append((-a, l, c))
			injections.extend(injections_buffer)

		#
		# Done with this file.
		#

		database.connection.close()
		dbtables.discard_connection_filename(filename, working_filename, verbose = verbose)

	#
	# Dump scatter plot data.
	#

	if verbose:
		print >>sys.stderr, "writing scatter plot data ..."

	f = file("lalapps_excesspowerfinal_background_scatter.dat", "w")
	for a, l, c in background:
		print >>f, "%.16g %.16g" % (l, c)

	f = file("lalapps_excesspowerfinal_zero_lag_scatter.dat", "w")
	for a, l, c in zero_lag:
		print >>f, "%.16g %.16g" % (l, c)

	f = file("lalapps_excesspowerfinal_injections_scatter.dat", "w")
	for a, l, c in injections:
		print >>f, "%.16g %.16g" % (l, c)

	if verbose:
		print >>sys.stderr, "done."


#
# Phase 2:  generate scatter plot from data files.
#


def plot_confidence_likelihood_scatter_data(slope, verbose = False):
	#
	# Start plot
	#

	fig = SnglBurstUtils.figure.Figure()
	SnglBurstUtils.FigureCanvas(fig)
	fig.set_size_inches(6.5, 6.5 / ((1 + math.sqrt(5)) / 2))
	axes = fig.gca()
	axes.grid(True)
	axes.loglog()
	axes.set_xlabel(r"Confidence")
	axes.set_ylabel(r"Likelihood Ratio")
	axes.set_title(r"Likelihood Ratio vs.\ Confidence Scatter Plot")

	#
	# Plot scatter data
	#

	def read_and_plot(filename, colour, verbose = False):
		if verbose:
			print >>sys.stderr, "reading '%s' ..." % filename
		X = []
		Y = []
		for line in file(filename):
			if line[0] == "#":
				continue
			y, x = map(float, line.strip().split())
			X.append(x)
			Y.append(y)
		if verbose:
			print >>sys.stderr, "plotting ..."
		return axes.plot(X, Y, colour)

	set1 = read_and_plot("lalapps_excesspowerfinal_injections_scatter.dat", "r+", verbose = verbose)
	set2 = read_and_plot("lalapps_excesspowerfinal_background_scatter.dat", "k+", verbose = verbose)
	#set3 = read_and_plot("lalapps_excesspowerfinal_zero_lag_scatter.dat", "b+", verbose = verbose)

	#axes.legend((set1, set2, set3), (r"Injections", r"Background (time slides)", r"Zero lag"), loc = "lower right")
	axes.legend((set1, set2), (r"Injections", r"Background (time slides)"), loc = "lower right")

	#
	# Plot threshold contours
	#

	def fit(x, lnb):
		return numpy.exp(slope * numpy.log(x) + lnb)

	if verbose:
		print >>sys.stderr, "plotting contours ..."
	ymin, ymax = axes.get_ylim()
	for lnb in range(10, 110, 10):
		x = 10**numpy.arange(0.85, 3.0, 0.01)
		y = fit(x, lnb)
		x = numpy.compress((ymin <= y) & (y <= ymax), x)
		y = numpy.compress((ymin <= y) & (y <= ymax), y)
		axes.plot(x, y, "g-")

	#
	# Adjust plot range and add likelihood = 1.0 marker
	#

	axes.plot(axes.get_xlim(), (1, 1), "k-")
	axes.set_xlim(10**0.85, 10**3)

	#
	# Write plot
	#

	if verbose:
		print >>sys.stderr, "writing 'lalapps_excesspowerfinal_scatter.png' ..."
	fig.savefig("lalapps_excesspowerfinal_scatter.png")

	#
	# Done
	#

	if verbose:
		print >>sys.stderr, "done."


#
# =============================================================================
#
#                              Rate vs. Threshold
#
# =============================================================================
#


#
# Book-keeping class
#


class RateVsThresholdData(object):
	def __init__(self):
		self.background_amplitudes = iterutils.Highest(max = 1e7)
		self.zero_lag_amplitudes = []
		self.zero_lag_live_time = 0.0
		self.background_live_time = 0.0


	def update_from(self, contents, filename = "", verbose = False):
		#
		# Add the live times from this file to the totals.
		#

		if verbose:
			print >>sys.stderr, "measuring live time ..."
		zero_lag_time_slides, background_time_slides = SnglBurstUtils.get_time_slides(contents.connection)
		self.zero_lag_live_time += SnglBurstUtils.time_slides_livetime(contents.seglists, zero_lag_time_slides.values(), verbose = verbose)
		self.background_live_time += SnglBurstUtils.time_slides_livetime(contents.seglists, background_time_slides.values(), verbose = verbose)

		#
		# Iterate over burst<-->burst coincidences.  Assume there
		# are no injections in this file.  Accumulate confidences
		# in a separate buffer to avoid slow calls to the Highest
		# class' append() method.
		#

		if verbose:
			print >>sys.stderr, "retrieving sngl_burst<-->sngl_burst coincs ..."
		buffer = []
		for id, likelihood, confidence, is_background in bb_id_likelihood_confidence_background(contents):
			if is_background:
				buffer.append(coinc_detection_statistic(likelihood, confidence))
			else:
				self.zero_lag_amplitudes.append((coinc_detection_statistic(likelihood, confidence), filename, id))
		if verbose:
			print >>sys.stderr, "sorting background coincs ..."
		self.background_amplitudes.extend(buffer)


	def finish(self, zero_lag_survivors, open_box = False, verbose = False):
		# sort the zero lag amplitudes from highest to lowest.

		if not open_box:
			# just forget them for safety
			self.zero_lag_amplitudes = []
		else:
			self.zero_lag_amplitudes.sort(reverse = True)

		# compute the mean event rate vs. amplitude threshold as
		# observed in the time slides.  xcoords are the background
		# event amplitudes, ycoords are the are the corresponding
		# rates.  have to not trick array() with the Highest class'
		# fake element count.  this is \mu_{0}(x) in the Brady et
		# al. paper, but the actual value at threshold is extracted
		# from a smoothed approximation below.  note that because
		# the Highest class's contents are sorted in decreasing
		# order, the x co-ordinates come out in decreasing order,
		# but it's more convenient for them to be in increasing
		# order so we have to reverse the two arrays.

		self.background_rate_x = numpy.fromiter(self.background_amplitudes, dtype = "double", count = list.__len__(self.background_amplitudes))
		self.background_rate_y = numpy.arange(1, len(self.background_rate_x) + 1, dtype = "double") / self.background_live_time
		self.background_rate_x = self.background_rate_x[::-1]
		self.background_rate_y = self.background_rate_y[::-1]

		# compute the 1-sigma vertical error bars.  ycoords is the
		# expected event rate measured from the background, times
		# zero_lag_live_time is the expected event count at zero
		# lag, the square root of which is the 1 sigma Poisson
		# fluctuation in the expected zero lag event count, divided
		# by zero_lag_live_time gives the 1 sigma Poisson
		# fluctuation in the expected zero lag event rate.

		self.background_rate_dy = numpy.sqrt(self.background_rate_y * self.zero_lag_live_time) / self.zero_lag_live_time

		# given the desired number of survivors at zero lag, find
		# the "amplitude" threshold to cut coincs on.  The
		# interpretation is that only coincidences with an
		# "amplitude" greater than (not equal to) this value are to
		# be retained.

		if open_box:
			# obtain threshold from loudest zero-lag trigger
			self.amplitude_threshold = self.zero_lag_amplitudes[zero_lag_survivors][0]
		else:
			# obtain threshold from background rate
			self.amplitude_threshold = self.background_rate_x[-1 - bisect.bisect(self.background_rate_y[::-1], 1.0 / self.zero_lag_live_time)]

		# compute cubic spline approximation of background event
		# rate curve in the vacinity of the amplitude threshold.
		# this is used to obtain a smoothed approximation of the
		# background event rate and its first derivative, which are
		# \mu_{0}(x) and \mu_{0}'(x) in the Brady et al. paper.
		# because of the shape of the curve, it is more natural to
		# approximate ln rate(threshold) instead of
		# rate(threshold).  the curve samples used to construct the
		# spline are chosen by identifying the sample corresponding
		# to the amplitude threshold, and then taking the samples
		# on either side of that corresponding to +/- 1 order of
		# magnitude in rate (a symmetric interval in ln(rate)).  to
		# obtain \mu_{0}'(x) from the spline, if
		#
		#	y = ln rate(threshold)
		#
		# then
		#
		#	dy/dthreshold = 1/rate drate/dthreshold
		#
		# so
		#
		#	drate/dthreshold = rate * dy/dthreshold
		#
		# the uncertainty in \mu_{0}(x) is computed the same way as
		# background_rate_dy above, except using \mu_{0}(x) as the
		# rate and using the total background live time instead of
		# the zero lag live time.

		index = bisect.bisect(self.background_rate_x, self.amplitude_threshold)
		lo = len(self.background_rate_x) - (len(self.background_rate_x) - index) * 10
		hi = len(self.background_rate_x) - (len(self.background_rate_x) - index) / 10
		spline = interpolate.UnivariateSpline(self.background_rate_x[lo:hi], numpy.log(self.background_rate_y[lo:hi]), k = 3)
		self.background_rate_approximant = lambda x: numpy.exp(spline(x))

		self.mu_0 = math.exp(spline(self.amplitude_threshold))
		self.dmu_0 = numpy.sqrt(self.mu_0 * self.background_live_time) / self.background_live_time
		self.mu_0primed = self.mu_0 * spline(self.amplitude_threshold, nu = 1)

		# compute the mean event rate vs. amplitude threshold as
		# observed at zero-lag.  to make the plots nicer, don't
		# include zero lag events whose amplitudes are below the
		# lowest amplitude recorded for a background event.  again,
		# the x co-ordinates come out in decreasing order (the were
		# sorted that way above), but for convenience we reverse
		# the arrays so that the x co-ordinates are in increasing
		# order.

		self.zero_lag_rate_x = numpy.array([zero_lag_amplitude[0] for zero_lag_amplitude in self.zero_lag_amplitudes], dtype = "double")
		self.zero_lag_rate_y = numpy.arange(1, len(self.zero_lag_rate_x) + 1, dtype = "double") / self.zero_lag_live_time
		self.zero_lag_rate_x = self.zero_lag_rate_x[::-1]
		self.zero_lag_rate_y = self.zero_lag_rate_y[::-1]
		index = bisect.bisect(self.zero_lag_rate_x, self.background_rate_x[0])
		self.zero_lag_rate_x = self.zero_lag_rate_x[index:]
		self.zero_lag_rate_y = self.zero_lag_rate_y[index:]


#
# Measure threshold
#


def measure_threshold(filenames, n_survivors, live_time_program = "lalapps_power", tmp_path = None, open_box = False, verbose = False):
	#
	# Initialize the book-keeping object.
	#

	rate_vs_threshold_data = RateVsThresholdData()

	#
	# Iterate over non-injection files.
	#

	for n, filename in enumerate(filenames):
		#
		# Open the database file.
		#

		if verbose:
			print >>sys.stderr, "%d/%d: %s" % (n + 1, len(filenames), filename)
		working_filename = dbtables.get_connection_filename(filename, tmp_path = tmp_path, verbose = verbose)
		database = SnglBurstUtils.CoincDatabase(sqlite3.connect(working_filename), live_time_program)
		if verbose:
			SnglBurstUtils.summarize_coinc_database(database)

		#
		# Process database contents.
		#

		rate_vs_threshold_data.update_from(database, filename = filename, verbose = verbose)

		#
		# Done with this file.
		#

		database.connection.close()
		dbtables.discard_connection_filename(filename, working_filename, verbose = verbose)

	#
	# Determine likelihood threshold.
	#

	if verbose:
		print >>sys.stderr, "finishing rate vs. threshold measurement ..."
	rate_vs_threshold_data.finish(n_survivors, open_box = open_box, verbose = verbose)

	#
	# Done.
	#

	return rate_vs_threshold_data


#
# Print results.
#


def print_rate_vs_threshold_data(rate_vs_threshold_data, confidence_contour_slope):
	print >>sys.stderr
	print >>sys.stderr, "=== Threshold Summary ==="
	print >>sys.stderr
	print >>sys.stderr, "threshold definition:  ln likelihood > %.16g ln confidence + %.16g" % (confidence_contour_slope, rate_vs_threshold_data.amplitude_threshold)
	print >>sys.stderr, "total live time in background = %.16g s" % rate_vs_threshold_data.background_live_time
	print >>sys.stderr, "total live time at zero lag = %.16g s" % rate_vs_threshold_data.zero_lag_live_time
	print >>sys.stderr, "number of coincs in background = %d" % len(rate_vs_threshold_data.background_amplitudes)
	print >>sys.stderr, "average number of background coincs per zero lag live time = %.16g" % (len(rate_vs_threshold_data.background_amplitudes) / rate_vs_threshold_data.background_live_time * rate_vs_threshold_data.zero_lag_live_time)
	print >>sys.stderr, "number of coincs at zero lag = %d" % len(rate_vs_threshold_data.zero_lag_amplitudes)
	print >>sys.stderr, "at threshold, \\mu_{0} = %.16g Hz +/- %.16g Hz" % (rate_vs_threshold_data.mu_0, rate_vs_threshold_data.dmu_0)
	print >>sys.stderr, "at threshold, \\mu_{0}' = %.16g Hz / unit of amplitude" % rate_vs_threshold_data.mu_0primed
	print >>sys.stderr
	print >>sys.stderr, "100 Highest-Ranked Zero Lag Events"
	print >>sys.stderr, "----------------------------------"
	print >>sys.stderr
	print >>sys.stderr, "Detection Statistic\tFilename\tID"
	for amplitude, filename, id in rate_vs_threshold_data.zero_lag_amplitudes[:100]:
		print >>sys.stderr, "%.16g\t%s\t%s" % (amplitude, filename, id)
	print >>sys.stderr


#
# Plot
#


def plot_rate_vs_threshold(data):
	"""
	Input is a RateVsThresholdData instance.  Output is a matplotlib
	Figure instance.
	"""
	# start the rate vs. threshold plot


	print >>sys.stderr
	print >>sys.stderr, "plotting event rate ..."

	fig, axes = SnglBurstUtils.make_burst_plot(r"Detection Statistic Threshold", r"Mean Event Rate (Hz)")
	axes.semilogy()

	# draw the background rate 1-sigma error bars as a shaded polygon
	# clipped to the plot.  warning:  the error bar polygon is not
	# *really* clipped to the axes' bounding box, the result will be
	# incorrect if the number of sample points is small.

	ymin = 1e-9
	ymax = 1e0
	poly_x = numpy.concatenate((data.background_rate_x, data.background_rate_x[::-1]))
	poly_y = numpy.concatenate((data.background_rate_y + 1 * data.background_rate_dy, (data.background_rate_y - 1 * data.background_rate_dy)[::-1]))
	axes.add_patch(patches.Polygon(numpy.column_stack((poly_x, numpy.clip(poly_y, ymin, ymax))), edgecolor = "k", facecolor = "k", alpha = 0.3))

	# draw the mean background event rate

	line1 = axes.plot(data.background_rate_x, data.background_rate_y, "k-")

	# draw the smoothed approximation to the mean background event rate

	#index = bisect.bisect(data.background_rate_x, data.amplitude_threshold)
	#lo = len(data.background_rate_x) - (len(data.background_rate_x) - index) * 10
	#hi = len(data.background_rate_x) - (len(data.background_rate_x) - index) / 10
	#x = data.background_rate_x[lo:hi]
	#line2 = axes.plot(x, data.background_rate_approximant(x), "b-")

	# draw the mean zero-lag event rate.  if the box is closed, these
	# arrays will be empty.

	line3 = axes.plot(data.zero_lag_rate_x, data.zero_lag_rate_y, "r-")

	# draw a vertical marker indicating the threshold

	axes.axvline(data.amplitude_threshold, color = "k")

	#axes.set_ylim((ymin, ymax))
	#axes.xaxis.grid(True, which="minor")
	#axes.yaxis.grid(True, which="minor")

	# add a legend and a title

	axes.legend((line1, line3), (r"Background (time slides)", r"Zero lag"))
	axes.set_title(r"Mean Event Rate vs.\ Detection Statistic Threshold")

	# done rate vs. threshold plot

	#print >>sys.stderr, "writing lalapps_excesspowerfinal_rate_vs_threshold.pdf ..."
	#fig.savefig("lalapps_excesspowerfinal_rate_vs_threshold.pdf")
	print >>sys.stderr, "writing lalapps_excesspowerfinal_rate_vs_threshold.png ..."
	fig.savefig("lalapps_excesspowerfinal_rate_vs_threshold.png")

	# start rate vs. threshold residual plot

	fig, axes = SnglBurstUtils.make_burst_plot(r"Detection Statistic Threshold", r"Mean Event Rate Residual (standard deviations)")

	# construct a linear interpolator for the backround rate data so
	# that we can evaluate it at the x co-ordinates of the zero-lag
	# curve samples

	interp_y = interpolate.interpolate.interp1d(data.background_rate_x, data.background_rate_y, kind = "linear", bounds_error = False)
	interp_dy = interpolate.interpolate.interp1d(data.background_rate_x, data.background_rate_dy, kind = "linear", bounds_error = False)

	# plot the residual

	axes.plot(data.zero_lag_rate_x, (data.zero_lag_rate_y - interp_y(data.zero_lag_rate_x)) / interp_dy(data.zero_lag_rate_x), "k-")

	# add a title

	axes.set_title(r"Mean Event Rate Residual vs.\ Detection Statistic Threshold")

	# done rate vs. threshold residual plot

	#print >>sys.stderr, "writing lalapps_excesspowerfinal_rate_vs_threshold_residual.pdf ..."
	#fig.savefig("lalapps_excesspowerfinal_rate_vs_threshold_residual.pdf")
	print >>sys.stderr, "writing lalapps_excesspowerfinal_rate_vs_threshold_residual.png ..."
	fig.savefig("lalapps_excesspowerfinal_rate_vs_threshold_residual.png")

	# done


#
# =============================================================================
#
#                                  Efficiency
#
# =============================================================================
#


#
# Diagnostic plots
#


def diagnostic_plot(z, bins, title, ylabel, filename):
	print >>sys.stderr, "generating %s ..." % filename
	fig, axes = SnglBurstUtils.make_burst_plot("Centre Frequency (Hz)", ylabel)
	axes.loglog()
	xcoords, ycoords = bins.centres()
	cset = axes.contour(xcoords, ycoords, numpy.transpose(z), 9)
	cset.clabel(inline = True, fontsize = 5, fmt = r"$%g$", colors = "k")
	axes.set_title(title)
	fig.savefig(filename)


#
# Book-keeping class
#


class EfficiencyData(SimBurstUtils.Efficiency_hrss_vs_freq):
	def __init__(self, *args):
		SimBurstUtils.Efficiency_hrss_vs_freq.__init__(self, *args)
		self.recovered_xyz = []


	def add_contents(self, contents, threshold):
		if self.instruments != contents.instruments:
			raise ValueError("this document contains instruments %s, but we want %s" % ("+".join(contents.instruments), "+".join(self.instruments)))

		# iterate over all sims

		offsetvectors = contents.time_slide_table.as_dict()
		create_sim_coinc_map_view(contents.connection)
		for values in contents.connection.cursor().execute("""
SELECT
	sim_burst.*,
	(
		SELECT
			MAX(sim_coinc_map.burst_coinc_amplitude)
		FROM
			sim_coinc_map
		WHERE
			sim_coinc_map.simulation_id == sim_burst.simulation_id
			AND sim_coinc_map.sim_coinc_def_id == ?
	)
FROM
	sim_burst
		""", (contents.scn_definer_id,)):
			sim = contents.sim_burst_table.row_from_cols(values)
			detection_statistic = values[-1]
			amplitude = self.amplitude_func(sim, None)
			if SimBurstUtils.on_instruments(sim, contents.seglists, offsetvectors[sim.time_slide_id]) == self.instruments:
				# injection was made
				self.injected_x.append(sim.frequency)
				self.injected_y.append(amplitude)
				if detection_statistic is not None:
					# injection was recovered (found at
					# all)
					self.recovered_xyz.append((sim.frequency, amplitude, detection_statistic))
					if detection_statistic > threshold:
						# injection was found above
						# threshold
						self.found_x.append(sim.frequency)
						self.found_y.append(amplitude)
			elif (detection_statistic is not None) and (detection_statistic > threshold):
				# injection was found, and found above
				# threshold, but not "injected" because it
				# lies outside of the output segments.
				# this is not unusual, in particular
				# because we are only looking for coincs
				# that are "nearby" an injection which can
				# mean several seconds.
				print >>sys.stderr, "odd, injection %s was found but not injected ..." % sim.simulation_id


	def finish(self, threshold):
		# bin the injections

		self._bin_events()

		# use the same binning for the found injection density as
		# was constructed for the efficiency

		self.found_density = rate.BinnedArray(self.efficiency.denominator.bins)

		# construct the amplitude weighting function

		amplitude_weight = rate.BinnedArray(rate.NDBins((rate.LinearBins(threshold - 100, threshold + 100, 10001),)))

		# gaussian window's width is the number of bins
		# corresponding to 10 units of amplitude, which is computed
		# by dividing 10 by the "volume" of the bin corresponding
		# to threshold.  index is the index of the element in
		# amplitude_weight corresponding to the threshold.

		index, = amplitude_weight.bins[threshold,]
		window = rate.gaussian_window(10.0 / amplitude_weight.bins.volumes()[index])
		window /= 10 * window[(len(window) - 1) / 2]

		# set the window data into the BinnedArray object.  the
		# Gaussian peaks on the middle element of the window, which
		# we want to place at index in the amplitude_weight array.

		lo = index - (len(window) - 1) / 2
		hi = lo + len(window)
		if lo < 0 or hi > len(amplitude_weight.array):
			raise ValueError("amplitude weighting window too large")
		amplitude_weight.array[lo:hi] = window

		# store the recovered injections in the found density bins
		# weighted by amplitude

		for x, y, z in self.recovered_xyz:
			try:
				weight = amplitude_weight[z,]
			except IndexError:
				# beyond the edge of the window
				weight = 0.0
			self.found_density[x, y] += weight

		# the efficiency is only valid up to the highest energy
		# that has been injected.  this creates problems later on
		# so, instead, for each frequency, identify the highest
		# energy that has been measured and copy the values for
		# that bin's numerator and denominator into the bins for
		# all higher energies.  do the same for the counts in the
		# found injection density bins.
		#
		# numpy.indices() returns two arrays array, the first of
		# which has each element set equal to its x index, the
		# second has each element set equal to its y index, we keep
		# the latter.  meanwhile numpy.roll() cyclically permutes
		# the efficiency bins down one along the y (energy) axis.
		# from this, the conditional finds bins where the
		# efficiency is greater than 0.9 but <= 0.9 in the bin
		# immediately above it in energy.  we select the elements
		# from the y index array where the conditional is true, and
		# then use numpy.max() along the y direction to return the
		# highest such y index for each x index, which is a 1-D
		# array.  finally, enumerate() is used to iterate over x
		# index and corresponding y index, and if the y index is
		# not negative (was found) the value from that x-y bin is
		# copied to all higher bins.

		n = self.efficiency.numerator.array
		d = self.efficiency.denominator.array
		f = self.found_density.array
		bady = -1
		for x, y in enumerate(numpy.max(numpy.where((d > 0) & (numpy.roll(d, -1, axis = 1) <= 0), numpy.indices(d.shape)[1], bady), axis = 1)):
			if y != bady:
				n[x, y + 1:] = n[x, y]
				d[x, y + 1:] = d[x, y]
				f[x, y + 1:] = f[x, y]

		# now do the same for the bins at energies below those that
		# have been measured.

		bady = d.shape[1]
		for x, y in enumerate(numpy.min(numpy.where((d > 0) & (numpy.roll(d, 1, axis = 1) <= 0), numpy.indices(d.shape)[1], bady), axis = 1)):
			if y != bady:
				n[x, 0:y] = n[x, y]
				d[x, 0:y] = d[x, y]
				f[x, 0:y] = f[x, y]

		diagnostic_plot(self.efficiency.numerator.array, self.efficiency.denominator.bins, r"Efficiency Numerator (Before Averaging)", self.amplitude_lbl, "lalapps_excesspowerfinal_efficiency_numerator_before.png")
		diagnostic_plot(self.efficiency.denominator.array, self.efficiency.denominator.bins, r"Efficiency Denominator (Before Averaging)", self.amplitude_lbl, "lalapps_excesspowerfinal_efficiency_denominator_before.png")
		diagnostic_plot(self.found_density.array, self.efficiency.denominator.bins, r"Injections Lost / Unit of Threshold (Before Averaging)", self.amplitude_lbl, "lalapps_excesspowerfinal_found_density_before.png")

		# smooth the efficiency bins and the found injection
		# density bins using the same 2D window.

		window = rate.gaussian_window(self.window_size_x, self.window_size_y)
		window /= window[tuple((numpy.array(window.shape, dtype = "double") - 1) / 2)]
		rate.filter_binned_ratios(self.efficiency, window)
		rate.filter_array(self.found_density.array, window)

		diagnostic_plot(self.efficiency.numerator.array, self.efficiency.denominator.bins, r"Efficiency Numerator (After Averaging)", self.amplitude_lbl, "lalapps_excesspowerfinal_efficiency_numerator_after.png")
		diagnostic_plot(self.efficiency.denominator.array, self.efficiency.denominator.bins, r"Efficiency Denominator (After Averaging)", self.amplitude_lbl, "lalapps_excesspowerfinal_efficiency_denominator_after.png")
		diagnostic_plot(self.found_density.array, self.efficiency.denominator.bins, r"Injections Lost / Unit of Threshold (After Averaging)", self.amplitude_lbl, "lalapps_excesspowerfinal_found_density_after.png")

		# compute the uncertainties in the efficiency and its
		# derivative by assuming these to be the binomial counting
		# fluctuations in the numerators.

		p = self.efficiency.ratio()
		self.defficiency = numpy.sqrt(p * (1 - p) / self.efficiency.denominator.array)
		p = self.found_density.array / self.efficiency.denominator.array
		self.dfound_density = numpy.sqrt(p * (1 - p) / self.efficiency.denominator.array)


#
# Measure efficiency
#


def measure_efficiency(filenames, threshold, live_time_program = "lalapps_power", upper_limit_scale = "E", tmp_path = None, verbose = False):
	# FIXME:  instruments are hard-coded.  bad bad bad.  sigh...
	if upper_limit_scale == "E":
		efficiency = EfficiencyData(("H1", "H2", "L1"), (lambda sim, instrument: sim.egw_over_rsquared), r"Equivalent Isotropic Energy ($M_{\odot} / \mathrm{pc}^{2}$)", 0.1)
	elif upper_limit_scale == "hrss":
		efficiency = EfficiencyData(("H1", "H2", "L1"), (lambda sim, instrument: sim.hrss), r"$h_{\mathrm{rss}}$", 0.1)
	else:
		raise ValueError("bad upper_limit_scale %s" % repr(upper_limit_scale))

	#
	# Iterate over injection files.
	#

	for n, filename in enumerate(filenames):
		#
		# Open the database file.
		#

		if verbose:
			print >>sys.stderr, "%d/%d: %s" % (n + 1, len(filenames), filename)
		working_filename = dbtables.get_connection_filename(filename, tmp_path = tmp_path, verbose = verbose)
		connection = sqlite3.connect(working_filename)
		connection.create_function("coinc_detection_statistic", 2, coinc_detection_statistic)
		database = SnglBurstUtils.CoincDatabase(connection, live_time_program)
		if verbose:
			SnglBurstUtils.summarize_coinc_database(database)

		#
		# Process database contents.
		#

		efficiency.add_contents(database, threshold)

		#
		# Done with this file.
		#

		database.connection.close()
		dbtables.discard_connection_filename(filename, working_filename, verbose = verbose)

	#
	# Compute efficiency from the data that have been collected
	#

	if verbose:
		print >>sys.stderr, "binning and smoothnig efficiency data ..."
	efficiency.finish(threshold)

	#
	# Done
	#

	return efficiency


#
# Plot efficiency data
#


def plot_efficiency_data(efficiency_data):
	print >>sys.stderr, "plotting efficiency curves ..."

	# use the stock plotting routing in SimBurstUtils for the
	# efficiency contour plot

	fig = SimBurstUtils.plot_Efficiency_hrss_vs_freq(efficiency_data)
	#fig.gca().set_ylim((3e-17, 3e-10))

	# done

	#print >>sys.stderr, "writing lalapps_excesspowerfinal_efficiency.pdf ..."
	#fig.savefig("lalapps_excesspowerfinal_efficiency.pdf")
	print >>sys.stderr, "writing lalapps_excesspowerfinal_efficiency.png ..."
	fig.savefig("lalapps_excesspowerfinal_efficiency.png")


#
# =============================================================================
#
#                               Rate Upper Limit
#
# =============================================================================
#


#
# Non-linear root finding routine to solve for the rate
#


def mu_p_epsilon(xi, p):
	"""
	Given the background correction factor $\\xi$, and the upper limit
	confidence p, solve

	$p = 1 - \\ee^{-\\mu_{p} \\epsilon} (1 + \\xi \\mu_{p} \\epsilon)$

	for ($\\mu_{p} \\epsilon$), the product of the foreground rate
	upper limit times the efficiency.  Note that this product is always
	well defined by this relation (i.e., even when the efficiency is
	0).
	"""
	def f(mue):
		return math.exp(-mue) * (1 + xi * mue) - (1 - p)
	def fprime(mue):
		return -math.exp(-mue) * (1 + xi * mue - xi)
	return optimize.newton(f, 1, fprime = fprime, tol = 1e-8)


#
# Compute the rate upper limit array
#


class RateData(object):
	def __init__(self, efficiency_data, confidence):
		self.rate_array = rate.BinnedArray(efficiency_data.efficiency.denominator.bins)
		self.amplitude_lbl = efficiency_data.amplitude_lbl
		self.confidence = confidence


def rate_upper_limit(efficiency_data, mu_0primed, zero_lag_live_time, p):
	print >>sys.stderr, "computing rate upper limit ..."

	# initialize the rate upper limit array, giving it the same binning
	# as the efficiency array.

	rate_data = RateData(efficiency_data, p)

	# efficiency

	e = efficiency_data.efficiency.ratio()

	# the found density is the count of injections recovered with an
	# "amplitude" falling within a range of 1 unit centred on the
	# amplitude threshold, and is therefore the number of events lost
	# from the efficiency's numerator as the threshold is increased by
	# 1 unit.  if t is the threshold, f is the number found, m is the
	# number made, and e = f / m is the efficiency, then
 	#
	#	e / (de/dt) = 1 / (d ln e / dt)
	#
	# and
	#
	#	d ln e / dt = d/dt (ln f - ln m)
	#	            = d ln f / dt
	#	            = (1 / f) df/dt
	#
	# so
	#
	#	e / (de/dt) = f / (df/dt)
	#
	# f are the counts in the efficiency's numerator bins, and because
	# the found density's bins are the count of events lost as the
	# threshold is increased they are -df/dt.

	e_over_eprimed = efficiency_data.efficiency.numerator.array / -efficiency_data.found_density.array

	# background rate correction factor, \xi in Brady et al, also a
	# function of frequency and energy.

	xi = 1.0 / (1.0 - e_over_eprimed * abs(mu_0primed * zero_lag_live_time))
	xi = numpy.where(numpy.isnan(e_over_eprimed), 1.0, xi)

	diagnostic_plot(xi, rate_data.rate_array.bins, r"Background Correction Factor $\xi$", rate_data.amplitude_lbl, "lalapps_excesspowerfinal_xi.png")

	# compute the rate upper limit

	for xy in iterutils.MultiIter(*map(xrange, rate_data.rate_array.array.shape)):
		rate_data.rate_array.array[xy] = mu_p_epsilon(xi[xy], p)
	diagnostic_plot(rate_data.rate_array.array, rate_data.rate_array.bins, r"$\mu_{%.2g} \epsilon$" % p, rate_data.amplitude_lbl, "lalapps_excesspowerfinal_mu_p_epsilon.png")
	rate_data.rate_array.array /= zero_lag_live_time * e

	# done

	return rate_data


#
# Plot the rate upper limit array
#


def plot_rate_upper_limit(rate_data):
	print >>sys.stderr, "plotting rate upper limit ..."

	#
	# contour plot in frequency-energy plane
	#

	fig, axes = SnglBurstUtils.make_burst_plot("Centre Frequency (Hz)", rate_data.amplitude_lbl)
	axes.loglog()

	xcoords, ycoords = rate_data.rate_array.centres()
	# zvals = log_{10}(R_{0.9} / 1 Hz)
	zvals = numpy.log10(numpy.clip(rate_data.rate_array.array, 0, 1) / 1.0)
	cset = axes.contour(xcoords, ycoords, numpy.transpose(zvals), 19)
	cset.clabel(inline = True, fontsize = 5, fmt = r"$%g$", colors = "k")

	#axes.set_ylim((3e-17, 3e-10))

	axes.set_title(r"%g\%% Confidence Rate Upper Limit ($\log_{10} R_{%g} / 1\,\mathrm{Hz}$)" % (100 * rate_data.confidence, rate_data.confidence))

	#print >>sys.stderr, "writing lalapps_excesspowerfinal_upper_limit_1.pdf ..."
	#fig.savefig("lalapps_excesspowerfinal_upper_limit_1.pdf")
	print >>sys.stderr, "writing lalapps_excesspowerfinal_upper_limit_1.png ..."
	fig.savefig("lalapps_excesspowerfinal_upper_limit_1.png")

	#
	# rate vs. energy curve at a sample frequency
	#

	fig, axes = SnglBurstUtils.make_burst_plot(rate_data.amplitude_lbl, r"Rate Upper Limit (Hz)")
	axes.loglog()

	zvals = rate_data.rate_array[110, :]
	max = min(zvals) * 1e4
	ycoords = numpy.compress(zvals <= max, ycoords)
	zvals = numpy.compress(zvals <= max, zvals)
	max = min(ycoords) * 1e6
	zvals = numpy.compress(ycoords <= max, zvals)
	ycoords = numpy.compress(ycoords <= max, ycoords)
	axes.plot(ycoords, zvals, "k-")

	axes.set_title(r"%g\%% Confidence Rate Upper Limit at $110\,\mathrm{Hz}$" % (100 * rate_data.confidence))

	#print >>sys.stderr, "writing lalapps_excesspowerfinal_upper_limit_2.pdf ..."
	#fig.savefig("lalapps_excesspowerfinal_upper_limit_2.pdf")
	print >>sys.stderr, "writing lalapps_excesspowerfinal_upper_limit_2.png ..."
	fig.savefig("lalapps_excesspowerfinal_upper_limit_2.png")


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
# The "amplitude" of a coinc, the statistic on which we threshold.  This
# function is defined here, like this, so that the slope of the
# confidence-likelihood isodensity contours doesn't have to be carried
# around.  For example, we can export this function to sqlite3 as a
# two-argument function so that we do not have to pass the slope parameter
# into SQL queries.
#


def coinc_detection_statistic(likelihood, confidence, m = options.confidence_contour_slope):
	# In the 2-D confidence--likelihood parameter space, the background
	# isodensity contours for high confidence, high likelihood,
	# coincident n-tuples are found to be approximated by the family
	# of curves given by
	#
	#	likelihood = b * confidence^m
	#
	# or
	#
	#	ln likelihood = m ln confidence + ln b
	#
	# where m is a constant and b parameterizes the family of curves.
	# Injections are found to have high "b" values, and noise low "b"
	# values.  We use b (actually ln b) as the final, combined,
	# statistic measuring the "goodness" of a coincident n-tuple.
	# Given the likelihood ratio and confidence of a coincident
	# n-tuple, this function computes and returns ln b, which is
	# sometimes referred to as the "amplitude" of a coinc in this code
	# following the language of Brady et al.

	if likelihood <= 0:
		# log() doesn't like 0, so we handle this case separately.
		return NegInf

	if 1.0 / likelihood <= 0:
		# this time likelihood == +inf, which can happen when there
		# are regions of parameter space where no noise is ever,
		# *ever*, seen.
		return PosInf

	return math.log(likelihood) - m * math.log(confidence)


#
# Preparatory work?
#


if options.dump_scatter_data:
	print >>sys.stderr, "=== Confidence-Likelihood Scatter Dump ==="
	print >>sys.stderr
	dump_confidence_likelihood_scatter_data(options.background_glob + options.injections_glob, live_time_program = options.live_time_program, tmp_path = options.tmp_space, verbose = options.verbose)
	print >>sys.stderr
	print >>sys.stderr, "=== Done ==="
	sys.exit(0)


if options.plot_scatter_data:
	print >>sys.stderr, "=== Confidence-Likelihood Scatter Plot ==="
	print >>sys.stderr
	plot_confidence_likelihood_scatter_data(options.confidence_contour_slope, verbose = options.verbose)
	print >>sys.stderr
	print >>sys.stderr, "=== Done ==="
	sys.exit(0)


#
# Accumulate the statistics required to extract rate vs. threshold
# information, and measure the amplitude of the n_survivors+1 loudest
# event.
#


print >>sys.stderr, "=== Threshold ==="
print >>sys.stderr
print >>sys.stderr, "\t\tBOX IS %s" % (options.open_box and "OPEN!!" or "CLOSED!!")
print >>sys.stderr

if options.verbose:
	print >>sys.stderr, "building file list ..."
filenames = sorted(filename for g in options.background_glob for filename in glob.glob(g))
if not filenames:
	raise ValueError("no background/zero lag files found")

if True:
	rate_vs_threshold_data = measure_threshold(filenames, options.zero_lag_survivors, live_time_program = options.live_time_program, tmp_path = options.tmp_space, open_box = options.open_box, verbose = options.verbose)
	print_rate_vs_threshold_data(rate_vs_threshold_data, options.confidence_contour_slope)
	plot_rate_vs_threshold(rate_vs_threshold_data)
else:
	# fake the rate-vs-threshold data to skip straight to injections
	rate_vs_threshold_data = RateVsThresholdData()
	rate_vs_threshold_data.zero_lag_live_time = 1163382.623777472
	rate_vs_threshold_data.amplitude_threshold = 49.3975937091782
	rate_vs_threshold_data.mu_0primed = -2.000347702238217e-07

print >>sys.stderr, "done."
print >>sys.stderr


#
# Efficiency
#


print >>sys.stderr
print >>sys.stderr, "=== Efficiency =="
print >>sys.stderr

if options.verbose:
	print >>sys.stderr, "building file list ..."
filenames = sorted(filename for g in options.injections_glob for filename in glob.glob(g))
if not filenames:
	raise ValueError("no injection files found")

efficiency_data = measure_efficiency(filenames, rate_vs_threshold_data.amplitude_threshold, live_time_program = options.live_time_program, upper_limit_scale = options.upper_limit_scale, tmp_path = options.tmp_space, verbose = options.verbose)

print >>sys.stderr
print >>sys.stderr, "=== Efficiency Summary ==="
print >>sys.stderr

plot_efficiency_data(efficiency_data)

print >>sys.stderr, "done."
print >>sys.stderr


#
# Rate upper limit
#


print >>sys.stderr
print >>sys.stderr, "=== Rate Upper Limit ==="
print >>sys.stderr

rate_data = rate_upper_limit(efficiency_data, rate_vs_threshold_data.mu_0primed, rate_vs_threshold_data.zero_lag_live_time, options.upper_limit_confidence)
plot_rate_upper_limit(rate_data)

print >>sys.stderr, "done."
print >>sys.stderr


#
# Done
#


print >>sys.stderr
print >>sys.stderr, "=== Done ==="
print >>sys.stderr
