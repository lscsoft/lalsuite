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
from optparse import OptionParser
import numpy
import sys


from lal.utils import CacheEntry


from glue import segments
from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils
from glue.ligolw.utils import search_summary as ligolw_search_summary
from lalburst import git_version
from pylal import rate
from lalburst import SnglBurstUtils


lsctables.use_in(ligolw.LIGOLWContentHandler)


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
	parser.add_option("-b", "--base", metavar = "base", default = "plotburst_", help = "set the prefix for output filenames (default = plotburst_)")
	parser.add_option("-f", "--format", metavar = "format", default = "png", help = "set the output image format (default = png)")
	parser.add_option("-i", "--input-cache", metavar = "filename", default = None, help = "get file names from this LAL cache file")
	parser.add_option("--plot", metavar = "number", action = "append", default = None, help = "only generate the given plot number")
	parser.add_option("--frequency-range", metavar = "low,high", default = "0,2200", help = "set the peak frequency range for the plots (default = 0,2200)")
	parser.add_option("--livetime-program", metavar = "name", default = "lalapps_power", help = "set the name of the program whose search_summary rows will define the livetime (default = \"lalapps_power\").")
	parser.add_option("-v", "--verbose", action = "store_true", help = "be verbose")
	options, filenames = parser.parse_args()

	if options.input_cache:
		filenames.extend([c.path for c in map(CacheEntry, file(options.input_cache))])

	if options.plot:
		options.plot = map(int, options.plot)
	else:
		options.plot = range(8)

	options.frequency_range = map(float, options.frequency_range.split(","))

	return options, (filenames or [None])


options, filenames = parse_command_line()


#
# =============================================================================
#
#                                    Input
#
# =============================================================================
#


class Summary(object):
	def __init__(self):
		self.nevents = 0
		self.start_time = []
		self.duration = []
		self.peak_time = []
		self.peak_freq = []
		self.bandwidth = []
		self.lo_freq = []
		self.snr = []
		self.confidence = []


def snglburst_append(self, row):
	global summary
	if row.ifo not in summary:
		summary[row.ifo] = Summary()
	summary[row.ifo].nevents += 1
	summary[row.ifo].start_time.append(lsctables.LIGOTimeGPS(row.start_time, row.start_time_ns))
	summary[row.ifo].duration.append(row.duration)
	summary[row.ifo].peak_time.append(lsctables.LIGOTimeGPS(row.peak_time, row.peak_time_ns))
	summary[row.ifo].peak_freq.append(row.peak_frequency)
	summary[row.ifo].bandwidth.append(row.bandwidth)
	summary[row.ifo].lo_freq.append(row.central_freq - row.bandwidth / 2.0)
	summary[row.ifo].snr.append(row.snr)
	summary[row.ifo].confidence.append(row.confidence)


lsctables.SnglBurstTable.append = snglburst_append


#
# =============================================================================
#
#                             Confidence vs. Time
#
# =============================================================================
#


class ConfidenceVsTime(object):
	def __init__(self, ifo):
		self.fig, self.axes = SnglBurstUtils.make_burst_plot("GPS Time (s)", "Confidence")
		self.ifo = ifo
		self.nevents = 0
		self.x = []
		self.y = []
		self.seglist = segments.segmentlist()
		self.axes.semilogy()

	def add_contents(self, contents, seglists):
		self.nevents += contents[self.ifo].nevents
		self.x.extend(map(float, contents[self.ifo].peak_time))
		self.y.extend(contents[self.ifo].confidence)
		self.seglist |= seglists[self.ifo]

	def finish(self):
		self.axes.set_title("Trigger Confidence vs. Time\n(%d Triggers)" % self.nevents)
		self.axes.plot(self.x, self.y, "k+")
		for seg in ~self.seglist & segments.segmentlist([segments.segment(self.axes.get_xlim())]):
			self.axes.axvspan(float(seg[0]), float(seg[1]), facecolor = "k", alpha = 0.2)


#
# =============================================================================
#
#                        Confidence vs. Peak Frequency
#
# =============================================================================
#


class ConfidenceVsFrequencyScatter(object):
	def __init__(self, ifo):
		self.fig, self.axes = SnglBurstUtils.make_burst_plot("Peak Frequency (Hz)", "Confidence")
		self.ifo = ifo
		self.nevents = 0
		self.x = []
		self.y = []
		self.axes.semilogy()

	def add_contents(self, contents, seglists):
		self.nevents += contents[self.ifo].nevents
		self.x.extend(contents[self.ifo].peak_freq)
		self.y.extend(contents[self.ifo].confidence)

	def finish(self):
		self.axes.set_title("Trigger Confidence vs. Peak Frequency\n(%d Triggers)" % self.nevents)
		self.axes.plot(self.x, self.y, "k+")
		self.axes.set_xlim((min(self.x), max(self.x)))


#
# =============================================================================
#
#                           Rate vs. Peak Frequency
#
# =============================================================================
#


class RateVsPeakFreq(object):
	def __init__(self, ifo, interval, width):
		self.fig, self.axes = SnglBurstUtils.make_burst_plot("Peak Frequency (Hz)", "Trigger Rate Spectral Density (triggers / s / Hz)")
		self.ifo = ifo
		self.nevents = 0
		# 21 bins per filter width
		bins = int(float(abs(interval)) / width) * 21
		binning = rate.NDBins((rate.LinearBins(interval[0], interval[1], bins),))
		self.rate = rate.BinnedArray(binning)

	def add_contents(self, contents, seglists):
		self.nevents += contents[self.ifo].nevents
		for f in contents[self.ifo].peak_freq:
			try:
				self.rate[f,] += 1.0
			except IndexError:
				raise ValueError("trigger peak frequency %g Hz outside plot range [%g Hz, %g Hz]" % (f, self.rate.bins[0].min, self.rate.bins[0].max))

	def finish(self):
		self.axes.set_title("Trigger Rate vs. Peak Frequency\n(%d Triggers)" % self.nevents)
		# 21 bins per filter width
		rate.to_moving_mean_density(self.rate, rate.gaussian_window(21))
		xvals = self.rate.centres()[0]
		self.axes.plot(xvals, self.rate.array, "k")
		self.axes.semilogy()
		self.axes.set_xlim((min(xvals), max(xvals)))


#
# =============================================================================
#
#                          Trigger Duration Histogram
#
# =============================================================================
#


class Durations(object):
	def __init__(self, ifo):
		self.fig, self.axes = SnglBurstUtils.make_burst_plot("Duration (s)", "Trigger Count")
		self.ifo = ifo
		self.nevents = 0
		self.bins = {}

	def add_contents(self, contents, seglists):
		self.nevents += contents[self.ifo].nevents
		for dt in contents[self.ifo].duration:
			if dt not in self.bins:
				self.bins[dt] = 0
			self.bins[dt] += 1

	def finish(self):
		self.axes.set_title("Trigger Durations\n(%d Triggers)" % self.nevents)
		data = self.bins.items()
		data.sort()
		self.axes.plot([d[0] for d in data], [d[1] for d in data], "ko-")


#
# =============================================================================
#
#                       Time Between Triggers Histogram
#
# =============================================================================
#


class Delays(object):
	def __init__(self, ifo, width, max):
		self.fig, self.axes = SnglBurstUtils.make_burst_plot("Delay (s)", "Count / Delay")
		self.ifo = ifo
		self.nevents = 0
		# 21 bins per filter width
		interval = segments.segment(0, max + 2)
		self.bins = rate.BinnedArray(rate.NDBins((rate.LinearBins(interval[0], interval[1], int(float(abs(interval)) / width) * 21),)))
		self.axes.semilogy()

	def add_contents(self, contents, seglists):
		self.nevents += contents[self.ifo].nevents
		peaks = list(contents[self.ifo].peak_time)
		peaks.sort()
		for i in xrange(1, len(peaks)):
			dt = float(peaks[i] - peaks[i - 1])
			try:
				self.bins[dt,] += 1
			except IndexError:
				# out of bounds
				pass

	def finish(self):
		self.axes.set_title("Time Between Triggers\n(%d Triggers)" % self.nevents)

		xvals = self.bins.centres()[0]
		rate.to_moving_mean_density(self.bins, rate.gaussian_window(21))
		self.axes.plot(xvals, self.bins.array, "k")

		self.axes.set_xlim((0, xvals[-1]))
		self.axes.set_ylim((1, 10.0**(int(math.log10(max(self.bins.array))) + 1)))


#
# =============================================================================
#
#                                 Rate vs. SNR
#
# =============================================================================
#


class RateVsSNR(object):
	def __init__(self, ifo):
		self.fig, self.axes = SnglBurstUtils.make_burst_plot("SNR", "Trigger Rate (Hz)")
		self.ifo = ifo
		self.nevents = 0
		self.x = []
		self.seglist = segments.segmentlist()
		self.axes.loglog()

	def add_contents(self, contents, seglists):
		self.nevents += contents[self.ifo].nevents
		self.x.extend(contents[self.ifo].snr)
		self.seglist |= seglists[self.ifo]

	def finish(self):
		self.axes.set_title("Cummulative Trigger Rate vs. SNR\n(%d Triggers)" % self.nevents)
		self.x.sort()
		self.y = numpy.arange(len(self.x), 0.0, -1.0) / float(abs(self.seglist))
		self.axes.plot(self.x, self.y, "ko-")


#
# =============================================================================
#
#                             Rate vs. Confidence
#
# =============================================================================
#


class RateVsConfidence(object):
	def __init__(self, ifo):
		self.fig, self.axes = SnglBurstUtils.make_burst_plot("Confidence", "Trigger Rate (Hz)")
		self.ifo = ifo
		self.nevents = 0
		self.x = []
		self.seglist = segments.segmentlist()
		self.axes.loglog()

	def add_contents(self, contents, seglists):
		self.nevents += contents[self.ifo].nevents
		self.x.extend(contents[self.ifo].confidence)
		self.seglist |= seglists[self.ifo]

	def finish(self):
		self.axes.set_title("Cummulative Trigger Rate vs. Confidence\n(%d Triggers)" % self.nevents)
		self.x.sort()
		self.y = numpy.arange(len(self.x), 0.0, -1.0, "Float64") / float(abs(self.seglist))
		self.axes.plot(self.x, self.y, "ko-")


#
# =============================================================================
#
#                             Time-Frequency Plane
#
# =============================================================================
#


# moved from viz.py in pylal, but looks like it was ripped off of somebody
# else at some point.  maybe rectfill() in matplotlib?
def tfplot(*args, **kwargs):
	"""
	tfplot(x, y, s=20, c='b', marker='o', cmap=None, norm=None,
	vmin=None, vmax=None, alpha=1.0)

	Supported function signatures:

	TFPLOT(x, y)  : make a scatter plot of x vs y

	TFPLOT(x, y, s)  : make a scatter plot of x vs y with size in area
	given by s

	TFPLOT(x, y, s, c) : make a scatter plot of x vs y with size in area
	given by s and colors given by c

	TFPLOT(x, y, s, c, **kwargs) : control colormapping and scaling
	with keyword args; see below

	Make a scatter plot of x versus y.  s is a size in points^2 a scalar
	or an array of the same length as x or y.  c is a color and can be a
	"""
	shading = kwargs.get('shading', 'faceted')
	cmap = kwargs.get('cmap', cm.get_cmap())
	norm = kwargs.get('norm', normalize())
	alpha = kwargs.get('alpha', 1.0)
	vmin = kwargs.get('vmin', None)
	vmax = kwargs.get('vmax', None)
	a = kwargs.get('axes', gca())

	try:
		X, dX, Y, dY, C = args
	except ValueError:
		raise TypeError('Illegal arguments to rectfill; see help(rectfill)')

	Nx, = X.shape
	verts = [(
		(X[i,],        Y[i,]       ),
		(X[i,]+dX[i,], Y[i,]       ),
		(X[i,]+dX[i,], Y[i,]+dY[i,]),
		(X[i,],        Y[i,]+dY[i,])) for i in range(Nx-1)]
	C = array([C[i,] for i in range(Nx-1)])

	if shading == 'faceted':
		edgecolors = (0, 0, 0, 1),
	else:
		edgecolors = 'None'

	collection = PolyCollection(verts, edgecolors = edgecolors, antialiaseds = (0,), linewidths = (0.25,))
	collection.set_alpha(alpha)
	collection.set_array(C)
	if norm is not None:
		assert isinstance(norm, normalize)
	if cmap is not None:
		assert isinstance(cmap, Colormap)
	collection.set_cmap(cmap)
	collection.set_norm(norm)
	if norm is not None:
		collection.set_clim(vmin, vmax)
	minx = amin(X)
	maxx = amax(X)
	miny = amin(Y)
	maxy = amax(Y)
	corners = (minx, miny), (maxx, maxy)
	a.update_datalim( corners )
	a.autoscale_view()
	# add the collection last
	a.add_collection(collection)
	return collection


class TimeFrequencyPlane(object):
	def __init__(self, ifo):
		self.fig, self.axes = SnglBurstUtils.make_burst_plot("GPS Time (s)", "Frequency (Hz)")
		self.ifo = ifo
		self.nevents = 0
		self.seglist = segments.segmentlist()

	def add_contents(self, contents, seglists):
		self.nevents += contents[self.ifo].nevents
		tfplot(numpy.array(map(float, contents[self.ifo].start_time)), numpy.array(contents[self.ifo].duration), numpy.array(contents[self.ifo].lo_freq), numpy.array(contents[self.ifo].bandwidth), numpy.log(numpy.array(contents[self.ifo].confidence)), axes = self.axes)
		self.seglist |= seglists[self.ifo]

	def finish(self):
		self.axes.set_title("Time-Frequency Plane\n(%d Triggers)" % self.nevents)
		for seg in ~self.seglist & segments.segmentlist([segments.segment(self.axes.get_xlim())]):
			self.axes.axvspan(float(seg[0]), float(seg[1]), facecolor = "k", alpha = 0.2)


#
# =============================================================================
#
#                                  Load Data
#
# =============================================================================
#


summary = {}
seglists = segments.segmentlistdict()


for n, filename in enumerate(ligolw_utils.sort_files_by_size(filenames, options.verbose, reverse = True)):
	if options.verbose:
		print >>sys.stderr, "%d/%d:" % (n + 1, len(filenames)),
	xmldoc = ligolw_utils.load_filename(filename, verbose = options.verbose, contenthandler = ligolw.LIGOLWContentHandler)
	seglists |= ligolw_search_summary.segmentlistdict_fromsearchsummary(xmldoc, options.livetime_program).coalesce()
	xmldoc.unlink()


#
# =============================================================================
#
#                                     Plot
#
# =============================================================================
#


def new_plots(ifo, plots):
	l = (
		RateVsPeakFreq(ifo, segments.segment(options.frequency_range), 4),
		Durations(ifo),
		Delays(ifo, 0.25, 20),
		RateVsSNR(ifo),
		RateVsConfidence(ifo),
		ConfidenceVsTime(ifo),
		ConfidenceVsFrequencyScatter(ifo),
		TimeFrequencyPlane(ifo)
	)
	return [l[i] for i in plots]


for ifo in summary.keys():
	format = "%%s%s_%%0%dd.%%s" % (ifo, int(math.log10(max(options.plot) or 1)) + 1)
	for plotnum, plot in zip(options.plot, new_plots(ifo, options.plot)):
		filename = format % (options.base, plotnum, options.format)
		if options.verbose:
			print >>sys.stderr, "adding to %s plot %d ..." % (ifo, plotnum)
		plot.add_contents(summary, seglists)
		if options.verbose:
			print >>sys.stderr, "finishing %s plot %d ..." % (ifo, plotnum)
		plot.finish()
		if options.verbose:
			print >>sys.stderr, "writing %s ..." % filename
		plot.fig.savefig(filename)
