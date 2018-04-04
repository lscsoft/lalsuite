#!/usr/bin/python

import math
import matplotlib
matplotlib.use("Agg")
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import numarray
from numarray import nd_image

from glue import segments
from pylal import rate
from pylal.date import LIGOTimeGPS

import webplot


#
# Confidence rate ratio plot description
#

class Plot(webplot.PlotDescription):
	def __init__(self):
		webplot.PlotDescription.__init__(self)
		self.widthratio = 3.0
		self.glitchthreshold = math.sqrt(self.widthratio)

	def trig_segment(self):
		s = self.segment.protract(5.0 * self.widthratio * self.ratewidth)
		# hack to work around problems in numarray
		return segments.segment(float(s[0]), float(s[1]))


#
# How to make a confidence rate ratio plot
#

def glitchsegments(xvals, yvals, threshold):
	# find starts and ends of segments above threshold
	if yvals[0] >= threshold:
		starts = [LIGOTimeGPS(xvals[0])]
	else:
		starts = []
	ends = []
	for i in range(0, len(yvals) - 1):
		if (yvals[i] < threshold) and (yvals[i + 1] >= threshold):
			starts.append(LIGOTimeGPS(xvals[i]))
	for i in range(1, len(yvals)):
		if (yvals[i] < threshold) and (yvals[i - 1] >= threshold):
			ends.append(LIGOTimeGPS(xvals[i]))
	if yvals[-1] >= threshold:
		ends.append(LIGOTimeGPS(xvals[-1]))

	# turn start/end pairs into segments
	return segments.segmentlist(map(segments.segment, starts, ends)).coalesce()


def makeplot(desc, table):
	fig = figure.Figure()
	canvas = FigureCanvasAgg(fig)
	fig.set_figsize_inches(16,8)
	axes = fig.gca()

	# extract peak times and confidences
	peaktime = [float(row.get_peak()) for row in table]
	confidence = numarray.log(-table.getColumnByName("confidence").asarray())

	# construct short time scale average confidence rate, and a long
	# time scale average confidence rate.
	bins = rate.Rate(desc.trig_segment(), desc.ratewidth)
	bins2 = rate.Rate(desc.trig_segment(), desc.widthratio * desc.ratewidth)
	for i in xrange(len(peaktime)):
		try:
			bins[peaktime[i]] += confidence[i]
			bins2[peaktime[i]] += confidence[i]
		except IndexError:
			# trigger lies outside bounds of plot
			pass
	yvals = bins.array = bins.filter()
	yvals2 = bins2.array = bins2.filter()

	# resample to match sample rate of short time scale.
	bins2.array = nd_image.zoom(bins2.array, float(len(yvals)) / float(len(yvals2)))

	# compute ratio, setting 0/0 equal to 0
	yvals = numarray.where(yvals2 > 0.0, yvals, 0.0) / numarray.where(yvals2 > 0.0, yvals2, 1.0)

	# determine segments where ratio is above threshold
	glitchsegs = glitchsegments(bins.xvals(), yvals, desc.glitchthreshold)

	# subtract segments near boundaries of data
	glitchsegs -= (~desc.seglist).protract(desc.widthratio * desc.ratewidth)

	# plot ratio vs time
	axes.plot(bins.xvals(), yvals, "k")

	# tinker with graph
	axes.axhline(desc.glitchthreshold, color = "r")
	axes.set_xlim(list(desc.segment))
	axes.set_ylim([0, desc.widthratio])
	axes.grid(True)

	for seg in ~desc.seglist & segments.segmentlist([desc.segment]):
		axes.axvspan(seg[0], seg[1], facecolor = "k", alpha = 0.2)
	for seg in glitchsegs & segments.segmentlist([desc.segment]):
		axes.axvspan(seg[0], seg[1], facecolor = "r", alpha = 0.2)

	axes.set_title(desc.instrument + " Excess Power %g s Confidence Rate to %g s Confidence Rate Ratio vs. Time\n(%d Triggers)" % (desc.ratewidth, desc.widthratio * desc.ratewidth, len(table)))
	axes.set_xlabel("GPS Time (s)")
	axes.set_ylabel("Confidence Rate Ratio")

	fig.savefig(desc.filename)


#
# Make a plot and send to client
#

description = Plot().parse_form()

makeplot(description, webplot.gettriggers(description)[0])

webplot.SendImage(description)
