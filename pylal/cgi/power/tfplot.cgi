#!/usr/bin/python

import matplotlib
matplotlib.use("Agg")
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import numarray

from glue import segments
from pylal import viz

import webplot


#
# TF plot description
#

class Plot(webplot.PlotDescription):
	def trig_segment(self):
		return self.segment.protract(2)


#
# How to make a time-frequency plane plot
#

def makeplot(desc, table):
	fig = figure.Figure()
	canvas = FigureCanvasAgg(fig)
	fig.set_figsize_inches(16,8)
	axes = fig.gca()

	bandwidth = table.getColumnByName("bandwidth").asarray()
	lo_freq = table.getColumnByName("central_freq").asarray() - 0.5 * bandwidth
	start_time = numarray.asarray([float(row.get_start()) for row in table])

	if len(table):
		viz.tfplot(start_time, table.getColumnByName("duration").asarray(), lo_freq, bandwidth, numarray.log(-table.getColumnByName("confidence").asarray()), alpha=0.3, axes = axes)

	axes.set_xlim(list(desc.segment), ylim = list(desc.band))

	for seg in ~desc.seglist & segments.segmentlist([desc.segment]):
		axes.axvspan(seg[0], seg[1], facecolor = "k", alpha = 0.2)
	axes.set_title(desc.instrument + " Excess Power Time-Frequency Plane\n(%d Triggers)" % (len(table)))
	axes.set_xlabel("GPS Time (s)")
	axes.set_ylabel("Frequency (Hz)")

	fig.savefig(desc.filename)


#
# Make a plot and send to client
#

description = Plot().parse_form()

makeplot(description, webplot.gettriggers(description)[0])

webplot.SendImage(description)
