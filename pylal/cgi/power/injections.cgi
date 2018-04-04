#!/usr/bin/python

import matplotlib
matplotlib.use("Agg")
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import numarray

from glue import segments

import webplot


#
# Injection plot description
#

class Plot(webplot.PlotDescription):
	def trig_segment(self):
		return self.segment.protract(2)


#
# How to make a time-frequency plane plot of injections
#


def makeplot(desc, table):
	fig = figure.Figure()
	canvas = FigureCanvasAgg(fig)
	fig.set_figsize_inches(16,8)
	axes = fig.gca()

	if desc.instrument[:2] in ["H1", "H2"]:
		time = numarray.asarray([float(row.get_h_peak()) for row in table])
	elif desc.instrument[:2] in ["L1"]:
		time = numarray.asarray([float(row.get_l_peak()) for row in table])

	freq = table.getColumnByName("freq").asarray()

	axes.plot(time, freq, "k+")

	for seg in ~desc.seglist & segments.segmentlist([desc.segment]):
		axes.axvspan(seg[0], seg[1], facecolor = "k", alpha = 0.2)

	axes.set_xlim(list(desc.segment), ylim = list(desc.band))
	axes.grid(True)

	axes.set_title(desc.instrument + " Injection Locations\n(%d Injections)" % (len(table)))
	axes.set_xlabel("GPS Time (s)")
	axes.set_ylabel("Frequency (Hz)")

	fig.savefig(desc.filename)


#
# Make a plot and send to client
#

description = Plot().parse_form()

makeplot(description, webplot.gettriggers(description)[1])

webplot.SendImage(description)
