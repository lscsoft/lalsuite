#!/usr/bin/python

import matplotlib
matplotlib.use("Agg")
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import numarray

from glue import segments

import webplot


#
# Confidence vs. Time plot description
#

class Plot(webplot.PlotDescription):
	pass


#
# How to make a confidence vs. time plot
#

def makeplot(desc, table):
	fig = figure.Figure()
	canvas = FigureCanvasAgg(fig)
	fig.set_figsize_inches(16,8)
	axes = fig.gca()

	confidence = -table.getColumnByName("confidence").asarray()
	peak_time = numarray.asarray([float(row.get_peak()) for row in table])

	axes.semilogy(peak_time, confidence, "k+")

	for seg in ~desc.seglist & segments.segmentlist([desc.segment]):
		axes.axvspan(seg[0], seg[1], facecolor = "k", alpha = 0.2)

	axes.set_xlim(list(desc.segment))
	axes.grid(True)

	axes.set_title(desc.instrument + " Excess Power Trigger Confidence vs. Time\n(%d Triggers)" % (len(table)))
	axes.set_xlabel("GPS Time (s)")
	axes.set_ylabel("|Confidence|")

	fig.savefig(desc.filename)


#
# Make a plot and send to client
#

description = Plot().parse_form()

makeplot(description, webplot.gettriggers(description)[0])

webplot.SendImage(description)
