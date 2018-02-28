#!/usr/bin/python

import matplotlib
matplotlib.use("Agg")
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import numarray

from glue import segments

import webplot


#
# Cummulative rate vs. confidence plot description
#

class Plot(webplot.PlotDescription):
	pass


#
# How to make a cummulative rate vs. confidence plot
#

def makeplot(desc, table):
	duration = float(abs(desc.seglist & segments.segmentlist([desc.segment])))
	
	fig = figure.Figure()
	canvas = FigureCanvasAgg(fig)
	fig.set_figsize_inches(16,8)
	axes = fig.gca()

	confidence = numarray.sort(-table.getColumnByName("confidence").asarray())
	yvals = numarray.arrayrange(len(confidence), 0.0, -1.0) / duration

	axes.loglog(confidence, yvals, "ko-")
	axes.grid(True)

	axes.set_title(desc.instrument + " Excess Power Cummulative Trigger Rate vs. Confidence\n(GPS Times %s ... %s, %d Triggers)" % (desc.segment[0], desc.segment[1], len(table)))
	axes.set_xlabel("|Confidence|")
	axes.set_ylabel("Rate (Hz)")

	fig.savefig(desc.filename)


#
# Make a plot and send to client
#

description = Plot().parse_form()

makeplot(description, webplot.gettriggers(description)[0])

webplot.SendImage(description)
