#!/usr/bin/python

import matplotlib
matplotlib.use("Agg")
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import numarray

import webplot


#
# Confidence vs. Frequency plot description
#

class Plot(webplot.PlotDescription):
	pass


#
# How to make a confidence vs. frequency plot
#

def makeplot(desc, table):
	fig = figure.Figure()
	canvas = FigureCanvasAgg(fig)
	fig.set_figsize_inches(16,8)
	axes = fig.gca()

	confidence = -table.getColumnByName("confidence").asarray()
	central_freq = table.getColumnByName("central_freq").asarray()

	axes.semilogy(central_freq, confidence, "k+")

	axes.set_xlim(list(desc.band))
	axes.set_xticks(numarray.arange(desc.band[0], desc.band[1], 100))
	axes.grid(True)

	axes.set_title(desc.instrument + " Excess Power Trigger Confidence vs. Central Frequency\n(GPS Times %s ... %s, %d Triggers)" % (desc.segment[0], desc.segment[1], len(table)))
	axes.set_xlabel("Central Frequency (Hz)")
	axes.set_ylabel("|Confidence|")

	fig.savefig(desc.filename)


#
# Make a plot and send to client
#

description = Plot().parse_form()

makeplot(description, webplot.gettriggers(description)[0])

webplot.SendImage(description)
