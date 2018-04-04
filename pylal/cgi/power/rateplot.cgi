#!/usr/bin/python

import matplotlib
matplotlib.use("Agg")
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

from glue import segments
from pylal import rate

import webplot


#
# Rate plot description
#

class Plot(webplot.PlotDescription):
	def trig_segment(self):
		return self.segment.protract(5.0 * self.ratewidth)


#
# How to make a rate plot
#

def makeplot(desc, table):
	fig = figure.Figure()
	canvas = FigureCanvasAgg(fig)
	fig.set_figsize_inches(16,8)
	axes = fig.gca()

	bins = rate.Rate(segments.segment(float(desc.trig_segment()[0]), float(desc.trig_segment()[1])), desc.ratewidth)
	for row in table:
		try:
			bins[float(row.get_peak())] += 1
		except IndexError:
			# trigger lies outside the bounds of the plot
			pass

	axes.plot(bins.xvals(), bins.filter(), "k")

	axes.set_xlim(list(desc.segment))
	axes.grid(True)

	for seg in ~desc.seglist & segments.segmentlist([desc.segment]):
		axes.axvspan(seg[0], seg[1], facecolor = "k", alpha = 0.2)

	axes.set_title(desc.instrument + " Excess Power Trigger Rate vs. Time\n(%d Triggers, %g s Average)" % (len(table), desc.ratewidth))
	axes.set_xlabel("GPS Time (s)")
	axes.set_ylabel("Trigger Rate (Hz)")

	fig.savefig(desc.filename)


#
# Make a plot and send to client
#

description = Plot().parse_form()

makeplot(description, webplot.gettriggers(description)[0])

webplot.SendImage(description)
