#!/usr/bin/python

import matplotlib
matplotlib.use("Agg")
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import numarray
import os
import shutil
import sys
import tempfile
import time

from glue import segments
from pylal.date import LIGOTimeGPS, XLALUTCToGPS, XLALGPSToUTC, UTCMidnights

import eventdisplay


#
# Description
#

handle, filename = tempfile.mkstemp("." + "png", "webplot_")
os.close(handle)

now = XLALUTCToGPS(time.gmtime())
s5length = 365.25 * 24.0 * 60.0 * 60.0	# 1 year

trigsegs = eventdisplay.TrigSegs()
seglist = trigsegs.H1 & trigsegs.H2 & trigsegs.L1
del trigsegs


def end_date(t):
	l = seglist - segments.segmentlist([segments.segment(t, segments.infinity())])
	livetime = float(abs(l))
	if livetime > 0.0:
		return now + (s5length/livetime - 1) * (l[-1][1] - l[0][0])
	else:
		return now


#
# How to make a plot of projected end date vs. time.
#

def make_xticks(a, b):
	locs = list(UTCMidnights(a, b))
	labels = []
	for tm in map(lambda t: time.struct_time(XLALGPSToUTC(t)), locs):
		if tm.tm_wday == 1:     # tuesday
			labels.append(time.strftime("%H h, %a %b %d, %Y", tm))
		else:
			labels.append("")
	return map(float, locs), labels

def make_yticks(a, b):
	tms = map(lambda t: time.struct_time(XLALGPSToUTC(t)), list(UTCMidnights(a, b)))
	tms = filter(lambda tm: tm.tm_mday == 1, tms)
	labels = map(lambda tm: time.strftime("%H h, %a %b %d, %Y", tm), tms)
	locs = map(lambda tm: float(XLALUTCToGPS(tm)), tms)
	return locs, labels


def makeplot():
	fig = figure.Figure()
	canvas = FigureCanvasAgg(fig)
	fig.set_figsize_inches(16,8)
	axes = fig.gca()
	increment = 12 * 3600.0

	xvals = numarray.arange(float(seglist[0][0]) + increment, float(now), increment)
	yvals = map(end_date, xvals)

	axes.plot(map(float, xvals), map(float, yvals), "k")

	axes.set_xlim([seglist[0][0], now])
	axes.set_ylim([851644814, 946339214])
	axes.grid(True)

	ticks = make_xticks(seglist[0][0], now)
	axes.set_xticks(ticks[0])
	axes.set_xticklabels(ticks[1], horizontalalignment="right", fontsize=9, rotation=10)
	ticks = make_yticks(LIGOTimeGPS(851644814), LIGOTimeGPS(946339214))
	axes.set_yticks(ticks[0])
	axes.set_yticklabels(ticks[1], horizontalalignment="right", fontsize=9, rotation=10)

	axes.set_title("Projected S5 End Time vs. Time")
	axes.set_xlabel("Time (UTC)")
	axes.set_ylabel("Projected S5 End Time (UTC)")

	fig.savefig(filename)


#
# Make a plot and send to client
#

makeplot()

#webplot.SendImage(filename)
print >>sys.stdout, "Content-Type: image/png\n"
shutil.copyfileobj(file(filename), sys.stdout)
os.remove(filename)
