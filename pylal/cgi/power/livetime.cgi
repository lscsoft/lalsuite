#!/usr/bin/python

import time

from pylal.date import XLALGPSToUTC, XLALUTCToGPS

import eventdisplay


#
# Live time info
#

def live_time_row_markup(inst, seglist):
	livetime = float(abs(seglist))
	if livetime > 0.0:
		rate = livetime / float(abs(seglist.extent()))
	else:
		rate = 0.0
	if len(seglist):
		averageseg = livetime / len(seglist)
		durations = [float(abs(seg)) for seg in seglist]
		longestseg = max(durations)
		shortestseg = min(durations)
		del durations
	else:
		averageseg = 0.0
		longestseg = 0.0
		shortestseg = 0.0
	s = """<tr>\n"""
	s += """	<td align="left">%s</td>\n""" % inst
	s += """	<td align="right">%.2f h</td>\n""" % (livetime / 60.0 / 60.0,)
	s += """	<td align="center">%.3f</td>\n""" % (rate,)
	s += """	<td align="center">%d</td>\n""" % (len(seglist),)
	s += """	<td align="right">%.2f h</td>\n""" % (averageseg / 60.0 / 60.0,)
	s += """	<td align="right">%.2f h</td>\n""" % (longestseg / 60.0 / 60.0,)
	s += """	<td align="right">%.2f s</td>\n""" % (shortestseg,)
	s += """</tr>\n"""
	return s

def live_time_markup(trigsegs):
	G1H1 = trigsegs.G1 & trigsegs.H1
	G1H2 = trigsegs.G1 & trigsegs.H2
	G1L1 = trigsegs.G1 & trigsegs.L1
	H1H2 = trigsegs.H1 & trigsegs.H2
	H2L1 = trigsegs.H2 & trigsegs.L1
	s = """<table frame="box" rules="all">\n"""
	s += """<thead><tr>\n"""
	s += """	<th>Instruments</th>\n"""
	s += """	<th>Live Time</th>\n"""
	s += """	<th>Rate</th>\n"""
	s += """	<th>Number of<br>Segments</th>\n"""
	s += """	<th>Average<br>Length</th>\n"""
	s += """	<th>Longest</th>\n"""
	s += """	<th>Shortest</th>\n"""
	s += """</tr></thead>\n"""
	s += """<tbody>\n"""
	s += live_time_row_markup("G1", trigsegs.G1)
	s += live_time_row_markup("H1", trigsegs.H1)
	s += live_time_row_markup("H2", trigsegs.H2)
	s += live_time_row_markup("L1", trigsegs.L1)
	s += live_time_row_markup("G1 &cap; H1", G1H1)
	s += live_time_row_markup("G1 &cap; H2", G1H2)
	s += live_time_row_markup("G1 &cap; L1", G1L1)
	s += live_time_row_markup("H1 &cap; H2", H1H2)
	s += live_time_row_markup("H1 &cap; L1", trigsegs.H1 & trigsegs.L1)
	s += live_time_row_markup("H2 &cap; L1", H2L1)
	s += live_time_row_markup("G1 &cap; H1 &cap; H2", G1H1 & trigsegs.H2)
	s += live_time_row_markup("G1 &cap; H1 &cap; L1", G1H1 & trigsegs.L1)
	s += live_time_row_markup("G1 &cap; H2 &cap; L1", G1H2 & trigsegs.L1)
	s += live_time_row_markup("H1 &cap; H2 &cap; L1", H1H2 & trigsegs.L1)
	s += live_time_row_markup("G1 &cap; H1 &cap; H2 &cap; L1", G1H1 & H2L1)
	s += """</tbody>\n"""
	s += """</table>"""
	return s

def s5_live_time_summary(now, seglist):
	s5length = 1.0 * 365.25 * 24.0 * 60.0 * 60.0	# 1 year
	livetime = float(abs(seglist))
	if livetime > 0.0:
		rate = livetime / float(abs(seglist.extent()))
	else:
		rate = 0.0
	s5end = now + (s5length - livetime) / rate
	s5end.nanoseconds = 0	# convert to UTC requires integer seconds
	s = """<table>\n"""
	s += """<tr>\n"""
	s += """	<td>H1 &cap; H2 &cap; L1 hours required for S5</td>\n"""
	s += """	<td>%.2f h</td>\n""" % (s5length / 60.0 / 60.0,)
	s += """</tr>\n"""
	s += """<tr>\n"""
	s += """	<td>Fraction of S5 completed</td>\n"""
	s += """	<td>%.1f%%</td>\n""" % (100.0 * livetime / s5length,)
	s += """</tr>\n"""
	s += """<tr>\n"""
	s += """	<td>Estimated S5 termination (<a href="enddateplot.cgi">history</a>)</td>\n"""
	s += """	<td>%s UTC (GPS %d s)</td>\n""" % (time.asctime(XLALGPSToUTC(s5end)), s5end)
	s += """</tr>\n"""
	s += """</table>"""
	return s


#
# Generate output.
#

trigsegs = eventdisplay.TrigSegs()

print "Content-Type: text/html\n"

print "<html>"
print "<p><center>"
print live_time_markup(trigsegs)
print "</center></p>"
print "<p><center>"
print s5_live_time_summary(XLALUTCToGPS(time.gmtime()), trigsegs.H1 & trigsegs.H2 & trigsegs.L1)
print "</center></p>"
print "</html>"
