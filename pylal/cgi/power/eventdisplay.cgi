#!/usr/bin/python

import cgi
import cgitb ; cgitb.enable()
import time

from glue import segments
from pylal.date import XLALUTCToGPS, LIGOTimeGPS

#
# Init
#

now = XLALUTCToGPS(time.gmtime())

default_start = "now"
default_duration = str(-1 * 3600)


#
# Parse display request
#

class Query(object):
	segment = None
	ratewidth = None
	freqwidth = None
	band = None
	instrument = None
	cluster = None

form = cgi.FieldStorage()

query = Query()

start = form.getfirst("start", default_start).lower()
duration = LIGOTimeGPS(form.getfirst("dur", default_duration))
if start == "now":
	query.segment = segments.segment(now, now + duration)
	refresh = """<meta http-equiv="refresh" content="%d"></meta>""" % (abs(duration) / 3600 * 180 + 180)
else:
	query.segment = segments.segment(LIGOTimeGPS(start), LIGOTimeGPS(start) + duration)
	refresh = ""

query.ratewidth = float(form.getfirst("ratewidth", "60"))
query.freqwidth = float(form.getfirst("freqwidth", "16"))
query.band = segments.segment(float(form.getfirst("lofreq", "0")), float(form.getfirst("hifreq", "2500")))
query.instrument = form.getfirst("inst", "H1")
query.cluster = int(form.getfirst("cluster", "0"))


#
# An error message
#

def errormsg(msg):
	return """<center>Error: %s</center>""" % msg


#
# Full-run trigger rate plots
#

def fullrunrates():
	s = """<table width="100%">"""
	s += "<tr>"
	s += """	<th>H1</th>"""
	s += """	<td width="90%"><object data="http://www.lsc-group.phys.uwm.edu/~kipp/S5/static/H1.png" type="image/png" width="100%"></td>"""
	s += """	<td align="center"><a href="http://www.lsc-group.phys.uwm.edu/~kipp/S5/static/H1.png">PNG</a><br><a href="http://www.lsc-group.phys.uwm.edu/~kipp/S5/static/H1.svg">SVG</a><br><a href="http://www.lsc-group.phys.uwm.edu/~kipp/S5/static/H1.eps">EPS</a></td>"""
	s += "</tr>"
	s += "<tr>"
	s += """	<th>H2</th>"""
	s += """	<td width="90%"><object data="http://www.lsc-group.phys.uwm.edu/~kipp/S5/static/H2.png" type="image/png" width="100%"></td>"""
	s += """	<td align="center"><a href="http://www.lsc-group.phys.uwm.edu/~kipp/S5/static/H2.png">PNG</a><br><a href="http://www.lsc-group.phys.uwm.edu/~kipp/S5/static/H2.svg">SVG</a><br><a href="http://www.lsc-group.phys.uwm.edu/~kipp/S5/static/H2.eps">EPS</a></td>"""
	s += "</tr>"
	s += "<tr>"
	s += """	<th>L1</th>"""
	s += """	<td width="90%"><object data="http://www.lsc-group.phys.uwm.edu/~kipp/S5/static/L1.png" type="image/png" width="100%"></td>"""
	s += """	<td align="center"><a href="http://www.lsc-group.phys.uwm.edu/~kipp/S5/static/L1.png">PNG</a><br><a href="http://www.lsc-group.phys.uwm.edu/~kipp/S5/static/L1.svg">SVG</a><br><a href="http://www.lsc-group.phys.uwm.edu/~kipp/S5/static/L1.eps">EPS</a></td>"""
	s += "</tr>"
	s += "</table>"
	return s


#
# Plot markup
#

def _imgsrc(name, query):
	s = "%s?inst=%s&start=%s&dur=%s&ratewidth=%s&freqwidth=%s&lofreq=%s&hifreq=%s" % (name, query.instrument, query.segment[0], abs(query.segment), query.ratewidth, query.freqwidth, query.band[0], query.band[1])
	if query.cluster:
		s += "&cluster=1"
	return s

def plot_pnglink(name, query):
	src = _imgsrc(name, query) + "&format=png"
	return """<a href="%s">PNG</a>""" % src

def plot_pngthumbnail(name, query):
	src = _imgsrc(name, query) + "&format=png"
	return """<a href="%s"><object data="%s" type="image/png" width="800" standby="Generating image.  Please wait.">Failure loading image (timeout?).</object></a>""" % (src, src)

def plot_epslink(name, query):
	src = _imgsrc(name, query) + "&format=eps"
	return """<a href="%s">EPS</a>""" % src

def plot_svglink(name, query):
	src = _imgsrc(name, query) + "&format=svg"
	return """<a href="%s">SVG</a>""" % src

def plot_table_row(name, query):
	s = "<tr>\n"
	s += """\t<td align="center">""" + plot_pngthumbnail(name, query) + "</td>\n"
	s += """\t<td align="center">""" + plot_pnglink(name, query) + "<br>" + plot_svglink(name, query) + "<br>" + plot_epslink(name, query) + "</td>\n"
	s += "</tr>"
	return s


#
# Trigger download markup
#

def triggerlink(name, query):
	src = _imgsrc(name, query) + "&format=xml"
	return """<a href="%s">Download These Triggers</a>""" % src


#
# Form markup
#

def formmarkup(query):
	def instrumentlist(default):
		s = """<select name="inst"><option>""" + default + """</option>"""
		for inst in [inst for inst in ["G1", "G1 Injections", "H1", "H1 Injections", "H2", "H2 Injections", "L1", "L1 Injections"] if inst != default]:
			s += "<option>" + inst + "</option>"
		return s + "</select>"

	s = """<form action="eventdisplay.cgi" method="get">\n"""
	s += """<table>\n"""
	s += """<tr>\n"""
	s += """	<td><label for="inst">Instrument:</label></td>\n"""
	s += """	<td>""" + instrumentlist(query.instrument) + """</td>\n"""
	s += """</tr>\n"""
	s += """<tr>\n"""
	s += """	<td>GPS Time Range:</td>\n"""
	s += """	<td><label for="start">start=</label><input type="text" name="start" value=\"""" + form.getfirst("start", start) + """\"> s, <label for="dur">duration=</label><input type="text" name="dur" value=\"""" + form.getfirst("dur", "") + """\"> s</td>\n"""
	s += """</tr>\n"""
	s += """<tr>\n"""
	s += """	<td><label for="ratewidth">Triggers per Second Window:</label></td>\n"""
	s += """	<td><input type="text" name="ratewidth" value=\"""" + str(query.ratewidth) + """\"> s</td>\n"""
	s += """</tr>\n"""
	s += """<tr>\n"""
	s += """	<td><label for="freqwidth">Triggers per Hz Window:</label></td>\n"""
	s += """	<td><input type="text" name="freqwidth" value=\"""" + str(query.freqwidth) + """\"> Hz</td>\n"""
	s += """</tr>\n"""
	s += """<tr>\n"""
	s += """	<td><label for="lofreq">Frequency Band:</label></td>\n"""
	s += """	<td><input type="text" name="lofreq" value=\"""" + str(query.band[0]) + """\"> Hz to <input type="text" name="hifreq" value=\"""" + str(query.band[1]) + """\"> Hz</td>\n"""
	s += """</tr>\n"""
	s += """<tr>\n"""
	s += """	<td><label for="cluster">(Re)Cluster Triggers:</label></td>\n"""
	if query.cluster:
		s += """	<td><input type="checkbox" name="cluster" value="1" checked></td>\n"""
	else:
		s += """	<td><input type="checkbox" name="cluster" value="1"></td>\n"""
	s += """</tr>\n"""
	s += """</table>\n"""
	s += """<center><input type="submit" value="Submit"></center>
</form>"""
	return s


def formnotes():
	s = "Notes:\n"
	s += "<ul>\n"
	s += "	<li>The triggers were clustered by the job that created them, but no clustering has been done at the boundaries between analysis jobs.&nbsp; If the (Re)Cluster option is selected above, then the triggers are clustered before generating the plots.&nbsp; This removes duplicates at boundaries between analysis jobs, but is very slow because each plot clusters the triggers itself.</li>\n"
	s += "	<li>The frequency band selected only affects the limits on the axes of the plots, it does not limit the plots to those triggers alone.&nbsp; For example, the trigger rate plot shows the rate for all triggers, not only the ones lieing in the requested frequency band.</li>\n"
	s += "	<li>An explanation of the confidence rate ratio plot.&nbsp; Each trigger has a confidence assigned to it, and the average confidence rate in some interval is the number of triggers seen in that interval times their average confidence divided by the length of the interval.&nbsp; The confidence rate ratio plot computes this quantity over two time scales, and plots the ratio of the short time scale rate to the long time scale rate &mdash; how many times faster confidence is accumulating in the short time scale compared to the long time scale.&nbsp; A momentary cluster of many high confidence triggers causes a larger fluctuation in the short time average than in the long time average, and so is seen as an increase in their ratio.&nbsp; The plot identifies regions of the graph where the ratio excedes a threshold, which is what we would call a glitch.</li>\n"
	s += "</ul>"
	return s


#
# Excess Power Event Interface
#

print "Content-Type: text/html\n"

print """<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">"""
print """<html>"""

print """<head>"""
print refresh
print """</head>"""

print """<body>"""
print "<center><h1>Excess Power Event Interface of Patience</h1></center>"
print "<center>(Patience Required)</center>"
print """<p>You can <a href="http://www.lsc-group.phys.uwm.edu/cgi-bin/cvs/viewcvs.cgi/pylal/cgi/power/?cvsroot=lscsoft">browse the source code for these web pages</a> and see what takes so long.</p>"""

print """<hr width="90%">"""
print "<center>"
print "<h2>S5 Excess Power Live Times To Date</h2>"
print "</center>"
print """<object data="livetime.cgi" type="text/html" width="100%" standby="Computing live time chart.  Please wait.">Failure loading live time chart (timeout?).</object>"""
print """<hr width="90%">"""
print "<center>"
print "<h2>Excess Power Trigger Rates for S5</h2>"
print "</center>"
print "<p><center>"
print fullrunrates()
print "</center></p>"
print """<hr width="90%">"""

if abs(query.segment) > 24 * 3600:
	# Time interval too long error
	print errormsg("Requested segment is too long (24 hour max)")
else:
	# Table of plots
	print "<center>"
	print "<h2>%s s Starting At %s</h2>" % (duration, start.title())
	print "</center>"
	print "<p>"
	print formmarkup(query)
	print "</p>"
	if True:
		print "<p><center>"
		print triggerlink("triggers.cgi", query)
		print "</center></p>"
		print "<p>"
		print "<center>"
		print "<table>"
		print plot_table_row("rateplot.cgi", query)
		print plot_table_row("conf_vs_time.cgi", query)
		print plot_table_row("glitchmetric.cgi", query)
		print plot_table_row("tfplot.cgi", query)
		if "Injections" in query.instrument:
			print plot_table_row("injections.cgi", query)
		print plot_table_row("rate_vs_freq.cgi", query)
		print plot_table_row("conf_vs_freq.cgi", query)
		print plot_table_row("rate_vs_conf.cgi", query)
		print plot_table_row("rate_vs_snr.cgi", query)
		print "</table>"
		print "</center>"
		print "</p>"
	else:
		print "<p><center>"
		print "<b>This code undergoing maintenance.  Please check back later.</b>"
		print "</center></p>"
	print "<p>"
	print formnotes()
	print "</p>"

print """</body>"""
print """</html>"""
