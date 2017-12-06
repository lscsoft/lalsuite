#
# Copyright (C) 2006,2012  Kipp Cannon
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


from optparse import OptionParser
import matplotlib
matplotlib.rcParams.update({
	"text.usetex": False	# SVG backend doesn't support LaTeX
})
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import sys
import time


from glue import segments
from glue.ligolw import ligolw
from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw.utils import ligolw_add
from pylal import rate
import lal
from lal.utils import CacheEntry
from lalburst import date
from lalburst import git_version


bins_per_filterwidth = 21


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                 Command line
#
# =============================================================================
#


def parse_command_line():
	parser = OptionParser(
		version = "Name: %%prog\n%s" % git_version.verbose_msg,
		usage = """usage: %prog [options] cachename ...

Generate long time scale trigger rate plot, getting trigger file names from LAL
cache files."""
	)
	parser.add_option("-s", "--gps-start-time", metavar = "seconds", help = "Set start time of plot in GPS seconds (required).")
	parser.add_option("-e", "--gps-end-time", metavar = "seconds", help = "Set end time of plot in GPS seconds (required).")
	parser.add_option("-w", "--window", metavar = "seconds", type = "float", default = 3600.0, help = "Set width of averaging window in seconds (default = 3600.0).")
	parser.add_option("-i", "--instrument", metavar = "name", help = "Set instrument name (required).")
	parser.add_option("-o", "--output-base", metavar = "base", help = "Set base (no extension) of output file name (required).")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, cache_names = parser.parse_args()

	# check for required options
	required_options = ["gps_start_time", "gps_end_time", "output_base", "instrument"]
	missing_options = [option for option in required_options if getattr(options, option) is None]
	if missing_options:
		raise ValueError("missing required option(s) %s" % ", ".join("--%s" % option.replace("_", "-") for option in missing_options))

	# parse trigger cache files
	if not cache_names:
		raise ValueError("no cache files named on command line")
	cache = [CacheEntry(l) for name in cache_names for l in file(name)]

	# set segment
	options.gps_start_time = int(lal.LIGOTimeGPS(options.gps_start_time))
	options.gps_end_time = int(lal.LIGOTimeGPS(options.gps_end_time))
	options.segment = segments.segment(options.gps_start_time, options.gps_end_time)
	options.read_segment = options.segment.protract(5.0 * options.window)

	# filter cache entries and sort
	return options, [c.path for c in sorted(c for c in cache if options.read_segment.intersects(c.segment))]


options, filenames = parse_command_line()


#
# =============================================================================
#
#   Custom SnglBurstTable append() method to put triggers directly into bins
#
# =============================================================================
#


nbins = int(float(abs(options.read_segment)) / options.window) * bins_per_filterwidth
binning = rate.NDBins((rate.LinearBins(options.read_segment[0], options.read_segment[1], nbins),))
trigger_rate = rate.BinnedDensity(binning)

num_triggers = 0


def snglburst_append(self, row, verbose = options.verbose):
	global num_triggers, rate
	t = row.peak
	if t in options.read_segment:
		trigger_rate.count[t,] += 1.0
	num_triggers += 1
	if verbose and not (num_triggers % 125):
		print >>sys.stderr, "sngl_burst rows read:  %d\r" % num_triggers,


lsctables.SnglBurstTable.append = snglburst_append


#
# =============================================================================
#
#                                    Input
#
# =============================================================================
#


#
# a content handler that reads only sngl_burst and search_summary tables
#


def element_filter(name, attrs):
	return name == ligolw.Table.tagName and table.Table.TableName(attrs["Name"]) in (lsctables.SnglBurstTable.tableName, lsctables.SearchSummaryTable.tableName)


@lsctables.use_in
class ContentHandler(ligolw.PartialLIGOLWContentHandler):
	def __init__(self, doc):
		ligolw.PartialLIGOLWContentHandler.__init__(self, doc, element_filter)


#
# use ligolw_add module to load documents, and extract search_summary
# table's "in" segment list.
#


seglist = lsctables.SearchSummaryTable.get_table(ligolw_add.ligolw_add(ligolw.Document(), filenames, verbose = options.verbose, contenthandler = ContentHandler)).get_inlist().coalesce()


#
# =============================================================================
#
#                        How to generate X axis labels
#
# =============================================================================
#


def make_xticks(segment):
	# generate tick locations and labels
	values = list(date.UTCMidnights(*(lal.LIGOTimeGPS(t) for t in segment)))
	labels = []
	for tm in map(lambda t: time.struct_time(lal.GPSToUTC(t)), values):
		if tm.tm_wday == 1:	# tuesday
			labels.append(time.strftime("%H h, %a %b %d, %Y", tm))
		else:
			labels.append("")
	return map(float, values), labels


#
# =============================================================================
#
#                                    Figure
#
# =============================================================================
#


#
# build a figure whose axes are 3" wide per week, and whose height is the
# width of a US letter page (minus some typical printer margins).  if there
# isn't enough data, make sure plot is at least 10" wide.
#

def newfig(segment):
	weeks = float(abs(segment)) / 86400.0 / 7.0	# FIXME: leep seconds?
	border = (0.5, 0.75, 0.125, 0.625)	# inches
	width = max(10.0, weeks * 3.0 + border[0] + border[2])	# inches
	height = 8.0	# inches
	fig = figure.Figure()
	canvas = FigureCanvas(fig)
	fig.set_size_inches(width, height)
	fig.gca().set_position([border[0] / width, border[1] / height, (width - border[0] - border[2]) / width, (height - border[1] - border[3]) / height])
	return fig


fig = newfig(options.segment)
axes = fig.gca()


rate.filter_array(trigger_rate.array, rate.gaussian_window(bins_per_filterwidth))


axes.plot(trigger_rate.centres()[0], trigger_rate.at_centres())


axes.set_xlim(list(options.segment))
axes.grid(True)


for seg in ~seglist & segments.segmentlist([options.segment]):
	axes.axvspan(seg[0], seg[1], facecolor = "k", alpha = 0.2)


axes.set_title("%s Excess Power Trigger Rate vs. Time\n(%d Triggers, %g s Moving Average)" % (options.instrument, num_triggers, options.window))


ticks = make_xticks(options.segment)
axes.set_xticks(ticks[0])
axes.set_xticklabels(ticks[1], horizontalalignment = "right", fontsize = 10, rotation = 10)
axes.set_xlabel("UTC")
#axes.yticks(fontsize = 10)
axes.set_ylabel("Trigger Rate (Hz)")


for extension in (".pdf", ".png", ".svg"):
	filename = options.output_base + extension
	if options.verbose:
		print >>sys.stderr, "writing %s ..." % filename
	fig.savefig(filename)
