#
# Copyright (C) 2006,2013  Kipp Cannon
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
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import sys

from glue.ligolw import ligolw
from glue.ligolw import lsctables
from glue.ligolw import utils as ligolw_utils
from lalburst import git_version


lsctables.use_in(ligolw.LIGOLWContentHandler)


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                 Command Line
#
# =============================================================================
#


def parse_command_line():
	parser = OptionParser(
		version = "Name: %%prog\n%s" % git_version.verbose_msg,
		description = "Plot projections of the time slides contained in the input files onto a two-instrument cross section of the delay space.  One plot is generated for each input file."
	)
	parser.add_option("-o", "--output", metavar = "filename", help = "Set output file name (required).  If the string \"%n\" occurs in the filename, it will be replaced with the plot number starting from 0.  If the output filename does not contain a \"%n\" in it and more than one plot is generated then the plots will over-write one another.")
	parser.add_option("-x", "--x-instrument", metavar = "instrument", help = "Plot this instrument's offsets along the x axis (required).")
	parser.add_option("-y", "--y-instrument", metavar = "instrument", help = "Plot this instrument's offsets along the y axis (required).")
	parser.add_option("-n", "--require-instruments", metavar = "instrument[,instrument...]", help = "Plot only time slides involving exactly these instruments.")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()

	if not options.output:
		raise ValueError("no output file specified")

	if not options.x_instrument or not options.y_instrument:
		raise ValueError("must set instruments for x and y axes")

	if options.require_instruments is not None:
		options.require_instruments = set(map(str.strip, options.require_instruments.split(",")))

	return options, (filenames or [None])


#
# =============================================================================
#
#                           Time Slide Table Helper
#
# =============================================================================
#


def time_slide_append(self, row):
	if row.time_slide_id not in self.dict:
		self.dict[row.time_slide_id] = {}
	self.dict[row.time_slide_id][row.instrument] = row.offset


def time_slide__end_of_columns(self):
	self.dict = {}


lsctables.TimeSlideTable.append = time_slide_append
lsctables.TimeSlideTable._end_of_columns = time_slide__end_of_columns


#
# =============================================================================
#
#                                     Plot
#
# =============================================================================
#


class Plot(object):
	def __init__(self, xmldoc, x_instrument, y_instrument, require_instruments = None):
		self.fig = figure.Figure()
		self.canvas = FigureCanvas(self.fig)
		self.axes = self.fig.gca()

		self.axes.grid(True)

		self.axes.set_xlabel("%s Offset (s)" % x_instrument)
		self.axes.set_ylabel("%s Offset (s)" % y_instrument)

		tisitable = lsctables.TimeSlideTable.get_table(xmldoc)
		# grab an offset from one of the diciontaries of offsets in
		# the time slide table.
		max_offset = min_offset = tisitable.dict.itervalues().next().itervalues().next()
		x = []
		y = []
		for offsets in tisitable.dict.itervalues():
			if require_instruments is None or require_instruments == set(offsets.keys()):
				x.append(offsets[x_instrument])
				y.append(offsets[y_instrument])
		min_offset = min(x + y)
		max_offset = max(x + y)

		self.axes.plot(x, y, "k+")
		self.axes.set_xlim([min_offset, max_offset])
		self.axes.set_ylim([min_offset, max_offset])

		if require_instruments is not None:
			self.axes.set_title("%d %s Time Slides" % (len(x), "+".join(sorted(require_instruments))))
		else:
			self.axes.set_title("%d Time Slides" % len(x))


#
# =============================================================================
#
#                                     Main
#
# =============================================================================
#


options, filenames = parse_command_line()


for n, filename in enumerate(filenames):
	xmldoc = ligolw_utils.load_filename(filename, verbose = options.verbose, contenthandler = ligolw.LIGOLWContentHandler)

	if options.verbose:
		print >>sys.stderr, "plotting ..."
	plot = Plot(xmldoc, options.x_instrument, options.y_instrument, require_instruments = options.require_instruments)

	output = options.output.replace("%n", "%d" % n)
	if options.verbose:
		print >>sys.stderr, "writing %s ..." % output
	plot.fig.savefig(output)
