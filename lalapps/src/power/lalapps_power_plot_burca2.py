#
# Copyright (C) 2006  Kipp Cannon
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


try:
	from fpconst import PosInf, NegInf
except ImportError:
	# fpconst is not part of the standard library and might not be
	# available
	PosInf = float("+inf")
	NegInf = float("-inf")
import math
import numpy
from optparse import OptionParser
import sys


from lalburst import git_version
from lalburst import burca_tailor
from lalburst import SnglBurstUtils


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
		version = "Name: %%prog\n%s" % git_version.verbose_msg
	)
	parser.add_option("-b", "--base", metavar = "base", default = "plotburca2", help = "Set the prefix for output filenames (default = \"plotburca2\").")
	parser.add_option("-f", "--format", metavar = "format", default = "png", help = "Set the output image format (default = \"png\").")
	parser.add_option("-l", "--live-time-program", metavar = "program", default = "lalapps_power", help = "Set the name, as it appears in the process table, of the program whose search summary entries define the search live time (default = \"lalapps_power\").")
	parser.add_option("--with-zero-lag", action = "store_true", help = "Also plot the parameter distributions for zero lag events.")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()

	return options, (filenames or [None])


#
# =============================================================================
#
#                   Coincidence Parameter Distribution Plot
#
# =============================================================================
#


def add_to_plot(axes, red, black, blue, gmst, plottype = "LR"):
	if plottype == "P":
		if blue is not None:
			y = blue[:,gmst]
			l, = axes.plot(numpy.compress(y > 0, blue.centres()[0]), numpy.compress(y > 0, y), "b")
		y = black[:,gmst]
		l, = axes.plot(numpy.compress(y > 0, black.centres()[0]), numpy.compress(y > 0, y), "k")
		#l.set_dashes([8, 4])
		y = red[:,gmst]
		l, = axes.plot(numpy.compress(y > 0, red.centres()[0]), numpy.compress(y > 0, y), "r")
	elif plottype == "LR":
		y1 = red[:,gmst]
		y2 = black[:,gmst]
		l, = axes.plot(numpy.compress((y1 > 0) & (y2 > 0), black.centres()[0]), numpy.compress((y1 > 0) & (y2 > 0), y1 / y2), "k")
	else:
		raise ValueError(plottype)

	ymin, ymax = axes.get_ylim()
	if ymin < 1e-5 * ymax:
		ymin = 1e-5 * ymax
		axes.set_ylim((ymin, ymax))
	axes.set_xlim((-2.0, +2.0))


def add_to_dtplot(axes, red, black, blue, gmst, plottype = "LR"):
	x = numpy.arange(len(black[:,gmst]), dtype = "double")
	if plottype == "P":
		if blue is not None:
			y = blue[:,gmst]
			l, = axes.plot(numpy.compress(y > 0, x), numpy.compress(y > 0, y), "b")
		y = black[:,gmst]
		l, = axes.plot(numpy.compress(y > 0, x), numpy.compress(y > 0, y), "k")
		#l.set_dashes([8, 4])
		y = red[:,gmst]
		l, = axes.plot(numpy.compress(y > 0, x), numpy.compress(y > 0, y), "r")
	elif plottype == "LR":
		y1 = red[:,gmst]
		y2 = black[:,gmst]
		l, = axes.plot(numpy.compress((y1 > 0) & (y2 > 0), x), numpy.compress((y1 > 0) & (y2 > 0), y1 / y2), "k")
	else:
		raise ValueError(plottype)

	ymin, ymax = axes.get_ylim()
	if ymin < 1e-5 * ymax:
		ymin = 1e-5 * ymax
		axes.set_ylim((ymin, ymax))
	axes.set_xlim((0, len(x) - 1))

	if black.bins[0].scale > 150:
		xticks = (NegInf, -0.02, -0.01, -0.005, -0.002, 0, 0.002, 0.005, 0.01, 0.02, PosInf)
	else:
		xticks = (NegInf, -0.05, -0.02, -0.01, -0.005, 0, 0.005, 0.01, 0.02, 0.05, PosInf)
	xticklabels = [r"$-\infty$"] + [r"$%.3g$" % (1000 * x) for x in xticks[1:-1]] + [r"$\infty$"]
	axes.set_xticks([black.bins[x,:][0] for x in xticks])
	axes.set_xticklabels(xticklabels)


#
# distributions is a burca_tailor.BurcaCoincParamsDistributions instance
#


def plot_coinc_params(distributions, gmst, plottype = "LR", with_zero_lag = False):
	#
	# Create a figure.
	#

	fig = SnglBurstUtils.figure.Figure()
	SnglBurstUtils.FigureCanvas(fig)

	#
	# How many instrument pairs are there?
	#

	pairs = set(tuple(name.split("_")[:2]) for name in distributions.background_pdf.keys() if "_" in name)

	#
	# How many plots in each row?
	#

	n_horiz = len(pairs)

	#
	# How many rows?
	#

	n_vert = 5

	#
	# Each of the first len(pairs) sub plot's aspect ratio is the golden
	# ratio.
	#

	size = 4.0
	fig.set_size_inches(n_horiz * size, n_vert / ((1 + math.sqrt(5)) / 2) * size)

	#
	# Plot layout.
	#

	vlabel_allowance = (n_vert * .05) / fig.figheight.get()
	hlabel_allowance = (n_horiz * .1) / fig.figwidth.get()
	border = .07 / fig.figwidth.get()
	width = 1.0 / n_horiz - hlabel_allowance - 2 * border
	height = 1.0 / n_vert

	#
	# Iterate over instrument pairs.
	#

	for i, pair in enumerate(pairs):
		#
		# Construct the axes for this instrument pair.
		#

		left = float(i) / n_horiz + hlabel_allowance + border

		dt_axes = fig.add_axes((left, 0 * height + vlabel_allowance + border, width, height - vlabel_allowance - 2 * border))
		df_axes = fig.add_axes((left, 1 * height + vlabel_allowance + border, width, height - vlabel_allowance - 2 * border))
		dh_axes = fig.add_axes((left, 2 * height + vlabel_allowance + border, width, height - vlabel_allowance - 2 * border))
		dband_axes = fig.add_axes((left, 3 * height + vlabel_allowance + border, width, height - vlabel_allowance - 2 * border))
		ddur_axes = fig.add_axes((left, 4 * height + vlabel_allowance + border, width, height - vlabel_allowance - 2 * border))

		dt_axes.semilogy()
		df_axes.semilogy()
		dh_axes.semilogy()
		dband_axes.semilogy()
		ddur_axes.semilogy()

		dt_axes.set_xlabel(r"$t_{\mathrm{%s}} - t_{\mathrm{%s}}$ (ms)" % pair)
		df_axes.set_xlabel(r"$(f_{\mathrm{%s}} - f_{\mathrm{%s}}) / \left< f \right>$" % pair)
		dh_axes.set_xlabel(r"$({h_{\mathrm{rss}}}_{\mathrm{%s}} - {h_{\mathrm{rss}}}_{\mathrm{%s}}) / \left< h_{\mathrm{rss}} \right>$" % pair)
		dband_axes.set_xlabel(r"$(\Delta f_{\mathrm{%s}} - \Delta f_{\mathrm{%s}}) / \left< \Delta f \right>$" % pair)
		ddur_axes.set_xlabel(r"$(\Delta t_{\mathrm{%s}} - \Delta t_{\mathrm{%s}}) / \left< \Delta t \right>$" % pair)

		#
		# Plot the data on them.
		#

		prefix = "%s_%s_" % pair

		for axes, suffix in zip((df_axes, dh_axes, dband_axes, ddur_axes), ("df", "dh", "dband", "ddur")):
			red = distributions.injection_pdf[prefix + suffix]
			black = distributions.background_pdf[prefix + suffix]
			if with_zero_lag:
				blue = distributions.zero_lag_pdf[prefix + suffix]
			else:
				blue = None
			add_to_plot(axes, red, black, blue, gmst, plottype = plottype)

		suffix = "dt"
		red = distributions.injection_pdf[prefix + suffix]
		black = distributions.background_pdf[prefix + suffix]
		if with_zero_lag:
			blue = distributions.zero_lag_pdf[prefix + suffix]
		else:
			blue = None
		add_to_dtplot(dt_axes, red, black, blue, gmst, plottype = plottype)

	#
	# Done.
	#

	return fig


#
# =============================================================================
#
#                                     Main
#
# =============================================================================
#


options, filenames = parse_command_line()
filenames.sort()


distributions, seglists = burca_tailor.EPGalacticCoreCoincParamsDistributions.from_filenames(filenames, u"lalapps_burca_tailor", verbose = options.verbose)
if options.verbose:
	print >>sys.stderr, "applying filters ..."
distributions.finish()


for gmst in numpy.arange(0.0, 2 * math.pi, 2 * math.pi / 10):
	filename = "%s-%%s-%d-%d-%s.%s" % (options.base, int(seglists.extent_all()[0]), int(abs(seglists.extent_all())), ("%.2f" % gmst).replace(".", "_"), options.format)
	if options.verbose:
		print >>sys.stderr, "writing %s ..." % (filename % "P")
	plot_coinc_params(distributions, gmst, plottype = "P", with_zero_lag = options.with_zero_lag).savefig(filename % "P")
	if options.verbose:
		print >>sys.stderr, "writing %s ..." % (filename % "LR")
	plot_coinc_params(distributions, gmst, plottype = "LR").savefig(filename % "LR")


if options.verbose:
	print >>sys.stderr, "done."
