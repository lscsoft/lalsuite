#
# Copyright (C) 2009  Kipp Cannon
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


import math
import matplotlib
matplotlib.rcParams.update({
	"font.size": 8.0,
	"axes.titlesize": 10.0,
	"axes.labelsize": 10.0,
	"xtick.labelsize": 8.0,
	"ytick.labelsize": 8.0,
	"legend.fontsize": 8.0,
	"figure.dpi": 300,
	"savefig.dpi": 300,
	"text.usetex": True,
	"grid.linestyle": "-",
	"grid.linewidth": 0.25
})
from matplotlib import figure
from matplotlib import patches
from matplotlib import cm
from matplotlib import colorbar
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import numpy
from optparse import OptionParser
import sys


from lalburst import git_version
from pylal import snglcoinc
from lalburst import stringutils


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


golden_ratio = (1 + math.sqrt(5)) / 2


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
	parser.add_option("-f", "--format", metavar = "extension", action = "append", default = [], help = "Set the image output format by setting the filename extension (default = \"png\").  Can be given multiple times to generate plots in multiple formats.")
	parser.add_option("--no-filter", action = "store_true", help = "Do not apply smoothing/normalization filters to data, plot raw bin values.")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()

	if not options.format:
		options.format = ["png"]

	return options, (filenames or [None])


#
# =============================================================================
#
#                                     Blah
#
# =============================================================================
#


def clip_binned_array_1d(binnedarray, xlim):
	imin, = binnedarray.bins[xlim[0],]
	imax, = binnedarray.bins[xlim[1],]
	coords, = binnedarray.bins.centres()
	return coords[imin:imax], binnedarray.array[imin:imax]


def clip_binned_array_2d(binnedarray, xlim, ylim):
	imin, jmin = binnedarray.bins[xlim[0], ylim[0]]
	imax, jmax = binnedarray.bins[xlim[1], ylim[1]]
	xcoords, ycoords = binnedarray.bins.centres()
	return xcoords[imin:imax], ycoords[jmin:jmax], binnedarray.array[imin:imax,jmin:jmax]


def snr2_chi2_plot(key, background_xcoords, background_ycoords, background_data, injections_xcoords, injections_ycoords, injections_data, ncontours = 49):
	fig = figure.Figure(figsize=(3, 3))
	FigureCanvas(fig)
	axes = fig.add_axes((.15, .15, .95 - .15, .90 - .15))
	axes.loglog()

	background_data = numpy.nan_to_num(numpy.transpose(background_data))
	injections_data = numpy.nan_to_num(numpy.transpose(injections_data))

	if numpy.any(background_data) or numpy.any(injections_data):
		hi = math.log10(max(background_data.max(), injections_data.max()))
		contours = numpy.arange(hi - 10, hi, 0.5, dtype = "double")
		injections_cset = axes.contour(injections_xcoords, injections_ycoords, numpy.log10(injections_data), contours, cmap = cm.Reds)
		background_cset = axes.contour(background_xcoords, background_ycoords, numpy.log10(background_data), contours, cmap = cm.Greys)
	axes.set_xlim([min(background_xcoords[0], injections_xcoords[0]), max(background_xcoords[-1], injections_xcoords[-1])])
	axes.set_ylim([min(background_ycoords[0], injections_ycoords[0]), max(background_ycoords[-1], injections_ycoords[-1])])
	#cbar = fig.add_axes((.75,.15,.1,.75))
	#colorbar.Colorbar(cbar, cset)
	instrument = key.split("-")[0]
	axes.set_title(r"%s Event Density in $\rho^{2}$--$\chi^{2}$ Plane" % instrument)
	axes.set_ylabel(r"$\chi^{2} / \mathrm{DOF}$")
	axes.set_xlabel(r"$\rho^{2}$")
	axes.xaxis.grid(True, which = "major,minor")
	axes.yaxis.grid(True, which = "major,minor")
	return fig


def dA_plot(key, background_coords, background_data, injections_coords, injections_data):
	fig = figure.Figure(figsize=(4, 4 / golden_ratio))
	FigureCanvas(fig)
	axes = fig.add_axes((.14, .15, .98 - .14, .90 - .15))

	background_data = numpy.nan_to_num(background_data)
	injections_data = numpy.nan_to_num(injections_data)

	if numpy.any(background_data):
		axes.semilogy(background_coords, background_data, "k-", label = "Background")
	if numpy.any(injections_data):
		axes.semilogy(injections_coords, injections_data, "r-", label = "Injections")
	ymin, ymax = axes.get_ylim()
	axes.set_ylim((max(ymin, ymax * 1e-6), ymax))
	axes.legend(loc = "upper left")

	instrument1, instrument2 = key.split("-")[:2]
	axes.set_title(r"%s--%s Amplitude Ratio Distribution" % (instrument1, instrument2))
	axes.set_ylabel("Probability Density")
	axes.set_xlabel(r"$\log_{10} \left|A_{\mathrm{%s}} / A_{\mathrm{%s}}\right|$" % (instrument1, instrument2))
	axes.xaxis.grid(True, which = "major,minor")
	axes.yaxis.grid(True, which = "major,minor")
	return fig

def df_plot(key, background_coords, background_data, injections_coords, injections_data):
	fig = figure.Figure(figsize=(4, 4 / golden_ratio))
	FigureCanvas(fig)
	axes = fig.add_axes((.14, .15, .98 - .14, .90 - .15))

	background_data = numpy.nan_to_num(background_data)
	injections_data = numpy.nan_to_num(injections_data)

	if numpy.any(background_data):
		axes.semilogy(background_coords, background_data, "k-", label = "Background")
	if numpy.any(injections_data):
		axes.semilogy(injections_coords, injections_data, "r-", label = "Injections")
	ymin, ymax = axes.get_ylim()
	axes.set_ylim((max(ymin, ymax * 1e-6), ymax))
	axes.legend(loc = "upper left")

	instrument1, instrument2 = key.split("-")[:2]
	axes.set_title(r"%s--%s Frequency Cutoff Asymmetry Distribution" % (instrument1, instrument2))
	axes.set_ylabel("Probability Density")
	axes.set_xlabel(r"$\left({f_{\mathrm{cut}}}_{\mathrm{%s}} - {f_{\mathrm{cut}}}_{\mathrm{%s}}\right) / \frac{1}{2} \left({f_{\mathrm{cut}}}_{\mathrm{%s}} + {f_{\mathrm{cut}}}_{\mathrm{%s}}\right)$" % (instrument1, instrument2,instrument1, instrument2))
	axes.xaxis.grid(True, which = "major,minor")
	axes.yaxis.grid(True, which = "major,minor")
	return fig


def dt_plot(key, background_coords, background_data, injections_coords, injections_data):
	fig = figure.Figure(figsize=(4, 4 / golden_ratio))
	FigureCanvas(fig)
	axes = fig.add_axes((.14, .15, .98 - .14, .90 - .15))

	background_data = numpy.nan_to_num(background_data)
	injections_data = numpy.nan_to_num(injections_data)

	if numpy.any(background_data):
		axes.semilogy(background_coords, background_data, "k-", label = "Background")
	if numpy.any(injections_data):
		axes.semilogy(injections_coords, injections_data, "r-", label = "Injections")
	ymin, ymax = axes.get_ylim()
	axes.set_ylim((max(ymin, ymax * 1e-6), ymax))
	axes.legend(loc = "upper left")

	instrument1, instrument2 = key.split("-")[:2]
	axes.set_title(r"%s--%s Arrival Time Difference Distribution" % (instrument1, instrument2))
	axes.set_ylabel("Probability Density")
	axes.set_xlabel(r"$t_{\mathrm{%s}} - t_{\mathrm{%s}}$ (s)" % (instrument1, instrument2))
	axes.xaxis.grid(True, which = "major,minor")
	axes.yaxis.grid(True, which = "major,minor")
	return fig


def nevents_plot(key, xcoords, background_data, injections_data):
	fig = figure.Figure(figsize = (4, 4 / golden_ratio))
	FigureCanvas(fig)
	axes = fig.add_axes((.14, .15, .98 - .14, .90 - .15))

	background_data = numpy.nan_to_num(background_data)
	injections_data = numpy.nan_to_num(injections_data)

	axes.plot(xcoords, background_data, "ko-", label = "Background")
	axes.plot(xcoords, injections_data, "ro-", label = "Injections")
	axes.legend(loc = "lower left")

	axes.set_title("Number of Events Found in Coincidence")
	axes.set_xlabel("Number of Events $N$")
	axes.set_ylabel("$P(N)$")
	axes.xaxis.grid(True, which = "major,minor")
	axes.yaxis.grid(True, which = "major,minor")

	return fig

def instrumentgroup_plot(key, xcoords, background_data, injections_data):
	fig = figure.Figure(figsize = (4, 4 / golden_ratio))
	FigureCanvas(fig)
	axes = fig.add_axes((.14, .15, .98 - .14, .90 - .15))

	background_data = numpy.nan_to_num(background_data)
	injections_data = numpy.nan_to_num(injections_data)

	axes.semilogy(xcoords, background_data, "ko-", label = "Background")
	axes.semilogy(xcoords, injections_data, "ro-", label = "Injections")
	axes.set_ticklabels([",".join(sorted(stringutils.category_to_instruments(i))) for i in xcoords])
	axes.legend(loc = "upper right")

	axes.set_title("Instrument Combinations")
	axes.set_xlabel("Instrument Combination")
	axes.set_ylabel(r"$P(\mathrm{Instruments})$")
	axes.xaxis.grid(True, which = "major,minor")
	axes.yaxis.grid(True, which = "major,minor")

	return fig

def instrumentgroup_timingresidual_plot(key, background_xcoords, background_ycoords, background_data, injections_xcoords, injections_ycoords, injections_data):
	fig = figure.Figure(figsize = (4, 4 / golden_ratio))
	FigureCanvas(fig)
	axes = fig.add_axes((.14, .15, .98 - .14, .90 - .15))

	background_data = numpy.nan_to_num(background_data)
	injections_data = numpy.nan_to_num(injections_data)

	for i in range(len(background_xcoords)):
		if numpy.any(background_data[i,:]) or numpy.any(injections_data[i,:]):
			axes.semilogy(background_ycoords, background_data[i,:], label = "Background (%s)" % ", ".join(sorted(stringutils.category_to_instruments(i + 1))))
			axes.semilogy(injections_ycoords, injections_data[i,:], label = "Injections (%s)" % ", ".join(sorted(stringutils.category_to_instruments(i + 1))))
	ymin, ymax = axes.get_ylim()
	# find the smallest global maximum of all curves
	ymin = min(max(background_data[i,:].max(), injections_data[i,:].max()) for i in range(len(background_xcoords)))
	axes.set_ylim((ymin * 1e-6, ymax))
	axes.legend(loc = "upper right")

	axes.set_title("RSS Timing Residual")
	axes.set_xlabel("RSS Timing Residual (s)")
	axes.set_ylabel("Probability Density")
	axes.xaxis.grid(True, which = "major,minor")
	axes.yaxis.grid(True, which = "major,minor")

	return fig


#
# =============================================================================
#
#                                     Main
#
# =============================================================================
#


options, filenames = parse_command_line()


coincparamsdistributions, seglists = stringutils.load_likelihood_data(filenames, verbose = options.verbose)

if options.verbose:
	print >>sys.stderr, "computing event densities ..."
if options.no_filter:
	source_dicts = (coincparamsdistributions.background_rates, coincparamsdistributions.injection_rates, coincparamsdistributions.zero_lag_rates)
else:
	coincparamsdistributions.finish(verbose = options.verbose)
	source_dicts = (coincparamsdistributions.background_pdf, coincparamsdistributions.injection_pdf, coincparamsdistributions.zero_lag_pdf)

for (background_name, background_binnedarray), (injections_name, injections_binnedarray), (zero_lag_name, zero_lag_binnedarray) in zip(sorted(source_dicts[0].items()), sorted(source_dicts[1].items()), sorted(source_dicts[2].items())):
	assert injections_name == background_name
	assert zero_lag_name == background_name
	name = background_name
	instruments = set(name.split("_")) & set(stringutils.StringCoincParamsDistributions.instrument_categories)
	if name.endswith("_snr2_chi2"):
		if options.verbose:
			print >>sys.stderr, "generating plots for %s ..." % name
		background_xcoords, background_ycoords, background_data = clip_binned_array_2d(background_binnedarray, [10, 1e6], [.01, 1e4])
		injections_xcoords, injections_ycoords, injections_data = clip_binned_array_2d(injections_binnedarray, [10, 1e6], [.01, 1e4])
		zero_lag_xcoords, zero_lag_ycoords, zero_lag_data = clip_binned_array_2d(zero_lag_binnedarray, [10, 1e6], [.01, 1e4])
		fig = snr2_chi2_plot("%s" % name.replace("_", "-"), background_xcoords, background_ycoords, background_data, injections_xcoords, injections_ycoords, injections_data)
		for extension in options.format:
			outname = "%s.%s" % (name, extension)
			if options.verbose:
				print >>sys.stderr, "\twriting %s ..." % outname
			fig.savefig(outname)
	elif name.endswith("_dt"):
		if options.verbose:
			print >>sys.stderr, "generating plots for %s ..." % name
		dt = .010 + snglcoinc.light_travel_time(*instruments)
		background_coords, background_data = clip_binned_array_1d(background_binnedarray, (-dt, +dt))
		injections_coords, injections_data = clip_binned_array_1d(injections_binnedarray, (-dt, +dt))
		zero_lag_coords, zero_lag_data = clip_binned_array_1d(zero_lag_binnedarray, (-dt, +dt))
		fig = dt_plot("%s" % name.replace("_", "-"), background_coords, background_data, injections_coords, injections_data)
		for extension in options.format:
			outname = "%s.%s" % (name, extension)
			if options.verbose:
				print >>sys.stderr, "\twriting %s ..." % outname
			fig.savefig(outname)
	elif name.endswith("_dA"):
		if options.verbose:
			print >>sys.stderr, "generating plots for %s ..." % name
		background_coords, background_data = clip_binned_array_1d(background_binnedarray, (-2, +2))
		injections_coords, injections_data = clip_binned_array_1d(injections_binnedarray, (-2, +2))
		zero_lag_coords, zero_lag_data = clip_binned_array_1d(zero_lag_binnedarray, (-2, +2))
		fig = dA_plot("%s" % name.replace("_", "-"), background_coords, background_data, injections_coords, injections_data)
		for extension in options.format:
			outname = "%s.%s" % (name, extension)
			if options.verbose:
				print >>sys.stderr, "\twriting %s ..." % outname
			fig.savefig(outname)
	elif name.endswith("_df"):
		if options.verbose:
			print >>sys.stderr, "generating plots for %s ..." % name
		background_coords, background_data = clip_binned_array_1d(background_binnedarray, (-0.6, +0.6))
		injections_coords, injections_data = clip_binned_array_1d(injections_binnedarray, (-0.6, +0.6))
		zero_lag_coords, zero_lag_data = clip_binned_array_1d(zero_lag_binnedarray, (-0.6, +0.6))
		fig = df_plot("%s" % name.replace("_", "-"), background_coords, background_data, injections_coords, injections_data)
		for extension in options.format:
			outname = "%s.%s" % (name, extension)
			if options.verbose:
				print >>sys.stderr, "\twriting %s ..." % outname
			fig.savefig(outname)
	elif name == "nevents":
		if options.verbose:
			print >>sys.stderr, "generating plots for %s ..." % name
		fig = nevents_plot("%s" % name.replace("_", "-"), background_binnedarray.centres()[0], background_binnedarray.array, injections_binnedarray.array)
		for extension in options.format:
			outname = "%s.%s" % (name, extension)
			if options.verbose:
				print >>sys.stderr, "\twriting %s ..." % outname
			fig.savefig(outname)
	elif name == "instrumentgroup":
		if options.verbose:
			print >>sys.stderr, "generating plots for %s ..." % name
		fig = instrumentgroup_plot("%s" % name.replace("_", "-"), background_binnedarray.centres()[0], background_binnedarray.array, injections_binnedarray.array)
		for extension in options.format:
			outname = "%s.%s" % (name, extension)
			if options.verbose:
				print >>sys.stderr, "\twriting %s ..." % outname
			fig.savefig(outname)
	elif name == "instrumentgroup,rss_timing_residual":
		if options.verbose:
			print >>sys.stderr, "generating plots for %s ..." % name
		background_xcoords, background_ycoords, background_data = clip_binned_array_2d(background_binnedarray, [1, background_binnedarray.array.shape[0] + 0.5], [0, .02])
		injections_xcoords, injections_ycoords, injections_data = clip_binned_array_2d(injections_binnedarray, [1, background_binnedarray.array.shape[0] + 0.5], [0, .02])
		fig = instrumentgroup_timingresidual_plot("%s" % name.replace("_", "-").replace(",", "--"), background_xcoords, background_ycoords, background_data, injections_xcoords, injections_ycoords, injections_data)
		for extension in options.format:
			outname = "%s.%s" % (name.replace(",", "_"), extension)
			if options.verbose:
				print >>sys.stderr, "\twriting %s ..." % outname
			fig.savefig(outname)
