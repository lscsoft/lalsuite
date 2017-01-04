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


"""
Program for plotting LAL time and frequency series data from a LIGO Light
Weight XML file.
"""


import matplotlib
matplotlib.use("Agg")
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
matplotlib.rcParams.update({
	"text.usetex": True
})
import math
import numpy
from optparse import OptionParser
import re

import lal

from glue.ligolw import ligolw
from glue.ligolw import array
from glue.ligolw import utils
from lalburst import git_version


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
	parser.add_option("-a", "--axes", metavar = "[linlin|linlog|loglin|loglog]", default = "loglog", help = "Set the axes types (default = \"loglog\").")
	parser.add_option("-v", "--verbose", action = "store_true", help = "Be verbose.")
	options, filenames = parser.parse_args()

	if options.axes not in ("linlin", "linlog", "loglin", "loglog"):
		raise ValueError, "unrecognized --axes %s" % options.axes

	return options, (filenames or [None])


#
# =============================================================================
#
#                      LAL Metadata Extraction for Arrays
#
# =============================================================================
#


def LALmetadata(array_elem):
	"""
	Given an Array element in a LIGO Light Weight document tree,
	extract the LAL time or frequency series metadata for the contents
	of the array.
	"""
	lw_elem = array_elem.parentNode

	comment_elem = lw_elem.getElementsByTagName(ligolw.Comment.tagName)
	if comment_elem:
		comment = comment_elem[0].pcdata
	else:
		comment = None

	time_elem = lw_elem.getElementsByTagName(ligolw.Time.tagName)
	# FIXME:  LIGOTimeGPS needs to be able to handle unicode input
	epoch = lal.LIGOTimeGPS(str(time_elem[0].pcdata))

	param_elem = lw_elem.getElementsByTagName(ligolw.Param.tagName)
	f0 = float(param_elem[0].pcdata)

	dim_elem = array_elem.getElementsByTagName(ligolw.Dim.tagName)
	xaxis = dim_elem[0].Name
	xunit = dim_elem[0].Unit
	y_is_complex = dim_elem[1].pcdata == "3"

	name = array_elem.Name
	yunit = array_elem.Unit

	return {"comment": comment, "epoch": epoch, "f0": f0, "xaxis": xaxis, "xunit": xunit, "y_is_complex": y_is_complex, "name": name, "yunit": yunit}


#
# =============================================================================
#
#                                 TeX Helpers
#
# =============================================================================
#


def texifyunit(unit):
	"""
	Put braces around exponents, replace spaces with "\," (this space),
	then wrap in "$\mathrm{}$".
	"""
	# put braces around exponents
	unit = re.sub(r"(\S+)\^(\S+)", r"\1^{\2}", unit)

	# replace spaces with \, (thin space)
	unit = re.sub(r"\s+", r" \, ", unit)

	# wrap in $\mathrm{}$
	return "$\mathrm{" + unit + "}$"


def texifychannel(channel):
	"""
	Replace "_" with "\_".
	"""
	return re.sub(r"_", r"\_", channel)


#
# =============================================================================
#
#                                Plotting Loop
#
# =============================================================================
#


options, filenames = parse_command_line()


for filename in filenames:
	#
	# Load file
	#

	xmldoc = utils.load_filename(filename, verbose = options.verbose)
	if filename.endswith(".gz"):
		basename = filename[:-7]
	else:
		basename = filename[:-4]

	#
	# Iterate over arrays in file
	#

	for n, a in enumerate(xmldoc.getElementsByTagName(ligolw.Array.tagName)):
		if options.verbose:
			print "\tfound %s array \"%s\" ..." % ("x".join(map(str, a.array.shape)), a.Name)

		metadata = LALmetadata(a)

		fig = figure.Figure()
		FigureCanvasAgg(fig)
		# 6.5" wide, golden ratio high
		fig.set_size_inches(6.5, 6.5 / ((1 + math.sqrt(5)) / 2))
		axes = fig.gca()
		axes.grid(True)
		if options.axes == "linlin":
			pass
		elif options.axes == "linlog":
			axes.semilogy()
		elif options.axes == "loglin":
			axes.semilogx()
		else:
			# "loglog"
			axes.loglog()
		for i in range(1, a.array.shape[0]):
			if options.axes == "linlin":
				axes.plot(a.array[0], a.array[i])
			elif options.axes == "linlog":
				axes.plot(a.array[0], numpy.fabs(a.array[i]))
			elif options.axes == "loglin":
				axes.plot(numpy.fabs(a.array[0]), a.array[i])
			else:
				# "loglog"
				axes.plot(numpy.fabs(a.array[0]), numpy.fabs(a.array[i]))

		axes.set_ylabel("%s (%s)" % (texifychannel(metadata["name"]), texifyunit(metadata["yunit"])))
		axes.set_xlabel("%s (%s)" % (metadata["xaxis"], texifyunit(metadata["xunit"])))

		title = "%s at GPS %s s" % (texifychannel(metadata["name"]), str(metadata["epoch"]))
		if metadata["comment"]:
			title += " (%s)" % metadata["comment"]
		axes.set_title(title)

		plotname = "%s_%03d.png" % (basename, n)
		fig.savefig(plotname)
		if options.verbose:
			print "\t\tsaved as %s" % plotname
		del fig, axes

