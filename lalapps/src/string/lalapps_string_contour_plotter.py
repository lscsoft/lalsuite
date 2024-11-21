##python
#
# Copyright (C) 2010  Andrew Mergl
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

#
# This is a Python replacement for the contour_plotter.m MATLAB code
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
	"text.usetex": True
})
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

import numpy
from optparse import OptionParser

from lalburst import git_version

__author__ = "Andrew Mergl <>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


def parse_command_line():
    parser = OptionParser(
        version = "Name: %%prog\n%s" % git_version.verbose_msg
    )
    parser.add_option("-v", "--verbose", action = "store_true",
                      help = "Be verbose.")
    parser.add_option("-t", "--live-time", dest="livetime",
                      type = "float",
                      help = "The total amount of live time in the run")

    options, filenames = parser.parse_args()
    if options.livetime is None:
        raise ValueError("No live time specified. Use -t or --live-time.")
    if filenames is None:
        raise ValueError("No data file specified.")
    return options, filenames


def cosmo_contour(x, y, avg, min, max, p, neventsUL, filename):
    fig = Figure()
    FigureCanvas(fig)
    fig.set_size_inches(4.5, 4.5 / ((1 + math.sqrt(5)) / 2))
    axes = fig.add_axes((.12, .15, .98 - .12, .90 - .15))
    axes.loglog()
    #axes.contourf(x, y, numpy.log(avg/p), 100, cmap = cm.bone)
    #axes.contourf(x, y, numpy.log(min/p), 100, cmap = cm.bone)
    #axes.contourf(x, y, numpy.log(max/p), 100, cmap = cm.bone)
    axes.contour(x, y, avg/p, [neventsUL], linestyles='solid', colors='r')
    axes.contour(x, y, min/p, [neventsUL], linestyles='dashed', colors='r')
    axes.contour(x, y, max/p, [neventsUL], linestyles='dashed', colors='r')
    axes.set_xlim([x.min(),x.max()])
    axes.set_ylim([y.min(),y.max()])
    axes.set_title(r"$P = %g$" % p)
    axes.set_xlabel(r"$G\mu$")
    axes.set_ylabel(r"$\epsilon$")
    fig.savefig(filename)

#main

options, filenames = parse_command_line()

#read cs_gamma file output
#the columns should be p,index,epsilon,Gmu,gammaAverage,gammaMin,gammaMax
p, n, eps, gmu, avg, min, max = [], [], [], [], [], [], []
for line in open(filenames[0], 'r'):
    line = line.strip()
    if line[0] in "%#":
        continue
    for arr, val in zip((p, n, eps, gmu, avg, min, max), line.split()):
        arr.append(float(val))

epsindex = dict((b, a) for a, b in enumerate(sorted(set(eps))))
gmuindex = dict((b, a) for a, b in enumerate(sorted(set(gmu))))

# check for rounding errors producing a wrong count of x and y values
assert len(epsindex) * len(gmuindex) == len(gmu)
# program only works if all values in the data file are for the same p and that p is 1
assert set(p) == set([1])

avgarr = numpy.zeros([len(epsindex), len(gmuindex)], dtype = "double")
minarr = numpy.zeros([len(epsindex), len(gmuindex)], dtype = "double")
maxarr = numpy.zeros([len(epsindex), len(gmuindex)], dtype = "double")

for e, g, val in zip(eps, gmu, avg):
	avgarr[epsindex[e], gmuindex[g]] = val
for e, g, val in zip(eps, gmu, min):
	minarr[epsindex[e], gmuindex[g]] = val
for e, g, val in zip(eps, gmu, max):
	maxarr[epsindex[e], gmuindex[g]] = val

avgarr *= options.livetime
minarr *= options.livetime
maxarr *= options.livetime

numberofeventsUL = -numpy.log(1-0.90)

x = numpy.array(sorted(gmuindex.keys()))
y = numpy.array(sorted(epsindex.keys()))

cosmo_contour(x, y, avgarr, minarr, maxarr, 1.000, numberofeventsUL, "gmu-eps-p1e+0")
cosmo_contour(x, y, avgarr, minarr, maxarr, 0.100, numberofeventsUL, "gmu-eps-p1e-1")
cosmo_contour(x, y, avgarr, minarr, maxarr, 0.010, numberofeventsUL, "gmu-eps-p1e-2")
cosmo_contour(x, y, avgarr, minarr, maxarr, 0.001, numberofeventsUL, "gmu-eps-p1e-3")
