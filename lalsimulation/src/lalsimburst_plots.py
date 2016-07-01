#!/usr/bin/env python
#
# Copyright (C) 2015  Kipp Cannon
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


"""
Script to regenerate plots for the LALSimBurst.c documentation
"""


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
	"path.simplify": True,
	"figure.dpi": 100,
	"savefig.dpi": 100,
	"text.usetex": True
})
from matplotlib import figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import numpy

import lalsimulation


#
# =============================================================================
#
#                                Sine-Gaussians
#
# =============================================================================
#


def sinegaussian_example(axes1, axes2, e, pol):
	axes1.grid(True)
	axes2.grid(True)

	hp, hc = lalsimulation.SimBurstSineGaussian(8.89, 250., 1., e, pol, 1./16384)
	print "measured hrss = %.17g" % lalsimulation.MeasureHrss(hp, hc)
	t = float(hp.epoch) + numpy.arange(0., len(hp.data.data)) * hp.deltaT
	axes1.plot(t, hp.data.data, label = r"$h_{+}$")
	t = float(hc.epoch) + numpy.arange(0., len(hc.data.data)) * hc.deltaT
	axes1.plot(t, hc.data.data, label = r"$h_{\times}$")
	axes1.set_ylim((-15, +15))
	axes1.set_xlabel(r"Time (seconds)")
	axes1.set_ylabel(r"Strain")
	axes1.set_title(r"$f_{0} = 250\,\mathrm{Hz}$, $Q = 8.89$, $h_{\mathrm{rss}} = 1$, $\epsilon = %g$, $\phi = %g \pi$" % (e, pol / math.pi))
	axes1.legend()
	axes2.plot(hc.data.data, hp.data.data, color = "k")
	axes2.set_xlim((-15, +15))
	axes2.set_ylim((-15, +15))
	axes2.set_xlabel(r"$h_{\times}$")
	axes2.set_ylabel(r"$h_{+}$")
	axes2.set_title(r"$f_{0} = 250\,\mathrm{Hz}$, $Q = 8.89$, $h_{\mathrm{rss}} = 1$, $\epsilon = %g$, $\phi = %g \pi$" % (e, pol / math.pi))
	axes2.set_aspect(1., adjustable = "box")

if 1:
	fig = figure.Figure()
	FigureCanvas(fig)
	width = 200.
	fig.set_size_inches(1.2 * width / 25.4, 7. * width / 25.4 / ((1 + math.sqrt(5)) / 2))
	sinegaussian_example(fig.add_subplot(9, 2,  1), fig.add_subplot(9, 2,  2), 0., 0.)
	sinegaussian_example(fig.add_subplot(9, 2,  3), fig.add_subplot(9, 2,  4), .8, 0.)
	sinegaussian_example(fig.add_subplot(9, 2,  5), fig.add_subplot(9, 2,  6), 1., 0.)
	sinegaussian_example(fig.add_subplot(9, 2,  7), fig.add_subplot(9, 2,  8), 0., math.pi/4.)
	sinegaussian_example(fig.add_subplot(9, 2,  9), fig.add_subplot(9, 2, 10), .8, math.pi/4.)
	sinegaussian_example(fig.add_subplot(9, 2, 11), fig.add_subplot(9, 2, 12), 1., math.pi/4.)
	sinegaussian_example(fig.add_subplot(9, 2, 13), fig.add_subplot(9, 2, 14), 0., math.pi/2.)
	sinegaussian_example(fig.add_subplot(9, 2, 15), fig.add_subplot(9, 2, 16), .8, math.pi/2.)
	sinegaussian_example(fig.add_subplot(9, 2, 17), fig.add_subplot(9, 2, 18), 1., math.pi/2.)
	fig.tight_layout()
	fig.savefig("lalsimburst_sinegaussianexamples.svg")


#
# =============================================================================
#
#                             String Cusp Examples
#
# =============================================================================
#


def stringcusp_example(axes, fhigh):
	axes.grid(True)
	hp, hc = lalsimulation.GenerateStringCusp(1., fhigh, 1./16384)
	t = float(hp.epoch) + numpy.arange(0., len(hp.data.data)) * hp.deltaT
	axes.plot(t, hp.data.data, label = r"$h_{+}$")
	t = float(hc.epoch) + numpy.arange(0., len(hc.data.data)) * hc.deltaT
	axes.plot(t, hc.data.data, label = r"$h_{\times}$")
	axes.set_xlabel(r"Time (seconds)")
	axes.set_ylabel(r"Strain")
	axes.set_title(r"$A = 1$, $f_{\mathrm{high}} = %g\,\mathrm{Hz}$" % fhigh)
	axes.legend()

if 1:
	fig = figure.Figure()
	FigureCanvas(fig)
	width = 200.
	fig.set_size_inches(width / 25.4, .75 * width / 25.4 / ((1 + math.sqrt(5)) / 2))
	stringcusp_example(fig.gca(), 128.)
	fig.tight_layout()
	fig.savefig("lalsimburst_stringcuspexamples.svg")


#
# =============================================================================
#
#                               BTLWNB Examples
#
# =============================================================================
#


def btlwnb_example(axes1, axes2, f0, dof, duration, e, pol):
	axes1.grid(True)
	axes2.grid(True)

	# dof/2 = dof in each polarization
	bandwidth = dof/2 * (2./math.pi) / duration

	import lal
	import lalburst
	simburst_row = lalburst.CreateSimBurst()
	simburst_row.waveform = "BTLWNB"
	simburst_row.egw_over_rsquared = 1. / (lal.GMSUN_SI * 4 / lal.C_SI / lal.PC_SI / lal.PC_SI)
	simburst_row.waveform_number = 0
	simburst_row.duration = duration
	simburst_row.frequency = f0
	simburst_row.bandwidth = bandwidth
	simburst_row.pol_ellipse_e = e
	simburst_row.pol_ellipse_angle = pol
	hp, hc = lalburst.GenerateSimBurst(simburst_row, 1./16384)

	t = float(hp.epoch) + numpy.arange(0., len(hp.data.data)) * hp.deltaT
	axes1.plot(t, hp.data.data, label = r"$h_{+}$")
	t = float(hc.epoch) + numpy.arange(0., len(hc.data.data)) * hc.deltaT
	axes1.plot(t, hc.data.data, label = r"$h_{\times}$")
	axes1.set_ylim((-0.015, +0.015))
	axes1.set_xlabel(r"Time (seconds)")
	axes1.set_ylabel(r"Strain")
	axes1.set_title(r"$\int(\dot{h}_{+}^{2} + \dot{h}_{\times}^{2})\,\mathrm{d}t = 1$, $f_{0} = %g\,\mathrm{Hz}$, $\Delta t = %g\,\mathrm{s}$, $%g\,\mathrm{DOF}$ ($\Delta f = %g\,\mathrm{Hz}$), $\epsilon = %g$, $\phi = %g \pi$" % (f0, duration, dof, bandwidth, e, pol / math.pi))
	axes1.legend()
	axes2.plot(hc.data.data, hp.data.data, color = "k")
	axes2.set_xlim((-0.015, +0.015))
	axes2.set_ylim((-0.015, +0.015))
	axes2.set_xlabel(r"$h_{\times}$")
	axes2.set_ylabel(r"$h_{+}$")
	axes2.set_aspect(1., adjustable = "box")

if 1:
	fig = figure.Figure()
	FigureCanvas(fig)
	width = 200.
	fig.set_size_inches(1.2 * width / 25.4, 7. * width / 25.4 / ((1 + math.sqrt(5)) / 2))
	btlwnb_example(fig.add_subplot(10, 2,  1), fig.add_subplot(10, 2,  2), 250., 16., 0.016, 0., 0.)
	btlwnb_example(fig.add_subplot(10, 2,  3), fig.add_subplot(10, 2,  4), 250.,  4., 0.008, 0., 0.)
	btlwnb_example(fig.add_subplot(10, 2,  5), fig.add_subplot(10, 2,  6), 250.,  4., 0.008, 0.85, 0.)
	btlwnb_example(fig.add_subplot(10, 2,  7), fig.add_subplot(10, 2,  8), 250.,  4., 0.008, 1., 0.)
	btlwnb_example(fig.add_subplot(10, 2,  9), fig.add_subplot(10, 2, 10), 250.,  2., 0.008, 0., 0.)
	btlwnb_example(fig.add_subplot(10, 2, 11), fig.add_subplot(10, 2, 12), 250.,  2., 0.008, 0., math.pi / 2.)
	btlwnb_example(fig.add_subplot(10, 2, 13), fig.add_subplot(10, 2, 14), 250.,  2., 0.008, 0.85, 0.)
	btlwnb_example(fig.add_subplot(10, 2, 15), fig.add_subplot(10, 2, 16), 250.,  2., 0.008, 1., 0.)
	btlwnb_example(fig.add_subplot(10, 2, 17), fig.add_subplot(10, 2, 18),   0.,  2., 0.008, 0., 0.)
	btlwnb_example(fig.add_subplot(10, 2, 19), fig.add_subplot(10, 2, 20),   0.,  2., 0.008, 1., 0.)
	fig.tight_layout()
	fig.savefig("lalsimburst_btlwnbexamples.svg")

