#
# Copyright (C) 2015 Chris Pankow
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
"""Inspiral template bank plotter"""

import sys
import os
import argparse
from math import sqrt

import numpy

import matplotlib
matplotlib.use("agg")
from matplotlib import pyplot, cm

from glue.ligolw import utils, ligolw, lsctables
lsctables.use_in(ligolw.LIGOLWContentHandler)

from lalinference.rapid_pe import amrlib

__author__ = 'Chris Pankow <chris.pankow@ligo.org>'

def get_intr_params(tmplt_bank, template_id=None, params=("mass1", "mass2")):
    for t in tmplt_bank:
        if int(t.event_id) == template_id:
            return tuple(getattr(t, p) for p in params)
    if template_id is not None:
        raise ValueError("Did not find template ID %d in bank" % template_id)

    return numpy.array([tuple(getattr(t, p) for t in tmplt_bank) for p in params])

# Stolen shamelessly from:
# https://nipunbatra.github.io/2014/08/latexify/
def latexify(fig_width=None, fig_height=None, columns=1, pdict={}):
    """Set up matplotlib's RC params for LaTeX plotting.
    Call this before plotting a figure.

    Parameters
    ----------
    fig_width : float, optional, inches
    fig_height : float, optional, inches
    columns : {1, 2}
    """

    # code adapted from http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples

    # Width and max height in inches for IEEE journals taken from
    # computer.org/cms/Computer.org/Journal%20templates/transactions_art_guide.pdf

    assert(columns in [1,2])

    if fig_width is None:
        fig_width = 3.39 if columns==1 else 6.9 # width in inches

    if fig_height is None:
        golden_mean = (sqrt(5)-1.0)/2.0    # Aesthetic ratio
        fig_height = fig_width*golden_mean # height in inches

    MAX_HEIGHT_INCHES = 8.0
    if fig_height > MAX_HEIGHT_INCHES:
        print("WARNING: fig_height too large:" + fig_height + 
              "so will reduce to" + MAX_HEIGHT_INCHES + "inches.")
        fig_height = MAX_HEIGHT_INCHES

    params = {'backend': 'ps',
              'text.latex.preamble': ['\usepackage{gensymb}'],
              'axes.labelsize': 4, # fontsize for x and y labels (was 10)
              'axes.titlesize': 6,
              'text.fontsize': 6, # was 10
              'legend.fontsize': 6, # was 10
              'xtick.labelsize': 4,
              'ytick.labelsize': 4,
              'text.usetex': True,
              'figure.figsize': [fig_width,fig_height],
              'figure.dpi': 600,
              'savefig.dpi': 600,
              'font.family': 'serif'
    }
    # User specified params
    params.update(pdict)

    matplotlib.rcParams.update(params)

latexify()

argp = argparse.ArgumentParser()
argp.add_argument("-t", "--template-bank", help="Template bank XML file.")
#argp.add_argument("-I", "--template-id", action="append", help="Mark this template in the bank. Provide multiple times to indicate many.")
argp.add_argument("-I", "--template-id", type=int, help="Mark this template in the bank.")
argp.add_argument("-f", "--plot-name-bank", help="Name of output bank plot.")
argp.add_argument("-d", "--plot-name-bank-density", help="Name of output bank density plot.")
argp.add_argument("-v", "--verbose", action="store_true", help="Be verbose.")
argp.add_argument("-T", "--no-title", action="store_true", help="Omit title on plot.")
args = argp.parse_args()

if args.template_bank:
    xmldoc = utils.load_filename(args.template_bank, contenthandler=ligolw.LIGOLWContentHandler)
    tmplt_bank = lsctables.SnglInspiralTable.get_table(xmldoc)
else:
    exit("Template bank file must be supplied on the command line.")

if args.template_id:
    m1_in, m2_in = get_intr_params(tmplt_bank, args.template_id)
    s1z_in, s2z_in = get_intr_params(tmplt_bank, args.template_id, ("spin1z", "spin2z"))
m1, m2 = get_intr_params(tmplt_bank)
s1z, s2z = get_intr_params(tmplt_bank, params=("spin1z", "spin2z"))
chi = (s1z*m1 + s2z*m2) / (m1 + m2)

#
# Mass based plots
#
pyplot.subplot(2,2,1)
pyplot.grid(True, color='0.5', linewidth=0.125, linestyle='-')
pyplot.scatter(m1, m2, c=chi, s=1, edgecolor='none', cmap=cm.gist_earth)
if args.template_id:
    pyplot.scatter([m1_in], [m2_in], s=5, marker='x', color='m')
pyplot.xlabel("mass 1 ($M_{\odot}$)")
pyplot.ylabel("mass 2 ($M_{\odot}$)")
pyplot.xlim(0.95 * m1.min(), 1.05 * m1.max())
pyplot.ylim(0.95 * m2.min(), 1.05 * m2.max())
pyplot.locator_params(nbins=4)
pyplot.gca().set_axisbelow(True)

#
# Transform to mchirp / eta
#
mc, eta = amrlib.apply_transform(numpy.array((m1, m2)).T, ("mass1", "mass2"), "mchirp_eta").T

pyplot.subplot(2,2,2)
pyplot.grid(True, color='0.5', linewidth=0.125, linestyle='-')
pyplot.scatter(mc, eta, c=chi, s=1, edgecolor='none', cmap=cm.gist_earth)
#pyplot.scatter(mc, eta, s=1, edgecolor='none', cmap=cm.gist_earth)
if args.template_id:
    mc_in, eta_in = amrlib.apply_transform(numpy.array(([m1_in], [m2_in])).T, ("mass1", "mass2"), "mchirp_eta").T
    pyplot.scatter([mc_in], [eta_in], s=5, marker='x', color='m')
pyplot.xlabel("$\mathcal{M}_c$ ($M_{\odot}$)")
pyplot.ylabel("$\eta$")
pyplot.xlim(0.95 * mc.min(), 1.05 * mc.max())
pyplot.ylim(0.95 * eta.min(), 1.05 * eta.max())
pyplot.locator_params(nbins=4)
pyplot.gca().set_axisbelow(True)

#
# Transform to tau0 / tau3
#
# FIXME: default flow is 40
tau0, tau3 = amrlib.apply_transform(numpy.array((m1, m2)).T, ("mass1", "mass2"), "tau0_tau3").T
pyplot.subplot(2,2,3)
pyplot.grid(True, color='0.5', linewidth=0.125, linestyle='-')
pyplot.scatter(tau0, tau3, c=chi, s=1, edgecolor='none', cmap=cm.gist_earth)
#pyplot.scatter(tau0, tau3, s=1, edgecolor='none', cmap=cm.gist_earth)
if args.template_id:
    tau0_in, tau3_in = amrlib.apply_transform(numpy.array(([m1_in], [m2_in])).T, ("mass1", "mass2"), "tau0_tau3").T
    pyplot.scatter([tau0_in], [tau3_in], s=5, marker='x', color='m')
pyplot.xlabel(r"$\tau_0$")
pyplot.ylabel(r"$\tau_3$")
pyplot.xlim(0.95 * tau0.min(), 1.05 * tau0.max())
pyplot.ylim(0.95 * tau3.min(), 1.05 * tau3.max())
pyplot.locator_params(nbins=4)
pyplot.gca().set_axisbelow(True)

#
# Spin based plots
#
pyplot.subplot(2,2,4)
pyplot.grid(True, color='0.5', linewidth=0.125, linestyle='-')
res = pyplot.scatter(s1z, s2z, c=chi, s=1, edgecolor='none', cmap=cm.gist_earth, vmin=0.0, vmax=1.0)
if args.template_id:
    pyplot.scatter([s1z_in], [s2z_in], s=5, marker='x', color='m')
pyplot.xlabel(r"$s_{1,z}$")
pyplot.ylabel(r"$s_{2,z}$")
pyplot.xlim(-1, 1)
pyplot.ylim(-1, 1)
pyplot.locator_params(nbins=4)
pyplot.gca().set_axisbelow(True)

if args.no_title is not None and not args.no_title:
    pyplot.suptitle("Template bank with spin parameters")
#pyplot.tight_layout()

fig = pyplot.gcf()
fig.subplots_adjust(wspace=0.3, hspace=0.5)
ax = fig.add_axes([0.92, 0.15, 0.025, 0.7])
cb = fig.colorbar(res, cax=ax)
cb.ax.set_xlabel("$\chi$")

pyplot.savefig(args.plot_name_bank or "tmplt_bank.png")

pyplot.figure()

#
# Density - Mass based plots
#
pyplot.subplot(2,2,1)

h, bx, by = numpy.histogram2d(m1, m2, bins=(20, 20))
h = numpy.ma.masked_equal(h.T, 0)
pyplot.pcolor(bx, by, h, cmap=cm.gist_earth)

if args.template_id:
    pyplot.scatter([m1_in], [m2_in], s=5, marker='x', color='m')
pyplot.xlabel("mass 1 ($M_{\odot}$)")
pyplot.ylabel("mass 2 ($M_{\odot}$)")
pyplot.xlim(0.95 * m1.min(), 1.05 * m1.max())
pyplot.ylim(0.95 * m2.min(), 1.05 * m2.max())
pyplot.locator_params(nbins=4)
pyplot.gca().set_axisbelow(True)
pyplot.grid(True, color='0.5', linewidth=0.125, linestyle='-')

pyplot.subplot(2,2,2)

h, bx, by = numpy.histogram2d(mc, eta, bins=(20, 20))
h = numpy.ma.masked_equal(h.T, 0)
pyplot.pcolor(bx, by, h, cmap=cm.gist_earth)

if args.template_id:
    pyplot.scatter([mc_in], [eta_in], s=5, marker='x', color='m')
pyplot.xlabel("$\mathcal{M}_c$ ($M_{\odot}$)")
pyplot.ylabel("$\eta$")
pyplot.xlim(0.95 * mc.min(), 1.05 * mc.max())
pyplot.ylim(0.95 * eta.min(), 1.05 * eta.max())
pyplot.locator_params(nbins=4)
pyplot.grid(True, color='0.5', linewidth=0.125, linestyle='-')
pyplot.gca().set_axisbelow(True)

pyplot.subplot(2,2,3)

h, bx, by = numpy.histogram2d(tau0, tau3, bins=(20, 20))
h = numpy.ma.masked_equal(h.T, 0)
pyplot.pcolor(bx, by, h, cmap=cm.gist_earth)

if args.template_id:
    pyplot.scatter([tau0_in], [tau3_in], s=5, marker='x', color='m')
pyplot.xlabel(r"$\tau_0$")
pyplot.ylabel(r"$\tau_3$")
pyplot.xlim(0.95 * tau0.min(), 1.05 * tau0.max())
pyplot.ylim(0.95 * tau3.min(), 1.05 * tau3.max())
pyplot.locator_params(nbins=4)
pyplot.grid(True, color='0.5', linewidth=0.125, linestyle='-')
pyplot.gca().set_axisbelow(True)

#
# Density Spin based plots
#
pyplot.subplot(2,2,4)

h, bx, by = numpy.histogram2d(s1z, s2z, bins=(20, 20))
h = numpy.ma.masked_equal(h.T, 0)
pyplot.pcolor(bx, by, h, cmap=cm.gist_earth)

if args.template_id:
    pyplot.scatter([s1z_in], [s2z_in], s=5, marker='x', color='m')
pyplot.xlabel(r"$s_{1,z}$")
pyplot.ylabel(r"$s_{2,z}$")
pyplot.xlim(-1, 1)
pyplot.ylim(-1, 1)
pyplot.locator_params(nbins=4)
pyplot.grid(True, color='0.5', linewidth=0.125, linestyle='-')
pyplot.gca().set_axisbelow(True)

if args.no_title is not None and not args.no_title:
    pyplot.suptitle("Template bank density")
#pyplot.tight_layout()

fig = pyplot.gcf()
fig.subplots_adjust(wspace=0.3, hspace=0.5)
pyplot.savefig(args.plot_name_bank_density or "tmplt_bank_density.png")
