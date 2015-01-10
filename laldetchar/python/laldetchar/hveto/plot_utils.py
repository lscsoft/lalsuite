# Copyright (C) 2013 Chris Pankow
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation,.

## \addtogroup laldetchar_py_hveto
"""Some plotting utilities for HVeto based analyses."""
# \author Chris Pankow (<chris.pankow@ligo.org>)
# ### Synopsis ###
# ~~~
# from laldetchar.hveto import plot_utils
# ~~~

import sys

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot
# FIXME: gridspec...
#import breakout

# SWIG bindings
import laldetchar
import lal
# FIXME: Will this be needed with Karl's bindings?
# Yes... lalburst doesn't import the wrapped version of SnglBurst, so trying to
# get its attributes will fail until they are defined with this
# But that's been fixed, so we should try it without
import lalmetaio

from laldetchar import git_version as version
__author__ = "Chris Pankow <chris.pankow@ligo.org>"
__version__ = version.id
__date__ = version.date

## \addtogroup laldetchar_py_hveto_plot_utils
#@{

def plot_scatter(trigs, ref_chan, saveto, runseg):
	"""
	Make a colorized scatter plot of triggers from trig from channel ref_chan. Save this to file name saveto. Optionally, reference times to st.
	"""

	st, end = runseg[0], runseg[1]
	st_h, end_h = 0, (end-st)/3600.0

	pyplot.figure()
	t, f, s = zip(*[ [float(sb.get_peak()), sb.central_freq, math.log10(sb.confidence)] for sb in trigs if sb.channel == ref_chan ])
	t = [ (ti-st)/3600.0 for ti in t ]
	pyplot.xlabel("Hours since %d" % st)
	pyplot.ylabel("Central Frequency (Hz)")
	pyplot.scatter( t, f, c=s )
	pyplot.xlim([st_h, end_h])
	cbar = pyplot.colorbar()
	cbar.set_label("log10 Significance")
	pyplot.grid()
	pyplot.savefig(saveto)

def plot_sigdrop(chansig_prev, chansig_cur, rnd):
	"""
	Make a significance 'drop plot', where the significance from several channels is compared in difference across two rounds.
	"""
	fig = pyplot.figure(figsize=(30,10))
	ax = fig.add_subplot(111)
	nchan = 1
	for chan, sig in sorted(chansig_prev.iteritems()):
		try:
			sprev, snew = sig, chansig_cur[chan]
		except KeyError:
			chansig_cur[chan] = 0
			sprev, snew = sig, 0
		pyplot.arrow( nchan, sprev, 0, (snew-sprev), fc="k", head_width=1, head_length=0.1 )
		nchan += 1

	pyplot.title("Significance drop plot round %d" % (rnd-1))
	pyplot.ylabel("Significance")
	pyplot.xlim(0, nchan+1)
	pyplot.ylim(0, 1.1*max(chansig_prev.values()))
	pyplot.grid()

	ax.set_xticks(range(1, nchan))
	ax.set_xticklabels(sorted(chansig_prev.keys()), rotation=270, fontsize=6 , weight="light")
	fig.savefig("drop_round_%d.png" % (rnd-1))

# close doxygen
##@}
