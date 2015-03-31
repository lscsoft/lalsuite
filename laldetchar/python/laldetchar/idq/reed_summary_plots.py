# Copyright (C) 2013 Reed Essick, Ruslan Vaulin
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your# option
# ) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

## \addtogroup laldetchar_py_idq

## Synopsis
# ~~~
# from laldetchar.idq import idq_summary_plots
# ~~~
# \author Reed Essick (<reed.essick@ligo.org>), Ruslan Vaulin (<ruslan.vaulin@ligo.org>)

#=================================================

import sys

import numpy as np

from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.cm
from matplotlib import mpl
from matplotlib import pyplot as plt

plt.rcParams.update( {'text.usetex':True} )

from laldetchar.idq import pdf_estimation as pdf_e
from laldetchar.idq import ovl
from laldetchar.idq import reed

from laldetchar import git_version

__author__ = \
    'Reed Essick (<reed.essick@ligo.org>), Ruslan Vaulin <ruslan.vaulin@ligo.org>'
__version__ = git_version.id
__date__ = git_version.date

#===================================================================================================

description = """Module for plotting routines for idq pipeline."""

#===================================================================================================
# defaults
#===================================================================================================

default_figsize = [10, 10]
default_axpos = [0.15, 0.15, 0.8, 0.8]

#========================
# default names
#========================

def rocfig(directory, classifier, ifo, tag, t, stride, figtype="png"):
    return "%s/%s_%s%s_ROC-%d-%d.%s"%(directory, ifo, classifier, tag, t, stride, figtype)

def histfig(directory, classifier, ifo, tag, t, stride, figtype="png"):
    return "%s/%s_%s%s_HIST-%d-%d.%s"%(directory, ifo, classifier, tag, t, stride, figtype)

def kdefig(directory, classifier, ifo, tag, t, stride, figtype="png"):
    return "%s/%s_%s%s_KDE-%d-%d.%s"%(directory, ifo, classifier, tag, t, stride, figtype)

def kdename(directory, classifier, ifo, tag, t, stride):
    return "%s/%s_%s%s_KDE-%d-%d.npy.gz"%(directory, ifo, classifier, tag, t, stride)

def bitfig(directory, ifo, tag, t, stride, figtype="png"):
    return "%s/%s%s_BITWORD-%d-%d.%s"%(directory, ifo, tag, t, stride, figtype)

#===================================================================================================
# utility
#===================================================================================================

def close(figure):
    """ 
    straight delegation to pyplot.close 
    """
    plt.close(figure)

def kde_pwg(eval, r, ds, scale=0.1, s=0.01):
    """
    delegates to pdf_e.point_wise_gaussian_kde after constructing the correct input values
    """
    observ = []
    for R, dS in zip(r, ds):
        observ += [R]*dS

    if np.sum(ds):
        return pdf_e.point_wise_gaussian_kde(eval, observ, scale=scale, s=s)
    else:
        return np.ones_like(eval)


def bin_by_rankthr(rankthr, output, columns=None):

    if columns==None:
        columns = sorted(output.keys())

    gch = []
    cln = []
    for i in xrange(len(output['i'])):
        if output['rank'][i] >= rankthr:
            d = dict( (column, output[column][i]) for column in columns )
            if output['i'][i]:
                gch.append( d )
            else:
                cln.append( d )

    return cln, gch

#===================================================================================================
# single-stride plotting functions
#===================================================================================================

def rcg_to_rocFig(c, g, color='b', label=None, figax=None):

    if figax==None:
        fig = plt.figure(figsize=default_figsize)
        ax = plt.subplot(1,1,1)
    else:
        fig, ax = figax

    n_c = c[-1]
    n_g = g[-1]

    fap = np.array(c, dtype="float")/c[-1]
    eff = np.array(g, dtype="float")/g[-1]

    if color:
        ax.step(fap, eff, color=color, where='post', label=label) ### plot steps after increase in eff
        ax.plot(fap, eff, markersize=2, marker='o', markerfacecolor=color, markeredgecolor=color, linestyle='none') ### plot steps after increase in eff
    else:
        ax.step(fap, eff, where='post', label=label) ### plot steps after increase in eff
        ax.plot(fap, eff, markersize=2, marker='o', markerfacecolor=color, markeredgecolor=color, linestyle='none') ### plot steps after increase in eff

    ### plot error bars at sample points
    for _c, _g in zip(c, g):
        cln_CR = reed.binomialCR( _c, n_c, conf=0.68 )
        gch_CR = reed.binomialCR( _g, n_g, conf=0.68 )
        if color:
            ax.fill_between( cln_CR, 2*[gch_CR[0]], 2*[gch_CR[1]], color=color, alpha=0.05)
        else:
            ax.fill_between( cln_CR, 2*[gch_CR[0]], 2*[gch_CR[1]], alpha=0.05)

    ax.grid(True, which='both')

    ax.set_xscale('log')
    ax.set_yscale('linear')

    ax.set_xlim(xmin=1e-5, xmax=1)
    ax.set_ylim(ymin=0, ymax=1)

    ax.set_xlabel('False Alarm Probability')
    ax.set_ylabel('Glitch Detection Efficiency')

    return fig, ax

def stacked_hist(r, s, color='b', linestyle='solid', histtype='step', label=None, figax=None, nperbin=10, bmin=0, bmax=1):

    if figax==None:
        fig = plt.figure(figsize=default_figsize)
        axh = fig.add_axes([0.1, 0.1, 0.8, 0.4])
        axc = fig.add_axes([0.1, 0.5, 0.8, 0.4])
    else:
        fig, axh, axc = figax

    nbins = int(max(5, np.sum(s)/nperbin))
    bins = np.linspace(bmin, bmax, nbins+1)

    if color:
        n, b, p = axh.hist( r, bins, weights=s, color=color, label=label , linestyle=linestyle, histtype=histtype)
    else:
        n, b, p = axh.hist( r, bins, weights=s, label=label , linestyle=linestyle, histtype=histtype)

    nsamples = 501
    rsamples = np.linspace(min(r), max(1, max(r)), nsamples)
    cum = np.zeros(nsamples)
    for _r, _s in zip(r, s):
        cum[rsamples>=_r] += _s
    cum = cum[-1] - cum ### switch so it is cumulative starting at the right
    if color:
        axc.plot( rsamples, cum, color=color, label=label , linestyle=linestyle )
    else:
        axc.plot( rsamples, cum, label=label , linestyle=linestyle )

    axh.grid(True, which="both")
    axc.grid(True, which="both")

    axh.set_ylim(ymin=0, ymax=max(n)*1.1)
    axc.set_ylim(ymin=0, ymax=np.sum(s))

    plt.setp(axc.get_xticklabels(), visible=False)
    axc.set_ylabel('cumulative count')

    axh.set_xlabel('rank')
    axh.set_ylabel('count')

    return fig, axh, axc

def stacked_kde(r, s, color='b', linestyle='solid', label=None, figax=None):

    if figax==None:
        fig = plt.figure(figsize=default_figsize)
        axh = fig.add_axes([0.1, 0.1, 0.8, 0.4])
        axc = fig.add_axes([0.1, 0.5, 0.8, 0.4])
    else:
        fig, axh, axc = figax

    axh.plot( r, s, color=color, linestyle=linestyle, label=label)

    c = np.zeros_like(r)
    _s = 0.0
    cum = 0.0
    for i, S in enumerate(s[::-1]): ### integrate by hand, and start at the right
        cum += 0.5*(S+_s)
        c[i] = cum
        _s = S
    c /= c[-1]
    c = c[::-1] ### reverse order to match r
    axc.plot( r, c, color=color, linestyle=linestyle, label=label)

    axh.grid(True, which="both")
    axc.grid(True, which="both")

    axh.set_xlim(xmin=0, xmax=1)
    axc.set_xlim(xmin=0, xmax=1)

    axh.set_ylim(ymin=0, ymax=max(axh.get_ylim()[1], max(s)*1.1) )
    axc.set_ylim(ymin=0, ymax=1.0)

    plt.setp(axc.get_xticklabels(), visible=False)
    axc.set_ylabel('cumulative fraction of events')

    axh.set_xlabel('rank')
    axh.set_ylabel('p(rank)')

    return fig, axh, axc

def stacked_params_hist( samples, figax=None, nperbin=10 , labels=None , xlabel='parameter'):

    if figax==None:
        fig = plt.figure(figsize=default_figsize)
        ax = plt.subplot(1,1,1)
    else:
        fig, ax = figax

    nbins = int(max( 5, min([len(s)/nperbin for s in samples])))
    bins = np.linspace(np.min(samples), np.max(samples), nbins+1)

    ax.hist(samples, bins=bins, histtype='step', stacked=True, label=labels)

    ax.grid(True, which="both")

    ax.set_ylabel('count')
    ax.set_xlabel(xlabel)

    return fig, ax

def bitword( samples, classifiers, figax=None, label=None):

    if figax==None:
        fig = plt.figure(figsize=default_figsize)
        ax = plt.subplot(1,1,1)
    else:
        fig, ax = figax
   
    ### make associations by gps times
    gps = defaultdict( int ) ### get all gps times
    ticks = [0]
    ticklabels = [""]
    for ind, classifier in enumerate(classifiers):
        dint = 2**ind ### value for by which we increment gps's int to signify detection by this classifier

        ticks += [t+dint for t in ticks] ### set up ticks and tick labels
#        ticklabels = ["%s\n!%s"%(l, classifier) for l in ticklabels] + ["%s\n%s"%(l, classifier) for l in ticklabels]
        ticklabels += ["%s\n%s"%(l, classifier) for l in ticklabels]

        for output in samples[classifier]:
            gps[ output['GPS'] ] += dint
    
#    bins = np.arange(2**len(classifiers)+1)-0.5
    bins = np.arange(2**len(classifiers))+0.5
    values = gps.values()
    if values: 
#        ax.hist( values, bins, histtype='step', label=label )
        N = len(values)
        ax.hist( values, bins, weights=np.ones(N, dtype=float)/N, histtype="step", label=label )

#    ax.set_ylabel('count')
    ax.set_ylabel('fraction of identified events')

    ax.xaxis.set_ticks(ticks)
    ax.xaxis.set_ticklabels([l.strip("\n") for l in ticklabels])

    ax.set_xlim(xmin=bins[0], xmax=bins[-1])
    ax.set_ylim(ymin=0, ymax=1)

    return fig, ax 

#===================================================================================================
# trending (multi-stride) plotting functions
#===================================================================================================

def rates( ranges, values, color='b', label=None, figax=None):

    if figax==None:
        fig = plt.figure()
        ax = plt.subplot(1,1,1)
    else:
        fig, ax = figax

    for v, s_e in zip(values, ranges):
        if color:
            ax.plot(s_e, [v]*2, color=color, label=label)
        else:
            ax.plot(s_e, [v]*2, label=label)
        label=None

    ax.grid(True, which="both")

    return fig, ax

