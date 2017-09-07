# Copyright (C) 2015 Reed Essick
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

## \defgroup laldetchar_py_idq_idq_summary_plots Plotting Routines
## \ingroup laldetchar_py_idq
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

plt.rcParams.update({
    'text.usetex':True,
    'font.family': 'serif',
    'font.serif': ['Computer Modern'],
    "axes.grid": True,
    "axes.axisbelow": False,
    "axes.formatter.limits": (-3, 4),
    "axes.labelsize": 20,
    "axes.titlesize": 22,
    "figure.subplot.bottom": 0.17,
    "figure.subplot.left": 0.1,
    "figure.subplot.right": 0.9,
    "figure.subplot.top": 0.90,
    "grid.color": 'gray', 
    "image.aspect": 'auto',
    "image.interpolation": 'nearest',
    "image.origin": 'lower',
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
})


from laldetchar.idq import ovl
from laldetchar.idq import idq
#from laldetchar.idq import reed as idq

from laldetchar import git_version

__author__ = \
    'Reed Essick (<reed.essick@ligo.org>)'
__version__ = git_version.id
__date__ = git_version.date

#===================================================================================================

description = """Module for plotting routines for idq pipeline."""

## \addtogroup laldetchar_py_idq_idq_summary_plots
# @{

#===================================================================================================
# defaults
#===================================================================================================

default_figsize = [10, 10]
default_axpos = [0.15, 0.15, 0.8, 0.8]

stacked_axhpos = [0.1, 0.1, 0.8, 0.4]
stacked_axcpos = [0.1, 0.5, 0.8, 0.4]

chan_perf_figsize = [20, 10]
chan_perf_axpos = [0.30, 0.15, 0.45, 0.80]

rank_timeseries_figsize = [10, 4]
#rank_timeseries_axpos = [0.10, 0.15, 0.8, 0.8]
rank_timeseries_axpos = [plt.rcParams['figure.subplot.left'], plt.rcParams['figure.subplot.bottom'], plt.rcParams['figure.subplot.right']-plt.rcParams['figure.subplot.left'], plt.rcParams['figure.subplot.top']-plt.rcParams['figure.subplot.bottom']]

rank_splittimeseries_figsize = [10, 8]
rank_splittimeseries_axpos = [plt.rcParams['figure.subplot.left'], plt.rcParams['figure.subplot.bottom'], plt.rcParams['figure.subplot.right']-plt.rcParams['figure.subplot.left'], 0.5-plt.rcParams['figure.subplot.bottom']]
fap_splittimeseries_axpos = [plt.rcParams['figure.subplot.left'], 0.52, plt.rcParams['figure.subplot.right']-plt.rcParams['figure.subplot.left'], plt.rcParams['figure.subplot.top']-0.52]

#========================
# default names
#========================

def rocfig(directory, classifier, ifo, tag, t, stride, figtype="png"):
    return "%s/%s_%s%s_ROC-%d-%d.%s"%(directory, ifo, classifier, tag, t, stride, figtype)

def histfig(directory, classifier, ifo, tag, t, stride, figtype="png"):
    return "%s/%s_%s%s_HIST-%d-%d.%s"%(directory, ifo, classifier, tag, t, stride, figtype)

def kdefig(directory, classifier, ifo, tag, t, stride, figtype="png"):
    return "%s/%s_%s%s_KDE-%d-%d.%s"%(directory, ifo, classifier, tag, t, stride, figtype)

def Lfig(directory, classifier, ifo, tag, t, stride, figtype="png"):
    return "%s/%s_%s%s_L-%d-%d.%s"%(directory, ifo, classifier, tag, t, stride, figtype)

def bitfig(directory, ifo, tag, t, stride, figtype="png"):
    return "%s/%s%s_BITWORD-%d-%d.%s"%(directory, ifo, tag, t, stride, figtype)

def ratefig(directory, ifo, tag, t, stride, figtype="png"):
    return "%s/%s%s_RATE-%d-%d.%s"%(directory, ifo, tag, t, stride, figtype)

def chanfig(directory, ifo, classifier, metric, tag, t, stride, figtype="png"):
    return "%s/%s_%s%s_chan-perf_%s-%d-%d.%s"%(directory, ifo, classifier, tag, metric, t, stride, figtype)

def calibfig(directory, ifo, classifier, tag, t, stride, figtype="png"):
    return "%s/%s_%s%s_calib-%d-%d.%s"%(directory, ifo, classifier, tag, t, stride, figtype)

def timeseriesfig(directory, ifo, classifier, tag, t, stride, figtype="png"):
    return "%s/%s_%s%s_timeseries-%d-%d.%s"%(directory, ifo, classifier, tag, t, stride, figtype)

def ovlstripchart(directory, ifo, classifier, tag, t, stride, figtype="png"):
    return "%s/%s_%s%s_chanstrip-%d-%d.%s"%(directory, ifo, classifier, tag, t, stride, figtype)

#===================================================================================================
# utility
#===================================================================================================

def close(figure):
    """ 
    straight delegation to pyplot.close 
    """
    plt.close(figure)

#===================================================================================================
# single-stride plotting functions
#===================================================================================================

def rcg_to_rocFig(c, g, color='b', label=None, figax=None):

    if figax==None:
        fig = plt.figure(figsize=default_figsize)
        ax = fig.add_axes(default_axpos)
    else:
        fig, ax = figax

    n_c = c[-1]
    n_g = g[-1]

    fap = np.array(c, dtype="float")/n_c
    eff = np.array(g, dtype="float")/n_g

    if color:
        if label:
            ax.step(fap, eff, color=color, where='post', label=label) ### plot steps after increase in eff
        else:
            ax.step(fap, eff, color=color, where='post') ### plot steps after increase in eff
        ax.plot(fap, eff, markersize=2, marker='o', markerfacecolor=color, markeredgecolor=color, linestyle='none') ### plot steps after increase in eff
    else:
        if label:
            ax.step(fap, eff, where='post', label=label) ### plot steps after increase in eff
        else:
            ax.step(fap, eff, where='post') ### plot steps after increase in eff
        ax.plot(fap, eff, markersize=2, marker='o', markerfacecolor=color, markeredgecolor=color, linestyle='none') ### plot steps after increase in eff

    ### plot error bars at sample points
    for _c, _g in zip(c, g):
        cln_CR = idq.binomialCR( _c, n_c, conf=0.68 )
        gch_CR = idq.binomialCR( _g, n_g, conf=0.68 )
        if color:
            ax.fill_between( cln_CR, 2*[gch_CR[0]], 2*[gch_CR[1]], color=color, alpha=0.05)
        else:
            ax.fill_between( cln_CR, 2*[gch_CR[0]], 2*[gch_CR[1]], alpha=0.05)

#    ### plot 90% UL on fap and eff
#    fap_UL = idq.binomialUL( c, n_c, conf=0.90 )
#    eff_UL = idq.binomialUL( g, n_g, conf=0.90 )
#    if color:
#        ax.plot( fap, eff_UL, color=color, alpha=0.25 )
#        ax.plot( fap_UL, eff, color=color, alpha=0.25 )
#    else:
#        ax.plot( fap, eff_UL, alpha=0.25 )
#        ax.plot( fap_UL, eff, alpha=0.25 )

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
        axh = fig.add_axes(stacked_axhpos)
        axc = fig.add_axes(stacked_axcpos)
    else:
        fig, axh, axc = figax

    nbins = int(max(5, np.sum(s)/nperbin))
    bins = np.linspace(bmin, bmax, nbins+1)

    if color:
        if label:
            n, b, p = axh.hist( r, bins, weights=s, color=color, label=label, linestyle=linestyle, histtype=histtype)
        else:
            n, b, p = axh.hist( r, bins, weights=s, color=color, linestyle=linestyle, histtype=histtype)
    else:
        if label:
            n, b, p = axh.hist( r, bins, weights=s, label=label, linestyle=linestyle, histtype=histtype)
        else:
            n, b, p = axh.hist( r, bins, weights=s, linestyle=linestyle, histtype=histtype)

    nsamples = 501
    rsamples = np.linspace(min(r), max(1, max(r)), nsamples)
    cum = np.zeros(nsamples)
    for _r, _s in zip(r, s):
        cum[rsamples>=_r] += _s
    cum = cum[-1] - cum ### switch so it is cumulative starting at the right
    if color:
        if label:
            axc.plot( rsamples, cum, color=color, label=label, linestyle=linestyle )
        else:
            axc.plot( rsamples, cum, color=color, linestyle=linestyle )
    else:
        if label:
            axc.plot( rsamples, cum, label=label , linestyle=linestyle )
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
        axh = fig.add_axes(stacked_axhpos)
        axc = fig.add_axes(stacked_axcpos)
    else:
        fig, axh, axc = figax

    if label:
        axh.plot( r, s, color=color, linestyle=linestyle, label=label)
    else:
        axh.plot( r, s, color=color, linestyle=linestyle)

    c = idq.kde_to_ckde( s )
    if label:
        axc.plot( r, c, color=color, linestyle=linestyle, label=label)
    else:
        axc.plot( r, c, color=color, linestyle=linestyle)

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

def stacked_L(r, c, g, color='b', linestyle='solid', label=None, figax=None):

    if not isinstance(c, np.ndarray):
        c = np.array(c)
    if not isinstance(g, np.ndarray):
        g = np.array(g)

    if figax==None:
        fig = plt.figure(figsize=default_figsize)
        axh = fig.add_axes(stacked_axhpos)
        axc = fig.add_axes(stacked_axcpos)
    else:
        fig, axh, axc = figax

    truth = c > 0
    if label:
        axh.plot( r[truth], g[truth]/c[truth], color=color, linestyle=linestyle, label=label)
    else:
        axh.plot( r[truth], g[truth]/c[truth], color=color, linestyle=linestyle)

    G = idq.kde_to_ckde( g )
    C = idq.kde_to_ckde( c )

    Truth = C > 0
    if label:
        axc.plot( r[Truth], G[Truth]/C[Truth], color=color, linestyle=linestyle, label=label)
    else:
        axc.plot( r[Truth], G[Truth]/C[Truth], color=color, linestyle=linestyle)

    axh.grid(True, which='both')
    axc.grid(True, which='both')

    axh.set_xlim(xmin=0, xmax=1)
    axc.set_xlim(xmin=0, xmax=1)

    axh.set_yscale('log')
    axc.set_yscale('log')

    if np.any(truth):
        axh.set_ylim(ymin=0, ymax=max(axh.get_ylim()[1], max(g[truth]/c[truth])*1.1))
    if np.any(Truth):
        axc.set_ylim(ymin=0, ymax=max(axc.get_ylim()[1], max(G[Truth]/C[Truth])*1.1))

    plt.setp(axc.get_xticklabels(), visible=False)
    axc.set_ylabel('cdf(g)/cdf(c)')

    axh.set_xlabel('rank')
    axh.set_ylabel('pdf(g)/pdf(c)')

    return fig, axh, axc

def stacked_params_hist( samples, figax=None, nperbin=10 , labels=None , xlabel='parameter'):

    if figax==None:
        fig = plt.figure(figsize=default_figsize)
        ax = fig.add_axes(default_axpos)
    else:
        fig, ax = figax

    nbins = int(max( 5, min([len(s)/nperbin for s in samples])))
    bins = np.linspace(np.min(samples), np.max(samples), nbins+1)

    if labels:
        ax.hist(samples, bins=bins, histtype='step', stacked=True, label=labels)
    else:
        ax.hist(samples, bins=bins, histtype='step', stacked=True)

    ax.grid(True, which="both")

    ax.set_ylabel('count')
    ax.set_xlabel(xlabel)

    return fig, ax

def bitword( samples, classifiers, figax=None, label=None):

    if figax==None:
        fig = plt.figure(figsize=default_figsize)
        ax = fig.add_axes(default_axpos)
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
        N = len(values)
        if label:
            ax.hist( values, bins, weights=np.ones(N, dtype=float)/N, histtype="step", label=label)
        else:
            ax.hist( values, bins, weights=np.ones(N, dtype=float)/N, histtype="step")

#    ax.set_ylabel('count')
    ax.set_ylabel('fraction of identified events')

    ax.xaxis.set_ticks(ticks)
    ax.xaxis.set_ticklabels([l.strip("\n") for l in ticklabels])

    ax.set_xlim(xmin=bins[0], xmax=bins[-1])
    ax.set_ylim(ymin=0, ymax=1)

    return fig, ax 

def calibration_scatter( deadtimes, statedFAPs, color='b', label=None, marker='o', figax=None):

    if figax==None:
        fig = plt.figure(figsize=default_figsize)
        ax = fig.add_axes(default_axpos)
    else:
        fig, ax = figax

    if color:
        if label:
            ax.plot( statedFAPs, deadtimes, label=label, linestyle='none', marker=marker, markeredgecolor=color, markerfacecolor='none' ) 
        else:
            ax.plot( statedFAPs, deadtimes, linestyle='none', marker=marker, markeredgecolor=color, markerfacecolor='none' ) 
    else:
        if label:
            ax.plot( statedFAPs, deadtimes, label=label, linestyle='none', marker=marker, markerfacecolor='none' ) 
        else:
            ax.plot( statedFAPs, deadtimes, linestyle='none', marker=marker, markerfacecolor='none' ) 

    ax.set_xlabel( 'stated False Alarm Probability' )
    ax.set_ylabel( 'observed deadtime' )

    ax.grid(True, which='both')

#    ax.set_xscale('log')
#    ax.set_yscale('log')
#    ax.set_xlim(xmin=1e-3, xmax=1)
#    ax.set_ylim(ymin=1e-3, ymax=1)

    ax.set_xscale('linear')
    ax.set_yscale('linear')
    ax.set_xlim(xmin=0, xmax=1)
    ax.set_ylim(ymin=0, ymax=1)

    return fig, ax

#===================================================================================================
# trending (multi-stride) plotting functions
#===================================================================================================

def rates( ranges, values, color='b', label=None, figax=None, linestyle='solid'):

    if figax==None:
        fig = plt.figure(figsize=default_figsize)
        ax = fig.add_axes(default_axpos)
    else:
        fig, ax = figax

    S = int(np.max(ranges))

    for v, (s, e) in zip(values, ranges):
        if color:
            ax.plot((s-S, e-S), [v]*2, color=color, linestyle=linestyle, alpha=0.25)
        else:
            ax.plot((s-S, e-S), [v]*2, linestyle=linestyle, alpha=0.25)

    if color:
         if label:
             ax.plot( [0.5*(s+e)-S for s, e in ranges], values, color=color, label=label, linestyle=linestyle, marker='o', markerfacecolor=color, markeredgecolor=color, markersize=4)
         else:
             ax.plot( [0.5*(s+e)-S for s, e in ranges], values, color=color, linestyle=linestyle, marker='o', markerfacecolor=color, markeredgecolor=color, markersize=4)
    else:
         if label:
             ax.plot( [0.5*(s+e)-S for s, e in ranges], values, label=label, linestyle=linestyle, marker='o', markersize=4)
         else:
             ax.plot( [0.5*(s+e)-S for s, e in ranges], values, linestyle=linestyle, marker='o', markersize=4)

    ax.grid(True, which="both")

    ax.set_xlabel('sec relative to %d'%S)

    return fig, ax

def channel_performance( vchans, ranges, performances, num_gchs, num_clns, figax=None, cmap=None):
    
    if figax == None:
        figrank = plt.figure(figsize=chan_perf_figsize)
        axrank = figrank.add_axes(chan_perf_axpos)

        figeff = plt.figure(figsize=chan_perf_figsize)
        axeff = figeff.add_axes(chan_perf_axpos)

        figfap = plt.figure(figsize=chan_perf_figsize)
        axfap = figfap.add_axes(chan_perf_axpos)
    else:
        (figrank, axrank), (figeff, axeff), (figfap, axfap) = figax

    if cmap==None:
        cmap = matplotlib.cm.get_cmap()

    S = int(np.max(ranges))
    ind = 0
    for ind, vchan in enumerate(vchans):
        for (s, e), performance, num_gch, num_cln in zip(ranges, performances, num_gchs, num_clns):
            if performance.has_key(vchan):
                gch = performance[vchan]['gch']
                cln = performance[vchan]['cln']

                if num_gch:
                    eff = 1.0*gch/num_gch
                else:
                    eff = 0

                if num_cln:
                    fap = 1.0*cln/num_cln
                else:
                    fap = 0

                if fap > 0:
                    rank = ovl.effbydt_to_rank( eff/fap )
                elif eff > 0:
                    rank = 1
                else:
                    rank = 0
                
                axrank.fill_between( [s-S, e-S], [ind-0.5]*2, [ind+0.5]*2, color=cmap(rank) )
                axeff.fill_between( [s-S, e-S], [ind-0.5]*2, [ind+0.5]*2, color=cmap(eff) )
                axfap.fill_between( [s-S, e-S], [ind-0.5]*2, [ind+0.5]*2, color=cmap(fap) )

    ticks = range(ind+1)
    for ax in [axrank, axeff, axfap]:
        ax.yaxis.set_ticks( ticks )
        ax.yaxis.set_ticklabels( [vchan.replace("_","\_") for vchan in vchans] )
       
        ax.set_xlabel('sec relative to %d'%S)
       
        ax.set_ylim(ymin=-0.5, ymax=ind+0.5)
 
    return (figrank, axrank), (figeff, axeff), (figfap, axfap)

#===================================================================================================
# Plots for GraceDB follow-up

def rank_timeseries( times, timeseries, start=None, end=None, to=None, shade_alpha=0.25, color='r', linestyle='solid', figax=None, zorder=1, alpha=1.0, linewidth=1.0 ):
    """
    makes a timeseries plot for GraceDB follow-up annotations
    """
    if figax==None:
        fig = plt.figure(figsize=rank_timeseries_figsize)
        ax = fig.add_axes(rank_timeseries_axpos)
    else:
        fig, ax = figax

    if start==None:
        start = times[0][0]
    if end==None:
        end = times[-1][-1]
    if to==None:
        to = start

    this_start = start
    for t, ts in zip(times, timeseries):
        if shade_alpha and (t[0]!=this_start):
            ax.fill_between(np.array([this_start, t[0]])-to, [1e-5,1e-5], [1,1], color='k', alpha=shade_alpha, linewidth=0.0 )
        ax.plot( np.array(t)-to, ts, color=color, linestyle=linestyle, zorder=zorder, alpha=alpha, linewidth=linewidth )
        this_start = t[-1] + (t[1]-t[0])

    return fig, ax

def rank_timeseries_filled( times, timeseries, baseline=1, start=None, end=None, to=None, shade_alpha=0.25, color='r', linestyle='solid', figax=None, zorder=1, alpha=1.0, linewidth=0.0 ):
    """
    makes a timeseries plot for GraceDB follow-up annotations
    """
    if figax==None:
        fig = plt.figure(figsize=rank_timeseries_figsize)
        ax = fig.add_axes(rank_timeseries_axpos)
    else:
        fig, ax = figax

    if start==None:
        start = times[0][0]
    if end==None:
        end = times[-1][-1]
    if to==None:
        to = start

    this_start = start
    for t, ts in zip(times, timeseries):
        if shade_alpha and (t[0]!=this_start):
            ax.fill_between(np.array([this_start, t[0]])-to, [1e-5,1e-5], [1,1], color='k', alpha=shade_alpha, linewidth=0.0 )

        ts[ts!=ts] = 1 ### get rid of Nan, in case they exist
                       ### this is hacky, and we shoudl remove it if possible
                       ### Nan is likely produced within idq.BinomialUL, and may be fixed but we'll keep this here just in case

        ax.fill_between( np.array(t)-to, ts, np.ones(np.shape(ts))*baseline, color=color, alpha=alpha, edgecolor=color, linewidth=linewidth )
        this_start = t[-1] + (t[1]-t[0])

    return fig, ax

def ovl_vconfig_stripchart( channel_ordered_vconfigs, to=None, cmap='jet', mapIt=dict( (k,i) for i,k in enumerate("vchan vthr vwin metric metric_exp rank".split()) ) ):
    """
    builds a strip chart based on the info in channel_ordered_vconfigs

    WARNING: 
        you may get gaps in your strip chart corresponding to edges of frames.
        this is a plotting artefact, and nothing more
    """
    if to==None:
        to = min([values[0][0] for value in channel_ordered_vconfigs.values()])

    ### find configs contained
    configs = set()
    for c in channel_ordered_vconfigs.keys():
        for C in c:
            configs.add( C )
    configs = list( configs )

    try:
        configs.remove( (None, None, None, None, 0, 0) )
    except ValueError as e: ### None is not present...
        pass

    configs.sort( key=lambda l: l[mapIt['vwin']], reverse=True )
    configs.sort( key=lambda l: l[mapIt['vthr']], reverse=True )
    configs.sort( key=lambda l: l[mapIt['vchan']] )

    configs.insert(0, (None, None, None, None, 0, 0) )
    
    N = len(configs)

    ### figure out the size of the axis object and the figure
    chan_stripchart_left = 0.3
    chan_stripchart_right = 0.90
    chan_stripchart_figwidth = 10

    chan_stripchart_topmargin = 0.4
    chan_stripchart_botmargin = 0.68
    chan_stripchart_chanheight = 0.95

    chan_stripchart_figheight = chan_stripchart_topmargin + N*chan_stripchart_chanheight + chan_stripchart_botmargin

    chan_stripchart_axpos = [chan_stripchart_left, chan_stripchart_botmargin/chan_stripchart_figheight, chan_stripchart_right-chan_stripchart_left, 1-chan_stripchart_topmargin/chan_stripchart_figheight - chan_stripchart_botmargin/chan_stripchart_figheight]
#    chan_stripchart_axpos = [chan_stripchart_left, plt.rcParams['figure.subplot.bottom'], plt.rcParams['figure.subplot.right']-chan_stripchart_left, plt.rcParams['figure.subplot.top']-plt.rcParams['figure.subplot.bottom']]

    chan_stripchart_cb_left = chan_stripchart_right + 0.01
    chan_stripchart_cb_right = 0.93

    chan_stripchart_cbpos = [chan_stripchart_cb_left, chan_stripchart_axpos[1], chan_stripchart_cb_right-chan_stripchart_cb_left, chan_stripchart_axpos[3]]

    fig = plt.figure(figsize=[chan_stripchart_figwidth, chan_stripchart_figheight])
    ax = fig.add_axes(chan_stripchart_axpos)

    cbax = fig.add_axes(chan_stripchart_cbpos)

    ymap = dict( (config, ind) for ind, config in enumerate(configs) )

    cmap = plt.get_cmap(cmap)

    ### plot segments for each vconfig
    alpha = 1.0
    for cnfgs, segs in channel_ordered_vconfigs.items():
        for config in cnfgs:
            y = ymap[config]

            metric_exp = config[mapIt['metric_exp']]
     
            color = cmap( config[mapIt['rank']] )

            for start, end in segs:
                ax.fill_between( [start-to, end-to], [y+0.02]*2, [y+0.98]*2, edgecolor=color, color=color, linewidth=0.5, alpha=alpha ) 

    ### plot the colorbar
    cmap_norm=matplotlib.colors.Normalize( vmin=0.0, vmax=1.0 )
    cb = matplotlib.colorbar.ColorbarBase( cbax, cmap=cmap, norm=cmap_norm, orientation='vertical', alpha=alpha )

    ### decorate the plot
    ax.set_xlabel('Time [sec] since %.3f'%to)
    ax.set_yticks( [ind+0.5 for ind in range(len(configs))])
    yticklabels = ["None"]
    for c in configs[1:]:
        vchan = c[mapIt['vchan']]
        try:
            vchan = vchan.split("_")
            flow, fhigh = [float(l) for l in vchan[-2:]]
            vchan = "%s:%s"%(vchan[0], "\_".join(vchan[1:-2]))
            yticklabels.append( "%s\nmin freq : %.1f Hz\nmax freq : %.1f Hz\nthreshold : %.1f\nwindow : %.3f sec"%(vchan, flow, fhigh, c[mapIt['vthr']], c[mapIt['vwin']]) )
        except:
            yticklabels.append( "%s\nthreshold : %.1f\nwindow : %.3f sec"%(vchan.replace("_","\_"), c[mapIt['vthr']], c[mapIt['vwin']]) )
    ax.set_yticklabels( yticklabels, fontsize=10)

    ax.set_ylim(ymin=0, ymax=N)

    for tic in ax.yaxis.get_major_ticks():
         tic.tick1On = tic.tick2On = False
#         tic.label1On = tic.label2On = False

    ax.grid(False, which="both")

    ### decorate colorbar
    if N < 3:
        cb.set_ticks([0.0, 0.5, 1.0])
    elif N < 5:
        cb.set_ticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

    return fig, ax, cbax

##@}
