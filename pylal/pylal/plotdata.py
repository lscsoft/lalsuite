#!/usr/bin/env python

# Copyright (C) 2011 Duncan Macleod
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
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
This module provides some wrapper functions for plotting TimeSeries and FrequencySeries using classes defined in pylal.plotutils.
"""

# =============================================================================
# Preamble
# =============================================================================

from __future__ import division

import datetime
import re
import itertools
import numpy
import numbers

from scipy import interpolate

import lal

from pylal import plotutils
from pylal import seriesutils
from pylal import git_version
from pylal import rate

__author__  = "Duncan M. Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__    = git_version.date

# =============================================================================
# Plot TimeSeries
# =============================================================================

def plottimeseries(series, outfile, t0=0, zeroline=False, **kwargs):
    """
    Plot a REALXTimeSeries.
    """
    # construct list of series
    if hasattr(series, "__contains__"):
        serieslist = series
    else:
        serieslist = [series]

    # get limits
    xlim = kwargs.pop("xlim", None)
    ylim = kwargs.pop("ylim", None)
    if xlim:
        start,end = xlim
    else:
        start = min(float(s.epoch) for s in serieslist)
        end   = max(float(s.epoch + s.data.length*s.deltaT) for s in serieslist)

    # get axis scales
    logx = kwargs.pop("logx", False)
    logy = kwargs.pop("logy", False)

    # get legend loc
    loc = kwargs.pop("loc", 0)
    alpha = kwargs.pop("alpha", 0.8)

    # get colorbar options
    hidden_colorbar = kwargs.pop("hidden_colorbar", False)

    # get savefig option
    bbox_inches = kwargs.pop("bbox_inches", None)

    #
    # get labels
    #

    xlabel = kwargs.pop("xlabel", None)
    if xlabel:
        unit = 1
    if not xlabel:
        unit, timestr = plotutils.time_axis_unit(end-start)
        if not t0:
            t0 = start
        t0 = lal.LIGOTimeGPS(t0)
        if int(t0.gpsNanoSeconds) == 0:
            xlabel = datetime.datetime(*lal.GPSToUTC(int(t0))[:6])\
                         .strftime("%B %d %Y, %H:%M:%S %ZUTC")
            xlabel = "Time (%s) since %s (%s)" % (timestr, xlabel, int(t0))
        else:
            xlabel = datetime.datetime(*lal.GPSToUTC(t0.gpsSeconds)[:6])\
                          .strftime("%B %d %Y, %H:%M:%S %ZUTC")
            xlabel = "Time (%s) since %s (%s)"\
                     % (timestr,
                        xlabel.replace(" UTC",".%.3s UTC" % t0.gpsNanoSeconds),\
                        t0)
        t0 = float(t0)
    ylabel   = kwargs.pop("ylabel", "Amplitude (counts)")
    title    = kwargs.pop("title", "")
    subtitle = kwargs.pop("subtitle", "")

    #
    # set default line params
    #

    color = kwargs.pop("color", None)
    if isinstance(color, str):
        color = color.split(',')
        if len(color) == 1:
            color = [color[0]]*len(serieslist)
    kwargs2 = dict()
    kwargs2.update(kwargs)
    if kwargs.has_key("marker"):
        kwargs.setdefault("markersize", 5)
        kwargs2.setdefault("markersize", 2)
    else:
        kwargs.setdefault("linestyle", "-")
        kwargs.setdefault("linewidth", "1")
        kwargs2.setdefault("linewidth", "0.1")

    #
    # make plot
    #

    allnames    = [s.name for s in serieslist]
    namedseries = [s for s in serieslist if not re.search("(min|max)\Z",s.name)]

    plot = plotutils.SimplePlot(xlabel, ylabel, title, subtitle)
    for i,(series,c) in enumerate(itertools.izip(namedseries,\
                                                 plotutils.default_colors())):
        if color:
            c = color[serieslist.index(series)]
        x = numpy.arange(series.data.length) * series.deltaT +\
            float(series.epoch) - float(t0)
        x = x.astype(float)
        x /= unit
        d = series.data.data
        if logy and ylim:
            numpy.putmask(d, d==0, ylim[0]-abs(ylim[0])*0.01)
        plot.add_content(x, d, color=c,\
                         label=plotutils.display_name(series.name), **kwargs)
        # find min/max and plot
        for i,name in enumerate(allnames):
            for ext in ["min", "max"]:
                if re.match("%s[- _]%s" % (re.escape(series.name), ext), name):
                    if color:
                        c = color[i]
                    series2 = serieslist[i]
                    x2 = numpy.arange(series2.data.length) * series2.deltaT\
                         + float(series2.epoch) - float(t0)
                    x2 /= unit
                    d2 = series2.data.data
                    if logy and ylim:
                        numpy.putmask(d2, d2==0, ylim[0]-abs(ylim[0])*0.01)
                    plot.ax.plot(x2, d2, color=c, **kwargs2)
                    plot.ax.fill_between(x2, d, d2, alpha=0.1,\
                                         color=c)

    # finalize
    plot.finalize(loc=loc, alpha=alpha)
    if hidden_colorbar:
        plotutils.add_colorbar(plot.ax, visible=False)
    
    # add zero line
    axis_lims = plot.ax.get_ylim()
    if zeroline:
        plot.ax.plot([0, 0], [axis_lims[0], axis_lims[1]], 'r--', linewidth=2)
        plot.ax.set_ylim([ axis_lims[0], axis_lims[1] ])

    # set logscale
    if logx:
        plot.ax.xaxis.set_scale("log")
    if logy:
        plot.ax.yaxis.set_scale("log")
    plot.ax._update_transScale()

    # format axes
    if xlim:
        xlim = (numpy.asarray(xlim).astype(float)-t0)/unit
        plot.ax.set_xlim(xlim)
    if ylim:
        plot.ax.set_ylim(ylim)

    plotutils.set_time_ticks(plot.ax)
    plotutils.set_minor_ticks(plot.ax, x=False)

    # save and close
    plot.savefig(outfile, bbox_inches=bbox_inches,\
                 bbox_extra_artists=plot.ax.texts)
    plot.close()

# =============================================================================
# Plot FrequencySeries
# =============================================================================

def plotfrequencyseries(series, outfile, **kwargs):
    """
    Plot a (swig)LAL FrequencySeries.
    """
    # construct list of series
    if hasattr(series, "__contains__"):
        serieslist = series
    else:
        serieslist = [series]

    # get limits
    xlim = kwargs.pop("xlim", None)
    ylim = kwargs.pop("ylim", None)
    if not xlim:
        fmin = min(float(s.f0) for s in serieslist)
        fmax = max(float(s.f0 + s.data.length*s.deltaF) for s in serieslist)
        xlim = [fmin, fmax]
        
    # get axis scales
    logx = kwargs.pop("logx", False)
    logy = kwargs.pop("logy", False)

    # get legend loc
    loc = kwargs.pop("loc", 0)
    alpha = kwargs.pop("alpha", 0.8)

    # get colorbar options
    hidden_colorbar = kwargs.pop("hidden_colorbar", False)

    # get savefig option
    bbox_inches = kwargs.pop("bbox_inches", None)

    #
    # get labels
    #

    xlabel = kwargs.pop("xlabel", "Frequency (Hz)")
    ylabel   = kwargs.pop("ylabel", "Amplitude")
    title    = kwargs.pop("title", "")
    subtitle = kwargs.pop("subtitle", "")

    #
    # set default line params
    #

    color = kwargs.pop("color", None)
    if isinstance(color, str):
        color = [color]*len(serieslist)
    kwargs2 = dict()
    kwargs2.update(kwargs)
    if kwargs.has_key("marker"):
        kwargs.setdefault("markersize", 5)
        kwargs2.setdefault("markersize", 2)
    else:
        kwargs.setdefault("linestyle", "-")
        kwargs.setdefault("linewidth", "1")
        kwargs2.setdefault("linewidth", "0.1")

    #
    # make plot
    #

    allnames    = [s.name for s in serieslist]
    namedseries = [s for s in serieslist if not re.search("(min|max)\Z",s.name)]

    plot = plotutils.SimplePlot(xlabel, ylabel, title, subtitle)
    for i,(series,c) in enumerate(itertools.izip(namedseries,\
                                                 plotutils.default_colors())):
        if color:
            c = color[serieslist.index(series)]
        if series.f_array is not None:
            x = series.f_array
        else:
            x = numpy.arange(series.data.length) * series.deltaF + series.f0 
        x = x.astype(float)
        plot.add_content(x, series.data.data, color=c,
                         label=plotutils.display_name(series.name), **kwargs)
        # find min/max and plot
        for i,name in enumerate(allnames):
            for ext in ["min", "max"]:
                if re.match("%s[- _]%s" % (re.escape(series.name), ext), name):
                    if color:
                        c = color[i]
                    series2 = serieslist[i]
                    if series2.f_array is not None:
                        x2 = series2.f_array
                    else:
                        x2 = numpy.arange(series2.data.length) * series2.deltaF\
                            + series2.f0 
                    plot.ax.plot(x2, series2.data.data, color=c, **kwargs2)
                    # sanity check for malformed inputs
                    if series.data.data.shape == series2.data.data.shape:
                        plot.ax.fill_between(x2, series.data.data,\
                                                 series2.data.data, alpha=0.1,\
                                                 color=c)

    # finalize
    plot.finalize(loc=loc, alpha=alpha)
    if hidden_colorbar:
        plotutils.add_colorbar(plot.ax, visible=False)
    
    # set logscale
    if logx:
        plot.ax.xaxis.set_scale("log")
    if logy:
        plot.ax.yaxis.set_scale("log")
    plot.ax._update_transScale()

    # format axes
    if xlim:
        plot.ax.set_xlim(xlim)
    if ylim:
        plot.ax.set_ylim(ylim)

    plot.ax.grid(True, which="both")
    plotutils.set_minor_ticks(plot.ax)

    # save and close
    plot.savefig(outfile, bbox_inches=bbox_inches,\
                 bbox_extra_artists=plot.ax.texts)
    plot.close()

# =============================================================================
# Plot spectrogram
# =============================================================================

def plotspectrogram(sequencelist, outfile, epoch=0, deltaT=1, f0=0, deltaF=1,\
                    t0=0, ydata=None, **kwargs):
    """
    Plots a list of REAL?VectorSequences on a time-frequency-amplitude colour
    map. The epochand deltaT arguments define the spacing in the x direction
    for the given VectorSequence, and similarly f0 and deltaF in the
    y-direction. If a list of VectorSequences is given, any of these arguments
    can be in list form, one for each sequence.

    ydata can be given as to explicitly define the frequencies at which the
    sequences are sampled.
    """
    # construct list of series
    if not hasattr(sequencelist, "__contains__"):
        sequencelist = [sequencelist]
    numseq = len(sequencelist)

    # format variables
    if isinstance(epoch, numbers.Number) or isinstance(epoch, lal.LIGOTimeGPS):
        epoch = [epoch]*numseq
        epoch = map(float, epoch)
    if not len(epoch) == numseq:
        raise ValueError("Wrong number of epoch arguments given.")
    if isinstance(deltaT, numbers.Number):
        deltaT = [deltaT]*numseq
        deltaT = map(float, deltaT)
    if not len(deltaT) == numseq:
        raise ValueError("Wrong number of deltaT arguments given.")
    if not ydata is None:
        if isinstance(f0, numbers.Number):
            f0 = [f0]*numseq
            f0 = map(float, f0)
        if not len(f0) == numseq:
            raise ValueError("Wrong number of f0 arguments given.")
        if isinstance(deltaF, numbers.Number):
            deltaF = [deltaF]*numseq
            deltaF = map(float, deltaF)
        if not len(deltaF) == numseq:
            raise ValueError("Wrong number of deltaF arguments given.")

    # get limits
    xlim = kwargs.pop("xlim", None)
    ylim = kwargs.pop("ylim", None)
    colorlim = kwargs.pop("colorlim", None)
    if xlim:
        start,end = xlim
    else:
        start = min(epoch)
        end   = max(e + l.length * dt\
                    for e,l,dt in zip(epoch, sequencelist, deltaT))
    if not ydata is None and not ylim:
        ylim = [ydata.min(), ydata.max()]

    # get axis scales
    logx = kwargs.pop("logx", False)
    logy = kwargs.pop("logy", False)
    logcolor = kwargs.pop("logcolor", False)

    # get legend loc
    loc = kwargs.pop("loc", 0)
    alpha = kwargs.pop("alpha", 0.8)

    # get colorbar options
    hidden_colorbar = kwargs.pop("hidden_colorbar", False)

    # get savefig option
    bbox_inches = kwargs.pop("bbox_inches", None)
    
    #
    # get labels
    #

    xlabel = kwargs.pop("xlabel", None)
    if xlabel:
        unit = 1
    if not xlabel:
        unit, timestr = plotutils.time_axis_unit(end-start)
        if not t0:
            t0 = start
        t0 = lal.LIGOTimeGPS(t0)
        if int(t0.gpsNanoSeconds) == 0:
            xlabel = datetime.datetime(*lal.GPSToUTC(int(t0))[:6])\
                         .strftime("%B %d %Y, %H:%M:%S %ZUTC")
            xlabel = "Time (%s) since %s (%s)" % (timestr, xlabel, int(t0))
        else:
            xlabel = datetime.datetime(*lal.GPSToUTC(t0.gpsSeconds)[:6])\
                          .strftime("%B %d %Y, %H:%M:%S %ZUTC")
            xlabel = "Time (%s) since %s (%s)"\
                     % (timestr,
                        xlabel.replace(" UTC",".%.3s UTC" % t0.gpsNanoSeconds),\
                        t0)
        t0 = float(t0)
    ylabel     = kwargs.pop("ylabel", "Frequency (Hz)")
    colorlabel = kwargs.pop("colorlabel", "Amplitude")
    title      = kwargs.pop("title", "")
    subtitle   = kwargs.pop("subtitle", "")

    #
    # restrict data to the correct limits for plotting
    #

    interpolate = logy and ydata is None

    for i,sequence in enumerate(sequencelist):
        if interpolate:
           # interpolate the data onto a log-scale
           sequence,ydata = loginterpolate(sequence, f0[i], deltaF[i])
        if logy and ylim:
           plotted = (ydata > ylim[0]) & (ydata <= ylim[1])
           newVectorLength = int(plotted.sum())
           newsequence = lal.CreateREAL8VectorSequence(sequence.length,\
                                                           newVectorLength)
           for j in range(sequence.length):
               newsequence.data[j,:] = sequence.data[j,:][plotted]
           del sequence
           sequencelist[i] = newsequence
    if len(sequencelist) and logy and ylim:
        ydata = ydata[plotted]

    #
    # format bins
    #

    xbins = []
    for i in range(numseq):
        xmin = epoch[i]
        xmax = epoch[i] + sequencelist[i].length * deltaT[i]
        xbins.append(rate.LinearBins(float(xmin-t0)/unit, float(xmax-t0)/unit,\
                                     2))

    ybins = []
    for i in range(numseq):
        if ydata is not None:
            ydata = numpy.asarray(ydata)
            ymin = ydata.min()
            ymax = ydata.max()
        else:
            ymin = f0[i]
            ymax = f0[i] + sequencelist[i].vectorLength * deltaF[i]
        if logy:
            if ymin == 0:
                ymin = deltaF[i]
            ybins.append(rate.LogarithmicBins(ymin, ymax, 2))
        else:
            ybins.append(rate.LinearBins(ymin, ymax, 2))

    #
    # plot
    # 

    kwargs.setdefault("interpolation", "kaiser")

    plot = plotutils.ImagePlot(xlabel=xlabel, ylabel=ylabel, title=title,\
                               subtitle=subtitle, colorlabel=colorlabel)

    for sequence,x,y in zip(sequencelist, xbins, ybins):
        data = numpy.ma.masked_where(numpy.isnan(sequence.data), sequence.data,\
                                     copy=False)
        plot.add_content(data.T, x, y, **kwargs)

    # finalize
    plot.finalize(colorbar=True, logcolor=logcolor, minorticks=True,\
                  clim=colorlim)
    if hidden_colorbar:
        plotutils.add_colorbar(plot.ax, visible=False)

    # set logscale
    if logx:
        plot.ax.xaxis.set_scale("log")
    if logy:
        plot.ax.yaxis.set_scale("log")
    plot.ax._update_transScale()

    # format axes
    if xlim:
        xlim = (numpy.asarray(xlim).astype(float)-t0)/unit
        plot.ax.set_xlim(xlim)
    if ylim:
        plot.ax.set_ylim(ylim)

    # set grid and ticks
    plot.ax.grid(True, which="both")
    plotutils.set_time_ticks(plot.ax)
    plotutils.set_minor_ticks(plot.ax)

    # save and close
    plot.savefig(outfile, bbox_inches=bbox_inches,\
                 bbox_extra_artists=plot.ax.texts)
    plot.close()

def loginterpolate(sequence, y0, deltaY, N=None):
    """
    Interpolate a REAL?VectorSequence into logarithmic spacing in the second
    dimension (y-axis, or frequency).
    """
    # work out y-axis limits and new y-range
    if not N:
        N = sequence.vectorLength
    ylin = numpy.arange(sequence.vectorLength)*deltaY + y0
    ymin = y0 and y0 or deltaY
    ymax = y0 + sequence.vectorLength * deltaY
    ylog = numpy.logspace(numpy.log10(ymin*1.01), numpy.log10(ymax*0.99), num=N,\
                           endpoint=False)

    # make new sequence
    datatype = seriesutils.typecode(type(sequence))
    TYPESTR  = seriesutils._typestr[datatype]
    func     = getattr(lal, "Create%sVectorSequence" % TYPESTR)
    out      = func(sequence.length, N)

    # interpolate columnwise
    for i in range(sequence.length):
        intplt   = interpolate.interp1d(ylin, sequence.data[i,:])
        out.data[i,:] = intplt(ylog)

    return out,ylog

# =============================================================================
# Plot histogram of series
# =============================================================================

def plothistogram(serieslist, outfile, nonzero=False, num_bins=100,\
                  cumulative=False, bar=False, **kwargs):
    """
    Plots a line histogram of the data contained in the series, or list of
    series'.
    """
    # construct list of series
    if not hasattr(serieslist, "__contains__"):
        serieslist = [serieslist]

    # get limits
    xlim = kwargs.pop("xlim", None)
    ylim = kwargs.pop("ylim", None)

    # set labels
    xlabel   = kwargs.pop("xlabel", "Amplitude")
    ylabel   = kwargs.pop("ylabel", "Fraction of data")
    title    = kwargs.pop("title", "")
    subtitle = kwargs.pop("subtitle","")

    # get axis scales
    logx     = kwargs.pop("logx", False)
    logy     = kwargs.pop("logy", False)

    # get colorbar options
    hidden_colorbar = kwargs.pop("hidden_colorbar", False)

    # get savefig option
    bbox_inches = kwargs.pop("bbox_inches", "tight")

    # get fill
    fill = kwargs.pop("fill", False)

    # get legend loc
    loc = kwargs.pop("loc", 0)
    alpha = kwargs.pop("alpha", 0.8)

    # get colors to use
    color = kwargs.pop("color", None)
    if isinstance(color, str):
        color = color.split(',')
        if len(color) == 1:
            color = [color[0]]*len(serieslist)


    #
    # plot
    #

    plot = plotutils.LineHistogram(xlabel, ylabel, title, subtitle)

    for i, series in enumerate(serieslist):
        data = series.data.data
        if nonzero:
            data = data[data!=0]
        if color:
            kwargs["color"] = color[i]
        plot.add_content(data, normalize=data.size, label=series.name, **kwargs)
    plot.finalize(loc=loc, alpha=alpha, logx=logx, logy=logy, bar=bar,\
                  xlim=xlim, num_bins=num_bins, cumulative=cumulative)
    if hidden_colorbar:
        plotutils.add_colorbar(plot.ax, visible=False)

    plot.ax.autoscale_view(tight=True, scalex=True, scaley=True)

    if xlim:
        plot.ax.set_xlim(xlim)
    if ylim:
        plot.ax.set_ylim(ylim)

    plotutils.set_minor_ticks(plot.ax)

    # save and close
    plot.savefig(outfile, bbox_inches=bbox_inches,\
                 bbox_extra_artists=plot.ax.texts)
    plot.close()

