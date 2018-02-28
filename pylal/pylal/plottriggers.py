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
This modules provides some wrapper functions for the classes defined in pylal.plotutils.
"""

# =============================================================================
# Preamble
# =============================================================================

from __future__ import division
import re
import copy
import numpy
import datetime
from pylal import git_version
from pylal import date
from pylal import plotutils
from pylal.datatypes import LIGOTimeGPS
from glue.ligolw import lsctables
from glue.ligolw import table

pylab      = plotutils.pylab
matplotlib = pylab.matplotlib

__author__  = "Duncan M. Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__    = git_version.date


jet = copy.deepcopy(matplotlib.cm.jet)
# remap red scale
jet._segmentdata['red'] = ([x == jet._segmentdata["red"][-1] and (1,0.8,0.8)\
                            or x for x in jet._segmentdata['red'] ])
jet._segmentdata['blue'] = ([x == jet._segmentdata["blue"][0] and (0,0.7,0.7)\
                            or x for x in jet._segmentdata['blue'] ])

# =============================================================================
# Get column from a table
# =============================================================================

def get_column(lsctable, column):
    """
    Extract column from the given glue.ligolw.table lsctable as numpy array.
    Tries to use a 'get_col() function if available, otherwise uses
    getColumnByName(), treating 'time' as a special case for the known
    Burst/Inspiral/Ringdown tables.
    """
    # format column
    column = str(column).lower()
    obj_type = str(type(lsctable))

    # if there's a 'get_' function, use it
    if hasattr(lsctable, 'get_%s' % column):
        return numpy.asarray(getattr(lsctable, 'get_%s' % column)())

    # treat 'time' as a special case
    if column == 'time':
        # find time column
        if re.match("sim_inspiral", lsctable.tableName, re.I):
            timecolumn = "geocent_end_time"
        elif re.match("sim_burst", lsctable.tableName, re.I):
            timecolumn = "time_geocent_gps"
        elif re.match("sim_ringdown", lsctable.tableName, re.I):
            timecolumn = "geocent_start_time"
        elif re.search("inspiral", lsctable.tableName, re.I):
            timecolumn = "end_time"
        elif re.search("burst", lsctable.tableName, re.I):
            timecolumn = "peak_time"
        elif re.search("(ringdown|stochastic)", lsctable.tableName, re.I):
            timecolumn = "start_time"

        # return the correctly combined columns
        return numpy.asarray(lsctable.getColumnByName(timecolumn)) + \
               numpy.asarray(lsctable.getColumnByName("%s_ns"%timecolumn))*1e-9

    return numpy.asarray(lsctable.getColumnByName(column))

# =============================================================================
# Plot LIGOLW trigger table
# =============================================================================

def plottable(lsctable, outfile, xcolumn="time", ycolumn="snr",\
              colorcolumn=None, rankcolumn=None, t0=None, **kwargs):
    """
    Plot any column of a valid LigoLW table against any other, coloured by any
    other, or plot all the same if you're that way inclined. Multiple tables
    can be given for any non-coloured plot.

    "time" as a column means the output of get_peak for Burst tables, get_end
    for Inspiral tables, and get_start for ringdown tables

    Arguments:

        tables : [ dict | glue.ligolw.table.Table ]
            dict of ("name", table) pairs, or single LigoLW tables
        outfile : string
            string path for output plot

    Keyword arguments:

        xcolumn : string
            valid column of triggers table to plot on x-axis
        ycolumn : string
            valid column of triggers table to plot on y-axis
        zcolumn : string
            valid column of triggers table to use for colorbar (optional).
    """
    # work out dictionary or table
    if isinstance(lsctable, table.Table):
        tables = {"_":lsctable}
    else:
        tables = lsctable
    tablenames = tables.keys()
    tables     = [tables[n] for n in tablenames]

    # get axis limits
    xlim = kwargs.pop('xlim', None)
    ylim = kwargs.pop('ylim', None)
    zlim = kwargs.pop('zlim', None)
    colorlim = kwargs.pop('colorlim', None)
    if zlim and not colorlim:
        colorlim = zlim

    # set up columns
    columns = list(map(str.lower, [xcolumn, ycolumn]))
    if colorcolumn:
        columns.append(colorcolumn.lower())
        if rankcolumn:
            columns.append(rankcolumn.lower())
        else:
            columns.append(colorcolumn.lower())

    # set up limits
    limits = [xlim, ylim, zlim, None]

    # get zero
    if "time" in columns and not t0:
        if xlim:
            t0 = float(xlim[0])
        else:
            t0 = numpy.infty
            for t in tables:
                timedata = get_column(t, "time").astype(float)
                if len(timedata):
                    t0 = min(t0, timedata.min())
            if numpy.isinf(t0):
                t0 = 0
            t0 = int(t0)

    #
    # get data
    #

    # extract columns
    data = list()
    for i,name in enumerate(tablenames):
        data.append(list())
        for c in columns:
            c = c.lower()
            if ((c not in tables[i].columnnames and c != "time") and not
                    hasattr(tables[i], "get_%s" % c.lower())):
                raise RuntimeError("Column '%s' not found in table '%s'."\
                                   % (c,name))
            data[-1].append(get_column(tables[i], c.lower()).astype(float))
            if not len(data[-1][-1].shape)\
            or data[-1][-1].shape[0] != len(tables[i]):
                raise AttributeError("No data loaded for column \"%s\" and "\
                                     "table \"%s\"" % (c, name))

    # add our own limits
    for i in range(len(columns)):
        if not limits[i]:
            mins = [data[j][i].min() for j in range(len(tablenames))\
                    if len(data[j][i].shape) and data[j][i].shape[0] != 0]
            if len(mins):
                lmin = min(mins)
                lmax = max(data[j][i].max() for j in range(len(tablenames))\
                           if len(data[j][i]))
                limits[i] = [lmin, lmax]

    # get time unit
    if "time" in columns:
        idx = columns.index("time")
        if limits[idx]:
            unit, timestr = plotutils.time_axis_unit(limits[idx][1]\
                                                     - limits[idx][0])
        else:
            unit, timestr = plotutils.time_axis_unit(1)

    # format data
    for i,name in enumerate(tablenames):
        for j,column in enumerate(columns):
            if column == "time":
                data[i][j] = (data[i][j] - t0) / unit
                if i==0 and limits[j]:
                    limits[j] = [float(limits[j][0] - t0) / unit,\
                                 float(limits[j][1] - t0) / unit]

        # format a condition to reject triggers outside the plot range
        plotted = True
        for j,column in enumerate(columns):
            if limits[j]:
                plotted = plotted & (limits[j][0] <= data[i][j])\
                                  & (limits[j][1] >=  data[i][j])

        # apply the limits
        if not isinstance(plotted, bool):
            for j,d in enumerate(data[i]):
                data[i][j] = d[plotted]

    #
    # find loudest event
    #

    loudest = None
    if len(columns) == 4 and len(tablenames) == 1 and len(data[0][0]) > 0:
        idx = data[0][3].argmax()   
        loudest = [data[0][j][idx] for j in range(len(columns))]

    #
    # get or set the labels
    #

    label = {}
    for i,(label,column) in enumerate(zip(["xlabel", "ylabel", "colorlabel"],\
                                          columns)):
        # get given label
        l = kwargs.pop(label, None)
        if l is None:
            # format time string
            if column == "time" and limits[i]:
                zerostr = datetime.datetime(\
                              *date.XLALGPSToUTC(LIGOTimeGPS(t0))[:6])\
                                       .strftime("%B %d %Y, %H:%M:%S %ZUTC")
                l = "Time (%s) since %s (%s)" % (timestr, zerostr, t0)
            # format any other column
            else:
                l = plotutils.display_name(column)
        if label == "xlabel":
            xlabel = l
        elif label == "ylabel":
            ylabel = l
        else:
            colorlabel = l

    title = kwargs.pop("title", "")
    subtitle = kwargs.pop("subtitle", None)
    if subtitle is None and loudest:
        subtitle = "Loudest event by %s:" % plotutils.display_name(columns[-1])
        for j,c in enumerate(columns):
            if j!= 0 and c == columns[j-1]: continue
            lstr = loudest[j]
            if c == "time":
                lstr = lstr * unit + t0
            subtitle += " %s=%.2f" % (plotutils.display_name(c), lstr)
    else:
        loudest = None

    #
    # get parameters
    #

    # get axis scales
    logx     = kwargs.pop("logx", False)
    logy     = kwargs.pop("logy", False)
    logcolor = kwargs.pop("logcolor", False)

    # get legend loc
    loc  = kwargs.pop("loc", 0)

    # get colorbar options
    hidden_colorbar = kwargs.pop("hidden_colorbar", False)

    # get savefig option
    bbox_inches = kwargs.pop("bbox_inches", "tight")

    # get detchar plot params
    dqstyle  = (kwargs.pop("detchar", False) or
                kwargs.pop('detchar_style', False))
    dqthresh = (kwargs.pop('dcthreshold', False) or
                kwargs.pop('detchar_style_threshold', 10))

    # get greyscale param
    greyscale = kwargs.pop("greyscale", False)
    if greyscale and not kwargs.has_key("cmap"):
        kwargs["cmap"] =\
            pylab.matplotlib.colors.LinearSegmentedColormap("clrs",\
                                        pylab.matplotlib.cm.hot_r._segmentdata)
    elif not kwargs.has_key("cmap"):
        kwargs["cmap"] = jet

    #
    # make the plot
    #

    tablenames = map(plotutils.display_name, tablenames)

    if len(columns) == 2:
        plot = plotutils.ScatterPlot(xlabel, ylabel,title, subtitle)
        for i in range(len(tablenames)):
            plot.add_content(data[i][0], data[i][1], label=tablenames[i],\
                             **kwargs)
        plot.finalize(loc=loc)
        if hidden_colorbar:
            plotutils.add_colorbar(plot.ax, visible=False)
    else:
        if dqstyle:
            plot = plotutils.DQScatterPlot(xlabel, ylabel, colorlabel,\
                                           title, subtitle)
        else:
            plot = plotutils.ColorbarScatterPlot(xlabel, ylabel, colorlabel,\
                                                 title, subtitle)
        plot.add_content(data[0][0], data[0][1], data[0][2],\
                         label=tablenames[0], **kwargs)
        kwargs.pop("cmap", None)
        if dqstyle:
            plot.finalize(logcolor=logcolor, clim=colorlim, loc=loc,\
                          threshold=dqthresh)
        else:
            plot.finalize(logcolor=logcolor, clim=colorlim, loc=loc)

    # plot loudest
    if loudest:
        plot.ax.plot([loudest[0]], [loudest[1]], marker="*", color="gold",\
                     markersize=15)

    # set axes
    if logx:
        plot.ax.set_xscale("log")
    if logy:
        plot.ax.set_yscale("log")
    plot.ax.autoscale_view(tight=True, scalex=True, scaley=True)

    if limits[0]:
        plot.ax.set_xlim(limits[0])
    if limits[1]:
        plot.ax.set_ylim(limits[1])

    plot.ax.grid(True, which="both")
    if xcolumn == "time":
        plotutils.set_time_ticks(plot.ax)
    plotutils.set_minor_ticks(plot.ax)

    # save and close
    if greyscale:
        plot.ax.patch.set_facecolor("#E8E8E8")
    plot.savefig(outfile, bbox_inches=bbox_inches,\
                 bbox_extra_artists=plot.ax.texts)
    plot.close()
    
# =============================================================================
# Plot a histogram of any column
# =============================================================================

def plothistogram(lsctable, outfile, column="snr", numbins=100,\
                  colorcolumn=None, colorbins=10, normalize=None,\
                  cumulative=None, bar=False, **kwargs):
    """
    Plot any column of a valid LigoLW table against any other, coloured by any
    other, or plot all the same if you're that way inclined. Multiple tables
    can be given for any non-coloured plot.

    "time" as a column means the output of get_peak for Burst tables, get_end
    for Inspiral tables, and get_start for ringdown tables

    Arguments:

        table : [ dict | glue.ligolw.table.Table ]
            dict of ("name", table) pairs, or single LigoLW table
        outfile : string
            string path for output plot

    Keyword arguments:
        xcolumn : string
            valid column of triggers table to plot on x-axis
        ycolumn : string
            valid column of triggers table to plot on y-axis
        zcolumn : string
            valid column of triggers table to use for colorbar (optional).
    """
    # work out dictionary or table
    if isinstance(lsctable, table.Table):
        tables = {"_":lsctable}
    else:
        tables = lsctable
    tablenames = sorted(tables.keys())
    tables     = [tables[n] for n in tablenames]

    # get axis limits
    xlim = kwargs.pop("xlim", None)
    ylim = kwargs.pop("ylim", None)
    zlim = kwargs.pop("zlim", None)
    clim = kwargs.pop("clim", zlim)

    # get axis scales
    logx     = kwargs.pop("logx", False)
    logy     = kwargs.pop("logy", False)
    logcolor = kwargs.pop("logcolor", False)

    # get colorbar options
    hidden_colorbar = kwargs.pop("hidden_colorbar", False)

    # get savefig option
    bbox_inches = kwargs.pop("bbox_inches", "tight")

    # get fill
    fill = kwargs.pop("fill", False)

    # get extras
    rate = kwargs.pop("rate", False)
    if normalize:
        rate = True
        if not hasattr(normalize, "__iter__"):
           normalize = [normalize]*len(tablenames)
    else:
        normalize = [1]*len(tablenames)
    loc        = kwargs.pop("loc", 0)

    # set labels
    xlabel = kwargs.pop("xlabel", plotutils.display_name(column))
    if rate and cumulative:
        ylabel = kwargs.pop("ylabel", "Cumulative rate (Hz)")
    elif rate:
        ylabel = kwargs.pop("ylabel", "Rate (Hz)")
    elif not rate and cumulative:
        ylabel = kwargs.pop("ylabel", "Cumulative number")
    else:
        ylabel = kwargs.pop("ylabel", "Number")

    title = kwargs.pop("title", "")
    subtitle = kwargs.pop("subtitle","")
    colorlabel  = kwargs.pop("colorlabel", colorcolumn\
                                  and plotutils.display_name(colorcolumn) or "")

    #
    # get column data
    #

    if colorcolumn:
        colordata = get_column(tables[0], colorcolumn).astype(float)
        if not clim and colordata.size != 0:
            clim = [colordata.min(), colordata.max()]
        elif not clim:
            clim = [1, 10]
        if isinstance(colorbins, int):
            numcolorbins = colorbins
            if logcolor:
                colorbins = numpy.logspace(*numpy.log10(clim), num=numcolorbins)
            else:
                colorbins = numpy.linspace(*clim, num=numcolorbins)
        else:
            numcolorbins = len(colorbins)

    data = list()
    for i,lsctable in enumerate(tables):
        d = get_column(lsctable, column).astype(float)
        if i==0 and colorcolumn:
            data.append(list())
            for j in range(numcolorbins):
                data[i].append(d[d>=colorbins[j]])
        else:
            data.append(d)

    #
    # generate plot and make
    #

    if colorcolumn:
        raise NotImplementedError("This has not been implemented in "+\
                                  "plotutils yet, here's your chance!")
        plot = plotutils.ColorbarLineHistogram(xlabel, ylabel, colorlabel,\
                                               title, subtitle)
        for dataset,colorbin in zip(data[0], colorbins):
            plot.add_content(dataset, normalize=normalize[0],\
                             colorvalue=colorbin, label=tablenames[0], **kwargs)
        for i,name in enumerate(tablenames[1:]):
            plot.add_content(data[i+1], normalize=normalize[i+1],\
                             label=name, **kwargs)
        plot.finalize(logcolor=logcolor, colorlabel=colorlabel, loc=loc)
    else:
        color = kwargs.pop("color", None)
        if isinstance(color, str):
            color = color.split(',')
            if len(color) == 1:
                color = [color[0]]*len(serieslist)
        plot = plotutils.LineHistogram(xlabel, ylabel, title, subtitle)
        for i,name in enumerate(tablenames):
            if color:
                kwargs["color"] = color[i]
            plot.add_content(data[i], normalize=normalize[i], label=name,\
                             **kwargs)
        plot.finalize(loc=loc, logx=logx, logy=logy, bar=bar, num_bins=numbins,\
                      cumulative=cumulative)

    plot.ax.autoscale_view(tight=True, scalex=True, scaley=True)

    if xlim:
        plot.ax.set_xlim(xlim)
    if ylim:
        plot.ax.set_ylim(ylim)

    # save and close
    plot.savefig(outfile, bbox_inches=bbox_inches,\
                 bbox_extra_artists=plot.ax.texts)
    plot.close()

# =============================================================================
# Plot trigger rate
# =============================================================================

def plotrate(lsctable, outfile, stride=60, column="peak_frequency", bins=[],\
             t0=None, **kwargs):
    """
    Plot rate versus time for the given ligolw table triggers, binned by the
    given bincolumn using the bins list.

    Arguments:

        lsctable : glue.ligolw.table.Table
            LIGOLW table containing a list of triggers
        outfile : string
            string path for output plot
    """
    # get column data
    timedata = get_column(lsctable, "time")
    ratedata = get_column(lsctable, column)
    # sort in time
    if timedata.size:
        timedata, ratedata = map(numpy.asarray,\
                                 zip(*sorted(zip(timedata,ratedata),\
                                             key=lambda (t,r): t)))
    

    # get x-axis limits
    xlim = kwargs.pop("xlim", None)
    if xlim:
        start,end = xlim
    elif timedata.size:
        start = timedata.min()
        end   = timedata.max()
        xlim  = start, end
    else:
        start = 0
        end   = 1

    # get other params
    logx            = kwargs.pop("logx", False)
    logy            = kwargs.pop("logy", False)
    ylim            = kwargs.pop("ylim", False)
    hidden_colorbar = kwargs.pop("hidden_colorbar", False)
    bbox_inches     = kwargs.pop("bbox_inches", "tight")
    loc             = kwargs.pop("loc", 0)

    # get bins
    xbins = numpy.arange(float(start), float(end), stride)
    if not bins:
        bins = [[0, float("inf")]]
    ybins = [map(float, bin_) for bin_ in bins]

    # remove all triggers not in any bins
    xextent = (xbins[0], xbins[-1]+stride)
    yextent = (ybins[0][0], ybins[-1][-1])
    inanybin = (timedata > xextent[0]) & (timedata <= xextent[1])\
               & (ratedata > yextent[0]) & (ratedata <= yextent[1])
    timedata = timedata[inanybin]
    ratedata = ratedata[inanybin]

    ydata = numpy.zeros((len(xbins), len(ybins)))
    for i,xbin in enumerate(xbins):
        # work out index of first trigger outside x bin
        try:
            idx = (timedata >= xbin+stride).nonzero()[0][0]
            # get triggers in this time bin
            ratebin = numpy.asarray(sorted(ratedata[:idx]))
            break_ = False
            # get triggers in this time bin
            ratebin = numpy.asarray(sorted(ratedata[:idx]))
        except IndexError:
            idx = None
            ratebin = numpy.asarray(sorted(ratedata))
            break_ = True
        # loop over ybins
        for j,ybin in enumerate(ybins):
            # work out index of first trigger outside this y bin
            try:
                yidx = (ratebin >= ybin[1]).nonzero()[0][0]
            except IndexError:
                ydata[i][j] = len(ratebin)/stride
                break
            ydata[i][j] = len(ratebin[:yidx])/stride
            ratebin = ratebin[yidx:]
        if break_: break
        timedata = timedata[idx:]
    xdata = xbins + stride/2
 
    # set xlabel and renomalize time
    if t0:
        xdata -= t0
        start -= t0
        end -= t0
    xlabel = kwargs.pop("xlabel", None)
    if not xlabel:
        if not t0:
            t0 = start
            start -= t0
            end   -= t0
            xdata -= t0
        unit, timestr = plotutils.time_axis_unit(end-start)
        start /= unit
        end   /= unit
        xdata = xdata/unit
        t0 = LIGOTimeGPS(t0)
        if t0.nanoseconds==0:
            xlabel = datetime.datetime(*date.XLALGPSToUTC(LIGOTimeGPS(t0))[:6])\
                         .strftime("%B %d %Y, %H:%M:%S %ZUTC")
        else:
            xlabel = datetime.datetime(*date.XLALGPSToUTC(\
                                                  LIGOTimeGPS(t0.seconds))[:6])\
                          .strftime("%B %d %Y, %H:%M:%S %ZUTC")
            xlabel = "Time (%s) since %s (%s)"\
                     % (timestr,
                        xlabel.replace(" UTC", ".%.3s UTC" % t0.nanoseconds),\
                        t0)
    ylabel   = kwargs.pop("ylabel",   "Rate (Hz)")
    title    = kwargs.pop("title",    "")
    subtitle = kwargs.pop("subtitle", "")

    #
    # plot object
    #

    plot = plotutils.ScatterPlot(xlabel, ylabel, title, subtitle)

    # mask zeros for log plot
    if logy:
        ydata = numpy.ma.masked_where(ydata==0, ydata, copy=False)

    for i,ybin in enumerate(ybins):
        if list(ybin) == [0, numpy.inf]:
            label = "_"
        else:
            label = "%s-%s" % (ybin[0], ybin[1])
        plot.add_content(xdata, ydata[:,i], label=label, **kwargs)

    plot.finalize(loc=loc)
    if hidden_colorbar:
        plotutils.add_colorbar(plot.ax, visible=False)
    if logx:
        plot.ax.set_xscale("log")
    if logy:
        plot.ax.set_yscale("log")
    plot.ax.autoscale_view(tight=True, scalex=True, scaley=True)
    plot.ax.set_xlim(start, end)
    if ylim:
        plot.ax.set_ylim(ylim)

    # set grid and ticks
    plot.ax.grid(True, which="both")
    plotutils.set_time_ticks(plot.ax)
    plotutils.set_minor_ticks(plot.ax)

    # save
    plot.savefig(outfile, bbox_inches=bbox_inches,\
                 bbox_extra_artists=plot.ax.texts)
    plot.close()
