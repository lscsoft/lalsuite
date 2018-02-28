# encoding: utf-8
#
# Copyright (C) 2008  Nickolas V Fotopoulos
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

__author__ = "Nickolas Fotopoulos <nvf@gravity.phys.uwm.edu>"

import datetime, re
import numpy

from warnings import warn
from glue import segments
from pylal import plotutils
from pylal import date
from pylal.datatypes import LIGOTimeGPS
import pylab

##############################################################################
# Plotting class

class PlotSegmentsPlot(object):
  """
  Represent a plotsegments plot.  To use it, one must instantiate a
  PlotSegmentsPlot object, call add_contents(), then call finalize().
  To save, call the savefig() method.  To display to the screen, call
  pylab.show().
  
  Developer note: There is the transformation function _time_transform.
  Follow the convention of applying it before storing internal data so that
  there is no ambiguity about whether data has been transformed or not.
  """
  color_code = {'H1':'r', 'H2':'b', 'L1':'g', 'V1':'m', 'G1':'k'}
  
  def __init__(self, t0=0):
    """
    Create a fresh plot.  Provide t0 to provide a reference time to use as
    zero.
    """
    warn("This class is now deprecated by pylal.plotutils.SegmentPlot",\
         DeprecationWarning)
    self.fig = pylab.figure()
    self.ax = self.fig.add_subplot(111)
    self.savefig = self.fig.savefig
    
    self.window = None
    self.ifos = []
    self._time_transform = lambda t: t - t0
    
    self.ax.set_xlabel("time (s)")
    self.ax.set_ylabel("IFO")
  
  def add_contents(self, segdict, ifos=None):
    """
    Add the contents of segdict to the plot.  Provide the list ifos, if you
    wish to specify the order in which they appear or wish to plot a subset of
    the segmentlists in segdict.
    """
    if ifos is None:
      ifos = segdict.keys()
      ifos.sort()
    self.ifos.extend(ifos[::-1])
    for row, ifo in enumerate(ifos[::-1]):
      color = self.color_code[ifo]
      for seg in segdict[ifo]:
        a = self._time_transform(seg[0])
        b = self._time_transform(seg[1])
        self.ax.fill([a, b, b, a, a], [row, row, row+1, row+1, row], color)
  
  def set_window(self, window_seg, padding=0):
    """
    Define a window of interest by setting the x-limits of the plot
    appropriately.  If padding is also present, protract the x-limits by
    that quantity and mark the unpadded window with solid black lines.
    """
    a = self._time_transform(window_seg[0])
    b = self._time_transform(window_seg[1])
    self.window = segments.segment((a - padding, b + padding))
    
    if padding > 0:
      self.ax.axvline(a, color='k', linewidth=2)
      self.ax.axvline(b, color='k', linewidth=2)
  
  def highlight_segment(self, seg):
    """
    Highlight a particular segment with dashed lines.
    """
    self.ax.axvline(self._time_transform(seg[0]), color='k', linestyle='--')
    self.ax.axvline(self._time_transform(seg[1]), color='k', linestyle='--')
  
  def finalize(self):
    """
    Make final changes to the plot with the guarantee that no additional data
    will be added.
    """
    ticks = pylab.arange(len(self.ifos)) + 0.5
    self.ax.set_yticks(ticks)
    self.ax.set_yticklabels(self.ifos)
    
    if self.window is not None:
      self.ax.set_xlim(self.window)
    self.ax.set_ylim((0, len(self.ifos)))

  def close(self):
    pylab.close(self.fig)
  
  def __del__(self):
    self.close()

# =============================================================================
# Plot segments
# =============================================================================

def plotsegmentlistdict(segdict, outfile, keys=None, t0=None,\
                        highlight_segments=None, insetlabels=None,\
                        **kwargs):
    """
    Plots a glue.segments.segmentlistdict using the PlotSegmentsPlot class from
    pylal.plotutils.

    Arguments:

        segdict : glue.segments.segmentlistdict
            dictionary of (name, segmentlist) pairs to plot.
        outfile : str
            filepath to write to

    Keyword arguments:

        keys : list
            ordered list of keys to use in this order on the plot
        t0 : float
            GPS time around which to zero plot
        highlight_segments : glue.segments.segmentlistdict
            list of segments to highlight with vertical red dashed lines
        insetlabels : [ True | False ]
            write labels inside the plot axis

    Unnamed keyword arguments:

        xlabel : string
            label for x-axis
        ylabel : string
            label for y-axis
        title : string
            title for plot
        subtitle : string
            subtitle for plot
        bbox_inches : str
            use "tight" to get a bounding box tight to the axis.
    """
    # get time limits
    xlim = kwargs.pop("xlim", None)
    if xlim is None:
        try:
            extents = [seg.extent() for seg in segdict.values()]
            start = min(s[0] for s in extents)
            end   = max(s[1] for s in extents)
        except ValueError:
            start = 0
            end   = 1 
        xlim  = start,end
    else:
        start,end = xlim

    # get unit for plot
    unit, timestr = plotutils.time_axis_unit(end-start)

    # set xlabel and renomalize time
    xlabel = kwargs.pop("xlabel", None)
    if not xlabel:
        if not t0:
            t0 = start
        unit, timestr = plotutils.time_axis_unit(end-start)
        t0 = LIGOTimeGPS(float(t0))
        if t0.nanoseconds==0:
            xlabel = datetime.datetime(*date.XLALGPSToUTC(t0)[:6])\
                         .strftime("%B %d %Y, %H:%M:%S %ZUTC")
            xlabel = "Time (%s) since %s (%s)" % (timestr, xlabel, int(t0))
        else:
            xlabel = datetime.datetime(*date.XLALGPSToUTC(\
                                                  LIGOTimeGPS(t0.seconds))[:6])\
                          .strftime("%B %d %Y, %H:%M:%S %ZUTC")
            xlabel = "Time (%s) since %s (%s)"\
                     % (timestr,
                        xlabel.replace(" UTC", ".%.3s UTC" % t0.nanoseconds),\
                        t0)
        t0 = float(t0)
    ylabel   = kwargs.pop("ylabel",   "")
    title    = kwargs.pop("title",    "")
    subtitle = kwargs.pop("subtitle", "")

    # get other parameters
    labels_inset    = kwargs.pop("labels_inset", False)
    bbox_inches     = kwargs.pop("bbox_inches", "tight")
    hidden_colorbar = kwargs.pop("hidden_colorbar", False)

    # escape underscores for latex text
    if not keys:
        keys = segdict.keys()
    if pylab.rcParams["text.usetex"]:
        newdict = segments.segmentlistdict()
        for i,key in enumerate(keys):
           newkey = re.sub('(?<!\\\\)_', '\_', key)
           keys[i] = newkey
           newdict[newkey] = segdict[key]
        segdict = newdict

    # generate plot
    plot = plotutils.PlotSegmentsPlot(xlabel, ylabel, title, subtitle,\
                                       t0=t0, dt=unit)
    plot.add_content(segdict, keys, **kwargs)
    plot.finalize(labels_inset=insetlabels)

    # highlight segments
    if highlight_segments:
        for seg in highlight_segments:
            plot.highlight_segment(seg)
 
    # add colorbar
    if hidden_colorbar:
        plotutils.add_colorbar(plot.ax, visible=False)

    # set axis limits
    xlim = [float(start-t0)/unit, float(end-t0)/unit]
    plot.ax.set_xlim(*xlim)

    # set grid
    plot.ax.grid(True, which="both")
    plotutils.set_time_ticks(plot.ax)
    plotutils.set_minor_ticks(plot.ax, x=False)

    # save
    plot.savefig(outfile, bbox_inches=bbox_inches,\
                 bbox_extra_artists=plot.ax.texts)
    plot.close()

# =============================================================================
# Plot segment histogram
# =============================================================================

def plothistogram(segdict, outfile, keys=None, num_bins=100, **kwargs):
    """
    Plots a histogram of segment duration for each entry in the
    glue.segments.segmentlistdict segdict.

    Arguments:

      segdict : glue.segments.segmentlistdict
          list of segments with which to veto triggers, use dict for multiple
          datasets
      outfile : string 
          string path for output plot

    Keyword arguments:

        keys : [ str | list]
            ordered list of keys to use in this order on the plot
        num_bins : int
            number of bins.

    Unnamed keyword arguments:

        logx : [ True | False ]
            display x-axis in log scale
        logy : [ True | False ]
            display y-axis in log scale
        xlim : tuple
            (xmin,xmax) limits for x-axis 
        ylim : tuple
            (ymin,ymax) limits for y-axis 
        xlabel : string
            label for x-axis
        ylabel : string
            label for y-axis
        title : string
            title for plot
        subtitle : string
            subtitle for plot
        bbox_inches : str
            use "tight" to get a bounding box tight to the axis.
    """
    # get limits
    xlim = kwargs.pop('xlim', None)
    ylim = kwargs.pop('ylim', None)

    # get labels
    xlabel = kwargs.pop('xlabel', 'Length of segment (seconds)')
    ylabel = kwargs.pop('ylabel', 'Number of segments')
    title  = kwargs.pop('title',  'Segment Duration Histogram')
    subtitle = kwargs.pop('subtitle', "")

    # get axis scale
    logx = kwargs.pop('logx', False)
    logy = kwargs.pop('logy', False)

    # get savefig option
    bbox_inches = kwargs.pop('bbox_inches', 'tight')
    hidden_colorbar = kwargs.pop("hidden_colorbar", False)
 
    # escape underscores for latex text
    if not keys:
        keys = segdict.keys()
    if pylab.rcParams["text.usetex"]:
        newdict = segments.segmentlistdict()
        for i,key in enumerate(keys):
           newkey = re.sub('(?<!\\\\)_', '\_', key)
           keys[i] = newkey
           newdict[newkey] = segdict[key]
        segdict = newdict

    # generate plot object
    plot = plotutils.VerticalBarHistogram(xlabel, ylabel, title, subtitle)

    # add each segmentlist
    for flag in keys:
        plot.add_content([float(abs(seg)) for seg in segdict[flag]],\
                          label=flag, **kwargs)

    # finalize plot with histograms
    plot.finalize(num_bins=num_bins, logx=logx, logy=logy)

    # set limits
    plot.ax.autoscale_view(tight=True, scalex=True, scaley=True)
    if ylim:
      plot.ax.set_ylim(map(float, ylim))
    if xlim:
      plot.ax.set_xlim(map(float, xlim))

    # save figure
    plot.savefig(outfile, bbox_inches=bbox_inches,\
                 bbox_extra_artists=plot.ax.texts)
    plot.close()

# =============================================================================
# Plot duty cycle
# =============================================================================

def plotdutycycle(segdict, outfile, binlength=3600, keys=None, t0=None,\
                  showmean=False, **kwargs):
    """
    Plot the percentage duty cycle each flag in the given
    glue.segments.segmentlistdict, binned over the given duration.
    """
    # get time limits
    xlim = kwargs.pop("xlim", None)
    if xlim is None:
        try:
            extents = [seg.extent() for seg in segdict.values()]
            start = min(s[0] for s in extents)
            end   = max(s[1] for s in extents)
        except ValueError:
            start = 0
            end   = 1 
        xlim  = start,end
    else:
        start,end = xlim

    # get unit for plot
    unit, timestr = plotutils.time_axis_unit(end-start)

    # set xlabel and renomalize time
    xlabel = kwargs.pop("xlabel", None)
    if not xlabel:
        if not t0:
            t0 = start
        unit, timestr = plotutils.time_axis_unit(end-start)
        t0 = LIGOTimeGPS(float(t0))
        if t0.nanoseconds==0:
            xlabel = datetime.datetime(*date.XLALGPSToUTC(t0)[:6])\
                         .strftime("%B %d %Y, %H:%M:%S %ZUTC")
            xlabel = "Time (%s) since %s (%s)" % (timestr, xlabel, int(t0))
        else:
            xlabel = datetime.datetime(*date.XLALGPSToUTC(\
                                                  LIGOTimeGPS(t0.seconds))[:6])\
                          .strftime("%B %d %Y, %H:%M:%S %ZUTC")
            xlabel = "Time (%s) since %s (%s)"\
                     % (timestr,
                        xlabel.replace(" UTC", ".%.3s UTC" % t0.nanoseconds),\
                        t0)
        t0 = float(t0)
        xlim[0] = (start - t0)/unit
        xlim[1] = (end - t0)/unit
    ylabel   = kwargs.pop("ylabel",   "")
    title    = kwargs.pop("title",    "")
    subtitle = kwargs.pop("subtitle", "")

    # get other parameters
    loc = kwargs.pop("loc", 0)
    legalpha = kwargs.pop("alpha", 0.8)
    labels_inset    = kwargs.pop("labels_inset", False)
    bbox_inches     = kwargs.pop("bbox_inches", "tight")
    hidden_colorbar = kwargs.pop("hidden_colorbar", False)

    # escape underscores for latex text
    if not keys:
        keys = segdict.keys()
    if pylab.rcParams["text.usetex"]:
        newdict = segments.segmentlistdict()
        for i,key in enumerate(keys):
           newkey = re.sub('(?<!\\\\)_', '\_', key)
           keys[i] = newkey
           newdict[newkey] = segdict[key]
        segdict = newdict

    #
    # generate duty cycle info
    #

    # generate bins
    binlength = float(binlength)
    if int(end-start) % binlength == 0:
        numbins = int(end-start)/binlength
    else:
        numbins = float(end-start)//binlength+1
    bins = numpy.arange(float(start), float(end), binlength) + binlength/2
    duty = dict((key, numpy.zeros(numbins)) for key in keys)

    bs = float(start)
    for i in range(numbins):
        be = float(bs + binlength)
        seg = segments.segmentlist([segments.segment(bs, be)])
        for key in keys:
            duty[key][i] = float(abs(segdict[key] & seg))/abs(seg) * 100
        bs += binlength

    if showmean:
      mean = dict((key, numpy.zeros(numbins)) for key in keys)
      for key in keys:
        mean[key] = [duty[key][:i+1].mean() for i in range(numbins)]

    #
    # generate plot
    #

    bins = (bins-t0)/unit

    plot = plotutils.BarPlot(xlabel, ylabel, title, subtitle)
    for i,key in enumerate(keys):
      if showmean:
        thislabel = plotutils.display_name(key) + ' (%.2f\%%)' % (mean[key][-1])
      else:
        thislabel = plotutils.display_name(key)
      plot.add_content(bins, duty[key], label=thislabel,\
                         alpha=0.8, width=binlength/unit)
    plot.finalize(loc=loc, alpha=legalpha)

    # add running mean
    if showmean:
      for i,key in enumerate(keys):
        print i, key
        plot.ax.plot(bins, mean[key], linestyle = '--')
      plot.ax.get_legend().get_frame().set_alpha(0.5)

    # add colorbar
    if hidden_colorbar:
        plotutils.add_colorbar(plot.ax, visible=False)

    # set limits
    plot.ax.autoscale_view(tight=True, scalex=True)
    if xlim:
        plot.ax.set_xlim(map(float, xlim))
    plot.ax.set_ylim(0, 100)

    # set grid
    plot.ax.grid(True, which="both")
    plotutils.set_time_ticks(plot.ax)
    plotutils.set_minor_ticks(plot.ax, x=False)

    # save figure
    plot.savefig(outfile, bbox_inches=bbox_inches,\
                 bbox_extra_artists=plot.ax.texts)
