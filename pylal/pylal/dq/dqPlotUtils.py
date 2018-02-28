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

# =============================================================================
# Preamble
# =============================================================================

from __future__ import division
import math,re,numpy,itertools,copy,matplotlib,sys,warnings

# test matplotlib backend and reset if needed
from os import getenv
_display = getenv('DISPLAY','')
_backend_warn = """No display detected, moving to 'Agg' backend in matplotlib."""
if not _display and not matplotlib.get_backend().lower() == 'agg':
  warnings.warn(_backend_warn)
  matplotlib.use('Agg', warn=False)
import pylab

try: from mpl_toolkits.basemap import Basemap
except ImportError,e: warnings.warn(str(e))
try: from mpl_toolkits import axes_grid
except ImportError,e: warnings.warn(str(e))

from datetime import datetime
from glue import segments,git_version
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
from pylal import date,plotutils
from pylal.dq.dqTriggerUtils import def_get_time,get_column

__author__  = "Duncan Macleod <duncan.macleod@astro.cf.ac.uk>"
__version__ = "git id %s" % git_version.id
__date__    = git_version.date

"""
This module provides plotting routines for use in data quality investigations. All routines are written to work in as general a way as possible with ligolw tables and lsctables compatible columns, and to plot in a similar pythonic way to pylal.plotutils. 
"""

# =============================================================================
# Set plot parameters aux helper functions
# =============================================================================

def set_rcParams():

  """
    Customise the figure parameters.
  """

  # customise plot appearance
  pylab.rcParams.update({"text.usetex": True,
                         "text.verticalalignment": "center",
                         "lines.linewidth": 2,
                         "xtick.labelsize": 16,
                         "ytick.labelsize": 16,
                         "axes.titlesize": 22,
                         "axes.labelsize": 16,
                         "axes.linewidth": 1,
                         "grid.linewidth": 1,
                         "legend.fontsize": 16,
                         "legend.loc": "best",
                         "figure.figsize": [12,6],
                         "figure.dpi": 80,
                         "image.origin": 'lower',
                         "axes.grid": True,
                         "axes.axisbelow": False })

def set_ticks(ax, calendar_time=False):

  """
    Format the x- and y-axis ticks to ensure minor ticks appear when needed
    and the x-axis is set for spaces of 4 rather than 5.

    Arguments:

      ax : matplotlib.axes.AxesSubplot
        Axes object to format
  """

  #
  # make sure we get minor ticks if there are no major ticks in the range
  #

  # xticks
  ticks = list(ax.get_xticks())
  xlim  = ax.get_xlim()
  for i,tick in enumerate(ticks[::-1]):
    if not xlim[0] <= tick <= xlim[1]: ticks.pop(-1)
  if len(ticks)<=1:
    ax.xaxis.set_minor_formatter(pylab.matplotlib.ticker.ScalarFormatter())

  # yticks
  ticks = list(ax.get_yticks())
  ylim  = ax.get_ylim()
  for i,tick in enumerate(ticks[::-1]): 
    if not ylim[0] <= tick <= ylim[1]: ticks.pop(-1)
  if len(ticks)<=1:
    ax.yaxis.set_minor_formatter(pylab.matplotlib.ticker.ScalarFormatter())

  # set xticks in time format, python2.5 is not new enough for
  # flexibility, recoding part of AutoDateFormatter to get it
  if calendar_time:
    dateLocator = pylab.matplotlib.dates.AutoDateLocator()
    dateLocator.set_axis(ax.xaxis)
    dateLocator.refresh()
    scale = float( dateLocator._get_unit() )
    if ( scale == 365.0 ):
      dateFormatter = pylab.matplotlib.dates.DateFormatter("$%Y$")
    elif ( scale == 30.0 ):
      dateFormatter = pylab.matplotlib.dates.DateFormatter("$%y$/$%m$")
    elif ( (scale == 1.0) or (scale == 7.0) ):
      dateFormatter = pylab.matplotlib.dates.DateFormatter("$%m$/$%d$")
    elif ( scale == (1.0/24.0) ):
      dateFormatter = pylab.matplotlib.dates.DateFormatter("$%d$-$%H$")
    elif ( scale == (1.0/(24*60)) ):
      dateFormatter = pylab.matplotlib.dates.DateFormatter("$%H$:$%M$")
    elif ( scale == (1.0/(24*3600)) ):
      dateFormatter = pylab.matplotlib.dates.DateFormatter("$%H$:$%M$")
  
    ax.xaxis.set_major_locator(dateLocator)
    ax.xaxis.set_major_formatter(dateFormatter)

  else:
    # set xticks for 4 hours rather than 5
    xticks = ax.get_xticks()
    if len(xticks)>1 and xticks[1]-xticks[0]==5:
      ax.xaxis.set_major_locator(\
          pylab.matplotlib.ticker.MultipleLocator(base=2))

  
  # set tick linewidth
  for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markersize(10)
    line.set_markeredgewidth(1)

def display_name(columnName):

  """
    Format the string columnName (e.g. xml table column) into latex format for 
    an axis label.

    Examples:

    >>> display_name('snr')
    'SNR'

    >>> display_name('bank_chisq_dof')
    'Bank $\\chi^2$ DOF'

    Arguments:

      columnName : string
        string to format
  """

  acro  = ['snr', 'ra','dof', 'id', 'ms', 'far']
  greek = ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta',\
           'theta', 'iota', 'kappa', 'lamda', 'mu', 'nu', 'xi', 'omicron',\
           'pi', 'rho', 'sigma', 'tau', 'upsilon', 'phi', 'chi', 'psi', 'omega']
  unit  = ['ns']
  sub   = ['flow', 'fhigh', 'hrss', 'mtotal', 'mchirp']

  words = []
  for w in re.split('\s', columnName):
    if w.isupper(): words.append(w)
    else:           words.extend(re.split('_', w))

  # parse words
  for i,w in enumerate(words):
    # get acronym in lower case
    if w in acro:
      words[i] = w.upper()
    # get numerical unit
    elif w in unit:
      words[i] = '$(%s)$' % w
    # get character with subscript text
    elif w in sub:
      words[i] = '%s$_{\mbox{\\small %s}}$' % (w[0], w[1:])
    # get greek word
    elif w in greek:
      words[i] = '$\%s$' % w
    # get starting with greek word
    elif re.match('(%s)' % '|'.join(greek), w):
      if w[-1].isdigit():
        words[i] = '$\%s_{%s}$''' % tuple(re.findall(r"[a-zA-Z]+|\d+",w))
      elif w.endswith('sq'):
        words[i] = '$\%s^2$' % w.rstrip('sq')
    # get everything else
    else:
      if w.isupper():
        words[i] = w
      else:
        words[i] = w.title()
      # escape underscore
      words[i] = re.sub('(?<!\\\\)_', '\_', words[i])

  return ' '.join(words) 

def gps2datenum(gpstime):

  """
    Convert GPS time into floating in standard python time format
    (days since Jan 01 0000), don't seem to use properly leap seconds
  """
  
  # set time of 0 GPS in datenum units
  zeroGPS = 722820.0
  ## correct for leap seconds assuming that all time stamps are within
  # a range not including a leap
  # select the first time stamp
  if isinstance(gpstime,float) or isinstance(gpstime,int):
    repTime = gpstime
  else:
    if len(gpstime)>0:
      repTime = gpstime[0]
    else:
      return gpstime

  zeroGPS = zeroGPS  + float(date.XLALLeapSeconds(LIGOTimeGPS(0)) -\
                                 date.XLALLeapSeconds(LIGOTimeGPS(repTime)))/86400
  # convert to datenum (days)
  datenum = gpstime/86400 + zeroGPS

  return datenum

def calendar_time_unit(duration):

  if duration > 63072000:
    return "year"
  elif duration > 5184000:
    return "year/month"
  elif duration > 604800:
    return "month/day"
  elif duration > 7200:
    return "day-hour"
  elif duration > 600:
    return "hour:minute"
  else: 
    return "hour:minute:second"
  
def time_unit(duration):

  """
    Work out renormalisation for the time axis, makes the label more
    appropriate. Returns unit (in seconds) and string descriptor

    Example:

    >>> time_unit(100)
    (1, 'seconds')

    >>> time_unit(604800)
    (86400, 'days')

    Arguments:

      duration : float
        plot duration to normalise
  """

  # set plot time unit whether it's used or not
  if (duration) < 1000:
    unit = 1
  elif (duration) < 20000:
    unit = 60
  elif (duration) >= 20000 and (duration) < 604800:
    unit = 3600
  else:
    unit = 86400
  timestring = {1:'seconds', 60:'minutes',3600:'hours',86400:'days'}

  return unit, timestring[unit] 

def getTrigAttribute(trig, col):

  """
    Returns the value of trig.col or trig.get_col() for the given string col,
    and the object trig. If col='time' is given, trig.get_peak() is returned for
    *Burst* objects, trig.get_end() for *Inspiral* objects and trig.get_start()
    for *Ringdown* objects. Raises KeyError if cannot execute.

    Arguments:

      trig : [ lsctables.SnglBurst | lscatbles.SnglInspiral |
               lsctables.SnglRingdown ]
        xml table entry from which to extract parameter
      col : string
        glue.ligolw.table column name to extract
  """

  # return time
  if col=='time':
    if re.search('Burst', str(trig)):
      return trig.get_peak()
    elif re.search('Inspiral', str(trig)):
      return trig.get_end()
    elif re.search('Ringdown', str(trig)):
      return trig.get_start()
    else:
      return trig.time

  # return simple column 
  if col in trig.__slots__:
    return getattr(trig, col)

  # return get_XXX() parameter
  try:
    return eval('trig.get_%s()', col)

  # if we get here, panic
  except:
    raise KeyError, "Column '%s' not found in %s." % (col, type(trig))

def add_colorbar(self, mappable, cax=None, clim=None, log=False, base=10,\
                 label=None, **kwargs):

  # set axes for colorbar
  #if not cax:
  #  div = axes_grid.make_axes_locatable(self.ax)
  #  cax = div.new_horizontal(size="2%", pad=0.05, visible=True)
    #cax = div.append_axes("right", size="2%", pad=0.05)

  cmin, cmax = clim

  # construct colorbar tick formatter, using logic not supported before python
  # 2.5
  colorticksize = None
  if log and numpy.power(base,cmax)-numpy.power(base,cmin) > 4:
    formatter = pylab.matplotlib.ticker.FuncFormatter(lambda x,pos:\
                    numpy.power(base,x)>=1 and\
                        "$%d$" % round(numpy.power(base, x)) or\
                    numpy.power(base,x)>=0.01 and\
                        "$%.2f$" % numpy.power(base, x) or\
                    numpy.power(base,x)<0.01 and "$%f$")
  elif log and numpy.power(base,cmax)-numpy.power(base,cmin) > 0.4:
    formatter = pylab.matplotlib.ticker.FuncFormatter(lambda x,pos: "$%.2f$"\
                                                      % numpy.power(base, x))
  elif log:
    colorticksize = 'small'
    formatter = pylab.matplotlib.ticker.FuncFormatter(lambda x,pos: '$%.2e$'\
                                                      % numpy.power(base, x))
  elif not log and cmax-cmin > 4:
    formatter = pylab.matplotlib.ticker.FuncFormatter(lambda x,pos:\
                    x>=1 and "$%d$" % round(x) or\
                    x>=0.01 and "$%.2f$" % x or\
                    x==0.0 and "$0$" or\
                    x and "$%f$")
  else:
    formatter = None

  if clim:
    colorticks = numpy.linspace(cmin, cmax, 4)
  else:
    colorticks = None

  self.colorbar = self.ax.figure.colorbar(mappable, cax=cax, format=formatter,\
                                          ticks=colorticks, pad=0.005,\
                                          fraction=0.02, aspect=40, **kwargs)
  if colorticksize is not None:
    for l in self.colorbar.ax.yaxis.get_ticklabels():
      l.set_size(colorticksize)
#  colorbar = pylab.colorbar(mappable, cax=cax, format=formatter,\
                            #ticks=colorticks, **kwargs)
  self.colorbar.set_label(label)
  self.colorbar.draw_all()

  #try:
  #  self.colorbar = colorbar
  #except:
  #  return colorbar

def add_hidden_colorbar(self):

  # set axes for colorbar
  div = axes_grid.make_axes_locatable(self.ax)
  cax = div.new_horizontal(size="2%", pad=0.05, visible=False)

# =============================================================================
# Abstract classes for plots
# =============================================================================

# =============================================================================
# Class for segment plot 

class PlotSegmentsPlot(plotutils.BasicPlot):

  """
    Horizontal bar segment plot. Based originally on PlotSegmentsPlot class in
    pylal/bin/plotsegments.py
  """

  color_code = {'H1':'r', 'H2':'b', 'L1':'g', 'V1':'m', 'G1':'k'}

  def __init__(self, xlabel="", ylabel="", title="", subtitle="", t0=0, unit=1,\
               calendar_time=False):
    """
    Create a fresh plot.  Provide t0 to provide a reference time to use as
    zero.
    """
    plotutils.BasicPlot.__init__(self, xlabel, ylabel, title)
    self.ax.set_title(title, x=0.5, y=1.025)
    self.ax.text(0.5, 1.035, subtitle, horizontalalignment='center',
                 transform=self.ax.transAxes, verticalalignment='top')
    self.segdict = segments.segmentlistdict()
    self.keys = []
    if calendar_time:
      self._time_transform =\
          lambda seg: segments.segment(gps2datenum(float(seg[0])),\
                                       gps2datenum(float(seg[1])))
    else:
      self._time_transform =\
          lambda seg: segments.segment(float(seg[0]-t0)/unit,\
                                       float(seg[1]-t0)/unit)



  def add_content(self, segdict, keys=None, t0=0, unit=1):
    if not keys:
      keys = sorted(segdict.keys())
    for key in keys:
      self.segdict[key] = segdict[key]
    self.keys.extend(keys)

  def highlight_segment(self, seg, **plot_args):
    """
    Highlight a particular segment with dashed lines.
    """
    a,b = self._time_transform(seg)
    plot_args.setdefault('linestyle', '--')
    plot_args.setdefault('color','r')
    self.ax.axvline(a, **plot_args)
    self.ax.axvline(b, **plot_args)

  @plotutils.method_callable_once
  def finalize(self, labels_inset=False, hidden_colorbar=False):

    for row,key in enumerate(self.keys):
      if self.color_code.has_key(key):
        col = self.color_code[key]
      else:
        col = 'b'
      for seg in self.segdict[key]:
        a,b = self._time_transform(seg)
        self.ax.fill([a, b, b, a, a],\
                     [row-0.4, row-0.4, row+0.4, row+0.4, row-0.4], 'b')
      if labels_inset:
        self.ax.text(0.01,(row+1)/(len(self.keys)+1), re.sub('\\+_+','\_',key),\
                     horizontalalignment='left', verticalalignment='center',\
                     transform=self.ax.transAxes, backgroundcolor='white',\
                     bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))

    ticks = pylab.arange(len(self.keys))
    self.ax.set_yticks(ticks)
    if labels_inset:
      self.ax.set_yticklabels(ticks, color='white')
    else:
      self.ax.set_yticklabels([re.sub(r'\\+_+', '\_', k)\
                               for k in self.keys], size='small')
    self.ax.set_ylim(-1, len(self.keys))

    # add hidden colorbar
    if hidden_colorbar:
      self.add_hidden_colorbar()

PlotSegmentsPlot.add_hidden_colorbar = add_hidden_colorbar

# =============================================================================
# Class for standard scatter plot

class ScatterPlot(plotutils.BasicPlot):

  """
    A simple scatter plot, taking x- and y-axis data.
  """

  def __init__(self, xlabel="", ylabel="", title="", subtitle=""):
    plotutils.BasicPlot.__init__(self, xlabel, ylabel, title)
    self.ax.set_title(title, x=0.5, y=1.025)
    self.ax.text(0.5, 1.035, subtitle, horizontalalignment='center',
                 transform=self.ax.transAxes, verticalalignment='top')
    self.x_data_sets = []
    self.y_data_sets = []
    self.kwarg_sets = []

  def add_content(self, x_data, y_data, **kwargs):
    self.x_data_sets.append(x_data)
    self.y_data_sets.append(y_data)
    self.kwarg_sets.append(kwargs)

  @plotutils.method_callable_once
  def finalize(self, loc='best', logx=False, logy=False, hidden_colorbar=False):
    # make plot
    for x_vals, y_vals, plot_kwargs, c in \
        itertools.izip(self.x_data_sets, self.y_data_sets, self.kwarg_sets,\
                       plotutils.default_colors()):
      plot_kwargs.setdefault("marker", "o")
      plot_kwargs.setdefault("s", 20)
      plot_kwargs.setdefault("linewidth", 1)

      self.ax.scatter(x_vals, y_vals, c=c, **plot_kwargs)

    # add hidden colorbar to make plots the right aspect
    if hidden_colorbar:
      self.add_hidden_colorbar()

    # add legend if there are any non-trivial labels
    self.add_legend_if_labels_exist(loc=loc)

    leg = self.ax.get_legend()
    # set transparent legend
    if leg:
      legfr = leg.get_frame()
      legfr.set_alpha(0.5)

    # set axes
    if logx:
      self.ax.set_xscale('log')
    if logy:
      self.ax.set_yscale('log')

ScatterPlot.add_hidden_colorbar = add_hidden_colorbar

# =============================================================================
# Class for scatter plot with colorbar

class ColorbarScatterPlot(plotutils.BasicPlot):

  """
    A scatter plot of x- versus y-data, coloured by z-data. A single colorbar
    is used, from the union of the ranges of z-data for each content set,
    unless strict limits are passed as clim=(cmin, cmax) to finalize().
  """

  def __init__(self, xlabel="", ylabel="", zlabel="", title="", subtitle=""):
    plotutils.BasicPlot.__init__(self, xlabel, ylabel, title)
    self.ax.set_title(title, x=0.5, y=1.025)
    self.ax.text(0.5, 1.035, subtitle, horizontalalignment='center',
                 transform=self.ax.transAxes, verticalalignment='top')

    self.x_data_sets = []
    self.y_data_sets = []
    self.z_data_sets = []
    self.rank_data_sets = []
    self.kwarg_sets = []
    self.color_label = zlabel

  def add_content(self, x_data, y_data, z_data, rank_data, **kwargs):
    self.x_data_sets.append(x_data)
    self.y_data_sets.append(y_data)
    self.z_data_sets.append(z_data)
    self.rank_data_sets.append(rank_data)
    self.kwarg_sets.append(kwargs)

  @plotutils.method_callable_once
  def finalize(self, loc='best', logx=False, logy=False, logz=False, clim=None,\
               base=10):
  
    # set up things we'll need later
    p = []

    # set colorbar limits
    if clim:
      cmin,cmax = clim
    else:
      try:
        cmin = min([min(z_vals) for z_vals in self.z_data_sets])*0.99
        cmax = max([max(z_vals) for z_vals in self.z_data_sets])*1.01
      # catch no triggers
      except ValueError:
        cmin = 1
        cmax = 10

    # reset logs
    if logz:
      cmin = numpy.math.log(cmin, base)
      cmax = numpy.math.log(cmax, base)

    # make plot
    for x_vals, y_vals, z_vals, rank_vals, plot_kwargs in\
        itertools.izip(self.x_data_sets, self.y_data_sets, self.z_data_sets,\
                       self.rank_data_sets, self.kwarg_sets):

      if logz:  z_vals = [numpy.math.log(z, base) for z in z_vals]

      plot_kwargs.setdefault("marker", "o")
      plot_kwargs.setdefault("vmin", cmin)
      plot_kwargs.setdefault("vmax", cmax)

      # sort data by z-value
      zipped = zip(x_vals, y_vals, z_vals, rank_vals)
      zipped.sort(key=lambda (x,y,z,r): r)
      x_vals, y_vals, z_vals, rank_vals = map(list, zip(*zipped))

      p.append(self.ax.scatter(x_vals, y_vals, c=z_vals, **plot_kwargs))

    if len(p)<1:
      p.append(self.ax.scatter([1], [1], c=[cmin], vmin=cmin,\
                               vmax=cmax, visible=False))


    # write colorbar
    #self.set_colorbar(p[-1], [cmin, cmax], logz, base)
    self.add_colorbar(p[-1], clim=[cmin, cmax], log=logz,\
                      label=self.color_label, base=10)

    # set axes
    if logx:
      self.ax.set_xscale('log')
    if logy:
      self.ax.set_yscale('log')

    self.add_legend_if_labels_exist(loc=loc)

    # set transparent legend
    leg = self.ax.get_legend()
    if leg:
      legfr = leg.get_frame()
      legfr.set_alpha(0.5)

  def set_colorbar(self, mappable, clim=None, log=True, base=10): 

    cmin, cmax = clim

    # construct colorbar tick formatter, using logic not supported before python
    # 2.5
    if log and numpy.power(base,cmax)-numpy.power(base,cmin) > 4:
      formatter = pylab.matplotlib.ticker.FuncFormatter(lambda x,pos:\
                      numpy.power(base,x)>=1 and\
                          "$%d$" % round(numpy.power(base, x)) or\
                      numpy.power(base,x)>=0.01 and\
                          "$%.2f$" % numpy.power(base, x) or\
                      numpy.power(base,x)<0.01 and "$%f$")
    elif log and numpy.power(base,cmax)-numpy.power(base,cmin) > 0.4:
      formatter = pylab.matplotlib.ticker.FuncFormatter(lambda x,pos: "$%.2f$"\
                                                        % numpy.power(base, x))
    elif log:
      formatter = pylab.matplotlib.ticker.FuncFormatter(lambda x,pos: "$%.4e$"\
                                                        % numpy.power(base, x))
    else:
      formatter = pylab.matplotlib.ticker.FuncFormatter(lambda x,pos: "$%.4e$"\
                                                        % x)

    if clim:
      colorticks = numpy.linspace(cmin, cmax, 4)
      self.colorbar = self.ax.figure.colorbar(mappable, format=formatter,\
                                   ticks=colorticks)
    else:
      self.colorbar = self.ax.figure.colorbar(mappable, format=formatter)

    self.colorbar.set_label(self.color_label)
    self.colorbar.draw_all()

ColorbarScatterPlot.add_colorbar = add_colorbar

# =============================================================================
# Extension of ColorbarScatterPlot to plot DetChar-style scatter plot

class DetCharScatterPlot(ColorbarScatterPlot):

  """
    A 'DetChar' style scatter plot, whereby those triggers under a threshold
    on the colour column are plotted much smaller than others, allowing line
    features to be shown easily. 
  """

  @plotutils.method_callable_once
  def finalize(self, loc='best', logx=False, logy=False, logz=True, base=10,\
               clim=None, rankthreshold=None):

    p = []
    # set colorbar limits
    if clim:
      cmin,cmax = clim
    else:
      try:
        cmin = min([min(z_vals) for z_vals in self.z_data_sets])*0.99
        cmax = max([max(z_vals) for z_vals in self.z_data_sets])*1.01
      # catch no triggers
      except ValueError:
        cmin = 1
        cmax = 10 

    # reset logs
    if logz:
      cmin = numpy.math.log(cmin, base)
      cmax = numpy.math.log(cmax, base)

    # make plot
    for x_vals, y_vals, z_vals, rank_vals, plot_kwargs in\
        itertools.izip(self.x_data_sets, self.y_data_sets, self.z_data_sets,\
                       self.rank_data_sets, self.kwarg_sets):
      plot_kwargs.setdefault("vmin", cmin)
      plot_kwargs.setdefault("vmax", cmax)
      plot_kwargs.setdefault("s", 15)
      plot_kwargs.setdefault("marker", "o")
      plot_kwargs.setdefault("edgecolor", "k")

      if logz:
        z_vals = [numpy.math.log(z, base) for z in z_vals]

      # sort data by z-value
      zipped = zip(x_vals, y_vals, z_vals, rank_vals)
      zipped.sort(key=lambda (x,y,z,r): r)
      x_vals, y_vals, z_vals, rank_vals = map(list, zip(*zipped))

      bins = ['low', 'high']
      x_bins = {}
      y_bins = {}
      z_bins = {}
      for bin in bins:
        x_bins[str(bin)] = []
        y_bins[str(bin)] = []
        z_bins[str(bin)] = []
      for i in xrange(len(x_vals)):
        if rankthreshold and rank_vals[i] < rankthreshold:
          x_bins[bins[0]].append(float(x_vals[i]))
          y_bins[bins[0]].append(float(y_vals[i]))
          z_bins[bins[0]].append(float(z_vals[i]))
        else:
          x_bins[bins[1]].append(float(x_vals[i]))
          y_bins[bins[1]].append(float(y_vals[i]))
          z_bins[bins[1]].append(float(z_vals[i]))

      # plot bins
      for i,bin in enumerate(bins):
        if bin == bins[0]:
          args = copy.deepcopy(plot_kwargs)
          args['s']/=4
          if not (args.has_key('marker') and args['marker']=='x'):
            args['edgecolor'] = 'none'
          if len(x_bins[bin])>=1:
            p.append(self.ax.scatter(x_bins[bin], y_bins[bin], c=z_bins[bin],\
                                     **args))
        else:
          if len(x_bins[bin])>=1:
            p.append(self.ax.scatter(x_bins[bin], y_bins[bin], c=z_bins[bin],\
                                     **plot_kwargs))
      args = plot_kwargs

    if len(p)<1:
      p.append(self.ax.scatter([1], [1], c=[cmin], vmin=cmin,\
                               vmax=cmax, visible=False))

    # write colorbar
    self.add_colorbar(p[-1], clim=[cmin, cmax], log=logz,\
                      label=self.color_label, base=10)
    #self.set_colorbar(p[-1], [cmin, cmax], logz, base)

    # set axes
    if logx:
      self.ax.set_xscale('log')
    if logy:
      self.ax.set_yscale('log')

    self.add_legend_if_labels_exist(loc=loc)

# =============================================================================
# Class for line histogram

class LineHistogram(ColorbarScatterPlot, plotutils.BasicPlot):

  """
    A simple line histogram plot. The values of each histogram bin are plotted
    using pylab.plot(), with points centred on the x values and height equal
    to the y values.

    Cumulative, and rate options can be passeed to the finalize() method to
    format each trace individually.

  """

  def __init__(self, xlabel="", ylabel="", title="", subtitle="",\
               colorlabel=""):
    plotutils.BasicPlot.__init__(self, xlabel, ylabel, title)
    self.ax.set_title(title, x=0.5, y=1.025)
    self.ax.text(0.5, 1.035, subtitle, horizontalalignment='center',
                 transform=self.ax.transAxes, verticalalignment='top')
    self.data_sets = []
    self.livetimes = []
    self.kwarg_sets = []
    self.color_label = colorlabel
    self.color_values = []

  def add_content(self, data, livetime=1, cval=None, **kwargs):
    self.data_sets.append(data)
    self.livetimes.append(livetime)
    self.kwarg_sets.append(kwargs)
    self.color_values.append(cval)

  @plotutils.method_callable_once
  def finalize(self, loc='best', num_bins=100, cumulative=False, rate=False,\
               logx=False, logy=False, fill=False, base=10, xlim=None,\
               colorbar=None, clim=None, hidden_colorbar=False):

    # determine binning
    if xlim:
      min_stat, max_stat = xlim
    else:
      min_stat, max_stat = plotutils.determine_common_bin_limits(self.data_sets)
    
    if min_stat!=max_stat:
      max_stat *= 1.001
      if logx:
        bins = numpy.logspace(numpy.math.log(min_stat, base),\
                              numpy.math.log(max_stat, base),\
                              num_bins+1, endpoint=True)
      else:
        bins = numpy.linspace(min_stat, max_stat, num_bins + 1, endpoint=True)
    else:
      bins = [min_stat]

    if logy:
      ymin = 1E-100
    else:
      ymin = 0

    # get colors
    if colorbar:
      if isinstance(colorbar, str):
        ctype = colorbar
        cmap = pylab.matplotlib.colors.LinearSegmentedColormap('clrs',\
                                           pylab.matplotlib.cm.jet._segmentdata)
      else:
        ctype, cmap = colorbar
      assert re.match('(lin|log)', ctype, re.I),\
             "colorbar must have type 'linear', or 'log'"
      try:
        colors = pylab.matplotlib.colors.makeMappingArray(100000,\
                                                          cmap)
      except IndexError:
        xind = numpy.linspace(0, 1, 100000)
        colors = numpy.clip(numpy.array(cmap(xind), dtype=numpy.float), 0, 1)
      
      if clim:
        cmin,cmax = clim
      else:
        try:
          cmin = min(self.color_values)
          cmax = max(self.color_values)
        except ValueError:
          cmin = 1
          cmax = 10

      if re.search('log', ctype, re.I):
        cmin = numpy.math.log(cmin, base)
        cmax = numpy.math.log(cmax, base)
        self.color_values = [numpy.math.log(x, base) for x in self.color_values]
      color_idx = lambda x: int((x-cmin)/(cmax-cmin)*(len(colors)-1))
      self.color_values = [colors[color_idx(x)] for x in self.color_values]

    # add hidden colorbar to make plots the right aspect
    elif hidden_colorbar:
      self.add_hidden_colorbar()
 
    p = []
    for i, (data_set, livetime, plot_kwargs, col) in\
        enumerate(itertools.izip(self.data_sets, self.livetimes,\
                  self.kwarg_sets, self.color_values)):

      #
      # make histogram
      #

      # get version
      v = [int(i) for i in numpy.version.version.split('.')]
      if v[1] < 1:
        y, x = numpy.histogram(data_set, bins=bins)
      elif v[1] < 3:
        y, x = numpy.histogram(data_set, bins=bins, new=True)
        x = x[:-1]
      else:
        y, x = numpy.histogram(data_set, bins=bins)
        x = x[:-1]

      # get cumulative sum
      if cumulative:
        y = y[::-1].cumsum()[::-1]

      # convert to rate
      if rate:
        if livetime > 0:
          y = y/livetime
          ymin /= livetime
        else:
          y = numpy.empty(y.shape)
          y.fill(numpy.nan)
          ymin = numpy.nan

      # reset zeros on logscale, tried with numpy, unreliable
      if logy:
        y = y.astype(float)
        numpy.putmask(y, y==0, ymin)

      # set color
      if colorbar:
        plot_kwargs['color'] = col

      # plot
      if fill:
        plot_kwargs.setdefault('linewidth', 1)
        plot_kwargs.setdefault('alpha', 0.8)
        p.append(self.ax.plot(x, y, **plot_kwargs))
        self.ax.fill_between(x, 1E-100, y, **plot_kwargs)
      else:
        p.append(self.ax.plot(x, y, **plot_kwargs))

    # set axes
    if logx:
      self.ax.set_xscale('log')
    if logy:
      self.ax.set_yscale('log')

    # set colorbar
    if colorbar:
      self.add_colorbar(self.ax.scatter([1], [1], c=1, cmap=cmap,\
                                          vmin=cmin, vmax=cmax, visible=False),\
                          clim=[cmin, cmax], label=self.color_label)
    elif hidden_colorbar:
      self.add_hidden_colorbar()

    # add legend if there are any non-trivial labels
    self.add_legend_if_labels_exist(loc=loc)

    # fix legend
    leg = self.ax.get_legend()
    if leg:
      for l in leg.get_lines():
        l.set_linewidth(4)

LineHistogram.add_hidden_colorbar = add_hidden_colorbar

# =============================================================================
# Extension of VerticalBarHistogram to include log axes

class VerticalBarHistogram(plotutils.VerticalBarHistogram):

  def __init__(self, xlabel="", ylabel="", title="", subtitle=""):
    plotutils.BasicPlot.__init__(self, xlabel, ylabel, title)
    self.ax.set_title(title, x=0.5, y=1.025)
    self.ax.text(0.5, 1.035, subtitle, horizontalalignment='center',
                 transform=self.ax.transAxes, verticalalignment='top')
    self.data_sets = []
    self.kwarg_sets = []

  @plotutils.method_callable_once
  def finalize(self, num_bins=20, normed=False, logx=False, logy=False,\
               base=10, hidden_colorbar=False, loc='best'):

    # determine binning
    min_stat, max_stat = plotutils.determine_common_bin_limits(self.data_sets)
    if logx and min_stat!=0 and max_stat!=0:
      bins = numpy.logspace(numpy.math.log(min_stat, base),\
                            numpy.math.log(max_stat, base),\
                            num_bins+1, endpoint=True)
    else:
      bins = numpy.linspace(min_stat, max_stat, num_bins + 1, endpoint=True)

    # determine bar width; gets silly for more than a few data sets
    if logx:
      width = [bins[i+1]-bins[i] for i in xrange(len(bins)-1)]
      width.append(width[-1])
    else:
      width = (1 - 0.1 * len(self.data_sets)) * (bins[1] - bins[0])

    width = numpy.asarray(width)/2

    # get version
    v = map(int, numpy.version.version.split('.'))
    if v[1] >= 1: width = width[:-1]

    # set base of plot in log scale
    if logy:
      ymin = (base**-1)*5
    else:
      ymin = 0

    # make plot
    legends = []
    plot_list = []
    for i, (data_set, plot_kwargs) in \
        enumerate(itertools.izip(self.data_sets, self.kwarg_sets)):
      # set default values
      plot_kwargs.setdefault("alpha", 0.6)
      plot_kwargs.setdefault("align", "center")
      plot_kwargs.setdefault("width", width)
      if logy:
        plot_kwargs.setdefault("bottom", ymin)
      else:
        plot_kwargs.setdefault("bottom", ymin)

      # make histogram
      if v[1] < 1:
        y, x = numpy.histogram(data_set, bins=bins, normed=normed)
      elif v[1] < 3:
        y, x = numpy.histogram(data_set, bins=bins, new=True, normed=normed)
        x = x[:-1]
      else:
        y, x = numpy.histogram(data_set, bins=bins, normed=normed)
        x = x[:-1]

      if logy:
        y = y-ymin

      # stagger bins for pure aesthetics
      #x += 0.1 * i * max_stat / num_bins

      # plot
      self.ax.bar(x, y, **plot_kwargs)

    # set axes
    if logx:
      self.ax.set_xscale('log')
    if logy:
      self.ax.set_yscale('log')

    # set logy minimum
    self.ax.set_ybound(lower=ymin)

    # add legend if there are any non-trivial labels
    self.add_legend_if_labels_exist(loc=loc)
    # set transparent legend
    leg = self.ax.get_legend()
    if leg: leg.get_frame().set_alpha(0.5)


    # add hidden colorbar
    if hidden_colorbar:
      self.add_hidden_colorbar()

VerticalBarHistogram.add_hidden_colorbar = add_hidden_colorbar

# =============================================================================
# Class for time series plot

class DataPlot(plotutils.BasicPlot):

  """
    Time-series data plot. Just a nice looking line plot.
  """

  def __init__(self, xlabel="", ylabel="", title="", subtitle=""):
    plotutils.BasicPlot.__init__(self, xlabel, ylabel, title)
    self.ax.set_title(title, x=0.5, y=1.025)
    self.ax.text(0.5, 1.035, subtitle, horizontalalignment='center',
                 transform=self.ax.transAxes, verticalalignment='top')
    self.x_data_sets = []
    self.y_data_sets = []
    self.kwarg_sets = []

  def add_content(self, x_data, y_data, **kwargs):
    self.x_data_sets.append(numpy.ma.asarray(x_data))
    self.y_data_sets.append(numpy.ma.asarray(y_data))
    self.kwarg_sets.append(kwargs)

  @plotutils.method_callable_once
  def finalize(self, loc='best', logx=False, logy=False, hidden_colorbar=False):
    # make plot
    plots = []
    markersizes = []
    markerscale = 4

    for x_vals, y_vals, plot_kwargs, c in \
        itertools.izip(self.x_data_sets, self.y_data_sets,\
                       self.kwarg_sets, plotutils.default_colors()):

      # magnify the markers on the legend
      plot_kwargs.setdefault('markersize', 5)
      markersizes.append(plot_kwargs['markersize'])
      plot_kwargs['markersize'] = min(20, plot_kwargs['markersize']*markerscale)
      # plot
      if logx:
        x_vals = x_vals.astype(float)
        numpy.putmask(x_vals, x_vals==0, 1E-100)
      if logy:
        y_vals = y_vals.astype(float)
        numpy.putmask(y_vals, y_vals==0, 1E-100)

      plots.append(self.ax.plot(x_vals, y_vals, **plot_kwargs))

    # add legend if there are any non-trivial labels
    self.add_legend_if_labels_exist(loc=loc)
    leg = self.ax.get_legend()
    # magnify the lines on the legend 
    if leg:
      try:
        for l in leg.get_lines():
          if l.get_linewidth():
            l.set_linewidth(4)
      except AttributeError:
        pass

    # reset markersizes on plot
    for i,plot in enumerate(plots):
      l = plot[0] 
      l.set_markersize(markersizes[i])

    # set transparent legend
    if leg:
      legfr = leg.get_frame()
      legfr.set_alpha(0.5)

    # add hidden colorbar
    if hidden_colorbar: self.add_hidden_colorbar()

    # set axes
    if logx:  self.ax.set_xscale('log')
    if logy:  self.ax.set_yscale('log')

DataPlot.add_hidden_colorbar = add_hidden_colorbar

# =============================================================================
# Class for color map plot

class ColorMap(ColorbarScatterPlot):

  def __init__(self, xlabel="", ylabel="", zlabel="", title="", subtitle=""):
    plotutils.BasicPlot.__init__(self, xlabel, ylabel, title)
    self.color_label = zlabel
    self.ax.set_title(title, x=0.5, y=1.025)
    self.ax.text(0.5, 1.035, subtitle, horizontalalignment='center',
                 transform=self.ax.transAxes, verticalalignment='top')
    self.data_sets   = []
    self.extent_sets = []
    self.kwarg_sets  = []

  def add_content(self, data, extent, **kwargs):
    self.data_sets.append(numpy.asarray(data))
    self.extent_sets.append(extent)
    self.kwarg_sets.append(kwargs)

  @plotutils.method_callable_once
  def finalize(self, loc='best', logx=False, logy=False, logz=False, clim=None,\
               origin=pylab.rcParams['image.origin'], base=10):

    # set colorbar limits
    if clim:
      cmin,cmax = clim
    else:
      cmin = min([z_vals.min() for z_vals in self.data_sets])*0.99
      cmax = max([z_vals.max() for z_vals in self.data_sets])*1.01

    # reset logs
    if logz:
      cmin = numpy.math.log(cmin, base)
      cmax = numpy.math.log(cmax, base)

    p = []

    for i, (data_set, extent, plot_kwargs) in \
        enumerate(itertools.izip(self.data_sets, self.extent_sets,\
                  self.kwarg_sets)):
 
      plot_kwargs.setdefault("vmin", cmin)
      plot_kwargs.setdefault("vmax", cmax)
      plot_kwargs.setdefault("norm", None)
      plot_kwargs.setdefault("interpolation", "kaiser")

      if logz:
        if base==10:
          data_set = numpy.log10(data_set)
        elif base==numpy.e:
          data_set = numpy.log(data_set)
        elif base==2:
          data_set = numpy.log2(data_set)
        else:
          raise AttributeError("Can only use base = 2, e, or 10.")

      p.append(self.ax.imshow(data_set, extent=extent, origin=origin,\
                              **plot_kwargs))

    if len(p)==0:
      p.append(self.ax.imshow([[1]], vmin=cmin,\
                               vmax=cmax, visible=False))

    # write colorbar
    self.add_colorbar(p[-1], clim=[cmin, cmax], log=logz,\
                      label=self.color_label, base=10)

    # set axes
    if logx:  self.ax.set_xscale('log')
    if logy:  self.ax.set_yscale('log')

# =============================================================================
# Projection for sky position plotting

class SkyPositionsPlot(ScatterPlot):

  """
    Plot of sky positions plotted onto mpl_toolkits basemap projection.
  """

  @plotutils.method_callable_once
  def finalize(self, loc='lower right', projection='moll', centre=None,\
               range=[(0, 0), (1, 1)]):

    """
      Finalize and draw plot. Can only be called once.

      Keyword arguments:

        toc : [ str | int ]
          compatible value for legend loc, default: 'lower right'
        projection : str
          valid projection name for mpl_toolkits.basemap.Basemap
        centre : list
          (ra, dec) in degrees for point to centre on plot, only works for
          some projections
        range : list
          [(xmin, ymin), (xmax, ymax)] pair of tuples for lower left and upper
          right corners of plot area. Full areai (entire sphere) is
          [(0,0), (1,1)]. Narrower range is zoomed in, wider range zoomed out.
    """

    # set up projection
    projection = projection.lower()
    if not centre:
      centre = (0,0)
    centre = [(-(centre[0]-180)+180) % 360, centre[1]]
    m1 = Basemap(projection=projection, lon_0=centre[0], lat_0=centre[1],\
                 resolution=None, ax=self.ax)

    # set up range
    width  = m1.urcrnrx
    height = m1.urcrnry
    range    = [list(range[0]), list(range[1])]
    maprange = [[-width/2, -height/2], [width/2, height/2]]
    if range[0][0] != None: maprange[0][0] = width/2*(range[0][0]-1)
    if range[0][1] != None: maprange[0][1] = height/2*(range[0][1]-1)
    if range[1][0] != None: maprange[1][0] = width/2*(range[1][0])
    if range[1][1] != None: maprange[1][1] = height/2*(range[1][1])

    # set projection
    m  = Basemap(projection=projection, lon_0=centre[0], lat_0=centre[1],\
                 llcrnrx=maprange[0][0], llcrnry=maprange[0][1],\
                 urcrnrx=maprange[1][0], urcrnry=maprange[1][1],\
                 resolution=None, ax=self.ax)

    # turn on 'fulldisk' property when using full disk
    if maprange == [[-width/2, -height/2], [width/2, height/2]]:
      m._fulldisk = True

    xrange = (m.llcrnrx, m.urcrnrx)
    yrange = (m.llcrnry, m.urcrnry) 

    # make plot
    for x_vals, y_vals, plot_kwargs, c in\
        itertools.izip(self.x_data_sets, self.y_data_sets, self.kwarg_sets,\
                       plotutils.default_colors()):

      plot_kwargs.setdefault("marker", "o")
      plot_kwargs.setdefault("edgecolor", "k")
      plot_kwargs.setdefault("facecolor", c)

      # project data
      convert = lambda x: (x>=180 and x-360) or (x<180 and x)
      x_vals  = [-convert(x) for x in x_vals]
      if projection in ['moll', 'hammer', 'orth']:
        convert = lambda x: (x>=0 and x) or (x<0 and x+360)
      x_vals, y_vals = m(x_vals, y_vals)
      m.scatter(x_vals, y_vals, **plot_kwargs)

    # finish projection
    m.drawmapboundary()
   
    # set labels
    if projection in ['ortho']:
      plabels = [0, 0, 0, 0]
      mlabels = [0, 0, 0, 0]
    else:
      plabels = [1, 0, 0, 0]
      mlabels = [0, 0, 0, 0]

    # draw parallels
    parallels = numpy.arange(-90., 120., 30.)
    m.drawparallels(parallels, labels=plabels,\
                    labelstyle='+/-', latmax=90)

    # draw meridians
    if projection in ['moll', 'hammer', 'ortho']:
      meridians = numpy.arange(0., 360., 45.)
    else:
      meridians = numpy.arange(-180, 181, 45)
    m.drawmeridians(meridians, labels=mlabels,\
                    latmax=90, labelstyle='+/-')

    # label parallels for certain projections
    if projection in ['ortho']:
      for lon,lat in zip([0.]*len(parallels), parallels):
        #if lonrange[0]<=lon<=lonrange[1] and latrange[0]<=lat<=latrange[1]:
        x, y = m(lon, lat)
        lon, lat = m1(x, y, inverse=True)
        if x<=10**20 and y<=10**20\
        and xrange[0]<x<xrange[1] and yrange[0]<=y<=yrange[1]:
          m.ax.text(x, y, r"$%0.0f^\circ$" % lat)

    # label meridians for all projections
    for lon,lat in zip(meridians, [0.]*len(meridians)):
      tlon = (-(lon-180)+180) % 360
      x, y = m(lon, lat)
      lon, lat = m1(x, y, inverse=True)

      if x<=10**20 and y<=10**20\
      and xrange[0]<x<xrange[1] and yrange[0]<=y<=yrange[1]:
        m.ax.text(x, y, r"$%0.0f^\circ$" % tlon)

    # set legend
    self.add_legend_if_labels_exist(loc=loc, scatterpoints=1)

class ColorbarSkyPositionsPlot(ColorbarScatterPlot):

  @plotutils.method_callable_once
  def finalize(self, loc='lower right', projection='moll', centre=None,\
               range=[(0, 0), (1, 1)], logz=False, clim=None,\
               base=10):

    """
      Finalize and draw plot. Can only be called once.

      Keyword arguments:

        toc : [ str | int ]
          compatible value for legend loc, default: 'lower right'
        projection : str
          valid projection name for mpl_toolkits.basemap.Basemap
        centre : list
          (ra, dec) in degrees for point to centre on plot, only works for
          some projections
        range : list
          [(xmin, ymin), (xmax, ymax)] pair of tuples for lower left and upper
          right corners of plot area. Full area (entire sphere) is
          [(0,0), (1,1)]. Narrower range is zoomed in, wider range zoomed out.
        logz : [ True | False ]
          plot z-axis in log scale
        clim : 
    """

    p = []

    # set up projection
    projection = projection.lower()
    if not centre:
      centre = (0,0)
    m1 = Basemap(projection=projection, lon_0=centre[0], lat_0=centre[1],\
                 resolution=None, ax=self.ax)

    # set up range
    width  = m1.urcrnrx
    height = m1.urcrnry
    range    = [list(range[0]), list(range[1])]
    maprange = [[-width/2, -height/2], [width/2, height/2]]
    if range[0][0] != None: maprange[0][0] = width/2*(range[0][0]-1)
    if range[0][1] != None: maprange[0][1] = height/2*(range[0][1]-1)
    if range[1][0] != None: maprange[1][0] = width/2*(range[1][0])
    if range[1][1] != None: maprange[1][1] = height/2*(range[1][1])

    # set projection
    m  = Basemap(projection=projection, lon_0=centre[0], lat_0=centre[1],\
                 llcrnrx=maprange[0][0], llcrnry=maprange[0][1],\
                 urcrnrx=maprange[1][0], urcrnry=maprange[1][1],\
                 resolution=None, ax=self.ax)

    # turn on 'fulldisk' property when using full disk
    if maprange == [[-width/2, -height/2], [width/2, height/2]]:
      m._fulldisk = True

    xrange = (m.llcrnrx, m.urcrnrx)
    yrange = (m.llcrnry, m.urcrnry)

    # set colorbar limits
    if clim:
      cmin,cmax = clim
    else:
      try:
        cmin = min([min(z_vals) for z_vals in self.z_data_sets])*0.99
        cmax = max([max(z_vals) for z_vals in self.z_data_sets])*1.01
      # catch no triggers
      except ValueError:
        cmin = 1
        cmax = 10

    # reset logs
    if logz:
      cmin = numpy.math.log(cmin, base)
      cmax = numpy.math.log(cmax, base)

    # make plot
    for x_vals, y_vals, z_vals, plot_kwargs in\
        itertools.izip(self.x_data_sets, self.y_data_sets, self.z_data_sets,\
                       self.kwarg_sets):

      if logz:  z_vals = [numpy.math.log(z, base) for z in z_vals]

      plot_kwargs.setdefault("marker", "o")
      plot_kwargs.setdefault("edgecolor", "k")
      plot_kwargs.setdefault("vmin", cmin)
      plot_kwargs.setdefault("vmax", cmax)

      # sort data by z-value
      zipped = zip(x_vals, y_vals, z_vals)
      zipped.sort(key=lambda (x,y,z): z)
      x_vals, y_vals, z_vals = map(list, zip(*zipped))

      # project data
      if projection in ['moll', 'hammer', 'orth']:
        convert = lambda x: (x>=0 and x) or (x<0 and x+360)
      else:
        convert = lambda x: (x>=180 and x-360) or (x<180 and x)
      x_vals  = [convert(x) for x in x_vals]
      x_vals, y_vals = m(x_vals, y_vals)

      # pull out z value if given
      if plot_kwargs.has_key('facecolor'):
        self.ax.scatter(x_vals, y_vals, **plot_kwargs)
      else:
        p.append(self.ax.scatter(x_vals, y_vals, c=z_vals, **plot_kwargs))

    # finish projection
    m.drawmapboundary()

    # set labels
    if projection in ['ortho']:
      plabels = [0, 0, 0, 0]
      mlabels = [0, 0, 0, 0]
    else:
      plabels = [1, 0, 0, 0]
      mlabels = [0, 0, 0, 1]

    # draw parallels
    parallels = numpy.arange(-90., 120., 30.)
    m.drawparallels(parallels, labels=plabels,\
                    labelstyle='+/-', latmax=90)

    # draw meridians
    if projection in ['moll', 'hammer', 'ortho']:
      meridians = numpy.arange(0., 360., 45.)
    else:
      meridians = numpy.arange(-180, 181, 45)
    m.drawmeridians(meridians, labels=mlabels,\
                    latmax=90, labelstyle='+/-')

    # label parallels for certain projections
    if projection in ['ortho']:
      for lon,lat in zip([0.]*len(parallels), parallels):
        x, y = m(lon, lat)
        if x<=10**20 and y<=10**20\
        and xrange[0]<x<xrange[1] and yrange[0]<=y<=yrange[1]:
          self.ax.text(x, y, r"$%0.0f^\circ$" % lat)
    # label meridians for certain projections
    if projection in ['moll', 'hammer', 'ortho']:
      for lon,lat in zip(meridians, [0.]*len(meridians)):
        tlon = lon
        while tlon < 0:  tlon += 360
        x, y = m(lon, lat)
        if x<=10**20 and y<=10**20\
        and xrange[0]<x<xrange[1] and yrange[0]<=y<=yrange[1]:
          self.ax.text(x, y, r"$%0.0f^\circ$" % tlon)

    # write colorbar
    self.set_colorbar(p[-1], [cmin, cmax], logz, base)

    # add legend if there are any non-trivial labels
    self.add_legend_if_labels_exist(loc=loc, scatterpoints=1)

# =============================================================================
# Wrappers to translate ligolw table into plot
# =============================================================================

# =============================================================================
# Plot time series of data set

def plot_data_series(data, outfile, x_format='time', zero=None, \
                     zeroindicator=False, **kwargs):

  """
    Plot the time series / spectrum of a given set (or given sets) of data.

    Arguments:

      data : list
        list of (ChannelName,x_data,y_data) tuples with channel name (or data 
        source) and time/freq, amplitude arrays for each channel. Channels are
        plotted in the order given.
      outfile : str
        output plot path

    Keyword Arguments:

      x_format : [ 'time' | 'frequency' ]
        type of data for x_axis, allows formatting of axes
      zero : [ float | int | LIGOTimeGPS ]
        time around which to centre time series plot
      zeroindicator : [ False | True ]
        indicate zero time with veritcal dashed line, default: False

    Unnamed keyword arguments:

      logx : [ True | False ]
        boolean option to display x-axis in log scale.
      logy : [ True | False ]
        boolean option to display y-axis in log scale.
      xlim : tuple
        (xmin, xmax) limits for x-axis
      ylim : tuple
        (ymin, ymax) limits for y-axis
      xlabel : string
        label for x-axis
      ylabel : string
        label for y-axis
      title : string
        title for plot
      subtitle : string
        subtitle for plot

    All other given arguments will be passed to matplotlib.axes.Axes.plot.
  """

  # format times
  if not kwargs.has_key('xlim') or not kwargs['xlim']:
    start = min([min(d[1]) for d in data])
    end   = max([max(d[1]) for d in data])
    kwargs['xlim'] = [start,end]

  if not zero:
    zero = kwargs['xlim'][0]

  # set plot time unit whether it's used or not
  if x_format=='time':
    unit, timestr = time_unit(kwargs['xlim'][1]-kwargs['xlim'][0])

  # set labels
  if x_format=='time':
    zero = LIGOTimeGPS('%.3f' % zero)
    if zero.nanoseconds==0:
      tlabel = datetime(*date.XLALGPSToUTC(LIGOTimeGPS(zero))[:6])\
                   .strftime("%B %d %Y, %H:%M:%S %ZUTC")
    else:
      tlabel = datetime(*date.XLALGPSToUTC(LIGOTimeGPS(zero.seconds))[:6])\
                    .strftime("%B %d %Y, %H:%M:%S %ZUTC")
      tlabel = tlabel.replace(' UTC', '.%.3s UTC' % zero.nanoseconds)
    xlabel   = kwargs.pop('xlabel',\
                          'Time (%s) since %s (%s)' % (timestr, tlabel, zero))
    title    = kwargs.pop('title', 'Time series')
  else:
    xlabel = kwargs.pop('xlabel', 'Frequency (Hz)')
    title  = kwargs.pop('title', 'Frequency spectrum')

  ylabel   = kwargs.pop('ylabel', 'Amplitude')
  subtitle = kwargs.pop('subtitle', '')

  # customise plot appearance
  set_rcParams()

  # get limits
  xlim = kwargs.pop('xlim', None)
  ylim = kwargs.pop('ylim', None)
  calendar_time = kwargs.pop('calendar_time', False)

  # get axis scales
  logx = kwargs.pop('logx', False)
  logy = kwargs.pop('logy', False)

  # get legend loc
  loc = kwargs.pop('loc', 'best')

  # get colorbar options
  hidden_colorbar = kwargs.pop('hidden_colorbar', False)

  # get savefig option
  bbox_inches = kwargs.pop('bbox_inches', 'tight')

  # generate plot object
  plot = DataPlot(xlabel, ylabel, title, subtitle)

  # set plot params
  style = kwargs.pop('style', '-')
  if style in ['-', '--', '-.', ':']:
    kwargs.setdefault('linestyle', style)
    kwargs.setdefault('linewidth', 2)
    kwargs.pop('marker', None)
  else:
    kwargs.setdefault('marker', style)
    kwargs.setdefault('markersize', 5)
    kwargs.setdefault('linewidth', 0)
    kwargs.pop('linestyle', ' ')

  # get uniq data sets that aren't errors
  allchannels = []
  channels    = []
  for i,(c,_,_) in enumerate(data):
    allchannels.append(c)
    if not re.search('(min|max)\Z', str(c)):
      channels.append((i,c))

  # add data
  for i,c in channels:
    x_data,y_data = data[i][1:]
    if x_format=='time':
      if calendar_time:
        x_data = gps2datenum(numpy.array(map(float, x_data)))
      else:
        x_data = (x_data.astype(float)-float(zero))/unit
    lab = str(c)
    if lab != '_': lab = lab.replace('_', '\_')
    plot.add_content(x_data, y_data, label=lab,**kwargs)

  # finalize plot
  plot.finalize(logx=logx, logy=logy, loc=loc, hidden_colorbar=hidden_colorbar)

  # plot errors
  l = None
  for (i,channel),c in itertools.izip(channels, plotutils.default_colors()):
    try:
      minidx = allchannels.index('%s_min' % str(channel))
      maxidx = allchannels.index('%s_max' % str(channel))
    except ValueError:
      continue
    y = []
    for idx in [minidx,maxidx]:
      x_data,y_data = data[idx][1:]
      y.append(y_data)
      if x_format=='time':
        x_data = (numpy.array(map(float, x_data))-float(zero))/unit
      l = l!=None and l or float(kwargs.pop('linewidth', 1))/4
      plot.ax.plot(x_data, y_data, color=c, linewidth=l, **kwargs)
    plot.ax.fill_between(x_data, y[1], y[0], color=c, alpha=0.25)

  # set axes
  plot.ax.autoscale_view(tight=True, scalex=True, scaley=True)
  if ylim:
    plot.ax.set_ylim(tuple(ylim))

  # FIXME add zero indicator
  if x_format=='time':
    axis_lims = plot.ax.get_ylim()
    if zeroindicator:
      plot.ax.plot([0, 0], [axis_lims[0], axis_lims[1]], 'r--', linewidth=2)
      plot.ax.set_ylim([ axis_lims[0], axis_lims[1] ])

    # set x axis
    if xlim and calendar_time:
      plot.ax.set_xlim([gps2datenum(float(xlim[0])),\
                        gps2datenum(float(xlim[1]))])
    elif xlim:
      plot.ax.set_xlim([ float(xlim[0]-zero)/unit, float(xlim[1]-zero)/unit ])
  else:
      # set global axis limits
    if xlim:
      plot.ax.set_xlim(xlim)
    if ylim:
      plot.ax.set_ylim(ylim)


  set_ticks(plot.ax, calendar_time=calendar_time)

  plot.savefig(outfile, bbox_inches=bbox_inches,\
               bbox_extra_artists=plot.ax.texts)
  pylab.close(plot.fig)

# =============================================================================
# Plot a histogram of any column

def plot_trigger_hist(triggers, outfile, column='snr', num_bins=1000,\
                      color_column=None, color_bins=10,\
                      seglist=None, flag='unknown', start=None, end=None,\
                      livetime=None, etg=None, **kwargs):

  """
    Wrapper for dqPlotUtils.LineHistogram to plot a histogram of the value in
    any column of the ligolw table triggers. If a glue.segments.segmentlist
    seglist is given, the histogram is presented before and after removal of
    triggers falling inside any segment in the list.

    Arguments:

      triggers : glue.ligolw.table.Table
        ligolw table containing triggers
      outfile : string
        string path for output plot

    Keyword arguments:

      column : string
        valid column of triggers table to plot as histrogram
      num_bins : int
        number of histogram bins to use
      seglist : glue.segments.segmentlist
        list of segments with which to veto triggers
      flag : string
        display name of segmentlist, normally the name of the DQ flag
      start : [ float | int | LIGOTimeGPS]
        GPS start time (exclude triggers and segments before this time)
      end : [ float | int | LIGOTimeGPS]
        GPS end time (exclude triggers and segments after this time)
      livetime : [ float | int | LIGOTimeGPS ]
        span of time from which triggers and segments are valid, used to
        display histogram counts in terms of rate (Hz) for easy comparisons
      etg : string
        display name of trigger generator, defaults based on triggers tableName

    Unnamed keyword arguments:

      cumulative : [ True | False ]
        plot cumulative histogram
      rate : [ True | False ]
        plot rate histogram (normalises with given or calculated livetime)
      fill : [ True | False ]
        fill below the histogram curves, default colors:
            red (vetoed), green (not vetoed).
      logx : [ True | False ]
        boolean option to display x-axis in log scale.
      logy : [ True | False ]
        boolean option to display y-axis in log scale.
      xlim : tuple
        (xmin, xmax) limits for x-axis
      ylim : tuple
        (ymin, ymax) limits for y-axis
      xlabel : string
        label for x-axis
      ylabel : string
        label for y-axis
      title : string
        title for plot
      subtitle : string
        subtitle for plot
      greyscale : [ True | False ]
        use (non-greyscale) colour scheme suitable for greyscale plots

    All other given arguments will be passed to matplotlib.axes.Axes.plot and
    matplotlib.axes.Axes.fill_between. 
  """

  get_time = def_get_time(triggers.tableName)

  # calculate livetime
  if not start or not end:
    times = [get_time(t) for t in triggers]
  if not start:
    start = min(times)
  if not end:
    end   = max(times)
  if not livetime:
    livetime = end-start
  livetime = float(livetime)

  # format seglist
  if seglist==None:
    seglist = segments.segmentlist()
  else:
    seglist = segments.segmentlist(seglist)

  # format colours
  logz = kwargs.pop('logz', False)
  clim = kwargs.pop('clim', None)
 
  if color_column:
    colData = get_column(triggers, color_column).astype(float)
    if not clim:
      if len(colData)>0:
        clim = [colData.min(), colData.max()]
      else:
        clim = [1,10]
    if isinstance(color_bins, int):
      num_color_bins = color_bins
      if logz:
        color_bins = numpy.logspace(numpy.math.log10(clim[0]),\
                                    numpy.math.log10(clim[1]),\
                                    num=num_color_bins)
      else:
        color_bins = numpy.linspace(0, clim[1],\
                                    num=num_color_bins)
    else:
      num_color_bins = len(color_bins)
  else:
    num_color_bins = 1

  # get data
  tdata    = get_column(triggers, 'time')
  preData = []
  postData = []
  for i in range(num_color_bins):
    if color_column:
      preData.append(get_column(triggers, column).astype(float)\
                                                 [colData>=color_bins[i]])
    else:
      preData.append(get_column(triggers, column).astype(float))
    postData.append([p for j,p in enumerate(preData[i])\
                     if tdata[j] not in seglist])

  # get veto livetime
  vetoLivetime = livetime-float(abs(seglist))

  # set some random plot parameters
  greyscale = kwargs.pop('greyscale', False)
  if seglist:
    color = ['r','g']
    label = ['Before vetoes', 'After vetoes']
  else:
    color = ['b']
    label = ['_']
  if greyscale:
    color = ['k','k']
    linestyle = ['-','--']
  else:
    linestyle = ['-','-']

  # fix names for latex
  if flag:
    flag = flag.replace('_','\_')
  else:
    flag = 'Unknown'

  # customise plot appearance
  set_rcParams()

  # get limits
  xlim = kwargs.pop('xlim', None)
  ylim = kwargs.pop('ylim', None)

  # get axis scales
  logx = kwargs.pop('logx', False)
  logy = kwargs.pop('logy', False)

  # get colorbar options
  hidden_colorbar = kwargs.pop('hidden_colorbar', False)

  # get savefig option
  bbox_inches = kwargs.pop('bbox_inches', 'tight')

  # get fill
  fill = kwargs.pop('fill', False)

  # get extras
  cumulative = kwargs.pop('cumulative', False)
  rate       = kwargs.pop('rate', False)

  # set labels
  xlabel = kwargs.pop('xlabel', display_name(column))
  if rate and cumulative:
    ylabel = kwargs.pop('ylabel', 'Cumulative rate (Hz)')
  elif rate:
    ylabel = kwargs.pop('ylabel', 'Rate (Hz)')
  elif not rate and cumulative:
    ylabel = kwargs.pop('ylabel', 'Cumulative number')
  elif not rate and not cumulative:
    ylabel = kwargs.pop('ylabel', 'Number')

  # get ETG
  if not etg:
    if re.search('burst', triggers.tableName.lower()):
      etg = 'Burst'
    elif re.search('inspiral', triggers.tableName.lower()):
      etg = 'Inspiral'
    elif re.search('ringdown', triggers.tableName.lower()):
      etg = 'Ringdown'
    else:
      etg = 'Unknown'

  title = '%s triggers' % (etg.replace('_','\_'))
  if seglist:
    title += ' and %s segments' % (flag)
  title = kwargs.pop('title', title)
  if start and end:
    subtitle = '%d-%d' % (start, math.ceil(end))
  else:
    subtitle = ""
  subtitle = kwargs.pop('subtitle', subtitle)
  zlabel  = kwargs.pop('zlabel',\
                       color_column and display_name(color_column) or "")

  # generate plot object
  plot = LineHistogram(xlabel, ylabel, title, subtitle, zlabel)

  # add each data set
  if color_column: 
    for i in range(len(preData)):
      plot.add_content(preData[i], livetime=livetime, cval=color_bins[i],\
                       linestyle=linestyle[0], label=label[0], **kwargs)
  else:
    plot.add_content(preData[0], livetime=livetime, color=color[0],\
                     linestyle=linestyle[0], label=label[0], **kwargs)
    if seglist:
      plot.add_content(postData[0], livetime=vetoLivetime, color=color[1],\
                       linestyle=linestyle[1], label=label[1], **kwargs)

  

  # finalize plot with histograms
  if not num_bins: num_bins=100
  if color_column:
    colorbar = logz and 'log' or 'linear'
  else:
    colorbar = None
  plot.finalize(num_bins=num_bins, logx=logx, logy=logy, cumulative=cumulative,\
                rate=rate, fill=fill, colorbar=colorbar,\
                hidden_colorbar=hidden_colorbar)
  plot.ax.autoscale_view(tight=True, scalex=True, scaley=True)

  # set lower y axis limit
  if rate:
    ymin = 1/livetime
  elif logy:
    ymin = 0.5
  else:
    ymin = plot.ax.get_ylim()[0]
  plot.ax.set_ybound(lower=ymin)

  # set global axis limits
  if xlim:
    plot.ax.set_xlim(xlim)
  if ylim:
    plot.ax.set_ylim(ylim)

  # set global ticks
  set_ticks(plot.ax)

  # save figure
  plot.savefig(outfile, bbox_inches=bbox_inches,\
               bbox_extra_artists=plot.ax.texts)

# =============================================================================
# Plot one column against another column coloured by any third column

def plot_triggers(triggers, outfile, reftriggers=None, xcolumn='time',\
                  ycolumn='snr', zcolumn=None, rankcolumn=None, etg=None, start=None, end=None,\
                  zero=None, seglist=None, flag=None, **kwargs):

  """
    Plots ycolumn against xcolumn for columns in given
    Sngl{Burst,Inspiral}Table object triggers, coloured by the zcolumn
    highlighting those entries falling inside one of the entries in the
    glue.segments.segmentlist object segments, if given. 

    'time' given as a column name is a special case, since s and ns times are
    stored separately in the SnglTable structures. In this case the
    trigger.get_xxxx() function is called.

    Arguments:

      triggers : glue.ligolw.table.Table
        ligolw table containing triggers
      outfile : string
        string path for output plot

    Keyword arguments:

      xcolumn : string
        valid column of triggers table to plot on x-axis
      ycolumn : string
        valid column of triggers table to plot on y-axis
      zcolumn : string
        valid column of triggers table to use for colorbar (optional).
      rankcolumn : string
        valid column of triggers table to use for ranking events (optional).
      etg : string
        display name of trigger generator, defaults based on triggers tableName 
      start : [ float | int | LIGOTimeGPS ]
        GPS start time of plot
      end : [ float | int | LIGOTimeGPS ]
        GPS end time of plot
      zero : [ float | int | LIGOTimeGPS ]
        time around which to centre plot
      seglist : glue.segments.segmentlist
        list of segments with which to veto triggers
      flag : string
        display name of segmentlist, normally the name of the DQ flag

    Unnamed keyword arguments:

      detchar : [ True | False ]
        use 'DetChar' style for scatter plot with colorbar, triggers below given
        dcthreshold are small with no edges, whilst other triggers are normal
      dcthreshold : float
        threshold below which scatter points are small with no edges when using
        DetChar plotting style
      logx : [ True | False ]
        boolean option to display x-axis in log scale.
      logy : [ True | False ]
        boolean option to display y-axis in log scale.
      logz : [ True | False ]
        boolean option to display z-axis in log scale.
      xlim : tuple
        (xmin, xmax) limits for x-axis. Triggers outside range are removed.
      ylim : tuple
        (ymin, ymax) limits for y-axis. Triggers outside range are removed.
      zlim : tuple
        (zmin, zmax) limits for z-axis. Triggers outside range are removed.
      clim : tuple
        (cmin, cmax) limits for color scale. Triggers outside range are moved
        onto boundary.
      xlabel : string
        label for x-axis
      ylabel : string
        label for y-axis
      zlabel : string
        label for z-axis
      title : string
        title for plot
      subtitle : string
        subtitle for plot
      greyscale : [ True | False ]
        use (non-greyscale) colour scheme suitable for greyscale plots

    All other given arguments will be passed to matplotlib.axes.Axes.scatter. 
  """

  from pylal import plotutils

  # test multiple tables
  if not len(triggers)==0 and \
     (isinstance(triggers[0], tuple) or isinstance(triggers[0], list)):
    assert not zcolumn,\
           "Can only plot single table when using colorbar plot"
    tables = [t[1] for t in triggers]
    tablelabel = [t[0] for t in triggers]
    for i,t in enumerate(tablelabel):
      if t!='_':
        tablelabel[i] = t.replace('_','\_')
  else:
    tables = [triggers]
    tablelabel = '_'

  # get time column
  get_time = []
  for t in tables:
    get_time.append(def_get_time(t.tableName))

  # set start and end time if needed
  if not start or not end:
    times = [get_time[i](t)  for i in xrange(len(tables)) for t in tables[i]]
  if not start and len(times)>=1:
    start = int(math.floor(min(times)))
  elif not start:
    start = 0
  if not end and len(times)>=1:
    end   = int(math.ceil(max(times)))
  elif not end:
    end   = 1

  # set zero time
  if not zero:
    zero = start

  # get time params
  unit,timestr = time_unit(end-start)

  # set up segmentlist
  if seglist:
    segs = segments.segmentlist(seglist)
  if not seglist:
    segs = segments.segmentlist()

  # get axis limits
  xlim = kwargs.pop('xlim', None)
  ylim = kwargs.pop('ylim', None)
  zlim = kwargs.pop('zlim', None)
  clim = kwargs.pop('clim', zlim)
  calendar_time = kwargs.pop('calendar_time', None)

  # set up columns
  columns = list(map(str.lower, [xcolumn, ycolumn]))
  if zcolumn: 
    columns.append(zcolumn.lower())
    if rankcolumn: 
      columns.append(rankcolumn.lower())
    else: columns.append(zcolumn.lower())

  # set up limits
  limits    = [xlim, ylim, zlim, None]
  for i,col in enumerate(columns):
    if re.search('time\Z', col) and not limits[i]:
      limits[i] = [start,end]

  # get veto info
  if seglist:
    tdata = get_column(triggers, 'time')

  # get all data
  vetoData  = []
  nvetoData = []
  for j,tab in enumerate(tables):
    vetoData.append({})
    nvetoData.append({})
    # get veto info
    if seglist:
      tdata = get_column(tab, 'time')
    # get data
    for i,col in enumerate(columns):
      nvetoData[j][col]  = get_column(tab, col).astype(float)
    # apply limits and vetoes
    condition = True
    for i,col in enumerate(columns):
      if limits[i]:
        condition = condition & (limits[i][0] <= nvetoData[j][col])\
                              & (nvetoData[j][col] <= limits[i][1])
    for col in nvetoData[j].keys():
      nvetoData[j][col] = nvetoData[j][col][condition]
      if seglist:
        vetoData[j][col] = numpy.asarray([d for i,d in\
                                          enumerate(nvetoData[j][col]) if\
                                          tdata[i] in seglist])
      else:
        vetoData[j][col] = numpy.array([])

  data = {}

      
    
  # normalize zcolumn by time-averaged value
  whitenedFlag = kwargs.pop('whitened', False)
  if zcolumn and whitenedFlag:
    # get ref data if provided
    refData = {}
    if reftriggers:
      for i,col in enumerate(columns):
        refData[col]  = get_column(reftriggers, col).astype(float)
    else:
      for i,col in enumerate(columns):
        refData[col]  = nvetoData[0][col]
      
    uniqYvalues = numpy.unique1d(refData[ycolumn])
    # building look back table by hand, is included in unique1d for numpy >= v1.3
    for yVal in uniqYvalues:
      backTable = numpy.where(yVal == nvetoData[0][ycolumn])
      zMedian =  numpy.median(refData[zcolumn][yVal == refData[ycolumn]])
      for  iTrig in backTable[0]:
        nvetoData[0][zcolumn][iTrig] /= zMedian

  # filter zcolumn by  provided poles/zeros filter as a function of ycolumn
  flatenedFlag = kwargs.pop('filter', False)
  if zcolumn and flatenedFlag:
    # get filter params
    polesList = kwargs.pop('poles', None)
    zerosList = kwargs.pop('zeros', None)
    amplitude = kwargs.pop('amplitude', 1)
    nvetoData[0][zcolumn] *= amplitude
    for filtPole in polesList:
      nvetoData[0][zcolumn] /= abs(nvetoData[0][ycolumn] - filtPole)
    for filtZero in zerosList:
      nvetoData[0][zcolumn] *= abs(nvetoData[0][ycolumn] - filtZero)
    nvetoData[0][zcolumn].astype(float)

  # flaten zcolumn by 1/sqrt of sum given rational fraction mononomes as a 
  # function of ycolumn
  flatenedFlag = kwargs.pop('flaten', False)
  if zcolumn and flatenedFlag:
    # get filter params
    expList = kwargs.pop('exponents', None)
    constList = kwargs.pop('constants', None)
    filter = numpy.zeros(len(nvetoData[0][zcolumn]))
    for iTerm, exponent in enumerate(expList):
      filter += pow(constList[iTerm]*numpy.power(nvetoData[0][ycolumn],expList[iTerm]),2)
    filter = numpy.sqrt(filter)
    nvetoData[0][zcolumn] /= filter
  
  # median/min/max of ycolumn binned by exact xcolumn values
  minmaxmedianFlag = kwargs.pop('minmaxmedian', False)
  if minmaxmedianFlag:
    uniqXvalues = numpy.unique1d(nvetoData[j][xcolumn])
    # building look back table by hand, is included in unique1d for numpy >= v1.3
    for xVal in uniqXvalues:
      backTable = numpy.where(xVal == nvetoData[j][xcolumn])
      if len(backTable[0]) > 3:
        nvetoData[j][ycolumn][backTable[0][0]] =\
            numpy.median(nvetoData[j][ycolumn][xVal == nvetoData[j][xcolumn]])
        nvetoData[j][ycolumn][backTable[0][1]] =\
            numpy.min(nvetoData[j][ycolumn][xVal == nvetoData[j][xcolumn]])
        nvetoData[j][ycolumn][0][backTable[0][2]] =\
            numpy.max(nvetoData[j][ycolumn][xVal == nvetoData[j][xcolumn]])
        for iTrig in backTable[0][3:]:
          nvetoData[j][ycolumn][iTrig] = numpy.nan

  # down-sample (xaxis) the triggers by plotting only median z-value over averaging window 
  medianDuration = float(kwargs.pop('medianduration', 0))
  stdDuration = float(kwargs.pop('stdduration', 0))
  if medianDuration or stdDuration:
    uniqYvalues = numpy.unique1d(nvetoData[0][ycolumn])
    if medianDuration:
      avDuration = medianDuration
    elif stdDuration:
      avDuration = stdDuration
    tedges = numpy.arange(float(start), float(end), avDuration)
    tedges = numpy.append(tedges, float(end))
    # array for repacking triggers into time-frequency bins
    repackTrig = {}
    for yVal in uniqYvalues:
      repackTrig[yVal] = {}
      for iBin in range(len(tedges)-1):
        repackTrig[yVal][iBin] = []
    # building new set of triggers
    newTrigs = {}
    newTrigs[xcolumn] = []
    newTrigs[ycolumn] = []
    newTrigs[zcolumn] = []
    for iTrig in range(len(nvetoData[0][zcolumn])):
      trigVal = nvetoData[0][ycolumn][iTrig]
      trigBin = math.floor((nvetoData[0][xcolumn][iTrig]-float(start))/avDuration)
      repackTrig[trigVal][trigBin].append(nvetoData[0][zcolumn][iTrig])
    for yVal in uniqYvalues:
      for iBin in range(len(tedges)-1):
        # keep only if at least a few elements, std doesn't make sens otherwise
        if medianDuration and len(repackTrig[yVal][iBin]) > 0:
          newTrigs[xcolumn].append( (tedges[iBin]+tedges[iBin+1])/2 )
          newTrigs[ycolumn].append(yVal)
          newTrigs[zcolumn].append(numpy.median(numpy.array(repackTrig[yVal][iBin])))
        elif stdDuration and len(repackTrig[yVal][iBin]) > 5:
          newTrigs[xcolumn].append( (tedges[iBin]+tedges[iBin+1])/2 )
          newTrigs[ycolumn].append(yVal)
          newTrigs[zcolumn].append(numpy.std(numpy.array(repackTrig[yVal][iBin]))/ \
                                     numpy.median(numpy.array(repackTrig[yVal][iBin])))
          if newTrigs[zcolumn][-1] == 0:
            newTrigs[zcolumn][-1] = numpy.nan
    for i,col in enumerate(columns):
      nvetoData[0][col] = numpy.array(newTrigs[col])

  # get limits
  for i,col in enumerate(columns):
    if not limits[i]:
      limits[i] = [0,0]
      for j in xrange(len(tables)):
        data[col] = numpy.concatenate((nvetoData[j][col], vetoData[j][col]))
        if len(data[col])>=1:
          limits[i][0] = min(data[col].min()*0.99, limits[i][0])
          limits[i][1] = max(data[col].max()*1.01, limits[i][1])

    # renormalise time and set time axis label unless given
    if re.search('time\Z', col):
      renormalise = True
      if kwargs.has_key('xlabel'):  renormalise = False
      if renormalise:
        for j in xrange(len(tables)):
          if calendar_time:
            vetoData[j][col] = gps2datenum(vetoData[j][col])
            nvetoData[j][col] = gps2datenum(nvetoData[j][col])
          else:
            vetoData[j][col] = (vetoData[j][col]-float(zero))/unit
            nvetoData[j][col] = (nvetoData[j][col]-float(zero))/unit
        if calendar_time:        
          limits[i] = [gps2datenum(float(limits[i][0])),\
                       gps2datenum(float(limits[i][1]))]
        else:
          limits[i] = [float(limits[i][0]-zero)/unit,\
                       float(limits[i][1]-zero)/unit]

  # set labels
  label = {}
  for i,col in enumerate(['xcolumn', 'ycolumn', 'zcolumn']):
    if i >= len(columns):  continue
    if re.search('time\Z', columns[i]):
      zerostr = datetime(*date.XLALGPSToUTC(LIGOTimeGPS(zero))[:6])\
                         .strftime("%B %d %Y, %H:%M:%S %ZUTC")
      lab = 'Time (%s) since %s (%s)' % (timestr, zerostr, zero)
      label[columns[i]] = kwargs.pop('%slabel' % col[0], lab)
    else:
      label[columns[i]] = kwargs.pop('%slabel' % col[0],\
                                     display_name(columns[i]))

  # find loudest event
  loudest = {}
  if len(columns)==4 and\
     len(nvetoData[0][columns[0]])+len(vetoData[0][columns[0]])>=1:
    # find loudest vetoed event
    vetomax = 0
    if len(vetoData[0][columns[3]])>=1:
      vetomax = vetoData[0][columns[3]].max()
    nvetomax = 0
    # find loudest unvetoed event
    if len(nvetoData[0][columns[3]])>=1:
      nvetomax = nvetoData[0][columns[3]].max()
    if vetomax == nvetomax == 0:
      pass
    # depending on which one is loudest, find loudest overall event
    elif vetomax > nvetomax:
      index = vetoData[0][columns[3]].argmax()
      for col in columns:
        loudest[col] = vetoData[0][col][index]
    else:
      index = nvetoData[0][columns[3]].argmax()
      for col in columns:
        loudest[col] = nvetoData[0][col][index]

  # fix flag for use with latex
  if flag:
    flag = flag.replace('_', '\_')
  else:
    flag = 'Unknown'

  # get ETG
  if not etg:
    if re.search('burst', tables[0].tableName.lower()):
      etg = 'Burst'
    elif re.search('inspiral', tables[0].tableName.lower()):
      etg = 'Inspiral'
    elif re.search('ringdown', tables[0].tableName.lower()):
      etg = 'Ringdown'
    else:
      etg = 'Unknown'
  else:
    etg = etg.replace('_', '\_')

  # customise plot appearance
  set_rcParams()

  # set title
  title = '%s triggers' % etg
  if seglist:
    title += ' and %s segments' % (flag)
  title = kwargs.pop('title', title)

  if len(columns)==4 and loudest:
    subtitle = "Loudest event by %s:" % display_name(columns[-1])
    for col in columns:
      maxcol = loudest[col]
      if re.search('time\Z', col) and renormalise:
        maxcol = maxcol*unit+zero
      loudstr = "%s=%.2f" % (display_name(col), maxcol)
      if not re.search(loudstr, subtitle):
        subtitle += ' %s' % loudstr
  elif start and end:
    subtitle = '%s-%s' % (start, end)
  else:
    subtitle = ''
  subtitle = kwargs.pop('subtitle', subtitle)

  # get axis scales
  logx = kwargs.pop('logx', False)
  logy = kwargs.pop('logy', False)
  logz = kwargs.pop('logz', False)

  # get colorbar options
  hidden_colorbar = kwargs.pop('hidden_colorbar', False)

  # get savefig option
  bbox_inches = kwargs.pop('bbox_inches', 'tight')

  # get detchar plot params
  detchar = kwargs.pop('detchar', False)
  dcthresh = float(kwargs.pop('dcthreshold', 10))

  # get greyscale param
  greyscale = kwargs.pop('greyscale', False)
  if greyscale and not kwargs.has_key('cmap'):
    kwargs['cmap'] = pylab.matplotlib.colors.LinearSegmentedColormap('clrs',\
                                           pylab.matplotlib.cm.hot._segmentdata)

  # initialise standard scatter plot
  if len(columns)==2:
    plotutils.default_colors = lambda: itertools.cycle(('b', 'r', 'g', 'c', 'm', 'y', 'k'))
    plot = ScatterPlot(label[columns[0]], label[columns[1]], title, subtitle)
    for j in xrange(len(tables)):
      if len(nvetoData[j][columns[0]])>=1:
        plot.add_content(nvetoData[j][columns[0]], nvetoData[j][columns[1]],\
                         label=tablelabel[j], **kwargs)
      # add veto triggers
      if len(vetoData[j][columns[0]])>=1:
        plot.add_content(vetoData[j][columns[0]], vetoData[j][columns[1]],\
                         label=tablelabel[j], marker='x', color='r')
    # finalise
    plot.finalize(logx=logx, logy=logy, hidden_colorbar=hidden_colorbar)
  # initialise scatter plot with colorbar
  elif len(columns)==4:
    # initialize color bar plot
    if detchar:
      plot = DetCharScatterPlot(label[columns[0]], label[columns[1]],\
                                label[columns[2]], title, subtitle)
    else:
      plot = ColorbarScatterPlot(label[columns[0]], label[columns[1]],\
                                 label[columns[2]], title, subtitle)

    # add non veto triggers
    if len(nvetoData[0][columns[0]])>=1:
      plot.add_content(nvetoData[0][columns[0]], nvetoData[0][columns[1]],\
                       nvetoData[0][columns[2]], nvetoData[0][columns[3]],\
                       **kwargs)
    # add veto triggers
    if len(vetoData[0][columns[0]])>=1:
      plot.add_content(vetoData[0][columns[0]], vetoData[0][columns[1]],\
                       vetoData[0][columns[2]], vetoData[0][columns[2]],\
                       marker='x', edgecolor='r',\
                       **kwargs)
    # finalise
    if detchar:
      plot.finalize(logx=logx, logy=logy, logz=logz, clim=clim,\
                    rankthreshold=dcthresh)
    else:
      plot.finalize(logx=logx, logy=logy, logz=logz, clim=clim)
    # add loudest event to plot
    if loudest:
      plot.ax.plot([loudest[columns[0]]], [loudest[columns[1]]],\
                   marker='*', color='gold', markersize=15)

  # set axes
  plot.ax.autoscale_view(tight=True, scalex=True, scaley=True)

  if limits[0]:
    plot.ax.set_xlim(limits[0])
  if limits[1]:
    plot.ax.set_ylim(limits[1])

  # reset ticks
  set_ticks(plot.ax, calendar_time=calendar_time)

  # get both major and minor grid lines
  plot.savefig(outfile, bbox_inches=bbox_inches,\
               bbox_extra_artists=plot.ax.texts)
  pylab.close(plot.fig)

# =============================================================================
# Plot a histogram of segment duration

def plot_segment_hist(segs, outfile, keys=None, num_bins=100, coltype=int,\
                      **kwargs):

  """
    segments.
    Plots a histogram of segment duration for the glue.segments.segmentlist

    Arguments:

      segs : [ glue.segments.segmentlist | glue.segments.segmentlistdict ]
        list of segments with which to veto triggers, use dict for multiple
        datasets
      outfile : string 
        string path for output plot
   
    Keyword arguments:

      flag : string
        display name for segments, normally the name of the DQ flag
      logx : [ True | False ]
        boolean option to display x-axis in log scale.
      logy : [ True | False ]
        boolean option to display y-axis in log scale.
  """

  # customise plot appearance
  set_rcParams()

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

  # get colorbar option
  hidden_colorbar = kwargs.pop('hidden_colorbar', False)

  # format mutltiple segments
  if isinstance(segs,list):
    if keys and isinstance(keys, 'str'):
      segs = segments.segmentlistdict({keys:segs})
      keys = [keys]
    elif keys:
      segs = segments.segmentlistdict({keys[0]:segs})
    else:
      segs = segments.segmentlistdict({'_':segs})
      keys = ['_']
  else:
    if not keys:
      keys = segs.keys()
    for i,flag in enumerate(keys):
      keys[i] = display_name(flag)
      if keys[i]!=flag:
        segs[keys[i]] = segs[flag]
        del segs[flag]

  if kwargs.has_key('color'):
    colors = kwargs.pop('color').split(',')
  else:
    colors = zip(*zip(range(len(keys)), plotutils.default_colors()))[-1]

  # generate plot object
  plot = VerticalBarHistogram(xlabel, ylabel, title, subtitle)

  # add each segmentlist
  for flag,c in zip(keys, colors):
    plot.add_content([float(abs(seg)) for seg in segs[flag]],\
                      label=flag, color=c, **kwargs)

  # finalize plot with histograms
  plot.finalize(num_bins=num_bins, logx=logx, logy=logy,\
                hidden_colorbar=hidden_colorbar)

  # set limits
  plot.ax.autoscale_view(tight=True, scalex=True, scaley=True)
  if ylim:
    ylim = map(float, ylim)
    plot.ax.set_ylim(ylim)
  if xlim:
    xlim = map(float, xlim)
    plot.ax.set_xlim(xlim)

  # save figure
  plot.savefig(outfile, bbox_inches=bbox_inches,\
               bbox_extra_artists=plot.ax.texts)
  pylab.close(plot.fig)

# =============================================================================
# Plot rate versus time in bins

def plot_trigger_rate(triggers, outfile, average=600, start=None, end=None,\
                      zero=None, bincolumn='peak_frequency', bins=[],\
                      etg='Unknown', **kwargs):

  """
    Plot rate versus time for the given ligolw table triggers, binned by the
    given bincolumn using the bins list.

    Arguments:

      triggers : glue.ligolw.table
        LIGOLW table containing a list of triggers
      outfile : string
        string path for output plot

    Keyword arguments:

      average : float
        Length (seconds) of rate segment
      start : [ float | int | LIGOTimeGPS ]
        GPS start time
      end : [ float | int | LIGOTimeGPS ]
        GPS end time
      zero : [ float | int | LIGOTimeGPS ]
        GPS time to use for 0 on time axis
      bincolumn : string
        valid column of the trigger table to use for binning
      bins : list
        list of tuples defining the rate bins
      etg : string
        display name of trigger generator
      logy : [ True | False ]
        boolean option to display y-axis in log scale
      ylim : tuple
        (ymin, ymax) limits for rate axis
  """

  tableName = triggers.tableName.lower()
  get_time = def_get_time(tableName)

  # set start and end times
  if not start and not end:
    times = [get_time(t) for t in triggers]
  if not start:
    start = min(times)
  if not end:
    end   = max(times)

  if not zero:
    zero = start

  # set plot time unit whether it's used or not
  unit, timestr = time_unit(end-start)

  # set ETG
  if not etg:
    if re.search('burst', tableName):
      etg = 'Burst'
    elif re.search('inspiral', tableName):
      etg = 'Inspiral'
    elif re.search('ringdown', tableName):
      etg = 'Ringdown'
    else:
      etg = 'Unknown'

  # get limits
  calendar_time = kwargs.pop('calendar_time', False)
  if calendar_time:
    xlim = kwargs.pop('xlim', [gps2datenum(float(start)),\
                               gps2datenum(float(end))])
  else:
    xlim = kwargs.pop('xlim', [float(start-zero)/unit, float(end-zero)/unit])
  ylim = kwargs.pop('ylim', None)

  # get axis scales
  logx = kwargs.pop('logx', False)
  logy = kwargs.pop('logy', False)

  # get averaging time
  average = kwargs.pop('average',average)

  # get colorbar options
  hidden_colorbar = kwargs.pop('hidden_colorbar', False)

  # get savefig option
  bbox_inches = kwargs.pop('bbox_inches', 'tight')

  # format ybins
  if not bins:
    bins  = [[0,float('inf')]]
  ybins   = [map(float, bin) for bin in bins]

  # bin data
  tbins   = {}
  rate = {}
  for bin in ybins:
    if calendar_time:
      tbins[bin[0]] = list(gps2datenum(numpy.arange(float(start), float(end),\
                                                    average)))
    else:
      tbins[bin[0]] = list(numpy.arange(0,float(end-start), average)/unit)
    rate[bin[0]] = list(numpy.zeros(len(tbins[bin[0]])))

  for trig in triggers:
    x = int(float(getTrigAttribute(trig, 'time')-start)//average)
    y = getTrigAttribute(trig, bincolumn)
    for bin in ybins:
      if bin[0] <= y < bin[1]:
        if x>=len(rate[bin[0]]):
          print "trigger after end time, something is wrong", x, len(rate[bin[0]])
          continue
        rate[bin[0]][x] += 1/average
        break

  # if logscale includes zeros, pylab.scatter will break, so remove zeros
  if logy:
    for bin in ybins:
      removes = 0
      numtbins = len(tbins[bin[0]])
      for rbin in xrange(0,numtbins):
        if rate[bin[0]][rbin-removes]==0:
          rate[bin[0]].pop(rbin-removes)
          tbins[bin[0]].pop(rbin-removes)
          removes+=1

  # set labels
  etg   = etg.replace('_', '\_')
  zero = LIGOTimeGPS('%.3f' % zero)
  if zero.nanoseconds==0:
    tlabel = datetime(*date.XLALGPSToUTC(LIGOTimeGPS(zero))[:6])\
                 .strftime("%B %d %Y, %H:%M:%S %ZUTC")
  else:
    tlabel = datetime(*date.XLALGPSToUTC(LIGOTimeGPS(zero.seconds))[:6])\
                  .strftime("%B %d %Y, %H:%M:%S %ZUTC")
    tlabel = tlabel.replace(' UTC', '.%.3s UTC' % zero.nanoseconds)
  xlabel = kwargs.pop('xlabel',\
                      'Time (%s) since %s (%s)' % (timestr, tlabel, zero))
  ylabel = kwargs.pop('ylabel', 'Rate (Hz)')
  title = kwargs.pop('title', '%s triggers binned by %s'\
                              % (etg, display_name(bincolumn)))
  if start and end:
    subtitle = '%s-%s' % (start, end)
  else:
    subtitle = " "
  subtitle = kwargs.pop('subtitle', subtitle)

  # customise plot appearance
  set_rcParams()

  # generate plot object
  plot = ScatterPlot(xlabel, ylabel, title, subtitle)

  # plot rates
  for bin in ybins:
    if logy:
      if len(rate[bin[0]])>0:
        plot.add_content(tbins[bin[0]], rate[bin[0]],\
                         label='-'.join(map(str, bin)), **kwargs)
      else:
        plot.add_content([1],[0.1], label='-'.join(map(str, bin)),\
                         visible=False)
    else:
      plot.add_content(tbins[bin[0]], rate[bin[0]], label='-'.join(map(str, bin)),\
                       **kwargs)

  # finalise plot
  plot.finalize(logx=logx, logy=logy, hidden_colorbar=hidden_colorbar)

  # set limits
  plot.ax.autoscale_view(tight=True, scalex=True, scaley=True)
  plot.ax.set_xlim(xlim)
  if ylim:
    plot.ax.set_ylim(ylim)

  # normalize ticks
  set_ticks(plot.ax, calendar_time)

  # save
  plot.savefig(outfile, bbox_inches=bbox_inches,\
               bbox_extra_artists=plot.ax.texts)
  pylab.close(plot.fig)

# =============================================================================
# Plot RMS versus time in bins

def plot_trigger_rms(triggers, outfile, average=600, start=None, end=None,\
                     zero=None, rmscolumn='snr', bincolumn='peak_frequency',\
                     bins=[], etg='Unknown', **kwargs):

  """
    Plot RMS versus time for the given ligolw table triggers, binned by the
    given bincolumn using the bins list.

    Arguments:

      triggers : glue.ligolw.table
        LIGOLW table containing a list of triggers
      outfile : string
        string path for output plot

    Keyword arguments:

      average : float
        Length (seconds) of RMS segment
      start : [ float | int | LIGOTimeGPS ]
        GPS start time
      end : [ float | int | LIGOTimeGPS ]
        GPS end time
      zero : [ float | int | LIGOTimeGPS ]
        GPS time to use for 0 on time axis
      rmscolumn : string
        valid column of the trigger table to RMS over
      bincolumn : string
        valid column of the trigger table to use for binning
      bins : list
        list of tuples defining the rate bins
      etg : string
        display name of trigger generator
      logy : [ True | False ]
        boolean option to display y-axis in log scale
      ylim : tuple
        (ymin, ymax) limits for rate axis
  """

  tableName = triggers.tableName.lower()
  get_time = def_get_time(tableName)

  # set start and end times
  if not start and not end:
    times = [get_time(t) for t in triggers]
  if not start:
    start = min(times)
  if not end:
    end   = max(times)

  if not zero:
    zero = start

  # set plot time unit whether it's used or not
  unit, timestr = time_unit(end-start)

  # set ETG
  if not etg:
    if re.search('burst', tableName):
      etg = 'Burst'
    elif re.search('inspiral', tableName):
      etg = 'Inspiral'
    elif re.search('ringdown', tableName):
      etg = 'Ringdown'
    else:
      etg = 'Unknown'

  # get limits
  calendar_time = kwargs.pop('calendar_time', False)
  if calendar_time:
    xlim = kwargs.pop('xlim', [gps2datenum(float(start)),\
                               gps2datenum(float(end))])
  else:
    xlim = kwargs.pop('xlim', [float(start-zero)/unit, float(end-zero)/unit])
  ylim = kwargs.pop('ylim', None)

  # get axis scales
  logx = kwargs.pop('logx', False)
  logy = kwargs.pop('logy', False)

  # get averaging time
  average = kwargs.pop('average',average)

  # get colorbar options
  hidden_colorbar = kwargs.pop('hidden_colorbar', False)

  # get savefig option
  bbox_inches = kwargs.pop('bbox_inches', 'tight')

  # format ybins
  if not bins:
    bins  = [[0,float('inf')]]
  ybins   = [map(float, bin) for bin in bins]

  # bin data
  tbins   = {}
  rate = {}
  rms = {}
  for bin in ybins:
    if calendar_time:
      tbins[bin[0]] = list(gps2datenum(numpy.arange(float(start), float(end),\
                                                    average)))
    else:
      tbins[bin[0]] = list(numpy.arange(0,float(end-start), average)/unit)
    rate[bin[0]] = list(numpy.zeros(len(tbins[bin[0]])))
    rms[bin[0]] = list(numpy.zeros(len(tbins[bin[0]])))

  for trig in triggers:
    x = int(float(getTrigAttribute(trig, 'time')-start)//average)
    y = getTrigAttribute(trig, bincolumn)
    z = getTrigAttribute(trig, rmscolumn)
    for bin in ybins:
      if bin[0] <= y < bin[1]:
        rms[bin[0]][x] += z*z
        rate[bin[0]][x] += 1
        break

  # Normalize the RMS to get the mean not the sum
  for bin in ybins:
    for x in range(len(tbins[bin[0]])):
      if rate[bin[0]][x] :
        rms[bin[0]][x] = math.sqrt(rms[bin[0]][x]/rate[bin[0]][x])

  # if logscale includes zeros, pylab.scatter will break, so remove zeros
  if logy:
    for bin in ybins:
      removes = 0
      numtbins = len(tbins[bin[0]])
      for rbin in xrange(0,numtbins):
        if rms[bin[0]][rbin-removes]==0:
          rms[bin[0]].pop(rbin-removes)
          tbins[bin[0]].pop(rbin-removes)
          removes+=1

  # set labels
  etg   = etg.replace('_', '\_')
  zero = LIGOTimeGPS('%.3f' % zero)
  if zero.nanoseconds==0:
    tlabel = datetime(*date.XLALGPSToUTC(LIGOTimeGPS(zero))[:6])\
                 .strftime("%B %d %Y, %H:%M:%S %ZUTC")
  else:
    tlabel = datetime(*date.XLALGPSToUTC(LIGOTimeGPS(zero.seconds))[:6])\
                  .strftime("%B %d %Y, %H:%M:%S %ZUTC")
    tlabel = tlabel.replace(' UTC', '.%.3s UTC' % zero.nanoseconds)
  xlabel = kwargs.pop('xlabel',\
                      'Time (%s) since %s (%s)' % (timestr, tlabel, zero))
  ylabel = kwargs.pop('ylabel', 'RMS')
  title = kwargs.pop('title', '%s triggers binned by %s'\
                              % (etg, display_name(bincolumn)))
  if start and end:
    subtitle = '%s-%s' % (start, end)
  else:
    subtitle = " "
  subtitle = kwargs.pop('subtitle', subtitle)

  # customise plot appearance
  set_rcParams()

  # generms plot object
  plot = ScatterPlot(xlabel, ylabel, title, subtitle)

  # plot rmss
  for bin in ybins:
    if logy:
      if len(rms[bin[0]])>0:
        plot.add_content(tbins[bin[0]], rms[bin[0]],\
                         label='-'.join(map(str, bin)), **kwargs)
      else:
        plot.add_content([1],[0.1], label='-'.join(map(str, bin)),\
                         visible=False)
    else:
      plot.add_content(tbins[bin[0]], rms[bin[0]], label='-'.join(map(str, bin)),\
                       **kwargs)

  # finalise plot
  plot.finalize(logx=logx, logy=logy, hidden_colorbar=hidden_colorbar)

  # set limits
  plot.ax.autoscale_view(tight=True, scalex=True, scaley=True)
  plot.ax.set_xlim(xlim)
  if ylim:
    plot.ax.set_ylim(ylim)

  # normalize ticks
  set_ticks(plot.ax, calendar_time=calendar_time)

  # save
  plot.savefig(outfile, bbox_inches=bbox_inches,\
               bbox_extra_artists=plot.ax.texts)
  pylab.close(plot.fig)

# =============================================================================
# Plot segments

def plot_segments(segdict, outfile, start=None, end=None, zero=None, 
                  keys=None, highlight_segments=None, **kwargs):

  """
    Plot the segments contained within the glue.segments.segmentlistdict
    segdict to the given path string outfile. The list keys can be given to
    guarantee the order of the segments on the y-axis. x-axis limits can be
    controlled using start, end and zero. The glue.segments.segmentlist object
    highlight_segments can be given to highlight a number of segments.

    Arguments:

      segdict : glue.segments.segmentlistdict

  """

  if not start:
    minstart = lambda seglist: min(segdict[key][0][0] for key in segdict.keys()\
                                   if segdict[key])
    start = min(minstart(segdict[key]) for key in segdict.keys())
  if not end:
    maxend = lambda seglist: max(segdict[key][0][1] for key in segdict.keys()\
                                 if segdict[key])
    end    = max(maxend(segdict[key]) for key in segdict.keys())
  if not zero:
    zero = start

  # set plot time unit whether it's used or not
  unit, timestr = time_unit(end-start)

  # set labels
  zero = LIGOTimeGPS('%.3f' % zero)
  if zero.nanoseconds==0:
    tlabel = datetime(*date.XLALGPSToUTC(LIGOTimeGPS(zero))[:6])\
                 .strftime("%B %d %Y, %H:%M:%S %ZUTC")
  else:
    tlabel = datetime(*date.XLALGPSToUTC(LIGOTimeGPS(zero.seconds))[:6])\
                  .strftime("%B %d %Y, %H:%M:%S %ZUTC")
    tlabel = tlabel.replace(' UTC', '.%.3s UTC' % zero.nanoseconds)
  xlabel   = kwargs.pop('xlabel',\
                        'Time (%s) since %s (%s)' % (timestr, tlabel, zero))
  ylabel   = kwargs.pop('ylabel', "")
  title    = kwargs.pop('title', '')
  subtitle = kwargs.pop('subtitle', '')

  # get axis limits
  calendar_time = kwargs.pop('calendar_time', False)
  if calendar_time:
    xlim = kwargs.pop('xlim', [gps2datenum(float(start)),\
                               gps2datenum(float(end))])
  else:
    xlim = kwargs.pop('xlim', [float(start-zero)/unit, float(end-zero)/unit])

  # get label param
  labels_inset = kwargs.pop('labels_inset', False)

  # get savefig option
  bbox_inches = kwargs.pop('bbox_inches', 'tight')

  # get colorbar option
  hidden_colorbar = kwargs.pop('hidden_colorbar', False)

  if keys:
    # escape underscore, but don't do it twice
    keys = [key.replace('_','\_').replace('\\_','\_') for key in keys]
  segkeys = segdict.keys()
  newdict = segments.segmentlistdict()
  for key in segkeys:
    newkey = key.replace('_','\_').replace('\\\_','\_')
    newdict[newkey] = segdict[key]
  segdict = newdict

  # set params
  set_rcParams()

  plot = PlotSegmentsPlot(xlabel, ylabel, title, subtitle, t0=zero, unit=unit,\
                          calendar_time=calendar_time)
  plot.add_content(segdict, keys, **kwargs)
  plot.finalize(labels_inset=labels_inset, hidden_colorbar=hidden_colorbar)

  # indicate last frame
  if highlight_segments:
    for seg in highlight_segments:
      plot.highlight_segment(seg)

  # set x axis
  plot.ax.set_xlim(xlim)

  plot.ax.grid(True,which='major')
  plot.ax.grid(True,which='majorminor')

  set_ticks(plot.ax, calendar_time)

  plot.savefig(outfile, bbox_inches=bbox_inches,\
               bbox_extra_artists=plot.ax.texts)
  pylab.close(plot.fig)

# =============================================================================
# Helper functions
# =============================================================================

def parse_plot_config(cp, section):

  """
    Parse ConfigParser.ConfigParser section for plot parameters. Sections should
    be name '[plot xcolumn-ycolumn-zcolumn]' e.g.
    '[plot time-peak_frequency-snr]'. Returns a pair of dicts with the
    following keys:

    columns:

      xcolumn : [ string | None ]
        column string to plot on x-axis
      ycolumn : [ string | None ]
        column string to plot on y-axis
      zcolumn : [ string | None ]
        column string to plot on z-axis

    params:

      xlim : list
        [xmin, xmax] pair for x-axis limits
      ylim : list
        [ymin, ymax] pair for y-axis limits
      zlim : list
        [zmin, zmax] pair for z-axis limits
      clim : list
        [cmin, cmax] pair for colorbar limits
      logx : bool
        True / False to plot log scale on x-axis
      logy : bool
        True / False to plot log scale on y-axis
      logz : bool
        True / False to plot log scale on z-axis
  """

  columns = {'xcolumn':None, 'ycolumn':None, 'zcolumn':None, 'rankcolumn':None}
  params  = {}

  plot = re.split('[\s-]', section)[1:]
  if len(plot)>=1:
    columns['xcolumn'] = plot[0]
  if len(plot)>=2:
    columns['ycolumn'] = plot[1]
  if len(plot)>=3:
    columns['zcolumn'] = plot[2]
  if len(plot)>=4:
    columns['rankcolumn'] = plot[3]

  limits   = ['xlim', 'ylim', 'zlim', 'clim', 'exponents', 'constants',\
              'color_bins']
  filters  = ['poles', 'zeros']
  bins     = ['bins']
  booleans = ['logx', 'logy', 'logz', 'cumulative', 'rate', 'detchar',\
              'greyscale', 'zeroindicator', 'normalized', 'include_downtime',\
              'calendar_time', 'fill', 'hidden_colorbar']
  values   = ['dcthresh','amplitude','num_bins','linewidth','average','s']

  # extract plot params as a dict
  params = {}
  for key,val in cp.items(section):
    val = val.rstrip('"').strip('"')
    if key in limits:
      params[key] = map(float, val.split(','))
    elif key in filters:
      params[key] = map(complex, val.split(','))
    elif key in bins:
       params[key] = map(lambda p: map(float, p.split(',')), val.split(';'))
    elif key in booleans:
      params[key] = cp.getboolean(section, key)
    elif key in values:
      params[key] = float(val)
    else:
      params[key] = val

  return columns, params

# ==============================================================================
# Plot color map

def plot_color_map(data, outfile, data_limits=None, x_format='time',\
                   y_format='frequency', z_format='amplitude', zero=None,\
                   x_range=None, y_range=None, **kwargs):

  """
    Plots data in a 2d color map.
  """
 
  if not (isinstance(data, list) or isinstance(data, tuple)):
    data_sets = [data]
    data_limit_sets = [data_limits]
  else:
    data_sets = data
    data_limit_sets = data_limits
    if not len(data_sets) == len(data_limit_sets):
      raise AttributeError("You have given %d data sets and %d limit sets! Eh?"\
                           % (len(data_sets), len(data_limit_sets)))

  # set data limits
  for i,data_limits in enumerate(data_limit_sets):
    if not data_limits:
      numrows, numcols = numpy.shape(data_sets[i])
      if kwargs.has_key('xlim'):
        xlim = kwargs['xlim']
      else:
        xlim = [0, numrows]
      if kwargs.has_key('ylim'):
        ylim = kwargs['ylim']
      elif pylab.rcParams['image.origin'] == 'upper':
        ylim = [numrows, 0]
      else:
        ylim = [0, numrows]

      data_limit_sets[i] = [xlim[0], xlim[1], ylim[0], ylim[1]]

  # get limits
  xlim = kwargs.pop('xlim', None)
  ylim = kwargs.pop('ylim', None)
  zlim = kwargs.pop('zlim', None)
  clim = kwargs.pop('clim', None)
  calendar_time = kwargs.pop('calendar_time', False)

  # get axis scales
  logx = kwargs.pop('logx', False)
  logy = kwargs.pop('logy', False)
  logz = kwargs.pop('logz', False)

  # get savefig option
  bbox_inches = kwargs.pop('bbox_inches', 'tight')

  # restrict data to meet limits if using log
  for i,(data, data_limits) in enumerate(zip(data_sets, data_limit_sets)):
    if logx and xlim:
      shape = numpy.shape(data)
      xmin, xmax = data_limits[0], data_limits[1]
      if x_range==None:
        x_range = numpy.logspace(numpy.log10(xmin), numpy.log10(xmax),\
                               num.shape[-1])
      condition = (x_range>xlim[0]) & (x_range<=xlim[1])
      newx = numpy.where(condition)[0]
      data2 = numpy.resize(data, (len(newx), shape[-2]))
      for j in xrange(newx):
        data2[:,j] = data[newx[j],:]
      data = data2.transpose()

    if logy and ylim:
      shape = numpy.shape(data)
      ymin, ymax = data_limits[2], data_limits[3]
      if y_range==None:
        y_range = numpy.logspace(numpy.log10(ymin), numpy.log10(ymax),\
                               num=shape[-2])
      condition = (y_range>ylim[0]) & (y_range<=ylim[1])
      newy  = numpy.where((condition))[0]
      data2 = numpy.resize(data, (len(newy), shape[-1]))
      for j in xrange(shape[-1]):
        data2[:,j] = data[:,j][condition]
      data = data2
      data_limits[2:] = ylim
    data_sets[i] = data

  # get columnar params
  columns = list(map(str.lower, [x_format, y_format, z_format]))
  limits  = [xlim, ylim, zlim]
  for i,lim in enumerate(limits):
    if lim:
      limits[i] = list(map(float, lim))
  labels = ["", "", "", "", ""]

  # get zero time for normalisation
  if 'time' in columns and not zero:
    i = columns.index('time')
    if limits[i]:
      zero = limits[i][0]
    else:
      zero = min(data_limits[0] for data_limits in data_limit_sets)

  # format labels
  for i,(col,c) in enumerate(zip(columns, ['x', 'y', 'z'])):
    if col.lower() == 'time' and limits[i]:
      unit, timestr = time_unit(float(limits[i][1]-limits[i][0]))
      zero = LIGOTimeGPS('%.3f' % zero)
      if zero.nanoseconds==0:
        tlabel = datetime(*date.XLALGPSToUTC(LIGOTimeGPS(zero))[:6])\
                     .strftime("%B %d %Y, %H:%M:%S %ZUTC")
      else:
        tlabel = datetime(*date.XLALGPSToUTC(LIGOTimeGPS(zero.seconds))[:6])\
                    .strftime("%B %d %Y, %H:%M:%S %ZUTC")
        tlabel = tlabel.replace(' UTC', '.%.3s UTC' % zero.nanoseconds)
      labels[i] = kwargs.pop('%slabel' % c, 'Time (%s) since %s (%s)'\
                                            % (timestr, tlabel, zero))
      if calendar_time:
        limits[i] = gps2datenum(numpy.asarray(limits[i]))
      else:
        limits[i] = (numpy.asarray(limits[i])-float(zero))/unit
      for i,data_limits in enumerate(data_limit_sets):
        if calendar_time:
          data_limit_sets[i][0] = gps2datenum(float(data_limits[0]))
          data_limit_sets[i][1] = gps2datenum(float(data_limits[1]))
        else:
          data_limit_sets[i][0] = float(data_limits[0]-zero)/unit
          data_limit_sets[i][1] = float(data_limits[1]-zero)/unit
    else:
      labels[i] = kwargs.pop('%slabel' % c, display_name(col))

  labels[3] = kwargs.pop('title', "")
  labels[4] = kwargs.pop('subtitle', "")

  # customise plot appearance
  set_rcParams()

  # generate plot object
  plot = ColorMap(*labels)

  # add data
  for i, (data, data_limits) in enumerate(zip(data_sets, data_limit_sets)):
    plot.add_content(data, data_limits, aspect='auto')

  # finalize
  plot.finalize(logx=logx, logy=logy, logz=logz, clim=clim, origin='lower')

  if limits[0] is not None and len(limits[0])==2:
    plot.ax.set_xlim(limits[0])
  if limits[1] is not None and len(limits[1])==2:
    plot.ax.set_ylim(limits[1])

  # set global ticks
  set_ticks(plot.ax, calendar_time=calendar_time)

  # save figure
  plot.savefig(outfile, bbox_inches=bbox_inches,\
               bbox_extra_artists=plot.ax.texts)
  pylab.close(plot.fig)

# =============================================================================
# Significance drop plot (HVeto style)

def plot_significance_drop(startsig, endsig, outfile, **kwargs):

  """
    Plot significance drop for each channel relative to the application of
    HVeto round veto segments.
  """

  # get channels
  channels = startsig.keys()
  for c in channels:
    if c not in endsig.keys():
      raise AttributeError("Significance lists do not match.")
  channels.sort()

  # find winner
  wch,wsig = max(startsig.items(), key=lambda x: x[1])

  # extract parameters
  kwargs.pop('xlim', None)
  ylim     = kwargs.pop('ylim', (0, wsig+1))
  xlabel   = kwargs.pop('xlabel',   "")
  ylabel   = kwargs.pop('ylabel',   "Significance")
  title    = kwargs.pop('title',    "Coincidence significance drop plot")
  subtitle = kwargs.pop('subtitle',\
                        "Winner: %s, significance: %s" % (wch, wsig))

  kwargs.setdefault('linestyle', '-')
  kwargs.setdefault('marker', 'o')
  color    = kwargs.pop('color', None)

  # customise plot appearance
  set_rcParams()
  pylab.rcParams.update({"figure.figsize":[24,6], "xtick.labelsize": 8})

  # get colorbar options
  hidden_colorbar = kwargs.pop('hidden_colorbar', False)

  # get savefig option
  bbox_inches = kwargs.pop('bbox_inches', 'tight')

  # generate plot object
  plot = DataPlot(xlabel, ylabel, title, subtitle)

  kwargs.setdefault('markersize', 10)

  # plot each channel's drop
  for i,c in enumerate(channels):
    s   = startsig[c]
    e   = endsig[c]
    col = color and color or s>e and 'b' or 'r'
    plot.add_content([i,i], [s,e], color=col, **kwargs)

  # finalise plot object
  plot.finalize(logx=False, logy=False, hidden_colorbar=hidden_colorbar)

  # set xticks to channel names and rotate
  plot.ax.set_xlim(-1, len(channels))
  plot.ax.set_xticks(numpy.arange(0,len(channels)))
  plot.ax.set_xticklabels([c.replace('_','\_') for c in channels])
  for i,t in enumerate(plot.ax.get_xticklabels()):
    t.set_rotation(315)
    #t.set_position(((i+1)/(len(channels)+1), 0.1))
    t.set_verticalalignment('top')
    t.set_horizontalalignment('left')
    #t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='none'))

  # set ylim
  plot.ax.set_ylim(ylim)

  # turn off x grid
  plot.ax.xaxis.grid(False)

  plot.fig.subplots_adjust(bottom=0.3)

  # save figure
  plot.savefig(outfile, bbox_inches=bbox_inches,\
               bbox_extra_artists=plot.ax.texts)
  pylab.close(plot.fig)

# =============================================================================
# Plot sky positions

def plot_sky_positions(skyTable, outfile, format='radians', zcolumn=None,\
                       projection=None, centre=None, detectors=[], lines={},\
                       range=[(0,0), (1,1)], **kwargs):

  """
    Plot latitude against longitude for the given ligolw table skyTable into the
    given outfile. Uses the mpl_toolkits basemap module to plot the sky sphere 
    in a variety of projections, or simply a scatter plot if projection=None.
    Can plot lines and detector positions on top if given.

    Arguments:

      skyTable : glue.ligolw.table.Table
        ligolw table containing triggers or SkyPosition objects
      outfile : string
        string path for output plot

    Keyword arguments:

      zcolumn : string
        valid column of ligolw table to use for colorbar (optional).
      format : [ 'radians' | 'degrees' ]
        str identifying format of longtiude/ra, latitude/dec colums in table.
        Plot is always drawn in degrees.
      projection : str
        type of spherical projection to use, if any. See matplotlib Basemap
        documentation for details, recommended: 'ortho', 'hammer'.
      centre : tuple
        (longitude, latitude) pair on which to centre plot. Latitude centring
        only works for certain projections.
      detectors : list
        list of detector prefixes to plot, e.g. ['H1', 'L1'].
      lines : dict
        dict of name:table pairs from which to plot lines (+ marker) on top of
        points
      range : tuple
        ((xmin, ymin), (xmax, ymax)) tuples for lower left and upper right
        corners, in range 0-1, e.g. (0,0) is full lower left, (1,1) full upper
        right corners for full spherical projection. Only works with certain
        projections.

    Unnamed keyword arguments:

      logx : [ True | False ]
        boolean option to display x-axis in log scale. Not applicable when using
        projection.
      logy : [ True | False ]
        boolean option to display y-axis in log scale. Not applicable when using
        projection.
      logz : [ True | False ]
        boolean option to display z-axis in log scale.
      xlim : tuple
        (xmin, xmax) limits for x-axis (longitude). Triggers outside range are
        removed.
      ylim : tuple
        (ymin, ymax) limits for y-axis (latitude). Triggers outside range are
        removed.
      zlim : tuple
        (zmin, zmax) limits for z-axis. Triggers outside range are removed.
      clim : tuple
        (cmin, cmax) limits for color scale. Triggers outside range are moved
        onto boundary.
      xlabel : string
        label for x-axis
      ylabel : string
        label for y-axis
      zlabel : string
        label for z-axis
      title : string
        title for plot
      subtitle : string
        subtitle for plot

    All other given arguments will be passed to matplotlib.axes.Axes.scatter. 
  """

  format = format.lower()

  # get labels
  if projection:
    xlabel = kwargs.pop("xlabel", "")
    ylabel = kwargs.pop("ylabel", "")
  else:
    xlabel = kwargs.pop("xlabel", "Longitude (degrees of arc)")
    ylabel = kwargs.pop("ylabel", "Latitude (degrees of arc)")
  if zcolumn: tmp = zcolumn
  else:       tmp = ""
  zlabel   = kwargs.pop("zlabel", display_name(tmp))

  # get ETG
  if re.search('burst', skyTable.tableName.lower()):
    etg = 'Burst'
  elif re.search('inspiral', skyTable.tableName.lower()):
    etg = 'Inspiral'
  elif re.search('ringdown', skyTable.tableName.lower()):
    etg = 'Ringdown'
  else:
    etg = 'Unknown'
  etg = etg.replace('_', '\_')

  # get title
  title    = kwargs.pop("title", "")
  subtitle = kwargs.pop("subtitle", "")

  # get limits
  xlim = kwargs.pop('xlim', None)
  ylim = kwargs.pop('ylim', None)
  zlim = kwargs.pop('zlim', None)
  clim = kwargs.pop('clim', None)

  # get logs
  logx = kwargs.pop('logx', False)
  logy = kwargs.pop('logy', False)
  logz = kwargs.pop('logz', False)

  # get savefig option
  bbox_inches = kwargs.pop('bbox_inches', 'tight')

  # get IFO data
  detData = []
  detCol  = []
  for ifo in detectors:
    d = XLALTools.cached_detector[inject.prefix_to_name[ifo]]
    if format=='radians':
      detData.append((numpy.degrees(d.vertexLongitudeRadians),\
                      numpy.degrees(d.vertexLatitudeRadians)))
    else:
      detData.append((d.vertexLongitudeRadians, d.vertexLatitudeRadians))
    if ifo in plotsegments.PlotSegmentsPlot.color_code.keys():
      detCol.append(plotsegments.PlotSegmentsPlot.color_code[ifo])
    else:
      detCol.append('cyan')

  # get line data

  lineData = []
  lineCol  = []
  for line in lines:
    if isinstance(lines, dict):
      label = line
      line  = lines[label]
    else:
      label = '_'
    # get column
    if re.search('multi_inspiral', line.tableName):
      loncol = 'ra'
      latcol = 'dec'
    else:
      loncol = 'longitude'
      latcol = 'latitude'
    # get data
    lon = []
    lat = []
    for p in line:
      if format=='radians':
        lon.append(numpy.degrees(getattr(p, loncol)))
        lat.append(numpy.degrees(getattr(p, latcol)))
      else:
        lon.append(getattr(p, loncol))
        lat.append(getattr(p, latcol))
    lineData.append((label, [lon, lat]))
  if len(lineData):
    for c in plotutils.default_colors():
      if c=='b': continue
      lineCol.append(c)
      if len(lineCol)==len(lineData):  break

  # get point data
  columns = ['longitude', 'latitude']
  if re.search('multi_inspiral', skyTable.tableName):
    columns[0] = 'ra'
    columns[1] = 'dec'
  if zcolumn:
    columns.append(zcolumn)

  data = {}
  for col in columns:
    data[col] = []
  for i,p in enumerate(skyTable):
    # test limits
    if xlim and not xlim[0] <= getTrigAttribute(p, columns[0]) <= xlim[1]:
      continue
    if ylim and not ylim[0] <= getTrigAttribute(p, columns[1]) <= ylim[1]:
      continue
    if zlim and not zlim[0] <= getTrigAttribute(p, columns[2]) <= zlim[1]:
      continue
    # get sky data
    for col in columns[0:2]:
      if format=='radians':
        data[col].append(numpy.degrees(getTrigAttribute(p, col)))
      else:
        data[col].append(getTrigAttribute(p, col))
    # get z-axis data
    if zcolumn:
      data[columns[2]].append(getTrigAttribute(p, columns[2]))

  # customise plot appearance
  set_rcParams()
  if projection not in [None, 'hammer', 'moll', 'robin', 'sinu', 'cyl', 'merc',\
                        'mill', 'gall']:
    pylab.rcParams.update({"figure.figsize":[8,6]})

  # generate plot
  if projection:
    # colorbar plot
    if zcolumn:
      plot = ColorbarSkyPositionsPlot(xlabel, ylabel, zlabel,\
                                      title, subtitle)
      plot.add_content(data[columns[0]], data[columns[1]], data[zcolumn],\
                       **kwargs)
      for p,ifo,c in zip(detData, detectors, detCol):
        lon, lat = p
        plot.add_content([lon], [lat], [1], label=ifo, facecolor=c,\
                         edgecolor=c, **kwargs)
      for line,c in zip(lineData, lineCol):
        label, data = line
        plot.add_content(data[0], data[1], label=label, marker='+',\
                         edgecolor=c, s=4)
      plot.finalize(projection=projection, centre=centre, range=range,\
                    logz=logz, clim=clim)
    # normal plot
    else:
      plot = SkyPositionsPlot(xlabel, ylabel, title, subtitle)
      plot.add_content(data[columns[0]], data[columns[1]], label='_', **kwargs)
      for p,ifo,c in zip(detData, detectors, detCol):
        lon, lat = p
        plot.add_content([lon], [lat], label=ifo, facecolor=c, edgecolor=c,\
                         **kwargs)
      for line,c in zip(lineData, lineCol):
        label, data = line
        plot.add_content(data[0], data[1], label=label, marker='+',\
                         edgecolor=c, s=4)
      plot.finalize(projection=projection, centre=centre, range=range)

  # generate square plot
  else:
    if zcolumn:
      plot = ColorbarScatterPlot(xlabel, ylabel, zlabel, title,\
                                             subtitle)
      plot.add_content(data[columns[0]], data[columns[1]], data[columns[2]],\
                       **kwargs)
      plot.finalize(logx=logx, logy=logy, logz=logz, clim=clim)
    else:
      plot = ScatterPlot(xlabel, ylabel, title, subtitle)
      plot.add_content(data[columns[0]], data[columns[1]], **kwargs)
      plot.finalize(logx=logx, logy=logy)

    # set axes
    plot.ax.autoscale_view(tight=True, scalex=True, scaley=True)
    plot.ax.set_xlim(xlim)
    plot.ax.set_ylim(ylim)

  # save
  plot.savefig(outfile, bbox_inches=bbox_inches,\
               bbox_extra_artists=plot.ax.texts)
  pylab.close(plot.fig)
