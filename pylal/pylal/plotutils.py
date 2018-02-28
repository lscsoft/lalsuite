# Copyright (C) 2008  Nickolas Fotopoulos
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
"""
This module is intended to store generic, reusable, sub-classable plot classes
to minimize formulaic copying and pasting.
"""

from __future__ import division

__author__ = "Nickolas Fotopoulos <nvf@gravity.phys.uwm.edu>"

import itertools
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid import make_axes_locatable
from mpl_toolkits.basemap import Basemap

import numpy
import pylab
import re
import copy
import ConfigParser

from glue import iterutils
from glue import segments

from pylal import viz

# general defaults
pylab.rc("lines", markersize=12)
pylab.rc("text", usetex=True)

# Utility function
def float_to_latex(x, format="%.2g"):
    """
    Convert a floating point number to a latex representation.  In particular,
    scientific notation is handled gracefully: e -> 10^
    """
    base_str = format % x
    if "e" not in base_str:
        return base_str
    mantissa, exponent = base_str.split("e")
    exponent = exponent.lstrip("0+")
    if mantissa == "1":
        return r"10^{%s}" % exponent
    else:
        return r"%s\times 10^{%s}" % (mantissa, exponent)

##############################################################################
# abstract classes

class BasicPlot(object):
    """
    A very default meta-class to almost any plot you might want to make.
    It provides basic initialization, a savefig method, and a close method.
    It is up to developers to subclass BasicPlot and fill in the add_content()
    and finalize() methods.
    """
    def __init__(self, xlabel="", ylabel="", title="", subtitle="", **kwargs):
        """
        Basic plot initialization.  A subclass can override __init__ and call
        this one (plotutils.BasicPlot.__init__(self, *args, **kwargs)) and
        then initialize variables to hold data to plot and labels.
        """
        self.fig = pylab.figure(**kwargs)
        self.ax = self.fig.add_subplot(111)

        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        if subtitle:
            self.ax.set_title(title, x=0.5, y=1.03)
            self.ax.text(0.5, 1.035, subtitle, horizontalalignment='center',
                         transform=self.ax.transAxes, verticalalignment='top')
        else:
            self.ax.set_title(title)
        self.ax.grid(True)

    def add_content(self, data, label="_nolabel_"):
        """
        Stub.  Replace with a method that appends values or lists of values
        to self.data_sets and appends labels to self.data_labels.  Feel free
        to accept complicated inputs, but try to store only the raw numbers
        that will enter the plot.
        """
        raise NotImplementedError

    def finalize(self):
        """
        Stub.  Replace with a function that creates and makes your plot
        pretty.  Do not do I/O here.
        """
        raise NotImplementedError

    def savefig(self, *args, **kwargs):
        self.fig.savefig(*args, **kwargs)

    def close(self):
        """
        Close the plot and release its memory.
        """
        pylab.close(self.fig)

    def add_legend_if_labels_exist(self, *args, **kwargs):
        """
        Create a legend if there are any non-trivial labels.
        """

        # extract useable parameters that don't get passed to ax.legend
        alpha      = kwargs.pop("alpha", None)
        linewidth  = kwargs.pop("linewidth", None)
        markersize = kwargs.pop("markersize", None)

        # make legend if any data set requires it
        for plot_kwargs in self.kwarg_sets:
            if "label" in plot_kwargs and \
                not plot_kwargs["label"].startswith("_"):
                # generate legend
                self.ax.legend(*args, **kwargs)
                # apply extra formatting
                leg   = self.ax.get_legend()
                frame = leg.get_frame()
                if alpha:
                    frame.set_alpha(alpha)
                if linewidth:
                    for l in leg.get_lines():
                        l.set_linewidth(linewidth)
                return

##############################################################################
# utility functions

def default_colors():
    """
    An infinite iterator of some default colors.
    """
    return itertools.cycle(('b', 'g', 'r', 'c', 'm', 'y', 'k'))

def default_symbols():
    """
    An infinite iterator of some default symbols.
    """
    return itertools.cycle(('x', '^', 'D', 'H', 'o', '1', '+'))

def determine_common_bin_limits(data_sets, default_min=0, default_max=0):
    """
    Given a some nested sequences (e.g. list of lists), determine the largest
    and smallest values over the data sets and determine a common binning.
    """
    max_stat = max(list(iterutils.flatten(data_sets)) + [-numpy.inf])
    min_stat = min(list(iterutils.flatten(data_sets)) + [numpy.inf])
    if numpy.isinf(-max_stat):
        max_stat = default_max
    if numpy.isinf(min_stat):
        min_stat = default_min
    return min_stat, max_stat

def method_callable_once(f):
    """
    Decorator to make a method complain if called more than once.
    """
    def _new(self, *args, **kwargs):
        attr = "_" + f.__name__ + "_already_called"
        if hasattr(self, attr) and getattr(self, attr):
            raise ValueError, f.__name__ + " can only be called once"
        setattr(self, attr, True)
        return f(self, *args, **kwargs)
    _new.__doc__ == f.__doc__
    _new.__name__ = f.__name__
    return _new

_dq_params = {"text.usetex": True,   "text.verticalalignment": "center",
              "lines.linewidth": 2,  "xtick.labelsize": 16,
              "ytick.labelsize": 16, "axes.titlesize": 22,
              "axes.labelsize": 16,  "axes.linewidth": 1,
              "grid.linewidth": 1,   "legend.fontsize": 16,
              "legend.loc": "best",  "figure.figsize": [12,6],
              "figure.dpi": 80,      "image.origin": 'lower',
              "axes.grid": True,     "axes.axisbelow": False}

def set_rcParams(params=_dq_params):
    """
    Update pylab plot parameters, defaulting to parameters for DQ-style trigger
    plots.
    """

    # customise plot appearance
    pylab.rcParams.update(params)

def display_name(columnName):
    """
    Format the string columnName (e.g. xml table column) into latex format for
    an axis label. Formats known acronyms, greek letters, units, subscripts, and
    some miscellaneous entries.

    Examples:

    >>> display_name('snr')
    'SNR'
    >>> display_name('bank_chisq_dof')
    'Bank $\\chi^2$ DOF'
    >>> display_name('hoft')
    '$h(t)$'

    Arguments:

      columnName : str
        string to format
    """

    # define known acronyms
    acro    = ['snr', 'ra','dof', 'id', 'ms', 'far']
    # define greek letters
    greek   = ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta',\
               'theta', 'iota', 'kappa', 'lamda', 'mu', 'nu', 'xi',\
               'pi', 'rho', 'sigma', 'tau', 'upsilon', 'phi', 'chi', 'psi',\
               'omega']
    # define known units
    unit    = {'ns':'ns', 'hz':'Hz'}
    # define known subscripts
    sub     = ['flow', 'fhigh', 'hrss', 'mtotal', 'mchirp']
    # define miscellaneous entries
    misc    = {'hoft': '$h(t)$'}

    if len(columnName)==1:
        return columnName

    # find all words, preserving acronyms in upper case
    words = []
    for w in re.split('\s', columnName):
        if w[:-1].isupper(): words.append(w)
        else:           words.extend(re.split('_', w))

    # parse words
    for i,w in enumerate(words):
        if w.startswith('\\'):
            pass
        wl = w.lower()
        # get miscellaneous definitions
        if wl in misc.keys():
            words[i] = misc[wl]
        # get acronym in lower case
        elif wl in acro:
            words[i] = w.upper()
        # get numerical unit
        elif wl in unit:
            words[i] = '(%s)' % unit[wl]
        # get character with subscript text
        elif wl in sub:
            words[i] = '%s$_{\mbox{\\small %s}}$' % (w[0], w[1:])
        # get greek word
        elif wl in greek:
            words[i] = '$\%s$' % w
        # get starting with greek word
        elif re.match('(%s)' % '|'.join(greek), w):
            if w[-1].isdigit():
                words[i] = '$\%s_{%s}$''' %tuple(re.findall(r"[a-zA-Z]+|\d+",w))
            elif wl.endswith('sq'):
                words[i] = '$\%s^2$' % w[:-2]
        # get everything else
        else:
            if w[:-1].isupper():
                words[i] = w
            else:
                words[i] = w.title()
            # escape underscore
            words[i] = re.sub('(?<!\\\\)_', '\_', words[i])

    return ' '.join(words)

def add_colorbar(ax, mappable=None, visible=True, log=False, clim=None,\
                 label=None, **kwargs):
    """
    Adds a figure colorbar to the given Axes object ax, based on the values
    found in the mappable object. If visible=True, returns the Colorbar object, 
    otherwise, no return.

    Arguments:

        ax : matplotlib.axes.AxesSubplot
            axes object beside which to draw colorbar

    Keyword arguments:

        mappable : [ matplotlib.image.Image | matplotlib.contour.ContourSet... ]
            image object from which to map colorbar values
        visible : [ True | False]
            add colorbar to figure, or simply resposition ax as if to draw one
        log : [ True | False ]
            use logarithmic scale for colorbar
        clim : tuple
            (vmin, vmax) pair for limits of colorbar
        label : str
            label string for colorbar

    All other keyword arguments will be passed to pylab.colorbar. Logarithmic
    colorbars can be created by plotting log10 of the data and setting log=True.
    """

    div = make_axes_locatable(ax)
    cax = div.new_horizontal("3%", pad="1%")
    if not visible:
        return
    else:
        div._fig.add_axes(cax)
    
    # set default tex formatting for colorbar
    if pylab.rcParams['text.usetex']:
        kwargs.setdefault('format',\
            pylab.matplotlib.ticker.FuncFormatter(lambda x,pos: "$%s$"\
                                                  % float_to_latex(x)))

    # set limits
    if not clim:
        cmin = min([c.get_array().min() for c in ax.collections+ax.images])
        cmax = max([c.get_array().max() for c in ax.collections+ax.images])
        clim = [cmin, cmax]
    if log:
        kwargs.setdefault("ticks", numpy.logspace(numpy.log10(clim[0]),\
                                                  numpy.log10(clim[1]), num=9,\
                                                  endpoint=True))
    else:
        kwargs.setdefault("ticks", numpy.linspace(clim[0], clim[1], num=9,\
                                                  endpoint=True))

    # find mappable with lowest maximum
    if not mappable:
        if len(ax.collections+ax.images) == 0:
            norm = log and pylab.matplotlib.colors.LogNorm() or none
            mappable = ax.scatter([1], [1], c=[clim[0]], vmin=clim[0],\
                                   vmax=clim[1], visible=False, norm=norm)
        else:
            minindex = numpy.asarray([c.get_array().min() for c in\
                                      ax.collections+ax.images]).argmin()
            mappable = (ax.collections+ax.images)[minindex]

    # make sure the mappable has at least one element
    if mappable.get_array() is None:
        norm = log and pylab.matplotlib.colors.LogNorm() or none
        mappable = ax.scatter([1], [1], c=[clim[0]], vmin=clim[0],\
                               vmax=clim[1], visible=False, norm=norm)
    
    # generate colorbar
    colorbar = ax.figure.colorbar(mappable, cax=cax, **kwargs)
    if clim: colorbar.set_clim(clim)
    if label: colorbar.set_label(label)
    colorbar.draw_all()

    return colorbar

def parse_plot_config(cp, section):
    """
    Parser ConfigParser.ConfigParser section for plotting parameters. Returns
    a dict that can be passed to any plotutils.plot_xxx function in **kwargs
    form. Set ycolumn to 'hist' or 'rate' to generate those types of plots.

    Arguments:

        cp : ConfigParser.ConfigParser
            INI file object from which to read
        section : str
            section name to read for options

    Basic parseable options:

        xcolumn : str
            parameter to plot on x-axis    
        ycolumn : str
            parameter to plot on y-axis    
        zcolumn : str
            parameter to plot on z-axis    
        rank-by : str
            parameter by which to rank elements
        xlim : list
            [xmin, xmax] pair for x-axis limits
        ylim : list
            [ymin, ymax] pair for y-axis limits
        zlim : list
            [zmin, zmax] pair for z-axis limits
        clim : list
            [cmin, cmax] pair for colorbar limits
        logx : [ True | False ]
            plot x-axis in log scale
        logy : [ True | False ]
            plot y-axis in log scale
        logz : [ True | False ]
            plot z-axis in log scale

    Trigger plot options:

        detchar-style : [ True | False ]
            use S6-style plotting: low snr triggers small with no edges
        detchar-style-theshold : float
            z-column threshold at below which to apply detchar-style

    Trigger rate plot options:

        bins : str
            semi-colon-separated list of comma-separated bins for rate plot

    Histogram options:

        cumulative : [ True | False ]
            plot cumulative counts in histogram
        rate : [ True | False ]
            plot histogram counts as rate
        num-bins : int
            number of bins for histogram
        fill : [ True | False ]
            plot solid colour underneath histogram curve
        color-bins : str
            semi-colon-separated list of comma-separated bins for colorbar
            histogram

    Data plot options:

        zero-indicator : [ True | False ]
            draw vertical dashed red line at t=0

    Other options:

        greyscale : [ True | False ]
            save plot in black-and-white
        bbox-inches : 'tight'
            save figure with tight bounding box around Axes
        calendar-time : [ True | False ]
            plot time axis with date and time instead of time from zero.
    """
    params = dict()

    # define option types
    pairs    = ['xlim', 'ylim', 'zlim', 'colorlim']
    pairlist = ['bins', 'color-bins']
    booleans = ['logx', 'logy', 'logz', 'cumulative', 'rate', 'detchar-style',\
                'greyscale', 'zero-indicator', 'normalized', 'fill',\
                'calendar-time', 'bar']
    floats   = ['detchar-style-threshold', 'dcthreshold']
    ints     = ['num-bins']

    # construct param dict
    for key,val in cp.items(section, raw=True):
        if val == None: continue
        # remove quotes
        val = val.rstrip('"').strip('"')
        # format key key
        hkey = re.sub('_', '-', key)
        ukey = re.sub('-', '_', key)
        # get limit pairs
        if hkey in pairs:
            params[ukey] = map(float, val.split(','))
        # get bins
        elif hkey in pairlist:
            params[ukey] = map(lambda p: map(float,p.split(',')),val.split(';'))
        # get booleans
        elif hkey in booleans:
            params[ukey] = cp.getboolean(section, key)
        # get float values
        elif hkey in floats:
            params[ukey] = float(val)
        # get float values
        elif hkey in ints:
            params[ukey] = int(val)
        # else construct strings
        else:
            params[ukey] = str(val)

    return params

def log_transform(lin_range):
    """
    Return the logarithmic ticks and labels corresponding to the
    input lin_range.
    """
    log_range = numpy.log10(lin_range)
    slope = (lin_range[1] - lin_range[0]) / (log_range[1] - log_range[0])
    inter = lin_range[0] - slope * log_range[0]
    tick_range = [tick for tick in range(int(log_range[0] - 1.0),\
                                         int(log_range[1] + 1.0))\
                  if tick >= log_range[0] and tick<=log_range[1]]
    ticks = [inter + slope * tick for tick in tick_range]
    labels = ["${10^{%d}}$" % tick for tick in tick_range]
    minorticks = []
    for i in range(len(ticks[:-1])):
        minorticks.extend(numpy.logspace(numpy.log10(ticks[i]),\
                                         numpy.log10(ticks[i+1]), num=10)[1:-1])
    return ticks, labels, minorticks

def time_axis_unit(duration):
    """
    Work out renormalisation for the time axis, makes the label more
    appropriate. Returns unit (in seconds) and string descriptor.

    Example:

    >>> time_axis_unit(100)
    (1, 'seconds')

    >>> time_axis_unit(604800)
    (86400, 'days')

    Arguments:

        duration : float
            plot duration to normalise
    """
    if (duration) < 1000:
        return 1,"seconds"
    elif (duration) < 20000:
        return 60,"minutes"
    elif (duration) >= 20000 and (duration) < 604800:
        return 3600,"hours"
    elif (duration) < 8640000:
        return 86400,"days"
    else:
        return 2592000,"months"

def set_time_ticks(ax):
    """
    Quick utility to set better formatting for ticks on a time axis.
    """
    xticks = ax.get_xticks()
    if len(xticks)>1 and xticks[1]-xticks[0]==5:
        ax.xaxis.set_major_locator(pylab.matplotlib.ticker.MultipleLocator(base=2))
    return

def set_minor_ticks(ax, x=True, y=True):
    """
    Labels first minor tick in the case that there is only a single major
    tick label visible.
    """

    def even(x, pos):
        if int(str(int(x*10**8))[0]) % 2:
            return ""
        elif pylab.rcParams["text.usetex"]:
            return "$%s$" % float_to_latex(x)
        else:
            return str(int(x))

    # xticks
    if x:
        ticks = list(ax.get_xticks())
        xlim  = ax.get_xlim()
        for i,tick in enumerate(ticks[::-1]):
            if not xlim[0] <= tick <= xlim[1]:
                ticks.pop(-1)
        if len(ticks) <= 1:
            ax.xaxis.set_minor_formatter(pylab.FuncFormatter(even))

    # yticks
    if y:
        ticks = list(ax.get_yticks())
        ylim  = ax.get_ylim()
        for i,tick in enumerate(ticks[::-1]):
            if not ylim[0] <= tick <= ylim[1]:
                ticks.pop(-1)
        if len(ticks)<=1:
            ax.yaxis.set_minor_formatter(pylab.FuncFormatter(even))

    return

##############################################################################
# generic, but usable classes

class SimplePlot(BasicPlot):
    """
    Exactly what you get by calling pylab.plot(), but with the handy extras
    of the BasicPlot class.
    """
    def __init__(self, *args, **kwargs):
        BasicPlot.__init__(self, *args, **kwargs)
        self.x_data_sets = []
        self.y_data_sets = []
        self.kwarg_sets = []

    def add_content(self, x_data, y_data, **kwargs):
        self.x_data_sets.append(x_data)
        self.y_data_sets.append(y_data)
        self.kwarg_sets.append(kwargs)

    @method_callable_once
    def finalize(self, loc=0, alpha=0.8):
        # make plot
        for x_vals, y_vals, plot_kwargs in \
            itertools.izip(self.x_data_sets, self.y_data_sets, self.kwarg_sets):
            self.ax.plot(x_vals, y_vals, **plot_kwargs)

        # add legend if there are any non-trivial labels
        self.add_legend_if_labels_exist(loc=loc, alpha=alpha)

        # decrement reference counts
        del self.x_data_sets
        del self.y_data_sets
        del self.kwarg_sets

class BarPlot(BasicPlot):
    """
    A simple vertical bar plot.  Bars are centered on the x values and have
    height equal to the y values.
    """
    def __init__(self, *args, **kwargs):
        BasicPlot.__init__(self, *args, **kwargs)
        self.x_data_sets = []
        self.y_data_sets = []
        self.kwarg_sets = []

    def add_content(self, x_data, y_data, **kwargs):
        self.x_data_sets.append(x_data)
        self.y_data_sets.append(y_data)
        self.kwarg_sets.append(kwargs)

    @method_callable_once
    def finalize(self, loc=0, orientation="vertical", alpha=0.8):
        # make plot
        for x_vals, y_vals, plot_kwargs, c in \
            itertools.izip(self.x_data_sets, self.y_data_sets,
                           self.kwarg_sets, default_colors()):
            plot_kwargs.setdefault("align", "center")
            plot_kwargs.setdefault("color", c)
            # FIXME: linewidth is not a valid kwarg in matplotlib 0.87.7
            # Reenable once clusters upgrade to CentOS 5.  Until then,
            # all bars have thick, black borders.
            #plot_kwargs.setdefault("linewidth", 0)
            plot_kwargs.setdefault("orientation", orientation)
            self.ax.bar(x_vals, y_vals, **plot_kwargs)

        # add legend if there are any non-trivial labels
        self.add_legend_if_labels_exist(loc=loc, alpha=alpha)

        # decrement reference counts
        del self.x_data_sets
        del self.y_data_sets
        del self.kwarg_sets

class VerticalBarHistogram(BasicPlot):
    """
    Histogram data sets with a common binning, then make a vertical bar plot.
    """
    def __init__(self, *args, **kwargs):
        BasicPlot.__init__(self, *args, **kwargs)
        self.data_sets = []
        self.kwarg_sets = []

    def add_content(self, data, **kwargs):
        self.data_sets.append(data)
        self.kwarg_sets.append(kwargs)

    @method_callable_once
    def finalize(self, loc=0, alpha=0.8, num_bins=20, normed=False,\
                 logx=False, logy=False):
        # determine binning
        min_stat, max_stat = determine_common_bin_limits(self.data_sets)
        if logx:
            bins = numpy.logspace(numpy.log10(min_stat), numpy.log10(max_stat),\
                                  num_bins + 1, endpoint=True)
        else:
            bins = numpy.linspace(min_stat, max_stat, num_bins+1, endpoint=True)

        # determine bar width; gets silly for more than a few data sets
        if logx:
            width = list(numpy.diff(bins))
        else:
            width = (1 - 0.1 * len(self.data_sets)) * (bins[1] - bins[0])

        # make plot
        for i, (data_set, plot_kwargs, c) in \
            enumerate(itertools.izip(self.data_sets, self.kwarg_sets,\
                                     default_colors())):
            # set default values
            plot_kwargs.setdefault("alpha", 0.6)
            plot_kwargs.setdefault("width", width)
            plot_kwargs.setdefault("color", c)

            # make histogram
            y, x = numpy.histogram(data_set, bins=bins, normed=normed)
            x = x[:-1]

            # mask zeros for logy
            if logy:
                y = numpy.ma.masked_where(y==0, y, copy=False)

            # plot
            plot_item = self.ax.bar(x, y, **plot_kwargs)

        # add legend if there are any non-trivial labels
        self.add_legend_if_labels_exist(loc=loc, alpha=alpha)

        if logx:
            self.ax.set_xscale("log")
        if logy:
            self.ax.set_yscale("log")

        # decrement reference counts
        del self.data_sets
        del self.kwarg_sets

class NumberVsBinBarPlot(BasicPlot):
    """
    Make a bar plot in which the width and placement of the bars are set
    by the given bins.
    """
    def __init__(self, *args, **kwargs):
        BasicPlot.__init__(self, *args, **kwargs)
        self.bin_sets = []
        self.value_sets = []
        self.kwarg_sets = []

    def add_content(self, bins, values, **kwargs):
        if len(bins) != len(values):
            raise ValueError, "length of bins and values do not match"

        self.bin_sets.append(bins)
        self.value_sets.append(values)
        self.kwarg_sets.append(kwargs)

    @method_callable_once
    def finalize(self, orientation="vertical"):
        for bins, values, plot_kwargs in itertools.izip(self.bin_sets,
            self.value_sets, self.kwarg_sets):
            x_vals = bins.centres()

            # prevent each bar from getting a separate legend entry
            label = "_nolegend_"
            if "label" in plot_kwargs:
                label = plot_kwargs["label"]
                del plot_kwargs["label"]

            # set default
            plot_kwargs.setdefault("align", "center")

            if orientation == "vertical":
                plot_kwargs.setdefault("width", bins.upper() - bins.lower())
                patches = self.ax.bar(x_vals, values, **plot_kwargs)
            elif orientation == "horizontal":
                plot_kwargs.setdefault("height", bins.upper() - bins.lower())
                patches = self.ax.barh(x_vals, values, **plot_kwargs)
            else:
                raise ValueError, orientation + " must be 'vertical' " \
                    "or 'horizontal'"

            # prevent each bar from getting a separate legend entry
            if len(patches) > 0:
                patches[0].set_label(label)

        pylab.axis('tight')

        # add legend if there are any non-trivial labels
        self.add_legend_if_labels_exist()

        # decrement reference counts
        del self.bin_sets
        del self.value_sets
        del self.kwarg_sets

class CumulativeHistogramPlot(BasicPlot):
    """
    Cumulative histogram of foreground that also has a shaded region,
    determined by the mean and standard deviation of the background
    population coincidence statistics.
    """
    def __init__(self, *args, **kwargs):
        BasicPlot.__init__(self, *args, **kwargs)
        self.fg_data_sets = []
        self.fg_kwarg_sets = []
        self.bg_data_sets = []
        self.bg_kwargs = {}

    def add_content(self, fg_data_set, **kwargs):
        self.fg_data_sets.append(fg_data_set)
        self.fg_kwarg_sets.append(kwargs)

    def add_background(self, bg_data_sets, **kwargs):
        self.bg_data_sets.extend(bg_data_sets)
        self.bg_kwargs = kwargs

    @method_callable_once
    def finalize(self, num_bins=20, normalization=1):
        epsilon = 1e-8

        # determine binning
        min_stat, max_stat = determine_common_bin_limits(\
            self.fg_data_sets + self.bg_data_sets)
        bins = numpy.linspace(min_stat, max_stat, num_bins + 1, endpoint=True)
        bins_bg = numpy.append(bins, float('Inf'))
        dx = bins[1] - bins[0]

        # plot foreground
        for data_set, plot_kwargs in \
            itertools.izip(self.fg_data_sets, self.fg_kwarg_sets):
            # make histogram
            y, x = numpy.histogram(data_set, bins=bins)
            y = y[::-1].cumsum()[::-1]
            x = x[:-1]

            # plot
            y = numpy.array(y, dtype=numpy.float32)
            y[y <= epsilon] = epsilon
            self.ax.plot(x + dx/2, y*normalization, **plot_kwargs)

        # shade background region
        if len(self.bg_data_sets) > 0:
            # histogram each background instance separately and take stats
            hist_sum = numpy.zeros(len(bins), dtype=float)
            sq_hist_sum = numpy.zeros(len(bins), dtype=float)
            for instance in self.bg_data_sets:
                # make histogram
	        y, x = numpy.histogram(instance, bins=bins_bg)
                x = numpy.delete(x, -1)
                y = y[::-1].cumsum()[::-1]
                hist_sum += y
                sq_hist_sum += y*y

            # get statistics
            N = len(self.bg_data_sets)
            means = hist_sum / N
            stds = numpy.sqrt((sq_hist_sum - hist_sum*means) / (N - 1))

            # plot mean
            means[means <= epsilon] = epsilon
            self.ax.plot(x + dx/2, means*normalization, 'r+', **self.bg_kwargs)

            # shade in the area
            if "label" in self.bg_kwargs:
                self.bg_kwargs["label"] = r"$\mu_\mathrm{%s}$" \
                    % self.bg_kwargs["label"]
            self.bg_kwargs.setdefault("alpha", 0.3)
            self.bg_kwargs.setdefault("facecolor", "y")
            upper = means + stds
            lower = means - stds
            lower[lower <= epsilon] = epsilon
            tmp_x, tmp_y = viz.makesteps(bins, upper, lower)
            self.ax.fill(tmp_x, tmp_y*normalization, **self.bg_kwargs)

        # make semilogy plot
        self.ax.set_yscale("log")

        # adjust plot range
        self.ax.set_xlim((0.9 * min_stat, 1.1 * max_stat))
        possible_ymins = [0.6]
        if len(self.bg_data_sets) > 0:
            possible_ymins.append(0.6 / N)
        else:
            possible_ymins.append(0.6 * normalization)
        self.ax.set_ylim(min(possible_ymins))

        # add legend if there are any non-trivial labels
        self.kwarg_sets = self.fg_kwarg_sets
        self.add_legend_if_labels_exist()

        # decrement reference counts
        del self.kwarg_sets
        del self.fg_data_sets
        del self.fg_kwarg_sets
        del self.bg_data_sets
        del self.bg_kwargs


class ImagePlot(BasicPlot):
    """
    The equivalent of pylab.imshow(), but with the BasicPlot niceties and
    some defaults that are more in tune with what a scientist wants --
    origin="lower", requiring x and y bins so that we can label axes
    correctly, and a colorbar.
    """
    def __init__(self, *args, **kwargs):
        colorlabel = kwargs.pop("colorlabel", None)
        BasicPlot.__init__(self, *args, **kwargs)
        self.image_sets  = []
        self.x_bins_sets = []
        self.y_bins_sets = []
        self.kwarg_sets  = []
        self.colorlabel  = colorlabel

    def add_content(self, image, x_bins, y_bins, **kwargs):
        """
        Add a given image to this plot.

        @param image: two-dimensional numpy array to plot
        @param x_bins: pylal.rate.Bins instance describing the x binning
        @param y_bins: pylal.rate.Bins instance describing the y binning
        """
        if image.ndim != 2:
            raise ValueError, "require 2-D array"
        self.image_sets.append(image)
        self.x_bins_sets.append(x_bins)
        self.y_bins_sets.append(y_bins)
        self.kwarg_sets.append(kwargs)

    @method_callable_once
    def finalize(self, colorbar=True, logcolor=False, minorticks=False,\
                 clim=None):
        from pylal import rate

        logx = False
        logy = False

        for image,x_bins,y_bins,plot_kwargs in\
            itertools.izip(self.image_sets, self.x_bins_sets, self.y_bins_sets,\
                           self.kwarg_sets):
            extent = [x_bins.lower()[0], x_bins.upper()[-1],
                      y_bins.lower()[0], y_bins.upper()[-1]]
            plot_kwargs.setdefault("origin", "lower")
            plot_kwargs.setdefault("interpolation", "nearest")
            if logcolor and clim:
                plot_kwargs.setdefault(\
                    "norm", pylab.matplotlib.colors.LogNorm(vmin=clim[0],\
                                                            vmax=clim[1]))
            elif logcolor:
                plot_kwargs.setdefault("norm",\
                                       pylab.matplotlib.colors.LogNorm())
            elif clim:
                plot_kwargs.setdefault(\
                    "norm", pylab.matplotlib.colors.Normalize(vmin=clim[0],\
                                                              vmax=clim[1]))
            if colorbar:
                plot_kwargs.setdefault("cmap", \
                     pylab.matplotlib.colors.LinearSegmentedColormap("clrs",\
                                       pylab.matplotlib.cm.jet._segmentdata))

            im = self.ax.imshow(image, extent=extent, **plot_kwargs)

            # set log axes, raise error if you try to use both log and lin on
            # the same plot
            if isinstance(x_bins, rate.LogarithmicBins):
                logx = True
            if logx and not isinstance(x_bins, rate.LogarithmicBins):
                raise ValueError("Cannot process both linear and logarithmic "+\
                                 "Images on the same Axis.")
            if isinstance(y_bins, rate.LogarithmicBins):
                logy = True
            if logy and not isinstance(y_bins, rate.LogarithmicBins):
                raise ValueError("Cannot process both linear and logarithmic "+\
                                 "Images on the same Axis.")

        pylab.axis('tight')

        if colorbar:
            add_colorbar(self.ax, log=logcolor, clim=clim,\
                         label=self.colorlabel)

        if logx:
            xticks, xlabels, xminorticks = log_transform(self.ax.get_xlim())
            self.ax.set_xticks(xticks)
            self.ax.set_xticklabels(xlabels)
            if minorticks: self.ax.set_xticks(xminorticks, minor=True)
        if logy:
            yticks, ylabels, yminorticks = log_transform(self.ax.get_ylim())
            self.ax.set_yticks(yticks)
            self.ax.set_yticklabels(ylabels)
            if minorticks: self.ax.set_yticks(yminorticks, minor=True)

class FillPlot(BasicPlot):
    """
    Given a list of vertices (passed by x-coords and y-coords), fill the
    regions (by default with successively darker gray shades).
    """
    def __init__(self, *args, **kwargs):
        BasicPlot.__init__(self, *args, **kwargs)
        self.x_coord_sets = []
        self.y_coord_sets = []
        self.kwarg_sets = []
        self.shades = []

    def add_content(self, x_coords, y_coords, shade=None, **kwargs):
        if len(x_coords) != len(y_coords):
            raise ValueError, "x and y coords have different length"
        if iterutils.any(s is None for s in self.shades) and shade is not None \
        or iterutils.any(s is not None for s in self.shades) and shade is None:
            raise ValueError, "cannot mix explicit and automatic shading"

        self.x_coord_sets.append(x_coords)
        self.y_coord_sets.append(y_coords)
        self.kwarg_sets.append(kwargs)
        self.shades.append(shade)

    @method_callable_once
    def finalize(self):
        # fill in default shades if necessary
        if iterutils.any(s is None for s in self.shades):
            n = len(self.shades)
            grays = numpy.linspace(0, 1, n, endpoint=False)[::-1]
            self.shades = numpy.vstack((grays, grays, grays)).T

        # plot
        for x, y, s, plot_kwargs in zip(self.x_coord_sets, self.y_coord_sets,
                              self.shades, self.kwarg_sets):
            self.ax.fill(x, y, facecolor=s, **plot_kwargs)

        # add legend if there are any non-trivial labels
        self.add_legend_if_labels_exist(loc=0)

class SixStripSeriesPlot(BasicPlot):
    """
    Given a time- or frequency-series, plot it across six horizontal axes,
    stacked on top of one another.  This is good for representing a
    high-resolution one-dimensional data set.  To support missing data,
    we require x and y coordinates for each data set.
    """
    def __init__(self, xlabel="", ylabel="", title=""):
        self.fig = pylab.figure(figsize=(5.54, 7.5))
        self.title = title

        # do not want to create axes yet
        self.xlabel = xlabel
        self.ylabel = ylabel

        self.x_coord_sets = []
        self.y_coord_sets = []
        self.kwarg_sets = []

    def add_content(self, x_coords, y_coords, **kwargs):
        if len(x_coords) != len(y_coords):
            raise ValueError, "x and y coords have different length"
        if len(self.kwarg_sets) and \
            (iterutils.any("color" in kw for kw in self.kwarg_sets)
             ^ ("color" in kwargs)):
            raise ValueError, "cannot mix explicit and automatic coloring"

        self.x_coord_sets.append(x_coords)
        self.y_coord_sets.append(y_coords)
        self.kwarg_sets.append(kwargs)

    @method_callable_once
    def finalize(self, yscale="linear"):
        for kw, c in zip(self.kwarg_sets, default_colors()):
            kw.setdefault("color", c)

        min_x, max_x = determine_common_bin_limits(self.x_coord_sets)

        numaxes = 6  # This is hardcoded below.  Change at your own risk.
        ticks_per_axis = 6 # Hardcoded, but you can change it if it looks bad.
        fperaxis = (max_x - min_x) / numaxes
        freq_boundaries = numpy.linspace(min_x, max_x, numaxes + 1)
        freq_ranges = zip(freq_boundaries[:-1], freq_boundaries[1:])

        # attempt to put ticks at every multiple of 10 unless there are too few
        # or too many
        tickspacing = int(numpy.ceil(fperaxis / ticks_per_axis / 10)) * 10
        if abs(fperaxis / tickspacing - ticks_per_axis) > 2:
            tickspacing = int(numpy.ceil(fperaxis / ticks_per_axis))

        # iterate over axes
        for j, freq_range in enumerate(freq_ranges):
            # create one of the 6 axes
            ax = self.fig.add_axes([.12, .84-.155*j, .86, 0.124])

            for x_coords, y_coords, kwarg_set in zip(self.x_coord_sets,
                self.y_coord_sets, self.kwarg_sets):
                # just look at the region relevant to our axis
                ind = (x_coords >= freq_range[0]) & (x_coords < freq_range[1])
                x = x_coords[ind]
                y = y_coords[ind]

                # add data to axes
                ax.plot(x, y, **kwarg_set)

                # fix the limits and ticks
                ax.set_xlim(freq_range)

                mintick = int(numpy.ceil(freq_range[0] / tickspacing)) \
                    * tickspacing
                maxtick = int(numpy.floor(freq_range[1] / tickspacing)) \
                    * tickspacing
                ax.set_xticks(xrange(mintick, maxtick+1, tickspacing))

                ax.set_ylabel(self.ylabel)
                ax.set_yscale(yscale)
                ax.grid(True)

        # label bottom row
        ax.set_xlabel(self.xlabel)

        # title top row
        self.fig.axes[0].set_title(self.title)

        # apply common y limits
        y_lims = [ax.get_ylim() for ax in self.fig.axes]
        new_y_lims = determine_common_bin_limits(y_lims)
        for ax in self.fig.axes:
            ax.set_ylim(new_y_lims)

class ROCPlot(BasicPlot):
    """
    Plot the receiver operating characteristic (ROC) based on the foreground
    and background values from given techniques.  For example, to compare
    SNR vs IFAR, do something like:

    plot = ROCPlot("FAP", "EFF", "ROC IFAR vs SNR")
    plot.add_content(ifar_bg, ifar_fg, label="IFAR")
    plot.add_content(snr_bg, snr_fg, label="SNR")
    plot.finalize()
    """
    def __init__(self, *args, **kwargs):
        BasicPlot.__init__(self, *args, **kwargs)
        self.bg_sets = []
        self.fg_sets = []
        self.eff_weight_sets = []
        self.kwarg_sets = []

    def add_content(self, bg, fg, eff_weight=None, **kwargs):
        """
        Enter a particular technique's background, foreground, and efficiency
        weights.  These should be one-dimensional arrays with the values of
        of your statistics for backgrounds and foregrounds.  Eff_weight
        are weights on the efficiency, useful if, say, you have a different
        prior than your injection set reflects.
        """
        if eff_weight is not None and len(fg) != len(eff_weight):
            raise ValueError, "efficiency weights and foreground values "\
                "must be the same length"
        self.bg_sets.append(bg)
        self.fg_sets.append(fg)
        self.eff_weight_sets.append(eff_weight)
        self.kwarg_sets.append(kwargs)

    @method_callable_once
    def finalize(self, loc=0):
        for bg, fg, weights, kwargs in itertools.izip(self.bg_sets,
            self.fg_sets, self.eff_weight_sets, self.kwarg_sets):
            # sort, keeping the weights and fg values together
            bg = numpy.array(bg)  # always copies
            bg.sort()
            fg = numpy.array(fg)  # always copies
            if weights is not None:
                weights = numpy.array(weights)[fg.argsort()]
            fg.sort()

            # calculate false alarm probability and efficiency
            FAP = 1 - numpy.arange(len(bg), dtype=float) / len(bg)
            if weights is not None:
                EFF = weights[::-1].cumsum()[::-1]
            else:
                EFF = 1 - numpy.arange(len(fg), dtype=float) / len(fg)

            # now find the efficiency *at* each false alarm probability
            # if louder than loudest event, then efficiency is 0
            EFF_ind = numpy.array([fg.searchsorted(x) for x in bg])
            fg_too_loud = EFF_ind == len(fg)
            EFF_ind[fg_too_loud] = 0  # dummy index
            EFF_by_FAP = EFF[EFF_ind]
            EFF_by_FAP[fg_too_loud] = 0.

            # plot!
            self.ax.plot(FAP, EFF_by_FAP, **kwargs)

        # make it pretty
        self.ax.grid(True)
        self.ax.set_xlim((0, 1))
        self.ax.set_ylim((0, 1))

        # resize figure to make axes square
        fig_side = min(self.fig.get_size_inches())
        self.fig.set_size_inches(fig_side, fig_side)

        # add legend if there are any non-trivial labels
        self.add_legend_if_labels_exist(loc=loc)

        # decrement reference counts
        del self.fg_sets
        del self.bg_sets
        del self.eff_weight_sets
        del self.kwarg_sets

class QQPlot(BasicPlot):
    """
    Plot the global rank versus the self rank, like a Q-Q plot, i.e.
    differences between the probability distribution of a
    statistical population to a test population.
    """
    def __init__(self, *args, **kwargs):
        BasicPlot.__init__(self, *args, **kwargs)
        self.fg_data = []
        self.bg_data = []
        self.fg_kwargs = []
        self.bg_kwargs = []
        self.kwargs = None

    def add_content(self, *args):
        raise NotImplementedError, "A function 'add_content' is not "\
              "implemented for QQPlot. "\
              "Please use the functions 'add_fg' and 'add_bg'."

    def _compute_quantiles(self, fg, bg):
        """
        Calculates the sorted quantiles to be plotted.
        """
        from scipy import stats

        N_FG = len(fg)
        N_BG = len(bg)
        N_tot = N_BG + N_FG

        # global rank vs self rank; very much like a Q-Q plot
        bg_self_quantile = stats.rankdata(bg).astype(float) / N_BG
        fg_self_quantile = stats.rankdata(fg).astype(float) / N_FG

        popAB = numpy.concatenate((fg, bg))
        popAB_ranked = stats.rankdata(popAB).astype(float)
        fg_global_quantile = popAB_ranked[:N_FG] / N_tot
        bg_global_quantile = popAB_ranked[N_FG:] / N_tot

        fg_self_quantile.sort()
        bg_self_quantile.sort()
        fg_global_quantile.sort()
        bg_global_quantile.sort()

        return fg_self_quantile, bg_self_quantile, fg_global_quantile, bg_global_quantile

    def add_fg(self, data, **kwargs):
        """
        Add the foreground data and the kwargs to a list for the
        """
        self.fg_data.append(data)
        self.fg_kwargs.append(kwargs)

    def add_bg(self, data, **kwargs):
        """
        Add the foreground data and the kwargs to a list for the
        """
        self.bg_data.append(data)
        self.bg_kwargs.append(kwargs)

    @method_callable_once
    def finalize(self, loc=0):

        if len(self.fg_data)!=len(self.bg_data):
            raise ValueError, "You have to add the same number of foreground data"\
                  " as background data. But you have specified %d sets of foreground data but %d"\
                  " sets of background data. " % ( len(self.fg_data), len(self.bg_data))

        for fg, bg, fg_kwargs, bg_kwargs in itertools.izip(self.fg_data, self.bg_data, \
                                                           self.fg_kwargs, self.bg_kwargs):

            fg_self, bg_self, fg_global, bg_global = \
                     self._compute_quantiles(fg, bg)

            self.ax.plot(bg_self, bg_global, **bg_kwargs)
            self.ax.plot(fg_self, fg_global, **fg_kwargs)

        # make it pretty
        self.ax.grid(True)

        # add legend if there are any non-trivial labels
        self.ax.legend(loc=4)

        # decrement reference counts
        del self.fg_data
        del self.bg_data
        del self.fg_kwargs
        del self.bg_kwargs
        del self.kwargs

class SimpleMapPlot(BasicPlot):
    """
    Class to create a clickable map html page.
    """
    def __init__(self, *args, **kwargs):
        BasicPlot.__init__(self, *args, **kwargs)
        self.x_data_sets = []
        self.y_data_sets = []
        self.link_list = []
        self.kwarg_sets = []

        self.click_size = 5

        self.click_x = []
        self.click_y = []
        self.click_link = []

    def set_click_size(self, size):
        """
        Sets the size of the area around a point that can be clicked at.
        """
        self.click_size = size

    def add_content(self, x_data, y_data, link, **kwargs):
        self.x_data_sets.append(x_data)
        self.y_data_sets.append(y_data)
        self.link_list.append(link)
        self.kwarg_sets.append(kwargs)


    @method_callable_once
    def finalize(self, loc=0):

        # make plot
        for x_vals, y_vals, plot_kwargs in \
            itertools.izip(self.x_data_sets, self.y_data_sets, self.kwarg_sets):
            self.ax.plot(x_vals, y_vals, **plot_kwargs)

        # add legend if there are any non-trivial labels
        self.add_legend_if_labels_exist(loc=loc)


    def rescale(self, dpi):
        """
        Calculate the rescaling to map the point coordinates
        to the actual coordinates on the image.
        """

        # set the dpi used for the savefig
        if dpi is None:
            dpi = pylab.rcParams['savefig.dpi']
        self.ax.get_figure().set_dpi(dpi)

        # get the full extend of the final image
        _, _, width, height = self.fig.bbox.extents

        # set the pixels for each point that can be clicked
        self.click_x = []
        self.click_y = []
        self.click_link = []
        for xvec, yvec, link in \
                itertools.izip(self.x_data_sets, self.y_data_sets, self.link_list):
            # skip if no link is associated with these point(s)
            if link is None:
                continue

            for point in zip(xvec, yvec):

                # transform the data coordinates into piel coordinates
                pixel = self.ax.transData.transform(point)

                # save the new coordinates. The y-coordinate is just
                # the other way around.
                self.click_x.append(pixel[0])
                self.click_y.append(height - pixel[1])
                self.click_link.append(link)


    def create_html(self, plotname, dpi = None):
        """
        Create the html file with the name of the plot
        as the only additional input information needed.
        If the actual plot is saved with a different
        dpi than the standard one, it must be specified here!!
        """

        # points need to be rescaled first
        self.rescale(dpi)

        # get the full extend of the final image
        _, _, width, height = self.fig.bbox.extents

        # create the map-page
        page = ''
        page += '<IMG src="%s" width=%dpx '\
            'usemap="#map">' % ( plotname, width)
        page +=  '<MAP name="map"> <P>'
        n=0
        for px, py, link in zip( self.click_x, self.click_y, self.click_link):
            n+=1
            page +=  '<area href="%s" shape="circle" '\
                'coords="%d, %d, 5"> Point%d</a>' %\
                ( link, px, py, n)
        page += '</P></MAP></OBJECT><br>'
        page += "<hr/>"

        # return the html code, to be saved to a file or embedded in a
        # larger html page
        return page



plot3D_head = """
  <style type="text/css">
    .rotationViewer {
      position:relative;
      width:640px;
      height:378px;
      border-style:solid;
      border-width:1px;
      margin:auto;
      margin-bottom:10px;
      cursor:pointer;
    }
    .rotationViewer img {
      position:absolute;
      visibility:hidden;
      left:0;
      top:0;
      width:100%;
      height:100%;
    }
  </style>
"""

plot3D_body = """
<script type="text/javascript" src="PATHTOUIZE/Uize.js"></script>

<div class="main">
  <div id="page_rotationViewer" class="rotationViewer insetBorderColor"></div>
</div>

<!-- JavaScript code to make the static bar HTML "come alive" -->
<!-- Based on http://www.uize.com/examples/3d-rotation-viewer.html -->
<script type="text/javascript">

Uize.module ({
  required:[
    'UizeDotCom.Page.Example.library',
    'UizeDotCom.Page.Example',
    'Uize.Widget.Drag',
    'Uize.Fade.xFactory',
    'Uize.Curve.Rubber',
    'Uize.Curve.Mod'
  ],
  builder:function () {
    /*** create the example page widget ***/
      var page = window.page = Uize.Widget.Page  ({evaluator:function (code) {eval (code)}});

    /*** configuration variables ***/
      var
        totalFrames = TOTALFRAMES,
        frameUrlTemplate = 'URLTEMPLATE'
      ;

    /*** state variables ***/
      var
        rotation = 0,
        lastFrameNo = -1,
        dragStartRotation
      ;

    /*** create the Uize.Widget.Drag instance ***/
      var rotationViewer = page.addChild (
        'rotationViewer',
        Uize.Widget.Drag,
        {
          cancelFade:{duration:5000,curve:Uize.Curve.Rubber.easeOutBounce ()},
          releaseTravel:function (speed) {
            var
              deceleration = 5000, // measured in pixels/s/s
              duration = speed / deceleration
            ;
            return {
              duration:duration,
              distance:Math.round (speed * duration / 2),
              curve:function (_value) {return 1 - (_value = 1 - _value) * _value}
            };
          },
          html:function (input) {
            var
              htmlChunks = [],
              frameNodeIdPrefix = input.idPrefix + '-frame'
            ;
            for (var frameNo = 0; ++frameNo <= totalFrames;) {
              htmlChunks.push (
                '<img' +
                  ' id="' + frameNodeIdPrefix + frameNo + '"' +
                  ' src="' + Uize.substituteInto (frameUrlTemplate,{frame:(frameNo < 10 ? '0' : '') + frameNo}) +'"' +
                '/>'
              );
            }
            return htmlChunks.join ('');
          },
          built:false
        }
      );

    /*** wire up the drag widget with events for updating rotation degree ***/
      function updateRotation (newRotation) {
        rotation = ((newRotation % 360) + 360) % 360;
        var frameNo = 1 + Math.round (rotation / 360 * (totalFrames - 1));
        if (frameNo != lastFrameNo) {
          rotationViewer.showNode ('frame'+ lastFrameNo,false);
          rotationViewer.showNode ('frame'+ (lastFrameNo = frameNo));
        }
      }
      rotationViewer.wire ({
        'Drag Start':function () {dragStartRotation = rotation},
        'Drag Update':function (e) {updateRotation (dragStartRotation - e.source.eventDeltaPos [0] / 2.5)}
      });

    /*** function for animating spin ***/
      function spin (degrees,duration,curve) {
        Uize.Fade.fade (updateRotation,rotation,rotation + degrees,duration,{quantization:1,curve:curve});
      }

    /*** initialization ***/
      Uize.Node.wire (window,'load',function () {spin (360,2700,Uize.Curve.easeInOutPow (4))});

    /*** wire up the page widget ***/
      page.wireUi ();
  }
});

</script>
"""


class Plot3D(BasicPlot):
    """
    Exactly what you get by calling pylab.plot(), but with the handy extras
    of the BasicPlot class.
    """
    def __init__(self, *args, **kwargs):


        BasicPlot.__init__(self)
        self.x_data_sets = []
        self.y_data_sets = []
        self.z_data_sets = []
        self.kwarg_sets = []

        # need the 3D axes
        self.ax =  Axes3D(self.fig)

        # set the labels on the new 3D axis
        self.ax.set_xlabel(args[0])
        self.ax.set_ylabel(args[1])
        self.ax.set_zlabel(args[2])
        self.ax.set_title(args[3])

        self.ax.grid(True)

    def add_content(self, x_data, y_data, z_data, **kwargs):
        self.x_data_sets.append(x_data)
        self.y_data_sets.append(y_data)
        self.z_data_sets.append(z_data)
        self.kwarg_sets.append(kwargs)


    @method_callable_once
    def finalize(self, plot_basename, n=50, dpi=50, loc=0, plot_dir = '.',\
                     js_dir = 'https://ldas-jobs.phys.uwm.edu/~dietz/js'):

        # check input argument
        if n>100:
            raise ValueError("The 3D code cannot handle more than 100 frames currently.")

        fig = pylab.figure()
        ax = Axes3D(fig)

        # create plot
        for x_vals, y_vals, z_vals, plot_kwargs in \
            itertools.izip(self.x_data_sets, self.y_data_sets, \
                           self.z_data_sets, self.kwarg_sets):
            self.ax.scatter(x_vals, y_vals, z_vals, **plot_kwargs)

        # add legend if there are any non-trivial labels
        self.add_legend_if_labels_exist(loc=loc)

        # loop over a full circle with n steps
        dphi = 360.0/float(n)
        for i, angle in enumerate(numpy.arange(0.0, 360.0, dphi)):

            # make the rotation of the image
            self.ax.view_init(30.0, angle)

            # create the rotated image in the ploting directory
            picname = '%s/%s-%02d'%(plot_dir, plot_basename, i)
            self.savefig(picname+'.png',dpi=dpi)


        # decrement reference counts
        del self.x_data_sets
        del self.y_data_sets
        del self.z_data_sets
        del self.kwarg_sets


        # produce the html output code
        plot_template = '%s/%s-[#frame].png'%(plot_dir, plot_basename)
        body = plot3D_body.replace('TOTALFRAMES',str(n))
        body = body.replace('URLTEMPLATE', plot_template)
        body = body.replace('PATHTOUIZE', js_dir)

        # return the html snippets
        return plot3D_head, body

class ScatterPlot(SimplePlot):
    """
    Exactly what you get from calling pylab.scatter(), but with the handy extras
    of the BasicPlot class.
    """

    @method_callable_once
    def finalize(self, loc=0, alpha=0.8):
        # make plot
        for x_vals, y_vals, plot_kwargs, color in \
            itertools.izip(self.x_data_sets, self.y_data_sets, self.kwarg_sets,\
                           default_colors()):
            plot_kwargs.setdefault("c", color)
            if (len(x_vals) and
                (isinstance(y_vals, numpy.ma.MaskedArray) and y_vals.count() or
                 True)):
                self.ax.scatter(x_vals, y_vals, **plot_kwargs)
            else:
                plot_kwargs["visible"] = False
                self.ax.scatter([1], [1], **plot_kwargs)

        # add legend if there are any non-trivial labels
        self.add_legend_if_labels_exist(loc=loc, alpha=0.8)

        # decrement reference counts
        del self.x_data_sets
        del self.y_data_sets
        del self.kwarg_sets

class ColorbarScatterPlot(BasicPlot):
    """
    Exactly what you get from calling pylab.scatter() when you want a colorbar,
    but with the handy extras of the BasicPlot class.
    """

    def __init__(self, xlabel="", ylabel="", clabel="", title="", subtitle=""):
        BasicPlot.__init__(self, xlabel=xlabel, ylabel=ylabel, title=title,\
                           subtitle=subtitle)
        self.clabel      = clabel
        self.x_data_sets = []
        self.y_data_sets = []
        self.c_data_sets = []
        self.kwarg_sets  = []

    def add_content(self, x_data, y_data, c_data, **kwargs):
        self.x_data_sets.append(x_data)
        self.y_data_sets.append(y_data)
        self.c_data_sets.append(c_data)
        self.kwarg_sets.append(kwargs)

    def finalize(self, loc=0, colorbar=True, logcolor=False, clim=None,\
                 alpha=0.8):
        # make plot
        p = None
        for x_vals, y_vals, c_vals, plot_kwargs in\
            itertools.izip(self.x_data_sets, self.y_data_sets, self.c_data_sets,
                           self.kwarg_sets):
            if len(x_vals):
                zipped = zip(x_vals, y_vals, c_vals)
                zipped.sort(key=lambda (x,y,c): c)
                x_vals, y_vals, c_vals = map(numpy.asarray, zip(*zipped))
            if logcolor:
                plot_kwargs.setdefault("norm",pylab.matplotlib.colors.LogNorm())
            try:
                p = self.ax.scatter(x_vals, y_vals, c=c_vals, **plot_kwargs)
            except ValueError:
                plot_kwargs['visible'] = False
                p = self.ax.scatter([1], [1], c=[1], **plot_kwargs)
            p = self.ax.scatter(x_vals, y_vals, c=c_vals, **plot_kwargs)

        if colorbar and p is not None:
            add_colorbar(self.ax, mappable=p, log=logcolor, label=self.clabel,\
                         clim=clim)

        # add legend if there are any non-trivial labels
        self.add_legend_if_labels_exist(loc=loc, alpha=alpha)

        # decrement reference counts
        del self.x_data_sets
        del self.y_data_sets
        del self.c_data_sets
        del self.kwarg_sets

class PlotSegmentsPlot(BasicPlot):
    """
    Horizontal bar segment plot.
    """
    color_code = {'H1':'r', 'H2':'b', 'L1':'g', 'V1':'m', 'G1':'k'}

    def __init__(self, xlabel="", ylabel="", title="", subtitle="", t0=0,\
                 dt=1):
        """
        Create a fresh plot. Provide t0 to provide reference time to use as
        zero, and dt to use a different number of seconds for each unit of the
        x-axis.
        """
        BasicPlot.__init__(self, xlabel, ylabel, title)
        self.segdict = segments.segmentlistdict()
        self.keys = []
        self._time_transform = lambda t: float(t - t0)/dt
        self.kwarg_sets  = []

    def add_content(self, segdict, keys=None, **kwargs):
        if not keys:
            keys = sorted(segdict.keys())
        for key in keys:
            if key in self.segdict.keys():
                self.segdict[key] += segdict[key]
            else:
                self.keys.append(key)
                self.segdict[key] = segdict[key]
            self.kwarg_sets.append(dict())
            self.kwarg_sets[-1].update(kwargs)
                
    def highlight_segment(self, seg, **plot_args):
        """
        Highlight a particular segment with dashed lines.
        """
        a,b = map(self._time_transform,seg)
        plot_args.setdefault('linestyle', '--')
        plot_args.setdefault('color','r')
        self.ax.axvline(a, **plot_args)
        self.ax.axvline(b, **plot_args)

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

    @method_callable_once
    def finalize(self, labels_inset=False, hidden_colorbar=False):

        for row,(key,plot_kwargs)\
        in enumerate(itertools.izip(self.keys, self.kwarg_sets)):
            plot_kwargs.setdefault("edgecolor", "k")
            if self.color_code.has_key(key):
                plot_kwargs.setdefault("color", self.color_code[key])
            else:
                plot_kwargs.setdefault("color", "b")
            for seg in self.segdict[key]:
                a,b = map(self._time_transform, seg)
                self.ax.fill([a, b, b, a, a],\
                             [row-0.4, row-0.4, row+0.4, row+0.4, row-0.4],\
                             **plot_kwargs)
            if labels_inset:
                self.ax.text(0.01, (row+1)/(len(self.keys)+1),\
                             re.sub(r'\\+_+','\_',key),\
                             horizontalalignment='left',\
                             verticalalignment='center',\
                             transform=self.ax.transAxes,\
                             backgroundcolor='white',\
                             bbox=dict(facecolor='white', alpha=0.5,\
                                       edgecolor='white'))
    
        ticks = pylab.arange(len(self.keys))
        self.ax.set_yticks(ticks)
        if labels_inset:
            self.ax.set_yticklabels(ticks, color='white')
        else:
            self.ax.set_yticklabels([re.sub(r'\\+_+', '\_', k)\
                                       for k in self.keys], size='small')
        self.ax.set_ylim(-1, len(self.keys))

        # dereference
        del self.keys
        del self.segdict
        del self.kwarg_sets

class DQScatterPlot(ColorbarScatterPlot):
    """
    A DQ-formatted scatter plot, with those triggers below some threshold on the
    colour axis get plotted tiny wee.
    """
    def finalize(self, loc=0, alpha=0.8, colorbar=True, logcolor=False,\
                 clim=None, threshold=0):
        # make plot
        p = None
        for x_vals, y_vals, c_vals, plot_kwargs in\
            itertools.izip(self.x_data_sets, self.y_data_sets, self.c_data_sets,
                           self.kwarg_sets):
            plot_kwargs.setdefault("s", 15)
            if logcolor:
                plot_kwargs.setdefault("norm",pylab.matplotlib.colors.LogNorm())
            if clim:
                plot_kwargs.setdefault("vmin", clim[0])
                plot_kwargs.setdefault("vmax", clim[1])

            if len(x_vals):
                zipped = zip(x_vals, y_vals, c_vals)
                zipped.sort(key=lambda (x,y,c): c)
                x_vals, y_vals, c_vals = map(numpy.asarray, zip(*zipped))

            t = (c_vals>=threshold).nonzero()[0]
            if len(t):
                idx = t[0]
                xlow,xhigh = numpy.split(x_vals, [idx])
                ylow,yhigh = numpy.split(y_vals, [idx])
                clow,chigh = numpy.split(c_vals, [idx])
            else:
                xlow  = x_vals
                ylow  = y_vals
                clow  = c_vals
                xhigh = numpy.array([])
                yhigh = numpy.array([])
                chigh = numpy.array([])

            if xlow.size > 0:
                lowargs = copy.deepcopy(plot_kwargs)
                lowargs.pop("label", None)
                lowargs["s"] /= 4
                if lowargs.get("marker", None) != "x":
                    lowargs["edgecolor"] = "none"
                self.ax.scatter(xlow, ylow, c=clow, **lowargs)
            if xhigh.size > 0:
                self.ax.scatter(xhigh, yhigh, c=chigh, **plot_kwargs)
            if not len(x_vals):
                plot_kwargs["visible"] = False
                self.ax.scatter([1], [1], c=1, **plot_kwargs)

        if colorbar:
            add_colorbar(self.ax, log=logcolor, label=self.clabel, clim=clim,\
                         cmap=plot_kwargs.get("cmap", None))

        # add legend if there are any non-trivial labels
        self.add_legend_if_labels_exist(loc=loc, alpha=alpha)

        # decrement reference counts
        del self.x_data_sets
        del self.y_data_sets
        del self.c_data_sets
        del self.kwarg_sets

class LineHistogram(BasicPlot):
    """
    A simple line histogram plot. The values of each histogram bin are plotted
    using pylab.plot(), with points centred on the x values and height equal
    to the y values.

    Cumulative, and rate options can be passeed to the finalize() method to
    format each trace individually.
    """
    def __init__(self, xlabel="", ylabel="", title="", subtitle=""):
        BasicPlot.__init__(self, xlabel, ylabel, title, subtitle)
        self.data_sets        = []
        self.normalize_values = []
        self.kwarg_sets       = []

    def add_content(self, data, normalize=1, **kwargs):
        self.data_sets.append(data)
        self.normalize_values.append(normalize)
        self.kwarg_sets.append(kwargs)

    @method_callable_once
    def finalize(self, loc=0, num_bins=100, cumulative=False, rate=False,\
                 xlim=None, bar=False, fill=False, logx=False, logy=False,\
                 alpha=0.8):
        # determine binning
        if xlim:
            min_stat, max_stat = xlim
        else:
            min_stat, max_stat = determine_common_bin_limits(self.data_sets)
        if logx:
            if min_stat == max_stat == 0:
                min_stat = 1
                max_stat = 10
            bins = numpy.logspace(numpy.log10(min_stat), numpy.log10(max_stat),\
                                  int(num_bins) + 1, endpoint=True)
        else:
            bins = numpy.linspace(min_stat, max_stat, num_bins+1, endpoint=True)

        # determine bar width; gets silly for more than a few data sets
        if logx:
            width = list(numpy.diff(bins))
        else:
            width = (1 - 0.1 * len(self.data_sets)) * (bins[1] - bins[0])
        width = numpy.asarray(width)

        # make plot
        ymin = None

        for data_set, norm, plot_kwargs in\
            itertools.izip(self.data_sets, self.normalize_values,\
                           self.kwarg_sets):
            # make histogram
            y, x = numpy.histogram(data_set, bins=bins, normed=False)
            y = y.astype(float)
            x = x[:-1]

            # set defaults
            if fill:
                plot_kwargs.setdefault("linewidth", 1)
                plot_kwargs.setdefault("alpha", 0.8)

            # get cumulative sum
            if cumulative:
                y = y[::-1].cumsum()[::-1]

            # normalize
            if norm != 1:
                y /= norm

            # get bars
            if bar:
                x = numpy.vstack((x-width/2, x+width/2)).reshape((-1,),\
                                                             order="F")
                y = numpy.vstack((y, y)).reshape((-1,), order="F")

            # plot
            if logy:
                numpy.putmask(y, y==0, 1e-100)
            self.ax.plot(x, y, **plot_kwargs)
            if fill:
                plot_kwargs.pop("label", None)
                self.ax.fill_between(x, 1e-100, y, **plot_kwargs)

        if logx:
            try:
                self.ax.set_xscale("log", nonposx='clip')
            except OverflowError:
                self.ax._autoscaleXon = False
                self.ax.set_xscale("log", nonposx='clip')

        if logy:
            try:
                self.ax.set_yscale("log", nonposy='clip')
            except OverflowError:
                self.ax._autoscaleYon = False
                self.ax.set_yscale("log", nonposy='clip')

        # add legend if there are any non-trivial labels
        self.add_legend_if_labels_exist(loc=loc, alpha=0.8)

        # decrement reference counts
        del self.data_sets
        del self.normalize_values
        del self.kwarg_sets

class SkyPositionsPlot(ScatterPlot):
    """
    A spherical projection plot.
    """
    @method_callable_once
    def finalize(self, loc='lower right', projection='moll', centre=(0,0),\
                 frame=[(0, 0), (1, 1)]):
        """
        Finalize and draw plot. Can only be called once.

        Keyword arguments:

            toc : [ str | int ]
                compatible value for legend loc, default: 'lower right'
            projection : str
                valid projection name for mpl_toolkits.basemap.Basemap
            centre : list
                (ra, dec) in degrees for point to centre on plot, only works
                for some projections
            frame : list
                [(xmin, ymin), (xmax, ymax)] pair of tuples for lower left
                and upper right corners of plot area. Full arei (entire
                sphere) is [(0,0), (1,1)]. Narrower range is span in,
                wider range zoomed out.
        """
        # set up projection
        projection = projection.lower()
        centre = [(-(centre[0]-180)+180) % 360, centre[1]]
        m1 = Basemap(projection=projection, lon_0=centre[0], lat_0=centre[1],\
                       resolution=None, ax=self.ax)

        # set up range
        width  = m1.urcrnrx
        height = m1.urcrnry
        frame = [list(frame[0]), list(frame[1])]
        mapframe = [[-width/2, -height/2], [width/2, height/2]]
        if frame[0][0] != None:
            mapframe[0][0] = width/2*(frame[0][0]-1)
        if frame[0][1] != None:
            mapframe[0][1] = height/2*(frame[0][1]-1)
        if frame[1][0] != None:
            mapframe[1][0] = width/2*(frame[1][0])
        if frame[1][1] != None:
            mapframe[1][1] = height/2*(frame[1][1])

        # set projection
        m = Basemap(projection=projection, lon_0=centre[0], lat_0=centre[1],\
                    llcrnrx=mapframe[0][0], llcrnry=mapframe[0][1],\
                    urcrnrx=mapframe[1][0], urcrnry=mapframe[1][1],\
                    resolution=None, ax=self.ax)

        # turn on 'fulldisk' property when using full disk
        if mapframe == [[-width/2, -height/2], [width/2, height/2]]:
             m._fulldisk = True

        xrange_ = (m.llcrnrx, m.urcrnrx)
        yrange = (m.llcrnry, m.urcrnry)

        # make plot
        for x_vals, y_vals, plot_kwargs, c in\
            itertools.izip(self.x_data_sets, self.y_data_sets,\
                           self.kwarg_sets, default_colors()):

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
                x, y = m(lon, lat)
                lon, lat = m1(x, y, inverse=True)
                if x<=10**20 and y<=10**20\
                and xrange_[0]<x<xrange_[1] and yrange[0]<=y<=yrange[1]:
                    m.ax.text(x, y, r"$%0.0f^\circ$" % lat)

        # label meridians for all projections
        for lon,lat in zip(meridians, [0.]*len(meridians)):
            tlon = (-(lon-180)+180) % 360
            x, y = m(lon, lat)
            lon, lat = m1(x, y, inverse=True)

            if x<=10**20 and y<=10**20\
            and xrange_[0]<x<xrange_[1] and yrange[0]<=y<=yrange[1]:
                m.ax.text(x, y, r"$%0.0f^\circ$" % tlon)

        # set legend
        self.add_legend_if_labels_exist(loc=loc, scatterpoints=1)


###################################################
## unittest section
###################################################

import unittest

# --------------------------------------------
class TestSimpleMapPlot(unittest.TestCase):


    def test_plot(self):
        # set a different dpi for testing purposes
        pylab.rcParams.update({"savefig.dpi": 120})

        # define the SimpleMapPlot
        plot = SimpleMapPlot(r"Range [km]", r"Take off weight [t]", "4 engine planes")

        # define some data.
        # Note: The third item consists of two points (with the same link)
        # while the last data point has no link at all (link=None)
        range = [ [14600], [14800], [13000, 14800], [2800], [4800], [14000]]
        weight = [[365], [560], [374, 442], [46], [392], [50]]
        links = ['http://de.wikipedia.org/wiki/Airbus_A340',\
                'http://de.wikipedia.org/wiki/Airbus_A380',\
                'http://de.wikipedia.org/wiki/Boeing_747',\
                'http://de.wikipedia.org/wiki/Avro_RJ85',\
                'http://de.wikipedia.org/wiki/Antonow_An-124',\
                 None]
        markers = ['ro','bD','yx','cs','g^','k>']

        # set the data to plotutils
        for x,y,link, mark in zip(range, weight, links, markers):
            plot.add_content(x, y, link, color=mark[0], marker=mark[1])

        # finalize the plot
        plot.finalize()

        # save the plot
        plotname = 'plotutils_TestSimpleMapPlot.png'
        plot.savefig(plotname)

        # and here create the html click map and save the html file
        html = plot.create_html(plotname)
        f = file('plotutils_TestSimpleMapPlot.html','w')
        f.write(html)
        f.close()

# --------------------------------------------
if __name__ == '__main__':


    unittest.main()
