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
Class definitions and helper functions for building a GW
interferometer or system noise budget.
Bugs/change requests: https://bugs.ligo.org/redmine/projects/detchar/issues/new
"""

# =============================================================================
# Preamble
# =============================================================================

from __future__ import division

import numpy
import re
import copy

from pylal import seriesutils
from pylal import plotutils
from pylal import git_version

# set metadata
__author__  = "Duncan M. Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__    = git_version.date

# =============================================================================
# NoiseBudget class
# =============================================================================

class NoiseTerm(object):
    """
    Object representing one term in a noise budget/projection. It holds a few
    parameters and has methods to read data and generate ASD spectrum borrowed
    from pylal.seriesutils, and simple plotting from pylal.plotutils
    """
    # set types 
    typemap = {\
               'channel':    str,\
               'frame_type': str,\
               'name':       str,\
               'sum':        bool,\
               'data':       numpy.array,\
               'deltaT':     float,\
               'spectrum':   numpy.array,\
               'deltaF':     float,\
               'f0':         float,\
               'epoch':      seriesutils.lal.LIGOTimeGPS,\
               'ref_spectrum': numpy.array,\
               'ref_deltaF': float,\
               'ref_f0':     float,\
               'ref_epoch':  seriesutils.lal.LIGOTimeGPS\
               }
    # set defaults
    name       = None
    channel    = None
    frame_type = None
    sum        = False
    data       = None
    deltaT     = 0
    spectrum   = None
    deltaF     = 0
    f0         = 0
    epoch      = seriesutils.lal.LIGOTimeGPS()
    ref_spectrum = None
    ref_deltaF = 0
    ref_f0     = 0
    ref_epoch  = seriesutils.lal.LIGOTimeGPS()

    def __init__(self, **kwargs):
        """"
    Initialise a NoiseTerm.

    Keyword arguments:

        name : str
            description of this NoiseTerm
        channel : str
            data channel from which this noise term is estimated
        frame_type : str
            LDR frame type for files containing data for this NoiseTerm
        """
        for arg,val in kwargs.items():
            setattr(self, arg, self.typemap[arg](val))

    #
    # define methods from seriesutils
    #

    def fromNDS(self, *args, **kwargs):
        self.timeseries = seriesutils.fromNDS(*args, **kwargs)
        self.data   = self.timeseries.data.data
        self.deltaT = self.timeseries.deltaT
        self.epoch  = self.timeseries.epoch
        self.f0     = self.timeseries.f0
    fromNDS.__doc__ = seriesutils.fromNDS.__doc__

    def fromcache(self, *args, **kwargs):
        self.timeseries = seriesutils.fromlalcache(*args, **kwargs)
        self.data = self.timeseries.data.data
        self.deltaT = self.timeseries.deltaT
        self.epoch  = self.timeseries.epoch
        self.f0     = self.timeseries.f0
    fromcache.__doc__ = seriesutils.fromlalcache.__doc__

    def fromframefile(self, *args, **kwargs):
        self.timeseries = seriesutils.fromframefile(*args, **kwargs)
        self.data = self.timeseries.data.data
        self.deltaT = self.timeseries.deltaT
        self.epoch  = self.timeseries.epoch
        self.f0     = self.timeseries.f0
    fromframefile.__doc__ = seriesutils.fromframefile.__doc__

    def compute_average_spectrum(self, *args, **kwargs):
        self.frequencyseries =\
            seriesutils.compute_average_spectrum(*args, **kwargs)
        self.spectrum = self.frequencyseries.data.data
        self.epoch    = self.frequencyseries.epoch
        self.f0       = self.frequencyseries.f0
        self.deltaF   = self.frequencyseries.deltaF
    compute_average_spectrum.__doc__=\
        seriesutils.compute_average_spectrum.__doc__

    #
    # transfer functions and calibration
    #

    def apply_calibration(self, func):
        """
    Apply a calibration function to the timeseries.

    Arguments:

        func : callable
            any callable function that accepts numpy.array as input
        """
        if not hasattr(self, "data"):
            raise RuntimeError("No data has been loaded for this channel. "+\
                               "Please use one of the data getting methods f"+\
                               "or this NoiseTerm before applying calibration.")
        self.data = func(self.data)

    def apply_spectrum_calibration(self, func):   
        """
    Apply the transfer function to the spectrum.

    Arguments:

        func : [ callable | numpy.array ]
            any callable function that accepts numpy.array as input or
            array by which to multiply spectrum.
        """
        if not hasattr(self, "spectrum"):
            raise RuntimeError("No spectrum has been loaded for this channel."+\
                               " Please use one of the spectrum getting "+\
                               "methods for this NoiseTerm before applying "+\
                               "calibration.")
        if hasattr(func, "shape"):
            self.spectrum *= func
        else:
            self.spectrum = func(self.spectrum)

    def loadreference(self, referencefile, epoch=seriesutils.lal.LIGOTimeGPS(),\
                      sampleUnits=seriesutils.lal.lalStrainUnit, fcol=0,acol=1):
        """
    Read the reference spectrum from a two-column file of frequency and ASD.

    Arguments:

        referencefile : str
            path to file with space-delimited columns of frequency and ASD

    Keyword arguments:

        epoch : lal.LIGOTimeGPS
            GPS epoch of this reference
        sampleUnits : lal.LALUnit
            amplitude unit for reference spectrum
        fcol : int
            number of column holding frequency array, default 0
        acol : int
            number of column holding amplitude array, default 1
        """
        f, data = numpy.loadtxt(referencefile, dtype=float, unpack=True,\
                                usecols=[fcol, acol])
        deltaF  = f[1]-f[0]
        f0      = f[0]
        self.ref_frequencyseries = seriesutils.fromarray(data, str(self.name),
                                                         epoch=epoch,\
                                                         deltaT=deltaF, f0=f0,\
                                                         frequencyseries=True)
        self.ref_spectrum = data
        self.ref_f0       = f0
        self.ref_deltaF   = deltaF

    #
    # restrict to frequency band
    #

    def apply_frequency_band(self, fmin, fmax):
        """
    Restrict the spectrum for this NoiseTerm to the given semiopen
    [fmin, fmax) interval.

    Arguments:

        fmin : float
            minimum frequency for this NoiseTerm
        fmax : float
            maximum frequency for this NoiseTerm
        """
        f = numpy.arange(len(self.spectrum))*self.deltaF + self.f0
        numpy.putmask(self.spectrum, ((f<flim[0])|(f>=flim[1])), 0)

    #
    # plotter
    #

    def plot(self, outfile, **params):
        """
        Plot the spectrum of this NoiseTerm.

        Arguments:

        outfile : str
            path to output file for this plot
        """
        # labels
        xlabel   = params.pop("xlabel", "Frequency (Hz)")
        ylabel   = params.pop("ylabel", "Strain amplitude spectral density "
                                        "$/\sqrt{\mbox{Hz}}$")
        title    = params.pop("title", "%s noise curve"\
                                       % plotutils.display_name(self.name))
        subtitle = params.pop("subtitle", "")

        # extract params
        bbox_inches     = params.pop("bbox_inches", "tight")
        hidden_colorbar = params.pop("hidden_colorbar", False)
        logx   = params.pop("logx", True)
        logy   = params.pop("logy", True)

        # build limits
        minx = self.deltaF
        maxx = self.spectrum.size*self.deltaF + self.f0
        xlim   = params.pop("xlim", [minx, maxx])
        ylim   = params.pop("ylim", None)

        # set plot object
        plot = plotutils.DataPlot(xlabel, ylabel, title, subtitle)

        f = numpy.arange(self.spectrum.size) * self.deltaF +\
                         self.f0
        plot.add_content(f, self.spectrum,\
                         label=plotutils.display_name(self.name),\
                         color=self.color, **params)
       
        if self.ref_spectrum is not None:
            params["linestyle"] = "--"
            params['linewidth'] = 0.3
            f = numpy.arange(self.ref_spectrum.size) *\
                self.ref_frequencyseries.deltaF + self.frequencyseries.f0
            plot.add_content(f, self.ref_spectrum,\
                         label=plotutils.display_name(self.name),\
                         color=self.color, **params)

        # finalize plot
        plot.finalize(logx=logx, logy=logy, loc='upper right',\
                      hidden_colorbar=hidden_colorbar)
        if xlim: plot.ax.set_xlim(xlim)
        if ylim: plot.ax.set_ylim(ylim)

        # save
        plot.savefig(outfile, bbox_inches=bbox_inches,\
                     bbox_extra_artists=plot.ax.texts)

# =============================================================================
# NoiseBudget class
# =============================================================================

class NoiseBudget(list):
    """
    Object representing a list of NoiseTerms comprising a noise budget
    estimation.
    """
    typemap = {\
               'name':      str,\
               'numTerms':  int,\
               'target':    NoiseTerm,\
               'noise_sum': NoiseTerm\
              }

    name      = None
    numTerms  = 0
    target    = None
    noise_sum = None

    # init
    def __init__(self, *args, **kwargs):
        """
    Initialise this NoiseBudget. Arguments should be NoiseTerms to add to the
    NoiseBudget.

    Keyword arguments:

        name : str
            name for this NoiseBudget
        """
        # init
        list.__init__([])

        # set arguments
        for arg in args:
            self.append(NoiseTerm(arg))

        # set keyword arguments
        for arg,val in kwargs.items():
            setattr(self, arg, self.typemap[arg](val))


    # append
    def append(self, *args, **kwargs):
        list.append(self, *args, **kwargs)
        self.numTerms += 1

    def sort(self, *args, **kwargs):
        kwargs.setdefault("key", lambda t: t.name)
        list.sort(self, *args, **kwargs)

    #
    # methods
    #

    # set target
    def set_target(self, noiseterm):
        """
    Set the target of this NoiseBudget to be the noiseterm. Argument should be
    a NoiseTerm object.
        """
        if not isinstance(noiseterm, NoiseTerm):
            raise TypeError("Target must be a NoiseTerm.")
        self.target = noiseterm
    def get_target(self):
        """
    Returns the target of this NoiseBudget. Returns NoneType if no target
    has been set.
        """
        return self.target
 
    # calculate sum of NoiseTerms
    def compute_noise_sum(self, name="Noise sum"):
        """
    Calculate the quadrature sum of terms in the NoiseBudget. All terms whose
    sum attribute evaluates to True will be included.
        """    
        # get metadata from first sum term
        sum_terms = [t for t in self if t.sum]
        if len(sum_terms):
            t = sum_terms[0]
            epoch = t.epoch
            deltaF = t.deltaF
            f0    = t.f0
        else:
            raise RuntimeError("No NoiseTerms were set to be included in the "+\
                               "budget sum.")
        self.noise_sum = NoiseTerm(name=name)
        data = (numpy.asarray([t.spectrum for t in\
                              sum_terms])**2).sum(axis=0)**(1/2)
        self.noise_sum.frequencyseries =\
            seriesutils.fromarray(data, str(name), epoch=epoch,\
                                  deltaT=deltaF, f0=f0, frequencyseries=True)
        self.noise_sum.fromarray

    # compute difference between noise target and noise_sum
    def compute_deficit(self, func=lambda (s,t): abs(1-(t/s)**2)**(1/2)):
        """
    Calculate the deficit of thise NoiseBudget as the normalised quadrature
    difference ratio the sum of the NoiseTerms and the target. Defaults to:

    deficit = $\sqrt{1 - \left(\frac{target}{sum of terms}\right)^2}

    Keyword arguments:

        func : callable
            any callable function that accepts noise sum and target arrays
        """
        if self.target is None:
            raise RuntimeError("NoiseBudget target has not been set. "+\
                               "Run NoiseBudget.set_target(mytarget).")
        if self.noise_sum is None:
            raise RuntimeError("NoiseBudget sum has not be computed. "+\
                               "Run NoiseBudget.compute_noise_sum().")
        deficit        = NoiseTerm('%d budget deficit' % self.name)
        deficit.epoch  = self.target.epoch
        deficit.deltaF = self.target.deltaF
        deficit.f0     = self.target.f0
        data           = func(self.noise_sum.spectrum, self.target.spectrum)
        deficit.frequencyseries = seriesutils.fromarray(data,\
                                                        name=deficit.name,\
                                                        epoch=deficit.epoch,\
                                                        deltaT=deficit.deltaF,\
                                                        f0=deficit.f0,\
                                                        frequencyseries=True)
        return deficit

    # compute difference between noise target and noise_sum
    def compute_ratio_deficit(self):
        """
    Calculate the deficit of thise NoiseBudget as the ratio of the sum of the 
    NoiseTerms and the target.

    ratio_deficit = $\frac{target}{sum of terms}
        """
        return self.compute_deficit(func=lambda (s,t): t/s)

    # plot
    def plot(self, outfile, **params):

        """
    Plot this NoiseBudget
        """

        # extract params
        bbox_inches     = params.pop("bbox_inches", "tight")
        hidden_colorbar = params.pop("hidden_colorbar", False)
        logx   = params.pop("logx", True)
        logy   = params.pop("logy", True)
        
        plot_component_ref = params.pop("plot_component_ref", False)
        if plot_component_ref == 'true':
          plot_component_ref = True

        # build limits
        minx = min(t.deltaF for t in self)
        maxx = max(t.spectrum.size*t.deltaF + t.f0 for t in self)
        xlim   = params.pop("xlim", [minx, maxx])
        ylim   = params.pop("ylim", None)

        # labels
        xlabel   = params.pop("xlabel", "Frequency (Hz)")
        ylabel   = params.pop("ylabel", "Strain amplitude spectral density "
                                        "$/\sqrt{\mbox{Hz}}$")
        title    = params.pop("title", plotutils.display_name(self.name))
        subtitle = params.pop("subtitle", "")

        # copy params for reference lines
        refparams = copy.deepcopy(params)
        refparams.setdefault('linestyle', '--')
        refparams['linewidth'] = 0.3

        # set plot object
        plot = plotutils.DataPlot(xlabel, ylabel, title, subtitle)
        reference_plotted = False

        # plot target and noise sum
        for term in ["target", "noise_sum"]:
            if hasattr(self, term) and getattr(self, term) is not None:
                term = getattr(self, term)
                f = numpy.arange(term.data.size) *term.frequencyseries.deltaF +\
                    term.frequencyseries.f0

                plot.add_content(f, term.data,\
                                 label=plotutils.display_name(term.name),\
                                 linecolor=term.color, **params)
 
                # plot reference
                if term.ref_spectrum is not None:
                    f = numpy.arange(term.ref_spectrum.size) *\
                        self.ref_frequencyseries.deltaF +\
                        self.ref_frequencyeries.f0
                    plot.add_content(f, term.ref_spectrum,\
                                     label=plotutils.display_name(self.name),\
                                     color=self.color, **refparams)
                    reference_plotted = True

        # plot terms in the budget
        for term in self:
            f = numpy.arange(term.spectrum.size) * term.deltaF +\
                term.f0
            plot.add_content(f, term.spectrum,\
                             label=plotutils.display_name(term.name),\
                             color=term.color, **params)
            if plot_component_ref and term.ref_spectrum is not None:
                    f = numpy.arange(term.ref_spectrum.size) *\
                        term.ref_frequencyseries.deltaF +\
                        term.ref_frequencyseries.f0
                    plot.add_content(f, term.ref_spectrum,\
                                     label=None, color=term.color, **refparams)
                    reference_plotted = True

        # finalize plot
        plot.finalize(logx=logx, logy=logy, loc='upper right',\
                      hidden_colorbar=hidden_colorbar)
        if xlim: plot.ax.set_xlim(xlim)
        if ylim: plot.ax.set_ylim(ylim)

        # add text box
        if reference_plotted:
            plot.ax.text(0.005, 0.99, 'Dashed line indicates reference',\
                         horizontalalignment='left', verticalalignment='top',\
                         transform=plot.ax.transAxes,\
                         backgroundcolor='white',\
                         bbox=dict(facecolor='white', alpha=0.6,\
                                   edgecolor='black'))

        # save
        plot.savefig(outfile, bbox_inches=bbox_inches,\
                     bbox_extra_artists=plot.ax.texts)

    def plot_deficit(self, outfile, **params):
        """
    Plot the difference between the noise budget sum and its target.
        """
        deficit = self.compute_deficit()
        deficit.plot(outfile, **params)

    def plot_ratio_deficit(self, outfile, **params):
        """
    Plot the difference between the noise budget sum and its target.
        """
        deficit = compute_deficit_ratio()
        deficit.plot(outfile, **params)
