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
import re,numpy,math,subprocess,scipy,sys

from glue.ligolw import ligolw,table,lsctables,utils
from glue.ligolw.utils import process as ligolw_process
from glue import segments

from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS
import lal as XLALConstants
from pylal.dq import dqTriggerUtils

from matplotlib import use
use('Agg')
from pylab import hanning

# Hey, scipy, shut up about your nose already.
import warnings
warnings.filterwarnings("ignore")
from scipy import signal as signal
from scipy import interpolate

from matplotlib import mlab
from glue import git_version

__author__  = "Duncan Macleod <duncan.macleod@astro.cf.ac.uk>"
__version__ = "git id %s" % git_version.id
__date__    = git_version.date

"""
This module provides a bank of useful functions for manipulating triggers and trigger files for data quality investigations.
"""

# =============================================================================
# Execute shell command and get output
# =============================================================================

def make_external_call(command):

  """
    Execute shell command and capture standard output and errors. 
    Returns tuple "(stdout,stderr)".
  """

  p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    shell=isinstance(command, str))
  out, err = p.communicate()

  return out,err

# =============================================================================
# Read injection files
# =============================================================================

def frominjectionfile(file, type, ifo=None, start=None, end=None):
  
  """
    Read generic injection file object file containing injections of the given
    type string. Returns an 'Sim' lsctable of the corresponding type.

    Arguments:
   
      file : file object
      type : [ "inspiral" | "burst" | "ringdown" ]

    Keyword arguments:

      ifo : [ "G1" | "H1" | "H2" | "L1" | "V1" ]
  """

  # read type
  type = type.lower()

  # read injection xml
  xml = re.compile('(xml$|xml.gz$)')
  if re.search(xml,file.name):
    xmldoc,digest = utils.load_fileobj(file)
    injtable = table.get_table(xmldoc,'sim_%s:table' % (type))

  # read injection txt
  else:
    cchar = re.compile('[#%<!()_\[\]{}:;\'\"]+')

    #== construct new Sim{Burst,Inspiral,Ringdown}Table
    injtable = lsctables.New(lsctables.__dict__['Sim%sTable' % (type.title())])
    if type=='inspiral':
      columns = ['geocent_end_time.geocent_end_time_ns',\
                 'h_end_time.h_end_time_ns',\
                 'l_end_time.l_end_time_ns',\
                 'v_end_time.v_end_time_ns',\
                 'distance'] 
      for line in file.readlines():
        if re.match(cchar,line):
          continue
        # set up siminspiral object
        inj = lsctables.SimInspiral()
        # split data
        sep = re.compile('[\s,=]+')
        data = sep.split(line)
        # set attributes
        inj.geocent_end_time    = int(data[0].split('.')[0])
        inj.geocent_end_time_ns = int(data[0].split('.')[1])
        inj.h_end_time          = int(data[1].split('.')[0])
        inj.h_end_time_ns       = int(data[1].split('.')[1])
        inj.l_end_time          = int(data[2].split('.')[0])
        inj.l_end_time_ns       = int(data[2].split('.')[1])
        inj.v_end_time          = int(data[3].split('.')[0])
        inj.v_end_time_ns       = int(data[3].split('.')[1])
        inj.distance            = float(data[4])

        injtable.append(inj)

    if type=='burst':
      if file.readlines()[0].startswith('filestart'):
        # if given parsed burst file
        file.seek(0)

        snrcol = { 'G1':23, 'H1':19, 'L1':21, 'V1':25 }

        for line in file.readlines():
          inj = lsctables.SimBurst()
          # split data
          sep = re.compile('[\s,=]+')
          data = sep.split(line)
          # set attributes

          # gps time
          if 'burstgps' in data:
            idx = data.index('burstgps')+1
            geocent = LIGOTimeGPS(data[idx])

            inj.time_geocent_gps    = geocent.seconds
            inj.time_geocent_gps_ns = geocent.nanoseconds
          else:
            continue


          #inj.waveform            = data[4]
          #inj.waveform_number     = int(data[5])

          # frequency
          if 'freq' in data:
            idx = data.index('freq')+1
            inj.frequency = float(data[idx])
          else:
            continue

          # SNR a.k.a. amplitude
          if ifo and 'snr%s' % ifo in data:
            idx = data.index('snr%s' % ifo)+1
            inj.amplitude = float(data[idx])
          elif 'rmsSNR' in data:
            idx = data.index('rmsSNR')+1
            inj.amplitude = float(data[idx])
          else:
            continue

          if 'phi' in data:
            idx = data.index('phi' )+1
            inj.ra = float(data[idx])*24/(2*math.pi)       

          if 'theta' in data:
            idx = data.index('theta' )+1 
            inj.ra = 90-(float(data[idx])*180/math.pi)

          if ifo and 'hrss%s' % ifo in data:
            idx = data.index('hrss%s' % ifo)+1
            inj.hrss = float(data[idx])
          elif 'hrss' in data:
            idx = data.index('hrss')+1
            inj.hrss = float(data[idx])

          # extra columns to be added when I know how
          #inj.q = 0
          #inj.q                   = float(data[11])
          #h_delay = LIGOTimeGPS(data[41])
          #inj.h_peak_time         = inj.time_geocent_gps+h_delay.seconds
          #inj.h_peak_time_ns      = inj.time_geocent_gps_ns+h_delay.nanoseconds
          #l_delay = LIGOTimeGPS(data[43])
          #inj.l_peak_time         = inj.time_geocent_gps+l_delay.seconds
          #inj.l_peak_time_ns      = inj.time_geocent_gps_ns+l_delay.nanoseconds
          #v_delay = LIGOTimeGPS(data[43])
          #inj.v_peak_time         = inj.time_geocent_gps+v_delay.seconds
          #inj.v_peak_time_ns      = inj.time_geocent_gps_ns+v_delay.nanoseconds

          injtable.append(inj)

      else:
        # if given parsed burst file
        file.seek(0)
        for line in file.readlines():
          inj = lsctables.SimBurst()
          # split data
          sep = re.compile('[\s,]+')
          data = sep.split(line)
          # set attributes
          geocent = LIGOTimeGPS(data[0])
          inj.time_geocent_gps    = geocent.seconds
          inj.time_geocent_gps_ns = geocent.nanoseconds

          injtable.append(inj)

  injections = table.new_from_template(injtable)
  if not start:  start = 0
  if not end:    end   = 9999999999
  span = segments.segmentlist([ segments.segment(start, end) ])
  get_time = dqTriggerUtils.def_get_time(injections.tableName)
  injections.extend(inj for inj in injtable if get_time(inj) in span)

  return injections

# =============================================================================
# Calculate band-limited root-mean-square
# =============================================================================

def blrms(data, sampling, average=1, band=None, ripple_db=50, width=2.0,\
          remove_mean=False, return_filter=False, verbose=False):

  """
    This function will calculate the band-limited root-mean-square of the given
    data, using averages of the given length in the given [fmin,fmax] band
    with a kaiser window.

    Options are included to offset the data, and weight frequencies given a 
    dict object of (frequency:weight) pairs.

    Arguments:

      data : numpy.ndarray
        array of data points
      sampling : int
        number of data points per second

    Keyword arguments:

      average : float
        length of rms in seconds
      band : tuple
        [fmin, fmax] for bandpass
      ripple_db : int
        Attenuation in the stop band, in dB
      width : float
        Desired width of the transition from pass to stop, in Hz
      remove_mean : boolean
      verbose : boolean
  """

  nyq = sampling/2

  # verify band variables
  if band==None:
    band=[0,sampling/2]
  fmin = float(band[0])
  fmax = float(band[1])

  if verbose:
    sys.stdout.write("Calculating BLRMS in band %s-%s Hz...\n" % (fmin, fmax))

  #
  # remove mean
  #

  if remove_mean:
    data = data-data.mean()
    if verbose: sys.stdout.write("Data mean removed.\n")

  #
  # Bandpass data
  #
  if return_filter:
    data, filter = bandpass(data, sampling, fmin, fmax, ripple_db=ripple_db,\
                            width=width, return_filter=True, verbose=verbose)
  else:
    data = bandpass(data, sampling, fmin, fmax, ripple_db=ripple_db,\
                    width=width, return_filter=False, verbose=verbose)


  #
  # calculate rms
  #

  # construct output array
  numsamp = int(average*sampling)
  numaverage = numpy.ceil(len(data)/sampling/average)
  output  = numpy.empty(numaverage)

  # loop over averages
  for i in xrange(len(output)):

    # get indices
    idxmin = i*sampling*average
    idxmax = idxmin + numsamp

    # get data chunk
    chunk = data[idxmin:idxmax]

    # get rms
    output[i] = (chunk**2).mean()**(1/2)

  if verbose: sys.stdout.write("RMS calculated for %d averages.\n"\
                               % len(output))

  if return_filter:
    return output, filter
  else:
    return output

# =============================================================================
# Bandpass
# =============================================================================

def bandpass(data, sampling, fmin, fmax, ripple_db=50, width=2.0,\
             return_filter=False, verbose=False):

  """
    This function will bandpass filter data in the given [fmin,fmax] band
    using a kaiser window.

    Arguments:

      data : numpy.ndarray
        array of data points
      sampling : int
        number of data points per second
      fmin : float
        frequency of lowpass
      fmax : float
        frequency of highpass

    Keyword arguments:

      ripple_db : int
        Attenuation in the stop band, in dB
      width : float
        Desired width of the transition from pass to stop, in Hz
      return_filter: boolean
        Return filter
      verbose : boolean
  """

  # construct filter
  order, beta = signal.kaiserord(ripple_db, width*2/sampling)

  lowpass = signal.firwin(order, fmin*2/sampling, window=('kaiser', beta))
  highpass = - signal.firwin(order, fmax*2/sampling, window=('kaiser', beta))
  highpass[order//2] = highpass[order//2] + 1

  bandpass = -(lowpass + highpass); bandpass[order//2] = bandpass[order//2] + 1

  # filter data forward then backward
  data = signal.lfilter(bandpass,1.0,data)
  data = data[::-1]
  data = signal.lfilter(bandpass,1.0,data)
  data = data[::-1]

  if verbose: sys.stdout.write("Bandpass filter applied to data.\n")

  if return_filter:
    return data, bandpass
  else:
    return data

# =============================================================================
# Lowpass
# =============================================================================

def lowpass(data, sampling, fmin, ripple_db=50, width=2.0,\
            return_filter=False, verbose=False):

  """
    This function will lowpass filter data in the given fmin band
    using a kaiser window.

    Arguments:

      data : numpy.ndarray
        array of data points
      sampling : int
        number of data points per second
      fmin : float
        frequency of lowpass

    Keyword arguments:

      ripple_db : int
        Attenuation in the stop band, in dB
      width : float
        Desired width of the transition from pass to stop, in Hz
      return_filter: boolean
        Return filter
      verbose : boolean
  """

  # construct filter
  order, beta = signal.kaiserord(ripple_db, width*2/sampling)

  lowpass = signal.firwin(order, fmin*2/sampling, window=('kaiser', beta))

  # filter data forward then backward
  data = signal.lfilter(lowpass,1.0,data)
  data = data[::-1]
  data = signal.lfilter(lowpass,1.0,data)
  data = data[::-1]

  if verbose: sys.stdout.write("Lowpass filter applied to data.\n")

  if return_filter:
    return data, lowpass
  else:
    return data

# =============================================================================
# Highpass
# =============================================================================

def highpass(data, sampling, fmax, ripple_db=50, width=2.0,\
             return_filter=False, verbose=False):

  """
    This function will highpass filter data in the given fmax band
    using a kaiser window.

    Arguments:

      data : numpy.ndarray
        array of data points
      sampling : int
        number of data points per second
      fmax : float
        frequency of highpass

    Keyword arguments:

      ripple_db : int
        Attenuation in the stop band, in dB
      width : float
        Desired width of the transition from pass to stop, in Hz
      return_filter: boolean
        Return filter
      verbose : boolean
  """

  # construct filter
  order, beta = signal.kaiserord(ripple_db, width*2/sampling)

  highpass = - signal.firwin(order, fmax*2/sampling, window=('kaiser', beta))
  highpass[order//2] = highpass[order//2] + 1

  # filter data forward then backward
  data = signal.lfilter(highpass,1.0,data)
  data = data[::-1]
  data = signal.lfilter(highpass,1.0,data)
  data = data[::-1]

  if verbose: sys.stdout.write("Highpass filter applied to data.\n")

  if return_filter:
    return data, highpass
  else:
    return data

# =============================================================================
# Calculate spectrum
# =============================================================================

def spectrum(data, sampling, NFFT=256, overlap=0.5,\
             window='hanning', detrender=mlab.detrend_linear,\
             sides='onesided', scale='PSD'):

  numpoints  = len(data)
  numoverlap = int(sampling * (1.0 - overlap))

  if isinstance(window,str):
    window=window.lower()

  win = signal.get_window(window, NFFT)

  # calculate PSD with given parameters
  spec,freq = mlab.psd(data, NFFT=NFFT, Fs=sampling, noverlap=numoverlap,\
                       window=win, sides=sides, detrend=detrender)

  # rescale data to meet user's request
  scale = scale.lower()
  if scale == 'asd':
    spec = numpy.sqrt(spec) * numpy.sqrt(2 / (sampling*sum(win**2)))
  elif scale == 'psd':
    spec *= 2/(sampling*sum(win**2))
  elif scale == 'as':
    spec = nump.sqrt(spec) * numpy.sqrt(2) / sum(win)
  elif scale == 'ps':
    spec = spec * 2 / (sum(win)**2)

  return freq, spec.flatten()

# =============================================================================
# Median Mean Spectrum
# =============================================================================

def AverageSpectrumMedianMean(data, fs, NFFT=256, overlap=128,\
                              window=('kaiser',24), sides='onesided',\
                              verbose=False, log=False, warn=True):

  """
    Computes power spectral density of a data series using the median-mean
    average method.
  """

  if sides!='onesided':
    raise NotImplementedError('Only one sided spectrum implemented for the momen')

  # cast data series to numpy array
  data = numpy.asarray(data)

  # number of segments (must be even)
  if overlap==0:
    numseg = int(len(data)/NFFT)
  else:
    numseg = 1 + int((len(data)-NFFT)/overlap)
  assert (numseg - 1)*overlap + NFFT == len(data),\
         "Data is wrong length to be covered completely, please resize"

  # construct window
  win = scipy.signal.get_window(window, NFFT)

  if verbose: sys.stdout.write("%s window constructed.\nConstructing "
                               "median-mean average spectrum "
                               "with %d segments...\n"\
                               % (window, numseg))

  #
  # construct PSD
  #

  # fft scaling factor for units of Hz^-1
  scaling_factor = 1 / (fs * NFFT)

  # construct frequency
  f = numpy.arange(NFFT//2 + 1) * (fs / NFFT)

  odd  = numpy.arange(0, numseg, 2)
  even = numpy.arange(1, numseg, 2)

  # if odd number of segments, ignore the first one (better suggestions welcome)
  if numseg == 1:
    odd = [0]
    even = []
  elif numseg % 2 == 1:
    odd = odd[:-1]
    numseg -= 1
    if warn:
      sys.stderr.write("WARNING: odd number of FFT segments, skipping last.\n")

  # get bias factor
  biasfac = MedianBias(numseg//2)
  # construct normalisation factor
  normfac = 1/(2*biasfac)

  # set data holder
  S = numpy.empty((numseg, len(f)))

  # loop over segments
  for i in xrange(numseg):

    # get data
    chunk = data[i*overlap:i*overlap+NFFT]
    # apply window
    wdata = WindowDataSeries(chunk, win)
    # FFT
    S[i]  = PowerSpectrum(wdata, sides) * scaling_factor

  if verbose: sys.stdout.write("Generated spectrum for each chunk.\n")

  # compute median-mean average
  if numseg > 1:
    S_odd = numpy.median([S[i] for i in odd],0)
    S_even = numpy.median([S[i] for i in even],0)
    S = (S_even  + S_odd) * normfac
  else:
    S = S.flatten()
  if verbose: sys.stdout.write("Calculated median-mean average.\n")

  if log:
    f_log = numpy.logspace(numpy.log10(f[1]), numpy.log10(f[-1]), num=len(f)/2,\
                           endpoint=False)
    I = interpolate.interp1d(f, S)
    S = I(f_log)
    return f_log, S

  return f, S

# =============================================================================
# Median bias factor
# =============================================================================

def MedianBias(nn):

  """
    Returns the median bias factor.
  """

  nmax = 1000;
  ans  = 1;
  n    = (nn - 1)//2;
  if nn >= nmax:
   return numpy.log(2)

  for i in xrange(1, n+1):
    ans -= 1.0/(2*i);
    ans += 1.0/(2*i + 1);

  return ans;

# =============================================================================
# Median average spectrum
# =============================================================================

def AverageSpectrumMedian(data, fs, NFFT=256, overlap=128,\
                          window='hanning', sides='onesided',\
                          verbose=False):

  """
    Construct power spectral density for given data set using the median
    average method.  
  """

  if sides!='onesided':
    raise NotImplementedError('Only one sided spectrum implemented for the momen')

  # cast data series to numpy array
  data = numpy.asarray(data)

  # number of segments (must be even)
  if overlap==0:
    numseg = int(len(data)/NFFT)
  else:
    numseg = 1 + int((len(data)-NFFT)/overlap)
  assert (numseg - 1)*overlap + NFFT == len(data),\
         "Data is wrong length to be covered completely, please resize"

  # construct window
  win = scipy.signal.get_window(window, NFFT)

  if verbose: sys.stdout.write("%s window constructed.\nConstructing "
                               "median average spectrum "
                               "with %d segments...\n"\
                               % (window.title(), numseg))

  #
  # construct PSD
  #

  # fft scaling factor for units of Hz^-1
  scaling_factor = 1 / (fs * NFFT)

  # construct frequency
  f = numpy.arange(NFFT//2 + 1) * (fs / NFFT)

  # get bias factor
  biasfac = MedianBias(numseg)

  # construct normalisation factor
  normfac = 1/(biasfac)

  # set data holder
  S = numpy.empty((numseg, len(f)))

  # loop over segments
  for i in xrange(numseg):

    # get data
    chunk = data[i*overlap:i*overlap+NFFT]
    # apply window
    wdata = WindowDataSeries(chunk, win)
    # FFT
    S[i]  = PowerSpectrum(wdata, sides) * scaling_factor

  if verbose: sys.stdout.write("Generated spectrum for each chunk.\n")

  # compute median-mean average
  if numseg > 1:
    S = scipy.median([S[i] for i in odd])*normfac
  else:
    S = S.flatten()
  if verbose: sys.stdout.write("Calculated median average.\n")

  return f, S 

# =============================================================================
# Apply window
# =============================================================================

def WindowDataSeries(series, window=None):

  """
    Apply window function to data set, defaults to Hanning window.
  """

  # generate default window
  if window == None:
    window = scipy.signal.hanning(len(series))

  # check dimensions
  assert len(series)==len(window), 'Window and data must be same shape'

  # get sum of squares
  sumofsquares = (window**2).sum()
  assert sumofsquares > 0, 'Sum of squares of window non-positive.'

  # generate norm
  norm = (len(window)/sumofsquares)**(1/2)

  # apply window
  return series * window * norm

# =============================================================================
# Power spectrum
# =============================================================================

def PowerSpectrum(series, sides='onesided'):

  """
    Calculate power spectum of given series
  """

  if sides!='onesided':
    raise NotImplementedError('Only one sided spectrum implemented for the moment')

  # apply FFT
  tmp = numpy.fft.fft(series, n=len(series))

  # construct spectrum
  if sides=='onesided':
    spec = numpy.empty(len(tmp)//2+1)
  elif sides=='twosided':
    spec = numpy.empty(len(tmp))

  # DC component
  spec[0] = tmp[0]**2

  # others
  s = (len(series)+1)//2
  spec[1:s] = 2 * (tmp[1:s].real**2 + tmp[1:s].imag**2)

  # Nyquist
  if len(series) % 2 == 0:
    spec[len(series)/2] = tmp[len(series)/2]**2

  return spec

# =============================================================================
# Inspiral range
# =============================================================================

def inspiral_range(f, S, rho=8, m1=1.4, m2=1.4, fmin=30, fmax=4096,\
                   horizon=False):

  """
    Calculate inspiral range for a given spectrum.
  """

  Mpc = 10**6 * XLALConstants.PC_SI

  # compute chirp mass and total mass (in units of solar masses)
  mtot = m1 + m2;
  reducedmass = m1*m2/mtot;
  mchirp = reducedmass**(3/5)*mtot**(2/5);

  # calculate prefactor in m^2
  mchirp *= XLALConstants.MSUN_SI * XLALConstants.G_SI /\
            XLALConstants.C_SI**2
  pre = (5 * XLALConstants.C_SI**(1/3) * mchirp**(5/3) * 1.77**2) /\
        (96 * numpy.pi ** (4/3) * rho**2)

  # include fisco
  fisco = XLALConstants.C_SI**3/XLALConstants.G_SI/XLALConstants.MSUN_SI/\
      6**1.5/numpy.pi/mtot

  # restrict to range, include fisco
  condition = (f >= fmin) & (f < min(fmax,fisco))
  S         = S[condition]
  f         = f[condition]

  # calculate integrand
  integrand = (f**(-7/3))/S

  # integrate
  result = scipy.integrate.trapz(integrand, f)

  R = (pre*result) ** 0.5 / Mpc

  if horizon: R *= 2.26

  return R

# =============================================================================
# Frequency dependent Burst range
# =============================================================================

def f_dependent_burst_range(f, S, rho=8, E=1e-2):
  """
    Calculate GRB-like or supernov-like burst range for a given spectrum    and background trigger SNR at a given time as a function of freqeucy.
  """

  Mpc = 10**6 * XLALConstants.PC_SI

  # generate frequency dependent range
  A = (((XLALConstants.G_SI * (E*XLALConstants.MSUN_SI) * 2/5)/(XLALConstants.PI**2 * XLALConstants.C_SI))**(1/2))/Mpc
  R = A/ (rho * S**(1/2) * f)

  return R

# =============================================================================
# Burst range
# =============================================================================

def burst_range(f, S, rho=8, E=1e-2, fmin=64, fmax=500):
  """
    Calculate GRB-like or supernova-like burst range for a given spectrum
    and background trigger SNR.
  """

  # restrict spectrum to given frequency range
  condition = (f>=fmin) & (f<fmax)
  S2 = S[condition]
  f2 = f[condition]

  # calculate integral
  FOM1 = scipy.integrate.trapz(f_dependent_burst_range(f2, S2, rho, E)**3, f2)
  FOM2 = FOM1/(fmax-fmin)

  return FOM2**(1/3)

def burst_sg_range(f, S, centralFreq, Q, rho=8, E=1e-2, fmin=64, fmax=500):
  """
    Calculate range for sine-Gaussians for a given spectrum
    and background trigger SNR, assuming isotropic GW emission (unphysical but simple)
  """

  # restrict spectrum to given frequency range
  condition = (f>=fmin) & (f<fmax)
  S2 = S[condition]
  f2 = f[condition]

  # generate frequency dependent range
  Mpc = 10**6 * XLALConstants.PC_SI
  1/centralFreq
  A = (((XLALConstants.G_SI * (E*XLALConstants.MSUN_SI) )/(XLALConstants.PI**2 * XLALConstants.C_SI))**(1/2))/Mpc/centralFreq
  sigmaSq = Q**2 / (4 * XLALConstants.PI**2 * centralFreq**2)
  sg = numpy.exp( - (f2 - centralFreq)**2 * sigmaSq / 2) 
  normSG = scipy.integrate.trapz(sg**2, f2)**(1/2)
  sg = sg/normSG
  R = A * sg / (rho * S2**(1/2) )

  # calculate integral
  FOM1 = scipy.integrate.trapz(R**2, f2)

  # factor 0.36, median antenna pattern factor for a single detector
  # and linearly polarized GWs
  return FOM1**(1/2)*0.36
