# -*- coding: utf-8 -*-
#
#       pulsarpputils.py
#
#       Copyright 2012
#       Matthew Pitkin <matthew.pitkin@ligo.org>
#
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

# known pulsar analysis post-processing utilities

# Many functions in this a taken from, or derived from equivalents available in
# the PRESTO pulsar software package http://www.cv.nrao.edu/~sransom/presto/

import sys
import math
import os
import numpy as np

import matplotlib
#matplotlib.use("Agg")

from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib.mlab import specgram, find, psd
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d

# pylal stuff
from pylal import date
from pylal.xlal import inject
from pylal.xlal import tools
from pylal.xlal.datatypes.ligotimegps import LIGOTimeGPS

from types import StringType, FloatType

# some common constants taken from psr_constants.py in PRESTO
ARCSECTORAD = float('4.8481368110953599358991410235794797595635330237270e-6')
RADTOARCSEC = float('206264.80624709635515647335733077861319665970087963')
SECTORAD    = float('7.2722052166430399038487115353692196393452995355905e-5')
RADTOSEC    = float('13750.987083139757010431557155385240879777313391975')
RADTODEG    = float('57.295779513082320876798154814105170332405472466564')
DEGTORAD    = float('1.7453292519943295769236907684886127134428718885417e-2')
RADTOHRS    = float('3.8197186342054880584532103209403446888270314977710')
HRSTORAD    = float('2.6179938779914943653855361527329190701643078328126e-1')
PI          = float('3.1415926535897932384626433832795028841971693993751')
TWOPI       = float('6.2831853071795864769252867665590057683943387987502')
PIBYTWO     = float('1.5707963267948966192313216916397514420985846996876')
SECPERDAY   = float('86400.0')
SECPERJULYR = float('31557600.0')
KMPERPC     = float('3.0856776e13')
KMPERKPC    = float('3.0856776e16')
Tsun        = float('4.925490947e-6') # sec
Msun        = float('1.9891e30')      # kg
Mjup        = float('1.8987e27')      # kg
Rsun        = float('6.9551e8')       # m
Rearth      = float('6.378e6')        # m
SOL         = float('299792458.0')    # m/s
MSUN        = float('1.989e+30')      # kg
G           = float('6.673e-11')      # m^3/s^2/kg 
C           = SOL
KPC         = float('3.0856776e19')   # kiloparsec in metres
I38         = float('1e38')           # moment of inertia kg m^2 

# some angle conversion functions taken from psr_utils.py in PRESTO
def rad_to_dms(rad):
  """
  rad_to_dms(rad):
     Convert radians to degrees, minutes, and seconds of arc.
  """
  if (rad < 0.0): sign = -1
  else: sign = 1
  arc = RADTODEG * np.fmod(np.fabs(rad), math.pi)
  d = int(arc)
  arc = (arc - d) * 60.0
  m = int(arc)
  s = (arc - m) * 60.0
  if sign==-1 and d==0:
    return (sign * d, sign * m, sign * s)
  else:
    return (sign * d, m, s)

def dms_to_rad(deg, min, sec):
  """
  dms_to_rad(deg, min, sec):
     Convert degrees, minutes, and seconds of arc to radians.
  """
  if (deg < 0.0):
    sign = -1
  elif (deg==0.0 and (min < 0.0 or sec < 0.0)):
    sign = -1
  else:
    sign = 1
  return sign * ARCSECTORAD * \
    (60.0 * (60.0 * np.fabs(deg) + np.fabs(min)) + np.fabs(sec))

def dms_to_deg(deg, min, sec):
  """
  dms_to_deg(deg, min, sec):
     Convert degrees, minutes, and seconds of arc to degrees.
  """
  return RADTODEG * dms_to_rad(deg, min, sec)

def rad_to_hms(rad):
  """
  rad_to_hms(rad):
     Convert radians to hours, minutes, and seconds of arc.
  """
  rad = np.fmod(rad, 2.*math.pi)
  if (rad < 0.0): rad = rad + 2.*math.pi
  arc = RADTOHRS * rad
  h = int(arc)
  arc = (arc - h) * 60.0
  m = int(arc)
  s = (arc - m) * 60.0
  return (h, m, s)
  
def hms_to_rad(hour, min, sec):
  """
  hms_to_rad(hour, min, sec):
     Convert hours, minutes, and seconds of arc to radians
  """
  if (hour < 0.0): sign = -1
  else: sign = 1
  return sign * SECTORAD * \
         (60.0 * (60.0 * np.fabs(hour) + np.fabs(min)) + np.fabs(sec))

def coord_to_string(h_or_d, m, s):
  """
  coord_to_string(h_or_d, m, s):
     Return a formatted string of RA or DEC values as
     'hh:mm:ss.ssss' if RA, or 'dd:mm:ss.ssss' if DEC.
  """
  retstr = ""
  if h_or_d < 0:
    retstr = "-"
  elif abs(h_or_d)==0:
    if (m < 0.0) or (s < 0.0):
      retstr = "-"
  
  h_or_d, m, s = abs(h_or_d), abs(m), abs(s)
  if (s >= 9.9995):
    return retstr+"%.2d:%.2d:%.4f" % (h_or_d, m, s)
  else:
    return retstr+"%.2d:%.2d:0%.4f" % (h_or_d, m, s)

def ra_to_rad(ra_string):
  """
  ra_to_rad(ar_string):
     Given a string containing RA information as
     'hh:mm:ss.ssss', return the equivalent decimal
     radians. Also deal with cases where input 
     string is just hh:mm, or hh.
  """
  hms = ra_string.split(":")
  if len(hms) == 3:
    return hms_to_rad(int(hms[0]), int(hms[1]), float(hms[2]))
  elif len(hms) == 2:
    return hms_to_rad(int(hms[0]), int(hms[1]), 0.0)
  elif len(hms) == 1:
    return hms_to_rad(float(hms[0]), 0.0, 0.0)
  else:
    print >> sys.stderr, "Problem parsing RA string %s" % ra_string
    sys.exit(1)

def dec_to_rad(dec_string):
  """
  dec_to_rad(dec_string):
     Given a string containing DEC information as
     'dd:mm:ss.ssss', return the equivalent decimal
     radians. Also deal with cases where input string 
     is just dd:mm or dd
  """
  dms = dec_string.split(":")
  if "-" in dms[0] and float(dms[0]) == 0.0:
    m = '-'
  else:
    m = ''
  
  if len(dms) == 3:
    return dms_to_rad(int(dms[0]), int(m+dms[1]), float(m+dms[2]))
  elif len(dms) == 2:
    return dms_to_rad(int(dms[0]), int(m+dms[1]), 0.0)
  elif len(dms) == 1:
    return dms_to_rad(float(dms[0]), 0.0, 0.0)
  else:
    print >> sys.stderr, "Problem parsing DEC string %s" % dec_string
    sys.exit(1)

def p_to_f(p, pd, pdd=None):
  """
  p_to_f(p, pd, pdd=None):
    Convert period, period derivative and period second
    derivative to the equivalent frequency counterparts.
    Will also convert from f to p.
  """
  f = 1.0 / p
  fd = -pd / (p * p)
  if (pdd==None):
    return [f, fd]
  else:
    if (pdd==0.0):
      fdd = 0.0
    else:
      fdd = 2.0 * pd * pd / (p**3.0) - pdd / (p * p)
     
    return [f, fd, fdd]

def pferrs(porf, porferr, pdorfd=None, pdorfderr=None):
  """
  pferrs(porf, porferr, pdorfd=None, pdorfderr=None):
     Calculate the period or frequency errors and
     the pdot or fdot errors from the opposite one.
  """
  if (pdorfd==None):
    return [1.0 / porf, porferr / porf**2.0]
  else:
    forperr = porferr / porf**2.0
    fdorpderr = np.sqrt((4.0 * pdorfd**2.0 * porferr**2.0) / porf**6.0 +
                          pdorfderr**2.0 / porf**4.0)
    [forp, fdorpd] = p_to_f(porf, pdorfd)
  
    return [forp, forperr, fdorpd, fdorpderr]

# class to read in a pulsar par file - this is heavily based on the function
# in parfile.py in PRESTO
float_keys = ["F", "F0", "F1", "F2", "F3", "F4", "F5", "F6",
              "P", "P0", "P1", "P2", "P3", "P4", "P5", "P6",
              "PEPOCH", "POSEPOCH", "DM", "START", "FINISH", "NTOA",
              "TRES", "TZRMJD", "TZRFRQ", "TZRSITE", "NITS",
              "A1", "XDOT", "E", "ECC", "EDOT", "T0", "PB", "PBDOT", "OM",
              "OMDOT", "EPS1", "EPS2", "EPS1DOT", "EPS2DOT", "TASC", "LAMBDA",
              "BETA", "RA_RAD", "DEC_RAD", "GAMMA", "SINI", "M2", "MTOT",
              "FB0", "FB1", "FB2", "ELAT", "ELONG", "PMRA", "PMDEC", "DIST",
              # GW PARAMETERS
              "H0", "COSIOTA", "PSI", "PHI0", "THETA", "I21", "I31"]
str_keys = ["FILE", "PSR", "PSRJ", "NAME", "RAJ", "DECJ", "RA", "DEC", "EPHEM",
            "CLK", "BINARY", "UNITS"]

class psr_par:
  def __init__(self, parfilenm):
    self.FILE = parfilenm
    pf = open(parfilenm)
    for line in pf.readlines():
      # ignore empty lines (i.e. containing only whitespace)
      if not line.strip():
        continue
      
      # Convert any 'D-' or 'D+' to 'E-' or 'E+'
      line = line.replace("D-", "E-")
      line = line.replace("D+", "E+")
      # also check for lower case
      line = line.replace("d-", "e-")
      line = line.replace("d+", "e+")
      splitline = line.split()
      
      # get all upper case version in case lower case in par file
      key = splitline[0].upper()
      
      if key in str_keys:
        setattr(self, key, splitline[1])
      elif key in float_keys:
        try:
          setattr(self, key, float(splitline[1]))
        except:
          continue
      
      if len(splitline)==3: # Some parfiles don't have flags, but do have errors
        if splitline[2] not in ['0', '1']:
          setattr(self, key+'_ERR', float(splitline[2]))
      
      if len(splitline)==4:
        setattr(self, key+'_ERR', float(splitline[3]))
     
    # sky position
    if hasattr(self, 'RAJ'):
      setattr(self, 'RA_RAD', ra_to_rad(self.RAJ))
      
      # set RA error in rads (rather than secs)
      if hasattr(self, 'RAJ_ERR'):
        setattr(self, 'RA_RAD_ERR', hms_to_rad(0, 0, self.RAJ_ERR))
    if hasattr(self, 'DECJ'):
      setattr(self, 'DEC_RAD', dec_to_rad(self.DECJ))
      
      # set DEC error in rads (rather than arcsecs)
      if hasattr(self, 'DECJ_ERR'):
        setattr(self, 'DEC_RAD_ERR', dms_to_rad(0, 0, self.DECJ_ERR))
      
    # periods and frequencies
    if hasattr(self, 'P'):
      setattr(self, 'P0', self.P)
    if hasattr(self, 'P0'):
      setattr(self, 'F0', 1.0/self.P0)
    if hasattr(self, 'F0'):
      setattr(self, 'P0', 1.0/self.F0)
    if hasattr(self, 'FB0'):
      setattr(self, 'PB', (1.0/self.FB0)/86400.0)
    if hasattr(self, 'P0_ERR'):
      if hasattr(self, 'P1_ERR'):
        f, ferr, fd, fderr = pferrs(self.P0, self.P0_ERR,
                                       self.P1, self.P1_ERR)
        setattr(self, 'F0_ERR', ferr) 
        setattr(self, 'F1', fd) 
        setattr(self, 'F1_ERR', fderr) 
      else:
        if hasattr(self, 'P1'):
          f, fd, = p_to_f(self.P0, self.P1)
          setattr(self, 'F0_ERR', self.P0_ERR/(self.P0*self.P0))
          setattr(self, 'F1', fd) 
    if hasattr(self, 'F0_ERR'):
      if hasattr(self, 'F1_ERR'):
        p, perr, pd, pderr = pferrs(self.F0, self.F0_ERR,
                                    self.F1, self.F1_ERR)
        setattr(self, 'P0_ERR', perr) 
        setattr(self, 'P1', pd) 
        setattr(self, 'P1_ERR', pderr) 
      else:
        if hasattr(self, 'F1'):
          p, pd, = p_to_f(self.F0, self.F1)
          setattr(self, 'P0_ERR', self.F0_ERR/(self.F0*self.F0))
          setattr(self, 'P1', pd) 
    
    # binary parameters
    if hasattr(self, 'EPS1') and hasattr(self, 'EPS2'):
      ecc = math.sqrt(self.EPS1 * self.EPS1 + self.EPS2 * self.EPS2)
      omega = math.atan2(self.EPS1, self.EPS2)
      setattr(self, 'E', ecc)
      setattr(self, 'OM', omega * RADTODEG)
      setattr(self, 'T0', self.TASC + self.PB * omega/TWOPI)
    if hasattr(self, 'PB') and hasattr(self, 'A1') and not \
       (hasattr(self, 'E') or hasattr(self, 'ECC')):
      setattr(self, 'E', 0.0)  
    if hasattr(self, 'T0') and not hasattr(self, 'TASC') and hasattr(self, 'OM') and hasattr(self, 'PB'):
      setattr(self, 'TASC', self.T0 - self.PB * self.OM/360.0)
        
    pf.close()
  
  def __getitem__(self, key):
    try:
      par = getattr(self, key)
    except:
      par = None
      
    return par
  
  def __str__(self):
    out = ""
    for k, v in self.__dict__.items():
      if k[:2]!="__":
        if type(self.__dict__[k]) is StringType:
          out += "%10s = '%s'\n" % (k, v)
        else:
          out += "%10s = %-20.15g\n" % (k, v)
      
    return out

# class to read in a nested sampling prior file   
class psr_prior:
  def __init__(self, priorfilenm):
    self.FILE = priorfilenm
    pf = open(priorfilenm)
    for line in pf.readlines():
      splitline = line.split()
      
      # get all upper case version in case lower case in par file
      key = splitline[0].upper()
      
      if key in str_keys:
        # everything in a prior files should be numeric
        setattr(self, key, [float(splitline[1]), float(splitline[2])])
      elif key in float_keys:
        setattr(self, key, [float(splitline[1]), float(splitline[2])])        
 
    # get sky positions in rads as strings 'dd/hh:mm:ss.s'  
    if hasattr(self, 'RA'):
      hl, ml, sl = rad_to_hms(self.RA[0])
      rastrl = coord_to_string(hl, ml, sl)
      hu, mu, su = rad_to_hms(self.RA[1])
      rastru = coord_to_string(hu, mu, su)
      setattr(self, 'RA_STR', [rastrl, rastru])
      
    if hasattr(self, 'DEC'):
      dl, ml, sl = rad_to_dms(self.DEC[0])
      decstrl = coord_to_string(dl, ml, sl)
      du, mu, su = rad_to_dms(self.DEC[1])
      decstru = coord_to_string(du, mu, su)
      setattr(self, 'DEC_STR', [decstrl, decstru])
    
    pf.close()
  
  def __getitem__(self, key):
    try:
      atr = getattr(self, key)
    except:
      atr = None
    
    return atr
  
  def __str__(self):
    out = ""
    for k, v in self.__dict__.items():
      if k[:2]!="__":
        if type(self.__dict__[k]) is StringType:
          out += "%10s = '%s'\n" % (k, v)
        else:
          out += "%10s = %-20.15g, %-20.15g\n" % (k, float(v[0]), float(v[1]))
      
    return out

# Function to return a pulsar's strain spin-down limit given its spin frequency
#(Hz), spin-down (Hz/s) and distance (kpc). The canonical value of moment of
# inertia of 1e38 kg m^2 is used
def spin_down_limit(freq, fdot, dist):
  hsd = math.sqrt((5./2.)*(G/C**3)*I38*math.fabs(fdot)/freq)/(dist*KPC)
  
  return hsd
  
# Function to convert a pulsar stain into ellipticity assuming the canonical
# moment of inertia
def h0_to_ellipticity(h0, freq, dist):
  ell = h0*C**4*dist*KPC/(16.*math.pi**2*G*I38*freq**2)
  
  return ell

# function to convert the psi' and phi0' coordinates used in nested sampling
# into the standard psi and phi0 coordinates (using vectors of those parameters
def phipsiconvert(phipchain, psipchain):
  chainlen=len(phipchain)
  
  phichain = []
  psichain = []
  
  theta = math.atan2(1,2);
  ct = math.cos(theta);
  st = math.sin(theta);
  
  for i in range(0,chainlen):  
    phi0 = (1/(2*st))*phipchain[i] - (1/(2*st))*psipchain[i];
    psi = (1/(2*ct))*phipchain[i] + (1/(2*ct))*psipchain[i];
    
    # put psi between +/-pi/4
    if math.fabs(psi) > math.pi/4.:
      # shift phi0 by pi
      phi0 = phi0 + math.pi;
      
      # wrap around psi
      if psi > math.pi/4.:
        psi = -(math.pi/4.) + math.fmod(psi+(math.pi/4.), math.pi/2.);
      else:
        psi = (math.pi/4.) - math.fmod((math.pi/4.)-psi, math.pi/2.);
  
    # get phi0 into 0 -> 2pi range
    if phi0 > 2.*math.pi: 
      phi0 = math.fmod(phi0, 2.*math.pi);
    else:
      phi0 = 2.*math.pi - math.fmod(2.*math.pi-phi0, 2.*math.pi);
    
    phichain.append(phi0)
    psichain.append(psi)

  return phichain, psichain
 
# function to create histogram plot of the 1D posterior (potentially for
# multiple IFOs) for a parameter (param). If an upper limit is given then
# that will be output
def plot_posterior_hist(poslist, param, ifos,
                        parambounds=[float("-inf"), float("inf")], 
                        nbins=50, upperlimit=0, overplot=False,
                        parfile=None, mplparams=False):
  # create list of figures
  myfigs = []
  
  # create a list of upper limits
  ulvals = []
  
  # set some matplotlib defaults for hist
  if not mplparams:
    mplparams = { \
      'backend': 'Agg',
      'text.usetex': True, # use LaTeX for all text
      'axes.linewidth': 0.5, # set axes linewidths to 0.5
      'axes.grid': True, # add a grid
      'grid.linwidth': 0.5,
      'font.family': 'serif',
      'font.size': 12 }
  
  matplotlib.rcParams.update(mplparams)
  
  # ifos line colour specs
  coldict = {'H1': 'r', 'H2': 'c', 'L1': 'g', 'V1': 'b', 'G1': 'm', 'Joint':'k'}
  
  # some parameter names for special LaTeX treatment in figures
  paramdict = {'H0': '$h_0$', 'COSIOTA': '$\cos{\iota}$', \
               'PSI': '$\psi$ (rad)', 'PHI0': '$\phi_0$ (rad)', \
               'RA': '$\\alpha$ (rad)', 'DEC': '$\delta$ (rad)', \
               'F0': '$f_0$ (Hz)', 'F1': '$\dot{f}$ (Hz/s)', 'F2':
               '$\\\\ddot{f}$ (Hz/s$^2$)', 'LOGL': '$\log{L}$', \
               'PMRA': 'proper motion $\\alpha$ (rad/s)', \
               'PMDEC': 'proper motion $\delta$ (rad/s)'}
               
  # param name for axis label
  try:
    paraxis = paramdict[param.upper()]
  except:
    paraxis = param
  
  ymax = []
  
  # if a par file object is given expect that we have an injection file
  # containing the injected value to overplot on the histgram
  parval = None
  if parfile:
    parval = parfile[param.upper()]
  
  # loop over ifos
  for idx, ifo in enumerate(ifos):
    # check whether to plot all figures on top of each other
    if overplot and idx == 0:
      myfig = plt.figure(figsize=(4,4),dpi=200)
      plt.hold(True)
    elif not overplot:
      myfig = plt.figure(figsize=(4,4),dpi=200)
    
    pos = poslist[idx]
    
    # check for cosiota
    if 'cosiota' in param:
      pos_samps = np.cos(pos['iota'].samples)
    else:
      pos_samps = pos[param].samples
 
    # get a normalised histogram for each
    n, bins = hist_norm_bounds( pos_samps, int(nbins), parambounds[0], \
                                parambounds[1] )
    
    # plot histogram
    #plt.plot(bins, n, color=coldict[ifo])
    plt.step(bins, n, color=coldict[ifo])
    
    if 'h0' not in param:
      plt.xlim(parambounds[0], parambounds[1])
   
    plt.xlabel(r''+paraxis, fontsize=14, fontweight=100)
    plt.ylabel(r'Probability Density', fontsize=14, fontweight=100)
    myfig.subplots_adjust(left=0.18, bottom=0.15) # adjust size
   
    if not overplot:
      plt.ylim(0, n.max()+0.1*n.max()) 
      #plt.legend(ifo)
      # set background colour of axes
      ax = plt.gca()
      ax.set_axis_bgcolor("#F2F1F0")
      myfigs.append(myfig)
    else:
      ymax.append(n.max()+0.1*n.max())
    
    # if upper limit is needed then integrate posterior using trapezium rule
    if upperlimit != 0:
      ct = cumtrapz(n, bins)
      
      # prepend a zero to ct
      ct = np.insert(ct, 0, 0)
      
      # use spline interpolation to find the value at 'upper limit'
      ctu, ui = np.unique(ct, return_index=True)
      intf = interp1d(ctu, bins[ui], kind='cubic')
      ulvals.append(intf(float(upperlimit)))
  
  # plot parameter values
  if parval:
    if not overplot:
      plt.hold(True)
      plt.plot([parval, parval], [0, n.max()+0.1*n.max()], 'k--', linewidth=1.5)
    else:
      plt.plot([parval, parval], [0, max(ymax)], 'k--', linewidth=1.5)
  
  if overplot:
    plt.ylim(0, max(ymax))
    #plt.legend(ifos)
    ax = plt.gca()
    ax.set_axis_bgcolor("#F2F1F0")
    myfigs.append(myfig)
  
  return myfigs, ulvals

# function to plot a posterior chain (be it MCMC chains or nested samples)
# the input should be a list of posteriors for each IFO, and the parameter
# required, the list of IFO. grr is a list of dictionaries giving
# the Gelman-Rubins statistics for the given parameter for each IFO.
# If withhist is set then it will also output a histgram, with withhist number
# of bins
def plot_posterior_chain(poslist, param, ifos, grr=None, withhist=0, \
                         mplparams=False):
  try:
    from matplotlib import gridspec
  except:
    return None
  
  if not mplparams:
    mplparams = { \
      'backend': 'Agg',
      'text.usetex': True, # use LaTeX for all text
      'axes.linewidth': 0.5, # set axes linewidths to 0.5
      'axes.grid': True, # add a grid
      'grid.linwidth': 0.5,
      'font.family': 'serif',
      'font.size': 14 }
  
  matplotlib.rcParams.update(mplparams)
  
  # some parameter names for special LaTeX treatment in figures
  paramdict = {'H0': '$h_0$', 'COSIOTA': '$\cos{\iota}$', \
               'PSI': '$\psi$ (rad)', 'PHI0': '$\phi_0$ (rad)', \
               'RA': '$\\alpha$ (rad)', 'DEC': '$\delta$ (rad)', \
               'F0': '$f_0$ (Hz)', 'F1': '$\dot{f}$ (Hz/s)', 'F2':
               '$\\\\ddot{f}$ (Hz/s$^2$)', 'LOGL': '$\log{L}$', \
               'PMRA': 'proper motion $\\alpha$ (rad/s)', \
               'PMDC': 'proper motion $\delta$ (rad/s)'}
 
  coldict = {'H1': 'r', 'H2': 'c', 'L1': 'g', 'V1': 'b', 'G1': 'm', \
             'Joint': 'k'}
  
  # param name for axis label    
  try:
    if param == 'iota':
      p = 'cosiota'
    else:
      p = param
    
    paryaxis = paramdict[p.upper()]
  except:
    paryaxis = param
  
  if grr:
    legendvals = []
  
  maxiter = 0
  maxn = 0
  minsamp = float('inf')
  maxsamp = -float('inf')
  
  for idx, ifo in enumerate(ifos):
    if idx == 0:
      myfig = plt.figure(figsize=(12,4),dpi=200)
      myfig.subplots_adjust(bottom=0.15)
      
      if withhist:
        gs = gridspec.GridSpec(1,4, wspace=0)
        ax1 = plt.subplot(gs[:-1])
        ax2 = plt.subplot(gs[-1])
        
    pos = poslist[idx]
    
    # check for cosiota
    if 'iota' == param:
      pos_samps = np.cos(pos['iota'].samples)
    else:
      pos_samps = pos[param].samples
    
    if np.min(pos_samps) < minsamp:
      minsamp = np.min(pos_samps)
    if np.max(pos_samps) > maxsamp:
      maxsamp = np.max(pos_samps)
    
    if withhist:
      ax1.hold(True)
      ax1.plot(pos_samps, '.', color=coldict[ifo], markersize=1)

      n, binedges = np.histogram( pos_samps, withhist )
      n = np.append(n, 0)
      ax2.hold(True)
      ax2.step(n, binedges, color=coldict[ifo])
      
      if np.max(n) > maxn:
        maxn = np.max(n)
    else:
      plt.plot(pos_samps, '.', color=coldict[ifo], markersize=1)
      plt.hold(True)

    if grr:
      if 'iota' == param:
        p = 'cosiota'
      else:
        p = param
      
      try:
        legendvals.append(r'$R = %.2f$' % grr[idx][p])
      except:
        legendval = []
    
    if len(pos_samps) > maxiter:
      maxiter = len(pos_samps)
  
  if not withhist:
    ax1 = plt.gca()
  
  bounds = [minsamp, maxsamp]
  
  ax1.set_ylabel(r''+paryaxis, fontsize=16, fontweight=100)
  ax1.set_xlabel(r'Iterations', fontsize=16, fontweight=100)
  
  ax1.set_xlim(0, maxiter)
  ax1.set_ylim(bounds[0], bounds[1])
  
  if withhist:
    ax2.set_ylim(bounds[0], bounds[1])
    ax2.set_xlim(0, maxn+0.1*maxn)
    ax2.set_xlabel(r'Count', fontsize=16, fontweight=100)
    ax2.set_yticklabels([])
    ax2.set_axis_bgcolor("#F2F1F0")
  
  # add gelman-rubins stat data
  if legendvals:
    ax1.legend(legendvals, title='Gelman-Rubins test')
  
  return myfig
  
# function to create a histogram plot of the 2D posterior
def plot_posterior_hist2D(poslist, params, ifos, bounds=None, nbins=[50,50], \
                          parfile=None, mplparams=False):
  if len(params) != 2:
    print >> sys.stderr, "Require 2 parameters"
    sys.exit(1)
  
  # set some matplotlib defaults for amplitude spectral density
  if not mplparams:
    mplparams = { \
      'backend': 'Agg',
      'text.usetex': True, # use LaTeX for all text
      'axes.linewidth': 0.5, # set axes linewidths to 0.5
      'axes.grid': True, # add a grid
      'grid.linwidth': 0.5,
      'font.family': 'serif',
      'font.size': 12 }
  
  matplotlib.rcParams.update(mplparams)
  
  # some parameter names for special LaTeX treatment in figures
  paramdict = {'H0': '$h_0$', 'COSIOTA': '$\cos{\iota}$', \
               'PSI': '$\psi$ (rad)', 'PHI0': '$\phi_0$ (rad)', \
               'RA': '$\\alpha$ (rad)', 'DEC': '$\delta$ (rad)', \
               'F0': '$f_0$ (Hz)', 'F1': '$\dot{f}$ (Hz/s)', 'F2':
               '$\\\\ddot{f}$ (Hz/s$^2$)', 'LOGL': '$\log{L}$', \
               'PMRA': 'proper motion $\\alpha$ (rad/s)', \
               'PMDC': 'proper motion $\delta$ (rad/s)'}
  
  myfigs = []
  
  # param name for axis label
  try:
    parxaxis = paramdict[params[0].upper()]
  except:
    parxaxis = params[0]
    
  try:
    paryaxis = paramdict[params[1].upper()]
  except:
    paryaxis = params[1]
  
  parval1 = None
  parval2 = None
  if parfile:
    parval1 = parfile[params[0].upper]
    parval2 = parfile[params[1].upper]
  
  for idx, ifo in enumerate(ifos):
    posterior = poslist[idx]
    
    if 'cosiota' in params[0]:
      a = np.squeeze(posterior['iota'].samples)
      a = np.cos(a)
      b = np.squeeze(posterior[params[1]].samples)
    elif 'cosiota' in params[1]:
      b = np.squeeze(posterior['iota'].samples)
      b = np.cos(b)
      a = np.squeeze(posterior[params[0]].samples)
    else:
      a = np.squeeze(posterior[params[0]].samples)
      b = np.squeeze(posterior[params[1]].samples)
    
    # Create 2D bin array
    par1pos_min = a.min()
    par2pos_min = b.min()

    par1pos_max = a.max()
    par2pos_max = b.max()

    myfig = plt.figure(figsize=(4,4),dpi=200)
   
    plt.xlabel(r''+parxaxis, fontsize=14, fontweight=100)
    plt.ylabel(r''+paryaxis, fontsize=14, fontweight=100, rotation=270)
    
    H, xedges, yedges = np.histogram2d(a, b, nbins, normed=True)

    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    plt.imshow(np.transpose(H), aspect='auto', extent=extent, \
interpolation='bicubic', cmap='gray_r')
    #plt.colorbar()
    if bounds:
      plt.xlim(bounds[0][0], bounds[0][1])
      plt.ylim(bounds[1][0], bounds[1][1])
    
    # plot injection values if given
    if parval1 and parval2:
      plt.hold(True)
      plt.plot(parval1, parval2, 'rx', markersize=2)
    
    myfig.subplots_adjust(left=0.18, bottom=0.15) # adjust size
    
    myfigs.append(myfig)
    
  return myfigs
  
# a function that creates and normalises a histograms of samples, with nbins
# between an upper and lower bound a upper and lower bound. The values at the
# bin points (with length nbins+2) will be returned as numpy arrays
def hist_norm_bounds(samples, nbins, low=float("-inf"), high=float("inf")):
  # get histogram
  n, binedges = np.histogram( samples, nbins )
  
  # get bin width
  binwidth = binedges[1] - binedges[0]
  
  # create bin centres
  bincentres = np.array([])
  for i in range(0, len(binedges)-1):
    bincentres = np.append(bincentres, binedges[i]+binwidth/2)
  
  # if histogram points are not close to boundaries (i.e. within a bin of the
  # boundaries) then add zeros to histrogram edges
  if bincentres[0] - binwidth > low:
    # prepend a zero to n
    n = np.insert(n, 0, 0);
    
    # prepend a new bin centre at bincentres[0] - binwidth
    bincentres = np.insert(bincentres, 0, bincentres[0] - binwidth)
  else:
    # we're  closer to the boundary edge than the bin width then, so set a new
    # bin on the boundary with a value linearly extrapolated from the
    # gradiant of the adjacent points
    dx = bincentres[0] - low;
    
    # prepend low value to bins
    bincentres = np.insert(bincentres, 0, low)
    
    dn = n[1]-n[0]
    
    nbound = n[0] - (dn/binwidth)*dx
    
    # prepend to n
    n = np.insert(n, 0, nbound)
  
  # now the other end!
  if bincentres[-1] + binwidth < high:
    # append a zero to n
    n = np.append(n, 0)
    
    # append a new bin centre at bincentres[end] + binwidth
    bincentres = np.append(bincentres, bincentres[-1] + binwidth)
  else:
    dx = high - bincentres[-1];
    
    # prepend low value to bins
    bincentres = np.append(bincentres, high)
    
    dn = n[-1]-n[-2]
    
    nbound = n[-1] + (dn/binwidth)*dx
    
    # prepend to n
    n = np.append(n, nbound)
  
  # now calculate area and normalise
  area = np.trapz(n, x=bincentres)
  
  ns = np.array([])
  for i in range(0, len(bincentres)):
    ns = np.append(ns, float(n[i])/area)
  
  return ns, bincentres

# create a Tukey window of length N
def tukey_window(N, alpha=0.5):
  # if alpha >= 1 just return a Hanning window
  if alpha >= 1:
    return np.hanning(N)
  
  # get x values at which to calculate window
  x = np.linspace(0, 1, N)
  
  # initial square window
  win = np.ones(x.shape)

  # get the left-hand side of the window  0 <= x < alpha/2  
  lhs = x<alpha/2
  win[lhs] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[lhs] - alpha/2) ))

  # get right hand side condition 1 - alpha / 2 <= x <= 1
  rhs = x>=(1 - alpha/2)
  win[rhs] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[rhs] - 1 + alpha/2)))
  
  return win

# create a function for plotting the absolute value of Bk data (read in from
# data files) and an averaged 1 day amplitude spectral density spectrogram for
# each IFO
def plot_Bks_ASDs( Bkdata, ifos, delt=86400, plotpsds=True,
                   plotfscan=False, removeoutlier=None, mplparams=False ):
  # create list of figures
  Bkfigs = []
  psdfigs = []
  fscanfigs = []
  asdlist = []
  
  # set some matplotlib defaults for amplitude spectral density
  if not mplparams:
    mplparams = { \
      'backend': 'Agg',
      'text.usetex': True, # use LaTeX for all text
      'axes.linewidth': 0.5, # set axes linewidths to 0.5
      'axes.grid': True, # add a grid
      'grid.linwidth': 0.5,
      'font.family': 'serif',
      'font.size': 12 } 
      #'text.latex.preamble': \usepackage{xfrac} }
  
  matplotlib.rcParams.update(mplparams)
  # xfrac causes problems when compiling an eps (option clashes with graphicx)
  # so use nicefrac package instead
  #matplotlib.rcParams['text.latex.preamble']=r'\usepackage{xfrac}'
  matplotlib.rcParams['text.latex.preamble']=r'\usepackage{nicefrac}'
  
  # ifos line colour specs
  coldict = {'H1': 'r', 'H2': 'c', 'L1': 'g', 'V1': 'b', 'G1': 'm'}
  colmapdic = {'H1': 'Reds', 'H2': 'PuBu', 'L1': 'Greens', \
    'V1': 'Blues', 'G1': 'PuRd'}
  
  # there should be data for each ifo
  for i, ifo in enumerate(ifos):
    # get data for given ifo
    try:
      Bk = np.loadtxt(Bkdata[i])
    except:
      print "Could not open file ", Bkdata[i]
      exit(-1)
   
    # should be three lines in file
    gpstime = []
    
    # remove outliers at Xsigma by working out sigma from the peak of the
    # distribution of the log absolute value
    if removeoutlier:
      n, binedges = np.histogram(np.log(np.fabs(np.concatenate((Bk[:,1], \
Bk[:,2])))), 50)
      j = n.argmax(0)
      
      # standard devaition estimate
      stdest = math.exp((binedges[j]+binedges[j+1])/2)

      # get values within +/-8 sigma
      # test real parts
      vals = np.where(np.fabs(Bk[:,1]) < removeoutlier*stdest)
      Bknew = Bk[vals,:][-1]
      # test imag parts
      vals = np.where(np.fabs(Bknew[:,2]) < removeoutlier*stdest)
      Bk = Bknew[vals,:][-1]
   
    gpstime = Bk[:,0]
    
    # minimum time step between points (should generally be 60 seconds)
    mindt = min(np.diff(gpstime))
    
    Bkabs = np.sqrt(Bk[:,1]**2 + Bk[:,2]**2)
          
    # plot the time series of the data
    Bkfig = plt.figure(figsize=(11,3.5), dpi=200)
    Bkfig.subplots_adjust(bottom=0.15, left=0.09, right=0.94)
    
    tms = map(lambda x: x-gpstime[0], gpstime)
    
    plt.plot(tms, Bkabs, '.', color=coldict[ifo], markersize=1)
    plt.xlabel(r'GPS - %d' % int(gpstime[0]), fontsize=14, fontweight=100)
    plt.ylabel(r'$|B_k|$', fontsize=14, fontweight=100)
    #plt.title(r'$B_k$s for ' + ifo.upper(), fontsize=14)
    plt.xlim(tms[0], tms[-1])
    
    Bkfigs.append(Bkfig)
    
    if plotpsds or plotfscan:
      # create PSD by splitting data into days, padding with zeros to give a
      # sample a second, getting the PSD for each day and combining them
      totlen = gpstime[-1] - gpstime[0] # total data length
      
      # check mindt is an integer and greater than 1
      if math.fmod(mindt, 1) != 0. or mindt < 1:
        print "Error time steps between data points must be integers"
        exit(-1)
          
      count = 0
      npsds = 0
      totalpsd = np.zeros(math.floor(delt/60)) # add sum of PSDs
      
      asdlist.append(totalpsd)
      
      # zero pad the data and bin each point in the nearest 60s bin
      datazeropad = np.zeros(math.ceil(totlen/60.)+1, dtype=complex)
      
      idx = map(lambda x: math.floor((x/60.)+0.5), tms)
      for i in range(0, len(idx)):
        datazeropad[idx[i]] = complex(Bk[i,1], Bk[i,2])
      
      win = tukey_window(math.floor(delt/60), alpha=0.25)
        
      Fs = 1./60. # sample rate in Hz
      
      fscan, freqs, t = specgram(datazeropad, NFFT=int(math.floor(delt/60)), \
Fs=Fs, window=win)
      
      if plotpsds:
        fshape = fscan.shape
      
        for i in range(0, fshape[1]):
          # add to total psd if the psd does not contain zero
          if fscan[0,i] != 0.:
            # add psd onto total value
            totalpsd = map(lambda x, y: x+y, totalpsd, fscan[:,i])
      
            # count number of psds
            npsds = npsds+1
      
        # average the PSD and convert to amplitude spectral density
        totalpsd = map(lambda x: math.sqrt(x/npsds), totalpsd)
    
        # plot PSD
        psdfig = plt.figure(figsize=(4,3.5), dpi=200)
        psdfig.subplots_adjust(left=0.18, bottom=0.15)
    
        plt.plot(freqs, totalpsd, color=coldict[ifo])
        plt.xlim(freqs[0], freqs[-1])
        plt.xlabel(r'Frequency (Hz)', fontsize=14, fontweight=100)
        plt.ylabel(r'$h/\sqrt{\rm Hz}$', fontsize=14, fontweight=100)
        
        # convert frequency labels to fractions
        ax = plt.gca()
        xt = [-Fs/2., -Fs/4., 0., Fs/2., Fs/4.]
        ax.set_xticks(xt)
        xl = []
        for item in xt:
          if item == 0:
            xl.append('0')
          else:
            if item < 0:
              xl.append(r'$-\nicefrac{1}{%d}$' % (-1./item))
            else:
              xl.append(r'$\nicefrac{1}{%d}$' % (1./item))
        ax.set_xticklabels(xl)
        #plt.setp(ax.get_xticklabels(), fontsize=16)  # increase font size
        plt.tick_params(axis='x', which='major', labelsize=14)
    
        psdfigs.append(psdfig)
      
      if plotfscan:
        fscanfig = plt.figure(figsize=(11,3.5), dpi=200)
        fscanfig.subplots_adjust(bottom=0.15, left=0.09, right=0.94)
        
        extent = [tms[0], tms[-1], freqs[0], freqs[-1]]
        plt.imshow(np.sqrt(fscan), aspect='auto', extent=extent,
          interpolation=None, cmap=colmapdic[ifo], norm=colors.Normalize())
        plt.ylabel(r'Frequency (Hz)', fontsize=14, fontweight=100)
        plt.xlabel(r'GPS - %d' % int(gpstime[0]), fontsize=14, fontweight=100)
      
        # convert frequency labels to fractions
        ax = plt.gca()
        yt = [-Fs/2., -Fs/4., 0., Fs/2., Fs/4.]
        ax.set_yticks(yt)
        yl = []
        for item in yt:
          if item == 0:
            yl.append('0')
          else:
            if item < 0:
              yl.append(r'$-\nicefrac{1}{%d}$' % (-1./item))
            else:
              yl.append(r'$\nicefrac{1}{%d}$' % (1./item))
        ax.set_yticklabels(yl)
        #plt.setp(ax.get_yticklabels(), fontsize=16) 
        plt.tick_params(axis='y', which='major', labelsize=14)
        
        fscanfigs.append(fscanfig)
      
  return Bkfigs, psdfigs, fscanfigs, asdlist

# a function to create the signal model for a heterodyned triaxial pulsar given
# a signal GPS start time, signal duration (in seconds), sample interval
# (seconds), detector, and the parameters h0, cos(iota), psi (rads), initial
# phase phi0 (rads), right ascension (rads) and declination (rads) in a
# dictionary. The list of time stamps, and the real and imaginary parts of the
# signal are returned
def heterodyned_triaxial_pulsar(starttime, duration, dt, detector, pardict):    
                   
  # create a list of times stamps
  ts = []
  tmpts = starttime
  
  # create real and imaginary parts of the signal
  sr = []
  si = []
  
  i = 0
  while tmpts < starttime + duration:
    ts.append(starttime+(dt*i))
    
    # get the antenna response
    fp, fc = antenna_response(ts[i], pardict['ra'], pardict['dec'], \
                              pardict['psi'], detector)
    
    Xplus = 0.25*(1.+pardict['cosiota']*pardict['cosiota'])*pardict['h0']
    Xcross = 0.5*pardict['cosiota']*pardict['h0']
    Xpsinphi = Xplus*np.sin(pardict['phi0'])
    Xcsinphi = Xcross*np.sin(pardict['phi0'])
    Xpcosphi = Xplus*np.cos(pardict['phi0'])
    Xccosphi = Xcross*np.cos(pardict['phi0'])
    
    # create real part of signal
    sr.append(fp*Xpcosphi + fc*Xcsinphi)
    si.append(fp*Xpsinphi - fc*Xccosphi)
    
    tmpts = ts[i]+dt
    
    i = i+1;
  
  return ts, sr, si

# a function to create the signal model for a heterodyned pinned superfluid
# model pulsar a signal GPS start time, signal duration (in seconds),
# sample interval (seconds), detector, and the parameters I21, I31,
# cos(theta), lambda, cos(iota), psi (rads), initial phase phi0 (rads), right
# ascension (rads), declination (rads), distance (kpc) and frequency (f0) in a
# dictionary. The list of time stamps, and the real and imaginary parts of the
# 1f and 2f signals are returned in an array
def heterodyned_pinsf_pulsar(starttime, duration, dt, detector, pardict):
  iota = np.arccos(pardict['cosiota'])
  theta = np.arccos(pardict['costheta'])
  siniota = math.sin(iota)
  sintheta = math.sin(theta)
  sin2theta = math.sin( 2.0*theta )
  sinlambda = math.sin(pardict['lambda'])
  coslambda = math.cos(pardict['lambda'])
  sin2lambda = math.sin( 2.0*pardict['lambda'] )
  sinphi = math.sin(pardict['phi0'])
  cosphi = math.sin(pardict['phi0'])
  sin2phi = math.sin(2.0*pardict['phi0'])
  cos2phi = math.sin(2.0*pardict['phi0'])
 
  f2_r = pardict['f0'] * pardict['f0'] / pardict['dist'];
  
  """
    This model is a complex heterodyned time series for a pinned superfluid
    neutron star emitting at its roation frequency and twice its rotation
    frequency (as defined in Jones 2009)
  """
  
  Xplusf = -( f2_r / 2.0 ) * siniota * pardict['cosiota'];
  Xcrossf = -( f2_r / 2.0 ) * siniota;
  Xplus2f = -f2_r * ( 1.0 + pardict['cosiota'] * pardict['cosiota'] );
  Xcross2f = -f2_r * 2.0 * pardict['cosiota'];
  
  A1 = ( pardict['I21'] * coslambda * coslambda - pardict['I31'] ) * sin2theta;
  A2 = pardict['I21'] * sin2lambda * sintheta;
  B1 = pardict['I21'] * ( coslambda * coslambda * pardict['costheta'] * \
    pardict['costheta'] - sinlambda * sinlambda ) + pardict['I31'] * \
    sintheta * sintheta;
  B2 = pardict['I21'] * sin2lambda * pardict['costheta'];
  
  # create a list of times stamps
  ts1 = []
  ts2 = []
  tmpts = starttime
  
  # create real and imaginary parts of the 1f signal
  sr1 = []
  si1 = []
  sr2 = []
  si2 = []
  
  i = 0
  while tmpts < starttime + duration:
    ts1.append(starttime+(dt*i))
    ts2.append(starttime+(dt*i))
    
    # get the antenna response
    fp, fc = antenna_response(ts1[i], pardict['ra'], pardict['dec'], \
                              pardict['psi'], detector)
    
    # create the complex signal amplitude model at 1f
    sr1.append( fp * Xplusf * ( A1 * cosphi - A2 * sinphi ) + \
                fc * Xcrossf * ( A2 * cosphi + A1 * sinphi ) )
    
    si1.append( fp * Xplusf * ( A2 * cosphi + A1 * sinphi ) + \
                fc * Xcrossf * ( A2 * sinphi - A1 * cosphi ) )
   
    
    # create the complex signal amplitude model at 2f
    sr2.append( fp * Xplus2f * ( B1 * cos2phi - B2 * sin2phi ) + \
                fc * Xcross2f * ( B2 * cos2phi + B1 * sin2phi ) )
    
    si2.append( fp * Xplus2f * ( B2 * cos2phi + B1 * sin2phi ) + \
                fc * Xcross2f * ( B2 * sin2phi - B1 * cos2phi ) )
    
    tmpts = ts1[i]+dt
    
    i = i+1;
 
  # combine data into 1 array
  ts = np.vstack([ts1, ts2])
  sr = np.vstack([sr1, sr2])
  si = np.vstack([si1, si2])
    
  return ts, sr, si
  
# function to get the antenna response for a given detector. This is based on
# the response function in pylal/antenna.py. It takes in a GPS time, right
# ascension (rads), declination (rads), polarisation angle (rads) and a
# detector name e.g. H1, L1, V1. The plus and cross polarisations are returned.
def antenna_response( gpsTime, ra, dec, psi, det ):
  gps = LIGOTimeGPS( gpsTime )
  gmst_rad = date.XLALGreenwichMeanSiderealTime(gps)

  # create detector-name map
  detMap = {'H1': 'LHO_4k', 'H2': 'LHO_2k', 'L1': 'LLO_4k',
            'G1': 'GEO_600', 'V1': 'VIRGO', 'T1': 'TAMA_300'}
  try:
    detector=detMap[det]
  except KeyError:
    raise ValueError, "ERROR. Key %s is not a valid detector name."\
          % (det)

  # get detector
  if detector not in tools.cached_detector.keys():
    raise ValueError, "%s is not a cached detector.  "\
          "Cached detectors are: %s" \
          % (det, tools.cached_detector.keys())

  # get the correct response data
  response = tools.cached_detector[detector].response

  # actual computation of antenna factors
  fp, fc = inject.XLALComputeDetAMResponse(response, ra, dec,
                                                    psi, gmst_rad)
  
  return fp, fc

# a function to inject a heterodyned triaxial pulsar signal into noise of a
# given level for detectors. it will return the signal + noise and the optimal
# SNR. If an snrscale value is passed to the function the signal will be scaled
# to that SNR. nsigs
def inject_pulsar_signal(starttime, duration, dt, detectors, pardict, \
                         model='triaxial', npsds=None, snrscale=None):
  if model == 'triaxial':
    freqfac = [2.0]
  elif model == 'pinsf':
    freqfac = [1.0, 2.0]
  else:
    print >> sys.stderr, "Model must be triaxial or pinsf"
    sys.exit(1)
  
  # if detectors is just a string (i.e. one detector) then make it a list
  if isinstance(detectors, basestring):
    detectors = [detectors]
  
  # if not noise sigma's are given then generate a noise level from the given
  # detector noise curves
  if npsds is None:
    npsds = []
    
    for det in detectors:
      for frf in freqfac:
        psd = detector_noise( det, frf*pardict['f0'] )
      
        # convert to time domain standard devaition shared between real and
        # imaginary signal
        ns = np.sqrt( (psd/2.0)/(2.0*dt) )
      
        npsds.append(ns)
  else:
    # convert input psds into time domain noise standard deviation
    tmpnpsds = []
    
    count = 0
    for j, det in enumerate(detectors):
      for frf in freqfac:
        if len(npsds) == 1:
          tmpnpsds.append( (npsds/2.0)/(2.0*dt) )
        else:
          tmpnpsds.append( (npsds[count]/2.0)/(2.0*dt) )
        
        count = count+1
        
    npsds = tmpnpsds
  
  if model == 'triaxial' and len(detectors) != len(npsds):
    raise ValueError, "Number of detectors %d not the same as number of "\
                      "noises %d" % (len(detectors), len(npsds))
  
  if model == 'pinsf' and 2*len(detectors) != len(npsds):
    raise ValueError, "Number of detectors %d not half the number of "\
                      "noises %d" % (len(detectors), len(npsds))
  
  tss = np.array([])
  srs = np.array([])
  sis = np.array([])
 
  snrtot = 0
  for j, det in enumerate(detectors):
    # create the pulsar signal
    if model == 'triaxial':
      ts, sr, si = heterodyned_triaxial_pulsar(starttime, duration, dt, det, \
                                               pardict)
    
      if j == 0:
        tss = np.append(tss, ts)
        srs = np.append(srs, sr)
        sis = np.append(sis, si)
      else:
        tss = np.vstack([tss, ts])
        srs = np.vstack([srs, sr])
        sis = np.vstack([sis, si])
    
      # get SNR
      snrtmp = get_optimal_snr( sr, si, npsds[j] )
    elif model == 'pinsf':
      ts, sr, si = heterodyned_pinsf_pulsar(starttime, duration, dt, det, \
                                            pardict)
      
      snrtmp = 0
      for k, frf in enumerate(freqfac):
        if j == 0 and k == 0:
          tss = np.append(tss, ts[k][:])
          srs = np.append(srs, sr[k][:])
          sis = np.append(sis, si[k][:])
        else:
          tss = np.vstack([tss, ts[k][:]])
          srs = np.vstack([srs, sr[k][:]])
          sis = np.vstack([sis, si[k][:]])
        
        snrtmp2 = get_optimal_snr( sr[k][:], si[k][:], npsds[2*j+k] )
        snrtmp = snrtmp + snrtmp2*snrtmp2
        
      snrtmp = np.sqrt(snrtmp)
    
    snrtot = snrtot + snrtmp*snrtmp
    
  # total multidetector/data stream snr
  snrtot = np.sqrt(snrtot)
  
  # add noise and rescale signals if necessary
  if snrscale is not None:
    snrscale = snrscale / snrtot
  else:
    snrscale = 1
  
  i = 0
  for det in detectors:
    # for triaxial model
    if len(freqfac) == 1 and model == 'triaxial':
      # generate random numbers
      rs = np.random.randn(len(ts), 2)
      
      for j, t in enumerate(ts):
        if len(tss.shape) == 1:
          srs[j] = snrscale*srs[j] + npsds[i]*rs[j][0]
          sis[j] = snrscale*sis[j] + npsds[i]*rs[j][1]
        else:
          srs[i][j] = snrscale*srs[i][j] + npsds[i]*rs[j][0]
          sis[i][j] = snrscale*sis[i][j] + npsds[i]*rs[j][1]

      i = i+1
    elif len(freqfac) == 2 and model == 'pinsf':
      # generate random numbers
      rs = np.random.randn(len(ts[0][:]), 4)
      
      for j, t in enumerate(ts[0][:]):
        srs[i][j] = snrscale*srs[i][j] + npsds[i]*rs[j][0]
        sis[i][j] = snrscale*sis[i][j] + npsds[i]*rs[j][1]
        srs[i+1][j] = snrscale*srs[i+1][j] + npsds[i+1]*rs[j][2]
        sis[i+1][j] = snrscale*sis[i+1][j] + npsds[i+1]*rs[j][3]
        
      i = i+2
    else:
      print >> sys.stderr, "Something wrong with injection"
      sys.exit(1)
  
  snrtot = snrtot*snrscale 
  
  return tss, srs, sis, snrtot
  
# function to create a time domain PSD from theoretical
# detector noise curves. It takes in the detector name and the frequency at
# which to generate the noise.
#
# The noise models are taken from those in lalsimulation/src/LALSimNoisePSD.c
def detector_noise( det, f ):
  if det == 'AV1': # Advanced Virgo
    x = np.log(f / 300.);
    x2 = x*x;

    asd = 1.259e-24 * ( 0.07*np.exp(-0.142 - 1.437*x + 0.407*x2)
                       + 3.1*np.exp(-0.466 - 1.043*x - 0.548*x2)
                       + 0.4*np.exp(-0.304 + 2.896*x - 0.293*x2)
                       + 0.09*np.exp(1.466 + 3.722*x - 0.984*x2) )

    return asd*asd;
  elif det == 'H1' or det == 'L1': # iLIGO SRD
    aseis = 1.57271;
    pseis = -14.0;
    athrm = 3.80591e-19;
    pthrm = -2.0;
    ashot = 1.12277e-23;
    fshot = 89.3676;
    seis = aseis * aseis * np.power(f, 2.0*pseis);
    thrm = athrm * athrm * np.power(f, 2.0*pthrm);
    shot = ashot * ashot * (1.0 + np.power(f / fshot, 2.0));
    
    return seis + thrm + shot;
  elif det == 'G1': # GEO_600
    x = f/150.
    seismic = np.power(10.,-16.) * np.power(x,-30.)
    thermal = 34. / x
    shot = 20. * (1 - np.power(x,2.) + 0.5 * np.power(x,4.)) / \
           (1. + 0.5 * np.power(x,2.))

    return 1e-46*(seismic + thermal + shot)
  elif det == 'V1': # initial Virgo
    x = f/500.;
    s0 = 10.2e-46;

    return s0*( np.power(7.87*x,-4.8) + 6./17./x + 1. + x*x)
  else:
    raise ValueError, "%s is not a recognised detector" % (det)

# function to calculate the optimal SNR of a heterodyned pulsar signal - it
# takes in a complex signal model and noise standard deviation
def get_optimal_snr( sr, si, sig ):  
  ss = 0
  # sum square of signal
  for i, s in enumerate(sr):
    ss = ss + sr[i]*sr[i] + si[i]*si[i]
   
  return np.sqrt( ss / (sig*sig) )

# use the Gelman-Rubins convergence test for MCMC chains, where chains is a list
# of MCMC numpy chain arrays - this copies the gelman_rubins function in
# bayespputils
def gelman_rubins(chains):  
  chainMeans = [np.mean(data) for data in chains]
  chainVars = [np.var(data) for data in chains]
  
  BoverN = np.var(chainMeans)
  W = np.mean(chainVars)
  
  sigmaHat2 = W + BoverN
  
  m = len(chains)
  
  VHat=sigmaHat2 + BoverN/m
  
  R = VHat/W
  
  return R