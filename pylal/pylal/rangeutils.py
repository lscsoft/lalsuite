# Copyright (C) 2012 Duncan M. Macleod
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
This module provides a bunch of user-friendly wrappers to the SWIG-bound REAL8TimeSeries and REAL8FrequencySeries objects and their associated functions.
"""

from __future__ import division

import numpy
from scipy import integrate
import lal
from pylal import git_version

__author__  = "Duncan M. Macleod <duncan.macleod@ligo.org>"
__version__ = git_version.id
__date__    = git_version.date

# =============================================================================
# Inspiral range
# =============================================================================

def inspiralrange(f, S, snr=8, m1=1.4, m2=1.4, fmin=10, fmax=None,\
                  horizon=False):
    """
    Calculate the sensitivity distance to an inspiral with the given
    masses, for the given signal-to-noise ratio. See the following
    reference for information on the integral performed:

    https://dcc.ligo.org/cgi-bin/private/DocDB/ShowDocument?docid=27267

    @param f: frequency array
    @type  f: C{numpy.array}
    @param S: power spectral density array
    @type  S: C{numpy.array}
    @param snr: signal-to-noise ratio at which to calculate range
    @type  snr: C{float}
    @param m1: mass (in solar masses) of first binary component,
        default: 1.4
    @type  m1: C{float}
    @param m2: mass (in solar masses) of second binary component,
        default: 1.4
    @type  m2: C{float}
    @param fmin: minimum frequency limit of integral, default: 10 Hz
    @type  fmin: C{float}
    @param fmax: maximum frequency limit of integral, default: ISCO
    @type  fmax: C{float}
    @param horizon: return horizon distance in stead of angle-averaged,
        default: False
    @type  horizon: C{bool}

    @return: sensitive distance to an inspiral (in solar masses) for the
        given PSD and parameters
    @rtype: C{float}
    """
    # compute chirp mass and symmetric mass ratio
    mtot   = m1+m2
    mchirp = (m1*m2)**(0.6) / mtot**(0.4)
  
    # get mchirp in solar masses and compute integral prefactor
    mchirp *= lal.LAL_MSUN_SI
    pre = (5 * lal.LAL_C_SI**(1/3) *\
           (mchirp * lal.LAL_G_SI / lal.LAL_C_SI**2)**(5/3) * 1.77**2) /\
          (96 * lal.LAL_PI ** (4/3) * snr**2)


    # compute ISCO
    fisco = lal.LAL_C_SI**3\
            / (lal.LAL_G_SI * lal.LAL_MSUN_SI * 6**1.5 * lal.LAL_PI * mtot)
    if not fmax:
        fmax = fisco
    elif fmax > fisco:
        warnings.warn("Upper frequency bound greater than %s-%s ISCO "\
                      "frequency of %.2g Hz, using ISCO" % (m1,m2,fisco))
        fmax = fisco

    # integrate
    condition = (f >= fmin) & (f < fmax)
    integrand = f[condition]**(-7/3)/S[condition]
    result = integrate.trapz(integrand, f[condition])

    d = (pre*result) ** (1/2) / (lal.LAL_PC_SI*1e6)
    return d

# =============================================================================
# Squeezing ratio
# =============================================================================

def squeezing(f, S, dc=None, opgain=None, cavpf=None, fmin=1000, fmax=2000):
    """
    Calculate squeezing factor based on observed noise in given band
    and predicted shot noise

    @param f: frequency array
    @type  f: C{numpy.array}
    @param S: power spectral density array
    @type  S: C{numpy.array}
    @param dc: 
    @type  dc: C{float}
    @param opgain: 
    @type  opgain: C{float}
    @param cavpf: 
    @type  cavpf: C{float}
    @param fmin: minimum frequency limit of integral, default: 1000 Hz
    @type  fmin: C{float}
    @param fmax: maximum frequency limit of integral, default: 2000 Hz
    @type  fmax: C{float}

    @return: squeezing ratio in dB
    @rtype: C{float}
    """

    # set up frequency band for squeezing estimation
    condition = (f >= fmin) & (f < fmax)

    # compute model shot noise spectrum
    model = abs(1+1j*f[condition]/cavpf) * (dc)**(1/2) /opgain

    # compare to actual noise spectrum
    d =  numpy.median(model /S[condition]**(1/2))
    
    # convert to dB
    d = 20 * numpy.log10(d)

    return d

# =============================================================================
# Burst range
# =============================================================================

def fdependent_burstrange(f, S, snr=8, E=1e-2):
    """
    Calculate the sensitive distance to a GW burst with the given intrinsic
    energy for the given signal-to-noise ratio snr, as a function of f.

    @param f: frequency of interest
    @type  f: C{float} or C{numpy.array}
    @param S: power spectral density
    @type  S: C{numpy.array}
    @param snr: signal-to-noise ratio of burst
    @type  snr: C{float}
    @param E: instrinsic energy of burst, default: grb-like 0.01
    @type  E: C{float}

    @return: sensitive distance in pc at which a GW burst with the given
        energy would be detected with the given SNR, as a function of
        it's frequency
    @rtype: C{float}
    """
    A = ((lal.LAL_G_SI * E * lal.LAL_MSUN_SI * 0.4)\
         / (lal.LAL_PI**2 * lal.LAL_C_SI))**(1/2) / lal.LAL_PC_SI 
    return A / (snr * S**(1/2) * f)

def burstrange(f, S, snr=8, E=1e-2, fmin=0, fmax=None, unit="Mpc"):
    """
    Calculate the sensitive distance to a GW burst with the given intrinsic
    energy for the given signal-to-noise ratio snr, integrated over frequency.

    @param f: frequency of interest
    @type  f: C{float} or C{numpy.array}
    @param S: power spectral density
    @type  S: C{numpy.array}
    @param snr: signal-to-noise ratio of burst, default: 8
    @type  snr: C{float}
    @param E: instrinsic energy of burst, default: grb-like 0.01
    @type  E: C{float}
    @param fmin: minimum frequency limit of integral, default: 10 Hz
    @type  fmin: C{float}
    @param fmax: maximum frequency limit of integral, default: ISCO
    @type  fmax: C{float}

    @return: sensitive distance at which a GW burst with the given
        energy would be detected with the given SNR, integrated over
        frequency.
    @rtype: C{float}
    """
    if not fmin: 
        fmin = f.min()
    if not fmax:
        fmax = f.max()

    # restrict to band
    condition = (f >= fmin) & (f < fmax)

    # integrate
    integrand = fdependent_burstrange(f[condition], S[condition], snr, E)**3
    result = integrate.trapz(integrand, f[condition])
    
    d = (result / (fmax-fmin))**(1/3)
    if unit == "Mpc":
        d = d/1e6
    elif unit == "kpc":
        d = d/1e3
    elif unit == "pc":
        d = d
    else:
        raise ValueError("Unrecognized unit: %s" % unit)

    return d
