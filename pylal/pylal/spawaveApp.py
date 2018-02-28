#spawaveApp.py
#Created by Satya Mohapatra on 2/26/10.
#Copyright (C) 2010 Satya Mohapatra
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
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

"""
This module is an enhancement based on the module spawaveform available in pylal. 
This module can do the following:
1) Find the expected SNR of a spin-aligned IMR waveform (see Ajith et al.) with respect to initial LIGO SRD.
  submodule:IMRsnr
2) Find the hrss of the waveform.
  submodule:IMRhrss
3) Find the peak amplitude of the waveform.
  submodule:IMRpeakAmp
4) Find the target burst frequency of the waveform:
 submodule:IMRtargetburstfreq
 
"""



from pylal import spawaveform
import math
from scipy import interpolate
import scipy
import numpy
import time


from pylal import git_version
__author__ = "Satya Mohapatra <satya@physics.umass.edu>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date

design = numpy.loadtxt('ligo4k_design.txt')
interpolatation = interpolate.interp1d(design[:,0],design[:,1]) 
#interpolatation = interpolate.splrep(design[:,0],design[:,1],s=0) 

def IMRsnr(m1,m2,spin1z,spin2z,d):

    """
    IMRsnr finds the SNR of the waveform with Initial LIGO SRD for a given source parameters and the source distance.
	usage: IMRsnr(m1,m2,spin1z,spin2z,distance)
	e.g.
	spawaveApp.IMRsnr(30,40,0.45,0.5,100)

	"""

    chi = spawaveform.computechi(m1, m2, spin1z, spin2z)   
    imrfFinal = spawaveform.imrffinal(m1, m2, chi, 'fcut')
    fLower = 10.0
    order = 7
    dur = 2**numpy.ceil(numpy.log2(spawaveform.chirptime(m1,m2,order,fLower)))
    sr = 2**numpy.ceil(numpy.log2(imrfFinal*2))
    deltaF = 1.0 / dur
    deltaT = 1.0 / sr
    s = numpy.empty(sr * dur, 'complex128')	
    spawaveform.imrwaveform(m1, m2, deltaF, fLower, s, spin1z, spin2z)
    S = numpy.abs(s)
    x = scipy.linspace(fLower, imrfFinal, numpy.size(S))
    N = interpolatation(x)
    #N = interpolate.splev(x,interpolatation)
    SNR = math.sqrt(numpy.sum(numpy.divide(numpy.square(S),numpy.square(N))))/d
    return SNR

def IMRhrss(m1,m2,spin1z,spin2z,d):

    """
	IMRhrss finds the SNR of the waveform for a given source parameters and the source distance.
	usage: IMRhrss(m1,m2,spin1z,spin2z,distance)
	e.g.
	spawaveApp.IMRhrss(30,40,0.45,0.5,100)

	"""

    chi = spawaveform.computechi(m1, m2, spin1z, spin2z)   
    imrfFinal = spawaveform.imrffinal(m1, m2, chi, 'fcut')
    fLower = 10.0
    order = 7
    dur = 2**numpy.ceil(numpy.log2(spawaveform.chirptime(m1,m2,order,fLower)))
    sr = 2**numpy.ceil(numpy.log2(imrfFinal*2))
    deltaF = 1.0 / dur
    deltaT = 1.0 / sr
    s = numpy.empty(sr * dur, 'complex128')	
    spawaveform.imrwaveform(m1, m2, deltaF, fLower, s, spin1z, spin2z)
    s = numpy.abs(s)
    s = numpy.square(s)
    hrss = math.sqrt(numpy.sum(s))/d
    return hrss

def IMRpeakAmp(m1,m2,spin1z,spin2z,d):

    """
	IMRpeakAmp finds the peak amplitude of the waveform for a given source parameters and the source distance.
	usage: IMRpeakAmp(m1,m2,spin1z,spin2z,distance)
	e.g.
	spawaveApp.IMRpeakAmp(30,40,0.45,0.5,100)

	"""

    chi = spawaveform.computechi(m1, m2, spin1z, spin2z)   
    imrfFinal = spawaveform.imrffinal(m1, m2, chi, 'fcut')
    fLower = 10.0
    order = 7
    dur = 2**numpy.ceil(numpy.log2(spawaveform.chirptime(m1,m2,order,fLower)))
    sr = 2**numpy.ceil(numpy.log2(imrfFinal*2))
    deltaF = 1.0 / dur
    deltaT = 1.0 / sr
    s = numpy.empty(sr * dur, 'complex128')	
    spawaveform.imrwaveform(m1, m2, deltaF, fLower, s, spin1z, spin2z)
    s = scipy.ifft(s)
    #s = numpy.abs(s)
    s = numpy.real(s)
    max = numpy.max(s)/d
    return max

def IMRtargetburstfreq(m1,m2,spin1z,spin2z):

    """
	IMRtargetburstfreq finds the peak amplitude of the waveform for a given source parameters and the source distance.
	usage: IMRtargetburstfreq(m1,m2,spin1z,spin2z)
	e.g.
	spawaveApp.IMRtargetburstfreq(30,40,0.45,0.5)

	"""
    chi = spawaveform.computechi(m1, m2, spin1z, spin2z)
    fFinal = spawaveform.ffinal(m1,m2,'schwarz_isco')   
    imrfFinal = spawaveform.imrffinal(m1, m2, chi, 'fcut')
    fLower = 40.0
    order = 7
    dur = 2**numpy.ceil(numpy.log2(spawaveform.chirptime(m1,m2,order,fFinal)))
    sr = 2**numpy.ceil(numpy.log2(imrfFinal*2))
    deltaF = 1.0 / dur
    deltaT = 1.0 / sr
    s = numpy.empty(sr * dur, 'complex128')	
    spawaveform.imrwaveform(m1, m2, deltaF, fFinal, s, spin1z, spin2z)
    #S = numpy.real(s)
    S = numpy.abs(s)
    x = scipy.linspace(fFinal, imrfFinal, numpy.size(S))
    N = interpolatation(x)
    #N = interpolate.splev(x,interpolatation)
    ratio = numpy.divide(numpy.square(S),numpy.square(N))
    #ratio = numpy.divide(S,N)
    maxindex = numpy.argmax(ratio)
    maxfreq = x[maxindex]
    return maxfreq


def IMRfinalspin(m1, m2,spin1z,spin2z):

    """
	The final spin calculation is based on Rezolla et al. 
	IMRfinalspin final spin of the end product black hole.	
	usage: IMRfinalspin(m1,m2,spin1z,spin2z)
	e.g.
	spawaveApp.IMRfinalspin(30.,40.,0.45,0.5)

	"""

    s4 = -0.129
    s5 = -0.384
    t0 = -2.686
    t2 = -3.454
    t3 = 2.353
    spin1x = 0.
    spin1y = 0.
    spin2x = 0.
    spin2y = 0.
    M = m1 + m2
    q = m2/m1
    eta = m1*m2/(m1+m2)**2
    a1 = math.sqrt(spin1x**2 + spin1y**2 + spin1z**2)
    a2 = math.sqrt(spin2x**2 + spin2y**2 + spin2z**2)
    if (a1 != 0) and (a2 != 0): cosa = (spin1x*spin2x + spin1y*spin2y + spin1z*spin2z)/(a1*a2)
    else:cosa = 0
    if a1 != 0: cosb = spin1z/a1
    else: cosb = 0
    if a2 != 0: cosc = spin2z/a2
    else: cosc = 0
    l = (s4/(1+q**2)**2 * (a1**2 + (a2**2)*(q**4) + 2*a1*a2*(q**2)*cosa) + (s5*eta + t0 + 2)/(1+q**2) * (a1*cosb + a2*(q**2)*cosc) + 2*math.sqrt(3.) + t2*eta + t3*(eta**2))
    afin = 1/(1+q)**2 * math.sqrt(a1**2 + (a2**2)*(q**4) + 2*a1*a2*(q**2)*cosa + 2*(a1*cosb + a2*(q**2)*cosc)*l*q + (l**2)*(q**2))
    return afin

