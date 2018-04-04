# Copyright (C) 2012,2013  Drew Keppel
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
# ============================================================================
#
#                                   Preamble
#
# ============================================================================
#


"""
This module provides instances for some of the common cases for use with the
coherent inspiral metric module.
"""

import pylab
import scipy
import numpy
import lalsimulation
from pylal import coherent_inspiral_metric as metric
from scipy import pi,sin,cos
#from pylal import inject
#from pylal.xlal.tools import cached_detector

__author__ = "Drew Keppel <drew.keppel@ligo.org>"


def f_for_fft(fLow, fNyq, deltaF):
	f = scipy.arange(2*(fNyq/deltaF))*deltaF
	f[fNyq/deltaF+1:] = -f[fNyq/deltaF-1:0:-1]
	return f

def f_PSD_from_file(filename, fLow, fNyq, deltaF):
	"""
	Read a detector ascii ASD file and return the PSD and frequency vector
	for use with ffts.
	"""
	f_in,S_in = numpy.loadtxt(filename, unpack=True)
	f = numpy.linspace(fLow,fNyq,scipy.ceil((fNyq-fLow)/deltaF)+1)
	S = pylab.interp(f, f_in, S_in)
	# packing is of the form:
	# [0 deltaF 2*deltaF ... fNyquist-deltaF fNyquist -fNyquist+deltaF ... -2*deltaF -deltaF]
	PSD = scipy.zeros(2*(fNyq/deltaF), dtype='float')+scipy.inf
	PSD[round(fLow/deltaF):fNyq/deltaF+1] = S**2
	if -round(fLow/deltaF) == 0:
		PSD[fNyq/deltaF+1:] = S[-2:0:-1]**2
	else:
		PSD[fNyq/deltaF+1:-round(fLow/deltaF)] = S[-2:0:-1]**2
	f = f_for_fft(fLow, fNyq, deltaF)
	nNyq = round(fNyq/deltaF)
	nLow = round(fLow/deltaF)
	PSD = scipy.zeros(2*nNyq)+scipy.inf
	PSD[nLow:nNyq+1] = S**2
	if -nLow == 0:
		PSD[nNyq+1:] = S[-2:0:-1]**2
	else:
		PSD[nNyq+1:-nLow] = S[-2:0:-1]**2
	return f,PSD

def f_PSD_from_func(psd_func, fLow, fNyq, deltaF):
	"""
	Construct a (frequency, PSD) vector for use with FFTs from the specified
	lalsimulation PSD functions.
	"""
	f = numpy.linspace(fLow,fNyq,scipy.ceil((fNyq-fLow)/deltaF)+1)
	S = numpy.array([psd_func(x) for x in f])
	# packing is of the form:
	# [0 deltaF 2*deltaF ... fNyquist-deltaF fNyquist -fNyquist+deltaF ... -2*deltaF -deltaF]
	f = f_for_fft(fLow, fNyq, deltaF)
	nNyq = round(fNyq/deltaF)
	nLow = round(fLow/deltaF)
	PSD = scipy.zeros(2*nNyq)+scipy.inf
	PSD[nLow:nNyq+1] = S
	if -nLow == 0:
		PSD[nNyq+1:] = S[-2:0:-1]
	else:
		PSD[nNyq+1:-nLow] = S[-2:0:-1]
	return f,PSD

def DegMinSec2Rad(sign,deg,minutes,seconds):
	"""
	A routine to convet from Degrees, Minutes, Seconds to radians.
	"""
	return sign*(deg + minutes/60. + seconds/3600.)*(pi/180.)

def lat_lon_2_vertex(lat,lon):
	"""
	A routine to return the location of a detector's vertex in 3D
	Cartesean Coordinates given the latitude and longitude of the
	detector. This routine approximates teh Earth as a sphere.
	"""
	return (metric.R_earth*cos(lon)*cos(lat), metric.R_earth*sin(lon)*cos(lat), metric.R_earth*sin(lat))

def lat_lon_ori_2_xarm(lat, lon, ori):
	"""
	A routine to return the x-arm vector of a detector in 3D
	Cartesean Coordinates given the latitude, longitude, and orientation
	of the detector. This routine approximates the Earth as a sphere.
	"""
	ori += pi/4.
	dx0_dRA,dy0_dRA,dz0_dRA = -sin(lon)*cos(lat), cos(lon)*cos(lat), 0
	dx0_ddec,dy0_ddec,dz0_ddec = -cos(lon)*sin(lat), -sin(lon)*sin(lat), cos(lat)
	x = dx0_dRA*sin(ori) + dx0_ddec*cos(ori)
	y = dy0_dRA*sin(ori) + dy0_ddec*cos(ori)
	z = dz0_dRA*sin(ori) + dz0_ddec*cos(ori)
	return (x,y,z)

def lat_lon_ori_2_yarm(lat, lon, ori):
	"""
	A routine to return the y-arm vector of a detector in 3D
	Cartesean Coordinates given the latitude, longitude, and orientation
	of the detector. This routine approximates the Earth as a sphere.
	"""
	ori -= pi/4.
	dx0_dRA,dy0_dRA,dz0_dRA = -sin(lon)*cos(lat), cos(lon)*cos(lat), 0
	dx0_ddec,dy0_ddec,dz0_ddec = -cos(lon)*sin(lat), -sin(lon)*sin(lat), cos(lat)
	x = dx0_dRA*sin(ori) + dx0_ddec*cos(ori)
	y = dy0_dRA*sin(ori) + dy0_ddec*cos(ori)
	z = dz0_dRA*sin(ori) + dz0_ddec*cos(ori)
	return (x,y,z)


# PSDs obtainable from:
# Adv. LIGO:
# https://dcc.ligo.org/cgi-bin/DocDB/ShowDocument?docid=2974

# Adv. Virgo:
# https://wwwcascina.virgo.infn.it/advirgo/
# https://wwwcascina.virgo.infn.it/advirgo/docs/AdV_refsens_100512.txt

# LCGT (detuned):
# http://gwcenter.icrr.u-tokyo.ac.jp/en/researcher/parameter
# https://granite.phys.s.u-tokyo.ac.jp/trac/LCGT/browser/trunk/sensitivity/spectrum/BW2009_VRSED.dat

# LCGT (broadband):
# https://granite.phys.s.u-tokyo.ac.jp/trac/LCGT/browser/trunk/sensitivity/spectrum/BW2009_VRSEB.dat

# detector geometries from LALDetectors.h
def make_LHO(fLow=10., fNyq=2048., deltaF=0.1, psd_filename=None):
	if psd_filename is None:
		f, psd = f_PSD_from_func(
			lalsimulation.SimNoisePSDaLIGOZeroDetHighPower, fLow, fNyq, deltaF)
	else:
		f, psd = f_PSD_from_file(psd_filename,fLow,fNyq,deltaF)
	detector = metric.Detector(
		'H',
		(-0.22389266154,0.79983062746,0.55690487831),
		(-0.91397818574,0.02609403989,-0.40492342125),
		(-2.16141492636e+06,-3.83469517889e+06,4.60035022664e+06),
		f=f,
		psd=psd
		)
	detector.set_required_moments()
	return detector

def make_LLO(fLow=10., fNyq=2048., deltaF=0.1, psd_filename=None):
	if psd_filename is None:
		f, psd = f_PSD_from_func(
			lalsimulation.SimNoisePSDaLIGOZeroDetHighPower, fLow, fNyq, deltaF)
	else:
		f, psd = f_PSD_from_file(psd_filename,fLow,fNyq,deltaF)
	detector = metric.Detector(
		'L',
		(-0.95457412153,-0.14158077340,-0.26218911324),
		(0.29774156894,-0.48791033647,-0.82054461286),
		(-7.42760447238e+04,-5.49628371971e+06,3.22425701744e+06),
		f=f,
		psd=psd
		)
	detector.set_required_moments()
	return detector

def make_Virgo(fLow=10., fNyq=2048., deltaF=0.1, psd_filename=None):
	if psd_filename is None:
		f, psd = f_PSD_from_func(
			lalsimulation.SimNoisePSDAdvVirgo, fLow, fNyq, deltaF)
	else:
		f, psd = f_PSD_from_file(psd_filename,fLow,fNyq,deltaF)
	detector = metric.Detector(
		'V',
		(-0.70045821479,0.20848948619,0.68256166277),
		(-0.05379255368,-0.96908180549,0.24080451708),
		(4.54637409900e+06,8.42989697626e+05,4.37857696241e+06),
		f=f,
		psd=psd
		)
	detector.set_required_moments()
	return detector

def make_GEO(fLow=10., fNyq=2048., deltaF=0.1, psd_filename=None):
	if psd_filename is None:
		# Adv. LIGO but 10x less sensitive
		f, psd = f_PSD_from_func(
			lalsimulation.SimNoisePSDaLIGOZeroDetHighPower, fLow, fNyq, deltaF)
		psd *= 100
	else:
		f, psd = f_PSD_from_file(psd_filename,fLow,fNyq,deltaF)
	detector = metric.Detector(
		'G',
		(-0.44530676905, 0.86651354130, 0.22551311312),
		(-0.62605756776, -0.55218609524, 0.55058372486),
		(3.85630994926e+06, 6.66598956317e+05, 5.01964141725e+06),
		f=f,
		psd=psd
		)
	detector.set_required_moments()
	return detector

# detector geometry from http://gwcenter.icrr.u-tokyo.ac.jp/en/
def make_KAGRA(fLow=10., fNyq=2048., deltaF=0.1, psd_filename=None):
	if psd_filename is None:
		f, psd = f_PSD_from_func(
			lalsimulation.SimNoisePSDKAGRA, fLow, fNyq, deltaF)
	else:
		f, psd = f_PSD_from_file(psd_filename,fLow,fNyq,deltaF)
	detector = metric.Detector(
		'K',
		(-0.3758971922940067, -0.83615832915098853, 0.3994252738835008),
		(0.71644138561445658, 0.011148317434329416, 0.6975582097554448),
		(-3.777336055e+06, 3.484898386e+06, 3.765313690e+06),
		f=f,
		psd=psd
		)
	detector.set_required_moments()
	return detector

# detector geometries from arxiv:1102.5421v2
def make_LCGT(fLow=10., fNyq=2048., deltaF=0.1, psd_filename=None):
	if psd_filename is None:
		f, psd = f_PSD_from_func(
			lalsimulation.SimNoisePSDKAGRA, fLow, fNyq, deltaF)
	else:
		f, psd = f_PSD_from_file(psd_filename,fLow,fNyq,deltaF)
	detector = metric.Detector(
		'C',
		lat_lon_ori_2_xarm(DegMinSec2Rad(1,36.,15.,0.),DegMinSec2Rad(1,137.,10.,48.),DegMinSec2Rad(1,20.,0.,0.)),
		lat_lon_ori_2_yarm(DegMinSec2Rad(1,36.,15.,0.),DegMinSec2Rad(1,137.,10.,48.),DegMinSec2Rad(1,20.,0.,0.)),
		lat_lon_2_vertex(DegMinSec2Rad(1,36.,15.,0.),DegMinSec2Rad(1,137.,10.,48.)),
		f=f,
		psd=psd
		)
	detector.set_required_moments()
	return detector

def make_IndIGO(fLow=10., fNyq=2048., deltaF=0.1, psd_filename=None):
	if psd_filename is None:
		f, psd = f_PSD_from_func(
			lalsimulation.SimNoisePSDaLIGOZeroDetHighPower, fLow, fNyq, deltaF)
	else:
		f, psd = f_PSD_from_file(psd_filename,fLow,fNyq,deltaF)
	detector = metric.Detector(
		'I',
		lat_lon_ori_2_xarm(DegMinSec2Rad(1,19.,5.,47.),DegMinSec2Rad(1,74.,2.,51.),DegMinSec2Rad(1,270.,0.,0.)),
		lat_lon_ori_2_yarm(DegMinSec2Rad(1,19.,5.,47.),DegMinSec2Rad(1,74.,2.,51.),DegMinSec2Rad(1,270.,0.,0.)),
		lat_lon_2_vertex(DegMinSec2Rad(1,19.,5.,47.),DegMinSec2Rad(1,74.,2.,51.)),
		f=f,
		psd=psd
		)
	detector.set_required_moments()
	return detector

