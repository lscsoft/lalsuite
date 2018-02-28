# Copyright (C) 2009  Kipp Cannon
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
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


"""
This module provides some utilities to assist with FFTs
"""


import lal


import git_version


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                  Utilities
#
# =============================================================================
#


def prepare_fseries_for_real8tseries(series):
	"""
	Construct a COMPLEX16FrequencySeries object suitable for storing
	the Fourier transform of a REAL8TimeSeries object.
	"""
	n = len(series.data.data)
	return lal.CreateCOMPLEX16FrequencySeries(
		name = series.name,
		epoch = series.epoch,
		f0 = series.f0,	# note: non-zero f0 not supported by LAL
		deltaF = 1.0 / (n * series.deltaT),
		sampleUnits = series.sampleUnits * lal.SecondUnit,
		length = n // 2 + 1
	)


def prepare_fseries_for_complex16tseries(series):
	"""
	Construct a COMPLEX16FrequencySeries object suitable for storing
	the Fourier transform of a COMPLEX16TimeSeries object.
	"""
	n = len(series.data.data)
	return lal.CreateCOMPLEX16FrequencySeries(
		name = series.name,
		epoch = series.epoch,
		f0 = series.f0,	# note: non-zero f0 not supported by LAL
		deltaF = 1.0 / (n * series.deltaT),
		sampleUnits = series.sampleUnits * lal.SecondUnit,
		length = n
	)
