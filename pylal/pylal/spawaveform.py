# Copyright (C) 2010  Kipp Cannon,  Drew Keppel
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
This module is a wrapper of the _spawaveform module, supplementing the C
code in that module with additional features that are more easily
implemented in Python.  It is recommended that you import this module
rather than importing _spawaveform directly.
"""


import math


import lal
from pylal import git_version
from pylal._spawaveform import *


__author__ = "Kipp Cannon <kipp.cannon@ligo.org>"
__version__ = "git id %s" % git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                    Empty
#
# =============================================================================
#


def imrchirptime(m1, m2, fLower, chi, a_hat = 0.98, e_folds = 10):
	"""
	An approximate IMR chirptime in seconds.

	FIXME this should be replaced by something better

	1) compute the SPA chirptime up to ringdown frequency at 3.5 PN, verify that it is nonnegative
	2) then add efolds worth of ringdown time

	Ringdown decay time forumla in solar masses is:

	tau = 2 * (m1+m2) * 5e-6 * (0.7 + 1.4187 * (1-a_hat)**-0.4990) / (1.5251 - 1.1568 * (1-a_hat)**0.1292)

	from (7-9) of LIGO-P1300156.  

	@param m1 Mass 1
	@param m2 Mass 2
	@param fLower the starting frequency
	@param chi the effective spin parameter from computechi()
	@param e_folds The number of efolds to use in the ringdown signal duration, default 10
	@param a_hat The dimensionless spin of the final black hole, default 0.98
	"""

	assert (a_hat < 0.9999999999999999) # demand spin less than 1 (or approximately the closest floating point representation of 1)
	fFinal = imrffinal(m1, m2, chi, 'ringdown')
	assert (fFinal > fLower) # demand that the low frequency comes before the ringdown frequency
	tau = 2 * (m1+m2) * 5e-6 * (0.7 + 1.4187 * (1-a_hat)**-0.4990) / (1.5251 - 1.1568 * (1-a_hat)**0.1292)
	inspiral_time = chirptime(m1, m2, 7, fLower, fFinal, chi)
	if inspiral_time < 0:
		raise ValueError("Inspiral time is negative: m1 = %e, m2 = %e, flow = %e, chi = %e" % (m1, m2, fLower, chi)) # demand positive inspiral times
	return inspiral_time + e_folds * tau


def eta(m1, m2):
	"""
	Compute the symmetric mass ratio, eta.
	"""
	return m1*m2/(m1+m2)**2.


def chirpmass(m1, m2):
	"""
	Compute the chirp mass in seconds.
	"""
	return lal.MTSUN_SI * (m1+m2) * eta(m1, m2)**.6


def ms2taus(m1, m2, f0 = 40.0):
	"""
	Solve for tau_0 and tau_3 from m1 and m2.
	"""
	tau0 = 5./256./(math.pi*f0)**(8./3.) * chirpmass(m1,m2)**(-5./3.)
	tau3 = math.pi/8./eta(m1,m2)**.6/(math.pi*f0)**(5./3.) * chirpmass(m1,m2)**(-2./3.)
	return tau0, tau3


def taus2ms(tau0, tau3, f0 = 40.0):
	"""
	Solve for m1 and m2 from tau_0 and tau_3.
	"""
	Mc = (5./256./(math.pi*f0)**(8./3.) / tau0)**(3./5.)
	eta = (math.pi/8./(math.pi*f0)**(5./3.) / tau3 / Mc**(2./3.))**(5./3.)

	M = Mc / eta**(3./5.)

	m1 = (1. + abs(1. - 4.*eta)**.5) * M / 2.
	m2 = M - m1

	return m1 / lal.MTSUN_SI, m2 / lal.MTSUN_SI
