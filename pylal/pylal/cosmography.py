# Copyright (C) 2011 Nickolas Fotopoulos
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
"""
Cosmography based on Hogg 1998 (astro-ph/9905116)
Cosmological parameters taken from WMAP 7 year BAO + H0 ML value (arXiv:1001.4538)
"""

from numpy import pi, sinh, arcsinh
from scipy.integrate import quad
from scipy.optimize import newton

from lal import PC_SI as LAL_PC_SI
from lal import C_SI as LAL_C_SI

h0 = 0.703  # H0 = h0 * 100 km / s / Mpc
H0_SI = h0 * 100000 / (1000000 * LAL_PC_SI)  # Hubble constant in inverse seconds
DH_SI = LAL_C_SI / H0_SI  # Hubble distance in meters
OmegaL = 0.729
OmegaM = 0.271 # sum of baryonic and cold dark matter Omega
OmegaR = 0 # the remainder is consistent with zero given the error bars

def Einv(z):
    """
    Integrand used in many cosmography integrals. 1 / E(z) in Hogg's notation.
    """
    return (OmegaM * (1 + z)**3 + OmegaR * (1 + z)**2 + OmegaL)**-0.5

def DM(z):
    """
    Comoving distance (transverse) at z. Hard-coded for OmegaR = 0.
    """
    return DH_SI * quad(Einv, 0, z)[0]

def DL(z):
    """
    Luminosity distance
    """
    return DM(z) * (1 + z)

def compute_redshift(DL0):
    """
    Compute the redshift from a luminosity distance.
    Use the Newton-Raphson refinement method on a first-order guess.
    """
    return newton(lambda z: DL(z) - DL0, DL0 / DH_SI)

def dVdz(z):
    """
    Different volume element per unit redshift. Hard-coded for OmegaR = 0.
    """
    return DH_SI * DM(z)**2 * Einv(z) * 4 * pi

def V(z):
    """
    Analytical integration of dVdz. Hard-coded for OmegaR = 0.
    Double precision craps out below about 100 kpc (z=2.3e-6).
    """
    tmpDM = DM(z)
    return 4 * pi / 3 * tmpDM * tmpDM * tmpDM

def V_from_DL_z(D, z=None):
    """
    Analytical integration of dVdz. Hard-coded for OmegaR = 0.
    Sped up for the case that you have D rather than (or in addition to) z.
    Double precision craps out below about 100 kpc (z=2.3e-6).
    """
    tmpDM = D / (1 + (z or compute_redshift(D)))
    return 4 * pi / 3 * tmpDM * tmpDM * tmpDM
