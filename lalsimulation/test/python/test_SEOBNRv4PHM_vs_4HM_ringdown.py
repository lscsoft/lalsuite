# -*- coding: utf-8 -*-
#
# Copyright (C) 2020 Serguei Ossokine
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http: //www.gnu.org/licenses/>.

"""Test that the aligned-spin limit of SEOBNRv4PHM is the same as SEOBNRv4HM
In particular we consider cases where the final spin is *negative* with respect to the inital
 orbital angular momentum. """

import sys
import pytest
from scipy.interpolate import InterpolatedUnivariateSpline
import lalsimulation as ls
import lal
import numpy as np

# -- test functions ---------------------


def get_SEOBNRv4HM_modes(q, M, chi1, chi2, f_start, distance, deltaT):
    """Generate SEOBNRv4HM modes"""
    m1SI = lal.MSUN_SI * q * M / (1.0 + q)
    m2SI = lal.MSUN_SI * M / (1.0 + q)
    nqcCoeffsInput = lal.CreateREAL8Vector(10)
    sphtseries, dyn, dynHI = ls.SimIMRSpinAlignedEOBModes(
        deltaT,
        m1SI,
        m2SI,
        f_start,
        distance,
        chi1,
        chi2,
        41,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
        1.0,
        nqcCoeffsInput,
        0,
    )

    # The minus sign in front of the modes takes into account the fact that the polarization basis in EOB
    # conventions is different wrt the one in LAL
    hI = {}
    modes = [(2, 2), (2, 1), (3, 3), (4, 4), (5, 5)]
    for lm in modes:
        hI[lm] = np.trim_zeros(
            -1 * ls.SphHarmTimeSeriesGetMode(sphtseries, lm[0], lm[1]).data.data, "b"
        )

    return hI


def get_SEOBNRv4PHM_modes(q, M, chi1, chi2, f_start, distance, deltaT):
    """Generate SEOBNRv4PHM modes"""

    m1SI = lal.MSUN_SI * q * M / (1.0 + q)
    m2SI = lal.MSUN_SI * M / (1.0 + q)
    approx = ls.SEOBNRv4PHM
    hlm = ls.SimInspiralChooseTDModes(0.,
        deltaT,
        m1SI,
        m2SI,
        chi1[0],
        chi1[1],
        chi1[2],
        chi2[0],
        chi2[1],
        chi2[2],
        f_start,
        f_start,
        distance,
        None,
        5,
        approx,
    )
    hI = {}
    modes = [(2, 2), (2, 1), (3, 3), (4, 4), (5, 5)]
    for lm in modes:
        hI[lm] = ls.SphHarmTimeSeriesGetMode(hlm, lm[0], lm[1]).data.data

    return hI


@pytest.mark.parametrize("q", [5.0, 10.0, 15.0, 20.0])
def test_aligned_spin_limit_ringdown(q):
    """ Test that the (2,2) mode frequency in the ringdown is within in 1Hz between v4HM and v4PHM.
    In particular, consider cases where the final spin is negative with respect to the initial
    orbital angular momentum. Of course considers only aligned-spin cases"""
    # Step 0: set parameters

    chi1_x, chi1_y, chi1_z = 0.0, 0.0, -0.95
    chi2_x, chi2_y, chi2_z = 0.0, 0.0, -0.95
    m_total = 150  # Feducial total mass in solar masses
    deltaT = 1.0 / 16384
    distance = 500 * 1e6 * lal.PC_SI  # 500 Mpc in m
    f_start22 = 20

    # Step 1: Get the SEOBNRv4PHM modes
    hlmP = get_SEOBNRv4PHM_modes(
        q,
        m_total,
        [chi1_x, chi2_y, chi1_z],
        [chi1_x, chi2_y, chi2_z],
        f_start22 * 0.5,
        distance,
        deltaT,
    )

    # Step 2: Get the SEOBNRv4HM modes
    hlmA = get_SEOBNRv4HM_modes(
        q, m_total, chi1_z, chi2_z, f_start22 * 0.5, distance, deltaT
    )

    time = deltaT * np.arange(len(hlmA[(2, 2)]))
    time2 = deltaT * np.arange(len(hlmP[(2, 2)]))
    # Compute the (2,2) mode frequencies
    mode = (2, 2)

    phaseA = np.unwrap(np.angle(hlmA[mode]))
    phaseP = np.unwrap(np.angle(hlmP[mode]))

    intrpA = InterpolatedUnivariateSpline(time, phaseA)
    intrpP = InterpolatedUnivariateSpline(time2, phaseP)

    omega_A = intrpA.derivative()
    omega_P = intrpP.derivative()

    # Consider the frequency in the ringdown
    amp = np.abs(hlmA[mode].data)
    am = np.argmax(amp)

    omega_A_rd = omega_A(time[am : am + 350])
    omega_P_rd = omega_P(time[am : am + 350])

    # Ensure the difference is small. Note here atol=1 corresponds to 1 Hz difference on frequncies of hundreds of Hz. These small
    # diffrences are expected due to some minor code differences.
    # Before a bug was fixed the difference would be huge, i.e. ~100s Hz
    assert np.allclose(
        omega_A_rd, omega_P_rd, atol=1
    ), "The frequencies of the (2,2) mode don't agree between SEOBNRv4PHM and SEOBNRv4HM!"


# -- run the tests ------------------------------

if __name__ == "__main__":
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-SEOBNRv4PHM_vs_4HM_ringdown.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
