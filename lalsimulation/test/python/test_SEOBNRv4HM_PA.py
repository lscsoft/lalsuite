# -*- coding: utf-8 -*-
#
# Copyright (C) 2022 Serguei Ossokine and Deyan Mihaylov
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


import sys
import lal
import numpy as np
import pytest
import lalsimulation as lalsim

np.set_printoptions(precision=17)


def get_stored_data():
    """
    Pre-computed data for different waveform modes,
    in geometric units.
    Used in the regression test of SEOBNRv4HM_PA

    Returns:
        Dict: dictionary containing the modes
    """
    data = {}

    data[(2, 2)] = np.array(
        [
            0.04850856975870988 - 0.00637294698993718j,
            -0.04912612080712583 + 0.00655602913322001j,
            -0.05075251795235988 + 0.01393045443088556j,
            -0.05844555878687606 - 0.0025119285616761j,
            0.03895349499320805 + 0.2828489920874551j,
            -0.02451033205173241 + 0.18361572401407864j,
        ]
    )
    data[(2, 1)] = np.array(
        [
            0.00015134315993136 + 1.7876173028986764e-03j,
            0.00182388997146254 - 1.5706182019653632e-04j,
            -0.00198685243641532 + 3.1167553798632072e-04j,
            0.00237408798978243 - 1.0622079026797424e-05j,
            -0.01226851299029523 - 3.7258992649276464e-02j,
            0.03233503550956496 - 6.0741982156654298e-04j,
        ]
    )
    data[(3, 3)] = np.array(
        [
            -0.00080496742532625 - 0.00373617566967526j,
            0.00380605241800963 - 0.00083355145728985j,
            -0.00389075919837588 + 0.00173785170779303j,
            0.00498372964164358 + 0.00021365074856033j,
            -0.05898484227154988 - 0.00809194919440455j,
            0.04526307520945221 + 0.006154020171714j,
        ]
    )
    data[(4, 4)] = np.array(
        [
            -0.00063470216562911 + 1.7400501804824798e-04j,
            -0.00065050544045144 + 1.8127127646207350e-04j,
            -0.00065177365819477 + 3.9311670099763233e-04j,
            -0.00093674464079264 - 7.3097859220779265e-05j,
            0.01740509349942739 - 1.4258588336993408e-02j,
            0.01100446574069062 - 9.3188777674405466e-03j,
        ]
    )
    data[(5, 5)] = np.array(
        [
            3.3052984301931995e-05 + 9.6445484843000264e-05j,
            9.9436079870226238e-05 - 3.4639610751927051e-05j,
            -9.5674406767707075e-05 + 7.6254715602813841e-05j,
            1.5846049146157098e-04 + 1.6718394390263267e-05j,
            -9.6093981617729769e-03 + 2.3745387625395417e-03j,
            6.4808842912091073e-03 + 5.2195718745625512e-06j,
        ]
    )
    return data


def test_SEOBNRv4HMPA():
    """
    Test that the waveform modes of SEOBNRv4HM_PA have not changed.
    This is done by comparing the modes directly at a set of pre-determined
    time points that cover the inspiral, merger and RD.
    """
    # Note second last element is (2,2) mode peak
    indices = [0, 10000, 50000, 100000, 187006, 187100]
    # Physical params
    m1 = 60.0
    m2 = 20.0
    chi1_x, chi1_y, chi1_z = 0.0, 0.0, 0.1
    chi2_x, chi2_y, chi2_z = 0.0, 0.0, 0.3
    modes = [(2, 2), (2, 1), (3, 3), (4, 4), (5, 5)]

    phi_c = 0.3
    f_start22 = 7.3  # Frequency of the 22 mode at which the signal starts
    distance = 1  # Irrelevant

    deltaT = 1.0 / 16384.0
    approx = lalsim.SEOBNRv4HM_PA
    lmax = 5
    lal_params = lal.CreateDict()
    # Generate the modes
    hlm = lalsim.SimInspiralChooseTDModes(
        phi_c,
        deltaT,
        m1 * lal.MSUN_SI,
        m2 * lal.MSUN_SI,
        chi1_x,
        chi1_y,
        chi1_z,
        chi2_x,
        chi2_y,
        chi2_z,
        f_start22,
        f_start22,
        distance,
        lal_params,
        lmax,
        approx,
    )
    time_hI = lalsim.SphHarmTimeSeriesGetMode(hlm, 2, 2).deltaT * np.arange(
        len(lalsim.SphHarmTimeSeriesGetMode(hlm, 2, 2).data.data)
    )
    hI = {}

    for lm in modes:
        hI[lm] = lalsim.SphHarmTimeSeriesGetMode(hlm, lm[0], lm[1]).data.data

    # Compare stored to current data
    data = get_stored_data()
    current_data = {}
    for mode in modes:
        current_data[mode] = hI[mode][indices] / (
            (m1 + m2) * lal.MRSUN_SI
        )  # rescale to geometric units
        # Make sure data is the same
        np.testing.assert_allclose(data[mode], current_data[mode], atol=0, rtol=1e-13)


if __name__ == "__main__":
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-SEOBNRv4HMPA.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
