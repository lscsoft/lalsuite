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
    data = {
        (2, 2): [
            [
                0.40689716589895747,
                0.43309591538465037,
                0.4721748343322336,
                0.5447370783113126,
                2.0023458851711906,
                1.8235613978215368,
            ],
            [
                -1.3969263518385833e-01,
                -1.5925496591368452e02,
                -3.3713851919667729e02,
                -5.4743337706169893e02,
                -8.7302404936234018e02,
                -8.7492784847230996e02,
            ],
        ],
        (2, 1): [
            [
                0.0313225894843847,
                0.0347815707120482,
                0.04022820240989489,
                0.05123839474459091,
                0.7108940042712576,
                0.8451123158705316,
            ],
            [
                1.4800894807662368,
                -78.07974926209056,
                -167.0249793662184,
                -272.1793580797458,
                -435.7330365096353,
                -437.18207987235945,
            ],
        ],
        (3, 3): [
            [
                0.05111250481525122,
                0.05608646181014196,
                0.0637759095572762,
                0.0788643477891161,
                0.5769500315793753,
                0.5669386351414297,
            ],
            [
                -1.7986372334622316,
                -240.47335889154334,
                -507.3014845029567,
                -822.7492254842681,
                -1311.313003009921,
                -1314.2554619594298,
            ],
        ],
        (4, 4): [
            [
                0.0094732179486535,
                0.01072105034372493,
                0.01272251602463651,
                0.01688049604193659,
                0.24798290965344766,
                0.23257308477746422,
            ],
            [
                2.84188245155912,
                -315.39054175943875,
                -671.1605084583525,
                -1091.7556790940084,
                -1743.009269177671,
                -1746.8685414830832,
            ],
        ],
        (5, 5): [
            [
                0.00189622618209567,
                0.00221381745387731,
                0.00274267447438163,
                0.00390737715569474,
                0.11379890992679738,
                0.10803857266725103,
            ],
            [
                1.1933141845843078e00,
                -3.9659776027206124e02,
                -8.4131116974703161e02,
                -1.3670574321591489e03,
                -2.1821066712187244e03,
                -2.1870969032147736e03,
            ],
        ],
    }
    return data


rtol_amp = {(2, 2): 1e-5, (2, 1): 1e-5, (3, 3): 1e-5, (4, 4): 1e-5, (5, 5): 1e-5}
rtol_phase = {(2, 2): 1e-5, (2, 1): 1e-2, (3, 3): 1e-5, (4, 4): 1e-3, (5, 5): 1e-3}


def test_SEOBNRv4HMPA():
    """
    Test that the waveform modes of SEOBNRv4HM_PA have not changed.
    This is done by comparing the modes directly at a set of pre-determined
    time points that cover the inspiral, merger and RD.
    """
    # Note second last element is (2,2) mode peak

    # Physical params

    chi1_x, chi1_y = 0.0, 0.0
    chi2_x, chi2_y = 0.0, 0.0
    m1 = 8.16594601e01
    m2 = 1.08125863e01
    chi1_z = -7.19166988e-01
    chi2_z = 5.67932765e-01
    modes = [(2, 2), (2, 1), (3, 3), (4, 4), (5, 5)]

    # Note second last element is (2,2) mode peak
    indices = [0, 54291, 108582, 162873, 217164, 217214]
    phi_c = 0.3
    f_start22 = 7.3  # Frequency of the 22 mode at which the signal starts
    distance = 10000  # Irrelevant

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

    for mode in modes:
        amp = np.abs(hI[mode])[indices]
        phase = np.unwrap(np.angle(hI[mode]))[indices]
        amp_stored, phase_stored = data[mode]
        np.testing.assert_allclose(amp, amp_stored, atol=0, rtol=rtol_amp[mode])
        np.testing.assert_allclose(phase, phase_stored, atol=0, rtol=rtol_phase[mode])


if __name__ == "__main__":
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-SEOBNRv4HMPA.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
