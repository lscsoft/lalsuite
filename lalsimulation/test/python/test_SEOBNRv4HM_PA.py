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
                # 0.40689716589895747,
                # 0.43309591538465037,
                # 0.4721748343322336,
                # 0.5447370783113126,
                # 2.0023458851711906,
                # 1.8235613978215368,
                0.406897143874587,
                0.43309588629164864,
                0.4721747898690362,
                0.5447369805430688,
                2.0023481970562735,
                1.8245536590453957,
            ],
            [
                # -1.3969263518385833e-01,
                # -1.5925496591368452e02,
                # -3.3713851919667729e02,
                # -5.4743337706169893e02,
                # -8.7302404936234018e02,
                # -8.7492784847230996e02,
                -1.3968438162982541e-01,
                -1.5925495580678407e+02,
                -3.3713850517954216e+02,
                -5.4743335258722220e+02,
                -8.7301946486216650e+02,
                -8.7492212355961533e+02,
            ],
        ],
        (2, 1): [
            [
                # 0.0313225894843847,
                # 0.0347815707120482,
                # 0.04022820240989489,
                # 0.05123839474459091,
                # 0.7108940042712576,
                # 0.8451123158705316,
                0.031322421516861723,
                0.034781337721511429,
                0.040227835734167786,
                0.051237617991819341,
                0.710259630940248132,
                0.844879838309898212,
            ],
            [
                # 1.4800894807662368,
                # -78.07974926209056,
                # -167.0249793662184,
                # -272.1793580797458,
                # -435.7330365096353,
                # -437.18207987235945,
                1.480096232806204837e+00,
                -7.807974106262993530e+01,
                -1.670249682972653318e+02,
                -2.721793395565951528e+02,
                -4.357288791143605522e+02,
                -4.371762545684100587e+02,
            ],
        ],
        (3, 3): [
            [
                # 0.05111250481525122,
                # 0.05608646181014196,
                # 0.0637759095572762,
                # 0.0788643477891161,
                # 0.5769500315793753,
                # 0.5669386351414297,
                0.051112501745226545,
                0.056086457627118551,
                0.063775902886042721,
                0.078864332055436248,
                0.576776486652450626,
                0.567171822083015198,
            ],
            [
                # -1.7986372334622316,
                # -240.47335889154334,
                # -507.3014845029567,
                # -822.7492254842681,
                # -1311.313003009921,
                # -1314.2554619594298,
                -1.798624864004767732e+00,
                -2.404733437441886679e+02,
                -5.073014634939491430e+02,
                -8.227491887980770571e+02,
                -1.311305924173796029e+03,
                -1.314246406233119387e+03
            ],
        ],
        (4, 4): [
            [
                # 0.0094732179486535,
                # 0.01072105034372493,
                # 0.01272251602463651,
                # 0.01688049604193659,
                # 0.24798290965344766,
                # 0.23257308477746422,
                0.009473217127067403,
                0.010721049156921609,
                0.012722513984742818,
                0.016880490682204544,
                0.247871444279401887,
                0.232676698120316383,
            ],
            [
                # 2.84188245155912,
                # -315.39054175943875,
                # -671.1605084583525,
                # -1091.7556790940084,
                # -1743.009269177671,
                # -1746.8685414830832,
                2.841897673213301445e+00,
                -3.153905230860205506e+02,
                -6.711604824120537387e+02,
                -1.091755633222126562e+03,
                -1.743000212109299355e+03,
                -1.746856853681002576e+03,
            ],
        ],
        (5, 5): [
            [
                # 0.00189622618209567,
                # 0.00221381745387731,
                # 0.00274267447438163,
                # 0.00390737715569474,
                # 0.11379890992679738,
                # 0.10803857266725103,
                0.001896205803077448,
                0.002213789075042580,
                0.002742629541378499,
                0.003907280821832254,
                0.113725569203534563,
                0.108064932141790415,
            ],
            [
                # 1.1933141845843078e00,
                # -3.9659776027206124e02,
                # -8.4131116974703161e02,
                # -1.3670574321591489e03,
                # -2.1821066712187244e03,
                # -2.1870969032147736e03,
                1.193337715232760221e+00,
                -3.965977315334420723e+02,
                -8.413111302239032057e+02,
                -1.367057364036911622e+03,
                -2.182095223636347328e+03,
                -2.187082823278491560e+03,
            ],
        ],
    }
    return data


rtol_amp = {(2, 2): 1e-5, (2, 1): 1e-5, (3, 3): 1e-5, (4, 4): 1e-5, (5, 5): 1e-5}
rtol_phase = {(2, 2): 1e-5, (2, 1): 1e-5, (3, 3): 1e-5, (4, 4): 1e-5, (5, 5): 1e-5}


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
        # print(',\n'.join([f"{x:.18f}" for x in amp]))
        # print(',\n'.join([f"{x:.18e}" for x in phase]))
        np.testing.assert_allclose(amp, amp_stored, atol=0, rtol=rtol_amp[mode])
        np.testing.assert_allclose(phase, phase_stored, atol=0, rtol=rtol_phase[mode])


if __name__ == "__main__":
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-SEOBNRv4HMPA.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
