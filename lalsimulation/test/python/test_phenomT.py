# -*- coding: utf-8 -*-
#
# Copyright (C) 2020 Hector Estelles (adapted from test_phenomX.py of Cecilio García Quiros)
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

#Simple test to see if the PhenomT family models have changed: IMRPhenomT, IMRPhenomTHM, IMRPhenomTP and IMRPhenomTPHM.

import sys
import pytest
import lal
import lalsimulation
import numpy as np

# -- utility functions ---------------------

def get_amp_phase(h):
    amp = np.abs(h)
    phase = np.unwrap(np.angle(h))
    return amp, phase

def sum_sqr_diff(x, y):
    return np.sqrt( np.sum( (x-y)**2 )  )

def gen_test_data(spin1x, approximant, mode_array, PV, FS):
    """
    compute the difference between two waveforms
    and compare to expected value
    """
    lalparams = lal.CreateDict()

    if(mode_array!=None):
        ModeArray = lalsimulation.SimInspiralCreateModeArray()
        for mode in mode_array:
            lalsimulation.SimInspiralModeArrayActivateMode(ModeArray, mode[0], mode[1])
        lalsimulation.SimInspiralWaveformParamsInsertModeArray(lalparams, ModeArray)

    if PV!=None:
        lalsimulation.SimInspiralWaveformParamsInsertPhenomXPrecVersion(lalparams, PV)
    if FS!=None:
        lalsimulation.SimInspiralWaveformParamsInsertPhenomXPFinalSpinMod(lalparams, FS)

    common_pars=dict(
    m1=50*lal.MSUN_SI,
    m2=30*lal.MSUN_SI,
    s1x=spin1x,
    s1y=0.,
    s1z=0.,
    s2x=0.,
    s2y=0.,
    s2z=0.,
    distance=1,
    inclination=np.pi/3.,
    phiRef=0.,
    longAscNodes=0.,
    eccentricity=0.,
    meanPerAno=0.,
    deltaT=1./4096.,
    f_min=30.,
    f_ref=30.,
    params=lalparams,
    approximant=approximant
    )

    pars1=common_pars.copy()

    pars2=common_pars.copy()
    pars2.update({"inclination":0.17})
    pars2.update({"phiRef":0.5})

    hp1, hc1 = lalsimulation.SimInspiralChooseTDWaveform(**pars1)
    hp2, hc2 = lalsimulation.SimInspiralChooseTDWaveform(**pars2)

    # compute amp and phase
    hp1_amp, hp1_phase = get_amp_phase(hp1.data.data)
    hc1_amp, hc1_phase = get_amp_phase(hc1.data.data)

    hp2_amp, hp2_phase = get_amp_phase(hp2.data.data)
    hc2_amp, hc2_phase = get_amp_phase(hc2.data.data)

    hp_amp_diff = sum_sqr_diff(hp1_amp, hp2_amp)
    hp_phase_diff = sum_sqr_diff(hp1_phase, hp2_phase)

    hc_amp_diff = sum_sqr_diff(hc1_amp, hc2_amp)
    hc_phase_diff = sum_sqr_diff(hc1_phase, hc2_phase)

    return hp_amp_diff, hp_phase_diff, hc_amp_diff, hc_phase_diff


# -- test functions ---------------------

def test_IMRPhenomT():
    """
    This test checks that IMRPhenomT hasn't changed.
    It does this by generating two PhenomT waveforms and computing
    their difference (according to their amplitude and phases)
    and compares them to pre-computed values.

    these pre-computed values were computed using the following line:

    `expected_result  =  np.array(gen_test_data(0., lalsimulation.IMRPhenomT))`

    """

    expected_result = np.array([1.93196330e+05, 6.75990151e+01, 2.02364224e+05, 6.65691548e+01])

    new_result  =  np.array(gen_test_data(0., lalsimulation.IMRPhenomT, None, None, None))

    np.testing.assert_allclose(new_result, expected_result, rtol=1e-6, err_msg="IMRPhenomT test failed")

def test_IMRPhenomTHM():
    """
    This test checks that IMRPhenomTHM hasn't changed.
    It does this by generating two PhenomTHM waveforms and computing
    their difference (according to their amplitude and phases)
    and compares them to pre-computed values.

    these pre-computed values were computed using the following line:

    `expected_result  =  np.array(gen_test_data(0., lalsimulation.IMRPhenomTHM))`

    """

    expected_result = np.array([1.93318947e+05, 6.75990151e+01, 2.02583328e+05, 6.64207284e+01])

    new_result  =  np.array(gen_test_data(0., lalsimulation.IMRPhenomTHM, None, None, None))

    np.testing.assert_allclose(new_result, expected_result, rtol=1e-6, err_msg="IMRPhenomTHM test failed")

def test_IMRPhenomTP():
    """
    This test checks that IMRPhenomTPHM hasn't changed.
    It does this by generating two PhenomTHM waveforms and computing
    their difference (according to their amplitude and phases)
    and compares them to pre-computed values.

    these pre-computed values were computed using the following line:

    `expected_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomTP, None, PV, FS))`

    """

    stored_resultsTP={None: {None: (182073.7553633626, 65.14545181720996, 177303.3661438229, 55.49159011183478)}}

    PVs = stored_resultsTP.keys()

    for PV in PVs:
        FSs = stored_resultsTP[PV].keys()

        for FS in FSs:
            expected_result = stored_resultsTP[PV][FS]

            new_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomTP, None, PV, FS))

            np.testing.assert_allclose(new_result, expected_result, rtol=2e-3, err_msg="IMRPhenomTP test failed")


def test_IMRPhenomTPHM():
    """
    This test checks that IMRPhenomTPHM hasn't changed.
    It does this by generating two PhenomTHM waveforms and computing
    their difference (according to their amplitude and phases)
    and compares them to pre-computed values.

    these pre-computed values were computed using the following line:

    `expected_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomTPHM, None, PV, FS))`

    """

    stored_resultsTPHM={None: {None: (182535.7160204654, 64.68935203000481, 178021.9819382646, 56.11063895508146)}}

    PVs = stored_resultsTPHM.keys()

    for PV in PVs:
        FSs = stored_resultsTPHM[PV].keys()

        for FS in FSs:
            expected_result = stored_resultsTPHM[PV][FS]

            new_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomTPHM, None, PV, FS))

            np.testing.assert_allclose(new_result, expected_result, rtol=1e-6, err_msg="IMRPhenomTPHM test failed")


if __name__ == '__main__':
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-phenomT.xml"]
    sys.exit(pytest.main(args=[__file__] + args))

