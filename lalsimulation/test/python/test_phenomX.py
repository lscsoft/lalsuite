# -*- coding: utf-8 -*-
#
# Copyright (C) 2020 Cecilio García Quirós (adapted from test_phenomPv3HM.py of Sebastian Khan)
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

"""Simple test to see if the PhenomX family models have changed: IMRPhenomXAS, IMRPhenomXHM, IMRPhenomXP, IMRPhenomXPHM, IMRPhenomXAS_NRTidalv2, IMRPhenomXP_NRTidalv2, and IMRPhenomXO4a.
"""

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


def gen_test_data(spin1x, approximant, mode_array,  lalparams = None, type='bbh', generic_spins=False,deltaF=1./4.,release=None):
    """
    compute the difference between two waveforms
    and compare to expected value

    type: should either be 'bbh' or 'bns'
    generic_spins: Use generic spins (with in-plane spins scaled from spin1x, which must be at most ~0.69, and aligned spins hard-coded) if True;
                   otherwise the only potentially nonzero spin component is spin1x

    """

    if(lalparams == None):
        lalparams = lal.CreateDict()

    if(mode_array!=None):
        ModeArray = lalsimulation.SimInspiralCreateModeArray()
        for mode in mode_array:
            lalsimulation.SimInspiralModeArrayActivateMode(ModeArray, mode[0], mode[1])
        lalsimulation.SimInspiralWaveformParamsInsertModeArray(lalparams, ModeArray)

    if release is not None:
        lalsimulation.SimInspiralWaveformParamsInsertPhenomXHMReleaseVersion(lalparams, release)

    if type == 'bbh':
        m1_sel = 50*lal.MSUN_SI
        m2_sel = 30*lal.MSUN_SI
        m2_sel_prime = 20.*lal.MSUN_SI
        f_max_sel = 512.
    elif type == 'bns':
        m1_sel = 1.6*lal.MSUN_SI
        m2_sel = 1.4*lal.MSUN_SI
        m2_sel_prime = 1.3*lal.MSUN_SI
        f_max_sel = 2048.

        lalsimulation.SimInspiralWaveformParamsInsertTidalLambda1(lalparams, 200.)
        lalsimulation.SimInspiralWaveformParamsInsertTidalLambda2(lalparams, 300.)
    else:
        raise ValueError("Unknown binary type")

    if generic_spins:
        spin1y = 0.5*spin1x
        spin1z = -0.3
        spin2x = -0.9*spin1x
        spin2y = 1.2*spin1x
        spin2z = 0.2
    else:
        spin1y = 0.
        spin1z = 0.
        spin2x = 0.
        spin2y = 0.
        spin2z = 0.

    common_pars=dict(
    m1=m1_sel,
    m2=m2_sel,
    S1x=spin1x,
    S1y=spin1y,
    S1z=spin1z,
    S2x=spin2x,
    S2y=spin2y,
    S2z=spin2z,
    distance=1,
    inclination=np.pi/3.,
    phiRef=0.,
    longAscNodes=0.,
    eccentricity=0.,
    meanPerAno=0.,
    deltaF=deltaF,
    f_min=30.,
    f_max=f_max_sel,
    f_ref=30.,
    LALpars=lalparams,
    approximant=approximant
    )

    pars1=common_pars.copy()

    pars2=common_pars.copy()
    pars2.update({"m2":m2_sel_prime})

    hp1, hc1 = lalsimulation.SimInspiralChooseFDWaveform(**pars1)
    hp2, hc2 = lalsimulation.SimInspiralChooseFDWaveform(**pars2)

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

def test_IMRPhenomXAS():
    """
    This test checks that IMRPhenomXAS hasn't changed.
    It does this by generating two PhenomXAS waveforms and computing
    their difference (according to their amplitude and phases)
    and compares them to pre-computed values.

    these pre-computed values were computed using the following line:

    `expected_result  =  np.array(gen_test_data(0., lalsimulation.IMRPhenomXAS, None, generic_spins=True))`

    """

    expected_result = np.array([787.00663452, 165.89210208, 629.60530761, 165.89210208])

    new_result  =  np.array(gen_test_data(0., lalsimulation.IMRPhenomXAS, None, generic_spins=True))

    np.testing.assert_allclose(new_result, expected_result, rtol=1e-6, err_msg="IMRPhenomXAS test failed")

def test_IMRPhenomXHM():
    """
    This test checks that IMRPhenomXHM hasn't changed.
    It does this by generating two PhenomXHM waveforms and computing
    their difference (according to their amplitude and phases)
    and compares them to pre-computed values.

    these pre-computed values were computed using the following line:

    `expected_result  =  np.array(gen_test_data(0., lalsimulation.IMRPhenomXHM, mode_array))`

    where mode_array is [[2,2],[2,-2],[2,1],[2,-1],[3,3],[3,-3],[4,4],[4,-4]] or [[3,2],[3,-2]]

    """


    expected_result = np.array([1005.16009183, 169.88197475, 768.18401876, 169.12711241])

    new_result  =  np.array(gen_test_data(0., lalsimulation.IMRPhenomXHM, [[2,2],[2,-2],[2,1],[2,-1],[3,3],[3,-3],[4,4],[4,-4]], release=122019))

    np.testing.assert_allclose(new_result, expected_result, rtol=1e-6, err_msg="IMRPhenomXHM (122019) no 32 mode test failed")

    expected_result = np.array([ 32.17818789, 216.01992794,   4.02227349, 215.97103911])

    new_result  =  np.array(gen_test_data(0., lalsimulation.IMRPhenomXHM, [[3,2],[3,-2]], release=122019))

    # rtol with 32 mode needs to b more lenient
    np.testing.assert_allclose(new_result, expected_result, rtol=3e-4, err_msg="IMRPhenomXHM (122019) 32 mode test failed")

    expected_result = np.array([1005.01319319, 169.88945372, 768.34648494, 169.13261004])

    new_result  =  np.array(gen_test_data(0., lalsimulation.IMRPhenomXHM, [[2,2],[2,-2],[2,1],[2,-1],[3,3],[3,-3],[4,4],[4,-4]]))

    np.testing.assert_allclose(new_result, expected_result, rtol=1e-6, err_msg="IMRPhenomXHM no 32 mode test failed")

    expected_result = np.array([34.62153262, 218.09073730, 4.32769157, 218.09073730])

    new_result  =  np.array(gen_test_data(0., lalsimulation.IMRPhenomXHM, [[3,2],[3,-2]]))

    # rtol with 32 mode needs to b more lenient
    np.testing.assert_allclose(new_result, expected_result, rtol=3.1e-4, err_msg="IMRPhenomXHM 32 mode test failed")

def test_IMRPhenomXP():
    """
    This test checks that IMRPhenomXP hasn't changed.
    It does this by generating two PhenomXP waveforms and computing
    their difference (according to their amplitude and phases)
    and compares them to pre-computed values.

    these pre-computed values were computed using the following line:

    `expected_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXP, None, generic_spins=False))`

    """
    lalDict = lal.CreateDict()
    lalsimulation.SimInspiralWaveformParamsInsertPhenomXPrecVersion(lalDict, 223)
    lalsimulation.SimInspiralWaveformParamsInsertPhenomXPFinalSpinMod(lalDict, 2)

    expected_result = np.array([1070.22089507,  271.62095671,  533.04876588,  268.4520146])

    new_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXP, None, generic_spins=False))

    np.testing.assert_allclose(new_result, expected_result, rtol=1e-6, err_msg="IMRPhenomXP test failed")


def test_IMRPhenomXP_NNLO():
    """
    This test checks that IMRPhenomXP with the NNLO precession option (version 102) hasn't changed.
    It does this by generating two PhenomXP waveforms and computing
    their difference (according to their amplitude and phases)
    and compares them to pre-computed values.

    these pre-computed values were computed using the following code:

    `lalDict = lal.CreateDict()`

    `lalsimulation.SimInspiralWaveformParamsInsertPhenomXPrecVersion(lalDict, 102)`

    `expected_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXP, None, lalparams=lalDict, generic_spins=False))`

    """

    lalDict = lal.CreateDict()
    lalsimulation.SimInspiralWaveformParamsInsertPhenomXPrecVersion(lalDict, 102)

    expected_result = np.array([1235.47048998, 226.22617618, 1049.4091208, 225.46870308])

    new_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXP, None, lalparams=lalDict, generic_spins=True))

    np.testing.assert_allclose(new_result, expected_result, rtol=1e-6, err_msg="IMRPhenomXP_NNLO test failed")



def test_IMRPhenomXP_MB():
    """
    This test checks that IMRPhenomXP hasn't changed.
    It does this by generating two PhenomXP waveforms and computing
    their difference (according to their amplitude and phases)
    and compares them to pre-computed values.

    these pre-computed values were computed using the following line:

    `expected_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXP, None, generic_spins=True))`

    FIXME: This description needs to be corrected if this is kept.

    """
    lalDict = lal.CreateDict()
    lalsimulation.SimInspiralWaveformParamsInsertPhenomXPrecVersion(lalDict, 223)

    expected_result = np.array([1468.09702243, 190.76614342, 972.51053189, 189.80404795])

    new_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXPHM, [[2,2]], lalparams=lalDict, generic_spins=True))

    np.testing.assert_allclose(new_result, expected_result, rtol=1e-6, err_msg="IMRPhenomXP_MB test failed")


def test_IMRPhenomXP_SpinTaylor():
    """
    This test checks that IMRPhenomXP with the SpinTaylor precession option (version 310) hasn't changed.
    It does this by generating two PhenomXP waveforms and computing
    their difference (according to their amplitude and phases)
    and compares them to pre-computed values.

    these pre-computed values were computed using the following code:

    `lalDict = lal.CreateDict()`

    `lalsimulation.SimInspiralWaveformParamsInsertPhenomXPrecVersion(lalDict, 310)`

    `expected_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXP, None, lalparams=lalDict, generic_spins=True))`

    """

    expected_result = np.array([1570.905974,  190.514064, 1107.696605,  195.697882])

    lalDict = lal.CreateDict()

    lalsimulation.SimInspiralWaveformParamsInsertPhenomXPrecVersion(lalDict, 310)

    new_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXP, None, lalparams=lalDict, generic_spins=True))

    np.testing.assert_allclose(new_result, expected_result, rtol=1e-6, err_msg="IMRPhenomXP test failed")


def test_IMRPhenomXPHM():
    """
    This test checks that IMRPhenomXPHM hasn't changed.
    It does this by generating two PhenomXPHM waveforms and computing
    their difference (according to their amplitude and phases)
    and compares them to pre-computed values.

    these pre-computed values were computed using the following line:

    `expected_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXPHM, mode_array, generic_spins=True))`

    where mode_array is [[2,2],[2,1],[3,3],[4,4]] or [[3,2]]

    """


    expected_result = np.array([1166.01091848, 334.5693217,  767.82099062, 326.09652364])

    new_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXPHM, [[2,2],[2,1],[3,3],[4,4]], release=122019))

    # rtol here needs to be more lenient to pass on builds with arm64 or MKL
    np.testing.assert_allclose(new_result, expected_result, rtol=1e-5, err_msg="IMRPhenomXPHM (122019) no 32 mode test failed")


    expected_result = np.array([68.9282789725476, 240.20999880535206, 25.111569754767224, 234.7465084316962])

    new_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXPHM, [[3,2]], release=122019))

    # rtol with 32 mode needs to b more lenient
    np.testing.assert_allclose(new_result, expected_result, rtol=3e-4, err_msg="IMRPhenomXPHM (122019) 32 mode test failed")

    expected_result = np.array([1166.77270896, 334.86014307, 768.93672645, 326.38518250])

    new_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXPHM, [[2,2],[2,1],[3,3],[4,4]]))

    # rtol here needs to be more lenient to pass on builds with arm64 or MKL
    np.testing.assert_allclose(new_result, expected_result, rtol=1.3e-5, err_msg="IMRPhenomXPHM no 32 mode test failed")

    expected_result = np.array([71.43504434, 242.82287296, 26.54528442, 237.35077401])

    new_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXPHM, [[3,2]]))

    # rtol with 32 mode needs to b more lenient
    np.testing.assert_allclose(new_result, expected_result, rtol=3e-4, err_msg="IMRPhenomXPHM 32 mode test failed")

def test_IMRPhenomXO4a():
    """
    This test checks that IMRPhenomXO4a hasn't changed.
    It does this by generating two IMRPhenomXO4a waveforms and computing

    `expected_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXO4a))`

    """

    expected_result = np.array([1335.4548650055917,  327.35598432690784,  958.1145964826108,  327.67901031324624])

    new_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXO4a, [[2,2],[2,1],[3,3],[4,4]]))

    # rtol here needs to be more lenient to pass on builds with arm64 or MKL
    np.testing.assert_allclose(new_result, expected_result, rtol=1e-5, err_msg="IMRPhenomXO4a no 32 mode test failed")


    expected_result = np.array([ 48.56333067501462, 236.2424809925183,  21.61677264763599, 232.45649243965286 ])

    new_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXO4a, [[3,2]]))

    # rtol with 32 mode needs to be more lenient
    np.testing.assert_allclose(new_result, expected_result, rtol=3e-3, err_msg="IMRPhenomXO4a 32 mode test failed")

def test_IMRPhenomXPHM_SpinTaylor():
    """
    This test checks that IMRPhenomXPHM with the SpinTaylor precession option (version 310) hasn't changed.
    It does this by generating two PhenomXPHM waveforms and computing
    their difference (according to their amplitude and phases)
    and compares them to pre-computed values.

    these pre-computed values were computed using the following code:

    `lalDict = lal.CreateDict()`

    `lalsimulation.SimInspiralWaveformParamsInsertPhenomXPrecVersion(lalDict, 310)`

    `expected_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXPHM, mode_array, lalparams=lalDict, generic_spins=True))`

    where mode_array is [[2,2],[2,1],[3,3],[4,4]] or [[3,2]]

    """

    expected_result = np.array([1688.30370786,  274.96949069, 1248.22149474,  279.10374629])

    lalDict = lal.CreateDict()

    lalsimulation.SimInspiralWaveformParamsInsertPhenomXPrecVersion(lalDict, 310)

    new_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXPHM, [[2,2],[2,1],[3,3],[4,4]], lalparams=lalDict, generic_spins=True,release=122019))

    # rtol here needs to be more lenient to pass on builds with arm64 or MKL
    np.testing.assert_allclose(new_result, expected_result, rtol=1e-5, err_msg="IMRPhenomXPHM no 32 mode test failed")


    expected_result = np.array([12.68276953, 170.44495875,  42.97960664, 142.53609984])

    new_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXPHM, [[3,2]], lalparams=lalDict, generic_spins=True,release=122019))

    # rtol with 32 mode needs to be more lenient
    np.testing.assert_allclose(new_result, expected_result, rtol=5e-4, err_msg="IMRPhenomXPHM 32 mode test failed")


def test_IMRPhenomXAS_NRTidalv2():
    """
    This test checks that IMRPhenomXAS_NRTidalv2 hasn't changed.
    It does this by generating two IMRPhenomXAS_NRTidalv2 waveforms and computing
    their difference (according to their amplitude and phases)
    and compares them to pre-computed values.

    these pre-computed values were computed using the following line:

    `expected_result  =  np.array(gen_test_data(0., lalsimulation.IMRPhenomXAS_NRTidalv2, None, type="bns", generic_spins=True))`

    """

    expected_result = np.array([10.04639297, 561.5921041 ,8.03711438, 561.62725166])

    new_result  =  np.array(gen_test_data(0., lalsimulation.IMRPhenomXAS_NRTidalv2, None, type='bns', generic_spins=True))

    np.testing.assert_allclose(new_result, expected_result, rtol=1e-6, err_msg="IMRPhenomXAS_NRTidalv2 test failed")


def test_IMRPhenomXP_NRTidalv2():
    """
    This test checks that IMRPhenomXP_NRTidalv2 hasn't changed.
    It does this by generating two IMRPhenomXP_NRTidalv2 waveforms and computing
    their difference (according to their amplitude and phases)
    and compares them to pre-computed values.

    these pre-computed values were computed using the following line:

    `expected_result  =  np.array(gen_test_data(0.2, lalsimulation.IMRPhenomXP_NRTidalv2, None, type="bns", generic_spins=True))`

    """

    lalDict = lal.CreateDict()

    expected_result = np.array([ 13.9202092 , 561.0238095 ,  19.05392711, 550.93840153])

    new_result  =  np.array(gen_test_data(0.2, lalsimulation.IMRPhenomXP_NRTidalv2, None, type='bns', lalparams=lalDict, generic_spins=True))
    np.testing.assert_allclose(new_result, expected_result, rtol=1e-6, err_msg="IMRPhenomXP_NRTidalv2 test failed")


def test_IMRPhenomXP_NRTidalv2_SpinTaylor():
    """
    This test checks that IMRPhenomXP_NRTidalv2 with the SpinTaylor precession (version 310) hasn't changed.
    It does this by generating two IMRPhenomXP_NRTidalv2 waveforms and computing
    their difference (according to their amplitude and phases)
    and compares them to pre-computed values.

    these pre-computed values were computed using the following code:

    `lalDict = lal.CreateDict()`

    `lalsimulation.SimInspiralWaveformParamsInsertPhenomXPrecVersion(lalDict, 310)`

    `expected_result  =  np.array(gen_test_data(0.2, lalsimulation.IMRPhenomXP_NRTidalv2, None, lalparams=lalDict, type="bns", generic_spins=True))`

    """

    expected_result = np.array([11.92543681, 730.77807026,  13.69906426, 552.53065783])

    lalDict = lal.CreateDict()

    lalsimulation.SimInspiralWaveformParamsInsertPhenomXPrecVersion(lalDict, 310)

    new_result  =  np.array(gen_test_data(0.2, lalsimulation.IMRPhenomXP_NRTidalv2, None, lalparams=lalDict, type='bns', generic_spins=True))

    np.testing.assert_allclose(new_result, expected_result, rtol=1e-6, err_msg="IMRPhenomXP_NRTidalv2 test failed")

# -- run the tests ------------------------------

if __name__ == '__main__':
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-phenomX.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
