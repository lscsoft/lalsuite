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

"""Simple test to see if the PhenomX family models have changed: IMRPhenomXAS, IMRPhenomXHM, IMRPhenomXP and IMRPhenomXPHM.
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

def gen_test_data(spin1x, approximant, mode_array):
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

    common_pars=dict(
    m1=50*lal.MSUN_SI,
    m2=30*lal.MSUN_SI,
    S1x=spin1x,
    S1y=0.,
    S1z=0.,
    S2x=0.,
    S2y=0.,
    S2z=0.,
    distance=1,
    inclination=np.pi/3.,
    phiRef=0.,
    longAscNodes=0.,
    eccentricity=0.,
    meanPerAno=0.,
    deltaF=1./4.,
    f_min=30.,
    f_max=512.,
    f_ref=30.,
    LALpars=lalparams,
    approximant=approximant
    )

    pars1=common_pars.copy()

    pars2=common_pars.copy()
    pars2.update({"m2":20.*lal.MSUN_SI})

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

    `expected_result  =  np.array(gen_test_data(0., lalsimulation.IMRPhenomXAS))`

    """

    expected_result = np.array([781.6954825979799, 240.30628831736627, 625.356386078384, 240.30628831736624])

    new_result  =  np.array(gen_test_data(0., lalsimulation.IMRPhenomXAS, None))

    np.testing.assert_allclose(new_result, expected_result, rtol=1e-6, err_msg="IMRPhenomXAS test failed")

def test_IMRPhenomXHM():
    """
    This test checks that IMRPhenomXHM hasn't changed.
    It does this by generating two PhenomXHM waveforms and computing
    their difference (according to their amplitude and phases)
    and compares them to pre-computed values.

    these pre-computed values were computed using the following line:

    `expected_result  =  np.array(gen_test_data(0., lalsimulation.IMRPhenomXHM))`

    """

    expected_result = np.array([1005.1602387413318, 170.00035596090245, 768.1841403827161, 169.12938008572488])

    new_result  =  np.array(gen_test_data(0., lalsimulation.IMRPhenomXHM, [[2,2],[2,-2],[2,1],[2,-1],[3,3],[3,-3],[4,4],[4,-4]]))

    np.testing.assert_allclose(new_result, expected_result, rtol=1e-6, err_msg="IMRPhenomXHM no 32 mode test failed")

    expected_result = np.array([32.17820078828559, 215.95386544772447, 4.022275098535691, 215.9538654477245])

    new_result  =  np.array(gen_test_data(0., lalsimulation.IMRPhenomXHM, [[3,2],[3,-2]]))

    # rtol with 32 mode needs to b more lenient
    np.testing.assert_allclose(new_result, expected_result, rtol=3e-4, err_msg="IMRPhenomXHM 32 mode test failed")


def test_IMRPhenomXP():
    """
    This test checks that IMRPhenomXP hasn't changed.
    It does this by generating two PhenomXP waveforms and computing
    their difference (according to their amplitude and phases)
    and compares them to pre-computed values.

    these pre-computed values were computed using the following line:

    `expected_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXP))`

    """

    expected_result = np.array([1070.2208950681218, 271.6209567082584, 533.0487658836101, 268.45201460415575])

    new_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXP, None))

    np.testing.assert_allclose(new_result, expected_result, rtol=1e-6, err_msg="IMRPhenomXP test failed")


def test_IMRPhenomXPHM():
    """
    This test checks that IMRPhenomXPHM hasn't changed.
    It does this by generating two PhenomXPHM waveforms and computing
    their difference (according to their amplitude and phases)
    and compares them to pre-computed values.

    these pre-computed values were computed using the following line:

    `expected_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXPHM))`

    """

    expected_result = np.array([1166.0110734324703, 334.37564200163365, 767.8213461946891, 326.0800374360311])

    new_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXPHM, [[2,2],[2,1],[3,3],[4,4]]))

    # rtol here needs to be more lenient to pass on builds with arm64 or MKL
    np.testing.assert_allclose(new_result, expected_result, rtol=1e-5, err_msg="IMRPhenomXPHM no 32 mode test failed")


    expected_result = np.array([68.09797391055203, 235.35545843385682, 24.39359538033152, 229.92386957103975])

    new_result  =  np.array(gen_test_data(0.5, lalsimulation.IMRPhenomXPHM, [[3,2]]))

    # rtol with 32 mode needs to b more lenient
    np.testing.assert_allclose(new_result, expected_result, rtol=3e-4, err_msg="IMRPhenomXPHM 32 mode test failed")

# -- run the tests ------------------------------

if __name__ == '__main__':
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-phenomX.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
