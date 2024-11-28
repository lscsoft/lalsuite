# -*- coding: utf-8 -*-
#
# Copyright (C) 2024 Adrian Abac
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

"""
Simple test to see if IMRPhenomXP_NRTidalv3 have changed
Adapted from test_SEOBNRv5HM_ROM.py.
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

def gen_test_data(approximant):
    """
    compute the difference between two waveforms
    and compare to expected value
    """

    LALparams = lal.CreateDict()
    lambda1 = 400.0
    lambda2 = 600.0
    lalsimulation.SimInspiralWaveformParamsInsertTidalLambda1(LALparams, lambda1)
    lalsimulation.SimInspiralWaveformParamsInsertTidalLambda2(LALparams, lambda2)
    common_pars=dict(
    m1=1.4*lal.MSUN_SI,
    m2=1.2*lal.MSUN_SI,
    S1x=0,
    S1y=0,
    S1z=-0.45,
    S2x=0,
    S2y=0,
    S2z=0.98,
    distance=1,
    inclination=0.,
    phiRef=0.,
    longAscNodes=0.,
    eccentricity=0.,
    meanPerAno=0.,
    deltaF=1./4.,
    f_min=30.,
    f_max=512.,
    f_ref=30.,
    LALpars=LALparams,
    approximant=approximant
    )

    pars1 = common_pars.copy()

    pars2 = common_pars.copy()
    pars2.update({"m2":1.0*lal.MSUN_SI})
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


def test_IMRPhenomXP_NRTidalv3():
    """
    This test checks that IMRPhenomXP_NRTidalv3 hasn't changed.
    It does this by generating two IMRPhenomXP_NRTidalv3 waveforms and computing
    their difference (according to their amplitude and phases)
    and compares them to pre-computed values.

    these pre-computed values were computed using the following line:

    `expected_result  =  np.array(gen_test_data(lalsimulation.IMRPhenomXP_NRTidalv3))`
    """

    expected_result = np.array([34.92128149, 1092.39559414,   34.92128149, 1092.55480085])
    new_result  =  np.array(gen_test_data(lalsimulation.IMRPhenomXP_NRTidalv3))
    np.testing.assert_allclose(new_result, expected_result, rtol=0.002, err_msg="IMRPhenomXP_NRTidalv3 test failed", verbose=True)

#-- run the tests ------------------------------

if __name__ == '__main__':
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-IMRPhenomXP_NRTidalv3.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
