# -*- coding: utf-8 -*-
#
# Copyright (C) 2020 Michalis Agathos
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

# define required parameters for time domain and frequency domain
# as well as their default values

"""Simple test to see if TEOBResumS has changed

Some functions copied from test_phenomPv3HM.py.
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

def gen_test_data():
    """
    compute the difference between two waveforms
    and compare to expected value
    """

    common_pars=dict(
    m1=50*lal.MSUN_SI,
    m2=30*lal.MSUN_SI,
    s1x=0,
    s1y=0,
    s1z=-0.45,
    s2x=0,
    s2y=0,
    s2z=0.98,
    distance=1,
    inclination=0.,
    phiRef=0.,
    longAscNodes=0.,
    eccentricity=0.,
    meanPerAno=0.,
    deltaT=1./4096.,
    f_min=30.,
    f_ref=30.,
    params=None,
    approximant=lalsimulation.TEOBResumS
    )

    pars1 = common_pars.copy()
    pars2 = common_pars.copy()
    pars2.update({"m2":20.*lal.MSUN_SI})

    hp1, hc1 = lalsimulation.SimInspiralChooseTDWaveform(**pars1)
    hp2, hc2 = lalsimulation.SimInspiralChooseTDWaveform(**pars2)

    # compute amp and phase
    hp1_amp, hp1_phase = get_amp_phase(hp1.data.data)
    hc1_amp, hc1_phase = get_amp_phase(hc1.data.data)

    hp2_amp, hp2_phase = get_amp_phase(hp2.data.data)
    hc2_amp, hc2_phase = get_amp_phase(hc2.data.data)

    npeak1 = np.argmax(hp1_amp)
    npeak2 = np.argmax(hp2_amp)
    ni = min(npeak1, npeak2)
    no = min(hp1.data.length - npeak1, hp2.data.length - npeak2)

    # Truncate and align waveforms in place w.r.t. peak
    hp1_amp = hp1_amp[npeak1-ni:npeak1+no]
    hc1_amp = hc1_amp[npeak1-ni:npeak1+no]
    hp2_amp = hp2_amp[npeak2-ni:npeak2+no]
    hc2_amp = hc2_amp[npeak2-ni:npeak2+no]
    hp1_phase = hp1_phase[npeak1-ni:npeak1+no]
    hc1_phase = hc1_phase[npeak1-ni:npeak1+no]
    hp2_phase = hp2_phase[npeak2-ni:npeak2+no]
    hc2_phase = hc2_phase[npeak2-ni:npeak2+no]

    hp_amp_diff   = sum_sqr_diff(hp1_amp, hp2_amp)
    hp_phase_diff = sum_sqr_diff(hp1_phase, hp2_phase)

    hc_amp_diff   = sum_sqr_diff(hc1_amp, hc2_amp)
    hc_phase_diff = sum_sqr_diff(hc1_phase, hc2_phase)

    # since we want to compare decimals, we return the above quantities as 0.xxxxx..
    # by dividing them for their order of magnitude +1
    return hp_amp_diff/1e6, hp_phase_diff/1e3, hc_amp_diff/1e6, hc_phase_diff/1e3



# -- test functions ---------------------

def test_TEOBResumS():
    """
    This test checks that TEOBResumS hasn't changed.
    It does this by generating two TEOBResumS waveforms and computing
    their difference (according to their amplitude and phases)
    and compares them to pre-computed values.

    these pre-computed values were computed using the following line:

    `expected_result  =  np.array(gen_test_data())`

    """


    expected_result = np.array([0.27385, 0.14397, 0.27356, 0.14324])
    new_result  =  np.array(gen_test_data())
    np.testing.assert_almost_equal(new_result, expected_result, 5, "TEOBResumS test failed")



# -- run the tests ------------------------------

if __name__ == '__main__':
    sys.exit(pytest.main(args=[__file__] + sys.argv[1:] + ['-v']))
