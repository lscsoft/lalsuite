# -*- coding: utf-8 -*-
#
# Copyright (C) 2022 Colm Talbot
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
Test to see if the frequency series returned by lalsimulation.SimInspiralFD
has the expected behaviour.
"""

import sys

import lal
import lalsimulation
import numpy as np
import pytest

_test_data = [(0.25, 1024.0), (0.25, 600.0), (0.25, 40.0), (0.0, 1024.0), (0.0, 600.0)]


def _generate_time_and_frequency_waveforms(deltaF, f_max):
    parameters = dict(
        m1=40 * lal.MSUN_SI,
        m2=40 * lal.MSUN_SI,
        deltaF=deltaF,
        f_min=20.0,
        f_max=f_max,
        distance=1e8 * lal.PC_SI,
        LALparams=lal.CreateDict(),
    )
    for key in [
        "S1x", "S1y", "S1z", "S2x", "S2y", "S2z",
        "inclination", "phiRef", "longAscNodes",
        "eccentricity", "meanPerAno", "f_ref",
    ]:
        parameters[key] = 0.0
    td_approximant = lalsimulation.IMRPhenomT
    fd_approximant = lalsimulation.IMRPhenomXAS
    hplus_td = lalsimulation.SimInspiralFD(**parameters, approximant=td_approximant)[0]
    hplus_fd = lalsimulation.SimInspiralFD(**parameters, approximant=fd_approximant)[0]
    return hplus_fd, hplus_td


@pytest.mark.parametrize("deltaF, f_max", _test_data)
def test_deltaf_and_length_matches_between_td_and_fd_for_siminspiralfd(deltaF, f_max):
    hplus_fd, hplus_td = _generate_time_and_frequency_waveforms(deltaF=deltaF, f_max=f_max)
    assert hplus_td.deltaF == hplus_fd.deltaF, "frequency spacing mismatch"
    assert hplus_td.data.length == hplus_fd.data.length, "output length mismatch"


@pytest.mark.parametrize("deltaF, f_max", _test_data)
def test_overlap_between_time_and_frequency_domain_approximants(deltaF, f_max):
    """
    Test that the overlap maximized over time and phase is close to 1 for
    comparable waveform models.
    """
    if f_max < lalsimulation.IMRPhenomDGetPeakFreq(40.0, 40.0, 0.0, 0.0):
        pytest.skip(
            "Skipping overlap test as frequency band doesn't contain full signal."
        )
    hplus_fd, hplus_td = _generate_time_and_frequency_waveforms(deltaF=deltaF, f_max=f_max)
    frequency_array = np.arange(len(hplus_fd.data.data)) * hplus_fd.deltaF
    mask = (frequency_array >= 20) & (frequency_array <= f_max)
    hplus_fd = hplus_fd.data.data * mask
    hplus_td = hplus_td.data.data * mask
    psd = np.ones(len(frequency_array))
    psd[mask] = np.array([
        lalsimulation.SimNoisePSDaLIGOZeroDetHighPower(frequency)
        for frequency in frequency_array[mask]]
    )
    max_inner_product = max(abs(np.fft.fft(hplus_fd * hplus_td.conjugate() / psd)))
    overlap = (
        max_inner_product / (
            sum(abs(hplus_fd)**2 / psd) ** 0.5
            * sum(abs(hplus_td)**2 / psd) ** 0.5
        )
    )
    assert abs(overlap - 1) < 0.01


if __name__ == '__main__':
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-siminspiralfd.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
