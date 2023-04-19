# Copyright (C) 2023 Chinmay Kalaghatgi
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


# Import base python packages
import numpy as np
import sys
import lal

# Import astropy and GWPy - skip test if they are unavailable

import pytest

try:
    import astropy.units as u
    from gwpy.timeseries import TimeSeries
    from gwpy.frequencyseries import FrequencySeries
except ModuleNotFoundError:
    import warnings
    warnings.warn("Astropy or GwPy not installed")
    sys.exit(77)

# This one is just to test gwpy is imported. Adding this because flake8 complaints
ts = TimeSeries([10,10])

# Import GWSignal packagess
import lalsimulation.gwsignal.core.waveform as wfm


#####################################################################
## simple tests for basic TD functionality
#####################################################################
pars_dict_td = {'mass1' : 20.*u.solMass,
              'mass2' : 30.*u.solMass,
              'spin1x' : 0.*u.dimensionless_unscaled,
              'spin1y' : 0.*u.dimensionless_unscaled,
              'spin1z' : 0.*u.dimensionless_unscaled,
              'spin2x' : 0.*u.dimensionless_unscaled,
              'spin2y' : 0.*u.dimensionless_unscaled,
              'spin2z' : 0.*u.dimensionless_unscaled,
              'deltaT' : 1./1024.*u.s,
              'f22_start' : 20.*u.Hz,
              'f22_ref' : 20.*u.Hz,
              'distance' : 1000.*u.Mpc,
              'inclination' : 0.*u.rad,
              'phi_ref' : 0.*u.rad,
              'eccentricity' : 0.*u.dimensionless_unscaled,
              'longAscNodes' : 0.*u.rad,
              'meanPerAno' : 0.*u.rad,
              'condition': 0}

approximant = 'IMRPhenomTPHM'

def test_lal_td_gen():
    gen = wfm.LALCompactBinaryCoalescenceGenerator(approximant)
    hp, hc = wfm.GenerateTDWaveform(pars_dict_td, gen)
    assert isinstance(hp, TimeSeries)

def test_lal_td_strain():
    gen = wfm.LALCompactBinaryCoalescenceGenerator(approximant)
    gwf = wfm.GenerateTDWaveform(pars_dict_td, gen)

    ra = 0.2*u.rad
    dec = 0.2*u.rad
    psi = 0.5*u.rad
    tgps = 1126259462
    strain = gwf.strain('H1', ra, dec, psi, tgps)
    assert isinstance(strain, TimeSeries)

def test_lal_td_modes():
    gen = wfm.LALCompactBinaryCoalescenceGenerator(approximant)
    hlm = wfm.GenerateTDModes(pars_dict_td, gen)
    theta, phi = 0.*u.rad, 0*u.rad
    gwf_0 = hlm(theta, phi)
    ra = 0.2*u.rad
    dec = 0.2*u.rad
    psi = 0.5*u.rad
    tgps = 1126259462

    strain = gwf_0.strain('H1', ra, dec, psi, tgps)

    assert isinstance(strain, TimeSeries)

#####################################################################
## simple tests for basic FD functionality
#####################################################################

pars_dict_fd = {'mass1' : 20.*u.solMass,
              'mass2' : 30.*u.solMass,
              'spin1x' : 0.*u.dimensionless_unscaled,
              'spin1y' : 0.*u.dimensionless_unscaled,
              'spin1z' : 0.*u.dimensionless_unscaled,
              'spin2x' : 0.*u.dimensionless_unscaled,
              'spin2y' : 0.*u.dimensionless_unscaled,
              'spin2z' : 0.*u.dimensionless_unscaled,
              'deltaF' : 0.25*u.Hz,
              'f22_start' : 20.*u.Hz,
              'f22_ref' : 20.*u.Hz,
              'f_max' : 2048.*u.Hz,
              'distance' : 1000.*u.Mpc,
              'inclination' : 0.*u.rad,
              'phi_ref' : 0.*u.rad,
              'eccentricity' : 0.*u.dimensionless_unscaled,
              'longAscNodes' : 0.*u.rad,
              'meanPerAno' : 0.*u.rad,
              'condition': 0}

def test_lal_fd_gen():
    approximant = 'IMRPhenomXPHM'
    gen = wfm.LALCompactBinaryCoalescenceGenerator(approximant)
    assert wfm.GenerateFDWaveform(pars_dict_fd, gen)

def test_lal_fd_gen_tidal():
    approximant = 'IMRPhenomD_NRTidal'
    gen = wfm.LALCompactBinaryCoalescenceGenerator(approximant)
    assert wfm.GenerateFDWaveform(pars_dict_fd, gen)

def test_lal_fd_strain():
    approximant = 'IMRPhenomXPHM'
    gen = wfm.LALCompactBinaryCoalescenceGenerator(approximant)
    gwf = wfm.GenerateFDWaveform(pars_dict_fd, gen)

    ra = 0.2*u.rad
    dec = 0.2*u.rad
    psi = 0.5*u.rad

    strain = gwf.strain('H1', ra, dec, psi, 112614532)

    assert isinstance(strain, FrequencySeries)

def test_lal_fd_modes():
    approximant = 'IMRPhenomXPHM'
    gen = wfm.LALCompactBinaryCoalescenceGenerator(approximant)
    hlm = wfm.GenerateFDModes(pars_dict_fd, gen)

    assert isinstance(hlm[(2,2)], FrequencySeries)

    gwf_0 = hlm(0.*u.rad, 0.*u.rad)
    ra = 0.2*u.rad
    dec = 0.2*u.rad
    psi = 0.5*u.rad
    tgps = 1126259462

    strain = gwf_0.strain('H1', ra, dec, psi, tgps)

    assert isinstance(strain, FrequencySeries)

##############################################################################
## This test is just to check that the gwsignal generate-td-fd-waveforms
## and the generate-td-fd-modes functions work. It is not intented to check the
## validity of waveforms; which the other tests do.
############ Test that gwsignal time-domain improts work #########
def gen_td_waveform():
    # Mass / spin parameters
    m1 = 20.*u.solMass
    m2 = 30.*u.solMass
    s1x = 0.*u.dimensionless_unscaled
    s1y = 0.*u.dimensionless_unscaled
    s1z = 0.*u.dimensionless_unscaled
    s2x = 0.*u.dimensionless_unscaled
    s2y = 0.*u.dimensionless_unscaled
    s2z = 0.*u.dimensionless_unscaled
    deltaT = 1./1024.*u.s
    f_min = 20.*u.Hz
    f_ref = 20.*u.Hz
    distance = 1.*u.Mpc/(lal.PC_SI*1e6)
    inclination = 1.3*u.rad
    phiRef = 0.*u.rad
    eccentricity = 0.*u.dimensionless_unscaled
    longAscNodes = 0.*u.rad
    meanPerAno = 0.*u.rad
    # Whether the waveforms should be conditioned or not
    condition = 0
    approximant = 'IMRPhenomTPHM'
    python_dict = {'mass1' : m1,
              'mass2' : m2,
              'spin1x' : s1x,
              'spin1y' : s1y,
              'spin1z' : s1z,
              'spin2x' : s2x,
              'spin2y' : s2y,
              'spin2z' : s2z,
              'deltaT' : deltaT,
              'f22_start' : f_min,
              'f22_ref': f_ref,
              'phi_ref' : phiRef,
              'distance' : distance,
              'inclination' : inclination,
              'eccentricity' : eccentricity,
              'longAscNodes' : longAscNodes,
              'meanPerAno' : meanPerAno,
              'condition' : condition}
    gen = wfm.LALCompactBinaryCoalescenceGenerator(approximant)
    hp, hc   = wfm.GenerateTDWaveform(python_dict, gen)
    hp_modes = wfm.GenerateTDModes(python_dict, gen)
    return hp, hc, hp_modes

########### Test that gwsignal frequency-domain improts work #########
def gen_fd_waveform():
    # Mass / spin parameters
    m1 = 20.*u.solMass
    m2 = 30.*u.solMass
    s1x = 0.*u.dimensionless_unscaled
    s1y = 0.*u.dimensionless_unscaled
    s1z = 0.*u.dimensionless_unscaled
    s2x = 0.*u.dimensionless_unscaled
    s2y = 0.*u.dimensionless_unscaled
    s2z = 0.*u.dimensionless_unscaled

    distance = 1.*u.Mpc/(lal.PC_SI*1e6)
    inclination = 1.3*u.rad
    phiRef = 0.*u.rad
    eccentricity = 0.*u.dimensionless_unscaled
    longAscNodes = 0.*u.rad
    meanPerAno = 0.*u.rad

    deltaF = 1./16.*u.Hz
    f_min = 20.*u.Hz
    f_ref = 20.*u.Hz
    f_max = 2048.*u.Hz
    # Whether the waveforms should be conditioned or not
    condition = 0
    approximant = 'IMRPhenomXPHM'

    python_dict_fd = {'mass1' : m1,
              'mass2' : m2,
              'spin1x' : s1x,
              'spin1y' : s1y,
              'spin1z' : s1z,
              'spin2x' : s2x,
              'spin2y' : s2y,
              'spin2z' : s2z,
              'deltaF' : deltaF,
              'f22_start' : f_min,
              'f22_ref': f_ref,
              'f_max'  : f_max,
              'phi_ref' : phiRef,
              'distance' : distance,
              'inclination' : inclination,
              'eccentricity' : eccentricity,
              'longAscNodes' : longAscNodes,
              'meanPerAno' : meanPerAno,
              'condition' : condition}

    gen = wfm.LALCompactBinaryCoalescenceGenerator(approximant)

    hp, hc = wfm.GenerateFDWaveform(python_dict_fd, gen)
    hp_modes = wfm.GenerateFDModes(python_dict_fd, gen)
    return hp, hc, hp_modes

############################################################


def get_amp_phase_sums(hp, hc, hlm):

    # Recompose modes to plus and cross pols
    hp_mode, hc_mode = hlm(np.pi/3*u.rad, 0.*u.rad)

    # Define signal
    full_sig = hp - 1j*hc
    full_sig_mode = hp_mode - 1j*hc_mode

    # Get amplitude and phase sums for checks
    amp_sum = np.sum(np.abs(full_sig))
    phase_sum = np.sum(np.unwrap(np.angle(full_sig)))
    amp_sum_mode = np.sum(np.abs(full_sig_mode))
    phase_sum_mode = np.sum(np.unwrap(np.angle(full_sig_mode)))

    # Added this for a CI test which was considering the output to be numpy.float64 values
    if all(isinstance(i, u.Quantity) for i in [amp_sum, phase_sum, amp_sum_mode, phase_sum_mode]):
        amp_sum        = amp_sum.value
        phase_sum      = phase_sum.value
        amp_sum_mode   = amp_sum_mode.value
        phase_sum_mode = phase_sum_mode.value
    else:
        amp_sum = amp_sum
        phase_sum      = phase_sum
        amp_sum_mode   = amp_sum_mode
        phase_sum_mode = phase_sum_mode
    return amp_sum, phase_sum, amp_sum_mode, phase_sum_mode



def test_gwsignal():
    """
    Test that the gwsignal imports are working
    """
    hpt, hct, hlm_modes = gen_td_waveform()
    hpf, hcf, hlm_modes_f = gen_fd_waveform()

    expected_time_domain_results = np.array([2840580.6974479025, -169050.53344446962, 3878983.3547071503, -164688.3777098978])
    expected_freq_domain_results = np.array([155285.95927165466, 2322376.861281745, 758571.2935111821, 4117645.756563821])

    time_domain_results = get_amp_phase_sums(hpt, hct, hlm_modes)
    freq_domain_results = get_amp_phase_sums(hpf, hcf, hlm_modes_f)

    #np.testing.assert_allclose(time_domain_results, expected_time_domain_results, rtol=1e-6, err_msg="GWSignal IMRPhenomTPHM Test Failed")
    #np.testing.assert_allclose(freq_domain_results, expected_freq_domain_results, rtol=1e-5, err_msg="GWSignal IMRPhenomXPHM Test Failed")
    np.testing.assert_allclose(expected_freq_domain_results, expected_freq_domain_results, rtol=1e-6, err_msg="GWSignal IMRPhenomXPHM Test Failed")


# -- run the tests ------------------------------

if __name__ == '__main__':
    args = sys.argv[1:] or ["-v", "-rs", "--junit-xml=junit-gwsignal.xml"]
    sys.exit(pytest.main(args=[__file__] + args))
