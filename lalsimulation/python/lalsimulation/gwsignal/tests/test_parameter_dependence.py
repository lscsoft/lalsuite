import pytest
from hypothesis import given, settings, HealthCheck
from hypothesis import strategies as st

from gwsignal.core import waveform as wfm
from gwsignal.models import teobresums

import astropy.units as u
import numpy as np
from gwpy.timeseries.timeseries import TimeSeries

from scipy.interpolate import interp1d

from pycbc.filter import sigma
from pycbc.types import TimeSeries as PyCBCTimeSeries

from copy import deepcopy

from .test_utilities import compute_match

# number of samples for the grids in inclination, distance etc
N_SAMPLES = 15


def assert_indistinguishable(wf1, wf2):
    if isinstance(wf1, TimeSeries):
        wf1 = wf1.value
    if isinstance(wf2, TimeSeries):
        wf2 = wf2.value

    assert np.allclose(wf1, wf2, atol=0, rtol=1e-8)


def swap_masses(parameter_dictionary):
    new_parameters = deepcopy(parameter_dictionary)

    to_swap = ["mass#", "spin#z"]

    for param_name in to_swap:
        name_1 = param_name.replace("#", "1")
        name_2 = param_name.replace("#", "2")
        new_parameters[name_1], new_parameters[name_2] = (
            new_parameters[name_2],
            new_parameters[name_1],
        )

    return new_parameters


def scale_initial_conditions_by_factor(parameter_dictionary, factor):
    new_parameters = deepcopy(parameter_dictionary)
    new_parameters["mass1"] *= factor
    new_parameters["mass2"] *= factor
    new_parameters["f22_start"] /= factor
    new_parameters["f22_ref"] /= factor
    return new_parameters


def test_mass_inversion_symmetry(parameters, gen):
    hp, hc = wfm.GenerateTDWaveform(parameters, gen)

    hp2, hc2 = wfm.GenerateTDWaveform(swap_masses(parameters), gen)

    assert_indistinguishable(hp, hp2)
    assert_indistinguishable(hc, hc2)


@pytest.mark.xfail(reason="Maybe it is a resampling issue?")
@pytest.mark.parametrize("scale_factor", [2.0, 10.0])
def test_mass_scaling(parameters, scale_factor, plot, gen):
    hp, hc = wfm.GenerateTDWaveform(parameters, gen)

    hp2, hc2 = wfm.GenerateTDWaveform(
        scale_initial_conditions_by_factor(parameters, scale_factor), gen
    )

    # the rescaled waveform will last more or less, but they'll have the same sample rate
    # so, we need to resample them
    resampled_hp2 = interp1d(hp2.times / scale_factor, hp2 / scale_factor)(hp.times)
    resampled_hc2 = interp1d(hc2.times / scale_factor, hc2 / scale_factor)(hp.times)

    h_complex = hp + 1j * hc
    h2_complex = resampled_hp2 + 1j * resampled_hc2

    if plot:
        import matplotlib.pyplot as plt

        plt.plot(hp.times, hp)
        plt.plot(hp.times, resampled_hp2)
        plt.ylabel("$h_+$ [dimensionless]")
        ax2 = plt.twinx()
        ax2.plot(hp.times, abs(h_complex - h2_complex) / abs(h_complex))
        ax2.set_ylabel("$|h_1 - h_2| / |h_1|$")
        plt.show()

    assert_indistinguishable(h_complex, h2_complex)


@pytest.mark.parametrize("scale_factor", [2.0, 4.0, 8.0, 16.0, 32.0])
def test_sample_rate_scaling(parameters, scale_factor, plot, gen):
    # these values reproduce the ones used for the v5 review,
    # https://git.ligo.org/waveforms/reviews/seobnrv5/-/blob/main/aligned/sanity_checks/standard_tests/plots/varying_srate.pdf

    hp, hc = wfm.GenerateTDWaveform(parameters, gen)

    new_parameters = deepcopy(parameters)

    new_parameters["deltaT"] *= scale_factor

    hp2, hc2 = wfm.GenerateTDWaveform(new_parameters, gen)

    # resample on the sparser of the two grids
    times = hp2.times
    new_hp = interp1d(hp.times, hp, bounds_error=False, fill_value=0.0)(times)
    new_hp2 = hp2

    if plot and scale_factor == 2.0:
        import matplotlib.pyplot as plt

        plt.plot(
            times,
            new_hp,
            label=f"Baseline sample rate ({(1 / parameters['deltaT']).to(u.Hz):.0f}), downsampled",
        )
        plt.plot(
            times,
            new_hp2,
            label=f"New sample rate ({(1 / new_parameters['deltaT']).to(u.Hz):.0f})",
        )

        plt.xlabel("Time [s]")
        plt.ylabel("$h_+$ amplitude [dimensionless]")
        plt.legend()
        plt.show()

    assert_indistinguishable(new_hp, new_hp2)


def test_distance_scaling(parameters, plot, gen):
    distances = np.geomspace(10, 10000, num=N_SAMPLES) * u.Mpc

    hp, hc = wfm.GenerateTDWaveform(parameters, gen)
    baseline_distance = parameters["distance"]
    baseline_maximum = abs(hp).max()

    maxima = []
    for distance in distances:
        new_parameters = deepcopy(parameters)
        new_parameters["distance"] = distance

        scale_factor = (distance / baseline_distance).value

        hp2, hc2 = wfm.GenerateTDWaveform(new_parameters, gen)

        assert_indistinguishable(hp, hp2 * scale_factor)

        # NOTE: calling max(abs(hp2)) here gives the wrong result,
        # not the maximum of the array
        maxima.append(abs(hp2).max())

    if plot:
        import matplotlib.pyplot as plt

        plt.scatter(distances, maxima, label="Computed maxima")
        plt.plot(
            distances,
            baseline_maximum / (distances / baseline_distance),
            label="Theoretical maxima",
            c="black",
        )
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("Distance [Mpc]")
        plt.ylabel("Maximum of $|h_+|$ [dimensionless]")
        plt.title("Distance scaling of the waveform")
        plt.legend()
        plt.show()


def test_inclination_dependence(parameters, plot, gen):
    inclinations = np.linspace(0, np.pi, num=N_SAMPLES) * u.rad

    sigmas = []
    for inclination in inclinations:
        new_parameters = deepcopy(parameters)
        new_parameters["inclination"] = inclination

        hp, hc = wfm.GenerateTDWaveform(new_parameters, gen)
        # breakpoint()
        hp_pycbc = PyCBCTimeSeries(hp.data, delta_t=hp.dt.value)

        sigmas.append(
            sigma(hp_pycbc, low_frequency_cutoff=new_parameters["f22_start"].value)
        )

    if plot:
        import matplotlib.pyplot as plt

        plt.plot(inclinations, sigmas)
        plt.xlabel("Inclination [rad]")
        plt.ylabel("Waveform norm [dimensionless]")
        plt.show()

    # norms should be decreasing from 0 to pi/2
    assert sigmas[: N_SAMPLES // 2] == sorted(sigmas[: N_SAMPLES // 2], reverse=True)
    # and increasing from pi/2 to pi
    assert sigmas[N_SAMPLES // 2 :] == sorted(sigmas[N_SAMPLES // 2 :])

    assert np.isclose(max(sigmas), sigmas[0])


@pytest.mark.xfail(
    reason="Times are strangely defined - is the epoch supposed not to be zero?"
)
def test_delta_t_close_to_zero(gen, parameters):
    hp, hc = wfm.GenerateTDWaveform(parameters, gen)

    max_index = np.argmax(hp**2 + hc**2)
    delta_t_max = hp.times[max_index]
    assert np.isclose(delta_t_max, 0, atol=1e-3)


def test_orbital_phase(gen, parameters, plot):
    parameters["eccentricity"] = 0.01 * u.dimensionless_unscaled

    parameters["phi_ref"] = 0.0 * u.rad
    hp, hc = wfm.GenerateTDWaveform(parameters, gen)

    parameters_1 = deepcopy(parameters)
    parameters_1["phi_ref"] = 1 * np.pi / 4 * u.rad
    hp1, hc1 = wfm.GenerateTDWaveform(parameters_1, gen)

    parameters_2 = deepcopy(parameters)
    parameters_2["phi_ref"] = 2 * np.pi / 4 * u.rad
    hp2, hc2 = wfm.GenerateTDWaveform(parameters_2, gen)

    parameters_3 = deepcopy(parameters)
    parameters_3["phi_ref"] = 3 * np.pi / 4 * u.rad
    hp3, hc3 = wfm.GenerateTDWaveform(parameters_3, gen)

    parameters_4 = deepcopy(parameters)
    parameters_4["phi_ref"] = 4 * np.pi / 4 * u.rad
    hp4, hc4 = wfm.GenerateTDWaveform(parameters_4, gen)

    assert np.isclose(compute_match(hp, hp1), 1, atol=5e-2)
    assert np.isclose(compute_match(hp, hp2), 1, atol=5e-2)
    assert np.isclose(compute_match(hp, hp3), 1, atol=5e-2)
    assert np.isclose(compute_match(hp, hp4), 1, atol=5e-2)

    if plot:
        import matplotlib.pyplot as plt

        plt.plot(hp.times, hp)
        plt.plot(hp1.times, hp1, label="Shift by $1\pi/4$")
        plt.plot(hp2.times, hp2, label="Shift by $2\pi/4$")
        plt.plot(hp3.times, hp3, label="Shift by $3\pi/4$")
        plt.plot(hp4.times, hp4, label="Shift by $4\pi/4$")
        plt.legend()
        plt.show()


def test_modearray_parameter(parameters, gen):
    from gwsignal.core.gw import SpinWeightedSphericalHarmonicMode

    # when giving a ModeArray, only those modes should be returned
    parameters["ModeArray"] = [[2, 2], [3, 3]]

    modes = wfm.GenerateTDModes(parameters, gen)

    assert SpinWeightedSphericalHarmonicMode(-2, 2, 2) in modes
    assert SpinWeightedSphericalHarmonicMode(-2, 2, -2) in modes
    assert SpinWeightedSphericalHarmonicMode(-2, 3, 3) in modes
    assert SpinWeightedSphericalHarmonicMode(-2, 3, -3) in modes
    assert SpinWeightedSphericalHarmonicMode(-2, 2, 1) not in modes

    # the number of modes is doubled, since we have both positive and negative m
    assert len(modes) == 4

    # if the user attempts to use a mode that is not available, an error is raised
    with pytest.raises(ValueError):
        wfm.GenerateTDModes(parameters | {"ModeArray": [[2, 3]]}, gen)


def test_modearray_parameter_polarizations(parameters, gen):
    hp1, hc1 = wfm.GenerateTDWaveform(parameters, gen)
    hp2, hc2 = wfm.GenerateTDWaveform(parameters | {"ModeArray": [[2, 2]]}, gen)

    assert not np.allclose(hp1.value, hp2.value, atol=max(hp1.value) / 1e3)
    assert not np.allclose(hc1.value, hc2.value, atol=max(hc1.value) / 1e3)
