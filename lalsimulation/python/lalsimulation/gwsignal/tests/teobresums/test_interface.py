import pytest

from gwsignal.core import waveform as wfm
from gwsignal.models import teobresums

import EOBRun_module

import astropy.units as u
import numpy as np
from gwpy.timeseries import TimeSeries



def test_modes_dictionary():
    assert len(teobresums.TEOB_DALI_MODES) == len(teobresums.TEOB_DALI_MODES_FROM_K)


def test_modes_dictionary_backreference():
    for mode in teobresums.TEOB_DALI_MODES:
        assert (
            teobresums.TEOB_DALI_MODES_FROM_K[teobresums.modes_to_k([mode])[0]] == mode
        )


def test_modes_gen(gen, parameters):
    hlm = wfm.GenerateTDModes(parameters, gen)

    assert isinstance(hlm[(2, 2)], TimeSeries)


@pytest.fixture
def hphc_native_vs_gwsignal(gen, parameters):
    if gen.__class__.__name__ != "TEOBResumSDALI":
        pytest.skip()

    hp_gwsignal, hc_gwsignal = wfm.GenerateTDWaveform(parameters, gen)

    unitless_parameters = gen._strip_units(parameters)
    teob_parameters = teobresums.convert_parameters_to_teob(unitless_parameters)

    teob_parameters["use_mode_lm"] = gen._available_modes_teob_convention

    t, hp, hc, htlm, dyn = EOBRun_module.EOBRunPy(teob_parameters)

    return t, hp, hp_gwsignal, hc, hc_gwsignal


def test_hphc_gen(hphc_native_vs_gwsignal):
    t, hp, hp_gwsignal, hc, hc_gwsignal = hphc_native_vs_gwsignal

    assert np.allclose(hp, hp_gwsignal, atol=1e-30)
    assert np.allclose(hc, hc_gwsignal, atol=1e-30)
    assert np.allclose(t, hp_gwsignal.times.value)


def test_plotting_hp(plot, hphc_native_vs_gwsignal):
    t, hp, hp_gwsignal, hc, hc_gwsignal = hphc_native_vs_gwsignal

    if plot:
        import matplotlib.pyplot as plt

        plt.plot(
            hp_gwsignal.times, hp_gwsignal.real, label="$h_+$ reconstructed by gwsignal"
        )
        plt.plot(t, hp, label="$h_+$ from TEOBResumS")
        plt.xlabel("$t$ [s]")
        plt.legend()
        plt.show()


def test_plotting_hc(plot, hphc_native_vs_gwsignal):
    t, hp, hp_gwsignal, hc, hc_gwsignal = hphc_native_vs_gwsignal

    if plot:
        import matplotlib.pyplot as plt

        plt.plot(
            hc_gwsignal.times,
            hc_gwsignal.real,
            label="$h_\\times$ reconstructed by gwsignal",
        )
        plt.plot(t, hc, label="$h_\\times$ from model")
        plt.xlabel("$t$ [s]")
        plt.legend()
        plt.show()


def test_plotting_modes(plot, gen, parameters):
    hlm = wfm.GenerateTDModes(parameters, gen)

    if plot:
        import matplotlib.pyplot as plt

        plt.plot(hlm[(2, 2)].times, hlm[(2, 2)].real, label="$\\Re(h_{22})$")
        plt.plot(hlm[(3, 3)].times, hlm[(3, 3)].real, label="$\\Re(h_{32})$")
        plt.legend()
        plt.xlabel("$t$ [s]")
        plt.show()
