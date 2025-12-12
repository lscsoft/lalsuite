from gwsignal.core import waveform as wfm
from gwsignal.core import eccentricity_utils as ecc
import astropy.units as u
from gw_eccentricity import measure_eccentricity
from pycbc.filter import match, optimized_match
from pycbc.types import TimeSeries as PyCBCTimeSeries

from copy import deepcopy
import pytest
import numpy as np

from .test_utilities import compute_match


@pytest.mark.parametrize("eccentricity", np.linspace(0, 0.9, num=10))
def test_circular_conversion(eccentricity):
    mean_anomaly = np.linspace(0, 2 * np.pi, num=100)

    recomputed_mean = ecc.mean_anomaly_from_true(
        ecc.true_anomaly_from_eccentric(
            ecc.eccentric_anomaly_from_mean(mean_anomaly, eccentricity), eccentricity
        ),
        eccentricity,
    )

    assert np.allclose(mean_anomaly, recomputed_mean)


@pytest.mark.xfail(
    reason="Measuring at the wrong frequency. WIP. "
    "For now, just a demo of gw_eccentricity usage"
)
def test_eccentricity_value(gen, parameters, plot):
    parameters["f22_start"] = 10.0 * u.Hz
    parameters["f22_ref"] = 10.0 * u.Hz

    modes = wfm.GenerateTDModes(parameters, gen)

    modes_dict_numpy = {ellm: mode.value for ellm, mode in modes.items()}
    times = modes[(2, 2)].times.value

    result_dict = measure_eccentricity(
        fref_in=10.0, dataDict={"t": times, "hlm": modes_dict_numpy}
    )
    # keys: fref_out, eccentricity, mean_anomaly, gwecc_object
    assert result_dict["eccentricity"] == 0.5

    if plot:
        import matplotlib.pyplot as plt

        gwecc_object = result_dict["gwecc_object"]
        # fig, ax = gwecc_object.make_diagnostic_plots()
        fig, ax = gwecc_object.plot_eccentricity()
        plt.show()
        plt.close()
        fig, ax = gwecc_object.plot_mean_anomaly()
        plt.show()
        plt.close()
        fig, ax = gwecc_object.plot_omega22()
        plt.show()
        plt.close()


@pytest.mark.xfail(
    reason="Waveforms generated change significantly in length for small eccentricity!"
)
def test_quasicircular_limit(gen, parameters, plot):
    eccentricities = np.geomspace(1e-6, 1e-3, num=20)

    zeroecc_params = deepcopy(parameters)
    zeroecc_params["eccentricity"] = 0.0 * u.dimensionless_unscaled

    zeroecc_hp, zeroecc_hc = wfm.GenerateTDWaveform(zeroecc_params, gen)

    mp_array = []
    mc_array = []

    for ecc in eccentricities:
        new_params = deepcopy(parameters)
        new_params["eccentricity"] = ecc * u.dimensionless_unscaled

        hp, hc = wfm.GenerateTDWaveform(new_params, gen)

        if plot:
            import matplotlib.pyplot as plt

            plt.plot(hp.times - hp.times[np.argmax(hp)], hp)
            plt.plot(
                zeroecc_hp.times - zeroecc_hp.times[np.argmax(zeroecc_hp)], zeroecc_hp
            )
            plt.show()

        mp, _ = compute_match(hp, zeroecc_hp)
        mc, _ = compute_match(hc, zeroecc_hc)
        mp_array.append(1 - mp)
        mc_array.append(1 - mc)

    assert sorted(mp_array) == mp_array
    assert sorted(mc_array) == mc_array

    if plot:
        import matplotlib.pyplot as plt

        plt.loglog(eccentricities, mp_array)
        plt.loglog(eccentricities, mc_array)
        plt.show()
