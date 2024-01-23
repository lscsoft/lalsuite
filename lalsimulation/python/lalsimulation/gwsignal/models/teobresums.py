"""Notes on TEOBResumS installation (TEMPORARY):

- clone it from https://bitbucket.org/eob_ihes/teobresums
- checkout the branch `dev/DALI-4PNrho22`
- within the `Python` folder, run `make`
- export the `PYTHONPATH` environment variable to the `Python` folder
    inside the TEOBResumS repo
"""

from __future__ import annotations
import warnings
from typing import Optional

try:
    import EOBRun_module
except ImportError as e:
    warnings.warn("The TEOBResumS package has failed to load")

from importlib_metadata import metadata
import numpy as np
from astropy.coordinates import Angle, SkyCoord
from gwpy.timeseries import TimeSeries
import astropy.units as u
import astropy.constants as ac

from gwsignal.core.waveform import CompactBinaryCoalescenceGenerator
import gwsignal.core.gw as gw
from gwsignal.core.utils import add_params_units


def modes_to_k(modes: list[tuple[int, int]]) -> list[int]:
    """TEOBResumS-specific function for converting
    a list of modes expressed as (l, m) into a list of single
    integers k.
    """
    return [int(x[0] * (x[0] - 1) / 2 + x[1] - 2) for x in modes]


TEOB_DALI_MODES = []
for l in range(2, 6):
    for m in range(1, l + 1):
        TEOB_DALI_MODES.append((l, m))

# reversed dictionary to back-reference the modes
TEOB_DALI_MODES_FROM_K = {modes_to_k([mode])[0]: mode for mode in TEOB_DALI_MODES}


class TEOBResumSDALI(CompactBinaryCoalescenceGenerator):

    """
    Implements a wrapper for the eccentric, time domain
    EOB model TEOBResumS-DALI.
    """

    def __init__(self, modes_to_use: Optional[list[tuple[int, int]]] = None, **kwargs):

        super().__init__()
        if modes_to_use is None:
            self.available_modes = TEOB_DALI_MODES
        else:
            self.available_modes = modes_to_use
        self._update_domains()

    @property
    def metadata(self):
        metadata = {
            "type": "eccentric",
            "f_ref_spin": False,
            "modes": True,
            "polarizations": True,
            "implemented_domain": "time",
            "approximant": "TEOBResumS",
            "implementation": "C",
            "condition": False,
            "conditioning_routines": None,
        }
        return metadata

    @property
    def _available_modes_teob_convention(self):
        return modes_to_k(self.available_modes)

    def generate_td_modes(self, **parameters):
        self.parameter_check(units_sys="Cosmo", **parameters)
        self.waveform_dict = self._strip_units(self.waveform_dict)

        parameters = convert_parameters_to_teob(self.waveform_dict)

        parameters["use_mode_lm"] = self._available_modes_teob_convention

        t, hp, hc, htlm, dyn = EOBRun_module.EOBRunPy(parameters)

        hlm_reindexed = self._to_gwpy_series(
            {
                TEOB_DALI_MODES_FROM_K[int(ind)]:
                mode_wf[0]
                * np.exp(-1j * mode_wf[1])
                for ind, mode_wf in htlm.items()
            },
            t,
        )

        return gw.GravitationalWaveModes(hlm_reindexed)

    def generate_td_waveform(self, **parameters):
        theta, phi = parameters["inclination"], parameters["phi_ref"]
        hlm = self.generate_td_modes(**parameters)
        hp, hc = hlm(theta, phi)
        
        nu = compute_symmetric_mass_ratio(self.waveform_dict)
        
        distance_rescaling = (
            (
                nu
                * (parameters["mass1"] + parameters["mass2"])
                / parameters["distance"]
                * ac.G
                / ac.c ** 2
            )
            .to(u.dimensionless_unscaled)
            .value
        )

        hp = hp * distance_rescaling
        hc = hc * distance_rescaling

        hp, hc = TimeSeries(hp, name="hplus"), TimeSeries(hc, name="hcross")
        return hp, hc

    def _to_gwpy_series(self, modes_dict, times):
        """
        Iterate over the dict and return a dict of gwpy TimeSeries objects.
        TEOBResumS only outputs the positive-m modes, but in order
        to correctly recover the signal we also need the negative-m ones,
        which can be obtained by symmetry.
        """
        gwpy_dict = {}
        for ellm, mode in modes_dict.items():
            l, m = ellm
            gwpy_dict[ellm] = TimeSeries(mode, times=times, name=f"h_{l}_{m}")
            gwpy_dict[(l, -m)] = TimeSeries(
                (-1) ** l * np.conj(mode), times=times, name=f"h_{l}_-{m}"
            )
        return gwpy_dict

    def _strip_units(self, waveform_dict):
        new_dc = {}
        for key in waveform_dict.keys():
            if isinstance(waveform_dict[key], u.Quantity):
                new_dc[key] = waveform_dict[key].value
            else:
                new_dc[key] = waveform_dict[key]
        return new_dc


def convert_parameters_to_teob(parameter_dict):

    fmin, dt = parameter_dict["f22_start"], parameter_dict["deltaT"]
    m1, m2 = parameter_dict["mass1"], parameter_dict["mass2"]

    q = m1 / m2
    if q < 1.0:
        q = 1 / q

    return {
        "M": (m1 + m2),
        "q": q,
        "chi1x": parameter_dict["spin1x"],
        "chi1y": parameter_dict["spin1y"],
        "chi1z": parameter_dict["spin1z"],
        "chi2x": parameter_dict["spin2x"],
        "chi2y": parameter_dict["spin2y"],
        "chi2z": parameter_dict["spin2z"],
        "ecc": parameter_dict["eccentricity"],
        "inclination": parameter_dict["inclination"],
        "coalescence_angle": np.pi / 2 - parameter_dict["phi_ref"],
        "srate_interp": 1 / dt,
        "initial_frequency": fmin,
        "distance": parameter_dict["distance"],
        "arg_out": "yes",
        "interp_uniform_grid": "yes",
        "use_geometric_units": "no",
    }


def compute_symmetric_mass_ratio(parameter_dict):
    m1, m2 = parameter_dict["mass1"], parameter_dict["mass2"]

    return m1 * m2 / (m1 + m2) ** 2
