try:
    import pyseobnr.generate_waveform as pyseob_wf
except ImportError:
    print("The pyseobnr package has failed to load, you will not be able to employ SEOBNRv5 approximants.")

from numpy import arange
from astropy.units import Mpc
from gwpy.timeseries import TimeSeries
from gwpy.frequencyseries import FrequencySeries
from lal import LIGOTimeGPS

from ..core.waveform import CompactBinaryCoalescenceGenerator


class SEOBNRv5EHM(CompactBinaryCoalescenceGenerator):

    """
    Implements a wrapper for SEOBNRv5EHM in the pyseob package

    Parameters
    ----------

        No parameters required for initialization.

    """

    # _domain = 'time'
    def __init__(self, **kwargs):

        super().__init__()

        self.seobnr = pyseob_wf
        self._domain = "time"
        self._implemented_domain = "time"
        self._generation_domain = None

    @property
    def metadata(self):
        metadata = {
            "type": "aligned_spin",
            "f_ref_spin": True,
            "modes": True,
            "polarizations": True,
            "implemented_domain": "time",
            "approximant": "SEOBNRv5EHM",
            "implementation": "",
            "conditioning_routines": "",
        }
        return metadata

    def _strip_units(self, waveform_dict):
        new_dc = {}
        for key in waveform_dict.keys():
            new_dc[key] = waveform_dict[key].value
        return new_dc

    def _generate_waveform_class(self, **parameters):

        parameters_extra = {}

        # Get the optional parameters for this approximant, if any
        optional_params = [
            "ModeArray",
            "condition",
            "lmax_nyquist",
            "lmax",
            "secular_bwd_int",
            "warning_secular_bwd_int",
            "t_backwards",
            "warning_bwd_int",
        ]

        for key in optional_params:
            val = parameters.pop(key, None)
            if val is not None:
                parameters_extra[key] = val

        # indicates we are running pyseobnr through GWSignal. Emits proper
        # warnings on non-reviewed model calls.
        parameters_extra["gwsignal_environment"] = True

        self.parameter_check(units_sys="Cosmo", **parameters)
        self.waveform_dict["distance"] = self.waveform_dict["distance"].to(Mpc)
        self.waveform_dict = self._strip_units(self.waveform_dict)
        self.waveform_dict["approximant"] = self.metadata["approximant"]
        self.waveform_dict["f_ref"] = self.waveform_dict["f22_ref"]
        self.waveform_dict["rel_anomaly"] = self.waveform_dict["meanPerAno"]

        # Make sure to update the waveform dictionary with extra parameters
        self.waveform_dict.update(**parameters_extra)

        #print(f"waveform_dict = {self.waveform_dict}")
        return self.seobnr.GenerateWaveform(self.waveform_dict)

    def _generate_modes(self, **parameters):

        gen_wf = self._generate_waveform_class(**parameters)
        times, hlm = gen_wf.generate_td_modes()
        epoch = LIGOTimeGPS(times[0])
        dt = self.waveform_dict["deltaT"]

        hlm_dict = {}
        for k, v in hlm.items():
            hlm_lal = TimeSeries(v, times=times,  name=k)
            hlm_dict[k] = hlm_lal

        return hlm_dict

    def _generate_td_polarizations(self, **parameters):

        gen_wf = self._generate_waveform_class(**parameters)
        if self.waveform_dict.get("condition"):
            hp, hc = gen_wf.generate_td_polarizations_conditioned_1()
        else:
            hp, hc = gen_wf.generate_td_polarizations()
        epoch = hp.epoch.gpsSeconds + hp.epoch.gpsNanoSeconds / 1e9
        times = hp.deltaT * arange(hp.data.length) + epoch
        return (
            TimeSeries(hp.data.data, times=times, name="hp"),
            TimeSeries(hc.data.data, times=times, name="hc"),
        )

    def _generate_fd_polarizations_from_td(self, **parameters):

        gen_wf = self._generate_waveform_class(**parameters)
        hptilde, hctilde = gen_wf.generate_fd_polarizations()
        frequencies = hptilde.deltaF * arange(hptilde.data.length)
        epoch = hptilde.epoch.gpsSeconds + hptilde.epoch.gpsNanoSeconds / 1e9
        return (
            FrequencySeries(
                hptilde.data.data, frequencies=frequencies, epoch=epoch, name="hp"
            ),
            FrequencySeries(
                hctilde.data.data, frequencies=frequencies, epoch=epoch, name="hc"
            ),
        )

    def _generate_polarizations(self, **parameters):

        if self._generation_domain == "time":
            return self._generate_td_polarizations(**parameters)
        elif self._generation_domain == "freq":
            return self._generate_fd_polarizations_from_td(**parameters)
        else:
            raise ValueError("Generation domain must be 'time' or 'freq'.")

    def generate_fd_waveform(self, **parameters):
        return self._generate_fd_polarizations_from_td(**parameters)

    def generate_td_waveform(self, **parameters):
        return self._generate_td_polarizations(**parameters)

    def generate_td_modes(self, **parameters):
        return self._generate_modes(**parameters)


class SEOBNRv5HM(CompactBinaryCoalescenceGenerator):

    """
    Implements a wrapper for SEOBNRv5 in the pyseob package

    Parameters
    ----------

        No parameters required for initialization.

    """

    # _domain = 'time'
    def __init__(self, **kwargs):

        # super().__init__()

        self.seobnr = pyseob_wf
        self._domain = "time"
        self._implemented_domain = "time"
        self._generation_domain = None

    @property
    def metadata(self):
        metadata = {
            "type": "aligned_spin",
            "f_ref_spin": True,
            "modes": True,
            "polarizations": True,
            "implemented_domain": "time",
            "approximant": "SEOBNRv5HM",
            "implementation": "",
            "conditioning_routines": "",
        }
        return metadata

    def _strip_units(self, waveform_dict):
        new_dc = {}
        for key in waveform_dict.keys():
            new_dc[key] = waveform_dict[key].value
        return new_dc

    def _generate_waveform_class(self, **parameters):

        parameters_extra = {}

        # Get the optional parameters for this approximant, if any
        optional_params = [
            "ModeArray",
            "postadiabatic",
            "postadiabatic_type",
            "condition",
            "lmax_nyquist",
            "lmax",
            "dA_dict",
            "dw_dict",
            "dTpeak",
            "domega_dict",
            "dtau_dict",
            "da6",
            "ddSO",
            "tol_PA",
            "rtol_ode",
            "atol_ode",
            "deltaT_sampling",
        ]
        for key in optional_params:
            val = parameters.pop(key, None)
            if val is not None:
                parameters_extra[key] = val

        # indicates we are running pyseobnr through GWSignal. Emits proper
        # warnings on non-reviewed model calls.
        parameters_extra["gwsignal_environment"] = True

        self.parameter_check(units_sys="Cosmo", **parameters)
        self.waveform_dict["distance"] = self.waveform_dict["distance"].to(Mpc)
        self.waveform_dict = self._strip_units(self.waveform_dict)
        self.waveform_dict["approximant"] = self.metadata["approximant"]
        self.waveform_dict["f_ref"] = self.waveform_dict["f22_ref"]

        # Make sure to update the waveform dictionary with extra parameters
        self.waveform_dict.update(**parameters_extra)

        return self.seobnr.GenerateWaveform(self.waveform_dict)

    def _generate_modes(self, **parameters):

        gen_wf = self._generate_waveform_class(**parameters)
        times, hlm = gen_wf.generate_td_modes()
        epoch = LIGOTimeGPS(times[0])
        dt = self.waveform_dict["deltaT"]

        hlm_dict = {}
        for k, v in hlm.items():
            hlm_lal = TimeSeries(v, times=times,  name=k)
            hlm_dict[k] = hlm_lal

        return hlm_dict

    def _generate_td_polarizations(self, **parameters):

        gen_wf = self._generate_waveform_class(**parameters)
        if self.waveform_dict.get("condition"):
            hp, hc = gen_wf.generate_td_polarizations_conditioned_2()
        else:
            hp, hc = gen_wf.generate_td_polarizations()
        epoch = hp.epoch.gpsSeconds + hp.epoch.gpsNanoSeconds / 1e9
        times = hp.deltaT * arange(hp.data.length) + epoch
        return (
            TimeSeries(hp.data.data, times=times, name="hp"),
            TimeSeries(hc.data.data, times=times, name="hc"),
        )

    def _generate_fd_polarizations_from_td(self, **parameters):

        gen_wf = self._generate_waveform_class(**parameters)
        hptilde, hctilde = gen_wf.generate_fd_polarizations()
        frequencies = hptilde.deltaF * arange(hptilde.data.length)
        epoch = hptilde.epoch.gpsSeconds + hptilde.epoch.gpsNanoSeconds / 1e9
        return (
            FrequencySeries(
                hptilde.data.data, frequencies=frequencies, epoch=epoch, name="hp"
            ),
            FrequencySeries(
                hctilde.data.data, frequencies=frequencies, epoch=epoch, name="hc"
            ),
        )

    def _generate_polarizations(self, **parameters):

        if self._generation_domain == "time":
            return self._generate_td_polarizations(**parameters)
        elif self._generation_domain == "freq":
            return self._generate_fd_polarizations_from_td(**parameters)
        else:
            raise ValueError("Generation domain must be 'time' or 'freq'.")

    def generate_fd_waveform(self, **parameters):
        return self._generate_fd_polarizations_from_td(**parameters)

    def generate_td_waveform(self, **parameters):
        return self._generate_td_polarizations(**parameters)

    def generate_td_modes(self, **parameters):
        return self._generate_modes(**parameters)


class SEOBNRv5PHM(CompactBinaryCoalescenceGenerator):

    """
    Implements a wrapper for SEOBNRv5 in the pyseob package

    Parameters
    ----------

        No parameters required for initialization.

    """

    # _domain = 'time'
    def __init__(self, **kwargs):

        # super().__init__()

        self.seobnr = pyseob_wf
        self._domain = "time"
        self._implemented_domain = "time"
        self._generation_domain = None

    @property
    def metadata(self):
        metadata = {
            "type": "precessing_spin",
            "f_ref_spin": True,
            "modes": True,
            "polarizations": True,
            "implemented_domain": "time",
            "approximant": "SEOBNRv5PHM",
            "implementation": "",
            "conditioning_routines": "",
        }
        return metadata

    def _strip_units(self, waveform_dict):
        new_dc = {}
        for key in waveform_dict.keys():
            new_dc[key] = waveform_dict[key].value
        return new_dc

    def _generate_waveform_class(self, **parameters):

        parameters_extra = {}

        # Get the optional parameters for this approximant, if any
        optional_params = [
            "ModeArray",
            "polarizations_from_coprec",
            "postadiabatic",
            "postadiabatic_type",
            "condition",
            "lmax_nyquist",
            "lmax",
            "dA_dict",
            "dw_dict",
            "dTpeak",
            "domega_dict",
            "dtau_dict",
            "da6",
            "ddSO",
            "deltaT_sampling",
            "omega_prec_deviation",
            "enable_antisymmetric_modes",
            "antisymmetric_modes_hm",
            "antisymmetric_modes",
        ]
        for key in optional_params:
            val = parameters.pop(key, None)
            if val is not None:
                parameters_extra[key] = val

        # indicates we are running pyseobnr through GWSignal. Emits proper
        # warnings on non-reviewed model calls.
        parameters_extra["gwsignal_environment"] = True

        self.parameter_check(units_sys="Cosmo", **parameters)
        self.waveform_dict["distance"] = self.waveform_dict["distance"].to(Mpc)
        self.waveform_dict = self._strip_units(self.waveform_dict)
        self.waveform_dict["approximant"] = self.metadata["approximant"]
        self.waveform_dict["f_ref"] = self.waveform_dict["f22_ref"]

        # Make sure to update the waveform dictionary with extra parameters
        self.waveform_dict.update(**parameters_extra)

        return self.seobnr.GenerateWaveform(self.waveform_dict)

    def _generate_modes(self, **parameters):

        gen_wf = self._generate_waveform_class(**parameters)
        times, hlm = gen_wf.generate_td_modes()
        epoch = LIGOTimeGPS(times[0])
        dt = self.waveform_dict["deltaT"]

        hlm_dict = {}
        for k, v in hlm.items():
            hlm_lal = TimeSeries(v, times=times,  name=k)
            hlm_dict[k] = hlm_lal

        return hlm_dict

    def _generate_td_polarizations(self, **parameters):

        gen_wf = self._generate_waveform_class(**parameters)
        if self.waveform_dict.get("condition"):
            hp, hc = gen_wf.generate_td_polarizations_conditioned_2()
        else:
            hp, hc = gen_wf.generate_td_polarizations()
        epoch = hp.epoch.gpsSeconds + hp.epoch.gpsNanoSeconds / 1e9
        times = hp.deltaT * arange(hp.data.length) + epoch
        return (
            TimeSeries(hp.data.data, times=times, name="hp"),
            TimeSeries(hc.data.data, times=times, name="hc"),
        )

    def _generate_fd_polarizations_from_td(self, **parameters):

        gen_wf = self._generate_waveform_class(**parameters)
        hptilde, hctilde = gen_wf.generate_fd_polarizations()
        frequencies = hptilde.deltaF * arange(hptilde.data.length)
        epoch = hptilde.epoch.gpsSeconds + hptilde.epoch.gpsNanoSeconds / 1e9
        return (
            FrequencySeries(
                hptilde.data.data, frequencies=frequencies, epoch=epoch, name="hp"
            ),
            FrequencySeries(
                hctilde.data.data, frequencies=frequencies, epoch=epoch, name="hc"
            ),
        )

    def _generate_polarizations(self, **parameters):

        if self._generation_domain == "time":
            return self._generate_td_polarizations(**parameters)
        elif self._generation_domain == "freq":
            return self._generate_fd_polarizations_from_td(**parameters)
        else:
            raise ValueError("Generation domain must be 'time' or 'freq'.")

    def generate_fd_waveform(self, **parameters):
        return self._generate_fd_polarizations_from_td(**parameters)

    def generate_td_waveform(self, **parameters):
        return self._generate_td_polarizations(**parameters)

    def generate_td_modes(self, **parameters):
        return self._generate_modes(**parameters)
