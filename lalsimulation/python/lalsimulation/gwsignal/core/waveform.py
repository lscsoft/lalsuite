import lal
import lalsimulation as lalsim
import numpy as np
from astropy import units as u
from gwpy.timeseries import TimeSeries
from gwpy.frequencyseries import FrequencySeries

from . import utils as ut
from . import parameter_conventions as pc
from . import gw as gw
from . import waveform_conditioning as wave_cond
from . import errors as err
import warnings
###################################################################################################


class GravitationalWaveGenerator(object):
    """
    This is the Parent class for all gravitational wave generator classes.
    This class implements the basic generate modes and polarizations definitions which can then be used by the
    consequent children classes.
    """

    def __init__(self):
        """
        Initialize the class and set base domain to time.

        Parameters
        ----------

            No parameters required for initialization.
        """

    def generate_td_modes(self, **kwargs):
        """
        Generate time domain modes given kwargs. The child classes will provide these routines.

        Parameters
        ----------

            kwargs required for generation of waveform modes

        Returns
        -------

            Waveform modes as implemented in the child class

        """
        raise NotImplementedError

    def generate_fd_modes(self, **kwargs):
        """
        Generate frequency domain modes given kwargs. The child classes will provide these routines.

        Parameters
        ----------

            kwargs required for generation of waveform modes

        Returns
        -------

            Waveform modes as implemented in the child class

        """
        raise NotImplementedError

    def generate_td_waveform(self, **kwargs):
        """
        Generate time domain polarizations given kwargs. The child classes will provide these routines.

        Parameters
        ----------

            kwargs required for generation of waveform polarizations

        Returns
        -------

            Waveform polarizations implemented in the child class

        """
        raise NotImplementedError

    def generate_fd_waveform(self, **kwargs):
        """
        Generate frequency domain polarizations given kwargs. The child classes will provide these routines.

        Parameters
        ----------

            kwargs required for generation of waveform polarizations

        Returns
        -------

            Waveform polarizations as implemented in the child class

        """
        raise NotImplementedError

    @property
    def domain(self):
        return self._generation_domain, self._implemented_domain

    def _update_domains(self):
        """
        Update the generation or implementation domain for a given generator / if needed.

        Parameters
        ----------
            None

        Returns
        -------
            Updated domain tags for polarization generation and waveform implementation.

        """
        self._implemented_domain = self.metadata['implemented_domain']
        return 0

    @property
    def metadata(self):
        metadata = {
            "type": '',
            "f_ref_spin": None,
            "modes": None,
            "polarizations": None,
            "implemented_domain": None,
            "generation_domain": None,
            "approximant" : None,
            "implementation" : '',
            "conditioning_routines" : ''
        }
        return metadata

class CompactBinaryCoalescenceGenerator(GravitationalWaveGenerator):
    """
        This is the parent generator class for compact binary coalescence waveforms (BBH, BNS, NSBH etc.).

        The main intended functionality of this class is to check the different parameters passed in the kwargs
        and to populate the parameter dictionaries with all the other parameters from the lalsuite parameter checks.
    """

    def __init__(self):
        """
        Initialize class

        Parameters
        ----------

            No parameters required for initilazation.

        """
        super(CompactBinaryCoalescenceGenerator, self).__init__()


    def parameter_check(self, units_sys='S.I.', extra_parameters=dict(), **parameters):
        """
        Perform checks on the various parameters and populate the different parameters not passed
        in the kwargs required for generating waveforms.

        Parameters
        ----------

            Python dictionary of parameters required for waveform generation
            of the form specified in `parameter_conventions.py`.

        Returns
        -------

            Populate self.waveform_dict with python dictionary and self.lal_dict with LALSuite dictionary structure.
        """
        default_dict = pc.default_dict.copy()
        # Need to add this line to take care of the extra parameters passed to ExternalPython LAL Generator
        ExternalPythonParameters=['object', 'module']

        for key, value in parameters.items():
            if key in ExternalPythonParameters:
                pass
            else:
                default_dict[key] = value

        if not 'deltaF' in default_dict:
            default_dict['deltaF'] = 1./16.*u.Hz
        if not 'deltaT' in default_dict:
            default_dict['deltaT'] = 1./512.*u.s
        if not 'f_max' in default_dict:
            default_dict['f_max'] = (0.5/default_dict['deltaT'].value)*u.Hz

        # Add units if indicated
        if units_sys is not None:
            self.waveform_dict = ut.add_params_units(default_dict, units_sys, generic_param_dict=extra_parameters)
        #This is a mandatory check that units are correctly used
        ut.check_dict_parameters(default_dict, generic_param_dict=extra_parameters)

        self.lal_dict  = ut.to_lal_dict(default_dict)
        return default_dict

########################################################################
#  LAL CBC Generator
########################################################################


class LALCompactBinaryCoalescenceGenerator(CompactBinaryCoalescenceGenerator):
    """
    Generator class for all CBC waveforms as implemented in LALSimulation.
    """

    def __init__(self, approximant, **kwargs):
        """
        Initialize class. Parent class is the "CompactBinaryCoalescenceGenerator" class

        Parameters
        ----------
        approximant : type 'str'
                    Name of the LAL approximant of which waveform to be generated

        Returns
        -------
        generator : Python GW generator class
                  Generator class for LALSimulation waveforms required to generate time and frequency domain polarizations and modes

        """
        super(LALCompactBinaryCoalescenceGenerator, self).__init__()
        warnings.warn("This code is currently UNREVIEWED, use with caution!")
        self.approximant = approximant
        self._get_approximant(self.approximant)
        self._generation_domain = None
        self._update_domains()

    @property
    def metadata(self):
        metadata = {
            "type": 'cbc_lalsimulation',
            "f_ref_spin": True,
            "modes": True,
            "polarizations": True,
            "implemented_domain": self._implemented_domain,
            "generation_domain": self._generation_domain,
            "approximant" : self._approx_name,
            "implementation" : "LALSimulation",
            "conditioning_routines" : 'lalsimulation'
        }
        return metadata

    def _update_domains(self):
        """
        Update the implemented domain of the LAL generator based on LALSimulations ImplementedTD/FDApproximant flag
        """
        # Check which is the impltemented domain for the waveform

        td_flag = lalsim.SimInspiralImplementedTDApproximants(self._approx)
        fd_flag = lalsim.SimInspiralImplementedFDApproximants(self._approx)
        check_array = [td_flag, fd_flag]

        if check_array==[1,0]:
            self._implemented_domain = 'time'
        elif check_array==[0,1]:
            self._implemented_domain = 'freq'
        elif check_array==[1,1]:
            self._implemented_domain = self._generation_domain

    def _get_approximant(self, approx):
        """Populate the approximant name and approximant int (as used by lalsimulation) from waveform dictionary.

        Parameters
        ----------
        approx : 'str'
            Name of approximant.

        Returns
        -------

            Populate self._approx_name (str) and self._approx (int) to be used by lalsimulation ChooseWaveformGenerator
            function

        Raises
        ------
            ValueError if approximant is not included in LALSimulation
        """

        if type(approx) is str:
            self._approx_name = approx
            self._approx = lalsim.SimInspiralGetApproximantFromString(approx)
        elif type(approx) is int:
            self._approx = approx
            self._approx_name = lalsim.GetStringFromApproximant(approx)
        else:
            raise ValueError('approximant not of recognized type')


    def generate_td_waveform(self, **parameters):
        """
        Perform parameter check, choose lalsimulation generator based on domain and conditioning subroutines.

        Parameters
        ----------
        parameter_dict : dictionary
            Dictionary of waveform parameters of the form specified
            in `parameter_conventions.py`.

        Returns
        -------
        hp, hc : LAL Time Series
            Plus and cross polarizations of a gravitational waveform (hp,hc) as LAL Data objects.

        Raises
        ------
            ValueError if domain ('time' or 'freq') for approximant is not specified
        """
        self.parameter_check(**parameters)
        self._generation_domain = 'time'
        self._update_domains()
        if self._implemented_domain=='time':
            self._pol_gen_function = lalsim.SimInspiralGenerateTDWaveform

        elif self._implemented_domain=='freq':
            if self._generation_domain=='time' and parameters['condition']==1:
                self._pol_gen_function = lalsim.SimInspiralGenerateTDWaveform
            elif self._generation_domain=='time' and parameters['condition']==0:
                raise ValueError("Generator requires conditioning to be turned on to generate time domain waveform")

        self._lal_generator = lalsim.SimInspiralChooseGenerator(self._approx, self.lal_dict)

        @err.mapexception
        def gen_pol_func(lal_dict, lal_generator):
            return self._pol_gen_function(lal_dict, lal_generator)

        hp, hc = gen_pol_func(self.lal_dict, self._lal_generator)
        hp, hc = to_gwpy_Series(hp, name='hplus'), to_gwpy_Series(hc, name='hcross')
        return hp, hc

    def generate_fd_waveform(self, **parameters):
        """
        Perform parameter check, choose lalsimulation generator based on domain and conditioning subroutines.

        Parameters
        ----------
        parameter_dict : dictionary
            Dictionary of waveform parameters of the form specified
            in `parameter_conventions.py`.

        Returns
        -------
        hp, hc : LAL Frequency Series
            Plus and cross polarizations of a gravitational waveform (hp,hc) as LAL Data objects.

        Raises
        ------
            ValueError if domain ('time' or 'freq') for approximant is not specified
        """
        self.parameter_check(**parameters)
        self._generation_domain = 'freq'
        self._update_domains()

        if self._implemented_domain=='time':
            if self._generation_domain=='freq' and parameters['condition']==1:
                self._pol_gen_function = lalsim.SimInspiralGenerateFDWaveform
            elif self._generation_domain=='freq' and parameters['condition']==0:
                raise ValueError("Generator requires conditioning to be turned on to generate frequency domain waveform")

        elif self._implemented_domain=='freq':
            self._pol_gen_function = lalsim.SimInspiralGenerateFDWaveform

        self._lal_generator = lalsim.SimInspiralChooseGenerator(self._approx, self.lal_dict)

        @err.mapexception
        def gen_pol_func(lal_dict, lal_generator):
            return self._pol_gen_function(lal_dict, lal_generator)

        hp, hc = gen_pol_func(self.lal_dict, self._lal_generator)
        hp, hc = to_gwpy_Series(hp, name='hplus', epoch=0.), to_gwpy_Series(hc, name='hcross', epoch=0.)
        return hp, hc


    def generate_td_modes(self, **parameters):
        """
        Perform parameter check, choose lalsimulation generator based on domain and attempt to generate
        the waveform modes.

        Parameters
        ----------
        parameter_dict : dictionary
            Dictionary of waveform parameters of the form specified
            in `parameter_conventions.py`.

        Returns
        -------
        hlm : Python dictionary of modes
            Modes of a gravitational waveform as python dictionary where each mode is LAL Time series
            objects.
        """

        self.parameter_check(**parameters)
        self._generation_domain = 'time'
        self._update_domains()

        if self._implemented_domain=='time':
            self._mode_gen_function = lalsim.SimInspiralGenerateTDModes
        elif self._implemented_domain=='freq':
            raise ValueError('Approximant does not have time domain mode generator')
        else:
            raise ValueError('Approximant domain unspecified')


        self._lal_generator = lalsim.SimInspiralChooseGenerator(self._approx, self.lal_dict)

        @err.mapexception
        def gen_modes(lal_dict, gen):
            return self._mode_gen_function(lal_dict, gen)


        hlm = gen_modes(self.lal_dict, self._lal_generator)


        # Put hlms in a dictionary
        # Define python dictionary and populate the modes
        hlm_dict = {}

        # Define 22 mode to get epoch and generate time data
        mode_22 = lalsim.SphHarmTimeSeriesGetMode(hlm, 2, 2)
        dt = self.waveform_dict['deltaT'].si.value
        times = np.float64(mode_22.epoch) + np.arange(0, len(mode_22.data.data))*dt

        hlm_dict['time_array'] = np.array(times)

        lmax = hlm.l
        m_max = hlm.m
        for l in np.arange(2, lmax+1):
            for m in np.arange(-l, l+1):
                hlm_dict[l,m] = lalsim.SphHarmTimeSeriesGetMode(hlm, int(l), int(m))

        hlm_out = to_gwpy_dict(hlm_dict)
        return hlm_out


    def generate_fd_modes(self, **parameters):
        """
        Perform parameter check, choose lalsimulation generator based on domain and attempt to generate
        the waveform modes.

        Parameters
        ----------
        parameter_dict : dictionary
            Dictionary of waveform parameters of the form specified
            in `parameter_conventions.py`.

        Returns
        -------
        hlm : Python dictionary of modes
            Modes of a gravitational waveform as python dictionary where each mode is LAL Frequency series
            objects.
        """

        self.parameter_check(**parameters)
        self._generation_domain = 'freq'
        self._update_domains()

        if self._implemented_domain=='time':
            raise ValueError('Approximant does not have frequency domain mode generator')
        elif self._implemented_domain=='freq':
            self._mode_gen_function = lalsim.SimInspiralGenerateFDModes
        else:
            raise ValueError('approximant domain unspecified')


        self._lal_generator = lalsim.SimInspiralChooseGenerator(self._approx, self.lal_dict)
        @err.mapexception
        def gen_modes(lal_dict, gen):
            return self._mode_gen_function(lal_dict, gen)

        hlm = gen_modes(self.lal_dict, self._lal_generator)


        # Put hlms in a dictionary
        # Define python dictionary and populate the modes
        hlm_dict = {}

        hlm_dict['frequency_array'] = hlm.fdata.data

        lmax = hlm.l
        m_max = hlm.m
        for l in np.arange(2, lmax+1):
            for m in np.arange(-l, l+1):
                hlm_dict[l,m] = lalsim.SphHarmFrequencySeriesGetMode(hlm, int(l), int(m))

        hlm_out = to_gwpy_dict(hlm_dict)

        return hlm_out

def conditioning_generator(generator):
    """
    Given a generator with contioning as 1, return function that generates conditioned waveforms

    Parameters
    ----------
    generator :  `GravitationalWaveGenerator`
        GravitationalWaveGenerator object

    Returns
    -------
        Conditioning function for the generator
    """

    if generator._implemented_domain=='time' and generator._generation_domain=='time':
        generator_func = wave_cond.generate_conditioned_td_waveform_from_td

    elif generator._implemented_domain=='time' and generator._generation_domain=='freq':
        generator_func = wave_cond.generate_conditioned_fd_waveform_from_td

    elif generator._implemented_domain=='freq' and generator._generation_domain=='freq':
        generator_func = wave_cond.generate_conditioned_fd_waveform_from_fd

    elif generator._implemented_domain=='freq' and generator._generation_domain=='time':
        generator_func = wave_cond.generate_conditioned_td_waveform_from_fd

    return generator_func

########################################################################


########################################################################

def to_gwpy_Series(h, f0=0., **kwargs):
    '''
    Function to convert a lal series to a gwpy series.

    Parameters
    ----------
    h : lal time or frequency series
        'lal.REAL8TimeSeries' or 'lal.COMPLEX16FrequencySeries'

    f0 : Starting frequency passed to gwpy.FrequencySeries

    Returns
    -------

    h : GWpy Time/Freq series
            Time/Freq domain waveform
    '''

    if isinstance(h, lal.REAL8TimeSeries) or isinstance(h, lal.COMPLEX16TimeSeries):
        return TimeSeries(h.data.data, dt = h.deltaT, t0 = h.epoch, **kwargs)

    elif isinstance(h, lal.COMPLEX16FrequencySeries):
        return FrequencySeries(h.data.data, df = h.deltaF, f0=f0, **kwargs)

    elif isinstance(h, TimeSeries) :
        return TimeSeries(h, **kwargs)

    elif isinstance(h, FrequencySeries):
        return FrequencySeries(h, **kwargs)

    else:
        print('Input type not recognized')

def to_gwpy_dict(mode_dict, **kwargs):
    '''
    Function to convert a mode dictionary of lal series to a dictionary of
    gwpy series.

    Parameters
    ----------
    mode_dict : python dictionary
        entries of the form {(l, m): lal-series} and {`time/frequency_array` : }

    Returns
    -------
    new_dict : python dictionary
            entries of the form {(l, m): gwpy-series} and {`time/frequency_array` : }
    '''

    new_dict = {}
    f0=0.
    for k, v in mode_dict.items():
        if k == 'time_array':
           new_dict.update({k: v})
        elif k == 'frequency_array':
           new_dict.update({k: v})
           f0 = v[0]
        elif v is None:
           pass
        else:
           new_dict.update({k: to_gwpy_Series(v, f0, name='h_%i_%i'%(k[0], k[1]),**kwargs)})

    return new_dict


########################################################################################################################
################################  Generate Waveform Polarizations and Modes functions ##################################
########################################################################################################################

def GenerateTDWaveform(parameter_dict, generator):
    """
    Function to generate time domain gravitational wave polarizations.

    Parameters
    ----------
    parameter_dict : dictionary
            Dictionary of intrinsic / extrinsic gravitational wave parameters
    generator : Python GW generator class
            Generator class for the waveform approximant (either LAL or external)

    Returns
    -------
        `GravitationalWavePolarizations` object: allows to extract the zero noise strain
        given a detector, sky-position, polarization and time of arrival values.
    """

    generator._generation_domain = 'time'
    generator._update_domains()

    conditioning = parameter_dict['condition']
    conditioning_routines = generator.metadata['conditioning_routines']

    # First check if generator generates waveform in generation domain and
    # if conditioning is on. If conditioning is off, raise ValueError
    if not conditioning and generator._generation_domain!=generator._implemented_domain:
        raise ValueError("Generator does not provide a method to generate time-domain waveforms. \n Please turn on conditioning to generate time-domain waveform")

    if conditioning and conditioning_routines=='gwsignal':
        warnings.warn("This code is currently UNREVIEWED, use with caution!")
        generator_func = conditioning_generator(generator)
        hp, hc = generator_func(parameter_dict, generator)
    else:
        hp, hc = generator.generate_td_waveform(**parameter_dict)

    return gw.GravitationalWavePolarizations(hp, hc)

def GenerateFDWaveform(parameter_dict, generator):
    """
    Function to generate frequency domain gravitational wave polarizations.

    Parameters
    ----------
    parameter_dict : dictionary
            Dictionary of intrinsic / extrinsic gravitational wave parameters
    generator : Python GW generator class
            Generator class for the waveform approximant (either LAL or external)

    Returns
    -------
        `GravitationalWavePolarizations` object: allows to extract the zero noise strain
        given a detector, sky-position, polarization and time of arrival values.
    """

    generator._generation_domain = 'freq'
    generator._update_domains()

    conditioning = parameter_dict['condition']
    conditioning_routines = generator.metadata['conditioning_routines']

    # First check if generator generates waveform in generation domain and
    # if conditioning is on. If conditioning is off, raise ValueError
    if not conditioning and generator._generation_domain!=generator._implemented_domain:
        raise ValueError("Generator does not provide a method to generate frequency-domain waveforms. \n Please turn on conditioning to generate time-domain waveform")

    if conditioning and conditioning_routines=='gwsignal':
        warnings.warn("This code is currently UNREVIEWED, use with caution!")
        generator_func = conditioning_generator(generator)
        hp, hc = generator_func(parameter_dict, generator)
    else:
        hp, hc = generator.generate_fd_waveform(**parameter_dict)

    return gw.GravitationalWavePolarizations(hp, hc)

def GenerateTDModes(parameter_dict, generator):
    """
    Function to generate time domain gravitational wave modes.

    Parameters
    ----------
    parameter_dict : dictionary
            Dictionary of intrinsic / extrinsic gravitational wave parameters
    generator : Python GW generator class
            Generator class for the waveform approximant (either LAL or external)

    Returns
    -------
    hlm    : Python dictionary
            Time domain modes returned as a dictionary where each mode is a GWpy Time series object
    """

    generator._generation_domain = 'time'
    generator._update_domains()

    hlm = generator.generate_td_modes(**parameter_dict)

    hlm_out = to_gwpy_dict(hlm)
    return gw.GravitationalWaveModes(hlm_out)


def GenerateFDModes(parameter_dict, generator):
    """
    Function to generate frequency domain gravitational wave modes.

    Parameters
    ----------
    parameter_dict : dictionary
            Dictionary of intrinsic / extrinsic gravitational wave parameters
    generator : Python GW generator class
            Generator class for the waveform approximant (either LAL or external)

    Returns
    -------
    hlm    : Python dictionary
            Frequency domain modes returned as a dictionary where each mode is a GWpy Frequency series object
    """

    generator._generation_domain = 'freq'
    generator._update_domains()

    hlm = generator.generate_fd_modes(**parameter_dict)
    hlm_out = to_gwpy_dict(hlm)
    return gw.GravitationalWaveModes(hlm_out)





########################################################################################################################


########################################################################################################################
