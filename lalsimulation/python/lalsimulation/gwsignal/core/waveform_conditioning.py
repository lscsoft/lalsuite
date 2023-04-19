import lal
import lalsimulation as lalsim
import numpy as np
from astropy import units as u
from . import waveform as wave
from . import conditioning_subroutines as cond
import warnings
import math

# This is basically the same conditioning as done in LALSimulation.
# Specifically, what is done in "LALSimInspiralGeneratorConditioning.c".

def fix_ref_frequency(parameter_dict, generator):
    """
    Function to fix the reference frequency to f22_start or generator specific reference frequency for conditioning routines.

    Parameters
    ----------
    parameter_dict : `dictionary`
        Dictionary of waveform parameters
    generator : `GravitationalWaveGenerator`
        Generator object of the GravitationalWaveGenerator class
    """
    if generator.metadata['implementation']=='LALSimulation':
        return 0
    else:
        cases_by_approximant = {
        'NRSurr' : parameter_dict['f22_start']
        }
    return cases_by_approximant[generator.metadata['approximant']]

########################################################################################

def generate_conditioned_td_waveform_from_td(parameter_dict, generator):
    """
    Function to generate conditioned time-domain waveform from time domain waveform model.

    Parameters
    ----------
    parameter_dict : `dictionary`
            Dictionary of intrinsic / extrinsic gravitational wave parameters
    generator : `GravitationalWaveGenerator`
        Generator object of the GravitationalWaveGenerator class

    Returns
    -------
    hp, hc : GWpy Time series
            Conditioned time domain polarizations
    """

    extra_time_fraction = 0.1   # fraction of waveform duration to add as extra time for tapering
    extra_cycles = 3.0          # more extra time measured in cycles at the starting frequency

    # Get some of the required parameters and get their values (so as not to have issues with units and LAL)
    f_min = parameter_dict['f22_start'].value
    f_ref = parameter_dict['f22_ref'].value
    s1z   = parameter_dict['spin1z'].value
    s2z   = parameter_dict['spin2z'].value

    m1    = parameter_dict['mass1'].si.value
    m2    = parameter_dict['mass2'].si.value

    # Fix reference frequency to fmin/ case by case for Python generators. It is set to zero for lalsim generators as LALSim has its
    # own checks.
    if np.isclose(f_ref, 0):
        f_ref = fix_ref_frequency(parameter_dict, generator)


    # If loweset freq. is below fisco, set fmin to fisco.
    fisco = 1.0 / (np.power(9.0, 1.5) * np.pi * (m1 + m2) * lal.MTSUN_SI / lal.MSUN_SI)
    if (f_min > fisco):
        f_min = fisco


    # Upper chrip time bound
    tchirp = lalsim.SimInspiralChirpTimeBound(f_min, m1, m2, s1z, s2z)

    # Upper bound on the final black hole spin
    s = lalsim.SimInspiralFinalBlackHoleSpinBound(s1z, s2z)

    # Upper bound on the final plunge, merger, and ringdown time
    tmerge = lalsim.SimInspiralMergeTimeBound(m1, m2) + lalsim.SimInspiralRingdownTimeBound(m1 + m2, s)

    # extra time to include for all waveforms to take care of situations
    # where the frequency is close to merger (and is sweeping rapidly):
    # this is a few cycles at the low frequency
    textra = extra_cycles / f_min

    # For conditioning, start waveform at a lower frequency than f_min and then apply tapers between new low freq and f_min.
    fstart = lalsim.SimInspiralChirpStartFrequencyBound((1.0 + extra_time_fraction) * tchirp + tmerge + textra, m1, m2)

    # generate the waveform in the time domain starting at fstart. Add astropy units
    new_parameters = parameter_dict.copy()
    new_parameters['f22_ref'] = f_ref*parameter_dict['f22_start'].unit
    new_parameters['f22_start'] = fstart*parameter_dict['f22_start'].unit


    # Generate the new waveform
    new_parameters['condition']=0
    hp, hc = wave.GenerateTDWaveform(new_parameters, generator)

    times = hp.times
    dt = hp.dt.value
    # Condition the time domain waveform by tapering in the extra time at the beginning
    # And perform the high-pass filtering
    hp, hc = cond.time_array_condition_stage1(hp, hc, dt, extra_time_fraction * tchirp + textra, parameter_dict['f22_start'].value)

    fisco = 1.0 / (np.power(6.0, 1.5) * np.pi * (m1 + m2) * lal.MTSUN_SI / lal.MSUN_SI)
    hp, hc = cond.time_array_condition_stage2(hp, hc, dt, f_min, fisco)

    return hp, hc

########################################################################################

def generate_conditioned_td_waveform_from_fd(parameter_dict, generator):
    """
    Function to generate conditioned time-domain waveform from frequency domain waveform model.

    Parameters
    ----------
    parameter_dict : dictionary
            Dictionary of intrinsic / extrinsic gravitational wave parameters
    generator : `GravitationalWaveGenerator`
        Generator object of the GravitationalWaveGenerator class

    Returns
    -------
    hp, hc : GWpy Time series
            Conditioned time domain polarizations
    """

    extra_time_fraction = 0.1   # fraction of waveform duration to add as extra time for tapering
    extra_cycles = 3.0          # more extra time measured in cycles at the starting frequency

    dt = parameter_dict['deltaT'].value
    original_f_min = f_min = parameter_dict['f22_start'].value
    f_ref = parameter_dict['f22_ref'].value
    f_max = 0.5/dt
    s1z   = parameter_dict['spin1z'].value
    s2z   = parameter_dict['spin2z'].value

    m1    = parameter_dict['mass1'].si.value
    m2    = parameter_dict['mass2'].si.value

    # Fix reference frequency to fmin/ case by case for Python generators. It is set to zero for lalsim generators as LALSim has its
    # own checks.
    if np.isclose(f_ref, 0):
        f_ref = fix_ref_frequency(parameter_dict, generator)

    # If loweset freq. is below fisco, set fmin to fisco.
    fisco = 1.0 / (np.power(9.0, 1.5) * np.pi * (m1 + m2) * lal.MTSUN_SI / lal.MSUN_SI)
    if (f_min > fisco):
        f_min = fisco

    # Upper chrip time bound
    tchirp = lalsim.SimInspiralChirpTimeBound(f_min, m1, m2, s1z, s2z)

    # Upper bound on the final black hole spin
    s = lalsim.SimInspiralFinalBlackHoleSpinBound(s1z, s2z)

    # Upper bound on the final plunge, merger, and ringdown time
    tmerge = lalsim.SimInspiralMergeTimeBound(m1, m2) + lalsim.SimInspiralRingdownTimeBound(m1 + m2, s)

    # extra time to include for all waveforms to take care of situations
    # where the frequency is close to merger (and is sweeping rapidly):
    # this is a few cycles at the low frequency
    textra = extra_cycles / f_min

    # Generate frequency domain conditioned waveform with deltaF = 0 to get
    # enough frequency resolution.

    new_parameters = parameter_dict.copy()
    new_parameters['f22_start'] = f_min*u.Hz
    new_parameters['f22_ref'] = f_ref*u.Hz
    new_parameters['f_max'] = f_max*u.Hz
    new_parameters['deltaF'] = 0.*u.Hz

    # Generate the new waveform
    new_parameters['condition']=1
    hp, hc = wave.GenerateFDWaveform(new_parameters, generator)


    # We want to make sure that this waveform will give something
    # sensible if it is later transformed into the time domain:
    # to avoid the end of the waveform wrapping around to the beginning,
    # we shift waveform backwards in time and compensate for this
    # shift by adjusting the epoch -- note that the conditioned
    # generate_fd_waveform method guarantees that there is the
    # extra padding to do this

    tshift = np.round(textra / dt) * dt # integer number of samples

    kvals = np.arange(0, len(hp))
    phase_shift = np.exp(2*np.pi*1j*hp.df.value*tshift*kvals)
    hp *= phase_shift
    hc *= phase_shift

    if hp.epoch==None or hc.epoch==None:
        hp.epoch=tshift*u.s
        hc.epoch=tshift*u.s
    else:
        hp.epoch+=tshift*u.s
        hc.epoch+=tshift*u.s

    # transform the waveform into the time domain
    hpt = hp.ifft()*2*hp.df
    hct = hc.ifft()*2*hp.df

    # High pass timeseries
    hpt = cond.high_pass_time_series(hpt, dt, original_f_min, 0.99, 8.)
    hct = cond.high_pass_time_series(hct, dt, original_f_min, 0.99, 8.)

    # Resize time series
    fstart = lalsim.SimInspiralChirpStartFrequencyBound((1.0 + extra_time_fraction) * tchirp, m1, m2);
    tchirp = lalsim.SimInspiralChirpTimeBound(fstart, m1, m2, s1z, s2z);

    #total expected chirp length includes merger

    chirplen = round((tchirp + tmerge) /dt)
    end = len(hpt) - np.round(tshift/dt)

    # snip off extra time at beginning and at the end
    hpt = cond.resize_gwpy_timeseries(hpt, end - chirplen, chirplen)
    hct = cond.resize_gwpy_timeseries(hct, end - chirplen, chirplen)

    # final tapering at the beginning and at the end to remove filter transients
    fisco = 1.0 / (np.power(6.0, 1.5) * np.pi * (m1 + m2) * lal.MTSUN_SI / lal.MSUN_SI)
    hpt, hct = cond.time_array_condition_stage2(hpt, hct, dt, f_min, fisco)

    return hpt, hct


########################################################################################

def check_pow_of_2(n):
    """
    Return the next highest power of 2 given number n.
    For eg: if n = 7; exponent will be 3 (2**3=8)
    """
    nv = int(n)
    out=False
    if (nv & (nv-1) == 0) and n!=0:
        out = True
    mant, exponent = math.frexp(n)

    return out, exponent

def generate_conditioned_fd_waveform_from_fd(parameter_dict, generator):
    """
    Function to generate conditioned frequency-domain waveform from frequency domain waveform model.

    Parameters
    ----------
    parameter_dict : dictionary
        Dictionary of intrinsic / extrinsic gravitational wave parameters
    generator : `GravitationalWaveGenerator`
        Generator object of the GravitationalWaveGenerator class

    Returns
    -------
        hp, hc : GWpy Time series
            Conditioned frequency domain polarizations
    """

    extra_time_fraction = 0.1   # fraction of waveform duration to add as extra time for tapering
    extra_cycles = 3.0          # more extra time measured in cycles at the starting frequency

    df    = parameter_dict['deltaF'].value
    f_min = parameter_dict['f22_start'].value
    f_ref = parameter_dict['f22_ref'].value
    f_max = parameter_dict['f_max'].value
    s1z   = parameter_dict['spin1z'].value
    s2z   = parameter_dict['spin2z'].value

    m1    = parameter_dict['mass1'].si.value
    m2    = parameter_dict['mass2'].si.value

    # Fix reference frequency to fmin/ case by case for Python generators. It is set to zero for lalsim generators as LALSim has its
    # own checks.
    if np.isclose(f_ref, 0):
        f_ref = fix_ref_frequency(parameter_dict, generator)

    # Apply condition that f_max rounds to the next power-of-two multiple of deltaF.
    # Round f_max / deltaF to next power of two.
    # Set f_max to the new Nyquist frequency.
    # The length of the chirp signal is then 2 * f_nyquist / deltaF.
    # The time spacing is 1 / (2 * f_nyquist)
    f_nyquist = f_max


    # Check if n is power of 2
    if df!=0:
        n = np.round(f_max/df)
        truth, exponent = check_pow_of_2(n)
        if not truth:
            f_nyquist = 2**(exponent)*df
    deltaT = 0.5/f_nyquist

    # If loweset freq. is below fisco, set fmin to fisco.
    fisco = 1.0 / (np.power(9.0, 1.5) * np.pi * (m1 + m2) * lal.MTSUN_SI / lal.MSUN_SI)
    if (f_min > fisco):
        f_min = fisco

    # Upper chirp-time bound
    tchirp = lalsim.SimInspiralChirpTimeBound(f_min, m1, m2, s1z, s2z)

    # Upper bound on the final black hole spin
    s = lalsim.SimInspiralFinalBlackHoleSpinBound(s1z, s2z)

    # Upper bound on the final plunge, merger, and ringdown time
    tmerge = lalsim.SimInspiralMergeTimeBound(m1, m2) + lalsim.SimInspiralRingdownTimeBound(m1 + m2, s)

    # To condition waveform at lower freq, add some extra early part by changing f_min.
    # New f_min determined by a fixed fraction of the chirp time.
    textra = extra_cycles / f_min
    fstart = lalsim.SimInspiralChirpStartFrequencyBound((1.0 + extra_time_fraction) * tchirp, m1, m2)

    # Revise chirp estimate from fstart
    tchirp = lalsim.SimInspiralChirpTimeBound(fstart, m1, m2, s1z, s2z)

    # Get length required for full waveform with padding upto power of 2
    chirplen = round((tchirp + tmerge + 2.0 * textra) / deltaT);
    truth, exponent = check_pow_of_2(chirplen)
    chirplen = 2**exponent

    if df==0:
        df = 1./(chirplen*deltaT)
    elif df>1./(chirplen*deltaT):
        warnings.warn("Specified frequency interval of %.2f Hz is too large for a chirp of duration %.2f s"%(df, chirplen*deltaT))

    # generate the waveform in the time domain starting at fstart. Add astropy units
    new_parameters = parameter_dict.copy()
    new_parameters['f22_ref'] = f_ref*u.Hz
    new_parameters['f22_start'] = fstart*u.Hz
    new_parameters['deltaF'] = df*u.Hz

    # Generate the new waveform
    new_parameters['condition']=0
    hp, hc = wave.GenerateFDWaveform(new_parameters, generator)

    # Tapering the FD data.
    # Get indices at fstart and fmin
    k0 = int(np.round(fstart/df))
    k1 = int(np.round(f_min/df))

    # Below fstart, ensure all elements are 0
    hp[:k0]=0.
    hc[:k0]=0.

    # From fstart to fmin, apply the taper
    kvals = np.arange(k0, k1)
    w = 0.5 - 0.5*np.cos(np.pi*(kvals-k0)/(k1-k0))
    hp[k0:k1] = w*hp[k0:k1]
    hc[k0:k1] = w*hc[k0:k1]

    # Ensure nyquist frequency elemnt is zero
    hp[len(hp)-1]=0.
    hc[len(hc)-1]=0.

    # Time shift the waveform to avoid waveform wrapping back over itself when converted to time-domain.
    tshift = np.round(tmerge / deltaT) * deltaT
    kvals = np.arange(0, len(hp))
    phase_shift = np.exp(1j*2*np.pi*df*tshift*kvals)

    hp *= phase_shift
    hc *= phase_shift

    # Add tshift to epoch
    hp.epoch = tshift*u.s
    hc.epoch = tshift*u.s

    return hp, hc

########################################################################################

def generate_conditioned_fd_waveform_from_td(parameter_dict, generator):
    """
    Function to generate conditioned frequency-domain waveform from time domain waveform model.

    Parameters
    ----------
    parameter_dict : dictionary
            Dictionary of intrinsic / extrinsic gravitational wave parameters
    generator : `GravitationalWaveGenerator`
        Generator object of the GravitationalWaveGenerator class

    Returns
    -------
    hp, hc : GWpy Time series
            Conditioned frequency domain polarizations
    """

    extra_time_fraction = 0.1   # fraction of waveform duration to add as extra time for tapering
    extra_cycles = 3.0          # more extra time measured in cycles at the starting frequency

    df    = parameter_dict['deltaF'].value
    f_min = parameter_dict['f22_start'].value
    f_ref = parameter_dict['f22_ref'].value
    f_max = parameter_dict['f_max'].value

    s1z   = parameter_dict['spin1z'].value
    s2z   = parameter_dict['spin2z'].value

    m1    = parameter_dict['mass1'].si.value
    m2    = parameter_dict['mass2'].si.value

    # Fix reference frequency to fmin/ case by case for Python generators. It is set to zero for lalsim generators as LALSim has its
    # own checks.
    if np.isclose(f_ref, 0):
        f_ref = fix_ref_frequency(parameter_dict, generator)

    # Apply condition that f_max rounds to the next power-of-two multiple of deltaF.
    # Round f_max / deltaF to next power of two.
    # Set f_max to the new Nyquist frequency.
    # The length of the chirp signal is then 2 * f_nyquist / deltaF.
    # The time spacing is 1 / (2 * f_nyquist)
    f_nyquist = f_max

    # Check if n is power of 2
    if df!=0:
        n = np.round(f_max/df)
        truth, exponent = check_pow_of_2(n)
        if not truth:
            f_nyquist = 2**(exponent)*df

    deltaT = 0.5/f_nyquist

    # generate the waveform in the time domain starting at fstart. Add astropy units
    new_parameters = parameter_dict.copy()
    new_parameters['f22_ref'] = f_ref*u.Hz
    new_parameters['deltaT'] = deltaT*u.s

    # Generate the new waveform
    new_parameters['condition']=1
    hp, hc = wave.GenerateTDWaveform(new_parameters, generator)

    if df==0:
        chirplen = len(hp)
        tt, chirplen_exp = check_pow_of_2(chirplen)
        chirplen = 2**(chirplen_exp)
        df = 1./(chirplen*hp.dt)
    else:
        chirplen=2*f_nyquist/df


    hp = cond.resize_gwpy_timeseries(hp, len(hp)-chirplen,chirplen)
    hc = cond.resize_gwpy_timeseries(hc,len(hc)-chirplen,chirplen)

    hpf = hp.fft()
    hcf = hc.fft()

    # Normalize to match lalsuite
    hpf = hpf/(2*hpf.df)
    hcf = hcf/(2*hpf.df)
    return hpf, hcf
