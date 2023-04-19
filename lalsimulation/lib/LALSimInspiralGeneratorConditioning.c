#include <math.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALDict.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimInspiralWaveformParams.h>
#include "fix_reference_frequency_macro.h"
#include "LALSimInspiralGenerator_private.h"

/* Helper struct storing generator and approximant */
struct internal_data {
    LALSimInspiralGenerator *generator;
    int approx; /* if this is a known named approximant */
};

/* Free memory */
static int finalize(LALSimInspiralGenerator * myself)
{
    struct internal_data *internal_data = myself->internal_data;
    if (internal_data->generator->finalize)
        internal_data->generator->finalize(internal_data->generator);
    LALFree(internal_data->generator);
    LALFree(internal_data);
    return 0;
}

/* this routine is used when reference frequency is the starting frequency */
static int generate_conditioned_td_waveform_from_td_fallback(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, LALDict *params, LALSimInspiralGenerator *myself)
{
    struct internal_data *internal_data = myself->internal_data;
    LALSimInspiralGenerator *internal_generator = internal_data->generator;
    LALSimInspiralApplyTaper taper = LAL_SIM_INSPIRAL_TAPER_START;

    if (internal_generator->generate_td_waveform(hplus, hcross, params, internal_generator) < 0)
        XLAL_ERROR(XLAL_EFUNC);

    /* taper the waveform */
    if (XLALSimInspiralREAL8WaveTaper((*hplus)->data, taper) == XLAL_FAILURE)
        XLAL_ERROR(XLAL_EFUNC);
    if (XLALSimInspiralREAL8WaveTaper((*hcross)->data, taper) == XLAL_FAILURE)
        XLAL_ERROR(XLAL_EFUNC);

    return 0;
}

/* Perform conditioning of TD waveform so that the Fourier transform afterwards is sane. Copy of code from XLALSimInspiralTDFromTD().
 * The redshift correction has been removed, now it is up to the user to apply the proper corrections depending on the meanining of the masses and distance they use.
 */
static int generate_conditioned_td_waveform_from_td(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, LALDict *params, LALSimInspiralGenerator *myself)
{
    struct internal_data *internal_data = myself->internal_data;
    LALSimInspiralGenerator *internal_generator = internal_data->generator;
    int approx = internal_data->approx;
    LALDict *new_params;
    const double extra_time_fraction = 0.1; /* fraction of waveform duration to add as extra time for tapering */
    const double extra_cycles = 3.0; /* more extra time measured in cycles at the starting frequency */
    double original_f_min; /* f_min might be overwritten below, so keep original value */
    double f_min;
    double f_ref;
    double tchirp, tmerge, textra;
    double fisco, fstart;
    double m1, m2, s1z, s2z, s;
    int retval;

    original_f_min = f_min = XLALSimInspiralWaveformParamsLookupF22Start(params);
    f_ref = XLALSimInspiralWaveformParamsLookupF22Ref(params);
    m1 = XLALSimInspiralWaveformParamsLookupMass1(params);
    m2 = XLALSimInspiralWaveformParamsLookupMass2(params);
    s1z = XLALSimInspiralWaveformParamsLookupSpin1z(params);
    s2z = XLALSimInspiralWaveformParamsLookupSpin2z(params);

    /* adjust the reference frequency for certain precessing approximants:
     * if that approximate interprets f_ref==0 to be f_min, set f_ref=f_min;
     * otherwise do nothing */
    FIX_REFERENCE_FREQUENCY(f_ref, f_min, approx);
    
    /* This option recovers the behaviour of SimInspiralTD in the old interface */
    if (XLALDictLookupINT4Value(params, "condition") == 2)
    {
      /* apply redshift correction to dimensionful source-frame quantities */
      REAL8 z=XLALSimInspiralWaveformParamsLookupRedshift(params);
      if (z != 0.0) {
          m1 *= (1.0 + z);
          m2 *= (1.0 + z);
          REAL8 distance = XLALSimInspiralWaveformParamsLookupDistance(params) * (1.0 + z);  /* change from comoving (transverse) distance to luminosity distance */
          XLALSimInspiralWaveformParamsInsertMass1(params, m1);
          XLALSimInspiralWaveformParamsInsertMass2(params, m2);
          XLALSimInspiralWaveformParamsInsertDistance(params, distance);
      }
      /* set redshift to zero so we don't accidentally apply it again later */
      z = 0.0;
      if (params)
        XLALSimInspiralWaveformParamsInsertRedshift(params,z);
    }

    /* if the requested low frequency is below the lowest Kerr ISCO
     * frequency then change it to that frequency */
    fisco = 1.0 / (pow(9.0, 1.5) * LAL_PI * (m1 + m2) * LAL_MTSUN_SI / LAL_MSUN_SI);
    if (f_min > fisco)
        f_min = fisco;

    /* upper bound on the chirp time starting at f_min */
    tchirp = XLALSimInspiralChirpTimeBound(f_min, m1, m2, s1z, s2z);

    /* upper bound on the final black hole spin */
    s = XLALSimInspiralFinalBlackHoleSpinBound(s1z, s2z);

    /* upper bound on the final plunge, merger, and ringdown time */
    tmerge = XLALSimInspiralMergeTimeBound(m1, m2) + XLALSimInspiralRingdownTimeBound(m1 + m2, s);

    /* extra time to include for all waveforms to take care of situations
     * where the frequency is close to merger (and is sweeping rapidly):
     * this is a few cycles at the low frequency */
    textra = extra_cycles / f_min;

    /* time domain approximant: condition by generating a waveform
     * with a lower starting frequency and apply tapers in the
     * region between that lower frequency and the requested
     * frequency f_min; here compute a new lower frequency */
    fstart = XLALSimInspiralChirpStartFrequencyBound((1.0 + extra_time_fraction) * tchirp + tmerge + textra, m1, m2);

    /* generate the waveform in the time domain starting at fstart */
    new_params = XLALDictDuplicate(params);
    XLALSimInspiralWaveformParamsInsertF22Ref(new_params, f_ref);
    XLALSimInspiralWaveformParamsInsertF22Start(new_params, fstart);
    retval = internal_generator->generate_td_waveform(hplus, hcross, new_params, internal_generator);
    XLALDestroyDict(new_params);
    if (retval < 0)
        XLAL_ERROR(XLAL_EFUNC);

    /* condition the time domain waveform by tapering in the extra time
     * at the beginning and high-pass filtering above original f_min */
    XLALSimInspiralTDConditionStage1(*hplus, *hcross, extra_time_fraction * tchirp + textra, original_f_min);

    /* final tapering at the beginning and at the end to remove filter transients */

    /* waveform should terminate at a frequency >= Schwarzschild ISCO
     * so taper one cycle at this frequency at the end; should not make
     * any difference to IMR waveforms */
    fisco = 1.0 / (pow(6.0, 1.5) * LAL_PI * (m1 + m2) * LAL_MTSUN_SI / LAL_MSUN_SI);
    XLALSimInspiralTDConditionStage2(*hplus, *hcross, f_min, fisco);

    return 0;
}

/* Conditioning of a FD waveform and transform it ot TD. Copy of code from XLALSimInspiralTDFromFD().
 * The redshift correction has been removed, now it is up to the user to apply the proper corrections depending on the meanining of the masses and distance they use.
 */
static int generate_conditioned_td_waveform_from_fd(REAL8TimeSeries **hplus, REAL8TimeSeries **hcross, LALDict *params, LALSimInspiralGenerator *myself)
{
    struct internal_data *internal_data = myself->internal_data;
    int approx = internal_data->approx;
    COMPLEX16FrequencySeries *hptilde = NULL;
    COMPLEX16FrequencySeries *hctilde = NULL;
    LALDict *new_params;
    REAL8FFTPlan *plan;
    size_t chirplen, end, k;
    double tshift;
    const double extra_time_fraction = 0.1; /* fraction of waveform duration to add as extra time for tapering */
    const double extra_cycles = 3.0; /* more extra time measured in cycles at the starting frequency */
    double original_f_min; /* f_min might be overwritten below, so keep original value */
    double f_min, f_max, f_ref;
    double deltaT;
    double tchirp, tmerge, textra;
    double fisco, fstart;
    double m1, m2, s1z, s2z, s;
    int retval;

    deltaT = XLALSimInspiralWaveformParamsLookupDeltaT(params);
    original_f_min = f_min = XLALSimInspiralWaveformParamsLookupF22Start(params);
    f_ref = XLALSimInspiralWaveformParamsLookupF22Ref(params);
    f_max = 0.5 / deltaT;
    m1 = XLALSimInspiralWaveformParamsLookupMass1(params);
    m2 = XLALSimInspiralWaveformParamsLookupMass2(params);
    s1z = XLALSimInspiralWaveformParamsLookupSpin1z(params);
    s2z = XLALSimInspiralWaveformParamsLookupSpin2z(params);

    /* adjust the reference frequency for certain precessing approximants:
     * if that approximate interprets f_ref==0 to be f_min, set f_ref=f_min;
     * otherwise do nothing */
    FIX_REFERENCE_FREQUENCY(f_ref, f_min, approx);
    
    /* This option recovers the behaviour of SimInspiralTD in the old interface */
    if (XLALDictLookupINT4Value(params, "condition") == 2)
    {
      /* apply redshift correction to dimensionful source-frame quantities */
      REAL8 z=XLALSimInspiralWaveformParamsLookupRedshift(params);
      if (z != 0.0) {
          m1 *= (1.0 + z);
          m2 *= (1.0 + z);
          REAL8 distance = XLALSimInspiralWaveformParamsLookupDistance(params) * (1.0 + z);  /* change from comoving (transverse) distance to luminosity distance */
          XLALSimInspiralWaveformParamsInsertMass1(params, m1);
          XLALSimInspiralWaveformParamsInsertMass2(params, m2);
          XLALSimInspiralWaveformParamsInsertDistance(params, distance);
      }
      /* set redshift to zero so we don't accidentally apply it again later */
      z = 0.0;
      if (params)
        XLALSimInspiralWaveformParamsInsertRedshift(params,z);
    }

    /* if the requested low frequency is below the lowest Kerr ISCO
     * frequency then change it to that frequency */
    fisco = 1.0 / (pow(9.0, 1.5) * LAL_PI * (m1 + m2) * LAL_MTSUN_SI / LAL_MSUN_SI);
    if (f_min > fisco)
        f_min = fisco;

    /* upper bound on the chirp time starting at f_min */
    tchirp = XLALSimInspiralChirpTimeBound(f_min, m1, m2, s1z, s2z);

    /* upper bound on the final black hole spin */
    s = XLALSimInspiralFinalBlackHoleSpinBound(s1z, s2z);

    /* upper bound on the final plunge, merger, and ringdown time */
    tmerge = XLALSimInspiralMergeTimeBound(m1, m2) + XLALSimInspiralRingdownTimeBound(m1 + m2, s);

    /* extra time to include for all waveforms to take care of situations
     * where the frequency is close to merger (and is sweeping rapidly):
     * this is a few cycles at the low frequency */
    textra = extra_cycles / f_min;

    /* generate the conditioned waveform in the frequency domain */
    /* set deltaF = 0 to get a small enough resolution */
    new_params = XLALDictDuplicate(params);
    XLALSimInspiralWaveformParamsInsertF22Ref(new_params, f_ref);
    XLALSimInspiralWaveformParamsInsertFMax(new_params, f_max);
    XLALSimInspiralWaveformParamsInsertDeltaF(new_params, 0.0);
    retval = myself->generate_fd_waveform(&hptilde, &hctilde, new_params, myself);
    XLALDestroyDict(new_params);
    if (retval < 0)
        XLAL_ERROR(XLAL_EFUNC);

    /* we want to make sure that this waveform will give something
     * sensible if it is later transformed into the time domain:
     * to avoid the end of the waveform wrapping around to the beginning,
     * we shift waveform backwards in time and compensate for this
     * shift by adjusting the epoch -- note that the conditioned
     * generate_fd_waveform method guarantees that there is the
     * extra padding to do this */
    tshift = round(textra / deltaT) * deltaT; /* integer number of samples */
    for (k = 0; k < hptilde->data->length; ++k) {
        double complex phasefac = cexp(2.0 * LAL_PI * I * k * hptilde->deltaF * tshift);
        hptilde->data->data[k] *= phasefac;
        hctilde->data->data[k] *= phasefac;
    }
    XLALGPSAdd(&hptilde->epoch, tshift);
    XLALGPSAdd(&hctilde->epoch, tshift);

    /* transform the waveform into the time domain */
    chirplen = 2 * (hptilde->data->length - 1);
    *hplus = XLALCreateREAL8TimeSeries("H_PLUS", &hptilde->epoch, 0.0, deltaT, &lalStrainUnit, chirplen);
    *hcross = XLALCreateREAL8TimeSeries("H_CROSS", &hctilde->epoch, 0.0, deltaT, &lalStrainUnit, chirplen);
    plan = XLALCreateReverseREAL8FFTPlan(chirplen, 0);
    if (!(*hplus) || !(*hcross) || !plan) {
        XLALDestroyCOMPLEX16FrequencySeries(hptilde);
        XLALDestroyCOMPLEX16FrequencySeries(hctilde);
        XLALDestroyREAL8TimeSeries(*hcross);
        XLALDestroyREAL8TimeSeries(*hplus);
        XLALDestroyREAL8FFTPlan(plan);
        XLAL_ERROR(XLAL_EFUNC);
    }
    XLALREAL8FreqTimeFFT(*hplus, hptilde, plan);
    XLALREAL8FreqTimeFFT(*hcross, hctilde, plan);

    /* apply time domain filter at original f_min */
    XLALHighPassREAL8TimeSeries(*hplus, original_f_min, 0.99, 8);
    XLALHighPassREAL8TimeSeries(*hcross, original_f_min, 0.99, 8);

    /* compute how long a chirp we should have */
    /* revised estimate of chirp length from new start frequency */
    fstart = XLALSimInspiralChirpStartFrequencyBound((1.0 + extra_time_fraction) * tchirp, m1, m2);
    tchirp = XLALSimInspiralChirpTimeBound(fstart, m1, m2, s1z, s2z);

    /* total expected chirp length includes merger */
    chirplen = round((tchirp + tmerge) / deltaT);

    /* amount to snip off at the end is tshift */
    end = (*hplus)->data->length - round(tshift / deltaT);

    /* snip off extra time at beginning and at the end */
    XLALResizeREAL8TimeSeries(*hplus, end - chirplen, chirplen);
    XLALResizeREAL8TimeSeries(*hcross, end - chirplen, chirplen);

    /* clean up */
    XLALDestroyREAL8FFTPlan(plan);
    XLALDestroyCOMPLEX16FrequencySeries(hptilde);
    XLALDestroyCOMPLEX16FrequencySeries(hctilde);

    /* final tapering at the beginning and at the end to remove filter transients */

    /* waveform should terminate at a frequency >= Schwarzschild ISCO
     * so taper one cycle at this frequency at the end; should not make
     * any difference to IMR waveforms */
    fisco = 1.0 / (pow(6.0, 1.5) * LAL_PI * (m1 + m2) * LAL_MTSUN_SI / LAL_MSUN_SI);
    XLALSimInspiralTDConditionStage2(*hplus, *hcross, f_min, fisco);

    return 0;
}

/* Conditioning a Fourier domain waveform to be properly transformed to the time domain. This code was taken from the original XLALSimInspiralFD() function, corresponding to the FD approximants part.
 * The redshift correction has been removed, now it is up to the user to apply the proper corrections depending on the meanining of the masses and distance they use.
 */
static int generate_conditioned_fd_waveform_from_fd(COMPLEX16FrequencySeries **hplus, COMPLEX16FrequencySeries **hcross, LALDict *params, LALSimInspiralGenerator *myself)
{
    struct internal_data *internal_data = myself->internal_data;
    LALSimInspiralGenerator *internal_generator = internal_data->generator;
    int approx = internal_data->approx;
    LALDict *new_params;
    const double extra_time_fraction = 0.1; /* fraction of waveform duration to add as extra time for tapering */
    const double extra_cycles = 3.0; /* more extra time measured in cycles at the starting frequency */
    double chirplen, deltaT, deltaF, f_nyquist;
    double f_min, f_max, f_ref;
    double tchirp, tmerge, textra, tshift;
    double fstart, fisco;
    double m1, m2, s1z, s2z, s;
    size_t k, k0, k1;
    int chirplen_exp;
    int retval;
    size_t n;

    deltaF = XLALSimInspiralWaveformParamsLookupDeltaF(params);
    f_min = XLALSimInspiralWaveformParamsLookupF22Start(params);
    f_max = XLALSimInspiralWaveformParamsLookupFMax(params);
    f_ref = XLALSimInspiralWaveformParamsLookupF22Ref(params);
    m1 = XLALSimInspiralWaveformParamsLookupMass1(params);
    m2 = XLALSimInspiralWaveformParamsLookupMass2(params);
    s1z = XLALSimInspiralWaveformParamsLookupSpin1z(params);
    s2z = XLALSimInspiralWaveformParamsLookupSpin2z(params);

    /* adjust the reference frequency for certain precessing approximants:
     * if that approximate interprets f_ref==0 to be f_min, set f_ref=f_min;
     * otherwise do nothing */
    FIX_REFERENCE_FREQUENCY(f_ref, f_min, approx);

    /* This option recovers the behaviour of SimInspiralFD in the old interface */
    if (XLALDictLookupINT4Value(params, "condition") == 2)
    {
      /* apply redshift correction to dimensionful source-frame quantities */
      REAL8 z=XLALSimInspiralWaveformParamsLookupRedshift(params);
      if (z != 0.0) {
          m1 *= (1.0 + z);
          m2 *= (1.0 + z);
          REAL8 distance = XLALSimInspiralWaveformParamsLookupDistance(params) * (1.0 + z);  /* change from comoving (transverse) distance to luminosity distance */
          XLALSimInspiralWaveformParamsInsertMass1(params, m1);
          XLALSimInspiralWaveformParamsInsertMass2(params, m2);
          XLALSimInspiralWaveformParamsInsertDistance(params, distance);
      }
      /* set redshift to zero so we don't accidentally apply it again later */
      z = 0.0;
      if (params)
        XLALSimInspiralWaveformParamsInsertRedshift(params,z);
    }
    
    /* Apply condition that f_max rounds to the next power-of-two multiple
     * of deltaF.
     * Round f_max / deltaF to next power of two.
     * Set f_max to the new Nyquist frequency.
     * The length of the chirp signal is then 2 * f_nyquist / deltaF.
     * The time spacing is 1 / (2 * f_nyquist) */
    f_nyquist = f_max;
    if (deltaF != 0) {
        n = round(f_max / deltaF);
        if ((n & (n - 1))) { /* not a power of 2 */
            frexp(n, &chirplen_exp);
            f_nyquist = ldexp(1.0, chirplen_exp) * deltaF;
            XLAL_PRINT_WARNING("f_max/deltaF = %g/%g = %g is not a power of two: changing f_max to %g", f_max, deltaF, f_max/deltaF, f_nyquist);
        }
    }
    deltaT = 0.5 / f_nyquist;

    /* generate a FD waveform and condition it by applying tapers at
     * frequencies between a frequency below the requested f_min and
     * f_min; also wind the waveform in phase in case it would wrap-
     * around at the merger time */

    /* if the requested low frequency is below the lowest Kerr ISCO
     * frequency then change it to that frequency */
    fisco = 1.0 / (pow(9.0, 1.5) * LAL_PI * (m1 + m2) * LAL_MTSUN_SI / LAL_MSUN_SI);
    if (f_min > fisco)
        f_min = fisco;

    /* upper bound on the chirp time starting at f_min */
    tchirp = XLALSimInspiralChirpTimeBound(f_min, m1, m2, s1z, s2z);

    /* upper bound on the final plunge, merger, and ringdown time */
    switch (approx) {
    case TaylorF2:
    case TaylorF2Ecc:
    case TaylorF2NLTides:
    case SpinTaylorF2:
    case TaylorF2RedSpin:
    case TaylorF2RedSpinTidal:
    case SpinTaylorT4Fourier:
        /* inspiral-only models: no merger time */
        tmerge = 0.0;
        break;
    default:
        /* IMR model: estimate plunge and merger time */
        /* sometimes these waveforms have phases that
         * cause them to wrap-around an amount equal to
         * the merger-ringodwn time, so we will undo
         * that here */
        s = XLALSimInspiralFinalBlackHoleSpinBound(s1z, s2z);
        tmerge = XLALSimInspiralMergeTimeBound(m1, m2) + XLALSimInspiralRingdownTimeBound(m1 + m2, s);
        break;
    }

    /* new lower frequency to start the waveform: add some extra early
     * part over which tapers may be applied, the extra amount being
     * a fixed fraction of the chirp time; add some additional padding
     * equal to a few extra cycles at the low frequency as well for
     * safety and for other routines to use */
    textra = extra_cycles / f_min;
    fstart = XLALSimInspiralChirpStartFrequencyBound((1.0 + extra_time_fraction) * tchirp, m1, m2);

    /* revise (over-)estimate of chirp from new start frequency */
    tchirp = XLALSimInspiralChirpTimeBound(fstart, m1, m2, s1z, s2z);

    /* need a long enough segment to hold a whole chirp with some padding */
    /* length of the chirp in samples */
    chirplen = round((tchirp + tmerge + 2.0 * textra) / deltaT);
    /* make chirplen next power of two */
    frexp(chirplen, &chirplen_exp);
    chirplen = ldexp(1.0, chirplen_exp);
    /* frequency resolution */
    if (deltaF == 0.0)
        deltaF = 1.0 / (chirplen * deltaT);
    else if (deltaF > 1.0 / (chirplen * deltaT))
        XLAL_PRINT_WARNING("Specified frequency interval of %g Hz is too large for a chirp of duration %g s", deltaF, chirplen * deltaT);

    /* generate the waveform in the frequency domain starting at fstart */
    new_params = XLALDictDuplicate(params);
    XLALSimInspiralWaveformParamsInsertF22Ref(new_params, f_ref);
    XLALSimInspiralWaveformParamsInsertF22Start(new_params, fstart);
    XLALSimInspiralWaveformParamsInsertDeltaF(new_params, deltaF);
    retval = internal_generator->generate_fd_waveform(hplus, hcross, new_params, internal_generator);
    XLALDestroyDict(new_params);
    if (retval < 0)
        XLAL_ERROR(XLAL_EFUNC);

    /* taper frequencies between fstart and f_min */
    k0 = round(fstart / (*hplus)->deltaF);
    k1 = round(f_min / (*hplus)->deltaF);
    /* make sure it is zero below fstart */
    for (k = 0; k < k0; ++k) {
        (*hplus)->data->data[k] = 0.0;
        (*hcross)->data->data[k] = 0.0;
    }
    /* taper between fstart and f_min */
    for ( ; k < k1; ++k) {
        double w = 0.5 - 0.5 * cos(LAL_PI * (k - k0) / (double)(k1 - k0));
        (*hplus)->data->data[k] *= w;
        (*hcross)->data->data[k] *= w;
    }
    /* make sure Nyquist frequency is zero */
    (*hplus)->data->data[(*hplus)->data->length - 1] = 0.0;
    (*hcross)->data->data[(*hcross)->data->length - 1] = 0.0;

    /* we want to make sure that this waveform will give something
     * sensible if it is later transformed into the time domain:
     * to avoid the end of the waveform wrapping around to the beginning,
     * we shift waveform backwards in time and compensate for this
     * shift by adjusting the epoch */
    tshift = round(tmerge / deltaT) * deltaT; /* integer number of time samples */
    for (k = 0; k < (*hplus)->data->length; ++k) {
        double complex phasefac = cexp(2.0 * LAL_PI * I * k * deltaF * tshift);
        (*hplus)->data->data[k] *= phasefac;
        (*hcross)->data->data[k] *= phasefac;
    }
    XLALGPSAdd(&(*hplus)->epoch, tshift);
    XLALGPSAdd(&(*hcross)->epoch, tshift);

    return 0;
}

/* Transform a Fourier domain waveform to the time domain. This code was taken from the original XLALSimInspiralFD() function, corresponding to the TD approximants part.
 * The redshift correction has been removed, now it is up to the user to apply the proper corrections depending on the meanining of the masses and distance they use.
 */
static int generate_conditioned_fd_waveform_from_td(COMPLEX16FrequencySeries **hplus, COMPLEX16FrequencySeries **hcross, LALDict *params, LALSimInspiralGenerator *myself)
{
    struct internal_data *internal_data = myself->internal_data;
    int approx = internal_data->approx;
    REAL8TimeSeries *hp = NULL;
    REAL8TimeSeries *hc = NULL;
    LALDict *new_params;
    REAL8FFTPlan *plan;
    double chirplen, deltaT, deltaF, f_nyquist;
    double f_min, f_max, f_ref;
    int chirplen_exp;
    int retval;
    size_t n;

    deltaF = XLALSimInspiralWaveformParamsLookupDeltaF(params);
    f_min = XLALSimInspiralWaveformParamsLookupF22Start(params);
    f_max = XLALSimInspiralWaveformParamsLookupFMax(params);
    f_ref = XLALSimInspiralWaveformParamsLookupF22Ref(params);

    /* adjust the reference frequency for certain precessing approximants:
     * if that approximate interprets f_ref==0 to be f_min, set f_ref=f_min;
     * otherwise do nothing */
    FIX_REFERENCE_FREQUENCY(f_ref, f_min, approx);
    
    /* This option recovers the behaviour of SimInspiralFD in the old interface */
    if (XLALDictLookupINT4Value(params, "condition") == 2)
    {
      /* apply redshift correction to dimensionful source-frame quantities */
      REAL8 z=XLALSimInspiralWaveformParamsLookupRedshift(params);
      if (z != 0.0) {
          REAL8 m1 = XLALSimInspiralWaveformParamsLookupMass1(params) * (1.0 + z);
          REAL8 m2 = XLALSimInspiralWaveformParamsLookupMass2(params) * (1.0 + z);
          REAL8 distance = XLALSimInspiralWaveformParamsLookupDistance(params) * (1.0 + z);  /* change from comoving (transverse) distance to luminosity distance */
          XLALSimInspiralWaveformParamsInsertMass1(params, m1);
          XLALSimInspiralWaveformParamsInsertMass2(params, m2);
          XLALSimInspiralWaveformParamsInsertDistance(params, distance);
      }
      /* set redshift to zero so we don't accidentally apply it again later */
      z = 0.0;
      if (params)
        XLALSimInspiralWaveformParamsInsertRedshift(params,z);
    }
    

    /* Apply condition that f_max rounds to the next power-of-two multiple
     * of deltaF.
     * Round f_max / deltaF to next power of two.
     * Set f_max to the new Nyquist frequency.
     * The length of the chirp signal is then 2 * f_nyquist / deltaF.
     * The time spacing is 1 / (2 * f_nyquist) */
    f_nyquist = f_max;
    if (deltaF != 0) {
        n = round(f_max / deltaF);
        if ((n & (n - 1))) { /* not a power of 2 */
            frexp(n, &chirplen_exp);
            f_nyquist = ldexp(1.0, chirplen_exp) * deltaF;
            XLAL_PRINT_WARNING("f_max/deltaF = %g/%g = %g is not a power of two: changing f_max to %g", f_max, deltaF, f_max/deltaF, f_nyquist);
        }
    }
    deltaT = 0.5 / f_nyquist;

    /* generate conditioned waveform in time domain */
    new_params = XLALDictDuplicate(params);
    XLALSimInspiralWaveformParamsInsertF22Ref(new_params, f_ref);
    XLALSimInspiralWaveformParamsInsertDeltaT(new_params, deltaT);
    retval = myself->generate_td_waveform(&hp, &hc, new_params, myself);
    XLALDestroyDict(new_params);
    if (retval < 0)
        XLAL_ERROR(XLAL_EFUNC);

    /* frequency resolution */
    if (deltaF == 0.0) {
        /* round length of time domain signal to next power of two */
        chirplen = hp->data->length;
        frexp(chirplen, &chirplen_exp);
        chirplen = ldexp(1.0, chirplen_exp);
        deltaF = 1.0 / (chirplen * hp->deltaT);
    } else {
        /* set chirp length using precomputed Nyquist */
        chirplen = 2 * f_nyquist / deltaF;
        if (chirplen < hp->data->length)
            XLAL_PRINT_WARNING("Specified frequency interval of %g Hz is too large for a chirp of duration %g s with Nyquist frequency %g Hz. The inspiral will be truncated.", deltaF, hp->data->length * deltaT, f_nyquist);
    }

    /* resize waveforms to the required length */
    XLALResizeREAL8TimeSeries(hp, hp->data->length - (size_t) chirplen, (size_t) chirplen);
    XLALResizeREAL8TimeSeries(hc, hc->data->length - (size_t) chirplen, (size_t) chirplen);

    /* put the waveform in the frequency domain */
    /* (the units will correct themselves) */
    *hplus = XLALCreateCOMPLEX16FrequencySeries("FD H_PLUS", &hp->epoch, 0.0, deltaF, &lalDimensionlessUnit, (size_t) chirplen / 2 + 1);
    *hcross = XLALCreateCOMPLEX16FrequencySeries("FD H_CROSS", &hc->epoch, 0.0, deltaF, &lalDimensionlessUnit, (size_t) chirplen / 2 + 1);
    plan = XLALCreateForwardREAL8FFTPlan((size_t) chirplen, 0);
    XLALREAL8TimeFreqFFT(*hcross, hc, plan);
    XLALREAL8TimeFreqFFT(*hplus, hp, plan);

    /* clean up */
    XLALDestroyREAL8FFTPlan(plan);
    XLALDestroyREAL8TimeSeries(hc);
    XLALDestroyREAL8TimeSeries(hp);

    return 0;
}

//static int generate_conditioned_td_modes(SphHarmTimeSeries **hlm, LALDict *params, LALSimInspiralGenerator *myself);

//static int generate_conditioned_fd_modes(SphHarmFrequencySeries **hlm, LALDict *params, LALSimInspiralGenerator *myself);

/* Function to assign the proper conditioning generator method for a given approximant. */
int XLALSimInspiralGeneratorAddConditioningForApproximant(LALSimInspiralGenerator *generator, int approximant)
{
    struct internal_data *internal_data;

    internal_data = LALMalloc(sizeof(*internal_data));
    internal_data->approx = approximant;
    internal_data->generator = LALMalloc(sizeof(*internal_data->generator));
    memcpy(internal_data->generator, generator, sizeof(*internal_data->generator));

    generator->internal_data = internal_data;
    generator->finalize = finalize;

    if (internal_data->generator->generate_td_waveform) {
        if (internal_data->approx == -1) {
            /* Not a recognized approximant:
             * assume we can use regular conditioning */
            generator->generate_td_waveform = generate_conditioned_td_waveform_from_td;
        } else {
          /* If using approximants for which reference frequency is the starting frequency
            * generate using XLALSimInspiralChooseTDWaveform and apply the
            * LAL Taper 'LAL_SIM_INSPIRAL_TAPER_START' instead of
            * XLALSimInspiralTDConditionStage1 and XLALSimInspiralTDConditionStage2
            * as is done in XLALSimInspiralTDFromTD.
            * This is because XLALSimInspiralTDFromTD modifies the start frequency
            * which is not always possible with NR_hdf5 waveforms.
            * Do the same (ChooseTDWaveform+LALTaper) if using approximants for
            * which a starting frequency of zero is allowed, as determined from
            * XLALSimInspiralGetAllowZeroMinFreqFromApproximant. This is because
            * XLALSimInspiralTDFromTD does not properly handle f_min=0. For models
            * that allow f_min=0, this (ChooseTDWaveform+LALTaper) is the behaviour
            * independent of what f_min is passed.
            */

            // Check whether for the given approximant reference frequency is the starting frequency
            SpinFreq spin_freq_flag = XLALSimInspiralGetSpinFreqFromApproximant(approximant);
            // Check whether for the given approximant, f_min=0 is allowed.
            AllowZeroMinFreq allow_zero_fmin_flag = XLALSimInspiralGetAllowZeroMinFreqFromApproximant(approximant);

            if (spin_freq_flag == LAL_SIM_INSPIRAL_SPINS_CASEBYCASE || spin_freq_flag == LAL_SIM_INSPIRAL_SPINS_FLOW || allow_zero_fmin_flag == LAL_SIM_INSPIRAL_ALLOW_ZERO_FMIN)
                generator->generate_td_waveform = generate_conditioned_td_waveform_from_td_fallback;
            else
                generator->generate_td_waveform = generate_conditioned_td_waveform_from_td;
       }
    } else if (internal_data->generator->generate_fd_waveform)
        generator->generate_td_waveform = generate_conditioned_td_waveform_from_fd;

    if (internal_data->generator->generate_fd_waveform)
        generator->generate_fd_waveform = generate_conditioned_fd_waveform_from_fd;
    else if (internal_data->generator->generate_td_waveform)
        generator->generate_fd_waveform = generate_conditioned_fd_waveform_from_td;

    /* FUTURE: implement routines for conditioning modes */
    // generator->generate_td_modes = generate_conditioned_td_modes;
    // generator->generate_fd_modes = generate_conditioned_fd_modes;


    return 0;
}

/* Function to assign standard conditioning. Currently this is `generate_conditioned_td_waveform_from_td`. */
int XLALSimInspiralGeneratorAddStandardConditioning(LALSimInspiralGenerator *generator)
{
    return XLALSimInspiralGeneratorAddConditioningForApproximant(generator, -1);
}
