/*
 * Copyright (C) 2008 J. Creighton, T. Creighton
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <math.h>
#include <gsl/gsl_sf_trig.h>
#include <lal/LALSimulation.h>
#include <lal/LALDetectors.h>
#include <lal/DetResponse.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/TimeDelay.h>
#include <lal/SkyCoordinates.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Window.h>
#include "check_series_macros.h"


/*
 * ============================================================================
 *
 *                                 Utilities
 *
 * ============================================================================
 */


/*
 * Note:  this function will round really really large numbers "up" to 0.
 */


static unsigned long round_up_to_power_of_two(unsigned long x)
{
	unsigned n;

	x--;
	/* if x shares bits with x + 1, then x + 1 is not a power of 2 */
	for(n = 1; n && (x & (x + 1)); n *= 2)
		x |= x >> n;

	return x + 1;
}


/**
 * Turn a detector prefix string into a LALDetector structure.  The first
 * two characters of the input string are used as the instrument name,
 * which allows channel names in the form "H1:LSC-STRAIN" to be used.  The
 * return value is a pointer into the lalCachedDetectors array, so
 * modifications to the contents are global.  Make a copy of the structure
 * if you want to modify it safely.
 */


const LALDetector *XLALDetectorPrefixToLALDetector(
	const char *string
)
{
	int i;

	for(i = 0; i < LAL_NUM_DETECTORS; i++)
		if(!strncmp(string, lalCachedDetectors[i].frDetector.prefix, 2))
			return &lalCachedDetectors[i];

	XLALPrintError("%s(): error: can't identify instrument from string \"%s\"\n", __func__, string);
	XLAL_ERROR_NULL(XLAL_EDATA);
}


/*
 * ============================================================================
 *
 *                            Injection Machinery
 *
 * ============================================================================
 */


#if 0
REAL8TimeSeries * XLALSimQuasiPeriodicInjectionREAL8TimeSeries( REAL8TimeSeries *aplus, REAL8TimeSeries *across, REAL8TimeSeries *frequency, REAL8TimeSeries *phi, REAL8TimeSeries *psi, SkyPosition *position, LALDetector *detector, EphemerisData *ephemerides, LIGOTimeGPS *start, REAL8 deltaT, UINT4 length, COMPLEX16FrequencySeries *response )
{
	REAL8 slowDeltaT = 600.0; /* 10 minutes */
	UINT4 slowLength;

	/* sanity checks */

	/* SET UP LOOK-UP TABLES */

	/* construct time-dependent time-delay look-up table */
	/* uses ephemeris */
	/* sampled slowDeltaT */
	/* starts 8 min before injection start, ends 8 min after injection end */
	slowLength = floor(0.5 + (length*deltaT + 16.0)/slowDeltaT);

	/* construct time-dependent polarization response look-up table */
	/* sampled at slowDeltaT */
	/* starts 8 min before injection start, ends 8 min after injection end */
	fplus;
	fcross;



	/* apply time-dependent time-delay */
	/* apply time-dependent polarization response */
	/* apply frequency-dependent detector response */
	for ( j = 0; j < injection->length; ++j ) {
		REAL8 dt;
		REAL8 fpl;
		REAL8 fcr;
		REAL8 apl;
		REAL8 acr;
		REAL8 freq;
		REAL8 phase;
		REAL8 pol;

		/* interpolate to get dt, fpl, fcr */

		/* interpolate to get apl, acr, freq, phase, pol */
		/* but not at detector time: detector time plus dt */

		/* rotate fpl and fcr using pol */

		/* adjust apl, acr, and phase by detector response at freq */

		injection->data->data[j] = fpl*apl*cos(phase) + fcr*acr*sin(phase); /* or something like this */

	}



}
#endif


/**
 * Input
 *
 * - h+ and hx time series for the injection with their epochs set to the
 * start of those time series at the geocentre (for simplicity the epochs
 * must be the same),
 *
 * - the right ascension and declination of the source in radians.
 *
 * - the orientation of the wave co-ordinate system, psi, in radians.
 *
 * - the detector into which the injection is destined to be injected,
 *
 * Output
 *
 * The strain time series as seen in the detector, with the epoch set to
 * the start of the time series at that detector.
 *
 * Notes
 *
 * Antenna response factors are computed sample-by-sample, but waveform is
 * not Doppler corrected or otherwise re-interpolated to account for
 * detector motion.
 *
 * The output time series units are the same as the two input time series
 * (which must both have the same sample units).
 */


REAL8TimeSeries *XLALSimDetectorStrainREAL8TimeSeries(
	const REAL8TimeSeries *hplus,
	const REAL8TimeSeries *hcross,
	REAL8 right_ascension,
	REAL8 declination,
	REAL8 psi,
	const LALDetector *detector
)
{
	char *name;
	REAL8TimeSeries *h = NULL;
	unsigned i;

	LAL_CHECK_VALID_SERIES(hplus, NULL);
	LAL_CHECK_VALID_SERIES(hcross, NULL);
	LAL_CHECK_CONSISTENT_TIME_SERIES(hplus, hcross, NULL);

	/* generate name */

	name = XLALMalloc(strlen(detector->frDetector.prefix) + 11);
	if(!name)
		goto error;
	sprintf(name, "%s injection", detector->frDetector.prefix);

	/* allocate output time series. */

	h = XLALCreateREAL8TimeSeries(name, &hplus->epoch, hplus->f0, hplus->deltaT, &hplus->sampleUnits, hplus->data->length);
	XLALFree(name);
	if(!h)
		goto error;

	/* add the detector's geometric delay.  after this, epoch = the
	 * time of the injection time series' first sample at the desired
	 * detector */

	XLALGPSAdd(&h->epoch, XLALTimeDelayFromEarthCenter(detector->location, right_ascension, declination, &h->epoch));

	/* project + and x time series onto detector */

	for(i = 0; i < h->data->length; i++) {
		LIGOTimeGPS t = h->epoch;
		double fplus, fcross;

		/* time of sample in detector */
		XLALGPSAdd(&t, i * h->deltaT);

		/* detector's response at that time */
		XLALComputeDetAMResponse(&fplus, &fcross, (const REAL4(*)[3])detector->response, right_ascension, declination, psi, XLALGreenwichMeanSiderealTime(&t));
		if(XLAL_IS_REAL8_FAIL_NAN(fplus) || XLAL_IS_REAL8_FAIL_NAN(fcross))
			goto error;

		h->data->data[i] = fplus * hplus->data->data[i] + fcross * hcross->data->data[i];
	}

	/* done */

	return h;

error:
	XLALDestroyREAL8TimeSeries(h);
	XLAL_ERROR_NULL(XLAL_EFUNC);
}


/**
 * Essentially a wrapper for XLALAddREAL8TimeSeries(), but performs
 * sub-sample re-interpolation to adjust the source time series epoch to
 * lie on an integer sample boundary in the target time series.  This
 * transformation is done in the frequency domain, so it is convenient to
 * allow a response function to be applied at the same time.  Passing NULL
 * for the response function turns this feature off (i.e., uses a unit
 * response).
 *
 * NOTE:  the source time series is modified in place by this function!
 *
 * This function accepts source and target time series whose units are not
 * the same, and allows the two time series to be herterodyned (although it
 * currently requires them to have the same heterodyne frequency).
 */


int XLALSimAddInjectionREAL8TimeSeries(
	REAL8TimeSeries *target,
	REAL8TimeSeries *h,
	const COMPLEX16FrequencySeries *response
)
{
	/* 1 ns is about 10^-5 samples at 16384 Hz */
	const double noop_threshold = 1e-4;	/* samples */
	/* the source time series is padded with at least this many 0's at
	 * the start and end before re-interpolation in an attempt to
	 * suppress aperiodicity artifacts, and 1/2 this many samples is
	 * clipped from the start and end afterwards */
	const unsigned aperiodicity_suppression_buffer = 32768;
	REAL8Window *window;
	unsigned i;
	double start_sample_int;
	double start_sample_frac;

	/* check input */

	/* FIXME:  since we route the source time series through the
	 * frequency domain, it might not be hard to adjust the heterodyne
	 * frequency and sample rate to match the target time series
	 * instead of making it an error if they don't match. */

	if(h->deltaT != target->deltaT || h->f0 != target->f0) {
		XLALPrintError("%s(): error: input sample rates or heterodyne frequencies do not match\n", __func__);
		XLAL_ERROR(XLAL_EINVAL);
	}

	/* extend the source time series by adding the "aperiodicity
	 * padding" to the start and end.  for efficiency's sake, make sure
	 * the new length is a power of two. */

	i = round_up_to_power_of_two(h->data->length + 2 * aperiodicity_suppression_buffer);
	if(i < h->data->length) {
		/* integer overflow */
		XLALPrintError("%s(): error: source time series too long\n", __func__);
		XLAL_ERROR(XLAL_EBADLEN);
	}
	i -= h->data->length;
	if(!XLALResizeREAL8TimeSeries(h, -(int) (i / 2), h->data->length + i))
		XLAL_ERROR(XLAL_EFUNC);

	/* compute the integer and fractional parts of the sample index in
	 * the target time series on which the source time series begins.
	 * modf() returns integer and fractional parts that have the same
	 * sign, e.g. -3.9 --> -3 + -0.9.  we adjust these so that the
	 * magnitude of the fractional part is not greater than 0.5, e.g.
	 * -3.9 --> -4 + 0.1, so that we never do more than 1/2 a sample of
	 * re-interpolation.  I don't know if really makes any difference,
	 * though */

	start_sample_frac = modf(XLALGPSDiff(&h->epoch, &target->epoch) / target->deltaT, &start_sample_int);
	if(start_sample_frac < -0.5) {
		start_sample_frac += 1.0;
		start_sample_int -= 1.0;
	} else if(start_sample_frac > +0.5) {
		start_sample_frac -= 1.0;
		start_sample_int += 1.0;
	}

	if(fabs(start_sample_frac) > noop_threshold) {
		COMPLEX16FrequencySeries *tilde_h;
		REAL8FFTPlan *plan;

		/* transform source time series to frequency domain.  the FFT
		 * function populates the frequency series' metadata with the
		 * appropriate values. */

		tilde_h = XLALCreateCOMPLEX16FrequencySeries(NULL, &h->epoch, 0, 0, &lalDimensionlessUnit, h->data->length / 2 + 1);
		plan = XLALCreateForwardREAL8FFTPlan(h->data->length, 0);
		if(!tilde_h || !plan) {
			XLALDestroyCOMPLEX16FrequencySeries(tilde_h);
			XLALDestroyREAL8FFTPlan(plan);
			XLAL_ERROR(XLAL_EFUNC);
		}
		i = XLALREAL8TimeFreqFFT(tilde_h, h, plan);
		XLALDestroyREAL8FFTPlan(plan);
		if(i) {
			XLALDestroyCOMPLEX16FrequencySeries(tilde_h);
			XLAL_ERROR(XLAL_EFUNC);
		}

		/* apply sub-sample time correction and optional response function
		 * */

		for(i = 0; i < tilde_h->data->length; i++) {
			const double f = tilde_h->f0 + i * tilde_h->deltaF;
			COMPLEX16 fac;

			/* phase for sub-sample time correction */

			fac = cexp(-I * LAL_TWOPI * f * start_sample_frac * target->deltaT);

			/* divide the source by the response function.  if a
			 * frequency is required that lies outside the domain of
			 * definition of the response function, then the response
			 * is assumed equal to its value at the nearest edge of the
			 * domain of definition.  within the domain of definition,
			 * frequencies are rounded to the nearest bin.  if the
			 * response function is zero in some bin, then the source
			 * data is zeroed in that bin (instead of dividing by 0).
			 * */

			/* FIXME:  should we use GSL to construct an interpolator
			 * for the modulus and phase as functions of frequency, and
			 * use that to evaluate the response?  instead of rounding
			 * to nearest bin? */

			if(response) {
				int j = floor((f - response->f0) / response->deltaF + 0.5);
				if(j < 0)
					j = 0;
				else if((unsigned) j > response->data->length - 1)
					j = response->data->length - 1;
				if(response->data->data[j] == 0.0)
					fac = 0.0;
				else
					fac /= response->data->data[j];
			}

			/* apply factor */

			tilde_h->data->data[i] *= fac;
		}

		/* adjust DC and Nyquist components.  the DC component must always
		 * be real-valued.  because we have adjusted the source time series
		 * to have a length that is an even integer (we've made it a power
		 * of 2) the Nyquist component must also be real valued. */

		if(response) {
			/* a response function has been provided.  zero the DC and
			 * Nyquist components */
			if(tilde_h->f0 == 0.0)
				tilde_h->data->data[0] = 0.0;
			tilde_h->data->data[tilde_h->data->length - 1] = 0.0;
		} else {
			/* no response has been provided.  set the phase of the DC
			 * component to 0, set the imaginary component of the
			 * Nyquist to 0 */
			if(tilde_h->f0 == 0.0)
				tilde_h->data->data[0] = cabs(tilde_h->data->data[0]);
			tilde_h->data->data[tilde_h->data->length - 1] = creal(tilde_h->data->data[tilde_h->data->length - 1]);
		}

		/* return to time domain */

		plan = XLALCreateReverseREAL8FFTPlan(h->data->length, 0);
		if(!plan) {
			XLALDestroyCOMPLEX16FrequencySeries(tilde_h);
			XLAL_ERROR(XLAL_EFUNC);
		}
		i = XLALREAL8FreqTimeFFT(h, tilde_h, plan);
		XLALDestroyREAL8FFTPlan(plan);
		XLALDestroyCOMPLEX16FrequencySeries(tilde_h);
		if(i)
			XLAL_ERROR(XLAL_EFUNC);

		/* the deltaT can get "corrupted" by floating point round-off
		 * during its trip through the frequency domain.  since this
		 * function starts by confirming that the sample rate of the source
		 * matches that of the target time series, we can use the target
		 * series' sample rate to reset the source's sample rate to its
		 * original value.  but we do a check to make sure we're not
		 * masking a real bug */

		if(fabs(h->deltaT - target->deltaT) / target->deltaT > 1e-12) {
			XLALPrintError("%s(): error: oops, internal sample rate mismatch\n", __func__);
			XLAL_ERROR(XLAL_EERR);
		}
		h->deltaT = target->deltaT;
	}

	/* set source epoch from target epoch and integer sample offset */

	h->epoch = target->epoch;
	XLALGPSAdd(&h->epoch, start_sample_int * target->deltaT);

	/* clip half of the "aperiodicity padding" from the start and end
	 * of the source time series in a continuing effort to suppress
	 * aperiodicity artifacts. */

	if(!XLALResizeREAL8TimeSeries(h, aperiodicity_suppression_buffer / 2, h->data->length - aperiodicity_suppression_buffer))
		XLAL_ERROR(XLAL_EFUNC);

	/* apply a Tukey window whose tapers lie within the remaining
	 * aperiodicity padding. leaving one sample of the aperiodicty
	 * padding untouched on each side of the original time series
	 * because the data might have been shifted into it */

	window = XLALCreateTukeyREAL8Window(h->data->length, (double) (aperiodicity_suppression_buffer - 2) / h->data->length);
	if(!window)
		XLAL_ERROR(XLAL_EFUNC);
	for(i = 0; i < h->data->length; i++)
		h->data->data[i] *= window->data->data[i];
	XLALDestroyREAL8Window(window);

	/* add source time series to target time series */

	if(!XLALAddREAL8TimeSeries(target, h))
		XLAL_ERROR(XLAL_EFUNC);

	/* done */

	return 0;
}



/**
 * Essentially a wrapper for XLALAddREAL4TimeSeries(), but performs
 * sub-sample re-interpolation to adjust the source time series epoch to
 * lie on an integer sample boundary in the target time series.  This
 * transformation is done in the frequency domain, so it is convenient to
 * allow a response function to be applied at the same time.  Passing NULL
 * for the response function turns this feature off (i.e., uses a unit
 * response).
 *
 * NOTE:  the source time series is modified in place by this function!
 *
 * This function accepts source and target time series whose units are not
 * the same, and allows the two time series to be herterodyned (although it
 * currently requires them to have the same heterodyne frequency).
 */


int XLALSimAddInjectionREAL4TimeSeries(
	REAL4TimeSeries *target,
	REAL4TimeSeries *h,
	const COMPLEX8FrequencySeries *response
)
{
	/* 1 ns is about 10^-5 samples at 16384 Hz */
	const double noop_threshold = 1e-4;	/* samples */
	REAL4Window *window;
	/* the source time series is padded with at least this many 0's at
	 * the start and end before re-interpolation in an attempt to
	 * suppress aperiodicity artifacts, and 1/2 this many samples is
	 * clipped from the start and end afterwards */
	const unsigned aperiodicity_suppression_buffer = 32768;
	unsigned i;
	double start_sample_int;
	double start_sample_frac;

	/* check input */

	/* FIXME:  since we route the source time series through the
	 * frequency domain, it might not be hard to adjust the heterodyne
	 * frequency and sample rate to match the target time series
	 * instead of making it an error if they don't match. */

	if(h->deltaT != target->deltaT || h->f0 != target->f0) {
		XLALPrintError("%s(): error: input sample rates or heterodyne frequencies do not match\n", __func__);
		XLAL_ERROR(XLAL_EINVAL);
	}

	/* extend the source time series by adding the "aperiodicity
	 * padding" to the start and end.  for efficiency's sake, make sure
	 * the new length is a power of two. */

	i = round_up_to_power_of_two(h->data->length + 2 * aperiodicity_suppression_buffer);
	if(i < h->data->length) {
		/* integer overflow */
		XLALPrintError("%s(): error: source time series too long\n", __func__);
		XLAL_ERROR(XLAL_EBADLEN);
	}
	i -= h->data->length;
	if(!XLALResizeREAL4TimeSeries(h, -(int) (i / 2), h->data->length + i))
		XLAL_ERROR(XLAL_EFUNC);

	/* compute the integer and fractional parts of the sample index in
	 * the target time series on which the source time series begins.
	 * modf() returns integer and fractional parts that have the same
	 * sign, e.g. -3.9 --> -3 + -0.9.  we adjust these so that the
	 * magnitude of the fractional part is not greater than 0.5, e.g.
	 * -3.9 --> -4 + 0.1, so that we never do more than 1/2 a sample of
	 * re-interpolation.  I don't know if really makes any difference,
	 * though */

	start_sample_frac = modf(XLALGPSDiff(&h->epoch, &target->epoch) / target->deltaT, &start_sample_int);
	if(start_sample_frac < -0.5) {
		start_sample_frac += 1.0;
		start_sample_int -= 1.0;
	} else if(start_sample_frac > +0.5) {
		start_sample_frac -= 1.0;
		start_sample_int += 1.0;
	}

	if(fabs(start_sample_frac) > noop_threshold) {
		COMPLEX8FrequencySeries *tilde_h;
		REAL4FFTPlan *plan;

		/* transform source time series to frequency domain.  the FFT
		 * function populates the frequency series' metadata with the
		 * appropriate values. */

		tilde_h = XLALCreateCOMPLEX8FrequencySeries(NULL, &h->epoch, 0, 0, &lalDimensionlessUnit, h->data->length / 2 + 1);
		plan = XLALCreateForwardREAL4FFTPlan(h->data->length, 0);
		if(!tilde_h || !plan) {
			XLALDestroyCOMPLEX8FrequencySeries(tilde_h);
			XLALDestroyREAL4FFTPlan(plan);
			XLAL_ERROR(XLAL_EFUNC);
		}
		i = XLALREAL4TimeFreqFFT(tilde_h, h, plan);
		XLALDestroyREAL4FFTPlan(plan);
		if(i) {
			XLALDestroyCOMPLEX8FrequencySeries(tilde_h);
			XLAL_ERROR(XLAL_EFUNC);
		}

		/* apply sub-sample time correction and optional response function
		 * */

		for(i = 0; i < tilde_h->data->length; i++) {
			const double f = tilde_h->f0 + i * tilde_h->deltaF;
			COMPLEX8 fac;

			/* phase for sub-sample time correction */

			fac = cexp(-I * LAL_TWOPI * f * start_sample_frac * target->deltaT);

			/* divide the source by the response function.  if a
			 * frequency is required that lies outside the domain of
			 * definition of the response function, then the response
			 * is assumed equal to its value at the nearest edge of the
			 * domain of definition.  within the domain of definition,
			 * frequencies are rounded to the nearest bin.  if the
			 * response function is zero in some bin, then the source
			 * data is zeroed in that bin (instead of dividing by 0).
			 * */

			/* FIXME:  should we use GSL to construct an interpolator
			 * for the modulus and phase as functions of frequency, and
			 * use that to evaluate the response?  instead of rounding
			 * to nearest bin? */

			if(response) {
				int j = floor((f - response->f0) / response->deltaF + 0.5);
				if(j < 0)
					j = 0;
				else if((unsigned) j > response->data->length - 1)
					j = response->data->length - 1;
				if(response->data->data[j] == 0.0)
					fac = 0.0;
				else
					fac /= response->data->data[j];
			}

			/* apply factor */

			tilde_h->data->data[i] *= fac;
		}

		/* adjust DC and Nyquist components.  the DC component must always
		 * be real-valued.  because we have adjusted the source time series
		 * to have a length that is an even integer (we've made it a power
		 * of 2) the Nyquist component must also be real valued. */

		if(response) {
			/* a response function has been provided.  zero the DC and
			 * Nyquist components */
			if(tilde_h->f0 == 0.0)
				tilde_h->data->data[0] = 0.0;
			tilde_h->data->data[tilde_h->data->length - 1] = 0.0;
		} else {
			/* no response has been provided.  set the phase of the DC
			 * component to 0, set the imaginary component of the
			 * Nyquist to 0 */
			if(tilde_h->f0 == 0.0)
				tilde_h->data->data[0] = cabsf(tilde_h->data->data[0]);
			tilde_h->data->data[tilde_h->data->length - 1] = crealf(tilde_h->data->data[tilde_h->data->length - 1]);
		}

		/* return to time domain */

		plan = XLALCreateReverseREAL4FFTPlan(h->data->length, 0);
		if(!plan) {
			XLALDestroyCOMPLEX8FrequencySeries(tilde_h);
			XLAL_ERROR(XLAL_EFUNC);
		}
		i = XLALREAL4FreqTimeFFT(h, tilde_h, plan);
		XLALDestroyREAL4FFTPlan(plan);
		XLALDestroyCOMPLEX8FrequencySeries(tilde_h);
		if(i)
			XLAL_ERROR(XLAL_EFUNC);

		/* the deltaT can get "corrupted" by floating point round-off
		 * during its trip through the frequency domain.  since this
		 * function starts by confirming that the sample rate of the source
		 * matches that of the target time series, we can use the target
		 * series' sample rate to reset the source's sample rate to its
		 * original value.  but we do a check to make sure we're not
		 * masking a real bug */

		if(fabs(h->deltaT - target->deltaT) / target->deltaT > 1e-12) {
			XLALPrintError("%s(): error: oops, internal sample rate mismatch\n", __func__);
			XLAL_ERROR(XLAL_EERR);
		}
		h->deltaT = target->deltaT;
	}

	/* set source epoch from target epoch and integer sample offset */

	h->epoch = target->epoch;
	XLALGPSAdd(&h->epoch, start_sample_int * target->deltaT);

	/* clip half of the "aperiodicity padding" from the start and end
	 * of the source time series in a continuing effort to suppress
	 * aperiodicity artifacts. */

	if(!XLALResizeREAL4TimeSeries(h, aperiodicity_suppression_buffer / 2, h->data->length - aperiodicity_suppression_buffer))
		XLAL_ERROR(XLAL_EFUNC);

	/* apply a Tukey window whose tapers lie within the remaining
	 * aperiodicity padding. leaving one sample of the aperiodicty
	 * padding untouched on each side of the original time series
	 * because the data might have been shifted into it */

	window = XLALCreateTukeyREAL4Window(h->data->length, (double) (aperiodicity_suppression_buffer - 2) / h->data->length);
	if(!window)
		XLAL_ERROR(XLAL_EFUNC);
	for(i = 0; i < h->data->length; i++)
		h->data->data[i] *= window->data->data[i];
	XLALDestroyREAL4Window(window);

	/* add source time series to target time series */

	if(!XLALAddREAL4TimeSeries(target, h))
		XLAL_ERROR(XLAL_EFUNC);

	/* done */

	return 0;
}



/* TODO: ROUTINES SHOULD GO INTO DETRESPONSE.[CH] IN LAL/PACKAGES/TOOLS */ 


/*
 *
 * beta = pi f L / c
 * mu = k . u
 *
 * @sa
 * John T. Whelan, "Higher-Frequency Corrections to Stochastic Formulae",
 * LIGO-T070172.
 * @sa
 * Louis J. Rubbo, Neil J. Cornish, and Olivier Poujade, "Forward modeling of
 * space-borne gravitational wave detectors", Phys. Rev. D 69, 082003 (2004);
 * arXiv:gr-qc/0311069.
 * @sa
 * Malik Rakhmanov, "Response of LIGO to Gravitational Waves at High
 * Frequencies and in the Vicinity of the FSR (37.5 kHz)", LIGO-T060237.
 */
static COMPLEX16 XLALComputeDetArmTransferFunction(double beta, double mu)
{
	COMPLEX16 ans;
	ans = cexp(I * beta * (1.0 - mu)) * gsl_sf_sinc(beta * (1.0 + mu));
	ans += cexp(-I * beta * (1.0 + mu)) * gsl_sf_sinc(beta * (1.0 - mu));
	ans *= 0.5;
	return ans;
}


static void getarm(double u[3], double alt, double azi, double lat, double lon)
{
	double cosalt = cos(alt);
	double sinalt = sin(alt);
	double cosazi = cos(azi);
	double sinazi = sin(azi);
	double coslat = cos(lat);
	double sinlat = sin(lat);
	double coslon = cos(lon);
	double sinlon = sin(lon);
	double uNorth = cosalt * cosazi;
	double uEast = cosalt * sinazi;
	double uUp = sinalt;
	double uRho = - sinlat * uNorth + coslat * uUp;
	u[0] = coslon * uRho - sinlon * uEast;
	u[1] = sinlon * uRho + coslon * uEast;
	u[2] = coslat * uNorth + sinlat * uUp;
	return;
}

static void XLALComputeDetAMResponseParts(double *armlen, double *xcos, double *ycos, double *fxplus, double *fyplus, double *fxcross, double *fycross, const LALDetector *detector, double ra, double dec, double psi, double gmst)
{
	double X[3];	/* wave frame x axis */
	double Y[3];	/* wave frame y axis */
	double Z[3];	/* wave frame z axis (propagation direction) */
	double U[3];	/* x arm unit vector */
	double V[3];	/* y arm unit vector */
	double DU[3][3];	/* single arm response tensor for x arm */
	double DV[3][3];	/* single arm response tensor for y arm */
	double gha = gmst - ra;	/* greenwich hour angle */
	double cosgha = cos(gha);
	double singha = sin(gha);
	double cosdec = cos(dec);
	double sindec = sin(dec);
	double cospsi = cos(psi);
	double sinpsi = sin(psi);
	int i, j;

	/* compute unit vectors specifying the wave frame x, y, and z axes */

	X[0] = -cospsi * singha - sinpsi * cosgha * sindec;
	X[1] = -cospsi * cosgha + sinpsi * singha * sindec;
	X[2] =  sinpsi * cosdec;
	Y[0] =  sinpsi * singha - cospsi * cosgha * sindec;
	Y[1] =  sinpsi * cosgha + cospsi * singha * sindec;
	Y[2] =  cospsi * cosdec;
	Z[0] = -cosgha * cosdec;
	Z[1] =  singha * cosdec;
	Z[2] = -sindec;

	switch (detector->type) {

	case LALDETECTORTYPE_IFOCOMM:
	case LALDETECTORTYPE_IFODIFF:
	
		/* FIXME: should compute the effect of non-equal arm lengths;
		 * but, for now, just use the mean arm length */

		*armlen = detector->frDetector.xArmMidpoint
			+ detector->frDetector.yArmMidpoint;

		/* get the unit vectors along the arms */

		getarm(U, detector->frDetector.xArmAltitudeRadians,
			detector->frDetector.xArmAzimuthRadians,
			detector->frDetector.vertexLatitudeRadians,
			detector->frDetector.vertexLongitudeRadians);

		getarm(V, detector->frDetector.yArmAltitudeRadians,
			detector->frDetector.yArmAzimuthRadians,
			detector->frDetector.vertexLatitudeRadians,
			detector->frDetector.vertexLongitudeRadians);

		/* compute direction cosines for the signal direction relative
         	 * to the x-arm and the y-arm */

		*xcos = *ycos = 0.0;
		for (i = 0; i < 3; ++i) {
			*xcos += U[i] * Z[i];
			*ycos += V[i] * Z[i];
		}

		/* compute the single arm response tensors for the x-arm and
		 * y-arm */

		for (i = 0; i < 3; ++i) {
			DU[i][i] = 0.5 * U[i] * U[i];
			DV[i][i] = 0.5 * V[i] * V[i];
			for (j = i + 1; j < 3; ++j) {
				DU[i][j] = DU[j][i] = 0.5 * U[i] * U[j];
				DV[i][j] = DV[j][i] = 0.5 * V[i] * V[j];
			}
		}

		/* compute the beam pattern partial responses for the x-arm and
		 * y-arm */

		*fxplus = *fxcross = 0.0;
		*fyplus = *fycross = 0.0;
		for (i = 0; i < 3; ++i) {
			double DUX = DU[i][0]*X[0]+DU[i][1]*X[1]+DU[i][2]*X[2];
			double DUY = DU[i][0]*Y[0]+DU[i][1]*Y[1]+DU[i][2]*Y[2];
			double DVX = DV[i][0]*X[0]+DV[i][1]*X[1]+DV[i][2]*X[2];
			double DVY = DV[i][0]*Y[0]+DV[i][1]*Y[1]+DV[i][2]*Y[2];
			*fxplus  += X[i] * DUX - Y[i] * DUY;
			*fxcross += X[i] * DUY + Y[i] * DUX;
			*fyplus  += X[i] * DVX - Y[i] * DVY;
			*fycross += X[i] * DVY + Y[i] * DVX;
		}

		/* differential interferometer: arm y is subtracted from
		 * arm x */
		if (detector->type == LALDETECTORTYPE_IFODIFF) {
			*fyplus *= -1;
			*fycross *= -1;
		}

		break;

	case LALDETECTORTYPE_IFOXARM:

		/* no y-arm */

		*armlen = 2.0 * detector->frDetector.xArmMidpoint;

		getarm(U, detector->frDetector.xArmAltitudeRadians,
			detector->frDetector.xArmAzimuthRadians,
			detector->frDetector.vertexLatitudeRadians,
			detector->frDetector.vertexLongitudeRadians);

		*xcos = *ycos = 0.0;
		for (i = 0; i < 3; ++i)
			*xcos += U[i] * Z[i];

		*fyplus = *fycross = 0.0;
		XLALComputeDetAMResponse(fxplus, fxcross,
			(const REAL4(*)[3])detector->response, ra, dec, psi,
			gmst);

		break;

	case LALDETECTORTYPE_IFOYARM:

		/* no x-arm */

		*armlen = 2.0 * detector->frDetector.yArmMidpoint;

		getarm(V, detector->frDetector.yArmAltitudeRadians,
			detector->frDetector.yArmAzimuthRadians,
			detector->frDetector.vertexLatitudeRadians,
			detector->frDetector.vertexLongitudeRadians);

		*xcos = *ycos = 0.0;
		for (i = 0; i < 3; ++i)
			*ycos += V[i] * Z[i];

		*fxplus = *fxcross = 0.0;
		XLALComputeDetAMResponse(fyplus, fycross,
			(const REAL4(*)[3])detector->response, ra, dec, psi,
			gmst);

		break;

	default:

		/* FIXME: could handle this situation properly; fur now, just
		 * ignore non long-wavelength-limit effects by setting armlen
		 * to zero; also, pretend that all of the response is
		 * associated with the x-arm */

		*armlen = *xcos = *ycos = 0.0;
		*fyplus = *fycross = 0.0;
		XLALComputeDetAMResponse(fxplus, fxcross,
			(const REAL4(*)[3])detector->response, ra, dec, psi,
			gmst);

		break;

	}

	return;
}


/* Helper routine that computes a segment of strain data with a single
 * time delay and beam pattern applied to the whole segment.  The duration
 * of the segment must therefore be reasonably short or else the movement
 * of the earth will invalidate the use of a single time shift and beam
 * pattern for the entire segment. */
static int XLALSimComputeStrainSegmentREAL8TimeSeries(
	REAL8TimeSeries *segment,
	const REAL8TimeSeries *hplus,
	const REAL8TimeSeries *hcross,
	COMPLEX16FrequencySeries *work1,
	COMPLEX16FrequencySeries *work2,
	REAL8FFTPlan *fwdplan,
	REAL8FFTPlan *revplan,
	REAL8Window *window,
	double ra,
	double dec,
	double psi,
	LALDetector *detector,
	const COMPLEX16FrequencySeries *response
)
{
	LIGOTimeGPS t;
	double gmst;
	double xcos;
	double ycos;
	double fxplus;
	double fyplus;
	double fxcross;
	double fycross;
	double armlen;
	double deltaT;
	double offint;
	double offrac;
	int offset;
	int j;
	size_t k;

	/* this routine assumes the segment has a length of N points where N is
	 * a power of two and that the workspace frequency series and the FFT
	 * plans are compatible with the size N; these assumptions are not
	 * checked: the calling routine must ensure that they are true */

	/* compute fplus, fcross, and time delay from earth's center at the
 	 * time corresponding to the middle of the segment */

	t = segment->epoch;
	XLALGPSAdd(&t, 0.5 * segment->data->length * segment->deltaT);
	gmst = XLALGreenwichMeanSiderealTime(&t);
	XLALComputeDetAMResponseParts(&armlen, &xcos, &ycos, &fxplus, &fyplus,
		&fxcross, &fycross, detector, ra, dec, psi, gmst);
	deltaT = XLALTimeDelayFromEarthCenter(detector->location, ra, dec, &t);

	/* add to the geometric delay the difference in time between the
	 * beginning of the injection timeseries and the beginning of the
	 * segment */

	deltaT += XLALGPSDiff(&hplus->epoch, &segment->epoch);

	/* compute the integer and fractional parts of the sample index in the
	 * segment on which the hplus and hcross time series begins: modf()
	 * returns integer and fractional parts that have the same sign, e.g.,
	 * -3.9 --> -3 + -0.9, and we adjust these so that magnitude of the
	 * fractional part is not greater than 0.5, e.g., -3.9 --> -4.0 + 0.1,
	 * so that we never do more than 1/2 a sample of re-interpolation */

	offrac = modf(deltaT / segment->deltaT, &offint);
	if (offrac < -0.5) {
		offrac += 1.0;
		offint -= 1.0;
	} else if (offrac > 0.5) {
		offrac -= 1.0;
		offint += 1.0;
	}
	offset = offint;

	/* now compute the sub-sample time shift that must be applied to the
	 * segment data */

	deltaT = offrac * segment->deltaT;

	/* window the date and put it in frequency domain */

	for (j = 0; j < (int)segment->data->length; ++j)
		if (j >= offset && j < (int)hplus->data->length + offset) {
			segment->data->data[j] = window->data->data[j]
				* hplus->data->data[j - offset];
		} else
			segment->data->data[j] = 0.0;

	if (XLALREAL8TimeFreqFFT(work1, segment, fwdplan) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	for (j = 0; j < (int)segment->data->length; ++j)
		if (j >= offset && j < (int)hcross->data->length + offset) {
			segment->data->data[j] = window->data->data[j]
				* hcross->data->data[j - offset];
		} else
			segment->data->data[j] = 0.0;

	if (XLALREAL8TimeFreqFFT(work2, segment, fwdplan) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	/* apply sub-sample time shift in frequency domain */

	for (k = 0; k < work1->data->length; ++k) {
		double f = work1->f0 + k * work1->deltaF;
		double beta = LAL_PI * f * armlen / LAL_C_SI;
		COMPLEX16 Tx, Ty; /* x- and y-arm transfer functions */
		COMPLEX16 gplus, gcross;
		COMPLEX16 fac;
		
		/* phase for sub-sample time correction */
		fac = cexp(-I * LAL_TWOPI * f * deltaT);
		if (response)
			fac /= response->data->data[k];

		Tx = XLALComputeDetArmTransferFunction(beta, xcos);
		Ty = XLALComputeDetArmTransferFunction(beta, ycos);
		gplus = Tx * fxplus + Ty * fyplus;
		gcross = Tx * fxcross + Ty * fycross;

		work1->data->data[k] *= gplus;
		work1->data->data[k] += gcross * work2->data->data[k];
		work1->data->data[k] *= fac;
	}

	/* adjust DC and Nyquist components: the DC component must always be
	 * real-valued; because the calling routine has made the time series
	 * have an even length, the Nyquist component must also be real-valued;
	 * also this routine makes the assumption that both the DC and the
	 * Nyquist components are zero */

	work1->data->data[0] = cabs(work1->data->data[0]);
	work1->data->data[work1->data->length - 1] =
		creal(work1->data->data[work1->data->length - 1]);
	
	/* return data to time domain */

	if (XLALREAL8FreqTimeFFT(segment, work1, revplan) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	return 0;
}

/* Helper routine that computes a segment of strain data with a single
 * time delay and beam pattern applied to the whole segment.  The duration
 * of the segment must therefore be reasonably short or else the movement
 * of the earth will invalidate the use of a single time shift and beam
 * pattern for the entire segment. */
static int XLALSimComputeStrainSegmentREAL4TimeSeries(
	REAL4TimeSeries *segment,
	const REAL4TimeSeries *hplus,
	const REAL4TimeSeries *hcross,
	COMPLEX8FrequencySeries *work1,
	COMPLEX8FrequencySeries *work2,
	REAL4FFTPlan *fwdplan,
	REAL4FFTPlan *revplan,
	REAL4Window *window,
	double ra,
	double dec,
	double psi,
	LALDetector *detector,
	const COMPLEX8FrequencySeries *response
)
{
	LIGOTimeGPS t;
	double gmst;
	double xcos;
	double ycos;
	double fxplus;
	double fyplus;
	double fxcross;
	double fycross;
	double armlen;
	double deltaT;
	double offint;
	double offrac;
	int offset;
	int j;
	size_t k;

	/* this routine assumes the segment has a length of N points where N is
	 * a power of two and that the workspace frequency series and the FFT
	 * plans are compatible with the size N; these assumptions are not
	 * checked: the calling routine must ensure that they are true */

	/* compute fplus, fcross, and time delay from earth's center at the
 	 * time corresponding to the middle of the segment */

	t = segment->epoch;
	XLALGPSAdd(&t, 0.5 * segment->data->length * segment->deltaT);
	gmst = XLALGreenwichMeanSiderealTime(&t);
	XLALComputeDetAMResponseParts(&armlen, &xcos, &ycos, &fxplus, &fyplus,
		&fxcross, &fycross, detector, ra, dec, psi, gmst);
	deltaT = XLALTimeDelayFromEarthCenter(detector->location, ra, dec, &t);

	/* add to the geometric delay the difference in time between the
	 * beginning of the injection timeseries and the beginning of the
	 * segment */

	deltaT += XLALGPSDiff(&hplus->epoch, &segment->epoch);

	/* compute the integer and fractional parts of the sample index in the
	 * segment on which the hplus and hcross time series begins: modf()
	 * returns integer and fractional parts that have the same sign, e.g.,
	 * -3.9 --> -3 + -0.9, and we adjust these so that magnitude of the
	 * fractional part is not greater than 0.5, e.g., -3.9 --> -4.0 + 0.1,
	 * so that we never do more than 1/2 a sample of re-interpolation */

	offrac = modf(deltaT / segment->deltaT, &offint);
	if (offrac < -0.5) {
		offrac += 1.0;
		offint -= 1.0;
	} else if (offrac > 0.5) {
		offrac -= 1.0;
		offint += 1.0;
	}
	offset = offint;

	/* now compute the sub-sample time shift that must be applied to the
	 * segment data */

	deltaT = offrac * segment->deltaT;

	/* window the date and put it in frequency domain */

	for (j = 0; j < (int)segment->data->length; ++j)
		if (j >= offset && j < (int)hplus->data->length + offset) {
			segment->data->data[j] = window->data->data[j]
				* hplus->data->data[j - offset];
		} else
			segment->data->data[j] = 0.0;

	if (XLALREAL4TimeFreqFFT(work1, segment, fwdplan) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	for (j = 0; j < (int)segment->data->length; ++j)
		if (j >= offset && j < (int)hcross->data->length + offset) {
			segment->data->data[j] = window->data->data[j]
				* hcross->data->data[j - offset];
		} else
			segment->data->data[j] = 0.0;

	if (XLALREAL4TimeFreqFFT(work2, segment, fwdplan) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	/* apply sub-sample time shift in frequency domain */

	for (k = 0; k < work1->data->length; ++k) {
		double f = work1->f0 + k * work1->deltaF;
		double beta = LAL_PI * f * armlen / LAL_C_SI;
		COMPLEX16 Tx, Ty; /* x- and y-arm transfer functions */
		COMPLEX16 gplus, gcross;
		COMPLEX16 fac;
		
		/* phase for sub-sample time correction */
		fac = cexp(-I * LAL_TWOPI * f * deltaT);
		if (response)
			fac /= response->data->data[k];

		Tx = XLALComputeDetArmTransferFunction(beta, xcos);
		Ty = XLALComputeDetArmTransferFunction(beta, ycos);
		gplus = Tx * fxplus + Ty * fyplus;
		gcross = Tx * fxcross + Ty * fycross;

		work1->data->data[k] *= gplus;
		work1->data->data[k] += gcross * work2->data->data[k];
		work1->data->data[k] *= fac;
	}

	/* adjust DC and Nyquist components: the DC component must always be
	 * real-valued; because the calling routine has made the time series
	 * have an even length, the Nyquist component must also be real-valued;
	 * also this routine makes the assumption that both the DC and the
	 * Nyquist components are zero */

	work1->data->data[0] = cabs(work1->data->data[0]);
	work1->data->data[work1->data->length - 1] =
		creal(work1->data->data[work1->data->length - 1]);
	
	/* return data to time domain */

	if (XLALREAL4FreqTimeFFT(segment, work1, revplan) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	return 0;
}

/**
 * @brief Computes strain for a detector and injects into target time series.
 * @details This routine takes care of the time-changing time delay from
 * the Earth's center and the time-changing antenna response pattern; it
 * also accounts for deviations from the long-wavelength limit at high
 * frequencies.  An optional calibration response function can be provided
 * if the output time series is not in strain units.
 * @param[in,out] target Time series to inject strain into.
 * @param[in] hplus Time series with plus-polarization gravitational waveform.
 * @param[in] hcross Time series with cross-polarization gravitational waveform.
 * @param[in] ra Right ascension of the source (radians).
 * @param[in] dec Declination of the source (radians).
 * @param[in] psi Polarization angle of the source (radians).
 * @param[in] detector Detector to use when computing strain.
 * @param[in] response Response function to use, or NULL if none.
 * @retval 0 Success.
 * @retval <0 Failure.
 */
int XLALSimInjectDetectorStrainREAL8TimeSeries(
	REAL8TimeSeries *target,
	const REAL8TimeSeries *hplus,
	const REAL8TimeSeries *hcross,
	double ra,
	double dec,
	double psi,
	LALDetector *detector,
	const COMPLEX16FrequencySeries *response
)
{
	const double nominal_segdur = 2.0; /* nominal segment duration = 2s */
	const double max_time_delay = 0.1; /* generous allowed time delay */
	const size_t strides_per_segment = 2; /* 2 strides in one segment */
	LIGOTimeGPS t0;
	LIGOTimeGPS t1;
	size_t length;		/* length in samples of interval t0 - t1 */
	size_t seglen;		/* length of segment in samples */
	size_t padlen;		/* padding at beginning and end of segment */
	size_t ovrlap;		/* overlapping data length */
	size_t stride;		/* stride of each step */
	size_t nsteps;		/* number of steps to take */
	REAL8TimeSeries *h = NULL; /* strain timeseries to inject into target */
	REAL8TimeSeries *segment = NULL;
	COMPLEX16FrequencySeries *work1 = NULL;
	COMPLEX16FrequencySeries *work2 = NULL;
	REAL8FFTPlan *fwdplan = NULL;
	REAL8FFTPlan *revplan = NULL;
	REAL8Window *window = NULL;
	size_t step;
	size_t j;
	int errnum = 0;

	/* check validity and compatibility of time series */

	LAL_CHECK_VALID_SERIES(target, XLAL_FAILURE);
	LAL_CHECK_VALID_SERIES(hplus, XLAL_FAILURE);
	LAL_CHECK_VALID_SERIES(hcross, XLAL_FAILURE);
	LAL_CHECK_CONSISTENT_TIME_SERIES(hplus, hcross, XLAL_FAILURE);
	if (response == NULL) {
		LAL_CHECK_COMPATIBLE_TIME_SERIES(target, hplus, XLAL_FAILURE);
	} else {
		/* units do no need to agree, but sample interval and
		 * start frequency do */
		if (fabs(target->deltaT - hplus->deltaT ) > LAL_REAL8_EPS)
			XLAL_ERROR(XLAL_ETIME);
		if (fabs(target->f0 - hplus->f0) > LAL_REAL8_EPS)
			XLAL_ERROR(XLAL_EFREQ);
	}
			

	/* constants describing the data segmentation: the length of the
	 * segment must be a power of two and the segment duration is at least
	 * nominal_segdur */

	seglen = round_up_to_power_of_two(nominal_segdur / target->deltaT);
	stride = seglen / strides_per_segment;
	padlen = max_time_delay / target->deltaT;
	ovrlap = seglen;
	ovrlap -= 2 * padlen;
	ovrlap -= stride;

	/* determine start and end time: the start time is the later of the
	 * start of the hplus/hcross time series and the target time series;
	 * the end time is the earlier of the end of the hplus/hcross time
	 * series and the target time series */

	t0 = hplus->epoch;
	t1 = target->epoch;
	XLALGPSAdd(&t0, hplus->data->length * hplus->deltaT);
	XLALGPSAdd(&t1, target->data->length * target->deltaT);
	t1 = XLALGPSCmp(&t1, &t0) < 0 ? t1 : t0;
	t0 = hplus->epoch;
	t0 = XLALGPSCmp(&t0, &target->epoch) > 0 ? t0 : target->epoch;

	/* add padding of 1 stride before and after these start and end times */

	XLALGPSAdd(&t0, -1.0 * stride * target->deltaT);
	XLALGPSAdd(&t1, stride * target->deltaT);

	/* determine if this is a disjoint set: if so, there is nothing to do */

	if (XLALGPSCmp(&t1, &t0) <= 0)
		return 0;

	/* create a segment that is seglen samples long */

	segment = XLALCreateREAL8TimeSeries(NULL, &t0, target->f0,
		target->deltaT, &target->sampleUnits, seglen);
	if (!segment) {
		errnum = XLAL_EFUNC;
		goto freereturn;
	}

	/* create a time series to hold the strain to inject into the target */

	length = XLALGPSDiff(&t1, &t0) / target->deltaT;
	h = XLALCreateREAL8TimeSeries(NULL, &t0, target->f0, target->deltaT,
		&target->sampleUnits, length);
	if (!h) {
		errnum = XLAL_EFUNC;
		goto freereturn;
	}
	memset(h->data->data, 0, h->data->length * sizeof(*h->data->data));

	/* determine number of steps it takes to go from t0 to t1 */

	nsteps = ((length%stride) ? (1 + length/stride) : (length/stride));


	/* create frequency-domain workspace; note that the FFT function
	 * populates the frequency series' metadata with the appropriate
	 * values */

	work1 = XLALCreateCOMPLEX16FrequencySeries(NULL, &t0, 0, 0,
		&lalDimensionlessUnit, seglen / 2 + 1);
	work2 = XLALCreateCOMPLEX16FrequencySeries(NULL, &t0, 0, 0,
		&lalDimensionlessUnit, seglen / 2 + 1);
	if (!work1 || !work2) {
		errnum = XLAL_EFUNC;
		goto freereturn;
	}

	/* create forward and reverse FFT plans */

	fwdplan = XLALCreateForwardREAL8FFTPlan(seglen, 0);
	revplan = XLALCreateReverseREAL8FFTPlan(seglen, 0);
	if (!fwdplan || !revplan) {
		errnum = XLAL_EFUNC;
		goto freereturn;
	}

	/* create a Tukey window with tapers entirely within the padding */

	window = XLALCreateTukeyREAL8Window(seglen, (double)padlen / seglen); 


	/* loop over steps, adding data from the current step to the strain */

	for (step = 0; step < nsteps; ++step) {

		int status;
		size_t offset;

		/* compute one segment of strain with time appropriate beam
 		 * pattern functions and time delays from earth's center */

		status = XLALSimComputeStrainSegmentREAL8TimeSeries(segment,				hplus, hcross, work1, work2, fwdplan, revplan, window,
			ra, dec, psi, detector, response);
		if (status < 0) {
			errnum = XLAL_EFUNC;
			goto freereturn;
		}

		/* compute the offset of this segment relative to the strain series */

		offset = XLALGPSDiff(&segment->epoch, &h->epoch) / h->deltaT;
		for (j = padlen; j < seglen - padlen; ++j)
			if ((j + offset) < h->data->length) {
				if (step && j - padlen < ovrlap) {
					/* feather overlapping data */
					double x = (double)(j - padlen) / ovrlap;
					h->data->data[j + offset] = x * segment->data->data[j]
						+ (1.0 - x) * h->data->data[j + offset];
				} else /* no feathering of remaining data */
					h->data->data[j + offset] = segment->data->data[j];
			}

		/* advance segment start time the next step */

		XLALGPSAdd(&segment->epoch, stride * segment->deltaT);
	}

	/* apply window to beginning and end of time series to reduce ringing */

	for (j = 0; j < stride - padlen; ++j)
		h->data->data[j] = h->data->data[h->data->length - 1 - j] = 0.0;
	for ( ; j < stride; ++j) {
		double fac = window->data->data[j - (stride - padlen)];
		h->data->data[j] *= fac;
		h->data->data[h->data->length - 1 - j] *= fac;
	}
	

	/* add computed strain to target time series */

	XLALAddREAL8TimeSeries(target, h);

freereturn:

	/* free all memory and return */

	XLALDestroyREAL8Window(window);
	XLALDestroyREAL8FFTPlan(revplan);
	XLALDestroyREAL8FFTPlan(fwdplan);
	XLALDestroyCOMPLEX16FrequencySeries(work2);
	XLALDestroyCOMPLEX16FrequencySeries(work1);
	XLALDestroyREAL8TimeSeries(h);
	XLALDestroyREAL8TimeSeries(segment);

	if (errnum)
		XLAL_ERROR(errnum);
	return 0;
}

/**
 * @brief Computes strain for a detector and injects into target time series.
 * @details This routine takes care of the time-changing time delay from
 * the Earth's center and the time-changing antenna response pattern; it
 * also accounts for deviations from the long-wavelength limit at high
 * frequencies.  An optional calibration response function can be provided
 * if the output time series is not in strain units.
 * @param[in,out] target Time series to inject strain into.
 * @param[in] hplus Time series with plus-polarization gravitational waveform.
 * @param[in] hcross Time series with cross-polarization gravitational waveform.
 * @param[in] ra Right ascension of the source (radians).
 * @param[in] dec Declination of the source (radians).
 * @param[in] psi Polarization angle of the source (radians).
 * @param[in] detector Detector to use when computing strain.
 * @param[in] response Response function to use, or NULL if none.
 * @retval 0 Success.
 * @retval <0 Failure.
 */
int XLALSimInjectDetectorStrainREAL4TimeSeries(
	REAL4TimeSeries *target,
	const REAL4TimeSeries *hplus,
	const REAL4TimeSeries *hcross,
	double ra,
	double dec,
	double psi,
	LALDetector *detector,
	const COMPLEX8FrequencySeries *response
)
{
	const double nominal_segdur = 2.0; /* nominal segment duration = 2s */
	const double max_time_delay = 0.1; /* generous allowed time delay */
	const size_t strides_per_segment = 2; /* 2 strides in one segment */
	LIGOTimeGPS t0;
	LIGOTimeGPS t1;
	size_t length;		/* length in samples of interval t0 - t1 */
	size_t seglen;		/* length of segment in samples */
	size_t padlen;		/* padding at beginning and end of segment */
	size_t ovrlap;		/* overlapping data length */
	size_t stride;		/* stride of each step */
	size_t nsteps;		/* number of steps to take */
	REAL4TimeSeries *h = NULL; /* strain timeseries to inject into target */
	REAL4TimeSeries *segment = NULL;
	COMPLEX8FrequencySeries *work1 = NULL;
	COMPLEX8FrequencySeries *work2 = NULL;
	REAL4FFTPlan *fwdplan = NULL;
	REAL4FFTPlan *revplan = NULL;
	REAL4Window *window = NULL;
	size_t step;
	size_t j;
	int errnum = 0;

	/* check validity and compatibility of time series */

	LAL_CHECK_VALID_SERIES(target, XLAL_FAILURE);
	LAL_CHECK_VALID_SERIES(hplus, XLAL_FAILURE);
	LAL_CHECK_VALID_SERIES(hcross, XLAL_FAILURE);
	LAL_CHECK_CONSISTENT_TIME_SERIES(hplus, hcross, XLAL_FAILURE);
	if (response == NULL) {
		LAL_CHECK_COMPATIBLE_TIME_SERIES(target, hplus, XLAL_FAILURE);
	} else {
		/* units do no need to agree, but sample interval and
		 * start frequency do */
		if (fabs(target->deltaT - hplus->deltaT ) > LAL_REAL8_EPS)
			XLAL_ERROR(XLAL_ETIME);
		if (fabs(target->f0 - hplus->f0) > LAL_REAL8_EPS)
			XLAL_ERROR(XLAL_EFREQ);
	}
			

	/* constants describing the data segmentation: the length of the
	 * segment must be a power of two and the segment duration is at least
	 * nominal_segdur */

	seglen = round_up_to_power_of_two(nominal_segdur / target->deltaT);
	stride = seglen / strides_per_segment;
	padlen = max_time_delay / target->deltaT;
	ovrlap = seglen;
	ovrlap -= 2 * padlen;
	ovrlap -= stride;

	/* determine start and end time: the start time is the later of the
	 * start of the hplus/hcross time series and the target time series;
	 * the end time is the earlier of the end of the hplus/hcross time
	 * series and the target time series */

	t0 = hplus->epoch;
	t1 = target->epoch;
	XLALGPSAdd(&t0, hplus->data->length * hplus->deltaT);
	XLALGPSAdd(&t1, target->data->length * target->deltaT);
	t1 = XLALGPSCmp(&t1, &t0) < 0 ? t1 : t0;
	t0 = hplus->epoch;
	t0 = XLALGPSCmp(&t0, &target->epoch) > 0 ? t0 : target->epoch;

	/* add padding of 1 stride before and after these start and end times */

	XLALGPSAdd(&t0, -1.0 * stride * target->deltaT);
	XLALGPSAdd(&t1, stride * target->deltaT);

	/* determine if this is a disjoint set: if so, there is nothing to do */

	if (XLALGPSCmp(&t1, &t0) <= 0)
		return 0;

	/* create a segment that is seglen samples long */

	segment = XLALCreateREAL4TimeSeries(NULL, &t0, target->f0,
		target->deltaT, &target->sampleUnits, seglen);
	if (!segment) {
		errnum = XLAL_EFUNC;
		goto freereturn;
	}

	/* create a time series to hold the strain to inject into the target */

	length = XLALGPSDiff(&t1, &t0) / target->deltaT;
	h = XLALCreateREAL4TimeSeries(NULL, &t0, target->f0, target->deltaT,
		&target->sampleUnits, length);
	if (!h) {
		errnum = XLAL_EFUNC;
		goto freereturn;
	}
	memset(h->data->data, 0, h->data->length * sizeof(*h->data->data));

	/* determine number of steps it takes to go from t0 to t1 */

	nsteps = ((length%stride) ? (1 + length/stride) : (length/stride));


	/* create frequency-domain workspace; note that the FFT function
	 * populates the frequency series' metadata with the appropriate
	 * values */

	work1 = XLALCreateCOMPLEX8FrequencySeries(NULL, &t0, 0, 0,
		&lalDimensionlessUnit, seglen / 2 + 1);
	work2 = XLALCreateCOMPLEX8FrequencySeries(NULL, &t0, 0, 0,
		&lalDimensionlessUnit, seglen / 2 + 1);
	if (!work1 || !work2) {
		errnum = XLAL_EFUNC;
		goto freereturn;
	}

	/* create forward and reverse FFT plans */

	fwdplan = XLALCreateForwardREAL4FFTPlan(seglen, 0);
	revplan = XLALCreateReverseREAL4FFTPlan(seglen, 0);
	if (!fwdplan || !revplan) {
		errnum = XLAL_EFUNC;
		goto freereturn;
	}

	/* create a Tukey window with tapers entirely within the padding */

	window = XLALCreateTukeyREAL4Window(seglen, (double)padlen / seglen); 


	/* loop over steps, adding data from the current step to the strain */

	for (step = 0; step < nsteps; ++step) {

		int status;
		size_t offset;

		/* compute one segment of strain with time appropriate beam
 		 * pattern functions and time delays from earth's center */

		status = XLALSimComputeStrainSegmentREAL4TimeSeries(segment,				hplus, hcross, work1, work2, fwdplan, revplan, window,
			ra, dec, psi, detector, response);
		if (status < 0) {
			errnum = XLAL_EFUNC;
			goto freereturn;
		}

		/* compute the offset of this segment relative to the strain series */

		offset = XLALGPSDiff(&segment->epoch, &h->epoch) / h->deltaT;
		for (j = padlen; j < seglen - padlen; ++j)
			if ((j + offset) < h->data->length) {
				if (step && j - padlen < ovrlap) {
					/* feather overlapping data */
					double x = (double)(j - padlen) / ovrlap;
					h->data->data[j + offset] = x * segment->data->data[j]
						+ (1.0 - x) * h->data->data[j + offset];
				} else /* no feathering of remaining data */
					h->data->data[j + offset] = segment->data->data[j];
			}

		/* advance segment start time the next step */

		XLALGPSAdd(&segment->epoch, stride * segment->deltaT);
	}

	/* apply window to beginning and end of time series to reduce ringing */

	for (j = 0; j < stride - padlen; ++j)
		h->data->data[j] = h->data->data[h->data->length - 1 - j] = 0.0;
	for ( ; j < stride; ++j) {
		double fac = window->data->data[j - (stride - padlen)];
		h->data->data[j] *= fac;
		h->data->data[h->data->length - 1 - j] *= fac;
	}
	

	/* add computed strain to target time series */

	XLALAddREAL4TimeSeries(target, h);

freereturn:

	/* free all memory and return */

	XLALDestroyREAL4Window(window);
	XLALDestroyREAL4FFTPlan(revplan);
	XLALDestroyREAL4FFTPlan(fwdplan);
	XLALDestroyCOMPLEX8FrequencySeries(work2);
	XLALDestroyCOMPLEX8FrequencySeries(work1);
	XLALDestroyREAL4TimeSeries(h);
	XLALDestroyREAL4TimeSeries(segment);

	if (errnum)
		XLAL_ERROR(errnum);
	return 0;
}


/*
 * The following routines are more computationally efficient but they
 * assume the long-wavelength limit is valid.
 */


/* Helper routine that computes a segment of strain data with a single
 * time delay and beam pattern applied to the whole segment.  The duration
 * of the segment must therefore be reasonably short or else the movement
 * of the earth will invalidate the use of a single time shift and beam
 * pattern for the entire segment. */
static int XLALSimComputeLWLStrainSegmentREAL8TimeSeries(
	REAL8TimeSeries *segment,
	const REAL8TimeSeries *hplus,
	const REAL8TimeSeries *hcross,
	COMPLEX16FrequencySeries *work,
	REAL8FFTPlan *fwdplan,
	REAL8FFTPlan *revplan,
	REAL8Window *window,
	double ra,
	double dec,
	double psi,
	LALDetector *detector,
	const COMPLEX16FrequencySeries *response
)
{
	LIGOTimeGPS t;
	double gmst;
	double fplus;
	double fcross;
	double deltaT;
	double offint;
	double offrac;
	int offset;
	int j;
	size_t k;

	/* this routine assumes the segment has a length of N points where N is
	 * a power of two and that the workspace frequency series and the FFT
	 * plans are compatible with the size N; these assumptions are not
	 * checked: the calling routine must ensure that they are true */

	/* compute fplus, fcross, and time delay from earth's center at the
 	 * time corresponding to the middle of the segment */

	t = segment->epoch;
	XLALGPSAdd(&t, 0.5 * segment->data->length * segment->deltaT);
	gmst = XLALGreenwichMeanSiderealTime(&t);
	XLALComputeDetAMResponse(&fplus, &fcross,
		(const REAL4(*)[3])detector->response, ra, dec, psi, gmst);
	deltaT = XLALTimeDelayFromEarthCenter(detector->location, ra, dec, &t);

	/* add to the geometric delay the difference in time between the
	 * beginning of the injection timeseries and the beginning of the
	 * segment */

	deltaT += XLALGPSDiff(&hplus->epoch, &segment->epoch);

	/* compute the integer and fractional parts of the sample index in the
	 * segment on which the hplus and hcross time series begins: modf()
	 * returns integer and fractional parts that have the same sign, e.g.,
	 * -3.9 --> -3 + -0.9, and we adjust these so that magnitude of the
	 * fractional part is not greater than 0.5, e.g., -3.9 --> -4.0 + 0.1,
	 * so that we never do more than 1/2 a sample of re-interpolation */

	offrac = modf(deltaT / segment->deltaT, &offint);
	if (offrac < -0.5) {
		offrac += 1.0;
		offint -= 1.0;
	} else if (offrac > 0.5) {
		offrac -= 1.0;
		offint += 1.0;
	}
	offset = offint;

	/* now compute the sub-sample time shift that must be applied to the
	 * segment data */

	deltaT = offrac * segment->deltaT;

	/* compute the strain from hplus and hcross and populate the segment
	 * with this windowed data; apply appropriate integer offset */

	for (j = 0; j < (int)segment->data->length; ++j)
		if (j >= offset && j < (int)hplus->data->length + offset)
			segment->data->data[j] = window->data->data[j]
				* ( fplus * hplus->data->data[j - offset]
				+ fcross * hcross->data->data[j - offset] );
		else
			segment->data->data[j] = 0.0;


	/* put segment in frequency domain */

	if (XLALREAL8TimeFreqFFT(work, segment, fwdplan) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	/* apply sub-sample time shift in frequency domain */

	for (k = 0; k < work->data->length; ++k) {
		double f = work->f0 + k * work->deltaF;
		COMPLEX16 fac;
		
		/* phase for sub-sample time correction */
		fac = cexp(-I * LAL_TWOPI * f * deltaT);
		if (response)
			fac /= response->data->data[k];

		work->data->data[k] *= fac;
	}

	/* adjust DC and Nyquist components: the DC component must always be
	 * real-valued; because the calling routine has made the time series
	 * have an even length, the Nyquist component must also be real-valued;
	 * also this routine makes the assumption that both the DC and the
	 * Nyquist components are zero */

	work->data->data[0] = cabs(work->data->data[0]);
	work->data->data[work->data->length - 1] =
		creal(work->data->data[work->data->length - 1]);
	
	/* return data to time domain */

	if (XLALREAL8FreqTimeFFT(segment, work, revplan) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	return 0;
}


/* Helper routine that computes a segment of strain data with a single
 * time delay and beam pattern applied to the whole segment.  The duration
 * of the segment must therefore be reasonably short or else the movement
 * of the earth will invalidate the use of a single time shift and beam
 * pattern for the entire segment. */
static int XLALSimComputeLWLStrainSegmentREAL4TimeSeries(
	REAL4TimeSeries *segment,
	const REAL4TimeSeries *hplus,
	const REAL4TimeSeries *hcross,
	COMPLEX8FrequencySeries *work,
	REAL4FFTPlan *fwdplan,
	REAL4FFTPlan *revplan,
	REAL4Window *window,
	double ra,
	double dec,
	double psi,
	LALDetector *detector,
	const COMPLEX8FrequencySeries *response
)
{
	LIGOTimeGPS t;
	double gmst;
	double fplus;
	double fcross;
	double deltaT;
	double offint;
	double offrac;
	int offset;
	int j;
	size_t k;

	/* this routine assumes the segment has a length of N points where N is
	 * a power of two and that the workspace frequency series and the FFT
	 * plans are compatible with the size N; these assumptions are not
	 * checked: the calling routine must ensure that they are true */

	/* compute fplus, fcross, and time delay from earth's center at the
 	 * time corresponding to the middle of the segment */

	t = segment->epoch;
	XLALGPSAdd(&t, 0.5 * segment->data->length * segment->deltaT);
	gmst = XLALGreenwichMeanSiderealTime(&t);
	XLALComputeDetAMResponse(&fplus, &fcross,
		(const REAL4(*)[3])detector->response, ra, dec, psi, gmst);
	deltaT = XLALTimeDelayFromEarthCenter(detector->location, ra, dec, &t);

	/* add to the geometric delay the difference in time between the
	 * beginning of the injection timeseries and the beginning of the
	 * segment */

	deltaT += XLALGPSDiff(&hplus->epoch, &segment->epoch);

	/* compute the integer and fractional parts of the sample index in the
	 * segment on which the hplus and hcross time series begins: modf()
	 * returns integer and fractional parts that have the same sign, e.g.,
	 * -3.9 --> -3 + -0.9, and we adjust these so that magnitude of the
	 * fractional part is not greater than 0.5, e.g., -3.9 --> -4.0 + 0.1,
	 * so that we never do more than 1/2 a sample of re-interpolation */

	offrac = modf(deltaT / segment->deltaT, &offint);
	if (offrac < -0.5) {
		offrac += 1.0;
		offint -= 1.0;
	} else if (offrac > 0.5) {
		offrac -= 1.0;
		offint += 1.0;
	}
	offset = offint;

	/* now compute the sub-sample time shift that must be applied to the
	 * segment data */

	deltaT = offrac * segment->deltaT;

	/* compute the strain from hplus and hcross and populate the segment
	 * with this windowed data; apply appropriate integer offset */

	for (j = 0; j < (int)segment->data->length; ++j)
		if (j >= offset && j < (int)hplus->data->length + offset)
			segment->data->data[j] = window->data->data[j]
				* ( fplus * hplus->data->data[j - offset]
				+ fcross * hcross->data->data[j - offset] );
		else
			segment->data->data[j] = 0.0;


	/* put segment in frequency domain */

	if (XLALREAL4TimeFreqFFT(work, segment, fwdplan) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	/* apply sub-sample time shift in frequency domain */

	for (k = 0; k < work->data->length; ++k) {
		double f = work->f0 + k * work->deltaF;
		COMPLEX8 fac;
		
		/* phase for sub-sample time correction */
		fac = cexp(-I * LAL_TWOPI * f * deltaT);
		if (response)
			fac /= response->data->data[k];

		work->data->data[k] *= fac;
	}

	/* adjust DC and Nyquist components: the DC component must always be
	 * real-valued; because the calling routine has made the time series
	 * have an even length, the Nyquist component must also be real-valued;
	 * also this routine makes the assumption that both the DC and the
	 * Nyquist components are zero */

	work->data->data[0] = cabs(work->data->data[0]);
	work->data->data[work->data->length - 1] =
		creal(work->data->data[work->data->length - 1]);
	
	/* return data to time domain */

	if (XLALREAL4FreqTimeFFT(segment, work, revplan) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	return 0;
}


/**
 * @brief Computes strain for a detector and injects into target time series.
 * @details This routine takes care of the time-changing time delay from
 * the Earth's center and the time-changing antenna response pattern; it
 * does NOT account for deviations from the long-wavelength limit at high
 * frequencies.  An optional calibration response function can be provided
 * if the output time series is not in strain units.
 * @param[in,out] target Time series to inject strain into.
 * @param[in] hplus Time series with plus-polarization gravitational waveform.
 * @param[in] hcross Time series with cross-polarization gravitational waveform.
 * @param[in] ra Right ascension of the source (radians).
 * @param[in] dec Declination of the source (radians).
 * @param[in] psi Polarization angle of the source (radians).
 * @param[in] detector Detector to use when computing strain.
 * @param[in] response Response function to use, or NULL if none.
 * @retval 0 Success.
 * @retval <0 Failure.
 * @warning This routine assumes the long-wavelength limit (LWL) is valid
 * when computing the detector strain; at high frequencies near the free
 * spectral range of an interferometric detector this approximation becomes
 * invalid.
 */
int XLALSimInjectLWLDetectorStrainREAL8TimeSeries(
	REAL8TimeSeries *target,
	const REAL8TimeSeries *hplus,
	const REAL8TimeSeries *hcross,
	double ra,
	double dec,
	double psi,
	LALDetector *detector,
	const COMPLEX16FrequencySeries *response
)
{
	const double nominal_segdur = 2.0; /* nominal segment duration = 2s */
	const double max_time_delay = 0.1; /* generous allowed time delay */
	const size_t strides_per_segment = 2; /* 2 strides in one segment */
	LIGOTimeGPS t0;
	LIGOTimeGPS t1;
	size_t length;		/* length in samples of interval t0 - t1 */
	size_t seglen;		/* length of segment in samples */
	size_t padlen;		/* padding at beginning and end of segment */
	size_t ovrlap;		/* overlapping data length */
	size_t stride;		/* stride of each step */
	size_t nsteps;		/* number of steps to take */
	REAL8TimeSeries *h = NULL; /* strain timeseries to inject into target */
	REAL8TimeSeries *segment = NULL;
	COMPLEX16FrequencySeries *work = NULL;
	REAL8FFTPlan *fwdplan = NULL;
	REAL8FFTPlan *revplan = NULL;
	REAL8Window *window = NULL;
	size_t step;
	size_t j;
	int errnum = 0;

	/* check validity and compatibility of time series */

	LAL_CHECK_VALID_SERIES(target, XLAL_FAILURE);
	LAL_CHECK_VALID_SERIES(hplus, XLAL_FAILURE);
	LAL_CHECK_VALID_SERIES(hcross, XLAL_FAILURE);
	LAL_CHECK_CONSISTENT_TIME_SERIES(hplus, hcross, XLAL_FAILURE);
	if (response == NULL) {
		LAL_CHECK_COMPATIBLE_TIME_SERIES(target, hplus, XLAL_FAILURE);
	} else {
		/* units do no need to agree, but sample interval and
		 * start frequency do */
		if (fabs(target->deltaT - hplus->deltaT ) > LAL_REAL8_EPS)
			XLAL_ERROR(XLAL_ETIME);
		if (fabs(target->f0 - hplus->f0) > LAL_REAL8_EPS)
			XLAL_ERROR(XLAL_EFREQ);
	}
			

	/* constants describing the data segmentation: the length of the
	 * segment must be a power of two and the segment duration is at least
	 * nominal_segdur */

	seglen = round_up_to_power_of_two(nominal_segdur / target->deltaT);
	stride = seglen / strides_per_segment;
	padlen = max_time_delay / target->deltaT;
	ovrlap = seglen;
	ovrlap -= 2 * padlen;
	ovrlap -= stride;

	/* determine start and end time: the start time is the later of the
	 * start of the hplus/hcross time series and the target time series;
	 * the end time is the earlier of the end of the hplus/hcross time
	 * series and the target time series */

	t0 = hplus->epoch;
	t1 = target->epoch;
	XLALGPSAdd(&t0, hplus->data->length * hplus->deltaT);
	XLALGPSAdd(&t1, target->data->length * target->deltaT);
	t1 = XLALGPSCmp(&t1, &t0) < 0 ? t1 : t0;
	t0 = hplus->epoch;
	t0 = XLALGPSCmp(&t0, &target->epoch) > 0 ? t0 : target->epoch;

	/* add padding of 1 stride before and after these start and end times */

	XLALGPSAdd(&t0, -1.0 * stride * target->deltaT);
	XLALGPSAdd(&t1, stride * target->deltaT);

	/* determine if this is a disjoint set: if so, there is nothing to do */

	if (XLALGPSCmp(&t1, &t0) <= 0)
		return 0;

	/* create a segment that is seglen samples long */

	segment = XLALCreateREAL8TimeSeries(NULL, &t0, target->f0,
		target->deltaT, &target->sampleUnits, seglen);
	if (!segment) {
		errnum = XLAL_EFUNC;
		goto freereturn;
	}

	/* create a time series to hold the strain to inject into the target */

	length = XLALGPSDiff(&t1, &t0) / target->deltaT;
	h = XLALCreateREAL8TimeSeries(NULL, &t0, target->f0, target->deltaT,
		&target->sampleUnits, length);
	if (!h) {
		errnum = XLAL_EFUNC;
		goto freereturn;
	}
	memset(h->data->data, 0, h->data->length * sizeof(*h->data->data));

	/* determine number of steps it takes to go from t0 to t1 */

	nsteps = ((length%stride) ? (1 + length/stride) : (length/stride));


	/* create frequency-domain workspace; note that the FFT function
	 * populates the frequency series' metadata with the appropriate
	 * values */

	work = XLALCreateCOMPLEX16FrequencySeries(NULL, &t0, 0, 0,
		&lalDimensionlessUnit, seglen / 2 + 1);
	if (!work) {
		errnum = XLAL_EFUNC;
		goto freereturn;
	}

	/* create forward and reverse FFT plans */

	fwdplan = XLALCreateForwardREAL8FFTPlan(seglen, 0);
	revplan = XLALCreateReverseREAL8FFTPlan(seglen, 0);
	if (!fwdplan || !revplan) {
		errnum = XLAL_EFUNC;
		goto freereturn;
	}

	/* create a Tukey window with tapers entirely within the padding */

	window = XLALCreateTukeyREAL8Window(seglen, (double)padlen / seglen); 


	/* loop over steps, adding data from the current step to the strain */

	for (step = 0; step < nsteps; ++step) {

		int status;
		size_t offset;

		/* compute one segment of strain with time appropriate beam
 		 * pattern functions and time delays from earth's center */

		status = XLALSimComputeLWLStrainSegmentREAL8TimeSeries(segment,
			hplus, hcross, work, fwdplan, revplan, window, ra, dec,
			psi, detector, response);
		if (status < 0) {
			errnum = XLAL_EFUNC;
			goto freereturn;
		}

		/* compute the offset of this segment relative to the strain series */

		offset = XLALGPSDiff(&segment->epoch, &h->epoch) / h->deltaT;
		for (j = padlen; j < seglen - padlen; ++j)
			if ((j + offset) < h->data->length) {
				if (step && j - padlen < ovrlap) {
					/* feather overlapping data */
					double x = (double)(j - padlen) / ovrlap;
					h->data->data[j + offset] = x * segment->data->data[j]
						+ (1.0 - x) * h->data->data[j + offset];
				} else /* no feathering of remaining data */
					h->data->data[j + offset] = segment->data->data[j];
			}

		/* advance segment start time the next step */

		XLALGPSAdd(&segment->epoch, stride * segment->deltaT);
	}

	/* apply window to beginning and end of time series to reduce ringing */

	for (j = 0; j < stride - padlen; ++j)
		h->data->data[j] = h->data->data[h->data->length - 1 - j] = 0.0;
	for ( ; j < stride; ++j) {
		double fac = window->data->data[j - (stride - padlen)];
		h->data->data[j] *= fac;
		h->data->data[h->data->length - 1 - j] *= fac;
	}
	

	/* add computed strain to target time series */

	XLALAddREAL8TimeSeries(target, h);

freereturn:

	/* free all memory and return */

	XLALDestroyREAL8Window(window);
	XLALDestroyREAL8FFTPlan(revplan);
	XLALDestroyREAL8FFTPlan(fwdplan);
	XLALDestroyCOMPLEX16FrequencySeries(work);
	XLALDestroyREAL8TimeSeries(h);
	XLALDestroyREAL8TimeSeries(segment);

	if (errnum)
		XLAL_ERROR(errnum);
	return 0;
}


/**
 * @brief Computes strain for a detector and injects into target time series.
 * @details This routine takes care of the time-changing time delay from
 * the Earth's center and the time-changing antenna response pattern; it
 * does NOT account for deviations from the long-wavelength limit at high
 * frequencies.  An optional calibration response function can be provided
 * if the output time series is not in strain units.
 * @param[in,out] target Time series to inject strain into.
 * @param[in] hplus Time series with plus-polarization gravitational waveform.
 * @param[in] hcross Time series with cross-polarization gravitational waveform.
 * @param[in] ra Right ascension of the source (radians).
 * @param[in] dec Declination of the source (radians).
 * @param[in] psi Polarization angle of the source (radians).
 * @param[in] detector Detector to use when computing strain.
 * @param[in] response Response function to use, or NULL if none.
 * @retval 0 Success.
 * @retval <0 Failure.
 * @warning This routine assumes the long-wavelength limit (LWL) is valid
 * when computing the detector strain; at high frequencies near the free
 * spectral range of an interferometric detector this approximation becomes
 * invalid.
 */
int XLALSimInjectLWLDetectorStrainREAL4TimeSeries(
	REAL4TimeSeries *target,
	const REAL4TimeSeries *hplus,
	const REAL4TimeSeries *hcross,
	double ra,
	double dec,
	double psi,
	LALDetector *detector,
	const COMPLEX8FrequencySeries *response
)
{
	const double nominal_segdur = 2.0; /* nominal segment duration = 2s */
	const double max_time_delay = 0.1; /* generous allowed time delay */
	const size_t strides_per_segment = 2; /* 2 strides in one segment */
	LIGOTimeGPS t0;
	LIGOTimeGPS t1;
	size_t length;		/* length in samples of interval t0 - t1 */
	size_t seglen;		/* length of segment in samples */
	size_t padlen;		/* padding at beginning and end of segment */
	size_t ovrlap;		/* overlapping data length */
	size_t stride;		/* stride of each step */
	size_t nsteps;		/* number of steps to take */
	REAL4TimeSeries *h = NULL; /* strain timeseries to inject into target */
	REAL4TimeSeries *segment = NULL;
	COMPLEX8FrequencySeries *work = NULL;
	REAL4FFTPlan *fwdplan = NULL;
	REAL4FFTPlan *revplan = NULL;
	REAL4Window *window = NULL;
	size_t step;
	size_t j;
	int errnum = 0;

	/* check validity and compatibility of time series */

	LAL_CHECK_VALID_SERIES(target, XLAL_FAILURE);
	LAL_CHECK_VALID_SERIES(hplus, XLAL_FAILURE);
	LAL_CHECK_VALID_SERIES(hcross, XLAL_FAILURE);
	LAL_CHECK_CONSISTENT_TIME_SERIES(hplus, hcross, XLAL_FAILURE);
	if (response == NULL) {
		LAL_CHECK_COMPATIBLE_TIME_SERIES(target, hplus, XLAL_FAILURE);
	} else {
		/* units do no need to agree, but sample interval and
		 * start frequency do */
		if (fabs(target->deltaT - hplus->deltaT ) > LAL_REAL8_EPS)
			XLAL_ERROR(XLAL_ETIME);
		if (fabs(target->f0 - hplus->f0) > LAL_REAL8_EPS)
			XLAL_ERROR(XLAL_EFREQ);
	}
			

	/* constants describing the data segmentation: the length of the
	 * segment must be a power of two and the segment duration is at least
	 * nominal_segdur */

	seglen = round_up_to_power_of_two(nominal_segdur / target->deltaT);
	stride = seglen / strides_per_segment;
	padlen = max_time_delay / target->deltaT;
	ovrlap = seglen;
	ovrlap -= 2 * padlen;
	ovrlap -= stride;

	/* determine start and end time: the start time is the later of the
	 * start of the hplus/hcross time series and the target time series;
	 * the end time is the earlier of the end of the hplus/hcross time
	 * series and the target time series */

	t0 = hplus->epoch;
	t1 = target->epoch;
	XLALGPSAdd(&t0, hplus->data->length * hplus->deltaT);
	XLALGPSAdd(&t1, target->data->length * target->deltaT);
	t1 = XLALGPSCmp(&t1, &t0) < 0 ? t1 : t0;
	t0 = hplus->epoch;
	t0 = XLALGPSCmp(&t0, &target->epoch) > 0 ? t0 : target->epoch;

	/* add padding of 1 stride before and after these start and end times */

	XLALGPSAdd(&t0, -1.0 * stride * target->deltaT);
	XLALGPSAdd(&t1, stride * target->deltaT);

	/* determine if this is a disjoint set: if so, there is nothing to do */

	if (XLALGPSCmp(&t1, &t0) <= 0)
		return 0;

	/* create a segment that is seglen samples long */

	segment = XLALCreateREAL4TimeSeries(NULL, &t0, target->f0,
		target->deltaT, &target->sampleUnits, seglen);
	if (!segment) {
		errnum = XLAL_EFUNC;
		goto freereturn;
	}

	/* create a time series to hold the strain to inject into the target */

	length = XLALGPSDiff(&t1, &t0) / target->deltaT;
	h = XLALCreateREAL4TimeSeries(NULL, &t0, target->f0, target->deltaT,
		&target->sampleUnits, length);
	if (!h) {
		errnum = XLAL_EFUNC;
		goto freereturn;
	}
	memset(h->data->data, 0, h->data->length * sizeof(*h->data->data));

	/* determine number of steps it takes to go from t0 to t1 */

	nsteps = ((length%stride) ? (1 + length/stride) : (length/stride));


	/* create frequency-domain workspace; note that the FFT function
	 * populates the frequency series' metadata with the appropriate
	 * values */

	work = XLALCreateCOMPLEX8FrequencySeries(NULL, &t0, 0, 0,
		&lalDimensionlessUnit, seglen / 2 + 1);
	if (!work) {
		errnum = XLAL_EFUNC;
		goto freereturn;
	}

	/* create forward and reverse FFT plans */

	fwdplan = XLALCreateForwardREAL4FFTPlan(seglen, 0);
	revplan = XLALCreateReverseREAL4FFTPlan(seglen, 0);
	if (!fwdplan || !revplan) {
		errnum = XLAL_EFUNC;
		goto freereturn;
	}

	/* create a Tukey window with tapers entirely within the padding */

	window = XLALCreateTukeyREAL4Window(seglen, (double)padlen / seglen); 


	/* loop over steps, adding data from the current step to the strain */

	for (step = 0; step < nsteps; ++step) {

		int status;
		size_t offset;

		/* compute one segment of strain with time appropriate beam
 		 * pattern functions and time delays from earth's center */

		status = XLALSimComputeLWLStrainSegmentREAL4TimeSeries(segment,
			hplus, hcross, work, fwdplan, revplan, window, ra, dec,
			psi, detector, response);
		if (status < 0) {
			errnum = XLAL_EFUNC;
			goto freereturn;
		}

		/* compute the offset of this segment relative to the strain series */

		offset = XLALGPSDiff(&segment->epoch, &h->epoch) / h->deltaT;
		for (j = padlen; j < seglen - padlen; ++j)
			if ((j + offset) < h->data->length) {
				if (step && j - padlen < ovrlap) {
					/* feather overlapping data */
					double x = (double)(j - padlen) / ovrlap;
					h->data->data[j + offset] = x * segment->data->data[j]
						+ (1.0 - x) * h->data->data[j + offset];
				} else /* no feathering of remaining data */
					h->data->data[j + offset] = segment->data->data[j];
			}

		/* advance segment start time the next step */

		XLALGPSAdd(&segment->epoch, stride * segment->deltaT);
	}

	/* apply window to beginning and end of time series to reduce ringing */

	for (j = 0; j < stride - padlen; ++j)
		h->data->data[j] = h->data->data[h->data->length - 1 - j] = 0.0;
	for ( ; j < stride; ++j) {
		double fac = window->data->data[j - (stride - padlen)];
		h->data->data[j] *= fac;
		h->data->data[h->data->length - 1 - j] *= fac;
	}
	

	/* add computed strain to target time series */

	XLALAddREAL4TimeSeries(target, h);

freereturn:

	/* free all memory and return */

	XLALDestroyREAL4Window(window);
	XLALDestroyREAL4FFTPlan(revplan);
	XLALDestroyREAL4FFTPlan(fwdplan);
	XLALDestroyCOMPLEX8FrequencySeries(work);
	XLALDestroyREAL4TimeSeries(h);
	XLALDestroyREAL4TimeSeries(segment);

	if (errnum)
		XLAL_ERROR(errnum);
	return 0;
}
