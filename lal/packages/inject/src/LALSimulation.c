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
#include <lal/LALSimulation.h>
#include <lal/LALDetectors.h>
#include <lal/LALComplex.h>
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


#include <lal/LALRCSID.h>
NRCSID(LALSIMULATIONC, "$Id:");


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
 * Turn an instrument name into a LALDetector structure.  The first two
 * characters of the input string are used as the instrument name, which
 * allows channel names in the form "H1:LSC-STRAIN" to be used.  The return
 * value is a pointer into the lalCachedDetectors array, so modifications
 * to the contents are global.  Make a copy of the structure if you want to
 * modify it safely.
 */


const LALDetector *XLALInstrumentNameToLALDetector(
	const char *string
)
{
	static const char func[] = "XLALInstrumentNameToLALDetector";
	int i;

	for(i = 0; i < LAL_NUM_DETECTORS; i++)
		if(!strncmp(string, lalCachedDetectors[i].frDetector.prefix, 2))
			return &lalCachedDetectors[i];

	XLALPrintError("%s(): error: can't identify instrument from string \"%s\"\n", func, string);
	XLAL_ERROR_NULL(func, XLAL_EDATA);
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
 *   start of those time series at the geocentre (for simplicity the epochs
 *   must be the same),
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
	LALDetector *detector
)
{
	static const char func[] = "XLALSimDetectorStrainREAL8TimeSeries";
	char name[13];	/* "?? injection" + terminator */
	REAL8TimeSeries *h;
	unsigned i;

	LAL_CHECK_VALID_SERIES(hplus, NULL);
	LAL_CHECK_VALID_SERIES(hcross, NULL);
	LAL_CHECK_CONSISTENT_TIME_SERIES(hplus, hcross, NULL);

	/* generate name */

	sprintf(name, "%2s injection", detector->frDetector.prefix);

	/* allocate output time series. */

	h = XLALCreateREAL8TimeSeries(name, &hplus->epoch, hplus->f0, hplus->deltaT, &hplus->sampleUnits, hplus->data->length);
	if(!h)
		XLAL_ERROR_NULL(func, XLAL_EFUNC);

	/* add the detector's geometric delay.  after this, epoch = the
	 * time of the injection time series' first sample at the desired
	 * detector */

	XLALGPSAdd(&h->epoch, XLALTimeDelayFromEarthCenter(detector->location, right_ascension, declination, &h->epoch));

	/* project + and x time series onto detector */

	for(i = 0; i < h->data->length; i++) {
		LIGOTimeGPS t = h->epoch;
		double fplus, fcross;

		XLALGPSAdd(&t, i * h->deltaT);
		XLALComputeDetAMResponse(&fplus, &fcross, detector->response, right_ascension, declination, psi, XLALGreenwichMeanSiderealTime(&t));

		h->data->data[i] = fplus * hplus->data->data[i] + fcross * hcross->data->data[i];
	}

	/* done */

	return h;
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
	static const char func[] = "XLALSimAddInjectionREAL8TimeSeries";
	COMPLEX16FrequencySeries *tilde_h;
	REAL8FFTPlan *plan;
	REAL8Window *window;
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
		XLALPrintError("%s(): error: input sample rates or heterodyne frequencies do not match\n", func);
		XLAL_ERROR(func, XLAL_EINVAL);
	}

	/* extend the source time series by adding the "aperiodicity
	 * padding" to the start and end.  for efficiency's sake, make sure
	 * the new length is a power of two. */

	i = round_up_to_power_of_two(h->data->length + 2 * aperiodicity_suppression_buffer);
	if(i < h->data->length) {
		/* integer overflow */
		XLALPrintError("%s(): error: source time series too long\n", func);
		XLAL_ERROR(func, XLAL_EBADLEN);
	}
	i -= h->data->length;
	if(!XLALResizeREAL8TimeSeries(h, -(int) (i / 2), h->data->length + i))
		XLAL_ERROR(func, XLAL_EFUNC);

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

	/* transform source time series to frequency domain.  the FFT
	 * function populates the frequency series' metadata with the
	 * appropriate values. */

	tilde_h = XLALCreateCOMPLEX16FrequencySeries(NULL, &h->epoch, 0, 0, &lalDimensionlessUnit, h->data->length / 2 + 1);
	plan = XLALCreateForwardREAL8FFTPlan(h->data->length, 0);
	if(!tilde_h || !plan) {
		XLALDestroyCOMPLEX16FrequencySeries(tilde_h);
		XLALDestroyREAL8FFTPlan(plan);
		XLAL_ERROR(func, XLAL_EFUNC);
	}
	i = XLALREAL8TimeFreqFFT(tilde_h, h, plan);
	XLALDestroyREAL8FFTPlan(plan);
	if(i) {
		XLALDestroyCOMPLEX16FrequencySeries(tilde_h);
		XLAL_ERROR(func, XLAL_EFUNC);
	}

	/* apply sub-sample time correction and optional response function
	 * */

	for(i = 0; i < tilde_h->data->length; i++) {
		const double f = tilde_h->f0 + i * tilde_h->deltaF;
		COMPLEX16 fac;

		/* phase for sub-sample time correction */

		fac = LAL_CEXP(LAL_CMUL_REAL(LAL_COMPLEX16_I, -LAL_TWOPI * f * start_sample_frac * target->deltaT));

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
			if(LAL_COMPLEX_EQ(response->data->data[j], LAL_COMPLEX16_ZERO))
				fac = LAL_COMPLEX16_ZERO;
			else
				fac = LAL_CDIV(fac, response->data->data[j]);
		}

		/* apply factor */

		tilde_h->data->data[i] = LAL_CMUL(tilde_h->data->data[i], fac);
	}

	/* adjust DC and Nyquist components.  the DC component must always
	 * be real-valued.  because we have adjusted the source time series
	 * to have a length that is an even integer (we've made it a power
	 * of 2) the Nyquist component must also be real valued. */

	if(response) {
		/* a response function has been provided.  zero the DC and
		 * Nyquist components */
		if(tilde_h->f0 == 0.0)
			tilde_h->data->data[0] = LAL_COMPLEX16_ZERO;
		tilde_h->data->data[tilde_h->data->length - 1] = LAL_COMPLEX16_ZERO;
	} else {
		/* no response has been provided.  set the phase of the DC
		 * component to 0, set the imaginary component of the
		 * Nyquist to 0 */
		if(tilde_h->f0 == 0.0)
			tilde_h->data->data[0] = XLALCOMPLEX16Rect(LAL_CABS(tilde_h->data->data[0]), 0.0);
		tilde_h->data->data[tilde_h->data->length - 1] = XLALCOMPLEX16Rect(LAL_REAL(tilde_h->data->data[tilde_h->data->length - 1]), 0.0);
	}

	/* return to time domain */

	plan = XLALCreateReverseREAL8FFTPlan(h->data->length, 0);
	if(!plan) {
		XLALDestroyCOMPLEX16FrequencySeries(tilde_h);
		XLAL_ERROR(func, XLAL_EFUNC);
	}
	i = XLALREAL8FreqTimeFFT(h, tilde_h, plan);
	XLALDestroyREAL8FFTPlan(plan);
	XLALDestroyCOMPLEX16FrequencySeries(tilde_h);
	if(i)
		XLAL_ERROR(func, XLAL_EFUNC);

	/* the deltaT can get "corrupted" by floating point round-off
	 * during its trip through the frequency domain.  since this
	 * function starts by confirming that the sample rate of the source
	 * matches that of the target time series, we can use the target
	 * series' sample rate to reset the source's sample rate to its
	 * original value.  but we do a check to make sure we're not
	 * masking a real bug */

	if(fabs(h->deltaT - target->deltaT) / target->deltaT > 1e-12) {
		XLALPrintError("%s(): error: oops, internal sample rate mismatch\n", func);
		XLAL_ERROR(func, XLAL_EERR);
	}
	h->deltaT = target->deltaT;

	/* set source epoch from target epoch and integer sample offset */

	h->epoch = target->epoch;
	XLALGPSAdd(&h->epoch, start_sample_int * target->deltaT);

	/* clip half of the "aperiodicity padding" from the start and end
	 * of the source time series in a continuing effort to suppress
	 * aperiodicity artifacts. */

	if(!XLALResizeREAL8TimeSeries(h, aperiodicity_suppression_buffer / 2, h->data->length - aperiodicity_suppression_buffer))
		XLAL_ERROR(func, XLAL_EFUNC);

	/* apply a Tukey window whose tapers lie within the remaining
	 * aperiodicity padding. leaving one sample of the aperiodicty
	 * padding untouched on each side of the original time series
	 * because the data might have been shifted into it */

	window = XLALCreateTukeyREAL8Window(h->data->length, (double) (aperiodicity_suppression_buffer - 2) / h->data->length);
	if(!window)
		XLAL_ERROR(func, XLAL_EFUNC);
	for(i = 0; i < h->data->length; i++)
		h->data->data[i] *= window->data->data[i];
	XLALDestroyREAL8Window(window);

	/* add source time series to target time series */

	if(!XLALAddREAL8TimeSeries(target, h))
		XLAL_ERROR(func, XLAL_EFUNC);

	/* done */

	return 0;
}
