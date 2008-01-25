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
#include "check_series_macros.h"


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

#define KIPPSPROPOSAL
#ifndef KIPPSPROPOSAL

REAL8TimeSeries * XLALSimDetectorStrainREAL8TimeSeries( REAL8TimeSeries *hplus, REAL8TimeSeries *hcross, REAL8 right_ascension, REAL8 declination, REAL8 psi, LALDetector *detector )
{
	static const char *func = "XLALSimDetectorStrainREAL8TimeSeries";
	REAL8TimeSeries *h;
	LIGOTimeGPS epoch;
	REAL8 fplus;
	REAL8 fcross;
	REAL8 gmst;
	REAL8 delay;
	UINT4 j;

	LAL_CHECK_VALID_SERIES(hplus, NULL);
	LAL_CHECK_VALID_SERIES(hcross, NULL);
	LAL_CHECK_CONSISTENT_TIME_SERIES(hplus, hcross, NULL);

	/* FIXME: sanity check on parameters */

	epoch = hplus->epoch;
	gmst = XLALGreenwichMeanSiderealTime( &epoch );

	XLALComputeDetAMResponse( &fplus, &fcross, detector->response, right_ascension, declination, psi, gmst );
	delay = XLALTimeDelayFromEarthCenter( detector->location, right_ascension, declination, &epoch );
	XLALGPSAdd( &epoch, delay );

	/* TODO: Encode detector in name */
	h = XLALCreateREAL8TimeSeries( "DETECTOR_STRAIN", &epoch, hplus->f0, hplus->deltaT, &lalStrainUnit, hplus->data->length );
	if ( ! h )
		XLAL_ERROR_NULL( func, XLAL_EFUNC );
	for ( j = 0; j < hplus->data->length; ++j )
		h->data->data[j] = fplus * hplus->data->data[j] + fcross * hcross->data->data[j];

	return h;
}


REAL8TimeSeries * XLALSimInjectionREAL8TimeSeries( REAL8TimeSeries *h, LIGOTimeGPS *start, REAL8 deltaT, UINT4 length, COMPLEX16FrequencySeries *response )
{
	static const char *func = "XLALSimInjectionREAL8TimeSeries";

	REAL8TimeSeries *injection;
	COMPLEX16FrequencySeries *injtilde;
	REAL8FFTPlan *plan;
	REAL8 injdur = length * deltaT;
	REAL8 sigdur = h->data->length * h->deltaT;
	REAL8 dt;
	UINT4 newlength;
	UINT4 j;
	UINT4 k;

	/* cannot presently support different sample rates */
	if ( deltaT != h->deltaT )
		XLAL_ERROR_NULL( func, XLAL_EINVAL );

	injection = XLALCreateREAL8TimeSeries( "INJECTION", start, h->f0, deltaT, &lalStrainUnit, length );
	memset( injection->data->data, 0, injection->data->length*sizeof(*injection->data->data) );

	dt = XLALGPSDiff( &h->epoch, start );
	/* if signal starts after the injection period or
	 * if the signal ends before the injection period
	 * just return an empty injection */
	if ( dt > injdur || dt < -sigdur ) {
		if ( response )
			XLALUnitMultiply( &injection->sampleUnits, &injection->sampleUnits, &response->sampleUnits );
		return injection;
	}

	/* minimum duration is twice signal duration plus a padding
	 * of one second to take care of possible wrap-around corruption */
	newlength = length;
	if ( injdur < sigdur )
		newlength = (2.0*sigdur + 1.0)/deltaT;
	else
		newlength = (injdur + sigdur + 1.0)/deltaT;
	/* for efficiency sake, make sure newlength is a power of two */
	/* newlength = pow( 2, ceil( log(newlength)/log(2) ) ); */
	if ( ((newlength - 1) & newlength) ) { /* if statment not necessary */
		/* note: this doesn't work for values greater than 2147483648 */
		--newlength;
		newlength |= (newlength >> 1);
		newlength |= (newlength >> 2);
		newlength |= (newlength >> 4);
		newlength |= (newlength >> 8);
		newlength |= (newlength >> 16);
		++newlength;
	}

	XLALResizeREAL8TimeSeries( injection, 0, newlength );

	/* copy signal to beginning of injection */
	for ( j = 0; j < h->data->length; ++j )
		injection->data->data[j] = h->data->data[j];

	/* put injection into frequency domain and apply delay correction
	 * and response function, if necessary */

	injtilde = XLALCreateCOMPLEX16FrequencySeries( "INJTILDE", start, 0.0, 1.0/(newlength * deltaT), &lalDimensionlessUnit, newlength/2 + 1 );

	plan = XLALCreateForwardREAL8FFTPlan( newlength, 0 );
	XLALREAL8TimeFreqFFT( injtilde, injection, plan );
	XLALDestroyREAL8FFTPlan( plan );


	/* eliminate DC and Nyquist components and apply phase correction for
	 * time delay as well as response function */
	injtilde->data->data[0] = LAL_COMPLEX16_ZERO;
	injtilde->data->data[injtilde->data->length-1] = LAL_COMPLEX16_ZERO;
	for ( k = 1; k < injtilde->data->length - 1; ++k ) {
		COMPLEX16 fac;
		REAL8 f = k * injtilde->deltaF;
		REAL8 phi = -LAL_TWOPI * f * dt;
		fac = LAL_CMUL_REAL( LAL_COMPLEX16_I, phi );
		fac = LAL_CEXP( fac );
		if ( response ) {
			UINT4 k2 = floor( 0.5 + f / response->deltaF );
			if ( k2 > response->data->length )
				k2 = response->data->length;
			if ( LAL_COMPLEX_EQ( response->data->data[k2], LAL_COMPLEX16_ZERO ) )
				fac = LAL_COMPLEX16_ZERO;
			else
				fac = LAL_CDIV( fac, response->data->data[k2] );
		}
		injtilde->data->data[k] = LAL_CMUL( injtilde->data->data[k], fac );
	}

	plan = XLALCreateReverseREAL8FFTPlan( newlength, 0 );
	XLALREAL8FreqTimeFFT( injection, injtilde, plan );
	XLALDestroyREAL8FFTPlan( plan );

	XLALDestroyCOMPLEX16FrequencySeries( injtilde );

	XLALResizeREAL8TimeSeries( injection, 0, length );

	/* TODO: for some reason, response doesn't have correct units.
	if ( response )
		XLALUnitDivide( &injection->sampleUnits, &injection->sampleUnits, &response->sampleUnits );
	*/

	return injection;
}


#else	/* KIPPSPROPOSAL */


/*
 * Input
 *
 * - h+ and hx time series for the injection with t = 0 interpreted to be
 *   the "time" of the injection,
 *
 * - the right ascension and declination of the source,
 *
 * - the orientation of the wave co-ordinate system,
 *
 * - the detector into which the injection is destined to be injected,
 *
 * - and the "time" of the injection as observed at the geocentre.
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
 */


REAL8TimeSeries *XLALSimDetectorStrainREAL8TimeSeries(const REAL8TimeSeries *hplus, const REAL8TimeSeries *hcross, REAL8 right_ascension, REAL8 declination, REAL8 psi, LALDetector *detector, const LIGOTimeGPS *injection_time_at_geocentre)
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

	/* allocate output time series.  epoch = injection "time" at
	 * geocentre */

	h = XLALCreateREAL8TimeSeries(name, injection_time_at_geocentre, hplus->f0, hplus->deltaT, &lalStrainUnit, hplus->data->length);
	if(!h)
		XLAL_ERROR_NULL(func, XLAL_EFUNC);

	/* add start time of input time series and the detector's geometric
	 * delay.  after this, epoch = the time of the injection time
	 * series' first sample at the desired detector */

	XLALGPSAdd(&h->epoch, XLALGPSGetREAL8(&hplus->epoch) + XLALTimeDelayFromEarthCenter(detector->location, right_ascension, declination, injection_time_at_geocentre));

	/* project + and x time series onto detector */

	for(i = 0; i < h->data->length; i++) {
		LIGOTimeGPS epoch = h->epoch;
		double fplus, fcross;

		XLALGPSAdd(&epoch, i * h->deltaT);
		XLALComputeDetAMResponse(&fplus, &fcross, detector->response, right_ascension, declination, psi, XLALGreenwichMeanSiderealTime(&epoch));

		h->data->data[i] = fplus * hplus->data->data[i] + fcross * hcross->data->data[i];
	}

	/* done */

	return h;
}


/*
 * Essentially a wrapper for XLALAddREAL8TimeSeries(), but performs
 * sub-sample re-interpolation to adjust the source time series epoch to
 * lie on an integer sample boundary in the target time series.  This
 * transformation is done in the frequency domain, so it is convenient to
 * allow a response function to be applied at the same time.  Note that the
 * source time series is modified in place by this function.
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


int XLALAddInjectionREAL8TimeSeries(REAL8TimeSeries *target, REAL8TimeSeries *h, const COMPLEX16FrequencySeries *response)
{
	static const char func[] = "XLALAddInjectionREAL8TimeSeries";
	COMPLEX16FrequencySeries *tilde_h;
	REAL8FFTPlan *plan;
	unsigned i;
	long start_sample_int;
	double start_sample_frac;

	/* check input */

	if(h->deltaT != target->deltaT || h->f0 != target->f0)
		XLAL_ERROR(func, XLAL_EINVAL);

	/* extend the injection time series by (at least) 1 second of 0s in
	 * both directions with the hope that this suppresses
	 * (a)periodicity artifacts sufficiently.  computing count of
	 * samples first ensures that the same number is added to the start
	 * and end.  for efficiency's sake, make sure the new length is a
	 * power of two */

	i = 1.0 / h->deltaT;
	i = round_up_to_power_of_two(h->data->length + 2 * i);
	if(!XLALResizeREAL8TimeSeries(h, -(int) (i / 2), i))
		XLAL_ERROR(func, XLAL_EFUNC);

	/* compute the integer and fractional parts of the sample in the
	 * target time series on which the injection time series begins.
	 * the fractional part is always positive, e.g. -3.5 --> -4 + 0.5.
	 * the integer part will be used to set a new epoch and the
	 * fractional part used to re-interpolate the injection */

	start_sample_frac = XLALGPSDiff(&h->epoch, &target->epoch) / target->deltaT;
	start_sample_int = floor(start_sample_frac);
	start_sample_frac -= start_sample_int;

	/* transform injection to frequency domain.  the FFT function
	 * populates the frequency series' metadata with the appropirate
	 * values. */

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

	/* apply delay correction and optional response function */

	for(i = 0; i < tilde_h->data->length; i++) {
		const double f = tilde_h->f0 + i * tilde_h->deltaF;
		COMPLEX16 fac;

		/* phase for time delay */

		fac = LAL_CEXP(LAL_CMUL_REAL(LAL_COMPLEX16_I, -LAL_TWOPI * f * start_sample_frac * target->deltaT));

		/* response function */

		if(response) {
			unsigned j = floor(f / response->deltaF + 0.5);
			if(j > response->data->length)
				j = response->data->length;
			if(LAL_COMPLEX_EQ(response->data->data[j], LAL_COMPLEX16_ZERO))
				fac = LAL_COMPLEX16_ZERO;
			else
				fac = LAL_CDIV(fac, response->data->data[j]);
		}

		/* apply factor */

		tilde_h->data->data[i] = LAL_CMUL(tilde_h->data->data[i], fac);
	}

	/* zero DC and Nyquist components */

	if(tilde_h->f0 == 0)
		tilde_h->data->data[0] = LAL_COMPLEX16_ZERO;
	tilde_h->data->data[tilde_h->data->length - 1] = LAL_COMPLEX16_ZERO;

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

	/* the deltaT can get corrupted by floating point round-off during
	 * its trip through the frequency domain.  since this function
	 * starts by confirming that the sample rate of the injection
	 * matches that of the target time series, we can use the target
	 * series' sample rate to reset the injection's sample rate to its
	 * original value.  but do a check to make sure we're not masking a
	 * real bug */

	if(fabs(h->deltaT - target->deltaT) / h->deltaT > 1e-12) {
		XLALPrintError("sample rate mismatch");
		XLAL_ERROR(func, XLAL_EERR);
	}
	h->deltaT = target->deltaT;

	/* set epoch from integer sample offset */

	h->epoch = target->epoch;
	XLALGPSAdd(&h->epoch, start_sample_int * h->deltaT);

	/* add to target time series */

	if(!XLALAddREAL8TimeSeries(target, h))
		XLAL_ERROR(func, XLAL_EFUNC);

	/* done */

	return 0;
}

#endif	/* KIPPSPROPOSAL */
