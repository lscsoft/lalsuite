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

	/*
	 * TODO: for some reason, response doesn't have correct units.
	if ( response )
		XLALUnitDivide( &injection->sampleUnits, &injection->sampleUnits, &response->sampleUnits );
		*/

	return injection;
}
