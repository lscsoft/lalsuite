#include <math.h>
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
#include "inject.h"

REAL8TimeSeries * XLALSignalDetectorStrainREAL8TimeSeries( REAL8TimeSeries *hplus, REAL8TimeSeries *hcross, SkyPosition *position, REAL8 psi, LALDetector *detector )
{
	static const char *func = "XLALDetectorStrain";
	REAL8TimeSeries *h;
	LIGOTimeGPS epoch;
	SkyPosition equatorial;
	REAL8 fplus;
	REAL8 fcross;
	REAL8 gmst;
	REAL8 delay;
	UINT4 j;

	/* FIXME: sanity check on parameters */

	/* FIXME: Convert any system to equatorial */
	if ( position->system != COORDINATESYSTEM_EQUATORIAL )
		XLAL_ERROR_NULL( func, XLAL_ETYPE );
	equatorial = *position;
	
	epoch = hplus->epoch;
	gmst = XLALGreenwichMeanSiderealTime( &epoch );

	XLALComputeDetAMResponse( &fplus, &fcross, detector->response, equatorial.longitude, equatorial.latitude, psi, gmst );
	delay = XLALTimeDelayFromEarthCenter( detector->location, equatorial.longitude, equatorial.latitude, &epoch );
	XLALGPSAdd( &epoch, delay );

	/* TODO: Encode detector in name */
	h = XLALCreateREAL8TimeSeries( "DETECTOR_STRAIN", &epoch, hplus->f0, hplus->deltaT, &lalStrainUnit, hplus->data->length );
	if ( ! h )
		XLAL_ERROR_NULL( func, XLAL_EFUNC );
	for ( j = 0; j < hplus->data->length; ++j )
		h->data->data[j] = fplus * hplus->data->data[j] + fcross * hcross->data->data[j];

	return h;
}

REAL8TimeSeries * XLALInjectionREAL8TimeSeries( REAL8TimeSeries *h, LIGOTimeGPS *start, REAL8 deltaT, UINT4 length, COMPLEX8FrequencySeries *response )
{
	static const char *func = "XLALInjectionREAL8TimeSeries";

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
	newlength = pow( 2, ceil( log(newlength)/log(2) ) );
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
		fac = XLALCOMPLEX16MulReal( LAL_COMPLEX16_I, phi );
		fac = XLALCOMPLEX16Exp( fac );
		if ( response ) {
			COMPLEX16 r16;
			COMPLEX8  r8;
			UINT4 k2;
			k2 = f/response->deltaF;
			if (k2 > response->data->length)
				k2 = response->data->length;
			r8 = response->data->data[k2];
		       	LAL_SET_COMPLEX( &r16, LAL_REAL(r8), LAL_IMAG(r8) );
			fac = XLALCOMPLEX16Div( fac, r16 );
		}
		injtilde->data->data[k] = XLALCOMPLEX16Mul( injtilde->data->data[k], fac );
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
