/*
*  Copyright (C) 2007 David Chin, Jolien Creighton, Kipp Cannon, Peter Shawhan
*  Copyright (C) 2012 Matthew Pitkin
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
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALError.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/DetResponse.h>
#include <lal/TimeSeries.h>
#include <lal/XLALError.h>

/**
 * An implementation of the detector response formulae in Anderson et al
 * PRD 63 042003 (2001) \cite ABCF2001.
 *
 * Computes F+ and Fx for a source at a specified sky position,
 * polarization angle, and sidereal time.  Also requires the detector's
 * response matrix which is defined by Eq. (B6) of [ABCF] using either
 * Table 1 of \cite ABCF2001 or Eqs. (B11)--(B17) to compute the arm
 * direction unit vectors.
 */
void XLALComputeDetAMResponse(
	double *fplus,		/**< Returned value of F+ */
	double *fcross,		/**< Returned value of Fx */
	const REAL4 D[3][3],	/**< Detector response 3x3 matrix */
	const double ra,	/**< Right ascention of source (radians) */
	const double dec,	/**< Declination of source (radians) */
	const double psi,	/**< Polarization angle of source (radians) */
	const double gmst	/**< Greenwich mean sidereal time (radians) */
)
{
	int i;
	double X[3];
	double Y[3];

	/* Greenwich hour angle of source (radians). */
	const double gha = gmst - ra;

	/* pre-compute trig functions */
	const double cosgha = cos(gha);
	const double singha = sin(gha);
	const double cosdec = cos(dec);
	const double sindec = sin(dec);
	const double cospsi = cos(psi);
	const double sinpsi = sin(psi);

	/* Eq. (B4) of [ABCF].  Note that dec = pi/2 - theta, and gha =
	 * -phi where theta and phi are the standard spherical coordinates
	 * used in that paper. */
	X[0] = -cospsi * singha - sinpsi * cosgha * sindec;
	X[1] = -cospsi * cosgha + sinpsi * singha * sindec;
	X[2] =  sinpsi * cosdec;

	/* Eq. (B5) of [ABCF].  Note that dec = pi/2 - theta, and gha =
	 * -phi where theta and phi are the standard spherical coordinates
	 * used in that paper. */
	Y[0] =  sinpsi * singha - cospsi * cosgha * sindec;
	Y[1] =  sinpsi * cosgha + cospsi * singha * sindec;
	Y[2] =  cospsi * cosdec;

	/* Now compute Eq. (B7) of [ABCF] for each polarization state, i.e.,
	 * with s+=1 and sx=0 to get F+, with s+=0 and sx=1 to get Fx */
	*fplus = *fcross = 0.0;
	for(i = 0; i < 3; i++) {
		const double DX = D[i][0] * X[0] + D[i][1] * X[1] + D[i][2] * X[2];
		const double DY = D[i][0] * Y[0] + D[i][1] * Y[1] + D[i][2] * Y[2];
		*fplus  += X[i] * DX - Y[i] * DY;
		*fcross += X[i] * DY + Y[i] * DX;
	}
}


/**
 *
 * An implementation of the detector response for all six tensor, vector and
 * scalar polarisation modes of general metric theories of gravity. We follow
 * the convention of \cite Blaut2012 (this is also equivalent to the
 * Equations in \cite ABCF2001 and  \cite Nishizawa2009 in albeit with a
 * rotated set of coordinates), but with \cite Blaut2012's \f$\theta = \pi/2 -
 * dec\f$ and \f$\phi = ra-GMST\f$ rather than the gha = gmst - ra used here.
 *
 * The values computed are the tensor mode response, F+ and Fx ("cross"), the
 * scalar breathing and longitudinal modes, Fb and Fl, and the vector "x" and
 * "y" modes, Fx ("x") and Fy, for a source at a specified sky position,
 * polarization angle, and sidereal time.  Also requires the detector's
 * response matrix which is defined by Eq. (B6) of [ABCF] using either
 * Table 1 of \cite ABCF2001 or Eqs. (B11)--(B17) to compute the arm
 * direction unit vectors.
 */
void XLALComputeDetAMResponseExtraModes(
	double *fplus,		/**< Returned value of F+ */
	double *fcross,		/**< Returned value of Fx (cross) */
	double *fb,		/**< Returned value of Fb (breathing mode) */
	double *fl,		/**< Returned value of Fl (scalar longitudinal */
	double *fx,		/**< Returned value of Fx ("x" vector mode) */
	double *fy,		/**< Returned value of Fy (y vector mode) */
	REAL4 D[3][3],		/**< Detector response 3x3 matrix */
	const double ra,	/**< Right ascention of source (radians) */
	const double dec,	/**< Declination of source (radians) */
	const double psi,	/**< Polarization angle of source (radians) */
	const double gmst	/**< Greenwich mean sidereal time (radians) */
)
{
	int i;
	double X[3];
	double Y[3];
	double Z[3];

	/* Greenwich hour angle of source (radians). */
	const double gha = gmst - ra;

	/* pre-compute trig functions */
	const double cosgha = cos(gha);
	const double singha = sin(gha);
	const double cosdec = cos(dec);
	const double sindec = sin(dec);
	const double cospsi = cos(psi);
	const double sinpsi = sin(psi);

	/* Eq. (B4) of [ABCF].  Note that dec = pi/2 - theta, and gha =
	 * -phi where theta and phi are the standard spherical coordinates
	 * used in that paper. */
	X[0] = -cospsi * singha - sinpsi * cosgha * sindec;
	X[1] = -cospsi * cosgha + sinpsi * singha * sindec;
	X[2] = sinpsi * cosdec;

	/* Eq. (B5) of [ABCF].  Note that dec = pi/2 - theta, and gha =
	 * -phi where theta and phi are the standard spherical coordinates
	 * used in that paper. */
	Y[0] = sinpsi * singha - cospsi * cosgha * sindec;
	Y[1] = sinpsi * cosgha + cospsi * singha * sindec;
	Y[2] = cospsi * cosdec;

	/* Eqns from [Blaut2012] - but converted from Blaut's theta = pi/2 - dec
	 * given cos(dec) = sin(theta) and cos(theta) = sin(dec), and Blaut's phi = ra
	 * - gmst, so cos(phi) = cos(gha) and sin(phi) = -sin(gha). This is consistent
	 * with the convention in XLALComputeDetAMResponse.
	 */
	Z[0] = -cosdec * cosgha;
	Z[1] = cosdec * singha;
	Z[2] = -sindec;

	/* Now compute Eq. (B7) of [ABCF] for each polarization state, i.e.,
	 * with s+=1 and sx=0 to get F+, with s+=0 and sx=1 to get Fx */
	*fplus = *fcross = *fb = *fl = *fx = *fy = 0.0;
	for(i = 0; i < 3; i++) {
		const double DX = D[i][0] * X[0] + D[i][1] * X[1] + D[i][2] * X[2];
		const double DY = D[i][0] * Y[0] + D[i][1] * Y[1] + D[i][2] * Y[2];
		const double DZ = D[i][0] * Z[0] + D[i][1] * Z[1] + D[i][2] * Z[2];

		*fplus += X[i] * DX - Y[i] * DY;
		*fcross += X[i] * DY + Y[i] * DX;

		/* scalar and vector modes from [Nishizawa2009] */
		*fb += X[i] * DX + Y[i] * DY;
		*fl += LAL_SQRT2 * Z[i] * DZ;
		*fx += X[i] * DZ + Z[i] * DX;
		*fy += Y[i] * DZ + Z[i] * DY;
	}
}

/**
 * \deprecated Use XLALComputeDetAMResponse() instead.
 */
void LALComputeDetAMResponse(LALStatus * status, LALDetAMResponse * pResponse, const LALDetAndSource * pDetAndSrc, const LIGOTimeGPS * gps)
{
	double fplus, fcross;

	INITSTATUS(status);
	ATTATCHSTATUSPTR(status);

	ASSERT(pResponse != (LALDetAMResponse *) NULL, status, DETRESPONSEH_ENULLOUTPUT, DETRESPONSEH_MSGENULLOUTPUT);

	ASSERT(pDetAndSrc != NULL, status, DETRESPONSEH_ENULLINPUT, DETRESPONSEH_MSGENULLINPUT);

	ASSERT(gps != (LIGOTimeGPS *) NULL, status, DETRESPONSEH_ENULLINPUT, DETRESPONSEH_MSGENULLINPUT);

	/* source coordinates must be in equatorial system */
	ASSERT(pDetAndSrc->pSource->equatorialCoords.system == COORDINATESYSTEM_EQUATORIAL, status, DETRESPONSEH_ESRCNOTEQUATORIAL, DETRESPONSEH_MSGESRCNOTEQUATORIAL);

	XLALComputeDetAMResponse(&fplus, &fcross, pDetAndSrc->pDetector->response, pDetAndSrc->pSource->equatorialCoords.longitude, pDetAndSrc->pSource->equatorialCoords.latitude, pDetAndSrc->pSource->orientation, XLALGreenwichMeanSiderealTime(gps));

	pResponse->plus = fplus;
	pResponse->cross = fcross;
	pResponse->scalar = 0.0;	/* not implemented */

	DETATCHSTATUSPTR(status);
	RETURN(status);
}


/**
 * Computes REAL4TimeSeries containing time series of response amplitudes.
 * \see XLALComputeDetAMResponse() for more details.
 */
int XLALComputeDetAMResponseSeries(REAL4TimeSeries ** fplus, REAL4TimeSeries ** fcross, const REAL4 D[3][3], const double ra, const double dec, const double psi, const LIGOTimeGPS * start, const double deltaT, const int n)
{
	LIGOTimeGPS t;
	double gmst;
	int i;
	double p, c;

	*fplus = XLALCreateREAL4TimeSeries("plus", start, 0.0, deltaT, &lalDimensionlessUnit, n);
	*fcross = XLALCreateREAL4TimeSeries("cross", start, 0.0, deltaT, &lalDimensionlessUnit, n);
	if(!*fplus || !*fcross) {
		XLALDestroyREAL4TimeSeries(*fplus);
		XLALDestroyREAL4TimeSeries(*fcross);
		*fplus = *fcross = NULL;
		XLAL_ERROR(XLAL_EFUNC);
	}

	for(i = 0; i < n; i++) {
		t = *start;
		gmst = XLALGreenwichMeanSiderealTime(XLALGPSAdd(&t, i * deltaT));
		if(XLAL_IS_REAL8_FAIL_NAN(gmst)) {
			XLALDestroyREAL4TimeSeries(*fplus);
			XLALDestroyREAL4TimeSeries(*fcross);
			*fplus = *fcross = NULL;
			XLAL_ERROR(XLAL_EFUNC);
		}
		XLALComputeDetAMResponse(&p, &c, D, ra, dec, psi, gmst);
		(*fplus)->data->data[i] = p;
		(*fcross)->data->data[i] = c;
	}

	return 0;
}

/**
 * Computes REAL4TimeSeries containing time series of the full general
 * metric theory of gravity response amplitudes.
 * \see XLALComputeDetAMResponseExtraModes() for more details.
 */
int XLALComputeDetAMResponseExtraModesSeries(REAL4TimeSeries ** fplus, REAL4TimeSeries ** fcross, REAL4TimeSeries ** fb, REAL4TimeSeries ** fl, REAL4TimeSeries ** fx, REAL4TimeSeries ** fy, REAL4 D[3][3], const double ra, const double dec, const double psi, const LIGOTimeGPS * start, const double deltaT, const int n)
{
	LIGOTimeGPS t;
	double gmst;
	int i;
	double p, c, b, l, x, y;

	*fplus = XLALCreateREAL4TimeSeries("plus", start, 0.0, deltaT, &lalDimensionlessUnit, n);
	*fcross = XLALCreateREAL4TimeSeries("cross", start, 0.0, deltaT, &lalDimensionlessUnit, n);
	*fb = XLALCreateREAL4TimeSeries("b", start, 0.0, deltaT, &lalDimensionlessUnit, n);
	*fl = XLALCreateREAL4TimeSeries("l", start, 0.0, deltaT, &lalDimensionlessUnit, n);
	*fx = XLALCreateREAL4TimeSeries("x", start, 0.0, deltaT, &lalDimensionlessUnit, n);
	*fy = XLALCreateREAL4TimeSeries("y", start, 0.0, deltaT, &lalDimensionlessUnit, n);

	if(!*fplus || !*fcross || !*fb || !*fl || !*fx || !*fy) {
		XLALDestroyREAL4TimeSeries(*fplus);
		XLALDestroyREAL4TimeSeries(*fcross);
		XLALDestroyREAL4TimeSeries(*fb);
		XLALDestroyREAL4TimeSeries(*fl);
		XLALDestroyREAL4TimeSeries(*fx);
		XLALDestroyREAL4TimeSeries(*fy);

		*fplus = *fcross = *fb = *fl = *fx = *fy = NULL;
		XLAL_ERROR(XLAL_EFUNC);
	}

	for(i = 0; i < n; i++) {
		t = *start;
		gmst = XLALGreenwichMeanSiderealTime(XLALGPSAdd(&t, i * deltaT));
		if(XLAL_IS_REAL8_FAIL_NAN(gmst)) {
			XLALDestroyREAL4TimeSeries(*fplus);
			XLALDestroyREAL4TimeSeries(*fcross);
			XLALDestroyREAL4TimeSeries(*fb);
			XLALDestroyREAL4TimeSeries(*fl);
			XLALDestroyREAL4TimeSeries(*fx);
			XLALDestroyREAL4TimeSeries(*fy);

			*fplus = *fcross = *fb = *fl = *fx = *fy = NULL;
			XLAL_ERROR(XLAL_EFUNC);
		}
		XLALComputeDetAMResponseExtraModes(&p, &c, &b, &l, &x, &y, D, ra, dec, psi, gmst);
		(*fplus)->data->data[i] = p;
		(*fcross)->data->data[i] = c;
		(*fb)->data->data[i] = b;
		(*fl)->data->data[i] = l;
		(*fx)->data->data[i] = x;
		(*fy)->data->data[i] = y;
	}

	return 0;
}

/**
 * Computes REAL4TimeSeries containing time series of response amplitudes.
 * \deprecated Use XLALComputeDetAMResponseSeries() instead.
 */
void LALComputeDetAMResponseSeries(LALStatus * status, LALDetAMResponseSeries * pResponseSeries, const LALDetAndSource * pDetAndSource, const LALTimeIntervalAndNSample * pTimeInfo)
{
	/* Want to loop over the time and call LALComputeDetAMResponse() */
	LALDetAMResponse instResponse;
	unsigned i;
	char infostr[128];

	INITSTATUS(status);
	ATTATCHSTATUSPTR(status);

	if(lalDebugLevel >= 8) {
		sprintf(infostr, "pResponseSeries->pPlus->data->length = %d\npTimeInfo->nSample = %d\n", pResponseSeries->pPlus->data->length, pTimeInfo->nSample);
		LALInfo(status, infostr);
	}

	/*
	 * Error-checking assertions
	 */
	ASSERT(pResponseSeries != (LALDetAMResponseSeries *) NULL, status, DETRESPONSEH_ENULLOUTPUT, DETRESPONSEH_MSGENULLOUTPUT);

	ASSERT(pDetAndSource != (LALDetAndSource *) NULL, status, DETRESPONSEH_ENULLINPUT, DETRESPONSEH_MSGENULLINPUT);

	ASSERT(pTimeInfo != (LALTimeIntervalAndNSample *) NULL, status, DETRESPONSEH_ENULLINPUT, DETRESPONSEH_MSGENULLINPUT);

	/*
	 * Set names
	 */
	pResponseSeries->pPlus->name[0] = '\0';
	pResponseSeries->pCross->name[0] = '\0';
	pResponseSeries->pScalar->name[0] = '\0';

	strncpy(pResponseSeries->pPlus->name, "plus", LALNameLength);
	strncpy(pResponseSeries->pCross->name, "cross", LALNameLength);
	strncpy(pResponseSeries->pScalar->name, "scalar", LALNameLength);

	/*
	 * Set sampling parameters
	 */
	pResponseSeries->pPlus->epoch = pTimeInfo->epoch;
	pResponseSeries->pPlus->deltaT = pTimeInfo->deltaT;
	pResponseSeries->pPlus->f0 = 0.;
	pResponseSeries->pPlus->sampleUnits = lalDimensionlessUnit;

	pResponseSeries->pCross->epoch = pTimeInfo->epoch;
	pResponseSeries->pCross->deltaT = pTimeInfo->deltaT;
	pResponseSeries->pCross->f0 = 0.;
	pResponseSeries->pCross->sampleUnits = lalDimensionlessUnit;

	pResponseSeries->pScalar->epoch = pTimeInfo->epoch;
	pResponseSeries->pScalar->deltaT = pTimeInfo->deltaT;
	pResponseSeries->pScalar->f0 = 0.;
	pResponseSeries->pScalar->sampleUnits = lalDimensionlessUnit;

	/*
	 * Ensure enough memory for requested vectors
	 */
	if(pResponseSeries->pPlus->data->length < pTimeInfo->nSample) {
		if(lalDebugLevel >= 8)
			LALInfo(status, "plus sequence too short -- reallocating");

		TRY(LALSDestroyVector(status->statusPtr, &(pResponseSeries->pPlus->data)), status);

		TRY(LALSCreateVector(status->statusPtr, &(pResponseSeries->pPlus->data), pTimeInfo->nSample), status);

		if(lalDebugLevel > 0)
			printf("pResponseSeries->pPlus->data->length = %d\n", pResponseSeries->pPlus->data->length);

	}

	if(pResponseSeries->pCross->data->length < pTimeInfo->nSample) {
		if(lalDebugLevel >= 8)
			LALInfo(status, "cross sequence too short -- reallocating");

		TRY(LALSDestroyVector(status->statusPtr, &(pResponseSeries->pCross->data)), status);

		TRY(LALSCreateVector(status->statusPtr, &(pResponseSeries->pCross->data), pTimeInfo->nSample), status);

	}

	if(pResponseSeries->pScalar->data->length < pTimeInfo->nSample) {
		if(lalDebugLevel & 0x08)
			LALInfo(status, "scalar sequence too short -- reallocating");

		TRY(LALSDestroyVector(status->statusPtr, &(pResponseSeries->pScalar->data)), status);

		TRY(LALSCreateVector(status->statusPtr, &(pResponseSeries->pScalar->data), pTimeInfo->nSample), status);
	}


	/*
	 * Loop to compute each element in time series.
	 */

	for(i = 0; i < pTimeInfo->nSample; ++i) {
		LIGOTimeGPS gps = pTimeInfo->epoch;
		XLALGPSAdd(&gps, i * pTimeInfo->deltaT);

		if(lalDebugLevel >= 8) {
			sprintf(infostr, "LALComputeDetAMResponseSeries: i = %d\n", i);
			LALInfo(status, infostr);
		}

		TRY(LALComputeDetAMResponse(status->statusPtr, &instResponse, pDetAndSource, &gps), status);

		pResponseSeries->pPlus->data->data[i] = instResponse.plus;
		pResponseSeries->pCross->data->data[i] = instResponse.cross;
		pResponseSeries->pScalar->data->data[i] = instResponse.scalar;
	}

	DETATCHSTATUSPTR(status);
	RETURN(status);
}
