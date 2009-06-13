/*
*  Copyright (C) 2007 David Chin, Jolien Creighton, Kipp Cannon, Peter Shawhan
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

/*<lalVerbatim file="DetResponseCV">

Author: David Chin <dwchin@umich.edu> +1-734-709-9119, Kipp Cannon <kipp@gravity.phys.uwm.edu>
$Id$

</lalVerbatim> */

/*
<lalLaTeX>

\subsection{Module \texttt{DetResponse.c}}
\label{ss:DetResponse.c}

Computes the response of a detector.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{DetResponseCP}
\idx{LALComputeDetAMResponse()}
\idx{XLALComputeDetAMResponse()}
\idx{LALComputeDetAMResponseSeries()}
\idx{XLALComputeDetAMResponseSeries()}

\subsubsection*{Description}

These routines compute the antenna beam pattern for all supported detector
types.  \texttt{XLALComputeDetAMResponse()} computes the response at one
instance in time, and \texttt{XLALComputeDetAMResponseSeries()} computes a
vector of response for some length of time.

\subsubsection*{Algorithm}

This code is a translation of the algorithm in the Maple worksheet by
Anderson, \textit{et al.}~\cite{tools:Anderson:2000}.  We compute the $h$-tensors for
$+$- and $\times$-polarized in the Earth-fixed frame, and then contract
them (take the scalar product) with the detector response tensors as
described in the \texttt{DetectorSite.h} section of the \texttt{tools}
package.

\texttt{DetectorSite.h} in the \texttt{tools} package  provides predefined
\texttt{LALDetector} structures representing most current detectors,
including LIGO (Hanford and Livingston), and GEO.

\subsubsection*{Uses}
\texttt{LALGPStoGMST1()}

\subsubsection*{Notes}

For examples of usage, please see the test programs in the \texttt{test}
directory.

\vfill{\footnotesize\input{DetResponseCV}}

</lalLaTeX>
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

NRCSID(DETRESPONSEC, "$Id$");


/** XLALComputeDetAMResponse
 *
 * An implementation of the detector response formulae in Anderson et al
 * PRD 63 042003 (2001) [ABCF].
 *
 * Computes F+ and Fx for a source at a specified sky position,
 * polarization angle, and sidereal time.  Also requires the detector's
 * response matrix which is defined by Eq. (B6) of [ABCF] using either
 * Table 1 of [ABCF] or Eqs. (B11)--(B17) of [ABCF] to compute the arm
 * direction unit vectors.
 */

/* <lalVerbatim file="DetResponseCP"> */
void XLALComputeDetAMResponse(
	double *fplus,		/**< Returned value of F+ */
	double *fcross, 	/**< Returned value of Fx */
	REAL4 D[3][3],		/**< Detector response 3x3 matrix */
	const double ra,	/**< Right ascention of source (radians) */
	const double dec,	/**< Declination of source (radians) */
	const double psi,	/**< Polarization angle of source (radians) */
	const double gmst	/**< Greenwich mean sidereal time (radians) */
)
{				/* </lalVerbatim> */
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


/* <lalVerbatim file="DetResponseCP"> */
void LALComputeDetAMResponse(LALStatus * status, LALDetAMResponse * pResponse, const LALDetAndSource * pDetAndSrc, const LALGPSandAcc * pGPSandAcc)
{				/* </lalVerbatim> */
	double fplus, fcross;

	INITSTATUS(status, "LALComputeDetAMResponse", DETRESPONSEC);
	ATTATCHSTATUSPTR(status);

	ASSERT(pResponse != (LALDetAMResponse *) NULL, status, DETRESPONSEH_ENULLOUTPUT, DETRESPONSEH_MSGENULLOUTPUT);

	ASSERT(pDetAndSrc != NULL, status, DETRESPONSEH_ENULLINPUT, DETRESPONSEH_MSGENULLINPUT);

	ASSERT(pGPSandAcc != (LALGPSandAcc *) NULL, status, DETRESPONSEH_ENULLINPUT, DETRESPONSEH_MSGENULLINPUT);

	/* source coordinates must be in equatorial system */
	ASSERT(pDetAndSrc->pSource->equatorialCoords.system == COORDINATESYSTEM_EQUATORIAL, status, DETRESPONSEH_ESRCNOTEQUATORIAL, DETRESPONSEH_MSGESRCNOTEQUATORIAL);

	XLALComputeDetAMResponse(&fplus, &fcross, pDetAndSrc->pDetector->response, pDetAndSrc->pSource->equatorialCoords.longitude, pDetAndSrc->pSource->equatorialCoords.latitude, pDetAndSrc->pSource->orientation, XLALGreenwichMeanSiderealTime(&pGPSandAcc->gps));

	pResponse->plus = fplus;
	pResponse->cross = fcross;
	pResponse->scalar = 0.0;	/* not implemented */

	DETATCHSTATUSPTR(status);
	RETURN(status);
}


/*
 * Computes REAL4TimeSeries containing time series of response amplitudes.
 */

int XLALComputeDetAMResponseSeries(
	REAL4TimeSeries **fplus,
	REAL4TimeSeries **fcross,
	REAL4 D[3][3],
	const double ra,
	const double dec,
	const double psi,
	const LIGOTimeGPS *start,
	const double deltaT,
	const int n
)
{
	static const char func[] = "XLALComputeDetAMResponseSeries";
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
		XLAL_ERROR(func, XLAL_EFUNC);
	}

	for(i = 0; i < n; i++) {
		t = *start;
		gmst = XLALGreenwichMeanSiderealTime(XLALGPSAdd(&t, i * deltaT));
		if(XLAL_IS_REAL8_FAIL_NAN(gmst)) {
			XLALDestroyREAL4TimeSeries(*fplus);
			XLALDestroyREAL4TimeSeries(*fcross);
			*fplus = *fcross = NULL;
			XLAL_ERROR(func, XLAL_EFUNC);
		}
		XLALComputeDetAMResponse(&p, &c, D, ra, dec, psi, gmst);
		(*fplus)->data->data[i] = p;
		(*fcross)->data->data[i] = c;
	}

	return 0;
}


/*
 * Computes REAL4TimeSeries containing time series of response amplitudes.
 */

/* <lalVerbatim file="DetResponseCP"> */
void LALComputeDetAMResponseSeries(LALStatus * status, LALDetAMResponseSeries * pResponseSeries, const LALDetAndSource * pDetAndSource, const LALTimeIntervalAndNSample * pTimeInfo)
{				/* </lalVerbatim> */
	/* Want to loop over the time and call LALComputeDetAMResponse() */
	LALDetAMResponse instResponse;
	LIGOTimeGPS gps;
	LIGOTimeGPS tmpgps;
	LALGPSandAcc gps_and_acc;
	LALTimeInterval dt;
	unsigned i;
	char infostr[128];

	INITSTATUS(status, "LALComputeDetAMResponseSeries", DETRESPONSEC);
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
	 * Loop to compute each element in time series; rint(3) is a std C
	 * function that rounds floating point numbers properly.
	 */
	gps = pTimeInfo->epoch;

	TRY(LALFloatToInterval(status->statusPtr, &dt, &(pTimeInfo->deltaT)), status);

	for(i = 0; i < pTimeInfo->nSample; ++i) {
		if(lalDebugLevel >= 8) {
			sprintf(infostr, "LALComputeDetAMResponseSeries: i = %d\n", i);
			LALInfo(status, infostr);
		}

		gps_and_acc.gps = gps;
		gps_and_acc.accuracy = pTimeInfo->accuracy;
		TRY(LALComputeDetAMResponse(status->statusPtr, &instResponse, pDetAndSource, &gps_and_acc), status);

		pResponseSeries->pPlus->data->data[i] = instResponse.plus;
		pResponseSeries->pCross->data->data[i] = instResponse.cross;
		pResponseSeries->pScalar->data->data[i] = instResponse.scalar;

		tmpgps = gps;

		TRY(LALIncrementGPS(status->statusPtr, &gps, &tmpgps, &dt), status);

	}

	DETATCHSTATUSPTR(status);
	RETURN(status);
}
