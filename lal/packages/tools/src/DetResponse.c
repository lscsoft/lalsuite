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


/*
 * Matrix operations.
 */

static double matrix_dot(REAL4 a[3][3], double b[3][3])
{
	return a[0][0] * b[0][0] + a[0][1] * b[0][1] + a[0][2] * b[0][2] + a[1][0] * b[1][0] + a[1][1] * b[1][1] + a[1][2] * b[1][2] + a[2][0] * b[2][0] + a[2][1] * b[2][1] + a[2][2] * b[2][2];
}


/* rotate a matrix by computing rot * M * rot^T where rot is a rotation
 * matrix.  assumes that M[][2] and M[2][] = {0,0,0} */
static void matrix_rot(double out[3][3], double rot[3][3], const double M[3][3])
{
	double t[3][2];
	int i, j;

	for(i = 0; i < 3; i++)
		for(j = 0; j < 2; j++)
			t[i][j] = rot[i][0] * M[0][j] + rot[i][1] * M[1][j];
	for(i = 0; i < 3; i++)
		for(j = 0; j < 3; j++)
			out[i][j] = t[i][0] * rot[j][0] + t[i][1] * rot[j][1];
}


/*
 * Compute detector response, making use of the detector response tensor
 * stored in LALDetector
 */

/* <lalVerbatim file="DetResponseCP"> */
LALDetAMResponse *XLALComputeDetAMResponse(LALDetAMResponse *output, REAL4 response[3][3], const double right_ascension, const double declination, const double orientation, const double gmst)
{				/* </lalVerbatim> */
	/* unit strains */
	const double e_plus[3][3] = {{1.0, 0.0, 0.0}, {0.0, -1.0, 0.0}, {0.0, 0.0, 0.0}};
	const double e_cros[3][3] = {{0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

	/* rotated to Earth-fixed frame */
	double e_plus_E[3][3];
	double e_cros_E[3][3];

	/* rotation matrix */
	double rot[3][3];

	/* Greenwich hour angle of source */
	const double gha = gmst - right_ascension;

	/* pre-computed trig functions */
	const double singha = sin(gha);
	const double cosgha = cos(gha);
	const double sinorient = sin(orientation);
	const double cosorient = cos(orientation);
	const double sindec = sin(declination);
	const double cosdec = cos(declination);

	/*
	 * Construct the rotation matrix to convert the source perterbation
	 * matrix into the Earth-fixed basis.  This rotation matrix equals
	 *
	 * R_z(gha + pi/2) * R_x(-declination - pi/2) * R_z(-orientation)
	 *
	 * where R_i is the matrix for a right-hand rotation about axis i.
	 * Some elements are set to 0 because they are not needed due to
	 * zeros in the e+ and ex matrices.
	 */

	rot[0][0] = -singha * cosorient - cosgha * sinorient * sindec;
	rot[0][1] =  singha * sinorient - cosgha * cosorient * sindec;
	rot[0][2] =  0.0;	/* -cosgha * cosdec */
	rot[1][0] = -cosgha * cosorient + singha * sinorient * sindec;
	rot[1][1] =  cosgha * sinorient + singha * cosorient * sindec;
	rot[1][2] =  0.0;	/* singha * cosdec */
	rot[2][0] =  sinorient * cosdec;
	rot[2][1] =  cosorient * cosdec;
	rot[2][2] =  0.0;	/* -sindec */

	/*
	 * Now, get the unit perturbation tensors in the Earth fixed frame.
	 *
	 *    e_plus_Earth = rot &* e_plus &* transpose(rot)
	 *    e_cros_Earth = rot &* e_cros &* transpose(rot)
	 */

	matrix_rot(e_plus_E, rot, e_plus);
	matrix_rot(e_cros_E, rot, e_cros);

	/*
	 * Then, F_plus = det_response &. e_plus_E, F_cros = det_response
	 * &. e_cros_E
	 */

	output->plus = matrix_dot(response, e_plus_E);
	output->cross = matrix_dot(response, e_cros_E);
	/* FIXME: scalar response not implemented, yet.  Will have to read
	 * Waggoner's paper to do this. */
	output->scalar = 0.0;

	return output;
}

/* <lalVerbatim file="DetResponseCP"> */
void LALComputeDetAMResponse(LALStatus * status, LALDetAMResponse * pResponse, const LALDetAndSource * pDetAndSrc, const LALGPSandAcc * pGPSandAcc)
{				/* </lalVerbatim> */
	INITSTATUS(status, "LALComputeDetAMResponse", DETRESPONSEC);
	ATTATCHSTATUSPTR(status);

	ASSERT(pResponse != (LALDetAMResponse *) NULL, status, DETRESPONSEH_ENULLOUTPUT, DETRESPONSEH_MSGENULLOUTPUT);

	ASSERT(pDetAndSrc != NULL, status, DETRESPONSEH_ENULLINPUT, DETRESPONSEH_MSGENULLINPUT);

	ASSERT(pGPSandAcc != (LALGPSandAcc *) NULL, status, DETRESPONSEH_ENULLINPUT, DETRESPONSEH_MSGENULLINPUT);

	/* source coordinates must be in equatorial system */
	ASSERT(pDetAndSrc->pSource->equatorialCoords.system == COORDINATESYSTEM_EQUATORIAL, status, DETRESPONSEH_ESRCNOTEQUATORIAL, DETRESPONSEH_MSGESRCNOTEQUATORIAL);

	XLALComputeDetAMResponse(pResponse, pDetAndSrc->pDetector->response, pDetAndSrc->pSource->equatorialCoords.longitude, pDetAndSrc->pSource->equatorialCoords.latitude, pDetAndSrc->pSource->orientation, XLALGreenwichMeanSiderealTime(&pGPSandAcc->gps));

	DETATCHSTATUSPTR(status);
	RETURN(status);
}


/*
 * Computes REAL4TimeSeries containing time series of response amplitudes.
 */

int XLALComputeDetAMResponseSeries(REAL4TimeSeries **plus, REAL4TimeSeries **cross, REAL4 response[3][3], const double right_ascension, const double declination, const double orientation, const LIGOTimeGPS *start, const double deltaT, const int n)
{
	static const char func[] = "XLALComputeDetAMResponseSeries";
	LIGOTimeGPS t;
	double gmst;
	int i;
	LALDetAMResponse resp;

	*plus = XLALCreateREAL4TimeSeries("plus", start, 0.0, deltaT, &lalDimensionlessUnit, n);
	*cross = XLALCreateREAL4TimeSeries("cross", start, 0.0, deltaT, &lalDimensionlessUnit, n);
	if(!*plus || !*cross) {
		XLALDestroyREAL4TimeSeries(*plus);
		XLALDestroyREAL4TimeSeries(*cross);
		*plus = *cross = NULL;
		XLAL_ERROR(func, XLAL_EFUNC);
	}

	for(i = 0; i < n; i++) {
		t = *start;
		gmst = XLALGreenwichMeanSiderealTime(XLALAddFloatToGPS(&t, i * deltaT));
		if(XLAL_IS_REAL8_FAIL_NAN(gmst)) {
			XLALDestroyREAL4TimeSeries(*plus);
			XLALDestroyREAL4TimeSeries(*cross);
			*plus = *cross = NULL;
			XLAL_ERROR(func, XLAL_EFUNC);
		}
		XLALComputeDetAMResponse(&resp, response, right_ascension, declination, orientation, gmst);
		(*plus)->data->data[i] = resp.plus;
		(*cross)->data->data[i] = resp.cross;
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
