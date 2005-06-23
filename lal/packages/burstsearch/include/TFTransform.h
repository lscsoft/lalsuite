/********************************** <lalVerbatim file="TFTransformHV">
Author: Flanagan, E
$Id$
**************************************************** </lalVerbatim> */

#ifndef _TFTRANSFORM_H
#define _TFTRANSFORM_H

#include <lal/LALDatatypes.h>
#include <lal/Window.h>
#include <lal/RealFFT.h>
#include <lal/LALRCSID.h>

#ifdef  __cplusplus		/* C++ protection. */
extern "C" {
#endif

NRCSID(TFTRANSFORMH, "$Id$");


typedef struct tagTFPlaneParams {
	INT4 timeBins;	/* Number of time bins in TF plane    */
	INT4 freqBins;	/* Number of freq bins in TF plane    */
	REAL8 deltaT;	/* time resolution of the plane     */
	REAL8 deltaF;	/* freq. resolution of the plane */
	REAL8 flow;	/* minimum frequency to search for */
	REAL8 fhigh;	/* maximum frequency to search for */
	REAL8 timeDuration;	/* length of data to be used to create a TF plane at a time */
} TFPlaneParams;


typedef struct tagCOMPLEX8TimeFrequencyPlane {
	CHAR *name;
	LIGOTimeGPS epoch;
	CHARVector *sampleUnits;
	TFPlaneParams params;
	COMPLEX8 *data;
	/*
	 * data[i*params->freqBins+j] is a complex number
	 * corresponding to a time t_i = epoch + i*(deltaT)
	 * and a frequency f_j = flow + j / (deltaT)
	 */
} COMPLEX8TimeFrequencyPlane;


int
XLALComputeFrequencySeries(
	COMPLEX8FrequencySeries *freqSeries,
	const REAL4TimeSeries *timeSeries,
	const REAL4Window *window,
	const REAL4FFTPlan *plan
);


COMPLEX8TimeFrequencyPlane *
XLALCreateTFPlane(
	TFPlaneParams *params,
	INT4 minFreqBins
);


void
XLALDestroyTFPlane(
	COMPLEX8TimeFrequencyPlane *plane
);


int
XLALFreqSeriesToTFPlane(
	COMPLEX8TimeFrequencyPlane *tfp,
	const COMPLEX8FrequencySeries *freqSeries,
	UINT4 windowShift,
	REAL4 *norm,
	const REAL4FrequencySeries *psd
);

#ifdef  __cplusplus
}
#endif				/* C++ protection. */
#endif
