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
#include <lal/Sequence.h>

#ifdef  __cplusplus		/* C++ protection. */
extern "C" {
#endif

NRCSID(TFTRANSFORMH, "$Id$");


typedef struct tagREAL4TimeFrequencyPlane {
	/* name of data from which this was computed */
	CHAR name[LALNameLength];
	/* epoch of data from which this was computed */
	LIGOTimeGPS epoch;
	/* time resolution of the plane */
	REAL8 deltaT;
	/* Number of frequency channels in TF plane */
	UINT4 channels;
	/* TF plane's frequency resolution (channel spacing) */
	REAL8 deltaF;
	/* low frequency boundary of TF plane */
	REAL8 flow;
	/* inner product of filters for neighbouring channels */
	REAL8 channel_overlap;
	/* predicted mean square for the time series in each channel */
	REAL4Sequence *channel_mean_square;
	/* channel data;  channel[j]->data[i] corresponds to time epoch + i
	 * * deltaT and frequency flow + j * deltaF */
	REAL4Sequence **channel;
} REAL4TimeFrequencyPlane;


REAL4TimeFrequencyPlane *
XLALCreateTFPlane(
	UINT4 timeBins,
	REAL8 deltaT,
	UINT4 channels,
	REAL8 deltaF,
	REAL8 flow
);


void
XLALDestroyTFPlane(
	REAL4TimeFrequencyPlane *plane
);


COMPLEX8FrequencySeries *XLALWindowedREAL4ForwardFFT(
	const REAL4TimeSeries *tseries,
	const REAL4Window *window,
	const REAL4FFTPlan *plan
);

COMPLEX16FrequencySeries *XLALWindowedREAL8ForwardFFT(
	const REAL8TimeSeries *tseries,
	const REAL8Window *window,
	const REAL8FFTPlan *plan
);


int
XLALFreqSeriesToTFPlane(
	REAL4TimeFrequencyPlane *tfplane,
	const COMPLEX8FrequencySeries *fseries,
	const REAL4FrequencySeries *psd,
	const REAL4FFTPlan *reverseplan
);


#ifdef  __cplusplus
}
#endif				/* C++ protection. */
#endif
