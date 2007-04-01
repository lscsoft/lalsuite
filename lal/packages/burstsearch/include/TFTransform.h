/********************************** <lalVerbatim file="TFTransformHV">
Author: Flanagan, E., Kipp Cannon
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


/*
 * A tiling of the time-frequency plane
 */


typedef struct tagTFTile {
	/* tile specification as indexes into time-frequency plane data */
	UINT4 channel0;
	UINT4 channels;
	UINT4 tstart;
	UINT4 tend;
	/* time-frequency plane parameters for reconstructing tile dimensions */
	REAL8 flow;
	REAL8 deltaT;
	REAL8 deltaF;
	/* computed tile properties */
	REAL8 excessPower;
	REAL8 hrss;
	REAL8 lnweight;
	/* -ln P(event | stationary Gaussian white noise) */
	REAL8 confidence;
} TFTile;


typedef struct tagTFTiling {
	/* array of tiles */
	TFTile *tile;
	size_t numtiles;
} TFTiling;


REAL8
XLALTFTileDegreesOfFreedom(
	const TFTile *tile
);


TFTiling *
XLALCreateTFTiling(
	UINT4 plane_length,
	REAL8 plane_deltaT,
	REAL8 plane_flow,
	REAL8 plane_deltaF,
	UINT4 plane_num_channels,
	UINT4 tiling_tstart,
	UINT4 inv_fractional_stride,
	REAL8 maxTileBandwidth,
	REAL8 maxTileDuration
);


void
XLALDestroyTFTiling(
	TFTiling *tiling
);


REAL8
XLALComputeLikelihood(
	TFTiling *tiling
);


/*
 * A time-frequency plane
 */


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
	/* inner product of filters for neighbouring channels;
	 * channel_overlap[0] is the inner product of the filters for
	 * channels 0 and 1, and so on */
	REAL8Sequence *channel_overlap;
	/* predicted root mean square for the time series in each channel */
	REAL8Sequence *channel_rms;
	/* channel data;  channel[j]->data[i] corresponds to time epoch + i
	 * * deltaT and frequency flow + j * deltaF */
	REAL4Sequence **channel;
	/* time-frequency plane's tiling */
	TFTiling *tiling;
} REAL4TimeFrequencyPlane;


REAL4TimeFrequencyPlane *
XLALCreateTFPlane(
	UINT4 tseries_length,
	REAL8 tseries_deltaT,
	UINT4 channels,
	REAL8 deltaF,
	REAL8 flow,
	UINT4 tiling_start,
	INT4 tiling_inv_fractional_stride,
	REAL8 tiling_max_bandwidth,
	REAL8 tiling_max_duration
);


void
XLALDestroyTFPlane(
	REAL4TimeFrequencyPlane *plane
);


int
XLALFreqSeriesToTFPlane(
	REAL4TimeFrequencyPlane *tfplane,
	const COMPLEX8FrequencySeries *fseries,
	const REAL4FrequencySeries *psd,
	const REAL4FFTPlan *reverseplan,
	INT4 enable_over_whitening
);


int
XLALComputeExcessPower(
	const REAL4TimeFrequencyPlane *plane
);


#ifdef  __cplusplus
}
#endif				/* C++ protection. */
#endif
