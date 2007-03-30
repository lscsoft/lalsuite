/******** <lalVerbatim file="CreateTFPlaneCV"> ********
Author: Flanagan, E. and Cannon, K.
$Id$
********* </lalVerbatim> ********/

#include <lal/LALRCSID.h>

NRCSID(CREATETFPLANEC, "$Id$");

#include <lal/LALMalloc.h>
#include <lal/LALStdlib.h>
#include <lal/TFTransform.h>
#include <lal/Sequence.h>
#include <lal/XLALError.h>


/******** <lalVerbatim file="CreateTFPlaneCP"> ********/
REAL4TimeFrequencyPlane *XLALCreateTFPlane(
	/* length of time series from which TF plane will be computed */
	UINT4 tseries_length,
	/* sample rate of time series */
	REAL8 tseries_deltaT,
	/* number of frequency channels in TF plane */
	UINT4 channels,
	/* frequency resolution of the plane */
	REAL8 deltaF,
	/* minimum frequency to search for */
	REAL8 flow,
	/* sample at which to start the tiling */
	UINT4 tiling_start,
	/* overlap of adjacent tiles */
	INT4 tiling_inv_fractional_stride,
	/* largest tile's bandwidth */
	REAL8 tiling_max_bandwidth,
	/* largest tile's duration */
	REAL8 tiling_max_duration
)
/******** </lalVerbatim> ********/
{
	static const char func[] = "XLALCreateTFPlane";
	REAL4TimeFrequencyPlane *plane;
	REAL8Sequence *channel_overlap;
	REAL8Sequence *channel_rms;
	REAL4Sequence **channel;
	TFTiling *tiling;
	unsigned i;

	/*
	 * Make sure that input parameters are reasonable
	 */

	if((flow < 0) ||
	   (deltaF <= 0.0) ||
	   (tseries_deltaT <= 0.0) ||
	   (channels < 1))
		XLAL_ERROR_NULL(func, XLAL_EINVAL);

	/*
	 * Allocate memory and construct the tiling.
	 */

	plane = LALMalloc(sizeof(*plane));
	channel_overlap = XLALCreateREAL8Sequence(channels - 1);
	channel_rms = XLALCreateREAL8Sequence(channels);
	channel = LALMalloc(channels * sizeof(*channel));
	tiling = XLALCreateTFTiling(tseries_length, tseries_deltaT, flow, deltaF, channels, tiling_start, tiling_inv_fractional_stride, tiling_max_bandwidth, tiling_max_duration);
	if(!plane || !channel_overlap || !channel_rms || !channel || !tiling) {
		LALFree(plane);
		XLALDestroyREAL8Sequence(channel_overlap);
		XLALDestroyREAL8Sequence(channel_rms);
		LALFree(channel);
		XLALDestroyTFTiling(tiling);
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);
	}
	for(i = 0; i < channels; i++) {
		channel[i] = XLALCreateREAL4Sequence(tseries_length);
		if(!channel[i]) {
			while(--i)
				XLALDestroyREAL4Sequence(channel[i]);
			LALFree(plane);
			XLALDestroyREAL8Sequence(channel_overlap);
			XLALDestroyREAL8Sequence(channel_rms);
			LALFree(channel);
			XLALDestroyTFTiling(tiling);
			XLAL_ERROR_NULL(func, XLAL_ENOMEM);
		}
	}

	/* 
	 * Initialize the structure
	 */

	plane->name[0] = '\0';
	plane->epoch.gpsSeconds = 0;
	plane->epoch.gpsNanoSeconds = 0;
	plane->deltaT = tseries_deltaT;
	plane->channels = channels;
	plane->deltaF = deltaF;
	plane->flow = flow;
	plane->channel_overlap = channel_overlap;
	plane->channel_rms = channel_rms;
	plane->channel = channel;
	plane->tiling = tiling;

	/*
	 * Success
	 */

	return plane;
}


/******** <lalVerbatim file="DestroyTFPlaneCP"> ********/
void
XLALDestroyTFPlane(
	REAL4TimeFrequencyPlane *plane
)
/******** </lalVerbatim> ********/
{
	unsigned i;

	if(plane) {
		XLALDestroyREAL8Sequence(plane->channel_overlap);
		XLALDestroyREAL8Sequence(plane->channel_rms);
		for(i = 0; i < plane->channels; i++)
			XLALDestroyREAL4Sequence(plane->channel[i]);
		LALFree(plane->channel);
		XLALDestroyTFTiling(plane->tiling);
	}
	LALFree(plane);
}
