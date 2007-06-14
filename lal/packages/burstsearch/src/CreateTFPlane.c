/*
 * $Id$
 *
 * Copyright (C) 2007  Kipp Cannon and Flanagan, E
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


#include <lal/LALRCSID.h>


NRCSID(CREATETFPLANEC, "$Id$");


#include <math.h>
#include <lal/LALMalloc.h>
#include <lal/LALStdlib.h>
#include <lal/TFTransform.h>
#include <lal/Sequence.h>
#include <lal/XLALError.h>


/*
 * Return TRUE if a is an integer multiple of b.
 */


static int double_is_int_multiple_of(double a, double b)
{
	const double epsilon = 0;
	int n = a / b;
	return fabs(1 - n * b / a) <= epsilon;
}


/******** <lalVerbatim file="CreateTFPlaneCP"> ********/
REAL4TimeFrequencyPlane *XLALCreateTFPlane(
	/* length of time series from which TF plane will be computed */
	UINT4 tseries_length,
	/* sample rate of time series */
	REAL8 tseries_deltaT,
	/* minimum frequency to search for */
	REAL8 flow,
	/* bandwidth of TF plane */
	REAL8 bandwidth,
	/* sample at which to start the tiling */
	UINT4 tiling_start,
	/* overlap of adjacent tiles */
	REAL8 tiling_fractional_stride,
	/* largest tile's bandwidth */
	REAL8 tiling_max_bandwidth,
	/* largest tile's duration */
	REAL8 tiling_max_duration
)
/******** </lalVerbatim> ********/
{
	static const char func[] = "XLALCreateTFPlane";
	REAL8 fseries_deltaF = 1.0 / (tseries_length * tseries_deltaT);
	REAL8 deltaF = 1 / tiling_max_duration * tiling_fractional_stride;
	INT4 channels = bandwidth / deltaF;
	REAL4TimeFrequencyPlane *plane;
	REAL8Sequence *channel_overlap;
	REAL8Sequence *channel_rms;
	REAL4Sequence **channel;
	TFTiling *tiling;
	int i;

	/*
	 * Make sure that input parameters are reasonable
	 */

	if((flow < 0) ||
	   (bandwidth < 0) ||
	   (deltaF <= 0.0) ||
	   (!double_is_int_multiple_of(deltaF, fseries_deltaF)) ||
	   (tseries_deltaT <= 0.0) ||
	   (channels * deltaF != bandwidth))
		XLAL_ERROR_NULL(func, XLAL_EINVAL);

	/*
	 * Allocate memory and construct the tiling.
	 */

	plane = XLALMalloc(sizeof(*plane));
	channel_overlap = XLALCreateREAL8Sequence(channels - 1);
	channel_rms = XLALCreateREAL8Sequence(channels);
	channel = XLALMalloc(channels * sizeof(*channel));
	tiling = XLALCreateTFTiling(tiling_start, tseries_length / 2, channels, tseries_deltaT, deltaF, tiling_fractional_stride, tiling_max_bandwidth, tiling_max_duration);
	if(!plane || !channel_overlap || !channel_rms || !channel || !tiling) {
		XLALFree(plane);
		XLALDestroyREAL8Sequence(channel_overlap);
		XLALDestroyREAL8Sequence(channel_rms);
		XLALFree(channel);
		XLALDestroyTFTiling(tiling);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}
	for(i = 0; i < channels; i++) {
		channel[i] = XLALCreateREAL4Sequence(tseries_length);
		if(!channel[i]) {
			while(--i)
				XLALDestroyREAL4Sequence(channel[i]);
			XLALFree(plane);
			XLALDestroyREAL8Sequence(channel_overlap);
			XLALDestroyREAL8Sequence(channel_rms);
			XLALFree(channel);
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
void XLALDestroyTFPlane(
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
		XLALFree(plane->channel);
		XLALDestroyTFTiling(plane->tiling);
	}
	XLALFree(plane);
}
