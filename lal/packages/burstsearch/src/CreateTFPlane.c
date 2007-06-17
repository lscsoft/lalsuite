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
 * Round target_length down so that an integer number of intervals of
 * length segment_length, each shifted by segment_shift with respect to the
 * interval preceding it, fits into the result.
 */


INT4 XLALOverlappedSegmentsCommensurate(
	INT4 target_length,
	INT4 segment_length,
	INT4 segment_shift
)
{
	static const char func[] = "XLALOverlappedSegmentsCommensurate";
	UINT4 segments;

	if((segment_length < 1) || (segment_shift < 1))
		XLAL_ERROR(func, XLAL_EINVAL);

	/* trivial case */
	if(target_length < segment_length)
		return 0;

	segments = (target_length - segment_length) / segment_shift;

	return segments * segment_shift + segment_length;
}


/*
 * Compute and return the timing parameters for an excess power analysis.
 *
 * Input:
 * 	window_length:
 * 		number of samples in a window used for the time-frequency
 * 		plane
 *
 * 	max_tile_length:
 * 		number of samples in the tile of longest duration
 *
 * 	fractional_tile_stride:
 * 		fractional amount by which the start of a tile is shifted
 * 		from the tile preceding it
 *
 * 	psd_length:
 * 		user's desired number of samples to use in computing a PSD
 * 		estimate
 *
 * Output:
 * 	psd_length:
 * 		actual number of samples to use in computing a PSD estimate
 * 		(rounded down to be comensurate with the windowing)
 *
 *	psd_shift:
 *		number of samples by which the start of a PSD is shifted
 *		from the start of the PSD that preceded it
 *
 *	window_shift:
 *		number of samples by which the start of a time-frequency
 *		plane window is shifted from the window preceding it
 *
 *	window_pad:
 *		how many samples at the start and end of each window are
 *		treated as padding, and will not be covered by the tiling
 */


INT4 XLALEPGetTimingParameters(
	INT4 window_length,
	INT4 max_tile_length,
	REAL8 fractional_tile_stride,
	INT4 *psd_length,
	INT4 *psd_shift,
	INT4 *window_shift,
	INT4 *window_pad
)
{
	static const char func[] = "XLALEPGetTimingParameters";

	if((max_tile_length < 1) ||
	   (fractional_tile_stride < 0))
		XLAL_ERROR(func, XLAL_EINVAL);

	/* it is assumed that the middle 1/2 of the window is analyzed */
	*window_pad = window_length / 4;
	*window_shift = window_length / 2 - (1 - fractional_tile_stride) * max_tile_length;

	*psd_length = XLALOverlappedSegmentsCommensurate(*psd_length, window_length, *window_shift);
	if(*psd_length < 0)
		XLAL_ERROR(func, XLAL_EFUNC);

	*psd_shift = *psd_length - (window_length - *window_shift);

	return 0;
}


/*
 * Return TRUE if a is an integer multiple of b
 */


static int is_int_multiple_of(int a, int b)
{
	int n;
	if(a == 0 || b == 0)
		return 0;
	n = a / b;
	return n * b == a;
}


/*
 * Return TRUE if a is an integer multiple of b.
 */


static int double_is_int_multiple_of(double a, double b)
{
	const double epsilon = 0;
	int n = a / b;
	return fabs(1 - n * b / a) <= epsilon;
}


/*
 * Macro for looping over all tiles.  This is ugly but it ensures that the
 * initialization, increment, and terminate statements are the same in the two
 * places the loop is done.
 */


#define FOR_EACH_TILE \
	for(t_length = min_length; t_length <= max_length; t_length *= 2) \
		for(channels = min_channels; channels <= max_channels; channels *= 2) \
			for(t_start = tiling_t_start; t_start + t_length <= tiling_t_end; t_start += t_length / inv_fractional_stride) \
				for(channel_start = 0; channel_start + channels <= tiling_n_channels; channel_start += channels / inv_fractional_stride)



/*
 * Allocate and initialize a tiling of the time-frequency plane.
 */


TFTiling *XLALCreateTFTiling(
	UINT4 tiling_t_start,
	UINT4 tiling_t_length,
	UINT4 tiling_n_channels,
	REAL8 plane_deltaT,
	REAL8 plane_deltaF,
	REAL8 fractional_stride,
	REAL8 max_tile_bandwidth,
	REAL8 max_tile_duration
)
{
	const char func[] = "XLALCreateTFTiling";
	TFTiling *tiling;
	TFTile *tiles;
	int numtiles;

	/*
	 * stride
	 */

	const unsigned inv_fractional_stride = 1 / fractional_stride;

	/*
	 * coordinate limits
	 */

	const unsigned tiling_t_end = tiling_t_start + tiling_t_length;

	/*
	 * coordinates of a TF tile
	 */

	unsigned channel_start;
	unsigned channels;
	unsigned t_start;
	unsigned t_length;

	/*
	 * tile size limits
	 */

	const unsigned min_length = (1 / max_tile_bandwidth) / plane_deltaT;
	const unsigned max_length = max_tile_duration / plane_deltaT;
	const unsigned min_channels = inv_fractional_stride;
	const unsigned max_channels = max_tile_bandwidth / plane_deltaF;

	/*
	 * check the tile size limits.  note that because all tile
	 * durations are integer multiples of the smallest duration, if the
	 * largest duration fits an integer number of times in the tiling
	 * length, and the smallest duration does as well, then all tile
	 * sizes in between also fit an integer number of times so there's
	 * no need to test them all.  likewise for the bandwidths.
	 */

	if((inv_fractional_stride * fractional_stride != 1) ||
	   (min_length * plane_deltaT != (1 / max_tile_bandwidth)) ||
	   !is_int_multiple_of(min_length, inv_fractional_stride) ||
	   !is_int_multiple_of(tiling_t_length, min_length) ||
	   !is_int_multiple_of(tiling_t_length, max_length) ||
	   !is_int_multiple_of(tiling_n_channels, min_channels) ||
	   !is_int_multiple_of(tiling_n_channels, max_channels)) {
	   	XLALPrintError("unable to construct time-frequency tiling from input parameters\n");
		XLAL_ERROR_NULL(func, XLAL_EINVAL);
	}

	/*
	 * count the tiles
	 */

	numtiles = 0;
	FOR_EACH_TILE {
		numtiles++;
	}

	/*
	 * allocate memory
	 */

	tiling = XLALMalloc(sizeof(*tiling));
	tiles = XLALMalloc(numtiles * sizeof(*tiles));
	if(!tiling || !tiles) {
		XLALFree(tiling);
		XLALFree(tiles);
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);
	}
	tiling->tiles = tiles;
	tiling->numtiles = numtiles;

	/*
	 * initialize each tile
	 */

	FOR_EACH_TILE {
		/*
		 * co-ordinates
		 */

		tiles->channel0 = channel_start;
		tiles->channels = channels;
		tiles->tstart = t_start;
		tiles->tend = t_start + t_length;

		/*
		 * t_length * channels = # of pixels in tile
		 * 2 * deltaT * deltaF = degrees of freedom in 1 pixel
		 */

		tiles->dof = (t_length * channels) * 2 * plane_deltaT * plane_deltaF;

		/*
		 * for safety
		 */

		tiles->excess_power = XLAL_REAL8_FAIL_NAN;
		tiles->confidence = XLAL_REAL8_FAIL_NAN;
		tiles->h_rss = XLAL_REAL8_FAIL_NAN;

		/*
		 * Next
		 */

		tiles++;
	}

	return(tiling);
}


/*
 * Free a tiling of the time-frequency plane.
 */


void XLALDestroyTFTiling(
	TFTiling *tiling
)
{
	if(tiling)
		XLALFree(tiling->tiles);
	XLALFree(tiling);
}


/*
 * Create and initialize a time-frequency plane object.
 */


REAL4TimeFrequencyPlane *XLALCreateTFPlane(
	/* length of time series from which TF plane will be computed */
	UINT4 tseries_length,
	/* sample at which to start the tiling */
	UINT4 tiling_start,
	/* length of tiling */
	UINT4 tiling_length,
	/* sample rate of time series */
	REAL8 tseries_deltaT,
	/* minimum frequency to search for */
	REAL8 flow,
	/* bandwidth of TF plane */
	REAL8 bandwidth,
	/* overlap of adjacent tiles */
	REAL8 tiling_fractional_stride,
	/* largest tile's bandwidth */
	REAL8 tiling_max_bandwidth,
	/* largest tile's duration */
	REAL8 tiling_max_duration
)
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
	   (bandwidth <= 0) ||
	   (deltaF <= 0) ||
	   (tiling_start + tiling_length > tseries_length) ||
	   (!double_is_int_multiple_of(deltaF, fseries_deltaF)) ||
	   (tseries_deltaT <= 0) ||
	   (channels * deltaF != bandwidth))
		XLAL_ERROR_NULL(func, XLAL_EINVAL);

	/*
	 * Allocate memory and construct the tiling.
	 */

	plane = XLALMalloc(sizeof(*plane));
	channel_overlap = XLALCreateREAL8Sequence(channels - 1);
	channel_rms = XLALCreateREAL8Sequence(channels);
	channel = XLALMalloc(channels * sizeof(*channel));
	tiling = XLALCreateTFTiling(tiling_start, tiling_length, channels, tseries_deltaT, deltaF, tiling_fractional_stride, tiling_max_bandwidth, tiling_max_duration);
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


/*
 * Free a time-frequency plane object.
 */


void XLALDestroyTFPlane(
	REAL4TimeFrequencyPlane *plane
)
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
