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

	/*
	 * check input
	 */

	if(segment_length < 1) {
		XLALPrintError("segment_length < 1");
		XLAL_ERROR(func, XLAL_EINVAL);
	}
	if(segment_shift < 1) {
		XLALPrintError("segment_shift < 1");
		XLAL_ERROR(func, XLAL_EINVAL);
	}

	/*
	 * trivial case
	 */

	if(target_length < segment_length)
		return 0;

	/*
	 * do the arithmetic
	 */

	segments = (target_length - segment_length) / segment_shift;

	return segments * segment_shift + segment_length;
}


/*
 * Compute and return the timing parameters for an excess power analysis.
 * Pass NULL for any optional pointer to not compute and return that
 * parameter.
 *
 * Input:
 * 	window_length:
 * 		number of samples in a window used for the time-frequency
 * 		plane
 *
 * 	max_tile_length:
 * 		number of samples in the tile of longest duration
 *
 * 	fractional_tile_shift:
 * 		number of samples by which the start of the longest tile is
 * 		shifted from the start of the tile preceding it, as a
 * 		fraction of its length
 *
 * 	psd_length (optional, required for psd_shift):
 * 		user's desired number of samples to use in computing a PSD
 * 		estimate
 *
 * Output:
 * 	psd_length (optional):
 * 		actual number of samples to use in computing a PSD estimate
 * 		(rounded down to be comensurate with the windowing)
 *
 *	psd_shift (optional):
 *		number of samples by which the start of a PSD is shifted
 *		from the start of the PSD that preceded it
 *
 *	window_shift (optional):
 *		number of samples by which the start of a time-frequency
 *		plane window is shifted from the window preceding it
 *
 *	window_pad (optional):
 *		how many samples at the start and end of each window are
 *		treated as padding, and will not be covered by the tiling
 *
 *	tiling_length (options):
 *		how many samples will be covered by the tiling
 */


INT4 XLALEPGetTimingParameters(
	INT4 window_length,
	INT4 max_tile_length,
	REAL8 fractional_tile_shift,
	INT4 *psd_length,
	INT4 *psd_shift,
	INT4 *window_shift,
	INT4 *window_pad,
	INT4 *tiling_length
)
{
	static const char func[] = "XLALEPGetTimingParameters";
	int max_tile_shift = fractional_tile_shift * max_tile_length;
	int wpad;
	int tlength;
	int wshift;

	/*
	 * check input parameters
	 */

	if(window_length % 4 != 0) {
		XLALPrintError("window_length is not a multiple of 4");
		XLAL_ERROR(func, XLAL_EINVAL);
	}
	if(max_tile_length < 1) {
		XLALPrintError("max_tile_length < 1");
		XLAL_ERROR(func, XLAL_EINVAL);
	}
	if(fractional_tile_shift <= 0) {
		XLALPrintError("fractional_tile_shift <= 0");
		XLAL_ERROR(func, XLAL_EINVAL);
	}
	if(fmod(fractional_tile_shift * max_tile_length, 1) != 0) {
		XLALPrintError("fractional_tile_shift * max_tile_length not an integer");
		XLAL_ERROR(func, XLAL_EINVAL);
	}
	if(max_tile_shift < 1) {
		XLALPrintError("fractional_tile_shift * max_tile_length < 1");
		XLAL_ERROR(func, XLAL_EINVAL);
	}

	/*
	 * discard first and last 1/4 of the window
	 */

	wpad = window_length / 4;

	/*
	 * tiling covers the remainder, rounded down to fit an integer
	 * number of tiles
	 */

	tlength = window_length - 2 * wpad;
	tlength = XLALOverlappedSegmentsCommensurate(tlength, max_tile_length, max_tile_shift);
	if(tlength < 1) {
		XLALPrintError("tiling_length < 1");
		XLAL_ERROR(func, XLAL_EINVAL);
	}
	if(tiling_length)
		*tiling_length = tlength;

	/*
	 * now re-compute window_pad from tiling_length
	 */

	wpad = (window_length - tlength) / 2;
	if(tlength + 2 * wpad != window_length) {
		XLALPrintError("cannot find window parameters consistent with tiling parameters");
		XLAL_ERROR(func, XLAL_EINVAL);
	}
	if(window_pad)
		*window_pad = wpad;

	/*
	 * adjacent tilings overlap so that their largest tiles overlap the
	 * same as within each tiling
	 */

	wshift = tlength - (max_tile_length - max_tile_shift);
	if(wshift < 1) {
		XLALPrintError("window_shift < 1");
		XLAL_ERROR(func, XLAL_EINVAL);
	}
	if(window_shift)
		*window_shift = wshift;

	/*
	 * compute the adjusted PSD length if desired
	 */

	if(psd_length) {
		*psd_length = XLALOverlappedSegmentsCommensurate(*psd_length, window_length, wshift);
		if(*psd_length < 0)
			XLAL_ERROR(func, XLAL_EFUNC);

		if(psd_shift) {
			*psd_shift = *psd_length - (window_length - wshift);
			if(*psd_shift < 1) {
				XLALPrintError("psd_shift < 1");
				XLAL_ERROR(func, XLAL_EINVAL);
			}
		}
	} else if(psd_shift) {
		/* for safety */
		*psd_shift = -1;
		/* can't compute psd_shift without psd_length input */
		XLAL_ERROR(func, XLAL_EFAULT);
	}

	return 0;
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
	   (min_length % inv_fractional_stride != 0) ||
	   (tiling_t_length % min_length != 0) ||
	   (tiling_t_length % max_length != 0) ||
	   (tiling_n_channels % min_channels != 0) ||
	   (tiling_n_channels % max_channels != 0)) {
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
	/* sample on which tiling starts */
	int tiling_start;
	/* length of tiling */
	int tiling_length;
	/* window shift */
	int window_shift;
	/* resolution of FT of input time series */
	double fseries_deltaF = 1.0 / (tseries_length * tseries_deltaT);
	/* time-frequency plane's channel spacing */
	double deltaF = 1 / tiling_max_duration * tiling_fractional_stride;
	/* total number of channels */
	int channels = bandwidth / deltaF;
	REAL4TimeFrequencyPlane *plane;
	REAL8Sequence *channel_overlap;
	REAL8Sequence *channel_rms;
	REAL4Sequence **channel;
	REAL4Window *tukey;
	TFTiling *tiling;
	int i;

	/*
	 * Make sure that input parameters are reasonable
	 */

	if((flow < 0) ||
	   (bandwidth <= 0) ||
	   (deltaF <= 0) ||
	   (fmod(tiling_max_duration, tseries_deltaT) != 0) ||
	   (fmod(deltaF, fseries_deltaF) != 0) ||
	   (tseries_deltaT <= 0) ||
	   (channels * deltaF != bandwidth))
		XLAL_ERROR_NULL(func, XLAL_EINVAL);

	/*
	 * Compute timeing parameters
	 */

	if(XLALEPGetTimingParameters(tseries_length, tiling_max_duration / tseries_deltaT, tiling_fractional_stride, NULL, NULL, &window_shift, &tiling_start, &tiling_length) < 0)
		XLAL_ERROR_NULL(func, XLAL_EFUNC);

	/*
	 * Allocate memory and construct the tiling.
	 */

	plane = XLALMalloc(sizeof(*plane));
	channel_overlap = XLALCreateREAL8Sequence(channels - 1);
	channel_rms = XLALCreateREAL8Sequence(channels);
	channel = XLALMalloc(channels * sizeof(*channel));
	tiling = XLALCreateTFTiling(tiling_start, tiling_length, channels, tseries_deltaT, deltaF, tiling_fractional_stride, tiling_max_bandwidth, tiling_max_duration);
	tukey = XLALCreateTukeyREAL4Window(tseries_length, (tseries_length - tiling_length) / (double) tseries_length);
	if(!plane || !channel_overlap || !channel_rms || !channel || !tiling || !tukey) {
		XLALFree(plane);
		XLALDestroyREAL8Sequence(channel_overlap);
		XLALDestroyREAL8Sequence(channel_rms);
		XLALFree(channel);
		XLALDestroyTFTiling(tiling);
		XLALDestroyREAL4Window(tukey);
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
			XLALDestroyREAL4Window(tukey);
			XLAL_ERROR_NULL(func, XLAL_ENOMEM);
		}
	}

	/*
	 * Adjust the Tukey window's normalization so that it is
	 * RMS-preserving in the flat portion.  This is done by setting its
	 * normalization to be equal to that of a rectangular window
	 * (pretend the tapers aren't present).
	 */

	tukey->sumofsquares = tukey->data->length;

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
	plane->tukey = tukey;
	plane->window_shift = window_shift;

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
		XLALDestroyREAL4Window(plane->tukey);
	}
	XLALFree(plane);
}
