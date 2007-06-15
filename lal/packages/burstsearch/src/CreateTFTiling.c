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


#include <math.h>
#include <lal/LALRCSID.h>


NRCSID(CREATETFTILINGC, "$Id$");


#include <lal/LALStdlib.h>
#include <lal/TFTransform.h>
#include <lal/XLALError.h>


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
 * Macro for looping over all tiles.  This is ugly but it ensures that the
 * initialization, increment, and terminate statements are the same in the two
 * places the loop is done.
 */


#define FOR_EACH_TILE \
	for(t_length = min_length; t_length <= max_length; t_length *= 2) \
		for(channels = min_channels; channels <= max_channels; channels *= 2) \
			for(t_start = tiling_t_start; t_start + t_length <= tiling_t_end; t_start += t_length / inv_fractional_stride) \
				for(channel_start = 0; channel_start + channels <= tiling_n_channels; channel_start += channels / inv_fractional_stride)



/******** <lalVerbatim file="CreateTFTilingCP"> ********/
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
/******** </lalVerbatim> *********/
{
	const char func[] = "XLALCreateTFTiling";
	TFTiling *tiling;
	TFTile *tile;
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
	tile = XLALMalloc(numtiles * sizeof(*tile));
	if(!tiling || !tile) {
		XLALFree(tiling);
		XLALFree(tile);
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);
	}
	tiling->tile = tile;
	tiling->numtiles = numtiles;

	/*
	 * initialize each tile
	 */

	FOR_EACH_TILE {
		/*
		 * co-ordinates
		 */

		tile->channel0 = channel_start;
		tile->channels = channels;
		tile->tstart = t_start;
		tile->tend = t_start + t_length;

		/*
		 * t_length * channels = # of pixels in tile
		 * 2 * deltaT * deltaF = degrees of freedom in 1 pixel
		 */

		tile->dof = (t_length * channels) * 2 * plane_deltaT * plane_deltaF;

		/*
		 * for safety
		 */

		tile->excessPower = XLAL_REAL8_FAIL_NAN;
		tile->confidence = XLAL_REAL8_FAIL_NAN;
		tile++;
	}

	return(tiling);
}


/******** <lalVerbatim file="DestroyTFTilingCP"> ********/
void XLALDestroyTFTiling(
	TFTiling *tiling
)
/******** </lalVerbatim> ********/
{
	if(tiling)
		XLALFree(tiling->tile);
	XLALFree(tiling);
}
