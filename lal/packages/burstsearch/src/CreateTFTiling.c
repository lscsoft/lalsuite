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
 * Macro for looping over all tiles.  This is ugly but it ensures that the
 * initialization, increment, and terminate statements are the same in the two
 * places the loop is done.
 */


#define FOR_EACH_TILE \
	for(tbins = min_tbins; tbins <= max_tbins; tbins *= 2) \
		for(channels = 1 / (tbins * plane_deltaT * plane_deltaF); channels <= max_channels; channels *= 2) \
			for(tstart = tiling_t_start; tstart + tbins <= tmax; tstart += tbins / inv_fractional_stride) \
				for(channel0 = 0; channel0 + channels <= tiling_n_channels; channel0 += channels / inv_fractional_stride)



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

	/* coordinates of a TF tile */
	unsigned channel0;
	unsigned channels;
	unsigned tstart;
	unsigned tbins;

	/* stride */
	const unsigned inv_fractional_stride = 1 / fractional_stride;

	/* coordinate limits */
	const unsigned tmax = tiling_t_start + tiling_t_length;

	/* tile size limits */
	const unsigned min_tbins = (1 / max_tile_bandwidth) / plane_deltaT;
	const unsigned max_tbins = max_tile_duration / plane_deltaT;
	const unsigned min_channels = inv_fractional_stride;
	const unsigned max_channels = max_tile_bandwidth / plane_deltaF;

	/* FIXME:  move tiling parameter checks from lalapps_power into
	 * this function, so that any code that uses this function will
	 * have its input validated */

	/* check the tile size limits */
	if((min_tbins < inv_fractional_stride) ||
	   (tmax < tiling_t_start + max_tbins) ||
	   (max_tbins > tiling_t_length))
		XLAL_ERROR_NULL(func, XLAL_EINVAL);

	/* Count the tiles */
	numtiles = 0;
	FOR_EACH_TILE {
		numtiles++;
	}
	if(!numtiles)
		/* can't fit any tiles into the TF plane! */
		XLAL_ERROR_NULL(func, XLAL_EINVAL);

	/* allocate memory */
	tiling = XLALMalloc(sizeof(*tiling));
	tile = XLALMalloc(numtiles * sizeof(*tile));
	if(!tiling || !tile) {
		XLALFree(tiling);
		XLALFree(tile);
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);
	}
	tiling->tile = tile;
	tiling->numtiles = numtiles;

	/* initialize each tile */
	FOR_EACH_TILE {
		tile->channel0 = channel0;
		tile->channels = channels;
		tile->tstart = tstart;
		tile->tend = tstart + tbins;
		tile->dof = ((tile->tend - tile->tstart) * tile->channels) * 2 * plane_deltaT * plane_deltaF;
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
