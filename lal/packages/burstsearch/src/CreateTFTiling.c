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
 * The number of degrees of freedom in a tile.
 */


REAL8 XLALTFTileDegreesOfFreedom(const TFTile *tile)
{
	return((2 * (tile->tend - tile->tstart) * tile->channels) * tile->deltaT * tile->deltaF);
}


/*
 * Macro for looping over all tiles.  This is ugly but it ensures that the
 * initialization, increment, and terminate statements are the same in the two
 * places the loop is done.
 */


#define FOR_EACH_TILE \
	for(tbins = min_tbins; tbins <= max_tbins; tbins *= 2) \
		for(channels = 1 / (tbins * plane_deltaT * plane_deltaF); channels <= max_channels; channels *= 2) \
			for(tstart = tiling_tstart; tstart + tbins <= tmax; tstart += tbins / inv_fractional_stride) \
				for(channel0 = 0; channel0 + channels <= plane_num_channels; channel0 += channels / inv_fractional_stride)



/******** <lalVerbatim file="CreateTFTilingCP"> ********/
TFTiling *XLALCreateTFTiling(
	UINT4 plane_length,
	REAL8 plane_deltaT,
	REAL8 plane_flow,
	REAL8 plane_deltaF,
	UINT4 plane_num_channels,
	UINT4 tiling_tstart,
	UINT4 inv_fractional_stride,
	REAL8 maxTileBandwidth,
	REAL8 maxTileDuration
)
/******** </lalVerbatim> *********/
{
	const char func[] = "XLALCreateTFTiling";
	TFTiling *tiling;
	TFTile *tile;
	int *weight;
	int numtiles;

	/* coordinates of a TF tile */
	unsigned channel0;
	unsigned channels;
	unsigned tstart;
	unsigned tbins;

	/* coordinate limits */
	const unsigned tmax = plane_length - tiling_tstart;

	/* tile size limits */
	const unsigned min_tbins = 1.0 / maxTileBandwidth / plane_deltaT;
	const unsigned max_tbins = maxTileDuration / plane_deltaT;
	const unsigned min_channels = 1 / (max_tbins * plane_deltaT * plane_deltaF);
	const unsigned max_channels = maxTileBandwidth / plane_deltaF;
	const int maxDOF = (2 * max_tbins * max_channels) * plane_deltaT * plane_deltaF;

	/* check the tile size limits */
	if((min_tbins < inv_fractional_stride) ||
	   (min_channels < inv_fractional_stride) ||
	   (tmax < tiling_tstart + max_tbins) ||
	   (plane_num_channels < max_channels))
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
	weight = XLALCalloc(maxDOF + 1, sizeof(*weight));
	if(!tiling || !tile || !weight) {
		XLALFree(tiling);
		XLALFree(tile);
		XLALFree(weight);
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
		tile->flow = plane_flow;
		tile->deltaT = plane_deltaT;
		tile->deltaF = plane_deltaF;
		tile->excessPower = XLAL_REAL8_FAIL_NAN;
		tile->confidence = XLAL_REAL8_FAIL_NAN;
		weight[(int) XLALTFTileDegreesOfFreedom(tile)]++;
		tile++;
	}

	/* determine the weighting for each tile */
	for(tile = tiling->tile; numtiles; numtiles--, tile++)
		tile->lnweight = log(weight[(int) XLALTFTileDegreesOfFreedom(tile)]);
	XLALFree(weight);

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
