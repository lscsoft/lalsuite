/******** <lalVerbatim file="CreateTFTilingCV"> ********
Author: Eanna Flanagan, and Cannon, K.
$Id$
********* </lalVerbatim> **********/

#include <math.h>
#include <lal/LALRCSID.h>

NRCSID(CREATETFTILINGC, "$Id$");

#include <lal/ExcessPower.h>
#include <lal/LALStdlib.h>
#include <lal/TFTransform.h>
#include <lal/XLALError.h>


/*
 * The number of degrees of freedom in a tile.
 */

REAL8 XLALTFTileDegreesOfFreedom(const TFTile *tile)
{
	return((2 * tile->tbins * tile->fbins) * tile->deltaT * tile->deltaF);
}


/*
 * Macro for looping over all tiles.  This is ugly but it ensures that the
 * initialization, increment, and terminate statements are the same in the two
 * places the loop is done.
 */

#define FOR_EACH_TILE \
	for(tbins = 2; tbins <= max_tbins; tbins *= 2) \
		for(fbins = 1 / (tbins * plane->deltaT * plane->deltaF); fbins <= max_fbins; fbins *= 2) \
			for(tstart = 0; tstart + tbins <= plane->timeBins; tstart += tbins / inv_fractional_stride) \
				for(fstart = 0; fstart + fbins <= plane->freqBins; fstart += fbins / inv_fractional_stride)



/******** <lalVerbatim file="CreateTFTilingCP"> ********/
TFTiling *XLALCreateTFTiling(
	const REAL4TimeFrequencyPlane *plane,
	INT4 inv_fractional_stride,
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

	/* tile size limits */
	int max_tbins, max_fbins;
	int maxDOF;

	/* coordinates of a TF tile */
	INT4 fstart;
	INT4 fbins;
	INT4 tstart;
	INT4 tbins;

	if(!plane)
		XLAL_ERROR_NULL(func, XLAL_EFAULT);

	/* determine the tile size limits */
	/* FIXME: should the two conditionals in fact be errors? */
	max_tbins = maxTileDuration / plane->deltaT;
	if(plane->timeBins < max_tbins)
		max_tbins = plane->timeBins;
	max_fbins = maxTileBandwidth / plane->deltaF;
	if(plane->freqBins < max_fbins)
		max_fbins = plane->freqBins;
	maxDOF = (2 * max_tbins * max_fbins) * plane->deltaT * plane->deltaF;

	/* Count the tiles */
	numtiles = 0;
	FOR_EACH_TILE {
		numtiles++;
	}
	if(!numtiles)
		/* can't fit any tiles into the TF plane! */
		XLAL_ERROR_NULL(func, XLAL_EINVAL);

	/* allocate memory */
	tiling = LALMalloc(sizeof(*tiling));
	tile = LALMalloc(numtiles * sizeof(*tile));
	weight = LALCalloc(maxDOF + 1, sizeof(*weight));
	if(!tiling || !tile || !weight) {
		LALFree(tiling);
		LALFree(tile);
		LALFree(weight);
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);
	}
	tiling->numtiles = numtiles;
	tiling->tile = tile;

	/* initialize each tile */
	FOR_EACH_TILE {
		tile->fstart = fstart;
		tile->fbins = fbins;
		tile->tstart = tstart;
		tile->tbins = tbins;
		tile->flow = plane->flow;
		tile->deltaT = plane->deltaT;
		tile->deltaF = plane->deltaF;
		tile->excessPower = XLAL_REAL8_FAIL_NAN;
		tile->lnalpha = XLAL_REAL8_FAIL_NAN;
		weight[(int) XLALTFTileDegreesOfFreedom(tile)]++;
		tile++;
	}

	/* determine the weighting for each tile */
	for(tile = tiling->tile; numtiles; numtiles--, tile++)
		tile->lnweight = log(weight[(int) XLALTFTileDegreesOfFreedom(tile)]);
	LALFree(weight);

	return(tiling);
}


/******** <lalVerbatim file="DestroyTFTilingCP"> ********/
void
XLALDestroyTFTiling(
	TFTiling *tiling
)
/******** </lalVerbatim> ********/
{
	if(tiling)
		LALFree(tiling->tile);
	LALFree(tiling);
}
