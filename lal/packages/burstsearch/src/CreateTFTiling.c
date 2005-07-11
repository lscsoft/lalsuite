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
		for(fbins = 1 / (tbins * planeparams->deltaT * planeparams->deltaF); fbins <= max_fbins; fbins *= 2) \
			for(tstart = 0; tstart + tbins <= planeparams->timeBins; tstart += tbins / input->overlapFactor) \
				for(fstart = 0; fstart + fbins <= planeparams->freqBins; fstart += fbins / input->overlapFactor)



/******** <lalVerbatim file="CreateTFTilingCP"> ********/
TFTiling *XLALCreateTFTiling(
	const CreateTFTilingIn *input,
	const TFPlaneParams *planeparams
)
/******** </lalVerbatim> *********/
{
	static const char *func = "XLALCreateTFTiling";
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

	/* determine the tile size limits */
	max_tbins = input->maxTileDuration / planeparams->deltaT;
	if(planeparams->timeBins < max_tbins)
		max_tbins = planeparams->timeBins;
	max_fbins = input->maxTileBandwidth / planeparams->deltaF;
	if(planeparams->freqBins < max_fbins)
		max_fbins = planeparams->freqBins;
	maxDOF = (2 * max_tbins * max_fbins) * planeparams->deltaT * planeparams->deltaF;

	/* Count the tiles */
	numtiles = 0;
	FOR_EACH_TILE {
		numtiles++;
	}

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
		tile->deltaT = planeparams->deltaT;
		tile->deltaF = planeparams->deltaF;
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
