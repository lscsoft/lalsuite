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

#define TRUE 1


/*
 * The number of degrees of freedom in a tile.
 */

REAL8 XLALTFTileDegreesOfFreedom(TFTile *tile)
{
	return(2.0 * (tile->tend - tile->tstart) * tile->deltaT * (tile->fend - tile->fstart) * tile->deltaF);
}


/******** <lalVerbatim file="CreateTFTilingCP"> ********/
TFTiling *XLALCreateTFTiling(
	const CreateTFTilingIn *input,
	const TFPlaneParams *planeparams
)
/******** </lalVerbatim> *********/
{
	static const char *func = "XLALCreateTFTiling";
	TFTiling *tiling;
	int numtiles;
	TFTile *tile;
	const int maxDOF = (int) (2.0 * planeparams->timeBins * planeparams->deltaT * planeparams->freqBins * planeparams->deltaF);
	INT4 *weight;
	size_t i;

	/* coordinates of a TF tile */
	INT4 fstart;
	INT4 deltaf;
	INT4 tstart;
	INT4 deltat;

	/* Count the tiles */
	/* FIXME: this is stupid... figure out how to compute directly */
	numtiles = 0;
	for(deltat = input->minTimeBins; (deltat <= planeparams->timeBins) && (deltat * planeparams->deltaT <= input->maxTileDuration); deltat *= 2)
		for(tstart = 0; tstart <= planeparams->timeBins - deltat; tstart += deltat / input->overlapFactor)
			for(deltaf = 1 / (deltat * planeparams->deltaT * planeparams->deltaF); (deltaf <= planeparams->freqBins) && (deltaf * planeparams->deltaF <= input->maxTileBand); deltaf *= 2)
				for(fstart = 0; fstart <= planeparams->freqBins - deltaf; fstart += deltaf / input->overlapFactor)
					numtiles++;

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

	/* 
	 * Initialize each tile one-by-one.
	 */

	/* deltat should correspond to no. of time bins for a particular
	 * tile */
	for(deltat = input->minTimeBins; (deltat <= planeparams->timeBins) && (deltat * planeparams->deltaT <= input->maxTileDuration); deltat *= 2)
		for(tstart = 0; tstart <= planeparams->timeBins - deltat; tstart += deltat / input->overlapFactor)
			/* deltaf is set by the deltat, requiring that the
			 * TF vol is >= 1 */
			for(deltaf = 1 / (deltat * planeparams->deltaT * planeparams->deltaF); (deltaf <= planeparams->freqBins) && (deltaf * planeparams->deltaF <= input->maxTileBand); deltaf *= 2)
				for(fstart = 0; fstart <= planeparams->freqBins - deltaf; fstart += deltaf / input->overlapFactor) {
					tile->fstart = fstart;
					tile->fend = fstart + deltaf;
					tile->tstart = tstart;
					tile->tend = tstart + deltat;
					tile->deltaT = planeparams->deltaT;
					tile->deltaF = planeparams->deltaF;
					tile->excessPower = XLAL_REAL8_FAIL_NAN;
					tile->lnalpha = XLAL_REAL8_FAIL_NAN;
					tile->PassFirstCut = TRUE;

					weight[(int) XLALTFTileDegreesOfFreedom(tile)]++;

					tile++;
				}

	/*
	 * Determine the weighting for each tile.
	 */

	for(i = 0, tile = tiling->tile; i < tiling->numtiles; i++, tile++)
		tile->lnweight = log(weight[(int) XLALTFTileDegreesOfFreedom(tile)]);
	LALFree(weight);

	return(tiling);
}
