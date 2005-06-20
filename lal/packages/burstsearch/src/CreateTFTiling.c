/******** <lalVerbatim file="CreateTFTilingCV"> ********
Author: Eanna Flanagan, and Cannon, K.
$Id$
********* </lalVerbatim> **********/

#include <lal/LALRCSID.h>

NRCSID(CREATETFTILINGC, "$Id$");

#include <lal/ExcessPower.h>
#include <lal/LALStdlib.h>
#include <lal/TFTransform.h>
#include <lal/XLALError.h>

#define FALSE 0

static INT4 pow1(INT4 a, INT4 b)
{
	/* returns a^b */
	INT4 t = 1;
	INT4 i;
	for(i = 0; i < b; i++)
		t *= a;
	return(t);
}


static TFTile *CreateTile(INT4 f, INT4 df, INT4 t, INT4 dt, REAL8 deltaT, REAL8 deltaF)
{
	TFTile *tile;

	tile = LALMalloc(sizeof(*tile));
	if(tile) {
		tile->fstart = f;
		tile->fend = f + df;
		tile->tstart = t;
		tile->tend = t + dt;
		tile->deltaT = deltaT;
		tile->deltaF = deltaF;
		tile->excessPower = 0.0;
		tile->alpha = 0.0;
		tile->weight = 1.0;
		tile->firstCutFlag = FALSE;
		tile->nextTile = NULL;
	}

	return(tile);
}


/******** <lalVerbatim file="CreateTFTilingCP"> ********/
TFTile *XLALCreateTFTiling(
	const CreateTFTilingIn *input,
	const TFPlaneParams *planeparams
)
/******** </lalVerbatim> *********/
{
	static const char *func = "XLALCreateTFTiling";
	TFTile *list;
	TFTile **tile;

	LALWindowParams winParams;
	/* coordinates of a TF tile */
	INT4 fstart;
	INT4 deltaf;
	INT4 tstart;
	INT4 deltat;

	/* set DFT params */
	winParams.type = Rectangular;
	winParams.length = planeparams->timeBins;

	/* 
	 *  construct the linked list of Time Frequency Tiles
	 */

	tile = &list;
	*tile = NULL;

	/* deltat should correspond to no. of time bins for a particular
	 * tile */
	for(deltat = input->minTimeBins; (deltat <= planeparams->timeBins) && (deltat * planeparams->deltaT <= input->maxTileDuration); deltat *= 2)
		for(tstart = 0; tstart <= planeparams->timeBins - deltat; tstart += deltat / input->overlapFactor)
			/* deltaf is set by the deltat, requiring that the
			 * TF vol is >= 1 */
			for(deltaf = 1 / (deltat * planeparams->deltaT * planeparams->deltaF); (deltaf <= planeparams->freqBins) && (deltaf * planeparams->deltaF <= input->maxTileBand); deltaf *= 2)
				for(fstart = 0; fstart <= planeparams->freqBins - deltaf; fstart += deltaf / input->overlapFactor) {
					*tile = CreateTile(fstart, deltaf, tstart, deltat, planeparams->deltaT, planeparams->deltaF);
					if(!*tile)
						XLAL_ERROR_NULL(func, XLAL_ENOMEM);
					tile = &(*tile)->nextTile;
				}

	return(list);
}
