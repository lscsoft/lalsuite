/******** <lalVerbatim file="CreateTFTilingCV"> ********
Author: Eanna Flanagan
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
		tile->whichPlane = 0;
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
TFTiling *XLALCreateTFTiling(
	CreateTFTilingIn *input,
	TFPlaneParams *planeParams
)
/******** </lalVerbatim> *********/
{
	static const char *func = "XLALCreateTFTiling";
	TFTiling *tfTiling;
	COMPLEX8TimeFrequencyPlane *plane;
	TFTile **tile;

	REAL8 fhigh;		/* max frequency of the TF plane */
	REAL8 flow;		/* min frequency of the TF plane */
	REAL8 timeDuration;	/* time duration of the plane    */
	LALWindowParams winParams;
	/* coordinates of a TF tile */
	INT4 fstart;
	INT4 deltaf;
	INT4 tstart;
	INT4 deltat;

	/* lowest frequency to be used in the plane. */
	flow = planeParams->flow;

	/* highest frequency to be used in the plane  */
	fhigh = planeParams->fhigh;

	/* time duration to be searched for, this will determine the no. of
	 * time bins in the plane. */
	timeDuration = planeParams->timeDuration;

	/* set TF plane params */
	planeParams->timeBins = timeDuration / planeParams->deltaT;
	planeParams->freqBins = (fhigh - flow) / planeParams->deltaF;
	if(planeParams->freqBins < input->minFreqBins) {
		fprintf(stderr, "no of freqbins is less than the minimum allowed\n");
		XLAL_ERROR_NULL(func, XLAL_EINVAL);
	}

	/* set DFT params */
	winParams.type = Rectangular;
	winParams.length = planeParams->timeBins;

	/* Allocate memory for tfTiling  */
	tfTiling = LALMalloc(sizeof(*tfTiling));
	plane = XLALCreateTFPlane(planeParams);
	if(!tfTiling || !plane) {
		LALFree(tfTiling);
		XLALDestroyTFPlane(tfTiling->tfp);
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);
	}
	tfTiling->tfp = plane;

	/* set tiling parameters */
	tfTiling->numPlanes = 1;
	tfTiling->planesComputed = FALSE;
	tfTiling->excessPowerComputed = FALSE;
	tfTiling->tilesSorted = FALSE;

	/* 
	 *  construct the linked list of Time Frequency Tiles
	 */

	tfTiling->numTiles = 0;
	tile = &tfTiling->firstTile;
	*tile = NULL;

	/* deltat should correspond to no. of time bins for a particular
	 * tile */
	for(deltat = input->minTimeBins; (deltat <= plane->params->timeBins) && (deltat * plane->params->deltaT <= input->maxTileDuration); deltat *= 2)
		for(tstart = 0; tstart <= plane->params->timeBins - deltat; tstart += deltat / input->overlapFactor)
			/* deltaf is set by the deltat, requiring that the
			 * TF vol is >= 1 */
			for(deltaf = 1 / (deltat * plane->params->deltaT * plane->params->deltaF); (deltaf <= plane->params->freqBins) && (deltaf * plane->params->deltaF <= input->maxTileBand); deltaf *= 2)
				for(fstart = 0; fstart <= plane->params->freqBins - deltaf; fstart += deltaf / input->overlapFactor) {
					*tile = CreateTile(fstart, deltaf, tstart, deltat, plane->params->deltaT, plane->params->deltaF);
					if(!*tile)
						XLAL_ERROR_NULL(func, XLAL_ENOMEM);
					tfTiling->numTiles++;
					tile = &(*tile)->nextTile;
				}

	return(tfTiling);
}
