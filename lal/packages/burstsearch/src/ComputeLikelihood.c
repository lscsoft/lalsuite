/******** <lalVerbatim file="ComputeLikelihoodCV"> ********
Author: Flanagan, E. and Cannon, K.
$Id$
********* </lalVerbatim> ********/

#include <lal/LALRCSID.h>

NRCSID (COMPUTELIKELIHOODC, "$Id$");

#include <lal/ExcessPower.h>


/******** <lalVerbatim file="ComputeLikelihoodCP"> ********/
REAL8
XLALComputeLikelihood(
	TFTile *tile
)
/******** </lalVerbatim> ********/
{
	REAL8 avglambda = 0.0;
	REAL8 dof;
	REAL8 rho4;
	int numTiles;

	for(numTiles = 0; tile; tile = tile->nextTile)
		if(tile->firstCutFlag) {
			dof = 2.0 * (tile->tend - tile->tstart + 1) * (tile->fend - tile->fstart + 1);
			rho4 = tile->excessPower * tile->excessPower;
			avglambda += dof / (rho4 * tile->alpha) * tile->weight;
		numTiles++;
		}

	/* compute the likelihood averaged over TF tiles */
	avglambda /= numTiles;

	/* return value of statistic */
	return(avglambda);
}
