/******** <lalVerbatim file="ComputeLikelihoodCV"> ********
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/

#include <lal/LALRCSID.h>

NRCSID (COMPUTELIKELIHOODC, "$Id$");

#include <lal/ExcessPower.h>
#include <lal/XLALError.h>


/******** <lalVerbatim file="ComputeLikelihoodCP"> ********/
REAL8
XLALComputeLikelihood(
	TFTiling *tfTiling
)
/******** </lalVerbatim> ********/
{
	static const char *func = "XLALComputeLikelihood";
	REAL8 avglambda = 0.0;
	REAL8 dof;
	REAL8 rho4;
	TFTile *tile;

	/* make sure XLALComputeExcessPower() has been already called */
	if(!tfTiling->excessPowerComputed)
		XLAL_ERROR_REAL8(func, XLAL_EDATA);

	for(tile = tfTiling->firstTile; tile; tile = tile->nextTile)
		if(tile->firstCutFlag) {
			dof = 2.0 * (tile->tend - tile->tstart + 1) * (tile->fend - tile->fstart + 1);
			rho4 = tile->excessPower * tile->excessPower;
			avglambda += dof / (rho4 * tile->alpha) * tile->weight;
		}

	/* compute the likelihood averaged over TF tiles */
	avglambda /= tfTiling->numTiles;

	/* return value of statistic */
	return(avglambda);
}
