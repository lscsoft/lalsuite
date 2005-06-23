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
	TFTiling *tiling
)
/******** </lalVerbatim> ********/
{
	REAL8 avglambda = 0.0;
	REAL8 dof;
	REAL8 rho4;
	TFTile *tile;
	size_t i;

	for(i = 0, tile = tiling->tile; i < tiling->numtiles; i++, tile++)
		if(tile->firstCutFlag) {
			/* FIXME: should this be the XLAL function? */
			dof = 2.0 * (tile->tend - tile->tstart + 1) * (tile->fend - tile->fstart + 1);
			rho4 = tile->excessPower * tile->excessPower;
			avglambda += dof / (rho4 * tile->alpha) * tile->weight;
		}

	/* compute the likelihood averaged over TF tiles */
	avglambda /= tiling->numtiles;

	/* return value of statistic */
	return(avglambda);
}
