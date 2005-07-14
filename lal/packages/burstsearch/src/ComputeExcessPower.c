/******** <lalVerbatim file="ComputeExcessPowerCV"> ********
Author: Flanagan, E. and Cannon, K.
$Id$
********* </lalVerbatim> ********/

#include <math.h>
#include <lal/LALRCSID.h>

NRCSID (COMPUTEEXCESSPOWERC, "$Id$");

#include <lal/ExcessPower.h>
#include <lal/Sequence.h>
#include <lal/Thresholds.h>
#include <lal/XLALError.h>


/******** <lalVerbatim file="ComputeExcessPowerCP"> ********/
int
XLALComputeExcessPower(
	TFTiling *tiling,
	const REAL4TimeFrequencyPlane *plane,
	const REAL4 *norm
)
/******** </lalVerbatim> ********/
{
	static const char *func = "XLALComputeExcessPower";
	TFTile *tile = tiling->tile;
	REAL8 sum;
	REAL8 dof;
	INT4 bin;
	INT4 nf = plane->params.freqBins;
	INT4 nt = plane->params.timeBins;
	size_t i;

	/* check on some parameter values */
	if((nf <= 0) || (nt <= 0))
		XLAL_ERROR(func, XLAL_EDOM);

	for(i = 0; i < tiling->numtiles; i++, tile++) {
		if((tile->tstart < 0) || (tile->tstart + tile->tbins > nt) ||
		   (tile->fstart < 0) || (tile->fstart + tile->fbins > nf))
			XLAL_ERROR(func, XLAL_EDATA);

		dof = XLALTFTileDegreesOfFreedom(tile);

		sum = 0.0;
		for(bin = 0; bin < tile->tbins; bin += tile->tbins / dof)
			sum += pow(XLALREAL4Sum(&plane->data[(tile->tstart + bin) * nf], tile->fstart, tile->fbins), 2.0) / XLALREAL4SumSquares(norm, tile->fstart, tile->fbins);
		tile->excessPower = sum - dof;

		tile->lnalpha = XLALlnOneMinusChisqCdf(sum, dof);
		if(XLALIsREAL8FailNaN(tile->lnalpha))
			XLAL_ERROR(func, XLAL_EFUNC);
	}

	/* success */
	return(0);
}
