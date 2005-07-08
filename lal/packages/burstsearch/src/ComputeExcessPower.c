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

#define TRUE 1
#define FALSE 0


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
	INT4 offset;
	INT4 nf = plane->params.freqBins;
	INT4 nt = plane->params.timeBins;
	INT4 t1;
	INT4 t2;
	INT4 k1;
	INT4 k2;
	size_t i;

	/* check on some parameter values */
	if((nf <= 0) || (nt <= 0))
		XLAL_ERROR(func, XLAL_EDOM);

	for(i = 0; i < tiling->numtiles; i++, tile++) {
		t1 = tile->tstart;
		t2 = tile->tend;
		k1 = tile->fstart;
		k2 = tile->fend;

		if((t1 < 0) || (t1 > t2) || (t2 > nt) ||
		   (k1 < 0) || (k1 > k2) || (k2 > nf))
			XLAL_ERROR(func, XLAL_EDATA);

		dof = XLALTFTileDegreesOfFreedom(tile);

		sum = 0.0;
		for(offset = t1; offset < t2; offset += (t2 - t1) / dof)
			sum += pow(XLALREAL4Sum(&plane->data[offset * nf], k1, k2 - k1), 2.0) / XLALREAL4SumSquares(norm, k1, k2 - k1);
		tile->excessPower = sum - dof;

		tile->lnalpha = XLALlnOneMinusChisqCdf(sum, dof);
		if(XLALIsREAL8FailNaN(tile->lnalpha))
			XLAL_ERROR(func, XLAL_EFUNC);
	}

	/* success */
	return(0);
}
