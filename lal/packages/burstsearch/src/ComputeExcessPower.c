/******** <lalVerbatim file="ComputeExcessPowerCV"> ********
Author: Flanagan, E
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

static REAL8 square_complex8_sum(const COMPLEX8 *vec, int start, int length)
{
	COMPLEX8 sum = { 0.0, 0.0 };

	for(vec += start; length-- > 0; vec++) {
		sum.re += (*vec).re;
		sum.im += (*vec).im;
	}

	return(sum.re * sum.re + sum.im * sum.im);
}


/******** <lalVerbatim file="ComputeExcessPowerCP"> ********/
int
XLALComputeExcessPower(
	TFTiling *tfTiling,
	const ComputeExcessPowerIn *input,
	const REAL4 *norm
)
/******** </lalVerbatim> ********/
{
	static const char *func = "XLALComputeExcessPower";
	COMPLEX8TimeFrequencyPlane *plane = tfTiling->tfp;
	TFTile *tile;
	REAL8 sum;
	REAL8 dof;
	REAL8 numsigma;
	INT4 offset;
	INT4 nf = plane->params->freqBins;
	INT4 nt = plane->params->timeBins;
	INT4 t1;
	INT4 t2;
	INT4 k1;
	INT4 k2;

	/* check on some parameter values */
	if((nf <= 0) || (nt <= 0) ||
	   (input->numSigmaMin < 1.0) ||
	   (input->alphaDefault < 0.0) ||
	   (input->alphaDefault > 1.0))
		XLAL_ERROR(func, XLAL_EDOM);

	/* make sure TF planes have already been computed */
	if(!tfTiling->planesComputed)
		XLAL_ERROR(func, XLAL_EDATA);

	for(tile = tfTiling->firstTile; tile; tile = tile->nextTile) {
		t1 = tile->tstart;
		t2 = tile->tend;
		k1 = tile->fstart;
		k2 = tile->fend;

		if((t1 < 0) || (t1 > t2) || (t2 > nt) ||
		   (k1 < 0) || (k1 > k2) || (k2 > nf))
			XLAL_ERROR(func, XLAL_EDATA);

		/* Calculate the degrees of freedom of the TF tile */
		dof = 2.0 * (t2 - t1) * tile->deltaT * (k2 - k1) * tile->deltaF;

		sum = 0.0;
		for(offset = t1; offset < t2; offset += (t2 - t1) / dof)
			sum += square_complex8_sum(&plane->data[offset * nf], k1, k2 - k1) / XLALREAL4SumSquares(norm, k1, k2 - k1);

		tile->excessPower = sum - dof;
		tile->weight = 1.0;

		/* Need to compute an accurate value of likelihood only if
		 * excess power is greater than a few sigma */
		numsigma = (sum - dof) / sqrt(2.0 * dof);
		if(numsigma > input->numSigmaMin) {
			tile->firstCutFlag = TRUE;
			tile->alpha = XLALOneMinusChisqCdf(sum, dof);
			if(XLALIsREAL8FailNaN(tile->alpha))
				XLAL_ERROR(func, XLAL_EFUNC);
		} else {
			tile->firstCutFlag = FALSE;
			tile->alpha =  input->alphaDefault; /* default value */
		}
	}

	/* set flag saying alpha for each tile has been computed */
	tfTiling->excessPowerComputed = TRUE;

	/* success */
	return(0);
}
