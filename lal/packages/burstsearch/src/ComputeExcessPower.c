/******** <lalVerbatim file="ComputeExcessPowerCV"> ********
Author: Flanagan, E. and Cannon, K.
$Id$
********* </lalVerbatim> ********/

#include <math.h>
#include <lal/LALRCSID.h>

NRCSID (COMPUTEEXCESSPOWERC, "$Id$");

#include <lal/ExcessPower.h>
#include <lal/TFTransform.h>
#include <lal/Thresholds.h>
#include <lal/XLALError.h>


/******** <lalVerbatim file="ComputeExcessPowerCP"> ********/
int XLALComputeExcessPower(
	TFTiling *tiling,
	const REAL4TimeFrequencyPlane *plane
)
/******** </lalVerbatim> ********/
{
	const char func[] = "XLALComputeExcessPower";
	TFTile *tile = tiling->tile;
	size_t i;

	for(i = 0; i < tiling->numtiles; i++, tile++) {
		const double dof = XLALTFTileDegreesOfFreedom(tile);
		const unsigned tstep = tile->tbins / dof;
		double sum = 0.0;
		unsigned t, channel;

		for(t = tile->tstart; t < tile->tstart + tile->tbins; t += tstep) {
			double tmp = 0.0;
			for(channel = tile->channel0; channel < tile->channel0 + tile->channels; channel++)
				tmp += plane->channel[channel]->data[t];
			sum += pow(tmp, 2.0) / (tile->channels + (tile->channels - 1) * plane->channel_overlap);
		}
		tile->excessPower = sum - dof;
		tile->hrss = 0;
		tile->lnalpha = XLALlnOneMinusChisqCdf(sum, dof);
		if(XLALIsREAL8FailNaN(tile->lnalpha))
			XLAL_ERROR(func, XLAL_EFUNC);
	}

	/* success */
	return 0;
}


/******** <lalVerbatim file="ComputeLikelihoodCP"> ********/
REAL8 XLALComputeLikelihood(
	TFTiling *tiling
)
/******** </lalVerbatim> ********/
{
	REAL8 avglambda = 0.0;
	REAL8 rho4;
	TFTile *tile;
	size_t i;

	for(i = 0, tile = tiling->tile; i < tiling->numtiles; i++, tile++) {
		rho4 = tile->excessPower * tile->excessPower;
		avglambda += XLALTFTileDegreesOfFreedom(tile) / rho4 * exp(tile->lnweight - tile->lnalpha);
	}

	/* compute the likelihood averaged over TF tiles */
	avglambda /= tiling->numtiles;

	/* return value of statistic */
	return avglambda;
}
