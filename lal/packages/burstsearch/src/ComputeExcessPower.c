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
		double sumsquares = 0.0;
		double hsumsquares = 0.0;
		unsigned t;

		for(t = tile->tstart; t < tile->tstart + tile->tbins; t += tstep) {
			double sum = 0.0;
			double hsum = 0.0;
			unsigned channel;

			for(channel = tile->channel0; channel < tile->channel0 + tile->channels; channel++) {
				sum += plane->channel[channel]->data[t];
				hsum += plane->channel[channel]->data[t] * sqrt(plane->channel_mean_square->data[channel]);
			}

			sumsquares += sum * sum / (tile->channels + (tile->channels - 1) * plane->channel_overlap);
			hsumsquares += hsum * hsum / (tile->channels + (tile->channels - 1) * plane->channel_overlap);
		}
		tile->excessPower = sumsquares - dof;
		tile->hrss = sqrt(hsumsquares);
		tile->lnalpha = XLALlnOneMinusChisqCdf(sumsquares, dof);
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
