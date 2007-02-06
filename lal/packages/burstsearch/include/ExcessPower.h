/********************************** <lalVerbatim file="ExcessPowerHV">
Author: Flanagan, E
$Id$
**************************************************** </lalVerbatim> */

#ifndef _EXCESSPOWER_H
#define _EXCESSPOWER_H

#include <stdio.h>
#include <stdlib.h>
#include <lal/LALDatatypes.h>
#include <lal/TFTransform.h>
#include <lal/LALRCSID.h>
#include <lal/Random.h>

#ifdef  __cplusplus		/* C++ protection. */
extern "C" {
#endif

NRCSID(EXCESSPOWERH, "$Id$");

typedef struct tagTFTile {
	/* tile specification as indexes into time-frequency plane data */
	UINT4 channel0;
	UINT4 channels;
	UINT4 tstart;
	UINT4 tbins;
	/* time-frequency plane parameters for reconstructing tile dimensions */
	REAL8 flow;
	REAL8 deltaT;
	REAL8 deltaF;
	/* computed tile properties */
	REAL8 excessPower;
	REAL8 hrss;
	REAL8 lnalpha;
	REAL8 lnweight;
} TFTile;


typedef struct tagTFTiling {
	/* array of tiles */
	TFTile *tile;
	size_t numtiles;
} TFTiling;


REAL8
XLALTFTileDegreesOfFreedom(
	const TFTile *tile
);


TFTiling *
XLALCreateTFTiling(
	const REAL4TimeFrequencyPlane *plane,
	INT4 inv_fractional_stride,
	REAL8 maxTileBandwidth,
	REAL8 maxTileDuration
);


void
XLALDestroyTFTiling(
	TFTiling *tiling
);


int
XLALComputeExcessPower(
	TFTiling *tiling,
	const REAL4TimeFrequencyPlane *plane
);


REAL8
XLALComputeLikelihood(
	TFTiling *tiling
);


REAL4Sequence *XLALREAL4AddWhiteNoise(
	REAL4Sequence *sequence,
	REAL4 rms,
	RandomParams *params
);

COMPLEX8Sequence *XLALCOMPLEX8AddWhiteNoise(
	COMPLEX8Sequence *sequence,
	REAL8 rms,
	RandomParams *params
);


#ifdef  __cplusplus
}
#endif				/* C++ protection. */
#endif
