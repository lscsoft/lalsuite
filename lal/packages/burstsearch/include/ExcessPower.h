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

#ifdef  __cplusplus		/* C++ protection. */
extern "C" {
#endif

NRCSID(EXCESSPOWERH, "$Id$");

typedef struct tagTFTile {
	INT4 fstart;
	INT4 fend;
	INT4 tstart;
	INT4 tend;
	REAL8 deltaT;
	REAL8 deltaF;
	REAL8 excessPower;
	REAL8 lnalpha;
	REAL8 lnweight;
} TFTile;


typedef struct tagTFTiling {
	TFTile *tile;
	size_t numtiles;
} TFTiling;


typedef struct tagCreateTFTilingIn {
	INT4 overlapFactor;
	INT4 minTimeBins;
	REAL8 maxTileBand;
	REAL8 maxTileDuration;
} CreateTFTilingIn;


int
XLALAddWhiteNoise(
	COMPLEX8Vector *v,
	REAL8 noiseLevel
);


void
LALAddWhiteNoise(
	LALStatus *status,
	COMPLEX8Vector *v,
	REAL8 noiseLevel
);


REAL8
XLALTFTileDegreesOfFreedom(
	TFTile *tile
);


TFTiling *
XLALCreateTFTiling(
	const CreateTFTilingIn *input, 
	const TFPlaneParams *planeparams
);


void
XLALDestroyTFTiling(
	TFTiling *tiling
);


int
XLALComputeExcessPower(
	TFTiling *tiling,
	const COMPLEX8TimeFrequencyPlane *plane,
	const REAL4 *norm
);


void
XLALSortTFTilingByAlpha(
	TFTiling *tiling
);


REAL8
XLALComputeLikelihood(
	TFTiling *tiling
);


void
XLALPrintTFTileList(
	FILE *fp,
	const TFTiling *tiling,
	const COMPLEX8TimeFrequencyPlane *plane,
	size_t maxTiles
);

#ifdef  __cplusplus
}
#endif				/* C++ protection. */
#endif
