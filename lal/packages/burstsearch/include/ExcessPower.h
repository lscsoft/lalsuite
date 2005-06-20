/********************************** <lalVerbatim file="ExcessPowerHV">
Author: Flanagan, E
$Id$
**************************************************** </lalVerbatim> */

#ifndef _EXCESSPOWER_H
#define _EXCESSPOWER_H

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
	REAL8 alpha;
	REAL8 weight;
	BOOLEAN firstCutFlag;
	struct tagTFTile *nextTile;
} TFTile;


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


TFTile *
XLALCreateTFTiling(
	const CreateTFTilingIn *input, 
	const TFPlaneParams *planeparams
);


void
XLALDestroyTFTiling(
	TFTile *list
);


int
XLALComputeExcessPower(
	TFTile *list,
	const COMPLEX8TimeFrequencyPlane *plane,
	REAL8 numSigmaMin,
	REAL8 alphaDefault,
	const REAL4 *norm
);


int
XLALSortTFTiling(
	TFTile **list
);


REAL8
XLALComputeLikelihood(
	TFTile *list
);


void
XLALPrintTFTileList(
	FILE *fp,
	const TFTile *list,
	const COMPLEX8TimeFrequencyPlane *plane,
	INT4 maxTiles
);

#ifdef  __cplusplus
}
#endif				/* C++ protection. */
#endif
