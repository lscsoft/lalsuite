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
	/* tile specification as indeces into time-frequency plane data */
	INT4 fstart;
	INT4 fbins;
	INT4 tstart;
	INT4 tbins;
	/* time-frequency plan parameters */
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
	const REAL4TimeFrequencyPlane *plane,
	const REAL8 *hrssfactor,
	const REAL4 *norm
);


REAL8
XLALComputeLikelihood(
	TFTiling *tiling
);


void
XLALPrintTFTileList(
	FILE *fp,
	const TFTiling *tiling,
	const REAL4TimeFrequencyPlane *plane,
	size_t maxTiles
);


int
XLALAddWhiteNoise(
	COMPLEX8Sequence *v,
	REAL8 noiseLevel
);


#ifdef  __cplusplus
}
#endif				/* C++ protection. */
#endif
