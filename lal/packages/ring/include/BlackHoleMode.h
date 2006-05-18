#ifndef _BLACKHOLEMODE_H
#define _BLACKHOLEMODE_H

#include <lal/LALDatatypes.h>

NRCSID( BLACKHOLEMODEH, "$Id$" );


#ifdef  __cplusplus
extern "C" {
#pragma } /* to match the previous brace */
#endif


typedef struct tagBlackHoleMode {
  REAL8 a;
  int l;
  int m;
  int s;
  COMPLEX16 A;
  COMPLEX16 omega;
} BlackHoleMode;



typedef struct tagBlackHoleModeEigenvalues
{
  REAL8 a;
  COMPLEX16 A;
  COMPLEX16 omega;
}
BlackHoleModeEigenvalues;

typedef struct tagBlackHoleModeTable
{
  int l;
  int m; /* absolute value */
  int s;
  /* these describe the table of mode values */
  size_t eigenTableSize;
  struct tagBlackHoleModeEigenvalues *positiveEigenTable;
  struct tagBlackHoleModeEigenvalues *negativeEigenTable;
}
BlackHoleModeTable;


/* low-level routines: these use Leaver's conventions */
void XLALDestroyBlackHoleModeTable( BlackHoleModeTable *mode );
BlackHoleModeTable * XLALCreateBlackHoleModeTable( int l, int m, int s );
int XLALBlackHoleModeEigenvaluesLeaverT( COMPLEX16 *A, COMPLEX16 *omega, REAL8 a, BlackHoleModeTable *mode );
int XLALBlackHoleModeEigenvaluesLeaver( COMPLEX16 *A, COMPLEX16 *omega, REAL8 a, int l, int m, int s );


/* medium-level routines: the input values are in standard conventions but the
 * parameter structure has data in Leaver's conventions */
int XLALSetBlackHoleModeParams( struct tagBlackHoleMode *params, REAL8 a, int l, int m, int s );
int XLALSpheroidalWaveFunction1( COMPLEX16 *result, REAL8 mu, struct tagBlackHoleMode *mode );
int XLALSpheroidalWaveFunctionNorm( COMPLEX16 *norm, struct tagBlackHoleMode *params );
int XLALSpheroidalWaveFunction( COMPLEX16 *result, REAL8 mu, struct tagBlackHoleMode *params );


/* high-level routines: these use standard conventions */
int XLALSpinWeightedSpheroidalHarmonic( COMPLEX16 *result, REAL8 mu, REAL8 a, int l, int m, int s );

int XLALBlackHoleRingdownAmplitude(
    COMPLEX16 *amplitudePlus,
    COMPLEX16 *amplitudeCross,
    REAL8 massSolar,
    REAL8 fracMassLoss,
    REAL8 distanceMpc,
    REAL8 inclinationRad,
    REAL8 azimuthRad,
    BlackHoleMode *mode
    );
int XLALBlackHoleRingdownWaveform(
    REAL4TimeSeries *plus,
    REAL4TimeSeries *cross,
    REAL8 massSolar,
    REAL8 spinParam,
    REAL8 fracMassLoss,
    REAL8 distanceMpc,
    REAL8 inclinationRad,
    REAL8 azimuthRad,
    int l,
    int m
    );


#ifdef  __cplusplus
#pragma { /* to match the next brace */
}
#endif

#endif /* _BLACKHOLEMODE_H */
