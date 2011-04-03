/*
*  Copyright (C) 2007 Jolien Creighton
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#ifndef _BLACKHOLEMODE_H
#define _BLACKHOLEMODE_H

#include <lal/LALDatatypes.h>

NRCSID( BLACKHOLEMODEH, "$Id$" );


#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
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


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _BLACKHOLEMODE_H */
