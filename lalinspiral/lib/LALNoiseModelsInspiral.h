/*
 *  Copyright (C) 2007 Stas Babak, David Churches, Duncan Brown, Jolien Creighton, B.S. Sathyaprakash, Anand Sengupta, Thomas Cokelaer
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

#ifndef _LALNOISEMODELSINSPIRAL_H
#define _LALNOISEMODELSINSPIRAL_H

#include <lal/LALDatatypes.h>
#include <lal/LALInspiral.h>
#include <lal/RealFFT.h>

#ifdef  __cplusplus
extern "C" {
#endif

typedef struct
tagInspiralWaveCorrelateIn
{
  REAL8        df;
  REAL8        fCutoff;
  REAL8        samplingRate;
  REAL4Vector  signal1;
  REAL4Vector  signal2;
  REAL8Vector  psd;
  RealFFTPlan *revp;
}
InspiralWaveCorrelateIn;

typedef struct
tagInspiralWaveOverlapIn
{
  INT4             nBegin;
  INT4             nEnd;
  REAL4Vector      signal;
  REAL8Vector      psd;
  InspiralTemplate param;
  RealFFTPlan      *fwdp;
  RealFFTPlan      *revp;
  UINT2            ifExtOutput;       /* A flag which takes values 0 or 1 to denote
                                         if an extended output consisting of filter
                                         and xcorr vectors need to be filled out in
                                         the call to LALInspiralWaveOverlap ( )
                                         */
  UINT2             ifCorrelationOutput;/* a flag to fill the xcorr1 and xcorr2 outputs*/
}
InspiralWaveOverlapIn;

typedef struct
tagInspiralWaveOverlapOut
{
  REAL8        max, phase, alpha;
  INT4         bin;                /* bin at which max occurs */
  REAL4Vector  *filter1, *filter2; /* zero and pi/2 phase templates */
  REAL4Vector  *xcorr1, *xcorr2;   /* cross correlation against filter 1/2 */
}
InspiralWaveOverlapOut;

typedef struct
tagInspiralWaveNormaliseIn
{
  REAL8        df;
  REAL8        fCutoff;
  REAL8        samplingRate;
  REAL8Vector *psd;
}
InspiralWaveNormaliseIn;


/* Function prototypes */

void
LALInspiralWaveCorrelate
(
 LALStatus   *status,
 REAL4Vector *output,
 InspiralWaveCorrelateIn in
 );

void
LALInspiralWaveNormalise
(
 LALStatus   *status,
 REAL4Vector *dh,
 REAL8       *norm,
 REAL8Vector psd
 );

void
LALInspiralWaveNormaliseLSO
(
 LALStatus               *status,
 REAL4Vector             *filter,
 REAL8                   *norm,
 InspiralWaveNormaliseIn *in
 );

void
LALInspiralWaveOverlap
(
 LALStatus               *status,
 REAL4Vector             *output,
 InspiralWaveOverlapOut  *overlapout,
 InspiralWaveOverlapIn   *overlapin
 );

#ifdef  __cplusplus
}
#endif

#endif /* _LALNOISEMODELSINSPIRAL_H */
