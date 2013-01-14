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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALInspiral.h>
#include <lal/RealFFT.h>
#include <lal/LALNoiseModels.h>

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

typedef struct tagDiscoverInspiralEventsIn
{
  UINT4            currentGPSTime;
  INT4             nBegin, nEnd;
  INT4             chisqBins;
  REAL8            Threshold, ClusterThreshold;
  REAL8            dynRangeScalingFac;
  REAL4Vector      signal;
  REAL8Vector      psd;
  InspiralTemplate param;
  RealFFTPlan      *fwdp;
  RealFFTPlan      *revp;
  UINT2            displayCorrelationStats;
  UINT2            displayCorrelation;
}
DiscoverInspiralEventsIn;

typedef struct tagDiscoverInspiralEventsList
{
  UINT4               bin;
  INT4                endTime, endTimeNS, impulseTime, impulseTimeNS, chisqDOF;
  REAL8               amplitude, effDistance, effD8, phase, snr;
  REAL8               cmax1, cmax2, sigmasq;
  REAL8               chisq, chisq1, chisq2;
  REAL8               t0, t3, m1, m2;
}
DiscoverInspiralEventsList;

typedef struct tagInspiralChisqInput
{

  INT4                          chisqBins;
  REAL4Vector                   filter1, filter2;
  REAL4Vector                   rho1, rho2;
  REAL8                         flso;
  DiscoverInspiralEventsIn      *findEventsIn;
}
InspiralChisqInput;

typedef struct tagInspiralChisqOutput
{
  REAL8         *chisqZERO;
  REAL8         *chisqPIbyTWO;
  REAL8         *chisq;
}
InspiralChisqOutput;

typedef struct tagInvSpecTruncationParams
{
  INT4              n;
  REAL8             df;
  RealFFTPlan       *fwdp;
  RealFFTPlan       *revp;
  REAL8             psdTruncTime;
  UINT2             ifDebug;
}
InvSpecTruncationParams;

typedef struct
tagRandomInspiralSignalIn
{
  INT4 useed;       /**< Seed for the random number generator */
  INT4 type;        /**< Type of signal required to be generated */

  REAL8 mMin;       /**< smallest component mass allowed */
  REAL8 mMax;       /**< largest component mass allowed */
  /* OR */
  REAL8 MMax;       /**< largest total mass allowed */
  REAL8 MMin;       /**< largest total mass allowed */
  REAL8 SignalAmp;  /**< amplitude of the signal (relevant only when type=2) */
  REAL8 NoiseAmp;   /**< amplitude of noise (relevant only when type=2) */
  REAL8 etaMin;     /**< smallest value of the symmetric mass ratio */
  InspiralTemplate param;/**< parameter stuct; user to specify certain params. */
  REAL8Vector psd;  /**< power spectral density used for coloring the noise */
  RealFFTPlan *fwdp;/**< pre-computed fftw plan for forward fftw */

  /* Chirp times are needed only if param.massChoice is t02 or t03 */
  REAL8 t0Min;      /**< smallest Newtonian chirp time */
  REAL8 t0Max;      /**< largest Newtonian chirp time */
  REAL8 tnMin;      /**< smallest 1, 1.5 PN chirp time */
  REAL8 tnMax;      /**< largest 1, 1.5 PN chirp time  */
  /* min/max values of BCV parameters*/
  REAL8 psi0Min;      /**< smallest Newtonian psi-parameter */
  REAL8 psi0Max;      /**< largest Newtonian psi-parameter */
  REAL8 psi3Min;      /**< smallest 1.5 PN psi-parameter */
  REAL8 psi3Max;      /**< largest 1.5 PN psi-parameter */

  INT4  coalescenceTime ;/**< bin in which is maximum of the chirp (coalescence time)*/

  /** \name These are for spin Taylor waveforms */
  /*@{*/
  REAL8 minDistance, maxDistance;
  REAL8 spin1min, spin1max, spin2min, spin2max;
  REAL8 theta0min, theta0max, phi0min, phi0max;
  REAL8 polarisationAngleMin, polarisationAngleMax;
  REAL8 sourceThetaMin, sourceThetaMax, sourcePhiMin, sourcePhiMax;
  REAL8 inclinationMin, inclinationMax;
  /*@}*/
}
RandomInspiralSignalIn;

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
tagInspiralFindEventsIn
{
  UINT4            currentGPSTime;
  INT4             nBegin;
  INT4             nEnd;
  REAL8            Threshold;
  REAL8            ClusterThreshold;
  REAL4Vector      signal;
  REAL8Vector      psd;
  InspiralTemplate param;
  RealFFTPlan      *fwdp;
  RealFFTPlan      *revp;
  UINT2            displayData;
  UINT2            displayPSD;
  UINT2            displayTemplates;
  UINT2            displayCorrelation;
  UINT2            displayCorrelationStats;
}
InspiralFindEventsIn;

typedef struct
tagInspiralEventsList
{
  UINT4            bin;

  INT4             endTime;
  INT4             endTimeNS;
  INT4             impulseTime;
  INT4             impulseTimeNS;
  INT4             chisqDOF;

  REAL8            amplitude;
  REAL8            effDistance;
  REAL8            phase;
  REAL8            snr;
  REAL8            sigmasq;
  REAL8            chisq;

  InspiralTemplate param;
}
InspiralEventsList;

typedef struct
tagInspiralWaveNormaliseIn
{
  REAL8        df;
  REAL8        fCutoff;
  REAL8        samplingRate;
  REAL8Vector *psd;
}
InspiralWaveNormaliseIn;

typedef struct
tagInspiralChisqDataVec
{
  REAL4Vector *SNRIntegrand;
  REAL8Vector *psd;
}
InspiralChisqDataVec;

typedef struct
tagInspiralChisqParams
{
  INT4 nBins;   /* number of chi-squared bins to use */

  REAL8 totalMass;
  REAL8 fLower;
  REAL8 deltaT; /* sampling interval */
}
InspiralChisqParams;

typedef struct
tagInspiralSNRIntegrandParams
{
  INT4  lag;    /* the value of the lag which produced the largest correlation */

  REAL8 phase;  /* phase of the correlation where the max occurs */
  REAL8 deltaT; /* sampling interval */
}
InspiralSNRIntegrandParams;


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
LALRandomInspiralSignal
(
 LALStatus *status,
 REAL4Vector *signalvec,
 RandomInspiralSignalIn *randIn
 );

void
LALInspiralWaveOverlap
(
 LALStatus               *status,
 REAL4Vector             *output,
 InspiralWaveOverlapOut  *overlapout,
 InspiralWaveOverlapIn   *overlapin
 );

void LALInspiralGetOrthoNormalFilter(REAL4Vector *filter2, REAL4Vector *filter1);

void
LALInspiralFindEvents
(
 LALStatus   *status,
 INT4  *nEvents,
 InspiralEventsList   **findeventsout,
 InspiralFindEventsIn *findeventsin
 );

void
LALInspiralFindLoudestEvent
(
 LALStatus            *status,
 INT4                 *nEvents,
 InspiralEventsList   *eventlist,
 InspiralFindEventsIn *findeventsin
 );

void
LALInspiralFindEventsCluster
(
 LALStatus            *status,
 INT4                 *nEvents,
 InspiralEventsList   **findeventsout,
 InspiralFindEventsIn *findeventsin
 );

void
LALInspiralComputeChisq
(
 LALStatus *status,
 REAL4 *chisq,
 InspiralChisqDataVec *input,
 InspiralChisqParams *params
 );

void
LALInspiralComputeSNRIntegrand
(
 LALStatus *status,
 REAL4Vector *output,
 InspiralWaveCorrelateIn corrin,
 InspiralSNRIntegrandParams *params
 );

void LALDiscoverInspiralEvents
(
 LALStatus                     *status,
 INT4                          *nEvents,
 DiscoverInspiralEventsList    **eventlist,
 DiscoverInspiralEventsIn      *findeventsin
 );

void LALEstimateEffectiveDistance
(
 LALStatus          *status,
 InspiralTemplate    param,
 REAL8               df,
 REAL8Vector        *psd,
 REAL8               lal_nm_snr,
 REAL8              *effDistance
 );

void LALEvaluateInspiralChisqTest
(
 LALStatus             *status,
 InspiralChisqOutput   *chisqOut,
 InspiralChisqInput    *chisqIn
 );

void LALTruncateInvSpectrum
(
 LALStatus               *status,
 REAL8Vector             *inputVec,
 InvSpecTruncationParams *params
 );

void GenerateTimeDomainWaveformForInjection (
    LALStatus              *status,
    REAL4Vector            *buff,
    InspiralTemplate       *params
    );

#ifdef  __cplusplus
}
#endif

#endif /* _LALNOISEMODELSINSPIRAL_H */
