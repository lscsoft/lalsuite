/* <lalVerbatim file="LALNoiseModelsHV">

Author: Sathyaprakash, B.S.
$Id$

</lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALNoiseModels.h}}
\label{s:LALNoiseModels.h}

Header file for model noise generation codes.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALNoiseModels.h>
\end{verbatim}

\noindent This header file covers routines that are used in 
synthetic background noise  expected in various
detectors and signals with random parameters in background noise.

</lalLaTeX> */

#ifndef _LALNOISEMODELS_H
#define _LALNOISEMODELS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALInspiral.h>
#include <lal/RealFFT.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( LALNOISEMODELSH, "$Id$" );

/* <lalLaTeX>

\subsection*{Error codes}

</lalLaTeX>  */

/* <lalErrTable> */

#define LALNOISEMODELSH_ENULL 1
#define LALNOISEMODELSH_EMEM 2
#define LALNOISEMODELSH_ECHOICE 4
#define LALNOISEMODELSH_EDIV0 8
#define LALNOISEMODELSH_ESIZE 16
#define LALNOISEMODELSH_MSGENULL "Arguments contained an unexpected null pointer"
#define LALNOISEMODELSH_MSGEMEM "Memory allocation error"
#define LALNOISEMODELSH_MSGECHOICE "Invalid choice for an input parameter"
#define LALNOISEMODELSH_MSGEDIV0 "Division by zero"
#define LALNOISEMODELSH_MSGESIZE "Invalid input size"

/* </lalErrTable> */

/* <lalLaTeX>

\section*{Structures}
\input{LALNoiseModelsHS}
</lalLaTeX>  */

/*  <lalVerbatim file="LALNoiseModelsHS"> */

typedef enum
{
  geo, 
  ligo, 
  tama, 
  virgo
}
Detector;

/*  </lalVerbatim>  */

/*  <lalLaTeX> 
\idx[Type]{Detector} 
</lalLaTeX>  */

/*  <lalVerbatim file="LALNoiseModelsHS"> */
typedef struct
tagInspiralWaveCorrelateIn 
{
   REAL4Vector signal1;
   REAL4Vector signal2;
   REAL8Vector psd;
   RealFFTPlan *revp;
}
InspiralWaveCorrelateIn;
/* </lalVerbatim>  */

/*  <lalLaTeX>
\idx[Type]{InspiralWaveCorrelateIn}
</lalLaTeX>  */

/*  <lalVerbatim file="LALNoiseModelsHS"> */
typedef struct 
tagAddVectorsIn
{
   REAL4Vector *v1;
   REAL4Vector *v2;
   REAL8       a1;
   REAL8       a2;
} 
AddVectorsIn;
/*  </lalVerbatim>  */

/*  <lalLaTeX> 
\idx[Type]{AddVectorsIn} 
</lalLaTeX>  */

/*  <lalVerbatim file="LALNoiseModelsHS"> */
typedef struct 
tagRandomInspiralSignalIn
{
   INT4 useed;       /* Seed for the random number generator */
   INT4 type;        /* Type of signal required to be generated */

   REAL8 mMin;       /* smallest component mass allowed */
   REAL8 MMax;       /* largest total mass allowed */
   REAL8 SignalAmp;  /* amplitude of the signal (relevant only when type=2) */
   REAL8 NoiseAmp;   /* amplitude of noise (relevant only when type=2) */
   REAL8 etaMin;     /* smallest value of the symmetric mass ratio */
   InspiralTemplate 
	 param;      /* parameter stuct; user to specify certain params. */

   REAL8Vector psd;  /* power spectral density used for coloring the noise */
   RealFFTPlan *fwdp;/* pre-computed fftw plan for forward fftw */

   /* Chirp times are needed only if param.massChoice is t02 or t03 */
   REAL8 t0Min;      /* smallest Newtonian chirp time */
   REAL8 t0Max;      /* largest Newtonian chirp time */
   REAL8 tnMin;      /* smallest 1, 1.5 PN chirp time */
   REAL8 tnMax;      /* largest 1, 1.5 PN chirp time  */

} 
RandomInspiralSignalIn;
/*  </lalVerbatim>  */

/*  <lalLaTeX> 
\idx[Type]{RandomInspiralSignalIn} 
</lalLaTeX>  */

/*  <lalVerbatim file="LALNoiseModelsHS"> */
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
} 
InspiralWaveOverlapIn; 
/*  </lalVerbatim>  */

/*  <lalLaTeX> 
\idx[Type]{InspiralwaveOverlapIn} 
</lalLaTeX>  */

/*  <lalVerbatim file="LALNoiseModelsHS"> */
typedef struct
tagInspiralWaveOverlapOut 
{
   REAL8 max, phase;
   INT4 bin;
} 
InspiralWaveOverlapOut; 
/*  </lalVerbatim>  */

/*  <lalLaTeX> 
\idx[Type]{InspiralWaveOverlapOut} 
</lalLaTeX>  */


/*  <lalVerbatim file="LALNoiseModelsHS"> */
typedef struct
tagInspiralFindEventsIn 
{
   INT4             nBegin;
   INT4             nEnd;
   REAL8            Threshold;
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
/*  </lalVerbatim>  */

/*  <lalLaTeX> 
\index{\texttt{InspiralFindEventsIn}} 
</lalLaTeX>  */

/*  <lalVerbatim file="LALNoiseModelsHS"> */
typedef struct
tagInspiralEventsList 
{
   REAL8            max;
   REAL8            phase;
   REAL8            eventTime;
   REAL8            eventTime_ns;

   UINT4            bin;

   InspiralTemplate param;
} 
InspiralEventsList; 
/*  </lalVerbatim>  */

/*  <lalLaTeX> 
\index{\texttt{InspiralEventsList}} 
</lalLaTeX>  */


/*  <lalVerbatim file="LALNoiseModelsHS"> */
typedef struct
tagStatsREAL4VectorOut
{
  REAL8 mean;
  REAL8 var;
  REAL8 stddev;
  REAL8 min;
  REAL8 max;
} 
StatsREAL4VectorOut;
/*  </lalVerbatim>  */

/*  <lalLaTeX> 
\index{\texttt{StatsREAL4VectorOut}} 
</lalLaTeX>  */



/*  <lalLaTeX>
\vfill{\footnotesize\input{LALNoiseModelsHV}}
</lalLaTeX>  */

/* Function prototypes */

/* <lalLaTeX>
\newpage\input{LALNoiseSpectralDensityC}
</lalLaTeX>  */

void 
LALNoiseSpectralDensity 
   (
   LALStatus   *status, 
   REAL8Vector *psd, 
   void        (*NoisePsd)(LALStatus *status, REAL8 *shf, REAL8 f),
   REAL8       f
   );

/* <lalLaTeX>
\newpage\input{LALInspiralWaveCorrelateC}
</lalLaTeX>  */

void 
LALInspiralWaveCorrelate 
   (
   LALStatus   *status, 
   REAL4Vector *output, 
   InspiralWaveCorrelateIn in
   );

/* <lalLaTeX>
\newpage\input{LALInspiralWaveNormaliseC}
</lalLaTeX>  */

void 
LALInspiralWaveNormalise 
   (
   LALStatus   *status, 
   REAL4Vector *dh, 
   REAL8       *norm, 
   REAL8Vector psd
   );

/* <lalLaTeX>
\newpage\input{LALGEOPsdC}
</lalLaTeX>  */

void 
LALGEOPsd 
   (
   LALStatus *status, 
   REAL8     *shf, 
   REAL8     x
   );

/* <lalLaTeX>
\newpage\input{LALLIGOIPsdC}
</lalLaTeX>  */

void 
LALLIGOIPsd 
   (
   LALStatus *status, 
   REAL8     *shf, 
   REAL8     x
   );

/* <lalLaTeX>
\newpage\input{LALTAMAPsdC}
</lalLaTeX>  */

void 
LALTAMAPsd 
   (
   LALStatus *status, 
   REAL8     *shf, 
   REAL8     x
   );

/* <lalLaTeX>
\newpage\input{LALVIRGOPsdC}
</lalLaTeX>  */

void 
LALVIRGOPsd 
   (
   LALStatus *status, 
   REAL8     *shf, 
   REAL8     x
   );


/* <lalLaTeX>
\newpage\input{LALRandomInspiralSignalC}
</lalLaTeX>  */

void 
LALRandomInspiralSignal
   (
   LALStatus *status, 
   REAL4Vector *signal,
   RandomInspiralSignalIn *randIn
   );

/* <lalLaTeX>
\newpage\input{LALColoredNoiseC}
</lalLaTeX>  */

void 
LALColoredNoise 
   (
   LALStatus   *status,
   REAL4Vector *noisy, 
   REAL8Vector psd
   );

/* <lalLaTeX>
\newpage\input{LALAddVectorsC}
</lalLaTeX>  */

void 
LALAddVectors
   (
   LALStatus *status, 
   REAL4Vector *vector, 
   AddVectorsIn in);

/* <lalLaTeX>
\newpage\input{LALInspiralWaveOverlapC}
</lalLaTeX>  */

void 
LALInspiralWaveOverlap 
   (
   LALStatus   *status,
   REAL4Vector *output,
   InspiralWaveOverlapOut  *overlapout,
   InspiralWaveOverlapIn   *overlapin
   );

/* <lalLaTeX>
\newpage\input{LALInspiralFindEventsC}
</lalLaTeX>  */

void 
LALInspiralFindEvents 
   (
   LALStatus   *status,
   INT4  *nEvents,
   InspiralEventsList   **findeventsout,
   InspiralFindEventsIn *findeventsin
   );

/* <lalLaTeX>
\newpage\input{LALStatsREAL4VectorC}
</lalLaTeX>  */

void 
LALStatsREAL4Vector
   (
   LALStatus *status, 
   StatsREAL4VectorOut *out, 
   REAL4Vector *vector
   );

/* <lalLaTeX>
\newpage\input{FilterTestC}
</lalLaTeX> */

/* <lalLaTeX>
\newpage\input{RandomInspiralSignalTestC}
</lalLaTeX> */

/* <lalLaTeX>
\newpage\input{NoisePSDTestC}
</lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif /* _LALNOISEMODELS_H */
