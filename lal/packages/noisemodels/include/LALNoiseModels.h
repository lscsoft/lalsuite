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
#define LALNOISEMODELSH_ECHOICE 3
#define LALNOISEMODELSH_EDIV0 4
#define LALNOISEMODELSH_ESIZE 8
#define LALNOISEMODELSH_MSGENULL "Arguments contained an unexpected null pointer"
#define LALNOISEMODELSH_MSGEMEM "Memory allocation error"
#define LALNOISEMODELSH_MSGECHOICE "Invalid choice for an input parameter"
#define LALNOISEMODELSH_MSGEDIV0 "Division by zero"
#define LALNOISEMODELSH_MSGESIZE "Invalid input range"

/* </lalErrTable> */

/* <lalLaTeX>

\section*{Structures}
\input{LALNoiseModelsHS}
</lalLaTeX>  */

/*  <lalVerbatim file="LALNoiseModelsHS"> */

typedef enum
{
  geo, ligo, tama, virgo
}
Detector;

/*  </lalVerbatim>  */

/*  <lalLaTeX> 
\index{\texttt{Detector}} 
</lalLaTeX>  */

/*  <lalVerbatim file="LALNoiseModelsHS"> */
typedef struct
tagCorrelateIn {
   REAL4Vector signal1, signal2;
   REAL8Vector psd;
   RealFFTPlan *revp;
}
CorrelateIn;
/* </lalVerbatim>  */

/*  <lalLaTeX>
\index{\texttt{CorrelateIn}}
</lalLaTeX>  */

/*  <lalVerbatim file="LALNoiseModelsHS"> */
typedef struct 
tagAddVectorsIn
{
   REAL4Vector *v1, *v2;
   REAL8 a1, a2;
} 
AddVectorsIn;
/*  </lalVerbatim>  */

/*  <lalLaTeX> 
\index{\texttt{AddVectorsIn}} 
</lalLaTeX>  */

/*  <lalVerbatim file="LALNoiseModelsHS"> */
typedef struct 
tagRandomInspiralSignalIn
{
   InspiralTemplate param;
   UINT8 useed;
   REAL8 mMin, MMax, SignalAmp, NoiseAmp;
   INT4 type;
   REAL8Vector psd;
   RealFFTPlan *fwdp;
} 
RandomInspiralSignalIn;
/*  </lalVerbatim>  */

/*  <lalLaTeX> 
\index{\texttt{RandomInspiralSignalIn}} 
</lalLaTeX>  */

/*  <lalVerbatim file="LALNoiseModelsHS"> */
typedef struct
tagOverlapIn {
   REAL4Vector signal;
   REAL8Vector psd;
   InspiralTemplate param;
   RealFFTPlan *fwdp, *revp;
} 
OverlapIn; 
/*  </lalVerbatim>  */

/*  <lalLaTeX> 
\index{\texttt{OverlapIn}} 
</lalLaTeX>  */

/*  <lalVerbatim file="LALNoiseModelsHS"> */
typedef struct
tagOverlapOut {
   REAL8 max, phase;
   INT4 bin;
} 
OverlapOut; 
/*  </lalVerbatim>  */

/*  <lalLaTeX> 
\index{\texttt{OverlapOut}} 
</lalLaTeX>  */

/*  <lalLaTeX>
\vfill{\footnotesize\input{LALNoiseModelsHV}}
</lalLaTeX>  */

/* Function prototypes */

/* <lalLaTeX>
\newpage\input{LALNoiseSpectralDensityC}
</lalLaTeX>  */

void LALNoiseSpectralDensity (LALStatus   *status, 
                              REAL8Vector *psd, 
                              void (*NoisePsd)(LALStatus *status, REAL8 *shf, REAL8 f),
                              REAL8       f);

/* <lalLaTeX>
\newpage\input{LALCorrelateC}
</lalLaTeX>  */

void LALCorrelate (LALStatus   *status, 
                   REAL4Vector *output, 
                   CorrelateIn in);

/* <lalLaTeX>
\newpage\input{LALNormaliseC}
</lalLaTeX>  */

void LALNormalise (LALStatus   *status, 
                   REAL4Vector *dh, 
                   REAL8       *norm, 
                   REAL8Vector psd);

/* <lalLaTeX>
\newpage\input{LALGEOPsdC}
</lalLaTeX>  */

void LALGEOPsd (LALStatus *status, REAL8 *shf, REAL8 x);

/* <lalLaTeX>
\newpage\input{LALLIGOIPsdC}
</lalLaTeX>  */

void LALLIGOIPsd (LALStatus *status, REAL8 *shf, REAL8 x);

/* <lalLaTeX>
\newpage\input{LALTAMAPsdC}
</lalLaTeX>  */

void LALTAMAPsd (LALStatus *status, REAL8 *shf, REAL8 x);

/* <lalLaTeX>
\newpage\input{LALVIRGOPsdC}
</lalLaTeX>  */

void LALVIRGOPsd (LALStatus *status, REAL8 *shf, REAL8 x);


/* <lalLaTeX>
\newpage\input{LALRandomInspiralSignalC}
</lalLaTeX>  */

void LALRandomInspiralSignal(LALStatus *status, 
                             REAL4Vector *signal,
                             RandomInspiralSignalIn *randIn);

/* <lalLaTeX>
\newpage\input{LALGaussianNoiseC}
</lalLaTeX>  */

void LALGaussianNoise (LALStatus   *status,
                       REAL4Vector *noisy, 
                       UINT8       *seed);

/* <lalLaTeX>
\newpage\input{LALColoredNoiseC}
</lalLaTeX>  */

void LALColoredNoise (LALStatus   *status,
                      REAL4Vector *noisy, 
                      REAL8Vector psd);

/* <lalLaTeX>
\newpage\input{LALGaussianRandomNumberC}
</lalLaTeX>  */

void LALGaussianRandomNumber(LALStatus *status, 
                             REAL4 *randnum);

/* <lalLaTeX>
\newpage\input{LALAddVectorsC}
</lalLaTeX>  */

void LALAddVectors(LALStatus *status, 
                   REAL4Vector *vector, 
                   AddVectorsIn in);

/* <lalLaTeX>
\newpage\input{LALInspiralWaveOverlapC}
</lalLaTeX>  */

void LALInspiralWaveOverlap (LALStatus   *status,
                             REAL4Vector *output,
                             OverlapOut  *overlapout,
                             OverlapIn   *overlapin);

#ifdef  __cplusplus
}
#endif

#endif /* _LALNOISEMODELS_H */
