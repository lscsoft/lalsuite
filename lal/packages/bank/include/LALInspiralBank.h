/* <lalVerbatim file="LALInspiralBankHV">

Author: Churches, D. K and B.S. Sathyaprakash.
$Id$

</lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALInspiralBank.h}}
\label{s:LALInspiralBank.h}

Header file for the template placement codes.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALInspiralBank.h>
\end{verbatim}

\noindent This header file covers routines that are used in template placement.

</lalLaTeX> */

#ifndef _LALINSPIRALBANK_H
#define _LALINSPIRALBANK_H

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

NRCSID( LALINSPIRALBANKH, "$Id$" );

/* <lalLaTeX>

\subsection*{Error codes}

</lalLaTeX>  */

/* <lalErrTable> */

#define LALINSPIRALBANKH_ENULL 1
#define LALINSPIRALBANKH_EDIV0 2
#define LALINSPIRALBANKH_MSGENULL "Arguments contained an unexpected null pointer"
#define LALINSPIRALBANKH_MSGEDIV0 "Division by zero"

/* </lalErrTable> */

/* <lalLaTeX>

\section*{Structures}
\input{LALInspiralBankHS}
</lalLaTeX>  */

/*  <lalVerbatim file="LALInspiralBankHS"> */

typedef enum
{
  geo, ligo, tama, virgo
}
Detector;

typedef enum
{
  Tau0Tau2, Tau0Tau3
}
CoordinateSpace;


/*  </lalVerbatim>  */

/*  <lalLaTeX>
\index{\texttt{Detector}}
</lalLaTeX>  */

/*  <lalLaTeX>
\index{\texttt{LALInspiralMetric}}
</lalLaTeX>  */

/* metric and its dimension */

/*  <lalVerbatim file="LALInspiralBankHS"> */
typedef struct 
tagInspiralMetric 
{
   REAL8            g00;     /* 00-component of the diagonalised metric. */
   REAL8            g11;     /* 11-component of the diagonalised metric. */
   REAL8            theta;   /* Angle from t0 to x0 */
   Detector         detector;/* detector for noise psd */ 
   CoordinateSpace space;    /* Coordinate space in which metric is computed */
   INT4 iflso;
} 
InspiralMetric;
/* </lalVerbatim>  */


/* a grid of inspiral templates (i.e., a template list) */

/*  <lalVerbatim file="LALInspiralBankHS"> */
typedef struct
tagInspiralTemplateList
{
  INT4 ID;
  InspiralTemplate  params;       /* Pointer to parameter vectors */
  InspiralMetric    metric;       /* Pointer to metric at every point */
  struct tagInspiralTemplateList *next;  /* to create linked list */
}
InspiralTemplateList;
/* </lalVerbatim>  */
 
/*  <lalLaTeX>
\index{\texttt{InspiralTemplateList}}
</lalLaTeX>  */

/* Parameters needed in InspiralCreateCoarseBank */

/*  <lalVerbatim file="LALInspiralBankHS"> */
typedef struct
tagInspiralBankParams
{
   INT4 nparams;            /* for future use, presently 2-dimensional */
   REAL8 x0;                /* coordinates and increments at current location */
   REAL8 x1;
   REAL8 dx0;
   REAL8 dx1;
   REAL8 x0Min;             /* min and max values of parameters */
   REAL8 x0Max;
   REAL8 x1Min;
   REAL8 x1Max;
   InspiralMetric *metric;  /* metric at current location */
}
InspiralBankParams;
/* </lalVerbatim>  */
 
/*  <lalLaTeX>
\index{\texttt{InspiralBankParams}}
</lalLaTeX>  */

/* input for specifying a template bank */

/*  <lalVerbatim file="LALInspiralBankHS"> */
typedef struct
tagInspiralCoarseBankIn
{
  REAL8           mMin;           /* minimum mass of components to search for */
  REAL8           MMax;           /* maximum total mass of binary to search for */
  REAL8           mmCoarse;       /* Coarse grid minimal match */
  REAL8           mmFine;         /* Fine grid minimal match */
  REAL8           fLower;         /* Lower frequency cutoff */
  REAL8           fUpper;         /* Upper frequency cutoff */
  REAL8           tSampling;      /* Sampling rate */
  Detector        detector;       /* The detector */
  Method          method;         /* Method of waveform generation */
  INT4            order;          /* Post-Newtonian order of the waveform */
  Approximant     approximant;    /* Approximant of the waveform */
  Domain          domain;         /* Time- or Frequency-domain*/
  CoordinateSpace space;          /* which of t0-t2 or t0-t3 coordinates */
  REAL8 etamin;                   /* minimum value of eta in our search */
  INT4 iflso;                     /* flso will be used as an upper limit in
                                     moments integrals if iflso!=0; else 
                                     fUpper will be used */
}
InspiralCoarseBankIn;
/* </lalVerbatim>  */

/*  <lalLaTeX>
\index{\texttt{InspiralCoarseBankIn}}
</lalLaTeX>  */

typedef 
enum {ta, phi, t0, t1, t2, t3, t4, t5} 
DParamSelect;

typedef struct
tagInspiralWaveDerivativeIn 
{
   DParamSelect     dp;
   REAL8Vector      *psd;
   InspiralTemplate p;
   REAL8            eps;
}
InspiralWaveDerivativeIn;

/*  <lalVerbatim file="LALInspiralBankHS"> */
typedef struct {
   REAL8 xmin, xmax, ndx, norm;
   Detector detector;
} InspiralMomentsIn;
/* </lalVerbatim>  */

/*  <lalLaTeX>
\index{\texttt{InspiralMomentsIn}}
</lalLaTeX>  */


/*  <lalVerbatim file="LALInspiralBankHS"> */
typedef struct {
   REAL8 ndx;
   REAL8 (*NoisePsd)(REAL8 x);
} InspiralMomentsIntegrandIn;
/* </lalVerbatim>  */

/*  <lalLaTeX>
\index{\texttt{InspiralMomentsIntegrandIn}}
</lalLaTeX>  */

/*  <lalVerbatim file="LALInspiralBankHS"> */
typedef struct
tagCorrelateIn {
   REAL4Vector signal1, signal2;
   REAL8Vector psd;
   RealFFTPlan *revp;
}CorrelateIn;
/* </lalVerbatim>  */

/*  <lalLaTeX>
\index{\texttt{CorrelateIn}}
</lalLaTeX>  */

typedef struct
tagInspiralFineBankIn
{
   InspiralTemplateList templateList;
   InspiralCoarseBankIn coarseIn;
} InspiralFineBankIn;


/*  <lalLaTeX>

\vfill{\footnotesize\input{LALInspiralBankHV}}

</lalLaTeX>  */

/* Function prototypes */

/*  <lalLaTeX>
\newpage\input{LALInspiralComputeParamsC}
</lalLaTeX>  */

void 
LALInspiralComputeParams(
   LALStatus            *status,
   InspiralTemplate     *pars,
   InspiralBankParams   bankParams,
   InspiralCoarseBankIn coarseIn
);


void 
LALInspiralCopyBankParams(
   LALStatus             *status,
   InspiralBankParams *bankParams1,
   InspiralBankParams bankParams2
);


void 
LALInspiralFirstTemplateOnBoundary(
   LALStatus               *status,
   InspiralBankParams   *bankParams,
   InspiralCoarseBankIn bankIn
);

 
void 
LALInspiralLastTemplateOnBoundary(
   LALStatus               *status,
   InspiralBankParams   *bankParams,
   InspiralCoarseBankIn bankIn
);
 

/* <lalLaTeX>
\newpage\input{LALInspiralValidParamsC}
</lalLaTeX>  */

void 
LALInspiralValidParams(
   LALStatus            *status,
   INT4                 *valid,
   InspiralBankParams   bankParams,
   InspiralCoarseBankIn coarseIn
);

/*  <lalLaTeX>
\newpage\input{LALInspiralSetSearchLimitsC}
</lalLaTeX>  */

void 
LALInspiralSetSearchLimits(
   LALStatus               *status,
   InspiralBankParams   *bankParams,
   InspiralCoarseBankIn coarseIn
);

/* <lalLaTeX>
\newpage\input{LALInspiralCreateCoarseBankC}
</lalLaTeX>  */

void 
LALInspiralCreateCoarseBank(
   LALStatus              *status,
   InspiralTemplateList   *list,
   INT4                   *nlist,
   InspiralCoarseBankIn   bankIn
);

/* <lalLaTeX>
\newpage\input{LALInspiralCreateFineBankC}
</lalLaTeX>  */

void 
LALInspiralCreateFineBank(
   LALStatus              *status,
   InspiralTemplateList   *outlist,
   INT4                   *nlist,
   InspiralFineBankIn     fineIn
);

/* <lalLaTeX>
\newpage\input{LALInspiralComputeMetricC}
</lalLaTeX>  */

void 
LALInspiralComputeMetric(
   LALStatus           *status,
   InspiralMetric      *metric,
   InspiralTemplate    params,
   INT4                pass
);

/* <lalLaTeX>
\newpage\input{LALInspiralUpdateParamsC}
</lalLaTeX>  */

void 
LALInspiralUpdateParams(
   LALStatus            *status,
   InspiralBankParams   *bankParams,
   InspiralMetric       metric,
   REAL8 minimalMatch
);

/* <lalLaTeX>
%\newpage\input{LALInspiralFineGridSpacingC}
</lalLaTeX>  */

void 
LALInspiralFineGridSpacing(
   LALStatus            *status,
   InspiralBankParams   *bankParams,
   InspiralMetric       metric,
   InspiralCoarseBankIn coarseIn
);

void 
LALInspiralWaveLength(
   LALStatus           *status, 
   INT4             *length,
   InspiralTemplate p
);

void 
LALNoiseSpectralDensity (
   LALStatus   *status, 
   REAL8Vector *psd, 
   Detector    choice, 
   REAL8       f
);

void 
LALCorrelate (
   LALStatus   *status, 
   REAL4Vector *output, 
   CorrelateIn in
);

void 
LALNormalise (
   LALStatus   *status, 
   REAL4Vector *dh, 
   REAL8       *norm, 
   REAL8Vector psd
);

/* <lalLaTeX>
\newpage\input{LALGEOPsdC}
</lalLaTeX>  */

REAL8 LALGEOPsd (REAL8 x);

/* <lalLaTeX>
\newpage\input{LALLIGOIPsdC}
</lalLaTeX>  */

REAL8 LALLIGOIPsd (REAL8 x);

/* <lalLaTeX>
\newpage\input{LALTAMAPsdC}
</lalLaTeX>  */

REAL8 LALTAMAPsd (REAL8 x);

/* <lalLaTeX>
\newpage\input{LALVIRGOPsdC}
</lalLaTeX>  */

REAL8 LALVIRGOPsd (REAL8 x);


void InverseMatrix (
LALStatus *status,
   INT4   Dim, 
   REAL8  **mm3, 
   INT4   *buff2, 
   REAL8  **buff1);

/* <lalLaTeX>
\newpage\input{LALMatrixTransformC}
</lalLaTeX>  */

void 
LALMatrixTransform (
   LALStatus *status,
   INT4  Dim,
   REAL8 **trans,
   REAL8 **buff1,
   REAL8 **mm3);

/* <lalLaTeX>
\newpage\input{LALInspiralMomentsC}
</lalLaTeX>  */

void 
LALInspiralMoments(
   LALStatus         *status,
   REAL8             *moment,
   InspiralMomentsIn pars);

/* <lalLaTeX>
\newpage\input{LALInspiralMomentsIntegrandC}
</lalLaTeX>  */

void 
LALInspiralMomentsIntegrand(
   LALStatus *status,
   REAL8  *integrand,
   REAL8  f,
   void   *pars);

void LALDeterminant3(LALStatus *status, 
                     REAL8  *determinant, 
                     REAL8  **matrix) ;

/* <lalLaTeX>
\newpage\input{LALInverse3C}
</lalLaTeX>  */

void LALInverse3(LALStatus *status, 
                 REAL8     **inverse, 
                 REAL8     **matrix) ;

/* <lalLaTeX>
\newpage\input{LALInspiralSetParamsC}
</lalLaTeX>  */

void 
LALInspiralSetParams(
   LALStatus            *status, 
   InspiralTemplate     *tempPars,
   InspiralCoarseBankIn coarseIn);


/* <lalLaTeX>
\newpage\input{LALInspiralNextTemplateC}
</lalLaTeX>  */

void 
LALInspiralNextTemplate(
    LALStatus          *status, 
    InspiralBankParams *bankPars, 
    InspiralMetric      metric);            
      
/* <lalLaTeX>
\newpage\input{CoarseTestC}
</lalLaTeX> */

typedef struct 
tagAddVectorsIn
{
   REAL4Vector *v1, *v2;
   REAL8 a1, a2;
} AddVectorsIn;

typedef struct 
tagRandomInspiralSignalIn
{
   InspiralTemplate param;
   UINT8 useed;
   REAL8 mMin, MMax, SignalAmp, NoiseAmp;
   INT4 type;
   REAL8Vector psd;
   RealFFTPlan *fwdp;
} RandomInspiralSignalIn;

typedef struct
tagOverlapIn {
   REAL4Vector signal;
   REAL8Vector psd;
   InspiralTemplate param;
   RealFFTPlan *fwdp, *revp;
} OverlapIn; 

typedef struct
tagOverlapOut {
   REAL8 max, phase;
   INT4 bin;
} OverlapOut; 

void
LALRandomInspiralSignal(
   LALStatus *status, 
   REAL4Vector *signal,
   RandomInspiralSignalIn *randIn);

void
LALGaussianNoise (
   LALStatus   *status,
   REAL4Vector *noisy, 
   UINT8       *seed);

void
LALColoredNoise (
   LALStatus   *status,
   REAL4Vector *noisy, 
   REAL8Vector psd);

void
LALGaussianRandomNumber(
   LALStatus *status, 
   REAL4 *randnum);

void 
LALAddVectors(
   LALStatus *status, 
   REAL4Vector *vector, 
   AddVectorsIn in);

void
LALInspiralWaveOverlap (
   LALStatus   *status,
   REAL4Vector *output,
   OverlapOut  *overlapout,
   OverlapIn   *overlapin);

#ifdef  __cplusplus
}
#endif

#endif /* _LALINSPIRALBANK_H */
