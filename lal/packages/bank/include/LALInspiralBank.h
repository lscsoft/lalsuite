/* <lalVerbatim file="LALInspiralBankHV">

Author: Churches, D.K. and Sathyaprakash, B.S.
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
#include <lal/LALNoiseModels.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( LALINSPIRALBANKH, "$Id$" );

/* <lalLaTeX>

\subsection*{Error codes}

</lalLaTeX>  */

/* <lalErrTable> */

#define LALINSPIRALBANKH_ENULL 1
#define LALINSPIRALBANKH_EMEM  2
#define LALINSPIRALBANKH_ECHOICE 3
#define LALINSPIRALBANKH_EDIV0 4
#define LALINSPIRALBANKH_ESIZE 8
#define LALINSPIRALBANKH_MSGENULL "Arguments contained an unexpected null pointer"
#define LALINSPIRALBANKH_MSGEMEM "Memory allocation failure"
#define LALINSPIRALBANKH_MSGECHOICE "Invalid choice for an input parameter"
#define LALINSPIRALBANKH_MSGEDIV0 "Division by zero"
#define LALINSPIRALBANKH_MSGESIZE "Invalid input range"


/* </lalErrTable> */

/* <lalLaTeX>

\section*{Structures}
\input{LALInspiralBankHS}
</lalLaTeX>  */

/*  <lalLaTeX> 
\idx[Type]{Detector} 
</lalLaTeX>  */

/*  <lalVerbatim file="LALInspiralBankHS"> */
typedef enum
{
  Tau0Tau2, Tau0Tau3
}
CoordinateSpace;
/*  </lalVerbatim>  */

/*  <lalLaTeX> 
\idx[Type]{CoordinateSpace} 
</lalLaTeX>  */

/*  <lalLaTeX>
\idx[Type]{LALInspiralMetric}
</lalLaTeX>  */

/* Metric and its dimension */

/*  <lalVerbatim file="LALInspiralBankHS"> */
typedef struct 
tagInspiralMetric 
{
   REAL8            g00;     /* 00-component of the diagonalised metric. */
   REAL8            g11;     /* 11-component of the diagonalised metric. */
   REAL8            theta;   /* Angle from t0 to x0 */
   CoordinateSpace space;    /* Coordinate space in which metric is computed */
   INT4 iflso;
   REAL8FrequencySeries *shf; /* one sided strain power spectral density */
} 
InspiralMetric;
/* </lalVerbatim>  */

/*  <lalLaTeX> 
\idx[Type]{InspiralMetric} 
</lalLaTeX>  */

/*  <lalVerbatim file="LALInspiralBankHS"> */
/* a grid of inspiral templates (i.e., a template list) */

typedef struct
tagInspiralTemplateList
{
  INT4 ID;
  InspiralTemplate  params;              /* Pointer to parameter vectors */
  InspiralMetric    metric;              /* Pointer to metric at every point */
  struct tagInspiralTemplateList *next;  /* to create linked list */
}
InspiralTemplateList;
/* </lalVerbatim>  */
 
/*  <lalLaTeX>
\idx[Type]{InspiralTemplateList}
</lalLaTeX>  */


/*  <lalVerbatim file="LALInspiralBankHS"> */
/* Parameters needed in InspiralCreateCoarseBank */
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
\idx[Type]{InspiralBankParams}
</lalLaTeX>  */

/*  <lalVerbatim file="LALInspiralBankHS"> */
/* input for specifying a template bank */

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
  REAL8FrequencySeries shf;
  Order            order;          /* Post-Newtonian order of the waveform */
  Approximant     approximant;    /* Approximant of the waveform */
  CoordinateSpace space;          /* which of t0-t2 or t0-t3 coordinates */
  REAL8 etamin;                   /* minimum value of eta in our search */
  INT4 iflso;                     /* flso will be used as an upper limit in
                                     moments integrals if iflso!=0; else 
                                     fUpper will be used */
}
InspiralCoarseBankIn;
/* </lalVerbatim>  */

/*  <lalLaTeX>
\idx[Type]{InspiralCoarseBankIn}
</lalLaTeX>  */

/*  <lalVerbatim file="LALInspiralBankHS"> */
typedef struct {
   REAL8 xmin, xmax, ndx, norm;
   REAL8FrequencySeries *shf;
} InspiralMomentsIn;
/* </lalVerbatim>  */

/*  <lalLaTeX>
\idx[Type]{InspiralMomentsIn}
</lalLaTeX>  */

/*  <lalVerbatim file="LALInspiralBankHS"> */
typedef struct
tagInspiralFineBankIn
{
   InspiralTemplateList templateList;
   InspiralCoarseBankIn coarseIn;
} InspiralFineBankIn;
/*  </lalVerbatim>  */
/*  <lalLaTeX> 
\idx[Type]{InspiralFineBankIn} 
</lalLaTeX>  */

/*  <lalVerbatim file="LALInspiralBankHS"> */
typedef struct
tagRectangleIn 
   {REAL8 x0, y0, dx, dy, theta;}
RectangleIn;
/*  </lalVerbatim>  */

/*  <lalLaTeX>
\index{\texttt{RectangleIn}}
</lalLaTeX>  */

/*  <lalVerbatim file="LALInspiralBankHS"> */
typedef struct
tagRectangleOut 
   {REAL8 x1, y1, x2, y2, x3, y3, x4, y4, x5, y5;}
RectangleOut;
/*  </lalVerbatim>  */

/*  <lalLaTeX>
\index{\texttt{RectangleOut}}
</lalLaTeX>  */

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

/* <lalLaTeX>
\newpage\input{LALInspiralValidTemplateC}
</lalLaTeX>  */

void 
LALInspiralValidTemplate(
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
   InspiralTemplateList   **list,
   INT4                   *nlist,
   InspiralCoarseBankIn   bankIn
);

/* <lalLaTeX>
\newpage\input{LALInspiralCreateFineBankC}
</lalLaTeX>  */

void 
LALInspiralCreateFineBank(
   LALStatus              *status,
   InspiralTemplateList   **outlist,
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
% TAKEN OUT BY JOLIEN: THIS IS NOT IN BANK
%\newpage\input{LALInspiralWaveLengthC}
</lalLaTeX>  */

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

/* <lalLaTeX>
\newpage\input{LALDeterminant3C}
</lalLaTeX>  */

void
LALDeterminant3(LALStatus *status, 
                REAL8  *determinant, 
                REAL8  **matrix) ;

/* <lalLaTeX>
\newpage\input{LALInverse3C}
</lalLaTeX>  */

void 
LALInverse3(
            LALStatus *status, 
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
\newpage\input{LALRectangleVerticesC}
</lalLaTeX>  */

void 
LALRectangleVertices(
   LALStatus *status, 
   RectangleOut *out,
   RectangleIn *in
);
      
/* <lalLaTeX>
\newpage\input{CoarseTestC}
</lalLaTeX> */

/* <lalLaTeX>
\newpage\input{CoarseTest2C}
</lalLaTeX> */

/* <lalLaTeX>
\newpage\input{ChirpSpaceC}
</lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif /* _LALINSPIRALBANK_H */

