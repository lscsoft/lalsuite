/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpSP.h
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */


#if 0
<lalVerbatim file="FindChirpSPHV">
Author: Brown, D. A.
$Id$
</lalVerbatim> 

<lalLaTeX>

\section{Header \texttt{FindChirpSP.h}}
\label{s:FindChirpSP.h}

Provides routines to filter IFO data for binary inspiral chirps generated
using the stationary phase approximation. 
\subsection*{Synopsis}
\begin{verbatim}
#include "FindChirpSP.h"
\end{verbatim}

\noindent This header provides routines necessary to filter IFO 
data contained in a \verb|DataSegment| structure for binary inspiral 
chirps generated using the stationary phase approximation. 

In order to increase efficency, the filtering algorithm is divided 
into three parts:
\begin{itemize}
\item Those that are independent of the template and need to be calculated
      once per data segment.

\item Those that need to be done once per template.

\item Those that need to be done once for each data segment for each
      template.
\end{itemize}
The template independent functions are contained in the module
\verb|FindChirpSPData.c| and template dependent functions are in
the module \verb|FindChirpSPTemplate.c|.

\input{FindChirpSPHDoc}

</lalLaTeX>
#endif

#ifndef _FINDCHIRPSPH_H
#define _FINDCHIRPSPH_H

#include <lal/LALDatatypes.h>
#include <lal/RealFFT.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>

#ifdef  __cplusplus
extern "C" {
#endif

#if 0
<lalLaTeX> 
\subsection*{Error codes} 
</lalLaTeX>
#endif
/* <lalErrTable> */
#define FINDCHIRPSPH_ENULL 1
#define FINDCHIRPSPH_ENNUL 2
#define FINDCHIRPSPH_EALOC 3
#define FINDCHIRPSPH_ENUMZ 4
#define FINDCHIRPSPH_ESEGZ 5
#define FINDCHIRPSPH_EMISM 6
#define FINDCHIRPSPH_EDELT 7
#define FINDCHIRPSPH_EFLOW 8
#define FINDCHIRPSPH_EDYNR 9
#define FINDCHIRPSPH_MSGENULL "Null pointer"
#define FINDCHIRPSPH_MSGENNUL "Non-null pointer"
#define FINDCHIRPSPH_MSGEALOC "Memory allocation error"
#define FINDCHIRPSPH_MSGENUMZ "Invalid number of segments"
#define FINDCHIRPSPH_MSGESEGZ "Invalid number of points in segments"
#define FINDCHIRPSPH_MSGEMISM "Mismatch between number of points in segments"
#define FINDCHIRPSPH_MSGEDELT "deltaT is zero or negative"
#define FINDCHIRPSPH_MSGEFLOW "Low frequency cutoff is negative"
#define FINDCHIRPSPH_MSGEDYNR "Dynamic range scaling is zero or negative"
/* </lalErrTable> */



typedef struct
tagFindChirpSPDataParams
{
  REAL4Vector                  *ampVec;
  RealFFTPlan                  *fwdPlan;
  RealFFTPlan                  *invPlan;
  REAL4Vector                  *vVec;
  REAL4Vector                  *wVec;
  COMPLEX8Vector               *wtildeVec;
  REAL4Vector                  *tmpltPowerVec;
  REAL4                         deltaT;
  REAL4                         fLow;
  REAL4                         dynRange;
  UINT4                         invSpecTrunc;
}
FindChirpSPDataParams;

typedef struct
tagFindChirpSPTmpltParams
{
  REAL4                         deltaT;
  REAL4                         fLow;
  REAL4                         dynRange;
  REAL4Vector                  *xfacVec;
}
FindChirpSPTmpltParams;


void
LALFindChirpSPDataInit (
    LALStatus                  *status,
    FindChirpSPDataParams     **output,
    FindChirpInitParams        *params
    );

void
LALFindChirpSPData (
    LALStatus                  *status,
    FindChirpSegmentVector     *fcSegVec,
    DataSegmentVector          *dataSegVec,
    FindChirpSPDataParams      *params
    );

void
LALFindChirpSPDataFinalize (
    LALStatus                  *status,
    FindChirpSPDataParams     **output
    );


void
LALFindChirpSPTemplateInit (
    LALStatus                  *status,
    FindChirpSPTmpltParams    **output,
    FindChirpInitParams        *params
    );

void
LALFindChirpSPTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *tmplt,
    FindChirpSPTmpltParams     *params
    );

void
LALFindChirpSPTemplateFinalize (
    LALStatus                  *status,
    FindChirpSPTmpltParams    **output
    );


#ifdef  __cplusplus
}
#endif

/*
<lalLaTeX>
\vfill{\footnotesize\input{FindChirpSPHV}}
\newpage\input{FindChirpSPDataC}
</lalLaTeX> 
*/

#endif /* _FINDCHIRPSPH_H */
