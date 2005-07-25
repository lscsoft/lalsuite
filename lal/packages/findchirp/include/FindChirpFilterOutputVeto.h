/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpFilterOutputVeto.h
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpFilterOutputVetoHV">
Author: Brown, D. A.
$Id$
</lalVerbatim> 

<lalLaTeX>
\section{Header \texttt{FindChirpFilterOutputVeto.h}}
\label{s:FindChirpFilterOutputVeto.h}

\noindent Provides protypes to implement signal based vetoes other than the
$\chi^2$ veto. No vetoes are presently implemented.

\subsubsection*{Synopsis}

\begin{verbatim}
#include <lal/FindChirpFilterOutputVeto.h>
\end{verbatim}

</lalLaTeX>
#endif

#ifndef _FINDCHIRPFILTEROUTPUTVETOH_H
#define _FINDCHIRPFILTEROUTPUTVETOH_H

#include <lal/LALDatatypes.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/FindChirpDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif


NRCSID (FINDCHIRPFILTEROUTPUTVETOH, "$Id$");

#if 0
<lalLaTeX> 
\subsection*{Error codes} 
</lalLaTeX>
#endif
/* <lalErrTable> */
#define FINDCHIRPFILTEROUTPUTVETOH_ENULL 1 
#define FINDCHIRPFILTEROUTPUTVETOH_ENNUL 2 
#define FINDCHIRPFILTEROUTPUTVETOH_EALOC 3 
#define FINDCHIRPFILTEROUTPUTVETOH_ERSQT 4 
#define FINDCHIRPFILTEROUTPUTVETOH_EDTZO 5 
#define FINDCHIRPFILTEROUTPUTVETOH_ERSQW 6 

#define FINDCHIRPFILTEROUTPUTVETOH_MSGENULL "Null pointer" 
#define FINDCHIRPFILTEROUTPUTVETOH_MSGENNUL "Non-null pointer" 
#define FINDCHIRPFILTEROUTPUTVETOH_MSGEALOC "Memory allocation error" 
#define FINDCHIRPFILTEROUTPUTVETO_MSGERSQT "r^2 threshold is negative" 
#define FINDCHIRPFILTEROUTPUTVETO_MSGEDTZO "deltaT is zero or negative" 
#define FINDCHIRPFILTEROUTPUTVETO_MSGERSQW "r^2 veto window is zero or negative" 
/* </lalErrTable> */


#if 0
<lalLaTeX>
\subsection*{Types}

\subsubsection*{Structure \texttt{FindChirpFilterOutputVetoParams}}
\idx[Type]{FindChirpFilterOutputVetoParams}

\noindent This structure provides the parameters for the filter
output veto.

</lalLaTeX>
#endif
/* <lalVerbatim> */
typedef struct
tagFindChirpFilterOutputVetoParams
{
  INT4          rsqvetoWindow;
  INT4          rsqvetoThresh;
}
FindChirpFilterOutputVetoParams;
/* </lalVerbatim> */
#if 0
<lalLaTeX>

\begin{description}
\item[\texttt{UINT4 rsqvetoWindow}] Width of the veto window.

\item[\texttt{UINT4 rsqvetoThresh}] Threshold of the veto window.
\end{description}

\vfill{\footnotesize\input{FindChirpFilterOutputVetoHV}}
</lalLaTeX> 
#endif

#if 0
<lalLaTeX>
\newpage\input{FindChirpFilterOutputVetoC}
</lalLaTeX>
#endif

void LALFindChirpFilterOutputVeto( 
    LALStatus                          *status,
    SnglInspiralTable                 **eventList, 
    FindChirpSegment                   *segment,
    REAL4Vector                        *chisqVec,
    REAL8                              deltaT,
    FindChirpFilterOutputVetoParams    *params
    );

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _FINDCHIRPFILTEROUTPUTVETO_H */
