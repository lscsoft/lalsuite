/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpView.h
 *
 * Author: Brady, P.R, and Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpViewHV">
Author: Brady P., R., and Brown, D. A.
$Id$
</lalVerbatim> 

<lalLaTeX>
\section{Header \texttt{FindChirpView.h}}
\label{s:FindChirpView.h}

Provides protypes, structures and functions to allow visualisation of
the events generated \texttt{findchirp} and the \texttt{inspiral} shared
object.

\subsection*{Synopsis}

\begin{verbatim}
#include <lal/FindChirpView.h>
\end{verbatim}

</lalLaTeX>
#endif

#ifndef _FINDCHIRPVIEWH_H
#define _FINDCHIRPVIEWH_H

#include <lal/LALDatatypes.h>
#include <lal/LALInspiral.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif


NRCSID (FINDCHIRPVIEWH, "$Id$");

#if 0
<lalLaTeX> 
\newpage\subsection*{Error codes} 
</lalLaTeX>
#endif
/* <lalErrTable> */
#define FINDCHIRPVIEWH_ENULL 1
#define FINDCHIRPVIEWH_ENNUL 2
#define FINDCHIRPVIEWH_EALOC 3
#define FINDCHIRPVIEWH_MSGENULL "Null pointer"
#define FINDCHIRPVIEWH_MSGENNUL "Non-null pointer"
#define FINDCHIRPVIEWH_MSGEALOC "Memory allocation error"
/* </lalErrTable> */


/*
 *
 * typedefs of structures used by findchip view functions
 *
 */


#if 0
<lalLaTeX>
\subsection*{Types}
</lalLaTeX>
#endif

typedef struct
tagInspiralLDASJob
{
  INT4                          jobID;
  INT4                          numEvents;
  LIGOTimeGPS                   startGPS;
  LIGOTimeGPS                   stopGPS;
  REAL4                         mpiTime;
  REAL4                         totalTime;
  struct tagInspiralLDASJob    *next;
  struct tagInspiralLDASJob    *prev;
}
InspiralLDASJob;

#if 0
<lalLaTeX>
\vfill{\footnotesize\input{FindChirpViewHV}}
</lalLaTeX> 
#endif


/*
 *
 * function prototypes
 *
 */


LALFindChirpInspiralTimingGraph (
    LALStatus          *status,
    InspiralLDASJob    *jobList,
    const CHAR         *fileName
    );

#if 0
<lalLaTeX>
\newpage\input{FindChirpViewC}
</lalLaTeX> 
#endif

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _FINDCHIRPVIEWH_H */
