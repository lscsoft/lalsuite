/***************************** <lalVerbatim file="DestroyIIRFilterCV">
Author: Creighton, T. D.
$Id$
****************************** </lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{DestroyIIRFilter.c}}
\label{ss:DestroyIIRFilter.c}

Destroys IIR filter objects.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{DestroyIIRFilterCP}
\index{\verb&DestroyREAL4IIRFilter()&}
\index{\verb&DestroyREAL8IIRFilter()&}

\subsubsection*{Description}

These functions destroy an object \verb@**input@ of type
\texttt{REAL4IIRFilter} or \texttt{REAL8IIRFilter}, and set
\verb@*input@ to \verb@NULL@.

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
void LALFree()
void SDestroyVector()
void DDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{DestroyIIRFilterCV}}

</lalLaTeX> */

#include "LALStdlib.h"
#include "AVFactories.h"
#include "IIRFilter.h"

NRCSID(DESTROYIIRFILTERC,"$Id$");


/* <lalVerbatim file="DestroyIIRFilterCP"> */
void DestroyREAL4IIRFilter(Status         *stat,
			   REAL4IIRFilter **input)
{ /* </lalVerbatim> */
  INITSTATUS(stat,"DestroyREAL4IIRFilter",DESTROYIIRFILTERC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure handle is non-null, and points to a non-null pointer.
     (The routine SDestroyVector will check that the data fields are
     non-null.) */
  ASSERT(input,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(*input,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);

  /* Destroy the vector fields. */
  TRY(SDestroyVector(stat->statusPtr,&((*input)->directCoef)),stat);
  TRY(SDestroyVector(stat->statusPtr,&((*input)->recursCoef)),stat);
  TRY(SDestroyVector(stat->statusPtr,&((*input)->history)),stat);

  /* Free the filter, then point the handle to NULL. */
  LALFree(*input);
  *input=NULL;

  /* Normal exit */
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}


/* <lalVerbatim file="DestroyIIRFilterCP"> */
void DestroyREAL8IIRFilter(Status         *stat,
			   REAL8IIRFilter **input)
{ /* </lalVerbatim> */
  INITSTATUS(stat,"DestroyREAL8IIRFilter",DESTROYIIRFILTERC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure handle is non-null, and points to a non-null pointer.
     (The routine DDestroyVector will check that the data fields are
     non-null.) */
  ASSERT(input,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(*input,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);

  /* Destroy the vector fields. */
  TRY(DDestroyVector(stat->statusPtr,&((*input)->directCoef)),stat);
  TRY(DDestroyVector(stat->statusPtr,&((*input)->recursCoef)),stat);
  TRY(DDestroyVector(stat->statusPtr,&((*input)->history)),stat);

  /* Free the filter, then point the handle to NULL. */
  LALFree(*input);
  *input=NULL;

  /* Normal exit */
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}
