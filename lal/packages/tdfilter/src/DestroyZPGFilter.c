/***************************** <lalVerbatim file="DestroyZPGFilterCV">
Author: Creighton, T. D.
$Id$
****************************** </lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{DestroyZPGFilter.c}}
\label{ss:DestroyZPGFilter.c}

Destroys ZPG filter objects.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{DestroyZPGFilterCP}
\index{\verb&DestroyCOMPLEX8ZPGFilter()&}
\index{\verb&DestroyCOMPLEX16ZPGFilter()&}

\subsubsection*{Description}

These functions destroy an object \verb@**output@ of type
\verb@COMPLEX8ZPGFilter@ or \verb@COMPLEX16ZPGFilter@, and set
\verb@*output@ to \verb@NULL@.

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
LALFree()
CDestroyVector()
ZDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{DestroyZPGFilterCV}}

</lalLaTeX> */

#include "LALStdlib.h"
#include "AVFactories.h"
#include "ZPGFilter.h"

NRCSID(DESTROYZPGFILTERC,"$Id$");

/* <lalVerbatim file="DestroyZPGFilterCP"> */
void DestroyCOMPLEX8ZPGFilter(Status            *stat,
			      COMPLEX8ZPGFilter **input)
{ /* </lalVerbatim> */
  INITSTATUS(stat,"DestroyCOMPLEX8ZPGFilter",DESTROYZPGFILTERC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure handle is non-null, and points to a non-null pointer.
     (The routine CDestroyVector will check that the data fields are
     non-null.) */
  ASSERT(input,stat,ZPGFILTER_ENUL,ZPGFILTER_MSGENUL);
  ASSERT(*input,stat,ZPGFILTER_ENUL,ZPGFILTER_MSGENUL);

  /* Destroy the vector fields. */
  TRY(CDestroyVector(stat->statusPtr,&((*input)->zeros)),stat);
  TRY(CDestroyVector(stat->statusPtr,&((*input)->poles)),stat);

  /* Free the filter, then point the handle to NULL. */
  LALFree(*input);
  *input=NULL;

  /* Normal exit */
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}


/* <lalVerbatim file="DestroyZPGFilterCP"> */
void DestroyCOMPLEX16ZPGFilter(Status             *stat,
			       COMPLEX16ZPGFilter **input)
{ /* </lalVerbatim> */
  INITSTATUS(stat,"DestroyCOMPLEX16ZPGFilter",DESTROYZPGFILTERC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure handle is non-null, and points to a non-null pointer.
     (The routine ZDestroyVector will check that the data fields are
     non-null.) */
  ASSERT(input,stat,ZPGFILTER_ENUL,ZPGFILTER_MSGENUL);
  ASSERT(*input,stat,ZPGFILTER_ENUL,ZPGFILTER_MSGENUL);

  /* Destroy the vector fields. */
  TRY(ZDestroyVector(stat->statusPtr,&((*input)->zeros)),stat);
  TRY(ZDestroyVector(stat->statusPtr,&((*input)->poles)),stat);

  /* Free the filter, then point the handle to NULL. */
  LALFree(*input);
  *input=NULL;

  /* Normal exit */
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}
