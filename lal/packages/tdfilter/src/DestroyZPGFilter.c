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
\index{\verb&LALDestroyCOMPLEX8ZPGFilter()&}
\index{\verb&LALDestroyCOMPLEX16ZPGFilter()&}

\subsubsection*{Description}

These functions destroy an object \verb@**output@ of type
\verb@COMPLEX8ZPGFilter@ or \verb@COMPLEX16ZPGFilter@, and set
\verb@*output@ to \verb@NULL@.

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
LALFree()
LALCDestroyVector()
LALZDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{DestroyZPGFilterCV}}

</lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/ZPGFilter.h>

NRCSID(DESTROYZPGFILTERC,"$Id$");

/* <lalVerbatim file="DestroyZPGFilterCP"> */
void LALDestroyCOMPLEX8ZPGFilter(LALStatus            *stat,
			      COMPLEX8ZPGFilter **input)
{ /* </lalVerbatim> */
  INITSTATUS(stat,"LALDestroyCOMPLEX8ZPGFilter",DESTROYZPGFILTERC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure handle is non-null, and points to a non-null pointer.
     (The routine LALCDestroyVector will check that the data fields are
     non-null.) */
  ASSERT(input,stat,ZPGFILTER_ENUL,ZPGFILTER_MSGENUL);
  ASSERT(*input,stat,ZPGFILTER_ENUL,ZPGFILTER_MSGENUL);

  /* Destroy the vector fields. */
  TRY(LALCDestroyVector(stat->statusPtr,&((*input)->zeros)),stat);
  TRY(LALCDestroyVector(stat->statusPtr,&((*input)->poles)),stat);

  /* Free the filter, then point the handle to NULL. */
  LALFree(*input);
  *input=NULL;

  /* Normal exit */
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}


/* <lalVerbatim file="DestroyZPGFilterCP"> */
void LALDestroyCOMPLEX16ZPGFilter(LALStatus             *stat,
			       COMPLEX16ZPGFilter **input)
{ /* </lalVerbatim> */
  INITSTATUS(stat,"LALDestroyCOMPLEX16ZPGFilter",DESTROYZPGFILTERC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure handle is non-null, and points to a non-null pointer.
     (The routine LALZDestroyVector will check that the data fields are
     non-null.) */
  ASSERT(input,stat,ZPGFILTER_ENUL,ZPGFILTER_MSGENUL);
  ASSERT(*input,stat,ZPGFILTER_ENUL,ZPGFILTER_MSGENUL);

  /* Destroy the vector fields. */
  TRY(LALZDestroyVector(stat->statusPtr,&((*input)->zeros)),stat);
  TRY(LALZDestroyVector(stat->statusPtr,&((*input)->poles)),stat);

  /* Free the filter, then point the handle to NULL. */
  LALFree(*input);
  *input=NULL;

  /* Normal exit */
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}
