/***************************** <lalVerbatim file="DestroyZPGFilterCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{DestroyZPGFilter.c}}
\label{ss:DestroyZPGFilter.c}

Destroys ZPG filter objects.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{DestroyZPGFilterCP}
\index{\texttt{LALDestroyCOMPLEX8ZPGFilter()}}
\index{\texttt{LALDestroyCOMPLEX16ZPGFilter()}}

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

******************************************************* </lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/ZPGFilter.h>

NRCSID(DESTROYZPGFILTERC,"$Id$");

/* <lalVerbatim file="DestroyZPGFilterCP"> */
void
LALDestroyCOMPLEX8ZPGFilter( LALStatus         *stat,
			     COMPLEX8ZPGFilter **input )
{ /* </lalVerbatim> */
  INITSTATUS(stat,"LALDestroyCOMPLEX8ZPGFilter",DESTROYZPGFILTERC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure handle is non-null, and points to a non-null
     pointer. */
  ASSERT(input,stat,ZPGFILTERH_ENUL,ZPGFILTERH_MSGENUL);
  ASSERT(*input,stat,ZPGFILTERH_ENUL,ZPGFILTERH_MSGENUL);

  /* Destroy the vector fields (if they exist). */
  if((*input)->zeros) {
    TRY(LALCDestroyVector(stat->statusPtr,&((*input)->zeros)),stat);
  }
  if((*input)->poles) {
    TRY(LALCDestroyVector(stat->statusPtr,&((*input)->poles)),stat);
  }

  /* Free the filter, then point the handle to NULL. */
  LALFree(*input);
  *input=NULL;

  /* Normal exit */
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}


/* <lalVerbatim file="DestroyZPGFilterCP"> */
void
LALDestroyCOMPLEX16ZPGFilter( LALStatus          *stat,
			      COMPLEX16ZPGFilter **input )
{ /* </lalVerbatim> */
  INITSTATUS(stat,"LALDestroyCOMPLEX16ZPGFilter",DESTROYZPGFILTERC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure handle is non-null, and points to a non-null
     pointer. */
  ASSERT(input,stat,ZPGFILTERH_ENUL,ZPGFILTERH_MSGENUL);
  ASSERT(*input,stat,ZPGFILTERH_ENUL,ZPGFILTERH_MSGENUL);

  /* Destroy the vector fields. */
  if((*input)->zeros) {
    TRY(LALZDestroyVector(stat->statusPtr,&((*input)->zeros)),stat);
  }
  if((*input)->poles) {
    TRY(LALZDestroyVector(stat->statusPtr,&((*input)->poles)),stat);
  }

  /* Free the filter, then point the handle to NULL. */
  LALFree(*input);
  *input=NULL;

  /* Normal exit */
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}
