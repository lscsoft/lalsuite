/************************* <lalVerbatim file="DestroyResampleRulesCV">
Author: Creighton, T. D.
Revision: $Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{DestroyResampleRules.c}}
\label{ss:DestroyResampleRules.c}

Destroys an object of type \verb@ResampleRules@.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{DestroyResampleRulesCP}
\index{\texttt{LALDestroyResampleRules()}}

\subsubsection*{Description}

This function destroys an object \verb@**rules@ of type
\texttt{ResampleRules}, and sets \verb@*rules@ to \verb@NULL@.

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
void LALFree()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{DestroyResampleRulesCV}}

******************************************************* </lalLaTeX> */

#include "LALStdlib.h"
#include "AVFactories.h"
#include "Resample.h"

NRCSID(DESTROYRESAMPLERULESC,"$Id$");


/* <lalVerbatim file="DestroyResampleRulesCP"> */
void
LALDestroyResampleRules( LALStatus     *stat,
			 ResampleRules **rules )
{ /* </lalVerbatim> */
  INITSTATUS(stat,"LALDestroyResampleRules",DESTROYRESAMPLERULESC);

  /* Make sure that the handle is non-null, that it points to a
     non-null pointer, and that the interval and shift fields are
     non-null.) */
  ASSERT(rules,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(*rules,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT((*rules)->interval,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT((*rules)->shift,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);

  /* Free all the memory. */
  LALFree((*rules)->interval);
  LALFree((*rules)->shift);
  LALFree(*rules);

  /* Point the handle to NULL, then get out of here. */
  *rules=NULL;
  RETURN(stat);
}
