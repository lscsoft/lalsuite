/*  <lalVerbatim file="LALAddVectorsCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/* <lalLaTeX>
\subsection{Module \texttt{LALAddVectors.c}}
Module to add two vectors with weights as in \texttt{AddVectorsIn}
and return the resulting vector in \texttt{vector}.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALAddVectorsCP}
\index{\verb&LALAddVectors()&}

\subsubsection*{Description}
\subsubsection*{Algorithm}
\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALAddVectorsCV}}
</lalLaTeX>  */
#include <lal/LALNoiseModels.h>
NRCSID (LALADDVECTORSC, "$Id$");

/*  <lalVerbatim file="LALAddVectorsCP"> */
void 
LALAddVectors(
   LALStatus *status, 
   REAL4Vector *vector, 
   AddVectorsIn in)
{  /*  </lalVerbatim>  */
   INT4 i;

   INITSTATUS (status, "LALAddVectors", LALADDVECTORSC);
   ATTATCHSTATUSPTR(status);
   ASSERT (vector,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (vector->length >= 2, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
   ASSERT (vector->data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (vector->data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (in.v1->length == in.v2->length, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
   ASSERT (in.v1->length == vector->length, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);

   i=vector->length;
   while (i--)
      vector->data[i] = in.a1 * in.v1->data[i] + in.a2 * in.v2->data[i];

   DETATCHSTATUSPTR(status);
   RETURN (status);
}
