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
\index{\texttt{LALAddVectors()}}

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

   i=vector->length;
   if ( (in.v1->length != in.v2->length) ||
        (in.v1->length != vector->length) ) {
      fprintf(stderr, "LALAddVectors: Incompatible lengths\n");
      exit(1);
   }
   while (i--)
      vector->data[i] = in.a1 * in.v1->data[i] + in.a2 * in.v2->data[i];
   RETURN(status);
}
