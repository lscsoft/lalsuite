/*  <lalVerbatim file="LALStatsREAL4VectorCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */
/* <lalLaTeX>
\subsection{Module \texttt{LALStatsREAL4Vector.c}}
Module to compute the mean, rms, minimum and maximum
of a \texttt{REAL4Vector}.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALStatsREAL4VectorCP}
\index{\verb&LALStatsREAL4Vector()&}

\subsubsection*{Description}
\subsubsection*{Algorithm}
\subsubsection*{Uses}
None
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALStatsREAL4VectorCV}}
</lalLaTeX>  */
#include <lal/LALNoiseModels.h>

NRCSID (LALSTATSREAL4VECTORC, "$Id$");

/*  <lalVerbatim file="LALStatsREAL4VectorCP"> */ 
void
LALStatsREAL4Vector(
   LALStatus *status,
   StatsREAL4VectorOut *out,
   REAL4Vector *vector)

{  /*  </lalVerbatim>  */

   INT4 n;
   REAL8 x;

   INITSTATUS (status, "LALStatsREAL4Vector", LALSTATSREAL4VECTORC);
   ATTATCHSTATUSPTR(status);

   ASSERT (vector->data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (vector->length > 0,  status, LALNOISEMODELSH_ECHOICE, LALNOISEMODELSH_MSGECHOICE);
   ASSERT (out,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);

   out->mean = 0.;
   n = vector->length;
   while (--n)
   {
      out->mean+=vector->data[n];
   }
   out->mean/=(REAL8) vector->length;
   ASSERT (out->mean > 0,  status, LALNOISEMODELSH_ECHOICE, LALNOISEMODELSH_MSGECHOICE);
   
   n = vector->length;
   while (--n)
   {
      x = vector->data[n]-out->mean;
      out->var+=x*x;
   }
   out->var /=(REAL8) vector->length;
   out->stddev =sqrt(out->var);
   DETATCHSTATUSPTR(status);
   RETURN(status);
}
