/*  <lalVerbatim file="LALInspiralWaveNormaliseCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */
/* <lalLaTeX>
\subsection{Module \texttt{LALInspiralWaveNormalise.c}}
Module to find the norm of a signal and to return a normalised 
array. The original signal is left untouched.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralWaveNormaliseCP}
\index{\verb&LALInspiralWaveNormalise()&}

\subsubsection*{Description}
\subsubsection*{Algorithm}
\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralWaveNormaliseCV}}
</lalLaTeX>  */
#include <lal/LALNoiseModels.h>

NRCSID (LALINSPIRALWAVENORMALISEC, "$Id$");

/*  <lalVerbatim file="LALInspiralWaveNormaliseCP"> */
void
LALInspiralWaveNormalise (
   LALStatus *status, 
   REAL4Vector *in, 
   REAL8 *norm,
   REAL8Vector psd) 
{  /*  </lalVerbatim>  */

  INT4 i, n, nby2, k;
  REAL8 psdvalue;

  INITSTATUS (status, "LALInspiralWaveNormalise", LALINSPIRALWAVENORMALISEC);
  ATTATCHSTATUSPTR(status);

  ASSERT (in->data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
  ASSERT (psd.data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
  ASSERT (psd.length == in->length/2+1, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
  n = in->length;
  *norm = 0;
  nby2 = n/2;
  for (i=1; i<nby2; i++) {
     k = n-i;
     psdvalue = psd.data[i+1];
     if (psdvalue) *norm += (pow(in->data[i],2.)+pow(in->data[k],2.))/psdvalue;
  }
  *norm = sqrt(*norm);
  
  i=n;
  while (i--) *(in->data+i) = *(in->data+i) / *norm;

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

