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
Given the positive frequency Fourier components 
$H_k,$ $k=0,\ldots,n-1,$ of a vector 
and the noise PSD $S_m,$ $m=0,\ldots,n/2,$ 
this module first computes the norm $H$ of the vector treating
$S_m$ as the measure: 
(note that in {\em fftw} notation, the zeroth frequency 
component is $H_0,$ Nyquist 
is $H_{n/2},$ $H_k,$ $k \ne 0,n/2,$ ($H_{n-k})$ is the real (imaginary)
part of the $k$th harmonic) 
\begin{equation}
H = \sum_{k=1}^{n/2-1} \frac{H_k^2 + H^2_{n-k}}{S_k}.
\label{eq:inspiralnorm}
\end{equation}
(Note that the zeroth and Nyquist components are ignored in the
computation of the norm.)
It then replaces the original vector $H_k$ with normalized
vector using: 
\begin{equation}
\widehat H_k = \frac {H_k}{\sqrt H},\ \ k=0,\ldots n-1.
\end{equation}
Note that the norm of $\widehat H_k,$ as defined in Eq.~(\ref{eq:inspiralnorm}), 
is unity.
\subsubsection*{Algorithm}
\subsubsection*{Uses}
\begin{verbatim}
none.
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralWaveNormaliseCV}}
</lalLaTeX>  */
#include <lal/LALNoiseModels.h>

NRCSID (LALINSPIRALWAVENORMALISEC, "$Id$");

/*  <lalVerbatim file="LALInspiralWaveNormaliseCP"> */
void
LALInspiralWaveNormalise 
   (
   LALStatus   *status, 
   REAL4Vector *in, 
   REAL8       *norm,
   REAL8Vector psd
   ) 
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

  for (i=1; i<nby2; i++) 
  {
     k = n-i;
     psdvalue = psd.data[i];
     if (psdvalue) 
     {
	*norm += (pow(in->data[i], 2.) + pow(in->data[k], 2.))/psdvalue;
     }
  }
  *norm += pow(in->data[0],2.) + pow(in->data[nby2],2.);
  *norm = sqrt(*norm);
  
  for (i=0; i<n; i++)
  { 
	  *(in->data+i) /= *norm;
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

