/*  <lalVerbatim file="LALInspiralWaveNormaliseLSOCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */
/* <lalLaTeX>
\subsection{Module \texttt{LALInspiralWaveNormaliseLSO.c}}
Module to find the norm of a signal and to return a normaliseLSOd 
array. The original signal is left untouched.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralWaveNormaliseLSOCP}
\index{\verb&LALInspiralWaveNormaliseLSO()&}

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

\vfill{\footnotesize\input{LALInspiralWaveNormaliseLSOCV}}
</lalLaTeX>  */
#include <lal/LALNoiseModels.h>

NRCSID (LALINSPIRALWAVENORMALISEC, "$Id$");

/*  <lalVerbatim file="LALInspiralWaveNormaliseLSOCP"> */
void
LALInspiralWaveNormaliseLSO 
   (
   LALStatus   *status, 
   REAL4Vector *filter,
   REAL8       *norm,
   InspiralWaveNormaliseIn *in 
   ) 
{  /*  </lalVerbatim>  */

  INT4 i, n, nby2;
  REAL8 psd, f;

  INITSTATUS (status, "LALInspiralWaveNormaliseLSO", LALINSPIRALWAVENORMALISEC);
  ATTATCHSTATUSPTR(status);

  ASSERT (filter->data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
  ASSERT (in->psd->data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
  ASSERT (in->psd->length == filter->length/2+1, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);

  n = filter->length;
  *norm = 0;
  nby2 = n/2;

  for (i=1; i<nby2; i++) 
  {
     f = i*in->df;
     if (f>in->fCutoff) break;
     psd = in->psd->data[i];
     if (psd) 
     {
	*norm += (pow(filter->data[i], 2.) + pow(filter->data[n-i], 2.))/(psd*0.5);
     }
  }

  /* If we want to add the negative frequency components as well now is
   * the time to do that before 0th and Nyquist contributions are added
   */

  (*norm) *= 2.;

  if (in->psd->data[0]) *norm += pow(filter->data[0], 2.)/(0.5*in->psd->data[0]);
  f = nby2*in->df;
  if (f<in->fCutoff) if (in->psd->data[nby2]) *norm += pow(filter->data[nby2], 2.)/(0.5*in->psd->data[nby2]);

  *norm /= ((double) n * in->samplingRate);
  *norm = sqrt(*norm);
  
  for (i=0; i<n; i++)
  { 
	  *(filter->data+i) /= *norm;
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

