/*
*  Copyright (C) 2007 B.S. Sathyaprakash
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

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
{\bf The above sum is limited to frequency} {\tt in->fCutoff.}
Also, note that the zeroth and Nyquist frequency components
are ignored in the computation of the norm.
Moreover, {\bf array elements of} {\tt filter} corresponding
to frequencies greater than {\tt in->fCutoff} are {\bf set to zero}.
That is, the code replaces the original vector $H_k$ with {\it normalized
vector} using:
\begin{eqnarray}
\widehat H_k & = & \frac {H_k}{\sqrt H},
\ \ \ {\tt k \times in\rightarrow df} \le {\tt in\rightarrow fCutoff},\nonumber \\
& = & 0, \ \ \ {\tt k \times in\rightarrow df} > {\tt in\rightarrow fCutoff}.
\end{eqnarray}
In addition, the 0th and Nyquist frequency components are also set to zero.
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
	  if (f>in->fCutoff)
	  {
		  /* Since the normalisation is done using power only up to fCutoff
		   * it is better to terminate the frequency-domain waveform beyond
		   * fCutoff so that the signal doesn't contribute to the correlation
		   * outside this frequency band
		   */
		  filter->data[i] = filter->data[n-i] = 0.;
	  }
	  else
	  {
		  psd = in->psd->data[i];
		  if (psd)
		  {
			  *norm += (filter->data[i]*filter->data[i] + filter->data[n-i]*filter->data[n-i])/(psd*0.5);
		  }
	  }
  }

  /* If we want to add the negative frequency components as well now is
   * the time to do that before 0th and Nyquist contributions are added
   */

  (*norm) *= 2.;

  /* Set the 0th and Nyquist frequency bins to be zero. */
  filter->data[0] = filter->data[nby2] = 0.;

  *norm /= ((double) n * in->samplingRate);
  *norm = sqrt(*norm);

  for (i=0; i<n; i++)
  {
	  *(filter->data+i) /= *norm;
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

