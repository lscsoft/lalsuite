/*  <lalVerbatim file="LALColoredNoiseCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/* <lalLaTeX>
\subsection{Module \texttt{LALColoredNoise.c}}
This module colors a given white noise input into a colored noise
of power spectral density \texttt{psd}.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALColoredNoiseCP}
\index{\verb&LALColoredNoise()&}

\subsubsection*{Description}
Given the Fourier transform $N(f)$ of  white noise, the
Fourier transform of noise of power spectral density $S(f)$ is
given by ${\cal N}(f) = N(f) \times \sqrt{S(f)}.$
In the discrete version there is an additional normalisation:
$${\cal N}_k = N_k \times \sqrt{\frac{2 S_k}{n}},\ \ 
  {\cal N}_{n-k} = N_{n-k} \times \sqrt{\frac{2 S_k}{n}},\ \ 
  k=1, \ldots, \frac{n}{2}.$$

\subsubsection*{Algorithm}
\subsubsection*{Uses}
\begin{verbatim}
none
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALColoredNoiseCV}}
</lalLaTeX>  */

#include <lal/LALNoiseModels.h>
NRCSID (LALCOLOREDNOISEC, "$Id$");

/*  <lalVerbatim file="LALColoredNoiseCP"> */
void 
LALColoredNoise 
   (
   LALStatus   *status,
   REAL4Vector *noisy, 
   REAL8Vector  psd
   ) 
{  /*  </lalVerbatim>  */

   INT4 i, j, n, nby2;
   REAL8 x, length;
   

   INITSTATUS (status, "LALColoredNoise", LALCOLOREDNOISEC);
   ATTATCHSTATUSPTR(status);

   ASSERT (noisy,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (noisy->data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (psd.data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (psd.length == noisy->length/2+1, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);

   n = length = noisy->length;
   nby2 = n/2;

   for (i=1; i<nby2; i++) 
   {
      j = n-i;
      // Since fftw requires n and NOT n/2, I presume we
      // don't need the factor 2 in the normalisation
      // x = sqrt(2. * psd.data[i] / length);
      x = sqrt(psd.data[i] / length);
      noisy->data[i] *= x;
      noisy->data[j] *= x;
   }
   x = sqrt(psd.data[0] / length);
   noisy->data[0] *= x;
   x = sqrt(psd.data[nby2] / length);
   noisy->data[nby2] *= x;

   DETATCHSTATUSPTR(status);
   RETURN (status);
}
