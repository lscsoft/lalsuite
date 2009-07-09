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

/*  <lalVerbatim file="LALInspiralComputeSNRIntegrandCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/* <lalLaTeX>
\subsection{Module \texttt{LALInspiralComputeSNRIntegrand.c}}
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralComputeSNRIntegrandCP}
\index{\verb&LALInspiralComputeSNRIntegrand()&}

\subsubsection*{Description}
\subsubsection*{Algorithm}
\subsubsection*{Uses}
\begin{verbatim}
LALREAL4VectorFFT
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralComputeSNRIntegrandCV}}
</lalLaTeX>  */

#include <lal/LALNoiseModels.h>

NRCSID (LALINSPIRALWAVECORRELATEC, "$Id$");

/*  <lalVerbatim file="LALInspiralComputeSNRIntegrandCP"> */
void
LALInspiralComputeSNRIntegrand
   (
   LALStatus                *status,
   REAL4Vector              *output,
   InspiralWaveCorrelateIn  corrin,
   InspiralSNRIntegrandParams *params
   )
{  /*  </lalVerbatim>  */
  INT4 n, nby2, i, k;
  REAL8 psd, r1, r2, i1, i2, twoPiLagByN, cr, ci;
  REAL8 rShift, iShift;


  INITSTATUS (status, "LALInspiralComputeSNRIntegrand", LALINSPIRALWAVECORRELATEC);
  ATTATCHSTATUSPTR(status);

  ASSERT (output,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
  ASSERT (output->data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
  ASSERT (corrin.signal1.data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
  ASSERT (corrin.signal2.data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
  ASSERT (corrin.psd.data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
  ASSERT (corrin.signal1.length == corrin.signal2.length, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);
  ASSERT (corrin.psd.length == corrin.signal1.length/2+1, status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);

  n = corrin.signal1.length;
  twoPiLagByN = 2.L * LAL_PI * (double) params->lag /(double)n;


  nby2 = n/2;

  for (i=1; i<nby2; i++)
  {
     k=n-i;
     psd = corrin.psd.data[i];

     if (psd) {

        r1 = corrin.signal1.data[i];
        r2 = corrin.signal2.data[i];
        i1 = corrin.signal1.data[k];
        i2 = corrin.signal2.data[k];
	rShift = cos((double)i * twoPiLagByN - params->phase);  /* cos(2 \pi f_k \tau) */
	iShift = sin((double)i * twoPiLagByN - params->phase);  /* sin(2 \pi f_k \tau) */

	cr = r1*r2 + i1*i2;
	ci = i1*r2 - r1*i2;

	output->data[i] = (float) ((rShift*cr - iShift*ci)/psd);
	output->data[k] = (float) ((rShift*ci + iShift*cr)/psd);

     } else {

        output->data[i] = 0.L;
        output->data[k] = 0.L;
     }
  }
  psd = corrin.psd.data[0];
  if (psd)
  {
     r1 = corrin.signal1.data[0];
     r2 = corrin.signal2.data[0];
     output->data[0] = (float)(r1*r2 / psd);
     /*
      * printf("%d %e %e\n", i, output->data[0], output->data[0]);
      */
  }
  else
  {
     output->data[0] = 0.L;
  }

  psd = corrin.psd.data[nby2];
  if (psd)
  {
     r1 = corrin.signal1.data[nby2];
     r2 = corrin.signal2.data[nby2];
     rShift = cos(LAL_PI*(double) params->lag - params->phase);  /* cos(2 \pi f_k \tau) */
     output->data[nby2] = (float)(rShift*r1*r2 / psd);
  }
  else
  {
     output->data[nby2] = 0.L;
  }

  /*
   * printf("&\n");
   */

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
