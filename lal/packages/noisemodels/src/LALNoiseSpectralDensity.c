/*  <lalVerbatim file="LALNoiseSpectralDensityCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */
/* <lalLaTeX>
\subsection{Module \texttt{LALNoiseSpectralDensity.c}}
Module to create NoiseSpectralDensity.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALNoiseSpectralDensityCP}
\index{\verb&LALNoiseSpectralDensity()&}

\subsubsection*{Description}
\subsubsection*{Algorithm}
\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALNoiseSpectralDensityCV}}
</lalLaTeX>  */
#include <lal/LALNoiseModels.h>

NRCSID (LALNOISESPECTRALDENSITYC, "$Id$");

/*  <lalVerbatim file="LALNoiseSpectralDensityCP"> */
void 
LALNoiseSpectralDensity (
   LALStatus   *status, 
   REAL8Vector *psd, 
   void        (*NoisePsd)(LALStatus *status, REAL8 *shf, REAL8 f),
   REAL8       df) 
{  /*  </lalVerbatim>  */

    REAL8 shf, fs, xs, s0, dx, x;
    int i;

   INITSTATUS(status, "LALNoiseSpectralDensity", LALNOISESPECTRALDENSITYC);
   ATTATCHSTATUSPTR(status);

   ASSERT (psd->data,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (NoisePsd, status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (df > 0., status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);

   if (NoisePsd == LALGEOPsd) {
           s0 = 1.e-46;
           fs = 40;
   } else if(NoisePsd == LALLIGOIPsd) {
           s0 = 9.0e-46;
           fs = 40;
   } else if(NoisePsd == LALTAMAPsd) {
           s0 = 75.e-46;
           fs = 75;
   } else if(NoisePsd == LALVIRGOPsd) {
           s0 = 3.24e-46;
           fs = 20;
   } else {
           s0 = 1.;
           fs = 1.;
   }
   dx = df;
   xs = fs;
   psd->data[0] = psd->data[1] = 0.;
   for (i=2; i<(int)psd->length; i++) {
      x = (i-1)*dx;
      if (x>xs ) {
        (*NoisePsd)(status->statusPtr, &shf, x);
        psd->data[i] = s0 * shf;
      } else
        psd->data[i] = 0.;
   }
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

