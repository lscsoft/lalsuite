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
\index{\texttt{LALNoiseSpectralDensity()}}

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
void 
/*  <lalVerbatim file="LALNoiseSpectralDensityCP"> */
LALNoiseSpectralDensity (
   LALStatus   *status, 
   REAL8Vector *psd, 
   Detector    choice, 
   REAL8       df) 
{  /*  </lalVerbatim>  */
    REAL8 f0, fs, xs, s0, dx, x, (*NoisePsd)(REAL8);
    int i;

   INITSTATUS (status, "LALNoiseSpectralDensity", LALNOISESPECTRALDENSITYC);
   switch (choice) {
	case 0:
	   NoisePsd = &LALGEOPsd;
           s0 = 1.e-46;
           f0 = 150;
           fs = 40;
	   break;
	case 1:
	   NoisePsd = &LALLIGOIPsd;
           s0 = 9.0e-46;
           f0 = 150;
           fs = 40;
	   break;
	case 2:
	   NoisePsd = &LALTAMAPsd;
           s0 = 75.e-46;
           f0 = 400;
           fs = 75;
	   break;
	case 3:
	   NoisePsd = &LALVIRGOPsd;
           s0 = 3.24e-46;
           f0 = 500;
           fs = 20;
	   break;
        default:
           fprintf(stderr, "LALNoiseSpectralDensity: Invalid Detector\n");
           exit(1);
           break;
   }
   dx = df/f0;
   xs = fs/f0;
   psd->data[0] = psd->data[1] = 0.;
/*
   s0 = 1.;
*/
   for (i=2; i<(int)psd->length; i++) {
      x = (i-1)*dx;
      if (x>xs ) {
        psd->data[i] = s0 * (*NoisePsd)(x);
      } else
        psd->data[i] = 0.;
   }
  RETURN(status);
}

