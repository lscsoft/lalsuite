/*  <lalVerbatim file="NoisePSDTestCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */
/* <lalLaTeX>
\subsection{Module \texttt{NoisePSDTest.c}}
Module to create NoiseSpectralDensity.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{NoisePSDTestCP}
\idx{NoisePSDTest()}

\subsubsection*{Description}
\subsubsection*{Algorithm}
\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{NoisePSDTestCV}}
</lalLaTeX>  */
#include <lal/LALNoiseModels.h>

INT4 lalDebugLevel=1;
/*  <lalVerbatim file="NoisePSDTestCP"> */
int main () 
{  /*  </lalVerbatim>  */

   static LALStatus status;
   REAL8Vector psd;
   REAL8       df;
   INT4 i;

   df = 1.0;
   psd.length = 8193;

   psd.data = LALMalloc(sizeof(REAL8) * psd.length);

   LALNoiseSpectralDensity(&status, &psd, &LALGEOPsd, df); 

   for (i=2; i<(INT4)psd.length; i++) {
        if (psd.data[i]) printf ("%d %e\n", i, psd.data[i]);
   }
   printf("&\n");
   LALNoiseSpectralDensity(&status, &psd, &LALLIGOIPsd, df); 
   for (i=2; i<(INT4)psd.length; i++) {
        if (psd.data[i]) printf ("%d %e\n", i, psd.data[i]);
   }
   printf("&\n");
   LALNoiseSpectralDensity(&status, &psd, &LALVIRGOPsd, df); 
   for (i=2; i<(INT4)psd.length; i++) {
        if (psd.data[i]) printf ("%d %e\n", i, psd.data[i]);
   }
   printf("&\n");
   LALNoiseSpectralDensity(&status, &psd, &LALTAMAPsd, df); 
   for (i=2; i<(INT4)psd.length; i++) {
        if (psd.data[i]) printf ("%d %e\n", i, psd.data[i]);
   }
  return 0;
}

