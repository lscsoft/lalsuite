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
   FILE *NoisePsdFile;

   fprintf(stderr, "This test code computes the amplitude spectrum of 
    GEO, LIGO, TAMA and VIRGO and writes them in NoisePSDTest.out in a
    format suitable for display with xmgr/xgrace\n");

   if ( (NoisePsdFile = fopen("NoisePSDTest.out", "w")) == NULL) 
   {
      fprintf(stderr, "Can't open output file\n");
      exit(0);
   }
   df = 1.0;
   psd.length = 8193;

   psd.data = LALMalloc(sizeof(REAL8) * psd.length);

   LALNoiseSpectralDensity(&status, &psd, &LALGEOPsd, df); 

   for (i=2; i<(INT4)psd.length; i++) {
        if (psd.data[i]) fprintf (NoisePsdFile, "%d %e\n", i, sqrt(psd.data[i]));
   }
   fprintf(NoisePsdFile, "&\n");
   LALNoiseSpectralDensity(&status, &psd, &LALLIGOIPsd, df); 
   for (i=2; i<(INT4)psd.length; i++) {
        if (psd.data[i]) fprintf (NoisePsdFile, "%d %e\n", i, sqrt(psd.data[i]));
   }
   fprintf(NoisePsdFile, "&\n");
   LALNoiseSpectralDensity(&status, &psd, &LALVIRGOPsd, df); 
   for (i=2; i<(INT4)psd.length; i++) {
        if (psd.data[i]) fprintf (NoisePsdFile, "%d %e\n", i, sqrt(psd.data[i]));
   }
   fprintf(NoisePsdFile, "&\n");
   LALNoiseSpectralDensity(&status, &psd, &LALTAMAPsd, df); 
   for (i=2; i<(INT4)psd.length; i++) {
        if (psd.data[i]) fprintf (NoisePsdFile, "%d %e\n", i, sqrt(psd.data[i]));
   }
   LALFree(psd.data);
   LALCheckMemoryLeaks();
  return 0;
}

