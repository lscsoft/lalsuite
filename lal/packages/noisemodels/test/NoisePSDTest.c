/*  <lalVerbatim file="NoisePSDTestCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */
/* <lalLaTeX>
\subsection{Program \texttt{NoisePSDTest.c}}
\label{ss:NoisePSDTest.c}

This program can be used generate expected noise
NoiseSpectralDensity in various interferometers.
See the beginning of the NoiseModels module to see details on how
this test program works.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{NoisePSDTestCP}
\idx{NoisePSDTest()}

\subsubsection*{Description}
\subsubsection*{Algorithm}
\subsubsection*{Uses}
\begin{verbatim}
LALDCreateVector
LALNoiseSpectralDensity
LALGEOPsd
LALLIGOIPsd
LALTAMAPsd
LALVIRGOPsd
LALAdvLIGOPsd
LALDDestroyVector
LALCheckMemoryLeaks
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{NoisePSDTestCV}}
</lalLaTeX>  */
#include <lal/AVFactories.h>
#include <lal/LALNoiseModels.h>

INT4 lalDebugLevel=1;
/*  <lalVerbatim file="NoisePSDTestCP"> */
int main ( void ) 
{  /*  </lalVerbatim>  */

   static LALStatus status;
   REAL8Vector *psd=NULL;
   REAL8       df;
   INT4 i, length;
   FILE *NoisePsdFile;

   fprintf(stderr, "This test code computes the amplitude spectrum \n");
   fprintf(stderr, "of GEO, LIGO, VIRGO, TAMA and AdvLIGO and writes them in \n");
   fprintf(stderr, "NoisePSDTest.out in a format suitable for \n");
   fprintf(stderr, "display with xmgr/xgrace\n");

   if ( (NoisePsdFile = fopen("NoisePSDTest.out", "w")) == NULL) 
   {
      fprintf(stderr, "Can't open output file\n");
      exit(0);
   }
   df = 1.0;
   length = 8193;

   LALDCreateVector(&status, &psd, length);
   fprintf(stderr, "Length of vector=%d\n", psd->length);

   LALNoiseSpectralDensity(&status, psd, &LALGEOPsd, df); 

   for (i=2; i<length; i++) {
        if (psd->data[i]) fprintf (NoisePsdFile, "%d %e\n", i, sqrt(psd->data[i]));
   }
   fprintf(NoisePsdFile, "&\n");
   LALNoiseSpectralDensity(&status, psd, &LALLIGOIPsd, df); 
   for (i=2; i<length; i++) {
        if (psd->data[i]) fprintf (NoisePsdFile, "%d %e\n", i, sqrt(psd->data[i]));
   }
   fprintf(NoisePsdFile, "&\n");
   LALNoiseSpectralDensity(&status, psd, &LALVIRGOPsd, df); 
   for (i=2; i<length; i++) {
        if (psd->data[i]) fprintf (NoisePsdFile, "%d %e\n", i, sqrt(psd->data[i]));
   }
   fprintf(NoisePsdFile, "&\n");
   LALNoiseSpectralDensity(&status, psd, &LALTAMAPsd, df); 
   for (i=2; i<length; i++) {
        if (psd->data[i]) fprintf (NoisePsdFile, "%d %e\n", i, sqrt(psd->data[i]));
   }
   fprintf(NoisePsdFile, "&\n");
   LALNoiseSpectralDensity(&status, psd, &LALAdvLIGOPsd, df); 
   for (i=2; i<length; i++) {
        if (psd->data[i]) fprintf (NoisePsdFile, "%d %e\n", i, sqrt(psd->data[i]));
   }
   LALDDestroyVector(&status, &psd);
   LALCheckMemoryLeaks();
   return 0;
}

