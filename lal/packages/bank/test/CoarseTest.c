/******************************** <lalVerbatim file="LALCoarseTestCV">
Author: Churches, D. K. and Sathyaprakash, B. S.
$Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{LALCoarseTest.c}}
\label{ss:LALCoarseTest.c}

Test code for the inspiral modules.

\subsubsection*{Usage}
\begin{verbatim}
LALCoarseTest
\end{verbatim}

\subsubsection*{Description}

This test code gives an example of how one calls the functions which generate inspiral
waveforms.
Note that one must calculate the length of the waveform and allocate memory for it
\emph{before} calling
\texttt{InspiralWave}. The length of the waveform can be calculated by calling the function
\texttt{InspiralWaveLength} beforehand, as shown.

There are only two functions which one can call to generate waveforms. These are
\texttt{InspiralWave},
which will return a \emph{single} waveform, and \texttt{InspiralWaveTemplates}, which
returns a \emph{pair}
of waveforms which have phases which differ by $\pi/2$.

\subsubsection*{Exit codes}
\input{LALCoarseTestCE}

\subsubsection*{Uses}
\begin{verbatim}
laldebuglevel
InspiralWaveLength
InspiralWave
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALCoarseTestCV}}
******************************************************* </lalLaTeX> */


/***************************** <lalErrTable file="LALCoarseTestCE"> */
/***************************** </lalErrTable> */





#include <stdio.h>
#include <lal/LALInspiralBank.h>

INT4 lalDebugLevel=1;
/* Nlist = expected number of coarse grid templates; if not sufficient increase */
     
int 
main ( void )
{
   InspiralTemplateList *list;
   InspiralTemplateList *list2;
   static LALStatus status;
   InspiralCoarseBankIn coarseIn;
   InspiralFineBankIn   fineIn;
   INT4 Nlist, Nlist2, i, j, nlist, flist;


   Nlist = 50000;
/* Nlist2 = expected number of fine grid templates around a given coarse 
grid point */
   Nlist2 = 100;
   list = (InspiralTemplateList *) LALMalloc(sizeof(InspiralTemplateList) * Nlist);
   list2 = (InspiralTemplateList *) LALMalloc(sizeof(InspiralTemplateList) * Nlist2);
   coarseIn.mMin = 1.0;
   coarseIn.MMax = 40.0;
   coarseIn.mmCoarse = 0.80;
   coarseIn.mmFine = 0.97;
   coarseIn.fLower = 40.;
   coarseIn.fUpper = 2000;
   coarseIn.iflso = 0;
   coarseIn.tSampling = 4000.;
   coarseIn.detector = ligo;
   coarseIn.method = one;
   coarseIn.order = twoPN;
   coarseIn.approximant = taylor;
   coarseIn.domain = TimeDomain;
   coarseIn.space = Tau0Tau3;
/* minimum value of eta */
   coarseIn.etamin = coarseIn.mMin * ( coarseIn.MMax - coarseIn.mMin) /
      pow(coarseIn.MMax,2.);

   LALInspiralCreateCoarseBank(&status, list, &nlist, coarseIn);
   fprintf(stderr, "nlist=%d\n",nlist);
   for (i=0; i<nlist; i++) {
      printf("%e %e %e %e %e %e %e\n", 
         list[i].params.t0, 
         list[i].params.t2, 
         list[i].params.t3, 
         list[i].params.totalMass,
         list[i].params.eta, 
         list[i].params.mass1, 
         list[i].params.mass2);

   }

  printf("&\n");

  fineIn.coarseIn = coarseIn;
  for (j=0; j<nlist; j+=48) {
     fineIn.templateList = list[j];
     LALInspiralCreateFineBank(&status, list2, &flist, fineIn);
     fprintf(stderr, "flist=%d\n",flist);
      for (i=0; i<flist; i++) {
         printf("%e %e %e %e %e %e %e\n", 
         list2[i].params.t0, 
         list2[i].params.t2, 
         list2[i].params.t3, 
         list2[i].params.totalMass,
         list2[i].params.eta, 
         list2[i].params.mass1, 
         list2[i].params.mass1); 
      }
   }

   return(0);
}

