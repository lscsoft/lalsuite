/* <lalVerbatim file="LALCoarseTestCV">
Author: Churches, D. K. and Sathyaprakash, B. S.
$Id$
</lalVerbatim> */

/* <lalLaTeX>
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
lalDebugLevel
InspiralWaveLength
InspiralWave
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALCoarseTestCV}}
</lalLaTeX> */

/* <lalErrTable file="LALCoarseTestCE"> */
/* </lalErrTable> */

#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALNoiseModels.h>

INT4 lalDebugLevel=0;

/* Nlist = expected number of coarse grid templates; if not sufficient increase */
     
int 
main ( void )
{
   InspiralTemplateList *list=NULL;
   InspiralTemplateList *list2=NULL;
   static LALStatus status;
   static InspiralCoarseBankIn coarseIn;
   static InspiralFineBankIn   fineIn;
   static INT4 i, j, nlist, flist;

   coarseIn.mMin = 1.0;
   coarseIn.MMax = 40.0;
   coarseIn.mmCoarse = 0.80;
   coarseIn.mmFine = 0.97;
   coarseIn.fLower = 40.;
   coarseIn.fUpper = 2000;
   coarseIn.iflso = 0;
   coarseIn.tSampling = 4096.;
   coarseIn.NoisePsd = LALLIGOIPsd;
   coarseIn.order = twoPN;
   coarseIn.space = Tau0Tau2;
   coarseIn.method = one;
   coarseIn.approximant = taylor;
   coarseIn.domain = TimeDomain;
/* minimum value of eta */
   coarseIn.etamin = coarseIn.mMin * ( coarseIn.MMax - coarseIn.mMin) /
      pow(coarseIn.MMax,2.);

   LALInspiralCreateCoarseBank(&status, &list, &nlist, coarseIn);

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
     LALInspiralCreateFineBank(&status, &list2, &flist, fineIn);
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
  if (list!=NULL) LALFree(list);
  if (list2!=NULL) LALFree(list2);
  LALCheckMemoryLeaks();

  return(0);
}
