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
   InspiralTemplateList *coarseList=NULL;
   InspiralTemplateList *fineList=NULL;
   LALStatus* status = LALCalloc(1, sizeof(*status));
   InspiralCoarseBankIn *coarseIn=NULL;
   InspiralFineBankIn   *fineIn=NULL;
   INT4 i, j, clist, flist;

   coarseIn = (InspiralCoarseBankIn *)LALMalloc(sizeof(InspiralCoarseBankIn));
   fineIn = (InspiralFineBankIn *)LALMalloc(sizeof(InspiralFineBankIn));

   coarseIn->mMin = 1.0;
   coarseIn->MMax = 40.0;
   coarseIn->mmCoarse = 0.80;
   coarseIn->mmFine = 0.97;
   coarseIn->fLower = 40.;
   coarseIn->fUpper = 1024L;
   coarseIn->iflso = 0;
   coarseIn->tSampling = 4000L;
   coarseIn->NoisePsd = LALLIGOIPsd;
   coarseIn->order = twoPN;
   coarseIn->approximant = TaylorT1;
   coarseIn->space = Tau0Tau3;
/* minimum value of eta */
   coarseIn->etamin = coarseIn->mMin * ( coarseIn->MMax - coarseIn->mMin) /
      pow(coarseIn->MMax,2.);

   coarseIn->iflso = 0.;
   LALInspiralCreateCoarseBank(status, &coarseList, &clist, *coarseIn);

   fprintf(stderr, "clist=%d\n",clist);
   for (i=0; i<clist; i++) 
   {
      printf("%e %e %e %e %e %e %e\n", 
         coarseList[i].params.t0, 
         coarseList[i].params.t3, 
         coarseList[i].params.t2, 
         coarseList[i].params.totalMass,
         coarseList[i].params.eta, 
         coarseList[i].params.mass1, 
         coarseList[i].params.mass2);
   }

  printf("&\n");

  fineIn->coarseIn = *coarseIn;
  for (j=0; j<clist; j+=48) 
  {
     fineIn->templateList = coarseList[j];
     LALInspiralCreateFineBank(status, &fineList, &flist, *fineIn);
     fprintf(stderr, "flist=%d\n",flist);
     for (i=0; i<flist; i++) {
        printf("%e %e %e %e %e %e %e\n", 
        fineList[i].params.t0, 
        fineList[i].params.t3, 
        fineList[i].params.t2, 
        fineList[i].params.totalMass,
        fineList[i].params.eta, 
        fineList[i].params.mass1, 
        fineList[i].params.mass1); 
     }
     if (fineList!=NULL) LALFree(fineList);
     fineList = NULL;
     flist = 0;
  }
  if (fineList!=NULL) LALFree(fineList);
  if (coarseList!=NULL) LALFree(coarseList);
  if (coarseIn!=NULL) LALFree(coarseIn);
  if (fineIn!=NULL) LALFree(fineIn);
  LALFree(status);
  LALCheckMemoryLeaks();

  return(0);
}
