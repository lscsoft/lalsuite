/* <lalVerbatim file="CoarseTestCV">
Author: Churches, D. K. and Sathyaprakash, B. S.
$Id$
</lalVerbatim> */

/* <lalLaTeX>
\subsection{Program \texttt{CoarseTest.c}}
\label{ss:CoarseTest.c}

Test code for the inspiral modules.

\subsubsection*{Usage}
\begin{verbatim}
CoarseTest
\end{verbatim}

\subsubsection*{Description}

This test code gives an example of how one calls \texttt{LALInspiralCreateCoarseBank}
and \texttt{LALInspiralCreateFineBank} modules.

\subsubsection*{Exit codes}
\input{CoarseTestCE}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALInspiralCreateCoarseBank
LALInspiralCreateFineBank
LALFree
LALCheckMemoryLeaks
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{CoarseTestCV}}
</lalLaTeX> */

/* <lalErrTable file="CoarseTestCE"> */
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
   UINT4 numPSDpts = 1048576;
   void *noisemodel = LALLIGOIPsd;

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
   coarseIn->order = twoPN;
   coarseIn->approximant = TaylorT1;
   coarseIn->space = Tau0Tau3;
/* minimum value of eta */
   coarseIn->etamin = coarseIn->mMin * ( coarseIn->MMax - coarseIn->mMin) /
      pow(coarseIn->MMax,2.);
   /* fill the psd */
   memset( &(coarseIn->shf), 0, sizeof(REAL8FrequencySeries) );
   coarseIn->shf.f0 = 0;
   LALDCreateVector( status, &(coarseIn->shf.data), numPSDpts );
   coarseIn->shf.deltaF = coarseIn->tSampling / (REAL8) coarseIn->shf.data->length;
   LALNoiseSpectralDensity (status, coarseIn->shf.data, noisemodel, coarseIn->shf.deltaF );
   

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
  LALDDestroyVector( status, &(coarseIn->shf.data) );
  if (coarseList!=NULL) LALFree(coarseList);
  if (coarseIn!=NULL) LALFree(coarseIn);
  if (fineIn!=NULL) LALFree(fineIn);
  LALFree(status);
  LALCheckMemoryLeaks();

  return(0);
}
