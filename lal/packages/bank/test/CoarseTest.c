/* <lalVerbatim file="CoarseTestCV">
Author: Churches, D. K. and Sathyaprakash, B. S.
$Id$
</lalVerbatim> */

/* <lalLaTeX>
\subsection{Program \texttt{CoarseTest.c}}
\label{ss:CoarseTest.c}

Test code for the inspiral modules. A template bank can be created either using
a full range for component masses of the binary $m_1$ and $m_2$, that is
\texttt{(mMin,mMax),} OR minimum value of the component masses \texttt{mMin} and 
maximum value of the total mass \texttt{MMax.} In the first case chirptimes 
satisfying the constraint \texttt{mMin}~$\le m_1, m_2 \le$~\texttt{mMax} are accepted
as valid systems. In the second case chirptimes satisfying the
constraint \texttt{mMin}~$\le m_1, m_2,$ and \texttt{MMax}$\le m=m_1+m_2,$
are treated as valid.  Users are expected to provide both \texttt{mMax} and \texttt{MMax}.

For \texttt{LALLIGOIPsd} the choice \texttt{mMin}$=1M_\odot$ \texttt{mMax}$=20 M_\odot$
gives 2292 templates, while the same \texttt{mMin} but choosing \texttt{MMax}$=40 M_\odot$
instead gives 2512 templates -- about 10\% incrase. However, the extra templates are ALL
short-lived templates and therefore potentially trigger a large number of false alarms,
as noted by Brown in E7 analysis.

\subsubsection*{Usage}
Input the following values of the InspiralCoarseBankIn structure to
create a template bank:

\begin{itemize}
   \item Minimum component mass in solar masses.
   \texttt{mMin = 1.0;}

   \item Maximum component mass in solar masses.
   \texttt{mMax = 20.0;}

   \item Maximum total mass. {\bf This should be specified independently of
   whether or not mMax is specified.} It is used in setting up a 
   rectangular area in the space of chirptimes where the templates will
   be laid.
   \texttt{MMax = 40.0;  }

   \item The type of search space. 
   \texttt{massRange = MinComponentMassMaxTotalMass;} OR
   \texttt{massRange = MinMaxComponentMasses;}

   \item Coarse bank minimal match
   \texttt{coarseIn->mmCoarse = 0.80;}

   \item Fine bank minimal match
   \texttt{mmFine = 0.97;}

   \item Lower frequency cutoff
   \texttt{fLower = 40.;}

   \item Upper frequency cutoff
   \texttt{fUpper = 1024L;}

   \item Whether or not lso should be used as an upper frequency cutoff
   (Currently not used; so please specify \texttt{fUpper}.
   \texttt{coarseIn->iflso = 0;}

   \item Sampling rate
   \texttt{tSampling = 4000L;}

   \item Space in which templates should be created: whether $(\tau_0,\tau_2)$
   or $(\tau_0, \tau_3).$
   \texttt{coarseIn->space = Tau0Tau2;} OR
   \texttt{coarseIn->space = Tau0Tau3;} OR

   \item Order and type of the approximant to be used in template generation.
   These members will NOT be used in creating a template bank but in 
   filling up the \texttt{InspiralTemplate} structure created by the bank. 
   
   \texttt{order = twoPN;}
   \texttt{coarseIn->approximant = TaylorT1;}

   \item minimum value of eta 
   \texttt{etamin = mMin * ( MMax - mMin)/pow(MMax,2.);}

   \item Finally, fill the psd structure (see test code for an example).
   This involves allocating memory to vector \texttt{shf.data} and filling it with 
   noise PSD as also choosing the starting frequency and the frequency resolution. 
   
\end{itemize}

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
#include <lal/AVFactories.h>

INT4 lalDebugLevel=1;

/* Nlist = expected number of coarse grid templates; if not sufficient increase */
     
int 
main ( void )
{
   InspiralTemplateList *coarseList=NULL;
   InspiralTemplateList *fineList=NULL;
   LALStatus* status = LALCalloc(1, sizeof(*status));
   InspiralCoarseBankIn *coarseIn=NULL;
   InspiralFineBankIn   *fineIn=NULL;
   INT4 i, j, clist=0, flist=0;
   UINT4 numPSDpts = 262144;
   void *noisemodel = LALLIGOIPsd;

   coarseIn = (InspiralCoarseBankIn *)LALMalloc(sizeof(InspiralCoarseBankIn));
   fineIn = (InspiralFineBankIn *)LALMalloc(sizeof(InspiralFineBankIn));

   coarseIn->mMin = 5.0;
   coarseIn->mMax = 20.0;
   coarseIn->MMax = coarseIn->mMax * 2.0;
   coarseIn->massRange = MinComponentMassMaxTotalMass;
   /* coarseIn->massRange = MinMaxComponentMass; */

   coarseIn->mmCoarse = 0.97;
   coarseIn->mmFine = 0.97;
   coarseIn->fLower = 40.;
   coarseIn->fUpper = 400;
   coarseIn->iflso = 0;
   coarseIn->tSampling = 2048.L;
   coarseIn->order = twoPN;
   coarseIn->approximant = TaylorT1;
   coarseIn->space = Tau0Tau3;

   /* minimum value of eta  */
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
