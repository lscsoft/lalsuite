/*
*  Copyright (C) 2007 Bernd Machenschalk, David Churches, Duncan Brown, Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/* <lalVerbatim file="CoarseTestCV">
Author: Churches, D. K. and Sathyaprakash, B. S., Cokelaer, T
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

This test code creates a template bank and stores it into CoarseTest.out . Then, it creates
a finer template bank around a sub-set of the original template bank and stores it in the
same output file.

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

#include <lal/LALRCSID.h>
NRCSID (COARSETESTC,"$Id$");

INT4 lalDebugLevel=0;

/* Nlist = expected number of coarse grid templates; if not sufficient increase */

int
main ( void )
{
   InspiralTemplateList *coarseList=NULL;
   InspiralTemplateList *fineList=NULL;
   static LALStatus status;
   InspiralCoarseBankIn *coarseIn=NULL;
   InspiralFineBankIn   *fineIn=NULL;
   INT4 i, j, clist=0, flist=0;
   UINT4 numPSDpts = 262144;
   void (*noisemodel)(LALStatus*,REAL8*,REAL8) = LALLIGOIPsd;
   FILE *fpr;

   /* output filename*/
   fpr = fopen("CoarseTest.out", "w");

   /* define the two template banks*/
   coarseIn = (InspiralCoarseBankIn *)LALMalloc(sizeof(InspiralCoarseBankIn));
   fineIn = (InspiralFineBankIn *)LALMalloc(sizeof(InspiralFineBankIn));

   coarseIn->gridSpacing = Hexagonal;
   /* fill the coarseBankIn structure */
   coarseIn->mMin = 1.0;
   coarseIn->mMax = 30.0;
   coarseIn->MMax = coarseIn->mMax * 2.0;
   coarseIn->massRange = MinComponentMassMaxTotalMass;
   /* coarseIn->massRange = MinMaxComponentMass; */
   coarseIn->mmCoarse = 0.95;
   coarseIn->mmFine = 0.97;
   coarseIn->fLower = 40.;
   coarseIn->fUpper = 2000;
   coarseIn->iflso = 0;
   coarseIn->tSampling = 4096.L;
   coarseIn->order = twoPN;
   coarseIn->approximant = TaylorT3;
   coarseIn->space = Tau0Tau3;

   /* minimum value of eta  */
   coarseIn->etamin = coarseIn->mMin * ( coarseIn->MMax - coarseIn->mMin) /
      pow(coarseIn->MMax,2.);

   /* fill the psd */
   memset( &(coarseIn->shf), 0, sizeof(REAL8FrequencySeries) );
   coarseIn->shf.f0 = 0;
   LALDCreateVector(&status, &(coarseIn->shf.data), numPSDpts );
   coarseIn->shf.deltaF = coarseIn->tSampling / ( 2.L * (REAL8) coarseIn->shf.data->length + 1.L);
   LALNoiseSpectralDensity (&status, coarseIn->shf.data, noisemodel, coarseIn->shf.deltaF );

   /* create and save the coarse bank */
   coarseIn->iflso = 0.;
   LALInspiralCreateCoarseBank(&status, &coarseList, &clist, *coarseIn);
   fprintf(fpr, "#clist=%d\n",clist);
   for (i=0; i<clist; i++)
   {
      fprintf(fpr, "%e %e %e %e %e %e %e\n",
         coarseList[i].params.t0,
         coarseList[i].params.t3,
         coarseList[i].params.t2,
         coarseList[i].params.totalMass,
         coarseList[i].params.eta,
         coarseList[i].params.mass1,
         coarseList[i].params.mass2);
   }

  fprintf(fpr, "&\n");

  /* Then creates the finer bank around some of the original template bank*/
  fineIn->coarseIn = *coarseIn;
  for (j=0; j<clist; j+=48)
  {
     fineIn->templateList = coarseList[j];
     LALInspiralCreateFineBank(&status, &fineList, &flist, *fineIn);
     fprintf(fpr, "#flist=%d\n",flist);
     for (i=0; i<flist; i++) {
        fprintf(fpr, "%e %e %e %e %e %e %e\n",
        fineList[i].params.t0,
        fineList[i].params.t3,
        fineList[i].params.t2,
        fineList[i].params.totalMass,
        fineList[i].params.eta,
        fineList[i].params.mass1,
        fineList[i].params.mass2);
     }
     if (fineList!=NULL) LALFree(fineList);
     fineList = NULL;
     flist = 0;
  }
  fclose(fpr);
  LALDDestroyVector(&status, &(coarseIn->shf.data) );
  if (fineList!=NULL) LALFree(fineList);
  if (coarseList!=NULL) LALFree(coarseList);
  if (coarseIn!=NULL) LALFree(coarseIn);
  if (fineIn!=NULL) LALFree(fineIn);
  LALCheckMemoryLeaks();

  return(0);
}
