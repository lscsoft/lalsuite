/* <lalVerbatim file="CoarseTest2CV">
Author: Churches, D. K. and Sathyaprakash, B. S.
$Id$
</lalVerbatim> */

/* <lalLaTeX>
\subsection{Program \texttt{CoarseTest2.c}}
\label{ss:CoarseTest2.c}

Test code for the \texttt{bank} modules.

\subsubsection*{Usage}
\begin{verbatim}
CoarseTest2
\end{verbatim}

\subsubsection*{Description}

This test code gives an example of how to generate a template bank and
generates vertices of the ambiguity 'rectangle' around each lattice point
suitable for plotting with xmgr or xgrace.

\subsubsection*{Exit codes}
\input{CoarseTest2CE}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALRectangleVertices
LALInspiralCreateCoarseBank
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{CoarseTest2CV}}
</lalLaTeX> */

/* <lalErrTable file="CoarseTest2CE"> */
/* </lalErrTable> */

#include <stdio.h>
#include <lal/AVFactories.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALNoiseModels.h>

INT4 lalDebugLevel=0;

int 
main ( void )
{
   InspiralTemplateList *coarseList=NULL;
   static LALStatus status;
   static InspiralCoarseBankIn coarseIn;
   static InspiralFineBankIn   fineIn;
   static INT4 i, j, nlist;
   static RectangleIn RectIn;
   static RectangleOut RectOut;
   UINT4 numPSDpts = 131073;
   void *noisemodel = LALLIGOIPsd;
   static REAL8 f;

   coarseIn.mmCoarse = 0.95;
   coarseIn.mmFine = 0.99;
   coarseIn.fLower = 40.L;
   coarseIn.fUpper = 2000.L;
   coarseIn.iflso = 0.0L;
   coarseIn.tSampling = 4096.L;
   coarseIn.order = twoPN;
   coarseIn.space = Tau0Tau3;
   coarseIn.approximant = TaylorT1;

   coarseIn.mMin = 1.0;
   coarseIn.mMax = 20.0;
   coarseIn.MMax = coarseIn.mMax * 2.;

   coarseIn.massRange = MinMaxComponentMass; 
   /* coarseIn.massRange = MinComponentMassMaxTotalMass;*/

/* minimum value of eta */
   coarseIn.etamin = coarseIn.mMin * ( coarseIn.MMax - coarseIn.mMin) / pow(coarseIn.MMax,2.);
   /* fill the psd */
   memset( &(coarseIn.shf), 0, sizeof(REAL8FrequencySeries) );
   LALDCreateVector( &status, &(coarseIn.shf.data), numPSDpts );
   coarseIn.shf.f0 = 0.;
   coarseIn.shf.deltaF = coarseIn.tSampling / (REAL8) (2*(coarseIn.shf.data->length-1));

   /*
   for(i=1; i<=numPSDpts; i++)
   {
	   scanf("%le %le", &f, &(coarseIn.shf.data->data[i]));
	   printf("%e %e\n", f, coarseIn.shf.data->data[i]);
   }
   */
   LALNoiseSpectralDensity (&status, coarseIn.shf.data, noisemodel, coarseIn.shf.deltaF );

   LALInspiralCreateCoarseBank(&status, &coarseList, &nlist, coarseIn);

   fprintf(stderr, "nlist=%d\n",nlist);
   for (i=0; i<nlist; i++) 
   {
	 printf("%e %e %e %e %e %e %e\n", 
         coarseList[i].params.t0, 
         coarseList[i].params.t3, 
         coarseList[i].params.t2, 
         coarseList[i].params.mass1, 
         coarseList[i].params.mass2,
         coarseList[i].params.totalMass,
         coarseList[i].params.eta
	 );
   }


  printf("&\n");
  fineIn.coarseIn = coarseIn;
  for (j=0; j<nlist; j++) {
     
     RectIn.dx = sqrt(2. * (1. - coarseIn.mmCoarse)/coarseList[j].metric.g00 );
     RectIn.dy = sqrt(2. * (1. - coarseIn.mmCoarse)/coarseList[j].metric.g11 );
     RectIn.theta = fabs(coarseList[j].metric.theta);
     RectIn.x0 = coarseList[j].params.t0;

     switch (coarseIn.space) 
     {
        case Tau0Tau2: 
           RectIn.y0 = coarseList[j].params.t2;
        break;
        case Tau0Tau3: 
           RectIn.y0 = coarseList[j].params.t3;
        break;
        default: 
           exit(0);
     }

     LALRectangleVertices(&status, &RectOut, &RectIn);
     printf("%e %e\n%e %e\n%e %e\n%e %e\n%e %e\n", 
        RectOut.x1, RectOut.y1, 
        RectOut.x2, RectOut.y2, 
        RectOut.x3, RectOut.y3, 
        RectOut.x4, RectOut.y4, 
        RectOut.x5, RectOut.y5);
  
     /*
     printf("&\n");
     */
  }
  if (coarseList!=NULL) LALFree(coarseList);

  LALDDestroyVector( &status, &(coarseIn.shf.data) );
  LALCheckMemoryLeaks();
  return(0);
}
