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

INT4 lalDebugLevel=1;

int
main(int argc, char **argv)
{
  INT4 arg;
  /* top-level status structure */
  static LALStatus status;     
  /* Structure specifying the nature of the bank needed */
  static InspiralCoarseBankIn coarseIn;
/* Template bank lists */
  static InspiralTemplateList *list1, *list2;
/* Number of templates in list1 and list2 */
  INT4 nlist1=0, nlist2=0;

  void *noisemodel = LALLIGOIPsd;
  UINT4   j, numPSDpts=262144;
  FILE *fpr;


  fpr = fopen("CoarseTest2.out", "w");
  coarseIn.LowGM = 3;
  coarseIn.HighGM= 6;
  coarseIn.fLower = 40.L;
  coarseIn.fUpper = 2000.L;
  coarseIn.tSampling = 4096.L;
  coarseIn.order = twoPN;
  coarseIn.space = Tau0Tau3;
  coarseIn.mmCoarse = 0.95;
  coarseIn.mmFine = 0.97;
  coarseIn.iflso = 0.0L;
  coarseIn.mMin = 3.0;
  coarseIn.mMax = 20.0;
  coarseIn.MMax = coarseIn.mMax * 2.;
  coarseIn.massRange = MinMaxComponentMass; 
  /* coarseIn.massRange = MinComponentMassMaxTotalMass;*/
  /* minimum value of eta */
  coarseIn.etamin = coarseIn.mMin * ( coarseIn.MMax - coarseIn.mMin) / pow(coarseIn.MMax,2.);
  coarseIn.psi0Min = 1.e0;
  coarseIn.psi0Max = 2.5e4;
  coarseIn.psi3Min = -2.2e3;
  coarseIn.psi3Max = 8.e2;
  coarseIn.alpha = 0.L;
  coarseIn.numFcutTemplates = 4;

  memset( &(coarseIn.shf), 0, sizeof(REAL8FrequencySeries) );
  coarseIn.shf.f0 = 0;
  LALDCreateVector( &status, &(coarseIn.shf.data), numPSDpts );
  coarseIn.shf.deltaF = coarseIn.tSampling / (2.*(REAL8) coarseIn.shf.data->length + 1.L);
  LALNoiseSpectralDensity (&status, coarseIn.shf.data, noisemodel, coarseIn.shf.deltaF );

  coarseIn.approximant = BCV;
  coarseIn.space 	= Psi0Psi3;
  LALInspiralCreateCoarseBank(&status, &list1, &nlist1, coarseIn);
  for (j=0; j<nlist1; j++)
  {
	  fprintf(fpr, "%e %e %e %e\n", 
			  list1[j].params.psi0, 
			  list1[j].params.psi3, 
			  list1[j].params.totalMass, 
			  list1[j].params.fFinal);
  }
		
  fprintf(fpr, "&\n");
  coarseIn.approximant = TaylorT1;
  coarseIn.space 	= Tau0Tau3;
  LALInspiralCreateCoarseBank(&status, &list2, &nlist2, coarseIn);
    
  for (j=0; j<nlist2; j++)
  {
	  fprintf(fpr, "%e %e %e %e\n", 
			  list2[j].params.t0, 
			  list2[j].params.t3,
			  list2[j].params.mass1, 
			  list2[j].params.mass2 
			  );
  }
  fprintf(fpr, "&\n");

  coarseIn.approximant = EOB;
  LALInspiralCreateCoarseBank(&status, &list2, &nlist2, coarseIn);
  
  {
    UINT4 j;
    UINT4 valid;
  
    static RectangleIn RectIn;
    static RectangleOut RectOut;

  
    RectIn.dx = sqrt(2.0 * (1. - coarseIn.mmCoarse)/list1[0].metric.g00 );
    RectIn.dy = sqrt(2.0 * (1. - coarseIn.mmCoarse)/list1[0].metric.g11 );
    RectIn.theta = list1[0].metric.theta;
    
    /* Print out the template parameters */
    for (j=0; j<nlist1; j++)
    {
	/*
	Retain only those templates that have meaningful masses:
	*/
	RectIn.x0 = (REAL8) list1[j].params.psi0;
	RectIn.y0 = (REAL8) list1[j].params.psi3;
	/*
	LALInspiralValidParams(&status, &valid, bankParams, coarseIn);
	*/
	valid = 1;
        if (valid) 
	{
		LALRectangleVertices(&status, &RectOut, &RectIn);
		fprintf(fpr, "%e %e\n%e %e\n%e %e\n%e %e\n%e %e\n", 
				RectOut.x1, RectOut.y1, 
				RectOut.x2, RectOut.y2, 
				RectOut.x3, RectOut.y3, 
				RectOut.x4, RectOut.y4, 
				RectOut.x5, RectOut.y5);
		fprintf(fpr, "&\n");
	}
    }
  }
  /* Free the list, and exit. */
  if (list1 != NULL) LALFree (list1);
  if (list2 != NULL) LALFree (list2);
  LALDDestroyVector( &status, &(coarseIn.shf.data) );
  LALCheckMemoryLeaks();
  return(0);
}
