/********************************* <lalVerbatim file="BCVSpinTemplatesCV">
Author: B.S. Sathyaprakash
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Program \texttt{BCVSpinTemplates.c}}
\label{ss:BCVSpinTemplates.c}

Creates a template mesh for BCVSpin using the mismatch metric.

\subsubsection*{Usage}

\subsubsection*{Description}

\subsubsection*{Exit codes}
****************************************** </lalLaTeX><lalErrTable> */
/******************************************** </lalErrTable><lalLaTeX>

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{BCVSpinTemplatesCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <stdlib.h>
#include <lal/LALInspiralBank.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>

NRCSID(FLATMESHTESTC,"$Id$");

/* Default parameter settings. */
INT4 lalDebugLevel = 0;

int
main(int argc, char **argv)
{
  /* top-level status structure */
  UINT4   numPSDpts=65537;
  static LALStatus status;     
  static InspiralCoarseBankIn coarseIn;
  void *noisemodel = LALLIGOIPsd;
  double beta;
  SnglInspiralTable *tiles=NULL;
  FILE *fpr;

/* Number of templates is nlist */
  INT4 nlist1, j;


  fpr = fopen("BCVSpinTemplates.out", "w");
  nlist1 = 0;
  coarseIn.HighGM = 6.; 
  coarseIn.LowGM = -4.;
  coarseIn.fLower = 40.L;
  coarseIn.fUpper = 400.L;
  coarseIn.tSampling = 4096.L;
  coarseIn.order = twoPN;
  coarseIn.mmCoarse = 0.90;
  coarseIn.mmFine = 0.97;
  coarseIn.iflso = 0.0L;
  coarseIn.mMin = 3.0;
  coarseIn.mMax = 20.0;
  coarseIn.MMax = coarseIn.mMax * 2.;
  coarseIn.massRange = MinMaxComponentMass; 
  /* coarseIn.massRange = MinComponentMassMaxTotalMass;*/
  /* minimum value of eta */
  coarseIn.etamin = coarseIn.mMin * ( coarseIn.MMax - coarseIn.mMin) / pow(coarseIn.MMax,2.);
  coarseIn.psi0Min = 1.0e4;
  coarseIn.psi0Max = 6.0e5;
  coarseIn.psi3Min = -3.0e3;
  coarseIn.psi3Max = -1.0e1;
  coarseIn.alpha = 0.L;
  coarseIn.numFcutTemplates = 4;
  coarseIn.betaMin = 0.0;
  coarseIn.betaMax = 800.;
  coarseIn.spinBank = 1;
  coarseIn.gridSpacing = Hexagonal;

  memset( &(coarseIn.shf), 0, sizeof(REAL8FrequencySeries) );
  coarseIn.shf.f0 = 0;
  LALDCreateVector( &status, &(coarseIn.shf.data), numPSDpts );
  coarseIn.shf.deltaF = coarseIn.tSampling / (2.*(REAL8) coarseIn.shf.data->length + 1.L);
  LALNoiseSpectralDensity (&status, coarseIn.shf.data, noisemodel, coarseIn.shf.deltaF );

  coarseIn.approximant = BCVSpin;
  coarseIn.space       = Psi0Psi3;
  
  // LALInspiralBCVSpinBank (&status, &tiles, &nlist1, &coarseIn);
  LALInspiralBankGeneration(&status, &coarseIn, &tiles, &nlist1);
  fprintf (fpr, "#numtemplaes=%d\n", nlist1);
  beta = tiles->beta;
  for (j=0; j<nlist1; j++)
  {
	  fprintf(fpr, "%e %e %e\n", tiles->psi0, tiles->psi3, beta);
	  tiles = tiles->next;
	  if (tiles != NULL && beta != tiles->beta)
	  {
		  beta = tiles->beta;
		  fprintf(fpr, "&\n");
	  }
  }

  fclose(fpr);
  LALDDestroyVector( &status, &(coarseIn.shf.data) );
  LALFree(tiles);
  LALCheckMemoryLeaks();
  return 0;
}
