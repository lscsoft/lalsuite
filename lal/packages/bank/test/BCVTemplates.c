/********************************* <lalVerbatim file="BCVTemplatesCV">
Author: B.S. Sathyaprakash
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Program \texttt{BCVTemplates.c}}
\label{ss:BCVTemplates.c}

Creates a template mesh for BCV (or, alternatively, for SPA but
assuing a constant metric) using the mismatch metric.

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

\vfill{\footnotesize\input{BCVTemplatesCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <stdlib.h>
#include <lal/LALInspiralBank.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>

NRCSID(FLATMESHTESTC,"$Id$");

/* Default parameter settings. */
int lalDebugLevel = 0;

void  LALInspiralCreateFlatBank(LALStatus *status, REAL4VectorSequence *list, InspiralBankParams *bankParams);
static void
GetInspiralMoments (
		LALStatus            *status,
		InspiralMomentsEtc   *moments,
		REAL8FrequencySeries *psd,
		InspiralTemplate     *params );

int
main(int argc, char **argv)
{
  INT4 arg;
  static LALStatus status;     /* top-level status structure */
  static InspiralTemplate params;
  UINT4   numPSDpts=262144;
  INT4 nlist;
  REAL8FrequencySeries shf;
  REAL8 samplingRate, dx0, dx1;
  void *noisemodel = LALLIGOIPsd;
  static InspiralCoarseBankIn coarseIn;
  REAL4VectorSequence *list=NULL;
  static InspiralMetric metric;
  static InspiralMomentsEtc moments;
  static InspiralBankParams   bankParams; 
  static CreateVectorSequenceIn in; 

/* Number of templates is nlist */

  nlist = 0;

  params.OmegaS = 0.;
  params.Theta = 0.;
  params.ieta=1; 
  params.mass1=10.; 
  params.mass2=10.; 
  params.startTime=0.0; 
  params.startPhase=0.0;
  params.fLower=40.0; 
  params.fCutoff=2000.00;
  params.tSampling=4096.0;
  params.order=4;
  params.approximant=TaylorT3;
  params.signalAmplitude=1.0;
  params.nStartPad=0;
  params.nEndPad=1000;
  params.massChoice=m1Andm2;
  params.distance = 1.e8 * LAL_PC_SI/LAL_C_SI;
  LALInspiralParameterCalc(&status, &params);
    
  /* InspiralCoarseBankIn is not necessary except to get the ambiguity rectangles */
  coarseIn.fLower = params.fLower;
  coarseIn.fUpper = params.fCutoff;
  coarseIn.tSampling = params.tSampling;
  coarseIn.order = params.order;
  coarseIn.space = Tau0Tau3;
  coarseIn.approximant = params.approximant;
  coarseIn.mmCoarse = 0.95;
  coarseIn.mmFine = 0.97;
  coarseIn.iflso = 0.0L;
  coarseIn.mMin = 1.0;
  coarseIn.mMax = 20.0;
  coarseIn.MMax = coarseIn.mMax * 2.;
  coarseIn.massRange = MinMaxComponentMass; 
  /* coarseIn.massRange = MinComponentMassMaxTotalMass;*/
  /* minimum value of eta */
  coarseIn.etamin = coarseIn.mMin * ( coarseIn.MMax - coarseIn.mMin) / pow(coarseIn.MMax,2.);

  params.psi0 = 132250.;
  params.psi3 = -1314.2;
  params.alpha = 0.L;
  params.fendBCV = 868.7;

  samplingRate = params.tSampling;
  memset( &(shf), 0, sizeof(REAL8FrequencySeries) );
  shf.f0 = 0;
  LALDCreateVector( &status, &(shf.data), numPSDpts );
  shf.deltaF = samplingRate / (2.*(REAL8) shf.data->length + 1.L);
  LALNoiseSpectralDensity (&status, shf.data, noisemodel, shf.deltaF );

  /* compute the metric at this point, update bankPars and add the params to the list */
	  
  LALInspiralComputeMetricBCV(&status, &metric, &shf, &params);
  dx0 = sqrt(2.L * (1.L-coarseIn.mmCoarse)/metric.g00);
  dx1 = sqrt(2.L * (1.L-coarseIn.mmCoarse)/metric.g11);

  fprintf(stderr, "%e %e %e\n", metric.G00, metric.G01, metric.G11);
  fprintf(stderr, "%e %e %e\n", metric.g00, metric.g11, metric.theta);
  fprintf(stderr, "dp0=%e dp1=%e\n", dx0, dx1);

  bankParams.metric = &metric;
  bankParams.minimalMatch = coarseIn.mmCoarse;
  bankParams.x0Min = 1.0;
  bankParams.x0Max = 2.5e4;
  bankParams.x1Min = -2.2e3;
  bankParams.x1Max = 8.e2;

  in.length = 1;
  in.vectorLength = 2;
  LALSCreateVectorSequence(&status, &list, &in);

  list->vectorLength = 2;
  LALInspiralCreateFlatBank(&status, list, &bankParams);

  nlist = list->length;
  /* Prepare to print result. */
  {
    UINT4 j, numFcutTemplates=4;
    static InspiralTemplateList *tempParams;
    /* Print out the template parameters */
	    
    tempParams = (InspiralTemplateList *) LALMalloc (sizeof (InspiralTemplateList) * nlist);
    for (j=0; j<nlist; j++)
    {
	/*
	Retain only those templates that have meaningful chirptimes:
	*/
	    tempParams[j].params.psi0 = (REAL8) list->data[2*j];
	    tempParams[j].params.psi3 = (REAL8) list->data[2*j+1];
	    tempParams[j].params.fLower = params.fLower;
	    tempParams[j].params.nStartPad = params.nStartPad;
	    tempParams[j].params.nEndPad = params.nEndPad;
	    tempParams[j].params.tSampling= params.tSampling;
	    tempParams[j].params.distance =  params.distance;
	    tempParams[j].params.signalAmplitude= params.signalAmplitude;
	    printf("%e %e\n", tempParams[j].params.psi0, tempParams[j].params.psi3);
    }

    printf("#nlist before=%d\n&\n", nlist);
    LALInspiralBCVFcutBank( &status, &tempParams, &nlist, numFcutTemplates) ;
    printf("#nlist after=%d\n", nlist);
    for (j=0; j<nlist; j++)
    printf("%e %e %e %e\n", tempParams[j].params.psi0, tempParams[j].params.psi3, tempParams[j].params.totalMass, tempParams[j].params.fendBCV);
  }

  exit(0);
  {
    UINT4 j;
    UINT4 valid;
  
    static RectangleIn RectIn;
    static RectangleOut RectOut;

  
    RectIn.dx = sqrt(2.0 * (1. - coarseIn.mmCoarse)/metric.g00 );
    RectIn.dy = sqrt(2.0 * (1. - coarseIn.mmCoarse)/metric.g11 );
    RectIn.theta = metric.theta;
    
    params.massChoice=t03;
    /* Print out the template parameters */
    for (j=0; j<nlist; j++)
    {
	/*
	Retain only those templates that have meaningful masses:
	*/
	RectIn.x0 = bankParams.x0 = (REAL8) list->data[2*j];
	RectIn.y0 = bankParams.x1 = (REAL8) list->data[2*j+1];
	LALInspiralValidParams(&status, &valid, bankParams, coarseIn);
	valid = 1;
        if (valid) 
	{
		LALRectangleVertices(&status, &RectOut, &RectIn);
		printf("%e %e\n%e %e\n%e %e\n%e %e\n%e %e\n", 
				RectOut.x1, RectOut.y1, 
				RectOut.x2, RectOut.y2, 
				RectOut.x3, RectOut.y3, 
				RectOut.x4, RectOut.y4, 
				RectOut.x5, RectOut.y5);
		printf("&\n");
	}
    }
  }
  /* Free the list, and exit. */
  /*
  LALFree (list->data);
  LALFree (list);
  LALCheckMemoryLeaks();
  */
}


static void
GetInspiralMoments (
		LALStatus            *status,
		InspiralMomentsEtc   *moments,
		REAL8FrequencySeries *psd,
		InspiralTemplate     *params )
{

   UINT4 k;
   InspiralMomentsIn in;

   INITSTATUS (status, "GetInspiralMoments", FLATMESHTESTC);
   ATTATCHSTATUSPTR(status);
  
   ASSERT (params, status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   ASSERT (params->fLower>0, status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   ASSERT (moments, status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   ASSERT (psd, status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);

   moments->a01 = 3.L/5.L;
   moments->a21 = 11.L * LAL_PI/12.L;
   moments->a22 = 743.L/2016.L * pow(25.L/(2.L*LAL_PI*LAL_PI), 1.L/3.L);
   moments->a31 = -3.L/2.L;
   moments->a41 = 617.L * LAL_PI * LAL_PI / 384.L;
   moments->a42 = 5429.L/5376.L * pow ( 25.L * LAL_PI/2.L, 1.L/3.L);
   moments->a43 = 1.5293365L/1.0838016L * pow(5.L/(4.L*pow(LAL_PI,4.L)), 1.L/3.L);
   
   /* setup the input structure needed in the computation of the moments */

   in.shf = psd;
   in.shf->f0 /= params->fLower;
   in.shf->deltaF /= params->fLower;
   in.xmin = params->fLower/params->fLower;
   in.xmax = params->fCutoff/params->fLower;
	   
   /* First compute the norm */

   in.norm = 1.L;
   in.ndx = 7.L/3.L; 
   LALInspiralMoments(status->statusPtr, &moments->j[7], in); 
   CHECKSTATUSPTR(status);
   in.norm = moments->j[7];

   if (lalDebugLevel & LALINFO)
   {
	   fprintf (stderr, "a01=%e a21=%e a22=%e a31=%e a41=%e a42=%e a43=%e \n", 
			   moments->a01, moments->a21, moments->a22, moments->a31, 
			   moments->a41, moments->a42, moments->a43);
   
	   fprintf(stderr, "j7=%e\n", moments->j[7]);
   }

   /* Normalised moments of the noise PSD from 1/3 to 17/3. */

   for (k=1; k<=17; k++)
   {
	   in.ndx = (REAL8) k /3.L; 
	   LALInspiralMoments(status->statusPtr,&moments->j[k],in);  
	   CHECKSTATUSPTR(status);
	   if (lalDebugLevel==1) fprintf(stderr, "j%1i=%e\n", k,moments->j[k]);
   }
   in.shf->deltaF *= params->fLower;
   in.shf->f0 *= params->fLower;
  
   DETATCHSTATUSPTR(status);
   RETURN (status);
}
