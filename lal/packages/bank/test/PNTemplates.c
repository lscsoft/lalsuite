/*
*  Copyright (C) 2007 Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
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

/********************************* <lalVerbatim file="PNTemplatesCV">
Author: B.S. Sathyaprakash
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Program \texttt{PNTemplates.c}}
\label{ss:PNTemplates.c}

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

\vfill{\footnotesize\input{PNTemplatesCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <stdlib.h>
#include <lal/LALInspiralBank.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>

NRCSID(PNTEMPLATESC,"$Id$");

/* Default parameter settings. */
int lalDebugLevel = 0;

/*static void PSItoMasses (LALStatus *status, InspiralTemplate *params, UINT4 *valid, REAL4 psi0, REAL4 psi3);*/
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
  static LALStatus status;     /* top-level status structure */
  static InspiralTemplate params;
  UINT4   numPSDpts=262144;
  INT4 nlist;
  REAL8FrequencySeries shf;
  REAL8 samplingRate, dx0, dx1;
  void (*noisemodel)(LALStatus*,REAL8*,REAL8) = LALLIGOIPsd;
  static InspiralCoarseBankIn coarseIn;
  REAL4VectorSequence *list=NULL;
  static InspiralMetric metric;
  static InspiralMomentsEtc moments;
  static InspiralBankParams   bankParams;
  static CreateVectorSequenceIn in;
  INT4 valid;
  FILE *fpr;
  fpr = fopen("PNTemplates.out", "w");

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

  coarseIn.fLower = params.fLower;
  coarseIn.fUpper = params.fCutoff;
  coarseIn.tSampling = params.tSampling;
  coarseIn.order = params.order;
  coarseIn.space = Tau0Tau3;
  coarseIn.approximant = params.approximant;
  coarseIn.mmCoarse = 0.70;
  coarseIn.mmFine = 0.97;
  coarseIn.iflso = 0.0L;
  coarseIn.mMin = 1.0;
  coarseIn.mMax = 20.0;
  coarseIn.MMax = coarseIn.mMax * 2.;
  coarseIn.massRange = MinMaxComponentMass;
  /* coarseIn.massRange = MinComponentMassMaxTotalMass;*/
  /* minimum value of eta */
  coarseIn.etamin = coarseIn.mMin * ( coarseIn.MMax - coarseIn.mMin) / pow(coarseIn.MMax,2.);

  metric.space = Tau0Tau3;

  samplingRate = params.tSampling;
  memset( &(shf), 0, sizeof(REAL8FrequencySeries) );
  shf.f0 = 0;
  LALDCreateVector( &status, &(shf.data), numPSDpts );
  shf.deltaF = samplingRate / (2.*(REAL8) shf.data->length + 1.L);
  LALNoiseSpectralDensity (&status, shf.data, noisemodel, shf.deltaF );

  /* compute the metric */

  GetInspiralMoments (&status, &moments, &shf, &params);
  LALInspiralComputeMetric(&status, &metric, &params, &moments);
  dx0 = sqrt(2.L * (1.L-coarseIn.mmCoarse)/metric.g00);
  dx1 = sqrt(2.L * (1.L-coarseIn.mmCoarse)/metric.g11);

  fprintf(fpr, "%e %e %e\n", metric.G00, metric.G01, metric.G11);
  fprintf(fpr, "%e %e %e\n", metric.g00, metric.g11, metric.theta);
  fprintf(fpr, "dp0=%e dp1=%e\n", dx0, dx1);

  bankParams.metric = &metric;
  bankParams.minimalMatch = coarseIn.mmCoarse;
  bankParams.x0Min = 3.0;
  bankParams.x0Max = 10.00;
  bankParams.x1Min = 0.15;
  bankParams.x1Max = 1.25;

  in.length = 1;
  in.vectorLength = 2;
  LALSCreateVectorSequence(&status, &list, &in);

  list->vectorLength = 2;
  LALInspiralCreateFlatBank(&status, list, &bankParams);
  nlist = list->length;

  fprintf(fpr, "Number of templates=%d dx0=%e dx1=%e\n", nlist, bankParams.dx0, bankParams.dx1);


  /* Prepare to print result. */
  {
    INT4 j;
    /* Print out the template parameters */
    for (j=0; j<nlist; j++)
    {
	/*
	Retain only those templates that have meaningful masses:
	*/
	    bankParams.x0 = (REAL8) list->data[2*j];
	    bankParams.x1 = (REAL8) list->data[2*j+1];
	    LALInspiralValidParams(&status, &valid, bankParams, coarseIn);
	    if (valid) fprintf(fpr, "%10.4f %10.4f %10.3f %10.3f\n",
			    bankParams.x0, bankParams.x1);
    }
  }
  {
    UINT4 j;
    INT4 valid;

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
  fclose(fpr);
  /* Free the list, and exit. */
  if (list != NULL) LALFree (list);
  LALDDestroyVector(&status, &(shf.data) );
  LALCheckMemoryLeaks();
  return 0;
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

   INITSTATUS (status, "GetInspiralMoments", PNTEMPLATESC);
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

