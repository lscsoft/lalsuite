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

/********************************* <lalVerbatim file="MetricTestCV">
Author: B.S. Sathyaprakash
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Program \texttt{MetricTest.c}}
\label{ss:MetricTest.c}

Creates a template mesh for BCV (or, alternatively, for SPA but
assuing a constant metric) using the mismatch metric.

\subsubsection*{Usage}

\subsubsection*{Description}

\subsubsection*{Exit codes}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{MetricTestCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <stdlib.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALStdio.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>

NRCSID(METRICTESTC,"$Id$");

int lalDebugLevel = 0;

static void
GetInspiralMoments (
		LALStatus            *status,
		InspiralMomentsEtc   *moments,
		REAL8FrequencySeries *psd,
		InspiralTemplate     *params );

void
LALInspiralComputeBCVMetric(
   LALStatus            *status,
   InspiralMetric       *metric,
   REAL8FrequencySeries *shf,
   InspiralTemplate     *params
);

int
main()
{
  UINT4 dim;                 /* dimension of parameter space */
  static LALStatus status;     /* top-level status structure */

  static InspiralMetric metric;
  static InspiralTemplate params;
  UINT4  nlist, numPSDpts=65536;
  REAL8FrequencySeries shf;
  REAL8 samplingRate;
  void (*noisemodel)(LALStatus*,REAL8*,REAL8) = LALLIGOIPsd;
  InspiralMomentsEtc moments;
  REAL8 mismatch;
  FILE *fpr;

  fpr = fopen("MetricTest.out", "w");

  mismatch = 0.10L;

/* Number of templates is nlist */

  dim = 2;
  nlist = 0;

  params.OmegaS = 0.;
  params.Theta = 0.;
  params.ieta=1;
  params.mass1=1.;
  params.mass2=1.;
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

  params.psi0 = 132250.;
  params.psi3 = -1314.2;
  params.alpha = 0.528;
  params.fFinal = 868.7;
  metric.space = Tau0Tau3;

  samplingRate = params.tSampling;
  memset( &(shf), 0, sizeof(REAL8FrequencySeries) );
  shf.f0 = 0;
  LALDCreateVector( &status, &(shf.data), numPSDpts );
  shf.deltaF = samplingRate / (2.*(REAL8) shf.data->length + 1.L);
  LALNoiseSpectralDensity (&status, shf.data, noisemodel, shf.deltaF );

  /* compute the metric at this point, update bankPars and add the params to the list */

  GetInspiralMoments (&status, &moments, &shf, &params);
  LALInspiralComputeMetric(&status, &metric, &params, &moments);
  fprintf(fpr, "#%e %e %e\n", metric.G00, metric.G01, metric.G11);
  fprintf(fpr, "#%e %e %e\n", metric.g00, metric.g11, metric.theta);
  fprintf(fpr, "#dP0=%e dP1=%e\n", sqrt (mismatch/metric.G00), sqrt (mismatch/metric.G11));
  fprintf(fpr, "#dp0=%e dp1=%e\n", sqrt (mismatch/metric.g00), sqrt (mismatch/metric.g11));
  {
  double MM;
  double dp0, dp1;
  long n=100;
  double dp0min=-0.1;
  double dp0max=0.1;
  double dp1min=-0.1;
  double dp1max=0.1;
  double d0=(dp0max-dp0min)/(double)n;
  double d1=(dp1max-dp1min)/(double)n;
  for ( dp0= dp0min; dp0<=dp0max+d0 ; dp0+=d0)
    {
      for ( dp1= dp1min; dp1<=dp1max+d1 ; dp1+=d1)
	{

	  MM = 1. - (metric.G00 * dp0 * dp0 +  metric.G01 * dp0 * dp1
		  +  metric.G01 * dp1 * dp0 +  metric.G11 * dp1 * dp1);
	  fprintf(fpr, "%f %f %f\n", dp0, dp1, MM);
	}
              fprintf(fpr,"\n");
    }
  }
  fclose(fpr);
  LALDDestroyVector(&status, &(shf.data));
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

   INITSTATUS (status, "GetInspiralMoments", METRICTESTC);
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

