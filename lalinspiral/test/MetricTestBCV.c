/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
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

/********************************* <lalVerbatim file="MetricTestBCVCV">
Author: B.S. Sathyaprakash
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Program \texttt{MetricTestBCV.c}}
\label{ss:MetricTestBCV.c}

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

\vfill{\footnotesize\input{MetricTestBCVCV}}

******************************************************* </lalLaTeX> */

#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALNoiseModels.h>
#include <lal/AVFactories.h>

#include <lal/LALRCSID.h>
NRCSID (METRICTESTBCVC,"$Id$");

INT4 lalDebugLevel=0;

int
main ( void )
{
  InspiralMetric metric;
  static LALStatus status;
  InspiralTemplate     *params;
  REAL8FrequencySeries    psd;
  void (*noisemodel)(LALStatus*,REAL8*,REAL8) = LALLIGOIPsd;
  UINT4 numPSDpts = 65536;
  REAL8 tSampling;
  REAL8 mismatch;
  FILE *fpr;
  fpr = fopen("MetricTestBCV.out", "w");

  mismatch = 0.03;
  params = (InspiralTemplate *)LALMalloc(sizeof(InspiralTemplate));

  params->alpha = 0.L;
  params->fLower = 30;
  params->fCutoff = 400;

  tSampling = 4096.L;

  memset( &(psd), 0, sizeof(REAL8FrequencySeries) );
  psd.f0 = 0;
  LALDCreateVector(&status, &(psd.data), numPSDpts );
  psd.deltaF = tSampling / (2.L*(REAL8) psd.data->length + 1.L);
  LALNoiseSpectralDensity (&status, psd.data, noisemodel, psd.deltaF );

  LALInspiralComputeMetricBCV(&status, &metric, &psd, params);

  fprintf(fpr, "#%e %e %e\n", metric.G00, metric.G01, metric.G11);
  fprintf(fpr, "#%e %e %e\n", metric.g00, metric.g11, metric.theta);
  fprintf(fpr, "#dp0=%e dp1=%e\n", sqrt (mismatch/metric.G00), sqrt (mismatch/metric.G11));
  fprintf(fpr, "#dP0=%e dP1=%e\n", sqrt (mismatch/metric.g00), sqrt (mismatch/metric.g11));


  {
  double MM;
  double dp0, dp1;
  long n=100;
  double dp0min=-5750;
  double dp0max=5750;
  double dp1min=-220;
  double dp1max=220;
  double d0=(dp0max-dp0min)/(double)n;
  double d1=(dp1max-dp1min)/(double)n;
  for ( dp0= dp0min; dp0<=dp0max ; dp0+=d0)
  {
      for ( dp1= dp1min; dp1<=dp1max ; dp1+=d1)
      {
	  MM = 1. - (metric.G00 * dp0 * dp0 +  metric.G01 * dp0 * dp1
		  +  metric.G01 * dp1 * dp0 +  metric.G11 * dp1 * dp1);
	  fprintf(fpr, "%f %f %f\n", dp0, dp1, MM);
      }
      fprintf(fpr, "\n");
  }
  }
  fclose(fpr);
  LALFree(params);
  LALDDestroyVector(&status, &(psd.data) );
  LALCheckMemoryLeaks();
  return 0;

}
