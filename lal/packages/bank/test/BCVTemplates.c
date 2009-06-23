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
INT4 lalDebugLevel = 0;

int
main(int argc, char **argv)
{
  /* top-level status structure */
  static LALStatus status;
  static InspiralCoarseBankIn coarseIn;
  static InspiralTemplateList *list1, *list2;

  static RectangleIn RectIn;
  static RectangleOut RectOut;
  void (*noisemodel)(LALStatus*,REAL8*,REAL8) = LALLIGOIPsd;

  INT4   j, valid, numPSDpts=262144;
  FILE *fpr;
  INT4 nlist1, nlist2;
/* Number of templates is nlist */


  fpr = fopen("BCVTemplates.out", "w");
  nlist1 = 0;
  nlist2 = 0;
  coarseIn.HighGM =6.;
  coarseIn.LowGM = 3.;
  coarseIn.fLower = 40.L;
  coarseIn.fUpper = 2000.L;
  coarseIn.tSampling = 4096.L;
  coarseIn.order = LAL_PNORDER_TWO;
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
  coarseIn.space       = Psi0Psi3;

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
  coarseIn.approximant  = TaylorT1;
  coarseIn.space	= Tau0Tau3;

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

  /* Print rectagles*/

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

  /* Free the list, and exit. */
  if (list1 != NULL) LALFree (list1);
  if (list2 != NULL) LALFree (list2);
  LALDDestroyVector( &status, &(coarseIn.shf.data) );
  LALCheckMemoryLeaks();
  fclose(fpr);
  return 0;
}
