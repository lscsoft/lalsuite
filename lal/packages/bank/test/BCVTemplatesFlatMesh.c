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
#define FLATMESHTESTC_ENORM 0
#define FLATMESHTESTC_ESUB  1
#define FLATMESHTESTC_EARG  2
#define FLATMESHTESTC_EMEM  3
#define FLATMESHTESTC_EDIM  4
#define FLATMESHTESTC_ELEN  5
#define FLATMESHTESTC_EFILE 6

#define FLATMESHTESTC_MSGENORM "Normal exit"
#define FLATMESHTESTC_MSGESUB  "Subroutine failed"
#define FLATMESHTESTC_MSGEARG  "Error parsing arguments"
#define FLATMESHTESTC_MSGEMEM  "Memory allocation error"
#define FLATMESHTESTC_MSGEDIM  "Inconsistent parameter space dimension"
#define FLATMESHTESTC_MSGELEN  "Too few points specified"
#define FLATMESHTESTC_MSGEFILE "Could not open file"
/******************************************** </lalErrTable><lalLaTeX>

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALPrintError()                 LALCheckMemoryLeaks()
LALCalloc()                     LALFree()
LALCreateFlatMesh()             LALSReadVectorSequence()
LALSCreateVectorSequence()      LALSDestroyVectorSequence()
LALSCreateVector()              LALSDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{BCVTemplatesCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <stdlib.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/FlatMesh.h>
#include <lal/StreamInput.h>

NRCSID(FLATMESHTESTC,"$Id$");

/* Default parameter settings. */
int lalDebugLevel = 0;
#define MISMATCH 0.3
#define DIM 2

/* Usage format string. */
#define USAGE "Usage: %s [-o outfile] [-d debuglevel] [-m mismatch]\n\
                    [eigenvectorfile inversefile rangefile]\n"

/* Macros for printing errors and testing subroutines. */
#define ERROR( code, msg, statement )                                \
if ( lalDebugLevel & LALERROR )                                      \
{                                                                    \
  LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n"   \
		 "        %s %s\n", (code), *argv, __FILE__,         \
		 __LINE__, FLATMESHTESTC, statement ? statement : \
                 "", (msg) );                                        \
}                                                                    \
else (void)(0)

#define INFO( statement )                                            \
if ( lalDebugLevel & LALINFO )                                       \
{                                                                    \
  LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"       \
		 "        %s\n", *argv, __FILE__, __LINE__,          \
		 FLATMESHTESTC, (statement) );                    \
}                                                                    \
else (void)(0)

#define SUB( func, statusptr )                                       \
if ( (func), (statusptr)->statusCode )                               \
{                                                                    \
  ERROR( FLATMESHTESTC_ESUB, FLATMESHTESTC_MSGESUB,            \
         "Function call \"" #func "\" failed:" );                    \
  return FLATMESHTESTC_ESUB;                                      \
}                                                                    \
else (void)(0)

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

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
main(int argc, char **argv)
{
  UINT4 dim;                 /* dimension of parameter space */
  static LALStatus status;     /* top-level status structure */
  REAL4 mismatch = MISMATCH; /* maximum mismatch level */
  REAL4VectorSequence *matrix = NULL;    /* tranformation matrix */
  REAL4VectorSequence *matrixInv = NULL; /* inverse tranformation */
  REAL4VectorSequence *corners = NULL;   /* corners of serach region */
  REAL4VectorSequence *mesh = NULL;      /* mesh of parameter values */
  FILE *fp;                  /* input/output file pointer */

  static InspiralMetric metric;
  static InspiralTemplate params;
  UINT4   nlist, numPSDpts=262144;
  REAL8FrequencySeries shf;
  REAL8 samplingRate;
  void (*noisemodel)(LALStatus*,REAL8*,REAL8) = LALLIGOIPsd;
  InspiralMomentsEtc moments;

  argc = 0;

/* Number of templates is nlist */

  fp = fopen("BCVTemplatesFlatMesh.out", "w");
  dim = DIM;
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
  /*
  params.alpha = 0.528;
  */
  params.alpha = 0.L;
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
  /*
  LALInspiralComputeMetricBCV(&status, &metric, &shf, &params);
  */

  fprintf(fp, "%e %e %e\n", metric.G00, metric.G01, metric.G11);
  fprintf(fp, "%e %e %e\n", metric.g00, metric.g11, metric.theta);
  fprintf(fp, "dp0=%e dp1=%e\n", sqrt (mismatch/metric.G00), sqrt (mismatch/metric.G11));
  fprintf(fp, "dP0=%e dP1=%e\n", sqrt (mismatch/metric.g00), sqrt (mismatch/metric.g11));
  {
    CreateVectorSequenceIn in;
    in.length = dim;
    in.vectorLength = dim;
    SUB( LALSCreateVectorSequence( &status, &matrix, &in ), &status );
    SUB( LALSCreateVectorSequence( &status, &matrixInv, &in ), &status );
    in.length = 2;
    SUB( LALSCreateVectorSequence( &status, &corners, &in ), &status );

    corners->data[0] = 0.3;
    corners->data[1] = 0.15;
    corners->data[2] = 43.;
    corners->data[3] = 1.25;

    /*
    corners->data[0] = 1.e5;
    corners->data[1] = -1.e3;
    corners->data[2] = 2.0e5;
    corners->data[3] = 1.e3;
    */

  }
  {
	  REAL4 det;
	  UINT4 i;



	  matrix->data[0] = -metric.G01/sqrt(pow(metric.G01,2.) + pow(metric.G00-metric.g00,2.))
	  /sqrt(metric.g11);
	  matrix->data[1] = -(metric.G00-metric.g00)/sqrt(pow(metric.G01,2.) + pow(metric.G00-metric.g00,2.))
	  /sqrt(metric.g00);
	  matrix->data[2] = -(metric.G11-metric.g11)/sqrt(pow(metric.G01,2.) + pow(metric.G11-metric.g11,2.))
	  /sqrt(metric.g11);
	  matrix->data[3] = -metric.G01/sqrt(pow(metric.G01,2.) + pow(metric.G00-metric.g00,2.))
	  /sqrt(metric.g00);

	  det = matrix->data[0]*matrix->data[3] - matrix->data[1]*matrix->data[2];

	  matrixInv->data[0] = matrix->data[3]/det;
	  matrixInv->data[1] = -matrix->data[1]/det;
	  matrixInv->data[2] = -matrix->data[2]/det;
	  matrixInv->data[3] = matrix->data[0]/det;

          for (i=0; i<2*dim; i+=2) fprintf(fp, "%e\t%e\n", matrix->data[i], matrix->data[i+1]);
          fprintf(fp, "Det=%e\n", det);
          for (i=0; i<2*dim; i+=2) fprintf(fp, "%e\t%e\n", matrixInv->data[i], matrixInv->data[i+1]);
  }


  /* Apply mismatch threshold to the transformation matrices. */
  {

	  UINT4 i;
	  REAL4 adjust = 2.0*mismatch/sqrt( (REAL4)(dim) );
	  REAL4 *data;  /* pointer to matrix data */

	  i = matrix->length*matrix->vectorLength;
	  data = matrix->data;
	  while ( i-- )
		  *(data++) *= adjust;

	  adjust = 1.0/adjust;
	  i = matrixInv->length*matrixInv->vectorLength;
	  data = matrixInv->data;
	  while ( i-- )
		  *(data++) *= adjust;
  }

  /* Extend the range boundary to ensure edge coverage. */
  {
	  UINT4 i, j;    /* indecies */
	  INT2 direct;  /* sign of direction from first corner to second */
	  REAL4 *data;  /* pointer to matrix data */
	  REAL4 *width; /* maximum width of a patch in each dimension */

	  /* Allocate local memory. */
	  width = (REAL4 *)LALCalloc( dim, sizeof(REAL4) );
	  if ( !width ) {
		  ERROR( FLATMESHTESTC_EMEM, FLATMESHTESTC_MSGEMEM, 0 );
		  return FLATMESHTESTC_EMEM;
	  }

	  /* Determine patch width. */
	  for ( data = matrix->data, i = 0; i < dim; i++ )
		  for ( j = 0; j < dim; j++, data++ )
			  width[j] += fabs( *data );

	  /* Extend each corner by 0.5*width in the appropriate
	     direction. */
	  for ( data = corners->data, i = 0; i < dim; i++, data++ ) {
		  direct = ( data[0] < data[dim] ) ? -1 : 1;
		  data[0] += 0.5*direct*width[i];
		  data[dim] -= 0.5*direct*width[i];
	  }

	  /* Free local memory. */
	  LALFree( width );
  }

  /* Generate mesh using LALFlatMesh() and LALRectIntersect(). */
  {
    /* Set up parameter structure for LALFlatMesh. */
    static FlatMeshParamStruc flatmesh;
    flatmesh.matrix = matrix;
    flatmesh.matrixInv = matrixInv;
    flatmesh.controlPoints = corners;
    /*
    flatmesh.intersection = LALRectIntersect;
    */
    flatmesh.intersection = NULL;
    SUB( LALSCreateVector( &status, &(flatmesh.xMin), dim ), &status );
    SUB( LALSCreateVector( &status, &(flatmesh.xMax), dim ), &status );

    flatmesh.xMin->data[0] = 0.3;
    flatmesh.xMin->data[1] = 0.15;
    flatmesh.xMax->data[0] = 43.;
    flatmesh.xMax->data[1] = 1.25;
    /*
    flatmesh.xMin->data[0] = 1.e5;
    flatmesh.xMin->data[1] = -1.e3;
    flatmesh.xMax->data[0] = 2.0e5;
    flatmesh.xMax->data[1] = 1.e3;
    */


    /* Compute the mesh, and clean up local memory. */
    SUB( LALCreateFlatMesh( &status, &mesh, &flatmesh ), &status );
    SUB( LALSDestroyVector( &status, &(flatmesh.xMin) ), &status );
    SUB( LALSDestroyVector( &status, &(flatmesh.xMax) ), &status );
    SUB( LALSDestroyVectorSequence( &status, &matrix ), &status );
    SUB( LALSDestroyVectorSequence( &status, &matrixInv ), &status );
    SUB( LALSDestroyVectorSequence( &status, &corners ), &status );
  }

  /* Prepare to print result. */
  {
    UINT4 i;
    InspiralBankParams   bankParams;
    InspiralCoarseBankIn coarseIn;
    INT4 valid;
    UINT4 k=0;
    REAL4 *data;

    coarseIn.mmCoarse = 0.70;
    coarseIn.mmFine = 0.97;
    coarseIn.fLower = 40.L;
    coarseIn.fUpper = 2000.L;
    coarseIn.iflso = 0.0L;
    coarseIn.tSampling = 4096.L;
    coarseIn.order = LAL_PNORDER_TWO;
    coarseIn.space = Tau0Tau3;
    coarseIn.approximant = TaylorT1;

    coarseIn.mMin = 1.0;
    coarseIn.mMax = 20.0;
    coarseIn.MMax = coarseIn.mMax * 2.;

    coarseIn.massRange = MinMaxComponentMass;
    /* coarseIn.massRange = MinComponentMassMaxTotalMass;*/

    /* minimum value of eta */
    coarseIn.etamin = coarseIn.mMin * ( coarseIn.MMax - coarseIn.mMin) / pow(coarseIn.MMax,2.);

    /* Print out the template parameters */
    i = mesh->length;
    data = mesh->data;
    while ( i--) {
	/*
	Retain only those templates that have meaningful masses:
	*/
	bankParams.x0 = (REAL8) data[k++];
	bankParams.x1 = (REAL8) data[k++];
	LALInspiralValidParams(&status, &valid, bankParams, coarseIn);
        if (valid) fprintf(fp, "%10.3e %10.3e\n", bankParams.x0, bankParams.x1);
    }
  }


  fclose(fp);
  /* Free the mesh, and exit. */
  SUB( LALSDestroyVectorSequence( &status, &mesh ), &status );
  LALDDestroyVector(&status, &shf.data);
  LALCheckMemoryLeaks();
  INFO( FLATMESHTESTC_MSGENORM );
  return FLATMESHTESTC_ENORM;
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
