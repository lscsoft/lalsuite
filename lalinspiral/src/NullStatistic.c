/*
*  Copyright (C) 2007 Eirini Messaritaki
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

/*-----------------------------------------------------------------------
 *
 * File Name: NullStatistic.c
 *
 * Author: Messaritaki, E.
 *
 * Revision: $Id$
 *-----------------------------------------------------------------------
 */

#include <math.h>
#include <string.h>

#include <lal/LALRCSID.h>
#include <lal/LALConfig.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/DetectorSite.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/SkyCoordinates.h>
#include <lal/Date.h>
#include <lal/CoherentInspiral.h>
#include <lal/NullStatistic.h>

#define rint(x) (floor((x)+0.5))
double modf( double value, double *integerPart );

NRCSID (NULLSTATISTICC, "$Id$");


static REAL4 cartesianInnerProduct(REAL4 x[3], REAL4 y[3])
{
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}


int
XLALNullStatisticInputInit (
   NullStatInputParams    **input,
   NullStatInitParams      *init
   )
{
  static const char *func = "XLALNullStatisticInputInit";

  CVector               *cVecPtr = NULL;
  NullStatInputParams   *inputPtr = NULL;
  int m;

  /* check that the number of vectors is not zero and that it is 2 */
  if ( init->numDetectors != 2 )
  {
    XLALPrintError("number of detectors must be 2; only H1-H2 allowed");
    XLAL_ERROR(func, XLAL_EINVAL);
  }

  /* check that the number of segments is positive */
  if ( init->numSegments <= 0 )
  {
    XLALPrintError("number of segments must be positive");
    XLAL_ERROR(func, XLAL_EINVAL);
  }

  /* check that the number of points is positive */
  if ( init->numPoints <= 0 )
  {
    XLALPrintError("number of points must be positive");
    XLAL_ERROR(func, XLAL_EINVAL);
  }

  /* allocate memory for the NullStatInputParams structure */
  inputPtr = *input = (NullStatInputParams *)
    LALCalloc(1,sizeof(NullStatInputParams) );
  if ( !inputPtr )
  {
    XLALPrintError("could not allocate memory for NullStatInputParams");
    XLAL_ERROR(func, XLAL_ENOMEM);
  }

  /* allocate memory for the CVector */
  cVecPtr = (*input)->CData = (CVector *) LALCalloc(1, sizeof(CVector));
  if ( !cVecPtr )
  {
    XLALPrintError("could not allocate memory for CVector");
    XLAL_ERROR(func, XLAL_ENOMEM);
  }

  cVecPtr->numDetectors = init->numDetectors;

  for (m=0; m<LAL_NUM_IFO; m++)
  {
    cVecPtr->cData[m] = (COMPLEX8TimeSeries *)
      LALCalloc(1, sizeof(COMPLEX8TimeSeries) );

    if ( !(cVecPtr->cData) )
    {
      XLALPrintError("could not allocate memory for cData vector");
      XLAL_ERROR(func, XLAL_ENOMEM);
    }
  }

  for (m=0; m<LAL_NUM_IFO; m++)
  {
    cVecPtr->cData[m]->data = (COMPLEX8Sequence *)
      LALCalloc(1, sizeof(COMPLEX8Sequence) );

    if ( !(&(cVecPtr->cData[m]->data)) )
    {
      XLALPrintError("could not allocate memory for cData->data vector");
      XLAL_ERROR(func, XLAL_ENOMEM);
    }
  }


  for (m=0; m<LAL_NUM_IFO; m++)
  {
    cVecPtr->cData[m]->data->data = (COMPLEX8 *)
         LALCalloc(1, sizeof(COMPLEX8) );

    if ( !(&(cVecPtr->cData[m]->data->data)) )
    {
      XLALPrintError("could not allocate memory for cData->data->data vector");
      XLAL_ERROR(func, XLAL_ENOMEM);
    }
  }

  for (m=0; m<LAL_NUM_IFO; m++)
  {
    cVecPtr->cData[m]->data->length = init->numPoints;
  }


  return( 0 );
}


int
XLALNullStatisticParamsInit(
   NullStatParams      **output,
   NullStatInitParams   *init
   )
{
  static const char *func = "XLALNullStatisticParamsInit";

  NullStatParams  *outputPtr;


  /* check that the number of vectors is not zero and that it is 2 */
  if ( init->numDetectors != 2 )
  {
    XLALPrintError("number of detectors must be 2; only H1-H2 allowed");
    XLAL_ERROR(func, XLAL_EINVAL);
  }

  /* check that the number of segments is positive */
  if ( init->numSegments <= 0 )
  {
    XLALPrintError("number of segments must be positive");
    XLAL_ERROR(func, XLAL_EINVAL);
  }

  /* check that the number of points is positive */
  if ( init->numPoints <= 0 )
  {
    XLALPrintError("number of points must be positive");
    XLAL_ERROR(func, XLAL_EINVAL);
  }

  /* allocate memory for the NullStatParams structure */
  outputPtr = *output = (NullStatParams *) LALCalloc(1,sizeof(NullStatParams) );  if ( !outputPtr )
  {
    XLALPrintError("could not allocate memory for NullStatParams");
    XLAL_ERROR(func, XLAL_ENOMEM);
  }

  /* if the null statistic is to be output, create the nullStatVec */
  if ( outputPtr->nullStatOut )
  {
    outputPtr->nullStatVec = (REAL4TimeSeries *)
      LALCalloc(1, sizeof(REAL4TimeSeries) );
    if ( !(outputPtr->nullStatVec) )
    {
      XLALPrintError("could not allocate memory for nullStatVec.");
      XLAL_ERROR(func, XLAL_ENOMEM);
    }
#if 0
    outputPtr->nullStatVec->data =
          XLALCreateVector( (init->numPoints * sizeof(REAL4Sequence)) );
    outputPtr->nullStatVec->data->data =
          XLALCreateVector( (init->numPoints * sizeof(REAL4)) );
#endif

/*    outputPtr->nullStatVec->data = (REAL4Sequence *)
       LALCalloc(1, sizeof(REAL4Sequence) ); */
    outputPtr->nullStatVec->data->data = (REAL4 *)
       LALCalloc(1, init->numPoints*sizeof(REAL4) );
  }

  /* populate the output structure */
  outputPtr->numDetectors = init->numDetectors;
  outputPtr->numSegments = init->numSegments;
  outputPtr->numPoints = init->numPoints;

  return (0);
}


int
XLALComputeNullStatistic (
  MultiInspiralTable    **eventList,
  NullStatInputParams    *input,
  NullStatParams         *params
  )
{
  UINT4                 numDetectors          = 0;
  UINT4                 numSegments           = 0;
  UINT4                 numPoints             = 0;
  INT4                  segmentLength         = 0;
  REAL4                 nullStatThresh        = 0.0;
  UINT4                 nullStatOut           = 0;
  REAL4                 nullStatRe            = 0.0 ;
  REAL4                 nullStatIm            = 0.0;
  REAL4                 eventNullStat         = 0.0;
  LALDetector          *detector              = NULL;
  COMPLEX8TimeSeries   *cData[LAL_NUM_IFO];
  MultiInspiralTable   *thisEvent             = NULL;
  REAL4TimeSeries      *nullStatVec           = NULL;
  REAL4                 sigmasq[6];

  INT4    h1idx, h2idx, m, n, j, idx;
  REAL4   norm;

  static const char *func = "XLALComputeNullStatistic";

  /* check that the arguments are reasonable */

  /* assign the parameters to local variables */
  numDetectors   = params->numDetectors;
  numPoints      = params->numPoints;
  numSegments    = params->numSegments;
  nullStatOut    = params->nullStatOut;
  segmentLength  = params->segmentLength;
  nullStatThresh = params->nullStatThresh;

  /* if the null statistic time series is required: */
  if ( nullStatOut)
  {
    memset( params->nullStatVec->data->data, 0, numPoints*sizeof(REAL4));
  }

  h1idx = LAL_IFO_H1;
  h2idx = LAL_IFO_H2;

  norm = ( 1.0 / ( params->sigmasq[h1idx] * params->sigmasq[h1idx] ) ) +
         ( 1.0 / ( params->sigmasq[h2idx] * params->sigmasq[h2idx] ) );

  /* read in the c-data for multiple detectors */
  for (n=0; n<LAL_NUM_IFO; n++)
  {
    cData[n] = input->CData->cData[n];
  }

  /*
   * the time for which we want the null statistic is the one at the
   * center of the 4-s time interval
   */
   idx = (INT4) (numPoints/2);

  /* calculate the null statistic for H1-H2              */
  /* the c-data already contain one division by sigma-sq */
  if (nullStatOut)
  {
    for (m=0; m<(INT4)numPoints; m++)
    {
      nullStatRe =
          (cData[LAL_IFO_H1]->data->data[m].re / params->sigmasq[LAL_IFO_H1])
        - (cData[LAL_IFO_H2]->data->data[m].re / params->sigmasq[LAL_IFO_H2]);
      nullStatIm =
          (cData[LAL_IFO_H1]->data->data[m].im / params->sigmasq[LAL_IFO_H1])
        - (cData[LAL_IFO_H2]->data->data[m].im / params->sigmasq[LAL_IFO_H2]);
      nullStatVec->data->data[m] =
         ( nullStatRe*nullStatRe + nullStatIm*nullStatIm ) / norm ;
    }
    eventNullStat = nullStatVec->data->data[idx];
  }
  else
  {
#if 0
    nullStatRe =
        (input->CData->cData[1]->data->data[idx].re / params->sigmasq[1])
      - (input->CData->cData[2]->data->data[idx].re / params->sigmasq[2]);
#endif
    nullStatRe =
        (cData[LAL_IFO_H1]->data->data[idx].re / params->sigmasq[LAL_IFO_H1])
      - (cData[LAL_IFO_H2]->data->data[idx].re / params->sigmasq[LAL_IFO_H2]);
    nullStatIm =
        (cData[LAL_IFO_H1]->data->data[idx].im / params->sigmasq[LAL_IFO_H1])
      - (cData[LAL_IFO_H2]->data->data[idx].im / params->sigmasq[LAL_IFO_H2]);
    eventNullStat = ( nullStatRe*nullStatRe + nullStatIm*nullStatIm ) / norm ;
  }


  thisEvent->mass1 = input->tmplt->mass1;
  thisEvent->mass2 = input->tmplt->mass2;
  thisEvent->mchirp = input->tmplt->totalMass * pow(input->tmplt->eta, 3.0/5.0);
  thisEvent->eta = input->tmplt->eta;
  thisEvent->null_statistic = eventNullStat;

  return (0);
}


int
XLALNullStatisticParamsFinal(
   NullStatParams      **output
   )
{
  static const char *func = "XLALNullStatisticParamsFinal";

  NullStatParams  *outputPtr;

  outputPtr = *output ;

  /* destroy detector vector */
  LALFree( outputPtr->detVector->detector );
  outputPtr->detVector->detector = NULL;
  LALFree( outputPtr->detVector );
  outputPtr->detVector = NULL;

  /* destroy null statistic vector, if it exists */
  if ( outputPtr->nullStatOut )
  {
    XLALDestroyVector( outputPtr->nullStatVec->data );
    LALFree( outputPtr->nullStatVec );
  }

  LALFree( outputPtr );
  *output = NULL;

  return (0);
}

int
XLALNullStatisticInputFinal (
   NullStatInputParams    **input
   )
{
  static const char *func = "XLALNullStatisticInputFinal";

  CVector               *cVecPtr;
  NullStatInputParams   *inputPtr;
  INT4                   l;

  inputPtr = *input;

  /* free the template memory */
  LALFree( inputPtr->tmplt );
  inputPtr->tmplt = NULL;

  /* destroy the cVector if necessary */
  if ( inputPtr->CData )
  {
    cVecPtr = (*input)->CData ;

    for ( l=0; l < cVecPtr->numDetectors; l++)
    {
      if ( cVecPtr->cData[l]->data != NULL )
      {
        XLALDestroyVector( cVecPtr->cData[l]->data );
      }
    }

    LALFree( cVecPtr->cData );
    LALFree (cVecPtr);
    (*input)->CData = NULL;

  }

  LALFree( inputPtr );
  *input = NULL;

  return( 0 );
}

