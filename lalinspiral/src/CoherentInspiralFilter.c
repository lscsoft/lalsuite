/*
*  Copyright (C) 2007 Sukanta Bose, Sean Seader
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
 * File Name: CoherentInspiralFilter.c
 *
 * Author: Bose, S., Seader, S. E.
 *
 *-----------------------------------------------------------------------
 */

#include <math.h>
#include <string.h>

#include <lal/LALConfig.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/DetectorSite.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/SkyCoordinates.h>
#include <lal/Date.h>
#include <lal/TimeDelay.h>
#include <lal/DetResponse.h>
#include <lal/CoherentInspiral.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#define rint(x) (floor((x)+0.5))

double modf( double value, double *integerPart );
int compare( const void* a, const void* b );
void XLALAssignEventQuads(INT4 caseID[6], COMPLEX8 quadTemp[6], MultiInspiralTable *thisEvent);
void XLALSaveStatsInEventTable(REAL4 cohSnrSqLocal,REAL4 autoCorrCohSq,REAL4 crossCorrCohSq,REAL4 autoCorrNullSq,REAL4 crossCorrNullSq,REAL8 ampMetricEigenVal1,REAL8 ampMetricEigenVal2,REAL4 VVPlus[4],REAL4 VVMinus[4],MultiInspiralTable *thisEvent);
REAL4 XLALComputeCohStatCase3b(COMPLEX8TimeSeries *cData[4], REAL4 VVPlus[4], REAL4 VVMinus[4], INT4 timeP, INT4 slidePoints[3], UINT4 numDetectors, REAL8 MM1, REAL8 MM2, REAL4 chisq[4], REAL4 chisq_dof[4], REAL4 eff_snr_denom_fac, REAL4 *eff_coh_snr);
REAL4 XLALComputeCrossCohStatCase3b(COMPLEX8TimeSeries *cData[4], REAL4 VVPlus[4], REAL4 VVMinus[4], INT4 timeP, INT4 slidePoints[3], UINT4 numDetectors, REAL8 MM1, REAL8 MM2, REAL4 *autoCorrCohSq, REAL4 chisq[4], REAL4 chisq_dof[4], REAL4 eff_snr_denom_fac, REAL4 *crossCorrCohSqEff);
double XLALComputeIncohStatCase3b(INT4 caseID[6], double fplus[4], double fcross[4], REAL8 *sigmasq, INT4 timeP, INT4 slidePoints[3], COMPLEX8TimeSeries *cData[4], REAL4 nullStatistic, REAL4 *crossCorrNullSq, REAL4 chisq[4], REAL4 chisq_dof[4], REAL4 eff_snr_denom_fac);
double XLALComputeNonQuadNullStatCase3b(INT4 caseID[6], double fplus[4], double fcross[4], REAL8 *sigmasq, INT4 timeP, INT4 slidePoints[3], COMPLEX8TimeSeries *cData[4]);
double XLALComputeNullStatCase3b(INT4 caseID[6],double fplus[4], double fcross[4], REAL8 *sigmasq, MultiInspiralTable *thisEvent);
double XLALComputeNullTimeSeriesCase3b(INT4 caseID[6], double fplus[4], double fcross[4], REAL8 *sigmasq, COMPLEX8 quadTemp[6]);
REAL4 XLALComputeCohStatCase4a(COMPLEX8TimeSeries *cData[4], REAL4 VVPlus[4], REAL4 VVMinus[4], INT4 timeP, INT4 slidePoints4D[4], UINT4 numDetectors, REAL8 MM1, REAL8 MM2, REAL4 chisq[4], REAL4 chisq_dof[4], REAL4 eff_snr_denom_fac, REAL4 *eff_coh_snr);
REAL4 XLALComputeCrossCohStatCase4a(COMPLEX8TimeSeries *cData[4], REAL4 VVPlus[4], REAL4 VVMinus[4], INT4 timeP, INT4 slidePoints4D[4], UINT4 numDetectors, REAL8 MM1, REAL8 MM2, REAL4 *autoCorrCohSq, REAL4 chisq[4], REAL4 chisq_dof[4], REAL4 eff_snr_denom_fac, REAL4 *crossCorrCohSqEff);
double XLALComputeIncohStatCase4a(INT4 caseID[6], double fplus[4], double fcross[4], REAL8 *sigmasq, INT4 timeP, INT4 slidePoints4D[4], COMPLEX8TimeSeries *cData[4], REAL4 nullStatistic, REAL4 *crossCorrNullSq, REAL4 chisq[4], REAL4 chisq_dof[4], REAL4 eff_snr_denom_fac);
double XLALComputeNonQuadNullStatCase4a(double fplus[4], double fcross[4], REAL8 *sigmasq, INT4 timeP, INT4 slidePoints4D[4], COMPLEX8TimeSeries *cData[4]);
double XLALComputeNullStatCase4a(double fplus[4], double fcross[4], REAL8 *sigmasq, MultiInspiralTable *thisEvent);
double XLALComputeNullTimeSeriesCase4a(INT4 caseID[6], double fplus[4], double fcross[4], REAL8 *sigmasq, COMPLEX8 quadTemp[6]);
void XLALCoherentCBCEstimateDistanceCase2a(double C_Real0, double C_Im0, double C_Real1, double C_Im1, REAL8 sigmasq4DArray[4], REAL4 chisq4DArray[4], MultiInspiralTable *thisEvent);
void XLALCoherentCBCEstimateDistanceCase2b(INT4 caseID[6], double C_Real0, double C_Im0, double C_Real1, double C_Im1, REAL8 sigmasq4DArray[4], REAL4 chisq4DArray[4], MultiInspiralTable *thisEvent);
void XLALCoherentCBCEstimateDistanceCase3a(INT4 caseID[6], double C_Real0, double C_Im0, double C_Real1,double C_Im1, double C_Real2, double C_Im2, REAL8 sigmasq[4], REAL4 chisq[4], MultiInspiralTable *thisEvent);
void XLALCoherentCBCSigmasqFor3Sites(INT4 caseID[6], REAL8 sigmasq4DArray[4], REAL4 chisq4DArray[4], MultiInspiralTable *thisEvent, UINT4 case3b);
double XLALCoherentCBCParamEstim( double *psi_est, double *iota_est, double *coa_phase_est, double a1, double a2, double a3, double a4, double amplitudeConst, MultiInspiralTable *thisEvent, double chirpTime, double C_Real0, double C_Im0, double C_Real1,double C_Im1, double C_Real2, double C_Im2, double C_Real3, double C_Im3, REAL8 sigmasq[4], UINT4 case3b, INT4 caseID[6] );

static REAL4 cartesianInnerProduct(REAL4 x[3], REAL4 y[3])
{
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}


void
LALCoherentInspiralFilterInputInit (
    LALStatus                       *status,
    CoherentInspiralFilterInput    **input,
    CoherentInspiralInitParams      *params
    )
{
  UINT4                            i, l;
  UINT4                            length = 4; /*length of thetaPhiVs vectors*/
  CoherentInspiralFilterInput     *inputPtr = NULL;
  CoherentInspiralBeamVector      *beamVecPtr = NULL;
  DetectorBeamArray               *detBeamArray = NULL;
  CoherentInspiralCVector         *cVecPtr = NULL;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  /*
   *
   * ensure that the arguments are reasonable
   *
   */

  /* check that the output handle exists, but points to a null pointer */
  ASSERT( input, status, COHERENTINSPIRALH_ENULL, COHERENTINSPIRALH_MSGENULL );
  ASSERT( !*input, status, COHERENTINSPIRALH_ENNUL, COHERENTINSPIRALH_MSGENNUL );

  /* check that the parameter structure exists */
  ASSERT( params, status, COHERENTINSPIRALH_ENULL, COHERENTINSPIRALH_MSGENULL );

  /* check that the number of detectors is greater than 1 */
  ASSERT( params->numDetectors > 1, status,
      COHERENTINSPIRALH_EZDET, COHERENTINSPIRALH_MSGEZDET );

  /* check that the number of detectors is less than 5 */
  ASSERT( params->numDetectors < 5, status,
	  COHERENTINSPIRALH_EZDET, COHERENTINSPIRALH_MSGEZDET );

  /* check that the number of data segments is positive */
  ASSERT( params->numSegments > 0, status,
      COHERENTINSPIRALH_ESEGZ, COHERENTINSPIRALH_MSGESEGZ );

  /* check that the number of points is positive */
  ASSERT( params->numPoints > 0, status,
	  COHERENTINSPIRALH_ENUMZ, COHERENTINSPIRALH_MSGENUMZ );


  inputPtr= *input = (CoherentInspiralFilterInput *)
    LALCalloc(1, sizeof(CoherentInspiralFilterInput) );
  if ( !inputPtr )
    {
      ABORT( status, COHERENTINSPIRALH_EALOC, COHERENTINSPIRALH_MSGEALOC);
    }

  /*
   *
   * create the CoherentInspiralBeamVector structure if needed
   *
   */


  if( params->numBeamPoints != 0 )
    {
      /* allocate memory for an array  */
      beamVecPtr = inputPtr->beamVec = (CoherentInspiralBeamVector *)
	LALCalloc(1, sizeof(CoherentInspiralBeamVector));

      if ( !beamVecPtr )
	{
	  LALFree( inputPtr);
	  *input = NULL;
	  ABORT( status, COHERENTINSPIRALH_EALOC, COHERENTINSPIRALH_MSGEALOC);
	}

      /* set the number of detectors in the vector */
      beamVecPtr->numDetectors = params->numDetectors;

      detBeamArray = beamVecPtr->detBeamArray = (DetectorBeamArray *)
	LALCalloc(1, beamVecPtr->numDetectors*sizeof(DetectorBeamArray));
      if ( !detBeamArray )
	{
	  LALFree( beamVecPtr );
	  beamVecPtr = NULL;
	  inputPtr->beamVec = NULL;
	  LALFree( inputPtr );
	  *input = NULL;
	  ABORT( status, COHERENTINSPIRALH_EALOC, COHERENTINSPIRALH_MSGEALOC);
	}


      /* set the number of beamPoints in the vector */
      for ( l = 0 ; l < beamVecPtr->numDetectors ; l++) {
	detBeamArray[l].numBeamPoints = params->numBeamPoints;

	detBeamArray[l].thetaPhiVs = (REAL4TimeSeries *)
	  LALCalloc(1, params->numBeamPoints*sizeof(REAL4TimeSeries));

	BEGINFAIL(status) {
	  LALFree( detBeamArray );
	  beamVecPtr->detBeamArray = NULL;
	  LALFree( beamVecPtr );
	  inputPtr->beamVec = NULL;
	  LALFree( inputPtr );
	  *input = NULL;
	}
	ENDFAIL( status );
      }

      /* set the number of fields in the thetaPhiVs to 4: theta,phi,v+,and v-*/
      for ( l = 0 ; l < beamVecPtr->numDetectors ; l++) {
	for ( i = 0 ; i < detBeamArray[l].numBeamPoints ; i++ ) {
	  LALCreateVector(status->statusPtr, &(detBeamArray[l].thetaPhiVs[i].data), length);

	  BEGINFAIL( status ) {
	    for( l = 0; l < beamVecPtr->numDetectors; l++) {
	      LALFree( detBeamArray[l].thetaPhiVs );
	    }
	    LALFree( detBeamArray );
	    beamVecPtr->detBeamArray = NULL;
	    LALFree( beamVecPtr );
	    inputPtr->beamVec = NULL;
	    LALFree( inputPtr );
	    *input = NULL;
	  }
	  ENDFAIL( status );
	}
      }

    }

  /*
   *
   * create the CoherentInspiralCVector structure
   *
   */

  /* allocate memory for the cData vector  */

  cVecPtr = (*input)->multiCData = (CoherentInspiralCVector *)
    LALCalloc(1, sizeof(CoherentInspiralCVector));
  if ( !cVecPtr )
    {
      LALFree( inputPtr );
      *input = NULL;
      ABORT( status, COHERENTINSPIRALH_EALOC, COHERENTINSPIRALH_MSGEALOC);
    }

  BEGINFAIL( status ) {
    LALFree( cVecPtr );
    (*input)->multiCData = NULL;
    LALFree( inputPtr );
    *input = NULL;
  }
  ENDFAIL( status );

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

void
LALCoherentInspiralFilterInputFinalize (
    LALStatus                       *status,
    CoherentInspiralFilterInput    **input
    )
{
  UINT4                            i,l;
  CoherentInspiralBeamVector      *beamVecPtr;
  DetectorBeamArray               *detBeamArray;
  CoherentInspiralFilterInput     *inputPtr;
  CoherentInspiralCVector         *cVecPtr;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );


  /*
   * check that the arguments are reasonable
   *
   */

  ASSERT( input, status, COHERENTINSPIRALH_ENULL, COHERENTINSPIRALH_MSGENULL );
  ASSERT( *input, status, COHERENTINSPIRALH_ENULL, COHERENTINSPIRALH_MSGENULL );


  inputPtr= *input;

  /*
   *
   * destroy the contents of the CoherentInspiralCVector structure if necessary
   *
   */

  if( inputPtr->multiCData )
    {
      cVecPtr = (*input)->multiCData;
      LALFree( cVecPtr );
      (*input)->multiCData = NULL;

    }

  /*
   *
   * destroy the contents of the CoherentInspiralBeamVector structure if needed
   *
   */

  if( inputPtr->beamVec )
    {
      beamVecPtr = inputPtr->beamVec;
      detBeamArray = beamVecPtr->detBeamArray;

      /* destroy vector  */
      if ( beamVecPtr != NULL ) {
	for ( l = 0 ; l < beamVecPtr->numDetectors ; l++ ) {
	  for  ( i = 0 ; i < detBeamArray[l].numBeamPoints ; i++) {
	    if (inputPtr->beamVec->detBeamArray[l].thetaPhiVs != NULL ) {
	      LALDestroyVector (status->statusPtr, &(detBeamArray[l].thetaPhiVs[i].data));
	      CHECKSTATUSPTR( status );
	    }
	  }
	  LALFree(detBeamArray[l].thetaPhiVs);
	  detBeamArray[l].thetaPhiVs = NULL;
	}
      }

      LALFree(detBeamArray);
      beamVecPtr->detBeamArray = NULL;
      LALFree(beamVecPtr);
      inputPtr->beamVec = NULL;
    }

  LALFree( inputPtr );
  *input = NULL;

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );

}


void
LALCoherentInspiralFilterParamsInit (
    LALStatus                       *status,
    CoherentInspiralFilterParams   **output,
    CoherentInspiralInitParams      *params
    )
{
  CoherentInspiralFilterParams     *outputPtr;
  INT4                              networkLength = LAL_NUM_IFO;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );



  /*
   *
   * ensure that the arguments are reasonable
   *
   */


  /* check that the output handle exists, but points to a null pointer */
  ASSERT( output, status, COHERENTINSPIRALH_ENULL, COHERENTINSPIRALH_MSGENULL );
  ASSERT( !*output, status, COHERENTINSPIRALH_ENNUL, COHERENTINSPIRALH_MSGENNUL );

  /* check that the parameter structure exists */
  ASSERT( params, status, COHERENTINSPIRALH_ENULL, COHERENTINSPIRALH_MSGENULL );
 /* check that the number of detectors is greater than 1 */
  ASSERT( params->numDetectors > 1, status,
      COHERENTINSPIRALH_EZDET, COHERENTINSPIRALH_MSGEZDET );

  /* check that the number of detectors is less than 5 */
  ASSERT( params->numDetectors < 5, status,
	  COHERENTINSPIRALH_EZDET, COHERENTINSPIRALH_MSGEZDET );

  /* check that the number of data segments is positive */
  ASSERT( params->numSegments > 0, status,
      COHERENTINSPIRALH_ESEGZ, COHERENTINSPIRALH_MSGESEGZ );

  /* check that the number of points is positive */
  ASSERT( params->numPoints > 0, status,
      COHERENTINSPIRALH_ENUMZ, COHERENTINSPIRALH_MSGENUMZ );


  /*
   *
   * allocate memory for the FindChirpFilterParams
   *
   */


  /* create the output structure */
  outputPtr= *output = (CoherentInspiralFilterParams *)
    LALCalloc(1, sizeof(CoherentInspiralFilterParams) );
  if (!outputPtr )
    {
      ABORT( status, COHERENTINSPIRALH_EALOC, COHERENTINSPIRALH_MSGEALOC);
    }

  outputPtr->numDetectors = params->numDetectors;
  outputPtr->numSegments = params->numSegments;
  outputPtr->numPoints = params->numPoints;
  outputPtr->numBeamPoints = params->numBeamPoints;

  LALU2CreateVector(status->statusPtr, &(outputPtr->detIDVec),
		     networkLength);
  BEGINFAIL( status ) {
    LALFree( outputPtr );
    *output = NULL;
  }
  ENDFAIL( status );



  /*
   *
   * create vector to store coherent SNR if outputting it
   *
   */

  if( params->cohSNROut )
    {
        if (params->threeSiteCase) {
	  outputPtr->cohSNRVec3Sites = (REAL4TimeSeries *)
	    LALCalloc( 1, sizeof(REAL4TimeSeries) );
	  LALCreateVector (status->statusPtr, &(outputPtr->cohSNRVec3Sites->data),
			   params->numBeamPoints);
	}
	else {
	  outputPtr->cohSNRVec = (REAL4TimeSeries *)
	    LALCalloc( 1, sizeof(REAL4TimeSeries) );
	  LALCreateVector (status->statusPtr, &(outputPtr->cohSNRVec->data),
			   params->numPoints);
	}
        BEGINFAIL( status ) {
	  if( outputPtr->detIDVec )
	    {
	      TRY( LALU2DestroyVector (status->statusPtr, &(outputPtr->detIDVec)),
		   status);
	      LALFree( outputPtr->detIDVec );
	      outputPtr->detIDVec = NULL;
	    }
          if (params->threeSiteCase) {
	    LALFree( outputPtr->cohSNRVec3Sites );
	    outputPtr->cohSNRVec3Sites = NULL;
	  }
	  else {
	    LALFree( outputPtr->cohSNRVec );
	    outputPtr->cohSNRVec = NULL;
	  }
          LALFree( outputPtr->detectorVec );
          outputPtr->detectorVec = NULL;
          LALFree( outputPtr );
          *output = NULL;
        }
        ENDFAIL( status );
    }

  /*
   *
   * create vector to store H1-H2 coherent-snr if that pair is in
   * and numDetectors > 2 and if the user wants it
   *
   */

  if( params->cohH1H2SNROut )
    {
      outputPtr->cohH1H2SNRVec = (REAL4TimeSeries *)
	LALCalloc( 1, sizeof(REAL4TimeSeries) );
      LALCreateVector (status->statusPtr, &(outputPtr->cohH1H2SNRVec->data),
		       params->numPoints);
      BEGINFAIL( status ) {
	if( params->cohSNROut )
	  {
	    if (params->threeSiteCase) {
	      TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec3Sites->data)),
		   status);
	      LALFree( outputPtr->cohSNRVec3Sites );
	      outputPtr->cohSNRVec3Sites = NULL;
	    }
	    else {
	      TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec->data)),
		   status);
	      LALFree( outputPtr->cohSNRVec );
	      outputPtr->cohSNRVec = NULL;
	    }
	  }
	if( outputPtr->detIDVec )
	  {
	    TRY( LALU2DestroyVector (status->statusPtr, &(outputPtr->detIDVec)),
		 status);
	    LALFree( outputPtr->detIDVec );
	    outputPtr->detIDVec = NULL;
	  }
	LALFree( outputPtr->cohH1H2SNRVec );
	outputPtr->cohH1H2SNRVec = NULL;
	LALFree( outputPtr->detectorVec );
	outputPtr->detectorVec = NULL;
	LALFree( outputPtr );
	*output = NULL;
      }
      ENDFAIL( status );
    }

  /*
   *
   * create vector to store H1-H2 null statistic if outputting it
   *
   */

  if( params->nullStatOut )
    {
      if ( params->threeSiteCase ) {
        outputPtr->nullStatVec3Sites = (REAL4TimeSeries *)
          LALCalloc( 1, sizeof(REAL4TimeSeries) );
        LALCreateVector (status->statusPtr, &(outputPtr->nullStatVec3Sites->data),
			  params->numBeamPoints);
      }
      else {
        outputPtr->nullStatVec = (REAL4TimeSeries *)
          LALCalloc( 1, sizeof(REAL4TimeSeries) );
        LALCreateVector (status->statusPtr, &(outputPtr->nullStatVec->data),
		         params->numPoints);
      }
      BEGINFAIL( status ) {
	  if( params->cohH1H2SNROut )
	    {
	      TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohH1H2SNRVec->data)),
		   status);
	      LALFree( outputPtr->cohH1H2SNRVec );
	      outputPtr->cohH1H2SNRVec = NULL;
	    }
	  if( params->cohSNROut )
	    {
	      if (params->threeSiteCase) {
		TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec3Sites->data)),
		     status);
		LALFree( outputPtr->cohSNRVec3Sites );
		outputPtr->cohSNRVec3Sites = NULL;
	      }
	      else {
		TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec->data)),
		     status);
		LALFree( outputPtr->cohSNRVec );
		outputPtr->cohSNRVec = NULL;
	      }
	    }
	  if( outputPtr->detIDVec )
	    {
	      TRY( LALU2DestroyVector (status->statusPtr, &(outputPtr->detIDVec)),
		   status);
	      LALFree( outputPtr->detIDVec );
	      outputPtr->detIDVec = NULL;
	    }
	  if ( params->threeSiteCase ) {
	    LALFree( outputPtr->nullStatVec3Sites );
	    outputPtr->nullStatVec3Sites = NULL;
	  }
	  else {
	    LALFree( outputPtr->nullStatVec );
	    outputPtr->nullStatVec = NULL;
	  }
          LALFree( outputPtr->detectorVec );
          outputPtr->detectorVec = NULL;
          LALFree( outputPtr );
          *output = NULL;
      }
      ENDFAIL( status );
    }

  /*
   *
   * create detector vector to store detector site IDs
   *
   */


  outputPtr->detectorVec = (DetectorVector *)
    LALCalloc( 1, sizeof(DetectorVector) );

  if ( !(outputPtr->detectorVec) ) {
    ABORT( status, COHERENTINSPIRALH_EALOC, COHERENTINSPIRALH_MSGEALOC);
    BEGINFAIL( status ) {
      if( params->nullStatOut )
	{
	  if ( params->threeSiteCase ) {
	    TRY( LALDestroyVector( status->statusPtr, &(outputPtr->nullStatVec3Sites->data)),
		 status);
	    LALFree( outputPtr->nullStatVec3Sites );
	    outputPtr->nullStatVec3Sites = NULL;

	  }
	  else {
	    TRY( LALDestroyVector( status->statusPtr, &(outputPtr->nullStatVec->data)),
		 status);
	    LALFree( outputPtr->nullStatVec );
	    outputPtr->nullStatVec = NULL;
	  }
	}
      if( params->cohH1H2SNROut )
	{
	  TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohH1H2SNRVec->data)),
	       status);
	  LALFree( outputPtr->cohH1H2SNRVec );
	  outputPtr->cohH1H2SNRVec = NULL;
	}
      if( params->cohSNROut )
	{
	  if (params->threeSiteCase) {
	    TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec3Sites->data)),
		 status);
	    LALFree( outputPtr->cohSNRVec3Sites );
	    outputPtr->cohSNRVec3Sites = NULL;
	  }
	  else {
	    TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec->data)),
		 status);
	    LALFree( outputPtr->cohSNRVec );
	    outputPtr->cohSNRVec = NULL;
	  }
	}
      if( outputPtr->detIDVec )
	{
	  TRY( LALU2DestroyVector (status->statusPtr, &(outputPtr->detIDVec)),
	       status);
	  LALFree( outputPtr->detIDVec );
	  outputPtr->detIDVec = NULL;
	}
      LALFree( outputPtr->detectorVec );
      outputPtr->detectorVec = NULL;
      LALFree( outputPtr );
      *output = NULL;
    }
    ENDFAIL( status );
  }

  outputPtr->detectorVec->detector = (LALDetector *)
    LALCalloc( 1, outputPtr->numDetectors*sizeof(LALDetector));

  if ( !(outputPtr->detectorVec->detector ) ) {
    ABORT( status, COHERENTINSPIRALH_EALOC, COHERENTINSPIRALH_MSGEALOC);
    BEGINFAIL( status ) {
      LALFree( outputPtr->detectorVec );
      if( params->nullStatOut )
	{
	  if ( params->threeSiteCase ) {
	    TRY( LALDestroyVector( status->statusPtr, &(outputPtr->nullStatVec3Sites->data)),
		 status);
	    LALFree( outputPtr->nullStatVec3Sites );
	    outputPtr->nullStatVec3Sites = NULL;
	  }
	  else {
	    TRY( LALDestroyVector( status->statusPtr, &(outputPtr->nullStatVec->data)),
		 status);
	    LALFree( outputPtr->nullStatVec );
	    outputPtr->nullStatVec = NULL;
	  }
	}
      if( params->cohH1H2SNROut )
	{
	  TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohH1H2SNRVec->data)),
	       status);
	  LALFree( outputPtr->cohH1H2SNRVec );
	  outputPtr->cohH1H2SNRVec = NULL;
	}
      if( params->cohSNROut )
	{
	  if (params->threeSiteCase) {
	    TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec3Sites->data)),
		 status);
	    LALFree( outputPtr->cohSNRVec3Sites );
	    outputPtr->cohSNRVec3Sites = NULL;
	  }
	  else {
	    TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec->data)),
		 status);
	    LALFree( outputPtr->cohSNRVec );
	    outputPtr->cohSNRVec = NULL;
	  }
	}
      if( outputPtr->detIDVec )
	{
	  TRY( LALU2DestroyVector (status->statusPtr, &(outputPtr->detIDVec)),
	       status);
	  LALFree( outputPtr->detIDVec );
	  outputPtr->detIDVec = NULL;
	}
      LALFree( outputPtr->detectorVec );
      outputPtr->detectorVec = NULL;
      LALFree( outputPtr );
      *output = NULL;
    }
    ENDFAIL( status );
  }

  /* Fixed the maximum number of participating ifos to be LAL_NUM_IFO */
  LALDCreateVector(status->statusPtr, &(outputPtr->sigmasqVec), LAL_NUM_IFO );
  if ( !(outputPtr->sigmasqVec ) ) {
    ABORT( status, COHERENTINSPIRALH_EALOC, COHERENTINSPIRALH_MSGEALOC);
    BEGINFAIL( status ) {
      if( params->nullStatOut )
	{
	  if ( params->threeSiteCase ) {
	    TRY( LALDestroyVector( status->statusPtr, &(outputPtr->nullStatVec3Sites->data)),
		 status);
	    LALFree( outputPtr->nullStatVec3Sites );
	    outputPtr->nullStatVec3Sites = NULL;
	  }
	  else {
	    TRY( LALDestroyVector( status->statusPtr, &(outputPtr->nullStatVec->data)),
		 status);
	    LALFree( outputPtr->nullStatVec );
	    outputPtr->nullStatVec = NULL;
	  }
	}
      if( params->cohH1H2SNROut )
	{
	  TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohH1H2SNRVec->data)),
	       status);
	  LALFree( outputPtr->cohH1H2SNRVec );
	  outputPtr->cohH1H2SNRVec = NULL;
	}
      if( params->cohSNROut )
	{
	  if (params->threeSiteCase) {
	    TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec3Sites->data)),
		 status);
	    LALFree( outputPtr->cohSNRVec3Sites );
	    outputPtr->cohSNRVec3Sites = NULL;
	  }
	  else {
	    TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec->data)),
		 status);
	    LALFree( outputPtr->cohSNRVec );
	    outputPtr->cohSNRVec = NULL;
	  }
	}
      if( outputPtr->detIDVec )
	{
	  TRY( LALU2DestroyVector (status->statusPtr, &(outputPtr->detIDVec)),
	       status);
	  LALFree( outputPtr->detIDVec );
	  outputPtr->detIDVec = NULL;
	}
      LALFree( outputPtr->detectorVec );
      outputPtr->detectorVec = NULL;
      LALFree( outputPtr );
      *output = NULL;
    }
    ENDFAIL( status );
  }

  /* Fixed the maximum number of participating ifos to be LAL_NUM_IFO */
  LALSCreateVector(status->statusPtr, &(outputPtr->chisqVec), LAL_NUM_IFO );
  if ( !(outputPtr->chisqVec ) ) {
    ABORT( status, COHERENTINSPIRALH_EALOC, COHERENTINSPIRALH_MSGEALOC);
    BEGINFAIL( status ) {
      if( params->nullStatOut )
	{
	  if ( params->threeSiteCase ) {
	    TRY( LALDestroyVector( status->statusPtr, &(outputPtr->nullStatVec3Sites->data)),
		 status);
	    LALFree( outputPtr->nullStatVec3Sites );
	    outputPtr->nullStatVec3Sites = NULL;
	  }
	  else {
	    TRY( LALDestroyVector( status->statusPtr, &(outputPtr->nullStatVec->data)),
		 status);
	    LALFree( outputPtr->nullStatVec );
	    outputPtr->nullStatVec = NULL;
	  }
	}
      if( params->cohH1H2SNROut )
	{
	  TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohH1H2SNRVec->data)),
	       status);
	  LALFree( outputPtr->cohH1H2SNRVec );
	  outputPtr->cohH1H2SNRVec = NULL;
	}
      if( params->cohSNROut )
	{
	  if (params->threeSiteCase) {
	    TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec3Sites->data)),
		 status);
	    LALFree( outputPtr->cohSNRVec3Sites );
	    outputPtr->cohSNRVec3Sites = NULL;
	  }
	  else {
	    TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec->data)),
		 status);
	    LALFree( outputPtr->cohSNRVec );
	    outputPtr->cohSNRVec = NULL;
	  }
	}
      if( outputPtr->detIDVec )
	{
	  TRY( LALU2DestroyVector (status->statusPtr, &(outputPtr->detIDVec)),
	       status);
	  LALFree( outputPtr->detIDVec );
	  outputPtr->detIDVec = NULL;
	}
      if( outputPtr->sigmasqVec )
	{
	  TRY( LALDDestroyVector( status->statusPtr, &(outputPtr->sigmasqVec)),
	       status);
	  LALFree( outputPtr->sigmasqVec );
	  outputPtr->sigmasqVec = NULL;
	}
      LALFree( outputPtr->detectorVec );
      outputPtr->detectorVec = NULL;
      LALFree( outputPtr );
      *output = NULL;
    }
    ENDFAIL( status );
  }

  /* Allocate memory for H1H2 null-statistic, if required */
  if( params->nullStatH1H2Out )
    {
      outputPtr->nullStatH1H2Vec = (REAL4TimeSeries *)
	LALCalloc( 1, sizeof(REAL4TimeSeries) );
      LALCreateVector (status->statusPtr, &(outputPtr->nullStatH1H2Vec->data),
		       params->numPoints);
      BEGINFAIL( status ) {
       if( params->nullStatOut )
       {
	 if ( params->threeSiteCase ) {
	   TRY( LALDestroyVector( status->statusPtr, &(outputPtr->nullStatVec3Sites->data)),
		status);
	   LALFree( outputPtr->nullStatVec3Sites );
	   outputPtr->nullStatVec3Sites = NULL;
	 }
	 else {
	   TRY( LALDestroyVector( status->statusPtr, &(outputPtr->nullStatVec->data)),
		status);
	   LALFree( outputPtr->nullStatVec );
	   outputPtr->nullStatVec = NULL;
	 }
       }
       if( params->cohH1H2SNROut )
       {
	 TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohH1H2SNRVec->data)),
	      status);
	 LALFree( outputPtr->cohH1H2SNRVec );
	 outputPtr->cohH1H2SNRVec = NULL;
       }
       if( params->cohSNROut )
       {
	 if (params->threeSiteCase) {
	   TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec3Sites->data)),
		status);
	   LALFree( outputPtr->cohSNRVec3Sites );
	   outputPtr->cohSNRVec3Sites = NULL;
	 }
	 else {
	   TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec->data)),
		status);
	   LALFree( outputPtr->cohSNRVec );
	   outputPtr->cohSNRVec = NULL;
	 }
       }
       if( outputPtr->detIDVec )
       {
	 TRY( LALU2DestroyVector (status->statusPtr, &(outputPtr->detIDVec)),
	      status);
	 LALFree( outputPtr->detIDVec );
	 outputPtr->detIDVec = NULL;
       }
       if( outputPtr->sigmasqVec )
       {
	 TRY( LALDDestroyVector( status->statusPtr, &(outputPtr->sigmasqVec)),
	      status);
	 LALFree( outputPtr->sigmasqVec );
	 outputPtr->sigmasqVec = NULL;
       }
       if( outputPtr->chisqVec )
       {
	 TRY( LALSDestroyVector( status->statusPtr, &(outputPtr->chisqVec)),
	      status);
	 LALFree( outputPtr->chisqVec );
	 outputPtr->chisqVec = NULL;
       }
       LALFree( outputPtr->detectorVec );
       outputPtr->detectorVec = NULL;
       LALFree( outputPtr );
       *output = NULL;
      }
      ENDFAIL( status );
    }

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


void
LALCoherentInspiralFilterParamsFinalize (
    LALStatus                       *status,
    CoherentInspiralFilterParams   **output
    )
{
  CoherentInspiralFilterParams     *outputPtr;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  /*
   *
   * ensure that the arguments are reasonable
   *
   */


  /* check that the output handle exists, but points to a non-null pointer */
  ASSERT( output, status, COHERENTINSPIRALH_ENULL, COHERENTINSPIRALH_MSGENULL );
  ASSERT( *output, status, COHERENTINSPIRALH_ENULL, COHERENTINSPIRALH_MSGENULL );

  /* local pointer to output structure */
  outputPtr = *output;


  /* destroy detector vector */
  LALFree( outputPtr->detectorVec->detector );
  outputPtr->detectorVec->detector = NULL;
  LALFree( outputPtr->detectorVec );
  outputPtr->detectorVec = NULL;

  LALU2DestroyVector(status->statusPtr, &(outputPtr->detIDVec) );
  CHECKSTATUSPTR( status );

  LALDDestroyVector(status->statusPtr, &(outputPtr->sigmasqVec) );
  CHECKSTATUSPTR( status );

  LALSDestroyVector(status->statusPtr, &(outputPtr->chisqVec) );
  CHECKSTATUSPTR( status );

  /*
   *
   * destroy coherent SNR vector, if it exists
   *
   */
  if ( outputPtr->cohSNRVec3Sites ) {
      LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec3Sites->data) );
      CHECKSTATUSPTR( status );
      LALFree( outputPtr->cohSNRVec3Sites );
      outputPtr->cohSNRVec3Sites = NULL;
  }
  if ( outputPtr->cohSNRVec ) {
      LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec->data) );
      CHECKSTATUSPTR( status );
      LALFree( outputPtr->cohSNRVec );
      outputPtr->cohSNRVec = NULL;
  }

  /*
   *
   * destroy H1-H2 coherent phase difference vector, if it exists
   *
   */
  if ( outputPtr->cohH1H2SNRVec ) {
    LALDestroyVector( status->statusPtr, &(outputPtr->cohH1H2SNRVec->data) );
    CHECKSTATUSPTR( status );
    LALFree( outputPtr->cohH1H2SNRVec );
    outputPtr->cohH1H2SNRVec = NULL;
  }

  /*
   *
   * destroy network null statistic vector, if it exists
   *
   */
  if ( outputPtr->nullStatVec3Sites ) {
      LALDestroyVector( status->statusPtr, &(outputPtr->nullStatVec3Sites->data) );
      CHECKSTATUSPTR( status );
      LALFree( outputPtr->nullStatVec3Sites );
      outputPtr->nullStatVec3Sites = NULL;
  }
  if ( outputPtr->nullStatVec ) {
      LALDestroyVector( status->statusPtr, &(outputPtr->nullStatVec->data) );
      CHECKSTATUSPTR( status );
      LALFree( outputPtr->nullStatVec );
      outputPtr->nullStatVec = NULL;
  }

  /*
   *
   * destroy H1-H2 null statistic vector, if it exists
   *
   */
  if ( outputPtr->nullStatH1H2Vec ) {
    LALDestroyVector( status->statusPtr, &(outputPtr->nullStatH1H2Vec->data) );
    CHECKSTATUSPTR( status );
    LALFree( outputPtr->nullStatH1H2Vec );
    outputPtr->nullStatH1H2Vec = NULL;
  }

  LALFree( outputPtr );
  *output = NULL;


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

void
LALCoherentInspiralEstimatePsiEpsilonCoaPhase (
    LALStatus                             *status,
    INT4                                   caseID[LAL_NUM_IFO],
    REAL8                                 *sigmasq,
    REAL4                                  theta,
    REAL4                                  phi,
    COMPLEX8                               cData[4],
    REAL4                                 *inclination,
    REAL4                                 *polarization,
    REAL4                                 *coaPhase
    )
{
  INT4                   i = 0;
  INT4                   j = 0;
  INT4                   k = 0;

  REAL4                  alphaBetaGamma[6][3]; /* H1 L V G T H2 */
  REAL4                  alphaBetaGammaH[3] = { 38.11, 256.35, 107.43 };
  REAL4                  alphaBetaGammaL[3] = { 38.09, 283.54, 196.88 };
  REAL4                  alphaBetaGammaV[3] = { 320.34, 275.92, 159.02 };
  REAL4                  alphaBetaGammaG[3] = { 326.41, 269.86, 110.4 };
  REAL4                  alphaBetaGammaT[3] = { 16.67, 188.47, 326.3 };
  REAL4                  gelFandABGPlusRe[6][5];
  REAL4                  gelFandABGMinusRe[6][5];
  REAL4                  gelFandPTZPlusRe[5];
  REAL4                  gelFandPTZMinusRe[5];
  REAL4                  gelFandABGPlusIm[6][5];
  REAL4                  gelFandABGMinusIm[6][5];
  REAL4                  gelFandPTZPlusIm[5];
  REAL4                  gelFandPTZMinusIm[5];
  REAL4                  gelFandPEZPlusRe = 0.0;
  REAL4                  gelFandPEZPlusIm = 0.0;
  REAL4                  gelFandPEZMinusRe = 0.0;
  REAL4                  gelFandPEZMinusIm = 0.0;
  REAL4                  mPTZRe[3] = { 0.0, 0.0, 0.0 };
  REAL4                  mPTZIm[3] = { 0.0, 0.0, 0.0 };
  REAL4                  mABGRe[6][3];
  REAL4                  mABGIm[6][3];
  REAL4                  sphereHarmonicRe[5][3][3]; /* -2 -1 0 1 2 */
  REAL4                  sphereHarmonicIm[5][3][3];
  REAL4                  sphereHarmonicReZero[3][3] = {{1,0,0},{0,-1,0},{0,0,0}};
  REAL4                  sphereHarmonicImZero[3][3] = {{0,-1,0},{-1,0,0},{0,0,0}};
  REAL4                  sphereHarmonicReOne[3][3] = {{0,0,1},{0,0,0},{1,0,0}};
  REAL4                  sphereHarmonicImOne[3][3] = {{0,0,0},{0,0,-1},{0,-1,0}};
  REAL4                  sphereHarmonicReTwo[3][3] = {{-1,0,0},{0,-1,0},{0,0,2}};
  REAL4                  sphereHarmonicImTwo[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
  REAL4                  sphereHarmonicReThree[3][3] = {{0,0,1},{0,0,0},{1,0,0}};
  REAL4                  sphereHarmonicImThree[3][3] = {{0,0,0},{0,0,1},{0,1,0}};
  REAL4                  sphereHarmonicReFour[3][3] = {{1,0,0},{0,-1,0},{0,0,0}};
  REAL4                  sphereHarmonicImFour[3][3] = {{0,1,0},{1,0,0},{0,0,0}};
  REAL4                  dVectorPlusRe[6]; /* H1 L V G T H2 */
  REAL4                  dVectorPlusIm[6];
  REAL4                  dVectorMinusRe[6];
  REAL4                  dVectorMinusIm[6];
  REAL4                  eVectorRe[6];
  REAL4                  eVectorIm[6];
  REAL4                  eVectorNorm = 0.0;
  REAL4                  qVectorRe[6];
  REAL4                  qVectorIm[6];
  REAL4                  cDotqRe = 0.0;
  REAL4                  cDotqIm = 0.0;
  COMPLEX8               cPlus;
  COMPLEX8               cMinus;
  COMPLEX8               cRatio;  /* cMinus/cPlus */

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  /* Must have 3 sites to estimate psi and epsilon */
  ASSERT( sigmasq[0] && sigmasq[1] && sigmasq[2], status,
	  COHERENTINSPIRALH_ENUMZ, COHERENTINSPIRALH_MSGENUMZ );
  ASSERT( crealf(cData[0]) && crealf(cData[1]) && crealf(cData[2]), status,
	  COHERENTINSPIRALH_ENUMZ, COHERENTINSPIRALH_MSGENUMZ );

  /* Initialize the arrays */
  for( i = 0; i < 6; i++ )
    {
      for( j = 0; j < 3; j++ )
	{
	  alphaBetaGamma[i][j] = 0.0;
	  mABGRe[i][j] = 0.0;
	  mABGIm[i][j] = 0.0;
	}
    }
  for( i = 0; i < 6; i++ )
    {
      for( j = 0; j < 5; j++ )
	{
	  gelFandABGPlusRe[i][j] = 0.0;
	  gelFandABGMinusRe[i][j] = 0.0;
	  gelFandABGPlusIm[i][j] = 0.0;
	  gelFandABGMinusIm[i][j] = 0.0;
	}
    }
  for( i = 0; i < 5; i++ )
    {
      gelFandPTZPlusRe[i] = 0.0;
      gelFandPTZMinusRe[i] = 0.0;
      gelFandPTZPlusIm[i] = 0.0;
      gelFandPTZMinusIm[i] = 0.0;
    }
  for( i = 0; i < 5; i++ )
    {
      for( j = 0; j < 3; j++ )
	{
	  for( k = 0; k < 3; k++ )
	    {
	      sphereHarmonicRe[i][j][k] = 0.0;
	      sphereHarmonicIm[i][j][k] = 0.0;
	    }
	}
    }
  for( i = 0; i < 6; i++ )
    {
      dVectorPlusRe[i] = 0.0;
      dVectorPlusIm[i] = 0.0;
      dVectorMinusRe[i] = 0.0;
      dVectorMinusIm[i] = 0.0;
      eVectorRe[i] = 0.0;
      eVectorIm[i] = 0.0;
      qVectorRe[i] = 0.0;
      qVectorIm[i] = 0.0;
    }

  /* Set the euler angles for the various sites */
  for( i = 0; i < 3; i++ )
    {
      alphaBetaGamma[0][i] = LAL_PI_180 * alphaBetaGammaH[i];
      alphaBetaGamma[1][i] = LAL_PI_180 * alphaBetaGammaL[i];
      alphaBetaGamma[2][i] = LAL_PI_180 * alphaBetaGammaV[i];
      alphaBetaGamma[3][i] = LAL_PI_180 * alphaBetaGammaG[i];
      alphaBetaGamma[4][i] = LAL_PI_180 * alphaBetaGammaT[i];
      alphaBetaGamma[5][i] = LAL_PI_180 * alphaBetaGammaH[i];
    }

  /* Compute the various Gel'Fand Functions */
  mPTZRe[0] = LAL_SQRT1_2 * cos(phi);
  mPTZRe[1] = LAL_SQRT1_2 * sin(phi);
  mPTZRe[2] = 0;
  mPTZIm[0] = LAL_SQRT1_2 * -cos(theta)*sin(phi);
  mPTZIm[1] = LAL_SQRT1_2 * cos(theta)*cos(phi);
  mPTZIm[2] = LAL_SQRT1_2 * sin(theta);

  for( i = 0; i < 6; i++ )
    {
      mABGRe[i][0] = LAL_SQRT1_2 * (cos(alphaBetaGamma[i][0]) *
            cos(alphaBetaGamma[i][2]) - sin(alphaBetaGamma[i][0])*
            cos(alphaBetaGamma[i][1])*sin(alphaBetaGamma[i][2]));
      mABGRe[i][1] = LAL_SQRT1_2 * (sin(alphaBetaGamma[i][0]) *
            cos(alphaBetaGamma[i][2]) + cos(alphaBetaGamma[i][0])*
            cos(alphaBetaGamma[i][1])*sin(alphaBetaGamma[i][2]));
      mABGRe[i][2] = LAL_SQRT1_2 * sin(alphaBetaGamma[i][1])*
            sin(alphaBetaGamma[i][2]);
      mABGIm[i][0] = LAL_SQRT1_2 * -(cos(alphaBetaGamma[i][0]) *
            sin(alphaBetaGamma[i][2]) + sin(alphaBetaGamma[i][0])*
            cos(alphaBetaGamma[i][1])*cos(alphaBetaGamma[i][2]));
      mABGIm[i][1] = LAL_SQRT1_2 * (-sin(alphaBetaGamma[i][0]) *
            sin(alphaBetaGamma[i][2]) + cos(alphaBetaGamma[i][0])*
            cos(alphaBetaGamma[i][1])*cos(alphaBetaGamma[i][2]));
      mABGIm[i][2] = LAL_SQRT1_2 * sin(alphaBetaGamma[i][1])*
            cos(alphaBetaGamma[i][2]);

    }

  for( i = 0; i < 3; i++ )
    {
      for( j = 0;j < 3; j++ )
	{
	  sphereHarmonicRe[0][i][j] = 0.25*sqrt(15/(2*LAL_PI)) * sphereHarmonicReZero[i][j];
	  sphereHarmonicIm[0][i][j] = 0.25*sqrt(15/(2*LAL_PI)) * sphereHarmonicImZero[i][j];
	  sphereHarmonicRe[1][i][j] = 0.25*sqrt(15/(2*LAL_PI)) * sphereHarmonicReOne[i][j];
	  sphereHarmonicIm[1][i][j] = 0.25*sqrt(15/(2*LAL_PI)) * sphereHarmonicImOne[i][j];
	  sphereHarmonicRe[2][i][j] = 0.5*sqrt(5/(4*LAL_PI)) * sphereHarmonicReTwo[i][j];
	  sphereHarmonicIm[2][i][j] = 0.5*sqrt(5/(4*LAL_PI)) * sphereHarmonicImTwo[i][j];
	  sphereHarmonicRe[3][i][j] = -0.25*sqrt(15/(2*LAL_PI)) * sphereHarmonicReThree[i][j];
	  sphereHarmonicIm[3][i][j] = -0.25*sqrt(15/(2*LAL_PI)) * sphereHarmonicImThree[i][j];
	  sphereHarmonicRe[4][i][j] = 0.25*sqrt(15/(2*LAL_PI)) * sphereHarmonicReFour[i][j];
	  sphereHarmonicIm[4][i][j] = 0.25*sqrt(15/(2*LAL_PI)) * sphereHarmonicImFour[i][j];
	}
    }

  for( i = 0; i < 3; i++ )
    {
      for( j = 0; j < 3; j++ )
	{
	  gelFandPTZPlusRe[0] += sphereHarmonicRe[4][i][j]*mPTZRe[i]*mPTZRe[j]
                - sphereHarmonicRe[4][i][j]*mPTZIm[i]*mPTZIm[j]
                - sphereHarmonicIm[4][i][j]*mPTZRe[i]*mPTZIm[j]
                - sphereHarmonicIm[4][i][j]*mPTZIm[i]*mPTZRe[j];

	  gelFandPTZPlusIm[0] += sphereHarmonicRe[4][i][j]*mPTZRe[i]*mPTZIm[j]
                + sphereHarmonicRe[4][i][j]*mPTZIm[i]*mPTZRe[j]
                + sphereHarmonicIm[4][i][j]*mPTZRe[i]*mPTZRe[j]
                - sphereHarmonicIm[4][i][j]*mPTZIm[i]*mPTZIm[j];

	  gelFandPTZPlusRe[1] += -sphereHarmonicRe[3][i][j]*mPTZRe[i]*mPTZRe[j]
                + sphereHarmonicRe[3][i][j]*mPTZIm[i]*mPTZIm[j]
                + sphereHarmonicIm[3][i][j]*mPTZRe[i]*mPTZIm[j]
                + sphereHarmonicIm[3][i][j]*mPTZIm[i]*mPTZRe[j];

	  gelFandPTZPlusIm[1] += -sphereHarmonicRe[3][i][j]*mPTZRe[i]*mPTZIm[j]
                - sphereHarmonicRe[3][i][j]*mPTZIm[i]*mPTZRe[j]
                - sphereHarmonicIm[3][i][j]*mPTZRe[i]*mPTZRe[j]
                + sphereHarmonicIm[3][i][j]*mPTZIm[i]*mPTZIm[j];

	  gelFandPTZPlusRe[2] += sphereHarmonicRe[2][i][j]*mPTZRe[i]*mPTZRe[j]
                - sphereHarmonicRe[2][i][j]*mPTZIm[i]*mPTZIm[j]
                - sphereHarmonicIm[2][i][j]*mPTZRe[i]*mPTZIm[j]
                - sphereHarmonicIm[2][i][j]*mPTZIm[i]*mPTZRe[j];

	  gelFandPTZPlusIm[2] += sphereHarmonicRe[2][i][j]*mPTZRe[i]*mPTZIm[j]
                + sphereHarmonicRe[2][i][j]*mPTZIm[i]*mPTZRe[j]
                + sphereHarmonicIm[2][i][j]*mPTZRe[i]*mPTZRe[j]
                - sphereHarmonicIm[2][i][j]*mPTZIm[i]*mPTZIm[j];

	  gelFandPTZPlusRe[3] += -sphereHarmonicRe[1][i][j]*mPTZRe[i]*mPTZRe[j]
                + sphereHarmonicRe[1][i][j]*mPTZIm[i]*mPTZIm[j]
                + sphereHarmonicIm[1][i][j]*mPTZRe[i]*mPTZIm[j]
                + sphereHarmonicIm[1][i][j]*mPTZIm[i]*mPTZRe[j];

	  gelFandPTZPlusIm[3] += -sphereHarmonicRe[1][i][j]*mPTZRe[i]*mPTZIm[j]
                - sphereHarmonicRe[1][i][j]*mPTZIm[i]*mPTZRe[j]
                - sphereHarmonicIm[1][i][j]*mPTZRe[i]*mPTZRe[j]
                + sphereHarmonicIm[1][i][j]*mPTZIm[i]*mPTZIm[j];

	  gelFandPTZPlusRe[4] += sphereHarmonicRe[0][i][j]*mPTZRe[i]*mPTZRe[j]
                - sphereHarmonicRe[0][i][j]*mPTZIm[i]*mPTZIm[j]
                - sphereHarmonicIm[0][i][j]*mPTZRe[i]*mPTZIm[j]
                - sphereHarmonicIm[0][i][j]*mPTZIm[i]*mPTZRe[j];

	  gelFandPTZPlusIm[4] += sphereHarmonicRe[0][i][j]*mPTZRe[i]*mPTZIm[j]
                + sphereHarmonicRe[0][i][j]*mPTZIm[i]*mPTZRe[j]
                + sphereHarmonicIm[0][i][j]*mPTZRe[i]*mPTZRe[j]
                - sphereHarmonicIm[0][i][j]*mPTZIm[i]*mPTZIm[j];

	  gelFandPTZMinusRe[0] += sphereHarmonicRe[4][i][j]*mPTZRe[i]*mPTZRe[j]
                - sphereHarmonicRe[4][i][j]*mPTZIm[i]*mPTZIm[j]
                + sphereHarmonicIm[4][i][j]*mPTZRe[i]*mPTZIm[j]
                + sphereHarmonicIm[4][i][j]*mPTZIm[i]*mPTZRe[j];

	  gelFandPTZMinusIm[0] += -sphereHarmonicRe[4][i][j]*mPTZRe[i]*mPTZIm[j]
                - sphereHarmonicRe[4][i][j]*mPTZIm[i]*mPTZRe[j]
                + sphereHarmonicIm[4][i][j]*mPTZRe[i]*mPTZRe[j]
                - sphereHarmonicIm[4][i][j]*mPTZIm[i]*mPTZIm[j];

	  gelFandPTZMinusRe[1] += -sphereHarmonicRe[3][i][j]*mPTZRe[i]*mPTZRe[j]
                + sphereHarmonicRe[3][i][j]*mPTZIm[i]*mPTZIm[j]
                - sphereHarmonicIm[3][i][j]*mPTZRe[i]*mPTZIm[j]
                - sphereHarmonicIm[3][i][j]*mPTZIm[i]*mPTZRe[j];

	  gelFandPTZMinusIm[1] += sphereHarmonicRe[3][i][j]*mPTZRe[i]*mPTZIm[j]
                + sphereHarmonicRe[3][i][j]*mPTZIm[i]*mPTZRe[j]
                - sphereHarmonicIm[3][i][j]*mPTZRe[i]*mPTZRe[j]
                + sphereHarmonicIm[3][i][j]*mPTZIm[i]*mPTZIm[j];

	  gelFandPTZMinusRe[2] += sphereHarmonicRe[2][i][j]*mPTZRe[i]*mPTZRe[j]
                - sphereHarmonicRe[2][i][j]*mPTZIm[i]*mPTZIm[j]
                + sphereHarmonicIm[2][i][j]*mPTZRe[i]*mPTZIm[j]
                + sphereHarmonicIm[2][i][j]*mPTZIm[i]*mPTZRe[j];

	  gelFandPTZMinusIm[2] += -sphereHarmonicRe[2][i][j]*mPTZRe[i]*mPTZIm[j]
                - sphereHarmonicRe[2][i][j]*mPTZIm[i]*mPTZRe[j]
                + sphereHarmonicIm[2][i][j]*mPTZRe[i]*mPTZRe[j]
                - sphereHarmonicIm[2][i][j]*mPTZIm[i]*mPTZIm[j];

	  gelFandPTZMinusRe[3] += -sphereHarmonicRe[1][i][j]*mPTZRe[i]*mPTZRe[j]
                + sphereHarmonicRe[1][i][j]*mPTZIm[i]*mPTZIm[j]
                - sphereHarmonicIm[1][i][j]*mPTZRe[i]*mPTZIm[j]
                - sphereHarmonicIm[1][i][j]*mPTZIm[i]*mPTZRe[j];

	  gelFandPTZMinusIm[3] += sphereHarmonicRe[1][i][j]*mPTZRe[i]*mPTZIm[j]
                + sphereHarmonicRe[1][i][j]*mPTZIm[i]*mPTZRe[j]
                - sphereHarmonicIm[1][i][j]*mPTZRe[i]*mPTZRe[j]
                + sphereHarmonicIm[1][i][j]*mPTZIm[i]*mPTZIm[j];

	  gelFandPTZMinusRe[4] += sphereHarmonicRe[0][i][j]*mPTZRe[i]*mPTZRe[j]
                - sphereHarmonicRe[0][i][j]*mPTZIm[i]*mPTZIm[j]
                + sphereHarmonicIm[0][i][j]*mPTZRe[i]*mPTZIm[j]
                + sphereHarmonicIm[0][i][j]*mPTZIm[i]*mPTZRe[j];

	  gelFandPTZMinusIm[4] += -sphereHarmonicRe[0][i][j]*mPTZRe[i]*mPTZIm[j]
                - sphereHarmonicRe[0][i][j]*mPTZIm[i]*mPTZRe[j]
                + sphereHarmonicIm[0][i][j]*mPTZRe[i]*mPTZRe[j]
                - sphereHarmonicIm[0][i][j]*mPTZIm[i]*mPTZIm[j];
	}
    }

  for( k = 0; k < 5; k++)
    {
      for( i = 0; i < 3; i++ )
	{
	  for( j = 0; j < 3; j++ )
	    {
	      gelFandABGPlusRe[k][0] += sphereHarmonicRe[4][i][j]*mABGRe[k][i]*mABGRe[k][j]
                    - sphereHarmonicRe[4][i][j]*mABGIm[k][i]*mABGIm[k][j]
                    - sphereHarmonicIm[4][i][j]*mABGRe[k][i]*mABGIm[k][j]
	            - sphereHarmonicIm[4][i][j]*mABGIm[k][i]*mABGRe[k][j];

	      gelFandABGPlusIm[k][0] += sphereHarmonicRe[4][i][j]*mABGRe[k][i]*mABGIm[k][j]
                    + sphereHarmonicRe[4][i][j]*mABGIm[k][i]*mABGRe[k][j]
                    + sphereHarmonicIm[4][i][j]*mABGRe[k][i]*mABGRe[k][j]
                    - sphereHarmonicIm[4][i][j]*mABGIm[k][i]*mABGIm[k][j];

	      gelFandABGPlusRe[k][1] += -sphereHarmonicRe[3][i][j]*mABGRe[k][i]*mABGRe[k][j]
                    + sphereHarmonicRe[3][i][j]*mABGIm[k][i]*mABGIm[k][j]
                    + sphereHarmonicIm[3][i][j]*mABGRe[k][i]*mABGIm[k][j]
                    + sphereHarmonicIm[3][i][j]*mABGIm[k][i]*mABGRe[k][j];

	      gelFandABGPlusIm[k][1] += -sphereHarmonicRe[3][i][j]*mABGRe[k][i]*mABGIm[k][j]
                    - sphereHarmonicRe[3][i][j]*mABGIm[k][i]*mABGRe[k][j]
                    - sphereHarmonicIm[3][i][j]*mABGRe[k][i]*mABGRe[k][j]
                    + sphereHarmonicIm[3][i][j]*mABGIm[k][i]*mABGIm[k][j];

	      gelFandABGPlusRe[k][2] += sphereHarmonicRe[2][i][j]*mABGRe[k][i]*mABGRe[k][j]
                    - sphereHarmonicRe[2][i][j]*mABGIm[k][i]*mABGIm[k][j]
                    - sphereHarmonicIm[2][i][j]*mABGRe[k][i]*mABGIm[k][j]
                    - sphereHarmonicIm[2][i][j]*mABGIm[k][i]*mABGRe[k][j];

	      gelFandABGPlusIm[k][2] += sphereHarmonicRe[2][i][j]*mABGRe[k][i]*mABGIm[k][j]
                    + sphereHarmonicRe[2][i][j]*mABGIm[k][i]*mABGRe[k][j]
                    + sphereHarmonicIm[2][i][j]*mABGRe[k][i]*mABGRe[k][j]
                    - sphereHarmonicIm[2][i][j]*mABGIm[k][i]*mABGIm[k][j];

	      gelFandABGPlusRe[k][3] += -sphereHarmonicRe[1][i][j]*mABGRe[k][i]*mABGRe[k][j]
                    + sphereHarmonicRe[1][i][j]*mABGIm[k][i]*mABGIm[k][j]
                    + sphereHarmonicIm[1][i][j]*mABGRe[k][i]*mABGIm[k][j]
                    + sphereHarmonicIm[1][i][j]*mABGIm[k][i]*mABGRe[k][j];

	      gelFandABGPlusIm[k][3] += -sphereHarmonicRe[1][i][j]*mABGRe[k][i]*mABGIm[k][j]
                    - sphereHarmonicRe[1][i][j]*mABGIm[k][i]*mABGRe[k][j]
                    - sphereHarmonicIm[1][i][j]*mABGRe[k][i]*mABGRe[k][j]
                    + sphereHarmonicIm[1][i][j]*mABGIm[k][i]*mABGIm[k][j];

	      gelFandABGPlusRe[k][4] += sphereHarmonicRe[0][i][j]*mABGRe[k][i]*mABGRe[k][j]
                    - sphereHarmonicRe[0][i][j]*mABGIm[k][i]*mABGIm[k][j]
                    - sphereHarmonicIm[0][i][j]*mABGRe[k][i]*mABGIm[k][j]
                    - sphereHarmonicIm[0][i][j]*mABGIm[k][i]*mABGRe[k][j];

	      gelFandABGPlusIm[k][4] += sphereHarmonicRe[0][i][j]*mABGRe[k][i]*mABGIm[k][j]
                    + sphereHarmonicRe[0][i][j]*mABGIm[k][i]*mABGRe[k][j]
                    + sphereHarmonicIm[0][i][j]*mABGRe[k][i]*mABGRe[k][j]
                    - sphereHarmonicIm[0][i][j]*mABGIm[k][i]*mABGIm[k][j];

	      gelFandABGMinusRe[k][0] += sphereHarmonicRe[4][i][j]*mABGRe[k][i]*mABGRe[k][j]
                    - sphereHarmonicRe[4][i][j]*mABGIm[k][i]*mABGIm[k][j]
                    + sphereHarmonicIm[4][i][j]*mABGRe[k][i]*mABGIm[k][j]
                    + sphereHarmonicIm[4][i][j]*mABGIm[k][i]*mABGRe[k][j];

	      gelFandABGMinusIm[k][0] += -sphereHarmonicRe[4][i][j]*mABGRe[k][i]*mABGIm[k][j]
                    - sphereHarmonicRe[4][i][j]*mABGIm[k][i]*mABGRe[k][j]
                    + sphereHarmonicIm[4][i][j]*mABGRe[k][i]*mABGRe[k][j]
                    - sphereHarmonicIm[4][i][j]*mABGIm[k][i]*mABGIm[k][j];

	      gelFandABGMinusRe[k][1] += -sphereHarmonicRe[3][i][j]*mABGRe[k][i]*mABGRe[k][j]
                    + sphereHarmonicRe[3][i][j]*mABGIm[k][i]*mABGIm[k][j]
                    - sphereHarmonicIm[3][i][j]*mABGRe[k][i]*mABGIm[k][j]
                    - sphereHarmonicIm[3][i][j]*mABGIm[k][i]*mABGRe[k][j];

	      gelFandABGMinusIm[k][1] += sphereHarmonicRe[3][i][j]*mABGRe[k][i]*mABGIm[k][j]
                    + sphereHarmonicRe[3][i][j]*mABGIm[k][i]*mABGRe[k][j]
                    - sphereHarmonicIm[3][i][j]*mABGRe[k][i]*mABGRe[k][j]
                    + sphereHarmonicIm[3][i][j]*mABGIm[k][i]*mABGIm[k][j];

	      gelFandABGMinusRe[k][2] += sphereHarmonicRe[2][i][j]*mABGRe[k][i]*mABGRe[k][j]
                    - sphereHarmonicRe[2][i][j]*mABGIm[k][i]*mABGIm[k][j]
                    + sphereHarmonicIm[2][i][j]*mABGRe[k][i]*mABGIm[k][j]
                    + sphereHarmonicIm[2][i][j]*mABGIm[k][i]*mABGRe[k][j];

	      gelFandABGMinusIm[k][2] += -sphereHarmonicRe[2][i][j]*mABGRe[k][i]*mABGIm[k][j]
                    - sphereHarmonicRe[2][i][j]*mABGIm[k][i]*mABGRe[k][j]
                    + sphereHarmonicIm[2][i][j]*mABGRe[k][i]*mABGRe[k][j]
                    - sphereHarmonicIm[2][i][j]*mABGIm[k][i]*mABGIm[k][j];

	      gelFandABGMinusRe[k][3] += -sphereHarmonicRe[1][i][j]*mABGRe[k][i]*mABGRe[k][j]
                    + sphereHarmonicRe[1][i][j]*mABGIm[k][i]*mABGIm[k][j]
                    - sphereHarmonicIm[1][i][j]*mABGRe[k][i]*mABGIm[k][j]
                    - sphereHarmonicIm[1][i][j]*mABGIm[k][i]*mABGRe[k][j];

	      gelFandABGMinusIm[k][3] += sphereHarmonicRe[1][i][j]*mABGRe[k][i]*mABGIm[k][j]
                    + sphereHarmonicRe[1][i][j]*mABGIm[k][i]*mABGRe[k][j]
                    - sphereHarmonicIm[1][i][j]*mABGRe[k][i]*mABGRe[k][j]
                    + sphereHarmonicIm[1][i][j]*mABGIm[k][i]*mABGIm[k][j];

	      gelFandABGMinusRe[k][4] += sphereHarmonicRe[0][i][j]*mABGRe[k][i]*mABGRe[k][j]
                    - sphereHarmonicRe[0][i][j]*mABGIm[k][i]*mABGIm[k][j]
                    + sphereHarmonicIm[0][i][j]*mABGRe[k][i]*mABGIm[k][j]
                    + sphereHarmonicIm[0][i][j]*mABGIm[k][i]*mABGRe[k][j];

	      gelFandABGMinusIm[k][4] += -sphereHarmonicRe[0][i][j]*mABGRe[k][i]*mABGIm[k][j]
                    - sphereHarmonicRe[0][i][j]*mABGIm[k][i]*mABGRe[k][j]
                    + sphereHarmonicIm[0][i][j]*mABGRe[k][i]*mABGRe[k][j]
                    - sphereHarmonicIm[0][i][j]*mABGIm[k][i]*mABGIm[k][j];
	    }
	}
    }

  /* Dont forget the coefficients... */

  for( i = 0; i < 5; i++ )
    {
      gelFandPTZPlusRe[i] *= sqrt(8*LAL_PI/15);
      gelFandPTZPlusIm[i] *= sqrt(8*LAL_PI/15);
      gelFandPTZMinusRe[i] *= sqrt(8*LAL_PI/15);
      gelFandPTZMinusIm[i] *= sqrt(8*LAL_PI/15);
    }


  for( i = 0; i < 6; i++ )
    {
      for( j = 0; j < 5; j++ )
	{
	  gelFandABGPlusRe[i][j] *= sqrt(8*LAL_PI/15);
	  gelFandABGPlusIm[i][j] *= sqrt(8*LAL_PI/15);
	  gelFandABGMinusRe[i][j] *= sqrt(8*LAL_PI/15);
	  gelFandABGMinusIm[i][j] *= sqrt(8*LAL_PI/15);
	}
    }

  /* Armed with the gelFand functions, the dVectors can be constructed */

  for( i = 0; i < 6; i++ )
    {
      for( j = 0; j < 5; j++ )
	{
	  dVectorPlusRe[i] += gelFandPTZPlusIm[j] * ( gelFandABGPlusRe[i][j]
                 - gelFandABGMinusRe[i][j] ) + gelFandPTZPlusRe[j] *
                 ( gelFandABGMinusIm[i][j]  - gelFandABGPlusIm[i][j] );
	  dVectorMinusRe[i] += gelFandPTZMinusIm[j] * ( gelFandABGPlusRe[i][j]
                 - gelFandABGMinusRe[i][j] ) + gelFandPTZMinusRe[j] *
                 ( gelFandABGMinusIm[i][j]  - gelFandABGPlusIm[i][j] );
	  dVectorPlusIm[i] += gelFandPTZPlusIm[j] * ( gelFandABGMinusIm[i][j]
                 - gelFandABGPlusIm[i][j] ) - gelFandPTZPlusRe[j] *
                 ( gelFandABGPlusRe[i][j] - gelFandABGMinusRe[i][j] );
	  dVectorMinusIm[i] += gelFandPTZMinusIm[j] * ( gelFandABGMinusIm[i][j]
                 - gelFandABGPlusIm[i][j] ) - gelFandPTZMinusRe[j] *
                 ( gelFandABGPlusRe[i][j] - gelFandABGMinusRe[i][j] );
	}
    }

  /* Now form the C's ( +2 and -2 ) */

  cPlus = 0.0;
  cMinus = 0.0;
  cRatio = 0.0;

  k = 0;
  for( i = 0; i < 6; i++ )
    {
      if( caseID[i] )
	{
	  cPlus += crectf( sigmasq[k] * ( crealf(cData[k]) * dVectorPlusRe[i] - cimagf(cData[k]) * dVectorPlusIm[i] ), sigmasq[k] * ( crealf(cData[k]) * dVectorPlusIm[i] + cimagf(cData[k]) * dVectorPlusRe[i] ) );
	  cMinus += crectf( sigmasq[k] * ( crealf(cData[k]) * dVectorMinusRe[i] - cimagf(cData[k]) * dVectorMinusIm[i] ), sigmasq[k] * ( crealf(cData[k]) * dVectorMinusIm[i] + cimagf(cData[k]) * dVectorMinusRe[i] ) );
	  k++;

	}
    }

  cRatio = crectf( ( crealf(cPlus) * crealf(cMinus) + cimagf(cPlus) * cimagf(cMinus) ) / ( crealf(cPlus) * crealf(cPlus) + cimagf(cPlus) * cimagf(cPlus) ), ( crealf(cPlus) * cimagf(cMinus) - cimagf(cPlus) * crealf(cMinus) ) / ( crealf(cPlus) * crealf(cPlus) + cimagf(cPlus) * cimagf(cPlus) ) );

  /* Now the estimates can be computed */

  *inclination = acos(  ( 1 - sqrt( sqrt(crealf(cRatio)*crealf(cRatio) + cimagf(cRatio)*cimagf(cRatio)) ))/
	  ( 1 + sqrt( sqrt(crealf(cRatio)*crealf(cRatio) + cimagf(cRatio)*cimagf(cRatio)) )) );
  *polarization = 0.25 * atan( cimagf(cRatio) / crealf(cRatio) );

  /* To estimate the coaphase, I need to construct the eVectors */

  gelFandPEZMinusRe = 0.25 * (cos(*inclination) - 1)*(cos(*inclination) - 1) * cos(2 * (*polarization));
  gelFandPEZMinusIm = 0.25 * (cos(*inclination) - 1)*(cos(*inclination) - 1) * sin(2 * (*polarization));
  gelFandPEZPlusRe = 0.25 * (cos(*inclination) + 1)*(cos(*inclination) + 1) * cos(2 * (*polarization));
  gelFandPEZPlusIm = -0.25 * (cos(*inclination) + 1)*(cos(*inclination) + 1) * sin(2 * (*polarization));

  k = 0;
  for( i = 0; i < 6; i++ )
    {
      if( caseID[i] )
	{
	  eVectorRe[i] = sigmasq[k] * (gelFandPEZPlusRe*dVectorPlusRe[i] - gelFandPEZPlusIm*dVectorPlusIm[i] + gelFandPEZMinusRe*dVectorPlusRe[i] + gelFandPEZMinusIm*dVectorPlusIm[i]);
	  eVectorIm[i] = sigmasq[k] * (gelFandPEZPlusRe*dVectorPlusIm[i] + gelFandPEZPlusIm*dVectorPlusRe[i] + gelFandPEZMinusIm*dVectorPlusRe[i] - gelFandPEZMinusRe*dVectorPlusIm[i]);
	  k++;
	}
    }

  /* Now form the Q's */

  for( i = 0; i < 6; i++ )
    {
      eVectorNorm += eVectorRe[i]*eVectorRe[i] + eVectorIm[i]*eVectorIm[i];
    }

  eVectorNorm = sqrt( eVectorNorm );

  for( i = 0; i < 6; i++ )
    {
      qVectorRe[i] = eVectorRe[i] / eVectorNorm;
      qVectorIm[i] = eVectorIm[i] / eVectorNorm;
    }

  /* Now for the phase estimate */

  k = 0;
  for( i = 0; i < 6; i++ )
    {
      if( caseID[i] )
	{
	  cDotqRe += crealf(cData[k]) * qVectorRe[i] + cimagf(cData[k]) * qVectorIm[i];
	  cDotqIm += crealf(cData[k]) * qVectorIm[i] - cimagf(cData[k]) * qVectorRe[i];
	  k++;
	}
    }

  *coaPhase = atan( cDotqIm / cDotqRe );

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

void
LALCoherentInspiralEstimateDistance (
    LALStatus                             *status,
    REAL8                                 *sigmasq,
    REAL4                                  UNUSED templateNorm,
    REAL8                                  deltaT,
    INT4                                   segmentLength,  /* time pts */
    REAL4                                  coherentSNR,
    REAL4                                 *distance
    )
{
  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  /* This function computes the effective distance for the H1-H2 pair */
  *distance = sqrt( 0.5 * (sigmasq[0] + sigmasq[1]) * deltaT / segmentLength) / coherentSNR;

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

void
XLALCoherentInspiralFilterSegment (
    LALStatus                             *status,
    MultiInspiralTable                    **eventList,
    CoherentInspiralFilterInput           *input,
    CoherentInspiralFilterParams          *params,
    const SkyGrid                         *skyGridPtr,
    REAL4                                 chisq[4],
    REAL4                                 chisq_dof[4],
    REAL4                                 eff_snr_denom_fac,
    REAL4                                 nullStatRegul
    )
{
  UINT4                               detId = 0;
  UINT4                               cohSNROut = 0;
  UINT4                               nullStatOut = 0;
  UINT4                               nullStatH1H2Out = 0;
  UINT4                               UNUSED case2a = 0;
  UINT4                               UNUSED case2b = 0;
  UINT4                               UNUSED case3a = 0;
  UINT4                               case3b = 0;
  UINT4                               case4a = 0;
  INT4                                detIdSlidTimePt = 0;
  INT4                                caseID[6] = {0,0,0,0,0,0};
  INT4                                i,q,w,m,j,k,l;
  INT4                                indexarray[4] = {0,0,0,0};
  INT4                                numPoints = 0;
  INT4                                deltaEventIndex = 0;
  INT4                                eventStartIdx = 0;
  INT4                                slidePoints[3] = {0,0,0};
  INT4                                slidePoints4D[4] = {0,0,0,0};
  INT4                                sortedSlidePoints3D[3]= {0,0,0};
  INT4                                sortedSlidePoints4D[4]= {0,0,0,0};
  int                                 locIdx;
  REAL4                               buffer = 0.0; /*account for timing errors*/
  REAL4                               timingError = 0.0; /*0.00025;*/  /* allowed timing error of 2 ms */
  REAL4                               s[3]={0,0,0};/* up to 4 distances;in 3D space*/
  REAL4                               distance[4] = {0,0,0,0};
  REAL4                               chirpTime = 0.0;
  REAL4                               cohSNRThresh = 0.0;
  REAL4                               UNUSED cohSNRThreshSq = 0.0;
  double                              inclination = 0.0;
  double                              polarization = 0.0;
  /*REAL4                               distanceEstimate = 0.0;*/
  double                              coaPhase      = 0.0;
  REAL4                               cohSNR        = 0.0;
  REAL4                               cohSNRLocal   = 0.0;
  REAL4                               cohSNRLocalRe = 0.0;
  REAL4                               cohSNRLocalIm = 0.0;
  REAL4                               cohSnrSqLocal = 0.0;
  REAL4                               nullStatRe    = 0.0;
  REAL4                               nullStatIm    = 0.0;
  REAL4                               nullNorm      = 0.0;
  REAL4                               nullStatistic = 0.0;
  REAL4                               cohSnrRe      = 0.0;
  REAL4                               cohSnrIm      = 0.0;
#if 0
  REAL4                               nullNumerRe   = 0.0;
  REAL4                               nullNumerIm   = 0.0;
  REAL4                               chisq4DArray[4]={1.0,1.0,1.0,1.0};
#endif
  REAL8                              *sigmasq = NULL;
  REAL8                               sigmasq4DArray[4]={0.0,0.0,0.0,0.0};
  REAL8                               deltaT = 0.0;
  REAL8                               tempTime = 0.0;
  REAL8                               fracpart = 0.0;
  REAL8                               intpart = 0.0;
  double                              theta = 0.0;
  double                              phi = 0.0;
  double                              timeDelay[4]= {0.0,0.0,0.0,0.0};
  double                              sortedDelays3D[3]= {0.0,0.0,0.0};
  double                              sortedDelays4D[4]= {0.0,0.0,0.0,0.0};
  double                              sortedSigmasq3D[3]= {0.0,0.0,0.0};
  double                              sortedSigmasq4D[4]= {0.0,0.0,0.0,0.0};
  double                              fplus[4] = {0.0,0.0,0.0,0.0};
  double                              fcross[4] = {0.0,0.0,0.0,0.0};
  COMPLEX8                            quadTemp[6];
  LALDetector                         detectors[4];
  COMPLEX8TimeSeries                 *cData[4] = {NULL,NULL,NULL,NULL};
  MultiInspiralTable                 *thisEvent = NULL;
  CHAR                                idtag[6][3] = {"G1","H1","H2","L1","T1","V1"};
  CHAR                                caseStr[FILENAME_MAX];

  UINT4          UNUSED degenerateStat = 0;
  INT4           timePt[4] = {0,0,0,0};
  REAL4          AA=0.0;
  REAL4          BB=0.0;
  REAL4          CC=0.0;
  REAL4          zetaRe=0.0;
  REAL4          zetaIm=0.0;
  REAL4          zetaAbsSq=0.0;
  REAL4          zetaAbsSqMin=1.0;
  REAL4          discrimSqrt=0.0;
  REAL4          VVPlus[4]={0.0,0.0,0.0,0.0};
  REAL4          VVMinus[4]={0.0,0.0,0.0,0.0};
  REAL4          AAn[4]={0.0,0.0,0.0,0.0};
  REAL4          BBn[4]={0.0,0.0,0.0,0.0};
  REAL4          CCn[4]={0.0,0.0,0.0,0.0};
  REAL8          MM1=0.0;
  REAL8          MM2=0.0;
  REAL4          O11=0.0;
  REAL4          O12=0.0;
  REAL4          O21=0.0;
  REAL4          O22=0.0;
  REAL8          gmstInRadians=0.0;
  double         amplitudeConst=1.0,chirpMass=1.0;
  double         UNUSED DRe[4]={0.0,0.0,0.0,0.0};
  double         UNUSED DIm[4]={0.0,0.0,0.0,0.0};
  double         InvMMAA = 0.0, InvMMBB = 0.0, InvMMCC = 0.0;
  double         determinantMM=1.0;
  double         uSigma[4]={0.0,0.0,0.0,0.0};
  double         vSigma[4]={0.0,0.0,0.0,0.0};
  double         aa[4]={0.0,0.0,0.0,0.0};
  double         NN[4]={0.0,0.0,0.0,0.0};

  /* Variables for alternative statistics */
  REAL4 autoCorrCohSq=0.0;
  REAL4 autoCorrNullSq=0.0;
  REAL4 crossCorrCohSq=0.0;
  REAL4 crossCorrNullSq=0.0;
  REAL4 ratioStatLocal=0.0;
  REAL4 eff_coh_snr = 0.0;
  REAL4 crossCorrCohSqEff = 0.0;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* make sure the output handle exists, but points to a null pointer */
  ASSERT( eventList, status, COHERENTINSPIRALH_ENULL, COHERENTINSPIRALH_MSGENULL );
  ASSERT( !*eventList, status, COHERENTINSPIRALH_ENNUL, COHERENTINSPIRALH_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT( params, status, COHERENTINSPIRALH_ENULL, COHERENTINSPIRALH_MSGENULL );

  /* check that the filter parameters are reasonable */
  ASSERT( params->cohSNRThresh > 0, status,
	  COHERENTINSPIRALH_ERHOT, COHERENTINSPIRALH_MSGERHOT );

 /* check that the number of detectors is greater than 1 */
  ASSERT( params->numDetectors > 1, status,
      COHERENTINSPIRALH_EZDET, COHERENTINSPIRALH_MSGEZDET );

  /* check that the number of detectors is less than 5 */
  ASSERT( params->numDetectors < 5, status,
	  COHERENTINSPIRALH_EZDET, COHERENTINSPIRALH_MSGEZDET );

  /* check that the number of segments in positive */
  ASSERT( params->numSegments > 0, status,
      COHERENTINSPIRALH_ESEGZ, COHERENTINSPIRALH_MSGESEGZ );

  /* if a cohSNRVec vector has been created, check we can store data in it  */
  if ( params->cohSNROut ) {
    if ( params->threeSiteCase ) {
      ASSERT( params->cohSNRVec3Sites->data->data, status,
	    COHERENTINSPIRALH_ENULL, COHERENTINSPIRALH_MSGENULL );
      ASSERT( params->cohSNRVec3Sites->data, status,
	    COHERENTINSPIRALH_ENULL, COHERENTINSPIRALH_MSGENULL );
    }
    else {
      ASSERT( params->cohSNRVec->data->data, status,
            COHERENTINSPIRALH_ENULL, COHERENTINSPIRALH_MSGENULL );
      ASSERT( params->cohSNRVec->data, status,
            COHERENTINSPIRALH_ENULL, COHERENTINSPIRALH_MSGENULL );
    }
  }

  /* make sure that the input structure contains some input */
  ASSERT( input->tmplt, status,
	  COHERENTINSPIRALH_ENULL, COHERENTINSPIRALH_MSGENULL );


  /*** read in the parameters ***/
  numPoints = params->numPoints;
  cohSNRThreshSq = params->cohSNRThresh * params->cohSNRThresh;
  cohSNRThresh  = params->cohSNRThresh;
  cohSNROut = params->cohSNROut;
  nullStatOut = params->nullStatOut;
  nullStatH1H2Out = params->nullStatH1H2Out;
  deltaT = params->deltaT;
  /*chirpTime in seconds, just like deltaT*/
  chirpTime = params->chirpTime;
  deltaEventIndex = (UINT4) rint( (chirpTime / deltaT) + 1.0 );
  buffer = rint( (timingError/deltaT) + 1.0 );
  chirpMass = pow(input->tmplt->eta,3.0/5.0)*input->tmplt->totalMass;

  /* Prefactor (in Mpc) for distance estimation */
  amplitudeConst = 4 * pow(LAL_PI * params->fLow,2/3) * pow(LAL_MTSUN_SI*chirpMass,5/12)
    * pow(params->fLow,2/3) * pow(chirpTime/5,0.25);

  /* if the full coherent snr / null vector is required, set it to zero */
  if ( cohSNROut ) {
    if ( params->threeSiteCase ) {
      memset( params->cohSNRVec3Sites->data->data, 0, params->numBeamPoints * sizeof( REAL4 ));
    }
    else {
      memset( params->cohSNRVec->data->data, 0, numPoints * sizeof( REAL4 ));
    }
  }

  if ( nullStatOut ) {
    if ( params->threeSiteCase ) {
      memset( params->nullStatVec3Sites->data->data, 0, params->numBeamPoints * sizeof( REAL4 ));
    }
    else {
      memset( params->nullStatVec->data->data, 0, numPoints * sizeof( REAL4 ));
    }
  }

  /*CHECK: hardwired to 6 detectors for now */
  for (l=0 ;l< 6 ; l++)
    {
      caseID[l] = params->detIDVec->data[l];
    }

  /* Read in sigmasq and chisq values from filtering */
  sigmasq = params->sigmasqVec->data;
  i = 0;
  for (l=0 ;l<6 ;l++)
    {
      if( caseID[l] )
	{
	  indexarray[i] = l;
          sigmasq4DArray[i] = sigmasq[l];
	  i++;
	}
    }

  /*** get detector-site locations */
  /* FIXME: Replace with
     ~XLALReadIfo functions in the future.
     Also note that the ifo orders in InterferometerNumber and
     lalCachedDetectors are different:
     caseID[0,..,5]=(G1,H1,H2,L1,T1,V1), whereas
     lalCachedDetectors[0,...]=(T1(0),V1(1),G1(2),H2(3),H1(4),L1(5),...
     ...,LAL_NUM_DETECTORS=12)*/
  w=0;
  if( caseID[0] ) detectors[w++] = lalCachedDetectors[2];
  if( caseID[1] ) detectors[w++] = lalCachedDetectors[4];
  if( caseID[2] ) detectors[w++] = lalCachedDetectors[3];
  if( caseID[3] ) detectors[w++] = lalCachedDetectors[5];
  if( caseID[4] ) detectors[w++] = lalCachedDetectors[0];
  if( caseID[5] ) detectors[w++] = lalCachedDetectors[1];

  i = 0;
  j = 0;
  k = 0;
  l = 0;
  switch ( params->numDetectors )
    {
    case 2:
      i = indexarray[0];
      j = indexarray[1];
      snprintf( caseStr, FILENAME_MAX * sizeof(CHAR), "%s%s",idtag[i],idtag[j]);
      break;

    case 3:
      i=indexarray[0];
      j=indexarray[1];
      k=indexarray[2];
      snprintf( caseStr, FILENAME_MAX * sizeof(CHAR), "%s%s%s",idtag[i],idtag[j],idtag[k]);
      break;

    case 4:
      i=indexarray[0];
      j=indexarray[1];
      k=indexarray[2];
      l=indexarray[3];
      snprintf( caseStr, FILENAME_MAX * sizeof(CHAR), "%s%s%s%s",idtag[i],idtag[j],idtag[k],idtag[l]);
      break;
    }

  /* read in CData */
  for ( l=0 ; l < (INT4) params->numDetectors ; l++) {
    cData[l] = input->multiCData->cData[l];
  }

  /* Now construct the appropriate coherent SNR */
  switch (params->numDetectors)  {
  case 2:
    /* Network: H1 and H2*/
    if(caseID[1] && caseID[2]) {
      /* FIXME: This section not readied for using "newsnr" */
      case2a = 1;
      m = 0;
      for (k=0;k<(INT4)numPoints;k++) {
	cohSNR = 0.0;

	for (m=k-buffer; m<k+buffer; m++) {
	  if(m >=0 && m < (INT4) numPoints) {
	    REAL4          snrsq1 = 0.0;
	    REAL4          snrsq2 = 0.0;
	    REAL4          chisqFac1 = 1.0;
	    REAL4          chisqFac2 = 1.0;

	    snrsq1 = pow(crealf(cData[0]->data->data[k]),2) +
	      pow(cimagf(cData[0]->data->data[k]),2);
	    chisqFac1 = pow( (1 + snrsq1/eff_snr_denom_fac)*chisq[0]/
			     (2*chisq_dof[0] -2), 0.25);

	    snrsq2 = pow(crealf(cData[1]->data->data[m]),2) +
	      pow(cimagf(cData[1]->data->data[m]),2);
	    chisqFac2 = pow( (1 + snrsq2/eff_snr_denom_fac)*chisq[1]/
			     (2*chisq_dof[1] -2), 0.25);

	    cohSNRLocalRe = sqrt(sigmasq[1])*crealf(cData[0]->data->data[k])/chisqFac1
	      + sqrt(sigmasq[2])*crealf(cData[1]->data->data[m])/chisqFac2;
	    cohSNRLocalIm = sqrt(sigmasq[1])*cimagf(cData[0]->data->data[k])/chisqFac1
	      + sqrt(sigmasq[2])*cimagf(cData[1]->data->data[m])/chisqFac2;

	    cohSNRLocal = (cohSNRLocalRe*cohSNRLocalRe + cohSNRLocalIm*cohSNRLocalIm) /
	      (sigmasq[1] + sigmasq[2] );

	    if(cohSNRLocal > cohSNR) {
	      cohSNR = cohSNRLocal;
	      quadTemp[0] = cData[0]->data->data[k];
	      quadTemp[1] = cData[1]->data->data[m];
	    }
	  }
	  if( cohSNROut ) params->cohSNRVec->data->data[k]= cohSNR;

	  if ( cohSNR > cohSNRThresh ) {
	    if ( !*eventList ) {
	      /* store the start of the crossing */
	      eventStartIdx = k;

	      /* if this is the first event, start the list */
	      thisEvent = *eventList = (MultiInspiralTable *)
		LALCalloc( 1, sizeof(MultiInspiralTable) );

	      if ( !thisEvent ) {
		ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
	      }

	      /* record the data that we need for the clustering algorithm */
	      tempTime = cData[0]->epoch.gpsSeconds + 1.0e-9 * cData[0]->epoch.gpsNanoSeconds + k * deltaT;
	      fracpart = modf( tempTime, &intpart );
	      thisEvent->end_time.gpsSeconds = (INT4) intpart;
	      thisEvent->end_time.gpsNanoSeconds = (INT4) ( 1.0e9*fracpart );

	      /* Fill thisEvent with this trigger's quadrature values */
	      XLALAssignEventQuads(caseID,quadTemp,thisEvent);
	      thisEvent->snr = cohSNR;
	      strcpy(thisEvent->ifos,caseStr);
	      thisEvent->mass1 = input->tmplt->mass1;
	      thisEvent->mass2 = input->tmplt->mass2;
	      thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
	      thisEvent->eta = input->tmplt->eta;
              /* Compute null-statistic for H1-H2 at just trigger end-time */
              nullNorm = ( 1.0 / sigmasq[1]  + 1.0 /  sigmasq[2] );
              nullStatRe = crealf(thisEvent->h1quad) / sqrt(sigmasq[1])
                - crealf(thisEvent->h2quad) / sqrt(sigmasq[2]);
              nullStatIm = cimagf(thisEvent->h1quad) / sqrt(sigmasq[1])
                - cimagf(thisEvent->h2quad) / sqrt(sigmasq[2]);
              thisEvent->null_stat_h1h2 = ( nullStatRe*nullStatRe + nullStatIm*nullStatIm ) / nullNorm ;
              thisEvent->null_statistic = thisEvent->null_stat_h1h2;
	      /*Calculate distance/effective distance */
              XLALCoherentCBCEstimateDistanceCase2a(
		  (double) crealf(quadTemp[0]),(double) cimagf(quadTemp[0]),
		  (double) crealf(quadTemp[1]),(double) cimagf(quadTemp[1]),
		  sigmasq4DArray, chisq, thisEvent);

	      thisEvent->ligo_angle = -1001;
	      thisEvent->coa_phase = -1001;
              thisEvent->ligo_angle_sig = -1001;

	      tempTime = 0.0;
	      fracpart = 0.0;
	      intpart = 0.0;
	      fflush( stdout );

	    } /* done creating a new event; closes "if ( !*eventList )" */
	    else if (params->maximizeOverChirp && k <= (eventStartIdx + deltaEventIndex) && cohSNR > thisEvent->snr ) {
	      /* if this is the same event, update the maximum */
	      tempTime = cData[0]->epoch.gpsSeconds + 1.0e-9 * cData[0]->epoch.gpsNanoSeconds + k * deltaT;
	      fracpart = modf( tempTime, &intpart );
	      thisEvent->end_time.gpsSeconds = (INT4) intpart;
	      thisEvent->end_time.gpsNanoSeconds = (INT4) ( 1.0e9*fracpart );

	      /* Fill thisEvent with this trigger's quadrature values */
	      XLALAssignEventQuads(caseID,quadTemp,thisEvent);
	      thisEvent->snr = cohSNR;
	      strcpy(thisEvent->ifos, caseStr);
	      thisEvent->mass1 = input->tmplt->mass1;
	      thisEvent->mass2 = input->tmplt->mass2;
	      thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
	      thisEvent->eta = input->tmplt->eta;
              /* Compute null-statistic for H1-H2 at just trigger end-time */
              nullNorm = ( 1.0 / sigmasq[1]  + 1.0 /  sigmasq[2] );
              nullStatRe = crealf(thisEvent->h1quad) / sqrt(sigmasq[1])
                - crealf(thisEvent->h2quad) / sqrt(sigmasq[2]);
              nullStatIm = cimagf(thisEvent->h1quad) / sqrt(sigmasq[1])
                - cimagf(thisEvent->h2quad) / sqrt(sigmasq[2]);
              thisEvent->null_stat_h1h2 = ( nullStatRe*nullStatRe + nullStatIm*nullStatIm ) / nullNorm ;
              thisEvent->null_statistic = thisEvent->null_stat_h1h2;
	      /*Calculate distance/effective distance */
              XLALCoherentCBCEstimateDistanceCase2a(
                  (double) crealf(quadTemp[0]),(double) cimagf(quadTemp[0]),
		  (double) crealf(quadTemp[1]),(double) cimagf(quadTemp[1]),
                  sigmasq4DArray, chisq, thisEvent);

	      thisEvent->ligo_angle = -1001;
	      thisEvent->coa_phase = -1001;
              thisEvent->ligo_angle_sig = -1001;

	      tempTime = 0.0;
	      fracpart = 0.0;
	      intpart = 0.0;
	      fflush( stdout );

	    } /* done updating event; closes "else if ( params->maximizeOverChirp &&...) */
	    else if ( k > (eventStartIdx + deltaEventIndex) || !(params->maximizeOverChirp) ) {
	      /* clean up this event */
	      MultiInspiralTable      *lastEvent = NULL;

	      /* allocate memory for the newEvent */
	      lastEvent = thisEvent;
	      lastEvent->next = thisEvent = (MultiInspiralTable *)
		LALCalloc( 1, sizeof(MultiInspiralTable) );
	      if ( !(lastEvent->next) )
		{
		  ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
		}

	      /* stick minimal data into the event */

	      tempTime = cData[0]->epoch.gpsSeconds + 1.0e-9 * cData[0]->epoch.gpsNanoSeconds + k * deltaT;
	      fracpart = modf( tempTime, &intpart );
	      thisEvent->end_time.gpsSeconds = (INT4) intpart;
	      thisEvent->end_time.gpsNanoSeconds = (INT4) ( 1.0e9*fracpart );

	      /* Fill thisEvent with this trigger's quadrature values */
	      XLALAssignEventQuads(caseID,quadTemp,thisEvent);
	      thisEvent->snr = cohSNR;
	      strcpy(thisEvent->ifos,caseStr);
	      thisEvent->mass1 = input->tmplt->mass1;
	      thisEvent->mass2 = input->tmplt->mass2;
	      thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
	      thisEvent->eta = input->tmplt->eta;
              /* Compute null-statistic for H1-H2 at just trigger end-time */
              nullNorm = ( 1.0 / sigmasq[1]  + 1.0 /  sigmasq[2] );
              nullStatRe = crealf(thisEvent->h1quad) / sqrt(sigmasq[1])
                - crealf(thisEvent->h2quad) / sqrt(sigmasq[2]);
              nullStatIm = cimagf(thisEvent->h1quad) / sqrt(sigmasq[1])
                - cimagf(thisEvent->h2quad) / sqrt(sigmasq[2]);
              thisEvent->null_stat_h1h2 = ( nullStatRe*nullStatRe + nullStatIm*nullStatIm ) / nullNorm ;
              thisEvent->null_statistic = thisEvent->null_stat_h1h2;
	      /*Calculate distance/effective distance */
              XLALCoherentCBCEstimateDistanceCase2a(
                  (double) crealf(quadTemp[0]),(double) cimagf(quadTemp[0]),
                  (double) crealf(quadTemp[1]),(double) cimagf(quadTemp[1]),
                  sigmasq4DArray, chisq, thisEvent);

	      thisEvent->ligo_angle = -1001;
	      thisEvent->coa_phase = -1001;
              thisEvent->ligo_angle_sig = -1001;

	      /* Need to initialize the event start index to new value */
	      if( k > (eventStartIdx + deltaEventIndex) )
		{
		  eventStartIdx = k;
		}
	      tempTime = 0.0;
	      fracpart= 0.0;
	      intpart = 0.0;

	    }
	  } /* matches if (cohSNR > cohSNRThresh) */

	} /* matches for(m=k-buffer....) */
      }/* matches for(k=0;... */
    }
    else
      { /* Network: 2 detectors excluding either H1, H2, or both H1 and H2 */
	/*Here, the time delay looping must start */
	/* Now calculate the distance (in meters) */
	case2b = 1;
	REAL4          chisqFac1Sq = 1.0;
	REAL4          chisqFac2Sq = 1.0;
	REAL4          eff_snr_denom_fac_inv2 = 1.0;

	cohSNR = 0.0;
	eff_snr_denom_fac_inv2 = 2.0/eff_snr_denom_fac;

	if ( (chisq[0]/(2*chisq_dof[0] -2)) > 1.0 ) {
	  /* Compute the denominator in new-snr for this ifo */
	  chisqFac1Sq = pow( 0.5*(1 + pow(chisq[0]/(2*chisq_dof[0] -2), 0.5*eff_snr_denom_fac)), eff_snr_denom_fac_inv2);
	}

	if ( (chisq[1]/(2*chisq_dof[1] -2)) > 1.0 ) {
	  /* Compute the denominator in new-snr for this ifo */
	  chisqFac2Sq = pow( 0.5*(1 + pow(chisq[1]/(2*chisq_dof[1] -2), 0.5*eff_snr_denom_fac)), eff_snr_denom_fac_inv2);
	}

	for (i=0;i<3;i++) {
	  s[i] = (REAL4) ( detectors[1].location[i] - detectors[0].location[i]);
	}
	distance[1] = sqrt( cartesianInnerProduct(s,s) );
	timeDelay[1] = distance[1]/LAL_C_SI;
	slidePoints[1] = rint( (fabs(timeDelay[1])/deltaT) + 1.0 );

	k = 0;
	q = 0;
	w = 0;

	for(k=0 ; k<(INT4)numPoints ; k++)
	  {
            REAL4          snrsq1 = 0.0;
            REAL4          snrsq2 = 0.0;

            snrsq1 = pow(crealf(cData[0]->data->data[k]),2) +
              pow(cimagf(cData[0]->data->data[k]),2);

	    for (q = k-slidePoints[1]-buffer; q < k+slidePoints[1]+buffer; q++)
	      {
		if(q >= 0 && q < (INT4) numPoints)
		  {
		    snrsq2 = pow(crealf(cData[1]->data->data[q]),2) +
		      pow(cimagf(cData[1]->data->data[q]),2);

		    cohSNRLocal = snrsq1 + snrsq2;

		    if(cohSNRLocal > cohSNR)
		      {
			cohSNR = cohSNRLocal;

			quadTemp[0] = cData[0]->data->data[k];
			quadTemp[1] = cData[1]->data->data[q];
			w = q;

		      }
		  }
	      }
	    if( cohSNROut ) params->cohSNRVec->data->data[k] = cohSNR;

	    if ( cohSNR > cohSNRThresh ) {
	      if ( !*eventList ) {
		/* store the start of the crossing */
		eventStartIdx = k;

		/* if this is the first event, start the list */
		thisEvent = *eventList = (MultiInspiralTable *)
		  LALCalloc( 1, sizeof(MultiInspiralTable) );

		if ( !thisEvent )
		  {
		    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
		  }

		/* record the data that we need for the clustering algorithm */
		tempTime = cData[0]->epoch.gpsSeconds + 1.0e-9 * cData[0]->epoch.gpsNanoSeconds + k * deltaT;
		fracpart = modf( tempTime, &intpart );
		thisEvent->end_time.gpsSeconds = (INT4) intpart;
		thisEvent->end_time.gpsNanoSeconds = (INT4) ( 1.0e9*fracpart );

		/* Fill thisEvent with this trigger's quadrature values */
		XLALAssignEventQuads(caseID,quadTemp,thisEvent);
		thisEvent->snr = cohSNR;
		strcpy(thisEvent->ifos,caseStr);
		thisEvent->mass1 = input->tmplt->mass1;
		thisEvent->mass2 = input->tmplt->mass2;
		thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
		thisEvent->eta = input->tmplt->eta;
                /* With two non-coaligned ifo, the null-statistic is not meaningful */
                thisEvent->null_statistic = -1;
		thisEvent->kappa = nullStatRegul;
		thisEvent->tau5 = snrsq1/chisqFac1Sq + snrsq2/chisqFac2Sq;
	        /*Calculate distance/effective distance */
                /* NOTE: Using snrsq in sigmasq4DArray */
                XLALCoherentCBCEstimateDistanceCase2b( caseID,
		    (double) crealf(quadTemp[0]),(double) cimagf(quadTemp[0]),
		    (double) crealf(quadTemp[1]),(double) cimagf(quadTemp[1]),
                    sigmasq4DArray, chisq, thisEvent);

		thisEvent->coa_phase = -1001;

		thisEvent->ligo_angle = acos( LAL_C_SI * deltaT * abs(k-w) / distance[1] );
		if( (k-w) > 0 )
		  {
		    thisEvent->ligo_angle_sig = 1;
		  }
		else
		  {
		    thisEvent->ligo_angle_sig = -1;
		  }

		tempTime = 0.0;
		fracpart = 0.0;
		intpart = 0.0;
		fflush( stdout );

	      } /* done creating a new event */
	      else if (params->maximizeOverChirp && k <= (eventStartIdx + deltaEventIndex) && cohSNR > thisEvent->snr ) {
		/* if this is the same event, update the maximum */
		tempTime = cData[0]->epoch.gpsSeconds + 1.0e-9 * cData[0]->epoch.gpsNanoSeconds + k * deltaT;
		fracpart = modf( tempTime, &intpart );
		thisEvent->end_time.gpsSeconds = (INT4) intpart;
		thisEvent->end_time.gpsNanoSeconds = (INT4) ( 1.0e9*fracpart );
		/*	obtainNetworkPhase(caseID, quadTemp, thisEvent);*/

		/* Fill thisEvent with this trigger's quadrature values */
		XLALAssignEventQuads(caseID,quadTemp,thisEvent);
		thisEvent->snr = cohSNR;
		strcpy(thisEvent->ifos, caseStr);
		thisEvent->mass1 = input->tmplt->mass1;
		thisEvent->mass2 = input->tmplt->mass2;
		thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
		thisEvent->eta = input->tmplt->eta;
                /* With two non-coaligned ifo, the null-statistic is not meaningful */
                thisEvent->null_statistic = -1;
		thisEvent->kappa = nullStatRegul;
		thisEvent->tau5 = snrsq1/chisqFac1Sq + snrsq2/chisqFac2Sq;
	        /*Calculate effective distance */
                XLALCoherentCBCEstimateDistanceCase2b( caseID,
		    (double) crealf(quadTemp[0]),(double) cimagf(quadTemp[0]),
		    (double) crealf(quadTemp[1]),(double) cimagf(quadTemp[1]),
                    sigmasq4DArray, chisq, thisEvent);

		thisEvent->coa_phase = -1001;
		thisEvent->ligo_angle = acos( LAL_C_SI * deltaT * abs(k-w) / distance[1] );
		if( (k-w) > 0 )
		  {
		    thisEvent->ligo_angle_sig = 1;
		  }
		else
		  {
		    thisEvent->ligo_angle_sig = -1;
		  }

		tempTime = 0.0;
		fracpart = 0.0;
		intpart = 0.0;
		fflush( stdout );

	         }
	      else if ( k > (eventStartIdx + deltaEventIndex) || !(params->maximizeOverChirp) ) {
		/* clean up this event */
		MultiInspiralTable      *lastEvent = NULL;

		/* allocate memory for the newEvent */
		lastEvent = thisEvent;
		lastEvent->next = thisEvent = (MultiInspiralTable *)
		  LALCalloc( 1, sizeof(MultiInspiralTable) );
		if ( !(lastEvent->next) )
		  {
		    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
		  }

		/* stick minimal data into the event */

		tempTime = cData[0]->epoch.gpsSeconds + 1.0e-9 * cData[0]->epoch.gpsNanoSeconds + k * deltaT;
		fracpart = modf( tempTime, &intpart );
		thisEvent->end_time.gpsSeconds = (INT4) intpart;
		thisEvent->end_time.gpsNanoSeconds = (INT4) ( 1.0e9*fracpart );
		/*	obtainNetworkPhase(caseID, quadTemp, thisEvent);*/

		/* Fill thisEvent with this trigger's quadrature values */
		XLALAssignEventQuads(caseID,quadTemp,thisEvent);
		thisEvent->snr = cohSNR;
		strcpy(thisEvent->ifos,caseStr);
		thisEvent->mass1 = input->tmplt->mass1;
		thisEvent->mass2 = input->tmplt->mass2;
		thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
		thisEvent->eta = input->tmplt->eta;
                /* With two non-coaligned ifo, the null-statistic is not meaningful */
                thisEvent->null_statistic = -1;
		thisEvent->kappa = nullStatRegul;
		thisEvent->tau5 = snrsq1/chisqFac1Sq + snrsq2/chisqFac2Sq;
	        /*Calculate effective distance */
                XLALCoherentCBCEstimateDistanceCase2b( caseID,
		    (double) crealf(quadTemp[0]),(double) cimagf(quadTemp[0]),
		    (double) crealf(quadTemp[1]),(double) cimagf(quadTemp[1]),
                    sigmasq4DArray, chisq, thisEvent);

		thisEvent->coa_phase = -1001;
                thisEvent->ligo_angle = acos( LAL_C_SI * deltaT * abs(k-w) / distance[1] );
		if( (k-w) > 0 )
		  {
		    thisEvent->ligo_angle_sig = 1;
		  }
		else
		  {
		    thisEvent->ligo_angle_sig = -1;
		  }

		/* Need to initialize the event start index to new value */
		if( k > (eventStartIdx + deltaEventIndex) )
		  {
		    eventStartIdx = k;
		  }

		tempTime = 0.0;
		fracpart= 0.0;
		intpart = 0.0;

	      }
	    } /* matches if (cohSNR > cohSNRThresh) */
	  }
      }
    break;
  case 3: /* Network: 3 detectors including both H1 and H2 */
    if(caseID[1] && caseID[2])
      {
	case3a = 1;
	/* Now calculate the distance (in meters) between LHO and 2nd site*/
	for (i=0;i<3;i++) {
	  s[i] = (REAL4) ( detectors[2].location[i] - detectors[0].location[i]);
	}
	distance[1] = sqrt( cartesianInnerProduct( s,s) );
	timeDelay[1] = distance[1]/LAL_C_SI;
	slidePoints[1] = rint( (fabs(timeDelay[1])/deltaT) + 1.0 );

	k = 0;
	q = 0;
	m = 0;
	w = 0;

	for(k=0;k<(INT4)numPoints;k++)
	  {
	    cohSNR = 0.0;
	    for(m=k-buffer;m<k+buffer;m++)
	      {
		if(m >=0 && m < (INT4) numPoints)
		  {
		    for (q = m-slidePoints[1]-buffer;q < m+slidePoints[1]+buffer;q++)
		      {
			if(q >= 0 && q < (INT4) numPoints)
			  {
			    /*FIXME: This will NOT work if G1 is present!
                              because it assumes that the "0" det is H1 and
                              "1" det is H2! Rectify in next rev. */
			    REAL4          snrsq1 = 0.0;
			    REAL4          snrsq2 = 0.0;
			    REAL4          snrsq3 = 0.0;
			    REAL4          chisqFac1 = 1.0;
			    REAL4          chisqFac2 = 1.0;
			    REAL4          chisqFac3 = 1.0;

			    snrsq1 = pow(crealf(cData[0]->data->data[k]),2) +
			      pow(cimagf(cData[0]->data->data[k]),2);
			    chisqFac1 = pow( (1 + snrsq1/eff_snr_denom_fac)*chisq[0]/
					     (2*chisq_dof[0] -2), 0.25);

			    snrsq2 = pow(crealf(cData[1]->data->data[m]),2) +
			      pow(cimagf(cData[1]->data->data[m]),2);
			    chisqFac2 = pow( (1 + snrsq2/eff_snr_denom_fac)*chisq[1]/
					     (2*chisq_dof[1] -2), 0.25);

			    snrsq3 = pow(crealf(cData[2]->data->data[q]),2) +
			      pow(cimagf(cData[2]->data->data[q]),2);
			    chisqFac3 = pow( (1 + snrsq3/eff_snr_denom_fac)*chisq[2]/
					     (2*chisq_dof[2] -2), 0.5);

			    cohSNRLocalRe = sqrt(sigmasq[1])*crealf(cData[0]->data->data[k])/
			      chisqFac1
			      + sqrt(sigmasq[2])*crealf(cData[1]->data->data[m])/chisqFac2;
			    cohSNRLocalIm = sqrt(sigmasq[1])*cimagf(cData[0]->data->data[k])/
			      chisqFac1
			      + sqrt(sigmasq[2])*cimagf(cData[1]->data->data[m])/chisqFac2;

			    cohSNRLocal = (cohSNRLocalRe*cohSNRLocalRe + cohSNRLocalIm*cohSNRLocalIm) /
			      (sigmasq[1] + sigmasq[2]) ;

			    cohSNRLocal += (crealf(cData[2]->data->data[q])
					    *crealf(cData[2]->data->data[q]) +
					    cimagf(cData[2]->data->data[q])
					    *cimagf(cData[2]->data->data[q]))/chisqFac3;

			  if(cohSNRLocal > cohSNR)
			      {
				cohSNR = cohSNRLocal;

				quadTemp[0] = cData[0]->data->data[k];
				quadTemp[1] = cData[1]->data->data[m];
				quadTemp[2] = cData[2]->data->data[q];
				w = q;

			      }
			  }
		      }
		  }
		if( cohSNROut ) params->cohSNRVec->data->data[k] = cohSNR;

		if ( cohSNR > cohSNRThresh ) {
		  if ( !*eventList ) {
		    /* store the start of the crossing */
		    eventStartIdx = k;

		    /* if this is the first event, start the list */
		    thisEvent = *eventList = (MultiInspiralTable *)
		      LALCalloc( 1, sizeof(MultiInspiralTable) );

		    if ( !thisEvent )
		      {
			ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
		      }

		    /* record the data that we need for the clustering algorithm */
		    tempTime = cData[0]->epoch.gpsSeconds + 1.0e-9 * cData[0]->epoch.gpsNanoSeconds + k * deltaT;
		    fracpart = modf( tempTime, &intpart );
		    thisEvent->end_time.gpsSeconds = (INT4) intpart;
		    thisEvent->end_time.gpsNanoSeconds = (INT4) ( 1.0e9*fracpart );

		    /* Fill thisEvent with this trigger's quadrature values */
		    XLALAssignEventQuads(caseID,quadTemp,thisEvent);
		    thisEvent->snr = cohSNR;
                    cohSNRLocalRe = sqrt(sigmasq[1])*crealf(thisEvent->h1quad)
                                    + sqrt(sigmasq[2])*crealf(thisEvent->h2quad);
                    cohSNRLocalIm = sqrt(sigmasq[1])*cimagf(thisEvent->h1quad)
                                    + sqrt(sigmasq[2])*cimagf(thisEvent->h2quad);
                    thisEvent->coh_snr_h1h2 = (cohSNRLocalRe*cohSNRLocalRe +
		       cohSNRLocalIm*cohSNRLocalIm) /(sigmasq[1] + sigmasq[2] );
		    strcpy(thisEvent->ifos,caseStr);
		    thisEvent->mass1 = input->tmplt->mass1;
		    thisEvent->mass2 = input->tmplt->mass2;
		    thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
		    thisEvent->eta = input->tmplt->eta;
                    /* Compute null-statistic for H1-H2 at just trigger end-time */
                    nullNorm = ( 1.0 / sigmasq[1]  + 1.0 /  sigmasq[2] );
                    nullStatRe = crealf(thisEvent->h1quad) / sqrt(sigmasq[1])
                      - crealf(thisEvent->h2quad) / sqrt(sigmasq[2]);
                    nullStatIm = cimagf(thisEvent->h1quad) / sqrt(sigmasq[1])
                      - cimagf(thisEvent->h2quad) / sqrt(sigmasq[2]);
                    thisEvent->null_stat_h1h2 = ( nullStatRe*nullStatRe + nullStatIm*nullStatIm ) / nullNorm ;
		    /*Calculate distance/effective distance */
		    XLALCoherentCBCEstimateDistanceCase3a( caseID,
                       (double) crealf(quadTemp[0]),(double) cimagf(quadTemp[0]),
                       (double) crealf(quadTemp[1]),(double) cimagf(quadTemp[1]),
		       (double) crealf(quadTemp[2]),(double) cimagf(quadTemp[2]),
		       sigmasq4DArray, chisq, thisEvent);

		    thisEvent->null_statistic = -1; /* store network null-statistic for numDetectors >2*/
		    thisEvent->kappa = nullStatRegul;

		    thisEvent->ligo_angle = acos( LAL_C_SI * deltaT * abs(k-w) / distance[1] );
                    thisEvent->coa_phase = -1001;

		    if( (k-w) > 0 )
		      {
			thisEvent->ligo_angle_sig = 1;
		      }
		    else
		      {
			thisEvent->ligo_angle_sig = -1;
		      }

		    tempTime = 0.0;
		    fracpart = 0.0;
		    intpart = 0.0;
		    fflush( stdout );

		  } /* done creating a new event */
		  else if (params->maximizeOverChirp && k <= (eventStartIdx + deltaEventIndex) && cohSNR > thisEvent->snr ) {
		    /* if this is the same event, update the maximum */
		    tempTime = cData[0]->epoch.gpsSeconds + 1.0e-9 * cData[0]->epoch.gpsNanoSeconds + k * deltaT;
		    fracpart = modf( tempTime, &intpart );
		    thisEvent->end_time.gpsSeconds = (INT4) intpart;
		    thisEvent->end_time.gpsNanoSeconds = (INT4) ( 1.0e9*fracpart );

		    /* Fill thisEvent with this trigger's quadrature values */
		    XLALAssignEventQuads(caseID,quadTemp,thisEvent);
		    thisEvent->snr = cohSNR;
                    cohSNRLocalRe = sqrt(sigmasq[1])*crealf(thisEvent->h1quad)
                                    + sqrt(sigmasq[2])*crealf(thisEvent->h2quad);
                    cohSNRLocalIm = sqrt(sigmasq[1])*cimagf(thisEvent->h1quad)
                                    + sqrt(sigmasq[2])*cimagf(thisEvent->h2quad);
                    thisEvent->coh_snr_h1h2 = (cohSNRLocalRe*cohSNRLocalRe +
                      cohSNRLocalIm*cohSNRLocalIm) /(sigmasq[1] + sigmasq[2] );
		    strcpy(thisEvent->ifos, caseStr);
		    thisEvent->mass1 = input->tmplt->mass1;
		    thisEvent->mass2 = input->tmplt->mass2;
		    thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
		    thisEvent->eta = input->tmplt->eta;
                    /* Compute null-statistic for H1-H2 at just trigger end-time */
                    nullNorm = ( 1.0 / sigmasq[1]  + 1.0 /  sigmasq[2] );
                    nullStatRe = crealf(thisEvent->h1quad) / sqrt(sigmasq[1])
                      - crealf(thisEvent->h2quad) / sqrt(sigmasq[2]);
                    nullStatIm = cimagf(thisEvent->h1quad) / sqrt(sigmasq[1])
                      - cimagf(thisEvent->h2quad) / sqrt(sigmasq[2]);
                    thisEvent->null_stat_h1h2 = ( nullStatRe*nullStatRe + nullStatIm*nullStatIm ) / nullNorm ;
		    /*Calculate distance/effective distance */
		    XLALCoherentCBCEstimateDistanceCase3a( caseID,
                       (double) crealf(quadTemp[0]),(double) cimagf(quadTemp[0]),
                       (double) crealf(quadTemp[1]),(double) cimagf(quadTemp[1]),
		       (double) crealf(quadTemp[2]),(double) cimagf(quadTemp[2]),
                       sigmasq4DArray, chisq, thisEvent);

		    thisEvent->null_statistic = -1; /* store network null-statistic for numDetectors >2*/
		    thisEvent->kappa = nullStatRegul;
		    thisEvent->ligo_angle = acos( LAL_C_SI * deltaT * abs(k-w) / distance[1] );
                    thisEvent->coa_phase = -1001;

		    if( (k-w) > 0 )
		      {
			thisEvent->ligo_angle_sig = 1;
		      }
		    else
		      {
			thisEvent->ligo_angle_sig = -1;
		      }

		    tempTime = 0.0;
		    fracpart = 0.0;
		    intpart = 0.0;
		    fflush( stdout );

		  }
		  else if ( k > (eventStartIdx + deltaEventIndex) || !(params->maximizeOverChirp) ) {
		    /* clean up this event */
		    MultiInspiralTable      *lastEvent = NULL;

		    /* allocate memory for the newEvent */
		    lastEvent = thisEvent;
		    lastEvent->next = thisEvent = (MultiInspiralTable *)
		      LALCalloc( 1, sizeof(MultiInspiralTable) );
		    if ( !(lastEvent->next) )
		      {
			ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
		      }

		    /* stick minimal data into the event */

		    tempTime = cData[0]->epoch.gpsSeconds + 1.0e-9 * cData[0]->epoch.gpsNanoSeconds + k * deltaT;
		    fracpart = modf( tempTime, &intpart );
		    thisEvent->end_time.gpsSeconds = (INT4) intpart;
		    thisEvent->end_time.gpsNanoSeconds = (INT4) ( 1.0e9*fracpart );

		    /* Fill thisEvent with this trigger's quadrature values */
		    XLALAssignEventQuads(caseID,quadTemp,thisEvent);
		    thisEvent->snr = cohSNR;
                    cohSNRLocalRe = sqrt(sigmasq[1])*crealf(thisEvent->h1quad)
                                    + sqrt(sigmasq[2])*crealf(thisEvent->h2quad);
                    cohSNRLocalIm = sqrt(sigmasq[1])*cimagf(thisEvent->h1quad)
                                    + sqrt(sigmasq[2])*cimagf(thisEvent->h2quad);
                    thisEvent->coh_snr_h1h2 = (cohSNRLocalRe*cohSNRLocalRe +
                      cohSNRLocalIm*cohSNRLocalIm) /(sigmasq[1] + sigmasq[2] );
		    strcpy(thisEvent->ifos,caseStr);
		    thisEvent->mass1 = input->tmplt->mass1;
		    thisEvent->mass2 = input->tmplt->mass2;
		    thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
		    thisEvent->eta = input->tmplt->eta;
                    /* Compute null-statistic for H1-H2 at just trigger end-time */
                    nullNorm = ( 1.0 / sigmasq[1]  + 1.0 /  sigmasq[2] );
                    nullStatRe = crealf(thisEvent->h1quad) / sqrt(sigmasq[1])
                      - crealf(thisEvent->h2quad) / sqrt(sigmasq[2]);
                    nullStatIm = cimagf(thisEvent->h1quad) / sqrt(sigmasq[1])
                      - cimagf(thisEvent->h2quad) / sqrt(sigmasq[2]);
                    thisEvent->null_stat_h1h2 = ( nullStatRe*nullStatRe + nullStatIm*nullStatIm ) / nullNorm ;
		    /*Calculate distance/effective distance */
		    XLALCoherentCBCEstimateDistanceCase3a( caseID,
                       (double) crealf(quadTemp[0]),(double) cimagf(quadTemp[0]),
                       (double) crealf(quadTemp[1]),(double) cimagf(quadTemp[1]),
		       (double) crealf(quadTemp[2]),(double) cimagf(quadTemp[2]),
                       sigmasq4DArray, chisq, thisEvent);

		    thisEvent->null_statistic = -1; /* store network null-statistic for numDetectors >2*/
		    thisEvent->kappa = nullStatRegul;
		    thisEvent->ligo_angle = acos( LAL_C_SI * deltaT * abs(k-w) / distance[1] );
                    thisEvent->coa_phase = -1001;

		    if( (k-w) > 0 )
		      {
			thisEvent->ligo_angle_sig = 1;
		      }
		    else
		      {
			thisEvent->ligo_angle_sig = -1;
		      }

		    /* Need to initialize the event start index to new value */
		    if( k > (eventStartIdx + deltaEventIndex) )
		      {
			eventStartIdx = k;
		      }
		    tempTime = 0.0;
		    fracpart= 0.0;
		    intpart = 0.0;

		  }
		} /* matches if (cohSNR > cohSNRThresh) */
	      }
	  }
      }
    else
      { /* Network: 3 detectors excluding either H1, H2, or both (H1 && H2)*/
	/*Now the last 3 cases will involve the looping over the coefficients*/
	LIGOTimeGPS 	triggerGPSEndTime;/* Needed to calculate time-delays */
	double          psiInRadians = 0.0;
	double          detRefLocation[3];
	double          detNextLocation[3];
        const SkyGrid  *skyGrid;
	UINT4 skyGridIdx=0;

        /* This is case 3b, which pertains to a 3D network
           with 3 ifos at 3 different sites */
        case3b = 1;

	triggerGPSEndTime.gpsSeconds = cData[0]->epoch.gpsSeconds;
	triggerGPSEndTime.gpsNanoSeconds = cData[0]->epoch.gpsNanoSeconds;
	/* Convert GPS time of trigger to GMST time in radians
	   for computing F+, Fx */
	gmstInRadians = XLALGreenwichMeanSiderealTime(&triggerGPSEndTime);

	/* Following needed because XLALArrivalTimeDiff() uses doubles */
	for ( locIdx=0 ; locIdx<3 ; locIdx++ ) {
	  detRefLocation[locIdx] = (double) detectors[0].location[locIdx];
	}

	/* Loop over points in the sky-position grid */
	for ( skyGrid = skyGridPtr ; skyGrid ; skyGrid = skyGrid->next ) {
	  phi = skyGrid->Alpha;
	  theta = skyGrid->Delta;

	  /* Loop over detectors computing their F+, Fx and t_c's */
	  detId = 0;
	  for( j=0; j<LAL_NUM_IFO; j++ ) {
	    /* Compute time-delays if caseID[j] != 0 */
	    if ( !(params->detIDVec->data[j] == 0 )) {
	      for ( locIdx=0 ; locIdx<3 ; locIdx++ ) {
		detNextLocation[locIdx]
		  = (double) detectors[detId].location[locIdx];
	      }

	      /*Arrival time delays in detectors relative to the 1st;
		esp., dt21 = t2 - t1, dt31 = t3 - t1, etc. are computed */
	      /*	timeDelay[detId] = XLALArrivalTimeDiff(
			detectors[detId].location,
			detectors[0].location,
			raDec.longitude,raDec.latitude,&triggerGPSEndTime);
	      */
	      timeDelay[detId] = XLALArrivalTimeDiff(detNextLocation,
				   detRefLocation,phi,theta,&triggerGPSEndTime);

	      /* save unsorted time-delays for sorting after
		 this loop over detectors is exited */
	      sortedDelays3D[detId] = timeDelay[detId];

	      /* round off to nearest integer */
	      slidePoints[detId] = rint( timeDelay[detId]/deltaT );

	      detId++;
	    }
	  }
	  /* Sort the time-delays now */
	  qsort( sortedDelays3D, (int) params->numDetectors, sizeof(double), compare );

	  for ( detId=0 ; detId<params->numDetectors ; detId++ ) {
	    /* Compute the slide boundaries;
	       round off sorted time delays to nearest integer */
	    sortedSlidePoints3D[detId]
	      = rint( sortedDelays3D[detId]/deltaT );
	  }

	  /* Compute antenna-patterns and coherent SNR */
	  /* Loop over detectors computing the theta-, phi-dependent
	     pieces of F+, Fx */
	  detId = 0;
	  for( j=0; j<LAL_NUM_IFO; j++ ) {
	    /* Compute antenna-patterns if caseID[j] != 0 */
	    if ( !(params->detIDVec->data[j] == 0 )) {
	      XLALComputeDetAMResponse(&fplus[detId], &fcross[detId],
				       (const REAL4(*)[3])detectors[detId].response, phi, theta,
				       psiInRadians, (double) gmstInRadians);

	      /* Compute antenna-pattern factors */
	      AAn[detId] = ( fplus[detId]) * ( fplus[detId]);
	      BBn[detId] = ( fplus[detId]) * ( fcross[detId]);
	      CCn[detId] = ( fcross[detId]) * ( fcross[detId]);


	      /* Compute the elements of the helicity-plane projection matrix */
	      AAn[detId] *= sigmasq[j];
	      BBn[detId] *= sigmasq[j];
	      CCn[detId] *= sigmasq[j];

	      /* Calculate factors necessary for parameter estimation*/
	      uSigma[detId] = (double)fplus[detId] * sqrt((double)sigmasq[j]);
	      vSigma[detId] = (double)fcross[detId] * sqrt((double)sigmasq[j]);

	      /* Sort sigmasq for degeneracy resoln. */
	      sortedSigmasq3D[detId] = sigmasq[j]/3;

	      detId++;
	    }
	  }
	  qsort( sortedSigmasq3D, (int) params->numDetectors, sizeof(double), compare );

	  /* Construct network terms and factors required for
	     computing the coherent statistics */
	  AA = AAn[0] + AAn[1] +AAn[2];
	  BB = BBn[0] + BBn[1] +BBn[2];
	  CC = CCn[0] + CCn[1] +CCn[2];

	  discrimSqrt = (AA*AA + 4*BB*BB
			     - 2*AA*CC + CC*CC);
          if ( (discrimSqrt>0.0) ) {
            discrimSqrt = sqrt(discrimSqrt);
          }
          else {
            discrimSqrt = 0.0;
          }

	  O22 = ( -AA + CC + discrimSqrt);
	  O22 /= (BB * sqrt( 4 + pow( AA - CC - discrimSqrt, 2)/( BB*BB ) ) );

	  O21 = 1 / sqrt( 1 + pow(-AA+CC+discrimSqrt, 2)/ ( 4*BB*BB ) );

	  O12 = ( -AA + CC - discrimSqrt);
	  O12 /= ( BB * sqrt( 4 + pow(-AA+CC-discrimSqrt, 2) / ( BB*BB ) ) );

	  O11 = 1 / sqrt( 1 + pow(-AA+CC-discrimSqrt, 2)/ ( 4*BB*BB ) );

	  detId = 0;
	  for( j=0; j<LAL_NUM_IFO; j++ ) {
	    if ( !(params->detIDVec->data[j] == 0 )) {
	      VVPlus[detId] = O11 * ( fplus[detId])
		+ O12 * ( fcross[detId]);
	      VVMinus[detId] = O21 * ( fplus[detId])
		+ O22 * ( fcross[detId]);

	      VVPlus[detId] *= sqrt((REAL4) sigmasq[j]);
	      VVMinus[detId] *= sqrt((REAL4) sigmasq[j]);

	      detId++;
	    }
	  }

	  MM1 = 2*(AA*CC - BB*BB)/(AA+CC+discrimSqrt);
	  MM2 = 2*(AA*CC - BB*BB)/(AA+CC-discrimSqrt);

	  /* Factors for determining & computing degenerate case/statistics*/
	  zetaRe = (AA-CC)/(AA+CC);
	  zetaIm = 2*BB/(AA+CC);
	  zetaAbsSq = zetaRe*zetaRe + zetaIm*zetaIm;
	  degenerateStat = 0;
          if ( (pow(zetaAbsSq - zetaAbsSqMin,2) < 0.1) ) {
	    degenerateStat = 1;
	  }
	  for ( detId=0 ; detId<params->numDetectors ; detId++ ) {
	    DRe[detId] = (uSigma[detId] + zetaRe*uSigma[detId] +
			  zetaIm*vSigma[detId]) / sqrt(AA+CC) ;
	    DIm[detId] = (vSigma[detId] - zetaRe*vSigma[detId] +
			  zetaIm*uSigma[detId]) / sqrt(AA+CC) ;
	  }

	  /* Loop over time-points in the reference detector */
	  for( timePt[0]=0 ; timePt[0]<(INT4)numPoints ; timePt[0]++) {
	    /* Reset cohSNR to zero so that it can only be ratcheted upward by
	       cohSNRLocal computed below for every point in sky-position grid*/
	    cohSNR = 0.0;

	    /* Loop over time-points in reference detector, after
	       accounting for the rounded-off sortedSlidePoints */
	    if( ( timePt[0] < (0-sortedSlidePoints3D[0]) )
		||  ( timePt[0]> (numPoints - 1 - sortedSlidePoints3D[2]) ) ) {
	      cohSNR = 0.0;
	      nullStatistic = 0.0;
	      autoCorrNullSq = 0.0;
	      ratioStatLocal = 0.0;
	      if( cohSNROut ) {
                if ( !params->threeSiteCase ) {
                  params->cohSNRVec->data->data[timePt[0]] = cohSNR;
                }
                else {
                  params->cohSNRVec3Sites->data->data[skyGridIdx] = cohSNR;
                }
              }
	      if( nullStatOut ) {
                if ( !params->threeSiteCase ) {
                  params->nullStatVec->data->data[timePt[0]] = nullStatistic;
                }
                else {
                  params->nullStatVec3Sites->data->data[skyGridIdx] = nullStatistic;
                }
	      }
            }
	    else {
	      if ( (MM1 == 0.0) || (MM2 == 0.0 ) ) {
		cohSnrSqLocal = 0.0;
		ratioStatLocal = 0.0;
	      }
	      else {
		cohSnrSqLocal = XLALComputeCohStatCase3b(cData,VVPlus,VVMinus,timePt[0],slidePoints,params->numDetectors,MM1,MM2,chisq,chisq_dof,eff_snr_denom_fac,&eff_coh_snr);
		/* Note that unlike above the next 3 are *squared* statistics */
		nullStatistic = (REAL4) XLALComputeNonQuadNullStatCase3b(caseID,fplus,fcross,sigmasq,timePt[0],slidePoints,cData);

		autoCorrNullSq = (REAL4) XLALComputeIncohStatCase3b(caseID,fplus,fcross,sigmasq,timePt[0],slidePoints,cData,nullStatistic,&crossCorrNullSq,chisq,chisq_dof,eff_snr_denom_fac);

		ratioStatLocal = cohSnrSqLocal;
	      }
	      /* Updating RatioStat here */
	      if(ratioStatLocal > cohSNR) {
		cohSNR = ratioStatLocal;
		for ( detId=0 ; detId<params->numDetectors ; detId++ ) {
		  detIdSlidTimePt = timePt[0]+slidePoints[detId];
		  quadTemp[detId] = cData[detId]->data->data[detIdSlidTimePt];
		}
		/* Coh-stat and null-stream sky-maps */
                if( cohSNROut ) {
                  if ( !params->threeSiteCase ) {
                    params->cohSNRVec->data->data[timePt[0]] = cohSNR;
                  }
                  else {
                    params->cohSNRVec3Sites->data->data[skyGridIdx] = cohSNR;
                  }
                }
                if( nullStatOut ) {
                  if ( !params->threeSiteCase ) {
                    params->nullStatVec->data->data[timePt[0]] = nullStatistic;
                  }
                  else {
                    params->nullStatVec3Sites->data->data[skyGridIdx] = nullStatistic;
                  }
                }
	      }
	      /* Threshold on RatioStat here; cohSNR was set to RatioStat above */
	      if ( cohSNR > cohSNRThresh ) {
		/* Initialize CData factor for parameter-estimation */
		if ( !*eventList ) {
		  /* store the start of the crossing */
		  eventStartIdx = timePt[0];

		  /* if this is the first event, start the list */
		  thisEvent = *eventList = (MultiInspiralTable *)
		    LALCalloc( 1, sizeof(MultiInspiralTable) );

		  if ( !thisEvent )
		    {
		      ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
		    }
		  /* record the data that we need for the clustering algorithm */
		  tempTime = cData[0]->epoch.gpsSeconds + 1.0e-9 * cData[0]->epoch.gpsNanoSeconds + timePt[0] * deltaT;
		  fracpart = modf( tempTime, &intpart );
		  thisEvent->end_time.gpsSeconds = (INT4) intpart;
		  thisEvent->end_time.gpsNanoSeconds = (INT4) ( 1.0e9*fracpart );

		  /* Fill thisEvent with this trigger's quadrature values */
		  for ( detId=0 ; detId<params->numDetectors ; detId++ ) {
		    detIdSlidTimePt = timePt[0]+slidePoints[detId];
		    quadTemp[detId] = cData[detId]->data->data[detIdSlidTimePt];
		  }
		  XLALAssignEventQuads(caseID,quadTemp,thisEvent);
		  thisEvent->snr = cohSNR;
		  strcpy(thisEvent->ifos,caseStr);
		  thisEvent->mass1 = input->tmplt->mass1;
		  thisEvent->mass2 = input->tmplt->mass2;
		  thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
		  thisEvent->eta = input->tmplt->eta;
		  /* Save a set of alternative statistics */
		  crossCorrCohSq = XLALComputeCrossCohStatCase3b(cData,VVPlus,VVMinus,timePt[0],slidePoints,params->numDetectors,MM1,MM2,&autoCorrCohSq,chisq,chisq_dof,eff_snr_denom_fac,&crossCorrCohSqEff);
		  thisEvent->tau5 = eff_coh_snr;
                  thisEvent->tau0 = crossCorrCohSqEff;
                  thisEvent->kappa = nullStatRegul;
		  XLALSaveStatsInEventTable(cohSnrSqLocal,autoCorrCohSq,crossCorrCohSq,autoCorrNullSq,crossCorrNullSq, MM1, MM2, VVPlus, VVMinus, thisEvent);

		  /* Since not both H1 and H2 are present, the H1H2 null-stat is not meaningful*/
                  if ( (MM1 == 0.0) || (MM2 == 0.0 ) ) {
                    thisEvent->null_stat_h1h2 = -2;
                  }
                  else {
                    thisEvent->null_stat_h1h2 = -1;
                  }
		  /* Compute network null-statistic at just trigger end-time */
                  if ( (MM1 == 0.0) || (MM2 == 0.0 ) ) {
		    thisEvent->null_stat_degen = -100.0;
		    thisEvent->null_statistic = thisEvent->null_stat_degen;
		  }
		  else {
		    thisEvent->null_statistic = (REAL4) XLALComputeNullStatCase3b(caseID,fplus,fcross,sigmasq,thisEvent);
		  }
		  /* Store the sigmasq's into event tables */
		  XLALCoherentCBCSigmasqFor3Sites(caseID,sigmasq4DArray,
				    chisq,thisEvent,case3b);

		  if ( params->estimParams || params->followup ) {
		    /* Parameter estimation: Distance
		       for ( detId=0 ; detId<params->numDetectors ; detId++ ) {*/
		    for ( i=0 ; i < 4 ; i++ ) {
		      NN[i] = 0.0;
		    }
		    for ( detId=0 ; detId < params->numDetectors ; detId++ ) {
		      NN[0] += uSigma[detId] * (double)crealf(quadTemp[detId]);
		      NN[1] += vSigma[detId] * (double)crealf(quadTemp[detId]);
		      NN[2] += uSigma[detId] * (double)cimagf(quadTemp[detId]);
		      NN[3] += vSigma[detId] * (double)cimagf(quadTemp[detId]);
		    }
		    determinantMM = (double)AA*(double)CC-(double)BB*(double)BB;
		    if ( ( determinantMM*determinantMM < 1e-40 ) )
		      determinantMM = 1e-20; /* Note: Saving with positive sign */
		    InvMMAA = (double)CC / determinantMM;
		    InvMMBB = -(double)BB / determinantMM;
		    InvMMCC = (double)AA / determinantMM;
		    aa[0] = InvMMAA*NN[0] + InvMMBB*NN[1];
		    aa[1] = InvMMBB*NN[0] + InvMMCC*NN[1];
		    aa[2] = InvMMAA*NN[2] + InvMMBB*NN[3];
		    aa[3] = InvMMBB*NN[2] + InvMMCC*NN[3];

		    thisEvent->distance = XLALCoherentCBCParamEstim( &polarization,
		      &inclination, &coaPhase, aa[0], aa[1], aa[2], aa[3],
		      amplitudeConst, thisEvent,
		      (double) chirpTime,(double) crealf(quadTemp[0]),(double) cimagf(quadTemp[0]),
		      (double) crealf(quadTemp[1]),(double) cimagf(quadTemp[1]),(double) crealf(quadTemp[2]),
                      (double) cimagf(quadTemp[2]), (double) crealf(quadTemp[3]),
                      (double) cimagf(quadTemp[3]),sigmasq4DArray, case3b,caseID);

		    thisEvent->inclination = (REAL4) inclination;
		    thisEvent->polarization = (REAL4) polarization;
		    thisEvent->coa_phase = (REAL4) coaPhase;
		    thisEvent->ra = (REAL4) phi;
		    thisEvent->dec = (REAL4) theta;
		  }/*ends if condition on parameter-estimation or followup */

		  tempTime = 0.0;
		  fracpart = 0.0;
		  intpart = 0.0;
		  fflush( stdout );

		} /* done creating a new event */
		else if (params->maximizeOverChirp && timePt[0] <= (eventStartIdx + deltaEventIndex) && cohSNR > thisEvent->snr ) {
		  /* if this is the same event, update the maximum */
		  tempTime = cData[0]->epoch.gpsSeconds + 1.0e-9 * cData[0]->epoch.gpsNanoSeconds + timePt[0] * deltaT;
		  fracpart = modf( tempTime, &intpart );
		  thisEvent->end_time.gpsSeconds = (INT4) intpart;
		  thisEvent->end_time.gpsNanoSeconds = (INT4) ( 1.0e9*fracpart );

		  /* Fill thisEvent with this trigger's quadrature values */
		  for ( detId=0 ; detId<params->numDetectors ; detId++ ) {
		    detIdSlidTimePt = timePt[0]+slidePoints[detId];
		    quadTemp[detId] = cData[detId]->data->data[detIdSlidTimePt];
		  }
		  XLALAssignEventQuads(caseID,quadTemp,thisEvent);
		  thisEvent->snr = cohSNR;
		  strcpy(thisEvent->ifos, caseStr);
		  thisEvent->mass1 = input->tmplt->mass1;
		  thisEvent->mass2 = input->tmplt->mass2;
		  thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
		  thisEvent->eta = input->tmplt->eta;
		  /* Save a set of alternative statistics */
		  crossCorrCohSq = XLALComputeCrossCohStatCase3b(cData,VVPlus,VVMinus,timePt[0],slidePoints,params->numDetectors,MM1,MM2,&autoCorrCohSq,chisq,chisq_dof,eff_snr_denom_fac,&crossCorrCohSqEff);
		  XLALSaveStatsInEventTable(cohSnrSqLocal,autoCorrCohSq,crossCorrCohSq,autoCorrNullSq,crossCorrNullSq, MM1, MM2, VVPlus, VVMinus, thisEvent);

		  /* Since not both H1 and H2 are present, the H1H2 null-stat is not meaningful*/
                  if ( (MM1 == 0.0) || (MM2 == 0.0 ) ) {
		    thisEvent->tau5 = eff_coh_snr;
                    thisEvent->tau0 = crossCorrCohSqEff;
                    thisEvent->kappa = nullStatRegul;
                    thisEvent->null_stat_h1h2 = -2;
                  }
                  else {
                    thisEvent->null_stat_h1h2 = -1;
                  }

		  /* Compute network null-statistic at just trigger end-time */
                  if ( (MM1 == 0.0) || (MM2 == 0.0 ) ) {
		    thisEvent->null_stat_degen = -100.0;
		    thisEvent->null_statistic = thisEvent->null_stat_degen;
		  }
		  else {
		    thisEvent->null_statistic = (REAL4) XLALComputeNullStatCase3b(caseID,fplus,fcross,sigmasq,thisEvent);
		  }
		  /* Store the sigmasq's into event tables */
		  XLALCoherentCBCSigmasqFor3Sites(caseID,sigmasq4DArray,
				    chisq,thisEvent,case3b);
		  if ( params->estimParams || params->followup ) {
		    /* Parameter estimation: Distance
		       for ( detId=0 ; detId<params->numDetectors ; detId++ ) {*/
		    for ( i=0 ; i < 4 ; i++ ) {
		      NN[i] = 0.0;
		    }
		    for ( detId=0 ; detId < params->numDetectors ; detId++ ) {
		      NN[0] += uSigma[detId] * (double)crealf(quadTemp[detId]);
		      NN[1] += vSigma[detId] * (double)crealf(quadTemp[detId]);
		      NN[2] += uSigma[detId] * (double)cimagf(quadTemp[detId]);
		      NN[3] += vSigma[detId] * (double)cimagf(quadTemp[detId]);
		    }
		    determinantMM = (double)AA*(double)CC-(double)BB*(double)BB;
		    if ( ( determinantMM*determinantMM < 1e-40 ) )
		      determinantMM = 1e-20; /* Note: Saving with positive sign */
		    InvMMAA = (double)CC / determinantMM;
		    InvMMBB = -(double)BB / determinantMM;
		    InvMMCC = (double)AA / determinantMM;
		    aa[0] = InvMMAA*NN[0] + InvMMBB*NN[1];
		    aa[1] = InvMMBB*NN[0] + InvMMCC*NN[1];
		    aa[2] = InvMMAA*NN[2] + InvMMBB*NN[3];
		    aa[3] = InvMMBB*NN[2] + InvMMCC*NN[3];

		    thisEvent->distance = XLALCoherentCBCParamEstim( &polarization,
		      &inclination, &coaPhase, aa[0], aa[1], aa[2], aa[3],
		      amplitudeConst, thisEvent,
		      (double) chirpTime,(double) crealf(quadTemp[0]),(double) cimagf(quadTemp[0]),
		      (double) crealf(quadTemp[1]),(double) cimagf(quadTemp[1]),(double) crealf(quadTemp[2]),
		      (double) cimagf(quadTemp[2]), (double) crealf(quadTemp[3]),
                      (double) cimagf(quadTemp[3]),sigmasq4DArray, case3b, caseID);

		    thisEvent->inclination = (REAL4) inclination;
		    thisEvent->polarization = (REAL4) polarization;
		    thisEvent->coa_phase = (REAL4) coaPhase;
		    thisEvent->ra = (REAL4) phi;
		    thisEvent->dec = (REAL4) theta;
		  }/*ends if condition on parameter-estimation or followup */

		  tempTime = 0.0;
		  fracpart = 0.0;
		  intpart = 0.0;
		  fflush( stdout );
		}
		else if ( timePt[0] > (eventStartIdx + deltaEventIndex) || !(params->maximizeOverChirp) ) {
		  /* clean up this event */
		  MultiInspiralTable      *lastEvent = NULL;

		  /* allocate memory for the newEvent */
		  lastEvent = thisEvent;
		  lastEvent->next = thisEvent = (MultiInspiralTable *)
		    LALCalloc( 1, sizeof(MultiInspiralTable) );
		  if ( !(lastEvent->next) )
		    {
		      ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
		    }

		  /* stick minimal data into the event */

		  tempTime = cData[0]->epoch.gpsSeconds + 1.0e-9 * cData[0]->epoch.gpsNanoSeconds + timePt[0] * deltaT;
		  fracpart = modf( tempTime, &intpart );
		  thisEvent->end_time.gpsSeconds = (INT4) intpart;
		  thisEvent->end_time.gpsNanoSeconds = (INT4) ( 1.0e9*fracpart );

		  /* Fill thisEvent with this trigger's quadrature values */
		  for ( detId=0 ; detId<params->numDetectors ; detId++ ) {
		    detIdSlidTimePt = timePt[0]+slidePoints[detId];
		    quadTemp[detId] = cData[detId]->data->data[detIdSlidTimePt];
		  }
		  XLALAssignEventQuads(caseID,quadTemp,thisEvent);
		  thisEvent->snr = cohSNR;
		  strcpy(thisEvent->ifos,caseStr);
		  thisEvent->mass1 = input->tmplt->mass1;
		  thisEvent->mass2 = input->tmplt->mass2;
		  thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
		  thisEvent->eta = input->tmplt->eta;
		  /* Save a set of alternative statistics */
		  crossCorrCohSq = XLALComputeCrossCohStatCase3b(cData,VVPlus,VVMinus,timePt[0],slidePoints,params->numDetectors,MM1,MM2,&autoCorrCohSq,chisq,chisq_dof,eff_snr_denom_fac,&crossCorrCohSqEff);
		  XLALSaveStatsInEventTable(cohSnrSqLocal,autoCorrCohSq,crossCorrCohSq,autoCorrNullSq,crossCorrNullSq, MM1, MM2, VVPlus, VVMinus, thisEvent);

		  /* Since not both H1 and H2 are present, the H1H2 null-stat is not meaningful*/
		  thisEvent->tau5 = eff_coh_snr;
                  thisEvent->tau0 = crossCorrCohSqEff;
                  /*Note: Temporarily save combined eff-snr here */
                  thisEvent->kappa = nullStatRegul;
                  if ( (MM1 == 0.0) || (MM2 == 0.0 ) ) {
                    thisEvent->null_stat_h1h2 = -2;
                  }
                  else {
                    thisEvent->null_stat_h1h2 = -1;
                  }

		  /* Compute network null-statistic at just trigger end-time */
                  if ( (MM1 == 0.0) || (MM2 == 0.0 ) ) {
		    thisEvent->null_stat_degen = -100.0;
		    thisEvent->null_statistic = thisEvent->null_stat_degen;
		  }
		  else {
		    thisEvent->null_statistic = (REAL4) XLALComputeNullStatCase3b(caseID,fplus,fcross,sigmasq,thisEvent);
		  }
		  /* Store the sigmasq's into event tables */
		  XLALCoherentCBCSigmasqFor3Sites(caseID,sigmasq4DArray,
				    chisq,thisEvent,case3b);
		  if ( params->estimParams || params->followup ) {
		    /* Parameter estimation: Distance
		       for ( detId=0 ; detId<params->numDetectors ; detId++ ) {*/
		    for ( i=0 ; i < 4 ; i++ ) {
		      NN[i] = 0.0;
		    }
		    for ( detId=0 ; detId < params->numDetectors ; detId++ ) {
		      NN[0] += uSigma[detId] * (double)crealf(quadTemp[detId]);
		      NN[1] += vSigma[detId] * (double)crealf(quadTemp[detId]);
		      NN[2] += uSigma[detId] * (double)cimagf(quadTemp[detId]);
		      NN[3] += vSigma[detId] * (double)cimagf(quadTemp[detId]);
		    }
		    determinantMM = (double)AA*(double)CC-(double)BB*(double)BB;
		    if ( ( determinantMM*determinantMM < 1e-40 ) )
		      determinantMM = 1e-20; /* Note: Saving with positive sign */
		    InvMMAA = (double)CC / determinantMM;
		    InvMMBB = -(double)BB / determinantMM;
		    InvMMCC = (double)AA / determinantMM;
		    aa[0] = InvMMAA*NN[0] + InvMMBB*NN[1];
		    aa[1] = InvMMBB*NN[0] + InvMMCC*NN[1];
		    aa[2] = InvMMAA*NN[2] + InvMMBB*NN[3];
		    aa[3] = InvMMBB*NN[2] + InvMMCC*NN[3];

		    thisEvent->distance = XLALCoherentCBCParamEstim( &polarization,
		      &inclination, &coaPhase, aa[0], aa[1], aa[2], aa[3],
		      amplitudeConst, thisEvent,
		      (double) chirpTime,(double) crealf(quadTemp[0]),(double) cimagf(quadTemp[0]),
		      (double) crealf(quadTemp[1]),(double) cimagf(quadTemp[1]),(double) crealf(quadTemp[2]),
                      (double) cimagf(quadTemp[2]), (double) crealf(quadTemp[3]),
                      (double) cimagf(quadTemp[3]),sigmasq4DArray, case3b, caseID);

		    thisEvent->inclination = (REAL4) inclination;
		    thisEvent->polarization = (REAL4) polarization;
		    thisEvent->coa_phase = (REAL4) coaPhase;
		    thisEvent->ra = (REAL4) phi;
		    thisEvent->dec = (REAL4) theta;
		  }/*ends if condition on parameter-estimation or followup */

		  /* Need to initialize the event start index to new value */
		  if( timePt[0] > (eventStartIdx + deltaEventIndex) )
		    {
		      eventStartIdx = timePt[0];
		    }
		  tempTime = 0.0;

		  fracpart= 0.0;
		  intpart = 0.0;
		} /* end of elseif statement */
	      } /* matches if (cohSNR > cohSNRThresh) */
	    } /* ends the else condition on timePt[0] */
	  } /* ends the "for loop" over the reference detector's timePt[0] */

	  skyGridIdx++;

	} /* ends loop over sky-position indices */
      } /* ends loop over case 3 */
    break;
  case 4: /* Network: 4 detectors */
    {
      LIGOTimeGPS     triggerGPSEndTime;/* Needed to calculate time-delays */
      double          psiInRadians = 0.0;
      double          detRefLocation[3];
      double          detNextLocation[3];
      const SkyGrid  *skyGrid;
      UINT4 skyGridIdx=0;

      /* This is case "4a", which pertains to a 4D network
	 with 4 ifos distributed at 3 different sites and with two
	 ifos sharing one of the sites, a la H1 and H2 */
      /* CHECK: This section not readied for using "newsnr" */
      /* However, the Case4a functions called here have been prepared for that */
      case4a = 1;

      triggerGPSEndTime.gpsSeconds = cData[0]->epoch.gpsSeconds;
      triggerGPSEndTime.gpsNanoSeconds = cData[0]->epoch.gpsNanoSeconds;
      /* Convert GPS time of trigger to GMST time in radians
	 for computing F+, Fx */
      gmstInRadians = XLALGreenwichMeanSiderealTime(&triggerGPSEndTime);

      /* Following needed because XLALArrivalTimeDiff() uses doubles */
      for ( locIdx=0 ; locIdx<3 ; locIdx++ ) {
	detRefLocation[locIdx] = (double) detectors[0].location[locIdx];
      }

      /* Loop over points in the sky-position grid */
      for ( skyGrid = skyGridPtr ; skyGrid ; skyGrid = skyGrid->next ) {
	phi = skyGrid->Alpha;
	theta = skyGrid->Delta;

	/* Loop over detectors computing their F+, Fx and t_c's */
	detId = 0;
	for( j=0; j<LAL_NUM_IFO; j++ ) {
	  /* Compute time-delays if caseID[j] != 0 */
	  if ( !(params->detIDVec->data[j] == 0 )) {
	    for ( locIdx=0 ; locIdx<3 ; locIdx++ ) {
	      detNextLocation[locIdx]
		= (double) detectors[detId].location[locIdx];
	    }

	    /*Arrival time delays in detectors relative to the 1st;
	      esp., dt21 = t2 - t1, dt31 = t3 - t1, etc. are computed */
	    /*	timeDelay[detId] = XLALArrivalTimeDiff(
		detectors[detId].location,
		detectors[0].location,
		raDec.longitude,raDec.latitude,&triggerGPSEndTime);
	    */
	    timeDelay[detId] = XLALArrivalTimeDiff(detNextLocation,
				  detRefLocation,phi,theta,&triggerGPSEndTime);

	    /* save unsorted time-delays for sorting after
	       this loop over detectors is exited */
	    sortedDelays4D[detId] = timeDelay[detId];

	    /* round off to nearest integer */
	    slidePoints4D[detId] = rint( timeDelay[detId]/deltaT );

	    detId++;
	  }
	}
	/* Sort the time-delays now */
	qsort( sortedDelays4D, (int) params->numDetectors, sizeof(double), compare );

	for ( detId=0 ; detId<params->numDetectors ; detId++ ) {
	  /* Compute the slide boundaries;
	     round off sorted time delays to nearest integer */
	  sortedSlidePoints4D[detId]
	    = rint( sortedDelays4D[detId]/deltaT );
	}

	/* Compute antenna-patterns and coherent SNR */
	/* Loop over detectors computing the theta-, phi-dependent
	   pieces of F+, Fx */
	detId = 0;
	for( j=0; j<LAL_NUM_IFO; j++ ) {
	  /* Compute antenna-patterns if caseID[j] != 0 */
	  if ( !(params->detIDVec->data[j] == 0 )) {
	    XLALComputeDetAMResponse(&fplus[detId], &fcross[detId],
				     (const REAL4(*)[3])detectors[detId].response, phi, theta,
				     psiInRadians, (double) gmstInRadians);

	    /* Compute antenna-pattern factors */
	    AAn[detId] = ( fplus[detId]) * ( fplus[detId]);
	    BBn[detId] = ( fplus[detId]) * ( fcross[detId]);
	    CCn[detId] = ( fcross[detId]) * ( fcross[detId]);

            /* Compute the elements of the helicity-plane projection matrix */
            AAn[detId] *= sigmasq[j];
            BBn[detId] *= sigmasq[j];
            CCn[detId] *= sigmasq[j];

	    /* Calculate factors necessary for parameter estimation*/
	    uSigma[detId] = (double)fplus[detId] * sqrt((double)sigmasq[j]);
	    vSigma[detId] = (double)fcross[detId] * sqrt((double)sigmasq[j]);

	    /* Sort sigmasq for degeneracy resoln. */
	    sortedSigmasq4D[detId] = sigmasq[j]/3;

	    detId++;
	  }
	}
	qsort( sortedSigmasq4D, (int) params->numDetectors, sizeof(double), compare );

	/* Construct network terms and factors required for
	   computing the coherent statistics */
	AA = AAn[0] + AAn[1] +AAn[2] + AAn[3];
	BB = BBn[0] + BBn[1] +BBn[2] + BBn[3];
	CC = CCn[0] + CCn[1] +CCn[2] + CCn[3];

	discrimSqrt = (AA*AA + 4*BB*BB
		       - 2*AA*CC + CC*CC);
	if ( (discrimSqrt>0.0) ) {
	  discrimSqrt = sqrt(discrimSqrt);
	}
	else {
	  discrimSqrt = 0.0;
	}

        O22 = ( -AA + CC + discrimSqrt);
        O22 /= (BB * sqrt( 4 + pow( AA - CC - discrimSqrt, 2)/( BB*BB ) ) );

        O21 = 1 / sqrt( 1 + pow(-AA+CC+discrimSqrt, 2)/ ( 4*BB*BB ) );

        O12 = ( -AA + CC - discrimSqrt);
        O12 /= ( BB * sqrt( 4 + pow(-AA+CC-discrimSqrt, 2) / ( BB*BB ) ) );

        O11 = 1 / sqrt( 1 + pow(-AA+CC-discrimSqrt, 2)/ ( 4*BB*BB ) );

        detId = 0;
        for( j=0; j<LAL_NUM_IFO; j++ ) {
          if ( !(params->detIDVec->data[j] == 0 )) {
            VVPlus[detId] = O11 * ( fplus[detId])
              + O12 * ( fcross[detId]);
            VVMinus[detId] = O21 * ( fplus[detId])
              + O22 * ( fcross[detId]);

            VVPlus[detId] *= sqrt((REAL4) sigmasq[j]);
            VVMinus[detId] *= sqrt((REAL4) sigmasq[j]);

            detId++;
          }
        }

        MM1 = 2*(AA*CC - BB*BB)/(AA+CC+discrimSqrt);
        MM2 = 2*(AA*CC - BB*BB)/(AA+CC-discrimSqrt);

	/* Factors for determining & computing degenerate case/statistics*/
	zetaRe = (AA-CC)/(AA+CC);
	zetaIm = 2*BB/(AA+CC);
	zetaAbsSq = zetaRe*zetaRe + zetaIm*zetaIm;
        degenerateStat = 0;
	if ( (pow(zetaAbsSq - zetaAbsSqMin,2) < 0.1) ) {
	  degenerateStat = 1;
	}
	for ( detId=0 ; detId<params->numDetectors ; detId++ ) {
	  DRe[detId] = (uSigma[detId] + zetaRe*uSigma[detId] +
			zetaIm*vSigma[detId]) / sqrt(AA+CC) ;
	  DIm[detId] = (vSigma[detId] - zetaRe*vSigma[detId] +
			zetaIm*uSigma[detId]) / sqrt(AA+CC) ;
	}

	/* Loop over time-points in the reference detector */
	for( timePt[0]=0 ; timePt[0]<(INT4)numPoints ; timePt[0]++) {
	  /* Reset cohSNR to zero so that it can only be ratcheted upward by
	     cohSnrSqLocal computed below for every point in sky-position grid*/
	  cohSNR = 0.0;

	  /* Loop over time-points in reference detector, after
	     accounting for the rounded-off sortedSlidePoints */
	  if( ( timePt[0] < (0-sortedSlidePoints4D[0]) )
	      ||  ( timePt[0]> (numPoints - 1 - sortedSlidePoints4D[3]) ) ) {
	    cohSNR = 0.0;
	    nullStatistic = 0.0;
	    autoCorrNullSq = 0.0;
	    ratioStatLocal = 0.0;
            if( cohSNROut ) {
              if ( !params->threeSiteCase ) {
                params->cohSNRVec->data->data[timePt[0]] = cohSNR;
              }
              else {
                params->cohSNRVec3Sites->data->data[skyGridIdx] = cohSNR;
              }
            }
            if( nullStatOut ) {
              if ( !params->threeSiteCase ) {
                params->nullStatVec->data->data[timePt[0]] = nullStatistic;
              }
              else {
                params->nullStatVec3Sites->data->data[skyGridIdx] = nullStatistic;
              }
            }
	  }
	  else {
            if ( (MM1 == 0.0) || (MM2 == 0.0 ) ) {
	      cohSnrSqLocal = 0.0;
	      ratioStatLocal = 0.0;
	    }
	    else {
              /* Note: This computation can be pushed to within threshcross loop */
	      cohSnrSqLocal = XLALComputeCohStatCase4a(cData,VVPlus,VVMinus,timePt[0],slidePoints4D,params->numDetectors,MM1,MM2,chisq,chisq_dof,eff_snr_denom_fac,&eff_coh_snr);
	      /* Note that unlike above, the next 3 are *squared* statistics */
	      nullStatistic = (REAL4) XLALComputeNonQuadNullStatCase4a(fplus,fcross,sigmasq,timePt[0],slidePoints4D,cData);

	      autoCorrNullSq = (REAL4) XLALComputeIncohStatCase4a(caseID,fplus,fcross,sigmasq,timePt[0],slidePoints4D,cData,nullStatistic,&crossCorrNullSq,chisq,chisq_dof,eff_snr_denom_fac);

	      ratioStatLocal = cohSnrSqLocal;

	    }
	    /* Updating RatioStat here */
	    if(ratioStatLocal > cohSNR) {
	      cohSNR = ratioStatLocal;
	      for ( detId=0 ; detId<params->numDetectors ; detId++ ) {
		detIdSlidTimePt = timePt[0]+slidePoints4D[detId];
		quadTemp[detId] = cData[detId]->data->data[detIdSlidTimePt];
	      }
              if( cohSNROut ) {
                if ( !params->threeSiteCase ) {
                  params->cohSNRVec->data->data[timePt[0]] = cohSNR;
                }
                else {
                  params->cohSNRVec3Sites->data->data[skyGridIdx] = cohSNR;
                }
              }
              if( nullStatOut ) {
                if ( !params->threeSiteCase ) {
                  params->nullStatVec->data->data[timePt[0]] = nullStatistic;
                }
                else {
                  params->nullStatVec3Sites->data->data[skyGridIdx] = nullStatistic;
                }
              }
	    }
	    /* Threshold on RatioStat here; cohSNR was set to RatioStat above */
	    if ( cohSNR > cohSNRThresh ) {
	      if ( !*eventList ) {
		/* store the start of the crossing */
		eventStartIdx = timePt[0];
		/* if this is the first event, start the list */
		thisEvent = *eventList = (MultiInspiralTable *)
		  LALCalloc( 1, sizeof(MultiInspiralTable) );

		if ( !thisEvent )
		  {
		    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
		  }
		/* record the data that we need for the clustering algorithm */
		tempTime = cData[0]->epoch.gpsSeconds + 1.0e-9 * cData[0]->epoch.gpsNanoSeconds + timePt[0] * deltaT;
		fracpart = modf( tempTime, &intpart );
		thisEvent->end_time.gpsSeconds = (INT4) intpart;
		thisEvent->end_time.gpsNanoSeconds = (INT4) ( 1.0e9*fracpart );

		/* Fill thisEvent with this trigger's quadrature values */
		for ( detId=0 ; detId<params->numDetectors ; detId++ ) {
		  detIdSlidTimePt = timePt[0]+slidePoints4D[detId];
		  quadTemp[detId] = cData[detId]->data->data[detIdSlidTimePt];
		}
		XLALAssignEventQuads(caseID,quadTemp,thisEvent);
		thisEvent->snr = cohSNR;
		cohSNRLocalRe = sqrt(sigmasq[1])*crealf(thisEvent->h1quad)
		  + sqrt(sigmasq[2])*crealf(thisEvent->h2quad);
		cohSNRLocalIm = sqrt(sigmasq[1])*cimagf(thisEvent->h1quad)
		  + sqrt(sigmasq[2])*cimagf(thisEvent->h2quad);
		thisEvent->coh_snr_h1h2 = (cohSNRLocalRe*cohSNRLocalRe +
			cohSNRLocalIm*cohSNRLocalIm) /(sigmasq[1] + sigmasq[2] );
		strcpy(thisEvent->ifos,caseStr);
		thisEvent->mass1 = input->tmplt->mass1;
		thisEvent->mass2 = input->tmplt->mass2;
		thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
		thisEvent->eta = input->tmplt->eta;
		/* Save a set of statistics */
		crossCorrCohSq = XLALComputeCrossCohStatCase4a(cData,VVPlus,VVMinus,timePt[0],slidePoints4D,params->numDetectors,MM1,MM2,&autoCorrCohSq,chisq,chisq_dof,eff_snr_denom_fac,&crossCorrCohSqEff);
		XLALSaveStatsInEventTable(cohSnrSqLocal,autoCorrCohSq,crossCorrCohSq,autoCorrNullSq,crossCorrNullSq, MM1, MM2, VVPlus, VVMinus, thisEvent);

		/* Compute null-statistic for H1-H2 at just trigger end-time */
		nullNorm = ( 1.0 / sigmasq[1]  + 1.0 /  sigmasq[2] );
		nullStatRe = crealf(thisEvent->h1quad) / sqrt(sigmasq[1])
		  - crealf(thisEvent->h2quad) / sqrt(sigmasq[2]);
		nullStatIm = cimagf(thisEvent->h1quad) / sqrt(sigmasq[1])
		  - cimagf(thisEvent->h2quad) / sqrt(sigmasq[2]);
		thisEvent->null_stat_h1h2 = ( nullStatRe*nullStatRe + nullStatIm*nullStatIm ) / nullNorm ;
		/* Compute network null-statistic at just trigger end-time */
		thisEvent->tau5 = eff_coh_snr;
                thisEvent->tau0 = crossCorrCohSqEff;
                /*Note: Temporarily save combined eff-snr here */
                thisEvent->kappa = nullStatRegul;
                if ( (MM1 == 0.0) || (MM2 == 0.0 ) ) {
		  thisEvent->null_stat_degen = -100.0;
		  thisEvent->null_statistic = thisEvent->null_stat_degen;
		}
		else {
		  if ( nullStatistic == 0.0 ) {
		    thisEvent->null_statistic = (REAL4) XLALComputeNullStatCase4a(fplus,fcross,sigmasq,thisEvent);
		  }
		  else {
		    thisEvent->null_statistic = (REAL4) nullStatistic;
		  }
		}
		/* Store the sigmasq's into event tables */
		XLALCoherentCBCSigmasqFor3Sites(caseID,sigmasq4DArray,
				    chisq,thisEvent,case3b);
		if ( params->estimParams || params->followup ) {
		  /* Parameter estimation: Distance
		     for ( detId=0 ; detId<params->numDetectors ; detId++ ) {*/
		  for ( i=0 ; i < 4 ; i++ ) {
		    NN[i] = 0.0;
		  }
		  for ( detId=0 ; detId < params->numDetectors ; detId++ ) {
		    NN[0] += uSigma[detId] * (double)crealf(quadTemp[detId]);
		    NN[1] += vSigma[detId] * (double)crealf(quadTemp[detId]);
		    NN[2] += uSigma[detId] * (double)cimagf(quadTemp[detId]);
		    NN[3] += vSigma[detId] * (double)cimagf(quadTemp[detId]);
		  }
		  determinantMM = (double)AA*(double)CC-(double)BB*(double)BB;
		  if ( ( determinantMM*determinantMM < 1e-40 ) )
		    determinantMM = 1e-20; /* Note: Saving with positive sign */
		  InvMMAA = (double)CC / determinantMM;
		  InvMMBB = -(double)BB / determinantMM;
		  InvMMCC = (double)AA / determinantMM;
		  aa[0] = InvMMAA*NN[0] + InvMMBB*NN[1];
		  aa[1] = InvMMBB*NN[0] + InvMMCC*NN[1];
		  aa[2] = InvMMAA*NN[2] + InvMMBB*NN[3];
		  aa[3] = InvMMBB*NN[2] + InvMMCC*NN[3];

		  thisEvent->distance = XLALCoherentCBCParamEstim( &polarization,
		    &inclination, &coaPhase, aa[0], aa[1], aa[2], aa[3],
		    amplitudeConst, thisEvent,
		    (double) chirpTime,(double) crealf(quadTemp[0]),(double) cimagf(quadTemp[0]),
		    (double) crealf(quadTemp[1]),(double) cimagf(quadTemp[1]),(double) crealf(quadTemp[2]),
		    (double) cimagf(quadTemp[2]), (double) crealf(quadTemp[3]),
                    (double) cimagf(quadTemp[3]), sigmasq4DArray, case3b, caseID);

		  thisEvent->inclination = (REAL4) inclination;
		  thisEvent->polarization = (REAL4) polarization;
		  thisEvent->coa_phase = (REAL4) coaPhase;
		  thisEvent->ra = (REAL4) phi;
		  thisEvent->dec = (REAL4) theta;
		}/*ends if condition on parameter-estimation or followup */

		tempTime = 0.0;
		fracpart = 0.0;
		intpart = 0.0;
		fflush( stdout );

	      } /* done creating a new event */
	      else if (params->maximizeOverChirp && timePt[0] <= (eventStartIdx + deltaEventIndex) && cohSNR > thisEvent->snr ) {
		/* if this is the same event, update the maximum */
		tempTime = cData[0]->epoch.gpsSeconds + 1.0e-9 * cData[0]->epoch.gpsNanoSeconds + timePt[0] * deltaT;
		fracpart = modf( tempTime, &intpart );
		thisEvent->end_time.gpsSeconds = (INT4) intpart;
		thisEvent->end_time.gpsNanoSeconds = (INT4) ( 1.0e9*fracpart );

		/* Fill thisEvent with this trigger's quadrature values */
		for ( detId=0 ; detId<params->numDetectors ; detId++ ) {
		  detIdSlidTimePt = timePt[0]+slidePoints4D[detId];
		  quadTemp[detId] = cData[detId]->data->data[detIdSlidTimePt];
		}
		XLALAssignEventQuads(caseID,quadTemp,thisEvent);
		thisEvent->snr = cohSNR;
		cohSNRLocalRe = sqrt(sigmasq[1])*crealf(thisEvent->h1quad)
		  + sqrt(sigmasq[2])*crealf(thisEvent->h2quad);
		cohSNRLocalIm = sqrt(sigmasq[1])*cimagf(thisEvent->h1quad)
		  + sqrt(sigmasq[2])*cimagf(thisEvent->h2quad);
		thisEvent->coh_snr_h1h2 = (cohSNRLocalRe*cohSNRLocalRe +
		  cohSNRLocalIm*cohSNRLocalIm) /(sigmasq[1] + sigmasq[2] );
		strcpy(thisEvent->ifos, caseStr);
		thisEvent->mass1 = input->tmplt->mass1;
		thisEvent->mass2 = input->tmplt->mass2;
		thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
		thisEvent->eta = input->tmplt->eta;
		/* Save a set of statistics */
		crossCorrCohSq = XLALComputeCrossCohStatCase4a(cData,VVPlus,VVMinus,timePt[0],slidePoints4D,params->numDetectors,MM1,MM2,&autoCorrCohSq,chisq,chisq_dof,eff_snr_denom_fac,&crossCorrCohSqEff);
		XLALSaveStatsInEventTable(cohSnrSqLocal,autoCorrCohSq,crossCorrCohSq,autoCorrNullSq,crossCorrNullSq, MM1, MM2, VVPlus, VVMinus, thisEvent);

		/* Compute null-statistic for H1-H2 at just trigger end-time */
		nullNorm = ( 1.0 / sigmasq[1]  + 1.0 /  sigmasq[2] );
		nullStatRe = crealf(thisEvent->h1quad) / sqrt(sigmasq[1])
		  - crealf(thisEvent->h2quad) / sqrt(sigmasq[2]);
		nullStatIm = cimagf(thisEvent->h1quad) / sqrt(sigmasq[1])
		  - cimagf(thisEvent->h2quad) / sqrt(sigmasq[2]);
		thisEvent->null_stat_h1h2 = ( nullStatRe*nullStatRe + nullStatIm*nullStatIm ) / nullNorm ;
		/* Compute network null-statistic at just trigger end-time */
		thisEvent->tau5 = eff_coh_snr;
                thisEvent->tau0 = crossCorrCohSqEff;
                /*Note: Temporarily save combined eff-snr here */
                thisEvent->kappa = nullStatRegul;
                if ( (MM1 == 0.0) || (MM2 == 0.0 ) ) {
		  thisEvent->null_stat_degen = -100.0;
		  thisEvent->null_statistic = thisEvent->null_stat_degen;
		}
		else {
		  if ( nullStatistic == 0.0 ) {
		    thisEvent->null_statistic = (REAL4) XLALComputeNullStatCase4a(fplus,fcross,sigmasq,thisEvent);
		  }
		  else {
		    thisEvent->null_statistic = (REAL4) nullStatistic;
		  }
		}
		/* Store the sigmasq's into event tables */
		XLALCoherentCBCSigmasqFor3Sites(caseID,sigmasq4DArray,
				    chisq,thisEvent,case3b);
		if ( params->estimParams || params->followup ) {
		  /* Parameter estimation: Distance
		     for ( detId=0 ; detId<params->numDetectors ; detId++ ) {*/
		  for ( i=0 ; i < 4 ; i++ ) {
		    NN[i] = 0.0;
		  }
		  for ( detId=0 ; detId < params->numDetectors ; detId++ ) {
		    NN[0] += uSigma[detId] * (double)crealf(quadTemp[detId]);
		    NN[1] += vSigma[detId] * (double)crealf(quadTemp[detId]);
		    NN[2] += uSigma[detId] * (double)cimagf(quadTemp[detId]);
		    NN[3] += vSigma[detId] * (double)cimagf(quadTemp[detId]);
		  }
		  determinantMM = (double)AA*(double)CC-(double)BB*(double)BB;
		  if ( ( determinantMM*determinantMM < 1e-40 ) )
		    determinantMM = 1e-20; /* Note: Saving with positive sign */
		  InvMMAA = (double)CC / determinantMM;
		  InvMMBB = -(double)BB / determinantMM;
		  InvMMCC = (double)AA / determinantMM;
		  aa[0] = InvMMAA*NN[0] + InvMMBB*NN[1];
		  aa[1] = InvMMBB*NN[0] + InvMMCC*NN[1];
		  aa[2] = InvMMAA*NN[2] + InvMMBB*NN[3];
		  aa[3] = InvMMBB*NN[2] + InvMMCC*NN[3];

		  thisEvent->distance = XLALCoherentCBCParamEstim( &polarization,
                    &inclination, &coaPhase, aa[0], aa[1], aa[2], aa[3],
                    amplitudeConst, thisEvent,
                    (double) chirpTime,(double) crealf(quadTemp[0]),(double) cimagf(quadTemp[0]),
                    (double) crealf(quadTemp[1]),(double) cimagf(quadTemp[1]),(double) crealf(quadTemp[2]),
		    (double) cimagf(quadTemp[2]), (double) crealf(quadTemp[3]),
                    (double) cimagf(quadTemp[3]), sigmasq4DArray, case3b, caseID);

		  thisEvent->inclination = (REAL4) inclination;
		  thisEvent->polarization = (REAL4) polarization;
		  thisEvent->coa_phase = (REAL4) coaPhase;
		  thisEvent->ra = (REAL4) phi;
		  thisEvent->dec = (REAL4) theta;
		}/*ends if condition on parameter-estimation or followup */

		tempTime = 0.0;
		fracpart = 0.0;
		intpart = 0.0;
		fflush( stdout );

	      }
	      else if ( timePt[0] > (eventStartIdx + deltaEventIndex) || !(params->maximizeOverChirp) ) {
		/* clean up this event */
		MultiInspiralTable      *lastEvent = NULL;

		/* allocate memory for the newEvent */
		lastEvent = thisEvent;
		lastEvent->next = thisEvent = (MultiInspiralTable *)
		  LALCalloc( 1, sizeof(MultiInspiralTable) );
		if ( !(lastEvent->next) )
		  {
		    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
		  }

		/* stick minimal data into the event */

		tempTime = cData[0]->epoch.gpsSeconds + 1.0e-9 * cData[0]->epoch.gpsNanoSeconds + timePt[0] * deltaT;
		fracpart = modf( tempTime, &intpart );
		thisEvent->end_time.gpsSeconds = (INT4) intpart;
		thisEvent->end_time.gpsNanoSeconds = (INT4) ( 1.0e9*fracpart );

		/* Fill thisEvent with this trigger's quadrature values */
		for ( detId=0 ; detId<params->numDetectors ; detId++ ) {
		  detIdSlidTimePt = timePt[0]+slidePoints4D[detId];
		  quadTemp[detId] = cData[detId]->data->data[detIdSlidTimePt];
		}
		XLALAssignEventQuads(caseID,quadTemp,thisEvent);
		thisEvent->snr = cohSNR;
		cohSNRLocalRe = sqrt(sigmasq[1])*crealf(thisEvent->h1quad)
		  + sqrt(sigmasq[2])*crealf(thisEvent->h2quad);
		cohSNRLocalIm = sqrt(sigmasq[1])*cimagf(thisEvent->h1quad)
		  + sqrt(sigmasq[2])*cimagf(thisEvent->h2quad);
		thisEvent->coh_snr_h1h2 = (cohSNRLocalRe*cohSNRLocalRe +
		  cohSNRLocalIm*cohSNRLocalIm) /(sigmasq[1] + sigmasq[2] );
		strcpy(thisEvent->ifos,caseStr);
		thisEvent->mass1 = input->tmplt->mass1;
		thisEvent->mass2 = input->tmplt->mass2;
		thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
		thisEvent->eta = input->tmplt->eta;
		/* Save a set of statistics */
		crossCorrCohSq = XLALComputeCrossCohStatCase4a(cData,VVPlus,VVMinus,timePt[0],slidePoints4D,params->numDetectors,MM1,MM2,&autoCorrCohSq,chisq,chisq_dof,eff_snr_denom_fac,&crossCorrCohSqEff);
		XLALSaveStatsInEventTable(cohSnrSqLocal,autoCorrCohSq,crossCorrCohSq,autoCorrNullSq,crossCorrNullSq, MM1, MM2, VVPlus, VVMinus, thisEvent);

		/* Compute null-statistic for H1-H2 at just trigger end-time */
		nullNorm = ( 1.0 / sigmasq[1]  + 1.0 /  sigmasq[2] );
		nullStatRe = crealf(thisEvent->h1quad) / sqrt(sigmasq[1])
		  - crealf(thisEvent->h2quad) / sqrt(sigmasq[2]);
		nullStatIm = cimagf(thisEvent->h1quad) / sqrt(sigmasq[1])
		  - cimagf(thisEvent->h2quad) / sqrt(sigmasq[2]);
		thisEvent->null_stat_h1h2 = ( nullStatRe*nullStatRe + nullStatIm*nullStatIm ) / nullNorm ;
		/* Compute network null-statistic at just trigger end-time */
		thisEvent->tau5 = eff_coh_snr;
                thisEvent->tau0 = crossCorrCohSqEff;
                /*Note: Temporarily save combined eff-snr here */
                thisEvent->kappa = nullStatRegul;
                if ( (MM1 == 0.0) || (MM2 == 0.0 ) ) {
		  thisEvent->null_stat_degen = -100.0;
		  thisEvent->null_statistic = thisEvent->null_stat_degen;
		}
		else {
		  if ( nullStatistic == 0.0 ) {
		    thisEvent->null_statistic = (REAL4) XLALComputeNullStatCase4a(fplus,fcross,sigmasq,thisEvent);
		  }
		  else {
		    thisEvent->null_statistic = (REAL4) nullStatistic;
		  }
		}
		inclination = 0.0;
		polarization = 0.0;
		/* Store the sigmasq's into event tables */
		XLALCoherentCBCSigmasqFor3Sites(caseID,sigmasq4DArray,
				    chisq,thisEvent,case3b);
		if ( params->estimParams || params->followup ) {
		  /* Parameter estimation: Distance
		     for ( detId=0 ; detId<params->numDetectors ; detId++ ) {*/
		  for ( i=0 ; i < 4 ; i++ ) {
		    NN[i] = 0.0;
		  }
		  for ( detId=0 ; detId < params->numDetectors ; detId++ ) {
		    NN[0] += uSigma[detId] * (double)crealf(quadTemp[detId]);
		    NN[1] += vSigma[detId] * (double)crealf(quadTemp[detId]);
		    NN[2] += uSigma[detId] * (double)cimagf(quadTemp[detId]);
		    NN[3] += vSigma[detId] * (double)cimagf(quadTemp[detId]);
		  }
		  determinantMM = (double)AA*(double)CC-(double)BB*(double)BB;
		  if ( ( determinantMM*determinantMM < 1e-40 ) )
		    determinantMM = 1e-20; /* Note: Saving with positive sign */
		  InvMMAA = (double)CC / determinantMM;
		  InvMMBB = -(double)BB / determinantMM;
		  InvMMCC = (double)AA / determinantMM;
		  aa[0] = InvMMAA*NN[0] + InvMMBB*NN[1];
		  aa[1] = InvMMBB*NN[0] + InvMMCC*NN[1];
		  aa[2] = InvMMAA*NN[2] + InvMMBB*NN[3];
		  aa[3] = InvMMBB*NN[2] + InvMMCC*NN[3];

		  thisEvent->distance = XLALCoherentCBCParamEstim( &polarization,
		    &inclination, &coaPhase, aa[0], aa[1], aa[2], aa[3],
		    amplitudeConst, thisEvent,
		    (double) chirpTime,(double) crealf(quadTemp[0]),(double) cimagf(quadTemp[0]),
		    (double) crealf(quadTemp[1]),(double) cimagf(quadTemp[1]),(double) crealf(quadTemp[2]),
		    (double) cimagf(quadTemp[2]), (double) crealf(quadTemp[3]),
                    (double) cimagf(quadTemp[3]), sigmasq4DArray, case3b, caseID);

		  thisEvent->inclination = (REAL4) inclination;
		  thisEvent->polarization = (REAL4) polarization;
		  thisEvent->coa_phase = (REAL4) coaPhase;
		  thisEvent->ra = (REAL4) phi;
		  thisEvent->dec = (REAL4) theta;
		}/*ends if condition on parameter-estimation or followup */

		/* Need to initialize the event start index to new value */
		if( timePt[0] > (eventStartIdx + deltaEventIndex) )
		  {
		    eventStartIdx = timePt[0];
		  }
		tempTime = 0.0;
		fracpart= 0.0;
		intpart = 0.0;
	      } /* end of elseif statement */
	    } /* matches if (cohSNR > cohSNRThresh) */
	  } /* ends the else condition on timePt[0] */
	} /* ends the "for loop" over the reference detector's timePt[0] */

        skyGridIdx++;

      } /* ends loop over sky-position indices */
    }/* case of 4 detectors */
  } /* closes  switch(params->numDetectors) */

  /* Compute time-series of the null-statistic of H1-H2 and of cohSNRH1H2
     if( thisEvent && params->nullStatOut && params->cohH1H2SNROut
      && !(case3b || case4a) ) {
  */
  if( thisEvent && params->nullStatH1H2Out && params->cohH1H2SNROut ) {
    /* Prepare norm for null statistic */
    nullNorm = ( 1.0 / sigmasq[1]  + 1.0 /  sigmasq[2] );

    /* Allocate memory for null statistic */
    memset( params->nullStatH1H2Vec->data->data, 0, numPoints*sizeof(REAL4));

    /* Allocate memory for cohSNRH1H2Vec if that SNR has
       not been computed above already*/
    memset( params->cohH1H2SNRVec->data->data, 0, numPoints*sizeof(REAL4));

    /*Note: Will not give intended result if first det is "G1", since it
      assumes that cdata[0] is H1 and cdata[1] is H1; rectify this in next rev. */
    for (k=0;k<(INT4)numPoints;k++) {

      /*Compute cohH1H2 snr, with index m replaced by k */
      cohSnrRe = sqrt(sigmasq[1])*crealf(cData[0]->data->data[k])
        + sqrt(sigmasq[2])*crealf(cData[1]->data->data[k]);
      cohSnrIm = sqrt(sigmasq[1])*cimagf(cData[0]->data->data[k])
        + sqrt(sigmasq[2])*cimagf(cData[1]->data->data[k]);

      params->cohH1H2SNRVec->data->data[k]
        = sqrt( (cohSnrRe*cohSnrRe + cohSnrIm*cohSnrIm) /
                (sigmasq[1] + sigmasq[2] ) );

      /* Compute null-stream statistic;
         in next rev. report re and im parts separately */
      nullStatRe = crealf(cData[0]->data->data[k]) / sqrt(sigmasq[1])
        - crealf(cData[1]->data->data[k]) / sqrt(sigmasq[2]);
      nullStatIm = cimagf(cData[0]->data->data[k]) / sqrt(sigmasq[1])
        - cimagf(cData[1]->data->data[k]) / sqrt(sigmasq[2]);
      params->nullStatH1H2Vec->data->data[k] =
        ( nullStatRe*nullStatRe + nullStatIm*nullStatIm ) / nullNorm ;
    }
  }

  /* Compute H1H2 null-statistic time-series and NOT cohH1H2SNR */
  if( thisEvent && params->nullStatH1H2Out && !(params->cohH1H2SNROut)
      && !(case3b || case4a) ) {

    /* Prepare norm for null statistic */
    nullNorm = ( 1.0 / sigmasq[1]  + 1.0 /  sigmasq[2] );

    /* Allocate memory for null statistic */
    memset( params->nullStatH1H2Vec->data->data, 0, numPoints*sizeof(REAL4));

    /*Note: Will not give intended result if first det is "G1", since it
      assumes that cdata[0] is H1 and cdata[1] is H2; rectify this in next rev. */
    for (k=0;k<(INT4)numPoints;k++) {
      /* Compute null-stream statistic;
         in next rev. report re and im parts separately */
      nullStatRe = crealf(cData[0]->data->data[k]) / sqrt(sigmasq[1])
        - crealf(cData[1]->data->data[k]) / sqrt(sigmasq[2]);
      nullStatIm = cimagf(cData[0]->data->data[k]) / sqrt(sigmasq[1])
        - cimagf(cData[1]->data->data[k]) / sqrt(sigmasq[2]);
      params->nullStatH1H2Vec->data->data[k] =
        ( nullStatRe*nullStatRe + nullStatIm*nullStatIm ) / nullNorm ;
    }
  }

  /* Compute cohSNRH1H2 ONLY, if not computed above already,
     but NOT the full  null-statistic time-series */
  if( thisEvent && params->cohH1H2SNROut && !nullStatH1H2Out ) {
    /* Allocate memory for cohSNRH1H2Vec if that SNR has
       not been computed above already*/
    memset( params->cohH1H2SNRVec->data->data, 0, numPoints*sizeof(REAL4));

    /*Note: Will not give intended result if first det is "G1", since it
      assumes that cdata[0] is H1 and cdata[1] is H1; rectify this in next rev. */
    for (k=0;k<(INT4)numPoints;k++) {

      /*Compute cohH1H2 snr, with index m replaced by k */
      cohSnrRe = sqrt(sigmasq[1])*crealf(cData[0]->data->data[k])
        + sqrt(sigmasq[2])*crealf(cData[1]->data->data[k]);
      cohSnrIm = sqrt(sigmasq[1])*cimagf(cData[0]->data->data[k])
        + sqrt(sigmasq[2])*cimagf(cData[1]->data->data[k]);

      params->cohH1H2SNRVec->data->data[k]
        = sqrt( (cohSnrRe*cohSnrRe + cohSnrIm*cohSnrIm) /
                (sigmasq[1] + sigmasq[2] ) );
    }
  }
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


void XLALAssignEventQuads(INT4 caseID[6],COMPLEX8 quadTemp[6],MultiInspiralTable *thisEvent) {
  INT4 found = 0;
  if(caseID[0])
    {
      thisEvent->g1quad = quadTemp[found];
      found += 1;
    }
  else
    {
      thisEvent->g1quad = 0.0;
    }
  if(caseID[1])
    {
      thisEvent->h1quad = quadTemp[found];
      found += 1;
    }
  else
    {
      thisEvent->h1quad = 0.0;
    }
  if(caseID[2])
    {
      thisEvent->h2quad = quadTemp[found];
      found += 1;
    }
  else
    {
      thisEvent->h2quad = 0.0;
    }
  if(caseID[3])
    {
      thisEvent->l1quad = quadTemp[found];
      found += 1;
    }
  else
    {
      thisEvent->l1quad = 0.0;
    }
  if(caseID[4])
    {
      thisEvent->t1quad = quadTemp[found];
      found += 1;
    }
  else
    {
      thisEvent->t1quad = 0.0;
    }
  if(caseID[5])
    {
      thisEvent->v1quad = quadTemp[found];
      found += 1;
    }
  else
    {
      thisEvent->v1quad = 0.0;
    }
}

void XLALSaveStatsInEventTable(REAL4 cohSnrSqLocal,REAL4 autoCorrCohSq,REAL4 crossCorrCohSq,REAL4 autoCorrNullSq,REAL4 crossCorrNullSq,REAL8 ampMetricEigenVal1,REAL8 ampMetricEigenVal2,REAL4 VVPlus[4], REAL4 VVMinus[4], MultiInspiralTable *thisEvent) {
  thisEvent->cohSnrSqLocal = cohSnrSqLocal;
  thisEvent->autoCorrCohSq = autoCorrCohSq;
  thisEvent->crossCorrCohSq = crossCorrCohSq;
  thisEvent->autoCorrNullSq = autoCorrNullSq;
  thisEvent->crossCorrNullSq = crossCorrNullSq;
  thisEvent->ampMetricEigenVal1 = ampMetricEigenVal1;
  thisEvent->ampMetricEigenVal2 = ampMetricEigenVal2;
  /* Note: Storage, only for H1L1V1 */
  thisEvent->tau2 = VVPlus[0];
  thisEvent->tau3 = VVPlus[1];
  thisEvent->tau4 = VVPlus[2];
  thisEvent->ttotal = VVMinus[0];
  thisEvent->bank_chisq = VVMinus[1];
  thisEvent->cont_chisq = VVMinus[2];
}

/* Function for computing 3-site-3-ifo coh-statistic at trigger end-time*/
REAL4 XLALComputeCohStatCase3b(COMPLEX8TimeSeries *cData[4], REAL4 VVPlus[4], REAL4 VVMinus[4], INT4 timeP, INT4 slidePoints[3], UINT4 numDetectors, REAL8 MM1, REAL8 MM2, REAL4 chisq[4], REAL4 chisq_dof[4], REAL4 eff_snr_denom_fac, REAL4 *eff_coh_snr) {

  /*Initialize cohSNR components and time stamps */
  UINT4          detId = 0;
  INT4           detIdSlidTimePt = 0;
  REAL4          CRePlus=0.0;
  REAL4          CImPlus=0.0;
  REAL4          CReMinus=0.0;
  REAL4          CImMinus=0.0;
  REAL4          CRePlusEff=0.0;
  REAL4          CImPlusEff=0.0;
  REAL4          cohSnrSqLocal=0.0;
  REAL4          eff_snr_denom_fac_inv = 1.0;

  eff_snr_denom_fac_inv = 1.0/eff_snr_denom_fac;

  /* Compute components of the coherent SNR */
  for ( detId=0 ; detId< (UINT4) numDetectors ; detId++ ) {
    REAL4          chisqFac = 1.0;
    REAL4          chisq_r  = 0.0;

    detIdSlidTimePt = timeP+slidePoints[detId];

    chisq_r = chisq[detId]/(2*chisq_dof[detId] -2);

    if ( chisq_r > 1.0 ) {
      /* Compute the denominator in new-snr for this ifo */
      chisqFac = pow( 0.5*(1 + pow(chisq_r, 0.5*eff_snr_denom_fac)), eff_snr_denom_fac_inv);
    }
    else {
      chisqFac = 1.0;
    }

    CRePlus += VVPlus[detId] *
      crealf(cData[detId]->data->data[detIdSlidTimePt]);
    CImPlus += VVPlus[detId] *
      cimagf(cData[detId]->data->data[detIdSlidTimePt]);
    CReMinus += VVMinus[detId] *
      crealf(cData[detId]->data->data[detIdSlidTimePt]);
    CImMinus += VVMinus[detId] *
      cimagf(cData[detId]->data->data[detIdSlidTimePt]);
    CRePlusEff += crealf(cData[detId]->data->data[detIdSlidTimePt])
      * crealf(cData[detId]->data->data[detIdSlidTimePt]) / pow(chisqFac,2);

    CImPlusEff += cimagf(cData[detId]->data->data[detIdSlidTimePt])
      * cimagf(cData[detId]->data->data[detIdSlidTimePt]) / pow(chisqFac,2);
  }

  /* Compute coherent SNR */
  cohSnrSqLocal = CRePlus*CRePlus/MM1 +
                  CImPlus*CImPlus/MM1 +
                  CReMinus*CReMinus/MM2 +
		  CImMinus*CImMinus/MM2;

  *eff_coh_snr = CRePlusEff + CImPlusEff;

  return cohSnrSqLocal;
}

/* Function for computing 3-site-3-ifo coh-statistic at trigger end-time*/
REAL4 XLALComputeCrossCohStatCase3b(COMPLEX8TimeSeries *cData[4], REAL4 VVPlus[4], REAL4 VVMinus[4], INT4 timeP, INT4 slidePoints[3], UINT4 numDetectors, REAL8 MM1, REAL8 MM2, REAL4 *autoCorrCohSq, REAL4 chisq[4], REAL4 chisq_dof[4], REAL4 eff_snr_denom_fac, REAL4 *crossCorrCohSqEff) {

  /*Initialize cohSNR components and time stamps */
  UINT4          detId = 0;
  INT4           detIdSlidTimePt = 0;
  REAL4          CRePlus=0.0;
  REAL4          CImPlus=0.0;
  REAL4          CReMinus=0.0;
  REAL4          CImMinus=0.0;
  REAL4          crossCorrCohSq=0.0;
  REAL4          Incohsq1=0.0;
  REAL4          Incohsq2=0.0;
  REAL4          Incohsq3=0.0;
  REAL4          Incohsq4=0.0;
  REAL4          autoCorrCohSqLocal = 0.0;
  REAL4          Incohsq1Eff=0.0;
  REAL4          Incohsq2Eff=0.0;
  REAL4          Incohsq3Eff=0.0;
  REAL4          Incohsq4Eff=0.0;
  REAL4          CRePlusEff=0.0;
  REAL4          CImPlusEff=0.0;
  REAL4          CReMinusEff=0.0;
  REAL4          CImMinusEff=0.0;
  REAL4          eff_snr_denom_fac_inv2 = 1.0;

  eff_snr_denom_fac_inv2 = 2.0/eff_snr_denom_fac;

  /* Compute components of the coherent SNR */
  for ( detId=0 ; detId< (UINT4) numDetectors ; detId++ ) {
    REAL4          chisqFacSq = 1.0;
    REAL4          chisq_r  = 0.0;

    detIdSlidTimePt = timeP+slidePoints[detId];

    chisq_r = chisq[detId]/(2*chisq_dof[detId] -2);

    if ( chisq_r > 1.0 ) {
      /* Compute the denominator in new-snr for this ifo */
      chisqFacSq = pow( 0.5*(1 + pow(chisq_r, 0.5*eff_snr_denom_fac)), eff_snr_denom_fac_inv2);
    }
    else {
      chisqFacSq = 1.0;
    }

    CRePlus += VVPlus[detId] *
      crealf(cData[detId]->data->data[detIdSlidTimePt]);
    CImPlus += VVPlus[detId] *
      cimagf(cData[detId]->data->data[detIdSlidTimePt]);
    CReMinus += VVMinus[detId] *
      crealf(cData[detId]->data->data[detIdSlidTimePt]);
    CImMinus += VVMinus[detId] *
      cimagf(cData[detId]->data->data[detIdSlidTimePt]);

    /* New lines for computing cross terms in coh stat */
    Incohsq1 += pow(VVPlus[detId] *
		   crealf(cData[detId]->data->data[detIdSlidTimePt]), 2);
    Incohsq2 += pow(VVPlus[detId] *
		    cimagf(cData[detId]->data->data[detIdSlidTimePt]), 2);
    Incohsq3 += pow(VVMinus[detId] *
		    crealf(cData[detId]->data->data[detIdSlidTimePt]), 2);
    Incohsq4 += pow(VVMinus[detId] *
		    cimagf(cData[detId]->data->data[detIdSlidTimePt]), 2);

    Incohsq1Eff += pow(VVPlus[detId] *
                    crealf(cData[detId]->data->data[detIdSlidTimePt]), 2) / chisqFacSq;
    Incohsq2Eff += pow(VVPlus[detId] *
                    cimagf(cData[detId]->data->data[detIdSlidTimePt]), 2) / chisqFacSq;
    Incohsq3Eff += pow(VVMinus[detId] *
                    crealf(cData[detId]->data->data[detIdSlidTimePt]), 2) / chisqFacSq;
    Incohsq4Eff += pow(VVMinus[detId] *
                    cimagf(cData[detId]->data->data[detIdSlidTimePt]), 2) / chisqFacSq;

    CRePlusEff += VVPlus[detId] *
      crealf(cData[detId]->data->data[detIdSlidTimePt]) / sqrt(chisqFacSq);
    CImPlusEff += VVPlus[detId] *
      cimagf(cData[detId]->data->data[detIdSlidTimePt]) / sqrt(chisqFacSq);
    CReMinusEff += VVMinus[detId] *
      crealf(cData[detId]->data->data[detIdSlidTimePt]) / sqrt(chisqFacSq);
    CImMinusEff += VVMinus[detId] *
      cimagf(cData[detId]->data->data[detIdSlidTimePt]) / sqrt(chisqFacSq);
  }

  /* Compute coherent SNR */
  crossCorrCohSq = ( CRePlus*CRePlus/MM1 +
			CImPlus*CImPlus/MM1 +
			CReMinus*CReMinus/MM2 +
			CImMinus*CImMinus/MM2 );
  /* Note that we are not taking the sqrt here, lest the quantity is negative*/
  autoCorrCohSqLocal = (Incohsq1 + Incohsq2)/MM1 + (Incohsq3 + Incohsq4)/MM2;
  crossCorrCohSq -= autoCorrCohSqLocal;

  *autoCorrCohSq = (Incohsq1Eff + Incohsq2Eff)/MM1 + (Incohsq3Eff + Incohsq4Eff)/MM2;

  *crossCorrCohSqEff = ( CRePlusEff*CRePlusEff/MM1 +
                        CImPlusEff*CImPlusEff/MM1 +
                        CReMinusEff*CReMinusEff/MM2 +
                        CImMinusEff*CImMinusEff/MM2 );

  *crossCorrCohSqEff -=  *autoCorrCohSq;

  return crossCorrCohSq;
}


/* Function for computing 3-site-3-ifo coh-statistic at trigger end-time*/
double XLALComputeNonQuadNullStatCase3b(INT4 caseID[6], double fplus[4], double fcross[4], REAL8 *sigmasq, INT4 timeP, INT4 slidePoints[3], COMPLEX8TimeSeries *cData[4]) {

  /* This trigger is from either H1 or H2 but not both */
  double sigmasqH = 0.0;
  double nullNorm8 = 0.0;
  double nullNumerRe8 = 0.0;
  double nullNumerIm8 = 0.0;
  double nullStatistic = 0.0;
  REAL4 Hquadre=0.0;
  REAL4 Hquadim=0.0;
  REAL4 l1quadre=0.0;
  REAL4 l1quadim=0.0;
  REAL4 v1quadre=0.0;
  REAL4 v1quadim=0.0;

  /* Compute components of the coherent SNR */
  Hquadre = crealf(cData[0]->data->data[timeP+slidePoints[0]]);
  Hquadim = cimagf(cData[0]->data->data[timeP+slidePoints[0]]);

  l1quadre = crealf(cData[1]->data->data[timeP+slidePoints[1]]);
  l1quadim = cimagf(cData[1]->data->data[timeP+slidePoints[1]]);
  v1quadre = crealf(cData[2]->data->data[timeP+slidePoints[2]]);
  v1quadim = cimagf(cData[2]->data->data[timeP+slidePoints[2]]);

  if ( caseID[1] == 0 ) {
      /* This is a H2 trigger */
      sigmasqH = sigmasq[2];

      nullNumerRe8 = (fplus[1]*fcross[2]-fplus[2]*fcross[1])*Hquadre/ sqrt(sigmasqH) +
        (fplus[2]*fcross[0]-fplus[0]*fcross[2])*l1quadre / sqrt(sigmasq[3]) +
        (fplus[0]*fcross[1]-fplus[1]*fcross[0])*v1quadre / sqrt(sigmasq[5]);

      nullNumerIm8 = (fplus[1]*fcross[2]-fplus[2]*fcross[1])*Hquadim / sqrt(sigmasqH) +
        (fplus[2]*fcross[0]-fplus[0]*fcross[2])*l1quadim / sqrt(sigmasq[3]) +
        (fplus[0]*fcross[1]-fplus[1]*fcross[0])*v1quadim / sqrt(sigmasq[5]);
    }
    else {
      sigmasqH = sigmasq[1];

      nullNumerRe8 = (fplus[1]*fcross[2]-fplus[2]*fcross[1])*Hquadre/ sqrt(sigmasqH) +
        (fplus[2]*fcross[0]-fplus[0]*fcross[2])*l1quadre / sqrt(sigmasq[3]) +
        (fplus[0]*fcross[1]-fplus[1]*fcross[0])*v1quadre / sqrt(sigmasq[5]);

      nullNumerIm8 = (fplus[1]*fcross[2]-fplus[2]*fcross[1])*Hquadim / sqrt(sigmasqH) +
        (fplus[2]*fcross[0]-fplus[0]*fcross[2])*l1quadim / sqrt(sigmasq[3]) +
        (fplus[0]*fcross[1]-fplus[1]*fcross[0])*v1quadim / sqrt(sigmasq[5]);
    }
    /* Prepare norm for null statistic */
    nullNorm8 = pow(fplus[1]*fcross[2]-fplus[2]*fcross[1],2)/ sigmasqH +
      pow(fplus[2]*fcross[0]-fplus[0]*fcross[2],2)/ sigmasq[3] +
      pow(fplus[0]*fcross[1]-fplus[1]*fcross[0],2)/ sigmasq[5] ;

    nullStatistic = ( nullNumerRe8*nullNumerRe8
                         + nullNumerIm8*nullNumerIm8)  / nullNorm8;
    return nullStatistic;
}


/* Function for computing 3-site-3-ifo coh-statistic at trigger end-time*/
double XLALComputeIncohStatCase3b(INT4 caseID[6], double fplus[4], double fcross[4], REAL8 *sigmasq, INT4 timeP, INT4 slidePoints[3], COMPLEX8TimeSeries *cData[4], REAL4 nullStatistic, REAL4 *crossCorrNullSq, REAL4 chisq[4], REAL4 chisq_dof[4], REAL4 eff_snr_denom_fac) {

  /* This trigger is from either H1 or H2 but not both */
  double sigmasqH = 0.0;
  double nullNorm8 = 0.0;
  double nullNumerSqRe8 = 0.0;
  double nullNumerSqIm8 = 0.0;
  double incohStat = 0.0;
  REAL4 Hquadre=0.0;
  REAL4 Hquadim=0.0;
  REAL4 l1quadre=0.0;
  REAL4 l1quadim=0.0;
  REAL4 v1quadre=0.0;
  REAL4 v1quadim=0.0;

  REAL4 HchisqFacSq = 1.0;
  REAL4 l1chisqFacSq = 1.0;
  REAL4 v1chisqFacSq = 1.0;

  REAL4 chisqH = 1.0;
  REAL4 chisq_dofH = 1.0;
  REAL4 eff_snr_denom_fac_inv2 = 1.0;

  eff_snr_denom_fac_inv2 = 2.0/eff_snr_denom_fac;

  /* Compute components of the coherent SNR */
  Hquadre = crealf(cData[0]->data->data[timeP+slidePoints[0]]);
  Hquadim = cimagf(cData[0]->data->data[timeP+slidePoints[0]]);

  l1quadre = crealf(cData[1]->data->data[timeP+slidePoints[1]]);
  l1quadim = cimagf(cData[1]->data->data[timeP+slidePoints[1]]);
  v1quadre = crealf(cData[2]->data->data[timeP+slidePoints[2]]);
  v1quadim = cimagf(cData[2]->data->data[timeP+slidePoints[2]]);

  if ( caseID[1] == 0 ) {
    /* This is a H2 trigger */
    sigmasqH = sigmasq[2];
    chisqH = chisq[0];/*Note: Correct for G1*/
    chisq_dofH = chisq_dof[0];
  }
  else {
    sigmasqH = sigmasq[1];
    chisqH = chisq[0];
    chisq_dofH = chisq_dof[0];
  }

  if ( (chisqH != 0.0) && ( chisqH/(2*chisq_dofH -2) > 1 ) ) {
    HchisqFacSq = pow( 0.5*(1 + pow(chisqH/(2*chisq_dofH -2), 0.5*eff_snr_denom_fac)), eff_snr_denom_fac_inv2);
  }
  else HchisqFacSq  = 1.0;

  if ( (chisq[1] != 0.0) && ( chisq[1]/(2*chisq_dof[1] -2) > 1 ) ) {
    l1chisqFacSq = pow( 0.5*(1 + pow(chisq[1]/(2*chisq_dof[1] -2), 0.5*eff_snr_denom_fac)), eff_snr_denom_fac_inv2);
  }
  else l1chisqFacSq = 1.0;

  if ( (chisq[2] != 0.0) && ( chisq[2]/(2*chisq_dof[2] -2) > 1 ) ) {
    v1chisqFacSq = pow( 0.5*(1 + pow(chisq[2]/(2*chisq_dof[2] -2), 0.5*eff_snr_denom_fac)), eff_snr_denom_fac_inv2);
  }
  else v1chisqFacSq = 1.0;

  nullNumerSqRe8 = pow((fplus[1]*fcross[2]-fplus[2]*fcross[1])*Hquadre/ sqrt(sigmasqH),2) +
    pow((fplus[2]*fcross[0]-fplus[0]*fcross[2])*l1quadre / sqrt(sigmasq[3]),2) +
    pow((fplus[0]*fcross[1]-fplus[1]*fcross[0])*v1quadre / sqrt(sigmasq[5]),2);

  nullNumerSqIm8 = pow((fplus[1]*fcross[2]-fplus[2]*fcross[1])*Hquadim / sqrt(sigmasqH),2) +
    pow((fplus[2]*fcross[0]-fplus[0]*fcross[2])*l1quadim / sqrt(sigmasq[3]),2) +
    pow((fplus[0]*fcross[1]-fplus[1]*fcross[0])*v1quadim / sqrt(sigmasq[5]),2);

  /* Prepare norm for null statistic */
  nullNorm8 = pow(fplus[1]*fcross[2]-fplus[2]*fcross[1],2)/ sigmasqH +
    pow(fplus[2]*fcross[0]-fplus[0]*fcross[2],2)/ sigmasq[3] +
    pow(fplus[0]*fcross[1]-fplus[1]*fcross[0],2)/ sigmasq[5] ;

  incohStat = ( nullNumerSqRe8 + nullNumerSqIm8)  / nullNorm8;

  *crossCorrNullSq = nullStatistic - incohStat;

  nullNumerSqRe8 = pow((fplus[1]*fcross[2]-fplus[2]*fcross[1])*Hquadre/ sqrt(sigmasqH),4)/HchisqFacSq +
    pow((fplus[2]*fcross[0]-fplus[0]*fcross[2])*l1quadre / sqrt(sigmasq[3]),4)/l1chisqFacSq +
    pow((fplus[0]*fcross[1]-fplus[1]*fcross[0])*v1quadre / sqrt(sigmasq[5]),4)/v1chisqFacSq;

  nullNumerSqIm8 = pow((fplus[1]*fcross[2]-fplus[2]*fcross[1])*Hquadim / sqrt(sigmasqH),4)/HchisqFacSq +
    pow((fplus[2]*fcross[0]-fplus[0]*fcross[2])*l1quadim / sqrt(sigmasq[3]),4)/l1chisqFacSq +
    pow((fplus[0]*fcross[1]-fplus[1]*fcross[0])*v1quadim / sqrt(sigmasq[5]),4)/v1chisqFacSq;

  nullNorm8 = pow(fplus[1]*fcross[2]-fplus[2]*fcross[1],4)/ pow(sigmasqH*HchisqFacSq*HchisqFacSq,2) +
    pow(fplus[2]*fcross[0]-fplus[0]*fcross[2],4)/ pow(sigmasq[3]*l1chisqFacSq*l1chisqFacSq,2) +
    pow(fplus[0]*fcross[1]-fplus[1]*fcross[0],4)/ pow(sigmasq[5]*v1chisqFacSq*v1chisqFacSq,2) ;

  incohStat = ( nullNumerSqRe8 + nullNumerSqIm8)  / nullNorm8;

  return incohStat;
}


/* Function for computing 3-site-3-ifo null-statistic at trigger end-time*/
double XLALComputeNullStatCase3b(INT4 caseID[6], double fplus[4], double fcross[4], REAL8 *sigmasq, MultiInspiralTable *thisEvent) {
    /* This trigger is from either H1 or H2 but not both */
    double sigmasqH = 0.0;
    double nullNorm8 = 0.0;
    double nullNumerRe8 = 0.0;
    double nullNumerIm8 = 0.0;
    double nullStatistic = 0.0;

    if ( caseID[1] == 0 ) {
      /* This is a H2 trigger */
      sigmasqH = sigmasq[2];

      nullNumerRe8 = (fplus[1]*fcross[2]-fplus[2]*fcross[1])*crealf(thisEvent->h2quad)/ sqrt(sigmasqH) +
        (fplus[2]*fcross[0]-fplus[0]*fcross[2])*crealf(thisEvent->l1quad) / sqrt(sigmasq[3]) +
        (fplus[0]*fcross[1]-fplus[1]*fcross[0])*crealf(thisEvent->v1quad) / sqrt(sigmasq[5]);

      nullNumerIm8 = (fplus[1]*fcross[2]-fplus[2]*fcross[1])*cimagf(thisEvent->h2quad) / sqrt(sigmasqH) +
        (fplus[2]*fcross[0]-fplus[0]*fcross[2])*cimagf(thisEvent->l1quad) / sqrt(sigmasq[3]) +
        (fplus[0]*fcross[1]-fplus[1]*fcross[0])*cimagf(thisEvent->v1quad) / sqrt(sigmasq[5]);
    }
    else {
      sigmasqH = sigmasq[1];

      nullNumerRe8 = (fplus[1]*fcross[2]-fplus[2]*fcross[1])*crealf(thisEvent->h1quad)/ sqrt(sigmasqH) +
        (fplus[2]*fcross[0]-fplus[0]*fcross[2])*crealf(thisEvent->l1quad) / sqrt(sigmasq[3]) +
        (fplus[0]*fcross[1]-fplus[1]*fcross[0])*crealf(thisEvent->v1quad) / sqrt(sigmasq[5]);

      nullNumerIm8 = (fplus[1]*fcross[2]-fplus[2]*fcross[1])*cimagf(thisEvent->h1quad) / sqrt(sigmasqH) +
        (fplus[2]*fcross[0]-fplus[0]*fcross[2])*cimagf(thisEvent->l1quad) / sqrt(sigmasq[3]) +
        (fplus[0]*fcross[1]-fplus[1]*fcross[0])*cimagf(thisEvent->v1quad) / sqrt(sigmasq[5]);
    }
    /* Prepare norm for null statistic */
    nullNorm8 = pow(fplus[1]*fcross[2]-fplus[2]*fcross[1],2)/ sigmasqH +
      pow(fplus[2]*fcross[0]-fplus[0]*fcross[2],2)/ sigmasq[3] +
      pow(fplus[0]*fcross[1]-fplus[1]*fcross[0],2)/ sigmasq[5] ;

    nullStatistic = ( nullNumerRe8*nullNumerRe8
                         + nullNumerIm8*nullNumerIm8)  / nullNorm8;
    return nullStatistic;
}


double XLALComputeNullTimeSeriesCase3b(INT4 caseID[6], double fplus[4], double fcross[4], REAL8 *sigmasq, COMPLEX8 quadTemp[6]) {
    /* This trigger is from either H1 or H2 but not both */
    double sigmasqH = 0.0;
    double nullNorm8 = 0.0;
    double nullNumerRe8 = 0.0;
    double nullNumerIm8 = 0.0;
    double nullStatistic = 0.0;

    if ( caseID[1] == 0 ) {
      /* This is a H2 trigger */
      sigmasqH = sigmasq[2];

      nullNorm8 = pow(fplus[1]*fcross[2]-fplus[2]*fcross[1],2)/ sigmasqH +
	pow(fplus[2]*fcross[0]-fplus[0]*fcross[2],2)/ sigmasq[3] +
	pow(fplus[0]*fcross[1]-fplus[1]*fcross[0],2)/ sigmasq[5] ;

      nullNumerRe8 = (fplus[1]*fcross[2]-fplus[2]*fcross[1])*crealf(quadTemp[0])/ sqrt(sigmasqH) +
        (fplus[2]*fcross[0]-fplus[0]*fcross[2])*crealf(quadTemp[1]) / sqrt(sigmasq[3]) +
        (fplus[0]*fcross[1]-fplus[1]*fcross[0])*crealf(quadTemp[2]) / sqrt(sigmasq[5]);

      nullNumerIm8 = (fplus[1]*fcross[2]-fplus[2]*fcross[1])*cimagf(quadTemp[0]) / sqrt(sigmasqH) +
        (fplus[2]*fcross[0]-fplus[0]*fcross[2])*cimagf(quadTemp[1]) / sqrt(sigmasq[3]) +
        (fplus[0]*fcross[1]-fplus[1]*fcross[0])*cimagf(quadTemp[2]) / sqrt(sigmasq[5]) ;
    }
    else {
      sigmasqH = sigmasq[1];

      nullNumerRe8 = (fplus[1]*fcross[2]-fplus[2]*fcross[1])*crealf(quadTemp[0])/ sqrt(sigmasqH) +
        (fplus[2]*fcross[0]-fplus[0]*fcross[2])*crealf(quadTemp[1]) / sqrt(sigmasq[3]) +
        (fplus[0]*fcross[1]-fplus[1]*fcross[0])*crealf(quadTemp[2]) / sqrt(sigmasq[5]);

      nullNumerIm8 = (fplus[1]*fcross[2]-fplus[2]*fcross[1])*cimagf(quadTemp[0]) / sqrt(sigmasqH) +
        (fplus[2]*fcross[0]-fplus[0]*fcross[2])*cimagf(quadTemp[1]) / sqrt(sigmasq[3]) +
        (fplus[0]*fcross[1]-fplus[1]*fcross[0])*cimagf(quadTemp[2]) / sqrt(sigmasq[5]) ;
    }
    /* Prepare norm for null statistic */
    nullNorm8 = pow(fplus[1]*fcross[2]-fplus[2]*fcross[1],2)/ sigmasqH +
      pow(fplus[2]*fcross[0]-fplus[0]*fcross[2],2)/ sigmasq[3] +
      pow(fplus[0]*fcross[1]-fplus[1]*fcross[0],2)/ sigmasq[5] ;

    nullStatistic = ( nullNumerRe8*nullNumerRe8
                         + nullNumerIm8*nullNumerIm8)  / nullNorm8;
    return nullStatistic;
}


/* Function for computing 3-site-3-ifo coh-statistic at trigger end-time*/
REAL4 XLALComputeCohStatCase4a(COMPLEX8TimeSeries *cData[4], REAL4 VVPlus[4], REAL4 VVMinus[4], INT4 timeP, INT4 slidePoints4D[4], UINT4 numDetectors, REAL8 MM1, REAL8 MM2, REAL4 chisq[4], REAL4 chisq_dof[4], REAL4 eff_snr_denom_fac, REAL4 *eff_coh_snr) {

  /*Initialize cohSNR components and time stamps */
  UINT4          detId = 0;
  INT4           detIdSlidTimePt = 0;
  REAL4          CRePlus=0.0;
  REAL4          CImPlus=0.0;
  REAL4          CReMinus=0.0;
  REAL4          CImMinus=0.0;
  REAL4          CRePlusEff=0.0;
  REAL4          CImPlusEff=0.0;
  REAL4          cohSnrSqLocal=0.0;
  REAL4          eff_snr_denom_fac_inv = 1.0;

  eff_snr_denom_fac_inv = 1.0/eff_snr_denom_fac;

  /* Compute components of the coherent SNR */
  for ( detId=0 ; detId<numDetectors ; detId++ ) {
    REAL4          chisqFac = 1.0;
    REAL4          chisq_r  = 0.0;

    detIdSlidTimePt = timeP+slidePoints4D[detId];

    chisq_r = chisq[detId]/(2*chisq_dof[detId] -2);

    if ( chisq_r > 1.0 ) {
      /* Compute the denominator in new-snr for this ifo */
      chisqFac = pow( 0.5*(1 + pow(chisq_r, 0.5*eff_snr_denom_fac)), eff_snr_denom_fac_inv);
    }
    else {
      chisqFac = 1.0;
    }

    CRePlus += VVPlus[detId] *
      crealf(cData[detId]->data->data[detIdSlidTimePt]);
    CImPlus += VVPlus[detId] *
      cimagf(cData[detId]->data->data[detIdSlidTimePt]);
    CReMinus += VVMinus[detId] *
      crealf(cData[detId]->data->data[detIdSlidTimePt]);
    CImMinus += VVMinus[detId] *
      cimagf(cData[detId]->data->data[detIdSlidTimePt]);

    CRePlusEff += crealf(cData[detId]->data->data[detIdSlidTimePt])
      * crealf(cData[detId]->data->data[detIdSlidTimePt]) / pow(chisqFac,2);

    CImPlusEff += cimagf(cData[detId]->data->data[detIdSlidTimePt])
      * cimagf(cData[detId]->data->data[detIdSlidTimePt]) / pow(chisqFac,2);
  }

  /* Compute coherent SNR */
  cohSnrSqLocal = CRePlus*CRePlus/MM1 +
                  CImPlus*CImPlus/MM1 +
                  CReMinus*CReMinus/MM2 +
		  CImMinus*CImMinus/MM2;

  *eff_coh_snr = CRePlusEff + CImPlusEff;

  return cohSnrSqLocal;
}


/* Function for computing 3-site-3-ifo coh-statistic at trigger end-time*/
REAL4 XLALComputeCrossCohStatCase4a(COMPLEX8TimeSeries *cData[4], REAL4 VVPlus[4], REAL4 VVMinus[4], INT4 timeP, INT4 slidePoints4D[4], UINT4 numDetectors, REAL8 MM1, REAL8 MM2, REAL4 *autoCorrCohSq, REAL4 chisq[4], REAL4 chisq_dof[4], REAL4 eff_snr_denom_fac, REAL4 *crossCorrCohSqEff) {

  /*Initialize cohSNR components and time stamps */
  UINT4          detId = 0;
  INT4           detIdSlidTimePt = 0;
  REAL4          CRePlus=0.0;
  REAL4          CImPlus=0.0;
  REAL4          CReMinus=0.0;
  REAL4          CImMinus=0.0;
  REAL4          crossCorrCohSq=0.0;
  REAL4          Incohsq1=0.0;
  REAL4          Incohsq2=0.0;
  REAL4          Incohsq3=0.0;
  REAL4          Incohsq4=0.0;
  REAL4          autoCorrCohSqLocal = 0.0;
  REAL4          Incohsq1Eff=0.0;
  REAL4          Incohsq2Eff=0.0;
  REAL4          Incohsq3Eff=0.0;
  REAL4          Incohsq4Eff=0.0;
  REAL4          CRePlusEff=0.0;
  REAL4          CImPlusEff=0.0;
  REAL4          CReMinusEff=0.0;
  REAL4          CImMinusEff=0.0;
  REAL4          eff_snr_denom_fac_inv2 = 1.0;

  eff_snr_denom_fac_inv2 = 2.0/eff_snr_denom_fac;

  /* Compute components of the coherent SNR */
  for ( detId=0 ; detId< (UINT4) numDetectors ; detId++ ) {
    REAL4          chisqFacSq = 1.0;
    REAL4          chisq_r  = 0.0;

    detIdSlidTimePt = timeP+slidePoints4D[detId];

    chisq_r = chisq[detId]/(2*chisq_dof[detId] -2);

    if ( chisq_r > 1.0 ) {
      /* Compute the denominator in new-snr for this ifo */
      chisqFacSq = pow( 0.5*(1 + pow(chisq_r, 0.5*eff_snr_denom_fac)), eff_snr_denom_fac_inv2);
    }
    else {
      chisqFacSq = 1.0;
    }

    CRePlus += VVPlus[detId] *
      crealf(cData[detId]->data->data[detIdSlidTimePt]);
    CImPlus += VVPlus[detId] *
      cimagf(cData[detId]->data->data[detIdSlidTimePt]);
    CReMinus += VVMinus[detId] *
      crealf(cData[detId]->data->data[detIdSlidTimePt]);
    CImMinus += VVMinus[detId] *
      cimagf(cData[detId]->data->data[detIdSlidTimePt]);

    /* New lines for computing cross terms in coh stat */
    Incohsq1 += pow(VVPlus[detId] *
		   crealf(cData[detId]->data->data[detIdSlidTimePt]), 2);
    Incohsq2 += pow(VVPlus[detId] *
		    cimagf(cData[detId]->data->data[detIdSlidTimePt]), 2);
    Incohsq3 += pow(VVMinus[detId] *
		    crealf(cData[detId]->data->data[detIdSlidTimePt]), 2);
    Incohsq4 += pow(VVMinus[detId] *
		    cimagf(cData[detId]->data->data[detIdSlidTimePt]), 2);

    Incohsq1Eff += pow(VVPlus[detId] *
                    crealf(cData[detId]->data->data[detIdSlidTimePt]), 2) / chisqFacSq;
    Incohsq2Eff += pow(VVPlus[detId] *
                    cimagf(cData[detId]->data->data[detIdSlidTimePt]), 2) / chisqFacSq;
    Incohsq3Eff += pow(VVMinus[detId] *
                    crealf(cData[detId]->data->data[detIdSlidTimePt]), 2) / chisqFacSq;
    Incohsq4Eff += pow(VVMinus[detId] *
                    cimagf(cData[detId]->data->data[detIdSlidTimePt]), 2) / chisqFacSq;

    CRePlusEff += VVPlus[detId] *
      crealf(cData[detId]->data->data[detIdSlidTimePt]) / sqrt(chisqFacSq);
    CImPlusEff += VVPlus[detId] *
      cimagf(cData[detId]->data->data[detIdSlidTimePt]) / sqrt(chisqFacSq);
    CReMinusEff += VVMinus[detId] *
      crealf(cData[detId]->data->data[detIdSlidTimePt]) / sqrt(chisqFacSq);
    CImMinusEff += VVMinus[detId] *
      cimagf(cData[detId]->data->data[detIdSlidTimePt]) / sqrt(chisqFacSq);
  }

  /* Compute coherent SNR */
  crossCorrCohSq = ( CRePlus*CRePlus/MM1 +
			CImPlus*CImPlus/MM1 +
			CReMinus*CReMinus/MM2 +
			CImMinus*CImMinus/MM2 );
  /*Note that we are not taking the sqrt here, lest the quantity is negative*/
  autoCorrCohSqLocal = (Incohsq1 + Incohsq2)/MM1 + (Incohsq3 + Incohsq4)/MM2;
  crossCorrCohSq -= autoCorrCohSqLocal;

  *autoCorrCohSq = (Incohsq1Eff + Incohsq2Eff)/MM1 + (Incohsq3Eff + Incohsq4Eff)/MM2;

  *crossCorrCohSqEff = ( CRePlusEff*CRePlusEff/MM1 +
                        CImPlusEff*CImPlusEff/MM1 +
                        CReMinusEff*CReMinusEff/MM2 +
                        CImMinusEff*CImMinusEff/MM2 );

  *crossCorrCohSqEff -=  *autoCorrCohSq;

  return crossCorrCohSq;
}


/* Function for computing 3-site-3-ifo coh-statistic at trigger end-time*/
double XLALComputeNonQuadNullStatCase4a(double fplus[4], double fcross[4], REAL8 *sigmasq, INT4 timeP, INT4 slidePoints4D[4], COMPLEX8TimeSeries *cData[4]) {

  /* This trigger is from either H1 or H2 but not both */
  double nullNorm8 = 0.0;
  double nullNumerRe8 = 0.0;
  double nullNumerIm8 = 0.0;
  double nullStatistic = 0.0;
  REAL4 Hquadre=0.0;
  REAL4 Hquadim=0.0;
  REAL4 l1quadre=0.0;
  REAL4 l1quadim=0.0;
  REAL4 v1quadre=0.0;
  REAL4 v1quadim=0.0;

  /* Compute components of the coherent SNR */
  Hquadre = crealf(cData[0]->data->data[timeP+slidePoints4D[0]]);
  Hquadim = cimagf(cData[0]->data->data[timeP+slidePoints4D[0]]);

  l1quadre = crealf(cData[2]->data->data[timeP+slidePoints4D[2]]);
  l1quadim = cimagf(cData[2]->data->data[timeP+slidePoints4D[2]]);
  v1quadre = crealf(cData[3]->data->data[timeP+slidePoints4D[3]]);
  v1quadim = cimagf(cData[3]->data->data[timeP+slidePoints4D[3]]);

  nullNumerRe8 = (fplus[1]*fcross[2]-fplus[2]*fcross[1])*Hquadre/ sqrt(sigmasq[1]) +
    (fplus[2]*fcross[0]-fplus[0]*fcross[2])*l1quadre / sqrt(sigmasq[3]) +
    (fplus[0]*fcross[1]-fplus[1]*fcross[0])*v1quadre / sqrt(sigmasq[5]);

  nullNumerIm8 = (fplus[1]*fcross[2]-fplus[2]*fcross[1])*Hquadim / sqrt(sigmasq[1]) +
    (fplus[2]*fcross[0]-fplus[0]*fcross[2])*l1quadim / sqrt(sigmasq[3]) +
    (fplus[0]*fcross[1]-fplus[1]*fcross[0])*v1quadim / sqrt(sigmasq[5]);

  /* Prepare norm for null statistic */
  nullNorm8 = pow(fplus[1]*fcross[2]-fplus[2]*fcross[1],2)/ sigmasq[1] +
    pow(fplus[2]*fcross[0]-fplus[0]*fcross[2],2)/ sigmasq[3] +
    pow(fplus[0]*fcross[1]-fplus[1]*fcross[0],2)/ sigmasq[5] ;

  nullStatistic = ( nullNumerRe8*nullNumerRe8
		    + nullNumerIm8*nullNumerIm8)  / nullNorm8;
  return nullStatistic;
}


/* Function for computing 3-site-3-ifo coh-statistic at trigger end-time*/
double XLALComputeIncohStatCase4a(INT4 UNUSED caseID[6], double fplus[4], double fcross[4], REAL8 *sigmasq, INT4 timeP, INT4 slidePoints4D[4], COMPLEX8TimeSeries *cData[4], REAL4 nullStatistic, REAL4 *crossCorrNullSq, REAL4 chisq[4], REAL4 chisq_dof[4], REAL4 eff_snr_denom_fac) {
  /* caseID is unused in this function */

  /* This trigger is from either H1 or H2 but not both */
  double nullNorm8 = 0.0;
  double nullNumerSqRe8 = 0.0;
  double nullNumerSqIm8 = 0.0;
  double incohStat = 0.0;
  REAL4 Hquadre=0.0;
  REAL4 Hquadim=0.0;
  REAL4 l1quadre=0.0;
  REAL4 l1quadim=0.0;
  REAL4 v1quadre=0.0;
  REAL4 v1quadim=0.0;

  REAL4 HchisqFacSq = 1.0;
  REAL4 l1chisqFacSq = 1.0;
  REAL4 v1chisqFacSq = 1.0;
  REAL4 eff_snr_denom_fac_inv2 = 1.0;

  eff_snr_denom_fac_inv2 = 2.0/eff_snr_denom_fac;

  Hquadre = crealf(cData[0]->data->data[timeP+slidePoints4D[0]]);
  Hquadim = cimagf(cData[0]->data->data[timeP+slidePoints4D[0]]);

  l1quadre = crealf(cData[2]->data->data[timeP+slidePoints4D[2]]);
  l1quadim = cimagf(cData[2]->data->data[timeP+slidePoints4D[2]]);
  v1quadre = crealf(cData[3]->data->data[timeP+slidePoints4D[3]]);
  v1quadim = cimagf(cData[3]->data->data[timeP+slidePoints4D[3]]);

  if ( (chisq[0] != 0.0) && ( chisq[0]/(2*chisq_dof[0] -2) > 1 ) ) {
    HchisqFacSq = pow( 0.5*(1 + pow(chisq[0]/(2*chisq_dof[0] -2), 0.5*eff_snr_denom_fac)), eff_snr_denom_fac_inv2);
  }
  else HchisqFacSq  = 1.0;

  if ( (chisq[2] != 0.0) && ( chisq[2]/(2*chisq_dof[2] -2) > 1 ) ) {
    l1chisqFacSq = pow( 0.5*(1 + pow(chisq[2]/(2*chisq_dof[2] -2), 0.5*eff_snr_denom_fac)), eff_snr_denom_fac_inv2);
  }
  else l1chisqFacSq = 1.0;

  if ( (chisq[3] != 0.0) && ( chisq[3]/(2*chisq_dof[3] -2) > 1 ) ) {
    v1chisqFacSq = pow( 0.5*(1 + pow(chisq[3]/(2*chisq_dof[3] -2), 0.5*eff_snr_denom_fac)), eff_snr_denom_fac_inv2);
  }
  else v1chisqFacSq = 1.0;

  /* Compute components of the coherent SNR */
  nullNumerSqRe8 = pow((fplus[1]*fcross[2]-fplus[2]*fcross[1])*Hquadre/ sqrt(sigmasq[1]),2) +
    pow((fplus[2]*fcross[0]-fplus[0]*fcross[2])*l1quadre / sqrt(sigmasq[3]),2) +
    pow((fplus[0]*fcross[1]-fplus[1]*fcross[0])*v1quadre / sqrt(sigmasq[5]),2);

  nullNumerSqIm8 = pow((fplus[1]*fcross[2]-fplus[2]*fcross[1])*Hquadim / sqrt(sigmasq[1]),2) +
    pow((fplus[2]*fcross[0]-fplus[0]*fcross[2])*l1quadim / sqrt(sigmasq[3]),2) +
    pow((fplus[0]*fcross[1]-fplus[1]*fcross[0])*v1quadim / sqrt(sigmasq[5]),2);

  /* Prepare norm for null statistic */
  nullNorm8 = pow(fplus[1]*fcross[2]-fplus[2]*fcross[1],2)/ sigmasq[1] +
    pow(fplus[2]*fcross[0]-fplus[0]*fcross[2],2)/ sigmasq[3] +
    pow(fplus[0]*fcross[1]-fplus[1]*fcross[0],2)/ sigmasq[5] ;

  incohStat = ( nullNumerSqRe8 + nullNumerSqIm8)  / nullNorm8;

  *crossCorrNullSq = nullStatistic - incohStat;

  nullNumerSqRe8 = pow((fplus[1]*fcross[2]-fplus[2]*fcross[1])*Hquadre/ sqrt(sigmasq[1]),4)/HchisqFacSq +
    pow((fplus[2]*fcross[0]-fplus[0]*fcross[2])*l1quadre / sqrt(sigmasq[3]),4)/l1chisqFacSq +
    pow((fplus[0]*fcross[1]-fplus[1]*fcross[0])*v1quadre / sqrt(sigmasq[5]),4)/v1chisqFacSq;

  nullNumerSqIm8 = pow((fplus[1]*fcross[2]-fplus[2]*fcross[1])*Hquadim / sqrt(sigmasq[1]),4)/HchisqFacSq +
    pow((fplus[2]*fcross[0]-fplus[0]*fcross[2])*l1quadim / sqrt(sigmasq[3]),4)/l1chisqFacSq +
    pow((fplus[0]*fcross[1]-fplus[1]*fcross[0])*v1quadim / sqrt(sigmasq[5]),4)/v1chisqFacSq;

  /* Prepare norm for null statistic */
  nullNorm8 = pow(fplus[1]*fcross[2]-fplus[2]*fcross[1],4)/ pow(sigmasq[1]*HchisqFacSq*HchisqFacSq,2) +
    pow(fplus[2]*fcross[0]-fplus[0]*fcross[2],4)/ pow(sigmasq[3]*l1chisqFacSq*l1chisqFacSq,2) +
    pow(fplus[0]*fcross[1]-fplus[1]*fcross[0],4)/ pow(sigmasq[5]*v1chisqFacSq*v1chisqFacSq,2);

  incohStat = ( nullNumerSqRe8 + nullNumerSqIm8)  / nullNorm8;

  return incohStat;
}


/* Function for computing 3-site-4-ifo null-statistic at trigger end-time*/
double XLALComputeNullStatCase4a(double fplus[4], double fcross[4], REAL8 *sigmasq, MultiInspiralTable *thisEvent) {
    /* This trigger is from both H1 and H2;
     but using H1 and not H2 for now*/
    double nullNorm8 = 0.0;
    double nullNumerRe8 = 0.0;
    double nullNumerIm8 = 0.0;
    double nullStatistic = 0.0;

    /* Prepare norm for null statistic */
    nullNorm8 = pow(fplus[1]*fcross[2]-fplus[2]*fcross[1],2)/ sigmasq[1] +
      pow(fplus[2]*fcross[0]-fplus[0]*fcross[2],2)/ sigmasq[3] +
      pow(fplus[0]*fcross[1]-fplus[1]*fcross[0],2)/ sigmasq[5] ;

    nullNumerRe8 = (fplus[1]*fcross[2]-fplus[2]*fcross[1])*crealf(thisEvent->h1quad) / sqrt(sigmasq[1]) +
      (fplus[2]*fcross[0]-fplus[0]*fcross[2])*crealf(thisEvent->l1quad) / sqrt(sigmasq[3]) +
      (fplus[0]*fcross[1]-fplus[1]*fcross[0])*crealf(thisEvent->v1quad) / sqrt(sigmasq[5]);

    nullNumerIm8 = (fplus[1]*fcross[2]-fplus[2]*fcross[1])*cimagf(thisEvent->h1quad) / sqrt(sigmasq[1]) +
      (fplus[2]*fcross[0]-fplus[0]*fcross[2])*cimagf(thisEvent->l1quad) / sqrt(sigmasq[3]) +
      (fplus[0]*fcross[1]-fplus[1]*fcross[0])*cimagf(thisEvent->v1quad) / sqrt(sigmasq[5]) ;

    nullStatistic = ( nullNumerRe8*nullNumerRe8
                         + nullNumerIm8*nullNumerIm8)  / nullNorm8;

    return nullStatistic;
}

double XLALComputeNullTimeSeriesCase4a(INT4 UNUSED caseID[6], double fplus[4], double fcross[4], REAL8 *sigmasq, COMPLEX8 quadTemp[6]) {
    /* This trigger is from both H1 and H2;
     but using H1 and not H2 for now*/
    double nullNorm8 = 0.0;
    double nullNumerRe8 = 0.0;
    double nullNumerIm8 = 0.0;
    double nullStatistic = 0.0;

    /* Prepare norm for null statistic */
    nullNorm8 = pow(fplus[1]*fcross[2]-fplus[2]*fcross[1],2)/ sigmasq[1] +
      pow(fplus[2]*fcross[0]-fplus[0]*fcross[2],2)/ sigmasq[3] +
      pow(fplus[0]*fcross[1]-fplus[1]*fcross[0],2)/ sigmasq[5] ;

    nullNumerRe8 = (fplus[1]*fcross[2]-fplus[2]*fcross[1])*crealf(quadTemp[0]) / sqrt(sigmasq[1]) +
      (fplus[2]*fcross[0]-fplus[0]*fcross[2])*crealf(quadTemp[2]) / sqrt(sigmasq[3]) +
      (fplus[0]*fcross[1]-fplus[1]*fcross[0])*crealf(quadTemp[3]) / sqrt(sigmasq[5]);

    nullNumerIm8 = (fplus[1]*fcross[2]-fplus[2]*fcross[1])*cimagf(quadTemp[0]) / sqrt(sigmasq[1]) +
      (fplus[2]*fcross[0]-fplus[0]*fcross[2])*cimagf(quadTemp[2]) / sqrt(sigmasq[3]) +
      (fplus[0]*fcross[1]-fplus[1]*fcross[0])*cimagf(quadTemp[3]) / sqrt(sigmasq[5]) ;

    nullStatistic = ( nullNumerRe8*nullNumerRe8
                         + nullNumerIm8*nullNumerIm8)  / nullNorm8;


    return nullStatistic;
}


int compare( const void* a, const void* b ) {
  const double* arg1 = (const double*) a;
  const double* arg2 = (const double*) b;
  if( *arg1 < *arg2 ) return -1;
  else if( *arg1 == *arg2 ) return 0;
  else return 1;
}

void XLALCoherentCBCEstimateDistanceCase2a(double C_Real0, double C_Im0, double C_Real1, double C_Im1, REAL8 sigmasq[4], REAL4 chisq[4], MultiInspiralTable *thisEvent) {
 /* Sigmasq includes chirpTime! */
 thisEvent->eff_dist_h1 = sqrt(sigmasq[0])/pow((C_Real0*C_Real0+C_Im0*C_Im0),.5);
 thisEvent->eff_dist_h2 = sqrt(sigmasq[1])/pow((C_Real1*C_Real1+C_Im1*C_Im1),.5);

 /*Note: Need to update the code here so that the C's correspond to H1 and H2 */
 thisEvent->eff_dist_h1h2 = (sigmasq[0]+
     sigmasq[1])/sqrt(sigmasq[0]*(C_Real0*C_Real0+C_Im0*C_Im0)+
     sigmasq[1]*(C_Real1*C_Real1+C_Im1*C_Im1)+
     2*sqrt(sigmasq[0]*sigmasq[1])*(C_Real0*C_Real1+C_Im0*C_Im1));

 thisEvent->sigmasq_h1 = sigmasq[0];
 thisEvent->sigmasq_h2 = sigmasq[1];

 thisEvent->chisq_h1 = chisq[0];
 thisEvent->chisq_h2 = chisq[1];

 thisEvent->distance = -1;
 thisEvent->inclination = -1001;
 thisEvent->polarization = -1001;
 thisEvent->ra = -1001;
 thisEvent->dec = -1001;
}

void XLALCoherentCBCEstimateDistanceCase2b(INT4 caseID[6], double C_Real0, double C_Im0, double C_Real1, double C_Im1, REAL8 sigmasq[4], REAL4 chisq[4], MultiInspiralTable *thisEvent) {
  UINT4 det1Found = 0;

  thisEvent->chisq_g = chisq[2];
  thisEvent->chisq_t = chisq[3];

  thisEvent->distance = -1;
  thisEvent->eff_dist_g = -1;
  thisEvent->eff_dist_h1 = -1;
  thisEvent->eff_dist_h2 = -1;
  thisEvent->eff_dist_l = -1;
  thisEvent->eff_dist_t = -1;
  thisEvent->eff_dist_v = -1;
  thisEvent->eff_dist_h1h2 = -1;
  thisEvent->inclination = -1001;
  thisEvent->polarization = -1001;
  thisEvent->ra = -1001;
  thisEvent->dec = -1001;
  /* Sigmasq includes chirpTime! */
  if( caseID[0] ) {
    thisEvent->eff_dist_g = sqrt(sigmasq[0])/pow((C_Real0*C_Real0+C_Im0*C_Im0),.5);
    thisEvent->sigmasq_g = sigmasq[0];
    thisEvent->chisq_g = chisq[0];
    det1Found =1;
  }
  if( caseID[1] ) {
    if ( det1Found ) {
      thisEvent->eff_dist_h1 = sqrt(sigmasq[1])/pow((C_Real1*C_Real1+C_Im1*C_Im1),.5);
      thisEvent->sigmasq_h1 = sigmasq[1];
      thisEvent->chisq_h1 = chisq[1];
    }
    else {
      thisEvent->eff_dist_h1 = sqrt(sigmasq[0])/pow((C_Real0*C_Real0+C_Im0*C_Im0),.5);
      thisEvent->sigmasq_h1 = sigmasq[0];
      thisEvent->chisq_h1 = chisq[0];
      det1Found =1;
    }
  }
  if( caseID[2] ) {
    if ( det1Found ) {
      thisEvent->eff_dist_h2 = sqrt(sigmasq[1])/pow((C_Real1*C_Real1+C_Im1*C_Im1),.5);
      thisEvent->sigmasq_h2 = sigmasq[1];
      thisEvent->chisq_h2 = chisq[1];
    }
    else {
      thisEvent->eff_dist_h2 = sqrt(sigmasq[0])/pow((C_Real0*C_Real0+C_Im0*C_Im0),.5);
      thisEvent->sigmasq_h2 = sigmasq[0];
      thisEvent->chisq_h2 = chisq[0];
      det1Found =1;
    }
  }
  if( caseID[3] ) {
    if ( det1Found ) {
      thisEvent->eff_dist_l = sqrt(sigmasq[1])/pow((C_Real1*C_Real1+C_Im1*C_Im1),.5);
      thisEvent->sigmasq_l = sigmasq[1];
      thisEvent->chisq_l = chisq[1];
    }
    else {
      thisEvent->eff_dist_l = sqrt(sigmasq[0])/pow((C_Real0*C_Real0+C_Im0*C_Im0),.5);
      thisEvent->sigmasq_l = sigmasq[0];
      thisEvent->chisq_l = chisq[0];
      det1Found =1;
    }
  }
  if( caseID[4] ) {
    if ( det1Found ) {
      thisEvent->eff_dist_t = sqrt(sigmasq[1])/pow((C_Real1*C_Real1+C_Im1*C_Im1),.5);
      thisEvent->sigmasq_t = sigmasq[1];
      thisEvent->chisq_t = chisq[1];
    }
    else {
      thisEvent->eff_dist_t = sqrt(sigmasq[0])/pow((C_Real0*C_Real0+C_Im0*C_Im0),.5);
      thisEvent->sigmasq_t = sigmasq[0];
      thisEvent->chisq_t = chisq[0];
      det1Found =1;
    }
  }
  if( caseID[5] ) {
    /*If it has come this far, v1 must be detector 2 */
    thisEvent->eff_dist_v = sqrt(sigmasq[1])/pow((C_Real1*C_Real1+C_Im1*C_Im1),.5);
    thisEvent->sigmasq_v = sigmasq[1];
    thisEvent->chisq_v = chisq[1];
  }
}

void XLALCoherentCBCEstimateDistanceCase3a(INT4 caseID[6], double C_Real0, double C_Im0, double C_Real1,double C_Im1, double C_Real2, double C_Im2, REAL8 sigmasq[4], REAL4 chisq[4], MultiInspiralTable *thisEvent) {
  thisEvent->distance = -1;
  thisEvent->eff_dist_g = -1;
  thisEvent->eff_dist_l = -1;
  thisEvent->eff_dist_t = -1;
  thisEvent->eff_dist_v = -1;
  thisEvent->inclination = -1001;
  thisEvent->polarization = -1001;
  thisEvent->ra = -1001;
  thisEvent->dec = -1001;

 /* Sigmasq includes chirpTime! */
  if( caseID[0] ) {
    thisEvent->eff_dist_g = sqrt(sigmasq[0])/pow((C_Real0*C_Real0+C_Im0*C_Im0),.5);
    thisEvent->eff_dist_h1 = sqrt(sigmasq[1])/pow((C_Real1*C_Real1+C_Im1*C_Im1),.5);
    thisEvent->eff_dist_h2 = sqrt(sigmasq[2])/pow((C_Real2*C_Real2+C_Im2*C_Im2),.5);
    thisEvent->eff_dist_h1h2 = (sigmasq[1]+
     sigmasq[2])/sqrt(sigmasq[1]*(C_Real1*C_Real1+C_Im1*C_Im1)+
     sigmasq[2]*(C_Real2*C_Real2+C_Im2*C_Im2)+
     2*sqrt(sigmasq[1]*sigmasq[2])*(C_Real1*C_Real2+C_Im1*C_Im2));

    thisEvent->sigmasq_g = sigmasq[0];
    thisEvent->sigmasq_h1 = sigmasq[1];
    thisEvent->sigmasq_h2 = sigmasq[2];
    thisEvent->chisq_g = chisq[0];
    thisEvent->chisq_h1 = chisq[1];
    thisEvent->chisq_h2 = chisq[2];
  }
  else {
    thisEvent->eff_dist_h1 = sqrt(sigmasq[0])/pow((C_Real0*C_Real0+C_Im0*C_Im0),.5);
    thisEvent->eff_dist_h2 = sqrt(sigmasq[1])/pow((C_Real1*C_Real1+C_Im1*C_Im1),.5);
    thisEvent->eff_dist_h1h2 = (sigmasq[0]+
     sigmasq[1])/sqrt(sigmasq[0]*(C_Real0*C_Real0+C_Im0*C_Im0)+
     sigmasq[1]*(C_Real1*C_Real1+C_Im1*C_Im1)+
     2*sqrt(sigmasq[0]*sigmasq[1])*(C_Real0*C_Real1+C_Im0*C_Im1));

    thisEvent->sigmasq_h1 = sigmasq[0];
    thisEvent->sigmasq_h2 = sigmasq[1];

    thisEvent->chisq_h1 = chisq[0];
    thisEvent->chisq_h2 = chisq[1];

    if( caseID[3] ) {
      thisEvent->eff_dist_l = sqrt(sigmasq[2])/pow((C_Real2*C_Real2+C_Im2*C_Im2),.5);
      thisEvent->sigmasq_l = sigmasq[2];
      thisEvent->chisq_l = chisq[2];
    }
    else if( caseID[4] ) {
      thisEvent->eff_dist_t = sqrt(sigmasq[2])/pow((C_Real2*C_Real2+C_Im2*C_Im2),.5);
      thisEvent->sigmasq_t = sigmasq[2];
      thisEvent->chisq_t = chisq[2];
    }
    else if( caseID[5] ) {
      thisEvent->eff_dist_v = sqrt(sigmasq[2])/pow((C_Real2*C_Real2+C_Im2*C_Im2),.5);
      thisEvent->sigmasq_v = sigmasq[2];
      thisEvent->chisq_v = chisq[2];
    }
  }
}


void XLALCoherentCBCSigmasqFor3Sites(INT4 caseID[6], REAL8 sigmasq[4], REAL4 chisq[4], MultiInspiralTable *thisEvent, UINT4 case3b) {
  UINT4   det1Found = 0;
  UINT4   det2Found = 0;

  if ( case3b ) {
    /* Sigmasq includes chirpTime! */
    if( caseID[0] ) {
      thisEvent->sigmasq_g = sigmasq[0];
      thisEvent->chisq_g = chisq[0];
      det1Found =1;
    }
    if( caseID[1] ) {
      if ( det1Found ) {
	thisEvent->sigmasq_h1 = sigmasq[1];
	thisEvent->chisq_h1 = chisq[1];
	det2Found = 1;
      }
      else {
	thisEvent->sigmasq_h1 = sigmasq[0];
	thisEvent->chisq_h1 = chisq[0];
	det1Found =1;
      }
    }
    if( caseID[2] ) {
      if ( det1Found ) {
	thisEvent->sigmasq_h2 = sigmasq[1];
	thisEvent->chisq_h2 = chisq[1];
	det2Found = 1;
      }
      else {
	thisEvent->sigmasq_h2 = sigmasq[0];
	thisEvent->chisq_h2 = chisq[0];
	det1Found =1;
      }
    }
    if( caseID[3] ) {
      if ( det1Found && !det2Found ) {
	thisEvent->sigmasq_l = sigmasq[1];
	thisEvent->chisq_l = chisq[1];
	det2Found = 1;
      }
      else if( det1Found && det2Found) {
	thisEvent->sigmasq_l = sigmasq[2];
	thisEvent->chisq_l = chisq[2];
      }
      else {
	thisEvent->sigmasq_l = sigmasq[0];
	thisEvent->chisq_l = chisq[0];
	det1Found =1;
      }
    }
    if( caseID[4] ) {
      if ( det1Found && !det2Found  ) {
	/*If it has come this far, either this is detector 2 or... */
	thisEvent->sigmasq_t = sigmasq[1];
	thisEvent->chisq_t = chisq[1];
	det2Found = 1;
      }
      else {
	/* ... this is detector 3 */
	thisEvent->sigmasq_t = sigmasq[2];
	thisEvent->chisq_t = chisq[2];
      }
    }
    else if( caseID[5] ) {
      /*If it has come this far, v1 must be detector 3 */
      thisEvent->sigmasq_v = sigmasq[2];
      thisEvent->chisq_v = chisq[2];
    }
  }
  else {
    /* This is for "case4a = 1", and hence has both h1 and h2. */
    if( caseID[0] ) {
      det1Found = 1; /* One of the two detectors other than h1 or h2 found */
      thisEvent->sigmasq_g = sigmasq[0];
      thisEvent->sigmasq_h1 = sigmasq[1];
      thisEvent->sigmasq_h2 = sigmasq[2];
      thisEvent->chisq_g = chisq[0];
      thisEvent->chisq_h1 = chisq[1];
      thisEvent->chisq_h2 = chisq[2];
    }
    else {
      thisEvent->sigmasq_h1 = sigmasq[0];
      thisEvent->sigmasq_h2 = sigmasq[1];
      thisEvent->chisq_h1 = chisq[0];
      thisEvent->chisq_h2 = chisq[1];
    }
    if( caseID[3] ) {
      if( det1Found ) {
	thisEvent->sigmasq_l = sigmasq[3];
	thisEvent->chisq_l = chisq[3];
      }
      else {
	thisEvent->sigmasq_l = sigmasq[2];
	thisEvent->chisq_l = chisq[2];
	det1Found = 1;
      }
    }
    if( caseID[4] ) {
      if( det1Found ) {
	thisEvent->sigmasq_t = sigmasq[3];
	thisEvent->chisq_t = chisq[3];
      }
      else {
	thisEvent->sigmasq_t = sigmasq[2];
	thisEvent->chisq_t = chisq[2];
      }
    }
    if( caseID[5] ) {
      /*If it has come this far, v1 must be detector 3 */
      thisEvent->sigmasq_v = sigmasq[3];
      thisEvent->chisq_v = chisq[3];
    }
  }
}

double XLALCoherentCBCParamEstim( double *psi_est, double *iota_est, double *coa_phase_est, double a1, double a2, double a3, double a4, double amplitudeConst, MultiInspiralTable *thisEvent, double UNUSED chirpTime, double C_Real0, double C_Im0, double C_Real1,double C_Im1, double C_Real2, double C_Im2, double C_Real3, double C_Im3, REAL8 sigmasq[4], UINT4 case3b,INT4 caseID[6] ) {

  double lum_dist_est,f_a,f_a_sq,g_a,h_a,sine_psi=0.0,sine_coa_phase=0.0,p=0.0;
  UINT4 det1Found = 0;
  UINT4 det2Found = 0;

  thisEvent->eff_dist_g = -1;
  thisEvent->eff_dist_h1 = -1;
  thisEvent->eff_dist_h2 = -1;
  thisEvent->eff_dist_l = -1;
  thisEvent->eff_dist_t = -1;
  thisEvent->eff_dist_v = -1;
  thisEvent->eff_dist_h1h2 = -1;

  p = a1*a4-a2*a3;
  f_a=2*p/(a1*a1+a2*a2+a3*a3+a4*a4);
  f_a_sq=f_a*f_a;

  if(f_a_sq>=1.0)
    {
      f_a_sq = 1.0;
    }

  g_a = f_a/(1+pow((1-f_a_sq),.5));


  h_a=(1/g_a)-pow(((1/(g_a*g_a))-1),.5);

  if(h_a<-1.0)
    {
      h_a=(1/g_a)+pow(((1/(g_a*g_a)-1)),.5);
    }


 *iota_est=acos(h_a);

 /*   LUMINOSITY DISTANCE  */
 lum_dist_est = pow((1+6*h_a*h_a+h_a*h_a*h_a*h_a)/(a1*a1+a2*a2+a3*a3+a4*a4),0.5);
 lum_dist_est *= amplitudeConst;

 /*   Phase Angle, coa_phase_est  */


 sine_coa_phase = 2.0*lum_dist_est*lum_dist_est*(a1*a3+a2*a4)/((1+h_a*h_a)*(4*h_a*h_a/pow((1+h_a*h_a),2)-1));

 if(sine_coa_phase>1.0||sine_coa_phase<-1.0)
   {
     *coa_phase_est = -50;
   }
 else
   {
     *coa_phase_est=0.5*asin(2.0*lum_dist_est*lum_dist_est*(a1*a3+a2*a4)/((1+h_a*h_a)*(4*h_a*h_a/pow((1+h_a*h_a),2)-1)));

     if(*coa_phase_est<0)
       {
         *coa_phase_est=2*3.1415+*coa_phase_est;
       }

   }

 /*   Polarization Angle    */

 if((((float)a1==(float)a4)&&((float)a2==-(float)a3))||(((float)a1==-(float)a4)&&((float)a2==(float)a3)))
   {
     *psi_est = -50.;
   }
 else
   {
     sine_psi = 2.*(a1*a3+a2*a4)/pow((((a1-a4)*(a1-a4)+(a2+a3)*(a2+a3))*((a1+a4)*(a1+a4)+(a2-a3)*(a2-a3))),0.5);
     if(sine_psi>1.0||sine_psi<-1.0)
       {
         *psi_est = -50.;
       }
     else
       {
         *psi_est = 0.25*asin(sine_psi);
         if(*psi_est<0)
           {
             *psi_est=2*3.1415+*coa_phase_est;
           }

       }
   }

 if ( case3b ) {
   /* Sigmasq includes chirpTime! */
   if( caseID[0] ) {
     thisEvent->eff_dist_g = sqrt(sigmasq[0])/pow((C_Real0*C_Real0+C_Im0*C_Im0),.5);
     det1Found =1;
   }
   if( caseID[1] ) {
     if ( det1Found ) {
       thisEvent->eff_dist_h1 = sqrt(sigmasq[1])/pow((C_Real1*C_Real1+C_Im1*C_Im1),.5);
       det2Found = 1;
     }
     else {
       thisEvent->eff_dist_h1 = sqrt(sigmasq[0])/pow((C_Real0*C_Real0+C_Im0*C_Im0),.5);
       det1Found =1;
     }
   }
   if( caseID[2] ) {
     if ( det1Found ) {
       thisEvent->eff_dist_h2 = sqrt(sigmasq[1])/pow((C_Real1*C_Real1+C_Im1*C_Im1),.5);
       det2Found = 1;
     }
     else {
       thisEvent->eff_dist_h2 = sqrt(sigmasq[0])/pow((C_Real0*C_Real0+C_Im0*C_Im0),.5);
       det1Found =1;
     }
   }
   if( caseID[3] ) {
     if ( det1Found && !det2Found ) {
       thisEvent->eff_dist_l = sqrt(sigmasq[1])/pow((C_Real1*C_Real1+C_Im1*C_Im1),.5);
       det2Found = 1;
     }
     else if( det1Found && det2Found) {
       thisEvent->eff_dist_l = sqrt(sigmasq[2])/pow((C_Real2*C_Real2+C_Im2*C_Im2),.5);
     }
     else {
       thisEvent->eff_dist_l = sqrt(sigmasq[0])/pow((C_Real0*C_Real0+C_Im0*C_Im0),.5);
       det1Found =1;
     }
   }
   if( caseID[4] ) {
     if ( det1Found && !det2Found  ) {
       /*If it has come this far, either this is detector 2 or... */
       thisEvent->eff_dist_t = sqrt(sigmasq[1])/pow((C_Real1*C_Real1+C_Im1*C_Im1),.5);
       det2Found = 1;
     }
     else {
       /* ... this is detector 3 */
       thisEvent->eff_dist_t = sqrt(sigmasq[2])/pow((C_Real2*C_Real2+C_Im2*C_Im2),.5);
     }
   }
   else if( caseID[5] ) {
     /*If it has come this far, v1 must be detector 3 */
     thisEvent->eff_dist_v = sqrt(sigmasq[2])/pow((C_Real2*C_Real2+C_Im2*C_Im2),.5);
   }
 }
 else {
   /* This is for "case4a = 1", and hence has both h1 and h2. */
   if( caseID[0] ) {
     det1Found = 1; /* One of the two detectors other than h1 or h2 found */
     thisEvent->eff_dist_g = sqrt(sigmasq[0])/pow((C_Real0*C_Real0+C_Im0*C_Im0),.5);

     thisEvent->eff_dist_h1 = sqrt(sigmasq[1])/pow((C_Real1*C_Real1+C_Im1*C_Im1),.5);
     thisEvent->eff_dist_h2 = sqrt(sigmasq[2])/pow((C_Real2*C_Real2+C_Im2*C_Im2),.5);
     thisEvent->eff_dist_h1h2 = (sigmasq[1]+
       sigmasq[2])/sqrt(sigmasq[1]*(C_Real1*C_Real1+C_Im1*C_Im1)+
       sigmasq[2]*(C_Real2*C_Real2+C_Im2*C_Im2)+
       2*sqrt(sigmasq[1]*sigmasq[2])*(C_Real1*C_Real2+C_Im1*C_Im2));
   }
   else {
     thisEvent->eff_dist_h1 = sqrt(sigmasq[0])/pow((C_Real0*C_Real0+C_Im0*C_Im0),.5);
     thisEvent->eff_dist_h2 = sqrt(sigmasq[1])/pow((C_Real1*C_Real1+C_Im1*C_Im1),.5);
     thisEvent->eff_dist_h1h2 = (sigmasq[0]+
       sigmasq[1])/sqrt(sigmasq[0]*(C_Real0*C_Real0+C_Im0*C_Im0)+
       sigmasq[1]*(C_Real1*C_Real1+C_Im1*C_Im1)+
       2*sqrt(sigmasq[0]*sigmasq[1])*(C_Real0*C_Real1+C_Im0*C_Im1));
   }
   if( caseID[3] ) {
     if( det1Found ) {
      thisEvent->eff_dist_l = sqrt(sigmasq[3])/pow((C_Real3*C_Real3+C_Im3*C_Im3),.5);
     }
     else {
       thisEvent->eff_dist_l = sqrt(sigmasq[2])/pow((C_Real2*C_Real2+C_Im2*C_Im2),.5);
       det1Found = 1;
     }
   }
   if( caseID[4] ) {
     if( det1Found ) {
       thisEvent->eff_dist_t = sqrt(sigmasq[3])/pow((C_Real3*C_Real3+C_Im3*C_Im3),.5);
     }
     else {
       thisEvent->eff_dist_t = sqrt(sigmasq[2])/pow((C_Real2*C_Real2+C_Im2*C_Im2),.5);
     }
   }
   if( caseID[5] ) {
     /*If it has come this far, v1 must be detector 3 */
     thisEvent->eff_dist_v = sqrt(sigmasq[3])/pow((C_Real3*C_Real3+C_Im3*C_Im3),.5);
   }
 }

 return lum_dist_est;
}
