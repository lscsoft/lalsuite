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
 * Revision: $Id$
 * 
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
#include <lal/DataBuffer.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/DetectorSite.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/SkyCoordinates.h>
#include <lal/Date.h>
#include <lal/TimeDelay.h>
#include <lal/DetResponse.h>
#include <lal/CoherentInspiral.h>

#define rint(x) (floor((x)+0.5))
double modf( double value, double *integerPart );
int compare( const void* a, const void* b );

NRCSID (COHERENTINSPIRALFILTERC, "$Id$");


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
  
  INITSTATUS( status, "LALCoherentInspiralFilterInputInit",
	      COHERENTINSPIRALFILTERC );
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

  INITSTATUS( status, "LALCoherentInspiralFilterInputFinalize",
	      COHERENTINSPIRALFILTERC );
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

  INITSTATUS( status, "LALCoherentInspiralFilterParamsInit",
	      COHERENTINSPIRALFILTERC );
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
        outputPtr->cohSNRVec = (REAL4TimeSeries *) 
          LALCalloc( 1, sizeof(REAL4TimeSeries) );
        LALCreateVector (status->statusPtr, &(outputPtr->cohSNRVec->data), 
		         params->numPoints);
        BEGINFAIL( status ) {
	  if( outputPtr->detIDVec )
	    {
	      TRY( LALU2DestroyVector (status->statusPtr, &(outputPtr->detIDVec)),
		   status);
	      LALFree( outputPtr->detIDVec );
	      outputPtr->detIDVec = NULL;
	    }
          LALFree( outputPtr->cohSNRVec );
          outputPtr->cohSNRVec = NULL;
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
	    TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec->data)),
		 status); 
	    LALFree( outputPtr->cohSNRVec );
	    outputPtr->cohSNRVec = NULL;
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
        outputPtr->nullStatVec = (REAL4TimeSeries *) 
          LALCalloc( 1, sizeof(REAL4TimeSeries) );
        LALCreateVector (status->statusPtr, &(outputPtr->nullStatVec->data), 
		         params->numPoints);
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
	      TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec->data)),
		   status); 
	      LALFree( outputPtr->cohSNRVec );
	      outputPtr->cohSNRVec = NULL;
	    } 
	  if( outputPtr->detIDVec )
	    {
	      TRY( LALU2DestroyVector (status->statusPtr, &(outputPtr->detIDVec)),
		   status);
	      LALFree( outputPtr->detIDVec );
	      outputPtr->detIDVec = NULL;
	    }
          LALFree( outputPtr->nullStatVec );
          outputPtr->nullStatVec = NULL;
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
	  TRY( LALDestroyVector( status->statusPtr, &(outputPtr->nullStatVec->data)),
	       status); 
	  LALFree( outputPtr->nullStatVec );
	  outputPtr->nullStatVec = NULL;
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
          TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec->data)),
	       status); 
          LALFree( outputPtr->cohSNRVec );
          outputPtr->cohSNRVec = NULL;
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
	  TRY( LALDestroyVector( status->statusPtr, &(outputPtr->nullStatVec->data)),
	       status); 
	  LALFree( outputPtr->nullStatVec );
	  outputPtr->nullStatVec = NULL;
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
          TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec->data)),
	       status);   
          LALFree( outputPtr->cohSNRVec );
          outputPtr->cohSNRVec = NULL;
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

  /* CHECK: Fixed the maximum number of participating ifos to be LAL_NUM_IFO */
  LALDCreateVector(status->statusPtr, &(outputPtr->sigmasqVec), LAL_NUM_IFO );
  if ( !(outputPtr->sigmasqVec ) ) {
    ABORT( status, COHERENTINSPIRALH_EALOC, COHERENTINSPIRALH_MSGEALOC);
    BEGINFAIL( status ) {
      if( params->nullStatOut )
	{
	  TRY( LALDestroyVector( status->statusPtr, &(outputPtr->nullStatVec->data)),
	       status); 
	  LALFree( outputPtr->nullStatVec );
	  outputPtr->nullStatVec = NULL;
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
          TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec->data)),
	       status);   
          LALFree( outputPtr->cohSNRVec );
          outputPtr->cohSNRVec = NULL;
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

  INITSTATUS( status, "LALCoherentInspiralFilterParamsInit",
	      COHERENTINSPIRALFILTERC );
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

  /*
   *
   * destroy coherent SNR vector, if it exists
   *
   */
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
   * destroy H1-H2 null statistic vector, if it exists
   *
   */
  if ( outputPtr->nullStatVec ) {
    LALDestroyVector( status->statusPtr, &(outputPtr->nullStatVec->data) );
    CHECKSTATUSPTR( status );
    LALFree( outputPtr->nullStatVec );
    outputPtr->nullStatVec = NULL;
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

  INITSTATUS( status, "LALCoherentInspiralEstimatePsiEpsilon", 
	      COHERENTINSPIRALFILTERC );
  ATTATCHSTATUSPTR( status );

  /* Must have 3 sites to estimate psi and epsilon */
  ASSERT( sigmasq[0] && sigmasq[1] && sigmasq[2], status,
	  COHERENTINSPIRALH_ENUMZ, COHERENTINSPIRALH_MSGENUMZ );
  ASSERT( cData[0].re && cData[1].re && cData[2].re, status,
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

  cPlus.re = 0.0;
  cPlus.im = 0.0;
  cMinus.re = 0.0;
  cMinus.im = 0.0;
  cRatio.re = 0.0;
  cRatio.im = 0.0;

  k = 0;
  for( i = 0; i < 6; i++ )
    {
      if( caseID[i] )
	{
	  cPlus.re += sigmasq[k] * ( cData[k].re * dVectorPlusRe[i] 
                - cData[k].im * dVectorPlusIm[i] );
	  cPlus.im += sigmasq[k] * ( cData[k].re * dVectorPlusIm[i] 
		+ cData[k].im * dVectorPlusRe[i] );
	  cMinus.re += sigmasq[k] * ( cData[k].re * dVectorMinusRe[i] 
                - cData[k].im * dVectorMinusIm[i] );
	  cMinus.im += sigmasq[k] * ( cData[k].re * dVectorMinusIm[i] 
		+ cData[k].im * dVectorMinusRe[i] );
	  k++;

	} 
    }

  cRatio.re = ( cPlus.re * cMinus.re + cPlus.im * cMinus.im ) / ( cPlus.re * 
          cPlus.re + cPlus.im * cPlus.im );
  cRatio.im = ( cPlus.re * cMinus.im - cPlus.im * cMinus.re ) / ( cPlus.re *
	  cPlus.re + cPlus.im * cPlus.im );

  /* Now the estimates can be computed */

  *inclination = acos(  ( 1 - sqrt( sqrt(cRatio.re*cRatio.re + cRatio.im*cRatio.im) ))/
	  ( 1 + sqrt( sqrt(cRatio.re*cRatio.re + cRatio.im*cRatio.im) )) );
  *polarization = 0.25 * atan( cRatio.im / cRatio.re );

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
	  cDotqRe += cData[k].re * qVectorRe[i] + cData[k].im * qVectorIm[i];
	  cDotqIm += cData[k].re * qVectorIm[i] - cData[k].im * qVectorRe[i];
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
    REAL4                                  templateNorm,
    REAL8                                  deltaT,
    INT4                                   segmentLength,  /* time pts */
    REAL4                                  coherentSNR,
    REAL4                                 *distance
    )
{
  INT4                 i = 0;
  REAL4                sigmaCoherent = 0.0;

  INITSTATUS( status, "LALCoherentInspiralEstimateDistance", 
	      COHERENTINSPIRALFILTERC );
  ATTATCHSTATUSPTR( status );


  /* CHECK: Assume that sigma[0] (sigma[1]) are of H1 (H2) and calculate 
     the effective distance for the H1-H2 pair */
  *distance = sqrt( 0.5 * (sigmasq[0] + sigmasq[1]) * deltaT / segmentLength) / coherentSNR;

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
 
void
LALCoherentInspiralFilterSegment (
    LALStatus                             *status,
    MultiInspiralTable                    **eventList,
    CoherentInspiralFilterInput           *input,
    CoherentInspiralFilterParams          *params
    )
{
  UINT4                               numDetectors = 0;
  UINT4                               numSegments = 0;
  UINT4                               numBeamPoints = 0;
  UINT4                               cohSNROut = 0;
  UINT4                               nullStatOut = 0;
  INT4                                caseID[6] = {0,0,0,0,0,0};
  INT4                                i,q,w,m,j,k,l;
  INT4                                found = 0;
  INT4                                detId = 0;
  INT4                                indexarray[4] = {0,0,0,0};
  INT4                                numPoints = 0;
  INT4                                deltaEventIndex = 0;
  INT4                                eventStartIdx = 0;
  INT4                                slidePoints[3] = {0,0,0};
  INT4                                segmentLength = 0;
  INT4                                sortedSlidePoints3D[3]= {0,0,0};
  INT4                                sortedSlidePoints4D[4]= {0,0,0,0};
  INT4                                detIdSlidTimePt = 0;
  int                                 locIdx,decIdx,raIdx,decMax,raMax;
  REAL4                               buffer = 0.0; /*account for timing errors*/
  REAL4                               timingError = 0.0; /*0.00025;*/  /* allowed timing error of 2 ms */
  REAL4                               s[3]={0,0,0};/* up to 4 distances;in 3D space*/
  REAL4                               nHatVect[3] = {0,0,0};
  REAL4                               distance[4] = {0,0,0,0};
  REAL4                               chirpTime = 0.0;
  REAL4                               cohSNRThresh = 0.0;
  REAL4                               cohSNRThreshSq = 0.0;
  REAL4                               inclination = 0.0;
  REAL4                               polarization = 0.0;
  REAL4                               distanceEstimate = 0.0;
  REAL4                               coaPhase = 0.0;
  REAL4                               cohSNRLocal = 0.0; 
  REAL4                               cohSNR        = 0.0;
  REAL4                               cohSNRSq = 0.0;
  REAL4                               cohSNRLocalSq = 0.0;
  REAL4                               cohSNRLocalRe = 0.0;
  REAL4                               cohSNRLocalIm = 0.0;
  REAL4                               cohSNRLocalSq1 = 0.0;
  REAL4                               cohSNRLocalSq2 = 0.0;
  REAL4                               nullStatRe    = 0.0 ;
  REAL4                               nullStatIm    = 0.0;
  REAL4                               nullNorm      = 0.0;
  REAL4                               cohSnrRe      = 0.0;
  REAL4                               cohSnrIm      = 0.0;
  REAL8                              *sigmasq = NULL;
  REAL8                               deltaT = 0.0;
  REAL8                               tempTime = 0.0;
  REAL8                               fracpart = 0.0;
  REAL8                               intpart = 0.0;
  double                              decStep = 0.0;
  double                              raStep = 0.0;
  double                              theta = 0.0;
  double                              phi = 0.0;
  double                              timeDelay[4]= {0.0,0.0,0.0,0.0};
  double                              sortedDelays3D[3]= {0.0,0.0,0.0};
  double                              sortedDelays4D[4]= {0.0,0.0,0.0,0.0};
  double                              fplus[4] = {0.0,0.0,0.0,0.0};
  double                              fcross[4] = {0.0,0.0,0.0,0.0};
  COMPLEX8                            cDataTemp[4];
  COMPLEX8                            quadTemp[6];
  LALDetector                         detectors[4];
  COMPLEX8TimeSeries                 *cData[4] = {NULL,NULL,NULL,NULL};
  MultiInspiralTable                 *thisEvent = NULL; 
  CoherentInspiralBeamVector         *beamVec = NULL;
  CHAR                                idtag[6][3] = {"G1","H1","H2","L1","T1","V1"};
  CHAR                                caseStr[FILENAME_MAX];

	    
  INT4           timePt[4] = {0,0,0,0};
  INT4           timePtTemp[4] = {0,0,0,0};
  REAL4          AA=0.0;
  REAL4          BB=0.0;
  REAL4          CC=0.0;
  REAL4          discrimSqrt=0.0;
  REAL4          VVPlus[4]={0.0,0.0,0.0,0.0};
  REAL4          VVMinus[4]={0.0,0.0,0.0,0.0};
  REAL4          CRePlus=0.0;
  REAL4          CImPlus=0.0;
  REAL4          CReMinus=0.0;
  REAL4          CImMinus=0.0;
  REAL4          AAn[4]={0.0,0.0,0.0,0.0};
  REAL4          BBn[4]={0.0,0.0,0.0,0.0};
  REAL4          CCn[4]={0.0,0.0,0.0,0.0};
  REAL4          discrimSqrtn[4]={0.0,0.0,0.0,0.0};  
  REAL4          MM1=0.0;
  REAL4          MM2=0.0;
  REAL4          O11=0.0;
  REAL4          O12=0.0;
  REAL4          O21=0.0;
  REAL4          O22=0.0;

  INITSTATUS( status, "LALCoherentInspiralFilterSegment", 
	      COHERENTINSPIRALFILTERC );
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
    ASSERT( params->cohSNRVec->data->data, status, 
	    COHERENTINSPIRALH_ENULL, COHERENTINSPIRALH_MSGENULL );
    ASSERT( params->cohSNRVec->data, status, 
	    COHERENTINSPIRALH_ENULL, COHERENTINSPIRALH_MSGENULL );
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
  deltaT = params->deltaT;
  segmentLength = params->segmentLength;
  raStep = params->raStep; 
  decStep = params->decStep; 
  chirpTime = params->chirpTime;
  deltaEventIndex = (UINT4) rint( (chirpTime / deltaT) + 1.0 );
  buffer = rint( (timingError/deltaT) + 1.0 );

  /* if the full coherent snr vector is required, set it to zero */
  if ( cohSNROut ) {
    memset( params->cohSNRVec->data->data, 0, numPoints * sizeof( REAL4 ));
  }

  /*CHECK: hardwired to 6 detectors for now */
  for (l=0 ;l< 6 ; l++)
    {
      caseID[l] = params->detIDVec->data[l];
    }
  
  i = 0;  
  for (l=0 ;l<6 ;l++)
    {
      if( caseID[l] )
	{
	  indexarray[i] = l;
	  i++;
	}
    }

  /* CHECK: Note that this may be replaced with 
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
      LALSnprintf( caseStr, FILENAME_MAX * sizeof(CHAR), "%s%s",idtag[i],idtag[j]);
      break;

    case 3:
      i=indexarray[0];
      j=indexarray[1];
      k=indexarray[2];
      LALSnprintf( caseStr, FILENAME_MAX * sizeof(CHAR), "%s%s%s",idtag[i],idtag[j],idtag[k]);
      break;

    case 4:
      i=indexarray[0];
      j=indexarray[1];
      k=indexarray[2];
      l=indexarray[3];
      LALSnprintf( caseStr, FILENAME_MAX * sizeof(CHAR), "%s%s%s%s",idtag[i],idtag[j],idtag[k],idtag[l]);
      break;
    }
    
  /*** get detector-site locations */
  sigmasq = params->sigmasqVec->data;

  /* read in CData */  
  for ( l=0 ; l < (INT4) params->numDetectors ; l++) {
    cData[l] = input->multiCData->cData[l];
  }
  
  /* Now construct the appropriate coherent SNR */ 
  switch (params->numDetectors)  {
  case 2:
    /* Network: H1 and H2*/
    if(caseID[1] && caseID[2]) {
      m = 0;
      for (k=0;k<(INT4)numPoints;k++) {
	REAL4 cohSNR = 0.0;

	for (m=k-buffer; m<k+buffer; m++) {
	  if(m >=0 && m < (INT4) numPoints) {
	    cohSNRLocalRe = sqrt(sigmasq[1])*cData[0]->data->data[k].re
	      + sqrt(sigmasq[2])*cData[1]->data->data[m].re;
	    cohSNRLocalIm = sqrt(sigmasq[1])*cData[0]->data->data[k].im
	      + sqrt(sigmasq[2])*cData[1]->data->data[m].im;
	    
	    cohSNRLocal = sqrt( (cohSNRLocalRe*cohSNRLocalRe + cohSNRLocalIm*cohSNRLocalIm) / 
				(sigmasq[1] + sigmasq[2] ) );
	    
	    if(cohSNRLocal > cohSNR) {
	      cohSNR = cohSNRLocal;
	      quadTemp[0].re=cData[0]->data->data[k].re;
	      quadTemp[0].im=cData[0]->data->data[k].im;
	      quadTemp[1].re=cData[1]->data->data[m].re;
	      quadTemp[1].im=cData[1]->data->data[m].im;
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
	      
	      found=0; /* ifo counter */	
	      if(caseID[0]) {
		thisEvent->g1quad.re=quadTemp[found].re;	    
		thisEvent->g1quad.im=quadTemp[found].im;
		found=found+1;
	      }
	      else {
		thisEvent->g1quad.re=0;	    
		thisEvent->g1quad.im=0;
	      }
	      if(caseID[1]) {
		thisEvent->h1quad.re=quadTemp[found].re;	    
		thisEvent->h1quad.im=quadTemp[found].im;
		found=found+1;
	      }
	      else {
		thisEvent->h1quad.re=0;	    
		thisEvent->h1quad.im=0;
	      }
	      if(caseID[2]) {
                thisEvent->h2quad.re=quadTemp[found].re;	    
                thisEvent->h2quad.im=quadTemp[found].im;
                found=found+1;
	      }
	      else {
                thisEvent->h2quad.re=0;	    
                thisEvent->h2quad.im=0;
	      }
	      if(caseID[3]) {
                thisEvent->l1quad.re=quadTemp[found].re;	    
                thisEvent->l1quad.im=quadTemp[found].im;
                found=found+1;		   
	      }
	      else {
                thisEvent->l1quad.re=0;	    
                thisEvent->l1quad.im=0;
	      }
	      if(caseID[4]) {
                thisEvent->t1quad.re=quadTemp[found].re;	    
                thisEvent->t1quad.im=quadTemp[found].im;
                found=found+1;
	      }
	      else {
                thisEvent->t1quad.re=0;	    
                thisEvent->t1quad.im=0;
	      }
	      if(caseID[5]) {
                thisEvent->v1quad.re=quadTemp[found].re;	    
                thisEvent->v1quad.im=quadTemp[found].im;
                found=found+1;
	      }
	      else {
                thisEvent->v1quad.re=0;	    
                thisEvent->v1quad.im=0;
	      }		
	      
	      
	      thisEvent->snr = cohSNR;
	      strcpy(thisEvent->ifos,caseStr);
	      thisEvent->mass1 = input->tmplt->mass1;
	      thisEvent->mass2 = input->tmplt->mass2;
	      thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
	      thisEvent->eta = input->tmplt->eta;		       
	      /*Calculate distance/effective distance */
	      LALCoherentInspiralEstimateDistance( status->statusPtr, sigmasq, params->templateNorm, deltaT, segmentLength, cohSNR, &distanceEstimate );
	      thisEvent->eff_distance = distanceEstimate;
	      
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
	      
	      found=0;			
	      if(caseID[0])
                {
		  thisEvent->g1quad.re=quadTemp[found].re;	    
		  thisEvent->g1quad.im=quadTemp[found].im;
		  found=found+1;
                }
	      else
                {
		  thisEvent->g1quad.re=0;	    
		  thisEvent->g1quad.im=0;
                }
	      if(caseID[1])
                {
		  thisEvent->h1quad.re=quadTemp[found].re;	    
		  thisEvent->h1quad.im=quadTemp[found].im;
		  found=found+1;
                }
	      else
                {
		  thisEvent->h1quad.re=0;	    
		  thisEvent->h1quad.im=0;
                }
	      if(caseID[2])
                {
		  thisEvent->h2quad.re=quadTemp[found].re;	    
		  thisEvent->h2quad.im=quadTemp[found].im;
		  found=found+1;
                }
	      else
                { 
		  thisEvent->h2quad.re=0;	    
		  thisEvent->h2quad.im=0;
                }
	      if(caseID[3])
                {
		  thisEvent->l1quad.re=quadTemp[found].re;	    
		  thisEvent->l1quad.im=quadTemp[found].im;
		  found=found+1;		   
                }
	      else
                {
		  thisEvent->l1quad.re=0;	    
		  thisEvent->l1quad.im=0;
                }
	      if(caseID[4])
                {
		  thisEvent->t1quad.re=quadTemp[found].re;	    
		  thisEvent->t1quad.im=quadTemp[found].im;
		  found=found+1;
                }
	      else
                {
		  thisEvent->t1quad.re=0;	    
		  thisEvent->t1quad.im=0;
                }
	      if(caseID[5])
                {
		  thisEvent->v1quad.re=quadTemp[found].re;	    
		  thisEvent->v1quad.im=quadTemp[found].im;
		  found=found+1;
                }
	      else
                {
		  thisEvent->v1quad.re=0;	    
		  thisEvent->v1quad.im=0;
                }		
	      thisEvent->snr = cohSNR;
	      strcpy(thisEvent->ifos, caseStr);
	      thisEvent->mass1 = input->tmplt->mass1;
	      thisEvent->mass2 = input->tmplt->mass2;
	      thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
	      thisEvent->eta = input->tmplt->eta;
	      /*Calculate distance/effective distance */
	      LALCoherentInspiralEstimateDistance( status->statusPtr, sigmasq, params->templateNorm, deltaT, segmentLength, cohSNR, &distanceEstimate );
	      thisEvent->eff_distance = distanceEstimate;
	      
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
	      
	      found=0;			
	      if(caseID[0])
                {
		  thisEvent->g1quad.re=quadTemp[found].re;	    
		  thisEvent->g1quad.im=quadTemp[found].im;
		  found=found+1;
                }
	      else
                {
		  thisEvent->g1quad.re=0;	    
		  thisEvent->g1quad.im=0;
                }
	      if(caseID[1])
                {
		  thisEvent->h1quad.re=quadTemp[found].re;	    
		  thisEvent->h1quad.im=quadTemp[found].im;
		  found=found+1;
                }
	      else
                {
		  thisEvent->h1quad.re=0;	    
		  thisEvent->h1quad.im=0;
                }
	      if(caseID[2])
                {
		  thisEvent->h2quad.re=quadTemp[found].re;	    
		  thisEvent->h2quad.im=quadTemp[found].im;
		  found=found+1;
                }
	      else
                { 
		  thisEvent->h2quad.re=0;	    
		  thisEvent->h2quad.im=0;
                }
	      if(caseID[3])
                {
		  thisEvent->l1quad.re=quadTemp[found].re;	    
		  thisEvent->l1quad.im=quadTemp[found].im;
		  found=found+1;		   
                }
	      else
                {
		  thisEvent->l1quad.re=0;	    
		  thisEvent->l1quad.im=0;
                }
	      if(caseID[4])
                {
		  thisEvent->t1quad.re=quadTemp[found].re;	    
		  thisEvent->t1quad.im=quadTemp[found].im;
		  found=found+1;
                }
	      else
                {
		  thisEvent->t1quad.re=0;	    
		  thisEvent->t1quad.im=0;
                }
	      if(caseID[5])
                {
		  thisEvent->v1quad.re=quadTemp[found].re;	    
		  thisEvent->v1quad.im=quadTemp[found].im;
		  found=found+1;
                }
	      else
                {
		  thisEvent->v1quad.re=0;	    
		  thisEvent->v1quad.im=0;
                }		
	      
	      thisEvent->snr = cohSNR;
	      strcpy(thisEvent->ifos,caseStr);
	      thisEvent->mass1 = input->tmplt->mass1;
	      thisEvent->mass2 = input->tmplt->mass2;
	      thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
	      thisEvent->eta = input->tmplt->eta;
	      /*Calculate distance/effective distance */
	      LALCoherentInspiralEstimateDistance( status->statusPtr, sigmasq, params->templateNorm, deltaT, segmentLength, cohSNR, &distanceEstimate );
	      thisEvent->eff_distance = distanceEstimate;
	      
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
	for (i=0;i<3;i++) {
	  s[i] = (REAL4) ( detectors[1].location[i] - detectors[0].location[i]);
	}
	distance[1] = sqrt( cartesianInnerProduct(s,s) ); 
	timeDelay[1] = distance[1]/LAL_C_SI;
	slidePoints[1] = rint( (fabs(timeDelay[1])/deltaT) + 1.0 );
	
	k = 0;
	q = 0;
	w = 0;

	for(k=0;k<(INT4)numPoints;k++)
	  {
	    REAL4 cohSNR = 0.0;
	    
	    for (q = k-slidePoints[1]-buffer; q < k+slidePoints[1]+buffer; q++)
	      {
		if(q >= 0 && q < (INT4) numPoints)
		  {
		    cohSNRLocal = sqrt(cData[0]->data->data[k].re * 
			cData[0]->data->data[k].re + cData[1]->data->data[q].re * 
			cData[1]->data->data[q].re + cData[0]->data->data[k].im * 
			cData[0]->data->data[k].im + cData[1]->data->data[q].im * 
			cData[1]->data->data[q].im);
		    
		    
		    if(cohSNRLocal > cohSNR)
		      {
			cohSNR = cohSNRLocal; 
			
			quadTemp[0].re=cData[0]->data->data[k].re;
			quadTemp[0].im=cData[0]->data->data[k].im;
		       	quadTemp[1].re=cData[1]->data->data[q].re;
			quadTemp[1].im=cData[1]->data->data[q].im;
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

	        found=0;			
		if(caseID[0])
		  {
		    thisEvent->g1quad.re=quadTemp[found].re;	    
		    thisEvent->g1quad.im=quadTemp[found].im;
		    found=found+1;
		  }
		else
		  {
		    thisEvent->g1quad.re=0;	    
		    thisEvent->g1quad.im=0;
		  }
		if(caseID[1])
		  {
		    thisEvent->h1quad.re=quadTemp[found].re;	    
		    thisEvent->h1quad.im=quadTemp[found].im;
        	    found=found+1;
		  }
		else
		  {
		    thisEvent->h1quad.re=0;	    
		    thisEvent->h1quad.im=0;
		  }
		if(caseID[2])
		  {
		    thisEvent->h2quad.re=quadTemp[found].re;	    
		    thisEvent->h2quad.im=quadTemp[found].im;
		    found=found+1;
		  }
		else
		  { 
		    thisEvent->h2quad.re=0;	    
		    thisEvent->h2quad.im=0;
		  }
		if(caseID[3])
		  {
		    thisEvent->l1quad.re=quadTemp[found].re;	    
		    thisEvent->l1quad.im=quadTemp[found].im;
		    found=found+1;		   
		  }
		else
		  {
		    thisEvent->l1quad.re=0;	    
		    thisEvent->l1quad.im=0;
		  }
		if(caseID[4])
		  {
		    thisEvent->t1quad.re=quadTemp[found].re;	    
		    thisEvent->t1quad.im=quadTemp[found].im;
		    found=found+1;
		  }
		else
		  {
		    thisEvent->t1quad.re=0;	    
		    thisEvent->t1quad.im=0;
		  }
		if(caseID[5])
		  {
		    thisEvent->v1quad.re=quadTemp[found].re;	    
		    thisEvent->v1quad.im=quadTemp[found].im;
		    found=found+1;
		  }
		else
		  {
		    thisEvent->v1quad.re=0;	    
		    thisEvent->v1quad.im=0;
		  }		
		
		thisEvent->snr = cohSNR;
		strcpy(thisEvent->ifos,caseStr);
		thisEvent->mass1 = input->tmplt->mass1;
		thisEvent->mass2 = input->tmplt->mass2;
		thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
		thisEvent->eta = input->tmplt->eta;
		/* CHECK: LALCoherentInspiralEstimateDistance( status->statusPtr, sigmasq, params->templateNorm, deltaT, segmentLength, cohSNR, &distanceEstimate );
		   calculates a valid effective distance for H1-H2. Since not both 
		   are present here, set the effective distance to zero
		*/
		thisEvent->eff_distance = 0.0;
		
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

		found=0;			
		if(caseID[0])
		  {
		    thisEvent->g1quad.re=quadTemp[found].re;	    
		    thisEvent->g1quad.im=quadTemp[found].im;
		    found=found+1;
		  }
		else
		  {
		    thisEvent->g1quad.re=0;	    
		    thisEvent->g1quad.im=0;
		  }
		if(caseID[1])
		  {
		    thisEvent->h1quad.re=quadTemp[found].re;	    
		    thisEvent->h1quad.im=quadTemp[found].im;
        	    found=found+1;
		  }
		else
		  {
		    thisEvent->h1quad.re=0;	    
		    thisEvent->h1quad.im=0;
		  }
		if(caseID[2])
		  {
		    thisEvent->h2quad.re=quadTemp[found].re;	    
		    thisEvent->h2quad.im=quadTemp[found].im;
		    found=found+1;
		  }
		else
		  { 
		    thisEvent->h2quad.re=0;	    
		    thisEvent->h2quad.im=0;
		  }
		if(caseID[3])
		  {
		    thisEvent->l1quad.re=quadTemp[found].re;	    
		    thisEvent->l1quad.im=quadTemp[found].im;
		    found=found+1;		   
		  }
		else
		  {
		    thisEvent->l1quad.re=0;	    
		    thisEvent->l1quad.im=0;
		  }
		if(caseID[4])
		  {
		    thisEvent->t1quad.re=quadTemp[found].re;	    
		    thisEvent->t1quad.im=quadTemp[found].im;
		    found=found+1;
		  }
		else
		  {
		    thisEvent->t1quad.re=0;	    
		    thisEvent->t1quad.im=0;
		  }
		if(caseID[5])
		  {
		    thisEvent->v1quad.re=quadTemp[found].re;	    
		    thisEvent->v1quad.im=quadTemp[found].im;
		    found=found+1;
		  }
		else
		  {
		    thisEvent->v1quad.re=0;	    
		    thisEvent->v1quad.im=0;
		  }		
		
		thisEvent->snr = cohSNR;
		strcpy(thisEvent->ifos, caseStr);
		thisEvent->mass1 = input->tmplt->mass1;
		thisEvent->mass2 = input->tmplt->mass2;
		thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
		thisEvent->eta = input->tmplt->eta;
		LALCoherentInspiralEstimateDistance( status->statusPtr, sigmasq, params->templateNorm, deltaT, segmentLength, cohSNR, &distanceEstimate );
		thisEvent->eff_distance = distanceEstimate;
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

		found=0;			
		if(caseID[0])
		  {
		    thisEvent->g1quad.re=quadTemp[found].re;	    
		    thisEvent->g1quad.im=quadTemp[found].im;
		    found=found+1;
		  }
		else
		  {
		    thisEvent->g1quad.re=0;	    
		    thisEvent->g1quad.im=0;
		  }
		if(caseID[1])
		  {
		    thisEvent->h1quad.re=quadTemp[found].re;	    
		    thisEvent->h1quad.im=quadTemp[found].im;
        	    found=found+1;
		  }
		else
		  {
		    thisEvent->h1quad.re=0;	    
		    thisEvent->h1quad.im=0;
		  }
		if(caseID[2])
		  {
		    thisEvent->h2quad.re=quadTemp[found].re;	    
		    thisEvent->h2quad.im=quadTemp[found].im;
		    found=found+1;
		  }
		else
		  { 
		    thisEvent->h2quad.re=0;	    
		    thisEvent->h2quad.im=0;
		  }
		if(caseID[3])
		  {
		    thisEvent->l1quad.re=quadTemp[found].re;	    
		    thisEvent->l1quad.im=quadTemp[found].im;
		    found=found+1;		   
		  }
		else
		  {
		    thisEvent->l1quad.re=0;	    
		    thisEvent->l1quad.im=0;
		  }
		if(caseID[4])
		  {
		    thisEvent->t1quad.re=quadTemp[found].re;	    
		    thisEvent->t1quad.im=quadTemp[found].im;
		    found=found+1;
		  }
		else
		  {
		    thisEvent->t1quad.re=0;	    
		    thisEvent->t1quad.im=0;
		  }
		if(caseID[5])
		  {
		    thisEvent->v1quad.re=quadTemp[found].re;	    
		    thisEvent->v1quad.im=quadTemp[found].im;
		    found=found+1;
		  }
		else
		  {
		    thisEvent->v1quad.re=0;	    
		    thisEvent->v1quad.im=0;
		  }		
		
		thisEvent->snr = cohSNR;
		strcpy(thisEvent->ifos,caseStr);
		thisEvent->mass1 = input->tmplt->mass1;
		thisEvent->mass2 = input->tmplt->mass2;
		thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
		thisEvent->eta = input->tmplt->eta;
		LALCoherentInspiralEstimateDistance( status->statusPtr, sigmasq, params->templateNorm, deltaT, segmentLength, cohSNR, &distanceEstimate );
		thisEvent->eff_distance = distanceEstimate;
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
	    REAL4 cohSNR = 0.0;
	    for(m=k-buffer;m<k+buffer;m++)
	      {
		if(m >=0 && m < (INT4) numPoints)
		  {
		    for (q = m-slidePoints[1]-buffer;q < m+slidePoints[1]+buffer;q++)
		      {
			if(q >= 0 && q < (INT4) numPoints)
			  {
			    /*CHECK: This will NOT work if G1 is present! because it assumes that the "0" det is H1 and "1" det is H2! Rectify in next rev. */
			    cohSNRLocalRe = sqrt(sigmasq[1])*cData[0]->data->data[k].re
			      + sqrt(sigmasq[2])*cData[1]->data->data[m].re;
			    cohSNRLocalIm = sqrt(sigmasq[1])*cData[0]->data->data[k].im
			      + sqrt(sigmasq[2])*cData[1]->data->data[m].im;
			    
			    cohSNRLocal = (cohSNRLocalRe*cohSNRLocalRe + cohSNRLocalIm*cohSNRLocalIm) / 
			      (sigmasq[1] + sigmasq[2]) ;
			    
			    cohSNRLocal += cData[2]->data->data[q].re*cData[2]->data->data[q].re + 
			      cData[2]->data->data[q].im*cData[2]->data->data[q].im;

			    cohSNRLocal = sqrt( cohSNRLocal );
			    
			  if(cohSNRLocal > cohSNR)
			      {
				cohSNR = cohSNRLocal;
	    		  
				quadTemp[0].re=cData[0]->data->data[k].re;
				quadTemp[0].im=cData[0]->data->data[k].im;
				quadTemp[1].re=cData[1]->data->data[m].re;
				quadTemp[1].im=cData[1]->data->data[m].im; 
				quadTemp[2].re=cData[2]->data->data[q].re;
				quadTemp[2].im=cData[2]->data->data[q].im; 
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
		    
		    found=0;			
		    if(caseID[0])
		      {
			thisEvent->g1quad.re=quadTemp[found].re;	    
			thisEvent->g1quad.im=quadTemp[found].im;
			found=found+1;
		      }
		    else
		      {
			thisEvent->g1quad.re=0;	    
			thisEvent->g1quad.im=0;
		      }
		    if(caseID[1])
		      {
			thisEvent->h1quad.re=quadTemp[found].re;	    
			thisEvent->h1quad.im=quadTemp[found].im;
			found=found+1;
		      }
		    else
		      {
			thisEvent->h1quad.re=0;	    
			thisEvent->h1quad.im=0;
		      }
		    if(caseID[2])
		      {
			thisEvent->h2quad.re=quadTemp[found].re;	    
			thisEvent->h2quad.im=quadTemp[found].im;
			found=found+1;
		      }
		    else
		      { 
			thisEvent->h2quad.re=0;	    
			thisEvent->h2quad.im=0;
		      }
		    if(caseID[3])
		      {
			thisEvent->l1quad.re=quadTemp[found].re;	    
			thisEvent->l1quad.im=quadTemp[found].im;
			found=found+1;		   
		      }
		    else
		      {
			thisEvent->l1quad.re=0;	    
			thisEvent->l1quad.im=0;
		      }
		    if(caseID[4])
		      {
			thisEvent->t1quad.re=quadTemp[found].re;	    
			thisEvent->t1quad.im=quadTemp[found].im;
			found=found+1;
		      }
		    else
		      {
			thisEvent->t1quad.re=0;	    
			thisEvent->t1quad.im=0;
		      }
		    if(caseID[5])
		      {
			thisEvent->v1quad.re=quadTemp[found].re;	    
			thisEvent->v1quad.im=quadTemp[found].im;
			found=found+1;
		      }
		    else
		      {
			thisEvent->v1quad.re=0;	    
			thisEvent->v1quad.im=0;
		      }		
		    
		    thisEvent->snr = cohSNR;
		    strcpy(thisEvent->ifos,caseStr);
		    thisEvent->mass1 = input->tmplt->mass1;
		    thisEvent->mass2 = input->tmplt->mass2;
		    thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
		    thisEvent->eta = input->tmplt->eta;
		    LALCoherentInspiralEstimateDistance( status->statusPtr, sigmasq, params->templateNorm, deltaT, segmentLength, cohSNR, &distanceEstimate );
		    thisEvent->eff_distance = distanceEstimate;
		    
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
		    
		    found=0;			
		    if(caseID[0])
		      {
			thisEvent->g1quad.re=quadTemp[found].re;	    
			thisEvent->g1quad.im=quadTemp[found].im;
			found=found+1;
		      }
		    else
		      {
			thisEvent->g1quad.re=0;	    
			thisEvent->g1quad.im=0;
		      }
		    if(caseID[1])
		      {
			thisEvent->h1quad.re=quadTemp[found].re;	    
			thisEvent->h1quad.im=quadTemp[found].im;
			found=found+1;
		      }
		    else
		      {
			thisEvent->h1quad.re=0;	    
			thisEvent->h1quad.im=0;
		      }
		    if(caseID[2])
		      {
			thisEvent->h2quad.re=quadTemp[found].re;	    
			thisEvent->h2quad.im=quadTemp[found].im;
			found=found+1;
		      }
		    else
		      { 
			thisEvent->h2quad.re=0;	    
			thisEvent->h2quad.im=0;
		      }
		    if(caseID[3])
		      {
			thisEvent->l1quad.re=quadTemp[found].re;	    
			thisEvent->l1quad.im=quadTemp[found].im;
			found=found+1;		   
		      }
		    else
		      {
			thisEvent->l1quad.re=0;	    
			thisEvent->l1quad.im=0;
		      }
		    if(caseID[4])
		      {
			thisEvent->t1quad.re=quadTemp[found].re;	    
			thisEvent->t1quad.im=quadTemp[found].im;
			found=found+1;
		      }
		    else
		      {
			thisEvent->t1quad.re=0;	    
			thisEvent->t1quad.im=0;
		      }
		    if(caseID[5])
		      {
			thisEvent->v1quad.re=quadTemp[found].re;	    
			thisEvent->v1quad.im=quadTemp[found].im;
			found=found+1;
		      }
		    else
		      {
			thisEvent->v1quad.re=0;	    
			thisEvent->v1quad.im=0;
		      }		
		    
		    thisEvent->snr = cohSNR;
		    strcpy(thisEvent->ifos, caseStr);
		    thisEvent->mass1 = input->tmplt->mass1;
		    thisEvent->mass2 = input->tmplt->mass2;
		    thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
		    thisEvent->eta = input->tmplt->eta;
		    LALCoherentInspiralEstimateDistance( status->statusPtr, sigmasq, params->templateNorm, deltaT, segmentLength, cohSNR, &distanceEstimate );
		    thisEvent->eff_distance = distanceEstimate;
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
		    
		    found=0;			
		    if(caseID[0])
		      {
			thisEvent->g1quad.re=quadTemp[found].re;	    
			thisEvent->g1quad.im=quadTemp[found].im;
			found=found+1;
		      }
		    else
		      {
			thisEvent->g1quad.re=0;	    
			thisEvent->g1quad.im=0;
		      }
		    if(caseID[1])
		      {
			thisEvent->h1quad.re=quadTemp[found].re;	    
			thisEvent->h1quad.im=quadTemp[found].im;
			found=found+1;
		      }
		    else
		      {
			thisEvent->h1quad.re=0;	    
			thisEvent->h1quad.im=0;
		      }
		    if(caseID[2])
		      {
			thisEvent->h2quad.re=quadTemp[found].re;	    
			thisEvent->h2quad.im=quadTemp[found].im;
			found=found+1;
		      }
		    else
		      { 
			thisEvent->h2quad.re=0;	    
			thisEvent->h2quad.im=0;
		      }
		    if(caseID[3])
		      {
			thisEvent->l1quad.re=quadTemp[found].re;	    
			thisEvent->l1quad.im=quadTemp[found].im;
			found=found+1;		   
		      }
		    else
		      {
			thisEvent->l1quad.re=0;	    
			thisEvent->l1quad.im=0;
		      }
		    if(caseID[4])
		      {
			thisEvent->t1quad.re=quadTemp[found].re;	    
			thisEvent->t1quad.im=quadTemp[found].im;
			found=found+1;
		      }
		    else
		      {
			thisEvent->t1quad.re=0;	    
			thisEvent->t1quad.im=0;
		      }
		    if(caseID[5])
		      {
			thisEvent->v1quad.re=quadTemp[found].re;	    
			thisEvent->v1quad.im=quadTemp[found].im;
			found=found+1;
		      }
		    else
		      {
			thisEvent->v1quad.re=0;	    
			thisEvent->v1quad.im=0;
		      }		
		    
		    thisEvent->snr = cohSNR;
		    strcpy(thisEvent->ifos,caseStr);
		    thisEvent->mass1 = input->tmplt->mass1;
		    thisEvent->mass2 = input->tmplt->mass2;
		    thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
		    thisEvent->eta = input->tmplt->eta;
		    LALCoherentInspiralEstimateDistance( status->statusPtr, sigmasq, params->templateNorm, deltaT, segmentLength, cohSNR, &distanceEstimate );
		    thisEvent->eff_distance = distanceEstimate;
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
      }
    else
      { /* Network: 3 detectors excluding either H1, H2, or both (H1 && H2)*/
	/*Now the last 3 cases will involve the looping over the coefficients*/
	LIGOTimeGPS 	triggerGPSEndTime;/* Needed to calculate time-delays */
	LALMSTUnitsAndAcc pUnitsAndAcc;
	REAL8           gmstInRadians=0.0;
	double          psiInRadians = 0.0;
	double          detRefLocation[3];
	double          detNextLocation[3];

	triggerGPSEndTime.gpsSeconds = cData[0]->epoch.gpsSeconds;
	triggerGPSEndTime.gpsNanoSeconds = cData[0]->epoch.gpsNanoSeconds;
	/* Convert GPS time of trigger to GMST time in radians 
	   for computing F+, Fx */
	pUnitsAndAcc.units = MST_RAD;
	pUnitsAndAcc.accuracy = LALLEAPSEC_LOOSE;
	LALGPStoGMST1(status->statusPtr,&gmstInRadians,&triggerGPSEndTime,&pUnitsAndAcc);
	
	/* Following needed because XLALArrivalTimeDiff() uses doubles */
	for ( locIdx=0 ; locIdx<3 ; locIdx++ ) {
	  detRefLocation[locIdx] = (double) detectors[0].location[locIdx];
	}
	
	/* Convert to radians */
	raMax = floor(360/raStep);
	decMax = floor(180/decStep);
	
	raStep *= LAL_PI_180;
	decStep *= LAL_PI_180;
	
	/* Loop over time-points in the reference detector */	
	for( timePt[0]=0 ; timePt[0]<(INT4)numPoints ; timePt[0]++) {
	  /* Reset cohSNR to zero so that it can only be ratched upward by
	     cohSNRLocal computed below for every point in sky-position grid*/
	  cohSNR = 0.0;

	  /* Loop over points in the sky-position grid */
	  for (decIdx=1; decIdx<decMax; decIdx++) {
	    for ( raIdx=0; raIdx<raMax; raIdx++) {
	      phi = raIdx * raStep;
	      theta = decIdx * decStep - 90*LAL_PI_180;
	      
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

	      /* Loop over time-points in reference detector, after
		 accounting for the rounded-off sortedSlidePoints */
	      if( ( timePt[0] < (0-sortedSlidePoints3D[0]) )
		  ||  ( timePt[0]> (numPoints-sortedSlidePoints3D[2]) ) ) {
		cohSNR = 0.0;
		if( cohSNROut ) params->cohSNRVec->data->data[timePt[0]] = cohSNR;
	      }
	      else {
		/* Compute antenna-patterns and coherent SNR */
		/* Loop over detectors computing the theta-, phi-dependent
		   pieces of F+, Fx */
		detId = 0;
		for( j=0; j<LAL_NUM_IFO; j++ ) {
		  /* Compute antenna-patterns if caseID[j] != 0 */		  		
		  if ( !(params->detIDVec->data[j] == 0 )) { 
		    XLALComputeDetAMResponse(&fplus[detId], &fcross[detId],
				    detectors[detId].response, phi, theta,
				    psiInRadians, (double) gmstInRadians);
		    
		    /* Compute antenna-pattern factors */
		    AAn[detId] = ( fplus[detId]) * ( fplus[detId]);
		    BBn[detId] = ( fplus[detId]) * ( fcross[detId]);
		    CCn[detId] = ( fcross[detId]) * ( fcross[detId]);
		    
		    discrimSqrtn[detId] = sqrt(AAn[detId]*AAn[detId] 
					       + 4*BBn[detId]*BBn[detId]
					       - 2*AAn[detId]*CCn[detId] 
					       + CCn[detId]*CCn[detId]);
		    
		    O11 = ( AAn[detId] - CCn[detId] - discrimSqrtn[detId]);
		    O11 /= (BBn[detId]
			    *sqrt( 4 + ( AAn[detId] - CCn[detId] - discrimSqrtn[detId])
				   *( AAn[detId] - CCn[detId] - discrimSqrtn[detId])
				   / ( BBn[detId]*BBn[detId] ) ) );
		    O12 = 1 / sqrt( 1 + (-AAn[detId]+CCn[detId]+discrimSqrtn[detId])
				    *(-AAn[detId]+CCn[detId]+discrimSqrtn[detId]) 
				    / ( 4*BBn[detId]*BBn[detId] ) );
		    O21 = ( AAn[detId] - CCn[detId] + discrimSqrtn[detId]);
		    O21 /= ( BBn[detId]
			     *sqrt( 4 + (-AAn[detId]+CCn[detId]-discrimSqrtn[detId])
				    *(-AAn[detId]+CCn[detId]-discrimSqrtn[detId]) 
				    / ( BBn[detId]*BBn[detId] ) ) );
		    O22 = 1 / sqrt( 1 + (-AAn[detId]+CCn[detId]-discrimSqrtn[detId])
				    *(-AAn[detId]+CCn[detId]-discrimSqrtn[detId]) 
				    / ( 4*BBn[detId]*BBn[detId] ) );
		    
		    VVPlus[detId] = O11 * ( fplus[detId]) 
		      + O12 * ( fcross[detId]);		
		    VVMinus[detId] = O21 * ( fplus[detId]) 
		      + O22 * ( fcross[detId]);		
		    
		    /* CHECK: If sigmasq should be in the denominator!
		       Compute the elements of the helicity-plane projection matrix */
		    AAn[detId] *= sigmasq[j];
		    BBn[detId] *= sigmasq[j];
		    CCn[detId] *= sigmasq[j];
		    VVPlus[detId] *= sqrt(sigmasq[j]);
		    VVMinus[detId] *= sqrt(sigmasq[j]);
		    
		    detId++;
		  }
		}
		/* Construct network terms and factors required for
		   computing the coherent statistics */
		AA = AAn[0] + AAn[1] +AAn[2];
		BB = BBn[0] + BBn[1] +BBn[2];
		CC = CCn[0] + CCn[1] +CCn[2];
		
		discrimSqrt = sqrt(AA*AA + 4*BB*BB
				   - 2*AA*CC + CC*CC);
		MM1 = 2 * ( AA*AA*AA 
			    - 2*BB*BB*discrimSqrt
			    - AA*AA*( 2*CC + discrimSqrt )
			    + AA*( 4*BB*BB 
				   + CC * (CC + discrimSqrt ) ) );
		MM1 /= ( 4*BB*BB + (-AA+CC+discrimSqrt)*(-AA+CC+discrimSqrt) );
		MM2 = 2 * ( AA*AA*AA 
			    + 2*BB*BB*discrimSqrt
			    + AA*AA*(-2*CC + discrimSqrt )
			    + AA*( 4*BB*BB 
				   + CC * (CC - discrimSqrt ) ) );
		MM2 /= ( 4*BB*BB + (AA-CC+discrimSqrt)*(AA-CC+discrimSqrt) );
		
		/* Regularize */
		if ( (MM1<1.0e-4 ) ) 
		  MM1=1.0e-4;
		if ( (MM2<1.0e-4 ) ) 
		  MM2=1.0e-4;

		/*Initialize cohSNR components and time stamps */
		CRePlus = 0.0;
		CImPlus = 0.0;
		CReMinus = 0.0;
		CImMinus = 0.0;
		timePtTemp[1] = 0;
		timePtTemp[2] = 0;

		/* Compute components of the coherent SNR */		
		for ( detId=0 ; detId<params->numDetectors ; detId++ ) {
		  detIdSlidTimePt = timePt[0]+slidePoints[detId];
		  
		  CRePlus += VVPlus[detId] * 
		    cData[detId]->data->data[detIdSlidTimePt].re;
		  CImPlus += VVPlus[detId] * 
		    cData[detId]->data->data[detIdSlidTimePt].im;
		  CReMinus += VVMinus[detId] * 
		    cData[detId]->data->data[detIdSlidTimePt].re;
		  CImMinus += VVMinus[detId] *
		    cData[detId]->data->data[detIdSlidTimePt].im;
		}
		
		/* Compute coherent SNR */
		cohSNRLocal = sqrt( ( CRePlus*CRePlus/MM1 + 
				      CImPlus*CImPlus/MM1 + 
				      CReMinus*CReMinus/MM2 + 
				      CImMinus*CImMinus/MM2 ) );
		
		if(cohSNRLocal > cohSNR) {
		  cohSNR = cohSNRLocal;
		  for ( detId=0 ; detId<params->numDetectors ; detId++ ) {
		    detIdSlidTimePt = timePt[0]+slidePoints[detId];		    
		    quadTemp[detId].re=cData[detId]->data->data[detIdSlidTimePt].re;
		    quadTemp[detId].im=cData[detId]->data->data[detIdSlidTimePt].im;
		    timePtTemp[detId] = detIdSlidTimePt;
		  }
		}
		
		if( cohSNROut ) params->cohSNRVec->data->data[timePt[0]] = cohSNR;
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
		      
		      found=0;			
		      if(caseID[0])
			{
			  thisEvent->g1quad.re=quadTemp[found].re;	    
			  thisEvent->g1quad.im=quadTemp[found].im;
			  found=found+1;
			}
		      else
			{
			  thisEvent->g1quad.re=0;	    
			  thisEvent->g1quad.im=0;
			}
		      if(caseID[1])
			{
			  thisEvent->h1quad.re=quadTemp[found].re;	    
			  thisEvent->h1quad.im=quadTemp[found].im;
			  found=found+1;
			}
		      else
			{
			  thisEvent->h1quad.re=0;	    
			  thisEvent->h1quad.im=0;
			}
		      if(caseID[2])
			{
			  thisEvent->h2quad.re=quadTemp[found].re;	    
			  thisEvent->h2quad.im=quadTemp[found].im;
			  found=found+1;
			}
		      else
			{ 
			  thisEvent->h2quad.re=0;	    
			  thisEvent->h2quad.im=0;
			}
		      if(caseID[3])
			{
			  thisEvent->l1quad.re=quadTemp[found].re;	    
			  thisEvent->l1quad.im=quadTemp[found].im;
			  found=found+1;		   
			}
		      else
			{
			  thisEvent->l1quad.re=0;	    
			  thisEvent->l1quad.im=0;
			}
		      if(caseID[4])
			{
			  thisEvent->t1quad.re=quadTemp[found].re;	    
			  thisEvent->t1quad.im=quadTemp[found].im;
			  found=found+1;
			}
		      else
			{
			  thisEvent->t1quad.re=0;	    
			  thisEvent->t1quad.im=0;
			}
		      if(caseID[5])
			{
			  thisEvent->v1quad.re=quadTemp[found].re;	    
			  thisEvent->v1quad.im=quadTemp[found].im;
			  found=found+1;
			}
		      else
			{
			  thisEvent->v1quad.re=0;	    
			  thisEvent->v1quad.im=0;
			}		
		      
		      thisEvent->snr = cohSNR;
		      strcpy(thisEvent->ifos,caseStr);
		      thisEvent->mass1 = input->tmplt->mass1;
		      thisEvent->mass2 = input->tmplt->mass2;
		      thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
		      thisEvent->eta = input->tmplt->eta;
		      
		      inclination = 0.0;
		      polarization = 0.0;
		      /* CHECK: Move this to post-processing to save time
		      for ( detId=0 ; detId<params->numDetectors ; detId++ ) {
			cDataTemp[detId].re=cData[detId]->data->data[timePtTemp[detId]].re;
			cDataTemp[detId].im=cData[detId]->data->data[timePtTemp[detId]].im;
		      }
		      LALCoherentInspiralEstimatePsiEpsilonCoaPhase( status->statusPtr, caseID, sigmasq, (REAL4) theta, (REAL4) phi, cDataTemp, &inclination, &polarization, &coaPhase ); 
		      */
		      thisEvent->inclination = inclination;
		      thisEvent->polarization = polarization;
		      thisEvent->coa_phase = coaPhase;
		      /* CHECK: LALCoherentInspiralEstimateDistance( status->statusPtr, sigmasq, params->templateNorm, deltaT, segmentLength, cohSNR, &distanceEstimate );
			 calculates a valid effective distance for H1-H2. Since not both 
			 are present here, set the effective distance to zero */
		      thisEvent->eff_distance = 0.0;
		      
		      thisEvent->ligo_axis_ra = (REAL4) phi;
		      thisEvent->ligo_axis_dec = (REAL4) theta;
		      
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
		      
		      found=0;			
		      if(caseID[0])
			 {
			   thisEvent->g1quad.re=quadTemp[found].re;	    
			   thisEvent->g1quad.im=quadTemp[found].im;
			   found=found+1;
			 }
		       else
			 {
			   thisEvent->g1quad.re=0;	    
			   thisEvent->g1quad.im=0;
			 }
		       if(caseID[1])
			 {
			   thisEvent->h1quad.re=quadTemp[found].re;	    
			   thisEvent->h1quad.im=quadTemp[found].im;
			   found=found+1;
			 }
		       else
			 {
			   thisEvent->h1quad.re=0;	    
			   thisEvent->h1quad.im=0;
			 }
		       if(caseID[2])
			 {
			   thisEvent->h2quad.re=quadTemp[found].re;	    
			   thisEvent->h2quad.im=quadTemp[found].im;
			   found=found+1;
			 }
		       else
			 { 
			   thisEvent->h2quad.re=0;	    
			   thisEvent->h2quad.im=0;
			 }
		       if(caseID[3])
			 {
			   thisEvent->l1quad.re=quadTemp[found].re;	    
			   thisEvent->l1quad.im=quadTemp[found].im;
			   found=found+1;		   
			 }
		       else
			 {
			   thisEvent->l1quad.re=0;	    
			   thisEvent->l1quad.im=0;
			 }
		       if(caseID[4])
			 {
			   thisEvent->t1quad.re=quadTemp[found].re;	    
			   thisEvent->t1quad.im=quadTemp[found].im;
			   found=found+1;
			 }
		       else
			 {
			   thisEvent->t1quad.re=0;	    
			   thisEvent->t1quad.im=0;
			 }
		       if(caseID[5])
			 {
			   thisEvent->v1quad.re=quadTemp[found].re;	    
			   thisEvent->v1quad.im=quadTemp[found].im;
			   found=found+1;
			 }
		       else
			 {
			   thisEvent->v1quad.re=0;	    
			   thisEvent->v1quad.im=0;
			 }		
		       
		       thisEvent->snr = cohSNR;
		       strcpy(thisEvent->ifos, caseStr);
		       thisEvent->mass1 = input->tmplt->mass1;
		       thisEvent->mass2 = input->tmplt->mass2;
		       thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
		       thisEvent->eta = input->tmplt->eta;

		       inclination = 0.0;
		       polarization = 0.0;
		       /* CHECK: Move this to post-processing to save time
		       for ( detId=0 ; detId<params->numDetectors ; detId++ ) {
			 cDataTemp[detId].re=cData[detId]->data->data[timePtTemp[detId]].re;
			 cDataTemp[detId].im=cData[detId]->data->data[timePtTemp[detId]].im;
		       }
		       LALCoherentInspiralEstimatePsiEpsilonCoaPhase( status->statusPtr, caseID, sigmasq, (REAL4) theta, (REAL4) phi, cDataTemp, &inclination, &polarization, &coaPhase ); 
		      */
		       thisEvent->inclination = inclination;
		       thisEvent->polarization = polarization;
		       thisEvent->coa_phase = coaPhase;
		       /* CHECK: LALCoherentInspiralEstimateDistance( status->statusPtr, sigmasq, params->templateNorm, deltaT, segmentLength, cohSNR, &distanceEstimate );
			  calculates a valid effective distance for H1-H2. Since not both 
			  are present here, set the effective distance to zero */
		       thisEvent->eff_distance = 0.0;
		       thisEvent->ligo_axis_ra = (REAL4) phi;
		       thisEvent->ligo_axis_dec = (REAL4) theta;

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
		       
		       found=0;			
		       if(caseID[0])
			 {
			   thisEvent->g1quad.re=quadTemp[found].re;	    
			   thisEvent->g1quad.im=quadTemp[found].im;
			   found=found+1;
			 }
		       else
			 {
			   thisEvent->g1quad.re=0;	    
			   thisEvent->g1quad.im=0;
			 }
		       if(caseID[1])
			 {
			   thisEvent->h1quad.re=quadTemp[found].re;	    
			   thisEvent->h1quad.im=quadTemp[found].im;
			   found=found+1;
			 }
		       else
			 {
			   thisEvent->h1quad.re=0;	    
			   thisEvent->h1quad.im=0;
			 }
		       if(caseID[2])
			 {
			   thisEvent->h2quad.re=quadTemp[found].re;	    
			   thisEvent->h2quad.im=quadTemp[found].im;
			   found=found+1;
			 }
		       else
			 { 
			   thisEvent->h2quad.re=0;	    
			   thisEvent->h2quad.im=0;
			 }
		       if(caseID[3])
			 {
			   thisEvent->l1quad.re=quadTemp[found].re;	    
			   thisEvent->l1quad.im=quadTemp[found].im;
			   found=found+1;		   
			 }
		       else
			 {
			   thisEvent->l1quad.re=0;	    
			   thisEvent->l1quad.im=0;
			 }
		       if(caseID[4])
			 {
			   thisEvent->t1quad.re=quadTemp[found].re;	    
			   thisEvent->t1quad.im=quadTemp[found].im;
			   found=found+1;
			 }
		       else
			 {
			   thisEvent->t1quad.re=0;	    
			   thisEvent->t1quad.im=0;
			 }
		       if(caseID[5])
			 {
			   thisEvent->v1quad.re=quadTemp[found].re;	    
			   thisEvent->v1quad.im=quadTemp[found].im;
			   found=found+1;
			 }
		       else
			 {
			   thisEvent->v1quad.re=0;	    
			   thisEvent->v1quad.im=0;
			 }		
		       
		       thisEvent->snr = cohSNR;
		       strcpy(thisEvent->ifos,caseStr);
		       thisEvent->mass1 = input->tmplt->mass1;
		       thisEvent->mass2 = input->tmplt->mass2;
		       thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
		       thisEvent->eta = input->tmplt->eta;

		       inclination = 0.0;
		       polarization = 0.0;
		       /* CHECK: Move this to post-processing to save time
		       for ( detId=0 ; detId<params->numDetectors ; detId++ ) {
			 cDataTemp[detId].re=cData[detId]->data->data[timePtTemp[detId]].re;
			 cDataTemp[detId].im=cData[detId]->data->data[timePtTemp[detId]].im;
		       }
;		       LALCoherentInspiralEstimatePsiEpsilonCoaPhase( status->statusPtr, caseID, sigmasq, (REAL4) theta, (REAL4) phi, cDataTemp, &inclination, &polarization, &coaPhase ); 
		      */
		       thisEvent->inclination = inclination;
		       thisEvent->polarization = polarization;
		       thisEvent->coa_phase = coaPhase;
		       /* CHECK: LALCoherentInspiralEstimateDistance( status->statusPtr, sigmasq, params->templateNorm, deltaT, segmentLength, cohSNR, &distanceEstimate );
			  calculates a valid effective distance for H1-H2. Since not both 
			  are present here, set the effective distance to zero */
		       thisEvent->eff_distance = 0.0;
		       thisEvent->ligo_axis_ra = (REAL4) phi;
		       thisEvent->ligo_axis_dec = (REAL4) theta;

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
	      } /* ends loop over right-ascension index "raIdx" */
	    } /* ends loop over declination index "decIdx" */	      }
	} /* ends the "for loop" over the reference detector's timePt[0] */
      } /* ends loop over case 3 */
    break;
  case 4: /* Network: 4 detectors */
    {
	LIGOTimeGPS 	triggerGPSEndTime;/* Needed to calculate time-delays */
	LALMSTUnitsAndAcc pUnitsAndAcc;
	REAL8           gmstInRadians=0.0;
	double          psiInRadians = 0.0;
	double          detRefLocation[3];
	double          detNextLocation[3];
	
	triggerGPSEndTime.gpsSeconds = cData[0]->epoch.gpsSeconds;
	triggerGPSEndTime.gpsNanoSeconds = cData[0]->epoch.gpsNanoSeconds;
	/* Convert GPS time of trigger to GMST time in radians 
	   for computing F+, Fx */
	pUnitsAndAcc.units = MST_RAD;
	pUnitsAndAcc.accuracy = LALLEAPSEC_LOOSE;
	LALGPStoGMST1(status->statusPtr,&gmstInRadians,&triggerGPSEndTime,&pUnitsAndAcc);
	
	/* Following needed because XLALArrivalTimeDiff() uses doubles */
	for ( locIdx=0 ; locIdx<3 ; locIdx++ ) {
	  detRefLocation[locIdx] = (double) detectors[0].location[locIdx];
	}
	
	/* Convert to radians */
	raMax = floor(360/raStep);
	decMax = floor(180/decStep);
	
	raStep *= LAL_PI_180;
	decStep *= LAL_PI_180;
	
	/* Loop over time-points in the reference detector */	
	for( timePt[0]=0 ; timePt[0]<(INT4)numPoints ; timePt[0]++) {
	  /* Reset cohSNR to zero so that it can only be ratched upward by
	     cohSNRLocal computed below for every point in sky-position grid*/
	  cohSNR = 0.0;

	  /* Loop over points in the sky-position grid */
	  for (decIdx=1; decIdx<decMax; decIdx++) {
	    for ( raIdx=0; raIdx<raMax; raIdx++) {
	      phi = raIdx * raStep;
	      theta = decIdx * decStep - 90*LAL_PI_180;
	      
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
		  slidePoints[detId] = rint( timeDelay[detId]/deltaT );

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

	      /* Loop over time-points in reference detector, after
		 accounting for the rounded-off sortedSlidePoints */
	      if( ( timePt[0] < (0-sortedSlidePoints4D[0]) )
		  ||  ( timePt[0]> (numPoints-sortedSlidePoints4D[2]) ) ) {
		cohSNR = 0.0;
		if( cohSNROut ) params->cohSNRVec->data->data[timePt[0]] = cohSNR;
	      }
	      else {
		/* Compute antenna-patterns and coherent SNR */
		/* Loop over detectors computing the theta-, phi-dependent
		   pieces of F+, Fx */
		detId = 0;
		for( j=0; j<LAL_NUM_IFO; j++ ) {
		  /* Compute antenna-patterns if caseID[j] != 0 */		  		
		  if ( !(params->detIDVec->data[j] == 0 )) { 
		    XLALComputeDetAMResponse(&fplus[detId], &fcross[detId],
				    detectors[detId].response, phi, theta,
				    psiInRadians, (double) gmstInRadians);
		    
		    /* Compute antenna-pattern factors */
		    AAn[detId] = ( fplus[detId]) * ( fplus[detId]);
		    BBn[detId] = ( fplus[detId]) * ( fcross[detId]);
		    CCn[detId] = ( fcross[detId]) * ( fcross[detId]);
		    
		    discrimSqrtn[detId] = sqrt(AAn[detId]*AAn[detId] 
					       + 4*BBn[detId]*BBn[detId]
					       - 2*AAn[detId]*CCn[detId] 
					       + CCn[detId]*CCn[detId]);
		    
		    O11 = ( AAn[detId] - CCn[detId] - discrimSqrtn[detId]);
		    O11 /= (BBn[detId]
			    *sqrt( 4 + ( AAn[detId] - CCn[detId] - discrimSqrtn[detId])
				   *( AAn[detId] - CCn[detId] - discrimSqrtn[detId])
				   / ( BBn[detId]*BBn[detId] ) ) );
		    O12 = 1 / sqrt( 1 + (-AAn[detId]+CCn[detId]+discrimSqrtn[detId])
				    *(-AAn[detId]+CCn[detId]+discrimSqrtn[detId]) 
				    / ( 4*BBn[detId]*BBn[detId] ) );
		    O21 = ( AAn[detId] - CCn[detId] + discrimSqrtn[detId]);
		    O21 /= ( BBn[detId]
			     *sqrt( 4 + (-AAn[detId]+CCn[detId]-discrimSqrtn[detId])
				    *(-AAn[detId]+CCn[detId]-discrimSqrtn[detId]) 
				    / ( BBn[detId]*BBn[detId] ) ) );
		    O22 = 1 / sqrt( 1 + (-AAn[detId]+CCn[detId]-discrimSqrtn[detId])
				    *(-AAn[detId]+CCn[detId]-discrimSqrtn[detId]) 
				    / ( 4*BBn[detId]*BBn[detId] ) );
		    
		    VVPlus[detId] = O11 * ( fplus[detId]) 
		      + O12 * ( fcross[detId]);		
		    VVMinus[detId] = O21 * ( fplus[detId]) 
		      + O22 * ( fcross[detId]);		
		    
		    /* CHECK: If sigmasq should be in the denominator!
		       Compute the elements of the helicity-plane projection matrix */
		    AAn[detId] *= sigmasq[j];
		    BBn[detId] *= sigmasq[j];
		    CCn[detId] *= sigmasq[j];
		    VVPlus[detId] *= sqrt(sigmasq[j]);
		    VVMinus[detId] *= sqrt(sigmasq[j]);
		    
		    detId++;
		  }
		}
		/* Construct network terms and factors required for
		   computing the coherent statistics */
		AA = AAn[0] + AAn[1] +AAn[2] + AAn[3];
		BB = BBn[0] + BBn[1] +BBn[2] + BBn[3];
		CC = CCn[0] + CCn[1] +CCn[2] + CCn[3];
		
		discrimSqrt = sqrt(AA*AA + 4*BB*BB
				   - 2*AA*CC + CC*CC);
		MM1 = 2 * ( AA*AA*AA 
			    - 2*BB*BB*discrimSqrt
			    - AA*AA*( 2*CC + discrimSqrt )
			    + AA*( 4*BB*BB 
				   + CC * (CC + discrimSqrt ) ) );
		MM1 /= ( 4*BB*BB + (-AA+CC+discrimSqrt)*(-AA+CC+discrimSqrt) );
		MM2 = 2 * ( AA*AA*AA 
			    + 2*BB*BB*discrimSqrt
			    + AA*AA*(-2*CC + discrimSqrt )
			    + AA*( 4*BB*BB 
				   + CC * (CC - discrimSqrt ) ) );
		MM2 /= ( 4*BB*BB + (AA-CC+discrimSqrt)*(AA-CC+discrimSqrt) );
		
		/* Regularize */
		if ( (MM1<1.0e-4 ) ) 
		  MM1=1.0e-4;
		if ( (MM2<1.0e-4 ) ) 
		  MM2=1.0e-4;

		/*Initialize cohSNR components and time stamps */
		CRePlus = 0.0;
		CImPlus = 0.0;
		CReMinus = 0.0;
		CImMinus = 0.0;
		timePtTemp[1] = 0;
		timePtTemp[2] = 0;
		timePtTemp[3] = 0;

		/* Compute components of the coherent SNR */		
		for ( detId=0 ; detId<params->numDetectors ; detId++ ) {
		  detIdSlidTimePt = timePt[0]+slidePoints[detId];
		  
		  CRePlus += VVPlus[detId] * 
		    cData[detId]->data->data[detIdSlidTimePt].re;
		  CImPlus += VVPlus[detId] * 
		    cData[detId]->data->data[detIdSlidTimePt].im;
		  CReMinus += VVMinus[detId] * 
		    cData[detId]->data->data[detIdSlidTimePt].re;
		  CImMinus += VVMinus[detId] *
		    cData[detId]->data->data[detIdSlidTimePt].im;
		}
		
		/* Compute coherent SNR */
		cohSNRLocal = sqrt( ( CRePlus*CRePlus/MM1 + 
				      CImPlus*CImPlus/MM1 + 
				      CReMinus*CReMinus/MM2 + 
				      CImMinus*CImMinus/MM2 ) );
		
		if(cohSNRLocal > cohSNR) {
		  cohSNR = cohSNRLocal;
		  for ( detId=0 ; detId<params->numDetectors ; detId++ ) {
		    detIdSlidTimePt = timePt[0]+slidePoints[detId];		    
		    quadTemp[detId].re=cData[detId]->data->data[detIdSlidTimePt].re;
		    quadTemp[detId].im=cData[detId]->data->data[detIdSlidTimePt].im;
		    timePtTemp[detId] = detIdSlidTimePt;
		  }
		}
		
		if( cohSNROut ) params->cohSNRVec->data->data[timePt[0]] = cohSNR;
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
			
			found=0;			
			if(caseID[0])
			  {
			    thisEvent->g1quad.re=quadTemp[found].re;	    
			    thisEvent->g1quad.im=quadTemp[found].im;
			    found=found+1;
			  }
			else
			  {
			    thisEvent->g1quad.re=0;	    
			    thisEvent->g1quad.im=0;
			  }
			if(caseID[1])
			  {
			    thisEvent->h1quad.re=quadTemp[found].re;	    
			    thisEvent->h1quad.im=quadTemp[found].im;
			    found=found+1;
			  }
			else
			  {
			    thisEvent->h1quad.re=0;	    
			    thisEvent->h1quad.im=0;
			  }
			if(caseID[2])
			  {
			    thisEvent->h2quad.re=quadTemp[found].re;	    
			    thisEvent->h2quad.im=quadTemp[found].im;
			    found=found+1;
			  }
			else
			  { 
			    thisEvent->h2quad.re=0;	    
			    thisEvent->h2quad.im=0;
			  }
			if(caseID[3])
			  {
			    thisEvent->l1quad.re=quadTemp[found].re;	    
			    thisEvent->l1quad.im=quadTemp[found].im;
			    found=found+1;		   
			  }
			else
			  {
			    thisEvent->l1quad.re=0;	    
			    thisEvent->l1quad.im=0;
			  }
			if(caseID[4])
			  {
			    thisEvent->t1quad.re=quadTemp[found].re;	    
			    thisEvent->t1quad.im=quadTemp[found].im;
			    found=found+1;
			  }
			else
			  {
			    thisEvent->t1quad.re=0;	    
			    thisEvent->t1quad.im=0;
			  }
			if(caseID[5])
			  {
			    thisEvent->v1quad.re=quadTemp[found].re;	    
			    thisEvent->v1quad.im=quadTemp[found].im;
			    found=found+1;
			  }
			else
			  {
			    thisEvent->v1quad.re=0;	    
			    thisEvent->v1quad.im=0;
			  }		
			
		        thisEvent->snr = cohSNR;
		        strcpy(thisEvent->ifos,caseStr);
		        thisEvent->mass1 = input->tmplt->mass1;
		        thisEvent->mass2 = input->tmplt->mass2;
			thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
			thisEvent->eta = input->tmplt->eta;

		        inclination = 0.0;
		        polarization = 0.0;
			/* CHECK: Move this to post-processing to save time
			for ( detId=0 ; detId<params->numDetectors ; detId++ ) {
			  cDataTemp[detId].re=cData[detId]->data->data[timePtTemp[detId]].re;
			  cDataTemp[detId].im=cData[detId]->data->data[timePtTemp[detId]].im;
			}
		        LALCoherentInspiralEstimatePsiEpsilonCoaPhase( status->statusPtr, caseID, sigmasq, (REAL4) theta, (REAL4) phi, cDataTemp, &inclination, &polarization, &coaPhase ); 
			*/
		        thisEvent->inclination = inclination;
		        thisEvent->polarization = polarization;
			thisEvent->coa_phase = coaPhase;
			/* CHECK: LALCoherentInspiralEstimateDistance( status->statusPtr, sigmasq, params->templateNorm, deltaT, segmentLength, cohSNR, &distanceEstimate );
			   calculates a valid effective distance for H1-H2. Since not both 
			   are present here, set the effective distance to zero */
		        thisEvent->eff_distance = 0.0;
			
   	                thisEvent->ligo_axis_ra = (REAL4) phi;
		        thisEvent->ligo_axis_dec = (REAL4) theta;
			
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
			
			found=0;			
			if(caseID[0])
			  {
			    thisEvent->g1quad.re=quadTemp[found].re;	    
			    thisEvent->g1quad.im=quadTemp[found].im;
			    found=found+1;
			  }
			else
			  {
			    thisEvent->g1quad.re=0;	    
			    thisEvent->g1quad.im=0;
			  }
			if(caseID[1])
			  {
			    thisEvent->h1quad.re=quadTemp[found].re;	    
			    thisEvent->h1quad.im=quadTemp[found].im;
			    found=found+1;
			  }
			else
			  {
			    thisEvent->h1quad.re=0;	    
			    thisEvent->h1quad.im=0;
			  }
			if(caseID[2])
			  {
			    thisEvent->h2quad.re=quadTemp[found].re;	    
			    thisEvent->h2quad.im=quadTemp[found].im;
			    found=found+1;
			  }
			else
			  { 
			    thisEvent->h2quad.re=0;	    
			    thisEvent->h2quad.im=0;
			  }
			if(caseID[3])
			  {
			    thisEvent->l1quad.re=quadTemp[found].re;	    
			    thisEvent->l1quad.im=quadTemp[found].im;
			    found=found+1;		   
			  }
			else
			  {
			    thisEvent->l1quad.re=0;	    
			    thisEvent->l1quad.im=0;
			  }
			if(caseID[4])
			  {
			    thisEvent->t1quad.re=quadTemp[found].re;	    
			    thisEvent->t1quad.im=quadTemp[found].im;
			    found=found+1;
			  }
			else
			  {
			    thisEvent->t1quad.re=0;	    
			    thisEvent->t1quad.im=0;
			  }
			if(caseID[5])
			  {
			    thisEvent->v1quad.re=quadTemp[found].re;	    
			    thisEvent->v1quad.im=quadTemp[found].im;
			    found=found+1;
			  }
			else
			  {
			    thisEvent->v1quad.re=0;	    
			    thisEvent->v1quad.im=0;
			  }		
			
		        thisEvent->snr = cohSNR;
		        strcpy(thisEvent->ifos, caseStr);
		        thisEvent->mass1 = input->tmplt->mass1;
		        thisEvent->mass2 = input->tmplt->mass2;
			thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
			thisEvent->eta = input->tmplt->eta;

		        inclination = 0.0;
		        polarization = 0.0;
			/* CHECK: Move this to post-processing to save time
			for ( detId=0 ; detId<params->numDetectors ; detId++ ) {
			  cDataTemp[detId].re=cData[detId]->data->data[timePtTemp[detId]].re;
			  cDataTemp[detId].im=cData[detId]->data->data[timePtTemp[detId]].im;
			}
		        LALCoherentInspiralEstimatePsiEpsilonCoaPhase( status->statusPtr, caseID, sigmasq, (REAL4) theta, (REAL4) phi, cDataTemp, &inclination, &polarization, &coaPhase ); 
			*/
		        thisEvent->inclination = inclination;
		        thisEvent->polarization = polarization;
			thisEvent->coa_phase = coaPhase;
			/* CHECK: LALCoherentInspiralEstimateDistance( status->statusPtr, sigmasq, params->templateNorm, deltaT, segmentLength, cohSNR, &distanceEstimate );
			   calculates a valid effective distance for H1-H2. Since not both 
			   are present here, set the effective distance to zero */
		        thisEvent->eff_distance = 0.0;
		        thisEvent->ligo_axis_ra = (REAL4) phi;
		        thisEvent->ligo_axis_dec = (REAL4) theta;

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
			
			found=0;			
			if(caseID[0])
			  {
			    thisEvent->g1quad.re=quadTemp[found].re;	    
			    thisEvent->g1quad.im=quadTemp[found].im;
			    found=found+1;
			  }
			else
			  {
			    thisEvent->g1quad.re=0;	    
			    thisEvent->g1quad.im=0;
			  }
			if(caseID[1])
			  {
			    thisEvent->h1quad.re=quadTemp[found].re;	    
			    thisEvent->h1quad.im=quadTemp[found].im;
			    found=found+1;
			  }
			else
			  {
			    thisEvent->h1quad.re=0;	    
			    thisEvent->h1quad.im=0;
			  }
			if(caseID[2])
			  {
			    thisEvent->h2quad.re=quadTemp[found].re;	    
			    thisEvent->h2quad.im=quadTemp[found].im;
			    found=found+1;
			  }
			else
			  { 
			    thisEvent->h2quad.re=0;	    
			    thisEvent->h2quad.im=0;
			  }
			if(caseID[3])
			  {
			    thisEvent->l1quad.re=quadTemp[found].re;	    
			    thisEvent->l1quad.im=quadTemp[found].im;
			    found=found+1;		   
			  }
			else
			  {
			    thisEvent->l1quad.re=0;	    
			    thisEvent->l1quad.im=0;
			  }
			if(caseID[4])
			  {
			    thisEvent->t1quad.re=quadTemp[found].re;	    
			    thisEvent->t1quad.im=quadTemp[found].im;
			    found=found+1;
			  }
			else
			  {
			    thisEvent->t1quad.re=0;	    
			    thisEvent->t1quad.im=0;
			  }
			if(caseID[5])
			  {
			    thisEvent->v1quad.re=quadTemp[found].re;	    
			    thisEvent->v1quad.im=quadTemp[found].im;
			    found=found+1;
			  }
			else
			  {
			    thisEvent->v1quad.re=0;	    
			    thisEvent->v1quad.im=0;
			  }		
			
		        thisEvent->snr = cohSNR;
		        strcpy(thisEvent->ifos,caseStr);
		        thisEvent->mass1 = input->tmplt->mass1;
		        thisEvent->mass2 = input->tmplt->mass2;
			thisEvent->mchirp = input->tmplt->totalMass * pow( input->tmplt->eta, 3.0/5.0 );
			thisEvent->eta = input->tmplt->eta;

		        inclination = 0.0;
		        polarization = 0.0;
			/* CHECK: Move this to post-processing to save time
			for ( detId=0 ; detId<params->numDetectors ; detId++ ) {
			  cDataTemp[detId].re=cData[detId]->data->data[timePtTemp[detId]].re;
			  cDataTemp[detId].im=cData[detId]->data->data[timePtTemp[detId]].im;
			}
		        LALCoherentInspiralEstimatePsiEpsilonCoaPhase( status->statusPtr, caseID, sigmasq, (REAL4) theta, (REAL4) phi, cDataTemp, &inclination, &polarization, &coaPhase ); 
			*/
		        thisEvent->inclination = inclination;
		        thisEvent->polarization = polarization;
			thisEvent->coa_phase = coaPhase;
		        LALCoherentInspiralEstimateDistance( status->statusPtr, sigmasq, params->templateNorm, deltaT, segmentLength, cohSNR, &distanceEstimate );
		        thisEvent->eff_distance = distanceEstimate;
		        thisEvent->ligo_axis_ra = (REAL4) phi;
		        thisEvent->ligo_axis_dec = (REAL4) theta;

		        /* Need to initialize the event start index to new value */
		        if( timePt[0] > (eventStartIdx + deltaEventIndex) )
		          {
			    eventStartIdx = timePt[0];
		          }
		        tempTime = 0.0;
		        fracpart= 0.0;
		        intpart = 0.0;
		      }
		} /* matches if (cohSNR > cohSNRThresh) */
	      }/* close if/else loop on time-points */
	    }/* end loop over right-ascension, raIdx */
	  }/* end loop over declination, decIdx */
	} /* end loop over timePts in the reference detector*/
    }/* case of 4 detectors */
  } /* closes  switch(params->numDetectors) */
  
  /* Compute null-statistic, just for H1-H2 as of now,
     and cohSNRH1H2, if not computed above already */
  if( thisEvent && params->nullStatOut && params->cohH1H2SNROut ) {
    
    /* Prepare norm for null statistic */
    nullNorm = ( 1.0 / sigmasq[1]  + 1.0 /  sigmasq[2] );

    /* Allocate memory for null statistic */
    memset( params->nullStatVec->data->data, 0, numPoints*sizeof(REAL4));

    /* Allocate memory for cohSNRH1H2Vec if that SNR has 
       not been computed above already*/
    memset( params->cohH1H2SNRVec->data->data, 0, numPoints*sizeof(REAL4));

    /*CHECK: Will not give intended result if first det is "G1", since it 
      assumes that cdata[0] is H1 and cdata[1] is H1; rectify this in next rev. */
    for (k=0;k<(INT4)numPoints;k++) {

      /*Compute cohH1H2 snr, with index m replaced by k */
      cohSnrRe = sqrt(sigmasq[1])*cData[0]->data->data[k].re
        + sqrt(sigmasq[2])*cData[1]->data->data[k].re;
      cohSnrIm = sqrt(sigmasq[1])*cData[0]->data->data[k].im
        + sqrt(sigmasq[2])*cData[1]->data->data[k].im;

      params->cohH1H2SNRVec->data->data[k]
        = sqrt( (cohSnrRe*cohSnrRe + cohSnrIm*cohSnrIm) /
                (sigmasq[1] + sigmasq[2] ) );

      /* Compute null-stream statistic; 
         in next rev. report re and im parts separately */
      nullStatRe = cData[0]->data->data[k].re / sqrt(sigmasq[1])
        - cData[1]->data->data[k].re / sqrt(sigmasq[2]);
      nullStatIm = cData[0]->data->data[k].im / sqrt(sigmasq[1])
        - cData[1]->data->data[k].im / sqrt(sigmasq[2]);
      params->nullStatVec->data->data[k] = 
        ( nullStatRe*nullStatRe + nullStatIm*nullStatIm ) / nullNorm ;
    }
    thisEvent->null_statistic = params->nullStatVec->data->data[(INT4)(numPoints/2)];
  }

  /* Compute null-statistic ONLY, just for H1-H2 as of now, and NOT cohH1H2SNR */
  if( thisEvent && params->nullStatOut && !(params->cohH1H2SNROut) ) {

    /* Prepare norm for null statistic */
    nullNorm = ( 1.0 / sigmasq[1]  + 1.0 /  sigmasq[2] );

    /* Allocate memory for null statistic */
    memset( params->nullStatVec->data->data, 0, numPoints*sizeof(REAL4));

    /*CHECK: Will not give intended result if first det is "G1", since it 
      assumes that cdata[0] is H1 and cdata[1] is H2; rectify this in next rev. */
    for (k=0;k<(INT4)numPoints;k++) {
      /* Compute null-stream statistic; 
         in next rev. report re and im parts separately */
      nullStatRe = cData[0]->data->data[k].re / sqrt(sigmasq[1])
        - cData[1]->data->data[k].re / sqrt(sigmasq[2]);
      nullStatIm = cData[0]->data->data[k].im / sqrt(sigmasq[1])
        - cData[1]->data->data[k].im / sqrt(sigmasq[2]);
      params->nullStatVec->data->data[k] = 
        ( nullStatRe*nullStatRe + nullStatIm*nullStatIm ) / nullNorm ;
    }
    thisEvent->null_statistic = params->nullStatVec->data->data[(INT4)(numPoints/2)];
  }

  /* Compute cohSNRH1H2 ONLY, if not computed above already, and NOT null-stat */
  if( params->cohH1H2SNROut && !nullStatOut ) {
    /* Allocate memory for cohSNRH1H2Vec if that SNR has 
       not been computed above already*/
    memset( params->cohH1H2SNRVec->data->data, 0, numPoints*sizeof(REAL4));

    /*CHECK: Will not give intended result if first det is "G1", since it 
      assumes that cdata[0] is H1 and cdata[1] is H1; rectify this in next rev. */
    for (k=0;k<(INT4)numPoints;k++) {

      /*Compute cohH1H2 snr, with index m replaced by k */
      cohSnrRe = sqrt(sigmasq[1])*cData[0]->data->data[k].re
        + sqrt(sigmasq[2])*cData[1]->data->data[k].re;
      cohSnrIm = sqrt(sigmasq[1])*cData[0]->data->data[k].im
        + sqrt(sigmasq[2])*cData[1]->data->data[k].im;

      params->cohH1H2SNRVec->data->data[k]
        = sqrt( (cohSnrRe*cohSnrRe + cohSnrIm*cohSnrIm) /
                (sigmasq[1] + sigmasq[2] ) );
    }
  }

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

int compare( const void* a, const void* b ) {
  double* arg1 = (double*) a;
  double* arg2 = (double*) b;
  if( *arg1 < *arg2 ) return -1;
  else if( *arg1 == *arg2 ) return 0;
  else return 1;
}

