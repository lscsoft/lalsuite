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
#include <lal/CoherentInspiral.h>


double rint(double x);

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
	  if (params->numBeamPoints != 0)
	    {
	      for ( l = 0 ; l < beamVecPtr->numDetectors ; l++) {  
		for ( i = 0 ; i < detBeamArray[l].numBeamPoints ; i++ ) {
		  TRY( LALDestroyVector(status->statusPtr, &(detBeamArray[l].thetaPhiVs[i].data)), status);
		}
	      }
	      for( l = 0; l < beamVecPtr->numDetectors; l++) {
		LALFree( detBeamArray[l].thetaPhiVs );
	      }
	      LALFree( detBeamArray );
	      beamVecPtr->detBeamArray = NULL;
	      LALFree( beamVecPtr );
	      inputPtr->beamVec = NULL;
	    }
          LALFree( inputPtr );
          *input = NULL;
          ABORT( status, COHERENTINSPIRALH_EALOC, COHERENTINSPIRALH_MSGEALOC);
        }
  
      /* set the number of detectors in the vector */
      cVecPtr->numDetectors = params->numDetectors;
  
      cVecPtr->cData = (COMPLEX8TimeSeries *)
        LALCalloc(1, cVecPtr->numDetectors*sizeof(COMPLEX8TimeSeries));
      if ( !(cVecPtr->cData) )
        {
          LALFree( cVecPtr );
          (*input)->multiCData = NULL;
          if( params->numBeamPoints != 0)
	    {
	      for ( l = 0 ; l < beamVecPtr->numDetectors ; l++) {  
	        for ( i = 0 ; i < detBeamArray[l].numBeamPoints ; i++ ) {
	          TRY( LALDestroyVector(status->statusPtr, &(detBeamArray[l].thetaPhiVs[i].data)), status);
	        }
	      }
	      for( l = 0; l < beamVecPtr->numDetectors; l++) {
	        LALFree( detBeamArray[l].thetaPhiVs );
	      }
	      LALFree( detBeamArray );
	      beamVecPtr->detBeamArray = NULL;
	      LALFree( beamVecPtr );
	      inputPtr->beamVec = NULL;
	    }
          LALFree( inputPtr );
          *input = NULL;
          ABORT( status, COHERENTINSPIRALH_EALOC, COHERENTINSPIRALH_MSGEALOC);
        }
  
      for ( i = 0 ; i < cVecPtr->numDetectors ; i++ ) {
        LALCCreateVector(status->statusPtr, &(cVecPtr->cData[i].data), 
			 params->numPoints);
        BEGINFAIL( status ) {
          LALFree( cVecPtr->cData );
          LALFree( cVecPtr );
          (*input)->multiCData = NULL;
          if( params->numBeamPoints != 0)
	    {
	      for ( l = 0 ; l < beamVecPtr->numDetectors ; l++) {  
	        for ( i = 0 ; i < detBeamArray[l].numBeamPoints ; i++ ) {
	          TRY( LALDestroyVector(status->statusPtr, &(detBeamArray[l].thetaPhiVs[i].data)), status);
	        }
	      }
	      for( l = 0; l < beamVecPtr->numDetectors; l++) {
	        LALFree( detBeamArray[l].thetaPhiVs );
	      }
	      LALFree( detBeamArray );
	      beamVecPtr->detBeamArray = NULL;
	      LALFree( beamVecPtr );
	      inputPtr->beamVec = NULL;
	    }
          LALFree( inputPtr );
          *input = NULL;
        }
        ENDFAIL( status );
      }


  BEGINFAIL( status ) {
    for ( i = 0 ; i < cVecPtr->numDetectors ; i++ ) {
    TRY( LALCDestroyVector(status->statusPtr, &(cVecPtr->cData[i].data)),
	 status);
    }
    LALFree( cVecPtr->cData );
    LALFree( cVecPtr );
    (*input)->multiCData = NULL;
    if( params->numBeamPoints != 0)
      {
	for ( l = 0 ; l < beamVecPtr->numDetectors ; l++) {  
	  for ( i = 0 ; i < detBeamArray[l].numBeamPoints ; i++ ) {
	    TRY( LALDestroyVector(status->statusPtr, &(detBeamArray[l].thetaPhiVs[i].data)), status);
	  }
	}
	for( l = 0; l < beamVecPtr->numDetectors; l++) {
	  LALFree( detBeamArray[l].thetaPhiVs );
	}
	LALFree( detBeamArray );
	beamVecPtr->detBeamArray = NULL;
	LALFree( beamVecPtr );
	inputPtr->beamVec = NULL;
      }
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
  
      for ( l = 0 ; l < cVecPtr->numDetectors ; l++ ) {
	if (cVecPtr->cData[l].data != NULL ) {
	  LALCDestroyVector( status->statusPtr, &(cVecPtr->cData[l].data) );
	  CHECKSTATUSPTR( status );
	}
      }

      LALFree( cVecPtr->cData );
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
  INT4                              networkLength = 6; /*CHECK: hardwired 
							 to 6 for now*/

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
          TRY( LALU2DestroyVector (status->statusPtr, &(outputPtr->detIDVec)),
	       status);
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
   * create detector vector to store detector site IDs
   *
   */
  
  
  outputPtr->detectorVec = (DetectorVector *) 
    LALCalloc( 1, sizeof(DetectorVector) );

  if ( !(outputPtr->detectorVec) ) {
    ABORT( status, COHERENTINSPIRALH_EALOC, COHERENTINSPIRALH_MSGEALOC);
    BEGINFAIL( status ) {
      if( params->cohSNROut )
	{
          TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec->data)),
	       status); 
          LALFree( outputPtr->cohSNRVec );
          outputPtr->cohSNRVec = NULL;
	} 
      TRY( LALU2DestroyVector( status->statusPtr, &(outputPtr->detIDVec)),
	   status);
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
      if( params->cohSNROut )
	{
          TRY( LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec->data)),
	       status);   
          LALFree( outputPtr->cohSNRVec );
          outputPtr->cohSNRVec = NULL;
	}
      TRY( LALU2DestroyVector( status->statusPtr, &(outputPtr->detIDVec)),
	   status);
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
  
  LALFree( outputPtr );
  *output = NULL;


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
  INT4                                caseID[6] = {0,0,0,0,0,0};
  INT4                                i,q,w,m,j,k,l;
  CHAR                                idtag[6][3] = {"H1","L","V","G","T","H2"};
  INT4                                indexarray[4] = {0,0,0,0};
  CHAR                                caseStr[FILENAME_MAX];
  INT8                                chirpTimeNS = 0; 
  UINT4                               numDetectors = 0;
  UINT4                               numSegments = 0;
  UINT4                               numPoints = 0;
  UINT4                               numBeamPoints = 0;
  UINT4                               deltaEventIndex = 0;
  UINT4                               eventStartIdx = 0;
  INT4                                slidePoints[3] = {0,0,0};
  REAL4                               buffer = 0; /*account for timing errors*/
  REAL4                               timingError = 0.00025; /* allowed timing error of 2 ms */
  REAL4                               s[4][3];/*up to 4 distances;in 3D space*/
  REAL4                               deltaT = 0.0;
  REAL4                               nHatVect[3] = {0,0,0};
  REAL4                               distance[4] = {0,0,0,0};
  REAL4                               timeDelay[4] = {0,0,0,0};
  REAL4                               chirpTime = 0;
  REAL4                               cohSNRThresh = 0;
  UINT4                               cohSNROut;
  REAL4                               cohSNRLocal = 0;
  LALDetector                        *detector = NULL;
  COMPLEX8TimeSeries                 *cData = NULL;
  MultiInspiralTable                 *thisEvent = NULL; 
  CoherentInspiralBeamVector         *beamVec = NULL;
  
  
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
  numDetectors = params->numDetectors;
  numPoints = params->numPoints;
  numSegments = params->numSegments;
  numBeamPoints = params->numBeamPoints;
  cohSNRThresh = params->cohSNRThresh;
  cohSNROut = params->cohSNROut;
  deltaT = params->deltaT;
  

  /* if the full coherent snr vector is required, set it to zero */
  if ( cohSNROut ) {
    memset( params->cohSNRVec->data->data, 0, numPoints * sizeof( REAL4 ));
  }

  /*CHECK: hardwired to 6 detectors for now */
  for (l=0 ;l< 6 ; l++)
    {
      caseID[l] = params->detIDVec->data[l];
    }
  
  /* Now generate the network string for MultiInspiralTables */

  /* INT4 indexarray[4];
  for (j=0; j<numDetectors; j++)
    {
      indexarray[j] = 0;
      }*/
  i = 0;  
  for (l=0 ;l<6 ;l++)
    {
      if( caseID[l] )
	{
	  indexarray[i] = l;
	  i++;
	}
    }
  i = 0;
  j = 0;
  k = 0;
  l = 0;
  switch ( numDetectors )
    {
    case 2:
      i = indexarray[0];
      j = indexarray[1];
      LALSnprintf( caseStr, FILENAME_MAX * sizeof(CHAR), "%s-%s",idtag[i],idtag[j]);
      break;

    case 3:
      i=indexarray[0];
      j=indexarray[1];
      k=indexarray[2];
      LALSnprintf( caseStr, FILENAME_MAX * sizeof(CHAR), "%s-%s-%s",idtag[i],idtag[j],idtag[k]);
      break;

    case 4:
      i=indexarray[0];
      j=indexarray[1];
      k=indexarray[2];
      l=indexarray[3];
      LALSnprintf( caseStr, FILENAME_MAX * sizeof(CHAR), "%s-%s-%s-%s",idtag[i],idtag[j],idtag[k],idtag[l]);
      break;
    }
    
  /*** get detector beam-pattern information if we have 3 sites ***/

  if( (numDetectors == 3 && !(caseID[0] && caseID[5])) || numDetectors == 4 )
	{ 
	  beamVec = input->beamVec;
	}
  
  /*** get detector c outputs ***/
  
  /* read in the c-data for multiple detectors */
  cData = input->multiCData->cData;

  /*** get detector-site locations */
  detector = params->detectorVec->detector;
  

  /*Now compute the position vector of all detectors relative to first detector*/

  for ( l=1 ; l < (INT4) numDetectors ; l++) {
    for (i=0;i<3;i++)
      {
	s[l][i] = (REAL4) ( detector[l].location[i] - detector[0].location[i]);
      }
  }
  
  /* calculate the length of the chirp for clustering over chirp-length*/
  {  
    REAL4 eta = input->tmplt->eta;
    REAL4 m1 = input->tmplt->mass1;
    REAL4 m2 = input->tmplt->mass2;
    REAL4 fmin = params->fLow;
    REAL4 mtotal = m1 + m2;
    REAL4 c0 = 5*mtotal*LAL_MTSUN_SI/(256*eta);
    REAL4 c2 = 743.0/252.0 + eta*11.0/3.0;
    REAL4 c3 = -32*LAL_PI/3;
    REAL4 c4 = 3058673.0/508032.0 + eta*(5429.0/504.0 + eta*617.0/72.0);
    REAL4 x  = pow(LAL_PI*mtotal*LAL_MTSUN_SI*fmin, 1.0/3.0);
    REAL4 x2 = x*x;
    REAL4 x3 = x*x2;
    REAL4 x4 = x2*x2;
    REAL4 x8 = x4*x4;
  
    chirpTime = c0*(1 + c2*x2 + c3*x3 + c4*x4)/x8;

    deltaEventIndex = (UINT4) rint( (chirpTime / deltaT) + 1.0 );
    chirpTimeNS = (INT8) (1e9 * chirpTime);
  }

  
  buffer = rint( (timingError/deltaT) + 1.0 );

  /* Now construct the appropriate coherent SNR */ 
  switch (numDetectors)  {
  case 2:
    if(caseID[0] && caseID[5])  /* Network: H1 and H2*/
      {
	m = 0;
	for (k=0;k<(INT4)numPoints;k++)
	  {
	    REAL4 cohSNR = 0.0;

	    for (m=k-buffer; m<k+buffer; m++)
	      {
		if(m >=0 && m < (INT4) numPoints)
		  {
		    cohSNRLocal = sqrt( (cData[0].data->data[k].re + cData[1].data->data[m].re)*(cData[0].data->data[k].re + cData[1].data->data[m].re) + (cData[0].data->data[k].im + cData[1].data->data[m].im)*(cData[0].data->data[k].im + cData[1].data->data[m].im));

		    if(cohSNRLocal > cohSNR)
		      {
			cohSNR = cohSNRLocal;
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
		
		    if ( !thisEvent )
		      {
			ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
		      }
		
		    /* record the data that we need for the clustering algorithm */          
		    thisEvent->end_time.gpsSeconds = k;
		    thisEvent->snr = cohSNR;
		    fflush( stdout ); 

		  } /* done creating a new event */
		  else if (params->maximizeOverChirp && k <= (INT4)(thisEvent->end_time.gpsSeconds + deltaEventIndex) && cohSNR > thisEvent->snr ) {
		    /* if this is the same event, update the maximum */
		    thisEvent->end_time.gpsSeconds = k;
		    thisEvent->snr = cohSNR;
		    strcpy(thisEvent->ifos, caseStr);
		    thisEvent->mass1 = input->tmplt->mass1;
		    thisEvent->mass2 = input->tmplt->mass2;

		    fflush( stdout ); 

		  }
		  else if ( k > (INT4) (thisEvent->end_time.gpsSeconds + deltaEventIndex) || !(params->maximizeOverChirp) ) {
		    /* clean up this event */
		    MultiInspiralTable      *lastEvent = NULL;
		    INT8                    timeNS1;
		    INT4                    timeIndex = thisEvent->end_time.gpsSeconds;
		    timeNS1 = (INT8) (1e9 * timeIndex * deltaT);
		    thisEvent->end_time.gpsSeconds = (INT4) (timeNS1/1000000000L);
		    thisEvent->end_time.gpsNanoSeconds = (INT4) (timeNS1%1000000000L);
		
		    /* store the start of the crossing */
		    eventStartIdx = k;
		
		    /* allocate memory for the newEvent */
		    lastEvent = thisEvent;
		
		    lastEvent->next = thisEvent = (MultiInspiralTable *) 
		      LALCalloc( 1, sizeof(MultiInspiralTable) );
		    if ( !(lastEvent->next) )
		      {
			ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
		      }
		
		    /* stick minimal data into the event */
		    thisEvent->end_time.gpsSeconds = k;
		    thisEvent->snr = cohSNR;
		    strcpy(thisEvent->ifos,caseStr);
		    thisEvent->mass1 = input->tmplt->mass1;
		    thisEvent->mass2 = input->tmplt->mass2;
		   
		  }
		} /* matches if (cohSNR > cohSNRThresh) */

	      }
	  }
	
	/* 
	 *
	 * clean up the last event if there is one
	 * 
	 */
	
	
	if ( thisEvent )
	  {
	    INT8                    timeNS1;
	    INT4                    timeIndex = thisEvent->end_time.gpsSeconds;
	    timeNS1 = (INT8) (1e9 * timeIndex * deltaT);
	    thisEvent->end_time.gpsSeconds = (INT4) (timeNS1/1000000000L);
	    thisEvent->end_time.gpsNanoSeconds = (INT4) (timeNS1%1000000000L);
	    thisEvent->mass1 = input->tmplt->mass1;
	     thisEvent->mass2 = input->tmplt->mass2;
	    
	    fflush( stdout ); 
	  }
      }
    else 
      { /* Network: 2 detectors excluding either H1, H2, or both H1 and H2 */
	/*Here, the time delay looping must start */
	
	/* Now calculate the distance (in meters) */
	distance[1] = sqrt( cartesianInnerProduct(s[1],s[1]) ); 
	timeDelay[1] = distance[1]/LAL_C_SI;
	slidePoints[1] = rint( (fabs(timeDelay[1])/deltaT) + 1.0 );
	
	k = 0;
	q = 0;

	for(k=0;k<(INT4)numPoints;k++)
	  {
	    REAL4 cohSNR = 0.0;
	    
	    for (q = k-slidePoints[1]-buffer; q < k+slidePoints[1]+buffer; q++)
	      {
		if(q >= 0 && q < (INT4) numPoints)
		  {
		    cohSNRLocal = sqrt(cData[0].data->data[k].re*cData[0].data->data[k].re + cData[1].data->data[q].re*cData[1].data->data[q].re + cData[0].data->data[k].im*cData[0].data->data[k].im + cData[1].data->data[q].im*cData[1].data->data[q].im);
		    if(cohSNRLocal > cohSNR)
		      {
			cohSNR = cohSNRLocal; 
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
		thisEvent->end_time.gpsSeconds = k;
		thisEvent->snr = cohSNR;
		fflush( stdout ); 

	      } /* done creating a new event */
	      else if (params->maximizeOverChirp && k <= (INT4) (thisEvent->end_time.gpsSeconds +deltaEventIndex) && cohSNR > thisEvent->snr ) {
		/* if this is the same event, update the maximum */
		thisEvent->end_time.gpsSeconds = k;
		thisEvent->snr = cohSNR;
		strcpy(thisEvent->ifos, caseStr);
		thisEvent->mass1 = input->tmplt->mass1;
		thisEvent->mass2 = input->tmplt->mass2;
		fflush( stdout ); 
	      }
	      else if ( k > (INT4) (thisEvent->end_time.gpsSeconds + deltaEventIndex) || !(params->maximizeOverChirp) ) {
		/* clean up this event */
		MultiInspiralTable  *lastEvent = NULL;
		INT8                    timeNS1;
		INT4                    timeIndex = thisEvent->end_time.gpsSeconds;
		timeNS1 = (INT8) (1e9 * timeIndex * deltaT);
		thisEvent->end_time.gpsSeconds = (INT4) (timeNS1/1000000000L);
		thisEvent->end_time.gpsNanoSeconds = (INT4) (timeNS1%1000000000L);
		
		/* store the start of the crossing */
		eventStartIdx = k;
		
		/* allocate memory for the newEvent */
		lastEvent = thisEvent;		
		lastEvent->next = thisEvent = (MultiInspiralTable *) 
		  LALCalloc( 1, sizeof(MultiInspiralTable) );
		if ( !(lastEvent->next) )
		  {
		    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
		  }
		
		/* stick minimal data into the event */
		thisEvent->end_time.gpsSeconds = k;
		thisEvent->snr = cohSNR;
		strcpy(thisEvent->ifos,caseStr);
		thisEvent->mass1 = input->tmplt->mass1;
		thisEvent->mass2 = input->tmplt->mass2;
		   
	      }
	    } /* matches if (cohSNR > cohSNRThresh) */
	  }

	if ( thisEvent )
	  {
	    INT8                    timeNS1;
	    INT4                    timeIndex = thisEvent->end_time.gpsSeconds;
	    timeNS1 = (INT8) (1e9 * timeIndex * deltaT);	    
	    thisEvent->end_time.gpsSeconds = (INT4) (timeNS1/1000000000L);
	    thisEvent->end_time.gpsNanoSeconds = (INT4) (timeNS1%1000000000L);
	    thisEvent->mass1 = input->tmplt->mass1;
	    thisEvent->mass2 = input->tmplt->mass2;
	    
	    fflush( stdout ); 
	  }
      }
    break;
  case 3: /* Network: 3 detectors including both H1 and H2 */
    if(caseID[0] && caseID[5])
      {    
	/* Now calculate the distance (in meters) */
	distance[1] = sqrt( cartesianInnerProduct( s[1],s[1]) );
	timeDelay[1] = distance[1]/LAL_C_SI;
	slidePoints[1] = rint( (fabs(timeDelay[1])/deltaT) + 1.0 );

	k = 0;
	q = 0;
	m = 0;

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
			    cohSNRLocal = sqrt( ((cData[0].data->data[k].re + cData[2].data->data[m].re)*(cData[0].data->data[k].re + cData[2].data->data[m].re) + (cData[0].data->data[k].im + cData[2].data->data[m].im)*(cData[0].data->data[k].im + cData[2].data->data[m].im)) + cData[1].data->data[q].re*cData[1].data->data[q].re + cData[1].data->data[q].im*cData[1].data->data[q].im);
		    		    		  
			    if(cohSNRLocal > cohSNR)
			      {
				cohSNR = cohSNRLocal;
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
		    thisEvent->end_time.gpsSeconds = k;
		    thisEvent->snr = cohSNR;
		    fflush( stdout ); 

		  } /* done creating a new event */
		  else if (params->maximizeOverChirp && k <= (INT4) (thisEvent->end_time.gpsSeconds +deltaEventIndex) && cohSNR > thisEvent->snr ) {
		    /* if this is the same event, update the maximum */
		    thisEvent->end_time.gpsSeconds = k;
		    thisEvent->snr = cohSNR;
		    strcpy(thisEvent->ifos, caseStr);
		    thisEvent->mass1 = input->tmplt->mass1;
		    thisEvent->mass2 = input->tmplt->mass2;

		    fflush( stdout ); 

		  }
		  else if ( k > (INT4) (thisEvent->end_time.gpsSeconds + deltaEventIndex) || !(params->maximizeOverChirp) ) {
		    /* clean up this event */
		    MultiInspiralTable  *lastEvent = NULL;
		    INT8                    timeNS1;
		    INT4                    timeIndex = thisEvent->end_time.gpsSeconds;
		    timeNS1 = (INT8) (1e9 * timeIndex * deltaT);
		    thisEvent->end_time.gpsSeconds = (INT4) (timeNS1/1000000000L);
		    thisEvent->end_time.gpsNanoSeconds = (INT4) (timeNS1%1000000000L);
		
		    /* store the start of the crossing */
		    eventStartIdx = k;
		
		    /* allocate memory for the newEvent */
		    lastEvent = thisEvent;
		
		    lastEvent->next = thisEvent = (MultiInspiralTable *) 
		      LALCalloc( 1, sizeof(MultiInspiralTable) );
		    if ( !(lastEvent->next) )
		      {
			ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
		      }
		
		    /* stick minimal data into the event */
		    thisEvent->end_time.gpsSeconds = k;
		    thisEvent->snr = cohSNR;
		    strcpy(thisEvent->ifos,caseStr);
		    thisEvent->mass1 = input->tmplt->mass1;
		    thisEvent->mass2 = input->tmplt->mass2;
		   
		  }
		} /* matches if (cohSNR > cohSNRThresh) */
	      }
	  }

	if ( thisEvent )
	  {
	    INT8                    timeNS1;
	    INT4                    timeIndex = thisEvent->end_time.gpsSeconds;
	    timeNS1 = (INT8) (1e9 * timeIndex * deltaT);
	    thisEvent->end_time.gpsSeconds = (INT4) (timeNS1/1000000000L);
	    thisEvent->end_time.gpsNanoSeconds = (INT4) (timeNS1%1000000000L);
	    thisEvent->mass1 = input->tmplt->mass1;
	    thisEvent->mass2 = input->tmplt->mass2;
	    
	    fflush( stdout ); 
	  }
      }
    else
      { /* Network: 3 detectors excluding either H1, H2, or both (H1 && H2)*/
	/*Now the last 3 cases will involve the looping over the coefficients*/

	l = 0;
	k = 0;
	q = 0;
	w = 0;

	for (l=0;l < (INT4) numBeamPoints; l++)
	  {
	    /* position vector to source relative to first detector */
	    nHatVect[0] = cos(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[0] * (REAL4) LAL_PI_180)*sin(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[1] * (REAL4) LAL_PI_180);
	    nHatVect[1] = sin(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[0] * LAL_PI_180)*sin(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[1] * LAL_PI_180);
	    nHatVect[2] = cos(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[1] * LAL_PI_180); 
	    
	    /* Now calculate the distance (in meters) projected along sky-position */
	    distance[1] = cartesianInnerProduct(s[1],nHatVect);
	    distance[2] = cartesianInnerProduct(s[2],nHatVect);
	    timeDelay[1] = distance[1]/LAL_C_SI;
	    timeDelay[2] = distance[2]/LAL_C_SI;
	    slidePoints[1] = rint( (fabs(timeDelay[1])/deltaT) + 1.0 );
	    slidePoints[2] = rint( (fabs(timeDelay[2])/deltaT) + 1.0 );	    
	    
	    for(k=0;k<(INT4)numPoints;k++)
	      {
		REAL4 cohSNR = 0.0;
		for (q = k-slidePoints[1]-buffer;q < k+slidePoints[1]+buffer;q++)
		  {
		    if(q >= 0 && q < (INT4) numPoints)
		      {			
			for (w = q-slidePoints[2]-buffer; w < q+slidePoints[2]+buffer;w++)
			  {
			    if (w >= 0 && w < (INT4) numPoints)
			      {
				cohSNRLocal = sqrt( ( beamVec->detBeamArray[0].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[0].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[0].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[0].thetaPhiVs[l].data->data[3] )*( cData[0].data->data[k].re*cData[0].data->data[k].re + cData[0].data->data[k].im*cData[0].data->data[k].im ) + ( beamVec->detBeamArray[1].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[1].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[1].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[1].thetaPhiVs[l].data->data[3] )*( cData[1].data->data[q].re*cData[1].data->data[q].re + cData[1].data->data[q].im*cData[1].data->data[q].im ) + ( beamVec->detBeamArray[2].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[2].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[3] )*( cData[2].data->data[w].re*cData[2].data->data[w].re + cData[2].data->data[w].im*cData[2].data->data[w].im ) + 2*(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[1].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[0].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[1].thetaPhiVs[l].data->data[3])*(cData[0].data->data[k].re * cData[1].data->data[q].re + cData[0].data->data[k].im * cData[1].data->data[q].im) + 2*(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[0].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[3])*(cData[0].data->data[k].re*cData[2].data->data[w].re + cData[0].data->data[k].im * cData[2].data->data[w].im) + 2*(beamVec->detBeamArray[1].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[1].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[3])*(cData[1].data->data[q].re * cData[2].data->data[w].re + cData[1].data->data[q].im *cData[2].data->data[w].im));
			      
				if(cohSNRLocal > cohSNR)
				  {
				    cohSNR = cohSNRLocal;
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
			thisEvent->end_time.gpsSeconds = k;
		        thisEvent->snr = cohSNR;
		        fflush( stdout ); 

		      } /* done creating a new event */
		      else if (params->maximizeOverChirp && k <= (INT4) (thisEvent->end_time.gpsSeconds +deltaEventIndex) && cohSNR > thisEvent->snr ) {
			/* if this is the same event, update the maximum */
			thisEvent->end_time.gpsSeconds = k;
		        thisEvent->snr = cohSNR;
		        strcpy(thisEvent->ifos, caseStr);
		        thisEvent->mass1 = input->tmplt->mass1;
		        thisEvent->mass2 = input->tmplt->mass2;
			thisEvent->ligo_axis_ra = beamVec->detBeamArray[0].thetaPhiVs[l].data->data[0];
			thisEvent->ligo_axis_dec = beamVec->detBeamArray[0].thetaPhiVs[l].data->data[1];

		        fflush( stdout ); 

		      }
		      else if ( k > (INT4) (thisEvent->end_time.gpsSeconds + deltaEventIndex) || !(params->maximizeOverChirp) ) {
		        /* clean up this event */
		        MultiInspiralTable  *lastEvent = NULL;
		        INT8                    timeNS1;
		        INT4                    timeIndex = thisEvent->end_time.gpsSeconds;
		        timeNS1 = (INT8) (1e9 * timeIndex * deltaT);
		        thisEvent->end_time.gpsSeconds = (INT4) (timeNS1/1000000000L);
		        thisEvent->end_time.gpsNanoSeconds = (INT4) (timeNS1%1000000000L);
		
		        /* store the start of the crossing */
		        eventStartIdx = k;
		
	    	        /* allocate memory for the newEvent */
		        lastEvent = thisEvent;
		
		        lastEvent->next = thisEvent = (MultiInspiralTable *) 
		          LALCalloc( 1, sizeof(MultiInspiralTable) );
		        if ( !(lastEvent->next) )
		          {
			    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
		          }
		
		        /* stick minimal data into the event */
		        thisEvent->end_time.gpsSeconds = k;
		        thisEvent->snr = cohSNR;
		        strcpy(thisEvent->ifos,caseStr);
		        thisEvent->mass1 = input->tmplt->mass1;
		        thisEvent->mass2 = input->tmplt->mass2;
			thisEvent->ligo_axis_ra = beamVec->detBeamArray[0].thetaPhiVs[l].data->data[0];
			thisEvent->ligo_axis_dec = beamVec->detBeamArray[0].thetaPhiVs[l].data->data[1];
		   
		      }
		    } /* matches if (cohSNR > cohSNRThresh) */		    
		  }
	      }

	    if ( thisEvent )
	      {
		INT8                    timeNS1;
	        INT4                    timeIndex = thisEvent->end_time.gpsSeconds;
	        timeNS1 = (INT8) (1e9 * timeIndex * deltaT);
	        thisEvent->end_time.gpsSeconds = (INT4) (timeNS1/1000000000L);
	        thisEvent->end_time.gpsNanoSeconds = (INT4) (timeNS1%1000000000L);
	        thisEvent->mass1 = input->tmplt->mass1;
	        thisEvent->mass2 = input->tmplt->mass2;
	        thisEvent->ligo_axis_ra = beamVec->detBeamArray[0].thetaPhiVs[l].data->data[0];
	        thisEvent->ligo_axis_dec = beamVec->detBeamArray[0].thetaPhiVs[l].data->data[1];
	    
	        fflush( stdout ); 
	      }
	    
	  }/* outer loop end */


      } /* else statement end */
    break;
  case 4: /* Network: 4 detectors including both H1 and H2 */
    if(caseID[0] && caseID[5])
      {
	/*start search looping */
	
	for (l=0;l < (INT4) numBeamPoints; l++) 
	  {
	    /* position vector to source relative to first detector */
	    nHatVect[0] = cos(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[0] * LAL_PI_180)*sin(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[1] * LAL_PI_180);
	    nHatVect[1] = sin(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[0] * LAL_PI_180)*sin(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[1] * LAL_PI_180);
	    nHatVect[2] = cos(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[1] * LAL_PI_180); 
	    
	    /* Now calculate the distance (in meters) projected along sky-position */
	    distance[1] = cartesianInnerProduct(s[1],nHatVect);
	    distance[2] = cartesianInnerProduct(s[2],nHatVect);
	    timeDelay[1] = distance[1]/LAL_C_SI;
	    timeDelay[2] = distance[2]/LAL_C_SI;
	    slidePoints[1] = rint( (fabs(timeDelay[1])/deltaT) + 1.0 );
	    slidePoints[2] = rint( (fabs(timeDelay[2])/deltaT) + 1.0 );
	    	    
	    k = 0;
	    q = 0;
	    w = 0;
	    m = 0;

	    for(k=0;k<(INT4)numPoints;k++)
	      {
		REAL4 cohSNR = 0.0;
		for(m=k-buffer;m<k+buffer;m++)
		  {
		    if(m >= 0 && m < (INT4) numPoints)
		      {
			for (q = m-slidePoints[1]-buffer;q < m+slidePoints[1]+buffer;q++)
			  {
			    if(q >= 0 && q < (INT4) numPoints)
			      {			
				for (w = q-slidePoints[2]-buffer; w < q+slidePoints[2]+buffer;w++)
				  {
				    if (w >= 0 && w < (INT4) numPoints)
				      {
					cohSNRLocal = sqrt( ( beamVec->detBeamArray[0].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[0].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[0].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[0].thetaPhiVs[l].data->data[3] )*( (cData[0].data->data[k].re + cData[3].data->data[m].re)*(cData[0].data->data[k].re + cData[3].data->data[m].re) + (cData[0].data->data[k].im + cData[3].data->data[m].im)*(cData[0].data->data[k].im + cData[3].data->data[m].im) ) + ( beamVec->detBeamArray[1].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[1].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[1].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[1].thetaPhiVs[l].data->data[3] )*( cData[1].data->data[q].re*cData[1].data->data[q].re + cData[1].data->data[q].im*cData[1].data->data[q].im ) + ( beamVec->detBeamArray[2].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[2]  + beamVec->detBeamArray[2].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[3] )*( cData[2].data->data[w].re*cData[2].data->data[w].re + cData[2].data->data[w].im*cData[2].data->data[w].im ) + 2*(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[1].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[0].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[1].thetaPhiVs[l].data->data[3])*(cData[0].data->data[k].re*cData[1].data->data[q].re + cData[0].data->data[k].im*cData[1].data->data[q].im + cData[3].data->data[m].re*cData[1].data->data[q].re + cData[3].data->data[m].im*cData[1].data->data[q].im) + 2*(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[0].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[3])*(cData[0].data->data[k].re*cData[2].data->data[w].re + cData[0].data->data[k].im*cData[2].data->data[w].im + cData[3].data->data[m].re*cData[2].data->data[w].re + cData[3].data->data[m].im*cData[2].data->data[w].im) + 2*(beamVec->detBeamArray[1].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[1].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[3])*(cData[1].data->data[q].re*cData[2].data->data[w].re + cData[1].data->data[q].im*cData[2].data->data[w].im));
			      
					if(cohSNRLocal > cohSNR)
					  {
					    cohSNR = cohSNRLocal;
					  }
				      }
			    
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
			thisEvent->end_time.gpsSeconds = k;
		        thisEvent->snr = cohSNR;
		        fflush( stdout ); 

		      } /* done creating a new event */
		      else if (params->maximizeOverChirp && k <= (INT4) (thisEvent->end_time.gpsSeconds +deltaEventIndex) && cohSNR > thisEvent->snr ) {
			/* if this is the same event, update the maximum */
			thisEvent->end_time.gpsSeconds = k;
		        thisEvent->snr = cohSNR;
		        strcpy(thisEvent->ifos, caseStr);
		        thisEvent->mass1 = input->tmplt->mass1;
		        thisEvent->mass2 = input->tmplt->mass2;
			thisEvent->ligo_axis_ra = beamVec->detBeamArray[0].thetaPhiVs[l].data->data[0];
			thisEvent->ligo_axis_dec = beamVec->detBeamArray[0].thetaPhiVs[l].data->data[1];

		        fflush( stdout ); 

		      }
		      else if ( k > (INT4) (thisEvent->end_time.gpsSeconds + deltaEventIndex) || !(params->maximizeOverChirp) ) {
		        /* clean up this event */
		        MultiInspiralTable  *lastEvent = NULL;
		        INT8                    timeNS1;
		        INT4                    timeIndex = thisEvent->end_time.gpsSeconds;
		        timeNS1 = (INT8) (1e9 * timeIndex * deltaT);
		        thisEvent->end_time.gpsSeconds = (INT4) (timeNS1/1000000000L);
		        thisEvent->end_time.gpsNanoSeconds = (INT4) (timeNS1%1000000000L);
		
		        /* store the start of the crossing */
		        eventStartIdx = k;
		
	    	        /* allocate memory for the newEvent */
		        lastEvent = thisEvent;
		
		        lastEvent->next = thisEvent = (MultiInspiralTable *) 
		          LALCalloc( 1, sizeof(MultiInspiralTable) );
		        if ( !(lastEvent->next) )
		          {
			    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
		          }
		
		        /* stick minimal data into the event */
		        thisEvent->end_time.gpsSeconds = k;
		        thisEvent->snr = cohSNR;
		        strcpy(thisEvent->ifos,caseStr);
		        thisEvent->mass1 = input->tmplt->mass1;
		        thisEvent->mass2 = input->tmplt->mass2;
			thisEvent->ligo_axis_ra = beamVec->detBeamArray[0].thetaPhiVs[l].data->data[0];
			thisEvent->ligo_axis_dec = beamVec->detBeamArray[0].thetaPhiVs[l].data->data[1];
		   
		      }
		    } /* matches if (cohSNR > cohSNRThresh) */		    
		  }
	      }
	    if ( thisEvent )
	      {
	        INT8                    timeNS1;
	        INT4                    timeIndex = thisEvent->end_time.gpsSeconds;
	        timeNS1 = (INT8) (1e9 * timeIndex * deltaT);
	        thisEvent->end_time.gpsSeconds = (INT4) (timeNS1/1000000000L);
	        thisEvent->end_time.gpsNanoSeconds = (INT4) (timeNS1%1000000000L);
	        thisEvent->mass1 = input->tmplt->mass1;
	        thisEvent->mass2 = input->tmplt->mass2;
	        thisEvent->ligo_axis_ra = beamVec->detBeamArray[0].thetaPhiVs[l].data->data[0];
	        thisEvent->ligo_axis_dec = beamVec->detBeamArray[0].thetaPhiVs[l].data->data[1];
	    
	        fflush( stdout ); 
	      }		    
	    
	  } /* end for statement prior to computing distances*/
	

      } /* end outer if statement in case4*/
    else
      { /* Network: 4 detectors excluding either H1, H2, or both H1 and H2 */
	/* there will be one extra loop over the last case since there are 3 nonzero time delays*/
	
	/*start search looping */
	
	for (l=0;l < (INT4) numBeamPoints; l++)
	  {
	    /* position vector to source relative to first detector */
	    nHatVect[0] = cos(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[0] * LAL_PI_180)*sin(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[1] * LAL_PI_180);
	    nHatVect[1] = sin(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[0] * LAL_PI_180)*sin(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[1] * LAL_PI_180);
	    nHatVect[2] = cos(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[1] * LAL_PI_180); 
	    
	    /* Now calculate the distance (in meters) projected along sky-position */
	    distance[1] = cartesianInnerProduct(s[1],nHatVect);
	    distance[2] = cartesianInnerProduct(s[2],nHatVect);
	    distance[3] = cartesianInnerProduct(s[3],nHatVect);
	    timeDelay[1] = distance[1]/LAL_C_SI;
	    timeDelay[2] = distance[2]/LAL_C_SI;
	    timeDelay[3] = distance[3]/LAL_C_SI;
	    slidePoints[1] = rint( (fabs(timeDelay[1])/deltaT) + 1.0 );
	    slidePoints[2] = rint( (fabs(timeDelay[2])/deltaT) + 1.0 );
	    slidePoints[3] = rint( (fabs(timeDelay[3])/deltaT) + 1.0 );	    
	    
	    k = 0;
	    q = 0;
	    w = 0;
	    j = 0;

	    for(k=0;k<(INT4)numPoints;k++)
	      {
		REAL4 cohSNR = 0.0;
		for (q = k-slidePoints[1]-buffer;q < k+slidePoints[1]+buffer;q++)
		  {
		    if(q >= 0 && q < (INT4) numPoints)
		      {
			for (w = q-slidePoints[2]-buffer; w < q+slidePoints[2]+buffer;w++)
			  {
			    if (w >= 0 && w < (INT4) numPoints)
			      {
				for(j = w-slidePoints[3]-buffer; j < w+slidePoints[3]+buffer; j++)
				  {
				    if(j >=0 && j < (INT4) numPoints)
				      {
					cohSNRLocal = sqrt( ( beamVec->detBeamArray[0].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[0].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[0].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[0].thetaPhiVs[l].data->data[3] )*( cData[0].data->data[k].re*cData[0].data->data[k].re + cData[0].data->data[k].im*cData[0].data->data[k].im ) + ( beamVec->detBeamArray[1].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[1].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[1].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[1].thetaPhiVs[l].data->data[3] )*( cData[1].data->data[q].re*cData[1].data->data[q].re + cData[1].data->data[q].im*cData[1].data->data[q].im ) + ( beamVec->detBeamArray[2].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[2].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[3] )*( cData[2].data->data[w].re*cData[2].data->data[w].re + cData[2].data->data[w].im*cData[2].data->data[w].im ) + ( beamVec->detBeamArray[3].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[3].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[3].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[3].thetaPhiVs[l].data->data[3] )*( cData[3].data->data[j].re*cData[3].data->data[j].re + cData[3].data->data[j].im*cData[3].data->data[j].im ) + 2*(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[1].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[0].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[1].thetaPhiVs[l].data->data[3])*(cData[0].data->data[k].re * cData[1].data->data[q].re + cData[0].data->data[k].im * cData[1].data->data[q].im) + 2*(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[0].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[3])*(cData[0].data->data[k].re*cData[2].data->data[w].re + cData[0].data->data[k].im * cData[2].data->data[w].im) + 2*(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[3].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[0].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[3].thetaPhiVs[l].data->data[3])*(cData[0].data->data[k].re*cData[3].data->data[j].re + cData[0].data->data[k].im*cData[3].data->data[j].im) + 2*(beamVec->detBeamArray[1].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[1].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[3])*(cData[1].data->data[q].re * cData[2].data->data[w].re + cData[1].data->data[q].im*cData[2].data->data[w].im) + 2*(beamVec->detBeamArray[1].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[3].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[1].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[3].thetaPhiVs[l].data->data[3])*(cData[1].data->data[q].re*cData[3].data->data[j].re + cData[1].data->data[q].im*cData[3].data->data[j].im) + 2*(beamVec->detBeamArray[2].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[3].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[2].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[3].thetaPhiVs[l].data->data[3])*(cData[2].data->data[w].re*cData[3].data->data[j].re + cData[2].data->data[w].im * cData[3].data->data[j].im));
				      
					if(cohSNRLocal > cohSNR)
					  {
					    cohSNR = cohSNRLocal;
					  }
				      }
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
			thisEvent->end_time.gpsSeconds = k;
		        thisEvent->snr = cohSNR;
		        fflush( stdout ); 

		      } /* done creating a new event */
		      else if (params->maximizeOverChirp && k <= (INT4) (thisEvent->end_time.gpsSeconds +deltaEventIndex) && cohSNR > thisEvent->snr ) {
			/* if this is the same event, update the maximum */
			thisEvent->end_time.gpsSeconds = k;
		        thisEvent->snr = cohSNR;
		        strcpy(thisEvent->ifos, caseStr);
		        thisEvent->mass1 = input->tmplt->mass1;
		        thisEvent->mass2 = input->tmplt->mass2;
			thisEvent->ligo_axis_ra = beamVec->detBeamArray[0].thetaPhiVs[l].data->data[0];
			thisEvent->ligo_axis_dec = beamVec->detBeamArray[0].thetaPhiVs[l].data->data[1];

		        fflush( stdout ); 

		      }
		      else if ( k > (INT4) (thisEvent->end_time.gpsSeconds + deltaEventIndex) || !(params->maximizeOverChirp) ) {
		        /* clean up this event */
		        MultiInspiralTable  *lastEvent = NULL;
		        INT8                    timeNS1;
		        INT4                    timeIndex = thisEvent->end_time.gpsSeconds;
		        timeNS1 = (INT8) (1e9 * timeIndex * deltaT);        
		        thisEvent->end_time.gpsSeconds = (INT4) (timeNS1/1000000000L);
		        thisEvent->end_time.gpsNanoSeconds = (INT4) (timeNS1%1000000000L);
		
		        /* store the start of the crossing */
		        eventStartIdx = k;
		
	    	        /* allocate memory for the newEvent */
		        lastEvent = thisEvent;
		
		        lastEvent->next = thisEvent = (MultiInspiralTable *) 
		          LALCalloc( 1, sizeof(MultiInspiralTable) );
		        if ( !(lastEvent->next) )
		          {
			    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
		          }
		
		        /* stick minimal data into the event */
		        thisEvent->end_time.gpsSeconds = k;
		        thisEvent->snr = cohSNR;
		        strcpy(thisEvent->ifos,caseStr);
		        thisEvent->mass1 = input->tmplt->mass1;
		        thisEvent->mass2 = input->tmplt->mass2;
			thisEvent->ligo_axis_ra = beamVec->detBeamArray[0].thetaPhiVs[l].data->data[0];
			thisEvent->ligo_axis_dec = beamVec->detBeamArray[0].thetaPhiVs[l].data->data[1];
		   
		      }
		    } /* matches if (cohSNR > cohSNRThresh) */
		  }
	      }
	    
	    if ( thisEvent )
	      {
	        INT8                    timeNS1;
	        INT4                    timeIndex = thisEvent->end_time.gpsSeconds;
	        timeNS1 = (INT8) (1e9 * timeIndex * deltaT);
	        thisEvent->end_time.gpsSeconds = (INT4) (timeNS1/1000000000L);
	        thisEvent->end_time.gpsNanoSeconds = (INT4) (timeNS1%1000000000L);
	        thisEvent->mass1 = input->tmplt->mass1;
	        thisEvent->mass2 = input->tmplt->mass2;
	        thisEvent->ligo_axis_ra = beamVec->detBeamArray[0].thetaPhiVs[l].data->data[0];
	        thisEvent->ligo_axis_dec = beamVec->detBeamArray[0].thetaPhiVs[l].data->data[1];
	    
	        fflush( stdout ); 
	      }

	  } /*end the outermost for statement prior to computing distances*/
	
      } /*end else statement */
  } /* end case statement */
  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


