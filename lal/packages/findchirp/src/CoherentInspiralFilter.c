/*----------------------------------------------------------------------- 
 * 
 * File Name: Coherent.c
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
#include <lal/TwoInterfFindChirp.h>
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
  UINT4                            i,j,l;
  UINT4                            length=4; /*length of each thetaPhiVs vector*/
  CoherentInspiralFilterInput     *inputPtr;
  CoherentInspiralBeamVector      *beamVecPtr;
  DetectorBeamArray               *detBeamArray;
  CoherentInspiralZVector         *zVecPtr;
  COMPLEX8TimeSeries              *zData;
  
  LALDetector                     *detector;
  REAL4TimeSeries                 *thetaPhiVs;

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

  /* check that the number of detectors is positive */
  ASSERT( params->numDetectors > 0, status, 
      COHERENTINSPIRALH_ENUMZ, COHERENTINSPIRALH_MSGENUMZ ); 

  /* check that the number of data segments is positive */
  ASSERT( params->numSegments > 0, status, 
      COHERENTINSPIRALH_ENUMZ, COHERENTINSPIRALH_MSGENUMZ );

  /* check that the number of points is positive */
  ASSERT( params->numPoints > 0, status, 
      COHERENTINSPIRALH_ENUMZ, COHERENTINSPIRALH_MSGENUMZ );

  /* check that the number of theta-phi "template" points 
     in the beam-pattern functions file is positive */
  ASSERT( params->numBeamPoints > 0, status, 
      COHERENTINSPIRALH_ENUMZ, COHERENTINSPIRALH_MSGENUMZ );

  inputPtr= *input = (CoherentInspiralFilterInput *)
    LALCalloc(1, sizeof(CoherentInspiralFilterInput) );
  if (! inputPtr )
    {
      ABORT( status, COHERENTINSPIRALH_EALOC, COHERENTINSPIRALH_MSGEALOC);
    }
  

  /*
   *
   * create the CoherentInspiralBeamVector structure
   *
   */
  
  /* allocate memory for an array  */  
  beamVecPtr = inputPtr->beamVec = (CoherentInspiralBeamVector *)
    LALCalloc (1, sizeof(CoherentInspiralBeamVector));
  
  if ( !beamVecPtr )
    {
      LALFree( inputPtr);
      *input = NULL;
      ABORT( status, COHERENTINSPIRALH_EALOC, COHERENTINSPIRALH_MSGEALOC);
    }
  
  /* set the number of detectors in the vector */
  beamVecPtr->numDetectors = params->numDetectors;
  
  detBeamArray = beamVecPtr->detBeamArray = (DetectorBeamArray *)
    LALCalloc (1, inputPtr->beamVec->numDetectors*sizeof(DetectorBeamArray));
  if ( !detBeamArray )
    {
      LALFree( beamVecPtr );
      beamVecPtr = NULL;
      LALFree( inputPtr );
      *input = NULL;
      ABORT( status, COHERENTINSPIRALH_EALOC, COHERENTINSPIRALH_MSGEALOC);
    }
  
  
  /* set the number of beamPoints in the vector */
  for ( l = 0 ; l < beamVecPtr->numDetectors ; l++) {  
    detBeamArray[l].numBeamPoints = params->numBeamPoints;
    
    detBeamArray[l].thetaPhiVs = (REAL4TimeSeries *)
      LALCalloc (1, params->numBeamPoints*sizeof(REAL4TimeSeries));
  }
  
  /* set the number of fields in the thetaPhiVs to 4: theta,phi,v+, and v- */
  for ( l = 0 ; l < beamVecPtr->numDetectors ; l++) {  
    for ( i = 0 ; i < detBeamArray[l].numBeamPoints ; i++ ) {
      LALCreateVector (status->statusPtr, &(detBeamArray[l].thetaPhiVs[i].data), length);
      BEGINFAIL( status ) {
	LALFree( detBeamArray[l].thetaPhiVs );
	detBeamArray[l].thetaPhiVs = NULL;
	LALFree( detBeamArray );
	detBeamArray = NULL;
	LALFree( beamVecPtr );
	beamVecPtr = NULL;
	LALFree( inputPtr );
	*input = NULL;
      }
      ENDFAIL( status ); 
    }
  }
  
  
  /*
   *
   * create the CoherentInspiralZVector structure
   *
   */
  
  /* allocate memory for the zData vector  */  
  zVecPtr = (*input)->multiZData = (CoherentInspiralZVector *)
    LALCalloc (1, sizeof(CoherentInspiralZVector));
  if ( !zVecPtr )
    {
      ABORT( status, COHERENTINSPIRALH_EALOC, COHERENTINSPIRALH_MSGEALOC);
    }
  
  /* set the number of detectors in the vector */
  zVecPtr->numDetectors = params->numDetectors;
  
  zVecPtr->zData = (COMPLEX8TimeSeries *)
    LALCalloc (1, zVecPtr->numDetectors*sizeof(COMPLEX8TimeSeries));
  if ( !zData )
    {
      LALFree( zVecPtr );
      zVecPtr = NULL;
      LALFree( inputPtr );
      *input = NULL;
      ABORT( status, COHERENTINSPIRALH_EALOC, COHERENTINSPIRALH_MSGEALOC);
    }
  
  for ( i = 0 ; i < zVecPtr->numDetectors ; i++ ) {
    LALCCreateVector (status->statusPtr, &(zVecPtr->zData[i].data), 
		      params->numPoints);
    BEGINFAIL( status ) {
      LALFree( zVecPtr );
      zVecPtr = NULL;
      LALFree( inputPtr );
      *input = NULL;
    }
    ENDFAIL( status );
  }
  
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
  UINT4                            i,j,l;
  UINT4                            length=4;
  CoherentInspiralBeamVector      *beamVecPtr;
  DetectorBeamArray               *detBeamArray;
  CoherentInspiralFilterInput     *inputPtr;
  CoherentInspiralZVector         *zVecPtr;
  COMPLEX8TimeSeries              *zData;

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
   * destroy the contents of the CoherentInspiralZVector structure
   *
   */
  
  zVecPtr = (*input)->multiZData;
  
  for ( l = 0 ; l < zVecPtr->numDetectors ; l++ ) {
    if (zVecPtr->zData[l].data != NULL ) {
      LALCDestroyVector( status->statusPtr, &(zVecPtr->zData[l].data) );
      CHECKSTATUSPTR( status );
    }
  }
  LALFree( zVecPtr );
  zVecPtr = NULL;  

  /*
   *
   * destroy the contents of the CoherentInspiralBeamVector structure
   *
   */
  
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
  detBeamArray = NULL;
  LALFree(beamVecPtr);
  beamVecPtr = NULL;
  
  LALFree( inputPtr );
  inputPtr = NULL;

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

  /* check that the number of detectors is positive */
  ASSERT( params->numDetectors > 0, status, 
      COHERENTINSPIRALH_ENUMZ, COHERENTINSPIRALH_MSGENUMZ ); 

  /* check that the number of data segments is positive */
  ASSERT( params->numSegments > 0, status, 
      COHERENTINSPIRALH_ENUMZ, COHERENTINSPIRALH_MSGENUMZ );

  /* check that the number of points is positive */
  ASSERT( params->numPoints > 0, status, 
      COHERENTINSPIRALH_ENUMZ, COHERENTINSPIRALH_MSGENUMZ );

  /* check that the number of theta-phi "template" points 
     in the beam-pattern functions file is positive */
  ASSERT( params->numBeamPoints > 0, status, 
      COHERENTINSPIRALH_ENUMZ, COHERENTINSPIRALH_MSGENUMZ );

 

  /*
   *
   * allocate memory for the FindChirpFilterParams
   *
   */

  
  /* create the output structure */
  outputPtr= *output = (CoherentInspiralFilterParams *)
    LALCalloc(1, sizeof(CoherentInspiralFilterParams) );
  if (! outputPtr )
    {
      ABORT( status, COHERENTINSPIRALH_EALOC, COHERENTINSPIRALH_MSGEALOC);
    }

  LALU2CreateVector (status->statusPtr, &(outputPtr->detIDVec), 
		     networkLength);
  BEGINFAIL( status ) {
    LALFree( outputPtr );
    *output = NULL;
  }
  ENDFAIL( status );
  
  
  
  /*
   *
   * create vector to store coherent SNR, if required
   *
   */

  if ( params->cohSNROut ) {
    outputPtr->cohSNRVec = (REAL4TimeSeries *) 
      LALCalloc( 1, sizeof(REAL4TimeSeries) );
    LALCreateVector (status->statusPtr, &(outputPtr->cohSNRVec->data), 
		     params->numPoints);
    BEGINFAIL( status ) {
      LALFree( outputPtr->detectorVec );
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

  if ( !(outputPtr->detectorVec ) ) {
    LALFree( outputPtr );
    *output = NULL;
  }

  outputPtr->detectorVec->detector = (LALDetector *) 
    LALCalloc( 1, outputPtr->detectorVec->numDetectors*sizeof(LALDetector));
  
  if ( !(outputPtr->detectorVec->detector ) ) {
    LALFree( outputPtr->detectorVec );
    outputPtr->detectorVec = NULL;
    LALFree( outputPtr );
    *output = NULL;
  }
  
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
  LALFree( outputPtr->detectorVec );
  outputPtr->detectorVec = NULL;
  
  LALU2DestroyVector (status->statusPtr, &(outputPtr->detIDVec) ); 
  CHECKSTATUSPTR( status );

  /*
   *
   * destroy coherent SNR vector, if it exists
   *
   */
  if ( outputPtr->cohSNRVec ) {
    LALDestroyVector( status->statusPtr, &(outputPtr->cohSNRVec->data) );
    CHECKSTATUSPTR( status );
  }
  
  LALFree( outputPtr );
  *output = NULL;

}


void
LALCoherentInspiralFilterSegment (
    LALStatus                             *status,
    CoherentInspiralEvent                **eventList,
    CoherentInspiralFilterInput           *input,
    CoherentInspiralFilterParams          *params
    )
{
  INT4                                siteID[6] = {0,1,2,3,4,0}; /*H1 L1 G V T H2*/
  INT4                                caseID[6] = {0,0,0,0,0,0};
  INT4                                i,j,k,l,p,q,w;
  INT4                                dsites[4] = {0,0,0,0};
  UINT4                               numDetectors;
  UINT4                               numSegments;
  UINT4                               numPoints;
  UINT4                               numBeamPoints;
  UINT4                               deltaEventIndex;
  UINT4                               eventStartIdx = 0;
  INT4                                slidePoints[3] = {0,0,0};
  REAL4                               n;
  REAL4                               s[4][3];/* up to 4 distances; in 3D space */
  REAL4                               deltaT = 0.0;
  REAL4                               nHatVect[3] = {0,0,0};
  REAL4                               distance[4] = {0,0,0,0};
  REAL4                               timeDelay[4] = {0,0,0,0};
  REAL4                               chirpTime = 0;
  REAL4                               cohSNRThresh = 0;
  REAL4                               cohSNRLocal = 0;
  REAL4                               cohSNR = 0; /*CHECK remove this */
  BOOLEAN                             cohSNROut;
  LALDetector                        *detector = NULL;
  COMPLEX8TimeSeries                 *zData = NULL;
  CoherentInspiralEvent              *thisEvent = NULL;
  CoherentInspiralZVector            *multiZData = NULL;
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

  /* if a cohSNRVec vector has been created, check we can store data in it  */
  if ( params->cohSNRVec ) {
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
  cohSNROut = params->cohSNROut;
  cohSNRThresh = params->cohSNRThresh;
  deltaT = params->deltaT;
  

  /* if the full coherent snr vector is required, set it to zero */
  if ( params->cohSNRVec ) {
    memset( params->cohSNRVec->data->data, 0, numPoints * sizeof( REAL4 ));
  }

  /*CHECK: hardwired to 6 detectors for now */
  for (l=0 ;l< 6 ; l++)
    {
      caseID[l] = params->detIDVec->data[l];
    }
  
  fprintf(stdout, "You have specified %2d  detector(s).\n",numDetectors);
  fprintf(stdout, "The caseID is: %d %d %d %d %d %d (H1,L1,VIRGO,GEO,TAMA,H2) \n",caseID[0],caseID[1],caseID[2],caseID[3],caseID[4],caseID[5]);
  
  if (numDetectors > 4)
    {
      ABORT( status, COHERENTINSPIRALH_ENDET, COHERENTINSPIRALH_MSGENDET );
      
    }
  if (numDetectors == 0)
    {
      ABORT( status, COHERENTINSPIRALH_ENDET, COHERENTINSPIRALH_MSGENDET );
    }
  
  
  /* read in the coefficients for each specified detector.  This should be moved to the case statements below where the coherent snr is computed.  It should be put in each of the last three cases: 3b,4a,4b.  This saves some unnecessary reading in of data since the coefficients are only used in the last 3 cases. Note that for networks with H1 and H2, there is some redundancy because they have the same coefficients.  currently it reads in a separate file for each of them. */
  
  /*** get detector beam-pattern information ***/
  
  beamVec = input->beamVec;
  
  
  /*** get detector z outputs ***/
  
  /* read in the z-data for multiple detectors */
  zData = input->multiZData->zData;

  
  /*** get detector-site locations */
  detector = params->detectorVec->detector;
  

  /*Now compute the position vector of all detectors relative to first detector*/
  if (numDetectors == 1) {
    ABORT( status, COHERENTINSPIRALH_ENDET, COHERENTINSPIRALH_MSGENDET );
  }
  else {
    for ( l=1 ; l < numDetectors ; l++) {
      for (i=0;i<3;i++)
	{
	  s[l][i] = (REAL4) ( detector[l].location[i] - detector[0].location[i]);
	}
    }
  }
  
  /* calculate the length of the chirp for clustering over chirp-length*/
  {
    REAL4 eta = input->tmplt->eta;
    REAL4 m1 = input->tmplt->mass1;
    REAL4 m2 = input->tmplt->mass2;
    REAL4 fmin = params->fLow;
    REAL4 m = m1 + m2;
    REAL4 c0 = 5*m*LAL_MTSUN_SI/(256*eta);
    REAL4 c2 = 743.0/252.0 + eta*11.0/3.0;
    REAL4 c3 = -32*LAL_PI/3;
    REAL4 c4 = 3058673.0/508032.0 + eta*(5429.0/504.0 + eta*617.0/72.0);
    REAL4 x  = pow(LAL_PI*m*LAL_MTSUN_SI*fmin, 1.0/3.0);
    REAL4 x2 = x*x;
    REAL4 x3 = x*x2;
    REAL4 x4 = x2*x2;
    REAL4 x8 = x4*x4;
    chirpTime = c0*(1 + c2*x2 + c3*x3 + c4*x4)/x8;

    deltaEventIndex = (UINT4) rint( (chirpTime / deltaT) + 1.0 );
  }

  
  /* Now construct the appropriate coherent SNR */
  
  switch (numDetectors)  {
  case 1:
    ABORT( status, COHERENTINSPIRALH_ENDET, COHERENTINSPIRALH_MSGENDET );
    break;
  case 2:
    if(caseID[0] && caseID[5])
      {
	fprintf(stdout,"Network of 2 detectors - includes H1 and H2 \n");
	for (k=0;k<numPoints;k++)
	  {
	    REAL4 cohSNR = (1 / sqrt(2)) * sqrt( pow(zData[0].data->data[k].re + zData[1].data->data[k].re,2) + pow( zData[0].data->data[k].im + zData[1].data->data[k].im,2));
	    params->cohSNRVec->data->data[k] = cohSNR;
	    if ( cohSNR > cohSNRThresh ) {
	      if ( ! *eventList ) {
		/* if this is the first event, start the list */
		thisEvent = *eventList = (CoherentInspiralEvent *) 
		  LALCalloc( 1, sizeof(CoherentInspiralEvent) );
		if ( ! thisEvent )
		  {
		    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
		  }
		thisEvent->timeIndex = k;
		thisEvent->cohSNR = cohSNR;
	      }
	      else if (params->maximiseOverChirp && 
		       k <= thisEvent->timeIndex +deltaEventIndex ) {
		/* if this is the same event, update the maximum */
		thisEvent->timeIndex = k;
		thisEvent->cohSNR = cohSNR;
	      }
	      else if ( k >thisEvent->timeIndex  + deltaEventIndex ||
			! params->maximiseOverChirp )
		{
		  /* clean up this event */
		  CoherentInspiralEvent *lastEvent;
		  
		  thisEvent->timeIndex = k;
		  
		  /* copy the template into the event */
		  thisEvent->mass1  = (REAL4) input->tmplt->mass1;
		  thisEvent->mass2  = (REAL4) input->tmplt->mass2;
		  thisEvent->cohSNR = cohSNR;

		  /* store the start of the crossing */
		  eventStartIdx = k;
		  
		  /* allocate memory for the newEvent */
		  lastEvent = thisEvent;
		  
		  lastEvent->next = thisEvent = (CoherentInspiralEvent *) 
		    LALCalloc( 1, sizeof(CoherentInspiralEvent) );
		  if ( ! lastEvent->next )
		    {
		      ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
		    }
		  
		  /* stick minimal data into the event */
		  thisEvent->timeIndex = k;
		  thisEvent->cohSNR = cohSNR;
		}
	    }
	  } 
      }
    else 
      {
	/*Here, the time delay looping must start */
	fprintf(stdout,"Network of 2 detectors \n");
	
	p=0;
	for (l=0;l<6;l++)
	  {
	    if(caseID[l])
	      {
		dsites[p++] = siteID[l];
	      }
	  }
	
	fprintf(stdout,"%d %d \n",dsites[0],dsites[1]);
	
	for(k=0;k<numPoints;k++)
	  {
	    REAL4 cohSNR = 0.0;

	    /* Now calculate the distance (in meters) */
	    distance[1] = sqrt( cartesianInnerProduct(s[1],s[1]) );
	    
	    timeDelay[1] = distance[1]/LAL_C_SI;
	    slidePoints[1] = ceilf( fabsf(timeDelay[1])/deltaT );
	    
	    for (n = k-slidePoints[1];n < k+slidePoints[1];n++)
	      {
		if(n >= 0 && n < numPoints)
		  {
		    cohSNRLocal = sqrt(pow(zData[0].data->data[k].re,2) + pow(zData[1].data->data[(INT4) n].re,2) + pow(zData[0].data->data[k].im,2) + pow(zData[1].data->data[(INT4) n].im,2));
		    if(cohSNRLocal > cohSNR)
		      {
			cohSNR = cohSNRLocal; 
			/*CHECK: updated timedelays should be stored here too */
		      }
		  }
	      }
	    if ( cohSNR > cohSNRThresh ) {
	      if ( ! *eventList ) {
		/* if this is the first event, start the list */
		thisEvent = *eventList = (CoherentInspiralEvent *) 
		  LALCalloc( 1, sizeof(CoherentInspiralEvent) );
		if ( ! thisEvent )
		  {
		    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
		  }
		thisEvent->timeIndex = k;
		thisEvent->cohSNR = cohSNR;
	      }
	      else if (params->maximiseOverChirp && 
		       k <= thisEvent->timeIndex +deltaEventIndex ) {
		/* if this is the same event, update the maximum */
		thisEvent->timeIndex = k;
		thisEvent->cohSNR = cohSNR;
	      }
	      else if ( k >thisEvent->timeIndex  + deltaEventIndex ||
			! params->maximiseOverChirp )
		{
		  /* clean up this event */
		  CoherentInspiralEvent *lastEvent;
		  
		  thisEvent->timeIndex = k;
		  
		  /* copy the template into the event */
		  thisEvent->mass1  = (REAL4) input->tmplt->mass1;
		  thisEvent->mass2  = (REAL4) input->tmplt->mass2;
		  thisEvent->cohSNR = cohSNR;

		  /* store the start of the crossing */
		  eventStartIdx = j;
		  
		  /* allocate memory for the newEvent */
		  lastEvent = thisEvent;
		  
		  lastEvent->next = thisEvent = (CoherentInspiralEvent *) 
		    LALCalloc( 1, sizeof(CoherentInspiralEvent) );
		  if ( ! lastEvent->next )
		    {
		      ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
		    }
		  
		  /* stick minimal data into the event */
		  thisEvent->timeIndex = k;
		  thisEvent->cohSNR = cohSNR;
		}
	    }
	  }
      }
    break;
  case 3:
    if(caseID[0] && caseID[5])
      {
	fprintf(stdout,"Network of 3 detectors - includes H1 and H2 \n");
	/*The SNR for the H1-H2 pair is first computed, then the looping is done with one time delay */
	
	p=0;
	for (l=0;l<6;l++)
	  {
	    if(caseID[l])
	      {
		dsites[p++] = siteID[l];
	      }
	  }
	
	fprintf(stdout,"%d %d %d \n",dsites[0],dsites[1],dsites[2]);
	
	for(l=0;l<numPoints;l++)
	  {
	    for (n = l-slidePoints[1];n < l+slidePoints[1];n++)
	      {
		if(n >= 0 && n < numPoints)
		  {
		    cohSNRLocal = sqrt( 0.5 * (pow(zData[0].data->data[l].re + zData[1].data->data[l].re,2) + pow(zData[0].data->data[l].im + zData[1].data->data[l].im,2)) + pow(zData[2].data->data[(INT4) n].re,2) + pow(zData[2].data->data[(INT4) n].im,2));
		    
		    
		  }
		if(cohSNRLocal > cohSNR)
		  {
		    cohSNR = cohSNRLocal;
		  }
	      }
	  }
	
      }
    else
      {
	/* Now the last 3 cases will involve the looping over the coefficients */
	p=0;
	for (l=0;l<6;l++)
	  {
	    if(caseID[l])
	      {
		dsites[p++] = siteID[l];
	      }
	  }
	fprintf(stdout,"%d %d %d \n",dsites[0],dsites[1],dsites[2]);
	for (l=0;l < numBeamPoints; l++)
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
	    slidePoints[1] = ceilf( fabsf(timeDelay[1])/deltaT );
	    slidePoints[2] = ceilf( fabsf(timeDelay[2])/deltaT );
	    
	    fprintf(stdout,"%d %d\n",slidePoints[1],slidePoints[2]);
	    
	    for(p=0;p<numPoints;p++)
	      {
		for (n = p-slidePoints[1];n < p+slidePoints[1];n++)
		  {
		    if(n >= 0 && n < numPoints)
		      {
			
			for (q = n-slidePoints[2]; q < n+slidePoints[2];q++)
			  {
			    if (q >= 0 && q < numPoints)
			      {
				cohSNRLocal = sqrt( (pow(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[2],2) + pow(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[3],2))*(pow(zData[0].data->data[p].re,2) + pow(zData[0].data->data[p].im,2)) + (pow(beamVec->detBeamArray[1].thetaPhiVs[l].data->data[2],2) + pow(beamVec->detBeamArray[1].thetaPhiVs[l].data->data[3],2))*(pow(zData[1].data->data[(INT4) n].re,2) + pow(zData[1].data->data[(INT4) n].im,2)) + (pow(beamVec->detBeamArray[2].thetaPhiVs[l].data->data[2],2) + pow(beamVec->detBeamArray[2].thetaPhiVs[l].data->data[3],2))*(pow(zData[2].data->data[(INT4) q].re,2) + pow(zData[2].data->data[(INT4) q].im,2)) + 2*(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[1].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[0].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[1].thetaPhiVs[l].data->data[3])*(zData[0].data->data[p].re * zData[1].data->data[(INT4) n].re + zData[0].data->data[p].im * zData[1].data->data[(INT4) n].im) + 2*(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[0].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[3])*(zData[0].data->data[p].re*zData[2].data->data[(INT4) n].re + zData[0].data->data[p].im * zData[2].data->data[(INT4) n].im) + 2*(beamVec->detBeamArray[1].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[1].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[3])*(zData[1].data->data[(INT4) n].re * zData[2].data->data[(INT4) q].re + zData[1].data->data[(INT4) n].im *zData[2].data->data[(INT4) q].im));
				
				fprintf(stdout,"%f \n", cohSNRLocal);
			      }
			    if(cohSNRLocal > cohSNR)
			      {
				cohSNR = cohSNRLocal; /* CHECK: cohSNR[p]*/
			      }
			    
			  }
			
		      }
		    
		  }
	      }
	    
	  }/* outer loop end */
	
      } /* else statement end */
    break;
  case 4:
    if(caseID[0] && caseID[5])
      {
	fprintf(stdout,"Network of 4 detectors - includes H1 and H2 \n");
	
	p=0;
	for (l=0;l<6;l++)
	  {
	    if(caseID[l])
	      {
		dsites[p++] = siteID[l];
	      }
	  }
	fprintf(stdout,"%d %d %d %d\n",dsites[0],dsites[1],dsites[2],dsites[3]);
	
	/*start search looping */
	
	for (l=0;l < numBeamPoints; l++) 
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
	    slidePoints[1] = ceilf( fabsf(timeDelay[1])/deltaT );
	    slidePoints[2] = ceilf( fabsf(timeDelay[2])/deltaT );
	    slidePoints[3] = ceilf( fabsf(timeDelay[3])/deltaT );
	    
	    fprintf(stdout,"%d %d %d\n",slidePoints[1],slidePoints[2],slidePoints[3]);
	    
	    for(p=0;p<numPoints;p++)
	      {
		for (n = p-slidePoints[2];n < p+slidePoints[2];n++)
		  {
		    if(n >= 0 && n < numPoints)
		      {
			
			for (q = n-slidePoints[3]; q < n+slidePoints[3];q++)
			  {
			    if (q >= 0 && q < numPoints)
			      {
				cohSNRLocal = sqrt( 0.5*(pow(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[2],2) + pow(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[3],2))*(pow(zData[0].data->data[p].re,2) + pow(zData[0].data->data[p].im,2) + pow(zData[1].data->data[p].re,2) + pow(zData[1].data->data[p].im,2) + 2*(zData[0].data->data[p].re*zData[1].data->data[p].re + zData[0].data->data[p].im*zData[1].data->data[p].im)) + (pow(beamVec->detBeamArray[2].thetaPhiVs[l].data->data[2],2) + pow(beamVec->detBeamArray[2].thetaPhiVs[l].data->data[3],2))*(pow(zData[2].data->data[(INT4) n].re,2) + pow(zData[2].data->data[(INT4) n].im,2)) + (pow(beamVec->detBeamArray[3].thetaPhiVs[l].data->data[2],2) + pow(beamVec->detBeamArray[3].thetaPhiVs[l].data->data[3],2))*(pow(zData[3].data->data[(INT4) q].re,2) + pow(zData[3].data->data[(INT4) q].im,2)) + sqrt(2)*(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[0].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[3])*(zData[0].data->data[p].re*zData[2].data->data[(INT4) n].re + zData[0].data->data[p].im*zData[2].data->data[(INT4) n].im + zData[1].data->data[p].re*zData[2].data->data[(INT4) n].re + zData[1].data->data[p].im*zData[2].data->data[(INT4) n].im) + sqrt(2)*(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[3].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[0].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[3].thetaPhiVs[l].data->data[3])*(zData[0].data->data[p].re*zData[3].data->data[(INT4) q].re + zData[0].data->data[p].im*zData[3].data->data[(INT4) q].im + zData[1].data->data[p].re*zData[3].data->data[(INT4) q].re + zData[1].data->data[p].im*zData[3].data->data[(INT4) q].im) + 2*(beamVec->detBeamArray[2].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[3].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[2].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[3].thetaPhiVs[l].data->data[3])*(zData[2].data->data[(INT4) n].re*zData[3].data->data[(INT4) q].re + zData[2].data->data[(INT4) n].im*zData[3].data->data[(INT4) q].im));
				
				fprintf(stdout,"%f \n", cohSNRLocal);
			      }
			    if(cohSNRLocal > cohSNR)
			      {
				cohSNR = cohSNRLocal; /*CHECK: p*/
			      }
			    
			  }
			
		      }
		    
		  }
	      }
	    
	    
	  } /* end for statement prior to computing distances*/
	
	
	
      } /* end outer if statement in case4*/
    else
      {
	fprintf(stdout,"Network of four detctors \n");
	/* there will be one extra loop over the last case since there are 3 nonzero time delays*/
	
	p=0;
	for (l=0;l<6;l++)
	  {
	    if(caseID[l])
	      {
		dsites[p++] = siteID[l];
	      }
	  }
	fprintf(stdout,"%d %d %d %d\n",dsites[0],dsites[1],dsites[2],dsites[3]);
	
	/*start search looping */
	
	for (l=0;l < numBeamPoints; l++)
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
	    slidePoints[1] = ceilf( fabsf(timeDelay[1])/deltaT );
	    slidePoints[2] = ceilf( fabsf(timeDelay[2])/deltaT );
	    slidePoints[3] = ceilf( fabsf(timeDelay[3])/deltaT );
	    
	    fprintf(stdout,"%d %d %d\n",slidePoints[1],slidePoints[2],slidePoints[3]);
	    
	    for(p=0;p<numPoints;p++)
	      {
		for (n = p-slidePoints[1];n < p+slidePoints[1];n++)
		  {
		    if(n >= 0 && n < numPoints)
		      {
			for (q = n-slidePoints[2]; q < n+slidePoints[2];q++)
			  {
			    if (q >= 0 && q < numPoints)
			      {
				for(w = q-slidePoints[3]; w < q+slidePoints[3]; w++)
				  {
				    if(w >=0 && w < numPoints)
				      {
					cohSNRLocal = sqrt( (pow(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[2],2) + pow(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[3],2))*(pow(zData[0].data->data[p].re,2) + pow(zData[0].data->data[p].im,2)) + (pow(beamVec->detBeamArray[1].thetaPhiVs[l].data->data[2],2) + pow(beamVec->detBeamArray[1].thetaPhiVs[l].data->data[3],2))*(pow(zData[1].data->data[(INT4) n].re,2) + pow(zData[1].data->data[(INT4) n].im,2)) + (pow(beamVec->detBeamArray[2].thetaPhiVs[l].data->data[2],2) + pow(beamVec->detBeamArray[2].thetaPhiVs[l].data->data[3],2))*(pow(zData[2].data->data[(INT4) q].re,2) + pow(zData[2].data->data[(INT4) q].im,2)) + (pow(beamVec->detBeamArray[3].thetaPhiVs[l].data->data[2],2) + pow(beamVec->detBeamArray[3].thetaPhiVs[l].data->data[3],2))*(pow(zData[3].data->data[(INT4) w].re,2) + pow(zData[3].data->data[(INT4) w].im,2)) + 2*(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[1].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[0].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[1].thetaPhiVs[l].data->data[3])*(zData[0].data->data[p].re * zData[1].data->data[(INT4) n].re + zData[0].data->data[p].im * zData[1].data->data[(INT4) n].im) + 2*(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[0].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[3])*(zData[0].data->data[p].re*zData[2].data->data[(INT4) n].re + zData[0].data->data[p].im * zData[2].data->data[(INT4) n].im) + 2*(beamVec->detBeamArray[0].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[3].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[0].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[3].thetaPhiVs[l].data->data[3])*(zData[0].data->data[p].re*zData[3].data->data[(INT4) w].re + zData[0].data->data[p].im*zData[3].data->data[(INT4) w].im) + 2*(beamVec->detBeamArray[1].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[1].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[2].thetaPhiVs[l].data->data[3])*(zData[1].data->data[(INT4) n].re * zData[2].data->data[(INT4) q].re + zData[1].data->data[(INT4) n].im *zData[2].data->data[(INT4) q].im) + 2*(beamVec->detBeamArray[1].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[3].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[1].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[3].thetaPhiVs[l].data->data[3])*(zData[1].data->data[(INT4) n].re*zData[3].data->data[(INT4) w].re + zData[1].data->data[(INT4) n].im*zData[3].data->data[(INT4) w].im) + 2*(beamVec->detBeamArray[2].thetaPhiVs[l].data->data[2]*beamVec->detBeamArray[3].thetaPhiVs[l].data->data[2] + beamVec->detBeamArray[2].thetaPhiVs[l].data->data[3]*beamVec->detBeamArray[3].thetaPhiVs[l].data->data[3])*(zData[2].data->data[(INT4) q].re*zData[3].data->data[(INT4) w].re + zData[2].data->data[(INT4) q].im * zData[3].data->data[(INT4) w].im));
					
					
					fprintf(stdout,"%f \n", cohSNRLocal);
				      }
				    if(cohSNRLocal > cohSNR)
				      {
					cohSNR = cohSNRLocal;/*CHECK: p*/
				      }
				  }
			      }
			  }
		      }
		  }
	      }
	    
	  } /*end the outermost for statement prior to computing distances*/
	
      } /*end else statement */
    /*CHECK:      break; */
  } /* end case statement */
  

  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


