/*----------------------------------------------------------------------- 
 * 
 * File Name: TwoInterfFindChirpFilter.c
 *
 * Author: Bose, S., Brown, D. A., and Noel, J. S.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <math.h>
#include <string.h>

#include <lal/LALRCSID.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/DetectorSite.h> 
#include <lal/TwoInterfFindChirp.h>
#include <lal/PrintVector.h>
#include <lal/PrintFTSeries.h>
#include <lal/SkyCoordinates.h>

double rint(double x);

NRCSID (TWOINTERFFINDCHIRPFILTERC, "$Id$");

static REAL4 cartesianInnerProduct(REAL4 x[3], REAL4 y[3])
{
  return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

static void
FindBaseLine( 
	     LALStatus       *status,
             SkyPosition     *baseLine,
	     const LALDetectorPair *detectors,
	     LIGOTimeGPS     *gpsTime
	     )
{
  EarthPosition    earthBaseLine;

  INITSTATUS( status, "FindBaseLine", TWOINTERFFINDCHIRPFILTERC );
  ATTATCHSTATUSPTR( status );
  
  ASSERT( baseLine, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );  
  ASSERT( detectors, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );  
  ASSERT( gpsTime,   status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );  
  
  memset(&earthBaseLine, 0, sizeof(EarthPosition));
  memset(baseLine,       0, sizeof(SkyPosition));
  
  earthBaseLine.x = detectors->detectorTwo.location[0] 
    - detectors->detectorOne.location[0];
  earthBaseLine.y = detectors->detectorTwo.location[1] 
    - detectors->detectorOne.location[1];
  earthBaseLine.z = detectors->detectorTwo.location[2] 
    - detectors->detectorOne.location[2];
  
  baseLine->system  = COORDINATESYSTEM_EQUATORIAL;
  
  TRY( LALGeocentricToGeodetic  ( status->statusPtr, &earthBaseLine ), status );
  TRY( LALGeographicToEquatorial( status->statusPtr, baseLine, 
				  &earthBaseLine.geodetic, gpsTime), status);
  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

void
LALCreateTwoInterfFindChirpInputVector (
					LALStatus                             *status,
					TwoInterfFindChirpFilterInputVector  **vector,
					TwoInterfFindChirpInitParams          *params
					)
{
  UINT4                                    i;
  TwoInterfFindChirpFilterInputVector     *vectorPtr;
  FindChirpFilterInput                    *filterInputPtr;

  INITSTATUS( status, "LALCreateTwoInterfFindChirpFilterInputVector", 
	      TWOINTERFFINDCHIRPFILTERC );
  ATTATCHSTATUSPTR( status );

  /*
   *
   * ensure that the arguments are reasonable
   *
   */
  

  /* check that the output handle exists, but points to a null pointer */
  ASSERT( vector, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );  
  ASSERT( !*vector, status, TWOINTERFFINDCHIRPH_ENNUL, TWOINTERFFINDCHIRPH_MSGENNUL );

  /* check that the parameter structure exists */
  ASSERT( params, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );

  /* check that the number of detectors is positive */
  ASSERT( params->numDetectors > 0, status, 
      TWOINTERFFINDCHIRPH_ENUMZ, TWOINTERFFINDCHIRPH_MSGENUMZ );

  /* check that the number of points is positive */
  ASSERT( params->numPoints > 0, status, 
      TWOINTERFFINDCHIRPH_ENUMZ, TWOINTERFFINDCHIRPH_MSGENUMZ );

  
  /*
   *
   * create a vector of FindChirpFilterInputs
   *
   */

  vectorPtr= *vector = (TwoInterfFindChirpFilterInputVector *)
    LALCalloc(1, sizeof(TwoInterfFindChirpFilterInputVector) );
  if (! vectorPtr )
    {
      ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC);
    }

  /* set the number of detectors in the vector */
  vectorPtr->length = params->numDetectors;
  
  /* allocate memory for an array  */  
  filterInputPtr = vectorPtr->filterInput = (FindChirpFilterInput *)
    LALCalloc (1, vectorPtr->length*sizeof(FindChirpFilterInput));
  if ( !filterInputPtr )
    {
      LALFree( vectorPtr );
      vectorPtr = NULL;
      ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC);
    }

  /*
   *
   * create the twointerffindchirp filter input structure
   *
   */

  for ( i = 0 ; i < vectorPtr->length ; ++i )
    {
      /* create memory for the chirp template structure */
      filterInputPtr[i].fcTmplt = (FindChirpTemplate *)
	LALCalloc( 1, sizeof(FindChirpTemplate) );
      if ( !filterInputPtr[i].fcTmplt )
	{
	  LALFree( *vector );
	  *vector = NULL;
	  ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC );
	}
      
      /* create memory for the chirp template data */
      LALCCreateVector (status->statusPtr, &(filterInputPtr[i].fcTmplt->data), 
			params->numPoints);
      BEGINFAIL( status )
	{
	  LALFree( filterInputPtr[i].fcTmplt );
	  filterInputPtr[i].fcTmplt = NULL;
	  LALFree( *vector );
	  *vector = NULL;
	}
      ENDFAIL( status );
    }
  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

void
LALDestroyTwoInterfFindChirpInputVector (
    LALStatus                            *status,
    TwoInterfFindChirpFilterInputVector **vector
    )
{
  UINT4                                    i;
  TwoInterfFindChirpFilterInputVector     *vectorPtr;
  FindChirpFilterInput                    *filterInput;
  
  INITSTATUS( status, "LALDestroyTwoInterfFindChirpFilterInputVector", TWOINTERFFINDCHIRPFILTERC );
  ATTATCHSTATUSPTR( status );
  
  
  /* 
   * check that the arguments are reasonable
   *
   */
  
  ASSERT( vector, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( *vector, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  
  vectorPtr= *vector;
  
  filterInput = (*vector)->filterInput;

  /*
   *
   * destroy the contents of the array
   *
   */
  
  for ( i = 0 ; i < (*vector)->length ; ++i )
    {
      
      /*
       *
       * destroy the findchirp input structure
       *
       */
      
      if ( filterInput[i].fcTmplt->data != NULL )
	{
	  /* destroy the chirp template data storage */
	  LALCDestroyVector( status->statusPtr, &(filterInput[i].fcTmplt->data) );
	  CHECKSTATUSPTR( status );
	}
      if ( filterInput[i].fcTmplt != NULL )
	{
	  /* destroy the chirp template structure */
	  LALFree(filterInput[i].fcTmplt );
	  filterInput[i].fcTmplt = NULL;
	}
    }
  
  /* destroy the filter input structure */
  LALFree( filterInput );
  filterInput = NULL;
  
  /* free the vector structure */
  LALFree( vectorPtr );
  *vector = NULL;
  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

void
LALTwoInterfFindChirpFilterInit (
    LALStatus                                      *status,
    TwoInterfFindChirpFilterParams                **output,
    TwoInterfFindChirpInitParams                   *params
    )
{
  UINT4                                    i;
  TwoInterfFindChirpFilterParams          *outputPtr;
  TwoInterfFindChirpFilterParamsVector    *vectorPtr;
  FindChirpFilterParams                   *filterParamsPtr;

  INITSTATUS( status, "LALTwoInterfFindChirpFilterInit", TWOINTERFFINDCHIRPFILTERC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */
  
  
  /* check that the output handle exists, but points to a null pointer */
  ASSERT( output, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( !*output, status, TWOINTERFFINDCHIRPH_ENNUL, TWOINTERFFINDCHIRPH_MSGENNUL );
  
  /* check that the parameter structure exists */
  ASSERT( params, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );

  /* check that the number of points is positive */
  ASSERT( params->numPoints > 0,  status, 
	  TWOINTERFFINDCHIRPH_ENUMZ, TWOINTERFFINDCHIRPH_MSGENUMZ );
  
  /* create the output structure */
  outputPtr = *output = (TwoInterfFindChirpFilterParams *)
    LALCalloc( 1, sizeof(TwoInterfFindChirpFilterParams) );
  if ( ! outputPtr )
    {
      ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC );
    }
  
  vectorPtr= (*output)->paramsVec = (TwoInterfFindChirpFilterParamsVector *)
    LALCalloc(1, sizeof(TwoInterfFindChirpFilterParamsVector) );
  if (! vectorPtr )
    {
      ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC);
    }

  /* set the number of detectors in the vector */
  vectorPtr->length = params->numDetectors;
  /* allocate memory for an array  */  
  filterParamsPtr = vectorPtr->filterParams = (FindChirpFilterParams *)
    LALCalloc (1, vectorPtr->length*sizeof(FindChirpFilterParams));
  if ( !filterParamsPtr )
    {
      LALFree( vectorPtr );
      LALFree( outputPtr );
      *output = NULL;
      ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC);
    }

  /*
   *
   * create the twointerffindchirp filter input structure
   *
   */

  for ( i = 0 ; i < vectorPtr->length ; ++i )
    {  
      
      /* create memory for the chisq parameters */
      filterParamsPtr[i].chisqParams = (FindChirpChisqParams *)
	LALCalloc( 1, sizeof(FindChirpChisqParams) );
      if ( !filterParamsPtr[i].chisqParams )
	{
	  LALFree( outputPtr->twoInterfRhosqVec );
	  LALFree( vectorPtr );
	  LALFree( outputPtr );
	  *output = NULL;
	  ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC );
	}
      
      /* create memory for the chisq input */
      filterParamsPtr[i].chisqInput = (FindChirpChisqInput *)
	LALCalloc( 1, sizeof(FindChirpChisqInput) );
      if ( !filterParamsPtr[i].chisqInput )
	{
	  LALFree( outputPtr->twoInterfRhosqVec );
	  LALFree( filterParamsPtr[i].chisqParams );
	  LALFree( vectorPtr );
	  LALFree( outputPtr );
	  *output = NULL;
	  ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC );
	}
      
      
      /*
       *
       * create fft plans and workspace vectors
       *
       */
      
      
      /* create plan for optimal filter */
      
      LALCreateReverseComplexFFTPlan( status->statusPtr,
				      &(filterParamsPtr[i].invPlan), params->numPoints, 0 );
      BEGINFAIL( status )
	{
	  LALFree( outputPtr->twoInterfRhosqVec );
	  LALFree( filterParamsPtr[i].chisqInput );
	  LALFree( filterParamsPtr[i].chisqParams );
	  LALFree( vectorPtr );
	  LALFree( outputPtr );
	  *output = NULL;
	}
      ENDFAIL( status );
      
      /* create workspace vector for optimal filter: time domain */
      LALCCreateVector( status->statusPtr, &(filterParamsPtr[i].qVec), 
			params->numPoints );
      BEGINFAIL( status )
	{
	  TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
					 &(filterParamsPtr[i].invPlan) ), status );
	  
	  LALFree( outputPtr->twoInterfRhosqVec );
	  LALFree( filterParamsPtr[i].chisqInput );
	  LALFree( filterParamsPtr[i].chisqParams );
	  LALFree( vectorPtr );
	  LALFree( outputPtr );
	  *output = NULL;
	}
      ENDFAIL( status );
      
      /* create workspace vector for optimal filter: freq domain */
      LALCCreateVector( status->statusPtr, &(filterParamsPtr[i].qtildeVec), 
			params->numPoints );
      BEGINFAIL( status )
	{
	  TRY( LALCDestroyVector( status->statusPtr, &(filterParamsPtr[i].qVec) ), 
	       status );
	  
	  TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
					 &(filterParamsPtr[i].invPlan) ), status );
	  LALFree( outputPtr->twoInterfRhosqVec );
	  LALFree( filterParamsPtr[i].chisqInput );
	  LALFree( filterParamsPtr[i].chisqParams );
	  LALFree( vectorPtr );
	  LALFree( outputPtr );
	  *output = NULL;
	}
      ENDFAIL( status );
      
      /* create workspace vector for chisq filter */
      LALCreateVector (status->statusPtr, &(filterParamsPtr[i].chisqVec), 
		       params->numPoints);
      BEGINFAIL( status )
	{
	  TRY( LALCDestroyVector( status->statusPtr, &(filterParamsPtr[i].qtildeVec) ), 
	       status );
	  TRY( LALCDestroyVector( status->statusPtr, &(filterParamsPtr[i].qVec) ), 
	       status );
	  
	  TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
					 &(filterParamsPtr[i].invPlan) ), status );
	  LALFree( outputPtr->twoInterfRhosqVec );
	  LALFree( filterParamsPtr[i].chisqInput );
	  LALFree( filterParamsPtr[i].chisqParams );
	  LALFree( vectorPtr );
	  LALFree( outputPtr );
	  *output = NULL;
	}
      ENDFAIL( status );
  
      /*
       *
       * create vector to store snrsq, if required
       *
       */
      
      if ( params->createRhosqVec )
	{
	  filterParamsPtr[i].rhosqVec = (REAL4TimeSeries *) 
	    LALCalloc( 1, sizeof(REAL4TimeSeries) );
	  LALCreateVector (status->statusPtr, &(filterParamsPtr[i].rhosqVec->data), 
			   params->numPoints);
	  BEGINFAIL( status )
	    {
	      TRY( LALDestroyVector (status->statusPtr, &(filterParamsPtr[i].chisqVec) ), 
		   status ); 
	      TRY( LALCDestroyVector( status->statusPtr, &(filterParamsPtr[i].qtildeVec) ), 
		   status );
	      TRY( LALCDestroyVector( status->statusPtr, &(filterParamsPtr[i].qVec) ), 
		   status );
	      
	      TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
					     &(filterParamsPtr[i].invPlan) ), status );
	      
	      LALFree( outputPtr->twoInterfRhosqVec );
	      LALFree( filterParamsPtr[i].chisqInput );
	      LALFree( filterParamsPtr[i].chisqParams );
	      LALFree( vectorPtr );
	      LALFree( outputPtr );
	      *output = NULL;
	    }
	  ENDFAIL( status );
	}
    }/* end loop over detectors */
  
  
  /*
   * create vector to store network snrsq, if required
   *
   */
  
  
  if ( params->createTwoInterfRhosqVec )
    {
      outputPtr->twoInterfRhosqVec = (REAL4TimeSeries *) 
	LALCalloc( 1, sizeof(REAL4TimeSeries) );
      LALCreateVector (status->statusPtr, &(outputPtr->twoInterfRhosqVec->data), 
		       params->numPoints);
      BEGINFAIL( status )
	{
	  for ( i = 0 ; i < vectorPtr->length ; ++i )
	    {  
	      
	      TRY( LALDestroyVector( status->statusPtr, 
				     &(filterParamsPtr[i].rhosqVec->data)), status ); 
	      TRY( LALDestroyVector (status->statusPtr, &(filterParamsPtr[i].chisqVec) ), 
		   status ); 
	      TRY( LALCDestroyVector( status->statusPtr, &(filterParamsPtr[i].qtildeVec) ), 
		   status );
	      TRY( LALCDestroyVector( status->statusPtr, &(filterParamsPtr[i].qVec) ), 
		   status );
	      
	      TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
					     &(filterParamsPtr[i].invPlan) ), status );
	      LALFree( filterParamsPtr[i].chisqInput );
	      LALFree( filterParamsPtr[i].chisqParams );
	      
	      LALFree( vectorPtr );
	      LALFree( outputPtr );
	      *output = NULL;
	    }
	}/* end loop over detectors */
      
      ENDFAIL( status );
    }    
  
  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

void
LALTwoInterfFindChirpFilterFinalize (
    LALStatus                                   *status,
    TwoInterfFindChirpFilterParams             **output
    )
{
  UINT4                                    i;
  TwoInterfFindChirpFilterParams          *outputPtr;
  FindChirpFilterParams                   *paramsPtr;
  
  INITSTATUS( status, "LALTwoInterfFindChirpFilterFinalize", 
	      TWOINTERFFINDCHIRPFILTERC );
  ATTATCHSTATUSPTR( status );
  
  
  /*
   *
   * check that the arguments are reasonable
   *
   */

  
  ASSERT( output, status, TWOINTERFFINDCHIRPH_ENULL, 
	  TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( *output, status, TWOINTERFFINDCHIRPH_ENULL, 
	  TWOINTERFFINDCHIRPH_MSGENULL );
   
  /*local pointer to output structure */
  outputPtr = *output;
  
  /* destroy vector to store network snrsq, if it exists */
  if ( outputPtr->twoInterfRhosqVec )
    {
      LALDestroyVector (status->statusPtr, 
			&(outputPtr->twoInterfRhosqVec->data)); 
      CHECKSTATUSPTR( status );
      
      LALFree(outputPtr->twoInterfRhosqVec);
    }  
  
  /*
   *
   * destroy the network FilterParams Vector
   *
   */
  
  paramsPtr= (*output)->paramsVec->filterParams;
  
  for ( i = 0 ; i < (*output)->paramsVec->length ; ++i )
    {
      /*
       * free chisq structures 
       */      

      /* parameter structure */
      LALFree( paramsPtr[i].chisqParams );

      /* input structure */
      LALFree( paramsPtr[i].chisqInput );
      
      /* 
       * destroy fft plans and workspace vectors
       */

      LALDestroyComplexFFTPlan( status->statusPtr, &(paramsPtr[i].invPlan));
      CHECKSTATUSPTR( status );
      
      /* create workspace vector for optimal filter: time domain */
      LALCDestroyVector( status->statusPtr, &(paramsPtr[i].qVec));
      CHECKSTATUSPTR( status );
      
      /* create workspace vector for optimal filter: freq domain */
      LALCDestroyVector( status->statusPtr, &(paramsPtr[i].qtildeVec));
      CHECKSTATUSPTR( status );
      
      /* create workspace vector for chisq filter */
      LALDestroyVector (status->statusPtr, &(paramsPtr[i].chisqVec));
      CHECKSTATUSPTR( status );
      
      /*
       *
       * create vector to store snrsq, if required
       *
       */

      if ( paramsPtr[i].rhosqVec )
	{
	  LALDestroyVector (status->statusPtr, &(paramsPtr[i].rhosqVec->data)); 
	  CHECKSTATUSPTR( status );

	  LALFree(paramsPtr[i].rhosqVec);
	}    
    }
  
  /* free the vector structure */
  LALFree( paramsPtr );
  LALFree( outputPtr->paramsVec );
  LALFree( outputPtr );
  *output = NULL;

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

void
LALTwoInterfFindChirpFilterSegment (
    LALStatus                             *status,
    TwoInterfInspiralEvent               **eventList,
    TwoInterfFindChirpFilterInputVector   *input,
    TwoInterfFindChirpFilterParams        *params
    )
{
  /* UINT4                          j, m, n; */

  UINT4 j, maxLagPts, n;

  UINT4                          k=0;
  UINT4                          startIndex1, startIndex2;
  UINT4                          endIndex1, endIndex2;
  UINT4                          deltaEventIndex;
  UINT4                          twoInterfEventId = 0; 
  UINT4                          ignoreIndex[2];
  INT8                           chirpTimeNS;
  REAL4                          deltaT; 
  REAL4                          maxLag; /*lightTravel time between dets 1,2*/
  REAL4                          s[3];
  REAL4                          distance;
  REAL4                          norm1;
  REAL4                          norm2;
  REAL4                          modqsq1;
  REAL4                          modqsq2;
  REAL4                          rhosq = 0;
  REAL4                          twoInterfRhosqThresh;
  BOOLEAN                        haveChisq1    = 0;
  BOOLEAN                        haveChisq2    = 0;
  UINT4                          numPoints; 
  COMPLEX8                      *q1tilde      = NULL;
  COMPLEX8                      *q1           = NULL;
  COMPLEX8                      *inputData1   = NULL;
  COMPLEX8                      *q2tilde      = NULL;
  COMPLEX8                      *q2           = NULL;
  COMPLEX8                      *inputData2   = NULL;
  COMPLEX8                      *tmpltSignal1  = NULL;
  COMPLEX8                      *tmpltSignal2  = NULL;
  TwoInterfInspiralEvent        *thisEvent     = NULL; 
  BOOLEAN                        unresolvedEvent = 0;

  SkyPosition                   baseline;

  INITSTATUS( status, "LALTwoInterFindChirpFilter", 
	      TWOINTERFFINDCHIRPFILTERC );
  ATTATCHSTATUSPTR( status );
  
  
  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* check that the output handle exists, but points to a null pointer */
  ASSERT( eventList, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( !*eventList, status, 
	  TWOINTERFFINDCHIRPH_ENNUL, TWOINTERFFINDCHIRPH_MSGENNUL );
  
  /* check that the parameter structure exists */
  ASSERT( params, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  
  /* check that the parameter sub-structures exist */

  ASSERT( params->paramsVec, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  
  /* check that the filter parameters are reasonable */
  ASSERT( params->twoInterfRhosqThresh > 0, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );

  ASSERT( params->paramsVec->length > 0, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );

  ASSERT( params->twoInterfRhosqThresh > 0, status,
	  TWOINTERFFINDCHIRPH_ERHOT, TWOINTERFFINDCHIRPH_MSGERHOT );

  ASSERT( params->detectors != NULL, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );


  ASSERT( params->paramsVec->filterParams[0].deltaT > 0, status,
	  TWOINTERFFINDCHIRPH_EDTZO, TWOINTERFFINDCHIRPH_MSGEDTZO );
  ASSERT( params->paramsVec->filterParams[1].deltaT > 0, status,
	  TWOINTERFFINDCHIRPH_EDTZO, TWOINTERFFINDCHIRPH_MSGEDTZO );
  ASSERT( params->paramsVec->filterParams[0].chisqThresh > 0, status,
	  TWOINTERFFINDCHIRPH_ECHIT, TWOINTERFFINDCHIRPH_MSGECHIT ); 
  ASSERT( params->paramsVec->filterParams[1].chisqThresh > 0, status,
	  TWOINTERFFINDCHIRPH_ECHIT, TWOINTERFFINDCHIRPH_MSGECHIT ); 
  
  

  /* check that the fft plan exists */
  ASSERT( params->paramsVec->filterParams[0].invPlan, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( params->paramsVec->filterParams[1].invPlan, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  
  /* check that the workspace vectors exist for detector 1 */
  ASSERT( params->paramsVec->filterParams[0].qVec, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( params->paramsVec->filterParams[0].qVec->data, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( params->paramsVec->filterParams[0].qtildeVec, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( params->paramsVec->filterParams[0].qtildeVec->data, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  
  /* check that the workspace vectors exist for detector 2 */
  ASSERT( params->paramsVec->filterParams[1].qVec, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( params->paramsVec->filterParams[1].qVec->data, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( params->paramsVec->filterParams[1].qtildeVec, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( params->paramsVec->filterParams[1].qtildeVec->data, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  
  /* check that the chisq parameter and input structures exist */
  ASSERT( params->paramsVec->filterParams[0].chisqParams, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( params->paramsVec->filterParams[0].chisqInput, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( params->paramsVec->filterParams[1].chisqParams, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( params->paramsVec->filterParams[1].chisqInput, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  
  /* if a rhosqVec vector has been created, check we can store data in it */
  if ( params->paramsVec->filterParams[0].rhosqVec ) 
    {
      ASSERT( params->paramsVec->filterParams[0].rhosqVec->data->data, status, 
	      TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
      ASSERT( params->paramsVec->filterParams[0].rhosqVec->data, status, 
	      TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
    }
  if ( params->paramsVec->filterParams[1].rhosqVec ) 
    {
      ASSERT( params->paramsVec->filterParams[1].rhosqVec->data->data, status, 
	      TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
      ASSERT( params->paramsVec->filterParams[1].rhosqVec->data, status, 
	      TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
    }
 

  /* if a chisqVec vector has been created, check we can store data in it */
  if ( params->paramsVec->filterParams[0].chisqVec ) 
    {
      ASSERT( params->paramsVec->filterParams[0].chisqVec->data, status,
	      TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
    }
  if ( params->paramsVec->filterParams[1].chisqVec ) 
    {
      ASSERT( params->paramsVec->filterParams[1].chisqVec->data, status,
	      TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
    }

  /* check that the input structure exists */
  ASSERT( input, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( input->length > 0, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  
  /* check that the input structure contains some input */
  ASSERT( input->filterInput[0].tmplt, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( input->filterInput[0].fcTmplt, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( input->filterInput[0].segment, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( input->filterInput[0].segment, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL ); /* introduced this line */
  
  ASSERT( input->filterInput[1].tmplt, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( input->filterInput[1].fcTmplt, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( input->filterInput[1].segment, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( input->filterInput[1].segment, status, 
	  TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL ); /* introduced this line */
  
  /*
   *
   * point local pointers to input and output pointers
   *
   */

  
  /* workspace vectors for detector 1 */
  q1 = params->paramsVec->filterParams[0].qVec->data;
  q1tilde = params->paramsVec->filterParams[0].qtildeVec->data;
  
  /* workspace vectors for detector 2 */
  q2 = params->paramsVec->filterParams[1].qVec->data;
  q2tilde = params->paramsVec->filterParams[1].qtildeVec->data;
  
  /* template and data */
  inputData1 = input->filterInput[0].segment->data->data->data; 
  inputData2 = input->filterInput[1].segment->data->data->data;
  tmpltSignal1 = input->filterInput[0].fcTmplt->data->data;
  tmpltSignal2 = input->filterInput[1].fcTmplt->data->data;
  
  /* sampling rate of detector 1 data should match that in detector 2 */
  if ( !(params->paramsVec->filterParams[0].deltaT ==
	 params->paramsVec->filterParams[0].deltaT)) 
    {
      ABORT( status, TWOINTERFFINDCHIRPH_EDTZO, TWOINTERFFINDCHIRPH_MSGEDTZO );
    }
  deltaT = params->paramsVec->filterParams[0].deltaT; 

  /* calculate separation vector between sites */
  for ( j=0; j<3; j++) {
    s[j] = (REAL4) ( params->detectors->detectorOne.location[j] -
		     params->detectors->detectorTwo.location[j]);
  }
  /* calculate distance between sites (in meters) */
  distance = sqrt( cartesianInnerProduct(s,s) );
  
  /* light-travel time between sites (in seconds) */
  maxLag = distance / LAL_C_SI;
  maxLagPts = (UINT4) ( maxLag / deltaT ) ;
 
  /* number of points in detector 1 segment should match that in detector 2 */
  if ( !(params->paramsVec->filterParams[0].qVec->length == 
	 params->paramsVec->filterParams[1].qVec->length)) 
    {
      ABORT( status, TWOINTERFFINDCHIRPH_ENUMZ, TWOINTERFFINDCHIRPH_MSGENUMZ );
    }
  numPoints = params->paramsVec->filterParams[0].qVec->length;
  
  /*
   *
   * compute viable search regions in the snrsq vector
   *
   */


  {
    /* calculate the length of the chirp */
    REAL4 eta = input->filterInput[0].tmplt->eta;
    REAL4 m1 = input->filterInput[0].tmplt->mass1;
    REAL4 m2 = input->filterInput[0].tmplt->mass2;
    REAL4 fmin = input->filterInput[0].segment->fLow;
    REAL4 m =  m1+m2;
    REAL4 c0 = 5*m*LAL_MTSUN_SI/(256*eta);
    REAL4 c2 = 743.0/252.0 + eta*11.0/3.0;
    REAL4 c3 = -32*LAL_PI/3;
    REAL4 c4 = 3058673.0/508032.0 + eta*(5429.0/504.0 + eta*617.0/72.0);
    REAL4 x  = pow(LAL_PI*m*LAL_MTSUN_SI*fmin, 1.0/3.0);
    REAL4 x2 = x*x;
    REAL4 x3 = x*x2;
    REAL4 x4 = x2*x2;
    REAL4 x8 = x4*x4;
    REAL4 chirpTime = c0*(1 + c2*x2 + c3*x3 + c4*x4)/x8;

    deltaEventIndex = (UINT4) rint( (chirpTime / deltaT) + 1.0 );
    chirpTimeNS = (INT8) (1e9 * chirpTime);

    /* ignore corrupted data at start and end */
    for ( n = 0 ; n < 2 ; ++n )
      {
	ignoreIndex[n] = ( input->filterInput[n].segment->invSpecTrunc / 2 ) + deltaEventIndex;
      }
    
#if 0
    fprintf( stdout, "m1 = %e m2 = %e => %e seconds => %d points\n"
	     "invSpecTrunc = %d => ignoreIndex in detector 1 = %d\n" 
	     "ignoreIndex in detector 2 = %d\n", 
	     m1, m2, chirpTime, deltaEventIndex, 
	     input->filterInput[0].segment->invSpecTrunc, ignoreIndex[0]
	     input->filterInput[1].segment->invSpecTrunc, ignoreIndex[1] );
    fflush( stdout );
#endif
    
    /* XXX check that we are not filtering corrupted data XXX */
    /* XXX this is hardwired to 1/4 segment length        XXX */
    for ( n = 0 ; n < 2 ; ++n )
      {
	if ( ignoreIndex[n] > numPoints / 4 )
	  {
	    ABORT( status, TWOINTERFFINDCHIRPH_ECRUP, TWOINTERFFINDCHIRPH_MSGECRUP );
	  }
	/* XXX reset ignoreIndex to one quarter of a segment XXX */
	ignoreIndex[n] = numPoints / 4;
      }
  }
  
  /*timeIndex ranges of data being filtered*/
  fprintf( stdout, "filtering detector 1 data from timeIndex %d to %d\n",
	   startIndex1 = ignoreIndex[0], endIndex1 = numPoints - ignoreIndex[0]);
  fprintf( stdout, "filtering detector 2 data from timeIndex %d to %d\n",
	   startIndex2 = ignoreIndex[1], endIndex2 = numPoints - ignoreIndex[1]);
  fflush( stdout );
  
  /*
   *
   * compute q1tilde and q1
   *
   */
  
  
  memset( q1tilde, 0, numPoints * sizeof(COMPLEX8) );
  
  /* q1tilde positive frequency, not DC or nyquist */
  for ( k = 1; k < numPoints/2; ++k )
    {
      REAL4 r1 = inputData1[k].re;
      REAL4 s1 = inputData1[k].im;
      REAL4 x = tmpltSignal1[k].re;
      REAL4 y = 0 - tmpltSignal1[k].im;       /* note complex conjugate */
      
      q1tilde[k].re = r1*x - s1*y;
      q1tilde[k].im = r1*y + s1*x;
    }
  
  /* q1tilde negative frequency only: not DC or nyquist */
  if ( params->paramsVec->filterParams[0].computeNegFreq )
    {
      for ( k = numPoints/2 + 2; k < numPoints - 1; ++k )
	{
	  REAL4 r1 = inputData1[k].re;
	  REAL4 s1 = inputData1[k].im;
	  REAL4 x = tmpltSignal1[k].re;
	  REAL4 y = 0 - tmpltSignal1[k].im;     /* note complex conjugate */
	  
	  q1tilde[k].re = r1*x - s1*y;
	  q1tilde[k].im = r1*y + s1*x;
	}
    }
  LALCOMPLEX8VectorFFT( status->statusPtr, params->paramsVec->filterParams[0].qVec, 
			params->paramsVec->filterParams[0].qtildeVec, 
			params->paramsVec->filterParams[0].invPlan );
  CHECKSTATUSPTR( status );
  
  /* inverse fft to get q1 */
  
  /*
   *
   * compute q2tilde and q2
   *
   */
  
  
  memset( q2tilde, 0, numPoints * sizeof(COMPLEX8) );
  
  /* q2tilde positive frequency, not DC or nyquist */
  for ( k = 1; k < numPoints/2; ++k )
    {
      REAL4 r2 = inputData2[k].re;
      REAL4 s2 = inputData2[k].im;
      REAL4 x = tmpltSignal2[k].re;
      REAL4 y = 0 - tmpltSignal2[k].im;       /* note complex conjugate */
      
      q2tilde[k].re = r2*x - s2*y;
      q2tilde[k].im = r2*y + s2*x;
    }
  
  /* q2tilde negative frequency only: not DC or nyquist */
  if ( params->paramsVec->filterParams[1].computeNegFreq )
    {
      for ( k = numPoints/2 + 2; k < numPoints - 1; ++k )
	{
	  REAL4 r2 = inputData2[k].re;
	  REAL4 s2 = inputData2[k].im;
	  REAL4 x = tmpltSignal2[k].re;
	  REAL4 y = 0 - tmpltSignal2[k].im;     /* note complex conjugate */
	  
	  q2tilde[k].re = r2*x - s2*y;
	  q2tilde[k].im = r2*y + s2*x;
	}
    }
  
  /* inverse fft to get q2 */
  LALCOMPLEX8VectorFFT( status->statusPtr, params->paramsVec->filterParams[1].qVec, 
			params->paramsVec->filterParams[1].qtildeVec, 
			params->paramsVec->filterParams[1].invPlan );
  CHECKSTATUSPTR( status );
 
  /* 
   *
   * calculate signal to noise squared 
   *
   */
  
  /* if full snrsq vector for the network is required, set it to zero */
  if ( params->twoInterfRhosqVec )
    { 
      memset( params->twoInterfRhosqVec->data->data, 0, 
	      numPoints * sizeof( REAL4 ) );
    }
  
  /* if full snrsq vector for each detector is required, set it to zero */
  for ( n = 0 ; n < 2 ; ++n )
    {
      if ( params->paramsVec->filterParams[n].rhosqVec )
	memset( params->paramsVec->filterParams[n].rhosqVec->data->data, 0, 
		numPoints * sizeof( REAL4 ) );
    }/*end loop over detectors*/
  
  /* normalization */
  norm1 = 4.0 * (deltaT / (REAL4)numPoints) / input->filterInput[0].segment->segNorm;
  norm2 = 4.0 * (deltaT / (REAL4)numPoints) / input->filterInput[1].segment->segNorm;
  
 
  /* if full snrsq vector is required, store the snrsq */
  if ( params->paramsVec->filterParams[0].rhosqVec ) 
    {
      LALSnprintf( params->paramsVec->filterParams[0].rhosqVec->name, LALNameLength * sizeof(CHAR),
		"%s:rhosq:output", input->filterInput[0].segment->data->name );
      memcpy( &(params->paramsVec->filterParams[0].rhosqVec->epoch), &(input->filterInput[0].segment->data->epoch), 
	      sizeof(LIGOTimeGPS) );
      params->paramsVec->filterParams[0].rhosqVec->deltaT = input->filterInput[0].segment->deltaT;
      
      for ( j = 0; j < numPoints; ++j )
	{
	  REAL4 modqsq = q1[j].re * q1[j].re + q1[j].im * q1[j].im;
	  params->paramsVec->filterParams[0].rhosqVec->data->data[j] = norm1 * modqsq;
	}
    }
  
  if ( params->paramsVec->filterParams[1].rhosqVec ) 
    {
      snprintf( params->paramsVec->filterParams[1].rhosqVec->name, LALNameLength * sizeof(CHAR),
		"%s:rhosq:output", input->filterInput[1].segment->data->name );
      memcpy( &(params->paramsVec->filterParams[1].rhosqVec->epoch), &(input->filterInput[1].segment->data->epoch), 
	      sizeof(LIGOTimeGPS) );
      params->paramsVec->filterParams[1].rhosqVec->deltaT = input->filterInput[1].segment->deltaT;
      
      for ( j = 0; j < numPoints; ++j )
	{
	  REAL4 modqsq = q2[j].re * q2[j].re + q2[j].im * q2[j].im;
	  params->paramsVec->filterParams[1].rhosqVec->data->data[j] = norm2 * modqsq;
	}
    }
  
  twoInterfRhosqThresh = params->twoInterfRhosqThresh;

/* look for an events in the filter output */
  for ( j = startIndex1; j < endIndex1; ++j ) 
    {
      for ( k = ((j < startIndex2+maxLagPts) ? startIndex2 : (j-maxLagPts)) ; k < (( (j + maxLagPts+1) < endIndex2) ? (j + maxLagPts+1) : endIndex2) ; ++k )
	{
	  REAL4                  rhosq1;
	  REAL4                  rhosq2;
	  
	  modqsq1 = ( q1[j].re * q1[j].re + q1[j].im * q1[j].im );
	  modqsq2 = ( q2[k].re * q2[k].re + q2[k].im * q2[k].im );
	  
	  rhosq1  = norm1 * modqsq1;
	  rhosq2  = norm2 * modqsq2;
	  
	  rhosq = ( rhosq1 + rhosq2 )/sqrt(2);

	  /* if full snrsq vector is required, store the snrsq */
	  
	  if ( params->twoInterfRhosqVec ) 
	    {
	      snprintf( params->twoInterfRhosqVec->name, 
			LALNameLength * sizeof(CHAR),
			"%s:twoInterfRhosq:output", 
			input->filterInput[0].segment->data->name );
	      params->twoInterfRhosqVec->deltaT 
		= input->filterInput[0].segment->deltaT;
	      if ( rhosq > params->twoInterfRhosqVec->data->data[j] )
		params->twoInterfRhosqVec->data->data[j] = rhosq;
	      if ( ignoreIndex[0] > ignoreIndex[1] )
		{
		  memcpy( &(params->twoInterfRhosqVec->epoch), 
			  &(input->filterInput[0].segment->data->epoch), 
			  sizeof(LIGOTimeGPS) );
		}
	      else
		{
		   memcpy( &(params->twoInterfRhosqVec->epoch), 
			   &(input->filterInput[1].segment->data->epoch), 
			   sizeof(LIGOTimeGPS) );
		}
	    }
	  
	  /* if snrsq exceeds threshold at any point */
	  if ( rhosq > twoInterfRhosqThresh
	       && rhosq1 > params->paramsVec->filterParams[0].rhosqThresh
	       && rhosq2 > params->paramsVec->filterParams[1].rhosqThresh
	       )
	    {
	      /* compute chisq vector for detector 1 if it does not exist */
	      if ( ! haveChisq1   && input->filterInput[0].segment->chisqBinVec->length )
		{
		  memset( params->paramsVec->filterParams[0].chisqVec->data, 0, 
			  params->paramsVec->filterParams[0].chisqVec->length 
			  * sizeof(REAL4) );

		  /* pointers to chisq input for detector 1 */
		  params->paramsVec->filterParams[0].chisqInput->qtildeVec 
		    = params->paramsVec->filterParams[0].qtildeVec;
		  params->paramsVec->filterParams[0].chisqInput->qVec      
		    = params->paramsVec->filterParams[0].qVec;
		  
		  /* pointer to the chisq bin vector in segment of detector1 */
		  params->paramsVec->filterParams[0].chisqParams->chisqBinVec 
		    = input->filterInput[0].segment->chisqBinVec;
		  params->paramsVec->filterParams[0].chisqParams->norm  
		    = norm1;
		  /* XXX this should be passed in from the bank XXX */;
#if 0
		  params->paramsVec->filterParams[0].chisqParams->bankMatch=0.97;
#endif	  
		  /* compute the chisq threshold for detector 1 */

		  LALTwoInterfFindChirpChisqVeto ( status->statusPtr, 
					  params->paramsVec->filterParams[0].chisqVec, 
					  params->paramsVec->filterParams[0].chisqInput, 
					  params->paramsVec->filterParams[0].chisqParams );
		  CHECKSTATUSPTR (status); 
		  
		  haveChisq1 =1;
		}		  

	      /* compute chisq vector for detector 2 if it does not exist*/
	      if ( ! haveChisq2   && input->filterInput[1].segment->chisqBinVec->length )
		{
		  memset( params->paramsVec->filterParams[1].chisqVec->data, 0, 
			  params->paramsVec->filterParams[1].chisqVec->length 
			  * sizeof(REAL4) );
		  
		  /* pointers to chisq input for detector 2 */
		  params->paramsVec->filterParams[1].chisqInput->qtildeVec 
		    = params->paramsVec->filterParams[1].qtildeVec;
		  params->paramsVec->filterParams[1].chisqInput->qVec      
		    = params->paramsVec->filterParams[1].qVec;
		  
		  /* pointer to the chisq bin vector in segment of detector2 */
		  params->paramsVec->filterParams[1].chisqParams->chisqBinVec 
		    = input->filterInput[1].segment->chisqBinVec;
		  params->paramsVec->filterParams[1].chisqParams->norm   
		    = norm2;
		  /* XXX this should be passed in from the bank XXX */;
#if 0
		  params->paramsVec->filterParams[1].chisqParams->bankMatch=0.97;
#endif		  
		  /* compute the chisq threshold for detector 2 */
		  
		  LALTwoInterfFindChirpChisqVeto ( status->statusPtr, 
					  params->paramsVec->filterParams[1].chisqVec, 
					  params->paramsVec->filterParams[1].chisqInput, 
					  params->paramsVec->filterParams[1].chisqParams );
		  CHECKSTATUSPTR (status); 
		  
		  haveChisq2 =1;
		}
	      
	      /* if we don't have a chisq or the chisq drops below 
		 threshold, then start processing events */
	      if ( (! input->filterInput[0].segment->chisqBinVec->length ||
		    params->paramsVec->filterParams[0].chisqVec->data[j] 
		    < params->paramsVec->filterParams[0].chisqThresh ) 
		   && (! input->filterInput[1].segment->chisqBinVec->length ||
		       params->paramsVec->filterParams[1].chisqVec->data[k] 
		       < params->paramsVec->filterParams[1].chisqThresh ) )
		{

		  if ( unresolvedEvent ) {
		    /*  current rhosq peak falls within boundaries of last chirp, so it is treated as being
		     *  the same event 
                     */

		    if (rhosq > thisEvent->snrsq) {
		      /* if this is the same event, update the maximum */
		      thisEvent->eventIn1->timeIndex = j;
		      thisEvent->eventIn1->snrsq     = rhosq1;
		      thisEvent->eventIn2->timeIndex = k;
		      thisEvent->eventIn2->snrsq     = rhosq2;
		      
		      thisEvent->timeIndex = j;		      
		      thisEvent->snrsq = rhosq;
		    }
		    /* if the rhosq value was less than previous values for this event, it is ignored */
		  }

		  else { /* no events left unresolved => this is a new event */

		    if ( ! *eventList )
		    {
		      /* if this is the first event, start the list */
		      thisEvent = *eventList = (TwoInterfInspiralEvent *) 
			LALCalloc( 1, sizeof(TwoInterfInspiralEvent) );
		      if ( ! thisEvent )
			{
			  ABORT( status, TWOINTERFFINDCHIRPH_EALOC, 
				 TWOINTERFFINDCHIRPH_MSGEALOC );
			}

		    }
		    else { /* there are previous events (though they are resolved) */
		     
		      /* advance the list pointer to a new event structure */
		      thisEvent->next = (TwoInterfInspiralEvent *) 
			LALCalloc( 1, sizeof(TwoInterfInspiralEvent) );
		      if ( ! thisEvent->next )
			{
			  ABORT( status, TWOINTERFFINDCHIRPH_EALOC, 
				 TWOINTERFFINDCHIRPH_MSGEALOC );
			}
		      
		      thisEvent = thisEvent->next;
		    }

		    /* record event number */
		    thisEvent->twoInterfId = twoInterfEventId++;

		    /* allocate memory for associated single detector events */
		    thisEvent->eventIn1 = (InspiralEvent *) 
		      LALCalloc( 1, sizeof(InspiralEvent) );

		    thisEvent->eventIn2 = (InspiralEvent *) 
		      LALCalloc( 1, sizeof(InspiralEvent) );

		    if ( ! ( thisEvent->eventIn1 && thisEvent->eventIn2) )
		      {
			ABORT( status, TWOINTERFFINDCHIRPH_EALOC, 
			       TWOINTERFFINDCHIRPH_MSGEALOC );
		      }

		    /* record the segment id that the event was found in */
		    thisEvent->eventIn1->segmentNumber 
		      = input->filterInput[0].segment->number;
		    thisEvent->eventIn2->segmentNumber 
		      = input->filterInput[1].segment->number;
		    
		    thisEvent->segmentNumber
		      = input->filterInput[0].segment->number;
		    
		    /* stick minimal data into the event */
		    thisEvent->eventIn1->timeIndex = j;
		    thisEvent->eventIn1->snrsq     = rhosq1;
		    thisEvent->eventIn2->timeIndex = k;
		    thisEvent->eventIn2->snrsq     = rhosq2;
		    
		    thisEvent->timeIndex = j;
		    thisEvent->snrsq = rhosq;

		    /* set flag to finalize this chirp later */
		    unresolvedEvent = 1;

		  } /* done creating a new event */

		} /* done processing event which passed chissq test */

	    } /* done processing event which passed rhosq test */

	  if ( unresolvedEvent ) {
	    
	    /* determine whether we are done processing the current chirp */
	    if ( ! (params->paramsVec->filterParams[0].maximiseOverChirp 
		    && params->paramsVec->filterParams[1].maximiseOverChirp)
		 ||
		 
		 ( (j == thisEvent->eventIn1->timeIndex + deltaEventIndex ||
		    j == endIndex1 - 1 /* was j == endIndex1 */ )
		   && 
		   ( 
		    /* k == thisEvent->eventIn2->timeIndex + deltaEventIndex || */
		    k == ( ((j + maxLagPts+1) < endIndex2) ? (j + maxLagPts) : endIndex2-1)
		    /* was k == endIndex2 */
		    )
		   )
		 )
	      {
		/* finalize the event fields for the current chirp */

		INT8                    timeNS1;
		INT8                    timeNS2;
		
		/* set the event LIGO GPS time of the event 
		   for detector 1 */
		timeNS1 = 1000000000L * 
		  (INT8) (input->filterInput[0].segment->data->epoch.gpsSeconds);
		timeNS1 += 
		  (INT8) (input->filterInput[0].segment->data->epoch.gpsNanoSeconds);
		timeNS1 += 
		  (INT8) (1e9 * (thisEvent->eventIn1->timeIndex)*deltaT);
		thisEvent->eventIn1->time.gpsSeconds = 
		  (INT4) (timeNS1/1000000000L);
		thisEvent->eventIn1->time.gpsNanoSeconds = 
		  (INT4) (timeNS1%1000000000L);
		
		thisEvent->time.gpsSeconds 
		  = (INT4) (timeNS1/1000000000L);
		thisEvent->time.gpsNanoSeconds 
		  = (INT4) (timeNS1%1000000000L);		      
		
		/* set the event LIGO GPS time of the event 
		   for detector 2 */
		timeNS2 = 1000000000L * 
		  (INT8) (input->filterInput[1].segment->data->epoch.gpsSeconds);
		timeNS2 += 
		  (INT8) (input->filterInput[1].segment->data->epoch.gpsNanoSeconds);
		timeNS2 += 
		  (INT8) (1e9 * (thisEvent->eventIn2->timeIndex)*deltaT);
		thisEvent->eventIn2->time.gpsSeconds = 
		  (INT4) (timeNS2/1000000000L);
		thisEvent->eventIn2->time.gpsNanoSeconds = 
		  (INT4) (timeNS2%1000000000L);
		
		/* set the impulse time for the event */
		timeNS1 -= chirpTimeNS;
		thisEvent->eventIn1->impulseTime.gpsSeconds = (INT4) (timeNS1/1000000000L);
		thisEvent->eventIn1->impulseTime.gpsNanoSeconds = (INT4) (timeNS1%1000000000L);
		
		thisEvent->impulseTime.gpsSeconds = (INT4) (timeNS1/1000000000L);
		thisEvent->impulseTime.gpsNanoSeconds = (INT4) (timeNS1%1000000000L);

		timeNS2 -= chirpTimeNS;
		thisEvent->eventIn2->impulseTime.gpsSeconds = (INT4) (timeNS2/1000000000L);
		thisEvent->eventIn2->impulseTime.gpsNanoSeconds = (INT4) (timeNS2%1000000000L);
		
		/* record the ifo and channel name for the event */
		strncpy( thisEvent->eventIn1->ifoName, input->filterInput[0].segment->data->name, 
			 2 * sizeof(CHAR) );
		strncpy( thisEvent->eventIn1->channel, input->filterInput[0].segment->data->name + 3,
			 (LALNameLength - 3) * sizeof(CHAR) );
		
		strncpy( thisEvent->eventIn2->ifoName, input->filterInput[1].segment->data->name, 
			 2 * sizeof(CHAR) );
		strncpy( thisEvent->eventIn2->channel, input->filterInput[1].segment->data->name + 3,
			 (LALNameLength - 3) * sizeof(CHAR) );
		
		/* copy the template into the event */
		memcpy( &(thisEvent->tmplt), input->filterInput[0].tmplt, 
			sizeof(InspiralTemplate) );
		thisEvent->tmplt.next = NULL;
		thisEvent->tmplt.fine = NULL;
		
		/* set snrsq, chisq, sigma and effDist for this event */
		if ( input->filterInput[0].segment->chisqBinVec->length )
		  {
		    thisEvent->eventIn1->chisq =  
		      params->paramsVec->filterParams[0].chisqVec->data[thisEvent->eventIn1->timeIndex];
		    thisEvent->eventIn1->numChisqBins =  
		      input->filterInput[0].segment->chisqBinVec->length - 1;
		  }
		else
		  {
		    thisEvent->eventIn1->chisq = 0;
		    thisEvent->eventIn1->numChisqBins = 0;
		  }

		if ( input->filterInput[1].segment->chisqBinVec->length)
		  
		  {
		    thisEvent->eventIn2->chisq = 
		      params->paramsVec->filterParams[1].chisqVec->data[thisEvent->eventIn2->timeIndex];  
		    thisEvent->eventIn2->numChisqBins = 
		      input->filterInput[1].segment->chisqBinVec->length - 1;
		  }
		else
		  {
		    thisEvent->eventIn2->chisq = 0;
		    thisEvent->eventIn2->numChisqBins = 0;
		  }
		thisEvent->eventIn1->sigma= sqrt(norm1* 
						 input->filterInput[0].segment->segNorm*
						 input->filterInput[0].segment->segNorm*
						 input->filterInput[0].fcTmplt->tmpltNorm );
		thisEvent->eventIn2->sigma= sqrt(norm2* 
						 input->filterInput[1].segment->segNorm*
						 input->filterInput[1].segment->segNorm*
						 input->filterInput[1].fcTmplt->tmpltNorm );
		thisEvent->eventIn1->effDist = 
		  (input->filterInput[0].fcTmplt->tmpltNorm * input->filterInput[0].segment->segNorm * 
		   input->filterInput[0].segment->segNorm) * norm1 / thisEvent->eventIn1->snrsq;
		thisEvent->eventIn1->effDist = sqrt( thisEvent->eventIn1->effDist );
		
		thisEvent->eventIn2->effDist = 
		  (input->filterInput[1].fcTmplt->tmpltNorm * input->filterInput[1].segment->segNorm * 
		   input->filterInput[1].segment->segNorm) * norm2/ thisEvent->eventIn2->snrsq;
		thisEvent->eventIn2->effDist = sqrt( thisEvent->eventIn2->effDist );
		
		/* Compute orientation of the detector-pair's baseline */
		
		TRY (FindBaseLine( status->statusPtr, &baseline,
				   params->detectors, &thisEvent->time), status);

		thisEvent->twoInterfAxisRa  = (REAL4) baseline.longitude;
		thisEvent->twoInterfAxisDec = (REAL4) baseline.latitude;

		thisEvent->sigma =thisEvent->eventIn1->sigma;
		thisEvent->effDist =thisEvent->eventIn1->effDist;
		thisEvent->timeIndex2 =thisEvent->eventIn2->timeIndex;
		
		thisEvent->chisq1 = thisEvent->eventIn1->chisq;
		thisEvent->chisq2 = thisEvent->eventIn2->chisq;
	
		{
		  INT4  delayPts = (INT4) thisEvent->timeIndex2 - (INT4)thisEvent->timeIndex;
		  REAL4 delay = delayPts * deltaT;

		  thisEvent->twoInterfAngle   = acos(LAL_C_SI * delay / distance);
		    
		}

		unresolvedEvent = 0;

	      } /* done finalizing/resolving the event */

	  } /* done with case of having an unresolved event */

	} /* done looping over k (second interferometer time index) */

    } /* done looping over j (first interferometer time index) */

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
  
}
