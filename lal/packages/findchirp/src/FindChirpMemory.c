/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpMemory.c
 *
 * Author: Brown D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <lal/LALStdlib.h>
#include <lal/DataBuffer.h>
#include <lal/FindChirp.h>

NRCSID (FINDCHIRPMEMORYC, "$Id$");

void
LALCreateDataSegmentVector (
    LALStatus                  *status,
    DataSegmentVector         **vector,
    FindChirpInitParams        *params
    )
{
  UINT4                 i;
  DataSegmentVector    *vectorPtr;
  DataSegment          *segPtr;

  INITSTATUS( status, "LALCreateDataSegmentVector", FINDCHIRPMEMORYC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */
  

  /* make sure the output handle exists, but points to a null pointer */
  ASSERT( vector, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );
  ASSERT( !*vector, status, FINDCHIRP_ENNUL, FINDCHIRP_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT( params, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );

  /* make sure that the number of segments and number of points in a */
  /* segment are positive                                            */
  ASSERT( params->numSegments > 0, status, 
      FINDCHIRP_ESEGZ, FINDCHIRP_MSGESEGZ );
  ASSERT( params->numPoints > 0, status, 
      FINDCHIRP_ENUMZ, FINDCHIRP_MSGENUMZ );


  /*
   *
   * create a vector of DataSegments
   *
   */
  

  /* create the ouput structure */
  vectorPtr = *vector = (DataSegmentVector *) 
    LALMalloc( sizeof(DataSegmentVector) );
  ASSERT( vectorPtr, status, FINDCHIRP_EALOC, FINDCHIRP_MSGEALOC );
  memset( vectorPtr, 0, sizeof(DataSegmentVector) );

  /* set the number of segments in the vector */
  vectorPtr->length = params->numSegments;

  /* allocate memory for an array of data segments */
  segPtr = vectorPtr->data = (DataSegment *) 
    LALMalloc( vectorPtr->length * sizeof(DataSegment) );
  if ( !segPtr )
  {
    LALFree( vectorPtr );
    vectorPtr = NULL;
    ABORT( status, FINDCHIRP_EALOC, FINDCHIRP_MSGEALOC );
  }
  memset( segPtr, 0, vectorPtr->length * sizeof(DataSegment) );

  /* create the individual DataSegments in the vector */
  /* just give up if anything fails                   */
  /* this should be cleaned up to prevent leaks       */
  for (i = 0; i < vectorPtr->length; ++i)
  {
    /* ifodmro */
    segPtr[i].data = (INT2TimeSeries *) 
      LALMalloc( sizeof(INT2TimeSeries) );
    ASSERT( segPtr[i].data, status, FINDCHIRP_EALOC, FINDCHIRP_MSGEALOC );
    memset( segPtr[i].data, 0, sizeof(INT2TimeSeries) );

    LALI2CreateVector (status->statusPtr, 
        &segPtr[i].data->data, params->numPoints);
    CHECKSTATUSPTR (status);

    /* power spectrum */
    segPtr[i].spec = (REAL4FrequencySeries *) 
      LALMalloc( sizeof(REAL4FrequencySeries) );
    ASSERT( segPtr[i].spec, status, FINDCHIRP_EALOC, FINDCHIRP_MSGEALOC );
    memset( segPtr[i].spec, 0, sizeof(REAL4FrequencySeries) );

    LALCreateVector (status->statusPtr, 
        &segPtr[i].spec->data, params->numPoints/2 + 1);
    CHECKSTATUSPTR (status);
    
    /* response function */
    segPtr[i].resp = (COMPLEX8FrequencySeries *) 
      LALMalloc( sizeof(COMPLEX8FrequencySeries) );
    ASSERT( segPtr[i].resp, status, FINDCHIRP_EALOC, FINDCHIRP_MSGEALOC );
    memset( segPtr[i].resp, 0, sizeof(COMPLEX8FrequencySeries) );

    LALCCreateVector (status->statusPtr, 
        &segPtr[i].resp->data, params->numPoints/2 + 1);
    CHECKSTATUSPTR (status);
  }


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}



void
LALDestroyDataSegmentVector (
    LALStatus                  *status,
    DataSegmentVector         **vector
    )
{
  UINT4                 i;
  DataSegment          *segPtr;

  INITSTATUS( status, "LALDestroyDataSegmentVector", FINDCHIRPMEMORYC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */
  

  /* make sure handle is non-null and points to a non-null pointer */
  ASSERT( vector, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );
  ASSERT( *vector, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );


  /*
   *
   * destroy the data segment vector
   *
   */
  
  
  /* local pointer to the segment array */
  segPtr = (*vector)->data;

  /* destroy the contents of the array */
  for (i = 0; i < (*vector)->length; ++i)
  {
    /* ifodmro */
    LALI2DestroyVector (status->statusPtr, &segPtr[i].data->data);
    CHECKSTATUSPTR (status);

    LALFree (segPtr[i].data);

    /* power spectrum */
    LALDestroyVector (status->statusPtr, &segPtr[i].spec->data);
    CHECKSTATUSPTR (status);

    LALFree (segPtr[i].spec);

    /* response function */
    LALCDestroyVector (status->statusPtr, &segPtr[i].resp->data);
    CHECKSTATUSPTR (status);

    LALFree (segPtr[i].resp);
  }

  /* free the array */
  LALFree( segPtr );

  /* free the vector structure */
  LALFree( *vector );
  *vector = NULL;


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}



void
LALCreateFindChirpSegmentVector (
    LALStatus                  *status,
    FindChirpSegmentVector    **vector,
    FindChirpInitParams        *params
    )
{
  UINT4                         i;
  FindChirpSegmentVector       *vectorPtr;
  FindChirpSegment             *segPtr;

  INITSTATUS( status, "LALCreateFindChirpSegmentVector", FINDCHIRPMEMORYC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* make sure the output handle exists, but points to a null pointer */
  ASSERT( vector, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );
  ASSERT( !*vector, status, FINDCHIRP_ENNUL, FINDCHIRP_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT( params, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );

  /* make sure that the number of segments, number of points in a */
  /* segment and number of chisquared bins are positive           */
  ASSERT( params->numSegments > 0, 
      status, FINDCHIRP_ESEGZ, FINDCHIRP_MSGESEGZ );
  ASSERT( params->numPoints > 0, 
      status, FINDCHIRP_ENUMZ, FINDCHIRP_MSGENUMZ );
  ASSERT( params->numChisqBins > 0, status, 
      FINDCHIRP_ECHIZ, FINDCHIRP_MSGECHIZ );


  /*
   *
   * create a vector of findchirp segments 
   *
   */


  /* create the output structure */
  vectorPtr = *vector = (FindChirpSegmentVector *) 
    LALMalloc( sizeof(FindChirpSegmentVector) );
  ASSERT( vectorPtr, status, FINDCHIRP_EALOC, FINDCHIRP_MSGEALOC );
  memset( vectorPtr, 0, sizeof(FindChirpSegmentVector) );

  /* set the number of segments in the vector */
  vectorPtr->length = params->numSegments;

  /* allocate memory for an array of findchirp segments */
  segPtr = vectorPtr->data = (FindChirpSegment *) 
    LALMalloc( vectorPtr->length * sizeof(FindChirpSegment) );
  if ( !segPtr )
  {
    LALFree( vectorPtr );
    vectorPtr = NULL;
    ABORT( status, FINDCHIRP_EALOC, FINDCHIRP_MSGEALOC );
  }
  memset( segPtr, 0, vectorPtr->length * sizeof(FindChirpSegment) );

  /* create the individual FindChirpSegments in the vector */
  /* just give up if anything fails                        */
  /* this should be cleaned up to prevent leaks            */
  for (i = 0; i < vectorPtr->length; ++i)
  {
    /* template independent part of filter */
    segPtr[i].data = (COMPLEX8FrequencySeries *)
      LALMalloc( sizeof(COMPLEX8FrequencySeries));
    memset( segPtr[i].data, 0, sizeof(COMPLEX8FrequencySeries) );

    LALCCreateVector (status->statusPtr, 
        &segPtr[i].data->data, params->numPoints/2 + 1);
    CHECKSTATUSPTR (status);

    /* chi squared frequency bins */
    segPtr[i].chisqBinVec = NULL;

    LALU4CreateVector (status->statusPtr, 
        &segPtr[i].chisqBinVec, params->numChisqBins + 1);
    CHECKSTATUSPTR (status);

    /* segment dependent part of normalisation */
    segPtr[i].segNorm = 0.0;

    /* segment id number (invalid) */
    segPtr[i].number = -1;
  }


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}



void
LALDestroyFindChirpSegmentVector (
    LALStatus                  *status,
    FindChirpSegmentVector    **vector
    )
{
  UINT4                         i;
  FindChirpSegment             *segPtr;

  INITSTATUS( status, "LALDestroyFindChirpSegmentVector", FINDCHIRPMEMORYC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */

  
  /* make sure handle is non-null and points to a non-null pointer */
  ASSERT( vector, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );
  ASSERT( *vector, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );

  /* local pointer to the segment array */
  segPtr = (*vector)->data;

  /* destroy the contents of the array */
  for (i = 0; i < (*vector)->length; ++i)
  {
    /* template independent part of stationary phase filter */
    LALCDestroyVector (status->statusPtr, &segPtr[i].data->data);
    CHECKSTATUSPTR (status);

    /* chi squared frequency bins */
    LALU4DestroyVector (status->statusPtr, &segPtr[i].chisqBinVec);
    CHECKSTATUSPTR (status);

    /* frequency series pointer */
    LALFree (segPtr[i].data);
  }

  /* free the array */
  LALFree( segPtr );

  /* free the vector structure */
  LALFree( *vector );
  *vector = NULL;


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
