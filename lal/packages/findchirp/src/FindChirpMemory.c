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
  ASSERT( vector, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( !*vector, status, FINDCHIRPH_ENNUL, FINDCHIRPH_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure that the number of segments and number of points in a */
  /* segment are positive                                            */
  ASSERT( params->numSegments > 0, status, 
      FINDCHIRPH_ESEGZ, FINDCHIRPH_MSGESEGZ );
  ASSERT( params->numPoints > 0, status, 
      FINDCHIRPH_ENUMZ, FINDCHIRPH_MSGENUMZ );


  /*
   *
   * create a vector of DataSegments
   *
   */
  

  /* create the ouput structure */
  vectorPtr = *vector = (DataSegmentVector *) 
    LALCalloc( 1, sizeof(DataSegmentVector) );
  if ( ! vectorPtr )
  {
    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
  }

  /* set the number of segments in the vector */
  vectorPtr->length = params->numSegments;

  /* allocate memory for an array of data segments */
  segPtr = vectorPtr->data = (DataSegment *) 
    LALCalloc( 1, vectorPtr->length * sizeof(DataSegment) );
  if ( !segPtr )
  {
    LALFree( vectorPtr );
    vectorPtr = NULL;
    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
  }

  /* create the individual DataSegments in the vector */
  /* XXX just give up if anything fails               */
  /* XXX this should be cleaned up to prevent leaks   */
  for (i = 0; i < vectorPtr->length; ++i)
  {
    /* ifodmro */
    segPtr[i].data = (INT2TimeSeries *) 
      LALCalloc( 1, sizeof(INT2TimeSeries) );
    if ( ! segPtr[i].data )
    {
      ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
    }

    LALI2CreateVector (status->statusPtr, 
        &segPtr[i].data->data, params->numPoints);
    CHECKSTATUSPTR (status);

    /* as_q */
    segPtr[i].real4Data = (REAL4TimeSeries *) 
      LALCalloc( 1, sizeof(REAL4TimeSeries) );
    if ( ! segPtr[i].real4Data )
    {
      ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
    }

    LALCreateVector (status->statusPtr, 
        &segPtr[i].real4Data->data, params->numPoints);
    CHECKSTATUSPTR (status);

    /* power spectrum */
    segPtr[i].spec = (REAL4FrequencySeries *) 
      LALCalloc( 1, sizeof(REAL4FrequencySeries) );
    if ( ! segPtr[i].spec )
    {
      ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
    }

    LALCreateVector (status->statusPtr, 
        &segPtr[i].spec->data, params->numPoints/2 + 1);
    CHECKSTATUSPTR (status);
    
    /* response function */
    segPtr[i].resp = (COMPLEX8FrequencySeries *) 
      LALCalloc( 1, sizeof(COMPLEX8FrequencySeries) );
    if ( ! segPtr[i].resp )
    {
      ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
    }

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
  ASSERT( vector, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( *vector, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );


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

    /* as_q */
    LALDestroyVector (status->statusPtr, &segPtr[i].real4Data->data);
    CHECKSTATUSPTR (status);

    LALFree (segPtr[i].real4Data);

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
  ASSERT( vector, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( !*vector, status, FINDCHIRPH_ENNUL, FINDCHIRPH_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure that the number of segments, number of points in a */
  /* segment and number of chisquared bins are positive           */
  ASSERT( params->numSegments > 0, 
      status, FINDCHIRPH_ESEGZ, FINDCHIRPH_MSGESEGZ );
  ASSERT( params->numPoints > 0, 
      status, FINDCHIRPH_ENUMZ, FINDCHIRPH_MSGENUMZ );
  ASSERT( params->numChisqBins >= 0, status, 
      FINDCHIRPH_ECHIZ, FINDCHIRPH_MSGECHIZ );


  /*
   *
   * create a vector of findchirp segments 
   *
   */


  /* create the output structure */
  vectorPtr = *vector = (FindChirpSegmentVector *) 
    LALCalloc( 1, sizeof(FindChirpSegmentVector) );
  if ( ! vectorPtr )
  {
    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
  }

  /* set the number of segments in the vector */
  vectorPtr->length = params->numSegments;

  /* allocate memory for an array of findchirp segments */
  segPtr = vectorPtr->data = (FindChirpSegment *) 
    LALCalloc( 1, vectorPtr->length * sizeof(FindChirpSegment) );
  if ( !segPtr )
  {
    LALFree( vectorPtr );
    vectorPtr = NULL;
    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
  }

  /* create the individual FindChirpSegments in the vector */
  /* XXX just give up if anything fails                    */
  /* XXX this should be cleaned up to prevent leaks        */
  for (i = 0; i < vectorPtr->length; ++i)
  {
    /* template independent part of filter */
    segPtr[i].data = (COMPLEX8FrequencySeries *)
      LALCalloc( 1, sizeof(COMPLEX8FrequencySeries));
    if ( ! segPtr[i].data )
    {
      ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
    }

    LALCCreateVector (status->statusPtr, 
        &segPtr[i].data->data, params->numPoints/2 + 1);
    CHECKSTATUSPTR (status);

    /* chi squared frequency bins */
    if ( params->numChisqBins )
    {
      segPtr[i].chisqBinVec = NULL;

      LALU4CreateVector (status->statusPtr, 
          &segPtr[i].chisqBinVec, params->numChisqBins + 1);
      CHECKSTATUSPTR (status);
    }
    else
    {
      segPtr[i].chisqBinVec = (UINT4Vector *) 
        LALMalloc( sizeof(UINT4Vector) );

      segPtr[i].chisqBinVec->length = 0;
      segPtr[i].chisqBinVec->data = NULL;
    }

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
  ASSERT( vector, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( *vector, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* local pointer to the segment array */
  segPtr = (*vector)->data;

  /* destroy the contents of the array */
  for (i = 0; i < (*vector)->length; ++i)
  {
    /* template independent part of stationary phase filter */
    LALCDestroyVector (status->statusPtr, &segPtr[i].data->data);
    CHECKSTATUSPTR (status);

    /* chi squared frequency bins */
    if ( segPtr[i].chisqBinVec->length )
    {
      LALU4DestroyVector (status->statusPtr, &segPtr[i].chisqBinVec);
      CHECKSTATUSPTR (status);
    }
    else
    {
      LALFree( segPtr[i].chisqBinVec );
    }

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
