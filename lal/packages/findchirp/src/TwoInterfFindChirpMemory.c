/*----------------------------------------------------------------------- 
 * 
 * File Name: TwoInterfFindChirpMemory.c
 *
 * Author: Bose, S.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */
#include <lal/LALRCSID.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/DataBuffer.h>
#include <lal/FindChirp.h>
#include <lal/PrintVector.h>
#include <lal/PrintFTSeries.h>
#include <lal/TwoInterfFindChirp.h>

NRCSID (TWOINTERFFINDCHIRPMEMORYC, "$Id$");

void
LALCreateTwoInterfDataSegmentVector (
    LALStatus                           *status,
    TwoInterfDataSegmentVector         **vector,
    TwoInterfFindChirpInitParams        *params
    )
{
  UINT4                          i;
  UINT4                          j;
  TwoInterfDataSegmentVector    *vectorPtr;
  DataSegmentVector             *segVecPtr;
  DataSegment                   *segPtr;
  
  INITSTATUS( status, "LALCreateTwoInterfDataSegmentVector", TWOINTERFFINDCHIRPMEMORYC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */
  

  /* make sure the output handle exists, but points to a null pointer */
  ASSERT( vector, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( !*vector, status, TWOINTERFFINDCHIRPH_ENNUL, TWOINTERFFINDCHIRPH_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT( params, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );

  /* make sure that the number of segments and number of points in a */
  /* segment are positive                                            */
  ASSERT( params->numDetectors > 0, status, 
      TWOINTERFFINDCHIRPH_ESEGZ, TWOINTERFFINDCHIRPH_MSGESEGZ );
  ASSERT( params->numSegments > 0, status, 
      TWOINTERFFINDCHIRPH_ESEGZ, TWOINTERFFINDCHIRPH_MSGESEGZ );
  ASSERT( params->numPoints > 0, status, 
      TWOINTERFFINDCHIRPH_ENUMZ, TWOINTERFFINDCHIRPH_MSGENUMZ );


  /*
   *
   * create a vector of DataSegments
   *
   */
  

  /* create the ouput structure */
  vectorPtr = *vector = (TwoInterfDataSegmentVector *) 
    LALCalloc( 1, sizeof(TwoInterfDataSegmentVector) );
  if ( ! vectorPtr )
  {
    ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC );
  }

  /* set the number of detectors in the vector */
  vectorPtr->length = params->numDetectors;

  /* allocate memory for an array  */
  segVecPtr = vectorPtr->data = (DataSegmentVector *)
    LALCalloc (1, vectorPtr->length*sizeof(DataSegmentVector));
  if ( !segVecPtr )
    {
      LALFree( vectorPtr );
      vectorPtr = NULL;
      ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC);
    }
  /* create the individual DataSegments in the vector */
  /* XXX just give up if anything fails               */
  /* XXX this should be cleaned up to prevent leaks   */
  for (j = 0; j < vectorPtr->length; ++j)
    {
      segVecPtr[j].length = params->numSegments;
      
      /* allocate memory for an array of data segments */
      segPtr = segVecPtr[j].data = (DataSegment *) 
	LALCalloc( 1, segVecPtr[j].length * sizeof(DataSegment) );
      if ( !segPtr )
	{
	  LALFree( segVecPtr );
	  segVecPtr = NULL;
	  ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC );
	}
      
      /* create the individual DataSegments in the vector */
      /* XXX just give up if anything fails               */
      /* XXX this should be cleaned up to prevent leaks   */
      for (i = 0; i < segVecPtr[j].length; ++i)
	{
	  /* channel */
	  segPtr[i].chan = (REAL4TimeSeries *) 
	    LALCalloc( 1, sizeof(REAL4TimeSeries) );
	  if ( ! segPtr[i].chan )
	    {
	      ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC );
	    }
	  
	  LALCreateVector (status->statusPtr, 
			   &segPtr[i].chan->data, params->numPoints);
	  CHECKSTATUSPTR (status);
	  
	  /* power spectrum */
	  segPtr[i].spec = (REAL4FrequencySeries *) 
	    LALCalloc( 1, sizeof(REAL4FrequencySeries) );
	  if ( ! segPtr[i].spec )
	    {
	      ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC );
	    }
	  
	  LALCreateVector (status->statusPtr, 
			   &segPtr[i].spec->data, params->numPoints/2 + 1);
	  CHECKSTATUSPTR (status);
	  
	  /* response function */
	  segPtr[i].resp = (COMPLEX8FrequencySeries *) 
	    LALCalloc( 1, sizeof(COMPLEX8FrequencySeries) );
	  if ( ! segPtr[i].resp )
	    {
	      ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC );
	    }
	  
	  LALCCreateVector (status->statusPtr, 
			    &segPtr[i].resp->data, params->numPoints/2 + 1);
	  CHECKSTATUSPTR (status);
	}
    }
  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

#pragma <lalVerbatim file="TwoInterfFindChirpMemoryCP">
void
LALDestroyTwoInterfDataSegmentVector (
    LALStatus                           *status,
    TwoInterfDataSegmentVector         **vector
    )
#pragma </lalVerbatim>
{
  UINT4                          i;
  UINT4                          j;
  DataSegmentVector             *segVecPtr;
  DataSegment                   *segPtr;
  
  INITSTATUS( status, "LALDestroyTwoInterfDataSegmentVector", TWOINTERFFINDCHIRPMEMORYC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */
  

  ASSERT( vector, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( *vector, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );


  /*
   *
   * destroy the data segment vector
   *
   */
  
  segVecPtr = (*vector)->data;

  /* destroy the contents of the array */
  for (j = 0; j < (*vector)->length; ++j)
    {
      segPtr = segVecPtr[j].data;
 
      /* destroy the contents of the array */
      for (i = 0; i < segVecPtr[j].length; ++i)
	{
	  /* data channel */
	  LALDestroyVector (status->statusPtr, &segPtr[i].chan->data);
	  CHECKSTATUSPTR (status);
	  
	  LALFree (segPtr[i].chan);
	  
	  /* power spectrum */
	  LALDestroyVector (status->statusPtr, &segPtr[i].spec->data);
	  CHECKSTATUSPTR (status);
	  
	  LALFree (segPtr[i].spec);
	  
	  /* response function */
	  LALCDestroyVector (status->statusPtr, &segPtr[i].resp->data);
	  CHECKSTATUSPTR (status);
	  
	  LALFree (segPtr[i].resp);
	}
      LALFree( segPtr );
    }
  /* free the array */
  LALFree( segVecPtr );
  
  /* free the vector structure */
  LALFree( *vector );
  *vector = NULL;
  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


#pragma <lalVerbatim file="TWOINTERFFindChirpMemoryCP">
void
LALCreateTwoInterfFindChirpSegmentVector (
    LALStatus                           *status,
    TwoInterfFindChirpSegmentVector    **vector,
    TwoInterfFindChirpInitParams        *params
    )
#pragma </lalVerbatim>
{
  UINT4                                  i;
  UINT4                                  j;
  TwoInterfFindChirpSegmentVector       *vectorPtr;
  FindChirpSegmentVector                *segVecPtr;
  FindChirpSegment                      *segPtr;

  INITSTATUS( status, "LALCreateTwoInterfFindChirpSegmentVector", TWOINTERFFINDCHIRPMEMORYC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* make sure the output handle exists, but points to a null pointer */
  ASSERT( vector, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( !*vector, status, TWOINTERFFINDCHIRPH_ENNUL, TWOINTERFFINDCHIRPH_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT( params, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );

  /* make sure that the number of segments, number of points in a */
  /* segment and number of chisquared bins are positive           */
  ASSERT( params->numSegments > 0, 
      status, TWOINTERFFINDCHIRPH_ESEGZ, TWOINTERFFINDCHIRPH_MSGESEGZ );
  ASSERT( params->numPoints > 0, 
      status, TWOINTERFFINDCHIRPH_ENUMZ, TWOINTERFFINDCHIRPH_MSGENUMZ );


  /*
   *
   * create a vector of findchirp segments 
   *
   */


  /* create the output structure */
  vectorPtr = *vector = (TwoInterfFindChirpSegmentVector *) 
    LALCalloc( 1, sizeof(TwoInterfFindChirpSegmentVector) );
  if ( ! vectorPtr )
    {
      ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC );
    }
  /* set the number of detectors in the vector */
  vectorPtr->length = params->numDetectors;
  
  /* allocate memory for an array  */
  segVecPtr = vectorPtr->data = (FindChirpSegmentVector *)
    LALCalloc (1, vectorPtr->length*sizeof(FindChirpSegmentVector));
  if ( !segVecPtr )
    {
      LALFree( vectorPtr );
      vectorPtr = NULL;
      ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC);
    }

  for (j = 0; j < vectorPtr->length; ++j)
    {
      segVecPtr[j].length = params->numSegments;
      
      /* allocate memory for an array of data segments */
      segPtr = segVecPtr[j].data = (FindChirpSegment *) 
	LALCalloc( 1, segVecPtr[j].length * sizeof(FindChirpSegment) );
      if ( !segPtr )
	{
	  LALFree( segVecPtr );
	  segVecPtr = NULL;
	  ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC );
	}
      
      /* create the individual FindChirpSegments in the vector */
      /* XXX just give up if anything fails                    */
      /* XXX this should be cleaned up to prevent leaks        */
      for (i = 0; i < segVecPtr[j].length; ++i)
	{
	  /* template independent part of filter */
	  segPtr[i].data = (COMPLEX8FrequencySeries *)
	    LALCalloc( 1, sizeof(COMPLEX8FrequencySeries));
	  if ( ! segPtr[i].data )
	    {
	      ABORT( status, TWOINTERFFINDCHIRPH_EALOC, TWOINTERFFINDCHIRPH_MSGEALOC );
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
		LALCalloc( 1, sizeof(UINT4Vector) );
	    }
	  
	  /* segment dependent part of normalisation */
	  segPtr[i].segNorm = 0.0;
	  
	  /* segment id number (invalid) */
	  segPtr[i].number = -1;
	}
    } 
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

#pragma <lalVerbatim file="TwoInterfFindChirpMemoryCP">
void
LALDestroyTwoInterfFindChirpSegmentVector (
    LALStatus                           *status,
    TwoInterfFindChirpSegmentVector    **vector
    )
#pragma </lalVerbatim>
{
  UINT4                               i;
  UINT4                               j;
  FindChirpSegmentVector             *segVecPtr;
  FindChirpSegment                   *segPtr;

  INITSTATUS( status, "LALDestroyTwoInterfFindChirpSegmentVector", TWOINTERFFINDCHIRPMEMORYC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */

  ASSERT( vector, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );
  ASSERT( *vector, status, TWOINTERFFINDCHIRPH_ENULL, TWOINTERFFINDCHIRPH_MSGENULL );

  /*
   *
   * destroy the data segment vector
   *
   */
  
  segVecPtr = (*vector)->data;
  
  /* destroy the contents of the array */
  for (j = 0; j < (*vector)->length; ++j)
    {
      segPtr = segVecPtr[j].data;
      
      /* destroy the contents of the array */
      for (i = 0; i < segVecPtr[j].length; ++i)
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
      LALFree (segPtr);
    }
  /* free the array */
  LALFree( segVecPtr );
  
  /* free the vector structure */
  LALFree( *vector );
  *vector = NULL;

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
