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

#if 0
<lalVerbatim file="FindChirpMemoryCV">
Author: Brown D. A.
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{FindChirpMemory.c}}
\label{ss:FindChirpMemory.c}

Memory management functions for creating and destroying input data and 
workspace memory for findchirp.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpMemoryCP}
\idx{LALInitializeDataSegmentVector()}
\idx{LALFinalizeDataSegmentVector()}
\idx{LALCreateDataSegmentVector()}
\idx{LALDestroyDataSegmentVector()}
\idx{LALCreateFindChirpSegmentVector()}
\idx{LALDestroyFindChirpSegmentVector()}

\subsubsection*{Description}

The function \texttt{LALInitializeDataSegmentVector()} creates a vector of
\texttt{DataSegment} structures of length and dimension specified in the
\texttt{FindChirpInitParams} structure. The data segments are created using
the input data in \texttt{chan}, \texttt{spec} and \texttt{resp}. No storage
is allocated for the actual data in the data segments, instead the data
segment \texttt{chan}, \texttt{spec} and \texttt{resp} pointers are pointed at 
the input data appropriately. This means that the input \texttt{chan},
\texttt{spec} and \texttt{resp} time and frequency series must persist for the
duration of the filtering process.

The function \texttt{LALFinalizeDataSegmentVector()} frees any memory
allocated by the function \texttt{LALInitializeDataSegmentVector()}.

The function \texttt{LALCreateDataSegmentVector()} creates a vector of
\texttt{DataSegment} structures of length and dimension specified in the
\texttt{FindChirpInitParams} structure.  \texttt{**vector} must point to NULL
on entry and contains a handle to the address of the created
\texttt{DataSegmentVector} on exit.

The function \texttt{LALDestroyDataSegmentVector()} frees the memory of the
\texttt{DataSegmentVector} at address \texttt{*vector}.

The function \texttt{LALCreateFindChirpSegmentVector()} creates a vector of
\texttt{FindChirpSegment} structures of length and dimension specified in the
\texttt{FindChirpInitParams} structure.  \texttt{**vector} must point to NULL
on entry and contains a handle to the address of the created
\texttt{DataSegmentVector} on exit.

The function \texttt{LALDestroyFindChirpSegmentVector()} frees the memory of
the \texttt{FindChirpSegmentVector} at address \texttt{*vector}.

\subsubsection*{Algorithm}

None.

\subsubsection*{Uses}
\begin{verbatim}
LALCalloc()
LALCreateVector()
LALCCreateVector()
LALU4CreateVector()
LALFree()
LALDestroyVector()
LALCDestroyVector()
LALU4DestroyVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FindChirpMemoryCV}}
</lalLaTeX>
#endif

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/DataBuffer.h>
#include <lal/FindChirp.h>

NRCSID (FINDCHIRPMEMORYC, "$Id$");

/* <lalVerbatim file="FindChirpMemoryCP"> */
void
LALInitializeDataSegmentVector (
    LALStatus                  *status,
    DataSegmentVector         **dataSegVecPtr,
    REAL4TimeSeries            *chan,
    REAL4FrequencySeries       *spec,
    COMPLEX8FrequencySeries    *resp,
    FindChirpInitParams        *params
    )
/* </lalVerbatim> */
{
  INT8  chanStartTime;
  UINT4 i;
  UINT4 inputLength;
  REAL4 *dataPtr;
  DataSegmentVector *dataSegVec = NULL;

  INITSTATUS( status, "LALPopulateDataSegmentVector", FINDCHIRPMEMORYC );
  ATTATCHSTATUSPTR( status );

  ASSERT( dataSegVecPtr, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( ! *dataSegVecPtr, status, FINDCHIRPH_ENNUL, FINDCHIRPH_MSGENNUL );
  ASSERT( chan, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( chan->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( chan->data->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( spec, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( resp, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );


  /* check that there is enough data to populate the data seg vec */
  inputLength = params->numPoints * params->numSegments - 
    ( params->numSegments - 1 ) * params->ovrlap;
  if ( inputLength > chan->data->length )
  {
    ABORT( status, FINDCHIRPH_ESMSM, FINDCHIRPH_MSGESMSM );
  }

  /* allocate memory for the data segment vector */
  if ( ! (dataSegVec = *dataSegVecPtr = (DataSegmentVector *) 
        LALCalloc( 1, sizeof(DataSegmentVector) )) )
  {
    ABORT (status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC);
  }
  dataSegVec->length = params->numSegments;

  if ( ! ( dataSegVec->data =  (DataSegment *) 
        LALCalloc( dataSegVec->length, sizeof(DataSegment) )) )
  {
    ABORT (status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC);
  }

  /* allocate the DataSegments in the DataSegmentVector */
  for ( i = 0; i < dataSegVec->length; ++i )
  {
    /* point to current segment */
    DataSegment *currentSegment = dataSegVec->data + i;

    /* channel */
    if ( ! (currentSegment->chan = (REAL4TimeSeries *) 
          LALCalloc( 1, sizeof(REAL4TimeSeries) )) )
    {
      ABORT (status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC);
    }

    if ( ! (currentSegment->chan->data = (REAL4Sequence *) 
          LALCalloc( 1, sizeof(REAL4Sequence) )) )
    {
      ABORT (status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC);
    }
  }

  /* store the start of the input channel and its gps start time */
  dataPtr = chan->data->data;
  LALGPStoINT8( status->statusPtr, &chanStartTime, &(chan->epoch) );
  CHECKSTATUSPTR( status );

  for ( i = 0; i < dataSegVec->length; ++i )
  {
    /* point to current segment */
    DataSegment      *currentSegment = dataSegVec->data + i;
    REAL4Sequence    *origDataPtr;

    /* this should be set to a unique number for each segment   */
    currentSegment->number = i;

    /* copy the REAL4TimeSeries information into the data segment */
    origDataPtr = currentSegment->chan->data;
    memcpy( currentSegment->chan, chan, sizeof(REAL4TimeSeries) );
    currentSegment->chan->data = origDataPtr;

    /* set the correct GPS time for the current data segment */
    {
      INT8 currentSegmentTime = chanStartTime + (INT8) 
        ( 1.0e9 * 
          ((REAL8) params->numPoints - (REAL8) params->ovrlap) *
          (REAL8) i * chan->deltaT
        );

      LALINT8toGPS( status->statusPtr, &(currentSegment->chan->epoch), 
          &currentSegmentTime );
      CHECKSTATUSPTR( status );
    }

    /* set up the REAL4Sequence to contain the correct data */
    currentSegment->chan->data->length = params->numPoints;
    currentSegment->chan->data->data = dataPtr;

    /* advance the dataPtr to the next segment */
    dataPtr += currentSegment->chan->data->length - params->ovrlap;

    /* power spectrum */
    currentSegment->spec = spec;

    /* response function */
    currentSegment->resp = resp;
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="FindChirpMemoryCP"> */
void
LALFinalizeDataSegmentVector (
    LALStatus                  *status,
    DataSegmentVector         **vector
    )
/* </lalVerbatim> */
{
  DataSegmentVector *dataSegVec;
  UINT4 i;

  INITSTATUS( status, "LALFinalizeDataSegmentVector", FINDCHIRPMEMORYC );

  ASSERT( vector, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( *vector, status, FINDCHIRPH_ENNUL, FINDCHIRPH_MSGENNUL );

  dataSegVec = *vector;

  /* free the DataSegments in the DataSegmentVector */
  for ( i = 0; i < dataSegVec->length; ++i )
  {
    /* point to current segment */
    DataSegment *currentSegment = dataSegVec->data + i;
    LALFree( currentSegment->chan->data );
    LALFree( currentSegment->chan );
  }

  /* free the dataSegVec to hold the input data for FindChirpData */
  LALFree( dataSegVec->data );
  LALFree( dataSegVec );
  dataSegVec = NULL;

  RETURN( status );
}


/* <lalVerbatim file="FindChirpMemoryCP"> */
void
LALCreateDataSegmentVector (
    LALStatus                  *status,
    DataSegmentVector         **vector,
    FindChirpInitParams        *params
    )
/* </lalVerbatim> */
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
    /* channel */
    segPtr[i].chan = (REAL4TimeSeries *) 
      LALCalloc( 1, sizeof(REAL4TimeSeries) );
    if ( ! segPtr[i].chan )
    {
      ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
    }

    LALCreateVector (status->statusPtr, 
        &segPtr[i].chan->data, params->numPoints);
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


/* <lalVerbatim file="FindChirpMemoryCP"> */
void
LALDestroyDataSegmentVector (
    LALStatus                  *status,
    DataSegmentVector         **vector
    )
/* </lalVerbatim> */
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

  /* free the array */
  LALFree( segPtr );

  /* free the vector structure */
  LALFree( *vector );
  *vector = NULL;


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="FindChirpMemoryCP"> */
void
LALCreateFindChirpSegmentVector (
    LALStatus                  *status,
    FindChirpSegmentVector    **vector,
    FindChirpInitParams        *params
    )
/* </lalVerbatim> */
{
  UINT4                         i/*,k*/;
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

  /* check that the approximant is of a known type */
  switch ( params->approximant )
  {
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case TaylorF2:
    case GeneratePPN:
    case PadeT1:
    case EOB:
    case FindChirpSP:
    case BCV:
    case BCVSpin:
      break;
    default:
      ABORT( status, FINDCHIRPH_EUAPX, FINDCHIRPH_MSGEUAPX );
      break;
  }


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

    if ( params->approximant == BCV )
    {
      segPtr[i].dataBCV = (COMPLEX8FrequencySeries *)
        LALCalloc( 1, sizeof(COMPLEX8FrequencySeries));
      if ( ! segPtr[i].dataBCV )
      {
        ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
      }

      LALCCreateVector (status->statusPtr,
          &segPtr[i].dataBCV->data, params->numPoints/2 + 1);
      CHECKSTATUSPTR (status);
    }

    /* segment dependent part of the normalization, for the BCV templates */
    if ( params->approximant == BCV )
    {
      LALCreateVector (status->statusPtr,
          &(segPtr[i].a1), params->numPoints/2 + 1);
      CHECKSTATUSPTR (status);

      LALCreateVector (status->statusPtr,
          &(segPtr[i].b1), params->numPoints/2 + 1);
      CHECKSTATUSPTR (status);

      LALCreateVector (status->statusPtr,
          &(segPtr[i].b2), params->numPoints/2 + 1);
      CHECKSTATUSPTR (status);

#if 0
      LALCreateVector (status->statusPtr,
          &(segPtr[i].tmpltPowerVec), params->numPoints/2 + 1);
      CHECKSTATUSPTR (status);
#endif

      LALCreateVector (status->statusPtr,
          &(segPtr[i].tmpltPowerVecBCV),params->numPoints/2 + 1);  
      CHECKSTATUSPTR (status);
    }

    /* allocate space for the chi squared frequency bin boundaries */
    segPtr[i].chisqBinVec = (UINT4Vector *) 
      LALCalloc( 1, sizeof(UINT4Vector) );
    if ( ! segPtr[i].chisqBinVec )
    {
      ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
    }
    segPtr[i].chisqBinVec->length = 
      params->numChisqBins ? params->numChisqBins + 1 : 0;

    /* additional chisq frequency bins for the BCV templates */
    if ( params->approximant == BCV )
    { 
      if ( params->numChisqBins )
      {
        segPtr[i].chisqBinVecBCV = NULL;
        LALU4CreateVector (status->statusPtr, 
            &segPtr[i].chisqBinVecBCV, params->numChisqBins + 1);
        CHECKSTATUSPTR (status);
      }
      else
      {
        segPtr[i].chisqBinVecBCV = (UINT4Vector *) 
          LALCalloc( 1, sizeof(UINT4Vector) );
      }
    }

    /* segment dependent part of normalisation, for the SP templates */
    LALCreateVector( status->statusPtr, &(segPtr[i].segNorm), 
        params->numPoints/2 + 1 );
    CHECKSTATUSPTR( status );

    /* segment id number (invalid) */
    segPtr[i].number = -1;
  }


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="FindChirpMemoryCP"> */
void
LALDestroyFindChirpSegmentVector (
    LALStatus                  *status,
    FindChirpSegmentVector    **vector
    )
/* </lalVerbatim> */
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
    if ( segPtr[i].dataBCV && segPtr[i].dataBCV->data )
    {
      LALCDestroyVector (status->statusPtr, &segPtr[i].dataBCV->data);
      CHECKSTATUSPTR (status);
    }

    /* chi squared frequency bins */
    if ( segPtr[i].chisqBinVec->data )
    {
      LALFree( segPtr[i].chisqBinVec->data );
      segPtr[i].chisqBinVec->data = NULL;
    }
    LALFree( segPtr[i].chisqBinVec );

    if ( segPtr[i].chisqBinVecBCV )
    {
      if ( segPtr[i].chisqBinVecBCV->length )
      {
        LALU4DestroyVector (status->statusPtr, &segPtr[i].chisqBinVecBCV);
        CHECKSTATUSPTR (status);
      }
      else
      {
        LALFree( segPtr[i].chisqBinVecBCV );
      }
    }

    /* frequency series pointer */
    LALFree (segPtr[i].data);
    if ( segPtr[i].dataBCV )
    {
      LALFree( segPtr[i].dataBCV );
    }

    if ( segPtr[i].segNorm )
    {
      LALDestroyVector( status->statusPtr, &(segPtr[i].segNorm) );
      CHECKSTATUSPTR( status );
    }

    if ( segPtr[i].a1 )
    {
      LALDestroyVector( status->statusPtr, &(segPtr[i].a1) );
      CHECKSTATUSPTR( status );
    }
    if ( segPtr[i].b1 )
    {
      LALDestroyVector( status->statusPtr, &(segPtr[i].b1) );
      CHECKSTATUSPTR( status );
    }
    if ( segPtr[i].b2 )
    {
      LALDestroyVector( status->statusPtr, &(segPtr[i].b2) );
      CHECKSTATUSPTR( status );
    }
#if 0
    if ( segPtr[i].tmpltPowerVec )
    {
      LALDestroyVector( status->statusPtr, &(segPtr[i].tmpltPowerVec) );
      CHECKSTATUSPTR( status );
    }
#endif
    if ( segPtr[i].tmpltPowerVecBCV )
    {
      LALDestroyVector( status->statusPtr, &(segPtr[i].tmpltPowerVecBCV) );
      CHECKSTATUSPTR( status );
    }
  } /* end "for (i = 0; i < (*vector)->length; ++i)" */

  /* free the array */
  LALFree( segPtr );

  /* free the vector structure */
  LALFree( *vector );
  *vector = NULL;


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
