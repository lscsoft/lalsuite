/**** <lalVerbatim file="StringSearchConditionDataCV">
 * Author: Jocelyn Read
 * $Id$
 **** </lalVerbatim> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Units.h>
#include <lal/TimeFreqFFT.h>
#include <lal/VectorOps.h>
#include <lal/StringSearch.h>
#include <lal/String.h>

/**** <lalLaTeX>
 *
 * \subsection{Module \texttt{RingSearchConditionData.c}}
 *
 * Routine to condition data prior to a ring search.
 *
 * \subsubsection*{Prototypes}
 * \input{RingSearchConditionDataCP}
 * \idx{LALRingSearchConditionData()}
 * 
 * \subsubsection*{Description}
 *
 * The function \verb+LALRingSearchConditionData()+ takes an input channel
 * data, power spectrum data, and response data, and computes the strain data
 * segments (in frequency domain) and the (truncated) inverse strain noise
 * spectrum.  These are stored in the \verb+RingSearchParams+ structure.
 *
 * \vfill{\footnotesize\input{RingSearchConditionDataCV}}
 * 
 **** </lalLaTeX> */ 

NRCSID( STRINGSEARCHCONDITIONDATAC, "$Id$" );

/*
static const CHAR ifoNames[][3] =
    { "H0", "H1", "H2", "L0", "L1", "P0", "P1", "P2" };
    */

/* <lalVerbatim file="StringSearchConditionDataCP"> */
void
LALStringSearchConditionData(
    LALStatus               *status,
    StringSearchParams        *params,
    StringSearchData          *data
    )
{ /* </lalVerbatim> */
  RAT4 minusOne = { -1, 0 };
  LALUnitPair unitPair;
  LALUnit     unit;
  REAL4 norm = 1;
  UINT4 cut;
  UINT4 i;
  UINT4 k;

  INITSTATUS( status, "StringSearchConditionData", STRINGSEARCHCONDITIONDATAC );
  ATTATCHSTATUSPTR( status );

  ASSERT( data, status, STRINGSEARCHH_ENULL, STRINGSEARCHH_MSGENULL );
  ASSERT( params, status, STRINGSEARCHH_ENULL, STRINGSEARCHH_MSGENULL );

  params->sampleRate = 1 / data->channel->deltaT;


  /*
   *
   * Scale response function by dynamic range factor.
   *
   */

  if ( data->response->data->length != params->segmentSize / 2 + 1 )
  {
    ABORT( status, STRINGSEARCHH_ESZMM, STRINGSEARCHH_MSGESZMM );
  }
  for ( k = 0; k < data->response->data->length; ++k )
  {
    data->response->data->data[k].re *= params->dynRangeFac;
    data->response->data->data[k].im *= params->dynRangeFac;
  }


  /*
   *
   * Segment the data and compute the power spectrum.
   *
   */

  if ( ( 2 * data->channel->data->length ) % params->segmentSize )
  {
    ABORT( status, STRINGSEARCHH_ENSEG, STRINGSEARCHH_MSGENSEG );
  }
  params->numSegments = 2 * data->channel->data->length / params->segmentSize
    - 1;

  if ( data->spectrum ) /* copy given spectrum to inverse spectrum */
  {
    memcpy( params->invSpectrum->data->data, data->spectrum->data->data,
        params->invSpectrum->data->length
        * sizeof( *params->invSpectrum->data->data ) );
    params->invSpectrum->sampleUnits = data->spectrum->sampleUnits;
    params->invSpectrum->epoch       = data->spectrum->epoch;
    params->invSpectrum->deltaF      = data->spectrum->deltaF;
    params->invSpectrum->f0          = data->spectrum->f0;
  }
  else
  {
    params->invSpectrum->epoch  = data->channel->epoch;
    params->invSpectrum->deltaF = params->sampleRate / params->segmentSize;
    params->invSpectrum->f0     = 0;
    memset( params->invSpectrum->data->data, 0,
        params->invSpectrum->data->length
        * sizeof( *params->invSpectrum->data->data ) );
  }
  strncpy( params->invSpectrum->name, "inverse spectrum",
      sizeof( params->invSpectrum->name ) );

  params->dataSegment = LALCalloc( params->numSegments,
      sizeof( *params->dataSegment ) );
  if ( ! params->dataSegment )
  {
    ABORT( status, STRINGSEARCHH_EALOC, STRINGSEARCHH_MSGEALOC );
  }

  for ( i = 0; i < params->numSegments; ++i )
  {
    REAL4TimeSeries ser = *data->channel;
    REAL4Vector     vec;
    REAL8 dt;
    INT8 timeNS;

    /* create a dummy time series for this segment of data */
    vec.length = params->segmentSize;
    vec.data   = data->channel->data->data + i * params->segmentSize / 2;
    ser.data   = &vec;
    dt = 0.5 * i * params->segmentSize * ser.deltaT;
    timeNS  = (INT8)( 1000000000 ) * (INT8)( ser.epoch.gpsSeconds );
    timeNS += (INT8)( ser.epoch.gpsNanoSeconds );
    timeNS += (INT8)( 1e9 * dt );
    ser.epoch.gpsSeconds     = timeNS / (INT8)( 1000000000 );
    ser.epoch.gpsNanoSeconds = timeNS % (INT8)( 1000000000 );

    /* create this segment's frequency series */
    LALCCreateVector( status->statusPtr, &params->dataSegment[i].data,
        params->segmentSize / 2 + 1 );
    CHECKSTATUSPTR( status );

    /* TODO: NAME */

    /* fft the dummy time series */
    LALTimeFreqRealFFT( status->statusPtr, params->dataSegment + i, &ser,
        params->forwardPlan );
    CHECKSTATUSPTR( status );


    if ( ! data->spectrum ) /* compute the spectrum from the data */
    {
      for ( k = 0; k < params->invSpectrum->data->length; ++k )
      {
        REAL4 re = params->dataSegment[i].data->data[k].re;
        REAL4 im = params->dataSegment[i].data->data[k].im;
        params->invSpectrum->data->data[k] += re * re + im * im;
      }
    }

    /* multiply by the response function */
    LALCCVectorMultiply( status->statusPtr, params->dataSegment[i].data,
        params->dataSegment[i].data, data->response->data );
    CHECKSTATUSPTR( status );

    /* get the units right */
    unitPair.unitOne = &params->dataSegment[i].sampleUnits;
    unitPair.unitTwo = &data->response->sampleUnits;
    LALUnitMultiply( status->statusPtr, &params->dataSegment[i].sampleUnits,
        &unitPair );
    CHECKSTATUSPTR( status );
  }

  /* correct power spectrum normalization prior to inversion */
  if ( ! data->spectrum )
  {
    for ( k = 1; k < params->invSpectrum->data->length; ++k )
    {
      params->invSpectrum->data->data[k] *= 2 * params->invSpectrum->deltaF;
      params->invSpectrum->data->data[k] /= params->numSegments;
    }
    /* power spectrum units prior to inversion */
    unitPair.unitOne = &data->channel->sampleUnits;
    unitPair.unitTwo = &data->channel->sampleUnits;
    LALUnitMultiply( status->statusPtr, &params->invSpectrum->sampleUnits,
        &unitPair );
    CHECKSTATUSPTR( status );
    unitPair.unitOne = &params->invSpectrum->sampleUnits;
    unitPair.unitTwo = &lalSecondUnit;
    LALUnitMultiply( status->statusPtr, &params->invSpectrum->sampleUnits,
        &unitPair );
    CHECKSTATUSPTR( status );
  }



  /*
   *
   * Compute inverse power spectrum perhaps doing truncation.
   *
   */

  cut = params->lowFrequency / params->invSpectrum->deltaF;

  if ( params->invSpecTrunc ) /* truncate inverse spectrum */
  {
    COMPLEX8Vector *vtilde = NULL;
    REAL4Vector    *vector = NULL;

    LALSCreateVector( status->statusPtr, &vector, params->segmentSize );
    CHECKSTATUSPTR( status );
    LALCCreateVector( status->statusPtr, &vtilde, params->segmentSize / 2 + 1 );
    CHECKSTATUSPTR( status );

    memset( vtilde->data, 0, cut * sizeof( *vtilde->data ) );
    for ( k = cut; k < vtilde->length - 1; ++k )
    {
      vtilde->data[k].re = 1 / sqrt( params->invSpectrum->data->data[k] );
      vtilde->data[k].im = 0;
    }
    vtilde->data[vtilde->length - 1].re = 0;
    vtilde->data[vtilde->length - 1].im = 0;

    LALReverseRealFFT( status->statusPtr, vector, vtilde, params->reversePlan );
    CHECKSTATUSPTR( status );

    memset( vector->data + params->invSpecTrunc / 2, 0,
        ( vector->length - params->invSpecTrunc ) * sizeof( *vector->data ) );

    LALRealPowerSpectrum( status->statusPtr, params->invSpectrum->data, vector,
        params->forwardPlan );
    CHECKSTATUSPTR( status );
    /* CHANGE BACK TO ORIGINAL NORMALIZATION -- JC */
    {
      REAL4Vector *myvector = params->invSpectrum->data;
      UINT4 mybin;
      for ( mybin = 1; mybin < myvector->length - 1; ++mybin )
        myvector->data[mybin] *= 0.5;
    }

    /* adjust normalization to account for reverse/forward ffting */
    norm /= (REAL4)( vector->length ) * (REAL4)( vector->length );

    LALCDestroyVector( status->statusPtr, &vtilde );
    CHECKSTATUSPTR( status );
    LALSDestroyVector( status->statusPtr, &vector );
    CHECKSTATUSPTR( status );
  }
  else /* otherwise just invert */
  {
    for ( k = cut; k < params->invSpectrum->data->length - 1; ++k )
    {
      params->invSpectrum->data->data[k] =
        1 / params->invSpectrum->data->data[k];
    }
  }

  /* invert spectrum units */
  LALUnitRaise( status->statusPtr, &params->invSpectrum->sampleUnits,
      &params->invSpectrum->sampleUnits, &minusOne );
  CHECKSTATUSPTR( status );
  
  /* multiply by response function */
  memset( params->invSpectrum->data->data, 0, cut *
      sizeof( *params->invSpectrum->data->data ) );
  for ( k = cut; k < params->invSpectrum->data->length - 1; ++k )
  {
    REAL4 re = data->response->data->data[k].re;
    REAL4 im = data->response->data->data[k].im;
    params->invSpectrum->data->data[k] *= norm / ( re * re + im * im );
  }
  params->invSpectrum->data->data[params->invSpectrum->data->length - 1] = 0;
  LALUnitRaise( status->statusPtr, &unit, &data->response->sampleUnits,
      &minusOne );
  CHECKSTATUSPTR( status );
  unitPair.unitOne = &unit;
  unitPair.unitTwo = &unit;
  LALUnitMultiply( status->statusPtr, &unit, &unitPair );
  CHECKSTATUSPTR( status );
  unitPair.unitOne = &params->invSpectrum->sampleUnits;
  unitPair.unitTwo = &unit;
  LALUnitMultiply( status->statusPtr, &params->invSpectrum->sampleUnits,
      &unitPair );
  CHECKSTATUSPTR( status );


  /*
   *
   * Multiply data segments by inverse spectrum.
   *
   */

  for ( i = 0; i < params->numSegments; ++i )
  {
    LALSCVectorMultiply( status->statusPtr, params->dataSegment[i].data,
        params->invSpectrum->data, params->dataSegment[i].data );
    CHECKSTATUSPTR( status );
    unitPair.unitOne = &params->dataSegment[i].sampleUnits;
    unitPair.unitTwo = &params->invSpectrum->sampleUnits;
    LALUnitMultiply( status->statusPtr, &params->dataSegment[i].sampleUnits,
        &unitPair );
    CHECKSTATUSPTR( status );
  }


  DETATCHSTATUSPTR( status );
  RETURN( status );
}
