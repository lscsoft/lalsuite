/**** <lalVerbatim file="RingSearchConditionDataCV">
 * Author: Jolien Creighton
 * $Id$
 **** </lalVerbatim> */

#include <math.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/Units.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/VectorOps.h>
#include <lal/RingSearch.h>

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

NRCSID( RINGSEARCHCONDITIONDATAC, "$Id$" );

/* static const CHAR ifoNames[][3] =
    { "H0", "H1", "H2", "L0", "L1", "P0", "P1", "P2" }; */

/* <lalVerbatim file="RingSearchConditionDataCP"> */
void
LALRingSearchConditionData(
    LALStatus               *status,
    RingSearchParams        *params,
    RingSearchData          *data
    )
{ /* </lalVerbatim> */
  RAT4 minusOne = { -1, 0 };
  LALUnitPair unitPair;
  LALUnit     unit;
  REAL4 norm = 1;
  UINT4 cut;
  UINT4 i;
  UINT4 k;

  INITSTATUS( status, "LALRingSearchConditionData", RINGSEARCHCONDITIONDATAC );
  ATTATCHSTATUSPTR( status );

  ASSERT( data, status, RINGSEARCHH_ENULL, RINGSEARCHH_MSGENULL );
  ASSERT( params, status, RINGSEARCHH_ENULL, RINGSEARCHH_MSGENULL );

  params->sampleRate = 1 / data->channel->deltaT;

  /* check that the amount of data is correct */
  if ( ( 2 * data->channel->data->length ) % params->segmentSize )
  {
    ABORT( status, RINGSEARCHH_ENSEG, RINGSEARCHH_MSGENSEG );
  }
  params->numSegments = 2 * data->channel->data->length / params->segmentSize - 1;


  /*
   *
   * High-pass filter the data.
   *
   */

  if ( params->highpassFrequency > params->lowFrequency )
  {
    LALWarning( status,
        "highpass frequency should be less than low frequency" );
  }
  if ( params->highpassFrequency == 0 )
  {
    LALInfo( status, "highpass frequency set to low frequency\n"
        "	(use a negative highpass frequency to disable)" );
    params->highpassFrequency = params->lowFrequency;
  }
  if ( params->highpassFrequency > 0 )
  {
    PassBandParamStruc highpassParams;
    highpassParams.nMax =  8;
    highpassParams.f1   = -1;
    highpassParams.a1   = -1;
    highpassParams.f2   = params->highpassFrequency;
    highpassParams.a2   = 0.9; /* this means 10% attenuation at f2 */
    LALDButterworthREAL4TimeSeries( status->statusPtr, data->channel,
        &highpassParams );
    CHECKSTATUSPTR( status );
  }


  /*
   *
   * Scale response function by dynamic range factor.
   *
   */

  if ( data->response->data->length != params->segmentSize / 2 + 1 )
  {
    ABORT( status, RINGSEARCHH_ESZMM, RINGSEARCHH_MSGESZMM );
  }
  for ( k = 0; k < data->response->data->length; ++k )
  {
    data->response->data->data[k].re *= params->dynRangeFac;
    data->response->data->data[k].im *= params->dynRangeFac;
  }


  /*
   *
   * Compute the average power spectrum.  Store it in invSpectrum for now.
   *
   */

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
  else /* compute the average spectrum */
  {
    LALWindowParams       windowParams;
    AverageSpectrumParams avgSpecParams;
    windowParams.type     = Hann;
    windowParams.length   = params->segmentSize;
    avgSpecParams.window  = NULL;
    avgSpecParams.plan    = params->forwardPlan;
    avgSpecParams.method  = params->avgSpecMeth;
    avgSpecParams.overlap = params->segmentSize / 2;
    LALCreateREAL4Window( status->statusPtr, &avgSpecParams.window,
        &windowParams );
    CHECKSTATUSPTR( status );
    LALREAL4AverageSpectrum( status->statusPtr, params->invSpectrum,
        data->channel, &avgSpecParams );
    BEGINFAIL( status )
    {
      LALDestroyREAL4Window( status->statusPtr, &avgSpecParams.window );
    }
    ENDFAIL( status );
    LALDestroyREAL4Window( status->statusPtr, &avgSpecParams.window );
    CHECKSTATUSPTR( status );
    /* correct spectrum normalization for unit spectra */
    if ( params->avgSpecMeth == useUnity && params->avgSpecNorm > 0 )
    {
      for ( k = 0; k < params->invSpectrum->data->length; ++k )
        params->invSpectrum->data->data[k] *= params->avgSpecNorm;
    }
  }


  /*
   *
   * Compute inverse power spectrum perhaps doing truncation.
   *
   */

  strncpy( params->invSpectrum->name, "inverse spectrum",
      sizeof( params->invSpectrum->name ) );

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
   * Code for internal testing purposes: zero the data and/or inject
   * a simple ringdown signal.
   *
   */

  if ( params->testZeroData ) /* set data to zero for analysis */
  {
    memset( data->channel->data->data, 0,
        data->channel->data->length * sizeof( *data->channel->data->data ) );
  }
  if ( params->testInject ) /* very simple internal injection routine */
  {
    INT8 tstart;        /* start time of the data (ns) */
    INT8 tinject;       /* start time of the ring (ns) */
    REAL8 dt;           /* time of ring after start of data (s) */
    
    tstart   = (INT8)1000000000 * (INT8)data->channel->epoch.gpsSeconds;
    tstart  += (INT8)data->channel->epoch.gpsNanoSeconds;
    tinject  = (INT8)1000000000 * (INT8)params->testInjectTime.gpsSeconds;
    tinject += (INT8)params->testInjectTime.gpsNanoSeconds;
    dt       = 1e-9 * (REAL8)( tinject - tstart );
    for ( i = 0; i < data->channel->data->length; ++i )
    {
      REAL8 t; /* time since ringdown start time of this data sample */
      t = i * data->channel->deltaT - dt;
      if ( t >= 0 ) /* ringdown has started: do injection */
      {
        REAL4 phase;
        REAL4 exponent;
        phase = LAL_TWOPI * params->testInjectFreq * t;
        exponent = 0.5 * phase / params->testInjectQual;
        if ( exponent < 60.0 )
        {
          REAL4 val;
          val  = exp( - exponent) * cos( phase + params->testInjectPhase );
          val *= params->testInjectAmpl;
          data->channel->data->data[i] += val;
        }
        else /* ringdown has decayed enough to ignore */
        {
          break;
        }
      }
    }
  }


  /*
   *
   * Segment the data.  Multiply segments by response and invSpectrum.
   *
   */

  params->dataSegment = LALCalloc( params->numSegments,
      sizeof( *params->dataSegment ) );
  if ( ! params->dataSegment )
  {
    ABORT( status, RINGSEARCHH_EALOC, RINGSEARCHH_MSGEALOC );
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

    /* keep the same name as the original data so it can be extracted later */
    strncpy( params->dataSegment[i].name, ser.name,
        sizeof( params->dataSegment[i].name ) );

    /* fft the dummy time series */
    LALTimeFreqRealFFT( status->statusPtr, params->dataSegment + i, &ser,
        params->forwardPlan );
    CHECKSTATUSPTR( status );

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

    /* now multiply by the inverse spectrum */
    LALSCVectorMultiply( status->statusPtr, params->dataSegment[i].data,
        params->invSpectrum->data, params->dataSegment[i].data );
    CHECKSTATUSPTR( status );

    /* get units right again after this multiplication */
    unitPair.unitOne = &params->dataSegment[i].sampleUnits;
    unitPair.unitTwo = &params->invSpectrum->sampleUnits;
    LALUnitMultiply( status->statusPtr, &params->dataSegment[i].sampleUnits,
        &unitPair );
    CHECKSTATUSPTR( status );
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
