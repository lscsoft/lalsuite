/**** <lalVerbatim file="RingSearchCV">
 * Author: Jolien Creighton
 * $Id$
 **** </lalVerbatim> */

#include <string.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/TimeFreqFFT.h>
#include <lal/VectorOps.h>
#include <lal/Comm.h>
#include <lal/RingSearch.h>

/**** <lalLaTeX>
 *
 * \subsection{Module \texttt{RingSearch.c}}
 *
 * Routine to perform a ring search.
 *
 * \subsubsection*{Prototypes}
 * \input{RingSearchCP}
 * \idx{LALRingSearch()}
 *
 * \subsubsection*{Description}
 * 
 * The function \verb+LALRingSearch()+ performs a ring search over a specified
 * range of templates in the bank and returns a linked list of events where
 * the filter crosses a specified threshold.
 *
 * \vfill{\footnotesize\input{RingCV}}
 *
 **** </lalLaTeX> */



NRCSID( RINGSEARCHC, "$Id$" );



/* <lalVerbatim file="RingSearchCP"> */
void
LALRingSearch(
    LALStatus         *status,
    SnglBurstTable   **output,
    RingSearchInput   *input,
    RingSearchParams  *params
    )
{ /* </lalVerbatim> */
  SnglBurstTable          *thisEvent = NULL;
  REAL4TimeSeries          signal;
  REAL4TimeSeries          result;
  COMPLEX8FrequencySeries  stilde;
  COMPLEX8FrequencySeries  rtilde;
  UINT4                    ntmplt;
  UINT4                    tmplt;

  INITSTATUS( status, "LALRingSearch", RINGSEARCHC );
  ATTATCHSTATUSPTR( status );

  ASSERT(     input, status, RINGSEARCHH_ENULL, RINGSEARCHH_MSGENULL );
  ASSERT(    params, status, RINGSEARCHH_ENULL, RINGSEARCHH_MSGENULL );
  ASSERT(    output, status, RINGSEARCHH_ENULL, RINGSEARCHH_MSGENULL );
  ASSERT( ! *output, status, RINGSEARCHH_ENNUL, RINGSEARCHH_MSGENNUL );

  ntmplt = input->templatesToDo;
  tmplt  = input->startTemplate;

  if ( params->keepResults )
  {
    params->numResults = ntmplt * params->numSegments;
    params->result = LALCalloc( params->numResults,
        sizeof( *params->result ) );
    if ( ! params->result )
    {
      ABORT( status, RINGSEARCHH_EALOC, RINGSEARCHH_MSGEALOC );
    }
  }

  strncpy( signal.name, "ringdown signal", sizeof( signal.name ) );
  signal.f0          = 0;
  signal.deltaT      = 1 / params->sampleRate;
  signal.sampleUnits = lalStrainUnit;
  signal.data        = NULL;
  LALSCreateVector( status->statusPtr, &signal.data, params->segmentSize );
  CHECKSTATUSPTR( status );

  strncpy( stilde.name, "fft of ringdown signal", sizeof( stilde.name ) );
  stilde.f0   = 0;
  stilde.data = NULL;
  LALCCreateVector( status->statusPtr, &stilde.data,
      params->segmentSize / 2 + 1 );
  CHECKSTATUSPTR( status );

  strncpy( rtilde.name, "fft of filter results", sizeof( rtilde.name ) );
  rtilde.f0     = 0;
  rtilde.deltaF = params->sampleRate / params->segmentSize;
  rtilde.data   = NULL;
  LALCCreateVector( status->statusPtr, &rtilde.data,
      params->segmentSize / 2 + 1 );
  CHECKSTATUSPTR( status );

  strncpy( result.name, "filter results", sizeof( result.name ) );
  result.f0   = 0;
  result.data = NULL;
  if ( ! params->keepResults )
  {
    LALSCreateVector( status->statusPtr, &result.data, params->segmentSize );
    CHECKSTATUSPTR( status );
  }

  while ( ntmplt-- > 0 )
  {
    const REAL4 efolds = 10;
    REAL4 threshold;
    REAL4 duration;
    REAL4 sigma;
    REAL4 ssq;
    UINT4 seg;
    UINT4 k;

    LALComputeRingTemplate( status->statusPtr, &signal,
        params->templateBank->tmplt + tmplt );
    CHECKSTATUSPTR( status );

    duration  = efolds * params->templateBank->tmplt[tmplt].quality;
    duration /= LAL_PI * params->templateBank->tmplt[tmplt].frequency;
    duration += params->invSpecTrunc * signal.deltaT;

    LALTimeFreqRealFFT( status->statusPtr, &stilde, &signal,
        params->forwardPlan );
    CHECKSTATUSPTR( status );

    ssq = 0;
    for ( k = 0; k < stilde.data->length; ++k )
    {
      REAL4 re = stilde.data->data[k].re;
      REAL4 im = stilde.data->data[k].im;
      ssq += params->invSpectrum->data->data[k] * ( re * re + im * im );
    }
    sigma = 2 * params->dynRangeFac * sqrt( ssq * stilde.deltaF );
    threshold = 0.5 * sigma * params->threshold / params->dynRangeFac;

    for ( seg = 0; seg < params->numSegments; ++seg )
    {
      INT8        lastTimeNS = -1e11; /* before any likely time */
      INT8        gapTimeNS;
      LALUnitPair unitPair;
      UINT4 j;

      if ( params->maximizeEvents < 0 )
      {
        gapTimeNS = (INT8)( 1e9 * duration );
      }
      else if ( params->maximizeEvents > 0 )
      {
        gapTimeNS = (INT8)( 1e9 * params->maximizeEvents / params->sampleRate );
      }
      else
      {
        gapTimeNS = 0;
      }
            
      LALCCVectorMultiplyConjugate( status->statusPtr, rtilde.data,
          params->dataSegment[seg].data, stilde.data );
      CHECKSTATUSPTR( status );

      unitPair.unitOne = &params->dataSegment[seg].sampleUnits;
      unitPair.unitTwo = &stilde.sampleUnits;
      LALUnitMultiply( status->statusPtr, &rtilde.sampleUnits, &unitPair );
      CHECKSTATUSPTR( status );

      rtilde.epoch = params->dataSegment[seg].epoch;

      if ( params->keepResults )
      {
        result.data = NULL;
        LALSCreateVector( status->statusPtr, &result.data,
            params->segmentSize );
        CHECKSTATUSPTR( status );
      }

      LALFreqTimeRealFFT( status->statusPtr, &result, &rtilde,
          params->reversePlan );
      CHECKSTATUSPTR( status );

      /* search for threshold crossing in middle of segment */
      for ( j = result.data->length / 4; j < 3 * result.data->length / 4; ++j )
      {
        if ( fabs( result.data->data[j] ) > threshold )
        {
          REAL4 snr;
          INT8  timeNS;
          snr     = 2 * params->dynRangeFac * fabs( result.data->data[j] );
          snr    /= sigma;
          timeNS  = (INT8)( 1000000000 ) * (INT8)( result.epoch.gpsSeconds );
          timeNS += (INT8)( result.epoch.gpsNanoSeconds );
          timeNS += (INT8)( 1e9 * result.deltaT * j );
          if ( timeNS > lastTimeNS + gapTimeNS ) /* new event */
          {
            ++params->numEvents;
            if ( *output ) /* create a new event */
            {
              thisEvent->next = LALCalloc( 1, sizeof( *thisEvent->next ) );
              thisEvent = thisEvent->next;
            }
            else /* create the list */
            {
              thisEvent = *output = LALCalloc( 1, sizeof( **output ) );
            }
            if ( ! thisEvent ) /* allocation error */
            {
              ABORT( status, RINGSEARCHH_EALOC, RINGSEARCHH_MSGEALOC );
            }
            strncpy( thisEvent->ifo, params->ifoName,
                sizeof( thisEvent->ifo ) );
            strncpy( thisEvent->search, "ring",
                sizeof( thisEvent->search ) );
            strncpy( thisEvent->channel, params->dataSegment[seg].name,
                sizeof( thisEvent->channel ) );
            thisEvent->start_time.gpsSeconds     = timeNS / 1000000000;
            thisEvent->start_time.gpsNanoSeconds = timeNS % 1000000000;
            thisEvent->peak_time    = thisEvent->start_time;
            thisEvent->duration     = duration;
            thisEvent->central_freq =
              params->templateBank->tmplt[tmplt].frequency;
            thisEvent->bandwidth    = 
              params->templateBank->tmplt[tmplt].frequency /
              params->templateBank->tmplt[tmplt].quality;
            thisEvent->snr          = snr;
            thisEvent->amplitude    = snr / sigma;
            thisEvent->confidence   = 0; /* FIXME */
          }
          else /* maximize within existing event */
          {
            if ( snr > thisEvent->snr )
            {
              thisEvent->start_time.gpsSeconds     = timeNS / 1000000000;
              thisEvent->start_time.gpsNanoSeconds = timeNS % 1000000000;
              thisEvent->peak_time  = thisEvent->start_time;
              thisEvent->snr        = snr;
              thisEvent->amplitude  = snr / sigma;
              thisEvent->confidence = 0; /* FIXME */
            }
          }
          /* update last event time */
          lastTimeNS = timeNS;
        }
      }

      if ( params->keepResults )
      {
        sprintf( result.name, "snr-%d.%03d", tmplt, seg );
        for ( j = 0; j < result.data->length; ++j )
        {
          result.data->data[j] *= 2 * params->dynRangeFac / sigma;
        }
        params->result[tmplt * params->numSegments + seg] = result;
      }

    }

    ++tmplt;
  }

  if ( ! params->keepResults )
  {
    LALSDestroyVector( status->statusPtr, &result.data );
    CHECKSTATUSPTR( status );
  }

  LALCDestroyVector( status->statusPtr, &rtilde.data );
  CHECKSTATUSPTR( status );

  LALCDestroyVector( status->statusPtr, &stilde.data );
  CHECKSTATUSPTR( status );

  LALSDestroyVector( status->statusPtr, &signal.data );
  CHECKSTATUSPTR( status );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
