#include <lal/String.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/TimeFreqFFT.h>
#include <lal/VectorOps.h>
#include <lal/Comm.h>
#include <lal/StringSearch.h>

NRCSID( STRINGSEARCHC, "$Id$" );

void
LALStringSearch(
    LALStatus         *status,
    StringEventList    **output,
    StringSearchInput   *input,
    StringSearchParams  *params
    )
{
  StringEventList           *thisEvent = NULL;
  REAL4TimeSeries          result;
  COMPLEX8FrequencySeries  stilde;
  COMPLEX8FrequencySeries  rtilde;
  UINT4                    ntmplt;
  UINT4                    tmplt;

  INITSTATUS( status, "StringSearch", STRINGSEARCHC );
  ATTATCHSTATUSPTR( status );

  ASSERT(     input, status, STRINGSEARCHH_ENULL, STRINGSEARCHH_MSGENULL );
  ASSERT(    params, status, STRINGSEARCHH_ENULL, STRINGSEARCHH_MSGENULL );
  ASSERT(    output, status, STRINGSEARCHH_ENULL, STRINGSEARCHH_MSGENULL );
  ASSERT( ! *output, status, STRINGSEARCHH_ENNUL, STRINGSEARCHH_MSGENNUL );

  ntmplt = input->templatesToDo;
  tmplt  = input->startTemplate;

  if ( params->keepResults )
  {
    params->numResults = ntmplt * params->numSegments;
    params->result = LALCalloc( params->numResults,
        sizeof( *params->result ) );
    if ( ! params->result )
    {
      ABORT( status, STRINGSEARCHH_EALOC, STRINGSEARCHH_MSGEALOC );
    }
  }
  
  strncpy( stilde.name, "freq domain string signal", sizeof( stilde.name ) );
  stilde.f0   = 0;
  stilde.deltaF = params->sampleRate / params->segmentSize;
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

  LALComputeFreqStringTemplate( status->statusPtr, &stilde, params->freqPower);
    CHECKSTATUSPTR( status );


  while ( ntmplt-- > 0 )
  {
    REAL4 cutoff = params->templateBank->tmplt[tmplt].frequency 
      / stilde.deltaF;
    REAL4 threshold;
    REAL4 ssq = 0;
    REAL4 sigma;
    UINT4 seg;
    UINT4 k;

    for ( k= 0; k < cutoff; ++k)
    {
      ssq+= params->invSpectrum->data->data[k] * stilde.data->data[k].re
	* stilde.data->data[k].re;
    }
    sigma = 2 * params->dynRangeFac * sqrt(ssq*stilde.deltaF);
    threshold = 0.5 * sigma * params->threshold / params->dynRangeFac;

    for (seg = 0; seg<params->numSegments; ++seg )
    {
      INT8        lastTimeNS = -1e11; /* before any likely time */
      INT8        gapTimeNS;
      /* LALUnitPair unitPair; */
      UINT4 j;
      UINT4 n;
      
      if ( params->maximizeEvents < 0 )
      {
        gapTimeNS = 0;
      }
      else if ( params->maximizeEvents > 0 )
      {
        gapTimeNS = (INT8)( 1e9 * params->maximizeEvents / params->sampleRate );
      }
      else
      {
        gapTimeNS = 0;
      }

      /*conjugate multiplication with cutoff*/
      for( n=0; n<cutoff;++n)
      { /*   rtilde = stilde * conj(params->dataSegment[seg].data)   */
	REAL4 ar = stilde.data->data[n].re;
	REAL4 ai = stilde.data->data[n].im;
	REAL4 br = params->dataSegment[seg].data->data[n].re;
	REAL4 bi = params->dataSegment[seg].data->data[n].im;
	rtilde.data->data[n].re = ai*bi* + ar*br;
	rtilde.data->data[n].im = ai*br - ar*bi;
      }
      for( n=cutoff; n<rtilde.data->length; ++n)
      {
	rtilde.data->data[n].re = 0;
	rtilde.data->data[n].im = 0;
      }/*???*/
/*
      unitPair.unitOne = &params->dataSegment[seg].sampleUnits;
      unitPair.unitTwo = &stilde.sampleUnits;
      LALUnitMultiply( status->statusPtr, &rtilde.sampleUnits, &unitPair );
  */  
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
      
      /* search for threshold crossing XXX? in middle of segment XXX?*/
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
              ABORT( status, STRINGSEARCHH_EALOC, STRINGSEARCHH_MSGEALOC );
            }
            thisEvent->startTimeNS = timeNS;
            thisEvent->snr         = snr;
            thisEvent->amplitude   = snr / sigma;
            thisEvent->confidence  = 0; /* FIXME */
            thisEvent->duration    = 0;
            thisEvent->frequency   = 
	      params->templateBank->tmplt[tmplt].frequency;
            thisEvent->mass        = 0; /* FIXME */
            strncpy( thisEvent->ifoName, params->ifoName,
                sizeof( thisEvent->ifoName ) );
          }
          else /* maximize within existing event */
          {
            if ( snr > thisEvent->snr )
            {
              thisEvent->startTimeNS = timeNS;
              thisEvent->snr         = snr;
              thisEvent->amplitude   = snr / sigma;
              thisEvent->confidence  = 0; /* FIXME */
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

  DETATCHSTATUSPTR( status );
  RETURN( status );
}




    
    

	

      


