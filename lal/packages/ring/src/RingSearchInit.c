#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/RealFFT.h>
#include <lal/Ring.h>
#include <lal/RingSearch.h>

NRCSID( RINGSEARCHINITC, "$Id$" );

void LALRingSearchInit(
    LALStatus         *status,
    RingSearchParams **searchParams,
    const CHAR       **argv,
    INT4               argc
    )
{
  RingSearchParams      *params;
  RingTemplateBankInput  bankin;

  INITSTATUS( status, "LALRingSearchInit", RINGSEARCHINITC );
  ATTATCHSTATUSPTR( status );

  ASSERT( argc < 1 ||  argv, status, RINGSEARCHH_ENULL, RINGSEARCHH_MSGENULL );
  ASSERT( argc < 1 || *argv, status, RINGSEARCHH_ENULL, RINGSEARCHH_MSGENULL );
  ASSERT( searchParams, status, RINGSEARCHH_ENULL, RINGSEARCHH_MSGENULL );
  ASSERT( ! *searchParams, status, RINGSEARCHH_ENNUL, RINGSEARCHH_MSGENNUL );

  params = *searchParams = LALCalloc( 1, sizeof( *params ) );
  if ( ! params )
  {
    ABORT( status, RINGSEARCHH_EALOC, RINGSEARCHH_MSGEALOC );
  }

  params->numSlaves      = -1;
  params->myProcNumber   = -1;
  params->dynRangeFac    =  1;
  params->maximizeEvents = -1; /* natural amount: duration of filter */

  while ( --argc > 0 )
  {
    CHAR *arg;
    ++argv;
    if ( ( arg = strstr( *argv, "segsz:" ) ) )
    {
      params->segmentSize = atoi( strchr( arg, ':' ) + 1 );
    }
    else if ( ( arg = strstr( *argv, "scale:" ) ) )
    {
      params->dynRangeFac = atof( strchr( arg, ':' ) + 1 );
    }
    else if ( ( arg = strstr( *argv, "speclen:" ) ) )
    {
      params->invSpecTrunc = atoi( strchr( arg, ':' ) + 1 );
    }
    else if ( ( arg = strstr( *argv, "flow:" ) ) )
    {
      params->lowFrequency = atof( strchr( arg, ':' ) + 1 );
    }
    else if ( ( arg = strstr( *argv, "fmin:" ) ) )
    {
      params->minFrequency = atof( strchr( arg, ':' ) + 1 );
    }
    else if ( ( arg = strstr( *argv, "fmax:" ) ) )
    {
      params->maxFrequency = atof( strchr( arg, ':' ) + 1 );
    }
    else if ( ( arg = strstr( *argv, "qmin:" ) ) )
    {
      params->minQuality = atof( strchr( arg, ':' ) + 1 );
    }
    else if ( ( arg = strstr( *argv, "qmax:" ) ) )
    {
      params->maxQuality = atof( strchr( arg, ':' ) + 1 );
    }
    else if ( ( arg = strstr( *argv, "maxmm:" ) ) )
    {
      params->maxMismatch = atof( strchr( arg, ':' ) + 1 );
    }
    else if ( ( arg = strstr( *argv, "thresh:" ) ) )
    {
      params->threshold = atof( strchr( arg, ':' ) + 1 );
    }
    else
    {
      ABORT( status, RINGSEARCHH_EIOPT, RINGSEARCHH_MSGEIOPT );
    }
  }

  if ( ! params->segmentSize )
  {
    ABORT( status, RINGSEARCHH_ESIZE, RINGSEARCHH_MSGESIZE );
  }

  if ( ! ( params->lowFrequency > 0 ) )
  {
    ABORT( status, RINGSEARCHH_EFLOW, RINGSEARCHH_MSGEFLOW );
  }

  if ( ! ( params->minFrequency > 0
        && params->maxFrequency > params->minFrequency ) )
  {
    ABORT( status, RINGSEARCHH_EFREQ, RINGSEARCHH_MSGEFREQ );
  }

  if ( ! ( params->minQuality > 0
        && params->maxQuality > params->minQuality ) )
  {
    ABORT( status, RINGSEARCHH_EQUAL, RINGSEARCHH_MSGEQUAL );
  }


  /* memory for invSpectrum */

  params->invSpectrum = LALCalloc( 1, sizeof( *params->invSpectrum ) );
  if ( ! params->invSpectrum )
  {
    ABORT( status, RINGSEARCHH_EALOC, RINGSEARCHH_MSGEALOC );
  }

  LALSCreateVector( status->statusPtr, &params->invSpectrum->data,
      params->segmentSize / 2 + 1 );
  CHECKSTATUSPTR( status );

  LALCreateForwardRealFFTPlan( status->statusPtr, &params->forwardPlan,
      params->segmentSize, 0 );
  CHECKSTATUSPTR( status );
  
  LALCreateReverseRealFFTPlan( status->statusPtr, &params->reversePlan,
      params->segmentSize, 0 );
  CHECKSTATUSPTR( status );

  bankin.minQuality   = params->minQuality;
  bankin.maxQuality   = params->maxQuality;
  bankin.minFrequency = params->minFrequency;
  bankin.maxFrequency = params->maxFrequency;
  bankin.maxMismatch  = params->maxMismatch;
  LALCreateRingTemplateBank( status->statusPtr, &params->templateBank,
      &bankin );
  CHECKSTATUSPTR( status );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


void
LALRingSearchFini(
    LALStatus         *status,
    RingSearchParams **searchParams
    )
{
  RingSearchParams *params;

  INITSTATUS( status, "LALRingSearchFini", RINGSEARCHINITC );
  ATTATCHSTATUSPTR( status );

  ASSERT(  searchParams, status, RINGSEARCHH_ENULL, RINGSEARCHH_MSGENULL );
  ASSERT( *searchParams, status, RINGSEARCHH_ENULL, RINGSEARCHH_MSGENULL );

  params = *searchParams;

  while ( params->numResults )
  {
    --params->numResults;
    if ( params->result[params->numResults].data )
    {
      LALSDestroyVector( status->statusPtr,
          &params->result[params->numResults].data );
      CHECKSTATUSPTR( status );
    }
  }
  if ( params->result )
  {
    LALFree( params->result );
    params->result = NULL;
  }

  while ( params->numSegments )
  {
    --params->numSegments;
    if ( params->dataSegment[params->numSegments].data )
    {
      LALCDestroyVector( status->statusPtr,
          &params->dataSegment[params->numSegments].data );
      CHECKSTATUSPTR( status );
    }
  }
  if ( params->dataSegment )
  {
    LALFree( params->dataSegment );
    params->dataSegment = NULL;
  }

  if ( params->templateBank )
  {
    LALDestroyRingTemplateBank( status->statusPtr, &params->templateBank );
    CHECKSTATUSPTR( status );
  }

  if ( params->reversePlan )
  {
    LALDestroyRealFFTPlan( status->statusPtr, &params->reversePlan );
    CHECKSTATUSPTR( status );
  }

  if ( params->forwardPlan )
  {
    LALDestroyRealFFTPlan( status->statusPtr, &params->forwardPlan );
    CHECKSTATUSPTR( status );
  }
  
  if ( params->invSpectrum )
  {
    if ( params->invSpectrum->data )
    {
      LALSDestroyVector( status->statusPtr, &params->invSpectrum->data );
      CHECKSTATUSPTR( status );
    }
    LALFree( params->invSpectrum );
    params->invSpectrum = NULL;
  }

  LALFree( params );
  params = NULL;

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
