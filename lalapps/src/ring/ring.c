#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include <lal/LALStdlib.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/AVFactories.h>
#include <lal/RealFFT.h>
#include <lal/Ring.h>
#include <lal/FrameStream.h>

#include "lalapps.h"
#include "getdata.h"
#include "injsgnl.h"
#include "getresp.h"
#include "spectrm.h"
#include "segment.h"
#include "errutil.h"

#include "ring.h"

RCSID( "$Id$" );
#define PROGRAM_NAME "lalapps_ring"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE   "$Source$"
#define CVS_DATE     "$Date$"

static struct ring_params *ring_get_params( int argc, char **argv );
static REAL4FFTPlan *ring_get_fft_fwdplan( struct ring_params *params );
static REAL4FFTPlan *ring_get_fft_revplan( struct ring_params *params );
static REAL4TimeSeries *ring_get_data( struct ring_params *params );
static COMPLEX8FrequencySeries *ring_get_response( struct ring_params *params );
static REAL4FrequencySeries *ring_get_invspec(
    REAL4TimeSeries         *channel,
    COMPLEX8FrequencySeries *response,
    REAL4FFTPlan            *fwdplan,
    REAL4FFTPlan            *revplan,
    struct ring_params      *params
    );
static RingTemplateBank *ring_get_bank( struct ring_params *params );
static UINT4 *ring_get_tmpltnums(
    UINT4              *numtmplts,
    RingTemplateBank   *bank,
    struct ring_params *params
    );
static UINT4 *ring_get_sgmntnums(
    UINT4              *numsgmnts,
    struct ring_params *params
    );
static int is_in_list( int i, const char *list );
static void ring_cleanup(
    ProcessParamsTable      *procpar,
    REAL4FFTPlan            *fwdplan,
    REAL4FFTPlan            *revplan,
    REAL4TimeSeries         *channel,
    COMPLEX8FrequencySeries *response,
    REAL4FrequencySeries    *invspec,
    RingTemplateBank        *bank,
    UINT4                   *tmpltnums,
    UINT4                   *sgmntnums,
    UINT4                    numsgmnts,
    COMPLEX8FrequencySeries *segments,
    SnglBurstTable          *events
    );


int main( int argc, char **argv )
{
  struct ring_params      *params    = NULL;
  ProcessParamsTable      *procpar   = NULL;
  REAL4FFTPlan            *fwdplan   = NULL;
  REAL4FFTPlan            *revplan   = NULL;
  REAL4TimeSeries         *channel   = NULL;
  COMPLEX8FrequencySeries *response  = NULL;
  REAL4FrequencySeries    *invspec   = NULL;
  RingTemplateBank        *bank      = NULL;
  UINT4                    numsgmnts;
  UINT4                   *sgmntnums = NULL;
  UINT4                    numtmplts;
  UINT4                   *tmpltnums = NULL;
  COMPLEX8FrequencySeries *segments  = NULL;
  SnglBurstTable          *events    = NULL;

  /* set error handlers to abort on error */
  set_abrt_on_error();

  /* options are parsed and debug level is set here... */
  /* no lal mallocs before this! */
  params = ring_get_params( argc, argv );

  /* create process params */
  procpar = create_process_params( argc, argv, PROGRAM_NAME );

  /* create forward and reverse fft plans */
  fwdplan = ring_get_fft_fwdplan( params );
  revplan = ring_get_fft_revplan( params );

  /* get the data */
  channel = ring_get_data( params );

  /* get the response */
  response = ring_get_response( params );

  /* compute the spectrum */
  invspec = ring_get_invspec( channel, response, fwdplan, revplan, params );

  /* create the template bank */
  bank = ring_get_bank( params );

  /* get template and segment numbers to do */
  tmpltnums = ring_get_tmpltnums( &numtmplts, bank, params );
  sgmntnums = ring_get_sgmntnums( &numsgmnts, params );

  /* create the segments to do */
  segments = compute_data_segments( numsgmnts, sgmntnums, channel, invspec,
      response, params->segmentDuration, params->strideDuration, fwdplan );

  /* filter the data against the bank of templates */
  events = ring_filter( segments, numsgmnts, bank, numtmplts, tmpltnums,
      invspec, fwdplan, revplan, params );

  /* output the results */
  ring_output_events( events, procpar, params );

  /* cleanup */
  ring_cleanup( procpar, fwdplan, revplan, channel, response, invspec, bank,
      tmpltnums, sgmntnums, numsgmnts, segments, events );
  LALCheckMemoryLeaks();

  return 0;
}


/* warning: returns a pointer to a static variable... not reenterant */
/* only call this routine once to initialize params! */
/* also do not attempt to free this pointer! */
static struct ring_params *ring_get_params( int argc, char **argv )
{
  static struct ring_params params;
  static char programName[] = PROGRAM_NAME;
  static char cvsRevision[] = CVS_REVISION;
  static char cvsSource[]   = CVS_SOURCE;
  static char cvsDate[]     = CVS_DATE;
  ring_parse_options( &params, argc, argv );
  ring_params_sanity_check( &params ); /* this also sets various params */
  params.programName = programName;
  params.cvsRevision = cvsRevision;
  params.cvsSource   = cvsSource;
  params.cvsDate     = cvsDate;
  return &params;
}


static REAL4FFTPlan *ring_get_fft_fwdplan( struct ring_params *params )
{
  REAL4FFTPlan *plan = NULL;
  if ( params->segmentDuration > 0.0 )
  {
    UINT4 segmentLength;
    segmentLength = floor( params->segmentDuration * params->sampleRate + 0.5 );
    plan = XLALCreateForwardREAL4FFTPlan( segmentLength, 0 );
  }
  return plan;
}


static REAL4FFTPlan *ring_get_fft_revplan( struct ring_params *params )
{
  REAL4FFTPlan *plan = NULL;
  if ( params->segmentDuration > 0.0 )
  {
    UINT4 segmentLength;
    segmentLength = floor( params->segmentDuration * params->sampleRate + 0.5 );
    plan = XLALCreateReverseREAL4FFTPlan( segmentLength, 0 );
  }
  return plan;
}


static REAL4TimeSeries *ring_get_data( struct ring_params *params )
{
  REAL4TimeSeries *channel = NULL;

  if ( params->getData )
  {
    if ( params->simData )
      channel = get_simulated_data( params->channel, &params->startTime,
          params->duration, params->strainData, params->sampleRate, 
          params->randomSeed, 1.0 );
    else if ( params->geoData )
      channel = get_frame_data_dbl_convert( params->dataCache, params->channel,
          LAL_ADC_CHAN, &params->startTime, params->duration,
          params->strainData, params->geoHighpassFrequency, params->geoScale );
    else
      channel = get_frame_data( params->dataCache, params->channel,
          LAL_ADC_CHAN, &params->startTime, params->duration,
          params->strainData );
    if ( params->writeRawData ) /* write raw data */
      write_REAL4TimeSeries( channel );

    /* inject burst signals */
    if ( params->injectFile )
      inject_signal( channel, burst_inject, params->injectFile,
          params->calibCache, params->dynRangeFac );

    /* condition the data: resample and highpass */
    highpass_REAL4TimeSeries( channel, params->highpassFrequency );
    resample_REAL4TimeSeries( channel, params->sampleRate );

    if ( params->writeProcessedData ) /* write processed data */
      write_REAL4TimeSeries( channel );
  }
  
  return channel;
}


static COMPLEX8FrequencySeries *ring_get_response( struct ring_params *params )
{
  COMPLEX8FrequencySeries *response = NULL;
  if ( params->getResponse && ! params->strainData )
  {
    response = get_response( params->calibCache, params->ifoName,
        &params->startTime, params->segmentDuration, params->sampleRate,
        params->dynRangeFac, params->strainData );
    if ( params->writeResponse ) /* write response */
      write_COMPLEX8FrequencySeries( response );
  }
  return response;
}


static REAL4FrequencySeries *ring_get_invspec(
    REAL4TimeSeries         *channel,
    COMPLEX8FrequencySeries *response,
    REAL4FFTPlan            *fwdplan,
    REAL4FFTPlan            *revplan,
    struct ring_params      *params
    )
{
  REAL4FrequencySeries *invspec = NULL;
  if ( params->getSpectrum )
  {
    /* compute raw average spectrum; store spectrum in invspec for now */
    invspec = compute_average_spectrum( channel, params->segmentDuration,
        params->strideDuration, fwdplan, params->whiteSpectrum );
    if ( params->writeSpectrum ) /* write raw spectrum */
      write_REAL4FrequencySeries( invspec );

    /* invert spectrum */
    invert_spectrum( invspec, params->sampleRate, params->strideDuration,
        params->truncateDuration, params->lowCutoffFrequency, fwdplan,
        revplan );

    /* calibrate spectrum if necessary */
    if ( response )
      calibrate_spectrum( invspec, response, params->lowCutoffFrequency, 1 );

    if ( params->writeInvSpectrum ) /* write inverse calibrated spectrum */
      write_REAL4FrequencySeries( invspec );
  }

  return invspec;
}


static RingTemplateBank *ring_get_bank( struct ring_params *params )
{
  RingTemplateBank *bank = NULL;
  if ( params->getBank )
  {
    LALStatus status = blank_status;
    LAL_CALL( LALCreateRingTemplateBank( &status, &bank, &params->bankParams ),
        &status );
    if ( params->writeBank ) /* write template bank */
      write_bank( bank );
  }
  return bank;
}


static UINT4 *ring_get_tmpltnums(
    UINT4              *numtmplts,
    RingTemplateBank   *bank,
    struct ring_params *params
    )
{
  UINT4 *tmpltnums = NULL;
  UINT4  count = 0;
  UINT4  tmplt;

  /* if todo list is empty then do them all */
  if ( ! params->templatesToDoList || ! strlen( params->templatesToDoList ) )
  {
    *numtmplts = bank->numTmplt;
    return NULL;
  }

  /* first count the number of templates to do */
  for ( tmplt = 0; tmplt < bank->numTmplt; ++tmplt )
    if ( is_in_list( tmplt, params->templatesToDoList ) )
      ++count;
    else
      continue; /* skip this template: it is not in todo list */
  
  *numtmplts = count;

  /* now get the template numbers */
  if ( count )
  {
    UINT4 *thistmpltnum;
    thistmpltnum = tmpltnums = LALCalloc( count, sizeof( *tmpltnums ) );

    for ( tmplt = 0; tmplt < bank->numTmplt; ++tmplt )
      if ( is_in_list( tmplt, params->templatesToDoList ) )
        *thistmpltnum++ = tmplt;
  }

  return tmpltnums;
}


static UINT4 *ring_get_sgmntnums(
    UINT4              *numsgmnts,
    struct ring_params *params
    )
{
  UINT4 *sgmntnums = NULL;
  UINT4  count = 0;
  UINT4  sgmnt;

  /* TODO: trig start/end time condition */

  /* if todo list is empty then do them all */
  if ( ! params->segmentsToDoList || ! strlen( params->segmentsToDoList ) )
  {
    *numsgmnts = params->numOverlapSegments;
    return NULL;
  }

  /* first count the number of segments to do */
  for ( sgmnt = 0; sgmnt < params->numOverlapSegments; ++sgmnt )
    if ( is_in_list( sgmnt, params->segmentsToDoList ) )
      ++count;
    else
      continue; /* skip this segment: it is not in todo list */
  
  *numsgmnts = count;

  /* now get the segment numbers */
  if ( count )
  {
    UINT4 *thissgmntnum;
    thissgmntnum = sgmntnums = LALCalloc( count, sizeof( *sgmntnums ) );

    for ( sgmnt = 0; sgmnt < params->numOverlapSegments; ++sgmnt )
      if ( is_in_list( sgmnt, params->segmentsToDoList ) )
        *thissgmntnum++ = sgmnt;
  }

  return sgmntnums;
}


static void ring_cleanup(
    ProcessParamsTable      *procpar,
    REAL4FFTPlan            *fwdplan,
    REAL4FFTPlan            *revplan,
    REAL4TimeSeries         *channel,
    COMPLEX8FrequencySeries *response,
    REAL4FrequencySeries    *invspec,
    RingTemplateBank        *bank,
    UINT4                   *tmpltnums,
    UINT4                   *sgmntnums,
    UINT4                    numsgmnts,
    COMPLEX8FrequencySeries *segments,
    SnglBurstTable          *events
    )
{
  LALStatus status = blank_status;
  while ( events )
  {
    SnglBurstTable *thisEvent;
    thisEvent = events;
    events = events->next;
    LALFree( thisEvent );
  }
  if ( segments )
  {
    UINT4 sgmnt;
    for ( sgmnt = 0; sgmnt < numsgmnts; ++sgmnt )
      LAL_CALL( LALCDestroyVector( &status, &segments[sgmnt].data ), &status );
    LALFree( segments );
  }
  if ( sgmntnums )
    LALFree( sgmntnums );
  if ( tmpltnums )
    LALFree( tmpltnums );
  if ( bank )
  {
    LAL_CALL( LALDestroyRingTemplateBank( &status, &bank ), &status );
  }
  if ( invspec )
  {
    LAL_CALL( LALSDestroyVector( &status, &invspec->data ), &status );
    LALFree( invspec );
  }
  if ( response )
  {
    LAL_CALL( LALCDestroyVector( &status, &response->data ), &status );
    LALFree( response );
  }
  if ( channel )
  {
    LAL_CALL( LALSDestroyVector( &status, &channel->data ), &status );
    LALFree( channel );
  }
  if ( revplan )
    XLALDestroyREAL4FFTPlan( revplan );
  if ( fwdplan )
    XLALDestroyREAL4FFTPlan( fwdplan );
  while ( procpar )
  {
    ProcessParamsTable *thisParam;
    thisParam = procpar;
    procpar = procpar->next;
    LALFree( thisParam );
  }
  return;
}


#define BUFFER_SIZE 256
static int is_in_list( int i, const char *list )
{
  char  buffer[BUFFER_SIZE];
  char *str = buffer;
  int   ans = 0;

  strncpy( buffer, list, sizeof( buffer ) - 1 );

  while ( str )
  {
    char *tok;  /* token in a delimited list */
    char *tok2; /* second part of token if it is a range */
    tok = str;

    if ( ( str = strchr( str, ',' ) ) ) /* look for next delimiter */
      *str++ = 0; /* nul terminate current token; str is remaining string */

    /* now see if this token is a range */
    if ( ( tok2 = strchr( tok, '-' ) ) )
      *tok2++ = 0; /* nul terminate first part of token; tok2 is second part */

    if ( tok2 ) /* range */
    {
      int n1, n2;
      if ( strcmp( tok, "^" ) == 0 ) 
        n1 = INT_MIN;
      else
        n1 = atoi( tok );
      if ( strcmp( tok2, "$" ) == 0 )
        n2 = INT_MAX;
      else
        n2 = atoi( tok2 );
      if ( i >= n1 && i <= n2 ) /* see if i is in the range */
        ans = 1;
    }
    else if ( i == atoi( tok ) )
      ans = 1;

    if ( ans ) /* i is in the list */
      break;
  }

  return ans;
}
