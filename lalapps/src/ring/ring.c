/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Lisa M. Goggin, Stephen Fairhurst, Tania Regimbau
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include <lal/LALStdlib.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataRingdownUtils.h>
#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/RealFFT.h>
#include <lal/RingUtils.h>
#include <lal/LALFrStream.h>

#include "lalapps.h"
#include "getdata.h"
#include "injsgnl.h"
#include "getresp.h"
#include "spectrm.h"
#include "segment.h"
#include "errutil.h"

#include "ring.h"

#define PROGRAM_NAME "lalapps_ring"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE   "$Source$"
#define CVS_DATE     "$Date$"

/* a bunch of routines that are ring search specific */
/* such routines are prefixed with "ring_" */
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
static RingDataSegments *ring_get_segments(
    REAL4TimeSeries         *channel,
    COMPLEX8FrequencySeries *response,
    REAL4FrequencySeries    *invspec,
    REAL4FFTPlan            *fwdplan,
    struct ring_params      *params
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
    RingDataSegments        *segments,
    SnglRingdownTable       *events
    );


int main( int argc, char **argv )
{
  static LALStatus      status;
  struct ring_params      *params    = NULL;
  ProcessParamsTable      *procpar   = NULL;
  REAL4FFTPlan            *fwdplan   = NULL;
  REAL4FFTPlan            *revplan   = NULL;
  REAL4TimeSeries         *channel   = NULL;
  COMPLEX8FrequencySeries *response  = NULL;
  REAL4FrequencySeries    *invspec   = NULL;
  RingTemplateBank        *bank      = NULL;
  RingDataSegments        *segments  = NULL;
  SnglRingdownTable       *events    = NULL;
  SnglRingdownTable *tmpEventHead = NULL;
  SnglRingdownTable *checkEvents = NULL;
  SnglRingdownTable *lastEvent = NULL;
  UINT4 numEvents = 0;

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
  if ( params->bankFile[0] ) /* write out the bank */
    ring_output_events_xml( params->bankFile, bank->tmplt, procpar, params );

  /* create the segments to do */
  segments = ring_get_segments( channel, response, invspec, fwdplan, params );

  /* filter the data against the bank of templates */
  events = ring_filter( segments, bank, invspec, fwdplan, revplan, params );

  /* time sort the triggers */
  if ( vrbflg ) fprintf( stdout, "Sorting triggers\n" );
  LAL_CALL( LALSortSnglRingdown( &status, &(events),
        LALCompareSnglRingdownByTime ), &status );
  
  /* discard any triggers outside the trig start/end time window */
  
  if ( vrbflg ) fprintf( stdout, 
      "  discarding triggers outside trig start/end time... \n" );

  checkEvents = events;
 
  while ( checkEvents )
  {
    INT8 trigTimeNS;
    trigTimeNS = XLALGPSToINT8NS( &(checkEvents->start_time) );
    if ( trigTimeNS &&  ((params->trigStartTimeNS && 
            (trigTimeNS < params->trigStartTimeNS)) ||
          (params->trigEndTimeNS && (trigTimeNS >= params->trigEndTimeNS))) )
    {
      /* throw this trigger away */
      SnglRingdownTable *tmpEvent = checkEvents;

      if ( lastEvent )
      {
        lastEvent->next = checkEvents->next;
      }

      /* increment the linked list by one and free the event */
      checkEvents = checkEvents->next;
      LALFree( tmpEvent );
    }
    else
    {
      /* store the first event as the head of the new linked list */
      if ( ! tmpEventHead ) tmpEventHead = checkEvents;

      /* save the last event and increment the linked list by one */
      lastEvent = checkEvents;
     checkEvents = checkEvents->next;
     numEvents++;
    }
  }

  if ( vrbflg ) fprintf( stdout, "%u triggers remaining\n", numEvents );
  events = tmpEventHead;

  /* output the results */
  ring_output_events_xml( params->outputFile, events, procpar, params );

  /* cleanup */
  ring_cleanup( procpar, fwdplan, revplan, channel, response, invspec, bank,
      segments, events );
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


/* gets the forward fft plan */
static REAL4FFTPlan *ring_get_fft_fwdplan( struct ring_params *params )
{
  REAL4FFTPlan *plan = NULL;
  if ( params->segmentDuration > 0.0 )
  {
    UINT4 segmentLength;
    segmentLength = floor( params->segmentDuration * params->sampleRate + 0.5 );
    plan = XLALCreateForwardREAL4FFTPlan( segmentLength, 1 );
  }
  return plan;
}


/* gets the reverse fft plan */
static REAL4FFTPlan *ring_get_fft_revplan( struct ring_params *params )
{
  REAL4FFTPlan *plan = NULL;
  if ( params->segmentDuration > 0.0 )
  {
    UINT4 segmentLength;
    segmentLength = floor( params->segmentDuration * params->sampleRate + 0.5 );
    plan = XLALCreateReverseREAL4FFTPlan( segmentLength, 1 );
  }
  return plan;
}


/* gets the data, performs any injections, and conditions the data */
static REAL4TimeSeries *ring_get_data( struct ring_params *params )
{
  int stripPad = 0;
  REAL4TimeSeries *channel = NULL;
  UINT4 j;

  /* compute the start and duration needed to pad data */
  params->frameDataStartTime = params->startTime;
  XLALGPSAdd( &params->frameDataStartTime, -1.0 * params->padData );
  params->frameDataDuration = params->duration + 2.0 * params->padData;

  /* Needed for simulate data */
  REAL8FrequencySeries *spectrum = NULL;

  if ( params->getData )
  {
    switch ( params->dataType )
    {
      case LALRINGDOWN_DATATYPE_SIM:
        /* Value of 1E-20 is used for white spectrum. It sets the PSD scale
           to give SNRs that are kind of comparable to iLIGO PSD */
        spectrum = generate_theoretical_psd(1./params->sampleRate,
            params->duration,params->simDataType, 1E-20);
        channel = get_simulated_data_new( params->channel, &params->startTime,
            params->duration, params->sampleRate,
            params->randomSeed, spectrum );
        XLALDestroyREAL8FrequencySeries(spectrum);
      break;
      case LALRINGDOWN_DATATYPE_ZERO:
        channel = get_zero_data( params->channel, &params->startTime,
            params->duration, params->dataType, params->sampleRate );
      break;
      case LALRINGDOWN_DATATYPE_HT_REAL8:
        channel = get_frame_data_dbl_convert( params->dataCache, params->channel,
            &params->frameDataStartTime, params->frameDataDuration,
          params->dataType, params->highpassFrequency);
        stripPad = 1;
      break;
      case LALRINGDOWN_DATATYPE_UNCAL: case LALRINGDOWN_DATATYPE_HT_REAL4:
        channel = ring_get_frame_data( params->dataCache, params->channel,
            &params->frameDataStartTime, params->frameDataDuration,
            params->dataType );
        stripPad = 1;
      break;
    }
    if ( params->writeRawData ) /* write raw data */
      write_REAL4TimeSeries( channel );
    
    /* inject ring signals */
    if ( params->injectFile ) 
    {
      ring_inject_signal( channel, params->injectType, params->injectFile,
          params->calibCache, 1.0, params->channel ); 
      if ( params->writeRawData )
        write_REAL4TimeSeries( channel );
    }  

    if ( params->dataType == LALRINGDOWN_DATATYPE_HT_REAL4 ||
         params->dataType == LALRINGDOWN_DATATYPE_HT_REAL8)
    {  
      for (j=0; j<channel->data->length; j++)
      {
        channel->data->data[j] *= params->dynRangeFac;
      }
    }
    
    /* condition the data: resample and highpass */
    resample_REAL4TimeSeries( channel, params->sampleRate );
    if ( params->writeProcessedData ) /* write processed data */
      write_REAL4TimeSeries( channel );

    highpass_REAL4TimeSeries( channel, params->highpassFrequency );
    if ( params->writeProcessedData ) /* write processed data */
      write_REAL4TimeSeries( channel );

    if ( stripPad )
    {
      trimpad_REAL4TimeSeries( channel, params->padData );
      if ( params->writeProcessedData ) /* write data with padding removed */
        write_REAL4TimeSeries( channel );  
    }
  }
  
  return channel;
}


/* gets the up-to-date response function */
static COMPLEX8FrequencySeries *ring_get_response( struct ring_params *params )
{
  COMPLEX8FrequencySeries *response = NULL;
  if ( params->getResponse && (params->dataType != LALRINGDOWN_DATATYPE_HT_REAL4
         && params->dataType != LALRINGDOWN_DATATYPE_HT_REAL8 ))
  {
  response = get_response( params->calibCache, params->ifoName,
        &params->startTime, params->segmentDuration, params->sampleRate,
        params->dynRangeFac, params->dataType, params->channel ) ; 
    if ( params->writeResponse ) /* write response */
      write_COMPLEX8FrequencySeries( response );
  }
  return response;
}


/* computes the inverse power spectrum */
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
    if ( ! params->whiteSpectrum )
    {
      /* compute raw average spectrum; store spectrum in invspec for now */
      invspec = compute_average_spectrum( channel, params->spectrumType,
          params->segmentDuration, params->strideDuration,
          fwdplan, params->whiteSpectrum );
    }
    else
    {
      UINT4 k;
      REAL8FrequencySeries *spectrum;
      spectrum = generate_theoretical_psd(1./params->sampleRate,
          params->segmentDuration,params->simDataType, 1E-20);
      
      /* Need to convert to a REAL4 FrequencySeries */
      UINT4 segmentLength = floor( params->segmentDuration/channel->deltaT\
                                    + 0.5 );
      invspec = XLALCreateREAL4FrequencySeries("TEMP",&(channel->epoch),0,
          1.0/params->segmentDuration,&lalDimensionlessUnit,
          segmentLength/2 + 1);
      snprintf( invspec->name, sizeof( invspec->name),
          "%s_SPEC", channel->name);
      for ( k = 0; k < spectrum->data->length; ++k )
      {
        invspec->data->data[k] = (REAL4)(spectrum->data->data[k]*params->dynRangeFac*params->dynRangeFac);
      }
      XLALDestroyREAL8FrequencySeries(spectrum);
    }

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


/* creates the bank of requested ringdown templates (those in the list to do) */
static RingTemplateBank *ring_get_bank( struct ring_params *params )
{
  RingTemplateBank *bank = NULL;

  if ( params->getBank )
  {
    /* get the complete bank */
    bank = XLALCreateRingTemplateBank( &params->bankParams );

    /* if todo list empty then use the full bank */
    if ( params->templatesToDoList && strlen( params->templatesToDoList ) )
    {
      UINT4 count = 0;
      UINT4 tmplt;

      /* loop over templates in full bank and re-insert the desired ones
       * at the beginning of the bank; simultaneously count the number */
      for ( tmplt = 0; tmplt < bank->numTmplt; ++tmplt )
        if ( is_in_list( tmplt, params->templatesToDoList ) )
          bank->tmplt[count++] = bank->tmplt[tmplt];
      
      /* reallocate memory to the (possibly) smaller size */
      if ( count )
      {
        bank->numTmplt = count;
        bank->tmplt = LALRealloc( bank->tmplt,
            bank->numTmplt * sizeof( *bank->tmplt ) );
        for ( tmplt = 0; tmplt < count - 1; ++tmplt )
          bank->tmplt[tmplt].next = bank->tmplt + tmplt + 1;
        bank->tmplt[count - 1].next = NULL;
      }
      else /* no templates: return an NULL bank */
      {
        LALFree( bank->tmplt );
        LALFree( bank );
        bank = NULL;
      }
    }
  }

  if ( bank )
    verbose( "number of templates in bank = %d\n", bank->numTmplt );

  return bank;
}


/* creates the requested data segments (those in the list of segments to do) */
static RingDataSegments *ring_get_segments(
    REAL4TimeSeries         *channel,
    COMPLEX8FrequencySeries *response,
    REAL4FrequencySeries    *invspec,
    REAL4FFTPlan            *fwdplan,
    struct ring_params      *params
    )
{
  RingDataSegments *segments = NULL;
  UINT4  sgmnt;

  segments = LALCalloc( 1, sizeof( *segments ) );

  /* TODO: trig start/end time condition */

  /* if todo list is empty then do them all */
  if ( ! params->segmentsToDoList || ! strlen( params->segmentsToDoList ) )
  {
    segments->numSgmnt = params->numOverlapSegments;
    segments->sgmnt = LALCalloc( segments->numSgmnt, sizeof(*segments->sgmnt) );
    for ( sgmnt = 0; sgmnt < params->numOverlapSegments; ++sgmnt )
      compute_data_segment( &segments->sgmnt[sgmnt], sgmnt, channel, invspec,
          response, params->segmentDuration, params->strideDuration, fwdplan );
  }
  else  /* only do the segments in the todo list */
  {
    UINT4 count;
    
    /* first count the number of segments to do */
    count = 0;
    for ( sgmnt = 0; sgmnt < params->numOverlapSegments; ++sgmnt )
      if ( is_in_list( sgmnt, params->segmentsToDoList ) )
        ++count;
      else
        continue; /* skip this segment: it is not in todo list */

    if ( ! count ) /* no segments to do */
      return NULL;

    segments->numSgmnt = count;
    segments->sgmnt = LALCalloc( segments->numSgmnt, sizeof(*segments->sgmnt) );
  
    count = 0;
    for ( sgmnt = 0; sgmnt < params->numOverlapSegments; ++sgmnt )
      if ( is_in_list( sgmnt, params->segmentsToDoList ) )
        compute_data_segment( &segments->sgmnt[count++], sgmnt, channel,
            invspec, response, params->segmentDuration, params->strideDuration,
            fwdplan );

    if ( params->writeSegment) /* write data segment */
            write_REAL4TimeSeries( channel );
    
  }

  return segments;
}


/* frees all memory that these routines have created */
static void ring_cleanup(
    ProcessParamsTable      *procpar,
    REAL4FFTPlan            *fwdplan,
    REAL4FFTPlan            *revplan,
    REAL4TimeSeries         *channel,
    COMPLEX8FrequencySeries *response,
    REAL4FrequencySeries    *invspec,
    RingTemplateBank        *bank,
    RingDataSegments        *segments,
    SnglRingdownTable          *events
    )
{
  while ( events )
  {
    SnglRingdownTable *thisEvent;
    thisEvent = events;
    events = events->next;
    LALFree( thisEvent );
  }
  if ( segments )
  {
    UINT4 sgmnt;
    for ( sgmnt = 0; sgmnt < segments->numSgmnt; ++sgmnt )
      XLALDestroyCOMPLEX8Vector( segments->sgmnt[sgmnt].data );
    LALFree( segments->sgmnt );
    LALFree( segments );
  }
  if ( bank )
  {
    XLALDestroyRingTemplateBank( bank );
  }
  if ( invspec )
  {
    XLALDestroyREAL4Vector( invspec->data );
    LALFree( invspec );
  }
  if ( response )
  {
    XLALDestroyCOMPLEX8Vector( response->data );
    LALFree( response );
  }
  if ( channel )
  {
    XLALDestroyREAL4Vector( channel->data );
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


/* routine to see if integer i is in a list of integers to do */
/* e.g., 2, 7, and 222 are in the list "1-3,5,7-" but 4 is not */
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
