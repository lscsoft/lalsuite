#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/Date.h>
#include <lal/LIGOLwXML.h>
#include <lal/PrintFTSeries.h>

#include "processtable.h"
#include "lalapps.h"
#include "errutil.h"
#include "gpstime.h"
#include "ring.h"

RCSID( "$Id$" );


/*
 *
 * Routines to output event triggers.
 *
 */


/* macro is_option() (and friends): determine if string is an option */
/* i.e., does it start with "-[a-zA-Z]" or "--[a-zA-Z]" */
#define is_long_option(s) \
  ( strlen(s) > 2 && (s)[0] == '-' && (s)[1] == '-' && isalpha( s[2] ) )
#define is_short_option(s) \
  ( strlen(s) > 1 && (s)[0] == '-' && isalpha( s[1] ) )
#define is_option(s) ( is_long_option(s) || is_short_option(s) )

/* routine to output events as ascii */
static int ring_output_events_asc( 
    SnglBurstTable     *events,
    struct ring_params *params
    );

/* routine to output events as LIGOLw XML */
static int ring_output_events_xml( 
    SnglBurstTable     *events,
    ProcessParamsTable *processParamsTable,
    struct ring_params *params
    );


/* creates a process table */
static ProcessTable * ring_create_process_table( struct ring_params *params );


/* creates a search-summary table */
static SearchSummaryTable *ring_create_search_summary( struct ring_params *params );


/* creates a process params table from command line arguments */
ProcessParamsTable * create_process_params( int argc, char **argv,
    const char *program )
{
  ProcessParamsTable *processParamsTable = NULL;
  int c;
  for ( c = 1; c < argc; ++c )
  {
    char *arg = argv[c];
    char *opt = NULL;
    char *val = NULL;

    /* duplicate arg in opt */
    opt = LALMalloc( strlen( arg ) + 1 );
    strcpy( opt, arg );

    if ( ! is_option( opt ) )
      error( "%s is not an option\n", opt );

    /* search for equals-sign to identify val */
    if ( ( val = strchr( opt, '=' ) ) )
      *val++ = 0; /* nul terminate opt and set val to value */
    else if ( c < argc - 1 && ! is_option( argv[c+1] ) ) /* next arg is value */
      val = argv[++c]; /* set value and increment counter */
    else /* no value for this option */
      val = NULL;

    /* now write the option and value */
    if ( opt )
    {
      ProcessParamsTable *thisParam;
      thisParam = LALCalloc( 1, sizeof( *thisParam ) );
      thisParam->next = processParamsTable;
      processParamsTable = thisParam;

      strncpy( thisParam->program, program, LIGOMETA_PROGRAM_MAX - 1 );
      strncpy( thisParam->param, opt, LIGOMETA_PARAM_MAX - 1 );
      if ( val )
      {
        strncpy( thisParam->type, "string", LIGOMETA_TYPE_MAX - 1 );
        strncpy( thisParam->value, val, LIGOMETA_VALUE_MAX - 1 );
      }
    }

    LALFree( opt );
  }

  return processParamsTable;
}


/* routine to output events */
int ring_output_events( 
    SnglBurstTable     *events,
    ProcessParamsTable *processParamsTable,
    struct ring_params *params
    )
{
  if ( ! strlen( params->outputFile ) )
    return 0;
  if ( params->outputFormat == output_ligolw )
    ring_output_events_xml( events, processParamsTable, params );
  if ( params->outputFormat == output_ascii )
    ring_output_events_asc( events, params );
  return 0;
}


/* routine to output events as an ascii file */
static int ring_output_events_asc( 
    SnglBurstTable     *events,
    struct ring_params *params
    )
{
  SnglBurstTable *thisEvent;
  FILE *fp;
  verbose( "output events to ascii file %s\n", params->outputFile );
  fp = fopen( params->outputFile, "w" );
  fprintf( fp, "# gps start time\tsignal/noise\tamplitude\tfrequency\tbandwidth\n" );
  thisEvent = events;
  while ( thisEvent )
  {
    fprintf( fp, "%9d.%09d\t%e\t%e\t%e\t%e\n",
        (int) thisEvent->start_time.gpsSeconds,
        (int) thisEvent->start_time.gpsNanoSeconds,
        thisEvent->snr,
        thisEvent->amplitude,
        thisEvent->central_freq,
        thisEvent->bandwidth );
    thisEvent = thisEvent->next;
  }
  fclose( fp );
  return 0;
}


/* routine to output events as LIGOLw XML file */
static int ring_output_events_xml( 
    SnglBurstTable     *events,
    ProcessParamsTable *processParamsTable,
    struct ring_params *params
    )
{
  LALStatus status = blank_status;
  MetadataTable   process;
  MetadataTable   processParams;
  MetadataTable   searchSummary;
  MetadataTable   ringEvents;
  LIGOLwXMLStream results;

  verbose( "output events to LIGOLw XML file %s\n", params->outputFile );

  memset( &process, 0, sizeof( process ) );
  memset( &processParams, 0, sizeof( processParams ) );
  memset( &searchSummary, 0, sizeof( searchSummary ) );
  memset( &ringEvents, 0, sizeof( ringEvents ) );
  memset( &results, 0, sizeof( results ) );

  /* create process table and search summary tables */
  process.processTable = ring_create_process_table( params );
  processParams.processParamsTable = processParamsTable;
  searchSummary.searchSummaryTable = ring_create_search_summary( params );
  ringEvents.snglBurstTable = events;

  /* open results xml file */
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &results, params->outputFile ), &status );

  /* output the process table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, process, process_table ), &status );

  LAL_CALL( LALEndLIGOLwXMLTable( &status, &results ), &status );

  /* output process params table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_params_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, processParams, process_params_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable( &status, &results ), &status );

  /* output search summary table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, search_summary_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, searchSummary, search_summary_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable( &status, &results ), &status );

  /* output the events */
  if ( ringEvents.snglBurstTable )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, sngl_burst_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, ringEvents, sngl_burst_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &results ), &status );
  }

  /* close the xml file */
  LAL_CALL( LALCloseLIGOLwXMLFile( &status, &results ), &status );

  LALFree( searchSummary.searchSummaryTable );
  LALFree( process.processTable );

  return 0;
}


/* routine to create process table */
ProcessTable *ring_create_process_table( struct ring_params *params )
{
  LALStatus status = blank_status;
  LALLeapSecAccuracy accuracy = LALLEAPSEC_LOOSE;
  ProcessTable *processTable = NULL;

  processTable = LALCalloc( 1, sizeof( *processTable ) );

  /* call lalapps routine to populate the process table */
  LAL_CALL( populate_process_table( &status, processTable, params->programName,
        params->cvsRevision, params->cvsSource, params->cvsDate ), &status );
  strncpy( processTable->comment, " ", LIGOMETA_COMMENT_MAX );
  strncpy( processTable->ifos, params->ifoName, LIGOMETA_IFOS_MAX );
  LAL_CALL( LALGPSTimeNow( &status, &processTable->end_time, &accuracy ), &status );

  return processTable;
}


/* routine to create search summary table */
/* FIXME: need to handle trigstarttime and trigendtime */
static SearchSummaryTable *ring_create_search_summary( struct ring_params *params )
{
  SearchSummaryTable *searchSummary = NULL;
  LIGOTimeGPS outStartTime;
  LIGOTimeGPS outEndTime;
  INT8 outStartTimeNS;
  INT8 outEndTimeNS;

  /* setup search summary table */
  searchSummary = LALCalloc( 1, sizeof( *searchSummary ) );
  strncpy( searchSummary->comment, params->programName, LIGOMETA_COMMENT_MAX );
  searchSummary->nnodes = 1;

  /* compute the start and end times of data analyzed */
  outStartTimeNS  = epoch_to_ns( &params->startTime ) 
    + sec_to_ns( 0.5 * params->strideDuration );
  outEndTimeNS    = outStartTimeNS 
    + sec_to_ns( params->strideDuration * params->numOverlapSegments );
  ns_to_epoch( &outStartTime, outStartTimeNS );
  ns_to_epoch( &outEndTime, outEndTimeNS );

  /* store input start time and end time of raw data in search summary */
  searchSummary->in_start_time  = params->startTime;
  searchSummary->in_end_time    = params->endTime;
  searchSummary->out_start_time = outStartTime;
  searchSummary->out_end_time   = outEndTime;
  searchSummary->nevents        = params->numEvents;

  return searchSummary;
}


/*
 *
 * Routines to write intermediate results (time/frequency series and bank).
 *
 */


/* routine to construct an appropriately-formatted filename from series name */
static int generate_file_name( char *fname, size_t size,
    const char *sname, int t, int dt );
#define FILENAME_SIZE 256


/* routine to write a time series */
int write_REAL4TimeSeries( REAL4TimeSeries *series )
{
  char fname[FILENAME_SIZE];
  int t, dt;
  t  = series->epoch.gpsSeconds;
  dt = ceil(series->epoch.gpsNanoSeconds + series->data->length*series->deltaT);
  generate_file_name( fname, sizeof( fname ), series->name, t, dt );
  verbose( "writing series %s to file %s\n", series->name, fname );
  LALSPrintTimeSeries( series, fname );
  return 0;
}


/* routine to write a real frequency series */
int write_REAL4FrequencySeries( REAL4FrequencySeries *series )
{
  char fname[FILENAME_SIZE];
  int t, dt;
  t  = series->epoch.gpsSeconds;
  dt = ceil(series->epoch.gpsNanoSeconds + 1.0/series->deltaF);
  generate_file_name( fname, sizeof( fname ), series->name, t, dt );
  verbose( "writing series %s to file %s\n", series->name, fname );
  LALSPrintFrequencySeries( series, fname );
  return 0;
}


/* routine to write a complex frequency series */
int write_COMPLEX8FrequencySeries( COMPLEX8FrequencySeries *series )
{
  char fname[FILENAME_SIZE];
  int t, dt;
  t  = series->epoch.gpsSeconds;
  dt = ceil(series->epoch.gpsNanoSeconds + 1.0/series->deltaF);
  generate_file_name( fname, sizeof( fname ), series->name, t, dt );
  verbose( "writing series %s to file %s\n", series->name, fname );
  LALCPrintFrequencySeries( series, fname );
  return 0;
}


/* routine to write a ringdown template bank */
int write_bank( RingTemplateBank *bank )
{
  const char *fname = "RING_BANK.dat";
  FILE *fp;
  UINT4 tmplt;
  verbose( "writing template bank to file %s\n", fname );
  fp = fopen( fname, "w" );
  fprintf( fp, "# template bank\n" );
  fprintf( fp, "# template phase=%e\n", bank->tmplt->phase );
  fprintf( fp, "# frequency (Hz)\tquality\n" );
  for ( tmplt = 0; tmplt < bank->numTmplt; ++tmplt )
    fprintf( fp, "%e\t%e\n", bank->tmplt[tmplt].frequency,
        bank->tmplt[tmplt].quality );
  fclose( fp );
  return 0;
}


/* routine to construct an appropriately-formatted filename from series name */
static int generate_file_name( char *fname, size_t size,
    const char *sname, int t, int dt )
{
  char *c;

  strncpy( fname, sname, size - 1 );
  fname[size-1] = 0;

  /* slashes are not allowed */
  if ( strchr( fname, '/' ) )
    error( "slashes are not allowed in output file name %s\n", fname );

  /* convert hyphens to underscores */
  while ( ( c = strchr( fname, '-' ) ) )
    *c = '_';

  /* convert colons to hypens */
  while ( ( c = strchr( fname, ':' ) ) )
    *c = '-';

  /* convert spaces to underscores */
  while ( ( c = strchr( fname, ' ' ) ) )
    *c = '_';

  LALSnprintf( fname, size, "%s-%d-%d.dat", fname, t, dt );

  return 0;
}
