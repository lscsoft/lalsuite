/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Lisa M. Goggin
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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

#include "config.h"
#include "coh_PTF.h"
#include <lal/LIGOLwXML.h>
#include "LIGOLwXMLlegacy.h"

/*
 *
 * Routines to output event triggers.
 * Note that a lot of functions should be merged with ring_output
 *
 */


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
      strncpy( thisParam->type, "string", LIGOMETA_TYPE_MAX - 1 );
      if ( val )
        strncpy( thisParam->value, val, LIGOMETA_VALUE_MAX - 1 );
    }

    LALFree( opt );
  }

  return processParamsTable;
}

static int XLALWriteLIGOLwXMLTimeSlideSegmentMapTable(
	LIGOLwXMLStream *xml,
	const TimeSlideSegmentMapTable *time_slide_seg_map
)
{
	const char *row_head = "\n\t\t\t";

	/* table header */

	XLALClearErrno();
	XLALFilePuts("\t<Table Name=\"time_slide_segment_map:table\">\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"time_slide_segment_map:segment_def_id\" Type=\"ilwd:char\"/>\n", xml->fp);
	XLALFilePuts("\t\t<Column Name=\"time_slide_segment_map:time_slide_id\" Type=\"ilwd:char\"/>\n", xml->fp);
        XLALFilePuts("\t\t<Stream Name=\"time_slide_segment_map:table\" Type=\"Local\" Delimiter=\",\">", xml->fp);
	if(XLALGetBaseErrno())
		XLAL_ERROR(XLAL_EFUNC);

	/* rows */

	for(; time_slide_seg_map; time_slide_seg_map = time_slide_seg_map->next) {
		if(XLALFilePrintf(xml->fp, "%s\"segment_def:segment_def_id:%ld\",\"time_slide:time_slide_id:%ld\"",
			row_head,
			time_slide_seg_map->segment_def_id,
			time_slide_seg_map->time_slide_id
		) < 0)
			XLAL_ERROR(XLAL_EFUNC);
		row_head = ",\n\t\t\t";
	}

	/* table footer */

	if(XLALFilePuts("\n\t\t</Stream>\n\t</Table>\n", xml->fp) < 0)
		XLAL_ERROR(XLAL_EFUNC);

	/* done */

	return 0;
}


/* routine to output events as LIGOLw XML file */
int coh_PTF_output_events_xml(
    char               *outputFile,
    MultiInspiralTable  *events,
    SnglInspiralTable *snglEvents,
    SimInspiralTable *injections,
    ProcessParamsTable *processParamsTable,
    TimeSlide          *time_slide_head,
    TimeSlideSegmentMapTable *time_slide_map_head,
    SegmentTable       *segment_table_head,
    struct coh_PTF_params *params
    )
{
  LALStatus XLAL_INIT_DECL(status);
  LIGOLwXMLStream *results;

  verbose( "output events to LIGOLw XML file %s\n", outputFile );

  /* open results xml file */
  results = XLALOpenLIGOLwXMLFile(outputFile);

  /* output the process table */
  ProcessTable* processTable = coh_PTF_create_process_table(params);
  XLALWriteLIGOLwXMLProcessTable(results,processTable);
  LALFree(processTable);

  /* output process params table */
  XLALWriteLIGOLwXMLProcessParamsTable(results, processParamsTable);

  /* output search summary table */
  SearchSummaryTable* searchSummTable = coh_PTF_create_search_summary(params);
  XLALWriteLIGOLwXMLSearchSummaryTable(results,searchSummTable);
  LALFree(searchSummTable);

  /* write the signals injected in a template bank simulation */
  if ( injections )
    XLALWriteLIGOLwXMLSimInspiralTable( results, injections );

  /* output time slide table */
  XLALWriteLIGOLwXMLTimeSlideTable( results, time_slide_head);

  /* output time slide map table */
  XLALWriteLIGOLwXMLTimeSlideSegmentMapTable( results, time_slide_map_head);

  /* output segment list */
  XLALWriteLIGOLwXMLSegmentTable( results, segment_table_head);

  /* output the events */
  if (! params->writeSnglInspiralTable)
  {
    MetadataTable   ringEvents;
    memset( &ringEvents, 0, sizeof( ringEvents ) );
    ringEvents.multiInspiralTable = events;
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, results, multi_inspiral_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, results, ringEvents,multi_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, results ), &status );
  }
  else
  {
    (void) events;	/* silence unused variable warning */
    XLALWriteLIGOLwXMLSnglInspiralTable( results, snglEvents );
  }

  /* close the xml file */
  XLALCloseLIGOLwXMLFile(results);

  return 0;
}

/* routine to output template bank as LIGOLw XML file */
int coh_PTF_output_tmpltbank(
    char               *outputFile,
    SnglInspiralTable  *tmplts,
    ProcessParamsTable *processParamsTable,
    struct coh_PTF_params *params
    )
{
  ProcessTable *process;
  ProcessParamsTable *processParams;
  SearchSummaryTable *searchSummary;
  LIGOLwXMLStream *results;

  verbose( "output template bank to LIGOLw XML file %s\n", outputFile );

  /* create process table and search summary tables */
  process = coh_PTF_create_process_table( params );
  processParams = processParamsTable;
  searchSummary = coh_PTF_create_search_summary( params );

  /* open results xml file */
  results = XLALOpenLIGOLwXMLFile( outputFile );

  /* output the process table */
  XLALWriteLIGOLwXMLProcessTable( results, process );

  /* output process params table */
  XLALWriteLIGOLwXMLProcessParamsTable( results, processParams );

  /* output search summary table */
  XLALWriteLIGOLwXMLSearchSummaryTable( results, searchSummary );

  /* output the events */
  if ( tmplts )
    XLALWriteLIGOLwXMLSnglInspiralTable( results, tmplts );

  /* close the xml file */
  XLALCloseLIGOLwXMLFile( results );

  XLALDestroyProcessTable( process );
  XLALDestroyProcessParamsTable( processParams );
  XLALDestroySearchSummaryTable( searchSummary );

  return 0;
}



/* routine to create process table */
ProcessTable *coh_PTF_create_process_table( struct coh_PTF_params *params )
{
  ProcessTable *processTable = NULL;

  processTable = LALCalloc( 1, sizeof( *processTable ) );

  /* call lalapps routine to populate the process table */

  XLALPopulateProcessTable(processTable, params->programName,
      lalAppsVCSIdentInfo.vcsId,lalAppsVCSIdentInfo.vcsStatus,lalAppsVCSIdentInfo.vcsDate,0);

  strncpy( processTable->comment, " ", LIGOMETA_COMMENT_MAX );

  /* store ifos */
  if ( params->numIFO == 1 )
  {
    snprintf( processTable->ifos, LIGOMETA_IFOS_MAX,\
              "%s", params->ifoName[0] );
  }
  else if( params->numIFO == 2 )
  {
    XLALStringPrint( processTable->ifos, LIGOMETA_IFOS_MAX,\
              "%s%s", params->ifoName[0], params->ifoName[1] );
  }
  else if ( params->numIFO == 3 )
  {
    XLALStringPrint( processTable->ifos, LIGOMETA_IFOS_MAX,\
              "%s%s%s", params->ifoName[0], params->ifoName[1],
        params->ifoName[2] );
  }
  else if ( params->numIFO == 4 )
  {
    XLALStringPrint( processTable->ifos, LIGOMETA_IFOS_MAX,\
              "%s%s%s%s", params->ifoName[0], params->ifoName[1],
              params->ifoName[2], params->ifoName[3]);
  }

  processTable->start_time = params->jobStartTime;
  XLALGPSTimeNow(&processTable->end_time);

  return processTable;
}


/* routine to create search summary table */
SearchSummaryTable *coh_PTF_create_search_summary( struct coh_PTF_params *params )
{
  SearchSummaryTable *searchSummary = NULL;
  LIGOTimeGPS outStartTime;
  LIGOTimeGPS outEndTime;
  INT8 outStartTimeNS;
  INT8 outEndTimeNS;

  /* setup search summary table */
  searchSummary = LALCalloc( 1, sizeof( *searchSummary ) );
  strncpy( searchSummary->comment, params->programName, LIGOMETA_COMMENT_MAX-1 );
  searchSummary->nnodes = 1;

  /* compute the start and end times of data analyzed */
  outStartTimeNS  = epoch_to_ns( &params->startTime )
    + sec_to_ns( params->analStartPoint / (REAL4)params->sampleRate );
  outEndTimeNS    = epoch_to_ns( &params->endTime )
    - sec_to_ns( (params->numTimePoints - params->analEndPoint)
                 / (REAL4)params->sampleRate );
  if( params->trigStartTimeNS && (params->trigStartTimeNS > outStartTimeNS))
    outStartTimeNS = (outEndTimeNS < params->trigStartTimeNS) ?
      outEndTimeNS : params->trigStartTimeNS;
  if( params->trigEndTimeNS && (params->trigEndTimeNS < outEndTimeNS))
    outEndTimeNS = (outStartTimeNS > params->trigEndTimeNS) ?
      outStartTimeNS : params->trigEndTimeNS;
  ns_to_epoch( &outStartTime, outStartTimeNS );
  ns_to_epoch( &outEndTime, outEndTimeNS );

  /* store input start time and end time of raw data in search summary */
  searchSummary->in_start_time  = params->startTime;
  searchSummary->in_end_time    = params->endTime;
  searchSummary->out_start_time = outStartTime;
  /*XLALGPSAdd( &searchSummary->out_start_time, -1.0 * params->padData ); */
  searchSummary->out_end_time   = outEndTime;
  /*XLALGPSAdd( &searchSummary->out_end_time, 1.0 * params->padData); */
  searchSummary->nevents        = params->numEvents;

  /* store ifos */
  if ( params->numIFO == 1 )
  {
    snprintf( searchSummary->ifos, LIGOMETA_IFOS_MAX,\
              "%s", params->ifoName[0] );
  }
  else if( params->numIFO == 2 )
  {
    XLALStringPrint( searchSummary->ifos, LIGOMETA_IFOS_MAX,\
              "%s%s", params->ifoName[0], params->ifoName[1] );
  }
  else if ( params->numIFO == 3 )
  {
    XLALStringPrint( searchSummary->ifos, LIGOMETA_IFOS_MAX,\
              "%s%s%s", params->ifoName[0], params->ifoName[1],
        params->ifoName[2] );
  }
  else if ( params->numIFO == 4 )
  {
    XLALStringPrint( searchSummary->ifos, LIGOMETA_IFOS_MAX,\
              "%s%s%s%s", params->ifoName[0], params->ifoName[1],
              params->ifoName[2], params->ifoName[3]);
  }

  return searchSummary;
}


/*
 *
 * Routines to write intermediate results (time/frequency series and bank).
 *
 */

/* routine to write a time series */
int write_REAL4TimeSeries( REAL4TimeSeries *series )
{
  char fname[FILENAME_SIZE];
  int t, dt;
  t  = series->epoch.gpsSeconds;
  dt = ceil(1e-9*series->epoch.gpsNanoSeconds + series->data->length*series->deltaT);
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
  dt = ceil(1e-9*series->epoch.gpsNanoSeconds + 1.0/series->deltaF);
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
  dt = ceil(1e-9*series->epoch.gpsNanoSeconds + 1.0/series->deltaF);
  generate_file_name( fname, sizeof( fname ), series->name, t, dt );
  verbose( "writing series %s to file %s\n", series->name, fname );
  LALCPrintFrequencySeries( series, fname );
  return 0;
}

/* routine to construct an appropriately-formatted filename from series name */
int generate_file_name( char *fname, size_t size,
    const char *sname, int t, int dt )
{
  char *c;
  char *tmp_name;

  tmp_name = (char *) LALCalloc( size, sizeof(char) );

  strncpy( tmp_name, sname, size - 1 );
  tmp_name[size-1] = 0;

  /* slashes are not allowed */
  if ( strchr( tmp_name, '/' ) )
    error( "slashes are not allowed in output file name %s\n", tmp_name );

  /* convert hyphens to underscores */
  while ( ( c = strchr( tmp_name, '-' ) ) )
    *c = '_';

  /* convert colons to hypens */
  while ( ( c = strchr( tmp_name, ':' ) ) )
    *c = '-';

  /* convert spaces to underscores */
  while ( ( c = strchr( tmp_name, ' ' ) ) )
    *c = '_';

  snprintf( fname, size, "%s-%d-%d.dat", tmp_name, t, dt );

  LALFree( tmp_name );

  return 0;
}
