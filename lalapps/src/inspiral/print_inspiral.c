/*----------------------------------------------------------------------- 
 * 
 * File Name: trig2tmplt.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <regex.h>
#include <time.h>
#include <glob.h>
#include <lal/Date.h>
#include <lal/LIGOLwXML.h>
#include <lalapps.h>
#include <processtable.h>
#include "ligolwbank.h"

RCSID( "$Id$" );

#define PROGRAM_NAME "trig2tmplt"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"


/*
 *
 * compare function used by qsort
 *
 */


int compareTmpltsByTime ( const void *a, const void *b )
{
  SnglInspiralTable *aPtr = *((SnglInspiralTable **)a);
  SnglInspiralTable *bPtr = *((SnglInspiralTable **)b);
  INT8 ta, tb;
  LALStatus status = blank_status;

  LAL_CALL( LALGPStoINT8( &status, &ta, &(aPtr->end_time) ), &status );
  LAL_CALL( LALGPStoINT8( &status, &tb, &(bPtr->end_time) ), &status );

  if ( ta > tb )
  {
    return 1;
  }
  else if ( ta < tb )
  {
    return -1;
  }
  else
  {
    return 0;
  }
}


#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
  LALCalloc( 1, sizeof(ProcessParamsTable) ); \
  LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
   PROGRAM_NAME ); \
   LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
     long_options[option_index].name ); \
     LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
     LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

int main ( int argc, char *argv[] )
{
  extern int vrbflg;
  int playground = 0;
  LALStatus status = blank_status;
  LALLeapSecAccuracy accuracy = LALLEAPSEC_LOOSE;
  char *inputGlob = NULL;
  char *outputFileName = NULL;
  glob_t globbedFiles;
  CHAR comment[LIGOMETA_COMMENT_MAX];
  struct option long_options[] =
  {
    /* these options set a flag */
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"playground",              no_argument,       &playground,       1 },
    {"comment",                 required_argument, 0,                'c'},
    {"input",                   required_argument, 0,                'i'},
    {"output",                  required_argument, 0,                'o'},
    {"debug-level",             required_argument, 0,                'z'},
    {0, 0, 0, 0}
  };
  int i, numEvents = 0;
  size_t j;
  SnglInspiralTable    *eventHead = NULL;
  SnglInspiralTable    *thisEvent = NULL;
  SnglInspiralTable    *prevEvent = NULL;
  SnglInspiralTable   **eventHandle = NULL;
  ProcessParamsTable   *this_proc_param;
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         outputTable;
  LIGOLwXMLStream       results;
  

  /*
   *
   * initialization
   *
   */


  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "1" );

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) 
    LALCalloc( 1, sizeof(ProcessTable) );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->start_time),
        &accuracy ), &status );
  LAL_CALL( populate_process_table( &status, proctable.processTable, 
        PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE ), &status );
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
    LALCalloc( 1, sizeof(ProcessParamsTable) );
  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );

  /* parse the command line arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    size_t namelen;

    i = getopt_long_only( argc, argv, 
        "c:i:o:m:z:", long_options, &option_index );

    /* detect the end of the options */
    if ( i == - 1 )
    {
      break;
    }

    switch ( i )
    {
      case 0:
        /* if this option set a flag, do nothing else now */
        if ( long_options[option_index].flag != 0 )
        {
          break;
        }
        else
        {
          fprintf( stderr, "error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        break;

      case 'c':
        if ( strlen( optarg ) > LIGOMETA_COMMENT_MAX - 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "comment must be less than %d characters\n",
              long_options[option_index].name, LIGOMETA_COMMENT_MAX );
          exit( 1 );
        }
        else
        {
          LALSnprintf( comment, LIGOMETA_COMMENT_MAX, "%s", optarg);
        }
        break;

      case 'i':
        {
          /* create storage for the input file name name */
          namelen = strlen( optarg ) + 1;
          inputGlob = (CHAR *) LALCalloc( namelen, sizeof(CHAR));
          memcpy( inputGlob, optarg, namelen );
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        break;

      case 'o':
        {
          /* create storage for the output file name name */
          namelen = strlen( optarg ) + 1;
          outputFileName = (CHAR *) LALCalloc( namelen, sizeof(CHAR));
          memcpy( outputFileName, optarg, namelen );
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        break;

      case 'z':
        set_debug_level( optarg );
        break;

      case '?':
        exit( 1 );
        break;

      default:
        fprintf( stderr, "unknown error while parsing options\n" );
        exit( 1 );
    }
  }

  if ( optind < argc )
  {
    fprintf( stderr, "extraneous command line arguments:\n" );
    while ( optind < argc )
    {
      fprintf ( stderr, "%s\n", argv[optind++] );
    }
    exit( 1 );
  }

  /* fill the comment, if a user has specified it, or leave it blank */
  if ( ! *comment )
  {
    LALSnprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
  } else {
    LALSnprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
  }
  
  /* check that the input and output file names and match have been specified */
  if ( ! inputGlob )
  {
    fprintf( stderr, "--input must be specified\n" );
    exit( 1 );
  }
  if ( ! outputFileName )
  {
    fprintf( stderr, "--output must be specified\n" );
    exit( 1 );
  }


  /*
   *
   * find all the input files and read them in
   *
   */


  if ( glob( inputGlob, GLOB_ERR, NULL, &globbedFiles ) )
  {
    fprintf( stderr, "error globbing files from %s\n", inputGlob );
    perror( "error:" );
    exit( 1 );
  }

  if ( vrbflg )
  {
    fprintf( stdout, "globbed input files:\n" );
    for ( j = 0; j < globbedFiles.gl_pathc; ++j )
    {
      fprintf( stdout, "%s\n", globbedFiles.gl_pathv[j] );
    }
  }

  eventHandle = &eventHead;

  for ( j = 0; j < globbedFiles.gl_pathc; ++j )
  {
    int thisNumEvents;
    char *inputFileName = globbedFiles.gl_pathv[j];

    if ( vrbflg ) 
      fprintf( stdout, "reading sngl_inspiral table from %s\n", inputFileName );

    thisNumEvents = 
      SnglInspiralTableFromLIGOLw( eventHandle, inputFileName, 0, -1 );

    if ( thisNumEvents < 0 )
    {
      fprintf( stderr, "error: unable to read sngl_inspiral table from %s\n", 
          inputFileName );
      exit( 1 );
    }

    if ( vrbflg ) fprintf( stdout, "parsed %d sngl_inspiral rows from %s\n", 
        thisNumEvents, inputFileName );

    if ( thisNumEvents )
    {
      numEvents += thisNumEvents;

      thisEvent = *eventHandle;
      while ( thisEvent->next )
      {
        thisEvent = thisEvent->next;
      }
      eventHandle = &(thisEvent->next);
    }
  }

  if ( vrbflg ) fprintf( stdout, "got %d sngl_inspiral events\n", numEvents );

  globfree( &globbedFiles );

  
  /*
   *
   * order the templates by time
   *
   */


  /* create an array of ptrs, qsort it and turn it back into a linked list */
  eventHandle = (SnglInspiralTable **) 
    LALCalloc( numEvents, sizeof(SnglInspiralTable *) );

  for ( i = 0, thisEvent = eventHead; i < numEvents; 
      ++i, thisEvent = thisEvent->next )
  {
    eventHandle[i] = thisEvent;
  }

  if ( vrbflg ) fprintf( stdout, "sorting events... " );
  qsort( eventHandle, numEvents, sizeof(eventHandle[0]), compareTmpltsByTime );
  if ( vrbflg ) fprintf( stdout, "done\n" );

  /* create a linked list of sorted templates */
  if ( vrbflg ) 
  {
    fprintf( stdout, "reordering linked list " );
  }
  if ( vrbflg && playground )
  {
    fprintf( stdout, "for playground events only... " );
  }
  else if ( vrbflg )
  {
    fprintf( stdout, "... " );
  }
  
  if ( playground )
  {
    /* get only the playground events */
    int numPlayEvents = 0;
    eventHead = NULL;
    for ( i = 0; i < numEvents; ++i )
    {
      if ( ((eventHandle[i]->end_time.gpsSeconds - 729273613) % 6370) < 600 )
      {
        ++numPlayEvents;
        if ( ! eventHead )
        {
          prevEvent = eventHead = eventHandle[i];
        }
        else
        {
          prevEvent = prevEvent->next = eventHandle[i];
        }
      }
    }
    numEvents = numPlayEvents;
  }
  else
  {
    /* add all events to the linked list */
    prevEvent = eventHead = eventHandle[0];
    for ( i = 1; i < numEvents; ++i )
    {
      prevEvent = prevEvent->next = eventHandle[i];
    }
  }

  if ( vrbflg ) fprintf( stdout, "done\n" );

  if ( eventHead )
  {
    prevEvent->next = NULL;
  }
  else
  {
    if ( vrbflg ) fprintf( stdout, "no events left\n" );
  }
  

  /*
   *
   * write the sorted, unique list out to an XML file
   *
   */


  if ( vrbflg ) fprintf( stdout, "writing output file... " );

  memset( &results, 0, sizeof(LIGOLwXMLStream) );
  outputTable.snglInspiralTable = eventHead;
  
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &results, outputFileName),
      &status );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->end_time),
        &accuracy ), &status );
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, proctable, 
        process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  LALFree( proctable.processTable );

  /* erase the first empty process params entry */
  {
    ProcessParamsTable *emptyPPtable = procparams.processParamsTable;
    procparams.processParamsTable = procparams.processParamsTable->next;
    LALFree( emptyPPtable );
  }

  /* write the process params table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_params_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, procparams, 
        process_params_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  while( procparams.processParamsTable )
  {
    this_proc_param = procparams.processParamsTable;
    procparams.processParamsTable = this_proc_param->next;
    LALFree( this_proc_param );
  }
  
  if ( eventHead )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, sngl_inspiral_table ), 
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, outputTable, 
          sngl_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  }

  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &results ), &status );

  if ( vrbflg ) fprintf( stdout, "wrote %d sngl_inspiral rows to %s\n", 
      numEvents, outputFileName );


  /*
   *
   * free memory and exit
   *
   */


  if ( vrbflg ) fprintf( stdout, "freeing memory... ");
  fflush( stdout );
  if ( eventHead )
  {
    LALFree( eventHandle );
    while ( eventHead )
    {
      thisEvent = eventHead;
      eventHead = eventHead->next;
      LALFree( thisEvent );
    }
  }
  LALFree( inputGlob );
  LALFree( outputFileName );
  LALCheckMemoryLeaks();
  if ( vrbflg ) fprintf( stdout, "done\n");
  
  return 0;
}
