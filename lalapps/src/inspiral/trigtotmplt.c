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
#include <lal/Date.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lalapps.h>
#include <processtable.h>

RCSID( "$Id$" );

#define PROGRAM_NAME "trig2tmplt"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"

#define USAGE \
"Usage: lalapps_print_inspiral [options]\n\n" \
"  --help                       display this message\n"\
"  --verbose                    print progress information\n"\
"  --debug-level LEVEL          set the LAL debug level to LEVEL\n"\
"\n"\
"  --minimal-match M            set output bank minimal match with M\n"\
"  --snr-threshold RHO          discard all triggers with snr less than RHO\n"\
"\n"\
"  --input FILE                 read triggers from FILE\n"\
"  --output FILE                write template bank to FILE\n"\
"\n"


/*
 *
 * compare function used by qsort
 *
 */


int compareTmpltsByMass ( const void *a, const void *b )
{
  SnglInspiralTable *aPtr = *((SnglInspiralTable **)a);
  SnglInspiralTable *bPtr = *((SnglInspiralTable **)b);

  if ( aPtr->mass1 > bPtr->mass1 )
  {
    return 1;
  }
  else if ( aPtr->mass1 < bPtr->mass1 )
  {
    return -1;
  }
  else if ( aPtr->mass2 > bPtr->mass2 )
  {
    return 1;
  }
  else if ( aPtr->mass2 < bPtr->mass2 )
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
  calloc( 1, sizeof(ProcessParamsTable) ); \
  LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
   PROGRAM_NAME ); \
   LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
     long_options[option_index].name ); \
     LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
     LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

int main ( int argc, char *argv[] )
{
  extern int vrbflg;
  LALStatus status = blank_status;
  LALLeapSecAccuracy accuracy = LALLEAPSEC_LOOSE;
  char *inputFileName = NULL;
  char *outputFileName = NULL;
  CHAR comment[LIGOMETA_COMMENT_MAX];
  REAL4 snr = 0;
  struct option long_options[] =
  {
    /* these options set a flag */
    {"help",                    no_argument,       0,                'h'},
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"debug-level",             required_argument, 0,                'z'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    {"comment",                 required_argument, 0,                'c'},
    {"minimal-match",           required_argument, 0,                'm'},
    {"snr-threshold",           required_argument, 0,                's'},
    {"input",                   required_argument, 0,                'i'},
    {"output",                  required_argument, 0,                'o'},
    {0, 0, 0, 0}
  };
  int i, numEvents, numUniq = 0;
  REAL4 minMatch = -1;
  CHAR                 *userTag = NULL;
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
    calloc( 1, sizeof(ProcessTable) );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->start_time),
        &accuracy ), &status );
  LAL_CALL( populate_process_table( &status, proctable.processTable, 
        PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE ), &status );
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
    calloc( 1, sizeof(ProcessParamsTable) );
  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );

  /* parse the command line arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    size_t optarg_len;

    i = getopt_long_only( argc, argv, 
        "c:i:o:m:z:Z:", long_options, &option_index );

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

      case 's':
        snr = (REAL4) atof( optarg );
        if ( snr <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "threshold must be > 0: "
              "(%f specified)\n",
              long_options[option_index].name, snr );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", snr );
        break;

      case 'i':
          /* create storage for the input file name name */
          optarg_len = strlen( optarg ) + 1;
          inputFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
          memcpy( inputFileName, optarg, optarg_len );
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'o':
        {
          /* create storage for the output file name name */
          optarg_len = strlen( optarg ) + 1;
          outputFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
          memcpy( outputFileName, optarg, optarg_len );
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        break;

      case 'm':
        minMatch = (REAL4) atof( optarg );
        if ( minMatch <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "minimal match of bank must be > 0: "
              "(%f specified)\n",
              long_options[option_index].name, minMatch );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", minMatch );
        break;

      case 'h':
        fprintf( stdout, USAGE );
        exit( 0 );
        break;

      case 'z':
        set_debug_level( optarg );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'Z':
        /* create storage for the usertag */
        optarg_len = strlen( optarg ) + 1;
        userTag = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( userTag, optarg, optarg_len );

        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", 
            PROGRAM_NAME );
        LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "-userTag" );
        LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
            optarg );
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

  /* can use LALCalloc() / LALMalloc() from here */

  /* fill the comment, if a user has specified it, or leave it blank */
  if ( ! *comment )
  {
    LALSnprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
  } else {
    LALSnprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
  }
  
  /* check that the input and output file names and match have been specified */
  if ( ! inputFileName )
  {
    fprintf( stderr, "--input-file must be specified\n" );
    exit( 1 );
  }
  if ( ! outputFileName )
  {
    fprintf( stderr, "--output-file must be specified\n" );
    exit( 1 );
  }
  if ( minMatch < 0 )
  {
    fprintf( stderr, "--minimal-match must be specified\n" );
    exit( 1 );
  }


  /*
   *
   * read the entire sngl_inspiral table from the input file
   *
   */

  
  if ( vrbflg ) 
    fprintf( stdout, "reading sngl_inspiral table from %s\n", inputFileName );

  numEvents = SnglInspiralTableFromLIGOLw( &eventHead, inputFileName, 0, -1 );

  if ( numEvents < 0 )
  {
    fprintf( stderr, "error: unable to read sngl_inspiral table from %s\n", 
        inputFileName );
    exit( 1 );
  }
  else if ( numEvents == 0 )
  {
    goto cleanexit;
  }

  if ( vrbflg ) fprintf( stdout, "parsed %d sngl_inspiral rows from %s\n", 
      numEvents, inputFileName );


  /*
   *
   * discard template that are less than the snr threshold
   *
   */


  if ( eventHead )
  {
    SnglInspiralTable *tmpEventHead;
    prevEvent = tmpEventHead = NULL;
    if ( snr )
    {
      int numSnrEvents = 0;
      if ( vrbflg ) fprintf( stdout, "discarding events with snr < %e... ",
          snr );
      thisEvent = eventHead;
      while ( thisEvent )
      {
        if ( thisEvent->snr >= snr )
        {
          ++numSnrEvents;
          if ( ! tmpEventHead )
          {
            prevEvent = tmpEventHead = thisEvent;
          }
          else
          {
            prevEvent = prevEvent->next = thisEvent;
          }
          thisEvent = thisEvent->next;
        }
        else
        {
          SnglInspiralTable *discard = thisEvent;
          thisEvent = thisEvent->next;
          LALFree( discard );
        }
      }
      if ( (eventHead = tmpEventHead) ) prevEvent->next = NULL;
      numEvents = numSnrEvents;
      if ( vrbflg ) fprintf( stdout, "done\n" );
    }
  }
  
  if ( ! numEvents ) goto cleanexit;

  
  /*
   *
   * order the templates and discard duplicates
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
  qsort( eventHandle, numEvents, sizeof(eventHandle[0]), compareTmpltsByMass );
  if ( vrbflg ) fprintf( stdout, "done\n" );

  /* create a linked list of sorted templates */
  prevEvent = eventHead = eventHandle[0];
  numUniq = 1;
  for ( i = 1; i < numEvents; ++i )
  {
    if ( (prevEvent->mass1 == eventHandle[i]->mass1)  &&
        (prevEvent->mass2 == eventHandle[i]->mass2) ) 
    {
      /* discard the event as it is a duplicate */
      LALFree( eventHandle[i] );
    }
    else
    {
      /* add the event to the linked list */
      prevEvent = prevEvent->next = eventHandle[i];
      ++numUniq;
    }
  }
  prevEvent->next = NULL;
  

  /*
   *
   * write the sorted, unique list out to an XML file
   *
   */

cleanexit:

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
  free( proctable.processTable );

  /* erase the first empty process params entry */
  {
    ProcessParamsTable *emptyPPtable = procparams.processParamsTable;
    procparams.processParamsTable = procparams.processParamsTable->next;
    free( emptyPPtable );
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
    free( this_proc_param );
  }
  
  if ( eventHead )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, sngl_inspiral_table ), 
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, outputTable, 
          sngl_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );

    LALFree( eventHandle );
    while ( eventHead )
    {
      thisEvent = eventHead;
      eventHead = eventHead->next;
      LALFree( thisEvent );
    }
  }

  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &results ), &status );

  if ( vrbflg ) fprintf( stdout, "wrote %d sngl_inspiral rows to %s\n", 
      numUniq, outputFileName );


  /*
   *
   * free memory and exit
   *
   */


  if ( userTag ) free ( userTag );
  free( inputFileName );
  free( outputFileName );
  LALCheckMemoryLeaks();
  
  return 0;
}
