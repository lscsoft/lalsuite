/*----------------------------------------------------------------------- 
 * 
 * File Name: inca.c
 *
 * Author: Brady, P. R.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataUtils.h>
#include <lalapps.h>
#include <processtable.h>

RCSID("$Id$");

#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "inca"

#define INCA_EARG   1
#define INCA_EROW   2
#define INCA_EFILE  3

#define INCA_MSGEARG   "Error parsing arguments"
#define INCA_MSGROW    "Error reading row from XML table"
#define INCA_MSGEFILE  "Could not open file"

#define MAXIFO 2
#define MSEC (1000000LL)

/* Usage format string. */
#define USAGE \
"Usage: %s [options] [LIGOLW XML input files]\n\n"\
"  --help                       display this message\n"\
"  --verbose                    print progress information\n"\
"  --debug-level LEVEL          set the LAL debug level to LEVEL\n"\
"  --user-tag STRING            set the process_params usertag to STRING\n"\
"  --comment STRING             set the process table comment to STRING\n"\
"\n"\
"  --gps-start-time SEC         GPS second of data start time\n"\
"  --gps-start-time-ns NS       GPS nanosecond of data start time\n"\
"\n"\
"  --ifo-a ifo_name             name of first ifo (e.g. L1, H1 or H2)\n"\
"  --ifo-b ifo_name             name of second ifo (e.g. L1, H1 or H2)\n"\
"\n"\
"  --drhoplus snr               positive signal to noise window (default 0)\n"\
"  --drhominus snr              negative signal to noise windoe (default 0)\n"\
"  --dm mass                    mass coincidence window (default 0)\n"\
"  --dt time                    time coincidence window (milliseconds)\n"\
"\n"\
"  --no-playground              do not select triggers from playground\n"\
"  --playground-only            only use triggers that are in playground\n"\
"  --write-uniq-triggers        make sure triggers from IFO A are unique\n" \
"\n"\
"[LIGOLW XML input files] list of the input trigger files.\n"\
"\n"

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
  calloc( 1, sizeof(ProcessParamsTable) ); \
  LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
      PROGRAM_NAME ); \
      LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
          long_options[option_index].name ); \
          LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
          LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

int main( int argc, char *argv[] )
{
  static LALStatus      status;
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;

  INT4  writeUniqTrigs = 0;
  INT4  usePlayground = 1;
  INT4  verbose = 0;
  INT4  startCoincidence = -1;
  INT4  endCoincidence = -1;
  CHAR  ifoName[MAXIFO][LIGOMETA_IFO_MAX];
  CHAR  comment[LIGOMETA_COMMENT_MAX];
  CHAR *userTag = NULL;

  CHAR  fileName[FILENAME_MAX];

  INT4  have_ifo_a_trigger = 0;
  INT4  isPlay = 0;
  INT8  ta, tb;

  SnglInspiralTable    *inspiralEventList[MAXIFO];
  SnglInspiralTable    *currentTrigger[MAXIFO];

  SnglInspiralTable    *currentEvent = NULL;
  SnglInspiralTable    *outEvent[MAXIFO];
  SnglInspiralTable    *coincidentEvents[MAXIFO];
  SnglInspiralAccuracy  errorParams;

  MetadataTable         proctable;
  MetadataTable         processParamsTable;
  MetadataTable         inspiralTable;
  ProcessParamsTable   *this_proc_param;
  LIGOLwXMLStream       xmlStream;
  
  INT4                  i, j;

  /* getopt arguments */
  struct option long_options[] =
  {
    {"verbose",                 no_argument,       &verbose,          1 },
    {"no-playground",           no_argument,       &usePlayground,    0 },
    {"playground-only",         no_argument,       &usePlayground,    1 },
    {"write-uniq-triggers",     no_argument,       &writeUniqTrigs,   1 },
    {"ifo-a",                   required_argument, 0,                'a'},
    {"ifo-b",                   required_argument, 0,                'b'},
    {"drhoplus",                required_argument, 0,                'c'},
    {"drhominus",               required_argument, 0,                'd'},
    {"dm",                      required_argument, 0,                'm'},
    {"dt",                      required_argument, 0,                't'},
    {"gps-start-time",          required_argument, 0,                'q'},
    {"gps-end-time",            required_argument, 0,                'r'},
    {"comment",                 required_argument, 0,                's'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    {"help",                    no_argument,       0,                'h'}, 
    {"debug-level",             required_argument, 0,                'z'},
    {0, 0, 0, 0}
  };
  int c;


  /*
   * 
   * initialize things
   *
   */


  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "1" );
  setvbuf( stdout, NULL, _IONBF, 0 );

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->start_time),
        &accuracy ), &status );
  LAL_CALL( populate_process_table( &status, proctable.processTable, 
        PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE ), &status );
  this_proc_param = processParamsTable.processParamsTable = 
    (ProcessParamsTable *) calloc( 1, sizeof(ProcessParamsTable) );
  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );

  memset( &errorParams, 0, sizeof(SnglInspiralAccuracy) );
  memset( inspiralEventList, 0, MAXIFO * sizeof(SnglInspiralTable *) );
  memset( currentTrigger, 0, MAXIFO * sizeof(SnglInspiralTable *) );
  memset( coincidentEvents, 0, MAXIFO * sizeof(SnglInspiralTable *) );
  memset( outEvent, 0, MAXIFO * sizeof(SnglInspiralTable *) );

  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    long int gpstime;
    size_t optarg_len;

    c = getopt_long_only( argc, argv, 
        "a:b:c:d:m:t:q:r:s:hz:Z:", long_options, &option_index );

    /* detect the end of the options */
    if ( c == -1 )
    {
      break;
    }
    
    switch ( c )
    {
      case 0:
        /* if this option set a flag, do nothing else now */
        if ( long_options[option_index].flag != 0 )
        {
          break;
        }
        else
        {
          fprintf( stderr, "Error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        break;

      case 'a':
        /* name of interferometer a */
        strncpy( ifoName[0], optarg, LIGOMETA_IFO_MAX * sizeof(CHAR) );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'b':
        /* name of interferometer b */
        strncpy( ifoName[1], optarg, LIGOMETA_IFO_MAX * sizeof(CHAR) );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'c':
        /* SNR error upward */
        errorParams.dRhoPlus = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'd':
        /* SNR error downward */
        errorParams.dRhoMinus = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'm':
        /* mass errors allowed */
        errorParams.dm = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 't':
        /* time coincidence window */
        errorParams.dtime = atof(optarg) * MSEC;
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'q':
        /* time coincidence window */
        gpstime = atol( optarg );
        if ( gpstime < 441417609 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is prior to " 
              "Jan 01, 1994  00:00:00 UTC:\n"
              "(%ld specified)\n",
              long_options[option_index].name, gpstime );
          exit( 1 );
        }
        if ( gpstime > 999999999 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is after " 
              "Sep 14, 2011  01:46:26 UTC:\n"
              "(%ld specified)\n", 
              long_options[option_index].name, gpstime );
          exit( 1 );
        }
        startCoincidence = (INT4) gpstime;
        ADD_PROCESS_PARAM( "int", "%ld", startCoincidence );
        break;

      case 'r':
        /* time coincidence window */
        gpstime = atol( optarg );
        if ( gpstime < 441417609 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is prior to " 
              "Jan 01, 1994  00:00:00 UTC:\n"
              "(%ld specified)\n",
              long_options[option_index].name, gpstime );
          exit( 1 );
        }
        if ( gpstime > 999999999 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is after " 
              "Sep 14, 2011  01:46:26 UTC:\n"
              "(%ld specified)\n", 
              long_options[option_index].name, gpstime );
          exit( 1 );
        }
        endCoincidence = (INT4) gpstime;
        ADD_PROCESS_PARAM( "int", "%ld", endCoincidence );
        break;

      case 's':
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

      case 'h':
        /* help message */
        fprintf( stderr, USAGE , argv[0]);
        exit( 1 );
        break;

      case 'z':
        set_debug_level( optarg );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'Z':
        /* create storage for the usertag */
        optarg_len = strlen(optarg) + 1;
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
        fprintf( stderr, USAGE , argv[0]);
        exit( 1 );
        break;

      default:
        fprintf( stderr, "Error: Unknown error while parsing options\n" );
        fprintf( stderr, USAGE, argv[0] );
        exit( 1 );
    }
  }
  
  /* check the values of the arguments */
  if ( startCoincidence < 0 )
  {
    fprintf( stderr, "Error: --gps-start-time must be specified\n" );
    exit( 1 );
  }

  if ( endCoincidence < 0 )
  {
    fprintf( stderr, "Error: --gps-end-time must be specified\n" );
    exit( 1 );
  }
    
  /* fill the comment, if a user has specified on, or leave it blank */
  if ( ! *comment )
  {
    LALSnprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
  } 
  else 
  {
    LALSnprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
  }

  /* store the write all trigs option */
  if ( writeUniqTrigs )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--write-uniq-triggers" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }

  /* store the playground argument in the process_params */
  LALSnprintf( processParamsTable.processParamsTable->program, 
      LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
  LALSnprintf( processParamsTable.processParamsTable->type, 
      LIGOMETA_TYPE_MAX, "string" );
  LALSnprintf( processParamsTable.processParamsTable->value, 
      LIGOMETA_TYPE_MAX, " " );
  if ( usePlayground )
  {
    LALSnprintf( processParamsTable.processParamsTable->param, 
        LIGOMETA_PARAM_MAX, "--playground-only" );
  }
  else
  {
    LALSnprintf( processParamsTable.processParamsTable->param, 
        LIGOMETA_PARAM_MAX, "--no-playground" );
  }

    

  /*
   *
   * read in the input data from the rest of the arguments
   *
   */


  if ( optind < argc )
  {
    for( i = optind; i < argc; ++i )
    {
      INT4 numTriggers = 0;
      SnglInspiralTable *inputData = NULL;

      if ( verbose ) 
        fprintf( stdout, "reading triggers from file: %s\n", argv[i] );

      numTriggers = 
        LALSnglInspiralTableFromLIGOLw( &inputData, argv[i], 0, -1 );

      if ( numTriggers < 0 )
      {
        fprintf( stderr, "error: unable to read sngl_inspiral table from %s\n", 
            argv[i] );
        exit( 1 );
      }
      else if ( numTriggers > 0 )
      {
        INT4 knownIFO = 0;

        if ( verbose ) 
          fprintf( stdout, "got %d sngl_inspiral rows from %s for ifo %s\n", 
            numTriggers, argv[i], inputData->ifo );

        /* locate the ifo associated with these triggers and store them */
        for ( j = 0; j < MAXIFO; ++j )
        {
          if ( ! strncmp( ifoName[j], inputData->ifo, LIGOMETA_IFO_MAX ) )
          {
            knownIFO = 1;

            if ( ! inspiralEventList[j] )
            {
              /* store the head of the linked list */
              inspiralEventList[j] = currentTrigger[j] = inputData;
            }
            else
            {
              /* append to the end of the linked list */
              currentTrigger[j]->next = inputData;
            }
            while ( currentTrigger[j]->next )
            {
              /* spin on to the end of the linked list */
              currentTrigger[j] = currentTrigger[j]->next;
            }

            if ( verbose ) fprintf( stdout, "added triggers to list\n" );
            break;
          }
        }

        if ( ! knownIFO )
        {
          /* catch an unknown ifo name among the input files */
          fprintf( stderr, "Error: unknown interferometer %s\n", 
              inputData->ifo );
          exit( 1 );
        }
      }
      else
      {
        if ( verbose ) 
          fprintf( stdout, "%s contains no triggers, skipping\n", argv[i] );
      }
    }
  }
  else
  {
    fprintf( stderr, "Error: No trigger files specified.\n" );
    exit( 1 );
  }


  for ( j = 0; j < MAXIFO; ++j )
  {
    if ( ! currentTrigger[j] )
    {
      fprintf( stdout, "No triggers read in for interferometer %d\n", j );
      goto cleanexit;
    }
  }


  /*
   *
   * sort the input data by time
   *
   */


  for ( j = 0; j < MAXIFO; ++j )
  {
    if ( verbose ) fprintf( stdout, "Sorting triggers from ifo %d\n", j );
    LAL_CALL( LALSortSnglInspiral( &status, &(inspiralEventList[j]),
          LALCompareSnglInspiralByTime ), &status );
  }


  /*
   * 
   * find the first trigger after coincidence start time for ifo A
   *
   */


  if ( verbose ) fprintf( stdout, "Moving to first trigger in window\n" );

  for ( j = 0; j < MAXIFO; ++j )
  {
    currentTrigger[j] = inspiralEventList[j];
  }

  while ( currentTrigger[0] && 
      ( currentTrigger[0]->end_time.gpsSeconds < startCoincidence ) )
  {
    currentTrigger[0] = currentTrigger[0]->next;
  }

  if ( ! currentTrigger[0] )
  {
    fprintf( stdout, "No triggers found in coincidence window\n" );
    goto cleanexit;
  }


  /*
   * 
   * outer loop over triggers from interferometer A
   *
   */


  if ( verbose ) fprintf( stdout, "start loop over ifo A\n" );

  while ( (currentTrigger[0] ) && 
      (currentTrigger[0]->end_time.gpsSeconds < endCoincidence) )
  {
    if ( verbose ) fprintf( stdout, "  using IFO A trigger at %d + %10.10f\n",
        currentTrigger[0]->end_time.gpsSeconds, 
        ((REAL4) currentTrigger[0]->end_time.gpsNanoSeconds * 1e-9) );

    LAL_CALL( LALGPStoINT8( &status, &ta, &(currentTrigger[0]->end_time) ), 
        &status );

    LAL_CALL( LALINT8NanoSecIsPlayground( &status, &isPlay, &ta ), &status );

    if ( verbose )
    {
      if ( isPlay )
      {
        fprintf( stdout, "  trigger is playground\n" );
      } 
      else
      {
        fprintf( stdout, "  trigger is not playground\n" );
      }
    }

    /* spin ifo b until the current trigger is within the coinicdence */
    /* window of the current ifo a trigger                            */
    while ( currentTrigger[1] )
    {
      LAL_CALL( LALGPStoINT8( &status, &tb, &(currentTrigger[1]->end_time) ), 
          &status );

      if ( tb > ta - errorParams.dtime )
      {
        /* we have reached the time coinicidence window */
        break;
      }

      currentTrigger[1] = currentTrigger[1]->next;
    }

    /* if we are playground only and the trigger is in playground or we are */
    /* not using playground and the trigger is not in the playground...     */
    if ( ( usePlayground && isPlay ) || ( ! usePlayground && ! isPlay) )
    {
      if ( verbose ) 
        fprintf( stdout, "  start loop over IFO B trigger at %d + %10.10f\n",
            currentTrigger[1]->end_time.gpsSeconds, 
            ((REAL4)currentTrigger[1]->end_time.gpsNanoSeconds * 1e-9) );

      /* look for coincident events in B within the time window */
      currentEvent = currentTrigger[1];

      while ( currentTrigger[1] )
      {
        LAL_CALL( LALGPStoINT8( &status, &tb, &(currentTrigger[1]->end_time) ), 
            &status );

        if (tb > ta + errorParams.dtime )
        {
          /* we are outside the time coincidence so move to the next event */
          break;
        }
        else
        {
          /* call the LAL function which compares events parameters */
          if ( verbose ) 
            fprintf( stdout, "    comparing IFO B trigger at %d + %10.10f\n",
                currentTrigger[1]->end_time.gpsSeconds, 
                ((REAL4)currentTrigger[1]->end_time.gpsNanoSeconds * 1e-9) );

          LAL_CALL( LALCompareSnglInspiral( &status, currentTrigger[0],
                currentTrigger[1], &errorParams ), &status );
        }

        if ( errorParams.match )
        {
          /* store this event for output */
          if ( verbose )
            fprintf( stdout, "    >>> found coincidence <<<\n" );

          for ( j = 0; j < MAXIFO; ++j )
          {
            /* only record the triggers from the primary ifo once */
            if ( ! writeUniqTrigs || j || ( ! j && ! have_ifo_a_trigger ) )
            {
              if ( ! coincidentEvents[j] )
              {
                coincidentEvents[j] = outEvent[j] = (SnglInspiralTable *) 
                  LALMalloc( sizeof(SnglInspiralTable) );
              }
              else
              {
                outEvent[j] = outEvent[j]->next = (SnglInspiralTable *) 
                  LALMalloc( sizeof(SnglInspiralTable) );
              }

              memcpy( outEvent[j], currentTrigger[j], 
                  sizeof(SnglInspiralTable) );
              outEvent[j]->next = NULL;

              have_ifo_a_trigger = 1;
            }
          }

        }

        currentTrigger[1] = currentTrigger[1]->next;

      } /* end loop over current events */

      /* go back to saved current IFO B trigger */
      currentTrigger[1] = currentEvent;
      have_ifo_a_trigger = 0;

    } /* end if ( IFO A is playground ) */

    /* go to the next ifo a trigger */
    currentTrigger[0] = currentTrigger[0]->next;

  } /* end loop over ifo A events */


  /*
   *
   * write the output xml file
   *
   */


cleanexit:

  if ( verbose ) fprintf( stdout, "writing output file... " );

  if ( userTag )
  {
    LALSnprintf( fileName, FILENAME_MAX, "%s-INCA_%s-%d-%d.xml", ifoName[0],
        userTag, startCoincidence, endCoincidence - startCoincidence );
  }
  else
  {
    LALSnprintf( fileName, FILENAME_MAX, "%s-INCA-%d-%d.xml", ifoName[0],
        startCoincidence, endCoincidence - startCoincidence );
  }

  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status , &xmlStream, fileName), &status );

  /* write process table */
  LALSnprintf( proctable.processTable->ifos, LIGOMETA_IFOS_MAX, "%s%s", 
      ifoName[0], ifoName[1] );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->end_time),
        &accuracy ), &status );
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, process_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, proctable, 
        process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );
  free( proctable.processTable );

  /* write process_params table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, process_params_table ),
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, processParamsTable, 
        process_params_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );

  /* write the sngl_inspiral table using events from IFO A */
  if ( coincidentEvents[0] )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status ,&xmlStream, sngl_inspiral_table),
        &status );
    inspiralTable.snglInspiralTable = coincidentEvents[0];
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, inspiralTable,
          sngl_inspiral_table), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream), &status );
  }

  LAL_CALL( LALCloseLIGOLwXMLFile( &status, &xmlStream), &status );

  /* write secondary ifo events to another file */
  if ( userTag )
  {
    LALSnprintf( fileName, FILENAME_MAX, "%s-INCA_%s-%d-%d.xml", ifoName[1],
        userTag, startCoincidence, endCoincidence - startCoincidence );
  }
  else
  {
    LALSnprintf( fileName, FILENAME_MAX, "%s-INCA-%d-%d.xml", ifoName[1],
        startCoincidence, endCoincidence - startCoincidence );
  }

  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status , &xmlStream, fileName), &status );
  
  /* write the sngl_inspiral table using events from IFO B */
  if ( coincidentEvents[1] )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status ,&xmlStream, sngl_inspiral_table),
        &status );
    inspiralTable.snglInspiralTable = coincidentEvents[1];
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, inspiralTable,
          sngl_inspiral_table), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream), &status );
  }

  LAL_CALL( LALCloseLIGOLwXMLFile( &status, &xmlStream), &status );

  if ( verbose ) fprintf( stdout, "done\n" );


  /*
   *
   * clean up the memory that has been allocated 
   *
   */


  if ( verbose ) fprintf( stdout, "freeing memory... " );

  while( processParamsTable.processParamsTable )
  {
    this_proc_param = processParamsTable.processParamsTable;
    processParamsTable.processParamsTable = this_proc_param->next;
    free( this_proc_param );
  }

  for( j = 0; j < MAXIFO ; ++j )
  {

    while ( coincidentEvents[j] )
    {
      currentEvent = coincidentEvents[j];
      coincidentEvents[j] = coincidentEvents[j]->next;
      LALFree( currentEvent );
    }

    while ( inspiralEventList[j] )
    {
      currentEvent = inspiralEventList[j];
      inspiralEventList[j] = inspiralEventList[j]->next;
      LALFree( currentEvent );
    }
  }

  if ( userTag ) free( userTag );

  if ( verbose ) fprintf( stdout, "done\n" );

  LALCheckMemoryLeaks();

  exit( 0 );
}
