/*----------------------------------------------------------------------- 
 * 
 * File Name: inca.c
 *
 * Author: Brady, P. R., Brown, D. A. and Fairhurst, S.
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
#include <sys/types.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataUtils.h>
#include <lalapps.h>
#include <processtable.h>

RCSID("$Id$");

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
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

/* Usage format string. */
#define USAGE \
  "Usage: %s [options] [LIGOLW XML input files]\n\n"\
"  --help                    display this message\n"\
"  --verbose                 print progress information\n"\
"  --version                 print version information and exit\n"\
"  --debug-level LEVEL       set the LAL debug level to LEVEL\n"\
"  --user-tag STRING         set the process_params usertag to STRING\n"\
"  --ifo-tag STRING          set the ifo-tag to STRING - for file naming\n"\
"  --comment STRING          set the process table comment to STRING\n"\
"\n"\
"  --ligo-coinc              use the LIGO (inca) coincidence routine\n"\
"  --tama-coinc CLUSTER      use the TAMA coincidence routine, with CLUSTER\n"\
"                            one of [ snr_and_chisq | snrsq_over_chisq | snr ]\n"\
"  --gps-start-time SEC      GPS second of data start time\n"\
"  --gps-end-time SEC        GPS second of data end time\n"\
"\n"\
"  --silde-time SEC          slide all triggers of IFOB by SEC\n"\
"  --slide-time-ns NS        slide all triggers of IFOB by NS\n"\
"\n"\
"  --ifo-a IFOA              name of first ifo (e.g. L1, H1 or H2)\n"\
"  --ifo-b IFOB              name of second ifo (e.g. L1, H1 or H2)\n"\
"\n"\
"  --parameter-test TEST    set the desired parameters to test coincidence\n"\
"                           (m1_and_m2|psi0_and_psi3|mchirp_and_eta)\n"\
"  --dm Dm                   mass coincidence window (default 0)\n"\
"  --dpsi0 Dpsi0             psi0 coincidence window\n"\
"  --dpsi3 Dpsi3             psi3 coincidence window\n"\
"  --dmchirp Dmchirp         mchirp coincidence window\n"\
"  --deta  Deta              eta coincidence window\n"\
"  --dt Dt                   time coincidence window (milliseconds)\n"\
"\n"\
"  --no-playground           do not select triggers from playground\n"\
"  --playground-only         only use triggers that are in playground\n"\
"  --all-data                use all triggers\n"\
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

enum { undefined, ligo_coinc, tama_coinc } coincTest = undefined;

int main( int argc, char *argv[] )
{
  static LALStatus      status;
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;

  extern int vrbflg;

  LALPlaygroundDataMask dataType;

  INT4  startCoincidence = -1;
  INT4  endCoincidence = -1;
  LIGOTimeGPS startCoinc = {0,0};
  LIGOTimeGPS endCoinc = {0,0};
  CHAR  ifoName[MAXIFO][LIGOMETA_IFO_MAX];
  CHAR  comment[LIGOMETA_COMMENT_MAX];
  CHAR *userTag = NULL;
  CHAR *ifoTag = NULL;

  CHAR  fileName[FILENAME_MAX];
  CHAR *xmlFileName;

  LIGOTimeGPS slideData = {0,0};
  INT8  slideDataNS = 0;
  INT4  numIFO;
  INT4  numEvents = 0;

  SnglInspiralTable    *inspiralEventList[MAXIFO];
  SnglInspiralTable    *currentTrigger[MAXIFO];

  SnglInspiralTable    *currentEvent = NULL;
  SnglInspiralTable    *coincidentEvents[MAXIFO];
  SnglInspiralAccuracy  errorParams;

  SearchSummvarsTable  *inputFiles = NULL;
  SearchSummvarsTable  *thisInputFile = NULL;

  SnglInspiralClusterChoice   clusterchoice = none;
  
  MetadataTable         proctable;
  MetadataTable         processParamsTable;
  MetadataTable         searchsumm;
  MetadataTable         searchSummvarsTable;
  MetadataTable         inspiralTable;
  ProcessParamsTable   *this_proc_param = NULL;
  LIGOLwXMLStream       xmlStream;

  INT4                  i, j;

  /* getopt arguments */
  struct option long_options[] =
  {
    {"verbose",                 no_argument,     &vrbflg,             1 },
    {"ligo-coinc",              no_argument,     &coincTest, ligo_coinc },
    {"no-playground",           no_argument,       0,                'Q'},
    {"playground-only",         no_argument,       0,                'R'},
    {"all-data",                no_argument,       0,                'D'},
    {"tama-coinc",              required_argument, 0,                'T'},
    {"ifo-a",                   required_argument, 0,                'a'},
    {"ifo-b",                   required_argument, 0,                'b'},
    {"dm",                      required_argument, 0,                'm'},
    {"parameter-test",          required_argument, 0,                'A'},
    {"dpsi0",                   required_argument, 0,                'p'},
    {"dpsi3",                   required_argument, 0,                'P'},
    {"dmchirp",                 required_argument, 0,                'c'},
    {"deta",                    required_argument, 0,                'n'},
    {"dt",                      required_argument, 0,                't'},
    {"gps-start-time",          required_argument, 0,                'q'},
    {"gps-end-time",            required_argument, 0,                'r'},
    {"comment",                 required_argument, 0,                's'},
    {"slide-time",              required_argument, 0,                'X'},
    {"slide-time-ns",           required_argument, 0,                'Y'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    {"ifo-tag",                 required_argument, 0,                'I'},
    {"help",                    no_argument,       0,                'h'}, 
    {"debug-level",             required_argument, 0,                'z'},
    {"version",                 no_argument,       0,                'V'},
    {0, 0, 0, 0}
  };
  int c;
  INT4 haveTest = 0;


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

  /* create the search summary and zero out the summvars table */
  searchsumm.searchSummaryTable = (SearchSummaryTable *)
    calloc( 1, sizeof(SearchSummaryTable) );

  memset( &errorParams, 0, sizeof(SnglInspiralAccuracy) );
  memset( inspiralEventList, 0, MAXIFO * sizeof(SnglInspiralTable *) );
  memset( currentTrigger, 0, MAXIFO * sizeof(SnglInspiralTable *) );
  memset( coincidentEvents, 0, MAXIFO * sizeof(SnglInspiralTable *) );


  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    long int gpstime;
    size_t optarg_len;

    c = getopt_long_only( argc, argv, 
        "a:b:e:k:A:m:p:P:t:q:r:s:hz:I:Z:M:T:S:c:n:QRDT:", long_options, 
        &option_index );

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
        LALSnprintf( ifoName[0], LIGOMETA_IFO_MAX, "%s", optarg );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'b':
        /* name of interferometer b */
        LALSnprintf( ifoName[1], LIGOMETA_IFO_MAX, "%s", optarg );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;


      case 'A':
        if ( ! strcmp( "m1_and_m2", optarg ) )
        {
          errorParams.test = m1_and_m2;
        }
        else if ( ! strcmp( "psi0_and_psi3", optarg ) )
        {
          errorParams.test = psi0_and_psi3;
        }
        else if ( ! strcmp( "mchirp_and_eta", optarg ) )
        {
          errorParams.test = mchirp_and_eta;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown test specified: "
              "%s (must be m1_and_m2, psi0_and_psi3 or mchirp_and_eta)\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        haveTest = 1;
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;


      case 'm':
        /* mass errors allowed */
        errorParams.dm = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'p':
        /* psi0 errors allowed */
        errorParams.dpsi0 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'P':
        /* psi3 errors allowed */
        errorParams.dpsi3 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'c':
        /* mass errors allowed */
        errorParams.dmchirp = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'n':
        /* mass errors allowed */
        errorParams.deta = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 't':
        /* time coincidence window, argument is in milliseconds */
        errorParams.dt = atof(optarg) * 1000000LL;
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
        startCoinc.gpsSeconds = startCoincidence;
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
        endCoinc.gpsSeconds = endCoincidence;
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

      case 'Q':
        dataType = exclude_play;
        break;

      case 'R':
        dataType = playground_only;
        break;

      case 'D':
        dataType = all_data;
        break;

      case 'T':
        /* use the TAMA coincidence routine */
        coincTest = tama_coinc;
        
        /* set the clustering algorithm */
        {        
          if ( ! strcmp( "snr_and_chisq", optarg ) )
          {
            clusterchoice = snr_and_chisq;
          }
          else if ( ! strcmp( "snrsq_over_chisq", optarg) )
          {
            clusterchoice = snrsq_over_chisq;
          }
          else if ( ! strcmp( "snr", optarg) )
          {
            clusterchoice = snr;
          }        
          else
          {
            fprintf( stderr, "invalid argument to  --%s:\n"
                "unknown clustering specified:\n "
                "%s (must be one of: snr_and_chisq, \n"
                "   snrsq_over_chisq or snr)\n",
                long_options[option_index].name, optarg);
            exit( 1 );
          }
          ADD_PROCESS_PARAM( "string", "%s", optarg );
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

      case 'I':
        /* create storage for the ifo-tag */
        optarg_len = strlen(optarg) + 1;
        ifoTag = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( ifoTag, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'X':
        slideData.gpsSeconds = (INT4) atoi( optarg );
        ADD_PROCESS_PARAM( "int", "%d", slideData.gpsSeconds );
        break;

      case 'Y':
        slideData.gpsNanoSeconds = (INT4) atoi( optarg );
        ADD_PROCESS_PARAM( "int", "%d", slideData.gpsNanoSeconds );
        break;

      case 'V':
        /* print version information and exit */
        fprintf( stdout, "Inspiral Coincidence and Triggered Bank Generator\n" 
            "Patrick Brady, Duncan Brown and Steve Fairhurst\n"
            "CVS Version: " CVS_ID_STRING "\n"
            "CVS Tag: " CVS_NAME_STRING "\n" );
        exit( 0 );
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

  if ( ! haveTest )
  {
    fprintf( stderr, "--parameter-test must be specified\n" );
    exit( 1 );
  }

  /* decide how many ifos we have based on what we are doing */
  numIFO = 2;

  /* calculate the slide time in nanoseconds */
  LAL_CALL( LALGPStoINT8( &status, &slideDataNS, &slideData ), &status );

  if (  slideDataNS )
  {
    /* check that a dataType option is not specified if */
    /* doing a slide                        */
    if ( dataType )
    {
      fprintf( stderr, "Should not specify --playground-only, --no-playground"
          "or --all-data for a time slide\n" );
      exit( 1 );
    }
    else
    {
      dataType = all_data;
    }
  }
    
  /* fill the comment, if a user has specified one, or leave it blank */
  if ( ! *comment )
  {
    LALSnprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
    LALSnprintf( searchsumm.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX, 
        " " );
  } 
  else 
  {
    LALSnprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
    LALSnprintf( searchsumm.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
  }

  if ( slideDataNS )
  {
    /* erase the first empty process params */
    ProcessParamsTable *tmpProc = processParamsTable.processParamsTable;
    processParamsTable.processParamsTable = 
      processParamsTable.processParamsTable->next;
    free( tmpProc );
  }
  else
  {
    /* store the playground argument in the process_params */
    LALSnprintf( processParamsTable.processParamsTable->program, 
        LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
    LALSnprintf( processParamsTable.processParamsTable->type, 
        LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( processParamsTable.processParamsTable->value, 
        LIGOMETA_TYPE_MAX, " " );
    if ( dataType == playground_only )
    {
      LALSnprintf( processParamsTable.processParamsTable->param, 
          LIGOMETA_PARAM_MAX, "--playground-only" );
    }
    else if ( dataType == exclude_play  )
    {
      LALSnprintf( processParamsTable.processParamsTable->param, 
          LIGOMETA_PARAM_MAX, "--no-playground" );
    }
    else
    {
      LALSnprintf( processParamsTable.processParamsTable->param, 
          LIGOMETA_PARAM_MAX, "--all-data" );
    }       
  }


  if ( coincTest == undefined )
  {
    fprintf( stderr, "Must specify either the LIGO or TAMA coincidence\n"
        "--ligo-coinc or --tama-coinc\n");
    exit( 1 );
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
      INT4 numFileTriggers = 0;
      SnglInspiralTable *inputData = NULL;     

      if ( vrbflg ) fprintf( stdout, 
          "storing input file name %s in search summvars table\n", argv[i] );

      if ( ! inputFiles )
      {
        inputFiles = thisInputFile = (SearchSummvarsTable *)
          LALCalloc( 1, sizeof(SearchSummvarsTable) );
      }
      else
      {
        thisInputFile = thisInputFile->next = (SearchSummvarsTable *)
          LALCalloc( 1, sizeof(SearchSummvarsTable) );
      }
      LALSnprintf( thisInputFile->name, LIGOMETA_NAME_MAX, 
          "input_file" );
      LALSnprintf( thisInputFile->string, LIGOMETA_NAME_MAX, 
          "%s", argv[i] );      


      /* read in triggers */

      if ( vrbflg ) 
        fprintf( stdout, "reading triggers from file: %s\n", argv[i] );

      numFileTriggers = 
        LALSnglInspiralTableFromLIGOLw( &inputData, argv[i], 0, -1 );

      if ( numFileTriggers < 0 )
      {
        fprintf( stderr, "error: unable to read sngl_inspiral table from %s\n", 
            argv[i] );
        exit( 1 );
      }
      else if ( numFileTriggers > 0 )
      {
        INT4 knownIFO = 0;

        if ( vrbflg ) 
          fprintf( stdout, "got %d sngl_inspiral rows from %s for ifo %s\n", 
              numFileTriggers, argv[i], inputData->ifo );

        /* locate the ifo associated with these triggers and store them */
        for ( j = 0; j < numIFO ; ++j )
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
              /* append to the end of the linked list and set current    */
              /* trigger to the first trigger of the list being appended */
              currentTrigger[j] = currentTrigger[j]->next = inputData;
            }

            if ( slideDataNS && j == 1 && vrbflg)  fprintf( stdout, 
                "Doing a time slide of %d sec %d nanosec on IFOB triggers\n",
                slideData.gpsSeconds, slideData.gpsNanoSeconds );       

            while ( currentTrigger[j]->next )
            {
              /* spin on to the end of the linked list */
              /* doing time slides if necessary */
              if ( slideDataNS && j == 1 )
              {
                INT8 trigTimeNS = 0;
                LAL_CALL( LALGPStoINT8( &status, &trigTimeNS, 
                      &(currentTrigger[j]->end_time) ), &status );
                trigTimeNS += slideDataNS;
                LAL_CALL( LALINT8toGPS( &status, 
                      &(currentTrigger[j]->end_time), &trigTimeNS ), 
                    &status );
              }     
              currentTrigger[j] = currentTrigger[j]->next;
            }

            /* slide the last trigger */
            if ( slideDataNS && j == 1)
            {
              INT8 trigTimeNS = 0;
              LAL_CALL( LALGPStoINT8( &status, &trigTimeNS, 
                    &(currentTrigger[j]->end_time) ), &status );
              trigTimeNS += slideDataNS;
              LAL_CALL( LALINT8toGPS( &status, &(currentTrigger[j]->end_time), 
                    &trigTimeNS ), &status );
            }

            if ( vrbflg ) fprintf( stdout, "added triggers to list\n" );
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
        if ( vrbflg ) 
          fprintf( stdout, "%s contains no triggers, skipping\n", argv[i] );
      }
    }
  }
  else
  {
    fprintf( stderr, "Error: No trigger files specified.\n" );
    exit( 1 );
  }


  for ( j = 0; j < numIFO; ++j )
  {
    if ( ! inspiralEventList[j] )
    {
      /* no triggers in this ifo so no coincidences can be found */

      fprintf( stdout, "No triggers read in for interferometer %d\n", j );
      goto cleanexit;
    }
  }


  /*
   *
   * Preprocess the input data
   *
   */
  for ( j = 0; j < numIFO; ++j )
  {
    /* time sort the triggers */
    if ( vrbflg ) fprintf( stdout, "Sorting triggers from ifo %d\n", j );
    LAL_CALL( LALSortSnglInspiral( &status, &(inspiralEventList[j]),
          LALCompareSnglInspiralByTime ), &status );

    /* keep only triggers in coinc window */
    if ( vrbflg ) fprintf( stdout, 
        "Discarding triggers outside requested interval\n" );
    LAL_CALL( LALTimeCutSingleInspiral( &status, &inspiralEventList[j],
          &startCoinc, &endCoinc), &status );

    /* keep play/non-play/all triggers */
    if ( dataType == playground_only && vrbflg ) fprintf( stdout, 
        "Keeping only playground triggers\n" );
    else if ( dataType == exclude_play && vrbflg ) fprintf( stdout, 
        "Keeping only non-playground triggers\n" );
    else if ( dataType == all_data && vrbflg ) fprintf( stdout, 
        "Keeping all triggers\n" );
    LAL_CALL( LALPlayTestSingleInspiral( &status, &inspiralEventList[j],
          &dataType ), &status );
  }


  /* Test for coincidence */
  if ( coincTest == ligo_coinc )
  {
    LALIncaCoincidenceTest( &status, &coincidentEvents[0], 
        &coincidentEvents[1], inspiralEventList[0], inspiralEventList[1], 
        &errorParams ); 
  }
  else if (coincTest == tama_coinc )
  {
    LALTamaCoincidenceTest( &status, &coincidentEvents[1], 
        &coincidentEvents[0], inspiralEventList[1], inspiralEventList[0], 
        &errorParams, clusterchoice );
  }

  /* count number of coincident events */
  for( currentEvent = coincidentEvents[0], numEvents = 0; currentEvent;
      currentEvent = currentEvent->next, ++numEvents);
  
cleanexit:

  /* search summary entries: nevents is from primary ifo */
  searchsumm.searchSummaryTable->in_start_time = startCoinc;
  searchsumm.searchSummaryTable->in_end_time = endCoinc;
  searchsumm.searchSummaryTable->out_start_time = startCoinc;
  searchsumm.searchSummaryTable->out_end_time = endCoinc;
  searchsumm.searchSummaryTable->nevents = numEvents;
  searchsumm.searchSummaryTable->nnodes = 1;

  if ( vrbflg ) fprintf( stdout, "writing output file... " );

  for ( j = 0; j < numIFO; ++j )
  {
    /* set the file name correctly */
    if ( userTag && ifoTag )
    {
      LALSnprintf( fileName, FILENAME_MAX, "%s-INCA_%s_%s-%d-%d.xml", 
          ifoName[j], ifoTag, userTag, startCoincidence, 
          endCoincidence - startCoincidence );
    }
    else if ( userTag && !ifoTag )
    {
      LALSnprintf( fileName, FILENAME_MAX, "%s-INCA_%s-%d-%d.xml", 
          ifoName[j], userTag, startCoincidence, 
          endCoincidence - startCoincidence );
    }
    else if ( !userTag && ifoTag )
    {
      LALSnprintf( fileName, FILENAME_MAX, "%s-INCA_%s-%d-%d.xml", 
          ifoName[j], ifoTag, startCoincidence, 
          endCoincidence - startCoincidence );
    }
    else
    {
      LALSnprintf( fileName, FILENAME_MAX, "%s-INCA-%d-%d.xml", ifoName[j],
          startCoincidence, endCoincidence - startCoincidence );
    }

    xmlFileName = fileName;

    memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
    LAL_CALL( LALOpenLIGOLwXMLFile( &status , &xmlStream, xmlFileName), 
        &status );

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

    /* write process_params table */
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, 
          process_params_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, processParamsTable, 
          process_params_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );

    /* write search_summary table */
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, 
          search_summary_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, searchsumm, 
          search_summary_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );

    /* write the search_summvars tabls */
    LAL_CALL( LALBeginLIGOLwXMLTable( &status ,&xmlStream, 
          search_summvars_table), &status );
    searchSummvarsTable.searchSummvarsTable = inputFiles;
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, searchSummvarsTable,
          search_summvars_table), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream), &status );

    /* write the sngl_inspiral table using events from ifoName[j] */
    if ( coincidentEvents[j] )
    {
      LAL_CALL( LALBeginLIGOLwXMLTable( &status ,&xmlStream, 
            sngl_inspiral_table), &status );
      inspiralTable.snglInspiralTable = coincidentEvents[j];
      LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, inspiralTable,
            sngl_inspiral_table), &status );
      LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream), &status );
    }

    LAL_CALL( LALCloseLIGOLwXMLFile( &status, &xmlStream), &status );
  }

  if ( vrbflg ) fprintf( stdout, "done\n" );


  /*
   *
   * clean up the memory that has been allocated 
   *
   */


  if ( vrbflg ) fprintf( stdout, "freeing memory... " );

  free( proctable.processTable );
  free( searchsumm.searchSummaryTable );

  while ( processParamsTable.processParamsTable )
  {
    this_proc_param = processParamsTable.processParamsTable;
    processParamsTable.processParamsTable = this_proc_param->next;
    free( this_proc_param );
  }

  while ( inputFiles )
  {
    thisInputFile = inputFiles;
    inputFiles = thisInputFile->next;
    LALFree( thisInputFile );
  }

  for( j = 0; j < numIFO; ++j )
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
  if ( ifoTag ) free( ifoTag );

  if ( vrbflg ) fprintf( stdout, "done\n" );

  LALCheckMemoryLeaks();

  exit( 0 );
}
