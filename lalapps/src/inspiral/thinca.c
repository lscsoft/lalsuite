/*----------------------------------------------------------------------- 
 * 
 * File Name: thinca.c
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
#include <sys/stat.h>
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

#define MAXIFO 4

/* Usage format string. */
#define USAGE \
  "Usage: %s [options] [LIGOLW XML input files]\n\n"\
"  --help                    display this message\n"\
"  --verbose                 print progress information\n"\
"  --version                 print version information and exit\n"\
"  --debug-level LEVEL       set the LAL debug level to LEVEL\n"\
"  --user-tag STRING         set the process_params usertag to STRING\n"\
"  --comment STRING          set the process table comment to STRING\n"\
"\n"\
"  --gps-start-time SEC      GPS second of data start time\n"\
"  --gps-end-time SEC        GPS second of data end time\n"\
"\n"\
"  --g1-triggers             input triggers from G1\n"\
"  --h1-triggers             input triggers from H1\n"\
"  --h2-triggers             input triggers from H2\n"\
"  --l1-triggers             input triggers from L1\n"\
"  --t1-triggers             input triggers from T1\n"\
"  --v1-triggers             input triggers from V1\n"\
"\n"\
"  --parameter-test TEST    set the desired parameters to test coincidence\n"\
"                           for inca: (m1_and_m2|psi0_and_psi3|mchirp_and_eta)\n"\
"  --dm Dm                   mass coincidence window (default 0)\n"\
"  --dpsi0 Dpsi0             psi0 coincidence window\n"\
"  --dpsi3 Dpsi3             psi3 coincidence window\n"\
"  --dmchirp Dmchirp         mchirp coincidence window\n"\
"  --deta  Deta              eta coincidence window\n"\
"  --dt Dt                   time coincidence window (milliseconds)\n"\
"\n"\
"  --data-type DATA_TYPE     specify the data type, must be one of\n"\
"                            (playground_only|exclude_play|all_data)\n"\
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

int g1_trig = 0;
int h1_trig = 0;
int h2_trig = 0;
int l1_trig = 0;
int t1_trig = 0;
int v1_trig = 0;


int main( int argc, char *argv[] )
{
  static LALStatus      status;
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;

  extern int vrbflg;

  LALPlaygroundDataMask dataType;
  INT4  startCoincidence = -1;
  LIGOTimeGPS startCoinc = {0,0};
  INT4  endCoincidence = -1;
  LIGOTimeGPS endCoinc = {0,0};
  CHAR  ifoName[MAXIFO][LIGOMETA_IFO_MAX];
  CHAR  comment[LIGOMETA_COMMENT_MAX];
  CHAR *userTag = NULL;

  CHAR  fileName[FILENAME_MAX];

  INT8  currentTriggerNS[MAXIFO];
  INT4  haveCoinc[MAXIFO];
  INT4  numIFO;
  INT4  numTriggers = 0;
  INT4  inStartTime = -1;
  INT4  inEndTime = -1;

  INT8  maxTC = 0;

  SnglInspiralTable    *inspiralEventList = NULL;
  SnglInspiralTable    *thisInspiralTrigger = NULL;
  SnglInspiralTable    *currentTrigger[MAXIFO];
  CoincInspiralTable   *currentCoincCand[MAXIFO];
  INT4  haveTest = 0;
  INT4  haveDataType = 0;

  CoincInspiralTable   *coincInspiralList = NULL;
  CoincInspiralTable   *prevCoincInspiral = NULL;

  SnglInspiralAccuracy  errorParams;

  SearchSummvarsTable  *inputFiles = NULL;
  SearchSummvarsTable  *thisInputFile = NULL;

  SearchSummaryTable   *searchSummList = NULL;
  SearchSummaryTable   *thisSearchSumm = NULL;

  SummValueTable       *summValueList = NULL;
  SummValueTable       *thisSummValue = NULL;

  MetadataTable         proctable;
  MetadataTable         processParamsTable;
  MetadataTable         searchsumm;
  MetadataTable   searchSummvarsTable;
  MetadataTable   summValueTable;
  MetadataTable         inspiralTable;
  ProcessParamsTable   *this_proc_param = NULL;
  LIGOLwXMLStream       xmlStream;

  INT4                  i;

  /* getopt arguments */
  struct option long_options[] =
  {
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"g1-triggers",             no_argument,       &g1_trig,          1 },
    {"h1-triggers",             no_argument,       &h1_trig,          1 },
    {"h2-triggers",             no_argument,       &h2_trig,          1 },
    {"l1-triggers",             no_argument,       &l1_trig,          1 },
    {"t1-triggers",             no_argument,       &t1_trig,          1 },
    {"v1-triggers",             no_argument,       &v1_trig,          1 },
    {"parameter-test",          required_argument, 0,                'A'},
    {"dm",                      required_argument, 0,                'm'},
    {"dpsi0",                   required_argument, 0,                'p'},
    {"dpsi3",                   required_argument, 0,                'P'},
    {"dmchirp",                 required_argument, 0,                'c'},
    {"deta",                    required_argument, 0,                'n'},
    {"dt",                      required_argument, 0,                't'},
    {"gps-start-time",          required_argument, 0,                'q'},
    {"gps-end-time",            required_argument, 0,                'r'},
    {"data-type",               required_argument, 0,                'D'},
    {"comment",                 required_argument, 0,                's'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    {"help",                    no_argument,       0,                'h'}, 
    {"debug-level",             required_argument, 0,                'z'},
    {"version",                 no_argument,       0,                'V'},
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

  /* create the search summary and zero out the summvars table */
  searchsumm.searchSummaryTable = (SearchSummaryTable *)
    calloc( 1, sizeof(SearchSummaryTable) );

  memset( &errorParams, 0, sizeof(SnglInspiralAccuracy) );
  memset( currentTrigger, 0, MAXIFO * sizeof(SnglInspiralTable *) );
  memset( currentCoincCand, 0, MAXIFO * sizeof(CoincInspiralTable *) );
  memset( currentTriggerNS, 0, MAXIFO * sizeof(INT8) );
  memset( haveCoinc, 0, MAXIFO * sizeof(INT4) );

  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    long int gpstime;
    size_t optarg_len;

    c = getopt_long_only( argc, argv, 
        "c:hm:n:p:q:r:s:t:z:A:D:P:V:Z:", long_options, 
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

      case 'A':
        /* set the parameter test */
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
        /* chirp mass errors allowed */
        errorParams.dmchirp = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'n':
        /* mass ratio, eta, errors allowed */
        errorParams.deta = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 't':
        /* time coincidence window, argument is in milliseconds */
        errorParams.dt = atof(optarg) * 1000000LL;
        maxTC = errorParams.dt; 
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'q':
        /* start time coincidence window */
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
        /* end time coincidence window */
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
        /* comment */
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

      case 'D':
        /* type of data to analyze */
        if ( ! strcmp( "playground_only", optarg ) )
        {
          dataType = playground_only;
        }
        else if ( ! strcmp( "exclude_play", optarg ) )
        {
          dataType = exclude_play;
        }
        else if ( ! strcmp( "all_data", optarg ) )
        {
          dataType = all_data;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown data type, %s, specified: "
              "(must be playground_only, exclude_play or all_data)\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        haveDataType = 1;
        ADD_PROCESS_PARAM( "string", "%s", optarg );
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

  /* set the gps times startCoinc and endCoinc */
  startCoinc.gpsSeconds = startCoincidence;
  endCoinc.gpsSeconds = endCoincidence;

  if ( ! haveTest )
  {
    fprintf( stderr, "Error: --parameter-test must be specified\n" );
    exit( 1 );
  }


  if ( ! haveDataType )
  {
    fprintf( stderr, "Error: --data-type must be specified\n");
    exit(1);
  }

  if ( ! maxTC )
  {
    fprintf( stderr, "Error: --dt must be specified\n");
    exit(1);
  }

  numIFO = 0;
  if ( g1_trig )
  {
    LALSnprintf( ifoName[numIFO], LIGOMETA_IFO_MAX, "G1" ); 
    numIFO++;
  }
  if ( h1_trig )
  {
    LALSnprintf( ifoName[numIFO], LIGOMETA_IFO_MAX, "H1" ); 
    numIFO++;
  }
  if ( h2_trig )
  {
    LALSnprintf( ifoName[numIFO], LIGOMETA_IFO_MAX, "H2" ); 
    numIFO++;
  }
  if ( l1_trig )
  {
    LALSnprintf( ifoName[numIFO], LIGOMETA_IFO_MAX, "L1" ); 
    numIFO++;
  }
  if ( t1_trig )
  {
    LALSnprintf( ifoName[numIFO], LIGOMETA_IFO_MAX, "T1" ); 
    numIFO++;
  }
  if ( v1_trig )
  {
    LALSnprintf( ifoName[numIFO], LIGOMETA_IFO_MAX, "V1" ); 
    numIFO++;
  }

  /* check that we have at least two IFOs specified, or can't do coincidence */
  if ( numIFO < 2 )
  {
    fprintf( stderr, "Must specify at least two IFOs to do coincidence\n"
        "%d specified\n", numIFO );
    exit ( 1 );
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

  /*
   *
   * read in the input data from the rest of the arguments
   *
   */


  if ( optind < argc )
  {
    for( i = optind; i < argc; ++i )
    {
      struct stat infileStatus;
      INT4 haveSearchSum = 0;
      INT4 numFileTriggers = 0;
      INT4 haveSummValue = 0;
      SummValueTable   *inputSummValue = NULL;
      SnglInspiralTable     *inputData = NULL;
      SearchSummaryTable *inputSummary = NULL;

      /* if the named input file does not exist, exit with an error */
      if ( stat( argv[i], &infileStatus ) == -1 )
      {
        fprintf( stderr, "Error opening input file %s\n", argv[i] );
        perror( "failed to stat() file" );
        exit( 1 );
      }

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


      /* read in the search summary and store */ 
      if ( vrbflg ) fprintf( stdout, 
          "reading search_summary table from file: %s\n", argv[i] );

      haveSearchSum = SearchSummaryTableFromLIGOLw( &inputSummary, argv[i] );

      if ( haveSearchSum < 1 || ! inputSummary )
      {
        if ( vrbflg ) 
          fprintf( stdout, "no valid search_summary table, continuing\n" );
      }
      else
      {
        /* store the search summary table in searchSummList list */
        if ( !searchSummList )
        {
          searchSummList = thisSearchSumm = inputSummary;
        }
        else
        {
          thisSearchSumm = thisSearchSumm->next = inputSummary;
        }
        inputSummary = NULL;
      }


      /* read in the summ_value table and store */
      if ( vrbflg ) fprintf( stdout, 
          "reading summ_value table from file: %s\n", argv[i] );

      haveSummValue = SummValueTableFromLIGOLw( &inputSummValue, argv[i] );

      if ( haveSummValue < 1 || ! inputSummValue )
      {
        if ( vrbflg ) fprintf( stdout, 
            "Unable to read summ_value table from %s\n", argv[i] );
      }
      else
      {
        /* store the summ value table in summValueList list */
        if ( !summValueList )
        {
          summValueList = thisSummValue = inputSummValue;
        }
        else
        {
          thisSummValue = thisSummValue->next = inputSummValue;
        }
        inputSummValue = NULL;

        /* scroll to the end of the linked list of summValues */
        for ( ; thisSummValue->next; thisSummValue = thisSummValue->next );
      }



      /* read in the triggers */
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

        if ( vrbflg ) 
          fprintf( stdout, "got %d sngl_inspiral rows from %s\n", 
              numFileTriggers, argv[i] );

        /* store them */
        if ( ! inspiralEventList )
        {
          /* store the head of the linked list */
          inspiralEventList = thisInspiralTrigger = inputData;
        }
        else
        {
          /* append to the end of the linked list and set current    */
          /* trigger to the first trigger of the list being appended */
          thisInspiralTrigger = thisInspiralTrigger->next = inputData;
        }

        /* scroll to the end of the linked list of triggers */
        for ( ; thisInspiralTrigger->next; thisInspiralTrigger = 
            thisInspiralTrigger->next );

        if ( vrbflg ) fprintf( stdout, "added %d triggers to list\n",
            numFileTriggers );
        numTriggers += numFileTriggers;
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

  if ( vrbflg ) fprintf( stdout, "Read in a total of %d triggers.\n",
      numTriggers );


  /* check that we have read in data for all the requested times
     in all the requested instruments */
  /* for ( j = 0; j < numIFO ; ++j )
     {

     LAL_CALL( LALCheckOutTimeFromSearchSummary ( &status, searchSummList, 
     ifoName[j], &startCoinc, &endCoinc ), &status);
     }
   */
  if ( ! inspiralEventList )
  {
    /* no triggers, so no coincidences can be found */
    if ( vrbflg ) fprintf( stdout,
        "No triggers read in so no coincidences can befound\n" );

    goto cleanexit;
  }

  /* time sort the triggers */
  if ( vrbflg ) fprintf( stdout, "Sorting triggers\n" );
  LAL_CALL( LALSortSnglInspiral( &status, &(inspiralEventList),
        LALCompareSnglInspiralByTime ), &status );

  /* keep only triggers within the requested interval */
  if ( vrbflg ) fprintf( stdout, 
      "Discarding triggers outside requested interval\n" );
  LAL_CALL( LALTimeCutSingleInspiral( &status, &inspiralEventList,
        &startCoinc, &endCoinc), &status );


  /* keep play/non-play/all triggers */
  if ( dataType == playground_only && vrbflg ) fprintf( stdout, 
      "Keeping only playground triggers\n" );
  else if ( dataType == exclude_play && vrbflg ) fprintf( stdout, 
      "Keeping only non-playground triggers\n" );
  else if ( dataType == all_data && vrbflg ) fprintf( stdout, 
      "Keeping all triggers\n" );
  LAL_CALL( LALPlayTestSingleInspiral( &status, &inspiralEventList,
        &dataType ), &status );

  /* scroll to the end of the linked list of triggers, counting triggers */
  thisInspiralTrigger = inspiralEventList;
  for (numTriggers = 0 ; thisInspiralTrigger->next; ++numTriggers,
      thisInspiralTrigger = thisInspiralTrigger->next );
  if ( vrbflg ) fprintf( stdout, 
      "%d remaining triggers after time and data type cut.\n", numTriggers );


  /* 
   *  
   * check for two IFO coincidence
   *
   */

  LAL_CALL( LALCreateTwoIFOCoincList(&status, &coincInspiralList,
        inspiralEventList, &errorParams ), &status); 




  /*
   *
   * write the output xml file
   *
   */


cleanexit:

  /* search summary entries: nevents is from primary ifo */
  if ( inStartTime > 0 && inEndTime > 0 )
  {
    searchsumm.searchSummaryTable->in_start_time.gpsSeconds = inStartTime;
    searchsumm.searchSummaryTable->in_end_time.gpsSeconds = inEndTime;
  }
  searchsumm.searchSummaryTable->out_start_time.gpsSeconds = 
    inStartTime > startCoincidence ? inStartTime : startCoincidence;
  searchsumm.searchSummaryTable->out_end_time.gpsSeconds = 
    inEndTime < endCoincidence ? inEndTime : endCoincidence;
  searchsumm.searchSummaryTable->nnodes = 1;

  if ( vrbflg ) fprintf( stdout, "writing output file... " );

  if ( userTag )
  {
    LALSnprintf( fileName, FILENAME_MAX, "%s-THINCA_%s-%d-%d.xml", 
        ifoName[0], userTag, startCoincidence, 
        endCoincidence - startCoincidence );
  }
  else
  {
    LALSnprintf( fileName, FILENAME_MAX, "%s-THINCA-%d-%d.xml", ifoName[0],
        startCoincidence, endCoincidence - startCoincidence );
  }
  searchsumm.searchSummaryTable->nevents = numTriggers;

  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status , &xmlStream, fileName), 
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

  LAL_CALL( LALBeginLIGOLwXMLTable( &status ,&xmlStream, 
        summ_value_table), &status );
  summValueTable.summValueTable = summValueList;
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, summValueTable,
        summ_value_table), &status );
  LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream), &status );

  /* write the sngl_inspiral table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status ,&xmlStream, 
        sngl_inspiral_table), &status );
  inspiralTable.snglInspiralTable = inspiralEventList;
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, inspiralTable,
        sngl_inspiral_table), &status );
  LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream), &status );

  LAL_CALL( LALCloseLIGOLwXMLFile( &status, &xmlStream), &status );

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

  while ( summValueList )
  {
    thisSummValue = summValueList;
    summValueList = summValueList->next;
    LALFree( thisSummValue );
  }

  while ( inspiralEventList )
  {
    thisInspiralTrigger = inspiralEventList;
    inspiralEventList = inspiralEventList->next;
    LALFree( thisInspiralTrigger );
  }

  while ( coincInspiralList )
  {
    prevCoincInspiral = coincInspiralList;
    coincInspiralList = coincInspiralList->next;
    LALFree( prevCoincInspiral );
  }


  if ( userTag ) free( userTag );

  if ( vrbflg ) fprintf( stdout, "done\n" );

  LALCheckMemoryLeaks();

  exit( 0 );
}
