/*----------------------------------------------------------------------- 
 * 
 * File Name: thinca.c
 *
 * Author: Fairhurst, S.
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
#define PROGRAM_NAME "thinca"

#define INCA_EARG   1
#define INCA_EROW   2
#define INCA_EFILE  3

#define INCA_MSGEARG   "Error parsing arguments"
#define INCA_MSGROW    "Error reading row from XML table"
#define INCA_MSGEFILE  "Could not open file"

#define MAXIFO 4
#define NUMIFO 7

#define KAPPA 1000
#define EPSILON 2

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
  this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
calloc( 1, sizeof(ProcessParamsTable) ); \
LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
    PROGRAM_NAME ); \
LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
    long_options[option_index].name ); \
LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

int haveTrig[NUMIFO];
int checkTimes = 0;


/*
 * 
 * USAGE
 *
 */
static void print_usage(char *program)
{
  fprintf(stderr,
      "Usage:  %s [options] [LIGOLW XML input files]\n" \
      "The following options are recognized.  Options not surrounded in [] are\n" \
      "required.\n" \
      "  [--help]                      display this message\n"\
      "  [--verbose]                   print progress information\n"\
      "  [--version]                   print version information and exit\n"\
      "  [--debug-level]   level       set the LAL debug level to LEVEL\n"\
      "  [--user-tag]      usertag     set the process_params usertag\n"\
      "  [--comment]       string      set the process table comment to STRING\n"\
      "\n"\
      "   --gps-start-time start_time  GPS second of data start time\n"\
      "   --gps-end-time   end_time    GPS second of data end time\n"\
      "  [--check-times]               Check that all times were analyzed\n"\   
      "\n"\
      "  [--g1-triggers]               input triggers from G1\n"\
      "  [--h1-triggers]               input triggers from H1\n"\
      "  [--h2-triggers]               input triggers from H2\n"\
      "  [--l1-triggers]               input triggers from L1\n"\
      "  [--t1-triggers]               input triggers from T1\n"\
      "  [--v1-triggers]               input triggers from V1\n"\
      "\n"\
      "   --parameter-test test        set parameters with which to test coincidence:\n"\
      "                                (m1_and_m2|psi0_and_psi3|mchirp_and_eta)\n"\
      "  [--dm Dm]                     mass coincidence window (default 0)\n"\
      "  [--dpsi0] Dpsi0               psi0 coincidence window\n"\
      "  [--dpsi3] Dpsi3               psi3 coincidence window\n"\
      "  [--dmchirp] Dmchirp           mchirp coincidence window\n"\
      "  [--deta]  Deta                eta coincidence window\n"\
      "   --dt Dt                      time coincidence window (milliseconds)\n"\
      "\n"\
      "   --data-type DATA_TYPE        specify the data type, must be one of\n"\
      "                                (playground_only|exclude_play|all_data)\n"\
      "\n"\
      "[LIGOLW XML input files] list of the input trigger files.\n"\
      "\n", program);
}


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
  CHAR  ifos[LIGOMETA_IFOS_MAX];
  CHAR  comment[LIGOMETA_COMMENT_MAX];
  CHAR *userTag = NULL;

  CHAR  fileName[FILENAME_MAX];

  UINT4  numIFO = 0;
  UINT4  numTrigIFO = 0;
  UINT4  numTriggers = 0;
  UINT4  numCoinc = 0;
  UINT4  numTrigs[NUMIFO];

  SnglInspiralTable    *inspiralEventList = NULL;
  SnglInspiralTable    *thisInspiralTrigger = NULL;
  SnglInspiralTable    *snglOutput;

  EventIDColumn        *eventId;

  CoincInspiralTable   *coincInspiralList = NULL;
  CoincInspiralTable   *thisCoinc = NULL;

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
  MetadataTable         searchSummvarsTable;
  MetadataTable         summValueTable;
  MetadataTable         inspiralTable;
  ProcessParamsTable   *this_proc_param = NULL;
  LIGOLwXMLStream       xmlStream;

  UINT4                  j;
  INT4                   i;

  const CHAR                   ifoList[NUMIFO][LIGOMETA_IFO_MAX] = 
                                   {"??","G1", "H1", "H2", "L1", "T1", "V1"};
  const CHAR                  *ifoArg[NUMIFO] = 
                                   {"??","g1-triggers", "h1-triggers", 
                                         "h2-triggers", "l1-triggers", 
                                         "t1-triggers", "v1-triggers"};


  /* getopt arguments */
  struct option long_options[] =
  {
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"g1-triggers",             no_argument,       &(haveTrig[g1]),   1 },
    {"h1-triggers",             no_argument,       &(haveTrig[h1]),   1 },
    {"h2-triggers",             no_argument,       &(haveTrig[h2]),   1 },
    {"l1-triggers",             no_argument,       &(haveTrig[l1]),   1 },
    {"t1-triggers",             no_argument,       &(haveTrig[t1]),   1 },
    {"v1-triggers",             no_argument,       &(haveTrig[v1]),   1 },
    {"check-times",             no_argument,       &checkTimes,       1 },
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
  /* XXX set kappa and epsilon to default values of 0.5, 10000 to effectively
   * disable the effective distance cut XXX*/
  errorParams.kappa = KAPPA;
  errorParams.epsilon = EPSILON;
  
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
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;


      case 'h':
        /* help message */
        print_usage(argv[0]);
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
        fprintf( stdout, "Inspiral Coincidence\n" 
            "Steve Fairhurst\n"
            "CVS Version: " CVS_ID_STRING "\n"
            "CVS Tag: " CVS_NAME_STRING "\n" );
        exit( 0 );
        break;

      case '?':
        print_usage(argv[0]);
        exit( 1 );
        break;

      default:
        fprintf( stderr, "Error: Unknown error while parsing options\n" );
        print_usage(argv[0]);
        exit( 1 );
    }
  }

  /*
   *
   * check the values of the arguments
   *
   */


  /* Start and End times  */

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


  /* Parameter Test */
  if ( errorParams.test == no_test )
  {
    fprintf( stderr, "Error: --parameter-test must be specified\n" );
    exit( 1 );
  }


  /* Data Type */
  if ( dataType == unspecified_data_type )
  {
    fprintf( stderr, "Error: --data-type must be specified\n");
    exit(1);
  }


  /* Time coincidence window, dt */
  if ( ! errorParams.dt )
  {
    fprintf( stderr, "Error: --dt must be specified\n");
    exit(1);
  }

  /* Store the IFOs we expect triggers from */
  for( j=1; j< NUMIFO; j++)
  {
    if ( haveTrig[j] )
    {
      /* write ifo name in ifoName list */
      LALSnprintf( ifoName[numIFO], LIGOMETA_IFO_MAX, ifoList[j] );
      numIFO++;

      /* store the argument in the process_params table */
      this_proc_param = this_proc_param->next = (ProcessParamsTable *)
        calloc( 1, sizeof(ProcessParamsTable) );
      LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
          "%s", PROGRAM_NAME );
      LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", 
          ifoArg[j]);
      LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
      LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
    }
  }

  
  /* check that we have at least two IFOs specified, or can't do coincidence */
  if ( numIFO < 2 )
  {
    fprintf( stderr, "Must specify at least two IFOs to do coincidence\n"
        "%d specified\n", numIFO );
    exit ( 1 );
  }

  if ( numIFO > 2 )
  {
    fprintf( stdout, 
        "At present, only two IFO coincidence code is available.\n"
        "Will find all double coincidences in %d IFO time.\n",
        numIFO);
  }
  
  /* set ifos to be the alphabetical list of the ifos with triggers */
  if( numIFO == 2 )
  {
    LALSnprintf( ifos, LIGOMETA_IFOS_MAX, "%s%s", ifoName[0], ifoName[1] );
  }
  else if ( numIFO == 3 )
  {
    LALSnprintf( ifos, LIGOMETA_IFOS_MAX, "%s%s%s", ifoName[0], ifoName[1],
        ifoName[2] );
  }
  else if ( numIFO == 4 )
  {
    LALSnprintf( ifos, LIGOMETA_IFOS_MAX, "%s%s%s%s", ifoName[0], ifoName[1],
        ifoName[2], ifoName[3]);
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


  /* store the check-times in the process_params table */
  if ( checkTimes )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--check-times", 
        ifoArg[j]);
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }

  /* delete the first, empty process_params entry */
  this_proc_param = processParamsTable.processParamsTable;
  processParamsTable.processParamsTable = 
    processParamsTable.processParamsTable->next;
  free( this_proc_param );

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
  if ( checkTimes )
  {
    if ( vrbflg ) fprintf( stdout, 
        "Checking that we have data for all times from all IFOs\n");
    for ( j = 0; j < numIFO ; ++j )
    {
      LAL_CALL( LALCheckOutTimeFromSearchSummary ( &status, searchSummList, 
            ifoName[j], &startCoinc, &endCoinc ), &status);
    }
  }

  if ( ! inspiralEventList )
  {
    /* no triggers, so no coincidences can be found */
    if ( vrbflg ) fprintf( stdout,
        "No triggers read in so no coincidences can be found\n" );

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
  for (numTriggers = 0 ; thisInspiralTrigger; ++numTriggers,
      thisInspiralTrigger = thisInspiralTrigger->next );
  if ( vrbflg ) fprintf( stdout, 
      "%d remaining triggers after time and data type cut.\n", numTriggers );


  if ( ! inspiralEventList )
  {
    /* no triggers remaining, so no coincidences can be found */
    if ( vrbflg ) fprintf( stdout,
        "No triggers remain after time/playground cuts\n"
        "No coincidences can be found\n" );

    goto cleanexit;
  }

  /* count the number of triggers for each IFO */
  for( j = 1; j <7; j++)
  {
    LALIfoCountSingleInspiral(&status, &numTrigs[j], inspiralEventList, j);
    if ( vrbflg ) fprintf( stdout, 
        "Have %d triggers from %s.\n", numTrigs[j], ifoList[j] );
    if ( numTrigs[j] && !haveTrig[j] )
    {
      fprintf( stderr, "Read in triggers from %s, none expected.\n",
          ifoList[j]);
    }
    if ( haveTrig[j] && numTrigs[j] )
    {
      ++numTrigIFO;
    }
  }

  if ( !numTrigIFO )
  {
    if ( vrbflg ) fprintf( stdout, "Have no triggers from any IFOs\n"
        "Cannot be coincidences, so exiting without looking.\n");
    goto cleanexit;
  }
  else if ( numTrigIFO==1 )
  {
    if ( vrbflg ) fprintf( stdout, "Have triggers from only one IFO\n"
        "Cannot be coincidences, so exiting without looking.\n");
    goto cleanexit;
  }

  /* 
   *  
   * check for two IFO coincidence
   *
   */

  LAL_CALL( LALCreateTwoIFOCoincList(&status, &coincInspiralList,
        inspiralEventList, &errorParams ), &status); 


  /* count the coincs */
  if( coincInspiralList )
  {  
    for (numCoinc = 1, thisCoinc = coincInspiralList; 
        thisCoinc->next; ++numCoinc, thisCoinc = thisCoinc->next );
  }

  if ( vrbflg ) fprintf( stdout,
      "%d coincident triggers found.\n", numCoinc);


  /*
   *
   * write the output xml file
   *
   */


  /* since we don't yet write coinc inspiral tables, we must make a list of
   * sngl_inspiral tables with the eventId's appropriately poplulated */
  LAL_CALL( LALExtractCoincSngls( &status, &snglOutput, coincInspiralList, 
        &startCoinc), &status );


cleanexit:

  searchsumm.searchSummaryTable->in_start_time = startCoinc;
  searchsumm.searchSummaryTable->in_end_time = endCoinc;
  searchsumm.searchSummaryTable->out_start_time = startCoinc;
  searchsumm.searchSummaryTable->out_end_time = endCoinc;
  searchsumm.searchSummaryTable->nnodes = 1;

  if ( vrbflg ) fprintf( stdout, "writing output file... " );

  if ( userTag )
  {
    LALSnprintf( fileName, FILENAME_MAX, "%s-THINCA_%s-%d-%d.xml", 
        ifos, userTag, startCoincidence, endCoincidence - startCoincidence );
  }
  else
  {
    LALSnprintf( fileName, FILENAME_MAX, "%s-THINCA-%d-%d.xml", ifos,
        startCoincidence, endCoincidence - startCoincidence );
  }
  searchsumm.searchSummaryTable->nevents = numCoinc;

  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status , &xmlStream, fileName), 
      &status );

  /* write process table */

  LALSnprintf( proctable.processTable->ifos, LIGOMETA_IFOS_MAX, ifos );

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
  LALSnprintf( searchsumm.searchSummaryTable->ifos, LIGOMETA_IFOS_MAX, ifos );

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
  if( snglOutput )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status ,&xmlStream, 
          sngl_inspiral_table), &status );
    inspiralTable.snglInspiralTable = snglOutput;
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, inspiralTable,
          sngl_inspiral_table), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream), &status );
  }

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

  while ( searchSummList )
  {
    thisSearchSumm = searchSummList;
    searchSummList = searchSummList->next;
    LALFree( thisSearchSumm );
  }

  while ( summValueList )
  {
    thisSummValue = summValueList;
    summValueList = summValueList->next;
    LALFree( thisSummValue );
  }

  /* free the snglInspirals */
  while ( inspiralEventList )
  {
    thisInspiralTrigger = inspiralEventList;
    inspiralEventList = inspiralEventList->next;
    while ( thisInspiralTrigger->event_id )
    {
      /* free any associated event_id's */
      eventId = thisInspiralTrigger->event_id;
      thisInspiralTrigger->event_id = thisInspiralTrigger->event_id->next;
      LALFree( eventId );
    }
    LALFree( thisInspiralTrigger );
  }

  while ( snglOutput )
  {
    thisInspiralTrigger = snglOutput;
    snglOutput = snglOutput->next;
    while ( thisInspiralTrigger->event_id )
    {
      /* free any associated event_id's */
      eventId = thisInspiralTrigger->event_id;
      thisInspiralTrigger->event_id = thisInspiralTrigger->event_id->next;
      LALFree( eventId );
    }
    LALFree( thisInspiralTrigger );
  }

  while ( coincInspiralList )
  {
    thisCoinc = coincInspiralList;
    coincInspiralList = coincInspiralList->next;
    LALFree( thisCoinc );
  }


  if ( userTag ) free( userTag );

  if ( vrbflg ) fprintf( stdout, "done\n" );

  LALCheckMemoryLeaks();

  exit( 0 );
}
