/*----------------------------------------------------------------------- 
 * 
 * File Name: inca.c
 *
 * Author: Brady, P. R. and Brown, D. A.
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

#define CVS_ID_STRING "$Id$"
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
#define IFOB_SNRTHRESH 6
#define KAPPA 0.01
#define EPSILON 2

/* Usage format string. */
#define USAGE \
"Usage: %s [options] [LIGOLW XML input files]\n\n"\
"  --help                    display this message\n"\
"  --verbose                 print progress information\n"\
"  --version                    print version information and exit\n"\
"  --debug-level LEVEL       set the LAL debug level to LEVEL\n"\
"  --user-tag STRING         set the process_params usertag to STRING\n"\
"  --comment STRING          set the process table comment to STRING\n"\
"\n"\
"  --gps-start-time SEC      GPS second of data start time\n"\
"  --gps-end-time SEC        GPS second of data end time\n"\
"\n"\
"  --silde-time SEC          slide triggers by SEC when determining playground\n"\
"  --slide-time-ns NS        slide triggers by NS when determining playground\n"\
"\n"\
"  --ifo-a ifo_name          name of first ifo (e.g. L1, H1 or H2)\n"\
"  --ifo-b ifo_name          name of second ifo (e.g. L1, H1 or H2)\n"\
"\n"\
"  --triggered-bank FILE     write a triggered bank insted of doing inca\n"\
"  --minimal-match M         set minimal match of triggered bank to M\n"\
"\n"\
"  --epsilon ERROR           set effective distance test epsilon (default 2)\n"\
"  --kappa ERROR             set effective distance test kappa (default 0.01)\n"\
"  --ifo-b-snr-threshold SNR set minimum snr in IFO B (default 6)\n"\
"  --ifo-b-range-cut         test range of IFO B to see if sensitive to trigger\n"\
"  --dm mass                 mass coincidence window (default 0)\n"\
"  --dt time                 time coincidence window (milliseconds)\n"\
"\n"\
"  --no-playground           do not select triggers from playground\n"\
"  --playground-only         only use triggers that are in playground\n"\
"  --write-uniq-triggers     make sure triggers from IFO A are unique\n" \
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


int compareTmpltsByPsi ( const void *a, const void *b )
{
  SnglInspiralTable *aPtr = *((SnglInspiralTable **)a);
  SnglInspiralTable *bPtr = *((SnglInspiralTable **)b);

  if ( aPtr->psi0 > bPtr->psi0 )
  {
    return 1;
  }
  else if ( aPtr->psi0 < bPtr->psi0 )
  {
    return -1;
  }
  else if ( aPtr->psi3 > bPtr->psi3 )
  {
    return 1;
  }
  else if ( aPtr->psi3 < bPtr->psi3 )
  {
    return -1;
  }
  else
  {
    return 0;
  }
}


int main( int argc, char *argv[] )
{
  static LALStatus      status;
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;

  extern int vrbflg;
  static INT4  writeUniqTrigs = 0;
  static INT4  usePlayground = 1;
  INT4  startCoincidence = -1;
  INT4  endCoincidence = -1;
  CHAR  ifoName[MAXIFO][LIGOMETA_IFO_MAX];
  CHAR  comment[LIGOMETA_COMMENT_MAX];
  CHAR *userTag = NULL;

  CHAR  fileName[FILENAME_MAX];
  CHAR *trigBankFile = NULL;
  CHAR *xmlFileName;

  LIGOTimeGPS slideData = {0,0};
  INT8  slideDataNS;
  INT4  numEvents = 0;
  INT4  numIFO;
  INT4  have_ifo_a_trigger = 0;
  INT4	keep_a_trig = 0;
  INT4	dont_search_b = 0;
  INT4  isPlay = 0;
  INT8  ta, tb, tp;
  INT4  numTriggers[MAXIFO];
  INT4  inStartTime = -1;
  INT4  inEndTime = -1;
  REAL4 minMatch = -1;
  INT4  useRangeCut = 0;
  REAL4 ifob_snrthresh = IFOB_SNRTHRESH;
  REAL4	d_range[MAXIFO];

  SnglInspiralTable    *inspiralEventList[MAXIFO];
  SnglInspiralTable    *currentTrigger[MAXIFO];

  SnglInspiralTable    *currentEvent = NULL;
  SnglInspiralTable    *outEvent[MAXIFO];
  SnglInspiralTable    *coincidentEvents[MAXIFO];
  SnglInspiralAccuracy  errorParams;

  SummValueTable       *inspEffRange[MAXIFO];
  SummValueTable       *currentEffRange[MAXIFO];
  
  MetadataTable         proctable;
  MetadataTable         processParamsTable;
  MetadataTable         searchsumm;
  MetadataTable		summValueTable;
  MetadataTable         inspiralTable;
  ProcessParamsTable   *this_proc_param = NULL;
  LIGOLwXMLStream       xmlStream;
  
  INT4                  i, j;

  /* getopt arguments */
  struct option long_options[] =
  {
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"no-playground",           no_argument,       &usePlayground,    0 },
    {"playground-only",         no_argument,       &usePlayground,    1 },
    {"write-uniq-triggers",     no_argument,       &writeUniqTrigs,   1 },
    {"ifo-b-range-cut",		no_argument,       &useRangeCut,      1 },
    {"ifo-a",                   required_argument, 0,                'a'},
    {"ifo-b",                   required_argument, 0,                'b'},
    {"epsilon",                 required_argument, 0,                'e'},
    {"triggered-bank",          required_argument, 0,                'T'},
    {"minimal-match",           required_argument, 0,                'M'},
    {"kappa",                   required_argument, 0,                'k'},
    {"ifo-b-snr-threshold",     required_argument, 0,                'S'},
    {"dm",                      required_argument, 0,                'm'},
    {"dt",                      required_argument, 0,                't'},
    {"gps-start-time",          required_argument, 0,                'q'},
    {"gps-end-time",            required_argument, 0,                'r'},
    {"comment",                 required_argument, 0,                's'},
    {"slide-time",              required_argument, 0,                'X'},
    {"slide-time-ns",           required_argument, 0,                'Y'},
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
  memset( inspiralEventList, 0, MAXIFO * sizeof(SnglInspiralTable *) );
  memset( currentTrigger, 0, MAXIFO * sizeof(SnglInspiralTable *) );
  memset( coincidentEvents, 0, MAXIFO * sizeof(SnglInspiralTable *) );
  memset( outEvent, 0, MAXIFO * sizeof(SnglInspiralTable *) );
  memset( numTriggers, 0, MAXIFO * sizeof(INT4) );
  memset( inspEffRange, 0 , MAXIFO * sizeof(SummValueTable *) ); 

  /* default values */
  errorParams.epsilon = EPSILON;
  errorParams.kappa = KAPPA;

  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    long int gpstime;
    size_t optarg_len;

    c = getopt_long_only( argc, argv, 
        "a:b:e:k:m:t:q:r:s:hz:Z:M:T:S:", long_options, &option_index );

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

      case 'e':
        /* epsilon */
        errorParams.epsilon = atof(optarg);
        if ( errorParams.epsilon < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "epsilon must be non-negative: "
              "(%s given)\n", 
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'k':
        /* kappa */
        errorParams.kappa = atof(optarg);
        if ( errorParams.kappa < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "epsilon must be non-negative: "
              "(%s given)\n", 
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'S':
        /* set the snr threshold in ifo b.  Used when deciding if ifo b
	 * could have seen the triggers. */
        ifob_snrthresh = atof(optarg);
        if ( ifob_snrthresh < 0.0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "IFO B snr threshold must be positive"
              "(%s given)\n", 
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'm':
        /* mass errors allowed */
        errorParams.dm = atof(optarg);
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

      case 'T':
        optarg_len = strlen( optarg ) + 1;
        trigBankFile = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( trigBankFile, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'M':
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
            "Duncan Brown <duncan@gravity.phys.uwm.edu>\n"
            "CVS Version: " CVS_ID_STRING "\n" );
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

  /* check for minimal match when doing a triggered bank */
  if ( trigBankFile && minMatch < 0 )
  {
    fprintf( stderr, "--minimal-match must be specified\n" );
    exit( 1 );
  }
  
  /* fill the comment, if a user has specified on, or leave it blank */
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

  /* store the ifo b range cut option */
  if ( useRangeCut )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--ifo-b-range-cut" );
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


  /* decide how many ifos we have based on what we are doing */
  if ( trigBankFile )
  {
    numIFO = 1;
  }
  else
  {
    numIFO = 2;
  }
    
  /* calculate the slide time in nanoseconds */
  LAL_CALL( LALGPStoINT8( &status, &slideDataNS, &slideData ), &status );


  /*
   *
   * read in the input data from the rest of the arguments
   *
   */


  if ( optind < argc )
  {
    for( i = optind; i < argc; ++i )
    {
      INT4 haveSearchSum = 0;
      INT4 numFileTriggers = 0;
      SnglInspiralTable *inputData = NULL;
      SearchSummaryTable *inputSummary = NULL;

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
        if ( inStartTime < 0 || 
            inputSummary->out_start_time.gpsSeconds < inStartTime )
        {
          inStartTime = inputSummary->out_start_time.gpsSeconds;
        }

        if ( inEndTime < 0 ||
            inputSummary->out_end_time.gpsSeconds > inEndTime )
        {
          inEndTime = inputSummary->out_end_time.gpsSeconds;
        }

        LALFree( inputSummary );
        inputSummary = NULL;
      }
      
      if ( numIFO == 2 )
      {
	INT4 haveSummValue = 0;
	SummValueTable *thisSummValue = NULL;

	if ( vrbflg ) fprintf( stdout, 
	      "reading summ_value table from file: %s\n", argv[i] );
	
	haveSummValue = SummValueTableFromLIGOLw( &thisSummValue, argv[i] );
	
	if ( haveSummValue < 1 || ! thisSummValue )
	{
	  if ( vrbflg ) fprintf( stdout, 
	      "Unable to read summ_value table from %s\n", argv[i] );
	}
	else
	{
	  INT4 knownIFO = 0;
	  SummValueTable *tempSummValue = NULL;

	  if ( vrbflg ) fprintf( stdout, 
	     "checking summ_value table for inspiral effective distance\n" );
 
	  while ( thisSummValue )
          {  
	    if ( strncmp( thisSummValue->name, "inspiral_effective_distance",
		  LIGOMETA_SUMMVALUE_NAME_MAX ) )
	    {
	      /* not an effective distance -- discard */
	      tempSummValue = thisSummValue;
	      thisSummValue = thisSummValue->next;
	      LALFree( tempSummValue );
	    }
	    else
	    {
	      /* check that effective distance was calculated using 
	       * 1.4_1.4 solar mass inspiral and snr = 8 */
	      if ( strncmp( thisSummValue->comment, "1.4_1.4_8",
		  LIGOMETA_SUMMVALUE_COMM_MAX ) )
	      {
		fprintf( stdout, "effective distance not calculated\n");
		fprintf( stdout, "using 1.4-1.4 solar mass, snr = 8\n");
		fprintf( stdout, "comment was %s\n", thisSummValue->comment );
		tempSummValue = thisSummValue;
		thisSummValue = thisSummValue->next;
		LALFree( tempSummValue );
	      }
	      else
	      {
		if ( vrbflg )
		{  
		  fprintf( stdout, "got inspiral effective distance of %f ",
		     thisSummValue->value );
		  fprintf( stdout, "between %d and %d GPS secs for ifo %s\n",
		    thisSummValue->start_time.gpsSeconds, 
		    thisSummValue->end_time.gpsSeconds, thisSummValue->ifo );
		}
		/* locate the ifo associated to this summ_value and store it */
		for ( j = 0; j < numIFO ; ++j )
		{
		  if ( ! strncmp( ifoName[j], thisSummValue->ifo,
		      LIGOMETA_IFO_MAX ) )
		  {
		    knownIFO = 1;

		    if ( ! inspEffRange[j] )
		    {
		      /* store the head of the linked list */
		      inspEffRange[j] = currentEffRange[j] = thisSummValue;
		    }
		    else
		    {
		      /* append to the end of the linked list */
		      currentEffRange[j] = currentEffRange[j]->next = 
			  thisSummValue;
		    }
		    thisSummValue = thisSummValue->next;
		    currentEffRange[j]->next = NULL;
		    break;
		  }
		}
		if ( ! knownIFO )
		{
		  /* catch an unknown ifo name among the input files */
		  if ( vrbflg ) fprintf( stdout, 
		      "Unknown interferometer %s\n, discarding", 
		    thisSummValue->ifo );
		  tempSummValue = thisSummValue;
		  thisSummValue = thisSummValue->next;
		  LALFree( tempSummValue );
		}
	      } /* close for ( j = 0; j < numIFO ; ++j ) */
	    }
          } /* close while ( thisSummValue ) */
	}
      } /* close if( useRangeCut ) */

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
              /* append to the end of the linked list */
              currentTrigger[j]->next = inputData;
            }
            while ( currentTrigger[j]->next )
            {
              /* spin on to the end of the linked list */
              currentTrigger[j] = currentTrigger[j]->next;
            }

            /* store number of triggers from ifo a for trigtotmplt algorithm */
            if ( j == 0 ) 
            {
              numEvents += numFileTriggers;
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
      fprintf( stdout, "No triggers read in for interferometer %d\n", j );
      goto cleanexit;
    }
  }


  /*
   *
   * code for generating a triggered bank
   *
   */


  if ( trigBankFile )
  {
    SnglInspiralTable   **eventHandle = NULL;
    SnglInspiralTable    *thisEvent = NULL;
    SnglInspiralTable    *prevEvent = NULL;
    
    eventHandle = (SnglInspiralTable **) 
      LALCalloc( numEvents, sizeof(SnglInspiralTable *) );

    for ( i = 0, thisEvent = inspiralEventList[0]; i < numEvents; 
        ++i, thisEvent = thisEvent->next )
    {
      eventHandle[i] = thisEvent;
    }
    if ( vrbflg ) fprintf( stdout, "sorting events by mass... " );
    qsort( eventHandle, numEvents, sizeof(eventHandle[0]), 
        compareTmpltsByMass );
    if ( vrbflg ) fprintf( stdout, "done\n" );

    /* create a linked list of sorted templates */
    coincidentEvents[0] = prevEvent = eventHandle[0];
    numTriggers[0] = 1;
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
        ++numTriggers[0];
      }
    }
    prevEvent->next = NULL;

    if ( vrbflg ) fprintf( stdout, "found %d sngl_inspiral rows for bank %s\n", 
        numTriggers[0], trigBankFile );

    LALFree( eventHandle );

    /* skip the inca code and write out the bank */
    goto cleanexit;
  }


  /*
   *
   * sort the input data by time
   *
   */


  for ( j = 0; j < numIFO; ++j )
  {
    if ( vrbflg ) fprintf( stdout, "Sorting triggers from ifo %d\n", j );
    LAL_CALL( LALSortSnglInspiral( &status, &(inspiralEventList[j]),
          LALCompareSnglInspiralByTime ), &status );
  }


  /*
   * 
   * find the first trigger after coincidence start time for ifo A
   *
   */


  if ( vrbflg ) fprintf( stdout, "Moving to first trigger in window\n" );

  for ( j = 0; j < numIFO; ++j )
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


  if ( vrbflg ) fprintf( stdout, "start loop over ifo A\n" );

  while ( (currentTrigger[0] ) && 
      (currentTrigger[0]->end_time.gpsSeconds < endCoincidence) )
  {
    if ( vrbflg ) fprintf( stdout, "  using IFO A trigger at %d + %10.10f\n",
        currentTrigger[0]->end_time.gpsSeconds, 
        ((REAL4) currentTrigger[0]->end_time.gpsNanoSeconds * 1e-9) );

    LAL_CALL( LALGPStoINT8( &status, &ta, &(currentTrigger[0]->end_time) ), 
        &status );

    /* use the proper trigger times to determine playground */
    tp = ta - slideDataNS;

    LAL_CALL( LALINT8NanoSecIsPlayground( &status, &isPlay, &tp ), &status );

    if ( vrbflg )
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

      if ( tb > ta - errorParams.dt )
      {
        /* we have reached the time coinicidence window */
        break;
      }

      currentTrigger[1] = currentTrigger[1]->next;
    }

    /* if we are playground only and the trigger is in playground or we are */
    /* not using playground and the trigger is not in the playground...     */
    if ( ( usePlayground && isPlay ) || ( ! usePlayground && ! isPlay)) 
    {
 
      /* determine whether we should expect to see a trigger in ifo b */
      if ( useRangeCut )
      {
	REAL4 lower_limit = 0;
	REAL4 upper_limit = 0;
	/* get the relevant values of inspiral_effective_distance from */
	/* the summ_value table */
	for ( j = 0; j < numIFO; ++j )
	{
	  currentEffRange[j]=inspEffRange[j];
	  d_range[j] = 0;
	  while ( currentEffRange[j] )
	  {
	    INT8 ts, te;
	    LAL_CALL( LALGPStoINT8( &status, &ts, 
		&(currentEffRange[j]->start_time) ), &status );
	    LAL_CALL( LALGPStoINT8( &status, &te, 
		&(currentEffRange[j]->end_time) ), &status );
	  
	    if ( (ts < ta) && (ta < te) )
	    {
	      /* use this value of inspiral_effective_distance */
	      d_range[j] = currentEffRange[j]->value;
	      if( vrbflg ) fprintf( stdout,
		"range for %s is %f Mpc\n",  ifoName[j], d_range[j]);
	      break;
	    }
	    currentEffRange[j] = currentEffRange[j]->next;
	  }
	  if ( d_range[j] <= 0 )
	  {
	    fprintf( stderr, "error: unable to find range for %s\n", 
		ifoName[j]);
	    exit( 1 );
	  }
	}

	/* test whether we expect to be able to see anything in IFO B */
	/* calculate lowest and highest allowed SNRs in IFO B */
	lower_limit = ( ( d_range[1] / d_range[0] ) * currentTrigger[0]->snr 
	    - errorParams.epsilon) / ( 1 + errorParams.kappa);
	if ( errorParams.kappa < 1 )
	{
	  upper_limit = ( ( d_range[1] / d_range[0] ) * currentTrigger[0]->snr 
	  + errorParams.epsilon) / ( 1 - errorParams.kappa);
	}
	else 
	{
	  upper_limit = 0;
	}
      
	if ( vrbflg ) 
	{
	  fprintf( stdout, 
	    "trigger in IFO B expected to have SNR between %f and %f\n",
	    lower_limit, upper_limit );
	  fprintf( stdout, "SNR threshold in IFO B is %f\n", ifob_snrthresh );
	}
	if ( ifob_snrthresh <= lower_limit )
	{
	  if ( vrbflg ) fprintf( stdout, 
	      "looking for a coincident trigger in IFO B\n" );
	  keep_a_trig = 0;
	  dont_search_b = 0;
	}
	else if ( upper_limit  && ( ifob_snrthresh > upper_limit ) )
	{
	  if ( vrbflg ) fprintf( stdout, 
	      "don't expect a trigger in IFO B, keep IFO A trigger\n" );
	  keep_a_trig = 1;
	  dont_search_b = 1;
	}
	else
	{
	  if ( vrbflg ) fprintf( stdout, 
	      "we may see a trigger in IFO B, keep IFO A regardless\n" );
	  keep_a_trig = 1;
	  dont_search_b = 0;
	}
      } /* closes if ( useRangeCut ) */

      if ( ! dont_search_b )
      {
	if ( vrbflg &&  currentTrigger[1] ) fprintf( stdout, 
	    "  start loop over IFO B trigger at %d + %10.10f\n",
            currentTrigger[1]->end_time.gpsSeconds, 
            ((REAL4)currentTrigger[1]->end_time.gpsNanoSeconds * 1e-9) );
	  
	/* look for coincident events in B within the time window */
	currentEvent = currentTrigger[1];

	while ( currentTrigger[1] )
	{
	  LAL_CALL( LALGPStoINT8( &status, &tb, 
		&(currentTrigger[1]->end_time) ), &status );

	  if (tb > ta + errorParams.dt )
	  {
	    /* we are outside the time coincidence so move to the next event */
	    break;
	  }
	  else
	  {
	    /* call the LAL function which compares events parameters */
	    if ( vrbflg ) fprintf( stdout, 
		"    comparing IFO B trigger at %d + %10.10f\n",
                currentTrigger[1]->end_time.gpsSeconds, 
                ((REAL4)currentTrigger[1]->end_time.gpsNanoSeconds * 1e-9) );

	    LAL_CALL( LALCompareSnglInspiral( &status, currentTrigger[0],
                currentTrigger[1], &errorParams ), &status );
	  }

	  if ( errorParams.match )
	  {
	    /* store this event for output */
	    if ( vrbflg )
	      fprintf( stdout, "    >>> found coincidence <<<\n" );

	    for ( j = 0; j < numIFO; ++j )
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

		++numTriggers[j];
		have_ifo_a_trigger = 1;
	      }
	    }
	  }

	  currentTrigger[1] = currentTrigger[1]->next;

	} /* end loop over current events */

	/* go back to saved current IFO B trigger */
	currentTrigger[1] = currentEvent;

      } /* end loop over ! dont_search_b */
        
      
      if ( keep_a_trig && ! have_ifo_a_trigger )
      {
	if ( vrbflg ) fprintf( stdout, 
	  "kept trigger from IFO A although no coincident trigger in IFO B\n");
	if ( ! coincidentEvents[0] )
	{
	  coincidentEvents[0] = outEvent[0] = (SnglInspiralTable *) 
	      LALMalloc( sizeof(SnglInspiralTable) );
	}
	else
	{
	  outEvent[0] = outEvent[0]->next = (SnglInspiralTable *) 
	      LALMalloc( sizeof(SnglInspiralTable) );         
	}

	memcpy( outEvent[0], currentTrigger[0],  sizeof(SnglInspiralTable) );
	outEvent[0]->next = NULL;

	++numTriggers[0];
      }

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

  for ( j = 0; j < numIFO; ++j )
  {

    /* set the file name correctly */
    if ( trigBankFile )
    {
      xmlFileName = trigBankFile;
    }
    else
    {
      if ( userTag )
      {
        LALSnprintf( fileName, FILENAME_MAX, "%s-INCA_%s-%d-%d.xml", ifoName[j],
            userTag, startCoincidence, endCoincidence - startCoincidence );
      }
      else
      {
        LALSnprintf( fileName, FILENAME_MAX, "%s-INCA-%d-%d.xml", ifoName[j],
            startCoincidence, endCoincidence - startCoincidence );
      }

      xmlFileName = fileName;
    }

    searchsumm.searchSummaryTable->nevents = numTriggers[j];

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

    /* write the summ_value table for ifoName[j] */
    if ( inspEffRange[j] )
    {
      LAL_CALL( LALBeginLIGOLwXMLTable( &status ,&xmlStream, 
            summ_value_table), &status );
      summValueTable.summValueTable = inspEffRange[j];
      LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, summValueTable,
            summ_value_table), &status );
      LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream), &status );
    }

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

  while( processParamsTable.processParamsTable )
  {
    this_proc_param = processParamsTable.processParamsTable;
    processParamsTable.processParamsTable = this_proc_param->next;
    free( this_proc_param );
  }

  for( j = 0; j < numIFO; ++j )
  {
    
    while ( inspEffRange[j] )
    {
      currentEffRange[j] = inspEffRange[j];
      inspEffRange[j] = inspEffRange[j]->next;
      LALFree( currentEffRange[j] );
    }
      
    while ( coincidentEvents[j] )
    {
      currentEvent = coincidentEvents[j];
      coincidentEvents[j] = coincidentEvents[j]->next;
      LALFree( currentEvent );
    }

    if ( ! trigBankFile )
    {
      while ( inspiralEventList[j] )
      {
        currentEvent = inspiralEventList[j];
        inspiralEventList[j] = inspiralEventList[j]->next;
        LALFree( currentEvent );
      }
    }
  }

  if ( userTag ) free( userTag );

  if ( vrbflg ) fprintf( stdout, "done\n" );

  LALCheckMemoryLeaks();

  exit( 0 );
}
