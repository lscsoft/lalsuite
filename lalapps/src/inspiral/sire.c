/*----------------------------------------------------------------------- 
 * 
 * File Name: sire.c
 *
 * Author: Brady, P. R, Brown, D. A., and Fairhurst, S
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
#include <fcntl.h>
#include <regex.h>
#include <time.h>
#include <glob.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/Date.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXMLRead.h>
#include <lalapps.h>
#include <processtable.h>

RCSID("$Id$");

#define PROGRAM_NAME "sire"
#define CVS_ID_STRING "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"

#define USAGE \
"Usage: lalapps_sire [options]\n"\
"\n"\
"  --help                       display this message\n"\
"  --verbose                    print progress information\n"\
"  --debug-level LEVEL          set the LAL debug level to LEVEL\n"\
"  --user-tag STRING            set the process_params usertag to STRING\n"\
"  --comment STRING             set the process table comment to STRING\n"\
"  --version                    print the CVS version string\n"\
"\n"\
"Input data source:\n"\
"  --glob GLOB                  use pattern GLOB to determine the input files\n"\
"  --input FILE                 read list of input XML files from FILE\n"\
"\n"\
"Output data destination:\n"\
"  --output FILE                write output data to FILE\n"\
"  --tama-output FILE           write out text triggers for tama\n"\
"  --summary-file FILE          write trigger analysis summary to FILE\n"\
"\n"\
"Playground data:\n"\
"  --playground-only            write triggers that are in playground\n"\
"  --exclude-playground         write triggers that are NOT in playground\n"\
"  --all-data                   write triggers for all times read in\n"\
"\n"\
"Clustering and Sorting:\n"\
"  --sort-triggers              time sort the inspiral triggers\n"\
"  --snr-threshold RHO          discard all triggers with snr less than RHO\n"\
"  --cluster-algorithm CHOICE   use trigger clustering algorithm CHOICE\n"\
"                               [ snr_and_chisq | snrsq_over_chisq | snr ]\n"\
"  --cluster-time T             cluster triggers with T ms window\n"\
"\n"\
"Injection analysis:\n"\
"  --injection-file FILE        read injection parameters from FILE\n"\
"  --injection-coincidence T    trigger and injection coincidence window (ms)\n"\
"  --missed-injections FILE     write sim_inspiral for missed injections to FILE\n"\
"  --hardware-injections GPS    assume hardware injections starting at GPS\n"\
"\n"\
"Maintainer flags:\n"\
"  --disable-trig-start-time    don't modify the search summary table\n"

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
  calloc( 1, sizeof(ProcessParamsTable) ); \
  LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
   PROGRAM_NAME ); \
   LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
     long_options[option_index].name ); \
     LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
     LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

#define MAX_PATH 4096

/* function to read the next line of data from the input file list */
char *get_next_line( char *line, size_t size, FILE *fp )
{
  char *s;
  do
    s = fgets( line, size, fp );
  while ( ( line[0] == '#' || line[0] == '%' ) && s );
  return s;
}

int sortTriggers = 0;
LALPlaygroundDataMask dataType;

int main( int argc, char *argv[] )
{
  /* lal initialization variables */
  LALStatus stat = blank_status;
  LALLeapSecAccuracy accuracy = LALLEAPSEC_LOOSE;

  /*  program option variables */
  extern int vrbflg;
  CHAR *userTag = NULL;
  CHAR comment[LIGOMETA_COMMENT_MAX];
  char *inputGlob = NULL;
  char *inputFileName = NULL;
  char *outputFileName = NULL;
  char *tamaFileName = NULL;
  char *summFileName = NULL;
  REAL4 snrStar = -1;
  SnglInspiralClusterChoice clusterchoice = none;
  INT8 cluster_dt = -1;
  char *injectFileName = NULL;
  INT8 inject_dt = -1;
  char *missedFileName = NULL;
  INT4 hardware = 0;
  int  enableTrigStartTime = 1;
  int j;
  FILE *fp = NULL;
  glob_t globbedFiles;
  int numInFiles = 0;
  char **inFileNameList;
  char line[MAX_PATH];

  UINT8 triggerInputTimeNS = 0;

  MetadataTable         proctable;
  MetadataTable         procparams;
  ProcessParamsTable   *this_proc_param;

  int                   numSimEvents = 0;
  int                   numSimInData = 0;
  int                   numSimFound  = 0;
  int                   numSimMissed = 0;
  int                   numSimDiscard = 0;
  int                   numSimProcessed = 0;

  SimInspiralTable     *simEventHead = NULL;
  SimInspiralTable     *thisSimEvent = NULL;
  SimInspiralTable     *missedSimHead = NULL;
  SimInspiralTable     *thisMissedSim = NULL;
  SimInspiralTable     *tmpSimEvent = NULL;
  SimInspiralTable     *prevSimEvent = NULL;

  SearchSummaryTable   *searchSummaryTable = NULL;

  int                   numEvents = 0;
  int                   numEventsInXML = 0;
  int                   numEventsKept = 0;
  int                   numEventsCoinc = 0;
  int                   numEventsDiscard = 0;
  int                   numEventsProcessed = 0;
  int                   numClusteredEvents = 0;

  SnglInspiralTable   **eventHandle = NULL;      
  SnglInspiralTable    *eventHead = NULL;
  SnglInspiralTable    *thisEvent = NULL;
  SnglInspiralTable    *tmpEvent = NULL;
  SnglInspiralTable    *prevEvent = NULL;

  LIGOLwXMLStream       xmlStream;
  MetadataTable         outputTable;


  /*
   *
   * initialization
   *
   */


  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "33" );

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) 
    calloc( 1, sizeof(ProcessTable) );
  LAL_CALL( LALGPSTimeNow ( &stat, &(proctable.processTable->start_time),
        &accuracy ), &stat );
  LAL_CALL( populate_process_table( &stat, proctable.processTable, 
        PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE ), &stat );
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
    calloc( 1, sizeof(ProcessParamsTable) );
  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );


  /*
   *
   * parse command line arguments
   *
   */


  while (1)
  {
    /* getopt arguments */
    static struct option long_options[] = 
    {
      {"verbose",             no_argument,           &vrbflg,              1 },
      {"sort-triggers",       no_argument,     &sortTriggers,              1 },
      {"playground-only",     no_argument, (int *) &dataType, playground_only},
      {"exclude-playground",  no_argument, (int *) &dataType,    exclude_play},
      {"all-data",            no_argument, (int *) &dataType,        all_data},
      {"help",                    no_argument,            0,              'h'},
      {"debug-level",             required_argument,      0,              'z'},
      {"user-tag",                required_argument,      0,              'Z'},
      {"userTag",                 required_argument,      0,              'Z'},
      {"comment",                 required_argument,      0,              'c'},
      {"version",                 no_argument,            0,              'V'},
      {"glob",                    required_argument,      0,              'g'},
      {"input",                   required_argument,      0,              'i'},
      {"output",                  required_argument,      0,              'o'},
      {"tama-output",             required_argument,      0,              'j'},
      {"summary-file",            required_argument,      0,              'S'},
      {"snr-threshold",           required_argument,      0,              's'},
      {"cluster-algorithm",       required_argument,      0,              'C'},
      {"cluster-time",            required_argument,      0,              't'},
      {"injection-file",          required_argument,      0,              'I'},
      {"injection-coincidence",   required_argument,      0,              'T'},
      {"missed-injections",       required_argument,      0,              'm'},
      {"hardware-injections",     required_argument,      0,              'H'},
      {"disable-trig-start-time", no_argument,            0,              'D'},
      {0, 0, 0, 0}
    };
    int c;

    /* getopt_long stores the option index here. */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only ( argc, argv, "hzZ:c:g:i:o:j:S:s:C:Vt:I:T:m:H:D", 
        long_options, &option_index );

    /* detect the end of the options */
    if ( c == - 1 )
      break;

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
          fprintf( stderr, "error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
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

      case 'V':
        fprintf( stdout, "Single Inspiral Reader and Injection Analysis\n"
            "Patrick Brady, Duncan Brown and Steve Fairhurst\n"
            "CVS Version: " CVS_ID_STRING "\n" );
        exit( 0 );
        break;

      case 'g':
        /* create storage for the input file glob */
        optarg_len = strlen( optarg ) + 1;
        inputGlob = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( inputGlob, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "'%s'", optarg );
        break;

      case 'i':
        /* create storage for the input file name */
        optarg_len = strlen( optarg ) + 1;
        inputFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( inputFileName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'o':
        /* create storage for the output file name */
        optarg_len = strlen( optarg ) + 1;
        outputFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( outputFileName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'j':
	/* create storage of the TAMA file name */
	optarg_len = strlen( optarg ) + 1;
        tamaFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( tamaFileName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'S':
        /* create storage for the summ file name */
        optarg_len = strlen( optarg ) + 1;
        summFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( summFileName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 's':
        snrStar = (REAL4) atof( optarg );
        if ( snrStar < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "threshold must be >= 0: "
              "(%f specified)\n",
              long_options[option_index].name, snrStar );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", snrStar );
        break;

      case 'C':
        /* choose the clustering algorithm */
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

      case 't':
        /* cluster time is specified on command line in ms */
        cluster_dt = (INT8) atoi( optarg );
        if ( cluster_dt <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "custer window must be > 0: "
              "(%lld specified)\n",
              long_options[option_index].name, cluster_dt );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%lld", cluster_dt );
        /* convert cluster time from ms to ns */
        cluster_dt *= 1000000LL;
        break;

      case 'I':
        /* create storage for the injection file name */
        optarg_len = strlen( optarg ) + 1;
        injectFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( injectFileName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'T':
        /* injection coincidence time is specified on command line in ms */
        inject_dt = (INT8) atoi( optarg );
        if ( inject_dt < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "injection coincidence window must be >= 0: "
              "(%lld specified)\n",
              long_options[option_index].name, inject_dt );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%lld", inject_dt );
        /* convert inject time from ms to ns */
        inject_dt *= 1000000LL;
        break;

      case 'm':
        /* create storage for the missed injection file name */
        optarg_len = strlen( optarg ) + 1;
        missedFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( missedFileName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'H':
        hardware = (INT4) atoi( optarg );
        if ( hardware <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "GPS start time of hardware injections must be > 0: "
              "(%d specified)\n",
              long_options[option_index].name, hardware );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%ld", hardware );
        break;

      case 'D':
        enableTrigStartTime = 0;
        ADD_PROCESS_PARAM( "string", "%s", " " );
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

  
  /*
   *
   * can use LALCalloc() / LALMalloc() from here
   *
   */


  /* don't buffer stdout if we are in verbose mode */
  if ( vrbflg ) setvbuf( stdout, NULL, _IONBF, 0 );

  /* fill the comment, if a user has specified it, or leave it blank */
  if ( ! *comment )
  {
    LALSnprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
  }
  else
  {
    LALSnprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
  }
  
  /* check that the input and output file names have been specified */
  if ( (! inputGlob && ! inputFileName) || (inputGlob && inputFileName) )
  {
    fprintf( stderr, "exactly one of --glob or --input must be specified\n" );
    exit( 1 );
  }
  if ( ! outputFileName )
  {
    fprintf( stderr, "--output must be specified\n" );
    exit( 1 );
  }

  /* check that if clustering is being done that we have all the options */
  if ( clusterchoice && cluster_dt < 0 )
  {
    fprintf( stderr, "--cluster-time must be specified if --cluster-algorithm "
        "is given\n" );
    exit( 1 );
  }
  else if ( ! clusterchoice && cluster_dt >= 0 )
  {
    fprintf( stderr, "--cluster-algorithm must be specified if --cluster-time "
        "is given\n" );
    exit( 1 );
  }

  /* check that we have all the options to do injections */
  if ( injectFileName && inject_dt < 0 )
  {
    fprintf( stderr, "--injection-coincidence must be specified if "
        "--injection-file is given\n" );
    exit( 1 );
  }
  else if ( ! injectFileName && inject_dt >= 0 )
  {
    fprintf( stderr, "--injection-file must be specified if "
        "--injection-coincidence is given\n" );
    exit( 1 );
  }

  /* save the sort triggers flag */
  if ( sortTriggers )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
      calloc( 1, sizeof(ProcessParamsTable) ); \
      LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
          PROGRAM_NAME ); \
      LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
          "--sort-triggers" );
      LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" ); \
      LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, " " );
  }

  switch ( dataType )
  {
    case playground_only:
      if ( vrbflg )
        fprintf( stdout, "using data from playground times only\n" );
      LALSnprintf( procparams.processParamsTable->program, 
          LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
      LALSnprintf( procparams.processParamsTable->param,
          LIGOMETA_PARAM_MAX, "--playground-only" );
      LALSnprintf( procparams.processParamsTable->type, 
          LIGOMETA_TYPE_MAX, "string" );
      LALSnprintf( procparams.processParamsTable->value, 
          LIGOMETA_TYPE_MAX, " " );
      break;
      
    case exclude_play:
      if ( vrbflg )
        fprintf( stdout, "excluding all triggers in playground times\n" );
      LALSnprintf( procparams.processParamsTable->program, 
          LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
      LALSnprintf( procparams.processParamsTable->param,
          LIGOMETA_PARAM_MAX, "--exclude-play" );
      LALSnprintf( procparams.processParamsTable->type, 
          LIGOMETA_TYPE_MAX, "string" );
      LALSnprintf( procparams.processParamsTable->value, 
          LIGOMETA_TYPE_MAX, " " );
      break;
      
    case all_data:
      if ( vrbflg )
        fprintf( stdout, "using all input data\n" );
      LALSnprintf( procparams.processParamsTable->program, 
          LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
      LALSnprintf( procparams.processParamsTable->param,
          LIGOMETA_PARAM_MAX, "--all-data" );
      LALSnprintf( procparams.processParamsTable->type, 
          LIGOMETA_TYPE_MAX, "string" );
      LALSnprintf( procparams.processParamsTable->value, 
          LIGOMETA_TYPE_MAX, " " );
      break;

    default:
      fprintf( stderr, "data set not defined\n" );
      exit( 1 );
  }


  /*
   *
   * read in the injection XML file, if we are doing an injection analysis
   *
   */


  if ( injectFileName )
  {
    if ( vrbflg ) 
      fprintf( stdout, "reading injections from %s... ", injectFileName );

    numSimEvents = SimInspiralTableFromLIGOLw( &simEventHead, 
        injectFileName, 0, 0 );

    if ( vrbflg ) fprintf( stdout, "got %d injections\n", numSimEvents );

    if ( numSimEvents < 0 )
    {
      fprintf( stderr, "error: unable to read sim_inspiral table from %s\n", 
          injectFileName );
      exit( 1 );
    }

    /* if we are doing hardware injections, increment all the start times */
    if ( hardware )
    {
      if ( vrbflg ) fprintf( stdout, 
          "incrementing GPS times of injections by %d seconds\n", hardware );
      
      for ( thisSimEvent = simEventHead; 
          thisSimEvent; thisSimEvent = thisSimEvent->next )
      {
        thisSimEvent->geocent_end_time.gpsSeconds += hardware;
        thisSimEvent->h_end_time.gpsSeconds       += hardware;
        thisSimEvent->l_end_time.gpsSeconds       += hardware;
      }
    }

    /* discard all injection events that are not in the data we want */
    if ( dataType != all_data )
    {
      numSimDiscard = 0;

      thisSimEvent = simEventHead;
      simEventHead = NULL;
      prevSimEvent = NULL;

      if ( vrbflg ) fprintf( stdout, "discarding injections not in data\n" );

      while ( thisSimEvent )
      {
        INT4 isPlayground;
        LAL_CALL( LALGPSIsPlayground( &stat, &isPlayground, 
              &(thisSimEvent->geocent_end_time)), &stat );

        if ( (dataType == playground_only && isPlayground) || 
            (dataType == exclude_play && ! isPlayground) )
        {
          /* store the head of the linked list */
          if ( ! simEventHead ) simEventHead = thisSimEvent;

          /* keep this event */
          prevSimEvent = thisSimEvent;
          thisSimEvent = thisSimEvent->next;
          ++numSimInData;
          if ( vrbflg ) fprintf( stdout, "+" );
        }
        else
        {
          /* throw this event away */
          tmpSimEvent = thisSimEvent;
          if ( prevSimEvent ) prevSimEvent->next = thisSimEvent->next;
          thisSimEvent = thisSimEvent->next;
          LALFree( tmpSimEvent );
          ++numSimDiscard;
          if ( vrbflg ) fprintf( stdout, "-" );
        }
      }

      if ( vrbflg ) 
        fprintf( stdout, "\nusing %d (discarded %d) of %d injections\n",
            numSimInData, numSimDiscard, numSimEvents );
    }
    else
    {
      if ( vrbflg ) 
        fprintf( stdout, "using all %d injections\n", numSimInData );
      numSimInData = numSimEvents;
    }
  }


  /*
   *
   * read in the input triggers from the xml files
   *
   */


  if ( inputGlob )
  {
    /* use glob() to get a list of the input file names */
    if ( glob( inputGlob, GLOB_ERR, NULL, &globbedFiles ) )
    {
      fprintf( stderr, "error globbing files from %s\n", inputGlob );
      perror( "error:" );
      exit( 1 );
    }

    numInFiles = globbedFiles.gl_pathc;
    inFileNameList = (char **) LALCalloc( numInFiles, sizeof(char *) );

    for ( j = 0; j < numInFiles; ++j )
    {
      inFileNameList[j] = globbedFiles.gl_pathv[j];
    }
  }
  else if ( inputFileName )
  {
    /* read the list of input filenames from a file */
    fp = fopen( inputFileName, "r" );
    if ( ! fp )
    {
      fprintf( stderr, "could not open file containing list of xml files\n" );
      perror( "error:" );
      exit( 1 );
    }

    /* count the number of lines in the file */
    while ( get_next_line( line, sizeof(line), fp ) )
    {
      ++numInFiles;
    }
    rewind( fp );

    /* allocate memory to store the input file names */
    inFileNameList = (char **) LALCalloc( numInFiles, sizeof(char *) );

    /* read in the input file names */
    for ( j = 0; j < numInFiles; ++j )
    {
      inFileNameList[j] = (char *) LALCalloc( MAX_PATH, sizeof(char) );
      get_next_line( line, sizeof(line), fp );
      strncpy( inFileNameList[j], line, strlen(line) - 1);
    }
    
    fclose( fp );
  }
  else
  {
    fprintf( stderr, "no input file mechanism specified\n" );
    exit( 1 );
  }
    
  if ( vrbflg )
  {
    fprintf( stdout, "reading input triggers from:\n" );
    for ( j = 0; j < numInFiles; ++j )
    {
      fprintf( stdout, "%s\n", inFileNameList[j] );
    }
  }


  /*
   *
   * read in the triggers from the input xml files
   *
   */


  if ( injectFileName )
  {
    thisSimEvent = simEventHead;
    simEventHead = NULL;
    prevSimEvent = NULL;
    numSimDiscard = 0;
    numSimInData = 0;
    
    if ( vrbflg ) 
      fprintf( stdout, "discarding injections not in input data\n" );
  }

  for ( j = 0; j < numInFiles; ++j )
  {
    LIGOTimeGPS inPlay, outPlay;
    UINT8 outPlayNS, outStartNS, outEndNS, triggerTimeNS;
    INT4 trigStartTimeArg = 0;

    if ( SearchSummaryTableFromLIGOLw( &searchSummaryTable, 
          inFileNameList[j] ) != 1 )
    {
      fprintf( stderr, 
          "error: zero or multiple search_summary tables in %s\n",
          inFileNameList[j] );
      exit( 1 );
    }

    if ( enableTrigStartTime )
    {
      /* override the value of out_start_time if there is a non-zero */
      /* --trig-start-time option in the process_params table        */
      /* this is necessary to get round a bug in early versions of   */
      /* the inspiral code                                           */

      int mioStatus;
      int pParParam;
      int pParValue;
      struct MetaioParseEnvironment parseEnv;
      const  MetaioParseEnv env = &parseEnv;

      /* open the procress_params table from the input file */
      mioStatus = MetaioOpenTable( env, inFileNameList[j], "process_params" );
      if ( mioStatus )
      {
        fprintf( stderr, "error opening process_params table from file %s\n", 
            inFileNameList[j] );
        exit( 1 );
      }

      /* figure out where the param and value columns are */
      if ( (pParParam = MetaioFindColumn( env, "param" )) < 0 )
      {
        fprintf( stderr, "unable to find column param in process_params\n" );
        MetaioClose(env);
        exit( 1 );
      }
      if ( (pParValue = MetaioFindColumn( env, "value" )) < 0 )
      {
        fprintf( stderr, "unable to find column value in process_params\n" );
        MetaioClose(env);
        exit( 1 );
      }

      /* get the trigger start time from the process params */
      while ( (mioStatus = MetaioGetRow(env)) == 1 )
      {
        if ( ! strcmp( env->ligo_lw.table.elt[pParParam].data.lstring.data, 
              "--trig-start-time" ) )
        {
          trigStartTimeArg = (INT4) 
            atoi( env->ligo_lw.table.elt[pParValue].data.lstring.data );
        }
      }

      MetaioClose( env );

      if ( trigStartTimeArg )
      {
        searchSummaryTable->out_start_time.gpsSeconds = trigStartTimeArg;
        searchSummaryTable->out_start_time.gpsNanoSeconds = 0;
        if ( vrbflg ) fprintf( stdout, "file %s has --trig-start-time %d\n",
            inFileNameList[j], trigStartTimeArg );
      }
    }

    /* compute the out time from the search summary table */
    LAL_CALL( LALGPStoINT8( &stat, 
          &outStartNS, &(searchSummaryTable->out_start_time) ), &stat );
    LAL_CALL( LALGPStoINT8( &stat, 
          &outEndNS, &(searchSummaryTable->out_end_time) ), &stat );
    triggerTimeNS = outEndNS - outStartNS;

    /* check for events and playground */
    if ( dataType != all_data )
    {
      LAL_CALL( LALPlaygroundInSearchSummary( &stat, searchSummaryTable,
            &inPlay, &outPlay ), &stat );
      LAL_CALL( LALGPStoINT8( &stat, &outPlayNS, &outPlay ), &stat );

      if ( dataType == playground_only )
      {
        if ( outPlayNS )
        {
          /* increment the total trigger time by the amount of playground */
          triggerInputTimeNS += outPlayNS;
        }
        else
        {
          /* skip this file as it does not contain any playground data */
          if ( vrbflg )
          {
            fprintf( stdout, "file %s not in playground, continuing\n", 
                inFileNameList[j] );
          }
          LALFree( searchSummaryTable );
          searchSummaryTable = NULL;
          continue;
        }
      }
      else if ( dataType == exclude_play )
      {
        /* increment the total trigger time by the out time minus */
        /* the time that is in the playground                     */
        triggerInputTimeNS += triggerTimeNS - outPlayNS;
      }
    }
    else
    {
      /* increment the total trigger time by the out time minus */
      triggerInputTimeNS += triggerTimeNS;
    }

    if ( injectFileName )
    {
      if ( vrbflg ) fprintf( stdout, "discarding injections not in file: " );

      /* throw away injections that are outside analyzed times */
      while ( thisSimEvent && thisSimEvent->geocent_end_time.gpsSeconds < 
          searchSummaryTable->out_end_time.gpsSeconds )
      {
        /* check if injection is before file start time */
        if ( thisSimEvent->geocent_end_time.gpsSeconds < 
            searchSummaryTable->out_start_time.gpsSeconds )
        {
          /* discard the current injection */
          if ( prevSimEvent ) prevSimEvent->next = thisSimEvent->next;
          tmpSimEvent = thisSimEvent;
          thisSimEvent = thisSimEvent->next;
          LALFree( tmpSimEvent );
          ++numSimDiscard;
          if ( vrbflg ) fprintf( stdout, "-" );
        }
        else
        {
          /* store the head of the linked list */
          if ( ! simEventHead ) simEventHead = thisSimEvent;

          /* keep this injection */
          prevSimEvent = thisSimEvent;
          thisSimEvent = thisSimEvent->next;
          ++numSimInData;
          if ( vrbflg ) fprintf( stdout, "+" );
        }
      }
      if ( vrbflg ) fprintf( stdout, "\n" );
    }

    
    /*
     *
     * if there are any events in the file, read them in
     *
     */
    

    if ( searchSummaryTable->nevents )
    {
      INT4 isPlay;

      if ( vrbflg ) fprintf( stdout, "file %s contains %d events, processing\n",
          inFileNameList[j], searchSummaryTable->nevents );
      
      if ( ! prevEvent )
      {
        eventHandle = &thisEvent;
      }
      else
      {
        eventHandle = &(prevEvent->next);
      }

      /* read the events from the file into a temporary list */
      if ( (numEventsInXML = LALSnglInspiralTableFromLIGOLw( eventHandle, 
              inFileNameList[j], 0, -1 )) < 0 )
      {
        fprintf( stderr, "error: unable to read sngl_inspiral table from %s\n", 
           inFileNameList[j] );
        exit( 1 );
      }
      numEvents += numEventsInXML;

      /* only keep triggers from the data that we want to analyze */
      thisEvent = *eventHandle;
      while ( thisEvent )
      {
        LAL_CALL( LALGPSIsPlayground( &stat, &isPlay, &(thisEvent->end_time) ),
            &stat );

        if ( (dataType == all_data || 
            (dataType == playground_only && isPlay) ||
            (dataType == exclude_play && ! isPlay))
            && ( snrStar < 0 || thisEvent->snr > snrStar) )
        {
          /* keep the trigger and increment the count of triggers */
          if ( ! eventHead ) eventHead = thisEvent;
          prevEvent = thisEvent;
          thisEvent = thisEvent->next;
          ++numEventsKept;
        }
        else
        {
          /* discard the trigger and move to the next one */
          if ( prevEvent ) prevEvent->next = thisEvent->next;
          tmpEvent = thisEvent;
          thisEvent = thisEvent->next;
          LALFree( tmpEvent );
        }
      }

      /* make sure that the linked list is properly terminated */
      if ( prevEvent && prevEvent->next ) prevEvent->next->next = NULL;
    }
    else
    {
      if ( vrbflg ) fprintf( stdout, "file %s contains no events, skipping\n",
          inFileNameList[j] );
    }

    LALFree( searchSummaryTable );
    searchSummaryTable = NULL;
  }

  /* discard the remaining injections which occured after the last file */
  if ( injectFileName )
  {
    if ( vrbflg ) fprintf( stdout, "kept %d injections, discarded %d\n",
        numSimInData, numSimDiscard );

    if ( prevSimEvent ) prevSimEvent->next = NULL;

    numSimDiscard = 0;
    while ( thisSimEvent )
    {
      tmpSimEvent = thisSimEvent;
      thisSimEvent = thisSimEvent->next;
      LALFree( tmpSimEvent );
      ++numSimDiscard;
      if ( vrbflg ) fprintf( stdout, "-" );
    }

    if ( vrbflg ) fprintf( stdout, "\ndiscarded %d injections at end of list\n",
        numSimDiscard );
  }


  /*
   *
   * sort the inspiral events by time
   *
   */


  if ( injectFileName || sortTriggers )
  {
    if ( vrbflg ) fprintf( stdout, "sorting inspiral trigger list..." );
    LAL_CALL( LALSortSnglInspiral( &stat, &eventHead, 
          *LALCompareSnglInspiralByTime ), &stat );
    if ( vrbflg ) fprintf( stdout, "done\n" );
  }


  /*
   *
   * check for events that are coincident with injections
   *
   */


  if ( injectFileName )
  {
    int coincidence = 0;
    UINT8 simTime, inspiralTime;

    if ( vrbflg ) fprintf( stdout, 
        "checking for events that are coincident with injections\n" );

    /* Note: we are assuming that both the inspiral and */
    /* injection events are time sorted                 */
    thisSimEvent = simEventHead;
    thisEvent    = eventHead;
    
    simEventHead = NULL;
    eventHead    = NULL;
    prevSimEvent = NULL;
    prevEvent    = NULL;
    
    numSimFound      = 0;
    numSimDiscard    = 0;
    numEventsDiscard = 0;
    numEventsCoinc   = 0;

    if ( ! thisEvent )
    {
      /* no triggers in the input data, so all injections are missed */
      if ( vrbflg ) fprintf( stdout, "no triggers in input data\n" );

      thisMissedSim = missedSimHead = thisSimEvent;

      while ( thisMissedSim )
      {
        /* count the number of injections just stuck in the missed list */
        if ( vrbflg ) fprintf( stdout, "M" );
        ++numSimMissed;
        ++numSimProcessed;
        thisMissedSim = thisMissedSim->next;
      }
    }
    else
    {
      /* begin loop over the sim_inspiral events */
      while ( thisSimEvent )
      {
        /* compute the end time in nanosec for the injection */
        /* at the relevant detector                          */
        if ( ! strcmp( "L1", thisEvent->ifo ) )
        {
          LAL_CALL( LALGPStoINT8( &stat, &simTime, 
                &(thisSimEvent->l_end_time) ), &stat );
        }
        else if ( ! strcmp( "H1", thisEvent->ifo ) || 
            ! strcmp( "H2", thisEvent->ifo ) )
        {
          LAL_CALL( LALGPStoINT8( &stat, &simTime, 
                &(thisSimEvent->h_end_time) ), &stat );
        }
	else if ( ! strcmp( "G1", thisEvent->ifo ) )
        {
          LAL_CALL( LALGPStoINT8( &stat, &simTime, 
                &(thisSimEvent->g_end_time) ), &stat );
        }
	else if ( ! strcmp( "T1", thisEvent->ifo ) )
        {
          LAL_CALL( LALGPStoINT8( &stat, &simTime, 
                &(thisSimEvent->t_end_time) ), &stat );
        }
	else if ( ! strcmp( "V1", thisEvent->ifo ) )
        {
          LAL_CALL( LALGPStoINT8( &stat, &simTime, 
                &(thisSimEvent->v_end_time) ), &stat );
        }
        else
        {
          fprintf( stderr, "unknown detector found in event list: %s\n", 
              thisEvent->ifo );
          fprintf( stderr, "Detector must be one of (G1|H1|H2|L1|T1|V1)\n");
          exit( 1 );
        }

        /* find the first inspiral event after the current sim event */
        while ( thisEvent )
        {
          coincidence = 0;

          /* compute the time in nanosec for the inspiral */
          LAL_CALL( LALGPStoINT8( &stat, &inspiralTime, 
                &(thisEvent->end_time) ), &stat );

          if ( inspiralTime < (simTime - inject_dt) )
          {
            /* discard this event and move on to the next one */
            if ( prevEvent ) prevEvent->next = thisEvent->next;
            tmpEvent = thisEvent;
            thisEvent = thisEvent->next;
            LALFree( tmpEvent );
            ++numEventsProcessed;
            ++numEventsDiscard;
            if ( vrbflg ) fprintf( stdout, "-" );
          }
          else
          {
            /* we have reached the negative coincincidence window */
            break;
          }
        }

        while ( thisEvent )
        {
          /* compute the time in nanosec for the inspiral */
          LAL_CALL( LALGPStoINT8( &stat, &inspiralTime, 
                &(thisEvent->end_time) ), &stat );

          if ( inspiralTime < (simTime + inject_dt) )
          {
            /* this event is within the coincidence window  */
            /* store this event and move on to the next one */
            if ( ! eventHead ) eventHead = thisEvent;
            prevEvent = thisEvent;
            thisEvent = thisEvent->next;
            coincidence = 1;
            ++numEventsProcessed;
            ++numEventsCoinc;
            if ( vrbflg ) fprintf( stdout, "+" );
          }
          else
          {
            /* we have reached the end of the positive coincincidence window */
            break;
          }
        }

        if ( coincidence )
        {
          /* keep this event in the list and move to the next sim event */
          if ( ! simEventHead ) simEventHead = thisSimEvent;
          prevSimEvent = thisSimEvent;
          ++numSimFound;
          ++numSimProcessed;
          thisSimEvent = thisSimEvent->next;
          if ( vrbflg ) fprintf( stdout, "F" );
        }
        else
        {
          /* save this sim event in the list of missed events... */
          if ( ! missedSimHead )
          {
            missedSimHead = thisMissedSim = thisSimEvent;
          }
          else
          {
            thisMissedSim = thisMissedSim->next = thisSimEvent;
          }

          /* ...and remove it from the list of found events */
          if ( prevSimEvent ) prevSimEvent->next = thisSimEvent->next;
          ++numSimMissed;
          if ( vrbflg ) fprintf( stdout, "M" );

          /* move to the next sim in the list */
          ++numSimProcessed;
          thisSimEvent = thisSimEvent->next;

          /* make sure the missed sim list is terminated */
          thisMissedSim->next = NULL;
        }

        if ( ! thisEvent )
        {
          /* these are no more events to process so all the rest of the */
          /* injections must be put in the missed injections list       */
          if ( ! missedSimHead )
          {
            /* this and any subsequent events are in the missed sim list */
            if ( thisSimEvent ) thisMissedSim = missedSimHead = thisSimEvent;
          }
          else
          {
            if ( thisSimEvent )
            {
              /* append the rest of the list to the list of missed injections */
              thisMissedSim = thisMissedSim->next = thisSimEvent;
            }
            else
            {
              /* there are no injections after this one */
              thisMissedSim = thisMissedSim->next = NULL;
            }
          }

          /* terminate the list of found injections correctly */
          if ( prevSimEvent ) prevSimEvent->next = NULL;

          while ( thisMissedSim )
          {
            /* count the number of injections just stuck in the missed list */
            if ( vrbflg ) fprintf( stdout, "M" );
            ++numSimMissed;
            ++numSimProcessed;
            thisMissedSim = thisMissedSim->next;
          }
          thisSimEvent = NULL;
          break;
        }
      }

      if ( thisEvent )
      {
        /* discard any remaining inspiral triggers -- including thisEvent */
        /* as we have run out of injections */
        tmpEvent = thisEvent;
	if ( prevEvent ) prevEvent->next = NULL;
        while ( tmpEvent )
        {
          thisEvent = tmpEvent;
          tmpEvent = tmpEvent->next;
          LALFree( thisEvent );
          ++numEventsDiscard;
          ++numEventsProcessed;
          if ( vrbflg ) fprintf( stdout, "-" );
        }
      }
    }

    if ( vrbflg )
    {
      fprintf( stdout, "\nfound %d injections, missed %d injections "
          "(%d injections processed)\n",
          numSimFound, numSimMissed, numSimProcessed );

      fprintf( stdout, "found %d coincident events, %d events discarded "
          "(%d events processed)\n",
          numEventsCoinc, numEventsDiscard, numEventsProcessed );
    }

  } /* end if ( injectFileName ) */


  /*
   *
   * cluster the remaining events
   *
   */


  if ( eventHead && clusterchoice )
  {
    if ( vrbflg ) fprintf( stdout, "clustering remaining triggers... " );
    LAL_CALL( LALClusterSnglInspiralTable( &stat, eventHead,
          cluster_dt, clusterchoice ), &stat );
    if ( vrbflg ) fprintf( stdout, "done\n" );
    
    /* count the number of triggers surviving the clustering */
    thisEvent = eventHead;
    numClusteredEvents = 0;
    while ( thisEvent )
    {
      ++numClusteredEvents;
      thisEvent = thisEvent->next;
    }
  }


  /*
   *
   * write output data
   *
   */


  /* write the main output file containing found injections */
  if ( vrbflg ) fprintf( stdout, "writing output xml files... " );
  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
  LAL_CALL( LALOpenLIGOLwXMLFile( &stat, &xmlStream, outputFileName ), &stat );

  /* write out the process and process params tables */
  if ( vrbflg ) fprintf( stdout, "process... " );
  LAL_CALL( LALGPSTimeNow ( &stat, &(proctable.processTable->end_time),
        &accuracy ), &stat );
  LAL_CALL( LALBeginLIGOLwXMLTable( &stat, &xmlStream, process_table ), 
      &stat );
  LAL_CALL( LALWriteLIGOLwXMLTable( &stat, &xmlStream, proctable, 
        process_table ), &stat );
  LAL_CALL( LALEndLIGOLwXMLTable ( &stat, &xmlStream ), &stat );
  free( proctable.processTable );

  /* write the process params table */
  if ( vrbflg ) fprintf( stdout, "process_params... " );
  LAL_CALL( LALBeginLIGOLwXMLTable( &stat, &xmlStream, 
        process_params_table ),	&stat );
  LAL_CALL( LALWriteLIGOLwXMLTable( &stat, &xmlStream, procparams, 
        process_params_table ), &stat );
  LAL_CALL( LALEndLIGOLwXMLTable ( &stat, &xmlStream ), &stat );

  /* Write the found injections to the sim table */
  if ( simEventHead )
  {
    if ( vrbflg ) fprintf( stdout, "sim_inspiral... " );
    outputTable.simInspiralTable = simEventHead;
    LAL_CALL( LALBeginLIGOLwXMLTable( &stat, &xmlStream, 
          sim_inspiral_table ), &stat );
    LAL_CALL( LALWriteLIGOLwXMLTable( &stat, &xmlStream, outputTable, 
          sim_inspiral_table ), &stat );
    LAL_CALL( LALEndLIGOLwXMLTable( &stat, &xmlStream ), &stat );
  }
  
  /* Write the results to the inspiral table */
  if ( eventHead )
  {
    if ( vrbflg ) fprintf( stdout, "sngl_inspiral... " );
    outputTable.snglInspiralTable = eventHead;
    LAL_CALL( LALBeginLIGOLwXMLTable( &stat, &xmlStream, 
          sngl_inspiral_table ), &stat );
    LAL_CALL( LALWriteLIGOLwXMLTable( &stat, &xmlStream, outputTable, 
          sngl_inspiral_table ), &stat );
    LAL_CALL( LALEndLIGOLwXMLTable( &stat, &xmlStream ), &stat);
  }

  /* close the output file */
  LAL_CALL( LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);
  if ( vrbflg ) fprintf( stdout, "done\n" );

  /* write out the TAMA file if it is requested */
  if ( tamaFileName )
  {
    REAL8 trigtime;
    REAL4 mtotal;

    fp = fopen( tamaFileName, "w" );
    if ( ! fp )
    {
      perror( "TAMA file" );
      exit( 1 );
    }

    fprintf( fp, "IFO   trigger time       snr         chisq       "
        " total mass     eta       eff dist (kpc)\n" );

    for ( thisEvent = eventHead; thisEvent; thisEvent = thisEvent->next )
    {
      LAL_CALL( LALGPStoFloat( &stat, &trigtime, &(thisEvent->end_time) ),
          &stat );

      if ( thisEvent->eta <= 0.25 )
      {
        mtotal = thisEvent->mass1 + thisEvent->mass2;
      }
      else
      {
        mtotal = thisEvent->mchirp / pow( thisEvent->eta, 0.6 );
      }

      fprintf( fp, "%s %20.9f %12.6e %12.6e %12.6e %12.6e %12.6e\n", 
          thisEvent->ifo, trigtime, thisEvent->snr, thisEvent->chisq, 
          mtotal, thisEvent->eta, 1.0e+03 * thisEvent->eff_distance );
    }

    fclose( fp );
  }

  if ( missedFileName )
  {
    /* open the missed injections file and write the missed injections to it */
    if ( vrbflg ) fprintf( stdout, "writing missed injections... " );
    memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
    LAL_CALL( LALOpenLIGOLwXMLFile( &stat, &xmlStream, missedFileName ), 
        &stat );

    if ( missedSimHead )
    {
      outputTable.simInspiralTable = missedSimHead;
      LAL_CALL( LALBeginLIGOLwXMLTable( &stat, &xmlStream, sim_inspiral_table ),
          &stat );
      LAL_CALL( LALWriteLIGOLwXMLTable( &stat, &xmlStream, outputTable, 
            sim_inspiral_table ), &stat );
      LAL_CALL( LALEndLIGOLwXMLTable( &stat, &xmlStream ), &stat );
    }

    LAL_CALL( LALCloseLIGOLwXMLFile( &stat, &xmlStream ), &stat );
    if ( vrbflg ) fprintf( stdout, "done\n" );
  }

  if ( summFileName )
  {
    LIGOTimeGPS triggerTime;

    /* write out a summary file */
    fp = fopen( summFileName, "w" );

    switch ( dataType )
    {
      case playground_only:
        fprintf( fp, "using data from playground times only\n" );
        break;
      case exclude_play:
        fprintf( fp, "excluding all triggers in playground times\n" );
        break;
      case all_data:
        fprintf( fp, "using all input data\n" );
        break;
      default:
        fprintf( stderr, "data set not defined\n" );
        exit( 1 );
    }

    fprintf( fp, "read triggers from %d files\n", numInFiles );
    fprintf( fp, "number of triggers in input files: %d \n", numEvents );
    if ( snrStar >= 0 )
    {
      fprintf( fp, "number of triggers in input data with snr above %f: %d \n",
          snrStar, numEventsKept );
    }
    else
    {
      fprintf( fp, "number of triggers in input data %d \n", numEventsKept );
    }

    LAL_CALL( LALINT8toGPS( &stat, &triggerTime, &triggerInputTimeNS ), &stat );
    fprintf( fp, "amount of time analysed for triggers %d sec %d ns\n", 
        triggerTime.gpsSeconds, triggerTime.gpsNanoSeconds );
    
    if ( injectFileName )
    {
      fprintf( fp, "read %d injections from file %s\n", 
          numSimEvents, injectFileName );
      
      fprintf( fp, "number of injections in input data: %d\n", numSimInData );
      fprintf( fp, "number of injections found in input data: %d\n", 
          numSimFound );
      fprintf( fp, 
          "number of triggers found within %lld msec of injection: %d\n",
          (inject_dt / 1000000LL), numEventsCoinc );
      
      fprintf( fp, "efficiency: %f \n", 
          (REAL4) numSimFound / (REAL4) numSimInData );
    }

    if ( clusterchoice )
    {
      fprintf( fp, "number of event clusters with %lld msec window: %d\n",
	  cluster_dt/ 1000000LL, numClusteredEvents ); 
    }
    
    fclose( fp ); 
  }

   
  /*
   *
   * free memory and exit
   *
   */


  /* free the inspiral events we saved */
  while ( eventHead )
  {
    thisEvent = eventHead;
    eventHead = eventHead->next;
    LALFree( thisEvent );
  }

  /* free the process params */
  while( procparams.processParamsTable )
  {
    this_proc_param = procparams.processParamsTable;
    procparams.processParamsTable = this_proc_param->next;
    free( this_proc_param );
  }

  /* free the found injections */
  while ( simEventHead )
  {
    thisSimEvent = simEventHead;
    simEventHead = simEventHead->next;
    LALFree( thisSimEvent );
  }

  /* free the temporary memory containing the missed injections */
  while ( missedSimHead )
  {
    tmpSimEvent = missedSimHead;
    missedSimHead = missedSimHead->next;
    LALFree( tmpSimEvent );
  }
  
  /* free the input file name data */
  if ( inputGlob )
  {
    LALFree( inFileNameList ); 
    globfree( &globbedFiles );
  }
  else
  {
    for ( j = 0; j < numInFiles; ++j )
    {
      LALFree( inFileNameList[j] );
    }
    LALFree( inFileNameList );
  }

  if ( vrbflg ) fprintf( stdout, "checking memory leaks and exiting\n" );
  LALCheckMemoryLeaks();
  exit( 0 );
}
