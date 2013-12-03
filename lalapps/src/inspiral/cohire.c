/*-----------------------------------------------------------------------
 *
 * File Name: cohire.c
 *
 * Author: Bose, S.
 * 
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
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/Date.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataInspiralUtils.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/Segments.h>
#include <lal/SegmentsIO.h>
#include <LALAppsVCSInfo.h>
#include <lalapps.h>
#include <processtable.h>

#define PROGRAM_NAME "sire"
#define CVS_ID_STRING "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
  this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
calloc( 1, sizeof(ProcessParamsTable) ); \
snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
    PROGRAM_NAME ); \
snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
    long_options[option_index].name ); \
snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

#define MAX_PATH 4096

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/*
 *
 * USAGE
 *
 */


static void print_usage(char *program)
{
  fprintf(stderr,
      "Usage: %s [options] [LIGOLW XML input files]\n"\
      "The following options are recognized.  Options not surrounded in []\n"\
      "are required.\n", program );
  fprintf(stderr,
      " [--help                       display this message\n"\
      " [--verbose                    print progress information\n"\
      " [--user-tag]      usertag     set the process_params usertag\n"\
      " [--comment]       string      set the process table comment to string\n"\
      " [--version]                   print the CVS version string\n"\
      "\n"\
      "Output data destination:\n"\
      "  --output         output_file write output data to output_file\n"\
      " [--summary-file]  summ_file   write trigger analysis summary to summ_file\n"\
      "\n"\
      "Playground data:\n"\
      "  --data-type      datatype    specify the data type, must be one of\n"\
      "                               (playground_only|exclude_play|all_data)\n"\
      "\n"\
      "Cuts and Vetos:\n"\
      " [--ifo-cut]       ifo         only keep triggers from specified ifo\n"\
      " [--coinc-cut]         ifos     only keep triggers from IFOS\n"\
      " [--extract-slide]     slide    only keep triggers from specified slide\n"\
      " [--num-slides]        slides   number of time slides performed \n"\
      " [--snr-threshold] snr_star    discard all triggers with snr less than snr_star\n"\
      " [--rsq-threshold] rsq_thresh  discard all triggers whose rsqveto_duration\n"\
      "                               exceeds rsq_thresh\n"\
      " [--rsq-max-snr]   rsq_max_snr apply rsq on triggers with snr < rsq_max_snr\n"\
      "                               exceeds rsq_thresh\n"\
      " [--rsq-coeff]     rsq_coeff   apply rsq on triggers with snr > rsq_max_snr\n"\
      "                               exceeds rsq_coeff * snr ^ rsq_power\n"\
      " [--rsq-power]     rsq_power   apply rsq on triggers with snr > rsq_max_snr\n"\
      "                               exceeds rsq_coeff * snr ^ rsq_power\n"\
      " [--veto-file]     veto_file   discard all triggers which occur during times\n"\
      "                               contained in the segments of the veto_file\n"\
      "\n"\
      "Sorting and Clustering:\n"\
      " [--sort-triggers]             time sort the inspiral triggers\n"\
      " [--cluster-time]   clust_time cluster triggers with clust_time ms window\n"\
      " [--cluster-algorithm] clust   use trigger clustering algorithm clust\n"\
      "                               [ cohsnr | effCohSnr | nullstat | snrByNullstat\n"\
      "                               | autoCorrCohSqByNullstat | crossCorrCohSqByNullstat\n"\
      "                               | autoCorrNullSqByNullstat| crossCorrNullSqByNullstat ]\n"\
      "\n"\
      "Injection analysis:\n"\
      " [--injection-file]   inj_file read injection parameters from inj_file\n"\
      " [--injection-window] inj_win  trigger and injection coincidence window (ms)\n"\
      " [--missed-injections] missed  write sim_inspiral for missed injections to FILE\n");
  fprintf( stderr, "\n");
  fprintf( stderr, "[LIGOLW XML input files] list of the input trigger files.\n");
}

int sortTriggers = 0;
LALPlaygroundDataMask dataType;
extern int vrbflg;

int main( int argc, char *argv[] )
{
  /* lal initialization variables */
  LALStatus status = blank_status;

  /*  program option variables */
  CHAR *userTag = NULL;
  CHAR comment[LIGOMETA_COMMENT_MAX];
  char *ifos = NULL;
  char *ifoName = NULL;
  char *outputFileName = NULL;
  char *summFileName = NULL;
  char *injectFileName = NULL;
  char *vetoFileName = NULL;
  char *missedFileName = NULL;
  REAL4 snrStar = -1;
  REAL4 rsqVetoThresh = -1;
  REAL4 rsqMaxSnr     = -1;
  REAL4 rsqAboveSnrCoeff = -1;
  REAL4 rsqAboveSnrPow     = -1;
  LALSegList vetoSegs;
  MultiInspiralClusterChoice clusterchoice = no_statistic;
  INT8 cluster_dt = -1;
  INT8 injectWindowNS = -1;
  int j;
  FILE *fp = NULL;
  int numInFiles = 0;

  UINT8 triggerInputTimeNS = 0;

  MetadataTable         proctable;
  MetadataTable         procparams;
  ProcessParamsTable   *this_proc_param;

  SimInspiralTable     *simEventHead = NULL;
  SimInspiralTable     *thisSimEvent = NULL;
  SimInspiralTable     *missedSimHead = NULL;
  SimInspiralTable     *tmpSimEvent = NULL;

  SearchSummvarsTable  *inputFiles = NULL;

  SearchSummaryTable   *searchSummList = NULL;
  SearchSummaryTable   *thisSearchSumm = NULL;
  SummValueTable       *summValueList = NULL;

  int                   extractSlide = 0;
  int                   numSlides = 0;
  int                   numEvents = 0;
  int                   numEventsKept = 0;
  int                   numEventsInIFO = 0;
  int                   numEventsAboveSNRThresh = 0;
  int                   numEventsBelowRsqThresh = 0;
  int                   numEventsSurvivingVeto = 0;
  int                   numClusteredEvents = 0;
  int                   numEventsInIfos = 0;

  int                   numSimEvents = 0;
  int                   numSimInData = 0;
  int                   numSimFound  = 0;
  int                   numMultiFound  = 0;

  MultiInspiralTable   *missedHead = NULL;
  MultiInspiralTable   *thisEvent = NULL;
  MultiInspiralTable   *thisInspiralTrigger = NULL;
  MultiInspiralTable   *inspiralEventList = NULL;
  MultiInspiralTable   *slideEvent = NULL;

  LIGOLwXMLStream       xmlStream;
  MetadataTable         outputTable;
  MetadataTable         UNUSED savedEvents;
  MetadataTable         searchSummvarsTable;

  /*
   *
   * initialization
   *
   */


  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *)
    calloc( 1, sizeof(ProcessTable) );
  XLALGPSTimeNow(&(proctable.processTable->start_time));
  XLALPopulateProcessTable(proctable.processTable, PROGRAM_NAME, LALAPPS_VCS_IDENT_ID,
      LALAPPS_VCS_IDENT_STATUS, LALAPPS_VCS_IDENT_DATE, 0);
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *)
    calloc( 1, sizeof(ProcessParamsTable) );
  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );

  savedEvents.multiInspiralTable = NULL;


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
      {"help",                    no_argument,            0,              'h'},
      {"user-tag",                required_argument,      0,              'Z'},
      {"userTag",                 required_argument,      0,              'Z'},
      {"comment",                 required_argument,      0,              'c'},
      {"version",                 no_argument,            0,              'V'},
      {"data-type",               required_argument,      0,              'k'},
      {"output",                  required_argument,      0,              'o'},
      {"summary-file",            required_argument,      0,              'S'},
      {"extract-slide",           required_argument,      0,              'e'},
      {"num-slides",              required_argument,      0,              'N'},
      {"snr-threshold",           required_argument,      0,              's'},
      {"rsq-threshold",           required_argument,      0,              'r'},
      {"rsq-max-snr",             required_argument,      0,              'R'},
      {"rsq-coeff",               required_argument,      0,              'p'},
      {"rsq-power",               required_argument,      0,              'P'},
      {"cluster-algorithm",       required_argument,      0,              'C'},
      {"cluster-time",            required_argument,      0,              't'},
      {"ifo-cut",                 required_argument,      0,              'd'},
      {"coinc-cut",               required_argument,      0,              'D'},
      {"veto-file",               required_argument,      0,              'v'},
      {"injection-file",          required_argument,      0,              'I'},
      {"injection-window",        required_argument,      0,              'T'},
      {"missed-injections",       required_argument,      0,              'm'},
      {0, 0, 0, 0}
    };
    int c;

    /* getopt_long stores the option index here. */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only ( argc, argv,
        "c:d:D:hj:k:m:o:r:s:t:v:C:DH:I:R:ST:VZ:",
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
        print_usage(argv[0]);
        exit( 0 );
        break;

      case 'Z':
        /* create storage for the usertag */
        optarg_len = strlen( optarg ) + 1;
        userTag = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( userTag, optarg, optarg_len );

        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s",
            PROGRAM_NAME );
        snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "-userTag" );
        snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
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
          snprintf( comment, LIGOMETA_COMMENT_MAX, "%s", optarg);
        }
        break;

      case 'V':
        fprintf( stdout, "Coherent Inspiral Reader and Injection Analysis\n"
            "Sukanta Bose\n");
        XLALOutputVersionString(stderr, 0);
        exit( 0 );
        break;

      case 'o':
        /* create storage for the output file name */
        optarg_len = strlen( optarg ) + 1;
        outputFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( outputFileName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'e':
        /* store the number of slides */
        extractSlide = atoi( optarg );
        if ( extractSlide == 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "extractSlide must be non-zero: "
              "(%d specified)\n",
              long_options[option_index].name, extractSlide );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", extractSlide );
        break;

      case 'N':
        /* store the number of slides */
        numSlides = atoi( optarg );
        if ( numSlides < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "numSlides >= 0: "
              "(%d specified)\n",
              long_options[option_index].name, numSlides );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", numSlides );
        break;

      case 'S':
        /* create storage for the summ file name */
        optarg_len = strlen( optarg ) + 1;
        summFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( summFileName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'k':
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

      case 'r':
        rsqVetoThresh = (REAL4) atof( optarg );
        if ( rsqVetoThresh < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "threshold must be >= 0: "
              "(%f specified)\n",
              long_options[option_index].name, rsqVetoThresh );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", rsqVetoThresh );
        break;

      case 'R':
        rsqMaxSnr = (REAL4) atof( optarg );
        if ( rsqMaxSnr < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "threshold must be >= 0: "
              "(%f specified)\n",
              long_options[option_index].name, rsqMaxSnr );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", rsqMaxSnr );
        break;

      case 'p':
        rsqAboveSnrCoeff = (REAL4) atof( optarg );
        if ( rsqAboveSnrCoeff < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "coefficient must be >= 0: "
              "(%f specified)\n",
              long_options[option_index].name, rsqAboveSnrCoeff );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", rsqAboveSnrCoeff );
        break;

      case 'P':
        rsqAboveSnrPow = (REAL4) atof( optarg );
        if ( rsqAboveSnrPow < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "power must be >= 0: "
              "(%f specified)\n",
              long_options[option_index].name, rsqAboveSnrPow );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", rsqAboveSnrPow );
        break;

      case 'C':
        /* choose the clustering algorithm */
        {
          if ( ! strcmp( "nullstat", optarg) )
          {
            clusterchoice = nullstat;
          }
          else if ( ! strcmp( "cohsnr", optarg) )
          {
            clusterchoice = cohsnr;
          }
          else if ( ! strcmp( "effCohSnr", optarg) )
          {
            clusterchoice = effCohSnr;
          }
          else if ( ! strcmp( "snrByNullstat", optarg) )
          {
            clusterchoice = snrByNullstat;
          }
          else if ( ! strcmp( "autoCorrCohSqByNullstat", optarg) )
          {
            clusterchoice = autoCorrCohSqByNullstat;
          }
          else if ( ! strcmp( "crossCorrCohSqByNullstat", optarg) )
          {
            clusterchoice = autoCorrCohSqByNullstat;
          }
          else if ( ! strcmp( "autoCorrNullSqByNullstat", optarg) )
          {
            clusterchoice = autoCorrCohSqByNullstat;
          }
          else if ( ! strcmp( "crossCorrNullSqByNullstat", optarg) )
          {
            clusterchoice = crossCorrCohSqByNullstat;
          }
          else
          {
            fprintf( stderr, "invalid argument to  --%s:\n"
                "unknown clustering specified:\n "
                "%s (must be one of: cohsnr, effCohSnr, nullstat, snrByNullstat, autoCorrCohSqByNullstat, \n"
                "crossCorrCohSqByNullstat, autoCorrNullSqByNullstat, or crossCorrNullSqByNullstat)\n",
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
              "cluster window must be > 0: "
              "(%" LAL_INT8_FORMAT " specified)\n",
              long_options[option_index].name, cluster_dt );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%" LAL_INT8_FORMAT, cluster_dt );
        /* convert cluster time from ms to ns */
        cluster_dt *= 1000000LL;
        break;

      case 'v':
        /* create storage for the injection file name */
        optarg_len = strlen( optarg ) + 1;
        vetoFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( vetoFileName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'I':
        /* create storage for the injection file name */
        optarg_len = strlen( optarg ) + 1;
        injectFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( injectFileName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'd':
        optarg_len = strlen( optarg ) + 1;
        ifoName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( ifoName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'D':
        /* keep only coincs found in ifos */
        optarg_len = strlen( optarg ) + 1;
        ifos = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( ifos, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'T':
        /* injection coincidence time is specified on command line in ms */
        injectWindowNS = (INT8) atoi( optarg );
        if ( injectWindowNS < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "injection coincidence window must be >= 0: "
              "(%" LAL_INT8_FORMAT " specified)\n",
              long_options[option_index].name, injectWindowNS );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%" LAL_INT8_FORMAT, injectWindowNS );
        /* convert inject time from ms to ns */
        injectWindowNS *= 1000000LL;
        break;

      case 'm':
        /* create storage for the missed injection file name */
        optarg_len = strlen( optarg ) + 1;
        missedFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( missedFileName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case '?':
        exit( 1 );
        break;

      default:
        fprintf( stderr, "unknown error while parsing options\n" );
        exit( 1 );
    }
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
    snprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
  }
  else
  {
    snprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
  }

  /* check that the output file name has been specified */
  if ( ! outputFileName )
  {
    fprintf( stderr, "--output must be specified\n" );
    exit( 1 );
  }

  /* check that Data Type has been specified */
  if ( dataType == unspecified_data_type )
  {
    fprintf( stderr, "Error: --data-type must be specified\n");
    exit(1);
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

  /* check that if the rsq veto is being preformed,
                         we have the required options */
  if ( ( (rsqVetoThresh > 0) || (rsqMaxSnr > 0) ) && ( (rsqVetoThresh < 0)
    || (rsqMaxSnr < 0) ) )
  {
    fprintf( stderr, "--rsq-threshold and --rsq-max-snr and must be "
      "specified together" );
    exit( 1 );
  }
  else if ( (rsqAboveSnrCoeff > 0) && ( (rsqMaxSnr < 0) || (rsqVetoThresh < 0)
    || (rsqAboveSnrPow < 0) ) )
  {
    fprintf( stderr, "--rsq-max-snr --rsq-threshold and --rsq-power "
      "must be specified if --rsq-coeff is given\n" );
    exit( 1 );
  }
  else if ( (rsqAboveSnrPow > 0) && ( (rsqMaxSnr < 0) || (rsqVetoThresh < 0)
    || (rsqAboveSnrCoeff < 0) ) )
  {
    fprintf( stderr, "--rsq-max-snr --rsq-threshold and --rsq-coeff "
      "must be specified if --rsq-power is given\n" );
    exit( 1 );
  }

  /* check that we have all the options to do injections */
  if ( injectFileName && injectWindowNS < 0 )
  {
    fprintf( stderr, "--injection-coincidence must be specified if "
        "--injection-file is given\n" );
    exit( 1 );
  }
  else if ( ! injectFileName && injectWindowNS >= 0 )
  {
    fprintf( stderr, "--injection-file must be specified if "
        "--injection-coincidence is given\n" );
    exit( 1 );
  }

  if ( numSlides && extractSlide )
  {
    fprintf( stderr, "--num-slides and --extract-slide both specified\n"
        "this doesn't make sense\n" );
    exit( 1 );
  }

  /* save the sort triggers flag */
  if ( sortTriggers )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s",
        PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--sort-triggers" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, " " );
  }

  /* read in the veto file (if specified */

  if ( vetoFileName )
  {
    XLALSegListInit( &vetoSegs );
    LAL_CALL( LALSegListRead( &status, &vetoSegs, vetoFileName, NULL ),
        &status );
    XLALSegListCoalesce( &vetoSegs );
  }


  /*
   *
   * read in the input triggers from the xml files
   *
   */


  /* if we have run out of arguments on the command line, throw an error */
  if ( ! (optind < argc) )
  {
    fprintf( stderr, "Error: No input trigger files specified.\n" );
    exit( 1 );
  }

  /* read in the triggers */
  for( j = optind; j < argc; ++j )
  {
    INT4 numFileTriggers = 0;
    MultiInspiralTable   *inspiralFileList = NULL;
    MultiInspiralTable   *thisFileTrigger  = NULL;

    numInFiles++;

    numFileTriggers = XLALReadMultiInspiralTriggerFile( &inspiralFileList,
        &thisFileTrigger, &searchSummList, &inputFiles, argv[j] );
    numEvents += numFileTriggers;

    if (numFileTriggers < 0)
    {
      fprintf(stderr, "Error reading triggers from file %s\n",
          argv[j]);
      exit( 1 );
    }
    else
    {
      if ( vrbflg )
      {
        fprintf(stdout, "Read %d reading triggers from file %s\n",
            numFileTriggers, argv[j]);
      }
    }

    /* read the summ value table as well. */
    XLALReadSummValueFile(&summValueList, argv[j]);

    /*
     *
     *  keep only relevant triggers
     *
     */

    if( ifos )
    {
      numFileTriggers = XLALMultiInspiralIfosCut( &inspiralFileList, ifos );
      if ( vrbflg ) fprintf( stdout,
          "Kept %d coincs from %s instruments\n", numFileTriggers, ifos );
      numEventsInIfos += numFileTriggers;
    }

    /* Do playground_only or exclude_play cut */
    if ( dataType != all_data )
    {
      inspiralFileList = XLALPlayTestMultiInspiral( inspiralFileList,
          &dataType );
      /* count the triggers */
      numFileTriggers = XLALCountMultiInspiralTable( inspiralFileList );

      if ( dataType == playground_only && vrbflg ) fprintf( stdout,
          "Have %d playground triggers\n", numFileTriggers );
      else if ( dataType == exclude_play && vrbflg ) fprintf( stdout,
          "Have %d non-playground triggers\n", numFileTriggers );
    }
    numEventsKept += numFileTriggers;

    /*  Do snr cut */
    if ( snrStar > 0 )
    {
      inspiralFileList = XLALSNRCutMultiInspiral( inspiralFileList,
          snrStar );
      /* count the triggers  */
      numFileTriggers = XLALCountMultiInspiral( inspiralFileList );

      if ( vrbflg ) fprintf( stdout, "Have %d triggers after snr cut\n",
          numFileTriggers );
      numEventsAboveSNRThresh += numFileTriggers;
    }

    /* NOTE: Add vetoing:
       if ( vetoFileName )
       {
       inspiralFileList = XLALVetoMultiInspiral( inspiralFileList, &vetoSegs , ifoName);
       count the triggers
       numFileTriggers = XLALCountMultiInspiral( inspiralFileList );
       if ( vrbflg ) fprintf( stdout, "Have %d triggers after applying veto\n",
       numFileTriggers );
       numEventsSurvivingVeto += numFileTriggers;

       }
     */

    /* If there are any remaining triggers ... */
    if ( inspiralFileList )
    {
      /* add inspirals to list */
      if ( thisInspiralTrigger )
      {
        thisInspiralTrigger->next = inspiralFileList;
      }
      else
      {
        inspiralEventList = thisInspiralTrigger = inspiralFileList;
      }
      for( ; thisInspiralTrigger->next;
          thisInspiralTrigger = thisInspiralTrigger->next);
    }
  }

  for ( thisSearchSumm = searchSummList; thisSearchSumm;
      thisSearchSumm = thisSearchSumm->next )
  {
    UINT8 outPlayNS, outStartNS, outEndNS, triggerTimeNS;
    LIGOTimeGPS inPlay, outPlay;
    outStartNS = XLALGPSToINT8NS( &(thisSearchSumm->out_start_time) );
    outEndNS = XLALGPSToINT8NS( &(thisSearchSumm->out_end_time) );
    triggerTimeNS = outEndNS - outStartNS;

    /* check for events and playground */
    if ( dataType != all_data )
    {
      XLALPlaygroundInSearchSummary( thisSearchSumm, &inPlay, &outPlay );
      outPlayNS = XLALGPSToINT8NS( &outPlay );

      if ( dataType == playground_only )
      {
        /* increment the total trigger time by the amount of playground */
        triggerInputTimeNS += outPlayNS;
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
  }


  /*
   *
   * sort the inspiral events by time
   *
   */


  if ( injectFileName || sortTriggers )
  {
    inspiralEventList = XLALSortMultiInspiral( inspiralEventList,
        *LALCompareMultiInspiralByTime );
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

    /* keep play/non-play/all injections */
    if ( dataType == playground_only && vrbflg ) fprintf( stdout,
        "Keeping only playground injections\n" );
    else if ( dataType == exclude_play && vrbflg ) fprintf( stdout,
        "Keeping only non-playground injections\n" );
    else if ( dataType == all_data && vrbflg ) fprintf( stdout,
        "Keeping all injections\n" );
    XLALPlayTestSimInspiral( &simEventHead, &dataType );

    /* keep only injections in times analyzed */
    numSimInData = XLALSimInspiralInSearchedData( &simEventHead,
        &searchSummList );

    if ( vrbflg ) fprintf( stdout, "%d injections in analyzed data\n",
        numSimInData );


    /* check for events that are coincident with injections */
    numSimFound = XLALMultiSimInspiralTest( &simEventHead,
        &inspiralEventList, &missedSimHead, &missedHead, injectWindowNS );

    if ( vrbflg ) fprintf( stdout, "%d injections found in the ifos\n",
        numSimFound );

    if ( numSimFound )
    {
      for ( thisEvent = inspiralEventList; thisEvent;
          thisEvent = thisEvent->next, numMultiFound++ );
      if ( vrbflg ) fprintf( stdout,
          "%d triggers found at times of injection\n", numMultiFound );
    }

    /* free the missed singles  */
    while ( missedHead )
    {
      thisEvent = missedHead;
      missedHead = missedHead->next;
      XLALFreeMultiInspiral( &thisEvent );
    }
  }


  /*
   *
   * extract specified slide
   *
   */

  if ( extractSlide )
  {
    slideEvent = XLALMultiInspiralSlideCut( &inspiralEventList, extractSlide );
    /* free events from other slides */
    while ( inspiralEventList )
    {
      thisEvent = inspiralEventList;
      inspiralEventList = inspiralEventList->next;
      XLALFreeMultiInspiral( &thisEvent );
    }

    /* move events to inspiralEventList */
    inspiralEventList = slideEvent;
    slideEvent = NULL;
  }


  /*
   *
   * cluster the remaining events
   *
   */


  if ( inspiralEventList && clusterchoice )
  {
    if ( vrbflg ) fprintf( stdout, "clustering remaining triggers... " );

    if ( !numSlides ) {
      numClusteredEvents = XLALClusterMultiInspiralTable( &inspiralEventList,
        cluster_dt, clusterchoice );
    }
    else
    {
      int slide = 0;
      int numClusteredSlide = 0;
      MultiInspiralTable *tmp_slideEvent = NULL;
      MultiInspiralTable *slideClust = NULL;

      if ( vrbflg ) fprintf( stdout, "splitting events by slide\n" );

      for( slide = -numSlides; slide < (numSlides + 1); slide++)
      {
        if ( vrbflg ) fprintf( stdout, "slide number %d; ", slide );
        /* extract the slide */
        tmp_slideEvent = XLALMultiInspiralSlideCut( &inspiralEventList, slide );
        /* run clustering */
        numClusteredSlide = XLALClusterMultiInspiralTable( &tmp_slideEvent,
          cluster_dt, clusterchoice);

        if ( vrbflg ) fprintf( stdout, "%d clustered events \n",
          numClusteredSlide );
        numClusteredEvents += numClusteredSlide;

        /* add clustered triggers */
        if( tmp_slideEvent )
        {
          if( slideClust )
          {
            thisEvent = thisEvent->next = tmp_slideEvent;
          }
          else
          {
            slideClust = thisEvent = tmp_slideEvent;
          }
          /* scroll to end of list */
          for( ; thisEvent->next; thisEvent = thisEvent->next);
        }
      }

      /* free inspiralEventList -- although we expect it to be empty */
      while ( inspiralEventList )
      {
        thisEvent = inspiralEventList;
        inspiralEventList = inspiralEventList->next;
        XLALFreeMultiInspiral( &thisEvent );
      }

      /* move events to coincHead */
      inspiralEventList = slideClust;
      slideClust = NULL;
    }

    if ( vrbflg ) fprintf( stdout, "done\n" );
    if ( vrbflg ) fprintf( stdout, "%d clustered events \n",
        numClusteredEvents );
  }


  /*
   *
   * update search_summary->nevents with an authoritative count of triggers
   *
   */

  searchSummList->nevents = 0;
  thisEvent = inspiralEventList;
  while (thisEvent) {
    searchSummList->nevents += 1;
    thisEvent = thisEvent->next;
  }

  /*
   *
   * write output data
   *
   */


  /* write the main output file containing found injections */
  if ( vrbflg ) fprintf( stdout, "writing output xml files... " );
  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlStream, outputFileName ), &status );

  /* write out the process and process params tables */
  if ( vrbflg ) fprintf( stdout, "process... " );
  XLALGPSTimeNow(&(proctable.processTable->end_time));
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, process_table ),
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, proctable,
        process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );
  free( proctable.processTable );

  /* erase the first empty process params entry */
  {
    ProcessParamsTable *emptyPPtable = procparams.processParamsTable;
    procparams.processParamsTable = procparams.processParamsTable->next;
    free( emptyPPtable );
  }

  /* write the process params table */
  if ( vrbflg ) fprintf( stdout, "process_params... " );
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream,
        process_params_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, procparams,
        process_params_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );

  /* write search_summary table */
  if ( vrbflg ) fprintf( stdout, "search_summary... " );
  outputTable.searchSummaryTable = searchSummList;
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream,
        search_summary_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable,
        search_summary_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );

  /* write the search_summvars table */
  if ( vrbflg ) fprintf( stdout, "search_summvars... " );
  LAL_CALL( LALBeginLIGOLwXMLTable( &status ,&xmlStream,
        search_summvars_table), &status );
  searchSummvarsTable.searchSummvarsTable = inputFiles;
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, searchSummvarsTable,
        search_summvars_table), &status );
  LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream), &status );

  /* write summ_value table */
  if ( summValueList )
  {
    if ( vrbflg ) fprintf( stdout, "search_summary... " );
    outputTable.summValueTable = summValueList;
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream,
          summ_value_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable,
          summ_value_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );
  }

  /* Write the found injections to the sim table */
  if ( simEventHead )
  {
    if ( vrbflg ) fprintf( stdout, "sim_inspiral... " );
    outputTable.simInspiralTable = simEventHead;
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream,
          sim_inspiral_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable,
          sim_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream ), &status );
  }

  /* Write the results to the inspiral table */
  if ( inspiralEventList )
  {
    if ( vrbflg ) fprintf( stdout, "multi_inspiral... " );
    outputTable.multiInspiralTable = inspiralEventList;
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream,
          multi_inspiral_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable,
          multi_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream ), &status);
  }

  /* close the output file */
  LAL_CALL( LALCloseLIGOLwXMLFile(&status, &xmlStream), &status);
  if ( vrbflg ) fprintf( stdout, "done\n" );


  if ( missedFileName )
  {
    /* open the missed injections file and write the missed injections to it */
    if ( vrbflg ) fprintf( stdout, "writing missed injections... " );
    memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
    LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlStream, missedFileName ),
        &status );

    if ( missedSimHead )
    {
      outputTable.simInspiralTable = missedSimHead;
      LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, sim_inspiral_table ),
          &status );
      LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable,
            sim_inspiral_table ), &status );
      LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream ), &status );
    }

    LAL_CALL( LALCloseLIGOLwXMLFile( &status, &xmlStream ), &status );
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
    fprintf( fp, "number of triggers in input data %d \n", numEventsKept );
    if ( ifoName )
    {
      fprintf( fp, "number of triggers from %s ifo %d \n", ifoName,
          numEventsInIFO );
    }


    if ( snrStar > 0 )
    {
      fprintf( fp, "number of triggers in input data with snr above %f: %d \n",
          snrStar, numEventsAboveSNRThresh );
    }

    if ( rsqVetoThresh > 0 )
    {
      fprintf( fp, "performed R-squared veto on triggers with snr < %f\n",
          rsqMaxSnr);
      fprintf( fp, "with rsqveto_duration below %f\n",
          rsqVetoThresh);
      if ( (rsqAboveSnrCoeff > 0) && (rsqAboveSnrPow > 0) )
      {
        fprintf( fp, "and on triggers with snr > %f\n",
            rsqMaxSnr);
        fprintf( fp, "with rsqveto_duration above %f * snr ^ %f\n",
            rsqAboveSnrCoeff, rsqAboveSnrPow );
      }
      fprintf( fp, "the number of triggers below the R-squared veto are: %d \n",
          numEventsBelowRsqThresh);
    }

    if ( vetoFileName )
    {
      fprintf( fp, "number of triggers not vetoed by %s: %d \n",
          vetoFileName, numEventsSurvivingVeto );
    }

    XLALINT8NSToGPS( &triggerTime, triggerInputTimeNS );
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
          (injectWindowNS / 1000000LL), numMultiFound );

      fprintf( fp, "efficiency: %f \n",
          (REAL4) numSimFound / (REAL4) numSimInData );
    }

    if ( extractSlide )
    {
      fprintf( fp, "kept only triggers from slide %d\n", extractSlide );
    }

    if ( clusterchoice )
    {
      if ( numSlides )
      {
        fprintf( fp, "clustering triggers from %d slides separately\n",
            numSlides );
      }
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
  while ( inspiralEventList )
  {
    thisEvent = inspiralEventList;
    inspiralEventList = inspiralEventList->next;
    LAL_CALL ( LALFreeMultiInspiral ( &status, &thisEvent ), &status);
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

  /* free search summaries read in */
  while ( searchSummList )
  {
    thisSearchSumm = searchSummList;
    searchSummList = searchSummList->next;
    LALFree( thisSearchSumm );
  }

  while ( summValueList )
  {
    SummValueTable *thisSummValue;
    thisSummValue = summValueList;
    summValueList = summValueList->next;
    LALFree( thisSummValue );
  }

  if ( vetoFileName )
  {
    XLALSegListClear( &vetoSegs );
  }


  if ( vrbflg ) fprintf( stdout, "checking memory leaks and exiting\n" );
  LALCheckMemoryLeaks();
  exit( 0 );
}
