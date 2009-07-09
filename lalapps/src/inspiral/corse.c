/*----------------------------------------------------------------------- 
 * 
 * File Name: corse.c
 *
 * Author: Keppel, D
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
#include <lal/lalGitID.h>
#include <lalappsGitID.h>
#include <lalapps.h>
#include <processtable.h>

RCSID("$Id$");

#define PROGRAM_NAME "corse"
#define CVS_ID_STRING "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define CVS_NAME_STRING "$Name$"

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
      " [--help]                       display this message\n"\
      " [--verbose]                    print progress information\n"\
      " [--version]                    print version information and exit\n"\
      " [--debug-level]       level    set the LAL debug level to LEVEL\n"\
      " [--user-tag]          usertag  set the process_params usertag\n"\
      " [--comment]           string   set the process table comment\n"\
      " [--mass-tag]          string   specify what mass bin in summfile\n"\
      "\n"\
      " [--glob-zero]         glob     use pattern glob to determine the input files\n"\
      " [--input-zero]        input    read list of input XML files from input\n"\
      "\n"\
      " [--glob-slide]        glob     use pattern glob to determine the input files\n"\
      " [--input-slide]       input    read list of input XML files from input\n"\
      "\n"\
      "  --output             output   write output data to file: output\n"\
      " [--summary-file]      summ     write trigger analysis summary to summ\n"\
      " [--loudest]           file     write file only containing loudest triggers\n"\
      "\n"\
      "  --data-type          datatype data type of zero-lag, must be one of\n"\
      "                                (playground_only|exclude_play|all_data)\n"\
      "\n"\
      " [--coinc-cut]         ifos     only keep triggers from IFOS\n"\
      "\n"\
      "  --num-slides         slides   number of time slides performed \n"\
      "                                (slides * time = bkg time analyzed)\n"\
      "  --time-analyzed-file file     file containing the amount of time\n"\
      "                                analyzed per time slide\n"\
      "                                for the background estimation\n"\
      "                                (used in rate calculation)\n"\
      " [--background-modifier]\n"\
      "                       mod      the background modifier divides\n"\
      "                                the background time analyzed by mod\n"\
      " [--fit-num]           n        use an exponential fit to n loudest\n"\
      "                                background triggers to calculate\n"\
      "                                the FAR for triggers in tail\n"\
      " [--sort-triggers]              time sort the coincident triggers\n"\
      "  --coinc-stat         stat     use coinc statistic for cluster/cut\n"\
      "                       (snrsq|effective_snrsq|s3_snr_chi_stat|bitten_l)\n"\
      " [--stat-threshold]    thresh   discard all triggers with stat less than thresh\n"\
      " [--rate-threshold]    rate     discard all triggers with rate greater than thresh\n"\
      " [--eff-snr-denom-fac] number   parameter for clustering effective snr denominator (traditionally 250) \n"\
      " [--h1-bittenl-a]      bitten   paramater a for clustering\n"\
      " [--h1-bittenl-b]      bitten   paramater b for clustering\n"\
      " [--h2-bittenl-a]      bitten   paramater a for clustering\n"\
      " [--h2-bittenl-b]      bitten   paramater b for clustering\n"\
      " [--l1-bittenl-a]      bitten   paramater a for clustering\n"\
      " [--l1-bittenl-b]      bitten   paramater b for clustering\n"\
      "\n"\
      " [--injection-file]    inj_file read injection parameters from inj_file\n"\
      " [--injection-window]  inj_win  trigger and injection coincidence window (ms)\n"\
      " [--missed-injections] missed   write missed injections to file missed\n"\
      "\n");
}

/* function to read the next line of data from the input file list */
static char *get_next_line( char *line, size_t size, FILE *fp )
{
  char *s;
  do
    s = fgets( line, size, fp );
  while ( ( line[0] == '#' || line[0] == '%' ) && s );
  return s;
}

int sortTriggers = 0;
LALPlaygroundDataMask dataType;
extern int vrbflg;

int main( int argc, char *argv[] )
{
  /* lal initialization variables */
  LALLeapSecAccuracy accuracy = LALLEAPSEC_LOOSE;
  LALStatus status = blank_status ;

  /*  program option variables */
  CHAR *userTag = NULL;
  CHAR comment[LIGOMETA_COMMENT_MAX];
  CHAR *massTag = NULL;
  char *ifos = NULL;
  char *inputGlobZero = NULL;
  char *inputFileNameZero = NULL;
  char *inputGlobSlide = NULL;
  char *inputFileNameSlide = NULL;
  char *outputFileName = NULL;
  char *timeAnalyzedFileName = NULL;
  char *summFileName = NULL;
  char *loudestFileName = NULL;
  CoincInspiralStatistic coincstat = no_stat;
  int   fitNum = 0;
  REAL4 fitStat = -1;
  REAL4 fitA = 0;
  REAL4 fitB = 0;
  REAL4 statThreshold = -1;
  REAL4 rateThreshold = -1;
  REAL4 loudestRate = 0;
  REAL4 timeAnalyzed = 1;
  REAL4 timeModifier = 1;
  REAL4 bkgtimeAnalyzed = 0;
  char *coincifos = NULL;
  char *injectFileName = NULL;
  INT8 injectWindowNS = -1;
  char *missedFileName = NULL;
  int j;
  FILE *fp = NULL;
  glob_t globbedZeroFiles;
  glob_t globbedSlideFiles;
  int numInZeroFiles = 0;
  int numInSlideFiles = 0;
  char **inZeroFileNameList;
  char **inSlideFileNameList;
  char line[MAX_PATH];

  UINT8 triggerInputTimeNS = 0;

  MetadataTable         proctable;
  MetadataTable         procparams;
  ProcessParamsTable   *this_proc_param;

  int                   numSimEvents = 0;
  int                   numSimInData = 0;

  SearchSummvarsTable  *inputZeroFiles = NULL;
  SearchSummvarsTable  *inputSlideFiles = NULL;
  SearchSummvarsTable  *thisInputFile = NULL;
  SearchSummaryTable   *searchSummList = NULL;
  SearchSummaryTable   *searchSummSlideList = NULL;
  SearchSummaryTable   *thisSearchSumm = NULL;
  SummValueTable       *summValueList = NULL;
  SimInspiralTable     *simEventHead = NULL;
  SimInspiralTable     *thisSimEvent = NULL;
  SimInspiralTable     *missedSimHead = NULL;
  SimInspiralTable     *missedSimCoincHead = NULL;
  SimInspiralTable     *tmpSimEvent = NULL;

  int                   numSlides = -1;
  int                   numZeroTriggers = 0;
  int                   numSlideTriggers = 0;
  int                   numZeroCoincs = 0;
  int                   numSlideCoincs = 0;
  int                   numEventsPlayTest = 0;
  int                   numEventsAboveThresh = 0;
  int                   numEventsCoinc = 0;
  int                   numSnglFound = 0;
  int                   numCoincFound = 0;

  SnglInspiralTable    *inspiralZeroEventList = NULL;
  SnglInspiralTable    *inspiralSlideEventList = NULL;
  SnglInspiralTable    *thisSngl = NULL;
  SnglInspiralTable    *missedSnglHead = NULL;
  SnglInspiralTable    *thisInspiralTrigger = NULL;
  SnglInspiralTable    *snglOutput = NULL;

  CoincInspiralTable   *coincZeroHead = NULL;
  CoincInspiralTable   *coincSlideHead = NULL;
  CoincInspiralTable   *thisCoinc = NULL;
  CoincInspiralTable   *missedCoincHead = NULL;

  CoincInspiralStatParams    bittenLParams;

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
  LAL_CALL(LALGPSTimeNow(&status, &(proctable.processTable->start_time), 
        &accuracy), &status);
  if (strcmp(CVS_REVISION,"$Revi" "sion$"))
    {
      LAL_CALL( populate_process_table( &status, proctable.processTable, 
					PROGRAM_NAME, CVS_REVISION,
					CVS_SOURCE, CVS_DATE ), &status );
    }
  else
    {
      LAL_CALL( populate_process_table( &status, proctable.processTable, 
					PROGRAM_NAME, lalappsGitCommitID,
					lalappsGitGitStatus,
					lalappsGitCommitDate ), &status );
    }
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
    calloc( 1, sizeof(ProcessParamsTable) );
  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );

  memset( &bittenLParams, 0, sizeof(CoincInspiralStatParams   ) );
  /* Assign a default effective snr denominator factor */
  bittenLParams.eff_snr_denom_fac = 250.0;

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
      {"verbose",                 no_argument,       &vrbflg,              1 },
      {"sort-triggers",           no_argument,  &sortTriggers,             1 },
      {"help",                    no_argument,            0,              'h'},
      {"debug-level",             required_argument,      0,              'z'},
      {"user-tag",                required_argument,      0,              'Z'},
      {"userTag",                 required_argument,      0,              'Z'},
      {"comment",                 required_argument,      0,              'c'},
      {"mass-tag",                required_argument,      0,              'M'},
      {"version",                 no_argument,            0,              'V'},
      {"data-type",               required_argument,      0,              'k'},
      {"glob-zero",               required_argument,      0,              'g'},
      {"glob-slide",              required_argument,      0,              'G'},
      {"input-zero",              required_argument,      0,              'i'},
      {"input-slide",             required_argument,      0,              'I'},
      {"output",                  required_argument,      0,              'o'},
      {"summary-file",            required_argument,      0,              'S'},
      {"loudest",                 required_argument,      0,              'L'},
      {"num-slides",              required_argument,      0,              'N'},
      {"background-modifier",     required_argument,      0,              't'},
      {"fit-num",                 required_argument,      0,              's'},
      {"time-analyzed-file",      required_argument,      0,              'A'},
      {"coinc-stat",              required_argument,      0,              'C'},
      {"stat-threshold",          required_argument,      0,              'E'},
      {"rate-threshold",          required_argument,      0,              'R'},
      {"coinc-cut",               required_argument,      0,              'D'},
      {"injection-file",          required_argument,      0,              'f'},
      {"injection-window",        required_argument,      0,              'T'},
      {"eff-snr-denom-fac",    required_argument,      0,              'a'},
      {"h1-bittenl-a",            required_argument,      0,              'q'},
      {"h1-bittenl-b",            required_argument,      0,              'r'},
      {"h2-bittenl-a",            required_argument,      0,              'j'},
      {"h2-bittenl-b",            required_argument,      0,              'n'},
      {"l1-bittenl-a",            required_argument,      0,              'l'},
      {"l1-bittenl-b",            required_argument,      0,              'p'},
      {"missed-injections",       required_argument,      0,              'm'},

      {0, 0, 0, 0}
    };
    int c;

    /* getopt_long stores the option index here. */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only ( argc, argv, "a:b:c:f:g:hi:j:k:l:m:n:o:p:q:r:s:t:z:"
                                       "A:C:D:E:G:I:L:N:R:S:T:VZ",
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

     case 'a':
        bittenLParams.eff_snr_denom_fac = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg);
        break;
 
      case 'q':
        bittenLParams.param_a[LAL_IFO_H1] = atof(optarg);
      ADD_PROCESS_PARAM( "float", "%s", optarg);
        break;
      
      case 'r':
        bittenLParams.param_b[LAL_IFO_H1] = atof(optarg);
      ADD_PROCESS_PARAM( "float", "%s", optarg);
        break;

      case 'j':
        bittenLParams.param_a[LAL_IFO_H2] = atof(optarg);
      ADD_PROCESS_PARAM( "float", "%s", optarg);
        break;

      case 'n':
        bittenLParams.param_b[LAL_IFO_H2] = atof(optarg);
      ADD_PROCESS_PARAM( "float", "%s", optarg);
        break;

      case 'l':
        bittenLParams.param_a[LAL_IFO_L1] = atof(optarg);
      ADD_PROCESS_PARAM( "float", "%s", optarg);
        break;

      case 'p':
        bittenLParams.param_b[LAL_IFO_L1] = atof(optarg);
      ADD_PROCESS_PARAM( "float", "%s", optarg);
        break;

      case 'h':
        print_usage(argv[0]);
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

      case 'M':
        /* create storage for and save massTag string; this is not written to
         * the process params table; just used for summfile */
        massTag = (CHAR *) calloc( strlen(optarg), sizeof(CHAR) );
        memcpy( massTag, optarg, strlen(optarg) );
        break;

      case 'V':
        fprintf( stdout, "COincidence Rate-Statistic Estimator\n"
            "Drew Keppel\n"
            "CVS Version: " CVS_ID_STRING "\n" 
            "CVS Tag: " CVS_NAME_STRING "\n" );
	fprintf( stdout, lalappsGitID );
        exit( 0 );
        break;

      case 'g':
        /* create storage for the input file glob */
        optarg_len = strlen( optarg ) + 1;
        inputGlobZero = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( inputGlobZero, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "'%s'", optarg );
        break;

      case 'G':
        /* create storage for the input file glob */
        optarg_len = strlen( optarg ) + 1;
        inputGlobSlide = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( inputGlobSlide, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "'%s'", optarg );
        break;

      case 'i':
        /* create storage for the input file name */
        optarg_len = strlen( optarg ) + 1;
        inputFileNameZero = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( inputFileNameZero, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'I':
        /* create storage for the input file name */
        optarg_len = strlen( optarg ) + 1;
        inputFileNameSlide = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( inputFileNameSlide, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'o':
        /* create storage for the output file name */
        optarg_len = strlen( optarg ) + 1;
        outputFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( outputFileName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 't':
        /* store the time modifier */
        timeModifier = atof( optarg );
        if ( timeModifier <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "timeModifier > 0: "
              "(%e specified)\n",
              long_options[option_index].name, timeModifier );
          exit( 1 );    
        }
        ADD_PROCESS_PARAM( "float", "%e", timeModifier );
        break;

      case 's':
        /* store the number to fit above */
        fitNum = atoi( optarg );
        if ( fitNum <= 1 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "fitNum > 1: "
              "(%d specified)\n",
              long_options[option_index].name, fitNum );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", fitNum );
        break;

      case 'A':
        /* create storage for the output file name */
        optarg_len = strlen( optarg ) + 1;
        timeAnalyzedFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( timeAnalyzedFileName, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'N':
        /* store the number of slides */
        numSlides = atoi( optarg );
        if ( numSlides < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "numSlides > 0: "
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

      case 'L':
        /* create storage for the loudest file name */
        optarg_len = strlen( optarg ) + 1;
        loudestFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( loudestFileName, optarg, optarg_len );
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

      case 'C':
        /* choose the coinc statistic */
        {        
          if ( ! strcmp( "snrsq", optarg ) )
          {
            coincstat = snrsq;
          }
          else if ( ! strcmp( "bitten_l", optarg ) )
          {
            coincstat = bitten_l;
          }
          else if ( ! strcmp( "bitten_lsq", optarg ) )
          {
            coincstat = bitten_lsq;
          }
          else if ( ! strcmp( "effective_snrsq", optarg) )
          {
            coincstat = effective_snrsq;
          }
          else
          {
            fprintf( stderr, "invalid argument to  --%s:\n"
                "unknown coinc statistic:\n "
                "%s (must be one of:\n"
                "snrsq, effective_snrsq, bitten_l)\n",
                long_options[option_index].name, optarg);
            exit( 1 );
          }
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        break;

      case 'E':
        /* store the stat threshold for a cut */
        statThreshold = atof( optarg );
        if ( statThreshold < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "statThreshold must be >= 0: (%f specified)\n",
              long_options[option_index].name, statThreshold );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%f", statThreshold );
        break;

      case 'R':        /* store the rate threshold for a cut */
        rateThreshold = atof( optarg );
        if ( rateThreshold < 0 ) 
        {
          fprintf( stdout, "invalid argument to --%s:\n"  
              "rateThreshold must be <= 0: (%f specified)\n",
              long_options[option_index].name, rateThreshold );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%f", rateThreshold );
        break;
 
      case 'f':
        /* create storage for the injection file name */
        optarg_len = strlen( optarg ) + 1;
        injectFileName = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( injectFileName, optarg, optarg_len );
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
              "(%ld specified)\n",
              long_options[option_index].name, injectWindowNS );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%lld", injectWindowNS );
        /* convert inject time from ms to ns */
        injectWindowNS *= LAL_INT8_C(1000000);
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
  if ( (! inputGlobZero && ! inputFileNameZero) ||
       (inputGlobZero && inputFileNameZero) )
  {
    fprintf( stderr, "exactly one of --glob-zero or --input-zero must be "
        "specified\n" );
    exit( 1 );
  }
  if ( (! inputGlobSlide && ! inputFileNameSlide) ||
       (inputGlobSlide && inputFileNameSlide) )
  {
    fprintf( stderr, "exactly one of --glob-slide or --input-slide must be "
        "specified\n" );
    exit( 1 );
  }
  if ( ! outputFileName )
  {
    fprintf( stderr, "--output must be specified");
    exit(1);
  }

  /* check that Data Type has been specified */
  if ( dataType == unspecified_data_type )
  {
    fprintf( stderr, "Error: --data-type must be specified\n");
    exit(1);
  }

  /* check that if clustering is being done that we have all the options */
  if ( coincstat == no_stat )
  {
    fprintf( stderr, 
        "--coinc-stat must be specified\n" );
    exit( 1 );
  }

  /* check that we have all the options to do injections */
  if ( injectFileName && injectWindowNS < 0 )
  {
    fprintf( stderr, "--injection-window must be specified if "
        "--injection-file is given\n" );
    exit( 1 );
  }
  else if ( ! injectFileName && injectWindowNS >= 0 )
  {
    fprintf( stderr, "--injection-file must be specified if "
        "--injection-window is given\n" );
    exit( 1 );
  }

  if ( numSlides < 0 )
  {
    fprintf( stderr, "--num-slides must be specified\n" );
    exit( 1 );
  }

  if ( ! timeAnalyzedFileName )
  {
    fprintf( stderr, "--time-analyzed-file must be specified \n" );
    exit( 1 );
  }

  if ( ! timeModifier )
  {     
    fprintf( stderr, "--background-modifier must be non-zero\n" );
    exit( 1 );
  } 
 
  /* save the sort triggers flag */
  if ( sortTriggers )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *) 
      calloc( 1, sizeof(ProcessParamsTable) ); 
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s",
        PROGRAM_NAME ); 
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--sort-triggers" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" ); 
    LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, " " );
  }

  /*
   *
   * read in the input zero-lag triggers from the xml files
   *
   */


  if ( inputGlobZero )
  {
    /* use glob() to get a list of the input file names */
    if ( glob( inputGlobZero, GLOB_ERR, NULL, &globbedZeroFiles ) )
    {
      perror( "error:" );
      fprintf( stderr, "error globbing files from %s\n", inputGlobZero );
      exit( 1 );
    }

    numInZeroFiles = globbedZeroFiles.gl_pathc;
    inZeroFileNameList =
        (char **) LALCalloc( numInZeroFiles, sizeof(char *) );

    for ( j = 0; j < numInZeroFiles; ++j )
    {
      inZeroFileNameList[j] = globbedZeroFiles.gl_pathv[j];
    }
  }
  else if ( inputFileNameZero )
  {
    /* read the list of input filenames from a file */
    fp = fopen( inputFileNameZero, "r" );
    if ( ! fp )
    {
      perror( "error:" );
      fprintf( stderr,
          "could not open file containing list of zero-lag xml files\n" );
      exit( 1 );
    }

    /* count the number of lines in the file */
    while ( get_next_line( line, sizeof(line), fp ) )
    {
      ++numInZeroFiles;
    }
    rewind( fp );

    /* allocate memory to store the input file names */
    inZeroFileNameList =
        (char **) LALCalloc( numInZeroFiles, sizeof(char *) );

    /* read in the input file names */
    for ( j = 0; j < numInZeroFiles; ++j )
    {
      inZeroFileNameList[j] = (char *) LALCalloc( MAX_PATH, sizeof(char) );
      get_next_line( line, sizeof(line), fp );
      strncpy( inZeroFileNameList[j], line, strlen(line) - 1);
    }

    fclose( fp );
  }
  else
  {
    fprintf( stderr, "no zero-lag input file mechanism specified\n" );
    exit( 1 );
  }


  /*
   *
   * read in the input slide triggers from the xml files
   *
   */


  if ( inputGlobSlide )
  {
    /* use glob() to get a list of the input file names */
    if ( glob( inputGlobSlide, GLOB_ERR, NULL, &globbedSlideFiles ) )
    {
      perror( "error:" );
      fprintf( stderr, "error globbing files from %s\n", inputGlobSlide );
      exit( 1 );
    }

    numInSlideFiles = globbedSlideFiles.gl_pathc;
    inSlideFileNameList =
        (char **) LALCalloc( numInSlideFiles, sizeof(char *) );

    for ( j = 0; j < numInSlideFiles; ++j )
    {
      inSlideFileNameList[j] = globbedSlideFiles.gl_pathv[j];
    }
  }
  else if ( inputFileNameSlide )
  {
    /* read the list of input filenames from a file */
    fp = fopen( inputFileNameSlide, "r" );
    if ( ! fp )
    {
      perror( "error:" );
      fprintf( stderr,
          "could not open file containing list of slide xml files\n" );
      exit( 1 );
    }

    /* count the number of lines in the file */
    while ( get_next_line( line, sizeof(line), fp ) )
    {
      ++numInSlideFiles;
    }
    rewind( fp );

    /* allocate memory to store the input file names */
    inSlideFileNameList =
        (char **) LALCalloc( numInSlideFiles, sizeof(char *) );

    /* read in the input file names */
    for ( j = 0; j < numInSlideFiles; ++j )
    {
      inSlideFileNameList[j] = (char *) LALCalloc( MAX_PATH, sizeof(char) );
      get_next_line( line, sizeof(line), fp );
      strncpy( inSlideFileNameList[j], line, strlen(line) - 1);
    }

    fclose( fp );
  }
  else
  {
    fprintf( stderr, "no slide input file mechanism specified\n" );
    exit( 1 );
  }

  /* read in the zero-lag triggers */
  for( j = 0; j < numInZeroFiles; ++j )
  {
    INT4 numFileTriggers = 0;
    INT4 numFileCoincs   = 0;
    SnglInspiralTable   *inspiralFileList = NULL;
    SnglInspiralTable   *thisFileTrigger  = NULL;
    CoincInspiralTable  *coincFileHead    = NULL;
    
    numFileTriggers = XLALReadInspiralTriggerFile( &inspiralFileList,
        &thisFileTrigger, &searchSummList, &inputZeroFiles,
        inZeroFileNameList[j] );
    if (numFileTriggers < 0)
    {
      fprintf(stderr, "Error reading triggers from file %s\n",
          inZeroFileNameList[j]);
      exit( 1 );
    }
    else
    {
      if ( vrbflg )
      {
        fprintf(stdout, "Read %d reading triggers from file %s\n",
            numFileTriggers, inZeroFileNameList[j]);
      }
    }

    /* read the summ value table as well. */
    XLALReadSummValueFile(&summValueList, inZeroFileNameList[j]);

    
    /* reconstruct the coincs */
    numFileCoincs = XLALRecreateCoincFromSngls( &coincFileHead, 
        &inspiralFileList );
    if( numFileCoincs < 0 )
    {
      fprintf(stderr, 
          "Unable to reconstruct coincs from single ifo triggers");
      exit( 1 );
    }
    
    if ( vrbflg )
    {
      fprintf( stdout,
          "Recreated %d coincs from the %d triggers in file %s\n", 
          numFileCoincs, numFileTriggers, inZeroFileNameList[j] );
    }
    numZeroCoincs += numFileCoincs;

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
        inspiralZeroEventList = thisInspiralTrigger = inspiralFileList;
      }
      for( ; thisInspiralTrigger->next; 
          thisInspiralTrigger = thisInspiralTrigger->next);
      numZeroTriggers += numFileTriggers;
    }

    /* Do playground_only or exclude_play cut */
    if ( dataType != all_data )
    {
      coincFileHead = XLALPlayTestCoincInspiral( coincFileHead, 
          &dataType );
      /* count the triggers, scroll to end of list */
      numFileCoincs = XLALCountCoincInspiral( coincFileHead );

      if ( dataType == playground_only && vrbflg ) fprintf( stdout, 
        "Have %d playground triggers\n", numFileCoincs );
      else if ( dataType == exclude_play && vrbflg ) fprintf( stdout, 
        "Have %d non-playground triggers\n", numFileCoincs );
      numEventsPlayTest += numFileCoincs;
    }

    /* add coincs to list */
    if( numFileCoincs )
    {
      if ( thisCoinc )
      {
        thisCoinc->next = coincFileHead;
      }
      else
      {
        coincZeroHead = thisCoinc = coincFileHead;
      }
      for ( ; thisCoinc->next; thisCoinc = thisCoinc->next );
    }
  }

  thisCoinc = NULL;

  /* read in the slide triggers */
  for( j = 0; j < numInSlideFiles; ++j )
  {
    INT4 numFileTriggers = 0;
    INT4 numFileCoincs   = 0;
    SnglInspiralTable   *inspiralFileList = NULL;
    SnglInspiralTable   *thisFileTrigger  = NULL;
    CoincInspiralTable  *coincFileHead    = NULL;

    numFileTriggers = XLALReadInspiralTriggerFile( &inspiralFileList,
        &thisFileTrigger, &searchSummSlideList, &inputSlideFiles,
        inSlideFileNameList[j] );
    if (numFileTriggers < 0)
    {
      fprintf(stderr, "Error reading triggers from file %s\n",
          inSlideFileNameList[j]);
      exit( 1 );
    }
    else
    {
      if ( vrbflg )
      {
        fprintf(stdout, "Read %d reading triggers from file %s\n",
            numFileTriggers, inSlideFileNameList[j]);
      }
    }

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
        inspiralSlideEventList = thisInspiralTrigger = inspiralFileList;
      }
      for( ; thisInspiralTrigger->next;
          thisInspiralTrigger = thisInspiralTrigger->next);
      numSlideTriggers += numFileTriggers;
    }

    /* reconstruct the coincs */
    numFileCoincs = XLALRecreateCoincFromSngls( &coincFileHead,
        &inspiralFileList );
    if( numFileCoincs < 0 )
    {
      fprintf(stderr,
          "Unable to reconstruct coincs from single ifo triggers");
      exit( 1 );
    }

    if ( vrbflg )
    {
      fprintf( stdout,
          "Recreated %d coincs from the %d triggers in file %s\n",
          numFileCoincs, numFileTriggers, inSlideFileNameList[j] );
    }
    numSlideCoincs += numFileCoincs;

    /* add coincs to list */
    if( numFileCoincs )
    {
      if ( thisCoinc )
      {
        thisCoinc->next = coincFileHead;
      }
      else
      {
        coincSlideHead = thisCoinc = coincFileHead;
      }
      for ( ; thisCoinc->next; thisCoinc = thisCoinc->next );
    }

  }

  thisCoinc = NULL;

  /* sort triggers by statistic */
  if ( coincstat == effective_snrsq )
  {
    if ( vrbflg ) fprintf( stdout,
        "sorting zero-lag coinc inspiral trigger list by effective snr..." );
    coincZeroHead = XLALSortCoincInspiralByStat( coincZeroHead,
        *XLALCompareCoincInspiralByStat, &bittenLParams, &coincstat );
    if ( vrbflg ) fprintf( stdout, "done\n" );

    if ( vrbflg ) fprintf( stdout,
        "sorting slide coinc inspiral trigger list by effective snr..." );
    coincSlideHead = XLALSortCoincInspiralByStat( coincSlideHead,
        *XLALCompareCoincInspiralByStat, &bittenLParams, &coincstat );
    if ( vrbflg ) fprintf( stdout, "done\n" );
  }
  else
  {
    exit( 1 );
  }

  if ( fitNum )
  {
    if ( vrbflg ) fprintf( stdout,
        "Calculating background tail fit with %d loudest triggers...\n",
        fitNum );
    XLALCalcExpFitNLoudestBackground( coincSlideHead, fitNum, coincstat,
        &bittenLParams, &fitStat, &fitA, &fitB );
    if ( vrbflg ) fprintf( stdout,
        "Got a fit of %1.2e * exp(%1.2e * stat**2) above a stat**2 of %6.2f\n",
        fitA, fitB, fitStat );
  }

  {
    CoincInspiralSlideTable *slideHeads = NULL;
    CoincInspiralSlideTable *thisSlideHead = NULL;
    CoincInspiralTable      *thisEvent = NULL;
    CoincInspiralTable      *prevEvent = NULL;

    XLALCreateCoincSlideTable( &slideHeads, numSlides );

    timeAnalyzed = XLALSetupCoincSlideTable( slideHeads, coincSlideHead,
        timeAnalyzedFileName, timeModifier, numSlides );
    if( timeAnalyzed < 0 )
    {
      fprintf(stderr,
          "Unable to setup CoincSlideTable from slide coincs");
      exit( 1 );
    }

    thisSlideHead = slideHeads;

    /* calculating the FAR for the coincs */
    loudestRate = XLALRateErrorCalcCoincInspiral( coincZeroHead,
        thisSlideHead, coincstat, &bittenLParams, numSlides, timeAnalyzed,
        fitStat, fitA, fitB );
    if( loudestRate < 0 && coincZeroHead )
    {
      fprintf(stderr,
          "Error in calculating the FAR");
      exit( 1 );
    }
    /* for summary file, get calculate background time analyzed and get the
     * type of ifo coincidence (this is retrieved from the background in case
     * there are no foreground triggers */
    if ( summFileName )
    {
      thisSlideHead = slideHeads;
      while ( thisSlideHead )
      {
        bkgtimeAnalyzed = bkgtimeAnalyzed + thisSlideHead->slideTimeAnalyzed;
        /* get coinc ifo; this is probably an arcane way to do it */
        if ( !coincifos && thisSlideHead->coincInspiral )
        {
        /* allocate memory for coincifos string; needed size deteremined using
         * the numIfos element in the coincInspiral table */
          coincifos = (CHAR *) calloc(thisSlideHead->coincInspiral->numIfos * 2 + 1, sizeof(CHAR));
          for( j = 0; j < 6; ++j )
          {
            if ( thisSlideHead->coincInspiral->snglInspiral[j] )
            {
              strcat(coincifos, thisSlideHead->coincInspiral->snglInspiral[j]->ifo);
            }
          }
         }
         thisSlideHead = thisSlideHead->next;
       }
      }

    numEventsAboveThresh = XLALCountCoincInspiral( coincZeroHead );
    if ( vrbflg ) fprintf( stdout,
        "Loudest zero-lag coinc has a rate of %6.2f\n", loudestRate );

    if ( loudestFileName )
    {
      FILE *out;
      CoincInspiralTable   *tmpCoinc = NULL;

      out = fopen( loudestFileName, "w" );
      fprintf( out, "[corse]\n" );
      fprintf( out, "rate-threshold = %f\n", loudestRate );
      if ( coincZeroHead )
      {
        fprintf( out, "stat-threshold = %f\n",
            XLALCoincInspiralStat(coincZeroHead, coincstat, &bittenLParams) );

        /* free all but the loudest coinc inspirals */
        thisCoinc = coincZeroHead->next;
        while ( thisCoinc )
        {
          tmpCoinc = thisCoinc;
          thisCoinc= thisCoinc->next;
          XLALFreeCoincInspiral( &tmpCoinc );
        }
        coincZeroHead->next = NULL;
      }
      fclose( out );
    }

    /* free the CoincInspiralSlideTables from slideHeads */
    while ( slideHeads )
    {
      thisSlideHead = slideHeads;

      /* free all the coinc inspirals */
      coincSlideHead = thisSlideHead->coincInspiral;
      while ( coincSlideHead )
      {
        thisCoinc = coincSlideHead;
        coincSlideHead = thisCoinc->next;
        XLALFreeCoincInspiral( &thisCoinc );
      }

      slideHeads = slideHeads->next;
      LALFree( thisSlideHead );
    }
  }

  /* perform the statistic cut */
  if( rateThreshold >= 0 && statThreshold >= 0 )
  {
    coincZeroHead = XLALRateStatCutCoincInspiral ( coincZeroHead,
        coincstat, &bittenLParams, statThreshold, rateThreshold );
    numEventsAboveThresh = XLALCountCoincInspiral( coincZeroHead );
    if ( vrbflg ) fprintf( stdout,
        "Kept %d coincs below a rate threshold of %6.2f\n"
        "and above a statistic threshold of %6.2f\n", numEventsAboveThresh,
        rateThreshold, statThreshold );
  }

  if ( vrbflg )
  {
    fprintf( stdout, "Read in %d zero-lag triggers\n", numZeroTriggers );
    fprintf( stdout, "Recreated %d zero-lag coincs\n", numZeroCoincs );
    fprintf( stdout, "Read in %d slide triggers\n", numSlideTriggers );
    fprintf( stdout, "Recreated %d slide coincs\n", numSlideCoincs );
    if ( dataType != all_data ) 
    {
      fprintf( stdout, 
          "Have %d zero-lag coincs after play test\n", numEventsPlayTest);
    }
    if ( rateThreshold >= 0 )
    {
      fprintf( stdout,
          "Have %d zero-lag coincs below a rate of %6.2f\n"
          "and above a statistic threshold of %6.2f\n", numEventsAboveThresh,
          rateThreshold, statThreshold);
    }
    else
    {
      fprintf( stdout,
          "Have %d zero-lag coincs above a statistic threshold of %6.2f\n",
          numEventsAboveThresh, statThreshold);
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
  if ( sortTriggers || injectFileName )
  {
    if ( vrbflg ) fprintf( stdout,
        "sorting zero-lag coinc inspiral trigger list..." );
    coincZeroHead = XLALSortCoincInspiral( coincZeroHead, 
        *XLALCompareCoincInspiralByTime );
    if ( vrbflg ) fprintf( stdout, "done\n" );
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

    
    if ( vrbflg ) fprintf( stdout, 
        "Sorting single inspiral triggers before injection coinc test\n" );
    inspiralZeroEventList = XLALSortSnglInspiral( inspiralZeroEventList, 
        *LALCompareSnglInspiralByTime );
    
    /* first find singles which are coincident with injections */
    numSnglFound = XLALSnglSimInspiralTest( &simEventHead, 
        &inspiralZeroEventList, &missedSimHead, &missedSnglHead,
        injectWindowNS );

    if ( vrbflg ) fprintf( stdout, "%d injections found in single ifo\n", 
        numSnglFound );

    /* then check for coincs coincident with injections */
    numCoincFound = XLALCoincSimInspiralTest ( &simEventHead,  &coincZeroHead, 
        &missedSimCoincHead, &missedCoincHead );

    if ( vrbflg ) fprintf( stdout, "%d injections found in coincidence\n", 
        numCoincFound );

    if ( numCoincFound )
    {
      for ( thisCoinc = coincZeroHead; thisCoinc; thisCoinc = thisCoinc->next,
          numEventsCoinc++ );
      if ( vrbflg ) fprintf( stdout, "%d coincs found at times of injection\n",
          numEventsCoinc );
    }
    
    if ( missedSimCoincHead )
    {
      /* add these to the list of missed Sim's */
      if ( missedSimHead )
      {
        for (thisSimEvent = missedSimHead; thisSimEvent->next; 
            thisSimEvent = thisSimEvent->next );
        thisSimEvent->next = missedSimCoincHead;
      }
      else
      {
        missedSimHead = missedSimCoincHead;
      }
    }

    /* free the missed singles and coincs */
    while ( missedCoincHead )
    {
      thisCoinc = missedCoincHead;
      missedCoincHead = missedCoincHead->next;
      XLALFreeCoincInspiral( &thisCoinc );
    }

    while ( missedSnglHead )
    {
      thisSngl = missedSnglHead;
      missedSnglHead = missedSnglHead->next;
      XLALFreeSnglInspiral( &thisSngl );
    }

  } 


  /*
   *
   * write output data
   *
   */


  /* write out all coincs as singles with event IDs */
  snglOutput = XLALExtractSnglInspiralFromCoinc( coincZeroHead, 
      NULL, 0);


  /* write the main output file containing found injections */
  if ( vrbflg ) fprintf( stdout, "writing output xml files... " );
  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlStream, outputFileName ), 
      &status );

  /* write out the process and process params tables */
  if ( vrbflg ) fprintf( stdout, "process... " );
  LAL_CALL(LALGPSTimeNow(&status, &(proctable.processTable->end_time), 
        &accuracy), &status);
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, process_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, proctable, 
        process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );

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
  if ( snglOutput )
  {
    if ( vrbflg ) fprintf( stdout, "sngl_inspiral... " );
    outputTable.snglInspiralTable = snglOutput;
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, 
          sngl_inspiral_table ), &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable, 
          sngl_inspiral_table ), &status );
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

    /* write out the process and process params tables */
    if ( vrbflg ) fprintf( stdout, "process... " );
    LAL_CALL(LALGPSTimeNow(&status, &(proctable.processTable->end_time),
          &accuracy), &status);
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, process_table ),
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, proctable,
          process_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlStream ), &status );

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

    if ( missedSimHead )
    {
      outputTable.simInspiralTable = missedSimHead;
      LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlStream, 
            sim_inspiral_table ), &status );
      LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, outputTable, 
            sim_inspiral_table ), &status );
      LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream ), &status );
    }

    LAL_CALL( LALCloseLIGOLwXMLFile( &status, &xmlStream ), &status );
    if ( vrbflg ) fprintf( stdout, "done\n" );
  }
  
  /* free process proctable data after found and missed files have been written */  
  free( proctable.processTable );

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

    if ( !coincifos )
    {
      coincifos = (CHAR *) calloc(14, sizeof(CHAR));
      strcpy(coincifos, "no_background");
    }
    fprintf( fp, "coincident ifos: %s\n", coincifos );
    fprintf( fp, "read zero-lag triggers from %d files\n", numInZeroFiles );
    fprintf( fp, "number of zero-lag triggers in input files: %d \n",
        numZeroTriggers );
    fprintf( fp, "number of reconstructed zero-lag coincidences: %d \n",
        numZeroCoincs );

    if ( dataType != all_data )
    {
      fprintf( fp, "number of zero-lag triggers after playground test: %d \n", 
          numEventsPlayTest);
    }

    fprintf( fp, "read slide triggers from %d files\n", numInSlideFiles );
    fprintf( fp, "number of slide triggers in input files: %d \n",
        numSlideTriggers );
    fprintf( fp, "number of reconstructed slide coincidences: %d \n",
        numSlideCoincs );

    if ( rateThreshold >= 0 )
    {
      fprintf( fp, "number of zero-lag triggers with rate below %6.2f\n"
          "and statistic above %6.2f is: %d \n", rateThreshold, statThreshold,
          numEventsAboveThresh);
    }
    else
    {
      fprintf( fp, "number of zero-lag triggers with statistic above %6.2f"
          " is: %d \n", statThreshold, numEventsAboveThresh);
    }

    XLALINT8NSToGPS( &triggerTime, triggerInputTimeNS );
    fprintf( fp, "amount of time analysed for triggers %d sec %d ns\n", 
        triggerTime.gpsSeconds, triggerTime.gpsNanoSeconds );

    fprintf( fp, "amount of background time analyzed %f sec\n", bkgtimeAnalyzed );


    if ( massTag )
    {
      fprintf( fp, "mass-bin: %s\n", massTag );
    }

    if ( injectFileName )
    {
      fprintf( fp, "read %d injections from file %s\n", 
          numSimEvents, injectFileName );

      fprintf( fp, "number of injections in input data: %d\n", numSimInData );
      fprintf( fp, "number of injections found in input data: %d\n", 
          numCoincFound );
      fprintf( fp, 
          "number of triggers found within %lld msec of injection: %d\n",
          (injectWindowNS / 1000000LL), numEventsCoinc );

      fprintf( fp, "efficiency: %f \n", 
          (REAL4) numCoincFound / (REAL4) numSimInData );
    }

    fclose( fp ); 
  }


  /*
   *
   * free memory and exit
   *
   */


  /* free the coinc inspirals */
  while ( coincZeroHead )
  {
    thisCoinc = coincZeroHead;
    coincZeroHead = thisCoinc->next;
    XLALFreeCoincInspiral( &thisCoinc );
  }

  while ( coincSlideHead )
  {   
    thisCoinc = coincSlideHead;
    coincSlideHead = thisCoinc->next;
    XLALFreeCoincInspiral( &thisCoinc );
  }

  /* free the inspiral events we saved */
  while ( inspiralZeroEventList )
  {
    thisSngl = inspiralZeroEventList;
    inspiralZeroEventList = thisSngl->next;
    XLALFreeSnglInspiral( &thisSngl );
  }

  while ( inspiralSlideEventList )
  {     
    thisSngl = inspiralSlideEventList;
    inspiralSlideEventList = thisSngl->next;
    XLALFreeSnglInspiral( &thisSngl );
  }    

  while ( snglOutput )
  {
    thisInspiralTrigger = snglOutput;
    snglOutput = snglOutput->next;
    XLALFreeSnglInspiral( &thisInspiralTrigger );
  }

  /* free the process params */
  while( procparams.processParamsTable )
  {
    this_proc_param = procparams.processParamsTable;
    procparams.processParamsTable = this_proc_param->next;
    free( this_proc_param );
  }

  while ( summValueList )
  {
    SummValueTable *thisSummValue;
    thisSummValue = summValueList;
    summValueList = summValueList->next;
    LALFree( thisSummValue );
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

  /* free the zero-lag input file name data */
  if ( inputGlobZero )
  {
    LALFree( inZeroFileNameList ); 
    globfree( &globbedZeroFiles );
  }
  else
  {
    for ( j = 0; j < numInZeroFiles; ++j )
    {
      LALFree( inZeroFileNameList[j] );
    }
    LALFree( inZeroFileNameList );
  }

  /* free the slide input file name data */
  if ( inputGlobSlide )
  {
    LALFree( inSlideFileNameList );
    globfree( &globbedSlideFiles );
  }
  else
  {
    for ( j = 0; j < numInSlideFiles; ++j )
    {
      LALFree( inSlideFileNameList[j] );
    }
    LALFree( inSlideFileNameList );
  }

  /* free zero-lag input files list */
  while ( inputZeroFiles )
  {
    thisInputFile = inputZeroFiles;
    inputZeroFiles = thisInputFile->next;
    LALFree( thisInputFile );
  }

  /* free slide input files list */
  while ( inputSlideFiles )
  {
    thisInputFile = inputSlideFiles;
    inputSlideFiles = thisInputFile->next;
    LALFree( thisInputFile );
  }

  /* free search summaries read in */
  while ( searchSummList )
  {
    thisSearchSumm = searchSummList;
    searchSummList = thisSearchSumm->next;
    LALFree( thisSearchSumm );
  }

  /* free search summaries read in */
  while ( searchSummSlideList )
  {
    thisSearchSumm = searchSummSlideList;
    searchSummSlideList = thisSearchSumm->next;
    LALFree( thisSearchSumm );
  }

  if ( vrbflg ) fprintf( stdout, "checking memory leaks and exiting\n" );
  LALCheckMemoryLeaks();
  exit( 0 );
}
