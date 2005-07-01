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
#include <lal/TimeDelay.h>
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

int haveTrig[LAL_NUM_IFO];
int checkTimes = 0;
int multiIfoCoinc = 0;


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
      "  [--ifo-tag]       ifotag      set the ifo-tag - for file naming\n"\
      "  [--comment]       string      set the process table comment to STRING\n"\
      "\n"\
      "   --gps-start-time start_time  GPS second of data start time\n"\
      "   --gps-end-time   end_time    GPS second of data end time\n"\
      "  [--check-times]               Check that all times were analyzed\n"\
      "  [--multi-ifo-coinc]           Look for triple/quadruple ifo coincidence\n"\
      "  [--maximization-interval] max_dt set length of maximization interval in ms\n"\
      "\n"\
      "  [--g1-slide]      g1_slide    Slide G1 data by multiples of g1_slide\n"\
      "  [--h1-slide]      h1_slide    Slide H1 data by multiples of h1_slide\n"\
      "  [--h2-slide]      h2_slide    Slide H2 data by multiples of h2_slide\n"\
      "  [--l1-slide]      l1_slide    Slide L1 data by multiples of l1_slide\n"\
      "  [--t1-slide]      t1_slide    Slide T1 data by multiples of t1_slide\n"\
      "  [--v1-slide]      v1_slide    Slide V1 data by multiples of v1_slide\n"\
      "  [--num-slides]    num_slides  The number of time slides to perform\n"\
      "\n"\
      "  [--g1-triggers]               input triggers from G1\n"\
      "  [--h1-triggers]               input triggers from H1\n"\
      "  [--h2-triggers]               input triggers from H2\n"\
      "  [--l1-triggers]               input triggers from L1\n"\
      "  [--t1-triggers]               input triggers from T1\n"\
      "  [--v1-triggers]               input triggers from V1\n"\
      "\n"\
      "   --parameter-test     test    set parameters with which to test coincidence:\n"\
      "                                (m1_and_m2|mchirp_and_eta|psi0_and_psi3)\n"\
      "  [--g1-time-accuracy]  g1_dt   specify the timing accuracy of G1 in ms\n"\
      "  [--h1-time-accuracy]  h1_dt   specify the timing accuracy of H1 in ms\n"\
      "  [--h2-time-accuracy]  h2_dt   specify the timing accuracy of H2 in ms\n"\
      "  [--l1-time-accuracy]  l1_dt   specify the timing accuracy of L1 in ms\n"\
      "  [--t1-time-accuracy]  t1_dt   specify the timing accuracy of T1 in ms\n"\
      "  [--v1-time-accuracy]  v1_dt   specify the timing accuracy of V1 in ms\n"\
      "\n"\
      "  [--g1-mass-accuracy]  g1_dm   specify the mass accuracy of G1\n"\
      "  [--h1-mass-accuracy]  h1_dm   specify the mass accuracy of H1\n"\
      "  [--h2-mass-accuracy]  h2_dm   specify the mass accuracy of H2\n"\
      "  [--l1-mass-accuracy]  l1_dm   specify the mass accuracy of L1\n"\
      "  [--t1-mass-accuracy]  t1_dm   specify the mass accuracy of T1\n"\
      "  [--v1-mass-accuracy]  v1_dm   specify the mass accuracy of V1\n"\
      "\n"\
      "  [--g1-mchirp-accuracy] g1_dmchirp  specify the mchirp accuracy of G1\n"\
      "  [--h1-mchirp-accuracy] h1_dmchirp  specify the mchirp accuracy of H1\n"\
      "  [--h2-mchirp-accuracy] h2_dmchirp  specify the mchirp accuracy of H2\n"\
      "  [--l1-mchirp-accuracy] l1_dmchirp  specify the mchirp accuracy of L1\n"\
      "  [--t1-mchirp-accuracy] t1_dmchirp  specify the mchirp accuracy of T1\n"\
      "  [--v1-mchirp-accuracy] v1_dmchirp  specify the mchirp accuracy of V1\n"\
      "\n"\
      "  [--g1-eta-accuracy] g1_deta   specify the eta accuracy of G1\n"\
      "  [--h1-eta-accuracy] h1_deta   specify the eta accuracy of H1\n"\
      "  [--h2-eta-accuracy] h2_deta   specify the eta accuracy of H2\n"\
      "  [--l1-eta-accuracy] l1_deta   specify the eta accuracy of L1\n"\
      "  [--t1-eta-accuracy] t1_deta   specify the eta accuracy of T1\n"\
      "  [--v1-eta-accuracy] v1_deta   specify the eta accuracy of V1\n"\
      "\n"\
      "  [--g1-psi0-accuracy]  g1_dpsi0   specify the psi0 accuracy of G1\n"\
      "  [--h1-psi0-accuracy]  h1_dpsi0   specify the psi0 accuracy of H1\n"\
      "  [--h2-psi0-accuracy]  h2_dpsi0   specify the psi0 accuracy of H2\n"\
      "  [--l1-psi0-accuracy]  l1_dpsi0   specify the psi0 accuracy of L1\n"\
      "  [--t1-psi0-accuracy]  t1_dpsi0   specify the psi0 accuracy of T1\n"\
      "  [--v1-psi0-accuracy]  v1_dpsi0   specify the psi0 accuracy of V1\n"\
      "\n"\
      "  [--g1-psi3-accuracy]  g1_dpsi3   specify the psi3 accuracy of G1\n"\
      "  [--h1-psi3-accuracy]  h1_dpsi3   specify the psi3 accuracy of H1\n"\
      "  [--h2-psi3-accuracy]  h2_dpsi3   specify the psi3 accuracy of H2\n"\
      "  [--l1-psi3-accuracy]  l1_dpsi3   specify the psi3 accuracy of L1\n"\
      "  [--t1-psi3-accuracy]  t1_dpsi3   specify the psi3 accuracy of T1\n"\
      "  [--v1-psi3-accuracy]  v1_dpsi3   specify the psi3 accuracy of V1\n"\
      "\n"\
      "   --data-type        data_type specify the data type, must be one of\n"\
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

  LALPlaygroundDataMask dataType = unspecified_data_type;
  INT4  startCoincidence = -1;
  LIGOTimeGPS startCoinc = {0,0};
  INT4  endCoincidence = -1;
  LIGOTimeGPS endCoinc = {0,0};

  INT4         slideStep[LAL_NUM_IFO];
  LIGOTimeGPS  slideTimes[LAL_NUM_IFO];
  LIGOTimeGPS  slideReset[LAL_NUM_IFO];
  INT4         numSlides = 0;
  INT4         slideNum  = 0;

  CHAR  ifoName[MAXIFO][LIGOMETA_IFO_MAX];
  CHAR  ifos[LIGOMETA_IFOS_MAX];
  CHAR  ifoA[LIGOMETA_IFO_MAX];
  CHAR  ifoB[LIGOMETA_IFO_MAX];
  
  CHAR  comment[LIGOMETA_COMMENT_MAX];
  CHAR *userTag = NULL;
  CHAR *ifoTag = NULL;

  CHAR  fileName[FILENAME_MAX];
  CHAR  fileSlide[FILENAME_MAX];

  UINT4  numIFO = 0;
  UINT4  numTrigIFO = 0;
  UINT4  numTriggers = 0;
  UINT4  numCoinc = 0;
  UINT4  numDoubles = 0;
  UINT4  numTriples = 0;
  UINT4  numQuadruples = 0;
  UINT4  numTrigs[LAL_NUM_IFO];
  UINT4  N = 0;

  LALDetector          aDet;
  LALDetector          bDet;

  SnglInspiralTable    *inspiralEventList = NULL;
  SnglInspiralTable    *thisInspiralTrigger = NULL;
  SnglInspiralTable    *snglOutput = NULL;
  CoincInspiralTable   *coincInspiralList = NULL;
  CoincInspiralTable   *thisCoinc = NULL;

  InspiralAccuracyList  accuracyParams;

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

  InterferometerNumber  ifoNumber = LAL_UNKNOWN_IFO;
  InterferometerNumber  ifoTwo    = LAL_UNKNOWN_IFO;
  INT4                  i;
  INT4                  maximizationInterval = 0;

  const CHAR                   ifoList[LAL_NUM_IFO][LIGOMETA_IFO_MAX] = 
                                   {"G1", "H1", "H2", "L1", "T1", "V1"};
  const CHAR                  *ifoArg[LAL_NUM_IFO] = 
                                   {"g1-triggers", "h1-triggers", 
                                    "h2-triggers", "l1-triggers", 
                                    "t1-triggers", "v1-triggers"};


  /* getopt arguments */
  struct option long_options[] =
  {
    {"verbose",             no_argument,   &vrbflg,                   1 },
    {"g1-triggers",         no_argument,   &(haveTrig[LAL_IFO_G1]),   1 },
    {"h1-triggers",         no_argument,   &(haveTrig[LAL_IFO_H1]),   1 },
    {"h2-triggers",         no_argument,   &(haveTrig[LAL_IFO_H2]),   1 },
    {"l1-triggers",         no_argument,   &(haveTrig[LAL_IFO_L1]),   1 },
    {"t1-triggers",         no_argument,   &(haveTrig[LAL_IFO_T1]),   1 },
    {"v1-triggers",         no_argument,   &(haveTrig[LAL_IFO_V1]),   1 },
    {"check-times",         no_argument,   &checkTimes,               1 },
    {"multi-ifo-coinc",     no_argument,   &multiIfoCoinc,            1 },
    {"g1-slide",            required_argument, 0,                    'b'},
    {"h1-slide",            required_argument, 0,                    'c'},
    {"h2-slide",            required_argument, 0,                    'd'},
    {"l1-slide",            required_argument, 0,                    'e'},
    {"t1-slide",            required_argument, 0,                    'f'},
    {"v1-slide",            required_argument, 0,                    'g'},
    {"num-slides",          required_argument, 0,                    'T'},
    {"g1-time-accuracy",    required_argument, 0,                    'A'},
    {"h1-time-accuracy",    required_argument, 0,                    'B'}, 
    {"h2-time-accuracy",    required_argument, 0,                    'C'}, 
    {"l1-time-accuracy",    required_argument, 0,                    'D'},
    {"t1-time-accuracy",    required_argument, 0,                    'E'}, 
    {"v1-time-accuracy",    required_argument, 0,                    'F'}, 
    {"g1-mass-accuracy",    required_argument, 0,                    'G'},
    {"h1-mass-accuracy",    required_argument, 0,                    'H'},
    {"h2-mass-accuracy",    required_argument, 0,                    'I'},
    {"l1-mass-accuracy",    required_argument, 0,                    'J'},
    {"t1-mass-accuracy",    required_argument, 0,                    'K'},
    {"v1-mass-accuracy",    required_argument, 0,                    'L'},
    {"g1-mchirp-accuracy",  required_argument, 0,                    'M'},
    {"h1-mchirp-accuracy",  required_argument, 0,                    'N'},
    {"h2-mchirp-accuracy",  required_argument, 0,                    'O'},
    {"l1-mchirp-accuracy",  required_argument, 0,                    'P'},
    {"t1-mchirp-accuracy",  required_argument, 0,                    'Q'},
    {"v1-mchirp-accuracy",  required_argument, 0,                    'R'},
    {"g1-eta-accuracy",     required_argument, 0,                    'm'},
    {"h1-eta-accuracy",     required_argument, 0,                    'n'},
    {"h2-eta-accuracy",     required_argument, 0,                    'o'},
    {"l1-eta-accuracy",     required_argument, 0,                    'p'},
    {"t1-eta-accuracy",     required_argument, 0,                    'q'},
    {"v1-eta-accuracy",     required_argument, 0,                    'r'},
    {"g1-psi0-accuracy",    required_argument, 0,                    '2'},
    {"h1-psi0-accuracy",    required_argument, 0,                    '3'},
    {"h2-psi0-accuracy",    required_argument, 0,                    '4'},
    {"l1-psi0-accuracy",    required_argument, 0,                    '5'},
    {"t1-psi0-accuracy",    required_argument, 0,                    '6'},
    {"v1-psi0-accuracy",    required_argument, 0,                    '7'},
    {"g1-psi3-accuracy",    required_argument, 0,                    '8'},
    {"h1-psi3-accuracy",    required_argument, 0,                    '9'},
    {"h2-psi3-accuracy",    required_argument, 0,                    '!'},
    {"l1-psi3-accuracy",    required_argument, 0,                    '-'},
    {"t1-psi3-accuracy",    required_argument, 0,                    '+'},
    {"v1-psi3-accuracy",    required_argument, 0,                    '='},
    {"parameter-test",      required_argument, 0,                    'a'},
    {"gps-start-time",      required_argument, 0,                    's'},
    {"gps-end-time",        required_argument, 0,                    't'},
    {"maximization-interval",required_argument, 0,                   '@'},    
    {"data-type",           required_argument, 0,                    'k'},
    {"comment",             required_argument, 0,                    'x'},
    {"user-tag",            required_argument, 0,                    'Z'},
    {"userTag",             required_argument, 0,                    'Z'},
    {"ifo-tag",             required_argument, 0,                    'i'},
    {"help",                no_argument,       0,                    'h'}, 
    {"debug-level",         required_argument, 0,                    'z'},
    {"version",             no_argument,       0,                    'V'},
    {0, 0, 0, 0}
  };
  int c;


  /*
   * 
   * initialize things
   *
   */

  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "33" );
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

  memset( &accuracyParams, 0, sizeof(InspiralAccuracyList) );
  memset( &aDet, 0, sizeof(LALDetector) );

  /* set the time slide data to zero */
  memset( &slideStep, 0, LAL_NUM_IFO * sizeof(INT4) );
  memset( &slideTimes, 0, LAL_NUM_IFO * sizeof(LIGOTimeGPS) );
  memset( &slideReset, 0, LAL_NUM_IFO * sizeof(LIGOTimeGPS) );
  memset( &haveTrig, 0, LAL_NUM_IFO * sizeof(int) );
  
  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    long int gpstime;
    size_t optarg_len;

    c = getopt_long_only( argc, argv, 
        "A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:P:Q:R:S:T:U:VZ:"
        "a:b:c:d:e:f:g:hi:k:m:n:o:p:q:r:s:t:x:z:"
        "2:3:4:5:6:7:8:9:!:-:+:=:@:", 
        long_options, &option_index );

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
        /* set the parameter test */
        if ( ! strcmp( "m1_and_m2", optarg ) )
        {
          accuracyParams.test = m1_and_m2;
        }
        else if ( ! strcmp( "psi0_and_psi3", optarg ) )
        {
          accuracyParams.test = psi0_and_psi3;
        }
        else if ( ! strcmp( "mchirp_and_eta", optarg ) )
        {
          accuracyParams.test = mchirp_and_eta;
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

      case 'A':
        /* time accuracy G1, argument is in milliseconds */
        accuracyParams.ifoAccuracy[LAL_IFO_G1].dt = atof(optarg) * 1000000LL;
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
      
      case 'B':
        /* time accuracy H1, argument is in milliseconds */
        accuracyParams.ifoAccuracy[LAL_IFO_H1].dt = atof(optarg) * 1000000LL;
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'C':
        /* time accuracy H2, argument is in milliseconds */
        accuracyParams.ifoAccuracy[LAL_IFO_H2].dt = atof(optarg) * 1000000LL;
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'D':
        /* time accuracy L1, argument is in milliseconds */
        accuracyParams.ifoAccuracy[LAL_IFO_L1].dt = atof(optarg) * 1000000LL;
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'E':
        /* time accuracy T1, argument is in milliseconds */
        accuracyParams.ifoAccuracy[LAL_IFO_T1].dt = atof(optarg) * 1000000LL;
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'F':
        /* time accuracy V1, argument is in milliseconds */
        accuracyParams.ifoAccuracy[LAL_IFO_V1].dt = atof(optarg) * 1000000LL;
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'G':
        /* mass accuracy G1, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_G1].dm = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
      
      case 'H':
        /* mass accuracy H1, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_H1].dm = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'I':
        /* mass accuracy H2, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_H2].dm = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'J':
        /* mass accuracy L1, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_L1].dm = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'K':
        /* mass accuracy T1, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_T1].dm = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'L':
        /* mass accuracy V1, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_V1].dm = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'M':
        /* chirp mass accuracy G1, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_G1].dmchirp = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
      
      case 'N':
        /* chirp mass accuracy H1, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_H1].dmchirp = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'O':
        /* chirp mass accuracy H2, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_H2].dmchirp = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'P':
        /* chirp mass accuracy L1, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_L1].dmchirp = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'Q':
        /* chirp mass accuracy T1, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_T1].dmchirp = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'R':
        /* chirp mass accuracy V1, argument is in solar masses */
        accuracyParams.ifoAccuracy[LAL_IFO_V1].dmchirp = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
       
      case 'b':
        /* slide time for G1 */
        slideStep[LAL_IFO_G1] = atoi( optarg );
        if ( slideStep[LAL_IFO_G1] < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "The slideStep must be positive\n"
              "(%d specified)\n",
              long_options[option_index].name, slideStep[LAL_IFO_G1] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%ld", slideStep[LAL_IFO_G1] );
        break;

      case 'c':
        /* slide time for H1 */
        slideStep[LAL_IFO_H1] = atoi( optarg );
        if ( slideStep[LAL_IFO_H1] < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "The slideStep must be positive\n"
              "(%d specified)\n",
              long_options[option_index].name, slideStep[LAL_IFO_H1] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%ld", slideStep[LAL_IFO_H1] );
        break;
        
      case 'd':
        /* slide time for H2 */
        slideStep[LAL_IFO_H2] = atoi( optarg );
        if ( slideStep[LAL_IFO_H2] < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "The slideStep must be positive\n"
              "(%d specified)\n",
              long_options[option_index].name, slideStep[LAL_IFO_H2] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%ld", slideStep[LAL_IFO_H2] );
        break;
        
      case 'e':
        /* slide time for L1 */
        slideStep[LAL_IFO_L1] = atoi( optarg );
        if ( slideStep[LAL_IFO_L1] < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "The slideStep must be positive\n"
              "(%d specified)\n",
              long_options[option_index].name, slideStep[LAL_IFO_L1] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%ld", slideStep[LAL_IFO_L1] );
        break;
        
      case 'f':
        /* slide time for T1 */
        slideStep[LAL_IFO_T1] = atoi( optarg );
        if ( slideStep[LAL_IFO_T1] < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "The slideStep must be positive\n"
              "(%d specified)\n",
              long_options[option_index].name, slideStep[LAL_IFO_T1] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%ld", slideStep[LAL_IFO_T1] );
        break;
        
      case 'g':
        /* slide time for V1 */
        slideStep[LAL_IFO_V1] = atoi( optarg );
        if ( slideStep[LAL_IFO_V1] < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "The slideStep must be positive\n"
              "(%d specified)\n",
              long_options[option_index].name, slideStep[LAL_IFO_V1] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%ld", slideStep[LAL_IFO_V1] );
        break;

      case 'T':
        /* num slides*/
        numSlides = atoi( optarg );
        if ( numSlides < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "The number of time slides must be positive\n"
              "(%d specified)\n",
              long_options[option_index].name, numSlides );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%ld", numSlides );
        break;

      case 'm':
        /* eta accuracy G1, argument is dimensionless */
        accuracyParams.ifoAccuracy[LAL_IFO_G1].deta = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
      
      case 'n':
        /* eta accuracy H1, argument is dimensionless */
        accuracyParams.ifoAccuracy[LAL_IFO_H1].deta = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'o':
        /* eta accuracy H2, argument is dimensionless */
        accuracyParams.ifoAccuracy[LAL_IFO_H2].deta = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'p':
        /* eta accuracy L1, argument is dimensionless */
        accuracyParams.ifoAccuracy[LAL_IFO_L1].deta = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'q':
        /* eta accuracy T1, argument is dimensionless */
        accuracyParams.ifoAccuracy[LAL_IFO_T1].deta = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'r':
        /* eta accuracy V1, argument is dimensionless */
        accuracyParams.ifoAccuracy[LAL_IFO_V1].deta = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 's':
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

      case 't':
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

      case 'x':
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
        
      case 'i':
        /* create storage for the ifotag */
        optarg_len = strlen(optarg) + 1;
        ifoTag = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( ifoTag, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;


      case 'V':
        /* print version information and exit */
        fprintf( stdout, "The Hierarchical INspiral Coincidence Analysis\n" 
            "Steve Fairhurst\n"
            "CVS Version: " CVS_ID_STRING "\n"
            "CVS Tag: " CVS_NAME_STRING "\n" );
        exit( 0 );
        break;

      case '?':
        print_usage(argv[0]);
        exit( 1 );
        break;

      case '2':
     /* psi0 mass accuracy G1  */
        accuracyParams.ifoAccuracy[LAL_IFO_G1].dpsi0 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '3':
     /* psi0 mass accuracy H1  */
        accuracyParams.ifoAccuracy[LAL_IFO_H1].dpsi0 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '4':
     /* psi0 mass accuracy H2  */
        accuracyParams.ifoAccuracy[LAL_IFO_H2].dpsi0 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
     
      case '5':
     /* psi0 mass accuracy L1  */
        accuracyParams.ifoAccuracy[LAL_IFO_L1].dpsi0 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '6':
     /* psi0 mass accuracy T1  */
        accuracyParams.ifoAccuracy[LAL_IFO_T1].dpsi0 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '7':
     /* psi0 mass accuracy V1  */
        accuracyParams.ifoAccuracy[LAL_IFO_V1].dpsi0 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '8':
     /* psi3 mass accuracy G1  */
        accuracyParams.ifoAccuracy[LAL_IFO_G1].dpsi3 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '9':
     /* psi3 mass accuracy H1  */
        accuracyParams.ifoAccuracy[LAL_IFO_H1].dpsi3 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '!':
     /* psi3 mass accuracy H2  */
        accuracyParams.ifoAccuracy[LAL_IFO_H2].dpsi3 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '-':
     /* psi3 mass accuracy L1  */
        accuracyParams.ifoAccuracy[LAL_IFO_L1].dpsi3 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '+':
     /* psi3 mass accuracy T1  */
        accuracyParams.ifoAccuracy[LAL_IFO_T1].dpsi3 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '=':
     /* psi3 mass accuracy V1  */
        accuracyParams.ifoAccuracy[LAL_IFO_V1].dpsi3 = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case '@':
        /* set the maximization window */
        maximizationInterval = atoi( optarg );
        if ( maximizationInterval < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "maximization interval must be positive:\n "
              "(%d ms specified)\n",
              long_options[option_index].name, maximizationInterval );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%ld",  maximizationInterval );
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
  if ( accuracyParams.test == no_test )
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


  /* Store the IFOs we expect triggers from */
  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    if ( haveTrig[ifoNumber] )
    {
      /* write ifo name in ifoName list */
      LALSnprintf( ifoName[numIFO], LIGOMETA_IFO_MAX, ifoList[ifoNumber] );
      numIFO++;

      /* store the argument in the process_params table */
      this_proc_param = this_proc_param->next = (ProcessParamsTable *)
        calloc( 1, sizeof(ProcessParamsTable) );
      LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
          "%s", PROGRAM_NAME );
      LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", 
          ifoArg[ifoNumber]);
      LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
      LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );

      /* check that a non-zero timing accuracy was specified */
      if ( ! accuracyParams.ifoAccuracy[ifoNumber].dt )
      {
        fprintf( stderr, "Error: --dt must be specified for %s\n", 
            ifoName[ifoNumber]);
        exit(1);
      }
    }
  }

  
  /* check that we have at least two IFOs specified, or can't do coincidence */
  if ( numIFO < 2 )
  {
    fprintf( stderr, "Must specify at least two IFOs to do coincidence\n"
        "%d specified\n", numIFO );
    exit ( 1 );
  }

  if ( numIFO > 2 && vrbflg)
  {
    if ( !multiIfoCoinc )
    {
      fprintf( stdout, 
          "Finding all double coincidences in %d IFO time.\n"
          "If you want triples/quadruples please specify --multi-ifo-coinc.\n",
          numIFO);
    }
    else
    {
      fprintf( stdout, 
          "Finding all double/triple/quadruple coincidences in %d IFO time.\n"
          numIFO);

  
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

  /* if numSlides is set, check that the slide times are different for
   * different ifos (for which we have triggers */
  if( numSlides )
  {
    for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
    {
      XLALReturnIFO( ifoA, ifoNumber );
     
      if( vrbflg && haveTrig[ifoNumber] ) fprintf( stdout, 
          "Performing a slide of multiples of %d seconds on %s\n", 
          slideStep[ifoNumber], ifoA);

      for( ifoTwo = ifoNumber + 1; ifoTwo < LAL_NUM_IFO; ifoTwo++ )
      {
        if ( haveTrig[ifoTwo] && haveTrig[ifoNumber] &&
            slideStep[ifoTwo] == slideStep[ifoNumber] )
        {
          XLALReturnIFO( ifoB, ifoTwo );
           
          fprintf( stderr,
            "The time slide specified for ifo %s is %d\n"
            "The time slide specified for ifo %s is also %d\n"
            "Must specify unique time slides for all instruments\n",
            ifoA, slideStep[ifoNumber], ifoB, slideStep[ifoTwo]);

          exit( 1 );
        }
      }
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


  /* store the check-times in the process_params table */
  if ( checkTimes )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--check-times" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }

  /* store the multi-ifo-coinc in the process_params table */
  if ( multiIfoCoinc )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--multi-ifo-coinc" );
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

        /* maximize over a given interval */
        if ( maximizationInterval )
        {
          if (vrbflg)
          {
            fprintf( stdout, "Clustering triggers for over %d ms window\n",
                maximizationInterval);
          }
          XLALMaxSnglInspiralOverIntervals( &inputData, 
              (1.0e6 * maximizationInterval) );
        }

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
    for ( ifoNumber = 0; ifoNumber < numIFO; ++ifoNumber )
    {
      LAL_CALL( LALCheckOutTimeFromSearchSummary ( &status, searchSummList, 
            ifoName[ifoNumber], &startCoinc, &endCoinc ), &status);
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
  for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO ; ++ifoNumber )
  {
    LALIfoCountSingleInspiral(&status, &numTrigs[ifoNumber], 
        inspiralEventList, ifoNumber);
    if ( vrbflg ) fprintf( stdout, 
        "Have %d triggers from %s.\n", numTrigs[ifoNumber], 
        ifoList[ifoNumber] );
    if ( numTrigs[ifoNumber] && !haveTrig[ifoNumber] )
    {
      fprintf( stderr, "Read in triggers from %s, none expected.\n",
          ifoList[ifoNumber]);
      exit( 1 );
    }
    if ( haveTrig[ifoNumber] && numTrigs[ifoNumber] )
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

  /* Populate the lightTravel matrix */
  for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
  {
    XLALReturnDetector( &aDet, ifoNumber );
    
    for ( ifoTwo = 0; ifoTwo < LAL_NUM_IFO; ifoTwo++)
    {
      XLALReturnDetector( &bDet, ifoTwo );
      accuracyParams.lightTravelTime[ ifoNumber][ ifoTwo ] = 
        XLALLightTravelTime( &aDet, &bDet );
    }
  }

  /* 
   *  
   * check for two IFO coincidence
   *
   */

  if ( !numSlides )
  {
    LAL_CALL( LALCreateTwoIFOCoincList( &status, &coincInspiralList,
        inspiralEventList, &accuracyParams ), &status );
  
    if ( multiIfoCoinc )
    {
      for( N = 3; N <= numIFO; N++)
      {
        LAL_CALL( LALCreateNIFOCoincList( &status, &coincInspiralList, 
              &accuracyParams, N ), &status );
      }

      LAL_CALL( LALRemoveRepeatedCoincs( &status, &coincInspiralList ), 
          &status );
    }

    /* count the coincs */
    if( coincInspiralList )
    {
      for (numCoinc = 0, thisCoinc = coincInspiralList;
            thisCoinc; ++numCoinc, thisCoinc = thisCoinc->next )
      {
        if ( thisCoinc->numIfos == 2 )
        {
          ++numDoubles;
        }
        else if ( thisCoinc->numIfos == 3 )
        {
          ++numTriples;
        }
        else if ( thisCoinc->numIfos == 4 )
        {
          ++numQuadruples;
        }
      }     
    }
   
    if ( vrbflg ) 
    {
      fprintf( stdout, "%d coincident triggers found.\n", numCoinc );
      fprintf( stdout, "%d double coincident triggers\n"
                       "%d triple coincident triggers\n"
                       "%d quadruple coincident triggers\n",
                       numDoubles, numTriples, numQuadruples );
    }
    
    /* write out all coincs as singles with event IDs */
    LAL_CALL( LALExtractSnglInspiralFromCoinc( &status, &snglOutput, 
          coincInspiralList, &startCoinc, slideNum), &status );
  }
  else
  {
    /* perform the time slides */
    for( slideNum = -numSlides; slideNum <= numSlides; slideNum++ )
    {
      SnglInspiralTable    *slideOutput = NULL;
      INT4                  numCoincInSlide = 0;

      coincInspiralList = NULL;
      
      for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
      {
        if( slideNum == -numSlides )
        {
          INT4 tmpSlide = (- numSlides * slideStep[ifoNumber]);
          slideTimes[ifoNumber].gpsSeconds = tmpSlide;
          slideReset[ifoNumber].gpsSeconds = (-tmpSlide);
        }
        else
        {
          slideTimes[ifoNumber].gpsSeconds = slideStep[ifoNumber];
          slideReset[ifoNumber].gpsSeconds -= slideStep[ifoNumber];
        }
      }
    
      if ( vrbflg ) fprintf(stdout,
          "Performing time slide %d\n", slideNum );

      /* slide the data */
      LAL_CALL( LALTimeSlideSingleInspiral( &status, inspiralEventList,
            &startCoinc, &endCoinc, slideTimes), &status) ;
      LAL_CALL( LALSortSnglInspiral( &status, &(inspiralEventList),
            LALCompareSnglInspiralByTime ), &status );
      
      /* look for coincidences, except in the zero lag */
      if( slideNum )
      {
        LAL_CALL( LALCreateTwoIFOCoincList(&status, &coincInspiralList,
          inspiralEventList, &accuracyParams ), &status);

        if ( multiIfoCoinc )
        {
          for( N = 3; N <= numIFO; N++)
          {
            LAL_CALL( LALCreateNIFOCoincList( &status, &coincInspiralList, 
                  &accuracyParams, N ), &status );
          }

          LAL_CALL( LALRemoveRepeatedCoincs( &status, &coincInspiralList ), 
              &status );
        }

        /* count the coincs, scroll to end of list */
        if( coincInspiralList )
        {  
          for (numCoincInSlide = 1, thisCoinc = coincInspiralList; 
              thisCoinc->next; ++numCoincInSlide, thisCoinc = thisCoinc->next );

          if ( vrbflg ) fprintf( stdout,
              "%d coincident triggers found in slide.\n", numCoincInSlide );

          numCoinc += numCoincInSlide;


          /* write out all coincs as singles with event IDs */
          LAL_CALL( LALExtractSnglInspiralFromCoinc( &status, &slideOutput, 
                coincInspiralList, &startCoinc, slideNum), &status );

          /* the output triggers should be slid back to original time */
          LAL_CALL( LALTimeSlideSingleInspiral( &status, slideOutput,
                &startCoinc, &endCoinc, slideReset), &status) ;
          LAL_CALL( LALSortSnglInspiral( &status, &(slideOutput),
                LALCompareSnglInspiralByTime ), &status );

          while ( coincInspiralList )
          {
            thisCoinc = coincInspiralList;
            coincInspiralList = coincInspiralList->next;
            LALFree( thisCoinc );
          }

        }
      }
      
      if ( snglOutput )
      {
        thisInspiralTrigger->next = slideOutput;
      }
      else
      {
        snglOutput = slideOutput;
      }

      /* scroll to the end of the list */
      if ( slideOutput )
      {
        for( thisInspiralTrigger = slideOutput; thisInspiralTrigger->next; 
          thisInspiralTrigger = thisInspiralTrigger->next);
      }
    } 
  }


  /*
   *
   * write the output xml file
   *
   */


  /* since we don't yet write coinc inspiral tables, we must make a list of
   * sngl_inspiral tables with the eventId's appropriately poplulated */
  
   
cleanexit:

  searchsumm.searchSummaryTable->in_start_time = startCoinc;
  searchsumm.searchSummaryTable->in_end_time = endCoinc;
  searchsumm.searchSummaryTable->out_start_time = startCoinc;
  searchsumm.searchSummaryTable->out_end_time = endCoinc;
  searchsumm.searchSummaryTable->nnodes = 1;

  if ( vrbflg ) fprintf( stdout, "writing output file... " );

  if ( userTag && ifoTag)
  {
    LALSnprintf( fileName, FILENAME_MAX, "%s-THINCA_%s_%s-%d-%d.xml", 
        ifos, ifoTag, userTag, startCoincidence, 
        endCoincidence - startCoincidence );
    LALSnprintf( fileSlide, FILENAME_MAX, "%s-THINCA_SLIDE_%s_%s-%d-%d.xml", 
        ifos, ifoTag, userTag, startCoincidence, 
        endCoincidence - startCoincidence );
  }
  else if ( ifoTag )
  {
    LALSnprintf( fileName, FILENAME_MAX, "%s-THINCA_%s-%d-%d.xml", ifos,
        ifoTag, startCoincidence, endCoincidence - startCoincidence );
    LALSnprintf( fileSlide, FILENAME_MAX, "%s-THINCA_SLIDE_%s-%d-%d.xml", ifos,
        ifoTag, startCoincidence, endCoincidence - startCoincidence );
  }
  else if ( userTag )
  {
    LALSnprintf( fileName, FILENAME_MAX, "%s-THINCA_%s-%d-%d.xml", 
        ifos, userTag, startCoincidence, endCoincidence - startCoincidence );
    LALSnprintf( fileSlide, FILENAME_MAX, "%s-THINCA_SLIDE_%s-%d-%d.xml", 
        ifos, userTag, startCoincidence, endCoincidence - startCoincidence );
  }
  else
  {
    LALSnprintf( fileName, FILENAME_MAX, "%s-THINCA-%d-%d.xml", ifos,
        startCoincidence, endCoincidence - startCoincidence );
    LALSnprintf( fileSlide, FILENAME_MAX, "%s-THINCA_SLIDE-%d-%d.xml", ifos,
        startCoincidence, endCoincidence - startCoincidence );
  }
  searchsumm.searchSummaryTable->nevents = numCoinc;

  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );

  if ( !numSlides )
  {
    LAL_CALL( LALOpenLIGOLwXMLFile( &status , &xmlStream, fileName), 
        &status );
  }
  else
  {
    LAL_CALL( LALOpenLIGOLwXMLFile( &status , &xmlStream, fileSlide), 
        &status );
  }
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
    LAL_CALL( LALFreeSnglInspiral( &status, &thisInspiralTrigger ), &status );
  }

  while ( snglOutput )
  {
    thisInspiralTrigger = snglOutput;
    snglOutput = snglOutput->next;
    LAL_CALL( LALFreeSnglInspiral( &status, &thisInspiralTrigger ), &status );
  }
 
  while ( coincInspiralList )
  {
    thisCoinc = coincInspiralList;
    coincInspiralList = coincInspiralList->next;
    LALFree( thisCoinc );
  }


  if ( userTag ) free( userTag );
  if ( ifoTag ) free( ifoTag );

  if ( vrbflg ) fprintf( stdout, "done\n" );

  LALCheckMemoryLeaks();

  exit( 0 );
}
