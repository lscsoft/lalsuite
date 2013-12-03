/*
*  Copyright (C) 2007 Alexander Dietz, Stephen Fairhurst, Sean Seader
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/*-----------------------------------------------------------------------
 *
 * File Name: coherent_bank.c
 *
 * Author: Fairhust, S. and Seader, S.E.
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
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/LIGOMetadataInspiralUtils.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/Segments.h>
#include <lal/SegmentsIO.h>
#include <lalapps.h>
#include <processtable.h>
#include <LALAppsVCSInfo.h>

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "lalapps_coherentbank"

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
  this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
calloc( 1, sizeof(ProcessParamsTable) ); \
snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
    PROGRAM_NAME ); \
snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
    long_options[option_index].name ); \
snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

extern int vrbflg;
int allIFO = -1;
int doVeto = 0;

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
      " [--help]                      display this message\n"\
      " [--verbose]                   print progress information\n"\
      " [--version]                   print version information and exit\n"\
      " [--user-tag]      usertag     set the process_params usertag\n"\
      " [--comment]       string      set the process table comment\n"\
      " [--write-compress]            write a compressed xml file\n"\
      "  --ifos           ifos        list of ifos for which we have data\n"\
      "  --gps-start-time start_time  start time of the job\n"\
      "  --gps-end-time   end_time    end time of the job\n"\
      "  --enable-all-ifo             generate bank with templates for all ifos\n"\
      "  --disable-all-ifo            only generate bank for triggers in coinc\n"\
      " [--coinc-stat]        stat     use coinc statistic for cluster/cut\n"\
      " [--sngl-stat]   clusterchoice use single-ifo statistic for cluster/cut\n"\
      "                     [ none (default) | snr_and_chisq | snrsq_over_chisq | snr ]\n"\
      "                     [ snrsq | effective_snrsq ]\n"\
      " [--run-type]      runType     create trigger bank on coincs or\n"\
      "                               inspiral-coherent triggers\n"\
      "                               [ cohbank (default) | cohinspbank ]\n"\
      " [--stat-threshold]    thresh   discard all triggers with stat less than thresh\n"\
      " [--cluster-time]      time     cluster triggers with time ms window\n" \
      "  [--g1-slide]      g1_slide    Slide G1 data by multiples of g1_slide\n" \
      "  [--h1-slide]      h1_slide    Slide H1 data by multiples of h1_slide\n"\
      "  [--h2-slide]      h2_slide    Slide H2 data by multiples of h2_slide\n"\
      "  [--l1-slide]      l1_slide    Slide L1 data by multiples of l1_slide\n"\
      "  [--t1-slide]      t1_slide    Slide T1 data by multiples of t1_slide\n"\
      "  [--v1-slide]      v1_slide    Slide V1 data by multiples of v1_slide\n"\
      "  [--do-veto]                   Perform a veto on single IFO triggers\n"\
      "                                at the times specified in the veto files below\n"\
      "  [--g1-veto-file]              Veto file for G1\n"\
      "  [--h1-veto-file]              Veto file for H1\n"\
      "  [--h2-veto-file]              Veto file for H2\n"\
      "  [--l1-veto-file]              Veto file for L1\n"\
      "  [--t1-veto-file]              Veto file for T1\n"\
      "  [--v1-veto-file]              Veto file for V1\n"\
      "\n");
}

int main( int argc, char *argv[] )
{
  static LALStatus      status;

  INT4 i;
  INT4 numTriggers = 0;
  INT4 numCoincs = 0;
  INT4 numTmplts = 0;
  INT4 numCoincSegCutTrigs = 0;
  CoincInspiralStatistic coincstat = no_stat;
  SnglInspiralClusterChoice clusterchoice = none;
  CohbankRunType runType = cohbank;
  REAL4 statThreshold = 0;
  INT8 cluster_dt = 0;
  int  numSlides = 0;
  REAL8  slideStep[LAL_NUM_IFO] = {0.0,0.0,0.0,0.0,0.0,0.0};
  LIGOTimeGPS slideTimes[LAL_NUM_IFO];
  CoincInspiralStatParams    bittenLParams;
  INT4        startTime = -1;
  LIGOTimeGPS startTimeGPS = {0,0};
  INT4        endTime = -1;
  LIGOTimeGPS endTimeGPS = {0,0};
  INT8        startTimeNS = 0;
  INT8        endTimeNS = 0;
  CHAR  ifos[LIGOMETA_IFOS_MAX];
  CHAR *vetoFileName[LAL_NUM_IFO] = {NULL, NULL, NULL, NULL, NULL, NULL};
  CHAR  comment[LIGOMETA_COMMENT_MAX];
  CHAR *userTag = NULL;

  CHAR  fileName[FILENAME_MAX];
  CHAR  ifo[LIGOMETA_IFO_MAX];
  LALSegList vetoSegs[LAL_NUM_IFO];
  SnglInspiralTable    *inspiralEventList=NULL;
  SnglInspiralTable    *currentTrigger = NULL;
  SnglInspiralTable    *newEventList = NULL;

  CoincInspiralTable   *coincHead = NULL;
  CoincInspiralTable   *thisCoinc = NULL;

  SearchSummvarsTable  *inputFiles = NULL;
  SearchSummvarsTable  *thisInputFile = NULL;

  SearchSummaryTable   *searchSummList = NULL;
  SearchSummaryTable   *thisSearchSumm = NULL;

  MetadataTable         proctable;
  MetadataTable         processParamsTable;
  MetadataTable         searchsumm;
  MetadataTable         inspiralTable;
  ProcessParamsTable   *this_proc_param = NULL;
  LIGOLwXMLStream       xmlStream;
  INT4                  outCompress = 0;
  InterferometerNumber  ifoNumber = LAL_UNKNOWN_IFO;

  /* getopt arguments */
  struct option long_options[] =
  {
    {"verbose",                no_argument,     &vrbflg,                  1 },
    {"enable-all-ifo",         no_argument,     &allIFO,                  1 },
    {"disable-all-ifo",        no_argument,     &allIFO,                  0 },
    {"write-compress",         no_argument,     &outCompress,             1 },
    {"comment",                required_argument,     0,                 'x'},
    {"user-tag",               required_argument,     0,                 'Z'},
    {"help",                   no_argument,           0,                 'h'},
    {"version",                no_argument,           0,                 'V'},
    {"gps-start-time",         required_argument,     0,                 's'},
    {"gps-end-time",           required_argument,     0,                 't'},
    {"ifos",                   required_argument,     0,                 'i'},
    {"coinc-stat",             required_argument,     0,                 'C'},
    {"sngl-stat",              required_argument,     0,                 'S'},
    {"run-type",               required_argument,     0,                 'r'},
    {"stat-threshold",         required_argument,     0,                 'E'},
    {"cluster-time",           required_argument,     0,                 'T'},
    {"num-slides",             required_argument,     0,                 'N'},
    {"g1-slide",               required_argument,     0,                 'g'},
    {"h1-slide",               required_argument,     0,                 'W'},
    {"h2-slide",               required_argument,     0,                 'X'},
    {"l1-slide",               required_argument,     0,                 'Y'},
    {"t1-slide",               required_argument,     0,                 'U'},
    {"v1-slide",               required_argument,     0,                 'v'},
    {"do-veto",                no_argument,     &doVeto,                  1 },
    {"h1-veto-file",           required_argument,     0,                 '('},
    {"h2-veto-file",           required_argument,     0,                 ')'},
    {"g1-veto-file",           required_argument,     0,                 '{'},
    {"l1-veto-file",           required_argument,     0,                 '}'},
    {"t1-veto-file",           required_argument,     0,                 '['},
    {"v1-veto-file",           required_argument,     0,                 ']'},
    {"eff-snr-denom-fac",      required_argument,     0,                 'A'},
    {0, 0, 0, 0}
  };
  int c;

  /*
   *
   * initialize things
   *
   */

  lal_errhandler = LAL_ERR_EXIT;
  /*  setvbuf( stdout, NULL, _IONBF, 0 );*/

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
  XLALGPSTimeNow(&(proctable.processTable->start_time));
  XLALPopulateProcessTable(proctable.processTable, PROGRAM_NAME, LALAPPS_VCS_IDENT_ID,
      LALAPPS_VCS_IDENT_STATUS, LALAPPS_VCS_IDENT_DATE, 0);
  this_proc_param = processParamsTable.processParamsTable =
    (ProcessParamsTable *) calloc( 1, sizeof(ProcessParamsTable) );

  /* initialize variables */
  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );
  memset( ifos, 0, LIGOMETA_IFOS_MAX * sizeof(CHAR) );
  memset( &slideTimes, 0, LAL_NUM_IFO * sizeof(LIGOTimeGPS) );

  /* create the search summary and zero out the summvars table */
  searchsumm.searchSummaryTable = (SearchSummaryTable *)
    calloc( 1, sizeof(SearchSummaryTable) );

  memset( &bittenLParams, 0, sizeof(CoincInspiralStatParams   ) );
  /* FIXME: hard-wired default value from traditional effective snr formula */
  bittenLParams.eff_snr_denom_fac = 250.0;

  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    size_t optarg_len;
    long int gpstime;

    c = getopt_long_only( argc, argv,
        "hi:r:s:t:x:C:E:T:N:A:W:X:Y:U:v:(:):{:}:[:]:VZ:", long_options,
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

      case 'i':
        /* set ifos */
        strncpy( ifos, optarg, LIGOMETA_IFOS_MAX );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
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
        startTime = (INT4) gpstime;
        startTimeGPS.gpsSeconds = startTime;
        ADD_PROCESS_PARAM( "int", "%" LAL_INT4_FORMAT, startTime );
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
        endTime = (INT4) gpstime;
        endTimeGPS.gpsSeconds = endTime;
        ADD_PROCESS_PARAM( "int", "%" LAL_INT4_FORMAT, endTime );
        break;

      case 'x':
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

      case 'h':
        /* help message */
        print_usage(argv[0]);
        exit( 1 );
        break;

      case 'C':
        /* choose the coinc statistic */
        {
          if ( ! strcmp( "snrsq", optarg ) )
          {
            coincstat = snrsq;
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
                "snrsq, effective_snrsq\n",
                long_options[option_index].name, optarg);
            exit( 1 );
          }
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        break;

      case 'S':
        /* choose the single-ifo cluster statistic */
        {
          if ( ! strcmp( "none", optarg ) )
          {
            clusterchoice = none;
          }
          else if ( ! strcmp( "snr", optarg ) )
          {
            clusterchoice = snr;
          }
          else if ( ! strcmp( "snr_and_chisq", optarg) )
          {
            clusterchoice = snr_and_chisq;
          }
          else if ( ! strcmp( "snrsq_over_chisq", optarg) )
          {
            clusterchoice = snrsq_over_chisq;
          }
          else
          {
            fprintf( stderr, "invalid argument to  --%s:\n"
                "unknown coinc statistic:\n "
                "%s (must be one of:\n"
		"none, snr, snr_and_chisq, or snrsq_over_chisq\n",
                long_options[option_index].name, optarg);
            exit( 1 );
          }
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        break;

      case 'r':
        /* choose the run-type for constructing a bank on either
           coinc or inspiral-coherent triggers */
        {
          if ( ! strcmp( "cohbank", optarg ) )
          {
            runType = cohbank;
          }
          else if ( ! strcmp( "cohinspbank", optarg) )
          {
            runType = cohinspbank;
          }
          else
          {
            fprintf( stderr, "invalid argument to  --%s:\n"
                "unknown run-type:\n "
                "%s (must be one of:\n"
                "cohbank, cohinspbank\n",
                long_options[option_index].name, optarg);
            exit( 1 );
          }
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        break;

      case 'E':
        /* store the stat threshold for a cut */
        statThreshold = atof( optarg );
        if ( statThreshold <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "statThreshold must be positive: (%f specified)\n",
              long_options[option_index].name, statThreshold );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%f", statThreshold );
        break;


      case 'T':
        /* cluster time is specified on command line in ms */
        cluster_dt = (INT8) atoi( optarg );
        if ( cluster_dt <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "custer window must be > 0: "
              "(%" LAL_INT8_FORMAT " specified)\n",
              long_options[option_index].name, cluster_dt );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%" LAL_INT8_FORMAT, cluster_dt );
        /* convert cluster time from ms to ns */
        cluster_dt *= 1000000LL;
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

      case 'A':
        bittenLParams.eff_snr_denom_fac = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg);
        break;

      /* Read in time-slide steps for all detectors */
      /* Read in time-slide step for G1 */
      case 'g':
        slideStep[LAL_IFO_G1] = (REAL8) atof(optarg);
        ADD_PROCESS_PARAM("float", "%e", slideStep[LAL_IFO_G1] );
        break;

      /* Read in time-slide step for H1 */
      case 'W':
        slideStep[LAL_IFO_H1] = (REAL8) atof(optarg);
        ADD_PROCESS_PARAM("float", "%e", slideStep[LAL_IFO_H1]);
        break;

      /* Read in time-slide step for H2 */
      case 'X':
        slideStep[LAL_IFO_H2] = (REAL8) atof(optarg);
        ADD_PROCESS_PARAM("float", "%e", slideStep[LAL_IFO_H2]);
        break;

      /* Read in time-slide step for L1 */
      case 'Y':
        slideStep[LAL_IFO_L1] = (REAL8) atof(optarg);
        ADD_PROCESS_PARAM("float", "%e", slideStep[LAL_IFO_L1]);
        break;

      /* Read in time-slide step for T1 */
      case 'U':
        slideStep[LAL_IFO_T1] = (REAL8) atof(optarg);
        ADD_PROCESS_PARAM("float", "%e", slideStep[LAL_IFO_T1]);
        break;

      /* Read in time-slide step for V1 */
      case 'v':
        slideStep[LAL_IFO_V1] = (REAL8) atof(optarg);
        ADD_PROCESS_PARAM("float", "%e", slideStep[LAL_IFO_V1]);
        break;

      case 'Z':
        /* create storage for the usertag */
        optarg_len = strlen(optarg) + 1;
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

      case '(':
        /* veto filename */
        optarg_len = strlen( optarg ) + 1;
        vetoFileName[LAL_IFO_H1] = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( vetoFileName[LAL_IFO_H1], optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case ')':
        /* veto filename */
        optarg_len = strlen( optarg ) + 1;
        vetoFileName[LAL_IFO_H2] = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( vetoFileName[LAL_IFO_H2], optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case '}':
        /* veto filename */
        optarg_len = strlen( optarg ) + 1;
        vetoFileName[LAL_IFO_L1] = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( vetoFileName[LAL_IFO_L1], optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case '{':
        /* veto filename */
        optarg_len = strlen( optarg ) + 1;
        vetoFileName[LAL_IFO_G1] = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( vetoFileName[LAL_IFO_G1], optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case '[':
        /* veto filename */
        optarg_len = strlen( optarg ) + 1;
        vetoFileName[LAL_IFO_T1] = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( vetoFileName[LAL_IFO_T1], optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case ']':
        /* veto filename */
        optarg_len = strlen( optarg ) + 1;
        vetoFileName[LAL_IFO_V1] = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( vetoFileName[LAL_IFO_V1], optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'V':
        /* print version information and exit */
        fprintf( stdout, "Coherent Bank Generator\n"
            "Steve Fairhurst and Shawn Seader\n");
        XLALOutputVersionString(stderr, 0);
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
  if ( ! *comment )
  {
    snprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
    snprintf( searchsumm.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX,
        " " );
  }
  else
  {
    snprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
    snprintf( searchsumm.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX,
        "%s", comment );
  }

  /* enable/disable-all-ifo is stored in the first process param row */
  if ( allIFO == 1 )
  {
    snprintf( processParamsTable.processParamsTable->program,
        LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
    snprintf( processParamsTable.processParamsTable->param,
        LIGOMETA_PARAM_MAX, "--enable-all-ifo" );
    snprintf( processParamsTable.processParamsTable->type,
        LIGOMETA_TYPE_MAX, "string" );
    snprintf( processParamsTable.processParamsTable->value,
        LIGOMETA_TYPE_MAX, " " );
  }
  else if ( allIFO == 0 )
  {
    snprintf( processParamsTable.processParamsTable->program,
        LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
    snprintf( processParamsTable.processParamsTable->param,
        LIGOMETA_PARAM_MAX, "--disable-all-ifo" );
    snprintf( processParamsTable.processParamsTable->type,
        LIGOMETA_TYPE_MAX, "string" );
    snprintf( processParamsTable.processParamsTable->value,
        LIGOMETA_TYPE_MAX, " " );
  }
  else
  {
    fprintf( stderr, "--enable-all-ifo or --disable-all-ifo "
        "argument must be specified\n" );
    exit( 1 );
  }

  /* store the veto option */
  if ( doVeto )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--do-veto" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }

  /*
   *
   * check the values of the arguments
   *
   */

  startTimeNS = XLALGPSToINT8NS( &startTimeGPS );
  endTimeNS = XLALGPSToINT8NS( &endTimeGPS );

  if ( startTime < 0 )
  {
    fprintf( stderr, "Error: --gps-start-time must be specified\n" );
    exit( 1 );
  }

  if ( endTime < 0 )
  {
    fprintf( stderr, "Error: --gps-end-time must be specified\n" );
    exit( 1 );
  }

  if ( !strlen(ifos) )
  {
    fprintf(stderr,"You must specify a list of ifos with --ifos. Exiting.\n");
    exit(1);
  }

  /* check that if clustering is being done that we have all the options */
  if ( cluster_dt && (coincstat == no_stat) )
  {
    fprintf( stderr,
        "--coinc-stat must be specified if --cluster-time is given\n" );
    exit( 1 );
  }

  /*
   *
   * read in the input data from the rest of the arguments
   *
   */

  if ( optind < argc )
  {
    if ( !(runType == cohinspbank ) ) {
      for( i = optind; i < argc; ++i )
      {
	INT4 numFileTriggers = 0;
	numFileTriggers = XLALReadInspiralTriggerFile( &inspiralEventList,
	     &currentTrigger, &searchSummList, &inputFiles, argv[i] );
	if (numFileTriggers < 0)
	{
	  fprintf(stderr, "Error reading triggers from file %s",
            argv[i]);
	  exit( 1 );
	}

	numTriggers += numFileTriggers;
      }
    }
    else {
      InterferometerNumber  ifoNumberTmp = LAL_UNKNOWN_IFO;
      INT4                  numCoincTrigs = 0;
      for ( ifoNumberTmp = 0; ifoNumberTmp< LAL_NUM_IFO; ifoNumberTmp++) {
	INT4                  numIfoTriggers = 0;

	for( i = optind; i < argc; ++i )
	{
	  INT4 numFileTriggers = 0;
	  SearchSummaryTable *inputSummary = NULL;

	  /* read in the search summary and store */
	  XLALPrintInfo(
		"XLALReadInspiralTriggerFile(): Reading search_summary table\n");

	  inputSummary = XLALSearchSummaryTableFromLIGOLw(argv[i]);

	  if ( ! inputSummary )
	  {
	    fprintf(stderr,"No valid search_summary table in %s, exiting\n",
		    argv[i] );
	    exit( 1 );
	  }
	  else
	  {
	    fprintf(stdout,"IFOs in input summary table is %s\n",inputSummary->ifos);
	    if (ifoNumberTmp == XLALIFONumber(inputSummary->ifos)){
	      fprintf(stdout,"Reading triggers from IFO %s\n",inputSummary->ifos);
	      numFileTriggers = XLALReadInspiralTriggerFile( &inspiralEventList,
			&currentTrigger, &searchSummList, &inputFiles, argv[i] );
	      if (numFileTriggers < 0)
	      {
		fprintf(stderr, "Error reading triggers from file %s",
			argv[i]);
		exit( 1 );
	      }

	      numIfoTriggers += numFileTriggers;

	    }
	  }
	}/* Loop to read all trigger files from a single ifo */

      	if( vrbflg )
	  {
	    fprintf( stdout,
		     "Number of triggers found in ifo number %d is %d\n",
		     ifoNumberTmp , numIfoTriggers);
	  }

        numCoincTrigs += numIfoTriggers;
        if( vrbflg )
        {
          fprintf( stdout,
                   "Number of unclustered triggers found in all ifos in this coinc-segment is %d\n",
                  numCoincTrigs);
        }
      }/* Closes for loop over ifoNumberTmp */

      /* keep only triggers within the requested interval */
      if ( vrbflg ) fprintf( stdout,
			     "Discarding triggers outside requested interval\n" );
      LAL_CALL( LALTimeCutSingleInspiral( &status, &inspiralEventList,
				  &startTimeGPS, &endTimeGPS), &status );

      if ( vrbflg ) fprintf( stdout,
             "Removing triggers with event-ids from outside this coinc-segment\n" );
      if ( vrbflg ) fprintf( stdout, "GPS start time of this coinc-segment is %d\n",
                   startTime);

      numCoincSegCutTrigs = XLALCoincSegCutSnglInspiral( startTime,
                                        endTime, &inspiralEventList);

      if ( vrbflg ) fprintf( stdout,
                             "Sorting triggers within requested interval\n" );
      /* sort single inspiral trigger list according to event_id */
      inspiralEventList = XLALSortSnglInspiral( inspiralEventList,
                                        LALCompareSnglInspiralByID);

      if ( vrbflg ) fprintf( stdout,
                             "Clustering triggers in event-id\n" );
      if ( ! ( clusterchoice == none ) ) {
        numTriggers = XLALClusterInEventID(&inspiralEventList,clusterchoice);
      }
      else {
        numTriggers = numCoincSegCutTrigs;
      }

      if( vrbflg )
      {
	fprintf( stdout,
		   "Number of CLUSTERED triggers found in coinc-segment is %d\n",
		   numTriggers);
      }
    } /*Closes if runType is not cohinspbank */
  } /* Closes if optind < argc */
  else
  {
    fprintf( stderr, "Error: No trigger files specified.\n" );
    exit( 1 );
  }

  if ( numTriggers == 0 )
  {
    if( vrbflg )
    {
      fprintf( stdout,
         "No triggers found - the coherent bank will be empty.\n");
    }
  }
  else
  {
    if( vrbflg )
    {
      fprintf( stdout,
          "Read in a total of %d triggers.\n", numTriggers);
    }

    /* reconstruct the coincs */
    numCoincs = XLALRecreateCoincFromSngls( &coincHead, &inspiralEventList );
    if( numCoincs < 0 )
    {
      fprintf(stderr, "Unable to reconstruct coincs from single ifo triggers");
      exit( 1 );
    }
    else if ( vrbflg )
    {
      fprintf( stdout,
          "Recreated %d coincs from the %d triggers\n", numCoincs,
          numTriggers );
    }

    /*
     *
     * sort the inspiral events by time
     *
     */


    if ( coincHead && cluster_dt )
    {
      if ( vrbflg ) fprintf( stdout, "sorting coinc inspiral trigger list..." );
      coincHead = XLALSortCoincInspiral( coincHead,
                          *XLALCompareCoincInspiralByTime );
      if ( vrbflg ) fprintf( stdout, "done\n" );

      if ( vrbflg ) fprintf( stdout, "clustering remaining triggers... " );

      if ( numSlides )
      {
        int slide = 0;
        CoincInspiralTable *slideCoinc = NULL;
        CoincInspiralTable *slideClust = NULL;

        if ( vrbflg ) fprintf( stdout, "splitting events by slide\n" );

        for( slide = -numSlides; slide < (numSlides + 1); slide++)
        {
          if ( vrbflg ) fprintf( stdout, "slide number %d; ", slide );
          /* extract the slide */
          slideCoinc = XLALCoincInspiralSlideCut( &coincHead, slide );

          /* add clustered triggers */
          if( slideCoinc )
          {
            if( slideClust )
            {
              thisCoinc = thisCoinc->next = slideCoinc;
            }
            else
            {
              slideClust = thisCoinc = slideCoinc;
            }
            /* scroll to end of list */
            for( ; thisCoinc->next; thisCoinc = thisCoinc->next);
          }
        }

        /* free coincHead -- although we expect it to be empty */
        while ( coincHead )
        {
          thisCoinc = coincHead;
          coincHead = coincHead->next;
          XLALFreeCoincInspiral( &thisCoinc );
        }

        /* move events to coincHead */
        coincHead = slideClust;
        slideClust = NULL;
      }

      if ( vrbflg ) fprintf( stdout, "done\n" );

    }

    /*
     *
     *  Create the coherent bank
     *
     */

    numTmplts = XLALGenerateCoherentBank( &newEventList, coincHead, runType,
                    startTimeNS, endTimeNS, numSlides, slideStep, statThreshold, ifos);

    if ( numTmplts < 0 )
    {
      fprintf(stderr, "Unable to generate coherent bank\n");
      exit( 1 );
    }

    /* do a veto on the new sngls */
    for ( ifoNumber = 0; ifoNumber< LAL_NUM_IFO; ifoNumber++)
    {
      if (doVeto && vetoFileName[ifoNumber])
      {

	/* CHECK: if this is required.
	 * If the veto segment list wasn't initialized, then don't try to slide
	 * it.  The rest of the code, except possibly the H1H2 consistency test,
	 * does not try to use vetoSegs if it wasn't loaded / initialized. */
	/* if ( vetoSegs[ifoNumber].initMagic == SEGMENTSH_INITMAGICVAL )
	   {
	   XLALTimeSlideSegList( &vetoSegs[ifoNumber], &startCoinc, &endCoinc,
	   &slideTimes[ifoNumber] );
	   }*/

	XLALReturnIFO(ifo,ifoNumber);
	XLALSegListInit( &(vetoSegs[ifoNumber]) );
	LAL_CALL( LALSegListRead( &status, &(vetoSegs[ifoNumber]),
				  vetoFileName[ifoNumber], NULL),&status);
	XLALSegListCoalesce( &(vetoSegs[ifoNumber]) );

	/* keep only the segments that lie within the data-segment part */
	XLALSegListKeep(  &(vetoSegs[ifoNumber]), &startTimeGPS, &endTimeGPS );

	if ( vrbflg ) fprintf( stdout,
            "Applying veto segment (%s) list on ifo  %s \n ",
            vetoFileName[ifoNumber], ifo );
	newEventList = XLALVetoSingleInspiral( newEventList,
					       &(vetoSegs[ifoNumber]), ifo );
      }
    }

    /* count remaining singles */
    numTmplts =  XLALCountSnglInspiral( newEventList );
    if ( vrbflg )
    {
      fprintf(stdout, "Generated a coherent bank with %d templates\n",
	      numTmplts);
    }
  }

  /*
   *
   * write the output xml file
   *
   */

  /* search summary entries: */
  searchsumm.searchSummaryTable->in_start_time = startTimeGPS;
  searchsumm.searchSummaryTable->in_end_time = endTimeGPS;
  searchsumm.searchSummaryTable->out_start_time = startTimeGPS;
  searchsumm.searchSummaryTable->out_end_time = endTimeGPS;
  searchsumm.searchSummaryTable->nevents = numTmplts;

  if ( vrbflg ) fprintf( stdout, "writing output file... " );

  /* set the file name correctly */
  if ( userTag && !outCompress )
  {
    if ( runType == cohinspbank ) {
      snprintf( fileName, FILENAME_MAX, "%s-COHINSPBANK_%s-%d-%d.xml",
        ifos, userTag, startTime, endTime - startTime );
    }
    else {
      snprintf( fileName, FILENAME_MAX, "%s-COHBANK_%s-%d-%d.xml",
        ifos, userTag, startTime, endTime - startTime );
    }
  }
  else if ( userTag && outCompress )
  {
    if ( runType == cohinspbank ) {
      snprintf( fileName, FILENAME_MAX, "%s-COHINSPBANK_%s-%d-%d.xml.gz",
        ifos, userTag, startTime, endTime - startTime );
    }
    else {
      snprintf( fileName, FILENAME_MAX, "%s-COHBANK_%s-%d-%d.xml.gz",
        ifos, userTag, startTime, endTime - startTime );
    }
  }
  else if ( !userTag && outCompress )
  {
    if ( runType == cohinspbank ) {
      snprintf( fileName, FILENAME_MAX, "%s-COHINSPBANK-%d-%d.xml.gz",
        ifos, startTime, endTime - startTime );
    }
    else {
      snprintf( fileName, FILENAME_MAX, "%s-COHBANK-%d-%d.xml.gz",
        ifos, startTime, endTime - startTime );
    }
  }
  else
  {
    if ( runType == cohinspbank ) {
      snprintf( fileName, FILENAME_MAX, "%s-COHINSPBANK-%d-%d.xml",
        ifos, startTime, endTime - startTime );
    }
    else {
      snprintf( fileName, FILENAME_MAX, "%s-COHBANK-%d-%d.xml",
        ifos, startTime, endTime - startTime );
    }
  }
  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status , &xmlStream, fileName ),
      &status );

  /* write process table */
  XLALGPSTimeNow(&(proctable.processTable->end_time));
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
  /* XXX not necessary as bank file specified in arguements XXX
  LAL_CALL( LALBeginLIGOLwXMLTable( &status ,&xmlStream,
        search_summvars_table), &status );
  searchSummvarsTable.searchSummvarsTable = inputFiles;
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, searchSummvarsTable,
        search_summvars_table), &status );
  LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream), &status ); */

  /* write the sngl_inspiral table if we have one*/
  if ( newEventList )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status ,&xmlStream,
          sngl_inspiral_table), &status );
    inspiralTable.snglInspiralTable = newEventList;
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

  /* free the veto segment list. */
  for (ifoNumber=0; ifoNumber<LAL_NUM_IFO; ifoNumber++)
  {
    if ( vetoFileName[ifoNumber] )
    {
      free( vetoFileName[ifoNumber] );
    }

   if (vetoSegs[ifoNumber].initMagic == SEGMENTSH_INITMAGICVAL )
   {
        XLALSegListClear( &vetoSegs[ifoNumber] );
    }
  }

  while ( inspiralEventList )
  {
    currentTrigger = inspiralEventList;
    inspiralEventList = inspiralEventList->next;
    LAL_CALL( LALFreeSnglInspiral( &status, &currentTrigger ), &status );
  }

  while ( newEventList )
  {
    currentTrigger = newEventList;
    newEventList = newEventList->next;
    LAL_CALL( LALFreeSnglInspiral( &status, &currentTrigger ), &status );
  }

  while ( coincHead )
  {
    thisCoinc = coincHead;
    coincHead = coincHead->next;
    LALFree( thisCoinc );
  }

  if ( userTag ) free( userTag );

  if ( vrbflg ) fprintf( stdout, "done\n" );

  LALCheckMemoryLeaks();

  exit( 0 );
}
