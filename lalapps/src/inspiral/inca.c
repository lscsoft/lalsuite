/*
*  Copyright (C) 2007 Duncan Brown, Eirini Messaritaki, Kipp Cannon, Stephen Fairhurst
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
#include <lal/lalGitID.h>
#include <lalappsGitID.h>
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
#define IFOB_SNRTHRESH 6
#define KAPPA 0.01
#define EPSILON 2

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
"  --gps-start-time SEC      GPS second of data start time\n"\
"  --gps-end-time SEC        GPS second of data end time\n"\
"\n"\
"  --silde-time SEC          slide all triggers of IFOB by SEC\n"\
"  --slide-time-ns NS        slide all triggers of IFOB by NS\n"\
"\n"\
"  --ifo-a IFOA              name of first ifo (e.g. L1, H1 or H2)\n"\
"  --ifo-b IFOB              name of second ifo (e.g. L1, H1 or H2)\n"\
"\n"\
"  --single-ifo              input triggers from only one IFO\n"\
"  --single-summ-value       only write the first summ_value table found\n"\
"  --triggered-bank FILE     write a triggered bank insted of doing inca\n"\
"  --minimal-match M         set minimal match of triggered bank to M\n"\
"\n"\
"  --epsilon ERROR           set effective distance test epsilon (default 2)\n"\
"  --kappa ERROR             set effective distance test kappa (default 0.01)\n"\
"  --ifo-b-snr-threshold SNR set minimum snr in IFO B (default 6)\n"\
"  --ifo-b-range-cut         test range of IFO B to see if sensitive to trigger\n"\
"  --parameter-test TEST    set the desired parameters to test coincidence\n"\
"                            for inca: (m1_and_m2|psi0_and_psi3|mchirp_and_eta)\n"\
"                            for triggered bank (m1_and_m2|psi0_and_psi3)\n"\
"  --dm Dm                   mass coincidence window (default 0)\n"\
"  --dpsi0 Dpsi0             psi0 coincidence window\n"\
"  --dpsi3 Dpsi3             psi3 coincidence window\n"\
"  --alphaf-cut AlphaFCut    ignore BCV trigs with alphaF > AlphaFCut for\n"\
"                            the coincidence step\n"\
"  --dmchirp Dmchirp         mchirp coincidence window\n"\
"  --deta  Deta              eta coincidence window\n"\
"  --dt Dt                   time coincidence window (milliseconds)\n"\
"\n"\
"  --no-playground           do not select triggers from playground\n"\
"  --playground-only         only use triggers that are in playground\n"\
"  --all-data                use all triggers\n"\
"  --write-uniq-triggers     make sure triggers from IFO A are unique\n" \
"  --write-compress          write a compressed xml file\n" \
"\n"\
"[LIGOLW XML input files] list of the input trigger files.\n"\
"\n"

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
  this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
  calloc( 1, sizeof(ProcessParamsTable) ); \
  snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
      PROGRAM_NAME ); \
  snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
      long_options[option_index].name ); \
  snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
  snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

int main( int argc, char *argv[] )
{
  static LALStatus      status;

  extern int vrbflg;
  static INT4  writeUniqTrigs = 0;
  static INT4  usePlayground = 1;
  static INT4  allData = 0;

  INT4  havePlgOpt = 0;
  INT4  startCoincidence = -1;
  INT4  endCoincidence = -1;
  CHAR  ifoName[MAXIFO][LIGOMETA_IFO_MAX];
  CHAR  comment[LIGOMETA_COMMENT_MAX];
  CHAR *userTag = NULL;
  CHAR *ifoTag = NULL;

  CHAR  fileName[FILENAME_MAX];
  CHAR *trigBankFile = NULL;
  CHAR *xmlFileName;
  UINT4 outCompress = 0;

  LIGOTimeGPS slideData = {0,0};
  INT8  slideDataNS = 0;
  INT4  numEvents = 0;
  INT4  numIFO;
  INT4  have_ifo_a_trigger = 0;
  INT4  keep_a_trig = 0;
  INT4  dont_search_b = 0;
  INT4  isPlay = 0;
  INT8  ta, tb;
  INT4  numTriggers[MAXIFO];
  INT4  inStartTime = -1;
  INT4  inEndTime = -1;
  REAL4 minMatch = -1;
  INT4  useRangeCut = 0;
  INT4  singleIfo = 0;
  INT4  singleSummValue = 0;
  REAL4 ifob_snrthresh = IFOB_SNRTHRESH;
  REAL4 d_range[MAXIFO];

  SnglInspiralTable    *inspiralEventList[MAXIFO];
  SnglInspiralTable    *currentTrigger[MAXIFO];

  SnglInspiralTable    *currentEvent = NULL;
  SnglInspiralTable    *outEvent[MAXIFO];
  SnglInspiralTable    *coincidentEvents[MAXIFO];
  SnglInspiralAccuracy  errorParams;
  REAL4                 alphaFcut = -1.0;

  SummValueTable       *inspEffRange[MAXIFO];
  SummValueTable       *currentEffRange[MAXIFO];

  SearchSummvarsTable  *inputFiles = NULL;
  SearchSummvarsTable  *thisInputFile = NULL;

  MetadataTable         proctable;
  MetadataTable         processParamsTable;
  MetadataTable         searchsumm;
  MetadataTable         searchSummvarsTable;
  MetadataTable         summValueTable;
  MetadataTable         inspiralTable;
  ProcessParamsTable   *this_proc_param = NULL;
  LIGOLwXMLStream       xmlStream;

  INT4                  i, j;

  /* getopt arguments */
  struct option long_options[] =
  {
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"write-compress",          no_argument,       &outCompress,      1 },
    {"write-uniq-triggers",     no_argument,       &writeUniqTrigs,   1 },
    {"ifo-b-range-cut",         no_argument,       &useRangeCut,      1 },
    {"single-ifo",              no_argument,       &singleIfo,        1 },
    {"single-summ-value",       no_argument,       &singleSummValue,  1 },
    {"no-playground",           no_argument,       0,                'Q'},
    {"playground-only",         no_argument,       0,                'R'},
    {"all-data",                no_argument,       0,                'D'},
    {"ifo-a",                   required_argument, 0,                'a'},
    {"ifo-b",                   required_argument, 0,                'b'},
    {"epsilon",                 required_argument, 0,                'e'},
    {"triggered-bank",          required_argument, 0,                'T'},
    {"minimal-match",           required_argument, 0,                'M'},
    {"kappa",                   required_argument, 0,                'k'},
    {"ifo-b-snr-threshold",     required_argument, 0,                'S'},
    {"dm",                      required_argument, 0,                'm'},
    {"parameter-test",          required_argument, 0,                'A'},
    {"dpsi0",                   required_argument, 0,                'p'},
    {"dpsi3",                   required_argument, 0,                'P'},
    {"alphaf-cut",              required_argument, 0,                'F'},
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
  XLALGPSTimeNow(&(proctable.processTable->start_time));
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
        "a:b:e:k:A:m:p:P:F:t:q:r:s:hz:I:Z:M:T:S:c:n:QRD", long_options, 
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
        snprintf( ifoName[0], LIGOMETA_IFO_MAX, "%s", optarg );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'b':
        /* name of interferometer b */
        snprintf( ifoName[1], LIGOMETA_IFO_MAX, "%s", optarg );
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

      case 'F':
        /* alphaF cut for coincidence */
        alphaFcut = atof(optarg);
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
          snprintf( comment, LIGOMETA_COMMENT_MAX, "%s", optarg);
        }
        break;

      case 'Q':
        usePlayground = 0;
        havePlgOpt = 1;
        break;

      case 'R':
        usePlayground = 1;
        havePlgOpt = 1;
        break;

      case 'D':
        allData = 1;
        usePlayground = 0;
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
        snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", 
            PROGRAM_NAME );
        snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "-userTag" );
        snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
            optarg );
        break;

      case 'I':
        /* create storage for the ifo-tag */
        optarg_len = strlen(optarg) + 1;
        ifoTag = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( ifoTag, optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
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
          fprintf( stderr, "invalid argument to --%s:\n"
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
            "Patrick Brady, Duncan Brown and Steve Fairhurst\n"
            "CVS Version: " CVS_ID_STRING "\n"
            "CVS Tag: " CVS_NAME_STRING "\n" );
        fprintf( stdout, lalappsGitID );
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

  if ( ! haveTest & ! singleIfo )
  {
    fprintf( stderr, "--parameter-test must be specified\n" );
    exit( 1 );
  }

  if ( haveTest & singleIfo )
  {
    fprintf( stderr, 
        "--parameter-test must not be specified in single IFO mode\n" );
    exit( 1 );
  }

  if ( (errorParams.test == psi0_and_psi3) && (alphaFcut < 0.0) )
  {
    fprintf( stderr, 
        "alphaf-cut MUST be specified for BCV triggers!\n");
    exit( 1 );
  }

  /* check for minimal match when doing a triggered bank */
  if ( trigBankFile && minMatch < 0 )
  {
    fprintf( stderr, "--minimal-match must be specified\n" );
    exit( 1 );
  }

  /* check that a playground option is not specified if */
  /* doing a slide or a trig bank                       */
  if ( ( trigBankFile || slideDataNS ) && havePlgOpt )
  {
    fprintf( stderr, "--playground-only or --no-playground should not "
        "be specified for a time slide\n" );
    exit( 1 );
  }

  /* fill the comment, if a user has specified one, or leave it blank */
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

  /* store the write all trigs option */
  if ( writeUniqTrigs )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--write-uniq-triggers" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }

  /* store the ifo b range cut option */
  if ( useRangeCut )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--ifo-b-range-cut" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }

  /* store the single ifo option */
  if ( singleIfo )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--single-ifo" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }

  if ( trigBankFile || slideDataNS )
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
    snprintf( processParamsTable.processParamsTable->program, 
        LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
    snprintf( processParamsTable.processParamsTable->type, 
        LIGOMETA_TYPE_MAX, "string" );
    snprintf( processParamsTable.processParamsTable->value, 
        LIGOMETA_TYPE_MAX, " " );
    if ( usePlayground )
    {
      snprintf( processParamsTable.processParamsTable->param, 
          LIGOMETA_PARAM_MAX, "--playground-only" );
    }
    else if ( !usePlayground && !allData )
    {
      snprintf( processParamsTable.processParamsTable->param, 
          LIGOMETA_PARAM_MAX, "--no-playground" );
    }
    else
    {
      snprintf( processParamsTable.processParamsTable->param, 
          LIGOMETA_PARAM_MAX, "--all-data" );
    }       
  }


  /* decide how many ifos we have based on what we are doing */
  if ( trigBankFile || singleIfo )
  {
    numIFO = 1;
  }
  else
  {
    numIFO = 2;
  }

  /* calculate the slide time in nanoseconds */
  slideDataNS = XLALGPSToINT8NS( &slideData );


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
      SearchSummaryTable *inputSummary = NULL;

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
      snprintf( thisInputFile->name, LIGOMETA_NAME_MAX, 
          "input_file" );
      snprintf( thisInputFile->string, LIGOMETA_NAME_MAX, 
          "%s", argv[i] );      


      if ( vrbflg ) fprintf( stdout, 
          "reading search_summary table from file: %s\n", argv[i] );

      inputSummary = XLALSearchSummaryTableFromLIGOLw( argv[i] );

      if ( ! inputSummary )
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

      if ( ! trigBankFile )
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
#if 0
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
#endif
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
                    if ( slideDataNS && j == 1)
                    {
                      INT8 startNS = 0;
                      INT8 endNS = 0;

                      if ( vrbflg ) fprintf( stdout, 
                          "Doing a time slide of %d sec %d nanosec on IFOB\n",
                          slideData.gpsSeconds, slideData.gpsNanoSeconds );

                      startNS = XLALGPSToINT8NS( &(thisSummValue->start_time) );
                      startNS += slideDataNS;
                      XLALINT8NSToGPS( &(thisSummValue->start_time), startNS );

                      endNS = XLALGPSToINT8NS( &(thisSummValue->end_time) );
                      endNS += slideDataNS;
                      XLALINT8NSToGPS( &(thisSummValue->end_time), endNS );
                      if ( vrbflg ) 
                      {
                        fprintf( stdout, "inspiral effective distance of %f ", 
                            thisSummValue->value );
                        fprintf( stdout, "now valid from %d sec %d nanosec\n",
                            thisSummValue->start_time.gpsSeconds, 
                            thisSummValue->start_time.gpsNanoSeconds);
                        fprintf( stdout, "to %d sec %d nanosec\n",
                            thisSummValue->end_time.gpsSeconds, 
                            thisSummValue->end_time.gpsNanoSeconds);
                      }         
                    } /* close if ( slideDataNS && j == 1) */

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
                  } /* close if ( ! strncmp( ifoName[j], thisSummValue->ifo,
                                               LIGOMETA_IFO_MAX ) ) */
                } /*close for ( j = 0; j < numIFO ; ++j ) */
                if ( ! knownIFO )
                {
                  /* catch an unknown ifo name among the input files */
                  if ( vrbflg ) fprintf( stdout, 
                      "Unknown interferometer %s, discarding\n", 
                      thisSummValue->ifo );
                  tempSummValue = thisSummValue;
                  thisSummValue = thisSummValue->next;
                  LALFree( tempSummValue );
                }
#if 0
              } /* close for second else */
#endif
            } /* close for first else */
          } /* close while ( thisSummValue ) */
        }
      } /* close if( ! trigBankFile ) */

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
                trigTimeNS = XLALGPSToINT8NS( &(currentTrigger[j]->end_time) );
                trigTimeNS += slideDataNS;
                XLALINT8NSToGPS( &(currentTrigger[j]->end_time), trigTimeNS );
              }     
              currentTrigger[j] = currentTrigger[j]->next;
            }

            /* slide the last trigger */
            if ( slideDataNS && j == 1)
            {
              INT8 trigTimeNS = 0;
              trigTimeNS = XLALGPSToINT8NS( &(currentTrigger[j]->end_time) );
              trigTimeNS += slideDataNS;
              XLALINT8NSToGPS( &(currentTrigger[j]->end_time), trigTimeNS );
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
      /* no triggers in this ifo so no coincidences can be found */
     
      fprintf( stdout, "No triggers read in for interferometer %d\n", j );

      if ( j && useRangeCut )
      {
        if (vrbflg)
          fprintf( stdout, "Still need to keep those triggers in ifo 0\n"
              "Which were too distant to be observable in ifo %d\n", j);
      }
      else
      {
        goto cleanexit;
      }
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

    /* sort the templates by time and remove all the tmplts   */
    /* before and after the requested gps start and end times */
    numTriggers[0] = 0;
    thisEvent = inspiralEventList[0];
    inspiralEventList[0] = NULL;

    if ( vrbflg ) fprintf( stdout, 
        "discarding triggers outside start/end times: " );
    while ( thisEvent )
    {
      SnglInspiralTable *tmpEvent = thisEvent;
      thisEvent = thisEvent->next;

      if ( tmpEvent->end_time.gpsSeconds >= startCoincidence &&
          tmpEvent->end_time.gpsSeconds < endCoincidence )
      {
        /* keep this template */
        if ( ! inspiralEventList[0] )
        {
          inspiralEventList[0] = tmpEvent;
        }
        else
        {
          prevEvent->next = tmpEvent;
        }
        tmpEvent->next = NULL;
        prevEvent = tmpEvent;
        ++numTriggers[0];
        if ( vrbflg ) fprintf( stdout, "+" );
      }
      else
      {
        /* discard this template */
        LALFree( tmpEvent );
        if ( vrbflg ) fprintf( stdout, "-" );
      }
    }

    if ( vrbflg ) fprintf( stdout, " done\nkept %d templates\n", 
        numTriggers[0] );

    numEvents = numTriggers[0];

    eventHandle = (SnglInspiralTable **) 
      LALCalloc( numEvents, sizeof(SnglInspiralTable *) );

    for ( i = 0, thisEvent = inspiralEventList[0]; i < numEvents; 
        ++i, thisEvent = thisEvent->next )
    {
      eventHandle[i] = thisEvent;
    }

    if ( errorParams.test == m1_and_m2 )
    {       
      if ( vrbflg ) fprintf( stdout, "sorting events by mass... " );
      qsort( eventHandle, numEvents, sizeof(eventHandle[0]), 
          LALCompareSnglInspiralByMass );
      if ( vrbflg ) fprintf( stdout, "done\n" );
    }
    else if ( errorParams.test == psi0_and_psi3 )
    { 
      if ( vrbflg ) fprintf( stdout, "sorting events by psi... " );
      qsort( eventHandle, numEvents, sizeof(eventHandle[0]),
          LALCompareSnglInspiralByPsi );
      if ( vrbflg ) fprintf( stdout, "done\n" );
    }
    else
    {
      fprintf( stderr, 
          "error: unknown test for sorting events \n" );
      exit( 1 );
    }

    /* create a linked list of sorted templates */
    if ( vrbflg ) fprintf( stdout, 
        "discarding template with duplicate masses: " );

    numTriggers[0] = 0;
    coincidentEvents[0] = prevEvent = eventHandle[0];
    if ( coincidentEvents[0] ) numTriggers[0] = 1;

    for ( i = 1; i < numEvents; ++i )
    {
      if ( errorParams.test == m1_and_m2 )
      {
        if ( (prevEvent->mass1 == eventHandle[i]->mass1)  &&
            (prevEvent->mass2 == eventHandle[i]->mass2) ) 
        {
          /* discard the event as it is a duplicate */
          LALFree( eventHandle[i] );
          if ( vrbflg ) fprintf( stdout, "-" );
        }
        else
        {
          /* add the event to the linked list */
          prevEvent = prevEvent->next = eventHandle[i];
          ++numTriggers[0];
          if ( vrbflg ) fprintf( stdout, "+" );
        }
      }
      else if ( errorParams.test == psi0_and_psi3 )
      {
        if ( (prevEvent->psi0 == eventHandle[i]->psi0)  &&
            (prevEvent->psi3 == eventHandle[i]->psi3) )
        {
          /* discard the event as it is a duplicate */
          LALFree( eventHandle[i] );
          if ( vrbflg ) fprintf( stdout, "-" );
        }
        else
        {
          /* add the event to the linked list */
          prevEvent = prevEvent->next = eventHandle[i];
          ++numTriggers[0];
          if ( vrbflg ) fprintf( stdout, "+" );
        }
      }
      else
      {
        fprintf( stderr, "error: unknown parameter test\n" );
        exit( 1 );
      }
    }

    /* if the list is non-emnpty, make sure it is terminated */
    if ( prevEvent ) prevEvent->next = NULL;

    if ( vrbflg ) 
    {
      fprintf( stdout, " done\n" );
      fprintf( stdout, "found %d sngl_inspiral rows for bank %s\n", 
          numTriggers[0], trigBankFile );
    }

    LALFree( eventHandle );

    /* skip the rest of the inca code and write out the bank */
    goto cleanexit;
  }


  /*
   *
   * for the case of BCV triggers, discard the triggers that
   * have alphaF greater than the alphaFcut specified
   *
   */

  if (errorParams.test == psi0_and_psi3) 
  {
    for ( j = 0; j < numIFO; ++j )
    {
      if ( vrbflg ) fprintf( stdout, 
          "Discarding triggers with alphaF > %f from ifo %d\n", alphaFcut, j );
      LAL_CALL( LALalphaFCutSingleInspiral( &status, &(inspiralEventList[j]),
        alphaFcut,0.0), &status );
    }  
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
    if ( vrbflg )
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

    ta = XLALGPSToINT8NS( &(currentTrigger[0]->end_time) );

    LAL_CALL( LALINT8NanoSecIsPlayground( &status, &isPlay, &ta ), &status );

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

    if( singleIfo )
    {
      if ( ( usePlayground && isPlay ) || ( ! usePlayground && ! isPlay) 
        || (allData) )
      {
        /* record the triggers */
        for ( j = 0; j < numIFO; ++j )
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
        }
      }
    }
    else
    {
      /* spin ifo b until the current trigger is within the coinicdence */
      /* window of the current ifo a trigger                            */
      while ( currentTrigger[1] )
      {
        tb = XLALGPSToINT8NS( &(currentTrigger[1]->end_time) );

        if ( tb > ta - errorParams.dt )
        {
          /* we have reached the time coinicidence window */
          break;
        }

        currentTrigger[1] = currentTrigger[1]->next;
      }

      /* if we are playground only and the trigger is in playground or    */
      /* we are not using playground and the trigger is not in the        */
      /* playground or we have a non-zero time-slide...                   */
      if ( ( usePlayground && isPlay ) || ( ! usePlayground && ! isPlay) ||
          (allData) || (slideDataNS) )
      {

        /* determine whether we should expect to see a trigger in ifo b  */
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
              ts = XLALGPSToINT8NS( &(currentEffRange[j]->start_time) );
              te = XLALGPSToINT8NS( &(currentEffRange[j]->end_time) );

              if ( (ts <= ta) && (ta < te) )
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
            upper_limit =( ( d_range[1] / d_range[0] ) * currentTrigger[0]->snr 
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
            fprintf( stdout, "SNR threshold in IFO B is %f\n", ifob_snrthresh);
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
          if ( vrbflg && currentTrigger[1] ) fprintf( stdout, 
              "  start loop over IFO B trigger at %d + %10.10f\n",
              currentTrigger[1]->end_time.gpsSeconds, 
              ((REAL4)currentTrigger[1]->end_time.gpsNanoSeconds * 1e-9) );

          /* look for coincident events in B within the time window */
          currentEvent = currentTrigger[1];

          while ( currentTrigger[1] )
          {
            tb = XLALGPSToINT8NS( &(currentTrigger[1]->end_time) );

            if (tb > ta + errorParams.dt )
            {
              /* we are outside the time coincidence so move to next event */
              if ( vrbflg ) 
                fprintf( stdout, "outside the time coincidence window\n" );
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
              "kept IFO A trigger although no coincident trigger in IFO B\n");
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
      }

    } 

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
      if ( userTag && ifoTag && outCompress )
      {
        snprintf( fileName, FILENAME_MAX, "%s-INCA_%s_%s-%d-%d.xml.gz",
            ifoName[j], ifoTag, userTag, startCoincidence,
            endCoincidence - startCoincidence );
      }
      else if ( userTag && !ifoTag && outCompress )
      {
        snprintf( fileName, FILENAME_MAX, "%s-INCA_%s-%d-%d.xml.gz",
            ifoName[j], userTag, startCoincidence,
            endCoincidence - startCoincidence );
      }
      else if ( !userTag && ifoTag && outCompress )
      {
        snprintf( fileName, FILENAME_MAX, "%s-INCA_%s-%d-%d.xml.gz",
            ifoName[j], ifoTag, startCoincidence,
            endCoincidence - startCoincidence );
      }
      else if ( !userTag && !ifoTag && outCompress )
      {
        snprintf( fileName, FILENAME_MAX, "%s-INCA-%d-%d.xml.gz",
            ifoName[j], startCoincidence, endCoincidence - startCoincidence );
      }
      else if ( userTag && ifoTag && !outCompress )
      {
        snprintf( fileName, FILENAME_MAX, "%s-INCA_%s_%s-%d-%d.xml", 
            ifoName[j], ifoTag, userTag, startCoincidence, 
            endCoincidence - startCoincidence );
      }
      else if ( userTag && !ifoTag && !outCompress )
      {
        snprintf( fileName, FILENAME_MAX, "%s-INCA_%s-%d-%d.xml", 
            ifoName[j], userTag, startCoincidence, 
            endCoincidence - startCoincidence );
      }
      else if ( !userTag && ifoTag && !outCompress )
      {
        snprintf( fileName, FILENAME_MAX, "%s-INCA_%s-%d-%d.xml", 
            ifoName[j], ifoTag, startCoincidence, 
            endCoincidence - startCoincidence );
      }
      else
      {
        snprintf( fileName, FILENAME_MAX, "%s-INCA-%d-%d.xml", ifoName[j],
            startCoincidence, endCoincidence - startCoincidence );
      }

      xmlFileName = fileName;
    }

    searchsumm.searchSummaryTable->nevents = numTriggers[j];

    memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
    LAL_CALL( LALOpenLIGOLwXMLFile( &status , &xmlStream, xmlFileName ), 
        &status );

    /* write process table */
    if ( trigBankFile || singleIfo )
    {
      snprintf( proctable.processTable->ifos, LIGOMETA_IFOS_MAX, "%s", 
          ifoName[0] );
    }
    else
    {
      snprintf( proctable.processTable->ifos, LIGOMETA_IFOS_MAX, "%s%s", 
          ifoName[0], ifoName[1] );
    }
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
    LAL_CALL( LALBeginLIGOLwXMLTable( &status ,&xmlStream, 
          search_summvars_table), &status );
    searchSummvarsTable.searchSummvarsTable = inputFiles;
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, searchSummvarsTable,
          search_summvars_table), &status );
    LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream), &status );

    /* write the summ_value table for ifoName[j] */
    if ( inspEffRange[j] )
    {
      SummValueTable tmpSummValueTable;
      if ( singleSummValue )
      {
        memcpy( &tmpSummValueTable, inspEffRange[j], sizeof(SummValueTable) );
        tmpSummValueTable.next = NULL;
        summValueTable.summValueTable = &tmpSummValueTable;
      }
      else
      {
        summValueTable.summValueTable = inspEffRange[j];
      }

      LAL_CALL( LALBeginLIGOLwXMLTable( &status ,&xmlStream, 
            summ_value_table), &status );
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
  if ( ifoTag ) free( ifoTag );

  if ( vrbflg ) fprintf( stdout, "done\n" );

  LALCheckMemoryLeaks();

  /* print a success message to stdout for parsing by exitcode */
  fprintf( stdout, "%s: EXITCODE0\n", argv[0] );
  exit( 0 );
}
