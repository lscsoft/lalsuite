/*
*  Copyright (C) 2007 Duncan Brown, Lisa M. Goggin
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
 * File Name: rinca.c
 *
 * Author: Goggin, L. M. based on thinca.c by Fairhurst, S.
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
#include <lal/TimeDelay.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRingdownRead.h>
#include <lal/LIGOMetadataRingdownUtils.h>
#include <lalapps.h>
#include <processtable.h>
#include <lal/SegmentsIO.h>

#include <LALAppsVCSInfo.h>

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "rinca"

#define INCA_EARG   1
#define INCA_EROW   2
#define INCA_EFILE  3

#define INCA_MSGEARG   "Error parsing arguments"
#define INCA_MSGROW    "Error reading row from XML table"
#define INCA_MSGEFILE  "Could not open file"

#define MAXIFO 4

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
  this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
calloc( 1, sizeof(ProcessParamsTable) ); \
snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
    PROGRAM_NAME ); \
snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
    long_options[option_index].name ); \
snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

int haveTrig[LAL_NUM_IFO];
int checkTimes = 0;
int multiIfoCoinc = 0;
int distCut = 0;
int h1h2Consistency = 0;
int doVeto = 0;
int completeCoincs = 0;
extern int vrbflg;

/*
 * 
 * USAGE
 *
 */
static void print_usage(char *program)
{
  fprintf(stderr,   "Usage:  %s [options] [LIGOLW XML input files]\n" ,  program                );
  fprintf(stderr,   "The following options are recognized.  Options not surrounded in [] are\n" );
  fprintf(stderr,   "required.\n"                                                               );

  fprintf(stderr,     "  [--help]                      display this message\n"                                 );
  fprintf(stderr,     "  [--verbose]                   print progress information\n"                           );
  fprintf(stderr,     "  [--version]                   print version information and exit\n"                   );
  fprintf(stderr,     "  [--debug-level]   level       set the LAL debug level to LEVEL\n"                     );
  fprintf(stderr,     "  [--user-tag]      usertag     set the process_params usertag\n"                       );
  fprintf(stderr,     "  [--ifo-tag]       ifotag      set the ifo-tag - for file naming\n"                    );
  fprintf(stderr,     "  [--comment]       string      set the process table comment to STRING\n"              );
  fprintf(stderr,     "  [--write-compress]            write a compressed xml file\n"                          );
  fprintf(stderr,     "\n"                                                                                     );
  fprintf(stderr,     "   --gps-start-time    start_time    GPS second of data start time\n"                   );
  fprintf(stderr,     "   --gps-end-time      end_time      GPS second of data end time\n"                     );
  fprintf(stderr,     "  [--check-times]                    Check that all times were analyzed\n"              );
  fprintf(stderr,     "  [--multi-ifo-coinc]                Look for triple/quadruple ifo coincidence\n"       );
  fprintf(stderr,     "  [--maximization-interval]  max_dt  set length of maximization interval in ms\n"       );
  fprintf(stderr,     "\n"                                                                                     );
  fprintf(stderr,     "  [--h1-slide]      h1_slide    Slide H1 data by multiples of h1_slide\n"               );
  fprintf(stderr,     "  [--h2-slide]      h2_slide    Slide H2 data by multiples of h2_slide\n"               );
  fprintf(stderr,     "  [--l1-slide]      l1_slide    Slide L1 data by multiples of l1_slide\n"               );
  fprintf(stderr,     "  [--v1-slide]      v1_slide    Slide V1 data by multiples of v1_slide\n"               );
  fprintf(stderr,     "  [--num-slides]    num_slides  The number of time slides to perform\n"                 );
  fprintf(stderr,     "\n"                                                                                     );
  fprintf(stderr,     "  [--h1-triggers]               input triggers from H1\n"                               );
  fprintf(stderr,     "  [--h2-triggers]               input triggers from H2\n"                               );
  fprintf(stderr,     "  [--l1-triggers]               input triggers from L1\n"                               );
  fprintf(stderr,     "  [--v1-triggers]               input triggers from L1\n"                               );
  fprintf(stderr,     "\n"                                                                                     );
  fprintf(stderr,     "   --parameter-test     test    set parameters with which to test coincidence:\n"       );
  fprintf(stderr,     "                                (f_and_Q, ds_sq, ds_sq_fQt)\n"                          );
  fprintf(stderr,     "\n"                                                                                     );
  fprintf(stderr,     "  [--h1-time-accuracy]  h1_dt   specify the timing accuracy of H1 in ms\n"              );
  fprintf(stderr,     "  [--h2-time-accuracy]  h2_dt   specify the timing accuracy of H2 in ms\n"              );
  fprintf(stderr,     "  [--l1-time-accuracy]  l1_dt   specify the timing accuracy of L1 in ms\n"              );
  fprintf(stderr,     "  [--v1-time-accuracy]  v1_dt   specify the timing accuracy of V1 in ms\n"              );
  fprintf(stderr,     "\n"                                                                                     );
  fprintf(stderr,     "  [--h1-freq-accuracy]  h1_df   specify the freq accuracy of H1\n"                      );
  fprintf(stderr,     "  [--h2-freq-accuracy]  h2_df   specify the freq accuracy of H2\n"                      );
  fprintf(stderr,     "  [--l1-freq-accuracy]  l1_df   specify the freq accuracy of L1\n"                      );
  fprintf(stderr,     "  [--v1-freq-accuracy]  v1_df   specify the freq accuracy of V1\n"                      );
  fprintf(stderr,     "\n"                                                                                     );
  fprintf(stderr,     "  [--h1-quality-accuracy]  h1_dq  specify the quality accuracy of H1\n"                 );
  fprintf(stderr,     "  [--h2-quality-accuracy]  h2_dq  specify the quality accuracy of H2\n"                 );
  fprintf(stderr,     "  [--l1-quality-accuracy]  l1_dq  specify the quality accuracy of L1\n"                 );
  fprintf(stderr,     "  [--v1-quality-accuracy]  v1_dq  specify the quality accuracy of V1\n"                 );
  fprintf(stderr,     "\n"                                                                                     );
  fprintf(stderr,     "  [--h1-distance-accuracy]  h1_ddeff   specify the effective distance accuracy of H1\n" );
  fprintf(stderr,     "  [--h2-distance-accuracy]  h2_ddeff   specify the effective distance accuracy of H2\n" );
  fprintf(stderr,     "  [--l1-distance-accuracy]  l1_ddeff   specify the effective distance accuracy of L1\n" );
  fprintf(stderr,     "  [--v1-distance-accuracy]  v1_ddeff   specify the effective distance accuracy of V1\n" );
  fprintf(stderr,     "\n"                                                                                     );
  fprintf(stderr,     "  [--h1-ds_sq-accuracy]  ds_sq     specify the ds squared accuracy\n"                   );
  fprintf(stderr,     "  [--h2-ds_sq-accuracy]  ds_sq     specify the ds squared accuracy\n"                   );
  fprintf(stderr,     "  [--l1-ds_sq-accuracy]  ds_sq     specify the ds squared accuracy\n"                   );
  fprintf(stderr,     "  [--v1-ds_sq-accuracy]  ds_sq     specify the ds squared accuracy\n"                   );
  fprintf(stderr,     "\n"                                                                                     );
  fprintf(stderr,     "   --data-type      data_type    specify the data type, must be one of\n"               );
  fprintf(stderr,     "                                 (playground_only|exclude_play|all_data)\n"             );
  fprintf(stderr,     "\n"                                                                                     );
  fprintf(stderr,     "  [--h1-h2-consistency]            perform H1-H2 consistency cut\n"                     );
  fprintf(stderr,     "  [--h1-snr-cut]         snr       reject h1 triggers below this snr when\n"            );
  fprintf(stderr,     "                                   there are no h2 triggers present\n"                  );
  fprintf(stderr,     "                                   needed when --h1-h2-consistency is given\n"          );
  fprintf(stderr,     "\n"                                                                                     );
  fprintf(stderr,     "  [--complete-coincs]                write out triggers from all non-vetoed ifos\n"     );
  fprintf(stderr,     "  [--do-veto]         do_veto        veto cetain segments\n"                            );
  fprintf(stderr,     "  [--h1-veto-file]    h1_veto_file   specify H1 triggers to be vetoed\n"                );
  fprintf(stderr,     "  [--h2-veto-file]    h2_veto_file   specify H2 triggers to be vetoed\n"                );
  fprintf(stderr,     "  [--l1-veto-file]    l1_veto_file   specify L1 triggers to be vetoed\n"                );
  fprintf(stderr,     "  [--v1-veto-file]    v1_veto_file   specify V1 triggers to be vetoed\n"                );
  fprintf(stderr,     "\n"                                                                                     );
  fprintf(stderr,     "  [LIGOLW XML input files]           list of the input trigger files.\n"                );
  fprintf(stderr,     "\n"                                                                                     );
}


int main( int argc, char *argv[] )
{
  static LALStatus      status;

  LALPlaygroundDataMask dataType = unspecified_data_type;
  INT4  startCoincidence = -1;
  LIGOTimeGPS startCoinc = {0,0};
  INT4  endCoincidence = -1;
  LIGOTimeGPS endCoinc = {0,0};

  REAL8        slideStep[LAL_NUM_IFO];
  LIGOTimeGPS  slideTimes[LAL_NUM_IFO];
  LIGOTimeGPS  slideReset[LAL_NUM_IFO];
  INT4         numSlides = 0;
  INT4         slideNum  = 0;

  CHAR  ifoName[MAXIFO][LIGOMETA_IFO_MAX];
  CHAR  ifos[LIGOMETA_IFOS_MAX];
  CHAR  ifoA[LIGOMETA_IFO_MAX];
  CHAR  ifoB[LIGOMETA_IFO_MAX];
  CHAR  ifo[LIGOMETA_IFO_MAX];

  CHAR  comment[LIGOMETA_COMMENT_MAX];
  CHAR *userTag = NULL;
  CHAR *ifoTag = NULL;

  CHAR  fileName[FILENAME_MAX];
  CHAR  fileSlide[FILENAME_MAX];
  CHAR *vetoFileName[LAL_NUM_IFO] = {NULL, NULL, NULL, NULL, NULL, NULL};

  LALSegList vetoSegs[LAL_NUM_IFO];

  UINT4   numIFO = 0;
  UINT4  numTrigIFO = 0;
  UINT4  numTriggers = 0;
  UINT4  numCoinc = 0;
  UINT4  numDoubles = 0;
  UINT4  numTriples = 0;
  UINT4  numTrigs[LAL_NUM_IFO];
  UINT4  N = 0;
  INT4 outCompress = 0;
  UINT4  slideH1H2Together = 0;

  LALDetector          aDet;
  LALDetector          bDet;

  SnglRingdownTable    *ringdownEventList = NULL;
  SnglRingdownTable    *thisRingdownTrigger = NULL;
  SnglRingdownTable    *snglOutput = NULL;
  SnglRingdownTable    *completionSngls = NULL;
  CoincRingdownTable   *coincRingdownList = NULL;
  CoincRingdownTable   *thisCoinc = NULL;

  EventIDColumn        *eventId;

  RingdownAccuracyList  accuracyParams;

  SearchSummvarsTable  *inputFiles = NULL;
  SearchSummvarsTable  *thisInputFile = NULL;

  SearchSummaryTable   *searchSummList = NULL;
  SearchSummaryTable   *thisSearchSumm = NULL;

  MetadataTable         proctable;
  MetadataTable         processParamsTable;
  MetadataTable         searchsumm;
  MetadataTable         searchSummvarsTable;
  MetadataTable         ringdownTable;
  ProcessParamsTable   *this_proc_param = NULL;
  LIGOLwXMLStream       xmlStream;

  InterferometerNumber  ifoNumber = LAL_UNKNOWN_IFO;
  InterferometerNumber  ifoTwo    = LAL_UNKNOWN_IFO;
  INT4                  i;
  REAL4                 maximizationInterval = 0;

  REAL4                 h1snrCut = 0;

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
    {"h1-triggers",         no_argument,   &(haveTrig[LAL_IFO_H1]),   1 },
    {"h2-triggers",         no_argument,   &(haveTrig[LAL_IFO_H2]),   1 },
    {"l1-triggers",         no_argument,   &(haveTrig[LAL_IFO_L1]),   1 },
    {"v1-triggers",         no_argument,   &(haveTrig[LAL_IFO_V1]),   1 },
    {"check-times",         no_argument,   &checkTimes,               1 },
    {"multi-ifo-coinc",     no_argument,   &multiIfoCoinc,            1 },
    {"h1-h2-distance-cut",  no_argument,   &distCut,                  1 },
    {"h1-h2-consistency",   no_argument,   &h1h2Consistency,          1 },
    {"do-veto",             no_argument,   &doVeto,                   1 },
    {"complete-coincs",     no_argument,   &completeCoincs,           1 },
    {"write-compress",      no_argument,   &outCompress,              1 },
    {"h1-slide",            required_argument, 0,                    'c'},
    {"h2-slide",            required_argument, 0,                    'd'},
    {"l1-slide",            required_argument, 0,                    'e'},
    {"v1-slide",            required_argument, 0,                    'f'},
    {"num-slides",          required_argument, 0,                    'T'},
    {"h1-time-accuracy",    required_argument, 0,                    'B'},
    {"h2-time-accuracy",    required_argument, 0,                    'C'},
    {"l1-time-accuracy",    required_argument, 0,                    'D'},
    {"v1-time-accuracy",    required_argument, 0,                    'A'},
    {"h1-freq-accuracy",    required_argument, 0,                    'H'},
    {"h2-freq-accuracy",    required_argument, 0,                    'I'},
    {"l1-freq-accuracy",    required_argument, 0,                    'J'},
    {"v1-freq-accuracy",    required_argument, 0,                    'K'},
    {"h1-quality-accuracy", required_argument, 0,                    'N'},
    {"h2-quality-accuracy", required_argument, 0,                    'O'},
    {"l1-quality-accuracy", required_argument, 0,                    'P'},
    {"v1-quality-accuracy", required_argument, 0,                    'Q'},
    {"h1-dist-accuracy",    required_argument, 0,                    'n'},
    {"h2-dist-accuracy",    required_argument, 0,                    'o'},
    {"l1-dist-accuracy",    required_argument, 0,                    'p'},
    {"v1-dist-accuracy",    required_argument, 0,                    'q'},
    {"h1-ds_sq-accuracy",   required_argument, 0,                    'E'},
    {"h2-ds_sq-accuracy",   required_argument, 0,                    'F'},
    {"l1-ds_sq-accuracy",   required_argument, 0,                    'G'},
    {"v1-ds_sq-accuracy",   required_argument, 0,                    'L'},
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
    {"high-mass",           required_argument, 0,                    '&'},
    {"h1-snr-cut",          required_argument, 0,                    '*'},
    {"h1-veto-file",        required_argument, 0,                    '('},
    {"h2-veto-file",        required_argument, 0,                    ')'},
    {"l1-veto-file",        required_argument, 0,                    '}'},
    {"v1-veto-file",        required_argument, 0,                    '{'},
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
  XLALGPSTimeNow(&(proctable.processTable->start_time));

  XLALPopulateProcessTable(proctable.processTable, PROGRAM_NAME,
      LALAPPS_VCS_IDENT_ID, LALAPPS_VCS_IDENT_STATUS, LALAPPS_VCS_IDENT_DATE, 0);

  this_proc_param = processParamsTable.processParamsTable = 
    (ProcessParamsTable *) calloc( 1, sizeof(ProcessParamsTable) );
  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );

  /* create the search summary and zero out the summvars table */
  searchsumm.searchSummaryTable = (SearchSummaryTable *)
    calloc( 1, sizeof(SearchSummaryTable) );

  memset( &accuracyParams, 0, sizeof(RingdownAccuracyList) );
  memset( &aDet, 0, sizeof(LALDetector) );

  /* set the time slide data to zero */
  memset( &slideStep, 0, LAL_NUM_IFO * sizeof(REAL8) );
  memset( &slideTimes, 0, LAL_NUM_IFO * sizeof(LIGOTimeGPS) );
  memset( &slideReset, 0, LAL_NUM_IFO * sizeof(LIGOTimeGPS) );
  memset( &haveTrig, 0, LAL_NUM_IFO * sizeof(int) );

  /* initialize array first */
  for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO ; ++ifoNumber )
    {      
      numTrigs[ifoNumber]=0;
    }

  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    long int gpstime;
    size_t optarg_len;

    c = getopt_long_only( argc, argv, 
        "B:C:D:E:F:G:H:I:J:K:N:O:P:T:V:Z:"
        "a:c:d:e:h:i:k:n:o:p:s:t:x:z:"
        "@:&:(:):}", 
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
        if ( ! strcmp( "f_and_Q", optarg ) )
          {
            accuracyParams.test = LALRINGDOWN_F_AND_Q;
          }
        else if ( ! strcmp( "ds_sq", optarg ) )
        {
          accuracyParams.test = LALRINGDOWN_DS_SQ;
        }
        else if ( ! strcmp( "ds_sq_fQt", optarg ) )
        {
          accuracyParams.test = LALRINGDOWN_DS_SQ_FQT;
        }
        else
          {
            fprintf( stderr, "invalid argument to --%s:\n"
                "unknown test specified: "
                "%s (must be f_and_Q, ds_sq, or ds_sq_fQt)\n",
                long_options[option_index].name, optarg );
            exit( 1 );
            }
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;
        
      case 'B':
        /* time accuracy H1, argument is in milliseconds */
        accuracyParams.ifoAccuracy[LAL_IFO_H1].dt = atof(optarg) * LAL_INT8_C(1000000);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'C':
        /* time accuracy H2, argument is in milliseconds */
        accuracyParams.ifoAccuracy[LAL_IFO_H2].dt = atof(optarg) * LAL_INT8_C(1000000);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'D':
        /* time accuracy L1, argument is in milliseconds */
        accuracyParams.ifoAccuracy[LAL_IFO_L1].dt = atof(optarg) * LAL_INT8_C(1000000);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
    
      case 'A':
        /* time accuracy V1, argument is in milliseconds */
        accuracyParams.ifoAccuracy[LAL_IFO_V1].dt = atof(optarg) * LAL_INT8_C(1000000);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
    
      case 'E':
        /* ds^2 accuracy H1*/
        accuracyParams.ifoAccuracy[LAL_IFO_H1].ds_sq = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'F':
        /* ds^2 accuracy H2*/
        accuracyParams.ifoAccuracy[LAL_IFO_H2].ds_sq = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

     case 'G':
        /* ds^2 accuracy L1*/
        accuracyParams.ifoAccuracy[LAL_IFO_L1].ds_sq = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

     case 'L':
        /* ds^2 accuracy V1*/
        accuracyParams.ifoAccuracy[LAL_IFO_V1].ds_sq = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'H':
        /* frequency accuracy H1, argument is in Hz */
        accuracyParams.ifoAccuracy[LAL_IFO_H1].df = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'I':
        /* frequency accuracy H2, argument is in Hz */
        accuracyParams.ifoAccuracy[LAL_IFO_H2].df = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'J':
        /* frequency accuracy L1, argument is in Hz */
        accuracyParams.ifoAccuracy[LAL_IFO_L1].df = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'K':
        /* frequency accuracy V1, argument is in Hz */
        accuracyParams.ifoAccuracy[LAL_IFO_V1].df = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'N':
        /* quality factor accuracy H1 */
        accuracyParams.ifoAccuracy[LAL_IFO_H1].dQ = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'O':
        /* quality factor accuracy H2 */
        accuracyParams.ifoAccuracy[LAL_IFO_H2].dQ = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'P':
        /* quality factor accuracy L1 */
        accuracyParams.ifoAccuracy[LAL_IFO_L1].dQ = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'Q':
        /* quality factor accuracy V1 */
        accuracyParams.ifoAccuracy[LAL_IFO_V1].dQ = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'c':
        /* slide time for H1 */
        slideStep[LAL_IFO_H1] = atof( optarg );
        if ( slideStep[LAL_IFO_H1] < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "The slideStep must be positive\n"
              "(%f specified)\n",
              long_options[option_index].name, slideStep[LAL_IFO_H1] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "double", "%lf", slideStep[LAL_IFO_H1] );
        break;
        
      case 'd':
        /* slide time for H2 */
        slideStep[LAL_IFO_H2] = atof( optarg );
        if ( slideStep[LAL_IFO_H2] < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "The slideStep must be positive\n"
              "(%f specified)\n",
              long_options[option_index].name, slideStep[LAL_IFO_H2] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "double", "%lf", slideStep[LAL_IFO_H2] );
        break;
        
      case 'e':
        /* slide time for L1 */
        slideStep[LAL_IFO_L1] = atof( optarg );
        if ( slideStep[LAL_IFO_L1] < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "The slideStep must be positive\n"
              "(%f specified)\n",
              long_options[option_index].name, slideStep[LAL_IFO_L1] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "double", "%lf", slideStep[LAL_IFO_L1] );
        break;
        
      case 'f':
        /* slide time for V1 */
        slideStep[LAL_IFO_V1] = atof( optarg );
        if ( slideStep[LAL_IFO_V1] < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "The slideStep must be positive\n"
              "(%f specified)\n",
              long_options[option_index].name, slideStep[LAL_IFO_V1] );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "double", "%lf", slideStep[LAL_IFO_V1] );
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
        ADD_PROCESS_PARAM( "int", "%" LAL_INT4_FORMAT, numSlides );
        break;

      case 'n':
        /* effective distance for H1, argument is in Mpc */
        accuracyParams.ifoAccuracy[LAL_IFO_H1].ddeff = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
        
      case 'o':
        /* effective distance H2, argument is in Mpc  */
        accuracyParams.ifoAccuracy[LAL_IFO_H2].ddeff = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;

      case 'p':
        /* effective distance L1, argument is in Mpc  */
        accuracyParams.ifoAccuracy[LAL_IFO_L1].ddeff = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
                                
      case 'q':
        /* effective distance V1, argument is in Mpc  */
        accuracyParams.ifoAccuracy[LAL_IFO_V1].ddeff = atof(optarg);
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
        startCoincidence = (INT4) gpstime;
        ADD_PROCESS_PARAM( "int", "%" LAL_INT4_FORMAT, startCoincidence );
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
        endCoincidence = (INT4) gpstime;
        ADD_PROCESS_PARAM( "int", "%" LAL_INT4_FORMAT, endCoincidence );
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
          snprintf( comment, LIGOMETA_COMMENT_MAX, "%s", optarg);
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
        snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", 
            PROGRAM_NAME );
        snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "-userTag" );
        snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
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
        fprintf( stdout, "RINgdown Coincidence Analysis\n" 
            "Lisa Goggin based on thinca.c by Steve Fairhurst\n");
        XLALOutputVersionString(stderr, 0);
        exit( 0 );
        break;

      case '?':
        print_usage(argv[0]);
        exit( 1 );
        break;

      case '@':
        /* set the maximization window */
        maximizationInterval = atof( optarg );
        if ( maximizationInterval < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "maximization interval must be positive:\n "
              "(%" LAL_REAL4_FORMAT " ms specified)\n",
              long_options[option_index].name, maximizationInterval );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%" LAL_REAL4_FORMAT,  maximizationInterval );
        break;

      case '*':
        /* snr cut */
        h1snrCut = atof(optarg);
        ADD_PROCESS_PARAM( "float", "%s", optarg );
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
        vetoFileName[LAL_IFO_V1] = (CHAR *) calloc( optarg_len, sizeof(CHAR));
        memcpy( vetoFileName[LAL_IFO_V1], optarg, optarg_len );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
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
  if ( accuracyParams.test == LALRINGDOWN_UNKNOWN_TEST )
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
      snprintf( ifoName[numIFO], LIGOMETA_IFO_MAX, "%s", ifoList[ifoNumber] );
      numIFO++;

      /* store the argument in the process_params table */
      this_proc_param = this_proc_param->next = (ProcessParamsTable *)
        calloc( 1, sizeof(ProcessParamsTable) );
      snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
          "%s", PROGRAM_NAME );
      snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", 
          ifoArg[ifoNumber]);
      snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
      snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );

      /* check that a non-zero timing accuracy was specified */
      if ( (accuracyParams.test == LALRINGDOWN_F_AND_Q || accuracyParams.test == LALRINGDOWN_DS_SQ)
           && ! accuracyParams.ifoAccuracy[ifoNumber].dt )
      {
        fprintf( stderr, "Error: --dt must be specified for %s\n",
            ifoName[ifoNumber]);
        exit(1);
      }
    }
  }

if ( vrbflg)
    {
      fprintf( stderr, "number of ifos=%d\n",numIFO);
    }


  /* check that we have at least two IFOs specified, or can't do coincidence */
  if ( numIFO < 2 )
  {
    fprintf( stderr, "Must specify at least two IFOs to do coincidence\n"
        "%d specified\n", numIFO );
    exit ( 1 );
  }
 
  if ( numIFO > 2 && vrbflg )
  {
    if ( !multiIfoCoinc )
    {
      fprintf( stdout, 
          "Finding all double coincidences in %d IFO time.\n"
          "If you want triples please specify --multi-ifo-coinc.\n",
          numIFO);
    }
    else
    {
       fprintf( stdout, 
           "Finding all double/triple coincidences in %d IFO time.\n",
           numIFO);
    }
  }
  
  /* set ifos to be the alphabetical list of the ifos with triggers */
  if( numIFO == 2 )
  {
    snprintf( ifos, LIGOMETA_IFOS_MAX, "%s%s", ifoName[0], ifoName[1] );
  }
  else if ( numIFO == 3 )
  {
    snprintf( ifos, LIGOMETA_IFOS_MAX, "%s%s%s", ifoName[0], ifoName[1], 
        ifoName[2] );
  }

  /* if numSlides is set, check that the slide times are different for
   * different ifos (for which we have triggers, H1 and H2 are allowed
   * to slide together */
  if( numSlides )
  {
    for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
    {
      XLALReturnIFO( ifoA, ifoNumber );
     
      if( vrbflg && haveTrig[ifoNumber] ) fprintf( stdout, 
          "Performing a slide of multiples of %f seconds on %s\n", 
          slideStep[ifoNumber], ifoA);
          
      for( ifoTwo = ifoNumber + 1; ifoTwo < LAL_NUM_IFO; ifoTwo++ )
      {
        if ( haveTrig[ifoTwo] && haveTrig[ifoNumber] &&
            slideStep[ifoTwo] == slideStep[ifoNumber] )
        {
          XLALReturnIFO( ifoB, ifoTwo );

          if ( ( ! strcmp(ifoA,"H1") && ! strcmp(ifoB,"H2") ) ||
               ( ! strcmp(ifoA,"H2") && ! strcmp(ifoB,"H1") ) )
          {
            if( vrbflg ) fprintf( stdout,
              "The time slide specified for ifo %s is %f\n"
              "The time slide specified for ifo %s is also %f\n",
              ifoA, slideStep[ifoNumber], ifoB, slideStep[ifoTwo]);
            slideH1H2Together = 1;
          }
          else
          {
            fprintf( stderr,
              "The time slide specified for ifo %s is %f\n"
              "The time slide specified for ifo %s is also %f\n"
              "Must specify unique time slides for all instruments\n",
              ifoA, slideStep[ifoNumber], ifoB, slideStep[ifoTwo]);

            exit( 1 );
          }
        }
      }
    }
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

  if ( outCompress )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--write-compress" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }


  /* store the check-times in the process_params table */
  if ( checkTimes )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--check-times" );
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
  }

  /* store the h1h2 consistency option */ 
  if ( h1h2Consistency )
  {
    if ( !( h1snrCut > 0 ) )
    {
      fprintf( stderr,"The --h1-snr-cut option must be specified and \n"
               "greater than zero when --h1-h2-consistency is given\n" );
      exit( 1 );
    }
     
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
   snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--h1-h2-consistency" );
   snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
   snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
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

  /* store the complete-coincs option */
  if (completeCoincs)
  {
    if (!multiIfoCoinc)
    {
      fprintf( stderr,
          "--multi-ifo-coinc must be specified when --complete-coincs is.\n" );
      exit( 1 );
    }

    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX,
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX,
        "--complete-coincs");
    snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, " " );
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
      INT4 numFileTriggers = 0;
      SnglRingdownTable   *ringdownFileList = NULL;
      SnglRingdownTable   *thisFileTrigger  = NULL;

      numFileTriggers = XLALReadRingdownTriggerFile( &ringdownFileList,
          &thisFileTrigger, &searchSummList, &inputFiles, argv[i] );
      if (numFileTriggers < 0)
      {
        fprintf(stderr, "Error reading triggers from file %s",
            argv[i]);
        exit( 1 );
      }
      
       /* maximize over a given interval */
      if ( maximizationInterval )
      {
        if (vrbflg)
        {
          fprintf( stdout, "Clustering triggers for over %" LAL_REAL4_FORMAT " ms window\n",
              maximizationInterval);
        }
        XLALMaxSnglRingdownOverIntervals( &ringdownFileList,
            (INT8) round(1.0e6 * maximizationInterval) );
        numFileTriggers = XLALCountSnglRingdown( ringdownFileList );
      }
      
      numTriggers += numFileTriggers;
      /* If there are any remaining triggers ... */
      if ( ringdownFileList )
      {
        /* add ringdowns to list */
        if ( thisRingdownTrigger )
        {
          thisRingdownTrigger->next = ringdownFileList;
        }
         else
        {
          ringdownEventList = thisRingdownTrigger = ringdownFileList;
        }
         for( ; thisRingdownTrigger->next;
             thisRingdownTrigger = thisRingdownTrigger->next);
      }
      
    }
  
  }
  else
  {
    fprintf( stderr, "Error: No trigger files specified.\n" );
    exit( 1 );
  }

  if ( vrbflg && !maximizationInterval ) 
    fprintf( stdout, "Read in a total of %d triggers.\n",  numTriggers );
  else if ( vrbflg && maximizationInterval )
    fprintf( stdout, "After maximization have a total of %d triggers.\n",  numTriggers );
   
   /*    we initialise the veto segment list needed either by the h1h2 
        consistency check or the veto option itself. */
   if ( h1h2Consistency || doVeto )
   {
     for ( ifoNumber = 0; ifoNumber< LAL_NUM_IFO; ifoNumber++)
     {
       /* to perform the veto, we  need the filename. */
       if ( vetoFileName[ifoNumber] && haveTrig[ifoNumber])
       {
         XLALSegListInit( &(vetoSegs[ifoNumber]) );
         LAL_CALL( LALSegListRead( &status, &(vetoSegs[ifoNumber]),
               vetoFileName[ifoNumber], NULL),&status);
         XLALSegListCoalesce( &(vetoSegs[ifoNumber]) );

         /* keep only the segments that lie within the data-segment part */
         XLALSegListKeep(  &(vetoSegs[ifoNumber]), &startCoinc, &endCoinc );

         /* if the veto option is set, we remove single ringdown triggers 
            inside the list provided but we need to loop over the different 
            ifo name. */
         if (doVeto)
         {
           XLALReturnIFO(ifo,ifoNumber);
           if ( vrbflg ) fprintf( stdout,
               "Applying veto segment (%s) list on ifo  %s \n ",
               vetoFileName[ifoNumber], ifo );
           ringdownEventList = XLALVetoSingleRingdown( ringdownEventList,
               &(vetoSegs[ifoNumber]), ifo );
           /* count the triggers  */
           numTriggers = XLALCountSnglRingdown( ringdownEventList );
           if ( vrbflg ) fprintf( stdout,
               " --> %d remaining triggers after veto segment list applied.\n",
               numTriggers );
         }
       }
     }
   }

  /* check that we have read in data for all the requested times
     in all the requested instruments */
  if ( checkTimes )
  {
    if ( vrbflg ) fprintf( stdout, 
        "Checking that we have data for all times from all IFOs\n");
    for ( ifoNumber = 0; ifoNumber < (int)numIFO; ++ifoNumber )
    {
      LAL_CALL( LALCheckOutTimeFromSearchSummary ( &status, searchSummList, 
            ifoName[ifoNumber], &startCoinc, &endCoinc ), &status);
    }
  }

  if ( ! ringdownEventList )
  {
    /* no triggers, so no coincidences can be found */
    if ( vrbflg ) fprintf( stdout,
        "No triggers read in so no coincidences can be found\n" );

    goto cleanexit;
  }

  /* time sort the triggers */
  if ( vrbflg ) fprintf( stdout, "Sorting triggers\n" );
  LAL_CALL( LALSortSnglRingdown( &status, &(ringdownEventList),
        LALCompareSnglRingdownByTime ), &status );

  /* keep only triggers within the requested interval */
  if ( vrbflg ) fprintf( stdout, 
      "Discarding triggers outside requested interval\n" );
  LAL_CALL( LALTimeCutSingleRingdown( &status, &ringdownEventList,
        &startCoinc, &endCoinc), &status );


  /* keep play/non-play/all triggers */
  if ( dataType == playground_only && vrbflg ) fprintf( stdout, 
      "Keeping only playground triggers\n" );
  else if ( dataType == exclude_play && vrbflg ) fprintf( stdout, 
      "Keeping only non-playground triggers\n" );
  else if ( dataType == all_data && vrbflg ) fprintf( stdout, 
      "Keeping all triggers\n" );
  LAL_CALL( LALPlayTestSingleRingdown( &status, &ringdownEventList,
        &dataType ), &status );

  /* scroll to the end of the linked list of triggers, counting triggers and
   * freeing event_id's that were read in */
  thisRingdownTrigger = ringdownEventList;
  for (numTriggers = 0 ; thisRingdownTrigger; ++numTriggers,
      thisRingdownTrigger = thisRingdownTrigger->next )
  {  
    while ( thisRingdownTrigger->event_id )
    {
      eventId = (thisRingdownTrigger)->event_id;
      (thisRingdownTrigger)->event_id = (thisRingdownTrigger)->event_id->next;
      LALFree( eventId );
    }
  }

  if ( vrbflg ) fprintf( stdout, 
      "%d remaining triggers after time and data type cut.\n", numTriggers );


  if ( ! ringdownEventList )
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
    LALIfoCountSingleRingdown(&status, &numTrigs[ifoNumber], 
        ringdownEventList, ifoNumber);
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

  /* perform the time slides */
  for( slideNum = -numSlides; slideNum <= numSlides; slideNum++ )
  {
    SnglRingdownTable    *slideOutput = NULL;
    INT4                  numCoincInSlide = 0;

    coincRingdownList = NULL;

    if ( numSlides )
    {
      for( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++ )
      {
        if( slideNum == -numSlides )
        {
          /* Initialize the slide-times with the initial value, 
          which is the negative extreme value */
          REAL8 tmpInt;
  	  INT4 tmpSlide, tmpSlideNS;
  	  tmpSlideNS = (INT4) (-1000000000*
  	    modf(numSlides * slideStep[ifoNumber], &tmpInt) );
  	  tmpSlide = (INT4) (-tmpInt);
          slideTimes[ifoNumber].gpsSeconds = tmpSlide;
          slideTimes[ifoNumber].gpsNanoSeconds = tmpSlideNS;
          slideReset[ifoNumber].gpsSeconds = (-tmpSlide);
          slideReset[ifoNumber].gpsNanoSeconds = (-tmpSlideNS);
        }
        else
        {
          /* Set the slide-times to the constant slide-step, 
             since this times refers to 'ringdownEventList', 
             which triggers are shifted in each slide by a constant amount.
             The reset-time however refers to the coincidence list, so it must
             be decreased for each slide. */
          REAL8 tmpInt;
          INT4 tmpSlide, tmpSlideNS;
  	  tmpSlideNS = (INT4) (-1000000000*
  	    modf(slideStep[ifoNumber], &tmpInt) );
  	  tmpSlide = (INT4) tmpInt;
 
  	  slideTimes[ifoNumber].gpsSeconds = tmpSlide;
  	  slideTimes[ifoNumber].gpsNanoSeconds = tmpSlideNS;
  	  slideReset[ifoNumber].gpsSeconds -= tmpSlide;
 	  slideReset[ifoNumber].gpsNanoSeconds -= tmpSlideNS;
        }
      }
    }

    if ( vrbflg ) fprintf(stdout, "Performing time slide %d\n", slideNum );

    /* slide the data */
    LAL_CALL( LALTimeSlideSingleRingdown( &status, ringdownEventList,
              &startCoinc, &endCoinc, slideTimes), &status) ;
    LAL_CALL( LALSortSnglRingdown( &status, &(ringdownEventList),
              LALCompareSnglRingdownByTime ), &status );

    for ( ifoNumber = 0; ifoNumber< LAL_NUM_IFO; ifoNumber++)
    {
      if ( vetoSegs[ifoNumber].initMagic != SEGMENTSH_INITMAGICVAL )
        /* no veto segs */
        continue;
      if ( XLALTimeSlideSegList( &vetoSegs[ifoNumber], &startCoinc, &endCoinc,
                                 &slideTimes[ifoNumber] ) < 0 )
        exit(1);
    }
 
    /* don't analyze zero-lag if numSlides>0 */
    if ( numSlides && !slideNum ) continue;

    /* look for coincidences */
    LAL_CALL( LALCreateTwoIFORingdownCoincList(&status, &coincRingdownList,
          ringdownEventList, &accuracyParams ), &status);
    
    /* count the zero-lag coincidences */ 
    if ( !numSlides)
    {
      if( coincRingdownList )
      {
        for (numCoinc = 0, thisCoinc = coincRingdownList;
               thisCoinc; ++numCoinc, thisCoinc = thisCoinc->next )
        {
        }
        if ( vrbflg ) fprintf( stdout,
            "%d coincident triggers found before coincidence cuts..\n",
            numCoinc);
      }
    }

    if ( multiIfoCoinc )
    {
      for( N = 3; N <= numIFO; N++)
      {
        LAL_CALL( LALCreateNIFORingdownCoincList( &status, &coincRingdownList,
               &accuracyParams, N ), &status );
      }
 
      LAL_CALL( LALRemoveRepeatedRingdownCoincs( &status, &coincRingdownList ),
           &status );
    }

    /* perform the h1h2-consistency check */
    if ( h1h2Consistency && haveTrig[LAL_IFO_H1] && haveTrig[LAL_IFO_H2] )
    {
      if(vrbflg) 
      {
        if (vetoFileName[LAL_IFO_H1] && vetoFileName[LAL_IFO_H2])
        {
          fprintf(stdout, 
              "Using h1-h2-consistency with veto segment list %s and %s\n", 
              vetoFileName[LAL_IFO_H1], vetoFileName[LAL_IFO_H2]);
        }
        else
        { 
          fprintf(stdout, 
              "Using h1-h2-consistency without veto segment list. NOT RECOMMENDED\n");
        }
      }
      LAL_CALL( LALRingdownH1H2Consistency(&status,  &coincRingdownList,
           h1snrCut, &vetoSegs[LAL_IFO_H1], &vetoSegs[LAL_IFO_H2]), &status);
      if ( vrbflg ) fprintf( stdout, 
          "%d remaining coincident triggers after h1-h2-consisteny .\n", 
          XLALCountCoincInspiral((CoincInspiralTable*)coincRingdownList));
   }


    /* no time-slide */
    if ( !slideNum )
    {
      /* count the coincs */
      if( coincRingdownList )
      {
        for (numCoinc = 0, thisCoinc = coincRingdownList;
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
        }
      }

      if ( vrbflg )
      {
        fprintf( stdout, "%d coincident triggers found.\n", numCoinc );
        fprintf( stdout, "%d double coincident triggers\n"
                         "%d triple coincident triggers\n",
                         numDoubles, numTriples );
      }
    }

    if ( coincRingdownList )
    {
      /* count the coincs, scroll to end of list */
      if( slideNum )
      {
        
        /* keep only the requested coincs */
        if( slideH1H2Together )
        {
          char slide_ifos[] = "H1H2";
          if ( vrbflg ) fprintf( stdout,
              "Throwing out slide coincs found only as H1H2 doubles.\n" );
          numCoincInSlide = XLALCoincRingdownIfosDiscard(
              &coincRingdownList, slide_ifos );
          if ( vrbflg ) fprintf( stdout,
              "Kept %d non-H1H2 coincs in slide.\n", numCoincInSlide );
        }

        numCoincInSlide = 0;
        thisCoinc = coincRingdownList;
  	while( thisCoinc )
  	{
  	  ++numCoincInSlide;
  	  thisCoinc = thisCoinc->next;
  	}
 	 
        if ( vrbflg ) fprintf( stdout,
            "%d coincident triggers found in slide.\n", numCoincInSlide );
        numCoinc += numCoincInSlide;
      }

      /* complete the coincs */
      if ( completeCoincs )
      {
        completionSngls = XLALCompleteCoincRingdown ( coincRingdownList, 
            haveTrig);

        /* do a veto on the new sngls */
        for ( ifoNumber = 0; ifoNumber< LAL_NUM_IFO; ifoNumber++)
        {
        
          if (doVeto && vetoFileName[ifoNumber] && haveTrig[ifoNumber])
          {
            XLALReturnIFO(ifo,ifoNumber);
            if ( vrbflg ) fprintf( stdout, 
                "Applying veto list %s on completion sngls for ifo %s \n",
                vetoFileName[ifoNumber], ifo );
            completionSngls = XLALVetoSingleRingdown( completionSngls,
               &(vetoSegs[ifoNumber]), ifo );
          }
        }
      }

      /* write out all coincs as singles with event IDs */
      LAL_CALL( LALExtractSnglRingdownFromCoinc( &status, &slideOutput,
            coincRingdownList, &startCoinc, slideNum), &status );

      if ( numSlides && coincRingdownList )
      {

        /* the output triggers should be slid back to original time */
        LAL_CALL( LALTimeSlideSingleRingdown( &status, slideOutput,
              &startCoinc, &endCoinc, slideReset), &status) ;
        LAL_CALL( LALSortSnglRingdown( &status, &(slideOutput),
              LALCompareSnglRingdownByTime ), &status );
      }

      while ( coincRingdownList )
      {
        thisCoinc = coincRingdownList;
        coincRingdownList = coincRingdownList->next;
        XLALFreeCoincRingdown( &thisCoinc );
      }

      while ( completionSngls )
      {
        SnglRingdownTable *thisSngl = NULL;
        thisSngl = completionSngls;
        completionSngls = completionSngls->next;
        XLALFreeSnglRingdown( &thisSngl );
      }

    }

    if ( snglOutput )
    {
       thisRingdownTrigger->next = slideOutput;
    }
    else
    {
      snglOutput = slideOutput;
    }

    if ( numSlides )
    {
      /* scroll to the end of the list */
      if ( slideOutput )
      {
        for( thisRingdownTrigger = slideOutput; thisRingdownTrigger->next;
             thisRingdownTrigger = thisRingdownTrigger->next);
      }
    }
  }


  /*
   *
   * write the output xml file
   *
   */


  /* since we don't yet write coinc ringdown tables, we must make a list of
   * sngl_ringdown tables with the eventId's appropriately poplulated */
  
   
cleanexit:

  searchsumm.searchSummaryTable->in_start_time = startCoinc;
  searchsumm.searchSummaryTable->in_end_time = endCoinc;
  searchsumm.searchSummaryTable->out_start_time = startCoinc;
  searchsumm.searchSummaryTable->out_end_time = endCoinc;
  searchsumm.searchSummaryTable->nnodes = 1;

  if ( vrbflg ) fprintf( stdout, "writing output file... " );

  if ( userTag && ifoTag && !outCompress)
  {
    snprintf( fileName, FILENAME_MAX, "%s-RINCA_%s_%s-%d-%d.xml", 
        ifos, ifoTag, userTag, startCoincidence, 
        endCoincidence - startCoincidence );
    snprintf( fileSlide, FILENAME_MAX, "%s-RINCA_SLIDE_%s_%s-%d-%d.xml", 
        ifos, ifoTag, userTag, startCoincidence, 
        endCoincidence - startCoincidence );
  }
  else if  ( !userTag && ifoTag && !outCompress ) 
  {
    snprintf( fileName, FILENAME_MAX, "%s-RINCA_%s-%d-%d.xml", ifos,
        ifoTag, startCoincidence, endCoincidence - startCoincidence );
    snprintf( fileSlide, FILENAME_MAX, "%s-RINCA_SLIDE_%s-%d-%d.xml", ifos,
        ifoTag, startCoincidence, endCoincidence - startCoincidence );
  }
  else if ( userTag && !ifoTag && !outCompress )
  {
    snprintf( fileName, FILENAME_MAX, "%s-RINCA_%s-%d-%d.xml", 
        ifos, userTag, startCoincidence, endCoincidence - startCoincidence );
    snprintf( fileSlide, FILENAME_MAX, "%s-RINCA_SLIDE_%s-%d-%d.xml", 
        ifos, userTag, startCoincidence, endCoincidence - startCoincidence );
  }
  else if ( userTag && ifoTag && outCompress)
  {
    snprintf( fileName, FILENAME_MAX, "%s-RINCA_%s_%s-%d-%d.xml.gz",   
        ifos, ifoTag, userTag, startCoincidence,
        endCoincidence - startCoincidence );
    snprintf( fileSlide, FILENAME_MAX, "%s-RINCA_SLIDE_%s_%s-%d-%d.xml.gz",   
        ifos, ifoTag, userTag, startCoincidence,
        endCoincidence - startCoincidence );
  }
  else if  ( !userTag && ifoTag && outCompress ) 
  {
    snprintf( fileName, FILENAME_MAX, "%s-RINCA_%s-%d-%d.xml.gz", ifos,
        ifoTag, startCoincidence, endCoincidence - startCoincidence );
    snprintf( fileSlide, FILENAME_MAX, "%s-RINCA_SLIDE_%s-%d-%d.xml.gz", ifos,
        ifoTag, startCoincidence, endCoincidence - startCoincidence );
  }
  else if ( userTag && !ifoTag && outCompress )
  {
    snprintf( fileName, FILENAME_MAX, "%s-RINCA_%s-%d-%d.xml.gz",
        ifos, userTag, startCoincidence, endCoincidence - startCoincidence );
    snprintf( fileSlide, FILENAME_MAX, "%s-RINCA_SLIDE_%s-%d-%d.xml.gz",
        ifos, userTag, startCoincidence, endCoincidence - startCoincidence );
  }
  else if ( !userTag && !ifoTag && outCompress )
  {
    snprintf( fileName, FILENAME_MAX, "%s-RINCA-%d-%d.xml.gz",
        ifos, startCoincidence, endCoincidence - startCoincidence );
    snprintf( fileSlide, FILENAME_MAX, "%s-RINCA_SLIDE-%d-%d.xml.gz",
        ifos, startCoincidence, endCoincidence - startCoincidence );
  }
  else
  {
    snprintf( fileName, FILENAME_MAX, "%s-RINCA-%d-%d.xml", ifos,
        startCoincidence, endCoincidence - startCoincidence );
    snprintf( fileSlide, FILENAME_MAX, "%s-RINCA_SLIDE-%d-%d.xml", ifos,
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

  snprintf( proctable.processTable->ifos, LIGOMETA_IFOS_MAX, "%s", ifos );

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
  snprintf( searchsumm.searchSummaryTable->ifos, LIGOMETA_IFOS_MAX, "%s", ifos );

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

  /* write the sngl_ringdown table */
  if( snglOutput )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status ,&xmlStream, 
          sngl_ringdown_table), &status );
    ringdownTable.snglRingdownTable = snglOutput;
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, ringdownTable,
          sngl_ringdown_table), &status );
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
    if ( vetoFileName[ifoNumber]  && haveTrig[ifoNumber])
    {
      free( vetoFileName[ifoNumber] );
    }

   if (vetoSegs[ifoNumber].initMagic == SEGMENTSH_INITMAGICVAL )
   {
        XLALSegListClear( &vetoSegs[ifoNumber] );
    }
  }

  /* free the snglRingdown */
  while ( ringdownEventList )
  {
    thisRingdownTrigger = ringdownEventList;
    ringdownEventList = ringdownEventList->next;
    LAL_CALL( LALFreeSnglRingdown( &status, &thisRingdownTrigger ), &status );
  }

  while ( snglOutput )
  {
    thisRingdownTrigger = snglOutput;
    snglOutput = snglOutput->next;
    LAL_CALL( LALFreeSnglRingdown( &status, &thisRingdownTrigger ), &status );
  }
 
  while ( coincRingdownList )
  {
    thisCoinc = coincRingdownList;
    coincRingdownList = coincRingdownList->next;
    LALFree( thisCoinc );
  }


  if ( userTag ) free( userTag );
  if ( ifoTag ) free( ifoTag );

  if ( vrbflg ) fprintf( stdout, "done\n" );

  LALCheckMemoryLeaks();

  exit( 0 );
}
