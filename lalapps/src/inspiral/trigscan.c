/*
*  Copyright (C) 2007 Alexander Dietz, Duncan Brown, Patrick Brady, Stephen Fairhurst, Sean Seader
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
 * File Name: trigbank.c
 *
 * Author: Brady, P. R., Brown, D. A., Fairhurst, S. and Pekowsky, L.
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
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LALTrigScanCluster.h>
#include <lalapps.h>
#include <processtable.h>
#include <LALAppsVCSInfo.h>

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "trigbank"

#define TRIGSCAN_EARG   1
#define TRIGSCAN_EROW   2
#define TRIGSCAN_EFILE  3

#define TRIGSCAN_MSGEARG   "Error parsing arguments"
#define TRIGSCAN_MSGROW    "Error reading row from XML table"
#define TRIGSCAN_MSGEFILE  "Could not open file"

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
  this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
calloc( 1, sizeof(ProcessParamsTable) ); \
snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
    PROGRAM_NAME ); \
snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
    long_options[option_index].name ); \
snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

extern int vrbflg;

/*
 *
 * USAGE
 *
 */

static void print_usage(char *program)
{
  fprintf(stderr,
      "Usage: %s [options] [LIGOLW XML input files]\n"\
      "The following options are recognized.  Options not surrounded in [] are\n" \
      "required.\n" \
      "  [--help]                      display this message\n"\
      "  [--verbose]                   print progress information\n"\
      "  [--version]                   print version information and exit\n"\
      "  [--user-tag]      usertag     set the process_params usertag\n"\
      "  [--comment]       string      set the process table comment\n"\
      "  [--write-compress]            write a compressed xml file\n"\
      "  --gps-start-time start_time   GPS second of data start time\n"\
      "  --gps-end-time   end_time     GPS second of data end time\n"\
      "  --ifo             ifo         set the ifo - for file naming\n"\
      "  --ts-cluster      MTHD        max over template and end time MTHD \n"\
      "                                (T0T3Tc|T0T3TcAS|Psi0Psi3Tc|Psi0Psi3TcAS) \n"\
      "  --ts-metric-scaling fac      scale the metric which defines the ellipsoids for TrigScan\n"\
      "                               Scaling must be > 0\n"\
      "\n"\
      "[LIGOLW XML input files] list of the input trigger files.\n"\
      "\n", program);
}


int main( int argc, char *argv[] )
{
  static LALStatus      status;

  INT4  startTime = -1;
  LIGOTimeGPS startTimeGPS = {0,0};
  INT4  endTime = -1;
  LIGOTimeGPS endTimeGPS = {0,0};
  CHAR  inputIFO[LIGOMETA_IFO_MAX];
  CHAR  outputIFO[LIGOMETA_IFO_MAX];
  CHAR  comment[LIGOMETA_COMMENT_MAX];
  CHAR *userTag = NULL;
  CHAR *ifoTag = NULL;

  CHAR  fileName[FILENAME_MAX];

  INT4  numTriggers = 0;

  SnglInspiralTable    *inspiralEventList=NULL;
  SnglInspiralTable    *currentTrigger = NULL;

  SearchSummvarsTable  *inputFiles = NULL;
  SearchSummvarsTable  *thisInputFile = NULL;

  SearchSummaryTable   *searchSummList = NULL;
  SearchSummaryTable   *thisSearchSumm = NULL;

  MetadataTable         proctable;
  MetadataTable         processParamsTable;
  MetadataTable         searchsumm;
  MetadataTable         searchSummvarsTable;
  MetadataTable         inspiralTable;
  ProcessParamsTable   *this_proc_param = NULL;
  LIGOLwXMLStream       xmlStream;
  INT4                  outCompress = 0;

  long int              gpstime;

  trigScanType trigScanMethod = trigScanNone;
                                          /* Switch for clustering        */
                                          /* triggers in template         */
                                          /* parameters and end time      */
  REAL8  trigScanMetricScalingFac = -1.0;
  /* Use this scaling factor for the volume spanned by a trigger in the   */
  /* parameter space. When set to x, the volume is taken to be that of the*/
  /* ambiguity ellipsoid at a 'minimal match' of (1.0-x).                 */
                                          /* original bank entered at the */
                                          /* command line                 */
  INT2  trigScanAppendStragglers = -1;    /* Switch to append cluster     */
                                          /* out-liers (stragglers)       */

  INT4                  i;

  /* getopt arguments */
  struct option long_options[] =
  {
    {"verbose",                no_argument,           &vrbflg,            1 },
    {"write-compress",         no_argument,           &outCompress,       1 },
    {"gps-start-time",         required_argument,     0,                 'q'},
    {"gps-end-time",           required_argument,     0,                 'r'},
    {"help",                   no_argument,           0,                 'h'}, 
    {"version",                no_argument,           0,                 'V'},
    {"user-tag",               required_argument,     0,                 'Z'},
    {"ifo",                    required_argument,     0,                 'I'},
    {"ts-cluster",             required_argument,     0,                 '*'},
    {"ts-metric-scaling",      required_argument,     0,                 '>'},
    {0, 0, 0, 0}
  };
  int c;

  /*
   * 
   * initialize things
   *
   */

  lal_errhandler = LAL_ERR_EXIT;
  setvbuf( stdout, NULL, _IONBF, 0 );

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
  XLALGPSTimeNow(&(proctable.processTable->start_time));
  XLALPopulateProcessTable(proctable.processTable, PROGRAM_NAME, lalAppsVCSIdentId,
      lalAppsVCSIdentStatus, lalAppsVCSIdentDate, 0);
  this_proc_param = processParamsTable.processParamsTable = 
    (ProcessParamsTable *) calloc( 1, sizeof(ProcessParamsTable) );
  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );

  /* create the search summary and zero out the summvars table */
  searchsumm.searchSummaryTable = (SearchSummaryTable *)
    calloc( 1, sizeof(SearchSummaryTable) );


  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only( argc, argv, 
        "a:b:hq:r:s:A:I:VZ:", long_options, 
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

      case 'h':
        /* help message */
        print_usage(argv[0]);
        exit( 1 );
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

        snprintf(inputIFO,  LIGOMETA_IFO_MAX, "%s", optarg);
        snprintf(outputIFO, LIGOMETA_IFO_MAX, "%s", optarg);
        break;
      case 'V':
        /* print version information and exit */
        fprintf( stdout, "TrigScan Cluster \n" 
            "Larne Pekowsky\n"
            "Based on trigbank and inspiral by Patrick Brady, Duncan Brown and Steve Fairhurst\n");
        XLALOutputVersionString(stderr, 0);
        exit( 0 );
        break;

      case '?':
        print_usage(argv[0]);
        exit( 1 );
        break;
      case 'q':
        /* start time */
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

      case 'r':
        /* end time  */
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

      case '*':
        /* store trigSanClustering method */
        if ( ! strcmp( "T0T3Tc", optarg ) )
        {
            trigScanMethod = T0T3Tc;
            trigScanAppendStragglers = 0;
        }
        else if ( ! strcmp( "T0T3TcAS", optarg ) )
        {
            trigScanMethod = T0T3Tc;
            trigScanAppendStragglers = 1;
        }
        else if ( ! strcmp( "Psi0Psi3Tc", optarg ) )
        {
            trigScanMethod = Psi0Psi3Tc;
            trigScanAppendStragglers = 0;
        }
        else if ( ! strcmp( "Psi0Psi3TcAS", optarg ) )
        {
            trigScanMethod = Psi0Psi3Tc;
            trigScanAppendStragglers = 1;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown scan method specified: %s\n"
              "(Must be one of T0T3Tc, T0T3TcAS, Psi0Psi3Tc, Psi0Psi3TcAS)\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case '>':
        /* TrigScan Template Metric Scaling Factor */
        trigScanMetricScalingFac = atof( optarg );
        if ( trigScanMetricScalingFac <= 0.0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "ts-volume-safety must be > 0.0 : "
              "(%f specified)\n",
              long_options[option_index].name, trigScanMetricScalingFac );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%s", optarg );
        break;
      default:
        fprintf( stderr, "Error: Unknown error while parsing options\n" );
        print_usage(argv[0]);
        exit( 1 );
    }
  }


  /* check the values of the arguments */
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

  /* Check the trigScan input parameters */
  if ( ! trigScanMethod )
  {
      fprintf ( stderr, "You must specify --ts-method\n" );
      exit(1);
  }

  if ( trigScanMetricScalingFac <= 0.0 )
  {
    fprintf ( stderr, "You must specify --ts-metric-scaling\n" );
    exit(1);
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
  else
  {
    fprintf( stderr, "Error: No trigger files specified.\n" );
    exit( 1 );
  }

  if ( vrbflg ) fprintf( stdout, "Read in a total of %d triggers.\n",
      numTriggers );


  if ( ! inspiralEventList )
  {
    /* no triggers read in so triggered bank will be empty */
    fprintf( stdout, "No triggers read in\n");
    exit( 0 );
  }


  /* trigScanClustering */ 
  /* Call the clustering routine */ 
  if (XLALTrigScanClusterTriggers( &(inspiralEventList),
                                 trigScanMethod,
                                 trigScanMetricScalingFac,
                                 trigScanAppendStragglers ) == XLAL_FAILURE )
  {
    fprintf( stderr, "New trig scan has failed!!\n" );
    exit(1);
  }

  /* time sort the triggers */
  if ( vrbflg ) fprintf( stdout, "Sorting triggers\n" );
  LAL_CALL( LALSortSnglInspiral( &status, &inspiralEventList,
        LALCompareSnglInspiralByTime ), &status );

  if( !inspiralEventList )
  {
    if ( vrbflg ) fprintf( stdout, 
        "No triggers remain after time and playground cuts.\n" );

    /* set numTriggers after cuts were applied */
    numTriggers = 0;
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
  searchsumm.searchSummaryTable->nevents = numTriggers;

  if ( vrbflg ) fprintf( stdout, "writing output file... " );

  /* set the file name correctly */
  if ( userTag && !outCompress )
 {
    snprintf( fileName, FILENAME_MAX, "%s-TRIGSCAN_%s-%d-%d.xml", 
        outputIFO, userTag, startTime, endTime - startTime );
  }
  else if ( !userTag && !outCompress )
  {
    snprintf( fileName, FILENAME_MAX, "%s-TRIGSCAN_%d-%d.xml", 
        outputIFO, startTime, endTime - startTime );
  }
  else if ( userTag && outCompress )
  {
    snprintf( fileName, FILENAME_MAX, "%s-TRIGSCAN_%s-%d-%d.xml.gz",
        outputIFO, userTag, startTime, endTime - startTime );
  }
  else 
  {
    snprintf( fileName, FILENAME_MAX, "%s-TRIGSCAN_%d-%d.xml.gz",
        outputIFO, startTime, endTime - startTime );
  }

  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status , &xmlStream, fileName ), 
      &status );

  /* write process table */
  snprintf( proctable.processTable->ifos, LIGOMETA_IFOS_MAX, "%s", 
      inputIFO );
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

  /* write the sngl_inspiral table */
  if ( inspiralEventList )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status ,&xmlStream, 
          sngl_inspiral_table), &status );
    inspiralTable.snglInspiralTable = inspiralEventList;
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


  while ( inspiralEventList )
  {
    currentTrigger = inspiralEventList;
    inspiralEventList = inspiralEventList->next;
    LAL_CALL( LALFreeSnglInspiral( &status, &currentTrigger ), &status );
  }

  if ( userTag ) free( userTag );
  if ( ifoTag ) free( ifoTag );

  if ( vrbflg ) fprintf( stdout, "done\n" );

  LALCheckMemoryLeaks();

  exit( 0 );
}
