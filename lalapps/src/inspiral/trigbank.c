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
 * Author: Brady, P. R., Brown, D. A. and Fairhurst, S.
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
#include <lal/LIGOMetadataInspiralUtils.h>
#include <lalapps.h>
#include <processtable.h>
#include <LALAppsVCSInfo.h>

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "trigbank"

#define TRIGBANK_EARG   1
#define TRIGBANK_EROW   2
#define TRIGBANK_EFILE  3

#define TRIGBANK_MSGEARG   "Error parsing arguments"
#define TRIGBANK_MSGROW    "Error reading row from XML table"
#define TRIGBANK_MSGEFILE  "Could not open file"

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
  this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
calloc( 1, sizeof(ProcessParamsTable) ); \
snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
    PROGRAM_NAME ); \
snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
    long_options[option_index].name ); \
snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
snprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

int checkTimes = 0;
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
      "  [--ifo-tag]       ifotag      set the ifo-tag - for file naming\n"\
      "  [--comment]       string      set the process table comment\n"\
      "  [--write-compress]            write a compressed xml file\n"\
      "   --gps-start-time start_time  GPS second of data start time\n"\
      "   --gps-end-time   end_time    GPS second of data end time\n"\
      "  [--check-times]               Check that all times were analyzed\n"\
      "   --input-ifo      input_ifo   the name of the input IFO triggers\n"\
      "   --output-ifo     output_ifo  the name of the IFO for which to create the bank\n"\
      "   --parameter-test test        set parameters with which to test coincidence:\n"\
      "                                (m1_and_m2|psi0_and_psi3|no_test)\n"\
      "  --data-type DATA_TYPE         specify the data type, must be one of\n"\
      "                                (playground_only|exclude_play|all_data)\n"\
      "  --coherent-run    coherentRun      run in coherent stage\n"\
      "  --coherent-buffer coherentBuffer   time-buffer (sec) for vetoing triggers\n"\
      "                                     in corrupted data\n"\
      "\n"\
      "[LIGOLW XML input files] list of the input trigger files.\n"\
      "\n", program);
}


int main( int argc, char *argv[] )
{
  static LALStatus      status;

  LALPlaygroundDataMask  dataType = unspecified_data_type;
  SnglInspiralParameterTest  test = unspecified_test;

  INT4  startTime = -1;
  LIGOTimeGPS startTimeGPS = {0,0};
  INT4  endTime = -1;
  LIGOTimeGPS endTimeGPS = {0,0};
  INT4  startChunkTime = -1;
  LIGOTimeGPS startChunkTimeGPS = {0,0};
  INT4  endChunkTime = -1;
  LIGOTimeGPS endChunkTimeGPS = {0,0};
  int   coherentRun = 0;
  INT4  coherentBuffer = 64;
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

  INT4                  i;

  /* getopt arguments */
  struct option long_options[] =
  {
    {"verbose",                no_argument,     &vrbflg,                  1 },
    {"check-times",            no_argument, &checkTimes,                  1 },
    {"write-compress",         no_argument,           &outCompress,       1 },
    {"input-ifo",              required_argument,     0,                 'a'},
    {"output-ifo",             required_argument,     0,                 'b'},
    {"parameter-test",         required_argument,     0,                 'A'},
    {"data-type",              required_argument,     0,                 'D'},
    {"gps-start-time",         required_argument,     0,                 'q'},
    {"gps-end-time",           required_argument,     0,                 'r'},
    {"comment",                required_argument,     0,                 's'},
    {"coherent-run",           no_argument,           &coherentRun,       1 },
    {"coherent-buffer",        required_argument,     &coherentBuffer,   't'},
    {"user-tag",               required_argument,     0,                 'Z'},
    {"userTag",                required_argument,     0,                 'Z'},
    {"ifo-tag",                required_argument,     0,                 'I'},
    {"help",                   no_argument,           0,                 'h'},
    {"version",                no_argument,           0,                 'V'},
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
    long int gpstime;
    size_t optarg_len;

    c = getopt_long_only( argc, argv,
        "a:b:hq:r:s:t:A:I:VZ:", long_options,
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
        /* name of input ifo*/
        strncpy( inputIFO, optarg, LIGOMETA_IFO_MAX );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'b':
        /* name of output ifo */
        strncpy( outputIFO, optarg, LIGOMETA_IFO_MAX );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case 'A':
        /* comparison used to test for uniqueness of triggers */
        if ( ! strcmp( "m1_and_m2", optarg ) )
        {
          test = m1_and_m2;
        }
        else if ( ! strcmp( "psi0_and_psi3", optarg ) )
        {
          test = psi0_and_psi3;
        }
        else if ( ! strcmp( "mchirp_and_eta", optarg ) )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "mchirp_and_eta test specified, not implemented for trigbank: "
              "%s (must be m1_and_m2, psi0_and_psi3, no_test)\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        else if ( ! strcmp( "no_test", optarg ) )
        {
          test = no_test;
        }
        else
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "unknown test specified: "
              "%s (must be m1_and_m2, psi0_and_psi3,no_test, or mchirp_and_eta)\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "string", "%s", optarg );
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

      case 't':
        coherentBuffer = (INT4) atoi( optarg );
        if ( coherentBuffer < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "coherent-buffer duration (sec) must be positive or zero: "
              "(%d specified)\n",
              long_options[option_index].name, coherentBuffer );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", coherentBuffer );
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
        break;

      case 'V':
        /* print version information and exit */
        fprintf( stdout, "Inspiral Triggered Bank Generator\n" 
            "Patrick Brady, Duncan Brown and Steve Fairhurst\n");
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

  if ( dataType == unspecified_data_type )
  {
    fprintf( stderr, "Error: --data-type must be specified\n");
    exit(1);
  }
  
  if ( test == unspecified_test )
  {
      fprintf( stderr, "Error: --parameter-test must be specified\n");
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

  /* store the check-times in the process_params table */
  if ( checkTimes )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
        "%s", PROGRAM_NAME );
    snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--check-times");
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

  /* check that we have read in data for all the requested time */
  if ( checkTimes )
  {
    if ( vrbflg ) fprintf( stdout, 
        "Checking that we have data for all times from %s\n",
        inputIFO);
    LAL_CALL( LALCheckOutTimeFromSearchSummary ( &status, searchSummList, 
          inputIFO, &startTimeGPS, &endTimeGPS ), &status);
  }


  if ( ! inspiralEventList )
  {
    /* no triggers read in so triggered bank will be empty */
    fprintf( stdout, "No triggers read in\n");

    /* set numTriggers in case any cuts are made in future */
    numTriggers = 0;
    goto cleanexit;
  }


  /* keep only triggers from input ifo */
  LAL_CALL( LALIfoCutSingleInspiral( &status, &inspiralEventList, inputIFO ), 
      &status );

  /* time sort the triggers */
  if ( vrbflg ) fprintf( stdout, "Sorting triggers\n" );
  LAL_CALL( LALSortSnglInspiral( &status, &inspiralEventList,
        LALCompareSnglInspiralByTime ), &status );

  /* keep only triggers within the requested interval */
  if ( vrbflg ) fprintf( stdout, 
      "Discarding triggers outside requested interval\n" );
  if ( coherentRun ) {
    startChunkTime = startTimeGPS.gpsSeconds + coherentBuffer;
    startChunkTimeGPS.gpsSeconds = startChunkTime;
    endChunkTime = endTimeGPS.gpsSeconds - coherentBuffer;
    endChunkTimeGPS.gpsSeconds = endChunkTime;
    LAL_CALL( LALTimeCutSingleInspiral( &status, &inspiralEventList,
         &startChunkTimeGPS, &endChunkTimeGPS), &status );
  }
  else {
    LAL_CALL( LALTimeCutSingleInspiral( &status, &inspiralEventList,
         &startTimeGPS, &endTimeGPS), &status );
  }

  /* keep play/non-play/all triggers */
  if ( dataType == playground_only && vrbflg ) fprintf( stdout, 
      "Keeping only playground triggers\n" );
  else if ( dataType == exclude_play && vrbflg ) fprintf( stdout, 
      "Keeping only non-playground triggers\n" );
  else if ( dataType == all_data && vrbflg ) fprintf( stdout, 
      "Keeping all triggers\n" );
  LAL_CALL( LALPlayTestSingleInspiral( &status, &inspiralEventList,
        &dataType ), &status );

  if( !inspiralEventList )
  {
    if ( vrbflg ) fprintf( stdout, 
        "No triggers remain after time and playground cuts.\n" );

    /* set numTriggers after cuts were applied */
    numTriggers = 0;
    goto cleanexit;
  }

  /* Generate the triggered bank */
  if( test != no_test )
  {
    LAL_CALL( LALCreateTrigBank( &status, &inspiralEventList, &test ), 
        &status );
  }

  /* count the number of triggers  */
  for( currentTrigger = inspiralEventList, numTriggers = 0; currentTrigger; 
      currentTrigger = currentTrigger->next, ++numTriggers);

  if ( vrbflg ) fprintf( stdout, "%d triggers to be written to trigbank.\n",
      numTriggers );

  /*
   *
   * write the output xml file
   *
   */


cleanexit:

  /* search summary entries: */
  searchsumm.searchSummaryTable->in_start_time = startTimeGPS;
  searchsumm.searchSummaryTable->in_end_time = endTimeGPS;
  searchsumm.searchSummaryTable->out_start_time = startTimeGPS;
  searchsumm.searchSummaryTable->out_end_time = endTimeGPS;
  searchsumm.searchSummaryTable->nevents = numTriggers;

  if ( vrbflg ) fprintf( stdout, "writing output file... " );

  /* set the file name correctly */
  if ( userTag && ifoTag && !outCompress )
  {
    snprintf( fileName, FILENAME_MAX, "%s-TRIGBANK_%s_%s-%d-%d.xml", 
        outputIFO, ifoTag, userTag, startTime, endTime - startTime );
  }
  else if ( userTag && !ifoTag && !outCompress )
  {
    snprintf( fileName, FILENAME_MAX, "%s-TRIGBANK_%s-%d-%d.xml", 
        outputIFO, userTag, startTime, endTime - startTime );
  }
  else if ( !userTag && ifoTag && !outCompress )
  {
    snprintf( fileName, FILENAME_MAX, "%s-TRIGBANK_%s-%d-%d.xml", 
        outputIFO, ifoTag, startTime, endTime - startTime );
  }
  else if ( userTag && ifoTag && outCompress )
  {
    snprintf( fileName, FILENAME_MAX, "%s-TRIGBANK_%s_%s-%d-%d.xml.gz",
        outputIFO, ifoTag, userTag, startTime, endTime - startTime );
  }
  else if ( userTag && !ifoTag && outCompress )
  {
    snprintf( fileName, FILENAME_MAX, "%s-TRIGBANK_%s-%d-%d.xml.gz",
        outputIFO, userTag, startTime, endTime - startTime );
  }
  else if ( !userTag && ifoTag && outCompress )
  {
    snprintf( fileName, FILENAME_MAX, "%s-TRIGBANK_%s-%d-%d.xml.gz",
        outputIFO, ifoTag, startTime, endTime - startTime );
  }
  else if ( !userTag && !ifoTag && outCompress )
  {
    snprintf( fileName, FILENAME_MAX, "%s-TRIGBANK-%d-%d.xml.gz",
        outputIFO, startTime, endTime - startTime );
  }
  else 
  {
    snprintf( fileName, FILENAME_MAX, "%s-TRIGBANK-%d-%d.xml", 
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
