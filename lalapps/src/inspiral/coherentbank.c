/*----------------------------------------------------------------------- 
 * 
 * File Name: coherent_bank.c
 *
 * Author: Seader, S.E.
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
#include <lalapps.h>
#include <processtable.h>

RCSID("$Id$");

#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "lalapps_coherentbank"

#define COHBANK_EARG   1
#define COHBANK_EROW   2
#define COHBANK_EFILE  3

#define COHBANK_MSGEARG   "Error parsing arguments"
#define COHBANK_MSGROW    "Error reading row from XML table"
#define COHBANK_MSGEFILE  "Could not open file"

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
      "  [--debug-level]   level       set the LAL debug level to LEVEL\n"\
      "  [--user-tag]      usertag     set the process_params usertag\n"\
      "  [--comment]       string      set the process table comment\n"\
      "  [--bank-file]     file        the input trigger file.\n"\
      "\n", program);
}


int main( int argc, char *argv[] )
{
  static LALStatus      status;
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;

  extern int vrbflg;

  INT8 eventID = 0;
  INT4 j = 0;
  INT4 i = 0;
  INT4 k = 0;
  INT4 inputFileNameFlag = 0;
  INT4 inputFileNameLength = 0;
  REAL4 snrTemp = 0.0;
  REAL4 mass1 = 0.0;
  REAL4 mass2 = 0.0;

  CHAR tempStr[FILENAME_MAX];

  CHAR  comment[LIGOMETA_NAME_MAX];
  CHAR *userTag = NULL;

  CHAR  fileName[FILENAME_MAX];
  CHAR  inputFileName[FILENAME_MAX];

  SnglInspiralTable    *inspiralEventList=NULL;
  SnglInspiralTable    *currentTrigger = NULL;
  SnglInspiralTable    *tempTable = NULL;
  SnglInspiralTable    *newEventList = NULL;

  SearchSummvarsTable  *inputFiles = NULL;
  SearchSummvarsTable  *thisInputFile = NULL;

  SearchSummaryTable   *searchSummList = NULL;
  SearchSummaryTable   *thisSearchSumm = NULL;

  MetadataTable         proctable;
  MetadataTable         processParamsTable;
  MetadataTable         searchsumm;
  MetadataTable		searchSummvarsTable;
  MetadataTable         inspiralTable;
  ProcessParamsTable   *this_proc_param = NULL;
  LIGOLwXMLStream       xmlStream;


  /* getopt arguments */
  struct option long_options[] =
  {
    {"verbose",                no_argument,     &vrbflg,                  1 },
    {"comment",                required_argument,     0,                 's'},
    {"user-tag",               required_argument,     0,                 'Z'},
    {"help",                   no_argument,           0,                 'h'}, 
    {"debug-level",            required_argument,     0,                 'z'},
    {"version",                no_argument,           0,                 'V'},
    {"bank-file",              required_argument,     0,                 'b'},
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
  /*  setvbuf( stdout, NULL, _IONBF, 0 );*/

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


  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only( argc, argv, 
        "hs:z:b:Z:V", long_options, 
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

      case 'b':
	strcpy(inputFileName, optarg);
	inputFileNameFlag = 1;
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

      case 'V':
        /* print version information and exit */
        fprintf( stdout, "Inspiral Triggered Bank Generator\n" 
            "Patrick Brady, Duncan Brown and Steve Fairhurst\n"
            "CVS Version: " CVS_ID_STRING "\n"
            "CVS Tag: " CVS_NAME_STRING "\n" );
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

  if (optind < argc)
    {
      fprintf( stderr, "extraneous command line arguments:\n" );
      while ( optind < argc )
	{
          fprintf ( stderr, "%s\n", argv[optind++] );
        }
      exit( 1 );      
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

  if( !inputFileNameFlag )
    {
      fprintf(stderr,"You must specify a bank file. Exiting.\n");
      exit(1);
    }

    INT4 haveSearchSum = 0;
    INT4 numFileTriggers = 0;
    SnglInspiralTable  *inputData = NULL;
    SearchSummaryTable *inputSummary = NULL;


    /* store the file name in search summvars */
    if ( vrbflg ) fprintf( stdout, 
        "storing input file name %s in search summvars table\n", inputFileName );

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
        "%s", inputFileName );      


    /* read in the search summary and store */ 
    if ( vrbflg ) fprintf( stdout, 
        "reading search_summary table from file: %s\n", inputFileName );

    haveSearchSum = SearchSummaryTableFromLIGOLw( &inputSummary, inputFileName );

    if ( haveSearchSum < 1 || ! inputSummary )
    {
      if ( vrbflg ) 
        fprintf( stdout, "no valid search_summary table, exiting\n" );
      exit( 1 );
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

    /* read in the triggers */
    if ( vrbflg ) 
      fprintf( stdout, "reading triggers from file: %s\n", inputFileName );

    numFileTriggers = 
      LALSnglInspiralTableFromLIGOLw( &inputData, inputFileName, 0, -1 );
      
    if ( numFileTriggers < 0 )
    {
      fprintf( stderr, "error: unable to read sngl_inspiral table from %s\n", 
          inputFileName );
      exit( 1 );
    }
    else if ( numFileTriggers > 0 )
    {
      if ( vrbflg ) 
        fprintf( stdout, "got %d sngl_inspiral rows from %s\n" , 
            numFileTriggers, inputFileName );

        /* store the triggers */
	inspiralEventList = currentTrigger = inputData;
    


        /* store the values of masses corresponding to the loudest trigger */
	/* for each unique event ID */

        REAL4 massMatrix[numFileTriggers][2];
	for( i = 0; i < numFileTriggers; i++)
	  {
	    /* get the values from first snglInspiralTable */
	    if( i == 0 )
	      {
		eventID = inspiralEventList->event_id->id;
		snrTemp = inspiralEventList->snr;
		mass1 = inspiralEventList->mass1;
		mass2 = inspiralEventList->mass2;
		j++;
		inspiralEventList = inspiralEventList->next;
	      }
	    else if( inspiralEventList->event_id->id == (UINT8)eventID)
	      {
		/* update values if greater snr for same event ID */ 
		if( inspiralEventList->snr > snrTemp )
		  {
		    snrTemp = inspiralEventList->snr;
		    mass1 = inspiralEventList->mass1;
		    mass2 = inspiralEventList->mass2;
		    j++;
		    inspiralEventList = inspiralEventList->next;
		  }
		else
		  {
		    /* do nothing if snr is not greater for same event ID */
		    inspiralEventList = inspiralEventList->next;
		    j++;
		  }
		
		if( i == numFileTriggers - 1 )
		  {
		    /* store masses for the last event */
		    for( k = i-j+1; k < i+1; k++ )
		      {
			massMatrix[k][0] = mass1;
			massMatrix[k][1] = mass2;
		      }
		  }

	      }
	    else
	      {
		/* store masses in mass matrix */
		for( k = i-j; k < i; k++ )
		  {
		    massMatrix[k][0] = mass1;
		    massMatrix[k][1] = mass2;
		  }
		j = 1;
		k = 0;
		mass1 = inspiralEventList->mass1;
		mass2 = inspiralEventList->mass2;
		snrTemp = inspiralEventList->snr;
		eventID = inspiralEventList->event_id->id;
		inspiralEventList = inspiralEventList->next;
	       
	      }
	  }

	/* Now generate an event table with the new masses */
	
	newEventList = tempTable = inputData;
	for( i = 0; i < numFileTriggers; i++ )
	  {
	    tempTable->mass1 = massMatrix[i][0];
	    tempTable->mass2 = massMatrix[i][1];
	    tempTable = tempTable->next;
	  }
      }
      else
      { 
	fprintf( stderr, "%s contains no triggers - the cohernt bank will be empty\n",inputFileName );
      }
  

  goto cleanexit;

  /*
   *
   * write the output xml file
   *
   */


cleanexit:

  /* search summary entries: nevents is from primary ifo */
 
  if ( vrbflg ) fprintf( stdout, "writing output file... " );

  /* set the file name correctly */
  inputFileNameLength = strlen( inputFileName );
  memcpy(tempStr, inputFileName, inputFileNameLength -4);
 
  if ( userTag )
  {
    LALSnprintf( fileName, FILENAME_MAX, "%s-COHBANK_%s.xml", 
        tempStr, userTag );
  }
  else
  {
    LALSnprintf( fileName, FILENAME_MAX, "%s-COHBANK.xml", 
        tempStr);
  }
  memset( &xmlStream, 0, sizeof(LIGOLwXMLStream) );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status , &xmlStream, fileName), 
      &status );

  /* write process table */
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

  /* write the search_summvars tabls */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status ,&xmlStream, 
        search_summvars_table), &status );
  searchSummvarsTable.searchSummvarsTable = inputFiles;
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlStream, searchSummvarsTable,
        search_summvars_table), &status );
  LAL_CALL( LALEndLIGOLwXMLTable( &status, &xmlStream), &status );

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


  while ( inspiralEventList )
  {
    currentTrigger = inspiralEventList;
    inspiralEventList = inspiralEventList->next;
    LALFree( currentTrigger );
  }
  
 while ( newEventList )
  {
    tempTable = newEventList;
    newEventList = newEventList->next;
    LALFree( tempTable );
  }

  if ( userTag ) free( userTag );

  if ( vrbflg ) fprintf( stdout, "done\n" );

  LALCheckMemoryLeaks();

  exit( 0 );
}
