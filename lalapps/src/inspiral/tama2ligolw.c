/*----------------------------------------------------------------------- 
 * 
 * File Name: tama2ligolw.c
 *
 * Author:  Fairhurst, S
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
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/Date.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/StreamInput.h>
#include <lalapps.h>
#include <processtable.h>

RCSID("$Id$");

#define PROGRAM_NAME "tama2ligolw"
#define CVS_ID_STRING "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"

#define USAGE \
"Usage: lalapps_tama2ligolw [options]\n"\
"\n"\
"  --help                       display this message\n"\
"  --verbose                    print progress information\n"\
"  --debug-level LEVEL          set the LAL debug level to LEVEL\n"\
"  --user-tag STRING            set the process_params usertag to STRING\n"\
"  --comment STRING             set the process table comment to STRING\n"\
"  --version                    print the CVS version string\n"\
"\n"\
"  --gps-start-time SEC         start time (default to S2 start, 729273613)\n"\
"  --gps-start-time-ns NS       start time nanoseconds\n"\
"  --gps-end-time SEC           end time (default to S2 end, 734367613)\n"\
"  --gps-end-time-ns NS         end time nanoseconds\n"\
"  --input FILE                 TAMA file to be converted to xml\n"\
"  --output FILE                write output data to FILE\n"\

#define S2START 729273613
#define S2END 734367613

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
  this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
calloc( 1, sizeof(ProcessParamsTable) ); \
LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
    PROGRAM_NAME ); \
LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
    long_options[option_index].name ); \
LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

int main( int argc, char *argv[] )
{
  /* lal initialization variables */
  LALStatus stat = blank_status;
  LALLeapSecAccuracy accuracy = LALLEAPSEC_LOOSE;
  LALMSTUnitsAndAcc     gmstUnits;

  extern int vrbflg;
  CHAR *userTag = NULL;
  CHAR  comment[LIGOMETA_COMMENT_MAX];
  CHAR  line[LIGOMETA_STRING_MAX];

  char *inputFileName = NULL;
  char *outputFileName = NULL;

  FILE *fp;
  INT4			num_trigs = 0;
  MetadataTable         proctable;
  MetadataTable         procparams;
  ProcessParamsTable   *this_proc_param;
  LIGOLwXMLStream       xmlStream;
  MetadataTable         outputTable;
  MetadataTable         searchsumm;
  SnglInspiralTable    *eventHead = NULL;
  SnglInspiralTable    *thisEvent = NULL;
  REAL8                 gpsStartTimeFloat = 0;
  REAL8                 gpsEndTimeFloat = 0;
  LIGOTimeGPS           gpsStartTime = {0,0};
  LIGOTimeGPS           gpsEndTime = {0,0};
  /*
   *
   * initialization
   *
   */


  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "33" );

  /* set the gmst units and strictness */
  gmstUnits.units = MST_HRS;
  gmstUnits.accuracy = LALLEAPSEC_STRICT;

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
      {"verbose",                 no_argument,            &vrbflg,         1 },
      {"help",                    no_argument,            0,              'h'},
      {"debug-level",             required_argument,      0,              'z'},
      {"user-tag",                required_argument,      0,              'Z'},
      {"userTag",                 required_argument,      0,              'Z'},
      {"comment",                 required_argument,      0,              'c'},
      {"version",                 no_argument,            0,              'V'},
      {"input",                   required_argument,      0,              'i'},
      {"output",                  required_argument,      0,              'o'},
      {"gps-start-time",          required_argument,      0,              's'},
      {"gps-start-time-ns",       required_argument,      0,              'S'}, 
      {"gps-end-time",            required_argument,      0,              'e'},
      {"gps-end-time-ns",         required_argument,      0,              'E'},
      {0, 0, 0, 0}
    };
    int c;

    /* getopt_long stores the option index here. */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long ( argc, argv, "c:e:E:hi:o:s:S:VzZ:", 
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
	fprintf( stdout, "Tama text to LIGO Lw XML converter\n"
	    "Steve Fairhurst\n"
	    "CVS Version: " CVS_ID_STRING "\n" );
	exit( 0 );
	break;

      case 's':
	{
	  long int gstartt = atol( optarg );
	  if ( gstartt < 441417609 )
	  {
	    fprintf( stderr, "invalid argument to --%s:\n"
		"GPS start time is prior to " 
		"Jan 01, 1994  00:00:00 UTC:\n"
		"(%ld specified)\n",
		long_options[option_index].name, gstartt );
	    exit( 1 );
	  }
	  if ( gstartt > 999999999 )
	  {
	    fprintf( stderr, "invalid argument to --%s:\n"
		"GPS start time is after " 
		"Sep 14, 2011  01:46:26 UTC:\n"
		"(%ld specified)\n", 
		long_options[option_index].name, gstartt );
	    exit( 1 );
	  }
	  gpsStartTimeFloat += (REAL8) gstartt;
	  ADD_PROCESS_PARAM( "int", "%ld", gstartt );
	}
	break;

      case 'S':
	{
	  long int gstarttns = atol( optarg );
	  if ( gstarttns < 0 )
	  {
	    fprintf( stderr, "invalid argument to --%s:\n"
		"GPS start time nanoseconds is negative\n",
		long_options[option_index].name );
	    exit( 1 );
	  }
	  if ( gstarttns > 999999999 )
	  {
	    fprintf( stderr, "invalid argument to --%s:\n"
		"GPS start time nanoseconds is greater than unity:\n" 
		"Must be <= 999999999 (%ld specified)\n", 
		long_options[option_index].name, gstarttns );
	    exit( 1 );
	  }
	  gpsStartTimeFloat += (REAL8) 1.0e-09 * gstarttns;
	  ADD_PROCESS_PARAM( "int", "%ld", gstarttns );
	}
	break;

      case 'e':
	{
	  long int gendt = atol( optarg );
	  if ( gendt > 999999999 )
	  {
	    fprintf( stderr, "invalid argument to --%s:\n"
		"GPS end time is after " 
		"Sep 14, 2011  01:46:26 UTC:\n"
		"(%ld specified)\n", 
		long_options[option_index].name, gendt );
	    exit( 1 );
	  }
	  else if ( gendt < 441417609 )
	  {
	    fprintf( stderr, "invalid argument to --%s:\n"
		"GPS end time is prior to " 
		"Jan 01, 1994  00:00:00 UTC:\n"
		"(%ld specified)\n", 
		long_options[option_index].name, gendt );
	    exit( 1 );
	  }            
	  gpsEndTimeFloat += (REAL8) gendt;
	  ADD_PROCESS_PARAM( "int", "%ld", gendt );
	}
	break;

      case 'E':
	{
	  long int gendtns = atol( optarg );
	  if ( gendtns < 0 )
	  {
	    fprintf( stderr, "invalid argument to --%s:\n"
		"GPS end time nanoseconds is negative\n",
		long_options[option_index].name );
	    exit( 1 );
	  }
	  else if ( gendtns > 999999999 )
	  {
	    fprintf( stderr, "invalid argument to --%s:\n"
		"GPS end time nanoseconds is greater than unity:\n" 
		"Must be <= 999999999:\n"
		"(%ld specified)\n", 
		long_options[option_index].name, gendtns );
	    exit( 1 );
	  }            
	  gpsEndTimeFloat += (REAL8) 1.0e-09 * gendtns;
	  ADD_PROCESS_PARAM( "int", "%ld", gendtns );
	}
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




  /* check whether a start and or end time have been given.  If not, then set 
   * them to the default values */
  if( gpsStartTimeFloat == 0 )
  {
    gpsStartTimeFloat = S2START;
  }
  if( gpsEndTimeFloat == 0 )
  {
    gpsEndTimeFloat = S2END;
  }

  LAL_CALL( LALFloatToGPS( &stat, &gpsStartTime, &gpsStartTimeFloat ),
      &stat );
  LAL_CALL( LALFloatToGPS( &stat, &gpsEndTime, &gpsEndTimeFloat ),
      &stat );



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
  if ( ! inputFileName )
  {
    fprintf( stderr, "--input must be specified\n" );
    exit( 1 );
  }
  if ( ! outputFileName )
  {
    fprintf( stderr, "--output must be specified\n" );
    exit( 1 );
  }

  /*
   *
   * read in the tama trigger file and store as a linked list
   *
   */

  if ( vrbflg ) fprintf( stdout, 
      "reading TAMA triggers from %s... \n", inputFileName );
  fp = fopen( inputFileName, "r" );
  if ( fp == NULL )
  {
    fprintf( stderr, "Could not open file %s\n", inputFileName );
    exit ( 1 );
  }
  while( fgets( line, sizeof( line ), fp ) )
  {
    REAL8 trig_time = 0;
    REAL4 snr = 0;
    REAL4 chisq = 0;
    REAL4 mtot = 0;
    REAL4 eta = 0;
    REAL4 deff = 0;
    REAL4 col_six = 0;
    REAL4 col_seven = 0;
    REAL4 col_eight = 0;
    REAL4 col_nine = 0;
    REAL4 col_ten = 0;
    INT4 num_columns = 0;

    
    /* effective distance not required */
    num_columns =  sscanf( line, "%lf %f %f %f %f %f %f %f %f %f\n", 
	&trig_time, &snr, &chisq, &mtot, &eta, 
	&col_six, &col_seven, &col_eight, &col_nine, &col_ten ); 
    
    
    if( num_columns >= 5 )
    {
      if ( vrbflg )
      {
	fprintf( stdout, "obtained trigger at time %f\n", trig_time);
      }
    }
    else
    {
      fprintf( stderr, "Invalid line format\n" );
      exit ( 1 );
    }

    if( num_columns == 6 || num_columns == 10)
    {
      if ( vrbflg )
      {
	fprintf( stdout, "column six contains effective distance of %f", 
	    col_six);
      }
      deff = col_six;
    }
    
    /* check that trigger is within our analyzed time */
    if( trig_time >= gpsStartTimeFloat && trig_time <= gpsEndTimeFloat )
    {
      if( thisEvent )
      {
	thisEvent = thisEvent->next = (SnglInspiralTable *) 
	  LALCalloc( 1, sizeof(SnglInspiralTable) );
      }
      else
      {
	eventHead = thisEvent = (SnglInspiralTable *) 
	  LALCalloc( 1, sizeof(SnglInspiralTable) );
      }
      thisEvent->snr = snr;
      thisEvent->chisq = chisq;
      thisEvent->chisq_dof = 30;
      thisEvent->eta = eta;

      /* if eta is less than 0.25, store the masses of the components */
      if( eta <= 0.25 )
      {
	thisEvent->mass1 = mtot * ( 1 + sqrt( 1 - 4 * eta ) ) / 2;
	thisEvent->mass2 = mtot * ( 1 - sqrt( 1 - 4 * eta ) ) / 2;
      }
      else
      {
	thisEvent->mass1 = 0;
	thisEvent->mass2 = 0;
      } 
      thisEvent->mchirp = pow( eta, 0.6 ) * mtot;
      thisEvent->eff_distance = 1.0e-03 * deff; /* TAMA eff dist in kpc */

      LALSnprintf( thisEvent->ifo, LIGOMETA_IFO_MAX * sizeof(CHAR), "T1" );
      LAL_CALL( LALFloatToGPS( &stat, &(thisEvent->end_time), &trig_time ),
	  &stat );
      LAL_CALL( LALGPStoGMST1( &stat, &(thisEvent->end_time_gmst),
	    &(thisEvent->end_time),  &gmstUnits ), &stat );	  
      ++num_trigs;
    }
    else if( vrbflg )
    {
      fprintf( stdout, "trigger outside of analyzed time");
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
  if ( vrbflg )  fprintf( stdout, "process_params... " );
  LAL_CALL( LALBeginLIGOLwXMLTable( &stat, &xmlStream, 
	process_params_table ),	&stat );
  LAL_CALL( LALWriteLIGOLwXMLTable( &stat, &xmlStream, procparams, 
	process_params_table ), &stat );
  LAL_CALL( LALEndLIGOLwXMLTable ( &stat, &xmlStream ), &stat );

  while( procparams.processParamsTable )
  {
    this_proc_param = procparams.processParamsTable;
    procparams.processParamsTable = this_proc_param->next;
    free( this_proc_param );
  }

  /* create and write out the search summary table */
  /* create the search summary and zero out the summvars table */
  searchsumm.searchSummaryTable = (SearchSummaryTable *)
    calloc( 1, sizeof(SearchSummaryTable) );
  LALSnprintf( searchsumm.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX, 
      " " );
  searchsumm.searchSummaryTable->out_start_time = 
    searchsumm.searchSummaryTable->in_start_time = gpsStartTime;
  searchsumm.searchSummaryTable->out_end_time = 
    searchsumm.searchSummaryTable->in_end_time = gpsEndTime;
  searchsumm.searchSummaryTable->nevents = num_trigs;

  /* write the search summary table */
  if ( vrbflg ) fprintf( stdout, "  search_summary table...\n" );
  LAL_CALL( LALBeginLIGOLwXMLTable( &stat, &xmlStream, 
	search_summary_table ), &stat );
  LAL_CALL( LALWriteLIGOLwXMLTable( &stat, &xmlStream, searchsumm, 
	search_summary_table ), &stat );
  LAL_CALL( LALEndLIGOLwXMLTable ( &stat, &xmlStream ), &stat );

  /* free the search summary table */
  free( searchsumm.searchSummaryTable );


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

    /* free the temporary memory containing the events */
    while ( eventHead )
    {
      thisEvent = eventHead;
      eventHead = eventHead->next;
      LALFree( thisEvent );
    }
  }

  /* close the output file */
  LAL_CALL( LALCloseLIGOLwXMLFile(&stat, &xmlStream), &stat);
  if ( vrbflg ) fprintf( stdout, "done\n" );

  /*
   *
   * free memory and exit
   *
   */

  if ( vrbflg ) fprintf( stdout, "checking memory leaks and exiting\n" );
  LALCheckMemoryLeaks();
  exit( 0 );
}     
