/*----------------------------------------------------------------------- 
 * 
 * File Name: splitbank.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <regex.h>
#include <time.h>

#include <lalapps.h>
#include <processtable.h>

#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/Date.h>

RCSID( "$Id$" );
#define CVS_ID_STRING "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "splitbank"

#define USAGE \
"Usage: %s [options] [LIGOLW XML input files]\n\n"\
"  --help                    display this message\n"\
"  --verbose                 print progress information\n"\
"  --version                 print version information\n"\
"  --debug-level LEVEL       set the LAL debug level to LEVEL\n"\
"  --user-tag STRING         set the process_params usertag to STRING\n"\
"  --comment STRING          set the process table comment to STRING\n"\
"\n"\
"  --input-bank FILE         read templates from FILE\n"\
"  --number-of-banks N       split template bank into N files\n"


/*
 *
 * variables that control program behaviour
 *
 */


/* debugging */
extern int vrbflg;                      /* verbocity of lal function    */

/* template bank generation parameters */
CHAR    inputBankFile[4096];            /* name of the input template bank */
INT4    numOutBanks = 0;                /* number of banks to split into   */

int main ( int argc, char *argv[] )
{
  /* lal function variables */
  LALStatus             status = blank_status;
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;

  /* output data */
  /* MetadataTable         inputBank; */
  /* MetadataTable         outputBank; */
  MetadataTable         proctable;
  MetadataTable         procparams;
  ProcessParamsTable   *this_proc_param = NULL;
  /* LIGOLwXMLStream       results; */

  /* counters and other variables */
  /* INT4 i; */
  /* char fname[4096]; */
  CHAR comment[LIGOMETA_COMMENT_MAX];  
  CHAR *userTag = NULL;


  /* getopt arguments */
  struct option long_options[] =
  {
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"version",                 no_argument,       0,                'V'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    {"comment",                 required_argument, 0,                's'},    
    {"help",                    no_argument,       0,                'h'}, 
    {"debug-level",             required_argument, 0,                'z'},
    {"input-bank",              required_argument, 0,                'i'},
    {"number-of-banks",         required_argument, 0,                'n'},
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
  setvbuf( stdout, NULL, _IONBF, 0 );

  memset( inputBankFile, 0, sizeof(inputBankFile) * sizeof(*inputBankFile) );

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) calloc( 1, sizeof(ProcessTable) );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->start_time),
        &accuracy ), &status );
  LAL_CALL( populate_process_table( &status, proctable.processTable, 
        PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE ), &status );
  this_proc_param = procparams.processParamsTable = 
    (ProcessParamsTable *) calloc( 1, sizeof(ProcessParamsTable) );
  memset( comment, 0, LIGOMETA_COMMENT_MAX * sizeof(CHAR) );


  /*
   *
   * parse command line arguments
   *
   */


  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    size_t optarg_len;

    c = getopt_long_only( argc, argv, 
        "i:n:VZ:hz:s:", 
        long_options, &option_index );

    /* detect the end of the options */
    if ( c == - 1 )
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
          fprintf( stderr, "error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        break;

      case 'i':
        strncpy( inputBankFile, 
            optarg, sizeof(inputBankFile) * sizeof(*inputBankFile) );
        LALSnprintf( procparams.processParamsTable->program, 
            LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
        LALSnprintf( procparams.processParamsTable->type, 
            LIGOMETA_TYPE_MAX, "string" );
        LALSnprintf( procparams.processParamsTable->param, 
            LIGOMETA_PARAM_MAX, "--%s", long_options[option_index].name );
        LALSnprintf( procparams.processParamsTable->value, 
            LIGOMETA_TYPE_MAX, "%s", optarg );
        break;

      case 'n':
        numOutBanks = (INT4) atoi( optarg );
        if ( numOutBanks < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "Number of output banks must be greater than zero:" 
              "(%d specified)\n",
              long_options[option_index].name, numOutBanks );
          exit( 1 );
        }
        else if ( numOutBanks > 20 )
        {
          fprintf( stderr, 
              "Warning: generating more than 20 banks is not reccomended!\n" );
        }
        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
            "%s", PROGRAM_NAME );
        LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
            "--%s", long_options[option_index].name );
        LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "int" );
        LALSnprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, "%d", 
            numOutBanks );
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
          LALSnprintf( comment, LIGOMETA_COMMENT_MAX, "%s", optarg );
        }
        break;

      case 'z':
        set_debug_level( optarg );
        this_proc_param = this_proc_param->next = (ProcessParamsTable *)
          calloc( 1, sizeof(ProcessParamsTable) );
        LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", 
            PROGRAM_NAME );
        LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
            "--debug-level" );
        LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
        LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, "%s",
            optarg );
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

      case 'V':
        /* print version information and exit */
        fprintf( stdout, "Template Bank Splitter\n" 
            "Duncan Brown <duncan@gravity.phys.uwm.edu>\n"
            "CVS Version: " CVS_ID_STRING "\n" );
        exit( 0 );
        break;

      case '?':
        fprintf( stderr, USAGE, argv[0] );
        exit( 1 );
        break;

      default:
        fprintf( stderr, "unknown error while parsing options\n" );
        fprintf( stderr, USAGE, argv[0] );
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


  /* check the values of the arguments */
  if ( ! inputBankFile[0] )
  {
    fprintf( stderr, "Error: --input-bank must be specified\n" );
    exit( 1 );
  }

  if ( ! numOutBanks )
  {
    fprintf( stderr, "Error: --number-of-banks must be specified\n" );
    exit( 1 );
  }


  return 0;
}
