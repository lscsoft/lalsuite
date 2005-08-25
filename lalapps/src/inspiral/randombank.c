/*----------------------------------------------------------------------- 
 * 
 * File Name: randombank.c
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
#include <series.h>
#include <processtable.h>

#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/LALInspiral.h>
#include <lal/LALInspiralBank.h>
#include "inspiral.h"

RCSID( "$Id$" );
#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "tmpltbank"

int arg_parse_check( int argc, char *argv[], MetadataTable procparams );

/*
 *
 * variables that control program behaviour
 *
 */

/* debugging */
extern int vrbflg;                      /* verbocity of lal function    */

/* template bank generation parameters */
LIGOTimeGPS gpsStartTime = { 0, 0 };    /* input data GPS start time    */
LIGOTimeGPS gpsEndTime = { 0, 0 };      /* input data GPS end time      */
REAL4   minMass         = -1;           /* minimum component mass       */
REAL4   maxMass         = -1;           /* maximum component mass       */
REAL4  fLow             = -1;           /* low frequency cutoff         */
INT4    numTmplts       = -1;           /* number of template to create */
enum { unset, urandom, user } randSeedType = unset;    /* sim seed type */
INT4  randomSeed        = 0;            /* value of sim rand seed       */

/* output parameters */
CHAR  *userTag          = NULL;

int main ( int argc, char *argv[] )
{
  /* lal function variables */
  LALStatus             status = blank_status;
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;

  /* templates */
  RandomParams         *randParams = NULL;
  InspiralTemplate      newTmplt;
  SnglInspiralTable    *thisTmplt  = NULL;

  /* output data */
  MetadataTable         templateBank;
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         searchsummvars;
  SearchSummvarsTable  *this_search_summvar = NULL;
  ProcessParamsTable   *this_proc_param = NULL;
  LIGOLwXMLStream       results;

  /* counters and other variables */
  INT4 i;
  CHAR  fname[256];


  /*
   *
   * initialization
   *
   */


  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "1" );

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) 
    calloc( 1, sizeof(ProcessTable) );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->start_time),
        &accuracy ), &status );
  LAL_CALL( populate_process_table( &status, proctable.processTable, 
        PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE ), &status );
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
    calloc( 1, sizeof(ProcessParamsTable) );

  /* call the argument parse and check function */
  arg_parse_check( argc, argv, procparams );

  /* can use LALMalloc() / LALCalloc() from here */


  /*
   *
   * create the radom number seed
   *
   */


  /* store the seed in the search summvars table */
  this_search_summvar = searchsummvars.searchSummvarsTable = 
    (SearchSummvarsTable *) LALCalloc( 1, sizeof(SearchSummvarsTable) );
  LALSnprintf( this_search_summvar->name, 
      LIGOMETA_NAME_MAX * sizeof(CHAR), "template bank simulation seed" );

  if ( randSeedType == urandom )
  {
    FILE   *fpRand = NULL;
    INT4    randByte;

    if ( vrbflg ) 
      fprintf( stdout, "obtaining random seed from /dev/urandom: " );

    randomSeed = 0;
    fpRand = fopen( "/dev/urandom", "r" );
    if ( fpRand )
    {
      for ( randByte = 0; randByte < 4 ; ++randByte )
      {
        INT4 tmpSeed = (INT4) fgetc( fpRand );
        randomSeed += tmpSeed << ( randByte * 8 );
      }
      fclose( fpRand );
    }
    else
    {
      perror( "error obtaining random seed from /dev/urandom" );
      exit( 1 );
    }
  }
  else if ( randSeedType == user )
  {
    if ( vrbflg ) 
      fprintf( stdout, "using user specified random seed: " );
  }

  this_search_summvar->value = randomSeed;
  LALSnprintf( this_search_summvar->string, LIGOMETA_STRING_MAX * sizeof(CHAR),
      "%d", randomSeed );
  if ( vrbflg ) fprintf( stdout, "%d\n", randomSeed );

  /* create the tmplt bank random parameter structure */
  LAL_CALL( LALCreateRandomParams( &status, &randParams, randomSeed ),
      &status );


  /*
   *
   * create a random template bank
   *
   */


  /* make sure the pointer to the first template is null */
  templateBank.snglInspiralTable = NULL;

  for ( i = 0; i < numTmplts; ++i )
  {
    memset( &newTmplt, 0, sizeof(InspiralTemplate) );
    newTmplt.massChoice = m1Andm2;
    newTmplt.order = twoPN;
    newTmplt.fLower = fLow;

    /* set up the injection masses */
    if ( maxMass == minMass )
    {
      newTmplt.mass1 = (REAL8) maxMass;
      newTmplt.mass2 = (REAL8) maxMass;
    }
    else
    {
      REAL4 mass;
      /* generate random parameters for the injection */
      LAL_CALL( LALUniformDeviate( &status, &mass, randParams ), &status );
      newTmplt.mass1 = (maxMass - minMass) * mass;
      newTmplt.mass1 += minMass;

      LAL_CALL( LALUniformDeviate( &status, &mass, randParams ), &status );
      newTmplt.mass2 = (maxMass - minMass) * mass;
      newTmplt.mass2 += minMass;
    }

    LAL_CALL( LALInspiralParameterCalc( &status, &newTmplt ), &status );

    if ( ! templateBank.snglInspiralTable )
    {
      thisTmplt = templateBank.snglInspiralTable =
        (SnglInspiralTable *) LALCalloc(1, sizeof(SnglInspiralTable));
    }
    else
    {
      thisTmplt = thisTmplt->next = 
        (SnglInspiralTable *) LALCalloc(1, sizeof(SnglInspiralTable));
    }

    thisTmplt->mass1 = newTmplt.mass1;
    thisTmplt->mass2 = newTmplt.mass2;
    thisTmplt->mchirp = newTmplt.chirpMass;
    thisTmplt->eta = newTmplt.eta;
    thisTmplt->tau0 = newTmplt.t0;
    thisTmplt->tau2 = newTmplt.t2;
    thisTmplt->tau3 = newTmplt.t3;
    thisTmplt->tau4 = newTmplt.t4;
    thisTmplt->tau5 = newTmplt.t5;
    thisTmplt->ttotal = newTmplt.tC;
    thisTmplt->psi0 = newTmplt.psi0;
    thisTmplt->psi3 = newTmplt.psi3;
    thisTmplt->f_final = newTmplt.fFinal;
    thisTmplt->eta = newTmplt.eta;
    thisTmplt->beta = newTmplt.beta;
    LALSnprintf( thisTmplt->ifo, LIGOMETA_IFO_MAX * sizeof(CHAR), "P1" );
    LALSnprintf( thisTmplt->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
        "randombank" );
    LALSnprintf( thisTmplt->channel, LIGOMETA_CHANNEL_MAX * sizeof(CHAR),
        "SIM-BANK" );
  }


  /*
   *
   * write the output data
   *
   */

  /* open the output xml file */
  memset( &results, 0, sizeof(LIGOLwXMLStream) );
  if ( userTag )
  {
    LALSnprintf( fname, sizeof(fname), "P1-TMPLTBANK_%s-%d-%d.xml",
        userTag, gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  else
  {
    LALSnprintf( fname, sizeof(fname), "P1-TMPLTBANK-%d-%d.xml",
        gpsStartTime.gpsSeconds,
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &results, fname), &status );

  /* write the process table */
  LALSnprintf( proctable.processTable->ifos, LIGOMETA_IFO_MAX, "P1" );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->end_time),
        &accuracy ), &status );
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, proctable, 
        process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  free( proctable.processTable );

  /* erase the first empty process params entry */
  {
    ProcessParamsTable *emptyPPtable = procparams.processParamsTable;
    procparams.processParamsTable = procparams.processParamsTable->next;
    free( emptyPPtable );
  }

  /* write the process params table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, process_params_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, procparams, 
        process_params_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  while( procparams.processParamsTable )
  {
    this_proc_param = procparams.processParamsTable;
    procparams.processParamsTable = this_proc_param->next;
    free( this_proc_param );
  }

  /* write the search summvars table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, 
        search_summvars_table ), &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, searchsummvars, 
        search_summvars_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  while( searchsummvars.searchSummvarsTable )
  {
    this_search_summvar = searchsummvars.searchSummvarsTable;
    searchsummvars.searchSummvarsTable = this_search_summvar->next;
    LALFree( this_search_summvar );
  }

  /* write the template bank to the file */
  if ( templateBank.snglInspiralTable )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &results, sngl_inspiral_table ), 
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &results, templateBank, 
          sngl_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &results ), &status );
  }
  while ( templateBank.snglInspiralTable )
  {
    thisTmplt = templateBank.snglInspiralTable;
    templateBank.snglInspiralTable = templateBank.snglInspiralTable->next;
    LALFree( thisTmplt );
  }

  /* close the output xml file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &results ), &status );

  /* free the rest of the memory, check for memory leaks and exit */
  LAL_CALL( LALDestroyRandomParams( &status, &randParams ), &status );
  LALCheckMemoryLeaks();
  exit( 0 );
}


/* ------------------------------------------------------------------------- */

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
  calloc( 1, sizeof(ProcessParamsTable) ); \
  LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
      PROGRAM_NAME ); \
      LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
          long_options[option_index].name ); \
          LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
          LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

#define USAGE \
"  --help                       display this message\n"\
"  --version                    print version information and exit\n"\
"  --debug-level LEVEL          set the LAL debug level to LEVEL\n"\
"  --user-tag STRING            set the process_params usertag to STRING\n"\
"\n"\
"  --gps-start-time SEC         GPS second of data start time\n"\
"  --gps-end-time SEC           GPS second of data end time\n"\
"  --low-frequency-cutoff F     Compute tau parameters from F Hz\n"\
"  --minimum-mass MASS          set minimum component mass of bank to MASS\n"\
"  --maximum-mass MASS          set maximum component mass of bank to MASS\n"\
"  --number-of-template N       create N random templatess of bank to MASS\n"\
"  --random-seed SEED           set random number seed for injections to SEED\n"\
  "                                 (urandom|integer)\n"\
  "\n"

int arg_parse_check( int argc, char *argv[], MetadataTable procparams )
{
  /* getopt arguments */
  struct option long_options[] =
  {
    /* these options set a flag */
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"help",                    no_argument,       0,                'h'},
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"low-frequency-cutoff",    required_argument, 0,                'i'},
    {"minimum-mass",            required_argument, 0,                'A'},
    {"maximum-mass",            required_argument, 0,                'B'},
    {"number-of-templates",     required_argument, 0,                'K'},
    {"random-seed",             required_argument, 0,                'J'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    {"debug-level",             required_argument, 0,                'z'},
    {0, 0, 0, 0}
  };
  int c;
  ProcessParamsTable *this_proc_param = procparams.processParamsTable;


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

    c = getopt_long_only( argc, argv, "hs:a:b:i:A:B:K:J:Z:z:", 
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

      case 'h':
        fprintf( stdout, USAGE );
        exit( 0 );
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

      case 'A':
        minMass = (REAL4) atof( optarg );
        if ( minMass <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "miniumum component mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, minMass );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", minMass );
        break;

      case 'B':
        maxMass = (REAL4) atof( optarg );
        if ( maxMass <= 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "maxiumum component mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, maxMass );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", maxMass );
        break;

      case 'J':
        if ( ! strcmp( "urandom", optarg ) )
        {
          randSeedType = urandom;
          ADD_PROCESS_PARAM( "string", "%s", optarg );
        }
        else
        {
          randSeedType = user;
          randomSeed = (INT4) atoi( optarg );
          ADD_PROCESS_PARAM( "int", "%d", randomSeed );
        }
        break;

      case 'K':
        numTmplts = (INT4) atoi( optarg );
        if ( numTmplts < 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "number of template bank simulations"
              "must be greater than 1: (%d specified)\n", 
              long_options[option_index].name, numTmplts );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "int", "%d", numTmplts );
        break;

      case 'a':
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
          gpsStartTime.gpsSeconds = (INT4) gstartt;
          gpsStartTime.gpsNanoSeconds = 0;
          ADD_PROCESS_PARAM( "int", "%ld", gstartt );
        }
        break;

      case 'b':
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
          gpsEndTime.gpsSeconds = (INT4) gendt;
          gpsEndTime.gpsNanoSeconds = 0;
          ADD_PROCESS_PARAM( "int", "%ld", gendt );
        }
        break;

      case 'i':
        fLow = (REAL4) atof( optarg );
        if ( fLow < 0 )
        {
          fprintf( stdout, "invalid argument to --%s:\n"
              "low frequency cutoff is less than 0 Hz: "
              "(%f Hz specified)\n",
              long_options[option_index].name, fLow );
          exit( 1 );
        }
        ADD_PROCESS_PARAM( "float", "%e", fLow );
        break;

      case 'z':
        set_debug_level( optarg );
        ADD_PROCESS_PARAM( "string", "%s", optarg );
        break;

      case '?':
        fprintf( stderr, USAGE );
        exit( 1 );
        break;

      default:
        fprintf( stderr, "unknown error while parsing options\n" );
        fprintf( stderr, USAGE );
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
   * check validity of arguments
   *
   */


  if ( minMass < 0 )
  {
    fprintf( stderr, "--minimum-mass must be specified\n" );
    exit( 1 );
  }
  if ( maxMass < 0 )
  {
    fprintf( stderr, "--maximum-mass must be specified\n" );
    exit( 1 );
  }

  /* check that the bank parameters have been specified */
  if ( numTmplts < 0 )
  {
    fprintf( stderr, "--number-of-templates must be specified\n" );
    exit( 1 );
  }

  return 0;
}

#undef ADD_PROCESS_PARAM
