/*----------------------------------------------------------------------- 
 * 
 * File Name: minj.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <config.h>
#include <lalapps.h>
#include <processtable.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXML.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>

RCSID( "$Id$" );
#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "minj"

#define USAGE \
"lalapps_inspinj [options]\n"\
  "\nDefaults are shown in brackets\n\n" \
  "  --help                   display this message\n"\
  "  --gps-start-time TIME    start injections at GPS time TIME (729273613)\n"\
  "  --gps-end-time TIME      end injections at GPS time TIME (734367613)\n"\
  "  --time-step STEP         space injections by ave of STEP sec (2630/PI)\n"\
  "  --time-interval TIME     distribute injections in interval TIME (0)\n"\
  "  --seed SEED              seed random number generator with SEED (1)\n"\
  "  --user-tag STRING        set the usertag to STRING\n"\
  "  --minimum-mass MIN       set the minimum componenet mass to MIN (0.2)\n"\
  "  --maximum-mass MAX       set the maximum componenet mass to MAX (1.0)\n"\
  "  --core-radius A          set the galactic core radius to A kpc (8.5)\n"\
  "  --flatten-halo Q         set the halo flattening parameter to Q (1)\n"\
  "  --halo-radius RMAX       set the maximum halo radius to RMAX kpc (50)\n"\
  "\n"

  ProcessParamsTable *next_process_param( const char *name, const char *type,
      const char *fmt, ... )
{
  ProcessParamsTable *pp;
  va_list ap;
  pp = calloc( 1, sizeof( *pp ) );
  if ( ! pp )
  {
    perror( "next_process_param" );
    exit( 1 );
  }
  strncpy( pp->program, PROGRAM_NAME, LIGOMETA_PROGRAM_MAX );
  LALSnprintf( pp->param, LIGOMETA_PARAM_MAX, "--%s", name );
  strncpy( pp->type, type, LIGOMETA_TYPE_MAX );
  va_start( ap, fmt );
  LALVsnprintf( pp->value, LIGOMETA_VALUE_MAX, fmt, ap );
  va_end( ap );
  return pp;
}

#define KPC ( 1e3 * LAL_PC_SI )
#define MPC ( 1e6 * LAL_PC_SI )
#define GPC ( 1e9 * LAL_PC_SI )

extern int vrbflg;

int main( int argc, char *argv[] )
{
  LALStatus             status = blank_status;
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;
  const INT4            S2StartTime = 729273613; /* Feb 14 2003 16:00:00 UTC */
  const INT4            S2StopTime  = 734367613; /* Apr 14 2003 15:00:00 UTC */

  /* command line options */
  LIGOTimeGPS   gpsStartTime = {S2StartTime, 0};
  LIGOTimeGPS   gpsEndTime   = {S2StopTime, 0};
  REAL8         meanTimeStep = 2630 / LAL_PI;
  REAL8         timeInterval = 0;
  UINT4         randSeed = 1;
  CHAR         *userTag = NULL;
  REAL4         minMass = 0.1;
  REAL4         maxMass = 1.0;
  REAL4         a = 8.5 * KPC;          /* core radius */
  REAL4         Rmax = 50.0 * KPC;      /* halo radius */
  REAL4         q = 1.0;                /* flatten halo */

#if 0
  /* parameters of the injection */
  GalacticInspiralParamStruc    galacticPar;
  PPNParamStruc                 injection;

  /* site end time and effective distance */
  LALPlaceAndGPS       *place_and_gps;
  LALDetector           lho = lalCachedDetectors[LALDetectorIndexLHODIFF];
  LALDetector           llo = lalCachedDetectors[LALDetectorIndexLLODIFF];
  SkyPosition	       *sky_pos;
  DetTimeAndASource    *det_time_and_source;
  REAL8			time_diff;
  REAL8                 site_time;
#endif

  /* xml output data */
  CHAR                  fname[256];
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         injections;
  ProcessParamsTable   *this_proc_param;
  SimInspiralTable     *this_inj;
  LIGOLwXMLStream       xmlfp;

  /* getopt arguments */
  struct option long_options[] =
  {
    {"help",                    no_argument,       0,                'h'},
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"time-step",               required_argument, 0,                't'},
    {"time-interval",		required_argument, 0,		     'i'},
    {"seed",                    required_argument, 0,                's'},
    {"minimum-mass",            required_argument, 0,                'A'},
    {"maximum-mass",            required_argument, 0,                'B'},
    {"core-radius",             required_argument, 0,                'p'},
    {"flatten-halo",            required_argument, 0,                'q'},
    {"halo-radius",             required_argument, 0,                'r'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    {0, 0, 0, 0}
  };
  int c;

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
  LALSnprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
    calloc( 1, sizeof(ProcessParamsTable) );


  /*
   *
   * parse command line arguments
   *
   */

     
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    long int gpsinput;
    size_t optarg_len;

    c = getopt_long_only( argc, argv, 
        "ha:b:t:i:s:A:B:p:q:r:Z:", long_options, &option_index );

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

      case 'a':
        gpsinput = atol( optarg );
        if ( gpsinput < 441417609 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is prior to " 
              "Jan 01, 1994  00:00:00 UTC:\n"
              "(%ld specified)\n",
              long_options[option_index].name, gpsinput );
          exit( 1 );
        }
        if ( gpsinput > 999999999 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is after " 
              "Sep 14, 2011  01:46:26 UTC:\n"
              "(%ld specified)\n", 
              long_options[option_index].name, gpsinput );
          exit( 1 );
        }
        gpsStartTime.gpsSeconds = gpsinput;

        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "int", 
              "%ld", gpsinput );
        break;

      case 'b':
        gpsinput = atol( optarg );
        if ( gpsinput < 441417609 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is prior to " 
              "Jan 01, 1994  00:00:00 UTC:\n"
              "(%ld specified)\n",
              long_options[option_index].name, gpsinput );
          exit( 1 );
        }
        if ( gpsinput > 999999999 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "GPS start time is after " 
              "Sep 14, 2011  01:46:26 UTC:\n"
              "(%ld specified)\n", 
              long_options[option_index].name, gpsinput );
          exit( 1 );
        }
        gpsEndTime.gpsSeconds = gpsinput;
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "int", 
              "%ld", gpsinput );
        break;

      case 's':
        randSeed = atoi( optarg );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "int", 
              "%d", randSeed );
        break;

      case 't':
        meanTimeStep = (REAL8) atof( optarg );
        if ( meanTimeStep <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "time step must be > 0: (%le seconds specified)\n",
              long_options[option_index].name, meanTimeStep );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, "float", 
              "%le", meanTimeStep );
        break;

      case 'i':
        timeInterval = 1000000000LL * atof( optarg );
        if ( timeInterval < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "time interval must be >= 0: (%le seconds specified)\n",
              long_options[option_index].name, meanTimeStep );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%le", timeInterval/1000000000 );
        break;

      case 'A':
        minMass = (REAL4) atof( optarg );
        if ( minMass <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "miniumum component mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, minMass );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%e", minMass );
        break;

      case 'B':
        maxMass = (REAL4) atof( optarg );
        if ( maxMass <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "maxiumum component mass must be > 0: "
              "(%f solar masses specified)\n",
              long_options[option_index].name, maxMass );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%e", maxMass );
        break;

      case 'p':
        /* core-radius */
        a = (REAL4) atof( optarg );
        if ( a <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "galactic core radius must be > 0: "
              "(%f kpc specified)\n",
              long_options[option_index].name, a );
          exit( 1 );
        }
        next_process_param( long_options[option_index].name, 
            "float", "%e", a );
        break;

      case 'q':
        /* flatten-halo */
        q = (REAL4) atof( optarg );
        if ( q <= 0 || q > 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "halo flattening parameter must be in range (0,1]: "
              "(%f specified)\n",
              long_options[option_index].name, q );
          exit( 1 );
        }
        next_process_param( long_options[option_index].name, 
            "float", "%e", q );
        break;

      case 'r':
        /* max halo radius */
        Rmax = (REAL4) atof( optarg );
        if ( Rmax <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "halo radius must be greater than 0: "
              "(%f kpc specified)\n",
              long_options[option_index].name, Rmax );
          exit( 1 );
        }
        next_process_param( long_options[option_index].name, 
            "float", "%e", Rmax );
        break;

      case 'Z':
        /* create storage for the usertag */
        optarg_len = strlen( optarg ) + 1;
        userTag = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( userTag, optarg, optarg_len );
        next_process_param( "-userTag", "string", "%s", optarg );
        break;

      case 'z':
        set_debug_level( optarg );
        next_process_param( long_options[option_index].name, 
            "string", "%s", optarg );
        break;

      case 'h':
        fprintf( stderr, USAGE );
        exit( 0 );
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
   * write output to LIGO_LW XML file
   *
   */


  /* create the file name */
  if ( userTag )
  {
    LALSnprintf( fname, sizeof(fname), "HL-INJECTIONS_%d_%s-%d-%d.xml", 
        randSeed, userTag, gpsStartTime, 
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  else
  {
    LALSnprintf( fname, sizeof(fname), "HL-INJECTIONS_%d-%d-%d.xml", 
        randSeed, gpsStartTime, 
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }

  /* open the xml file */
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlfp, fname), &status );

  /* write the process table */
  LALSnprintf( proctable.processTable->ifos, LIGOMETA_IFO_MAX, "H1H2L1" );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->end_time),
        &accuracy ), &status );
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, process_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, proctable, 
        process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );
  free( proctable.processTable );

  /* free the unused process param entry */
  this_proc_param = procparams.processParamsTable;
  procparams.processParamsTable = procparams.processParamsTable->next;
  free( this_proc_param );

  /* write the process params table */
  if ( procparams.processParamsTable )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, process_params_table ), 
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, procparams, 
          process_params_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );
    while( procparams.processParamsTable )
    {
      this_proc_param = procparams.processParamsTable;
      procparams.processParamsTable = this_proc_param->next;
      free( this_proc_param );
    }
  }

  /* write the sim_inspiral table */
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, sim_inspiral_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, injections, 
        sim_inspiral_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );
  while ( injections.simInspiralTable )
  {
    this_inj= injections.simInspiralTable;
    injections.simInspiralTable = injections.simInspiralTable->next;
    free( this_inj );
  }

  /* close the injection file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlfp ), &status );

  /* check for memory leaks and exit */
  LALCheckMemoryLeaks();
  return 0;
}


