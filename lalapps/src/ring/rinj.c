/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Lisa M. Goggin
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
 * File Name: rinj.c
 *
 * Author: Goggin, L. M., and  Brown, D. A.
 *
 * $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <config.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <regex.h>
#include <glob.h>
#include <lalapps.h>
#include <processtable.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXML.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/Ring.h>
#include <lal/lalGitID.h>
#include <lalappsGitID.h>


RCSID( "$Id$" );
#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "rinj"

#define USAGE \
"lalapps_rinj [options]\n"\
"\nDefaults are shown in brackets\n\n" \
"  --help                   display this message\n"\
"  --verbose                turn verbose flag on\n"\
"  --gps-start-time TIME    start injections at GPS time TIME (793130413)\n"\
"  --gps-end-time TIME      end injections at GPS time TIME (795679213)\n"\
"  --time-step STEP         space injections by STEP / pi seconds apart (2630)\n"\
"  --time-interval TIME     distribute injections in interval TIME (250)\n"\
"  --seed SEED              seed random number generator with SEED (1)\n"\
"  --user-tag STRING        set the usertag to STRING\n"\
"  --inj-distr INJDISTR     distribute injections uniformly in\n"\
"                           log_10(frequency) ( INJDISTR = 0), frequency (INJDISTR = 1)\n"\
"                           or mass (INJDISTR = 2) (default INJDISTR = 0)\n"\
"  --minimum-mass MIN       set the minimum componenet mass to MIN (13.8)\n"\
"  --maximum-mass MAX       set the maximum componenet mass to MAX (236.8)\n"\
"  --minimum-spin AMIN      set the minimum component of the dimensionless spin parameter (0)\n"\
"  --maximum-spin AMAX      set the maximum component of the dimensionless spin parameter (0.994)\n"\
"  --minimum-quality MIN    set minimum quality factor to MIN (2)\n"\
"  --maximum-quality MAX    set maximum quality factor to MAX (20)\n"\
"  --minimum-frequency MIN  set minimum frequency to MIN (50)\n"\
"  --maximum-frequency MAX  set maximum frequency to MAX (2000)\n"\
"  --minimum-distance DMIN  set the minimum distance to DMIN kpc (1)\n"\
"  --maximum-distance DMAX      set the maximum distance to DMAX kpc (200000)\n"\
"  --epsilon EPS            amount of energy radiated as gravitational waves (0.01)\n"\
"  --waveform WVF           set the injection waveform to WVF\n"\
"\n"

/* all units are in kpc since this is what GalacticInspiralParamStruc expects */

extern int vrbflg;
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
  if ( ! strcmp( name, "userTag" ) || ! strcmp( name, "user-tag" ) )
    LALSnprintf( pp->param, LIGOMETA_PARAM_MAX, "-userTag" );
  else
    LALSnprintf( pp->param, LIGOMETA_PARAM_MAX, "--%s", name );
  strncpy( pp->type, type, LIGOMETA_TYPE_MAX );
  va_start( ap, fmt );
  vsnprintf( pp->value, LIGOMETA_VALUE_MAX, fmt, ap );
  va_end( ap );
  return pp;
}

int main( int argc, char *argv[] )
{
  LALStatus             status = blank_status;
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;
  const INT4            S4StartTime = 793130413; /* Tues Feb 22, 2005 at 12:00 CST (10:00 PST) */
  const INT4            S4StopTime  = 795679213; /* Wed March 23, 2005 at 24:00 CST (22:00 PST) */

  /* command line options */
  LIGOTimeGPS   gpsStartTime;
  LIGOTimeGPS   gpsEndTime;
  REAL8         meanTimeStep = 7000 / LAL_PI;
  REAL8         timeInterval = 250;
  REAL8         tstep = 2630;
  UINT4         randSeed = 1;
  CHAR         *userTag = NULL;
  REAL4         minMass = 13.8;
  REAL4         maxMass = 236.8;
  REAL4         minSpin = 0;
  REAL4         maxSpin = 0.994;
  REAL4         minFreq = 50.0;
  REAL4         maxFreq = 2000.0;
  REAL4         minQuality = 2;
  REAL4         maxQuality = 20.0;
  REAL4         dmin = 1;
  REAL4         dmax = 200000;
  REAL4         epsilon = 0.010;
  UINT4         injdistr = 0;
  
  /* program variables */
  RandomParams *randParams = NULL;
  REAL4  u, exponent, expt;
  REAL4  deltaM;
  LALMSTUnitsAndAcc     gmstUnits = { MST_HRS, LALLEAPSEC_STRICT };
  LALGPSandAcc          gpsAndAcc;
  SkyPosition           skyPos;
  LALSource             source;
  LALPlaceAndGPS        placeAndGPS;
  DetTimeAndASource     detTimeAndSource;
  LALDetector           lho = lalCachedDetectors[LALDetectorIndexLHODIFF];
  LALDetector           llo = lalCachedDetectors[LALDetectorIndexLLODIFF];
  LALDetAndSource       detAndSource;
  LALDetAMResponse      resp;
  REAL8                 time_diff_ns;
  REAL4                 splus, scross, cosiota;
  
  /* waveform */
  CHAR waveform[LIGOMETA_WAVEFORM_MAX];
  CHAR coordinates[LIGOMETA_COORDINATES_MAX];  
  
  LALGPSCompareResult        compareGPS;

  /*  xml output data */
  CHAR                  fname[256];
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         injections;
  ProcessParamsTable   *this_proc_param;
  SimRingdownTable     *this_inj = NULL; 
  LIGOLwXMLStream       xmlfp;

  /* getopt arguments */
  struct option long_options[] =
  {
    {"help",                    no_argument,       0,                'h'},
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"playground",              no_argument,       0,                'x'},
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"time-step",               required_argument, 0,                't'},
    {"time-interval",           required_argument, 0,                'i'},
    {"seed",                    required_argument, 0,                's'},
    {"inj-distr",               required_argument, 0,                'G'},
    {"minimum-mass",            required_argument, 0,                'A'},
    {"maximum-mass",            required_argument, 0,                'B'},
    {"minimum-spin",            required_argument, 0,                'P'},
    {"maximum-spin",            required_argument, 0,                'Q'},
    {"minimum-frequency",       required_argument, 0,                'C'},
    {"maximum-frequency",       required_argument, 0,                'D'},
    {"minimum-quality",         required_argument, 0,                'E'},
    {"maximum-quality",         required_argument, 0,                'F'},
    {"minimum-distance",        required_argument, 0,                'V'},
    {"maximum-distance",        required_argument, 0,                'W'},
    {"epsilon",                 required_argument, 0,                'r'},
    {"debug-level",             required_argument, 0,                'z'},
    {"waveform",                required_argument, 0,                'w'},
    {"coordinates",             required_argument, 0,                'c'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    {0, 0, 0, 0}
  };
  int c;

  REAL4 lmin;
  REAL4 lmax;
  REAL4 deltaL;
  REAL4 deltaA;
  REAL4 fmin;
  REAL4 fmax;
  REAL4 deltaf;

  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "1" );

  gpsStartTime.gpsSeconds = S4StartTime;
  gpsEndTime.gpsSeconds   = S4StopTime;
  gpsStartTime.gpsNanoSeconds = 0;
  gpsEndTime.gpsNanoSeconds   = 0;

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) 
    calloc( 1, sizeof(ProcessTable) );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->start_time),
        &accuracy ), &status );
  if (strcmp(CVS_REVISION, "$Revi" "sion$"))
  {
    XLALPopulateProcessTable(proctable.processTable, PROGRAM_NAME,
        CVS_REVISION, CVS_SOURCE, CVS_DATE, 0);
  }
  else
  {
    XLALPopulateProcessTable(proctable.processTable, PROGRAM_NAME,
        lalappsGitCommitID, lalappsGitGitStatus, lalappsGitCommitDate, 0);
  }
  LALSnprintf( proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " " );
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
    calloc( 1, sizeof(ProcessParamsTable) );
  

  /* clear the waveform field */
  memset( waveform, 0, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR) );
  memset( coordinates, 0, LIGOMETA_COORDINATES_MAX * sizeof(CHAR) );
    
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
        "a:A:b:B:C:D:E:F:G:h:P:Q:r:s:t:V:W:vz:Z:", long_options, &option_index );

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
        tstep = (REAL8) atof( optarg );
        meanTimeStep = tstep / LAL_PI;
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
        timeInterval = atof( optarg );
        if ( timeInterval < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "time interval must be >= 0: (%le seconds specified)\n",
              long_options[option_index].name, meanTimeStep );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%le", timeInterval );
        break;
     
        case 'G':
        injdistr = (UINT4) atoi( optarg );
        if ( injdistr != 0 && injdistr != 1 && injdistr != 2 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "inj-distr must be either 0, 1 or 2\n",
              long_options[option_index].name);
          exit(1);
        }
          this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
             "int", "%d", injdistr );
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

      case 'P':
        /* minimum spin */
        minSpin = (REAL4) atof( optarg );
        if ( minSpin < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "the minimum spin must be > 0: "
              "(%f specified)\n",
              long_options[option_index].name, minSpin );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%e", minSpin );
        break;

      case 'Q':
        /* maximum spin */
        maxSpin = (REAL4) atof( optarg );
        if ( maxSpin > 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "spin parameter must be < 1: "
              "(%f specified)\n",
              long_options[option_index].name, maxSpin );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%e", maxSpin );
        break;

      case 'C':
        minFreq = (REAL4) atof( optarg );
        if ( minFreq <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "miniumum frequency must be > 0: "
              "(%f Hz specified)\n",
              long_options[option_index].name, minFreq );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%e", minFreq );
        break;

      case 'D':
        maxFreq = (REAL4) atof( optarg );
        if ( maxFreq <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "maxiumum frequency must be > 0: "
              "(%f Hz specified)\n",
              long_options[option_index].name, maxFreq );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%e", maxFreq );
        break;

      case 'E':
        /* minimum quality factor */
        minQuality = (REAL4) atof( optarg );
        if ( minQuality < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "the minimum quality factor must be > 0: "
              "(%f specified)\n",
              long_options[option_index].name, minQuality );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%e", minQuality );
        break;

      case 'F':
        /* maximum quality factor */
        maxQuality = (REAL4) atof( optarg );
        if ( maxQuality < 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "the maximum quality factor must be > 0: "
              "(%f specified)\n",
              long_options[option_index].name, maxQuality );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name,
              "float", "%e", maxQuality );
        break;

       
      case 'V':
        /* minimum distance from earth */
        dmin = (REAL4) atof( optarg );
        if ( dmin <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "minimum distance must be > 0: "
              "(%f kpc specified)\n",
              long_options[option_index].name, dmin );
              exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%e", dmin );
        break;

      case 'W':
        /* max distance from earth */
        dmax = (REAL4) atof( optarg );
        if ( dmax <= 0 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "maximum distance must be greater than 0: "
              "(%f kpc specified)\n",
              long_options[option_index].name, dmax );
              exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
        next_process_param( long_options[option_index].name, 
            "float", "%e", dmax );
        break;
        
      case 'r':
        /* epsilon */
        epsilon = (REAL4) atof( optarg );
        if ( epsilon <= 0 || epsilon > 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "fraction of radiated energy must be between 0 and 1: "
              "(%f specified)\n",
              long_options[option_index].name, epsilon );
          exit( 1 );
        }
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "float", "%e", epsilon );
        break;

      case 'Z':
       /* create storage for the usertag */
        optarg_len = strlen( optarg ) + 1;
        userTag = (CHAR *) calloc( optarg_len, sizeof(CHAR) );
        memcpy( userTag, optarg, optarg_len );
        this_proc_param = this_proc_param->next = 
          next_process_param( long_options[option_index].name, 
              "string", "%s", optarg );
        break;
      
      case 'w':
        LALSnprintf( waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "%s",
            optarg);
        this_proc_param = this_proc_param->next =
           next_process_param( long_options[option_index].name, "string",
              "%s", optarg);
        break;
      
      case 'c':
        LALSnprintf( coordinates, LIGOMETA_COORDINATES_MAX * sizeof(CHAR), "%s",
            optarg);
        this_proc_param = this_proc_param->next =
          next_process_param( long_options[option_index].name, "string",
              "%s", optarg);
        break;
                                
      case 'v':
        vrbflg = 1;
        break;
      
      case 'z':
        set_debug_level( optarg );
        this_proc_param = this_proc_param->next = 
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

  if ( !*waveform )
    {
      /* use Ringdown as the default waveform */
      LALSnprintf( waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR),
          "Ringdown");
      }
  
  if ( !*coordinates )
        {
          /* use equatorial as the default system */
          LALSnprintf( coordinates, LIGOMETA_COORDINATES_MAX * sizeof(CHAR),
                              "EQUATORIAL");
                      }
  
                
  /*
   *
   * initialization
   *
   */


  /* initialize the random number generator */
  LAL_CALL( LALCreateRandomParams( &status, &randParams, randSeed ), &status );

  /* distance range */
   lmin = log10(dmin);
   lmax = log10(dmax);
   deltaL = lmax - lmin;
           
  /* spin range */
   deltaA = maxSpin - minSpin;
  
  /* null out the head of the linked list */
   injections.simRingdownTable = NULL;

  /* create the output file name */
  if ( userTag )
  {
    LALSnprintf( fname, sizeof(fname), "HL-INJECTIONS_%d_%s-%d-%d.xml", 
        randSeed, userTag, gpsStartTime.gpsSeconds, 
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }
  else
  {
    LALSnprintf( fname, sizeof(fname), "HL-INJECTIONS_%d-%d-%d.xml", 
        randSeed, gpsStartTime.gpsSeconds, 
        gpsEndTime.gpsSeconds - gpsStartTime.gpsSeconds );
  }

  /* check that the start time is before the end time */
  LAL_CALL( LALCompareGPS( &status, &compareGPS, &gpsStartTime, &gpsEndTime ),
      &status );
  
  /*
   *
   * loop over duration of desired output times
   *
   */

  while ( compareGPS == LALGPS_EARLIER )
  {
    /* create the sim_ringdown table */
    if ( injections.simRingdownTable )
      {
        this_inj = this_inj->next = (SimRingdownTable *)
          LALCalloc( 1, sizeof(SimRingdownTable) );
        }
    else
      {
        injections.simRingdownTable = this_inj = (SimRingdownTable *)
          LALCalloc( 1, sizeof(SimRingdownTable) );
        }

    /* set the waveform and coordinates fields */
    memcpy( this_inj->waveform, waveform, LIGOMETA_WAVEFORM_MAX *
        sizeof(CHAR));
    memcpy( this_inj->coordinates, coordinates, LIGOMETA_COORDINATES_MAX *
                sizeof(CHAR));
    
    this_inj->epsilon = epsilon;
    
    /* set the geocentric start time of the injection */
    this_inj->geocent_start_time = gpsStartTime;
    if ( timeInterval )
    {
      LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
      LAL_CALL( LALAddFloatToGPS( &status, &(this_inj->geocent_start_time),
          &(this_inj->geocent_start_time), u * timeInterval ), &status );
    }    

 
    if ( injdistr == 0 )
    /* uniform in log frequency */
    {
      /* set frequency, f0, and quality factor Q */
      fmin = log10(minFreq);
      fmax = log10(maxFreq);
      deltaf = fmax - fmin;
      LAL_CALL(  LALUniformDeviate(&status,&u,randParams),&status );
      expt = fmin + deltaf * u;
      this_inj->frequency = pow(10.0,(REAL4) expt);
    }    
    else if ( injdistr == 1)
    /* uniform in frequency */
    {
      LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
      this_inj->frequency = minFreq + u * (maxFreq - minFreq);
    }

    if ( injdistr == 0 || injdistr == 1 )
    {
      LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
      this_inj->quality = minQuality + u * (maxQuality - minQuality);
      /* calculate M and a from f and Q */
      this_inj->spin = XLALBlackHoleRingSpin( this_inj->quality);
      this_inj->mass = XLALBlackHoleRingMass( this_inj->frequency, this_inj->quality);
    }

    if ( injdistr == 2 )
    {
      /* mass distribution */
      deltaM = maxMass - minMass;
      LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
      this_inj->mass = minMass + u * deltaM;
      /* generate random spin parameter */  
      LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
      this_inj->spin = minSpin + u * deltaA;
      /* calculate central frequency, f0, and quality factor Q */
      this_inj->frequency = XLALBlackHoleRingFrequency( this_inj->mass, this_inj->spin );
      this_inj->quality = XLALBlackHoleRingQuality( this_inj->spin );
    }


            
    /* spatial distribution */
    /* compute random longitude and latitude */ 
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    this_inj->longitude = LAL_TWOPI * u ;
    
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    this_inj->latitude = asin( 2.0 * u - 1.0 ) ;
   
    /* initial phase */
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    this_inj->phase = u * LAL_TWOPI;
    
    /* uniform distribution in log(distance) */
    LAL_CALL(  LALUniformDeviate(&status,&u,randParams),&status );
    exponent = lmin + deltaL * u;
    this_inj->distance = pow(10.0,(REAL4) exponent); 
    this_inj->distance = this_inj->distance / 1000.0; /*convert to Mpc */

    /* compute random inclination, polarization */
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    this_inj->inclination = acos( 2.0 * u - 1.0 );
    
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    this_inj->polarization = LAL_TWOPI * u ;
    
    gpsAndAcc.gps = this_inj->geocent_start_time;

    /* set gmst */
    LAL_CALL( LALGPStoGMST1( &status, &(this_inj->start_time_gmst),
          &(this_inj->geocent_start_time), &gmstUnits ), &status);
    
    memset( &skyPos, 0, sizeof(SkyPosition) );
    memset( &source, 0, sizeof(LALSource) );
    memset( &placeAndGPS, 0, sizeof(LALPlaceAndGPS) );
    memset( &detTimeAndSource, 0, sizeof(DetTimeAndASource) );
    memset( &detAndSource, 0, sizeof(LALDetAndSource) );

    skyPos.longitude = this_inj->longitude;
    skyPos.latitude  = this_inj->latitude;
    skyPos.system    = COORDINATESYSTEM_EQUATORIAL;

    source.equatorialCoords = skyPos;
    source.orientation      = this_inj->polarization;
    
    placeAndGPS.p_gps = &(this_inj->geocent_start_time);
    
    detTimeAndSource.p_det_and_time = &placeAndGPS;
    detTimeAndSource.p_source = &skyPos;
    detAndSource.pSource = &source;
                    
    gpsAndAcc.accuracy = LALLEAPSEC_STRICT;
    gpsAndAcc.gps = this_inj->geocent_start_time;
    
    /* calculate h0 */
    this_inj->amplitude = XLALBlackHoleRingAmplitude( this_inj->frequency,
        this_inj->quality, this_inj->distance, this_inj->epsilon );
      
    /* calculate hrss */
    this_inj->hrss = this_inj->amplitude * sqrt( 2 / LAL_PI / this_inj->frequency ) * 
      pow( ( 2.0 * pow( this_inj->quality, 3.0 ) + this_inj->quality ) / 
          ( 1.0 + 4.0 * pow ( this_inj->quality, 2 ) ) , 0.5);
      
    /* initialize end times with geocentric value */
    this_inj->h_start_time = this_inj->l_start_time = this_inj->geocent_start_time;
    
    /* initialize distances with real distance and compute splus and scross*/
    this_inj->eff_dist_h = this_inj->eff_dist_l = 2.0 * this_inj->distance;
    cosiota = cos( this_inj->inclination );
    splus = -( 1.0 + cosiota * cosiota );
    scross = -2.0 * cosiota; 
      
    /* lho */
    placeAndGPS.p_detector = &lho;
    LAL_CALL( LALTimeDelayFromEarthCenter( &status, &time_diff_ns,
          &detTimeAndSource ), &status );
    LAL_CALL( LALAddFloatToGPS( &status, &(this_inj->h_start_time),
          &(this_inj->h_start_time), time_diff_ns ), &status );

    /* compute the response of the LHO detectors */
    detAndSource.pDetector = &lho;
    LAL_CALL( LALComputeDetAMResponse( &status, &resp, &detAndSource,
          &gpsAndAcc ), &status );
    
    /* compute the effective distance for LHO */
    this_inj->eff_dist_h /= sqrt( splus*splus*resp.plus*resp.plus +
        scross*scross*resp.cross*resp.cross );

    /* compute hrss at LHO */ 
    this_inj->hrss_h = this_inj->amplitude * pow ( ( 
          (2*pow(this_inj->quality,3)+this_inj->quality ) * splus*splus*resp.plus*resp.plus +
          2*pow(this_inj->quality,2) * splus*scross*resp.plus*resp.cross +
          2*pow(this_inj->quality,3) * scross*scross*resp.cross*resp.cross )
        /  2.0 / LAL_PI / this_inj->frequency / ( 1.0 + 4.0 * pow ( this_inj->quality, 2 ) ) , 0.5 );
      
    /* llo */
    placeAndGPS.p_detector = &llo;
    LAL_CALL( LALTimeDelayFromEarthCenter( &status,  &time_diff_ns,
          &detTimeAndSource ), &status);
    LAL_CALL( LALAddFloatToGPS( &status,  &(this_inj->l_start_time),
          &(this_inj->l_start_time), time_diff_ns ), &status);

    /* compute the response of the LLO detector */
    detAndSource.pDetector = &llo;
    LAL_CALL( LALComputeDetAMResponse( &status, &resp, &detAndSource,
          &gpsAndAcc ), &status);
    
    /* compute the effective distance for LLO */
    this_inj->eff_dist_l /= sqrt( splus*splus*resp.plus*resp.plus 
        + scross*scross*resp.cross*resp.cross );
    
    /* compute hrss at LLO */
    this_inj->hrss_l = this_inj->amplitude * pow ( (
          (2*pow(this_inj->quality,3)+this_inj->quality ) * splus*splus*resp.plus*resp.plus +
          2*pow(this_inj->quality,2) * splus*scross*resp.plus*resp.cross +
          2*pow(this_inj->quality,3) * scross*scross*resp.cross*resp.cross )
          /  2.0 / LAL_PI / this_inj->frequency / ( 1.0 + 4.0 * pow ( this_inj->quality, 2 ) ) , 0.5 );
        
    /* increment the injection time */
    LAL_CALL( LALAddFloatToGPS( &status, &gpsStartTime, &gpsStartTime, 
          meanTimeStep ), &status );
    LAL_CALL( LALCompareGPS( &status, &compareGPS, &gpsStartTime, 
          &gpsEndTime ), &status );


  } /* end loop over injection times */

  
  /* destroy random parameters */
  LAL_CALL( LALDestroyRandomParams( &status, &randParams ), &status );
  
  


  
  
  /*
   *
   * write output to LIGO_LW XML file
   *
   */


  /* open the xml file */
  memset( &xmlfp, 0, sizeof(LIGOLwXMLStream) );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlfp, fname), &status );

  /* write the process table */
  LALSnprintf( proctable.processTable->ifos, LIGOMETA_IFOS_MAX, "H1H2L1" );
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

  /* write the sim_ringdown table */
  if ( injections.simRingdownTable )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, sim_ringdown_table ), 
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, injections, 
          sim_ringdown_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );
  }
  while ( injections.simRingdownTable )
  {
    this_inj = injections.simRingdownTable;
    injections.simRingdownTable = injections.simRingdownTable->next;
    LALFree( this_inj );
  }
  /* close the injection file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlfp ), &status );

  /* check for memory leaks and exit */
  LALCheckMemoryLeaks();
  return 0;

}


