/*----------------------------------------------------------------------- 
 * 
 * File Name: rinj.c
 *
 * Author: Goggin, L. M. based on minj.c by Brown, D. A. and bbhinj.c by Brown,
 * D. A., Messaritaki E
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

/* ??? */
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
"  --verbose                print mass and galactocentic cartesian coordinates\n"\
"  --gps-start-time TIME    start injections at GPS time TIME (793130413)\n"\
"  --gps-end-time TIME      end injections at GPS time TIME (795679213)\n"\
"  --time-step STEP         space injections by ave of STEP sec (2630/PI)\n"\
"  --seed SEED              seed random number generator with SEED (1)\n"\
"  --user-tag STRING        set the usertag to STRING\n"\
"  --minimum-mass MIN       set the minimum componenet mass to MIN (5)\n"\
"  --maximum-mass MAX       set the maximum componenet mass to MAX (560)\n"\
"  --minimum-spin AMIN      set the minimum component of the dimensionless spin parameter (0)\n"\
"  --maximum-spin AMAX      set the maximum component of the dimensionless spin parameter (0.994)\n"\
"  --minimum-distance DMIN  set the minimum distance to DMIN kpc (1)\n"\
"  --max-distance DMAX      set the maximum distance to DMAX kpc (200000)\n"\
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
  LALVsnprintf( pp->value, LIGOMETA_VALUE_MAX, fmt, ap );
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
  LIGOTimeGPS   gpsStartTime = {S4StartTime, 0};
  LIGOTimeGPS   gpsEndTime   = {S4StopTime, 0};
  REAL8         meanTimeStep = 2630 / LAL_PI;
  UINT4         randSeed = 1;
  CHAR         *userTag = NULL;
  REAL4         minMass = 13.8;
  REAL4         maxMass = 236.8;
  REAL4         minSpin = 0;
  REAL4         maxSpin = 0.994;
  REAL4         dmin = 1;
  REAL4         dmax = 200000;
  REAL4         epsilon = 0.01;
  
  /* program variables */
  RandomParams *randParams = NULL;
  REAL4  u, exponent;
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
  SimRingTable         *this_inj = NULL; 
  LIGOLwXMLStream       xmlfp;

  /* getopt arguments */
  struct option long_options[] =
  {
    {"help",                    no_argument,       0,                'h'},
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"time-step",               required_argument, 0,                't'},
    {"seed",                    required_argument, 0,                's'},
    {"minimum-mass",            required_argument, 0,                'A'},
    {"maximum-mass",            required_argument, 0,                'B'},
    {"minimum-spin",            required_argument, 0,                'P'},
    {"maximum-spin",            required_argument, 0,                'Q'},
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
        "a:A:b:B:h:P:Q:r:s:t:V:W:vz:Z:", long_options, &option_index );

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

  /* mass range */
  deltaM = maxMass - minMass;

  /* distance range */
   REAL4 lmin = log10(dmin);
   REAL4 lmax = log10(dmax);
   REAL4 deltaL = lmax - lmin;
           
  /* spin range */
   REAL4 deltaA = maxSpin - minSpin;
  
  /* null out the head of the linked list */
   injections.simRingTable = NULL;

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
    /* create the sim_ring table */
    if ( injections.simRingTable )
      {
        this_inj = this_inj->next = (SimRingTable *)
          LALCalloc( 1, sizeof(SimRingTable) );
        }
    else
      {
        injections.simRingTable = this_inj = (SimRingTable *)
          LALCalloc( 1, sizeof(SimRingTable) );
        }

    /* set the waveform and coordinates fields */
    memcpy( this_inj->waveform, waveform, LIGOMETA_WAVEFORM_MAX *
        sizeof(CHAR));
    memcpy( this_inj->coordinates, coordinates, LIGOMETA_COORDINATES_MAX *
                sizeof(CHAR));
    
    this_inj->epsilon = epsilon;
    /* set the geocentric start time of the injection */
    /* XXX CHECK XXX */
    this_inj->geocent_start_time = gpsStartTime;

    /* mass distribution */
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    this_inj->totalmass = minMass + u * deltaM;
    
    /* generate random spin parameter */  
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    this_inj->spin = minSpin + u * deltaA;

    /* calculate central frequency, f0, and quality factor Q */
    this_inj->centralfreq = 32000.0 * ( 1.0 - 0.63 * pow( ( 1.0 - this_inj->spin ), 0.3 ) ) / this_inj->totalmass;
    this_inj->quality = 2.0 * pow( ( 1.0 - this_inj->spin ), -0.45 );
            
    /* spatial distribution */
    /* compute random right ascension and declination */ 
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    this_inj->rightascension = LAL_TWOPI * u ;
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    this_inj->declination = asin( 2.0 * u - 1.0 ) ;
    
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

    skyPos.longitude = this_inj->rightascension;
    skyPos.latitude  = this_inj->declination;
    skyPos.system    = COORDINATESYSTEM_EQUATORIAL;

    source.equatorialCoords = skyPos;
    source.orientation      = this_inj->polarization;
    
    placeAndGPS.p_gps = &(this_inj->geocent_start_time);
    
    detTimeAndSource.p_det_and_time = &placeAndGPS;
    detTimeAndSource.p_source = &skyPos;
    detAndSource.pSource = &source;
                    
    gpsAndAcc.accuracy = LALLEAPSEC_STRICT;
    gpsAndAcc.gps = this_inj->geocent_start_time;
    
    /* calculate hpeak */
    this_inj->hpeak = sqrt( 5.0 / 2.0 * this_inj->epsilon )  
      * ( LAL_G_SI * this_inj->totalmass * LAL_MSUN_SI / pow( LAL_C_SI, 2) /
           this_inj->distance / LAL_PC_SI /1000000.0 ) 
       * pow( this_inj->quality, -0.5 ) * pow( 1.0 + 7.0 / 
           24.0 / pow( this_inj->quality, 2.0), -0.5 )
       * pow(  1.0 - 0.63 * pow( 1.0 - this_inj->spin,0.3 ), -0.5);
    
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

    /* initialize 'effective hpeaks' with real hpeak */
    this_inj->hpeak_h = this_inj->hpeak_l = this_inj->hpeak / 2.0;
    
    /* compute hpeak at LHO */ 
    this_inj->hpeak_h *= sqrt(splus*splus*resp.plus*resp.plus + 
        scross*scross*resp.cross*resp.cross);
    
    
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
    
    /* compute hpeak at LLO */
    this_inj->hpeak_l *= sqrt(splus*splus*resp.plus*resp.plus + 
        scross*scross*resp.cross*resp.cross); 

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

  /* write the sim_ring table */
  if ( injections.simRingTable )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, sim_ring_table ), 
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, injections, 
          sim_ring_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );
  }
  while ( injections.simRingTable )
  {
    this_inj = injections.simRingTable;
    injections.simRingTable = injections.simRingTable->next;
    LALFree( this_inj );
  }
  /* close the injection file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlfp ), &status );

  /* check for memory leaks and exit */
  LALCheckMemoryLeaks();
  return 0;

}
