/*----------------------------------------------------------------------- 
 * 
 * File Name: bbhinj.c
 *
 * Author: Brown, D. A., Messaritaki E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <config.h>

#if !defined HAVE_GSL_GSL_FFT_REAL_H || !defined HAVE_LIBGSL
int main( void ) { fprintf( stderr, "no gsl: disabled\n" ); return 77; }
#else

#include <math.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <lalapps.h>
#include <processtable.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXML.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/GenerateInspiral.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>

RCSID( "$Id$" );
#define CVS_ID_STRING "$Id$"
#define CVS_NAME_STRING "$Name$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "bbhinj"

#define USAGE \
"lalapps_bbhinj [options]\n"\
"\nDefaults are shown in brackets\n\n" \
"  --help                   display this message\n"\
"  --verbose                print mass and galactocentic cartesian coordinates\n"\
"  --gps-start-time TIME    start injections at GPS time TIME (729273613)\n"\
"  --gps-end-time TIME      end injections at GPS time TIME (734367613)\n"\
"  --time-step STEP         space injections by ave of STEP sec (2630/PI)\n"\
"  --time-interval TIME     distribute injections in interval TIME (0)\n"\
"  --seed SEED              seed random number generator with SEED (1)\n"\
"  --user-tag STRING        set the usertag to STRING\n"\
"  --min-mass MIN           set the minimum component mass to MIN (3.0)\n"\
"  --max-mass MAX           set the maximum component mass to MAX (20.0)\n"\
"  --min-distance DMIN      set the minimum distance to DMIN kpc (1)\n"\
"  --max-distance DMAX      set the maximum distance to DMAX kpc (20000)\n"\
"  --d-distr DDISTR         distribute injections uniformly over\n"\
"                           d (DDISTR = 0), or over log10(d) (DDISTR = 1)\n"\
"                           or over volume (DDISTR = 2)\n"\
"                           (default: DDISTR = 0)\n"\
"  --m-distr MDISTR         distribute injections uniformly over\n"\
"                           total mass (MDISTR = 0), or over mass1 and\n"\
"                           over mass2 (MDISTR = 1) (default: MDISTR=0)\n"\
"  --waveform WVF           set the injection waveform to WVF\n"\
"                           (EOB, TaylorT1, TaylorT3,PadeT1; followed by the\n"\
"                           order: newtonian, onePN, onePointFivePN, twoPN,\n"\
"                           twoPointFivePN, threePN) (default: EOBtwoPN)\n"\
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
  const INT4            S2StartTime = 729273613; /* Feb 14 2003 16:00:00 UTC */
  const INT4            S2StopTime  = 734367613; /* Apr 14 2003 15:00:00 UTC */

  /* command line options */
  LIGOTimeGPS   gpsStartTime = {S2StartTime, 0};
  LIGOTimeGPS   gpsEndTime   = {S2StopTime, 0};
  REAL8         meanTimeStep = 2630 / LAL_PI;
  REAL8         timeInterval = 0;
  UINT4         randSeed = 1;
  CHAR         *userTag = NULL;
  REAL4         minMass = 3.0;       /* minimum component mass */
  REAL4         maxMass = 20.0;      /* maximum component mass */
  REAL4         dmin = 1.0;          /* minimum distance from earth (kpc) */
  REAL4         dmax = 20000.0 ;     /* maximum distance from earth (kpc) */
 /* REAL4         Rcore = 0.0; */
  UINT4         ddistr = 0, mdistr=0;

  /* program variables */
  RandomParams *randParams = NULL;
  REAL4  u, exponent, d2;
  REAL4  deltaM, mtotal;
  /* XXX CHECK XXX */
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
  /* XXX END CHECK XXX */


  /* waveform */
  CHAR waveform[LIGOMETA_WAVEFORM_MAX];

#if 0
  int i, stat;

  double d, cosphi, sinphi;
#endif

/*  GalacticInspiralParamStruc galacticPar; */
  LALGPSCompareResult        compareGPS;

  /* xml output data */
  CHAR                  fname[256];
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         injections;
  ProcessParamsTable   *this_proc_param;
  SimInspiralTable     *this_inj = NULL;
  LIGOLwXMLStream       xmlfp;

  /* getopt arguments */
  struct option long_options[] =
  {
    {"help",                    no_argument,       0,                'h'},
    {"verbose",                 no_argument,       &vrbflg,           1 },
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"time-step",               required_argument, 0,                't'},
    {"time-interval",		required_argument, 0,		     'i'},
    {"seed",                    required_argument, 0,                's'},
    {"min-mass",                required_argument, 0,                'A'},
    {"max-mass",                required_argument, 0,                'B'},
    {"min-distance",            required_argument, 0,                'p'},
    {"max-distance",            required_argument, 0,                'r'},
    {"d-distr",                 required_argument, 0,                'd'},
    {"m-distr",                 required_argument, 0,                'm'},
    {"waveform",                required_argument, 0,                'w'},
    {"debug-level",             required_argument, 0,                'z'},
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
        "a:A:b:B:hi:p:q:r:s:t:vz:Z:w:", long_options, &option_index );

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

      case 'A':
        /* minimum component mass */
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
        /* maximum component mass */
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

      case 'r':
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

      case 'd':
        ddistr = (UINT4) atof( optarg );
        if ( ddistr != 0 && ddistr != 1 && ddistr != 2)
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "DDISTR must be either 0 or 1 or 2\n",
              long_options[option_index].name);
          exit(1);
        }
        break;

      case 'm':
        mdistr = (UINT4) atof( optarg );
        if ( mdistr != 0 && mdistr != 1 )
        {
          fprintf( stderr, "invalid argument to --%s:\n"
              "MDISTR must be either 0 or 1\n",
              long_options[option_index].name);
          exit(1);
        }
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
              "%le", optarg);

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
    /* use EOBtwoPN as the default waveform */
    LALSnprintf( waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR),
        "EOBtwoPN");
  }


  /*
   *
   * initialization
   *
   */


  /* initialize the random number generator */
  LAL_CALL( LALCreateRandomParams( &status, &randParams, randSeed ), &status );

  /* mass range, per component */
  deltaM = maxMass - minMass;

  /* null out the head of the linked list */
  injections.simInspiralTable = NULL;

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

    /* rho, z and lGal are the galactocentric galactic axial coordinates */
    /* r and phi are the geocentric galactic spherical coordinates       */
#if 0
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    galacticPar.lGal = LAL_TWOPI * u;
    
    galacticPar.z = r * sinphi ;
    galacticPar.rho = 0.0 - Rcore * cos(galacticPar.lGal) +
      sqrt( r * r - galacticPar.z * galacticPar.z - 
      Rcore * Rcore * sin(galacticPar.lGal) * sin(galacticPar.lGal) );
#endif

#if 0
    if ( vrbflg ) fprintf( stdout, "%e %e %e %e %e\n", 
        galacticPar.m1, galacticPar.m2,
        galacticPar.rho * cos( galacticPar.lGal ),
        galacticPar.rho * sin( galacticPar.lGal ),
        galacticPar.z );
#endif

    /* create the sim_inspiral table */
    if ( injections.simInspiralTable )
    {
      this_inj = this_inj->next = (SimInspiralTable *)
        LALCalloc( 1, sizeof(SimInspiralTable) );
    }
    else
    {
      injections.simInspiralTable = this_inj = (SimInspiralTable *)
        LALCalloc( 1, sizeof(SimInspiralTable) );
    }

    /* set the geocentric end time of the injection */
    /* XXX CHECK XXX */
    this_inj->geocent_end_time = gpsStartTime;
    if ( timeInterval )
    {
      LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
      LAL_CALL( LALAddFloatToGPS( &status, &(this_inj->geocent_end_time),
          &(this_inj->geocent_end_time), u * timeInterval ), &status );
    }
    gpsAndAcc.gps = this_inj->geocent_end_time;

    /* set gmst */
    LAL_CALL( LALGPStoGMST1( &status, &(this_inj->end_time_gmst),
          &(this_inj->geocent_end_time), &gmstUnits ), &status);
    /* XXX END CHECK XXX */

    /* populate the sim_inspiral table */

    if (mdistr == 1)
    /* uniformly distributed mass1 and uniformly distributed mass2 */
    {
      LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
      this_inj->mass1 = minMass + u * deltaM;
      LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
      this_inj->mass2 = minMass + u * deltaM;
      mtotal = this_inj->mass1 + this_inj->mass2 ;
      this_inj->eta = this_inj->mass1 * this_inj->mass2 / ( mtotal * mtotal );
    }
    else if (mdistr == 0)
    /*uniformly distributed total mass */
    {
      LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status);
      mtotal = 2.0 * minMass + u * 2.0 *deltaM ;
      LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
      this_inj->mass1 = minMass + u * deltaM;
      this_inj->mass2 = mtotal - this_inj->mass1;
      while (this_inj->mass1 >= mtotal || 
          this_inj->mass2 >= maxMass || this_inj->mass2 <= minMass )
      {
        LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
        this_inj->mass1 = minMass + u * deltaM;
        this_inj->mass2 = mtotal - this_inj->mass1;
      }
      this_inj->eta = this_inj->mass1 * this_inj->mass2 / ( mtotal * mtotal );
     }

     /* spatial distribution */

#if 0
     LAL_CALL( LALUniformDeviate( &status, &u, randParams ),
            &status );
     sinphi = 2.0 * u - 1.0;
     cosphi = sqrt( 1.0 - sinphi*sinphi );
#endif

     if (ddistr == 0)
     /* uniform distribution in distance */
     {
       REAL4 deltaD = dmax - dmin ;
       LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
       this_inj->distance = dmin + deltaD * u ;
      }
      else if (ddistr == 1)
      /* uniform distribution in log(distance) */
      {
        REAL4 lmin = log10(dmin);
        REAL4 lmax = log10(dmax);
        REAL4 deltaL = lmax - lmin;
        LAL_CALL(  LALUniformDeviate(&status,&u,randParams),&status );
        exponent = lmin + deltaL * u;
        this_inj->distance = pow(10.0,(REAL4) exponent);
      }
     else if (ddistr == 2)
     /* uniform volume distribution */
     {
       REAL4 d2min = dmin * dmin ;
       REAL4 d2max = dmax * dmax ;
       REAL4 deltad2 = d2max - d2min ;
       LAL_CALL(  LALUniformDeviate(&status,&u,randParams),&status );
       d2 = d2min + u * deltad2 ;
       this_inj->distance = sqrt(d2);
     }

     this_inj->distance = this_inj->distance / 1000.0; /*convert to Mpc */
       

      /* compute random longitude and latitude */
      LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
      this_inj->longitude = LAL_PI * u - LAL_PI_2 ;
      LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
      this_inj->latitude = LAL_TWOPI * u ;
     

#if 0
    LAL_CALL( LALGalacticInspiralParamsToSimInspiralTable( &status,
          this_inj, &galacticPar, randParams ), &status );
    if (vrbflg)
    { 
      fprintf( stdout, "%e\n",
      sqrt(galacticPar.z*galacticPar.z+galacticPar.rho*galacticPar.rho
          + Rcore*Rcore + 2.0*Rcore*galacticPar.rho*cos(galacticPar.lGal))-
      this_inj->distance*1000.0);
    }
#endif


    /* set the source and waveform fields */
    LALSnprintf( this_inj->source, LIGOMETA_SOURCE_MAX * sizeof(CHAR), "???" );
    memcpy( this_inj->waveform, waveform, LIGOMETA_WAVEFORM_MAX *
        sizeof(CHAR));

    /* XXX CHECK XXX */
    /* compute random inclination, polarization and coalescence phase */
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    this_inj->inclination = acos( 2.0 * u - 1.0 );
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    this_inj->polarization = LAL_TWOPI * u ;
    LAL_CALL( LALUniformDeviate( &status, &u, randParams ), &status );
    this_inj->coa_phase = LAL_TWOPI * u ;
    
    /* set up params for the site end times and detector response */
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

    placeAndGPS.p_gps = &(this_inj->geocent_end_time);

    detTimeAndSource.p_det_and_time = &placeAndGPS;
    detTimeAndSource.p_source = &skyPos;

    detAndSource.pSource = &source;

    gpsAndAcc.accuracy = LALLEAPSEC_STRICT;
    gpsAndAcc.gps = this_inj->geocent_end_time;

    /*
    * compute site end times
    * (copied from SnglInspiralUtils.c)
    */

    /* initialize end times with geocentric value */
    this_inj->h_end_time = this_inj->l_end_time = this_inj->geocent_end_time;

    /* lho */
    placeAndGPS.p_detector = &lho;
    LAL_CALL( LALTimeDelayFromEarthCenter( &status, &time_diff_ns,
          &detTimeAndSource ), &status );
    LAL_CALL( LALAddFloatToGPS( &status, &(this_inj->h_end_time),
          &(this_inj->h_end_time), time_diff_ns ), &status );

    /* llo */
    placeAndGPS.p_detector = &llo;
    LAL_CALL( LALTimeDelayFromEarthCenter( &status,  &time_diff_ns,
          &detTimeAndSource ), &status);
    LAL_CALL( LALAddFloatToGPS( &status,  &(this_inj->l_end_time),
          &(this_inj->l_end_time), time_diff_ns ), &status);

    /* temporarily, populate the fields for the */
    /* GEO, TAMA and VIRGO times                */

    memset( &(this_inj->g_end_time), 0,sizeof(&(this_inj->l_end_time)) ); 
    memset( &(this_inj->t_end_time), 0,sizeof(&(this_inj->l_end_time)) );
    memset( &(this_inj->v_end_time), 0,sizeof(&(this_inj->l_end_time)) );

     /*
     * compute the effective distance of the inspiral
     * (copied from SnglInspiralUtils.c)
     */
      
     /* initialize distances with real distance and compute splus and scross*/
     this_inj->eff_dist_h = this_inj->eff_dist_l = 2.0 * this_inj->distance;
     cosiota = cos( this_inj->inclination );
     splus = -( 1.0 + cosiota * cosiota );
     scross = -2.0 * cosiota;
              
     /* compute the response of the LHO detectors */
     detAndSource.pDetector = &lho;
     LAL_CALL( LALComputeDetAMResponse( &status, &resp, &detAndSource,
                  &gpsAndAcc ), &status );
                    
     /* compute the effective distance for LHO */
     this_inj->eff_dist_h /= sqrt( splus*splus*resp.plus*resp.plus +
                 scross*scross*resp.cross*resp.cross );
                      
     /* compute the response of the LLO detector */
     detAndSource.pDetector = &llo;
     LAL_CALL( LALComputeDetAMResponse( &status, &resp, &detAndSource,
                  &gpsAndAcc ), &status);
                            
     /* compute the effective distance for LLO */
     this_inj->eff_dist_l /= sqrt( splus*splus*resp.plus*resp.plus 
                 + scross*scross*resp.cross*resp.cross );

     /* temporarily, populate the fields for the */
     /* GEO, TAMA and VIRGO effective distances  */

     memset( &(this_inj->eff_dist_g), 0, sizeof(&(this_inj->eff_dist_g)) );
     memset( &(this_inj->eff_dist_t), 0, sizeof(&(this_inj->eff_dist_g)) );
     memset( &(this_inj->eff_dist_v), 0, sizeof(&(this_inj->eff_dist_g)) );

    /* increment the injection time */
    LAL_CALL( LALAddFloatToGPS( &status, &gpsStartTime, &gpsStartTime, 
          meanTimeStep ), &status );
    LAL_CALL( LALCompareGPS( &status, &compareGPS, &gpsStartTime, 
          &gpsEndTime ), &status );
    /* XXX END CHECK XXX */

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

  /* write the sim_inspiral table */
  if ( injections.simInspiralTable )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, sim_inspiral_table ), 
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, injections, 
          sim_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );
  }
  while ( injections.simInspiralTable )
  {
    this_inj = injections.simInspiralTable;
    injections.simInspiralTable = injections.simInspiralTable->next;
    LALFree( this_inj );
  }

  /* close the injection file */
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlfp ), &status );

  /* check for memory leaks and exit */
  LALCheckMemoryLeaks();
  return 0;
}
#endif
