/*----------------------------------------------------------------------- 
 * 
 * File Name: inspinj.c
 *
 * Author: Brady, P. R., Brown, D. A. and Crieghton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

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
#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXML.h>
#include <lal/Date.h>
#include <lal/Random.h>
#include <lal/TimeDelay.h>
#include <lalapps.h>
#include <processtable.h>

#define USAGE \
"lalapps_binj [options]\n"\
"\nDefaults are shown in brackets\n\n" \
"  --help                   display this message\n"\
"  --gps-start-time TIME    start injections at GPS time TIME (729273613)\n"\
"  --gps-end-time TIME      end injections at GPS time TIME (734367613)\n"\
"  --time-step STEP         space injections STEP / pi seconds appart (2630)\n"\
"  --coordinates COORDS     coordinate system to use for injections\n"\
"  --flow FLOW              first frequency of injection (150.0)\n"\
"  --fhigh FHIGH            only inject frequencies smaller than FHIGH (1000.0)\n"\
"  --deltaf DELF            spacing of injections frequencies (0.0)\n"\
"  --quality Q              quality factor for SG waveforms TAU=Q/(sqrt(2) pi F)\n"\
"  --tau TAU                duration of SG waveforms.  Q overrides TAU setting\n"\
"  --hpeak HPEAK            amplitude of SG injection in strain units\n"\
"  --seed SEED              seed random number generator with SEED (1)\n"\
"  --waveform NAME          set waveform type to NAME (SineGaussian)\n"\
"  --user-tag STRING        set the usertag to STRING\n"\
"\n"

RCSID( "$Id$" );

#define KPC ( 1e3 * LAL_PC_SI )
#define MPC ( 1e6 * LAL_PC_SI )
#define GPC ( 1e9 * LAL_PC_SI )

#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "binj"

#define ADD_PROCESS_PARAM( pptype, format, ppvalue ) \
this_proc_param = this_proc_param->next = (ProcessParamsTable *) \
  calloc( 1, sizeof(ProcessParamsTable) ); \
  LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", \
   PROGRAM_NAME ); \
   LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--%s", \
     long_options[option_index].name ); \
     LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "%s", pptype ); \
     LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, format, ppvalue );

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

const int randm = 2147483647;
struct { int i; int y; int v[32]; } randpar;


/*
 *
 * Random number generators based on Numerical Recipes.
 *
 */


int basic_random( int i )
{
  const int a = 16807;
  const int q = 127773;
  const int r = 2836;
  int k;
  k = i/q;
  i = a*(i - k*q) - r*k;
  if (i < 0)
    i += randm;
  return i;
}

void seed_random( int seed )
{
  int n;
  while ( seed == 0 )
    seed = time( NULL );

  seed = abs( seed );
  randpar.i = seed;
  for ( n = 0; n < 8; ++n )
    randpar.i = basic_random( randpar.i );
  for ( n = 0; n < (int)(sizeof(randpar.v)/sizeof(*randpar.v)); ++n )
    randpar.v[n] = randpar.i = basic_random( randpar.i );
  randpar.y = randpar.v[0];
  return;
}

int my_random( void )
{
  int ans;
  int ndiv;
  int n;

  ndiv = 1 + (randm-1)/(sizeof(randpar.v)/sizeof(*randpar.v));
  n = randpar.y/ndiv;
  ans = randpar.y = randpar.v[n];
  randpar.v[n] = randpar.i = basic_random( randpar.i );

  return ans;
}

double my_urandom( void )
{
  double u;
  int i;
  i = my_random();
  u = (double)(i) / (double)(randm + 1.0);
  return u;
}


/* 
 *
 * computes Greenwich mean sidereal time in radians (2pi rad per day) 
 *
 */


double greenwich_mean_sidereal_time( int gpssec, int gpsnan, int taiutc )
{
  /* cf. S. Aoki et al., A&A 105, 359 (1982) eqs. 13 & 19 */
  /* also cf. http://aa.usno.navy.mil */
  /* Note: 00h UT 01 Jan 2000 has JD=2451544.5 and GPS=630720013 */
  const double JD_12h_01_Jan_2000     = 2451545.0;
  const double JD_00h_01_Jan_2000     = 2451544.5;
  const double GPS_00h_01_Jan_2000    = 630720013;
  const double TAIUTC_00h_01_Jan_2000 = 32; /* leap seconds: TAI - UTC */

  double t;
  double dpU;
  double TpU;
  double gmst;

  /* compute number of seconds since 00h UT 01 Jan 2000 */
  t  = gpssec - GPS_00h_01_Jan_2000;
  t += 1e-9 * gpsnan;
  t += taiutc - TAIUTC_00h_01_Jan_2000;

  /* compute number of days since 12h UT 01 Jan 2000 */
  dpU  = floor( t / ( 24.0 * 3600.0 ) ); /* full days since 0h UT 01 Jan 2000 */
  dpU += JD_00h_01_Jan_2000 - JD_12h_01_Jan_2000; /* i.e., -0.5 */

  /* compute number of centuries since 12h UT 31 Dec 1899 */
  TpU = dpU / 36525.0;

  /* compute the gmst at 0h of the current day */
  gmst = 24110.54841
    + TpU * ( 8640184.812866
        + TpU * ( 0.093104
          - TpU * 6.2e-6 ) ); /* seconds */

  /* add the sidereal time since the start of the day */
  t = fmod( t, 24.0 * 3600.0 ); /* seconds since start of day */
  gmst += t * 1.002737909350795; /* corrections omitted */

  /* convert to fractions of a day and to radians */
  gmst = fmod( gmst / ( 24.0 * 3600.0 ), 1.0 ); /* fraction of day */
  gmst *= 2.0 * LAL_PI; /* radians */
  return gmst;
}

/* output format for LIGO-TAMA simulations */
int ligo_tama_output(FILE *fpout, SimBurstTable *simBursts)
{
  SimBurstTable *thisEvent=NULL;

  thisEvent = simBursts;
  fprintf(fpout,"# $I""d$\n");
  fprintf(fpout,"# %s\n",thisEvent->waveform);
  fprintf(fpout,"# geocent_peak_time\tnSec\tdtminus\t\tdtplus\t\tlongitude\tlatitude\tpolarization\tcoordinates\thrss\thpeak\tfreq\ttau\n");
  fflush(fpout);

  while ( thisEvent ){
    fprintf(fpout,"%0d\t%0d\t%f\t%f\t%f\t%f\t%f\t%s\t%e\t%e\t%f\t%f\n",
        thisEvent->geocent_peak_time.gpsSeconds,
        thisEvent->geocent_peak_time.gpsNanoSeconds,
        thisEvent->dtminus,
        thisEvent->dtplus,
        thisEvent->longitude,
        thisEvent->latitude,
        thisEvent->polarization,
        thisEvent->coordinates,
        thisEvent->hrss,
        thisEvent->hpeak,
        thisEvent->freq,
        thisEvent->tau
        );
    thisEvent = thisEvent->next;
  }
  fprintf(fpout,"# $I""d$\n");

  return 0;
}


int main( int argc, char *argv[] ){

  const long S2StartTime   = 729273613;  /* Feb 14 2003 16:00:00 UTC */
  const long S2StopTime    = 734367613;  /* Apr 14 2003 15:00:00 UTC */
  long gpsStartTime = S2StartTime;
  long gpsEndTime = S2StopTime;
  double meanTimeStep = 210 / LAL_PI; /* seconds between injections     */
  long long tinj = 1000000000LL * gpsStartTime;

  size_t ninj;
  int rand_seed = 1;
  RandomParams *randParams=NULL;
  REAL4 deviate=0.0;

  /* waveform */
  CHAR waveform[LIGOMETA_WAVEFORM_MAX];
  INT4   use_quality=0;
  REAL4  tau=0.1;
  REAL4  quality=0;
  REAL4  freq=150.0;
  REAL4  flow=150.0;
  REAL4  fhigh=1000.0;
  REAL4  deltaf=0.0;
  REAL4  hpeak=1.0e-20;

  /* site end time */
  CHAR                  coordinates[LIGOMETA_COORDINATES_MAX];
  LALPlaceAndGPS       *place_and_gps;
  LALDetector           lho = lalCachedDetectors[LALDetectorIndexLHODIFF];
  LALDetector           llo = lalCachedDetectors[LALDetectorIndexLLODIFF];
  SkyPosition	       *sky_pos;
  DetTimeAndASource    *det_time_and_source;
  REAL8			time_diff;
  REAL8                 site_time;

  /* xml output data */
  CHAR                  fname[256];
  CHAR                 *userTag = NULL;
  LALStatus             status = blank_status;
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         injections;
  ProcessParamsTable   *this_proc_param=NULL;
  LIGOLwXMLStream       xmlfp;

  SimBurstTable *this_sim_burst=NULL;


  /* getopt arguments */
  struct option long_options[] =
  {
    {"help",                          no_argument, 0,                'h'},
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"coordinates",             required_argument, 0,                'c'},
    {"flow",                    required_argument, 0,                'd'},
    {"fhigh",                   required_argument, 0,                'e'},
    {"deltaf",                  required_argument, 0,                'f'},
    {"quality",                 required_argument, 0,                'g'},
    {"seed",                    required_argument, 0,                's'},
    {"time-step",               required_argument, 0,                't'},
    {"waveform",                required_argument, 0,                'w'},
    {"tau",                     required_argument, 0,                'x'},
    {"freq",                    required_argument, 0,                'y'},
    {"hpeak",                   required_argument, 0,                'z'},
    {"user-tag",                required_argument, 0,                'Z'},
    {0, 0, 0, 0}
  };
  int c;

  /* set up inital debugging values */
  lal_errhandler = LAL_ERR_EXIT;
  set_debug_level( "LALMSGLVL2" );

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

  /* clear the waveform field and set the coordinates */
  memset( waveform, 0, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR) );
  memset( coordinates, 0, LIGOMETA_COORDINATES_MAX * sizeof(CHAR) );
  LALSnprintf( coordinates, LIGOMETA_COORDINATES_MAX * sizeof(CHAR), 
      "EQUATORIAL" );

  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    long int gpsinput;
    size_t optarg_len;

    c = getopt_long_only( argc, argv, 
        "ha:b:t:s:w:x:y:z:c:Z:", long_options, &option_index );

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
        gpsStartTime = gpsinput;
        tinj = 1000000000LL * gpsStartTime;
        this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "int", "%ld", gpsinput );
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
        gpsEndTime = gpsinput;
        this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "int", "%ld", gpsinput );
        break;

      case 'd':
        flow = atof( optarg );
        this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "float", "%e", flow );
        break;

      case 'e':
        fhigh = atof( optarg );
        this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "float", "%e", fhigh );
        break;

      case 'f':
        deltaf = atof( optarg );
        this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "float", "%e", deltaf );
        break;

      case 'g':
        quality = atof( optarg );
        use_quality = 1;
        this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "float", "%e", deltaf );
        break;

      case 's':
        rand_seed = atoi( optarg );
        this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "int", "%d", rand_seed );
        break;

      case 't':
        {
          double tstep = atof( optarg );
          meanTimeStep = tstep / LAL_PI;
          this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "float", "%le", tstep );
        }
        break;

      case 'w':
        LALSnprintf( waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "%s",
            optarg );
        this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "string", "%le", optarg );
        break;


      case 'x':
        {
          tau = atof( optarg );
          this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "float", "%e", tau );
        }
        break;

      case 'y':
        {
          freq = atof( optarg );
          this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "float", "%e", freq );
        }
        break;

      case 'z':
        {
          hpeak = atof( optarg );
          this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "float", "%e", hpeak );
        }
        break;

      case 'c':
        LALSnprintf( coordinates, LIGOMETA_COORDINATES_MAX * sizeof(CHAR), "%s",
            optarg );
        this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "string", "%le", optarg );
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


  /* initialize random number generator */
  LAL_CALL( LALCreateRandomParams( &status, &randParams, rand_seed ), &status);

  if ( ! *waveform )
  {
    /* default to SineGaussian */
    LALSnprintf( waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), 
        "SineGaussian" );
  }

  this_proc_param = procparams.processParamsTable;
    procparams.processParamsTable = procparams.processParamsTable->next;
    free( this_proc_param );

    /* create the detector time map structures */
    place_and_gps = (LALPlaceAndGPS *) calloc( 1, sizeof(LALPlaceAndGPS) );
    sky_pos =  (SkyPosition *) 	calloc( 1, sizeof(SkyPosition ) );
    det_time_and_source = (DetTimeAndASource *) calloc( 1, sizeof(DetTimeAndASource));
    sky_pos->system = COORDINATESYSTEM_EQUATORIAL;

    /* create the first injection */
    this_sim_burst = injections.simBurstTable = (SimBurstTable *)
      calloc( 1, sizeof(SimBurstTable) );

    /* make injection times at intervals of 100/pi seconds */
    ninj = 1;
    freq = flow;
    while ( 1 )
    {
      long tsec;
      long tnan;
      double gmst;

      /* compute tau if quality was specified */
      if (use_quality) 
        tau = quality / ( sqrt(2.0) * LAL_PI * freq);

      tinj += (long long)( 1e9 * meanTimeStep );
      if ( tinj > 1000000000LL * gpsEndTime )
      {
        break;
      }

      if ( ninj == 1 )
      {
        /* create the first injection */
        this_sim_burst = injections.simBurstTable = (SimBurstTable *)
          calloc( 1, sizeof(SimBurstTable) );
      }
      else
      {
        this_sim_burst = this_sim_burst->next = (SimBurstTable *)
          calloc( 1, sizeof(SimBurstTable) );
      }

      ++ninj;

      /* GPS time of burst */
      tsec = this_sim_burst->geocent_peak_time.gpsSeconds = 
        (long)( tinj / 1000000000LL );
      tnan = this_sim_burst->geocent_peak_time.gpsNanoSeconds = 
        (long)( tinj % 1000000000LL );

      /* get gmst (radians) */
      gmst =  greenwich_mean_sidereal_time( tsec, tnan, 32 );

      /* save gmst (hours) in sim_burst table */
      this_sim_burst->peak_time_gmst = gmst * 12.0 / LAL_PI;

      /* populate the sim burst table */
      memcpy( this_sim_burst->waveform, waveform, 
          sizeof(CHAR) * LIGOMETA_WAVEFORM_MAX );
      memcpy( this_sim_burst->coordinates, coordinates, 
          sizeof(CHAR) * LIGOMETA_COORDINATES_MAX );

      LAL_CALL( LALUniformDeviate ( &status, &deviate, randParams ), &status);
      this_sim_burst->longitude = 2.0 * LAL_PI * deviate;

      LAL_CALL( LALUniformDeviate ( &status, &deviate, randParams ), &status);
      this_sim_burst->latitude = LAL_PI/2.0 - acos(2.0*deviate-1.0) ;

      LAL_CALL( LALUniformDeviate ( &status, &deviate, randParams ), &status);
      this_sim_burst->polarization = 2.0 * LAL_PI * deviate;

      this_sim_burst->hrss = sqrt( sqrt(2.0 * LAL_PI) * tau / 4.0 ) * hpeak;
      this_sim_burst->hpeak = hpeak;
      this_sim_burst->freq = freq;
      this_sim_burst->tau = tau;
      this_sim_burst->dtplus = 4.0 * tau;
      this_sim_burst->dtminus = 4.0 * tau;
      this_sim_burst->zm_number = 0;

      place_and_gps->p_gps = &(this_sim_burst->geocent_peak_time);
      det_time_and_source->p_det_and_time = place_and_gps;
      sky_pos->longitude = this_sim_burst->longitude;
      sky_pos->latitude = this_sim_burst->latitude;
      det_time_and_source->p_source = sky_pos;

      /* compute site arrival time for lho */
      place_and_gps->p_detector = &lho;
      LAL_CALL( LALTimeDelayFromEarthCenter( &status, &time_diff, 
            det_time_and_source), &status );
      LAL_CALL( LALGPStoFloat( &status, &site_time, 
            &(this_sim_burst->geocent_peak_time) ), &status );
      site_time += time_diff;
      LAL_CALL( LALFloatToGPS( &status, &(this_sim_burst->h_peak_time),
            &site_time ), &status );

      /* compute site arrival time for llo */
      place_and_gps->p_detector = &llo;
      LAL_CALL( LALTimeDelayFromEarthCenter( &status, &time_diff, 
            det_time_and_source), &status );
      LAL_CALL( LALGPStoFloat( &status, &site_time, 
            &(this_sim_burst->geocent_peak_time) ), &status );
      site_time += time_diff;
      LAL_CALL( LALFloatToGPS( &status, &(this_sim_burst->l_peak_time),
            &site_time ), &status );

      /* increment to next frequency and test it's still in band */
      freq += deltaf;
      if ( freq > fhigh ) freq=flow;

    }

    memset( &xmlfp, 0, sizeof(LIGOLwXMLStream) );

    if ( userTag )
    {
      LALSnprintf( fname, sizeof(fname), "HL-INJECTIONS_%d_%s-%d-%d.xml", 
          rand_seed, userTag, gpsStartTime, gpsEndTime - gpsStartTime );
    }
    else
    {
      LALSnprintf( fname, sizeof(fname), "HL-INJECTIONS_%d-%d-%d.xml", 
          rand_seed, gpsStartTime, gpsEndTime - gpsStartTime );
    }

    LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlfp, fname), &status );

    LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->end_time),
          &accuracy ), &status );
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, process_table ), 
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, proctable, 
          process_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );

    if ( procparams.processParamsTable )
    {
      LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, process_params_table ), 
          &status );
      LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, procparams, 
            process_params_table ), &status );
      LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );
    }

    if ( injections.simBurstTable )
    {
      LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, sim_burst_table ), 
          &status );
      LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, injections, 
            sim_burst_table ), &status );
      LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );
    }

    LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlfp ), &status );

    {
      FILE *fpout=NULL;

      if ( userTag )
      {
        LALSnprintf( fname, sizeof(fname), "HLT-INJECTIONS_%d_%s-%d-%d.txt", 
            rand_seed, userTag, gpsStartTime, gpsEndTime - gpsStartTime );
      }
      else
      {
        LALSnprintf( fname, sizeof(fname), "HLT-INJECTIONS_%d-%d-%d.txt", 
            rand_seed, gpsStartTime, gpsEndTime - gpsStartTime );
      }
      fpout = fopen(fname,"w");
      ligo_tama_output(fpout,injections.simBurstTable);
      fclose(fpout);
    }

    return 0;
}

