/*----------------------------------------------------------------------- 
 * 
 * File Name: binj.c
 *
 * Author: Brady, P. R., Brown, D. A., Crieghton, J. D. E. and Ray Majumder S
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

int snprintf(char *str, size_t size, const  char  *format, ...);


#define USAGE \
"lalapps_binj [options]\n"\
"\nDefaults are shown in brackets\n\n" \
"  --help                   display this message\n"\
"  --gps-start-time TIME    start injections at GPS time TIME (729273613)\n"\
"  --gps-end-time TIME      end injections at GPS time TIME (734367613)\n"\
"  --time-step STEP         space injections STEP / pi seconds appart (210)\n"\
"  --coordinates COORDS     coordinate system to use for injections\n"\
"  --flow FLOW              first frequency of injection (150.0)\n"\
"  --fhigh FHIGH            only inject frequencies smaller than FHIGH (1000.0)\n"\
"  --deltaf DELF            spacing of injections frequencies (0.0)\n"\
"  --quality Q              quality factor for SG waveforms TAU=Q/(sqrt(2) pi F)\n"\
"  --tau TAU                duration of SG waveforms.  Q overrides TAU setting\n"\
"  --hpeak HPEAK            amplitude of SG injection in strain units\n"\
"  --log-hpeak-min LOGHMIN  min amplitude of SG injection in strain units\n"\
"  --log-hpeak-max LOGHMAX  max amplitude of SG injection in strain units\n"\
"  --min-distance           min distance of source in Kpc(default 100Kpc) \n"\
"  --max-distance           max distance of source in Kpc(default 10000Kpc) \n"\
"  --d-distr                distance distribution ( 0 = logarithmic(only one as of now!! ) \n"\
"  --simwaveform-duration   duration of th esimulated waveform (Warren/Ott/ZM)\n"\
"  --simwaveform-min-number min # of the simulated waveform \n"\
"  --simwaveform-max-number max # of the simulated waveform \n"\
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

#define TRUE      1
#define FALSE     0


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
  REAL4 iamp;
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
  INT4   useRandomStrain=0;
  REAL4  log_hpeakMin=0.0;
  REAL4  log_hpeakMax=0.0;
  REAL4  logAmpRange=0.0;
  REAL4  simwavedur = 0.0;
  INT4   minsimnumber = 0;
  INT4   maxsimnumber = 10;
  INT4   deltasim = 0;
  REAL4  dmin = 100;
  REAL4  dmax = 10000;
  REAL4  deltaL = 0;
  REAL4  logdmin = 0;
  int    ddistr = 0;

  /* Need to set this to true if one wants to use MDC signals*/
  INT4 mdcFlag =FALSE;
 
  /* site end time */
  CHAR                  coordinates[LIGOMETA_COORDINATES_MAX];
  INT4                  useZenith=0;
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
  MetadataTable         mdcinjections;
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
    {"log-hpeak-min",           required_argument, 0,                'j'},
    {"log-hpeak-max",           required_argument, 0,                'k'},
    {"max-distance",            required_argument, 0,                'm'},
    {"min-distance",            required_argument, 0,                'n'},
    {"d-distr",                 required_argument, 0,                'p'},
    {"seed",                    required_argument, 0,                's'},
    {"time-step",               required_argument, 0,                't'},
    {"waveform",                required_argument, 0,                'w'},
    {"simwaveform-duration",    required_argument, 0,                'i'},
    {"simwaveform-min-number",  required_argument, 0,                'l'},
    {"simwaveform-max-number",  required_argument, 0,                'q'},
    {"tau",                     required_argument, 0,                'x'},
    {"freq",                    required_argument, 0,                'y'},
    {"hpeak",                   required_argument, 0,                'z'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"mdcFlag",                 no_argument,        &mdcFlag,      TRUE },
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
        "ha:b:t:s:w:i:d:e:f:g:j:k:m:n:p:q:l:x:y:z:c:Z:", long_options, &option_index );

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
        this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "float", "%e", quality );
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
        this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "string", "%s", optarg );
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

      case 'j':
        {
          useRandomStrain+=1;
          log_hpeakMin = atof( optarg );
          this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "float", "%e", log_hpeakMin );
        }
        break;

      case 'k':
        {
          useRandomStrain+=1;
          log_hpeakMax = atof( optarg );
          this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "float", "%e", log_hpeakMax );
        }
        break;

      case 'm':
        {
          dmax = atof( optarg );
          this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "float", "%e", dmax );
        }
        break;

      case 'n':
        {
          dmin = atof( optarg );
          this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "float", "%e", dmin );
        }
        break;

      case 'p':
        {
          ddistr = atoi( optarg );
	  if ( ddistr != 0 )
	    {
	      fprintf( stderr, "invalid argument to --d-distr:\n"
		       "ddistr must be 0 as of now\n"
		       );
	      exit(1);
	    }
          this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "int", "%d", ddistr );
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

        if( !(strcmp( coordinates, "ZENITH" )) )
          useZenith=1;

        this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "string", "%s", optarg );
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

      case 'i':
        {
          simwavedur = atof( optarg );
          this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "float", "%e", simwavedur );
        }
        break;

      case 'l':
        {
          minsimnumber = atoi( optarg );
          this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "int", "%d", minsimnumber );
        }
        break;

      case 'q':
        {
          maxsimnumber = atoi( optarg );
          this_proc_param = this_proc_param->next = next_process_param( long_options[option_index].name, "int", "%d", maxsimnumber );
        }
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

  /* check some of the input parameters for consistency */
  if( (useRandomStrain==1) )
  {
        fprintf( stderr, "Must supply upper and lower limits when using"
            "random strain\n" );
        exit( 1 );
  }
  else if ( (useRandomStrain == 2) )
  {
    logAmpRange = (log_hpeakMax - log_hpeakMin);
  }

  if (ddistr == 0)
  {
    logdmin = log10(dmin);
    deltaL = log10(dmax) - logdmin;
  }

  deltasim = maxsimnumber - minsimnumber;

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
    iamp = 1.0;
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


      /* sky location and polarizatoin angle */
      if(useZenith)
      {
        this_sim_burst->longitude = 0.0;
        this_sim_burst->latitude = 0.0;
        this_sim_burst->polarization = 0.0;
      }
      else
      {
        LAL_CALL( LALUniformDeviate ( &status, &deviate, randParams ), &status);
        this_sim_burst->longitude = 2.0 * LAL_PI * deviate;

        LAL_CALL( LALUniformDeviate ( &status, &deviate, randParams ), &status);
        this_sim_burst->latitude = LAL_PI/2.0 - acos(2.0*deviate-1.0) ;

        LAL_CALL( LALUniformDeviate ( &status, &deviate, randParams ), &status);
        this_sim_burst->polarization = 2.0 * LAL_PI * deviate;
      }

      /* compute amplitude information */
      if(useRandomStrain)
      {
	LAL_CALL( LALUniformDeviate ( &status, &deviate, randParams ), &status);
	hpeak = pow(10, logAmpRange * deviate + log_hpeakMin);
	/*Uncomment the lines below and comment the above two line
	 *if want to inject uniformly spaced amplitude waveforms
	 */
	/*if (iamp == 1.0 ){
	hpeak = pow(10,-19);
	iamp++;
	}
	else if( iamp == 2.0 ){
	  hpeak = (REAL4)2.0*pow(10,-19);
	  iamp++;
	}
	else if (iamp == 3.0 ){
	  hpeak = (REAL4)3.0*pow(10,-19);
	  iamp++;
	}
	else if (iamp == 4.0 ){
	  hpeak = (REAL4)3.3*pow(10,-18);
	  iamp++;
	}
	else if (iamp == 5.0 ){
	  hpeak = (REAL4)5.3*pow(10,-17);
	  iamp = 1.0;
	  }*/
      }

      if (ddistr == 0)
      /* uniform distribution in log(distance) */
      {
	LAL_CALL(  LALUniformDeviate( &status, &deviate, randParams), &status );
	this_sim_burst->distance = pow(10.0, (REAL4)( logdmin + deltaL * deviate ) );
      }


      /* deal with the intrinsic signal parameters */
      this_sim_burst->hrss = sqrt( sqrt(2.0 * LAL_PI) * tau / 4.0 ) * hpeak;
      this_sim_burst->hpeak = hpeak;
      this_sim_burst->freq = freq;
      this_sim_burst->tau = tau;
      if (simwavedur){
	this_sim_burst->dtplus = simwavedur / 2.0;
	this_sim_burst->dtminus = simwavedur / 2.0 ;
      }
      else{
	this_sim_burst->dtplus = 4.0 * tau;
	this_sim_burst->dtminus = 4.0 * tau;
      }

      /* set the simulated wavenumber */
      LAL_CALL(  LALUniformDeviate( &status, &deviate, randParams), &status );
      this_sim_burst->zm_number = (INT4)(minsimnumber + deltasim * deviate);

      /* set up for site arrival time calculation */
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

    /* if using the MDC signal frames for injections:
     * then read in the ascii file of inj. parameters
     * and write a sim burst table out of that. However
     * currently all the fields are not filled.[20040430: SKRM]
     */
    if (mdcFlag )
      {
	INT4   n = 0;
	INT4   nn, x, rc;
	CHAR mdcSimfile[20];
	FILE *fp;
	float *fHp,*fHx,*fLp,*fLx,*theta,*phi,*psi,*hrss;
	float *gg,*hh,*ii,*jj,*kk,*ll,*mm,*rr;
	int *gps,*tH,*tL,*tauHL;
	int *cc,*dd,*ee,*ff;
	char **waveform;
	char *qq, *y;
	
	SimBurstTable *this_mdcsim_burst=NULL;
     
	/*read in the ascii file containing the injection parameters 
	 *the ascii file is assumed to be present in the working dir. 
	 */
	snprintf(mdcSimfile,20,"mdcsim_%d.dat",1);
	fp = fopen(mdcSimfile,"r");
	if ( fp == NULL ){
	  fprintf(stdout,"Error:Must have file mdcsim_1.dat in the working dir.\n");
	  exit ( 1 );
	}
	while (fscanf(fp,"%*d %*d %*d %*d %*e %*e %*e %*e %*e %*e %*e %*s %*e")!= EOF)
	  ++n;
	rewind(fp);

	gps = cc = (int*)malloc(n*sizeof(int));
	tL = dd = (int*)malloc(n*sizeof(int));
	tH = ee = (int*)malloc(n*sizeof(int));
	tauHL = ff = (int*)malloc(n*sizeof(int));
	fHp = gg = (float*)malloc(n*sizeof(float));
	fHx = hh = (float*)malloc(n*sizeof(float));
	fLp = ii = (float*)malloc(n*sizeof(float));
	fLx = jj = (float*)malloc(n*sizeof(float));
	theta = kk = (float*)malloc(n*sizeof(float));
	phi = ll = (float*)malloc(n*sizeof(float));
	psi = mm = (float*)malloc(n*sizeof(float));
	hrss = rr = (float*)malloc(n*sizeof(float));
	qq = (char*)malloc((n*10)*sizeof(char)); /*memory allocation may need to be 
						  *changed if the name of waveform
						  *changes in the mdc file 
						  */
	waveform = (char**)malloc(n*10*sizeof(char));

	/* reads in the diff. parameters from the text file */
	nn = 0;
	while ((rc = fscanf(fp,"%d %d %d %d %e %e %e %e %e %e %e %s %e", cc++, 
			    dd++, ee++, ff++, gg++, hh++, ii++, jj++, kk++, 
			    ll++, mm++, qq, rr++)) == 13) 
	  {
	    waveform[nn++] = qq; /* scans the name of the waveform */
	    qq = qq + 10;
	    continue;
	  }

	if(rc != EOF)
	  printf("ERROR READING FILE\n");

	fclose(fp);

	/* create the first injection */
	this_mdcsim_burst = mdcinjections.simBurstTable = (SimBurstTable *)
	  calloc( 1, sizeof(SimBurstTable) );
	/*create the corresponding simburst table */
	for(x=0;x<n;x++)
	  {
	    this_mdcsim_burst->geocent_peak_time.gpsSeconds=gps[x];/* wrong value*/
	    this_mdcsim_burst->geocent_peak_time.gpsNanoSeconds=0;/* wrong value*/
	    this_mdcsim_burst->h_peak_time.gpsSeconds=gps[x];
	    this_mdcsim_burst->h_peak_time.gpsNanoSeconds=tH[x] + (0.5 * 1000000000LL) ; 
	    this_mdcsim_burst->l_peak_time.gpsSeconds=gps[x];
	    this_mdcsim_burst->l_peak_time.gpsNanoSeconds=tL[x] + (0.5 * 1000000000LL) ;
	    this_mdcsim_burst->longitude=phi[x];
	    this_mdcsim_burst->latitude=theta[x];
	    this_mdcsim_burst->polarization=psi[x];
	    this_mdcsim_burst->hrss=hrss[x];

	    /*rip off the frequency from the name of the waveform */
	    waveform[x] = waveform[x] + 2;
	    y = (char*)malloc(10*sizeof(char));
	    strncpy(y,waveform[x],3);
	    this_mdcsim_burst->freq = atof(y); /*freq of the SineGaussian*/
	    free(y);
	    this_mdcsim_burst = this_mdcsim_burst->next = (SimBurstTable *)
          calloc( 1, sizeof(SimBurstTable) );
	  }


	free(gps);
	free(tH);
	free(tL);
	free(tauHL);
	free(fHp);
	free(fHx);
	free(fLp);
	free(fLx);
	free(theta);
	free(phi);
	free(psi);
	free(hrss);
      }



    memset( &xmlfp, 0, sizeof(LIGOLwXMLStream) );

    if(!mdcFlag){/* these tables are written when mdc signals are not used */
    if ( userTag )
    {
      LALSnprintf( fname, sizeof(fname), "HL-INJECTIONS_%s-%d-%d.xml", 
          userTag, gpsStartTime, gpsEndTime - gpsStartTime );
    }
    else
    {
      LALSnprintf( fname, sizeof(fname), "HL-INJECTIONS-%d-%d.xml", 
          gpsStartTime, gpsEndTime - gpsStartTime );
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
        LALSnprintf( fname, sizeof(fname), "HLT-INJECTIONS_%s-%d-%d.txt", 
            userTag, gpsStartTime, gpsEndTime - gpsStartTime );
      }
      else
      {
        LALSnprintf( fname, sizeof(fname), "HLT-INJECTIONS-%d-%d.txt", 
            gpsStartTime, gpsEndTime - gpsStartTime );
      }
      fpout = fopen(fname,"w");
      ligo_tama_output(fpout,injections.simBurstTable);
      fclose(fpout);
    }
    }
    /* this is used when mdc frames are used */
    if( mdcFlag )
      {
	LALSnprintf( fname, sizeof(fname), "HL-MDCSG10_%d.xml",1 );
	LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlfp, fname), &status );
	LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, sim_burst_table ), 
		  &status );
	LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, mdcinjections, 
					  sim_burst_table ), &status );
	LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );
	LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlfp ), &status );
      }

    return 0;
}

