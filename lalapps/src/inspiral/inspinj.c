/*
 * Program to produce ILWD & LIGOLw files containing injection times and
 * parameters for binary inspirals in a population specified by external files.
 *
 * This program produces output files injepochs.ilwd and injparams.ilwd
 * containing injection times and parameters (respectively) in ILWD format.
 * Parameters are also written as a plain textfile to logfile injlog.txt.
 * XML output is written to injections-XX.xml where XX is the input random
 * number seed.
 *
 * $Id$
 *
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
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXML.h>
#include <lal/Date.h>

#define USAGE \
"lalapps_inspinj --gps-start-time gpsStartTime --gps-end-time gpsEndTime \n"\
"    --time-step timeStep --seed seed [--datafile inspsrcs.dat] [--waveform string]\n"

int snprintf(char *str, size_t size, const  char  *format, ...);

RCSID( "$Id$" );

#ifndef M_PI
#define M_PI 3.14159265358979323846  /* pi */
#endif

#define PC  3.0856775807e16
#define KPC ( 1e3 * PC )
#define MPC ( 1e6 * PC )
#define GPC ( 1e9 * PC )

#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "inspinj"

enum { mTotElem, etaElem, distElem, incElem, phiElem, lonElem, latElem,
  psiElem, numElem };

FILE *fplog;
SimInspiralTable *this_sim_insp;


/*
 *
 * Random number generators based on Numerical Recipes.
 *
 */
const int randm = 2147483647;
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

struct { int i; int y; int v[32]; } randpar;

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

/* computes Greenwich mean sidereal time in radians (2pi rad per day) */
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
  gmst *= 2.0 * M_PI; /* radians */
  return gmst;
}

/* this is crude... may not get leap seconds exactly right */
struct tm * lho_time( int gpssec )
{
  time_t now = (time_t)gpssec + (time_t)315964811;
  putenv( "TZ=PST08PDT,M4.1.0,M10.5.0" );
  return localtime( &now );
}

/* this is crude... may not get leap seconds exactly right */
struct tm * llo_time( int gpssec )
{
  time_t now = (time_t)gpssec + (time_t)315964811;
  putenv( "TZ=CST06CDT,M4.1.0,M10.5.0" );
  return localtime( &now );
}

/* leap seconds may not be correct */
int site_time_of_day( double *lhotod, double *llotod, int gpssec, int gpsnan )
{
  struct tm *t;
  t = lho_time( gpssec );
  *lhotod  = t->tm_hour;
  *lhotod += t->tm_min / 60.0;
  *lhotod += t->tm_sec / 3600.0;
  *lhotod += gpsnan / 3600000000000.0;
  t = llo_time( gpssec );
  *llotod  = t->tm_hour;
  *llotod += t->tm_min / 60.0;
  *llotod += t->tm_sec / 3600.0;
  *llotod += gpsnan / 3600000000000.0;
  return 0;
}

double eff_dist( double nx[], double ny[], double *injPar, double gmst )
{
  double mu;
  double theta;
  double phi;
  double psi;
  double splus;
  double scross;
  double x[3];
  double y[3];
  double d[3][3];
  double eplus[3][3];
  double ecross[3][3];
  double fplus;
  double fcross;
  double deff;
  int i;
  int j;

  theta = 0.5 * M_PI - injPar[latElem];
  phi = injPar[lonElem] - gmst;
  psi = injPar[psiElem];
  mu = cos( injPar[incElem] );
  splus = -( 1.0 + mu * mu );
  scross = -2.0 * mu;

  x[0] = +( sin( phi ) * cos( psi ) - sin( psi ) * cos( phi ) * cos( theta ) );
  x[1] = -( cos( phi ) * cos( psi ) + sin( psi ) * sin( phi ) * cos( theta ) );
  x[2] = sin( psi ) * sin( theta );
  y[0] = -( sin( phi ) * sin( psi ) + cos( psi ) * cos( phi ) * cos( theta ) );
  y[1] = +( cos( phi ) * sin( psi ) - cos( psi ) * sin( phi ) * cos( theta ) );
  y[2] = cos( psi ) * sin( theta );

  fplus = 0;
  fcross = 0;
  for ( i = 0; i < 3; ++i )
  {
    for ( j = 0; j < 3; ++j )
    {
      d[i][j]  = 0.5 * ( nx[i] * nx[j] - ny[i] * ny[j] );
      eplus[i][j] = x[i] * x[j] - y[i] * y[j];
      ecross[i][j] = x[i] * y[j] + y[i] * x[j];
      fplus += d[i][j] * eplus[i][j];
      fcross += d[i][j] * ecross[i][j];
    }
  }

  deff  = 2.0 * injPar[distElem];
  deff /= sqrt( splus*splus*fplus*fplus + scross*scross*fcross*fcross );

  return deff;
}


/* convert galactic coords to equatorial coords (all angles in radians) */
int galactic_to_equatorial( double *alpha, double *delta, double l, double b )
{
  const double alphaNGP = 192.8594813 * M_PI / 180.0;
  const double deltaNGP = 27.1282511 * M_PI / 180.0;
  const double lascend  = 33.0 * M_PI / 180.0;
  double lm = l - lascend;
  *alpha  = atan2( cos(b)*cos(lm),
      sin(b)*cos(deltaNGP) - cos(b)*sin(deltaNGP)*sin(lm) );
  *alpha += alphaNGP;
  *delta  = asin( cos(b)*cos(deltaNGP)*sin(lm) + sin(b)*sin(deltaNGP) );
  return 0;
}

/* generate a sky position for a random Galactic inspiral */
int galactic_sky_position( double *dist, double *alpha, double *delta )
{
  const double h_scale = 1.5 * KPC;
  const double r_scale = 4.0 * KPC;
  const double r_sun   = 8.5 * KPC;
  double r, z, phi, rho2, l, b;
  double u;

  r = r_scale * sqrt( -2 * log( my_urandom() ) );
  u = my_urandom();
  if ( u > 0.5)
    z = -h_scale * log( 2 * ( 1 - u ) );
  else
    z = h_scale * log( 2 * u );
  phi = 2 * M_PI * my_urandom();

  rho2 = r_sun * r_sun + r * r - 2 * r_sun * r * cos( phi );
  *dist = sqrt( z * z + rho2 );
  l = atan2( r * sin( phi ), r_sun - r * cos( phi ) );
  b = asin( z / (*dist) );
  galactic_to_equatorial( alpha, delta, l, b );

  return 0;
}


const char *sourceFile = NULL;
int num_source;
struct {
  char   name[16];
  double ra;
  double dec;
  double dist;
  double lum;
  double fudge;
} *source_data;

int read_source_data( const char *datafile )
{
  const char *path;
  const char *basename = datafile ? datafile : "inspsrcs.dat";
  char fname[256];
  double *p1, *p2;
  char line[256];
  FILE *fp;
  int i;
  if ( *basename == '/' || *basename == '.' ) /* full path */
    sprintf( fname, "%s", basename );
  else
  {
    path = getenv( "LALAPPS_DATA_PATH" );
    sprintf( fname, "%s/%s", path ? path : PREFIX "/" PACKAGE "/share",
        basename );
  }
  fp = fopen( fname, "r" );
  if ( ! fp )
  {
    fprintf( stderr, "Could not find file %s\n", fname );
    fprintf( stderr, "Set environment LALAPPS_DATA_PATH to location of file %s\n", basename );
    exit( 1 );
  }
  num_source = 0;
  while ( fgets( line, sizeof( line ), fp ) )
    if ( line[0] == '#' )
      continue;
    else
      ++num_source;
  rewind( fp );
  source_data = calloc( num_source, sizeof( *source_data ) );
  if ( ! source_data )
  {
    fprintf( stderr, "alloc error\n" );
    exit( 1 );
  }
  i = 0;
  while ( fgets( line, sizeof( line ), fp ) )
    if ( line[0] == '#' )
      continue;
    else
    {
      char ra_sgn, dec_sgn;
      double ra_h, ra_m, dec_d, dec_m;
      int c;
      c = sscanf( line, "%s %c%le:%le %c%le:%le %le %le %le",
          source_data[i].name, &ra_sgn, &ra_h, &ra_m, &dec_sgn, &dec_d, &dec_m,
          &source_data[i].dist, &source_data[i].lum, &source_data[i].fudge );
      if ( c != 10 )
      {
        fprintf( stderr, "error parsing source datafile %s\n", fname );
        exit( 1 );
      }
      /* by convention, overall sign is carried only on hours/degrees entry */
      source_data[i].ra  = ( ra_h + ra_m / 60.0 ) * M_PI / 12.0;
      source_data[i].dec = ( dec_d + dec_m / 60.0 ) * M_PI / 180.0;
      if ( ra_sgn == '-' )
        source_data[i].ra *= -1;
      if ( dec_sgn == '-' )
        source_data[i].dec *= -1;
      source_data[i].dist *= KPC;
      ++i;
    }
  fclose( fp );
  return num_source;
}

/* generate a sky position for a random inspiral from the Galaxy, LMC, or SMC */
int sky_position( double *dist, double *alpha, double *delta )
{
  static double *ratio;
  static double *frac;
  static double  norm;
  static int     init;

  double u;
  int i;

  if ( ! init )
  {
    init = 1;
    norm = 1; /* milky way */
    num_source = read_source_data( sourceFile );
    ratio = calloc( num_source, sizeof( *ratio ) );
    if ( ! ratio )
    {
      fprintf( stderr, "alloc error\n" );
      exit( 1 );
    }
    frac  = calloc( num_source, sizeof( *frac  ) );
    if ( ! frac )
    {
      fprintf( stderr, "alloc error\n" );
      exit( 1 );
    }
    for ( i = 0; i < num_source; ++i )
      norm += ratio[i] = source_data[i].lum * source_data[i].fudge;
    frac[0] = ratio[0] / norm;
    for ( i = 1; i < num_source; ++i )
      frac[i] = frac[i-1] + ratio[i] / norm;
  }

  u = my_urandom();

  for ( i = 0; i < num_source; ++i )
    if ( u < frac[i] )
    {
      fprintf( fplog, "\t%s", source_data[i].name );
      *dist  = source_data[i].dist;
      *alpha = source_data[i].ra;
      *delta = source_data[i].dec;
      return 0;
    }

  /* galactic event */
  fprintf( fplog, "\tMW" );
  return galactic_sky_position( dist, alpha, delta );
}

/* generate all parameters (sky position and angles) for a random inspiral */
int inj_params( double *injPar )
{
  static double *m1arr;
  static double *m2arr;
  static size_t n;
  size_t i;
  double m1;
  double m2;
  double alpha;
  double delta;
  double dist;

  if ( ! n )
  {
    const char *path;
    const char *basename = "BNSMasses.dat";
    char fname[256];
    double *p1, *p2;
    char line[64];
    FILE *fp;
    path = getenv( "LALAPPS_DATA_PATH" );
    sprintf( fname, "%s/%s", path ? path : PREFIX "/" PACKAGE "/share",
        basename );
    fp = fopen( fname, "r" );

    if ( ! fp )
    {
      fprintf( stderr, "Could not find file %s\n", fname );
      fprintf( stderr, "Set environment LALAPPS_DATA_PATH to location of file %s\n", basename );
      exit( 1 );
    }
    while ( fgets( line, sizeof( line ), fp ) )
      ++n;
    p1 = m1arr = calloc( n, sizeof( *m1arr ) );
    if ( ! p1 )
    {
      fprintf( stderr, "alloc error\n" );
      exit( 1 );
    }
    p2 = m2arr = calloc( n, sizeof( *m2arr ) );
    if ( ! p2 )
    {
      fprintf( stderr, "alloc error\n" );
      exit( 1 );
    }
    rewind( fp );
    while ( fgets( line, sizeof( line ), fp ) )
      sscanf( line, "%le %le", p1++, p2++ );
    fclose( fp );
  }

  sky_position( &dist, &alpha, &delta );
  i = (size_t)( n * my_urandom() );

  m1 = m1arr[i];
  m2 = m2arr[i];
  injPar[mTotElem] = m1 + m2;
  injPar[etaElem]  = m1 * m2 / ( ( m1 + m2 ) * ( m1 + m2 ) );
  injPar[incElem]  = acos( -1.0 + 2.0 * my_urandom() );
  injPar[phiElem]  = 2 * M_PI * my_urandom();
  injPar[psiElem]  = 2 * M_PI * my_urandom();
  injPar[distElem] = dist;
  injPar[lonElem]  = alpha;
  injPar[latElem]  = delta;

  return 0;
}


#define UNITS "msun,none,m,rad,rad,rad,rad,rad"
struct time_list { long long tinj; struct time_list *next; };

int main( int argc, char *argv[] )
{
  double nxH[3] = { -0.2239, +0.7998, +0.5569 };
  double nyH[3] = { -0.9140, +0.0261, -0.4049 };
  double nxL[3] = { -0.9546, -0.1416, -0.2622 };
  double nyL[3] = { +0.2977, -0.4879, -0.8205 };
  const long S1StartTime   = 714150013;  /* Aug 23, 2002  08:00:00 PDT */
  const long S1StopTime    = 715618813;  /* Sep 09, 2002  08:00:00 PDT */
  long gpsStartTime = S1StartTime;
  long gpsStopTime = S1StopTime;
  double meanTimeStep = 2630 / M_PI; /* seconds between injections     */

  long long tinj              = 1000000000LL * gpsStartTime;
  struct time_list  tlisthead;
  struct time_list *tlistelem = &tlisthead;

  double injPar[numElem];
  size_t ninj;
  size_t inj;
  FILE *fp;
  int rand_seed = 1;

  /* waveform */
  CHAR waveform[LIGOMETA_WAVEFORM_MAX];

  /* xml output data */
  CHAR                  fname[256];
  LALStatus             status = blank_status;
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         injections;
  ProcessParamsTable   *this_proc_param;
  LIGOLwXMLStream       xmlfp;

  /* getopt arguments */
  struct option long_options[] =
  {
    {"help",                          no_argument, 0,                'h'},
    {"source-file",             required_argument, 0,                'f'},
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"time-step",               required_argument, 0,                't'},
    {"seed",                    required_argument, 0,                's'},
    {"waveform",                required_argument, 0,                'w'},
    {0, 0, 0, 0}
  };
  int c;

  lalDebugLevel = LALMSGLVL3;

  /* create the process and process params tables */
  proctable.processTable = (ProcessTable *) 
    LALCalloc( 1, sizeof(ProcessTable) );
  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->start_time),
        &accuracy ), &status );
  LAL_CALL( populate_process_table( &status, proctable.processTable, 
        PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE ), &status );
  this_proc_param = procparams.processParamsTable = (ProcessParamsTable *) 
    LALCalloc( 1, sizeof(ProcessParamsTable) );

  /* create the first injection */
  this_sim_insp = injections.simInspiralTable = (SimInspiralTable *)
    LALCalloc( 1, sizeof(SimInspiralTable) );

  /* clear the waveform field */
  memset( waveform, 0, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR) );

  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;

    c = getopt_long_only( argc, argv, 
        "ha:b:f:s:t:w:", long_options, &option_index );

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

      case 'f':
        sourceFile = optarg;
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
          gpsStartTime = gstartt;
          tinj              = 1000000000LL * gpsStartTime;
          /* ADD_PROCESS_PARAM( "int", "%ld", gstartt ); */
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
          gpsStopTime = gendt;
          /* ADD_PROCESS_PARAM( "int", "%ld", gendt ); */
        }
        break;

      case 's':
        {
            int s = atoi( optarg );
            rand_seed = s;
            snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
                    "%s", PROGRAM_NAME );
            snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
                    "--seed" );
            snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "int" );
            snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, "%d", rand_seed );
            /* ADD_PROCESS_PARAM( "int", "%ld", gendt ); */
        }
        break;

      case 't':
        {
            double tstep = atof( optarg );
            meanTimeStep = tstep / M_PI;
            fprintf(stdout,"Time step is %e/PI = %e\n",tstep, meanTimeStep);
            snprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, 
                    "%s", PROGRAM_NAME );
            snprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
                    "--time-step" );
            snprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "float" );
            snprintf( this_proc_param->value, LIGOMETA_TYPE_MAX, "%e", tstep );
            /* ADD_PROCESS_PARAM( "int", "%ld", gendt ); */
        }
        break;

      case 'w':
        LALSnprintf( waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "%s",
            optarg );

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

  seed_random( rand_seed );

  tlisthead.tinj = tinj;
  tlisthead.next = NULL;

  if ( ! *waveform )
  {
    /* default to Tev's GeneratePPNInspiral as used in */
    LALSnprintf( waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), 
        "GeneratePPNtwoPN" );
  }

  /* open logfile for injection parameters */
  fplog = fopen( "injlog.txt", "w" );

  /* make injection times at intervals of 100/pi seconds */
  ninj = 1;
  tlistelem = &tlisthead;
  while ( 1 )
  {
    tinj += (long long)( 1e9 * meanTimeStep );
    if ( tinj > 1000000000LL * gpsStopTime )
      break;
    tlistelem = tlistelem->next = calloc( 1, sizeof( *tlistelem ) );
    tlistelem->tinj = tinj;
    ++ninj;
  }

  /*
   *
   * First Sequence: injection epochs.
   *
   */

  fp = fopen( "injepochs.ilwd", "w" );
  fputs( "<?ilwd?>\n", fp );
  fputs( "<ilwd name='injepochs::sequence' size='7'>\n", fp );

  fprintf( fp, "\t<lstring name='real:domain' size='4'>TIME</lstring>\n" );
  fprintf( fp, "\t<int_4u name='gps_sec:start_time' units='sec'>%ld</int_4u>\n",
      gpsStartTime );
  fprintf( fp, "\t<int_4u name='gps_nan:start_time' units='nanosec'>0</int_4u>\n" );
  fprintf( fp, "\t<int_4u name='gps_sec:stop_time' units='sec'>%ld</int_4u>\n",
      gpsStopTime );
  fprintf( fp, "\t<int_4u name='gps_nan:stop_time' units='nanosec'>0</int_4u>\n");
  fprintf( fp, "\t\t<real_8 name='time:step_size' units='sec'>%e</real_8>\n",
      meanTimeStep );

  fprintf( fp, "\t<int_4u ndim='2' dims='2,%d' name='data' units='s,ns'>",
      ninj );
  fprintf( fp, "%ld 0", gpsStartTime );
  tlistelem = tlisthead.next;

  for ( inj = 1; inj < ninj; ++inj )
  {
    long tsec = (long)( tlistelem->tinj / 1000000000LL );
    long tnan = (long)( tlistelem->tinj % 1000000000LL );
    fprintf( fp, " %ld %ld", tsec, tnan );
    tlistelem = tlistelem->next;
  }
  fprintf( fp, "</int_4u>\n" );

  fputs( "</ilwd>\n", fp );
  fclose( fp );

  /*
   * 
   * Second sequence: injection parameters.
   * 
   */

  fp = fopen( "injparams.ilwd", "w" );
  fputs( "<?ilwd?>\n", fp );
  fputs( "<ilwd name='injparams::sequence' size='7'>\n", fp );

  fprintf( fp, "\t<lstring name='real:domain' size='4'>TIME</lstring>\n" );
  fprintf( fp, "\t<int_4u name='gps_sec:start_time' units='sec'>%ld</int_4u>\n",
      gpsStartTime );
  fprintf( fp, "\t<int_4u name='gps_nan:start_time' units='nanosec'>0</int_4u>\n" );
  fprintf( fp, "\t<int_4u name='gps_sec:stop_time' units='sec'>%ld</int_4u>\n",
      gpsStopTime );
  fprintf( fp, "\t<int_4u name='gps_nan:stop_time' units='nanosec'>0</int_4u>\n" );
  fprintf( fp, "\t<real_8 name='time:step_size' units='sec'>%e</real_8>\n",
      meanTimeStep );

  fprintf( fp, "\t<real_4 ndim='2' dims='%d,%d' name='data' units='" UNITS "'>",
      numElem, ninj );
    fprintf( fplog, "# GPS Time (s.ns)\tGMST (h)\tSource\tMtot (MSun)\tEta      \tDist (Mpc)\t"
              "Incl (rad)\tPhase (rad)\tRA (rad)\tDEC (rad)\tPsi (rad)\tLHO time (h)\tLLO time (h)\tDeffH (Mpc)\tDeffL (Mpc)\n" );
  tlistelem = &tlisthead;
  for ( inj = 0; inj < ninj; ++inj )
  {
    int elem;
    long tsec = this_sim_insp->geocent_end_time.gpsSeconds = 
      (long)( tlistelem->tinj / 1000000000LL );
    long tnan = this_sim_insp->geocent_end_time.gpsNanoSeconds = 
      (long)( tlistelem->tinj % 1000000000LL );
    double gmst;
    double deffH;
    double deffL;
    double todH;
    double todL;
    gmst =  greenwich_mean_sidereal_time( tsec, tnan, 32 );
    fprintf( fplog, "%ld.%09ld\t%e", tsec, tnan, 
        (this_sim_insp->end_time_gmst = gmst * 12.0 / M_PI) );
    tlistelem = tlistelem->next;
    inj_params( injPar );
    fprintf( fp, "%s%e", inj ? " " : "", injPar[0] );
    fprintf( fplog, "\t%e", injPar[0] );
    for ( elem = 1; elem < numElem; ++elem )
    {
      fprintf( fp, " %e", injPar[elem] );
      if ( elem == distElem )
        fprintf( fplog, "\t%e", injPar[elem] / MPC );
      else
        fprintf( fplog, "\t%e", injPar[elem] );
    }

    memcpy( this_sim_insp->waveform, waveform, 
        sizeof(CHAR) * LIGOMETA_WAVEFORM_MAX );
    this_sim_insp->mtotal = injPar[mTotElem];
    this_sim_insp->eta = injPar[etaElem];
    this_sim_insp->distance = injPar[distElem] / MPC;
    this_sim_insp->longitude = injPar[lonElem];
    this_sim_insp->latitude = injPar[latElem];
    this_sim_insp->inclination = injPar[incElem];
    this_sim_insp->coa_phase = injPar[phiElem];
    this_sim_insp->polarization = injPar[psiElem];

    this_sim_insp->eff_dist_h = deffH = eff_dist( nxH, nyH, injPar, gmst )/MPC;
    this_sim_insp->eff_dist_l = deffL = eff_dist( nxL, nyL, injPar, gmst )/MPC;
    site_time_of_day( &todH, &todL, tsec, tnan );
    fprintf( fplog, "\t%e\t%e\t%e\t%e\n", todH, todL, deffH, deffL );

    if ( inj < ninj - 1 )
    {
      this_sim_insp = this_sim_insp->next = (SimInspiralTable *)
        LALCalloc( 1, sizeof(SimInspiralTable) );
    }
  }

  fprintf( fp, "</real_4>\n" );

  fputs( "</ilwd>\n", fp );
  fclose( fp );
  fclose( fplog );

  memset( &xmlfp, 0, sizeof(LIGOLwXMLStream) );
  snprintf( fname, sizeof(fname), "injections-%d.xml", rand_seed );
  LAL_CALL( LALOpenLIGOLwXMLFile( &status, &xmlfp, fname), &status );

  LAL_CALL( LALGPSTimeNow ( &status, &(proctable.processTable->end_time),
        &accuracy ), &status );
  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, process_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, proctable, 
        process_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );

  LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, process_params_table ), 
      &status );
  LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, procparams, 
        process_params_table ), &status );
  LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );

  if ( injections.simInspiralTable )
  {
    LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, sim_inspiral_table ), 
        &status );
    LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, injections, 
          sim_inspiral_table ), &status );
    LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );
  }
  
  LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlfp ), &status );

  return 0;
}
