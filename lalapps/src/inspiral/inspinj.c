/*----------------------------------------------------------------------- 
 * 
 * File Name: inspinj.c
 *
 * Author: Brown, D. A. and Crieghton, J. D. E.
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
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXML.h>
#include <lal/Date.h>
#include <lal/TimeDelay.h>

#define USAGE \
  "lalapps_inspinj [options]\n"\
"\nDefaults are shown in brackets\n\n" \
"  --help                   display this message\n"\
"  --source-file FILE       read source parameters from FILE\n"\
"  --mass-file FILE         read population mass parameters from FILE\n"\
"  --gps-start-time TIME    start injections at GPS time TIME (729273613)\n"\
"  --gps-end-time TIME      end injections at GPS time TIME (734367613)\n"\
"  --time-step STEP         space injections by ave of STEP sec (2630 / PI)\n"\
"  --time-interval TIME     distribute injections in interval TIME (0)\n"\
"  --seed SEED              seed random number generator with SEED (1)\n"\
"  --waveform NAME          set waveform type to NAME (GeneratePPNtwoPN)\n"\
"  --user-tag STRING        set the usertag to STRING\n"\
"  --tama-output            generate a text file for tama\n"\
"  --write-eff-dist         output the effective distances in tama file\n"\
"  --ilwd                   generate an ILWD file for LDAS\n"\
"\n"

RCSID( "$Id$" );

#define KPC ( 1e3 * LAL_PC_SI )
#define MPC ( 1e6 * LAL_PC_SI )
#define GPC ( 1e9 * LAL_PC_SI )

#define UNITS "msun,none,m,rad,rad,rad,rad,rad"

#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "inspinj"

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

enum { mTotElem, etaElem, distElem, incElem, phiElem, lonElem, latElem,
  psiElem, m1Elem, m2Elem, numElem };

SimInspiralTable *this_sim_insp;

const int randm = 2147483647;
struct { int i; int y; int v[32]; } randpar;

char *massFileName = NULL;
char *sourceFileName = NULL;

int num_source;
struct {
  char   name[16];
  double ra;
  double dec;
  double dist;
  double lum;
  double fudge;
} *source_data;

struct time_list { long long tinj; struct time_list *next; };


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


/*
 *
 * compute effective distance of injection
 *
 */


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

  theta = 0.5 * LAL_PI - injPar[latElem];
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


/* 
 *
 * convert galactic coords to equatorial coords (all angles in radians) 
 *
 */


int galactic_to_equatorial( double *alpha, double *delta, double l, double b )
{
  const double alphaNGP = 192.8594813 * LAL_PI / 180.0;
  const double deltaNGP = 27.1282511 * LAL_PI / 180.0;
  const double lascend  = 33.0 * LAL_PI / 180.0;
  double lm = l - lascend;

  *alpha  = atan2( cos(b)*cos(lm),
      sin(b)*cos(deltaNGP) - cos(b)*sin(deltaNGP)*sin(lm) );
  *alpha += alphaNGP;
  *delta  = asin( cos(b)*cos(deltaNGP)*sin(lm) + sin(b)*sin(deltaNGP) );

  return 0;
}


/* 
 *
 * generate a sky position for a random Galactic inspiral 
 *
 */


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
  {
    z = -h_scale * log( 2 * ( 1 - u ) );
  }
  else
  {
    z = h_scale * log( 2 * u );
  }
  phi = 2 * LAL_PI * my_urandom();

  rho2 = r_sun * r_sun + r * r - 2 * r_sun * r * cos( phi );
  *dist = sqrt( z * z + rho2 );
  l = atan2( r * sin( phi ), r_sun - r * cos( phi ) );
  b = asin( z / (*dist) );
  galactic_to_equatorial( alpha, delta, l, b );

  return 0;
}


/*
 *
 * functions to read source masses and distribution
 *
 */

int read_source_mass_data( double **pm1, double **pm2 )
{
  const char *basename = massFileName ? massFileName : "BNSMasses.dat";
  char fname[256];
  char line[256];
  FILE   *fp;
  double *m1;
  double *m2;
  int n = 0;

  LALSnprintf( fname, sizeof( fname ), basename );
  if ( ! massFileName || ! ( fp = fopen( massFileName, "r" ) ) )
  {
    if ( *basename != '.' && *basename != '/' ) /* not abs or rel path */
    { /* prepend path from env variable */
      const char *path = getenv( "LALAPPS_DATA_PATH" );
      LALSnprintf( fname, sizeof( fname ), "%s/%s",
	  path ? path : PREFIX "/" PACKAGE "/share", basename );
    }
    fp = fopen( fname, "r" );
  }

  if ( ! fp )
  {
    perror( "read_source_mass_data" );
    fprintf( stderr, "Could not find file %s\n", fname );
    fprintf( stderr, "Set environment LALAPPS_DATA_PATH to location of file %s\n", basename );
    exit( 1 );
  }

  while ( fgets( line, sizeof( line ), fp ) )
    ++n;
  m1 = *pm1 = calloc( n, sizeof( *m1 ) );
  if ( ! m1 )
  {
    fprintf( stderr, "alloc error\n" );
    exit( 1 );
  }
  m2 = *pm2 = calloc( n, sizeof( *m2 ) );
  if ( ! m2 )
  {
    fprintf( stderr, "alloc error\n" );
    exit( 1 );
  }
  rewind( fp );
  while ( fgets( line, sizeof( line ), fp ) )
  {
    sscanf( line, "%le %le", m1++, m2++ );
  }
  fclose( fp );
  return n;
}

int read_source_data( void )
{
  const char *basename = sourceFileName ? sourceFileName : "inspsrcs.dat";
  char fname[256];
  char line[256];
  FILE *fp;
  int i;

  LALSnprintf( fname, sizeof( fname ), sourceFileName );
  if ( ! sourceFileName || ! ( fp = fopen( sourceFileName, "r" ) ) )
  {
    if ( *basename != '.' && *basename != '/' ) /* not abs or rel path */
    { /* prepend path from env variable */
      const char *path = getenv( "LALAPPS_DATA_PATH" );
      LALSnprintf( fname, sizeof( fname ), "%s/%s",
	  path ? path : PREFIX "/" PACKAGE "/share", basename );
    }
    fp = fopen( fname, "r" );
  }

  if ( ! fp )
  {
    perror( "read_source_data" );
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
	fprintf( stderr, "error parsing source datafile %s\n", sourceFileName );
	exit( 1 );
      }

      /* by convention, overall sign is carried only on hours/degrees entry */
      source_data[i].ra  = ( ra_h + ra_m / 60.0 ) * LAL_PI / 12.0;
      source_data[i].dec = ( dec_d + dec_m / 60.0 ) * LAL_PI / 180.0;

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


/* 
 *
 * generate a sky position for a random inspiral from the Galaxy, LMC, or SMC 
 *
 */


int sky_position( double *dist, double *alpha, double *delta, char *source )
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
    num_source = read_source_data();
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
      LALSnprintf( source, sizeof(source)/sizeof(*source), 
	  "%s", source_data[i].name );
      *dist  = source_data[i].dist;
      *alpha = source_data[i].ra;
      *delta = source_data[i].dec;
      return 0;
    }

  /* galactic event */
  LALSnprintf( source, sizeof(source)/sizeof(*source), "MW" );
  return galactic_sky_position( dist, alpha, delta );
}


/* 
 *
 * generate all parameters (sky position and angles) for a random inspiral 
 *
 */


int inj_params( double *injPar, char *source )
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
    n = read_source_mass_data( &m1arr, &m2arr );

  sky_position( &dist, &alpha, &delta, source );
  i = (size_t)( n * my_urandom() );

  injPar[m1Elem] = m1 = m1arr[i];
  injPar[m2Elem] = m2 = m2arr[i];
  injPar[mTotElem] = m1 + m2;
  injPar[etaElem]  = m1 * m2 / ( ( m1 + m2 ) * ( m1 + m2 ) );
  injPar[incElem]  = acos( -1.0 + 2.0 * my_urandom() );
  injPar[phiElem]  = 2 * LAL_PI * my_urandom();
  injPar[psiElem]  = 2 * LAL_PI * my_urandom();
  injPar[distElem] = dist;
  injPar[lonElem]  = alpha;
  injPar[latElem]  = delta;

  return 0;
}

int main( int argc, char *argv[] )
{
  double nxH[3] = { -0.2239, +0.7998, +0.5569 };
  double nyH[3] = { -0.9140, +0.0261, -0.4049 };
  double nxL[3] = { -0.9546, -0.1416, -0.2622 };
  double nyL[3] = { +0.2977, -0.4879, -0.8205 };
  double nxG[3] = { -0.6261, -0.5522, +0.5506 };
  double nyG[3] = { -0.4453, +0.8665, +0.2255 };
  double nxT[3] = { +0.6490, +0.7608, +0.0000 };
  double nyT[3] = { -0.4437, +0.3785, -0.8123 };
  double nxV[3] = { -0.7005, +0.2085, +0.6826 };
  double nyV[3] = { -0.0538, -0.9691, +0.2408 };
  const long S2StartTime   = 729273613;  /* Feb 14 2003 16:00:00 UTC */
  const long S2StopTime    = 734367613;  /* Apr 14 2003 15:00:00 UTC */
  long gpsStartTime = S2StartTime;
  long gpsEndTime = S2StopTime;
  double meanTimeStep = 2630 / LAL_PI; /* seconds between injections     */
  double timeInterval = 0;

  long long tinj = 1000000000LL * gpsStartTime;
  struct time_list  tlisthead;
  struct time_list *tlistelem = &tlisthead;

  double injPar[numElem];
  size_t ninj;
  size_t inj;
  FILE *fp = NULL;
  FILE *fq = NULL;
  int rand_seed = 1;
  static int ilwd = 0;
  int tamaOutput = 0;
  int writeDeff = 0;

  /* waveform */
  CHAR waveform[LIGOMETA_WAVEFORM_MAX];

  /* xml output data */
  CHAR                  fname[256];
  CHAR                 *userTag = NULL;
  LALStatus             status = blank_status;
  LALLeapSecAccuracy    accuracy = LALLEAPSEC_LOOSE;
  MetadataTable         proctable;
  MetadataTable         procparams;
  MetadataTable         injections;
  ProcessParamsTable   *this_proc_param;
  LIGOLwXMLStream       xmlfp;
  CHAR		       tamaFileName[256];

  /* getopt arguments */
  struct option long_options[] =
  {
    {"help",                          no_argument, 0,                'h'},
    {"source-file",             required_argument, 0,                'f'},
    {"mass-file",               required_argument, 0,                'm'},
    {"gps-start-time",          required_argument, 0,                'a'},
    {"gps-end-time",            required_argument, 0,                'b'},
    {"time-step",               required_argument, 0,                't'},
    {"time-interval",		required_argument, 0,		     'i'},
    {"seed",                    required_argument, 0,                's'},
    {"waveform",                required_argument, 0,                'w'},
    {"user-tag",                required_argument, 0,                'Z'},
    {"userTag",                 required_argument, 0,                'Z'},
    {"tama-output",		      no_argument, &tamaOutput,	      1 },
    {"write-eff-dist",                no_argument, &writeDeff,        1 },
    {"ilwd",                          no_argument, &ilwd,             1 },
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

  /* clear the waveform field */
  memset( waveform, 0, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR) );

  /* parse the arguments */
  while ( 1 )
  {
    /* getopt_long stores long option here */
    int option_index = 0;
    long int gpsinput;
    size_t optarg_len;

    c = getopt_long_only( argc, argv, 
	"hf:m:a:b:t:s:w:i:", long_options, &option_index );

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
	optarg_len = strlen( optarg ) + 1;
	sourceFileName = calloc( 1, optarg_len * sizeof(char) );
	memcpy( sourceFileName, optarg, optarg_len * sizeof(char) );
	this_proc_param = this_proc_param->next = 
	  next_process_param( long_options[option_index].name, "string", 
	      "%s", optarg );
	break;

      case 'm':
	optarg_len = strlen( optarg ) + 1;
	massFileName = calloc( 1, optarg_len * sizeof(char) );
	memcpy( massFileName, optarg, optarg_len * sizeof(char) );
	this_proc_param = this_proc_param->next = 
	  next_process_param( long_options[option_index].name, "string", 
	      "%s", optarg );
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
	gpsEndTime = gpsinput;
	this_proc_param = this_proc_param->next = 
	  next_process_param( long_options[option_index].name, "int", 
	      "%ld", gpsinput );
	break;

      case 's':
	rand_seed = atoi( optarg );
	this_proc_param = this_proc_param->next = 
	  next_process_param( long_options[option_index].name, "int", 
	      "%d", rand_seed );
	break;

      case 't':
	{
	  double tstep = atof( optarg );
	  meanTimeStep = tstep ;
	  this_proc_param = this_proc_param->next = 
	    next_process_param( long_options[option_index].name, "float", 
		"%le", tstep );
	}
	break;

      case 'i':
	timeInterval = 1000000000LL * atof( optarg );
	this_proc_param = this_proc_param->next = 
	  next_process_param( long_options[option_index].name, 
	      "float", "%le", timeInterval/1000000000 );
	break;

      case 'w':
	LALSnprintf( waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), "%s",
	    optarg );
	this_proc_param = this_proc_param->next = 
	  next_process_param( long_options[option_index].name, "string", 
	      "%s", optarg );

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

  seed_random( rand_seed );

  tlisthead.tinj = tinj + timeInterval * my_urandom();
  tlisthead.next = NULL;

  if ( ! *waveform )
  {
    /* default to Tev's GeneratePPNInspiral as used in */
    LALSnprintf( waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR), 
	"GeneratePPNtwoPN" );
  }

  /* store the tamaOutput argument */
  if ( tamaOutput )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", 
	PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, "--tama-output" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, " " );
  }

  /* store the writeDeff argument */
  if ( writeDeff )
  {
    this_proc_param = this_proc_param->next = (ProcessParamsTable *)
      calloc( 1, sizeof(ProcessParamsTable) );
    LALSnprintf( this_proc_param->program, LIGOMETA_PROGRAM_MAX, "%s", 
	PROGRAM_NAME );
    LALSnprintf( this_proc_param->param, LIGOMETA_PARAM_MAX, 
        "--write-eff-dist" );
    LALSnprintf( this_proc_param->type, LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( this_proc_param->value, LIGOMETA_VALUE_MAX, " " );
  }

  /* store the ilwd argument */
  if ( ilwd )
  {
    LALSnprintf( procparams.processParamsTable->program, 
	LIGOMETA_PROGRAM_MAX, "%s", PROGRAM_NAME );
    LALSnprintf( procparams.processParamsTable->param,
	LIGOMETA_PARAM_MAX, "--ilwd" );
    LALSnprintf( procparams.processParamsTable->type, 
	LIGOMETA_TYPE_MAX, "string" );
    LALSnprintf( procparams.processParamsTable->value, 
	LIGOMETA_TYPE_MAX, " " );
  }
  else
  {
    this_proc_param = procparams.processParamsTable;
    procparams.processParamsTable = procparams.processParamsTable->next;
    free( this_proc_param );
  }

  /* create the first injection */
  this_sim_insp = injections.simInspiralTable = (SimInspiralTable *)
    calloc( 1, sizeof(SimInspiralTable) );

  /* make injection times at intervals of TSTEP seconds, 
   * add random time between 0 and timeInterval */
  ninj = 1;
  tlistelem = &tlisthead;
  while ( 1 )
  {
    tinj += (long long)( 1e9 * meanTimeStep );
    if ( tinj > 1000000000LL * gpsEndTime )
      break;
    tlistelem = tlistelem->next = calloc( 1, sizeof( *tlistelem ) );
    tlistelem->tinj = tinj + timeInterval * my_urandom();
    ++ninj;
  }


  /*
   *
   * First Sequence: injection epochs.
   *
   */


  if ( ilwd )
  {
    fp = fopen( "injepochs.ilwd", "w" );
    fputs( "<?ilwd?>\n", fp );
    fputs( "<ilwd name='injepochs::sequence' size='7'>\n", fp );

    fprintf( fp, 
	"\t<lstring name='real:domain' size='4'>TIME</lstring>\n" );
    fprintf( fp, 
	"\t<int_4u name='gps_sec:start_time' units='sec'>%ld</int_4u>\n",
	gpsStartTime );
    fprintf( fp, 
	"\t<int_4u name='gps_nan:start_time' units='nanosec'>0</int_4u>\n" );
    fprintf( fp, 
	"\t<int_4u name='gps_sec:stop_time' units='sec'>%ld</int_4u>\n",
	gpsEndTime );
    fprintf( fp, 
	"\t<int_4u name='gps_nan:stop_time' units='nanosec'>0</int_4u>\n");
    fprintf( fp,
	"\t\t<real_8 name='time:step_size' units='sec'>%e</real_8>\n",
	meanTimeStep );
    fprintf( fp, "\t<int_4u ndim='2' dims='2,%d' name='data' units='s,ns'>",
	(int) ninj );
    fprintf( fp, "%ld 0", gpsStartTime );
  }

  tlistelem = tlisthead.next;

  for ( inj = 1; inj < ninj; ++inj )
  {
    long tsec = (long)( tlistelem->tinj / 1000000000LL );
    long tnan = (long)( tlistelem->tinj % 1000000000LL );
    if ( ilwd ) fprintf( fp, " %ld %ld", tsec, tnan );
    tlistelem = tlistelem->next;
  }

  if ( ilwd ) 
  {
    fprintf( fp, "</int_4u>\n" );
    fputs( "</ilwd>\n", fp );
    fclose( fp );
  }


  /*
   * 
   * Second sequence: injection parameters.
   * 
   */


  if ( ilwd )
  {
    fp = fopen( "injparams.ilwd", "w" );
    fputs( "<?ilwd?>\n", fp );
    fputs( "<ilwd name='injparams::sequence' size='7'>\n", fp );

    fprintf( fp, "\t<lstring name='real:domain' size='4'>TIME</lstring>\n" );
    fprintf( fp, 
	"\t<int_4u name='gps_sec:start_time' units='sec'>%ld</int_4u>\n",
	gpsStartTime );
    fprintf( fp, 
	"\t<int_4u name='gps_nan:start_time' units='nanosec'>0</int_4u>\n" );
    fprintf( fp, 
	"\t<int_4u name='gps_sec:stop_time' units='sec'>%ld</int_4u>\n",
	gpsEndTime );
    fprintf( fp, 
	"\t<int_4u name='gps_nan:stop_time' units='nanosec'>0</int_4u>\n" );
    fprintf( fp, 
	"\t<real_8 name='time:step_size' units='sec'>%e</real_8>\n",
	meanTimeStep );

    fprintf( fp, 
	"\t<real_4 ndim='2' dims='%d,%d' name='data' units='" UNITS "'>",
	numElem, (int) ninj );
  }

  tlistelem = &tlisthead;

  for ( inj = 0; inj < ninj; ++inj )
  {
    int elem;
    long tsec = this_sim_insp->geocent_end_time.gpsSeconds = 
      (long)( tlistelem->tinj / 1000000000LL );
    long tnan = this_sim_insp->geocent_end_time.gpsNanoSeconds = 
      (long)( tlistelem->tinj % 1000000000LL );
    double gmst;
    
    /* get gmst (radians) */
    gmst =  greenwich_mean_sidereal_time( tsec, tnan, 32 );

    /* save gmst (hours) in sim_inspiral table */
    this_sim_insp->end_time_gmst = gmst * 12.0 / LAL_PI;

    tlistelem = tlistelem->next;

    inj_params( injPar, this_sim_insp->source );

    if ( ilwd ) fprintf( fp, "%s%e", inj ? " " : "", injPar[0] );

    for ( elem = 1; elem < numElem; ++elem )
    {
      if ( ilwd )
      {
	fprintf( fp, " %e", injPar[elem] );
      }
    }

    memcpy( this_sim_insp->waveform, waveform, 
	sizeof(CHAR) * LIGOMETA_WAVEFORM_MAX );
    this_sim_insp->mass1 = injPar[m1Elem];
    this_sim_insp->mass2 = injPar[m2Elem];
    this_sim_insp->eta = injPar[etaElem];
    this_sim_insp->distance = injPar[distElem] / MPC;
    this_sim_insp->longitude = injPar[lonElem];
    this_sim_insp->latitude = injPar[latElem];
    this_sim_insp->inclination = injPar[incElem];
    this_sim_insp->coa_phase = injPar[phiElem];
    this_sim_insp->polarization = injPar[psiElem];

    /* populate the site specific information */
    LAL_CALL(LALPopulateSimInspiralSiteInfo( &status, this_sim_insp ), 
	&status);

    /* inspinj has an independent way of calculating the effective distances.
       Calculate them instead using this */
    /* XXX Need to implement check that inspinj answers agree with LAL XXX */
    
    /* Hanford */
    this_sim_insp->eff_dist_h = eff_dist( nxH, nyH, injPar, gmst )/MPC;

    /* Livingston */
    this_sim_insp->eff_dist_l = eff_dist( nxL, nyL, injPar, gmst )/MPC;

    /* GEO */
    this_sim_insp->eff_dist_g = eff_dist( nxG, nyG, injPar, gmst )/MPC;

    /* TAMA */
    this_sim_insp->eff_dist_t = eff_dist( nxT, nyT, injPar, gmst )/MPC;

    /* Virgo */
    this_sim_insp->eff_dist_v = eff_dist( nxV, nyV, injPar, gmst )/MPC;
    
    if ( inj < ninj - 1 )
    {
      this_sim_insp = this_sim_insp->next = (SimInspiralTable *)
	calloc( 1, sizeof(SimInspiralTable) );
    }
  }

  if ( ilwd )
  {
    fprintf( fp, "</real_4>\n" );
    fputs( "</ilwd>\n", fp );
    fclose( fp );
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

  if ( ! ilwd )
  {
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

    if ( injections.simInspiralTable )
    {
      LAL_CALL( LALBeginLIGOLwXMLTable( &status, &xmlfp, sim_inspiral_table ), 
	  &status );
      LAL_CALL( LALWriteLIGOLwXMLTable( &status, &xmlfp, injections, 
	    sim_inspiral_table ), &status );
      LAL_CALL( LALEndLIGOLwXMLTable ( &status, &xmlfp ), &status );

      /* write out text file for Tama */
      if ( tamaOutput )
      {
	if ( userTag )
	{
	  LALSnprintf( tamaFileName, sizeof(tamaFileName), 
	      "HLT-INJECTIONS_%d_%s-%d-%d.txt", rand_seed, userTag, 
	      gpsStartTime, gpsEndTime - gpsStartTime );
	}
	else
	{
	  LALSnprintf( tamaFileName, sizeof(tamaFileName), 
	      "HLT-INJECTIONS_%d-%d-%d.txt", rand_seed, gpsStartTime, 
	      gpsEndTime - gpsStartTime );
	}
	fq = fopen( tamaFileName, "w" );
	fprintf( fq, "Injections for joint LIGO-TAMA analysis \n" );
	fprintf( fq, "geocentric end time  hanford end time   " );
	fprintf( fq, "livingston end time   tama end time     " );
        if ( writeDeff )
        {
          fprintf( fq, "H d_eff(kpc) L d_eff(kpc) T d_eff(kpc) " );
	}
        fprintf( fq, "  mtotal        eta      distance(kpc) " );
	fprintf( fq, " longitude     latitude    " );
	fprintf( fq, "inclination   coa_phase   polarization\n");

	this_sim_insp = injections.simInspiralTable;
	while ( this_sim_insp )
	{
	  REAL8 injtime, ltime, htime, ttime;
	  LAL_CALL( LALGPStoFloat( &status, &injtime, 
		&(this_sim_insp->geocent_end_time) ),  &status );
	  LAL_CALL( LALGPStoFloat( &status, &htime, 
		&(this_sim_insp->h_end_time) ),  &status );
	  LAL_CALL( LALGPStoFloat( &status, &ltime, 
		&(this_sim_insp->l_end_time) ),  &status );
	  LAL_CALL( LALGPStoFloat( &status, &ttime, 
		&(this_sim_insp->t_end_time) ),  &status );

	  fprintf( fq, "%19.9f %19.9f %19.9f %19.9f ", injtime, htime, ltime,
	      ttime );
          if ( writeDeff )
          {
            fprintf( fq, "%12.6e %12.6e %12.6e ", 
                1.0e+03 * this_sim_insp->eff_dist_h,
                1.0e+03 * this_sim_insp->eff_dist_l, 
                1.0e+03 * this_sim_insp->eff_dist_t );
	  }
          fprintf( fq, "%12.6e %12.6e %12.6e %12.6e %13.6e ",
	      (this_sim_insp->mass1 + this_sim_insp->mass2), 
	      this_sim_insp->eta, 1.0e+03 * this_sim_insp->distance,
	      this_sim_insp->longitude, this_sim_insp->latitude );
	  fprintf( fq, "%12.6e %12.6e %12.6e\n", 
	      this_sim_insp->inclination, this_sim_insp->coa_phase,
	      this_sim_insp->polarization );
	  this_sim_insp = this_sim_insp->next;
	}
      }

    }

    LAL_CALL( LALCloseLIGOLwXMLFile ( &status, &xmlfp ), &status );
  }

  return 0;
}
