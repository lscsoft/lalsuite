#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lalapps.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Random.h>

RCSID( "$Id$" );

#define usgfmt \
  "Usage: %s [options]\n" \
  "Options [default in brackets]:\n" \
  "  -h            print this message\n" \
  "  -V            print version info\n" \
  "  -v            verbose\n" \
  "  -d dbglvl     set debug level to dbglvl [0]\n" \
  "  -m min-max    minimum and maximum binary masses (solar) [1.3-1.5]\n" \
  "  -r range      maximum range of binary (kpc) [1.0]\n" \
  "  -s seed       set random number generator seed [0 = from clock]\n" \
  "  -t start-end  inject events between GPS times start and end\n" \
  "                [600000000-600010240]\n" \
  "  -u step       use uniform timestep step (seconds) (Poisson if negative)\n"\
  "                [31.8309886183791]\n"

#define usage( program ) fprintf( stderr, usgfmt, program )
#define UNITS "msun,none,m,rad,rad,rad,rad,rad"
struct time_list { UINT8 tinj; struct time_list *next; };
enum { mTotElem, etaElem, distElem, incElem, phiElem, lonElem, latElem,
  psiElem, numElem };

extern char *optarg;
extern int optind, opterr, optopt;
extern int vrbflg;

int main( int argc, char *argv[] )
{
  const char *program = argv[0];
  const char *dbglvl  = NULL;
  struct time_list *tlisthead;
  struct time_list *tlistelem;
  static RandomParams *rpar;
  REAL4 injPar[numElem];
  LALStatus status   = blank_status;
  REAL4 meanTimeStep = 1e2 / LAL_PI;
  REAL4 maxMass      = 1.5;
  REAL4 minMass      = 1.3;
  REAL4 maxDist      = 1e3 * LAL_PC_SI;
  UINT4 seed         = 0;
  UINT4 duration     = 10240;
  UINT4 startTime    = 600000000;
  UINT4 stopTime     = startTime + duration;
  UINT8 tinj;
  UINT4 ninj;
  UINT4 inj;
  int   opt;
  FILE *fp;

  /* parse options */
  while ( 0 < ( opt = getopt( argc, argv, "hVvd:m:r:s:t:u:" ) ) )
  {
    switch ( opt )
    {
      case 'h':
        usage( program );
        return 0;
      case 'V':
        PRINT_VERSION( "hello" );
        return 0;
      case 'v':
        vrbflg = 1;
        break;
      case 'd':
        dbglvl = optarg;
        break;
      case 'm':
        minMass = atof( optarg );
        maxMass = atof( ( optarg = strchr( optarg, '-' ) ) ? ++optarg : "" );
        if ( maxMass < minMass || minMass <= 0 )
        {
          fprintf( stderr, "invalid mass range %f-%f\n", minMass, maxMass );
          return 1;
        }
        break;
      case 'r':
        maxDist = atof( optarg ) * 1e3 * LAL_PC_SI;
        break;
      case 's':
        seed = atoi( optarg );
        break;
      case 't':
        startTime = atoi( optarg );
        stopTime = atoi( ( optarg = strchr( optarg, '-' ) ) ? ++optarg : "" );
        if ( stopTime <= startTime )
          stopTime = startTime + duration;
        else
          duration = stopTime - startTime;
        break;
      case 'u':
        meanTimeStep = atof( optarg );
        break;
      default:
        usage( program );
        return 1;
    }
  }
  if ( optind < argc )
  {
    usage( program );
    return 1;
  }

  fprintf( stderr, "maximum distance: %g m\n", maxDist );
  fprintf( stderr, "mass range: [%g,%g] msun\n", minMass, maxMass );
  fprintf( stderr, "time range: [%d,%d] s\n", startTime, stopTime );
  fprintf( stderr, "mean time step: %g s\n", fabs( meanTimeStep ) );
  if ( meanTimeStep < 0 )
    fprintf( stderr, "poisson random time intervals\n" );
  else
    fprintf( stderr, "constant time intervals\n" );


  /* set debug level */
  set_debug_level( dbglvl );

  LAL_CALL( LALCreateRandomParams( &status, &rpar, seed ), &status );

  /* make list of random times */
  ninj = 1;
  tinj = (INT8)(1000000000) * (INT8)startTime;
  tlistelem = tlisthead = LALCalloc( 1, sizeof( *tlisthead ) );
  tlistelem->tinj = tinj;
  while ( 1 )
  {
    if ( meanTimeStep < 0 ) /* Poisson */
    {
      REAL4 udev;
      LAL_CALL( LALUniformDeviate( &status, &udev, rpar ), &status );
      tinj += (UINT8)( 1e9 * meanTimeStep * log( udev ) );
    }
    else /* uniform */
    {
      tinj += (UINT8)( 1e9 * meanTimeStep );
    }
    if ( tinj > (UINT8)(1000000000) * (UINT8)(stopTime) )
      break;
    tlistelem = tlistelem->next = LALCalloc( 1, sizeof( *tlistelem ) );
    tlistelem->tinj = tinj;
    ++ninj;
  }
  
  /* output first sequence */

  fp = fopen( "injepochs.ilwd", "w" );
  fputs( "<?ilwd?>\n", fp );
  fputs( "<ilwd name='injepochs::sequence' size='7'>\n", fp );

  fprintf( fp, "\t<lstring name='real:domain' size='4'>TIME</lstring>\n" );
  fprintf( fp, "\t<int_4u name='gps_sec:start_time' units='sec'>%d</int_4u>\n",
      startTime );
  fprintf( fp, "\t<int_4u name='gps_nan:start_time' units='nanosec'>0</int_4u>\n" );
  fprintf( fp, "\t<int_4u name='gps_sec:stop_time' units='sec'>%d</int_4u>\n",
      stopTime );
  fprintf( fp, "\t<int_4u name='gps_nan:stop_time' units='nanosec'>0</int_4u>\n" );
  fprintf( fp, "\t\t<real_8 name='time:step_size' units='sec'>%e</real_8>\n",
      meanTimeStep );

  fprintf( fp, "\t<int_4u name='data' ndim='2' dims='2,%d' units='s,ns'>", ninj );
  fprintf( fp, "%d 0", startTime );
  tlistelem = tlisthead->next;
  for ( inj = 1; inj < ninj; ++inj )
  {
    UINT4 tsec = (UINT4)( tlistelem->tinj / 1000000000 );
    UINT4 tnan = (UINT4)( tlistelem->tinj % 1000000000 );
    fprintf( fp, " %d %d", tsec, tnan );
    tlistelem = tlistelem->next;
  }
  fprintf( fp, "</int_4u>\n" );

  fputs( "</ilwd>\n", fp );
  fclose( fp );

  /* output second sequence */
  fp = fopen( "injparams.ilwd", "w" );
  fputs( "<?ilwd?>\n", fp );
  fputs( "<ilwd name='injparams::sequence' size='7'>\n", fp );

  fprintf( fp, "\t<lstring name='real:domain' size='4'>TIME</lstring>\n" );
  fprintf( fp, "\t<int_4u name='gps_sec:start_time' units='sec'>%d</int_4u>\n",
      startTime );
  fprintf( fp, "\t<int_4u name='gps_nan:start_time' units='nanosec'>0</int_4u>\n" );
  fprintf( fp, "\t<int_4u name='gps_sec:stop_time' units='sec'>%d</int_4u>\n",
      stopTime );
  fprintf( fp, "\t<int_4u name='gps_nan:stop_time' units='nanosec'>0</int_4u>\n" );
  fprintf( fp, "\t<real_8 name='time:step_size' units='sec'>%e</real_8>\n",
      meanTimeStep );

  fprintf( fp, "\t<real_4 name='data' ndim='2' dims='%d,%d' units='" UNITS "'>",
      numElem, ninj );
  for ( inj = 0; inj < ninj; ++inj )
  {
    REAL4 u1, u2, u3, u4, u5;
    REAL4 x, y, z, r, rsq;
    REAL4 m1, m2;
    INT4 elem;

    do
    {
      LAL_CALL( LALUniformDeviate( &status, &u1, rpar), &status );
      LAL_CALL( LALUniformDeviate( &status, &u2, rpar), &status );
      LAL_CALL( LALUniformDeviate( &status, &u3, rpar), &status );
      x = 2 * u1 - 1;
      y = 2 * u2 - 1;
      z = 2 * u3 - 1;
      rsq = x * x + y * y + z * z;
    }
    while ( rsq > 1 );
    r = sqrt( rsq );

    LAL_CALL( LALUniformDeviate( &status, &u1, rpar), &status );
    LAL_CALL( LALUniformDeviate( &status, &u2, rpar), &status );
    LAL_CALL( LALUniformDeviate( &status, &u3, rpar), &status );
    LAL_CALL( LALUniformDeviate( &status, &u4, rpar), &status );
    LAL_CALL( LALUniformDeviate( &status, &u5, rpar), &status );

    m1 = minMass + ( maxMass - minMass ) * u1;
    m2 = minMass + ( maxMass - minMass ) * u2;

    injPar[mTotElem] = m1 + m2;
    injPar[etaElem]  = m1 * m2 / ( ( m1 + m2 ) * ( m1 + m2 ) );
    injPar[distElem] = r * maxDist;
    injPar[incElem]  = -0.5 * LAL_PI + LAL_PI * u3; 
    injPar[phiElem]  = 2 * LAL_PI * u4;
    injPar[lonElem]  = atan2( y, x );
    injPar[latElem]  = atan( z / r );
    injPar[psiElem]  = 2 * LAL_PI * u5;
#if 0
    injPar[incElem]  = 0;
    injPar[phiElem]  = 0;
    injPar[lonElem]  = 0;
    injPar[latElem]  = 0;
    injPar[psiElem]  = 0;
#endif

    fprintf( fp, "%s%e", inj ? " " : "", injPar[0] );
    for ( elem = 1; elem < numElem; ++elem )
      fprintf( fp, " %e", injPar[elem] );
  }
  fprintf( fp, "</real_4>\n" );

  fputs( "</ilwd>\n", fp );
  fclose( fp );  

  /* free memory */
  while ( tlisthead )
  {
    tlistelem = tlisthead;
    tlisthead = tlistelem->next;
    LALFree( tlistelem );
  }
    
  LAL_CALL( LALDestroyRandomParams( &status, &rpar ), &status );

  LALCheckMemoryLeaks();
  return 0;
}
