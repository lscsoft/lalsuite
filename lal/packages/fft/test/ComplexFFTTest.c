/*----------------------------------------------------------------------- 
 * 
 * File Name: ComplexFFTTest.c
 * 
 * Author: Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LALConfig.h"

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include "LALStdlib.h"
#include "AVFactories.h"
#include "ComplexFFT.h"

#define CODES_(x) #x
#define CODES(x) CODES_(x)

NRCSID( MAIN, "$Id$" );

extern char *optarg;
extern int   optind;

int debuglevel = 0;
int verbose    = 0;

static void
Usage( const char *program, int exitflag );

static void
ParseOptions( int argc, char *argv[] );

static void
TestStatus( Status *status, const char *expectedCodes, int exitCode );

static void
ClearStatus( Status *status );

static void
CheckErrorCodes();

int
main( int argc, char *argv[] )
{
  const UINT4 n       = 17;
  const REAL4 epsilon = 1e-6;

  Status          stat = {};
  ComplexFFTPlan *pfwd = NULL;
  ComplexFFTPlan *pinv = NULL;
  COMPLEX8Vector *avec = NULL;
  COMPLEX8Vector *bvec = NULL;
  COMPLEX8Vector *cvec = NULL;

  FILE *fp;
  INT4  i;

  ParseOptions( argc, argv );

  fp = verbose ? stdout : NULL ;

  CCreateVector( &stat, &avec, n );
  TestStatus( &stat, CODES( 0 ), 1 );

  CCreateVector( &stat, &bvec, n );
  TestStatus( &stat, CODES( 0 ), 1 );

  CCreateVector( &stat, &cvec, n );
  TestStatus( &stat, CODES( 0 ), 1 );

  EstimateFwdComplexFFTPlan( &stat, &pfwd, n );
  TestStatus( &stat, CODES( 0 ), 1 );

  EstimateInvComplexFFTPlan( &stat, &pinv, n );
  TestStatus( &stat, CODES( 0 ), 1 );


  for ( i = 0; i < n; ++i )
  {
    avec->data[i].im = i % 5 - 2;
    avec->data[i].re = i % 3 - 1;
    fp ? fprintf( fp, "%+.0f\t%+.0f\n", avec->data[i].re, avec->data[i].im )
       : 0;
  }
  fp ? fprintf( fp, "\n" ) : 0;

  COMPLEX8VectorFFT( &stat, bvec, avec, pfwd );
  TestStatus( &stat, CODES( 0 ), 1 );

  for ( i = 0; i < n; ++i )
  {
    fp ? fprintf( fp, "%+f\t%+f\n", bvec->data[i].re, bvec->data[i].im ) : 0;
  }
  fp ? fprintf( fp, "\n" ) : 0;

  COMPLEX8VectorFFT( &stat, cvec, bvec, pinv );
  TestStatus( &stat, CODES( 0 ), 1 );

  for ( i = 0; i < n; ++i )
  {
    cvec->data[i].re /= n;
    cvec->data[i].im /= n;
    fp ? fprintf( fp, "%+.0f\t%+.0f\n", cvec->data[i].re, cvec->data[i].im )
       : 0;
    if ( fabs( avec->data[i].re - cvec->data[i].re ) > epsilon )
    {
      fprintf( stderr, "FAIL: IFFT( FFT( a[] ) ) not equal to a[].\n" );
      return 1;
    }
    if ( fabs( avec->data[i].im - cvec->data[i].im ) > epsilon )
    {
      fprintf( stderr, "FAIL: IFFT( FFT( a[] ) ) not equal to a[].\n" );
      return 1;
    }
  }

  DestroyComplexFFTPlan( &stat, &pinv );
  TestStatus( &stat, CODES( 0 ), 1 );

  DestroyComplexFFTPlan( &stat, &pfwd );
  TestStatus( &stat, CODES( 0 ), 1 );

  CDestroyVector( &stat, &cvec );
  TestStatus( &stat, CODES( 0 ), 1 );

  CDestroyVector( &stat, &bvec );
  TestStatus( &stat, CODES( 0 ), 1 );

  CDestroyVector( &stat, &avec );
  TestStatus( &stat, CODES( 0 ), 1 );

  fp ? fprintf( fp, "\nChecking error codes:\n\n" ) : 0;
  CheckErrorCodes();

  LALCheckMemoryLeaks();
  return 0;
}


static void
CheckErrorCodes()
{
  enum { Size = 19 };
  const UINT4 size = Size;
  Status          stat;
  ComplexFFTPlan *plan;
  ComplexFFTPlan  aplan;
  COMPLEX8Vector  avec;
  COMPLEX8Vector  bvec;
  COMPLEX8        adat[Size];
  COMPLEX8        bdat[Size];

  aplan.size = size;

  EstimateFwdComplexFFTPlan( &stat, NULL, size );
  TestStatus( &stat, CODES( COMPLEXFFT_ENULL ), 1 );
  EstimateInvComplexFFTPlan( &stat, NULL, size );
  TestStatus( &stat, CODES( COMPLEXFFT_ENULL ), 1 );
  MeasureFwdComplexFFTPlan( &stat, NULL, size );
  TestStatus( &stat, CODES( COMPLEXFFT_ENULL ), 1 );
  MeasureInvComplexFFTPlan( &stat, NULL, size );
  TestStatus( &stat, CODES( COMPLEXFFT_ENULL ), 1 );

  plan = &aplan;
  EstimateFwdComplexFFTPlan( &stat, &plan, size );
  TestStatus( &stat, CODES( COMPLEXFFT_ENNUL ), 1 );
  EstimateInvComplexFFTPlan( &stat, &plan, size );
  TestStatus( &stat, CODES( COMPLEXFFT_ENNUL ), 1 );
  MeasureFwdComplexFFTPlan( &stat, &plan, size );
  TestStatus( &stat, CODES( COMPLEXFFT_ENNUL ), 1 );
  MeasureInvComplexFFTPlan( &stat, &plan, size );
  TestStatus( &stat, CODES( COMPLEXFFT_ENNUL ), 1 );

  plan = NULL;
  EstimateFwdComplexFFTPlan( &stat, &plan, 0 );
  TestStatus( &stat, CODES( COMPLEXFFT_ESIZE ), 1 );
  EstimateInvComplexFFTPlan( &stat, &plan, 0 );
  TestStatus( &stat, CODES( COMPLEXFFT_ESIZE ), 1 );
  MeasureFwdComplexFFTPlan( &stat, &plan, 0 );
  TestStatus( &stat, CODES( COMPLEXFFT_ESIZE ), 1 );
  MeasureInvComplexFFTPlan( &stat, &plan, 0 );
  TestStatus( &stat, CODES( COMPLEXFFT_ESIZE ), 1 );

  DestroyComplexFFTPlan( &stat, &plan );
  TestStatus( &stat, CODES( COMPLEXFFT_ENULL ), 1 );
  DestroyComplexFFTPlan( &stat, NULL );
  TestStatus( &stat, CODES( COMPLEXFFT_ENULL ), 1 );

  plan        = &aplan;
  avec.length = size;
  bvec.length = size;
  avec.data   = adat;
  bvec.data   = bdat;
  COMPLEX8VectorFFT( &stat, NULL, &avec, plan );
  TestStatus( &stat, CODES( COMPLEXFFT_ENULL ), 1 );
  COMPLEX8VectorFFT( &stat, &bvec, NULL, plan );
  TestStatus( &stat, CODES( COMPLEXFFT_ENULL ), 1 );
  COMPLEX8VectorFFT( &stat, &bvec, &avec, NULL );
  TestStatus( &stat, CODES( COMPLEXFFT_ENULL ), 1 );

  avec.data = NULL;
  bvec.data = bdat;
  COMPLEX8VectorFFT( &stat, &bvec, &avec, plan );
  TestStatus( &stat, CODES( COMPLEXFFT_ENULL ), 1 );

  avec.data = adat;
  bvec.data = NULL;
  COMPLEX8VectorFFT( &stat, &bvec, &avec, plan );
  TestStatus( &stat, CODES( COMPLEXFFT_ENULL ), 1 );

  avec.data = adat;
  bvec.data = adat;
  COMPLEX8VectorFFT( &stat, &bvec, &avec, plan );
  TestStatus( &stat, CODES( COMPLEXFFT_ESAME ), 1 );

  aplan.size = 0;
  avec.data  = adat;
  avec.data  = bdat;
  COMPLEX8VectorFFT( &stat, &bvec, &avec, plan );
  TestStatus( &stat, CODES( COMPLEXFFT_ESIZE ), 1 );

  aplan.size  = size;
  avec.length = 0;
  bvec.length = size;
  COMPLEX8VectorFFT( &stat, &bvec, &avec, plan );
  TestStatus( &stat, CODES( COMPLEXFFT_ESZMM ), 1 );

  aplan.size  = size;
  avec.length = size;
  bvec.length = 0;
  COMPLEX8VectorFFT( &stat, &bvec, &avec, plan );
  TestStatus( &stat, CODES( COMPLEXFFT_ESZMM ), 1 );

  return;
}


/*
 * TestStatus()
 *
 * Routine to check that the status code status->statusCode agrees with one of
 * the codes specified in the space-delimited string ignored; if not,
 * exit to the system with code exitcode.
 *
 */
static void
TestStatus( Status *status, const char *ignored, int exitcode )
{
  char  str[64];
  char *tok;

  if ( verbose )
  {
    REPORTSTATUS( status );
  }

  if ( strncpy( str, ignored, sizeof( str ) ) )
  {
    if ( (tok = strtok( str, " " ) ) )
    {
      do
      {
        if ( status->statusCode == atoi( tok ) )
        {
          return;
        }
      }
      while ( ( tok = strtok( NULL, " " ) ) );
    }
    else
    {
      if ( status->statusCode == atoi( tok ) )
      {
        return;
      }
    }
  }

  fprintf( stderr, "\nExiting to system with code %d\n", exitcode );
  exit( exitcode );
}


/*
 *
 * ClearStatus()
 *
 * Recursively applies DETATCHSTATUSPTR() to status structure to destroy
 * linked list of statuses.
 *
 */
void
ClearStatus( Status *status )
{
  if ( status->statusPtr )
  {
    ClearStatus( status->statusPtr );
    DETATCHSTATUSPTR( status );
  }
}


/*
 *
 * Usage()
 *
 * Prints a usage message for program program and exits with code exitcode.
 *
 */
static void
Usage( const char *program, int exitcode )
{
  fprintf( stderr, "Usage: %s [options]\n", program );
  fprintf( stderr, "Options:\n" );
  fprintf( stderr, "  -h         print this message\n" );
  fprintf( stderr, "  -q         quiet: run silently\n" );
  fprintf( stderr, "  -v         verbose: print extra information\n" );
  fprintf( stderr, "  -d level   set debuglevel to level\n" );
  exit( exitcode );
}


/*
 * ParseOptions()
 *
 * Parses the argc - 1 option strings in argv[].
 *
 */
static void
ParseOptions( int argc, char *argv[] )
{
  while ( 1 )
  {
    int c = -1;

    c = getopt( argc, argv, "hqvd:" );
    if ( c == -1 )
    {
      break;
    }

    switch ( c )
    {
      case 'd': /* set debug level */
        debuglevel = atoi( optarg );
        break;

      case 'v': /* verbose */
        ++verbose;
        break;

      case 'q': /* quiet: run silently (ignore error messages) */
        freopen( "/dev/null", "w", stderr );
        freopen( "/dev/null", "w", stdout );
        break;

      case 'h':
        Usage( argv[0], 0 );
        break;

      default:
        Usage( argv[0], 1 );
    }

  }

  if ( optind < argc )
  {
    Usage( argv[0], 1 );
  }

  return;
}

