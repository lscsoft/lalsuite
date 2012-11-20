/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton
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

/**
 * \file
 * \ingroup ComplexFFT_h
 *
 * \brief Tests the routines in \ref ComplexFFT.h.
 *
 * \heading{Usage}
 * \code
 * ComplexFFTTest [options]
 * Options:
 *   -h         print this message
 *   -q         quiet: run silently
 *   -v         verbose: print extra information
 *   -d level   set lalDebugLevel to level
 * \endcode
 *
 * \heading{Exit codes}
 * <table><tr><th>Code</th><th>Explanation</th></tr>
 * <tr><td>0</td><td>Success, normal exit.</td></tr>
 * <tr><td>1</td><td>Subroutine failed.</td></tr>
 * </table>
 *
 * \heading{Uses}
 * \heading{Notes}
 *
 *
 *
*/

/** \cond DONT_DOXYGEN */
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
#include <config.h>

#define CODES_(x) #x
#define CODES(x) CODES_(x)

extern char *optarg;
extern int   optind;

extern int lalDebugLevel;
int verbose       = 0;

static void
Usage( const char *program, int exitflag );

static void
ParseOptions( int argc, char *argv[] );

static void
TestStatus( LALStatus *status, const char *expectedCodes, int exitCode );

static void
CheckErrorCodes( void );

int
main( int argc, char *argv[] )
{
  const UINT4 n   = 17;
#if LAL_CUDA_ENABLED
  const REAL4 eps = 1e-4;
#else
  const REAL4 eps = 1e-6;
#endif

  static LALStatus  status;
  ComplexFFTPlan   *pfwd = NULL;
  ComplexFFTPlan   *prev = NULL;
  COMPLEX8Vector   *avec = NULL;
  COMPLEX8Vector   *bvec = NULL;
  COMPLEX8Vector   *cvec = NULL;

  FILE *fp;
  UINT4 i;

  lalDebugLevel = 0;

  ParseOptions( argc, argv );

  fp = verbose ? stdout : NULL ;

  LALCCreateVector( &status, &avec, n );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCCreateVector( &status, &bvec, n );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCCreateVector( &status, &cvec, n );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCreateForwardComplexFFTPlan( &status, &pfwd, n, 0 );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCreateReverseComplexFFTPlan( &status, &prev, n, 0 );
  TestStatus( &status, CODES( 0 ), 1 );


  for ( i = 0; i < n; ++i )
  {
    avec->data[i].im = rand() % 5 - 2;
    avec->data[i].re = rand() % 3 - 1;
    fp ? fprintf( fp, "%+.0f\t%+.0f\n", avec->data[i].re, avec->data[i].im )
       : 0;
  }
  fp ? fprintf( fp, "\n" ) : 0;
  fflush( stdout );

  LALCOMPLEX8VectorFFT( &status, bvec, avec, pfwd );
  TestStatus( &status, CODES( 0 ), 1 );

  for ( i = 0; i < n; ++i )
  {
    fp ? fprintf( fp, "%+f\t%+f\n", bvec->data[i].re, bvec->data[i].im ) : 0;
  }
  fp ? fprintf( fp, "\n" ) : 0;
  fflush( stdout );

  LALCOMPLEX8VectorFFT( &status, cvec, bvec, prev );
  TestStatus( &status, CODES( 0 ), 1 );

  for ( i = 0; i < n; ++i )
  {
    cvec->data[i].re /= n;
    cvec->data[i].im /= n;
    fp ? fprintf( fp, "%+.0f\t%+.0f\n", cvec->data[i].re, cvec->data[i].im )
       : 0;
    fflush( stdout );
    if ( fabs( avec->data[i].re - cvec->data[i].re ) > eps )
    {
      fprintf( stderr, "FAIL: IFFT( FFT( a[] ) ) not equal to a[].\n" );
      fprintf( stderr, "avec->data[%d].re = %e\n", i, avec->data[i].re );
      fprintf( stderr, "cvec->data[%d].re = %e\n", i, cvec->data[i].re );
      return 1;
    }
    if ( fabs( avec->data[i].im - cvec->data[i].im ) > eps )
    {
      fprintf( stderr, "FAIL: IFFT( FFT( a[] ) ) not equal to a[].\n" );
      fprintf( stderr, "avec->data[%d].im = %e\n", i, avec->data[i].im );
      fprintf( stderr, "cvec->data[%d].im = %e\n", i, cvec->data[i].im );
      return 1;
    }
  }

  LALDestroyComplexFFTPlan( &status, &prev );
  TestStatus( &status, CODES( 0 ), 1 );

  LALDestroyComplexFFTPlan( &status, &pfwd );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCDestroyVector( &status, &cvec );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCDestroyVector( &status, &bvec );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCDestroyVector( &status, &avec );
  TestStatus( &status, CODES( 0 ), 1 );

  fp ? fprintf( fp, "\nChecking error codes:\n\n" ) : 0;
  CheckErrorCodes();

  LALCheckMemoryLeaks();
  return 0;
}

/* need to have definition of tagComplexFFTPlan to corrupt it */
struct
tagComplexFFTPlan
{
  INT4   sign;
  UINT4  size;
  void  *plan;
};

static void
CheckErrorCodes( void )
{
#ifndef LAL_NDEBUG
  enum { Size = 19 };
  const UINT4 size = Size;
  static LALStatus status;
  ComplexFFTPlan *plan;
  ComplexFFTPlan  aplan;
  COMPLEX8Vector  avec;
  COMPLEX8Vector  bvec;
  COMPLEX8        adat[Size];
  COMPLEX8        bdat[Size];

  if ( ! lalNoDebug )
  {
    aplan.size = size;

    LALCreateForwardComplexFFTPlan( &status, NULL, size, 0 );
    TestStatus( &status, CODES( COMPLEXFFTH_ENULL ), 1 );
    LALCreateReverseComplexFFTPlan( &status, NULL, size, 0 );
    TestStatus( &status, CODES( COMPLEXFFTH_ENULL ), 1 );

    plan = &aplan;
    LALCreateForwardComplexFFTPlan( &status, &plan, size, 0 );
    TestStatus( &status, CODES( COMPLEXFFTH_ENNUL ), 1 );
    LALCreateReverseComplexFFTPlan( &status, &plan, size, 0 );
    TestStatus( &status, CODES( COMPLEXFFTH_ENNUL ), 1 );

    plan = NULL;
    LALCreateForwardComplexFFTPlan( &status, &plan, 0, 0 );
    TestStatus( &status, CODES( COMPLEXFFTH_ESIZE ), 1 );
    LALCreateReverseComplexFFTPlan( &status, &plan, 0, 0 );
    TestStatus( &status, CODES( COMPLEXFFTH_ESIZE ), 1 );

    LALDestroyComplexFFTPlan( &status, &plan );
    TestStatus( &status, CODES( COMPLEXFFTH_ENULL ), 1 );
    LALDestroyComplexFFTPlan( &status, NULL );
    TestStatus( &status, CODES( COMPLEXFFTH_ENULL ), 1 );

    plan        = &aplan;
    avec.length = size;
    bvec.length = size;
    avec.data   = adat;
    bvec.data   = bdat;
    LALCOMPLEX8VectorFFT( &status, NULL, &avec, plan );
    TestStatus( &status, CODES( COMPLEXFFTH_ENULL ), 1 );
    LALCOMPLEX8VectorFFT( &status, &bvec, NULL, plan );
    TestStatus( &status, CODES( COMPLEXFFTH_ENULL ), 1 );
    LALCOMPLEX8VectorFFT( &status, &bvec, &avec, NULL );
    TestStatus( &status, CODES( COMPLEXFFTH_ENULL ), 1 );

    avec.data = NULL;
    bvec.data = bdat;
    LALCOMPLEX8VectorFFT( &status, &bvec, &avec, plan );
    TestStatus( &status, CODES( COMPLEXFFTH_ENULL ), 1 );

    avec.data = adat;
    bvec.data = NULL;
    LALCOMPLEX8VectorFFT( &status, &bvec, &avec, plan );
    TestStatus( &status, CODES( COMPLEXFFTH_ENULL ), 1 );

    avec.data = adat;
    bvec.data = adat;
    LALCOMPLEX8VectorFFT( &status, &bvec, &avec, plan );
    TestStatus( &status, CODES( COMPLEXFFTH_ESAME ), 1 );

    aplan.size = 0;
    avec.data  = adat;
    avec.data  = bdat;
    LALCOMPLEX8VectorFFT( &status, &bvec, &avec, plan );
    TestStatus( &status, CODES( COMPLEXFFTH_ESIZE ), 1 );

    aplan.size  = size;
    avec.length = 0;
    bvec.length = size;
    LALCOMPLEX8VectorFFT( &status, &bvec, &avec, plan );
    TestStatus( &status, CODES( COMPLEXFFTH_ESZMM ), 1 );

    aplan.size  = size;
    avec.length = size;
    bvec.length = 0;
    LALCOMPLEX8VectorFFT( &status, &bvec, &avec, plan );
    TestStatus( &status, CODES( COMPLEXFFTH_ESZMM ), 1 );
  }
#endif

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
TestStatus( LALStatus *status, const char *ignored, int exitcode )
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
  fprintf( stderr, "  -d level   set lalDebugLevel to level\n" );
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
  FILE *fp;

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
        lalDebugLevel = atoi( optarg );
        break;

      case 'v': /* verbose */
        ++verbose;
        break;

      case 'q': /* quiet: run silently (ignore error messages) */
        fp = freopen( "/dev/null", "w", stderr );
        if (fp == NULL)
        {
          fprintf(stderr, "Error: Unable to open /dev/null\n");
          exit(1);
        }
        fp = freopen( "/dev/null", "w", stdout );
        if (fp == NULL)
        {
          fprintf(stderr, "Error: Unable to open /dev/null\n");
          exit(1);
        }
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

/** \endcond */
