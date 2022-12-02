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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

/**
 * \file
 * \ingroup ComplexFFT_h
 *
 * \brief Tests the routines in \ref ComplexFFT.h.
 *
 * ### Usage ###
 *
 * \code
 * ComplexFFTTest [options]
 * Options:
 * -h         print this message
 * -q         quiet: run silently
 * -v         verbose: print extra information
 * -d level   set lalDebugLevel to level
 * \endcode
 *
 * ### Exit codes ###
 *
 * <table><tr><th>Code</th><th>Explanation</th></tr>
 * <tr><td>0</td><td>Success, normal exit.</td></tr>
 * <tr><td>1</td><td>Subroutine failed.</td></tr>
 * </table>
 *
 * ### Uses ###
 *
 *
 * ### Notes ###
 *
 */

/** \cond DONT_DOXYGEN */
#include <config.h>

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/LALgetopt.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
#include <lal/LALString.h>
#include <config.h>

#define CODES_(x) #x
#define CODES(x) CODES_(x)

int verbose       = 0;

static void
Usage( const char *program, int exitflag );

static void
ParseOptions( int argc, char *argv[] );

static void
TestStatus( LALStatus *status, const char *expectedCodes, int exitCode );

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


  ParseOptions( argc, argv );

  fp = verbose ? stdout : NULL ;

  LALCCreateVector( &status, &avec, n );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCCreateVector( &status, &bvec, n );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCCreateVector( &status, &cvec, n );
  TestStatus( &status, CODES( 0 ), 1 );

  pfwd = XLALCreateForwardCOMPLEX8FFTPlan( n, 0 );

  prev = XLALCreateReverseCOMPLEX8FFTPlan( n, 0 );


  for ( i = 0; i < n; ++i )
  {
    avec->data[i] = rand() % 5 - 2;
    avec->data[i] += I * (rand() % 3 - 1);
    fp ? fprintf( fp, "%+.0f\t%+.0f\n", crealf(avec->data[i]), cimagf(avec->data[i]) )
       : 0;
  }
  fp ? fprintf( fp, "\n" ) : 0;
  fflush( stdout );

  XLALCOMPLEX8VectorFFT( bvec, avec, pfwd );

  for ( i = 0; i < n; ++i )
  {
    fp ? fprintf( fp, "%+f\t%+f\n", crealf(bvec->data[i]), cimagf(bvec->data[i]) ) : 0;
  }
  fp ? fprintf( fp, "\n" ) : 0;
  fflush( stdout );

  XLALCOMPLEX8VectorFFT( cvec, bvec, prev );

  for ( i = 0; i < n; ++i )
  {
    cvec->data[i] /= n;
    fp ? fprintf( fp, "%+.0f\t%+.0f\n", crealf(cvec->data[i]), cimagf(cvec->data[i]) )
       : 0;
    fflush( stdout );
    if ( fabs( creal(avec->data[i] - cvec->data[i]) ) > eps )
    {
      fprintf( stderr, "FAIL: IFFT( FFT( a[] ) ) not equal to a[].\n" );
      fprintf( stderr, "Re(avec->data[%d]) = %e\n", i, crealf(avec->data[i]) );
      fprintf( stderr, "Re(cvec->data[%d]) = %e\n", i, crealf(cvec->data[i]) );
      return 1;
    }
    if ( fabs( cimag(avec->data[i] - cvec->data[i]) ) > eps )
    {
      fprintf( stderr, "FAIL: IFFT( FFT( a[] ) ) not equal to a[].\n" );
      fprintf( stderr, "Im(avec->data[%d]) = %e\n", i, cimagf(avec->data[i]) );
      fprintf( stderr, "Im(cvec->data[%d]) = %e\n", i, cimagf(cvec->data[i]) );
      return 1;
    }
  }

  XLALDestroyCOMPLEX8FFTPlan( prev );
  XLALDestroyCOMPLEX8FFTPlan( pfwd );

  /* Null pointers should be a no-op */
  prev = NULL;
  pfwd = NULL;
  XLALDestroyCOMPLEX8FFTPlan( prev );
  XLALDestroyCOMPLEX8FFTPlan( pfwd );

  LALCDestroyVector( &status, &cvec );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCDestroyVector( &status, &bvec );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCDestroyVector( &status, &avec );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCheckMemoryLeaks();
  return 0;
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

  if ( XLALStringCopy( str, ignored, sizeof( str ) ) )
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
      if ( status->statusCode == atoi( str ) )
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

    c = LALgetopt( argc, argv, "hqvd:" );
    if ( c == -1 )
    {
      break;
    }

    switch ( c )
    {
      case 'd': /* set debug level */
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

  if ( LALoptind < argc )
  {
    Usage( argv[0], 1 );
  }

  return;
}

/** \endcond */
