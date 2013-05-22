/*
*  Copyright (C) 2007 Jolien Creighton
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
   \file
   \ingroup VectorOps_h

   \brief Tests the routines in \ref VectorOps.h.  Exercises some of the error
   conditions and makes sure that they work.

   \heading{Usage}
   \code
   VectorOpsTest [options]
   Options:
   -h         print help
   -q         quiet: run silently
   -v         verbose: print extra information
   -d level   set lalDebugLevel to level
   \endcode

   \heading{Exit codes}
   <table><tr><th>Code</th><th>Explanation</th></tr>
   <tr><td>0</td><td>Success, normal exit.</td></tr>
   <tr><td>1</td><td>Subroutine failed.</td></tr>
   </table>

*/

/** \cond DONT_DOXYGEN */

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <config.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/VectorOps.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#define CODES_(x) #x
#define CODES(x) CODES_(x)

extern char *optarg;
extern int   optind;

int verbose    = 0;

static void
Usage( const char *program, int exitflag );

static void
ParseOptions( int argc, char *argv[] );

static void
TestStatus( LALStatus *status, const char *expectedCodes, int exitCode );

int
main ( int argc, char *argv[] )
{
  const int size = 8;
  COMPLEX8Vector *z1 = NULL;
  COMPLEX8Vector *z2 = NULL;
  COMPLEX8Vector *z3 = NULL;
  REAL4Vector    *x1 = NULL;
  REAL4Vector    *x2 = NULL;
  REAL4Vector    *x3 = NULL;
  REAL4Vector    *y_1 = NULL;
  REAL4Vector    *y2 = NULL;
  REAL4Vector    *y3 = NULL;
  static LALStatus   status;
  INT4            i;


  ParseOptions( argc, argv );

  LALCCreateVector(&status, &z1, size);
  TestStatus( &status, CODES(0), 1 );
  LALCCreateVector(&status, &z2, size);
  TestStatus( &status, CODES(0), 1 );
  LALCCreateVector(&status, &z3, size);
  TestStatus( &status, CODES(0), 1 );
  LALSCreateVector(&status, &x1, size);
  TestStatus( &status, CODES(0), 1 );
  LALSCreateVector(&status, &x2, size);
  TestStatus( &status, CODES(0), 1 );
  LALSCreateVector(&status, &x3, size);
  TestStatus( &status, CODES(0), 1 );
  LALSCreateVector(&status, &y_1, size/2);
  TestStatus( &status, CODES(0), 1 );
  y2         = (REAL4Vector *)LALMalloc(sizeof(REAL4Vector));
  y2->data   = NULL;
  y2->length = size;
  y3         = (REAL4Vector *)LALMalloc(sizeof(REAL4Vector));
  y3->data   = (REAL4 *)LALMalloc(size*sizeof(REAL4));
  y3->length = 0;

  for (i = 0; i < size; ++i)
  {
    z1->data[i] = 1 + i;
    z1->data[i] += I * (2 + i*i);
    z2->data[i] = 3 + i + i*i*i;
    z2->data[i] += I * (4 + i*i + i*i*i);
    x1->data[i]    = 5 + i + i*i;
    x2->data[i]    = 6 + i + i*i + i*i*i;
  }

  if (verbose) printf("\n");
  LALCCVectorMultiply(&status, z3, z1, z2);
  TestStatus( &status, CODES(0), 1 );
  for (i = 0; i < size; ++i)
    if (verbose) printf("(% 6.0f,% 6.0f) x (% 6.0f,% 6.0f) = (% 6.0f,% 6.0f)\n",
        crealf(z1->data[i]), cimagf(z1->data[i]),
        crealf(z2->data[i]), cimagf(z2->data[i]),
        crealf(z3->data[i]), cimagf(z3->data[i]));

  if (verbose) printf("\n");
  LALCCVectorMultiplyConjugate(&status, z3, z1, z2);
  TestStatus( &status, CODES(0), 1 );
  for (i = 0; i < size; ++i)
    if (verbose) printf("(% 6.0f,% 6.0f) x (% 6.0f,% 6.0f)* = (% 6.0f,% 6.0f)\n",
        crealf(z1->data[i]), cimagf(z1->data[i]),
        crealf(z2->data[i]), cimagf(z2->data[i]),
        crealf(z3->data[i]), cimagf(z3->data[i]));

  if (verbose) printf("\n");
  LALCCVectorDivide(&status, z3, z1, z2);
  TestStatus( &status, CODES(0), 1 );
  for (i = 0; i < size; ++i)
    if (verbose) printf("(% 6.0f,% 6.0f) / (% 6.0f,% 6.0f) = (% 9.6f,% 9.6f)\n",
        crealf(z1->data[i]), cimagf(z1->data[i]),
        crealf(z2->data[i]), cimagf(z2->data[i]),
        crealf(z3->data[i]), cimagf(z3->data[i]));

  if (verbose) printf("\n");
  LALSCVectorMultiply(&status, z3, x1, z1);
  TestStatus( &status, CODES(0), 1 );
  for (i = 0; i < size; ++i)
    if (verbose) printf("% 6.0f x (% 6.0f,% 6.0f) = (% 6.0f,% 6.0f)\n",
        crealf(x1->data[i]),
        crealf(z1->data[i]), cimagf(z1->data[i]),
        crealf(z3->data[i]), cimagf(z3->data[i]));

  if (verbose) printf("\n");
  LALSSVectorMultiply(&status, x3, x1, x2);
  TestStatus( &status, CODES(0), 1 );
  for (i = 0; i < size; ++i)
    if (verbose) printf("% 6.0f x % 6.0f = % 6.0f\n",
        x1->data[i], x2->data[i], x3->data[i]);

  if (verbose) printf("\n");
#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
    LALSSVectorMultiply(&status, x3, x1, NULL);
    TestStatus( &status, CODES(VECTOROPSH_ENULL), 1 );
    LALSSVectorMultiply(&status, x3, y2, x2);
    TestStatus( &status, CODES(VECTOROPSH_ENULL), 1 );
    LALSSVectorMultiply(&status, y3, x1, x2);
    TestStatus( &status, CODES(VECTOROPSH_ESIZE), 1 );
    LALSSVectorMultiply(&status, x3, x1, y_1);
    TestStatus( &status, CODES(VECTOROPSH_ESZMM), 1 );
  }
#endif

  LALCDestroyVector(&status, &z1);
  TestStatus( &status, CODES(0), 1 );
  LALCDestroyVector(&status, &z2);
  TestStatus( &status, CODES(0), 1 );
  LALCDestroyVector(&status, &z3);
  TestStatus( &status, CODES(0), 1 );
  LALSDestroyVector(&status, &x1);
  TestStatus( &status, CODES(0), 1 );
  LALSDestroyVector(&status, &x2);
  TestStatus( &status, CODES(0), 1 );
  LALSDestroyVector(&status, &x3);
  TestStatus( &status, CODES(0), 1 );
  LALSDestroyVector(&status, &y_1);
  TestStatus( &status, CODES(0), 1 );
  LALFree(y2);
  LALFree(y3->data);
  LALFree(y3);

  x1 = x2 = x3 = y_1 = y2 = y3 = NULL;
  z1 = z2 = z3 = NULL;


  LALCCreateVector(&status, &z1, size);
  TestStatus( &status, CODES(0), 1 );

  LALSCreateVector(&status, &x1, size);
  TestStatus( &status, CODES(0), 1 );
  LALSCreateVector(&status, &x2, size);
  TestStatus( &status, CODES(0), 1 );
  LALSCreateVector(&status, &x3, size);
  TestStatus( &status, CODES(0), 1 );


  for (i = 0; i < size; ++i)
  {
    z1->data[i] = (12.0 + i) *cos(LAL_PI/3.0*i);
    z1->data[i] += I * (12.0 + i )*sin(LAL_PI/3.0*i);
  }

  if (verbose) printf("\n");
  LALCVectorAbs(&status, x1, z1);
  TestStatus( &status, CODES(0), 1 );
  for (i = 0; i < size; ++i)
    if (verbose) printf(" Abs(% f,%f)  = %f \n",
        crealf(z1->data[i]), cimagf(z1->data[i]),
        x1->data[i]);

  LALCVectorAngle(&status, x2, z1);
  TestStatus( &status, CODES(0), 1 );
  for (i = 0; i < size; ++i)
    if (verbose) printf(" Angle(%f,%f)  = %f \n",
        crealf(z1->data[i]), cimagf(z1->data[i]),
        x2->data[i]);

  LALUnwrapREAL4Angle(&status, x3, x2);
  TestStatus( &status, CODES(0), 1 );
  for (i = 0; i < size; ++i)
    if (verbose) printf(" Unwrap Phase Angle ( %f )  = %f \n",
        x2->data[i],
        x3->data[i]);


  LALSCreateVector(&status, &y_1, size/2);
  TestStatus( &status, CODES(0), 1 );

  y2         = (REAL4Vector *)LALMalloc(sizeof(REAL4Vector));
  y2->data   = NULL;
  y2->length = size;

  y3         = (REAL4Vector *)LALMalloc(sizeof(REAL4Vector));
  y3->data   = (REAL4 *)LALMalloc(size*sizeof(REAL4));
  y3->length = 0;

  if (verbose) printf("\n");

#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
    LALCVectorAbs(&status, x1, NULL);
    TestStatus( &status, CODES(VECTOROPSH_ENULL), 1 );
    LALCVectorAbs(&status, NULL, z1);
    TestStatus( &status, CODES(VECTOROPSH_ENULL), 1 );
    LALCVectorAbs(&status, y_1, z1);
    TestStatus( &status, CODES(VECTOROPSH_ESZMM), 1 );
    LALCVectorAbs(&status, y2, z1);
    TestStatus( &status, CODES(VECTOROPSH_ENULL), 1 );
    LALCVectorAbs(&status, y3, z1);
    TestStatus( &status, CODES(VECTOROPSH_ESIZE), 1 );


    LALCVectorAngle(&status, x2, NULL);
    TestStatus( &status, CODES(VECTOROPSH_ENULL), 1 );
    LALCVectorAngle(&status, NULL, z1);
    TestStatus( &status, CODES(VECTOROPSH_ENULL), 1 );
    LALCVectorAngle(&status, y_1, z1);
    TestStatus( &status, CODES(VECTOROPSH_ESZMM), 1 );
    LALCVectorAngle(&status, y2, z1);
    TestStatus( &status, CODES(VECTOROPSH_ENULL), 1 );
    LALCVectorAngle(&status, y3, z1);
    TestStatus( &status, CODES(VECTOROPSH_ESIZE), 1 );

    LALUnwrapREAL4Angle(&status, x3, NULL);
    TestStatus( &status, CODES(VECTOROPSH_ENULL), 1 );
    LALUnwrapREAL4Angle(&status, NULL, x2);
    TestStatus( &status, CODES(VECTOROPSH_ENULL), 1 );
    LALUnwrapREAL4Angle(&status, y_1, x2);
    TestStatus( &status, CODES(VECTOROPSH_ESZMM), 1 );
    LALUnwrapREAL4Angle(&status, y2, x2);
    TestStatus( &status, CODES(VECTOROPSH_ENULL), 1 );
    LALUnwrapREAL4Angle(&status, y3, x2);
    TestStatus( &status, CODES(VECTOROPSH_ESIZE), 1 );
    LALUnwrapREAL4Angle(&status, x2, x2);
    TestStatus( &status, CODES(VECTOROPSH_ESAME), 1 );
  }
#endif


  LALCDestroyVector(&status, &z1);
  TestStatus( &status, CODES(0), 1 );

  LALSDestroyVector(&status, &x1);
  TestStatus( &status, CODES(0), 1 );
  LALSDestroyVector(&status, &x2);
  TestStatus( &status, CODES(0), 1 );
  LALSDestroyVector(&status, &x3);
  TestStatus( &status, CODES(0), 1 );

  LALSDestroyVector(&status, &y_1);
  TestStatus( &status, CODES(0), 1 );
  LALFree(y2);
  LALFree(y3->data);
  LALFree(y3);

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

  if ( strncpy( str, ignored, sizeof( str ) ) )
  {
    if ( ( tok = strtok( str, " " ) ) )
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
