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

/**** <lalVerbatim file="RealFFTTestCV">
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \subsection{Program \texttt{RealFFTTest.c}}
 * \label{ss:RealFFTTest.c}
 *
 * Tests the routines in \verb+RealFFT.h+.
 *
 * \subsection*{Usage}
 * \begin{verbatim}
 * RealFFTTest [options]
 * Options:
 *   -h         print this message
 *   -q         quiet: run silently
 *   -v         verbose: print extra information
 *   -d level   set lalDebugLevel to level
 *   -m trials  set number of random trials
 *   -n size    set size of FFTs
 * \end{verbatim}
 *
 * Use the \verb+-n+ option to specify the size of the test transform and
 * the \verb+-m+ option to specify the number of test transforms of that size.
 * (Default is to test transforms of size 1 to 128 in unit steps and then
 * powers of two up to 65536.)
 *
 * \subsubsection*{Description}
 * \subsubsection*{Exit codes}
 * \begin{tabular}{|c|l|}
 * \hline
 * Code & Explanation                   \\
 * \hline
 * \tt 0 & Success, normal exit.         \\
 * \tt 1 & Subroutine failed.            \\
 * \hline
 * \end{tabular}
 *
 * \subsubsection*{Uses}
 * \subsubsection*{Notes}
 *
 * \vfill{\footnotesize\input{RealFFTTestCV}}
 *
 **** </lalLaTeX> */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALConfig.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>
#include <lal/VectorOps.h>
#include <config.h>

#define CODES_(x) #x
#define CODES(x) CODES_(x)

NRCSID( MAIN, "$Id$" );

extern char *optarg;
extern int   optind;

int lalDebugLevel = 0;
int verbose       = 0;
UINT4 m_ = 1; /* number of random trials */
UINT4 n_ = 0; /* size of each transform  */

static void
Usage( const char *program, int exitflag );

static void
ParseOptions( int argc, char *argv[] );

static void
TestStatus( LALStatus *status, const char *expectedCodes, int exitCode );

static void
ClearStatus( LALStatus *status );

void LALForwardRealDFT(
    LALStatus      *status,
    COMPLEX8Vector *output,
    REAL4Vector    *input
    );

int main( int argc, char *argv[] )
{
  static LALStatus status;

  RealFFTPlan    *fwd = NULL;
  RealFFTPlan    *rev = NULL;
  REAL4Vector    *dat = NULL;
  REAL4Vector    *rfft = NULL;
  REAL4Vector    *ans = NULL;
  COMPLEX8Vector *dft = NULL;
  COMPLEX8Vector *fft = NULL;
#if LAL_CUDA_ENABLED
  /* The test itself should pass at 1e-4, but it might fail at
   * some rare cases where accuracy is bad for some numbers. */
  REAL8           eps = 3e-4;
#else
  /* very conservative floating point precision */
  REAL8           eps = 1e-6;
#endif
  REAL8           lbn;
  REAL8           ssq;
  REAL8           var;
  REAL8           tol;

  UINT4 nmax;
  UINT4 m;
  UINT4 n;
  UINT4 i;
  UINT4 j;
  UINT4 k;
  UINT4 s = 0;

  FILE *fp;

  ParseOptions( argc, argv );
  m = m_;
  n = n_;

  fp = verbose ? stdout : NULL ;

  if ( n == 0 )
  {
    nmax = 65536;
  }
  else
  {
    nmax = n--;
  }

  while ( n < nmax )
  {
    if ( n < 128 )
    {
      ++n;
    }
    else
    {
      n *= 2;
    }

    LALSCreateVector( &status, &dat, n );
    TestStatus( &status, CODES( 0 ), 1 );
    LALSCreateVector( &status, &rfft, n );
    TestStatus( &status, CODES( 0 ), 1 );
    LALSCreateVector( &status, &ans, n );
    TestStatus( &status, CODES( 0 ), 1 );
    LALCCreateVector( &status, &dft, n / 2 + 1 );
    TestStatus( &status, CODES( 0 ), 1 );
    LALCCreateVector( &status, &fft, n / 2 + 1 );
    TestStatus( &status, CODES( 0 ), 1 );
    LALCreateForwardRealFFTPlan( &status, &fwd, n, 0 );
    TestStatus( &status, CODES( 0 ), 1 );
    LALCreateReverseRealFFTPlan( &status, &rev, n, 0 );
    TestStatus( &status, CODES( 0 ), 1 );

    /*
     *
     * Do m trials of random data.
     *
     */
    for ( i = 0; i < m; ++i )
    {
      srand( s++ ); /* seed the random number generator */

      /*
       *
       * Create data and compute error tolerance.
       *
       * Reference: Kaneko and Liu,
       * "Accumulation of round-off error in fast fourier tranforms"
       * J. Asssoc. Comp. Mach, Vol 17 (No 4) 637-654, October 1970.
       *
       */
      srand( i ); /* seed the random number generator */
      ssq = 0;
      for ( j = 0; j < n; ++j )
      {
        dat->data[j] = 20.0 * rand() / (REAL4)( RAND_MAX + 1.0 ) - 10.0;
        ssq += dat->data[j] * dat->data[j];
        fp ? fprintf( fp, "%e\n", dat->data[j] ) : 0;
      }
      lbn = log( n ) / log( 2 );
      var = 2.5 * lbn * eps * eps * ssq / n;
      tol = 5 * sqrt( var ); /* up to 5 sigma excursions */
      fp ? fprintf( fp, "\neps = %e \ntol = %e\n", eps, tol ) : 0;

      /*
       *
       * Perform forward FFT and DFT (only if n < 100).
       *
       */
      LALForwardRealFFT( &status, fft, dat, fwd );
      TestStatus( &status, CODES( 0 ), 1 );
      LALREAL4VectorFFT( &status, rfft, dat, fwd );
      TestStatus( &status, CODES( 0 ), 1 );
      LALREAL4VectorFFT( &status, ans, rfft, rev );
      TestStatus( &status, CODES( 0 ), 1 );
      fp ?  fprintf( fp, "rfft()\t\trfft(rfft())\trfft(rfft())\n\n"  ) : 0;
      for ( j = 0; j < n; ++j )
      {
        fp ? fprintf( fp, "%e\t%e\t%e\n",
            rfft->data[j], ans->data[j], ans->data[j] / n ) : 0;
      }
      if ( n < 128 )
      {
        LALForwardRealDFT( &status, dft, dat );
        TestStatus( &status, CODES( 0 ), 1 );

        /*
         *
         * Check accuracy of FFT vs DFT.
         *
         */
        fp ? fprintf( fp, "\nfftre\t\tfftim\t\t" ) : 0;
        fp ? fprintf( fp, "dtfre\t\tdftim\n" ) : 0;
        for ( k = 0; k <= n / 2; ++k )
        {
          REAL8 fftre = fft->data[k].re;
          REAL8 fftim = fft->data[k].im;
          REAL8 dftre = dft->data[k].re;
          REAL8 dftim = dft->data[k].im;
          REAL8 errre = fabs( dftre - fftre );
          REAL8 errim = fabs( dftim - fftim );
          REAL8 avere = fabs( dftre + fftre ) / 2 + eps;
          REAL8 aveim = fabs( dftim + fftim ) / 2 + eps;
          REAL8 ferre = errre / avere;
          REAL8 ferim = errim / aveim;
          fp ? fprintf( fp, "%e\t%e\t", fftre, fftim ) : 0;
          fp ? fprintf( fp, "%e\t%e\n", dftre, dftim ) : 0;
          /* fp ? fprintf( fp, "%e\t%e\t", errre, errim ) : 0; */
          /* fp ? fprintf( fp, "%e\t%e\n", ferre, ferim ) : 0; */
          if ( ferre > eps && errre > tol )
          {
            fputs( "FAIL: Incorrect result from forward transform\n", stderr );
            fprintf( stderr, "\tdifference = %e\n", errre );
            fprintf( stderr, "\ttolerance  = %e\n", tol );
            fprintf( stderr, "\tfrac error = %e\n", ferre );
            fprintf( stderr, "\tprecision  = %e\n", eps );
            return 1;
          }
          if ( ferim > eps && errim > tol )
          {
            fputs( "FAIL: Incorrect result from forward transform\n", stderr );
            fprintf( stderr, "\tdifference = %e\n", errim );
            fprintf( stderr, "\ttolerance  = %e\n", tol );
            fprintf( stderr, "\tfrac error = %e\n", ferim );
            fprintf( stderr, "\tprecision  = %e\n", eps );
            return 1;
          }
        }
      }

      /*
       *
       * Perform reverse FFT and check accuracy vs original data.
       *
       */
      LALReverseRealFFT( &status, ans, fft, rev );
      TestStatus( &status, CODES( 0 ), 1 );
      fp ? fprintf( fp, "\ndat->data[j]\tans->data[j] / n\n" ) : 0;
      for ( j = 0; j < n; ++j )
      {
        REAL8 err = fabs( dat->data[j] - ans->data[j] / n );
        REAL8 ave = fabs( dat->data[j] + ans->data[j] / n ) / 2 + eps;
        REAL8 fer = err / ave;
        fp ? fprintf( fp, "%e\t%e\n", dat->data[j], ans->data[j] / n ) : 0;
        /* fp ? fprintf( fp, "%e\t%e\n", err, fer ) : 0; */
        if ( fer > eps && err > tol )
        {
          fputs( "FAIL: Incorrect result after reverse transform\n", stderr );
          fprintf( stderr, "\tdifference = %e\n", err );
          fprintf( stderr, "\ttolerance  = %e\n", tol );
          fprintf( stderr, "\tfrac error = %e\n", fer );
          fprintf( stderr, "\tprecision  = %e\n", eps );
          return 1;
        }
      }
    }

    LALSDestroyVector( &status, &dat );
    TestStatus( &status, CODES( 0 ), 1 );
    LALSDestroyVector( &status, &rfft );
    TestStatus( &status, CODES( 0 ), 1 );
    LALSDestroyVector( &status, &ans );
    TestStatus( &status, CODES( 0 ), 1 );
    LALCDestroyVector( &status, &dft );
    TestStatus( &status, CODES( 0 ), 1 );
    LALCDestroyVector( &status, &fft );
    TestStatus( &status, CODES( 0 ), 1 );
    LALDestroyRealFFTPlan( &status, &fwd );
    TestStatus( &status, CODES( 0 ), 1 );
    LALDestroyRealFFTPlan( &status, &rev );
    TestStatus( &status, CODES( 0 ), 1 );
  }

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
ClearStatus( LALStatus *status )
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
  fprintf( stderr, "  -d level   set lalDebugLevel to level\n" );
  fprintf( stderr, "  -m trials  set number of random trials\n" );
  fprintf( stderr, "  -n size    set size of FFTs\n" );
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

    c = getopt( argc, argv, "hqvd:m:n:" );
    if ( c == -1 )
    {
      break;
    }

    switch ( c )
    {
      case 'n': /* set FFT size */
        n_ = atoi( optarg );
        break;

      case 'm': /* set number of trials */
        m_ = atoi( optarg );
        break;

      case 'd': /* set debug level */
        lalDebugLevel = atoi( optarg );
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


static void
LALDFT(
    LALStatus      *status,
    COMPLEX8Vector *output,
    COMPLEX8Vector *input,
    INT4            sign
    )
{
  UINT4 n;
  UINT4 j;
  UINT4 k;

  INITSTATUS( status, "DFT", MAIN );

  n = output->length;

  for ( k = 0; k < n; ++k )
  {
    REAL8 sre = 0;
    REAL8 sim = 0;
    for ( j = 0; j < n; ++j )
    {
      REAL8 phi = sign * LAL_TWOPI * (REAL8)(j) * (REAL8)(k) / (REAL8)(n);
      REAL8 c   = cos( phi );
      REAL8 s   = sin( phi );
      REAL8 re  = input->data[j].re;
      REAL8 im  = input->data[j].im;
      sre += c * re - s * im;
      sim += c * im + s * re;
    }
    output->data[k].re = sre;
    output->data[k].im = sim;
  }

  RETURN( status );
}

void LALForwardRealDFT(
    LALStatus      *status,
    COMPLEX8Vector *output,
    REAL4Vector    *input
    )
{
  COMPLEX8Vector *a = NULL;
  COMPLEX8Vector *b = NULL;
  UINT4 n;
  UINT4 j;
  UINT4 k;

  INITSTATUS( status, "LALForwardRealDFT", MAIN );
  ATTATCHSTATUSPTR( status );

  n = input->length;

  TRY( LALCCreateVector( status->statusPtr, &a, n ), status );
  TRY( LALCCreateVector( status->statusPtr, &b, n ), status );

  for ( j = 0; j < n; ++j )
  {
    a->data[j].re = input->data[j];
    a->data[j].im = 0;
  }

  TRY( LALDFT( status->statusPtr, b, a, -1 ), status );

  for ( k = 0; k <= n / 2; ++k )
  {
    output->data[k].re = b->data[k].re;
    output->data[k].im = b->data[k].im;
  }

  TRY( LALCDestroyVector( status->statusPtr, &a ), status );
  TRY( LALCDestroyVector( status->statusPtr, &b ), status );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

