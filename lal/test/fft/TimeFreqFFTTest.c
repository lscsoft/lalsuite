/*
*  Copyright (C) 2007 Bernd Machenschalk, Duncan Brown, Jolien Creighton
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
 * \ingroup TimeFreqFFT_h
 *
 * \brief Tests the routines in \ref TimeFreqFFT.h.
 *
 * ### Usage ###
 *
 * \code
 * TimeFreqFFTTest [options]
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
 * <tr><td>2</td><td>PSD estimation tolerance exceeded</td></tr>
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
#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/LALgetopt.h>
#include <lal/LALConstants.h>
#include <lal/LALStdio.h>
#include <lal/AVFactories.h>
#include <lal/Units.h>
#include <lal/Random.h>
#include <lal/PrintFTSeries.h>
#include <lal/TimeFreqFFT.h>

#define CODES_(x) #x
#define CODES(x) CODES_(x)

int verbose       = 0;

static void
Usage( const char *program, int exitflag );

static void
ParseOptions( int argc, char *argv[] );

static void
TestStatus( LALStatus *status, const char *expectedCodes, int exitCode );

int main( int argc, char *argv[] )
{
  const UINT4 n  = 65536;
  const REAL4 dt = 1.0 / 16384.0;
  static LALStatus status;

  static REAL4TimeSeries         x;
  static COMPLEX8FrequencySeries X;

  static REAL4TimeSeries         y;
  static REAL4FrequencySeries    Y;

  static COMPLEX8TimeSeries      z;
  static COMPLEX8FrequencySeries Z;

  RealFFTPlan    *fwdRealPlan    = NULL;
  RealFFTPlan    *revRealPlan    = NULL;
  ComplexFFTPlan *fwdComplexPlan = NULL;
  ComplexFFTPlan *revComplexPlan = NULL;
  RandomParams   *randpar        = NULL;

  UINT4 srate[] = { 4096, 9000 };
  UINT4 npts[] = { 262144, 1048576 };
  REAL4 var[] = { 5, 16 };

  UINT4 j, sr, np, vr;


  /*CHAR fname[2048];*/

  ParseOptions( argc, argv );

  LALSCreateVector( &status, &x.data, n );
  TestStatus( &status, CODES( 0 ), 1 );
  LALCCreateVector( &status, &X.data, n / 2 + 1 );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCCreateVector( &status, &z.data, n );
  TestStatus( &status, CODES( 0 ), 1 );
  LALCCreateVector( &status, &Z.data, n );
  TestStatus( &status, CODES( 0 ), 1 );

  fwdRealPlan = XLALCreateForwardREAL4FFTPlan( n, 0 );
  revRealPlan = XLALCreateReverseREAL4FFTPlan( n, 0 );
  fwdComplexPlan = XLALCreateForwardCOMPLEX8FFTPlan( n, 0 );
  revComplexPlan = XLALCreateReverseCOMPLEX8FFTPlan( n, 0 );
  TestStatus( &status, CODES( 0 ), 1 );

  randpar = XLALCreateRandomParams( 100 );


  /*
   *
   * Try the real transform.
   *
   */


  x.f0 = 0;
  x.deltaT = dt;
  x.sampleUnits = lalMeterUnit;
  snprintf( x.name, sizeof( x.name ), "x" );
  XLALNormalDeviates( x.data, randpar );
  for ( j = 0; j < n; ++j ) /* add a 60 Hz line */
  {
    REAL4 t = j * dt;
    x.data->data[j] += 0.1 * cos( LAL_TWOPI * 60.0 * t );
  }
  LALSPrintTimeSeries( &x, "x.out" );

  snprintf( X.name, sizeof( X.name ), "X" );
  XLALREAL4TimeFreqFFT( &X, &x, fwdRealPlan );
  LALCPrintFrequencySeries( &X, "X.out" );

  XLALREAL4FreqTimeFFT( &x, &X, revRealPlan );
  LALSPrintTimeSeries( &x, "xx.out" );


  /*
   *
   * Try the average power spectum.
   *
   */


  for ( np = 0; np < XLAL_NUM_ELEM(npts) ; ++np )
  {
    REAL4Window *window = XLALCreateHannREAL4Window(npts[np]);
    RealFFTPlan *plan = XLALCreateForwardREAL4FFTPlan(npts[np], 0);
    /* length of time series for 7 segments, overlapped by 1/2 segment */
    UINT4 tsLength = npts[np] * 7 - 6 * npts[np] / 2;
    LALCreateVector( &status, &y.data, tsLength );
    TestStatus( &status, CODES( 0 ), 1 );
    LALCreateVector( &status, &Y.data, npts[np] / 2 + 1  );
    TestStatus( &status, CODES( 0 ), 1 );

    for ( sr = 0; sr < XLAL_NUM_ELEM(srate) ; ++sr )
    {
      /* set the sample rate of the time series */
      y.deltaT = 1.0 / (REAL8) srate[sr];
      for ( vr = 0; vr < XLAL_NUM_ELEM(var) ; ++vr )
      {
        REAL4 eps = 1e-6; /* very conservative fp precision */
        REAL4 Sfk = 2.0 * var[vr] * var[vr] * y.deltaT;
        REAL4 sfk;
        REAL4 lbn;
        REAL4 sig;
        REAL4 ssq;
        REAL4 tol;

        /* create the data */
        XLALNormalDeviates( y.data, randpar );
        ssq = 0;
        for ( j = 0; j < y.data->length; ++j )
        {
          y.data->data[j] *= var[vr];
          ssq += y.data->data[j] * y.data->data[j];
        }

        /* compute tolerance for comparison */
        lbn = log( y.data->length ) / log( 2 );
        sig = sqrt( 2.5 * lbn * eps * eps * ssq / y.data->length );
        tol = 5 * sig;

        /* compute the psd and find the average */
        XLALREAL4AverageSpectrumWelch(&Y, &y, npts[np], npts[np] / 2, window, plan);
        TestStatus( &status, CODES( 0 ), 1 );
        {
          unsigned k;
          for(sfk = 0., k = 0; k < Y.data->length; sfk += Y.data->data[k++])
            ;
          sfk /= Y.data->length;
        }

        /* check the result */
        if ( fabs(Sfk-sfk) > tol )
        {
          fprintf( stderr, "FAIL: PSD estimate appears incorrect\n");
          fprintf( stderr, "expected %e, got %e ", Sfk, sfk );
          fprintf( stderr, "(difference = %e, tolerance = %e)\n",
              fabs(Sfk-sfk), tol );
          exit(2);
        }

      }
    }

    /* destroy structures that need to be resized */
    XLALDestroyREAL4FFTPlan( plan );
    XLALDestroyREAL4Window( window );
    LALDestroyVector( &status, &y.data );
    TestStatus( &status, CODES( 0 ), 1 );
    LALDestroyVector( &status, &Y.data );
    TestStatus( &status, CODES( 0 ), 1 );
  }


  /*
   *
   * Try the complex transform.
   *
   */


  z.f0 = 0;
  z.deltaT = dt;
  z.sampleUnits = lalVoltUnit;
  snprintf( z.name, sizeof( z.name ), "z" );
  { /* dirty hack */
    REAL4Vector tmp;
    tmp.length = 2 * z.data->length;
    tmp.data   = (REAL4 *)z.data->data;
    XLALNormalDeviates( &tmp, randpar );
  }
  for ( j = 0; j < n; ++j ) /* add a 50 Hz line and a 500 Hz ringdown */
  {
    REAL4 t = j * dt;
    z.data->data[j] += 0.2 * cos( LAL_TWOPI * 50.0 * t );
    z.data->data[j] += I * exp( -t ) * sin( LAL_TWOPI * 500.0 * t );
  }
  LALCPrintTimeSeries( &z, "z.out" );
  TestStatus( &status, CODES( 0 ), 1 );

  snprintf( Z.name, sizeof( Z.name ), "Z" );
  XLALCOMPLEX8TimeFreqFFT( &Z, &z, fwdComplexPlan );
  LALCPrintFrequencySeries( &Z, "Z.out" );

  XLALCOMPLEX8FreqTimeFFT( &z, &Z, revComplexPlan );
  TestStatus( &status, CODES( 0 ), 1 );
  LALCPrintTimeSeries( &z, "zz.out" );

  XLALDestroyRandomParams( randpar );

  XLALDestroyREAL4FFTPlan( fwdRealPlan );
  XLALDestroyREAL4FFTPlan( revRealPlan );
  XLALDestroyCOMPLEX8FFTPlan( fwdComplexPlan );
  XLALDestroyCOMPLEX8FFTPlan( revComplexPlan );

  LALCDestroyVector( &status, &Z.data );
  TestStatus( &status, CODES( 0 ), 1 );
  LALCDestroyVector( &status, &z.data );
  TestStatus( &status, CODES( 0 ), 1 );

  LALCDestroyVector( &status, &X.data );
  TestStatus( &status, CODES( 0 ), 1 );
  LALSDestroyVector( &status, &x.data );
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
