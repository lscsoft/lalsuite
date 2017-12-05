/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton
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

#include <math.h>
#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Units.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Window.h>
#include <lal/Random.h>
#include <lal/ResampleTimeSeries.h>

#define TESTSTATUS( s ) \
  if ( (s)->statusCode ) { REPORTSTATUS( s ); exit( 1 ); } else ((void)0)

int main( void )
{
  UINT4 resamplefac = 4;
  UINT4 seglen = 64 * 1024; /* after resampling */
  UINT4 numovrlpseg = 10; /* after resampling */
  UINT4 stride = seglen / 2; /* after resampling */
  UINT4 reclen = numovrlpseg * stride + stride; /* after resampling */
  /* UINT4 numovrlpseg = 8; */ /* after resampling */
  /* UINT4 stride = seglen; */ /* after resampling */
  /* UINT4 reclen = numovrlpseg * stride; */ /* after resampling */

  static LALStatus status;

  static AverageSpectrumParams  specpar;
  static REAL4FrequencySeries   fseries;
  static REAL4TimeSeries        tseries;
  static RandomParams          *randpar;
  static ResampleTSParams       resamplepar;
  static PassBandParamStruc     highpasspar;
  REAL8 avg;
  UINT4 j;
  UINT4 k;
  FILE *fp;


  /* allocate memory for time and frequency series */
  tseries.deltaT = 0.1;
  LALCreateVector( &status, &tseries.data, reclen * resamplefac );
  TESTSTATUS( &status );
  LALCreateVector( &status, &fseries.data, seglen / 2 + 1 );
  TESTSTATUS( &status );

  for ( j = 0; j < tseries.data->length; ++j )
    tseries.data->data[j] = sin( 0.1 * j );
  /*
  memset( tseries.data->data, 0,
      tseries.data->length * sizeof( *tseries.data->data ) );
  tseries.data->data[seglen/2-1] = 1;
  tseries.data->data[seglen/2] = 1;
  tseries.data->data[seglen/2+1] = 1;
  */
  /* create white Gaussian noise */
  LALCreateRandomParams( &status, &randpar, 0 );
  TESTSTATUS( &status );
  LALNormalDeviates( &status, tseries.data, randpar );
  TESTSTATUS( &status );
  LALDestroyRandomParams( &status, &randpar );
  TESTSTATUS( &status );

  /* resample */
  resamplepar.deltaT = resamplefac * tseries.deltaT;
  resamplepar.filterType = defaultButterworth;
  LALResampleREAL4TimeSeries( &status, &tseries, &resamplepar );
  TESTSTATUS( &status );

  /* do some simple colouring: high-pass filter */
  highpasspar.nMax = 4;
  highpasspar.f1   = -1;
  highpasspar.a1   = -1;
  highpasspar.f2   = 0.1 / tseries.deltaT; /* ~20% of Nyquist */
  highpasspar.a2   = 0.9; /* this means 10% attenuation at f2 */
  LALDButterworthREAL4TimeSeries( &status, &tseries, &highpasspar );
  TESTSTATUS( &status );

  specpar.method  = useMean;
  specpar.overlap = seglen - stride;

  LALCreateForwardRealFFTPlan( &status, &specpar.plan, seglen, 0 );
  TESTSTATUS( &status );
  specpar.window = XLALCreateRectangularREAL4Window(seglen);

  /* compute spectrum */
  LALREAL4AverageSpectrum( &status, &fseries, &tseries, &specpar );
  TESTSTATUS( &status );

  /* output results */
  fp = fopen( "out1.dat", "w" );
  for ( k = 0; k < fseries.data->length; ++k )
    fprintf( fp, "%e\t%e\n", k * fseries.deltaF, fseries.data->data[k] );
  fclose( fp );

  /* compute average */
  avg = 0;
  for ( k = 1; k < fseries.data->length - 1; ++k )
    avg += fseries.data->data[k];
  avg /= ( fseries.data->length - 2 );
  printf( "lal mean:\t%g\n", avg );

  /* use the xlal function */
  XLALREAL4AverageSpectrumWelch( &fseries, &tseries, seglen, stride,
      specpar.window, specpar.plan );
  if ( xlalErrno )
  {
    XLAL_PERROR();
    exit( 1 );
  }

  /* output results */
  fp = fopen( "out2.dat", "w" );
  for ( k = 0; k < fseries.data->length; ++k )
    fprintf( fp, "%e\t%e\n", k * fseries.deltaF, fseries.data->data[k] );
  fclose( fp );

  /* compute average */
  avg = 0;
  for ( k = 1; k < fseries.data->length - 1; ++k )
    avg += fseries.data->data[k];
  avg /= ( fseries.data->length - 2 );
  printf( "xlal mean:\t%g\n", avg );

  /* median-mean method */
  XLALREAL4AverageSpectrumMedianMean( &fseries, &tseries, seglen, stride,
      specpar.window, specpar.plan );
  if ( xlalErrno )
  {
    XLAL_PERROR();
    exit( 1 );
  }

  /* output results */
  fp = fopen( "out3.dat", "w" );
  for ( k = 0; k < fseries.data->length; ++k )
    fprintf( fp, "%e\t%e\n", k * fseries.deltaF, fseries.data->data[k] );
  fclose( fp );

  /* compute average */
  avg = 0;
  for ( k = 1; k < fseries.data->length - 1; ++k )
    avg += fseries.data->data[k];
  avg /= ( fseries.data->length - 2 );
  printf( "med mean:\t%g\n", avg );

  /* median method */
  XLALREAL4AverageSpectrumMedian( &fseries, &tseries, seglen, stride,
      specpar.window, specpar.plan );
  if ( xlalErrno )
  {
    XLAL_PERROR();
    exit( 1 );
  }

  /* output results */
  fp = fopen( "out4.dat", "w" );
  for ( k = 0; k < fseries.data->length; ++k )
    fprintf( fp, "%e\t%e\n", k * fseries.deltaF, fseries.data->data[k] );
  fclose( fp );

  /* compute average */
  avg = 0;
  for ( k = 1; k < fseries.data->length - 1; ++k )
    avg += fseries.data->data[k];
  avg /= ( fseries.data->length - 2 );
  printf( "median:\t\t%g\n", avg );

  /* cleanup */
  XLALDestroyREAL4Window( specpar.window );
  LALDestroyRealFFTPlan( &status, &specpar.plan );
  TESTSTATUS( &status );
  LALDestroyVector( &status, &fseries.data );
  TESTSTATUS( &status );
  LALDestroyVector( &status, &tseries.data );
  TESTSTATUS( &status );

  /* exit */
  LALCheckMemoryLeaks();
  return 0;
}
