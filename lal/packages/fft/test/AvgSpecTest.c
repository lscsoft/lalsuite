#include <math.h>
#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Units.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Window.h>
#include <lal/Random.h>


#define TESTSTATUS( s ) \
  if ( (s)->statusCode ) { REPORTSTATUS( s ); exit( 1 ); } else ((void)0)

int main( void )
{
  UINT4 seglen = 1024 * 1024;
 /* UINT4 stride = seglen / 2; */
  UINT4 stride = seglen;
  UINT4 numovrlpseg = 8;
  /* UINT4 reclen = numovrlpseg * stride + stride; */
  UINT4 reclen = numovrlpseg * stride;

  static LALStatus status;

  static LALWindowParams        winpar;
  static AverageSpectrumParams  specpar;
  static REAL4FrequencySeries   fseries;
  static REAL4TimeSeries        tseries;
  static RandomParams          *randpar;
  REAL8 avg;
  UINT4 j;
  UINT4 k;
  FILE *fp;

  lalDebugLevel = 3;

  /* allocate memory for time and frequency series */
  tseries.deltaT = 0.1;
  LALCreateVector( &status, &tseries.data, reclen );
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

  winpar.length = seglen;
  winpar.type   = Welch;
  winpar.type   = Rectangular;

  specpar.method  = useMean;
  specpar.overlap = seglen - stride;

  LALCreateForwardRealFFTPlan( &status, &specpar.plan, seglen, 0 );
  TESTSTATUS( &status );
  LALCreateREAL4Window( &status, &specpar.window, &winpar );
  TESTSTATUS( &status );

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
  printf( "%g\n", avg );

  /* use the xlal function */
  XLALREAL4AverageSpectrumWelch( &fseries, &tseries, seglen, stride,
      specpar.window, specpar.plan );
  if ( xlalErrno )
  {
    XLAL_PERROR( "main" );
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
  printf( "%g\n", avg );

  /* median-mean method */
  XLALREAL4AverageSpectrumMedianMean( &fseries, &tseries, seglen, stride,
      specpar.window, specpar.plan );
  if ( xlalErrno )
  {
    XLAL_PERROR( "main" );
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
  printf( "%g\n", avg );

  /* median method */
  XLALREAL4AverageSpectrumMedian( &fseries, &tseries, seglen, stride,
      specpar.window, specpar.plan );
  if ( xlalErrno )
  {
    XLAL_PERROR( "main" );
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
  printf( "%g\n", avg );

  /* cleanup */
  LALDestroyREAL4Window( &status, &specpar.window );
  TESTSTATUS( &status );
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
