#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/TimeFreqFFT.h>
#include <lal/RealFFT.h>
#include <lal/Window.h>
#include <lal/Random.h>

#define TESTSTATUS( s ) \
  if ( (s)->statusCode ) { REPORTSTATUS( s ); exit( 1 ); } else \
((void)0)

int main( void )
{
  const UINT4 n = 65536;
  const UINT4 m = 8;
  static LALWindowParams winpar;
  static AverageSpectrumParams specpar;
  static REAL4FrequencySeries fseries;
  static REAL4TimeSeries tseries;
  static LALStatus status;
  static RandomParams *randpar;
  REAL8 ave;
  UINT4 i;

  lalDebugLevel = 3;

  /* allocate memory for time and frequency series */
  tseries.deltaT = 1;
  LALCreateVector( &status, &tseries.data, n * m );
  TESTSTATUS( &status );
  LALCreateVector( &status, &fseries.data, n / 2 + 1 );
  TESTSTATUS( &status );

  /* set time series data to be a unit impulse */
  /*
  memset( tseries.data->data, 0, tseries.data->length * sizeof( 
        *tseries.data->data ) );
  tseries.data->data[0] = 1;
  */
  LALCreateRandomParams( &status, &randpar, 1 );
  TESTSTATUS( &status );
  LALNormalDeviates( &status, tseries.data, randpar );
  TESTSTATUS( &status );
  LALDestroyRandomParams( &status, &randpar );
  TESTSTATUS( &status );

  /* prepare window paramers */
  winpar.length = n;
  winpar.type = Welch;
  /* winpar.type = Rectangular; */

  /* prepare average spectrum parameters */
  specpar.method  = useMedian;
  specpar.overlap = n / 2;
  /* specpar.overlap = 0; */
  LALCreateForwardRealFFTPlan( &status, &specpar.plan, n, 0 );
  TESTSTATUS( &status );
  LALCreateREAL4Window( &status, &specpar.window, &winpar );
  TESTSTATUS( &status );

  /* compute spectrum */
  LALREAL4AverageSpectrum( &status, &fseries, &tseries, &specpar );
  TESTSTATUS( &status );

  /* output results -- omit DC & Nyquist */
  /*
  for ( i = 1; i < fseries.data->length - 1; ++i )
    fprintf( stdout, "%e\t%e\n", i * fseries.deltaF, 
        fseries.data->data[i] );
   */

  /* average values of power spectrum (omit DC & Nyquist ) */
  ave = 0;
  for ( i = 1; i < fseries.data->length - 1; ++i )
    ave += fseries.data->data[i];
  ave /= fseries.data->length - 2;
  fprintf( stdout, "median:\t%e\terror:\t%f%%\n", ave, fabs( ave - 2.0 ) / 0.02 );

  /* now do the same for mean */
  specpar.method  = useMean;
  LALREAL4AverageSpectrum( &status, &fseries, &tseries, &specpar );
  TESTSTATUS( &status );
  /* average values of power spectrum (omit DC & Nyquist ) */
  ave = 0;
  for ( i = 1; i < fseries.data->length - 1; ++i )
    ave += fseries.data->data[i];
  ave /= fseries.data->length - 2;
  fprintf( stdout, "mean:\t%e\terror:\t%f%%\n", ave, fabs( ave - 2.0 ) / 0.02 );


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
