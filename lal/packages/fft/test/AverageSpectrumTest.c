#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/TimeFreqFFT.h>
#include <lal/RealFFT.h>
#include <lal/Window.h>

#define TESTSTATUS( s ) \
  if ( (s)->statusCode ) { REPORTSTATUS( s ); exit( 1 ); } else \
((void)0)

int main( void )
{
  const UINT4 n = 16;
  static LALWindowParams winpar;
  static AverageSpectrumParams specpar;
  static REAL4FrequencySeries fseries;
  static REAL4TimeSeries tseries;
  static LALStatus status;
  UINT4 i;

  /* allocate memory for time and frequency series */
  tseries.deltaT = 1;
  LALCreateVector( &status, &tseries.data, n );
  TESTSTATUS( &status );
  LALCreateVector( &status, &fseries.data, n / 2 + 1 );
  TESTSTATUS( &status );

  /* set time series data to be a unit impulse */
  memset( tseries.data->data, 0, tseries.data->length * sizeof( 
        *tseries.data->data ) );
  tseries.data->data[0] = 1;

  /* prepare window paramers */
  winpar.length = n;
  winpar.type = Rectangular;

  /* prepare average spectrum parameters */
  specpar.method  = useMean;
  specpar.overlap = 0;
  LALCreateForwardRealFFTPlan( &status, &specpar.plan, n, 0 );
  TESTSTATUS( &status );
  LALCreateREAL4Window( &status, &specpar.window, &winpar );
  TESTSTATUS( &status );

  /* compute spectrum */
  LALREAL4AverageSpectrum( &status, &fseries, &tseries, &specpar );
  TESTSTATUS( &status );

  /* output results */
  for ( i = 0; i < fseries.data->length; ++i )
    fprintf( stdout, "%e\t%e\n", i * fseries.deltaF, 
        fseries.data->data[i] );

  /* cleanup */
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
