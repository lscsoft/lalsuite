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
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/TimeFreqFFT.h>
#include <lal/RealFFT.h>
#include <lal/Window.h>
#include <lal/Random.h>

#include <lal/LALRCSID.h>
NRCSID (AVERAGESPECTRUMTESTC,"$Id$");

#define TESTSTATUS( s ) \
  if ( (s)->statusCode ) { REPORTSTATUS( s ); exit( 1 ); } else \
((void)0)

int main( void )
{
  const UINT4 n = 65536;
  const UINT4 m = 8;
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

  /* prepare average spectrum parameters */
  specpar.method  = useMedian;
  specpar.overlap = n / 2;
  /* specpar.overlap = 0; */
  LALCreateForwardRealFFTPlan( &status, &specpar.plan, n, 0 );
  TESTSTATUS( &status );
  specpar.window = XLALCreateWelchREAL4Window(n);

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
