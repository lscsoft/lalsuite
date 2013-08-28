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

/**
 * \author Jolien D. E. Creighton
 * \file
 *
 * ### Program \ref MakeFrames.c ###
 *
 * Make some frames with random Gaussian noise.
 *
 * ### Usage ###
 *
 * \code
 * MakeFrames
 * \endcode
 *
 * ### Description ###
 *
 * This program makes some frames with one ADC channel containing random
 * Gaussian noise.
 *
 */

#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Random.h>
#include <lal/Units.h>
#include <lal/LALFrStream.h>

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 1; } else ((void)0)

#ifndef CHANNEL
#define CHANNEL "H1:LSC-AS_Q"
#endif


int main( void )
{
  static LALStatus      status;
  static INT4TimeSeries series;
  const UINT4   duration = 60;      /* duration of frame file */
  const INT4    seed = 10;          /* random number seed     */
  const REAL4   rate = 16384;       /* sample rate (Hz)       */
  const UINT4   npts = duration * rate; /* number of points       */
  FrOutPar      opar = { "F", "TEST", ADCDataChannel, 6, 0, 0 };
  RandomParams *rpar = NULL;
  UINT4 count = 0;
  UINT4 i;

  /* initialize */

  LALCreateRandomParams( &status, &rpar, seed );
  TESTSTATUS( &status );

  strncpy( series.name, CHANNEL, sizeof( series.name ) );
  series.epoch.gpsSeconds = 600000000;
  series.sampleUnits = lalADCCountUnit;
  series.deltaT = 1 / rate;
  LALI4CreateVector( &status, &series.data, npts );
  TESTSTATUS( &status );

  /* generate first frame file worth of data and write it */

  /*
  LALNormalDeviates( &status, series.data, rpar );
  TESTSTATUS( &status );
  */
  for ( i = 0; i < series.data->length; ++i )
  {
    INT8 ns;
    ns  = (INT8)1000000000 * (INT8)series.epoch.gpsSeconds;
    ns += (INT8)series.epoch.gpsNanoSeconds;
    ns += (INT8)( 1e9 * i * series.deltaT );
    ns %= (INT8)1000000000;
    series.data->data[i] = count++;
  }

  LALFrWriteINT4TimeSeries( &status, &series, &opar );
  TESTSTATUS( &status );

  /* generate second frame file worth of data and write it */

  series.epoch.gpsSeconds += duration;
  /*
  LALNormalDeviates( &status, series.data, rpar );
  TESTSTATUS( &status );
  */
  for ( i = 0; i < series.data->length; ++i )
  {
    INT8 ns;
    ns  = (INT8)1000000000 * (INT8)series.epoch.gpsSeconds;
    ns += (INT8)series.epoch.gpsNanoSeconds;
    ns += (INT8)( 1e9 * i * series.deltaT );
    ns %= (INT8)1000000000;
    series.data->data[i] = count++;
  }

  LALFrWriteINT4TimeSeries( &status, &series, &opar );
  TESTSTATUS( &status );

  /* generate third frame file worth of data and write it */

  series.epoch.gpsSeconds += duration;
  /*
  LALNormalDeviates( &status, series.data, rpar );
  TESTSTATUS( &status );
  */
  for ( i = 0; i < series.data->length; ++i )
  {
    INT8 ns;
    ns  = (INT8)1000000000 * (INT8)series.epoch.gpsSeconds;
    ns += (INT8)series.epoch.gpsNanoSeconds;
    ns += (INT8)( 1e9 * i * series.deltaT );
    ns %= (INT8)1000000000;
    series.data->data[i] = count++;
  }

  LALFrWriteINT4TimeSeries( &status, &series, &opar );
  TESTSTATUS( &status );

  //if ( XLALFrWriteINT4TimeSeries( &series, 7 ) < 0 )
    //return 1;

  /* cleanup */

  LALI4DestroyVector( &status, &series.data );
  TESTSTATUS( &status );

  LALDestroyRandomParams( &status, &rpar );
  TESTSTATUS( &status );

  LALCheckMemoryLeaks();
  return 0;
}
