/*
 * Copyright (C) 2007 Gregory Mendell
 * Copyright (C) 2010, 2011, 2016 Bernd Machenschalk
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA
 */

#include "config.h"

#include <math.h>

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>

/* GLOBAL VARIABLES */
REAL8 winFncRMS = 1.0; /* 10/05/12 gam; global variable with the RMS of the window function; default value is 1.0 */
REAL8TimeSeries dataDouble;
REAL4TimeSeries dataSingle;

/* FUNCTION PROTOTYPES */
int WindowData( REAL8 r );
int WindowDataTukey2( void );
int WindowDataHann( void );

int WindowData( const REAL8 r )
{
  INT4 k, N, kl, kh;
  /* 10/05/12 gam */
  REAL8 win;
  /* This implementation of a Tukey window is describes
     in the Matlab tukey window documentation */

  /* initialize the RMS of the window function */
  winFncRMS = 0.0;

  N = dataDouble.data->length;
  kl = r / 2 * ( N - 1 ) + 1;
  kh = N - r / 2 * ( N - 1 ) + 1;
  for ( k = 1; k < kl; k++ ) {
    /* dataDouble.data->data[k-1]=dataDouble.data->data[k-1]*0.5*( 1.0 + cos(LAL_TWOPI/r*(k-1)/(N-1) - LAL_PI) ); */
    win = 0.5 * ( 1.0 + cos( LAL_TWOPI / r * ( k - 1 ) / ( N - 1 ) - LAL_PI ) );
    dataDouble.data->data[k - 1] *= win;
    winFncRMS += win * win;
  }
  for ( k = kh; k <= N; k++ ) {
    /* dataDouble.data->data[k-1]=dataDouble.data->data[k-1]*0.5*( 1.0 + cos(LAL_TWOPI/r - LAL_TWOPI/r*(k-1)/(N-1) - LAL_PI) ); */
    win = 0.5 * ( 1.0 + cos( LAL_TWOPI / r - LAL_TWOPI / r * ( k - 1 ) / ( N - 1 ) - LAL_PI ) );
    dataDouble.data->data[k - 1] *= win;
    winFncRMS += win * win;
  }

  /* Add to sum of squares of the window function the parts of window which are equal to 1, and then find RMS value*/
  winFncRMS += ( REAL8 )( kh - kl );
  winFncRMS = sqrt( ( winFncRMS / ( ( REAL8 ) N ) ) );

  return 0;
}

/* Same as window function given in lalapps/src/pulsar/make_sfts.c */
int WindowDataTukey2( void )
{
  /* Define the parameters to make the window */
  INT4 WINSTART = 4096;
  INT4 WINEND = 8192;
  INT4 WINLEN = ( WINEND - WINSTART );
  /* INT4 i; */
  INT4 i, N;
  REAL8 win;


  /* initialize the RMS of the window function */
  winFncRMS = 0.0;

  N = dataDouble.data->length;
  /* window data.  Off portion */
  for ( i = 0; i < WINSTART; i++ ) {
    dataDouble.data->data[i] = 0.0;
    dataDouble.data->data[dataDouble.data->length - 1 - i] = 0.0;
  }
  /* window data, smooth turn-on portion */
  for ( i = WINSTART; i < WINEND; i++ ) {
    win = ( ( sin( ( i - WINSTART ) * LAL_PI / ( WINLEN ) - LAL_PI_2 ) + 1.0 ) / 2.0 );
    dataDouble.data->data[i] *= win;
    dataDouble.data->data[dataDouble.data->length - 1 - i]  *= win;
    winFncRMS += 2.0 * win * win;
  }

  /* Add to sum of squares of the window function the parts of window which are equal to 1, and then find RMS value*/
  winFncRMS += ( REAL8 )( N - 2 * WINEND );
  winFncRMS = sqrt( ( winFncRMS / ( ( REAL8 ) N ) ) );

  return 0;
}

/* Hann window based on Matlab, but with C indexing: w[k] = 0.5*( 1 - cos(2*pi*k/(N-1)) ) k = 0, 1, 2,...N-1 */
int WindowDataHann( void )
{
  INT4 k;
  REAL8 win, N, Nm1;
  REAL8 real8TwoPi = 2.0 * ( ( REAL8 )( LAL_PI ) );

  /* initialize the RMS of the window function */
  winFncRMS = 0.0;

  N = ( ( REAL8 )dataDouble.data->length );
  Nm1 = N - 1;
  for ( k = 0; k < N; k++ ) {
    win = 0.5 * ( 1.0 - cos( real8TwoPi * ( ( REAL8 )( k ) ) / Nm1 ) );
    dataDouble.data->data[k] *= win;
    winFncRMS += win * win;
  }

  /* Find RMS value; note that N is REAL8 in this function */
  winFncRMS = sqrt( ( winFncRMS / N ) );

  return 0;
}
