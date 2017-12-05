/*
*  Copyright (C) 2007 Jolien Creighton, Teviet Creighton, John Whelan
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
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/Units.h>
#include <lal/Random.h>
#include <lal/Inject.h>

/**
 * \author Creighton, T. D.
 * \addtogroup InjectTimeSeries_c
 *
 * \brief Injects a time series of floating-point numbers into a time
 * series of floating-point numbers using nearest-neighbour interpolation.
 *
 * ### Algorithm ###
 *
 * Samples are added to the target time series.  If the timestamps of
 * samples in the two time series are not identical, the offsets in the
 * target time series to which samples are added are rounded to the nearest
 * integer --- no sub-sample interpolation is performed.
 */
/*@{*/


/** \see See documentation in \ref InjectTimeSeries_c */
void
LALSSInjectTimeSeries( LALStatus       *stat,
		       REAL4TimeSeries *output,
		       REAL4TimeSeries *signalvec )
{
  INT4 n;  /* 1 + highest index of output touched by the injection */
  INT4 i;  /* index over output data */
  REAL4 *outData; /* pointer to output->data->data */
  REAL8 dt;       /* output->deltaT in units of signalvec->deltaT */
  REAL8 offset;   /* the time from the start of *signalvec to the start
		     of *output, in units of signalvec->deltaT. */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* Make sure parameter structures and their fields exist. */
  ASSERT( signalvec, stat, INJECTH_ENUL, INJECTH_MSGENUL );
  ASSERT( signalvec->data, stat, INJECTH_ENUL, INJECTH_MSGENUL );
  ASSERT( signalvec->data->data, stat, INJECTH_ENUL, INJECTH_MSGENUL );
  ASSERT( output, stat, INJECTH_ENUL, INJECTH_MSGENUL );
  ASSERT( output->data, stat, INJECTH_ENUL, INJECTH_MSGENUL );
  ASSERT( output->data->data, stat, INJECTH_ENUL, INJECTH_MSGENUL );
  outData = output->data->data;

  /* Make sure we never divide by zero. */
  ASSERT( signalvec->deltaT != 0.0, stat, INJECTH_EBAD, INJECTH_MSGEBAD );
  dt = output->deltaT / signalvec->deltaT;
  ASSERT( dt != 0.0, stat, INJECTH_EBAD, INJECTH_MSGEBAD );

  /* Check dimensions. */
  {
    CHAR newName[LALNameLength];

    ASSERT( XLALUnitCompare( &(signalvec->sampleUnits), &lalADCCountUnit ) == 0, stat, INJECTH_EUNIT, INJECTH_MSGEUNIT );
    ASSERT( XLALUnitCompare( &(output->sampleUnits), &lalADCCountUnit ) == 0, stat, INJECTH_EUNIT, INJECTH_MSGEUNIT );
    if(snprintf( newName, LALNameLength, "%s plus %s", output->name,
		 signalvec->name ) >= LALNameLength)
      ABORT( stat, INJECTH_ENUL, INJECTH_MSGENUL );
    memcpy( output->name, newName, LALNameLength*sizeof(CHAR) );
  }

  /* Compute offset. */
  offset = ( output->epoch.gpsSeconds - signalvec->epoch.gpsSeconds ) /
    signalvec->deltaT;
  offset += ( output->epoch.gpsNanoSeconds -
	      signalvec->epoch.gpsNanoSeconds ) * 1.0e-9 /
    signalvec->deltaT;

  /* Compute initial value of i, and correct to ensure we will never
     index either array out of its bounds. */
  i = (INT4)( -offset / dt );
  if ( i < 0 )
    i = 0;
  while ( offset + i*dt < 0.0 )
    i++;
  if ( i >= (INT4)( output->data->length ) )
    LALWarning( stat, "Signal starts after the end of the output"
		" time series." );

  /* Compute final value of i+1, and correct to ensure we will never
     index either array out of its bounds. */
  n = (INT4)( ( signalvec->data->length - offset ) / dt );
  if ( n > (INT4)( output->data->length ) )
    n = output->data->length;
  while ( offset + n*dt > signalvec->data->length )
    n--;
  if ( n <= 0 )
    LALWarning( stat, "Signal ends before the start of the output"
		" time series." );

  /* Start injecting... */
  for ( ; i < n; i++ ) {

    /* Interpolate the signal.  ***REMOVED***
    REAL8 t = offset + i*dt;          interpolated signal index
    INT4 j = (INT4)floor( t );        signal index preceding t
    REAL4 frac = (REAL4)( t - j );    interpolation fraction
    REAL4 y = frac*(signalvec->data->data[j+1]) +    interpolated signal
      ( 1.0 - frac )*(signalvec->data->data[j]);     value */

    /* Extract the nearest signal sample. */
    INT4 j = (INT4)floor( offset + i*dt + 0.5 );
    REAL4 y = signalvec->data->data[j];

    /* Add the signal to the output. */
    outData[i] += y;
  }

  /* Exit. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}

/*@}*/
