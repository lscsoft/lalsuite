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

/***************************** <lalVerbatim file="InjectTimeSeriesCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{InjectTimeSeries.c}}
\label{ss:InjectTimeSeries.c}

Injects a time series of floating-point numbers into a time series of
integers, with dithering.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{InjectTimeSeriesCP}
\idx{LALSI2InjectTimeSeries()}
\idx{LALSSInjectTimeSeries()}

\subsubsection*{Description}

The function \verb@LALSI2InjectTimeSeries()@ (i.e.\ ``Single-precision
to \verb@INT2@'') dithers each sample in \verb@*output@, adds the
nearest time sample from \verb@*signalvec@, and rounds to the nearest
integer, storing the result back in \verb@*output@.  If desired, the
random parameters for the dithering can be created outside this
routine and passed in as \verb@*params@ (see \verb@Random.h@); if this
pointer is \verb@NULL@, the parameters will be generated internally.

The function \verb@LALSSInjectVector()@ (i.e.\ ``Single-precision to
single-precision'') simply takes each sample from \verb@*output@ and
adds the nearest corresponding time sample from \verb@*signalvec@,
without performing any dithering.

\subsubsection*{Algorithm}

The algorithm is as given in \verb@InjectVector.c@, with the following
additional considerations.  Since the two time series each carry their
own information about epoch and sampling interval, the value to be
injected at a given point in \verb@*output@ is found by taking the
nearest time sample in \verb@*signalvec@.  Injection is only performed
over the range in times that \verb@*output@ and \verb@*signalvec@
overlap; other values in \verb@*output@ are untouched.

Previous versions of this algorithm found the value to be injected by
interpolating the two nearest samples in \verb@*signalvec@, which reduces
high-frequency aliasing noise and ensures that the pre- and
post-injection signals agree in timing to within a fraction of a
sample.  However, this interpolation effectively convolved the signal
with a triangular function of width $2\Delta t$, where $\Delta t$ is
the sampling interval of the \emph{signal}.  This has the effect of a
low-pass filter with an attenuation factor of $\sim0.8$ at frequencies
$\sim1/4\Delta t$.  Since input signals are typically sampled at or
near their Nyquist frequencies, this would represent an unacceptable
level of attenuation.  For this reason, the current version of the
algorithm eliminates the interpolation procedure.


\subsubsection*{Uses}
\begin{verbatim}
LALCreateRandomParams()
LALDestroyRandomParams()
LALUniformDeviate()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{InjectTimeSeriesCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/Units.h>
#include <lal/Random.h>
#include <lal/Inject.h>

NRCSID( INJECTTIMESERIESC, "$Id$" );

/* <lalVerbatim file="InjectTimeSeriesCP"> */
void
LALSI2InjectTimeSeries( LALStatus       *stat,
			INT2TimeSeries  *output,
			REAL4TimeSeries *signalvec,
			RandomParams    *params )
{ /* </lalVerbatim> */
  INT4 n;  /* 1 + highest index of output touched by the injection */
  INT4 i;  /* index over output data */
  INT2 *outData; /* pointer to output->data->data */
  REAL8 dt;      /* output->deltaT in units of signalvec->deltaT */
  REAL8 offset;  /* the time from the start of *signalvec to the start
		     of *output, in units of signalvec->deltaT. */
  RandomParams *internal = NULL; /* internal random parameters */

  const INT2 max = 32767;  /* largest INT2 */
  const INT2 min = -32768; /* smallest INT2 */

  INITSTATUS( stat, "LALSI2InjectTimeSeries", INJECTTIMESERIESC );
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
    BOOLEAN unitsOK;
    LALUnitPair pair;

    pair.unitOne = &(signalvec->sampleUnits);
    pair.unitTwo = &lalADCCountUnit;
    TRY( LALUnitCompare( stat->statusPtr, &unitsOK, &pair ), stat );
    ASSERT( unitsOK, stat, INJECTH_EUNIT, INJECTH_MSGEUNIT );
    pair.unitOne = &(output->sampleUnits);
    TRY( LALUnitCompare( stat->statusPtr, &unitsOK, &pair ), stat );
    ASSERT( unitsOK, stat, INJECTH_EUNIT, INJECTH_MSGEUNIT );
    LALSnprintf( newName, LALNameLength, "%s plus %s", output->name,
		 signalvec->name );
    memcpy( output->name, newName, LALNameLength*sizeof(CHAR) );
  }

  /* If params = NULL, generate an internal set of parameters. */
  if ( !params )
    TRY( LALCreateRandomParams( stat->statusPtr, &internal, 0 ),
	 stat );
  else
    internal = params;

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
    REAL4 x = (REAL4)( outData[i] ); /* current output sample */
    REAL4 d;                         /* current dithering */

    /* Interpolate the signal.  ***REMOVED***
    REAL8 t = offset + i*dt;          interpolated signal index
    INT4 j = (INT4)floor( t );        signal index preceding t
    REAL4 frac = (REAL4)( t - j );    interpolation fraction
    REAL4 y = frac*(signalvec->data->data[j+1]) +    interpolated signal
      ( 1.0 - frac )*(signalvec->data->data[j]);     value */

    /* Extract the nearest signal sample. */
    INT4 j = (INT4)floor( offset + i*dt + 0.5 );
    REAL4 y = signalvec->data->data[j];

    /* Compute the dithering. */
    LALUniformDeviate( stat->statusPtr, &d, internal );
    BEGINFAIL( stat )
      if ( !params ) {
	TRY( LALDestroyRandomParams( stat->statusPtr, &internal ),
	     stat );
      }
    ENDFAIL( stat );

    /* Dither and inject. */
    x += d + y;
    if ( x > max )
      outData[i] = max;
    else if ( x < min )
      outData[i] = min;
    else
      outData[i] = (INT2)( floor( x ) );
  }

  /* Cleanup and exit. */
  if ( !params ) {
    TRY( LALDestroyRandomParams( stat->statusPtr, &internal ), stat );
  }
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="InjectTimeSeriesCP"> */
void
LALSSInjectTimeSeries( LALStatus       *stat,
		       REAL4TimeSeries *output,
		       REAL4TimeSeries *signalvec )
{ /* </lalVerbatim> */
  INT4 n;  /* 1 + highest index of output touched by the injection */
  INT4 i;  /* index over output data */
  REAL4 *outData; /* pointer to output->data->data */
  REAL8 dt;       /* output->deltaT in units of signalvec->deltaT */
  REAL8 offset;   /* the time from the start of *signalvec to the start
		     of *output, in units of signalvec->deltaT. */

  INITSTATUS( stat, "LALSSInjectTimeSeries", INJECTTIMESERIESC );
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
    BOOLEAN unitsOK;
    LALUnitPair pair;

    pair.unitOne = &(signalvec->sampleUnits);
    pair.unitTwo = &lalADCCountUnit;
    TRY( LALUnitCompare( stat->statusPtr, &unitsOK, &pair ), stat );
    ASSERT( unitsOK, stat, INJECTH_EUNIT, INJECTH_MSGEUNIT );
    pair.unitOne = &(output->sampleUnits);
    TRY( LALUnitCompare( stat->statusPtr, &unitsOK, &pair ), stat );
    ASSERT( unitsOK, stat, INJECTH_EUNIT, INJECTH_MSGEUNIT );
    LALSnprintf( newName, LALNameLength, "%s plus %s", output->name,
		 signalvec->name );
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
