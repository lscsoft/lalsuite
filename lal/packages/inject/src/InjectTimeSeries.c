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

\subsubsection*{Description}

This function dithers the contents of \verb@*output@, adds the
contents of \verb@*signal@, and rounds to the nearest integer, storing
the result back in \verb@*output@.  If desired, the random parameters
for the dithering can be created outside this routine and passed in as
\verb@*params@ (see \verb@Random.h@); if this pointer is \verb@NULL@,
the parameters will be generated internally.

The characters \verb@SI2@ refer to ``Single-precision (i.e.\
\verb@REAL4@) to \verb@INT2@''.  At present I see no need for other
injection routines, since LIGO data will be 2-byte integers, and since
the dynamic range of a \verb@REAL8@ is completely inappropriate for
injection into integer data.  However, the namespace convention does
allow for other routines to be added.

\subsubsection*{Algorithm}

The algorithm is as given in \verb@InjectVector.c@, with the following
additional considerations.  Since the two time series each carry their
own information about epoch and sampling interval, the value to be
injected at a given point in \verb@*output@ is found by interpolating
the two nearest time samples in \verb@*signal@.  Injection is only
performed over the range in times that \verb@*output@ and
\verb@*signal@ overlap; other values in \verb@*output@ are untouched.

An unfortunate side effect of the interpolation is that it effectively
convolves the signal with a triangular function of width $2\Delta t$,
where $\Delta t$ is the sampling interval of the signal.  This
resultis in a low-pass filter with an attenuation factor of $\sim0.8$
at frequencies $\sim1/4\Delta t$.  There are two ways around this:
\begin{enumerate}
\item Sample your waveform at a rate at least 8 times its highest
frequency (i.e.\ $4\times$ the Nyquist rate).  This is the preferred
method.
\item If your waveform is sampled at exactly the same sampling
interval as the data into which it will be injected, you should adjust
the epochs to ensure that the waveform samples are precisely aligned
with the data samples.  In this case, no interpolation will occur.
\end{enumerate}


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
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/Random.h>
#include <lal/Inject.h>

NRCSID( INJECTTIMESERIESC, "$Id$" );

/* <lalVerbatim file="InjectTimeSeriesCP"> */
void
LALSI2InjectTimeSeries( LALStatus       *stat,
			INT2TimeSeries  *output,
			REAL4TimeSeries *signal,
			RandomParams    *params )
{ /* </lalVerbatim> */
  INT4 n;  /* 1 + highest index of output touched by the injection */
  INT4 i;  /* index over output data */
  INT2 *outData; /* pointer to output->data->data */
  REAL8 dt;      /* output->deltaT in units of signal->deltaT */
  REAL8 offset;  /* the time from the start of *signal to the start
		     of *output, in units of signal->deltaT. */
  RandomParams *internal = NULL; /* internal random parameters */

  const INT2 max = 32767;  /* largest INT2 */
  const INT2 min = -32768; /* smallest INT2 */

  INITSTATUS( stat, "LALSI2InjectTimeSeries", INJECTTIMESERIESC );
  ATTATCHSTATUSPTR( stat );

  /* Make sure parameter structures and their fields exist. */
  ASSERT( signal, stat, INJECTH_ENUL, INJECTH_MSGENUL );
  ASSERT( signal->data, stat, INJECTH_ENUL, INJECTH_MSGENUL );
  ASSERT( signal->data->data, stat, INJECTH_ENUL, INJECTH_MSGENUL );
  ASSERT( output, stat, INJECTH_ENUL, INJECTH_MSGENUL );
  ASSERT( output->data, stat, INJECTH_ENUL, INJECTH_MSGENUL );
  ASSERT( output->data->data, stat, INJECTH_ENUL, INJECTH_MSGENUL );
  outData = output->data->data;

  /* Make sure we never divide by zero. */
  ASSERT( signal->deltaT != 0.0, stat, INJECTH_EBAD, INJECTH_MSGEBAD );
  dt = output->deltaT / signal->deltaT;
  ASSERT( dt != 0.0, stat, INJECTH_EBAD, INJECTH_MSGEBAD );

  /* If params = NULL, generate an internal set of parameters. */
  if ( !params )
    TRY( LALCreateRandomParams( stat->statusPtr, &internal, 0 ),
	 stat );
  else
    internal = params;

  /* Compute offset. */
  offset = ( output->epoch.gpsSeconds - signal->epoch.gpsSeconds ) /
    signal->deltaT;
  offset += ( output->epoch.gpsNanoSeconds -
	      signal->epoch.gpsNanoSeconds ) * 1.0e-9 /
    signal->deltaT;

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
  n = (INT4)( ( signal->data->length - offset ) / dt );
  if ( n > (INT4)( output->data->length ) )
    n = output->data->length;
  while ( offset + n*dt > signal->data->length )
    n--;
  if ( n <= 0 )
    LALWarning( stat, "Signal ends before the start of the output"
		" time series." );

  /* Start injecting... */
  for ( ; i < n; i++ ) {
    REAL4 x = (REAL4)( outData[i] ); /* current output sample */
    REAL4 d;                         /* current dithering */

    /* Interpolate the signal. */
    REAL8 t = offset + i*dt;       /* interpolated signal index */
    INT4 j = (INT4)floor( t );     /* signal index preceding t */
    REAL4 frac = (REAL4)( t - j ); /* interpolation fraction */
    REAL4 y = frac*(signal->data->data[j+1]) + /* interpolated signal */
      ( 1.0 - frac )*(signal->data->data[j]);  /* value */

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
