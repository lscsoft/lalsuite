/**** <lalVerbatim file="ResampleTimeSeriesCV">
 * Author: Brown, D. A.
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \subsection{Module \texttt{ResampleTimeSeries.c}}
 * \label{ss:ResampleTimeSeries.c}
 *
 * The module contains functions to resampe a time series.
 *
 * \textcolor{red}{Warning: these functions have not yet been tested}
 *
 * \subsection*{Prototypes}
 * \input{ResampleTimeSeriesCP}
 *
 * The routine \texttt{LALResampleREAL4TimeSertes()} resample the data in
 * the specifeid time series \texttt{data} to the sample rate requested in
 * \texttt{params}. Resampling is done in place.
 *
 * \vfill{\footnotesize\input{ResampleTimeSeriesCV}}
 *
 **** </lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/AVFactories.h>
#include <lal/ResampleTimeSeries.h>

NRCSID( RESAMPLETIMESERIESC, "$Id$" );

/* <lalVerbatim file="ResampleTimeSeriesCP"> */
void
LALResampleREAL4TimeSeries(
    LALStatus                  *status,
    REAL4TimeSeries            *data, 
    ResampleTimeSeriesParams   *params
    )
{ /* </lalVerbatim> */

  INITSTATUS( status, "LALResampleREAL4TimeSeries", RESAMPLETIMESERIESC );
  ATTATCHSTATUSPTR( status );

  ASSERT( data, status, 
      RESAMPLETIMESERIESH_ENULL, RESAMPLETIMESERIESH_MSGENULL );
  ASSERT( data->data, status, 
      RESAMPLETIMESERIESH_ENULL, RESAMPLETIMESERIESH_MSGENULL );
  ASSERT( data->data->length, status, 
      RESAMPLETIMESERIESH_EZERO, RESAMPLETIMESERIESH_MSGEZERO );
  ASSERT( data->data->data, status, 
      RESAMPLETIMESERIESH_ENULL, RESAMPLETIMESERIESH_MSGENULL );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
