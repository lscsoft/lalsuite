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
 * \textcolor{red}{Danger: In fact this function is very unsafe. It is only a
 * placeholder until I get new code from Isabel.}
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

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/AVFactories.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/ResampleTimeSeries.h>

NRCSID( RESAMPLETIMESERIESC, "$Id$" );

#define ONEBY32K 0.000030517578125

/* <lalVerbatim file="ResampleTimeSeriesCP"> */
void
LALResampleREAL4TimeSeries(
    LALStatus          *status,
    REAL4TimeSeries    *ts, 
    ResampleTSParams   *params
    )
{ /* </lalVerbatim> */
  UINT4  i;
  int    ratio;
  REAL8  nyquist, ratioDeltaT;
  REAL4  *dataPtr;
  const REAL8 epsilon = 1.0e-8;
  
  INITSTATUS( status, "LALResampleREAL4TimeSeries", RESAMPLETIMESERIESC );

  ASSERT( ts, status, 
      RESAMPLETIMESERIESH_ENULL, RESAMPLETIMESERIESH_MSGENULL );
  ASSERT( ts->deltaT > 0, status,
      RESAMPLETIMESERIESH_ERATE, RESAMPLETIMESERIESH_MSGERATE );
  ASSERT( ts->data, status, 
      RESAMPLETIMESERIESH_ENULL, RESAMPLETIMESERIESH_MSGENULL );
  ASSERT( ts->data->length, status, 
      RESAMPLETIMESERIESH_EZERO, RESAMPLETIMESERIESH_MSGEZERO );
  ASSERT( ts->data->data, status, 
      RESAMPLETIMESERIESH_ENULL, RESAMPLETIMESERIESH_MSGENULL );
  ASSERT( params, status,
      RESAMPLETIMESERIESH_ENULL, RESAMPLETIMESERIESH_MSGENULL );
  ASSERT( params->deltaT > 0, status,
      RESAMPLETIMESERIESH_ERATE, RESAMPLETIMESERIESH_MSGERATE );

  /* if the input and output rates are the same do nothing */
  if ( fabs( ts->deltaT - params->deltaT ) < epsilon )
  {
    LALWarning( status, "no resampling performed:"
        " input and output sample rates are identical" );
    RETURN( status );
  }

  ATTATCHSTATUSPTR( status );

  /* check that we are downsampling */
  if ( ts->deltaT > params->deltaT )
  {
    ABORT( status, RESAMPLETIMESERIESH_EUPSM, RESAMPLETIMESERIESH_MSGEUPSM );
  }

  /* check that the input sample rate is below 32kHz */
  if ( ts->deltaT < ONEBY32K )
  {
    ABORT( status, RESAMPLETIMESERIESH_EHIGH, RESAMPLETIMESERIESH_MSGEHIGH );
  }

  /* check that we are resampling by a power of two */
  ratioDeltaT = params->deltaT / ts->deltaT;
  ratio = (int) rint( ratioDeltaT );
  if ( ! ( ratio==0x1 || ratio==0x2 || ratio==0x4 || ratio==0x8 ||
        ratio==0x10 || ratio==0x20 || ratio==0x40 || ratio==0x80 ) )
  {
    ABORT( status, RESAMPLETIMESERIESH_ELOG2, RESAMPLETIMESERIESH_MSGELOG2 );
  }

  /* get the new nyquist frequency to use in the filter */
  nyquist = 1.0 / ( 2.0 * (REAL8) ratio * ts->deltaT );
  fprintf( stderr, "nyquist = %e\n", nyquist );

  if ( params->filterType == defaultButterworth )
  {
    /* a butterworth filter with some reasonable params */
    params->filterParams.butterworth.nMax = 4;
    params->filterParams.butterworth.f1 = nyquist - 100.0;
    params->filterParams.butterworth.a1 = 0.1;
    params->filterParams.butterworth.f2 = 0;
    params->filterParams.butterworth.a2 = 0;
    
    fprintf( stdout, "attn %e at %e\n", 
        params->filterParams.butterworth.a1,
        params->filterParams.butterworth.f1 );

    LALDButterworthREAL4TimeSeries( status->statusPtr, ts, 
        &(params->filterParams.butterworth) );
    CHECKSTATUSPTR( status );
  }
  else
  {
    ABORT( status, RESAMPLETIMESERIESH_EFILT, RESAMPLETIMESERIESH_MSGEFILT );
  }

  /* decimate the time series */
  ts->deltaT = params->deltaT;
  ts->data->length /= ratio;
  dataPtr = ts->data->data;
  for ( i = 0; i < ts->data->length; ++i )
  {
    ts->data->data[i] = *dataPtr;
    dataPtr += ratio;
  }
  LALRealloc( ts->data->data, ts->data->length * sizeof(REAL4) );

  DETATCHSTATUSPTR( status );
  RETURN( status );

}

#undef ONEBY32K
