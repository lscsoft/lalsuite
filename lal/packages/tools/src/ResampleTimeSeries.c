/**** <lalVerbatim file="ResampleTimeSeriesCV">
 * Author: Brown, D. A., Brady, P. R. and Fritschel, P. F.
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \subsection{Module \texttt{ResampleTimeSeries.c}}
 * \label{ss:ResampleTimeSeries.c}
 *
 * The module contains functions to resample time series. 
 *
 * \textcolor{red}{Warning: these functions have not yet been tested}
 *
 * \textcolor{red}{Danger: In fact this function is very unsafe. It is only a
 * placeholder until I get new code from Isabel.}
 *
 * \subsection*{Prototypes}
 * \input{ResampleTimeSeriesCP}
 *
 * The routine \texttt{LALDecimateREAL4TimeSertes()} resample the data in
 * the specifeid time series \texttt{data} to the sample rate requested in
 * \texttt{params}. Resampling is done in place.  The downsampling
 * factor must be a power of 2. 
 *
 * The code is based on decimate functions in the GDS source tree which was
 * written by Peter Fritschel in 1998.  The current implementation is
 * intended to be called only once.   The data is low-pass filtered
 * using FIR filters from the same source.
 *
 *
 * \vfill{\footnotesize\input{ResampleTimeSeriesCV}}
 *
 **** </lalLaTeX> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/ResampleTimeSeries.h>

NRCSID( RESAMPLETIMESERIESC, "$Id$" );

#define ONEBY32K 0.000030517578125

/**********************************************************************
 * 
 * Filters taken from the GDS code "decimate.c" written by P. Fritschel 
 *
 *********************************************************************/
   static const float firls1[11] = 
   {2.225549e-3, -3.914217e-3, 6.226653e-3,
    -9.330331e-3, 1.347478e-2, -1.907462e-2,
    2.690526e-2, -3.863524e-2, 5.863624e-2,
    -1.030298e-1, 3.172755e-1};
   static const float firPM1[11] = 
   {5.704943e-3, -5.292245e-3, 7.672972e-3,
    -1.083958e-2, 1.493768e-2, -2.047472e-2,
    2.813003e-2, -3.967826e-2, 5.939968e-2,
    -1.035220e-1, 3.174278e-1};
   static const float firls2[6] =
   {-1.353431e-2, 2.193691e-2, -3.448320e-2,
    5.550899e-2, -1.010866e-1, 3.166165e-1};
   static const float firls3[21] =
   {8.549310e-5, -1.882289e-4, 3.483209e-4,
    -5.835810e-4, 9.146788e-4, -1.365369e-3,
    1.962879e-3, -2.738601e-3, 3.729302e-3,
    -4.979237e-3, 6.543786e-3, -8.495763e-3,
    1.093660e-2, -1.401691e-2, 1.797646e-2,
    -2.322807e-2, 3.055373e-2, -4.163686e-2,
    6.087163e-2, -1.044086e-1, 3.177415e-1};

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
  fprintf( stderr, "nyquist = %e %i\n", nyquist, ratio );

  if ( params->filterType == defaultButterworth )
  {
    /* a butterworth filter with some reasonable params */
    params->filterParams.butterworth.nMax = 8;
    params->filterParams.butterworth.f1 = nyquist-10.0;
    params->filterParams.butterworth.a1 = 0.99;
    params->filterParams.butterworth.f2 = nyquist;
    params->filterParams.butterworth.a2 = 0.01;
    
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


/* <lalVerbatim file="ResampleTimeSeriesCP"> */
void
LALDecimateREAL4TimeSeries(
    LALStatus          *status,
    REAL4TimeSeries    *ts, 
    ResampleTSParams   *params
    )
{ /* </lalVerbatim> */
    UINT4  i;
    INT4   ratio;
    REAL8  nyquist, ratioDeltaT, delay;
    REAL4  *dataPtr;
    const REAL8         epsilon = 1.0e-8;
    const REAL4*	filt_coeff;
    REAL8  mantissa;
    REAL4  sum;
    INT4   filt_length;
    INT4   tempnpt;
    INT4   npt;
    INT4   filt_ord;
    INT4   j;
    INT4   k;
    INT4   npt_y;
    INT4   ndec;

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
  fprintf( stderr, "nyquist = %e %i\n", nyquist, ratio );

  /* select the FIR filter to use for low-pass */
  if ( params->filterType == firLSOne )
  {
      filt_coeff = firls1;
      filt_length = 11;
  }
  else if ( params->filterType == firLSTwo )
  {
      filt_coeff = firls2;
      filt_length = 6;
  }
  else if ( params->filterType == firLSThree )
  {
      filt_coeff = firls3;
      filt_length = 21;
  }
  else if ( params->filterType == firEqOne )
  {
      filt_coeff = firPM1;
      filt_length = 11;
  }
  else
  {
    ABORT( status, RESAMPLETIMESERIESH_EFILT, RESAMPLETIMESERIESH_MSGEFILT );
  }

  /* BEGIN code fragment modified from GDS tree                  */
  /*  check that dec_factor is a power of 2; return error if not */
  if ((mantissa = frexp((double) ratio, &ndec)) != 0.5)
  {
      ABORT( status, RESAMPLETIMESERIESH_ELOG2, RESAMPLETIMESERIESH_MSGELOG2 );
  }
  --ndec;			       /* this makes 2^(ndec) = dec_factor   */
  npt_y = (ts->data->length/ratio);    /* # pts in the output array   */
  filt_ord = 4*filt_length - 2;        /* filter order  */
  tempnpt = ndec*filt_ord;             /* #pts needed in temp array to */

  /* filter algorithm; half-band, linear phase FIR filter;	*/
  /* taking advantage of the symmetry of the filter		*/
  /* coefficients, the appropriate data pairs are added	        */
  /* together, then multiplied by the filter coefficients	*/

  npt = ts->data->length - tempnpt;
  dataPtr = ts->data->data;
  delay = 0.0;
  for (k = 1; k <= ndec; ++k) {
      delay += pow(2.0, k-1);
      npt /= 2;     /* decrease #pts by factor of 2 at each stage */
      for (i = 0; i < npt; ++i) {	
          sum = 0.0;
          for (j = 0; j < filt_length; j++) {
              sum += filt_coeff[j] *
                  (dataPtr[2*(i+j)] + dataPtr[filt_ord + 2*(i-j)]);
          }
          dataPtr[i] = sum + dataPtr[filt_ord/2 + 2*i]/2.0;
      }
  }

  delay = delay * 0.5 * filt_ord;

  printf("data shuffling done %e %e\n",ts->deltaT * delay, delay);
  ts->deltaT = params->deltaT;
  ts->data->length /= ratio;
  LALRealloc( ts->data->data, ts->data->length * sizeof(REAL4) );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}




#undef ONEBY32K
