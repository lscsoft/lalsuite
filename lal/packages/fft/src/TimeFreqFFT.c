/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Patrick Brady, Reinhard Prix, Tania Regimbau, John Whelan
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


#include <complex.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/TimeFreqFFT.h>
#include <lal/LALConstants.h>
#include <lal/RngMedBias.h>

/* ---------- see TimeFreqFFT.h for doxygen documentation */

/*
 *
 * XLAL REAL4 Time->Freq and Freq->Time FFT routines
 *
 */


int XLALREAL4TimeFreqFFT(
    COMPLEX8FrequencySeries *freq,
    const REAL4TimeSeries   *time,
    const REAL4FFTPlan      *plan
    )
{
  UINT4 k;

  if ( ! freq || ! time || ! plan )
    XLAL_ERROR( XLAL_EFAULT );
  if ( time->deltaT <= 0.0 )
    XLAL_ERROR( XLAL_EINVAL );

  /* perform the transform */
  if ( XLALREAL4ForwardFFT( freq->data, time->data, plan ) == XLAL_FAILURE )
    XLAL_ERROR( XLAL_EFUNC );

  /* adjust the units */
  if ( ! XLALUnitMultiply( &freq->sampleUnits, &time->sampleUnits, &lalSecondUnit ) )
    XLAL_ERROR( XLAL_EFUNC );

  /* remaining fields */
  if ( time->f0 )  /* TODO: need to figure out what to do here */
    XLALPrintWarning( "XLAL Warning - frequency series may have incorrect f0" );
  freq->f0     = 0.0; /* FIXME: what if heterodyned data? */
  freq->epoch  = time->epoch;
  freq->deltaF = 1.0 / ( time->deltaT * time->data->length );

  /* provide the correct scaling of the result */
  for ( k = 0; k < freq->data->length; ++k )
    freq->data->data[k] *= time->deltaT;

  return 0;
}


int XLALREAL4FreqTimeFFT(
    REAL4TimeSeries               *time,
    const COMPLEX8FrequencySeries *freq,
    const REAL4FFTPlan            *plan
    )
{
  UINT4 j;

  if ( ! freq || ! time || ! plan )
    XLAL_ERROR( XLAL_EFAULT );
  if ( freq->deltaF <= 0.0 )
    XLAL_ERROR( XLAL_EINVAL );

  /* perform the transform */
  if ( XLALREAL4ReverseFFT( time->data, freq->data, plan ) == XLAL_FAILURE )
    XLAL_ERROR( XLAL_EFUNC );

  /* adjust the units */
  if ( ! XLALUnitMultiply( &time->sampleUnits, &freq->sampleUnits, &lalHertzUnit ) )
    XLAL_ERROR( XLAL_EFUNC );

  /* remaining fields */
  if ( freq->f0 )  /* TODO: need to figure out what to do here */
    XLALPrintWarning( "XLAL Warning - time series may have incorrect f0" );
  time->f0     = 0.0; /* FIXME: what if heterodyned data? */
  time->epoch  = freq->epoch;
  time->deltaT = 1.0 / ( freq->deltaF * time->data->length );

  /* provide the correct scaling of the result */
  for ( j = 0; j < time->data->length; ++j )
    time->data->data[j] *= freq->deltaF;

  return 0;
}


/*
 *
 * XLAL REAL8 Time->Freq and Freq->Time FFT routines
 *
 */


int XLALREAL8TimeFreqFFT(
    COMPLEX16FrequencySeries *freq,
    const REAL8TimeSeries    *time,
    const REAL8FFTPlan       *plan
    )
{
  UINT4 k;

  if ( ! freq || ! time || ! plan )
    XLAL_ERROR( XLAL_EFAULT );
  if ( time->deltaT <= 0.0 )
    XLAL_ERROR( XLAL_EINVAL );

  /* perform the transform */
  if ( XLALREAL8ForwardFFT( freq->data, time->data, plan ) == XLAL_FAILURE )
    XLAL_ERROR( XLAL_EFUNC );

  /* adjust the units */
  if ( ! XLALUnitMultiply( &freq->sampleUnits, &time->sampleUnits, &lalSecondUnit ) )
    XLAL_ERROR( XLAL_EFUNC );

  /* remaining fields */
  if ( time->f0 )  /* TODO: need to figure out what to do here */
    XLALPrintWarning( "XLAL Warning - frequency series may have incorrect f0" );
  freq->f0     = 0.0; /* FIXME: what if heterodyned data? */
  freq->epoch  = time->epoch;
  freq->deltaF = 1.0 / ( time->deltaT * time->data->length );

  /* provide the correct scaling of the result */
  for ( k = 0; k < freq->data->length; ++k )
    freq->data->data[k] *= time->deltaT;

  return 0;
}


int XLALREAL8FreqTimeFFT(
    REAL8TimeSeries                *time,
    const COMPLEX16FrequencySeries *freq,
    const REAL8FFTPlan             *plan
    )
{
  UINT4 j;

  if ( ! freq || ! time || ! plan )
    XLAL_ERROR( XLAL_EFAULT );
  if ( freq->deltaF <= 0.0 )
    XLAL_ERROR( XLAL_EINVAL );

  /* perform the transform */
  if ( XLALREAL8ReverseFFT( time->data, freq->data, plan ) == XLAL_FAILURE )
    XLAL_ERROR( XLAL_EFUNC );

  /* adjust the units */
  if ( ! XLALUnitMultiply( &time->sampleUnits, &freq->sampleUnits, &lalHertzUnit ) )
    XLAL_ERROR( XLAL_EFUNC );

  /* remaining fields */
  if ( freq->f0 )  /* TODO: need to figure out what to do here */
    XLALPrintWarning( "XLAL Warning - time series may have incorrect f0" );
  time->f0     = 0.0; /* FIXME: what if heterodyned data? */
  time->epoch  = freq->epoch;
  time->deltaT = 1.0 / ( freq->deltaF * time->data->length );

  /* provide the correct scaling of the result */
  for ( j = 0; j < time->data->length; ++j )
    time->data->data[j] *= freq->deltaF;

  return 0;
}


/*
 *
 * XLAL COMPLEX8 Time->Freq and Freq->Time FFT routines
 *
 */

/* We need to define the COMPLEX8FFTPlan structure so that we can check the FFT
 * sign.  The plan is not really a void*, but it is a pointer so a void* is
 * good enough (we don't need it).  */
struct tagCOMPLEX8FFTPlan
{
  INT4  sign;
  UINT4 size;
  void *junk;
};


int XLALCOMPLEX8TimeFreqFFT(
    COMPLEX8FrequencySeries  *freq,
    const COMPLEX8TimeSeries *time,
    const COMPLEX8FFTPlan    *plan
    )
{
  COMPLEX8Vector *tmp;
  UINT4 k;

  if ( ! freq || ! time || ! plan )
    XLAL_ERROR( XLAL_EFAULT );
  if ( time->deltaT <= 0.0 )
    XLAL_ERROR( XLAL_EINVAL );
  if ( plan->sign != -1 )
    XLAL_ERROR( XLAL_EINVAL );

  /* create temporary workspace */
  tmp = XLALCreateCOMPLEX8Vector( time->data->length );
  if ( ! tmp )
    XLAL_ERROR( XLAL_EFUNC );

  /* perform transform */
  if ( XLALCOMPLEX8VectorFFT( tmp, time->data, plan ) == XLAL_FAILURE )
  {
    int saveErrno = xlalErrno;
    xlalErrno = 0;
    XLALDestroyCOMPLEX8Vector( tmp );
    xlalErrno = saveErrno;
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* unpack the frequency series and multiply by deltaT */
  for ( k = 0; k < time->data->length / 2; ++k )
    freq->data->data[k] = time->deltaT * tmp->data[k+(time->data->length+1)/2];
  for ( k = time->data->length / 2; k < time->data->length; ++k )
    freq->data->data[k] = time->deltaT * tmp->data[k-time->data->length/2];

  /* destroy temporary workspace */
  XLALDestroyCOMPLEX8Vector( tmp );
  if ( xlalErrno )
    XLAL_ERROR( XLAL_EFUNC );

  /* adjust the units */
  if ( ! XLALUnitMultiply( &freq->sampleUnits, &time->sampleUnits, &lalSecondUnit ) )
    XLAL_ERROR( XLAL_EFUNC );

  /* remaining fields */
  freq->epoch  = time->epoch;
  freq->deltaF = 1.0 / ( time->deltaT * time->data->length );
  freq->f0     = time->f0 - freq->deltaF * floor( time->data->length / 2 );

  return 0;
}


int XLALCOMPLEX8FreqTimeFFT(
    COMPLEX8TimeSeries            *time,
    const COMPLEX8FrequencySeries *freq,
    const COMPLEX8FFTPlan         *plan
    )
{
  COMPLEX8Vector *tmp;
  UINT4 k;

  if ( ! freq || ! time || ! plan )
    XLAL_ERROR( XLAL_EFAULT );
  if ( freq->deltaF <= 0.0 )
    XLAL_ERROR( XLAL_EINVAL );
  if ( plan->sign != 1 )
    XLAL_ERROR( XLAL_EINVAL );

  /* create temporary workspace */
  tmp = XLALCreateCOMPLEX8Vector( freq->data->length );
  if ( ! tmp )
    XLAL_ERROR( XLAL_EFUNC );

  /* pack the frequency series and multiply by deltaF */
  for ( k = 0; k < freq->data->length / 2; ++k )
    tmp->data[k+(freq->data->length+1)/2] = freq->deltaF * freq->data->data[k];
  for ( k = freq->data->length / 2; k < freq->data->length; ++k )
    tmp->data[k-freq->data->length/2] = freq->deltaF * freq->data->data[k];

  /* perform transform */
  if ( XLALCOMPLEX8VectorFFT( time->data, tmp, plan ) == XLAL_FAILURE )
  {
    int saveErrno = xlalErrno;
    xlalErrno = 0;
    XLALDestroyCOMPLEX8Vector( tmp );
    xlalErrno = saveErrno;
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* destroy temporary workspace */
  XLALDestroyCOMPLEX8Vector( tmp );
  if ( xlalErrno )
    XLAL_ERROR( XLAL_EFUNC );

  /* adjust the units */
  if ( ! XLALUnitMultiply( &time->sampleUnits, &freq->sampleUnits, &lalHertzUnit ) )
    XLAL_ERROR( XLAL_EFUNC );

  /* remaining fields */
  time->f0     = freq->f0 + freq->deltaF * floor( freq->data->length / 2 );
  time->epoch  = freq->epoch;
  time->deltaT = 1.0 / ( freq->deltaF * freq->data->length );

  return 0;
}


/*
 *
 * XLAL COMPLEX16 Time->Freq and Freq->Time FFT routines
 *
 */

/* We need to define the COMPLEX16FFTPlan structure so that we can check the FFT
 * sign.  The plan is not really a void*, but it is a pointer so a void* is
 * good enough (we don't need it).  */
struct tagCOMPLEX16FFTPlan
{
  INT4  sign;
  UINT4 size;
  void *junk;
};


int XLALCOMPLEX16TimeFreqFFT(
    COMPLEX16FrequencySeries  *freq,
    const COMPLEX16TimeSeries *time,
    const COMPLEX16FFTPlan    *plan
    )
{
  COMPLEX16Vector *tmp;
  UINT4 k;

  if ( ! freq || ! time || ! plan )
    XLAL_ERROR( XLAL_EFAULT );
  if ( time->deltaT <= 0.0 )
    XLAL_ERROR( XLAL_EINVAL );
  if ( plan->sign != -1 )
    XLAL_ERROR( XLAL_EINVAL );

  /* create temporary workspace */
  tmp = XLALCreateCOMPLEX16Vector( time->data->length );
  if ( ! tmp )
    XLAL_ERROR( XLAL_EFUNC );

  /* perform transform */
  if ( XLALCOMPLEX16VectorFFT( tmp, time->data, plan ) == XLAL_FAILURE )
  {
    int saveErrno = xlalErrno;
    xlalErrno = 0;
    XLALDestroyCOMPLEX16Vector( tmp );
    xlalErrno = saveErrno;
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* unpack the frequency series and multiply by deltaT */
  for ( k = 0; k < time->data->length / 2; ++k )
    freq->data->data[k] = time->deltaT * tmp->data[k+(time->data->length+1)/2];
  for ( k = time->data->length / 2; k < time->data->length; ++k )
    freq->data->data[k] = time->deltaT * tmp->data[k-time->data->length/2];

  /* destroy temporary workspace */
  XLALDestroyCOMPLEX16Vector( tmp );
  if ( xlalErrno )
    XLAL_ERROR( XLAL_EFUNC );

  /* adjust the units */
  if ( ! XLALUnitMultiply( &freq->sampleUnits, &time->sampleUnits, &lalSecondUnit ) )
    XLAL_ERROR( XLAL_EFUNC );

  /* remaining fields */
  freq->epoch  = time->epoch;
  freq->deltaF = 1.0 / ( time->deltaT * time->data->length );
  freq->f0     = time->f0 - freq->deltaF * floor( time->data->length / 2 );

  return 0;
}


int XLALCOMPLEX16FreqTimeFFT(
    COMPLEX16TimeSeries            *time,
    const COMPLEX16FrequencySeries *freq,
    const COMPLEX16FFTPlan         *plan
    )
{
  COMPLEX16Vector *tmp;
  UINT4 k;

  if ( ! freq || ! time || ! plan )
    XLAL_ERROR( XLAL_EFAULT );
  if ( freq->deltaF <= 0.0 )
    XLAL_ERROR( XLAL_EINVAL );
  if ( plan->sign != 1 )
    XLAL_ERROR( XLAL_EINVAL );

  /* create temporary workspace */
  tmp = XLALCreateCOMPLEX16Vector( freq->data->length );
  if ( ! tmp )
    XLAL_ERROR( XLAL_EFUNC );

  /* pack the frequency series and multiply by deltaF */
  for ( k = 0; k < freq->data->length / 2; ++k )
    tmp->data[k+(freq->data->length+1)/2] = freq->deltaF * freq->data->data[k];
  for ( k = freq->data->length / 2; k < freq->data->length; ++k )
    tmp->data[k-freq->data->length/2] = freq->deltaF * freq->data->data[k];

  /* perform transform */
  if ( XLALCOMPLEX16VectorFFT( time->data, tmp, plan ) == XLAL_FAILURE )
  {
    int saveErrno = xlalErrno;
    xlalErrno = 0;
    XLALDestroyCOMPLEX16Vector( tmp );
    xlalErrno = saveErrno;
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* destroy temporary workspace */
  XLALDestroyCOMPLEX16Vector( tmp );
  if ( xlalErrno )
    XLAL_ERROR( XLAL_EFUNC );

  /* adjust the units */
  if ( ! XLALUnitMultiply( &time->sampleUnits, &freq->sampleUnits, &lalHertzUnit ) )
    XLAL_ERROR( XLAL_EFUNC );

  /* remaining fields */
  time->f0     = freq->f0 + freq->deltaF * floor( freq->data->length / 2 );
  time->epoch  = freq->epoch;
  time->deltaT = 1.0 / ( freq->deltaF * freq->data->length );

  return 0;
}


/*
 *
 * LAL Real Time -> Freq and Freq -> Time routines
 * (single-precision only)
 *
 */



void
LALTimeFreqRealFFT(
    LALStatus               *status,
    COMPLEX8FrequencySeries *freq,
    REAL4TimeSeries         *time,
    RealFFTPlan             *plan
    )
{
  INITSTATUS(status);
  XLAL_PRINT_DEPRECATION_WARNING("XLALREAL4TimeFreqFFT");

  ASSERT( plan, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( freq, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( time, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( time->data, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( time->data->length, status,
      TIMEFREQFFTH_ESIZE, TIMEFREQFFTH_MSGESIZE );
  ASSERT( time->deltaT > 0, status, TIMEFREQFFTH_ERATE, TIMEFREQFFTH_MSGERATE );

  /* call the XLAL function */
  if ( XLALREAL4TimeFreqFFT( freq, time, plan ) == XLAL_FAILURE )
  {
    XLALClearErrno();
    ABORTXLAL( status );
  }

  RETURN( status );
}



void
LALFreqTimeRealFFT(
    LALStatus               *status,
    REAL4TimeSeries         *time,
    COMPLEX8FrequencySeries *freq,
    RealFFTPlan             *plan
    )
{
  INITSTATUS(status);
  XLAL_PRINT_DEPRECATION_WARNING("XLALREAL4FreqTimeFFT");
  ATTATCHSTATUSPTR( status );

  ASSERT( plan, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( freq, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( time, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( time->data, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( time->data->length, status,
      TIMEFREQFFTH_ESIZE, TIMEFREQFFTH_MSGESIZE );
  ASSERT( freq->deltaF > 0, status, TIMEFREQFFTH_ERATE, TIMEFREQFFTH_MSGERATE );

  /* call the XLAL function */
  if ( XLALREAL4FreqTimeFFT( time, freq, plan ) == XLAL_FAILURE )
  {
    XLALClearErrno();
    ABORTXLAL( status );
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/*
 *
 * median spectrum estimator based on mark's version
 *
 */


static REAL4
MedianSpec(
    REAL4      *p,
    REAL4      *s,
    UINT4       j,
    UINT4       flength,
    UINT4       numSegs
    )
{
  /* p points to array of power spectra data over time slices */
  /* s is a scratch array used to sort the values             */
  /* j is desired frequency offset into power spectra array   */
  /* flength is size of frequency series obtained from DFT    */
  /* numSegs is the number of time slices to be evaluated     */
  /* returns the median value, over time slice at given freq. */

  UINT4  outer  = 0;       /* local loop counter */
  UINT4  middle = 0;       /* local loop counter */
  UINT4  inner  = 0;       /* local loop counter */
  REAL4 returnVal = 0.0;  /* holder for return value */

  /* zero out the sort array */
  memset( s, 0, numSegs * sizeof(REAL4) );

  /* scan time slices for a given frequency */
  for ( outer = 0; outer < numSegs; ++outer )
  {
    /* insert power value into sort array */
    REAL4 tmp = p[outer * flength + j]; /* obtain value to insert */
    for ( middle = 0; middle < numSegs; ++middle )
    {
      if ( tmp > s[middle] )
      {
        /* insert taking place of s[middle] */
        for ( inner = numSegs - 1; inner > middle; --inner )
        {
          s[inner] = s [inner - 1];  /* move old values */
        }
        s[middle] = tmp;   /* insert new value */
        break;  /* terminate loop */
      }
    }
  }  /* done inserting into sort array */

  /* check for odd or even number of segments */
  if ( numSegs % 2 )
  {
    /* if odd number of segments, return median */
    returnVal = s[numSegs / 2];
  }
  else
  {
    /* if even number of segments, return average of two medians */
    returnVal = 0.5 * (s[numSegs/2] + s[(numSegs/2) - 1]);
  }

  return returnVal;
}


/*
 *
 * compute average power spectrum
 *
 */



void
LALREAL4AverageSpectrum (
    LALStatus                   *status,
    REAL4FrequencySeries        *fSeries,
    REAL4TimeSeries             *tSeries,
    AverageSpectrumParams       *params
    )

{
  UINT4                 i, j, k;          /* seg, ts and freq counters       */
  UINT4                 numSeg;           /* number of segments in average   */
  UINT4                 fLength;          /* length of requested power spec  */
  UINT4                 tLength;          /* length of time series segments  */
  REAL4Vector          *tSegment = NULL;  /* dummy time series segment       */
  COMPLEX8Vector       *fSegment = NULL;  /* dummy freq series segment       */
  REAL4                *tSeriesPtr;       /* pointer to the segment data     */
  REAL4                 psdNorm = 0;      /* factor to multiply windows data */
  REAL4                 fftRe, fftIm;     /* real and imag parts of fft      */
  REAL4                *s;                /* work space for computing mean   */
  REAL4                *psdSeg = NULL;    /* storage for individual specta   */
  LALUnit               unit;
  LALUnitPair           pair;
  /* RAT4                  negRootTwo = { -1, 1 }; */

  INITSTATUS(status);
  XLAL_PRINT_DEPRECATION_WARNING("XLALREAL4AverageSpectrumWelch");
  ATTATCHSTATUSPTR (status);

  /* check the input and output data pointers are non-null */
  ASSERT( fSeries, status,
      TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( fSeries->data, status,
      TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( fSeries->data->data, status,
      TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( tSeries, status,
      TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( tSeries->data, status,
      TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( tSeries->data->data, status,
      TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );

  /* check the contents of the parameter structure */
  ASSERT( params, status,
      TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( params->window, status,
      TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( params->window->data, status,
      TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( params->window->data->length > 0, status,
      TIMEFREQFFTH_EZSEG, TIMEFREQFFTH_MSGEZSEG );
  ASSERT( params->plan, status,
      TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  if ( !  ( params->method == useUnity || params->method == useMean ||
        params->method == useMedian ) )
  {
    ABORT( status, TIMEFREQFFTH_EUAVG, TIMEFREQFFTH_MSGEUAVG );
  }

  /* check that the window length and fft storage lengths agree */
  fLength = fSeries->data->length;
  tLength = params->window->data->length;
  if ( fLength != tLength / 2 + 1 )
  {
    ABORT( status, TIMEFREQFFTH_EMISM, TIMEFREQFFTH_MSGEMISM );
  }

  /* compute the number of segs, check that the length and overlap are valid */
  numSeg = (tSeries->data->length - params->overlap) / (tLength - params->overlap);
  if ( (tSeries->data->length - params->overlap) % (tLength - params->overlap) )
  {
    ABORT( status, TIMEFREQFFTH_EMISM, TIMEFREQFFTH_MSGEMISM );
  }

  /* clear the output spectrum and the workspace frequency series */
  memset( fSeries->data->data, 0, fLength * sizeof(REAL4) );

  /* compute the parameters of the output frequency series data */
  fSeries->epoch = tSeries->epoch;
  fSeries->f0 = tSeries->f0;
  fSeries->deltaF = 1.0 / ( (REAL8) tLength * tSeries->deltaT );
  pair.unitOne = &(tSeries->sampleUnits);
  pair.unitTwo = &(tSeries->sampleUnits);
  LALUnitMultiply( status->statusPtr, &unit, &pair );
  CHECKSTATUSPTR( status );
  pair.unitOne = &unit;
  pair.unitTwo = &lalSecondUnit;
  LALUnitMultiply( status->statusPtr, &(fSeries->sampleUnits), &pair );
  CHECKSTATUSPTR( status );

  /* if this is a unit spectrum, just set the conents to unity and return */
  if ( params->method == useUnity )
  {
    for ( k = 0; k < fLength; ++k )
    {
      fSeries->data->data[k] = 1.0;
    }
    DETATCHSTATUSPTR( status );
    RETURN( status );
  }

  /* create temporary storage for the dummy time domain segment */
  LALCreateVector( status->statusPtr, &tSegment, tLength );
  CHECKSTATUSPTR( status );

  /* create temporary storage for the individual ffts */
  LALCCreateVector( status->statusPtr, &fSegment, fLength );
  CHECKSTATUSPTR( status );

  if ( params->method == useMedian )
  {
    /* create enough storage for the indivdiual power spectra */
    psdSeg = XLALCalloc( numSeg, fLength * sizeof(REAL4) );
  }

  /* compute each of the power spectra used in the average */
  for ( i = 0, tSeriesPtr = tSeries->data->data; i < (UINT4) numSeg; ++i )
  {
    /* copy the time series data to the dummy segment */
    memcpy( tSegment->data, tSeriesPtr, tLength * sizeof(REAL4) );

    /* window the time series segment */
    for ( j = 0; j < tLength; ++j )
    {
      tSegment->data[j] *= params->window->data->data[j];
    }

    /* compute the fft of the data segment */
    LALForwardRealFFT( status->statusPtr, fSegment, tSegment, params->plan );
    CHECKSTATUSPTR (status);

    /* advance the segment data pointer to the start of the next segment */
    tSeriesPtr += tLength - params->overlap;

    /* compute the psd components */
    if ( params->method == useMean )
    {
      /* we can get away with less storage */
      for ( k = 0; k < fLength; ++k )
      {
        fftRe = crealf(fSegment->data[k]);
        fftIm = cimagf(fSegment->data[k]);
        fSeries->data->data[k] += fftRe * fftRe + fftIm * fftIm;
      }

      /* halve the DC and Nyquist components to be consistent with T010095 */
      fSeries->data->data[0] /= 2;
      fSeries->data->data[fLength - 1] /= 2;
    }
    else if ( params->method == useMedian )
    {
      /* we must store all the spectra */
      for ( k = 0; k < fLength; ++k )
      {
        fftRe = crealf(fSegment->data[k]);
        fftIm = cimagf(fSegment->data[k]);
        psdSeg[i * fLength + k] = fftRe * fftRe + fftIm * fftIm;
      }

      /* halve the DC and Nyquist components to be consistent with T010095 */
      psdSeg[i * fLength] /= 2;
      psdSeg[i * fLength + fLength - 1] /= 2;
    }
  }

  /* destroy the dummy time series segment and the fft scratch space */
  LALDestroyVector( status->statusPtr, &tSegment );
  CHECKSTATUSPTR( status );
  LALCDestroyVector( status->statusPtr, &fSegment );
  CHECKSTATUSPTR( status );

  /* compute the desired average of the spectra */
  if ( params->method == useMean )
  {
    /* normalization constant for the arithmentic mean */
    psdNorm = ( 2.0 * tSeries->deltaT ) /
      ( (REAL4) numSeg * params->window->sumofsquares );

    /* normalize the psd to it matches the conventions document */
    for ( k = 0; k < fLength; ++k )
    {
      fSeries->data->data[k] *= psdNorm;
    }
  }
  else if ( params->method == useMedian )
  {
    REAL8 bias;

    /* determine the running median bias */
    /* note: this is not the correct bias if the segments are overlapped */
    if ( params->overlap )
    {
      LALWarning( status, "Overlapping segments with median method causes a biased spectrum." );
    }
    TRY( LALRngMedBias( status->statusPtr, &bias, numSeg ), status );

    /* normalization constant for the median */
    psdNorm = ( 2.0 * tSeries->deltaT ) /
      ( bias * params->window->sumofsquares );

    /* allocate memory array for insert sort */
    s = XLALMalloc( numSeg * sizeof(REAL4) );
    if ( ! s )
    {
      ABORT( status, TIMEFREQFFTH_EMALLOC, TIMEFREQFFTH_MSGEMALLOC );
    }

    /* compute the median spectra and normalize to the conventions doc */
    for ( k = 0; k < fLength; ++k )
    {
      fSeries->data->data[k] = psdNorm *
        MedianSpec( psdSeg, s, k, fLength, numSeg );
    }

    /* free memory used for sort array */
    XLALFree( s );

    /* destroy the storage for the individual spectra */
    XLALFree( psdSeg );
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


void
LALCOMPLEX8AverageSpectrum (
    LALStatus                   *status,
    COMPLEX8FrequencySeries     *fSeries,
    REAL4TimeSeries             *tSeries0,
    REAL4TimeSeries             *tSeries1,
    AverageSpectrumParams       *params
    )

{
  UINT4                 i, j, k, l;          /* seg, ts and freq counters       */
  UINT4                 numSeg;           /* number of segments in average   */
  UINT4                 fLength;          /* length of requested power spec  */
  UINT4                 tLength;          /* length of time series segments  */
  REAL4Vector          *tSegment[2] = {NULL,NULL};
  COMPLEX8Vector       *fSegment[2] = {NULL,NULL};
  REAL4                *tSeriesPtr0,*tSeriesPtr1 ;
  REAL4                 psdNorm = 0;      /* factor to multiply windows data */
  REAL4                 fftRe0, fftIm0, fftRe1, fftIm1;
  LALUnit               unit;
  LALUnitPair           pair;
  /* RAT4                  negRootTwo = { -1, 1 }; */

  INITSTATUS(status);
  XLAL_PRINT_DEPRECATION_WARNING("XLALREAL8AverageSpectrumWelch");
  ATTATCHSTATUSPTR (status);

  /* check the input and output data pointers are non-null */
  ASSERT( fSeries, status,
      TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( fSeries->data, status,
      TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( fSeries->data->data, status,
      TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( tSeries0, status,
      TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( tSeries0->data, status,
      TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( tSeries0->data->data, status,
      TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
   ASSERT( tSeries1, status,
      TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( tSeries1->data, status,
      TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( tSeries1->data->data, status,
      TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );

  /* check the contents of the parameter structure */
  ASSERT( params, status,
      TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( params->window, status,
      TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( params->window->data, status,
      TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( params->window->data->length > 0, status,
      TIMEFREQFFTH_EZSEG, TIMEFREQFFTH_MSGEZSEG );
  ASSERT( params->plan, status,
      TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  if ( !  ( params->method == useUnity || params->method == useMean ||
        params->method == useMedian ) )
  {
    ABORT( status, TIMEFREQFFTH_EUAVG, TIMEFREQFFTH_MSGEUAVG );
  }

  /* check that the window length and fft storage lengths agree */
  fLength = fSeries->data->length;
  tLength = params->window->data->length;
  if ( fLength != tLength / 2 + 1 )
  {
    ABORT( status, TIMEFREQFFTH_EMISM, TIMEFREQFFTH_MSGEMISM );
  }

  /* compute the number of segs, check that the length and overlap are valid */
  numSeg = (tSeries0->data->length - params->overlap) / (tLength - params->overlap);
  if ( (tSeries0->data->length - params->overlap) % (tLength - params->overlap) )
  {
    ABORT( status, TIMEFREQFFTH_EMISM, TIMEFREQFFTH_MSGEMISM );
  }

  /* clear the output spectrum and the workspace frequency series */
  memset( fSeries->data->data, 0, fLength * sizeof(COMPLEX8) );

  /* compute the parameters of the output frequency series data */
  fSeries->epoch = tSeries0->epoch;
  fSeries->f0 = tSeries0->f0;
  fSeries->deltaF = 1.0 / ( (REAL8) tLength * tSeries0->deltaT );
  pair.unitOne = &(tSeries0->sampleUnits);
  pair.unitTwo = &(tSeries0->sampleUnits);
  LALUnitMultiply( status->statusPtr, &unit, &pair );
  CHECKSTATUSPTR( status );
  pair.unitOne = &unit;
  pair.unitTwo = &lalSecondUnit;
  LALUnitMultiply( status->statusPtr, &(fSeries->sampleUnits), &pair );
  CHECKSTATUSPTR( status );

  /* create temporary storage for the dummy time domain segment */
  for (l = 0; l < 2; l ++){
  LALCreateVector( status->statusPtr, &tSegment[l], tLength );
  CHECKSTATUSPTR( status );

  /* create temporary storage for the individual ffts */
  LALCCreateVector( status->statusPtr, &fSegment[l], fLength );
  CHECKSTATUSPTR( status );}


  /* compute each of the power spectra used in the average */
  for ( i = 0, tSeriesPtr0 = tSeries0->data->data, tSeriesPtr1 = tSeries1->data->data; i < (UINT4) numSeg; ++i )
  {
    /* copy the time series data to the dummy segment */
    memcpy( tSegment[0]->data, tSeriesPtr0, tLength * sizeof(REAL4) );
    memcpy( tSegment[1]->data, tSeriesPtr1, tLength * sizeof(REAL4) );
    /* window the time series segment */
    for ( j = 0; j < tLength; ++j )
    {
      tSegment[0]->data[j] *= params->window->data->data[j];
      tSegment[1]->data[j] *= params->window->data->data[j];
    }

    /* compute the fft of the data segment */
    LALForwardRealFFT( status->statusPtr, fSegment[0], tSegment[0], params->plan );
    CHECKSTATUSPTR (status);
    LALForwardRealFFT( status->statusPtr, fSegment[1], tSegment[1], params->plan );
    CHECKSTATUSPTR (status);

    /* advance the segment data pointer to the start of the next segment */
    tSeriesPtr0 += tLength - params->overlap;
    tSeriesPtr1 += tLength - params->overlap;

    /* compute the psd components */
    /*use mean method here*/
      /* we can get away with less storage */
      for ( k = 0; k < fLength; ++k )
      {
        fftRe0 = crealf(fSegment[0]->data[k]);
        fftIm0 = cimagf(fSegment[0]->data[k]);
        fftRe1 = crealf(fSegment[1]->data[k]);
        fftIm1 = cimagf(fSegment[1]->data[k]);
        fSeries->data->data[k] += fftRe0 * fftRe1 + fftIm0 * fftIm1;
        fSeries->data->data[k] += I * (- fftIm0 * fftRe1 + fftRe0 * fftIm1);

      }

      /* halve the DC and Nyquist components to be consistent with T010095 */
      fSeries->data->data[0] /= 2;
      fSeries->data->data[fLength - 1] /= 2;

      }

  /* destroy the dummy time series segment and the fft scratch space */
  for (l = 0; l < 2; l ++){
  LALDestroyVector( status->statusPtr, &tSegment[l] );
  CHECKSTATUSPTR( status );
  LALCDestroyVector( status->statusPtr, &fSegment[l] );
  CHECKSTATUSPTR( status );}

  /* compute the desired average of the spectra */

    /* normalization constant for the arithmentic mean */
    psdNorm = ( 2.0 * tSeries1->deltaT ) /
      ( (REAL4) numSeg * params->window->sumofsquares );

    /* normalize the psd to it matches the conventions document */
    for ( k = 0; k < fLength; ++k )
      fSeries->data->data[k] *= psdNorm;


  DETATCHSTATUSPTR( status );
  RETURN( status );
}




void
LALTimeFreqComplexFFT(
    LALStatus               *status,
    COMPLEX8FrequencySeries *freq,
    COMPLEX8TimeSeries      *time,
    ComplexFFTPlan          *plan
    )
{
  INITSTATUS(status);
  XLAL_PRINT_DEPRECATION_WARNING("XLALCOMPLEX8TimeFreqFFT");

  ASSERT( plan, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( freq, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( time, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( time->data, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( time->data->length, status, TIMEFREQFFTH_ESIZE, TIMEFREQFFTH_MSGESIZE );
  ASSERT( time->deltaT > 0, status, TIMEFREQFFTH_ERATE, TIMEFREQFFTH_MSGERATE );
  ASSERT( plan->sign == -1, status, TIMEFREQFFTH_ESIGN, TIMEFREQFFTH_MSGESIGN );

  if ( XLALCOMPLEX8TimeFreqFFT( freq, time, plan ) == XLAL_FAILURE )
  {
    XLALClearErrno();
    ABORTXLAL( status );
  }

  RETURN( status );
}



void
LALFreqTimeComplexFFT(
    LALStatus               *status,
    COMPLEX8TimeSeries      *time,
    COMPLEX8FrequencySeries *freq,
    ComplexFFTPlan          *plan
    )
{
  INITSTATUS(status);
  XLAL_PRINT_DEPRECATION_WARNING("XLALCOMPLEX8FreqTimeFFT");

  ASSERT( plan, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( time, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( freq, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( freq->data, status, TIMEFREQFFTH_ENULL, TIMEFREQFFTH_MSGENULL );
  ASSERT( freq->data->length, status, TIMEFREQFFTH_ESIZE, TIMEFREQFFTH_MSGESIZE );
  ASSERT( freq->deltaF > 0, status, TIMEFREQFFTH_ERATE, TIMEFREQFFTH_MSGERATE );
  ASSERT( plan->sign == 1, status, TIMEFREQFFTH_ESIGN, TIMEFREQFFTH_MSGESIGN );

  if ( XLALCOMPLEX8FreqTimeFFT( time, freq, plan ) == XLAL_FAILURE )
  {
    XLALClearErrno();
    ABORTXLAL( status );
  }

  RETURN( status );
}
