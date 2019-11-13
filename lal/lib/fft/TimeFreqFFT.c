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
