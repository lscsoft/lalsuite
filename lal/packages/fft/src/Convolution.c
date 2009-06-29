/*
*  Copyright (C) 2007 Bernd Machenschalk, Kipp Cannon, Patrick Brady, Saikat Ray-Majumder
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
#include <string.h>
#include <lal/AVFactories.h>
#include <lal/Calibration.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Units.h>
#include <lal/VectorOps.h>
#include <lal/Window.h>
#include <lal/FrequencySeries.h>

#include <lal/LALRCSID.h>
NRCSID (CONVOLUTIONC,"$Id$");

/**
 * \addtogroup fft
 * @{
 */

/**
 *
 * \brief Apply transfer function to time series
 *
 *
 * This function returns the convolution of the time series with
 * the frequency domain transfer function that has been supplied by
 * the user.  It zero pads the input data by a factor of two to
 * alleviate wraparound from the FFT.   This means that the transfer
 * function must have
 * \f[
 * \texttt{deltaF} =  1.0 / ( 2.0 \times \texttt{strain->data->length}
 * \times \texttt{strain->data->deltaT} )
 * \f]
 *
 */
REAL4TimeSeries *XLALRespFilt(
    REAL4TimeSeries             *strain,
    COMPLEX8FrequencySeries     *transfer
    )
{
  static const char *func = "XLALRespFilt";
  REAL4Vector *tmpWave=NULL;
  COMPLEX8Vector *tmpFFTWave=NULL;
  COMPLEX8FrequencySeries *tmpTransfer=NULL;
  REAL4FFTPlan *fwdPlan=NULL, *invPlan=NULL;
  UINT4 k;
  UINT4 inTimeLength = 0;
  UINT4 paddedTimeLength = 0;
  const LALUnit countPerStrain = {0,{0,0,0,0,0,-1,1},{0,0,0,0,0,0,0}};
  const CHAR *chname = "Temporary Transfer";

  if ( ! strain || ! transfer )
    XLAL_ERROR_NULL( func, XLAL_EFAULT );
  if ( ! strain->data || ! transfer->data )
    XLAL_ERROR_NULL( func, XLAL_EINVAL );
  if ( strain->deltaT <= 0.0 || transfer->deltaF <=0 )
    XLAL_ERROR_NULL( func, XLAL_EINVAL );

  inTimeLength = strain->data->length;
  paddedTimeLength = 2 * inTimeLength;

  tmpWave = XLALCreateREAL4Vector( paddedTimeLength );
  if ( ! tmpWave )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );

  tmpFFTWave = XLALCreateCOMPLEX8Vector( inTimeLength + 1 );
  if ( ! tmpFFTWave )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );

  fwdPlan = XLALCreateForwardREAL4FFTPlan( paddedTimeLength, 0 );
  if ( ! fwdPlan )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );

  invPlan = XLALCreateReverseREAL4FFTPlan( paddedTimeLength, 0 );
  if ( ! invPlan )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );

  /* copy the signal into zero padded data */
  memset( tmpWave->data, 0, paddedTimeLength * sizeof(REAL4) );
  for( k=0; k < inTimeLength; k++ )
  {
    tmpWave->data[k] = strain->data->data[k];
  }

  /* Fourier transform the signal and multiply by transfer fn */
  XLALREAL4ForwardFFT( tmpFFTWave, tmpWave, fwdPlan );

  /* make sure the transfer function has the right units and df */
  tmpTransfer = XLALCreateCOMPLEX8FrequencySeries(chname,
      &strain->epoch, 0.0, 1.0/(paddedTimeLength * strain->deltaT),
      &countPerStrain, tmpFFTWave->length);
  XLALResponseConvert( tmpTransfer, transfer );
  XLALCCVectorMultiply( tmpFFTWave, tmpFFTWave, tmpTransfer->data );
  XLALUnitMultiply( &strain->sampleUnits, &tmpTransfer->sampleUnits, &strain->sampleUnits);
  XLALDestroyCOMPLEX8FrequencySeries(tmpTransfer);

  /* Now make sure the DC term has zero real and imaginary parts and
   * that the nyquist term is real - (the nyquist term is the last
   * term).  This may in generaly not quite be true due to numerical
   * rounding and accuracies so I will manually set these to the
   * appropriate values */
  tmpFFTWave->data[0].im = 0.0;
  tmpFFTWave->data[0].re = 0.0;
  tmpFFTWave->data[inTimeLength].im = 0.0;

  memset( tmpWave->data, 0, 2*inTimeLength*sizeof(REAL4) );
  XLALREAL4ReverseFFT( tmpWave, tmpFFTWave, invPlan );

  for( k=0; k < inTimeLength ; k++)
  {
    strain->data->data[k] = tmpWave->data[k] / (2 * inTimeLength);
  }

  /* Make sure everything I allocated is deallocated below */
  XLALDestroyREAL4Vector( tmpWave );
  XLALDestroyCOMPLEX8Vector( tmpFFTWave );
  XLALDestroyREAL4FFTPlan( fwdPlan );
  XLALDestroyREAL4FFTPlan( invPlan );

  return  strain;
}

/**
 *
 * \brief SHOULD Convolve two time series, but doesn't
 *
 * This function does nothing yet
 *
 */
REAL4TimeSeries *XLALREAL4Convolution(
    REAL4TimeSeries             *strain,
    REAL4TimeSeries             *transfer
    )
{
  static const char *func = "XLALRespFilt";
  /*
  REAL4Vector *tmpWave;
  REAL4 normfac;
  UINT4 k;
  */

  if ( ! strain || ! transfer )
      XLAL_ERROR_NULL( func, XLAL_EFAULT );
  if ( ! strain->data || ! transfer->data )
      XLAL_ERROR_NULL( func, XLAL_EINVAL );
  if ( strain->deltaT <= 0.0 )
      XLAL_ERROR_NULL( func, XLAL_EINVAL );
  if ( transfer->data->length != strain->data->length )
      XLAL_ERROR_NULL( func, XLAL_EBADLEN );

  return strain;
}

/** @} */


