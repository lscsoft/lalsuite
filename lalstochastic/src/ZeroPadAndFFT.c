/*
*  Copyright (C) 2007 Jolien Creighton, Robert Adam Mercer, John Whelan
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

/**
 * \author UTB Relativity Group; contact whelan@phys.utb.edu
 * \addtogroup ZeroPadAndFFT_c
 *
 * \brief Routines for zero-padding and Fourier transforming a time series.
 *
 * As described in \ref StochasticCrossCorrelation_c, data
 * streams to be cross-correlated need to be zero-padded to the same
 * length as the optimal filter via
 *
 * \f{equation}{
 * \bar{h}[k]=\
 * \left\{ \begin{array}{cl}
 * w[k] h[k]  &    k = 0, \ldots, N-1 \\
 * 0     &    k = N \ldots, M-1
 * \end{array}
 * \right.
 * \f}
 *
 * (where \f$w[k]\f$ is a windowing function)
 * before being Fourier transformed via
 * \f{equation}{
 * \widetilde{h}[\ell] := \sum_{\ell=0}^{M-1}
 * \delta t\,h[k]\,e^{-i2\pi k\ell/M}
 * \ .
 * \f}
 *
 * <tt>LALSZeroPadAndFFT()</tt> performs this operaton on a
 * \c REAL4TimeSeries of length \f$N\f$, zero-padding it to length
 * \f$M\f$ and Fourier-transforming it into a
 * \c COMPLEX8FrequencySeries of length \f$[M/2]+1\f$.
 *
 * <tt>LALCZeroPadAndFFT()</tt> performs this operaton on a
 * \c COMPLEX8TimeSeries of length \f$N\f$, zero-padding it to length
 * \f$M\f$ and Fourier-transforming it into a
 * \c COMPLEX8FrequencySeries of length \f$M\f$.
 *
 * ### Algorithm ###
 *
 * <tt>LALSZeroPadAndFFT()</tt> constructs the sequence \f$\bar{h}[k]\f$, and
 * then applies a real-to-complex time-to-frequency discrete Fourier
 * transform from the \c fft package.
 *
 * <tt>LALCZeroPadAndFFT()</tt> constructs the sequence \f$\bar{h}[k]\f$, and
 * then applies a complex-to-complex time-to-frequency discrete Fourier
 * transform from the \c fft package.
 *
 * ### Uses ###
 *
 * <tt>LALSZeroPadAndFFT()\/</tt> calls:
 *
 * \code
 * LALSCreateVector()
 * LALSDestroyVector()
 * LALTimeFreqRealFFT()
 * memset()
 * strncpy()
 * \endcode
 *
 * <tt>LALSZeroPadAndFFT()\/</tt> calls:
 *
 * \code
 * LALCCreateVector()
 * LALCDestroyVector()
 * LALTimeFreqComplexFFT()
 * memset()
 * strncpy()
 * \endcode
 *
 * ### Notes ###
 *
 * <ul>
 *
 * <li> The Fourier transform is defined to be the discrete
 * approximation of a continuous Fourier transorm, which makes it \f$\delta
 * t\f$ times the discrete Fourier transform.</li>
 *
 * <li> The Fourier transform of a series of \f$M\f$ points is calculated
 * with the FFTW [\ref fj_1998] (via the interfaces in
 * the \c fft package), which is efficient for products of small
 * primes, so \f$M\f$ should be chosen to have this property.  The minimum
 * value, \f$2N-1\f$, is odd and can thus be at best a power of 3.
 * Additionally, if \f$2N-1\f$ is a convenient number, \f$N\f$ will likely not
 * be, which is one reason it might be convenient to work with \f$M=2N\f$
 * instead.</li>
 *
 * <li> <tt>LALCZeroPadAndFFT()</tt> inherits its behavior from
 * <tt>LALTimeFreqComplexFFT()</tt>, which currently does not use the
 * initial phase of the reference oscillator.  The calling routine must
 * therefore remove the effects of this phase explicitly in order to
 * obtain the band-limited FFT of the unheterodyned data.</li>
 *
 * <li> The output units are determined from the input units, but under
 * normal circumstances in the context of a stochastic background
 * search, we will have
 * \f{eqnarray}{
 * {} [h(t)] &=& \textrm{count}\\
 * {} [\widetilde{\bar{h}}(f)] &:=& [h(t)] \,\textrm{Hz}^{-1}
 * = \textrm{count}\,\textrm{Hz}^{-1}
 * \f}
 *
 * </li>
 * </ul>
 *
 * @{
 */

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/LALStdlib.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/AVFactories.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Units.h>
#include <string.h>

void
LALSZeroPadAndFFT(
    LALStatus                *status,
    COMPLEX8FrequencySeries  *output,
    const REAL4TimeSeries    *input,
    SZeroPadAndFFTParameters *parameters)

{
  UINT4 length, fullLength;
  REAL4TimeSeries  hBar;
  REAL4 *sPtr, *sStopPtr, *hBarPtr, *windowPtr;

  /* initialize status structure */
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  /* ERROR CHECKING --------------------------------------------------- */

  /* check that pointer to real timer series for input is non-null */
  ASSERT(input != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to data member of real time series for input is
   * non-null */
  ASSERT(input->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that length of data member of real time series for input is
   * not equal to zero */
  length = input->data->length;
  ASSERT(length != 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* check that pointer to data-data member of real time series for input
   * is non-null */
  ASSERT(input->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to parameter structure is non-null */
  ASSERT(parameters != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to FFT plan parameter is non-null */
  ASSERT(parameters->fftPlan != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* These checks are only relevant if a window function is specified */
  if (parameters->window != NULL)
  {
    /* check that pointer to data member of window function is non-null */
    ASSERT(parameters->window->data->data != NULL, status, \
        STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
        STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* check that window function is same length as input time series */
    if (parameters->window->data->length != length)
    {
      ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
          STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
    }
  }

  /* check that zero-padded output is not shorter than input */
  fullLength = parameters->length;
  if (fullLength < length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* check that pointer to complex frequency series for output is non-null */
  ASSERT(output != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to data member of complex frequency series for
   * output is non-null */
  ASSERT(output->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that length of complex frequency series for output
     is consistent with length of zero-padded input */
  if (fullLength/2 + 1 != output->data->length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* check that pointer to data-data member of complex frequency series
   * for output is non-null */
  ASSERT(output->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that heterodyning frequency of real frequency series for
   * the input is equal to zero */
  if (input->f0 != 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENONZEROHETERO, \
           STOCHASTICCROSSCORRELATIONH_MSGENONZEROHETERO);
  }

  /* check that frequency spacing is positive */
#ifndef LAL_NDEBUG
  REAL8 deltaT;
  deltaT = input->deltaT;
  ASSERT(deltaT > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAT, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAT);
#endif

  /* EVERYTHING OKAY HERE! -------------------------------------------- */

  /* replicate input for zero-padding */
  hBar = *input;
  hBar.data = NULL;

  /* allocate memory for zero-padded vector */
  TRY(LALSCreateVector(status->statusPtr, &(hBar.data), fullLength), status);

  if (parameters->window == NULL)
  {
    sPtr = memcpy(hBar.data->data, input->data->data, length * sizeof(REAL4));
    if (sPtr != hBar.data->data)
    {
      TRY(LALSDestroyVector(status->statusPtr, &(hBar.data)), status);
      ABORT(status, STOCHASTICCROSSCORRELATIONH_EMEMORY, \
          STOCHASTICCROSSCORRELATIONH_MSGEMEMORY );
    }
  }
  else
  {
    /* window data */
    sStopPtr = input->data->data + input->data->length;
    for (sPtr = input->data->data, hBarPtr = hBar.data->data, \
        windowPtr = parameters->window->data->data ; sPtr < sStopPtr ; \
        ++sPtr, ++hBarPtr, ++windowPtr)
    {
      *(hBarPtr) = *(sPtr) * *(windowPtr);
    }
  }

  /* zero pad */
  sPtr = memset(hBar.data->data + length, 0, \
      (fullLength - length) * sizeof(REAL4));

  if (sPtr != hBar.data->data + length)
  {
    TRY(LALSDestroyVector(status->statusPtr, &(hBar.data)), status);
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMEMORY, \
        STOCHASTICCROSSCORRELATIONH_MSGEMEMORY);
  }

  /* take DFT */
  LALTimeFreqRealFFT(status->statusPtr, output, &hBar, parameters->fftPlan);

  /* Can't use TRY because we have memory allocated */
  BEGINFAIL(status)
    TRY(LALSDestroyVector(status->statusPtr, &(hBar.data)), status);
  ENDFAIL(status);

  /* fill output parameters */
  strncpy(output->name, "Fourier Transform of Zero-Padded Time Series", \
      LALNameLength);

  /* clean up */
  TRY(LALSDestroyVector(status->statusPtr, &(hBar.data)), status);

  /* normal exit*/
  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* SZeroPadAndFFT() */


void
LALCZeroPadAndFFT(
    LALStatus                *status,
    COMPLEX8FrequencySeries  *output,
    const COMPLEX8TimeSeries *input,
    CZeroPadAndFFTParameters *parameters)

{
  UINT4 length, fullLength;
  COMPLEX8TimeSeries  hBar;
  COMPLEX8 *cPtr, *cStopPtr, *hBarPtr;
  REAL4 *windowPtr;

  /* initialize status structure */
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  /* ERROR CHECKING --------------------------------------------------- */

  /* check that pointer to real time series for input is non-null */
  ASSERT(input != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to data member of real time series for input is
   * non-null */
  ASSERT(input->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that length of data member of real time series for input is
   * not equal to zero */
  length = input->data->length;
  ASSERT(length != 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* check that pointer to data-data member of real time series for input
   * is non-null */
  ASSERT(input->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to parameter structure is non-null */
  ASSERT(parameters != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to FFT plan parameter is non-null */
  ASSERT(parameters->fftPlan != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* These checks are only relevant if a window function is specified */
  if (parameters->window != NULL)
  {
    /* check that pointer to data member of window function is non-null */
    ASSERT(parameters->window->data != NULL, status, \
        STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
        STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* check that window function is same length as input time series */
    if (parameters->window->data->length != length)
    {
      ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
          STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
    }
  }

  /* check that zero-padded output is not shorter than input */
  fullLength = parameters->length;
  if (fullLength < length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* check that pointer to complex frequency series for output is non-null */
  ASSERT(output != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to data member of complex frequency series for
   * output is non-null */
  ASSERT(output->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that length of complex frequency series for output
   * is consistent with length of zero-padded input */
  if (fullLength != output->data->length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* check that pointer to data-data member of complex frequency series
   * for output is non-null */
  ASSERT(output->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that frequency spacing is positive */
  ASSERT(input->deltaT > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAT, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAT);

  /* EVERYTHING OKAY HERE! -------------------------------------------- */

  /* replicate input for zero-padding */
  hBar = *input;
  hBar.data = NULL;

  /* allocate memory for zero-padded vector */
  TRY(LALCCreateVector(status->statusPtr, &(hBar.data), fullLength), status);

  if (parameters->window == NULL)
  {
    cPtr = memcpy(hBar.data->data, input->data->data, \
        length * sizeof(COMPLEX8));
    if (cPtr != hBar.data->data)
    {
      TRY(LALCDestroyVector(status->statusPtr, &(hBar.data)), status);
      ABORT(status, STOCHASTICCROSSCORRELATIONH_EMEMORY, \
          STOCHASTICCROSSCORRELATIONH_MSGEMEMORY );
    }
  }
  else
  {
    /* window data */
    cStopPtr = input->data->data + input->data->length;
    for (cPtr = input->data->data, hBarPtr = hBar.data->data, \
        windowPtr = parameters->window->data->data ; cPtr < cStopPtr ; \
        ++cPtr, ++hBarPtr, ++windowPtr )
    {
      hBarPtr->realf_FIXME = crealf(*cPtr) * *(windowPtr);
      hBarPtr->imagf_FIXME = cimagf(*cPtr) * *(windowPtr);
    }
  }

  /* zero pad */
  cPtr = memset(hBar.data->data + length, 0, \
      (fullLength - length) * sizeof(COMPLEX8));

  if (cPtr != hBar.data->data + length)
  {
    TRY(LALCDestroyVector(status->statusPtr, &(hBar.data)), status);
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMEMORY, \
        STOCHASTICCROSSCORRELATIONH_MSGEMEMORY );
  }

  /* take DFT */
  LALTimeFreqComplexFFT(status->statusPtr, output, &hBar, parameters->fftPlan);

  /* Can't use TRY because we have memory allocated */
  BEGINFAIL(status)
    TRY(LALCDestroyVector(status->statusPtr, &(hBar.data)), status);
  ENDFAIL(status);

  /* fill output parameters */
  strncpy(output->name, "Fourier Transform of Zero-Padded Time Series", \
      LALNameLength);

  /* clean up */
  TRY(LALCDestroyVector(status->statusPtr, &(hBar.data)), status);

  /* normal exit*/
  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* CZeroPadAndFFT() */

/** @} */
