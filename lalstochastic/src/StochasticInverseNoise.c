/*
*  Copyright (C) 2007 Jolien Creighton, Robert Adam Mercer, Tania Regimbau, John Whelan
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
 * \addtogroup StochasticInverseNoise_c
 *
 * \brief Calculates the values of the calibrated and half-calibrated inverse
 * noise power spectra from the uncalibrated noise power spectrum and the
 * frequency-domain instrument response function.
 *
 * As described in \ref StochasticOptimalFilter_c,
 * the most convenient combinations of the noise \f$P(f)\f$ (defined by \f$\langle h(f)h(f')^*\rangle=\delta(f-f')P(f)\f$) and
 * instrument response
 * \f$\widetilde{R}(f)=h(f)/h(f)\f$ to use in
 * constructing an optimal filter are the inverse half-calibrated power
 * spectral density
 * \f{equation}{
 * \label{stochastic_e_halfCalibratedPSD}
 * \frac{1}{P^{\mathrm{HC}}(f)}=\frac{1}{\widetilde{R}(f)\,P^{\mathrm{C}}(f)}
 * =\frac{\widetilde{R}(f)^*}{P(f)}
 * \f}
 * and the inverse calibrated PSD
 * \f{equation}{
 * \label{stochastic_e_calibratedPSD}
 * \frac{1}{P^{\mathrm{C}}(f)}
 * =\frac{|\widetilde{R}(f)|^2}{P(f)}
 * \f}
 * The function <tt>LALStochasticInverseNoise()</tt> takes in a
 * \c REAL4FrequencySeries describing the uncalibrated PSD
 * \f$P(f)\f$ along with a
 * \c COMPLEX8FrequencySeries describing the frequency-domain
 * response \f$\widetilde{R}(f)\f$, and outputs a
 * \c REAL4FrequencySeries describing the calibrated inverse PSD
 * \f$1/P^{\mathrm{C}}(f)\f$
 * along with a \c COMPLEX8FrequencySeries describing the
 * half-calibrated inverse PSD \f$1/P^{\mathrm{HC}}(f)\f$.
 *
 * ### Algorithm ###
 *
 * The output series are filled according to a straightforward
 * implemementation of
 * \eqref{stochastic_e_halfCalibratedPSD}--\eqref{stochastic_e_calibratedPSD}.
 * The DC components, if included in the series, are set to zero.
 *
 * ### Uses ###
 *
 * \code
 * XLALUnitRaiseRAT4()
 * XLALUnitMultiply()
 * strncpy()
 * \endcode
 *
 * ### Notes ###
 *
 * <ul>
 * <li> Note that although \f$P^{\mathrm{C}}(f)\f$
 * and \f$P(f)\f$
 * are real, \f$P^{\mathrm{HC}}(f)\f$ is \e complex.</li>
 * <li> The output units are constructed by combining the input units,
 * but under normal circumstances the units will be as follows:
 * \f{eqnarray}{
 * {} [P] &=& \mathrm{count}^{2}\, \mathrm{Hz}^{-1}\\
 * {} [\widetilde{R}] &=& 10^{18}\,\mathrm{strain}^{-1}\,\mathrm{count} \\
 * {} [1/P^{\mathrm{C}}] &:=& [\widetilde{R}]^2 [P] = 10^{36}\,\mathrm{Hz}\,\mathrm{strain}^{-2} \\
 * {} [1/P^{\mathrm{HC}}] &:=&  [\widetilde{R}] [P] = 10^{18}\,\mathrm{Hz}\,\mathrm{strain}^{-1}\,\mathrm{count}^{-1}
 * \f}</li>
 * </ul>
 *
 * @{
 */


#include <lal/LALStdlib.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/FrequencySeries.h>
#include <string.h>

/** \cond DONT_DOXYGEN */
#define invNoise output->calibratedInverseNoisePSD
#define hwInvNoise output->halfCalibratedInverseNoisePSD
#define wNoise input->unCalibratedNoisePSD
#define wFilter input->responseFunction

/** \endcond */

void
LALStochasticInverseNoiseCal(
    LALStatus                         *status,
    StochasticInverseNoiseCalOutput   *output,
    const StochasticInverseNoiseInput *input )
{
  REAL8 deltaF;
  REAL8 f0;
  UINT4 length;

  REAL4 *sPtrPW, *sPtrIP, *sStopPtr;
  COMPLEX8 *cPtrR, *cPtrIPHC;

  RAT4 power;
  LALUnit wInvNoiseUnits;

  COMPLEX8FrequencySeries *hcInvNoise;

  /* initialize status structure */
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  /*****************************************************************
   *                                                               *
   *                    Test validity of inputs                    *
   *                                                               *
   *****************************************************************/

  /* check that pointer to input structure is not null */
  ASSERT(input != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to output structure is not null */
  ASSERT(output != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointers to members of input structure are not null */
  ASSERT(wNoise != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(wFilter != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointers to members of output structure are not null */
  ASSERT(invNoise != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);


  /* check that pointers to data members of series are not null */
  ASSERT(wNoise->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(wFilter->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(invNoise->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointers to data-data members of series are not null */
  ASSERT(wNoise->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(wFilter->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(invNoise->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that length is not zero */
  length = wNoise->data->length;
  ASSERT(length > 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* check that lengths of all series match */
  if (wFilter->data->length != length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  if (invNoise->data->length != length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* check that frequency spacing is positive */
  deltaF = wNoise->deltaF;
  ASSERT(deltaF > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  /* check that frequency spacings of input series match */
  if (wFilter->deltaF != deltaF)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }

  /* set frequency spacing of output series */
  invNoise->deltaF = deltaF;

  /* check that initial frequency is non-negative */
  f0 = wNoise->f0;
  if (f0 < 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);
  }

  /* check that initial frequency of input series match */
  if (wFilter->f0 != f0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }

  /* set initial frequency of output series */
  invNoise->f0 = f0;

  /* set epochs */
  invNoise->epoch = wNoise->epoch;

  /*---------------Valid data here---------------*/

  strncpy(invNoise->name, "Calibrated invserse noise PSD", LALNameLength);

  /* allocate memory for half calibrated inverse noise */
  hcInvNoise = XLALCreateCOMPLEX8FrequencySeries(\
        "half-calibrated invserse noise PSD", &wNoise->epoch, f0, deltaF, \
        &lalDimensionlessUnit, length);
  if(!hcInvNoise)
    ABORTXLAL(status);

  /* unit structure manipulation */

  /* Find units of uncalibrated inverse power spectrum */
  power.numerator = -1;
  power.denominatorMinusOne = 0;
  if (XLALUnitRaiseRAT4(&wInvNoiseUnits, &(wNoise->sampleUnits), &power) == NULL) {
    ABORTXLAL(status);
  }

  /* multiply by response function units to get half-calibrated inv noise
   * units */
  if (XLALUnitMultiply(&(hcInvNoise->sampleUnits), &(wFilter->sampleUnits), &wInvNoiseUnits) == NULL) {
    ABORTXLAL(status);
  }

  /* multiply by response function units to get calibrated inv noise units */
  if (XLALUnitMultiply(&(invNoise->sampleUnits), &(wFilter->sampleUnits), &(hcInvNoise->sampleUnits)) == NULL) {
    ABORTXLAL(status);
  }


  sStopPtr = wNoise->data->data + length;

  if (f0 == 0)
  {
    /* set DC channel to zero */
    hcInvNoise->data->data[0] = 0;
    invNoise->data->data[0] = 0;

    /* initialize pointers */
    sPtrPW = wNoise->data->data + 1;
    cPtrR = wFilter->data->data + 1;
    sPtrIP = invNoise->data->data + 1;
    cPtrIPHC = hcInvNoise->data->data + 1;
  } /* if (f0 == 0) */
  else
  {
    /* initialize pointers */
    sPtrPW = wNoise->data->data;
    cPtrR = wFilter->data->data;
    sPtrIP = invNoise->data->data;
    cPtrIPHC = hcInvNoise->data->data;
  }

  for (; sPtrPW < sStopPtr ; ++sPtrPW, ++cPtrR, ++sPtrIP, ++cPtrIPHC)
  {
    *sPtrIP = (crealf(*cPtrR)*crealf(*cPtrR) + cimagf(*cPtrR)*cimagf(*cPtrR)) / *sPtrPW;
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);
} /* LALStochasticInverseNoiseCal() */


void
LALStochasticInverseNoise(
    LALStatus                         *status,
    StochasticInverseNoiseOutput      *output,
    const StochasticInverseNoiseInput *input )

{
  REAL8 deltaF;
  REAL8 f0;
  UINT4 length;

  REAL4 *sPtrPW, *sPtrIP, *sStopPtr;
  COMPLEX8 *cPtrR, *cPtrIPHC;

  RAT4 power;
  LALUnit wInvNoiseUnits;

  /* initialize status structure */
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  /*****************************************************************
   *                                                               *
   *                    Test validity of inputs                    *
   *                                                               *
   *****************************************************************/

  /* check that pointer to input structure is not null */
  ASSERT(input != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to output structure is not null */
  ASSERT(output != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointers to members of input structure are not null */
  ASSERT(wNoise != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(wFilter != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointers to members of output structure are not null */
  ASSERT(invNoise != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(hwInvNoise != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointers to data members of series are not null */
  ASSERT(wNoise->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(wFilter->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(invNoise->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(hwInvNoise->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointers to data-data members of series are not null */
  ASSERT(wNoise->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(wFilter->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(invNoise->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(hwInvNoise->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that length is not zero */
  length = wNoise->data->length;
  ASSERT(length > 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* check that lengths of all series match */
  if (wFilter->data->length != length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  if (invNoise->data->length != length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  if (hwInvNoise->data->length != length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* check that frequency spacing is positive */
  deltaF = wNoise->deltaF;
  ASSERT(deltaF > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  /* check that frequency spacings of input series match */
  if (wFilter->deltaF != deltaF)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }

  /* set frequency spacing of output series */
  invNoise->deltaF = deltaF;
  hwInvNoise->deltaF = deltaF;

  /* check that initial frequency is non-negative */
  f0 = wNoise->f0;
  if (f0 < 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);
  }

  /* check that initial frequency of input series match */
  if (wFilter->f0 != f0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }

  /* set initial frequency of output series */
  invNoise->f0 = f0;
  hwInvNoise->f0 = f0;

  /* set epochs */
  invNoise->epoch = hwInvNoise->epoch = wNoise->epoch;

  /*---------------Valid data here---------------*/

  strncpy(invNoise->name, "Calibrated invserse noise PSD", LALNameLength);
  strncpy(hwInvNoise->name, "half-calibrated invserse noise PSD", \
      LALNameLength);

  /* unit structure manipulation */

  /* Find units of uncalibrated inverse power spectrum */
  power.numerator = -1;
  power.denominatorMinusOne = 0;
  if (XLALUnitRaiseRAT4(&wInvNoiseUnits, &(wNoise->sampleUnits), &power) == NULL) {
    ABORTXLAL(status);
  }

  /* multiply by response function units to get half-calibrated inv noise
   * units */
  if (XLALUnitMultiply(&(hwInvNoise->sampleUnits), &(wFilter->sampleUnits), &wInvNoiseUnits) == NULL) {
    ABORTXLAL(status);
  }

  /* multiply by response function units to get calibrated inv noise units */
  if (XLALUnitMultiply(&(invNoise->sampleUnits), &(wFilter->sampleUnits), &(hwInvNoise->sampleUnits)) == NULL) {
    ABORTXLAL(status);
  }

  sStopPtr = wNoise->data->data + length;

  if (f0 == 0)
  {
    /* set DC channel to zero */
    hwInvNoise->data->data[0] = 0;
    invNoise->data->data[0] = 0;

    /* initialize pointers */
    sPtrPW = wNoise->data->data + 1;
    cPtrR = wFilter->data->data + 1;
    sPtrIP = invNoise->data->data + 1;
    cPtrIPHC = hwInvNoise->data->data + 1;
  } /* if (f0 == 0) */
  else
  {
    /* initialize pointers */
    sPtrPW = wNoise->data->data;
    cPtrR = wFilter->data->data;
    sPtrIP = invNoise->data->data;
    cPtrIPHC = hwInvNoise->data->data;
  }

  for (; sPtrPW < sStopPtr ; ++sPtrPW, ++cPtrR, ++sPtrIP, ++cPtrIPHC)
  {
    *sPtrIP = (crealf(*cPtrR)*crealf(*cPtrR) + cimagf(*cPtrR)*cimagf(*cPtrR)) / *sPtrPW;
    *cPtrIPHC = crectf( crealf(*cPtrR) / *sPtrPW,

    /* minus sign because of complex conjugate */
                        -cimagf(*cPtrR) / *sPtrPW );
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);
} /* LALStochasticInverseNoise() */

/** @} */
