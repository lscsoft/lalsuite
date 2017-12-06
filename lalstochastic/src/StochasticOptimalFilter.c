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
 * \author UTB Relativity Group; contact john.whelan@ligo.org
 * \addtogroup StochasticOptimalFilter_c
 *
 * \brief Calculates the values of the optimal filter function for the standard cross-correlation statistic.
 *
 * As described in
 * \cite Allen1997, \cite Allen1999, \cite Finn2001,
 * the optimal filter \f$\widetilde{Q}^{\mathrm{C}}(f)\f$ which maximizes the ratio of the
 * mean \f$\mu=\langle Y\rangle\f$ to the standard deviation
 * \f$\sigma=\sqrt{\langle (Y-\mu)^2\rangle}\f$ of the cross-correlation
 * statistic \eqref{stochastic_e_ymax} is
 *
 * \f{equation}{
 * \widetilde{Q}^{\mathrm{C}}(f)=\lambda\,
 * \frac{\gamma(f)\,\Omega_{\mathrm{GW}}(f)}
 * {|f|^3\,P^{\mathrm{C}}_1(f)\,P^{\mathrm{C}}_2(f)}
 * \f}
 *
 * where \f$\lambda\f$ is a normalization constant, \f$\gamma(f)\f$ is the
 * overlap reduction function (\e cf \ref OverlapReductionFunction_c) for the two
 * detectors, \f$\Omega_{\mathrm{GW}}(f)\f$ is the stochastic
 * gravitational wave background strength (\e cf \ref OverlapReductionFunction_c), and \f$P^{\mathrm{C}}_i(f)\f$ is
 * the power spectral density (\f$\langle
 * h^{\mathrm{C}}_i(f)h^{\mathrm{C}}_i(f')^*\rangle=\delta(f-f')P^{\mathrm{C}}_i(f)\f$)
 * for the \f$i\f$th detector.
 *
 * However, in practice, the data stream coming out of the \f$i\f$th detector
 * is not the strain \f$h^{\mathrm{C}}_i(t)=h_{ab}(t,\vec{x}_i)d^{ab}\f$, but that
 * convolved with an instrumental response function \f$R_i(\tau)\f$ to
 * produce an "uncalibrated" data stream
 * \f{equation}{
 * h_i(t) = \int_0^{\infty} d\tau\, R_i(\tau)\,
 * h^{\mathrm{C}}_i(t-\tau)
 * \f}
 * which has the simpler frequency-domain representation
 * \f{equation}{
 * \widetilde{h}_i(f)
 * =  \widetilde{R}_i(f)\, \widetilde{h}_i(f)
 * \f}
 * If we want to calculate the cross-correlation statistic \f$Y\f$ using the
 * uncalibrated detector output, the expression is
 * \f{equation}{
 * Y
 * = \int_{-\infty}^{\infty} df\,
 * \left(
 * \frac{\widetilde{\bar{h}}{}_{1}(f)}
 * {\widetilde{R}_{1}(f)}
 * \right)^*
 * \,
 * \widetilde{Q}^{\mathrm{C}}(f)\,
 * \left(
 * \frac{\widetilde{\bar{h}}{}_{2}(f)}
 * {\widetilde{R}_{2}(f)}
 * \right)
 * = \int_{-\infty}^{\infty} df\,
 * \widetilde{\bar{h}}{}_{1}(f) ^*
 * \,
 * \widetilde{Q}(f)\,
 * \widetilde{\bar{h}}{}_{2}(f)
 * \f}
 * where the "uncalibrated optimal filter" is
 * \f{eqnarray}{
 * \label{stochastic_e_QW}
 * \widetilde{Q}(f)
 * &=&\frac{\widetilde{Q}^{\mathrm{C}}(f)}{\widetilde{R}_1(f)^*\widetilde{R}_2(f)}
 * =\lambda\,\left(\frac{1}{\widetilde{R}_1(f)^*P^{\mathrm{C}}_1(f)}\right)
 * \frac{\gamma(f)\,\Omega_{\mathrm{GW}}(f)}
 * {|f|^3}\left(\frac{1}{\widetilde{R}_2(f)P^{\mathrm{C}}_2(f)}\right) \\
 * &=&\lambda\,
 * \frac{\gamma(f)\,\Omega_{\mathrm{GW}}(f)}
 * {|f|^3\,P^{\mathrm{HC}}_1(f)^*\,P^{\mathrm{HC}}_2(f)}
 * \ ,
 * \f}
 * where \f$P^{\mathrm{HC}}_i(f)=\widetilde{R}_i(f)\,P^{\mathrm{C}}_i(f)\f$ is the
 * "half-calibrated" PSD.  (The uncalibrated PSD is
 * \f$P_i(f)=|\widetilde{R}_i(f)|^2\,P^{\mathrm{C}}_i(f)\f$.)
 *
 * <tt>LALStochasticOptimalFilter()</tt> generates a complex frequency
 * series containing the uncalibrated optimal filter
 * \f$\widetilde{Q}(f)\f$, taking as inputs real
 * frequency series representing the overlap reduction function
 * \f$\gamma(f)\f$ and the stochastic gravitational wave background spectrum
 * \f${h_{100}}^2\Omega_{\mathrm{GW}}(f)\f$, as well as complex
 * frequency series representing the half-calibrated (inverse) PSDs
 * \f$\{1/P^{\mathrm{HC}}_i(f)|i=1,2\}\f$, and as a real parameter
 * the normalization constant \f$\lambda\f$.
 *
 * ### Algorithm ###
 *
 * The routine <tt>LALStochasticOptimalFilter()</tt> fills its output
 * series is filled with the values corresponding to the definition
 * \eqref{stochastic_e_QW}.
 *
 * ### Uses ###
 *
 * \code
 * XLALUnitMultiply()
 * XLALUnitRaiseRAT4()
 * XLALUnitCompare()
 * \endcode
 *
 * ### Notes ###
 *
 * <ul>
 * <li> If \f$f_0=0\f$, the DC element \f$Q(0)\f$ is set to zero, regardless of
 * the values of the inputs, because the \f$f^3\f$ term would make it
 * diverge otherwise, and because any conceivable realistic noise
 * spectrum will end up blowing up at zero frequency fast enough to
 * kill the optimal filter.</li>
 * <li> The implementation of the optimal filter function given here
 * assumes a large observation time continuum-limit approximation.  In
 * this limit, the Dirichlet kernels (which appear in an exact
 * expression for the standard cross-correlation statistic, when
 * evaluated in discrete time \cite Finn2001; see also
 * the documentation for the module Dirichlet.c in the utilities package)
 * may be replaced by Dirac delta functions.</li>
 * <li> Although \f$Q^{\mathrm{C}}(f)\f$ is real by construction, the uncalibrated optimal
 * filter \f$\widetilde{Q}(f)\f$ will in general be
 * complex because the response functions \f$\widetilde{R}_i(f)\f$ for the
 * two sites will be different.</li>
 * <li> The expected units for the inputs and output of this function
 * are as follows (although the actual output units will be constructed
 * from the input units):
 * \f{eqnarray}{
 * {} [\lambda] &=& 10^{-36}\,\textrm{s}^{-1}\\
 * {} [\gamma] &=& \textrm{strain}^{2} \\
 * {} [\Omega_{\mathrm{GW}}] &=& 1 \\
 * {} [1/P^{\mathrm{HC}}_{1,2}]
 * &=& 10^{18}\,\textrm{Hz}\,\textrm{strain}^{-1}\,\textrm{count}^{-1} \\
 * {} [\widetilde{Q}] &:=&
 * [\lambda] [\gamma][\Omega_{\mathrm{GW}}]
 * \left[\frac{1}{P^{\mathrm{HC}}_1}\right]
 * \left[\frac{1}{P^{\mathrm{HC}}_2}\right]
 * \,\textrm{s}^3
 * =
 * \textrm{count}^{-2}
 * \f}</li>
 * </ul>
 *
 * @{
 */

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <math.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/Units.h>

void
LALStochasticOptimalFilterCal(
    LALStatus                             *status,
    REAL4FrequencySeries                  *optimalFilter,
    const StochasticOptimalFilterCalInput *input,
    const REAL4WithUnits                  *lambda)
{
  REAL4 mygamma;
  REAL4 omega;
  REAL4 p1WInv;
  REAL4 p2WInv;

  REAL8 f;
  REAL8 f0;
  REAL8 f3;
  REAL8 deltaF;

  /* normalization factor */
  UINT4 i;
  REAL4 realFactor;

  UINT4 length;

  RAT4 power;
  LALUnit tmpUnit1, tmpUnit2, checkUnit;

  /* initialize status pointer */
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  /* ERROR CHECKING ----------------------------------------------------- */

  /* check for null pointers */
  /* input structure */
  ASSERT(input != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* output structure */
  ASSERT(optimalFilter != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* overlap member of input */
  ASSERT(input->overlapReductionFunction != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* omega member of input */
  ASSERT(input->omegaGW != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* calibrated inverse noise 1 of input */
  ASSERT(input->calibratedInverseNoisePSD1 != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* calibrated inverse noise 2 of input */
  ASSERT(input->calibratedInverseNoisePSD2 != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member of overlap */
  ASSERT(input->overlapReductionFunction->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member of omega */
  ASSERT(input->omegaGW->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member of calibrated inverse noise 1 */
  ASSERT(input->calibratedInverseNoisePSD1->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member of calibrated inverse noise 2 */
  ASSERT(input->calibratedInverseNoisePSD2->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member of output */
  ASSERT(optimalFilter->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member of overlap */
  ASSERT(input->overlapReductionFunction->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member of omega */
  ASSERT(input->omegaGW->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member of calibrated inverse noise 1 */
  ASSERT(input->calibratedInverseNoisePSD1->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member of calibrated inverse noise 2 */
  ASSERT(input->calibratedInverseNoisePSD2->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member of output structure */
  ASSERT(optimalFilter->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /*** done with null pointers ***/

  /* extract parameters from overlap */
  length = input->overlapReductionFunction->data->length;
  f0 = input->overlapReductionFunction->f0;
  deltaF = input->overlapReductionFunction->deltaF;

  /**** check for legality ****/
  /* length must be positive */
  ASSERT(length != 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* start frequency must not be negative */
  if (f0 < 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN );
  }

  /* frequency spacing must be positive */
  ASSERT(deltaF > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  /** check for mismatches **/
  /* length */
  if (input->omegaGW->data->length != length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  if (input->calibratedInverseNoisePSD1->data->length != length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  if (input->calibratedInverseNoisePSD2->data->length != length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  if (optimalFilter->data->length != length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* initial frequency */
  if (input->omegaGW->f0 != f0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }
  if (input->calibratedInverseNoisePSD1->f0 != f0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }
  if (input->calibratedInverseNoisePSD2->f0 != f0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }

  /* frequency spacing */
  if (input->omegaGW->deltaF != deltaF)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }
  if (input->calibratedInverseNoisePSD1->deltaF != deltaF)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }
  if (input->calibratedInverseNoisePSD2->deltaF != deltaF)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }

  /* EVERYHTING OKAY HERE! ---------------------------------------------- */

  /* assign parameters to optimalFilter */
  optimalFilter->f0 = f0;
  optimalFilter->deltaF = deltaF;
  optimalFilter->epoch.gpsSeconds = 0;
  optimalFilter->epoch.gpsNanoSeconds = 0;
  strncpy(optimalFilter->name, "Optimal filter for stochastic search", \
      LALNameLength);

  /* All the powers we use are integers, so we can do this once here */
  power.denominatorMinusOne = 0;

  /* Set tmpUnit1 to dims of Omega/H0^2 ******/

  /* First, set it to dims of H0 */
  tmpUnit1 = lalHertzUnit;

  /* Account for scaled units of Hubble constant */
  tmpUnit1.powerOfTen -= 18;

  /* Now set tmpUnit2 to dims of H0^-2 */
  power.numerator = -2;
  if (XLALUnitRaiseRAT4(&tmpUnit2, &tmpUnit1, &power) == NULL) {
    ABORTXLAL(status);
  }

  if (XLALUnitMultiply(&tmpUnit1, &(input->omegaGW->sampleUnits), &tmpUnit2) == NULL) {
    ABORTXLAL(status);
  }

  /* Now tmpUnit1 has units of Omega/H0^2 */

  /* Now we need to set the Optimal Filter Units equal to the units of
   * lambda*mygamma*Omega*f^-3*P1W^-1*P2W^-1) */

  if (XLALUnitMultiply(&tmpUnit1, &(input->calibratedInverseNoisePSD1->sampleUnits), &(input->calibratedInverseNoisePSD2->sampleUnits)) == NULL) {
    ABORTXLAL(status);
  }

  /* tmpUnit1 now holds the units of P1W^-1*P2W^-1 */

  power.numerator = -3;
  if (XLALUnitRaiseRAT4(&tmpUnit2, &lalHertzUnit, &power) == NULL) {
    ABORTXLAL(status);
  }

  /* tmpUnit2 now holds the units of f^-3 */

  if (XLALUnitMultiply(&checkUnit, &tmpUnit1, &tmpUnit2) == NULL) {
    ABORTXLAL(status);
  }

  /* checkUnit now holds the units of f^-3*P1HW^-1*P2HW^-1) */
  if (XLALUnitMultiply(&tmpUnit1, &checkUnit, &(input->omegaGW->sampleUnits)) == NULL) {
    ABORTXLAL(status);
  }

  /* tmpUnit1 now holds units of Omega*f^-3*P1W^-1*P2W^-1) */

  if (XLALUnitMultiply(&tmpUnit2, &tmpUnit1, &(input->overlapReductionFunction->sampleUnits)) == NULL) {
    ABORTXLAL(status);
  }

  /* tmpUnit2 now holds units of mygamma*Omega*f^-3*P1W^-1*P2W^-1) */

  if (XLALUnitMultiply(&(optimalFilter->sampleUnits), &(lambda->units), &tmpUnit2) == NULL) {
    ABORTXLAL(status);
  }

  /* Done with unit manipulation */

  optimalFilter->data->data[0] = 0;

  /* calculate optimal filter values */
  for (i = (f0 == 0 ? 1 : 0) ; i < length; ++i)
  {
    f = f0 + deltaF * (REAL8)i;

    f3 = f * f * f;

    omega = input->omegaGW->data->data[i];
    mygamma = input->overlapReductionFunction->data->data[i];
    p1WInv = input->calibratedInverseNoisePSD1->data->data[i];
    p2WInv = input->calibratedInverseNoisePSD2->data->data[i];

		realFactor = (mygamma * omega * lambda->value) / f3;

    optimalFilter->data->data[i] = realFactor * p1WInv * p2WInv;
	}

  DETATCHSTATUSPTR(status);
  RETURN(status);
} /* LALStochasticOptimalFilterCal() */



void
LALStochasticOptimalFilter(
    LALStatus                          *status,
    COMPLEX8FrequencySeries            *optimalFilter,
    const StochasticOptimalFilterInput *input,
    const REAL4WithUnits               *lambda)

{
  REAL4 mygamma;
  REAL4 omega;
  COMPLEX8 p1HWInv;
  COMPLEX8 p2HWInv;

  COMPLEX8 *cPtrOptimalFilter;

  REAL8 f;
  REAL8 f0;
  REAL8 f3;
  REAL8 deltaF;

  /* normalization factor */
  UINT4 i;
  REAL8 realFactor;

  UINT4 length;

  RAT4 power;
  LALUnit tmpUnit1, tmpUnit2, checkUnit;

  /* initialize status pointer */
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  /* ERROR CHECKING ----------------------------------------------------- */

  /***** check for null pointers *****/
  /* input structure */
  ASSERT(input != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* output structure */
  ASSERT(optimalFilter != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* overlap member of input */
  ASSERT(input->overlapReductionFunction != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* omega member of input */
  ASSERT(input->omegaGW != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* half-calibrated inverse noise 1 of input */
  ASSERT(input->halfCalibratedInverseNoisePSD1 != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* half-calibrated inverse noise 2 of input */
  ASSERT(input->halfCalibratedInverseNoisePSD2 != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member of overlap */
  ASSERT(input->overlapReductionFunction->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member of omega */
  ASSERT(input->omegaGW->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member of half-calibrated inverse noise 1 */
  ASSERT(input->halfCalibratedInverseNoisePSD1->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member of half-calibrated inverse noise 2 */
  ASSERT(input->halfCalibratedInverseNoisePSD2->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member of output */
  ASSERT(optimalFilter->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member of overlap */
  ASSERT(input->overlapReductionFunction->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member of omega */
  ASSERT(input->omegaGW->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member of half calibrated inverse noise 1 */
  ASSERT(input->halfCalibratedInverseNoisePSD1->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member of half calibrated inverse noise 2 */
  ASSERT(input->halfCalibratedInverseNoisePSD2->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member of output structure */
  ASSERT(optimalFilter->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /*** done with null pointers ***/

  /* extract parameters from overlap */
  length = input->overlapReductionFunction->data->length;
  f0 = input->overlapReductionFunction->f0;
  deltaF = input->overlapReductionFunction->deltaF;

  /**** check for legality ****/
  /* length must be positive */
  ASSERT(length != 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* start frequency must not be negative */
  if (f0 < 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN );
  }

  /* frequency spacing must be positive */
  ASSERT(deltaF > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  /** check for mismatches **/
  /* length */
  if (input->omegaGW->data->length != length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  if (input->halfCalibratedInverseNoisePSD1->data->length != length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  if (input->halfCalibratedInverseNoisePSD2->data->length != length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  if (optimalFilter->data->length != length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* initial frequency */
  if (input->omegaGW->f0 != f0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }
  if (input->halfCalibratedInverseNoisePSD1->f0 != f0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }
  if (input->halfCalibratedInverseNoisePSD2->f0 != f0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }

  /* frequency spacing */
  if (input->omegaGW->deltaF != deltaF)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }
  if (input->halfCalibratedInverseNoisePSD1->deltaF != deltaF)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }
  if (input->halfCalibratedInverseNoisePSD2->deltaF != deltaF)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }
  /* EVERYHTING OKAY HERE! ---------------------------------------------- */

  /* assign parameters to optimalFilter */
  optimalFilter->f0 = f0;
  optimalFilter->deltaF = deltaF;
  optimalFilter->epoch.gpsSeconds = 0;
  optimalFilter->epoch.gpsNanoSeconds = 0;
  strncpy(optimalFilter->name, "Optimal filter for stochastic search", \
           LALNameLength);

  /* All the powers we use are integers, so we can do this once here */
  power.denominatorMinusOne = 0;

  /* Set tmpUnit1 to dims of Omega/H0^2 ******/

  /* First, set it to dims of H0 */

  tmpUnit1 = lalHertzUnit;

  /* Account for scaled units of Hubble constant */
  tmpUnit1.powerOfTen -= 18;

  /* Now set tmpUnit2 to dims of H0^-2 */
  power.numerator = -2;
  if (XLALUnitRaiseRAT4(&tmpUnit2, &tmpUnit1, &power) == NULL) {
    ABORTXLAL(status);
  }

  if (XLALUnitMultiply(&tmpUnit1, &(input->omegaGW->sampleUnits), &tmpUnit2) == NULL) {
    ABORTXLAL(status);
  }

  /* Now tmpUnit1 has units of Omega/H0^2 */

  /* Now we need to set the Optimal Filter Units equal to the units of */
  /* lambda*mygamma*Omega*f^-3*P1HW^-1*P2HW^-1) */

  if (XLALUnitMultiply(&tmpUnit1, &(input->halfCalibratedInverseNoisePSD1->sampleUnits), &(input->halfCalibratedInverseNoisePSD2->sampleUnits)) == NULL) {
    ABORTXLAL(status);
  }

  /* tmpUnit1 now holds the units of P1HW^-1*P2HW^-1 */

  power.numerator = -3;
  if (XLALUnitRaiseRAT4(&tmpUnit2, &lalHertzUnit, &power) == NULL) {
    ABORTXLAL(status);
  }

  /* tmpUnit2 now holds the units of f^-3 */

  if (XLALUnitMultiply(&checkUnit, &tmpUnit1, &tmpUnit2) == NULL) {
    ABORTXLAL(status);
  }

  /* checkUnit now holds the units of f^-3*P1HW^-1*P2HW^-1) */
  if (XLALUnitMultiply(&tmpUnit1, &checkUnit, &(input->omegaGW->sampleUnits)) == NULL) {
    ABORTXLAL(status);
  }

  /* tmpUnit1 now holds units of Omega*f^-3*P1HW^-1*P2HW^-1) */

  if (XLALUnitMultiply(&tmpUnit2, &tmpUnit1, &(input->overlapReductionFunction->sampleUnits)) == NULL) {
    ABORTXLAL(status);
  }

  /* tmpUnit2 now holds units of mygamma*Omega*f^-3*P1HW^-1*P2HW^-1) */

  if (XLALUnitMultiply(&(optimalFilter->sampleUnits), &(lambda->units), &tmpUnit2) == NULL) {
    ABORTXLAL(status);
  }

  /* Done with unit manipulation */

  optimalFilter->data->data[0] = 0.0;

  /* calculate optimal filter values */
  for (i = (f0 == 0 ? 1 : 0) ; i < length; ++i)
  {
    f = f0 + deltaF * (REAL8)i;

    f3 = f * f * f;

    omega = input->omegaGW->data->data[i];
    mygamma = input->overlapReductionFunction->data->data[i];
    p1HWInv = input->halfCalibratedInverseNoisePSD1->data->data[i];
    p2HWInv = input->halfCalibratedInverseNoisePSD2->data->data[i];

    cPtrOptimalFilter = &(optimalFilter->data->data[i]);

    realFactor = (mygamma * omega * lambda->value) / f3;

    *(cPtrOptimalFilter) = crectf( realFactor * ((crealf(p1HWInv) * crealf(p2HWInv)) + (cimagf(p1HWInv) * cimagf(p2HWInv))), realFactor * ((crealf(p1HWInv) * cimagf(p2HWInv)) - (cimagf(p1HWInv) * crealf(p2HWInv))) );
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);
} /* LALStochasticOptimalFilter() */

/** @} */
