/*
*  Copyright (C) 2007 Jolien Creighton, Kaice T. Reilly, Robert Adam Mercer, John Whelan
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
 * \addtogroup StochasticOptimalFilterNormalization_c
 *
 * \brief Calculates the normalization factor for the optimal filter and the
 * expected variance per unit time of the standard cross-correlation
 * statistic.
 *
 * As described in \ref StochasticOptimalFilter_c,
 * the optimal filter for stochastic searches is defined as
 *
 * \f{equation}{
 * \label{stochastic_e_Q}
 * \widetilde{Q}{}^{\mathrm{C}}(f)=\lambda\,
 * \frac{\gamma(f)\,\Omega_{\mathrm{GW}}(f)}
 * {|f|^3\,P^{\mathrm{C}}_1(f)\,P^{\mathrm{C}}_2(f)}
 * \f}
 *
 * The normalization constant \f$\lambda\f$ is chosen so that the expected mean
 * value of the cross-correlation statistic is \cite Allen1999
 * \f{equation}{
 * \label{stochastic_e_mu}
 * \mu = \frac{3 {H_0}^2}{20\pi^2}\, T \,\overline{w_1w_2}
 * \int_{-\infty}^{\infty} df\, |f|^{-3}\,
 * \gamma(f)\,\Omega_{\mathrm{GW}}(f)
 * \widetilde{Q}{}^{\mathrm{C}}(f) = \Omega_{\mathrm{R}} T
 * \f}
 * where \f$T\f$ is the integration time
 * (cf. \eqref{stochastic_e_ymax}, \f$w_1\f$ and \f$w_2\f$ are the functions
 * used to window the data, and
 * \f$\Omega_{\mathrm{R}} =\Omega_{\mathrm{GW}}(f_{\mathrm{R}})\f$ is the overall strength of the
 * stochastic background (see \ref OverlapReductionFunction_c).  This sets the
 * value at
 * \f{equation}{
 * \label{stochastic_e_lambda}
 * \lambda = \frac{20\pi^2\, \Omega_{\mathrm{R}}}
 * {3\,{H_0}^2 \overline{w_1w_2}}
 * \left(
 * \int_{-\infty}^\infty \frac{df}{f^6}
 * \frac{[\gamma(f)\,\Omega_{\mathrm{GW}}(f)]^2}{P^{\mathrm{C}}_1(f)P^{\mathrm{C}}_2(f)}
 * \right)^{-1}
 * \f}
 *
 * The same integral used to calculate \f$\lambda\f$ also allows one to
 * calculate the expected variance per unit integration time of the
 * cross-correlation statistic, since
 * \f{eqnarray}{
 * \label{stochastic_e_variance}
 * \frac{\sigma^2}{T}
 * &=& \frac{\overline{(w_1w_2)^2}}{4T}\int_{-\infty}^{\infty} df
 * \, P^{\mathrm{C}}_1(f)\, P^{\mathrm{C}}_2(f)\,
 * \left(
 * \widetilde{Q}{}^{\mathrm{C}}(f)
 * \right)^2
 * = \frac{\lambda}{4T} \overline{(w_1w_2)^2}
 * \int_{-\infty}^{\infty} \frac{df}{|f|^3}\,\gamma(f)\,
 * \Omega_{\mathrm{GW}}(f)\,\widetilde{Q}{}^{\mathrm{C}}(f) \\
 * &=& \frac{5\pi^2}{3 {H_0}^2}
 * \,\frac{\overline{(w_1w_2)^2}}{\overline{w_1w_2}}
 * \,\Omega_{\mathrm{R}} \,\lambda
 * \f}
 * where we have used \eqref{stochastic_e_Q} to replace one of the two
 * factors of \f$\widetilde{Q}{}^{\mathrm{C}}(f)\f$ and \eqref{stochastic_e_mu} to replace
 * the integral.
 *
 * <tt>LALStochasticOptimalFilterNormalization()</tt> uses
 * \eqref{stochastic_e_lambda} to calculate the normalization constant
 * \f$\lambda\f$ and \eqref{stochastic_e_variance} to calculate the expected
 * variance per unit time \f$\sigma^2/T\f$ of the cross-correlation
 * statistic.
 *
 * ### Algorithm ###
 *
 * The routine <tt>LALStochasticOptimalFilterNormalization()</tt> first uses
 * \eqref{stochastic_e_lambda} to find the normalization constant
 * \f$\lambda\f$ (the amplitude \f${h_{100}}^2\Omega_{\mathrm{R}}\f$ is
 * found by logarithmic interpolation using the reference frequency
 * \f$f_{\mathrm{R}}\f$ specified in the parameter structure and the
 * input series representing \f${h_{100}}^2\Omega_{\mathrm{GW}}(f)\f$).
 *
 * The precise behavior of the normalization depends on the boolean
 * parameter <tt>parameters->heterodyned</tt>, which indicates whether the
 * filter is to be used on heterodyned data or not.  In the case of
 * heterodyned data, the integral is approximated by the sum
 * \f{eqnarray}{
 * \lambda &\approx& \frac{20\pi^2\, \Omega_{\mathrm{R}}}
 * {3\,{H_0}^2 \overline{w_1w_2}}
 * \left(
 * \delta f\sum_{k=0}^{N-1}
 * (f_0 + k\,\delta f)^{-6}
 * \frac{(\gamma[k]\,\Omega_{\mathrm{GW}}[k])^2}{P^{\mathrm{C}}_1[k]P^{\mathrm{C}}_2[k]}
 * \right)^{-1} \\
 * &\approx&
 * \frac{20\pi^2\, \Omega_{\mathrm{R}}}{3\,{H_0}^2}
 * \left(
 * \int_{f_0}^{f_0+N\delta f} \frac{df}{f^6}
 * \frac{[\gamma(f)\,\Omega_{\mathrm{GW}}(f)]^2}{P^{\mathrm{C}}_1(f)P^{\mathrm{C}}_2(f)}
 * \right)^{-1}
 * \f}
 * (Leaving out frequencies outside the band is equivalent to assuming
 * one or both of the noise PSDs \f$P^{\mathrm{C}}_{1,2}(f)\f$ blows up outside that
 * range.)
 *
 * In the case of non-heterodyned data with \f$f_0=0\f$, we calculate
 * \f{equation}{
 * \lambda \approx \frac{20\pi^2\, \Omega_{\mathrm{R}}}
 * {3\,{H_0}^2 \overline{w_1w_2}}
 * \left(
 * \delta f\, 2\ {\mathrm{Re}}  \sum_{k=0 \mbox{ or } 1}^{N-1}
 * (k\,\delta f)^{-6}
 * \frac{(\gamma[k]\,\Omega_{\mathrm{GW}}[k])^2}{P^{\mathrm{C}}_1[k]P^{\mathrm{C}}_2[k]}
 * \right)^{-1}
 * \f}
 * which includes negative frequencies as well.  The difference
 * between the two is because the cross-correlation statistic appearing
 * in the definition \eqref{stochastic_e_mu} is the one calculated by
 * <tt>StochasticHeterodynedCrossCorrelationStatistic()</tt> in the case
 * of heterodyned and <tt>StochasticCrossCorrelationStatistic()</tt> in
 * the case of non-heterodyned data.
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
 * <li> The reference frequency \f$f_{\mathrm{R}}\f$ must lie
 * safely enough in the frequency range of the inputs to allow the
 * value of \f${h_{100}}^2\Omega_{\mathrm{R}}\f$ to be determined
 * by interpolation.</li>
 * <li> The implementation of the optimal filter function given here
 * assumes a large observation time continuum-limit approximation.  In
 * this limit, the Dirichlet kernels (which appear in an exact
 * expression for the standard cross-correlation statistic, when
 * evaluated in discrete time \cite Finn2001; see also
 * the documentation for the module Dirichlet.c in the utilities package)
 * may be replaced by Dirac delta functions.</li>
 * <li> The units of the input series are checked for consistency; since
 * \cite Allen1999
 * \f{equation}{
 * \langle\widetilde{h}{}^{\mathrm{C}}_1(f)^*
 * \widetilde{h}{}^{\mathrm{C}}_2(f')\rangle
 * = \frac{3H_0^2}{20\pi^2}\delta(f-f')
 * \gamma(|f|)\Omega_{\mathrm{GW}}(|f|)
 * \f}
 * and
 * \f{equation}{
 * \langle\widetilde{h}_i(f)^*\widetilde{h}_i(f')\rangle
 * = \delta(f-f')P^{\mathrm{C}}_i(f)
 * \,
 * \f}
 * we demand that, up to powers of ten,
 * \f{equation}{
 * ([\gamma][\Omega_{\mathrm{GW}}])^2
 * =([\widetilde{h}_1][\widetilde{h}_2][T]^{-2})^2
 * =[P^{\mathrm{C}}_1][P^{\mathrm{C}}_2][T]^{-2}
 * \f}</li>
 * <li> This routine, like all those in the \c stochastic package,
 * uses with single precision arithmetic, consistent with the accuracy
 * of the continuous approximation for the intended number of data
 * points.  However, the limited dynamic range can pose a problem,
 * especially in this routine, which uses numbers like \f$H_0^2\f$ which
 * are far from unity when expressed in the standard SI units.  To
 * avoid this problem, the application of \eqref{stochastic_e_lambda}
 * takes the units \f$H_0\f$ of to be \f$10^{-18}\,\textrm{s}^{-1}\f$ rather
 * than \f$\textrm{s}^{-1}\f$.  In these units \f$h_{100}H_0\f$ has a numerical
 * value of \f$3.2407792903\f$, and inclusion of \f$(h_{100}H_0)^2\f$ in an
 * expression doesn't risk single-precision overflow.  When %'
 * constructing the unit structures of its outputs,
 * <tt>LALStochasticOptimalFilterNormalization()</tt> uses the
 * power-of-ten feature of the \c LALUnit structure to account
 * for the units of \f$H_0\f$.</li>
 * <li> The expected units for the inputs and outputs of this function
 * are as follows (although the actual output units will be constructed
 * from the input units):
 * \f{equation}{
 * ([\gamma][\Omega_{\mathrm{GW}}])^2
 * =([\widetilde{h}_1][\widetilde{h}_2][T]^{-2})^2
 * =[P^{\mathrm{C}}_1][P^{\mathrm{C}}_2][T]^{-2}
 * \f}
 * \f{eqnarray}{
 * {} [\gamma] &=& \textrm{strain}^{2} \\
 * {} [\Omega_{\mathrm{GW}}] &=& 1 \\
 * {} [1/P^{\mathrm{C}}_{1,2}] &=& 10^{36}\,\textrm{Hz}\,\textrm{strain}^{-2} \\
 * {} [\lambda] &=&
 * 10^{36}\, [P^{\mathrm{C}}_1]\,[P^{\mathrm{C}}_2]
 * \,[\gamma]^{-2}\,[\Omega_{\mathrm{GW}}]^{-1}\,\textrm{s}^{-3}
 * =  10^{-36}\,\textrm{s}^{-1}\\
 * {} [\sigma^2/T] &=&
 * 10^{36}\, [\lambda]\, [\Omega_{\mathrm{GW}}] \,\textrm{s}^2
 * = \textrm{s}
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
LALStochasticOptimalFilterNormalization(
    LALStatus                                            *status,
    StochasticOptimalFilterNormalizationOutput           *output,
    const StochasticOptimalFilterNormalizationInput      *input,
    const StochasticOptimalFilterNormalizationParameters *parameters)

{
  REAL4 omegaTimesGamma;
  REAL4 p1Inv;
  REAL4 p2Inv;

  REAL8 f;
  REAL8 f0;
  REAL8 f3;
  REAL8 deltaF;

  /* normalization factor */
  UINT4 i;
  UINT4 xRef;
  REAL4 omega1;
  REAL4 omega2;
  REAL8 freq1;
  REAL8 freq2;
  REAL4 exponent;
  REAL4 omegaRef;
  REAL8 lambdaInv;
  REAL8 f6;

  /* windowing */
  REAL4 w2, meanW2, meanW4;
  REAL4 *sPtrW1, *sPtrW2, *sStopPtr;
  UINT4 winLength = 0;

  UINT4 length;

  LALUnit tmpUnit1, tmpUnit2, checkUnit;
  REAL4WithUnits *lamPtr;

  /* initialize status pointer */
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  /* ERROR CHECKING ----------------------------------------------------- */

  /* check for null pointers *****/
  /* input structure */
  ASSERT(input != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* output structure */
  ASSERT(output != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* parameter structure */
  ASSERT(parameters != NULL, status, \
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

  /* inverse noise 1 of input */
  ASSERT(input->inverseNoisePSD1 != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* inverse noise 2 of input */
  ASSERT(input->inverseNoisePSD2 != NULL, status, \
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

  /* data member of inverse noise 1 */
  ASSERT(input->inverseNoisePSD1->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member of inverse noise 2 */
  ASSERT(input->inverseNoisePSD2->data != NULL, status, \
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

  /* data-data member of inverse noise 1 */
  ASSERT(input->inverseNoisePSD1->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member of inverse noise 2 */
  ASSERT(input->inverseNoisePSD2->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* Checks that only apply if windowing is specified */

  if (parameters->window1 || parameters->window2)
  {
    /* window 1 parameter */
    ASSERT(parameters->window1 != NULL, status, \
        STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
        STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* window 2 parameter */
    ASSERT(parameters->window2 != NULL, status, \
        STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
        STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* data member of window 1 */
    ASSERT(parameters->window1->data != NULL, status, \
        STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
        STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* data member of window 2 */
    ASSERT(parameters->window2->data != NULL, status, \
        STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
        STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    winLength = parameters->window1->length;

    ASSERT(winLength != 0, status, \
        STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

    /* Check that windows are the same length */
    if (parameters->window2->length != winLength)
    {
      ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
          STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
    }
  }

  /* done with null pointers ***/

  /* extract parameters from overlap */
  length = input->overlapReductionFunction->data->length;
  f0 = input->overlapReductionFunction->f0;
  deltaF = input->overlapReductionFunction->deltaF;

  /* check for legality ****/
  /* length must be positive */
  ASSERT(length != 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* start frequency must not be negative */
  if (f0 < 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);
  }

  /* frequency spacing must be positive */
  ASSERT(deltaF > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  /* check for mismatches **/
  /* length */
  if (input->omegaGW->data->length != length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  if (input->inverseNoisePSD1->data->length != length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  if (input->inverseNoisePSD2->data->length != length)
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
  if (input->inverseNoisePSD1->f0 != f0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }
  if (input->inverseNoisePSD2->f0 != f0)
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
  if (input->inverseNoisePSD1->deltaF != deltaF)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }
  if (input->inverseNoisePSD2->deltaF != deltaF)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }
  /* check for reference frequency lower and upper limits **/
  if (parameters->fRef < f0 + deltaF)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EOORFREF, \
        STOCHASTICCROSSCORRELATIONH_MSGEOORFREF);
  }
  if (parameters->fRef > f0 + ((REAL8)(length-1)*deltaF))
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EOORFREF, \
        STOCHASTICCROSSCORRELATIONH_MSGEOORFREF);
  }

  /* EVERYHTING OKAY HERE! ---------------------------------------------- */

  /* check that units of gamma, Omega, P1 and P2 are consistent
   * we must have (gamma*Omega)^2 with the same units as f^2*P1*P2
   * up to a power of ten. We check this by constructing
   * checkUnit = f^2*P1*P2*(gamma*Omega)^-2 */
  if (XLALUnitMultiply(&tmpUnit1, &(input->inverseNoisePSD1->sampleUnits), &(input->inverseNoisePSD2->sampleUnits)) == NULL) {
    ABORTXLAL(status);
  }

  /* tmpUnit1 holds units of 1/(P1*P2) */

  if (XLALUnitRaiseINT2(&tmpUnit2, &tmpUnit1, -1) == NULL) {
    ABORTXLAL(status);
  }

  /* tmpUnit2 holds units of P1*P2 */

  if (XLALUnitMultiply(&tmpUnit1, &(input->overlapReductionFunction->sampleUnits), &(input->omegaGW->sampleUnits)) == NULL) {
    ABORTXLAL(status);
  }

  /* tmpUnit1 holds units of Omega*Gamma */

  if (XLALUnitMultiply(&checkUnit, &tmpUnit1, &lalSecondUnit) == NULL) {
    ABORTXLAL(status);
  }

  /* checkUnit holds units of f^-1*Omega*Gamma */

  if (XLALUnitRaiseINT2(&tmpUnit1, &checkUnit, -2) == NULL) {
    ABORTXLAL(status);
  }

  /* tmpUnit1 holds units of f^2*(Omega*Gamma)^-2 */

  if (XLALUnitMultiply(&checkUnit, &tmpUnit1, &tmpUnit2) == NULL) {
    ABORTXLAL(status);
  }

  /* checkUnit holds units of f^2*P1*P2(Omega*Gamma)^-2 */

  /* Check that checkUnit is dimensionless up to a power of ten ***/

  for (i=0; i<LALNumUnits; ++i)
  {
    if (checkUnit.unitNumerator[i] || checkUnit.unitDenominatorMinusOne[i])
    {
      ABORT(status, STOCHASTICCROSSCORRELATIONH_EWRONGUNITS, \
          STOCHASTICCROSSCORRELATIONH_MSGEWRONGUNITS);
    }
  }

  /* Set tmpUnit1 to dims of Omega/H0^2 ******/

  /* First, set it to dims of H0 */

  tmpUnit1 = lalHertzUnit;

  /* Account for scaled units of Hubble constant */
  tmpUnit1.powerOfTen -= 18;

  /* Now set tmpUnit2 to dims of H0^-2 */
  if (XLALUnitRaiseINT2(&tmpUnit2, &tmpUnit1, -2) == NULL) {
    ABORTXLAL(status);
  }
  if (XLALUnitMultiply(&tmpUnit1, &(input->omegaGW->sampleUnits), &tmpUnit2) == NULL) {
    ABORTXLAL(status);
  }

  /* Now tmpUnit1 has units of Omega/H0^2 */

  /* assign correct units to normalization constant ********/
  /* These are Omega/H0^2*f^5*P1*P2*(gamma*Omega)^-2
   * which is the same as tmpUnit1*f^3*checkUnit */
  if (XLALUnitMultiply(&tmpUnit2, &tmpUnit1, &checkUnit) == NULL) {
    ABORTXLAL(status);
  }

  /* tmpUnit2 has units of Omega/H0^2*f^2*P1*P2*(gamma*Omega)^-2 */

  /* Now that we've used checkUnit, we don't need it any more and */
  /* can use it for temp storage (of f^3) */
  if (XLALUnitRaiseINT2(&checkUnit, &lalHertzUnit, 3) == NULL) {
    ABORTXLAL(status);
  }

  /* In case the normalization output was NULL, we need to allocate it
   * since we still use it as an intermediate; we do so here to
   * minimize the BEGINFAIL/ENDFAIL pairs needed */

  if (output->normalization != NULL)
  {
    lamPtr = output->normalization;
  }
  else
  {
    lamPtr = (REAL4WithUnits*)LALMalloc(sizeof(REAL4WithUnits));
  }

  if (XLALUnitMultiply(&(lamPtr->units), &tmpUnit2, &checkUnit) == NULL) {
    if (output->normalization == NULL) LALFree(lamPtr);
    ABORTXLAL(status);
  }

  if (output->variance != NULL)
  {
    /* assign correct units to variance per time of CC stat ********/
    if (XLALUnitMultiply(&(output->variance->units), &(lamPtr->units), &tmpUnit1) == NULL) {
      if (output->normalization == NULL) LALFree(lamPtr);
      ABORTXLAL(status);
    }
  }

  /* Calculate properties of windows */
  if (parameters->window1 == NULL)
  {
    meanW2 = meanW4 = 1.0;
  }
  else
  {
    meanW2 = meanW4 = 0.0;
    for (sPtrW1 = parameters->window1->data, \
        sPtrW2 = parameters->window2->data, sStopPtr = sPtrW1 + winLength; \
        sPtrW1 < sStopPtr; ++sPtrW1, ++sPtrW2)
    {
      w2 = *sPtrW1 * *sPtrW2;
      meanW2 += w2;
      meanW4 += w2 * w2;
    }
    meanW2 /= winLength;
    meanW4 /= winLength;

    if (meanW2 <= 0)
    {
      ABORT(status, STOCHASTICCROSSCORRELATIONH_ENONPOSWIN, \
            STOCHASTICCROSSCORRELATIONH_MSGENONPOSWIN);
    }
  }

  /* ************* calculate lambda ***********************/
  /* find omegaRef */
  xRef = (UINT4)((parameters->fRef - f0)/deltaF);
  if (((parameters->fRef - f0)/deltaF) == (REAL8)(xRef))
  {
    omegaRef = input->omegaGW->data->data[xRef];
  }
  else
  {
    omega1 = input->omegaGW->data->data[xRef];
    omega2 = input->omegaGW->data->data[xRef+1];
    freq1 = f0 + xRef * deltaF;
    freq2 = f0 + (xRef + 1) * deltaF;
    if (omega1 <= 0)
    {
      ABORT(status, STOCHASTICCROSSCORRELATIONH_ENONPOSOMEGA, \
          STOCHASTICCROSSCORRELATIONH_MSGENONPOSOMEGA);
    }
    else
    {
      exponent = ((log((parameters->fRef)/freq1))/(log(freq2/freq1)));
    }
    omegaRef = omega1*(pow((omega2/omega1),exponent));
  }

  /* calculate inverse lambda value */
  lambdaInv = 0.0;

  for (i = (f0 == 0 ? 1 : 0) ; i < length; ++i)
  {
    f = f0 + deltaF * (REAL8) i;

    f3 = f * f * f;
    f6 = f3 * f3;

    omegaTimesGamma = input->omegaGW->data->data[i] * \
                      input->overlapReductionFunction->data->data[i];
    p1Inv = input->inverseNoisePSD1->data->data[i];
    p2Inv = input->inverseNoisePSD2->data->data[i];
    lambdaInv += (omegaTimesGamma * omegaTimesGamma * p1Inv * p2Inv) / f6;
  }

  lambdaInv /= (omegaRef / deltaF) * ((20.0L * LAL_PI * LAL_PI) / \
      (3.0L * (LAL_H0FAC_SI*1e+18) * (LAL_H0FAC_SI*1e+18) * meanW2));

  if (!parameters->heterodyned)
    lambdaInv *= 2.0;

  lamPtr->value = 1 / lambdaInv;

  if (output->variance != NULL)
  {
    output->variance->value = (meanW4/meanW2) * ((5.0L * LAL_PI * LAL_PI) / \
        (3.0L * (LAL_H0FAC_SI*1e+18) * (LAL_H0FAC_SI*1e+18))) * \
                              omegaRef / lambdaInv;
  }
  if (output->normalization == NULL)
    LALFree(lamPtr);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/** @} */
