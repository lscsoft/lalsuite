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
 * \addtogroup StochasticOmegaGW_c
 *
 * \brief Generates a frequency series containing a simple power law spectrum.
 *
 * The strength of a stochastic gravitational wave background is
 * defined as
 * \f{equation}{
 * \Omega_{\mathrm{GW}}(f)
 * :=\frac{f}{\rho_{\mathrm{crit}}}
 * \frac{d\rho_{\mathrm{GW}}}{df}
 * \ ,
 * \f}
 * where
 * \f{equation}{
 * \rho_{\mathrm{crit}} = \frac{3 H_0^2 c^2}{8\pi G}
 * \f}
 * is the critical density needed to close the universe.  Since the value
 * of \f$\rho_{\mathrm{crit}}\f$ depends on the observed value of the
 * Hubble constant \f$H_0\f$, it is traditional to remove this experimental
 * uncertainty from the definition by working with
 * \f${h_{100}}^2\Omega_{\mathrm{GW}}\f$, where \f$h_{100}\f$ is the Hubble
 * constant divided by
 * \f$100\,\textrm{km}\,\textrm{s}^{-1}\,\textrm{Mpc}^{-1}\f$.
 *
 * <tt>LALStochasticOmegaGW()</tt> generates a simple power law spectrum
 * \f{equation}{
 * {h_{100}}^2\Omega_{\mathrm{GW}}(f)
 * =\Omega_{\mathrm{R}}
 * \left(
 * \frac{f}{f_{\mathrm{R}}}
 * \right)^\alpha
 * \f}
 * The parameter <tt>parameters.omegaRef</tt> specifies the amplitude
 * \f${h_{100}}^2\Omega_{\mathrm{R}}\f$ of the spectrum.  This is simply
 * defined as the value of \f${h_{100}}^2\Omega_{\mathrm{GW}}(f)\f$ at the
 * reference frequency \f$f_{\mathrm{R}}\f$ which is specified in
 * <tt>parameters.omegaRef</tt>.
 *
 * \heading{Algorithm}
 *
 * <tt>LALStochasticOmegaGW()</tt> treats the constant spectrum \f$\alpha=0\f$ as a
 * special case, and simply sets every element of the output series to
 * \f$\Omega_{\mathrm{R}}\f$.
 *
 * If \f$\alpha\ne 0\f$, the DC value <tt>output->data->data[0]</tt> is set
 * to 0 or \c LAL_REAL4_MAX, depending on whether \f$\alpha\f$ is
 * positive or negative, respectively.
 *
 * The output units are set to be dimensionless.
 *
 * \heading{Uses}
 * \code
 * LAL_REAL4_MAX
 * \endcode
 *
 * \heading{Notes}
 *
 * <ul>
 * <li> This routine will eventually be generalized to
 * include "broken" power law spectra
 * \f{equation}{
 * {h_{100}}^2\Omega_{\mathrm{GW}}
 * = \left\{
 * \begin{array}{cc}
 * {h_{100}}^2\Omega_1 f^{\alpha_1} & f\le f_c\\
 * {h_{100}}^2\Omega_2 f^{\alpha_2} & f\ge f_c
 * \end{array}
 * \right.
 * \f}</li>
 * </ul>
 *
 * @{
 */

#include <lal/LALStdlib.h>

#include <lal/LALConstants.h>
#include <math.h>
#include <lal/StochasticCrossCorrelation.h>

void
LALStochasticOmegaGW(
    LALStatus                         *status,
    REAL4FrequencySeries              *output,
    const StochasticOmegaGWParameters *parameters)

{
  REAL4* sPtr;
  REAL4* sStopPtr;

  REAL4 alpha;
  REAL8 f0;
  REAL8 deltaF;
  REAL4 x, deltaX, x0;   /* x = f / fRef */
  REAL8 fRef;
  REAL4 omegaRef;
  UINT4 i;

  UINT4 length;

  INITSTATUS(status);

  /* ERROR CHECKING -------------------------------------------------- */

  /* check that pointer to input parameters is non-null */
  ASSERT(parameters != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that frequency spacing is greater than zero */
  deltaF = parameters->deltaF;
  ASSERT(deltaF > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  /* check that length parameter is greater than zero */
  length = parameters->length;
  ASSERT(length > 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* check that heterodyning doesn't include negative physical frequencies */

  f0 = parameters->f0;
  if (f0 < 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);
  }

  /* check that pointer to real frequency series for output is non-null */
  ASSERT(output != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to data member of real frequency series for
   * output is non-null */
  ASSERT(output->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that length of data member of real frequency series for output
   * equals length specified in input parameters */
  if (output->data->length != length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* check that pointer to data-data member of real frequency series for
   * output is non-null */
  ASSERT(output->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that the fRef value is positive */
  fRef = parameters->fRef;
  if (fRef <= 0.0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EOORFREF, \
        STOCHASTICCROSSCORRELATIONH_MSGEOORFREF);
  }

  /* check that omegaRef is larger than zero */
  omegaRef = parameters->omegaRef;
  if (omegaRef <= 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENONPOSOMEGA, \
        STOCHASTICCROSSCORRELATIONH_MSGENONPOSOMEGA);
  }

  /* EVERYTHING OKAY HERE! -------------------------------------------- */

  alpha = parameters->alpha;

  /* assign output to be dimensionless */
  output->sampleUnits = lalDimensionlessUnit;

  strncpy(output->name, "Gravitational wave strength OmegaGW", LALNameLength);

  /* assign parameters to frequency series */
  output->epoch.gpsSeconds = 0;
  output->epoch.gpsNanoSeconds = 0;
  output->f0 = f0;
  output->deltaF = deltaF;

  deltaX = deltaF / fRef;
  x0 = f0 / fRef;

  /* assign pointers */
  sStopPtr = output->data->data + length;

  /* calculate output(f) values */
  if (alpha == 0){
    for (sPtr = output->data->data; sPtr < sStopPtr; ++sPtr)
    {
      *sPtr = omegaRef;
    }
  }
  else
  {
    if (f0 == 0)
    {
      output->data->data[0] = (alpha>0 ? 0 : LAL_REAL4_MAX);
      for (i=1 ; i < length ; ++i)
      {
        x = deltaX * (REAL4) i;
        output->data->data[i] = omegaRef * pow(x,alpha);
      }
    }
    else
    {
      for (i=0 ; i < length ; ++i)
      {
        x = x0 + deltaX * (REAL4) i;
        output->data->data[i] = omegaRef * pow(x,alpha);
      }
    }
  }

  RETURN(status);
}

/** @} */
