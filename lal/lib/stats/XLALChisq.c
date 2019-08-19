/*
 * Copyright (C) 2005,2007,2008,2010,2015  Kipp Cannon
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


#include <math.h>


#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf.h>


#include <lal/LALChisq.h>
#include <lal/LALConstants.h>
#include <lal/XLALGSL.h>
#include <lal/XLALError.h>


/**
 * Compute the natural logarithm of the complementary cumulative
 * probability function of the \f$\chi^{2}\f$ distribution.
 *
 * Returns the natural logarithm of the probability that \f$x_{1}^{2} +
 * \cdots + x_{\mathrm{dof}}^{2} \geq \chi^{2}\f$, where \f$x_{1}, \ldots,
 * x_{\mathrm{dof}}\f$ are independent zero mean unit variance Gaussian
 * random variables.  The integral expression is \f$\ln Q = \ln
 * \int_{\chi^{2}/2}^{\infty} x^{\frac{n}{2}-1} \mathrm{e}^{-x} /
 * \Gamma(n/2) \, \mathrm{d}x\f$, where \f$n = \mathrm{dof} =\f$ number of
 * degrees of freedom.
 *
 * Results agree with Mathematica to 15 digits or more where the two have
 * been compared and where Mathematica is able to compute a result.  For
 * example, it does not seem to be possible to obtain a numerical result
 * from Mathematica for the equivalent of XLALLogChisqCCDF(10000.0, 8.5)
 * using any number of digits in the intermediate calculations, though the
 * two implementations agree to better than 15 digits at
 * XLALLogChisqCCDF(10000.0, 8.0)
 */

double XLALLogChisqCCDF(
	double chi2,
	double dof
)
{
	/*
	 * Notes:
	 *
	 * Q(chi^2, dof) = Gamma(dof/2, chi^2/2) / Gamma(dof/2)
	 *
	 * The incomplete gamma function, Gamma(a, x), can by approximated by
	 * computing the first few terms of
	 *
	 * Gamma(a, x) = e^-x x^(a-1) [ 1 + (a-1)/x + (a-1)(a-2)/x^2 + ... ]
	 *
	 * See Abramowitz and Stegun, (6.5.32).  The natural logarithm of
	 * this is
	 *
	 * ln Gamma(a, x) = (a - 1) * ln(x) - x + ln[ 1 + (a-1)/x +
	 * (a-1)(a-2)/x^2 + ... ]
	 *
	 * and we can use the log1p() function to evaluate the log of the
	 * sum.
	 */

	double ln_prob;

	if((chi2 < 0.0) || (dof <= 0.0))
		XLAL_ERROR_REAL8(XLAL_EDOM);

	/* start with a technique for intermediate probabilities */
	XLAL_CALLGSL(ln_prob = log(gsl_cdf_chisq_Q(chi2, dof)));

	if(ln_prob > -LAL_LN2) {
		/* the result turned out to be large, so recompute using
		 * the CDF */
		XLAL_CALLGSL(ln_prob = log1p(-gsl_cdf_chisq_P(chi2, dof)));
	} else if(ln_prob < -600.) {
		/* the result turned out to be very small, so recompute
		 * using an asymptotic approximation */
		const double a = dof / 2.;
		const double x = chi2 / 2.;
		double sum;
		double term;
		double i;

		/* calculate the log of the numerator.  if this aymptotic
		 * form doesn't work out, we can choose to report that as
		 * an error or we can fall back to whatever result we got
		 * from GSL.  this is left as a compile-time choice.  if
		 * the choice is to report this as an error, it can either
		 * be reported as the failure of the series to converge or
		 * as a general loss of precision in the result (because
		 * the only option left is to fall back to the GSL result).
		 * this is also left as a compile-time choice. */
		for(i = 1., sum = 0., term = 1.; fabs(term/sum) > 1e-17; i++) {
			double factor = (a - i) / x;
#if 0
#if 1
			if(fabs(factor) >= 1.)
				XLAL_ERROR_REAL8(XLAL_EDIVERGE);
			if(i > 100.)
				XLAL_ERROR_REAL8(XLAL_EMAXITER);
#else
			if(fabs(factor) >= 1. || i > 100.)
				XLAL_ERROR_REAL8(XLAL_ELOSS);
#endif
#else
			if(fabs(factor) >= 1. || i > 100.)
				goto done;
#endif
			term *= factor;
			sum += term;
		}
		ln_prob = (a - 1.) * log(x) - x + log1p(sum);

		/* subtract the log of the denominator */
		XLAL_CALLGSL(ln_prob -= gsl_sf_lngamma(a));
	}

	/* check that the final answer is the log of a legal probability */
done:
	if(ln_prob > 0.0)
		XLAL_ERROR_REAL8(XLAL_ERANGE);

	return ln_prob;
}
