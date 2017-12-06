/*
 *
 * Copyright (C) 2007  Kipp Cannon and Flanagan, E.
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
#include <stdlib.h>


#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf.h>
#include <lal/FindRoot.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/XLALGSL.h>
#include <lal/Thresholds.h>
#include <lal/XLALError.h>

/**
 * Cumulative Probability Distribution for Chi Squared distribution.
 *
 * returns probability that x_1^2 + .. x_dof^2 <= chi2, where x_1, ..
 * x_dof are independent Gaussians of zero mean and unit variance.  The
 * integral expression is prob = int_0^{chi^2/2} dx x^((n/2)-1) e^(-x) /
 * Gamma(n/2), where n = dof = number of degrees of freedom.  note chi2 = 2
 * * cal E, calE = variable used in paper.
 */
REAL8 XLALChisqCdf(
	REAL8 chi2,
	REAL8 dof
)
{
	double prob;

	/* Arguments chi2 and dof must be non-negative and positive
	 * respectively */
	if((chi2 < 0.0) || (dof <= 0.0))
		XLAL_ERROR_REAL8(XLAL_EDOM);

	/* use GSL because our previous version sucked */
	XLAL_CALLGSL(prob = gsl_cdf_chisq_P(chi2, dof));

	/* Check that final answer is a legal probability.  */
	if((prob < 0.0) || (prob > 1.0))
		XLAL_ERROR_REAL8(XLAL_ERANGE);

	return prob;
}


/**
 * Cumulative Probability Distribution for Chi Squared distribution.
 * Alternative version which is more accurate for large rho.
 *
 * returns probability that \f$x_1^2 + .. x_dof^2 >= chi2\f$, where \f$x_1, ..\f$
 * x_dof are independent Gaussians of zero mean and unit variance.  The
 * integral expression is \f$prob = \int_{chi^2/2}^\infty dx  x^((n/2)-1)
 * e^(-x) / Gamma(n/2)\f$, where n = dof = number of degrees of freedom.  note
 * chi2 = 2 * cal E, calE = variable used in paper.
 */
REAL8 XLALOneMinusChisqCdf(
	REAL8 chi2,
	REAL8 dof
)
{
	double prob;

	if((chi2 < 0.0) || (dof <= 0.0))
		XLAL_ERROR_REAL8(XLAL_EDOM);

	/* Use GSL because our previous version sucked */
	XLAL_CALLGSL(prob = gsl_cdf_chisq_Q(chi2, dof));

	/* Check that final answer is a legal probability. */
	if((prob < 0.0) || (prob > 1.0))
		XLAL_ERROR_REAL8(XLAL_ERANGE);

	return prob;
}


/**
 * This function returns the natural logarithm of the result returned by
 * XLALOneMinusChisqCdf(), i.e. ln(Q(chi^2, dof))
 */
REAL8 XLALlnOneMinusChisqCdf(
	REAL8 chi2,
	REAL8 dof
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
	 * See Abramowitz and Stegun, (6.5.32).
	 */

	const REAL8 a = dof / 2;
	const REAL8 x = chi2 / 2;
	REAL8 ln_prob, term;
	int i;

	if((chi2 < 0.0) || (dof <= 0.0))
		XLAL_ERROR_REAL8(XLAL_EDOM);

	/* start with a high precision technique for large probabilities */
	XLAL_CALLGSL(ln_prob = log(gsl_cdf_chisq_Q(chi2, dof)));

	/* if the result turned out to be very small, then recompute it using
	 * an asymptotic approximation */
	if(ln_prob < -600) {
		/* calculate the log of the numerator */
		for(i = 1, ln_prob = term = 1; (i < a) && fabs(term/ln_prob) > 1e-17; i++) {
			term *= (a - i) / x;
			ln_prob += term;
		}
		ln_prob = (a - 1) * log(x) - x + log(ln_prob);

		/* subtract the log of the denominator */
		XLAL_CALLGSL(ln_prob -= gsl_sf_lngamma(a));
	}

	/* check that the final answer is the log of a legal probability */
	if(ln_prob > 0.0)
		XLAL_ERROR_REAL8(XLAL_ERANGE);

	return ln_prob;
}


/* returns n! */
static REAL8 Factorial(INT4 n)
{
	REAL8 nfactorial = 1.0;

	for(; n >= 2; n--)
		nfactorial *= n;

	return nfactorial;
}


/**
 * Cumulative distribution function for noncentral chi-squared distribution
 *
 * returns probability that \f$(x_1+rho)^2 + x_2^2 + .. x_dof^2 \le chi2\f$
 * where x_1, ..  x_dof are independent Gaussians of zero mean and unit
 * variance, and where nonCentral = rho^2 and dof = number of degrees of
 * freedom
 *
 * We use the series formula to evaluate the probability.  Each term in the
 * series involves a call to XLALChisqCdf().
 */
REAL8 XLALNoncChisqCdf(
	REAL8 chi2,
	REAL8 dof,
	REAL8 nonCentral
)
{
	const double epsilon = 1.0e-8;
	const int maxloop = 170;	/* Factorial() breaks down here */
	double term;
	double sum;
	int n;

	/* Arguments chi2, dof and nonCentral must be non-negative */
	if((dof <= 0.0) ||
	   (chi2 < 0.0) ||
	   (nonCentral < 0.0))
		XLAL_ERROR_REAL8(XLAL_EDOM);

	/* Add terms from the series until either sufficient accuracy is
	 * achieved, or we exceed the maximum allowed number of terms */

	sum = 0;
	n = 0;
	do {
		double P = XLALChisqCdf(chi2, dof + 2.0 * n);
		if(XLALIsREAL8FailNaN(P))
			XLAL_ERROR_REAL8(XLAL_EFUNC);
		sum += term = exp(-nonCentral / 2.0 + n * log(nonCentral / 2.0)) * P / Factorial(n);
		if(++n >= maxloop)
			XLAL_ERROR_REAL8(XLAL_EMAXITER);
	} while(fabs(term / sum) > epsilon);

	/* check that final answer is a legal probability. */
	if((sum < 0.0) || (sum > 1.0))
		XLAL_ERROR_REAL8(XLAL_ERANGE);

	return sum;
}


/**
 * Threshold for chi2.  Returns value of chi2 such that
 *
 * falseAlarm = 1 - chisqCdf(chi2,dof)
 */
REAL8 XLALChi2Threshold(
	REAL8 dof,
	REAL8 falseAlarm
)
{
	REAL8 chi2;

	/* Argument dof must be positive, and supplied false alarm probability
	 * must be between 0 and 1 */
	if((dof <= 0.0) ||
	   (falseAlarm <= 0.0) || (falseAlarm >= 1.0))
		XLAL_ERROR_REAL8(XLAL_EDOM);

	/* call GSL */
	XLAL_CALLGSL(chi2 = gsl_cdf_chisq_Qinv(falseAlarm, dof));

	return(chi2);
}


/* Wrap XLALNoncChisqCdf() for use with the root-finding functions in
 * XLALRhoThreshold() below. */


struct NoncChisqCdfParams {
	REAL8 chi2;
	REAL8 dof;
	REAL8 falseDismissal;
};


static REAL8 NoncChisqCdf(REAL8 lnrho, void *data)
{
	struct NoncChisqCdfParams *params = data;

	return XLALNoncChisqCdf(params->chi2, params->dof, exp(2.0 * lnrho)) - params->falseDismissal;
}


/**
 * Threshold for rho.  Returns value of rho such that
 *
 * falseAlarm = noncChisqCdf(chi2,dof,rho^2)
 *
 * note that rho^2 is the same as nonCentral.
 */
REAL8 XLALRhoThreshold(
	REAL8 chi2,
	REAL8 dof,
	REAL8 falseDismissal
)
{
	REAL8 xmin = -2.0;
	REAL8 xmax = +2.0;
	struct NoncChisqCdfParams params;

	/* Arguments dof and chi2 must be positive, supplied false
	 * dismissal probability must be between 0 and 1 */
	if((dof <= 0.0) ||
	   (chi2 < 0.0) ||
	   (falseDismissal <= 0.0) ||
	   (falseDismissal >= 1.0))
		XLAL_ERROR_REAL8(XLAL_EDOM);

	/* Setup NoncChisqCdf() parameters */
	params.chi2 = chi2;
	params.dof = dof;
	params.falseDismissal = falseDismissal;

	/* Now bracket and find the root */
	XLALDBracketRoot(NoncChisqCdf, &xmin, &xmax, &params);

	return exp(XLALDBisectionFindRoot(NoncChisqCdf, xmin, xmax, 1e-5, &params));
}
