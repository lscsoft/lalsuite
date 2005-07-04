/********************************** <lalVerbatim file="ThresholdsCV">
Author: Flanagan, E. and Cannon, K.
$Id$
**************************************************** </lalVerbatim> */

#include <math.h>
#include <stdlib.h>
#include <lal/LALRCSID.h>

NRCSID(THRESHOLDSC, "$Id$");

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf.h>
#include <lal/FindRoot.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/XLALGSL.h>
#include <lal/Thresholds.h>
#include <lal/XLALError.h>


/******** <lalVerbatim file="ChisqCdfP"> ********/
REAL8 XLALChisqCdf(
	REAL8 chi2,
	REAL8 dof
)
/******** </lalVerbatim> ********/
{
	/*
	 * Cumulative Probability Distribution for Chi Squared
	 * distribution.
	 *  
	 *  returns probability that x_1^2 + .. x_dof^2 <= chi2, where x_1,
	 *  ..  x_dof are independent Gaussians of zero mean and unit
	 *  variance.  The integral expression is prob = int_0^{chi^2/2} dx
	 *  x^((n/2)-1) e^(-x) / Gamma(n/2), where n = dof = number of
	 *  degrees of freedom.  note chi2 = 2 * cal E, calE = variable
	 *  used in paper.
	 */

	static const char *func = "XLALChisqCdf";
	REAL8 prob;

	/* Arguments chi2 and dof must be non-negative */
	if((chi2 < 0.0) || (dof <= 0.0))
		XLAL_ERROR_REAL8(func, XLAL_EDOM);

	/* use GSL because our previous version sucked */
	XLAL_CALLGSL(prob = gsl_cdf_chisq_P(chi2, dof));

	/* Check that final answer is a legal probability.  The third test
	 * is necessary since there are some numbers x for which (x>0.0)
	 * evaluates as TRUE but for which 1/x evaluates to inf */
	if((prob < 0.0) || (prob > 1.0) || (1.0 / prob > LAL_REAL8_MAX))
		XLAL_ERROR_REAL8(func, XLAL_ERANGE);

	return (prob);
}


/******** <lalVerbatim file="OneMinusChisqCdfP"> ********/
REAL8 XLALOneMinusChisqCdf(
	REAL8 chi2,
	REAL8 dof
)
/******** </lalVerbatim> ********/
{
	/*
	 *  Cumulative Probability Distribution for Chi Squared
	 *  distribution.  Alternative version which is more accurate for
	 *  large rho.
	 *  
	 *  returns probability that x_1^2 + .. x_dof^2 >= chi2, where x_1,
	 *  ..  x_dof are independent Gaussians of zero mean and unit
	 *  variance.  The integral expression is prob =
	 *  int_{chi^2/2}^\infty dx  x^((n/2)-1) e^(-x) / Gamma(n/2), where
	 *  n = dof = number of degrees of freedom.  note chi2 = 2 * cal E,
	 *  calE = variable used in paper.
	 *  
	 *  This function's code is the same as the function XLALChisqCdf()
	 *  except for the very end.
	 */

	static const char *func = "XLALOneMinusChisqCdf";
	REAL8 prob;

	if ((chi2 < 0.0) || (dof <= 0.0))
		XLAL_ERROR_REAL8(func, XLAL_EDOM);

	/* Use GSL because our previous version sucked */
	XLAL_CALLGSL(prob = gsl_cdf_chisq_Q(chi2, dof));

	/* Check that final answer is a legal probability. Third test is
	 * necessary since there are some numbers x for which (x>0.0)
	 * evaluates as TRUE but for which 1/x evaluates to inf */
	if ((prob < 0.0) || (prob > 1.0))
		XLAL_ERROR_REAL8(func, XLAL_ERANGE);
	if (1.0 / prob > LAL_REAL8_MAX)
		prob = exp(-700.0);

	return (prob);
}


#if 0
/******** <lalVerbatim file="OneMinusChisqCdfP"> ********/
REAL8 XLALlnOneMinusChisqCdf(
	REAL8 chi2,
	REAL8 dof
)
/******** </lalVerbatim> ********/
{
	/*
	 * This function returns the natural logarithm of the result returned
	 * by XLALOneMinusChisqCdf(), i.e. ln(Q(chi^2, dof))
	 *
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

	static const char *func = "XLALlnOneMinusChisqCdf";
	REAL8 a = dof/2;
	REAL8 x = chi2/2;
	REAL8 ln_prob, term;
	int i;

	if((chi2 < 0.0) || (dof <= 0.0))
		XLAL_ERROR_REAL8(func, XLAL_EDOM);

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
		XLAL_ERROR_REAL8(func, XLAL_ERANGE);

	return(ln_prob);
}
#else
REAL8 XLALlnOneMinusChisqCdf(
	REAL8 chi2,
	REAL8 dof
)
{
	/*
	 *  Cumulative Probability Distribution for Chi Squared
	 *  distribution.  Alternative version which is more accurate for
	 *  large rho.
	 *  
	 *  returns probability that x_1^2 + .. x_dof^2 >= chi2, where x_1,
	 *  ..  x_dof are independent Gaussians of zero mean and unit
	 *  variance.  The integral expression is prob =
	 *  int_{chi^2/2}^\infty dx  x^((n/2)-1) e^(-x) / Gamma(n/2), where
	 *  n = dof = number of degrees of freedom.  note chi2 = 2 * cal E,
	 *  calE = variable used in paper.
	 *  
	 *  This function's code is the same as the function XLALChisqCdf()
	 *  except for the very end.
	 */

	static const char *func = "XLALOneMinusChisqCdf";
	REAL8 prob;

	if ((chi2 < 0.0) || (dof <= 0.0))
		XLAL_ERROR_REAL8(func, XLAL_EDOM);

	/* Use GSL because our previous version sucked */
	XLAL_CALLGSL(prob = gsl_cdf_chisq_Q(chi2, dof));

	/* Check that final answer is a legal probability. Third test is
	 * necessary since there are some numbers x for which (x>0.0)
	 * evaluates as TRUE but for which 1/x evaluates to inf */
	if ((prob < 0.0) || (prob > 1.0))
		XLAL_ERROR_REAL8(func, XLAL_ERANGE);
	if (1.0 / prob > LAL_REAL8_MAX)
		prob = -700.0;
	else
		prob = log(prob);

	return (prob);
}
#endif


/* returns n! */
static REAL8 Factorial(INT4 n)
{
	REAL8 nfactorial = 1.0;

	for(; n >= 2; n--)
		nfactorial *= n;

	return(nfactorial);
}


/******** <lalVerbatim file="NoncChisqCdfP"> ********/
REAL8 XLALNoncChisqCdfNonSafe(
	REAL8 chi2,
	REAL8 dof,
	REAL8 nonCentral
)
/******** </lalVerbatim> ********/
{
	/*
	 *  Cumulative distribution function for noncentral chi-squared
	 *  distribution
	 *  
	 *  returns probability that (x_1+rho)^2 + x_2^2 + .. x_dof^2 \le
	 *  chi2, where x_1, ..  x_dof are independent Gaussians of zero
	 *  mean and unit variance, and where nonCentral = rho^2 and dof =
	 *  number of degrees of freedom
	 *
	 *  We use the series formula to evaluate the probability.  Each
	 *  term in the series involves a call to XLALChisqCdf().
	 */

	static const char *func = "XLALNoncChisqCdfNonSafe";
	const REAL8 fractionalAccuracy = 1.0e-5;
	const INT4 maxloop = 170;	/* Factorial() breaks down here */
	REAL8 temp;
	REAL8 current;
	REAL8 sum;
	INT4 n;

	/* Arguments chi2, dof and nonCentral must be non-negative */
	if((dof <= 0.0) ||
	   (chi2 < 0.0) ||
	   (nonCentral < 0.0))
		XLAL_ERROR_REAL8(func, XLAL_EDOM);

	/* Evaluate the first term in the series   */
	temp = XLALChisqCdf(chi2, dof);
	if(XLALIsREAL8FailNaN(temp))
		XLAL_ERROR_REAL8(func, XLAL_EFUNC);
	sum = temp * exp(-nonCentral / 2.0);

	/* Now add successive higher terms in the series until either
	 * sufficient accuracy is achieved, or we exceed the maximum allowed
	 * number of terms */

	n = 0;
	do {
		if(++n >= maxloop)
			XLAL_ERROR_REAL8(func, XLAL_EMAXITER);
		temp = XLALChisqCdf(chi2, dof + 2.0 * n);
		if(XLALIsREAL8FailNaN(temp))
			XLAL_ERROR_REAL8(func, XLAL_EFUNC);
		current = exp(-nonCentral / 2.0 + n * log(nonCentral / 2.0)) * temp / Factorial(n);
		sum += current;
	}
	while(fabs(current / sum) > fractionalAccuracy);

	return(sum);
}


/******** <lalVerbatim file="NoncChisqCdfP"> ********/
REAL8 XLALNoncChisqCdf(
	REAL8 chi2,
	REAL8 dof,
	REAL8 nonCentral
)
/******** </lalVerbatim> ********/
{
	static const char *func = "XLALNoncChisqCdf";
	REAL8 prob = XLALNoncChisqCdfNonSafe(chi2, dof, nonCentral);

	if(XLALIsREAL8FailNaN(prob))
		XLAL_ERROR_REAL8(func, XLAL_EFUNC);

	/* check that final answer is a legal probability.  third test is
	 * necessary since there are some numbers x for which (x > 0.0)
	 * evaluates as TRUE but for which 1/x evaluates to inf */

	if((prob <= 0.0) || (prob >= 1.0) || (1.0 / prob >= LAL_REAL8_MAX))
		XLAL_ERROR_REAL8(func, XLAL_ERANGE);
	return(prob);
}


/******** <lalVerbatim file="Chi2ThresholdP"> ********/
REAL8 XLALChi2Threshold(
	REAL8 dof,
	REAL8 falseAlarm
)
/******** </lalVerbatim> ********/
{
	/* threshold for chi2:  returns value of chi2 such that falseAlarm = 1
	 * - chisqCdf(chi2,dof) */
	static const char *func = "XLALChi2Threshold";
	REAL8 chi2;

	/* Argument dof must be positive, and supplied false alarm probability
	 * must be between 0 and 1 */
	if((dof <= 0.0) ||
	   (falseAlarm <= 0.0) || (falseAlarm >= 1.0))
		XLAL_ERROR_REAL8(func, XLAL_EDOM);

	/* call GSL */
	XLAL_CALLGSL(chi2 = gsl_cdf_chisq_Qinv(falseAlarm, dof));

	return(chi2);
}


struct NoncChisqCdfParams
{
	REAL8 chi2;
	REAL8 dof;
	REAL8 falseDismissal;
};


static REAL8 NoncChisqCdf(REAL8 lnrho, void *data)
{
	/* Wrap XLALNoncChisqCdfNonSafe() for use with the root-finding
	 * functions in XLALRhoThreshold() below. */
	struct NoncChisqCdfParams *params = data;
	REAL8 chi2 = params->chi2;
	REAL8 dof = params->dof;
	REAL8 nonCentral = exp(2.0 * lnrho);
	REAL8 falseDismissal = params->falseDismissal;

	/* use the "nonsafe" version because it's ok for the result to be 1 or
	 * 0 */
	return(XLALNoncChisqCdfNonSafe(chi2, dof, nonCentral) - falseDismissal);
}


/******** <lalVerbatim file="RhoThresholdP"> ********/
REAL8 XLALRhoThreshold(
	REAL8 chi2,
	REAL8 dof,
	REAL8 falseDismissal
)
/******** </lalVerbatim> ********/
{
	/*
	 * threshold for rho:  returns value of rho such that falseAlarm =
	 * noncChisqCdf(chi2,dof,rho^2) note that rho^2 is the same as
	 * nonCentral
	 */
	static const char *func = "XLALRhoThreshold";
	REAL8 xmin = -2.0;
	REAL8 xmax = +2.0;
	struct NoncChisqCdfParams params;

	/* Arguments dof and chi2 must be positive, supplied false
	 * dismissal probability must be between 0 and 1 */
	if((dof <= 0.0) ||
	   (chi2 < 0.0) ||
	   (falseDismissal <= 0.0) ||
	   (falseDismissal >= 1.0))
		XLAL_ERROR_REAL8(func, XLAL_EDOM);

	/* Setup NoncChisqCdf1() parameters */
	params.chi2 = chi2;
	params.dof = dof;
	params.falseDismissal = falseDismissal;

	/* Now bracket and find the root */
	XLALDBracketRoot(NoncChisqCdf, &xmin, &xmax, &params);

	return(exp(XLALDBisectionFindRoot(NoncChisqCdf, xmin, xmax, 1e-5, &params)));
}
