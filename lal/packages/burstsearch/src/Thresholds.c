/********************************** <lalVerbatim file="ThresholdsCV">
Author: Flanagan, E
$Id$
**************************************************** </lalVerbatim> */

#include <math.h>
#include <stdlib.h>

#include <lal/LALRCSID.h>

NRCSID(THRESHOLDSC, "$Id$");

#include <gsl/gsl_cdf.h>
#include <lal/FindRoot.h>
#include <lal/LALConstants.h>
#include <lal/LALErrno.h>
#include <lal/XLALGSL.h>
#include <lal/Thresholds.h>
#include <lal/XLALError.h>

extern int lalDebugLevel;

/*
 *
 * Private functions not accessible outside this file Thresholds.c 
 *
 */

/* returns n! */
static REAL8 Factorial(INT4 n)
{
	REAL8 nfactorial = 1.0;

	for(; n >= 2; n--)
		nfactorial *= n;

	return(nfactorial);
}


static void NoncChisqCdf1(LALStatus *status, REAL8 *prob, REAL8 lnrho, void *params)
{
	/* This is just the XLALNoncChisqCdf() function plus an added
	 * constant.  Its used in LALRhoThreshold() in the call to
	 * DFindRoot().  params is actually a pointer to a structure
	 * RhoThresholdIn It passed to NoncChisqCdf1 as a pointer to void
	 * so that NoncChisqCdf1 will be a REAL8LALFunction and be passable
	 * to DFindRoot(); see header file FindRoot.h */
	REAL8 chi2 = ((RhoThresholdIn *) params)->chi2;
	REAL8 dof = ((RhoThresholdIn *) params)->dof;
	REAL8 nonCentral = exp(2.0 * lnrho);
	REAL8 falseDismissal = ((RhoThresholdIn *) params)->falseDismissal;

	/* use the "nonsafe" version because it's ok for the result to be 1 or
	 * 0 */
	*prob = XLALNoncChisqCdfNonSafe(chi2, dof, nonCentral) - falseDismissal;
}



/*
 *  Public functions
 */


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


/******** <lalVerbatim file="ChisqCdfP"> ********/
void LALChisqCdf (
	LALStatus  *status,
	REAL8      *prob,
	ChisqCdfIn *input
)
/******** </lalVerbatim> ********/
{  
	INITSTATUS (status, "LALChisqCdf", THRESHOLDSC);
	ATTATCHSTATUSPTR (status);

	/* check that arguments are reasonable */
	ASSERT(prob, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(input, status, LAL_NULL_ERR, LAL_NULL_MSG);

	*prob = XLALChisqCdf(input->chi2, input->dof);

	/* raise a range error if the result is no good */
	if(XLALIsREAL8FailNaN(*prob)) {
		XLALClearErrno();
		ABORT(status, LAL_RANGE_ERR, LAL_RANGE_MSG);
	}

	DETATCHSTATUSPTR( status );
	RETURN (status);
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
	REAL8 temp;
	REAL8 current;
	REAL8 sum;
	INT4 n;
	const REAL8 fractionalAccuracy = 1.0e-10;
	const INT4 maxloop = 1000;

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
	while((current * current) / (sum * sum) > fractionalAccuracy);

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

	/* check that final answer is a legal probability third test is
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

/******** <lalVerbatim file="Chi2ThresholdP"> ********/
void LALChi2Threshold(
	LALStatus *status,
	REAL8 *chi2,
	Chi2ThresholdIn *input
)
/******** </lalVerbatim> ********/
{
	INITSTATUS(status, "LALChi2Threshold", THRESHOLDSC);
	ATTATCHSTATUSPTR(status);

	*chi2 = XLALChi2Threshold(input->dof, input->falseAlarm);

	DETATCHSTATUSPTR(status);
	RETURN(status);
}


/******** <lalVerbatim file="RhoThresholdP"> ********/
void LALRhoThreshold(
	LALStatus *status,
	REAL8 *rho,
	RhoThresholdIn *input
)
/******** </lalVerbatim> ********/
{
	/*
	 * threshold for rho:  returns value of rho such that falseAlarm =
	 * noncChisqCdf(chi2,dof,rho^2) note that rho^2 is the same as
	 * nonCentral
	 */

	REAL8 lnrhoAns;
	DFindRootIn frInput;

	INITSTATUS(status, "LALRhoThreshold", THRESHOLDSC);
	ATTATCHSTATUSPTR(status);

	/* check that pointers are not null */
	ASSERT(rho, status, LAL_NULL_ERR, LAL_NULL_MSG);
	ASSERT(input, status, LAL_NULL_ERR, LAL_NULL_MSG);

	/* Arguments dof and chi2 must be positive */
	ASSERT((input->dof > 0.0) && (input->chi2 >= 0.0), status, LAL_RANGE_ERR, LAL_RANGE_MSG);

	/* Supplied false dismissal probability must be between 0 and 1 */
	ASSERT((input->falseDismissal > 0.0) && (input->falseDismissal < 1.0), status, LAL_RANGE_ERR, LAL_RANGE_MSG);


	/* Initialize input structure for DFindRoot() */
	frInput.function = NoncChisqCdf1;
	frInput.xmin = -2.0;
	frInput.xmax = 2.0;
	frInput.xacc = 1e-5;


	/* Now bracket and find the root */
	LALDBracketRoot(status->statusPtr, &frInput, input);
	CHECKSTATUSPTR(status);

	LALDBisectionFindRoot(status->statusPtr, &lnrhoAns, &frInput, input);
	CHECKSTATUSPTR(status);

	*rho = exp(lnrhoAns);

	DETATCHSTATUSPTR(status);
	RETURN(status);

}
