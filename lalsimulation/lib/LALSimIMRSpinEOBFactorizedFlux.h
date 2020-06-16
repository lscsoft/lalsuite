/*
 * LALSimIMRSpinEOBFactorizedFlux.h
 */

extern int UsePrec;

REAL8 XLALInspiralSpinFactorizedFlux (REAL8Vector * values,
					     EOBNonQCCoeffs * nqcCoeffs,
					     const REAL8 omega,
					     SpinEOBParams * ak,
					     const REAL8 H,
					     const UINT4 lMax,
					     const UINT4
					     SpinAlignedEOBversion);