/* 
 *  LALInferencePriorTest.c: Testing the Nested Sampling routines in LALInferencePrior.c
 *
 *  Copyright (C) 2011 Ben Aylott, Ilya Mandel, Chiara Mingarelli, Vivien Raymond, Christian Roever, Marc van der Sluys, John Veitch, Will Vousden
 *
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

#include <stdio.h>
#include <lal/LALInferencePrior.h>
#include <lal/LALInference.h>
#include <lal/XLALError.h>
#include <math.h>
#include <assert.h>
#include "LALInferenceTest.h"

#define EPSILON 1e-13

int computePriorMassNormTest(void);
int LALInferenceRotateInitialPhaseTest(void);
int LALInferenceCyclicReflectiveBoundTest(void);
int LALInferenceDrawFromPriorTest(void);
int LALInferenceInspiralPriorTest(void);

//Old tests
REAL8 BasicUniformLALPrior(LALInferenceRunState *runState, LALInferenceVariables *params);
REAL8 ASinOmegaTPrior(LALInferenceRunState *runState, LALInferenceVariables *params);

int main(void)
{
	int failureCount = 0;

	failureCount += computePriorMassNormTest();
	printf("\n");
	failureCount += LALInferenceRotateInitialPhaseTest();
	printf("\n");
	failureCount += LALInferenceCyclicReflectiveBoundTest();
	printf("\n");
	failureCount += LALInferenceDrawFromPriorTest();
	printf("\n");
	failureCount += LALInferenceInspiralPriorTest();
	printf("\n");

	printf("Test results: %i failure(s).\n", failureCount);
	return failureCount;
}

int computePriorMassNormTest(void)
{
	TEST_HEADER();

	int errnum;
	double result;

	char massRatioName[VARNAME_MAX];
	REAL8 MMin;
	REAL8 MMax;
	REAL8 MTotMax;
	REAL8 McMin;
	REAL8 McMax;
	REAL8 massRatioMin;
	REAL8 massRatioMax;

	MMin = 1;
	MMax = 10;
	MTotMax = 100;
	McMin = 1;
	McMax = 2;
	massRatioMin = 1;
	massRatioMax = 10;
	XLAL_TRY(result = LALInferenceComputePriorMassNorm(MMin, MMax, MTotMax, McMin, McMax, massRatioMin, massRatioMax, NULL), errnum);
	if (!XLAL_IS_REAL8_FAIL_NAN(result) || errnum != XLAL_EFAULT)
		TEST_FAIL("Null reference check failed.");

	strcpy(massRatioName, "foo");
	XLAL_TRY(result = LALInferenceComputePriorMassNorm(MMin, MMax, MTotMax, McMin, McMax, massRatioMin, massRatioMax, massRatioName), errnum);
	if (!XLAL_IS_REAL8_FAIL_NAN(result) || errnum != XLAL_ENAME)
		TEST_FAIL("Invalid mass ratio name specified but appropriate error not generated.");

	strcpy(massRatioName, "asym_massratio");
	MMin = 1;
	MMax = -10;
	MTotMax = -100;
	McMin = -1;
	McMax = -2;
	massRatioMin = 1;
	massRatioMax = 10;
	XLAL_TRY(result = LALInferenceComputePriorMassNorm(MMin, MMax, MTotMax, McMin, McMax, massRatioMin, massRatioMax, massRatioName), errnum);
	if (!XLAL_IS_REAL8_FAIL_NAN(result) || errnum == XLAL_SUCCESS)
		TEST_FAIL("Unphysical masses given but appropriate error not generated.");

	TEST_FOOTER();
}

int LALInferenceRotateInitialPhaseTest(void)
{
	TEST_HEADER();

	int errnum;

	// A basic null reference check.
	XLAL_TRY(LALInferenceRotateInitialPhase(NULL), errnum);
	if (errnum != XLAL_EFAULT)
		TEST_FAIL("Null reference check failed.");

	// Construct a variable list containing phi0 and psi.
	REAL8 psi = 0;
	REAL8 phi0 = 0;
	REAL8 phi0_2 = 0;
	LALInferenceVariables *variables = XLALCalloc(1, sizeof(LALInferenceVariables));
	LALInferenceAddVariable(variables, "psi", &psi, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
	LALInferenceAddVariable(variables, "phi0", &phi0, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);

	// Check that if psi is in [0,2pi], phi0 remains unchanged.
	psi = LAL_PI;
	phi0 = 0;
	LALInferenceSetVariable(variables, "psi", &psi);
	LALInferenceSetVariable(variables, "phi0", &phi0);
	LALInferenceRotateInitialPhase(variables);
	phi0_2 = *(REAL8 *)LALInferenceGetVariable(variables, "phi0");
	if (!compareFloats(phi0_2, phi0, EPSILON))
		TEST_FAIL("Psi in [0,2pi] but phi0 has changed!");

	// Check that if psi is outside [0,2pi], phi0 is rotated.
	psi = -2 * LAL_TWOPI;
	phi0 = 0;
	LALInferenceSetVariable(variables, "psi", &psi);
	LALInferenceSetVariable(variables, "phi0", &phi0);
	LALInferenceRotateInitialPhase(variables);
	phi0_2 = *(REAL8 *)LALInferenceGetVariable(variables, "phi0");
	if (!compareFloats(phi0_2, phi0 + LAL_PI, EPSILON))
		TEST_FAIL("Psi outside [0,2pi] but phi0 not rotated by 2pi!");
	psi = 2 * LAL_TWOPI;
	phi0 = 0;
	LALInferenceSetVariable(variables, "psi", &psi);
	LALInferenceSetVariable(variables, "phi0", &phi0);
	LALInferenceRotateInitialPhase(variables);
	phi0_2 = *(REAL8 *)LALInferenceGetVariable(variables, "phi0");
	if (!compareFloats(phi0_2, phi0 + LAL_PI, EPSILON))
		TEST_FAIL("Psi outside [0,2pi] but phi0 not rotated by 2pi!");

	// Check boundary cases: if psi=0 or psi=2pi, phi0 shouldn't be rotated.
	psi = 0;
	phi0 = 0;
	LALInferenceSetVariable(variables, "psi", &psi);
	LALInferenceSetVariable(variables, "phi0", &phi0);
	LALInferenceRotateInitialPhase(variables);
	phi0_2 = *(REAL8 *)LALInferenceGetVariable(variables, "phi0");
	if (phi0 != phi0_2) // Should be exact.)
		TEST_FAIL("Psi on boundary of [0,2pi] but phi0 has changed!");
	psi = LAL_TWOPI;
	phi0 = 0;
	LALInferenceSetVariable(variables, "psi", &psi);
	LALInferenceSetVariable(variables, "phi0", &phi0);
	LALInferenceRotateInitialPhase(variables);
	phi0_2 = *(REAL8 *)LALInferenceGetVariable(variables, "phi0");
	if (phi0 != phi0_2) // Should be exact.)
		TEST_FAIL("Psi on boundary of [0,2pi] but phi0 has changed!");

	XLALFree(variables);
	TEST_FOOTER();
}

int LALInferenceCyclicReflectiveBoundTest(void)
{
	TEST_HEADER();

	int errnum;
	int outcome;

	const LALInferenceVariableType type = LALINFERENCE_REAL8_t;
	REAL8 a;
	REAL8 b;
	REAL8 a_2;
	REAL8 b_2;
	REAL8 a_min;
	REAL8 a_max;
	REAL8 a_delta;
	REAL8 b_min;
	REAL8 b_max;
	REAL8 b_delta;
	LALInferenceVariables *parameters = XLALCalloc(1, sizeof(LALInferenceVariables));
	LALInferenceVariables *priorArgs = XLALCalloc(1, sizeof(LALInferenceVariables));

	// A basic null reference check.
	outcome = 1;
	XLAL_TRY(LALInferenceCyclicReflectiveBound(NULL, priorArgs), errnum);
	outcome &= errnum == XLAL_EFAULT;
	XLAL_TRY(LALInferenceCyclicReflectiveBound(parameters, NULL), errnum);
	outcome &= errnum == XLAL_EFAULT;
	if (!outcome)
		TEST_FAIL("Null reference check failed.");
	
	// Check some (meaningful) minima/maxima.
	a_min = -LAL_PI;
	a_max = LAL_PI;
	a_delta = a_max - a_min;
	b_min = -1;
	b_max = 1;
	b_delta = b_max - b_min;
	LALInferenceRemoveMinMaxPrior(priorArgs, "a");
	LALInferenceRemoveMinMaxPrior(priorArgs, "b");
	LALInferenceAddMinMaxPrior(priorArgs, "a", &a_min, &a_max, type);
	LALInferenceAddMinMaxPrior(priorArgs, "b", &b_min, &b_max, type);

	// Variables within [min,max]: should remain unchanged.
	a = a_min + (a_max - a_min) / 2;
	b = b_min + (b_max - b_min) / 2;
	LALInferenceAddVariable(parameters, "a", &a, type, LALINFERENCE_PARAM_CIRCULAR);
	LALInferenceAddVariable(parameters, "b", &b, type, LALINFERENCE_PARAM_LINEAR);
	LALInferenceCyclicReflectiveBound(parameters, priorArgs);
	a_2 = *(REAL8 *)LALInferenceGetVariable(parameters, "a");
	b_2 = *(REAL8 *)LALInferenceGetVariable(parameters, "b");
	if (!compareFloats(a, a_2, EPSILON) || !compareFloats(b, b_2, EPSILON))
		TEST_FAIL("Values within bounds should remain unchanged.");

	// Boundary cases (circular): variables on [min, max] boundaries should be equal modulo period.
	outcome = 1;
	a = a_min;
	LALInferenceAddVariable(parameters, "a", &a, type, LALINFERENCE_PARAM_CIRCULAR);
	LALInferenceCyclicReflectiveBound(parameters, priorArgs);
	a_2 = *(REAL8 *)LALInferenceGetVariable(parameters, "a");
	outcome &= compareFloats(fmod(a - a_2, a_delta), 0, EPSILON);
	a = a_max;
	LALInferenceAddVariable(parameters, "a", &a, type, LALINFERENCE_PARAM_CIRCULAR);
	LALInferenceCyclicReflectiveBound(parameters, priorArgs);
	a_2 = *(REAL8 *)LALInferenceGetVariable(parameters, "a");
	outcome &= compareFloats(fmod(a - a_2, a_delta), 0, EPSILON);
	if (!outcome)
		TEST_FAIL("Circular boundary values should remain equal modulo their period.");

	// Boundary cases (linear): variables on [min, max] boundaries should be equal modulo period.
	outcome = 1;
	b = b_min;
	LALInferenceAddVariable(parameters, "b", &b, type, LALINFERENCE_PARAM_LINEAR);
	LALInferenceCyclicReflectiveBound(parameters, priorArgs);
	b_2 = *(REAL8 *)LALInferenceGetVariable(parameters, "b");
	outcome &= compareFloats(b, b_2, EPSILON);
	b = b_max;
	LALInferenceAddVariable(parameters, "b", &b, type, LALINFERENCE_PARAM_LINEAR);
	LALInferenceCyclicReflectiveBound(parameters, priorArgs);
	b_2 = *(REAL8 *)LALInferenceGetVariable(parameters, "b");
	outcome &= compareFloats(b, b_2, EPSILON);
	if (!outcome)
		TEST_FAIL("Linear boundary values should remain unchanged.");

	// Outside range (circular).
	outcome = 1;
	a = a_min - a_delta / 3;
	LALInferenceAddVariable(parameters, "a", &a, type, LALINFERENCE_PARAM_CIRCULAR);
	LALInferenceCyclicReflectiveBound(parameters, priorArgs);
	a_2 = *(REAL8 *)LALInferenceGetVariable(parameters, "a");
	outcome &= compareFloats(a_2, a_max - a_delta / 3, EPSILON);
	a = a_max + a_delta / 3;
	LALInferenceAddVariable(parameters, "a", &a, type, LALINFERENCE_PARAM_CIRCULAR);
	LALInferenceCyclicReflectiveBound(parameters, priorArgs);
	a_2 = *(REAL8 *)LALInferenceGetVariable(parameters, "a");
	outcome &= compareFloats(a_2, a_min + a_delta / 3, EPSILON);
	if (!outcome)
		TEST_FAIL("Circular values outside range should be correctly modded into range.");
	
	// Outside range (linear).
	outcome = 1;
	b = b_min - 10 * b_delta / 3;
	LALInferenceAddVariable(parameters, "b", &b, type, LALINFERENCE_PARAM_LINEAR);
	LALInferenceCyclicReflectiveBound(parameters, priorArgs);
	b_2 = *(REAL8 *)LALInferenceGetVariable(parameters, "b");
	outcome &= compareFloats(b_2, b_max - b_delta / 3, EPSILON);
	b = b_max + 7 * b_delta / 5;
	LALInferenceAddVariable(parameters, "b", &b, type, LALINFERENCE_PARAM_LINEAR);
	LALInferenceCyclicReflectiveBound(parameters, priorArgs);
	b_2 = *(REAL8 *)LALInferenceGetVariable(parameters, "b");
	outcome &= compareFloats(b_2, b_min + 2 * b_delta / 5, EPSILON);
	if (!outcome)
		TEST_FAIL("Linear values outside range should be correctly reflected into range.");

	XLALFree(parameters);
	XLALFree(priorArgs);
	TEST_FOOTER();
}

int LALInferenceDrawFromPriorTest(void)
{
	TEST_HEADER();

	int errnum;
	const char *name;
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(rng, 0);

	LALInferenceVariables *output = XLALCalloc(1, sizeof(LALInferenceVariables));
	LALInferenceVariables *priorArgs = XLALCalloc(1, sizeof(LALInferenceVariables));

	// Null reference checks.
	int outcome = 1;
	XLAL_TRY(LALInferenceDrawFromPrior(NULL, priorArgs, rng), errnum);
	outcome &= errnum == XLAL_EFAULT;
	XLAL_TRY(LALInferenceDrawFromPrior(output, NULL, rng), errnum);
	outcome &= errnum == XLAL_EFAULT;
	XLAL_TRY(LALInferenceDrawFromPrior(output, priorArgs, NULL), errnum);
	outcome &= errnum == XLAL_EFAULT;
	if (!outcome)
		TEST_FAIL("Null reference check failed.");

	int i;
	const char *varyName=NULL;
	char caseTag[VARNAME_MAX];
	LALInferenceVariableType type = LALINFERENCE_REAL8_t;
	LALInferenceParamVaryType vary=-1;
	for (i = 0; i < 2; i++)
	{
		switch (i)
		{
			case 0:
				vary = LALINFERENCE_PARAM_LINEAR;
				varyName = "linear";
				break;
			case 1:
				vary = LALINFERENCE_PARAM_CIRCULAR;
				varyName = "circular";
				break;
		}
		sprintf(caseTag, "[%s] ", varyName);

		// Try and generate some normally distributed variables for various mu and sigma.
		REAL8 gaussian = 0;
		REAL8 mu;
		REAL8 sigma;
		name = "gaussian";
		LALInferenceAddVariable(output, name, &gaussian, type, vary);

		// Zero standard deviation; should always equal mean.
		mu = -50;
		sigma = 0;
		LALInferenceRemoveGaussianPrior(priorArgs, name);
		LALInferenceAddGaussianPrior(priorArgs, name, &mu, &sigma, type);
		XLAL_TRY(LALInferenceDrawFromPrior(output, priorArgs, rng), errnum);
		gaussian = *(REAL8*)LALInferenceGetVariable(output, name);
		if (errnum != XLAL_SUCCESS)
		{
			TEST_FAIL("%sFailed to generate Gaussian variable; XLAL error: %s.", caseTag, XLALErrorString(errnum));
		}
		else if (!compareFloats(gaussian, mu, EPSILON))
		{
			TEST_FAIL("%sGaussian variable with zero standard deviation did not match the mean; X = %f, mu = %f.", caseTag, gaussian, mu);
		}

		LALInferenceRemoveVariable(output, name);
		LALInferenceRemoveGaussianPrior(priorArgs, name);

		// Try a uniform variable!
		REAL8 uniform = 0;
		REAL8 min;
		REAL8 max;
		name = "uniform";
		LALInferenceAddVariable(output, name, &uniform, type, vary);

		min = -1;
		max = 1;
		LALInferenceRemoveMinMaxPrior(priorArgs, name);
		LALInferenceAddMinMaxPrior(priorArgs, name, &min, &max, type);
		XLAL_TRY(LALInferenceDrawFromPrior(output, priorArgs, rng), errnum);
		if (errnum != XLAL_SUCCESS)
			TEST_FAIL("%sFailed to generate uniform variable; XLAL error: %s.", caseTag, XLALErrorString(errnum));

		LALInferenceRemoveVariable(output, name);
		LALInferenceRemoveMinMaxPrior(priorArgs, name);

		// Try a correlated variable!
		REAL8 correlated = 0;
		UINT4 idx = 0;
		name = "correlated";
		gsl_matrix *covariance = gsl_matrix_calloc(3, 3);
		LALInferenceAddVariable(output, name, &correlated, type, vary);
		LALInferenceRemoveCorrelatedPrior(priorArgs, name);
		LALInferenceAddCorrelatedPrior(priorArgs, name, &covariance, &idx);

		// See what happens when we try to generate correlated values from a non-positive-definite
		// covariance matrix.
		gsl_matrix_set(covariance, 0, 0, -1);
		XLAL_TRY(LALInferenceDrawFromPrior(output, priorArgs, rng), errnum);
		if (errnum == XLAL_SUCCESS)
			TEST_FAIL("%sNon-positive-definite covariance matrix was not rejected.", caseTag);

		// Now try a positive-semi-definite matrix; this should be accepted (need only update matrix, not add it afresh).
		gsl_matrix_set(covariance, 0, 0, 1);
		XLAL_TRY(LALInferenceDrawFromPrior(output, priorArgs, rng), errnum);
		if (errnum != XLAL_SUCCESS)
			TEST_FAIL("%sCould not generate correlated variable from positive-semi-definite matrix; XLAL error: %s.", caseTag, XLALErrorString(errnum));
		
		// Try a legitimate positive-definite covariance matrix.
		gsl_matrix_set(covariance, 0, 0, 2);
		gsl_matrix_set(covariance, 0, 1, 1);
		gsl_matrix_set(covariance, 0, 2, 0);
		gsl_matrix_set(covariance, 1, 0, 1);
		gsl_matrix_set(covariance, 1, 1, 5);
		gsl_matrix_set(covariance, 1, 2, 1);
		gsl_matrix_set(covariance, 2, 0, 0);
		gsl_matrix_set(covariance, 2, 1, 1);
		gsl_matrix_set(covariance, 2, 2, 1);
		XLAL_TRY(LALInferenceDrawFromPrior(output, priorArgs, rng), errnum);
		if (errnum != XLAL_SUCCESS)
			TEST_FAIL("%sCould not generate correlated variable from positive-definite matrix; XLAL error: %s.", caseTag, XLALErrorString(errnum));

		LALInferenceRemoveVariable(output, name);
		LALInferenceRemoveCorrelatedPrior(priorArgs, name);

		gsl_matrix_free(covariance);
		LALInferenceRemoveVariable(output, "gaussian");
		LALInferenceRemoveVariable(output, "uniform");
		LALInferenceRemoveVariable(output, "correlated");
	}

	XLALFree(output);
	XLALFree(priorArgs);
	TEST_FOOTER();
}

int LALInferenceInspiralPriorTest(void)
{
	TEST_HEADER();

	int errnum;

	REAL8 result;
	LALInferenceRunState *runState = XLALCalloc(1, sizeof(LALInferenceRunState));
	LALInferenceVariables *params = XLALCalloc(1, sizeof(LALInferenceVariables));
	LALInferenceVariables *priorArgs = XLALCalloc(1, sizeof(LALInferenceVariables));

	// Standard null reference check.
	int failed = 1;
	runState->priorArgs = NULL;
	XLAL_TRY(result = LALInferenceInspiralPrior(runState, params), errnum);
	failed &= !XLAL_IS_REAL8_FAIL_NAN(result) || errnum != XLAL_EFAULT;
	runState->priorArgs = priorArgs;
	XLAL_TRY(result = LALInferenceInspiralPrior(NULL, params), errnum);
	failed &= !XLAL_IS_REAL8_FAIL_NAN(result) || errnum != XLAL_EFAULT;
	XLAL_TRY(result = LALInferenceInspiralPrior(runState, NULL), errnum);
	failed &= !XLAL_IS_REAL8_FAIL_NAN(result) || errnum != XLAL_EFAULT;
	if (failed)
		TEST_FAIL("Null reference check failed.");

	// Set up parameters.
	REAL8 value = 0;
	LALInferenceAddVariable(params, "logdistance", &value, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	LALInferenceAddVariable(params, "distance", &value, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	LALInferenceAddVariable(params, "inclination", &value, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
	LALInferenceAddVariable(params, "rightascension", &value, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
	LALInferenceAddVariable(params, "declination", &value, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
	LALInferenceAddVariable(params, "theta_spin1", &value, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
	LALInferenceAddVariable(params, "theta_spin2", &value, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_CIRCULAR);
	LALInferenceAddVariable(params, "logmc", &value, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	LALInferenceAddVariable(params, "chirpmass", &value, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	/*LALInferenceAddVariable(params, "asym_massratio", &value, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);*/
	LALInferenceAddVariable(params, "massratio", &value, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_LINEAR);
	
	REAL8 min, max;
	min = 1.0; max = 100.0;
	LALInferenceAddMinMaxPrior(priorArgs, "distance", &min, &max, LALINFERENCE_REAL8_t);
	min = log(min); max = log(max);
	LALInferenceAddMinMaxPrior(priorArgs, "logdistance", &min, &max, LALINFERENCE_REAL8_t);
	min = 1.0; max = 20.5;
	LALInferenceAddMinMaxPrior(priorArgs, "chirpmass", &min, &max, LALINFERENCE_REAL8_t);
	min = log(min); max = log(max);
	LALInferenceAddMinMaxPrior(priorArgs, "logmc", &min, &max, LALINFERENCE_REAL8_t);
	min = 1.0; max = 30.0;
	LALInferenceAddMinMaxPrior(priorArgs, "component", &min, &max, LALINFERENCE_REAL8_t);
	/*max *= 2;*/
	/*LALInferenceAddVariable(priorArgs, "MTotMax", &max, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);*/
	min = -LAL_PI; max = LAL_PI;
	LALInferenceAddMinMaxPrior(priorArgs, "inclination", &min, &max, LALINFERENCE_REAL8_t);
	min = 0; max = LAL_TWOPI;
	LALInferenceAddMinMaxPrior(priorArgs, "rightascension", &min, &max, LALINFERENCE_REAL8_t);
	min = -LAL_PI / 2.0; max = LAL_PI / 2.0;
	LALInferenceAddMinMaxPrior(priorArgs, "declination", &min, &max, LALINFERENCE_REAL8_t);
	min = -LAL_PI / 2.0; max = LAL_PI / 2.0;
	LALInferenceAddMinMaxPrior(priorArgs, "theta_spin1", &min, &max, LALINFERENCE_REAL8_t);
	min = -LAL_PI / 2.0; max = LAL_PI / 2.0;
	LALInferenceAddMinMaxPrior(priorArgs, "theta_spin2", &min, &max, LALINFERENCE_REAL8_t);
	min = 0.01; max = 0.25;
	LALInferenceAddMinMaxPrior(priorArgs, "massratio", &min, &max, LALINFERENCE_REAL8_t);

	// Pick a random point in the non-zero region of the parameter space.
	gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(rng, 0);
	LALInferenceDrawFromPrior(params, priorArgs, rng);

	// Check that we get a finite log prior.
	XLAL_TRY(result = LALInferenceInspiralPrior(runState, params), errnum);
	if (XLAL_IS_REAL8_FAIL_NAN(result) || errnum != XLAL_SUCCESS)
	{
		TEST_FAIL("Could not generate inspiral prior; XLAL error: %s", XLALErrorString(errnum));
	}
	else if (result == -DBL_MAX)
	{
		TEST_FAIL("Parameter configuration within specified min/max bounds for each parameter gave zero prior.");
	}

	// Now set a parameter outside its bounds and see what happens.
	LALInferenceGetMinMaxPrior(priorArgs, "distance", &min, &max);
	value = max + (max - min) / 2;
	LALInferenceSetVariable(params, "distance", &value);
	XLAL_TRY(result = LALInferenceInspiralPrior(runState, params), errnum);
	if (XLAL_IS_REAL8_FAIL_NAN(result) || errnum != XLAL_SUCCESS)
	{
		TEST_FAIL("Could not generate inspiral prior; XLAL error: %s", XLALErrorString(errnum));
	}
	else if (result != -DBL_MAX)
	{
		TEST_FAIL("Distance %f is outside [%f,%f] but prior is non-zero.", value, min, max);
	}

	// Try another configuration; this time set m1 and m2 such that one is *outside* its bounds,
	// but the chirp mass and symmetric mass ratio are still OK; this should be picked up and a
	// zero prior returned.
	LALInferenceDrawFromPrior(params, priorArgs, rng);
	LALInferenceGetMinMaxPrior(priorArgs, "component", &min, &max);
	REAL8 m2 = 0.5;
	REAL8 m1 = 3.82;
	REAL8 eta = m1 * m2 / pow(m1 + m2, 2);
	LALInferenceSetVariable(params, "massratio", &eta);
	REAL8 Mc = pow(m1 * m2, 3.0 / 5.0) / pow(m1 + m2, 1.0 / 5.0);
	LALInferenceSetVariable(params, "chirpmass", &Mc);
	REAL8 logMc = log(Mc);
	LALInferenceSetVariable(params, "logmc", &logMc);
	XLAL_TRY(result = LALInferenceInspiralPrior(runState, params), errnum);
	if (XLAL_IS_REAL8_FAIL_NAN(result) || errnum != XLAL_SUCCESS)
	{
		TEST_FAIL("Could not generate inspiral prior; XLAL error: %s", XLALErrorString(errnum));
	}
	else if (result != -DBL_MAX)
	{
		TEST_FAIL("Mass ratio %f and chirp mass %f define masses outside bounds [%f,%f], but prior is non-zero.", eta, Mc, min, max);
	}

	TEST_FOOTER();
}

/******************************************
 * 
 * Old tests
 * 
 ******************************************/

//Test LALInferencePriorFunction
REAL8 BasicUniformLALPrior(LALInferenceRunState *runState, LALInferenceVariables *params)
/****************************************/
/* Returns unnormalized (!),            */
/* logarithmic (!) prior density.      	*/
/****************************************/
/* Assumes the following parameters	*/
/* exist (e.g., for TaylorT1):		*/
/* chirpmass, massratio, inclination,	*/
/* phase, time, rightascension,		*/
/* desclination, polarisation, distance.*/
/* Prior is flat if within range	*/
/****************************************/
{
  (void) runState; /* avoid warning about unused parameter */
  REAL8 eta, iota, phi, ra, dec, psi;
  REAL8 logdensity;
  
  // UNUSED!!: REAL8 mc   = *(REAL8*) LALInferenceGetVariable(params, "chirpmass");		/* solar masses*/
  eta  = *(REAL8*) LALInferenceGetVariable(params, "massratio");		/* dim-less    */
  iota = *(REAL8*) LALInferenceGetVariable(params, "inclination");		/* radian      */
  // UNUSED!!: REAL8 tc   = *(REAL8*) LALInferenceGetVariable(params, "time");			/* GPS seconds */
  phi  = *(REAL8*) LALInferenceGetVariable(params, "phase");		/* radian      */
  ra   = *(REAL8*) LALInferenceGetVariable(params, "rightascension");	/* radian      */
  dec  = *(REAL8*) LALInferenceGetVariable(params, "declination");		/* radian      */
  psi  = *(REAL8*) LALInferenceGetVariable(params, "polarisation"); 	/* radian      */
  // UNUSED!!: REAL8  dist = *(REAL8*) LALInferenceGetVariable(params, "distance");		/* Mpc         */

  if(eta>0.0 && eta<=0.25 && iota>=0.0 && iota<=LAL_PI && phi>=0.0 && phi<=LAL_TWOPI 
     && ra>=0.0 && ra<=LAL_TWOPI && dec>=-LAL_PI_2 && dec<=LAL_PI_2 && psi>=0.0 && psi<=LAL_PI)	
    logdensity = 0.0;
  else
    logdensity = -HUGE_VAL;
  //TODO: should be properly normalized; pass in range via priorArgs?	

  return(logdensity);
}

//Test LALPriorFunction
REAL8 ASinOmegaTPrior(LALInferenceRunState *runState, LALInferenceVariables *params)
/****************************************/
/* Prior for two-parameter				*/
/* waveform family ASinOmegaT			*/
/* Assumes the following parameters		*/
/* exist:	A, Omega					*/
/* Prior is flat if within range		*/
/****************************************/
{
  (void) runState; /* avoid warning about unused parameter */
  REAL8 A, Omega;
  REAL8 logdensity;
  
  A     = *(REAL8*) LALInferenceGetVariable(params, "A");				/* dim-less	   */
  Omega = *(REAL8*) LALInferenceGetVariable(params, "Omega");			/* rad/sec     */
  
  if ((A>0.0) & (Omega>0))
    logdensity = 0.0;
  else
    logdensity = -HUGE_VAL;

  return logdensity;
}
