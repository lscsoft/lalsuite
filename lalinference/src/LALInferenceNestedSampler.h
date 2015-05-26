/*
 *
 *  LALInferenceNestedSampler: Nested sampler written for LALInference
 *  LALInferenceNestedSampler.h:   main header file
 *
 *  Copyright (C) 2011 John Veitch
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
#ifndef LALInferenceNestedSampler_h
#define LALInferenceNestedSampler_h

#include <lal/LALInference.h>

#ifdef SWIG // SWIG interface directives
SWIGLAL(
	FUNCTION_POINTER(
			 LALInferenceNestedSamplingAlgorithm,
			 LALInferenceNestedSamplingOneStep
			 )
);
#endif



/**
 * \defgroup LALInferenceNestedSampler_h Header LALInferenceNestedSampler.h
 * \ingroup lalinference_general
 * \brief Nested sampler written for LALInference. Independent of model.
 *
 * Nested Sampling algorithm defined using the LALInference
 * Infrastructure. This code should be independent of choice
 * of model. Provided are a LALAlgorithm function and a
 * LALEvolveOneStepFunction which implement the evidence
 * calculation
 */
/*@{*/

/* logadd(a,b) = log(exp(a) + exp(b)) using Stirling's approximation */
/* double logadd(double a,double b); */

/**
 * Initialise the nested sampling algorithm by reading from the commandLine
 * and setting up algorithmParams
 */
void LALInferenceNestedSamplingAlgorithmInit(LALInferenceRunState *runState);

/**
 * NestedSamplingAlgorithm implements the nested sampling algorithm,
 * see e.g. Sivia "Data Analysis: A Bayesian Tutorial, 2nd edition
 */
void LALInferenceNestedSamplingAlgorithm(LALInferenceRunState *runState);

/** A single iteration of the NS algorithm */
INT4 LALInferenceNestedSamplingOneStep(LALInferenceRunState *runState);

/** Compute the autocorrelation length from the sampler at the current global iteration */
LALInferenceVariables *LALInferenceComputeAutoCorrelation(LALInferenceRunState *runState, UINT4 max_iterations, LALInferenceEvolveOneStepFunction evolve);

/** Perform one MCMC iteration on runState->currentParams. Return 1 if accepted or 0 if not */
UINT4 LALInferenceMCMCSamplePrior(LALInferenceRunState *runState);

/** Sample the prior N times, returns number of acceptances */
UINT4 LALInferenceMCMCSamplePriorNTimes(LALInferenceRunState *runState, UINT4 N);

/**
 * Sample the limited prior distribution using the MCMC method as usual, but
 * run a sub-chain of x iterations which doesn't check the likelihood bound.
 * x=LALInferenceGetVariable(runState->algorithmParams,"sloppyratio")
 */
INT4 LALInferenceNestedSamplingSloppySample(LALInferenceRunState *runState);

/* Cache wrapper around the sloppy sampler */
INT4 LALInferenceNestedSamplingCachedSampler(LALInferenceRunState *runState);

/**
 * Setup the live points by calling runState->initVariables on each of them
 * if it is specified. Otherwise clones runState->currentParams (legacy)
 */
void LALInferenceSetupLivePointsArray(LALInferenceRunState *runState);

/** Project the sample in params onto the eigenvectors given in eigenvectors. */
void LALInferenceProjectSampleOntoEigenvectors(LALInferenceVariables *params, gsl_matrix *eigenvectors, REAL8Vector **projection);

/*@}*/

#endif
