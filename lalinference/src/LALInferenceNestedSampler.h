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

/**
 * \file LALInferenceNestedSampler.h
 * \brief Nested sampler written for LALInference. Independent of model.
 */

#ifndef LALInferenceNestedSampler_h
#define LALInferenceNestedSampler_h


/** Nested Sampling algorithm defined using the LALInference
 * Infrastructure. This code should be independent of choice
 * of model. Provided are a LALAlgorithm function and a 
 * LALEvolveOneStepFunction which implement the evidence
 * calculation
 */

#include <lal/LALInference.h>

/* logadd(a,b) = log(exp(a) + exp(b)) using Stirling's approximation */
/* double logadd(double a,double b); */

/* Append sample to an array maintained elsewhere */
void LALInferenceLogSampleToArray(LALInferenceRunState *state, LALInferenceVariables *vars);

/** NestedSamplingAlgorithm implements the nested sampling algorithm,
 see e.g. Sivia "Data Analysis: A Bayesian Tutorial, 2nd edition */
void LALInferenceNestedSamplingAlgorithm(LALInferenceRunState *runState);

/** A differential evolution proposal for nested sampling algorithm */
void LALInferenceProposalDifferentialEvolution(LALInferenceRunState *runState, LALInferenceVariables *parameter);

/** Calculate covariance matrix from a collection of live points */
void LALInferenceNScalcCVM(gsl_matrix **cvm, LALInferenceVariables **Live, UINT4 Nlive);
/** This should be moved */
/* double logadd(double a,double b); */

/** A single iteration of the NS algorithm */
void LALInferenceNestedSamplingOneStep(LALInferenceRunState *runState);

/** Compute the autocorrelation length from the sampler at the current global iteration */
LALInferenceVariables *LALInferenceComputeAutoCorrelation(LALInferenceRunState *runState, UINT4 max_iterations, LALInferenceEvolveOneStepFunction *evolve);

/** Perform one MCMC iteration on runState->currentParams. Return 1 if accepted or 0 if not */
UINT4 LALInferenceMCMCSamplePrior(LALInferenceRunState *runState);

/** Sample the prior N times, returns number of acceptances */
UINT4 LALInferenceMCMCSamplePriorNTimes(LALInferenceRunState *runState, UINT4 N);

/** Sample the limited prior distribution using the MCMC method as usual, but
   run a sub-chain of x iterations which doesn't check the likelihood bound.
   x=LALInferenceGetVariable(runState->algorithmParams,"sloppyratio")
*/
void LALInferenceNestedSamplingSloppySample(LALInferenceRunState *runState);

/* REAL8 mean(REAL8 *array,int N); */
REAL8 LALInferenceNSSample_logt(int Nlive,gsl_rng *RNG);
REAL8 LALInferenceAngularDistance(REAL8 a1, REAL8 a2);
REAL8 LALInferenceAngularVariance(LALInferenceVariables **list,const char *pname, int N);

/* Functions for proposal distributions */
void LALInferenceProposalNS(LALInferenceRunState *runState, LALInferenceVariables *parameter);
void LALInferenceProposalMultiStudentT(LALInferenceRunState *runState, LALInferenceVariables *parameter);
//Declared in LALInferencePrior.h instead:
//void LALInferenceCyclicReflectiveBound(LALInferenceVariables *parameter, LALInferenceVariables *priorArgs);
INT4 LALInferenceReflectDetPlane(
							 LALInferenceRunState *state,
							 LALInferenceVariables *parameter
								 );
void LALInferenceRotateSky(
						   LALInferenceRunState *state,
						   LALInferenceVariables *parameter
						   );

/*void crossProduct(REAL8 out[3],REAL8 x[3],REAL8 y[3]);
void CartesianToSkyPos(REAL8 pos[3],REAL8 *longitude, REAL8 *latitude);
void GetCartesianPos(REAL8 vec[3],REAL8 longitude, REAL8 latitude);*/



/* Returns a REAL4vector drawn from the multivariate normal distribution
 with covariance matrix matrix, mean zero, of dimension dim. randParam is
 an already initialised RandomParams structure */

void XLALMultiNormalDeviates(REAL4Vector *vector,
							 gsl_matrix *matrix,
							 UINT4 dim,
							 RandomParams *randParam);

/* Returns a REAL4vector drawn from the multivariate student-t distribution
 with n degrees of freedom, with covariance matrix matrix, mean zero,
 of dimension dim. randParam is an already initialised RandomParams structure */

void
XLALMultiStudentDeviates(REAL4Vector  *vector,
						 gsl_matrix   *matrix,
						 UINT4         dim,
						 UINT4         n,
						 RandomParams *randParam);


/* Check that the gsl_matrix is positive definite. dim = number of dimensions */
UINT4 LALInferenceCheckPositiveDefinite(gsl_matrix *matrix, UINT4 dim);
void LALInferenceSetupLivePointsArray(LALInferenceRunState *runState);

/** Setup a k-D tree from the current set of nested sampling live points for use
    as a proposal distribution. */
void LALInferenceSetupkDTreeNSLivePoints( LALInferenceRunState *runState );

#endif

