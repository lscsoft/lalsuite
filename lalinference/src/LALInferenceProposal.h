/*
 * 
 *  LALInferenceProposal:             LALInference Proposal functions       
 *  include/LALInferenceProposal.h:   main header file
 *
 *  Copyright (C) 2009 Ilya Mandel, Vivien Raymond, Christian Roever, Marc van der Sluys and John Veitch
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
 * \file LALInferenceProposal.h
 * \brief LALInference Proposal functions for use in algorithms.
 * \ingroup LALInference
 * 
 * LALInference is a Bayesian analysis toolkit for use with LAL. It contains
 * common requirements for Bayesian codes such as Likelihood functions, data
 * handling routines, MCMC and Nested Sampling algorithms and a template generation
 * interface to the LALInspiral package.
 * 
 * This file contains jump proposals for use with LALInference. Supported jump types
 * include:
 * covariance-matrix, adaptive, differential evolution, distance-inclination, etc
 * 
 * 
 */




/* PTMCMC algorithm defined using the LALInference
 * Infrastructure. This code should be independent of choice
 * of model. Provided are a LALAlgorithm function and a 
 * LALEvolveOneStepFunction which implement the evidence
 * calculation
 */

#ifndef LALInferenceProposal_h
#define LALInferenceProposal_h

#include <lal/LALInference.h>

#define COVMATRIXNAME "covarianceMatrix"
#define UNCORRSAMPNAME "uncorrelatedSample"
#define SIGMAVECTORNAME "sigmaJump"

/** Wrapper for PTMCMCLALProposal that sets certains variables used in the MCMC code that aren't in the NS code */
void NSWrapMCMCLALProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);

void PTMCMCLALProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);
void PTMCMCLALBlockProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);
void PTMCMCLALSingleProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);

/* Makes correlated jumps if the square root (Cholesky decomp) of a
   covariance matrix has been specified; if not, falls back to PTMCMCLALBlockProposal. */
void PTMCMCLALBlockCorrelatedProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);
void PTMCMCLALSingleAdaptProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);

/* Flip inclination about observing plane (iota -> Pi - iota) */
void PTMCMCLALInferenceInclinationFlip(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);

/* Rotate one or both spins about L. */
void PTMCMCLALInferenceRotateSpins(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);

/* Increment orbital phase by Pi */
void PTMCMCLALInferenceOrbitalPhaseJump(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);

/* Choose a random covariance matrix eigenvector to jump along. */
void PTMCMCLALInferenceCovarianceEigenvectorJump(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);

/* Jump around by 0.01 radians in angle on the sky */
void PTMCMCLALInferenceSkyLocWanderJump(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);

void PTMCMCLALAdaptationProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);
void PTMCMCLALAdaptationSingleProposal(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);

void PTMCMCLALInferenceRotateSky(LALInferenceRunState *state,LALInferenceVariables *parameter);
INT4 PTMCMCLALInferenceReflectDetPlane(LALInferenceRunState *state,LALInferenceVariables *parameter);

/* Jump in iota and distance so that approximate magnitude remains
   constant in one of the detectors. */
void PTMCMCLALInferenceInclinationDistanceConstAmplitudeJump(LALInferenceRunState *state, LALInferenceVariables *proposedParams);

/* Differential evolution */
void PTMCMCLALInferenceDifferentialEvolutionFull(LALInferenceRunState *state, LALInferenceVariables *proposedParams);
void PTMCMCLALInferenceDifferentialEvolutionNames(LALInferenceRunState *state, LALInferenceVariables *proposedParams, const char *names[]);
void PTMCMCLALInferenceDifferentialEvolutionMasses(LALInferenceRunState *state, LALInferenceVariables *proposedParams);
void PTMCMCLALInferenceDifferentialEvolutionAmp(LALInferenceRunState *state, LALInferenceVariables *proposedParams);
void PTMCMCLALInferenceDifferentialEvolutionSpins(LALInferenceRunState *state, LALInferenceVariables *proposedParams);
void PTMCMCLALInferenceDifferentialEvolutionSky(LALInferenceRunState *state, LALInferenceVariables *proposedParams);

/*draws a value from the prior, uniformly in individual parameters used for jumps.*/
void PTMCMCLALInferenceDrawFromPrior(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);
/*draws a value from the prior, using Von Neumann rejection sampling.*/
void PTMCMCLALInferenceDrawUniformlyFromPrior(LALInferenceRunState *runState, LALInferenceVariables *proposedParams);
// Von Neumann rejection sampler for the prior !!
//void VNRPriorOneStep(LALInferenceRunState *runState);

#endif

