/*
 *
 *  LALInferenceMCMCSampler:    Markov-Chain Monte Carlo sampler for LALInference        
 *  LALInferenceMCMCSampler.h:  main header file
 *
 *  Copyright (C) 2011 Vivien Raymond, Ben Farr, Will Farr, Ilya Mandel, Christian Roever, Marc van der Sluys and John Veitch
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
 * \file LALInferenceMCMCSampler.h
 * \brief Markov-Chain Monte Carlo sampler written for LALInference. Independent of model.
 * \ingroup LALInference
 * 
 * Markov-Chain Monte Carlo sampler incorporating parallel tempering using MPI and 
 * the possibility of adaptative jumps.
 * 
 * Provided are a LALAlgorithm function and a 
 * LALEvolveOneStepFunction which implement a single step forward.
 * 
 */


#include <lal/LALInference.h>

#define COVMATRIXNAME "covarianceMatrix"
#define UNCORRSAMPNAME "uncorrelatedSample"
#define SIGMAVECTORNAME "sigmaJump"

/** Implements the parallel tempered MCMC algorithm. Designes to use PTMCMCOneStep() as the runstate->evolve function */
void PTMCMCAlgorithm(struct tagLALInferenceRunState *runState);
/** Implements one MCMC step forward, updating the sigma values for the jump proposals if required.*/
void PTMCMCOneStep(LALInferenceRunState *runState);
