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
 * \ingroup lalapps_inspiral
 * \brief Markov-Chain Monte Carlo sampler written for LALInference. Independent of model.
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
void mcmc_step(LALInferenceRunState *runState, LALInferenceThreadState *thread);

/* MPI communications */
typedef enum {
    PT_COM,          /** Parallel tempering communications */
    LADDER_UPDATE_COM,    /** Update positions across the ladder */
    RUN_PHASE_COM,   /** runPhase passing */
    RUN_COMPLETE       /** Run complete */
} LALInferenceMPIcomm;

/* Standard parallel temperature swap proposal function */
void LALInferencePTswap(LALInferenceRunState *runState, INT4 i, FILE *swapfile);

/* Metropolis-coupled MCMC swap proposal, when the likelihood is not identical between chains */
UINT4 LALInferenceMCMCMCswap(LALInferenceRunState *runState, REAL8 *ladder, INT4 i, FILE *swapfile);

/* Functions for controlling adaptation */
void acknowledgePhase(LALInferenceRunState *runState);
void LALInferenceAdaptation(LALInferenceThreadState *thread, INT4 cycle);
void LALInferenceAdaptationRestart(LALInferenceThreadState *thread, INT4 cycle);
REAL8 LALInferenceAdaptationEnvelope(INT4 cycle, INT4 start, INT4 length, INT4 tau, INT4 reset);
void LALInferenceShutdownLadder(void);
void LALInferenceFlushPTswap(void);
void LALInferenceLadderUpdate(LALInferenceRunState *runState, INT4 sourceChainFlag, INT4 cycle);

/* Data IO routines */
FILE** LALInferencePrintPTMCMCHeadersOrResume(LALInferenceRunState *runState);
void LALInferencePrintPTMCMCHeaderFiles(LALInferenceRunState *runState, FILE **files);
void LALInferencePrintPTMCMCInjectionSample(LALInferenceRunState *runState);
void LALInferenceDataDump(LALInferenceIFOData *data, LALInferenceModel *model);
void LALInferencePrintMCMCSample(LALInferenceThreadState *thread, LALInferenceIFOData *data, INT4 iteration, REAL8 timestamp, FILE *threadoutput);

/** Reads final parameter values from the given output file, and
    stores them in the current params to try to continue the run. */
INT4 LALInferenceMCMCResumeRead(LALInferenceThreadState *thread, FILE *resumeFile);

void init_mpi_randomstate(LALInferenceRunState *run_state);
