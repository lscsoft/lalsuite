/*
 *  LALInferenceMCMC.c:  Bayesian Followup function testing site
 *
 *  Copyright (C) 2011 Ilya Mandel, Vivien Raymond, Christian Roever,
 *  Marc van der Sluys, John Veitch and Will M. Farr
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
#include <lal/Date.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALInference.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/StringInput.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/TimeSeries.h>
#include "LALInferenceEnsembleSampler.h"
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferenceProposal.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceReadData.h>
#include <lal/LALInferenceInit.h>
#include <lalapps.h>

#include <mpi.h>


LALInferenceRunState *initialize(ProcessParamsTable *commandLine);
void initializeMCMC(LALInferenceRunState *runState);


/* Read samples from file as initial state of ensemble */
void LALInferenceInitEnsemble(LALInferenceRunState *state);
void LALInferenceInitEnsemble(LALInferenceRunState *state) {
    LALInferenceVariableItem *item;
    LALInferenceIFOData *headData;
    ProcessParamsTable *ppt;
    REAL8 *sampleArray = NULL;
    UINT4 i=0;

    INT4 walker;
    MPI_Comm_rank(MPI_COMM_WORLD, &walker);

    INT4 ndim = LALInferenceGetVariableDimensionNonFixed(state->currentParams);
    INT4 *col_order = XLALMalloc(ndim *sizeof(INT4));

    ppt = LALInferenceGetProcParamVal(state->commandLine, "--init-samples");
    if (ppt) {
        if (walker == 0) {
            char *infile = ppt->value;
            FILE *input = fopen(infile, "r");

            char params[128][128];
            UINT4 ncols, nsamps;

            /* Parse parameter names */
            LALInferenceReadAsciiHeader(input, params, &ncols);

            LALInferenceVariables *backward_params = XLALCalloc(1, sizeof(LALInferenceVariables));

            /* Only cluster parameters that are being sampled */
            UINT4 nvalid_cols=0, j=0;
            UINT4 *valid_cols = XLALMalloc(ncols * sizeof(UINT4));
            for (j=0; j<ncols; j++)
                valid_cols[j] = 0;

            UINT4 logl_idx;
            for (j=0; j<ncols; j++) {
                char* internal_param_name = XLALMalloc(512*sizeof(char));
                LALInferenceTranslateExternalToInternalParamName(internal_param_name, params[j]);

                i=0;
                for (item = state->currentParams->head; item; item = item->next) {
                    if (LALInferenceCheckVariableNonFixed(state->currentParams, item->name)) {
                        if (!strcmp(item->name, internal_param_name)) {
                            col_order[i] = nvalid_cols;
                            nvalid_cols++;
                            valid_cols[j] = 1;
                            LALInferenceAddVariable(backward_params, item->name, item->value, item->type, item->vary);
                            break;
                        }
                        i++;
                    }
                }
            }

            /* Double check dimensions */
            if (nvalid_cols != ndim) {
                fprintf(stderr, "Inconsistent dimensions for starting state!\n");
                fprintf(stderr, "Sampling in %i dimensions, %i read from file!\n", ndim, nvalid_cols);
                exit(1);
            }

            /* LALInferenceAddVariable() builds the array backwards, so reverse it. */
            LALInferenceVariables *input_params = XLALCalloc(1, sizeof(LALInferenceVariables));

            for (item = backward_params->head; item; item = item->next)
                LALInferenceAddVariable(input_params, item->name, item->value, item->type, item->vary);

            sampleArray = LALInferenceParseDelimitedAscii(input, ncols, valid_cols, &nsamps);

            /* Choose unique samples to intialize ensemble */
            gsl_ran_shuffle(state->GSLrandom, sampleArray, nsamps, ndim*sizeof(REAL8));

            /* Cleanup */
            LALInferenceClearVariables(backward_params);
            XLALFree(backward_params);
        }

        /* Broadcast parameter order */
        MPI_Bcast(col_order, ndim, MPI_INT, 0, MPI_COMM_WORLD);

        REAL8Vector *parameters = XLALCreateREAL8Vector(ndim);
        MPI_Scatter(sampleArray, ndim, MPI_DOUBLE, parameters->data, ndim, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        REAL8Vector *reordered_parameters = XLALCreateREAL8Vector(ndim);
        for (i=0; i<ndim; i++)
            reordered_parameters->data[i] = parameters->data[col_order[i]];

        LALInferenceCopyArrayToVariables(reordered_parameters, state->currentParams);

        XLALDestroyREAL8Vector(parameters);
        XLALDestroyREAL8Vector(reordered_parameters);
        if (walker == 0)
            XLALFree(sampleArray);
    } else {
        LALInferenceDrawApproxPrior(state, state->currentParams);
        while (state->prior(state, state->currentParams) <= -DBL_MAX)
            LALInferenceDrawApproxPrior(state, state->currentParams);
    }

    /* Determine null loglikelihood that will be subtracted from printed likelihoods */
    REAL8 null_likelihood = 0.0;
    if (state->likelihood==&LALInferenceUndecomposedFreqDomainLogLikelihood ||
            state->likelihood==&LALInferenceFreqDomainLogLikelihood){
        null_likelihood = LALInferenceNullLogLikelihood(state->data);

    /* If no simple null likelihood method exists, scale signal into nothingness */
    } else if (state->likelihood == &LALInferenceFreqDomainStudentTLogLikelihood || 
                (state->likelihood == &LALInferenceMarginalisedTimeLogLikelihood &&
                (!LALInferenceGetProcParamVal(state->commandLine, "--malmquistPrior") ||
                !(LALInferenceGetProcParamVal(state->commandLine,"--psdFit") ||
                  LALInferenceGetProcParamVal(state->commandLine,"--glitchFit"))))) {
        headData = state->data;
        REAL8 d = *(REAL8 *)LALInferenceGetVariable(state->currentParams, "distance");
        REAL8 bigD = INFINITY;

        /* Don't store to cache, since distance scaling won't work */
        LALSimInspiralWaveformCache *cache = headData->waveformCache;
        while (headData != NULL) {
            headData->waveformCache = NULL;
            headData = headData->next;
        }
        headData = state->data;

        LALInferenceSetVariable(state->currentParams, "distance", &bigD);
        null_likelihood = state->likelihood(state->currentParams, state->data, state->templt);

        /* Restore cache to data structure */
        while (headData != NULL) {
            headData->nullloglikelihood = headData->loglikelihood;
            headData->waveformCache = cache;
            headData = headData->next;
        }

        /* Replace finite distance */
        LALInferenceSetVariable(state->currentParams, "distance", &d);
    }
    LALInferenceAddVariable(state->proposalArgs, "nullLikelihood", &null_likelihood,
                            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);

    /* Initialize starting likelihood and prior */
    state->currentLikelihood = state->likelihood(state->currentParams,
                                                 state->data, state->templt);
    headData = state->data;
    while (headData != NULL) {
        headData->acceptedloglikelihood = headData->loglikelihood;
        headData = headData->next;
    }

    state->currentPrior = state->prior(state, state->currentParams);
}

LALInferenceRunState *initialize(ProcessParamsTable *commandLine)
/* calls the "ReadData()" function to gather data & PSD from files, */
/* and initializes other variables accordingly.                     */
{
    LALInferenceRunState *irs=NULL;
    LALInferenceIFOData *ifoPtr, *ifoListStart;
  
    /* read data from files: */
    irs = XLALCalloc(1, sizeof(LALInferenceRunState));
    irs->commandLine=commandLine;
    LALInferenceCheckOptionsConsistency(commandLine);
    irs->data = LALInferenceReadData(commandLine);
  
    if (irs->data != NULL) {
        LALInferenceInjectInspiralSignal(irs->data,commandLine);
  
        ifoPtr = irs->data;
        ifoListStart = irs->data;
        while (ifoPtr != NULL) {
            /*If two IFOs have the same sampling rate, they should have the same timeModelh*,
              freqModelh*, and modelParams variables to avoid excess computation
              in model waveform generation in the future*/
            LALInferenceIFOData * ifoPtrCompare=ifoListStart;
            int foundIFOwithSameSampleRate=0;
            while (ifoPtrCompare != NULL && ifoPtrCompare!=ifoPtr) {
                if(ifoPtrCompare->timeData->deltaT == ifoPtr->timeData->deltaT){
                    ifoPtr->timeModelhPlus=ifoPtrCompare->timeModelhPlus;
                    ifoPtr->freqModelhPlus=ifoPtrCompare->freqModelhPlus;
                    ifoPtr->timeModelhCross=ifoPtrCompare->timeModelhCross;
                    ifoPtr->freqModelhCross=ifoPtrCompare->freqModelhCross;
                    ifoPtr->modelParams=ifoPtrCompare->modelParams;
                    foundIFOwithSameSampleRate=1;
                    break;
                }
                ifoPtrCompare = ifoPtrCompare->next;
            }
            if(!foundIFOwithSameSampleRate){
                ifoPtr->timeModelhPlus  = XLALCreateREAL8TimeSeries("timeModelhPlus",
                                                                    &(ifoPtr->timeData->epoch),
                                                                    0.0,
                                                                    ifoPtr->timeData->deltaT,
                                                                    &lalDimensionlessUnit,
                                                                    ifoPtr->timeData->data->length);
                ifoPtr->timeModelhCross = XLALCreateREAL8TimeSeries("timeModelhCross",
                                                                    &(ifoPtr->timeData->epoch),
                                                                    0.0,
                                                                    ifoPtr->timeData->deltaT,
                                                                    &lalDimensionlessUnit,
                                                                    ifoPtr->timeData->data->length);
                ifoPtr->freqModelhPlus = XLALCreateCOMPLEX16FrequencySeries("freqModelhPlus",
                                                                            &(ifoPtr->freqData->epoch),
                                                                            0.0,
                                                                            ifoPtr->freqData->deltaF,
                                                                            &lalDimensionlessUnit,
                                                                            ifoPtr->freqData->data->length);
                ifoPtr->freqModelhCross = XLALCreateCOMPLEX16FrequencySeries("freqModelhCross",
                                                                             &(ifoPtr->freqData->epoch),
                                                                             0.0,
                                                                             ifoPtr->freqData->deltaF,
                                                                             &lalDimensionlessUnit,
                                                                             ifoPtr->freqData->data->length);
                ifoPtr->modelParams = XLALCalloc(1, sizeof(LALInferenceVariables));
            }
            ifoPtr = ifoPtr->next;
        }
        irs->currentLikelihood=LALInferenceNullLogLikelihood(irs->data);
        printf("Injection Null Log Likelihood: %g\n", irs->currentLikelihood);
    } else {
        fprintf(stdout, " initialize(): no data read.\n");
        irs = NULL;
        return(irs);
    }
  
    /* Turn off differential evolution */
    irs->differentialPoints = NULL;
    irs->differentialPointsLength = 0;
    irs->differentialPointsSize = 0;
    return(irs);
}

/********** Initialise MCMC structures *********/

/************************************************/
void initializeMCMC(LALInferenceRunState *runState) {
    char help[]="\
                ---------------------------------------------------------------------------------------------------\n\
                 --- General Algorithm Parameters ------------------------------------------------------------------\n\
                 ---------------------------------------------------------------------------------------------------\n\
                 (--nsteps n)                     Total number of steps for all walkers to make (100000).\n\
                 (--skip n)                       Number of steps between writes to file (100).\n\
                 (--update-interval n)            Number of steps between ensemble updates (1000).\n\
                 (--randomseed seed)              Random seed of sampling distribution (random).\n\
                 \n\
                 ---------------------------------------------------------------------------------------------------\n\
                 --- Likelihood Functions --------------------------------------------------------------------------\n\
                 ---------------------------------------------------------------------------------------------------\n\
                 (--zeroLogLike)                  Use flat, null likelihood.\n\
                 (--studentTLikelihood)           Use the Student-T Likelihood that marginalizes over noise.\n\
                 (--correlatedGaussianLikelihood) Use analytic, correlated Gaussian for Likelihood.\n\
                 (--bimodalGaussianLikelihood)    Use analytic, bimodal correlated Gaussian for Likelihood.\n\
                 (--rosenbrockLikelihood)         Use analytic, Rosenbrock banana for Likelihood.\n\
                 (--analyticnullprior)            Use analytic null prior.\n\
                 (--nullprior)                    Use null prior in the sampled parameters.\n\
                 (--noiseonly)                    Use signal-free log likelihood (noise model only).\n\
                 \n\
                 ---------------------------------------------------------------------------------------------------\n\
                 --- Noise Model -----------------------------------------------------------------------------------\n\
                 ---------------------------------------------------------------------------------------------------\n\
                 (--psdFit)                       Run with PSD fitting\n\
                 (--psdNblock)                    Number of noise parameters per IFO channel (8)\n\
                 (--psdFlatPrior)                 Use flat prior on psd parameters (Gaussian)\n\
                 (--removeLines)                  Do include persistent PSD lines in fourier-domain integration\n\
                 (--KSlines)                      Run with the KS test line removal\n\
                 (--KSlinesWidth)                 Width of the lines removed by the KS test (deltaF)\n\
                 (--chisquaredlines)              Run with the Chi squared test line removal\n\
                 (--chisquaredlinesWidth)         Width of the lines removed by the Chi squared test (deltaF)\n\
                 (--powerlawlines)                Run with the power law line removal\n\
                 (--powerlawlinesWidth)           Width of the lines removed by the power law test (deltaF)\n\
                 (--xcorrbands)                   Run PSD fitting with correlated frequency bands\n\
                 \n\
                 ---------------------------------------------------------------------------------------------------\n\
                 --- Output ----------------------------------------------------------------------------------------\n\
                 ---------------------------------------------------------------------------------------------------\n\
                 (--data-dump)                    Output waveforms to file.\n\
                 (--outfile file)                 Write output files <file>.<chain_number> (ensemble.output.<random_seed>.<walker_number>).\n";

    INT4 walker, i;
    MPI_Comm_rank(MPI_COMM_WORLD, &walker);

    /* Print command line arguments if runState was not allocated */
    if(runState==NULL) {
        fprintf(stdout,"%s",help);
        return;
    }

    INT4 tmpi=0;
    unsigned int randomseed=0;
    ProcessParamsTable *commandLine = runState->commandLine;
    ProcessParamsTable *ppt = NULL;
    FILE *devrandom;
    struct timeval tv;

    /* Print command line arguments if help requested */
    if (LALInferenceGetProcParamVal(runState->commandLine,"--help")) {
        fprintf(stdout,"%s",help);
        return;
    }

    /* Initialise parameters structure */
    runState->algorithmParams = XLALCalloc(1,sizeof(LALInferenceVariables));
    runState->priorArgs = XLALCalloc(1,sizeof(LALInferenceVariables));
    runState->proposalArgs = XLALCalloc(1,sizeof(LALInferenceVariables));

    /* Set up the appropriate functions for the MCMC algorithm */
    runState->algorithm = &ensemble_sampler;
    runState->proposal = &LALInferenceClusteredKDEProposal;

    /* Choose the template generator for inspiral signals */
    LALInferenceInitCBCTemplate(runState);

    UINT4 malmquist = 0;
    if(LALInferenceGetProcParamVal(commandLine,"--skyLocPrior")){
        runState->prior=&LALInferenceInspiralSkyLocPrior;
    } else if (LALInferenceGetProcParamVal(commandLine, "--correlatedGaussianLikelihood") || 
               LALInferenceGetProcParamVal(commandLine, "--bimodalGaussianLikelihood") ||
               LALInferenceGetProcParamVal(commandLine, "--rosenbrockLikelihood") ||
               LALInferenceGetProcParamVal(commandLine, "--analyticnullprior")) {
        runState->prior=&LALInferenceAnalyticNullPrior;
    } else if (LALInferenceGetProcParamVal(commandLine, "--nullprior")) {
        runState->prior=&LALInferenceNullPrior;
    } else if (LALInferenceGetProcParamVal(commandLine, "--malmquistprior")) {
        printf("Using malmquist prior.\n");
        malmquist = 1;
        LALInferenceAddVariable(runState->priorArgs, "malmquist", &malmquist, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
        runState->prior=&LALInferenceInspiralPrior;
    } else {
        runState->prior=&LALInferenceInspiralPriorNormalised;
    }

    if (malmquist) {
        REAL8 malmquist_loudest = 0.0;
        REAL8 malmquist_second_loudest = 5.0;
        REAL8 malmquist_network = 0.0;

        ppt = LALInferenceGetProcParamVal(commandLine,"--malmquist-loudest-snr");
        if (ppt)
            malmquist_loudest = atof(ppt->value);

        ppt = LALInferenceGetProcParamVal(commandLine,"--malmquist-second-loudest-snr");
        if (ppt)
            malmquist_second_loudest = atof(ppt->value);

        ppt = LALInferenceGetProcParamVal(commandLine,"--malmquist-network-snr");
        if (ppt)
            malmquist_network = atof(ppt->value);

        LALInferenceAddVariable(runState->priorArgs, "malmquist_loudest_snr",
                &malmquist_loudest, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        LALInferenceAddVariable(runState->priorArgs, "malmquist_second_loudest_snr",
                &malmquist_second_loudest, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        LALInferenceAddVariable(runState->priorArgs, "malmquist_network_snr",
                &malmquist_network, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
    }

    /* Print more stuff */
    UINT4 verbose=0;
    if (LALInferenceGetProcParamVal(commandLine, "--verbose"))
        verbose=1;

    /* Number of steps between ensemble updates */
    UINT4 nsteps = 100000;
    ppt = LALInferenceGetProcParamVal(commandLine, "--nsteps");
    if(ppt)
        nsteps = atoi(ppt->value);

    /* Print sample every skip iterations */
    UINT4 skip = 100;
    ppt = LALInferenceGetProcParamVal(commandLine, "--skip");
    if(ppt)
        skip = atoi(ppt->value);

    /* Update ensemble every *update_interval* iterations */
    UINT4 update_interval = 1000;
    ppt = LALInferenceGetProcParamVal(commandLine, "--update-interval");
    if(ppt)
        update_interval = atoi(ppt->value);

    /* Keep track of time if benchmarking */
    UINT4 benchmark = 0;
    if (LALInferenceGetProcParamVal(runState->commandLine, "--benchmark"))
        benchmark = 1;

    /* Initialize a random number generator. */
    gsl_rng_env_setup();
    runState->GSLrandom = gsl_rng_alloc(gsl_rng_mt19937);

    /* Use clocktime if seed isn't provided */
    ppt = LALInferenceGetProcParamVal(commandLine, "--randomseed");
    if (ppt)
        randomseed = atoi(ppt->value);
    else {
        if ((devrandom = fopen("/dev/urandom","r")) == NULL) {
            if (walker == 0) {
                gettimeofday(&tv, 0);
                randomseed = tv.tv_sec + tv.tv_usec;
            }
            MPI_Bcast(&randomseed, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        } else {
            if (walker == 0) {
                fread(&randomseed, sizeof(randomseed), 1, devrandom);
                fclose(devrandom);
            }
            MPI_Bcast(&randomseed, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        }
    }

    fprintf(stdout, " initialize(): random seed: %u\n", randomseed);

    gsl_rng_set(runState->GSLrandom, randomseed);

    /* Now make sure that everyone is running with un-correlated
       jumps!  We re-seed rank i process with the ith output of
       the RNG stream from the rank 0 process. Otherwise the
       random stream is the same across all processes. */
    for (i = 0; i < walker; i++)
        randomseed = gsl_rng_get(runState->GSLrandom);

    gsl_rng_set(runState->GSLrandom, randomseed);

    /* Store flags to keep from checking the command line all the time */
    LALInferenceAddVariable(runState->algorithmParams,"verbose", &verbose,
                            LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);

    LALInferenceAddVariable(runState->algorithmParams, "nsteps", &nsteps,
                            LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);

    LALInferenceAddVariable(runState->algorithmParams, "skip", &skip,
                            LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);

    LALInferenceAddVariable(runState->algorithmParams, "update_interval", &update_interval,
                            LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);

    LALInferenceAddVariable(runState->algorithmParams, "benchmark", &benchmark,
                            LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);

    LALInferenceAddVariable(runState->algorithmParams, "random_seed", &randomseed,
                            LALINFERENCE_UINT4_t,LALINFERENCE_PARAM_FIXED);

    return;
}


int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);

    INT4 walker;
    MPI_Comm_rank(MPI_COMM_WORLD, &walker);

    if (walker == 0) fprintf(stdout," ========== lalinference_ensemble ==========\n");

    /* Read command line and parse */
    ProcessParamsTable *procParams=LALInferenceParseCommandLine(argc,argv);

    /* initialise runstate based on command line */
    /* This includes reading in the data */
    /* And performing any injections specified */
    /* And allocating memory */
    LALInferenceRunState *runState = initialize(procParams);

    /* Set up structures for MCMC */
    initializeMCMC(runState);

    /* Set up currentParams with variables to be used */
    LALInferenceInitCBCVariables(runState);

    /* Choose the likelihood */
    LALInferenceInitLikelihood(runState);
 
    /* Call the extra code that was removed from previous function */
    LALInferenceInitEnsemble(runState);
 
    if(runState==NULL) {
        fprintf(stderr, "runState not allocated (%s, line %d).\n",
                __FILE__, __LINE__);
        exit(1);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    /* Call MCMC algorithm */
    runState->algorithm(runState);

    if (walker == 0) printf(" ========== main(): finished. ==========\n");

    MPI_Finalize();
    return 0;
}
