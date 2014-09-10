/*
 *  LALInferenceEnsemble.c:  Bayesian Followup function testing site
 *
 *  Copyright (C) 2014 Ben Farr
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
#include <lal/LALInference.h>
#include "LALInferenceEnsembleSampler.h"
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceProposal.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceReadData.h>
#include <lal/LALInferenceInit.h>

#include <mpi.h>


LALInferenceRunState *initialize(ProcessParamsTable *commandLine);
void initializeMCMC(LALInferenceRunState *runState);
void sample_prior(LALInferenceRunState *runState);


/* Read samples from file as initial state of ensemble */
void LALInferenceInitEnsemble(LALInferenceRunState *state);
void LALInferenceInitEnsemble(LALInferenceRunState *state) {
    LALInferenceVariableItem *item;
    LALInferenceIFOData *headData;
    ProcessParamsTable *ppt;
    REAL8 *sampleArray = NULL;
    UINT4 i=0;

    INT4 MPIrank, MPIsize, walker;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
    MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);
    INT4 nwalkers_per_thread = *(INT4 *) LALInferenceGetVariable(state->algorithmParams, "nwalkers_per_thread");

    INT4 ndim = LALInferenceGetVariableDimensionNonFixed(state->currentParamArray[0]);

    ppt = LALInferenceGetProcParamVal(state->commandLine, "--init-samples");
    if (ppt) {
        char *infile = ppt->value;
        FILE *input = fopen(infile, "r");

        char params[128][128];
        INT4 *col_order = XLALMalloc(ndim *sizeof(INT4));
        UINT4 ncols;

        /* Parse parameter names */
        LALInferenceReadAsciiHeader(input, params, &ncols);

        /* Only cluster parameters that are being sampled */
        UINT4 nvalid_cols=0, j=0;
        UINT4 *valid_cols = XLALMalloc(ncols * sizeof(UINT4));

        for (j = 0; j < ncols; j++) {
            char* internal_param_name = XLALMalloc(512*sizeof(char));
            LALInferenceTranslateExternalToInternalParamName(internal_param_name, params[j]);

            i=0;
            valid_cols[j] = 0;
            for (item = state->currentParamArray[0]->head; item; item = item->next) {
                if (LALInferenceCheckVariableNonFixed(state->currentParamArray[0], item->name)) {
                    if (!strcmp(item->name, internal_param_name)) {
                        col_order[i] = nvalid_cols;
                        nvalid_cols++;
                        valid_cols[j] = 1;
                        break;
                    }
                    i++;
                }
            }

            XLALFree(internal_param_name);
        }

        /* Double check dimensions */
        if (nvalid_cols != ndim) {
            fprintf(stderr, "Inconsistent dimensions for starting state!\n");
            fprintf(stderr, "Sampling in %i dimensions, %i read from file!\n", ndim, nvalid_cols);
            exit(1);
        }

        /* Give a different chunk of samples to each MPI thread */
        INT4 ch, nsamples = 0;
        while ( nsamples < MPIrank*nwalkers_per_thread &&
                (ch = getc(input)) != EOF) {
            if (ch=='\n')
                ++nsamples;
        }

        sampleArray = LALInferenceParseDelimitedAscii(input, ncols, valid_cols, &nwalkers_per_thread);

        REAL8 *parameters = XLALMalloc(ndim * sizeof(REAL8));
        for (walker = 0; walker < nwalkers_per_thread; walker++) {
            for (i = 0; i < ndim; i++)
                parameters[i] = sampleArray[walker*ndim + col_order[i]];
            LALInferenceCopyArrayToVariables(parameters, state->currentParamArray[walker]);
        }

        /* Cleanup */
        XLALFree(col_order);
        XLALFree(valid_cols);
        XLALFree(parameters);
        XLALFree(sampleArray);
    } else {
        #pragma omp parallel for
        for (walker = 0; walker < nwalkers_per_thread; walker++) {
            LALInferenceDrawApproxPrior(state, state->currentParamArray[walker], state->currentParamArray[walker]);
            while (state->prior(state, state->currentParamArray[walker], state->modelArray[walker]) <= -DBL_MAX)
                LALInferenceDrawApproxPrior(state, state->currentParamArray[walker], state->currentParamArray[walker]);
        }
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
        REAL8 d = *(REAL8 *)LALInferenceGetVariable(state->currentParamArray[0], "distance");
        REAL8 bigD = INFINITY;

        /* Don't store to cache, since distance scaling won't work */
        LALSimInspiralWaveformCache *cache = state->modelArray[0]->waveformCache;
        state->modelArray[0]->waveformCache = NULL;

        LALInferenceSetVariable(state->currentParamArray[0], "distance", &bigD);
        null_likelihood = state->likelihood(state->currentParamArray[0], state->data, state->modelArray[0]);

        /* Restore cache to data structure */
        while (headData != NULL) {
            headData->nullloglikelihood = state->modelArray[0]->loglikelihood;
            headData = headData->next;
        }
        state->modelArray[0]->waveformCache = cache;

        /* Replace finite distance */
        LALInferenceSetVariable(state->currentParamArray[0], "distance", &d);
    }
    LALInferenceAddVariable(state->proposalArgs, "nullLikelihood", &null_likelihood,
                            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);

    /* Initialize starting likelihood and prior */
    state->currentPriors = XLALMalloc(nwalkers_per_thread * sizeof(REAL8));
    state->currentLikelihoods = XLALMalloc(nwalkers_per_thread * sizeof(REAL8));
    for (walker = 0; walker < nwalkers_per_thread; walker++) {
        state->currentPriors[walker] = state->prior(state,
                                                    state->currentParamArray[walker],
                                                    state->modelArray[walker]);

        state->currentLikelihoods[walker] = 0.0;
    }

    /* Distribute ensemble according to prior when randomly initializing */
    if (!LALInferenceGetProcParamVal(state->commandLine, "--init-samples"))
        sample_prior(state);

    /* Set starting likelihood values (prior function hasn't changed) */
    #pragma omp parallel for
    for (walker = 0; walker < nwalkers_per_thread; walker++)
        state->currentLikelihoods[walker] = state->likelihood(state->currentParamArray[walker],
                                                                 state->data,
                                                                 state->modelArray[walker]);
}

LALInferenceRunState *initialize(ProcessParamsTable *commandLine)
/* calls the "ReadData()" function to gather data & PSD from files, */
/* and initializes other variables accordingly.                     */
{
    LALInferenceRunState *irs=NULL;
  
    /* read data from files: */
    irs = XLALCalloc(1, sizeof(LALInferenceRunState));
    irs->commandLine=commandLine;
    LALInferenceCheckOptionsConsistency(commandLine);
    irs->data = LALInferenceReadData(commandLine);
  
    if (irs->data != NULL) {
        LALInferenceInjectInspiralSignal(irs->data,commandLine);
        irs->currentLikelihood=LALInferenceNullLogLikelihood(irs->data);
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

    INT4 i, walker;
    INT4 MPIrank, MPIsize;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
    MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);

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
        LALInferenceAddVariable(runState->priorArgs, "malmquist", &malmquist, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);
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
                &malmquist_loudest, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_OUTPUT);
        LALInferenceAddVariable(runState->priorArgs, "malmquist_second_loudest_snr",
                &malmquist_second_loudest, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_OUTPUT);
        LALInferenceAddVariable(runState->priorArgs, "malmquist_network_snr",
                &malmquist_network, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_OUTPUT);
    }

    /* Print more stuff */
    UINT4 verbose=0;
    if (LALInferenceGetProcParamVal(commandLine, "--verbose"))
        verbose=1;

    /* Determine number of walkers */
    INT4 nwalkers = MPIsize;
    INT4 nwalkers_per_thread = 1;
    ppt = LALInferenceGetProcParamVal(commandLine, "--nwalkers");
    if (ppt) {
        INT4 requested_nwalkers = atoi(ppt->value);

        /* Round up to have consistent number of walkers across MPI threads */
        REAL8 requested_nwalkers_per_thread =  floor((REAL8)requested_nwalkers / (REAL8)MPIsize);
        if (requested_nwalkers % MPIsize != 0.0) {
            nwalkers_per_thread = (INT4)requested_nwalkers_per_thread + 1;
            nwalkers = MPIsize * nwalkers_per_thread;
            printf("Rounding up number of walkers to %i to provide \
                    consistent performance across the %i available \
                    MPI threads.\n", nwalkers, MPIsize);
        } else {
            nwalkers = requested_nwalkers;
            nwalkers_per_thread = (INT4) requested_nwalkers_per_thread;
        }
    }

    /* Number of steps between ensemble updates */
    UINT4 step = 0;
    UINT4 nsteps = 100000;
    ppt = LALInferenceGetProcParamVal(commandLine, "--nsteps");
    if (ppt)
        nsteps = atoi(ppt->value);

    /* Print sample every skip iterations */
    UINT4 skip = 100;
    ppt = LALInferenceGetProcParamVal(commandLine, "--skip");
    if (ppt)
        skip = atoi(ppt->value);

    /* Update ensemble every *update_interval* iterations */
    UINT4 update_interval = 1000;
    ppt = LALInferenceGetProcParamVal(commandLine, "--update-interval");
    if (ppt)
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
            if (MPIrank == 0) {
                gettimeofday(&tv, 0);
                randomseed = tv.tv_sec + tv.tv_usec;
            }
            MPI_Bcast(&randomseed, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        } else {
            if (MPIrank == 0) {
                fread(&randomseed, sizeof(randomseed), 1, devrandom);
                fclose(devrandom);
            }
            MPI_Bcast(&randomseed, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        }
    }

    fprintf(stdout, " initialize(): random seed: %u\n", randomseed);

    gsl_rng_set(runState->GSLrandom, randomseed);

     /* Set up CBC model and parameter array */
    runState->modelArray = XLALMalloc(nwalkers_per_thread * sizeof(LALInferenceModel*));
    runState->currentParamArray = XLALMalloc(nwalkers_per_thread * sizeof(LALInferenceVariables*));

    for (walker = 0; walker < nwalkers_per_thread; walker++) {
        runState->modelArray[walker] = LALInferenceInitCBCModel(runState);

        runState->currentParamArray[walker] = XLALMalloc(sizeof(LALInferenceVariables));
        memset(runState->currentParamArray[walker], 0, sizeof(LALInferenceVariables));
        LALInferenceCopyVariables(runState->modelArray[walker]->params, runState->currentParamArray[walker]);
    }

    /* Have currentParams in runState point to the first parameter set, since currentParams
     *   is often used to count dimensions, determine variable names, etc. */
    runState->currentParams = runState->currentParamArray[0];

    /* Store flags to keep from checking the command line all the time */
    LALInferenceAddVariable(runState->algorithmParams,"verbose", &verbose,
                            LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(runState->algorithmParams, "nwalkers_per_thread", &nwalkers_per_thread,
                            LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(runState->algorithmParams, "nwalkers", &nwalkers,
                            LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(runState->algorithmParams, "nsteps", &nsteps,
                            LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(runState->algorithmParams, "skip", &skip,
                            LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(runState->algorithmParams, "update_interval", &update_interval,
                            LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(runState->algorithmParams, "benchmark", &benchmark,
                            LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(runState->algorithmParams, "random_seed", &randomseed,
                            LALINFERENCE_UINT4_t,LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(runState->algorithmParams, "step", &step,
                            LALINFERENCE_UINT4_t,LALINFERENCE_PARAM_OUTPUT);

   /* Now make sure that everyone is running with un-correlated
       jumps!  We re-seed rank i process with the ith output of
       the RNG stream from the rank 0 process. Otherwise the
       random stream is the same across all processes. */
    for (i = 0; i < MPIrank; i++)
        randomseed = gsl_rng_get(runState->GSLrandom);

    gsl_rng_set(runState->GSLrandom, randomseed);

    return;
}


void sample_prior(LALInferenceRunState *runState) {
    INT4 update_interval, nprior_steps, nsteps;
    INT4 walker, nwalkers_per_thread;
    LALInferenceVariables *algorithmParams = runState->algorithmParams;

    /* Ensure ensemble distributed according to prior */
    nsteps = *(INT4 *) LALInferenceGetVariable(algorithmParams, "nsteps");
    update_interval = *(INT4 *) LALInferenceGetVariable(algorithmParams, "update_interval");
    nwalkers_per_thread = *(INT4 *) LALInferenceGetVariable(algorithmParams, "nwalkers_per_thread");

    nprior_steps = 2 * update_interval - 1;
    LALInferenceSetVariable(algorithmParams, "nsteps", &nprior_steps);

    runState->likelihood = &LALInferenceZeroLogLikelihood;

    printf("Distributing ensemble according to prior.\n");
    runState->algorithm(runState);

    /* Call MCMC algorithm */
    LALInferenceSetVariable(algorithmParams, "nsteps", &nsteps);

    /* Reset likelihood function */
    LALInferenceInitLikelihood(runState);

    printf("Beginning posterior sampling.\n");
}

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);

    INT4 MPIrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

    if (MPIrank == 0) fprintf(stdout," ========== lalinference_ensemble ==========\n");

    /* Read command line and parse */
    ProcessParamsTable *procParams=LALInferenceParseCommandLine(argc,argv);

    /* initialise runstate based on command line */
    /* This includes reading in the data */
    /* And performing any injections specified */
    /* And allocating memory */
    LALInferenceRunState *runState = initialize(procParams);

    /* Set up structures for MCMC */
    initializeMCMC(runState);

    /* Set up model struct and set currentVariables to match the initialized model params */
    runState->model = LALInferenceInitCBCModel(runState);
    runState->currentParams = XLALMalloc(sizeof(LALInferenceVariables));
    memset(runState->currentParams, 0, sizeof(LALInferenceVariables));
    LALInferenceCopyVariables(runState->model->params, runState->currentParams);
  
    /* Set template function in runState, since it's sometimes used */
    runState->templt = runState->model->templt;

    /* Choose the likelihood */
    LALInferenceInitLikelihood(runState);
 
    /* Call the extra code that was removed from previous function */
    LALInferenceInitEnsemble(runState);
 
    if(runState==NULL) {
        fprintf(stderr, "runState not allocated (%s, line %d).\n",
                __FILE__, __LINE__);
        exit(1);
    }

    /* Call MCMC algorithm */
    runState->algorithm(runState);

    if (MPIrank == 0) printf(" ========== main(): finished. ==========\n");

    MPI_Finalize();
    return 0;
}
