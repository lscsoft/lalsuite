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


LALInferenceRunState *initialize(ProcessParamsTable *command_line);
void init_sampler(LALInferenceRunState *run_state);
void init_ensemble(LALInferenceRunState *run_state);
void sample_prior(LALInferenceRunState *run_state);

/* Initialize ensemble randomly or from file */
void init_ensemble(LALInferenceRunState *run_state) {
    LALInferenceVariables *current_param;
    LALInferenceVariableItem *item;
    LALInferenceIFOData *headData;
    ProcessParamsTable *ppt;
    REAL8 *sampleArray = NULL;
    UINT4 i=0;

    /* Determine number of MPI threads, and the
     *   number of walkers run by each thread */
    INT4 mpi_rank, mpi_size, walker;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    INT4 nwalkers_per_thread =
        *(INT4 *) LALInferenceGetVariable(run_state->algorithmParams,
                                            "nwalkers_per_thread");

    current_param = run_state->currentParamArray[0];
    UINT4 ndim =
        LALInferenceGetVariableDimensionNonFixed(current_param);

    ProcessParamsTable *command_line = run_state->commandLine;
    ppt = LALInferenceGetProcParamVal(command_line, "--init-samples");
    if (ppt) {
        char *infile = ppt->value;
        FILE *input = fopen(infile, "r");

        char params[128][128];
        UINT4 *col_order = XLALMalloc(ndim *sizeof(INT4));
        UINT4 ncols;

        /* Parse parameter names */
        LALInferenceReadAsciiHeader(input, params, &ncols);

        /* Only cluster parameters that are being sampled */
        UINT4 nvalid_cols=0, j=0;
        UINT4 *valid_cols = XLALMalloc(ncols * sizeof(UINT4));

        for (j = 0; j < ncols; j++) {
            char* internal_pname = XLALMalloc(512*sizeof(char));
            LALInferenceTranslateExternalToInternalParamName(internal_pname,
                                                                params[j]);

            i=0;
            valid_cols[j] = 0;
            for (item = current_param->head; item; item = item->next) {
                if (LALInferenceCheckVariableNonFixed(current_param,
                                                        item->name)) {
                    if (!strcmp(item->name, internal_pname)) {
                        col_order[i] = nvalid_cols;
                        nvalid_cols++;
                        valid_cols[j] = 1;
                        break;
                    }
                    i++;
                }
            }

            XLALFree(internal_pname);
        }

        /* Double check dimensions */
        if (nvalid_cols != ndim) {
            fprintf(stderr, "Inconsistent dimensions for starting state!\n");
            fprintf(stderr, "Sampling in %i dimensions, %i read from file!\n",
                    ndim, nvalid_cols);
            exit(1);
        }

        /* Give a different chunk of samples to each MPI thread */
        INT4 ch, nsamples = 0;
        while ( nsamples < mpi_rank*nwalkers_per_thread &&
                (ch = getc(input)) != EOF) {
            if (ch=='\n')
                ++nsamples;
        }

        UINT4 nlines = (UINT4) nwalkers_per_thread;
        sampleArray = LALInferenceParseDelimitedAscii(input,
                                                        ncols, valid_cols,
                                                        &nlines);

        REAL8 *parameters = XLALMalloc(ndim * sizeof(REAL8));
        for (walker = 0; walker < nwalkers_per_thread; walker++) {
            for (i = 0; i < ndim; i++)
                parameters[i] = sampleArray[walker*ndim + col_order[i]];

            LALInferenceCopyArrayToVariables(parameters,
                                              run_state->currentParamArray[walker]);
        }

        /* Cleanup */
        XLALFree(col_order);
        XLALFree(valid_cols);
        XLALFree(parameters);
        XLALFree(sampleArray);
    } else {
        #pragma omp parallel for
        for (walker = 0; walker < nwalkers_per_thread; walker++) {
            LALInferenceDrawApproxPrior(run_state,
                                        run_state->currentParamArray[walker],
                                        run_state->currentParamArray[walker]);
            while (run_state->prior(run_state,
                                run_state->currentParamArray[walker],
                                run_state->modelArray[walker]) <= -DBL_MAX) {
                LALInferenceDrawApproxPrior(run_state,
                                            run_state->currentParamArray[walker],
                                            run_state->currentParamArray[walker]);
            }
        }
    }

    /* Determine null loglikelihood to be subtracted from printed likelihoods */
    REAL8 null_likelihood = 0.0;
    if (run_state->likelihood==&LALInferenceUndecomposedFreqDomainLogLikelihood ||
        run_state->likelihood==&LALInferenceFreqDomainLogLikelihood){

        null_likelihood = LALInferenceNullLogLikelihood(run_state->data);

    /* If no simple null likelihood method exists, scale signal away */
    } else if (run_state->likelihood == &LALInferenceFreqDomainStudentTLogLikelihood ||
                (run_state->likelihood == &LALInferenceMarginalisedTimeLogLikelihood &&
                !(LALInferenceGetProcParamVal(command_line,"--psdFit")))) {

        headData = run_state->data;
        REAL8 d = *(REAL8 *)LALInferenceGetVariable(run_state->currentParamArray[0],
                                                    "distance");
        REAL8 bigD = INFINITY;

        /* Don't store to cache, since distance scaling won't work */
        LALSimInspiralWaveformCache *cache =
            run_state->modelArray[0]->waveformCache;
        run_state->modelArray[0]->waveformCache = NULL;

        LALInferenceSetVariable(run_state->currentParamArray[0], "distance", &bigD);
        null_likelihood = run_state->likelihood(run_state->currentParamArray[0],
                                            run_state->data,
                                            run_state->modelArray[0]);

        /* Restore cache to data structure */
        while (headData != NULL) {
            headData->nullloglikelihood = run_state->modelArray[0]->loglikelihood;
            headData = headData->next;
        }
        run_state->modelArray[0]->waveformCache = cache;

        /* Replace finite distance */
        LALInferenceSetVariable(run_state->currentParamArray[0], "distance", &d);
    }

    LALInferenceAddVariable(run_state->proposalArgs, "nullLikelihood",
                            &null_likelihood, LALINFERENCE_REAL8_t,
                            LALINFERENCE_PARAM_FIXED);

    /* Initialize starting likelihood and prior */
    run_state->currentPriors = XLALMalloc(nwalkers_per_thread * sizeof(REAL8));
    run_state->currentLikelihoods = XLALMalloc(nwalkers_per_thread * sizeof(REAL8));
    for (walker = 0; walker < nwalkers_per_thread; walker++) {
        run_state->currentPriors[walker] =
            run_state->prior(run_state,
                            run_state->currentParamArray[walker],
                            run_state->modelArray[walker]);

        run_state->currentLikelihoods[walker] = 0.0;
    }

    /* Distribute ensemble according to prior when randomly initializing */
    if (!LALInferenceGetProcParamVal(command_line, "--init-samples") &&
            !LALInferenceGetProcParamVal(command_line, "--skip-prior")) {

        if (mpi_rank == 0)
            printf("Distributing ensemble according to prior.\n");

        sample_prior(run_state);

        if (mpi_rank == 0)
            printf("Completed prior sampling.\n");
    }

    /* Set starting likelihood values (prior function hasn't changed) */
    #pragma omp parallel for
    for (walker = 0; walker < nwalkers_per_thread; walker++)
        run_state->currentLikelihoods[walker] =
            run_state->likelihood(run_state->currentParamArray[walker],
                                run_state->data,
                                run_state->modelArray[walker]);
}

LALInferenceRunState *initialize(ProcessParamsTable *command_line)
/* calls the "ReadData()" function to gather data & PSD from files, */
/* and initializes other variables accordingly.                     */
{
    LALInferenceRunState *run_state = NULL;
    run_state = XLALCalloc(1, sizeof(LALInferenceRunState));

    /* Check that command line is consistent first */
    LALInferenceCheckOptionsConsistency(command_line);
    run_state->commandLine = command_line;

    /* Read data from files or generate fake data */
    run_state->data = LALInferenceReadData(command_line);
  
    /* Perform injections if data successful read or created */
    if (run_state->data != NULL) {
        LALInferenceInjectInspiralSignal(run_state->data,command_line);
        run_state->currentLikelihood=LALInferenceNullLogLikelihood(run_state->data);
    } else {
        run_state = NULL;
        return(run_state);
    }
  
    /* Turn off differential evolution */
    run_state->differentialPoints = NULL;
    run_state->differentialPointsLength = 0;
    run_state->differentialPointsSize = 0;

    return(run_state);
}

/********** Initialise MCMC structures *********/

/************************************************/
void init_sampler(LALInferenceRunState *run_state) {
    char help[]="\
                ---------------------------------------------------------------------------------------------------\n\
                 --- General Algorithm Parameters ------------------------------------------------------------------\n\
                 ---------------------------------------------------------------------------------------------------\n\
                 (--nsteps n)                     Total number of steps for all walkers to make (10000).\n\
                 (--skip n)                       Number of steps between writing samples to file (100).\n\
                 (--update-interval n)            Number of steps between ensemble updates (100).\n\
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
    INT4 mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    /* Print command line arguments if run_state was not allocated */
    if(run_state==NULL) {
        if (mpi_rank == 0)
            fprintf(stdout, "%s", help);
        return;
    }

    unsigned int randomseed=0;
    ProcessParamsTable *command_line = run_state->commandLine;
    ProcessParamsTable *ppt = NULL;
    FILE *devrandom;
    struct timeval tv;

    /* Print command line arguments if help requested */
    if (LALInferenceGetProcParamVal(command_line,"--help")) {
        if (mpi_rank == 0)
            fprintf(stdout, "%s", help);
        return;
    }

    /* Initialise parameters structure */
    run_state->algorithmParams = XLALCalloc(1,sizeof(LALInferenceVariables));
    run_state->priorArgs = XLALCalloc(1,sizeof(LALInferenceVariables));
    run_state->proposalArgs = XLALCalloc(1,sizeof(LALInferenceVariables));

    /* Set up the appropriate functions for the MCMC algorithm */
    run_state->algorithm = &ensemble_sampler;
    run_state->proposal = &LALInferenceClusteredKDEProposal;

    /* Choose the template generator for inspiral signals */
    LALInferenceInitCBCTemplate(run_state);

    UINT4 malmquist = 0;
    if(LALInferenceGetProcParamVal(command_line,"--skyLocPrior")){
        run_state->prior=&LALInferenceInspiralSkyLocPrior;
    } else if (LALInferenceGetProcParamVal(command_line, "--correlatedGaussianLikelihood") || 
               LALInferenceGetProcParamVal(command_line, "--bimodalGaussianLikelihood") ||
               LALInferenceGetProcParamVal(command_line, "--rosenbrockLikelihood") ||
               LALInferenceGetProcParamVal(command_line, "--analyticnullprior")) {
        run_state->prior=&LALInferenceAnalyticNullPrior;
    } else if (LALInferenceGetProcParamVal(command_line, "--nullprior")) {
        run_state->prior=&LALInferenceNullPrior;
    } else if (LALInferenceGetProcParamVal(command_line, "--malmquistprior")) {
        malmquist = 1;
        LALInferenceAddVariable(run_state->priorArgs,
                                "malmquist", &malmquist,
                                LALINFERENCE_UINT4_t,
                                LALINFERENCE_PARAM_OUTPUT);
        run_state->prior=&LALInferenceInspiralPrior;
    } else {
        run_state->prior=&LALInferenceInspiralPriorNormalised;
    }

    if (malmquist) {
        REAL8 malmquist_loudest = 0.0;
        REAL8 malmquist_second_loudest = 5.0;
        REAL8 malmquist_network = 0.0;

        ppt = LALInferenceGetProcParamVal(command_line,
                                            "--malmquist-loudest-snr");
        if (ppt)
            malmquist_loudest = atof(ppt->value);

        ppt = LALInferenceGetProcParamVal(command_line,
                                            "--malmquist-second-loudest-snr");
        if (ppt)
            malmquist_second_loudest = atof(ppt->value);

        ppt = LALInferenceGetProcParamVal(command_line,
                                            "--malmquist-network-snr");
        if (ppt)
            malmquist_network = atof(ppt->value);

        LALInferenceAddVariable(run_state->priorArgs,
                                "malmquist_loudest_snr",
                                &malmquist_loudest,
                                LALINFERENCE_REAL8_t,
                                LALINFERENCE_PARAM_OUTPUT);

        LALInferenceAddVariable(run_state->priorArgs,
                                "malmquist_second_loudest_snr",
                                &malmquist_second_loudest,
                                LALINFERENCE_REAL8_t,
                                LALINFERENCE_PARAM_OUTPUT);

        LALInferenceAddVariable(run_state->priorArgs,
                                "malmquist_network_snr",
                                &malmquist_network,
                                LALINFERENCE_REAL8_t,
                                LALINFERENCE_PARAM_OUTPUT);
    }

    /* Print more stuff */
    UINT4 verbose=0;
    if (LALInferenceGetProcParamVal(command_line, "--verbose"))
        verbose=1;

    /* Determine number of walkers */
    INT4 nwalkers = mpi_size;
    INT4 nwalkers_per_thread = 1;
    ppt = LALInferenceGetProcParamVal(command_line, "--nwalkers");
    if (ppt) {
        INT4 requested_nwalkers = atoi(ppt->value);

        /* Round up to have consistent number of walkers across MPI threads */
        REAL8 requested_nwalkers_per_thread =
            floor((REAL8)requested_nwalkers / (REAL8)mpi_size);

        if (requested_nwalkers % mpi_size != 0.0) {
            nwalkers_per_thread = (INT4)requested_nwalkers_per_thread + 1;
            nwalkers = mpi_size * nwalkers_per_thread;

            if (mpi_rank == 0)
                printf("Rounding up number of walkers to %i to provide \
                        consistent performance across the %i available \
                        MPI threads.\n", nwalkers, mpi_size);
        } else {
            nwalkers = requested_nwalkers;
            nwalkers_per_thread = (INT4) requested_nwalkers_per_thread;
        }
    }

    /* Number of steps between ensemble updates */
    UINT4 step = 0;
    UINT4 nsteps = 10000;
    ppt = LALInferenceGetProcParamVal(command_line, "--nsteps");
    if (ppt)
        nsteps = atoi(ppt->value);

    /* Print sample every skip iterations */
    UINT4 skip = 100;
    ppt = LALInferenceGetProcParamVal(command_line, "--skip");
    if (ppt)
        skip = atoi(ppt->value);

    /* Update ensemble every *update_interval* iterations */
    UINT4 update_interval = 1000;
    ppt = LALInferenceGetProcParamVal(command_line, "--update-interval");
    if (ppt)
        update_interval = atoi(ppt->value);

    /* Keep track of time if benchmarking */
    UINT4 benchmark = 0;
    if (LALInferenceGetProcParamVal(command_line, "--benchmark"))
        benchmark = 1;

    /* Initialize a random number generator. */
    gsl_rng_env_setup();
    run_state->GSLrandom = gsl_rng_alloc(gsl_rng_mt19937);

    /* Use clocktime if seed isn't provided */
    ppt = LALInferenceGetProcParamVal(command_line, "--randomseed");
    if (ppt)
        randomseed = atoi(ppt->value);
    else {
        if ((devrandom = fopen("/dev/urandom","r")) == NULL) {
            if (mpi_rank == 0) {
                gettimeofday(&tv, 0);
                randomseed = tv.tv_sec + tv.tv_usec;
            }
            MPI_Bcast(&randomseed, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        } else {
            if (mpi_rank == 0) {
                fread(&randomseed, sizeof(randomseed), 1, devrandom);
                fclose(devrandom);
            }
            MPI_Bcast(&randomseed, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
        }
    }

    if (mpi_rank == 0)
        fprintf(stdout, " initialize(): random seed: %u\n", randomseed);

    gsl_rng_set(run_state->GSLrandom, randomseed);

    /* Set up CBC model and parameter array */
    run_state->modelArray =
        XLALMalloc(nwalkers_per_thread * sizeof(LALInferenceModel*));

    run_state->currentParamArray =
        XLALMalloc(nwalkers_per_thread * sizeof(LALInferenceVariables*));

    run_state->currentPropDensityArray =
        XLALMalloc(nwalkers_per_thread * sizeof(REAL8));

    for (walker = 0; walker < nwalkers_per_thread; walker++) {
        run_state->currentPropDensityArray[walker] = -DBL_MAX;
        run_state->modelArray[walker] = LALInferenceInitCBCModel(run_state);

        run_state->currentParamArray[walker] =
            XLALMalloc(sizeof(LALInferenceVariables));

        memset(run_state->currentParamArray[walker], 0,
                sizeof(LALInferenceVariables));

        LALInferenceCopyVariables(run_state->modelArray[walker]->params,
                                    run_state->currentParamArray[walker]);
    }

    /* Have currentParams in run_state point to the first parameter set.
     *  This isn't used explicitly by the ensemble sampler, but it is often
     *  used to count dimensions, determine variable names, etc. */
    run_state->currentParams = run_state->currentParamArray[0];

    /* Store flags to keep from checking the command line all the time */
    LALInferenceVariables *algorithm_params = run_state->algorithmParams;

    LALInferenceAddVariable(algorithm_params, "verbose", &verbose,
                            LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(algorithm_params, "nwalkers_per_thread", &nwalkers_per_thread,
                            LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(algorithm_params, "nwalkers", &nwalkers,
                            LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(algorithm_params, "nsteps", &nsteps,
                            LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(algorithm_params, "skip", &skip,
                            LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(algorithm_params, "update_interval", &update_interval,
                            LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(algorithm_params, "benchmark", &benchmark,
                            LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(algorithm_params, "random_seed", &randomseed,
                            LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(algorithm_params, "step", &step,
                            LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_OUTPUT);

   /* Now make sure each MPI-thread is running with un-correlated
       jumps. Re-seed this process with the ith output of
       the RNG stream from the rank 0 thread. Otherwise the
       random stream is the same across all threads. */
    for (i = 0; i < mpi_rank; i++)
        randomseed = gsl_rng_get(run_state->GSLrandom);

    gsl_rng_set(run_state->GSLrandom, randomseed);

    return;
}


/* Sample the prior. Useful for defining an initial state for the ensemble */
void sample_prior(LALInferenceRunState *run_state) {
    INT4 update_interval, nprior_steps, nsteps;

    LALInferenceVariables *algorithm_params = run_state->algorithmParams;

    /* Store old algorithm parameters for later restoration */
    nsteps = *(INT4 *) LALInferenceGetVariable(algorithm_params, "nsteps");
    update_interval =
        *(INT4 *) LALInferenceGetVariable(algorithm_params, "update_interval");

    /* Sample prior for two update intervals */
    nprior_steps = 2 * update_interval - 1;
    LALInferenceSetVariable(algorithm_params, "nsteps", &nprior_steps);

    /* Use the "null" likelihood function in order to sample the prior */
    run_state->likelihood = &LALInferenceZeroLogLikelihood;

    /* Run the sampler to completion */
    run_state->algorithm(run_state);

    /* Restore algorithm parameters and likelihood function */
    LALInferenceSetVariable(algorithm_params, "nsteps", &nsteps);
    LALInferenceInitLikelihood(run_state);
}

int main(int argc, char *argv[]){
    INT4 mpi_rank;

    /* Initialize MPI parallelization */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    if (mpi_rank == 0)
        fprintf(stdout," ========== lalinference_ensemble ==========\n");

    /* Read command line and parse */
    ProcessParamsTable *proc_params = LALInferenceParseCommandLine(argc,argv);

    /* Initialise run_state based on command line. This includes allocating
     *   memory, reading data, and performing any injections specified. */
    LALInferenceRunState *run_state = initialize(proc_params);

    if(run_state==NULL) {
        fprintf(stderr, "run_state not allocated (%s, line %d).\n",
                __FILE__, __LINE__);
        exit(1);
    }

    /* Read arguments and set up sampler */
    init_sampler(run_state);

    /* Setup model struct and set current variables to match the
     *  initialized model params */
    run_state->model = LALInferenceInitCBCModel(run_state);
    run_state->currentParams = XLALMalloc(sizeof(LALInferenceVariables));
    memset(run_state->currentParams, 0, sizeof(LALInferenceVariables));
    LALInferenceCopyVariables(run_state->model->params, run_state->currentParams);
  
    /* Set template function in run_state.  This is already specified in the model,
     *   but it is sometimes looked for in the run_state itself */
    run_state->templt = run_state->model->templt;

    /* Choose the likelihood */
    LALInferenceInitLikelihood(run_state);
 
    /* Setup the initial state of the walkers */
    init_ensemble(run_state);
 
    /* Run the sampler to completion */
    run_state->algorithm(run_state);

    if (mpi_rank == 0)
        fprintf(stdout," ==========  sampling complete ==========\n");

    /* Close down MPI parallelization and return */
    MPI_Finalize();

    return 0;
}
