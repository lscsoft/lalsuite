/*
 *  LALInferenceMCMC.c:  Bayesian Followup function testing site
 *
 *  Copyright (C) 2011 Ilya Mandel, Vivien Raymond, Christian Roever,
 *  Marc van der Sluys, John Veitch, Will M. Farr, and Ben Farr
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
#include "LALInferenceMCMCSampler.h"
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferenceProposal.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceReadData.h>
#include <lal/LALInferenceInit.h>
#include <lalapps.h>
#include <lal/LALInferenceCalibrationErrors.h>

#include <mpi.h>


int MPIrank, MPIsize;

static INT4 readSquareMatrix(gsl_matrix *m, UINT4 N, FILE *inp) {
    UINT4 i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            REAL8 value;
            INT4 nread;

            nread = fscanf(inp, " %lg ", &value);

            if (nread != 1) {
                fprintf(stderr, "Cannot read from matrix file (in %s, line %d)\n",
                        __FILE__, __LINE__);
                exit(1);
            }

            gsl_matrix_set(m, i, j, value);
        }
    }

    return 0;
}


void initializeMCMC(LALInferenceRunState *runState);
void init_mpi_randomstate(LALInferenceRunState *run_state);
REAL8 **parseMCMCoutput(char ***params, UINT4 *nInPar, UINT4 *nInSamps, char *infilename, UINT4 burnin);

/* Set the starting seed of rank 0, and give the rest of the threads
    a seed based on it */
void init_mpi_randomstate(LALInferenceRunState *run_state) {
    ProcessParamsTable *ppt = NULL;
    INT4 i, randomseed;
    INT4 mpi_rank, mpi_size;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    /* Broadcast rank=0's randomseed to everyone */
    randomseed = LALInferenceGetINT4Variable(run_state->algorithmParams, "random_seed");
    MPI_Bcast(&randomseed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    LALInferenceSetVariable(run_state->algorithmParams, "random_seed", &randomseed);

    if (mpi_rank == 0)
        printf(" initialize(): random seed: %u\n", randomseed);

    /* Now make sure each MPI-thread is running with un-correlated
        jumps. Re-seed this process with the ith output of
        the RNG stream from the rank 0 thread. Otherwise the
        random stream is the same across all threads. */
     for (i = 0; i < mpi_rank; i++)
         randomseed = gsl_rng_get(run_state->GSLrandom);

     gsl_rng_set(run_state->GSLrandom, randomseed);

     return;
}

void LALInferenceDrawThreads(LALInferenceRunState *state) {
    INT4 c;

    /* If using a malmquist prior, force a strict prior window on distance for starting point, otherwise
     * the approximate prior draws are very unlikely to be within the malmquist prior */
    REAL8 dist_low, dist_high;
    REAL8 restricted_dist_low = 10.0;
    REAL8 restricted_dist_high = 100.0;
    INT4 changed_dist = 0;
    if (LALInferenceCheckVariable(state->priorArgs, "malmquist") && LALInferenceCheckVariableNonFixed(currentParams, "distance")) {
        changed_dist = 1;
        LALInferenceGetMinMaxPrior(state->priorArgs, "distance", &dist_low, &dist_high);
        LALInferenceRemoveMinMaxPrior(state->priorArgs, "distance");
        LALInferenceAddMinMaxPrior(state->priorArgs, "distance", &restricted_dist_low, &restricted_dist_high, LALINFERENCE_REAL8_t);
    }

    /* If the currentParams are not in the prior, overwrite and pick paramaters from the priors. OVERWRITE EVEN USER CHOICES.
     *     (necessary for complicated prior shapes where LALInferenceCyclicReflectiveBound() is not enough */
    #pragma omp parallel for
    for (c = 0; c < state->nthreads; c++) {
        LALInferenceDrawApproxPrior(state,
                                    state->threads[c]->currentParams,
                                    state->threads[c]->currentParams);
        while (run_state->prior(state,
                                state->threads[c]->currentParams,
                                state->threads[c]->model) <= -DBL_MAX) {
            LALInferenceDrawApproxPrior(state,
                                        state->threads[c]->currentParams,
                                        state->threads[c]->currentParams);
        }

        /* Make sure that our initial value is within the
        *     prior-supported volume. */
        LALInferenceCyclicReflectiveBound(run_state->threads[c]->currentParams, priorArgs);

        /* Initialize starting likelihood and prior */
        run_state->threads[c]->currentPrior =
            run_state->prior(run_state,
                             run_state->threads[c]->currentParams,
                             run_state->threads[c]->model);

        run_state->threads[c]->currentLikelihood =
            run_state->likelihood(run_state,
                                  run_state->threads[c]->currentParams,
                                  run_state->threads[c]->model);
    }

    /* Replace distance prior if changed for initial sample draw */
    if (changed_dist) {
        LALInferenceRemoveMinMaxPrior(state->priorArgs, "distance");
        LALInferenceAddMinMaxPrior(state->priorArgs, "distance", &dist_low, &dist_high, LALINFERENCE_REAL8_t);
    }
}

/********** Initialise MCMC structures *********/

/************************************************/
void init_ptmcmc(LALInferenceRunState *runState)
{
  char help[]="\
               ---------------------------------------------------------------------------------------------------\n\
               --- General Algorithm Parameters ------------------------------------------------------------------\n\
               ---------------------------------------------------------------------------------------------------\n\
               (--nsteps n)                     Maximum number of steps to take (1e7).\n\
               (--neff N)                       Number of independent samples to collect (nsteps).\n\
               (--skip n)                       Number of steps between writing samples to file (100).\n\
               (--adapt-tau)                    Adaptation decay power, results in adapt length of 10^tau (5).\n\
               (--no-adapt)                     Do not adapt run.\n\
               (--randomseed seed)              Random seed of sampling distribution (random).\n\
               \n\
               ---------------------------------------------------------------------------------------------------\n\
               --- Parallel Tempering Algorithm Parameters -------------------------------------------------------\n\
               ---------------------------------------------------------------------------------------------------\n\
               (--temp-skip N)                  Number of steps between temperature swap proposals (100).\n\
               (--tempKill N)                   Iteration number to stop temperature swapping (Niter).\n\
               (--ntemp N)                      Number of temperature chains in ladder (as many as needed).\n
               (--temp-min T)                   Lowest temperature for parallel tempering (1.0).\n\
               (--temp-max T)                   Highest temperature for parallel tempering (50.0).\n\
               (--anneal)                       Anneal hot temperature linearly to T=1.0.\n\
               (--annealStart N)                Iteration number to start annealing (5*10^5).\n\
               (--annealLength N)               Number of iterations to anneal all chains to T=1.0 (1*10^5).\n\
               \n\
               ---------------------------------------------------------------------------------------------------\n\
               --- Noise Model -----------------------------------------------------------------------------------\n\
               ---------------------------------------------------------------------------------------------------\n\
               (--psdFit)                       Run with PSD fitting\n\
               (--psdNblock)                    Number of noise parameters per IFO channel (8)\n\
               (--psdFlatPrior)                 Use flat prior on psd parameters (Gaussian)\n\
               (--glitchFit)                    Run with glitch fitting\n\
               (--glitchNmax)                   Maximum number of glitch basis functions per IFO (20)\n\
               \n\
               ---------------------------------------------------------------------------------------------------\n\
               --- Output ----------------------------------------------------------------------------------------\n\
               ---------------------------------------------------------------------------------------------------\n\
               (--data-dump)                    Output waveforms to file.\n\
               (--adapt-verbose)                Output parameter jump sizes and acceptance rate stats to file.\n\
               (--temp-verbose)                 Output temperature swapping stats to file.\n\
               (--prop-verbose)                 Output proposal stats to file.\n\
               (--outfile file)                 Write output files <file>.<chain_number> (PTMCMC.output.<random_seed>.<mpi_thread>).\n";
    INT4 i, t;
    INT4 ntemp_per_thread;
    ProcessParamsTable *command_line, *ppt = NULL;
    LALInferenceThreadState *thread;
    LALInferenceVariables *propArgs;

    /* Send help if runState was not allocated */
    if(runState == NULL || LALInferenceGetProcParamVal(runState->commandLine, "--help")) {
        fprintf(stdout, "%s", help);
        return;
    }

    /* Set up the appropriate functions for the MCMC algorithm */
    runState->algorithm = &PTMCMCAlgorithm;
    runState->evolve = &mcmc_step;
    runState->proposal = &LALInferenceCyclicProposal;

    /* Choose the appropriate swapping method */
    if (LALInferenceGetProcParamVal(command_line, "--varyFlow")) {
        /* Metropolis-coupled MCMC Swap (assumes likelihood function differs between chains).*/
        runState->parallelSwap = &LALInferenceMCMCMCswap;
    } else {
        /* Standard parallel tempering swap. */
        runState->parallelSwap = &LALInferencePTswap;
    }

    command_line = run_state->commandLine;

    /* Store flags to keep from checking the command line all the time */
    LALInferenceVariables *algorithm_params = run_state->algorithmParams;

    /* Print more stuff */
    INT4 verbose = 0;
    if (LALInferenceGetProcParamVal(command_line, "--verbose"))
        verbose = 1;

    /* Step counter */
    INT4 step = 0;

    /* Keep track of time if benchmarking */
    INT4 benchmark = 0;
    if (LALInferenceGetProcParamVal(command_line, "--benchmark"))
        benchmark = 1;

    /* Number of steps between ensemble updates */
    INT4 nsteps = 100000000;
    ppt = LALInferenceGetProcParamVal(command_line, "--nsteps");
    if (ppt)
        nsteps = atoi(ppt->value);

    /* Number of independent samples to collect.  Default behavior (Neff=0)
        is to not calculate Neff and just got *nsteps* iterations. */
    INT4 neff = nsteps;
    ppt = LALInferenceGetProcParamVal(command_line, "--neff");
    if (ppt)
        neff = atoi(ppt->value);

    /* Print sample every skip iterations */
    INT4 skip = 100;
    ppt = LALInferenceGetProcParamVal(command_line, "--skip");
    if (ppt)
        skip = atoi(ppt->value);

    /* Iterations between proposed temperature swaps */
    INT4 Tskip = 100;
    ppt = LALInferenceGetProcParamVal(command_line, "--temp-skip");
    if (ppt)
        Tskip = atoi(ppt->value);

    /* Allow user to restrict size of temperature ladder */
    INT4 ntemp = 0;
    ppt = LALInferenceGetProcParamVal(command_line, "--ntemp");
    if (ppt)
        ntemp = atoi(ppt->value);

    /* Starting temperature of the ladder */
    REAL8 tempMin = 1.0;
    ppt = LALInferenceGetProcParamVal(command_line, "--temp-min");
    if (ppt)
        tempMin = strtod(ppt->value, (char **)NULL);

    /* Starting temperature of the ladder */
    REAL8 tempMax = 50.0;
    ppt = LALInferenceGetProcParamVal(command_line, "--temp-max");
    if (ppt)
        tempMax = strtod(ppt->value, (char **)NULL);

    /* Limit the size of the differential evolution buffer */
    INT4 de_buffer_limit = 1000000;
    ppt = LALInferenceGetProcParamVal(command_line, "--differential-buffer-limit");
    if (ppt)
        de_buffer_limit = atoi(ppt);

    /* Network SNR of trigger */
    REAL8 trigSNR = 0.0;
    ppt = LALInferenceGetProcParamVal(command_line, "--trigger-snr");
    if (ppt)
        trigSNR = strtod(ppt->value, (char **)NULL);

    /* Variable for tracking autocorrelation length */
    INT4 acl = 0;

    /* Print out temperature swapping info */
    INT4 temp_verbose = 0;
    if (LALInferenceGetProcParamVal(command_line, "--temp-verbose"))
        temp_verbose = 1;

    /* Print out adaptation info */
    INT4 adapt_verbose = 0;
    if (LALInferenceGetProcParamVal(command_line, "--adapt-verbose"))
        adapt_verbose = 1;

    /* Save everything in the run state */
    LALInferenceAddVariable(algorithm_params, "verbose", &verbose,
                            LALINFERENCE_INT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(algorithm_params, "benchmark", &benchmark,
                            LALINFERENCE_INT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(algorithm_params, "step", &step,
                            LALINFERENCE_INT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(algorithm_params, "nsteps", &nsteps,
                            LALINFERENCE_INT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(algorithm_params, "skip", &skip,
                            LALINFERENCE_INT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(algorithm_params, "neff", &neff,
                            LALINFERENCE_INT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(algorithm_params, "tskip", &Tskip,
                            LALINFERENCE_INT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(algorithm_params, "ntemp", &ntemp,
                            LALINFERENCE_INT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(algorithm_params, "temp_min", &tempMin,
                            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(algorithm_params, "temp_max", &tempMax,
                            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(algorithm_params, "de_buffer_limit", &de_buffer_limit,
                            LALINFERENCE_INT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(algorithm_params, "trig_snr", &trigSNR,
                            LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(algorithm_params, "mpirank", &mpi_rank,
                            LALINFERENCE_INT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(algorithm_params, "acl", &acl,
                            LALINFERENCE_INT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(algorithm_params, "temp_verbose", &temp_verbose,
                            LALINFERENCE_INT4_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddVariable(algorithm_params, "adapt_verbose", &adapt_verbose,
                            LALINFERENCE_INT4_t, LALINFERENCE_PARAM_OUTPUT);

    /* Build the temperature ladder */
    INT4 ntemp_per_thread = LALInferenceBuildHybridTempLadder(runState);

    /* Initialize the walkers on this MPI thread */
    LALInferenceInitCBCThreads(runState, ntemp_per_thread);

    /* Establish the random state across MPI threads */
    init_mpi_randomstate(runState);

    /* Build the proposals and randomize */
    propArgs = LALInferenceParseProposalArgs(runState);
    for (i=0; i<runState->nthreads; i++) {
        thread = runState->threads[i];
        thread->cycle = LALInferenceSetupDefaultInspiralProposalCycle(propArgs);
        LALInferenceRandomizeProposalCycle(thread->cycle, thread->GSLrandom);
    }
}


INT4 LALInferenceBuildHybridTempLadder(LALInferenceRunState *runState) {
    REAL8 temp, tempMin, tempMax, tempDelta;
    REAL8 targetHotLike;
    REAL8 trigSNR, networkSNRsqrd=0.0;
    REAL8 *ladder=NULL;
    INT4 flexible_tempmax;
    INT4 ndim, ntemp, ntemp_per_thread;
    INT4 mpi_rank, mpi_size;
    INT4 t;
    LALInferenceIFOData *ifo;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    ndim = LALInferenceGetVariableDimensionNonFixed(runState->currentParams);
    ntemp = LALInferenceGetINT4Variable(runState->algorithmParams, "ntemp");
    tempMin = LALInferenceGetREAL8Variable(runState->algorithmParams, "temp_min");
    tempMax = LALInferenceGetREAL8Variable(runState->algorithmParams, "temp_max");

    /* Targeted max 'experienced' log(likelihood) of hottest chain */
    targetHotLike = ndim/2.;

    /* Set maximum temperature (command line value take precidence) */
    if (LALInferenceGetProcParamVal(runState->commandLine,"--temp-max")) {
        if (MPIrank==0)
            fprintf(stdout,"Using tempMax specified by commandline: %f.\n", ladderMax);

    } else if (LALInferenceGetProcParamVal(runState->commandLine,"--trigger-snr")) {        //--trigSNR given, choose ladderMax to get targetHotLike
        trigSNR = LALInferenceGetREAL8Variable(runState->algorithmParams, "trigger-snr");
        networkSNRsqrd = trigSNR * trigSNR;
        tempMax = networkSNRsqrd/(2*targetHotLike);

        if (MPIrank==0)
            fprintf(stdout,"Trigger SNR of %f specified, setting max temperature to %f.\n", trigSNR, tempMax);

    } else {
        /* Determine network SNR if injection was done */
        ifo = runState->data;
        while (ifo != NULL) {
            networkSNRsqrd += ifo->SNR * ifo->SNR;
            ifo = ifo->next;
        }

        if (networkSNRsqrd > 0.0) {
            tempMax = networkSNRsqrd/(2*targetHotLike);
            if (MPIrank==0)
                fprintf(stdout,"Injecting SNR of %f, setting tempMax to %f.\n", sqrt(networkSNRsqrd), tempMax);

        /* If all else fails, use the default */
        } else {
            if (MPIrank==0)
                fprintf(stdout, "No --trigger-snr or --temp-max specified, and \
                        not injecting a signal. Setting max temperature to default of %f.\n", tempMax);
        }
    }

    LALInferenceSetVariable(runState->algorithmParams, "temp_max", &tempMax);

    if (tempMin > tempMax) {
        fprintf(stdout,"WARNING: temp_min > temp_max.  Forcing temp_min=1.0.\n");
        tempMin = 1.0;

        LALInferenceSetVariable(runState->algorithmParams, "temp_min", &tempMin);
    }

    /* Construct temperature ladder */

    /* If ntemp wasn't specified, figure out how big the ladder should be
        to fill between the min and max with the calculated tempDelta */
    if (ntemp == 0) {
        flexible_tempmax = 1;
        tempDelta = 1. + sqrt(2./(REAL8)ndim);

        t = 1;
        temp = tempMin;
        while (temp < tempMax) {
            t++;
            temp = tempMin * pow(tempDelta, t)
        }

        ntemp = t;
    }

    ntemp_per_thread = ntemp / mpi_size;

    /* Increase ntemp to distribute evenly across MPI threads */
    if (ntemp % mpi_size != 0.0) {
        /* Round up to have consistent number of walkers across MPI threads */
        ntemp_per_thread = (INT4)ceil((REAL8)ntemp / (REAL8)mpi_size);
        ntemp = mpi_size * ntemp_per_thread;
    }

    LALInferenceSetVariable(runState->algorithmParams, "ntemp", &ntemp);

    if (flexible_tempmax)
        tempDelta = 1. + sqrt(2./(REAL8)ndim);
    else
        tempDelta = pow(tempMax/tempMin, 1.0/(REAL8)(ntemp-1));

    ladder = XLALCalloc(ntemp, sizeof(REAL8));
    for (t=0; t<ntemp; ++t)
        ladder[t] = tempMin * pow(tempDelta, t);

    for (c = 0; c < ntemp_per_thread; c++)
        runState->threads[c]->temperature = ladder[MPIrank*ntemp_per_thread + c];

    if (MPIrank == 0 && LALInferenceGetVariable(runState->algorithmParams, "verbose")) {
        printf("\nTemperature ladder:\n");
        for (t=0; t<ntemp; ++t)
            printf(" ladder[%d]=%f\n", t, ladder[t]);
    }

    return ntemp_per_thread;
}


REAL8 **parseMCMCoutput(char ***params, UINT4 *nInPar, UINT4 *nInSamps, char *infileName, UINT4 burnin) {
    char str[999];
    char header[999];
    char *word;
    UINT4 nread;
    UINT4 i=0, j=0, nCols=0, nPar=0, par=0, col=0;
    UINT4 cycle=0;
    REAL8 val=0;

    const char *non_params[] = {"cycle","logpost","logprior","logl","loglH1","loglL1","loglV1","",NULL};

    FILE *infile = fopen(infileName,"r");

    fgets(str, 999, infile);
    strcpy(header, str);
    word = strtok(header, " \t");
    // Find column headers
    while (strcmp(word,"cycle") && str != NULL) {
        fgets(str, 999, infile);
        strcpy(header, str);
        word = strtok(header, " \t");
    }

    if (str == NULL) {
        fprintf(stderr, "Couldn't find column headers in file %s\n",infileName);
        exit(1);
    }

    // Read in column names and check if they are parameters
    strcpy(header, str);
    word = strtok(header, " \t");
    while (word != NULL) {
        nCols++;
        word = strtok(NULL, " \t");
    }
    // FIXME Remove a false column due to trailing whitespace
    nCols--;

    UINT4 is_param[nCols];

    strcpy(header, str);
    word = strtok(header, " \t");
    for (i=0; i<nCols; i++) {
        j=0;
        is_param[i] = 1;
        nPar++;
        while (non_params[j] != NULL) {
            if (!strcmp(non_params[j],word)) {
                is_param[i] = 0;
                nPar--;
                break;
            }
            j++;
        }
        word = strtok(NULL, " \t");
    }

    char** in_params = XLALMalloc((nPar)*sizeof(char *));

    word = strtok(str, " \t");
    // Already assumed cycle is the first column, so skip it
    par=0;
    for (i=1; i<nCols; i++) {
        char *param_name = strtok(NULL, " \t");
        if (is_param[i]) {
            in_params[par] = param_name;
            par++;
        }
    }

    printf("Reading the following params from %s:\n", infileName);
    for (par=0; par<nPar; par++)
        printf("\t%s\n",in_params[par]);

    // Move past burnin
    INT4 ch;
    if (burnin > 0) {
        while (cycle <= burnin) {
            fscanf(infile, "%i", &cycle);
            for (j=1;j<nCols;j++)
                fscanf(infile, "%lg", &val);
        }

        // Make sure at end of line
        ch = getc(infile);
        while (ch != '\n') ch = getc(infile);
    }

    // Determine number of samples after burnin
    unsigned long startPostBurnin = ftell(infile);
    UINT4 nSamples=0;

    while ( (ch = getc(infile)) != EOF) {
        if (ch=='\n')
            ++nSamples;
    }
    fseek(infile,startPostBurnin,SEEK_SET);
    printf("%i samples read from %s.\n", nSamples, infileName);

    // Read in samples
    REAL8 **sampleArray;
    sampleArray = (REAL8**) XLALMalloc(nSamples * sizeof(REAL8*));

    for (i = 0; i < nSamples; i++) {
        sampleArray[i] = XLALMalloc(nPar * sizeof(REAL8));

        nread = fscanf(infile, "%i", &cycle);
        if (nread != 1) {
            fprintf(stderr, "Cannot read sample from file (in %s, line %d)\n",
            __FILE__, __LINE__);
            exit(1);
        }

        par=0;
        for (col = 1; col < nCols; col++) {
            nread = fscanf(infile, "%lg", &val);
            if (nread != 1) {
                fprintf(stderr, "Cannot read sample from file (in %s, line %d)\n",
                __FILE__, __LINE__);
                exit(1);
            }

            if (is_param[col]) {
                sampleArray[i][par] = val;
                par++;
            }
        }
    }

    *params = in_params;
    *nInPar = nPar;
    *nInSamps = nSamples;
    return sampleArray;
}

ProcessParamsTable *LALInferenceContinueMCMC(char *infileName) {
    INT4 n;
    INT4 fileargc=1;
    char *infileName, *pch;
    char *fileargv[99], str[999], buffer[99];
    FILE *infile;

    infile = fopen(infileName, "r");
    if (infile==NULL) {
        fprintf(stderr,"Cannot read %s/n", infileName);
        exit (1);
    }

    n = sprintf(buffer, "lalinference_mcmcmpi_from_file_%s", infileName);
    fileargv[0] = (char*)XLALCalloc(n+1, sizeof(char*));
    fileargv[0] = buffer;

    fgets(str, 999, infile);
    fgets(str, 999, infile);
    fclose(infile);

    pch = strtok (str," ");
    while (pch != NULL) {
        if (strcmp(pch, "Command") != 0 && strcmp(pch, "line:") != 0) {
            n = strlen(pch);
            fileargv[fileargc] = (char*)XLALCalloc(n+1, sizeof(char*));
            fileargv[fileargc] = pch;

            fileargc++;
            if(fileargc>=99) {
                fprintf(stderr,"Too many arguments in file %s\n",infileName);
                exit (1);
            }
        }
        pch = strtok (NULL, " ");
    }

    /* In order to get rid of the '\n' than fgets returns when reading the command line. */
    fileargv[fileargc-1][strlen(fileargv[fileargc-1])-1] = '\0';
    procParams = LALInferenceParseCommandLine(fileargc, fileargv);
}


int main(int argc, char *argv[]){
    INT4 mpiRank;
    ProcessParamsTable *procParams = NULL;
    LALInferenceRunState *runState = NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

    if (MPIrank == 0) fprintf(stdout," ========== LALInference_MCMC ==========\n");

    /* Read command line and parse */
    procParams = LALInferenceParseCommandLine(argc, argv);

    ppt = LALInferenceGetProcParamVal(procParams, "--continue-run");
    if (ppt)
        procParams = LALInferenceContinueMCMC(ppt->value);

    /* initialise runstate based on command line */
    /* This includes reading in the data */
    /* And performing any injections specified */
    /* And allocating memory */
    runState = LALInferenceInitCBCRunState(procParams);

    if (runState == NULL) {
        if (LALInferenceGetProcParamVal(proc_params, "--help")) {
            exit(0);
        else {
            fprintf(stderr, "run_state not allocated (%s, line %d).\n",
                    __FILE__, __LINE__);
            exit(1);
        }
    }

    /* Handle PTMCMC setup */
    init_ptmcmc(runState);

    /* Choose the prior */
    LALInferenceInitPrior(runState);

    /* Choose the likelihood */
    LALInferenceInitLikelihood(runState);

    /* Draw starting positions */
    LALInferenceDrawThreads(runState);

    /* Call MCMC algorithm */
    runState->algorithm(runState);

    if (MPIrank == 0) printf(" ========== main(): finished. ==========\n");
    MPI_Finalize();

    return 0;
}
