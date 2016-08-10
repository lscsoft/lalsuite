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


void init_mpi_randomstate(LALInferenceRunState *run_state);
void initializeMCMC(LALInferenceRunState *runState);
INT4 init_ptmcmc(LALInferenceRunState *runState);
REAL8 *LALInferenceBuildHybridTempLadder(LALInferenceRunState *runState, INT4 ndim);
ProcessParamsTable *LALInferenceContinueMCMC(char *infileName);

REAL8 **parseMCMCoutput(char ***params, UINT4 *nInPar, UINT4 *nInSamps, char *infilename, UINT4 burnin);

/********** Initialise MCMC structures *********/

/************************************************/
/* Set the starting seed of rank 0, and give the rest of the threads
    a seed based on it.  This is a cut-and-paste from kombine, which needs
    to be unified in LALInference when MPI-enabled.*/
void init_mpi_randomstate(LALInferenceRunState *run_state) {
    INT4 i, randomseed;
    INT4 mpi_rank;

    mpi_rank = LALInferenceGetINT4Variable(run_state->algorithmParams, "mpirank");

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

INT4 init_ptmcmc(LALInferenceRunState *runState) {
  char help[]="\
    ----------------------------------------------\n\
    --- MCMC Algorithm Parameters ----------------\n\
    ----------------------------------------------\n\
    (--nsteps n)        Maximum number of steps to take (1e7)\n\
    (--neff N)          Number of independent samples to collect (nsteps)\n\
    (--skip n)          Number of steps between writing samples to file (100)\n\
    (--adapt-tau)       Adaptation decay power, results in adapt length of 10^tau (5)\n\
    (--no-adapt)        Do not adapt run\n\
    (--randomseed seed) Random seed of sampling distribution (random)\n\
    \n\
    ----------------------------------------------\n\
    --- Parallel Tempering Algorithm Parameters --\n\
    ----------------------------------------------\n\
    (--temp-skip N)     Number of steps between temperature swap proposals (100)\n\
    (--tempKill N)      Iteration number to stop temperature swapping (Niter)\n\
    (--ntemps N)         Number of temperature chains in ladder (as many as needed)\n\
    (--temp-min T)      Lowest temperature for parallel tempering (1.0)\n\
    (--temp-max T)      Highest temperature for parallel tempering (50.0)\n\
    (--anneal)          Anneal hot temperature linearly to T=1.0\n\
    (--annealStart N)   Iteration number to start annealing (5*10^5)\n\
    (--annealLength N)  Number of iterations to anneal all chains to T=1.0 (1*10^5)\n\
    \n\
    ----------------------------------------------\n\
    --- Noise Model ------------------------------\n\
    ----------------------------------------------\n\
    (--psd-fit)         Run with PSD fitting\n\
    (--psdNblock)       Number of noise parameters per IFO channel (8)\n\
    (--psdFlatPrior)    Use flat prior on psd parameters (Gaussian)\n\
    (--glitch-fit)       Run with glitch fitting\n\
    (--glitchNmax)      Maximum number of glitch basis functions per IFO (20)\n\
    \n\
    ----------------------------------------------\n\
    --- Output -----------------------------------\n\
    ----------------------------------------------\n\
    (--data-dump)       Output waveforms to file\n\
    (--adapt-verbose)   Output parameter jump sizes and acceptance rate stats to file\n\
    (--temp-verbose)    Output temperature swapping stats to file\n\
    (--prop-verbose)    Output proposal stats to file\n\
    (--prop-track)      Output proposal parameters\n\
    (--outfile file)    Write output files <file>.<chain_number> \n\
                            (PTMCMC.output.<random_seed>.<mpi_thread>)\n\
    \n";
    INT4 mpi_rank, mpi_size;
    INT4 ntemps_per_thread;
    INT4 noAdapt;
    INT4 i, ndim;
    INT4 count_vectors = 0;
    ProcessParamsTable *command_line = NULL, *ppt = NULL;
    LALInferenceThreadState *thread = NULL;
    LALInferenceModel *model;
    REAL8 *ladder;

    /* Send help if runState was not allocated */
    if(runState == NULL || LALInferenceGetProcParamVal(runState->commandLine, "--help")) {
        fprintf(stdout, "%s", help);
        LALInferenceInitCBCThreads(runState,0);
        return XLAL_FAILURE;
    }

    command_line = runState->commandLine;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    /* Set up the appropriate functions for the MCMC algorithm */
    runState->algorithm = &PTMCMCAlgorithm;

    /* Choose the appropriate swapping method */
    if (LALInferenceGetProcParamVal(command_line, "--varyFlow")) {
        /* Metropolis-coupled MCMC Swap (assumes likelihood function differs between chains).*/
        //runState->parallelSwap = &LALInferenceMCMCMCswap;
        fprintf(stderr, "ERROR: MCMCMC sampling hasn't been brought up-to-date since restructuring.\n");
        return XLAL_FAILURE;
    } else {
        /* Standard parallel tempering swap. */
        runState->parallelSwap = &LALInferencePTswap;
    }

    /* Store flags to keep from checking the command line all the time */
    LALInferenceVariables *algorithm_params = runState->algorithmParams;

    /* Print more stuff */
    INT4 verbose = 0;
    if (LALInferenceGetProcParamVal(command_line, "--verbose"))
        verbose = 1;

    INT4 propVerbose = 0;
    if (LALInferenceGetProcParamVal(command_line, "--prop-verbose"))
        propVerbose = 1;

    INT4 propTrack = 0;
    if (LALInferenceGetProcParamVal(command_line, "--prop-track"))
        propTrack = 1;

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
    ppt = LALInferenceGetProcParamVal(command_line, "--Neff");
    if (ppt) {
        printf("WARNING: --Neff has been deprecated in favor of --neff\n");
    } else {
        ppt = LALInferenceGetProcParamVal(command_line, "--neff");
    }
    if (ppt)
        neff = atoi(ppt->value);

    /* Print sample every skip iterations */
    INT4 skip = 500;
    ppt = LALInferenceGetProcParamVal(command_line, "--skip");
    if (ppt)
        skip = atoi(ppt->value);

    /* Iterations between proposed temperature swaps */
    INT4 Tskip = 500;
    ppt = LALInferenceGetProcParamVal(command_line, "--temp-skip");
    if (ppt)
        Tskip = atoi(ppt->value);

    /* Counter for triggering PT swaps */
    INT4 nsteps_until_swap = Tskip;

    /* Allow user to restrict size of temperature ladder */
    INT4 ntemps = 0;
    ppt = LALInferenceGetProcParamVal(command_line, "--ntemp");
    if (ppt)
        printf("WARNING: --ntemp has been deprecated in favor of --ntemps\n");
    else
        ppt = LALInferenceGetProcParamVal(command_line, "--ntemps");
    if (ppt)
        ntemps = atoi(ppt->value);

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
        de_buffer_limit = atoi(ppt->value);

    /* Network SNR of trigger */
    REAL8 trigSNR = 0.0;
    ppt = LALInferenceGetProcParamVal(command_line, "--trigger-snr");
    if (ppt)
        trigSNR = strtod(ppt->value, (char **)NULL);

    /* Variable for tracking autocorrelation length */
    INT4 acl = INT_MAX;

    /* Print out temperature swapping info */
    INT4 temp_verbose = 0;
    if (LALInferenceGetProcParamVal(command_line, "--temp-verbose"))
        temp_verbose = 1;

    /* Print out adaptation info */
    INT4 adapt_verbose = 0;
    if (LALInferenceGetProcParamVal(command_line, "--adapt-verbose"))
        adapt_verbose = 1;

    /* Output SNRs */
    INT4 outputSNRs = 0;
    if (LALInferenceGetProcParamVal(command_line, "--output-snrs"))
        outputSNRs = 1;

    /* Save everything in the run state */
    LALInferenceAddINT4Variable(algorithm_params, "verbose", verbose, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddINT4Variable(algorithm_params, "prop_verbose", propVerbose, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddINT4Variable(algorithm_params, "prop_track", propTrack, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddINT4Variable(algorithm_params, "benchmark", benchmark, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddINT4Variable(algorithm_params, "nsteps", nsteps, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddINT4Variable(algorithm_params, "skip", skip, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddINT4Variable(algorithm_params, "neff", neff, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddINT4Variable(algorithm_params, "tskip", Tskip, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddINT4Variable(algorithm_params, "nsteps_until_swap", nsteps_until_swap, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddINT4Variable(algorithm_params, "mpirank", mpi_rank, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddINT4Variable(algorithm_params, "mpisize", mpi_size, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddINT4Variable(algorithm_params, "ntemps", ntemps, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddREAL8Variable(algorithm_params, "temp_min", tempMin, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddREAL8Variable(algorithm_params, "temp_max", tempMax, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddINT4Variable(algorithm_params, "de_buffer_limit", de_buffer_limit, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddREAL8Variable(algorithm_params, "trigger_snr", trigSNR, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddINT4Variable(algorithm_params, "temp_verbose", temp_verbose, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddINT4Variable(algorithm_params, "adapt_verbose", adapt_verbose, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddINT4Variable(algorithm_params, "output_snrs", outputSNRs, LALINFERENCE_PARAM_OUTPUT);

    /* Make a single model just to count dimensions */
    model = LALInferenceInitCBCModel(runState);
    ndim = LALInferenceGetVariableDimensionNonFixedChooseVectors(model->params, count_vectors);

    /* Build the temperature ladder */
    ladder = LALInferenceBuildHybridTempLadder(runState, ndim);
    ntemps_per_thread = LALInferenceGetINT4Variable(runState->algorithmParams, "ntemps_per_thread");

    /* Add some settings settings to runstate proposal args so their copied to threads */
    LALInferenceAddINT4Variable(runState->proposalArgs, "de_skip", skip, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddINT4Variable(runState->proposalArgs, "output_snrs", outputSNRs, LALINFERENCE_PARAM_OUTPUT);

    /* Parse proposal args for runSTate and initialize the walkers on this MPI thread */
    LALInferenceInitCBCThreads(runState, ntemps_per_thread);

    /* Establish the random state across MPI threads */
    init_mpi_randomstate(runState);

    /* Give a new set of proposal args to each thread */
    for (i=0; i<runState->nthreads; i++) {
        thread = runState->threads[i];

        thread->id = mpi_rank*ntemps_per_thread + i;
        thread->temperature = ladder[thread->id];

        thread->proposal = &LALInferenceCyclicProposal;
        thread->proposalArgs = LALInferenceParseProposalArgs(runState);
        thread->cycle = LALInferenceSetupDefaultInspiralProposalCycle(thread->proposalArgs);

        LALInferenceRandomizeProposalCycle(thread->cycle, thread->GSLrandom);

        LALInferenceAddINT4Variable(thread->proposalArgs, "acl",
                                    acl, LALINFERENCE_PARAM_OUTPUT);
    }

    /* Grab adaptation settings from the last thread and add to algorithm params */
    noAdapt = LALInferenceGetINT4Variable(thread->proposalArgs, "no_adapt");

    INT4 adaptTau = 0;
    if (LALInferenceCheckVariable(thread->proposalArgs, "adaptTau"))
        adaptTau = LALInferenceGetINT4Variable(thread->proposalArgs, "adaptTau");

    LALInferenceAddINT4Variable(algorithm_params, "adaptTau",
                                adaptTau, LALINFERENCE_PARAM_OUTPUT);

    INT4 adaptLength = 0;
    if (LALInferenceCheckVariable(thread->proposalArgs, "adaptLength"))
        adaptLength = LALInferenceGetINT4Variable(thread->proposalArgs, "adaptLength");

    LALInferenceAddINT4Variable(algorithm_params, "adaptLength",
                                adaptLength, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddINT4Variable(algorithm_params, "no_adapt",
                                noAdapt, LALINFERENCE_PARAM_OUTPUT);

    /* Set counters to negative numbers, indicating adaptation is happening */
    for (i=0; i<runState->nthreads; i++)
        runState->threads[i]->step = -adaptLength;

    XLALFree(ladder);
    return XLAL_SUCCESS;
}


REAL8 *LALInferenceBuildHybridTempLadder(LALInferenceRunState *runState, INT4 ndim) {
    REAL8 temp, tempMin, tempMax, tempDelta;
    REAL8 targetHotLike;
    REAL8 trigSNR, networkSNRsqrd=0.0;
    REAL8 *ladder=NULL;
    INT4 flexible_tempmax=0;
    INT4 ntemps, ntemps_per_thread;
    INT4 mpi_rank, mpi_size;
    INT4 t;
    LALInferenceIFOData *ifo;

    mpi_rank = LALInferenceGetINT4Variable(runState->algorithmParams, "mpirank");
    mpi_size = LALInferenceGetINT4Variable(runState->algorithmParams, "mpisize");

    ntemps = LALInferenceGetINT4Variable(runState->algorithmParams, "ntemps");
    tempMin = LALInferenceGetREAL8Variable(runState->algorithmParams, "temp_min");
    tempMax = LALInferenceGetREAL8Variable(runState->algorithmParams, "temp_max");
    trigSNR = LALInferenceGetREAL8Variable(runState->algorithmParams, "trigger_snr");

    /* Targeted max 'experienced' log(likelihood) of hottest chain */
    targetHotLike = ndim/2.;

    /* Set maximum temperature (command line value take precidence) */
    if (LALInferenceGetProcParamVal(runState->commandLine,"--temp-max")) {
        if (mpi_rank==0)
            fprintf(stdout,"Using tempMax specified by commandline: %f.\n", tempMax);
    } else if (trigSNR > 0) {
        networkSNRsqrd = trigSNR * trigSNR;
        tempMax = networkSNRsqrd/(2*targetHotLike);

        if (mpi_rank==0)
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
            if (mpi_rank==0)
                fprintf(stdout,"Injecting SNR of %f, setting tempMax to %f.\n", sqrt(networkSNRsqrd), tempMax);

        /* If all else fails, use the default */
        } else {
            if (mpi_rank==0) {
                fprintf(stdout, "No --trigger-snr or --temp-max specified, and ");
                fprintf(stdout, "not injecting a signal. Setting max temperature");
                fprintf(stdout, "to default of %f.\n", tempMax);
            }
        }
    }

    LALInferenceSetVariable(runState->algorithmParams, "temp_max", &tempMax);

    if (tempMin > tempMax) {
        fprintf(stdout,"WARNING: temp_min > temp_max.  Forcing temp_min=1.0.\n");
        tempMin = 1.0;

        LALInferenceSetVariable(runState->algorithmParams, "temp_min", &tempMin);
    }

    /* Construct temperature ladder */

    /* If ntemps wasn't specified, figure out how big the ladder should be
        to fill between the min and max with the calculated tempDelta */
    if (ntemps == 0) {
        flexible_tempmax = 1;
        tempDelta = 1. + sqrt(2./(REAL8)ndim);

        t = 1;
        temp = tempMin;
        while (temp < tempMax) {
            t++;
            temp = tempMin * pow(tempDelta, t);
        }

        ntemps = t;
    }

    ntemps_per_thread = ntemps / mpi_size;

    /* Increase ntemps to distribute evenly across MPI threads */
    if (ntemps % mpi_size != 0.0) {
        /* Round up to have consistent number of walkers across MPI threads */
        ntemps_per_thread = (INT4)ceil((REAL8)ntemps / (REAL8)mpi_size);
        ntemps = mpi_size * ntemps_per_thread;
    }

    LALInferenceAddINT4Variable(runState->algorithmParams, "ntemps",
                                ntemps, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddINT4Variable(runState->algorithmParams, "ntemps_per_thread",
                                ntemps_per_thread, LALINFERENCE_PARAM_OUTPUT);

    if (flexible_tempmax)
        tempDelta = 1. + sqrt(2./(REAL8)ndim);
    else
        tempDelta = pow(tempMax/tempMin, 1.0/(REAL8)(ntemps-1));

    ladder = XLALCalloc(ntemps, sizeof(REAL8));
    for (t=0; t<ntemps; ++t)
        ladder[t] = tempMin * pow(tempDelta, t);

    return ladder;
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
    while (strcmp(word,"cycle")) {
        fgets(str, 999, infile);
        strcpy(header, str);
        word = strtok(header, " \t");
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
    char *pch;
    ProcessParamsTable *procParams = NULL;
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

    return procParams;
}


int main(int argc, char *argv[]){
    INT4 mpirank;
    ProcessParamsTable *procParams = NULL, *ppt = NULL;
    LALInferenceRunState *runState = NULL;
    LALInferenceIFOData *data = NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);

    if (mpirank == 0) fprintf(stdout," ========== LALInference_MCMC ==========\n");

    /* Read command line and parse */
    procParams = LALInferenceParseCommandLine(argc, argv);

    ppt = LALInferenceGetProcParamVal(procParams, "--continue-run");
    if (ppt)
        procParams = LALInferenceContinueMCMC(ppt->value);

    /* initialise runstate based on command line */
    /* This includes reading in the data */
    /* And performing any injections specified */
    /* And allocating memory */
    runState = LALInferenceInitRunState(procParams);

    if (runState == NULL) {
        if (!LALInferenceGetProcParamVal(procParams, "--help")) {
            fprintf(stderr, "run_state not allocated (%s, line %d).\n",
                    __FILE__, __LINE__);
        }
    } else {
        data = runState->data;
    }

    /* Perform injections if data successful read or created */
    if (runState){
      LALInferenceInjectInspiralSignal(data, runState->commandLine);
    }

    /* Simulate calibration errors. 
     * NOTE: this must be called after both ReadData and (if relevant) 
     * injectInspiralTD/FD are called! */
    LALInferenceApplyCalibrationErrors(data, procParams);

    /* Handle PTMCMC setup */
    init_ptmcmc(runState);

    /* Choose the prior */
    LALInferenceInitCBCPrior(runState);

    /* Choose the likelihood */
    LALInferenceInitLikelihood(runState);

    /* Draw starting positions */
    LALInferenceDrawThreads(runState);

    if (runState == NULL)
        return XLAL_FAILURE;

    /* Call MCMC algorithm */
    if (mpirank == 0) printf("sampling...\n");
    runState->algorithm(runState);

    if (mpirank == 0) printf(" ========== main(): finished. ==========\n");
    MPI_Finalize();

    return XLAL_SUCCESS;
}
