/*
 *  LALInferenceKombine.c:  Bayesian Followup function testing site
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
#include "LALInferenceKombineSampler.h"
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceProposal.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceReadData.h>
#include <lal/LALInferenceInit.h>
#include <lal/LALInferenceCalibrationErrors.h>

#include <mpi.h>

#ifndef _OPENMP
#define omp ignore
#endif


INT4 init_mpi_randomstate(LALInferenceRunState *run_state);
INT4 on_your_marks(LALInferenceRunState *run_state);
INT4 sample_prior(LALInferenceRunState *run_state);
INT4 init_ensemble(LALInferenceRunState *run_state);

/* Set the starting seed of rank 0, and give the rest of the threads
    a seed based on it.  This is a cut-and-paste from PTMCMC, which needs
    to be unified in LALInference when MPI-enabled.*/
INT4 init_mpi_randomstate(LALInferenceRunState *run_state) {
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

     return XLAL_SUCCESS;
}

/* Initialize ensemble randomly or from file */
INT4 on_your_marks(LALInferenceRunState *run_state) {
    LALInferenceVariableItem *item;
    ProcessParamsTable *ppt = NULL;
    REAL8 *sampleArray = NULL;
    INT4 i=0;

    ProcessParamsTable *command_line = run_state->commandLine;

    LALInferenceThreadState *thread = run_state->threads[0];
    LALInferenceVariables *current_param = thread->currentParams;
    INT4 ndim = LALInferenceGetVariableDimensionNonFixed(thread->currentParams);

    /* Determine number of MPI threads, and the
     *   number of walkers run by each thread */
    INT4 nwalkers_per_thread = run_state->nthreads;

    INT4 mpi_rank, walker;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    ppt = LALInferenceGetProcParamVal(command_line, "--init-samples");
    if (ppt) {
        char *infile = ppt->value;
        FILE *input = fopen(infile, "r");

        char params[128][VARNAME_MAX];
        INT4 *col_order = XLALCalloc(ndim, sizeof(INT4));
        INT4 ncols;

        /* Parse parameter names */
        LALInferenceReadAsciiHeader(input, params, &ncols);

        /* Only cluster parameters that are being sampled */
        INT4 nvalid_cols=0, j=0;
        INT4 *valid_cols = XLALCalloc(ncols, sizeof(INT4));

        for (j = 0; j < ncols; j++) {
            char* internal_pname = XLALCalloc(512, sizeof(char));
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

        INT4 nlines = (INT4) nwalkers_per_thread;
        sampleArray = LALInferenceParseDelimitedAscii(input,
                                                        ncols, valid_cols,
                                                        &nlines);


        REAL8 *parameters = XLALCalloc(ndim, sizeof(REAL8));
        for (walker=0; walker<run_state->nthreads; walker++) {
            thread = run_state->threads[walker];

            for (i = 0; i < ndim; i++)
                parameters[i] = sampleArray[mpi_rank*nwalkers_per_thread*ndim + walker*ndim + col_order[i]];

            LALInferenceCopyArrayToVariables(parameters, thread->currentParams);
        }

        /* Cleanup */
        XLALFree(col_order);
        XLALFree(valid_cols);
        XLALFree(parameters);
        XLALFree(sampleArray);
    } else
        LALInferenceDrawThreads(run_state);

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
        thread = run_state->threads[walker];

        thread->currentLikelihood = run_state->likelihood(thread->currentParams,
                                    run_state->data,
                                    thread->model);

    return XLAL_SUCCESS;
}

/********** Initialise MCMC structures *********/

INT4 init_ensemble(LALInferenceRunState *run_state) {
    char help[]="\
    ----------------------------------------------\n\
    --- General Algorithm Parameters -------------\n\
    ----------------------------------------------\n\
    (--nwalkers n)        Number of MCMC walkers to sample with (1000).\n\
    (--nsteps n)          Total number of steps for all walkers to make (10000).\n\
    (--skip n)            Number of steps between writing samples to file (100).\n\
    (--update-interval n) Number of steps between ensemble updates (100).\n\
    (--randomseed seed)   Random seed of sampling distribution (random).\n\
    \n\
    ----------------------------------------------\n\
    --- Output -----------------------------------\n\
    ----------------------------------------------\n\
    (--data-dump)         Output waveforms to file.\n\
    (--outfile file)      Write output files <file>.<chain_number>\n\
                            (ensemble.output.<random_seed>.<mpi_thread>).\n";
    INT4 mpi_rank, mpi_size, i;
    ProcessParamsTable *command_line=NULL, *ppt=NULL;
    LALInferenceThreadState *thread;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    /* Print command line arguments if help requested */
    if (run_state == NULL || LALInferenceGetProcParamVal(run_state->commandLine, "--help")) {
        if (mpi_rank == 0)
            printf("%s", help);

        /* Print other help messages */
        LALInferenceInitCBCThreads(NULL, 0);

        return XLAL_FAILURE;
    }

    command_line = run_state->commandLine;

    /* Set up the appropriate functions for the MCMC algorithm */
    run_state->algorithm = &ensemble_sampler;

    /* Determine number of walkers */
    INT4 nwalkers = 1000;
    INT4 nwalkers_per_thread = nwalkers;
    ppt = LALInferenceGetProcParamVal(command_line, "--nwalkers");
    if (ppt) {
        nwalkers = atoi(ppt->value);
        nwalkers_per_thread = nwalkers / mpi_size;

        if (nwalkers % mpi_size != 0.0) {
            /* Round up to have consistent number of walkers across MPI threads */
            nwalkers_per_thread = (INT4)ceil((REAL8)nwalkers / (REAL8)mpi_size);
            nwalkers = mpi_size * nwalkers_per_thread;

            if (mpi_rank == 0)
                printf("Rounding up number of walkers to %i to provide \
                        consistent performance across the %i available \
                        MPI threads.\n", nwalkers, mpi_size);
        }
    }

    /* Step counter */
    INT4 step = 0;

    /* Number of steps between ensemble updates */
    INT4 nsteps = 10000;
    ppt = LALInferenceGetProcParamVal(command_line, "--nsteps");
    if (ppt)
        nsteps = atoi(ppt->value);

    /* Print sample every skip iterations */
    INT4 skip = 100;
    ppt = LALInferenceGetProcParamVal(command_line, "--skip");
    if (ppt)
        skip = atoi(ppt->value);

    /* Update ensemble every *update_interval* iterations */
    INT4 update_interval = 1000;
    ppt = LALInferenceGetProcParamVal(command_line, "--update-interval");
    if (ppt)
        update_interval = atoi(ppt->value);

    /* Impose cyclic/reflective bounds in KDE */
    INT4 cyclic_reflective = 0;
    if (LALInferenceGetProcParamVal(command_line, "--cyclic-reflective-kde"))
        cyclic_reflective = 1;

    /* Print more stuff */
    INT4 verbose = 0;
    if (LALInferenceGetProcParamVal(command_line, "--verbose"))
        verbose = 1;

    /* Keep track of time if benchmarking */
    INT4 benchmark = 0;
    if (LALInferenceGetProcParamVal(command_line, "--benchmark"))
        benchmark = 1;

    /* Save everything in the run state */
    LALInferenceVariables *algorithm_params = run_state->algorithmParams;

    LALInferenceAddINT4Variable(algorithm_params, "step",
                                step, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddINT4Variable(algorithm_params, "verbose",
                                verbose, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddINT4Variable(algorithm_params, "benchmark",
                                benchmark, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddINT4Variable(algorithm_params, "mpirank",
                                mpi_rank, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddINT4Variable(algorithm_params, "cyclic_reflective",
                                cyclic_reflective, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddINT4Variable(algorithm_params, "nwalkers_per_thread",
                                nwalkers_per_thread, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddINT4Variable(algorithm_params, "nwalkers",
                                nwalkers, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddINT4Variable(algorithm_params, "nsteps",
                                nsteps, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddINT4Variable(algorithm_params, "skip",
                                skip, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceAddINT4Variable(algorithm_params, "update_interval",
                                update_interval, LALINFERENCE_PARAM_OUTPUT);

    /* Initialize the walkers on this MPI thread */
    LALInferenceInitCBCThreads(run_state, nwalkers_per_thread);

    /* Establish the random state across MPI threads */
    init_mpi_randomstate(run_state);

    for (i=0; i<run_state->nthreads; i++) {
        thread = run_state->threads[i];

        thread->id = mpi_rank*nwalkers_per_thread + i;
        thread->proposalArgs = LALInferenceParseProposalArgs(run_state);
    }

    return XLAL_SUCCESS;
}


/* Sample the prior. Useful for defining an initial state for the ensemble */
INT4 sample_prior(LALInferenceRunState *run_state) {
    INT4 update_interval, nprior_steps, nsteps;

    LALInferenceVariables *algorithm_params = run_state->algorithmParams;

    /* Store old algorithm parameters for later restoration */
    nsteps = LALInferenceGetINT4Variable(algorithm_params, "nsteps");
    update_interval = LALInferenceGetINT4Variable(algorithm_params, "update_interval");

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

    return XLAL_SUCCESS;
}

int main(int argc, char *argv[]){
    INT4 mpi_rank;
    ProcessParamsTable *proc_params;
    LALInferenceRunState *run_state = NULL;
    LALInferenceIFOData *data = NULL;

    /* Initialize MPI parallelization */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    if (mpi_rank == 0)
        printf(" ========== lalinference_kombine ==========\n");

    /* Read command line and parse */
    proc_params = LALInferenceParseCommandLine(argc, argv);

    /* Initialise run_state based on command line. This includes allocating
     *   memory, reading data, and performing any injections specified. */
    run_state = LALInferenceInitRunState(proc_params);

    /* Build the ensemble based on command line args */
    init_ensemble(run_state);

    if (run_state == NULL) {
        if (!LALInferenceGetProcParamVal(proc_params, "--help")) {
            fprintf(stderr, "run_state not allocated (%s, line %d).\n",
                    __FILE__, __LINE__);
        }
    } else {
        data = run_state->data;
    }

    /* Perform injections if data successful read or created */
    if (run_state)
        LALInferenceInjectInspiralSignal(data, run_state->commandLine);

    /* Simulate calibration errors. 
    * NOTE: this must be called after both ReadData and (if relevant) 
    * injectInspiralTD/FD are called! */
    LALInferenceApplyCalibrationErrors(data, proc_params);

    /* Choose the prior */
    LALInferenceInitCBCPrior(run_state);

    /* Choose the likelihood */
    LALInferenceInitLikelihood(run_state);

    if (run_state == NULL)
        return XLAL_FAILURE;

    /* Setup the initial state of the walkers */
    on_your_marks(run_state);

    /* Run the sampler to completion */
    run_state->algorithm(run_state);

    if (mpi_rank == 0)
        printf(" ==========  sampling complete ==========\n");

    /* Close down MPI parallelization and return */
    MPI_Finalize();

    return XLAL_SUCCESS;
}
