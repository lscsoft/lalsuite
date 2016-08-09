/*
 *  LALInferenceMCMC.c:  Bayesian Followup, MCMC algorithm.
 *
 *  Copyright (C) 2009, 2012 Ilya Mandel, Vivien Raymond, Christian
 *  Roever, Marc van der Sluys, John Veitch, Will M. Farr, and Ben Farr
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

#include <config.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <lal/LALInspiral.h>
#include <lal/DetResponse.h>
#include <lal/SeqFactories.h>
#include <lal/Date.h>
#include <lal/VectorOps.h>
#include <lal/TimeFreqFFT.h>
#include <lal/GenerateInspiral.h>
#include <lal/TimeDelay.h>
#include <lalapps.h>
#include <mpi.h>
#include <lal/LALInference.h>
#include "LALInferenceMCMCSampler.h"
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferenceProposal.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/LALInferenceReadData.h>
#include <lal/LALInferenceHDF5.h>
#include <sys/time.h>

#include <LALAppsVCSInfo.h>
#include <lal/LALStdlib.h>

#define PROGRAM_NAME "LALInferenceMCMCSampler.c"
#define CVS_ID_STRING "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define CVS_NAME_STRING "$Name$"

const char *const clusteredKDEProposalName = "ClusteredKDEProposal";

static void
thinDifferentialEvolutionPoints(LALInferenceThreadState *thread) {
    size_t i;
    size_t newSize;

    /* Delete all the even-index points. */
    for (i = 0; i < thread->differentialPointsLength; i += 2) {
        LALInferenceClearVariables(thread->differentialPoints[i]);
        XLALFree(thread->differentialPoints[i]);
        thread->differentialPoints[i] = NULL;
    }

    /* Copy the odd points into the first part of the array. */
    for (i = 1; i < thread->differentialPointsLength; i += 2) {
        thread->differentialPoints[i/2] = thread->differentialPoints[i];
        thread->differentialPoints[i] = NULL;
    }

    newSize = thread->differentialPointsLength / 2;

    /* Now shrink the buffer down. */
    thread->differentialPoints = XLALRealloc(thread->differentialPoints, 2*newSize*sizeof(LALInferenceVariables *));
    thread->differentialPointsSize = 2*newSize;
    thread->differentialPointsLength = newSize;
    thread->differentialPointsSkip *= 2;
}

static void
accumulateDifferentialEvolutionSample(LALInferenceThreadState *thread, size_t buffer_limit) {
    if (thread->differentialPointsSize == thread->differentialPointsLength) {
        size_t newSize = thread->differentialPointsSize*2;

        if (buffer_limit < newSize) {
            /* Then thin, and record sample. */
            thinDifferentialEvolutionPoints(thread);
            return accumulateDifferentialEvolutionSample(thread, buffer_limit);
        } else {
            thread->differentialPoints = XLALRealloc(thread->differentialPoints, newSize*sizeof(LALInferenceVariables *));
            thread->differentialPointsSize = newSize;
        }
    }

    thread->differentialPoints[thread->differentialPointsLength] = XLALCalloc(1, sizeof(LALInferenceVariables));
    LALInferenceCopyVariables(thread->currentParams, thread->differentialPoints[thread->differentialPointsLength]);

    thread->differentialPointsLength += 1;
}

static void
resetDifferentialEvolutionBuffer(LALInferenceThreadState *thread) {
    size_t i;

    for (i = 0; i < thread->differentialPointsLength; i++) {
        LALInferenceClearVariables(thread->differentialPoints[i]);
        XLALFree(thread->differentialPoints[i]);
        thread->differentialPoints[i] = NULL;
    }

    thread->differentialPoints = XLALRealloc(thread->differentialPoints, sizeof(LALInferenceVariables *));
    thread->differentialPointsLength = 0;
    thread->differentialPointsSize = 1;
    thread->differentialPointsSkip = LALInferenceGetINT4Variable(thread->proposalArgs, "de_skip");
}


void PTMCMCAlgorithm(struct tagLALInferenceRunState *runState) {
    INT4 t=0; //indexes for for() loops
    INT4 runComplete = 0;
    REAL8 timestamp_epoch=0.0;
    INT4 MPIrank, MPIsize;
    LALStatus status;
    LALInferenceVariables *algorithm_params;
    LALInferenceThreadState *thread;
    INT4 write_interval = 1;
    INT4 step_last_written;

    memset(&status, 0, sizeof(status));

    MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

    algorithm_params = runState->algorithmParams;

    /* Get all algorithm params */
    INT4 n_local_threads = runState->nthreads;
    INT4 nPar = LALInferenceGetVariableDimensionNonFixed(runState->threads[0]->currentParams);
    INT4 Niter = LALInferenceGetINT4Variable(algorithm_params, "nsteps");
    INT4 Neff = LALInferenceGetINT4Variable(algorithm_params, "neff");
    INT4 Nskip = LALInferenceGetINT4Variable(algorithm_params, "skip");
    INT4 temp_skip = LALInferenceGetINT4Variable(algorithm_params, "tskip");
    INT4 de_buffer_limit = LALInferenceGetINT4Variable(algorithm_params, "de_buffer_limit");
    INT4 randomseed = LALInferenceGetINT4Variable(algorithm_params, "random_seed");

    INT4 verbose = LALInferenceGetINT4Variable(algorithm_params, "verbose");
    INT4 propVerbose = LALInferenceGetINT4Variable(algorithm_params, "prop_verbose");
    INT4 propTrack = LALInferenceGetINT4Variable(algorithm_params, "prop_track");
    INT4 tempVerbose = LALInferenceGetINT4Variable(algorithm_params, "temp_verbose");
    INT4 adaptVerbose = LALInferenceGetINT4Variable(algorithm_params, "adapt_verbose");
    INT4 benchmark = LALInferenceGetINT4Variable(algorithm_params, "benchmark");

    /* Clustered-KDE proposal updates */
    INT4 kde_update_start = 200;  // rough number of effective samples to start KDE updates

    /* proposal will be updated 5 times per decade, so this interval will change */
    INT4 *kde_update_interval = XLALCalloc(n_local_threads, sizeof(INT4));
    INT4 *last_kde_update = XLALCalloc(n_local_threads, sizeof(INT4)); // effective sample size at last KDE update

    INT4 diffEvo = 1;
    if (runState->threads[0]->differentialPoints == NULL)
        diffEvo = 0;

    /* Adaptation settings */
    INT4 no_adapt = LALInferenceGetINT4Variable(runState->algorithmParams, "no_adapt");
    INT4 adaptTau = LALInferenceGetINT4Variable(runState->algorithmParams, "adaptTau");
    INT4 adaptLength = LALInferenceGetINT4Variable(runState->algorithmParams, "adaptLength"); // Number of iterations to adapt before turning off

    /* File outputs */
    FILE *verbose_file = NULL;
    char verbose_filename[256];

    if (tempVerbose) {
        sprintf(verbose_filename, "PTMCMC.tempswaps.%u.%2.2d", randomseed, MPIrank);
        verbose_file = fopen(verbose_filename, "w");

        fprintf(verbose_file,
            "cycle\tlow_temp\thigh_temp\tlog(chain_swap)\tlow_temp_likelihood\thigh_temp_likelihood\tswap_accepted\n");

        fclose(verbose_file);
    }

    for (t = 0; t < n_local_threads; t++) {
        thread = runState->threads[t];

        if (adaptVerbose & !no_adapt) {
            sprintf(verbose_filename, "PTMCMC.statistics.%u.%2.2d",
                    randomseed, n_local_threads*MPIrank+t);
            verbose_file = fopen(verbose_filename, "w");

            fprintf(verbose_file, "cycle\ts_gamma");
            LALInferencePrintAdaptationHeader(verbose_file, thread);
            fclose(verbose_file);
        }

        if (propVerbose) {
            sprintf(verbose_filename, "PTMCMC.propstats.%u.%2.2d",
                    randomseed, n_local_threads*MPIrank+t);
            verbose_file = fopen(verbose_filename, "w");

            fprintf(verbose_file, "cycle\t");
            LALInferencePrintProposalStatsHeader(verbose_file, thread->cycle);
            fclose(verbose_file);
        }

        if (propTrack) {
            sprintf(verbose_filename, "PTMCMC.proptrack.%u.%2.2d",
                    randomseed, n_local_threads*MPIrank+t);
            verbose_file = fopen(verbose_filename, "w");

            fprintf(verbose_file, "cycle\t");
            LALInferencePrintProposalTrackingHeader(verbose_file, thread->currentParams);
            fclose(verbose_file);
        }
    }

    for (t = 0; t < n_local_threads; t++)
        record_likelihoods(runState->threads[t]);

    LALInferenceNameOutputs(runState);
    LALInferenceResumeMCMC(runState);
    if (MPIrank == 0)
        LALInferencePrintInjectionSample(runState);

    if (benchmark) {
        struct timeval start_tv;
        gettimeofday(&start_tv, NULL);
        timestamp_epoch = start_tv.tv_sec + start_tv.tv_usec/1E6;
        LALInferenceAddVariable(runState->algorithmParams,
                                "timestamp_epoch", &timestamp_epoch,
                                LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_OUTPUT);
    }

    /* Print run details */
    if (MPIrank == 0) {
        thread = runState->threads[0];
        if (verbose) {
            printf("\nParallel Behavior:\n");
            if (!no_adapt)
                printf(" Adapting with decay power %i for %i iterations after max log(L) increases by nParams/2 (%1.2f).\n", adaptTau, adaptLength, (double)nPar/2.0);
            else
                printf(" Adaptation off.\n");
            if (Neff != Niter)
                printf(" Collecting %i effective samples.\n", Neff);

            printf("\nPTMCMCAlgorithm(); starting parameter values:\n");
            LALInferencePrintVariables(thread->currentParams);
            printf(" MCMC iteration: 0\t");
            printf("%f\t", thread->currentLikelihood - thread->nullLikelihood);
            printf("\n");
        }

        /* Print to file the contents of model->freqhPlus. */
        if (LALInferenceGetProcParamVal(runState->commandLine, "--data-dump"))
            LALInferenceDataDump(runState->data, thread->model);
    }

    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    // iterate:
    step_last_written = runState->threads[0]->step;
    while (!runComplete) {
        #pragma omp parallel for private(thread)
        for (t = 0; t < n_local_threads; t++) {
            FILE *outfile = NULL;
            char outfilename[256];
            //struct timeval tv;
            //REAL8 timestamp=-1.0;
            INT4 i=0;

            thread = runState->threads[t];

            for (i=0; i<temp_skip; i++) {
                /* Increment iteration counter */
                thread->step += 1;

                INT4 adapting = LALInferenceGetINT4Variable(thread->proposalArgs, "adapting");

                if (!no_adapt)
                    LALInferenceAdaptation(thread);

                //ACL calculation
                thread->effective_sample_size = 0;
                if (thread->step % (100*Nskip) == 0) {
                    if (adapting)
                        thread->effective_sample_size = 0;
                    else {
                        thread->effective_sample_size = LALInferenceComputeEffectiveSampleSize(thread);
                        if (verbose && thread->temperature == 1.)
                            printf("Cold thread has collected %i samples.\n", thread->effective_sample_size);
                    }
                }

                if (MPIrank == 0 && t == 0 && thread->effective_sample_size > Neff) {
                    fprintf(stdout,"Thread %i has %i effective samples. Stopping...\n", MPIrank, thread->effective_sample_size);
                    runComplete = 1;          // Sampling is done!
                }

                mcmc_step(runState, thread); //evolve the chain at temperature ladder[t]
                record_likelihoods(thread);

                if (propVerbose)
                    LALInferenceTrackProposalAcceptance(thread);

                if ((thread->step % Nskip) == 0) {
                    /* Update clustered-KDE proposal every time the buffer is expanded */
                    if (LALInferenceGetProcParamVal(runState->commandLine, "--proposal-kde")
                        && (thread->effective_sample_size > kde_update_start)
                        && (((thread->effective_sample_size - last_kde_update[t]) > kde_update_interval[t]) ||
                          ((last_kde_update[t] - thread->effective_sample_size) > kde_update_interval[t]))) {
                        LALInferenceSetupClusteredKDEProposalFromDEBuffer(thread);

                        /* Update 5 times each decade.  This keeps hot chains (with lower ACLs) under control */
                        kde_update_interval[t] = 2 * ((INT4) pow(10.0, floor(log10((REAL8) thread->effective_sample_size))));

                        last_kde_update[t] = thread->effective_sample_size;
                    }

                    if (diffEvo && (thread->step % thread->differentialPointsSkip == 0))
                        accumulateDifferentialEvolutionSample(thread, de_buffer_limit);
                    /*
                    if (benchmark) {
                        gettimeofday(&tv, NULL);
                        timestamp = tv.tv_sec + tv.tv_usec/1E6 - timestamp_epoch;
                    }*/

                    //LALInferenceSaveSample(thread, resumeoutputs[t]);
                    //LALInferencePrintMCMCSample(thread, runState->data, thread->step, timestamp, threadoutputs[t]);
                    LALInferenceLogSampleToArray(thread->algorithmParams, thread->currentParams);

                    if (adaptVerbose && !no_adapt) {
                        sprintf(outfilename, "PTMCMC.statistics.%u.%2.2d",
                                randomseed, n_local_threads*MPIrank+t);
                        outfile = fopen(outfilename, "a");
                        fprintf(outfile, "%d\t", thread->step);
                        LALInferencePrintAdaptationSettings(outfile, thread);
                        fclose(outfile);
                    }

                    if (propVerbose){
                        sprintf(outfilename, "PTMCMC.propstats.%u.%2.2d", randomseed,
                                n_local_threads*MPIrank+t);
                        outfile = fopen(outfilename, "a");
                        fprintf(outfile, "%d\t", thread->step);
                        LALInferencePrintProposalStats(outfile, thread->cycle);
                        fclose(outfile);
                    }

                    if (propTrack) {
                        REAL8 logProposalRatio = LALInferenceGetREAL8Variable(thread->proposalArgs, "logProposalRatio");
                        sprintf(outfilename, "PTMCMC.proptrack.%u.%2.2d", randomseed,
                                n_local_threads*MPIrank+t);
                        outfile = fopen(outfilename, "w");
                        fprintf(outfile, "%d\t", thread->step);
                        LALInferencePrintProposalTracking(outfile, thread->cycle, thread->preProposalParams, thread->proposedParams, logProposalRatio, thread->accepted);
                        fclose(outfile);
                    }
                }
            }
        }

        write_interval = (INT4) pow(10, floor(log10(runState->threads[0]->step + adaptLength)) - 1);
        write_interval = write_interval > 1 ? write_interval : 1;
        if ((runState->threads[0]->step < step_last_written) ||
            (runState->threads[0]->step - step_last_written > write_interval)) {
            LALInferenceCheckpointMCMC(runState);
            LALInferenceWriteMCMCSamples(runState);
            step_last_written = runState->threads[0]->step;
        }

        /* Open swap file if going verbose */
        verbose_file = NULL;
        if (tempVerbose) {
            sprintf(verbose_filename, "PTMCMC.tempswaps.%u.%2.2d", randomseed, MPIrank);
            verbose_file = fopen(verbose_filename, "a");
        }

        /* Excute swap proposal. */
        runState->parallelSwap(runState, verbose_file);

        if (tempVerbose)
            fclose(verbose_file);

        /* Check if run should end */
        if (runState->threads[0]->step > Niter)
            runComplete=1;

        /* Broadcast the root's decision on run completion */
        MPI_Bcast(&runComplete, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }// while (!runComplete)
}

void record_likelihoods(LALInferenceThreadState *thread) {
    REAL8 deltalogl = thread->currentLikelihood - thread->nullLikelihood;
    REAL8 logpost = deltalogl + thread->currentPrior;

    LALInferenceAddVariable(thread->currentParams, "cycle", &(thread->step), LALINFERENCE_INT4_t, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddVariable(thread->currentParams, "logpost", &logpost, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddVariable(thread->currentParams, "logprior", &(thread->currentPrior), LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddVariable(thread->currentParams, "logl", &(thread->currentLikelihood), LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_OUTPUT);
    LALInferenceAddVariable(thread->currentParams, "deltalogl", &deltalogl, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceIFOData *headIFO = thread->parent->data;
    char name[256];
    INT4 ifo = 0;
    while (headIFO != NULL) {
        sprintf(name, "logl%s", headIFO->name);
        REAL8 ifo_logl = thread->currentIFOLikelihoods[ifo] - headIFO->nullloglikelihood;

        LALInferenceAddVariable(thread->currentParams, name, &ifo_logl, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_OUTPUT);
        headIFO = headIFO->next;
        ifo++;
    }
}

void mcmc_step(LALInferenceRunState *runState, LALInferenceThreadState *thread) {
    // Metropolis-Hastings sampler.
    INT4 MPIrank;
    REAL8 logPriorCurrent, logPriorProposed;
    REAL8 logLikelihoodCurrent, logLikelihoodProposed;
    REAL8 logProposalRatio = 0.0;  // = log(P(backward)/P(forward))
    REAL8 logAcceptanceProbability;
    REAL8 targetAcceptance = 0.234;

    MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

    INT4 outputSNRs = LALInferenceGetINT4Variable(runState->algorithmParams, "output_snrs");
    INT4 propTrack = LALInferenceGetINT4Variable(runState->algorithmParams, "prop_track");

    // current values:
    logPriorCurrent = thread->currentPrior;
    logLikelihoodCurrent = thread->currentLikelihood;

    // generate proposal:
    logProposalRatio = thread->proposal(thread, thread->currentParams, thread->proposedParams);

    // compute prior & likelihood:
    logPriorProposed = runState->prior(runState, thread->proposedParams, thread->model);
    if (logPriorProposed > -DBL_MAX)
        logLikelihoodProposed = runState->likelihood(thread->proposedParams, runState->data, thread->model);
    else
        logLikelihoodProposed = -DBL_MAX;

    if (propTrack)
        LALInferenceCopyVariables(thread->currentParams, thread->preProposalParams);

    // determine acceptance probability:
    logAcceptanceProbability = (1.0/thread->temperature)*(logLikelihoodProposed - logLikelihoodCurrent)
                                + (logPriorProposed - logPriorCurrent)
                                + logProposalRatio;

    // accept/reject:
    if ((logAcceptanceProbability > 0)
        || (log(gsl_rng_uniform(thread->GSLrandom)) < logAcceptanceProbability)) {   //accept
        LALInferenceCopyVariables(thread->proposedParams, thread->currentParams);
        thread->currentLikelihood = logLikelihoodProposed;
        thread->currentPrior = logPriorProposed;

        /* Calculate SNR if requested, and not already calculated by prior */
        LALInferenceIFOData *ifoPtr = runState->data;
        INT4 ifo = 0;
        if (outputSNRs) {
            if (thread->model->SNR == 0.0)
                LALInferenceNetworkSNR(thread->currentParams, runState->data, thread->model);
            thread->currentSNR = thread->model->SNR;
            while (ifoPtr) {
                thread->currentIFOSNRs[ifo] = thread->model->ifo_SNRs[ifo];
                ifo++;
                ifoPtr = ifoPtr->next;
            }
        }

        if (thread->model->ifo_loglikelihoods != NULL) {
            ifoPtr = runState->data;
            ifo = 0;
            while (ifoPtr) {
                thread->currentIFOLikelihoods[ifo] = thread->model->ifo_loglikelihoods[ifo];
                ifo++;
                ifoPtr = ifoPtr->next;
            }
        }

        thread->acceptanceCount++;
        thread->accepted = 1;
    } else {
        thread->accepted = 0;
    }

    LALInferenceUpdateAdaptiveJumps(thread, targetAcceptance);

    return;
}


//-----------------------------------------
// Swap routines:
//-----------------------------------------
void LALInferencePTswap(LALInferenceRunState *runState, FILE *swapfile) {
    INT4 MPIrank, MPIsize;
    MPI_Status MPIstatus;
    INT4 nPar, adjNPar, n_local_threads, ntemps;
    INT4 low;
    INT4 ind, cold_ind, hot_ind;
    INT4 cold_rank, hot_rank;
    INT4 swapAccepted;
    INT4 *cold_inds;
    REAL8 adjCurrentLikelihood, adjCurrentPrior;
    REAL8 logThreadSwap, temp_prior, temp_like, cold_temp;
    LALInferenceThreadState *cold_thread = runState->threads[0];
    LALInferenceThreadState *hot_thread;
    LALInferenceVariables *temp_params;

    MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
    MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);

    n_local_threads = runState->nthreads;
    ntemps = MPIsize*n_local_threads;

    /* Return if running with only a single temperature */
    if (ntemps == 1)
        return;

    cold_inds = XLALCalloc(ntemps-1, sizeof(INT4));

    /* Have the root process choose a random order of ladder swaps and share it */
    if (MPIrank == 0) {
        for (low = 0; low < ntemps-1; low++)
            cold_inds[low] = low;

        gsl_ran_shuffle(runState->GSLrandom, cold_inds, ntemps-1, sizeof(INT4));
    }

    MPI_Bcast(cold_inds, ntemps-1, MPI_INT, 0, MPI_COMM_WORLD);

    for (ind=0; ind < ntemps-1; ind++) {
        cold_ind = cold_inds[ind];
        hot_ind = cold_ind+1;

        /* Check if cold chain is handled by this thread */
        cold_rank = cold_ind/n_local_threads;
        hot_rank = hot_ind/n_local_threads;

        if (cold_rank == hot_rank) {
            if (MPIrank == cold_rank) {
                cold_thread = runState->threads[cold_ind % n_local_threads];
                hot_thread = runState->threads[hot_ind % n_local_threads];

                /* Determine if swap is accepted and tell the other chain */
                logThreadSwap = 1.0/cold_thread->temperature - 1.0/hot_thread->temperature;
                logThreadSwap *= hot_thread->currentLikelihood - cold_thread->currentLikelihood;

                if ((logThreadSwap > 0) || (log(gsl_rng_uniform(runState->GSLrandom)) < logThreadSwap ))
                    swapAccepted = 1;
                else
                    swapAccepted = 0;

                /* Print to file if verbose is chosen */
                if (swapfile != NULL)
                    fprintf(swapfile, "%d\t%f\t%f\t%f\t%f\t%f\t%i\n",
                            hot_thread->step, cold_thread->temperature, hot_thread->temperature,
                            logThreadSwap, cold_thread->currentLikelihood,
                            hot_thread->currentLikelihood, swapAccepted);

                if (swapAccepted) {
                    temp_params = hot_thread->currentParams;
                    temp_prior = hot_thread->currentPrior;
                    temp_like = hot_thread->currentLikelihood;

                    hot_thread->currentParams = cold_thread->currentParams;
                    hot_thread->currentPrior = cold_thread->currentPrior;
                    hot_thread->currentLikelihood = cold_thread->currentLikelihood;

                    cold_thread->currentParams = temp_params;
                    cold_thread->currentPrior = temp_prior;
                    cold_thread->currentLikelihood = temp_like;
                }
            }
        } else {
            if (MPIrank == cold_rank) {
                cold_thread = runState->threads[cold_ind % n_local_threads];

                /* Send current likelihood for swap proposal */
                MPI_Send(&(cold_thread->temperature), 1, MPI_DOUBLE, hot_rank, PT_COM, MPI_COMM_WORLD);
                MPI_Send(&(cold_thread->currentLikelihood), 1, MPI_DOUBLE, hot_rank, PT_COM, MPI_COMM_WORLD);

                /* Determine if swap was accepted */
                MPI_Recv(&swapAccepted, 1, MPI_INT, hot_rank, PT_COM, MPI_COMM_WORLD, &MPIstatus);

                /* Perform Swap */
                if (swapAccepted) {
                    /* Set new likelihood */
                    MPI_Recv(&adjCurrentLikelihood, 1, MPI_DOUBLE, hot_rank, PT_COM, MPI_COMM_WORLD, &MPIstatus);
                    cold_thread->currentLikelihood = adjCurrentLikelihood;

                    /* Exchange current prior values */
                    MPI_Send(&(cold_thread->currentPrior), 1, MPI_DOUBLE, hot_rank, PT_COM, MPI_COMM_WORLD);
                    MPI_Recv(&adjCurrentPrior, 1, MPI_DOUBLE, hot_rank, 0, MPI_COMM_WORLD, &MPIstatus);
                    cold_thread->currentPrior = adjCurrentPrior;

                    /* Package and send parameters */
                    nPar = LALInferenceGetVariableDimensionNonFixed(cold_thread->currentParams);
                    MPI_Send(&nPar, 1, MPI_INT, hot_rank, PT_COM, MPI_COMM_WORLD);
                    MPI_Recv(&adjNPar, 1, MPI_INT, hot_rank, PT_COM, MPI_COMM_WORLD, &MPIstatus);

                    REAL8 *parameters = XLALMalloc(nPar * sizeof(REAL8));
                    LALInferenceCopyVariablesToArray(cold_thread->currentParams, parameters);
                    MPI_Send(parameters, nPar, MPI_DOUBLE, hot_rank, PT_COM, MPI_COMM_WORLD);

                    /* Recieve and unpack parameters */
                    REAL8 *adjParameters = XLALMalloc(adjNPar * sizeof(REAL8));
                    MPI_Recv(adjParameters, adjNPar, MPI_DOUBLE, hot_rank, PT_COM, MPI_COMM_WORLD, &MPIstatus);
                    LALInferenceCopyArrayToVariables(adjParameters, cold_thread->currentParams);

                    XLALFree(parameters);
                    XLALFree(adjParameters);
                }
            } else if (MPIrank == hot_rank) {
                hot_thread = runState->threads[hot_ind % n_local_threads];

                MPI_Recv(&cold_temp, 1, MPI_DOUBLE, cold_rank, PT_COM, MPI_COMM_WORLD, &MPIstatus);

                /* Receive adjacent likelilhood */
                MPI_Recv(&adjCurrentLikelihood, 1, MPI_DOUBLE, cold_rank, PT_COM, MPI_COMM_WORLD, &MPIstatus);

                /* Determine if swap is accepted and tell the other chain */
                logThreadSwap = 1.0/cold_temp - 1.0/hot_thread->temperature;
                logThreadSwap *= hot_thread->currentLikelihood - adjCurrentLikelihood;
                if ((logThreadSwap > 0) || (log(gsl_rng_uniform(runState->GSLrandom)) < logThreadSwap ))
                    swapAccepted = 1;
                else
                    swapAccepted = 0;

                MPI_Send(&swapAccepted, 1, MPI_INT, cold_rank, PT_COM, MPI_COMM_WORLD);

                /* Print to file if verbose is chosen */
                if (swapfile != NULL) {
                    fprintf(swapfile, "%d%f\t%f\t\t%f\t%f\t%f\t%i\n",
                            hot_thread->step, cold_temp, hot_thread->temperature,
                            logThreadSwap, adjCurrentLikelihood,
                            hot_thread->currentLikelihood, swapAccepted);
                    fflush(swapfile);
                }

                /* Perform Swap */
                if (swapAccepted) {
                    /* Swap likelihoods */
                    MPI_Send(&(hot_thread->currentLikelihood), 1, MPI_DOUBLE, cold_rank, PT_COM, MPI_COMM_WORLD);
                    hot_thread->currentLikelihood = adjCurrentLikelihood;

                    /* Exchange current prior values */
                    MPI_Recv(&adjCurrentPrior, 1, MPI_DOUBLE, cold_rank, PT_COM, MPI_COMM_WORLD, &MPIstatus);
                    MPI_Send(&(hot_thread->currentPrior), 1, MPI_DOUBLE, cold_rank, PT_COM, MPI_COMM_WORLD);
                    hot_thread->currentPrior = adjCurrentPrior;

                    /* Package parameters */
                    nPar = LALInferenceGetVariableDimensionNonFixed(hot_thread->currentParams);
                    MPI_Recv(&adjNPar, 1, MPI_INT, cold_rank, PT_COM, MPI_COMM_WORLD, &MPIstatus);
                    MPI_Send(&nPar, 1, MPI_INT, cold_rank, PT_COM, MPI_COMM_WORLD);

                    REAL8 *parameters = XLALMalloc(nPar * sizeof(REAL8));
                    LALInferenceCopyVariablesToArray(hot_thread->currentParams, parameters);

                    /* Swap parameters */
                    REAL8 *adjParameters = XLALMalloc(adjNPar * sizeof(REAL8));
                    MPI_Recv(adjParameters, adjNPar, MPI_DOUBLE, cold_rank, PT_COM, MPI_COMM_WORLD, &MPIstatus);
                    MPI_Send(parameters, nPar, MPI_DOUBLE, cold_rank, PT_COM, MPI_COMM_WORLD);

                    /* Unpack parameters */
                    LALInferenceCopyArrayToVariables(adjParameters, hot_thread->currentParams);

                    XLALFree(parameters);
                    XLALFree(adjParameters);
                }
            }
        }
    }

    XLALFree(cold_inds);

    return;
}


// UINT4 LALInferenceMCMCMCswap(LALInferenceRunState *runState, REAL8 *ladder, INT4 i, FILE *swapfile) {
//     INT4 MPIrank;
//     MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);
//     MPI_Status MPIstatus;
//     INT4 nPar = LALInferenceGetVariableDimensionNonFixed(runState->currentParams);
//     INT4 nChain = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "nChain");
//     INT4 Tskip = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "tempSkip");
//     REAL8 adjCurrentPrior;
//     REAL8 logChainSwap;
//     INT4 readyToSwap = 0;
//     UINT4 swapProposed=0;
//     INT4 swapAccepted=0;
//
//     REAL8 lowLikeLowParams = 0;
//     REAL8 highLikeHighParams = 0;
//     REAL8 lowLikeHighParams = 0;
//     REAL8 highLikeLowParams = 0;
//
//     LALInferenceParamVaryType fLowVary = LALINFERENCE_PARAM_FIXED;
//     LALInferenceVariables *adjCurrentParams = NULL;
//     adjCurrentParams = (LALInferenceVariables *)XLALCalloc(sizeof(LALInferenceVariables), 1);
//
//     /* If Tskip reached, then block until next chain in ladder is prepared to accept swap proposal */
//     if (((i % Tskip) == 0) && MPIrank < nChain-1) {
//         swapProposed = 1;
//         /* Send current likelihood for swap proposal */
//         MPI_Send(&(runState->currentLikelihood), 1, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD);
//
//         /* Fix fLow if it is being varied, to avoid swapping it */
//         if(LALInferenceCheckVariable(runState->currentParams,"fLow")) {
//             fLowVary = LALInferenceGetVariableVaryType(runState->currentParams,"fLow");
//             if(fLowVary != LALINFERENCE_PARAM_FIXED) {
//                 LALInferenceSetParamVaryType(runState->currentParams,"fLow",LALINFERENCE_PARAM_FIXED);
//                 nPar -= 1;
//             }
//         }
//
//         /* Prepare Variables structure to recieve adjacent parameters */
//         LALInferenceCopyVariables(runState->currentParams, adjCurrentParams);
//
//         /* Package, swap, and unpack parameters */
//         REAL8 *parameters = XLALMalloc(nPar * sizeof(REAL8));
//         REAL8 *adjParameters = XLALMalloc(nPar * sizeof(REAL8));
//         LALInferenceCopyVariablesToArray(runState->currentParams, parameters);
//         MPI_Send(parameters, nPar, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD);
//         MPI_Recv(adjParameters, nPar, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
//         LALInferenceCopyArrayToVariables(adjParameters, adjCurrentParams);
//         XLALFree(parameters);
//         XLALFree(adjParameters);
//
//         /* Calculate likelihood at adjacent parameters and send */
//         lowLikeHighParams = runState->likelihood(adjCurrentParams, runState->data, runState->model);
//         MPI_Send(&lowLikeHighParams, 1, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD);
//
//         /* Determine if swap was accepted */
//         MPI_Recv(&swapAccepted, 1, MPI_INT, MPIrank+1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
//
//         if (swapAccepted) {
//             /* Set new likelihood */
//             runState->currentLikelihood = lowLikeHighParams;
//
//             MPI_Send(&(runState->currentPrior), 1, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD);
//             MPI_Recv(&adjCurrentPrior, 1, MPI_DOUBLE, MPIrank+1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
//             runState->currentPrior = adjCurrentPrior;
//
//             /* Set new parameters */
//             LALInferenceCopyVariables(adjCurrentParams,runState->currentParams);
//         }
//
//         /* Unfix fLow if it was originally unfixed */
//         if(fLowVary!=LALINFERENCE_PARAM_FIXED) {
//             LALInferenceSetParamVaryType(runState->currentParams,"fLow",fLowVary);
//             nPar += 1;
//         }
//
//         if (swapAccepted) {
//             /* Calculate prior at new values */
//             runState->currentPrior = runState->prior(runState, runState->currentParams, runState->model);
//         }
//     }
//
//     /* Check if next lower temperature is ready to swap */
//     else if (MPIrank > 0) {
//             MPI_Iprobe(MPIrank-1, PT_COM, MPI_COMM_WORLD, &readyToSwap, &MPIstatus);
//
//         /* Hotter chain decides acceptance */
//         if (readyToSwap) {
//             /* Receive adjacent likelilhood */
//             MPI_Recv(&lowLikeLowParams, 1, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
//             highLikeHighParams = runState->currentLikelihood;
//
//             /* Fix fLow if it is being varied, to avoid swapping it */
//             if(LALInferenceCheckVariable(runState->currentParams,"fLow")) {
//                 fLowVary = LALInferenceGetVariableVaryType(runState->currentParams,"fLow");
//             if(fLowVary != LALINFERENCE_PARAM_FIXED) {
//                 LALInferenceSetParamVaryType(runState->currentParams,"fLow",LALINFERENCE_PARAM_FIXED);
//                 nPar -= 1;
//             }
//         }
//
//         /* Prepare Variables structure to recieve adjacent parameters */
//         LALInferenceCopyVariables(runState->currentParams, adjCurrentParams);
//
//         /* Package, swap, and unpack parameters */
//         REAL8 *parameters = XLALMalloc(nPar * sizeof(REAL8));
//         REAL8 *adjParameters = XLALMalloc(nPar * sizeof(REAL8));
//         LALInferenceCopyVariablesToArray(runState->currentParams, parameters);
//         MPI_Recv(adjParameters, nPar, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
//         MPI_Send(parameters, nPar, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD);
//         LALInferenceCopyArrayToVariables(adjParameters, adjCurrentParams);
//         XLALFree(parameters);
//         XLALFree(adjParameters);
//
//         /* Calculate likelihood at adjacent parameters */
//         highLikeLowParams = runState->likelihood(adjCurrentParams, runState->data, runState->model);
//
//         /* Recieve likelihood from adjacent chain */
//         MPI_Recv(&lowLikeHighParams, 1, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
//
//         /* Propose swap */
//         logChainSwap = (1./ladder[MPIrank-1])*(lowLikeHighParams-lowLikeLowParams)+(1./ladder[MPIrank])*(highLikeLowParams-highLikeHighParams);
//
//         if ((logChainSwap > 0) || (log(gsl_rng_uniform(runState->GSLrandom)) < logChainSwap )) {
//             swapAccepted = 1;
//         } else {
//             swapAccepted = 0;
//         }
//         MPI_Send(&swapAccepted, 1, MPI_INT, MPIrank-1, PT_COM, MPI_COMM_WORLD);
//
//         /* Print to file if verbose is chosen */
//         if (swapfile != NULL) {
//             fprintf(swapfile,"%d\t%f\t%f\t%f\t%i\n",i,logChainSwap,lowLikeLowParams,highLikeHighParams,swapAccepted);
//             fflush(swapfile);
//         }
//
//         if (swapAccepted) {
//             /* Swap likelihoods */
//             runState->currentLikelihood = highLikeLowParams;
//
//             /* Exchange current prior values (assumes prior functions are identical) */
//             MPI_Recv(&adjCurrentPrior, 1, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD, &MPIstatus);
//             MPI_Send(&(runState->currentPrior), 1, MPI_DOUBLE, MPIrank-1, PT_COM, MPI_COMM_WORLD);
//             runState->currentPrior = adjCurrentPrior;
//
//             /* Set new parameters */
//             LALInferenceCopyVariables(adjCurrentParams,runState->currentParams);
//         }
//
//         /* Unfix fLow if it was originally unfixed */
//         if(fLowVary!=LALINFERENCE_PARAM_FIXED) {
//             LALInferenceSetParamVaryType(runState->currentParams,"fLow",fLowVary);
//             nPar += 1;
//         }
//
//         if (swapAccepted) {
//             /* Calculate prior at new values */
//             runState->currentPrior = runState->prior(runState, runState->currentParams, runState->model);
//             }
//         }
//     }
//
//     /* Return values for colder chain: 0=nothing happened; 1=swap proposed, not accepted; 2=swap proposed & accepted */
//     if (swapProposed && swapAccepted)
//         swapProposed++;
//
//     LALInferenceClearVariables(adjCurrentParams);
//     XLALFree(adjCurrentParams);
//     return swapProposed;
// }


//-----------------------------------------
// Adaptation:
//-----------------------------------------
void LALInferenceAdaptation(LALInferenceThreadState *thread) {
    REAL8 s_gamma;

    INT4 nPar = LALInferenceGetVariableDimensionNonFixed(thread->currentParams);

    /* if maximum logL has increased by more than nParam/2, restart it */
    INT4 adapting = LALInferenceGetINT4Variable(thread->proposalArgs, "adapting");
    INT4 adaptTau = LALInferenceGetINT4Variable(thread->proposalArgs, "adaptTau");
    INT4 length = LALInferenceGetINT4Variable(thread->proposalArgs, "adaptLength");
    INT4 adaptRestartBuffer = LALInferenceGetINT4Variable(thread->proposalArgs, "adaptResetBuffer");
    REAL8 logLAtAdaptStart = LALInferenceGetREAL8Variable(thread->proposalArgs, "logLAtAdaptStart");

    if (thread->currentLikelihood > logLAtAdaptStart+(REAL8)nPar/2) {
        LALInferenceAdaptationRestart(thread);
    } else if (adapting) {
        /* Turn off adaption after adaptLength steps without restarting */
        if (thread->step > 0) {
            adapting = 0;  //turn off adaptation
            LALInferenceSetVariable(thread->proposalArgs, "adapting", &adapting);
            LALInferenceRemoveVariable(thread->proposalArgs, "s_gamma");

            /* Clear differential buffer so that it contains only post-burnin samples */
            resetDifferentialEvolutionBuffer(thread);

        /* Else set adaptation envelope */
        } else {
            s_gamma = LALInferenceAdaptationEnvelope(thread->step, adaptTau, length, adaptRestartBuffer);
            LALInferenceSetVariable(thread->proposalArgs, "s_gamma", &s_gamma);
        }
    }
}


//-----------------------------------------
// Restart adaptation:
//-----------------------------------------
void LALInferenceAdaptationRestart(LALInferenceThreadState *thread) {
    INT4 adapting=1;
    LALInferenceVariableItem *item;
    INT4 bigACL = INT_MAX;

    for (item=thread->currentParams->head; item; item=item->next) {
        if (LALInferenceCheckVariableNonFixed(thread->currentParams, item->name) &&
                item->type == LALINFERENCE_REAL8_t) {
            char tmpname[MAX_STRLEN]="";

            sprintf(tmpname, "%s_%s", item->name, ACCEPTSUFFIX);
            REAL8 *accepted = (REAL8 *)LALInferenceGetVariable(thread->proposalArgs, tmpname);
            *accepted = 0;

            sprintf(tmpname,"%s_%s", item->name, PROPOSEDSUFFIX);
            REAL8 *proposed = (REAL8 *)LALInferenceGetVariable(thread->proposalArgs, tmpname);
            *proposed = 0;
        }
    }

    INT4 length = LALInferenceGetINT4Variable(thread->proposalArgs, "adaptLength");
    INT4 tau = LALInferenceGetINT4Variable(thread->proposalArgs, "adaptTau");
    INT4 adaptRestartBuffer = LALInferenceGetINT4Variable(thread->proposalArgs, "adaptResetBuffer");
    REAL8 s_gamma = LALInferenceAdaptationEnvelope(thread->step, tau, length, adaptRestartBuffer);
    LALInferenceAddVariable(thread->proposalArgs, "s_gamma", &s_gamma, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_OUTPUT);

    LALInferenceSetVariable(thread->proposalArgs, "adapting", &adapting);
    LALInferenceSetVariable(thread->proposalArgs, "logLAtAdaptStart", &(thread->currentLikelihood));
    LALInferenceSetVariable(thread->proposalArgs, "acl", &bigACL);

    thread->step = -length;
}

//-----------------------------------------
// Adaptation envelope function:
//-----------------------------------------
REAL8 LALInferenceAdaptationEnvelope(INT4 step, INT4 tau, INT4 length, INT4 fix_adapt_len) {
    REAL8 s_gamma = 0.0;

    if (step <= -length + fix_adapt_len) {
        s_gamma=((REAL8)(length + step)/(REAL8)fix_adapt_len);
        s_gamma *= s_gamma;
    } else if (step < 0) {
        s_gamma = 10.0*exp(-(1.0/tau)*log((REAL8)(length + step)))-1;
    } else {
        s_gamma = 0.0;
    }

    return s_gamma;
}


//-----------------------------------------
// file output routines:
//-----------------------------------------
/* Decide on output and resume file names and store to run state */
void LALInferenceNameOutputs(LALInferenceRunState *runState) {
    ProcessParamsTable *ppt;
    INT4 randomseed;
    INT4 MPIrank, t, n_local_threads;
    LALInferenceThreadState *thread;

    MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

    n_local_threads = runState->nthreads;
    for (t=0; t<n_local_threads; t++) {
        thread = runState->threads[t];

        if (n_local_threads * MPIrank + t == 0)
            sprintf(thread->name, "posterior_samples");
        else
            sprintf(thread->name, "chain_%2.2d", thread->id);
    }

    randomseed = LALInferenceGetINT4Variable(runState->algorithmParams,"random_seed");

    ppt = LALInferenceGetProcParamVal(runState->commandLine, "--outfile");
    if (ppt) {
        runState->outFileName = (char*)XLALCalloc(strlen(ppt->value)+255, sizeof(char*));
        runState->resumeOutFileName = (char*)XLALCalloc(strlen(ppt->value)+255, sizeof(char*));
        if (MPIrank == 0) {
            sprintf(runState->outFileName, "%s", ppt->value);
            sprintf(runState->resumeOutFileName, "%s.resume", ppt->value);
        } else {
            sprintf(runState->outFileName, "%s.%2.2d", ppt->value, MPIrank);
            sprintf(runState->resumeOutFileName, "%s.%2.2d.resume", ppt->value, MPIrank);
        }

    } else {
        runState->outFileName = (char*)XLALCalloc(255, sizeof(char*));
        runState->resumeOutFileName = (char*)XLALCalloc(255, sizeof(char*));
        if (MPIrank == 0) {
            sprintf(runState->outFileName, "PTMCMC.output.%u.h5", randomseed);
            sprintf(runState->resumeOutFileName, "PTMCMC.output.%u.resume.h5", randomseed);
        } else {
            sprintf(runState->outFileName, "PTMCMC.output.%u.%2.2d.h5", randomseed, MPIrank);
            sprintf(runState->resumeOutFileName, "PTMCMC.output.%u.%2.2d.resume.h5", randomseed, MPIrank);
        }
    }

    if((ppt=LALInferenceGetProcParamVal(runState->commandLine, "--runid")))
        snprintf(runState->runID, 255, "%s_%s", "lalinference_mcmc",ppt->value);
    else
        snprintf(runState->runID, 255, "lalinference_mcmc");
}


/* If running with --resume, check if the resume file exists.  If it does, restore the runstate */
void LALInferenceResumeMCMC(LALInferenceRunState *runState) {
    //INT4 t, n_local_threads;
    //LALInferenceThreadState *thread;

    if (LALInferenceGetProcParamVal(runState->commandLine, "--resume") &&
            access(runState->outFileName, R_OK) == 0 &&
            access(runState->resumeOutFileName, R_OK) == 0) {

        /* Then file already exists for reading, and we're going to resume
        from it, so don't write the header. */
        LALInferenceReadMCMCCheckpoint(runState);
    }

    return;
}

/* Store the MCMC run state to HDF5 for use by --resume */
void LALInferenceCheckpointMCMC(LALInferenceRunState *runState) {
    //ProcessParamsTable *ppt;
    INT4 t, n_local_threads;
    INT4 MPIrank;
    LALH5File *resume_file = NULL;
    LALInferenceThreadState *thread;

    MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

    resume_file = XLALH5FileOpen(runState->resumeOutFileName, "w");
    if(resume_file == NULL){
        XLALErrorHandler = XLALExitErrorHandler;
        XLALPrintError("Output file error. Please check that the specified path exists. (in %s, line %d)\n",__FILE__, __LINE__);
        XLAL_ERROR_VOID(XLAL_EIO);
    }

    LALH5File *group = LALInferenceH5CreateGroupStructure(resume_file, "lalinference", runState->runID);

    n_local_threads = runState->nthreads;
    for (t = 0; t < n_local_threads; t++) {
        thread = runState->threads[t];

        char chain_group_name[1024];
        snprintf(chain_group_name, sizeof(chain_group_name), "%s-checkpoint", thread->name);
        LALH5File *chain_group = XLALH5GroupOpen(group, chain_group_name);
        /*
         * FIXME: use GNUism asprintf or fall back to libiberty
         *

        char *chain_group_name = NULL;
        asprintf(&chain_group_name, "%s-checkpoint", thread->name);
        LALH5File *chain_group = XLALH5GroupOpen(group, chain_group_name);
        free(chain_group_name);
        */

        /* Create run identifier group */
        LALInferenceH5VariablesArrayToDataset(chain_group, thread->differentialPoints, thread->differentialPointsLength, "differential_points");
        LALInferenceH5VariablesArrayToDataset(chain_group, &(thread->proposalArgs), 1, "proposal_arguments");
        LALInferenceH5VariablesArrayToDataset(chain_group, &(thread->currentParams), 1, "current_parameters");
        XLALH5FileAddScalarAttribute(chain_group, "last_step", &(thread->step), LAL_D_TYPE_CODE);
        XLALH5FileAddScalarAttribute(chain_group, "effective_sample_size", &(thread->effective_sample_size), LAL_D_TYPE_CODE);
        XLALH5FileAddScalarAttribute(chain_group, "differential_point_skip", &(thread->differentialPointsSkip), LAL_D_TYPE_CODE);

        /* TODO: Write metadata */
        XLALH5FileClose(chain_group);
    }
    XLALH5FileClose(group);

    XLALH5FileClose(resume_file);

    return;
}

/* Read in and restore the run state from an MCMC checkpoint */
void LALInferenceReadMCMCCheckpoint(LALInferenceRunState *runState) {
    //ProcessParamsTable *ppt;
    INT4 t, n_local_threads;
    UINT4 n;
    LALH5File *resume_file = NULL;
    LALH5File *output = NULL;
    LALInferenceThreadState *thread;

    /* Read it the resume file, which stores info needed to restore the proposals (adaptation settings, etc.) */
    resume_file = XLALH5FileOpen(runState->resumeOutFileName, "r");
    if(resume_file == NULL){
        XLALErrorHandler = XLALExitErrorHandler;
        XLALPrintError("Output file error. Please check that the specified path exists. (in %s, line %d)\n",__FILE__, __LINE__);
        XLAL_ERROR_VOID(XLAL_EIO);
    }

    LALH5File *li_group = XLALH5GroupOpen(resume_file, "lalinference");
    LALH5File *group = XLALH5GroupOpen(li_group, runState->runID);

    n_local_threads = runState->nthreads;
    for (t = 0; t < n_local_threads; t++) {
        thread = runState->threads[t];

        char chain_group_name[1024];
        snprintf(chain_group_name, sizeof(chain_group_name), "%s-checkpoint", thread->name);
        LALH5File *chain_group = XLALH5GroupOpen(group, chain_group_name);
        /*
         * FIXME: use GNUism asprintf or fall back to libiberty
         *

        char *chain_group_name = NULL;
        asprintf(&chain_group_name, "%s-checkpoint", thread->name);
        LALH5File *chain_group = XLALH5GroupOpen(group, chain_group_name);
        free(chain_group_name);
        */

        /* Restore differential evolution buffer */
        LALH5Dataset *de_group = XLALH5DatasetRead(chain_group, "differential_points");
        LALInferenceH5DatasetToVariablesArray(de_group, &(thread->differentialPoints), (UINT4 *)&(thread->differentialPointsLength));
        thread->differentialPointsSize = thread->differentialPointsLength;

        /* Restore proposal arguments, most importantly adaptation settings */
        LALInferenceVariables **propArgs;
        LALH5Dataset *prop_arg_group = XLALH5DatasetRead(chain_group, "proposal_arguments");
        LALInferenceH5DatasetToVariablesArray(prop_arg_group, &propArgs, &n);
        LALInferenceCopyVariables(propArgs[0], thread->proposalArgs);

        /* We don't save strings, so we can't tell which parameter was last updated with adaptation.
         * So we'll treat the last step as non-adaptable */
        INT4 adaptable_step = 0;
        LALInferenceSetVariable(thread->proposalArgs, "adaptableStep", &adaptable_step);


        /* Restore the parameters of the last sample */
        LALInferenceVariables **currentParams;
        LALH5Dataset *current_param_group = XLALH5DatasetRead(chain_group, "current_parameters");
        LALInferenceH5DatasetToVariablesArray(current_param_group, &currentParams, &n);
        LALInferenceCopyVariables(currentParams[0], thread->currentParams);

        /* Recalculate the likelihood and prior for the restored parameters */
        thread->currentPrior = runState->prior(runState,
                                               thread->currentParams,
                                               thread->model);

        thread->currentLikelihood = runState ->likelihood(thread->currentParams,
                                                          runState->data,
                                                          thread->model);

        /* Recover the estimated effective sample size, iteration number, and DE buffer thinning multiplier */
        XLALH5FileQueryScalarAttributeValue(&(thread->effective_sample_size), chain_group, "effective_sample_size");
        XLALH5FileQueryScalarAttributeValue(&(thread->step), chain_group, "last_step");
        XLALH5FileQueryScalarAttributeValue(&(thread->differentialPointsSkip), chain_group, "differential_point_skip");

        /* TODO: Write metadata */
        XLALH5DatasetFree(current_param_group);
        XLALH5DatasetFree(prop_arg_group);
        XLALH5DatasetFree(de_group);
        XLALH5FileClose(chain_group);
    }
    XLALH5FileClose(group);
    XLALH5FileClose(li_group);
    XLALH5FileClose(resume_file);

    /* Read in samples collected so far */
    output = XLALH5FileOpen(runState->outFileName, "r");
    if(output == NULL){
        XLALErrorHandler = XLALExitErrorHandler;
        XLALPrintError("Output file error. Please check that the specified path exists. (in %s, line %d)\n",__FILE__, __LINE__);
        XLAL_ERROR_VOID(XLAL_EIO);
    }

    li_group = XLALH5GroupOpen(output, "lalinference");
    group = XLALH5GroupOpen(li_group, runState->runID);

    n_local_threads = runState->nthreads;
    for (t = 0; t < n_local_threads; t++) {
        thread = runState->threads[t];

        LALH5Dataset *chain_dataset = XLALH5DatasetRead(group, thread->name);

        LALInferenceVariables **input_array;
        UINT4 i,N;
        LALInferenceH5DatasetToVariablesArray(chain_dataset, &input_array, &N);
        for (i=0; i<N; i++)
            LALInferenceLogSampleToArray(thread->algorithmParams, input_array[i]);

        /* TODO: Write metadata */
        XLALH5DatasetFree(chain_dataset);
    }
    XLALH5FileClose(group);
    XLALH5FileClose(li_group);
    XLALH5FileClose(output);

    return;
}


void LALInferenceWriteMCMCSamples(LALInferenceRunState *runState) {
    //ProcessParamsTable *ppt;
    INT4 MPIrank;
    INT4 t, n_local_threads;
    LALH5File *output = NULL;
    LALInferenceThreadState *thread;

    MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

    output = XLALH5FileOpen(runState->outFileName, "w");
    if(output == NULL){
        XLALErrorHandler = XLALExitErrorHandler;
        XLALPrintError("Output file error. Please check that the specified path exists. (in %s, line %d)\n",__FILE__, __LINE__);
        XLAL_ERROR_VOID(XLAL_EIO);
    }

    LALH5File *group = LALInferenceH5CreateGroupStructure(output, "lalinference", runState->runID);

    n_local_threads = runState->nthreads;
    for (t = 0; t < n_local_threads; t++) {
        thread = runState->threads[t];

        LALInferenceVariables **output_array=NULL;
        UINT4 N_output_array=0;
        if(LALInferenceCheckVariable(thread->algorithmParams, "outputarray")
                && LALInferenceCheckVariable(thread->algorithmParams, "N_outputarray") ) {
            output_array=*(LALInferenceVariables ***)LALInferenceGetVariable(thread->algorithmParams,"outputarray");
            N_output_array=*(UINT4 *)LALInferenceGetVariable(thread->algorithmParams,"N_outputarray");
        }

        /* Create run identifier group */
        LALInferenceH5VariablesArrayToDataset(group, output_array, N_output_array, thread->name);

        /* TODO: Write metadata */
    }
    XLALH5FileClose(group);

    XLALH5FileClose(output);
    return;
}


void LALInferencePrintPTMCMCHeaderFile(LALInferenceRunState *runState, LALInferenceThreadState *thread, FILE *threadoutput) {
    INT4 MPIrank, nthreads;
    INT4 randomseed, nPar, Niter;
    LALInferenceIFOData *ifodata1;
    char *arg_str;

    MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);

    nthreads = LALInferenceGetINT4Variable(runState->algorithmParams, "ntemps");

    randomseed = LALInferenceGetINT4Variable(runState->algorithmParams, "random_seed");
    nPar = LALInferenceGetVariableDimensionNonFixed(thread->currentParams);
    Niter = LALInferenceGetINT4Variable(runState->algorithmParams, "nsteps");

    REAL8 f_ref = 0.0;
    if(LALInferenceCheckVariable(thread->currentParams, "f_ref"))
        f_ref = LALInferenceGetREAL8Variable(thread->currentParams, "f_ref");

    UINT4 nIFO = 0;
    ifodata1 = runState->data;
    while(ifodata1){
        nIFO++;
        ifodata1 = ifodata1->next;
    }

    UINT4 waveform = 0;
    if(LALInferenceCheckVariable(thread->currentParams, "LAL_APPROXIMANT"))
        waveform = LALInferenceGetUINT4Variable(thread->currentParams, "LAL_APPROXIMANT");

    REAL8 pnorder = 0.0;
    if(LALInferenceCheckVariable(thread->currentParams,"LAL_PNORDER"))
        pnorder = ((REAL8)LALInferenceGetINT4Variable(thread->currentParams, "LAL_PNORDER"))/2.0;

    arg_str = LALInferencePrintCommandLine(runState->commandLine);

    ifodata1 = runState->data;
    REAL8 networkSNR=0.0;
    while(ifodata1){
        networkSNR += ifodata1->SNR * ifodata1->SNR;
        ifodata1 = ifodata1->next;
    }
    networkSNR = sqrt(networkSNR);

    REAL8 SampleRate = 4096.0; //default value of the sample rate from LALInferenceReadData()
    if(LALInferenceGetProcParamVal(runState->commandLine, "--srate"))
        SampleRate = atof(LALInferenceGetProcParamVal(runState->commandLine, "--srate")->value);

    UINT4 benchmark=0;
    if(LALInferenceGetProcParamVal(runState->commandLine, "--benchmark"))
        benchmark = 1;

    /* Print version info */
    fprintf(threadoutput, "  LALInference version:%s,%s,%s,%s,%s\n", lalAppsVCSId,lalAppsVCSDate,lalAppsVCSBranch,lalAppsVCSAuthor,lalAppsVCSStatus);
    fprintf(threadoutput,"  %s\n", arg_str);

    /* Print algorithm parameters */
    fprintf(threadoutput, "%10s  %6s  %20s  %6s %6s  %10s  %12s  %9s  %9s %8s %8s\n",
            "nIter", "seed", "null_likelihood", "Ndet", "nTemps",
            "Tchain", "NetworkSNR", "Waveform", "pNorder", "Npar", "f_ref");
    fprintf(threadoutput, "%10d  %u  %20.10lf  %6d %6d %12.1f %14.6f  %9i  %12.1f  %8i %12.1f\n",
            Niter, randomseed, thread->nullLikelihood, nIFO, nthreads,
            thread->temperature, networkSNR, waveform, pnorder, nPar, f_ref);

    /* Print detector-specific settings */
    fprintf(threadoutput, "\n%16s  %16s  %10s  %10s  %20s  %15s  %12s\n",
            "Detector", "SNR", "f_low", "f_high",
            "Sample_start", "Sample_length", "Sample_rate");
    ifodata1=runState->data;
    while(ifodata1){
        fprintf(threadoutput, "%16s  %16.8lf  %10.2lf  %10.2lf  %20.8lf  %15.7lf  %.1f\n",
                ifodata1->detector->frDetector.name, ifodata1->SNR, ifodata1->fLow, ifodata1->fHigh,
                XLALGPSGetREAL8(&(ifodata1->epoch)), atof(LALInferenceGetProcParamVal(runState->commandLine,"--seglen")->value),
                SampleRate);
        ifodata1=ifodata1->next;
    }

    /* Print column header */
    fprintf(threadoutput, "\n\n%31s\n","");
    fprintf(threadoutput, "cycle\tlogpost\tlogprior\t");
    LALInferenceFprintParameterHeaders(threadoutput, thread->currentParams);

    /* Print the likelihood and SNR of each individual detector */
    fprintf(threadoutput, "logl\t");
    fprintf(threadoutput, "deltalogl\t");
    LALInferenceIFOData *headIFO = runState->data;
    while (headIFO != NULL) {
        fprintf(threadoutput, "logl");
        fprintf(threadoutput, "%s",headIFO->name);
        fprintf(threadoutput, "\t");
        headIFO = headIFO->next;
    }
    if (LALInferenceGetProcParamVal(runState->commandLine, "--output-SNRs")) {
        headIFO = runState->data;
        while (headIFO != NULL) {
            fprintf(threadoutput, "SNR");
            fprintf(threadoutput, "%s",headIFO->name);
            fprintf(threadoutput, "\t");
            headIFO = headIFO->next;
        }
        fprintf(threadoutput, "SNR\t");
    }

    /* Print the cumulative runtime at each sample */
    if (benchmark)
        fprintf(threadoutput, "timestamp\t");
    fprintf(threadoutput,"\n");
}

void LALInferencePrintAdaptationHeader(FILE *outfile, LALInferenceThreadState *thread) {
    LALInferenceVariableItem *item;

    for (item=thread->currentParams->head; item; item=item->next) {
        if (LALInferenceCheckVariableNonFixed(thread->currentParams, item->name) &&
                item->type == LALINFERENCE_REAL8_t) {
            fprintf(outfile, "\tsigma_%s",
                    LALInferenceTranslateInternalToExternalParamName(item->name));
        }
    }

    for (item=thread->currentParams->head; item; item=item->next) {
        if (LALInferenceCheckVariableNonFixed(thread->currentParams, item->name) &&
                item->type == LALINFERENCE_REAL8_t) {
            fprintf(outfile, "\tPaccept_%s",
                    LALInferenceTranslateInternalToExternalParamName(item->name));
        }
    }

    fprintf(outfile,"\n");

    return;
}

void LALInferencePrintPTMCMCInjectionSample(LALInferenceRunState *runState) {
    ProcessParamsTable *ppt;
    LALInferenceThreadState *thread = runState->threads[0];

    ppt = LALInferenceGetProcParamVal(runState->commandLine, "--inj");
    if (ppt) {
        ProcessParamsTable *ppt2 = LALInferenceGetProcParamVal(runState->commandLine, "--outfile");
        INT4 randomseed = LALInferenceGetINT4Variable(runState->algorithmParams, "random_seed");
        FILE *out = NULL;
        char *fname = NULL;
        LALInferenceVariables *saveParams = NULL;

        saveParams = (LALInferenceVariables *)XLALCalloc(sizeof(LALInferenceVariables), 1);

        if (ppt2) {
            fname = (char *) XLALCalloc((strlen(ppt2->value)+255)*sizeof(char), 1);
            sprintf(fname, "%s.injection", ppt2->value);
        } else {
            fname = (char *) XLALCalloc(255*sizeof(char), 1);
            sprintf(fname, "PTMCMC.output.%u.injection", randomseed);
        }
        out = fopen(fname, "w");

        LALInferenceCopyVariables(thread->currentParams, saveParams);

        SimInspiralTable *injTable = NULL;
        SimInspiralTable *theEventTable = NULL;

        SimInspiralTableFromLIGOLw(&injTable,ppt->value,0,0);

        ppt2 = LALInferenceGetProcParamVal(runState->commandLine, "--event");
        if (ppt2) {
            UINT4 event = atoi(ppt2->value);
            UINT4 i;
            theEventTable = injTable;
        for (i = 0; i < event; i++) {
            theEventTable = theEventTable->next;
        }
            theEventTable->next = NULL;
        } else {
            theEventTable=injTable;
            theEventTable->next = NULL;
        }

        REAL8 m1 = theEventTable->mass1;
        REAL8 m2 = theEventTable->mass2;
        REAL8 q = m2/m1;
        REAL8 eta = m1*m2/(m1+m2)/(m1+m2);

        if (q > 1.0) q = 1.0/q;

        REAL8 sx = theEventTable->spin1x;
        REAL8 sy = theEventTable->spin1y;
        REAL8 sz = theEventTable->spin1z;

        REAL8 a_spin1 = sqrt(sx*sx + sy*sy + sz*sz);

        REAL8 theta_spin1, phi_spin1;
        if (a_spin1 == 0.0) {
            theta_spin1 = 0.0;
            phi_spin1 = 0.0;
        } else {
            theta_spin1 = acos(sz / a_spin1);
            phi_spin1 = atan2(sy, sx);
            if (phi_spin1 < 0.0) phi_spin1 += 2.0*M_PI;
        }

        sx = theEventTable->spin2x;
        sy = theEventTable->spin2y;
        sz = theEventTable->spin2z;

        REAL8 a_spin2 = sqrt(sx*sx + sy*sy + sz*sz), theta_spin2, phi_spin2;
        if (a_spin2 == 0.0) {
            theta_spin2 = 0.0;
            phi_spin2 = 0.0;
        } else {
            theta_spin2 = acos(sz / a_spin2);
            phi_spin2 = atan2(sy, sx);
            if (phi_spin2 < 0.0) phi_spin2 += 2.0*M_PI;
        }

        REAL8 psi = theEventTable->polarization;
        if (psi>=M_PI) psi -= M_PI;

        REAL8 injGPSTime = XLALGPSGetREAL8(&(theEventTable->geocent_end_time));

        REAL8 chirpmass = theEventTable->mchirp;

        REAL8 dist = theEventTable->distance;
        REAL8 inclination = theEventTable->inclination;
        REAL8 phase = theEventTable->coa_phase;
        REAL8 dec = theEventTable->latitude;
        REAL8 ra = theEventTable->longitude;

        LALInferenceSetVariable(thread->currentParams, "chirpmass", &chirpmass);
        if (LALInferenceCheckVariable(thread->currentParams, "q")) {
            LALInferenceSetVariable(thread->currentParams, "q", &q);
        } else if (LALInferenceCheckVariable(thread->currentParams, "eta")) {
            LALInferenceSetVariable(thread->currentParams, "eta", &eta);
        } else {
            /* Restore state, cleanup, and throw error */
            LALInferenceClearVariables(thread->currentParams);
            LALInferenceCopyVariables(saveParams, thread->currentParams);
            XLALFree(fname);
            LALInferenceClearVariables(saveParams);
            XLALFree(saveParams);
            XLAL_ERROR_VOID(XLAL_EINVAL, "unknown mass ratio parameter name (allowed are 'eta' or 'q')");
        }

        UINT4 added_time_param = 0;
        if (!LALInferenceCheckVariable(thread->currentParams, "time")) {
            added_time_param = 1;
            LALInferenceAddVariable(thread->currentParams, "time", &injGPSTime, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        } else {
            LALInferenceSetVariable(thread->currentParams, "time", &injGPSTime);
        }

        UINT4 added_phase_param = 0;
        if (!LALInferenceCheckVariable(thread->currentParams, "phase")) {
            added_phase_param = 1;
            LALInferenceAddVariable(thread->currentParams, "phase", &phase, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        } else {
            LALInferenceSetVariable(thread->currentParams, "phase", &phase);
        }
        REAL8 cosinclination=cos(inclination);
        LALInferenceSetVariable(thread->currentParams, "distance", &dist);
        LALInferenceSetVariable(thread->currentParams, "costheta_jn", &cosinclination);
        LALInferenceSetVariable(thread->currentParams, "polarisation", &(psi));
        LALInferenceSetVariable(thread->currentParams, "declination", &dec);
        LALInferenceSetVariable(thread->currentParams, "rightascension", &ra);
        if (LALInferenceCheckVariable(thread->currentParams, "a_spin1")) {
            LALInferenceSetVariable(thread->currentParams, "a_spin1", &a_spin1);
        }
        if (LALInferenceCheckVariable(thread->currentParams, "theta_spin1")) {
            LALInferenceSetVariable(thread->currentParams, "theta_spin1", &theta_spin1);
        }
        if (LALInferenceCheckVariable(thread->currentParams, "phi_spin1")) {
            LALInferenceSetVariable(thread->currentParams, "phi_spin1", &phi_spin1);
        }
        if (LALInferenceCheckVariable(thread->currentParams, "a_spin2")) {
            LALInferenceSetVariable(thread->currentParams, "a_spin2", &a_spin2);
        }
        if (LALInferenceCheckVariable(thread->currentParams, "theta_spin2")) {
            LALInferenceSetVariable(thread->currentParams, "theta_spin2", &theta_spin2);
        }
        if (LALInferenceCheckVariable(thread->currentParams, "phi_spin2")) {
            LALInferenceSetVariable(thread->currentParams, "phi_spin2", &phi_spin2);
        }

        thread->currentLikelihood = runState->likelihood(thread->currentParams, runState->data, thread->model);
        thread->currentPrior = runState->prior(runState, thread->currentParams, thread->model);

        LALInferencePrintPTMCMCHeaderFile(runState, thread, out);

        fclose(out);

        if (added_time_param) {
            LALInferenceRemoveVariable(thread->currentParams, "time");
            LALInferenceRemoveMinMaxPrior(runState->priorArgs, "time");
        }

        if (added_phase_param) {
            LALInferenceRemoveVariable(thread->currentParams, "phase");
            LALInferenceRemoveMinMaxPrior(runState->priorArgs, "phase");
        }

        LALInferenceCopyVariables(saveParams, thread->currentParams);
        thread->currentLikelihood = runState->likelihood(thread->currentParams, runState->data, thread->model);
        thread->currentPrior = runState->prior(runState, thread->currentParams, thread->model);

        XLALFree(fname);
        LALInferenceClearVariables(saveParams);
        XLALFree(saveParams);
    }
}

void LALInferenceSaveSample(LALInferenceThreadState *thread, FILE *output) {
    LALInferenceWriteVariablesBinary(output, thread->currentParams);
    rewind(output);

    return;
}

void LALInferencePrintMCMCSample(LALInferenceThreadState *thread, LALInferenceIFOData *data, INT4 iteration, REAL8 timestamp, FILE *threadoutput) {
    UINT4 ifo = 0;
    LALInferenceIFOData *headIFO;

    fseek(threadoutput, 0L, SEEK_END);

    fprintf(threadoutput, "%d\t%f\t%f\t",
            iteration, (thread->currentLikelihood - thread->nullLikelihood) + thread->currentPrior, thread->currentPrior);

    LALInferencePrintSample(threadoutput, thread->currentParams);
    LALInferenceLogSampleToArray(thread->algorithmParams, thread->currentParams);

    fprintf(threadoutput,"%f\t", thread->currentLikelihood);
    fprintf(threadoutput,"%f\t", thread->currentLikelihood - thread->nullLikelihood);

    ifo = 0;
    headIFO = data;
    while (headIFO != NULL) {
        fprintf(threadoutput, "%f\t", thread->currentIFOLikelihoods[ifo] - headIFO->nullloglikelihood);
        ifo++;
        headIFO = headIFO->next;
    }

    if (LALInferenceGetINT4Variable(thread->proposalArgs, "output_snrs")) {
        headIFO = data;
        ifo = 0;
        while (headIFO != NULL) {
            fprintf(threadoutput, "%f\t", thread->currentIFOSNRs[ifo]);
            headIFO = headIFO->next;
            ifo++;
        }
        fprintf(threadoutput, "%f\t", thread->model->SNR);
    }

    if (timestamp > 0)
        fprintf(threadoutput, "%f\t", timestamp);

    fprintf(threadoutput,"\n");
    fflush(threadoutput);
    return;
}

void LALInferencePrintAdaptationSettings(FILE *outfile, LALInferenceThreadState *thread) {
    LALInferenceVariableItem *item;
    REAL8 s_gamma = 0.;
    REAL8 sigma, naccepted, nproposed;

    if (LALInferenceCheckVariable(thread->proposalArgs, "s_gamma"))
        s_gamma = LALInferenceGetREAL8Variable(thread->proposalArgs, "s_gamma");
    fprintf(outfile,"%f\t",s_gamma);

    for (item=thread->currentParams->head; item; item=item->next) {
        if (LALInferenceCheckVariableNonFixed(thread->currentParams, item->name) &&
                item->type == LALINFERENCE_REAL8_t) {
            char tmpname[MAX_STRLEN]="";
            sprintf(tmpname, "%s_%s", item->name, ADAPTSUFFIX);
            sigma = LALInferenceGetREAL8Variable(thread->proposalArgs, tmpname);
            fprintf(outfile, "%g\t", sigma);
        }
    }

    for (item=thread->currentParams->head; item; item=item->next) {
        if (LALInferenceCheckVariableNonFixed(thread->currentParams, item->name) &&
                item->type == LALINFERENCE_REAL8_t) {
            char tmpname[MAX_STRLEN]="";

            sprintf(tmpname, "%s_%s", item->name, ACCEPTSUFFIX);
            naccepted = LALInferenceGetREAL8Variable(thread->proposalArgs, tmpname);

            sprintf(tmpname, "%s_%s", item->name, PROPOSEDSUFFIX);
            nproposed = LALInferenceGetREAL8Variable(thread->proposalArgs, tmpname);

            fprintf(outfile, "%f\t", naccepted/( nproposed==0 ? 1.0 : nproposed ));
        }
    }

    fprintf(outfile,"\n");
    return;
}

void LALInferenceDataDump(LALInferenceIFOData *data, LALInferenceModel *model) {
    const UINT4 nameLength=256;
    char filename[nameLength];
    FILE *out;
    UINT4 ui;

    snprintf(filename, nameLength, "freqTemplatehPlus.dat");
    out = fopen(filename, "w");
    for (ui = 0; ui < model->freqhPlus->data->length; ui++) {
        REAL8 f = model->freqhPlus->deltaF * ui;
        COMPLEX16 d = model->freqhPlus->data->data[ui];

        fprintf(out, "%g %g %g\n", f, creal(d), cimag(d));
    }
    fclose(out);

    snprintf(filename, nameLength, "freqTemplatehCross.dat");
    out = fopen(filename, "w");
    for (ui = 0; ui < model->freqhCross->data->length; ui++) {
        REAL8 f = model->freqhCross->deltaF * ui;
        COMPLEX16 d = model->freqhCross->data->data[ui];

        fprintf(out, "%g %g %g\n", f, creal(d), cimag(d));
    }
    fclose(out);

    while (data != NULL) {
        snprintf(filename, nameLength, "%s-freqTemplateStrain.dat", data->name);
        out = fopen(filename, "w");
        for (ui = 0; ui < model->freqhCross->data->length; ui++) {
            REAL8 f = model->freqhCross->deltaF * ui;
            COMPLEX16 d;
            d = data->fPlus * model->freqhPlus->data->data[ui] +
            data->fCross * model->freqhCross->data->data[ui];

            fprintf(out, "%g %g %g\n", f, creal(d), cimag(d) );
        }
        fclose(out);

        snprintf(filename, nameLength, "%s-timeTemplateStrain.dat", data->name);
        out = fopen(filename, "w");
        for (ui = 0; ui < model->timehCross->data->length; ui++) {
            REAL8 tt = XLALGPSGetREAL8(&(model->timehCross->epoch)) +
            data->timeshift + ui*model->timehCross->deltaT;
            REAL8 d = data->fPlus*model->timehPlus->data->data[ui] +
            data->fCross*model->timehCross->data->data[ui];

            fprintf(out, "%.6f %g\n", tt, d);
        }
        fclose(out);

        snprintf(filename, nameLength, "%s-timeTemplatehPlus.dat", data->name);
        out = fopen(filename, "w");
        for (ui = 0; ui < model->timehPlus->data->length; ui++) {
            REAL8 tt = XLALGPSGetREAL8(&(model->timehPlus->epoch)) +
            ui * model->timehPlus->deltaT;
            REAL8 d = model->timehPlus->data->data[ui];

            fprintf(out, "%.6f %g\n", tt, d);
        }
        fclose(out);

        snprintf(filename, nameLength, "%s-timeTemplatehCross.dat", data->name);
        out = fopen(filename, "w");
        for (ui = 0; ui < model->timehCross->data->length; ui++) {
            REAL8 tt = XLALGPSGetREAL8(&(model->timehCross->epoch)) +
            ui * model->timehCross->deltaT;
            REAL8 d = model->timehCross->data->data[ui];

            fprintf(out, "%.6f %g\n", tt, d);
        }
        fclose(out);

        data = data->next;
    }
}

void LALInferenceMCMCResumeRead(LALInferenceThreadState *thread, FILE *resumeFile) {
    LALInferenceClearVariables(thread->currentParams);
    thread->currentParams = LALInferenceReadVariablesBinary(resumeFile);

    return;
}
