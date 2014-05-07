/*
 *  LALInferenceMCMC.c:  Bayesian Followup, MCMC algorithm.
 *
 *  Copyright (C) 2009, 2012 Ilya Mandel, Vivien Raymond, Christian
 *  Roever, Marc van der Sluys, John Veitch and Will M. Farr
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
#include <stdlib.h>
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
#include "LALInferenceEnsembleSampler.h"
#include "LALInferenceMCMCSampler.h"
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferenceProposal.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <sys/time.h>

#include <LALAppsVCSInfo.h>
#include <lal/LALStdlib.h>

#define PROGRAM_NAME "LALInferenceEnsembleSampler.c"
#define CVS_ID_STRING "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define CVS_NAME_STRING "$Name$"

void ensemble_sampler(struct tagLALInferenceRunState *runState) {
    INT4 walker, nwalkers;
    INT4 i,t,c;
    FILE *walker_output = NULL;

    /* Initialize LIGO status */
    LALStatus status;
    memset(&status,0,sizeof(status));

    MPI_Comm_rank(MPI_COMM_WORLD, &walker);      // This walker's index
    MPI_Comm_size(MPI_COMM_WORLD, &nwalkers);    // Size of the ensemble

    /* Parameters controlling output and ensemble update frequency */
    INT4 nsteps = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "nsteps");
    INT4 skip = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "skip");
    INT4 update_interval = *(INT4*) LALInferenceGetVariable(runState->algorithmParams, "update_interval");

    walker_output = LALInferenceInitializeEnsembleOutput(runState);

    /* Setup clustered-KDE proposal from the inital state of the ensemble */
    ensemble_update(runState);

    UINT4 step = 0;
    while (step < nsteps) {
        step++;

        /* Update the proposal from the current state of the ensemble */
        if ((step % update_interval) == 0)
            ensemble_update(runState);

        walker_step(runState); //evolve the walker

        /* Don't print every sample to file */
        if ((step % skip) == 0)
            LALInferencePrintEnsembleSample(runState, walker_output, step);
    }

    fclose(walker_output);
}

void walker_step(LALInferenceRunState *runState) {
    REAL8 log_prior_proposed, log_likelihood_proposed;
    REAL8 log_proposal_ratio = 0.0;
    REAL8 log_acceptance_probability;
    LALInferenceIFOData *headData;

    /* Propose a new sample */
    LALInferenceVariables proposedParams;
    proposedParams.head = NULL;
    proposedParams.dimension = 0;
    runState->proposal(runState, &proposedParams);

    /* Get the probability of proposing the reverse jump */
    log_proposal_ratio =
        *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "logProposalRatio");

    /* Only bother calculating likelihood if within prior boundaries */
    log_prior_proposed = runState->prior(runState, &proposedParams);
    if (log_prior_proposed > -DBL_MAX)
        log_likelihood_proposed =
            runState->likelihood(&proposedParams, runState->data, runState->templt);
    else
        log_likelihood_proposed = -INFINITY;

    /* Find jump acceptance probability */
    log_acceptance_probability = (log_prior_proposed + log_likelihood_proposed)
                                - (runState->currentPrior + runState->currentLikelihood)
                                + log_proposal_ratio;

    /* Accept the jump with the calculated probability */
    if (log_acceptance_probability > 0
            || (log(gsl_rng_uniform(runState->GSLrandom)) < log_acceptance_probability)) {
        LALInferenceCopyVariables(&proposedParams, runState->currentParams);
        runState->currentLikelihood = log_likelihood_proposed;
        runState->currentPrior = log_prior_proposed;

        headData = runState->data;
        while (headData != NULL) {
            headData->acceptedloglikelihood = headData->loglikelihood;
            headData->acceptedSNR = headData->currentSNR;
            headData = headData->next;
        }
    }

    LALInferenceClearVariables(&proposedParams);
}


void ensemble_update(LALInferenceRunState *runState) {
    INT4 walker, nwalkers;
    MPI_Comm_rank(MPI_COMM_WORLD, &walker);
    MPI_Comm_size(MPI_COMM_WORLD, &nwalkers);

    /* Get current dimension */
    INT4 ndim = LALInferenceGetVariableDimensionNonFixed(runState->currentParams);

    /* Prepare array to constain samples */
    REAL8 *samples = (REAL8*) XLALMalloc(nwalkers * ndim * sizeof(REAL8));

    /* Get this walkers current location */
    REAL8Vector *parameters = LALInferenceCopyVariablesToArray(runState->currentParams);

    /* Gather samples from the ensemble */
    MPI_Allgather(parameters->data, ndim, MPI_DOUBLE, samples, ndim, MPI_DOUBLE, MPI_COMM_WORLD);

    /* Update the KDE proposal */
    UINT4 ntrials = 50;  // Number of trials to optimize clustering at fixed-k
    LALInferenceSetupClusteredKDEProposalFromRun(runState, samples, nwalkers, ntrials);

    /* Cleanup */
    XLALDestroyREAL8Vector(parameters);
    XLALFree(samples);
}


//-----------------------------------------
// file output routines:
//-----------------------------------------
FILE *LALInferenceInitializeEnsembleOutput(LALInferenceRunState *runState) {
    ProcessParamsTable *ppt;
    char *outfile_name = NULL;
    FILE *walker_output = NULL;

    INT4 walker;
    MPI_Comm_rank(MPI_COMM_WORLD, &walker);

    /* Randomseed used to prevent overwriting when peforming multiple analyses */
    UINT4 randomseed = *(UINT4*) LALInferenceGetVariable(runState->algorithmParams,"random_seed");

    ppt = LALInferenceGetProcParamVal(runState->commandLine, "--outfile");
    if (ppt) {
        outfile_name = (char*) XLALCalloc(strlen(ppt->value)+255, sizeof(char*));
        sprintf(outfile_name, "%s.%2.2d", ppt->value, walker);
    } else {
        outfile_name = (char*) XLALCalloc(255, sizeof(char*));
        sprintf(outfile_name, "ensemble.output.%u.%2.2d", randomseed, walker);
    }

    walker_output = fopen(outfile_name, "w");
    if (walker_output == NULL){
        XLALErrorHandler = XLALExitErrorHandler;
        XLALPrintError("Output file error. Please check that the specified path exists. \
                        (in %s, line %d)\n",__FILE__, __LINE__);
        XLAL_ERROR_NULL(XLAL_EIO);
    }

    LALInferencePrintEnsembleHeader(runState, walker_output);
    fclose(walker_output);

    walker_output = fopen(outfile_name, "a");
    if (walker_output == NULL) {
        XLALErrorHandler = XLALExitErrorHandler;
        XLALPrintError("Output file error. Please check that the specified path exists. \
                        (in %s, line %d)\n",__FILE__, __LINE__);
        XLAL_ERROR_NULL(XLAL_EIO);
    }

    XLALFree(outfile_name);
    return walker_output;
}

void LALInferencePrintEnsembleSample(LALInferenceRunState *runState, FILE *walker_output, INT4 step) {
    fseek(walker_output, 0L, SEEK_END);

    /* Print step number, log(posterior), and log(prior) */
    REAL8 null_likelihood = *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "nullLikelihood");
    fprintf(walker_output, "%d\t", step);
    fprintf(walker_output, "%f\t", (runState->currentLikelihood-null_likelihood)+runState->currentPrior);
    fprintf(walker_output, "%f\t", runState->currentPrior);

    /* Print the non-fixed parameter values */
    LALInferencePrintSampleNonFixed(walker_output, runState->currentParams);

    /* Print network and single-IFO likelihoods */
    fprintf(walker_output, "%f\t", runState->currentLikelihood - null_likelihood);

    LALInferenceIFOData *ifo_data = runState->data;
    while (ifo_data != NULL) {
        fprintf(walker_output, "%f\t", ifo_data->acceptedloglikelihood - ifo_data->nullloglikelihood);
        ifo_data = ifo_data->next;
    }

    /* Print network and single-IFO SNRs */
    REAL8 networkSNR = 0.0;
    ifo_data = runState->data;
    while (ifo_data != NULL) {
        fprintf(walker_output, "%f\t", ifo_data->acceptedSNR);
        networkSNR += ifo_data->acceptedSNR * ifo_data->acceptedSNR;
        ifo_data = ifo_data->next;
    }
    networkSNR = sqrt(networkSNR);
    fprintf(walker_output, "%f\t", networkSNR);

    UINT4 benchmark = *(UINT4 *) LALInferenceGetVariable(runState->algorithmParams, "benchmark");
    if (benchmark) {
        struct timeval tv;
        REAL8 timestamp_epoch = *(REAL8 *) LALInferenceGetVariable(runState->algorithmParams,
                                                                   "timestamp_epoch");

        gettimeofday(&tv, NULL);
        REAL8 timestamp = tv.tv_sec + tv.tv_usec/1E6 - timestamp_epoch;
        fprintf(walker_output, "%f\t", timestamp);
    }
    fprintf(walker_output,"\n");
    fflush(walker_output);
}


void LALInferencePrintEnsembleHeader(LALInferenceRunState *runState, FILE *walker_output) {
    LALInferenceIFOData *ifo_data;
    REAL8TimeSeries *time_data;
    INT4 walker, ndim;
    INT4 waveform = 0;
    UINT4 nifo, randomseed, benchmark;
    REAL8 timestamp;
    REAL8 null_likelihood, pn_order;
    REAL8 network_snr, sampling_rate;
    REAL8 delta_t, f_ref = 0.0;
    struct timeval tv;
    char *cmd_str;

    /* This walker's index */
    MPI_Comm_rank(MPI_COMM_WORLD, &walker);

    /* Save the command line for repoducability */
    cmd_str = LALInferencePrintCommandLine(runState->commandLine);

    null_likelihood = *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "nullLikelihood");
    ndim = LALInferenceGetVariableDimensionNonFixed(runState->currentParams);
    randomseed = *(UINT4*) LALInferenceGetVariable(runState->algorithmParams, "random_seed");

    if (LALInferenceCheckVariable(runState->currentParams, "fRef"))
        f_ref = *(REAL8*) LALInferenceGetVariable(runState->currentParams, "fRef");

    /* Count number of detectors */
    ifo_data = runState->data;
    nifo = 0;
    while (ifo_data){
      nifo++;
      ifo_data = ifo_data->next;
    }

    /* Integer (from an enum) identfying the waveform family used */
    if (LALInferenceCheckVariable(runState->currentParams, "LAL_APPROXIMANT"))
        waveform = *(INT4 *) LALInferenceGetVariable(runState->currentParams, "LAL_APPROXIMANT");

    /* Determine post-Newtonian (pN) order (half of the integer stored in currentParams) */
    pn_order = 0.0;
    if (LALInferenceCheckVariable(runState->currentParams, "LAL_PNORDER"))
        pn_order = (REAL8)(*(INT4*) LALInferenceGetVariable(runState->currentParams, "LAL_PNORDER"));
    pn_order /= 2.0;

    /* Calculated the network signal-to-noise ratio if an injection was done */
    ifo_data = runState->data;
    network_snr = 0.0;
    while (ifo_data) {
        network_snr += ifo_data->SNR * ifo_data->SNR;
        ifo_data = ifo_data->next;
    }
    network_snr = sqrt(network_snr);

    /* Keep track of time if benchmarking */
    benchmark = *(UINT4 *) LALInferenceGetVariable(runState->algorithmParams, "benchmark");

    /* Write the header information to file */
    fprintf(walker_output,
            "  LALInference version:%s,%s,%s,%s,%s\n",
            LALAPPS_VCS_ID, LALAPPS_VCS_DATE, LALAPPS_VCS_BRANCH,
            LALAPPS_VCS_AUTHOR, LALAPPS_VCS_STATUS);

    fprintf(walker_output, "  %s\n", cmd_str);

    fprintf(walker_output, "%6s\t%20s\t%6s\t%12s\t%9s\t%9s\t%8s\t%8s\n",
        "seed", "null_likelihood", "ndet", "network_snr", "waveform", "pn_order", "ndim", "f_ref");

    fprintf(walker_output, "%u\t%20.10lf\t%6d\t%14.6f\t%9i\t%9.1f\t%8i\t%12.1f\n",
        randomseed, null_likelihood, nifo, network_snr, waveform, pn_order, ndim, f_ref);

    /* Use time step in time-domain data to determine sampling rate */
    fprintf(walker_output, "\n%16s\t%16s\t%10s\t%10s\t%10s\t%10s\t%20s\n",
        "Detector", "SNR", "f_low", "f_high", "start_time", "segment_length", "sampling_rate");
    ifo_data=runState->data;
    while(ifo_data){
        time_data = ifo_data->timeData;
        delta_t = time_data->deltaT;
        fprintf(walker_output,
                "%16s\t%16.8lf\t%10.2lf\t%10.2lf\t%15.7lf\t%12d\n",
                ifo_data->detector->frDetector.name,
                ifo_data->SNR, ifo_data->fLow, ifo_data->fHigh,
                XLALGPSGetREAL8(&(ifo_data->epoch)), time_data->data->length);
        ifo_data=ifo_data->next;
    }
    fprintf(walker_output, "\n\n\n");

    /* These are the actual column headers for the samples to be output */
    fprintf(walker_output, "cycle\tlogpost\tlogprior\t");

    LALInferenceFprintParameterNonFixedHeaders(walker_output, runState->currentParams);

    fprintf(walker_output, "logl\t");
    ifo_data = runState->data;
    while (ifo_data != NULL) {
        fprintf(walker_output, "logl");
        fprintf(walker_output, "%s", ifo_data->name);
        fprintf(walker_output, "\t");
        ifo_data = ifo_data->next;
    }
    ifo_data = runState->data;
    while (ifo_data != NULL) {
        fprintf(walker_output, "SNR");
        fprintf(walker_output, "%s", ifo_data->name);
        fprintf(walker_output, "\t");
        ifo_data = ifo_data->next;
    }
    fprintf(walker_output, "SNR\t");

    if (benchmark)
        fprintf(walker_output, "timestamp\t");
    fprintf(walker_output,"\n");

    /* Print starting values as 0th iteration */
    fprintf(walker_output, "%d\t%f\t%f\t",
            0, (runState->currentLikelihood - null_likelihood)+runState->currentPrior,
            runState->currentPrior);

    LALInferencePrintSampleNonFixed(walker_output, runState->currentParams);

    fprintf(walker_output, "%f\t", runState->currentLikelihood - null_likelihood);
    ifo_data = runState->data;
    while (ifo_data != NULL) {
        fprintf(walker_output, "%f\t", ifo_data->acceptedloglikelihood - ifo_data->nullloglikelihood);
        ifo_data = ifo_data->next;
    }

    if (benchmark) {
        gettimeofday(&tv, NULL);
        timestamp = tv.tv_sec + tv.tv_usec/1E6;
        LALInferenceAddVariable(runState->algorithmParams, "timestamp_epoch", &timestamp,
                                LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
        fprintf(walker_output, "%f\t", 0.0);
    }
    fprintf(walker_output,"\n");

    XLALFree(cmd_str);
}
