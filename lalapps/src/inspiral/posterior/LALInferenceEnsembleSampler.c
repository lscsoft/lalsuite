/*
 *  LALInferenceEnsembleSampler.c:  Bayesian Followup, ensemble-sampling algorithm.
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

#include <mpi.h>
#include <lal/LALInference.h>
#include "LALInferenceEnsembleSampler.h"
#include <lal/LALInferenceProposal.h>

#include <LALAppsVCSInfo.h>

#define PROGRAM_NAME "LALInferenceEnsembleSampler.c"
#define CVS_ID_STRING "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define CVS_NAME_STRING "$Name$"

void ensemble_sampler(struct tagLALInferenceRunState *runState) {
    INT4 MPIrank, MPIsize;
    INT4 walker, nwalkers_per_thread, nsteps;
    INT4 skip, update_interval, verbose;
    INT4 i,t,c;
    char **walker_output_names = NULL;

    /* Initialize LIGO status */
    LALStatus status;
    memset(&status, 0, sizeof(status));

    MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);    // This thread's index
    MPI_Comm_size(MPI_COMM_WORLD, &MPIsize);    // Number of MPI processes

    /* Parameters controlling output and ensemble update frequency */
    LALInferenceVariables *algorithmParams = runState->algorithmParams;
    nwalkers_per_thread = *(INT4*) LALInferenceGetVariable(algorithmParams, "nwalkers_per_thread");
    nsteps = *(INT4*) LALInferenceGetVariable(algorithmParams, "nsteps");
    skip = *(INT4*) LALInferenceGetVariable(algorithmParams, "skip");
    update_interval = *(INT4*) LALInferenceGetVariable(algorithmParams, "update_interval");
    verbose = *(INT4*) LALInferenceGetVariable(algorithmParams, "verbose");

    /* Initialize walker outputs, closing due to the large numbers that are likley */
    walker_output_names = XLALMalloc(nwalkers_per_thread * sizeof(char*));
    for (walker=0; walker<nwalkers_per_thread; walker++)
        walker_output_names[walker] =
            LALInferenceInitializeEnsembleOutput(runState, walker,
                                                    MPIrank*nwalkers_per_thread, verbose);

    /* Setup clustered-KDE proposal from the inital state of the ensemble */
    ensemble_update(runState);

    UINT4 step = 0;
    while (step < nsteps) {
        step++;

        /* Update the proposal from the current state of the ensemble */
        if ((step % update_interval) == 0)
            ensemble_update(runState);

        /* Update all walkers on this MPI-thread */
        #pragma omp parallel for
        for (walker=0; walker<nwalkers_per_thread; walker++) {
            INT4 accepted;
            LALInferenceVariables proposedParams;
            proposedParams.head = NULL;
            proposedParams.dimension = 0;

            accepted = walker_step(runState, runState->modelArray[walker],
                                    runState->currentParamArray[walker],
                                    &proposedParams,
                                    &(runState->currentPriors[walker]),
                                    &(runState->currentLikelihoods[walker]));

            /* Don't print every sample to file */
            if ((step % skip) == 0) {
                LALInferencePrintEnsembleSample(runState, walker_output_names, walker, step);
                if (verbose)
                    LALInferencePrintProposedSample(runState, &proposedParams, walker, accepted);
            }
            LALInferenceClearVariables(&proposedParams);
        }
    }

    for (walker=0; walker<nwalkers_per_thread; walker++)
        XLALFree(walker_output_names[walker]);
    XLALFree(walker_output_names);
}

INT4 walker_step(LALInferenceRunState *runState,
                    LALInferenceModel *model,
                    LALInferenceVariables *currentParams,
                    LALInferenceVariables *proposedParams,
                    REAL8 *currentPrior, REAL8 *currentLikelihood) {
    REAL8 log_prior_proposed, log_likelihood_proposed;
    REAL8 log_proposal_ratio = 0.0;
    REAL8 log_acceptance_probability;
    INT4 accepted = 0;
    LALInferenceIFOData *headData;

    /* Propose a new sample */
    LALInferenceClearVariables(proposedParams);

    /* Get the probability of proposing the reverse jump */
    log_proposal_ratio = runState->proposal(runState, currentParams, proposedParams);

    /* Only bother calculating likelihood if within prior boundaries */
    log_prior_proposed = runState->prior(runState, proposedParams, model);
    if (log_prior_proposed > -DBL_MAX)
        log_likelihood_proposed =
            runState->likelihood(proposedParams, runState->data, model);
    else
        log_likelihood_proposed = -INFINITY;

    /* Find jump acceptance probability */
    log_acceptance_probability = (log_prior_proposed + log_likelihood_proposed)
                                - (*currentPrior +
                                    *currentLikelihood)
                                + log_proposal_ratio;

    /* Accept the jump with the calculated probability */
    if (log_acceptance_probability > 0
            || (log(gsl_rng_uniform(runState->GSLrandom)) < log_acceptance_probability)) {
        LALInferenceCopyVariables(proposedParams, currentParams);
        *currentPrior = log_prior_proposed;
        *currentLikelihood = log_likelihood_proposed;

        headData = runState->data;
        while (headData != NULL) {
            headData->acceptedloglikelihood = headData->loglikelihood;
            headData->acceptedSNR = headData->currentSNR;
            headData = headData->next;
        }
        accepted = 1;
    }

    return accepted;
}


void ensemble_update(LALInferenceRunState *runState) {
    INT4 nwalkers, nwalkers_per_thread, walker, ndim;
    REAL8 *parameters, *samples, *param_array;

    LALInferenceVariables *algorithmParams = runState->algorithmParams;
    nwalkers = *(INT4 *) LALInferenceGetVariable(algorithmParams, "nwalkers");
    nwalkers_per_thread = *(INT4*) LALInferenceGetVariable(algorithmParams, "nwalkers_per_thread");

    /* Get current dimension */
    ndim = LALInferenceGetVariableDimensionNonFixed(runState->currentParamArray[0]);

    /* Prepare array to constain samples */
    samples = (REAL8*) XLALMalloc(nwalkers * ndim * sizeof(REAL8));

    /* Get this walkers current location */
    param_array = XLALMalloc(nwalkers_per_thread * ndim * sizeof(REAL8));
    for (walker=0; walker<nwalkers_per_thread; walker++) {
        parameters = &(param_array[ndim*walker]);
        LALInferenceCopyVariablesToArray(runState->currentParamArray[walker], parameters);
    }

    /* Gather samples from the ensemble */
    MPI_Allgather(param_array, nwalkers_per_thread*ndim,
                    MPI_DOUBLE, samples, nwalkers_per_thread*ndim,
                    MPI_DOUBLE, MPI_COMM_WORLD);

    /* Update the KDE proposal */
    UINT4 ntrials = 50;  // Number of trials to optimize clustering at fixed-k
    LALInferenceSetupClusteredKDEProposalFromRun(runState, samples, nwalkers, ntrials);

    /* Cleanup */
    XLALFree(param_array);
    XLALFree(samples);
}


//-----------------------------------------
// file output routines:
//-----------------------------------------
char *LALInferenceInitializeEnsembleOutput(LALInferenceRunState *runState,
                                            INT4 walker,
                                            INT4 walker_offset,
                                            INT4 verbose) {
    ProcessParamsTable *ppt;
    UINT4 randomseed;
    char *outfile_name = NULL;
    char *prop_name = NULL;
    FILE *walker_output = NULL;
    FILE *prop_output = NULL;

    /* Randomseed used to prevent overwriting when peforming multiple analyses */
    randomseed = *(UINT4*) LALInferenceGetVariable(runState->algorithmParams,"random_seed");

    ppt = LALInferenceGetProcParamVal(runState->commandLine, "--outfile");
    if (ppt) {
        outfile_name = (char*) XLALCalloc(strlen(ppt->value)+255, sizeof(char*));
        sprintf(outfile_name, "%s.%2.2d", ppt->value, walker);
    } else {
        outfile_name = (char*) XLALCalloc(255, sizeof(char*));
        sprintf(outfile_name, "ensemble.output.%u.%2.2d", randomseed, walker_offset+walker);
        if (verbose) {
            prop_name = (char*) XLALCalloc(255, sizeof(char*));
            sprintf(prop_name, "ensemble.proposed.%u.%2.2d", randomseed, walker_offset+walker);
        }
    }

    walker_output = fopen(outfile_name, "w");
    if (walker_output == NULL){
        XLALErrorHandler = XLALExitErrorHandler;
        XLALPrintError("Output file error in %s, line %d. %s.\n",
                        __FILE__, __LINE__, strerror(errno));
        XLAL_ERROR_NULL(XLAL_EIO);
    }

    LALInferencePrintEnsembleHeader(runState, walker_output, walker);
    fclose(walker_output);
    if (verbose) {
        FILE *prop_out = fopen(prop_name, "w");
        LALInferenceFprintParameterNonFixedHeaders(prop_out, runState->currentParamArray[0]);
        fprintf(prop_out, "accepted\n");
        fclose(prop_out);
        XLALFree(prop_name);
    }

    walker_output = fopen(outfile_name, "a");
    if (walker_output == NULL) {
        XLALErrorHandler = XLALExitErrorHandler;
        XLALPrintError("Output file error. Please check that the specified path exists. \
                        (in %s, line %d)\n",__FILE__, __LINE__);
        XLAL_ERROR_NULL(XLAL_EIO);
    }

    fclose(walker_output);
    return outfile_name;
}

void LALInferencePrintEnsembleSample(LALInferenceRunState *runState,
                                        char **walker_output_names,
                                        UINT4 walker,
                                        INT4 step) {
    REAL8 null_likelihood, timestamp, timestamp_epoch;
    REAL8 networkSNR, normed_logl;
    REAL8 *currentPriors, *currentLikelihoods;
    UINT4 benchmark;
    FILE *walker_output;

    walker_output = fopen(walker_output_names[walker], "a");
    currentPriors = runState->currentPriors;
    currentLikelihoods = runState->currentLikelihoods;

    /* Print step number, log(posterior), and log(prior) */
    null_likelihood = *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "nullLikelihood");
    fprintf(walker_output, "%d\t", step);
    fprintf(walker_output, "%f\t",
            (currentLikelihoods[walker] - null_likelihood) + currentPriors[walker]);
    fprintf(walker_output, "%f\t", currentPriors[walker]);

    /* Print the non-fixed parameter values */
    LALInferencePrintSampleNonFixed(walker_output, runState->currentParamArray[walker]);

    /* Print network and single-IFO likelihoods */
    fprintf(walker_output, "%f\t", currentLikelihoods[walker] - null_likelihood);

    LALInferenceIFOData *ifo_data = runState->data;
    while (ifo_data != NULL) {
        normed_logl = ifo_data->acceptedloglikelihood-ifo_data->nullloglikelihood;
        fprintf(walker_output, "%f\t", normed_logl);
        ifo_data = ifo_data->next;
    }

    /* Print network and single-IFO SNRs */
    networkSNR = 0.0;
    ifo_data = runState->data;
    while (ifo_data != NULL) {
        fprintf(walker_output, "%f\t", ifo_data->acceptedSNR);
        networkSNR += ifo_data->acceptedSNR * ifo_data->acceptedSNR;
        ifo_data = ifo_data->next;
    }
    networkSNR = sqrt(networkSNR);
    fprintf(walker_output, "%f\t", networkSNR);

    benchmark = *(UINT4 *) LALInferenceGetVariable(runState->algorithmParams, "benchmark");
    if (benchmark) {
        struct timeval tv;
        timestamp_epoch = *(REAL8 *) LALInferenceGetVariable(runState->algorithmParams,
                                                                   "timestamp_epoch");

        gettimeofday(&tv, NULL);
        timestamp = tv.tv_sec + tv.tv_usec/1E6 - timestamp_epoch;
        fprintf(walker_output, "%f\t", timestamp);
    }
    fprintf(walker_output,"\n");

    fclose(walker_output);
}


void LALInferencePrintProposedSample(LALInferenceRunState *runState,
                                        LALInferenceVariables *proposedParams,
                                        INT4 walker,
                                        INT4 accepted) {
    INT4 MPIrank;
    INT4 nwalkers_per_thread;
    UINT4 randomseed;
    FILE *output = NULL;
    char *outname = (char*) XLALCalloc(255, sizeof(char*));

    MPI_Comm_rank(MPI_COMM_WORLD, &MPIrank);    // This thread's index

    LALInferenceVariables *algorithmParams = runState->algorithmParams;
    nwalkers_per_thread = *(INT4*) LALInferenceGetVariable(algorithmParams, "nwalkers_per_thread");
    randomseed = *(UINT4*) LALInferenceGetVariable(algorithmParams,"random_seed");
    sprintf(outname, "ensemble.proposed.%u.%2.2d", randomseed, nwalkers_per_thread*MPIrank+walker);

    output = fopen(outname, "a");
    LALInferencePrintSampleNonFixed(output, proposedParams);
    fprintf(output, "%i\n", accepted);
    fclose(output);

    XLALFree(outname);
}

void LALInferencePrintEnsembleHeader(LALInferenceRunState *runState,
                                        FILE *walker_output,
                                        INT4 walker) {
    LALInferenceVariables **currentParamArray;
    LALInferenceIFOData *ifo_data;
    REAL8TimeSeries *time_data;
    INT4 ndim, int_pn_order, waveform = 0;
    UINT4 nifo, randomseed, benchmark;
    REAL8 null_likelihood, normed_logl, pn_order;
    REAL8 network_snr, sampling_rate;
    REAL8 delta_t, timestamp, f_ref = 0.0;
    struct timeval tv;
    char *cmd_str;

    /* Save the command line for repoducability */
    cmd_str = LALInferencePrintCommandLine(runState->commandLine);

    randomseed = *(UINT4*) LALInferenceGetVariable(runState->algorithmParams, "random_seed");
    null_likelihood = *(REAL8*) LALInferenceGetVariable(runState->proposalArgs, "nullLikelihood");
    currentParamArray = runState->currentParamArray;
    ndim = LALInferenceGetVariableDimensionNonFixed(currentParamArray[walker]);

    if (LALInferenceCheckVariable(currentParamArray[walker], "fRef"))
        f_ref = *(REAL8*) LALInferenceGetVariable(currentParamArray[walker], "fRef");

    /* Count number of detectors */
    ifo_data = runState->data;
    nifo = 0;
    while (ifo_data){
      nifo++;
      ifo_data = ifo_data->next;
    }

    /* Integer (from an enum) identfying the waveform family used */
    if (LALInferenceCheckVariable(currentParamArray[walker], "LAL_APPROXIMANT"))
        waveform = *(INT4 *) LALInferenceGetVariable(currentParamArray[walker], "LAL_APPROXIMANT");

    /* Determine post-Newtonian (pN) order (half of the integer stored in currentParams) */
    if (LALInferenceCheckVariable(currentParamArray[walker], "LAL_PNORDER"))
        int_pn_order = *(INT4*)LALInferenceGetVariable(currentParamArray[walker], "LAL_PNORDER");
    pn_order = int_pn_order/2.0;

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
            lalAppsVCSId, lalAppsVCSDate, lalAppsVCSBranch,
            lalAppsVCSAuthor, lalAppsVCSStatus);

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

    LALInferenceFprintParameterNonFixedHeaders(walker_output, currentParamArray[walker]);

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
    normed_logl = runState->currentLikelihoods[walker]-null_likelihood;
    fprintf(walker_output, "%d\t%f\t%f\t", 0,
            runState->currentPriors[walker] + normed_logl,
            runState->currentPriors[walker]);

    LALInferencePrintSampleNonFixed(walker_output, currentParamArray[walker]);

    fprintf(walker_output, "%f\t", runState->currentLikelihoods[walker] - null_likelihood);
    ifo_data = runState->data;
    while (ifo_data != NULL) {
        normed_logl = ifo_data->acceptedloglikelihood-ifo_data->nullloglikelihood;
        fprintf(walker_output, "%f\t", normed_logl);
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
