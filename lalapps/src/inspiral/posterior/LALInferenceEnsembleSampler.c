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
#include <lal/LALInferenceClusteredKDE.h>

#include <LALAppsVCSInfo.h>

#define PROGRAM_NAME "LALInferenceEnsembleSampler.c"
#define CVS_ID_STRING "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define CVS_NAME_STRING "$Name$"

void ensemble_sampler(struct tagLALInferenceRunState *run_state) {
    INT4 mpi_rank, mpi_size;
    INT4 walker, nwalkers_per_thread;
    INT4 nsteps, skip, update_interval, verbose;
    INT4 *step, i;
    INT4 *acceptance_buffer;
    REAL8 *acceptance_rates;
    REAL8 *prop_priors, *prop_likelihoods, *prop_densities;
    FILE *output = NULL;

    /* Initialize LIGO status */
    LALStatus status;
    memset(&status, 0, sizeof(status));

    /* Determine number of MPI threads, and this thread's rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    /* Parameters controlling output and ensemble update frequency */
    LALInferenceVariables *algorithm_params = run_state->algorithmParams;
    nsteps = LALInferenceGetINT4Variable(algorithm_params, "nsteps");
    skip = LALInferenceGetINT4Variable(algorithm_params, "skip");
    verbose = LALInferenceGetINT4Variable(algorithm_params, "verbose");

    step = (INT4*) LALInferenceGetVariable(algorithm_params, "step");

    nwalkers_per_thread =
        LALInferenceGetINT4Variable(algorithm_params, "nwalkers_per_thread");
    update_interval =
        LALInferenceGetINT4Variable(algorithm_params, "update_interval");

    /* Initialize walker acceptance rate tracking */
    acceptance_buffer = XLALCalloc(nwalkers_per_thread * update_interval, sizeof(INT4));
    acceptance_rates = XLALCalloc(nwalkers_per_thread, sizeof(REAL8));

    /* Create arrays for storing proposal stats that are used to calculate evidence */
    prop_priors = XLALCalloc(nwalkers_per_thread, sizeof(REAL8));
    prop_likelihoods = XLALCalloc(nwalkers_per_thread, sizeof(REAL8));
    prop_densities = XLALCalloc(nwalkers_per_thread, sizeof(REAL8));

    /* Open output and print header */
    output = init_ensemble_output(run_state, verbose, mpi_rank);

    /* Setup clustered-KDE proposal from the current state of the ensemble */
    ensemble_update(run_state);

    /* Main sampling loop */
    while (*step < nsteps) {

        /* Update step counters */
        (*step)++;

        /* Update the proposal from the current state of the ensemble */
        if ((*step % update_interval) == 0)
            ensemble_update(run_state);

        /* Update all walkers on this MPI-thread */
        #pragma omp parallel for
        for (walker=0; walker<nwalkers_per_thread; walker++) {
            INT4 accepted;

            LALInferenceVariables *proposed_params =
                XLALCalloc(1, sizeof(LALInferenceVariables));

            accepted = walker_step(run_state, run_state->modelArray[walker],
                                    run_state->currentParamArray[walker], proposed_params,
                                    &(run_state->currentPriors[walker]),
                                    &(run_state->currentLikelihoods[walker]),
                                    &(run_state->currentPropDensityArray[walker]),
                                    &(prop_priors[walker]),
                                    &(prop_likelihoods[walker]),
                                    &(prop_densities[walker]));

            /* Track acceptance rates */
            acceptance_buffer[walker + (*step % update_interval)] = accepted;
            acceptance_rates[walker] = 0.0;
            for (i = 0; i < update_interval; i++)
                acceptance_rates[walker] += acceptance_buffer[walker+i];
            acceptance_rates[walker] /= update_interval;

            if (verbose)
                print_proposed_sample(run_state, proposed_params, walker, accepted);

            LALInferenceClearVariables(proposed_params);
            XLALFree(proposed_params);
        }

        /* Output samples to file */
        if ((*step % skip) == 0)
            print_samples(run_state, output, prop_priors,
                            prop_likelihoods, prop_densities,
                            acceptance_rates, mpi_rank);
    }

    /* Sampling complete, so clean up and return */
    for (walker=0; walker<nwalkers_per_thread; walker++)
        run_state->currentPropDensityArray[walker] = -DBL_MAX;

    fclose(output);
    return;
}

INT4 walker_step(LALInferenceRunState *run_state,
                    LALInferenceModel *model,
                    LALInferenceVariables *current_params,
                    LALInferenceVariables *proposed_params,
                    REAL8 *current_prior, REAL8 *current_likelihood,
                    REAL8 *current_prop_density,
                    REAL8 *proposed_prior, REAL8 *proposed_likelihood,
                    REAL8 *proposed_prop_density) {
    REAL8 proposal_ratio, acceptance_probability;
    INT4 accepted = 0;

    *proposed_prior = -INFINITY;
    *proposed_likelihood = -INFINITY;

    /* Propose a new sample */
    LALInferenceClearVariables(proposed_params);

    /* Get the probability of proposing the reverse jump */
    *proposed_prop_density = *current_prop_density;
    proposal_ratio = LALInferenceStoredClusteredKDEProposal(run_state,
                                                            current_params,
                                                            proposed_params,
                                                            proposed_prop_density);

    /* Only bother calculating likelihood if within prior boundaries */
    *proposed_prior = run_state->prior(run_state, proposed_params, model);
    if (*proposed_prior > -INFINITY)
        *proposed_likelihood = run_state->likelihood(proposed_params, run_state->data, model);

    /* Find jump acceptance probability */
    acceptance_probability = (*proposed_prior + *proposed_likelihood)
                            - (*current_prior + *current_likelihood)
                            + proposal_ratio;

    /* Accept the jump with the calculated probability */
    if (acceptance_probability > 0
            || (log(gsl_rng_uniform(run_state->GSLrandom)) < acceptance_probability)) {
        LALInferenceCopyVariables(proposed_params, current_params);
        *current_prior = *proposed_prior;
        *current_likelihood = *proposed_likelihood;
        *current_prop_density = *proposed_prop_density;

        accepted = 1;
    }

    return accepted;
}


void ensemble_update(LALInferenceRunState *run_state) {
    INT4 nwalkers, nwalkers_per_thread, walker, ndim, cyclic_reflective;
    REAL8 *parameters, *samples, *param_array;

    LALInferenceVariables *algorithm_params = run_state->algorithmParams;
    nwalkers = LALInferenceGetINT4Variable(algorithm_params, "nwalkers");
    nwalkers_per_thread = LALInferenceGetINT4Variable(algorithm_params, "nwalkers_per_thread");
    cyclic_reflective = LALInferenceGetINT4Variable(algorithm_params, "cyclic_reflective");

    /* Prepare array to contain samples */
    ndim = LALInferenceGetVariableDimensionNonFixed(run_state->currentParamArray[0]);
    samples = XLALCalloc(nwalkers * ndim, sizeof(REAL8));

    /* Get this thread's walkers' locations */
    param_array = XLALCalloc(nwalkers_per_thread * ndim, sizeof(REAL8));
    for (walker = 0; walker < nwalkers_per_thread; walker++) {
        parameters = &(param_array[ndim*walker]);
        LALInferenceCopyVariablesToArray(run_state->currentParamArray[walker], parameters);
    }

    /* Send all walker locations to all MPI threads */
    MPI_Allgather(param_array, nwalkers_per_thread*ndim,
                    MPI_DOUBLE, samples, nwalkers_per_thread*ndim,
                    MPI_DOUBLE, MPI_COMM_WORLD);

    /* Update the KDE proposal */
    parallel_incremental_kmeans(run_state, samples, nwalkers, cyclic_reflective);

    /* Clean up */
    XLALFree(param_array);
    XLALFree(samples);
}


/* This is a temporary, messy solution for now.
TODO: When MPI is enables in lalinference, move this routine over and clean up */
void parallel_incremental_kmeans(LALInferenceRunState *run_state,
                                    REAL8 *samples,
                                    INT4 nwalkers,
                                    INT4 cyclic_reflective) {
    INT4 i, ndim;
    INT4 k = 0, kmax = 10;
    INT4 mpi_rank, mpi_size, best_rank;
    REAL8 bic = -INFINITY;
    REAL8 best_bic = -INFINITY;
    REAL8 *bics;

    LALInferenceKmeans *kmeans;
    LALInferenceKmeans *best_clustering = NULL;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    bics = XLALCalloc(mpi_size, sizeof(REAL8));

    ndim = LALInferenceGetVariableDimensionNonFixed(run_state->currentParamArray[0]);
    gsl_matrix_view mview = gsl_matrix_view_array(samples, nwalkers, ndim);

    /* Keep track of clustered parameter names */
    LALInferenceVariables *backward_params = XLALCalloc(1, sizeof(LALInferenceVariables));
    LALInferenceVariables *cluster_params = XLALCalloc(1, sizeof(LALInferenceVariables));
    LALInferenceVariableItem *item;
    for (item = run_state->currentParamArray[0]->head; item; item = item->next)
        if (LALInferenceCheckVariableNonFixed(run_state->currentParamArray[0], item->name))
            LALInferenceAddVariable(backward_params, item->name,
                                    item->value, item->type, item->vary);

    for (item = backward_params->head; item; item = item->next)
        LALInferenceAddVariable(cluster_params, item->name, item->value, item->type, item->vary);

    /* Have each MPI thread handle a fixed-k clustering */
    k = mpi_rank + 1;
    while (k < kmax) {
        kmeans = LALInferenceKmeansRunBestOf(k, &mview.matrix, 8, run_state->GSLrandom);
        bic = -INFINITY;
        if (kmeans)
            bic = LALInferenceKmeansBIC(kmeans);
        if (bic > best_bic) {
            if (best_clustering)
                LALInferenceKmeansDestroy(best_clustering);
            best_clustering = kmeans;
            best_bic = bic;
        } else {
            LALInferenceKmeansDestroy(kmeans);
            break;
        }

        k += mpi_size;
    }

    MPI_Gather(&best_bic, 1, MPI_DOUBLE, bics, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (mpi_rank == 0) {
        best_bic = -INFINITY;
        for (i = 0; i < mpi_size; i++) {
           if (bics[i] > best_bic) {
               best_bic = bics[i];
               best_rank = i;
           }
        }
    }

    /* Send the winning k-size to everyone */
    if (best_clustering)
        k = best_clustering->k;
    MPI_Bcast(&best_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&k, 1, MPI_INT, best_rank, MPI_COMM_WORLD);

    /* Create a kmeans instance with the winning k */
    if (mpi_rank != best_rank) {
        if (best_clustering)
            LALInferenceKmeansDestroy(best_clustering);
        best_clustering = LALInferenceKmeansRunBestOf(k, &mview.matrix, 8, run_state->GSLrandom);
    }

    /* Broadcast cluster assignments */
    MPI_Bcast(best_clustering->assignments, nwalkers, MPI_INT, best_rank, MPI_COMM_WORLD);

    /* Calculate centroids, the run to compute sizes and weights */
    LALInferenceKmeansUpdate(best_clustering);
    LALInferenceKmeansRun(best_clustering);

    LALInferenceClusteredKDE *proposal = XLALCalloc(1, sizeof(LALInferenceClusteredKDE));

    proposal->kmeans = best_clustering;

    LALInferenceInitClusteredKDEProposal(run_state, proposal,
                                            samples, nwalkers,
                                            cluster_params, "ClusteredKDEProposal",
                                            1.0, NULL, cyclic_reflective, 0);

    LALInferenceAddClusteredKDEProposalToSet(run_state, proposal);

    LALInferenceClearVariables(backward_params);
    XLALFree(backward_params);
    XLALFree(bics);
}


//-----------------------------------------
// file output routines:
//-----------------------------------------
FILE *init_ensemble_output(LALInferenceRunState *run_state,
                            INT4 verbose,
                            INT4 rank) {
    FILE *output;

    /* Print header */
    output = print_ensemble_header(run_state, rank);

    /* Extra outputs when running verbosely */
    if (verbose)
        print_proposal_header(run_state, rank);

    return output;
}


char* ensemble_output_name(LALInferenceRunState *run_state,
                            const char *out_type,
                            INT4 rank) {
    char *name = NULL;
    INT4 randomseed;

    randomseed = LALInferenceGetINT4Variable(run_state->algorithmParams, "random_seed");

    name = (char*) XLALCalloc(255, sizeof(char*));
    sprintf(name, "ensemble.%s.%u.%i", out_type, randomseed, rank);
    return name;
}


void print_samples(LALInferenceRunState *run_state,
                    FILE *output,
                    REAL8* prop_priors,
                    REAL8* prop_likelihoods,
                    REAL8* prop_densities,
                    REAL8* acceptance_rates,
                    INT4 rank) {
    REAL8 null_likelihood, timestamp, timestamp_epoch;
    REAL8 *current_priors, *current_likelihoods;
    REAL8 evidence_ratio;
    LALInferenceVariables **current_params;
    INT4 nwalkers_per_thread;
    INT4 walker, start_id, step;
    INT4 benchmark;
    struct timeval tv;

    LALInferenceVariables *algorithm_params = run_state->algorithmParams;
    nwalkers_per_thread = LALInferenceGetINT4Variable(algorithm_params, "nwalkers_per_thread");
    start_id = rank * nwalkers_per_thread;

    step = LALInferenceGetINT4Variable(run_state->algorithmParams, "step");

    null_likelihood = LALInferenceGetREAL8Variable(run_state->proposalArgs, "nullLikelihood");
    current_priors = run_state->currentPriors;
    current_likelihoods = run_state->currentLikelihoods;
    current_params = run_state->currentParamArray;

    /* Keep track of wall time if benchmarking */
    benchmark = LALInferenceGetINT4Variable(run_state->algorithmParams, "benchmark");
    if (benchmark) {
        timestamp_epoch = LALInferenceGetREAL8Variable(run_state->algorithmParams,
                                                        "timestamp_epoch");
        gettimeofday(&tv, NULL);
        timestamp = tv.tv_sec + tv.tv_usec/1E6 - timestamp_epoch;
    }


    for (walker = 0; walker < nwalkers_per_thread; walker++) {
        /* Print step number, log(posterior) */
        fprintf(output, "%d\t", step);
        fprintf(output, "%f\t",
                (current_likelihoods[walker] - null_likelihood) + current_priors[walker]);

        /* Print the non-fixed parameter values */
        LALInferencePrintSampleNonFixed(output, current_params[walker]);

        /* Print log(prior) and log(likelihood)  */
        fprintf(output, "%f\t", current_priors[walker]);
        fprintf(output, "%f\t", current_likelihoods[walker] - null_likelihood);

        if (benchmark)
            fprintf(output, "%f\t", timestamp);

        evidence_ratio = prop_priors[walker] + prop_likelihoods[walker] - prop_densities[walker];

        fprintf(output, "%f\t", evidence_ratio);

        fprintf(output, "%f\t", acceptance_rates[walker]);

        fprintf(output, "%i\t", start_id + walker);

        fprintf(output, "\n");
    }
}


void print_proposed_sample(LALInferenceRunState *run_state,
                            LALInferenceVariables *proposed_params,
                            INT4 walker,
                            INT4 accepted) {
    INT4 mpi_rank, nwalkers_per_thread;
    FILE *output = NULL;
    char *outname;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);    // This thread's index

    LALInferenceVariables *algorithm_params = run_state->algorithmParams;
    nwalkers_per_thread = LALInferenceGetINT4Variable(algorithm_params, "nwalkers_per_thread");

    outname = ensemble_output_name(run_state, "proposed", nwalkers_per_thread*mpi_rank+walker);
    output = fopen(outname, "a");

    LALInferencePrintSampleNonFixed(output, proposed_params);
    fprintf(output, "%i\n", accepted);

    fclose(output);
    XLALFree(outname);
}


void print_evidence(LALInferenceRunState *run_state,
                            FILE *output,
                            REAL8* logprior,
                            REAL8* loglike,
                            REAL8* prop_density) {
    INT4 walker, nwalkers_per_thread;
    REAL8 *ratios;
    REAL8 evidence, std = 0.0;

    LALInferenceVariables *algorithm_params = run_state->algorithmParams;
    nwalkers_per_thread = LALInferenceGetINT4Variable(algorithm_params, "nwalkers_per_thread");

    ratios = XLALCalloc(nwalkers_per_thread, sizeof(REAL8));
    for (walker = 0; walker < nwalkers_per_thread; walker++)
        ratios[walker] = logprior[walker] + loglike[walker] - prop_density[walker];

    evidence = log_add_exps(ratios, nwalkers_per_thread) - log((REAL8)nwalkers_per_thread);

    for (walker = 0; walker < nwalkers_per_thread; walker++)
        std += pow(logprior[walker] + loglike[walker] - prop_density[walker] - evidence, 2.0);

    std = sqrt(std)/nwalkers_per_thread;

    fprintf(output, "%g\t%g\t", evidence, std);

    /* Close file before returning */
    XLALFree(ratios);
}


FILE *print_ensemble_header(LALInferenceRunState *run_state, INT4 rank) {
    ProcessParamsTable *ppt;
    LALInferenceVariables **current_param_array;
    LALInferenceVariables *sample_params;
    LALInferenceIFOData *ifo_data;
    REAL8TimeSeries *time_data;
    INT4 walker, start_id;
    INT4 nwalkers_per_thread;
    INT4 ndim, int_pn_order, waveform = 0;
    INT4 nifo, randomseed, benchmark;
    REAL8 null_likelihood, normed_logl, pn_order=-1.0;
    REAL8 network_snr, sampling_rate;
    REAL8 f_ref = 0.0;
    char *outfile_name = NULL;
    FILE *output = NULL;
    char *cmd_str;

    /* Get ensemble size */
    LALInferenceVariables *algorithm_params = run_state->algorithmParams;
    nwalkers_per_thread = LALInferenceGetINT4Variable(algorithm_params, "nwalkers_per_thread");
    start_id = rank * nwalkers_per_thread;

    /* Decide on file name(s) and open */
    ppt = LALInferenceGetProcParamVal(run_state->commandLine, "--outfile");
    if (ppt) {
        outfile_name = (char*) XLALCalloc(255, sizeof(char*));
        sprintf(outfile_name, "%s.%i", ppt->value, rank);
    } else
        outfile_name = ensemble_output_name(run_state, "output", rank);

    output = fopen(outfile_name, "w");
    if (output == NULL){
        XLALErrorHandler = XLALExitErrorHandler;
        XLALPrintError("Output file error in %s, line %d. %s.\n",
                        __FILE__, __LINE__, strerror(errno));
        XLAL_ERROR_NULL(XLAL_EIO);
    }

    /* Save the command line for repoducability */
    cmd_str = LALInferencePrintCommandLine(run_state->commandLine);

    /* Use the first set of current params to get information */
    current_param_array = run_state->currentParamArray;
    sample_params = current_param_array[0];

    randomseed = LALInferenceGetINT4Variable(run_state->algorithmParams, "random_seed");
    null_likelihood = LALInferenceGetREAL8Variable(run_state->proposalArgs, "nullLikelihood");
    ndim = LALInferenceGetVariableDimensionNonFixed(sample_params);

    /* Reference frequency for evolving parameters */
    if (LALInferenceCheckVariable(sample_params, "fRef"))
        f_ref = LALInferenceGetREAL8Variable(sample_params, "fRef");

    /* Count number of detectors */
    ifo_data = run_state->data;
    nifo = 0;
    while (ifo_data){
      nifo++;
      ifo_data = ifo_data->next;
    }

    /* Integer (from an enum) identfying the waveform family used */
    if (LALInferenceCheckVariable(sample_params, "LAL_APPROXIMANT"))
        waveform = LALInferenceGetINT4Variable(sample_params, "LAL_APPROXIMANT");

    /* Determine post-Newtonian (pN) order (half of the integer stored in currentParams) */
    if (LALInferenceCheckVariable(sample_params, "LAL_PNORDER")) {
        int_pn_order = LALInferenceGetINT4Variable(sample_params, "LAL_PNORDER");
        pn_order = int_pn_order/2.0;
    }

    /* Calculated the network signal-to-noise ratio if an injection was done */
    ifo_data = run_state->data;
    network_snr = 0.0;
    while (ifo_data) {
        network_snr += ifo_data->SNR * ifo_data->SNR;
        ifo_data = ifo_data->next;
    }
    network_snr = sqrt(network_snr);

    /* Keep track of time if benchmarking */
    benchmark = LALInferenceGetINT4Variable(run_state->algorithmParams, "benchmark");

    /* Write the header information to file */
    fprintf(output,
            "  LALInference version:%s,%s,%s,%s,%s\n",
            lalAppsVCSId, lalAppsVCSDate, lalAppsVCSBranch,
            lalAppsVCSAuthor, lalAppsVCSStatus);

    fprintf(output, "  %s\n", cmd_str);

    fprintf(output, "%6s\t%20s\t%6s\t%12s\t%9s\t%9s\t%8s\t%8s\n",
        "seed", "null_likelihood", "ndet", "network_snr", "waveform", "pn_order", "ndim", "f_ref");

    fprintf(output, "%u\t%20.10lf\t%6d\t%14.6f\t%9i\t%9.1f\t%8i\t%12.1f\n",
        randomseed, null_likelihood, nifo, network_snr, waveform, pn_order, ndim, f_ref);

    /* Use time step in time-domain data to determine sampling rate */
    fprintf(output, "\n%16s\t%16s\t%10s\t%10s\t%10s\t%10s\t%20s\n",
        "Detector", "SNR", "f_low", "f_high", "start_time", "segment_length", "sampling_rate");
    ifo_data=run_state->data;
    while(ifo_data){
        time_data = ifo_data->timeData;
        sampling_rate = 1.0/time_data->deltaT;
        fprintf(output,
                "%16s\t%16.8lf\t%10.2lf\t%10.2lf\t%15.7lf\t%12d\t%10.2lf\n",
                ifo_data->detector->frDetector.name,
                ifo_data->SNR, ifo_data->fLow, ifo_data->fHigh,
                XLALGPSGetREAL8(&(ifo_data->epoch)), time_data->data->length, sampling_rate);
        ifo_data=ifo_data->next;
    }
    fprintf(output, "\n\n\n");

    /* These are the actual column headers for the samples to be output */
    fprintf(output, "cycle\tlogpost\t");

    LALInferenceFprintParameterNonFixedHeaders(output, sample_params);

    fprintf(output, "logprior\tlogl\t");

    if (benchmark)
        fprintf(output, "timestamp\t");

    fprintf(output, "\tevidence_ratio\tacceptance_rate\twalker\t");
    fprintf(output, "\n");

    /* Print starting values as 0th iteration */
    for (walker = 0; walker < nwalkers_per_thread; walker++) {
        normed_logl = run_state->currentLikelihoods[walker]-null_likelihood;
        fprintf(output, "%d\t%f\t", 0, run_state->currentPriors[walker] + normed_logl);

        LALInferencePrintSampleNonFixed(output, current_param_array[walker]);

        /* Starting prior and likelihood values */
        fprintf(output, "%f\t%f\t", run_state->currentPriors[walker],
                run_state->currentLikelihoods[walker] - null_likelihood);

        /* Elapsed dime in seconds */
        if (benchmark)
            fprintf(output, "%f\t", 0.0);

        /* Ratio used to calculate evidence */
        fprintf(output, "%f\t", 0.0);

        /* Jump proposal acceptance rate */
        fprintf(output, "%f\t", 0.0);

        /* Unique walker ID */
        fprintf(output, "%i\t", start_id + walker);

        fprintf(output,"\n");
    }

    XLALFree(outfile_name);
    XLALFree(cmd_str);
    return output;
}


void print_proposal_header(LALInferenceRunState *run_state, INT4 rank) {
    char *name = ensemble_output_name(run_state, "proposal", rank);
    FILE *output = fopen(name, "w");

    LALInferenceFprintParameterNonFixedHeaders(output, run_state->currentParamArray[0]);
    fprintf(output, "accepted\n");

    fclose(output);
    XLALFree(name);
}
