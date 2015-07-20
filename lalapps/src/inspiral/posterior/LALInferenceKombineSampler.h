/*
 *
 *  LALInferenceKombineSampler:    Ensemble Markov-Chain Monte Carlo sampler for LALInference        
 *  LALInferenceKombineSampler.h:  main header file
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

/**
 * \file LALInferenceKombineSampler.h
 * \ingroup lalapps_inspiral
 * \brief Ensemble Markov-Chain Monte Carlo sampler written for LALInference.
 *
 * Ensemble Markov-Chain Monte Carlo sampler.
 *
 */

#include <lal/LALInference.h>

/** The sampling algorithm */
void ensemble_sampler(LALInferenceRunState *run_state);


/** Evolve a walker a single step */
void walker_step(LALInferenceRunState *run_state, LALInferenceThreadState *thread,
                 REAL8 *proposed_prior, REAL8 *proposed_likelihood, REAL8 *proposed_prop_density);

/** Update the ensemble proposal from the ensemble's current state */
REAL8 get_acceptance_rate(LALInferenceRunState *run_state, REAL8 *local_acceptance_rates);

void ensemble_update(LALInferenceRunState *run_state);

void parallel_incremental_kmeans(LALInferenceRunState *run_state,
                                    REAL8 *samples,
                                    INT4 nwalkers,
                                    INT4 cyclic_reflective);

/* Data IO routines */
FILE* init_ensemble_output(LALInferenceRunState *run_state,
                            INT4 verbose,
                            INT4 rank);

void print_samples(LALInferenceRunState *run_state,
                    FILE *output,
                    REAL8* prop_priors,
                    REAL8* prop_likelihoods,
                    REAL8* prop_densities,
                    REAL8* acceptance_rates,
                    INT4 rank);

void print_evidence(LALInferenceRunState *run_state,
                            FILE *output,
                            REAL8* logprior,
                            REAL8* loglike,
                            REAL8* prop_density);

void print_proposed_sample(LALInferenceThreadState *thread);

char* ensemble_output_name(const char *out_type, INT4 rank);

FILE* print_ensemble_header(LALInferenceRunState *run_state,
                            INT4 rank);

void print_proposal_header(LALInferenceRunState *run_state,
                            INT4 rank);
