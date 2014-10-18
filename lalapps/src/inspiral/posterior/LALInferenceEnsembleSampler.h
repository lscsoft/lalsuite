/*
 *
 *  LALInferenceEnsembleSampler:    Ensemble Markov-Chain Monte Carlo sampler for LALInference        
 *  LALInferenceEnsembleSampler.h:  main header file
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
 * \file LALInferenceEnsembleSampler.h
 * \brief Ensemble Markov-Chain Monte Carlo sampler written for LALInference.
 * \ingroup LALInference
 *
 * Ensemble Markov-Chain Monte Carlo sampler.
 *
 */

#include <lal/LALInference.h>

/** The sampling algorithm */
void ensemble_sampler(LALInferenceRunState *run_state);


/** Evolve a walker a single step */
INT4 walker_step(LALInferenceRunState *run_state,
                    LALInferenceModel *model,
                    LALInferenceVariables *current_params,
                    LALInferenceVariables *proposed_params,
                    REAL8 *current_prior,
                    REAL8 *current_likelihood);


/** Update the ensemble proposal from the ensemble's current state */
void ensemble_update(LALInferenceRunState *run_state);


/* Data IO routines */
char *init_ensemble_output(LALInferenceRunState *run_state,
                            INT4 walker,
                            INT4 walker_offset,
                            INT4 verbose);

void print_ensemble_sample(LALInferenceRunState *run_state,
                            char **walker_output_names,
                            UINT4 walker);

void print_proposed_sample(LALInferenceRunState *run_state,
                            LALInferenceVariables *proposed_params,
                            INT4 walker,
                            INT4 accepted);

void print_ensemble_header(LALInferenceRunState *run_state,
                            FILE *walker_output,
                            INT4 walker);
