/*
*  Copyright (C) 2014 Matthew Pitkin
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
 * \file
 * \ingroup lalapps_pulsar_HeterodyneSearch
 * \author Matthew Pitkin, John Veitch, Colin Gill
 *
 * \brief Header file for the initialisation functions for the parameter estimation code for known pulsar
 * searches using the nested sampling algorithm.
 */

#ifndef _PPE_INIT_H
#define _PPE_INIT_H

#include "pulsar_parameter_estimation_nested.h"
#include "ppe_models.h"
#include "ppe_utils.h"
#include "ppe_readdata.h"

#ifdef __cplusplus
extern "C" {
#endif

void nested_sampling_algorithm_wrapper(LALInferenceRunState *runState);
void setup_live_points_array_wrapper( LALInferenceRunState *runState );
void initialise_algorithm( LALInferenceRunState *runState );
void setup_lookup_tables(LALInferenceRunState *runState, LALSource *source);
void add_initial_variables( LALInferenceVariables *ini, PulsarParameters *pars );
void initialise_prior( LALInferenceRunState *runState );
void initialise_proposal( LALInferenceRunState *runState );
void add_correlation_matrix( LALInferenceVariables *ini,
                             LALInferenceVariables *priors, REAL8Array *corMat,
                             LALStringVector *parMat );
void sum_data( LALInferenceRunState *runState );
void LogSampleToFile(LALInferenceRunState *state, LALInferenceVariables *vars);
void LogSampleToArray(LALInferenceRunState *state, LALInferenceVariables *vars);

#ifdef __cplusplus
}
#endif

#endif /* _PPE_INIT_H */
