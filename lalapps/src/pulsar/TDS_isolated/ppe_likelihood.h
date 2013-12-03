/*******************************************************************************
  Matt Pitkin, Colin Gill, John Veitch - 2011

  ppe_models.h

  Header file for ppe_likelihood.c

*******************************************************************************/

/*
  Author:
*/

/**
 * \file
 * \ingroup pulsarApps
 * \author Matthew Pitkin, John Veitch, Colin Gill
 *
 * \brief Header file for the likelihood and prior functions used in parameter
 * estimation code for known pulsar searches using the nested sampling
 * algorithm.
 */

#ifndef _PPE_LIKELIHOOD_H
#define _PPE_LIKELIHOOD_H

#include "pulsar_parameter_estimation_nested.h"

#ifdef __cplusplus
extern "C" {
#endif

/* likelihood function */
REAL8 pulsar_log_likelihood( LALInferenceVariables *vars, 
                             LALInferenceIFOData *data,
                             LALInferenceTemplateFunction get_pulsar_model );

/* prior function */
REAL8 priorFunction( LALInferenceRunState *runState, 
                     LALInferenceVariables *params );

/* check params are within prior range */
UINT4 in_range( LALInferenceVariables *priors, LALInferenceVariables *params );

/* convert nested samples to posterior samples */
void ns_to_posterior( LALInferenceRunState *runState );

/* create kd-tree prior */
void create_kdtree_prior( LALInferenceRunState *runState );

#ifdef __cplusplus
}
#endif

#endif /* _PPE_LIKELIHOOD_H */
