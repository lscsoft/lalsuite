/*******************************************************************************
  Matt Pitkin - 2014

  ppe_roq.h

  Header file for ppe_roq.c

*******************************************************************************/

/*
  Author:
*/

/**
 * \file
 * \ingroup lalapps_pulsar_HeterodyneSearch
 * \author Matthew Pitkin
 *
 * \brief Header file for the reduced order quadrature generation used in parameter
 * estimation code for known pulsar searches using the nested sampling
 * algorithm.
 */

#ifndef _PPE_ROQ_H
#define _PPE_ROQ_H

#include "pulsar_parameter_estimation_nested.h"
#include "ppe_models.h"
#include "ppe_utils.h"

#define ROQTOLERANCE 1e-11

#ifdef __cplusplus
extern "C" {
#endif

/* generate Chebyshev-Gauss-Lobatto nodes in frequency */
REAL8 *chebyshev_gauss_lobatto_nodes( REAL8 freqmin, REAL8 freqmax, UINT4 nnodes );

/* generate the interpolants */
void generate_interpolant( LALInferenceRunState *runState );

/* generate a training set */
gsl_matrix_complex *generate_training_set( LALInferenceRunState *rs,
                                           UINT4 n,
                                           UINT4 freqnodes );

#ifdef __cplusplus
}
#endif

#endif /* _PPE_ROQ_H */
