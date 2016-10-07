/*
*  Copyright (C) 2014, 2016 Matthew Pitkin
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
COMPLEX16Array *generate_training_set( LALInferenceRunState *rs,
                                       UINT4 n );

REAL8Array *generate_training_set_quad( LALInferenceRunState *rs,
                                        UINT4 n );

#ifdef __cplusplus
}
#endif

#endif /* _PPE_ROQ_H */
