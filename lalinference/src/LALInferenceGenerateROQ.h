/*
 *  LALInferenceROQ.h: Reduced order quadrature basis and interpolant generation
 *
 *  Copyright (C) 2014 Matthew Pitkin, Rory Smith
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
#ifndef LALINFERENCEROQ_H
#define LALINFERENCEROQ_H

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALInference.h>
#include <lal/XLALError.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>

#ifdef __cplusplus
extern "C" {
#endif


/** A structure to hold a real (double precision) interpolant matrix and interpolation node indices */
typedef struct tagLALInferenceREALROQInterpolant{
  REAL8Array *B;  /**< The interpolant matrix */
  UINT4 *nodes;   /**< The nodes (indices) for the interpolation */
}LALInferenceREALROQInterpolant;

/** A structure to hold a complex (double precision) interpolant matrix and interpolation node indices */
typedef struct tagLALInferenceCOMPLEXROQInterpolant{
  COMPLEX16Array *B;      /**< The interpolant matrix */
  UINT4 *nodes;           /**< The nodes (indices) for the interpolation */
}LALInferenceCOMPLEXROQInterpolant;

/* function to create or enrich a real orthonormal basis set from a training set of models */
REAL8 LALInferenceGenerateREAL8OrthonormalBasis(REAL8Array **RB,
                                                REAL8Vector *delta,
                                                REAL8 tolerance,
                                                REAL8Array *TS);

REAL8 LALInferenceGenerateCOMPLEX16OrthonormalBasis(COMPLEX16Array **RB,
                                                    REAL8Vector *delta,
                                                    REAL8 tolerance,
                                                    COMPLEX16Array *TS);


/* functions to test the basis */
INT4 LALInferenceTestREAL8OrthonormalBasis(REAL8Vector *delta,
                                           REAL8 tolerance,
                                           REAL8Array *RB,
                                           REAL8Array *testmodels);

INT4 LALInferenceTestCOMPLEX16OrthonormalBasis(REAL8Vector *delta,
                                               REAL8 tolerance,
                                               COMPLEX16Array *RB,
                                               COMPLEX16Array *testmodels);

/* functions to enrich the training model set and basis set */
REAL8 LALInferenceEnrichREAL8Basis(REAL8Vector *delta,
                                   REAL8 tolerance,
                                   REAL8Array **RB,
                                   REAL8Array **testmodels,
                                   REAL8Array *testmodelsnew);

REAL8 LALInferenceEnrichCOMPLEX16Basis(REAL8Vector *delta,
                                       REAL8 tolerance,
                                       COMPLEX16Array **RB,
                                       COMPLEX16Array **testmodels,
                                       COMPLEX16Array *testmodelsnew);

/* function to create the empirical interpolant */
LALInferenceREALROQInterpolant *LALInferenceGenerateREALROQInterpolant(REAL8Array *RB);

LALInferenceCOMPLEXROQInterpolant *LALInferenceGenerateCOMPLEXROQInterpolant(COMPLEX16Array *RB);

/* create ROQ weights for interpolant to calculate the linear <d|h> terms */
REAL8Vector *LALInferenceGenerateREAL8LinearWeights(REAL8Array *B, REAL8Vector *data, REAL8Vector *vars);

/* create ROQ weights for interpolant to calculate the quadratic model terms real(<h|h>) */
REAL8Vector *LALInferenceGenerateQuadraticWeights(REAL8Array *B, REAL8Vector *vars);

/* create ROQ weights for interpolant to calculate the data dot model terms */
COMPLEX16Vector *LALInferenceGenerateCOMPLEX16LinearWeights(COMPLEX16Array *B, COMPLEX16Vector *data, REAL8Vector *vars);

/* calculate ROQ version of the dot product (where the model vector is the model just computed at the interpolant points) */
REAL8 LALInferenceROQREAL8DotProduct(REAL8Vector *weights, REAL8Vector *model);
COMPLEX16 LALInferenceROQCOMPLEX16DotProduct(COMPLEX16Vector *weights, COMPLEX16Vector *model);

/* memory destruction */
void LALInferenceRemoveREALROQInterpolant( LALInferenceREALROQInterpolant *a );
void LALInferenceRemoveCOMPLEXROQInterpolant( LALInferenceCOMPLEXROQInterpolant *a );

#ifdef __cplusplus
}
#endif

#endif
