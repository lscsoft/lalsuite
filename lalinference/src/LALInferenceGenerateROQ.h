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
  gsl_matrix *B;  /**< The interpolant matrix */
  UINT4 *nodes;   /**< The nodes (indices) for the interpolation */
}LALInferenceREALROQInterpolant;

/** A structure to hold a complex (double precision) interpolant matrix and interpolation node indices */
typedef struct tagLALInferenceCOMPLEXROQInterpolant{
  gsl_matrix_complex *B;  /**< The interpolant matrix */
  UINT4 *nodes;           /**< The nodes (indices) for the interpolation */
}LALInferenceCOMPLEXROQInterpolant;

/* function to create a real orthonormal basis set from a training set of models */
REAL8 *LALInferenceGenerateREAL8OrthonormalBasis(gsl_vector *delta,
                                                 REAL8 tolerance,
                                                 gsl_matrix *TS,
                                                 size_t *nbases);

COMPLEX16 *LALInferenceGenerateCOMPLEX16OrthonormalBasis(gsl_vector *delta,
                                                         REAL8 tolerance,
                                                         gsl_matrix_complex *TS,
                                                         size_t *nbases);

/* functions to test the basis */
INT4 LALInferenceTestREAL8OrthonormalBasis(gsl_vector *delta,
                                           REAL8 tolerance,
                                           gsl_matrix *RB,
                                           gsl_matrix *testmodels);

INT4 LALInferenceTestCOMPLEX16OrthonormalBasis(gsl_vector *delta,
                                               REAL8 tolerance,
                                               gsl_matrix_complex *RB,
                                               gsl_matrix_complex *testmodels);

/* function to create the empirical interpolant */
LALInferenceREALROQInterpolant *LALInferenceGenerateREALROQInterpolant(gsl_matrix *RB);

LALInferenceCOMPLEXROQInterpolant *LALInferenceGenerateCOMPLEXROQInterpolant(gsl_matrix_complex *RB);

/* create ROQ weights for interpolant to calculate the data dot model terms */
gsl_vector *LALInferenceGenerateREAL8DataModelWeights(gsl_matrix *B, gsl_vector *data, gsl_vector *vars);

/* create ROQ weights for interpolant to calculate the model dot model terms */
gsl_matrix *LALInferenceGenerateREALModelModelWeights(gsl_matrix *B, gsl_vector *vars);

/* create ROQ weights for interpolant to calculate the data dot model terms */
gsl_vector_complex *LALInferenceGenerateCOMPLEX16DataModelWeights(gsl_matrix_complex *B, gsl_vector_complex *data, gsl_vector *vars);

/* create ROQ weights for interpolant to calculate the model dot model terms */
gsl_matrix_complex *LALInferenceGenerateCOMPLEXModelModelWeights(gsl_matrix_complex *B, gsl_vector *vars);

/* calculate ROQ version of the data model dot product (where the model
   vector is the model just computed at the interpolant points) */
REAL8 LALInferenceROQREAL8DataDotModel(gsl_vector *weights, gsl_vector *model);

/* calculate ROQ version of the model model dot product (where the model
   vector is the model just computed at the interpolant points) */
REAL8 LALInferenceROQREAL8ModelDotModel(gsl_matrix *weights, gsl_vector *model);

COMPLEX16 LALInferenceROQCOMPLEX16DataDotModel(gsl_vector_complex *weights, gsl_vector_complex *model);
COMPLEX16 LALInferenceROQCOMPLEX16ModelDotModel(gsl_matrix_complex *weights, gsl_vector_complex *model);

/* memory destruction */
void LALInferenceRemoveREALROQInterpolant( LALInferenceREALROQInterpolant *a );
void LALInferenceRemoveCOMPLEXROQInterpolant( LALInferenceCOMPLEXROQInterpolant *a );

#ifdef __cplusplus
}
#endif

#endif
