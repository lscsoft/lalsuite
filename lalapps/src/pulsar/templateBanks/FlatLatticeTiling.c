/*
 *  Copyright (C) 2007 Karl Wette
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
 * \author K. Wette
 * \file
 * \brief Flat lattice tiling over multi-dimensioned parameter spaces
 */

#include "FlatLatticeTiling.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_sort_vector.h>

#include <lal/LALRCSID.h>
#include <lal/LALStdlib.h>
#include <lal/LALMalloc.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>

static const BOOLEAN TRUE  = (1==1);
static const BOOLEAN FALSE = (1==0);

NRCSID(FLATLATTICETILINGC, "$Id$");

static void SquareParameterSpaceBounds(gsl_vector*, INT4, gsl_vector*, REAL8*, REAL8*);

static void GetCurrentBounds(FlatLatticeTiling*, INT4, REAL8*, REAL8*);

/**
 * Creates a structure to hold information about the flat lattice tiling
 */
FlatLatticeTiling *XLALCreateFlatLatticeTiling(
					       INT4 dimension /**< [in] Dimension of the parameter space */
					       )
{

  FlatLatticeTiling *tiling = NULL;

  /* Check input */
  if (dimension < 1) {
    XLALPrintError("%s\nERROR: Dimension must be non-zero\n", FLATLATTICETILINGC);
    XLAL_ERROR_NULL("XLALCreateFlatLatticeTiling", XLAL_EINVAL);
  }

  if ((tiling = (FlatLatticeTiling*) LALMalloc(sizeof(FlatLatticeTiling))) == NULL) {
    XLALPrintError("%s\nERROR: Could not allocate a FlatLatticeTiling\n", FLATLATTICETILINGC);
    XLAL_ERROR_NULL("XLALCreateFlatLatticeTiling", XLAL_ENOMEM);
  }

  tiling->dimension = dimension;

  tiling->metric = NULL;

  tiling->mismatch = 0.0;

  tiling->generator = NULL;

  tiling->bounds = NULL;
  tiling->bounds_args = NULL;
  tiling->bounds_xml_desc = NULL;

  tiling->padding = NULL;

  tiling->increment = NULL;

  tiling->current = NULL;
  tiling->upper = NULL;

  tiling->temp = NULL;

  tiling->started = FALSE;
  tiling->finished = FALSE;

  tiling->template_count = 0;

  tiling->counted = FALSE;

  tiling->total_templates = 0;

  return tiling;

}

/**
 * Destroy a flat lattice tiling structure
 */
void XLALDestroyFlatLatticeTiling(
				  FlatLatticeTiling *tiling /**< [in] Flat lattice tiling structure */
				  )
{

  if (tiling) {

    if (tiling->metric) gsl_matrix_free(tiling->metric);

    if (tiling->generator) gsl_matrix_free(tiling->generator);

    if (tiling->bounds_args) gsl_vector_free(tiling->bounds_args);
    if (tiling->bounds_xml_desc) LALFree(tiling->bounds_xml_desc);

    if (tiling->padding) gsl_vector_free(tiling->padding);

    if (tiling->increment) gsl_matrix_free(tiling->increment);

    if (tiling->current) gsl_vector_free(tiling->current);
    if (tiling->upper) gsl_vector_free(tiling->upper);

    if (tiling->temp) gsl_vector_free(tiling->temp);

    LALFree(tiling);
    tiling = NULL;

  }

}

/**
 * Orthonormalise the columns of a matrix with respect to a metric
 */
int XLALOrthonormaliseMatrixWRTMetric(
				      gsl_matrix *matrix, /**< [in] Matrix of columns to orthonormalise */
				      gsl_matrix *metric  /**< [in] Metric to orthonormalise with respect to */
				      )
{
  
  INT4 i, j;
  INT4 m = matrix->size1;
  INT4 n = matrix->size2;
  gsl_vector *temp = NULL;
  double inner_prod = 0.0;

  /* Check input */
  if (m != (INT4)metric->size1 || m != (INT4)metric->size2) {
    XLALPrintError("%s\ERROR: Metric matrix is not the correct size\n", FLATLATTICETILINGC);
    XLAL_ERROR("XLALOrthonormaliseMatrixWRTMetric", XLAL_ESIZE);
  }

  /* Allocate */
  temp = gsl_vector_alloc(m);

  /* Orthonormalise the columns of the matrix using numerically stabilised Gram-Schmidt */
  for (i = n - 1; i >= 0; --i) {

    /* Store view of ith column */
    gsl_vector_view col_i = gsl_matrix_column(matrix, i);

    for (j = n - 1; j > i; --j) {

      /* Store view of jth column */
      gsl_vector_view col_j = gsl_matrix_column(matrix, j);

      /* Compute inner product of jth and ith columns with the metric */
      gsl_blas_dgemv(CblasNoTrans, 1.0, metric, &col_j.vector, 0.0, temp);
      gsl_blas_ddot(&col_i.vector, temp, &inner_prod);

      /* Subtract component of jth column from ith column */
      gsl_vector_memcpy(temp, &col_j.vector);
      gsl_vector_scale(temp, inner_prod);
      gsl_vector_sub(&col_i.vector, temp);
      
    }

    /* Compute inner product of ith column with itself */
    gsl_blas_dgemv(CblasNoTrans, 1.0, metric, &col_i.vector, 0.0, temp);
    gsl_blas_ddot(&col_i.vector, temp, &inner_prod);

    /* Normalise ith column */
    gsl_vector_scale(&col_i.vector, 1.0 / sqrt(inner_prod));    

  }

  /* Print debugging */
  XLALPrintInfo("Inner product of matrix columns with metric\n");
  for (i = 0; i < n; ++i) {

    /* Store view of ith column */
    gsl_vector_view col_i = gsl_matrix_column(matrix, i);

    for (j = 0; j < n; ++j) {

      /* Store view of jth column */
      gsl_vector_view col_j = gsl_matrix_column(matrix, j);

      /* Compute inner product of jth and ith columns with the metric */
      gsl_blas_dgemv(CblasNoTrans, 1.0, metric, &col_j.vector, 0.0, temp);
      gsl_blas_ddot(&col_i.vector, temp, &inner_prod);

      XLALPrintInfo("  % 0.16e", inner_prod);

    }

    XLALPrintInfo(" ;\n");

  }

  /* Cleanup */
  gsl_vector_free(temp);

  return XLAL_SUCCESS;

}

/**
 * Transform a generator matrix into a lower triangular form, and with unit covering radius
 */
gsl_matrix *XLALLowerTriangularLatticeGenerator(
						gsl_matrix *generator /**< [in] Generator matrix */
						)
{

  INT4 i, j;
  INT4 m = generator->size1;
  INT4 n = generator->size2;
  gsl_matrix *QR_decomp = NULL;
  gsl_vector *QR_tau = NULL;
  gsl_matrix *Q = NULL;
  gsl_matrix *R = NULL;
  gsl_matrix *lower_tri = NULL;

  /* Check input */
  if (m < n) {
    XLALPrintError("%s\nERROR: Generator must have at least the same number of rows as columns\n", FLATLATTICETILINGC);
    XLAL_ERROR_NULL("XLALLowerTriangularLatticeGenerator", XLAL_ESIZE);
  }

  /* Find the QR decomposition of the generator */
  QR_decomp = gsl_matrix_alloc(m, n);
  QR_tau = gsl_vector_alloc(n);
  gsl_matrix_memcpy(QR_decomp, generator);
  gsl_linalg_QR_decomp(QR_decomp, QR_tau);

  /* Unpack to get the Q and R matrices */
  Q = gsl_matrix_alloc(m, m);
  R = gsl_matrix_alloc(m, n);
  gsl_linalg_QR_unpack(QR_decomp, QR_tau, Q, R);

  /* Make signs of diagonals of R positive */
  for (i = 0; i < n; ++i) {
    gsl_vector_view R_i = gsl_matrix_column(R, i);
    gsl_vector_scale(&R_i.vector, gsl_vector_get(&R_i.vector, i) < 0 ? -1.0 : 1.0);
  }

  /* Convert R from upper to lower triangular, equivalent to changing the ordering of the dimensions */
  lower_tri = gsl_matrix_alloc(n, n);
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      gsl_matrix_set(lower_tri, i, j, gsl_matrix_get(R, n - i - 1, n - j - 1));
    }
  }
  
  /* Cleanup */
  gsl_matrix_free(QR_decomp);
  gsl_vector_free(QR_tau);
  gsl_matrix_free(Q);
  gsl_matrix_free(R);

  return lower_tri;  
  
}

/**
 * Normalise a lattice generator matrix to have unit covering radius
 */
int XLALNormaliseLatticeGenerator(
				  gsl_matrix *generator, /**< [in/out] Generator matrix of lattice */
				  REAL8 norm_thickness   /**< [in] Normalised thickness of lattice */
				  )
{
  
  gsl_matrix *LU_decomp = NULL;
  gsl_permutation *LU_perm = NULL;
  int LU_sign = 0;
  double generator_determinant = 0.0;
  double covering_radius = 0.0;

  /* Check input */
  if (generator->size1 != generator->size2) {
    XLALPrintError("%s\nERROR: Generator must be square\n", FLATLATTICETILINGC);
    XLAL_ERROR("XLALNormaliseLatticeGenerator", XLAL_ESIZE);
  }

  /* Allocate */
  LU_decomp = gsl_matrix_alloc(generator->size1, generator->size2);
  LU_perm = gsl_permutation_alloc(generator->size1);

  /* Compute generator LU decomposition */
  gsl_matrix_memcpy(LU_decomp, generator);
  gsl_linalg_LU_decomp(LU_decomp, LU_perm, &LU_sign);

  /* Compute generator determinant */
  generator_determinant = gsl_linalg_LU_det(LU_decomp, LU_sign);

  /* Compute covering radius */
  covering_radius = pow(norm_thickness * generator_determinant, 1.0 / ((double) generator->size1));

  /* Divide generator by covering radius, so covering spheres now have unit radius */
  gsl_matrix_scale(generator, 1.0 / covering_radius);

  /* Cleanup */
  gsl_matrix_free(LU_decomp);
  gsl_permutation_free(LU_perm);

  return XLAL_SUCCESS;  
  
}

/**
 * Returns the generator matrix for a $Z_n$ or cubic lattice, normalised to unit covering radius.
 * See \ref CS99, page 106.
 */
gsl_matrix *XLALCubicLatticeGenerator(
				      INT4 dimension /**< [in] Dimension of the lattice */
				      )
{

  gsl_matrix *generator = NULL;
  double norm_thickness = 0.0;
  double n = dimension;

  /* Check input */
  if (dimension < 1) {
    XLALPrintError("%s\nERROR: Dimension must be non-zero\n", FLATLATTICETILINGC);
    XLAL_ERROR_NULL("XLALCubicLatticeGenerator", XLAL_EINVAL);
  }

  /* Generator */
  generator = gsl_matrix_alloc(dimension, dimension);
  gsl_matrix_set_identity(generator);

  /* Normalised thickness */
  norm_thickness = pow(sqrt(n)/2, n);

  /* Normalise generator */
  XLALNormaliseLatticeGenerator(generator, norm_thickness);

  return generator;

}

/**
 * Returns the generator matrix for an $A_n^*$ lattice, normalised to unit covering radius.
 * $A_3^*$ is the body-centered cubic lattice. See \ref CS99, page 115.
 */
gsl_matrix *XLALAnstarLatticeGenerator(
				       INT4 dimension /**< [in] Dimension of the lattice */
				       )
{

  INT4 i, j;
  gsl_matrix *generator = NULL;
  gsl_matrix *lower_tri = NULL;
  double norm_thickness = 0.0;
  double n = dimension;

  /* Check input */
  if (dimension < 1) {
    XLALPrintError("%s\nERROR: Dimension must be non-zero\n", FLATLATTICETILINGC);
    XLAL_ERROR_NULL("XLALAnstarLatticeGenerator", XLAL_EINVAL);
  }

  /* Generator in dimension+1 space */
  generator = gsl_matrix_alloc(dimension + 1, dimension);
  gsl_matrix_set_all(generator, 0.0);
  gsl_vector_view first_row = gsl_matrix_row(generator, 0);
  gsl_vector_set_all(&first_row.vector, 1.0);
  gsl_vector_view sub_diag = gsl_matrix_subdiagonal(generator, 1);
  gsl_vector_set_all(&sub_diag.vector, -1.0);
  gsl_vector_view last_col = gsl_matrix_column(generator, dimension-1);
  gsl_vector_set_all(&last_col.vector, 1.0 / (n + 1));
  gsl_vector_set(&last_col.vector, 0, -n / (n + 1));
  XLALPrintInfo("A_%i^* generator matrix in %i-dimensional space:\n", dimension, dimension + 1);
  for (i = 0; i < dimension + 1; ++i) {
    for (j = 0; j < dimension; ++j) {
      XLALPrintInfo("  % 0.16e", gsl_matrix_get(generator, i, j));
    }
    XLALPrintInfo(" ;\n");
  }
 

  /* Generator transformed to lower triangular */
  lower_tri = XLALLowerTriangularLatticeGenerator(generator);

  /* Normalised thickness */
  norm_thickness = sqrt(n + 1)*pow((n*(n + 2))/(12*(n + 1)), n/2);

  /* Normalise generator */
  XLALNormaliseLatticeGenerator(lower_tri, norm_thickness);

  /* Cleanup */
  gsl_matrix_free(generator);

  return lower_tri;

}

/**
 * Create a parameter described by vectors of starting points and widths of each dimension.
 */
int XLALSquareParameterSpace(
			     FlatLatticeTiling* tiling, /**< [in] Flat lattice tiling structure */
			     gsl_vector *start,         /**< [in] Starting point of dimension*/
			     gsl_vector *width          /**< [in] Width of dimension */
			     )
{

  INT4 i;
  INT4 n = tiling->dimension;
  double x  = 0.0;
  double dx = 0.0;
  CHAR *p = NULL;

  /* Lower and upper bounds */
  tiling->bounds_args = gsl_vector_alloc(2 * n);
  for (i = 0; i < n; ++i) {
    x  = gsl_vector_get(start, i);
    dx = gsl_vector_get(width, i);
    gsl_vector_set(tiling->bounds_args, 2*i  , GSL_MIN_DBL(x, x + dx));
    gsl_vector_set(tiling->bounds_args, 2*i+1, GSL_MAX_DBL(x, x + dx));
  }

  /* Bound function */
  tiling->bounds = &SquareParameterSpaceBounds;

  /* XML description */
  if ((tiling->bounds_xml_desc = (CHAR*) LALMalloc((300 + 100*n) * sizeof(CHAR))) == NULL) {
    XLALPrintError("%s\nERROR: Could not allocate a CHAR*\n", FLATLATTICETILINGC);
    XLAL_ERROR("XLALSquareParameterSpace", XLAL_ENOMEM);
  }
  sprintf(tiling->bounds_xml_desc, "<type>square</type><lower>");
  for (i = 0; i < n; ++i) {
    p = tiling->bounds_xml_desc + strlen(tiling->bounds_xml_desc);
    sprintf(p, "%0.16e ", gsl_vector_get(tiling->bounds_args, 2*i));
  }
  p = tiling->bounds_xml_desc + strlen(tiling->bounds_xml_desc);
  sprintf(p, ";</lower><upper>");
  for (i = 0; i < n; ++i) {
    p = tiling->bounds_xml_desc + strlen(tiling->bounds_xml_desc);
    sprintf(p, "%0.16e ", gsl_vector_get(tiling->bounds_args, 2*i+1));
  }
  p = tiling->bounds_xml_desc + strlen(tiling->bounds_xml_desc);
  sprintf(p, ";</upper>");

  return XLAL_SUCCESS;

}

/**
 * Returns the bounds for a square parameter
 */
void SquareParameterSpaceBounds(
				gsl_vector *current, /**< [in] Current point */
				INT4 index,          /**< [in] Dimension to set bound on */
				gsl_vector *args,    /**< [in] Other arguments */
				REAL8 *lower,        /**< [out] Lower bound */
				REAL8 *upper         /**< [out] Upper bound */
				) 
{
  
  *lower = gsl_vector_get(args, 2*index  );
  *upper = gsl_vector_get(args, 2*index+1);

}

/**
 * Perform setup tasks for the flat lattice tiling
 */
int XLALSetupFlatLatticeTiling(
			       FlatLatticeTiling* tiling /**< [in] Flat lattice tiling structure */
			       )
{

  INT4 i, j;
  INT4 n = tiling->dimension;
  gsl_matrix *directions = NULL;
  gsl_matrix *LU_decomp = NULL;
  gsl_permutation *LU_perm = NULL;
  gsl_matrix *inverse = NULL;
  int LU_sign = 0;

  /* Check metric */
  if ((INT4)tiling->metric->size1 != n || (INT4)tiling->metric->size2 != n) {
    XLALPrintError("%s\nERROR: Metric matrix is not the correct size\n", FLATLATTICETILINGC);
    XLAL_ERROR("XLALSetupFlatLatticeTiling", XLAL_ESIZE);
  }
  for (i = 0; i < n; ++i) {
    for (j = i + 1; j < n; ++j) {
      if (gsl_matrix_get(tiling->metric, i, j) != gsl_matrix_get(tiling->metric, j, i)) {
	XLALPrintError("%s\nERROR: Generator matrix is not symmetric\n", FLATLATTICETILINGC);
	XLAL_ERROR("XLALSetupFlatLatticeTiling", XLAL_EDATA);
      }
    }
  }
  XLALPrintInfo("Parameter space metric:\n");
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      XLALPrintInfo("  % 0.16e", gsl_matrix_get(tiling->metric, i, j));
    }
    XLALPrintInfo(" ;\n");
  }

  /* Check mismatch */
  if (tiling->mismatch <= 0.0) {
    XLALPrintError("%s\nERROR: Mismatch must be non-zero\n", FLATLATTICETILINGC);
    XLAL_ERROR("XLALSetupFlatLatticeTiling", XLAL_EINVAL);
  }

  /* Check generator */
  if ((INT4)tiling->generator->size1 != n || (INT4)tiling->generator->size2 != n) {
    XLALPrintError("%s\nERROR: Generator matrix is not the correct size\n", FLATLATTICETILINGC);
    XLAL_ERROR("XLALSetupFlatLatticeTiling", XLAL_ESIZE);
  }
  for (i = 0; i < n; ++i) {
    for (j = i + 1; j < n; ++j) {
      if (gsl_matrix_get(tiling->generator, i, j) != 0.0) {
	XLALPrintError("%s\nERROR: Generator matrix is not upper triangular\n", FLATLATTICETILINGC);
	XLAL_ERROR("XLALSetupFlatLatticeTiling", XLAL_EDATA);
      }
    }
  }
  XLALPrintInfo("Lattice generator:\n");
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      XLALPrintInfo("  % 0.16e", gsl_matrix_get(tiling->generator, i, j));
    }
    XLALPrintInfo(" ;\n");
  }

  /* Find orthogonal directions with respect to the metric, starting from 
     the identity matrix, and scale the orthnormal directions to the mismatch. */
  directions = gsl_matrix_alloc(n, n);
  gsl_matrix_set_identity(directions);
  XLALOrthonormaliseMatrixWRTMetric(directions, tiling->metric);
  gsl_matrix_scale(directions, tiling->mismatch);
  XLALPrintInfo("Orthogonal directions with respect to metric:\n");
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      XLALPrintInfo("  % 0.16e", gsl_matrix_get(directions, i, j));
    }
    XLALPrintInfo(" ;\n");
  }

  /* Compute the increment vectors between lattice points in the parameter space. */
  tiling->increment = gsl_matrix_alloc(n, n);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, directions, tiling->generator, 0.0, tiling->increment);
  XLALPrintInfo("Increment vectors between lattice points:\n");
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      XLALPrintInfo("  % 0.16e", gsl_matrix_get(tiling->increment, i, j));
    }
    XLALPrintInfo(" ;\n");
  }

  /* Compute the inverse of the metric, and scale its columns to the 
     mismatch to find the extent of the box bounding the metric ellipses. */
  LU_decomp = gsl_matrix_alloc(n, n);
  LU_perm = gsl_permutation_alloc(n);
  inverse = gsl_matrix_alloc(n, n);
  gsl_matrix_memcpy(LU_decomp, tiling->metric);
  gsl_linalg_LU_decomp(LU_decomp, LU_perm, &LU_sign);
  gsl_linalg_LU_invert(LU_decomp, LU_perm, inverse);
  for (i = 0; i < n; ++i) {
    gsl_vector_view col_i = gsl_matrix_column(inverse, i);
    gsl_vector_scale(&col_i.vector, tiling->mismatch / sqrt(gsl_matrix_get(inverse, i, i)));
  }
  XLALPrintInfo("Inverse of the metric, scaled to the mismatch:\n");
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      XLALPrintInfo("  % 0.16e", gsl_matrix_get(inverse, i, j));
    }
    XLALPrintInfo(" ;\n");
  }

  /* Set the padding at either end along each dimension to be the size of the 
     bounding box of the ellipses, to take care of possible gaps at the edges. */
  tiling->padding = gsl_vector_alloc(n);
  gsl_vector_view inverse_diag = gsl_matrix_diagonal(inverse);
  gsl_vector_memcpy(tiling->padding, &inverse_diag.vector);
  XLALPrintInfo("Padding at either end along each dimension:\n");
  for (i = 0; i < n; ++i) {
    XLALPrintInfo("  % 0.16e", gsl_vector_get(tiling->padding, i));
  }
  XLALPrintInfo(" ;\n");

  /* We haven't started generating templates yet, nor have we counted them */
  tiling->started = FALSE;
  tiling->finished = FALSE;
  tiling->template_count = 0;
  tiling->counted = FALSE;
  tiling->total_templates = 0;

  /* Cleanup */
  gsl_matrix_free(directions);
  gsl_matrix_free(LU_decomp);
  gsl_permutation_free(LU_perm);
  gsl_matrix_free(inverse);

  return XLAL_SUCCESS;

}

/**
 * Call the bounds function to get bounds on the current point
 */
void GetCurrentBounds(
		      FlatLatticeTiling* tiling, /**< [in] Flat lattice tiling structure */
		      INT4 index,                /**< [in] Dimension to set bound on */
		      REAL8 *lower,              /**< [out] Lower bound */
		      REAL8 *upper               /**< [out] Upper bound */
		      ) 
{

  /* Padding */
  double padding = gsl_vector_get(tiling->padding, index);

  /* Call bounds function, enforcing that the function can 
     only the value of the current point for higher dimensions */
  if (index == 0) {
    (tiling->bounds)(NULL,                        index, tiling->bounds_args, lower, upper);
  }
  else {
    gsl_vector_view current_up_to_index = gsl_vector_subvector(tiling->current, 0, index);
    (tiling->bounds)(&current_up_to_index.vector, index, tiling->bounds_args, lower, upper);
  }
						
  /* Add padding if dimension is non-singular */
  if (*lower < *upper) {
    *lower -= padding;
    *upper += padding;
  }
  
}

/**
 * Find the next point in the flat lattice tiling
 */
BOOLEAN XLALNextFlatLatticePoint(
				 FlatLatticeTiling* tiling /**< [in] Flat lattice tiling structure */
				 )
{
  
  INT4 i, j;
  INT4 n = tiling->dimension;
  REAL8 lower = 0.0;
  REAL8 upper = 0.0;
  REAL8 scale = 0.0;

  /* Have we already finished? */
  if (tiling->finished) {
    return FALSE;
  }

  /* If we haven't started generating templates yet ... */ 
 if (!tiling->started) {

    /* We've started, but we're not finished */
    tiling->started = TRUE;
    tiling->finished = FALSE;
    tiling->template_count = 0;

    /* Allocate */
    tiling->current = gsl_vector_alloc(n);
    tiling->upper = gsl_vector_alloc(n);
    tiling->temp = gsl_vector_alloc(n);

    /* Calculate starting point and upper bounds */
    for (i = 0; i < n; ++i) {
      
      /* Get bounds with padding */
      GetCurrentBounds(tiling, i, &lower, &upper);

      /* Set starting point to lower bound */
      gsl_vector_set(tiling->current, i, lower);

      /* Save upper bound */
      gsl_vector_set(tiling->upper, i, upper);
      
    }

  }

  /* Otherwise advance to next template */
  else {

    for (i = n - 1; i >= 0; --i) {

      /* Get increment vector for this dimension */
      gsl_vector_view increment_i = gsl_matrix_column(tiling->increment, i);

      /* Increment current point */
      gsl_vector_add(tiling->current, &increment_i.vector);

      /* Reset the higher dimensions to their lower bounds */
      for (j = i + 1; j < n; ++j) {

	/* Get increment vector for this dimension */
	gsl_vector_view increment_j = gsl_matrix_column(tiling->increment, j);
	gsl_vector_memcpy(tiling->temp, &increment_j.vector);
	
	/* Get bounds with padding */
	GetCurrentBounds(tiling, j, &lower, &upper);

	/* Move current point in this dimension back to its lower bound */
	scale = floor((gsl_vector_get(tiling->current, j) - lower) / gsl_matrix_get(tiling->increment, j, j));
	gsl_vector_scale(tiling->temp, scale);
	gsl_vector_sub(tiling->current, tiling->temp);

	/* Save upper bound */
	gsl_vector_set(tiling->upper, j, upper);
	
      }

      /* If it is still in bounds ... */
      if (gsl_vector_get(tiling->current, i) < gsl_vector_get(tiling->upper, i)) {

	/* This is a valid point */
	break;

      }

      /* Otherwise, if this is the first dimension ... */
      else if (i == 0) {

	/* There are no more points: we are done */
	tiling->finished = TRUE;
	return FALSE;

      }

      /* Otherwise continue to next dimension */

    }

  }
    
  /* Increment template count */
  ++tiling->template_count;

  return TRUE;

}

/**
 * Return the current point in the flat lattice tiling
 */
REAL8 XLALCurrentFlatLatticePoint(
				  FlatLatticeTiling* tiling, /**< [in] Flat lattice tiling structure */
				  INT4 index                 /**< [in] Index of the dimension */
				  )
{

  return gsl_vector_get(tiling->current, index);

}

/**
 * Return the number of points generated so far by the flat lattice tiling
 */
UINT8 XLALCurrentNumberOfFlatLatticePoints(
					   FlatLatticeTiling* tiling /**< [in] Flat lattice tiling structure */
					   )
{
  
  return tiling->template_count;

}

/**
 * Count and return the total number of points generated by the flat lattice tiling
 */
UINT8 XLALTotalNumberOfFlatLatticePoints(
					 FlatLatticeTiling* tiling /**< [in] Flat lattice tiling structure */
					 )
{

  /* If we haven't already counted them ... */
  if (!tiling->counted) {

    /* Count the number of templates */
    while (XLALNextFlatLatticePoint(tiling));

    /* Save number of points */
    tiling->total_templates = tiling->template_count;
    tiling->counted = TRUE;

    /* Reset tiling structure */
    tiling->started = FALSE;
    tiling->finished = FALSE;
    tiling->template_count = 0;

  }

  return tiling->total_templates;

}
