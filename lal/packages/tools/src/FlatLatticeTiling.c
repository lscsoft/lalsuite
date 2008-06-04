/*
 *  Copyright (C) 2007, 2008 Karl Wette
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
 * \author Karl Wette
 * \file
 * \brief Lattice-based template generation for flat metric parameter spaces
 */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_nan.h>

#include <lal/LALRCSID.h>
#include <lal/LALStdlib.h>
#include <lal/LALMalloc.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/FlatLatticeTiling.h>
#include <lal/FlatLatticeTilingSupport.h>

NRCSID(FLATLATTICETILINGC, "$Id$");

static const BOOLEAN TRUE  = (1==1);
static const BOOLEAN FALSE = (0==1);

/**
 * State of the tiling
 */
enum tagFlatLatticeTilingState {
  FLT_S_NotInitialised,
  FLT_S_NotStarted,
  FLT_S_InProgress,
  FLT_S_Finished
};

/**
 * Create a new flat lattice tiling structure
 */
FlatLatticeTiling *XLALCreateFlatLatticeTiling(
	 				       INT4 dimension /**< Dimension of the parameter space */
					       )
{

  FlatLatticeTiling *tiling = NULL;

  /* Check dimension */
  if (dimension < 1)
    XLAL_ERROR_NULL("'dimension' must be non-zero", XLAL_EINVAL);

  /* Allocate memory */
  if ((tiling = (FlatLatticeTiling*) LALMalloc(sizeof(FlatLatticeTiling))) == NULL)
    XLAL_ERROR_NULL("Could not allocate 'tiling'", XLAL_ENOMEM);

  /* Initialise structure */
  tiling->dimension = dimension;
  ALLOC_GSL_VECTOR_INT(tiling->bound_type, dimension, XLAL_ERROR_NULL);
  gsl_vector_int_set_all(tiling->bound_type, FLT_PSBT_Undefined);
  ALLOC_GSL_VECTOR(tiling->singular_bound, dimension, XLAL_ERROR_NULL);
  gsl_vector_set_zero(tiling->singular_bound);
  tiling->linear_bound_A = NULL;
  tiling->linear_bound_b = NULL;
  tiling->linear_tableau = NULL;
  tiling->linear_transl_b = NULL;
  tiling->linear_optimal = NULL;
  tiling->linear_temp = NULL;  
  tiling->quadratic_bound_Q = NULL;
  tiling->metric = NULL;
  tiling->max_mismatch = 0.0;
  tiling->norm_factors = NULL;
  tiling->reduced_dim = -1;
  tiling->dim_map = NULL;
  tiling->reduced_metric = NULL;
  tiling->bound_box = NULL;
  tiling->directions = NULL;
  tiling->generator = NULL;
  tiling->latt_to_param = NULL;
  tiling->state = FLT_S_NotInitialised;
  tiling->latt_current = NULL;
  tiling->param_current = NULL;
  tiling->param_offset = NULL;
  tiling->param_lower = NULL;
  tiling->param_upper = NULL;
  tiling->template_count = 0;

  return tiling;

}

/**
 * Free a flat lattice tiling structure
 */
void XLALFreeFlatLatticeTiling(
			       FlatLatticeTiling *tiling /**< Tiling structure */
			       )
{
  
  if (tiling) {

    FREE_GSL_VECTOR_INT(tiling->bound_type);
    FREE_GSL_VECTOR(tiling->singular_bound);
    FREE_GSL_MATRIX(tiling->linear_bound_A);
    FREE_GSL_VECTOR(tiling->linear_bound_b);
    XLALFreeSimplexMethodTableau(tiling->linear_tableau);
    FREE_GSL_VECTOR(tiling->linear_transl_b);
    FREE_GSL_VECTOR(tiling->linear_optimal);
    FREE_GSL_VECTOR(tiling->linear_temp);
    FREE_GSL_MATRIX(tiling->quadratic_bound_Q);
    FREE_GSL_MATRIX(tiling->metric);
    FREE_GSL_VECTOR(tiling->norm_factors);
    FREE_GSL_VECTOR_INT(tiling->dim_map);
    FREE_GSL_MATRIX(tiling->reduced_metric);
    FREE_GSL_VECTOR(tiling->bound_box);
    FREE_GSL_MATRIX(tiling->directions);
    FREE_GSL_MATRIX(tiling->generator);
    FREE_GSL_MATRIX(tiling->latt_to_param);
    FREE_GSL_VECTOR_INT(tiling->latt_current);
    FREE_GSL_VECTOR(tiling->param_current);
    FREE_GSL_VECTOR(tiling->param_offset);
    FREE_GSL_VECTOR(tiling->param_lower);
    FREE_GSL_VECTOR(tiling->param_upper);

    LALFree(tiling);
  
  }
  
}

/**
 * Set a singular bound on a dimension of the parameter space
 */
int XLALSetFlatLatticeTilingSingularBound(
					  FlatLatticeTiling *tiling, /**< Tiling structure */
					  INT4 index,                /**< Index of the singular dimension */
					  REAL8 value                /**< Value of the singular bound */
					  )
{

  const int n = tiling->dimension;

  /* Check index */
  if (index < 0 || n <= index)
    XLAL_ERROR("'index' is out of range", XLAL_EINVAL);
  switch (gsl_vector_int_get(tiling->bound_type, index)) {
  case FLT_PSBT_Undefined:
  case FLT_PSBT_Singular:
    break;
  default:
    XLAL_ERROR("bound type of 'index' has been set to a different type", XLAL_EINVAL);
  }
  
  /* Set singular bound */
  gsl_vector_int_set(tiling->bound_type, index, FLT_PSBT_Singular);
  gsl_vector_set(tiling->singular_bound, index, value);

  return XLAL_SUCCESS;

}

/**
 * Add a linear bound on a dimension of the parameter space
 */
int XLALAddFlatLatticeTilingLinearBound(
					FlatLatticeTiling *tiling, /**< Tiling structure */
					INT4 index,                /**< Index of the singular dimension */
					gsl_vector* bound_A,       /**< Linear bound matrix row */
					REAL8 bound_b              /**< Linear bound vector value */
					)
{

  const int n = tiling->dimension;

  int i, m;
  
  /* Check index */
  if (index < 0 || n <= index)
    XLAL_ERROR("'index' is out of range", XLAL_EINVAL);
  switch (gsl_vector_int_get(tiling->bound_type, index)) {
  case FLT_PSBT_Undefined:
  case FLT_PSBT_Linear:
    break;
  default:
    XLAL_ERROR("bound type of 'index' has been set to a different type", XLAL_EINVAL);
  }

  /* Set linear bound */
  if (gsl_vector_int_get(tiling->bound_type, index) == FLT_PSBT_Undefined)
    gsl_vector_int_set(tiling->bound_type, index, FLT_PSBT_Linear);

  /* Check that only linear bound dimensions are being bounded */
  for (i = 0; i < n; ++i)
    if (gsl_vector_int_get(tiling->bound_type, index) != FLT_PSBT_Linear)
      if (gsl_vector_get(bound_A, i) != 0.0)
	XLAL_ERROR("Bounds are being place on non-linear bound dimensions", XLAL_EINVAL);

  /* Deduce current size */
  m = 0;
  if (tiling->linear_bound_b)
    m = (int)tiling->linear_bound_b->size;

  /* Expand convex polytope matrix A and vector b */
  ++m;
  if (NULL == (tiling->linear_bound_A = XLALResizeGSLMatrix(tiling->linear_bound_A, m, n)))
    XLAL_ERROR("XLALResizeGSLMatrix failed", XLAL_EFAILED);
  if (NULL == (tiling->linear_bound_b = XLALResizeGSLVector(tiling->linear_bound_b, m)))
    XLAL_ERROR("XLALResizeGSLVector failed", XLAL_EFAILED);

  /* Set last row of A/b to new bounding hyperplane */
  {
    gsl_vector_view v = gsl_matrix_row(tiling->linear_bound_A, m - 1);
    gsl_vector_memcpy(&v.vector, bound_A);
  }
  gsl_vector_set(tiling->linear_bound_b, m - 1, bound_b);

  return XLAL_SUCCESS;    

}

/**
 * Set the flat lattice tiling metric and maximum mismatch
 */
int XLALSetFlatLatticeTilingMetric(
				   FlatLatticeTiling *tiling, /**< Tiling structure */
				   gsl_matrix *metric,        /**< Metric */
				   REAL8 max_mismatch,        /**< Maximum mismatch */
				   gsl_vector *norm_factors,  /**< Normalisation factors, may be NULL */
				   BOOLEAN is_metric_norm     /**< Is the given metric already normalised? */
				   )
{

  const int n = tiling->dimension;
  
  int i, j;
  gsl_matrix *norm_metric = NULL;
  
  /* Check tiling state */
  if (tiling->state != FLT_S_NotInitialised)
    XLAL_ERROR("'tiling' has already been initialised", XLAL_EFAILED);
  
  /* Check metric */
  if (!metric)
    XLAL_ERROR("'metric' must be allocated", XLAL_EINVAL);
  if (metric->size1 != metric->size2)
    XLAL_ERROR("'metric' must be square", XLAL_EINVAL);
  if (n != (int)metric->size1)
    XLAL_ERROR("'metric' size must match tiling dimension", XLAL_EINVAL);
  for (i = 0; i < n; ++i)
    for (j = 0; j < n; ++j)
      if (gsl_matrix_get(metric, i, j) != gsl_matrix_get(metric, j, i))
	XLAL_ERROR("'metric' is not symmetric", XLAL_EINVAL);
  if (tiling->metric)
    XLAL_ERROR("'tiling->metric' has already been set", XLAL_EFAILED);

  /* Check metric and maximum mismatch */
  if (max_mismatch < 0.0)
    XLAL_ERROR("'max_mismatch' must be strictly positive", XLAL_EINVAL);

  /* Allocate memory */
  ALLOC_GSL_MATRIX(tiling->metric, n, n, XLAL_ERROR);
  ALLOC_GSL_MATRIX(norm_metric,    n, n, XLAL_ERROR);
  ALLOC_GSL_VECTOR(tiling->param_current, n, XLAL_ERROR);
  ALLOC_GSL_VECTOR(tiling->param_offset,  n, XLAL_ERROR);

  /* Copy metric, mismatch */
  gsl_matrix_memcpy(tiling->metric, metric);
  gsl_matrix_memcpy(norm_metric, metric);
  tiling->max_mismatch = max_mismatch;

  /* Normalise metric and bounds if needed */
  if (norm_factors) {

    /* Allocate memory */
    ALLOC_GSL_VECTOR(tiling->norm_factors, n, XLAL_ERROR);
    
    /* Copy normalisation factors */
    gsl_vector_memcpy(tiling->norm_factors, norm_factors);

    /* Normalise or un-normalise metric */
    for (i = 0; i < n; ++i) {

      /* If metric is already normalised, un-normalise it */
      if (is_metric_norm) {
	gsl_vector_view r = gsl_matrix_row(tiling->metric, i);
	gsl_vector_view c = gsl_matrix_column(tiling->metric, i);
	gsl_vector_mul(&r.vector, tiling->norm_factors);
	gsl_vector_mul(&c.vector, tiling->norm_factors);
      }
    
      /* If metric is not normalised, normalise copy */
      else {
	gsl_vector_view r = gsl_matrix_row(norm_metric, i);
	gsl_vector_view c = gsl_matrix_column(norm_metric, i);
	gsl_vector_div(&r.vector, tiling->norm_factors);
	gsl_vector_div(&c.vector, tiling->norm_factors);
      }

    }

    /* Normalise singular bounds */
    gsl_vector_div(tiling->singular_bound, tiling->norm_factors);

    /* Normalise linear bounds so that b is +/- 1.0 */
    for (i = 0; i < (int)tiling->linear_bound_A->size1; ++i) {
      gsl_vector_view r = gsl_matrix_row(tiling->linear_bound_A, i);
      const double b = gsl_vector_get(tiling->linear_bound_b, i);
      gsl_vector_mul(&r.vector, tiling->norm_factors);
      if (b != 0.0) {
	gsl_vector_scale(&r.vector, 1.0 / fabs(b));
	gsl_vector_set(tiling->linear_bound_b, i, b < 0.0 ? -1.0 : 1.0);
      }
    }
    
  }

  /* Check all dimensions have been bounded
     and count number of non-singular dimensions */
  tiling->reduced_dim = tiling->dimension;
  for (i = 0; i < n; ++i)
    switch (gsl_vector_int_get(tiling->bound_type, i)) {
    case FLT_PSBT_Undefined:
      XLAL_ERROR("No all dimensions have been bounded", XLAL_EFAILED);
    case FLT_PSBT_Singular:
      --tiling->reduced_dim;
    }

  /* If there are non-singular dimensions */
  if (tiling->reduced_dim > 0) {

    const int r = tiling->reduced_dim;

    /* Allocate memory */
    ALLOC_GSL_VECTOR_INT(tiling->dim_map, r, XLAL_ERROR);
    ALLOC_GSL_MATRIX(tiling->reduced_metric, r, r, XLAL_ERROR);
    ALLOC_GSL_MATRIX(tiling->directions, r, r, XLAL_ERROR);
    ALLOC_GSL_MATRIX(tiling->latt_to_param, r, r, XLAL_ERROR);
    ALLOC_GSL_VECTOR_INT(tiling->latt_current, r, XLAL_ERROR);
    ALLOC_GSL_VECTOR(tiling->param_lower, r, XLAL_ERROR);
    ALLOC_GSL_VECTOR(tiling->param_upper, r, XLAL_ERROR);

    /* Make map from reduced to full dimensions */
    for (i = 0, j = 0; j < r; ++i, ++j) {
      while (gsl_vector_int_get(tiling->bound_type, i) == FLT_PSBT_Singular)
	++i;
      gsl_vector_int_set(tiling->dim_map, j, i);
    }

    /* Extract reduced metric from normalised metric */
    for (i = 0; i < r; ++i) {
      const int ii = gsl_vector_int_get(tiling->dim_map, i);
      for (j = 0; j < r; ++j) {
	const int jj = gsl_vector_int_get(tiling->dim_map, j);
	gsl_matrix_set(tiling->reduced_metric, i, j,
		       gsl_matrix_get(norm_metric, ii, jj));
      }
    }

    /* Find reduced metric ellipse bounding box */
    if (NULL == (tiling->bound_box = XLALMetricEllipseBoundingBox(tiling->reduced_metric, tiling->max_mismatch)))
      XLAL_ERROR("XLALMetricEllipseBoundingBox failed", XLAL_EFAILED);

    /* Find orthonormal directions of reduced metric */
    gsl_matrix_set_identity(tiling->directions);
    if (XLAL_SUCCESS != XLALOrthonormaliseWRTMetric(tiling->directions, tiling->reduced_metric))
      XLAL_ERROR("XLALOrthonormaliseWRTMetric failed", XLAL_EFAILED);

  }

  /* Initialise the linear boundary checking */
  if (tiling->linear_bound_A && tiling->linear_bound_b) {

    /* Create the simplex tableau */
    if (XLAL_SUCCESS != XLALQuadraticSimplexMethodTableau(-1, tiling->metric, NULL, tiling->linear_bound_A, NULL, &tiling->linear_tableau))
      XLAL_ERROR("XLALQuadraticSimplexMethodTableau failed", XLAL_EFAILED);
    
    /* Allocate memory */
    ALLOC_GSL_VECTOR(tiling->linear_transl_b, tiling->linear_tableau->m, XLAL_ERROR);
    ALLOC_GSL_VECTOR(tiling->linear_optimal, n, XLAL_ERROR);
    ALLOC_GSL_VECTOR(tiling->linear_temp, n, XLAL_ERROR);

  }
  
  /* Cleanup */
  FREE_GSL_MATRIX(norm_metric);

  return XLAL_SUCCESS;

}

/**
 * Set the tiling lattice generator
 */
static int SetGenerator(
			FlatLatticeTiling *tiling, /**< Tiling structure */
			gsl_matrix *generator,     /**< Lattice generator */
			REAL8 norm_thickness       /**< Normalised thickness of lattice */
			)
{
  
  const int r = tiling->reduced_dim;

  /* Check tiling state */
  if (tiling->state != FLT_S_NotInitialised)
    XLAL_ERROR("'tiling' has already been initialised", XLAL_EFAILED);
  
  /* If there are non-singular dimensions */
  if (r > 0) {

    /* Transform lattice generator to square lower triangular */
    if (NULL == (tiling->generator = XLALSquareLowerTriangularLatticeGenerator(generator)))
      XLAL_ERROR("XLALSquareLowerTriangularLatticeGenerator failed", XLAL_EFAILED);
    
    /* Normalise lattice generator so covering radius is sqrt(mismatch) */
    if (XLAL_SUCCESS != XLALNormaliseLatticeGenerator(tiling->generator, norm_thickness, sqrt(tiling->max_mismatch)))
      XLAL_ERROR("XLALNormaliseLatticeGenerator failed", XLAL_EFAILED);
    
    /* Compute transformation from lattice coordinates to parameter space */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tiling->directions, tiling->generator, 0.0, tiling->latt_to_param);
    
  }

  /* Tiling is now fully initialised */
  tiling->state = FLT_S_NotStarted;

  return XLAL_SUCCESS;

}

/**
 * Set the tiling to a cubic lattice generator
 */
int XLALSetCubicLatticeTiling(
			      FlatLatticeTiling *tiling /**< Tiling structure */
			      )
{

  const int r = tiling->reduced_dim;
  gsl_matrix *generator = NULL;
  REAL8 norm_thickness = 0.0;

  /* Allocate memory */
  ALLOC_GSL_MATRIX(generator, r, r, XLAL_ERROR);

  /* Create generator */
  gsl_matrix_set_identity(generator);

  /* Calculate normalised thickness */
  norm_thickness = pow(sqrt(r)/2, r);

  /* Set lattice generator */
  if (SetGenerator(tiling, generator, norm_thickness) != XLAL_SUCCESS)
    XLAL_ERROR("SetGenerator failed", XLAL_EFAILED);

  /* Cleanup */
  FREE_GSL_MATRIX(generator);

  return XLAL_SUCCESS;

}

/**
 * Set the tiling to an \fA_n^*\f lattice generator
 */
int XLALSetAnstarLatticeTiling(
			      FlatLatticeTiling *tiling /**< Tiling structure */
			      )
{

  const int r = tiling->reduced_dim;
  gsl_matrix *generator = NULL;
  REAL8 norm_thickness = 0.0;

  /* Allocate memory */
  ALLOC_GSL_MATRIX(generator, r + 1, r, XLAL_ERROR);

  /* Create generator in (r + 1) space */
  gsl_matrix_set_all(generator, 0.0);
  {
    gsl_vector_view first_row = gsl_matrix_row(generator, 0);
    gsl_vector_view sub_diag = gsl_matrix_subdiagonal(generator, 1);
    gsl_vector_view last_col = gsl_matrix_column(generator, r - 1);
    gsl_vector_set_all(&first_row.vector, 1.0);
    gsl_vector_set_all(&sub_diag.vector, -1.0);
    gsl_vector_set_all(&last_col.vector, 1.0 / (r + 1.0));
    gsl_vector_set(&last_col.vector, 0, -r * 1.0 / (r + 1.0));
  }

  /* Calculate normalised thickness */
  norm_thickness = sqrt(r + 1.0)*pow((1.0*r*(r + 2))/(12.0*(r + 1)), 0.5*r);

  /* Set lattice generator */
  if (SetGenerator(tiling, generator, norm_thickness) != XLAL_SUCCESS)
    XLAL_ERROR("SetGenerator failed", XLAL_EFAILED);

  /* Cleanup */
  FREE_GSL_MATRIX(generator);

  return XLAL_SUCCESS;

}

/**
 * Get bounds on the specified dimension
 */
static int GetBounds(
		     FlatLatticeTiling *tiling, /**< Tiling structure */
		     gsl_vector *current,       /**< Vector to set bounds on */
		     INT4 index,                /**< Index of the dimension to bound */
		     REAL8 *lower,              /**< Lower bound on dimension */
		     REAL8 *upper               /**< Upper bound on dimension */
		     )
{
  
  const int n = tiling->dimension;
  const int bound_type = gsl_vector_int_get(tiling->bound_type, index);

  /* Switch on type of bound */
  switch (bound_type) {
  case FLT_PSBT_Singular:

    *lower = *upper = gsl_vector_get(tiling->singular_bound, index);
    break;

  case FLT_PSBT_Linear:
    {

      int i, j;
      double ax_lt_index = 0.0;
      double limit = 0.0;

      /* Starting bounds */
      *upper = GSL_NEGINF;
      *lower = GSL_POSINF;

      /* Consider all bounds */
      for (i = 0; i < (int)tiling->linear_bound_A->size1; ++i) {

	/* Only use rows which are non-zero at index */
	const double aindex = gsl_matrix_get(tiling->linear_bound_A, i, index);
	if (aindex == 0.0)
	  continue;

	/* Only use rows which are all zero in higher (>index) dimensions */
	for (j = index + 1; j < n; ++j)
	  if (gsl_matrix_get(tiling->linear_bound_A, i, j) != 0.0) {
	    j = -1;
	    break;
	  }
	if (j < 0)
	  continue;
	
	/* Get dot product of current point with lower (<index) bounds */
	if (index > 0) {
	  gsl_vector_view row = gsl_matrix_row(tiling->linear_bound_A, i);
	  gsl_vector_view u = gsl_vector_subvector(&row.vector, 0, index);
	  gsl_vector_view v = gsl_vector_subvector(current, 0, index);
	  gsl_blas_ddot(&u.vector, &v.vector, &ax_lt_index);
	}
	else
	  ax_lt_index = 0.0;

	/* Get right hand side, subtract dot product and divide by value at index */
	limit = (gsl_vector_get(tiling->linear_bound_b, i) - ax_lt_index) / aindex;

	/* If coeffficient is positive, set upper bound */
	if (aindex > 0) {
	  if (limit > *upper)
	    *upper = limit;
	}
	
	/* If coefficient is negative, set lower bound */
	else {
	  if (limit < *lower)
	    *lower = limit;
	}

      }

    }
    break;

  default:
    XLAL_ERROR("Invalid boundary type", XLAL_EINVAL);

  }

  /* Check if bounds are finite and consistent */
  if (!gsl_finite(*lower) || !gsl_finite(*upper) || *lower > *upper)
    XLAL_ERROR("Boundaries are erroneous", XLAL_EINVAL);

  return XLAL_SUCCESS;

}

/**
 * Update the current point in parameter space from 
 * the current point in lattice coordinates
 */
static void UpdateParamCurrent(
			       FlatLatticeTiling *tiling /**< Tiling structure */
			       )
{
  
  const int r = tiling->reduced_dim;

  int i, j;
  
  /* Compute the matrix transform for non-singular dimensions */
  gsl_vector_memcpy(tiling->param_current, tiling->param_offset);
  for (i = 0; i < r; ++i) {
    const int ii = gsl_vector_int_get(tiling->dim_map, i);
    double x = gsl_vector_get(tiling->param_current, ii);
    for (j = 0; j < r; ++j) {
      x += gsl_matrix_get(tiling->latt_to_param, i, j) * 
	gsl_vector_int_get(tiling->latt_current, j);
    }
    gsl_vector_set(tiling->param_current, ii, x);
  }

}

/**
 * Check whether the metric ellipse about the
 * current point intersects the parameter space.
 */
static BOOLEAN DoesMetricEllipseIntersectParameterSpace(
							FlatLatticeTiling *tiling /**< Tiling structure */
							)
{

  const int r = tiling->reduced_dim;
  
  int i;
  
  /* If point is in bounds, metric ellipse
     definitely intersects parameter space */
  for (i = 0; i < r; ++i) {
    const int ii = gsl_vector_int_get(tiling->dim_map, i);
    const double lower = gsl_vector_get(tiling->param_lower, i);
    const double point = gsl_vector_get(tiling->param_current, ii);
    const double upper = gsl_vector_get(tiling->param_upper, i);
    if (point < lower || upper < point)
      break;
  }
  if (i >= r)
    return TRUE;

  /* Check the linear bounds:
     Solve the quadratic simplex problem to get the smallest
     distance from the current point to the parameter space. */
  if (tiling->linear_bound_A && tiling->linear_bound_b) {

    double dist = 0.0;

    /* Shift parameter space relative to current point */
    gsl_vector_memcpy(tiling->linear_transl_b, tiling->linear_bound_b);
    gsl_blas_dgemv(CblasNoTrans, -1.0, tiling->linear_bound_A, tiling->param_current, 1.0, tiling->linear_transl_b);

    /* Update the quadratic problem simplex */
    if (XLAL_SUCCESS != XLALQuadraticSimplexMethodTableau(0, NULL, NULL, NULL, tiling->linear_transl_b, &tiling->linear_tableau))
      XLAL_ERROR("XLALQuadraticSimplexMethodTableau failed", XLAL_EFAILED);
    
    /* Solve to find the smallest distance */
    if (XLAL_SUCCESS != XLALQuadraticSimplexMethod(tiling->linear_tableau, tiling->linear_optimal))
      return FALSE;

    /* Calculate the distance with respect to the metric */
    gsl_blas_dgemv(CblasNoTrans, 1.0, tiling->metric, tiling->linear_optimal, 0.0, tiling->linear_temp);
    gsl_blas_ddot(tiling->linear_optimal, tiling->linear_temp, &dist);

    /* If distance is greater than mismatch, no intersection */
    if (dist > tiling->max_mismatch)
      return FALSE;
    
  }

  return TRUE;

}

/**
 * Move to the next point in the flat lattice tiling
 */
int XLALNextFlatLatticeTile(
			    FlatLatticeTiling *tiling /**< Tiling structure */
			    )
{

  const int r = tiling->reduced_dim;
  REAL8 lower = 0.0;
  REAL8 upper = 0.0;
  REAL8 dist = 0.0;
  int dlatt = 0;

  int i, j;

  /* Switch on tiling state */
  switch (tiling->state) {

  case FLT_S_NotInitialised:
    /* Fail if uninitialised */
    XLAL_ERROR("'tiling' has not been fully initialised", XLAL_EINVAL);
    
  case FLT_S_Finished:
    /* Exit if finished */
    return XLAL_FAILURE;
    
  case FLT_S_NotStarted:
    /* If not started */
      
    /* If all dimensions singular, return the single lattice point */
    if (r == 0) {
      gsl_vector_memcpy(tiling->param_current, tiling->singular_bound);
      tiling->state = FLT_S_Finished;
      return XLAL_FAILURE;
    }
    
    /* Find lowest point to use as overall offset,
       and store starting bounds on current point. */
    gsl_vector_memcpy(tiling->param_offset, tiling->singular_bound);
    for (i = 0; i < r; ++i) {

      /* Get mapped index and metric ellipse bounding box */
      const int ii = gsl_vector_int_get(tiling->dim_map, i);
      const double bound_box = gsl_vector_get(tiling->bound_box, i);
      
      /* Get current bounds on the offset and choose the lowest */
      if (XLAL_SUCCESS != GetBounds(tiling, tiling->param_offset, ii, &lower, &upper))
	XLAL_ERROR("GetBounds failed", XLAL_EFAILED);
      gsl_vector_set(tiling->param_offset, ii, lower - bound_box);

      /* Store current bounds */
      gsl_vector_set(tiling->param_lower, i, lower);
      gsl_vector_set(tiling->param_upper, i, upper);

    }
    
    /* Initialise current point */
    gsl_vector_int_set_zero(tiling->latt_current);
    gsl_vector_int_set(tiling->latt_current, r - 1, -1);
    UpdateParamCurrent(tiling);
    
    /* Tiling has started */
    tiling->template_count = 0;
    tiling->state = FLT_S_InProgress;
    
  }
  
  /* Loop until found point in parameter space. This 
     will only happen at the parameter space boundaries. */
  do {

    /* Loop over dimensions starting from highest */
    for (i = r - 1; i >= 0; --i) {

      /* Get mapped index and metric ellipse bounding box */
      const int ii = gsl_vector_int_get(tiling->dim_map, i);
      const double bound_box_ii = gsl_vector_get(tiling->bound_box, i);

      /* Increment lattice point along this dimension */
      gsl_vector_int_set(tiling->latt_current, i, 1 +
			 gsl_vector_int_get(tiling->latt_current, i));
      UpdateParamCurrent(tiling);

      /* Move point back to lower bound in higher dimensions */
      for (j = i + 1; j < r; ++j) {

	/* Get mapped index and metric ellipse bounding box */
	const int jj = gsl_vector_int_get(tiling->dim_map, j);
	const double bound_box_jj = gsl_vector_get(tiling->bound_box, j);

	/* Get bounds on this dimension */
	if (XLAL_SUCCESS != GetBounds(tiling, tiling->param_current, jj, &lower, &upper))
	  XLAL_ERROR("GetBounds failed", XLAL_EFAILED);

	/* Store current bounds */
	gsl_vector_set(tiling->param_lower, j, lower);
	gsl_vector_set(tiling->param_upper, j, upper);

	/* Move the current point past the lower bound */
	dist = lower - bound_box_jj - gsl_vector_get(tiling->param_current, jj);
	dlatt = (int) ceil(dist / gsl_matrix_get(tiling->latt_to_param, j, j));
	gsl_vector_int_set(tiling->latt_current, j, dlatt +
			   gsl_vector_int_get(tiling->latt_current, j));
	UpdateParamCurrent(tiling);

      }

      /* If point is within bounds, stop moving point */
      if (gsl_vector_get(tiling->param_current, ii) <=
	  gsl_vector_get(tiling->param_upper, i) + bound_box_ii)
	break;

      /* If this is lowest dimension, we're done! */
      if (i == 0) {
	tiling->state = FLT_S_Finished;
	return XLAL_FAILURE;
      }

    }

  } while (!DoesMetricEllipseIntersectParameterSpace(tiling));

  /* Template found, increment count */
  ++tiling->template_count;

  /* Un-normalise template point */
  if (tiling->norm_factors)
    gsl_vector_mul(tiling->param_current, tiling->norm_factors);

  return XLAL_SUCCESS;

}
