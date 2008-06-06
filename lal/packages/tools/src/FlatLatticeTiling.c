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
  tiling->linear_vertices = NULL;
  tiling->linear_map = NULL;
  tiling->linear_min_vertex = NULL;
  tiling->linear_max_vertex = NULL;
  tiling->linear_tableau = NULL;
  tiling->linear_current = NULL;
  tiling->linear_transl_b = NULL;
  tiling->linear_optimal = NULL;
  tiling->linear_temp = NULL;  
  tiling->quadratic_bound_Q = NULL;
  tiling->quadratic_map = NULL;
  tiling->bound_inverse_map = NULL;
  ALLOC_GSL_VECTOR_INT(tiling->careful_linear, dimension, XLAL_ERROR_NULL);
  gsl_vector_int_set_all(tiling->careful_linear, TRUE);
  ALLOC_GSL_VECTOR_INT(tiling->careful_quadratic, dimension, XLAL_ERROR_NULL);
  gsl_vector_int_set_all(tiling->careful_quadratic, TRUE);
  tiling->norm_metric = NULL;
  tiling->norm_to_real_mul = NULL;
  tiling->norm_to_real_add = NULL;
  tiling->max_mismatch = 0.0;
  tiling->reduced_dim = -1;
  tiling->reduced_map = NULL;
  tiling->latt_to_norm = NULL;
  tiling->norm_padding = NULL;
  tiling->latt_current = NULL;
  tiling->norm_current = NULL;
  tiling->norm_lower = NULL;
  tiling->norm_upper = NULL;
  tiling->current_tile = NULL;
  tiling->tile_count = 0;
  tiling->state = FLT_S_NotInitialised;

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
    FREE_GSL_MATRIX(tiling->linear_vertices);
    FREE_GSL_VECTOR_INT(tiling->linear_map);
    FREE_GSL_VECTOR(tiling->linear_min_vertex);
    FREE_GSL_VECTOR(tiling->linear_max_vertex);
    XLALFreeSimplexMethodTableau(tiling->linear_tableau);
    FREE_GSL_VECTOR(tiling->linear_current);
    FREE_GSL_VECTOR(tiling->linear_transl_b);
    FREE_GSL_VECTOR(tiling->linear_optimal);
    FREE_GSL_VECTOR(tiling->linear_temp);
    FREE_GSL_MATRIX(tiling->quadratic_bound_Q);
    FREE_GSL_VECTOR_INT(tiling->quadratic_map);
    FREE_GSL_VECTOR_INT(tiling->bound_inverse_map);
    FREE_GSL_VECTOR_INT(tiling->careful_linear);
    FREE_GSL_VECTOR_INT(tiling->careful_quadratic);
    FREE_GSL_MATRIX(tiling->norm_metric);
    FREE_GSL_VECTOR(tiling->norm_to_real_mul);
    FREE_GSL_VECTOR(tiling->norm_to_real_add);
    FREE_GSL_VECTOR_INT(tiling->reduced_map);
    FREE_GSL_MATRIX(tiling->latt_to_norm);
    FREE_GSL_VECTOR(tiling->norm_padding);
    FREE_GSL_VECTOR_INT(tiling->latt_current);
    FREE_GSL_VECTOR(tiling->norm_current);
    FREE_GSL_VECTOR(tiling->norm_lower);
    FREE_GSL_VECTOR(tiling->norm_upper);

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
  gsl_vector_int_set(tiling->careful_linear, index, FALSE);
  gsl_vector_int_set(tiling->careful_quadratic, index, FALSE);

  return XLAL_SUCCESS;

}

/**
 * Add a linear bound on a dimension of the parameter space
 */
int XLALAddFlatLatticeTilingLinearBound(
					FlatLatticeTiling *tiling, /**< Tiling structure */
					INT4 index,                /**< Index of the singular dimension */
					gsl_vector* bound_A,       /**< Linear bound matrix row */
					REAL8 bound_b,             /**< Linear bound vector value */
					BOOLEAN careful            /**< Do careful boundary checking on this dimension? */
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

  /* Set careful boundary checking */
  gsl_vector_int_set(tiling->careful_linear, index, careful);

  return XLAL_SUCCESS;    

}

/**
 * Finalise the bounds
 */
static int FinaliseBounds(
			  FlatLatticeTiling *tiling /**< Tiling structure */
			  )
{

  const int n = tiling->dimension;

  /* Allocate memory */
  ALLOC_GSL_VECTOR_INT(tiling->bound_inverse_map, n, XLAL_ERROR);
  gsl_vector_int_set_all(tiling->bound_inverse_map, -1);

  /* Linear bounds: remove non-bound dimensions,
     create map and find minimum and maximum vertices */
  if (tiling->linear_bound_A && tiling->linear_bound_b) {

    int i, ln;
    gsl_matrix *old_linear_bound_A = NULL;

    /* Save old matrix */
    old_linear_bound_A = tiling->linear_bound_A;
    tiling->linear_bound_A = NULL;

    /* Count number of bound dimensions */
    ln = 0;
    for (i = 0; i < n; ++i)
      if (gsl_vector_int_get(tiling->bound_type, i) == FLT_PSBT_Linear)
	++ln;

    /* Allocate memory */
    ALLOC_GSL_MATRIX(tiling->linear_bound_A, old_linear_bound_A->size1, ln, XLAL_ERROR);
    ALLOC_GSL_VECTOR_INT(tiling->linear_map, ln, XLAL_ERROR);
    ALLOC_GSL_VECTOR(tiling->linear_min_vertex, ln, XLAL_ERROR);
    ALLOC_GSL_VECTOR(tiling->linear_max_vertex, ln, XLAL_ERROR);
    
    /* Create (inverse)map from bound to full dimensions */
    {
      int ii;
      for (i = 0, ii = 0; i < ln; ++i, ++ii) {
	while (gsl_vector_int_get(tiling->bound_type, ii) != FLT_PSBT_Linear)
	  ++ii;
	gsl_vector_int_set(tiling->linear_map, i, ii);
	if (gsl_vector_int_get(tiling->bound_inverse_map, ii) != -1)
	  XLAL_ERROR("Inverse lookup for bound index maps is corrupted", XLAL_EFAILED);
	gsl_vector_int_set(tiling->bound_inverse_map, ii, i);
      }
    }

    /* Copy columns across */
    for (i = 0; i < ln; ++i) {
      const int ii = gsl_vector_int_get(tiling->linear_map, i);
      gsl_vector_view old = gsl_matrix_column(old_linear_bound_A, ii);
      gsl_vector_view new = gsl_matrix_column(tiling->linear_bound_A, i);
      gsl_vector_memcpy(&new.vector, &old.vector);
    }

    /* Get vertices of polytope */
    if (XLAL_SUCCESS != XLALPolytopeVertices(tiling->linear_bound_A, tiling->linear_bound_b, &tiling->linear_vertices))
      XLAL_ERROR("XLALPolytopeVertices failed", XLAL_EFAILED);

    /* Find minimum and maximum on each dimension */
    for (i = 0; i < ln; ++i) {
      gsl_vector_view c = gsl_matrix_column(tiling->linear_vertices, i);
      gsl_vector_set(tiling->linear_min_vertex, i, gsl_vector_min(&c.vector));
      gsl_vector_set(tiling->linear_max_vertex, i, gsl_vector_max(&c.vector));
    }

    /* Cleanup */
    FREE_GSL_MATRIX(old_linear_bound_A);

  }

  return XLAL_SUCCESS;

}

/**
 * Normalise the bounds
 */
static int NormaliseBounds(
			   FlatLatticeTiling *tiling /**< Tiling structure */
			   )
{

  /* Linear bounds: normalise using the normalised to real conversions in reverse */
  if (tiling->linear_bound_A && tiling->linear_bound_b) {
    
    const int ln = (int)tiling->linear_map->size;

    int i, j;

    /* Subtract linear_bound_A*norm_to_real_add off linear_bound_b */
    for (i = 0; i < (int)tiling->linear_bound_b->size; ++i) {
      double b = gsl_vector_get(tiling->linear_bound_b, i);
      for (j = 0; j < ln; ++j) {
	const int jj = gsl_vector_int_get(tiling->linear_map, j);
	b -= gsl_matrix_get(tiling->linear_bound_A, i, j) * gsl_vector_get(tiling->norm_to_real_add, jj);
      }
      gsl_vector_set(tiling->linear_bound_b, i, b);
    }
    
    /* Normalise rows of linear_bound_A */
    for (i = 0; i < (int)tiling->linear_bound_A->size1; ++i) {
      for (j = 0; j < ln; ++j) {
	const int jj = gsl_vector_int_get(tiling->linear_map, j);
	gsl_matrix_set(tiling->linear_bound_A, i, j,
		       gsl_matrix_get(tiling->linear_bound_A, i, j) *
		       gsl_vector_get(tiling->norm_to_real_mul, jj));
      }
    }
    
    /* Normalise linear_bound_A/b so that |elements of linear_bound_b| = 1 */
    for (i = 0; i < (int)tiling->linear_bound_A->size1; ++i) {
      const double b = gsl_vector_get(tiling->linear_bound_b, i);
      gsl_vector_view r = gsl_matrix_row(tiling->linear_bound_A, i);
      if (b != 0.0) {
	gsl_vector_scale(&r.vector, 1.0 / fabs(b));
	gsl_vector_set(tiling->linear_bound_b, i, b < 0.0 ? -1.0 : 1.0);
      }
    }

  }

  return XLAL_SUCCESS;

}

/**
 * Initialise the bound checking routines
 */
static int InitialiseBoundChecking(
				   FlatLatticeTiling *tiling /**< Tiling structure */
				   )
{
  
  /* Linear bounds: initialise the simplex tableau */
  if (tiling->linear_bound_A && tiling->linear_bound_b) {

    const int ln = (int)tiling->linear_map->size;

    int i, j;

    /* Allocate memory */
    ALLOC_GSL_MATRIX(tiling->linear_metric, ln, ln, XLAL_ERROR);
    ALLOC_GSL_VECTOR(tiling->linear_current,    ln, XLAL_ERROR);

    /* Extract linear bound metric from normalised metric */
    for (i = 0; i < ln; ++i) {
      const int ii = gsl_vector_int_get(tiling->linear_map, i);
      for (j = 0; j < ln; ++j) {
	const int jj = gsl_vector_int_get(tiling->linear_map, j);
	gsl_matrix_set(tiling->linear_metric, i, j, gsl_matrix_get(tiling->norm_metric, ii, jj));
      }
    }

    /* Create the simplex tableau */
    if (XLAL_SUCCESS != XLALQuadraticSimplexMethodTableau(-1, tiling->linear_metric, NULL, tiling->linear_bound_A, NULL, &tiling->linear_tableau))
      XLAL_ERROR("XLALQuadraticSimplexMethodTableau failed", XLAL_EFAILED);
    
    /* Allocate memory */
    ALLOC_GSL_VECTOR(tiling->linear_transl_b, tiling->linear_tableau->m, XLAL_ERROR);
    ALLOC_GSL_VECTOR(tiling->linear_optimal,  ln, XLAL_ERROR);
    ALLOC_GSL_VECTOR(tiling->linear_temp,     ln, XLAL_ERROR);

  }

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
  
  const int bound_type = gsl_vector_int_get(tiling->bound_type, index);
  const int rindex = gsl_vector_int_get(tiling->bound_inverse_map, index);

  /* Switch on type of bound */
  switch (bound_type) {
  case FLT_PSBT_Singular:

    *lower = *upper = gsl_vector_get(tiling->singular_bound, index);
    break;

  case FLT_PSBT_Linear:
    {

      const int ln = (int)tiling->linear_map->size;

      int i, j;
      double ax_lt_rindex = 0.0;
      double limit = 0.0;

      /* Starting bounds are minimum/maximum vertices */
      *lower = GSL_POSINF; /* gsl_vector_get(tiling->linear_min_vertex, rindex); */
      *upper = GSL_NEGINF; /* gsl_vector_get(tiling->linear_max_vertex, rindex); */

      /* Consider all bounds */
      for (i = 0; i < (int)tiling->linear_bound_A->size1; ++i) {

	/* Only use rows which are non-zero at index */
	const double a_rindex = gsl_matrix_get(tiling->linear_bound_A, i, rindex);
	if (a_rindex == 0.0)
	  continue;

	/* Only use rows which are all zero in higher (>index) dimensions */
	for (j = rindex + 1; j < ln; ++j)
	  if (gsl_matrix_get(tiling->linear_bound_A, i, j) != 0.0)
	    break;
	if (j < ln)
	  continue;
	
	/* Get dot product of current point with lower (<index) bounds */
	ax_lt_rindex = 0.0;
	for (j = 0; j < rindex; ++j) {
	  const int jj = gsl_vector_int_get(tiling->linear_map, j);
	  ax_lt_rindex += gsl_vector_get(current, jj) * 
	    gsl_matrix_get(tiling->linear_bound_A, i, j);
	}
	
	/* Get right hand side, subtract dot product and divide by value at index */
	limit = (gsl_vector_get(tiling->linear_bound_b, i) - ax_lt_rindex) / a_rindex;

	/* If coeffficient is positive, set upper bound */
	if (a_rindex > 0) {
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
 * Set the flat lattice tiling metric and maximum mismatch, and perform other initialisation tasks
 */
int XLALSetFlatLatticeTilingMetric(
				   FlatLatticeTiling *tiling, /**< Tiling structure */
				   gsl_matrix *metric,        /**< Metric */
				   REAL8 max_mismatch,        /**< Maximum mismatch */
				   gsl_vector *norm_to_real   /**< Multiply to get real metric, may be NULL */
				   )
{

  const int n = tiling->dimension;
  
  int i, j;
  
  /* Check tiling state */
  if (tiling->state != FLT_S_NotInitialised)
    XLAL_ERROR("'tiling' has already been initialised", XLAL_EFAILED);

  /* Finalise bounds */
  if (XLAL_SUCCESS != FinaliseBounds(tiling))
    XLAL_ERROR("FinaliseBounds failed", XLAL_EFAILED);
  
  /* Check metric and mismatch */
  if (tiling->norm_metric)
    XLAL_ERROR("'tiling->norm_metric' has already been set", XLAL_EFAILED);
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
  if (max_mismatch < 0.0)
    XLAL_ERROR("'max_mismatch' must be strictly positive", XLAL_EINVAL);

  /* Allocate memory */
  ALLOC_GSL_MATRIX(tiling->norm_metric,      n, n, XLAL_ERROR);
  ALLOC_GSL_VECTOR(tiling->norm_to_real_mul,    n, XLAL_ERROR);
  ALLOC_GSL_VECTOR(tiling->norm_to_real_add,    n, XLAL_ERROR);
  ALLOC_GSL_VECTOR(tiling->norm_current,        n, XLAL_ERROR);

  /* Copy metric, normalised to real conversion, and mismatch */
  gsl_matrix_memcpy(tiling->norm_metric, metric);
  if (norm_to_real)
    gsl_vector_memcpy(tiling->norm_to_real_mul, norm_to_real);
  else
    gsl_vector_set_all(tiling->norm_to_real_mul, 1.0);
  tiling->max_mismatch = max_mismatch;

  /* Check all dimensions have been bounded
     and count number of non-singular dimensions */
  tiling->reduced_dim = tiling->dimension;
  for (i = 0; i < n; ++i) {
    switch (gsl_vector_int_get(tiling->bound_type, i)) {
    case FLT_PSBT_Undefined:
      XLAL_ERROR("No all dimensions have been bounded", XLAL_EFAILED);
    case FLT_PSBT_Singular:
      --tiling->reduced_dim;
    }
  }

  /* Set the normalised to real offset to be
     the midpoint of the bounds on each dimension */
  gsl_vector_set_zero(tiling->norm_to_real_add);
  for (i = 0; i < n; ++i) {
    double lower = 0.0;
    double upper = 0.0;
    if (XLAL_SUCCESS != GetBounds(tiling, tiling->norm_to_real_add, i, &lower, &upper))
      XLAL_ERROR("GetBounds failed", XLAL_EFAILED);
    gsl_vector_set(tiling->norm_to_real_add, i, 0.5*(lower + upper));
  }

  return XLAL_SUCCESS;

}

/**
 * Set the tiling lattice generator, and perform other initialisation tasks
 */
static int SetTilingLattice(
			    FlatLatticeTiling *tiling, /**< Tiling structure */
			    gsl_matrix *generator,     /**< Lattice generator */
			    REAL8 norm_thickness       /**< Normalised thickness of lattice */
			    )
{
  
  const int r = tiling->reduced_dim;

  /* Check tiling state */
  if (tiling->state != FLT_S_NotInitialised)
    XLAL_ERROR("'tiling' has already been initialised", XLAL_EFAILED);
  
  /* Normalise bounds */
  if (XLAL_SUCCESS != NormaliseBounds(tiling))
    XLAL_ERROR("NormaliseBounds failed", XLAL_EFAILED);
  
  /* Initialise bound checking */
  if (XLAL_SUCCESS != InitialiseBoundChecking(tiling))
    XLAL_ERROR("InitialiseBoundChecking failed", XLAL_EFAILED);

  /* If there are non-singular dimensions */
  if (r > 0) {

    int i, j;  
    gsl_matrix *reduced_metric = NULL;
    gsl_matrix *ortho_directions = NULL;
    gsl_matrix *sq_lwtri_generator = NULL;
    
    /* Allocate memory */
    ALLOC_GSL_VECTOR_INT(tiling->reduced_map,     r, XLAL_ERROR);
    ALLOC_GSL_MATRIX(tiling->latt_to_norm,     r, r, XLAL_ERROR);
    ALLOC_GSL_VECTOR_INT(tiling->latt_current,    r, XLAL_ERROR);
    ALLOC_GSL_VECTOR(tiling->norm_lower,          r, XLAL_ERROR);
    ALLOC_GSL_VECTOR(tiling->norm_upper,          r, XLAL_ERROR);
    ALLOC_GSL_MATRIX(reduced_metric,           r, r, XLAL_ERROR);
    ALLOC_GSL_MATRIX(ortho_directions,         r, r, XLAL_ERROR);

    /* Make map from reduced to full dimensions */
    {
      int ii;
      for (i = 0, ii = 0; i < r; ++i, ++ii) {
	while (gsl_vector_int_get(tiling->bound_type, ii) == FLT_PSBT_Singular)
	  ++ii;
	gsl_vector_int_set(tiling->reduced_map, i, ii);
      }
    }

    /* Extract reduced metric from normalised metric */
    for (i = 0; i < r; ++i) {
      const int ii = gsl_vector_int_get(tiling->reduced_map, i);
      for (j = 0; j < r; ++j) {
	const int jj = gsl_vector_int_get(tiling->reduced_map, j);
	gsl_matrix_set(reduced_metric, i, j, gsl_matrix_get(tiling->norm_metric, ii, jj));
      }
    }

    /* Use reduced metric ellipse bounding box as padding */
    if (NULL == (tiling->norm_padding = XLALMetricEllipseBoundingBox(reduced_metric, tiling->max_mismatch)))
      XLAL_ERROR("XLALMetricEllipseBoundingBox failed", XLAL_EFAILED);

    /* Find orthonormal directions of reduced metric */
    gsl_matrix_set_identity(ortho_directions);
    if (XLAL_SUCCESS != XLALOrthonormaliseWRTMetric(ortho_directions, reduced_metric))
      XLAL_ERROR("XLALOrthonormaliseWRTMetric failed", XLAL_EFAILED);

    /* Transform lattice generator to square lower triangular */
    if (NULL == (sq_lwtri_generator = XLALSquareLowerTriangularLatticeGenerator(generator)))
      XLAL_ERROR("XLALSquareLowerTriangularLatticeGenerator failed", XLAL_EFAILED);
    
    /* Normalise lattice generator so covering radius is sqrt(mismatch) */
    if (XLAL_SUCCESS != XLALNormaliseLatticeGenerator(sq_lwtri_generator, norm_thickness, sqrt(tiling->max_mismatch)))
      XLAL_ERROR("XLALNormaliseLatticeGenerator failed", XLAL_EFAILED);
    
    /* Compute transformation from lattice coordinates to parameter space */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, ortho_directions, sq_lwtri_generator, 0.0, tiling->latt_to_norm);
    
    /* Cleanup */
    FREE_GSL_MATRIX(reduced_metric);
    FREE_GSL_MATRIX(ortho_directions);
    FREE_GSL_MATRIX(sq_lwtri_generator);
    
  }

  /* Tiling is now fully initialised */
  tiling->state = FLT_S_NotStarted;

  return XLAL_SUCCESS;

}

/**
 * Set the tiling to a cubic lattice generator
 */
int XLALSetCubicTilingLattice(
			      FlatLatticeTiling *tiling /**< Tiling structure */
			      )
{

  const int r = tiling->reduced_dim;
  gsl_matrix *generator = NULL;
  REAL8 norm_thickness = 0.0;

  /* Check non-singular dimensions */
  if (r > 0) {

    /* Allocate memory */
    ALLOC_GSL_MATRIX(generator, r, r, XLAL_ERROR);
    
    /* Create generator */
    gsl_matrix_set_identity(generator);
    
    /* Calculate normalised thickness */
    norm_thickness = pow(sqrt(r)/2, r);

  }

  /* Set lattice generator */
  if (SetTilingLattice(tiling, generator, norm_thickness) != XLAL_SUCCESS)
    XLAL_ERROR("SetTilingLattice failed", XLAL_EFAILED);

  /* Cleanup */
  FREE_GSL_MATRIX(generator);

  return XLAL_SUCCESS;

}

/**
 * Set the tiling to an \fA_n^*\f lattice generator
 */
int XLALSetAnstarTilingLattice(
			       FlatLatticeTiling *tiling /**< Tiling structure */
			       )
{

  const int r = tiling->reduced_dim;
  gsl_matrix *generator = NULL;
  REAL8 norm_thickness = 0.0;

  /* Check non-singular dimensions */
  if (r > 0) {

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

  }

  /* Set lattice generator */
  if (SetTilingLattice(tiling, generator, norm_thickness) != XLAL_SUCCESS)
    XLAL_ERROR("SetTilingLattice failed", XLAL_EFAILED);

  /* Cleanup */
  FREE_GSL_MATRIX(generator);

  return XLAL_SUCCESS;

}

/**
 * Update the current point in normalised parameter space from 
 * the current point in lattice coordinates
 */
static void UpdateNormCurrent(
			      FlatLatticeTiling *tiling /**< Tiling structure */
			      )
{
  
  const int r = tiling->reduced_dim;

  int i, j;
  
  /* Compute the matrix transform for non-singular dimensions */
  gsl_vector_set_zero(tiling->norm_current);
  for (i = 0; i < r; ++i) {
    const int ii = gsl_vector_int_get(tiling->reduced_map, i);
    double x = gsl_vector_get(tiling->norm_current, ii);
    for (j = 0; j < r; ++j) {
      x += gsl_matrix_get(tiling->latt_to_norm, i, j) * 
	gsl_vector_int_get(tiling->latt_current, j);
    }
    gsl_vector_set(tiling->norm_current, ii, x);
  }

}

/**
 * Return to the first point in specified and higher dimensions
 */
static int ReturnToFirstPoint(
			      FlatLatticeTiling *tiling, /**< Tiling structure */
			      INT4 dimension             /**< Starting dimension */
			      )
{
  
  const int r = tiling->reduced_dim;

  int i;
  REAL8 lower = 0.0;
  REAL8 upper = 0.0;
  REAL8 dist = 0.0;
  int dlatt = 0;

  for (i = dimension; i < r; ++i) {
    
    /* Get mapped index and padding */
    const int ii = gsl_vector_int_get(tiling->reduced_map, i);
    const double norm_padding = gsl_vector_get(tiling->norm_padding, i);
	  
    /* Get bounds on this dimension */
    if (XLAL_SUCCESS != GetBounds(tiling, tiling->norm_current, ii, &lower, &upper))
      XLAL_ERROR("GetBounds failed", XLAL_EFAILED);
	  
    /* Store current bounds */
    gsl_vector_set(tiling->norm_lower, i, lower);
    gsl_vector_set(tiling->norm_upper, i, upper);
	  
    /* Move the current point past the lower bound */
    dist = lower - norm_padding - gsl_vector_get(tiling->norm_current, ii);
    dlatt = (int) ceil(dist / gsl_matrix_get(tiling->latt_to_norm, i, i));
    gsl_vector_int_set(tiling->latt_current, i, dlatt +
		       gsl_vector_int_get(tiling->latt_current, i));
    UpdateNormCurrent(tiling);
    
  }

  return XLAL_SUCCESS;

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
  BOOLEAN in_bounds = TRUE;
  BOOLEAN careful_linear = FALSE;
  BOOLEAN careful_quadratic = FALSE;
  
  /* If point is in bounds, metric ellipse
     definitely intersects parameter space */
  for (i = 0; i < r; ++i) {
    const int ii = gsl_vector_int_get(tiling->reduced_map, i);
    const double lower = gsl_vector_get(tiling->norm_lower, i);
    const double point = gsl_vector_get(tiling->norm_current, ii);
    const double upper = gsl_vector_get(tiling->norm_upper, i);
    if (point < lower || upper < point) {
      in_bounds = FALSE;
      careful_linear |= gsl_vector_int_get(tiling->careful_linear, ii);
      careful_quadratic |= gsl_vector_int_get(tiling->careful_quadratic, ii);
    }
  }
  if (in_bounds)
    return TRUE;

  /* Check the linear bounds:
     Solve the quadratic simplex problem to get the smallest
     distance from the current point to the parameter space. */
  if (careful_linear && tiling->linear_bound_A && tiling->linear_bound_b) {

    const int ln = (int)tiling->linear_map->size;

    double dist = 0.0;

    /* Copy current point onto linear bound space */
    for (i = 0; i < ln; ++i) {
      const int ii = gsl_vector_int_get(tiling->linear_map, i);
      gsl_vector_set(tiling->linear_current, i, gsl_vector_get(tiling->norm_current, ii));
    }

    /* Shift parameter space relative to current point */
    gsl_vector_memcpy(tiling->linear_transl_b, tiling->linear_bound_b);
    gsl_blas_dgemv(CblasNoTrans, -1.0, tiling->linear_bound_A, tiling->linear_current, 1.0, tiling->linear_transl_b);

    /* Update the quadratic problem simplex */
    if (XLAL_SUCCESS != XLALQuadraticSimplexMethodTableau(0, NULL, NULL, NULL, tiling->linear_transl_b, &tiling->linear_tableau))
      XLAL_ERROR("XLALQuadraticSimplexMethodTableau failed", XLAL_EFAILED);
    
    /* Solve to find the smallest distance */
    if (XLAL_SUCCESS != XLALQuadraticSimplexMethod(tiling->linear_tableau, tiling->linear_optimal))
      return FALSE;

    /* Calculate the distance with respect to the metric */
    gsl_blas_dgemv(CblasNoTrans, 1.0, tiling->linear_metric, tiling->linear_optimal, 0.0, tiling->linear_temp);
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

  int i;

  /* Switch on tiling state */
  switch (tiling->state) {

  case FLT_S_NotInitialised:

    /* Fail if uninitialised */
    XLAL_ERROR("'tiling' has not been fully initialised", XLAL_EINVAL);
    
  case FLT_S_Finished:

    /* Exit if finished */
    tiling->current_tile = NULL;
    return XLAL_FAILURE;
    
  case FLT_S_NotStarted:
    
    /* Tiling has started */
    tiling->tile_count = 0;
    tiling->state = FLT_S_InProgress;
    
    /* If all dimensions singular, return the single lattice point */
    if (r == 0) {
      gsl_vector_set_zero(tiling->norm_current);
      tiling->state = FLT_S_Finished;
    }

    /* Otherwise initialise current point */
    else {
      gsl_vector_int_set_zero(tiling->latt_current);
      UpdateNormCurrent(tiling);
      if (XLAL_SUCCESS != ReturnToFirstPoint(tiling, 0))
	XLAL_ERROR("ReturnToFirstPoint failed", XLAL_EFAILED);
      gsl_vector_int_set(tiling->latt_current, r - 1, 
			 gsl_vector_int_get(tiling->latt_current, r - 1)
			 - 1);
      UpdateNormCurrent(tiling);
    }

    /* Set pointer to current point */
    tiling->current_tile = tiling->norm_current;
    
  }
  
  /* Loop until found point in parameter space. This 
     will only happen at the parameter space boundaries. */
  if (tiling->state == FLT_S_InProgress)
    do {
      
      /* Loop over dimensions starting from highest */
      for (i = r - 1; i >= 0; --i) {
	
	/* Get mapped index and padding */
	const int ii = gsl_vector_int_get(tiling->reduced_map, i);
	const double norm_padding = gsl_vector_get(tiling->norm_padding, i);
	
	/* Increment lattice point along this dimension */
	gsl_vector_int_set(tiling->latt_current, i, 
			   gsl_vector_int_get(tiling->latt_current, i)
			   + 1);
	UpdateNormCurrent(tiling);
	
	/* Move point back to lower bound in higher dimensions */
	if (XLAL_SUCCESS != ReturnToFirstPoint(tiling, i + 1))
	  XLAL_ERROR("ReturnToFirstPoint failed", XLAL_EFAILED);
	
	/* If point is within bounds, stop moving point */
	if (gsl_vector_get(tiling->norm_current, ii) <=
	    gsl_vector_get(tiling->norm_upper, i) + norm_padding)
	  break;
	
	/* If this is lowest dimension, we're done! */
	if (i == 0) {
	  tiling->state = FLT_S_Finished;
	  return XLAL_FAILURE;
	}
	
      }
     
      /* Exit if metric ellipse intersects parameter space */
    } while (!DoesMetricEllipseIntersectParameterSpace(tiling));
  
  /* Template found, increment count */
  ++tiling->tile_count;
  
  /* Transform template point from normalised to real coordinates */
  gsl_vector_mul(tiling->norm_current, tiling->norm_to_real_mul);
  gsl_vector_add(tiling->norm_current, tiling->norm_to_real_add);

  return XLAL_SUCCESS;

}

/**
 * Return the count of the total number of flat lattice points
 */
REAL8 XLALTotalFlatLatticeTileCount(
				    FlatLatticeTiling *tiling /**< Tiling structure */
				    )
{
  
  /* Switch on tiling state */
  switch (tiling->state) {
    
  case FLT_S_NotInitialised:

    /* Fail if uninitialised */
    XLAL_ERROR("'tiling' has not been fully initialised", XLAL_EINVAL);
    
  case FLT_S_NotStarted:
    {
      
      int retn;

      /* Iterate through all templates */
      while ((retn = XLALNextFlatLatticeTile(tiling)) == XLAL_SUCCESS);
      if (retn != XLAL_FAILURE)
	XLAL_ERROR_REAL8("XLALNextFlatLatticeTile failed", XLAL_EFAILED);
      
      /* Reset tiling */
      tiling->state = FLT_S_NotStarted;
      
    }

  }
    
  /* Return the current/final number of templates */
  return (REAL8) tiling->tile_count;
  
}
