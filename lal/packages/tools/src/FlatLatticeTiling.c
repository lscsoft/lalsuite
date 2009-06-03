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
#include <time.h>

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_sf.h>

#include <lal/LALRCSID.h>
#include <lal/LALStdlib.h>
#include <lal/LALMalloc.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/GSLSupport.h>
#include <lal/BitField.h>
#include <lal/FlatLatticeTiling.h>

NRCSID(FLATLATTICETILINGC, "$Id$");

#define TRUE  (1==1)
#define FALSE (1==0)

/**
 * Create a new flat lattice tiling bound structure
 */
static FlatLatticeTilingBound *CreateFlatLatticeTilingBound(void)
{

  FlatLatticeTilingBound *bound = NULL;

  /* Allocate memory */
  if ((bound = (FlatLatticeTilingBound*)XLALMalloc(sizeof(FlatLatticeTilingBound))) == NULL)
    XLAL_ERROR_NULL("Could not allocate 'bound'", XLAL_ENOMEM);

  /* Initialise structure */
  bound->dimensions = 0;
  bound->is_bound = 0;
  bound->func = NULL;
  bound->data = NULL;
  bound->free = NULL;

  return bound;

}

/**
 * Free a flat lattice tiling bound structure
 */
static void FreeFlatLatticeTilingBound(
				       FlatLatticeTilingBound *bound /**< Tiling bound structure */
				       )
{

  if (bound) {

    if (bound->free)
      (bound->free)(bound->data);

    XLALFree(bound);

  }

}

/**
 * Create a new flat lattice tiling subspace structure
 */
static FlatLatticeTilingSubspace *CreateFlatLatticeTilingSubspace(void)
{

  FlatLatticeTilingSubspace *subspace = NULL;

  /* Allocate memory */
  if ((subspace = (FlatLatticeTilingSubspace*)XLALMalloc(sizeof(FlatLatticeTilingSubspace))) == NULL)
    XLAL_ERROR_NULL("Could not allocate 'subspace'", XLAL_ENOMEM);

  /* Initialise structure */
  subspace->is_tiled = 0;
  subspace->dimensions = 0;
  subspace->padding = NULL;
  subspace->increment = NULL;

  return subspace;

}

/**
 * Free a flat lattice tiling subspace structure
 */
static void FreeFlatLatticeTilingSubspace(
					  FlatLatticeTilingSubspace *subspace /**< Tiling subspace structure */
					  )
{

  if (subspace) {

    FREE_GSL_VECTOR(subspace->padding);
    FREE_GSL_MATRIX(subspace->increment);

    XLALFree(subspace);

  }

}

/**
 * Create a new flat lattice tiling structure
 */
FlatLatticeTiling *XLALCreateFlatLatticeTiling(
					       INT4 dimensions /**< Number of parameter space dimensions */
					       )
{

  FlatLatticeTiling *tiling = NULL;

  /* Check input */
  if (dimensions <= 0)
    XLAL_ERROR_NULL("'dimensions' must be strictly positive", XLAL_EINVAL);

  /* Allocate memory */
  if ((tiling = (FlatLatticeTiling*)XLALMalloc(sizeof(FlatLatticeTiling))) == NULL)
    XLAL_ERROR_NULL("Could not allocate 'tiling'", XLAL_ENOMEM);

  /* Initialise structure */
  tiling->dimensions = dimensions;
  tiling->num_bounds = 0;
  tiling->bounds = NULL;
  tiling->bound_map = NULL;
  tiling->bound_point = NULL;
  tiling->metric = NULL;
  tiling->real_scale = NULL;
  tiling->real_offset = NULL;
  tiling->max_mismatch = 0.0;
  tiling->generator = NULL;
  tiling->num_subspaces = 0;
  tiling->subspaces = NULL;
  tiling->scale_padding = 1.0;
  tiling->curr_is_tiled = 0;
  tiling->curr_subspace = NULL;
  tiling->curr_point = NULL;
  tiling->curr_lower = NULL;
  tiling->curr_upper = NULL;
  tiling->current = NULL;
  tiling->count = 0;
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

  INT4 k;

  if (tiling) {

    for (k = 0; k < tiling->num_bounds; ++k)
      FreeFlatLatticeTilingBound(tiling->bounds[k]);
    XLALFree(tiling->bounds);
    FREE_GSL_VECTOR_INT(tiling->bound_map);
    FREE_GSL_VECTOR(tiling->bound_point);
    FREE_GSL_MATRIX(tiling->metric);
    FREE_GSL_VECTOR(tiling->real_scale);
    FREE_GSL_VECTOR(tiling->real_offset);
    for (k = 0; k < tiling->num_subspaces; ++k)
      FreeFlatLatticeTilingSubspace(tiling->subspaces[k]);
    XLALFree(tiling->subspaces);
    FREE_GSL_VECTOR(tiling->curr_point);
    FREE_GSL_VECTOR(tiling->curr_lower);
    FREE_GSL_VECTOR(tiling->curr_upper);
    FREE_GSL_VECTOR(tiling->current);

    XLALFree(tiling);

  }

}

/**
 * Add a parameter space for the flat lattice tiling
 */
int XLALAddFlatLatticeTilingBound(
				  FlatLatticeTiling *tiling,             /**< Tiling structure */
				  UINT8 bound_dimensions,                /**< Bit field indicating the dimensions bound */
				  FlatLatticeTilingBoundFunc bound_func, /**< Parameter space bound function */
				  void *bound_data,                      /**< Arbitrary data describing parameter space */
				  FlatLatticeTilingBoundFree bound_free  /**< Cleanup function */
				  )
{

  const int n = tiling->dimensions;

  int i, k;
  FlatLatticeTilingBound *bound = NULL;

  /* Check tiling state */
  if (tiling->state != FLT_S_NotInitialised)
    XLAL_ERROR("'tiling' has already been initialised", XLAL_EFAILED);

  /* Check input */
  if (bound_dimensions == 0)
    XLAL_ERROR("'bound_dimensions' must be non-zero", XLAL_EINVAL);
  if (bound_dimensions & ~ALL_BITS(UINT8, tiling->dimensions))
    XLAL_ERROR("'bound_dimensions' has bits set outside the dimensions of this tiling", XLAL_EINVAL);
  for (k = 0; k < tiling->num_bounds; ++k)
    if (tiling->bounds[k]->is_bound & bound_dimensions)
      XLAL_ERROR("'bound_dimensions' has bits set which conflict with an already defined bound of this tiling", XLAL_EINVAL);

  /* (Re)Allocate memory */
  if (!tiling->bound_map) {
    ALLOC_GSL_VECTOR_INT(tiling->bound_map, n, XLAL_ERROR);
    gsl_vector_int_set_all(tiling->bound_map, -1);
  }
  if (!tiling->bound_point)
    ALLOC_GSL_VECTOR(tiling->bound_point, n, XLAL_ERROR);
  if (NULL == (tiling->bounds = (FlatLatticeTilingBound**)XLALRealloc(tiling->bounds, ++tiling->num_bounds * sizeof(FlatLatticeTilingBound*))))
    XLAL_ERROR("Could not (re)allocate 'tiling->bounds'", XLAL_ENOMEM);
  if (NULL == (bound = (tiling->bounds[tiling->num_bounds - 1] = CreateFlatLatticeTilingBound())))
    XLAL_ERROR("CreateFlatLatticeTilingBound failed", XLAL_EFAILED);

  /* Initialise structure */
  bound->dimensions = 0;
  bound->is_bound = bound_dimensions;
  bound->func = bound_func;
  bound->data = bound_data;
  bound->free = bound_free;

  /* Check bound map and count the number of bound dimensions */
  for (i = 0; i < n; ++i)
    if (GET_BIT(UINT8, bound_dimensions, i)) {

      /* Check bound map */
      if (gsl_vector_int_get(tiling->bound_map, i) >= 0)
	XLAL_ERROR("'bound_dimensions' has bits set which conflict with previously defined bounds", XLAL_EINVAL);

      /* Set bound map */
      gsl_vector_int_set(tiling->bound_map, i, tiling->num_bounds - 1);

      /* Increment number of bound dimensions */
      ++bound->dimensions;

    }

  return XLAL_SUCCESS;

}

/**
 * Get bounds of the specified dimension
 */
static void GetBounds(
		      FlatLatticeTiling *tiling, /**< Tiling structure */
		      INT4 dimension,            /**< Dimension on which bound applies */
		      gsl_vector *point,         /**< Point on which to find bounds */
		      REAL8 *lower,              /**< Lower bound on point in dimension */
		      REAL8 *upper,              /**< Upper bound on point in dimension */
		      UINT8 *is_tiled            /**< Bit field of tiled dimensions */
		      )
{

  int i, ii;
  FlatLatticeTilingBound *bound = NULL;
  double x;
  BOOLEAN retn;

  /* Get the appropriate bound dimension */
  bound = tiling->bounds[gsl_vector_int_get(tiling->bound_map, dimension)];

  /* Copy the relevant values of point */
  for (i = 0, ii = 0; i < dimension; ++i) {
    if (GET_BIT(UINT8, bound->is_bound, i)) {
      x = gsl_vector_get(point, i);
      x *= gsl_vector_get(tiling->real_scale, i);
      x += gsl_vector_get(tiling->real_offset, i);
      gsl_vector_set(tiling->bound_point, ii, x);
      ++ii;
    }
  }

  /* Call parameter space bounds function */
  if (ii == 0) {
    retn = (bound->func)(bound->data, ii, NULL, lower, upper);
  }
  else {
    gsl_vector_view v = gsl_vector_subvector(tiling->bound_point, 0, ii);
    retn = (bound->func)(bound->data, ii, &v.vector, lower, upper);
  }

  /* Normalise bounds */
  *lower -= gsl_vector_get(tiling->real_offset, dimension);
  *upper -= gsl_vector_get(tiling->real_offset, dimension);
  *lower /= gsl_vector_get(tiling->real_scale, dimension);
  *upper /= gsl_vector_get(tiling->real_scale, dimension);

  /* Update whether dimension is flat */
  if (is_tiled)
    SET_BIT(UINT8, *is_tiled, dimension, (*upper - *lower) > GSL_DBL_EPSILON);

}

/**
 * Set the flat lattice tiling metric and maximum mismatch, and perform other initialisations
 */
int XLALSetFlatLatticeTilingMetric(
				   FlatLatticeTiling *tiling, /**< Tiling structure */
				   gsl_matrix *metric,        /**< Metric */
				   REAL8 max_mismatch,        /**< Maximum mismatch */
				   gsl_vector *real_scale     /**< Multiply to get real metric, may be NULL */
				   )
{

  const int n = tiling->dimensions;

  int i, j;
  REAL8 lower, upper;

  /* Check tiling state */
  if (!tiling->bounds)
    XLAL_ERROR("'tiling->bounds' has not been created", XLAL_EINVAL);
  if (tiling->state != FLT_S_NotInitialised)
    XLAL_ERROR("'tiling' has already been initialised", XLAL_EFAILED);

  /* Check that all parameter space dimensions are bounded */
  if (!tiling->bound_map)
    XLAL_ERROR("No parameter space bounds have been set on 'tiling'", XLAL_EFAILED);
  for (i = 0; i < n; ++i)
    if (gsl_vector_int_get(tiling->bound_map, i) < 0)
      XLAL_ERROR("Some parameter space dimensions have not been bounded", XLAL_EFAILED);

  /* Check input */
  if (tiling->metric)
    XLAL_ERROR("'tiling->metric' has already been set", XLAL_EFAILED);
  if (!metric)
    XLAL_ERROR("'metric' must be allocated", XLAL_EINVAL);
  if (metric->size1 != metric->size2)
    XLAL_ERROR("'metric' must be square", XLAL_EINVAL);
  if (n != (int)metric->size1)
    XLAL_ERROR("'metric' size must match tiling dimension", XLAL_EINVAL);
  for (i = 0; i < n; ++i)
    for (j = 0; j < n; ++j)
      if (gsl_matrix_get(metric, i, j) != gsl_matrix_get(metric, j, i))
	XLAL_ERROR("'metric' must be symmetric", XLAL_EINVAL);
  if (max_mismatch <= 0.0)
    XLAL_ERROR("'max_mismatch' must be strictly positive", XLAL_EINVAL);
  if (real_scale) {
    if (n != (int)real_scale->size)
      XLAL_ERROR("'real_scale' is not the correct size", XLAL_EINVAL);
    for (i = 0; i < n; ++i)
      if (gsl_vector_get(real_scale, i) <= 0.0)
	XLAL_ERROR("'real_scale' must be strictly positive", XLAL_EINVAL);
  }

  /* Allocate memory */
  ALLOC_GSL_MATRIX(tiling->metric,          n, n, XLAL_ERROR);
  ALLOC_GSL_VECTOR(tiling->real_scale,         n, XLAL_ERROR);
  ALLOC_GSL_VECTOR(tiling->real_offset,        n, XLAL_ERROR);
  ALLOC_GSL_VECTOR(tiling->curr_point,         n, XLAL_ERROR);
  ALLOC_GSL_VECTOR(tiling->curr_lower,         n, XLAL_ERROR);
  ALLOC_GSL_VECTOR(tiling->curr_upper,         n, XLAL_ERROR);
  ALLOC_GSL_VECTOR(tiling->current,            n, XLAL_ERROR);

  /* Initialise normalised to real conversion */
  gsl_vector_set_all(tiling->real_scale, 1.0);
  gsl_vector_set_zero(tiling->real_offset);

  /* Find the real parameter space offset */
  for (i = 0; i < n; ++i) {
    GetBounds(tiling, i, tiling->real_offset, &lower, &upper, NULL);
    gsl_vector_set(tiling->real_offset, i, lower);
  }

  /* Copy metric, normalised to real conversion, and mismatch */
  gsl_matrix_memcpy(tiling->metric, metric);
  if (real_scale)
    gsl_vector_memcpy(tiling->real_scale, real_scale);
  tiling->max_mismatch = max_mismatch;

  return XLAL_SUCCESS;

}

/**
 * Set the tiling lattice generator, and perform other initialisations
 */
int XLALSetFlatTilingLattice(
			     FlatLatticeTiling *tiling,           /**< Tiling structure */
			     FlatTilingLatticeGenerator generator /**< Flat lattice tiling generator */
			     )
{

  /* Check tiling state */
  if (!tiling->bounds)
    XLAL_ERROR("'tiling->bounds' has not been created", XLAL_EINVAL);
  if (tiling->state != FLT_S_NotInitialised)
    XLAL_ERROR("'tiling' has already been initialised", XLAL_EFAILED);

  /* Set the flat lattice tiling generator */
  tiling->generator = generator;

  /* Tiling is now fully initialised */
  tiling->state = FLT_S_NotStarted;

  return XLAL_SUCCESS;

}

/**
 * Update the current flat lattice tiling subspace
 */
static int UpdateFlatLatticeTilingSubspace(
					   FlatLatticeTiling *tiling /**< Tiling structure */
					   )
{

  const int n = tiling->dimensions;
  int r;

  int i, j, ii, jj, k;
  gsl_matrix *metric = NULL;
  gsl_vector *padding = NULL;
  gsl_matrix *orth_directions = NULL;
  gsl_matrix *generator = NULL;
  gsl_matrix *sq_lwtri_generator = NULL;
  gsl_matrix *increment = NULL;
  REAL8 norm_thickness;

  /* Search for a previously-generated subspace */
  for (k = 0; k < tiling->num_subspaces; ++k) {
    tiling->curr_subspace = tiling->subspaces[k];
    if (tiling->curr_subspace->is_tiled == tiling->curr_is_tiled)
      return XLAL_SUCCESS;
  }

  /* (Re)Allocate memory */
  if (NULL == (tiling->subspaces = (FlatLatticeTilingSubspace**)XLALRealloc(tiling->subspaces, ++tiling->num_subspaces * sizeof(FlatLatticeTilingSubspace*))))
    XLAL_ERROR("Could not (re)allocate 'tiling->subspaces'", XLAL_ENOMEM);
  if (NULL == (tiling->curr_subspace = (tiling->subspaces[tiling->num_subspaces - 1] = CreateFlatLatticeTilingSubspace())))
    XLAL_ERROR("CreateFlatLatticeTilingSubspace failed", XLAL_EFAILED);

  /* Initialise structure */
  tiling->curr_subspace->dimensions = 0;
  tiling->curr_subspace->is_tiled = tiling->curr_is_tiled;
  tiling->curr_subspace->increment = NULL;

  /* Count the number of tiled dimensions */
  for (i = 0; i < n; ++i)
    if (GET_BIT(UINT8, tiling->curr_subspace->is_tiled, i))
      ++tiling->curr_subspace->dimensions;
  if ((r = tiling->curr_subspace->dimensions) > 0) {

    /* Allocate memory */
    ALLOC_GSL_MATRIX(metric,                           r, r, XLAL_ERROR);
    ALLOC_GSL_MATRIX(orth_directions,                  r, r, XLAL_ERROR);
    ALLOC_GSL_MATRIX(increment,                        r, r, XLAL_ERROR);
    ALLOC_GSL_VECTOR(tiling->curr_subspace->padding,      n, XLAL_ERROR);
    ALLOC_GSL_MATRIX(tiling->curr_subspace->increment, n, n, XLAL_ERROR);

    /* Copy tiled dimensions of the metric and orthogonal directions */
    for (i = 0, ii = 0; i < n; ++i) {
      if (GET_BIT(UINT8, tiling->curr_subspace->is_tiled, i)) {
	for (j = 0, jj = 0; j < n; ++j) {
	  if (GET_BIT(UINT8, tiling->curr_subspace->is_tiled, j)) {
	    gsl_matrix_set(metric, ii, jj, gsl_matrix_get(tiling->metric, i, j));
	    ++jj;
	  }
	}
	++ii;
      }
    }

    /* Use lengths of metric ellipse bounding box as padding along bounds */
    if (NULL == (padding = XLALMetricEllipseBoundingBox(metric, tiling->max_mismatch)))
      XLAL_ERROR("XLALMetricEllipseBoundingBox failed", XLAL_EFAILED);
    gsl_vector_scale(padding, tiling->scale_padding);

    /* Find orthonormalise directions with respect to subspace metric */
    gsl_matrix_set_identity(orth_directions);
    if (XLAL_SUCCESS != XLALOrthonormaliseWRTMetric(orth_directions, metric))
      XLAL_ERROR("XLALOrthonormaliseWRTMetric failed", XLAL_EFAILED);

    /* Get lattice generator */
    if (XLAL_SUCCESS != (tiling->generator)(r, &generator, &norm_thickness))
      XLAL_ERROR("(tiling->generator) failed", XLAL_EFAILED);

    /* Transform lattice generator to square lower triangular */
    if (NULL == (sq_lwtri_generator = XLALSquareLowerTriangularLatticeGenerator(generator)))
      XLAL_ERROR("XLALSquareLowerTriangularLatticeGenerator failed", XLAL_EFAILED);

    /* Normalise lattice generator so covering radius is sqrt(mismatch) */
    if (XLAL_SUCCESS != XLALNormaliseLatticeGenerator(sq_lwtri_generator, norm_thickness, sqrt(tiling->max_mismatch)))
      XLAL_ERROR("XLALNormaliseLatticeGenerator failed", XLAL_EFAILED);

    /* Compute the increment vectors of the lattice generator along the orthogonal directions */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, orth_directions, sq_lwtri_generator, 0.0, increment);

    /* Copy the increment vectors so that non-tiled dimensions are zero */
    gsl_vector_set_zero(tiling->curr_subspace->padding);
    gsl_matrix_set_zero(tiling->curr_subspace->increment);
    for (i = 0, ii = 0; i < n; ++i) {
      if (GET_BIT(UINT8, tiling->curr_subspace->is_tiled, i)) {
	gsl_vector_set(tiling->curr_subspace->padding, i, gsl_vector_get(padding, ii));
	for (j = 0, jj = 0; j < n; ++j) {
	  if (GET_BIT(UINT8, tiling->curr_subspace->is_tiled, j)) {
	    gsl_matrix_set(tiling->curr_subspace->increment, i, j, gsl_matrix_get(increment, ii, jj));
	    ++jj;
	  }
	}
	++ii;
      }
    }

    /* Cleanup */
    FREE_GSL_MATRIX(metric);
    FREE_GSL_VECTOR(padding);
    FREE_GSL_MATRIX(orth_directions);
    FREE_GSL_MATRIX(generator);
    FREE_GSL_MATRIX(sq_lwtri_generator);
    FREE_GSL_MATRIX(increment);

  }

  return XLAL_SUCCESS;

}

/**
 * Move to the next point in the flat lattice tiling
 */
int XLALNextFlatLatticePoint(
			     FlatLatticeTiling *tiling /**< Tiling structure */
			     )
{

  const int n = tiling->dimensions;

  int i, j;
  REAL8 lower, point, upper, padding, dist;

  /* Switch on tiling state */
  switch (tiling->state) {

  case FLT_S_NotInitialised:

    /* Fail if uninitialised */
    XLAL_ERROR("'tiling' has not been fully initialised", XLAL_EINVAL);

  case FLT_S_Finished:

    /* Exit if finished */
    return XLAL_FAILURE;

  case FLT_S_NotStarted:

    /* Initialise current point and bounds */
    for (i = 0; i < n; ++i) {

      /* Get and initialise bounds */
      GetBounds(tiling, i, tiling->curr_point, &lower, &upper, &tiling->curr_is_tiled);
      gsl_vector_set(tiling->curr_lower, i, lower);
      gsl_vector_set(tiling->curr_upper, i, upper);

      /* Initialise point */
      if (GET_BIT(UINT8, tiling->curr_is_tiled, i)) {
	gsl_vector_set(tiling->curr_point, i, lower);
      }
      else {
	gsl_vector_set(tiling->curr_point, i, 0.5*(lower + upper));
      }

    }

    /* Initialise subspace */
    if (XLAL_SUCCESS != UpdateFlatLatticeTilingSubspace(tiling))
      XLAL_ERROR("UpdateFlatLatticeTilingSubspace failed", XLAL_EFAILED);

    /* Add padding */
    gsl_vector_sub(tiling->curr_point, tiling->curr_subspace->padding);

    /* Initialise count */
    tiling->count = 0;

    break;

  default:
    break;

  }

  /* Do not advance first point */
  if (tiling->state == FLT_S_NotStarted)
    tiling->state = FLT_S_InProgress;
  else {

    /* Loop until a point is found */
    i = n;
    while (TRUE) {

      /* Decrease current dimension index */
      --i;

      /* If dimension index is less than zero, we're done! */
      if (i < 0)
	return XLAL_FAILURE;

      /* If dimension is not tiled, move to lower dimension */
      if (!GET_BIT(UINT8, tiling->curr_is_tiled, i))
	continue;

      /* Increment current point along index */
      {
	gsl_vector_view increment = gsl_matrix_column(tiling->curr_subspace->increment, i);
	gsl_vector_add(tiling->curr_point, &increment.vector);
      }

      /* Get current point and bounds */
      lower = gsl_vector_get(tiling->curr_lower, i);
      point = gsl_vector_get(tiling->curr_point, i);
      upper = gsl_vector_get(tiling->curr_upper, i);
      padding = gsl_vector_get(tiling->curr_subspace->padding, i);

      /* If point is out of bounds, move to lower dimension */
      if (point > upper + padding)
	continue;

      /* Return point to lower bound in higher dimensions */
      for (j = i + 1; j < n; ++j) {

	/* Get bounds */
	GetBounds(tiling, j, tiling->curr_point, &lower, &upper, &tiling->curr_is_tiled);

	/* Store bounds */
	gsl_vector_set(tiling->curr_lower, j, lower);
	gsl_vector_set(tiling->curr_upper, j, upper);

	/* Update subspace */
	if (tiling->curr_is_tiled != tiling->curr_subspace->is_tiled)
	  if (XLAL_SUCCESS != UpdateFlatLatticeTilingSubspace(tiling))
	    XLAL_ERROR("UpdateFlatLatticeTilingSubspace failed", XLAL_EFAILED);

	/* If dimension is tiled */
	if (GET_BIT(UINT8, tiling->curr_is_tiled, j)) {

	  /* Get increment vector */
	  gsl_vector_view increment = gsl_matrix_column(tiling->curr_subspace->increment, j);

	  /* Calculate the distance from current point to
	     the lower bound, in integer number of increments */
	  point = gsl_vector_get(tiling->curr_point, j);
	  padding = gsl_vector_get(tiling->curr_subspace->padding, j);
	  dist = floor((lower - padding - point) / gsl_vector_get(&increment.vector, j));

	  /* Move point back to lower bound */
	  gsl_blas_daxpy(dist, &increment.vector, tiling->curr_point);

	}

	/* Otherwise just centre point */
	else {
	  gsl_vector_set(tiling->curr_point, j, 0.5*(lower + upper));
	}

      }

      /* Found a template point */
      break;

    }

  }

  /* Template was found, so increase count */
  ++tiling->count;

  /* Convert template to real parameter space coordinates */
  gsl_vector_memcpy(tiling->current, tiling->curr_point);
  gsl_vector_mul(tiling->current, tiling->real_scale);
  gsl_vector_add(tiling->current, tiling->real_offset);

  return XLAL_SUCCESS;

}

/**
 * Return the count of the total number of flat lattice points
 */
UINT4 XLALTotalFlatLatticePointCount(
				     FlatLatticeTiling *tiling /**< Tiling structure */
				     )
{

  /* Switch on tiling state */
  switch (tiling->state) {

  case FLT_S_NotInitialised:

    /* Fail if uninitialised */
    XLAL_ERROR("'tiling' has not been fully initialised", -1);

  case FLT_S_NotStarted:
    {

      int retn;

      /* Iterate through all templates */
      while ((retn = XLALNextFlatLatticePoint(tiling)) == XLAL_SUCCESS);
      if (retn != XLAL_FAILURE)
	XLAL_ERROR("XLALNextFlatLatticePoint failed", -1);

      /* Reset tiling */
      tiling->state = FLT_S_NotStarted;

    }
    break;

  default:
    break;

  }

  /* Return the current/final number of templates */
  return tiling->count;

}

int XLALRandomPointInFlatLatticeParamSpace(
					   FlatLatticeTiling *tiling, /**< Tiling structure */
					   RandomParams *random,      /**< Random parameters for generating random point */
					   gsl_vector* random_point,  /**< Random point */
					   gsl_vector* point,         /**< Another point */
					   REAL8* metric_dist         /**< Distance from random point to other point w.r.t. metric */
					   )
{

  const int n = tiling->dimensions;

  int i, j;
  REAL8 random_number;
  REAL8 lower, upper;
  gsl_vector *diff;

  /* Create random point */
  gsl_vector_set_zero(random_point);
  for (i = 0; i < n; ++i) {

    /* Get bounds */
    GetBounds(tiling, i, random_point, &lower, &upper, NULL);

    /* Generate random number */
    random_number = XLALUniformDeviate(random);

    /* Generate random point */
    gsl_vector_set(random_point, i, lower + random_number*(upper - lower));

  }

  /* Convert random point to real parameter space coordinates */
  gsl_vector_mul(random_point, tiling->real_scale);
  gsl_vector_add(random_point, tiling->real_offset);

  /* Calculate distance from other point w.r.t metric */
  if (point && metric_dist) {
    *metric_dist = 0.0;

    /* Allocate memory */
    ALLOC_GSL_VECTOR(diff, n, XLAL_ERROR);

    /* Calculate difference between random and other point */
    gsl_vector_memcpy(diff, point);
    gsl_vector_sub(diff, random_point);
    gsl_vector_div(diff, tiling->real_scale);

    /* Calculate off-diagonal parts (metric is symmetric) */
    for (i = 0; i < n; ++i)
      if (gsl_vector_get(diff, i) != 0.0)
	for (j = i + 1; j < n; ++j)
	  if (gsl_vector_get(diff, j) != 0.0)
	    *metric_dist += gsl_matrix_get(tiling->metric, i, j) *
	      gsl_vector_get(diff, i) * gsl_vector_get(diff, j);
    *metric_dist *= 2.0;

    /* Calculate diagonal components */
    for (i = 0; i < n; ++i)
      if (gsl_vector_get(diff, i) != 0.0)
	*metric_dist += gsl_matrix_get(tiling->metric, i, i) *
	  gsl_vector_get(diff, i) * gsl_vector_get(diff, i);

    /* Cleanup */
    FREE_GSL_VECTOR(diff);

  }

  return XLAL_SUCCESS;

}

/**
 * Find the principal of the mismatch ellipses of a metric
 */
gsl_matrix* XLALMetricEllipsePrincipalAxes(
					   gsl_matrix *metric, /**< Metric to bound */
					   REAL8 max_mismatch  /**< Maximum mismatch w.r.t metric */
					   )
{

  const int n = metric->size1;

  int i;
  gsl_eigen_symmv_workspace *eig_wksp = NULL;
  gsl_matrix *temp_matrix = NULL;
  gsl_vector *temp_vector = NULL;
  gsl_matrix *eig_vec = NULL;
  double inner_prod;

  /* Check input */
  if (n != (int)metric->size1 || n != (int)metric->size2)
    XLAL_ERROR_NULL("'metric' is not square", XLAL_ESIZE);

  /* Allocate memory */
  ALLOC_GSL_1D(eigen_symmv, eig_wksp,    n, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(temp_matrix,       n, n, XLAL_ERROR_NULL);
  ALLOC_GSL_VECTOR(temp_vector,          n, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(eig_vec,           n, n, XLAL_ERROR_NULL);

  /* Calculate the eigenvector of the metric */
  gsl_matrix_memcpy(temp_matrix, metric);
  gsl_eigen_symmv(temp_matrix, temp_vector, eig_vec, eig_wksp);

  /* Normalise the eigenvectors to the mismatch */
  for (i = 0; i < n; ++i) {

    /* Get view of eigenvector */
    gsl_vector_view eig_vec_i = gsl_matrix_column(eig_vec, i);

    /* Calculate inner product of eigenvector with respect to the metric */
    gsl_blas_dgemv(CblasNoTrans, 1.0, metric, &eig_vec_i.vector, 0.0, temp_vector);
    gsl_blas_ddot(&eig_vec_i.vector, temp_vector, &inner_prod);

    /* Normalise the eigenvector */
    gsl_vector_scale(&eig_vec_i.vector, sqrt(max_mismatch / inner_prod));

  }

  /* Cleanup */
  FREE_GSL_1D(eigen_symmv, eig_wksp);
  FREE_GSL_MATRIX(temp_matrix);
  FREE_GSL_VECTOR(temp_vector);

  return eig_vec;

}

/**
 * Find the bounding box of the mismatch ellipses of a metric
 */
gsl_vector *XLALMetricEllipseBoundingBox(
					 gsl_matrix *metric, /**< Metric to bound */
					 REAL8 max_mismatch  /**< Maximum mismatch w.r.t metric */
					 )
{

  const int n = metric->size1;

  int i;
  gsl_matrix *LU_decomp = NULL;
  gsl_permutation *LU_perm = NULL;
  gsl_matrix *inverse = NULL;
  int LU_sign = 0;
  gsl_vector *bound_box = NULL;

  /* Check input */
  if (n != (int)metric->size1 || n != (int)metric->size2)
    XLAL_ERROR_NULL("'metric' is not square", XLAL_ESIZE);

  /* Allocate memory */
  ALLOC_GSL_MATRIX(LU_decomp,    n, n, XLAL_ERROR_NULL);
  ALLOC_GSL_PERMUTATION(LU_perm,    n, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(inverse,      n, n, XLAL_ERROR_NULL);
  ALLOC_GSL_VECTOR(bound_box,       n, XLAL_ERROR_NULL);

  /* Compute metric inverse */
  gsl_matrix_memcpy(LU_decomp, metric);
  gsl_linalg_LU_decomp(LU_decomp, LU_perm, &LU_sign);
  gsl_linalg_LU_invert(LU_decomp, LU_perm, inverse);

  /* Compute bounding box */
  for (i = 0; i < n; ++i) {
    gsl_vector_set(bound_box, i, sqrt(max_mismatch * gsl_matrix_get(inverse, i ,i)));
  }

  /* Cleanup */
  FREE_GSL_MATRIX(LU_decomp);
  FREE_GSL_PERMUTATION(LU_perm);
  FREE_GSL_MATRIX(inverse);

  return bound_box;

}

/**
 * Orthonormalise the columns of a matrix with respect to a metric (matrix is lower triangular)
 */
int XLALOrthonormaliseWRTMetric(
				gsl_matrix *matrix, /**< Matrix of columns to orthonormalise */
				gsl_matrix *metric  /**< Metric to orthonormalise with respect to */
				)
{

  const int n = metric->size1;

  int i, j;
  gsl_vector *temp = NULL;
  double inner_prod = 0.0;

  /* Check input */
  if (n != (int)metric->size1 || n != (int)metric->size2)
    XLAL_ERROR("'metric' is not square", XLAL_ESIZE);
  if (metric->size1 != matrix->size2 || metric->size2 != matrix->size2)
    XLAL_ERROR("'matrix' is not the same size as 'metric'", XLAL_ESIZE);

  /* Allocate */
  ALLOC_GSL_VECTOR(temp, n, XLAL_ERROR);

  /* Orthonormalise the columns of the matrix using numerically stabilised Gram-Schmidt */
  for (i = n - 1; i >= 0; --i) {
    gsl_vector_view col_i = gsl_matrix_column(matrix, i);
    for (j = n - 1; j > i; --j) {
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

  /* Cleanup */
  FREE_GSL_VECTOR(temp);

  return XLAL_SUCCESS;

}

/**
 * Transform a lattice generator to a square lower triangular form
 */
gsl_matrix *XLALSquareLowerTriangularLatticeGenerator(
						      gsl_matrix *generator /**< Generator matrix of lattice */
						      )
{

  const int m = generator->size1;
  const int n = generator->size2;

  int i, j;
  gsl_matrix *QR_decomp = NULL;
  gsl_vector *QR_tau = NULL;
  gsl_matrix *Q = NULL;
  gsl_matrix *R = NULL;
  gsl_matrix *perm_sign = NULL;
  gsl_matrix *left = NULL;
  gsl_matrix *right = NULL;
  gsl_matrix *temp = NULL;
  gsl_matrix *result = NULL;
  double x = 0.0;

  /* Check input */
  if (m < n)
    XLAL_ERROR_NULL("'generator' must have number of rows >= number of columns", XLAL_ESIZE);

  /* Allocate memory */
  ALLOC_GSL_MATRIX(QR_decomp, m, n, XLAL_ERROR_NULL);
  ALLOC_GSL_VECTOR(QR_tau,       n, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(Q,         m, m, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(R,         m, n, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(perm_sign, n, m, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(left,      n, m, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(right,     n, n, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(temp,      m, n, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(result,    n, n, XLAL_ERROR_NULL);

  /* Find the QR decomposition of the generator */
  gsl_matrix_memcpy(QR_decomp, generator);
  gsl_linalg_QR_decomp(QR_decomp, QR_tau);
  gsl_linalg_QR_unpack(QR_decomp, QR_tau, Q, R);

  /* Build matrix to permute column order and make signs to diagonal positive */
  gsl_matrix_set_zero(perm_sign);
  for (i = 0; i < n; ++i) {
    for (j = 0; j < m; ++j) {
      if (i + j == n - 1) {
	x = gsl_matrix_get(R, j, j);
	gsl_matrix_set(perm_sign, i, j, x < 0 ? -1.0 : (x > 0 ? 1.0 : 0.0));
      }
    }
  }

  /* Calculate left side of transform (Q is transposed to get inverse) */
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, perm_sign, Q, 0.0, left);

  /* Build right side of transform */
  gsl_matrix_set_zero(right);
  for (i = 0; i < n; ++i)
    for (j = 0; j < n; ++j)
      if (i + j == n - 1)
	gsl_matrix_set(right, i, j, 1.0);

  /* Transform generator */
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, generator, right, 0.0, temp);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, left, temp, 0.0, result);

  /* Generator will be lower triangular, so zero out upper triangle */
  for (i = 0; i < n; ++i)
    for (j = i + 1; j < n; ++j)
      gsl_matrix_set(result, i, j, 0.0);

  /* Cleanup */
  FREE_GSL_MATRIX(QR_decomp);
  FREE_GSL_VECTOR(QR_tau);
  FREE_GSL_MATRIX(Q);
  FREE_GSL_MATRIX(R);
  FREE_GSL_MATRIX(perm_sign);
  FREE_GSL_MATRIX(left);
  FREE_GSL_MATRIX(right);
  FREE_GSL_MATRIX(temp);

  return result;

}

/**
 * Normalise a lattice generator matrix to have a specified covering radius
 */
int XLALNormaliseLatticeGenerator(
				  gsl_matrix *generator, /**< Generator matrix of lattice */
				  REAL8 norm_thickness,  /**< Normalised thickness of lattice */
				  REAL8 covering_radius  /**< Desired covering radius */
				  )
{

  const int n = generator->size1;

  gsl_matrix *LU_decomp = NULL;
  gsl_permutation *LU_perm = NULL;
  int LU_sign = 0;
  double generator_determinant = 0.0;
  double generator_covering_radius = 0.0;

  /* Check input */
  if (n != (int)generator->size1 || n != (int)generator->size2)
    XLAL_ERROR("'generator' is not square", XLAL_ESIZE);

  /* Allocate memory */
  ALLOC_GSL_MATRIX(LU_decomp, n, n, XLAL_ERROR);
  ALLOC_GSL_PERMUTATION(LU_perm, n, XLAL_ERROR);

  /* Compute generator LU decomposition */
  gsl_matrix_memcpy(LU_decomp, generator);
  gsl_linalg_LU_decomp(LU_decomp, LU_perm, &LU_sign);

  /* Compute generator determinant */
  generator_determinant = gsl_linalg_LU_det(LU_decomp, LU_sign);

  /* Compute covering radius */
  generator_covering_radius = pow(norm_thickness * generator_determinant, 1.0 / n);

  /* Normalise so covering spheres have specified covering radius */
  gsl_matrix_scale(generator, covering_radius / generator_covering_radius);

  /* Cleanup */
  FREE_GSL_MATRIX(LU_decomp);
  FREE_GSL_PERMUTATION(LU_perm);

  return XLAL_SUCCESS;

}

/**
 * Set the tiling to a cubic lattice generator
 */
static int FlatTilingCubicLatticeGenerator(INT4 dimensions, gsl_matrix** generator, REAL8* norm_thickness)
{

  const int r = dimensions;

  /* Allocate memory */
  ALLOC_GSL_MATRIX(*generator, r, r, XLAL_ERROR);

  /* Create generator */
  gsl_matrix_set_identity(*generator);

  /* Calculate normalised thickness */
  *norm_thickness = pow(sqrt(r)/2, r);

  return XLAL_SUCCESS;

}

int XLALSetFlatTilingCubicLattice(
				  FlatLatticeTiling *tiling /**< Tiling structure */
				  )
{

  if (XLAL_SUCCESS != XLALSetFlatTilingLattice(tiling, FlatTilingCubicLatticeGenerator))
    XLAL_ERROR("XLALSetFlatTilingLattice failed", XLAL_EFAILED);

  return XLAL_SUCCESS;

}

/**
 * Set the tiling to an \fA_n^*\f lattice generator
 */
static int FlatTilingAnstarLatticeGenerator(INT4 dimensions, gsl_matrix** generator, REAL8* norm_thickness)
{

  const int r = dimensions;

  /* Allocate memory */
  ALLOC_GSL_MATRIX(*generator, r + 1, r, XLAL_ERROR);

  /* Create generator in (r + 1) space */
  gsl_matrix_set_all(*generator, 0.0);
  {
    gsl_vector_view first_row = gsl_matrix_row(*generator, 0);
    gsl_vector_view sub_diag = gsl_matrix_subdiagonal(*generator, 1);
    gsl_vector_view last_col = gsl_matrix_column(*generator, r - 1);
    gsl_vector_set_all(&first_row.vector, 1.0);
    gsl_vector_set_all(&sub_diag.vector, -1.0);
    gsl_vector_set_all(&last_col.vector, 1.0 / (r + 1.0));
    gsl_vector_set(&last_col.vector, 0, -r * 1.0 / (r + 1.0));
  }

  /* Calculate normalised thickness */
  *norm_thickness = sqrt(r + 1.0)*pow((1.0*r*(r + 2))/(12.0*(r + 1)), 0.5*r);

  return XLAL_SUCCESS;

}

int XLALSetFlatTilingAnstarLattice(
				   FlatLatticeTiling *tiling /**< Tiling structure */
				   )
{

  if (XLAL_SUCCESS != XLALSetFlatTilingLattice(tiling, FlatTilingAnstarLatticeGenerator))
    XLAL_ERROR("XLALSetFlatTilingLattice failed", XLAL_EFAILED);

  return XLAL_SUCCESS;

}

/**
 * Set a flat lattice tiling to a square parameter space
 */
static BOOLEAN ConstantBound(void *data, INT4 dimension, gsl_vector *point, REAL8 *lower, REAL8 *upper)
{

  /* Set lower and upper bound */
  *lower = gsl_vector_get((gsl_vector*)data, 0);
  *upper = gsl_vector_get((gsl_vector*)data, 1);

  return TRUE;

}
static void ConstantFree(void *data)
{

  /* Cleanup */
  FREE_GSL_VECTOR((gsl_vector*)data);

}
int XLALAddFlatLatticeTilingConstantBound(
					  FlatLatticeTiling *tiling, /**< Tiling structure */
					  INT4 dimension,            /**< Dimension to bound */
					  REAL8 lower,               /**< Lower bound on dimension */
					  REAL8 upper                /**< Upper bound on dimension */
					  )
{

  gsl_vector *data;

  /* Check input */
  if (dimension < 0 || tiling->dimensions <= dimension)
    XLAL_ERROR("'dimension' is out of bounds", XLAL_EINVAL);
  if (lower > upper)
    XLAL_ERROR("'lower' must be less than or equal to 'upper'", XLAL_EINVAL);

  /* Allocate memory */
  ALLOC_GSL_VECTOR(data, 2, XLAL_ERROR);

  /* Set bounds data */
  gsl_vector_set(data, 0, lower);
  gsl_vector_set(data, 1, upper);

  /* Set parameter space */
  if (XLAL_SUCCESS != XLALAddFlatLatticeTilingBound(tiling, ((UINT8)(1)) << dimension,
						    ConstantBound, (void*)data, ConstantFree))
    XLAL_ERROR("XLALAddFlatLatticeTilingBound failed", XLAL_EFAILED);

  return XLAL_SUCCESS;

}
