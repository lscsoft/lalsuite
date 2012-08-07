/*
 *  Copyright (C) 2007, 2008, 2012 Karl Wette
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

#include <lal/LALStdlib.h>
#include <lal/LALMalloc.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/BitField.h>
#include <lal/FlatLatticeTiling.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/********** Structures and Enumerations **********/

/**
 * Flat lattice tiling bound
 */
typedef struct tagFlatLatticeTilingBound {

  /* Number of bound dimensions */
  size_t dimensions;

  /* Dimensions which are bound */
  uint64_t is_bound;

  /* Parameter space bound function */
  FlatLatticeTilingBoundFunc func;

  /* Arbitrary data describing parameter space */
  void* data;

  /* Cleanup function */
  FlatLatticeTilingBoundFree free;

} FlatLatticeTilingBound;

/**
 * Flat lattice tiling subspace
 */
typedef struct tagFlatLatticeTilingSubspace {

  /* Total number of tiled (non-flat) dimensions */
  size_t dimensions;

  /* Dimensions which are tiled (non-flat) */
  uint64_t is_tiled;

  /* Padding of bounds along each dimension */
  gsl_vector* padding;

  /* Increment vectors of the lattice tiling generator */
  gsl_matrix* increment;

} FlatLatticeTilingSubspace;

/**
 * State of the flat lattice tiling algorithm
 */
enum {
  FLT_S_NotInitialised,
  FLT_S_NotStarted,
  FLT_S_InProgress,
  FLT_S_Finished
};

/**
 * Flat lattice tiling state/input structure
 */
struct tagFlatLatticeTiling {

  /* Dimension of the parameter space */
  size_t dimensions;

  /* Parameter space bounds */
  size_t num_bounds;
  FlatLatticeTilingBound **bounds;
  gsl_vector_int *bound_map;
  gsl_vector* bound_point;

  /* Metric of the parameter space in normalised coordinates */
  gsl_matrix* metric;

  /* Normalised to real parameter coordinates scaling and offset */
  gsl_vector* real_scale;
  gsl_vector* real_offset;

  /* Maximum metric mismatch between the templates */
  double max_mismatch;

  /* Flat tiling lattice generator */
  FlatTilingLatticeGenerator generator;

  /* Cache of generated tiling subspaces */
  size_t num_subspaces;
  FlatLatticeTilingSubspace **subspaces;

  /* Scaling of the padding of bounds (for testing) */
  double scale_padding;

  /* Current dimensions which are tiled (non-flat) */
  uint64_t curr_is_tiled;

  /* Current tiling subspace */
  FlatLatticeTilingSubspace *curr_subspace;

  /* Current lattice point */
  gsl_vector* curr_point;

  /* Bounds on current point */
  gsl_vector* curr_lower;
  gsl_vector* curr_upper;

  /* Current template */
  gsl_vector* current;

  /* Total number of points generated so far */
  uint64_t count;

  /* State of the tiling */
  int state;

};

/********** Functions **********/

/**
 * Create a new flat lattice tiling bound structure
 */
static FlatLatticeTilingBound *CreateFlatLatticeTilingBound(void)
{

  /* Allocate memory */
  FlatLatticeTilingBound *bound = XLALMalloc(sizeof(FlatLatticeTilingBound));
  XLAL_CHECK_NULL(bound != NULL, XLAL_ENOMEM);

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

    if (bound->free) {
      (bound->free)(bound->data);
    }

    XLALFree(bound);

  }

}

/**
 * Create a new flat lattice tiling subspace structure
 */
static FlatLatticeTilingSubspace *CreateFlatLatticeTilingSubspace(void)
{

  /* Allocate memory */
  FlatLatticeTilingSubspace *subspace = XLALMalloc(sizeof(FlatLatticeTilingSubspace));
  XLAL_CHECK_NULL(subspace != NULL, XLAL_ENOMEM);

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

    gsl_vector_free(subspace->padding);
    gsl_matrix_free(subspace->increment);

    XLALFree(subspace);

  }

}

/**
 * Create a new flat lattice tiling structure
 */
FlatLatticeTiling* XLALCreateFlatLatticeTiling(
  size_t dimensions /**< Number of parameter space dimensions */
  )
{

  /* Check input */
  XLAL_CHECK_NULL(dimensions > 0, XLAL_EINVAL);

  /* Allocate memory */
  FlatLatticeTiling* tiling = XLALMalloc(sizeof(FlatLatticeTiling));
  XLAL_CHECK_NULL(tiling != NULL, XLAL_ENOMEM);

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

size_t XLALFlatLatticeTilingDimension(
  FlatLatticeTiling* tiling /**< Tiling structure */
  )
{
  return tiling->dimensions;
}

gsl_matrix* XLALFlatLatticeTilingMetric(
  FlatLatticeTiling* tiling /**< Tiling structure */
  )
{
  return tiling->metric;
}

/**
 * Free a flat lattice tiling structure
 */
void XLALDestroyFlatLatticeTiling(
  FlatLatticeTiling* tiling /**< Tiling structure */
  )
{

  if (tiling) {

    for (size_t k = 0; k < tiling->num_bounds; ++k) {
      FreeFlatLatticeTilingBound(tiling->bounds[k]);
    }
    XLALFree(tiling->bounds);
    gsl_vector_int_free(tiling->bound_map);
    gsl_vector_free(tiling->bound_point);
    gsl_matrix_free(tiling->metric);
    gsl_vector_free(tiling->real_scale);
    gsl_vector_free(tiling->real_offset);
    for (size_t k = 0; k < tiling->num_subspaces; ++k) {
      FreeFlatLatticeTilingSubspace(tiling->subspaces[k]);
    }
    XLALFree(tiling->subspaces);
    gsl_vector_free(tiling->curr_point);
    gsl_vector_free(tiling->curr_lower);
    gsl_vector_free(tiling->curr_upper);
    gsl_vector_free(tiling->current);

    XLALFree(tiling);

  }

}

/**
 * Add a parameter space for the flat lattice tiling
 */
int XLALAddFlatLatticeTilingBound(
  FlatLatticeTiling* tiling,             /**< Tiling structure */
  uint64_t bound_dimensions,                /**< Bit field indicating the dimensions bound */
  FlatLatticeTilingBoundFunc bound_func, /**< Parameter space bound function */
  void* bound_data,                      /**< Arbitrary data describing parameter space */
  FlatLatticeTilingBoundFree bound_free  /**< Cleanup function */
  )
{

  const size_t n = tiling->dimensions;

  FlatLatticeTilingBound *bound = NULL;

  /* Check tiling state */
  XLAL_CHECK(tiling->state == FLT_S_NotInitialised, XLAL_EFAILED);

  /* Check input */
  XLAL_CHECK(bound_dimensions != 0, XLAL_EINVAL);
  XLAL_CHECK(!(bound_dimensions & ~ALL_BITS(uint64_t, tiling->dimensions)), XLAL_EINVAL);
  for (size_t k = 0; k < tiling->num_bounds; ++k) {
    XLAL_CHECK(!(tiling->bounds[k]->is_bound & bound_dimensions), XLAL_EINVAL);
  }

  /* (Re)Allocate memory */
  if (!tiling->bound_map) {
    tiling->bound_map = gsl_vector_int_alloc(n);
    XLAL_CHECK(tiling->bound_map != NULL, XLAL_ENOMEM);
    gsl_vector_int_set_all(tiling->bound_map, -1);
  }
  if (!tiling->bound_point) {
    tiling->bound_point = gsl_vector_alloc(n);
    XLAL_CHECK(tiling->bound_point != NULL, XLAL_ENOMEM);
  }
  tiling->bounds = XLALRealloc(tiling->bounds, ++tiling->num_bounds * sizeof(FlatLatticeTilingBound*));
  XLAL_CHECK(tiling->bounds != NULL, XLAL_ENOMEM);
  bound = tiling->bounds[tiling->num_bounds - 1] = CreateFlatLatticeTilingBound();
  XLAL_CHECK(bound != NULL, XLAL_ENOMEM);

  /* Initialise structure */
  bound->dimensions = 0;
  bound->is_bound = bound_dimensions;
  bound->func = bound_func;
  bound->data = bound_data;
  bound->free = bound_free;

  /* Check bound map and count the number of bound dimensions */
  for (size_t i = 0; i < n; ++i) {
    if (GET_BIT(uint64_t, bound_dimensions, i)) {

      /* Check bound map */
      XLAL_CHECK(gsl_vector_int_get(tiling->bound_map, i) < 0, XLAL_EINVAL);

      /* Set bound map */
      gsl_vector_int_set(tiling->bound_map, i, tiling->num_bounds - 1);

      /* Increment number of bound dimensions */
      ++bound->dimensions;

    }
  }
  XLAL_CHECK(bound->dimensions == 1, XLAL_EINVAL);

  return XLAL_SUCCESS;

}

/**
 * Get bounds of the specified dimension
 */
static void GetBounds(
  FlatLatticeTiling* tiling, /**< Tiling structure */
  size_t dimension,            /**< Dimension on which bound applies */
  gsl_vector* point,         /**< Point on which to find bounds */
  double* lower,              /**< Lower bound on point in dimension */
  double* upper,              /**< Upper bound on point in dimension */
  uint64_t *is_tiled            /**< Bit field of tiled dimensions */
  )
{

  FlatLatticeTilingBound *bound = NULL;

  /* Get the appropriate bound dimension */
  bound = tiling->bounds[gsl_vector_int_get(tiling->bound_map, dimension)];

  /* Convert template to real parameter space coordinates */
  gsl_vector_memcpy(tiling->bound_point, point);
  gsl_vector_mul(tiling->bound_point, tiling->real_scale);
  gsl_vector_add(tiling->bound_point, tiling->real_offset);

  /* Call parameter space bounds function */
  if (dimension == 0) {
    (bound->func)(bound->data, dimension, NULL, lower, upper);
  }
  else {
    gsl_vector_view v = gsl_vector_subvector(tiling->bound_point, 0, dimension);
    (bound->func)(bound->data, dimension, &v.vector, lower, upper);
  }

  /* Normalise bounds */
  *lower -= gsl_vector_get(tiling->real_offset, dimension);
  *upper -= gsl_vector_get(tiling->real_offset, dimension);
  *lower /= gsl_vector_get(tiling->real_scale, dimension);
  *upper /= gsl_vector_get(tiling->real_scale, dimension);

  /* Update whether dimension is flat */
  if (is_tiled) {
    SET_BIT(uint64_t, *is_tiled, dimension, (*upper - *lower) > GSL_DBL_EPSILON);
  }

}

/**
 * Set the flat lattice tiling metric and maximum mismatch, and perform other initialisations
 */
int XLALSetFlatLatticeTilingMetric(
  FlatLatticeTiling* tiling, /**< Tiling structure */
  gsl_matrix* metric,        /**< Metric */
  double max_mismatch,        /**< Maximum mismatch */
  gsl_vector* real_scale     /**< Multiply to get real metric, may be NULL */
  )
{

  const size_t n = tiling->dimensions;

  double lower, upper;

  /* Check tiling state */
  XLAL_CHECK(tiling->bounds != NULL, XLAL_EINVAL);
  XLAL_CHECK(tiling->state == FLT_S_NotInitialised, XLAL_EFAILED);

  /* Check that all parameter space dimensions are bounded */
  XLAL_CHECK(tiling->bound_map != NULL, XLAL_EFAILED);
  for (size_t i = 0; i < n; ++i) {
    XLAL_CHECK(gsl_vector_int_get(tiling->bound_map, i) >=0, XLAL_EFAILED);
  }

  /* Check input */
  XLAL_CHECK(tiling->metric == NULL, XLAL_EFAILED);
  XLAL_CHECK(metric != NULL, XLAL_EFAILED);
  XLAL_CHECK(metric->size1 == metric->size2, XLAL_EINVAL);
  XLAL_CHECK(metric->size1 == n, XLAL_EINVAL);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = i+1; j < n; ++j) {
      XLAL_CHECK(gsl_matrix_get(metric, i, j) == gsl_matrix_get(metric, j, i), XLAL_EINVAL);
    }
  }
  XLAL_CHECK(max_mismatch > 0, XLAL_EINVAL);
  if (real_scale) {
    XLAL_CHECK(real_scale->size == n, XLAL_EINVAL);
    for (size_t i = 0; i < n; ++i) {
      XLAL_CHECK(gsl_vector_get(real_scale, i) > 0, XLAL_EINVAL);
    }
  }

  /* Allocate memory */
  tiling->metric = gsl_matrix_alloc(n, n);
  XLAL_CHECK(tiling->metric != NULL, XLAL_ENOMEM);
  tiling->real_scale = gsl_vector_alloc(n);
  XLAL_CHECK(tiling->real_scale != NULL, XLAL_ENOMEM);
  tiling->real_offset = gsl_vector_alloc(n);
  XLAL_CHECK(tiling->real_offset != NULL, XLAL_ENOMEM);
  tiling->curr_point = gsl_vector_alloc(n);
  XLAL_CHECK(tiling->curr_point != NULL, XLAL_ENOMEM);
  tiling->curr_lower = gsl_vector_alloc(n);
  XLAL_CHECK(tiling->curr_lower != NULL, XLAL_ENOMEM);
  tiling->curr_upper = gsl_vector_alloc(n);
  XLAL_CHECK(tiling->curr_upper != NULL, XLAL_ENOMEM);
  tiling->current = gsl_vector_alloc(n);
  XLAL_CHECK(tiling->current != NULL, XLAL_ENOMEM);

  /* Initialise normalised to real conversion */
  gsl_vector_set_all(tiling->real_scale, 1.0);
  gsl_vector_set_zero(tiling->real_offset);

  /* Find the real parameter space offset */
  for (size_t i = 0; i < n; ++i) {
    GetBounds(tiling, i, tiling->real_offset, &lower, &upper, NULL);
    gsl_vector_set(tiling->real_offset, i, lower);
  }

  /* Copy metric, normalised to real conversion, and mismatch */
  gsl_matrix_memcpy(tiling->metric, metric);
  if (real_scale) {
    gsl_vector_memcpy(tiling->real_scale, real_scale);
  }
  tiling->max_mismatch = max_mismatch;

  return XLAL_SUCCESS;

}

/**
 * Set the tiling lattice generator, and perform other initialisations
 */
int XLALSetFlatTilingLattice(
  FlatLatticeTiling* tiling,           /**< Tiling structure */
  FlatTilingLatticeGenerator generator /**< Flat lattice tiling generator */
  )
{

  /* Check tiling state */
  XLAL_CHECK(tiling->bounds != NULL, XLAL_EINVAL);
  XLAL_CHECK(tiling->state == FLT_S_NotInitialised, XLAL_EFAILED);

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
  FlatLatticeTiling* tiling /**< Tiling structure */
  )
{

  const size_t n = tiling->dimensions;

  /* Search for a previously-generated subspace */
  for (size_t k = 0; k < tiling->num_subspaces; ++k) {
    tiling->curr_subspace = tiling->subspaces[k];
    if (tiling->curr_subspace->is_tiled == tiling->curr_is_tiled) {
      return XLAL_SUCCESS;
    }
  }

  /* (Re)Allocate memory */
  tiling->subspaces = XLALRealloc(tiling->subspaces, ++tiling->num_subspaces * sizeof(FlatLatticeTilingSubspace*));
  XLAL_CHECK(tiling->subspaces != NULL, XLAL_ENOMEM);
  tiling->curr_subspace = tiling->subspaces[tiling->num_subspaces - 1] = CreateFlatLatticeTilingSubspace();
  XLAL_CHECK(tiling->curr_subspace != NULL, XLAL_EFAILED);

  /* Initialise structure */
  tiling->curr_subspace->dimensions = 0;
  tiling->curr_subspace->is_tiled = tiling->curr_is_tiled;
  tiling->curr_subspace->increment = NULL;

  /* Count the number of tiled dimensions */
  for (size_t i = 0; i < n; ++i) {
    if (GET_BIT(uint64_t, tiling->curr_subspace->is_tiled, i)) {
      ++tiling->curr_subspace->dimensions;
    }
  }
  size_t r = tiling->curr_subspace->dimensions;
  if (r > 0) {

    /* Allocate memory */
    gsl_matrix* metric = gsl_matrix_alloc(r, r);
    XLAL_CHECK(metric != NULL, XLAL_ENOMEM);
    gsl_matrix* orth_directions = gsl_matrix_alloc(r, r);
    XLAL_CHECK(orth_directions != NULL, XLAL_ENOMEM);
    gsl_matrix* increment = gsl_matrix_alloc(r, r);
    XLAL_CHECK(increment != NULL, XLAL_ENOMEM);
    tiling->curr_subspace->padding = gsl_vector_alloc(n);
    XLAL_CHECK(tiling->curr_subspace->padding != NULL, XLAL_ENOMEM);
    tiling->curr_subspace->increment = gsl_matrix_alloc(n, n);
    XLAL_CHECK(tiling->curr_subspace->increment != NULL, XLAL_ENOMEM);

    /* Copy tiled dimensions of the metric and orthogonal directions */
    for (size_t i = 0, ii = 0; i < n; ++i) {
      if (GET_BIT(uint64_t, tiling->curr_subspace->is_tiled, i)) {
        for (size_t j = 0, jj = 0; j < n; ++j) {
          if (GET_BIT(uint64_t, tiling->curr_subspace->is_tiled, j)) {
            gsl_matrix_set(metric, ii, jj, gsl_matrix_get(tiling->metric, i, j));
            ++jj;
          }
        }
        ++ii;
      }
    }

    /* Use lengths of metric ellipse bounding box as padding along bounds */
    gsl_vector* padding = XLALMetricEllipseBoundingBox(metric, tiling->max_mismatch);
    XLAL_CHECK(padding != NULL, XLAL_EFAILED);
    gsl_vector_scale(padding, tiling->scale_padding);

    /* Find orthonormalise directions with respect to subspace metric */
    gsl_matrix_set_identity(orth_directions);
    XLAL_CHECK(XLALOrthonormaliseWRTMetric(orth_directions, metric) == XLAL_SUCCESS, XLAL_EFAILED);

    /* Get lattice generator */
    gsl_matrix* generator = NULL;
    double norm_thickness = 0.0;
    XLAL_CHECK((tiling->generator)(r, &generator, &norm_thickness) == XLAL_SUCCESS, XLAL_EFAILED);

    /* Transform lattice generator to square lower triangular */
    gsl_matrix* sq_lwtri_generator = XLALSquareLowerTriangularLatticeGenerator(generator);
    XLAL_CHECK(sq_lwtri_generator != NULL, XLAL_EFAILED);

    /* Normalise lattice generator so covering radius is sqrt(mismatch) */
    XLAL_CHECK(XLALNormaliseLatticeGenerator(sq_lwtri_generator, norm_thickness, sqrt(tiling->max_mismatch)) == XLAL_SUCCESS, XLAL_EFAILED);

    /* Compute the increment vectors of the lattice generator along the orthogonal directions */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, orth_directions, sq_lwtri_generator, 0.0, increment);

    /* Copy the increment vectors so that non-tiled dimensions are zero */
    gsl_vector_set_zero(tiling->curr_subspace->padding);
    gsl_matrix_set_zero(tiling->curr_subspace->increment);
    for (size_t i = 0, ii = 0; i < n; ++i) {
      if (GET_BIT(uint64_t, tiling->curr_subspace->is_tiled, i)) {
        gsl_vector_set(tiling->curr_subspace->padding, i, gsl_vector_get(padding, ii));
        for (size_t j = 0, jj = 0; j < n; ++j) {
          if (GET_BIT(uint64_t, tiling->curr_subspace->is_tiled, j)) {
            gsl_matrix_set(tiling->curr_subspace->increment, i, j, gsl_matrix_get(increment, ii, jj));
            ++jj;
          }
        }
        ++ii;
      }
    }

    /* Cleanup */
    gsl_matrix_free(metric);
    gsl_vector_free(padding);
    gsl_matrix_free(orth_directions);
    gsl_matrix_free(generator);
    gsl_matrix_free(sq_lwtri_generator);
    gsl_matrix_free(increment);

  }

  return XLAL_SUCCESS;

}

/**
 * Move to the next point in the flat lattice tiling
 */
int XLALNextFlatLatticePoint(
  FlatLatticeTiling* tiling /**< Tiling structure */
  )
{

  const size_t n = tiling->dimensions;

  double lower, point, upper, padding, dist;

  /* Switch on tiling state */
  switch (tiling->state) {

  case FLT_S_NotInitialised:

    /* Fail if uninitialised */
    XLAL_ERROR(XLAL_EINVAL);

  case FLT_S_Finished:

    /* Exit if finished */
    return XLAL_FAILURE;

  case FLT_S_NotStarted:

    /* Initialise current point and bounds */
    for (size_t i = 0; i < n; ++i) {

      /* Get and initialise bounds */
      GetBounds(tiling, i, tiling->curr_point, &lower, &upper, &tiling->curr_is_tiled);
      gsl_vector_set(tiling->curr_lower, i, lower);
      gsl_vector_set(tiling->curr_upper, i, upper);

      /* Initialise point */
      if (GET_BIT(uint64_t, tiling->curr_is_tiled, i)) {
        gsl_vector_set(tiling->curr_point, i, lower);
      }
      else {
        gsl_vector_set(tiling->curr_point, i, 0.5*(lower + upper));
      }

    }

    /* Initialise subspace */
    XLAL_CHECK(UpdateFlatLatticeTilingSubspace(tiling) == XLAL_SUCCESS, XLAL_EFAILED);

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
    ssize_t i = n;
    while (1) {

      /* Decrease current dimension index */
      --i;

      /* If dimension index is less than zero, we're done! */
      if (i < 0) {
        return XLAL_FAILURE;
      }

      /* If dimension is not tiled, move to lower dimension */
      if (!GET_BIT(uint64_t, tiling->curr_is_tiled, i)) {
        continue;
      }

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
      if (point > upper + padding) {
        continue;
      }

      /* Return point to lower bound in higher dimensions */
      for (size_t j = i + 1; j < n; ++j) {

        /* Get bounds */
        GetBounds(tiling, j, tiling->curr_point, &lower, &upper, &tiling->curr_is_tiled);

        /* Store bounds */
        gsl_vector_set(tiling->curr_lower, j, lower);
        gsl_vector_set(tiling->curr_upper, j, upper);

        /* Update subspace */
        if (tiling->curr_is_tiled != tiling->curr_subspace->is_tiled) {
          XLAL_CHECK(UpdateFlatLatticeTilingSubspace(tiling) == XLAL_SUCCESS, XLAL_EFAILED);
        }

        /* If dimension is tiled */
        if (GET_BIT(uint64_t, tiling->curr_is_tiled, j)) {

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

gsl_vector*
XLALCurrentFlatLatticePoint(
  FlatLatticeTiling* tiling /**< Tiling structure */
  )
{
  return tiling->current;
}

/**
 * Return the count of the total number of flat lattice points
 */
uint64_t XLALTotalFlatLatticePointCount(
  FlatLatticeTiling* tiling /**< Tiling structure */
  )
{

  /* Switch on tiling state */
  switch (tiling->state) {

  case FLT_S_NotInitialised:

    /* Fail if uninitialised */
    XLAL_ERROR(0);

  case FLT_S_NotStarted:
  {

    int retn;

    /* Iterate through all templates */
    while ((retn = XLALNextFlatLatticePoint(tiling)) == XLAL_SUCCESS);
    if (retn != XLAL_FAILURE) {
      XLAL_ERROR(0);
    }

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
  FlatLatticeTiling* tiling,  /**< Tiling structure */
  RandomParams *randomParams, /**< Random parameters for generating random point */
  gsl_vector* random_point,   /**< Random point */
  gsl_vector* point,          /**< Another point */
  double* metric_dist          /**< Distance from random point to other point w.r.t. metric */
  )
{

  const size_t n = tiling->dimensions;

  double random_number;
  double lower, upper;
  gsl_vector* diff;

  /* Create random point */
  gsl_vector_set_zero(random_point);
  for (size_t i = 0; i < n; ++i) {

    /* Get bounds */
    GetBounds(tiling, i, random_point, &lower, &upper, NULL);

    /* Generate random number */
    random_number = XLALUniformDeviate(randomParams);

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
    diff = gsl_vector_alloc(n);
    XLAL_CHECK(diff != NULL, XLAL_ENOMEM);

    /* Calculate difference between random and other point */
    gsl_vector_memcpy(diff, point);
    gsl_vector_sub(diff, random_point);
    gsl_vector_div(diff, tiling->real_scale);

    /* Calculate off-diagonal parts (metric is symmetric) TODO USE GSL BLAS */
    for (size_t i = 0; i < n; ++i) {
      if (gsl_vector_get(diff, i) != 0.0) {
        for (size_t j = i + 1; j < n; ++j) {
          if (gsl_vector_get(diff, j) != 0.0) {
            *metric_dist += gsl_matrix_get(tiling->metric, i, j) * gsl_vector_get(diff, i) * gsl_vector_get(diff, j);
          }
        }
      }
    }
    *metric_dist *= 2.0;

    /* Calculate diagonal components TODO USE GSL BLAS */
    for (size_t i = 0; i < n; ++i) {
      if (gsl_vector_get(diff, i) != 0.0) {
        *metric_dist += gsl_matrix_get(tiling->metric, i, i) * gsl_vector_get(diff, i) * gsl_vector_get(diff, i);
      }
    }

    /* Cleanup */
    gsl_vector_free(diff);

  }

  return XLAL_SUCCESS;

}

/**
 * Find the bounding box of the mismatch ellipses of a metric
 */
gsl_vector* XLALMetricEllipseBoundingBox(
  gsl_matrix* metric, /**< Metric to bound */
  double max_mismatch  /**< Maximum mismatch w.r.t metric */
  )
{

  const size_t n = metric->size1;

  /* Check input */
  XLAL_CHECK_NULL(metric->size1 == metric->size2, XLAL_ESIZE);

  /* Allocate memory */
  gsl_matrix* LU_decomp = gsl_matrix_alloc(n, n);
  XLAL_CHECK_NULL(LU_decomp != NULL, XLAL_ENOMEM);
  gsl_permutation* LU_perm = gsl_permutation_alloc(n);
  XLAL_CHECK_NULL(LU_perm != NULL, XLAL_ENOMEM);
  gsl_matrix* inverse = gsl_matrix_alloc(n, n);
  XLAL_CHECK_NULL(inverse != NULL, XLAL_ENOMEM);
  gsl_vector* bound_box = gsl_vector_alloc(n);
  XLAL_CHECK_NULL(bound_box != NULL, XLAL_ENOMEM);

  /* Compute metric inverse */
  int LU_sign = 0;
  gsl_matrix_memcpy(LU_decomp, metric);
  gsl_linalg_LU_decomp(LU_decomp, LU_perm, &LU_sign);
  gsl_linalg_LU_invert(LU_decomp, LU_perm, inverse);

  /* Compute bounding box */
  for (size_t i = 0; i < n; ++i) {
    gsl_vector_set(bound_box, i, sqrt(max_mismatch * gsl_matrix_get(inverse, i ,i)));
  }

  /* Cleanup */
  gsl_matrix_free(LU_decomp);
  gsl_permutation_free(LU_perm);
  gsl_matrix_free(inverse);

  return bound_box;

}

/**
 * Orthonormalise the columns of a matrix with respect to a metric (matrix is lower triangular)
 */
int XLALOrthonormaliseWRTMetric(
  gsl_matrix* matrix, /**< Matrix of columns to orthonormalise */
  gsl_matrix* metric  /**< Metric to orthonormalise with respect to */
  )
{

  const size_t n = metric->size1;

  /* Check input */
  XLAL_CHECK(metric->size1 == metric->size2, XLAL_ESIZE);
  XLAL_CHECK(metric->size1 == matrix->size2 && metric->size2 == matrix->size2, XLAL_ESIZE);

  /* Allocate */
  gsl_vector* temp = gsl_vector_alloc(n);
  XLAL_CHECK(temp != NULL, XLAL_ENOMEM);

  /* Orthonormalise the columns of the matrix using numerically stabilised Gram-Schmidt */
  for (ssize_t i = n - 1; i >= 0; --i) {
    gsl_vector_view col_i = gsl_matrix_column(matrix, i);

    double inner_prod = 0.0;

    for (ssize_t j = n - 1; j > i; --j) {
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
  gsl_vector_free(temp);

  return XLAL_SUCCESS;

}

/**
 * Transform a lattice generator to a square lower triangular form
 */
gsl_matrix* XLALSquareLowerTriangularLatticeGenerator(
  gsl_matrix* generator /**< Generator matrix of lattice */
  )
{

  const size_t m = generator->size1;
  const size_t n = generator->size2;

  /* Check input */
  XLAL_CHECK_NULL(m >= n, XLAL_ESIZE);

  /* Allocate memory */
  gsl_matrix* QR_decomp = gsl_matrix_alloc(m, n);
  XLAL_CHECK_NULL(QR_decomp != NULL, XLAL_ENOMEM);
  gsl_vector* QR_tau = gsl_vector_alloc(n);
  XLAL_CHECK_NULL(QR_tau != NULL, XLAL_ENOMEM);
  gsl_matrix* Q = gsl_matrix_alloc(m, m);
  XLAL_CHECK_NULL(Q != NULL, XLAL_ENOMEM);
  gsl_matrix* R = gsl_matrix_alloc(m, n);
  XLAL_CHECK_NULL(R != NULL, XLAL_ENOMEM);
  gsl_matrix* perm_sign = gsl_matrix_alloc(n, m);
  XLAL_CHECK_NULL(perm_sign != NULL, XLAL_ENOMEM);
  gsl_matrix* left = gsl_matrix_alloc(n, m);
  XLAL_CHECK_NULL(left != NULL, XLAL_ENOMEM);
  gsl_matrix* right = gsl_matrix_alloc(n, n);
  XLAL_CHECK_NULL(right != NULL, XLAL_ENOMEM);
  gsl_matrix* temp = gsl_matrix_alloc(m, n);
  XLAL_CHECK_NULL(temp != NULL, XLAL_ENOMEM);
  gsl_matrix* result = gsl_matrix_alloc(n, n);
  XLAL_CHECK_NULL(result != NULL, XLAL_ENOMEM);

  /* Find the QR decomposition of the generator */
  gsl_matrix_memcpy(QR_decomp, generator);
  gsl_linalg_QR_decomp(QR_decomp, QR_tau);
  gsl_linalg_QR_unpack(QR_decomp, QR_tau, Q, R);

  /* Build matrix to permute column order and make signs to diagonal positive */
  gsl_matrix_set_zero(perm_sign);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      if (i + j == n - 1) {
        double x = gsl_matrix_get(R, j, j);
        gsl_matrix_set(perm_sign, i, j, x < 0 ? -1.0 : (x > 0 ? 1.0 : 0.0));
      }
    }
  }

  /* Calculate left side of transform (Q is transposed to get inverse) */
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, perm_sign, Q, 0.0, left);

  /* Build right side of transform */
  gsl_matrix_set_zero(right);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      if (i + j == n - 1) {
        gsl_matrix_set(right, i, j, 1.0);
      }
    }
  }

  /* Transform generator */
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, generator, right, 0.0, temp);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, left, temp, 0.0, result);

  /* Generator will be lower triangular, so zero out upper triangle */
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = i + 1; j < n; ++j) {
      gsl_matrix_set(result, i, j, 0.0);
    }
  }

  /* Cleanup */
  gsl_matrix_free(QR_decomp);
  gsl_vector_free(QR_tau);
  gsl_matrix_free(Q);
  gsl_matrix_free(R);
  gsl_matrix_free(perm_sign);
  gsl_matrix_free(left);
  gsl_matrix_free(right);
  gsl_matrix_free(temp);

  return result;

}

/**
 * Normalise a lattice generator matrix to have a specified covering radius
 */
int XLALNormaliseLatticeGenerator(
  gsl_matrix* generator, /**< Generator matrix of lattice */
  double norm_thickness,  /**< Normalised thickness of lattice */
  double covering_radius  /**< Desired covering radius */
  )
{

  const size_t n = generator->size1;

  gsl_matrix* LU_decomp = NULL;
  gsl_permutation *LU_perm = NULL;
  int LU_sign = 0;
  double generator_determinant = 0.0;
  double generator_covering_radius = 0.0;

  /* Check input */
  XLAL_CHECK(generator->size1 == generator->size2, XLAL_ESIZE);

  /* Allocate memory */
  LU_decomp = gsl_matrix_alloc(n, n);
  XLAL_CHECK(LU_decomp != NULL, XLAL_ENOMEM);
  LU_perm = gsl_permutation_alloc(n);
  XLAL_CHECK(LU_perm != NULL, XLAL_ENOMEM);

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
  gsl_matrix_free(LU_decomp);
  gsl_permutation_free(LU_perm);

  return XLAL_SUCCESS;

}

/**
 * Set the tiling to a cubic lattice generator
 */
static int FlatTilingCubicLatticeGenerator(size_t dimensions, gsl_matrix** generator, double* norm_thickness)
{

  const size_t r = dimensions;

  /* Allocate memory */
  *generator = gsl_matrix_alloc(r, r);
  XLAL_CHECK(*generator != NULL, XLAL_ENOMEM);

  /* Create generator */
  gsl_matrix_set_identity(*generator);

  /* Calculate normalised thickness */
  *norm_thickness = pow(sqrt(r)/2, r);

  return XLAL_SUCCESS;

}

int XLALSetFlatTilingCubicLattice(
  FlatLatticeTiling* tiling /**< Tiling structure */
  )
{

  XLAL_CHECK(XLALSetFlatTilingLattice(tiling, FlatTilingCubicLatticeGenerator) == XLAL_SUCCESS, XLAL_EFAILED);

  return XLAL_SUCCESS;

}

/**
 * Set the tiling to an \fA_n^*\f lattice generator
 */
static int FlatTilingAnstarLatticeGenerator(size_t dimensions, gsl_matrix** generator, double* norm_thickness)
{

  const size_t r = dimensions;

  /* Allocate memory */
  *generator = gsl_matrix_alloc(r + 1, r);
  XLAL_CHECK(*generator != NULL, XLAL_ENOMEM);

  /* Create generator in (r + 1) space */
  gsl_matrix_set_all(*generator, 0.0);
  {
    gsl_vector_view first_row = gsl_matrix_row(*generator, 0);
    gsl_vector_view sub_diag = gsl_matrix_subdiagonal(*generator, 1);
    gsl_vector_view last_col = gsl_matrix_column(*generator, r - 1);
    gsl_vector_set_all(&first_row.vector, 1.0);
    gsl_vector_set_all(&sub_diag.vector, -1.0);
    gsl_vector_set_all(&last_col.vector, 1.0 / (r + 1.0));
    gsl_vector_set(&last_col.vector, 0, -1.0 * r / (r + 1.0));
  }

  /* Calculate normalised thickness */
  *norm_thickness = sqrt(r + 1.0)*pow((1.0*r*(r + 2))/(12.0*(r + 1)), 0.5*r);

  return XLAL_SUCCESS;

}

int XLALSetFlatTilingAnstarLattice(
  FlatLatticeTiling* tiling /**< Tiling structure */
  )
{

  XLAL_CHECK(XLALSetFlatTilingLattice(tiling, FlatTilingAnstarLatticeGenerator) == XLAL_SUCCESS, XLAL_EFAILED);

  return XLAL_SUCCESS;

}

/**
 * Set a flat lattice tiling to a square parameter space
 */
static void ConstantBound(void* data, size_t dimension UNUSED, gsl_vector* point UNUSED, double* lower, double* upper)
{

  /* Set lower and upper bound */
  *lower = gsl_vector_get((gsl_vector*)data, 0);
  *upper = gsl_vector_get((gsl_vector*)data, 1);

}
static void ConstantFree(void* data)
{

  /* Cleanup */
  gsl_vector_free((gsl_vector*)data);

}
int XLALAddFlatLatticeTilingConstantBound(
  FlatLatticeTiling* tiling, /**< Tiling structure */
  size_t dimension,            /**< Dimension to bound */
  double lower,               /**< Lower bound on dimension */
  double upper                /**< Upper bound on dimension */
  )
{

  gsl_vector* data;

  /* Check input */
  XLAL_CHECK(dimension < tiling->dimensions, XLAL_EINVAL);
  XLAL_CHECK(lower <= upper, XLAL_EINVAL);

  /* Allocate memory */
  data = gsl_vector_alloc(2);
  XLAL_CHECK(data != NULL, XLAL_ENOMEM);

  /* Set bounds data */
  gsl_vector_set(data, 0, lower);
  gsl_vector_set(data, 1, upper);

  /* Set parameter space */
  XLAL_CHECK(XLALAddFlatLatticeTilingBound(tiling, ((uint64_t)(1)) << dimension, ConstantBound, (void*)data, ConstantFree)
             == XLAL_SUCCESS, XLAL_EFAILED);

  return XLAL_SUCCESS;

}
