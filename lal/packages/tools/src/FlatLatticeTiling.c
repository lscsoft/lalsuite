//
// Copyright (C) 2007, 2008, 2012 Karl Wette
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA  02111-1307  USA
//

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

///
/// Get bounds on the specified dimension
///
static void GetBounds(
  FlatLatticeTiling* tiling,	///< [in] Tiling state
  size_t dimension,		///< [in] Dimension on which bound applies
  gsl_vector* point,		///< [in] Point at which to find bounds
  double* lower,		///< [out] Lower bound on point
  double* upper,		///< [out] Upper bound on point
  uint64_t *is_tiled		///< [out] Bit field of tiled dimensions
  );

///
/// Find the bounding box of the mismatch ellipses of a metric
//
static gsl_vector* MetricEllipseBoundingBox(
  gsl_matrix* metric,		///< [in] Metric to bound
  double max_mismatch		///< [in] Maximum mismatch with respect to metric
  );

///
/// Orthonormalise the columns of a matrix with respect to a metric (matrix is lower triangular)
///
static int OrthonormaliseWRTMetric(
  gsl_matrix* matrix,		///< [in] Matrix of columns to orthonormalise
  gsl_matrix* metric		///< [in] Metric to orthonormalise with respect to
  );

///
/// Transform a lattice generator to a square lower triangular form
///
static gsl_matrix* SquareLowerTriangularLatticeGenerator(
  gsl_matrix* generator		///< [in] Generator matrix of lattice
  );

///
/// Normalise a lattice generator matrix to have a specified covering radius
///
static int NormaliseLatticeGenerator(
  gsl_matrix* generator,	///< [in] Generator matrix of lattice
  double norm_thickness,	///< [in] Normalised thickness of lattice
  double covering_radius	///< [in] Desired covering radius
  );

///
/// Update the current flat lattice tiling subspace
///
static int UpdateSubspace(
  FlatLatticeTiling* tiling	///< [in] Tiling state
  );

///
/// Flat lattice tiling bound info
///
typedef struct tagFLT_Bound {
  FlatLatticeTilingBound func;		///< Parameter space bound function
  gsl_vector* data;			///< Arbitrary data describing parameter space
} FLT_Bound;

///
/// Flat lattice tiling subspace info
///
typedef struct tagFLT_Subspace {
  size_t dimensions;			///< Total number of tiled (non-singular) dimensions
  uint64_t is_tiled;			///< Dimensions which are tiled (non-singular)
  gsl_vector* padding;			///< Padding of bounds along each dimension
  gsl_matrix* increment;		///< Increment vectors of the lattice tiling generator
} FLT_Subspace;

///
/// Flat lattice tiling status
///
typedef enum tagFLT_Status {
  FLT_S_INITIAL,
  FLT_S_STARTED,
  FLT_S_FINISHED
} FLT_Status;

///
/// Flat lattice tiling state structure
///
struct tagFlatLatticeTiling {

  size_t dimensions;			///< Dimension of the parameter space
  FLT_Status status;			///< Status of the tiling

  size_t num_bounds;			///< Number of set parameter space bounds
  FLT_Bound *bounds;			///< Array of parameter space bound info
  gsl_vector* bound_point;		///< Temporary vector used in computing bounds
  double scale_padding;			///< Scaling of the padding of bounds (for testing)

  gsl_matrix* metric;			///< Parameter space metric in normalised coordinates
  gsl_vector* metric_scale;		///< Normalised to physical metric scaling

  double max_mismatch;			///< Maximum prescribed metric mismatch between templates
  FlatTilingLatticeGenerator generator;	///< Flat tiling lattice generator function

  size_t num_subspaces;			///< Number of cached tiling subspaces
  FLT_Subspace *subspaces;		///< Array of cached tiling subspaces
  FLT_Subspace *curr_subspace;		///< Current tiling subspace
  uint64_t curr_is_tiled;		///< Current dimensions which are tiled (non-singular)

  gsl_vector* curr_point;		///< Current lattice point
  gsl_vector* curr_lower;		///< Current lower bound on parameter space
  gsl_vector* curr_upper;		///< Current upper bound on parameter space

  gsl_vector* phys_scale;		///< Normalised to physical coordinate scaling
  gsl_vector* phys_offset;		///< Normalised to physical coordinate offset
  gsl_vector* curr_phys_point;		///< Current physical parameter-space point
  uint64_t count;			///< Total number of points generated so far

};

FlatLatticeTiling* XLALCreateFlatLatticeTiling(
  size_t dimensions
  )
{

  // Check input
  XLAL_CHECK_NULL(1 <= dimensions && dimensions <= 64, XLAL_EINVAL);

  const size_t n = dimensions;

  // Allocate and initialise tiling structure
  FlatLatticeTiling* tiling = XLALCalloc(1, sizeof(FlatLatticeTiling));
  XLAL_CHECK_NULL(tiling != NULL, XLAL_ENOMEM);
  tiling->dimensions = n;
  tiling->status = FLT_S_INITIAL;
  tiling->scale_padding = 1.0;
  tiling->bound_point = gsl_vector_alloc(n);
  XLAL_CHECK_NULL(tiling->bound_point != NULL, XLAL_ENOMEM);
  tiling->curr_point = gsl_vector_alloc(n);
  XLAL_CHECK_NULL(tiling->curr_point != NULL, XLAL_ENOMEM);
  tiling->curr_lower = gsl_vector_alloc(n);
  XLAL_CHECK_NULL(tiling->curr_lower != NULL, XLAL_ENOMEM);
  tiling->curr_upper = gsl_vector_alloc(n);
  XLAL_CHECK_NULL(tiling->curr_upper != NULL, XLAL_ENOMEM);
  tiling->phys_scale = gsl_vector_alloc(n);
  XLAL_CHECK_NULL(tiling->phys_scale != NULL, XLAL_ENOMEM);
  tiling->phys_offset = gsl_vector_alloc(n);
  XLAL_CHECK_NULL(tiling->phys_offset != NULL, XLAL_ENOMEM);
  tiling->curr_phys_point = gsl_vector_alloc(n);
  XLAL_CHECK_NULL(tiling->curr_phys_point != NULL, XLAL_ENOMEM);

  // Allocate and initialise bounds array
  tiling->bounds = XLALCalloc(n, sizeof(FLT_Bound));
  XLAL_CHECK_NULL(tiling->bounds != NULL, XLAL_ENOMEM);

  return tiling;

}

void XLALDestroyFlatLatticeTiling(
  FlatLatticeTiling* tiling
  )
{

  if (tiling) {

    // Destroy bounds array
    for (size_t k = 0; k < tiling->num_bounds; ++k) {
      gsl_vector_free(tiling->bounds[k].data);
    }
    XLALFree(tiling->bounds);

    // Destroy subspaces cache
    for (size_t k = 0; k < tiling->num_subspaces; ++k) {
      gsl_vector_free(tiling->subspaces[k].padding);
      gsl_matrix_free(tiling->subspaces[k].increment);
    }
    XLALFree(tiling->subspaces);

    // Destroy tiling structure
    gsl_vector_free(tiling->bound_point);
    gsl_matrix_free(tiling->metric);
    gsl_vector_free(tiling->metric_scale);
    gsl_vector_free(tiling->curr_point);
    gsl_vector_free(tiling->curr_lower);
    gsl_vector_free(tiling->curr_upper);
    gsl_vector_free(tiling->phys_scale);
    gsl_vector_free(tiling->phys_offset);
    gsl_vector_free(tiling->curr_phys_point);
    XLALFree(tiling);

  }

}

size_t XLALGetFlatLatticeTilingDimensions(
  FlatLatticeTiling* tiling
  )
{
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);
  return tiling->dimensions;
}

gsl_matrix* XLALGetFlatLatticeTilingMetric(
  FlatLatticeTiling* tiling
  )
{
  XLAL_CHECK_NULL(tiling != NULL, XLAL_EFAULT);
  return tiling->metric;
}

uint64_t XLALGetFlatLatticePointCount(
  FlatLatticeTiling* tiling
  )
{
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);
  return tiling->count;
}

int XLALAddFlatLatticeTilingBound(
  FlatLatticeTiling* tiling,
  FlatLatticeTilingBound func,
  gsl_vector* data
  )
{

  // Check tiling
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status == FLT_S_INITIAL, XLAL_EFAILED);

  // Check that there are dimensions still to be bounded
  XLAL_CHECK(tiling->num_bounds < tiling->dimensions, XLAL_EFAILED);

  // Set the next parameter space bound
  tiling->bounds[tiling->num_bounds].func = func;
  tiling->bounds[tiling->num_bounds].data = data;
  ++tiling->num_bounds;

  return XLAL_SUCCESS;

}

int XLALSetFlatLatticeTilingMetric(
  FlatLatticeTiling* tiling,
  gsl_matrix* metric,
  double max_mismatch
  )
{

  const size_t n = tiling->dimensions;

  // Check tiling
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status == FLT_S_INITIAL, XLAL_EFAILED);

  // Check that metric has not already been set
  XLAL_CHECK(tiling->metric == NULL, XLAL_EFAULT);

  // Check input
  XLAL_CHECK(metric != NULL, XLAL_EFAILED);
  XLAL_CHECK(metric->size1 == n && metric->size2 == n, XLAL_EINVAL);
  XLAL_CHECK(max_mismatch > 0, XLAL_EINVAL);

  // Allocate memory
  tiling->metric = gsl_matrix_alloc(n, n);
  XLAL_CHECK(tiling->metric != NULL, XLAL_ENOMEM);
  tiling->metric_scale = gsl_vector_alloc(n);
  XLAL_CHECK(tiling->metric_scale != NULL, XLAL_ENOMEM);

  // Check diagonal elements and calculate metric diagonal normalisation
  for (size_t i = 0; i < n; ++i) {
    const double metric_i_i = gsl_matrix_get(metric, i, i);
    XLAL_CHECK(metric_i_i > 0, XLAL_EINVAL);
    gsl_vector_set(tiling->metric_scale, i, 1.0 / sqrt(metric_i_i));
  }

  // Check metric is symmetric, and initialise metric with diagonal normalisation
  for (size_t i = 0; i < n; ++i) {
    const double scale_i = gsl_vector_get(tiling->metric_scale, i);
    gsl_matrix_set(tiling->metric, i, i, 1.0);
    for (size_t j = i+1; j < n; ++j) {
      const double scale_j = gsl_vector_get(tiling->metric_scale, j);
      double metric_i_j = gsl_matrix_get(metric, i, j);
      XLAL_CHECK(metric_i_j == gsl_matrix_get(metric, j, i), XLAL_EINVAL);
      metric_i_j *= scale_i * scale_j;
      gsl_matrix_set(tiling->metric, i, j, metric_i_j);
      gsl_matrix_set(tiling->metric, j, i, metric_i_j);
    }
  }

  // Initialise mismatch
  tiling->max_mismatch = max_mismatch;

  return XLAL_SUCCESS;

}

int XLALSetFlatTilingLatticeGenerator(
  FlatLatticeTiling* tiling,
  FlatTilingLatticeGenerator generator
  )
{

  // Check tiling
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status == FLT_S_INITIAL, XLAL_EFAILED);

  // Set the flat lattice tiling generator
  tiling->generator = generator;

  return XLAL_SUCCESS;

}

gsl_vector* XLALNextFlatLatticePoint(
  FlatLatticeTiling* tiling
  )
{

  const size_t n = tiling->dimensions;

  double lower, point, upper, padding, dist;

  // Check tiling
  XLAL_CHECK_NULL(tiling != NULL, XLAL_EFAULT);

  // If finished status, nothing to be done
  if (tiling->status == FLT_S_FINISHED) {
    return NULL;
  }

  // If initial status, perform final initialisation
  if (tiling->status == FLT_S_INITIAL) {

    // Check that all parameter space dimensions are bounded
    XLAL_CHECK_NULL(tiling->num_bounds == tiling->dimensions, XLAL_EFAILED);

    // Check that parameter space metric has been set
    XLAL_CHECK_NULL(tiling->metric != NULL, XLAL_EFAILED);

    // Check that lattice generator function has been set
    XLAL_CHECK_NULL(tiling->generator != NULL, XLAL_EFAILED);

    // Set physical parameter space offset
    gsl_vector_set_all(tiling->phys_scale, 1.0);
    gsl_vector_set_zero(tiling->phys_offset);
    for (size_t i = 0; i < n; ++i) {
      GetBounds(tiling, i, tiling->phys_offset, &lower, &upper, NULL);
      gsl_vector_set(tiling->phys_offset, i, lower);
    }

    // Set physical parameter space scaling
    gsl_vector_memcpy(tiling->phys_scale, tiling->metric_scale);

    // Initialise current point and bounds
    for (size_t i = 0; i < n; ++i) {

      // Get and initialise bounds
      GetBounds(tiling, i, tiling->curr_point, &lower, &upper, &tiling->curr_is_tiled);
      gsl_vector_set(tiling->curr_lower, i, lower);
      gsl_vector_set(tiling->curr_upper, i, upper);

      // Initialise point
      if (GET_BIT(uint64_t, tiling->curr_is_tiled, i)) {
        gsl_vector_set(tiling->curr_point, i, lower);
      }
      else {
        gsl_vector_set(tiling->curr_point, i, 0.5*(lower + upper));
      }

    }

    // Initialise subspace
    XLAL_CHECK_NULL(UpdateSubspace(tiling) == XLAL_SUCCESS, XLAL_EFAILED);

    // Subtract padding
    gsl_vector_sub(tiling->curr_point, tiling->curr_subspace->padding);

    // Initialise count
    tiling->count = 0;

    // Tiling has been started.
    tiling->status = FLT_S_STARTED;

  } else if (tiling->status == FLT_S_STARTED) {

    // Loop until a point is found
    ssize_t i = n;
    while (1) {

      // Decrease current dimension index
      --i;

      // If dimension index is less than zero, we're done!
      if (i < 0) {
        tiling->status = FLT_S_FINISHED;
        return NULL;
      }

      // If dimension is not tiled, move to lower dimension
      if (!GET_BIT(uint64_t, tiling->curr_is_tiled, i)) {
        continue;
      }

      // Increment current point along index
      {
        gsl_vector_view increment = gsl_matrix_column(tiling->curr_subspace->increment, i);
        gsl_vector_add(tiling->curr_point, &increment.vector);
      }

      // Get current point and bounds
      lower = gsl_vector_get(tiling->curr_lower, i);
      point = gsl_vector_get(tiling->curr_point, i);
      upper = gsl_vector_get(tiling->curr_upper, i);
      padding = gsl_vector_get(tiling->curr_subspace->padding, i);

      // If point is out of bounds, move to lower dimension
      if (point > upper + padding) {
        continue;
      }

      // Return point to lower bound in higher dimensions
      for (size_t j = i + 1; j < n; ++j) {

        // Get bounds
        GetBounds(tiling, j, tiling->curr_point, &lower, &upper, &tiling->curr_is_tiled);

        // Store bounds
        gsl_vector_set(tiling->curr_lower, j, lower);
        gsl_vector_set(tiling->curr_upper, j, upper);

        // Update subspace
        if (tiling->curr_is_tiled != tiling->curr_subspace->is_tiled) {
          XLAL_CHECK_NULL(UpdateSubspace(tiling) == XLAL_SUCCESS, XLAL_EFAILED);
        }

        // If dimension is tiled
        if (GET_BIT(uint64_t, tiling->curr_is_tiled, j)) {

          // Get increment vector
          gsl_vector_view increment = gsl_matrix_column(tiling->curr_subspace->increment, j);

          // Calculate the distance from current point to
          // the lower bound, in integer number of increments
          point = gsl_vector_get(tiling->curr_point, j);
          padding = gsl_vector_get(tiling->curr_subspace->padding, j);
          dist = floor((lower - padding - point) / gsl_vector_get(&increment.vector, j));

          // Move point back to lower bound
          gsl_blas_daxpy(dist, &increment.vector, tiling->curr_point);

        }

        // Otherwise just centre point
        else {
          gsl_vector_set(tiling->curr_point, j, 0.5*(lower + upper));
        }

      }

      // Found a template point
      break;

    }

  } else {
    XLAL_ERROR_NULL(XLAL_EFAILED, "Invalid tiling status!");
  }

  // Template was found, so increase count
  ++tiling->count;

  // Convert template to physical parameter space coordinates
  gsl_vector_memcpy(tiling->curr_phys_point, tiling->curr_point);
  gsl_vector_mul(tiling->curr_phys_point, tiling->phys_scale);
  gsl_vector_add(tiling->curr_phys_point, tiling->phys_offset);

  return tiling->curr_phys_point;

}

uint64_t XLALCountTotalFlatLatticePoints(
  FlatLatticeTiling* tiling
  )
{

  // Check tiling
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_VAL(0, tiling->status == FLT_S_INITIAL, XLAL_EFAILED);

  // Iterate over all templates
  int errnum;
  XLAL_TRY(while (XLALNextFlatLatticePoint(tiling) != NULL), errnum);
  XLAL_CHECK_VAL(0, errnum == 0, XLAL_EFAILED);

  // Save the template count
  uint64_t count = tiling->count;

  // Reset tiling
  tiling->status = FLT_S_INITIAL;
  tiling->count = 0;

  // Return the template count
  return count;

}

int XLALRandomPointInFlatLatticeParamSpace(
  FlatLatticeTiling* tiling,  ///< Tiling state
  RandomParams *randomParams, ///< Random parameters for generating random point
  gsl_vector* random_point,   ///< Random point
  gsl_vector* point,          ///< Another point
  double* metric_dist          ///< Distance from random point to other point w.r.t. metric
  )
{

  const size_t n = tiling->dimensions;

  double random_number;
  double lower, upper;
  gsl_vector* diff;

  // Create random point
  gsl_vector_set_zero(random_point);
  for (size_t i = 0; i < n; ++i) {

    // Get bounds
    GetBounds(tiling, i, random_point, &lower, &upper, NULL);

    // Generate random number
    random_number = XLALUniformDeviate(randomParams);

    // Generate random point
    gsl_vector_set(random_point, i, lower + random_number*(upper - lower));

  }

  // Convert random point to real parameter space coordinates
  gsl_vector_mul(random_point, tiling->phys_scale);
  gsl_vector_add(random_point, tiling->phys_offset);

  // Calculate distance from other point w.r.t metric
  if (point && metric_dist) {
    *metric_dist = 0.0;

    // Allocate memory
    diff = gsl_vector_alloc(n);
    XLAL_CHECK(diff != NULL, XLAL_ENOMEM);

    // Calculate difference between random and other point
    gsl_vector_memcpy(diff, point);
    gsl_vector_sub(diff, random_point);
    gsl_vector_div(diff, tiling->phys_scale);

    // Calculate off-diagonal parts (metric is symmetric) TODO USE GSL BLAS
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

    // Calculate diagonal components TODO USE GSL BLAS
    for (size_t i = 0; i < n; ++i) {
      if (gsl_vector_get(diff, i) != 0.0) {
        *metric_dist += gsl_matrix_get(tiling->metric, i, i) * gsl_vector_get(diff, i) * gsl_vector_get(diff, i);
      }
    }

    // Cleanup
    gsl_vector_free(diff);

  }

  return XLAL_SUCCESS;

}

int XLALCubicLatticeGenerator(
  size_t dimensions,
  gsl_matrix** generator,
  double* norm_thickness
  )
{

  const size_t r = dimensions;

  // Check input
  XLAL_CHECK(generator != NULL, XLAL_EFAULT);
  XLAL_CHECK(norm_thickness != NULL, XLAL_EFAULT);

  // Allocate memory
  *generator = gsl_matrix_alloc(r, r);
  XLAL_CHECK(*generator != NULL, XLAL_ENOMEM);

  // Create generator
  gsl_matrix_set_identity(*generator);

  // Calculate normalised thickness
  *norm_thickness = pow(sqrt(r)/2, r);

  return XLAL_SUCCESS;

}

int XLALAnstarLatticeGenerator(
  size_t dimensions,
  gsl_matrix** generator,
  double* norm_thickness
  )
{

  const size_t r = dimensions;

  // Check input
  XLAL_CHECK(generator != NULL, XLAL_EFAULT);
  XLAL_CHECK(norm_thickness != NULL, XLAL_EFAULT);

  // Allocate memory
  *generator = gsl_matrix_alloc(r + 1, r);
  XLAL_CHECK(*generator != NULL, XLAL_ENOMEM);

  // Create generator in (r + 1) space
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

  // Calculate normalised thickness
  *norm_thickness = sqrt(r + 1.0)*pow((1.0*r*(r + 2))/(12.0*(r + 1)), 0.5*r);

  return XLAL_SUCCESS;

}

static void ConstantBound(double* lower, double* upper, gsl_vector* point UNUSED, gsl_vector* data)
{

  // Set constant lower and upper bounds
  *lower = gsl_vector_get(data, 0);
  *upper = gsl_vector_get(data, 1);

}
int XLALAddFlatLatticeTilingConstantBound(
  FlatLatticeTiling* tiling,
  double lower,
  double upper
  )
{

  gsl_vector* data;

  // Check tiling
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);

  // Check input
  XLAL_CHECK(lower <= upper, XLAL_EINVAL);

  // Allocate and set bounds data
  data = gsl_vector_alloc(2);
  XLAL_CHECK(data != NULL, XLAL_ENOMEM);
  gsl_vector_set(data, 0, lower);
  gsl_vector_set(data, 1, upper);

  // Set parameter space
  XLAL_CHECK(XLALAddFlatLatticeTilingBound(tiling, ConstantBound, data) == XLAL_SUCCESS, XLAL_EFAILED);

  return XLAL_SUCCESS;

}

static void GetBounds(
  FlatLatticeTiling* tiling,
  size_t dimension,
  gsl_vector* point,
  double* lower,
  double* upper,
  uint64_t *is_tiled
  )
{

  // Get the appropriate bound dimension
  FLT_Bound* bound = &tiling->bounds[dimension];

  // Convert template to physical parameter space coordinates
  gsl_vector_memcpy(tiling->bound_point, point);
  gsl_vector_mul(tiling->bound_point, tiling->phys_scale);
  gsl_vector_add(tiling->bound_point, tiling->phys_offset);

  // Call parameter space bounds function
  if (dimension == 0) {
    (bound->func)(lower, upper, NULL, bound->data);
  }
  else {
    gsl_vector_view v = gsl_vector_subvector(tiling->bound_point, 0, dimension);
    (bound->func)(lower, upper, &v.vector, bound->data);
  }

  // Normalise bounds
  *lower -= gsl_vector_get(tiling->phys_offset, dimension);
  *upper -= gsl_vector_get(tiling->phys_offset, dimension);
  *lower /= gsl_vector_get(tiling->phys_scale, dimension);
  *upper /= gsl_vector_get(tiling->phys_scale, dimension);

  // Update whether dimension is singular
  if (is_tiled) {
    SET_BIT(uint64_t, *is_tiled, dimension, (*upper - *lower) > GSL_DBL_EPSILON);
  }

}

static gsl_vector* MetricEllipseBoundingBox(
  gsl_matrix* metric,
  double max_mismatch
  )
{

  const size_t n = metric->size1;

  // Check input
  XLAL_CHECK_NULL(metric->size1 == metric->size2, XLAL_ESIZE);

  // Allocate memory
  gsl_matrix* LU_decomp = gsl_matrix_alloc(n, n);
  XLAL_CHECK_NULL(LU_decomp != NULL, XLAL_ENOMEM);
  gsl_permutation* LU_perm = gsl_permutation_alloc(n);
  XLAL_CHECK_NULL(LU_perm != NULL, XLAL_ENOMEM);
  gsl_matrix* inverse = gsl_matrix_alloc(n, n);
  XLAL_CHECK_NULL(inverse != NULL, XLAL_ENOMEM);
  gsl_vector* bound_box = gsl_vector_alloc(n);
  XLAL_CHECK_NULL(bound_box != NULL, XLAL_ENOMEM);

  // Compute metric inverse
  int LU_sign = 0;
  gsl_matrix_memcpy(LU_decomp, metric);
  gsl_linalg_LU_decomp(LU_decomp, LU_perm, &LU_sign);
  gsl_linalg_LU_invert(LU_decomp, LU_perm, inverse);

  // Compute bounding box
  for (size_t i = 0; i < n; ++i) {
    gsl_vector_set(bound_box, i, sqrt(max_mismatch * gsl_matrix_get(inverse, i ,i)));
  }

  // Cleanup
  gsl_matrix_free(LU_decomp);
  gsl_permutation_free(LU_perm);
  gsl_matrix_free(inverse);

  return bound_box;

}

static int OrthonormaliseWRTMetric(
  gsl_matrix* matrix,
  gsl_matrix* metric
  )
{

  const size_t n = metric->size1;

  // Check input
  XLAL_CHECK(metric->size1 == metric->size2, XLAL_ESIZE);
  XLAL_CHECK(metric->size1 == matrix->size2 && metric->size2 == matrix->size2, XLAL_ESIZE);

  // Allocate
  gsl_vector* temp = gsl_vector_alloc(n);
  XLAL_CHECK(temp != NULL, XLAL_ENOMEM);

  // Orthonormalise the columns of the matrix using numerically stabilised Gram-Schmidt
  for (ssize_t i = n - 1; i >= 0; --i) {
    gsl_vector_view col_i = gsl_matrix_column(matrix, i);

    double inner_prod = 0.0;

    for (ssize_t j = n - 1; j > i; --j) {
      gsl_vector_view col_j = gsl_matrix_column(matrix, j);

      // Compute inner product of jth and ith columns with the metric
      gsl_blas_dgemv(CblasNoTrans, 1.0, metric, &col_j.vector, 0.0, temp);
      gsl_blas_ddot(&col_i.vector, temp, &inner_prod);

      // Subtract component of jth column from ith column
      gsl_vector_memcpy(temp, &col_j.vector);
      gsl_vector_scale(temp, inner_prod);
      gsl_vector_sub(&col_i.vector, temp);

    }

    // Compute inner product of ith column with itself
    gsl_blas_dgemv(CblasNoTrans, 1.0, metric, &col_i.vector, 0.0, temp);
    gsl_blas_ddot(&col_i.vector, temp, &inner_prod);

    // Normalise ith column
    gsl_vector_scale(&col_i.vector, 1.0 / sqrt(inner_prod));

  }

  // Cleanup
  gsl_vector_free(temp);

  return XLAL_SUCCESS;

}

static gsl_matrix* SquareLowerTriangularLatticeGenerator(
  gsl_matrix* generator
  )
{

  const size_t m = generator->size1;
  const size_t n = generator->size2;

  // Check input
  XLAL_CHECK_NULL(m >= n, XLAL_ESIZE);

  // Allocate memory
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

  // Find the QR decomposition of the generator
  gsl_matrix_memcpy(QR_decomp, generator);
  gsl_linalg_QR_decomp(QR_decomp, QR_tau);
  gsl_linalg_QR_unpack(QR_decomp, QR_tau, Q, R);

  // Build matrix to permute column order and make signs to diagonal positive
  gsl_matrix_set_zero(perm_sign);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      if (i + j == n - 1) {
        double x = gsl_matrix_get(R, j, j);
        gsl_matrix_set(perm_sign, i, j, x < 0 ? -1.0 : (x > 0 ? 1.0 : 0.0));
      }
    }
  }

  // Calculate left side of transform (Q is transposed to get inverse)
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, perm_sign, Q, 0.0, left);

  // Build right side of transform
  gsl_matrix_set_zero(right);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      if (i + j == n - 1) {
        gsl_matrix_set(right, i, j, 1.0);
      }
    }
  }

  // Transform generator
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, generator, right, 0.0, temp);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, left, temp, 0.0, result);

  // Generator will be lower triangular, so zero out upper triangle
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = i + 1; j < n; ++j) {
      gsl_matrix_set(result, i, j, 0.0);
    }
  }

  // Cleanup
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

static int NormaliseLatticeGenerator(
  gsl_matrix* generator,
  double norm_thickness,
  double covering_radius
  )
{

  const size_t n = generator->size1;

  gsl_matrix* LU_decomp = NULL;
  gsl_permutation *LU_perm = NULL;
  int LU_sign = 0;
  double generator_determinant = 0.0;
  double generator_covering_radius = 0.0;

  // Check input
  XLAL_CHECK(generator->size1 == generator->size2, XLAL_ESIZE);

  // Allocate memory
  LU_decomp = gsl_matrix_alloc(n, n);
  XLAL_CHECK(LU_decomp != NULL, XLAL_ENOMEM);
  LU_perm = gsl_permutation_alloc(n);
  XLAL_CHECK(LU_perm != NULL, XLAL_ENOMEM);

  // Compute generator LU decomposition
  gsl_matrix_memcpy(LU_decomp, generator);
  gsl_linalg_LU_decomp(LU_decomp, LU_perm, &LU_sign);

  // Compute generator determinant
  generator_determinant = gsl_linalg_LU_det(LU_decomp, LU_sign);

  // Compute covering radius
  generator_covering_radius = pow(norm_thickness * generator_determinant, 1.0 / n);

  // Normalise so covering spheres have specified covering radius
  gsl_matrix_scale(generator, covering_radius / generator_covering_radius);

  // Cleanup
  gsl_matrix_free(LU_decomp);
  gsl_permutation_free(LU_perm);

  return XLAL_SUCCESS;

}

static int UpdateSubspace(
  FlatLatticeTiling* tiling
  )
{

  const size_t n = tiling->dimensions;

  // Search for a previously-generated subspace
  for (size_t k = 0; k < tiling->num_subspaces; ++k) {
    tiling->curr_subspace = &tiling->subspaces[k];
    if (tiling->curr_subspace->is_tiled == tiling->curr_is_tiled) {
      return XLAL_SUCCESS;
    }
  }

  // (Re)Allocate memory for subspaces
  tiling->subspaces = XLALRealloc(tiling->subspaces, ++tiling->num_subspaces * sizeof(FLT_Subspace));
  XLAL_CHECK(tiling->subspaces != NULL, XLAL_ENOMEM);
  tiling->curr_subspace = &tiling->subspaces[tiling->num_subspaces - 1];

  // Initialise current subspace
  memset(tiling->curr_subspace, 0, sizeof(FLT_Subspace));
  tiling->curr_subspace->is_tiled = tiling->curr_is_tiled;

  // Count the number of tiled dimensions
  for (size_t i = 0; i < n; ++i) {
    if (GET_BIT(uint64_t, tiling->curr_subspace->is_tiled, i)) {
      ++tiling->curr_subspace->dimensions;
    }
  }

  // If subspace contains tiled dimensions
  if (tiling->curr_subspace->dimensions > 0) {

    const size_t r = tiling->curr_subspace->dimensions;

    // Allocate memory
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

    // Copy tiled dimensions of the metric and orthogonal directions
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

    // Use lengths of metric ellipse bounding box as padding along bounds
    gsl_vector* padding = MetricEllipseBoundingBox(metric, tiling->max_mismatch);
    XLAL_CHECK(padding != NULL, XLAL_EFAILED);
    gsl_vector_scale(padding, tiling->scale_padding);

    // Find orthonormalise directions with respect to subspace metric
    gsl_matrix_set_identity(orth_directions);
    XLAL_CHECK(OrthonormaliseWRTMetric(orth_directions, metric) == XLAL_SUCCESS, XLAL_EFAILED);

    // Get lattice generator
    gsl_matrix* generator = NULL;
    double norm_thickness = 0.0;
    XLAL_CHECK((tiling->generator)(r, &generator, &norm_thickness) == XLAL_SUCCESS, XLAL_EFAILED);

    // Transform lattice generator to square lower triangular
    gsl_matrix* sq_lwtri_generator = SquareLowerTriangularLatticeGenerator(generator);
    XLAL_CHECK(sq_lwtri_generator != NULL, XLAL_EFAILED);

    // Normalise lattice generator so covering radius is sqrt(mismatch)
    XLAL_CHECK(NormaliseLatticeGenerator(sq_lwtri_generator, norm_thickness, sqrt(tiling->max_mismatch)) == XLAL_SUCCESS, XLAL_EFAILED);

    // Compute the increment vectors of the lattice generator along the orthogonal directions
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, orth_directions, sq_lwtri_generator, 0.0, increment);

    // Copy the increment vectors so that non-tiled dimensions are zero
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

    // Cleanup
    gsl_matrix_free(metric);
    gsl_vector_free(padding);
    gsl_matrix_free(orth_directions);
    gsl_matrix_free(generator);
    gsl_matrix_free(sq_lwtri_generator);
    gsl_matrix_free(increment);

  }

  return XLAL_SUCCESS;

}
