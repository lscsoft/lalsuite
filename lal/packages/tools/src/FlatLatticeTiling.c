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
#include <lal/FlatLatticeTiling.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

///
/// Get physical bounds on the specified dimension
///
static void GetPhysBounds(
  FlatLatticeTiling* tiling,		///< [in] Tiling state
  const size_t dimension,		///< [in] Dimension on which bound applies
  const gsl_vector_uint* curr_bound,	///< [in] Indices of current bounds
  const gsl_vector* phys_point,		///< [in] Physical point at which to find bounds
  gsl_vector* phys_lower,		///< [out] Physical lower bounds on point
  gsl_vector* phys_upper		///< [out] Physical upper bounds on point
  );

///
/// Get padding for parameter space bounds
///
static double GetPadding(
  FlatLatticeTiling* tiling,		///< [in] Tiling state
  const size_t dimension,		///< [in] Dimension to get padding for
  const double lower,			///< [in] Normalised lower bound
  const double upper			///< [in] Normalised upper bound
  );

///
/// Flat lattice tiling bound info
///
typedef struct tagFLT_Bound {
  bool tiled;				///< Is the bound tiled, i.e. non-singular?
  FlatLatticeBound func;		///< Parameter space bound function
  void* data;				///< Arbitrary data describing parameter space
} FLT_Bound;

///
/// Flat lattice tiling status
///
typedef enum tagFLT_Status {
  FLT_S_INCOMPLETE,
  FLT_S_INITIALISED,
  FLT_S_STARTED,
  FLT_S_FINISHED
} FLT_Status;

///
/// Maximum number of parameter space bounds per dimension
///
#define MAX_BOUNDS 4

///
/// Flat lattice tiling state structure
///
struct tagFlatLatticeTiling {
  size_t dimensions;			///< Dimension of the parameter space
  FLT_Status status;			///< Status of the tiling
  FLT_Bound *bounds;			///< Array of parameter space bound info for each dimension
  FlatLatticeGenerator generator;	///< Flat tiling lattice generator function
  gsl_vector* phys_scale;		///< Normalised to physical coordinate scaling
  gsl_vector* phys_offset;		///< Normalised to physical coordinate offset
  gsl_vector* bounding_box;		///< Extent of bounding box of metric ellipse
  gsl_matrix* increment;		///< Increment vectors of the lattice tiling generator
  gsl_vector* curr_point;		///< Current lattice point
  gsl_vector_uint* curr_bound;		///< Indices of current bound on parameter space
  gsl_matrix* curr_lower;		///< Current lower bound on parameter space
  gsl_matrix* curr_upper;		///< Current upper bound on parameter space
  gsl_vector* curr_phys_point;		///< Current physical parameter-space point
  uint64_t count;			///< Total number of points generated so far
};

FlatLatticeTiling* XLALCreateFlatLatticeTiling(
  const size_t dimensions
  )
{

  // Check input
  XLAL_CHECK_NULL(dimensions > 0, XLAL_EINVAL);

  const size_t n = dimensions;

  // Allocate and initialise tiling structure
  FlatLatticeTiling* tiling = XLALCalloc(1, sizeof(FlatLatticeTiling));
  XLAL_CHECK_NULL(tiling != NULL, XLAL_ENOMEM);
  tiling->dimensions = n;
  tiling->status = FLT_S_INCOMPLETE;
  tiling->generator = NULL;
  tiling->count = 0;

  // Allocate parameter space bounds info
  tiling->bounds = XLALCalloc(n, sizeof(FLT_Bound));
  XLAL_CHECK_NULL(tiling->bounds != NULL, XLAL_ENOMEM);

  // Allocate vectors and matrices
  tiling->phys_scale = gsl_vector_alloc(n);
  XLAL_CHECK_NULL(tiling->phys_scale != NULL, XLAL_ENOMEM);
  tiling->phys_offset = gsl_vector_alloc(n);
  XLAL_CHECK_NULL(tiling->phys_offset != NULL, XLAL_ENOMEM);
  tiling->bounding_box = gsl_vector_alloc(n);
  XLAL_CHECK_NULL(tiling->bounding_box != NULL, XLAL_ENOMEM);
  tiling->increment = gsl_matrix_alloc(n, n);
  XLAL_CHECK_NULL(tiling->increment != NULL, XLAL_ENOMEM);
  tiling->curr_point = gsl_vector_alloc(n);
  XLAL_CHECK_NULL(tiling->curr_point != NULL, XLAL_ENOMEM);
  tiling->curr_bound = gsl_vector_uint_alloc(n);
  XLAL_CHECK_NULL(tiling->curr_bound != NULL, XLAL_ENOMEM);
  tiling->curr_lower = gsl_matrix_alloc(n, MAX_BOUNDS);
  XLAL_CHECK_NULL(tiling->curr_lower != NULL, XLAL_ENOMEM);
  tiling->curr_upper = gsl_matrix_alloc(n, MAX_BOUNDS);
  XLAL_CHECK_NULL(tiling->curr_upper != NULL, XLAL_ENOMEM);
  tiling->curr_phys_point = gsl_vector_alloc(n);
  XLAL_CHECK_NULL(tiling->curr_phys_point != NULL, XLAL_ENOMEM);

  return tiling;

}

void XLALDestroyFlatLatticeTiling(
  FlatLatticeTiling* tiling
  )
{

  if (tiling) {

    const size_t n = tiling->dimensions;

    gsl_error_handler_t* old_handler = gsl_set_error_handler_off();

    // Free bounds data, allowing bounds to share the same memory
    for (size_t i = 0; i < n; ++i) {
      void* data = tiling->bounds[i].data;
      if (data != NULL) {
        for (size_t j = i; j < n; ++j) {
          if (tiling->bounds[j].data == data) {
            tiling->bounds[j].data = NULL;
          }
        }
        XLALFree(data);
      }
    }
    XLALFree(tiling->bounds);

    // Free vectors and matrices
    gsl_vector_free(tiling->phys_scale);
    gsl_vector_free(tiling->phys_offset);
    gsl_vector_free(tiling->bounding_box);
    gsl_matrix_free(tiling->increment);
    gsl_vector_free(tiling->curr_point);
    gsl_vector_uint_free(tiling->curr_bound);
    gsl_matrix_free(tiling->curr_lower);
    gsl_matrix_free(tiling->curr_upper);
    gsl_vector_free(tiling->curr_phys_point);

    // Free tiling structure
    XLALFree(tiling);

    gsl_set_error_handler(old_handler);

  }

}

size_t XLALGetFlatLatticeDimensions(
  FlatLatticeTiling* tiling
  )
{
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);
  return tiling->dimensions;
}

uint64_t XLALGetFlatLatticePointCount(
  FlatLatticeTiling* tiling
  )
{
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);
  return tiling->count;
}

int XLALSetFlatLatticeBound(
  FlatLatticeTiling* tiling,
  const size_t dimension,
  const bool singular,
  const FlatLatticeBound func,
  void* data
  )
{

  // Check tiling
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status == FLT_S_INCOMPLETE, XLAL_EFAILED);

  // Check input
  XLAL_CHECK(dimension < tiling->dimensions, XLAL_ESIZE);
  XLAL_CHECK(func != NULL, XLAL_EFAULT);

  // Set the next parameter space bound
  tiling->bounds[dimension].tiled = !singular;
  tiling->bounds[dimension].func = func;
  tiling->bounds[dimension].data = data;

  return XLAL_SUCCESS;

}

int XLALSetFlatLatticeGenerator(
  FlatLatticeTiling* tiling,
  const FlatLatticeGenerator generator
  )
{

  // Check tiling
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status == FLT_S_INCOMPLETE, XLAL_EFAILED);

  // Check input
  XLAL_CHECK(generator != NULL, XLAL_EFAILED);

  // Set the flat lattice tiling generator
  tiling->generator = generator;

  return XLAL_SUCCESS;

}

int XLALSetFlatLatticeMetric(
  FlatLatticeTiling* tiling,
  const gsl_matrix* metric,
  const double max_mismatch
  )
{

  const size_t n = tiling->dimensions;

  // Check tiling
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status == FLT_S_INCOMPLETE, XLAL_EFAILED);

  // Check input
  XLAL_CHECK(metric != NULL, XLAL_EFAULT);
  XLAL_CHECK(metric->size1 == n && metric->size2 == n, XLAL_EINVAL);
  XLAL_CHECK(max_mismatch > 0, XLAL_EINVAL);

  // Check that all parameter space dimensions are bounded,
  // and count number of tiles dimensions
  size_t tn = 0;
  for (size_t i = 0; i < tiling->dimensions; ++i) {
    XLAL_CHECK(tiling->bounds[i].func != NULL, XLAL_EFAILED, "Dimension #%i is unbounded", i);
    tn += tiling->bounds[i].tiled ? 1 : 0;
  }

  // Check that the flat lattice tiling generator has been set
  XLAL_CHECK(tiling->generator != NULL, XLAL_EFAILED);

  // Initialise parameter space bound indices
  gsl_vector_uint_set_zero(tiling->curr_bound);

  // Get physical parameter space offset
  for (size_t i = 0; i < n; ++i) {

    // Get physical bounds
    gsl_vector_view phys_lower = gsl_matrix_row(tiling->curr_lower, i);
    gsl_vector_view phys_upper = gsl_matrix_row(tiling->curr_upper, i);
    GetPhysBounds(tiling, i, tiling->curr_bound, tiling->phys_offset, &phys_lower.vector, &phys_upper.vector);

    // Set physical parameter space offset
    gsl_vector_set(tiling->phys_offset, i, gsl_vector_get(&phys_lower.vector, 0));

  }

  // Check diagonal elements are positive, and calculate
  // physical parameter space scaling from metric diagonal
  gsl_vector_set_all(tiling->phys_scale, 1.0);
  for (size_t i = 0; i < n; ++i) {
    if (tiling->bounds[i].tiled) {
      const double metric_i_i = gsl_matrix_get(metric, i, i);
      XLAL_CHECK(metric_i_i > 0, XLAL_EINVAL, "metric(%zu,%zu) <= 0", i, i);
      gsl_vector_set(tiling->phys_scale, i, 1.0 / sqrt(metric_i_i));
    }
  }

  // Allocate memory
  gsl_matrix* tmetric = gsl_matrix_alloc(tn, tn);
  XLAL_CHECK(tmetric != NULL, XLAL_ENOMEM);
  gsl_matrix* tdirections = gsl_matrix_alloc(tn, tn);
  XLAL_CHECK(tdirections != NULL, XLAL_ENOMEM);
  gsl_matrix* tincrement = gsl_matrix_alloc(tn, tn);
  XLAL_CHECK(tincrement != NULL, XLAL_ENOMEM);

  // Check metric is symmetric, and copy and rescale tiled dimensions of metric
  for (size_t i = 0, ti = 0; i < n; ++i) {
    if (tiling->bounds[i].tiled) {
      const double scale_i = gsl_vector_get(tiling->phys_scale, i);
      for (size_t j = 0, tj = 0; j < n; ++j) {
        if (tiling->bounds[j].tiled) {
          const double scale_j = gsl_vector_get(tiling->phys_scale, j);
          double metric_i_j = gsl_matrix_get(metric, i, j);
          XLAL_CHECK(metric_i_j == gsl_matrix_get(metric, j, i), XLAL_EINVAL, "metric(%zu,%zu) != metric(%zu,%zu)", i, j, j, i);
          metric_i_j *= scale_i * scale_j;
          gsl_matrix_set(tmetric, ti, tj, metric_i_j);
          ++tj;
        }
      }
      ++ti;
    }
  }

  // Calculate metric ellipse bounding box
  gsl_vector* tbounding_box = XLALMetricEllipseBoundingBox(tmetric, max_mismatch);
  XLAL_CHECK(tbounding_box != NULL, XLAL_EFAILED);

  // Find orthonormalise directions with respect to subspace metric
  gsl_matrix_set_identity(tdirections);
  XLAL_CHECK(XLALOrthonormaliseWRTMetric(tdirections, tmetric) == XLAL_SUCCESS, XLAL_EFAILED);

  // Get lattice generator
  gsl_matrix* tgenerator = NULL;
  double norm_thickness = 0.0;
  XLAL_CHECK((tiling->generator)(tn, &tgenerator, &norm_thickness) == XLAL_SUCCESS, XLAL_EFAILED);

  // Transform lattice generator to square lower triangular
  gsl_matrix* sq_lwtri_generator = XLALSquareLowerTriangularLatticeGenerator(tgenerator);
  XLAL_CHECK(sq_lwtri_generator != NULL, XLAL_EFAILED);

  // Normalise lattice generator so covering radius is sqrt(mismatch)
  XLAL_CHECK(XLALNormaliseLatticeGenerator(sq_lwtri_generator, norm_thickness, sqrt(max_mismatch)) == XLAL_SUCCESS, XLAL_EFAILED);

  // Compute the increment vectors of the lattice generator along the orthogonal directions
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tdirections, sq_lwtri_generator, 0.0, tincrement);

  // Copy increment vectors and bounding box so that non-tiled dimensions are zero
  gsl_vector_set_zero(tiling->bounding_box);
  gsl_matrix_set_zero(tiling->increment);
  for (size_t i = 0, ti = 0; i < n; ++i) {
    if (tiling->bounds[i].tiled) {
      gsl_vector_set(tiling->bounding_box, i, gsl_vector_get(tbounding_box, ti));
      for (size_t j = 0, tj = 0; j < n; ++j) {
        if (tiling->bounds[j].tiled) {
          gsl_matrix_set(tiling->increment, i, j, gsl_matrix_get(tincrement, ti, tj));
          ++tj;
        }
      }
      ++ti;
    }
  }

  // Cleanup
  gsl_matrix_free(tmetric);
  gsl_matrix_free(tdirections);
  gsl_matrix_free(tincrement);
  gsl_vector_free(tbounding_box);
  gsl_matrix_free(tgenerator);
  gsl_matrix_free(sq_lwtri_generator);

  // Tiling has been fully initialised
  tiling->status = FLT_S_INITIALISED;
  XLALRestartFlatLatticeTiling(tiling);

  return XLAL_SUCCESS;

}

gsl_vector* XLALNextFlatLatticePoint(
  FlatLatticeTiling* tiling
  )
{

  const size_t n = tiling->dimensions;

  // Check tiling
  XLAL_CHECK_NULL(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(tiling->status != FLT_S_INCOMPLETE, XLAL_EFAILED);

  // If finished status, nothing more to be done!
  if (tiling->status == FLT_S_FINISHED) {
    return NULL;
  }

  // If initialised status, set and return starting point
  if (tiling->status == FLT_S_INITIALISED) {

    // Initialise parameter space bound indices
    gsl_vector_uint_set_zero(tiling->curr_bound);

    // Set parameter space bounds and starting point
    bool any_tiled = false;
    for (size_t i = 0; i < n; ++i) {

      // Get physical scales and offsets
      const double phys_scale = gsl_vector_get(tiling->phys_scale, i);
      const double phys_offset = gsl_vector_get(tiling->phys_offset, i);

      // Get physical bounds
      gsl_vector_view phys_lower = gsl_matrix_row(tiling->curr_lower, i);
      gsl_vector_view phys_upper = gsl_matrix_row(tiling->curr_upper, i);
      GetPhysBounds(tiling, i, tiling->curr_bound, tiling->curr_phys_point, &phys_lower.vector, &phys_upper.vector);

      // Convert physical bounds to normalised bounds
      gsl_vector_add_constant(&phys_lower.vector, -phys_offset);
      gsl_vector_scale(&phys_lower.vector, 1.0/phys_scale);
      gsl_vector_add_constant(&phys_upper.vector, -phys_offset);
      gsl_vector_scale(&phys_upper.vector, 1.0/phys_scale);

      // Get bounds
      const size_t bound = gsl_vector_uint_get(tiling->curr_bound, i);
      const double lower = gsl_matrix_get(tiling->curr_lower, i, bound);

      // Determine whether any bounds are tiled
      any_tiled |= tiling->bounds[i].tiled;

      // Initialise current point
      const double point = lower;
      gsl_vector_set(tiling->curr_point, i, point);

      // Update current physical point
      const double phys_point = (point * phys_scale) + phys_offset;
      gsl_vector_set(tiling->curr_phys_point, i, phys_point);

    }

    // Subtract padding
    for (size_t i = 0; i < n; ++i) {
      const size_t bound = gsl_vector_uint_get(tiling->curr_bound, i);
      const double lower = gsl_matrix_get(tiling->curr_lower, i, bound);
      const double upper = gsl_matrix_get(tiling->curr_upper, i, bound);
      const double padding = GetPadding(tiling, i, lower, upper);
      const double point = gsl_vector_get(tiling->curr_point, i);
      gsl_vector_set(tiling->curr_point, i, point - padding);
    }

    // Update current physical point
    gsl_vector_memcpy(tiling->curr_phys_point, tiling->curr_point);
    gsl_vector_mul(tiling->curr_phys_point, tiling->phys_scale);
    gsl_vector_add(tiling->curr_phys_point, tiling->phys_offset);

    // Initialise count
    tiling->count = 1;

    // If no bounds are tiled, there is only one template point!
    // Otherwise, tiling has been started.
    tiling->status = any_tiled ? FLT_S_STARTED : FLT_S_FINISHED;

    return tiling->curr_phys_point;

  }

  // Otherwise started status: loop until the next point is found
  size_t i = n, ir;
  while (true) {

    // If dimension index is now zero, we're done!
    if (i == 0) {
      tiling->status = FLT_S_FINISHED;
      return NULL;
    }

    // Decrement current dimension index
    --i;

    // Return point to lower bound in higher dimensions
    ir = i + 1;

    // Get current bound index
    size_t bound = gsl_vector_uint_get(tiling->curr_bound, i);

    // If dimension is tiled...
    if (tiling->bounds[i].tiled) {

      // Get increment vector
      gsl_vector_view increment = gsl_matrix_column(tiling->increment, i);

      // Increment current point along index
      gsl_vector_add(tiling->curr_point, &increment.vector);

      // Update current physical point
      gsl_vector_memcpy(tiling->curr_phys_point, tiling->curr_point);
      gsl_vector_mul(tiling->curr_phys_point, tiling->phys_scale);
      gsl_vector_add(tiling->curr_phys_point, tiling->phys_offset);

      // Get current point, bounds and padding
      const double point = gsl_vector_get(tiling->curr_point, i);
      const double lower = gsl_matrix_get(tiling->curr_lower, i, bound);
      const double upper = gsl_matrix_get(tiling->curr_upper, i, bound);
      const double padding = GetPadding(tiling, i, lower, upper);

      // If point is not out of bounds, we have found a template point
      if (point <= upper + padding) {
        break;
      }

    }

    // Increment bound index
    ++bound;

    if (bound < MAX_BOUNDS) {

      // Get bounds
      const double lower = gsl_matrix_get(tiling->curr_lower, i, bound);
      const double upper = gsl_matrix_get(tiling->curr_upper, i, bound);

      if (gsl_isnan(lower) && gsl_isnan(upper)) {

        // If no more bounds, reset bound index in this and higher dimensions
        for (size_t j = i; j < n; ++j) {
          gsl_vector_uint_set(tiling->curr_bound, j, 0);
        }

      } else {

        // Set current bound index
        gsl_vector_uint_set(tiling->curr_bound, i, bound);

        // Return point to new lower bound in this dimension
        ir = i;

        // Found a template point
        break;

      }

    }

    // Move on to lower dimensions
    continue;

  }

  // Return point to lower bound in appropriate dimensions
  for (; ir < n; ++ir) {

    // Get current bound index
    const size_t bound = gsl_vector_uint_get(tiling->curr_bound, ir);

    // Get new physical bounds if required
    if (bound == 0) {

      // Get physical scales and offsets
      const double phys_scale = gsl_vector_get(tiling->phys_scale, ir);
      const double phys_offset = gsl_vector_get(tiling->phys_offset, ir);

      // Get physical bounds
      gsl_vector_view phys_lower = gsl_matrix_row(tiling->curr_lower, ir);
      gsl_vector_view phys_upper = gsl_matrix_row(tiling->curr_upper, ir);
      GetPhysBounds(tiling, ir, tiling->curr_bound, tiling->curr_phys_point, &phys_lower.vector, &phys_upper.vector);

      // Convert physical bounds to normalised bounds
      gsl_vector_add_constant(&phys_lower.vector, -phys_offset);
      gsl_vector_scale(&phys_lower.vector, 1.0/phys_scale);
      gsl_vector_add_constant(&phys_upper.vector, -phys_offset);
      gsl_vector_scale(&phys_upper.vector, 1.0/phys_scale);

    }

    // Get lower bound
    const double lower = gsl_matrix_get(tiling->curr_lower, ir, bound);

    // If dimension is tiled...
    if (tiling->bounds[ir].tiled) {

      // Get increment vector
      gsl_vector_view increment = gsl_matrix_column(tiling->increment, ir);

      // Get upper bound and padding
      const double upper = gsl_matrix_get(tiling->curr_upper, ir, bound);
      const double padding = GetPadding(tiling, ir, lower, upper);

      // Calculate the distance from current point to the lower bound, in integer number of increments
      const double point = gsl_vector_get(tiling->curr_point, ir);
      const double dist = floor((lower - padding - point) / gsl_vector_get(&increment.vector, ir));

      // Move point back to lower bound
      gsl_blas_daxpy(dist, &increment.vector, tiling->curr_point);

    } else {

      // Otherwise set point to lower bound
      gsl_vector_set(tiling->curr_point, ir, lower);

    }

    // Update current physical point
    gsl_vector_memcpy(tiling->curr_phys_point, tiling->curr_point);
    gsl_vector_mul(tiling->curr_phys_point, tiling->phys_scale);
    gsl_vector_add(tiling->curr_phys_point, tiling->phys_offset);

  }

  // Template was found, so increase count
  ++tiling->count;

  return tiling->curr_phys_point;

}

size_t XLALNextFlatLatticePoints(
  FlatLatticeTiling* tiling,
  gsl_matrix* points,
  const bool fill_last
  )
{

  const size_t n = tiling->dimensions;

  // Check tiling
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_VAL(0, tiling->status != FLT_S_INCOMPLETE, XLAL_EFAILED);

  // Check input
  XLAL_CHECK_VAL(0, points != NULL, XLAL_EFAULT);
  XLAL_CHECK_VAL(0, points->size1 == n, XLAL_EFAULT);

  // Fill 'points' matrix columns with flat lattice tiling points
  size_t i = 0;
  gsl_vector* point = NULL;
  while (i < points->size2) {

    // Get next tiling point, checking for errors
    point = XLALNextFlatLatticePoint(tiling);
    XLAL_CHECK_VAL(0, xlalErrno == 0, XLAL_EFAILED);

    // If no point return, no more points available
    if (point == NULL) {
      break;
    }

    // Copy the point to the next available column of 'points'
    gsl_vector_view p = gsl_matrix_column(points, i);
    gsl_vector_memcpy(&p.vector, point);
    ++i;

  }

  if (fill_last && i > 0) {

    // Copy last template point to remaining columns of 'points'
    gsl_vector_view q = gsl_matrix_column(points, i-1);
    for (size_t j = i; j < points->size2; ++j) {
      gsl_vector_view p = gsl_matrix_column(points, i);
      gsl_vector_memcpy(&p.vector, &q.vector);
    }

  }

  // Return the number of points filled in 'points'
  return i;

}

int XLALRestartFlatLatticeTiling(
  FlatLatticeTiling* tiling
  )
{

  // Check tiling
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status != FLT_S_INCOMPLETE, XLAL_EFAILED);

  // Restart tiling
  tiling->status = FLT_S_INITIALISED;
  tiling->count = 0;

  return XLAL_SUCCESS;

}

uint64_t XLALCountTotalFlatLatticePoints(
  FlatLatticeTiling* tiling
  )
{

  // Check tiling
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_VAL(0, tiling->status != FLT_S_INCOMPLETE, XLAL_EFAILED);

  // Iterate over all templates
  while (XLALNextFlatLatticePoint(tiling) != NULL);
  XLAL_CHECK_VAL(0, xlalErrno == 0, XLAL_EFAILED);

  // Save the template count
  uint64_t count = tiling->count;

  // Restart tiling
  XLALRestartFlatLatticeTiling(tiling);

  // Return the template count
  return count;

}

int XLALGenerateRandomFlatLatticePoints(
  FlatLatticeTiling* tiling,
  RandomParams* randpar,
  gsl_matrix* randpoints
  )
{

  const size_t n = tiling->dimensions;

  // Check tiling
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status != FLT_S_INCOMPLETE, XLAL_EFAILED);

  // Check input
  XLAL_CHECK(randpar != NULL, XLAL_EFAULT);
  XLAL_CHECK(randpoints != NULL, XLAL_EFAULT);
  XLAL_CHECK(randpoints->size1 == n, XLAL_ESIZE);

  // Allocate temporary bound index and physical bound vectors
  gsl_vector_uint* curr_bound = gsl_vector_uint_alloc(n);
  XLAL_CHECK(curr_bound != NULL, XLAL_ENOMEM);
  gsl_vector* phys_lower = gsl_vector_alloc(MAX_BOUNDS);
  XLAL_CHECK(phys_lower != NULL, XLAL_ENOMEM);
  gsl_vector* phys_width = gsl_vector_alloc(MAX_BOUNDS);
  XLAL_CHECK(phys_width != NULL, XLAL_ENOMEM);

  // Create random points in tiling parameter space
  for (size_t j = 0; j < randpoints->size2; ++j) {
    gsl_vector_view point = gsl_matrix_column(randpoints, j);
    for (size_t i = 0; i < n; ++i) {

      // Get physical bounds
      GetPhysBounds(tiling, i, curr_bound, &point.vector, phys_lower, phys_width);
      gsl_vector_sub(phys_width, phys_lower);

      // Get total bounds width
      double phys_total_width = 0;
      size_t max_bounds = 0;
      while (max_bounds < MAX_BOUNDS) {
        const double width = gsl_vector_get(phys_width, max_bounds);
        if (gsl_isnan(width)) {
          break;
        }
        phys_total_width += width;
        ++max_bounds;
      }

      // Generate random number
      const double u = XLALUniformDeviate(randpar);

      double p;
      size_t bound = 0;
      if (tiling->bounds[i].tiled) {

        // Generate random point within total bounds widths
        p = u * phys_total_width;

        // Convert point to be within parameter space bounds
        while (bound + 1 < max_bounds) {
          const double width = gsl_vector_get(phys_width, bound);
          if (p <= width) {
            break;
          }
          p -= width;
          ++bound;
        }
        p += gsl_vector_get(phys_lower, bound);

      } else {

        // Generate random bound index
        bound = (size_t)floor(u * max_bounds);

        // Get point from random bound
        p = gsl_vector_get(phys_lower, bound);

      }

      // Set parameter space point and bound index
      gsl_vector_set(&point.vector, i, p);
      gsl_vector_uint_set(curr_bound, i, bound);

    }

  }

  // Cleanup
  gsl_vector_uint_free(curr_bound);
  gsl_vector_free(phys_lower);
  gsl_vector_free(phys_width);

  return XLAL_SUCCESS;

}

int XLALCubicLatticeGenerator(
  const size_t dimensions,
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
  const size_t dimensions,
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

static void ConstantBound(
  const gsl_vector_uint* bound UNUSED,
  const gsl_vector* point UNUSED,
  const void* data,
  gsl_vector* lower,
  gsl_vector* upper
  )
{

  // Set constant lower and upper bounds
  gsl_vector_set(lower, 0, ((const double*)data)[0]);
  if (upper) {
    gsl_vector_set(upper, 0, ((const double*)data)[1]);
  }

}

int XLALSetFlatLatticeConstantBound(
  FlatLatticeTiling* tiling,
  size_t dimension,
  double lower,
  double upper
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(lower <= upper, XLAL_EINVAL);

  // Allocate and set bounds data
  double* data = XLALCalloc(2, sizeof(double));
  XLAL_CHECK(data != NULL, XLAL_ENOMEM);
  data[0] = lower;
  data[1] = upper;

  // Set parameter space bound
  XLAL_CHECK(XLALSetFlatLatticeBound(tiling, dimension, lower == upper, ConstantBound, (void*)data) == XLAL_SUCCESS, XLAL_EFAILED);

  return XLAL_SUCCESS;

}

static void GetPhysBounds(
  FlatLatticeTiling* tiling,
  const size_t dimension,
  const gsl_vector_uint* curr_bound,
  const gsl_vector* phys_point,
  gsl_vector* phys_lower,
  gsl_vector* phys_upper
  )
{

  // Get the appropriate bound dimension
  FLT_Bound* bound = &tiling->bounds[dimension];

  // Initialise bound vectors
  gsl_vector_set_all(phys_lower, GSL_NAN);
  gsl_vector_set_all(phys_upper, GSL_NAN);

  // Pass upper bound vector only for tiled bounds
  gsl_vector* phys_upper_tiled = bound->tiled ? phys_upper : NULL;

  // Call parameter space bounds function, with
  // view of physical point only in lower dimensions
  if (dimension == 0) {
    (bound->func)(NULL, NULL, bound->data, phys_lower, phys_upper_tiled);
  } else {
    gsl_vector_uint_const_view curr_bound_view = gsl_vector_uint_const_subvector(curr_bound, 0, dimension);
    gsl_vector_const_view phys_point_view = gsl_vector_const_subvector(phys_point, 0, dimension);
    (bound->func)(&curr_bound_view.vector, &phys_point_view.vector, bound->data, phys_lower, phys_upper_tiled);
  }

}

static double GetPadding(
  FlatLatticeTiling* tiling,
  const size_t dimension,
  const double lower UNUSED,
  const double upper UNUSED
  )
{
  if (!tiling->bounds[dimension].tiled) {
    return 0.0;
  }
  return gsl_vector_get(tiling->bounding_box, dimension);
}

gsl_vector* XLALMetricEllipseBoundingBox(
  gsl_matrix* metric,
  const double max_mismatch
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

int XLALOrthonormaliseWRTMetric(
  gsl_matrix* matrix,
  const gsl_matrix* metric
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

gsl_matrix* XLALSquareLowerTriangularLatticeGenerator(
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

int XLALNormaliseLatticeGenerator(
  gsl_matrix* generator,
  const double norm_thickness,
  const double covering_radius
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

///
/// Workspace for computing the nearest template to a set of injections
///
struct tagNearestTemplateWorkspace {
  gsl_matrix* metric;			///< Parameter space metric
  gsl_matrix* metric_templates;		///< Templates points, multiplied by metric
  gsl_vector* dot_templates;		///< Dot product of templates
  gsl_matrix* injections;		///< Injection points
  gsl_vector* dot_injections;		///< Dot product of injections
  gsl_matrix* cross_terms;		///< Cross terms in distance between templates and injections
  gsl_vector* distances;		///< Distances between templates and injections
};

NearestTemplateWorkspace* XLALCreateNearestTemplateWorkspace(
  const gsl_matrix* metric,
  const size_t num_templates,
  const size_t num_injections
  )
{

  // Check input
  XLAL_CHECK_NULL(metric != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(metric->size1 == metric->size2, XLAL_ESIZE);

  // Check metric is symmetric, that metric diagonal elements are positive or zero,
  // and that rows and columns passing through zero diagonals are also zero
  for (size_t i = 0; i < metric->size1; ++i) {
    XLAL_CHECK_NULL(gsl_matrix_get(metric, i, i) >= 0.0, XLAL_EINVAL);
    for (size_t j = 0; j < metric->size2; ++j) {
      if (gsl_matrix_get(metric, i, i) == 0.0 || gsl_matrix_get(metric, i, i) == 0.0) {
        XLAL_CHECK_NULL(gsl_matrix_get(metric, i, j) == 0.0, XLAL_EINVAL);
      } else {
        XLAL_CHECK_NULL(gsl_matrix_get(metric, i, j) == gsl_matrix_get(metric, j, i), XLAL_EINVAL);
      }
    }
  }

  // Allocate and initialise workspace
  NearestTemplateWorkspace* wksp = XLALCalloc(1, sizeof(NearestTemplateWorkspace));
  XLAL_CHECK_NULL(wksp != NULL, XLAL_ENOMEM);

  // Allocate workspace memory
  wksp->metric = gsl_matrix_alloc(metric->size1, metric->size2);
  XLAL_CHECK_NULL(wksp->metric != NULL, XLAL_ENOMEM);
  wksp->metric_templates = gsl_matrix_alloc(wksp->metric->size1, num_templates);
  XLAL_CHECK_NULL(wksp->metric_templates != NULL, XLAL_ENOMEM);
  wksp->dot_templates = gsl_vector_alloc(num_templates);
  XLAL_CHECK_NULL(wksp->dot_templates != NULL, XLAL_ENOMEM);
  wksp->injections = gsl_matrix_alloc(wksp->metric->size1, num_injections);
  XLAL_CHECK_NULL(wksp->injections != NULL, XLAL_ENOMEM);
  wksp->dot_injections = gsl_vector_alloc(num_injections);
  XLAL_CHECK_NULL(wksp->dot_injections != NULL, XLAL_ENOMEM);
  wksp->cross_terms = gsl_matrix_alloc(num_injections, num_templates);
  XLAL_CHECK_NULL(wksp->cross_terms != NULL, XLAL_ENOMEM);
  wksp->distances = gsl_vector_alloc(num_templates);
  XLAL_CHECK_NULL(wksp->distances != NULL, XLAL_ENOMEM);

  // Copy metric
  gsl_matrix_memcpy(wksp->metric, metric);

  return wksp;

}

void XLALDestroyNearestTemplateWorkspace(
  NearestTemplateWorkspace* wksp
  )
{

  if (wksp) {

    gsl_error_handler_t* old_handler = gsl_set_error_handler_off();

    // Cleanup
    gsl_matrix_free(wksp->metric);
    gsl_matrix_free(wksp->metric_templates);
    gsl_vector_free(wksp->dot_templates);
    gsl_matrix_free(wksp->injections);
    gsl_vector_free(wksp->dot_injections);
    gsl_matrix_free(wksp->cross_terms);
    gsl_vector_free(wksp->distances);
    XLALFree(wksp);

    gsl_set_error_handler(old_handler);

  }

}

int XLALUpdateWorkspaceTemplates(
  NearestTemplateWorkspace* wksp,
  const gsl_matrix* templates,
  gsl_vector_uint* nearest_template
  )
{

  // Check input
  XLAL_CHECK(wksp != NULL, XLAL_EFAULT);
  XLAL_CHECK(templates != NULL, XLAL_EFAULT);
  XLAL_CHECK(templates->size1 == wksp->metric->size1, XLAL_ESIZE);
  XLAL_CHECK(templates->size2 == wksp->dot_templates->size, XLAL_ESIZE);
  XLAL_CHECK(nearest_template != NULL, XLAL_EFAULT);
  XLAL_CHECK(nearest_template->size == wksp->dot_templates->size, XLAL_EFAULT);

  // Compute dot product of templates
  gsl_blas_dsymm(CblasLeft, CblasUpper, 1.0, wksp->metric, templates, 0.0, wksp->metric_templates);
  gsl_matrix_mul_elements(wksp->metric_templates, templates);
  gsl_vector_set_zero(wksp->dot_templates);
  for (size_t i = 0; i < wksp->metric->size1; ++i) {
    gsl_vector_view v = gsl_matrix_row(wksp->metric_templates, i);
    gsl_vector_add(wksp->dot_templates, &v.vector);
  }

  // Copy templates and multiply by metric
  gsl_blas_dsymm(CblasLeft, CblasUpper, 1.0, wksp->metric, templates, 0.0, wksp->metric_templates);

  // Initialise nearest template index
  gsl_vector_uint_set_all(nearest_template, -1);

  return XLAL_SUCCESS;

}

int XLALUpdateWorkspaceInjections(
  NearestTemplateWorkspace* wksp,
  const gsl_matrix* injections,
  gsl_vector* min_distance
  )
{

  // Check input
  XLAL_CHECK(wksp != NULL, XLAL_EFAULT);
  XLAL_CHECK(injections != NULL, XLAL_EFAULT);
  XLAL_CHECK(injections->size1 == wksp->metric->size1, XLAL_ESIZE);
  XLAL_CHECK(injections->size2 == wksp->dot_injections->size, XLAL_ESIZE);
  XLAL_CHECK(min_distance != NULL, XLAL_EFAULT);
  XLAL_CHECK(min_distance->size == wksp->dot_injections->size, XLAL_EFAULT);

  // Compute dot product of injections
  gsl_blas_dsymm(CblasLeft, CblasUpper, 1.0, wksp->metric, injections, 0.0, wksp->injections);
  gsl_matrix_mul_elements(wksp->injections, injections);
  gsl_vector_set_zero(wksp->dot_injections);
  for (size_t i = 0; i < wksp->metric->size1; ++i) {
    gsl_vector_view v = gsl_matrix_row(wksp->injections, i);
    gsl_vector_add(wksp->dot_injections, &v.vector);
  }

  // Copy injections
  gsl_matrix_memcpy(wksp->injections, injections);

  // Initialise minimum distances
  gsl_vector_set_all(min_distance, GSL_POSINF);

  return XLAL_SUCCESS;

}

int XLALUpdateNearestTemplateToInjections(
  NearestTemplateWorkspace* wksp,
  gsl_vector* min_distance,
  gsl_vector_uint* nearest_template
  )
{

  // Check input
  XLAL_CHECK(wksp != NULL, XLAL_EFAULT);
  XLAL_CHECK(min_distance != NULL, XLAL_EFAULT);
  XLAL_CHECK(min_distance->size == wksp->dot_injections->size, XLAL_EFAULT);
  XLAL_CHECK(nearest_template != NULL, XLAL_EFAULT);
  XLAL_CHECK(nearest_template->size == wksp->dot_templates->size, XLAL_EFAULT);

  // Compute cross terms in distance between templates and injections
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, -2.0, wksp->injections, wksp->metric_templates, 0.0, wksp->cross_terms);

  // Find closest template to each injection
  for (size_t i = 0; i < wksp->dot_injections->size; ++i) {

    // Start with template dot products
    gsl_vector_memcpy(wksp->distances, wksp->dot_templates);

    // Add cross terms for this injection
    gsl_vector_view v = gsl_matrix_row(wksp->cross_terms, i);
    gsl_vector_add(wksp->distances, &v.vector);

    // Find smallest vector element
    size_t i_min = gsl_vector_min_index(wksp->distances);
    double mu_min = gsl_vector_get(wksp->distances, i_min);

    // Compute minimum distance by add injection dot product
    mu_min += gsl_vector_get(wksp->dot_injections, i);

    // Update minimum distance vector
    if (mu_min < gsl_vector_get(min_distance, i)) {
      gsl_vector_set(min_distance, i, mu_min);
      gsl_vector_uint_set(nearest_template, i, i_min);
    }

  }

  return XLAL_SUCCESS;

}
