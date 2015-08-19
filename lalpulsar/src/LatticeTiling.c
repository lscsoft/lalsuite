//
// Copyright (C) 2007, 2008, 2012, 2014, 2015 Karl Wette
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307 USA
//

#include <config.h>
#include <fenv.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <lal/LatticeTiling.h>
#include <lal/LALStdio.h>
#include <lal/MetricUtils.h>

#include "GSLHelpers.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

///
/// Lattice tiling parameter-space bound for one dimension.
///
typedef struct tagLT_Bound {
  bool is_tiled;			///< True if the dimension is tiled, false if it is a single point
  LatticeTilingBound func;		///< Parameter space bound function
  size_t data_len;			///< Length of arbitrary data describing parameter-space bounds
  void *data_lower;			///< Arbitrary data describing lower parameter-space bound
  void *data_upper;			///< Arbitrary data describing upper parameter-space bound
  UINT4 pad_lower;			///< Multiple of lower parameter space bound padding to apply
  UINT4 pad_upper;			///< Multiple of upper parameter space bound padding to apply
} LT_Bound;

///
/// Lattice tiling index trie for one dimension.
///
typedef struct tagLT_IndexTrie LT_IndexTrie;
struct tagLT_IndexTrie {
  INT8 int_lower;			///< Lower integer point bound in this dimension
  INT8 int_upper;			///< Upper integer point bound in this dimension
  UINT8 seqidx;				///< Sequential lattice tiling index up to this dimension
  LT_IndexTrie *next;			///< Pointer to array of index tries for the next-highest dimension
};

struct tagLatticeTiling {
  size_t ndim;				///< Number of parameter-space dimensions
  LT_Bound *bounds;			///< Array of parameter-space bound info for each dimension
  size_t tiled_ndim;			///< Number of tiled parameter-space dimensions
  size_t *tiled_idx;			///< Index to tiled parameter-space dimensions
  TilingLattice lattice;		///< Type of lattice to generate tiling with
  gsl_vector *phys_bbox;		///< Metric ellipse bounding box
  gsl_vector *phys_origin;		///< Parameter-space origin in physical coordinates
  gsl_matrix *int_from_phys;		///< Transform to generating integers from physical coordinates
  gsl_matrix *phys_from_int;		///< Transform to physical coordinates from generating integers
  gsl_matrix *tiled_generator;		///< Lattice generator matrix in tiled dimensions
  LatticeTilingStats *stats;		///< Array of parameter-space statistics for each dimension
};

struct tagLatticeTilingIterator {
  LatticeTiling *tiling;		///< Lattice tiling
  size_t itr_ndim;			///< Number of parameter-space dimensions to iterate over
  size_t tiled_itr_ndim;		///< Number of tiled parameter-space dimensions to iterate over
  TilingIteratorFlags flags;		///< Flags which determine behaviour of iterator
  int state;				///< Iterator state: 0=initialised, 1=in progress, 2=finished
  gsl_vector *phys_point;		///< Current lattice point in physical coordinates
  INT8 *int_point;			///< Current lattice point in generating integers
  INT8 *int_lower;			///< Current lower parameter-space bound in generating integers
  INT8 *int_upper;			///< Current upper parameter-space bound in generating integers
  INT8 *direction;			///< Direction of iteration in each tiled parameter-space dimension
  UINT8 index;				///< Index of current lattice tiling point
};

struct tagLatticeTilingLocator {
  LatticeTiling *tiling;		///< Lattice tiling
  size_t ndim;				///< Number of parameter-space dimensions
  size_t tiled_ndim;			///< Number of tiled parameter-space dimensions
  LT_IndexTrie *index_trie;		///< Trie for locating unique index of nearest point
};

///
/// Zero out the strictly upper triangular part of the matrix \c A.
///
static void LT_ZeroStrictUpperTriangle(gsl_matrix *A)
{
  for (size_t i = 0; i < A->size1; ++i) {
    for (size_t j = i + 1; j < A->size2; ++j) {
      gsl_matrix_set(A, i, j, 0.0);
    }
  }
}

///
/// Reverse the order of both the rows and columns of the matrix \c A.
///
static void LT_ReverseOrderRowsCols(gsl_matrix *A)
{
  for (size_t i = 0; i < A->size1 / 2; ++i) {
    gsl_matrix_swap_rows(A, i, A->size1 - i - 1);
  }
  for (size_t j = 0; j < A->size2 / 2; ++j) {
    gsl_matrix_swap_columns(A, j, A->size2 - j - 1);
  }
}

///
/// Return the parameter-space bounds on a given dimension.
///
static void LT_GetBounds(
  const LatticeTiling *tiling,		///< [in] Lattice tiling
  const size_t dim,			///< [in] Dimension on which bound applies
  const gsl_vector *phys_point,		///< [in] Physical point at which to find bounds
  double *phys_lower,			///< [out] Lower parameter-space bound
  double *phys_upper			///< [out] Upper parameter-space bound
  )
{

  // Get bound information for this dimension
  const LT_Bound *bound = &tiling->bounds[dim];

  // Get view of first (dimension) dimensions of physical point
  gsl_vector_const_view phys_point_subv_view = gsl_vector_const_subvector(phys_point, 0, GSL_MAX(1, dim));
  const gsl_vector *phys_point_subv = (dim == 0) ? NULL : &phys_point_subv_view.vector;

  // Get lower parameter-space bound
  *phys_lower = (bound->func)(bound->data_lower, dim, phys_point_subv);

  if (bound->is_tiled) {

    // Get upper parameter-space bound
    *phys_upper = (bound->func)(bound->data_upper, dim, phys_point_subv);

  } else {

    // Set upper bound to lower bound
    *phys_upper = *phys_lower;

  }

}

///
/// Find the extrema of the parameter-space bounds, by sampling the bounds around the current point.
///
static void LT_FindBoundExtrema(
  const LatticeTiling *tiling,		///< [in] Lattice tiling
  const size_t i,			///< [in] Current dimension of parameter space
  const size_t dim,			///< [in] Dimension in which to record bound extrema
  gsl_vector *phys_point,		///< [out] Current physical point being bounded
  double *phys_min_lower,		///< [out] Minimum lower bound on parameter space
  double *phys_max_upper		///< [out] Maximum upper bound on parameter space
  )
{

  if (i < dim) {

    if (tiling->bounds[i].is_tiled) {

      // Get the vector pointing to neighbouring points of the current point in this dimension
      gsl_vector_const_view phys_from_int_i_view = gsl_matrix_const_row(tiling->phys_from_int, i);

      for (int d = -1; d <= 1; ++d) {

        // Move current point half-way towards neighbouring point in the direction 'd'
        gsl_blas_daxpy(0.5 * d, &phys_from_int_i_view.vector, phys_point);

        // Find bound extrema in higher dimensions
        LT_FindBoundExtrema(tiling, i + 1, dim, phys_point, phys_min_lower, phys_max_upper);

        // Reset current point back to previous value
        gsl_blas_daxpy(-0.5 * d, &phys_from_int_i_view.vector, phys_point);

      }

    } else {

      // Find bound extrema in higher dimensions
      LT_FindBoundExtrema(tiling, i + 1, dim, phys_point, phys_min_lower, phys_max_upper);

    }

  } else {

    // Get the physical bounds on the current dimension
    double phys_lower = 0.0, phys_upper = 0.0;
    LT_GetBounds(tiling, i, phys_point, &phys_lower, &phys_upper);

    // Record the minimum lower bound and maximum upper bound
    if (phys_lower < *phys_min_lower) {
      *phys_min_lower = phys_lower;
    }
    if (phys_upper > *phys_max_upper) {
      *phys_max_upper = phys_upper;
    }

  }

}

///
/// Return the extrema of the parameter-space bounds on a given dimension.
///
static void LT_GetExtremalBounds(
  const LatticeTiling *tiling,		///< [in] Lattice tiling
  const size_t dim,			///< [in] Dimension on which bound applies
  const gsl_vector *phys_point,		///< [in] Physical point at which to find bounds
  double *phys_lower,			///< [out] Lower parameter-space bound
  double *phys_upper			///< [out] Upper parameter-space bound
  )
{

  // Get the physical bounds on the current dimension
  LT_GetBounds(tiling, dim, phys_point, phys_lower, phys_upper);

  if (tiling->bounds[dim].is_tiled) {

    if (dim > 0) {

      // Create a copy of current physical point for use in finding bound extrema
      double test_phys_point_array[phys_point->size];
      gsl_vector_view test_phys_point_view = gsl_vector_view_array(test_phys_point_array, phys_point->size);
      gsl_vector_memcpy(&test_phys_point_view.vector, phys_point);

      // Find the extreme values of the physical bounds
      LT_FindBoundExtrema(tiling, 0, dim, &test_phys_point_view.vector, phys_lower, phys_upper);

    }

    // Add padding of (multiples of) half the metric ellipse bounding box in this dimension
    const double phys_bbox_dim = 0.5 * gsl_vector_get(tiling->phys_bbox, dim);
    *phys_lower -= tiling->bounds[dim].pad_lower * phys_bbox_dim;
    *phys_upper += tiling->bounds[dim].pad_upper * phys_bbox_dim;

  }

}

///
/// Fast-forward a lattice tiling iterator over its highest tiled dimension.
/// Return the number of points fast-forwarded over.
///
static INT8 LT_FastForwardIterator(
  LatticeTilingIterator *itr		///< [in] Lattice tiling iterator
  )
{

  // Get indexes of highest tiled dimension
  const size_t tj = itr->tiling->tiled_ndim - 1;
  const size_t j = itr->tiling->tiled_idx[tj];

  // Get current iteration direction in highest tiled dimension
  const INT8 direction = itr->direction[tj];

  // Get difference from integer point to upper/lower bound, depending on direction
  const INT8 d = ((direction > 0) ? itr->int_upper[tj] : itr->int_lower[tj]) - itr->int_point[tj];

  // Set point in highest tiled dimension to upper/lower bound, so that the next call
  // to XLALNextLatticeTilingPoint() will advance the next-highest tiled dimension
  itr->int_point[tj] += d;
  {
    const double phys_point_j = gsl_vector_get(itr->phys_point, j);
    const double phys_from_int_j_j = gsl_matrix_get(itr->tiling->phys_from_int, j, j);
    gsl_vector_set(itr->phys_point, j, phys_point_j + phys_from_int_j_j * d);
  }

  // Update iterator index
  itr->index += labs(d);

  return labs(d);

}

///
/// Free memory pointed to by an index trie. The trie itself should be freed by the caller.
///
static void LT_FreeIndexTrie(
  LT_IndexTrie *trie			///< [in] Pointer to array of index tries
  )
{
  LT_IndexTrie *next = trie->next;
  if (next != NULL) {
    for (INT8 i = trie->int_lower; i <= trie->int_upper; ++i) {
      LT_FreeIndexTrie(next++);
    }
    XLALFree(trie->next);
  }
}

///
/// Find the nearest point within the parameter-space bounds of the lattice tiling, by polling
/// the neighbours of an 'original' nearest point found by XLALNearestLatticeTilingPoints().
///
static void LT_PollIndexTrie(
  const LatticeTiling *tiling,		///< [in] Lattice tiling
  const LT_IndexTrie *trie,		///< [in] Lattice tiling index trie
  const size_t ti,			///< [in] Current depth of the trie
  const INT8 *original_int,		///< [in] Original nearest point
  INT8 *poll_int,			///< [in] Neighbouring point currently being polled
  double *poll_min_distance,		///< [in] Minimum distance to neighbouring point found so far
  INT8 *nearest_int			///< [in] New nearest point found by polling
  )
{

  const size_t tn = tiling->tiled_ndim;

  // Get integer lower and upper bounds
  const INT8 int_lower = trie->int_lower;
  const INT8 int_upper = trie->int_upper;

  // Poll points within 1 of original nearest point, but within bounds
  const INT8 poll_lower = GSL_MAX(int_lower, GSL_MIN(original_int[ti] - 1, int_upper));
  const INT8 poll_upper = GSL_MAX(int_lower, GSL_MIN(original_int[ti] + 1, int_upper));

  for (poll_int[ti] = poll_lower; poll_int[ti] <= poll_upper; ++poll_int[ti]) {

    // Continue polling in higher dimensions
    if (ti + 1 < tn) {
      const LT_IndexTrie *next = &trie->next[poll_int[ti] - trie->int_lower];
      LT_PollIndexTrie(tiling, next, ti + 1, original_int, poll_int, poll_min_distance, nearest_int);
      continue;
    }

    // Compute distance between original and poll point with respect to lattice generator
    double poll_distance = 0;
    for (size_t tj = 0; tj < tn; ++tj) {
      const double diff_j = original_int[tj] - poll_int[tj];
      for (size_t tk = 0; tk < tn; ++tk) {
        const double diff_k = original_int[tk] - poll_int[tk];
        const double generator_j_k = gsl_matrix_get(tiling->tiled_generator, tj, tk);
        poll_distance += generator_j_k * diff_j * diff_k;
      }
    }

    // If distance is smaller than minimum distance, record current poll point as nearest
    if (poll_distance < *poll_min_distance) {
      *poll_min_distance = poll_distance;
      memcpy(nearest_int, poll_int, tn * sizeof(nearest_int[0]));
    }

  }

}

///
/// Print one level of a lattice tiling index trie.
///
static void LT_PrintIndexTrie(
  const LatticeTiling *tiling,		///< [in] Lattice tiling
  const LT_IndexTrie *trie,		///< [in] Lattice tiling index trie
  const size_t ti,			///< [in] Current depth of the trie
  FILE *file,				///< [in] File pointer to print trie to
  INT8 int_lower[]			///< [in] Current integer lower bound
  )
{

  const size_t tn = tiling->tiled_ndim;

  // Print indentation
  for (size_t s = 0; s <= ti; ++s) {
    fprintf(file, "   ");
  }

  // Return if 'trie' is NULL (which should never happen)
  if (trie == NULL) {
    fprintf(file, "ERROR: 'trie' is NULL\n");
    return;
  }

  // Set 'i'th integer lower bound to 'trie' lower bound, then
  // transform to physical coordinates from generating integers
  int_lower[ti] = trie->int_lower;
  const size_t i = tiling->tiled_idx[ti];
  double phys_lower = gsl_vector_get(tiling->phys_origin, i);
  for (size_t tj = 0; tj < tn; ++tj) {
    const size_t j = tiling->tiled_idx[tj];
    const double phys_from_int_i_j = gsl_matrix_get(tiling->phys_from_int, i, j);
    phys_lower += phys_from_int_i_j * int_lower[tj];
  }

  // Calculate physical upper bound from physical lower_bound
  const double phys_from_int_i_i = gsl_matrix_get(tiling->phys_from_int, i, i);
  double phys_upper = phys_lower + phys_from_int_i_i * (trie->int_upper - trie->int_lower);

  // Print information on the current trie trie dimension
  fprintf(file, "dim: #%zu/%zu   int: [%+5" LAL_INT8_FORMAT ",%+5" LAL_INT8_FORMAT "]   phys: [%+10g,%+10g]   seqidx:%" LAL_UINT8_FORMAT "\n",
          ti + 1, tn, trie->int_lower, trie->int_upper, phys_lower, phys_upper, trie->seqidx);

  // If this is not the highest dimension, loop over this dimension
  LT_IndexTrie *next = trie->next;
  if (next != NULL) {
    for (int32_t point = trie->int_lower; point <= trie->int_upper; ++point, ++next) {

      // Set 'i'th integer lower bound to this point
      int_lower[ti] = point;

      // Print higher dimensions
      LT_PrintIndexTrie(tiling, next, ti + 1, file, int_lower);

    }
  }

}

LatticeTiling *XLALCreateLatticeTiling(
  const size_t ndim
  )
{

  // Check input
  XLAL_CHECK_NULL(ndim > 0, XLAL_EINVAL);

  // Allocate memory
  LatticeTiling *tiling = XLALCalloc(1, sizeof(*tiling));
  XLAL_CHECK_NULL(tiling != NULL, XLAL_ENOMEM);
  tiling->bounds = XLALCalloc(ndim, sizeof(*tiling->bounds));
  XLAL_CHECK_NULL(tiling->bounds != NULL, XLAL_ENOMEM);

  // Initialise fields
  tiling->ndim = ndim;
  tiling->lattice = TILING_LATTICE_MAX;

  return tiling;

}

void XLALDestroyLatticeTiling(
  LatticeTiling *tiling
  )
{
  if (tiling != NULL) {
    if (tiling->bounds != NULL) {
      for (size_t i = 0; i < tiling->ndim; ++i) {
        XLALFree(tiling->bounds[i].data_lower);
        XLALFree(tiling->bounds[i].data_upper);
      }
      XLALFree(tiling->bounds);
    }
    XLALFree(tiling->tiled_idx);
    XLALFree(tiling->stats);
    GFMAT(tiling->int_from_phys, tiling->phys_from_int, tiling->tiled_generator);
    GFVEC(tiling->phys_bbox, tiling->phys_origin);
    XLALFree(tiling);
  }
}

int XLALSetLatticeTilingBound(
  LatticeTiling *tiling,
  const size_t dim,
  const LatticeTilingBound func,
  const size_t data_len,
  void *data_lower,
  void *data_upper
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->lattice == TILING_LATTICE_MAX, XLAL_EINVAL);
  XLAL_CHECK(dim < tiling->ndim, XLAL_ESIZE);
  XLAL_CHECK(func != NULL, XLAL_EFAULT);
  XLAL_CHECK(data_len > 0, XLAL_EFAULT);
  XLAL_CHECK(data_lower != NULL, XLAL_EFAULT);
  XLAL_CHECK(data_upper != NULL, XLAL_EFAULT);

  // Check that bound has not already been set
  XLAL_CHECK(tiling->bounds[dim].func == NULL, XLAL_EINVAL, "Lattice tiling dimension #%zu is already bounded", dim);

  // Determine if this dimension is tiled
  const bool is_tiled = (memcmp(data_lower, data_upper, data_len) != 0);

  // Set the parameter-space bound
  tiling->bounds[dim].is_tiled = is_tiled;
  tiling->bounds[dim].func = func;
  tiling->bounds[dim].data_len = data_len;
  tiling->bounds[dim].data_lower = data_lower;
  tiling->bounds[dim].data_upper = data_upper;
  tiling->bounds[dim].pad_lower = 1;
  tiling->bounds[dim].pad_upper = 1;

  return XLAL_SUCCESS;

}

int XLALSetLatticeTilingBoundPadding(
  LatticeTiling *tiling,
  const size_t dim,
  UINT4 pad_lower,
  UINT4 pad_upper
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->lattice == TILING_LATTICE_MAX, XLAL_EINVAL);
  XLAL_CHECK(dim < tiling->ndim, XLAL_ESIZE);

  // Check that bound has been set and is tiled
  XLAL_CHECK(tiling->bounds[dim].func != NULL, XLAL_EINVAL, "Lattice tiling dimension #%zu has not been bounded", dim);
  XLAL_CHECK(tiling->bounds[dim].is_tiled, XLAL_EINVAL, "Lattice tiling dimension #%zu is not tiled, so has no padding", dim);

  // Set the parameter-space bound padding
  tiling->bounds[dim].pad_lower = pad_lower;
  tiling->bounds[dim].pad_upper = pad_upper;

  return XLAL_SUCCESS;

}

static double ConstantBound(
  const void *data,
  const size_t dim UNUSED,
  const gsl_vector *point UNUSED
  )
{

  // Return bound
  return * ((const double *) data);

}

int XLALSetLatticeTilingConstantBound(
  LatticeTiling *tiling,
  const size_t dim,
  const double bound1,
  const double bound2
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(isfinite(bound1), XLAL_EINVAL);
  XLAL_CHECK(isfinite(bound2), XLAL_EINVAL);

  // Allocate memory
  const size_t data_len = sizeof(double);
  double *data_lower = XLALMalloc(data_len);
  XLAL_CHECK(data_lower != NULL, XLAL_ENOMEM);
  double *data_upper = XLALMalloc(data_len);
  XLAL_CHECK(data_upper != NULL, XLAL_ENOMEM);

  // Set the parameter-space bound
  *data_lower = GSL_MIN(bound1, bound2);
  *data_upper = GSL_MAX(bound1, bound2);
  XLAL_CHECK(XLALSetLatticeTilingBound(tiling, dim, ConstantBound, data_len, data_lower, data_upper) == XLAL_SUCCESS, XLAL_EFUNC);

  return XLAL_SUCCESS;

}

int XLALSetTilingLatticeAndMetric(
  LatticeTiling *tiling,
  const TilingLattice lattice,
  const gsl_matrix *metric,
  const double max_mismatch
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->lattice == TILING_LATTICE_MAX, XLAL_EINVAL);
  XLAL_CHECK(lattice < TILING_LATTICE_MAX, XLAL_EINVAL);
  XLAL_CHECK(metric != NULL, XLAL_EFAULT);
  XLAL_CHECK(metric->size1 == tiling->ndim && metric->size2 == tiling->ndim, XLAL_EINVAL);
  XLAL_CHECK(max_mismatch > 0, XLAL_EINVAL);

  const size_t n = tiling->ndim;

  // Check that all parameter-space dimensions are bounded
  for (size_t i = 0; i < n; ++i) {
    XLAL_CHECK(tiling->bounds[i].func != NULL, XLAL_EFAILED, "Lattice tiling dimension #%zu is unbounded", i);
  }

  // Check metric is symmetric and has positive diagonal elements
  for (size_t i = 0; i < n; ++i) {
    XLAL_CHECK(gsl_matrix_get(metric, i, i) > 0, XLAL_EINVAL, "Parameter-space metric(%zu,%zu) <= 0", i, i);
    for (size_t j = i + 1; j < n; ++j) {
      XLAL_CHECK(gsl_matrix_get(metric, i, j) == gsl_matrix_get(metric, j, i), XLAL_EINVAL, "Parameter-space metric(%zu,%zu) != metric(%zu,%zu)", i, j, j, i);
    }
  }

  // Set fields
  tiling->lattice = lattice;

  // Count number of tiled dimensions
  tiling->tiled_ndim = 0;
  for (size_t i = 0; i < n; ++i) {
    if (tiling->bounds[i].is_tiled) {
      ++tiling->tiled_ndim;
    }
  }

  // Allocate and initialise vectors and matrices
  GAVEC(tiling->phys_bbox, n);
  GAVEC(tiling->phys_origin, n);
  GAMAT(tiling->int_from_phys, n, n);
  gsl_matrix_set_identity(tiling->int_from_phys);
  GAMAT(tiling->phys_from_int, n, n);
  gsl_matrix_set_identity(tiling->phys_from_int);

  // If no parameter-space dimensions are tiled, we're done
  if (tiling->tiled_ndim == 0) {
    return XLAL_SUCCESS;
  }

  const size_t tn = tiling->tiled_ndim;

  // Build index to tiled parameter-space dimensions
  tiling->tiled_idx = XLALCalloc(tn, sizeof(*tiling->tiled_idx));
  XLAL_CHECK(tiling->tiled_idx != NULL, XLAL_ENOMEM);
  for (size_t i = 0, ti = 0; i < n; ++i) {
    if (tiling->bounds[i].is_tiled) {
      tiling->tiled_idx[ti++] = i;
    }
  }

  // Calculate normalisation scale from metric diagonal elements
  gsl_vector *GAVEC(t_norm, tn);
  for (size_t ti = 0; ti < tn; ++ti) {
    const size_t i = tiling->tiled_idx[ti];
    const double metric_i_i = gsl_matrix_get(metric, i, i);
    gsl_vector_set(t_norm, ti, sqrt(metric_i_i));
  }

  // Copy and normalise tiled dimensions of metric
  gsl_matrix *GAMAT(t_metric, tn, tn);
  for (size_t ti = 0; ti < tn; ++ti) {
    const size_t i = tiling->tiled_idx[ti];
    const double t_norm_ti = gsl_vector_get(t_norm, ti);
    for (size_t tj = 0; tj < tn; ++tj) {
      const size_t j = tiling->tiled_idx[tj];
      const double t_norm_tj = gsl_vector_get(t_norm, tj);
      gsl_matrix_set(t_metric, ti, tj, gsl_matrix_get(metric, i, j) / t_norm_ti / t_norm_tj);
    }
  }

  // Compute metric ellipse bounding box
  gsl_vector *t_bbox = XLALMetricEllipseBoundingBox(t_metric, max_mismatch);
  XLAL_CHECK(t_bbox != NULL, XLAL_EFUNC);

  // Copy bounding box in physical coordinates to tiled dimensions
  for (size_t ti = 0; ti < tn; ++ti) {
    const size_t i = tiling->tiled_idx[ti];
    const double t_norm_ti = gsl_vector_get(t_norm, ti);
    gsl_vector_set(tiling->phys_bbox, i, gsl_vector_get(t_bbox, ti) / t_norm_ti);
  }

  // Set physical parameter-space origin to mid-point of parameter-space bounds
  for (size_t i = 0; i < n; ++i) {
    double phys_lower = 0.0, phys_upper = 0.0;
    LT_GetBounds(tiling, i, tiling->phys_origin, &phys_lower, &phys_upper);
    gsl_vector_set(tiling->phys_origin, i, 0.5* (phys_lower + phys_upper));
  }

  // Set non-tiled dimensions of physical parameter-space origin back to zero
  for (size_t i = 0; i < n; ++i) {
    if (!tiling->bounds[i].is_tiled) {
      gsl_vector_set(tiling->phys_origin, i, 0);
    }
  }

  // Compute a lower-triangular basis matrix whose columns are orthonormal with respect to the tiled metric
  gsl_matrix *GAMAT(t_basis, tn, tn);
  {
    // We want to find a lower-triangular basis such that:
    //   basis^T * metric * basis = I
    // This is rearranged to give:
    //   metric^-1 = basis * basis^T
    // Hence basis is the Cholesky decomposition of metric^-1
    gsl_matrix_memcpy(t_basis, t_metric);
    GCALL(gsl_linalg_cholesky_decomp(t_basis), "Parameter-space metric is not positive definite");
    GCALL(gsl_linalg_cholesky_invert(t_basis), "Parameter-space metric cannot be inverted");
    GCALL(gsl_linalg_cholesky_decomp(t_basis), "Inverse of parameter-space metric is not positive definite");

    // gsl_linalg_cholesky_decomp() stores both basis and basis^T
    // in the same matrix; zero out upper triangle to get basis
    LT_ZeroStrictUpperTriangle(t_basis);
  }

  // Compute a lower-triangular generator matrix for a given lattice type and mismatch
  GAMAT(tiling->tiled_generator, tn, tn);
  {

    // Compute lattice generator and normalised thickness
    double norm_thickness = 0.0;
    switch (lattice) {

    case TILING_LATTICE_CUBIC: {	// Cubic (\f$Z_n\f$) lattice

      // Zn lattice generator is the identity
      gsl_matrix_set_identity(tiling->tiled_generator);

      // Zn normalised thickness
      norm_thickness = pow(sqrt(tn)/2.0, tn);

    }
      break;

    case TILING_LATTICE_ANSTAR: {	// An-star (\f$A_n^*\f$) lattice

      // An* lattice generator in tn+1 dimensions, given in:
      //   McKilliam et.al., "A linear-time nearest point algorithm for the lattice An*"
      //   in "International Symposium on Information Theory and Its Applications", ISITA2008,
      //   Auckland, New Zealand, 7-10 Dec. 2008. DOI: 10.1109/ISITA.2008.4895596
      gsl_matrix *GAMAT(G, tn + 1, tn + 1);
      gsl_matrix_set_identity(G);
      gsl_matrix_add_constant(G, -1.0 / (tn + 1));

      // Find the QL decomposition of the generator matrix G, excluding 1st column,
      // which is linearly dependent on the remaining columns:
      //   G(:, 2:end) = Gp = Q * L
      // where Q is an orthogonal matrix and L an lower-triangular matrix.
      // This is found using the more commonly implemented QR decomposition by:
      // - reversing the order of the rows/columns of Gp
      // - decomposing Gp = Qp * Lp, where Lp is upper triangular
      // - reversing the order of the rows/columns of Qp to give Q
      // - reversing the order of the rows/columns of Lp to give L
      gsl_matrix_view Gp = gsl_matrix_submatrix(G, 0, 1, tn + 1, tn);
      LT_ReverseOrderRowsCols(&Gp.matrix);
      gsl_vector *GAVEC(tau, tn);
      GCALL(gsl_linalg_QR_decomp(&Gp.matrix, tau), "'G' cannot be QR-decomposed");
      gsl_matrix *GAMAT(Q, tn + 1, tn + 1);
      gsl_matrix *GAMAT(L, tn + 1, tn);
      gsl_linalg_QR_unpack(&Gp.matrix, tau, Q, L);
      LT_ReverseOrderRowsCols(Q);
      LT_ReverseOrderRowsCols(L);

      // Discard the first row of L, which is zero, to get the generator in tn dimensions
      gsl_matrix_view L_view = gsl_matrix_submatrix(L, 1, 0, tn, tn);
      gsl_matrix_memcpy(tiling->tiled_generator, &L_view.matrix);

      // Cleanup
      GFMAT(G, L, Q);
      GFVEC(tau);

      // An* normalised thickness
      norm_thickness = sqrt(tn+1.0) * pow((1.0*tn* (tn+2.0)) / (12.0* (tn+1.0)), 0.5*tn);

    }
      break;

    default:
      XLAL_ERROR(XLAL_EINVAL, "Invalid lattice=%u", lattice);
    }

    // Generator will be lower-triangular, so zero out upper triangle
    LT_ZeroStrictUpperTriangle(tiling->tiled_generator);

    // Ensure that the generator has positive diagonal elements, by
    // changing the sign of the columns of the matrix, if necessary
    for (size_t tj = 0; tj < tn; ++tj) {
      gsl_vector_view generator_col = gsl_matrix_column(tiling->tiled_generator, tj);
      XLAL_CHECK(gsl_vector_get(&generator_col.vector, tj) != 0, XLAL_ERANGE, "Generator matrix(%zu,%zu) == 0", tj, tj);
      if (gsl_vector_get(&generator_col.vector, tj) < 0) {
        gsl_vector_scale(&generator_col.vector, -1);
      }
    }

    // Compute generator LU decomposition
    gsl_matrix *GAMAT(LU_decomp, tn, tn);
    gsl_matrix_memcpy(LU_decomp, tiling->tiled_generator);
    gsl_permutation *GAPERM(LU_perm, tn);
    int LU_sign = 0;
    GCALL(gsl_linalg_LU_decomp(LU_decomp, LU_perm, &LU_sign), "Generator matrix cannot be LU-decomposed");

    // Compute generator determinant
    const double generator_determinant = XLALMetricDeterminant(tiling->tiled_generator);
    XLAL_CHECK(!XLAL_IS_REAL8_FAIL_NAN(generator_determinant), XLAL_EFUNC);

    // Compute generator covering radius
    const double generator_covering_radius = pow(norm_thickness * generator_determinant, 1.0 / tn);

    // Normalise so covering spheres have sqrt(max_mismatch) covering radii
    gsl_matrix_scale(tiling->tiled_generator, sqrt(max_mismatch) / generator_covering_radius);

    // Cleanup
    GFMAT(LU_decomp);
    GFPERM(LU_perm);

  }

  // Compute transform to normalised physical coordinates from generating integers
  gsl_matrix *GAMAT(t_norm_from_int, tn, tn);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, t_basis, tiling->tiled_generator, 0.0, t_norm_from_int);
  LT_ZeroStrictUpperTriangle(t_norm_from_int);

  // Compute transform to generating integers from normalised physical coordinates
  gsl_matrix *GAMAT(t_int_from_norm, tn, tn);
  gsl_matrix_set_identity(t_int_from_norm);
  gsl_blas_dtrsm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, t_norm_from_int, t_int_from_norm);
  LT_ZeroStrictUpperTriangle(t_int_from_norm);

  // Set tiled dimensions of transforms, and convert to unnormalised physical coordinates
  for (size_t ti = 0; ti < tn; ++ti) {
    const size_t i = tiling->tiled_idx[ti];
    const double t_norm_ti = gsl_vector_get(t_norm, ti);
    for (size_t tj = 0; tj < tn; ++tj) {
      const size_t j = tiling->tiled_idx[tj];
      const double t_norm_tj = gsl_vector_get(t_norm, tj);
      gsl_matrix_set(tiling->int_from_phys, i, j, gsl_matrix_get(t_int_from_norm, ti, tj) * t_norm_tj);
      gsl_matrix_set(tiling->phys_from_int, i, j, gsl_matrix_get(t_norm_from_int, ti, tj) / t_norm_ti);
    }
  }

  // Round tiled dimensions of physical parameter-space origin to nearest lattice step size, then
  // shift by half a step size. This ensures that the tiling will never place a lattice point at zero
  // in physical coordinates, since the physical coordinates may not be well-defined at zero.
  for (size_t ti = 0; ti < tn; ++ti) {
    const size_t i = tiling->tiled_idx[ti];
    const double int_from_phys_i_i = gsl_matrix_get(tiling->int_from_phys, i, i);
    const double phys_from_int_i_i = gsl_matrix_get(tiling->phys_from_int, i, i);
    double phys_origin_i = gsl_vector_get(tiling->phys_origin, i);
    phys_origin_i = (round(phys_origin_i * int_from_phys_i_i) + 0.5) * phys_from_int_i_i;
    gsl_vector_set(tiling->phys_origin, i, phys_origin_i);
  }

  // Cleanup
  GFMAT(t_metric, t_basis, t_norm_from_int, t_int_from_norm);
  GFVEC(t_norm, t_bbox);

  return XLAL_SUCCESS;

}

size_t XLALTotalLatticeTilingDimensions(
  const LatticeTiling *tiling
  )
{

  // Check input
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);

  return tiling->ndim;

}

size_t XLALTiledLatticeTilingDimensions(
  const LatticeTiling *tiling
  )
{

  // Check input
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_VAL(0, tiling->lattice < TILING_LATTICE_MAX, XLAL_EINVAL);

  return tiling->tiled_ndim;

}

REAL8 XLALLatticeTilingStepSizes(
  const LatticeTiling *tiling,
  const size_t dim
  )
{

  // Check input
  XLAL_CHECK_REAL8(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_REAL8(tiling->lattice < TILING_LATTICE_MAX, XLAL_EINVAL);
  XLAL_CHECK_REAL8(dim < tiling->ndim, XLAL_ESIZE);

  // Return 0 for non-tiled dimensions
  if (!tiling->bounds[dim].is_tiled) {
    return 0.0;
  }

  // Step size is the (dim)th diagonal element of 'phys_from_int'
  return gsl_matrix_get(tiling->phys_from_int, dim, dim);

}

const LatticeTilingStats *XLALLatticeTilingStatistics(
  LatticeTiling *tiling,
  const size_t dim
  )
{

  // Check input
  XLAL_CHECK_NULL(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(tiling->lattice < TILING_LATTICE_MAX, XLAL_EINVAL);
  XLAL_CHECK_NULL(dim < tiling->ndim, XLAL_ESIZE);

  // Ensure lattice tiling statistics have been computed
  if (tiling->stats == NULL) {

    const size_t n = tiling->ndim;
    const size_t tn = tiling->tiled_ndim;

    // Allocate memory
    tiling->stats = XLALCalloc(n, sizeof(*tiling->stats));
    XLAL_CHECK_NULL(tiling->stats != NULL, XLAL_ENOMEM);

    // Create iterator over tiling
    LatticeTilingIterator *itr = XLALCreateLatticeTilingIterator(tiling, n, TILING_ITR_DEFAULT);
    XLAL_CHECK_NULL(itr != NULL, XLAL_EFUNC);

    // Start iterator at first point
    XLAL_CHECK_NULL(XLALNextLatticeTilingPoint(itr, NULL) > 0, XLAL_EFUNC);

    // Initialise statistics
    for (size_t i = 0; i < n; ++i) {
      const double phys_point = gsl_vector_get(itr->phys_point, i);
      tiling->stats[i].total_points = 1;
      tiling->stats[i].min_points_pass = 1;
      tiling->stats[i].avg_points_pass = 1;
      tiling->stats[i].max_points_pass = 1;
      tiling->stats[i].min_value_pass = phys_point;
      tiling->stats[i].max_value_pass = phys_point;
    }

    // If lattice tiling contains more than a single point
    if (tn > 0) {

      // Allocate and initialise arrays for computing statistics
      UINT8 t_num_points[tn];
      INT8 t_num_passes[tn], t_min_points[tn], t_max_points[tn];
      double t_min_value[tn], t_max_value[tn];
      for (size_t tj = 0; tj < tn; ++tj) {
        const size_t j = itr->tiling->tiled_idx[tj];
        const INT8 int_lower = itr->int_lower[tj];
        const INT8 int_upper = itr->int_upper[tj];
        const INT8 num_points = int_upper - int_lower + 1;
        t_num_points[tj] = 1;
        t_num_passes[tj] = 1;
        t_min_points[tj] = num_points;
        t_max_points[tj] = num_points;
        t_min_value[tj] = tiling->stats[j].min_value_pass;
        t_max_value[tj] = tiling->stats[j].max_value_pass;
      }

      // Iterate over remaining points; XLALNextLatticeTilingPoint() returns the index
      // (offset from 1) of the lowest dimension where the current point has changed
      xlalErrno = 0;
      int ti_plus_1;
      while ((ti_plus_1 = XLALNextLatticeTilingPoint(itr, NULL)) > 0) {
        const size_t ti = ti_plus_1 - 1;

        // Compute statistics for each dimension which has changed
        t_num_points[ti] += 1;
        for (size_t tj = ti + 1; tj < tn; ++tj) {
          const size_t j = itr->tiling->tiled_idx[tj];
          const INT8 int_lower = itr->int_lower[tj];
          const INT8 int_upper = itr->int_upper[tj];
          const INT8 num_points = int_upper - int_lower + 1;
          const double phys_point = gsl_vector_get(itr->phys_point, j);
          t_num_points[tj] += 1;
          t_num_passes[tj] += 1;
          t_min_points[tj] = GSL_MIN(t_min_points[tj], num_points);
          t_max_points[tj] = GSL_MAX(t_max_points[tj], num_points);
          t_min_value[tj] = GSL_MIN(t_min_value[tj], phys_point);
          t_max_value[tj] = GSL_MAX(t_max_value[tj], phys_point);
        }

        // Fast-forward iterator over highest tiled dimension
        t_num_points[tn - 1] += LT_FastForwardIterator(itr);
        {
          const size_t j = itr->tiling->tiled_idx[tn - 1];
          const double phys_point = gsl_vector_get(itr->phys_point, j);
          t_min_value[tn - 1] = GSL_MIN(t_min_value[tn - 1], phys_point);
          t_max_value[tn - 1] = GSL_MAX(t_max_value[tn - 1], phys_point);
        }

      }

      // Store statistics
      for (size_t tj = 0; tj < tn; ++tj) {
        const size_t j = itr->tiling->tiled_idx[tj];
        for (size_t k = j; k < n; ++k) {
          // Non-tiled dimensions should inherit their total number of points from lower dimensions
          tiling->stats[k].total_points = t_num_points[tj];
        }
        tiling->stats[j].min_points_pass = t_min_points[tj];
        tiling->stats[j].avg_points_pass = ((double) t_num_points[tj]) / t_num_passes[tj];
        tiling->stats[j].max_points_pass = t_max_points[tj];
        tiling->stats[j].min_value_pass = t_min_value[tj];
        tiling->stats[j].max_value_pass = t_max_value[tj];
      }

    }

    // Cleanup
    XLALDestroyLatticeTilingIterator(itr);

  }

  return &tiling->stats[dim];

}

int XLALRandomLatticeTilingPoints(
  const LatticeTiling *tiling,
  const double scale,
  RandomParams *rng,
  gsl_matrix *random_points
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->lattice < TILING_LATTICE_MAX, XLAL_EINVAL);
  XLAL_CHECK(scale > -1.0, XLAL_EINVAL);
  XLAL_CHECK(rng != NULL, XLAL_EFAULT);
  XLAL_CHECK(random_points != NULL, XLAL_EFAULT);
  XLAL_CHECK(random_points->size1 == tiling->ndim, XLAL_ESIZE);

  const size_t n = tiling->ndim;

  // Generate random points in parameter space
  for (size_t k = 0; k < random_points->size2; ++k) {
    gsl_vector_view phys_point = gsl_matrix_column(random_points, k);
    for (size_t i = 0; i < n; ++i) {

      // Get the physical bounds on the current dimension
      double phys_lower = 0.0, phys_upper = 0.0;
      LT_GetBounds(tiling, i, &phys_point.vector, &phys_lower, &phys_upper);

      // Generate random number
      const double u = (1.0 + scale) * (XLALUniformDeviate(rng) - 0.5) + 0.5;

      // Set parameter-space point
      gsl_vector_set(&phys_point.vector, i, phys_lower + u* (phys_upper - phys_lower));

    }

  }

  return XLAL_SUCCESS;

}

LatticeTilingIterator *XLALCreateLatticeTilingIterator(
  LatticeTiling *tiling,
  const size_t itr_ndim,
  const TilingIteratorFlags flags
  )
{

  // Check input
  XLAL_CHECK_NULL(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(tiling->lattice < TILING_LATTICE_MAX, XLAL_EINVAL);
  XLAL_CHECK_NULL(itr_ndim <= tiling->ndim, XLAL_EINVAL);
  XLAL_CHECK_NULL(flags < TILING_ITR_MAX, XLAL_EINVAL);

  // Allocate memory
  LatticeTilingIterator *itr = XLALCalloc(1, sizeof(*itr));
  XLAL_CHECK_NULL(itr != NULL, XLAL_ENOMEM);

  // Store reference to lattice tiling
  itr->tiling = tiling;

  // Set fields
  itr->itr_ndim = itr_ndim;
  itr->flags = flags;

  // Determine the maximum tiled dimension to iterate over
  itr->tiled_itr_ndim = 0;
  for (size_t i = 0; i < itr_ndim; ++i) {
    if (itr->tiling->bounds[i].is_tiled) {
      ++itr->tiled_itr_ndim;
    }
  }

  const size_t n = itr->tiling->ndim;
  const size_t tn = itr->tiling->tiled_ndim;

  // Allocate and initialise vectors and matrices
  GAVEC_NULL(itr->phys_point, n);
  if (tn > 0) {
    itr->int_lower = XLALCalloc(tn, sizeof(*itr->int_lower));
    XLAL_CHECK_NULL(itr->int_lower != NULL, XLAL_EINVAL);
    itr->int_point = XLALCalloc(tn, sizeof(*itr->int_point));
    XLAL_CHECK_NULL(itr->int_point != NULL, XLAL_EINVAL);
    itr->int_upper = XLALCalloc(tn, sizeof(*itr->int_upper));
    XLAL_CHECK_NULL(itr->int_upper != NULL, XLAL_EINVAL);
    itr->direction = XLALCalloc(tn, sizeof(*itr->direction));
    XLAL_CHECK_NULL(itr->direction != NULL, XLAL_EINVAL);
  }

  // Set iterator to beginning of lattice tiling
  XLAL_CHECK_NULL(XLALResetLatticeTilingIterator(itr) == XLAL_SUCCESS, XLAL_EFUNC);

  return itr;

}

void XLALDestroyLatticeTilingIterator(
  LatticeTilingIterator *itr
  )
{
  if (itr) {
    GFVEC(itr->phys_point);
    XLALFree(itr->int_lower);
    XLALFree(itr->int_point);
    XLALFree(itr->int_upper);
    XLALFree(itr->direction);
    XLALFree(itr);
  }
}

int XLALResetLatticeTilingIterator(
  LatticeTilingIterator *itr
  )
{

  // Check input
  XLAL_CHECK(itr != NULL, XLAL_EFAULT);

  // Return iterator to initialised state
  itr->state = 0;

  return XLAL_SUCCESS;

}

int XLALNextLatticeTilingPoint(
  LatticeTilingIterator *itr,
  gsl_vector *point
  )
{

  // Check input
  XLAL_CHECK(itr != NULL, XLAL_EFAULT);
  XLAL_CHECK(point == NULL || point->size == itr->tiling->ndim, XLAL_EINVAL);

  const size_t n = itr->tiling->ndim;
  const size_t tn = itr->tiling->tiled_ndim;

  // If iterator is finished, we're done
  if (itr->state > 1) {
    return 0;
  }

  // Which dimensions have changed?
  size_t changed_ti;

  // Which dimensions need to be reset?
  size_t reset_ti;

  if (itr->state == 0) {	// Iterator has been initialised

    // Initialise lattice point
    gsl_vector_set_zero(itr->phys_point);
    for (size_t ti = 0; ti < tn; ++ti) {
      itr->int_point[ti] = 0;
    }

    // Initialise iteration order to 1 for negative order.
    for (size_t ti = 0; ti < tn; ++ti) {
      itr->direction[ti] = 1;
    }

    // Initialise index
    itr->index = 0;

    // All dimensions have changed
    changed_ti = 0;

    // All dimensions need to be reset
    reset_ti = 0;

  } else {			// Iterator is in progress

    // Start iterating from the maximum tiled dimension specified at iterator creation
    size_t ti = itr->tiled_itr_ndim;

    // Find the next lattice point
    while (true) {

      // If dimension index is now zero, we're done
      if (ti == 0) {

        // Iterator is now finished
        itr->state = 2;

        return 0;

      }

      // Decrement current dimension index
      --ti;

      // Increment integer point in this dimension, in current direction
      const INT8 direction = itr->direction[ti];
      const INT8 int_point_ti = itr->int_point[ti] + direction;
      itr->int_point[ti] = int_point_ti;

      // Increment physical point in this dimension, in current direction
      const size_t i = itr->tiling->tiled_idx[ti];
      gsl_vector_const_view phys_from_int_i = gsl_matrix_const_column(itr->tiling->phys_from_int, i);
      gsl_blas_daxpy(direction, &phys_from_int_i.vector, itr->phys_point);

      // If point is not out of bounds, we have found the next lattice point
      const INT8 int_lower_ti = itr->int_lower[ti];
      const INT8 int_upper_ti = itr->int_upper[ti];
      if ((direction > 0 && int_point_ti <= int_upper_ti) || (direction < 0 && int_point_ti >= int_lower_ti)) {
        break;
      }

      // Move on to lower dimensions
      continue;

    }

    // Point was found, so increase index
    ++itr->index;

    // This dimension and higher have changed
    changed_ti = ti;

    // Higher dimensions need to be reset
    reset_ti = ti + 1;

  }

  // Reset specified dimensions
  for (size_t i = 0, ti = 0; i < n; ++i) {

    // Get extremal physical bounds
    double phys_lower = 0, phys_upper = 0;
    LT_GetExtremalBounds(itr->tiling, i, itr->phys_point, &phys_lower, &phys_upper);

    // If not tiled, set current physical point to non-tiled parameter-space bound
    if (!itr->tiling->bounds[i].is_tiled) {
      gsl_vector_set(itr->phys_point, i, phys_lower);
      continue;
    }

    // If tiled dimension needs to be reset:
    if (ti >= reset_ti) {

      // Transform physical point in lower dimensions to generating integer offset
      const double phys_origin_i = gsl_vector_get(itr->tiling->phys_origin, i);
      double int_from_phys_point_i = 0;
      for (size_t j = 0; j < i; ++j) {
        const double int_from_phys_i_j = gsl_matrix_get(itr->tiling->int_from_phys, i, j);
        const double phys_point_j = gsl_vector_get(itr->phys_point, j);
        const double phys_origin_j = gsl_vector_get(itr->tiling->phys_origin, j);
        int_from_phys_point_i += int_from_phys_i_j * (phys_point_j - phys_origin_j);
      }

      {
        // Transform physical bounds to generating integers
        const double int_from_phys_i_i = gsl_matrix_get(itr->tiling->int_from_phys, i, i);
        const double dbl_int_lower_i = int_from_phys_point_i + int_from_phys_i_i * (phys_lower - phys_origin_i);
        const double dbl_int_upper_i = int_from_phys_point_i + int_from_phys_i_i * (phys_upper - phys_origin_i);

        // Compute integer lower/upper bounds, rounded up/down to avoid extra boundary points
        feclearexcept(FE_ALL_EXCEPT);
        const INT8 int_lower_i = lround(ceil(dbl_int_lower_i));
        const INT8 int_upper_i = lround(floor(dbl_int_upper_i));
        XLAL_CHECK(fetestexcept(FE_INVALID) == 0, XLAL_EFAILED, "Integer bounds on dimension #%zu are too large: %0.2e to %0.2e", i, dbl_int_lower_i, dbl_int_upper_i);

        // Set integer lower/upper bounds
        itr->int_lower[ti] = int_lower_i;
        itr->int_upper[ti] = int_upper_i;
      }
      const INT8 int_lower_i = itr->int_lower[ti];
      const INT8 int_upper_i = itr->int_upper[ti];

      // Switch iterator direction, depending on given order
      INT8 direction = itr->direction[ti];
      if (itr->flags & TILING_ITR_ALT_ORDER) {

        // Only switch direction:
        // - if iterator is in progress
        // - for iterated-over dimensions
        // - if there is more than one point in this dimension
        if ((itr->state > 0) && (ti < itr->tiled_itr_ndim) && (int_lower_i < int_upper_i)) {
          direction = -direction;
          itr->direction[ti] = direction;
        }

      }

      // Set integer point to lower or bound, depending on current direction
      itr->int_point[ti] = (direction > 0) ? int_lower_i : int_upper_i;

      // Set current physical point from integer point
      double phys_point_i = phys_origin_i;
      for (size_t tj = 0; tj < tn; ++tj) {
        const size_t j = itr->tiling->tiled_idx[tj];
        const double phys_from_int_i_j = gsl_matrix_get(itr->tiling->phys_from_int, i, j);
        const INT8 int_point_tj = itr->int_point[tj];
        phys_point_i += phys_from_int_i_j * int_point_tj;
      }
      gsl_vector_set(itr->phys_point, i, phys_point_i);

    }

    ++ti;

  }

  // Iterator is in progress
  itr->state = 1;

  // Optionally, copy current physical point
  if (point != NULL) {
    gsl_vector_memcpy(point, itr->phys_point);
  }

  // Return index of changed dimensions (offset from 1, since 0 is used to indicate no more points)
  return 1 + changed_ti;

}

int XLALNextLatticeTilingPoints(
  LatticeTilingIterator *itr,
  gsl_matrix **points
  )
{

  // Check input
  XLAL_CHECK(itr != NULL, XLAL_EFAULT);
  XLAL_CHECK(points != NULL && *points != NULL, XLAL_EFAULT);
  XLAL_CHECK((*points)->size1 == itr->tiling->ndim, XLAL_EINVAL);

  // Fill 'points' with points from XLALNextLatticeTilingPoint(), but stop if there are none left
  size_t j = 0;
  for (; j < (*points)->size2; ++j) {
    gsl_vector_view point_j = gsl_matrix_column(*points, j);
    int retn = XLALNextLatticeTilingPoint(itr, &point_j.vector);
    XLAL_CHECK(retn >= 0, XLAL_EFUNC, "XLALNextLatticeTilingPoint() failed at j=%zu", j);
    if (retn == 0) {
      break;
    }
  }

  // If there are fewer points than the size of 'points', resize 'points' to fit
  if (0 < j && j < (*points)->size2) {
    gsl_matrix *GAMAT(new_points, (*points)->size1, j);
    gsl_matrix_view points_view = gsl_matrix_submatrix(*points, 0, 0, (*points)->size1, j);
    gsl_matrix_memcpy(new_points, &points_view.matrix);
    GFMAT(*points);
    *points = new_points;
  }

  return j;

}

UINT8 XLALLatticeTilingPointsInPass(
  const LatticeTilingIterator *itr,
  const size_t dim
  )
{

  // Check input
  XLAL_CHECK_VAL(0, itr != NULL, XLAL_EFAULT);
  XLAL_CHECK_VAL(0, itr->state > 0, XLAL_EINVAL);
  XLAL_CHECK_VAL(0, dim < itr->tiling->ndim, XLAL_EINVAL);

  // If this is a tiled dimension, return the number of points in this pass
  for (size_t ti = 0; ti < itr->tiling->tiled_ndim; ++ti) {

    // If this tiled dimension matches 'dim'
    const size_t i = itr->tiling->tiled_idx[ti];
    if (i == dim) {

      // Return number of points in current pass
      return itr->int_upper[ti] - itr->int_lower[ti] + 1;

    }

  }

  // This is either a non-tiled or a non-iterated dimension
  return 1;

}

UINT8 XLALTotalLatticeTilingPoints(
  LatticeTilingIterator *itr
  )
{

  // Check input
  XLAL_CHECK_VAL(0, itr != NULL, XLAL_EFAULT);

  // Return 1 if iterator or tiling contains only a single point
  if (itr->tiled_itr_ndim == 0 || itr->tiling->tiled_ndim == 0) {
    return 1;
  }

  // Get lattice tiling statistics for highest iterated dimension
  const LatticeTilingStats *stats = XLALLatticeTilingStatistics(itr->tiling, itr->itr_ndim - 1);
  XLAL_CHECK_VAL(0, stats != NULL, XLAL_EFUNC);

  return stats->total_points;

}

UINT8 XLALCurrentLatticeTilingIndex(
  LatticeTilingIterator *itr
  )
{

  // Check input
  XLAL_CHECK_VAL(0, itr != NULL, XLAL_EFAULT);
  XLAL_CHECK_VAL(0, itr->state > 0, XLAL_EINVAL);

  return itr->index;

}

LatticeTilingLocator *XLALCreateLatticeTilingLocator(
  LatticeTiling *tiling
  )
{

  // Check input
  XLAL_CHECK_NULL(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(tiling->lattice < TILING_LATTICE_MAX, XLAL_EINVAL);

  // Allocate memory
  LatticeTilingLocator *loc = XLALCalloc(1, sizeof(*loc));
  XLAL_CHECK_NULL(loc != NULL, XLAL_ENOMEM);

  // Store reference to lattice tiling
  loc->tiling = tiling;

  // Set fields
  loc->ndim = tiling->ndim;
  loc->tiled_ndim = tiling->tiled_ndim;

  // Build index trie to enforce parameter-space bounds
  if (loc->tiled_ndim > 0) {

    // Create iterator over the bounded dimensions
    LatticeTilingIterator *itr = XLALCreateLatticeTilingIterator(tiling, tiling->ndim, TILING_ITR_DEFAULT);

    const size_t n = itr->tiling->ndim;
    const size_t tn = itr->tiling->tiled_ndim;

    // Allocate array of pointers to the next index trie in each dimension
    LT_IndexTrie *next[tn];
    memset(next, 0, sizeof(next));

    // Allocate array containing sequential indices for every dimension
    UINT8 seqidx[n];
    memset(seqidx, 0, sizeof(seqidx));

    // Iterate over all points; XLALNextLatticeTilingPoint() returns the index
    // (offset from 1) of the lowest dimension where the current point has changed
    xlalErrno = 0;
    int ti_plus_1;
    while ((ti_plus_1 = XLALNextLatticeTilingPoint(itr, NULL)) > 0) {
      const size_t ti = ti_plus_1 - 1;

      // Iterate over all dimensions where the current point has changed
      for (size_t tj = ti; tj < tn; ++tj) {

        // If next index trie pointer is NULL, it needs to be initialised
        if (next[tj] == NULL) {

          // Get a pointer to the index trie which needs to be built:
          // - if 'tj' is non-zero, we should use the struct pointed to by 'next' in the lower dimension
          // - otherwise, this is the first point of the tiling, so initialise the base index trie
          LT_IndexTrie *trie = NULL;
          if (tj > 0) {
            trie = next[tj - 1];
          } else {
            trie = loc->index_trie = XLALCalloc(1, sizeof(*trie));
            XLAL_CHECK_NULL(loc->index_trie != NULL, XLAL_ENOMEM);
          }

          // Save the lower and upper integer point bounds
          trie->int_lower = itr->int_lower[tj];
          trie->int_upper = itr->int_upper[tj];

          // Save the sequential index of the current point up to this dimension
          trie->seqidx = seqidx[tj];

          if (tj + 1 < tn) {

            // If we are below the highest dimension, allocate a new
            // array of index tries for the next highest dimension
            const size_t next_length = trie->int_upper - trie->int_lower + 1;
            trie->next = XLALCalloc(next_length, sizeof(*trie->next));
            XLAL_CHECK_NULL(trie->next != NULL, XLAL_ENOMEM);

            // Point 'next[tj]' to this array, for higher dimensions to use
            next[tj] = trie->next;

          }

        } else {

          // Otherwise, advance to the next index trie in the array
          ++next[tj];

        }

        // If we are below the highest dimension, set 'next' in the next highest
        // dimension to NULL, so that on the next loop a new array will be created
        if (tj + 1 < tn) {
          next[tj + 1] = NULL;
        }

      }

      // Increment sequential index in every higher dimension
      for (size_t i = itr->tiling->tiled_idx[ti]; i < n; ++i) {
        ++seqidx[i];
      }

      // Fast-forward iterator over highest tiled dimension.
      seqidx[itr->tiling->tiled_idx[tn - 1]] += LT_FastForwardIterator(itr);

    }
    XLAL_CHECK_NULL(xlalErrno == 0, XLAL_EFUNC);

    // Cleanup
    XLALDestroyLatticeTilingIterator(itr);

  }

  return loc;

}

void XLALDestroyLatticeTilingLocator(
  LatticeTilingLocator *loc
  )
{
  if (loc) {
    if (loc->index_trie != NULL) {
      LT_FreeIndexTrie(loc->index_trie);
      XLALFree(loc->index_trie);
    }
    XLALFree(loc);
  }
}

int XLALNearestLatticeTilingPoints(
  const LatticeTilingLocator *loc,
  const gsl_matrix *points,
  const size_t seqidx_ndim,
  gsl_matrix **nearest_points,
  UINT8VectorSequence **nearest_idxs
  )
{

  // Check input
  XLAL_CHECK(loc != NULL, XLAL_EFAULT);
  XLAL_CHECK(points != NULL, XLAL_EFAULT);
  XLAL_CHECK(points->size1 == loc->ndim, XLAL_ESIZE);
  XLAL_CHECK(seqidx_ndim <= loc->ndim, XLAL_EINVAL);
  XLAL_CHECK(nearest_points != NULL, XLAL_EFAULT);

  const size_t n = loc->ndim;
  const size_t tn = loc->tiled_ndim;
  const size_t npoints = points->size2;

  // Resize or allocate nearest point matrix, if required, and create view of correct size
  if (*nearest_points != NULL) {
    if ((*nearest_points)->size1 != n || (*nearest_points)->size2 < npoints) {
      GFMAT(*nearest_points);
      *nearest_points = NULL;
    }
  }
  if (*nearest_points == NULL) {
    GAMAT(*nearest_points, n, npoints);
  }
  gsl_matrix_view nearest_view = gsl_matrix_submatrix(*nearest_points, 0, 0, n, npoints);
  gsl_matrix *const nearest = &nearest_view.matrix;

  // Resize or allocate nearest index vector sequence, if required, and initialise to zero
  if (nearest_idxs != NULL) {
    if (*nearest_idxs != NULL) {
      if ((*nearest_idxs)->length != npoints || (*nearest_idxs)->vectorLength != n) {
        XLALDestroyUINT8VectorSequence(*nearest_idxs);
        *nearest_idxs = NULL;
      }
    }
    if (*nearest_idxs == NULL) {
      *nearest_idxs = XLALCreateUINT8VectorSequence(npoints, n);
      XLAL_CHECK(*nearest_idxs != NULL, XLAL_ENOMEM);
    }
    memset((*nearest_idxs)->data, 0, npoints * n * sizeof((*nearest_idxs)->data[0]));
  }

  if (tn == 0) {

    // Set all columns of 'nearest' to 'phys_origin', the sole point in the tiling
    for (size_t i = 0; i < n; ++i) {
      const double phys_origin = gsl_vector_get(loc->tiling->phys_origin, i);
      gsl_vector_view nearest_row = gsl_matrix_row(nearest, i);
      gsl_vector_set_all(&nearest_row.vector, phys_origin);
    }

    return XLAL_SUCCESS;

  }

  // Copy 'points' to 'nearest'
  gsl_matrix_memcpy(nearest, points);

  // Subtract physical origin from every point in 'nearest'
  for (size_t i = 0; i < n; ++i) {
    const double phys_origin = gsl_vector_get(loc->tiling->phys_origin, i);
    gsl_vector_view nearest_row = gsl_matrix_row(nearest, i);
    gsl_vector_add_constant(&nearest_row.vector, -phys_origin);
  }

  // Transform 'nearest' from physical coordinates to generating integers
  gsl_blas_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, loc->tiling->int_from_phys, nearest);

  // Find the nearest points in the lattice tiling to the points in 'nearest'
  for (size_t j = 0; j < npoints; ++j) {

    // Find the nearest point to 'nearest[:,j]', the tiled dimensions of which are generating integers
    INT8 nearest_int[tn];
    switch (loc->tiling->lattice) {

    case TILING_LATTICE_CUBIC: {	// Cubic (\f$Z_n\f$) lattice

      // Round each dimension of 'nearest[:,j]' to nearest integer to find the nearest point in Zn
      feclearexcept(FE_ALL_EXCEPT);
      for (size_t ti = 0; ti < tn; ++ti) {
        const size_t i = loc->tiling->tiled_idx[ti];
        nearest_int[ti] = lround(gsl_matrix_get(nearest, i, j));
      }
      if (fetestexcept(FE_INVALID) != 0) {
        XLALPrintError("Rounding failed while finding nearest point #%zu:", j);
        for (size_t ti = 0; ti < tn; ++ti) {
          const size_t i = loc->tiling->tiled_idx[ti];
          XLALPrintError(" %0.2e", gsl_matrix_get(nearest, i, j));
        }
        XLALPrintError("\n");
        XLAL_ERROR(XLAL_EFAILED);
      }

    }
      break;

    case TILING_LATTICE_ANSTAR: {	// An-star (\f$A_n^*\f$) lattice

      // The nearest point algorithm used below embeds the An* lattice in tn+1 dimensions,
      // however 'nearest[:,j]' has only 'tn' tiled dimensional. The algorithm is only
      // sensitive to the differences between the 'ti'th and 'ti+1'th dimension, so we can
      // freely set one of the dimensions to a constant value. We choose to set the 0th
      // dimension to zero, i.e. the (tn+1)-dimensional lattice point is
      //   y = (0, tiled dimensions of 'nearest[:,j]').
      double y[tn+1];
      y[0] = 0;
      for (size_t ti = 0; ti < tn; ++ti) {
        const size_t i = loc->tiling->tiled_idx[ti];
        y[ti+1] = gsl_matrix_get(nearest, i, j);
      }

      // Find the nearest point in An* to the point 'y', using the O(tn) Algorithm 2 given in:
      //   McKilliam et.al., "A linear-time nearest point algorithm for the lattice An*"
      //   in "International Symposium on Information Theory and Its Applications", ISITA2008,
      //   Auckland, New Zealand, 7-10 Dec. 2008. DOI: 10.1109/ISITA.2008.4895596
      // Notes:
      //   * Since Algorithm 2 uses 1-based arrays, we have to translate, e.g.:
      //       z_t in paper <---> z[tn-1] in C code
      //   * Line 6 in Algorithm 2 as written in the paper is in error, see correction below.
      //   * We are only interested in 'k', the generating integers of the nearest point
      //     'x = Q * k', therefore line 26 in Algorithm 2 is not included.
      INT8 k[tn+1];
      {

        // Lines 1--4, 20
        double z[tn+1], alpha = 0, beta = 0;
        size_t bucket[tn+1], link[tn+1];
        feclearexcept(FE_ALL_EXCEPT);
        for (size_t ti = 1; ti <= tn + 1; ++ti) {
          k[ti-1] = lround(y[ti-1]);   // Line 20, moved here to avoid duplicate round
          z[ti-1] = y[ti-1] - k[ti-1];
          alpha += z[ti-1];
          beta += z[ti-1]*z[ti-1];
          bucket[ti-1] = 0;
        }
        if (fetestexcept(FE_INVALID) != 0) {
          XLALPrintError("Rounding failed while finding nearest point #%zu:", j);
          for (size_t ti = 1; ti <= tn + 1; ++ti) {
            XLALPrintError(" %0.2e", y[ti-1]);
          }
          XLALPrintError("\n");
          XLAL_ERROR(XLAL_EFAILED);
        }

        // Lines 5--8
        // Notes:
        //   * Correction to line 6, as as written in McKilliam et.al.:
        //       ti = tn + 1 - (tn + 1)*floor(z_t + 0.5)
        //     should instead read
        //       ti = tn + 1 - floor((tn + 1)*(z_t + 0.5))
        //   * We also convert the floor() operation into an lround():
        //       ti = tn + 1 - lround((tn + 1)*(z_t + 0.5) - 0.5)
        //     to avoid a casting operation. Rewriting the line as:
        //       ti = lround((tn + 1)*(0.5 - z_t) + 0.5)
        //     appears to improve numerical robustness in some cases.
        //   * No floating-point exception checking needed for lround()
        //     here since its argument will be of order 'tn'.
        for (size_t tt = 1; tt <= tn + 1; ++tt) {
          const INT8 ti = lround((tn + 1)* (0.5 - z[tt-1]) + 0.5);
          link[tt-1] = bucket[ti-1];
          bucket[ti-1] = tt;
        }

        // Lines 9--10
        double D = beta - alpha*alpha / (tn + 1);
        size_t tm = 0;

        // Lines 11--19
        for (size_t ti = 1; ti <= tn + 1; ++ti) {
          size_t tt = bucket[ti-1];
          while (tt != 0) {
            alpha = alpha - 1;
            beta = beta - 2*z[tt-1] + 1;
            tt = link[tt-1];
          }
          double d = beta - alpha*alpha / (tn + 1);
          if (d < D) {
            D = d;
            tm = ti;
          }
        }

        // Lines 21--25
        for (size_t ti = 1; ti <= tm; ++ti) {
          size_t tt = bucket[ti-1];
          while (tt != 0) {
            k[tt-1] = k[tt-1] + 1;
            tt = link[tt-1];
          }
        }

      }

      // The nearest point in An* is the tn differences between k[1]...k[tn] and k[0]
      for (size_t ti = 0; ti < tn; ++ti) {
        nearest_int[ti] = k[ti+1] - k[0];
      }

    }
      break;

    default:
      XLAL_ERROR(XLAL_EINVAL, "Invalid lattice=%u", loc->tiling->lattice);
    }

    // Bound generating integers and return the (sequential) indexes of 'nearest[:,j]', if requested
    if (nearest_idxs != NULL) {

      // Bound generating integers
      {
        const LT_IndexTrie *trie = loc->index_trie;
        size_t ti = 0;
        while (ti < tn) {

          // If 'nearest_int[ti]' is outside parameter-space bounds
          if (nearest_int[ti] < trie->int_lower || nearest_int[ti] > trie->int_upper) {

            // Find the nearest point within the parameter-space bounds of the lattice tiling
            INT8 original_int[tn], poll_int[tn];
            memcpy(original_int, nearest_int, sizeof(original_int));
            double poll_min_distance = GSL_POSINF;
            XLALPrintInfo("%s: calling LT_PollIndexTrie()\n", __func__);
            LT_PollIndexTrie(loc->tiling, loc->index_trie, 0, original_int, poll_int, &poll_min_distance, nearest_int);

            // Reset 'trie', given that 'nearest_int' may have changed in any dimension
            trie = loc->index_trie;
            ti = 0;
            continue;

          }

          // If we are below the highest dimension, jump to the next dimension based on 'nearest_int[ti]'
          if (ti + 1 < tn) {
            trie = &trie->next[nearest_int[ti] - trie->int_lower];
          }

          ++ti;

        }
      }

      // Return the (sequential) indexes of 'nearest[:,j]'
      {
        const LT_IndexTrie *trie = loc->index_trie;
        for (size_t ti = 0; ti < tn; ++ti) {
          const size_t i = loc->tiling->tiled_idx[ti];

          if (i < seqidx_ndim) {

            // If below or at 'seqidx_ndim', return the sequential index up to this dimension
            (*nearest_idxs)->data[ n * j + i] = trie->seqidx;
            if (i + 1 == seqidx_ndim) {

              // If at 'seqidx_ndim', need to also add the index in this dimension
              (*nearest_idxs)->data[ n * j + i] += nearest_int[ti] - trie->int_lower;

            }

          } else {

            // If above 'seqidx_ndim', return the index in this dimension
            (*nearest_idxs)->data[ n * j + i] = nearest_int[ti] - trie->int_lower;

          }

          // If we are below the highest dimension, jump to the next dimension based on 'nearest_int[ti]'
          if (ti + 1 < tn) {
            trie = &trie->next[nearest_int[ti] - trie->int_lower];
          }

        }
      }

    }

    // Set 'nearest[:,j]' to generating integers
    for (size_t ti = 0; ti < tn; ++ti) {
      const size_t i = loc->tiling->tiled_idx[ti];
      gsl_matrix_set(nearest, i, j, nearest_int[ti]);
    }

  }

  // Transform 'nearest' from generating integers to physical coordinates
  gsl_blas_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, loc->tiling->phys_from_int, nearest);

  // Add physical origin to every point in 'nearest'
  for (size_t i = 0; i < n; ++i) {
    const double phys_origin = gsl_vector_get(loc->tiling->phys_origin, i);
    gsl_vector_view nearest_row = gsl_matrix_row(nearest, i);
    gsl_vector_add_constant(&nearest_row.vector, phys_origin);
  }

  // Set any non-tiled dimensions in 'nearest'
  for (size_t i = 0; i < n; ++i) {
    if (!loc->tiling->bounds[i].is_tiled) {
      for (size_t j = 0; j < npoints; ++j) {
        gsl_vector_view nearest_col = gsl_matrix_column(nearest, j);

        // Get the physical bounds on the current dimension
        double phys_lower = 0.0, phys_upper = 0.0;
        LT_GetBounds(loc->tiling, i, &nearest_col.vector, &phys_lower, &phys_upper);

        // Set point to non-tiled parameter-space bound
        gsl_vector_set(&nearest_col.vector, i, phys_lower);

      }
    }
  }

  return XLAL_SUCCESS;

}

int XLALPrintLatticeTilingIndexTrie(
  const LatticeTilingLocator *loc,
  FILE *file
  )
{

  // Check input
  XLAL_CHECK(loc != NULL, XLAL_EFAULT);
  XLAL_CHECK(file != NULL, XLAL_EFAULT);

  const size_t tn = loc->tiled_ndim;

  // Return if index trie is NULL
  if (tn == 0 || loc->index_trie == NULL) {
    fprintf(file, "WARNING: index trie is NULL\n");
    return XLAL_SUCCESS;
  }

  // Print index trie
  INT8 int_lower[tn];
  LT_PrintIndexTrie(loc->tiling, loc->index_trie, 0, file, int_lower);

  return XLAL_SUCCESS;

}
