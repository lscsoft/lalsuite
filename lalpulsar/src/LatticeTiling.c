//
// Copyright (C) 2007, 2008, 2012, 2014 Karl Wette
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
#include <inttypes.h>
#include <sys/types.h>
#include <limits.h>
#include <math.h>
#include <fenv.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <lal/LatticeTiling.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALMalloc.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/XLALGSL.h>

#include "GSLHelpers.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

///
/// Lattice tiling parameter-space bound for one dimension
///
typedef struct tagLT_Bound {
  bool is_tiled;			///< True if the dimension is tiled, false if it is a single point
  LatticeTilingBound func;		///< Parameter space bound function
  void* data_lower;			///< Arbitrary data describing lower parameter-space bound
  void* data_upper;			///< Arbitrary data describing upper parameter-space bound
} LT_Bound;

///
/// Lattice tiling parameter-space bound trie for one dimension
///
typedef struct tagLT_BoundTrie LT_BoundTrie;
struct tagLT_BoundTrie {
  unsigned short is_tiled;		///< True if the dimension is tiled, false if it is a single point
  unsigned short is_highest;		///< True if this dimension is this highest
  const LT_BoundTrie* prev_trie;	///< Pointer to array of trie structs in next-lowest dimension
  union {
    struct {
      int32_t int_lower;		///< Integer lower bound on a tiled dimension
      int32_t int_upper;		///< Integer upper bound on a tiled dimension
    } tiled;
    double phys_point;			///< Physical point in a non-tiled dimension
  } bounds;
  union {
    LT_BoundTrie* next_trie;		///< Pointer to array of trie structs in next-highest dimension
    UINT8 tiling_index;			///< Lattice tiling index in the highest dimension
  } stored;
};

struct tagLatticeTilingSpace {
  size_t ndim;				///< Number of parameter-space dimensions
  LT_Bound* bounds;			///< Array of parameter-space bound info for each dimension
};

struct tagLatticeTiling {
  size_t ndim;				///< Number of parameter-space dimensions
  TilingLattice lattice;		///< Type of lattice tiling was generated with
  size_t tiled_ndim;			///< Number of tiled parameter-space dimensions
  size_t* tiled_idx;			///< Index to tiled parameter-space dimensions
  gsl_vector* phys_origin;		///< Parameter-space origin in physical coordinates
  gsl_vector* phys_bbox;		///< Metric ellipse bounding box
  gsl_matrix* int_from_phys;		///< Transform to generating integers from physical coordinates
  gsl_matrix* phys_from_int;		///< Transform to physical coordinates from generating integers
  LT_BoundTrie trie_base;		///< Base of lattice tiling parameter-space bound trie
  UINT8* point_counts;			///< Number of points cumulatively in each dimension
  UINT8* itr_link_counts;		///< Number of links to bound tries for iteration in each dimension
  LT_BoundTrie*** itr_links;		///< Links to bound tries for iteration in each dimension
};

struct tagLatticeTilingIterator {
  const LatticeTiling* tiling;		///< Lattice tiling state structure
  size_t dim;				///< Dimension being iterated
  gsl_vector* point;			///< Current parameter-space point
  UINT8 link_index;			///< Index to links to bound tries in iterated dimension
  UINT8 highest_count;			///< Number of remaining points in highest iterated dimension
  double highest_step;			///< Step size between points in highest iterated dimension
};

///
/// Zero out the strictly upper triangular part of the matrix A
///
static void LT_ZeroStrictUpperTriangle(gsl_matrix* A) {
  for (size_t i = 0; i < A->size1; ++i) {
    for (size_t j = i + 1; j < A->size2; ++j) {
      gsl_matrix_set(A, i, j, 0.0);
    }
  }
}

///
/// Reverse the order of both the rows and columns of the matrix A
///
static void LT_ReverseOrderRowsCols(gsl_matrix* A) {
  for (size_t i = 0; i < A->size1 / 2; ++i) {
    gsl_matrix_swap_rows(A, i, A->size1 - i - 1);
  }
  for (size_t j = 0; j < A->size2 / 2; ++j) {
    gsl_matrix_swap_columns(A, j, A->size2 - j - 1);
  }
}

LatticeTilingSpace* XLALCreateLatticeTilingSpace(
  const size_t n
  )
{

  // Check input
  XLAL_CHECK_NULL(n > 0, XLAL_EINVAL);

  // Allocate memory
  LatticeTilingSpace* space = XLALCalloc(1, sizeof(*space));
  XLAL_CHECK_NULL(space != NULL, XLAL_ENOMEM);
  space->bounds = XLALCalloc(n, sizeof(*space->bounds));
  XLAL_CHECK_NULL(space->bounds != NULL, XLAL_ENOMEM);

  // Initialise fields
  space->ndim = n;

  return space;

}

void XLALDestroyLatticeTilingSpace(
  LatticeTilingSpace* space
  )
{
  if (space != NULL) {
    if (space->bounds != NULL) {
      for (size_t i = 0; i < space->ndim; ++i) {
        XLALFree(space->bounds[i].data_lower);
        XLALFree(space->bounds[i].data_upper);
      }
      XLALFree(space->bounds);
    }
    XLALFree(space);
  }
}

int XLALSetLatticeTilingBound(
  LatticeTilingSpace* space,
  const size_t dim,
  const LatticeTilingBound func,
  const size_t data_len,
  void* data_lower,
  void* data_upper
  )
{

  // Check input
  XLAL_CHECK(space != NULL, XLAL_EFAULT);
  XLAL_CHECK(dim < space->ndim, XLAL_ESIZE);
  XLAL_CHECK(func != NULL, XLAL_EFAULT);
  XLAL_CHECK(data_len > 0, XLAL_EFAULT);
  XLAL_CHECK(data_lower != NULL, XLAL_EFAULT);
  XLAL_CHECK(data_upper != NULL, XLAL_EFAULT);

  // Check that bound has not already been set
  XLAL_CHECK(space->bounds[dim].func == NULL, XLAL_EINVAL, "Dimension #%zu is already bounded", dim);

  // Determine if this dimension is tiled
  const bool is_tiled = (memcmp(data_lower, data_upper, data_len) != 0);

  // Set the parameter-space bound
  space->bounds[dim].is_tiled = is_tiled;
  space->bounds[dim].func = func;
  space->bounds[dim].data_lower = data_lower;
  space->bounds[dim].data_upper = data_upper;

  return XLAL_SUCCESS;

}

static double ConstantBound(
  const void* data,
  const size_t dim UNUSED,
  const gsl_vector* point UNUSED
  )
{

  // Return bound
  return *((const double*) data);

}

int XLALSetLatticeTilingConstantBound(
  LatticeTilingSpace* space,
  const size_t dim,
  const double bound1,
  const double bound2
  )
{

  // Check input
  XLAL_CHECK(space != NULL, XLAL_EFAULT);
  XLAL_CHECK(isfinite(bound1), XLAL_EINVAL);
  XLAL_CHECK(isfinite(bound2), XLAL_EINVAL);

  // Allocate memory
  const size_t data_len = sizeof(double);
  double* data_lower = XLALMalloc(data_len);
  XLAL_CHECK(data_lower != NULL, XLAL_ENOMEM);
  double* data_upper = XLALMalloc(data_len);
  XLAL_CHECK(data_upper != NULL, XLAL_ENOMEM);

  // Set the parameter-space bound
  *data_lower = GSL_MIN(bound1, bound2);
  *data_upper = GSL_MAX(bound1, bound2);
  XLAL_CHECK(XLALSetLatticeTilingBound(space, dim, ConstantBound, data_len, data_lower, data_upper)
             == XLAL_SUCCESS, XLAL_EFUNC);

  return XLAL_SUCCESS;

}

///
/// Return the lower and upper parameter-space bounds on a given dimension
///
static void LT_GetPhysBounds(
  const LatticeTilingSpace* space,	///< [in] Lattice tiling parameter space
  const size_t dim,			///< [in] Dimension on which bound applies
  const gsl_vector* phys_point,		///< [in] Physical point at which to find bounds
  double* phys_lower,			///< [in,out] Lower parameter-space bound
  double* phys_upper			///< [in,out] Upper parameter-space bound
  )
{

  // Get bound information for this dimension
  const LT_Bound* bound = &space->bounds[dim];

  // Get view of first (dimension) dimensions of physical point
  gsl_vector_const_view phys_point_subv_view = gsl_vector_const_subvector(phys_point, 0, GSL_MAX(1, dim));
  const gsl_vector* phys_point_subv = (dim == 0) ? NULL : &phys_point_subv_view.vector;

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

int XLALRandomLatticeTilingPoints(
  const LatticeTilingSpace* space,
  RandomParams* rng,
  gsl_matrix* random_points
  )
{

  // Check input
  XLAL_CHECK(space != NULL, XLAL_EFAULT);
  const size_t n = space->ndim;
  XLAL_CHECK(rng != NULL, XLAL_EFAULT);
  XLAL_CHECK(random_points != NULL, XLAL_EFAULT);
  XLAL_CHECK(random_points->size1 == n, XLAL_ESIZE);

  // Generate random points in lattice tiling parameter space
  for (size_t k = 0; k < random_points->size2; ++k) {
    gsl_vector_view phys_point = gsl_matrix_column(random_points, k);
    for (size_t i = 0; i < n; ++i) {

      // Get the physical bounds on the current dimension
      double phys_lower = 0.0, phys_upper = 0.0;
      LT_GetPhysBounds(space, i, &phys_point.vector, &phys_lower, &phys_upper);

      // Generate random number
      const double u = XLALUniformDeviate(rng);

      // Set parameter space point
      gsl_vector_set(&phys_point.vector, i, phys_lower + u*(phys_upper - phys_lower));

    }

  }

  return XLAL_SUCCESS;

}

///
/// Compute the extent of the bounding box of the mismatch ellipses of a metric
///
static gsl_vector* LT_MetricEllipseBoundingBox(
  const gsl_matrix* metric,		///< [in] Metric to bound
  const double max_mismatch		///< [in] Maximum mismatch with respect to metric
  )
{

  // Check input
  XLAL_CHECK_NULL(metric != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(metric->size1 == metric->size2, XLAL_ESIZE);
  const size_t n = metric->size1;

  // Allocate memory
  gsl_matrix* GAMAT_NULL(LU_decomp, n, n);
  gsl_permutation* GAPERM_NULL(LU_perm, n);
  gsl_matrix* GAMAT_NULL(inverse, n, n);
  gsl_vector* GAVEC_NULL(bounding_box, n);

  // Compute metric inverse
  int LU_sign = 0;
  gsl_matrix_memcpy(LU_decomp, metric);
  GCALL_NULL(gsl_linalg_LU_decomp(LU_decomp, LU_perm, &LU_sign), "'metric' cannot be LU-decomposed");
  GCALL_NULL(gsl_linalg_LU_invert(LU_decomp, LU_perm, inverse), "'metric' cannot be inverted");

  // Compute bounding box, and reverse diagonal scaling
  for (size_t i = 0; i < n; ++i) {
    const double bounding_box_i = 2.0 * sqrt(max_mismatch * gsl_matrix_get(inverse, i, i));
    gsl_vector_set(bounding_box, i, bounding_box_i);
  }

  // Cleanup
  GFMAT(LU_decomp, inverse);
  GFPERM(LU_perm);

  return bounding_box;

}

///
/// Compute a lower-triangular basis matrix whose columns are orthonormal with respect to a given metric
///
static gsl_matrix* LT_ComputeMetricOrthoBasis(
  const gsl_matrix* metric		///< [in] Metric to orthonormalise with respect to
  )
{

  // Check input
  XLAL_CHECK_NULL(metric->size1 == metric->size2, XLAL_ESIZE);
  const size_t n = metric->size1;

  // Allocate memory
  gsl_matrix* GAMAT_NULL(basis, n, n);

  // We want to find a lower-triangular basis such that:
  //   basis^T * metric * basis = I
  // This is rearranged to give:
  //   metric^-1 = basis * basis^T
  // Hence basis is the Cholesky decomposition of metric^-1
  gsl_matrix_memcpy(basis, metric);
  GCALL_NULL(gsl_linalg_cholesky_decomp(basis), "'metric' is not positive definite");
  GCALL_NULL(gsl_linalg_cholesky_invert(basis), "'metric' cannot be inverted");
  GCALL_NULL(gsl_linalg_cholesky_decomp(basis), "Inverse of 'metric' is not positive definite");

  // gsl_linalg_cholesky_decomp() stores both basis and basis^T
  // in the same matrix; zero out upper triangle to get basis
  LT_ZeroStrictUpperTriangle(basis);

  return basis;

}

///
/// Compute a lower-triangular generator matrix for a given lattice type and mismatch
///
static gsl_matrix* LT_ComputeLatticeGenerator(
  const size_t n,			///< [in] Number of dimensions
  const TilingLattice lattice,		///< [in] Type of lattice
  const double max_mismatch		///< [in] Maximum prescribed mismatch
  )
{

  // Check input
  XLAL_CHECK_NULL(n > 0, XLAL_EINVAL);
  XLAL_CHECK_NULL(max_mismatch > 0.0, XLAL_EINVAL);

  // Allocate memory
  gsl_matrix* GAMAT_NULL(generator, n, n);
  gsl_matrix* GAMAT_NULL(LU_decomp, n, n);
  gsl_permutation* GAPERM_NULL(LU_perm, n);

  // Compute lattice generator and normalised thickness
  double norm_thickness = 0.0;
  switch (lattice) {

  case TILING_LATTICE_CUBIC:		// Cubic (\f$Z_n\f$) lattice
  {

    // Zn lattice generator is the identity
    gsl_matrix_set_identity(generator);

    // Zn normalised thickness
    norm_thickness = pow(sqrt(n)/2.0, n);

  }
  break;

  case TILING_LATTICE_ANSTAR:		// An-star (\f$A_n^*\f$) lattice
  {

    // Allocate memory
    gsl_matrix* GAMAT_NULL(G, n + 1, n + 1);
    gsl_vector* GAVEC_NULL(tau, n);
    gsl_matrix* GAMAT_NULL(Q, n + 1, n + 1);
    gsl_matrix* GAMAT_NULL(L, n + 1, n);

    // An* lattice generator in n+1 dimensions, given in:
    //   McKilliam et.al., "A linear-time nearest point algorithm for the lattice An*"
    //   in "International Symposium on Information Theory and Its Applications", ISITA2008,
    //   Auckland, New Zealand, 7-10 Dec. 2008. DOI: 10.1109/ISITA.2008.4895596
    gsl_matrix_set_identity(G);
    gsl_matrix_add_constant(G, -1.0 / (n + 1));

    // Find the QL decomposition of the generator matrix G, excluding 1st column,
    // which is linearly dependent on the remaining columns:
    //   G(:, 2:end) = Gp = Q * L
    // where Q is an orthogonal matrix and L an lower-triangular matrix.
    // This is found using the more commonly implemented QR decomposition by:
    // - reversing the order of the rows/columns of Gp
    // - decomposing Gp = Qp * Lp, where Lp is upper triangular
    // - reversing the order of the rows/columns of Qp to give Q
    // - reversing the order of the rows/columns of Lp to give L
    gsl_matrix_view Gp = gsl_matrix_submatrix(G, 0, 1, n + 1, n);
    LT_ReverseOrderRowsCols(&Gp.matrix);
    GCALL_NULL(gsl_linalg_QR_decomp(&Gp.matrix, tau), "'G' cannot be QR-decomposed");
    gsl_linalg_QR_unpack(&Gp.matrix, tau, Q, L);
    LT_ReverseOrderRowsCols(Q);
    LT_ReverseOrderRowsCols(L);

    // Discard the first row of L, which is zero, to get the generator in n dimensions
    gsl_matrix_view L_view = gsl_matrix_submatrix(L, 1, 0, n, n);
    gsl_matrix_memcpy(generator, &L_view.matrix);

    // Cleanup
    GFMAT(G, L, Q);
    GFVEC(tau);

    // Ans normalised thickness
    norm_thickness = sqrt(n+1.0) * pow( (1.0*n*(n+2.0)) / (12.0*(n+1.0)), 0.5*n );

  }
  break;

  default:
    XLAL_ERROR_NULL(XLAL_EINVAL, "Invalid lattice %u", lattice);
  }

  // Generator will be lower-triangular, so zero out upper triangle
  LT_ZeroStrictUpperTriangle(generator);

  // Ensure that the generator has positive diagonal elements, by
  // changing the sign of the columns of the matrix, if necessary
  for (size_t j = 0; j < generator->size2; ++j) {
    gsl_vector_view generator_col = gsl_matrix_column(generator, j);
    XLAL_CHECK_NULL(gsl_vector_get(&generator_col.vector, j) != 0, XLAL_ERANGE,
                    "Generator matrix generator(%zu,%zu) == 0", j, j);
    if (gsl_vector_get(&generator_col.vector, j) < 0) {
      gsl_vector_scale(&generator_col.vector, -1);
    }
  }

  // Compute generator LU decomposition
  gsl_matrix_memcpy(LU_decomp, generator);
  int LU_sign = 0;
  GCALL_NULL(gsl_linalg_LU_decomp(LU_decomp, LU_perm, &LU_sign), "'generator' cannot be LU-decomposed");

  // Compute generator determinant
  const double generator_determinant = gsl_linalg_LU_det(LU_decomp, LU_sign);

  // Compute generator covering radius
  const double generator_covering_radius = pow(norm_thickness * generator_determinant, 1.0 / n);

  // Normalise so covering spheres have sqrt(max_mismatch) covering radii
  gsl_matrix_scale(generator, sqrt(max_mismatch) / generator_covering_radius);

  // Cleanup
  GFMAT(LU_decomp);
  GFPERM(LU_perm);

  return generator;

}

//
// Find the extreme values of the physical bounds, by sampling the bounds around the current point
//
static void LT_FindBoundExtrema(
  const size_t dim,			///< [in] Current dimension of parameter space
  const size_t max_dim,			///< [in] Dimension in which to record bound extrema
  const LatticeTilingSpace* space,	///< [in] Lattice tiling parameter space
  const gsl_matrix* phys_from_int,	///< [in] Transform to physical coordinates from generating integers
  gsl_vector* phys_point,		///< [in/out] Current physical point being bounded
  double* phys_min_lower,		///< [in/out] Minimum lower bound on parameter space
  double* phys_max_upper		///< [in/out] Maximum upper bound on parameter space
  )
{

  if (dim < max_dim) {

    if (space->bounds[dim].is_tiled) {

      // Get the vector pointing to neighbouring points of the current point in this dimension
      gsl_vector_const_view phys_from_int_dim_view = gsl_matrix_const_row(phys_from_int, dim);

      for (int d = -1; d <= 1; ++d) {

        // Move current point half-way towards neighbouring point in the direction 'd'
        gsl_blas_daxpy(0.5 * d, &phys_from_int_dim_view.vector, phys_point);

        // Find bound extrema in higher dimensions
        LT_FindBoundExtrema(dim + 1, max_dim, space, phys_from_int, phys_point,
                            phys_min_lower, phys_max_upper);

        // Reset current point back to previous value
        gsl_blas_daxpy(-0.5 * d, &phys_from_int_dim_view.vector, phys_point);

      }

    } else {

      // Find bound extrema in higher dimensions
      LT_FindBoundExtrema(dim + 1, max_dim, space, phys_from_int, phys_point,
                          phys_min_lower, phys_max_upper);

    }

  } else {

    // Get the physical bounds on the current dimension
    double phys_lower = 0.0, phys_upper = 0.0;
    LT_GetPhysBounds(space, dim, phys_point, &phys_lower, &phys_upper);

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
/// Recursively build the lattice tiling parameter-space bound trie
///
static int LT_BuildBoundTrie(
  const size_t dim,			///< [in] Current dimension of parameter space
  const LatticeTilingSpace* space,	///< [in] Lattice tiling parameter space
  const LatticeTiling* tiling,		///< [in] Lattice tiling state structure
  gsl_vector* phys_point,		///< [in/out] Current physical point being bounded
  UINT8* point_counts,			///< [in/out] Count of number of points cumulatively in each dimension
  LT_BoundTrie* curr_trie		///< [in/out] Pointer to current trie struct for this dimension
  )
{

  // Set whether current trie is tiled and/or in highest dimension
  curr_trie->is_tiled = space->bounds[dim].is_tiled;
  curr_trie->is_highest = (dim + 1 == space->ndim);

  // Calculate integer conversion offset from lower dimensions
  double int_from_phys_offset = 0;
  if (curr_trie->is_tiled) {
    for (size_t j = 0; j < dim; ++j) {
      const double int_from_phys_dim_j = gsl_matrix_get(tiling->int_from_phys, dim, j);
      const double phys_origin_j = gsl_vector_get(tiling->phys_origin, j);
      int_from_phys_offset += int_from_phys_dim_j * (gsl_vector_get(phys_point, j) - phys_origin_j);
    }
  }

  // Get the physical bounds on the current dimension
  double phys_lower = 0.0, phys_upper = 0.0;
  LT_GetPhysBounds(space, dim, phys_point, &phys_lower, &phys_upper);

  if (curr_trie->is_tiled) {

    if (dim > 0) {

      // Create a copy of current physical point for use in finding bound extrema
      double test_phys_point_array[phys_point->size];
      gsl_vector_view test_phys_point_view = gsl_vector_view_array(test_phys_point_array, phys_point->size);
      gsl_vector_memcpy(&test_phys_point_view.vector, phys_point);

      // Find the extreme values of the physical bounds
      LT_FindBoundExtrema(0, dim, space, tiling->phys_from_int, &test_phys_point_view.vector,
                          &phys_lower, &phys_upper);

    }

    // Add padding of half the metric ellipse bounding box in this dimension
    const double phys_bbox_dim = 0.5 * gsl_vector_get(tiling->phys_bbox, dim);
    phys_lower -= phys_bbox_dim;
    phys_upper += phys_bbox_dim;

    // Convert parameter-space bounds to generating integers
    const double int_from_phys_dim_dim = gsl_matrix_get(tiling->int_from_phys, dim, dim);
    const double phys_origin_dim = gsl_vector_get(tiling->phys_origin, dim);
    const double int_lower_dbl = int_from_phys_offset + int_from_phys_dim_dim * (phys_lower - phys_origin_dim);
    const double int_upper_dbl = int_from_phys_offset + int_from_phys_dim_dim * (phys_upper - phys_origin_dim);

    // Store integer lower/upper bounds on a tiled dimension, rounded up/down to avoid extra
    // boundary points, and checking that bounds are within the range of an 'int32_t'
    feclearexcept(FE_ALL_EXCEPT);
    const long int int_lower = lround(ceil(int_lower_dbl));
    const long int int_upper = lround(floor(int_upper_dbl));
    XLAL_CHECK(fetestexcept(FE_INVALID) == 0, XLAL_EFAILED,
               "Bounds on dimension #%zu are too large (%0.2e to %0.2e)", dim, int_lower_dbl, int_upper_dbl);
    XLAL_CHECK(INT32_MIN <= int_lower && int_lower <= INT32_MAX &&
               INT32_MIN <= int_upper && int_upper <= INT32_MAX, XLAL_EFAILED,
               "Bounds on dimension #%zu are too large (%li to %li)", dim, int_lower, int_upper);
    XLAL_CHECK(int_lower <= int_upper, XLAL_EFAILED,
               "Bounds on dimension #%zu are inverted! (%li !<= %li)", dim, int_lower, int_upper);
    curr_trie->bounds.tiled.int_lower = ((int32_t) int_lower);
    curr_trie->bounds.tiled.int_upper = ((int32_t) int_upper);

  } else {

    // Store physical point in a non-tiled dimension
    curr_trie->bounds.phys_point = phys_lower;

  }

  // Number of points in this dimension
  const size_t point_count = curr_trie->is_tiled ?
    (curr_trie->bounds.tiled.int_upper - curr_trie->bounds.tiled.int_lower + 1) : 1;

  if (curr_trie->is_highest) {

    // Store current tiling index
    curr_trie->stored.tiling_index = point_counts[dim];

  } else {

    // Allocate memory
    curr_trie->stored.next_trie = XLALCalloc(point_count, sizeof(*curr_trie->stored.next_trie));
    XLAL_CHECK(curr_trie->stored.next_trie != NULL, XLAL_ENOMEM);

    if (curr_trie->is_tiled) {

      // Convert integer lower bound to physical coordinates
      const double phys_from_int_dim_dim = gsl_matrix_get(tiling->phys_from_int, dim, dim);
      const double int_lower_dbl = ((double) curr_trie->bounds.tiled.int_lower);
      const double phys_origin_dim = gsl_vector_get(tiling->phys_origin, dim);
      double phys_point_dim = phys_origin_dim + phys_from_int_dim_dim * (int_lower_dbl - int_from_phys_offset);

      for (size_t i = 0; i < point_count; ++i) {

        // Set physical point in this dimension to current point
        gsl_vector_set(phys_point, dim, phys_point_dim);

        // Build bound tries in higher dimensions
        XLAL_CHECK(LT_BuildBoundTrie(dim + 1, space, tiling, phys_point,
                                     point_counts, &curr_trie->stored.next_trie[i])
                   == XLAL_SUCCESS, XLAL_EFUNC);

        // Increment physical point
        phys_point_dim += phys_from_int_dim_dim;

      }

    } else {

      // Set physical point in this dimension
      gsl_vector_set(phys_point, dim, curr_trie->bounds.phys_point);

      // Build bound trie in higher dimensions
      XLAL_CHECK(LT_BuildBoundTrie(dim + 1, space, tiling, phys_point,
                                   point_counts, curr_trie->stored.next_trie)
                 == XLAL_SUCCESS, XLAL_EFUNC);

    }

  }

  // Increment point count
  point_counts[dim] += point_count;

  return XLAL_SUCCESS;

}

///
/// Recursively build links to bound tries for iteration in each dimension
///
static int LT_BuildItrLinks(
  const size_t dim,			///< [in] Current dimension of parameter space
  const LT_BoundTrie* prev_trie,	///< [in] Pointer to trie struct in next-lowest dimension
  LT_BoundTrie* curr_trie,		///< [in] Pointer to current trie struct for this dimension
  LT_BoundTrie*** const itr_links,	///< [in] Links to bound tries cumulatively in each dimension
  UINT8* link_index			///< [in] Indices to links to bound tries
  )
{

  // Save link to previous trie in current trie
  curr_trie->prev_trie = prev_trie;

  // Save link to current trie in iteration links
  itr_links[dim][link_index[dim]++] = curr_trie;

  if (!curr_trie->is_highest) {

    // Number of points in this dimension
    const size_t point_count = curr_trie->is_tiled ?
      (curr_trie->bounds.tiled.int_upper - curr_trie->bounds.tiled.int_lower + 1) : 1;

    for (size_t i = 0; i < point_count; ++i) {

      // Build links in higher dimensions
      XLAL_CHECK(LT_BuildItrLinks(dim + 1, curr_trie, &curr_trie->stored.next_trie[i],
                                  itr_links, link_index)
                 == XLAL_SUCCESS, XLAL_EFUNC);

    }

  }

  return XLAL_SUCCESS;

}

///
/// Recursively free a lattice tiling parameter-space bound trie
///
static void LT_FreeBoundTrie(
  LT_BoundTrie* curr_trie		///< [in] Pointer to current trie struct for this dimension
  )
{

  if (!curr_trie->is_highest) {

    // Number of points in this dimension
    const size_t point_count = curr_trie->is_tiled ?
      (curr_trie->bounds.tiled.int_upper - curr_trie->bounds.tiled.int_lower + 1) : 1;

    for (size_t i = 0; i < point_count; ++i) {

      // Free trie in higher dimensions
      LT_FreeBoundTrie(&curr_trie->stored.next_trie[i]);

    }

    // Free memory in this dimension
    XLALFree(curr_trie->stored.next_trie);

  }

}

LatticeTiling* XLALCreateLatticeTiling(
  const LatticeTilingSpace* space,
  const TilingLattice lattice,
  const gsl_matrix* metric,
  const double max_mismatch
  )
{

  // Check input
  XLAL_CHECK_NULL(space != NULL, XLAL_EFAULT);
  const size_t n = space->ndim;
  XLAL_CHECK_NULL(lattice < TILING_LATTICE_MAX, XLAL_EINVAL);
  XLAL_CHECK_NULL(metric != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(metric->size1 == n && metric->size2 == n, XLAL_EINVAL);
  XLAL_CHECK_NULL(max_mismatch > 0, XLAL_EINVAL);

  // Check that all parameter-space dimensions are bounded
  for (size_t i = 0; i < n; ++i) {
    XLAL_CHECK_NULL(space->bounds[i].func != NULL, XLAL_EFAILED,
                    "Lattice tiling dimension #%zu is unbounded", i);
  }

  // Check metric is symmetric and has positive diagonal elements
  for (size_t i = 0; i < n; ++i) {
    XLAL_CHECK_NULL(gsl_matrix_get(metric, i, i) > 0, XLAL_EINVAL,
                    "Parameter-space metric(%zu,%zu) <= 0", i, i);
    for (size_t j = i + 1; j < n; ++j) {
      XLAL_CHECK_NULL(gsl_matrix_get(metric, i, j) == gsl_matrix_get(metric, j, i), XLAL_EINVAL,
                      "Parameter-space metric(%zu,%zu) != metric(%zu,%zu)", i, j, j, i);
    }
  }

  // Allocate and initialise memory
  LatticeTiling* tiling = XLALCalloc(1, sizeof(*tiling));
  XLAL_CHECK_NULL(tiling != NULL, XLAL_ENOMEM);
  tiling->point_counts = XLALCalloc(n, sizeof(*tiling->point_counts));
  XLAL_CHECK_NULL(tiling->point_counts != NULL, XLAL_ENOMEM);
  tiling->itr_link_counts = XLALCalloc(n, sizeof(*tiling->itr_link_counts));
  XLAL_CHECK_NULL(tiling->itr_link_counts != NULL, XLAL_ENOMEM);
  tiling->itr_links = XLALCalloc(n, sizeof(*tiling->itr_links));
  XLAL_CHECK_NULL(tiling->itr_links != NULL, XLAL_ENOMEM);
  GAVEC_NULL(tiling->phys_origin, n);
  GAVEC_NULL(tiling->phys_bbox, n);
  GAMAT_NULL(tiling->int_from_phys, n, n);
  GAMAT_NULL(tiling->phys_from_int, n, n);
  gsl_matrix_set_identity(tiling->int_from_phys);
  gsl_matrix_set_identity(tiling->phys_from_int);

  // Save fields
  tiling->ndim = n;
  tiling->lattice = lattice;

  // Count number of tiled dimensions
  tiling->tiled_ndim = 0;
  for (size_t i = 0; i < n; ++i) {
    if (space->bounds[i].is_tiled) {
      ++tiling->tiled_ndim;
    }
  }

  if (tiling->tiled_ndim > 0) {
    const size_t tn = tiling->tiled_ndim;

    // Allocate memory
    tiling->tiled_idx = XLALCalloc(tn, sizeof(*tiling->tiled_idx));
    XLAL_CHECK_NULL(tiling->tiled_idx != NULL, XLAL_ENOMEM);

    // Record indices of tiled dimensions
    for (size_t i = 0, ti = 0; i < n; ++i) {
      if (space->bounds[i].is_tiled) {
        tiling->tiled_idx[ti++] = i;
      }
    }

    // Calculate normalisation scale from metric diagonal elements
    gsl_vector* GAVEC_NULL(t_norm, tn);
    for (size_t ti = 0; ti < tn; ++ti) {
      const size_t i = tiling->tiled_idx[ti];
      const double metric_i_i = gsl_matrix_get(metric, i, i);
      gsl_vector_set(t_norm, ti, sqrt(metric_i_i));
    }

    // Copy and normalise tiled dimensions of metric
    gsl_matrix* GAMAT_NULL(t_metric, tn, tn);
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
    gsl_vector* t_bbox = LT_MetricEllipseBoundingBox(t_metric, max_mismatch);
    XLAL_CHECK_NULL(t_bbox != NULL, XLAL_EFUNC);

    // Copy bounding box in physical coordinates to tiled dimensions
    for (size_t ti = 0; ti < tn; ++ti) {
      const size_t i = tiling->tiled_idx[ti];
      const double t_norm_ti = gsl_vector_get(t_norm, ti);
      gsl_vector_set(tiling->phys_bbox, i, gsl_vector_get(t_bbox, ti) / t_norm_ti);
    }

    // Set physical parameter-space origin to mid-point of parameter-space bounds
    for (size_t i = 0; i < n; ++i) {
      double phys_lower = 0.0, phys_upper = 0.0;
      LT_GetPhysBounds(space, i, tiling->phys_origin, &phys_lower, &phys_upper);
      gsl_vector_set(tiling->phys_origin, i, 0.5*(phys_lower + phys_upper));
    }

    // Set non-tiled dimensions of physical parameter-space origin back to zero
    for (size_t i = 0; i < n; ++i) {
      if (!space->bounds[i].is_tiled) {
        gsl_vector_set(tiling->phys_origin, i, 0);
      }
    }

    // Compute a lower-triangular basis matrix whose columns are orthonormal with respect to the tiled metric
    gsl_matrix* t_basis = LT_ComputeMetricOrthoBasis(t_metric);
    XLAL_CHECK_NULL(t_basis != NULL, XLAL_EFUNC);

    // Compute a lower-triangular generator matrix for a given lattice type and mismatch
    gsl_matrix* t_gen = LT_ComputeLatticeGenerator(tn, lattice, max_mismatch);
    XLAL_CHECK_NULL(t_gen != NULL, XLAL_EFUNC);

    // Compute transform to normalised physical coordinates from generating integers
    gsl_matrix* GAMAT_NULL(t_norm_from_int, tn, tn);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, t_basis, t_gen, 0.0, t_norm_from_int);
    LT_ZeroStrictUpperTriangle(t_norm_from_int);

    // Compute transform to generating integers from normalised physical coordinates
    gsl_matrix* GAMAT_NULL(t_int_from_norm, tn, tn);
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
    GFMAT(t_metric, t_basis, t_gen, t_norm_from_int, t_int_from_norm);
    GFVEC(t_norm, t_bbox);

  }

  // Build parameter-space bound trie
  double phys_point_array[n];
  gsl_vector_view phys_point_view = gsl_vector_view_array(phys_point_array, n);
  XLAL_CHECK_NULL(LT_BuildBoundTrie(0, space, tiling, &phys_point_view.vector,
                                    tiling->point_counts, &tiling->trie_base)
                  == XLAL_SUCCESS, XLAL_EFUNC);

  // Initialise number of links to bound tries for iteration in each dimension
  tiling->itr_link_counts[0] = 1;
  for (size_t i = 0; i + 1< n; ++i) {
    tiling->itr_link_counts[i + 1] = tiling->point_counts[i];
  }

  // Allocate and initialise memory
  UINT8 link_index[n];
  for (size_t i = 0; i < n; ++i) {
    tiling->itr_links[i] = XLALCalloc(tiling->point_counts[i], sizeof(*tiling->itr_links[i]));
    XLAL_CHECK_NULL(tiling->itr_links[i] != NULL, XLAL_ENOMEM);
    link_index[i] = 0;
  }

  // Build links within bound trie
  XLAL_CHECK_NULL(LT_BuildItrLinks(0, NULL, &tiling->trie_base, tiling->itr_links, link_index)
                  == XLAL_SUCCESS, XLAL_EFUNC);
  for (size_t i = 0; i < n; ++i) {
    XLAL_CHECK_NULL(link_index[i] == tiling->itr_link_counts[i], XLAL_EFAILED);
  }

  return tiling;

}

void XLALDestroyLatticeTiling(
  LatticeTiling* tiling
  )
{
  if (tiling != NULL) {
    const size_t n = tiling->ndim;
    LT_FreeBoundTrie(&tiling->trie_base);
    GFMAT(tiling->int_from_phys, tiling->phys_from_int);
    GFVEC(tiling->phys_origin, tiling->phys_bbox);
    XLALFree(tiling->tiled_idx);
    XLALFree(tiling->point_counts);
    XLALFree(tiling->itr_link_counts);
    for (size_t i = 0; i < n; ++i) {
      XLALFree(tiling->itr_links[i]);
    }
    XLALFree(tiling->itr_links);
    XLALFree(tiling);
  }
}

size_t XLALTotalLatticeTilingDimensions(
  const LatticeTiling* tiling
  )
{

  // Check input
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);

  return tiling->ndim;

}

size_t XLALTiledLatticeTilingDimensions(
  const LatticeTiling* tiling
  )
{

  // Check input
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);

  return tiling->tiled_ndim;

}

UINT8 XLALLatticeTilingTotalPointCount(
  const LatticeTiling* tiling
  )
{

  // Check input
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);

  return tiling->point_counts[tiling->ndim - 1];

}

UINT8 XLALLatticeTilingUniquePointCount(
  const LatticeTiling* tiling,
  const size_t dim
  )
{

  // Check input
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_VAL(0, dim < tiling->ndim, XLAL_ESIZE);

  return tiling->point_counts[dim];

}

gsl_matrix* XLALLatticeTilingUniquePoints(
  const LatticeTiling* tiling,
  const size_t dim
  )
{

  // Check input
  XLAL_CHECK_NULL(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(dim < tiling->ndim, XLAL_ESIZE);

  // Allocate memory
  gsl_matrix* GAMAT_NULL(points, dim + 1, tiling->point_counts[dim]);

  // Iterate over points and add to matrix
  LatticeTilingIterator* itr = XLALCreateLatticeTilingIterator(tiling, dim);
  XLAL_CHECK_NULL(itr != NULL, XLAL_EFUNC);
  for (size_t i = 0; i < points->size2; ++i) {
    gsl_vector_view point_i = gsl_matrix_column(points, i);
    XLAL_CHECK_NULL(XLALNextLatticeTilingPoint(itr, &point_i.vector) > 0, XLAL_EFUNC);
  }
  XLAL_CHECK_NULL(XLALNextLatticeTilingPoint(itr, NULL) == 0, XLAL_EFUNC);

  // Cleanup
  XLALDestroyLatticeTilingIterator(itr);

  return points;

}

LatticeTilingIterator* XLALCreateLatticeTilingIterator(
  const LatticeTiling* tiling,
  const size_t dim
  )
{

  // Check input
  XLAL_CHECK_NULL(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(dim < tiling->ndim, XLAL_ESIZE);

  // Allocate memory
  LatticeTilingIterator* itr = XLALCalloc(1, sizeof(*itr));
  XLAL_CHECK_NULL(itr != NULL, XLAL_ENOMEM);
  GAVEC_NULL(itr->point, dim + 1);

  // Initialise fields
  itr->tiling = tiling;
  itr->dim = dim;
  itr->link_index = 0;
  itr->highest_count = 0;
  itr->highest_step = 0.0;

  return itr;

}

void XLALDestroyLatticeTilingIterator(
  LatticeTilingIterator* itr
  )
{
  if (itr != NULL) {
    GFVEC(itr->point);
    XLALFree(itr);
  }
}

int XLALNextLatticeTilingPoint(
  LatticeTilingIterator* itr,
  gsl_vector* point
  )
{

  // Check input
  XLAL_CHECK(itr != NULL, XLAL_EFAULT);
  XLAL_CHECK(itr->tiling != NULL, XLAL_EFAULT);
  const LatticeTiling* const tiling = itr->tiling;
  XLAL_CHECK(point == NULL || point->size == itr->dim + 1, XLAL_ESIZE);

  if (itr->highest_count > 0) {

    // Decrement the number of remaining points in highest iterated dimension
    --itr->highest_count;

    // Increment the current point by the step size in the highest iterated dimension
    *gsl_vector_ptr(itr->point, itr->dim) += itr->highest_step;

  } else {

    // If no more trie links to iterate over, iterator is finished
    if (itr->link_index >= tiling->itr_link_counts[itr->dim]) {
      return 0;
    }

    // Get link to bound trie this iterator currently points to
    const LT_BoundTrie* curr_trie = tiling->itr_links[itr->dim][itr->link_index];
    XLAL_CHECK(curr_trie != NULL, XLAL_EFAILED);

    // Advance bound trie link index
    ++itr->link_index;

    if (curr_trie->is_tiled) {

      // Set point to generating integer lower bound in a tiled dimension
      gsl_vector_set(itr->point, itr->dim, curr_trie->bounds.tiled.int_lower);

    } else {

      // Set point to physical point in a non-tiled dimension
      gsl_vector_set(itr->point, itr->dim, curr_trie->bounds.phys_point);

    }

    // Save number of remaining points and step size in highest iterated dimension
    itr->highest_count = curr_trie->is_tiled ?
      (curr_trie->bounds.tiled.int_upper - curr_trie->bounds.tiled.int_lower) : 0;
    itr->highest_step = curr_trie->is_tiled ?
      gsl_matrix_get(itr->tiling->phys_from_int, itr->dim, itr->dim) : 0.0;

    // Move through lower dimensions of bound trie until base is reached
    const LT_BoundTrie* prev_trie = curr_trie;
    curr_trie = curr_trie->prev_trie;
    size_t dim = itr->dim;
    while (curr_trie != NULL) {
      --dim;

      if (curr_trie->is_tiled) {

        // Work out generating integer offset from difference in trie pointers
        const ptrdiff_t trie_offset = prev_trie - curr_trie->stored.next_trie;

        // Set point to generating integer in a tiled dimension
        gsl_vector_set(itr->point, dim, curr_trie->bounds.tiled.int_lower + trie_offset);

      } else {

        // Set point to physical point in a non-tiled dimension
        gsl_vector_set(itr->point, dim, curr_trie->bounds.phys_point);

      }

      // Move to lower dimension of bound trie
      prev_trie = curr_trie;
      curr_trie = curr_trie->prev_trie;

    }

    // Transform 'itr->point' from generating integers to physical coordinates
    gsl_matrix_const_view phys_from_int_view =
      gsl_matrix_const_submatrix(tiling->phys_from_int, 0, 0, itr->point->size, itr->point->size);
    gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, &phys_from_int_view.matrix, itr->point);

    // Add physical parameter-space origin
    gsl_vector_const_view phys_origin_view =
      gsl_vector_const_subvector(tiling->phys_origin, 0, itr->point->size);
    gsl_vector_add(itr->point, &phys_origin_view.vector);

  }

  // Copy 'itr->point' to output 'point'
  if (point != NULL) {
    gsl_vector_memcpy(point, itr->point);
  }

  return 1;

}

///
/// Find the nearest point within the parameter-space bounds of the lattice tiling,
/// by polling neighbours of a nearest point found with LT_FindNearestPoint()
///
static void LT_PollNearestPoint(
  const size_t tdim,			///< [in] Current tiled dimension of parameter space
  const size_t* tiled_idx,		///< [in] Index to tiled parameter-space dimensions
  const LT_BoundTrie* curr_trie,	///< [in] Pointer to current trie struct for this dimension
  const gsl_vector* nearest_point,	///< [in] Point for which nearest point is sought
  const long int* nearest_int_point,	///< [in] Nearest point found by LT_FindNearestPoint()
  const double poll_distance,		///< [in] Distance to currently polled point from lower dimensions
  long int* poll_int_point,		///< [in] Point currently being polled
  long int* poll_nearest_int_point,	///< [in] Nearest point found by polling
  double* poll_min_distance		///< [in] Minimum distance to nearest point found by polling
  )
{

  if (curr_trie->is_tiled) {

    // Get integer lower and upper bounds
    const long int int_lower = curr_trie->bounds.tiled.int_lower;
    const long int int_upper = curr_trie->bounds.tiled.int_upper;

    // Poll points within 1 of nearest point found by LT_FindNearestPoint(), but within bounds
    const long int poll_lower = GSL_MAX(nearest_int_point[tdim] - 1, int_lower);
    const long int poll_upper = GSL_MIN(nearest_int_point[tdim] + 1, int_upper);

    for (poll_int_point[tdim] = poll_lower; poll_int_point[tdim] <= poll_upper; ++poll_int_point[tdim]) {

      // Compute offset between point for which nearest point is sought, and polled point
      const size_t dim = tiled_idx[tdim];
      const double d = gsl_vector_get(nearest_point, dim) - poll_int_point[tdim];

      // Add to distance between point for which nearest point is sought, and polled point
      const double poll_distance_i = poll_distance + d*d;

      if (!curr_trie->is_highest) {

        // Set index to trie in next higher dimension
        const size_t next_trie_index = poll_int_point[tdim] - int_lower;

        // Continue polling in higher dimensions
        LT_PollNearestPoint(tdim + 1, tiled_idx, &curr_trie->stored.next_trie[next_trie_index],
                            nearest_point, nearest_int_point, poll_distance_i,
                            poll_int_point, poll_nearest_int_point, poll_min_distance);

      } else if (poll_distance_i < *poll_min_distance) {

        // If distance is smaller than minimum distance, record current polled point as nearest
        *poll_min_distance = poll_distance_i;
        memcpy(poll_nearest_int_point, poll_int_point, (tdim+1)*sizeof(poll_nearest_int_point[0]));

      }

    }

  } else if (!curr_trie->is_highest) {

    // Continue polling in higher dimensions
    LT_PollNearestPoint(tdim, tiled_idx, &curr_trie->stored.next_trie[0],
                        nearest_point, nearest_int_point, poll_distance,
                        poll_int_point, poll_nearest_int_point, poll_min_distance);

  }

}

///
/// Find the nearest point, and its tiling index, in the lattice tiling
///
static void LT_FindNearestPoint(
  const LatticeTiling* tiling,		///< [in] Lattice tiling state structure
  gsl_vector* nearest_point,		///< [out] Point to replace with nearest point
  UINT8* nearest_index			///< [in] Tiling index of nearest point
  )
{

  const size_t tn = tiling->tiled_ndim;
  long int nearest_int_point[tn];

  if (tn > 0) {

    // Find the nearest point to 'nearest_point', the tiled dimensions of which are generating integers
    switch (tiling->lattice) {

    case TILING_LATTICE_CUBIC:		// Cubic (\f$Z_n\f$) lattice
    {

      // Round each dimension of 'nearest_point' to nearest integer to find the nearest point in Zn
      for (size_t ti = 0; ti < tn; ++ti) {
        const size_t i = tiling->tiled_idx[ti];
        nearest_int_point[ti] = lround(gsl_vector_get(nearest_point, i));
      }

    }
    break;

    case TILING_LATTICE_ANSTAR:		// An-star (\f$A_n^*\f$) lattice
    {

      // The nearest point algorithm used below embeds the An* lattice in tn+1 dimensions,
      // however 'nearest_point' has only 'tn' tiled dimensional. The algorithm is only
      // sensitive to the differences between the (ti)th and (ti+1)th dimension, so we can
      // freely set one of the dimensions to a constant value. We choose to set the 0th
      // dimension to zero, i.e. the (tn+1)-dimensional lattice point is
      //   y = (0, tiled dimensions of 'nearest_point').
      double y[tn+1];
      y[0] = 0;
      for (size_t ti = 0; ti < tn; ++ti) {
        const size_t i = tiling->tiled_idx[ti];
        y[ti+1] = gsl_vector_get(nearest_point, i);
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
      long int k[tn+1];
      {

        // Lines 1--4, 20
        double z[tn+1], alpha = 0, beta = 0;
        size_t bucket[tn+1], link[tn+1];
        for (size_t ti = 1; ti <= tn + 1; ++ti) {
          k[ti-1] = lround(y[ti-1]);   // Line 20, moved here to avoid duplicate round
          z[ti-1] = y[ti-1] - k[ti-1];
          alpha += z[ti-1];
          beta += z[ti-1]*z[ti-1];
          bucket[ti-1] = 0;
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
        for (size_t tt = 1; tt <= tn + 1; ++tt) {
          const size_t ti = lround((tn + 1)*(0.5 - z[tt-1]) + 0.5);
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
        nearest_int_point[ti] = k[ti+1] - k[0];
      }

    }
    break;

    default:
      XLAL_ERROR_VOID(XLAL_EINVAL, "Invalid lattice %u", tiling->lattice);
    }

  }

  // Check nearest point is within parameter-space bounds
  {
    size_t ti = 0;
    const LT_BoundTrie* curr_trie = &tiling->trie_base;
    while (curr_trie != NULL) {

      size_t next_trie_index = 0;

      if (curr_trie->is_tiled) {

        // Get integer lower and upper bounds
        const long int int_lower = curr_trie->bounds.tiled.int_lower;
        const long int int_upper = curr_trie->bounds.tiled.int_upper;

        // If nearest point is outside parameter-space bounds ...
        if (nearest_int_point[ti] < int_lower || nearest_int_point[ti] > int_upper) {

          /// Find the nearest point within the parameter-space bounds of the lattice tiling
          long int poll_int_point[tn], poll_nearest_int_point[tn];
          double poll_min_distance = GSL_POSINF;
          memcpy(poll_nearest_int_point, nearest_int_point, sizeof(poll_nearest_int_point));
          LT_PollNearestPoint(0, tiling->tiled_idx, &tiling->trie_base,
                              nearest_point, nearest_int_point, 0.0,
                              poll_int_point, poll_nearest_int_point, &poll_min_distance);
          memcpy(nearest_int_point, poll_nearest_int_point, sizeof(nearest_int_point));

          break;

        }

        // Set index to trie in next higher dimension
        next_trie_index = nearest_int_point[ti] - curr_trie->bounds.tiled.int_lower;

        ++ti;

      }

      if (curr_trie->is_highest) {

        curr_trie = NULL;

      } else {

        // Jump to trie in next higher dimension based in 'next_trie_index'
        curr_trie = &curr_trie->stored.next_trie[next_trie_index];

      }

    }
  }

  // Set nearest point
  {
    size_t i = 0, ti = 0;
    const LT_BoundTrie* curr_trie = &tiling->trie_base;
    while (curr_trie != NULL) {

      size_t next_trie_index = 0;

      if (curr_trie->is_tiled) {

        // Set 'nearest_point' to 'nearest_int_point'
        gsl_vector_set(nearest_point, i, nearest_int_point[ti]);

        // Set index to trie in next higher dimension
        next_trie_index = nearest_int_point[ti] - curr_trie->bounds.tiled.int_lower;

        ++i;
        ++ti;

      } else {

        // Set 'nearest_point' to physical point in this dimension
        gsl_vector_set(nearest_point, i, curr_trie->bounds.phys_point);

        ++i;

      }

      if (curr_trie->is_highest) {

        // Set tiling index of 'nearest_point'
        *nearest_index = curr_trie->stored.tiling_index + next_trie_index;

        curr_trie = NULL;

      } else {

        // Jump to trie in next higher dimension based in 'next_trie_index'
        curr_trie = &curr_trie->stored.next_trie[next_trie_index];

      }

    }
  }

}

int XLALNearestLatticeTilingPoints(
  const LatticeTiling* tiling,
  const gsl_matrix* points,
  gsl_matrix** nearest_points,
  UINT8Vector** nearest_indices
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  const size_t n = tiling->ndim;
  XLAL_CHECK(points != NULL, XLAL_EFAULT);
  XLAL_CHECK(points->size1 == n, XLAL_ESIZE);
  const size_t point_count = points->size2;
  XLAL_CHECK(nearest_points != NULL, XLAL_EFAULT);

  // (Re)Allocate nearest point matrix, and create view of correct size
  if (*nearest_points != NULL && ((*nearest_points)->size1 != n || (*nearest_points)->size2 < point_count)) {
    gsl_matrix_free(*nearest_points);
    *nearest_points = NULL;
  }
  if (*nearest_points == NULL) {
    GAMAT(*nearest_points, n, point_count);
  }
  gsl_matrix_view nearest_view = gsl_matrix_submatrix(*nearest_points, 0, 0, n, point_count);
  gsl_matrix *const nearest = &nearest_view.matrix;

  // (Re)Allocate nearest point tiling index vector, if supplied
  if (nearest_indices != NULL) {
    if (*nearest_indices != NULL && (*nearest_indices)->length < point_count) {
      XLALDestroyUINT8Vector(*nearest_indices);
      *nearest_indices = NULL;
    }
    if (*nearest_indices == NULL) {
      *nearest_indices = XLALCreateUINT8Vector(point_count);
      XLAL_CHECK(*nearest_indices != NULL, XLAL_ENOMEM);
    }
  }

  // Copy 'points' to 'nearest'
  gsl_matrix_memcpy(nearest, points);

  // Subtract physical origin from every point in 'nearest'
  for (size_t i = 0; i < n; ++i) {
    const double phys_origin = gsl_vector_get(tiling->phys_origin, i);
    gsl_vector_view nearest_row = gsl_matrix_row(nearest, i);
    gsl_vector_add_constant(&nearest_row.vector, -phys_origin);
  }

  // Transform 'nearest' from physical coordinates to generating integers
  gsl_blas_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, tiling->int_from_phys, nearest);

  // Find the nearest points in the lattice tiling to the points in 'nearest'
  xlalErrno = 0;
  for (size_t j = 0; j < point_count; ++j) {
    gsl_vector_view nearest_j = gsl_matrix_column(nearest, j);
    UINT8 nearest_index = 0;
    LT_FindNearestPoint(tiling, &nearest_j.vector, &nearest_index);
    if (nearest_indices != NULL) {
      (*nearest_indices)->data[j] = nearest_index;
    }
  }
  XLAL_CHECK(xlalErrno == 0, XLAL_EFAILED, "Unexpected failure in LT_FindNearestPoint()");

  // Transform 'nearest' from generating integers to physical coordinates
  gsl_blas_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, tiling->phys_from_int, nearest);

  // Add physical origin to every point in 'nearest'
  for (size_t i = 0; i < n; ++i) {
    const double phys_origin = gsl_vector_get(tiling->phys_origin, i);
    gsl_vector_view nearest_row = gsl_matrix_row(nearest, i);
    gsl_vector_add_constant(&nearest_row.vector, phys_origin);
  }

  return XLAL_SUCCESS;

}
