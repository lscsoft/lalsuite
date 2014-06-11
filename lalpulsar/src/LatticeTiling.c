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
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <lal/LatticeTiling.h>
#include <lal/LALStdlib.h>
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
/// Lattice tiling status
///
typedef enum tagLT_Status {
  LT_S_INCOMPLETE,
  LT_S_INITIALISED,
  LT_S_STARTED,
  LT_S_FINISHED
} LT_Status;

///
/// Lattice tiling bound info
///
typedef struct tagLT_Bound {
  bool tiled;					///< True if the bound is tiled, false if it is a single point
  LatticeBoundFunction func;			///< Parameter space bound function
  void* data_lower;				///< Arbitrary data describing lower parameter space bound
  void* data_upper;				///< Arbitrary data describing upper parameter space bound
} LT_Bound;

///
/// Lattice tiling nearest point index lookup trie
///
typedef struct tagLT_IndexLookup {
  int32_t min;					///< Minimum integer point value in this dimension
  int32_t max;					///< Maximum integer point value in this dimension
  union {
    struct tagLT_IndexLookup* next;		///< Pointer to array of trie structures for the next-highest dimension
    uint64_t index;				///< Lattice tiling index in the highest dimension
  };
} LT_IndexLookup;

///
/// Lattice tiling state structure
///
struct tagLatticeTiling {
  size_t dimensions;				///< Dimension of the parameter space
  LT_Status status;				///< Status of the tiling
  LT_Bound* bounds;				///< Array of parameter-space bound info for each dimension
  gsl_vector_uint* tiled_idx;			///< Indices of the tiled dimensions of the parameter space
  LatticeType lattice;				///< Type of lattice to generate tiling with
  gsl_vector* phys_bbox;			///< Metric ellipse bounding box in physical coordinates
  gsl_vector* phys_scale;			///< Normalised to physical coordinate scale
  gsl_vector* phys_offset;			///< Normalised to physical coordinate offset
  gsl_matrix* increment;			///< increment matrix of the lattice tiling generator
  gsl_matrix* inv_increment;			///< Inverse of increment matrix of the lattice tiling generator
  gsl_matrix* tiled_increment;			///< Increment matrix of tiled dimensions of the lattice generator
  gsl_vector* point;				///< Current lattice point
  gsl_vector* lower;				///< Current lower bound on parameter space
  gsl_vector* upper;				///< Current upper bound on parameter space
  uint64_t count;				///< Number of points generated so far
  uint64_t total_count;				///< Total number of points in parameter space
  LT_IndexLookup* lookup_base;			///< Lookup trie for finding index of nearest point
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
/// Exchange all rows and columns of the matrix A
///
static void LT_ExchangeRowsCols(gsl_matrix* A) {
  for (size_t i = 0; i < A->size1 / 2; ++i) {
    gsl_matrix_swap_rows(A, i, A->size1 - i - 1);
  }
  for (size_t j = 0; j < A->size2 / 2; ++j) {
    gsl_matrix_swap_columns(A, j, A->size2 - j - 1);
  }
}

///
/// Returns the lower and upper parameter-space bounds
///
static void LT_GetBounds(
  const LatticeTiling* tiling,			///< [in] Tiling state
  const size_t dimension,			///< [in] Dimension on which bound applies
  const gsl_vector* point,			///< [in] Point at which to find bounds
  const bool padded,				///< [in] Whether to add padding to parameter-space bounds
  double* lower,				///< [out] Lower bound on point
  double* upper					///< [out] Upper bound on point
  )
{
  const LT_Bound* bound = &tiling->bounds[dimension];

  // Convert point from normalised to physical coordinates, and get
  // a view of only the first (dimension) number of dimensions; if
  // dimension == 0, NULL is passed to the bound function instead
  double phys_point_array[dimension + 1];
  for (size_t i = 0; i < dimension; ++i) {
    const double phys_scale = gsl_vector_get(tiling->phys_scale, i);
    const double phys_offset = gsl_vector_get(tiling->phys_offset, i);
    phys_point_array[i] = phys_scale * gsl_vector_get(point, i) + phys_offset;
  }
  gsl_vector_view phys_point_view;
  gsl_vector* phys_point = NULL;
  if (dimension > 0) {
    phys_point_view = gsl_vector_view_array(phys_point_array, dimension);
    phys_point = &phys_point_view.vector;
  }

  // Get a view of only the first (dimension+1) number of dimensions
  // of the metric ellipse bounding box in physical coordinates
  gsl_vector_const_view phys_bbox_view = gsl_vector_const_subvector(tiling->phys_bbox, 0, dimension + 1);
  const gsl_vector* phys_bbox = &phys_bbox_view.vector;

  // Set default padding to metric ellipse bounding box in this dimension,
  // for tiled dimensions, otherwise zero for non-tiled dimensions
  double lower_pad = bound->tiled ? gsl_vector_get(tiling->phys_bbox, dimension) : 0.0;
  double upper_pad = lower_pad;

  // Compute lower parameter space bound
  (bound->func)(dimension, phys_point, phys_bbox, bound->data_lower, lower, &lower_pad);

  // Convert lower bound from physical back to normalised coordinates
  const double phys_scale = gsl_vector_get(tiling->phys_scale, dimension);
  const double phys_offset = gsl_vector_get(tiling->phys_offset, dimension);
  *lower = (*lower - phys_offset) / phys_scale;

  // If this dimension is non-tiled, we're done
  if (!bound->tiled) {
    *upper = *lower;
    return;
  }

  // Compute upper parameter space bound
  (bound->func)(dimension, phys_point, phys_bbox, bound->data_upper, upper, &upper_pad);

  // Convert upper bound from physical back to normalised coordinates
  *upper = (*upper - phys_offset) / phys_scale;

  // Optionally add padding (converted from physical back to normalised coordinates)
  if (padded) {
    *lower -= lower_pad / phys_scale;
    *upper += upper_pad / phys_scale;
  }

}

///
/// Destroy a lattice tiling nearest point index lookup trie
///
static void LT_DestroyLookup(
  const size_t ti,				///< [in] Current depth of the trie
  const size_t tn,				///< [in] Total depth of the trie
  LT_IndexLookup* lookup			///< [in] Pointer to array of trie structures
  )
{

  // If at least two dimensions below the highest dimension, call recursively to free memory
  // allocated in higher dimensions. No need to do this for second-highest dimension, since
  // highest dimension does not allocate memory ('index' is used instead of 'next').
  if (ti + 2 < tn) {
    LT_IndexLookup* next = lookup->next;
    for (int32_t i = lookup->min; i <= lookup->max; ++i) {
      LT_DestroyLookup(ti + 1, tn, next++);
    }
  }

  // If at least one dimension below the highest dimension, free memory allocated in this
  // dimension. This will not be called if tn == 1, since then 'next' is never used.
  if (ti + 1 < tn) {
    XLALFree(lookup->next);
  }

}

///
/// Print a lattice tiling nearest point index lookup trie
///
static void LT_PrintLookup(
  const LatticeTiling* tiling,			///< [in] Tiling state
  FILE* fp,					///< [in] File pointer to print trie to
  const size_t ti,				///< [in] Current depth of the trie
  const size_t tn,				///< [in] Total depth of the trie
  const LT_IndexLookup* lookup,			///< [in] Pointer to array of trie structures
  gsl_vector* int_point				///< [in] Temporary point used in printing
  )
{

  // Print indentation
  for (size_t s = 0; s < ti; ++s) {
    fprintf(fp, "   ");
  }

  // If lookup is NULL (which should never happen), print error and return
  if (lookup == NULL) {
    fprintf(fp, "ERROR: lookup is NULL\n");
    return;
  }

  // Set (i)th dimension of point to minimum, then transform from generating integer
  // space, by multiplying by the lattice increments matrix, to get lower bound
  // in normalised coordinates
  gsl_vector_set(int_point, ti, lookup->min);
  gsl_vector_view tiled_increment_row = gsl_matrix_row(tiling->tiled_increment, ti);
  double lower = 0;
  gsl_blas_ddot(int_point, &tiled_increment_row.vector, &lower);

  // Get upper bound by adding increment to lower bound
  const double spacing = gsl_matrix_get(tiling->tiled_increment, ti, ti);
  double upper = lower + (lookup->max - lookup->min + 1) * spacing;

  // Transform lower/upper bounds from normalised to physical coordinates
  double phys_lower = 0, phys_upper = 0;
  {
    const size_t i = gsl_vector_uint_get(tiling->tiled_idx, ti);
    const double phys_scale = gsl_vector_get(tiling->phys_scale, i);
    const double phys_offset = gsl_vector_get(tiling->phys_offset, i);
    phys_lower = lower * phys_scale + phys_offset;
    phys_upper = upper * phys_scale + phys_offset;
  }

  // Print information on the current lookup trie dimension
  fprintf(fp, "dim #%zu/%zu   min %+5i   max %+5i   lower %+10e   upper %+10e",
          ti + 1, tn, lookup->min, lookup->max, phys_lower, phys_upper);

  // If this is the highest dimension, print index and return
  if (ti + 1 == tn) {
    fprintf(fp, "   index %"PRIu64"\n", lookup->index);
  } else {

    // Otherwise, loop over this dimension
    fprintf(fp, "\n");
    LT_IndexLookup* next = lookup->next;
    for (int32_t p = lookup->min; p <= lookup->max; ++p, ++next) {

      // Set (i)th dimension of integer point
      gsl_vector_set(int_point, ti, p);

      // Print higher dimensions
      LT_PrintLookup(tiling, fp, ti + 1, tn, next, int_point);

    }

  }

}

///
/// Find the nearest point, and its index, in the lattice tiling
///
static void LT_FindNearestPoint(
  const size_t n,				///< [in] Number of tiled dimensions
  const LatticeType lattice,			///< [in] Lattice type
  const LT_IndexLookup* lookup_base,		///< [in] Lookup trie for finding index of nearest point
  const gsl_vector* int_point,			///< [in] Point in generating integer space
  long long* nearest_int_point,			///< [in/out] Nearest point in generating integer space
  UINT8* nearest_idx				///< [in/out] Index of nearest point
  )
{

  // Find the nearest point to 'int_point', which is in generating integer space
  switch (lattice) {

  case LATTICE_TYPE_CUBIC:   // Cubic (\f$Z_n\f$) lattice
  {

    // Since Zn lattice generator is the identity, 'int_point' is already in lattice space;
    // to find the nearest point, just need to round to nearest integer
    for (size_t i = 0; i < n; ++i) {
      nearest_int_point[i] = llround(gsl_vector_get(int_point, i));
    }

  }
  break;

  case LATTICE_TYPE_ANSTAR:   // An-star (\f$A_n^*\f$) lattice
  {

    // Transform 'int_point' from n-D generating integer space to (n+1)-D lattice space:
    // The An* lattice generator in n+1 dimensions [see XLALComputeLatticeGenerator()]
    // projects the hypercubic Z(n+1) lattice onto the plane perpendicular to the
    // (n+1)-D 1-vector 'w' = (1,1,...,1), which generates the An* lattice in that plane.
    // The generator matrix is therefore the projection operator
    //   G = Q = I - w*w^T/(n+1)
    // where 'I' is the (n+1)x(n+1) identity, and w*w^T denotes a (n+1)x(n+1) matrix of 1s.
    // This projection is not one-to-one; points in Zn which differ by the 1-vector 'w'
    // will be mapped to the same point in An*. We can therefore choose an n-D subset of
    // points in Zn to map to and from 'int_point'; we choose the subset where the 1st
    // dimension is zero, and therefore the points in Zn are (0,'int_point'). Therefore
    //   w*w^T/(n+1) * (0,'int_point') = w*sum('int_point')/(n+1)
    // and
    //  G * (0,'int_point') = (0,'int_point') - w*sum('int_point')/(n+1)
    //                      = (S, 'int_point' + S)
    // where S = -sum('int_point')/(n+1), and 'int_point' + S denotes adding S to each
    // element of 'int_point'. We first calculate S, store it in y[0], then add S to
    // each element of 'int_point', storing them in y[1]...y[n].
    double y[n+1];
    y[0] = 0;
    for (size_t i = 0; i < n; ++i) {
      y[0] += gsl_vector_get(int_point, i);
    }
    y[0] *= -1.0 / (n + 1);
    for (size_t i = 0; i < n; ++i) {
      y[i+1] = gsl_vector_get(int_point, i) + y[0];
    }

    // Find the nearest point in An* to the point 'y', using the O(n) Algorithm 2 given in:
    //   McKilliam et.al., "A linear-time nearest point algorithm for the lattice An*"
    //   in "International Symposium on Information Theory and Its Applications", ISITA2008,
    //   Auckland, New Zealand, 7-10 Dec. 2008. DOI: 10.1109/ISITA.2008.4895596
    // Notes:
    //   * Since Algorithm 2 uses 1-based arrays, we have to translate, e.g.:
    //       z_t in paper <---> z[tn-1] in C code
    //   * Line 6 in Algorithm 2 as written in the paper is in error, see correction below.
    //   * We are only interested in 'k', the generating integers of the nearest point
    //     'x = Q * k', therefore line 26 in Algorithm 2 is not included.
    long long k[n+1];
    {

      // Lines 1--4, 20
      double z[n+1], alpha = 0, beta = 0;
      size_t bucket[n+1], link[n+1];
      for (size_t i = 1; i <= n + 1; ++i) {
        k[i-1] = llround(y[i-1]);   // Line 20, moved here to avoid duplicate round
        z[i-1] = y[i-1] - k[i-1];
        alpha += z[i-1];
        beta += z[i-1]*z[i-1];
        bucket[i-1] = 0;
      }

      // Lines 5--8. Note correction to line 6:
      //   i = n + 1 - (n + 1)*floor(z_t + 0.5)
      // as written in McKilliam et.al. should instead read
      //   i = n + 1 - floor((n + 1)*(z_t + 0.5))
      for (size_t t = 1; t <= n + 1; ++t) {
        const size_t i = n + 1 - (size_t)floor((n + 1) * (z[t-1] + 0.5));
        link[t-1] = bucket[i-1];
        bucket[i-1] = t;
      }

      // Lines 9--10
      double D = beta - alpha*alpha / (n + 1);
      size_t m = 0;

      // Lines 11--19
      for (size_t i = 1; i <= n + 1; ++i) {
        size_t t = bucket[i-1];
        while (t != 0) {
          alpha = alpha - 1;
          beta = beta - 2*z[t-1] + 1;
          t = link[t-1];
        }
        double d = beta - alpha*alpha / (n + 1);
        if (d < D) {
          D = d;
          m = i;
        }
      }

      // Lines 21--25
      for (size_t i = 1; i <= m; ++i) {
        size_t t = bucket[i-1];
        while (t != 0) {
          k[t-1] = k[t-1] + 1;
          t = link[t-1];
        }
      }

    }

    // Transform 'k' from (n+1)-D lattice space to n-D generating integers 'nearest_int_point':
    // Since points in Zn which differ by the 1-vector 'w' are mapped to the same point in An*,
    // we can freely add or subtract 'w' from 'k' until the first dimension of 'k' is zero, which
    // is our chosen set of points in Zn. Clearly this is achieved by subtracting k[0]*w from 'k',
    // i.e. by substracting 'k[0]' from each element of 'k'.
    for (size_t i = 0; i < n; ++i) {
      nearest_int_point[i] = k[i+1] - k[0];
    }

  }
  break;

  default:
    XLAL_ERROR_VOID(XLAL_EINVAL, "Invalid 'lattice'=%u", lattice);
  }

  // Look up the index of 'int_point' index lookup trie: start at the base of the trie
  const LT_IndexLookup* lookup = lookup_base;

  // Iterate over tiled dimensions, except the highest dimension
  for (size_t i = 0; i + 1 < n; ++i) {

    // Make sure (i)th dimension of nearest integer point is within lookup bounds, then subtract minimum to get index
    const size_t idx_i = GSL_MAX( lookup->min, GSL_MIN( nearest_int_point[i], lookup->max ) ) - lookup->min;

    // Jump to the next level of the lookup trie, based in index
    lookup = &lookup->next[idx_i];

  }

  // Make sure (n-1)th dimension of nearest integer point is within lookup bounds, then subtract minimum to get index
  const size_t idx_nm1 = GSL_MAX( lookup->min, GSL_MIN( nearest_int_point[n-1], lookup->max ) ) - lookup->min;

  // Set index of nearest point
  *nearest_idx = lookup->index + idx_nm1;

}

gsl_vector* XLALMetricEllipseBoundingBox(
  const gsl_matrix* metric,
  const double max_mismatch
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

  // Copy metric, and ensure it is diagonally normalised
  for (size_t i = 0; i < n; ++i) {
    const double norm_i = gsl_matrix_get(metric, i, i);
    for (size_t j = 0; j < n; ++j) {
      const double norm_j = gsl_matrix_get(metric, j, j);
      const double metric_i_j = gsl_matrix_get(metric, i, j);
      gsl_matrix_set(LU_decomp, i, j, metric_i_j / sqrt(norm_i * norm_j));
    }
  }

  // Compute metric inverse
  int LU_sign = 0;
  GCALL_NULL(gsl_linalg_LU_decomp(LU_decomp, LU_perm, &LU_sign));
  GCALL_NULL(gsl_linalg_LU_invert(LU_decomp, LU_perm, inverse));

  // Compute bounding box, and reverse diagonal scaling
  for (size_t i = 0; i < n; ++i) {
    const double norm_i = gsl_matrix_get(metric, i, i);
    const double bounding_box_i = sqrt(norm_i * max_mismatch * gsl_matrix_get(inverse, i, i));
    gsl_vector_set(bounding_box, i, bounding_box_i);
  }

  // Cleanup
  GFMAT(LU_decomp, inverse);
  GFPERM(LU_perm);

  return bounding_box;

}

gsl_matrix* XLALComputeMetricOrthoBasis(
  const gsl_matrix* metric
  )
{

  // Check input
  XLAL_CHECK_NULL(metric->size1 == metric->size2, XLAL_ESIZE);
  const size_t n = metric->size1;

  // Allocate memory
  gsl_matrix* GAMAT_NULL(basis, n, n);
  gsl_vector* GAVEC_NULL(temp, n);

  // Initialise basis to the identity matrix
  gsl_matrix_set_identity(basis);

  // Orthonormalise the columns of 'basis' using numerically stabilised Gram-Schmidt
  for (ssize_t i = n - 1; i >= 0; --i) {
    gsl_vector_view col_i = gsl_matrix_column(basis, i);
    for (ssize_t j = n - 1; j > i; --j) {
      gsl_vector_view col_j = gsl_matrix_column(basis, j);

      // Compute inner product of (j)th and (i)th columns with the metric
      gsl_blas_dgemv(CblasNoTrans, 1.0, metric, &col_j.vector, 0.0, temp);
      double inner_prod = 0.0;
      gsl_blas_ddot(&col_i.vector, temp, &inner_prod);

      // Subtract component of (j)th column from (i)th column
      gsl_vector_memcpy(temp, &col_j.vector);
      gsl_vector_scale(temp, inner_prod);
      gsl_vector_sub(&col_i.vector, temp);

    }

    // Compute inner product of (i)th column with itself
    gsl_blas_dgemv(CblasNoTrans, 1.0, metric, &col_i.vector, 0.0, temp);
    double inner_prod = 0.0;
    gsl_blas_ddot(&col_i.vector, temp, &inner_prod);

    // Normalise (i)th column
    gsl_vector_scale(&col_i.vector, 1.0 / sqrt(inner_prod));

  }

  // Matrix will be lower triangular, so zero out upper triangle
  LT_ZeroStrictUpperTriangle(basis);

  // Cleanup
  GFVEC(temp);

  return basis;

}

gsl_matrix* XLALComputeLatticeGenerator(
  const size_t dimensions,
  const LatticeType lattice,
  const double max_mismatch
  )
{

  // Check input
  XLAL_CHECK_NULL(dimensions > 0, XLAL_EINVAL);
  XLAL_CHECK_NULL(max_mismatch > 0.0, XLAL_EINVAL);
  const size_t n = dimensions;

  // Allocate memory
  gsl_matrix* GAMAT_NULL(generator, n, n);
  gsl_matrix* GAMAT_NULL(LU_decomp, n, n);
  gsl_permutation* GAPERM_NULL(LU_perm, n);

  // Compute lattice generator and normalised thickness
  double norm_thickness = 0.0;
  switch (lattice) {

  case LATTICE_TYPE_CUBIC:   // Cubic (\f$Z_n\f$) lattice
  {

    // Zn lattice generator is the identity
    gsl_matrix_set_identity(generator);

    // Zn normalised thickness
    norm_thickness = pow(sqrt(n)/2.0, n);

  }
  break;

  case LATTICE_TYPE_ANSTAR:   // An-star (\f$A_n^*\f$) lattice
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
    // where Q is an orthogonal matrix and L an lower triangular matrix.
    // This is found using the more commonly implemented QR decomposition by:
    // - exchanging the rows/columns of Gp
    // - decomposing Gp = Qp * Lp, where Lp is upper triangular
    // - exchanging the rows/columns of Qp to give Q
    // - exchanging the rows/columns of Lp to give L
    gsl_matrix_view Gp = gsl_matrix_submatrix(G, 0, 1, n + 1, n);
    LT_ExchangeRowsCols(&Gp.matrix);
    GCALL_NULL(gsl_linalg_QR_decomp(&Gp.matrix, tau));
    GCALL_NULL(gsl_linalg_QR_unpack(&Gp.matrix, tau, Q, L));
    LT_ExchangeRowsCols(Q);
    LT_ExchangeRowsCols(L);

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
    XLAL_ERROR_NULL(XLAL_EINVAL, "Invalid 'lattice'=%u", lattice);
  }

  // Generator will be lower triangular, so zero out upper triangle
  LT_ZeroStrictUpperTriangle(generator);

  // Ensure that the generator has positive diagonal elements, by
  // changing the sign of the columns of the matrix, if necessary
  for (size_t j = 0; j < generator->size2; ++j) {
    gsl_vector_view generator_col = gsl_matrix_column(generator, j);
    XLAL_CHECK_NULL( gsl_vector_get(&generator_col.vector, j) != 0, XLAL_ERANGE, "generator(%zu,%zu) == 0", j ,j );
    if (gsl_vector_get(&generator_col.vector, j) < 0) {
      gsl_vector_scale(&generator_col.vector, -1);
    }
  }

  // Compute generator LU decomposition
  gsl_matrix_memcpy(LU_decomp, generator);
  int LU_sign = 0;
  GCALL_NULL(gsl_linalg_LU_decomp(LU_decomp, LU_perm, &LU_sign));

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

LatticeTiling* XLALCreateLatticeTiling(
  const size_t dimensions
  )
{

  // Check input
  XLAL_CHECK_NULL(dimensions > 0, XLAL_EINVAL);
  const size_t n = dimensions;

  // Allocate and initialise tiling structure
  LatticeTiling* tiling = XLALCalloc(1, sizeof(*tiling));
  XLAL_CHECK_NULL(tiling != NULL, XLAL_ENOMEM);
  tiling->dimensions = n;
  tiling->status = LT_S_INCOMPLETE;
  tiling->lattice = LATTICE_TYPE_MAX;

  // Allocate tiling structure memory
  tiling->bounds = XLALCalloc(n, sizeof(*tiling->bounds));
  XLAL_CHECK_NULL(tiling->bounds != NULL, XLAL_ENOMEM);
  GAVEC_NULL(tiling->phys_bbox, n);
  GAVEC_NULL(tiling->phys_scale, n);
  GAVEC_NULL(tiling->phys_offset, n);
  GAVEC_NULL(tiling->point, n);
  GAVEC_NULL(tiling->lower, n);
  GAVEC_NULL(tiling->upper, n);

  return tiling;

}

void XLALDestroyLatticeTiling(
  LatticeTiling* tiling
  )
{

  if (tiling) {

    const size_t n = tiling->dimensions;
    const size_t tn = (tiling->tiled_idx != NULL) ? tiling->tiled_idx->size : 0;

    // Free bounds data
    if (tiling->bounds != NULL) {
      for (size_t i = 0; i < n; ++i) {
        XLALFree(tiling->bounds[i].data_lower);
        XLALFree(tiling->bounds[i].data_upper);
      }
      XLALFree(tiling->bounds);
    }

    // Free lookup trie
    if (tiling->lookup_base != NULL) {
      LT_DestroyLookup(0, tn, tiling->lookup_base);
      XLALFree(tiling->lookup_base);
    }

    // Free vectors and matrices
    GFMAT(tiling->increment, tiling->inv_increment, tiling->tiled_increment);
    GFVEC(tiling->phys_bbox, tiling->phys_scale, tiling->phys_offset, tiling->point, tiling->lower, tiling->upper);
    GFVECU(tiling->tiled_idx);

    // Free tiling structure
    XLALFree(tiling);

  }

}

size_t XLALGetLatticeTotalDimensions(
  const LatticeTiling* tiling
  )
{

  // Check input
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);

  return tiling->dimensions;

}

size_t XLALGetLatticeTiledDimensions(
  const LatticeTiling* tiling
  )
{

  // Check input
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_VAL(0, tiling->status > LT_S_INCOMPLETE, XLAL_EINVAL);

  return (tiling->tiled_idx != NULL) ? tiling->tiled_idx->size : 0;

}

uint64_t XLALGetLatticePointCount(
  const LatticeTiling* tiling
  )
{

  // Check input
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_VAL(0, tiling->status >= LT_S_INITIALISED, XLAL_EINVAL);

  return tiling->count;

}

uint64_t XLALCountLatticePoints(
  LatticeTiling* tiling
  )
{

  // Check input
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_VAL(0, tiling->status > LT_S_INCOMPLETE, XLAL_EINVAL);

  // If points have not already been counted, count them
  if (tiling->total_count == 0) {

    // Iterate over all points
    while (XLALNextLatticePoint(tiling, NULL) >= 0) {
      XLAL_CHECK_VAL(0, XLALFastForwardLatticeTiling(tiling, NULL, NULL) == XLAL_SUCCESS, XLAL_EFUNC);
    }
    XLAL_CHECK_VAL(0, xlalErrno == 0, XLAL_EFUNC);
    XLAL_CHECK_VAL(0, tiling->total_count > 0, XLAL_EFAILED);

    // Restart tiling
    XLAL_CHECK_VAL(0, XLALRestartLatticeTiling(tiling) == XLAL_SUCCESS, XLAL_EFUNC);

  }

  // Return the point count
  return tiling->total_count;

}

int XLALSetLatticeBound(
  LatticeTiling* tiling,
  const size_t dimension,
  const LatticeBoundFunction func,
  const size_t data_len,
  void* data_lower,
  void* data_upper
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status == LT_S_INCOMPLETE, XLAL_EINVAL);
  XLAL_CHECK(dimension < tiling->dimensions, XLAL_ESIZE);
  XLAL_CHECK(func != NULL, XLAL_EFAULT);
  XLAL_CHECK(data_len > 0, XLAL_EFAULT);
  XLAL_CHECK(data_lower != NULL, XLAL_EFAULT);
  XLAL_CHECK(data_upper != NULL, XLAL_EFAULT);

  // Determine if this dimension is tiled
  const bool tiled = (memcmp(data_lower, data_upper, data_len) != 0);

  // Set the parameter-space bound
  tiling->bounds[dimension].tiled = tiled;
  tiling->bounds[dimension].func = func;
  tiling->bounds[dimension].data_lower = data_lower;
  tiling->bounds[dimension].data_upper = data_upper;

  return XLAL_SUCCESS;

}

static void ConstantBound(
  const size_t dimension UNUSED,
  const gsl_vector* point UNUSED,
  const gsl_vector* bbox UNUSED,
  const void* data UNUSED,
  double* bound,
  double* padding UNUSED
  )
{

  // Set bound
  *bound = *((const double*) data);

}

int XLALSetLatticeConstantBound(
  LatticeTiling* tiling,
  const size_t dimension,
  const double bound1,
  const double bound2
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status == LT_S_INCOMPLETE, XLAL_EINVAL);
  XLAL_CHECK(isfinite(bound1), XLAL_EINVAL);
  XLAL_CHECK(isfinite(bound2), XLAL_EINVAL);

  // Allocate memory
  double* data_lower = XLALMalloc(sizeof(*data_lower));
  XLAL_CHECK(data_lower != NULL, XLAL_ENOMEM);
  double* data_upper = XLALMalloc(sizeof(*data_lower));
  XLAL_CHECK(data_upper != NULL, XLAL_ENOMEM);

  // Set the parameter-space bound
  *data_lower = GSL_MIN(bound1, bound2);
  *data_upper = GSL_MAX(bound1, bound2);
  XLAL_CHECK(XLALSetLatticeBound(tiling, dimension, ConstantBound,
                                 sizeof(*data_lower), data_lower, data_upper) == XLAL_SUCCESS, XLAL_EFUNC);

  return XLAL_SUCCESS;

}

int XLALSetLatticeTypeAndMetric(
  LatticeTiling* tiling,
  const LatticeType lattice,
  const gsl_matrix* metric,
  const double max_mismatch
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status == LT_S_INCOMPLETE, XLAL_EINVAL);
  const size_t n = tiling->dimensions;
  XLAL_CHECK(metric != NULL, XLAL_EFAULT);
  XLAL_CHECK(metric->size1 == n && metric->size2 == n, XLAL_EINVAL);
  XLAL_CHECK(max_mismatch > 0, XLAL_EINVAL);

  // Save the type of lattice to generate tiling with
  tiling->lattice = lattice;

  // Check that all parameter-space dimensions are bounded, and record indices of tiled dimensions
  size_t tn = 0;
  for (size_t i = 0; i < tiling->dimensions; ++i) {
    XLAL_CHECK(tiling->bounds[i].func != NULL, XLAL_EFAILED, "Dimension #%i is unbounded", i);
    if (tiling->bounds[i].tiled) {
      ++tn;
    }
  }
  if (tn > 0) {
    GAVECU(tiling->tiled_idx, tn);
    for (unsigned int i = 0, ti = 0; i < tiling->dimensions; ++i) {
      if (tiling->bounds[i].tiled) {
        gsl_vector_uint_set(tiling->tiled_idx, ti, i);
        ++ti;
      }
    }
  }

  // Check metric is symmetric and has positive diagonal elements
  for (size_t i = 0; i < n; ++i) {
    XLAL_CHECK(gsl_matrix_get(metric, i, i) > 0, XLAL_EINVAL, "metric(%zu,%zu) <= 0", i, i);
    for (size_t j = i + 1; j < n; ++j) {
      XLAL_CHECK(gsl_matrix_get(metric, i, j) == gsl_matrix_get(metric, j, i), XLAL_EINVAL, "metric(%zu,%zu) != metric(%zu,%zu)", i, j, j, i);
    }
  }

  // Set physical parameter-space offset from parameter-space bounds
  gsl_vector_set_zero(tiling->phys_bbox);
  gsl_vector_set_all(tiling->phys_scale, 1.0);
  gsl_vector_set_zero(tiling->phys_offset);
  for (size_t i = 0; i < n; ++i) {
    double phys_lower_i = 0, phys_upper_i = 0;
    LT_GetBounds(tiling, i, tiling->phys_offset, false, &phys_lower_i, &phys_upper_i);
    gsl_vector_set(tiling->phys_offset, i, 0.5*(phys_lower_i + phys_upper_i));
  }

  // Calculate physical parameter-space scale from metric diagonal elements
  for (size_t ti = 0; ti < tn; ++ti) {
    const size_t i = gsl_vector_uint_get(tiling->tiled_idx, ti);
    const double metric_i_i = gsl_matrix_get(metric, i, i);
    gsl_vector_set(tiling->phys_scale, i, 1.0 / sqrt(metric_i_i));
  }

  // If there are tiled dimensions:
  if (tn > 0) {

    // Allocate memory
    GAMAT(tiling->increment, n, n);
    GAMAT(tiling->tiled_increment, tn, tn);
    GAMAT(tiling->inv_increment, n, n);
    gsl_matrix* GAMAT(tiled_metric, tn, tn);
    gsl_matrix* GAMAT(tiled_metric_copy, tn, tn);
    gsl_matrix* GAMAT(inv_tiled_increment, tn, tn);

    // Copy and normalise tiled dimensions of metric
    for (size_t ti = 0; ti < tn; ++ti) {
      const size_t i = gsl_vector_uint_get(tiling->tiled_idx, ti);
      const double phys_scale_i = gsl_vector_get(tiling->phys_scale, i);
      for (size_t tj = 0; tj < tn; ++tj) {
        const size_t j = gsl_vector_uint_get(tiling->tiled_idx, tj);
        const double phys_scale_j = gsl_vector_get(tiling->phys_scale, j);
        const double metric_i_j = gsl_matrix_get(metric, i, j);
        gsl_matrix_set(tiled_metric, ti, tj, metric_i_j * phys_scale_i * phys_scale_j);
      }
    }

    // Check tiled metric is positive definite, by trying to compute its Cholesky decomposition
    gsl_matrix_memcpy(tiled_metric_copy, tiled_metric);
    GCALL(gsl_linalg_cholesky_decomp(tiled_metric_copy), "tiled metric is not positive definite");

    // Compute metric ellipse bounding box
    gsl_vector* tiled_bbox = XLALMetricEllipseBoundingBox(tiled_metric, max_mismatch);
    XLAL_CHECK(tiled_bbox != NULL, XLAL_EFUNC);

    // Copy bounding box in physical coordinates to tiled dimensions
    for (size_t ti = 0; ti < tn; ++ti) {
      const size_t i = gsl_vector_uint_get(tiling->tiled_idx, ti);
      const double phys_scale_i = gsl_vector_get(tiling->phys_scale, i);
      gsl_vector_set(tiling->phys_bbox, i, gsl_vector_get(tiled_bbox, ti) * phys_scale_i);
    }

    // Compute a lower triangular basis matrix whose columns are orthonormal with respect to the tiled metric
    gsl_matrix* tiled_basis = XLALComputeMetricOrthoBasis(tiled_metric);
    XLAL_CHECK(tiled_basis != NULL, XLAL_EFUNC);

    // Compute a lower triangular generator matrix for a given lattice type and mismatch
    gsl_matrix* tiled_generator = XLALComputeLatticeGenerator(tn, lattice, max_mismatch);
    XLAL_CHECK(tiled_generator != NULL, XLAL_EFUNC);

    // Compute the increment matrix, which is the lattice generator expressed in the orthonormal metric basis
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tiled_basis, tiled_generator, 0.0, tiling->tiled_increment);

    // Increment matrix will be lower triangular, so zero out upper triangle
    LT_ZeroStrictUpperTriangle(tiling->tiled_increment);

    // Calculate inverse of lattice increment matrix; set 'inv_tiled_increment'
    // to identity matrix, then multiply by inverse of lattice increment matrix
    gsl_matrix_set_identity(inv_tiled_increment);
    gsl_blas_dtrsm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, tiling->tiled_increment, inv_tiled_increment);

    // Copy lower triangular part of increment matrix to tiled dimensions (strict upper triangle is zero)
    gsl_matrix_set_zero(tiling->increment);
    gsl_matrix_set_zero(tiling->inv_increment);
    for (size_t ti = 0; ti < tn; ++ti) {
      const size_t i = gsl_vector_uint_get(tiling->tiled_idx, ti);
      for (size_t tj = 0; tj <= ti; ++tj) {
        const size_t j = gsl_vector_uint_get(tiling->tiled_idx, tj);
        gsl_matrix_set(tiling->increment, i, j, gsl_matrix_get(tiling->tiled_increment, ti, tj));
        gsl_matrix_set(tiling->inv_increment, i, j, gsl_matrix_get(inv_tiled_increment, ti, tj));
      }
    }

    // Cleanup
    GFMAT(inv_tiled_increment, tiled_basis, tiled_generator, tiled_metric, tiled_metric_copy);
    GFVEC(tiled_bbox);

  }

  // Tiling has been fully initialised
  tiling->status = LT_S_INITIALISED;
  XLAL_CHECK(XLALRestartLatticeTiling(tiling) == XLAL_SUCCESS, XLAL_EFUNC);

  return XLAL_SUCCESS;

}

int XLALNextLatticePoint(
  LatticeTiling* tiling,
  gsl_vector* curr_point
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status > LT_S_INCOMPLETE, XLAL_EINVAL);
  const size_t n = tiling->dimensions;
  const size_t tn = (tiling->tiled_idx != NULL) ? tiling->tiled_idx->size : 0;
  XLAL_CHECK(curr_point == NULL || curr_point->size == n, XLAL_EINVAL);

  // If finished status, nothing more to be done!
  if (tiling->status == LT_S_FINISHED) {
    return -1;
  }

  size_t ti;

  // If initialised status, set and return starting point
  if (tiling->status == LT_S_INITIALISED) {

    // Set parameter-space bounds and starting point
    gsl_vector_set_zero(tiling->point);
    for (size_t i = 0; i < n; ++i) {

      // Get normalised bounds, with padding
      double lower_i = 0, upper_i = 0;
      LT_GetBounds(tiling, i, tiling->point, true, &lower_i, &upper_i);

      // Set parameter-space bounds
      gsl_vector_set(tiling->lower, i, lower_i);
      gsl_vector_set(tiling->upper, i, upper_i);

      // Set starting point to lower bound
      gsl_vector_set(tiling->point, i, lower_i);

      if (tiling->bounds[i].tiled) {

        // Transform point from normalised coordinates to generating integer space,
        // by multiplying by the inverse of the lattice increments matrix
        gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, tiling->inv_increment, tiling->point);

        // Round starting point up to the nearest integer
        gsl_vector_set(tiling->point, i, ceil(gsl_vector_get(tiling->point, i)));

        // Transform point from generating integer space back to normalised coordinates,
        // by multiplying by the lattice increments matrix
        gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, tiling->increment, tiling->point);

      }

    }

    // Initialise count
    tiling->count = 1;

    // Tiling has been started
    tiling->status = LT_S_STARTED;

    // All dimensions of point have changed
    ti = 0;

  } else {

    // Otherwise started status: loop until the next point is found
    ti = tn;
    while (true) {

      // If dimension index is now zero, we're done!
      if (ti == 0) {

        // Tiling is now finished
        tiling->status = LT_S_FINISHED;

        // Save total point count
        tiling->total_count = tiling->count;

        return -1;

      }

      // Decrement current dimension index
      --ti;

      // Get index of tiled dimension
      const size_t i = gsl_vector_uint_get(tiling->tiled_idx, ti);

      // Get increment vector
      gsl_vector_view increment = gsl_matrix_column(tiling->increment, i);

      // Increment current point along index
      gsl_vector_add(tiling->point, &increment.vector);

      // If point is not out of bounds, we have found a point
      const double point_i = gsl_vector_get(tiling->point, i);
      const double upper_i = gsl_vector_get(tiling->upper, i);
      if (point_i <= upper_i) {
        break;
      }

      // Move on to lower dimensions
      continue;

    }

    // Return point to lower bound in higher dimensions
    for (size_t tj = ti + 1; tj < tn; ++tj) {

      // Get index of tiled dimension
      const size_t j = gsl_vector_uint_get(tiling->tiled_idx, tj);

      // Get normalised bounds, with padding
      double lower_j = 0, upper_j = 0;
      LT_GetBounds(tiling, j, tiling->point, true, &lower_j, &upper_j);

      // Set parameter-space bounds
      gsl_vector_set(tiling->lower, j, lower_j);
      gsl_vector_set(tiling->upper, j, upper_j);

      // Get increment vector
      gsl_vector_view increment = gsl_matrix_column(tiling->increment, j);

      // Calculate the distance from current point to the lower bound, in integer number of increments
      const double point_j = gsl_vector_get(tiling->point, j);
      const double dist = ceil( (lower_j - point_j) / gsl_vector_get(&increment.vector, j) );

      // Move point back to lower bound
      gsl_blas_daxpy(dist, &increment.vector, tiling->point);

    }

    // Point was found, so increase count
    ++tiling->count;

  }

  // Optionally, copy current point and transform from normalised to physical coordinates
  if (curr_point != NULL) {
    gsl_vector_memcpy(curr_point, tiling->point);
    gsl_vector_mul(curr_point, tiling->phys_scale);
    gsl_vector_add(curr_point, tiling->phys_offset);
  }

  // Return lowest dimension where point has changed
  return ti;

}

int XLALFastForwardLatticeTiling(
  LatticeTiling* tiling,
  uint32_t* point_count,
  double* point_spacing
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status == LT_S_STARTED, XLAL_EINVAL);
  const size_t tn = (tiling->tiled_idx != NULL) ? tiling->tiled_idx->size : 0;

  // If no tiled dimensions, nothing to fast-forward
  uint32_t count_i = 0;
  double spacing_i = 0.0;

  // If there are tiled dimensions:
  if (tn > 0) {

    // Get index of highest tiled dimension
    const size_t i = gsl_vector_uint_get(tiling->tiled_idx, tn - 1);

    // Get point, spacing, and upper bound in this dimension
    const double point_i = gsl_vector_get(tiling->point, i);
    spacing_i = gsl_matrix_get(tiling->increment, i, i);
    const double upper_i = gsl_vector_get(tiling->upper, i);

    // Calculate number of points to fast-forward, so that then calling
    // XLALNextLatticePoint() will advance the next highest tiled dimension
    count_i = (uint32_t)floor((upper_i - point_i) / spacing_i);

    // Fast-forward over dimension
    gsl_vector_set(tiling->point, i, point_i + count_i * spacing_i);
    tiling->count += count_i;

  }

  // Return count and spacing
  if (point_count != NULL) {
    *point_count = count_i;
  }
  if (point_spacing != NULL) {
    *point_spacing = spacing_i;
  }

  return XLAL_SUCCESS;

}

int XLALRestartLatticeTiling(
  LatticeTiling* tiling
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status > LT_S_INCOMPLETE, XLAL_EINVAL);

  // Restart tiling
  tiling->status = LT_S_INITIALISED;
  tiling->count = 0;

  return XLAL_SUCCESS;

}

int XLALRandomLatticePoints(
  const LatticeTiling* tiling,
  RandomParams* rng,
  gsl_matrix* random_points
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status > LT_S_INCOMPLETE, XLAL_EFAILED);
  const size_t n = tiling->dimensions;

  // Check input
  XLAL_CHECK(rng != NULL, XLAL_EFAULT);
  XLAL_CHECK(random_points != NULL, XLAL_EFAULT);
  XLAL_CHECK(random_points->size1 == n, XLAL_ESIZE);

  // Create random points in lattice tiling parameter space, in normalised coordinates
  for (size_t k = 0; k < random_points->size2; ++k) {
    gsl_vector_view point = gsl_matrix_column(random_points, k);
    for (size_t i = 0; i < n; ++i) {

      // Get normalised bounds
      double lower_i = 0, upper_i = 0;
      LT_GetBounds(tiling, i, &point.vector, false, &lower_i, &upper_i);

      // Generate random number
      const double u = XLALUniformDeviate(rng);

      // Set parameter space point
      gsl_vector_set(&point.vector, i, lower_i + u*(upper_i - lower_i));

    }

  }

  // Transform given points from normalised to physical coordinates
  for (size_t i = 0; i < n; ++i) {
    const double phys_scale = gsl_vector_get(tiling->phys_scale, i);
    const double phys_offset = gsl_vector_get(tiling->phys_offset, i);
    gsl_vector_view row = gsl_matrix_row(random_points, i);
    gsl_vector_scale(&row.vector, phys_scale);
    gsl_vector_add_constant(&row.vector, phys_offset);
  }

  return XLAL_SUCCESS;

}

int XLALBuildLatticeIndexLookup(
  LatticeTiling* tiling
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status == LT_S_INITIALISED, XLAL_EFAILED);
  const size_t n = tiling->dimensions;
  const size_t tn = (tiling->tiled_idx != NULL) ? tiling->tiled_idx->size : 0;

  // If lookup has already been built, or there are no tiled dimensions, we're done!
  if (tiling->lookup_base != NULL || tn == 0) {
    return XLAL_SUCCESS;
  }

  // Allocate vector
  gsl_vector* GAVEC(int_point, n);

  // Initialise pointer to the base lookup trie struct
  LT_IndexLookup* lookup_base = NULL;

  // Allocate array of pointers to the next lookup trie struct
  LT_IndexLookup* next[tn];
  memset(next, 0, sizeof(next));

  // Ensure tiling is restarted
  XLAL_CHECK(XLALRestartLatticeTiling(tiling) == XLAL_SUCCESS, XLAL_EFUNC);

  // Iterate over all points; while tiling is in progress, XLALNextLatticePoint()
  // returns the index of the lowest dimension where the current point has changed
  int ti;
  while ( (ti = XLALNextLatticePoint(tiling, NULL)) >= 0 ) {

    // Transform current point from normalised coordinates to generating integer space,
    // by multiplying by the inverse of the lattice increments matrix.
    gsl_vector_memcpy(int_point, tiling->point);
    gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, tiling->inv_increment, int_point);

    // Iterate over all dimensions where the current point has changed
    for (size_t tj = ti; tj < tn; ++tj) {

      // Get index of tiled dimension
      const size_t j = gsl_vector_uint_get(tiling->tiled_idx, tj);

      // If next lookup struct pointer is NULL, it needs to be initialised
      if (next[tj] == NULL) {

        // Get a pointer to the lookup struct which needs to be built:
        // - if tj is non-zero, we should use the struct pointed to by 'next' in the lower dimension
        // - otherwise, this is the first point of the tiling, so initialise the base lookup struct
        LT_IndexLookup* lookup = NULL;
        if (tj > 0) {
          lookup = next[tj - 1];
        } else {
          lookup_base = lookup = XLALCalloc(1, sizeof(*lookup));
          XLAL_CHECK( lookup_base != NULL, XLAL_ENOMEM );
        }

        // Get point, spacing, and upper bound in (j)th dimension
        const double point_j = gsl_vector_get(tiling->point, j);
        const double spacing_j = gsl_matrix_get(tiling->increment, j, j);
        const double upper_j = gsl_vector_get(tiling->upper, j);

        // Calculate the number of points at this point in (j)th dimension
        const long long num_points_j = 1 + (long long)floor((upper_j - point_j) / spacing_j);

        // Set the minimum and maximum integer point values
        // - Note that the llround() here is important! Even though the elements of 'int_point'
        //   should be very close to integers, they are not guaranteed to be exactly integers,
        //   so llround() is needed to make this safe against numerical errors.
        const long long min_j = llround(gsl_vector_get(int_point, j));
        const long long max_j = min_j + num_points_j - 1;
        XLAL_CHECK( INT32_MIN < min_j && max_j < INT32_MAX, XLAL_EDOM,
                    "Integer point range [%i,%i] is outside int32_t range", min_j, max_j );
        lookup->min = min_j;
        lookup->max = max_j;

        // If we are below the highest dimension:
        if (tj + 1 < tn) {

          // Allocate a new array of lookup structs for the next highest dimension
          lookup->next = XLALCalloc(num_points_j, sizeof(*lookup->next));
          XLAL_CHECK( lookup->next != NULL, XLAL_ENOMEM );

          // Point 'next' to this array in this dimension, for higher dimensions to use
          next[tj] = lookup->next;

        } else {

          // In the highest dimension, we instead saves an index to the current point
          lookup->index = tiling->count - 1;

        }

      } else {

        // Otherwise, advance to the next lookup struct in the array
        ++next[tj];

      }

      // If we are below the highest dimension, set 'next' in the next highest
      // dimension to NULL, so that on the next loop a new array will be created
      if (tj + 1 < tn) {
        next[tj + 1] = NULL;
      }

    }

    // Fast-forward the lattice tiling through the highest dimension
    XLAL_CHECK(XLALFastForwardLatticeTiling(tiling, NULL, NULL) == XLAL_SUCCESS, XLAL_EFUNC);

  }

  // Restart tiling
  XLAL_CHECK(XLALRestartLatticeTiling(tiling) == XLAL_SUCCESS, XLAL_EFUNC);

  // Save lookup trie
  tiling->lookup_base = lookup_base;

  // Cleanup
  GFVEC(int_point);

  return XLAL_SUCCESS;

}

int XLALPrintLatticeIndexLookup(
  const LatticeTiling* tiling,
  FILE* file
  )
{


  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->lookup_base != NULL, XLAL_EFAULT);
  const size_t tn = (tiling->tiled_idx != NULL) ? tiling->tiled_idx->size : 0;

  // Check input
  XLAL_CHECK(file != NULL, XLAL_EFAULT);

  // Allocate vector
  gsl_vector* GAVEC(int_point, tn);

  // Print lookup trie
  LT_PrintLookup(tiling, file, 0, tn, tiling->lookup_base, int_point);

  // Cleanup
  GFVEC(int_point);

  return XLAL_SUCCESS;

}

int XLALNearestLatticePoints(
  const LatticeTiling* tiling,
  const gsl_matrix* points,
  gsl_matrix* nearest_points,
  UINT8Vector* nearest_indices
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status > LT_S_INCOMPLETE, XLAL_EFAILED);
  const size_t n = tiling->dimensions;
  const size_t tn = (tiling->tiled_idx != NULL) ? tiling->tiled_idx->size : 0;
  XLAL_CHECK(nearest_indices == NULL || tn == 0 || tiling->lookup_base != NULL, XLAL_EFAILED,
             "Need to first run XLALBuildLatticeIndexLookup() to build index lookup" );

  // Check input
  XLAL_CHECK(points != NULL, XLAL_EFAULT);
  XLAL_CHECK(points->size1 == n, XLAL_ESIZE);
  const size_t num_points = points->size2;
  XLAL_CHECK(nearest_points != NULL, XLAL_EFAULT);
  XLAL_CHECK(nearest_points->size1 == n, XLAL_ESIZE);
  XLAL_CHECK(nearest_points->size2 >= num_points, XLAL_ESIZE);
  XLAL_CHECK(nearest_indices == NULL || nearest_indices->length >= num_points, XLAL_ESIZE);
  XLAL_CHECK(nearest_indices == NULL || nearest_indices->data != NULL, XLAL_EFAULT);

  // If there are no tiled dimensions:
  if (tn == 0) {

    // Set all columns of 'nearest_points' to the sole point in the tiling,
    // which is just the physical offset given by 'phys_offset'
    for (size_t i = 0; i < n; ++i) {
      const double phys_offset = gsl_vector_get(tiling->phys_offset, i);
      gsl_vector_view nearest_points_row = gsl_matrix_row(nearest_points, i);
      gsl_vector_set_all(&nearest_points_row.vector, phys_offset);
    }

    // Set all elements of 'nearest_indices' to zero, since there is only one point
    if (nearest_indices != NULL) {
      memset(nearest_indices->data, 0, nearest_indices->length * sizeof(nearest_indices->data[0]));
    }

    return XLAL_SUCCESS;

  }

  // Create a workspace view of 'nearest_points' which has the same size as 'points'
  gsl_matrix_view wksp = gsl_matrix_submatrix(nearest_points, 0, 0, tn, num_points);

  // Copy tiled dimensions of 'points' to 'wksp'
  for (size_t ti = 0; ti < tn; ++ti) {
    const size_t i = gsl_vector_uint_get(tiling->tiled_idx, ti);
    gsl_vector_const_view points_row = gsl_matrix_const_row(points, i);
    gsl_vector_view wksp_row = gsl_matrix_row(&wksp.matrix, ti);
    gsl_vector_memcpy(&wksp_row.vector, &points_row.vector);
  }

  // Transform 'wksp' from physical to normalised coordinates
  for (size_t ti = 0; ti < tn; ++ti) {
    const size_t i = gsl_vector_uint_get(tiling->tiled_idx, ti);
    const double phys_scale = gsl_vector_get(tiling->phys_scale, i);
    const double phys_offset = gsl_vector_get(tiling->phys_offset, i);
    gsl_vector_view wksp_row = gsl_matrix_row(&wksp.matrix, ti);
    gsl_vector_add_constant(&wksp_row.vector, -phys_offset);
    gsl_vector_scale(&wksp_row.vector, 1.0/phys_scale);
  }

  // Transform 'wksp' points from normalised coordinates to generating integer space, by
  // multiplying by the inverse of the lattice increments matrix
  gsl_blas_dtrsm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, tiling->tiled_increment, &wksp.matrix);

  // Find the nearest points in the lattice tiling to the 'wksp' points
  long long nearest_int_point[tn];
  for (size_t j = 0; j < num_points; ++j) {
    gsl_vector_view p = gsl_matrix_column(&wksp.matrix, j);
    LT_FindNearestPoint(tn, tiling->lattice, tiling->lookup_base, &p.vector, nearest_int_point, &(nearest_indices->data[j]));
    for (size_t ti = 0; ti < tn; ++ti) {
      gsl_vector_set(&p.vector, ti, nearest_int_point[ti]);
    }
  }
  XLAL_CHECK( xlalErrno == 0, XLAL_EFAILED, "LT_FindNearestPoint() failed" );

  // Transform 'wksp' points from generating integer space back to normalised coordinates,
  // by multiplying by the lattice increments matrix
  gsl_blas_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, tiling->tiled_increment, &wksp.matrix);

  // Transform 'wksp' from normalised back to physical coordinates
  for (size_t ti = 0; ti < tn; ++ti) {
    const size_t i = gsl_vector_uint_get(tiling->tiled_idx, ti);
    const double phys_scale = gsl_vector_get(tiling->phys_scale, i);
    const double phys_offset = gsl_vector_get(tiling->phys_offset, i);
    gsl_vector_view wksp_row = gsl_matrix_row(&wksp.matrix, ti);
    gsl_vector_scale(&wksp_row.vector, phys_scale);
    gsl_vector_add_constant(&wksp_row.vector, phys_offset);
  }

  // Copy 'wksp' to tiled dimensions of 'nearest_points'
  // - because 'wksp' points to same memory as 'nearest_points', we don't need to copy
  //   if n == tn; if we do need to copy, we must iterate in reverse from (tn - 1) to 0,
  //   in order not to erase results already stored in 'nearest_points'
  if (tn != n) {
    for (ssize_t ti = tn - 1; ti >= 0; --ti) {
      const size_t i = gsl_vector_uint_get(tiling->tiled_idx, ti);
      gsl_vector_view wksp_row = gsl_matrix_row(&wksp.matrix, ti);
      gsl_vector_view nearest_points_row = gsl_matrix_row(nearest_points, i);
      gsl_vector_memcpy(&nearest_points_row.vector, &wksp_row.vector);
    }
  }

  // Set any non-tiled dimensions in 'nearest_points' to the sole point in that dimension
  for (size_t i = 0; i < n; ++i) {
    const double phys_scale = gsl_vector_get(tiling->phys_scale, i);
    const double phys_offset = gsl_vector_get(tiling->phys_offset, i);
    if (!tiling->bounds[i].tiled) {
      for (size_t j = 0; j < num_points; ++j) {

        // Get normalised bounds
        gsl_vector_view point_j = gsl_matrix_column(nearest_points, j);
        double lower_i = 0, upper_i = 0;
        LT_GetBounds(tiling, i, &point_j.vector, false, &lower_i, &upper_i);

        // Set point
        gsl_vector_set(&point_j.vector, i, lower_i*phys_scale + phys_offset);

      }
    }
  }

  return XLAL_SUCCESS;

}
