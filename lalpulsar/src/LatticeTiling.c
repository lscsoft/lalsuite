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
  int int_lower;				///< Lower integer point bound in this dimension
  int int_upper;				///< Upper integer point bound in this dimension
  union {
    struct tagLT_IndexLookup* next;		///< Pointer to array of trie structures for the next-highest dimension
    UINT8 index;				///< Lattice tiling index in the highest dimension
  };
} LT_IndexLookup;

///
/// Lattice tiling state structure
///
struct tagLatticeTiling {
  size_t dims;					///< Number of dimensions of the parameter space
  gsl_vector_uint* tiled;			///< Tiled dimensions of the parameter space
  LT_Status status;				///< Status of the tiling
  LT_Bound* bounds;				///< Array of parameter-space bound info for each dimension
  LatticeType lattice;				///< Type of lattice to generate tiling with
  gsl_vector* phys_bbox;			///< Metric ellipse bounding box in physical coordinates
  gsl_vector* phys_offset;			///< Physical coordinate offset
  gsl_matrix* int_from_phys;			///< Transform to generating integers from physical coordinates
  gsl_matrix* phys_from_int;			///< Transform to physical coordinates from generating integers
  gsl_vector_int* int_point;			///< Current lattice point in generating integers
  gsl_vector* phys_point;			///< Current lattice point in physical coordinates
  gsl_vector_int* int_lower;			///< Current lower parameter-space bound in generating integers
  gsl_vector_int* int_upper;			///< Current upper parameter-space bound in generating integers
  UINT8 count;					///< Number of points generated so far
  UINT8 total_count;				///< Total number of points in parameter space
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
static void LT_GetPhysBounds(
  const LatticeTiling* tiling,			///< [in] Tiling state
  const size_t dimension,			///< [in] Dimension on which bound applies
  const gsl_vector* phys_point,			///< [in] Physical point at which to find bounds
  const bool padded,				///< [in] Whether to add padding to parameter-space bounds (ignored if not tiled)
  double* phys_lower,				///< [out] Physical lower bound on parameter space
  double* phys_upper				///< [out] Physical upper bound on parameter space
  )
{
  const LT_Bound* bound = &tiling->bounds[dimension];

  // Get a view of only the first (dimension) number of dimensions;
  // if dimension == 0, NULL is passed to the bound function instead
  gsl_vector_const_view phys_point_dim_view = gsl_vector_const_subvector(phys_point, 0, (dimension == 0) ? 1 : dimension);
  const gsl_vector* phys_point_dim = (dimension == 0) ? NULL : &phys_point_dim_view.vector;

  // Get a view of only the first (dimension+1) number of dimensions
  // of the metric ellipse bounding box in physical coordinates
  gsl_vector_const_view phys_bbox_dim_view = gsl_vector_const_subvector(tiling->phys_bbox, 0, dimension + 1);
  const gsl_vector* phys_bbox_dim = &phys_bbox_dim_view.vector;

  // Set default padding to metric ellipse bounding box in this dimension (is zero for non-tiled dimensions)
  double phys_lower_pad = 0.5 * gsl_vector_get(tiling->phys_bbox, dimension);
  double phys_upper_pad = phys_lower_pad;

  // Compute lower parameter space bound
  (bound->func)(dimension, phys_point_dim, phys_bbox_dim, bound->data_lower, phys_lower, &phys_lower_pad);

  // If this dimension is non-tiled, we're done
  if (!bound->tiled) {
    *phys_upper = *phys_lower;
    return;
  }

  // Compute upper parameter space bound
  (bound->func)(dimension, phys_point_dim, phys_bbox_dim, bound->data_upper, phys_upper, &phys_upper_pad);

  // Optionally add padding
  if (padded) {
    *phys_lower -= phys_lower_pad;
    *phys_upper += phys_upper_pad;
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
    for (int32_t i = lookup->int_lower; i <= lookup->int_upper; ++i) {
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
  int int_lower[]				///< [in] Current integer lower bound
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

  // Set (i)th integer lower bound to lookup lower bound, then
  // transform to physical coordinates from generating integers
  int_lower[ti] = lookup->int_lower;
  const size_t i = gsl_vector_uint_get(tiling->tiled, ti);
  double phys_lower_i = gsl_vector_get(tiling->phys_offset, i);
  for (size_t tj = 0; tj < tn; ++tj) {
    const double phys_from_int_i_tj = gsl_matrix_get(tiling->phys_from_int, i, tj);
    phys_lower_i += phys_from_int_i_tj * int_lower[tj];
  }

  // Calculate physical upper bound from physical lower_bound
  const double phys_from_int_i_ti = gsl_matrix_get(tiling->phys_from_int, i, ti);
  double phys_upper_i = phys_lower_i + phys_from_int_i_ti * (lookup->int_upper - lookup->int_lower);

  // Print information on the current lookup trie dimension
  fprintf(fp, "dim #%zu/%zu   min %+5i   max %+5i   lower %+10e   upper %+10e",
          ti + 1, tn, lookup->int_lower, lookup->int_upper, phys_lower_i, phys_upper_i);

  // If this is the highest dimension, print index and return
  if (ti + 1 == tn) {
    fprintf(fp, "   index %" LAL_UINT8_FORMAT "\n", lookup->index);
  } else {

    // Otherwise, loop over this dimension
    fprintf(fp, "\n");
    LT_IndexLookup* next = lookup->next;
    for (int32_t point = lookup->int_lower; point <= lookup->int_upper; ++point, ++next) {

      // Set (i)th integer lower bound to this point
      int_lower[ti] = point;

      // Print higher dimensions
      LT_PrintLookup(tiling, fp, ti + 1, tn, next, int_lower);

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
  int* nearest_int_point,			///< [in/out] Nearest point in generating integer space
  UINT8* nearest_index				///< [in/out] Index of nearest point
  )
{

  // Find the nearest point to 'int_point', which is in generating integer space
  switch (lattice) {

  case LATTICE_TYPE_CUBIC:   // Cubic (\f$Z_n\f$) lattice
  {

    // Round each 'int_point' to nearest integer to find the nearest point in Zn
    for (size_t i = 0; i < n; ++i) {
      nearest_int_point[i] = lround(gsl_vector_get(int_point, i));
    }

  }
  break;

  case LATTICE_TYPE_ANSTAR:   // An-star (\f$A_n^*\f$) lattice
  {

    // The nearest point algorithm used below embeds the An* lattice in n+1 dimensions,
    // however 'int_point' is only n-dimensional. The algorithm is, however, sensitive
    // to the n differences between n dimensions and the (n+1)th dimension, so we can
    // freely set one of the dimensions to a constant value. We choose to set the 0th
    // dimension to zero, i.e. the (n+1)-dimensional lattice point is (0,'int_point').
    double y[n+1];
    y[0] = 0;
    for (size_t i = 0; i < n; ++i) {
      y[i+1] = gsl_vector_get(int_point, i);
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
    int k[n+1];
    {

      // Lines 1--4, 20
      double z[n+1], alpha = 0, beta = 0;
      size_t bucket[n+1], link[n+1];
      for (size_t i = 1; i <= n + 1; ++i) {
        k[i-1] = lround(y[i-1]);   // Line 20, moved here to avoid duplicate round
        z[i-1] = y[i-1] - k[i-1];
        alpha += z[i-1];
        beta += z[i-1]*z[i-1];
        bucket[i-1] = 0;
      }

      // Lines 5--8
      // Notes:
      //   * Correction to line 6, as as written in McKilliam et.al.:
      //       i = n + 1 - (n + 1)*floor(z_t + 0.5)
      //     should instead read
      //       i = n + 1 - floor((n + 1)*(z_t + 0.5))
      //   * We also convert the floor() operation into a round():
      //       floor(x) = round(x - 0.5)
      //     so that lattice points themselves will fall in the middle
      //     of the rounded interval, and will therefore be mapped back
      //     to themselves robustly (i.e. against numerical errors),
      //     while points furthest from the lattice bound are at the
      //     edges of the rounded interval, where it is acceptable for
      //     them to be rounded to different nearest lattice points.
      //     Line 6 then simplifies to:
      //       i = n + 1 - round((n + 1)*(z_t + 0.5) - 0.5)
      //         = round((n + 1)*(0.5 - z_t) + 0.5)
      for (size_t t = 1; t <= n + 1; ++t) {
        const size_t i = lround((n + 1)*(0.5 - z[t-1]) + 0.5);
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

    // The nearest point in An* is the n differences between k[1]...k[n] and k[0]
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

    // Make sure (i)th dimension of nearest integer point is within lookup bounds
    nearest_int_point[i] = GSL_MAX( lookup->int_lower, GSL_MIN( nearest_int_point[i], lookup->int_upper ) );

    // Jump to the next level of the lookup trie, based in index
    lookup = &lookup->next[nearest_int_point[i] - lookup->int_lower];

  }

  // Make sure (n-1)th dimension of nearest integer point is within lookup bounds
  nearest_int_point[n-1] = GSL_MAX( lookup->int_lower, GSL_MIN( nearest_int_point[n-1], lookup->int_upper ) );

  // Set index of nearest point
  *nearest_index = lookup->index + (nearest_int_point[n-1] - lookup->int_lower);

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
    const double bounding_box_i = 2.0 * sqrt(max_mismatch * gsl_matrix_get(inverse, i, i) / norm_i);
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
  tiling->dims = n;
  tiling->status = LT_S_INCOMPLETE;
  tiling->lattice = LATTICE_TYPE_MAX;

  // Allocate bound info array
  tiling->bounds = XLALCalloc(n, sizeof(*tiling->bounds));
  XLAL_CHECK_NULL(tiling->bounds != NULL, XLAL_ENOMEM);

  return tiling;

}

void XLALDestroyLatticeTiling(
  LatticeTiling* tiling
  )
{

  if (tiling) {

    const size_t n = tiling->dims;
    const size_t tn = (tiling->tiled != NULL) ? tiling->tiled->size : 0;

    // Free bound info array
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
    GFMAT(tiling->int_from_phys, tiling->phys_from_int);
    GFVEC(tiling->phys_bbox, tiling->phys_offset, tiling->phys_point);
    GFVECI(tiling->int_lower, tiling->int_point, tiling->int_upper);
    GFVECU(tiling->tiled);

    // Free tiling structure
    XLALFree(tiling);

  }

}

size_t XLALLatticeDimensions(
  const LatticeTiling* tiling
  )
{

  // Check input
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);

  return tiling->dims;

}

gsl_vector_uint* XLALLatticeTiledDimensions(
  const LatticeTiling* tiling
  )
{

  // Check input
  XLAL_CHECK_NULL(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(tiling->status > LT_S_INCOMPLETE, XLAL_EINVAL);

  // Return tiled dimensions of the parameter space
  return tiling->tiled;

}

UINT8 XLALLatticePointCount(
  const LatticeTiling* tiling
  )
{

  // Check input
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_VAL(0, tiling->status >= LT_S_INITIALISED, XLAL_EINVAL);

  return tiling->count;

}

UINT8 XLALCountLatticePoints(
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

gsl_matrix* XLALLatticeBasisVectors(
  const LatticeTiling* tiling
  )
{

  // Check input
  XLAL_CHECK_NULL(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(tiling->status > LT_S_INCOMPLETE, XLAL_EINVAL);

  // Return 'phys_from_int', whose columns are the basis vectors of the lattice
  return tiling->phys_from_int;

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
  XLAL_CHECK(dimension < tiling->dims, XLAL_ESIZE);
  XLAL_CHECK(func != NULL, XLAL_EFAULT);
  XLAL_CHECK(data_len > 0, XLAL_EFAULT);
  XLAL_CHECK(data_lower != NULL, XLAL_EFAULT);
  XLAL_CHECK(data_upper != NULL, XLAL_EFAULT);

  // Check that bound has not already been set
  XLAL_CHECK(tiling->bounds[dimension].func == NULL, XLAL_EINVAL,
             "Dimension #%i has already been bounded", dimension);

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
  XLAL_CHECK(XLALSetLatticeBound(tiling, dimension, ConstantBound, data_len, data_lower, data_upper) == XLAL_SUCCESS, XLAL_EFUNC);

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
  const size_t n = tiling->dims;
  XLAL_CHECK(metric != NULL, XLAL_EFAULT);
  XLAL_CHECK(metric->size1 == n && metric->size2 == n, XLAL_EINVAL);
  XLAL_CHECK(max_mismatch > 0, XLAL_EINVAL);

  // Save the type of lattice to generate tiling with
  tiling->lattice = lattice;

  // Check metric is symmetric and has positive diagonal elements
  for (size_t i = 0; i < n; ++i) {
    XLAL_CHECK(gsl_matrix_get(metric, i, i) > 0, XLAL_EINVAL, "metric(%zu,%zu) <= 0", i, i);
    for (size_t j = i + 1; j < n; ++j) {
      XLAL_CHECK(gsl_matrix_get(metric, i, j) == gsl_matrix_get(metric, j, i), XLAL_EINVAL, "metric(%zu,%zu) != metric(%zu,%zu)", i, j, j, i);
    }
  }

  // Check that all parameter-space dimensions are bounded, and record indices of tiled dimensions
  {
    size_t tn = 0;
    for (size_t i = 0; i < tiling->dims; ++i) {
      XLAL_CHECK(tiling->bounds[i].func != NULL, XLAL_EFAILED, "Dimension #%i is unbounded", i);
      if (tiling->bounds[i].tiled) {
        ++tn;
      }
    }
    if (tn > 0) {
      GAVECU(tiling->tiled, tn);
      for (size_t i = 0, ti = 0; i < tiling->dims; ++i) {
        if (tiling->bounds[i].tiled) {
          gsl_vector_uint_set(tiling->tiled, ti++, i);
        }
      }
    }
  }
  const size_t tn = (tiling->tiled != NULL) ? tiling->tiled->size : 0;

  // Allocate memory
  GAVEC(tiling->phys_bbox, n);
  GAVEC(tiling->phys_offset, n);
  GAVEC(tiling->phys_point, n);

  // Initialise bounding box to zero
  gsl_vector_set_zero(tiling->phys_bbox);

  // Set physical parameter-space offset from parameter-space bounds
  for (size_t i = 0; i < n; ++i) {
    double phys_lower_i = 0, phys_upper_i = 0;
    LT_GetPhysBounds(tiling, i, tiling->phys_offset, false, &phys_lower_i, &phys_upper_i);
    gsl_vector_set(tiling->phys_offset, i, 0.5*(phys_lower_i + phys_upper_i));
  }

  // If there are tiled dimensions:
  if (tn > 0) {

    // Allocate memory
    GAVECI(tiling->int_lower, tn);
    GAVECI(tiling->int_point, tn);
    GAVECI(tiling->int_upper, tn);

    // Calculate normalisation scale from metric diagonal elements
    gsl_vector* GAVEC(t_norm, tn);
    for (size_t ti = 0; ti < tn; ++ti) {
      const size_t i = gsl_vector_uint_get(tiling->tiled, ti);
      const double metric_i_i = gsl_matrix_get(metric, i, i);
      gsl_vector_set(t_norm, ti, sqrt(metric_i_i));
    }

    // Copy and normalise tiled dimensions of metric
    gsl_matrix* GAMAT(t_metric, tn, tn);
    for (size_t ti = 0; ti < tn; ++ti) {
      const size_t i = gsl_vector_uint_get(tiling->tiled, ti);
      const double t_norm_ti = gsl_vector_get(t_norm, ti);
      for (size_t tj = 0; tj < tn; ++tj) {
        const size_t j = gsl_vector_uint_get(tiling->tiled, tj);
        const double t_norm_tj = gsl_vector_get(t_norm, tj);
        const double metric_i_j = gsl_matrix_get(metric, i, j);
        gsl_matrix_set(t_metric, ti, tj, metric_i_j / t_norm_ti / t_norm_tj);
      }
    }

    // Check tiled metric is positive definite, by trying to compute its Cholesky decomposition
    gsl_matrix* GAMAT(t_metric_copy, tn, tn);
    gsl_matrix_memcpy(t_metric_copy, t_metric);
    GCALL(gsl_linalg_cholesky_decomp(t_metric_copy), "tiled metric is not positive definite");

    // Compute metric ellipse bounding box
    gsl_vector* t_bbox = XLALMetricEllipseBoundingBox(t_metric, max_mismatch);
    XLAL_CHECK(t_bbox != NULL, XLAL_EFUNC);

    // Copy bounding box in physical coordinates to tiled dimensions
    for (size_t ti = 0; ti < tn; ++ti) {
      const size_t i = gsl_vector_uint_get(tiling->tiled, ti);
      const double t_norm_ti = gsl_vector_get(t_norm, ti);
      gsl_vector_set(tiling->phys_bbox, i, gsl_vector_get(t_bbox, ti) / t_norm_ti);
    }

    // Compute a lower triangular basis matrix whose columns are orthonormal with respect to the tiled metric
    gsl_matrix* t_basis = XLALComputeMetricOrthoBasis(t_metric);
    XLAL_CHECK(t_basis != NULL, XLAL_EFUNC);

    // Compute a lower triangular generator matrix for a given lattice type and mismatch
    gsl_matrix* t_generator = XLALComputeLatticeGenerator(tn, lattice, max_mismatch);
    XLAL_CHECK(t_generator != NULL, XLAL_EFUNC);

    // Compute transform to normalised coordinates from generating integers,
    // which is the lattice generator expressed in the orthonormal metric basis
    gsl_matrix* GAMAT(t_norm_from_int, tn, tn);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, t_basis, t_generator, 0.0, t_norm_from_int);

    // Transform will be lower triangular, so zero out upper triangle
    LT_ZeroStrictUpperTriangle(t_norm_from_int);

    // Compute inverse transform to generating integers from normalised coordinates
    gsl_matrix* GAMAT(t_int_from_norm, tn, tn);
    gsl_matrix_set_identity(t_int_from_norm);
    gsl_blas_dtrsm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, t_norm_from_int, t_int_from_norm);

    // Inverse of transform will be lower triangular, so zero out upper triangle
    LT_ZeroStrictUpperTriangle(t_int_from_norm);

    // Copy transform to tiling structure, inserting zero rows/columns for non-tiled dimensions
    GAMAT(tiling->phys_from_int, n, tn);
    GAMAT(tiling->int_from_phys, tn, n);
    gsl_matrix_set_zero(tiling->phys_from_int);
    gsl_matrix_set_zero(tiling->int_from_phys);
    for (size_t ti = 0; ti < tn; ++ti) {
      const size_t i = gsl_vector_uint_get(tiling->tiled, ti);
      const double t_norm_ti = gsl_vector_get(t_norm, ti);
      for (size_t tj = 0; tj < tn; ++tj) {
        const size_t j = gsl_vector_uint_get(tiling->tiled, tj);
        const double t_norm_tj = gsl_vector_get(t_norm, tj);
        const double t_norm_from_int_ti_tj = gsl_matrix_get(t_norm_from_int, ti, tj);
        gsl_matrix_set(tiling->phys_from_int, i, tj, t_norm_from_int_ti_tj / t_norm_ti);
        const double t_int_from_norm_ti_tj = gsl_matrix_get(t_int_from_norm, ti, tj);
        gsl_matrix_set(tiling->int_from_phys, ti, j, t_int_from_norm_ti_tj * t_norm_tj);
      }
    }

    // Shift the physical offset by half of a step-size in each tiled dimension.
    // This ensures that the tiling will never place a lattice point at zero in
    // physical coordinates, but rather two points either side of zero. This is
    // because the physical coordinates may not be well-defined at zero.
    for (size_t ti = 0; ti < tn; ++ti) {
      const size_t i = gsl_vector_uint_get(tiling->tiled, ti);
      const double phys_from_int_i_ti = gsl_matrix_get(tiling->phys_from_int, i, ti);
      const double phys_offset_i = gsl_vector_get(tiling->phys_offset, i);
      gsl_vector_set(tiling->phys_offset, i, phys_offset_i + 0.5 * phys_from_int_i_ti);
    }

    // Cleanup
    GFMAT(t_basis, t_generator, t_int_from_norm, t_metric, t_metric_copy, t_norm_from_int);
    GFVEC(t_bbox, t_norm);

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
  const size_t n = tiling->dims;
  const size_t tn = (tiling->tiled != NULL) ? tiling->tiled->size : 0;
  XLAL_CHECK(curr_point == NULL || curr_point->size == n, XLAL_EINVAL);

  // If finished status, nothing more to be done!
  if (tiling->status == LT_S_FINISHED) {
    return -1;
  }

  // Which dimensions have changed?
  size_t changed_ti;

  // Which dimensions need to be reset to lower bounds?
  size_t reset_ti;

  if (tiling->status == LT_S_INITIALISED) {

    // Initialise lattice point
    gsl_vector_int_set_zero(tiling->int_point);
    gsl_vector_set_zero(tiling->phys_point);

    // Initialise count
    tiling->count = 1;

    // ALl dimensions have changed
    changed_ti = 0;

    // All dimensions need to be reset to lower bounds
    reset_ti = 0;

    // Tiling has been started
    tiling->status = LT_S_STARTED;

  } else {

    // Find the next lattice point
    size_t ti = tn;
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

      // Increment integer point in this dimension
      const int int_point_ti = gsl_vector_int_get(tiling->int_point, ti) + 1;
      gsl_vector_int_set(tiling->int_point, ti, int_point_ti);

      // Increment physical point in this dimension
      gsl_vector_const_view phys_from_int_ti = gsl_matrix_const_column(tiling->phys_from_int, ti);
      gsl_vector_add(tiling->phys_point, &phys_from_int_ti.vector);

      // If point is not out of bounds, we have found the next lattice point
      const int int_upper_ti = gsl_vector_int_get(tiling->int_upper, ti);
      if (int_point_ti <= int_upper_ti) {
        break;
      }

      // Move on to lower dimensions
      continue;

    }

    // Point was found, so increase count
    ++tiling->count;

    // This dimension and higher have changed
    changed_ti = ti;

    // Higher dimensions need to be reset to lower bounds
    reset_ti = ti + 1;

  }

  // Reset specified dimensions to lower bounds
  for (size_t i = 0, ti = 0; i < n; ++i) {

    // Get physical bounds with padding
    double phys_lower_i = 0, phys_upper_i = 0;
    LT_GetPhysBounds(tiling, i, tiling->phys_point, true, &phys_lower_i, &phys_upper_i);

    // If not tiled, set current physical point to non-tiled parameter-space bound
    if (!tiling->bounds[i].tiled) {
      gsl_vector_set(tiling->phys_point, i, phys_lower_i);
      continue;
    }

    // If tiled dimension needs to be reset to lower bound:
    if (ti >= reset_ti) {

      // Transform physical point in lower dimensions to generating integer offset
      const double phys_offset_i = gsl_vector_get(tiling->phys_offset, i);
      double int_from_phys_point_i = 0;
      for (size_t j = 0; j < i; ++j) {
        const double int_from_phys_ti_j = gsl_matrix_get(tiling->int_from_phys, ti, j);
        const double phys_point_j = gsl_vector_get(tiling->phys_point, j);
        const double phys_offset_j = gsl_vector_get(tiling->phys_offset, j);
        int_from_phys_point_i += int_from_phys_ti_j * (phys_point_j - phys_offset_j);
      }

      // Transform physical bounds to generating integers
      const double int_from_phys_ti_i = gsl_matrix_get(tiling->int_from_phys, ti, i);
      const double int_lower_i = int_from_phys_point_i + int_from_phys_ti_i * (phys_lower_i - phys_offset_i);
      const double int_upper_i = int_from_phys_point_i + int_from_phys_ti_i * (phys_upper_i - phys_offset_i);

      // Set integer lower/upper bounds, rounded up/down to avoid extra boundary points
      gsl_vector_int_set(tiling->int_lower, ti, lround(ceil(int_lower_i)));
      gsl_vector_int_set(tiling->int_upper, ti, lround(floor(int_upper_i)));

      // Set integer point to lower bound
      gsl_vector_int_set(tiling->int_point, ti, gsl_vector_int_get(tiling->int_lower, ti));

      // Set current physical point from integer point
      double phys_point_i = phys_offset_i;
      for (size_t tj = 0; tj < tn; ++tj) {
        const double phys_from_int_i_tj = gsl_matrix_get(tiling->phys_from_int, i, tj);
        const int int_point_tj = gsl_vector_int_get(tiling->int_point, tj);
        phys_point_i += phys_from_int_i_tj * int_point_tj;
      }
      gsl_vector_set(tiling->phys_point, i, phys_point_i);

    }

    ++ti;

  }

  // Optionally, copy current physical point
  if (curr_point != NULL) {
    gsl_vector_memcpy(curr_point, tiling->phys_point);
  }

  return changed_ti;

}

int XLALFastForwardLatticeTiling(
  LatticeTiling* tiling,
  size_t* point_count,
  double* point_spacing
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status == LT_S_STARTED, XLAL_EINVAL);
  const size_t tn = (tiling->tiled != NULL) ? tiling->tiled->size : 0;

  // If no tiled dimensions, nothing to fast-forward
  size_t count = 0;
  double spacing = 0.0;

  // If there are tiled dimensions:
  if (tn > 0) {

    // Get integer point and upper bound in highest tiled dimension
    const int int_point = gsl_vector_int_get(tiling->int_point, tn - 1);
    const int int_upper = gsl_vector_int_get(tiling->int_upper, tn - 1);

    // Set integer point in highest tiled dimension to upper bound, so that the next
    // call to XLALNextLatticePoint() will advance the next-highest tiled dimension
    gsl_vector_int_set(tiling->int_point, tn - 1, int_upper);
    count = int_upper - int_point;
    tiling->count += count;

    // Return physical spacing of lattice points in highest tiled dimension
    const size_t i = gsl_vector_uint_get(tiling->tiled, tn - 1);
    spacing = gsl_matrix_get(tiling->phys_from_int, i, tn - 1);

  }

  // Return count and spacing
  if (point_count != NULL) {
    *point_count = count;
  }
  if (point_spacing != NULL) {
    *point_spacing = spacing;
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

gsl_matrix* XLALLatticeTilingSubset(
  LatticeTiling* tiling,
  size_t subset_dimension
  )
{

  // Check input
  XLAL_CHECK_NULL(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(tiling->status > LT_S_INCOMPLETE, XLAL_EFAILED);
  const size_t n = tiling->dims;
  const size_t tn = (tiling->tiled != NULL) ? tiling->tiled->size : 0;
  XLAL_CHECK_NULL(0 < subset_dimension && subset_dimension <= n, XLAL_EINVAL);
  const size_t subset_n = subset_dimension;

  // Count how many tiled dimensions are in the subset dimensions
  size_t subset_tn = 0;
  for (size_t i = 0; i < subset_n; ++i) {
    if (tiling->bounds[i].tiled) {
      ++subset_tn;
    }
  }

  // Save the total tiling count, since it will be destroyed by what follows
  const UINT8 save_total_count = tiling->total_count;

  // Restart the tiling, and iterate over all points
  XLAL_CHECK_NULL(XLALRestartLatticeTiling(tiling) == XLAL_SUCCESS, XLAL_EFUNC);
  while (XLALNextLatticePoint(tiling, NULL) >= 0) {

    // In all dimensions higher than the subset dimension, skip to the very last point
    for (size_t ti = subset_tn; ti < tn; ++ti) {
      const int int_upper_ti = gsl_vector_int_get(tiling->int_upper, ti);
      gsl_vector_int_set(tiling->int_point, ti, int_upper_ti);
    }

  }

  // Allocate a matrix to hold the subset points; the tiling count
  // will now hold the total number of points in the subset
  gsl_matrix* GAMAT_NULL(subset_points, subset_n, tiling->total_count);

  // Restart the tiling, and iterate over all points
  XLAL_CHECK_NULL(XLALRestartLatticeTiling(tiling) == XLAL_SUCCESS, XLAL_EFUNC);
  while (XLALNextLatticePoint(tiling, NULL) >= 0) {

    // Copy the first 'subset_n' dimensions of each point to the subset matrix
    for (size_t i = 0; i < subset_n; ++i) {
      const double phys_point_i = gsl_vector_get(tiling->phys_point, i);
      gsl_matrix_set(subset_points, i, tiling->count - 1, phys_point_i);
    }

    // In all dimensions higher than the subset dimension, skip to the very last point
    for (size_t ti = subset_tn; ti < tn; ++ti) {
      const int int_upper_ti = gsl_vector_int_get(tiling->int_upper, ti);
      gsl_vector_int_set(tiling->int_point, ti, int_upper_ti);
    }

  }

  // Restart the tiling, and restore the total tiling count
  XLAL_CHECK_NULL(XLALRestartLatticeTiling(tiling) == XLAL_SUCCESS, XLAL_EFUNC);
  tiling->total_count = save_total_count;

  return subset_points;

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
  const size_t n = tiling->dims;

  // Check input
  XLAL_CHECK(rng != NULL, XLAL_EFAULT);
  XLAL_CHECK(random_points != NULL, XLAL_EFAULT);
  XLAL_CHECK(random_points->size1 == n, XLAL_ESIZE);

  // Generate random points in lattice tiling parameter space
  for (size_t k = 0; k < random_points->size2; ++k) {
    gsl_vector_view phys_point = gsl_matrix_column(random_points, k);
    for (size_t i = 0; i < n; ++i) {

      // Get physical bounds
      double phys_lower_i = 0, phys_upper_i = 0;
      LT_GetPhysBounds(tiling, i, &phys_point.vector, false, &phys_lower_i, &phys_upper_i);

      // Generate random number
      const double u = XLALUniformDeviate(rng);

      // Set parameter space point
      gsl_vector_set(&phys_point.vector, i, phys_lower_i + u*(phys_upper_i - phys_lower_i));

    }

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
  const size_t tn = (tiling->tiled != NULL) ? tiling->tiled->size : 0;

  // If lookup has already been built, or there are no tiled dimensions, we're done!
  if (tiling->lookup_base != NULL || tn == 0) {
    return XLAL_SUCCESS;
  }

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

    // Iterate over all dimensions where the current point has changed
    for (size_t tj = ti; tj < tn; ++tj) {

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

        // Save the lower and upper integer point bounds
        lookup->int_lower = gsl_vector_int_get(tiling->int_lower, tj);
        lookup->int_upper = gsl_vector_int_get(tiling->int_upper, tj);

        // If we are below the highest dimension:
        if (tj + 1 < tn) {

          // Allocate a new array of lookup structs for the next highest dimension
          const size_t next_length = lookup->int_upper - lookup->int_lower + 1;
          lookup->next = XLALCalloc(next_length, sizeof(*lookup->next));
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
  const size_t tn = (tiling->tiled != NULL) ? tiling->tiled->size : 0;
  XLAL_CHECK(file != NULL, XLAL_EFAULT);

  // Print lookup trie
  int int_point[tn];
  LT_PrintLookup(tiling, file, 0, tn, tiling->lookup_base, int_point);

  return XLAL_SUCCESS;

}

int XLALNearestLatticePoints(
  const LatticeTiling* tiling,
  const gsl_matrix* points,
  gsl_matrix* nearest_points,
  UINT8Vector* nearest_indices,
  gsl_matrix** workspace
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status > LT_S_INCOMPLETE, XLAL_EFAILED);
  const size_t n = tiling->dims;
  const size_t tn = (tiling->tiled != NULL) ? tiling->tiled->size : 0;
  XLAL_CHECK(points != NULL, XLAL_EFAULT);
  XLAL_CHECK(points->size1 == n, XLAL_ESIZE);
  const size_t num_points = points->size2;
  XLAL_CHECK(nearest_points == NULL || nearest_points->size1 == n, XLAL_ESIZE);
  XLAL_CHECK(nearest_points == NULL || nearest_points->size2 >= num_points, XLAL_ESIZE);
  XLAL_CHECK(nearest_indices == NULL || nearest_indices->length >= num_points, XLAL_ESIZE);
  XLAL_CHECK(nearest_indices == NULL || nearest_indices->data != NULL, XLAL_EFAULT);
  XLAL_CHECK(nearest_indices == NULL || tn == 0 || tiling->lookup_base != NULL, XLAL_EFAILED,
             "Need to first run XLALBuildLatticeIndexLookup() to build index lookup" );
  XLAL_CHECK(workspace != NULL, XLAL_EFAULT);

  // If there are no tiled dimensions:
  if (tn == 0) {

    // Set all columns of 'nearest_points' to the sole point in the tiling,
    // which is just the physical offset given by 'phys_offset'
    if (nearest_points != NULL) {
      for (size_t i = 0; i < n; ++i) {
        const double phys_offset = gsl_vector_get(tiling->phys_offset, i);
        gsl_vector_view nearest_points_row = gsl_matrix_row(nearest_points, i);
        gsl_vector_set_all(&nearest_points_row.vector, phys_offset);
      }
    }

    // Set all elements of 'nearest_indices' to zero, since there is only one point
    if (nearest_indices != NULL) {
      memset(nearest_indices->data, 0, nearest_indices->length * sizeof(nearest_indices->data[0]));
    }

    return XLAL_SUCCESS;

  }

  // (Re)Allocate workspace matrix
  const size_t wksp_rows = tn + ((nearest_points == NULL) ? n : 0);
  if (*workspace != NULL && ((*workspace)->size1 != wksp_rows || (*workspace)->size2 < num_points)) {
    gsl_matrix_free(*workspace);
    *workspace = NULL;
  }
  if (*workspace == NULL) {
    GAMAT(*workspace, wksp_rows, num_points);
  }

  // Create view of either 'nearest_points' or '*workspace', for storing nearest points
  gsl_matrix_view nearest_view = gsl_matrix_submatrix((nearest_points == NULL) ? *workspace : nearest_points,
                                                      (nearest_points == NULL) ? tn : 0,
                                                      0, n, num_points);
  gsl_matrix *const nearest = &nearest_view.matrix;

  // Create view of '*workspace' of the required size, for storing tiled dimensions of nearest points
  gsl_matrix_view t_nearest_view = gsl_matrix_submatrix(*workspace, 0, 0, tn, num_points);
  gsl_matrix *const t_nearest = &t_nearest_view.matrix;

  // Copy 'points' to 'nearest', then subtract physical offset from every point
  gsl_matrix_memcpy(nearest, points);
  for (size_t i = 0; i < n; ++i) {
    const double phys_offset = gsl_vector_get(tiling->phys_offset, i);
    gsl_vector_view nearest_row = gsl_matrix_row(nearest, i);
    gsl_vector_add_constant(&nearest_row.vector, -phys_offset);
  }

  // Transform 'nearest' points to generating integers from physical coordinates, storing result in 't_nearest'
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tiling->int_from_phys, nearest, 0.0, t_nearest);

  // Find the nearest points in the lattice tiling to the 't_nearest' points
  int nearest_int_point[tn];
  for (size_t j = 0; j < num_points; ++j) {
    gsl_vector_view t_nearest_col = gsl_matrix_column(t_nearest, j);
    UINT8 nearest_index = 0;
    LT_FindNearestPoint(tn, tiling->lattice, tiling->lookup_base, &t_nearest_col.vector, nearest_int_point, &nearest_index);
    if (nearest_indices != NULL) {
      nearest_indices->data[j] = nearest_index;
    }
    for (size_t ti = 0; ti < tn; ++ti) {
      gsl_vector_set(&t_nearest_col.vector, ti, nearest_int_point[ti]);
    }
  }
  XLAL_CHECK( xlalErrno == 0, XLAL_EFAILED, "LT_FindNearestPoint() failed" );

  if (nearest_points != NULL) {

    // Transform 't_nearest' points to physical coordinates from generating integers, storing result in 'nearest'
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tiling->phys_from_int, t_nearest, 0.0, nearest);

    // Add physical offset back to every point in 'nearest'
    for (size_t i = 0; i < n; ++i) {
      const double phys_offset = gsl_vector_get(tiling->phys_offset, i);
      gsl_vector_view nearest_row = gsl_matrix_row(nearest, i);
      gsl_vector_add_constant(&nearest_row.vector, phys_offset);
    }

    // Set any non-tiled dimensions in 'nearest_points'
    for (size_t i = 0; i < n; ++i) {
      if (!tiling->bounds[i].tiled) {
        for (size_t j = 0; j < num_points; ++j) {
          gsl_vector_view nearest_points_col = gsl_matrix_column(nearest_points, j);

          // Get physical bounds
          double phys_lower_i = 0, phys_upper_i = 0;
          LT_GetPhysBounds(tiling, i, &nearest_points_col.vector, false, &phys_lower_i, &phys_upper_i);

          // Set point to non-tiled parameter-space bound
          gsl_vector_set(&nearest_points_col.vector, i, phys_lower_i);

        }
      }
    }

  }

  return XLAL_SUCCESS;

}
