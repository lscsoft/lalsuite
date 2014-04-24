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
#include <sys/types.h>
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

#include <lal/FlatLatticeTiling.h>
#include <lal/LALStdlib.h>
#include <lal/LALMalloc.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/GSLSupport.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

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
/// Flat lattice tiling bound info
///
typedef struct tagFLT_Bound {
  bool tiled;					///< Is the bound tiled, i.e. not a single point?
  size_t N;					///< Number of bound coefficients (\f$N\f$)
  gsl_vector* a;				///< Vector of offsets (\f$a\f$)
  gsl_matrix* c_lower;				///< Matrix of coefficients (\f$c\f$) for the lower bound
  gsl_matrix* m_lower;				///< Matrix of exponents (\f$m\f$) for the lower bound
  gsl_matrix* c_upper;				///< Matrix of coefficients (\f$c\f$) for the upper bound
  gsl_matrix* m_upper;				///< Matrix of exponents (\f$m\f$) for the upper bound
} FLT_Bound;

///
/// Flat lattice tiling state structure
///
struct tagFlatLatticeTiling {
  size_t dimensions;				///< Dimension of the parameter space
  FLT_Status status;				///< Status of the tiling
  FLT_Bound *bounds;				///< Array of parameter-space bound info for each dimension
  gsl_vector_uint* tiled_idx;			///< Indices of the tiled dimensions of the parameter space
  FlatLatticeType lattice;			///< Type of lattice to generate flat tiling with
  gsl_vector* phys_scale;			///< Normalised to physical coordinate scale
  gsl_vector* phys_offset;			///< Normalised to physical coordinate offset
  gsl_matrix* metric;				///< Normalised parameter-space metric
  gsl_matrix* increment;			///< Increment vectors of the lattice tiling generator
  gsl_vector* padding;				///< Padding at edges of parameter-space bounds
  gsl_vector* point;				///< Current lattice point
  gsl_vector* lower;				///< Current lower bound on parameter space
  gsl_vector* upper;				///< Current upper bound on parameter space
  uint64_t count;				///< Number of points generated so far
  uint64_t total_count;				///< Total number of points in parameter space
};

///
/// Find the bounding box of the mismatch ellipses of a metric
///
static gsl_vector* FLT_MetricEllipseBoundingBox(
  const gsl_matrix* metric,			///< [in] Metric to bound
  const double max_mismatch			///< [in] Maximum mismatch with respect to metric
  )
{

  // Check input
  XLAL_CHECK_NULL(metric != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(metric->size1 == metric->size2, XLAL_ESIZE);
  const size_t n = metric->size1;

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

///
/// Orthonormalise the columns of a matrix with respect to a metric (matrix is lower triangular)
///
static int FLT_OrthonormaliseWRTMetric(
  gsl_matrix* matrix,				///< [in] Matrix of columns to orthonormalise
  const gsl_matrix* metric			///< [in] Metric to orthonormalise with respect to
  )
{

  // Check input
  XLAL_CHECK(matrix != NULL, XLAL_EFAULT);
  XLAL_CHECK(metric != NULL, XLAL_EFAULT);
  XLAL_CHECK(metric->size1 == metric->size2, XLAL_ESIZE);
  XLAL_CHECK(metric->size1 == matrix->size2 && metric->size2 == matrix->size2, XLAL_ESIZE);
  const size_t n = metric->size1;

  // Allocate temporary vector
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

///
/// Transform a lattice generator to a square lower triangular form
///
static gsl_matrix* FLT_SquareLowerTriangularLatticeGenerator(
  gsl_matrix* generator				///< [in] Generator matrix of lattice
  )
{

  // Check input
  XLAL_CHECK_NULL(generator != NULL, XLAL_EFAULT);
  const size_t m = generator->size1;
  const size_t n = generator->size2;
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

///
/// Normalise a lattice generator matrix to have a specified covering radius
///
static int FLT_NormaliseLatticeGenerator(
  gsl_matrix* generator,			///< [in] Generator matrix of lattice
  const double norm_thickness,			///< [in] Normalised thickness of lattice
  const double covering_radius			///< [in] Desired covering radius
  )
{

  // Check input
  XLAL_CHECK(generator != NULL, XLAL_EFAULT);
  XLAL_CHECK(generator->size1 == generator->size2, XLAL_ESIZE);
  const size_t n = generator->size1;

  // Allocate memory
  gsl_matrix* LU_decomp = gsl_matrix_alloc(n, n);
  XLAL_CHECK(LU_decomp != NULL, XLAL_ENOMEM);
  gsl_permutation* LU_perm = gsl_permutation_alloc(n);
  XLAL_CHECK(LU_perm != NULL, XLAL_ENOMEM);

  // Compute generator LU decomposition
  gsl_matrix_memcpy(LU_decomp, generator);
  int LU_sign = 0;
  gsl_linalg_LU_decomp(LU_decomp, LU_perm, &LU_sign);

  // Compute generator determinant
  const double generator_determinant = gsl_linalg_LU_det(LU_decomp, LU_sign);

  // Compute covering radius
  const double generator_covering_radius = pow(norm_thickness * generator_determinant, 1.0 / n);

  // Normalise so covering spheres have specified covering radius
  gsl_matrix_scale(generator, covering_radius / generator_covering_radius);

  // Cleanup
  gsl_matrix_free(LU_decomp);
  gsl_permutation_free(LU_perm);

  return XLAL_SUCCESS;

}

///
/// Find the lattice increment vectors for a given lattice, metric and mismatch
///
static gsl_matrix* FLT_MetricLatticeIncrements(
  const FlatLatticeType lattice,		///< [in] Lattice type
  const gsl_matrix* metric,			///< [in] parameter-space metric
  const double max_mismatch			///< [in] Maximum prescribed mismatch
  )
{

  // Check input
  XLAL_CHECK_NULL(metric != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(metric->size1 == metric->size2, XLAL_ESIZE);
  XLAL_CHECK_NULL(max_mismatch > 0.0, XLAL_EINVAL);
  const size_t n = metric->size1;

  // Allocate memory
  gsl_matrix* directions = gsl_matrix_alloc(metric->size1, metric->size2);
  XLAL_CHECK_NULL(directions != NULL, XLAL_ENOMEM);
  gsl_matrix* increment = gsl_matrix_alloc(metric->size1, metric->size2);
  XLAL_CHECK_NULL(increment != NULL, XLAL_ENOMEM);

  // Check metric is positive definite, by trying to compute its Cholesky decomposition
  gsl_matrix_memcpy(directions, metric);   // Make copy to preserve original
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
  int retn = gsl_linalg_cholesky_decomp(directions);
  gsl_set_error_handler(old_handler);
  XLAL_CHECK_NULL(retn == 0, XLAL_EFAILED, "metric is not positive definite");

  // Find orthonormalise directions with respect to tiling metric
  gsl_matrix_set_identity(directions);
  XLAL_CHECK_NULL(FLT_OrthonormaliseWRTMetric(directions, metric) == XLAL_SUCCESS, XLAL_EFUNC);

  // Get lattice generator
  gsl_matrix* gen_matrix = NULL;
  double norm_thickness = 0.0;
  switch (lattice) {

  case FLAT_LATTICE_TYPE_CUBIC:	  // Cubic (\f$Z_n\f$) lattice
  {

    // Allocate memory
    gen_matrix = gsl_matrix_alloc(n, n);
    XLAL_CHECK_NULL(gen_matrix != NULL, XLAL_ENOMEM);

    // Create generator
    gsl_matrix_set_identity(gen_matrix);

    // Calculate normalised thickness
    norm_thickness = pow(sqrt(n)/2, n);

  }
  break;

  case FLAT_LATTICE_TYPE_ANSTAR:   // An-star (\f$A_n^*\f$) lattice
  {

    // Allocate memory
    gen_matrix = gsl_matrix_alloc(n + 1, n);
    XLAL_CHECK_NULL(gen_matrix != NULL, XLAL_ENOMEM);

    // Create generator in (n + 1) space
    gsl_matrix_set_all(gen_matrix, 0.0);
    gsl_vector_view first_row = gsl_matrix_row(gen_matrix, 0);
    gsl_vector_view sub_diag = gsl_matrix_subdiagonal(gen_matrix, 1);
    gsl_vector_view last_col = gsl_matrix_column(gen_matrix, n - 1);
    gsl_vector_set_all(&first_row.vector, 1.0);
    gsl_vector_set_all(&sub_diag.vector, -1.0);
    gsl_vector_set_all(&last_col.vector, 1.0 / (n + 1.0));
    gsl_vector_set(&last_col.vector, 0, -1.0 * n / (n + 1.0));

    // Calculate normalised thickness
    norm_thickness = sqrt(n + 1.0)*pow((1.0*n*(n + 2))/(12.0*(n + 1)), 0.5*n);

  }
  break;

  default:
    XLAL_ERROR_NULL(XLAL_EINVAL, "Invalid 'lattice'=%u", lattice);
  }

  // Transform lattice generator to square lower triangular
  gsl_matrix* sqlwtr_gen_matrix = FLT_SquareLowerTriangularLatticeGenerator(gen_matrix);
  XLAL_CHECK_NULL(sqlwtr_gen_matrix != NULL, XLAL_EFAILED);

  // Normalise lattice generator so covering radius is sqrt(mismatch)
  XLAL_CHECK_NULL(FLT_NormaliseLatticeGenerator(sqlwtr_gen_matrix, norm_thickness, sqrt(max_mismatch)) == XLAL_SUCCESS, XLAL_EFUNC);

  // Compute the increment vectors of the lattice generator along the orthogonal directions
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, directions, sqlwtr_gen_matrix, 0.0, increment);

  // Cleanup
  gsl_matrix_free(directions);
  gsl_matrix_free(gen_matrix);
  gsl_matrix_free(sqlwtr_gen_matrix);

  return increment;

}

///
/// Special power function used when calculating bounds
///
static inline double FLT_BoundPow(double x, double y) {
  double z = pow(x, y);
  return isfinite(z) ? z : 0.0;
}

///
/// Returns the parameter-space bound \f$X_n\f$ on dimension \f$n\f$; see XLALSetFlatLatticeBound().
///
static double FLT_GetBound(
  const size_t n,				///< [in] Dimension on which bound applies (\f$n\f$)
  const gsl_vector* x,				///< [in] Point at which to find bounds (\f$x\f$)
  const size_t N,				///< [in] Number of coefficients (\f$N\f$)
  const gsl_vector* a,				///< [in] Vector of offsets (\f$a\f$)
  const gsl_matrix* c,				///< [in] Matrix of coefficients (\f$c\f$)
  const gsl_matrix* m				///< [in] Matrix of exponents (\f$m\f$)
  )
{
  double bound = 0;
  for (size_t k = 0; k < /* P */ c->size2; ++k) {
    double bound_k = 0;
    for (size_t j = 0; j < N; ++j) {
      double bound_j = gsl_matrix_get(c, j, k);
      if (n > 0) {
        for (size_t i = 0; i < n; ++i) {
          const double dx_i = gsl_vector_get(x, i) - gsl_vector_get(a, i);
          bound_j *= FLT_BoundPow(dx_i, gsl_matrix_get(m, n*j + i, k));
        }
      }
      bound_k += bound_j;
    }
    bound_k = FLT_BoundPow(bound_k, gsl_matrix_get(m, n*N, k));
    bound_k = gsl_matrix_get(c, N, k) * bound_k + gsl_vector_get(a, n);
    bound += bound_k;
  }
  return bound;
}

///
/// Transform the parameter space bounds on coordinates \f$x\f$, represented by \f$(a,c,m)\f$, to
/// bounds on coordinates \f$\bar{x}\f$, represented by \f$(\bar{a},\bar{c},\bar{m})\f$, where the
/// relationship between coordinates is \f$x = s \bar{x} + u\f$, where \f$s\f$ is a vector of scales,
/// and \f$u\f$ is a vector of offsets. Substituting \f$x = s \bar{x} + u\f$ into the equation given
/// in the documentation for XLALSetFlatLatticeBound() gives the transformation rules:
/// \f{eqnarray*}{
/// \bar{a} &=& \frac{a - u}{s} &\quad \bar{m} &=& m \\{}
/// \bar{c}_{j,k} &=& c_{j,k} \prod_{i=0}^{n-1} s_i^{m_{nj+i,k}} &\quad \bar{c}_{N,k} &=& \frac{c_{N,k}}{s_n}
/// \f}
///
static void FLT_TransformBound(
  const size_t n,				///< [in] Dimension on which bound applies (\f$n\f$)
  const gsl_vector* s,				///< [in] Vector of scales (\f$s\f$)
  const gsl_vector* u,				///< [in] Vector of offsets (\f$u\f$)
  const size_t N,				///< [in] Number of coefficients (\f$N\f$)
  gsl_vector* a,				///< [in/out] Vector of offsets (\f$a \rightarrow \bar{a}\f$)
  gsl_matrix* c,				///< [in/out] Matrix of coefficients (\f$c \rightarrow \bar{c}\f$)
  const gsl_matrix* m				///< [in] Matrix of exponents (\f$m = \bar{m}\f$)
  )
{
  if (a != NULL) {
    gsl_vector_sub(a, u);
    gsl_vector_div(a, s);
  }
  if (c != NULL && m != NULL) {
    for (size_t k = 0; k < /* P */ c->size2; ++k) {
      for (size_t j = 0; j < N; ++j) {
        double c_j = gsl_matrix_get(c, j, k);
        for (size_t i = 0; i < n; ++i) {
          c_j *= FLT_BoundPow(gsl_vector_get(s, i), gsl_matrix_get(m, n*j + i, k));
        }
        gsl_matrix_set(c, j, k, c_j);
      }
      double c_N = gsl_matrix_get(c, N, k);
      c_N /= gsl_vector_get(s, n);
      gsl_matrix_set(c, N, k, c_N);
    }
  }
}

///
/// Returns the lower and upper parameter-space bounds
///
static void FLT_GetBounds(
  const FlatLatticeTiling* tiling,		///< [in] Tiling state
  const size_t dimension,			///< [in] Dimension on which bound applies
  const gsl_vector* point,			///< [in] Point at which to find bounds
  double* lower,				///< [out] Lower bound on point
  double* upper					///< [out] Upper bound on point
  )
{

  // Calculate lower and upper bounds
  const FLT_Bound* bound = &tiling->bounds[dimension];
  if (dimension == 0) {
    *lower = FLT_GetBound(dimension, NULL, bound->N, bound->a, bound->c_lower, bound->m_lower);
    *upper = FLT_GetBound(dimension, NULL, bound->N, bound->a, bound->c_upper, bound->m_upper);
  } else {
    gsl_vector_const_view point_n = gsl_vector_const_subvector(point, 0, dimension);
    *lower = FLT_GetBound(dimension, &point_n.vector, bound->N, bound->a, bound->c_lower, bound->m_lower);
    *upper = FLT_GetBound(dimension, &point_n.vector, bound->N, bound->a, bound->c_upper, bound->m_upper);
  }

}

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
  tiling->lattice = FLAT_LATTICE_TYPE_MAX;

  // Allocate parameter-space bounds info
  tiling->bounds = XLALCalloc(n, sizeof(FLT_Bound));
  XLAL_CHECK_NULL(tiling->bounds != NULL, XLAL_ENOMEM);

  // Allocate vectors and matrices
  tiling->phys_scale = gsl_vector_alloc(n);
  XLAL_CHECK_NULL(tiling->phys_scale != NULL, XLAL_ENOMEM);
  tiling->phys_offset = gsl_vector_alloc(n);
  XLAL_CHECK_NULL(tiling->phys_offset != NULL, XLAL_ENOMEM);
  tiling->metric = gsl_matrix_alloc(n, n);
  XLAL_CHECK_NULL(tiling->metric != NULL, XLAL_ENOMEM);
  tiling->increment = gsl_matrix_alloc(n, n);
  XLAL_CHECK_NULL(tiling->increment != NULL, XLAL_ENOMEM);
  tiling->padding = gsl_vector_alloc(n);
  XLAL_CHECK_NULL(tiling->padding != NULL, XLAL_ENOMEM);
  tiling->point = gsl_vector_alloc(n);
  XLAL_CHECK_NULL(tiling->point != NULL, XLAL_ENOMEM);
  tiling->lower = gsl_vector_alloc(n);
  XLAL_CHECK_NULL(tiling->lower != NULL, XLAL_ENOMEM);
  tiling->upper = gsl_vector_alloc(n);
  XLAL_CHECK_NULL(tiling->upper != NULL, XLAL_ENOMEM);

  return tiling;

}

void XLALDestroyFlatLatticeTiling(
  FlatLatticeTiling* tiling
  )
{

  if (tiling) {

    const size_t n = tiling->dimensions;

    // Free bounds data
    for (size_t i = 0; i < n; ++i) {
      gsl_vector_free(tiling->bounds[i].a);
      gsl_matrix_free(tiling->bounds[i].c_lower);
      gsl_matrix_free(tiling->bounds[i].m_lower);
      gsl_matrix_free(tiling->bounds[i].c_upper);
      gsl_matrix_free(tiling->bounds[i].m_upper);
    }
    XLALFree(tiling->bounds);

    // Free vectors and matrices
    gsl_vector_uint_free(tiling->tiled_idx);
    gsl_vector_free(tiling->phys_scale);
    gsl_vector_free(tiling->phys_offset);
    gsl_matrix_free(tiling->metric);
    gsl_matrix_free(tiling->increment);
    gsl_vector_free(tiling->padding);
    gsl_vector_free(tiling->point);
    gsl_vector_free(tiling->lower);
    gsl_vector_free(tiling->upper);

    // Free tiling structure
    XLALFree(tiling);

  }

}

size_t XLALGetFlatLatticeDimensions(
  const FlatLatticeTiling* tiling
  )
{

  // Check tiling
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);

  return tiling->dimensions;

}

uint64_t XLALGetFlatLatticePointCount(
  const FlatLatticeTiling* tiling
  )
{

  // Check tiling
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_VAL(0, tiling->status >= FLT_S_INITIALISED, XLAL_EINVAL);

  return tiling->count;

}

uint64_t XLALCountFlatLatticePoints(
  FlatLatticeTiling* tiling
  )
{

  // Check tiling
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_VAL(0, tiling->status > FLT_S_INCOMPLETE, XLAL_EINVAL);

  // If points have not already been counted, count them
  if (tiling->total_count == 0) {

    // Iterate over all points
    while (XLALNextFlatLatticePoint(tiling, NULL) >= 0) {
      XLAL_CHECK_VAL(0, XLALFastForwardFlatLatticeTiling(tiling, NULL, NULL) == XLAL_SUCCESS, XLAL_EFUNC);
    }
    XLAL_CHECK_VAL(0, xlalErrno == 0, XLAL_EFUNC);
    XLAL_CHECK_VAL(0, tiling->total_count > 0, XLAL_EFAILED);

    // Restart tiling
    XLAL_CHECK_VAL(0, XLALRestartFlatLatticeTiling(tiling) == XLAL_SUCCESS, XLAL_EFUNC);

  }

  // Return the point count
  return tiling->total_count;

}

int XLALSetFlatLatticeBound(
  FlatLatticeTiling* tiling,
  const size_t dimension,
  const gsl_vector* a,
  const gsl_matrix* c_lower,
  const gsl_matrix* m_lower,
  const gsl_matrix* c_upper,
  const gsl_matrix* m_upper
  )
{

  // Check tiling
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status == FLT_S_INCOMPLETE, XLAL_EINVAL);

  // Check input
  XLAL_CHECK(dimension < tiling->dimensions, XLAL_ESIZE);
  XLAL_CHECK(a != NULL, XLAL_EFAULT);
  XLAL_CHECK(c_lower != NULL, XLAL_EFAULT);
  XLAL_CHECK(m_lower != NULL, XLAL_EFAULT);
  XLAL_CHECK(c_upper != NULL, XLAL_EFAULT);
  XLAL_CHECK(m_upper != NULL, XLAL_EFAULT);
  XLAL_CHECK(a->size == dimension + 1, XLAL_ESIZE,
             "'a' should have length %u", dimension + 1);
  XLAL_CHECK(c_lower->size1 == c_upper->size1, XLAL_ESIZE);
  XLAL_CHECK(c_lower->size2 == c_upper->size2, XLAL_ESIZE);
  XLAL_CHECK(m_lower->size1 == m_upper->size1, XLAL_ESIZE);
  XLAL_CHECK(m_lower->size2 == m_upper->size2, XLAL_ESIZE);
  XLAL_CHECK(c_lower->size1 == c_upper->size1, XLAL_ESIZE);
  XLAL_CHECK(c_lower->size2 == c_upper->size2, XLAL_ESIZE);
  XLAL_CHECK(m_lower->size1 == m_upper->size1, XLAL_ESIZE);
  XLAL_CHECK(m_lower->size2 == m_upper->size2, XLAL_ESIZE);
  XLAL_CHECK(c_lower->size2 == m_lower->size2, XLAL_ESIZE);
  XLAL_CHECK(c_upper->size2 == m_upper->size2, XLAL_ESIZE);
  const size_t N = c_lower->size1 - 1;
  XLAL_CHECK(m_lower->size1 == dimension * N + 1, XLAL_ESIZE,
             "'m_{lower|upper}' should have %u rows", dimension * N + 1);

  // Allocate memory
  FLT_Bound* bound = &tiling->bounds[dimension];
  bound->a = gsl_vector_alloc(a->size);
  XLAL_CHECK(bound->a != NULL, XLAL_EFAULT);
  bound->c_lower = gsl_matrix_alloc(c_lower->size1, c_lower->size2);
  XLAL_CHECK(bound->c_lower != NULL, XLAL_EFAULT);
  bound->m_lower = gsl_matrix_alloc(m_lower->size1, m_lower->size2);
  XLAL_CHECK(bound->m_lower != NULL, XLAL_EFAULT);
  bound->c_upper = gsl_matrix_alloc(c_upper->size1, c_upper->size2);
  XLAL_CHECK(bound->c_upper != NULL, XLAL_EFAULT);
  bound->m_upper = gsl_matrix_alloc(m_upper->size1, m_upper->size2);
  XLAL_CHECK(bound->m_upper != NULL, XLAL_EFAULT);

  // Determine if bound is tiled, i.e. not a single point
  bound->tiled = false;
  for (size_t i = 0; i < c_lower->size1; ++i) {
    for (size_t j = 0; j < c_lower->size2; ++j) {
      bound->tiled |= gsl_matrix_get(c_lower, i, j) != gsl_matrix_get(c_upper, i, j);
    }
  }
  for (size_t i = 0; i < m_lower->size1; ++i) {
    for (size_t j = 0; j < m_lower->size2; ++j) {
      bound->tiled |= gsl_matrix_get(m_lower, i, j) != gsl_matrix_get(m_upper, i, j);
    }
  }

  // Set the parameter-space bound
  bound->N = N;
  gsl_vector_memcpy(bound->a, a);
  gsl_matrix_memcpy(bound->c_lower, c_lower);
  gsl_matrix_memcpy(bound->m_lower, m_lower);
  gsl_matrix_memcpy(bound->c_upper, c_upper);
  gsl_matrix_memcpy(bound->m_upper, m_upper);

  return XLAL_SUCCESS;

}

int XLALSetFlatLatticeConstantBound(
  FlatLatticeTiling* tiling,
  const size_t dimension,
  const double bound1,
  const double bound2
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(isfinite(bound1), XLAL_EINVAL);
  XLAL_CHECK(isfinite(bound2), XLAL_EINVAL);

  // Allocate zeroed memory
  gsl_vector* a = gsl_vector_calloc(dimension + 1);
  XLAL_CHECK(a != NULL, XLAL_EFAULT);
  gsl_matrix* c_lower = gsl_matrix_calloc(1, 1);
  XLAL_CHECK(c_lower != NULL, XLAL_EFAULT);
  gsl_matrix* m_lower = gsl_matrix_calloc(1, 1);
  XLAL_CHECK(m_lower != NULL, XLAL_EFAULT);
  gsl_matrix* c_upper = gsl_matrix_calloc(c_lower->size1, c_lower->size2);
  XLAL_CHECK(c_upper != NULL, XLAL_EFAULT);
  gsl_matrix* m_upper = gsl_matrix_calloc(m_lower->size1, m_lower->size2);
  XLAL_CHECK(m_upper != NULL, XLAL_EFAULT);

  // Set bounds data: min(bound1,bound2) <= x <= max(bound1,bound2)
  gsl_matrix_set(c_lower, 0, 0, GSL_MIN(bound1, bound2));
  gsl_matrix_set(c_upper, 0, 0, GSL_MAX(bound1, bound2));

  // Set parameter-space bound
  XLAL_CHECK(XLALSetFlatLatticeBound(tiling, dimension, a, c_lower, m_lower, c_upper, m_upper) == XLAL_SUCCESS, XLAL_EFUNC);

  // Cleanup
  gsl_vector_free(a);
  gsl_matrix_free(c_lower);
  gsl_matrix_free(m_lower);
  gsl_matrix_free(c_upper);
  gsl_matrix_free(m_upper);

  return XLAL_SUCCESS;

}

int XLALSetFlatLatticeEllipticalBounds(
  FlatLatticeTiling* tiling,
  const size_t x_dimension,
  const double x_centre,
  const double y_centre,
  const double x_semi,
  const double y_semi
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(x_semi >= 0.0, XLAL_EINVAL);
  XLAL_CHECK(y_semi >= 0.0, XLAL_EINVAL);
  const size_t y_dimension = x_dimension + 1;

  // Set parameter-space x bound
  XLAL_CHECK(XLALSetFlatLatticeConstantBound(tiling, x_dimension, x_centre - x_semi, x_centre + x_semi) == XLAL_SUCCESS, XLAL_EFUNC);

  // Allocate zeroed memory
  gsl_vector* a = gsl_vector_calloc(y_dimension + 1);
  XLAL_CHECK(a != NULL, XLAL_EFAULT);
  gsl_matrix* c_lower = gsl_matrix_calloc(3, 1);
  XLAL_CHECK(c_lower != NULL, XLAL_EFAULT);
  gsl_matrix* m_lower = gsl_matrix_calloc(2*y_dimension + 1, 1);
  XLAL_CHECK(m_lower != NULL, XLAL_EFAULT);
  gsl_matrix* c_upper = gsl_matrix_calloc(c_lower->size1, c_lower->size2);
  XLAL_CHECK(c_upper != NULL, XLAL_EFAULT);
  gsl_matrix* m_upper = gsl_matrix_calloc(m_lower->size1, m_lower->size2);
  XLAL_CHECK(m_upper != NULL, XLAL_EFAULT);

  // Set bounds data: -sqrt(1 - (x-x_centre)^2) <= y-y_centre <= +sqrt(1 - (x-x_centre))
  gsl_vector_set(a, y_dimension - 1, x_centre);
  gsl_vector_set(a, y_dimension, y_centre);
  gsl_matrix_set(c_lower, 0, 0, 1);
  gsl_matrix_set(c_lower, 1, 0, -1);
  gsl_matrix_set(c_lower, 2, 0, -1);
  gsl_matrix_set(c_upper, 0, 0, 1);
  gsl_matrix_set(c_upper, 1, 0, -1);
  gsl_matrix_set(c_upper, 2, 0, 1);
  gsl_matrix_set(m_lower, 2*y_dimension - 1, 0, 2);
  gsl_matrix_set(m_lower, 2*y_dimension, 0, 0.5);
  gsl_matrix_memcpy(m_upper, m_lower);

  // Set parameter-space y bound
  XLAL_CHECK(XLALSetFlatLatticeBound(tiling, y_dimension, a, c_lower, m_lower, c_upper, m_upper) == XLAL_SUCCESS, XLAL_EFUNC);

  // Cleanup
  gsl_vector_free(a);
  gsl_matrix_free(c_lower);
  gsl_matrix_free(m_lower);
  gsl_matrix_free(c_upper);
  gsl_matrix_free(m_upper);

  return XLAL_SUCCESS;

}

int XLALSetFlatLatticeTypeAndMetric(
  FlatLatticeTiling* tiling,
  const FlatLatticeType lattice,
  const gsl_matrix* metric,
  const double max_mismatch
  )
{

  // Check tiling
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status == FLT_S_INCOMPLETE, XLAL_EINVAL);
  const size_t n = tiling->dimensions;

  // Check input
  XLAL_CHECK(lattice < FLAT_LATTICE_TYPE_MAX, XLAL_EINVAL, "'lattice'=%u must be in [0,%u)", lattice, FLAT_LATTICE_TYPE_MAX);
  XLAL_CHECK(metric != NULL, XLAL_EFAULT);
  XLAL_CHECK(metric->size1 == n && metric->size2 == n, XLAL_EINVAL);
  XLAL_CHECK(max_mismatch > 0, XLAL_EINVAL);

  // Check that all parameter-space dimensions are bounded, and record indices of tiled dimensions
  size_t tn = 0;
  for (size_t i = 0; i < tiling->dimensions; ++i) {
    XLAL_CHECK(tiling->bounds[i].a != NULL, XLAL_EFAILED, "Dimension #%i is unbounded", i);
    if (tiling->bounds[i].tiled) {
      ++tn;
    }
  }
  if (tn > 0) {
    tiling->tiled_idx = gsl_vector_uint_alloc(tn);
    for (size_t i = 0, ti = 0; i < tiling->dimensions; ++i) {
      if (tiling->bounds[i].tiled) {
        gsl_vector_uint_set(tiling->tiled_idx, ti, i);
        ++ti;
      }
    }
  }

  // Save the type of lattice to generate flat tiling with
  tiling->lattice = lattice;

  // Check diagonal elements of tiled dimensions are positive, and calculate
  // physical parameter-space scale from metric diagonal elements
  gsl_vector_set_all(tiling->phys_scale, 1.0);
  for (size_t ti = 0; ti < tn; ++ti) {
    const size_t i = gsl_vector_uint_get(tiling->tiled_idx, ti);
    const double metric_i_i = gsl_matrix_get(metric, i, i);
    XLAL_CHECK(metric_i_i > 0, XLAL_EINVAL, "metric(%zu,%zu) <= 0", i, i);
    gsl_vector_set(tiling->phys_scale, i, 1.0 / sqrt(metric_i_i));
  }

  // Set physical parameter-space offset
  for (size_t i = 0; i < n; ++i) {
    double phys_lower = 0, phys_upper = 0;
    FLT_GetBounds(tiling, i, tiling->phys_offset, &phys_lower, &phys_upper);
    gsl_vector_set(tiling->phys_offset, i, phys_lower);
  }

  // Transform parameter-space bounds from physical to normalised coordinates
  for (size_t i = 0; i < n; ++i) {
    const FLT_Bound* bound = &tiling->bounds[i];
    gsl_vector_view phys_scale = gsl_vector_subvector(tiling->phys_scale, 0, i+1);
    gsl_vector_view phys_offset = gsl_vector_subvector(tiling->phys_offset, 0, i+1);
    FLT_TransformBound(i, &phys_scale.vector, &phys_offset.vector, bound->N, bound->a, bound->c_lower, bound->m_lower);
    FLT_TransformBound(i, &phys_scale.vector, &phys_offset.vector, bound->N,     NULL, bound->c_upper, bound->m_upper);
  }

  // Check metric is symmetric, and copy rescaled metric
  for (size_t i = 0; i < n; ++i) {
    const double scale_i = gsl_vector_get(tiling->phys_scale, i);
    for (size_t j = 0; j < n; ++j) {
      const double scale_j = gsl_vector_get(tiling->phys_scale, j);
      double metric_i_j = gsl_matrix_get(metric, i, j);
      XLAL_CHECK(metric_i_j == gsl_matrix_get(metric, j, i), XLAL_EINVAL, "metric(%zu,%zu) != metric(%zu,%zu)", i, j, j, i);
      metric_i_j *= scale_i * scale_j;
      gsl_matrix_set(tiling->metric, i, j, metric_i_j);
    }
  }

  // Initialise for zero-dimensional parameter space
  gsl_matrix_set_zero(tiling->increment);
  gsl_vector_set_zero(tiling->padding);

  // If there are tiled dimensions:
  if (tn > 0) {

    // Allocate memory
    gsl_matrix* tiled_metric = gsl_matrix_alloc(tn, tn);
    XLAL_CHECK(tiled_metric != NULL, XLAL_ENOMEM);

    // Copy tiled dimensions of metric
    for (size_t ti = 0; ti < tn; ++ti) {
      const size_t i = gsl_vector_uint_get(tiling->tiled_idx, ti);
      for (size_t tj = 0; tj < tn; ++tj) {
        const size_t j = gsl_vector_uint_get(tiling->tiled_idx, tj);
        gsl_matrix_set(tiled_metric, ti, tj, gsl_matrix_get(tiling->metric, i, j));
      }
    }

    // Calculate metric lattice increment vectors
    gsl_matrix* increment = FLT_MetricLatticeIncrements(tiling->lattice, tiled_metric, max_mismatch);
    XLAL_CHECK(increment != NULL, XLAL_EFUNC);

    // Calculate metric ellipse bounding box
    gsl_vector* bounding_box = FLT_MetricEllipseBoundingBox(tiled_metric, max_mismatch);
    XLAL_CHECK(bounding_box != NULL, XLAL_EFUNC);

    // Copy increment vectors and padding so that non-tiled dimensions are zero
    for (size_t ti = 0; ti < tn; ++ti) {
      const size_t i = gsl_vector_uint_get(tiling->tiled_idx, ti);
      gsl_vector_set(tiling->padding, i, gsl_vector_get(bounding_box, ti));
      for (size_t tj = 0; tj < tn; ++tj) {
        const size_t j = gsl_vector_uint_get(tiling->tiled_idx, tj);
        gsl_matrix_set(tiling->increment, i, j, gsl_matrix_get(increment, ti, tj));
      }
    }

    // Cleanup
    gsl_matrix_free(tiled_metric);
    gsl_matrix_free(increment);
    gsl_vector_free(bounding_box);

  }

  // Tiling has been fully initialised
  tiling->status = FLT_S_INITIALISED;
  XLAL_CHECK(XLALRestartFlatLatticeTiling(tiling) == XLAL_SUCCESS, XLAL_EFUNC);

  return XLAL_SUCCESS;

}

int XLALNextFlatLatticePoint(
  FlatLatticeTiling* tiling,
  gsl_vector* curr_point
  )
{

  // Check tiling
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status > FLT_S_INCOMPLETE, XLAL_EINVAL);
  const size_t n = tiling->dimensions;
  const size_t tn = (tiling->tiled_idx != NULL) ? tiling->tiled_idx->size : 0;
  XLAL_CHECK(curr_point == NULL || curr_point->size == n, XLAL_EINVAL);

  // If finished status, nothing more to be done!
  if (tiling->status == FLT_S_FINISHED) {
    return -1;
  }

  size_t ti;

  // If initialised status, set and return starting point
  if (tiling->status == FLT_S_INITIALISED) {

    // Set parameter-space bounds and starting point
    for (size_t i = 0; i < n; ++i) {

      // Get normalised bounds
      double lower = 0, upper = 0;
      FLT_GetBounds(tiling, i, tiling->point, &lower, &upper);

      // Set parameter-space bounds
      gsl_vector_set(tiling->lower, i, lower);
      gsl_vector_set(tiling->upper, i, upper);

      // Initialise current point
      const double point = lower - gsl_vector_get(tiling->padding, i);
      gsl_vector_set(tiling->point, i, point);

    }

    // Initialise count
    tiling->count = 1;

    // Tiling has been started
    tiling->status = FLT_S_STARTED;

    // All dimensions of point have changed
    ti = 0;

  } else {

    // Otherwise started status: loop until the next point is found
    ti = tn;
    while (true) {

      // If dimension index is now zero, we're done!
      if (ti == 0) {

        // Tiling is now finished
        tiling->status = FLT_S_FINISHED;

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
      const double point = gsl_vector_get(tiling->point, i);
      const double upper = gsl_vector_get(tiling->upper, i);
      const double padding = gsl_vector_get(tiling->padding, i);
      if (point <= upper + padding) {
        break;
      }

      // Move on to lower dimensions
      continue;

    }

    // Return point to lower bound in higher dimensions
    for (size_t tj = ti + 1; tj < tn; ++tj) {

      // Get index of tiled dimension
      const size_t j = gsl_vector_uint_get(tiling->tiled_idx, tj);

      // Get normalised bounds
      double lower = 0, upper = 0;
      FLT_GetBounds(tiling, j, tiling->point, &lower, &upper);

      // Set parameter-space bounds
      gsl_vector_set(tiling->lower, j, lower);
      gsl_vector_set(tiling->upper, j, upper);

      // Get increment vector
      gsl_vector_view increment = gsl_matrix_column(tiling->increment, j);

      // Calculate the distance from current point to the lower bound, in integer number of increments
      const double padding = gsl_vector_get(tiling->padding, j);
      const double point = gsl_vector_get(tiling->point, j);
      const double dist = ceil( (lower - padding - point) / gsl_vector_get(&increment.vector, j) );

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

int XLALFastForwardFlatLatticeTiling(
  FlatLatticeTiling* tiling,			///< [in] Tiling state
  uint64_t *point_count,			///< [in/out] Count of points fast-forwarded over
  double *point_spacing				///< [in/out] Spacing between points fast-forwarded over
  )
{

  // Check tiling
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status == FLT_S_STARTED, XLAL_EINVAL);
  const size_t tn = (tiling->tiled_idx != NULL) ? tiling->tiled_idx->size : 0;

  // If no tiled dimensions, nothing to fast-forward
  uint64_t count = 0;
  double spacing = 0.0;

  // If there are tiled dimensions:
  if (tn > 0) {

    // Get index of highest tiled dimension
    const size_t i = gsl_vector_uint_get(tiling->tiled_idx, tn - 1);

    // Get current point, spacing, upper bound, and padding
    const double point = gsl_vector_get(tiling->point, i);
    spacing = gsl_matrix_get(tiling->increment, i, i);
    const double upper = gsl_vector_get(tiling->upper, i);
    const double padding = gsl_vector_get(tiling->padding, i);

    // Calculate number of points to fast-forward, so that then calling
    // XLALNextFlatticePoint() will advance the next highest tiled dimension
    count = floor((upper + padding - point) / spacing);

    // Fast-forward over dimension
    gsl_vector_set(tiling->point, i, point + count*spacing);
    tiling->count += count;

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

int XLALRestartFlatLatticeTiling(
  FlatLatticeTiling* tiling
  )
{

  // Check tiling
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status > FLT_S_INCOMPLETE, XLAL_EINVAL);

  // Restart tiling
  tiling->status = FLT_S_INITIALISED;
  tiling->count = 0;

  return XLAL_SUCCESS;

}
