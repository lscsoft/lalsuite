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
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA  02111-1307  USA
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
/// Flat lattice tiling bound info
///
typedef struct tagFLT_Bound {
  bool tiled;                                   ///< Is the bound tiled, i.e. not a single point?
  size_t num_bounds;                            ///< Number of bound pieces in this dimension
  size_t num_coeffs;                            ///< Number of bound coefficients (\f$N\f$)
  gsl_vector* a;                                ///< Vector of offsets (\f$a\f$)
  gsl_matrix* c_lower;                          ///< Column vectors of coefficients (\f$c\f$) for the lower bound of each bound pair
  gsl_matrix* m_lower;                          ///< Column vectors of exponents (\f$m\f$) for the lower bound of each bound pair
  gsl_matrix* c_upper;                          ///< Column vectors of coefficients (\f$c\f$) for the upper bound of each bound pair
  gsl_matrix* m_upper;                          ///< Column vectors of exponents (\f$m\f$) for the upper bound of each bound pair
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
/// Flat lattice tiling state structure
///
struct tagFlatLatticeTiling {
  size_t dimensions;                            ///< Dimension of the parameter space
  size_t tiled_dimensions;                      ///< Tiled dimension of the parameter space
  FLT_Status status;                            ///< Status of the tiling
  FLT_Bound *bounds;                            ///< Array of parameter-space bound info for each dimension
  FlatLatticeType lattice;                      ///< Type of lattice to generate flat tiling with
  gsl_vector* phys_scale;                       ///< Normalised to physical coordinate scaling
  gsl_vector* phys_offset;                      ///< Normalised to physical coordinate offset
  gsl_matrix* metric;                           ///< Normalised parameter-space metric
  gsl_matrix* increment;                        ///< Increment vectors of the lattice tiling generator
  gsl_vector* padding;                          ///< Padding at edges of parameter-space bounds
  gsl_vector* point;                            ///< Current lattice point
  gsl_vector_uint* bound_idx;                   ///< Indices of current bound on parameter space
  gsl_vector* lower;                            ///< Current lower bound on parameter space
  gsl_vector* upper;                            ///< Current upper bound on parameter space
  unsigned long count;                          ///< Number of points generated so far
  unsigned long total_count;                    ///< Total number of points in parameter space
};

///
/// Find the bounding box of the mismatch ellipses of a metric
///
static gsl_vector* FLT_MetricEllipseBoundingBox(
  const gsl_matrix* metric,                     ///< [in] Metric to bound
  const double max_mismatch                     ///< [in] Maximum mismatch with respect to metric
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

///
/// Orthonormalise the columns of a matrix with respect to a metric (matrix is lower triangular)
///
static int FLT_OrthonormaliseWRTMetric(
  gsl_matrix* matrix,                           ///< [in] Matrix of columns to orthonormalise
  const gsl_matrix* metric                      ///< [in] Metric to orthonormalise with respect to
  )
{

  // Check input
  XLAL_CHECK(matrix != NULL, XLAL_EFAULT);
  XLAL_CHECK(metric != NULL, XLAL_EFAULT);
  XLAL_CHECK(metric->size1 == metric->size2, XLAL_ESIZE);
  XLAL_CHECK(metric->size1 == matrix->size2 && metric->size2 == matrix->size2, XLAL_ESIZE);
  const size_t n = metric->size1;

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

///
/// Transform a lattice generator to a square lower triangular form
///
static gsl_matrix* FLT_SquareLowerTriangularLatticeGenerator(
  gsl_matrix* generator                         ///< [in] Generator matrix of lattice
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
  gsl_matrix* generator,                        ///< [in] Generator matrix of lattice
  const double norm_thickness,                  ///< [in] Normalised thickness of lattice
  const double covering_radius                  ///< [in] Desired covering radius
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
  const FlatLatticeType lattice,                ///< [in] Lattice type
  const gsl_matrix* metric,                     ///< [in] parameter-space metric
  const double max_mismatch                     ///< [in] Maximum prescribed mismatch
  )
{

  // Check input
  XLAL_CHECK_NULL(metric != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(metric->size1 == metric->size2, XLAL_ESIZE);
  XLAL_CHECK_NULL(max_mismatch > 0.0, XLAL_EINVAL);
  const size_t r = metric->size1;

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
  XLAL_CHECK_NULL(FLT_OrthonormaliseWRTMetric(directions, metric) == XLAL_SUCCESS, XLAL_EFAILED);

  // Get lattice generator
  gsl_matrix* gen_matrix = NULL;
  double norm_thickness = 0.0;
  switch (lattice) {

  case FLAT_LATTICE_TYPE_CUBIC:   // Cubic (\f$Z_n\f$) lattice

    // Allocate memory
    gen_matrix = gsl_matrix_alloc(r, r);
    XLAL_CHECK_NULL(gen_matrix != NULL, XLAL_ENOMEM);

    // Create generator
    gsl_matrix_set_identity(gen_matrix);

    // Calculate normalised thickness
    norm_thickness = pow(sqrt(r)/2, r);

    break;

  case FLAT_LATTICE_TYPE_ANSTAR:   // An-star (\f$A_n^*\f$) lattice

    // Allocate memory
    gen_matrix = gsl_matrix_alloc(r + 1, r);
    XLAL_CHECK_NULL(gen_matrix != NULL, XLAL_ENOMEM);

    // Create generator in (r + 1) space
    gsl_matrix_set_all(gen_matrix, 0.0);
    {
      gsl_vector_view first_row = gsl_matrix_row(gen_matrix, 0);
      gsl_vector_view sub_diag = gsl_matrix_subdiagonal(gen_matrix, 1);
      gsl_vector_view last_col = gsl_matrix_column(gen_matrix, r - 1);
      gsl_vector_set_all(&first_row.vector, 1.0);
      gsl_vector_set_all(&sub_diag.vector, -1.0);
      gsl_vector_set_all(&last_col.vector, 1.0 / (r + 1.0));
      gsl_vector_set(&last_col.vector, 0, -1.0 * r / (r + 1.0));
    }

    // Calculate normalised thickness
    norm_thickness = sqrt(r + 1.0)*pow((1.0*r*(r + 2))/(12.0*(r + 1)), 0.5*r);

    break;

  default:
    XLAL_ERROR_NULL(XLAL_EINVAL, "Invalid 'lattice'=%u", lattice);
  }

  // Transform lattice generator to square lower triangular
  gsl_matrix* sqlwtr_gen_matrix = FLT_SquareLowerTriangularLatticeGenerator(gen_matrix);
  XLAL_CHECK_NULL(sqlwtr_gen_matrix != NULL, XLAL_EFAILED);

  // Normalise lattice generator so covering radius is sqrt(mismatch)
  XLAL_CHECK_NULL(FLT_NormaliseLatticeGenerator(sqlwtr_gen_matrix, norm_thickness, sqrt(max_mismatch)) == XLAL_SUCCESS, XLAL_EFAILED);

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
static inline double FLT_pow(double x, double y) {
  return y == 0.0 ? 1.0 : creal(cpow(x, y));
}

///
/// Returns the parameter-space bound \f$b_n\f$ on dimension \f$n\f$, given the current point
/// \f$x = (x_0,\cdots,x_{n-1})\f$ in the lower dimensions.
///
static double FLT_GetBound(
  const size_t dimension,                       ///< [in] Dimension on which bound applies (\f$n\f$)
  const gsl_vector* point,                      ///< [in] Point at which to find bounds (\f$x\f$)
  const size_t num_coeffs,                      ///< [in] Number of coefficients (\f$N\f$)
  const gsl_vector* a,                          ///< [in] Vector of offsets (\f$a\f$)
  const gsl_vector* c,                          ///< [in] Vector of coefficients (\f$c\f$)
  const gsl_vector* m                           ///< [in] Vector of exponents (\f$m\f$)
  )
{
  double bound = 0;
  for (size_t i = 0; i < num_coeffs; ++i) {
    double bound_i = gsl_vector_get(c, i);
    if (dimension > 0) {
      for (size_t j = 0; j < dimension; ++j) {
        const double dx_j = gsl_vector_get(point, j) - gsl_vector_get(a, j);
        bound_i *= FLT_pow(dx_j, gsl_vector_get(m, dimension*i + j));
      }
    }
    bound += bound_i;
  }
  bound = FLT_pow(bound, gsl_vector_get(m, dimension*num_coeffs));
  bound = gsl_vector_get(c, num_coeffs) * bound + gsl_vector_get(a, dimension);
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
/// \bar{c}_i &=& c_i \prod_{j=0}^{n-1} s_j^{m_{ni+j}} &\quad \bar{c}_N &=& \frac{c_N}{s_n}
/// \f}
///
static void FLT_TransformBound(
  const size_t dimension,                       ///< [in] Dimension on which bound applies (\f$n\f$)
  const gsl_vector* s,                          ///< [in] Vector of scales (\f$s\f$)
  const gsl_vector* u,                          ///< [in] Vector of offsets (\f$u\f$)
  const size_t num_coeffs,                      ///< [in] Number of coefficients (\f$N\f$)
  gsl_vector* a,                                ///< [in/out] Vector of offsets (\f$a \rightarrow \bar{a}\f$)
  gsl_vector* c,                                ///< [in/out] Vector of coefficients (\f$c \rightarrow \bar{c}\f$)
  const gsl_vector* m                           ///< [in] Vector of exponents (\f$m = \bar{m}\f$)
  )
{
  if (a != NULL) {
    gsl_vector_sub(a, u);
    gsl_vector_div(a, s);
  }
  if (c != NULL && m != NULL) {
    for (size_t i = 0; i < num_coeffs; ++i) {
      double c_i = gsl_vector_get(c, i);
      for (size_t j = 0; j < dimension; ++j) {
        c_i *= FLT_pow(gsl_vector_get(s, j), gsl_vector_get(m, dimension*i + j));
      }
      gsl_vector_set(c, i, c_i);
    }
    double c_N = gsl_vector_get(c, num_coeffs);
    c_N /= gsl_vector_get(s, dimension);
    gsl_vector_set(c, num_coeffs, c_N);
  }
}

///
/// Returns the lower and upper parameter-space bounds
///
static void FLT_GetBounds(
  const FlatLatticeTiling* tiling,              ///< [in] Tiling state
  const size_t dimension,                       ///< [in] Dimension on which bound applies
  const size_t bound,                           ///< [in] Index of bound within this dimension
  const gsl_vector* point,                      ///< [in] Point at which to find bounds
  double* lower,                                ///< [out] Lower bound on point
  double* upper                                 ///< [out] Upper bound on point
  )
{

  // Get views of c and m vectors for this bound
  gsl_vector_const_view c_lower = gsl_matrix_const_column(tiling->bounds[dimension].c_lower, bound);
  gsl_vector_const_view m_lower = gsl_matrix_const_column(tiling->bounds[dimension].m_lower, bound);
  gsl_vector_const_view c_upper = gsl_matrix_const_column(tiling->bounds[dimension].c_upper, bound);
  gsl_vector_const_view m_upper = gsl_matrix_const_column(tiling->bounds[dimension].m_upper, bound);

  // Calculate lower and upper bounds
  const size_t num_coeffs = tiling->bounds[dimension].num_coeffs;
  const gsl_vector* a = tiling->bounds[dimension].a;
  if (dimension == 0) {
    *lower = FLT_GetBound(dimension, NULL, num_coeffs, a, &c_lower.vector, &m_lower.vector);
    *upper = FLT_GetBound(dimension, NULL, num_coeffs, a, &c_upper.vector, &m_upper.vector);
  } else {
    gsl_vector_const_view point_n = gsl_vector_const_subvector(point, 0, dimension);
    *lower = FLT_GetBound(dimension, &point_n.vector, num_coeffs, a, &c_lower.vector, &m_lower.vector);
    *upper = FLT_GetBound(dimension, &point_n.vector, num_coeffs, a, &c_upper.vector, &m_upper.vector);
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
  tiling->bound_idx = gsl_vector_uint_alloc(n);
  XLAL_CHECK_NULL(tiling->bound_idx != NULL, XLAL_ENOMEM);
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

    gsl_error_handler_t* old_handler = gsl_set_error_handler_off();

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
    gsl_vector_free(tiling->phys_scale);
    gsl_vector_free(tiling->phys_offset);
    gsl_matrix_free(tiling->metric);
    gsl_matrix_free(tiling->increment);
    gsl_vector_free(tiling->padding);
    gsl_vector_free(tiling->point);
    gsl_vector_uint_free(tiling->bound_idx);
    gsl_vector_free(tiling->lower);
    gsl_vector_free(tiling->upper);

    // Free tiling structure
    XLALFree(tiling);

    gsl_set_error_handler(old_handler);

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

gsl_vector* XLALGetFlatLatticePoint(
  const FlatLatticeTiling* tiling,
  gsl_vector* point
  )
{

  // Check tiling
  XLAL_CHECK_NULL(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(tiling->status > FLT_S_INITIALISED, XLAL_EFAILED);

  // Allocate physical point vector if needed
  if (point == NULL) {
    point = gsl_vector_alloc(tiling->dimensions);
    XLAL_CHECK_NULL(point != NULL, XLAL_ENOMEM);
  }

  // Copy current point and transform from normalised to physical coordinates
  gsl_vector_memcpy(point, tiling->point);
  gsl_vector_mul(point, tiling->phys_scale);
  gsl_vector_add(point, tiling->phys_offset);

  return point;

}

unsigned long XLALGetFlatLatticePointCount(
  const FlatLatticeTiling* tiling
  )
{

  // Check tiling
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_VAL(0, tiling->status >= FLT_S_INITIALISED, XLAL_EFAILED);

  return tiling->count;

}

unsigned long XLALCountFlatLatticePoints(
  FlatLatticeTiling* tiling
  )
{

  // Check tiling
  XLAL_CHECK_VAL(0, tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_VAL(0, tiling->status > FLT_S_INCOMPLETE, XLAL_EFAILED);

  // If templates have not already been counted, count them
  if (tiling->total_count == 0) {

    // Iterate over all templates
    while (XLALNextFlatLatticePoint(tiling) >= 0);
    XLAL_CHECK_VAL(0, xlalErrno == 0, XLAL_EFUNC);
    XLAL_CHECK_VAL(0, tiling->total_count > 0, XLAL_EFAILED);

    // Restart tiling
    XLALRestartFlatLatticeTiling(tiling);

  }

  // Return the template count
  return tiling->total_count;

}

gsl_matrix* XLALGetFlatLatticeIncrements(
  const FlatLatticeTiling* tiling
  )
{

  // Check tiling
  XLAL_CHECK_NULL(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(tiling->status >= FLT_S_INITIALISED, XLAL_EFAILED);

  // Allocate increment vector
  gsl_matrix* increment = gsl_matrix_alloc(tiling->increment->size1, tiling->increment->size2);
  XLAL_CHECK_NULL(tiling->increment != NULL, XLAL_ENOMEM);

  // Copy increments, rescaled to physical coordinates
  for (size_t i = 0; i < increment->size2; ++i) {
    gsl_vector_view tiling_increment_i = gsl_matrix_column(tiling->increment, i);
    gsl_vector_view increment_i = gsl_matrix_column(increment, i);
    gsl_vector_memcpy(&increment_i.vector, &tiling_increment_i.vector);
    gsl_vector_mul(&increment_i.vector, tiling->phys_scale);
  }

  return increment;

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
  XLAL_CHECK(tiling->status == FLT_S_INCOMPLETE, XLAL_EFAILED);

  // Check input
  XLAL_CHECK(dimension < tiling->dimensions, XLAL_ESIZE);
  XLAL_CHECK(a != NULL, XLAL_EFAULT);
  XLAL_CHECK(c_lower != NULL, XLAL_EFAULT);
  XLAL_CHECK(m_lower != NULL, XLAL_EFAULT);
  XLAL_CHECK(c_upper != NULL, XLAL_EFAULT);
  XLAL_CHECK(m_upper != NULL, XLAL_EFAULT);
  XLAL_CHECK(a->size == dimension + 1, XLAL_EFAULT);
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

  // Determine number of bounds, and number of coefficients and exponents per bound
  const size_t num_bounds = c_lower->size2;
  const size_t num_coeffs = c_lower->size1 - 1;
  const size_t num_exponents = dimension * num_coeffs + 1;
  XLAL_CHECK(num_exponents == m_lower->size1, XLAL_ESIZE,
             "m_{lower|upper} has %u rows, should have %u", m_lower->size1, num_exponents);

  // Determine if bound is tiled, i.e. not a single point
  bool tiled = false;
  for (size_t i = 0; i < c_lower->size1; ++i) {
    for (size_t j = 0; j < c_lower->size2; ++j) {
      tiled |= gsl_matrix_get(c_lower, i, j) != gsl_matrix_get(c_upper, i, j);
    }
  }
  for (size_t i = 0; i < m_lower->size1; ++i) {
    for (size_t j = 0; j < m_lower->size2; ++j) {
      tiled |= gsl_matrix_get(m_lower, i, j) != gsl_matrix_get(m_upper, i, j);
    }
  }

  // Allocate memory
  tiling->bounds[dimension].a = gsl_vector_alloc(a->size);
  XLAL_CHECK(tiling->bounds[dimension].a != NULL, XLAL_EFAULT);
  tiling->bounds[dimension].c_lower = gsl_matrix_alloc(c_lower->size1, c_lower->size2);
  XLAL_CHECK(tiling->bounds[dimension].c_lower != NULL, XLAL_EFAULT);
  tiling->bounds[dimension].m_lower = gsl_matrix_alloc(m_lower->size1, m_lower->size2);
  XLAL_CHECK(tiling->bounds[dimension].m_lower != NULL, XLAL_EFAULT);
  tiling->bounds[dimension].c_upper = gsl_matrix_alloc(c_upper->size1, c_upper->size2);
  XLAL_CHECK(tiling->bounds[dimension].c_upper != NULL, XLAL_EFAULT);
  tiling->bounds[dimension].m_upper = gsl_matrix_alloc(m_upper->size1, m_upper->size2);
  XLAL_CHECK(tiling->bounds[dimension].m_upper != NULL, XLAL_EFAULT);

  // Set the parameter-space bound
  tiling->bounds[dimension].tiled = tiled;
  tiling->bounds[dimension].num_bounds = num_bounds;
  tiling->bounds[dimension].num_coeffs = num_coeffs;
  gsl_vector_memcpy(tiling->bounds[dimension].a, a);
  gsl_matrix_memcpy(tiling->bounds[dimension].c_lower, c_lower);
  gsl_matrix_memcpy(tiling->bounds[dimension].m_lower, m_lower);
  gsl_matrix_memcpy(tiling->bounds[dimension].c_upper, c_upper);
  gsl_matrix_memcpy(tiling->bounds[dimension].m_upper, m_upper);

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
  XLAL_CHECK(XLALSetFlatLatticeBound(tiling, dimension, a, c_lower, m_lower, c_upper, m_upper) == XLAL_SUCCESS, XLAL_EFAILED);

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
  XLAL_CHECK(XLALSetFlatLatticeConstantBound(tiling, x_dimension, x_centre - x_semi, x_centre + x_semi) == XLAL_SUCCESS, XLAL_EFAILED);

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
  XLAL_CHECK(XLALSetFlatLatticeBound(tiling, y_dimension, a, c_lower, m_lower, c_upper, m_upper) == XLAL_SUCCESS, XLAL_EFAILED);

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
  XLAL_CHECK(tiling->status == FLT_S_INCOMPLETE, XLAL_EFAILED);
  const size_t n = tiling->dimensions;

  // Check input
  XLAL_CHECK(lattice < FLAT_LATTICE_TYPE_MAX, XLAL_EINVAL, "'lattice'=%u must be in [0,%u)", lattice, FLAT_LATTICE_TYPE_MAX);
  XLAL_CHECK(metric != NULL, XLAL_EFAULT);
  XLAL_CHECK(metric->size1 == n && metric->size2 == n, XLAL_EINVAL);
  XLAL_CHECK(max_mismatch > 0, XLAL_EINVAL);

  // Check that all parameter-space dimensions are bounded, and count number of tiles dimensions
  tiling->tiled_dimensions = 0;
  for (size_t i = 0; i < tiling->dimensions; ++i) {
    XLAL_CHECK(tiling->bounds[i].num_bounds > 0, XLAL_EFAILED, "Dimension #%i is unbounded", i);
    tiling->tiled_dimensions += tiling->bounds[i].tiled ? 1 : 0;
  }

  // Save the type of lattice to generate flat tiling with
  tiling->lattice = lattice;

  // Initialise parameter-space bound indices
  gsl_vector_uint_set_zero(tiling->bound_idx);

  // Get physical parameter-space offset
  for (size_t i = 0; i < n; ++i) {

    // Get physical bounds at lowest bound index
    double phys_lower = 0, phys_upper = 0;
    FLT_GetBounds(tiling, i, 0, tiling->phys_offset, &phys_lower, &phys_upper);

    // Set physical parameter-space offset
    gsl_vector_set(tiling->phys_offset, i, phys_lower);

  }

  // Check diagonal elements of tiled dimensions are positive, and calculate
  // physical parameter-space scaling from metric diagonal elements
  gsl_vector_set_all(tiling->phys_scale, 1.0);
  for (size_t i = 0; i < n; ++i) {
    if (tiling->bounds[i].tiled) {
      const double metric_i_i = gsl_matrix_get(metric, i, i);
      XLAL_CHECK(metric_i_i > 0, XLAL_EINVAL, "metric(%zu,%zu) <= 0", i, i);
      gsl_vector_set(tiling->phys_scale, i, 1.0 / sqrt(metric_i_i));
    }
  }

  // Transform parameter-space bounds from physical to normalised coordinates
  for (size_t i = 0; i < n; ++i) {
    gsl_vector_view phys_scale = gsl_vector_subvector(tiling->phys_scale, 0, i+1);
    gsl_vector_view phys_offset = gsl_vector_subvector(tiling->phys_offset, 0, i+1);
    const size_t num_coeffs = tiling->bounds[i].num_coeffs;
    FLT_TransformBound(i, &phys_scale.vector, &phys_offset.vector, num_coeffs, tiling->bounds[i].a, NULL, NULL);
    for (size_t bound = 0; bound < tiling->bounds[i].num_bounds; ++bound) {
      gsl_vector_view c_lower = gsl_matrix_column(tiling->bounds[i].c_lower, bound);
      gsl_vector_view m_lower = gsl_matrix_column(tiling->bounds[i].m_lower, bound);
      gsl_vector_view c_upper = gsl_matrix_column(tiling->bounds[i].c_upper, bound);
      gsl_vector_view m_upper = gsl_matrix_column(tiling->bounds[i].m_upper, bound);
      FLT_TransformBound(i, &phys_scale.vector, &phys_offset.vector, num_coeffs, NULL, &c_lower.vector, &m_lower.vector);
      FLT_TransformBound(i, &phys_scale.vector, &phys_offset.vector, num_coeffs, NULL, &c_upper.vector, &m_upper.vector);
    }
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

  if (tiling->tiled_dimensions > 0) {

    const size_t tn = tiling->tiled_dimensions;

    // Allocate memory
    gsl_matrix* tiled_metric = gsl_matrix_alloc(tn, tn);
    XLAL_CHECK(tiled_metric != NULL, XLAL_ENOMEM);

    // Copy tiled dimensions of metric
    for (size_t i = 0, ti = 0; i < n; ++i) {
      if (tiling->bounds[i].tiled) {
        for (size_t j = 0, tj = 0; j < n; ++j) {
          if (tiling->bounds[j].tiled) {
            gsl_matrix_set(tiled_metric, ti, tj, gsl_matrix_get(tiling->metric, i, j));
            ++tj;
          }
        }
        ++ti;
      }
    }

    // Calculate metric lattice increment vectors
    gsl_matrix* increment = FLT_MetricLatticeIncrements(tiling->lattice, tiled_metric, max_mismatch);
    XLAL_CHECK(increment != NULL, XLAL_EFAILED);

    // Calculate metric ellipse bounding box
    gsl_vector* bounding_box = FLT_MetricEllipseBoundingBox(tiled_metric, max_mismatch);
    XLAL_CHECK(bounding_box != NULL, XLAL_EFAILED);

    // Copy increment vectors and padding so that non-tiled dimensions are zero
    for (size_t i = 0, ti = 0; i < n; ++i) {
      if (tiling->bounds[i].tiled) {
        gsl_vector_set(tiling->padding, i, gsl_vector_get(bounding_box, ti));
        for (size_t j = 0, tj = 0; j < n; ++j) {
          if (tiling->bounds[j].tiled) {
            gsl_matrix_set(tiling->increment, i, j, gsl_matrix_get(increment, ti, tj));
            ++tj;
          }
        }
        ++ti;
      }
    }

    // Cleanup
    gsl_matrix_free(tiled_metric);
    gsl_matrix_free(increment);
    gsl_vector_free(bounding_box);

  }

  // Tiling has been fully initialised
  tiling->status = FLT_S_INITIALISED;
  XLALRestartFlatLatticeTiling(tiling);

  return XLAL_SUCCESS;

}

int XLALNextFlatLatticePoint(
  FlatLatticeTiling* tiling
  )
{

  // Check tiling
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status > FLT_S_INCOMPLETE, XLAL_EFAILED);
  const size_t n = tiling->dimensions;

  // If finished status, nothing more to be done!
  if (tiling->status == FLT_S_FINISHED) {
    return -1;
  }

  // If started status, but no tiled dimensions, we're finished!
  if (tiling->status == FLT_S_STARTED && tiling->tiled_dimensions == 0) {
    tiling->status = FLT_S_FINISHED;
    tiling->total_count = 1;   // Set total template count
    return -1;
  }

  // If initialised status, set and return starting point
  if (tiling->status == FLT_S_INITIALISED) {

    // Initialise parameter-space bound indices
    gsl_vector_uint_set_zero(tiling->bound_idx);

    // Set parameter-space bounds and starting point
    for (size_t i = 0; i < n; ++i) {

      // Get normalised bounds at bound index
      const size_t bound = gsl_vector_uint_get(tiling->bound_idx, i);
      double lower = 0, upper = 0;
      FLT_GetBounds(tiling, i, bound, tiling->phys_offset, &lower, &upper);

      // Set parameter-space bounds
      gsl_vector_set(tiling->lower, i, lower);
      gsl_vector_set(tiling->upper, i, upper);

      // Initialise current point
      const double point = lower - gsl_vector_get(tiling->padding, i);
      gsl_vector_set(tiling->point, i, point);

    }

    // Initialise count
    tiling->count = 1;

    // Tiling has been started.
    tiling->status = FLT_S_STARTED;

    // All dimensions of point have changed.
    return 0;

  }

  // Otherwise started status: loop until the next point is found
  size_t i = n, ir;
  while (true) {

    // If dimension index is now zero, we're done!
    if (i == 0) {
      tiling->status = FLT_S_FINISHED;
      tiling->total_count = tiling->count;   // Save total template count
      return -1;
    }

    // Decrement current dimension index
    --i;

    // Return point to lower bound in higher dimensions
    ir = i + 1;

    // Get current bound index
    size_t bound = gsl_vector_uint_get(tiling->bound_idx, i);

    // If dimension is tiled...
    if (tiling->bounds[i].tiled) {

      // Get increment vector
      gsl_vector_view increment = gsl_matrix_column(tiling->increment, i);

      // Increment current point along index
      gsl_vector_add(tiling->point, &increment.vector);

      // If point is not out of bounds, we have found a template point
      const double point = gsl_vector_get(tiling->point, i);
      const double upper = gsl_vector_get(tiling->upper, i);
      const double padding = gsl_vector_get(tiling->padding, i);
      if (point <= upper + padding) {
        break;
      }

    }

    // Increment bound index
    ++bound;

    if (bound < tiling->bounds[i].num_bounds) {

      // Set current bound index
      gsl_vector_uint_set(tiling->bound_idx, i, bound);

      // Return point to new lower bound in this dimension
      ir = i;

      // Found a template point
      break;

    } else {

      // If no more bounds, reset bound index in this and higher dimensions
      for (size_t j = i; j < n; ++j) {
        gsl_vector_uint_set(tiling->bound_idx, j, 0);
      }

    }

    // Move on to lower dimensions
    continue;

  }

  // Return point to lower bound in appropriate dimensions
  for (; ir < n; ++ir) {

    // Get normalised bounds at bound index
    const size_t bound = gsl_vector_uint_get(tiling->bound_idx, ir);
    double lower = 0, upper = 0;
    FLT_GetBounds(tiling, ir, bound, tiling->phys_offset, &lower, &upper);

    // Set parameter-space bounds
    gsl_vector_set(tiling->lower, ir, lower);
    gsl_vector_set(tiling->upper, ir, upper);

    // If dimension is tiled...
    if (tiling->bounds[ir].tiled) {

      // Get increment vector
      gsl_vector_view increment = gsl_matrix_column(tiling->increment, ir);

      // Calculate the distance from current point to the lower bound, in integer number of increments
      const double padding = gsl_vector_get(tiling->padding, ir);
      const double point = gsl_vector_get(tiling->point, ir);
      const double dist = ceil((lower - padding - point) / gsl_vector_get(&increment.vector, ir));

      // Move point back to lower bound
      gsl_blas_daxpy(dist, &increment.vector, tiling->point);

    } else {

      // Otherwise set point to lower bound
      gsl_vector_set(tiling->point, ir, lower);

    }

  }

  // Template was found, so increase count
  ++tiling->count;

  // Return lowest dimension where point has changed
  return i;

}

int XLALRestartFlatLatticeTiling(
  FlatLatticeTiling* tiling
  )
{

  // Check tiling
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling->status > FLT_S_INCOMPLETE, XLAL_EFAILED);

  // Restart tiling
  tiling->status = FLT_S_INITIALISED;
  tiling->count = 0;

  return XLAL_SUCCESS;

}

/* int XLALNearestFlatLatticePointToRandomPoints( */
/*   FlatLatticeTiling* tiling, */
/*   RandomParams* rng, */
/*   const size_t num_random_points, */
/*   gsl_matrix** random_points, */
/*   gsl_vector_ulong** nearest_indices, */
/*   gsl_vector** nearest_distances, */
/*   gsl_matrix** workspace */
/*   ) */
/* { */

/*   // Check tiling */
/*   XLAL_CHECK(tiling != NULL, XLAL_EFAULT); */
/*   XLAL_CHECK(tiling->status > FLT_S_INCOMPLETE, XLAL_EFAILED); */
/*   const size_t n = tiling->dimensions; */

/*   // Check input */
/*   XLAL_CHECK(rng != NULL, XLAL_EFAULT); */
/*   XLAL_CHECK(num_random_points > 0, XLAL_ESIZE); */
/*   XLAL_CHECK(random_points != NULL, XLAL_EFAULT); */
/*   XLAL_CHECK(nearest_indices != NULL, XLAL_EFAULT); */
/*   XLAL_CHECK(nearest_distances != NULL, XLAL_EFAULT); */
/*   XLAL_CHECK(workspace != NULL, XLAL_EFAULT); */

/*   // (Re)Allocate matrix of random points */
/*   if (*random_points != NULL && (*random_points)->size2 != num_random_points) { */
/*     gsl_matrix_free(*random_points); */
/*     *random_points = NULL; */
/*   } */
/*   if (*random_points == NULL) { */
/*     *random_points = gsl_matrix_alloc(n, num_random_points); */
/*     XLAL_CHECK(*random_points != NULL, XLAL_ENOMEM); */
/*   } */

/*   // (Re)Allocate vector of indices of nearest lattice point */
/*   if (*nearest_indices != NULL && (*nearest_indices)->size != num_random_points) { */
/*     gsl_vector_ulong_free(*nearest_indices); */
/*     *nearest_indices = NULL; */
/*   } */
/*   if (*nearest_indices == NULL) { */
/*     *nearest_indices = gsl_vector_ulong_alloc(num_random_points); */
/*     XLAL_CHECK(*nearest_indices != NULL, XLAL_ENOMEM); */
/*   } */

/*   // (Re)Allocate vector of distances to nearest lattice point */
/*   if (*nearest_distances != NULL && (*nearest_distances)->size != num_random_points) { */
/*     gsl_vector_free(*nearest_distances); */
/*     *nearest_distances = NULL; */
/*   } */
/*   if (*nearest_distances == NULL) { */
/*     *nearest_distances = gsl_vector_alloc(num_random_points); */
/*     XLAL_CHECK(*nearest_distances != NULL, XLAL_ENOMEM); */
/*   } */

/*   // (Re)Allocate workspace matrix for computing distances */
/*   if (*workspace != NULL && (*workspace)->size2 != num_random_points) { */
/*     gsl_matrix_free(*workspace); */
/*     *workspace = NULL; */
/*   } */
/*   if (*workspace == NULL) { */
/*     *workspace = gsl_matrix_alloc(3*n - 1, num_random_points); */
/*     XLAL_CHECK(*workspace != NULL, XLAL_ENOMEM); */
/*   } */

/*   // Create temporary bound index and physical bound vectors */
/*   gsl_vector_uint* bound_idx = gsl_vector_uint_alloc(n); */
/*   XLAL_CHECK(bound_idx != NULL, XLAL_ENOMEM); */
/*   gsl_vector* phys_lower = gsl_vector_alloc(tiling->max_num_bounds); */
/*   XLAL_CHECK(phys_lower != NULL, XLAL_ENOMEM); */
/*   gsl_vector* phys_width = gsl_vector_alloc(tiling->max_num_bounds); */
/*   XLAL_CHECK(phys_width != NULL, XLAL_ENOMEM); */

/*   // Create random points in flat lattice tiling parameter space */
/*   for (size_t k = 0; k < num_random_points; ++k) { */
/*     gsl_vector_view point = gsl_matrix_column(*random_points, k); */
/*     for (size_t i = 0; i < n; ++i) { */

/*       // Get physical bounds and padding */
/*       FLT_GetBounds(tiling, i, bound_idx, &point.vector, */
/*                         phys_lower, phys_width, NULL, NULL); */
/*       gsl_vector_sub(phys_width, phys_lower); */

/*       // Get total bounds width */
/*       double phys_total_width = 0; */
/*       size_t max_bounds = 0; */
/*       while (max_bounds < tiling->bounds[i].num_bounds) { */
/*         const double lower = gsl_vector_get(phys_lower, max_bounds); */
/*         const double width = gsl_vector_get(phys_width, max_bounds); */
/*         if (gsl_isnan(lower) && gsl_isnan(width)) { */
/*           break; */
/*         } */
/*         phys_total_width += width; */
/*         ++max_bounds; */
/*       } */

/*       // Generate random number */
/*       const double u = XLALUniformDeviate(rng); */

/*       double p; */
/*       size_t bound = 0; */
/*       if (tiling->bounds[i].tiled) { */

/*         // Generate random point within total bounds widths */
/*         p = u * phys_total_width; */

/*         // Convert point to be within parameter-space bounds */
/*         while (bound + 1 < max_bounds) { */
/*           const double width = gsl_vector_get(phys_width, bound); */
/*           if (p <= width) { */
/*             break; */
/*           } */
/*           p -= width; */
/*           ++bound; */
/*         } */
/*         p += gsl_vector_get(phys_lower, bound); */

/*       } else { */

/*         // Generate random bound index */
/*         bound = (size_t)floor(u * max_bounds); */

/*         // Get point from random bound */
/*         p = gsl_vector_get(phys_lower, bound); */

/*       } */

/*       // Set parameter-space point and bound index */
/*       gsl_vector_set(&point.vector, i, p); */
/*       gsl_vector_uint_set(bound_idx, i, bound); */

/*     } */

/*   } */

/*   // Create temporary matrices in workspace */
/*   gsl_matrix_view point_diffs = gsl_matrix_submatrix(*workspace, 0, 0, n, num_random_points); */
/*   gsl_matrix_view off_diag_terms = gsl_matrix_submatrix(*workspace, n, 0, n - 1, num_random_points); */
/*   gsl_matrix_view distances = gsl_matrix_submatrix(*workspace, 2*n - 1, 0, n, num_random_points); */

/*   // Initialise minimum distance vector */
/*   gsl_vector_set_all(*nearest_distances, GSL_POSINF); */

/*   // Iterate over all flat lattice points */
/*   XLALRestartFlatLatticeTiling(tiling); */
/*   while (true) { */

/*     // Advance to the next lattice point */
/*     const int ich = XLALNextFlatLatticePoint(tiling); */
/*     if (ich < 0) { */
/*       break; */
/*     } */
/*     const gsl_vector* lattice_point = XLALGetFlatLatticePoint(tiling); */
/*     const unsigned long nearest_index = tiling->count - 1; */

/*     // For dimensions where flat lattice point has changed (given by ich), */
/*     // copy random points to workspace, subtract flat lattice point from each, */
/*     // and normalise by physical scaling */
/*     for (size_t i = (size_t)ich; i < n; ++i) { */
/*       const double phys_scale = gsl_vector_get(tiling->phys_scale, i); */
/*       gsl_vector_view point_diffs_i = gsl_matrix_row(&point_diffs.matrix, i); */
/*       gsl_vector_view random_points_i = gsl_matrix_row(*random_points, i); */
/*       gsl_vector_memcpy(&point_diffs_i.vector, &random_points_i.vector); */
/*       gsl_vector_add_constant(&point_diffs_i.vector, -gsl_vector_get(lattice_point, i)); */
/*       gsl_vector_scale(&point_diffs_i.vector, 1.0/phys_scale); */
/*     } */

/*     // For dimensions where flat lattice point has changed (given by ich), */
/*     // re-compute the off-diagonal terms of the metric distance, which */
/*     // are multiplied by the ith coordinate difference */
/*     for (size_t i = (size_t)ich; i < n - 1; ++i) { */
/*       gsl_vector_view off_diag_terms_i = gsl_matrix_row(&off_diag_terms.matrix, i); */
/*       gsl_vector_set_zero(&off_diag_terms_i.vector); */
/*       for (size_t j = 0; j <= i; ++j) { */
/*         const double metric_off_diag = gsl_matrix_get(tiling->metric, i + 1, j); */
/*         gsl_vector_view point_diffs_j = gsl_matrix_row(&point_diffs.matrix, j); */
/*         gsl_blas_daxpy(2.0 * metric_off_diag, &point_diffs_j.vector, &off_diag_terms_i.vector); */
/*       } */
/*     } */

/*     // For dimensions where flat lattice point has changed (given by ich), */
/*     // re-compute terms in the distances from random points to the flat lattice */
/*     // point which involve the ith coordinate difference, and cumulatively sum */
/*     // together to get the distance in the last row */
/*     for (size_t i = (size_t)ich; i < n; ++i) { */

/*       gsl_vector_view point_diffs_i = gsl_matrix_row(&point_diffs.matrix, i); */
/*       gsl_vector_view distances_i = gsl_matrix_row(&distances.matrix, i); */

/*       // Compute the diagonal term of the metric distance, */
/*       // which are multiplied by the ith coordinate difference */
/*       const double metric_diag = gsl_matrix_get(tiling->metric, i, i); */
/*       gsl_vector_memcpy(&distances_i.vector, &point_diffs_i.vector); */
/*       gsl_vector_scale(&distances_i.vector, metric_diag); */

/*       // Add the pre-computed off-diagomal terms of the metric distance, */
/*       // which are multiplied by the ith coordinate difference */
/*       if (i > 0) { */
/*         gsl_vector_view off_diag_terms_iprev = gsl_matrix_row(&off_diag_terms.matrix, i - 1); */
/*         gsl_vector_add(&distances_i.vector, &off_diag_terms_iprev.vector); */
/*       } */

/*       // Multiply by the ith coordinate difference */
/*       gsl_vector_mul(&distances_i.vector, &point_diffs_i.vector); */

/*       // Add the distance computed for the lower dimensions thus far */
/*       if (i > 0) { */
/*         gsl_vector_view distances_iprev = gsl_matrix_row(&distances.matrix, i - 1); */
/*         gsl_vector_add(&distances_i.vector, &distances_iprev.vector); */
/*       } */

/*     } */

/*     // For each random point, if the distance to the flat lattice point is */
/*     // the smallest so far, record the flat lattice point, distance, and index */
/*     gsl_vector_view distance = gsl_matrix_row(&distances.matrix, n - 1); */
/*     for (size_t k = 0; k < num_random_points; ++k) { */
/*       const double distance_k = gsl_vector_get(&distance.vector, k); */
/*       if (distance_k < gsl_vector_get(*nearest_distances, k)) { */
/*         gsl_vector_ulong_set(*nearest_indices, k, nearest_index); */
/*         gsl_vector_set(*nearest_distances, k, distance_k); */
/*       } */
/*     } */

/*   } */
/*   XLAL_CHECK(xlalErrno == 0, XLAL_EFAILED, "XLALNextFlatLatticePoint() failed"); */

/*   // Cleanup */
/*   gsl_vector_uint_free(bound_idx); */
/*   gsl_vector_free(phys_lower); */
/*   gsl_vector_free(phys_width); */

/*   return XLAL_SUCCESS; */

/* } */
