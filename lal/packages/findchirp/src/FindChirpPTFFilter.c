/*
*  Copyright (C) 2007 Diego Fazi, Duncan Brown
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

/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpPTFFilter.c
 *
 * Author: Brown, D. A. and Fazi, D.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpPTFFilterCV">
Author: Brown, D. A. and Fazi, D.
$Id$
</lalVerbatim>

<lalLaTeX>
\input{FindChirpPTFFilterCDoc}

\vfill{\footnotesize\input{FindChirpPTFFilterCV}}
#endif

#include <config.h>
#include <stdlib.h>
#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/FindChirp.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirpPTF.h>
#include <lal/MatrixUtils.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_linalg.h>

#define FLOAT_RADIX         2.0
#define FLOAT_RADIX_SQ      (FLOAT_RADIX * FLOAT_RADIX)
#define GSL_FRANCIS_COEFF1  (0.75)
#define GSL_FRANCIS_COEFF2  (-0.4375)

NRCSID (FINDCHIRPPTFFILTERC, "$Id$");

/* <lalVerbatim file="FindChirpPTFFilterCP"> */

double rint(double x);

/* * * * * * * * * * * * * * * * * * * * * 
 *
 *  Beginning of pseudo-gsl code
 *
 * * * * * * * * * * * * * * * * * * * * *
 */

typedef struct {
  size_t size;           /* matrix size */
  size_t max_iterations; /* max iterations since last eigenvalue found */
  size_t n_iter;         /* number of iterations since last eigenvalue found */
  size_t n_evals;        /* number of eigenvalues found so far */

  gsl_vector *hv2;       /* temporary 2-by-1 _householder vector */
  gsl_vector *hv3;       /* temporary 3-by-1 householder vector */

  int compute_t;         /* compute Schur form T = Z^t A Z */

  gsl_matrix *H;         /* pointer to Hessenberg matrix */

  gsl_matrix *Z;         /* pointer to Schur vector matrix */
} _gsl_eigen_francis_workspace;


static        _gsl_eigen_francis_workspace * _gsl_eigen_francis_alloc (void);
static int    _gsl_eigen_francis_Z (gsl_matrix * H, gsl_vector_complex * eval,
                            gsl_matrix * Z, _gsl_eigen_francis_workspace * w);
static void   _gsl_eigen_francis_free (_gsl_eigen_francis_workspace * w);
static void   _gsl_eigen_francis_T (const int compute_t,
                             _gsl_eigen_francis_workspace * w);
static int    _gsl_eigen_francis (gsl_matrix * H, gsl_vector_complex * eval,
                                  _gsl_eigen_francis_workspace * w);
static int    _gsl_linalg_balance_matrix( gsl_matrix * A, gsl_vector * D );
static int    _gsl_linalg_hessenberg(gsl_matrix *A, gsl_vector *tau);
static int    _gsl_linalg_hessenberg_unpack(gsl_matrix * H, gsl_vector * tau,
                                  gsl_matrix * U);
static int    _gsl_linalg_hessenberg_unpack_accum(gsl_matrix * H, gsl_vector * tau,
                                        gsl_matrix * U);
static void   _gsl_schur_standardize(gsl_matrix *T, size_t row, gsl_complex *eval1,
                                  gsl_complex *eval2, int update_t, gsl_matrix *Z);


typedef struct {
  size_t size;                 /* size of matrices */
  gsl_vector *diag;            /* diagonal matrix elements from balancing */
  gsl_vector *tau;             /* Householder coefficients */
  gsl_matrix *Z;               /* pointer to Z matrix */
  int do_balance;              /* perform balancing transformation? */
  size_t n_evals;              /* number of eigenvalues found */

  _gsl_eigen_francis_workspace *francis_workspace_p;
} _gsl_eigen_nonsymm_workspace;

static        _gsl_eigen_nonsymm_workspace * _gsl_eigen_nonsymm_alloc (const size_t n);
static void   _gsl_eigen_nonsymm_free (_gsl_eigen_nonsymm_workspace * w);
static void   _gsl_eigen_nonsymm_params (const int compute_t, const int balance,
                                          _gsl_eigen_nonsymm_workspace *w);
static int    _gsl_eigen_nonsymm (gsl_matrix * A, gsl_vector_complex * eval,
                                  _gsl_eigen_nonsymm_workspace * w);

/* Balance a general matrix by scaling the rows and columns, so the
 * new row and column norms are the same order of magnitude.
 *
 * B =  D^-1 A D
 *
 * where D is a diagonal matrix
 * 
 * This is necessary for the unsymmetric eigenvalue problem since the
 * calculation can become numerically unstable for unbalanced
 * matrices.  
 *
 * See Golub & Van Loan, "Matrix Computations" (3rd ed), Section 7.5.7
 * and Wilkinson & Reinsch, "Handbook for Automatic Computation", II/11 p320.
 */



static int
_gsl_linalg_balance_matrix(gsl_matrix * A, gsl_vector * D)
{
  const size_t N = A->size1;

  if (N != D->size)
    {
      GSL_ERROR ("vector must match matrix size", GSL_EBADLEN);
    }
  else
    {
      double row_norm,
             col_norm;
      int not_converged;
      gsl_vector_view v;

      /* initialize D to the identity matrix */
      gsl_vector_set_all(D, 1.0);

      not_converged = 1;

      while (not_converged)
        {
          size_t i, j;
          double g, f, s;

          not_converged = 0;

          for (i = 0; i < N; ++i)
            {
              row_norm = 0.0;
              col_norm = 0.0;

              for (j = 0; j < N; ++j)
                {
                  if (j != i)
                    {
                      col_norm += fabs(gsl_matrix_get(A, j, i));
                      row_norm += fabs(gsl_matrix_get(A, i, j));
                    }
                }

              if ((col_norm == 0.0) || (row_norm == 0.0))
                {
                  continue;
                }

              g = row_norm / FLOAT_RADIX;
              f = 1.0;
              s = col_norm + row_norm;

              /*
               * find the integer power of the machine radix which
               * comes closest to balancing the matrix
               */
              while (col_norm < g)
                {
                  f *= FLOAT_RADIX;
                  col_norm *= FLOAT_RADIX_SQ;
                }

              g = row_norm * FLOAT_RADIX;

              while (col_norm > g)
                {
                  f /= FLOAT_RADIX;
                  col_norm /= FLOAT_RADIX_SQ;
                }

              if ((row_norm + col_norm) < 0.95 * s * f)
                {
                  not_converged = 1;

                  g = 1.0 / f;

                  /*
                   * apply similarity transformation D, where
                   * D_{ij} = f_i * delta_{ij}
                   */

                  /* multiply by D^{-1} on the left */
                  v = gsl_matrix_row(A, i);
                  gsl_blas_dscal(g, &v.vector);

                  /* multiply by D on the right */
                  v = gsl_matrix_column(A, i);
                  gsl_blas_dscal(f, &v.vector);

                  /* keep track of transformation */
                  gsl_vector_set(D, i, gsl_vector_get(D, i) * f);
                }
            }
        }

      return GSL_SUCCESS;
    }
} /* gsl_linalg_balance_matrix() */

/*
gsl_linalg_balance_accum()
  Accumulate a balancing transformation into a matrix.
This is used during the computation of Schur vectors since the
Schur vectors computed are the vectors for the balanced matrix.
We must at some point accumulate the balancing transformation into
the Schur vector matrix to get the vectors for the original matrix.

A -> D A

where D is the diagonal matrix

Inputs: A - matrix to transform
        D - vector containing diagonal elements of D
*/

static int
_gsl_linalg_balance_accum(gsl_matrix *A, gsl_vector *D)
{
  const size_t N = A->size1;

  if (N != D->size)
    {
      GSL_ERROR ("vector must match matrix size", GSL_EBADLEN);
    }
  else
    {
      size_t i;
      double s;
      gsl_vector_view r;

      for (i = 0; i < N; ++i)
        {
          s = gsl_vector_get(D, i);
          r = gsl_matrix_row(A, i);

          gsl_blas_dscal(s, &r.vector);
        }

      return GSL_SUCCESS;
    }
} /* gsl_linalg_balance_accum() */


/*
gsl_linalg_hessenberg()
  Compute the Householder reduction to Hessenberg form of a
square N-by-N matrix A.

H = U^t A U

See Golub & Van Loan, "Matrix Computations" (3rd ed), algorithm
7.4.2

Inputs: A   - matrix to reduce
        tau - where to store scalar factors in Householder
              matrices; this vector must be of length N,
              where N is the order of A

Return: GSL_SUCCESS unless error occurs

Notes: on output, the upper triangular portion of A (including
the diagaonal and subdiagonal) contains the Hessenberg matrix.
The lower triangular portion (below the subdiagonal) contains
the Householder vectors which can be used to construct
the similarity transform matrix U.

The matrix U is

U = U(1) U(2) ... U(n - 2)

where

U(i) = I - tau(i) * v(i) * v(i)^t

and the vector v(i) is stored in column i of the matrix A
underneath the subdiagonal. So the first element of v(i)
is stored in row i + 2, column i, the second element at
row i + 3, column i, and so on.

Also note that for the purposes of computing U(i),
v(1:i) = 0, v(i + 1) = 1, and v(i+2:n) is what is stored in
column i of A beneath the subdiagonal.
*/

static int
_gsl_linalg_hessenberg(gsl_matrix *A, gsl_vector *tau)
{
  const size_t N = A->size1;

  if (N != A->size2)
    {
      GSL_ERROR ("Hessenberg reduction requires square matrix",
                 GSL_ENOTSQR);
    }
  else if (N != tau->size)
    {
      GSL_ERROR ("tau vector must match matrix size", GSL_EBADLEN);
    }
  else if (N < 3)
    {
      /* nothing to do */
      return GSL_SUCCESS;
    }
  else
    {
      size_t i;           /* looping */
      gsl_vector_view c,  /* matrix column */
                      hv; /* householder vector */
      gsl_matrix_view m;
      double tau_i;       /* beta in algorithm 7.4.2 */

      for (i = 0; i < N - 2; ++i)
        {
          /*
           * make a copy of A(i + 1:n, i) and store it in the section
           * of 'tau' that we haven't stored coefficients in yet
           */

          c = gsl_matrix_column(A, i);
          c = gsl_vector_subvector(&c.vector, i + 1, N - (i + 1));

          hv = gsl_vector_subvector(tau, i + 1, N - (i + 1));
          gsl_vector_memcpy(&hv.vector, &c.vector);

          /* compute householder transformation of A(i+1:n,i) */
          tau_i = gsl_linalg_householder_transform(&hv.vector);

          /* apply left householder matrix (I - tau_i v v') to A */
          m = gsl_matrix_submatrix(A, i + 1, i, N - (i + 1), N - i);
          gsl_linalg_householder_hm(tau_i, &hv.vector, &m.matrix);

          /* apply right householder matrix (I - tau_i v v') to A */
          m = gsl_matrix_submatrix(A, 0, i + 1, N, N - (i + 1));
          gsl_linalg_householder_mh(tau_i, &hv.vector, &m.matrix);

          /* save Householder coefficient */
          gsl_vector_set(tau, i, tau_i);

          /*
           * store Householder vector below the subdiagonal in column
           * i of the matrix. hv(1) does not need to be stored since
           * it is always 1.
           */
          c = gsl_vector_subvector(&c.vector, 1, c.vector.size - 1);
          hv = gsl_vector_subvector(&hv.vector, 1, hv.vector.size - 1);
          gsl_vector_memcpy(&c.vector, &hv.vector);
        }

      return GSL_SUCCESS;
    }
} /* gsl_linalg_hessenberg() */

/*
gsl_linalg_hessenberg_unpack()
  Construct the matrix U which transforms a matrix A into
its upper Hessenberg form:

H = U^t A U

by unpacking the information stored in H from gsl_linalg_hessenberg().

U is a product of Householder matrices:

U = U(1) U(2) ... U(n - 2)

where

U(i) = I - tau(i) * v(i) * v(i)^t

The v(i) are stored in the lower triangular part of H by
gsl_linalg_hessenberg(). The tau(i) are stored in the vector tau.

Inputs: H       - Hessenberg matrix computed from
                  gsl_linalg_hessenberg()
        tau     - tau vector computed from gsl_linalg_hessenberg()
        U       - (output) where to store similarity matrix

Return: success or error
*/

static int
_gsl_linalg_hessenberg_unpack(gsl_matrix * H, gsl_vector * tau,
                              gsl_matrix * U)
{
  int s;

  gsl_matrix_set_identity(U);

  s = _gsl_linalg_hessenberg_unpack_accum(H, tau, U);

  return s;
} /* gsl_linalg_hessenberg_unpack() */

/*
gsl_linalg_hessenberg_unpack_accum()
  This routine is the same as gsl_linalg_hessenberg_unpack(), except
instead of storing the similarity matrix in U, it accumulates it,
so that

U -> U * [ U(1) U(2) ... U(n - 2) ]

instead of:

U -> U(1) U(2) ... U(n - 2)

Inputs: H       - Hessenberg matrix computed from
                  gsl_linalg_hessenberg()
        tau     - tau vector computed from gsl_linalg_hessenberg()
        V       - (input/output) where to accumulate similarity matrix

Return: success or error

Notes: 1) On input, V needs to be initialized. The Householder matrices
          are accumulated into V, so on output,

            V_out = V_in * U(1) * U(2) * ... * U(n - 2)

          so if you just want the product of the Householder matrices,
          initialize V to the identity matrix before calling this
          function.

       2) V does not have to be square, but must have the same
          number of columns as the order of H
*/

static int
_gsl_linalg_hessenberg_unpack_accum(gsl_matrix * H, gsl_vector * tau,
                                    gsl_matrix * V)
{
  const size_t N = H->size1;

  if (N != H->size2)
    {
      GSL_ERROR ("Hessenberg reduction requires square matrix",
                 GSL_ENOTSQR);
    }
  else if (N != tau->size)
    {
      GSL_ERROR ("tau vector must match matrix size", GSL_EBADLEN);
    }
  else if (N != V->size2)
    {
      GSL_ERROR ("V matrix has wrong dimension", GSL_EBADLEN);
    }
  else
    {
      size_t j;           /* looping */
      double tau_j;       /* householder coefficient */
      gsl_vector_view c,  /* matrix column */
                      hv; /* householder vector */
      gsl_matrix_view m;

      if (N < 3)
        {
          /* nothing to do */
          return GSL_SUCCESS;
        }

      for (j = 0; j < (N - 2); ++j)
        {
          c = gsl_matrix_column(H, j);

          tau_j = gsl_vector_get(tau, j);

          /*
           * get a view to the householder vector in column j, but
           * make sure hv(2) starts at the element below the
           * subdiagonal, since hv(1) was never stored and is always
           * 1
           */
          hv = gsl_vector_subvector(&c.vector, j + 1, N - (j + 1));

          /*
           * Only operate on part of the matrix since the first
           * j + 1 entries of the real householder vector are 0
           *
           * V -> V * U(j)
           *
           * Note here that V->size1 is not necessarily equal to N
           */
          m = gsl_matrix_submatrix(V, 0, j + 1, V->size1, N - (j + 1));

          /* apply right Householder matrix to V */
          gsl_linalg_householder_mh(tau_j, &hv.vector, &m.matrix);
        }

      return GSL_SUCCESS;
    }
} /* gsl_linalg_hessenberg_unpack_accum() */

/*
 * This module computes the eigenvalues of a real upper hessenberg
 * matrix, using the classical double shift Francis QR algorithm.
 * It will also optionally compute the full Schur form and matrix of
 * Schur vectors.
 *
 * See Golub & Van Loan, "Matrix Computations" (3rd ed), algorithm 7.5.2
 */

/* exceptional shift coefficients - these values are from LAPACK DLAHQR */


static void _francis_schur_decomp(gsl_matrix * H,
                                         gsl_vector_complex * eval,
                                         _gsl_eigen_francis_workspace * w);
static size_t _francis_search_subdiag_small_elements(gsl_matrix * A);
static int _francis_qrstep(gsl_matrix * H,
                                     _gsl_eigen_francis_workspace * w);
static void _francis_schur_standardize(gsl_matrix *A,
                                              gsl_complex *eval1,
                                              gsl_complex *eval2,
                                              _gsl_eigen_francis_workspace *w);
static void _francis_get_submatrix(gsl_matrix *A, gsl_matrix *B,
                                             size_t *top); 


/*    
 *    francis_schur_decomp()
 *      Compute the Schur decomposition of the matrix H
 *
 *      Inputs: H     - hessenberg matrix
 *              eval  - where to store eigenvalues
 *                      w     - workspace
 *                            
 *                            Return: none
 *                            */
      
static  void
_francis_schur_decomp(gsl_matrix * H, gsl_vector_complex * eval,
                       _gsl_eigen_francis_workspace * w)
{     
  gsl_matrix_view m;   /* active matrix we are working on */
  size_t N;            /* size of matrix */ 
  size_t q;            /* index of small subdiagonal element */      
  gsl_complex lambda1, /* eigenvalues */
              lambda2;

  N = H->size1;

  if (N == 1)
  {
    GSL_SET_COMPLEX(&lambda1, gsl_matrix_get(H, 0, 0), 0.0);
    gsl_vector_complex_set(eval, w->n_evals, lambda1);
    w->n_evals += 1;
    w->n_iter = 0;
    return;
  }
  else if (N == 2)
  {
    _francis_schur_standardize(H, &lambda1, &lambda2, w);
    gsl_vector_complex_set(eval, w->n_evals, lambda1);
    gsl_vector_complex_set(eval, w->n_evals + 1, lambda2);
    w->n_evals += 2;
    w->n_iter = 0;
    return;
  }
  m = gsl_matrix_submatrix(H, 0, 0, N, N);

  while ((N > 2) && ((w->n_iter)++ < w->max_iterations))
  {
    q = _francis_search_subdiag_small_elements(&m.matrix);

    if (q == 0)
    {
      /*
       * no small subdiagonal element found - perform a QR
       * sweep on the active reduced hessenberg matrix
       */
      _francis_qrstep(&m.matrix, w);
      continue;
    }

    /*
     * a small subdiagonal element was found - one or two eigenvalues
     * have converged or the matrix has split into two smaller matrices
     */

    if (q == (N - 1))
    {
    /*
     * the last subdiagonal element of the matrix is 0 -
     * m_{NN} is a real eigenvalue
     */
      GSL_SET_COMPLEX(&lambda1,
          gsl_matrix_get(&m.matrix, q, q), 0.0);
      gsl_vector_complex_set(eval, w->n_evals, lambda1);
      w->n_evals += 1;
      w->n_iter = 0;

      --N;
      m = gsl_matrix_submatrix(&m.matrix, 0, 0, N, N);
    }
    else if (q == (N - 2))
    {
      gsl_matrix_view v;

      /*
       * The bottom right 2-by-2 block of m is an eigenvalue system
       */

      v = gsl_matrix_submatrix(&m.matrix, q, q, 2, 2);
      _francis_schur_standardize(&v.matrix, &lambda1, &lambda2, w);
      gsl_vector_complex_set(eval, w->n_evals, lambda1);
      gsl_vector_complex_set(eval, w->n_evals + 1, lambda2);
      w->n_evals += 2;
      w->n_iter = 0;

      N -= 2;
      m = gsl_matrix_submatrix(&m.matrix, 0, 0, N, N);
    }
    else if (q == 1)
    {
      /* the first matrix element is an eigenvalue */
      GSL_SET_COMPLEX(&lambda1,
          gsl_matrix_get(&m.matrix, 0, 0), 0.0);
      gsl_vector_complex_set(eval, w->n_evals, lambda1);
      w->n_evals += 1;
      w->n_iter = 0;

      --N;
      m = gsl_matrix_submatrix(&m.matrix, 1, 1, N, N);
    }
    else if (q == 2)
    {
      gsl_matrix_view v;

      /* the upper left 2-by-2 block is an eigenvalue system */

      v = gsl_matrix_submatrix(&m.matrix, 0, 0, 2, 2);
      _francis_schur_standardize(&v.matrix, &lambda1, &lambda2, w);

      gsl_vector_complex_set(eval, w->n_evals, lambda1);
      gsl_vector_complex_set(eval, w->n_evals + 1, lambda2);
      w->n_evals += 2;
      w->n_iter = 0;

      N -= 2;
      m = gsl_matrix_submatrix(&m.matrix, 2, 2, N, N);
    }
    else
    {
      gsl_matrix_view v;
      
      /*
       * There is a zero element on the subdiagonal somewhere 
       * in the middle of the matrix - we can now operate
       * separately on the two submatrices split by this
       * element. q is the row index of the zero element.
       */
      
      /* operate on lower right (N - q)-by-(N - q) block first */
      v = gsl_matrix_submatrix(&m.matrix, q, q, N - q, N - q);
      _francis_schur_decomp(&v.matrix, eval, w);

      /* operate on upper left q-by-q block */
      v = gsl_matrix_submatrix(&m.matrix, 0, 0, q, q);
      _francis_schur_decomp(&v.matrix, eval, w);

      N = 0;
    }
  }
  if (N == 1)
  {
    GSL_SET_COMPLEX(&lambda1, gsl_matrix_get(&m.matrix, 0, 0), 0.0);
    gsl_vector_complex_set(eval, w->n_evals, lambda1);
    w->n_evals += 1;
    w->n_iter = 0;
  }
  else if (N == 2)
  {
    _francis_schur_standardize(&m.matrix, &lambda1, &lambda2, w);
    gsl_vector_complex_set(eval, w->n_evals, lambda1);
    gsl_vector_complex_set(eval, w->n_evals + 1, lambda2);
    w->n_evals += 2;
    w->n_iter = 0;
                                              }
} /* francis_schur_decomp() */

/*
francis_qrstep()
  Perform a Francis QR step.

See Golub & Van Loan, "Matrix Computations" (3rd ed),
algorithm 7.5.1

Inputs: H - upper Hessenberg matrix
        w - workspace

Notes: The matrix H must be "reduced", ie: have no tiny subdiagonal
       elements. When computing the first householder reflection,
       we divide by H_{21} so it is necessary that this element
       is not zero. When a subdiagonal element becomes negligible,
       the calling function should call this routine with the
       submatrices split by that element, so that we don't divide
       by zeros.
*/


static  int
_francis_qrstep(gsl_matrix * H, _gsl_eigen_francis_workspace * w)
{
  const size_t N = H->size1;
  double x, y, z;  /* householder vector elements */
  double scale;    /* scale factor to avoid overflow */
  size_t i;        /* looping */
  gsl_matrix_view m;
  double tau_i;    /* householder coefficient */
  size_t q, r;
  size_t top;      /* location of H in original matrix */
  double s,
         disc;
  double h_nn,     /* H(n,n) */
         h_nm1nm1, /* H(n-1,n-1) */
         h_cross,  /* H(n,n-1) * H(n-1,n) */
         h_tmp1,
         h_tmp2;

  if ((w->n_iter == 10) || (w->n_iter == 20))
    {
      /*
       * exceptional shifts: we have gone 10 or 20 iterations
       * without finding a new eigenvalue, try a new choice of shifts.
       * See Numerical Recipes in C, eq 11.6.27 and LAPACK routine
       * DLAHQR
       */
      s = fabs(gsl_matrix_get(H, N - 1, N - 2)) +
          fabs(gsl_matrix_get(H, N - 2, N - 3));
      h_nn = gsl_matrix_get(H, N - 1, N - 1) + GSL_FRANCIS_COEFF1 * s;
      h_nm1nm1 = h_nn;
      h_cross = GSL_FRANCIS_COEFF2 * s * s;
    }
  else
    {
      /*
       * normal shifts - compute Rayleigh quotient and use
       * Wilkinson shift if possible
       */

      h_nn = gsl_matrix_get(H, N - 1, N - 1);
      h_nm1nm1 = gsl_matrix_get(H, N - 2, N - 2);
      h_cross = gsl_matrix_get(H, N - 1, N - 2) *
                gsl_matrix_get(H, N - 2, N - 1);

      disc = 0.5 * (h_nm1nm1 - h_nn);
      disc = disc * disc + h_cross;
      if (disc > 0.0)
        {
          double ave;

          /* real roots - use Wilkinson's shift twice */
          disc = sqrt(disc);
          ave = 0.5 * (h_nm1nm1 + h_nn);
          if (fabs(h_nm1nm1) - fabs(h_nn) > 0.0)
            {
              h_nm1nm1 = h_nm1nm1 * h_nn - h_cross;
              h_nn = h_nm1nm1 / (disc * GSL_SIGN(ave) + ave);
            }
          else
            {
              h_nn = disc * GSL_SIGN(ave) + ave;
            }

          h_nm1nm1 = h_nn;
          h_cross = 0.0;
        }
    }

  h_tmp1 = h_nm1nm1 - gsl_matrix_get(H, 0, 0);
  h_tmp2 = h_nn - gsl_matrix_get(H, 0, 0);

  /*
   * These formulas are equivalent to those in Golub & Van Loan
   * for the normal shift case - the terms have been rearranged
   * to reduce possible roundoff error when subdiagonal elements
   * are small
   */

  x = (h_tmp1*h_tmp2 - h_cross) / gsl_matrix_get(H, 1, 0) +
      gsl_matrix_get(H, 0, 1);
  y = gsl_matrix_get(H, 1, 1) - gsl_matrix_get(H, 0, 0) - h_tmp1 - h_tmp2;
  z = gsl_matrix_get(H, 2, 1);

  scale = fabs(x) + fabs(y) + fabs(z);
  if (scale != 0.0)
    {
      /* scale to prevent overflow or underflow */
      x /= scale;
      y /= scale;
      z /= scale;
    }

  if (w->Z || w->compute_t)
    {
      /*
       * get absolute indices of this (sub)matrix relative to the
       * original Hessenberg matrix
       */
      _francis_get_submatrix(w->H, H, &top);
    }

  for (i = 0; i < N - 2; ++i)
    {
      gsl_vector_set(w->hv3, 0, x);
      gsl_vector_set(w->hv3, 1, y);
      gsl_vector_set(w->hv3, 2, z);
      tau_i = gsl_linalg_householder_transform(w->hv3);

      if (tau_i != 0.0)
        {
          /* q = max(1, i - 1) */
          q = (1 > ((int)i - 1)) ? 0 : (i - 1);

          /* r = min(i + 3, N - 1) */
          r = ((i + 3) < (N - 1)) ? (i + 3) : (N - 1);

          if (w->compute_t)
            {
              /*
               * We are computing the Schur form T, so we
               * need to transform the whole matrix H
               *
               * H -> P_k^t H P_k
               *
               * where P_k is the current Householder matrix
               */

              /* apply left householder matrix (I - tau_i v v') to H */
              m = gsl_matrix_submatrix(w->H,
                                       top + i,
                                       top + q,
                                       3,
                                       w->size - top - q);
              gsl_linalg_householder_hm(tau_i, w->hv3, &m.matrix);

              /* apply right householder matrix (I - tau_i v v') to H */
              m = gsl_matrix_submatrix(w->H,
                                       0,
                                       top + i,
                                       top + r + 1,
                                       3);
              gsl_linalg_householder_mh(tau_i, w->hv3, &m.matrix);
            }
          else
            {
              /*
               * We are not computing the Schur form T, so we
               * only need to transform the active block
               */

              /* apply left householder matrix (I - tau_i v v') to H */
              m = gsl_matrix_submatrix(H, i, q, 3, N - q);
              gsl_linalg_householder_hm(tau_i, w->hv3, &m.matrix);

              /* apply right householder matrix (I - tau_i v v') to H */
              m = gsl_matrix_submatrix(H, 0, i, r + 1, 3);
              gsl_linalg_householder_mh(tau_i, w->hv3, &m.matrix);
            }

          if (w->Z)
            {
              /* accumulate the similarity transformation into Z */
              m = gsl_matrix_submatrix(w->Z, 0, top + i, w->size, 3);
              gsl_linalg_householder_mh(tau_i, w->hv3, &m.matrix);
            }
        } /* if (tau_i != 0.0) */

      x = gsl_matrix_get(H, i + 1, i);
      y = gsl_matrix_get(H, i + 2, i);
      if (i < (N - 3))
        {
          z = gsl_matrix_get(H, i + 3, i);
        }

      scale = fabs(x) + fabs(y) + fabs(z);
      if (scale != 0.0)
        {
          /* scale to prevent overflow or underflow */
          x /= scale;
          y /= scale;
          z /= scale;
        }
    } /* for (i = 0; i < N - 2; ++i) */

  gsl_vector_set(w->hv2, 0, x);
  gsl_vector_set(w->hv2, 1, y);
  tau_i = gsl_linalg_householder_transform(w->hv2);

  if (w->compute_t)
    {
      m = gsl_matrix_submatrix(w->H,
                               top + N - 2,
                               top + N - 3,
                               2,
                               w->size - top - N + 3);
      gsl_linalg_householder_hm(tau_i, w->hv2, &m.matrix);

      m = gsl_matrix_submatrix(w->H,
                               0,
                               top + N - 2,
                               top + N,
                               2);
      gsl_linalg_householder_mh(tau_i, w->hv2, &m.matrix);
    }
  else
    {
      m = gsl_matrix_submatrix(H, N - 2, N - 3, 2, 3);
      gsl_linalg_householder_hm(tau_i, w->hv2, &m.matrix);

      m = gsl_matrix_submatrix(H, 0, N - 2, N, 2);
      gsl_linalg_householder_mh(tau_i, w->hv2, &m.matrix);
    }

  if (w->Z)
    {
      /* accumulate transformation into Z */
      m = gsl_matrix_submatrix(w->Z, 0, top + N - 2, w->size, 2);
      gsl_linalg_householder_mh(tau_i, w->hv2, &m.matrix);
    }

  return GSL_SUCCESS;
} /* francis_qrstep() */

/*
francis_search_subdiag_small_elements()
  Search for a small subdiagonal element starting from the bottom
of a matrix A. A small element is one that satisfies:

|A_{i,i-1}| <= eps * (|A_{i,i}| + |A_{i-1,i-1}|)

Inputs: A - matrix (must be at least 3-by-3)

Return: row index of small subdiagonal element or 0 if not found

Notes: the first small element that is found (starting from bottom)
       is set to zero
*/

static  size_t
_francis_search_subdiag_small_elements(gsl_matrix * A)
{
  const size_t N = A->size1;
  size_t i;
  double dpel = gsl_matrix_get(A, N - 2, N - 2);

  for (i = N - 1; i > 0; --i)
    {
      double sel = gsl_matrix_get(A, i, i - 1);
      double del = gsl_matrix_get(A, i, i);

      if ((sel == 0.0) ||
          (fabs(sel) < GSL_DBL_EPSILON * (fabs(del) + fabs(dpel))))
        {
          gsl_matrix_set(A, i, i - 1, 0.0);
          return (i);
        }

      dpel = del;
    }

  return (0);
} /* francis_search_subdiag_small_elements() */

/*
 * This module contains some routines related to manipulating the
 * Schur form of a matrix which are needed by the eigenvalue solvers
 *
 * This file contains routines based on original code from LAPACK
 * which is distributed under the modified BSD license. The LAPACK
 * routine used is DLANV2.
 */

static  void _schur_standard_form(gsl_matrix *A, gsl_complex *eval1,
                                       gsl_complex *eval2, double *cs,
                                       double *sn);

/*
gsl_schur_standardize()
  Wrapper function for schur_standard_form - convert a 2-by-2 eigenvalue
block to standard form and then update the Schur form and
Schur vectors.

Inputs: T        - Schur form
        row      - row of T of 2-by-2 block to be updated
        eval1    - where to store eigenvalue 1
        eval2    - where to store eigenvalue 2
        update_t - 1 = update the entire matrix T with the transformation
                   0 = do not update rest of T
        Z        - (optional) if non-null, accumulate transformation
*/

static void
_gsl_schur_standardize(gsl_matrix *T, size_t row, gsl_complex *eval1,
                      gsl_complex *eval2, int update_t, gsl_matrix *Z)
{
  const size_t N = T->size1;
  gsl_matrix_view m;
  double cs, sn;

  m = gsl_matrix_submatrix(T, row, row, 2, 2);
  _schur_standard_form(&m.matrix, eval1, eval2, &cs, &sn);

  if (update_t)
    {
      gsl_vector_view xv, yv, v;

      /*
       * The above call to schur_standard_form transformed a 2-by-2 block
       * of T into upper triangular form via the transformation
       *
       * U = [ CS -SN ]
       *     [ SN  CS ]
       *
       * The original matrix T was
       *
       * T = [ T_{11} | T_{12} | T_{13} ]
       *     [   0*   |   A    | T_{23} ]
       *     [   0    |   0*   | T_{33} ]
       *
       * where 0* indicates all zeros except for possibly
       * one subdiagonal element next to A.
       *
       * After schur_standard_form, T looks like this:
       *
       * T = [ T_{11} | T_{12}  | T_{13} ]
       *     [   0*   | U^t A U | T_{23} ]
       *     [   0    |    0*   | T_{33} ]
       *
       * since only the 2-by-2 block of A was changed. However,
       * in order to be able to back transform T at the end,
       * we need to apply the U transformation to the rest
       * of the matrix T since there is no way to apply a
       * similarity transformation to T and change only the
       * middle 2-by-2 block. In other words, let
       *
       * M = [ I 0 0 ]
       *     [ 0 U 0 ]
       *     [ 0 0 I ]
       *
       * and compute
       *
       * M^t T M = [ T_{11} | T_{12} U |   T_{13}   ]
       *           [ U^t 0* | U^t A U  | U^t T_{23} ]
       *           [   0    |   0* U   |   T_{33}   ]
       *
       * So basically we need to apply the transformation U
       * to the i x 2 matrix T_{12} and the 2 x (n - i + 2)
       * matrix T_{23}, where i is the index of the top of A
       * in T.
       *
       * The BLAS routine drot() is suited for this.
       */

      if (row < (N - 2))
        {
          /* transform the 2 rows of T_{23} */

          v = gsl_matrix_row(T, row);
          xv = gsl_vector_subvector(&v.vector,
                                    row + 2,
                                    N - row - 2);

          v = gsl_matrix_row(T, row + 1);
          yv = gsl_vector_subvector(&v.vector,
                                    row + 2,
                                    N - row - 2);

          gsl_blas_drot(&xv.vector, &yv.vector, cs, sn);
        }

      if (row > 0)
        {
          /* transform the 2 columns of T_{12} */

          v = gsl_matrix_column(T, row);
          xv = gsl_vector_subvector(&v.vector,
                                    0,
                                    row);

          v = gsl_matrix_column(T, row + 1);
          yv = gsl_vector_subvector(&v.vector,
                                    0,
                                    row);

          gsl_blas_drot(&xv.vector, &yv.vector, cs, sn);
        }
    } /* if (update_t) */

  if (Z)
    {
      gsl_vector_view xv, yv;

      /*
       * Accumulate the transformation in Z. Here, Z -> Z * M
       *
       * So:
       *
       * Z -> [ Z_{11} | Z_{12} U | Z_{13} ]
       *      [ Z_{21} | Z_{22} U | Z_{23} ]
       *      [ Z_{31} | Z_{32} U | Z_{33} ]
       *
       * So we just need to apply drot() to the 2 columns
       * starting at index 'row'
       */

      xv = gsl_matrix_column(Z, row);
      yv = gsl_matrix_column(Z, row + 1);

      gsl_blas_drot(&xv.vector, &yv.vector, cs, sn);
    } /* if (Z) */
} /* gsl_schur_standardize() */

/*******************************************************
 *            INTERNAL ROUTINES                        *
 *******************************************************/

/*
schur_standard_form()
  Compute the Schur factorization of a real 2-by-2 matrix in
standard form:

[ A B ] = [ CS -SN ] [ T11 T12 ] [ CS SN ]
[ C D ]   [ SN  CS ] [ T21 T22 ] [-SN CS ]

where either:
1) T21 = 0 so that T11 and T22 are real eigenvalues of the matrix, or
2) T11 = T22 and T21*T12 < 0, so that T11 +/- sqrt(|T21*T12|) are
   complex conjugate eigenvalues

Inputs: A     - 2-by-2 matrix
        eval1 - where to store eigenvalue 1
        eval2 - where to store eigenvalue 2
        cs    - where to store cosine parameter of rotation matrix
        sn    - where to store sine parameter of rotation matrix

Notes: based on LAPACK routine DLANV2
*/

static  void
_schur_standard_form(gsl_matrix *A, gsl_complex *eval1, gsl_complex *eval2,
                    double *cs, double *sn)
{
  double a, b, c, d; /* input matrix values */
  double tmp;
  double p, z;
  double bcmax, bcmis, scale;
  double tau, sigma;
  double cs1, sn1;
  double aa, bb, cc, dd;
  double sab, sac;

  a = gsl_matrix_get(A, 0, 0);
  b = gsl_matrix_get(A, 0, 1);
  c = gsl_matrix_get(A, 1, 0);
  d = gsl_matrix_get(A, 1, 1);

  if (c == 0.0)
    {
      /*
       * matrix is already upper triangular - set rotation matrix
       * to the identity
       */
      *cs = 1.0;
      *sn = 0.0;
    }
  else if (b == 0.0)
    {
      /* swap rows and columns to make it upper triangular */

      *cs = 0.0;
      *sn = 1.0;

      tmp = d;
      d = a;
      a = tmp;
      b = -c;
      c = 0.0;
    }
  else if (((a - d) == 0.0) && (GSL_SIGN(b) != GSL_SIGN(c)))
    {
      /* the matrix has complex eigenvalues with a == d */
      *cs = 1.0;
      *sn = 0.0;
    }
  else
    {
      tmp = a - d;
      p = 0.5 * tmp;
      bcmax = GSL_MAX(fabs(b), fabs(c));
      bcmis = GSL_MIN(fabs(b), fabs(c)) * GSL_SIGN(b) * GSL_SIGN(c);
      scale = GSL_MAX(fabs(p), bcmax);
      z = (p / scale) * p + (bcmax / scale) * bcmis;

      if (z >= 4.0 * GSL_DBL_EPSILON)
        {
          /* real eigenvalues, compute a and d */

          z = p + GSL_SIGN(p) * fabs(sqrt(scale) * sqrt(z));
          a = d + z;
          d -= (bcmax / z) * bcmis;

          /* compute b and the rotation matrix */

          tau = gsl_hypot(c, z);
          *cs = z / tau;
          *sn = c / tau;
          b -= c;
          c = 0.0;
        }
      else
        {
          /*
           * complex eigenvalues, or real (almost) equal eigenvalues -
           * make diagonal elements equal
           */

          sigma = b + c;
          tau = gsl_hypot(sigma, tmp);
          *cs = sqrt(0.5 * (1.0 + fabs(sigma) / tau));
          *sn = -(p / (tau * (*cs))) * GSL_SIGN(sigma);

          /*
           * Compute [ AA BB ] = [ A B ] [ CS -SN ]
           *         [ CC DD ]   [ C D ] [ SN  CS ]
           */
          aa = a * (*cs) + b * (*sn);
          bb = -a * (*sn) + b * (*cs);
          cc = c * (*cs) + d * (*sn);
          dd = -c * (*sn) + d * (*cs);

          /*
           * Compute [ A B ] = [ CS SN ] [ AA BB ]
           *         [ C D ]   [-SN CS ] [ CC DD ]
           */
          a = aa * (*cs) + cc * (*sn);
          b = bb * (*cs) + dd * (*sn);
          c = -aa * (*sn) + cc * (*cs);
          d = -bb * (*sn) + dd * (*cs);

          tmp = 0.5 * (a + d);
          a = d = tmp;

          if (c != 0.0)
            {
              if (b != 0.0)
                {
                  if (GSL_SIGN(b) == GSL_SIGN(c))
                    {
                      /*
                       * real eigenvalues: reduce to upper triangular
                       * form
                       */
                      sab = sqrt(fabs(b));
                      sac = sqrt(fabs(c));
                      p = GSL_SIGN(c) * fabs(sab * sac);
                      tau = 1.0 / sqrt(fabs(b + c));
                      a = tmp + p;
                      d = tmp - p;
                      b -= c;
                      c = 0.0;

                      cs1 = sab * tau;
                      sn1 = sac * tau;
                      tmp = (*cs) * cs1 - (*sn) * sn1;
                      *sn = (*cs) * sn1 + (*sn) * cs1;
                      *cs = tmp;
                    }
                }
              else
                {
                  b = -c;
                  c = 0.0;
                  tmp = *cs;
                  *cs = -(*sn);
                  *sn = tmp;
                }
            }
        }
    }

  /* set eigenvalues */

  GSL_SET_REAL(eval1, a);
  GSL_SET_REAL(eval2, d);
  if (c == 0.0)
    {
      GSL_SET_IMAG(eval1, 0.0);
      GSL_SET_IMAG(eval2, 0.0);
    }
  else
    {
      tmp = sqrt(fabs(b) * fabs(c));
      GSL_SET_IMAG(eval1, tmp);
      GSL_SET_IMAG(eval2, -tmp);
    }

  /* set new matrix elements */

  gsl_matrix_set(A, 0, 0, a);
  gsl_matrix_set(A, 0, 1, b);
  gsl_matrix_set(A, 1, 0, c);
  gsl_matrix_set(A, 1, 1, d);
} /* schur_standard_form() */


/*
francis_schur_standardize()
  Convert a 2-by-2 diagonal block in the Schur form to standard form
and update the rest of T and Z matrices if required.

Inputs: A     - 2-by-2 matrix
        eval1 - where to store eigenvalue 1
        eval2 - where to store eigenvalue 2
        w     - francis workspace
*/

static  void
_francis_schur_standardize(gsl_matrix *A, gsl_complex *eval1,
                           gsl_complex *eval2,
                           _gsl_eigen_francis_workspace *w)
{
  size_t top;

  /*
   * figure out where the submatrix A resides in the
   * original matrix H
   */
  _francis_get_submatrix(w->H, A, &top);

  /* convert A to standard form and store eigenvalues */
  _gsl_schur_standardize(w->H, top, eval1, eval2, w->compute_t, w->Z);
} /* francis_schur_standardize() */

/*
francis_get_submatrix()
  B is a submatrix of A. The goal of this function is to
compute the indices in A of where the matrix B resides
*/

static  void
_francis_get_submatrix(gsl_matrix *A, gsl_matrix *B, size_t *top)
{
  size_t diff;
  double ratio;

  diff = (size_t) (B->data - A->data);

  ratio = (double)diff / ((double) (A->tda + 1));

  *top = (size_t) floor(ratio);
} /* francis_get_submatrix() */


/*
gsl_eigen_francis_alloc()

Allocate a workspace for solving the nonsymmetric eigenvalue problem.
The size of this workspace is O(1)

Inputs: none

Return: pointer to workspace
*/

static _gsl_eigen_francis_workspace *
_gsl_eigen_francis_alloc(void)
{
  _gsl_eigen_francis_workspace *w;

  w = (_gsl_eigen_francis_workspace *)
      malloc (sizeof (_gsl_eigen_francis_workspace));

  if (w == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for workspace", GSL_ENOMEM);
    }

  /* these are filled in later */

  w->size = 0;
  w->max_iterations = 0;
  w->n_iter = 0;
  w->n_evals = 0;

  w->compute_t = 0;
  w->Z = NULL;
  w->H = NULL;

  w->hv2 = gsl_vector_alloc(2);
  w->hv3 = gsl_vector_alloc(3);

  if ((w->hv2 == 0) || (w->hv3 == 0))
    {
      GSL_ERROR_NULL ("failed to allocate space for householder vectors", GSL_ENOMEM);
    }

  return (w);
} /* gsl_eigen_francis_alloc() */

/*
gsl_eigen_francis_free()
  Free francis workspace w
*/

static void
_gsl_eigen_francis_free (_gsl_eigen_francis_workspace *w)
{
  gsl_vector_free(w->hv2);
  gsl_vector_free(w->hv3);

  free(w);
} /* gsl_eigen_francis_free() */

/*
gsl_eigen_francis_T()
  Called when we want to compute the Schur form T, or no longer
compute the Schur form T

Inputs: compute_t - 1 to compute T, 0 to not compute T
        w         - francis workspace
*/

static void
_gsl_eigen_francis_T (const int compute_t, _gsl_eigen_francis_workspace *w)
{
  w->compute_t = compute_t;
}

/*
gsl_eigen_francis()

Solve the nonsymmetric eigenvalue problem

H x = \lambda x

for the eigenvalues \lambda using algorithm 7.5.2 of
Golub & Van Loan, "Matrix Computations" (3rd ed)

Inputs: H    - upper hessenberg matrix
        eval - where to store eigenvalues
        w    - workspace

Return: success or error - if error code is returned,
        then the QR procedure did not converge in the
        allowed number of iterations. In the event of non-
        convergence, the number of eigenvalues found will
        still be stored in the beginning of eval,

Notes: On output, the diagonal of H contains 1-by-1 or 2-by-2
       blocks containing the eigenvalues. If T is desired,
       H will contain the full Schur form on output.
*/

static int
_gsl_eigen_francis (gsl_matrix * H, gsl_vector_complex * eval,
                    _gsl_eigen_francis_workspace * w)
{
  /* check matrix and vector sizes */

  if (H->size1 != H->size2)
    {
      GSL_ERROR ("matrix must be square to compute eigenvalues", GSL_ENOTSQR);
    }
  else if (eval->size != H->size1)
    {
      GSL_ERROR ("eigenvalue vector must match matrix size", GSL_EBADLEN);
    }
  else
    {
      const size_t N = H->size1;
      int j;

      /*
       * Set internal parameters which depend on matrix size.
       * The Francis solver can be called with any size matrix
       * since the workspace does not depend on N.
       * Furthermore, multishift solvers which call the Francis
       * solver may need to call it with different sized matrices
       */
      w->size = N;
      w->max_iterations = 30 * N;

      /*
       * save a pointer to original matrix since francis_schur_decomp
       * is recursive
       */
      w->H = H;

      w->n_iter = 0;
      w->n_evals = 0;

      /*
       * zero out the first two subdiagonals (below the main subdiagonal)
       * needed as scratch space by the QR sweep routine
       */
      for (j = 0; j < (int) N - 3; ++j)
        {
          gsl_matrix_set(H, (size_t) j + 2, (size_t) j, 0.0);
          gsl_matrix_set(H, (size_t) j + 3, (size_t) j, 0.0);
        }

      if (N > 2)
        gsl_matrix_set(H, N - 1, N - 3, 0.0);

      /*
       * compute Schur decomposition of H and store eigenvalues
       * into eval
       */
      _francis_schur_decomp(H, eval, w);

      if (w->n_evals != N)
               return GSL_EMAXITER;

      return GSL_SUCCESS;
    }
} /* gsl_eigen_francis() */

/*
gsl_eigen_francis_Z()

Solve the nonsymmetric eigenvalue problem for a Hessenberg
matrix

H x = \lambda x

for the eigenvalues \lambda using the Francis double-shift
method.

Here we compute the real Schur form

T = Q^t H Q

with the diagonal blocks of T giving us the eigenvalues.
Q is the matrix of Schur vectors.

Originally, H was obtained from a general nonsymmetric matrix
A via a transformation

H = U^t A U

so that

T = (UQ)^t A (UQ) = Z^t A Z

Z is the matrix of Schur vectors computed by this algorithm

Inputs: H    - upper hessenberg matrix
        eval - where to store eigenvalues
        Z    - where to store Schur vectors
        w    - workspace

Notes: 1) If T is computed, it is stored in H on output. Otherwise,
          the diagonal of H will contain 1-by-1 and 2-by-2 blocks
          containing the eigenvalues.

       2) The matrix Z must be initialized to the Hessenberg
          similarity matrix U. Or if you want the eigenvalues
          of H, initialize Z to the identity matrix.
*/

static int
_gsl_eigen_francis_Z (gsl_matrix * H, gsl_vector_complex * eval,
                      gsl_matrix * Z, _gsl_eigen_francis_workspace * w)
{
  int s;

  /* set internal Z pointer so we know to accumulate transformations */
  w->Z = Z;

  s = _gsl_eigen_francis(H, eval, w);

  w->Z = NULL;

  return s;
} /* gsl_eigen_francis_Z() */


/*
 * This module computes the eigenvalues of a real nonsymmetric
 * matrix, using the double shift Francis method.
 *
 * See the references in francis.c.
 *
 * This module gets the matrix ready by balancing it and
 * reducing it to Hessenberg form before passing it to the
 * francis module.
 */

/*
gsl_eigen_nonsymm_alloc()

Allocate a workspace for solving the nonsymmetric eigenvalue problem.
The size of this workspace is O(2n)

Inputs: n - size of matrix

Return: pointer to workspace
*/

static _gsl_eigen_nonsymm_workspace *
_gsl_eigen_nonsymm_alloc(const size_t n)
{
  _gsl_eigen_nonsymm_workspace *w;

  if (n == 0)
    {
      GSL_ERROR_NULL ("matrix dimension must be positive integer",
                      GSL_EINVAL);
    }

  w = (_gsl_eigen_nonsymm_workspace *)
      malloc (sizeof (_gsl_eigen_nonsymm_workspace));

  if (w == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for workspace", GSL_ENOMEM);
    }

  w->size = n;
  w->Z = NULL;
  w->do_balance = 0;

  w->diag = gsl_vector_alloc(n);

  if (w->diag == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for balancing vector", GSL_ENOMEM);
    }

  w->tau = gsl_vector_alloc(n);

  if (w->tau == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for hessenberg coefficients", GSL_ENOMEM);
    }

  w->francis_workspace_p = _gsl_eigen_francis_alloc();

  if (w->francis_workspace_p == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for francis workspace", GSL_ENOMEM);
    }

  return (w);
} /* gsl_eigen_nonsymm_alloc() */

/*
gsl_eigen_nonsymm_free()
  Free workspace w
*/

static void
_gsl_eigen_nonsymm_free (_gsl_eigen_nonsymm_workspace * w)
{
  gsl_vector_free(w->tau);

  gsl_vector_free(w->diag);

  _gsl_eigen_francis_free(w->francis_workspace_p);

  free(w);
} /* gsl_eigen_nonsymm_free() */

/*
gsl_eigen_nonsymm_params()
  Set some parameters which define how we solve the eigenvalue
problem.

Inputs: compute_t - 1 if we want to compute T, 0 if not
        balance   - 1 if we want to balance the matrix, 0 if not
        w         - nonsymm workspace
*/

static void
_gsl_eigen_nonsymm_params (const int compute_t, const int balance,
                           _gsl_eigen_nonsymm_workspace *w)
{
  _gsl_eigen_francis_T(compute_t, w->francis_workspace_p);
  w->do_balance = balance;
} /* gsl_eigen_nonsymm_params() */

/*
gsl_eigen_nonsymm()

Solve the nonsymmetric eigenvalue problem

A x = \lambda x

for the eigenvalues \lambda using the Francis method.

Here we compute the real Schur form

T = Z^t A Z

with the diagonal blocks of T giving us the eigenvalues.
Z is a matrix of Schur vectors which is not computed by
this algorithm. See gsl_eigen_nonsymm_Z().

Inputs: A    - general real matrix
        eval - where to store eigenvalues
        w    - workspace

Return: success or error

Notes: If T is computed, it is stored in A on output. Otherwise
       the diagonal of A contains the 1-by-1 and 2-by-2 eigenvalue
       blocks.
*/

static int
_gsl_eigen_nonsymm (gsl_matrix * A, gsl_vector_complex * eval,
                    _gsl_eigen_nonsymm_workspace * w)
{
  const size_t N = A->size1;

  /* check matrix and vector sizes */

  if (N != A->size2)
    {
      GSL_ERROR ("matrix must be square to compute eigenvalues", GSL_ENOTSQR);
    }
  else if (eval->size != N)
    {
      GSL_ERROR ("eigenvalue vector must match matrix size", GSL_EBADLEN);
    }
  else
    {
      int s;

      if (w->do_balance)
        {
          /* balance the matrix */
          _gsl_linalg_balance_matrix(A, w->diag);
        }

      /* compute the Hessenberg reduction of A */
      _gsl_linalg_hessenberg(A, w->tau);

      if (w->Z)
        {
          /*
           * initialize the matrix Z to U, which is the matrix used
           * to construct the Hessenberg reduction.
           */

          /* compute U and store it in Z */
          _gsl_linalg_hessenberg_unpack(A, w->tau, w->Z);

          /* find the eigenvalues and Schur vectors */
          s = _gsl_eigen_francis_Z(A, eval, w->Z, w->francis_workspace_p);

          if (w->do_balance)
            {
              /*
               * The Schur vectors in Z are the vectors for the balanced
               * matrix. We now must undo the balancing to get the
               * vectors for the original matrix A.
               */
              _gsl_linalg_balance_accum(w->Z, w->diag);
            }
        }
      else
        {
          /* find the eigenvalues only */
          s = _gsl_eigen_francis(A, eval, w->francis_workspace_p);
        }

      w->n_evals = w->francis_workspace_p->n_evals;

      return s;
    }
} /* gsl_eigen_nonsymm() */


/* * * * * * * * * * * * * * * * * * *
 *
 *
 * 
 *   End of the pseudo-gsl code
 *
 *
 *
 * * * * * * * * * * * * * * * * * * *
 */ 


void
LALFindChirpPTFFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    )
/* </lalVerbatim> */
{
  UINT4                 i, j, k, l, kmax, kmin;
  UINT4                 errcode = 0;
  UINT4                 numPoints;
  UINT4                 deltaEventIndex;
  UINT4                 ignoreIndex;
  REAL4                 deltaT, max_eigen, r, s, x, y;
  REAL4                 deltaF, fFinal, fmin, length;
  COMPLEX8             *snr            = NULL;
  COMPLEX8             *PTFQtilde, *qtilde, *PTFq, *inputData;
  COMPLEX8Vector        qVec;
  FindChirpBankVetoData clusterInput;

  /* Define variables and allocate memory needed for gsl function */
  _gsl_eigen_nonsymm_workspace *workspace = _gsl_eigen_nonsymm_alloc ( 5 );
  gsl_matrix                  *PTFMatrix = gsl_matrix_alloc( 5, 5 );
  gsl_vector_complex          *eigenvalues = gsl_vector_complex_alloc( 5 );
  gsl_complex                  eval;
  
  /*
   *
   * point local pointers to input and output pointers
   *
   */

  /* number of points in a segment */
  numPoints = params->PTFqVec->vectorLength;
  /* workspace vectors */
  snr = params->PTFsnrVec->data;
  qtilde = params->qtildeVec->data;
  PTFq   = params->PTFqVec->data;
  qVec.length = numPoints;

  /* template and data */
  inputData = input->segment->data->data->data;
  length    = input->segment->data->data->length;
  PTFQtilde = input->fcTmplt->PTFQtilde->data;

  /* number of points and frequency cutoffs */
  deltaT = (REAL4) params->deltaT;
  deltaF = 1.0 / ( deltaT * (REAL4) numPoints );
  fFinal = (REAL4) input->fcTmplt->tmplt.fFinal;
  fmin   = (REAL4) input->fcTmplt->tmplt.fLower;
  kmax =  fFinal / deltaF < numPoints/2 ? fFinal / deltaF : numPoints/2;
  kmin =  fmin / deltaF > 1.0 ? fmin/ deltaF : 1;
  
  /* Set parameters for the gsl function */
  _gsl_eigen_nonsymm_params( 0, 1, workspace);
  
  INITSTATUS( status, "LALFindChirpPTFFilter", FINDCHIRPPTFFILTERC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* make sure the output handle exists, but points to a null pointer */
  ASSERT( eventList, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( !*eventList, status, FINDCHIRPH_ENNUL, FINDCHIRPH_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the filter parameters are reasonable */
  ASSERT( params->deltaT > 0, status,
      FINDCHIRPH_EDTZO, FINDCHIRPH_MSGEDTZO );
  ASSERT( params->rhosqThresh >= 0, status,
      FINDCHIRPH_ERHOT, FINDCHIRPH_MSGERHOT );

  /* check that the fft plan exists */
  ASSERT( params->invPlan, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the workspace vectors exist */

  /* if a rhosqVec vector has been created, check we can store data in it */
  if ( params->rhosqVec ) 
  {
    ASSERT( params->rhosqVec->data->data, status, 
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
    ASSERT( params->rhosqVec->data, status, 
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  }

  /* make sure that the input structure exists */
  ASSERT( input, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure that the input structure contains some input */
  ASSERT( input->fcTmplt, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( input->segment, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure that the filter has been initialized for the correct */
  /* approximant                                                    */
  if ( params->approximant != FindChirpPTF )
  {
    ABORT( status, FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  }

  /* make sure the approximant in the tmplt and segment agree */
  ASSERT( input->fcTmplt->tmplt.approximant == FindChirpPTF, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  ASSERT( input->segment->approximant == FindChirpPTF, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );

  fprintf(stderr,"LALFindChirpPTFFilterSegment called\n");

  /*
   *
   * compute viable search regions in the snrsq vector
   *
   */


  if ( input->fcTmplt->tmplt.tC <= 0 )
  {
    ABORT( status, FINDCHIRPH_ECHTZ, FINDCHIRPH_MSGECHTZ );
  }

  deltaEventIndex = (UINT4) rint( (input->fcTmplt->tmplt.tC / deltaT) + 1.0 );

  /* ignore corrupted data at start and end */
  ignoreIndex = ( input->segment->invSpecTrunc / 2 ) + deltaEventIndex;

  if ( lalDebugLevel & LALINFO )
  {
    CHAR infomsg[256];

    LALSnprintf( infomsg, sizeof(infomsg) / sizeof(*infomsg),
        "m1 = %e, m2 = %e, chi = %e, kappa = %e "
        "=> %e seconds => %d points\n"
        "invSpecTrunc = %d => ignoreIndex = %d\n", 
        input->fcTmplt->tmplt.mass1, input->fcTmplt->tmplt.mass2, 
        input->fcTmplt->tmplt.chi, input->fcTmplt->tmplt.kappa, 
        input->fcTmplt->tmplt.tC , deltaEventIndex, 
        input->segment->invSpecTrunc, ignoreIndex );
    LALInfo( status, infomsg );
  }

  /* XXX check that we are not filtering corrupted data XXX */
  /* XXX this is hardwired to 1/4 segment length        XXX */
  if ( ignoreIndex > numPoints / 4 )
  {
    ABORT( status, FINDCHIRPH_ECRUP, FINDCHIRPH_MSGECRUP );
  }
  /* XXX reset ignoreIndex to one quarter of a segment XXX */
  ignoreIndex = numPoints / 4;

  if ( lalDebugLevel & LALINFO )
  {
    CHAR infomsg[256];

    LALSnprintf( infomsg, sizeof(infomsg) / sizeof(*infomsg), 
        "filtering from %d to %d\n",
        ignoreIndex, numPoints - ignoreIndex );
    LALInfo( status, infomsg );
  }

  /*
   *
   * compute the PTF filter statistic
   *
   */

  /* clear the snr output vector and workspace*/
  memset( params->PTFsnrVec->data, 0, numPoints * sizeof(REAL4) );
  memset( params->PTFqVec->data, 0, 5 * numPoints * sizeof(COMPLEX8) );
  
  for ( i = 0; i < 5; ++i )
  {

    /* compute qtilde using data and Qtilde */

    memset( params->qtildeVec->data, 0, 
        params->qtildeVec->length * sizeof(COMPLEX8) );

    /* qtilde positive frequency, not DC or nyquist */
    for ( k = kmin; k < kmax ; ++k )
    {
      r = inputData[k].re;
      s = inputData[k].im;
      x = PTFQtilde[i * (numPoints / 2 + 1) + k].re;
      y = 0 - PTFQtilde[i * (numPoints / 2 + 1) + k].im; /* cplx conj */

      qtilde[k].re = 2 * (r*x - s*y);
      qtilde[k].im = 2 * (r*y + s*x);
    }

    qVec.data = params->PTFqVec->data + (i * numPoints);

    /* inverse fft to get q */
    LALCOMPLEX8VectorFFT( status->statusPtr, &qVec, params->qtildeVec, 
        params->invPlan );
    CHECKSTATUSPTR( status );
  }

  /* now we have PTFqVec which contains <s|Q^I_0> + i <s|Q^I_\pi/2> */

  for ( j = 0; j < numPoints; ++j ) /* beginning of main loop over time */
  {  
    /* Set PTFMatrix elements to zero */
    memset(params->PTFA->data, 0 , 25 * sizeof(REAL4));
    gsl_matrix_set_zero( PTFMatrix );
    
    /* construct A */
    for ( i = 0; i < 5; ++i )
    {  
      for ( l = 0; l < i + 1; ++l )
      { 
        params->PTFA->data[5 * i + l] = PTFq[i * numPoints + j].re * 
                                        PTFq[l * numPoints + j].re +
                                        PTFq[i * numPoints + j].im * 
                                        PTFq[l * numPoints + j].im ;
        params->PTFA->data[5 * l + i] = params->PTFA->data[ 5 * i + l]; 
      }  
    } 
    
    /* multiply by PTFBinverse to obtain AB^(-1) */

    LALSMatrixMultiply(status->statusPtr, 
        params->PTFMatrix, params->PTFA, 
        input->fcTmplt->PTFBinverse);
    CHECKSTATUSPTR( status );

    /* Construct the gsl matrix to be used by gsl_eigen_nonsymm */
    for ( i = 0; i < 5; ++i ) 
    {
      for ( l = 0; l < 5; ++l )
      {
        gsl_matrix_set(PTFMatrix, i, l, params->PTFMatrix->data[5 * i + l]);
      }
    }
     
    /* find max eigenvalue and store it in snr vector */
    errcode = _gsl_eigen_nonsymm ( PTFMatrix, eigenvalues, workspace); 
    if ( errcode != GSL_SUCCESS )
    {
      fprintf(stderr,"GSL eigenvalue error %d at time step %d\n",errcode,j);
      ABORT( status, FINDCHIRPH_EIGEN, FINDCHIRPH_MSGEIGEN);
    }
    
    max_eigen = 0.0;
    for ( i = 0; i < 5; ++i )
    {
      eval = gsl_vector_complex_get( eigenvalues, i );
      if ( GSL_IMAG( eval ) == 0) 
      {  
        if ( GSL_REAL( eval ) > max_eigen ) max_eigen = GSL_REAL( eval );        
      } 
    }
    
    snr[j].re = 2.0 * sqrt(max_eigen) / (REAL4) numPoints;
    snr[j].im = 0;
    
  } /* End of main loop over time */


  /* 
   *
   * calculate signal to noise squared 
   *
   */


  /* if full snrsq vector is required, set it to zero */
  if ( params->rhosqVec )
    memset( params->rhosqVec->data->data, 0, numPoints * sizeof( REAL4 ) );

  /* if full snrsq vector is required, store the snrsq */
  if ( params->rhosqVec ) 
  {
    memcpy( params->rhosqVec->name, input->segment->data->name,
        LALNameLength * sizeof(CHAR) );
    memcpy( &(params->rhosqVec->epoch), &(input->segment->data->epoch), 
        sizeof(LIGOTimeGPS) );
    params->rhosqVec->deltaT = input->segment->deltaT;

    for ( j = 0; j < numPoints; ++j )
    {
      params->rhosqVec->data->data[j] =  snr[j].re * snr[j].re;
    }
  }


  /*
   *
   * look for and cluster events in the snr vector
   *
   */


  clusterInput.length       = 1; /* do not do the bank veto */
  clusterInput.qVecArray    = NULL;
  clusterInput.fcInputArray = NULL;
  clusterInput.ccMat        = NULL;
  clusterInput.normMat      = NULL;
  clusterInput.spec         = NULL;
  clusterInput.resp         = NULL;
  params->ignoreIndex       = ignoreIndex;

  /* cluster events searches the qVec for events */
  params->qVec = params->PTFsnrVec;
  input->fcTmplt->norm = 1.0;

  LALFindChirpClusterEvents( status->statusPtr, eventList, 
      input, params, &clusterInput, 0 );
  CHECKSTATUSPTR( status );

  params->qVec = NULL;

  /* Free gsl allocated memory */
  gsl_matrix_free( PTFMatrix );
  _gsl_eigen_nonsymm_free( workspace );
  gsl_vector_complex_free( eigenvalues );
  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
