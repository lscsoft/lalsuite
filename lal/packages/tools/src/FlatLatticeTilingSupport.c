/*
 *  Copyright (C) 2007, 2008 Karl Wette
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
 * \brief Support routines for FlatLatticeTiling
 */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_combination.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_nan.h>

#include <lal/LALRCSID.h>
#include <lal/LALStdlib.h>
#include <lal/LALMalloc.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/FlatLatticeTiling.h>
#include <lal/FlatLatticeTilingSupport.h>

NRCSID(FLATLATTICETILINGSUPPORTC, "$Id$");

static const BOOLEAN TRUE  = (1==1);
static const BOOLEAN FALSE = (0==1);

/** 
 * Convert variable argument list to a gsl_vector
 */
gsl_vector *XLALGSLVectorFromVAList(
				    INT4 nargs,  /**< Number of arguments */
				    REAL8 first, /**< First argument */
				    ...          /**< Additional arguments */
				    )
{

  int i;
  gsl_vector *result = NULL;
  va_list va;

  /* Allocate memory */
  ALLOC_GSL_VECTOR(result, nargs, XLAL_ERROR_NULL);

  /* Initialise vector */
  gsl_vector_set(result, 0, first);
  va_start(va, first);
  for (i = 1; i < nargs; ++i) {
    gsl_vector_set(result, i, va_arg(va, REAL8));
  }
  va_end(va);

  return result;

}

/** 
 * Convert variable argument list to a gsl_vector
 */
gsl_vector *XLALGSLVectorFromLALStringVector(
					     LALStringVector *args /* Arguments */
					     )
{

  int i;
  gsl_vector *result = NULL;
  double x = 0.0;

  /* Allocate memory */
  ALLOC_GSL_VECTOR(result, args->length, XLAL_ERROR_NULL);

  /* Initialise vector */
  for (i = 0; i < (INT4)args->length; ++i) {
    if (sscanf(args->data[i], "%le", &x) != 1) {
      XLALPrintError("'%s' is not numeric\n", args->data[i]);
      XLAL_ERROR_NULL("An element of 'args' is not numeric", XLAL_EINVAL);
    }
    gsl_vector_set(result, i, x);
  }

  return result;

}

/**
 * Perform a Gaussian pivot around a matrix element
 */
void XLALGaussianPivot(
		       gsl_matrix* m, /**< Matrix */
		       size_t ii,     /**< Pivot row */
		       size_t jj      /**< Pivot column */
		       )
{
  
  size_t i, j;
  double x = 0.0;
  
  /* Normalise row */
  x = gsl_matrix_get(m, ii, jj);
  for (j = 0; j < m->size2; ++j)
    gsl_matrix_set(m, ii, j, gsl_matrix_get(m, ii, j) / x);

  /* Subtract from all other rows */
  for (i = 0; i < m->size1; ++i) {
    if (i == ii)
      continue;
    x = gsl_matrix_get(m, i, jj);
    for (j = 0; j < m->size2; ++j)
      gsl_matrix_set(m, i, j, gsl_matrix_get(m, i, j) - x * gsl_matrix_get(m, ii, j));
  }

}

/**
 * Reduce a matrix to reduced row echelon form 
 * by Gaussian elimination with partial pivoting
 */
void XLALReducedRowEchelonForm(
			       gsl_matrix *m /**< Matrix */
			       )
{
  
  size_t ii, jj, i, maxi;
  double x = 0.0;
  
  /* Loop over all columns */
  for (ii = 0, jj = 0; ii < m->size1 && jj < m->size2; ++jj) {
    
    /* Find the pivot: row with largest value */
    maxi = ii;
    for (i = ii + 1; i < m->size1; ++i)
      if (fabs(gsl_matrix_get(m, i, jj)) > fabs(gsl_matrix_get(m, maxi, jj)))
	maxi = i;
    
    /* If non-zero, swap rows, pivot and move to next row */
    if (gsl_matrix_get(m, maxi, jj) != 0.0) {
      gsl_matrix_swap_rows(m, maxi, ii);
      XLALGaussianPivot(m, ii, jj);
      ++ii;
    }

  }
  
}

/**
 * Resize a gsl_matrix and copy over existing contents
 */
gsl_matrix *XLALResizeGSLMatrix(
				gsl_matrix *m, /**< Matrix */
				size_t size1,  /**< New number of rows */
				size_t size2   /**< New number of columns */
				)
{

  gsl_matrix *n = NULL;
  
  /* Allocate memory */
  ALLOC_GSL_MATRIX(n, size1, size2, XLAL_ERROR_NULL);
  
  /* Copy contents if any */
  if (m != NULL) {
    gsl_matrix_view old = gsl_matrix_submatrix(m, 0, 0, GSL_MIN(m->size1, size1), GSL_MIN(m->size2, size2));
    gsl_matrix_view new = gsl_matrix_submatrix(n, 0, 0, GSL_MIN(m->size1, size1), GSL_MIN(m->size2, size2));
    gsl_matrix_memcpy(&new.matrix, &old.matrix);
    FREE_GSL_MATRIX(m);
  }

  return n;
  
}

/**
 * Resize a gsl_vector and copy over existing contents
 */
gsl_vector *XLALResizeGSLVector(
				gsl_vector *u, /**< Vector */
				size_t size    /**< New number of elements */
				)
{
  
  gsl_vector *v = NULL;

  /* Allocate memory */
  ALLOC_GSL_VECTOR(v, size, XLAL_ERROR_NULL);
  
  /* Copy contents if any */
  if (u != NULL) {
    gsl_vector_view old = gsl_vector_subvector(u, 0, GSL_MIN(u->size, size));
    gsl_vector_view new = gsl_vector_subvector(v, 0, GSL_MIN(u->size, size));
    gsl_vector_memcpy(&new.vector, &old.vector);
    FREE_GSL_VECTOR(u);
  }
  
  return v;

}

/**
 * Find the vertices of a given polytope \f$A x \le b\f$
 */
int XLALPolytopeVertices(
			 gsl_matrix *polytope_A, /**< Polytope matrix A */
			 gsl_vector *polytope_b, /**< Polytope vector b */
			 gsl_matrix **vertices   /**< Polytope vertices */
			 ) 
{

  const int m = polytope_A->size1;
  const int n = polytope_A->size2;
  
  int i;
  gsl_combination *vertex_hyperplanes = NULL;
  gsl_matrix *vertex_eqns = NULL;
  gsl_permutation *vertex_perm = NULL;
  int vertex_sign = 0;
  gsl_vector *vertex_soln = NULL;
  gsl_vector *vertex_b = NULL;
  int nvertices = 0;
  BOOLEAN passed = FALSE;

  /* Allocate memory */
  ALLOC_GSL_COMBINATION(vertex_hyperplanes, m, n, XLAL_ERROR);
  ALLOC_GSL_MATRIX(vertex_eqns, n, n, XLAL_ERROR);
  ALLOC_GSL_PERMUTATION(vertex_perm, n, XLAL_ERROR);
  ALLOC_GSL_VECTOR(vertex_soln, n, XLAL_ERROR);
  ALLOC_GSL_VECTOR(vertex_b, m, XLAL_ERROR);

  /* Loop over all possible combinations of hyperplanes */
  nvertices = 0;
  *vertices = NULL;
  gsl_combination_init_first(vertex_hyperplanes);
  do {

    /* Initialise equation matrix for vertex */
    for (i = 0; i < n; ++i) {
      gsl_vector_view row = gsl_matrix_row(vertex_eqns, i);
      gsl_vector_view hyperplane = gsl_matrix_row(polytope_A, gsl_combination_get(vertex_hyperplanes, i));
      gsl_vector_memcpy(&row.vector, &hyperplane.vector);
      gsl_vector_set(vertex_soln, i, gsl_vector_get(polytope_b, gsl_combination_get(vertex_hyperplanes, i)));
    }

    /* Solve for vertex */
    gsl_linalg_LU_decomp(vertex_eqns, vertex_perm, &vertex_sign);
    gsl_linalg_LU_svx(vertex_eqns, vertex_perm, vertex_soln);

    /* If vertex is finite */
    passed = TRUE;
    for (i = 0; i < n; ++i)
      passed &= (gsl_finite(gsl_vector_get(vertex_soln, i)) == 1);
    if (passed) {

      /* If vertex satisfied all hyperplanes */
      gsl_blas_dgemv(CblasNoTrans, 1.0, polytope_A, vertex_soln, 0.0, vertex_b);
      passed = TRUE;
      for (i = 0; i < n; ++i)
	passed &= (gsl_vector_get(vertex_b, i) <= gsl_vector_get(polytope_b, i));
      if (passed) {

	/* Add vertex */
	++nvertices;
	if ((*vertices = XLALResizeGSLMatrix(*vertices, n, nvertices)) == NULL)
	  XLAL_ERROR("XLALResizeMatrix failed", XLAL_EFAILED);
	{
	  gsl_vector_view last_col = gsl_matrix_column(*vertices, nvertices - 1);
	  gsl_vector_memcpy(&last_col.vector, vertex_soln);
	}

      }

    }

  } while (gsl_combination_next(vertex_hyperplanes) == GSL_SUCCESS);

  /* Cleanup */
  FREE_GSL_COMBINATION(vertex_hyperplanes);
  FREE_GSL_MATRIX(vertex_eqns);
  FREE_GSL_PERMUTATION(vertex_perm);
  FREE_GSL_VECTOR(vertex_soln);
  FREE_GSL_VECTOR(vertex_b);

  return (vertices == NULL) ? XLAL_FAILURE : XLAL_SUCCESS;

}

/**
 * Given vertices, find the corresponding polytope \f$A x \le b\f$
 */
int XLALPolytopeFromVertices(
			     gsl_matrix *vertices,    /**< Vertices */
			     gsl_matrix **polytope_A, /**< Polytope matrix A */
			     gsl_vector **polytope_b  /**< Polytope vector b */
			     )
{

  const int m = vertices->size1;
  const int n = vertices->size2;
  
  int i, j;
  gsl_combination *hyperplane_vertices = NULL;
  gsl_matrix *hyperplane_eqns = NULL;
  gsl_vector *hyperplane_soln = NULL;
  gsl_matrix *hyperplane_A = NULL;
  double hyperplane_b = 0.0;
  gsl_matrix *vertices_test = NULL;
  double vertices_test_min = 0.0;
  double vertices_test_max = 0.0;
  int nhyperplanes = 0;
  BOOLEAN found = FALSE;

  /* Allocate memory */
  ALLOC_GSL_COMBINATION(hyperplane_vertices, n, m, XLAL_ERROR);
  ALLOC_GSL_MATRIX(hyperplane_eqns, m, m + 1, XLAL_ERROR);
  ALLOC_GSL_VECTOR(hyperplane_soln, m + 1, XLAL_ERROR);
  ALLOC_GSL_MATRIX(hyperplane_A, 1, m, XLAL_ERROR);
  ALLOC_GSL_MATRIX(vertices_test, 1, n, XLAL_ERROR);

  /* Loop over all possible combinations of vertices */
  nhyperplanes = 0;
  *polytope_A = NULL;
  *polytope_b = NULL;
  gsl_combination_init_first(hyperplane_vertices);
  do {

    /* Initialise equation matrix for hyperplane */
    for (i = 0; i < m; ++i) {
      for (j = 0; j < m; ++j)
	gsl_matrix_set(hyperplane_eqns, i, j, gsl_matrix_get(vertices, j, gsl_combination_get(hyperplane_vertices, i)));
      gsl_matrix_set(hyperplane_eqns, i, m, 1.0);
    }

    /* Reduce equations to row reduced echelon */
    XLALReducedRowEchelonForm(hyperplane_eqns);
    
    /* Solve for hyperplane solution coefficients */
    gsl_vector_set_all(hyperplane_soln, 1.0);
    for (i = m - 1; i >= 0; --i) {
      gsl_vector_view row = gsl_matrix_row(hyperplane_eqns, i);
      double x = 0.0;
      for (j = 0; gsl_vector_get(&row.vector, j) == 0.0; ++j);
      gsl_vector_set(hyperplane_soln, j, 0.0);
      gsl_blas_ddot(&row.vector, hyperplane_soln, &x);
      gsl_vector_set(hyperplane_soln, j, -1.0 * x / gsl_vector_get(&row.vector, j));
    }

    /* Separate variable coefficients and constant term */
    for (j = 0; j < m; ++j)
      gsl_matrix_set(hyperplane_A, 0, j, gsl_vector_get(hyperplane_soln, j));
    hyperplane_b = -1.0 * gsl_vector_get(hyperplane_soln, m);      

    /* Compute range of A*x for all vertices  */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, hyperplane_A, vertices, 0.0, vertices_test);
    gsl_matrix_minmax(vertices_test, &vertices_test_min, &vertices_test_max);

    /* If all vertices fall on one side of constant, vertices are within half-space */
    if (vertices_test_max <= hyperplane_b || hyperplane_b <= vertices_test_min) {

      /* If points are greater, reverse sign */
      if (hyperplane_b <= vertices_test_min) {
	gsl_matrix_scale(hyperplane_A, -1.0);
	hyperplane_b *= -1.0;
      }

      /* Check that half-space has not already been added */
      found = FALSE;
      if ((*polytope_A) != NULL && (*polytope_b) != NULL) {
	for (i = 0; i < (int)(*polytope_A)->size1; ++i) {
	  found = TRUE;
	  for (j = 0; j < (int)(*polytope_A)->size2; ++j)
	    found &= (gsl_matrix_get(*polytope_A, i, j) == gsl_matrix_get(hyperplane_A, 0, j));
	  found &= (gsl_vector_get(*polytope_b, i) == hyperplane_b);
	  if (found)
	    break;
	}
      }
      if (!found) {

	/* Increase size of polytope and add new half-space */
	++nhyperplanes;
	if ((*polytope_A = XLALResizeGSLMatrix(*polytope_A, nhyperplanes, m)) == NULL)
	  XLAL_ERROR("XLALResizeGSLMatrix failed", XLAL_EFAILED);
	if ((*polytope_b = XLALResizeGSLVector(*polytope_b, nhyperplanes)) == NULL)
	  XLAL_ERROR("XLALResizeGSLVector failed", XLAL_EFAILED);
	{
	  gsl_matrix_view last_row = gsl_matrix_submatrix(*polytope_A, nhyperplanes - 1, 0, 1, m);
	  gsl_matrix_memcpy(&last_row.matrix, hyperplane_A);
	}
	gsl_vector_set(*polytope_b, nhyperplanes - 1, hyperplane_b);

      }

    }

  } while (gsl_combination_next(hyperplane_vertices) == GSL_SUCCESS);

  /* Cleanup */
  FREE_GSL_COMBINATION(hyperplane_vertices);
  FREE_GSL_MATRIX(hyperplane_eqns);
  FREE_GSL_VECTOR(hyperplane_soln);
  FREE_GSL_MATRIX(hyperplane_A);
  FREE_GSL_MATRIX(vertices_test);
  
  return (polytope_A == NULL) ? XLAL_FAILURE : XLAL_SUCCESS;

}

/**
 * Create a simplex method tableau.
 */
static SimplexMethodTableau *CreateSimplexMethodTableau(
							int m, /**< Number of constraints */
							int n  /**< Dimension of solution space */
							)
{

  const int N = 2*n + m;

  SimplexMethodTableau *tableau = NULL;

  /* Check dimensions */
  if (m < 1 || n < 1)
    XLAL_ERROR_NULL("'m'/'n' must be non-zero", XLAL_EINVAL);

  /* Allocate memory */
  if ((tableau = (SimplexMethodTableau*) LALMalloc(sizeof(SimplexMethodTableau))) == NULL)
    XLAL_ERROR_NULL("Could not allocate 'tableau'", XLAL_ENOMEM);
  ALLOC_GSL_MATRIX(tableau->T, 1 + N, 1 + 3*N, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(tableau->work_T, tableau->T->size1, tableau->T->size2, XLAL_ERROR_NULL);
  ALLOC_GSL_VECTOR_INT(tableau->basic, tableau->T->size2, XLAL_ERROR_NULL);

  /* Initialise structure */
  tableau->m = m;
  tableau->n = n;
  gsl_matrix_set_zero(tableau->T);

  return tableau;  

}

/**
 * Construct the simplex method tableau for the problem:
 * Find \f$x\f$ that minimizes/maximizes the quadratic
 * \f$f(x) = x^T Q x + c^T x\f$
 * subject to the constraint \f$A x \le b\f$.
 */
int XLALQuadraticSimplexMethodTableau(
				      int optimize,                  /**< >0 for maximize, <0 for minimize, =0 to ignore */
				      gsl_matrix *Q,                 /**< nxn matrix, =NULL to ignore */
				      gsl_vector *c,                 /**< n vector, =NULL to ignore */
				      gsl_matrix *A,                 /**< mxn matrix, =NULL to ignore */
				      gsl_vector *b,                 /**< m vector, =NULL to ignore */   
				      SimplexMethodTableau **tableau /**< Simplex method tableau */
				      )
{

  const BOOLEAN initialised = (*tableau != NULL);
  
  /* Allocate tableau */
  if (!initialised) {
    if (!Q || !A)
      XLAL_ERROR("'Q' and 'A' must be supplied for initialisation of 'tableau->T'", XLAL_EINVAL);
    if ((*tableau = CreateSimplexMethodTableau((int)A->size1, (int)A->size2)) == NULL)
      XLAL_ERROR("CreateSimplexMethodTableau failed", XLAL_EFAILED);
  }
  {

    const int m = (*tableau)->m;
    const int n = (*tableau)->n;
    const int N = 2*n + m;

    int i, j;

    /* Initialise tableau matrix */
    if (!initialised) {
      
      /* Initialise coefficients of artificial 
	 variables in the objective (first) row */
      {
	gsl_vector_view u = gsl_matrix_row((*tableau)->T, 0);
	gsl_vector_view v = gsl_vector_subvector(&u.vector, 1, N);
	gsl_vector_set_all(&v.vector, 1.0);
      }
      
      /* Initialise coefficients of y and v */
      {
	gsl_matrix_view U = gsl_matrix_submatrix((*tableau)->T, 1, 1 + 2*N, N, N);
	gsl_matrix_view V = gsl_matrix_submatrix(&U.matrix, 0, 0, 2*n, 2*n);
	gsl_matrix_set_identity(&U.matrix);
	gsl_matrix_scale(&V.matrix, -1.0);
      }
    }

    /* Normalise optimize */
    if (optimize != 0)
      optimize = (optimize < 0 ? -1 : 1);
    
    /* Update c in (*tableau)->T */
    if (c) {
      if ((int)c->size != n)
	XLAL_ERROR("'c' is not of length 'n'", XLAL_EINVAL);
      if (optimize == 0)
	XLAL_ERROR("'optimize' must be non-zero if 'c' is supplied", XLAL_EINVAL);
      {
	gsl_vector_view u = gsl_matrix_column((*tableau)->T, 0);
	gsl_vector_view v = gsl_vector_subvector(&u.vector, 1, n);
	gsl_vector_view w = gsl_vector_subvector(&u.vector, n + 1, n);
	gsl_vector_memcpy(&v.vector, c);
	gsl_vector_memcpy(&w.vector, c);
	gsl_vector_scale(&v.vector,  1.0 * optimize);
	gsl_vector_scale(&w.vector, -1.0 * optimize);
      }
    }
    
    /* Update b in (*tableau)->T */
    if (b) {
      if ((int)b->size != m)
	XLAL_ERROR("'b' is not of length 'm'", XLAL_EINVAL);
      {
	gsl_vector_view u = gsl_matrix_column((*tableau)->T, 0);
	gsl_vector_view v = gsl_vector_subvector(&u.vector, 2*n + 1, m);
	gsl_vector_memcpy(&v.vector, b);
      }
    }
    
    /* Update coefficients of artificial variables */
    if (c || b) {
      gsl_matrix_view U = gsl_matrix_submatrix((*tableau)->T, 1, 1, N, N);
      for (i = 0; i < N; ++i) {
	gsl_matrix_set(&U.matrix, i, i, gsl_matrix_get((*tableau)->T, i + 1, 0) < 0 ? -1.0 : 1.0);
      }
    }
    
    /* Update Q */
    if (Q) {
      if ((int)Q->size1 != n || (int)Q->size2 != n)
	XLAL_ERROR("'Q' is not of size 'n'x'n'", XLAL_EINVAL);
      if (optimize == 0)
	XLAL_ERROR("'optimize' must be non-zero if 'Q' is supplied", XLAL_EINVAL);
      {
	gsl_matrix_view U = gsl_matrix_submatrix((*tableau)->T, 1, 1 + N, N, N);
	for (i = 0; i < 2; ++i) {
	  for (j = 0; j < 2; ++j) {
	    gsl_matrix_view V = gsl_matrix_submatrix(&U.matrix, i*n, j*n, n, n);
	    gsl_matrix_memcpy(&V.matrix, Q);
	    gsl_matrix_scale(&V.matrix, -2.0 * optimize * pow(-1.0, i - j));
	  }
	}
      }
    }
    
    /* Update A */
    if (A) {
      if ((int)A->size1 != m || (int)A->size2 != n)
	XLAL_ERROR("'A' is not of size 'm'x'n'", XLAL_EINVAL);
      {
	gsl_matrix_view U = gsl_matrix_submatrix((*tableau)->T, 1, 1 + N, N, N);
	for (j = 0; j < 2; ++j) {
	  gsl_matrix_view V = gsl_matrix_submatrix(&U.matrix, 2*n, j*n, m, n);
	  gsl_matrix_memcpy(&V.matrix, A);
	  gsl_matrix_scale(&V.matrix, pow(-1.0, j));
	}
	{
	  gsl_matrix_view V  = gsl_matrix_submatrix(&U.matrix, 2*n, 0, m, 2*n);
	  gsl_matrix_view Vt = gsl_matrix_submatrix(&U.matrix, 0, 2*n, 2*n, m);
	  for (i = 0; i < m; ++i)
	    for (j = 0; j < 2*n; ++j)
	      gsl_matrix_set(&Vt.matrix, j, i, gsl_matrix_get(&V.matrix, i, j));
	}
      }
    }

  }

  return XLAL_SUCCESS;

}

/**
 * If vector represents a basic variable, return its index
 * otherwise return -1
 */
static int BasicVariable(gsl_vector *v) {

  int i;
  double x;

  /* Get length of vector */
  gsl_blas_ddot(v, v, &x);

  /* If vector has unit length */
  if (gsl_fcmp(x, 1.0, GSL_DBL_EPSILON) == 0)
    /* Search for unity element */
    for (i = 0; i < (int)v->size; ++i)
      if (gsl_fcmp(gsl_vector_get(v, i), 1.0, GSL_DBL_EPSILON) == 0)
	/* Return index of basic variable */
	return i;

  /* Not a basic variable */
  return -1;

}

/**
 * Find \f$x\f$ that minimizes/maximizes the quadratic
 * \f$f(x) = x^T Q x + c^T x\f$
 * subject to the constraint \f$A x \le b\f$, with the given tableau.
 */
int XLALQuadraticSimplexMethod(
			       SimplexMethodTableau *tableau, /**< Simplex method tableau */
			       gsl_vector *optimal            /**< Optimal solution */
			       )
{

  const int m = tableau->m;
  const int n = tableau->n;
  const int N = 2*n + m;

  int i, j, s, ii, jj;
  double x = 0.0;

  /* Check dimensions */
  if (!optimal)
    XLAL_ERROR("'optimal' must be allocated", XLAL_EINVAL);
  if ((int)optimal->size != n)
    XLAL_ERROR("'optimal' must be of length 'n'", XLAL_EINVAL);

  /* Create working copy of tableau */
  gsl_matrix_memcpy(tableau->work_T, tableau->T);

  /* Reduce tableau to its proper form from Gaussian elimination
     by pivoting around the coefficients of the artificial variables */
  for (i = 0; i < N; ++i)
    XLALGaussianPivot(tableau->work_T, 1 + i, 1 + i);
  
  /* Begin iterations */
  gsl_vector_int_set_all(tableau->basic, 0);
  while (TRUE) {

    /* Find the entering basic variable (column) with the largest 
       strictly negative coefficient in the objective (first) row, if any.
       Exclude variables whose complementary variable is basic, to satisfy 
       the complementary constraint. If no such variable exists, the optimal
       solution has been found. */
    x = 0.0;
    jj = -1;
    for (j = 1 + N; j < (int)tableau->work_T->size2; ++j) {

      /* Column of the complementary variable */
      const int jc = 1 + N + ((j - 1) % (2*N));
      
      /* If neither variable nor its complementary are basic */
      if (gsl_vector_int_get(tableau->basic, j) == 0 && 
	  gsl_vector_int_get(tableau->basic, jc) == 0) {

	/* Get value and store smallest */
	double y = gsl_matrix_get(tableau->work_T, 0, j);
	if (y < x) {
	  x = y;
	  jj = j;
	}

      }
    }
    
    /* Optimal solution has been found */
    if (jj < 0)
      break;
      
    /* Find the leaving basic variable (row) which is the 
       smallest of T(i,0)/T(i,jj), excluding rows which are
       not basic variables and rows where T(i,jj) <= 0. 
       If no such row, no solution exists. */
    x = GSL_POSINF;
    ii = -1;
    for (i = 1; i < 1 + N; ++i) {
	
      /* Check pivot element */
      double y = gsl_matrix_get(tableau->work_T, i, jj);
      if (y > 0) {
	
	/* Get ratio and store smallest */
	y = gsl_matrix_get(tableau->work_T, i, 0) / y;
	if (y < x) {
	  x = y;
	  ii = i;
	}
	
      }
    }
    
    /* No solution exists */
    if (ii < 0)
      return XLAL_FAILURE;
    
    /* Record entering and leaving basic variables */
    for (j = 0; j < (int)tableau->work_T->size2; ++j)
      if (j == jj)
	gsl_vector_int_set(tableau->basic, j, ii);
      else if (gsl_vector_int_get(tableau->basic, j) == ii)
	gsl_vector_int_set(tableau->basic, j, 0);

    /* Perform Gaussian pivot to change basic variable from ii to jj */
    XLALGaussianPivot(tableau->work_T, ii, jj);
    
    /* Combat roundoff error */
    x = 0.0;
    for (i = 0; i < (int)tableau->work_T->size1; ++i)
      for (j = 0; j < (int)tableau->work_T->size2; ++j) {
	double y = fabs(gsl_matrix_get(tableau->work_T, i, j));
	if (y > x)
	  x = y;
      }
    for (i = 0; i < (int)tableau->work_T->size1; ++i)
      for (j = 0; j < (int)tableau->work_T->size2; ++j)
	if (fabs(gsl_matrix_get(tableau->work_T, i, j) / x) < GSL_DBL_EPSILON)
	  gsl_matrix_set(tableau->work_T, i, j, 0.0);
    
  }
  
  /* Find the optimal solution: take values from first
     column if is a basic variable, otherwise take as zero. */
  for (j = 0; j < n; ++j) {
    x = 0.0;
    for (s = 0; s < 2; ++s) {
      i = gsl_vector_int_get(tableau->basic, 1 + N + s*n + j);
      if (i > 0)
	x += (s == 0 ? 1.0 : -1.0) * gsl_matrix_get(tableau->work_T, i, 0);
    }
    gsl_vector_set(optimal, j, x);
  }

  return XLAL_SUCCESS;

}

/**
 * Free a simplex method tableau
 */
void XLALFreeSimplexMethodTableau(
				  SimplexMethodTableau *tableau /**< Simplex method tableau */
				  )
{
  
  if (tableau) {
    
    FREE_GSL_MATRIX(tableau->T);
    FREE_GSL_MATRIX(tableau->work_T);
    FREE_GSL_VECTOR_INT(tableau->basic);
    
    LALFree(tableau);

  }
  
}

/**
 * Find the bounding box of the mismatch ellipses of a metric
 */
gsl_vector *XLALMetricEllipseBoundingBox(
					 gsl_matrix *metric, /**< Metric to bound */
					 REAL8 max_mismatch  /**< Maximum mismatch w.r.t metric */
					 )
{

  const int n = metric->size1;

  int i = 0;
  gsl_matrix *LU_decomp = NULL;
  gsl_permutation *LU_perm = NULL;
  gsl_matrix *inverse = NULL;
  int LU_sign = 0;
  gsl_vector *bound_box = NULL;

  /* Check input */
  if ((int)metric->size1 != n || (int)metric->size2 != n)
    XLAL_ERROR_NULL("'metric' is not square", XLAL_ESIZE);

  /* Allocate memory */
  ALLOC_GSL_MATRIX(LU_decomp, n, n, XLAL_ERROR_NULL);
  ALLOC_GSL_PERMUTATION(LU_perm, n, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(inverse, n, n, XLAL_ERROR_NULL);
  ALLOC_GSL_VECTOR(bound_box, n, XLAL_ERROR_NULL);

  /* Compute metric inverse */
  gsl_matrix_memcpy(LU_decomp, metric);
  gsl_linalg_LU_decomp(LU_decomp, LU_perm, &LU_sign);
  gsl_linalg_LU_invert(LU_decomp, LU_perm, inverse);

  /* Compute bounding box */
  for (i = 0; i < n; ++i) {
    gsl_vector_set(bound_box, i, sqrt(max_mismatch * gsl_matrix_get(inverse, i ,i)));
  }

  /* Cleanup */
  FREE_GSL_MATRIX(LU_decomp);
  FREE_GSL_PERMUTATION(LU_perm);
  FREE_GSL_MATRIX(inverse);

  return bound_box;

}

/**
 * Orthonormalise the columns of a matrix with respect to a metric (matrix is lower triangular)
 */
int XLALOrthonormaliseWRTMetric(
				gsl_matrix *matrix, /**< Matrix of columns to orthonormalise */
				gsl_matrix *metric  /**< Metric to orthonormalise with respect to */
				)
{
  
  const int n = metric->size1;

  int i = 0;
  int j = 0;
  gsl_vector *temp = NULL;
  double inner_prod = 0.0;

  /* Check input */
  if ((int)metric->size1 != n || (int)metric->size2 != n)
    XLAL_ERROR("'metric' is not square", XLAL_ESIZE);
  if (metric->size1 != matrix->size2 || metric->size2 != matrix->size2)
    XLAL_ERROR("'matrix' is not the same size as 'metric'", XLAL_ESIZE);

  /* Allocate */
  ALLOC_GSL_VECTOR(temp, n, XLAL_ERROR);

  /* Orthonormalise the columns of the matrix using numerically stabilised Gram-Schmidt */
  for (i = n - 1; i >= 0; --i) {
    gsl_vector_view col_i = gsl_matrix_column(matrix, i);
    for (j = n - 1; j > i; --j) {
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
  FREE_GSL_VECTOR(temp);

  return XLAL_SUCCESS;

}

/**
 * Transform a lattice generator to a square lower triangular form
 */
gsl_matrix *XLALSquareLowerTriangularLatticeGenerator(
						      gsl_matrix *generator /**< Generator matrix of lattice */
						      )
{

  const int m = generator->size1;
  const int n = generator->size2;
  int i = 0;
  int j = 0;
  gsl_matrix *QR_decomp = NULL;
  gsl_vector *QR_tau = NULL;
  gsl_matrix *Q = NULL;
  gsl_matrix *R = NULL;
  gsl_matrix *perm_sign = NULL;
  gsl_matrix *left = NULL;
  gsl_matrix *right = NULL;
  gsl_matrix *temp = NULL;
  gsl_matrix *result = NULL;
  double x = 0.0;

  /* Check input */
  if (m < n)
    XLAL_ERROR_NULL("'generator' must have number of rows >= number of columns", XLAL_ESIZE);
  
  /* Allocate memory */
  ALLOC_GSL_MATRIX(QR_decomp, m, n, XLAL_ERROR_NULL);
  ALLOC_GSL_VECTOR(QR_tau,       n, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(Q,         m, m, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(R,         m, n, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(perm_sign, n, m, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(left,      n, m, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(right,     n, n, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(temp,      m, n, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(result,    n, n, XLAL_ERROR_NULL);

  /* Find the QR decomposition of the generator */
  gsl_matrix_memcpy(QR_decomp, generator);
  gsl_linalg_QR_decomp(QR_decomp, QR_tau);
  gsl_linalg_QR_unpack(QR_decomp, QR_tau, Q, R);

  /* Build matrix to permute column order and make signs to diagonal positive */
  gsl_matrix_set_zero(perm_sign);
  for (i = 0; i < n; ++i) {
    for (j = 0; j < m; ++j) {
      if (i + j == n - 1) {
	x = gsl_matrix_get(R, j, j);
	gsl_matrix_set(perm_sign, i, j, x < 0 ? -1.0 : (x > 0 ? 1.0 : 0.0));
      }
    }
  }

  /* Calculate left side of transform (Q is transposed to get inverse) */
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, perm_sign, Q, 0.0, left);

  /* Build right side of transform */
  gsl_matrix_set_zero(right);
  for (i = 0; i < n; ++i)
    for (j = 0; j < n; ++j)
      if (i + j == n - 1)
	gsl_matrix_set(right, i, j, 1.0);

  /* Transform generator */
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, generator, right, 0.0, temp);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, left, temp, 0.0, result);

  /* Generator will be lower triangular, so zero out upper triangle */
  for (i = 0; i < n; ++i)
    for (j = i + 1; j < n; ++j)
      gsl_matrix_set(result, i, j, 0.0);

  /* Cleanup */
  FREE_GSL_MATRIX(QR_decomp);
  FREE_GSL_VECTOR(QR_tau);
  FREE_GSL_MATRIX(Q);
  FREE_GSL_MATRIX(R);
  FREE_GSL_MATRIX(perm_sign);
  FREE_GSL_MATRIX(left);
  FREE_GSL_MATRIX(right);
  FREE_GSL_MATRIX(temp);
  
  return result;

}

/**
 * Normalise a lattice generator matrix to have a specified covering radius
 */
int XLALNormaliseLatticeGenerator(
				  gsl_matrix *generator, /**< Generator matrix of lattice */
				  REAL8 norm_thickness,  /**< Normalised thickness of lattice */
				  REAL8 covering_radius  /**< Desired covering radius */
				  )
{
  
  const int m = generator->size1;
  const int n = generator->size2;
  gsl_matrix *LU_decomp = NULL;
  gsl_permutation *LU_perm = NULL;
  int LU_sign = 0;
  double generator_determinant = 0.0;
  double generator_covering_radius = 0.0;

  /* Allocate memory */
  ALLOC_GSL_MATRIX(LU_decomp, m, n, XLAL_ERROR);
  ALLOC_GSL_PERMUTATION(LU_perm, m, XLAL_ERROR);

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
  FREE_GSL_MATRIX(LU_decomp);
  FREE_GSL_PERMUTATION(LU_perm);

  return XLAL_SUCCESS;
  
}
