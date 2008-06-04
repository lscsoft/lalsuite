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

#ifndef _FLATLATTICETILINGSUPPORT_H
#define _FLATLATTICETILINGSUPPORT_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <lal/LALRCSID.h>
#include <lal/LALDatatypes.h>

#ifdef __cplusplus
extern "C" {
#endif

  NRCSID(FLATLATTICETILINGSUPPORTH, "$Id$");

  /**
   * Structures
   */
  typedef struct tagTableau {

    /**< Number of constraints */
    int m;
    
    /**< Dimension of solution space */
    int n;
    
    /**< Tableau matrix */
    gsl_matrix *T;

    /**< Working copy of tableau matrix */
    gsl_matrix *work_T;
    
    /**< Basic variables */
    gsl_vector_int *basic;

  } SimplexMethodTableau;
  
  /**
   * Functions
   */
  gsl_vector *XLALGSLVectorFromVAList(INT4, REAL8, ...);
  gsl_vector *XLALGSLVectorFromLALStringVector(LALStringVector*);
  void XLALGaussianPivot(gsl_matrix*, size_t, size_t);
  void XLALReducedRowEchelonForm(gsl_matrix*);
  gsl_matrix *XLALResizeGSLMatrix(gsl_matrix*, size_t, size_t);
  gsl_vector *XLALResizeGSLVector(gsl_vector*, size_t);
  int XLALPolytopeVertices(gsl_matrix*, gsl_vector*, gsl_matrix**);
  int XLALPolytopeFromVertices(gsl_matrix*, gsl_matrix**, gsl_vector**);
  int XLALQuadraticSimplexMethodTableau(int, gsl_matrix*, gsl_vector*, gsl_matrix*, gsl_vector*, SimplexMethodTableau**);
  int XLALQuadraticSimplexMethod(SimplexMethodTableau*, gsl_vector*);
  void XLALFreeSimplexMethodTableau(SimplexMethodTableau*);
  gsl_vector *XLALMetricEllipseBoundingBox(gsl_matrix*, REAL8);
  int XLALOrthonormaliseWRTMetric(gsl_matrix*, gsl_matrix*);
  gsl_matrix *XLALSquareLowerTriangularLatticeGenerator(gsl_matrix*);
  int XLALNormaliseLatticeGenerator(gsl_matrix*, REAL8, REAL8);
  
#ifdef __cplusplus
}
#endif

#endif

/**
 * Macros
 */

/* Allocate and free gsl_matrix */
#define ALLOC_GSL_MATRIX(matrix, m, n, xlalerr)				\
  if ((matrix = gsl_matrix_alloc(m, n)) == NULL)			\
    xlalerr("Could not allocate '" #matrix "'", XLAL_ENOMEM);
#define FREE_GSL_MATRIX(matrix)						\
  if (matrix != NULL)							\
    gsl_matrix_free(matrix);

/* Allocate and free gsl_vector */
#define ALLOC_GSL_VECTOR(vector, n, xlalerr)				\
  if ((vector = gsl_vector_alloc(n)) == NULL)				\
    xlalerr("Could not allocate '" #vector "'", XLAL_ENOMEM);
#define FREE_GSL_VECTOR(vector)						\
  if (vector != NULL)							\
    gsl_vector_free(vector);

/* Allocate and free gsl_vector_int */
#define ALLOC_GSL_VECTOR_INT(vector_int, n, xlalerr)			\
  if ((vector_int = gsl_vector_int_alloc(n)) == NULL)			\
    xlalerr("Could not allocate '" #vector_int "'", XLAL_ENOMEM);
#define FREE_GSL_VECTOR_INT(vector_int)					\
  if (vector_int != NULL)						\
    gsl_vector_int_free(vector_int);

/* Allocate and free gsl_combination */
#define ALLOC_GSL_COMBINATION(combination, n, k, xlalerr)		\
  if ((combination = gsl_combination_alloc(n, k)) == NULL) 		\
    xlalerr("Could not allocate '" #combination "'", XLAL_ENOMEM);
#define FREE_GSL_COMBINATION(combination)				\
  if (combination != NULL) 						\
    gsl_combination_free(combination);

/* Allocate and free gsl_permutation */
#define ALLOC_GSL_PERMUTATION(permutation, n, xlalerr)			\
  if ((permutation = gsl_permutation_alloc(n)) == NULL)			\
    xlalerr("Could not allocate '" #permutation "'", XLAL_ENOMEM);
#define FREE_GSL_PERMUTATION(permutation)				\
  if (permutation != NULL)						\
    gsl_permutation_free(permutation);
