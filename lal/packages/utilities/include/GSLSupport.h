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
 * \brief Support routines for GSL
 */

#ifndef _GSLSUPPORT_H
#define _GSLSUPPORT_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <lal/LALRCSID.h>
#include <lal/LALDatatypes.h>

#ifdef __cplusplus
extern "C" {
#endif

  NRCSID(GSLSUPPORTH, "$Id$");

  /**
   * Functions
   */
  gsl_vector *XLALGSLVectorFromVAList(INT4, REAL8, ...);
  gsl_vector *XLALGSLVectorFromLALStringVector(LALStringVector*);
  gsl_matrix *XLALResizeGSLMatrix(gsl_matrix*, size_t, size_t, double);
  gsl_vector *XLALResizeGSLVector(gsl_vector*, size_t, double);
  gsl_vector_int *XLALResizeGSLVectorInt(gsl_vector_int*, size_t, int);

#ifdef __cplusplus
}
#endif

#endif

/**
 * Macros
 */

/* Allocate, print, and free 1D things */
#define ALLOC_GSL_1D(type, var, n, xlalerr)                                       \
  if ((var = gsl_ ## type ## _alloc(n)) == NULL)                                  \
    xlalerr("Could not allocate '" #var "'", XLAL_ENOMEM);
#define PRINT_GSL_1D(type, var, format)                                           \
{                                                                                 \
  size_t GSLSUPPORT_i;                                                            \
  printf("%s: ", #var);                                                           \
  if (var)                                                                        \
    for (GSLSUPPORT_i = 0; GSLSUPPORT_i < var->size; ++GSLSUPPORT_i)              \
      printf(format, gsl_ ## type ## _get(var, GSLSUPPORT_i));                    \
  else                                                                            \
    printf("NULL");                                                               \
  printf("\n");                                                                   \
}
#define FREE_GSL_1D(type, var)                                                    \
  if (var != NULL)                                                                \
    gsl_ ## type ## _free(var);

/* Allocate, print, and free 2D things */
#define ALLOC_GSL_2D(type, var, m, n, xlalerr)                                    \
  if ((var = gsl_ ## type ## _alloc(m, n)) == NULL)                               \
    xlalerr("Could not allocate '" #var "'", XLAL_ENOMEM);
#define PRINT_GSL_2D(type, var, format)                                           \
{                                                                                 \
  size_t GSLSUPPORT_i, GSLSUPPORT_j;                                              \
  printf("%s:\n", #var);                                                          \
  if (var)                                                                        \
    for (GSLSUPPORT_i = 0; GSLSUPPORT_i < var->size1; ++GSLSUPPORT_i) {           \
      for (GSLSUPPORT_j = 0; GSLSUPPORT_j < var->size2; ++GSLSUPPORT_j)           \
        printf(format, gsl_ ## type ## _get(var, GSLSUPPORT_i, GSLSUPPORT_j));    \
      printf("\n");                                                               \
    }                                                                             \
  else                                                                            \
    printf("NULL\n");                                                             \
}
#define FREE_GSL_2D(type, var)                                                    \
  if (var != NULL)                                                                \
    gsl_ ## type ## _free(var);

/* Allocate, print, and free gsl_vector */
#define ALLOC_GSL_VECTOR(var, n, xlalerr) ALLOC_GSL_1D(vector, var, n, xlalerr)
#define PRINT_GSL_VECTOR(var)             PRINT_GSL_1D(vector, var, "%g ")
#define FREE_GSL_VECTOR(var)              FREE_GSL_1D(vector, var)

/* Allocate, print, and free gsl_vector_int */
#define ALLOC_GSL_VECTOR_INT(var, n, xlalerr) ALLOC_GSL_1D(vector_int, var, n, xlalerr)
#define PRINT_GSL_VECTOR_INT(var)             PRINT_GSL_1D(vector_int, var, "%i ")
#define FREE_GSL_VECTOR_INT(var)              FREE_GSL_1D(vector_int, var)

/* Allocate, print, and free gsl_matrix */
#define ALLOC_GSL_MATRIX(var, m, n, xlalerr) ALLOC_GSL_2D(matrix, var, m, n, xlalerr)
#define PRINT_GSL_MATRIX(var)                PRINT_GSL_2D(matrix, var, "%g ")
#define FREE_GSL_MATRIX(var)                 FREE_GSL_2D(matrix, var)

/* Allocate, print, and free gsl_matrix_int */
#define ALLOC_GSL_MATRIX_INT(var, m, n, xlalerr) ALLOC_GSL_2D(matrix_int, var, m, n, xlalerr)
#define PRINT_GSL_MATRIX_INT(var)                PRINT_GSL_2D(matrix_int, var, "%i ")
#define FREE_GSL_MATRIX_INT(var)                 FREE_GSL_2D(matrix_int, var)

/* Allocate, print, and free gsl_combination */
#define ALLOC_GSL_COMBINATION(var, n, xlalerr) ALLOC_GSL_1D(combination, var, n, xlalerr)
#define PRINT_GSL_COMBINATION(var)             PRINT_GSL_1D(combination, var, "%i ")
#define FREE_GSL_COMBINATION(var)              FREE_GSL_1D(combination, var)

/* Allocate, print, and free gsl_permutation */
#define ALLOC_GSL_PERMUTATION(var, n, xlalerr) ALLOC_GSL_1D(permutation, var, n, xlalerr)
#define PRINT_GSL_PERMUTATION(var)             PRINT_GSL_1D(permutation, var, "%i ")
#define FREE_GSL_PERMUTATION(var)              FREE_GSL_1D(permutation, var)
