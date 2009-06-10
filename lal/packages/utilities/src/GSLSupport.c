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

#include <math.h>
#include <stdio.h>
#include <stdarg.h>

#include <gsl/gsl_combination.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_math.h>

#include <lal/LALRCSID.h>
#include <lal/LALStdlib.h>
#include <lal/LALMalloc.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/GSLSupport.h>

NRCSID(GSLSUPPORTC, "$Id$");

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
 * Resize a gsl_matrix and copy over existing contents
 */
gsl_matrix *XLALResizeGSLMatrix(
				gsl_matrix *m, /**< Matrix */
				size_t size1,  /**< New number of rows */
				size_t size2,  /**< New number of columns */
				double value   /**< Default value for new elements */
				)
{

  gsl_matrix *n = NULL;

  /* Allocate memory */
  ALLOC_GSL_MATRIX(n, size1, size2, XLAL_ERROR_NULL);
  gsl_matrix_set_all(n, value);

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
				size_t size,   /**< New number of elements */
				double value   /**< Default value for new elements */
				)
{

  gsl_vector *v = NULL;

  /* Allocate memory */
  ALLOC_GSL_VECTOR(v, size, XLAL_ERROR_NULL);
  gsl_vector_set_all(v, value);

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
 * Resize a gsl_vector_int and copy over existing contents
 */
gsl_vector_int *XLALResizeGSLVectorInt(
				       gsl_vector_int *u, /**< Vector */
				       size_t size,       /**< New number of elements */
				       int value          /**< Default value for new elements */
				       )
{

  gsl_vector_int *v = NULL;

  /* Allocate memory */
  ALLOC_GSL_VECTOR_INT(v, size, XLAL_ERROR_NULL);
  gsl_vector_int_set_all(v, value);

  /* Copy contents if any */
  if (u != NULL) {
    gsl_vector_int_view old = gsl_vector_int_subvector(u, 0, GSL_MIN(u->size, size));
    gsl_vector_int_view new = gsl_vector_int_subvector(v, 0, GSL_MIN(u->size, size));
    gsl_vector_int_memcpy(&new.vector, &old.vector);
    FREE_GSL_VECTOR_INT(u);
  }

  return v;

}
