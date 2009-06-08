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
 * \brief Pulsar-specific routines for FlatLatticeTiling
 */

#include <math.h>
#include <stdio.h>
#include <string.h>

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
#include <lal/GSLSupport.h>
#include <lal/FlatLatticeTilingPulsar.h>

NRCSID(FLATLATTICETILINGPULSARC, "$Id$");

#define TRUE  (1==1)
#define FALSE (1==0)

/**
 * Calculate the factorial of an integer
 */
static INT4 factorial(
		      INT4 i
		      )
{
  
  INT4 f = 1;
  
  while (i > 1)
    f *= i--;
  
  return f;
  
}

/**
 * Set a flat lattice tiling to the spindown Fstat metric
 * (so no sky position tiling). Components are in the order
 * \f$\omega_0,\alpha,\delta,\omega_1,\omega_2,...\f$
 * and will be converted on output to
 * \f$f_0,\alpha,\delta,f_1,f_2,...\f$
 * using the conversions
 * \f$f_k = \omega_k {(k+1)! \over {2\pi T^{k+1}}}f$
 */
int XLALSetFlatLatticeTilingSpindownFstatMetric(
						FlatLatticeTiling *tiling, /**< Tiling structure */
						REAL8 max_mismatch,        /**< Maximum mismatch */
						REAL8 Tspan                /**< Time span of the search */
						)
{

  const int n = tiling->dimensions;
  
  int i, j;
  gsl_matrix *norm_metric = NULL;
  gsl_vector *norm_to_real;

  /* Check input */
  if (Tspan <= 0.0)
    XLAL_ERROR("Tspan must be strictly positive", XLAL_EINVAL);

  /* Allocate memory */
  ALLOC_GSL_MATRIX(norm_metric, n, n, XLAL_ERROR);
  ALLOC_GSL_VECTOR(norm_to_real,   n, XLAL_ERROR);
  gsl_matrix_set_identity(norm_metric);
  gsl_vector_set_all(norm_to_real, 1.0);

  /* Calculate metric and conversion factors */
  for (i = 0; i < n - 2; ++i) {

    for (j = 0; j < n - 2; ++j)
      gsl_matrix_set(norm_metric, i + 2, j + 2, (1.0 * (i + 1) * (j + 1)) / ((i + 2) * (j + 2) * (i + j + 3)));

    gsl_vector_set(norm_to_real, i + 2, 1.0 * factorial(i + 1) / (LAL_TWOPI * pow(Tspan, i + 1)));

  }

  /* Swap rows/columns 0 and 2 to get right order */
  gsl_matrix_swap_rows(norm_metric, 0, 2);
  gsl_matrix_swap_columns(norm_metric, 0, 2);
  gsl_vector_swap_elements(norm_to_real, 0, 2);

  /* Set the metric of the flat lattice tiling */
  if (XLALSetFlatLatticeTilingMetric(tiling, norm_metric, max_mismatch, norm_to_real) != XLAL_SUCCESS)
    XLAL_ERROR("XLALSetFlatLatticeTilingMetric failed", XLAL_EFAILED);

  /* Cleanup */
  FREE_GSL_MATRIX(norm_metric);
  FREE_GSL_VECTOR(norm_to_real);

  return XLAL_SUCCESS;

}

/**
 * Set a flat lattice tiling to a parameter space defined by
 * the age and possible braking index range of an object
 */
static BOOLEAN AgeBrakingIndexBound(void *data, INT4 dimension, gsl_vector *point, REAL8 *lower, REAL8 *upper)
{
  const char *fn = "AgeBrakingIndexBound()";
  double x;

  /* Set constant based on dimension */
  switch (dimension) {
  case 0:
    x = 1.0;
    break;
  case 1:
    x = gsl_vector_get(point, 0);
    break;
  case 2:
    x = gsl_vector_get(point, 1);
    x *= x;
    x /= gsl_vector_get(point, 0);
    break;
  default:
    XLALPrintError ("%s: invalid dimension %d input, allowed are 0-2.\n", fn, dimension );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  /* Set lower and upper bound */
  *lower = x * gsl_matrix_get((gsl_matrix*)data, dimension, 0);
  *upper = x * gsl_matrix_get((gsl_matrix*)data, dimension, 1);
  
  return TRUE;
  
}
static void AgeBrakingIndexFree(void *data)
{

  /* Cleanup */
  FREE_GSL_MATRIX((gsl_matrix*)data);

}
int XLALAddFlatLatticeTilingAgeBrakingIndexBounds(
						  FlatLatticeTiling *tiling, /**< Tiling structure */
						  REAL8 freq,                /**< Starting frequency */
						  REAL8 freq_band,           /**< Frequency band */
						  REAL8 age,                 /**< Spindown age */
						  REAL8 min_braking,         /**< Minimum braking index */
						  REAL8 max_braking,         /**< Maximum braking index */
						  INT4 offset,               /**< Number of dimensions offset between first dimension and frequency */
						  INT4 gap                   /**< Number of dimensions gap netween frequency and first spindown */
						  )
{
  
  gsl_matrix *data = NULL;
  UINT8 bound;

  /* Check tiling dimension */
  if (tiling->dimensions < (3 + gap))
    XLAL_ERROR("'tiling->dimensions' is too small", XLAL_EINVAL);
  
  /* Allocate memory */
  ALLOC_GSL_MATRIX(data, tiling->dimensions, 2, XLAL_ERROR);
  gsl_matrix_set_zero(data);

  /* Set frequency bounds */
  gsl_matrix_set(data, 0, 0, freq);
  gsl_matrix_set(data, 0, 1, freq + freq_band);

  /* Set first spindown bounds */
  gsl_matrix_set(data, 1, 0, -1.0 / ((min_braking - 1.0) * age));
  gsl_matrix_set(data, 1, 1, -1.0 / ((max_braking - 1.0) * age));

  /* Set second spindown bounds */
  gsl_matrix_set(data, 2, 0, min_braking);
  gsl_matrix_set(data, 2, 1, max_braking);

  /* Set bound dimensions */
  bound = ((UINT8)1) | (((UINT8)6) << gap) << offset;

  /* Set parameter space */
  if (XLAL_SUCCESS != XLALAddFlatLatticeTilingBound(tiling, bound, AgeBrakingIndexBound,
						    (void*)data, AgeBrakingIndexFree))
    XLAL_ERROR("XLALAddFlatLatticeTilingBound failed", XLAL_EFAILED);
  
  return XLAL_SUCCESS;
  
}
