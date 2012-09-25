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

#include <lal/FlatLatticeTilingPulsar.h>
#include <lal/LALStdlib.h>
#include <lal/LALMalloc.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/Factorial.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/**
 * Set a flat lattice tiling to the spindown Fstat metric
 * (so no sky position tiling). Components are in the order
 * \f$\alpha,\delta,\omega_0,\omega_1,\omega_2,...\f$
 * and will be converted on output to
 * \f$\alpha,\delta,f_0,f_1,f_2,...\f$
 * using the conversions
 * \f$f_k = \omega_k \frac{(k+1)!}{{2\pi T^{k+1}}}\f$
 */
int XLALSetFlatLatticeTilingSpindownFstatMetric(
  FlatLatticeTiling *tiling, /**< Tiling structure */
  double max_mismatch,        /**< Maximum mismatch */
  double Tspan                /**< Time span of the search */
  )
{

  const size_t n = XLALGetFlatLatticeDimensions(tiling);

  /* Check input */
  XLAL_CHECK(Tspan > 0.0, XLAL_EINVAL);

  /* Allocate memory */
  gsl_matrix* norm_metric = gsl_matrix_alloc(n, n);
  XLAL_CHECK(norm_metric != NULL, XLAL_ENOMEM);
  gsl_vector* norm_to_real = gsl_vector_alloc(n);
  XLAL_CHECK(norm_to_real != NULL, XLAL_ENOMEM);
  gsl_matrix_set_identity(norm_metric);

  /* Calculate metric and conversion factors */
  for (size_t i = 0; i < n - 2; ++i) {
    for (size_t j = i; j < n - 2; ++j) {
      gsl_matrix_set(norm_metric, i + 2, j + 2, (
                       4.0 * LAL_PI * LAL_PI * pow(Tspan, i + j + 2) * (i + 1) * (j + 1)
                       ) / (
                         LAL_FACT[i + 1] * LAL_FACT[j + 1] * (i + 2) * (j + 2) * (i + j + 3)
                         ));
      gsl_matrix_set(norm_metric, j + 2, i + 2, gsl_matrix_get(norm_metric, i + 2, j + 2));
    }
  }

  /* Set the metric of the flat lattice tiling */
  XLAL_CHECK(XLALSetFlatLatticeMetric(tiling, norm_metric, max_mismatch) == XLAL_SUCCESS, XLAL_EFAILED);

  /* Cleanup */
  FREE_GSL_MATRIX(norm_metric);
  FREE_GSL_VECTOR(norm_to_real);

  return XLAL_SUCCESS;

}

/**
 * Set a flat lattice tiling to a parameter space defined by
 * the age and possible braking index range of an object
 */
static void AgeBraking1stSpindownBound(
  const gsl_vector_uint* bound UNUSED,
  const gsl_vector* point,
  const void* data,
  gsl_vector* lower,
  gsl_vector* upper
  )
{

  /* Set lower and upper bound */
  double x = gsl_vector_get(point, point->size - 1);
  gsl_vector_set(lower, 0, x * ((const double*)data)[0]);
  gsl_vector_set(upper, 0, x * ((const double*)data)[1]);

}
static void AgeBraking2ndSpindownBound(
  const gsl_vector_uint* bound UNUSED,
  const gsl_vector* point,
  const void* data,
  gsl_vector* lower,
  gsl_vector* upper
  )
{

  /* Set lower and upper bound */
  double x = gsl_vector_get(point, point->size - 1);
  x *= x;
  x /= gsl_vector_get(point, point->size - 2);
  gsl_vector_set(lower, 0, x * ((const double*)data)[0]);
  gsl_vector_set(upper, 0, x * ((const double*)data)[1]);

}
int XLALSetFlatLatticeTilingAgeBrakingIndexBounds(
  FlatLatticeTiling *tiling, /**< Tiling structure */
  double freq,                /**< Starting frequency */
  double freq_band,           /**< Frequency band */
  double age,                 /**< Spindown age */
  double min_braking,         /**< Minimum braking index */
  double max_braking          /**< Maximum braking index */
  )
{

  // Allocate memory
  double* f1dot_data = XLALCalloc(2, sizeof(double));
  XLAL_CHECK(f1dot_data != NULL, XLAL_ENOMEM);
  double* f2dot_data = XLALCalloc(2, sizeof(double));
  XLAL_CHECK(f2dot_data != NULL, XLAL_ENOMEM);

  // Set frequency bounds
  XLAL_CHECK(XLALSetFlatLatticeConstantBound(tiling, 2, freq, freq + freq_band) == XLAL_SUCCESS, XLAL_EFAILED);

  // Set first spindown bounds
  f1dot_data[0] = -1.0 / ((min_braking - 1.0) * age);
  f1dot_data[1] = -1.0 / ((max_braking - 1.0) * age);
  XLAL_CHECK(XLALSetFlatLatticeBound(tiling, 3, false, AgeBraking1stSpindownBound, (void*)f1dot_data) == XLAL_SUCCESS, XLAL_EFAILED);

  // Set second spindown bounds
  f2dot_data[0] = min_braking;
  f2dot_data[1] = max_braking;
  XLAL_CHECK(XLALSetFlatLatticeBound(tiling, 4, false, AgeBraking2ndSpindownBound, (void*)f2dot_data) == XLAL_SUCCESS, XLAL_EFAILED);

  return XLAL_SUCCESS;

}
