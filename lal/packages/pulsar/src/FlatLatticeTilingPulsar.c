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
#include "FlatLatticeTilingPulsar.h"
#include "FlatLatticeTilingSupport.h"

NRCSID(FLATLATTICETILINGPULSARC, "$Id$");

static const BOOLEAN TRUE  = (1==1);
static const BOOLEAN FALSE = (0==1);

/**
 * Setup a flat lattice tiling to use the spindown Fstat metric
 * (so no sky position tiling). Components are in the order
 * \f$\w_0,\alpha,\delta,\w_1,\w_2,...\f$
 * and will be converted on output to
 * \f$f_0,\alpha,\delta,f_1,f_2,...\f$
 */
int XLALSetFlatLatticeTilingSpindownFstatMetric(
						FlatLatticeTiling *tiling, /**< Tiling structure */
						REAL8 max_mismatch,        /**< Maximum mismatch */
						REAL8 Tspan                /**< Time span of the search */
						)
{
  
  const int n = tiling->dimension;
  
  int i, j;
  gsl_matrix *metric = NULL;
  gsl_vector *to_real;

  /* Check input */
  if (Tspan <= 0.0)
    XLAL_ERROR("Tspan must be strictly positive", XLAL_EINVAL);

  /* Allocate memory */
  ALLOC_GSL_MATRIX(metric, n, n, XLAL_ERROR);
  ALLOC_GSL_VECTOR(to_real, n, XLAL_ERROR);
  gsl_matrix_set_zero(metric);
  gsl_vector_set_all(to_real, 1.0);

  /* Calculate metric */
  for (i = 0; i < n - 2; ++i)
    for (j = 0; j < n - 2; ++j)
      gsl_matrix_set(metric, i + 2, j + 2, (1.0 * (i + 1) * (j + 1)) / ((i + 2) * (j + 2) * (i + j + 3)));

  /* Swap row/column 0 and row/column 2 to get right order */
  gsl_matrix_swap_rows(metric, 0, 2);
  gsl_matrix_swap_columns(metric, 0, 2);

  /* Calculate conversion factors for frequency and spindowns */
  {
    gsl_vector_view v = gsl_vector_subvector(to_real, 2, n - 2);
    double x = 1.0 / LAL_TWOPI;
    for (i = 0; i < n - 2; ++i) {
      x *= (i + 1) / Tspan;
      gsl_vector_set(&v.vector, i, x);
    }
  }
  
  /* Swap elements 0 and 2 to get the right order */
  gsl_vector_swap_elements(to_real, 0, 2);

  /* Set the metric of the flat lattice tiling */
  if (XLALSetFlatLatticeTilingMetric(tiling, metric, max_mismatch, to_real) != XLAL_SUCCESS)
    XLAL_ERROR("XLALSetFlatLatticeTilingMetric failed", XLAL_EFAILED);

  /* Cleanup */
  FREE_GSL_MATRIX(metric);
  FREE_GSL_VECTOR(to_real);

  return XLAL_SUCCESS;

}
