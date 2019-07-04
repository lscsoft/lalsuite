/*
 * Copyright (C) 2019 Andrew Matas, Jonathan Thompson, Edward Fauchon-Jones, Sebastian Khan
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

#include <math.h>
#include <lal/LALSimIMR.h>
#include <lal/Units.h>
#include <lal/LALConstants.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_poly.h>

#ifndef _OPENMP
#define omp ignore
#endif

/**
 * @author Andrew Matas, Jonathan Thompson, Edward Fauchon-Jones, Sebastian Khan
 * @addtogroup LALSimNSBHProperties_c Module LALSimNSBHProperties.c
 * @ingroup lalsimulation_general
 * @brief Provides routines for NSBH waveform models.
 * @{
 *
 * @name Kerr routines
 * @{
 */

/**
 * GW frequency for a particle on Kerr
 *
 * GW frequency (in units of inverse  M_BH) for a particle at a given separation
 * r from a Kerr BH with mass M and dimensionless spin parameter a.
 */
double XLALSimNSBH_fGWinKerr(
    const REAL8 r,   /**< Separtion */
    const REAL8 M,   /**< Kerr BH mass */
    const REAL8 a    /**< Dimensionless spin parameter */
) {
  return 1.0/(LAL_PI*(a*M + sqrt(pow(r, 3.0)/M)));
}

/**
 * Kerr BH ISCO radius
 *
 * Kerr BH ISCO radius for a BH with unit mass as a function of its
 * dimensionless spin parameter.
 */
double XLALSimNSBH_rKerrISCO(
    const REAL8 a    /**< Dimensionless spin parameter */
) {
  REAL8 Z1 = 1.0 + pow(1.0 - pow(a, 2.0), 1.0/3.0)*(pow(1.0+a, 1.0/3.0) + pow(1.0-a, 1.0/3.0));
  REAL8 Z2 = sqrt(3.0*pow(a,2.0) + pow(Z1, 2.0));
  if (a > 0.0) {
    return 3.0 + Z2 - sqrt((3.0 - Z1)*(3.0 + Z1 + 2.0*Z2));
  } else {
    return 3.0 + Z2 + sqrt((3.0 - Z1)*(3.0 + Z1 + 2.0*Z2));
  }
}

/**
 * @}
 *
 * @name NSBH routines
 * @{
 */

/**
 * Relativistic correction to orbital radius at mass-shedding
 *
 * Relativistic correction to the standard Newtonian estimate of the orbital
 * radius at mass-shedding. See Eq. (8) in https://arxiv.org/abs/1509.00512 for
 * a description of the implicit equation for `xi_tide`.
 */
double XLALSimNSBH_xi_tide(
    const REAL8 q,     /**< Mass ratio of the NSBH system M_BH/M_NS > 1 */
    const REAL8 a,     /**< The BH dimensionless spin parameter */
    const REAL8 mu     /**< The ratio of the BH and the NS radius M_BH/R_NS = q*C where C is the compactness C = M_NS/R_NS. */
) {
  REAL8 p[11] = {
    -3.0*q*pow(mu,2.0)*pow(a,2.0),
    0.0,
    6*q*mu,
    0.0,
    -3.0*q,
    0.0,
    0.0,
    2.0*a*pow(mu,1.5),
    -3*mu,
    0.0,
    1.0};
  double z[20];

  gsl_error_handler_t * error_handler = gsl_set_error_handler_off();
  gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc(11);
  if (!w)
    XLAL_ERROR(XLAL_EFAILED, "Cannot setup workspace to solve for xi_tide");
  int status = gsl_poly_complex_solve(p, 11, w, z);
  gsl_poly_complex_workspace_free(w);
  gsl_set_error_handler(error_handler);
  if (status)
    XLAL_ERROR(XLAL_EFAILED, "Cannot find solution for xi_tide");

  double root = 0.0;
  for (int i=0; i<10; i++) {
    if (fabs(z[2*i+1]) < 1.0E-5) {
      if (z[2*i] > 0.0 && pow(z[2*i],2.0) > root) {
        root = pow(z[2*i],2.0);
      }
    }
  }

  return root;
}

/**
 * Compactness as a function of tidal deformability
 *
 * For Lambda > 1 this computes the NS compactness from the tidal
 * deformability as described by the universal relation from N. Yagi
 * and N. Yunes, Phys. Rep. 681 (2017), Eq. (78).
 *
 * For Lambda <= 1 this extrapolates to the BH limit such that the
 * interpolation over [0, 1] is C^1 continuous as Lambda=1.
 */
double XLALSimNSBH_compactness_from_lambda(
    const REAL8 Lambda    /**< Dimensionless tidal deformability */
) {
  // Fitting coefficients
  const REAL8 a0 = 0.360;
  const REAL8 a1 = -0.0355;
  const REAL8 a2 = 0.000705;

  if (Lambda > 1) {
    const REAL8 log_lambda = log(Lambda);
    return a0 + a1*log_lambda + a2*pow(log_lambda,2.0);
  } else {
    const REAL8 L2 = Lambda*Lambda;
    const REAL8 L3 = L2*Lambda;
    return 0.5 + (3*a0-a1-1.5)*L2 + (-2*a0+a1+1)*L3;
  }
}

/**
 * Baryonic mass of the torus remnant of a BH-NS merger
 *
 * Baryonic mass of the torus remnant of a BH-NS merger in units of the NS
 * baryonic mass. See arXiv:1207.6304v1 for details on the "Newtonian" and the
 * "relativistic" fits.
 *
 * See Eq. (11) in https://arxiv.org/abs/1509.00512 for the definition of this
 * quantity.
 */
double XLALSimNSBH_torus_mass_fit(
  const REAL8 q,   /**< The binary mass ratio (assumes q>1) */
  const REAL8 a,   /**< The BH dimensionless spin parameter */
  const REAL8 C    /**< Neutron star compactness. */
) {

  REAL8 mu = q*C;
  REAL8 alpha = 0.296;
  REAL8 beta = 0.171;
  REAL8 xi = XLALSimNSBH_xi_tide(q, a, mu);
  REAL8 Mtorus = alpha * xi * (1.0-2.0*C) - beta * mu * XLALSimNSBH_rKerrISCO(a);

  if (Mtorus <0.0)
    return 0.0;
  else
    return Mtorus;
}

/** @} */

/** @} */
