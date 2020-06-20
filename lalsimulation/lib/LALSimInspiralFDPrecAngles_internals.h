#ifndef _LALSIM_INS_FD_PREC_ANGLES_INTERNALS
#define _LALSIM_INS_FD_PREC_ANGLES_INTERNALS

/*
 * Copyright (C) 2017 Katerina Chatziioannou, Sebastian Khan
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
#include <stdio.h>
#include <string.h>
#include "gsl/gsl_sf_elljac.h"
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_mode.h>
#include <lal/LALError.h>

// #include "LALSimInspiralFDPrecAngles.h"
#include <lal/LALSimInspiralFDPrecAngles.h>

static int InitializeSystem(sysq* system, /**< [out] Pointer to sysq struct  */
                             const double m1,  /**< Primary mass in SI (kg) */
                             const double m2,  /**< Secondary mass in SI (kg) */
                             const double mul, /**< Cosine of Polar angle of orbital angular momentum */
                             const double phl, /**< Azimuthal angle of orbital angular momentum  */
                             const double mu1, /**< Cosine of Polar angle of primary spin w.r.t. orbital angular momentum */
                             const double ph1, /**< Azimuthal angle of primary spin  */
                             const double ch1, /**< Dimensionless spin magnitude of primary spin */
                             const double mu2, /**< Cosine of Polar angle of secondary spin w.r.t. orbital angular momentum */
                             const double ph2, /**< Azimuthal angle of secondary spin  */
                             const double ch2, /**< Dimensionless spin magnitude of secondary spin */
                             const double f_0, /**< Reference Gravitational Wave frequency (Hz) */
                             const int ExpansionOrder /**< Keep terms upto ExpansionOrder in precession angles phi_z and zeta */
                         );

static double DotProd(const vector vec1, const vector vec2);
static double Norm(const vector vec1);
static vector CreateSphere(const double r, const double th, const double ph);
static vector ScalarProd(const double c, const vector vec);
static vector Sum(const vector vec1, const vector vec2);
static vector CrossProd(const vector vec1, const vector vec2);

static vector Roots(const double L_norm, const double J_norm, const sysq *system);
static vector BCDcoeff(const double L_norm, const double J_norm, const sysq *system);

static double beta(const double a, const double b, const sysq *system);
static double sigma(const double a, const double b, const sysq *system);
static double tau(const double a, const double b, const sysq *system);

static double J_norm_of_xi(const double L_norm, const sysq *system);
static double S_norm_of_xi(const double xi, const double xi_2, const vector roots, const sysq *system);
static double L_norm_2PN_NonSpinning_of_xi(const double xi_2, const double L_norm, const sysq *system);
static double L_norm_3PN_of_xi(const double xi, const double xi_2, const double L_norm, const sysq *system);

static vector c(const double xi, const double xi_2, const double J_norm, const vector roots, const sysq *system);
static vector d(const double L_norm, const double J_norm, const vector roots);

static vector compute_phiz_zeta_costhetaL2PNNonSpinning(const double xi, const sysq *system);
static vector compute_phiz_zeta_costhetaL3PN(const double xi, const sysq *system);
static vector compute_phiz_zeta_costhetaL(const double xi, const sysq *system);

static double costhetaL(const double J_norm, const double L_norm, const double S_norm);

static double u_of_xi(const double xi, const double xi_2, const sysq *system);
static double psidot(const double xi, const double xi_2, const vector roots, const sysq *system);

static vector computeMScorrections (const double xi, const double xi_2, const double L_norm, const double J_norm, const vector roots, const sysq *system);
static double phiz_of_xi(const double xi, const double xi_2, const double J_norm, const sysq *system);
static double zeta_of_xi(const double xi, const double xi_2, const sysq *system);

#endif	// of #ifndef _LALSIM_INS_FD_PREC_ANGLES_INTERNALS
