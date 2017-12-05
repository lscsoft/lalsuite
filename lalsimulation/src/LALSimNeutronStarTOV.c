/*
 * Copyright (C) 2013 J. Creighton, B. Lackey
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
 * @author Jolien Creighton, Benjamin Lackey
 * @addtogroup LALSimNeutronStarTOV_c
 * @brief Provides routines for solving the Tolman-Oppenheimer-Volkov equation.
 * @{
 */

#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALSimNeutronStar.h>

/** @cond */

/* Implements Eq. (50) of Damour & Nagar, Phys. Rev. D 80 084035 (2009).
 * See also Eq. (14) of Hinderer et al. Phys. Rev. D 81 123016 (2010). */
static double tidal_Love_number_k2(double c, double y)
{
    double num;
    double den;

    num = (8.0 / 5.0) * pow(1 - 2 * c, 2.0) * pow(c, 5)
        * (2 * c * (y - 1) - y + 2);
    den = 2 * c * (4 * (y + 1) * pow(c, 4) + (6 * y - 4) * pow(c, 3)
        + (26 - 22 * y) * c * c + 3 * (5 * y - 8) * c - 3 * y + 6);
    den -= 3 * pow(1 - 2 * c, 2) * (2 * c * (y - 1) - y + 2)
        * log(1.0 / (1 - 2 * c));

    return num / den;
}

/* For convenience, use a structure that provides a dictionary between
 * a vector of the ode variables and what they represent. */
struct tov_ode_vars {
    double r;   /* radial coordinate, m */
    double m;   /* mass within r in geometric units, m */
    double H;   /* stellar perturbation in arbitrary units */
    double b;   /* derivative of metric pertubation in arbitrary units */
};

#define TOV_ODE_VARS_DIM (sizeof(struct tov_ode_vars)/sizeof(double))

/* Casts an array of doubles to a structure with named parameter. */
static struct tov_ode_vars *tov_ode_vars_cast(const double *y)
{
    union {
        const double *y;
        struct tov_ode_vars *v;
    } u = {
    y};
    return u.v;
}

/* ODE integrand for TOV equations with pseudo-enthalpy independent variable.
 * Implements Eqs. (5) and (6) of Lindblom, Astrophys. J. 398, 569 (1992).
 * Also uses Eqs. (7) and (8) [ibid] for inner boundary data, and
 * Eqs. (18), (27), (28) of Damour & Nagar, Phys. Rev. D 80, 084035 (2009)
 * [See also: Eqs. (11) & (12) Hinderer et al. Phys. Rev. D 81 123016 (2010)]
 * for the metric perturbation used to obtain the Love number. */
static int tov_ode(double h, const double *y, double *dy, void *params)
{
    struct tov_ode_vars *vars = tov_ode_vars_cast(y);
    struct tov_ode_vars *derivs = tov_ode_vars_cast(dy);
    LALSimNeutronStarEOS *eos = params;

    double r = vars->r;
    double m = vars->m;
    double H = vars->H;
    double b = vars->b;
    double p =
        XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometerized(h, eos);
    double e =
        XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometerized(h,
        eos);
    double dedp =
        XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometerized(p, eos);
    /* Eq. (18) of Damour & Nagar PRD 80 084035 (2009). */
    double A = 1.0 / (1.0 - 2.0 * m / r);
    /* Eq. (28) of Damour & Nagar PRD 80 084035 (2009). */
    double C1 = 2.0 / r + A * (2.0 * m / (r * r) + 4.0 * LAL_PI * r * (p - e));
    /* Eq. (29) of Damour & Nagar PRD 80 084035 (2009). */
    double C0 =
        A * (-(2) * (2 + 1) / (r * r) + 4.0 * LAL_PI * (e + p) * dedp +
        4.0 * LAL_PI * (5.0 * e + 9.0 * p)) - pow(2.0 * (m +
            4.0 * LAL_PI * r * r * r * p) / (r * (r - 2.0 * m)), 2.0);
    double dr = -r * (r - 2.0 * m) / (m + 4.0 * LAL_PI * r * r * r * p);
    double dm = 4.0 * LAL_PI * r * r * e * dr;
    double dH = b * dr;
    double db = -(C0 * H + C1 * b) * dr;

    derivs->r = dr;
    derivs->m = dm;
    derivs->H = dH;
    derivs->b = db;
    return 0;
}

/** @endcond */

/**
 * @brief Integrates the Tolman-Oppenheimer-Volkov stellar structure equations.
 * @details
 * Solves the Tolman-Oppenheimer-Volkov stellar structure equations using the
 * pseudo-enthalpy formalism introduced in:
 * Lindblom (1992) "Determining the Nuclear Equation of State from Neutron-Star
 * Masses and Radii", Astrophys. J. 398 569.
 * @param[out] radius The radius of the star in m.
 * @param[out] mass The mass of the star in kg.
 * @param[out] love_number_k2 The k_2 tidal love number of the star.
 * @param[in] central_pressure_si The central pressure of the star in Pa.
 * @param eos Pointer to the Equation of State structure.
 * @retval 0 Success.
 * @retval <0 Failure.
 */
int XLALSimNeutronStarTOVODEIntegrate(double *radius, double *mass,
    double *love_number_k2, double central_pressure_si,
    LALSimNeutronStarEOS * eos)
{
    /* ode integration variables */
    const double epsabs = 0.0, epsrel = 1e-6;
    double y[TOV_ODE_VARS_DIM];
    double dy[TOV_ODE_VARS_DIM];
    struct tov_ode_vars *vars = tov_ode_vars_cast(y);
    gsl_odeiv_system sys = { tov_ode, NULL, TOV_ODE_VARS_DIM, eos };
    gsl_odeiv_step *step =
        gsl_odeiv_step_alloc(gsl_odeiv_step_rk8pd, TOV_ODE_VARS_DIM);
    gsl_odeiv_control *ctrl = gsl_odeiv_control_y_new(epsabs, epsrel);
    gsl_odeiv_evolve *evolv = gsl_odeiv_evolve_alloc(TOV_ODE_VARS_DIM);

    /* central values */
    /* note: will be updated with Lindblom's series expansion */
    /* geometrisized units for variables in length (m) */
    double pc = central_pressure_si * LAL_G_C4_SI;
    double ec =
        XLALSimNeutronStarEOSEnergyDensityOfPressureGeometerized(pc, eos);
    double hc =
        XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometerized(pc, eos);
    double dedp_c =
        XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometerized(pc, eos);
    double dhdp_c = 1.0 / (ec + pc);
    double dedh_c = dedp_c / dhdp_c;
    double dh = -1e-3 * hc;
    double h0 = hc + dh;
    double h1 = 0.0 - dh;
    double r0 = sqrt(-3.0 * dh / (2.0 * LAL_PI * (ec + 3.0 * pc)));
    double m0 = 4.0 * LAL_PI * r0 * r0 * r0 * ec / 3.0;
    double H0 = r0 * r0;
    double b0 = 2.0 * r0;

    double yy;
    double c;
    double h;
    size_t i;

    /* series expansion for the initial core */

    /* second factor of Eq. (7) of Lindblom (1992) */
    r0 *= 1.0 + 0.25 * dh * (ec - 3.0 * pc  - 0.6 * dedh_c) / (ec + 3.0 * pc);
    /* second factor of Eq. (8) of Lindblom (1992) */
    m0 *= 1.0 + 0.6 * dh * dedh_c / ec;

    /* perform integration */
    vars->r = r0;
    vars->m = m0;
    vars->H = H0;
    vars->b = b0;
    h = h0;
    while (h > h1) {
        int s =
            gsl_odeiv_evolve_apply(evolv, ctrl, step, &sys, &h, h1, &dh, y);
        if (s != GSL_SUCCESS)
            XLAL_ERROR(XLAL_EERR,
                "Error encountered in GSL's ODE integrator\n");
    }

    /* take one final Euler step to get to surface */
    tov_ode(h, y, dy, eos);
    for (i = 0; i < TOV_ODE_VARS_DIM; ++i)
        y[i] += dy[i] * (0.0 - h1);

    /* compute tidal Love number k2 */
    c = vars->m / vars->r;      /* compactness */
    yy = vars->r * vars->b / vars->H;

    /* convert from geometric units to SI units */
    *radius = vars->r;
    *mass = vars->m * LAL_MSUN_SI / LAL_MRSUN_SI;
    *love_number_k2 = tidal_Love_number_k2(c, yy);

    /* free ode memory */
    gsl_odeiv_evolve_free(evolv);
    gsl_odeiv_control_free(ctrl);
    gsl_odeiv_step_free(step);
    return 0;
}

/** @} */
