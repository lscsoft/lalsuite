/*
 * Copyright (C) 2013 J. Creighton, B. Lackey
 * Copyright (C) 2026 P. Davis, M. Oertel, L. Suleiman, L. Wade
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */
/**
 * @author Jolien Creighton, Benjamin Lackey, Philip Davis, Micaela Oertel, Lami Suleiman, Leslie Wade
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


/*
 * Implements Eq.(50) of Damour & Nagar, Phys. Rev. D 80 084035 (2009).
 * See also Eq.(14) of Hinderer et al. Phys. Rev. D 81 123016 (2010).
 */
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

/*
 * Implements Eq.(51) of Damour & Nagar, Phys. Rev. D 80 084035 (2009).
 */
static double tidal_Love_number_k3(double c, double y)
{
    double num;
    double den;
    double c2 = c*c;
    double c3 = c2*c;
    double c4 = c2*c2;
    double c5 = c4*c;
    double c7 = c4*c3;

    num = (8.0 / 7.0) * pow(1 - 2 * c, 2.0) * c7
        * (2 * (y - 1) * c2 - 3.0 * (y - 2.0) * c + y - 3.0);
    den = 2.0 * c * (  4*(y + 1)*c5 + 2*(9*y - 2)*c4
                        - 20*(7*y - 9)*c3 + 5*(37*y - 72)*c2
                        - 45*(2*y - 5)*c  + 15*(y - 3) ) ;

    den -= 15*pow(1 - 2 * c, 2.0) * ( 2*(y - 1)*c2 - 3*(y-2)*c + y - 3 ) * log(1.0 / (1 - 2 * c));

    return num / den;
}

/*
 * Implements Eq.(52) of Damour & Nagar, Phys. Rev. D 80 084035 (2009).
 */
static double tidal_Love_number_k4(double c, double y)
{
    double num;
    double den;
    double c2 = c*c;
    double c3 = c2*c;
    double c4 = c2*c2;
    double c5 = c4*c;
    double c6 = c3*c3;
    double c9 = c6*c3;

    num = (32.0 / 147.0) * pow(1 - 2 * c, 2.0) * c9
        * ( 12*(y - 1)*c3 - 34*(y - 2)*c2 + 28*(y - 3)*c - 7*(y - 4) );
    den = 2.0 * c * (  8*(y + 1)*c6 + (68*y - 8)*c5 + (1284 - 996*y)*c4
                        + 40*(55*y - 116)*c3 + (5360 - 1910*y)*c2
                        + 105*(7*y - 24)*c - 105*(y - 4) ) ;

    den -= 15*pow(1 - 2 * c, 2.0) * ( 12*(y - 1)*c3 - 34*(y - 2)*c2 + 28*(y - 3)*c - 7*(y - 4) ) * log(1.0 / (1 - 2 * c));

    return num / den;
}

/*
 * Implements Eq.(29) of Damour & Nagar PRD 80 084035 (2009).
 */
static double constant_zero_for_kl(double r, double m, double A, double e, double p, double dedp, int l){
    double C0 =
    A * (-(l) * (l + 1) / (r * r) + 4.0 * LAL_PI * (e + p) * dedp +
        4.0 * LAL_PI * (5.0 * e + 9.0 * p)) - pow(2.0 * (m +
            4.0 * LAL_PI * r * r * r * p) / (r * (r - 2.0 * m)), 2.0);

    return C0;
}


/*
 * For convenience, use a structure that provides a dictionary between
 * a vector of the ode variables and what they represent.
 */
struct tov_ode_vars {
    double r;   /* radial coordinate, m */
    double m;   /* mass within r in geometric units, m */
    double H;   /* stellar perturbation in arbitrary units */
    double b;   /* derivative of metric pertubation in arbitrary units */
};

#define TOV_ODE_VARS_DIM (sizeof(struct tov_ode_vars)/sizeof(double))

/*
 * Casts an array of doubles to a structure with named parameter.
 */
static struct tov_ode_vars *tov_ode_vars_cast(const double *y)
{
    union {
        const double *y;
        struct tov_ode_vars *v;
    } u = {
    y};
    return u.v;
}

/*
 * ODE integrand for TOV equations with pseudo-enthalpy independent variable.
 * Implements Eqs. (5) and (6) of Lindblom, Astrophys. J. 398, 569 (1992).
 * Also uses Eqs. (7) and (8) [ibid] for inner boundary data, and
 * Eqs. (18), (27), (28) of Damour & Nagar, Phys. Rev. D 80, 084035 (2009)
 * [See also: Eqs. (11) & (12) Hinderer et al. Phys. Rev. D 81 123016 (2010)]
 * for the metric perturbation used to obtain the Love number.
 */
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
        XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometrized(h, eos);
    double e =
        XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrized(h,
        eos);
    double dedp =
        XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrized(p, eos);
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

/**
 * @brief Integrates the stellar structure system of ordinary differential
 * equations to obtain neutron star's macroscopic parameters.
 * @details
 * Solves the Tolman-Oppenheimer-Volkoff stellar structure equations and
 * the gravito-electric tidal Love number equations, given an equation of
 * state and a central pressure. The solution to this set of ODEs provides
 * the neutron star's radius, gravitational mass and second order gravito-electric
 * tidal Love number.
 * It uses the pseudo-enthalpy formalism introduced in Lindblom (1992)
 * Astrophys. J. 398 569 "Determining the Nuclear Equation of State from
 * Neutron-Star Masses and Radii".
 * The surface of the star is defined by a minimum value of the pseudo-enthalpy
 * h1 = 1e-12 * hc (with hc the central pseudo-enthalpy) ; if the minimum
 * pseudo-enthalpy of the equation of state is larger, a polytrope
 * (degenerate non-relativistic gas) is used to extend the model beyond
 * the EoS minimal value down to h1.
 * @param[out] radius The radius of the star in m.
 * @param[out] mass The mass of the star in kg.
 * @param[out] love_number_k2 The k_2 tidal love number of the star.
 * @param[in] central_pressure_si The central pressure of the star in Pa.
 * @param eos Pointer to the Equation of State structure.
 * @param[in] epsrel The relative error for the TOV solver routine.
 * @retval 0 Success.
 * @retval <0 Failure.
 */
int XLALSimNeutronStarTOVODEIntegrateWithTolerance(double *radius, double *mass,
    double *love_number_k2, double central_pressure_si, LALSimNeutronStarEOS * eos, double epsrel)
{
    /* ode integration variables */
    const double epsabs = 0.0;
    double y[TOV_ODE_VARS_DIM] = {0.0};
    double dy[TOV_ODE_VARS_DIM] = {0.0};
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
        XLALSimNeutronStarEOSEnergyDensityOfPressureGeometrized(pc, eos);
    double hc =
        XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometrized(pc, eos);
    double dedp_c =
        XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrized(pc, eos);
    double dhdp_c = 1.0 / (ec + pc);
    double dedh_c = dedp_c / dhdp_c;
    double dh = -1e-12 * hc;
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
        int s = gsl_odeiv_evolve_apply(evolv, ctrl, step, &sys, &h, h1, &dh, y);
        if (s != GSL_SUCCESS) {
            gsl_odeiv_evolve_free(evolv);
            gsl_odeiv_control_free(ctrl);
            gsl_odeiv_step_free(step);
            XLAL_ERROR(XLAL_EERR,
                "Error encountered in GSL's ODE integrator\n");
        }
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


/**
 * @brief Integrates the stellar structure system of ordinary differential
 * equations to obtain neutron star's macroscopic parameters with a solver
 * relative precision 1e-6.
 * @details
 * Solves the Tolman-Oppenheimer-Volkoff stellar structure equations and
 * the gravito-electric tidal Love number equations, given an equation of
 * state and a central pressure. The solution to this set of ODEs provides
 * the neutron star's radius, gravitational mass and second order gravito-electric
 * tidal Love number.
 * It uses the pseudo-enthalpy formalism introduced in Lindblom (1992)
 * Astrophys. J. 398 569 "Determining the Nuclear Equation of State from
 * Neutron-Star Masses and Radii".
 * The surface of the star is defined by a minimum value of the pseudo-enthalpy
 * h1 = 1e-12 * hc (with hc the central pseudo-enthalpy) ; if the minimum
 * pseudo-enthalpy of the equation of state is larger, a polytrope
 * (degenerate non-relativistic gas) is used to extend the model beyond
 * the EoS minimal value down to h1.
 * @param[out] radius The radius of the star in m.
 * @param[out] mass The mass of the star in kg.
 * @param[out] love_number_k2 The k_2 tidal love number of the star.
 * @param[in] central_pressure_si The central pressure of the star in Pa.
 * @param eos Pointer to the Equation of State structure.
 * @retval 0 Success.
 * @retval <0 Failure.
 */
int XLALSimNeutronStarTOVODEIntegrate(double *radius, double *mass,
    double *love_number_k2, double central_pressure_si, LALSimNeutronStarEOS * eos)
{
    const double epsrel = 1e-6;
    return XLALSimNeutronStarTOVODEIntegrateWithTolerance(radius, mass, love_number_k2, central_pressure_si, eos, epsrel);
}

/*
 * For convenience, use a structure that provides a dictionary between
 * a vector of the TOV+virial ode variables and what they represent.
 */
struct tov_virial_ode_vars {
    double r;   /* radial coordinate, m */
    double m;   /* mass within r in geometric units, m */
    double H;   /* stellar perturbation in arbitrary units */
    double b;   /* derivative of metric pertubation in arbitrary units */
    double I1;  /* dependent variable for Virial ODEs */
    double I2;  /* dependent variable for Virial ODEs */
    double J1; /* dependent variable for Virial ODEs */
    double J2; /* dependent variable for Virial ODEs */

};

#define TOV_VIRIAL_ODE_VARS_DIM (sizeof(struct tov_virial_ode_vars)/sizeof(double))

/*
 * Casts an array of doubles to a structure with named parameter.
 */
static struct tov_virial_ode_vars *tov_virial_ode_vars_cast(const double *y)
{
    union {
        const double *y;
        struct tov_virial_ode_vars *v;
    } u = {
    y};
    return u.v;
}



/* ODE integrand for TOV equations and Virial equations with pseudo-enthalpy independent variable.
 * Implements Eqs. (5) and (6) of Lindblom, Astrophys. J. 398, 569 (1992).
 * Also uses Eqs. (7) and (8) [ibid] for inner boundary data, and
 * Eqs. (18), (27), (28) of Damour & Nagar, Phys. Rev. D 80, 084035 (2009)
 * [See also: Eqs. (11) & (12) Hinderer et al. Phys. Rev. D 81 123016 (2010)]
 * for the metric perturbation used to obtain the Love number.
 * For the Virial portion, implements equations provided by A. Nikolaidis, N. Stergioulas, H. Markakis. */
static int tov_virial_ode(double h, const double *y, double *dy, void *params)
{
    struct tov_virial_ode_vars *vars = tov_virial_ode_vars_cast(y);
    struct tov_virial_ode_vars *derivs = tov_virial_ode_vars_cast(dy);
    LALSimNeutronStarEOS *eos = params;

    double r = vars->r;
    double m = vars->m;
    double H = vars->H;
    double b = vars->b;
    // printf("Inside the tov virial ode: %.16e \n", b);
    double p =
        XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometrized(h, eos);
    double e =
        XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrized(h,
        eos);
    double dedp =
        XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrized(p, eos);
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
    /* Eq.(12) of Hinderer et al. Phys. Rev. D 81 123016 (2010). */
    double db = -(C0 * H + C1 * b) * dr;
    double alpha = 1.0 - 2.0 * m / r;
    double beta = (m + 4.0 * LAL_PI * r * r * r * p) / (r * r);

    double dI1 = 8.0 * LAL_PI * r * pow(alpha, (-1.0/2.0)) * p * dr;
    double dI2 = r * pow(alpha, (-1.5)) * pow(beta, (2.0)) * dr;
    double dJ1 = 4.0 * LAL_PI * r * r * pow(alpha, (-0.5)) * 3.0 * p * dr;
    double dJ2 = pow(alpha, (-0.5)) * (pow(alpha, (-1.0)) * pow((beta * r), (2.0)) - 0.5 * pow((sqrt(alpha) - 1.0), (2.0))) * dr;


    derivs->r = dr;
    derivs->m = dm;
    derivs->H = dH;
    derivs->b = db;
    derivs->I1 = dI1;
    derivs->I2 = dI2;
    derivs->J1 = dJ1;
    derivs->J2 = dJ2;
    return 0;
}

/**
 * @brief Integrates the stellar structure and virial system of ordinary differential
 * equations to obtain neutron star's macroscopic parameters and virial parameters.
 * @details
 * Solves the Tolman-Oppenheimer-Volkoff stellar structure equations,
 * the gravito-electric tidal Love number equations and the virial equations,
 * given an equation of state and a central pressure. The solution to this set
 * of ODEs provides the neutron star's radius, gravitational mass, second order
 * gravito-electric tidal Love number and the virial parameters.
 * It uses the pseudo-enthalpy formalism introduced in Lindblom (1992)
 * Astrophys. J. 398 569 "Determining the Nuclear Equation of State from
 * Neutron-Star Masses and Radii".
 * The surface of the star is defined by a minimum value of the pseudo-enthalpy
 * h1 = 1e-12 * hc (with hc the central pseudo-enthalpy) ; if the minimum
 * pseudo-enthalpy of the equation of state is larger, a polytrope
 * (degenerate non-relativistic gas) is used to extend the model beyond
 * the EoS minimal value down to h1.
 * @param[out] radius The radius of the star in m.
 * @param[out] mass The mass of the star in kg.
 * @param[out] int1 Virial parameter.
 * @param[out] int2 Virial parameter.
 * @param[out] int3 Virial parameter.
 * @param[out] int4 Virial parameter.
 * @param[out] int5 Virial parameter.
 * @param[out] int6 Virial parameter.
 * @param[out] love_number_k2 The k_2 tidal love number of the star.
 * @param[in] central_pressure_si The central pressure of the star in Pa.
 * @param eos Pointer to the Equation of State structure.
 * @param[in] epsrel The relative error in the TOV solver routine.
 * @retval 0 Success.
 * @retval <0 Failure.
 */
int XLALSimNeutronStarVirialODEIntegrateWithTolerance(double *radius, double *mass,
    double *int1, double *int2, double *int3, double *int4, double *int5, double *int6,
    double *love_number_k2, double central_pressure_si,
    LALSimNeutronStarEOS * eos, double epsrel)
{
    /* ode integration variables */
    const double epsabs = 0.0;
    double y[TOV_VIRIAL_ODE_VARS_DIM] = {0.0};
    double dy[TOV_VIRIAL_ODE_VARS_DIM] = {0.0};
    struct tov_virial_ode_vars *vars = tov_virial_ode_vars_cast(y);
    gsl_odeiv_system sys = { tov_virial_ode, NULL, TOV_VIRIAL_ODE_VARS_DIM, eos };
    gsl_odeiv_step *step =
        gsl_odeiv_step_alloc(gsl_odeiv_step_rk8pd, TOV_VIRIAL_ODE_VARS_DIM);
    gsl_odeiv_control *ctrl = gsl_odeiv_control_y_new(epsabs, epsrel);
    gsl_odeiv_evolve *evolv = gsl_odeiv_evolve_alloc(TOV_VIRIAL_ODE_VARS_DIM);

    /* central values */
    /* note: will be updated with Lindblom's series expansion */
    /* geometrisized units for variables in length (m) */
    double pc = central_pressure_si * LAL_G_C4_SI;
    double ec =
        XLALSimNeutronStarEOSEnergyDensityOfPressureGeometrized(pc, eos);
    double hc =
        XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometrized(pc, eos);
    double dedp_c =
        XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrized(pc, eos);
    double dhdp_c = 1.0 / (ec + pc);
    double dedh_c = dedp_c / dhdp_c;
    double dh = -1e-12 * hc;
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

    /* Virial ODEs starting points */

    double I1_0 = - 8.0 * LAL_PI * r0 * r0 * pc * dh - dh * dh;
    double I2_0 = - 16.0 * LAL_PI * LAL_PI * r0 * r0 * r0 * r0 * pc * pc * dh - 6.0 * LAL_PI * r0 * r0 * pc * dh * dh;
    double J1_0 = - 12.0 * LAL_PI * r0 * r0 * r0 * pc * dh - 3.0 * r0 * dh * dh * dh;
    double J2_0 = - 16.0 * LAL_PI * LAL_PI * r0 * r0 * r0 * r0 * r0 * pc * pc * dh + 8.0 * LAL_PI * LAL_PI * r0 * r0 * r0 * r0 * r0 * pc * pc * dh * dh;


    /* perform integration */
    vars->r = r0;
    vars->m = m0;
    vars->H = H0;
    vars->b = b0;
    vars->I1 = I1_0;
    vars->I2 = I2_0;
    vars->J1 = J1_0;
    vars->J2 = J2_0;

    h = h0;
    while (h > h1) {
        int s = gsl_odeiv_evolve_apply(evolv, ctrl, step, &sys, &h, h1, &dh, y);
        if (s != GSL_SUCCESS) {
            gsl_odeiv_evolve_free(evolv);
            gsl_odeiv_control_free(ctrl);
            gsl_odeiv_step_free(step);
            XLAL_ERROR(XLAL_EERR,
                "Error encountered in GSL's ODE integrator\n");
        }
    }

    /*take one final Euler step to get to surface*/
    for (int w = 0 ; w < 1 ; ++w){
        tov_virial_ode(h, y, dy, eos);
        for (i = 0; i < TOV_ODE_VARS_DIM; ++i)
            y[i] += dy[i] * (0.0 - h1);
    }


    /* compute tidal Love number k2 */
    c = vars->m / vars->r;      /* compactness */
    yy = vars->r * vars->b / vars->H;

    *int3 = (1.0 - vars->m / vars->r) * pow((1.0 - 2.0 * vars->m / vars->r), (-0.5)) - 1.0;
    *int6 = vars->r * (*int3);

    *int1 = vars->I1;
    *int2 = vars->I2;
    *int4 = vars->J1;
    *int5 = vars->J2;

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

/**
 * @brief Integrates the stellar structure and virial system of ordinary differential
 * equations to obtain neutron star's macroscopic parameters and virial parameters,
 * with a solver relative precision 1e-6.
 * @details
 * Solves the Tolman-Oppenheimer-Volkoff stellar structure equations,
 * the gravito-electric tidal Love number equations and the virial equations,
 * given an equation of state and a central pressure. The solution to this set
 * of ODEs provides the neutron star's radius, gravitational mass, second order
 * gravito-electric tidal Love number and the virial parameters.
 * It uses the pseudo-enthalpy formalism introduced in Lindblom (1992)
 * Astrophys. J. 398 569 "Determining the Nuclear Equation of State from
 * Neutron-Star Masses and Radii".
 * The surface of the star is defined by a minimum value of the pseudo-enthalpy
 * h1 = 1e-12 * hc (with hc the central pseudo-enthalpy) ; if the minimum
 * pseudo-enthalpy of the equation of state is larger, a polytrope
 * (degenerate non-relativistic gas) is used to extend the model beyond
 * the EoS minimal value down to h1.
 * @param[out] radius The radius of the star in m.
 * @param[out] mass The mass of the star in kg.
 * @param[out] int1 Virial parameter.
 * @param[out] int2 Virial parameter.
 * @param[out] int3 Virial parameter.
 * @param[out] int4 Virial parameter.
 * @param[out] int5 Virial parameter.
 * @param[out] int6 Virial parameter.
 * @param[out] love_number_k2 The k_2 tidal love number of the star.
 * @param[in] central_pressure_si The central pressure of the star in Pa.
 * @param eos Pointer to the Equation of State structure.
 * @retval 0 Success.
 * @retval <0 Failure.
 */
int XLALSimNeutronStarVirialODEIntegrate(double *radius, double *mass,
    double *int1, double *int2, double *int3, double *int4, double *int5, double *int6,
    double *love_number_k2, double central_pressure_si,
    LALSimNeutronStarEOS * eos)
{
    const double epsrel = 1e-6;
    return XLALSimNeutronStarVirialODEIntegrateWithTolerance(radius, mass,
    int1, int2, int3, int4, int5, int6, love_number_k2, central_pressure_si, eos, epsrel);
}


/*
 * For convenience, use a structure that provides a dictionary between
 * a vector of the extended TOV ode variables and what they represent.
 */
struct tov_ext_ode_vars {
    double r;   /* radial coordinate, m */
    double m;   /* mass within r in geometric units, m */
    double mb;  /* baryon mass, m */
    double H_k2;   /* stellar perturbation in arbitrary units */
    double b_k2;   /* derivative of metric perturbation in arbitrary units */
    double H_k3;   /* stellar perturbation in arbitrary units */
    double b_k3;   /* derivative of metric perturbation in arbitrary units */
    double H_k4;   /* stellar perturbation in arbitrary units */
    double b_k4;   /* derivative of metric perturbation in arbitrary units */
};

#define TOV_EXT_ODE_VARS_DIM (sizeof(struct tov_ext_ode_vars)/sizeof(double)) // Dimension of the ODE equation system

/*
 * Casts an array of doubles to a structure with named parameter.
 */
static struct tov_ext_ode_vars *tov_ext_ode_vars_cast(const double *y)
{
    union {
        const double *y;
        struct tov_ext_ode_vars *v;
    } u = {
    y};
    return u.v;
}

/*
 * Boundary condition for gravitational mass, baryonic mass,
 * radius and tidal deformability related parameters (H and beta for order 2, 3 and 4)
 * defined for the center of the star.
*/
static int tov_ext_initial_condition(double eps, double p, double dh, LALSimNeutronStarEOS * eos, struct tov_ext_ode_vars *variables){

    double dedp = XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrized(p, eos); // Central energy density derivative
    double dhdp = 1.0 / (eps + p);
    double dedh = dedp / dhdp;

    double rval = sqrt(-3.0 * dh / (2.0 * LAL_PI * (eps + 3.0 * p)));
    double mval = 4.0 * LAL_PI * rval * rval * rval * eps / 3.0;
    double Hval_k2 = rval * rval; // See section IV.A of Damour & Nagar
    double bval_k2 = 2.0 * rval;
    double Hval_k3 = rval * rval * rval;
    double bval_k3 = 3.0 * rval * rval ;
    double Hval_k4 = rval * rval * rval * rval;
    double bval_k4 = 4.0 * rval * rval * rval;

    /* series expansion for the initial core */
    /* second factor of Eq. (7) of Lindblom (1992) */
    rval *= 1.0 + 0.25 * dh * (eps - 3.0 * p  - 0.6 * dedh) / (eps + 3.0 * p);
    /* second factor of Eq. (8) of Lindblom (1992) */
    mval *= 1.0 + 0.6 * dh * dedh / eps;
    double mbval = mval;

    /* perform integration at first point */
    variables->r = rval;
    variables->m = mval;
    variables->mb = mbval;
    variables->H_k2 = Hval_k2;
    variables->b_k2 = bval_k2;
    variables->H_k3 = Hval_k3;
    variables->b_k3 = bval_k3;
    variables->H_k4 = Hval_k4;
    variables->b_k4 = bval_k4;

    return 0;
}

/*
 * ODE integrand for extended TOV and Love number (k2, k3 and k4) equations
 * with pseudo-enthalpy as the independent variable.
 */
static int tov_ext_ode(double h, const double *y, double *dy, void *params)
{
    struct tov_ext_ode_vars *vars = tov_ext_ode_vars_cast(y);
    struct tov_ext_ode_vars *derivs = tov_ext_ode_vars_cast(dy);

    LALSimNeutronStarEOS * eos = *(LALSimNeutronStarEOS **) params;

    double r = vars->r;
    double m = vars->m;
    double H_k2 = vars->H_k2;
    double b_k2 = vars->b_k2;
    double H_k3 = vars->H_k3;
    double b_k3 = vars->b_k3;
    double H_k4 = vars->H_k4;
    double b_k4 = vars->b_k4;

    double p = XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometrized(h, eos);
    double e = XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrized(h, eos);
    double rho = XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometrized(h, eos);
    double dedp = XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrized(p, eos);
    /* Eq. (18) of Damour & Nagar PRD 80 084035 (2009). */
    double A = 1.0 / (1.0 - 2.0 * m / r);
    /* Eq. (28) of Damour & Nagar PRD 80 084035 (2009). */
    double C1 = 2.0 / r + A * (2.0 * m / (r * r) + 4.0 * LAL_PI * r * (p - e));
    /* Eq. (29) of Damour & Nagar PRD 80 084035 (2009). */
    double C0_k2 = constant_zero_for_kl(r, m, A, e, p, dedp, 2);
    double C0_k3 = constant_zero_for_kl(r, m, A, e, p, dedp, 3);
    double C0_k4 = constant_zero_for_kl(r, m, A, e, p, dedp, 4);
    double dr = -r * (r - 2.0 * m) / (m + 4.0 * LAL_PI * r * r * r * p);
    double dm = 4.0 * LAL_PI * r * r * e * dr;
    double dmb = 4.0 * LAL_PI * r * r * rho * sqrt(A) * dr;

    double dH_k2 = b_k2 * dr;
    double db_k2 = -(C0_k2 * H_k2 + C1 * b_k2) * dr;
    double dH_k3 = b_k3 * dr;
    double db_k3 = -(C0_k3 * H_k3 + C1 * b_k3) * dr;
    double dH_k4 = b_k4 * dr;
    double db_k4 = -(C0_k4 * H_k4 + C1 * b_k4) * dr;

    derivs->r = dr;
    derivs->m = dm;
    derivs->mb = dmb;
    derivs->H_k2 = dH_k2;
    derivs->b_k2 = db_k2;
    derivs->H_k3 = dH_k3;
    derivs->b_k3 = db_k3;
    derivs->H_k4 = dH_k4;
    derivs->b_k4 = db_k4;
    return 0;
}

/*
 * Phase transition correction for the tidal love number.
 * Implements Eq.(41) of Pereira et al. 2020 ApJ 895 28
*/
static double correction_phase_transition(double r, double m, double b_kl, double H_kl,
                                          double h, LALSimNeutronStarEOS * eos0, LALSimNeutronStarEOS * eos1){

    double corrected_term = 0.0;
    double r3 = r*r*r;
    double pres_pt = XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometrized(h, eos1);
    double dpt_eps =  XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrized(h, eos1)
                    - XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrized(h, eos0) ;

    double rho_bar = (m + 4.0 * LAL_PI * r3 * pres_pt)/ (4.0 / 3.0 * LAL_PI * r3 ) ;
    corrected_term = b_kl - H_kl / r * dpt_eps * 3.0 / rho_bar ;
    return corrected_term;
}


/**
 * @brief Integrates the stellar structure system of ordinary differential
 * equations to obtain neutron star's macroscopic parameters; can accomodate
 * equations of state with phase transitions.
 * @details
 * Solves the Tolman-Oppenheimer-Volkoff stellar structure equations and
 * the gravito-electric tidal Love number equations, given an equation of
 * state and a central pressure. The solution to this set of ODEs provides
 * the neutron star's radius, gravitational mass, baryonic mass, and
 * second, third and fourth order gravito-electric tidal Love number;
 * it also accounts for corrections required for equations of state with
 * phase transitions.
 * It uses the pseudo-enthalpy formalism introduced in Lindblom (1992)
 * Astrophys. J. 398 569 "Determining the Nuclear Equation of State from
 * Neutron-Star Masses and Radii".
 * The surface of the star is defined by a minimum pseudo-enthalpy of the
 * input equation of state.
 * @param[out] radius The radius of the star in m.
 * @param[out] mass The gravitational mass of the star in kg.
 * @param[out] baryon_mass The baryon mass of the star in kg.
 * @param[out] love_number_k2 The k_2 tidal love number of the star.
 * @param[out] love_number_k3 The k_3 tidal love number of the star.
 * @param[out] love_number_k4 The k_4 tidal love number of the star.
 * @param[in] central_pressure_si The central pressure of the star in Pa.
 * @param eos Pointer to the multiple-parts Equation of State structure.
 * @param[in] epsrel The relative error for the numerical solver routine.
 */
void XLALSimNeutronStarTOVODEExtendedIntegrateWithTolerance(double *radius, double *mass, double *baryon_mass,
             double *love_number_k2, double *love_number_k3, double *love_number_k4,
             double central_pressure_si,
             EOSMultiParts *eos,
             double epsrel){

    /* ode integration variables */
    const double epsabs = 0.0;
    int s;
    LALSimNeutronStarEOS * params;
    double y[TOV_EXT_ODE_VARS_DIM];
    double dy[TOV_EXT_ODE_VARS_DIM];
    struct tov_ext_ode_vars *vars = tov_ext_ode_vars_cast(y);
    /* Set up the actual ODE system to solve with tov_virial_ode */
    gsl_odeiv_system sys = {tov_ext_ode, NULL, TOV_EXT_ODE_VARS_DIM, &params};
    /* Set up the iterative method to solve the ODE, and the dimension of the ODE set */
    gsl_odeiv_step *step = gsl_odeiv_step_alloc(gsl_odeiv_step_rk8pd, TOV_EXT_ODE_VARS_DIM);
    /* Set up the precision of the numerical method with relative error epsrel */
    gsl_odeiv_control *ctrl = gsl_odeiv_control_y_new(epsabs, epsrel);
    /* Set up evolution function to solve ODE with the dimension of the ODE set */
    gsl_odeiv_evolve *evolv = gsl_odeiv_evolve_alloc(TOV_EXT_ODE_VARS_DIM);


    /* Equations of state quantities at the center of the star */
    double pc = central_pressure_si * LAL_G_C4_SI; // convert the pressure from Pascal to geometrized units [/m^2]
    int number_eos = XLALSimNeutronStarEOSMultiPartsNumber(eos);
    int central_piece = 0; // default to first piece; handles floating-point edge case at pmin
    for(int j = 0; j < number_eos; j++){
        double pmin_piece = XLALSimNeutronStarEOSMultiPartsPieceMinPressureGeometrized(eos, j);
        double pmax_piece = XLALSimNeutronStarEOSMultiPartsPieceMaxPressureGeometrized(eos, j);
        if (pc >= pmin_piece && pc <= pmax_piece){
            central_piece = j;
            break;
        }
        if (pc > pmax_piece)
            central_piece = j; // track last piece whose range pc exceeded; handles floating-point edge case at pmax
    }
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, central_piece);
    double ec = XLALSimNeutronStarEOSEnergyDensityOfPressureGeometrized(pc, eos_part);
    double hc = XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometrized(pc, eos_part);

    /* Initial condition */
    double dh = -1e-12 * hc;
    double h0 = hc + dh;
    double hmin = 0.0;
    double h;
    tov_ext_initial_condition(ec, pc, dh, eos_part, vars);
    h = h0;

    for(int j = central_piece; j >= 0; j--) {
        LALSimNeutronStarEOS * eos_part_1 = XLALSimNeutronStarEOSPart(eos,j);
        params = eos_part_1;
        hmin = XLALSimNeutronStarEOSMinPseudoEnthalpy(eos_part_1);
        while (h > hmin) {
            s = gsl_odeiv_evolve_apply(evolv, ctrl, step, &sys, &h, hmin, &dh, y);
            if (s != GSL_SUCCESS) {
                gsl_odeiv_evolve_free(evolv);
                gsl_odeiv_control_free(ctrl);
                gsl_odeiv_step_free(step);
                XLAL_ERROR_VOID(XLAL_EERR, "Error encountered in GSL's ODE integrator\n");
            }
        }

        // Correction related to the phase transition
        if (j != 0){
            LALSimNeutronStarEOS * eos_part_0 = XLALSimNeutronStarEOSPart(eos,j-1);
            vars->b_k2 = correction_phase_transition(vars->r, vars->m, vars->b_k2, vars->H_k2,
                                          h,  eos_part_0,  eos_part_1);
            vars->b_k3 = correction_phase_transition(vars->r, vars->m, vars->b_k3, vars->H_k3,
                                          h,  eos_part_0,  eos_part_1);
            vars->b_k4 = correction_phase_transition(vars->r, vars->m, vars->b_k4, vars->H_k4,
                                          h,  eos_part_0,  eos_part_1);
        }

    }
    /*take one final Euler step to get to surface*/
    for (int w = 0 ; w < 1 ; ++w){
        eos_part = XLALSimNeutronStarEOSPart(eos,0);
        tov_ext_ode(h, y, dy, &eos_part);
        for (size_t i = 0; i < TOV_EXT_ODE_VARS_DIM; ++i) // No need for all Virial variables, only physically relevant ones
            y[i] += dy[i] * (0.0 - hmin);
    }

    /* compute tidal Love numbers and compactness at the surface of the star */
    double comp = vars->m / vars->r;      /* compactness */
    double yy_k2 = vars->r * vars->b_k2 / vars->H_k2; /* Eq. 13 of Hinderer et al. Phys. Rev. D 81 123016 */
    double yy_k3 = vars->r * vars->b_k3 / vars->H_k3;
    double yy_k4 = vars->r * vars->b_k4 / vars->H_k4;
    /* convert from geometric units to SI units */
    *radius = vars->r;
    *mass = vars->m * LAL_MSUN_SI / LAL_MRSUN_SI;
    *baryon_mass = vars->mb * LAL_MSUN_SI / LAL_MRSUN_SI;
    *love_number_k2 = tidal_Love_number_k2(comp, yy_k2);
    *love_number_k3 = tidal_Love_number_k3(comp, yy_k3);
    *love_number_k4 = tidal_Love_number_k4(comp, yy_k4);

    /* free ode memory */
    gsl_odeiv_evolve_free(evolv);
    gsl_odeiv_control_free(ctrl);
    gsl_odeiv_step_free(step);

    return;
}

/**
 * @brief Integrates the stellar structure system of ordinary differential
 * equations to obtain neutron star's macroscopic parameters, with a solver
 * relative precision 1e-6.; can accomodate equations of state with phase transitions.
 * @details
 * Solves the Tolman-Oppenheimer-Volkoff stellar structure equations and
 * the gravito-electric tidal Love number equations, given an equation of
 * state and a central pressure. The solution to this set of ODEs provides
 * the neutron star's radius, gravitational mass, baryonic mass, and
 * second, third and fourth order gravito-electric tidal Love number;
 * it also accounts for corrections required for equations of state with
 * phase transitions.
 * It uses the pseudo-enthalpy formalism introduced in Lindblom (1992)
 * Astrophys. J. 398 569 "Determining the Nuclear Equation of State from
 * Neutron-Star Masses and Radii".
 * The surface of the star is defined by a minimum pseudo-enthalpy of the
 * input equation of state.
 * @param[out] radius The radius of the star in m.
 * @param[out] mass The gravitational mass of the star in kg.
 * @param[out] baryon_mass The baryon mass of the star in kg.
 * @param[out] love_number_k2 The k_2 tidal love number of the star.
 * @param[out] love_number_k3 The k_3 tidal love number of the star.
 * @param[out] love_number_k4 The k_4 tidal love number of the star.
 * @param[in] central_pressure_si The central pressure of the star in Pa.
 * @param eos Pointer to the multiple-parts Equation of State structure.
 */
void XLALSimNeutronStarTOVODEExtendedIntegrate(double *radius, double *mass, double *baryon_mass,
             double *love_number_k2, double *love_number_k3, double *love_number_k4,
             double central_pressure_si,
             EOSMultiParts *eos)
{
    const double epsrel = 1e-6;
    XLALSimNeutronStarTOVODEExtendedIntegrateWithTolerance(radius, mass, baryon_mass, love_number_k2,
                                                           love_number_k3, love_number_k4, central_pressure_si,
                                                           eos, epsrel);
    return ;
}

/*
 * For convenience, use a structure that provides a dictionary between
 * a vector of the mini TOV ode variables and what they represent.
 */
struct tov_mini_ode_vars {
    double r;   /* radial coordinate, m */
    double m;   /* mass within r in geometric units, m */
    double H_k2;   /* stellar perturbation in arbitrary units */
    double b_k2;   /* derivative of metric perturbation in arbitrary units */
};

#define TOV_MINI_ODE_VARS_DIM (sizeof(struct tov_mini_ode_vars)/sizeof(double)) // Dimension of the ODE equation system

/*
 * Casts an array of doubles to a structure with named parameter.
 */
static struct tov_mini_ode_vars *tov_mini_ode_vars_cast(const double *y)
{
    union {
        const double *y;
        struct tov_mini_ode_vars *v;
    } u = {
    y};
    return u.v;
}


/*
 * Boundary condition for gravitational mass,
 * radius and tidal deformability related parameters (H and beta)
 * defined for the center of the star.
*/
static int tov_mini_initial_condition(double eps, double p, double dh, LALSimNeutronStarEOS * eos, struct tov_mini_ode_vars *variables){

    double dedp = XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrized(p, eos); // Central energy density derivative
    double dhdp = 1.0 / (eps + p);
    double dedh = dedp / dhdp;

    double rval = sqrt(-3.0 * dh / (2.0 * LAL_PI * (eps + 3.0 * p)));
    double mval = 4.0 * LAL_PI * rval * rval * rval * eps / 3.0;
    double Hval_k2 = rval * rval;
    double bval_k2 = 2.0 * rval;

    /* series expansion for the initial core */
    /* second factor of Eq. (7) of Lindblom (1992) */
    rval *= 1.0 + 0.25 * dh * (eps - 3.0 * p  - 0.6 * dedh) / (eps + 3.0 * p);
    /* second factor of Eq. (8) of Lindblom (1992) */
    mval *= 1.0 + 0.6 * dh * dedh / eps;

    /* perform integration at first point */
    variables->r = rval;
    variables->m = mval;
    variables->H_k2 = Hval_k2;
    variables->b_k2 = bval_k2;

    return 0;
}

/*
 * ODE integrand for the TOV and Love number (k2) equations with pseudo-enthalpy as
 * the independent variable (pertinent for the mini TOV solver).
 */
static int tov_mini_ode(double h, const double *y, double *dy, void *params)
{
    struct tov_mini_ode_vars *vars = tov_mini_ode_vars_cast(y);
    struct tov_mini_ode_vars *derivs = tov_mini_ode_vars_cast(dy);

    LALSimNeutronStarEOS * eos = *(LALSimNeutronStarEOS **) params;

    double r = vars->r;
    double m = vars->m;
    double H_k2 = vars->H_k2;
    double b_k2 = vars->b_k2;

    double p = XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometrized(h, eos);
    double e = XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometrized(h, eos);
    double dedp = XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometrized(p, eos);
    /* Eq. (18) of Damour & Nagar PRD 80 084035 (2009). */
    double A = 1.0 / (1.0 - 2.0 * m / r);
    /* Eq. (28) of Damour & Nagar PRD 80 084035 (2009). */
    double C1 = 2.0 / r + A * (2.0 * m / (r * r) + 4.0 * LAL_PI * r * (p - e));
    /* Eq. (29) of Damour & Nagar PRD 80 084035 (2009). */
    double C0_k2 = constant_zero_for_kl(r, m, A, e, p, dedp, 2);
    double dr = -r * (r - 2.0 * m) / (m + 4.0 * LAL_PI * r * r * r * p);
    double dm = 4.0 * LAL_PI * r * r * e * dr;

    double dH_k2 = b_k2 * dr;
    double db_k2 = -(C0_k2 * H_k2 + C1 * b_k2) * dr;

    derivs->r = dr;
    derivs->m = dm;
    derivs->H_k2 = dH_k2;
    derivs->b_k2 = db_k2;
    return 0;
}


/**
 * @brief Integrates the stellar structure system of ordinary differential
 * equations to obtain a minimum number of neutron star's macroscopic parameters;
 * can accomodate equations of state with phase transitions.
 * @details
 * Solves the Tolman-Oppenheimer-Volkoff stellar structure equations and
 * the gravito-electric tidal Love number equations, given an equation of
 * state and a central pressure. The solution to this set of ODEs provide
 * the neutron star's radius, gravitational mass and
 * second order gravito-electric tidal Love number; it also accounts for
 * corrections required for equations of state with phase transitions.
 * It uses the pseudo-enthalpy formalism introduced in Lindblom (1992)
 * Astrophys. J. 398 569 "Determining the Nuclear Equation of State from
 * Neutron-Star Masses and Radii".
 * The surface of the star is defined by a minimum pseudo-enthalpy of the
 * input equation of state.
 * @param[out] radius The radius of the star in m.
 * @param[out] mass The gravitationnal mass of the star in kg.
 * @param[out] love_number_k2 The k_2 tidal love number of the star.
 * @param[in] central_pressure_si The central pressure of the star in Pa.
 * @param eos Pointer to the multiple-parts Equation of State structure.
 * @param[in] epsrel The relative error for the numerical solver routine.
 */
void XLALSimNeutronStarTOVODEMiniIntegrateWithTolerance(double *radius, double *mass,
             double *love_number_k2, double central_pressure_si,
             EOSMultiParts *eos,
             double epsrel){

    /* ode integration variables */
    const double epsabs = 0.0;
    int s;
    LALSimNeutronStarEOS * params;
    double y[TOV_MINI_ODE_VARS_DIM]; // array for the variables of the ODE equations
    double dy[TOV_MINI_ODE_VARS_DIM];
    struct tov_mini_ode_vars *vars = tov_mini_ode_vars_cast(y);
    /* Set up the actual ODE system to solve with tov_virial_ode */
    gsl_odeiv_system sys = {tov_mini_ode, NULL, TOV_MINI_ODE_VARS_DIM, &params};
    /* Set up the iterative method to solve the ODE, and the dimension of the ODE set */
    gsl_odeiv_step *step = gsl_odeiv_step_alloc(gsl_odeiv_step_rk8pd, TOV_MINI_ODE_VARS_DIM);
    /* Set up the precision of the numerical method with relative error epsrel */
    gsl_odeiv_control *ctrl = gsl_odeiv_control_y_new(epsabs, epsrel);
    /* Set up evolution function to solve ODE with the dimension of the ODE set */
    gsl_odeiv_evolve *evolv = gsl_odeiv_evolve_alloc(TOV_MINI_ODE_VARS_DIM);

    /* Equations of state quantities at the center of the star */
    double pc = central_pressure_si * LAL_G_C4_SI; // convert the pressure from Pascal to geometrized units [/m^2]
    int number_eos = XLALSimNeutronStarEOSMultiPartsNumber(eos);
    int central_piece = 0; // default to first piece; handles floating-point edge case at pmin
    for(int j = 0; j < number_eos; j++){
        double pmin_piece = XLALSimNeutronStarEOSMultiPartsPieceMinPressureGeometrized(eos, j);
        double pmax_piece = XLALSimNeutronStarEOSMultiPartsPieceMaxPressureGeometrized(eos, j);
        if (pc >= pmin_piece && pc <= pmax_piece){
            central_piece = j;
            break;
        }
        if (pc > pmax_piece)
            central_piece = j; // track last piece whose range pc exceeded; handles floating-point edge case at pmax
    }
    LALSimNeutronStarEOS * eos_part = XLALSimNeutronStarEOSPart(eos, central_piece);
    double ec = XLALSimNeutronStarEOSEnergyDensityOfPressureGeometrized(pc, eos_part);
    double hc = XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometrized(pc, eos_part);


    /* Initial condition */
    double dh = -1e-12 * hc;
    double h0 = hc + dh;
    double hmin = 0.0;
    double h;
    tov_mini_initial_condition(ec, pc, dh, eos_part, vars);
    h = h0;

    for(int j = central_piece; j >= 0; j--) {
        LALSimNeutronStarEOS * eos_part_1 = XLALSimNeutronStarEOSPart(eos,j);
        params = eos_part_1;
        hmin = XLALSimNeutronStarEOSMinPseudoEnthalpy(eos_part_1);
        while (h > hmin) {
            s = gsl_odeiv_evolve_apply(evolv, ctrl, step, &sys, &h, hmin, &dh, y);
            if (s != GSL_SUCCESS) {
                gsl_odeiv_evolve_free(evolv);
                gsl_odeiv_control_free(ctrl);
                gsl_odeiv_step_free(step);
                XLAL_ERROR_VOID(XLAL_EERR, "Error encountered in GSL's ODE integrator\n");
            }
        }
        // Correction related to the phase transition
        if (j != 0) {
            LALSimNeutronStarEOS * eos_part_0 = XLALSimNeutronStarEOSPart(eos,j-1);
            vars->b_k2 = correction_phase_transition(vars->r, vars->m, vars->b_k2, vars->H_k2,
                                          h,  eos_part_0,  eos_part_1);
        }

    }
    /*take one final Euler step to get to surface*/
    for (int w = 0 ; w < 1 ; ++w){
        eos_part = XLALSimNeutronStarEOSPart(eos,0);
        tov_mini_ode(h, y, dy, &eos_part);
        for (size_t i = 0; i < TOV_MINI_ODE_VARS_DIM; ++i) // No need for all Virial variables, only physically relevant ones
            y[i] += dy[i] * (0.0 - hmin);
    }

    /* compute tidal Love numbers and compactness at the surface of the star */
    double comp = vars->m / vars->r;      /* compactness */
    double yy_k2 = vars->r * vars->b_k2 / vars->H_k2; /* Eq. 13 of Hinderer et al. Phys. Rev. D 81 123016 */
    /* convert from geometric units to SI units */
    *radius = vars->r;
    *mass = vars->m * LAL_MSUN_SI / LAL_MRSUN_SI;
    *love_number_k2 = tidal_Love_number_k2(comp, yy_k2);

    /* free ode memory */
    gsl_odeiv_evolve_free(evolv);
    gsl_odeiv_control_free(ctrl);
    gsl_odeiv_step_free(step);

    return ;
}


/**
 * @brief Integrates the stellar structure system of ordinary differential
 * equations to obtain a minimum number of neutron star's macroscopic parameters,
 * with a solver relative precision 1e-6.; can accomodate equations of state with
 * phase transitions.
 * @details
 * Solves the Tolman-Oppenheimer-Volkoff stellar structure equations and
 * the gravito-electric tidal Love number equations, given an equation of
 * state and a central pressure. The solution to this set of ODEs provide
 * the neutron star's radius, gravitational mass and
 * second order gravito-electric tidal Love number; it also accounts for
 * corrections required for equations of state with phase transitions.
 * It uses the pseudo-enthalpy formalism introduced in Lindblom (1992)
 * Astrophys. J. 398 569 "Determining the Nuclear Equation of State from
 * Neutron-Star Masses and Radii".
 * The surface of the star is defined by a minimum pseudo-enthalpy of the
 * input equation of state.
 * @param[out] radius The radius of the star in m.
 * @param[out] mass The gravitationnal mass of the star in kg.
 * @param[out] love_number_k2 The k_2 tidal love number of the star.
 * @param[in] central_pressure_si The central pressure of the star in Pa.
 * @param eos Pointer to the multiple-parts Equation of State structure.
 */
void XLALSimNeutronStarTOVODEMiniIntegrate(double *radius, double *mass, double *love_number_k2,
                                       double central_pressure_si, EOSMultiParts *eos)
{
    const double epsrel = 1e-6;
    XLALSimNeutronStarTOVODEMiniIntegrateWithTolerance(radius, mass, love_number_k2, central_pressure_si, eos, epsrel);
    return ;
}


/** @} */
