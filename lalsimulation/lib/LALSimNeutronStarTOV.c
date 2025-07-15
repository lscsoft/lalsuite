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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
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


//CUTER-dev // TODO why can I not name this LAL blabla ??
/* Structure containing phase transition information (hpt, ppt, delta_eps),
 * one EOS structure (eos_low) for the equation of state before the phase transition,
 * and a one EOS structure (eos_up) for the eos after the phase transition.
 */
struct EOSTwoPartsWithPTinfo{
  LALSimNeutronStarEOS * eos_low;
  LALSimNeutronStarEOS * eos_up;
  double hpt;
  double ppt;
  double delta_eps;
};


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

//CUTER-dev
/* Implements Eq.(51) of Damour & Nagar, Phys. Rev. D 80 084035 (2009).*/
static double tidal_Love_number_k3(double c, double y){
    double num;
    double den;
    double c2 = pow(c, 2.0);
    double c3 = pow(c, 3.0);
    double c4 = pow(c, 4.0);
    double c5 = pow(c, 5.0);
    double c7 = pow(c, 7.0);

    num = (8.0 / 7.0) * pow(1 - 2 * c, 2.0) * c7
        * (2 * (y - 1) * c2 - 3.0 * (y - 2.0) * c + y - 3.0);
    den = 2.0 * c * (  4*(y + 1)*c5 + 2*(9*y - 2)*c4
                        - 20*(7*y - 9)*c3 + 5*(37*y - 72)*c2
                        - 45*(2*y - 5)*c  + 15*(y - 3) ) ;

    den -= 15*pow(1 - 2 * c, 2.0) * ( 2*(y - 1)*c2 - 3*(y-2)*c + y - 3 ) * log(1.0 / (1 - 2 * c));

    return num / den;
}

//CUTER-dev
/* Implements Eq.(52) of Damour & Nagar, Phys. Rev. D 80 084035 (2009).*/
static double tidal_Love_number_k4(double c, double y){
    double num;
    double den;
    double c2 = pow(c, 2.0);
    double c3 = pow(c, 3.0);
    double c4 = pow(c, 4.0);
    double c5 = pow(c, 5.0);
    double c6 = pow(c, 6.0);
    double c9 = pow(c, 9.0);

    num = (32.0 / 147.0) * pow(1 - 2 * c, 2.0) * c9
        * ( 12*(y - 1)*c3 - 34*(y - 2)*c2 + 28*(y - 3)*c - 7*(y - 4) );
    den = 2.0 * c * (  8*(y + 1)*c6 + (68*y - 8)*c5 + (1284 - 996*y)*c4
                        + 40*(55*y - 116)*c3 + (5360 - 1910*y)*c2
                        + 105*(7*y - 24)*c - 105*(y - 4) ) ;

    den -= 15*pow(1 - 2 * c, 2.0) * ( 12*(y - 1)*c3 - 34*(y - 2)*c2 + 28*(y - 3)*c - 7*(y - 4) ) * log(1.0 / (1 - 2 * c));

    return num / den;
}



//CUTER-dev
static double constant_zero_for_kl(double r, double m, double A, double e, double p, double dedp, int l){
    /* Eq. (29) of Damour & Nagar PRD 80 084035 (2009). */
    double C0 =
    A * (-(l) * (l + 1) / (r * r) + 4.0 * LAL_PI * (e + p) * dedp +
        4.0 * LAL_PI * (5.0 * e + 9.0 * p)) - pow(2.0 * (m +
            4.0 * LAL_PI * r * r * r * p) / (r * (r - 2.0 * m)), 2.0);

    return C0;
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
 * @param[in] epsrel The relative error for the TOV solver routine
 * @retval 0 Success.
 * @retval <0 Failure.
 */
int XLALSimNeutronStarTOVODEIntegrateWithTolerance(double *radius, double *mass,
    double *love_number_k2, double central_pressure_si,
    LALSimNeutronStarEOS * eos, double epsrel)
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
        XLALSimNeutronStarEOSEnergyDensityOfPressureGeometerized(pc, eos);
    double hc =
        XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometerized(pc, eos);
    double dedp_c =
        XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometerized(pc, eos);
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
    // printf("Initial value of h = %.6e\n", h);
    while (h > h1) {
        // printf("Star integration h= %.16e \t M = %.6e \n", h, vars->m  / LAL_MRSUN_SI);
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
    double db = -(C0 * H + C1 * b) * dr; // CUTER-dev exactly Eq.(12) in Hindered 2008 https://arxiv.org/pdf/0911.3535

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
 * @brief Integrates the Tolman-Oppenheimer-Volkov stellar structure equations and the Virial Equations.
 * @details
 * Solves the Tolman-Oppenheimer-Volkov stellar structure equations using the
 * pseudo-enthalpy formalism introduced in:
 * Lindblom (1992) "Determining the Nuclear Equation of State from Neutron-Star
 * Masses and Radii", Astrophys. J. 398 569.
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
 * @param eos EoS structure (multiple part EOS)
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
    double ec = // TODO find a way to add the new value of Eps after PT.
        XLALSimNeutronStarEOSEnergyDensityOfPressureGeometerized(pc, eos);
    double hc =
        XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometerized(pc, eos);
    double dedp_c =
        XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometerized(pc, eos);
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

    // double Gamma_c = (ec + pc)/(pc * dedp_c);
    // double r1 = sqrt(3.0/(2.0 * LAL_PI * (ec + 3.0 * pc)));
    // double r3 = - (r1 / (4.0 * (ec + 3.0 * pc))) * (ec - 3.0 * pc - 3.0 * (ec + pc) * (ec + pc) /  (5.0 * pc * Gamma_c));
    // double m3 = 4.0 * LAL_PI * ec * r1 * r1 * r1 / 3.0;
    // double m5 = 4.0 * LAL_PI * r1 * r1 * r1 * (r3 * ec / r1 - (ec + pc) * (ec + pc) / (5.0 * pc * Gamma_c));

    // r0 = r1 * sqrt(fabs(dh)) + r3 * pow(sqrt(fabs(dh)), 3.0);
    // m0 = m3 * pow(sqrt(fabs(dh)), 3.0) + m5 * pow(sqrt(fabs(dh)), 5.0);

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
        // printf("Old TOV solver star integration h= %.16e \t M = %.6e \n", h, vars->m  / LAL_MRSUN_SI);
        int s =
            gsl_odeiv_evolve_apply(evolv, ctrl, step, &sys, &h, h1, &dh, y);
        if (s != GSL_SUCCESS)
            XLAL_ERROR(XLAL_EERR,
                "Error encountered in GSL's ODE integrator\n");

//         printf("\t\t Star integration h= %.16e \t dh= %.16e \t M = %.6e \t I2 = %.6e\n", h, dh, vars->m  / LAL_MRSUN_SI, vars->I2);
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


//CUTER-dev
/* Structure that contains the set of physical variables to solve with ODE solver */
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

static struct tov_ext_ode_vars *tov_ext_ode_vars_cast(const double *y)
{
    union {
        const double *y;
        struct tov_ext_ode_vars *v;
    } u = {
    y};
    return u.v;
}


// CUTER-dev
/* Boundary condition for gravitaitonal mass, baryonic mass,
 * radius and tidal deformability related parameters (H and beta)
 * defined for the center of the star.
*/
static int tov_ext_initial_condition(double eps, double p, double dh, LALSimNeutronStarEOS * eos, struct tov_ext_ode_vars *variables){

    double dedp = XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometerized(p, eos); // Central energy density derivative
    double dhdp = 1.0 / (eps + p);
    double dedh = dedp / dhdp;

    double rval = sqrt(-3.0 * dh / (2.0 * LAL_PI * (eps + 3.0 * p)));
    double mval = 4.0 * LAL_PI * rval * rval * rval * eps / 3.0;
    double Hval_k2 = rval * rval;
    double bval_k2 = 2.0 * rval;
    double Hval_k3 = rval * rval;
    double bval_k3 = 2.0 * rval;
    double Hval_k4 = rval * rval;
    double bval_k4 = 2.0 * rval; // TODO check the initial conditions ?

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


//CUTER-dev
/* ODE integrand for TOV (including the baryonic mass) and
 * Love number (k2, k3 and k4) equations with pseudo-enthalpy as
 * the independent variable.
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

    double p = XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometerized(h, eos);
    double e = XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometerized(h, eos);
    double rho = XLALSimNeutronStarEOSRestMassDensityOfPseudoEnthalpyGeometerized(h, eos);
    double dedp = XLALSimNeutronStarEOSEnergyDensityDerivOfPressureGeometerized(p, eos);
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

//     printf("\t\tinside the tov_ode : h = %.6e p = %.6e \n", h, p);

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


//CUTER-dev
/**
 * @brief Integrates the Tolman-Oppenheimer-Volkov stellar structure equations.
 * @details
 * Solves the Tolman-Oppenheimer-Volkov stellar structure equations using the
 * pseudo-enthalpy formalism introduced in:
 * Lindblom (1992) "Determining the Nuclear Equation of State from Neutron-Star
 * Masses and Radii", Astrophys. J. 398 569.
 * @param[out] radius The radius of the star in m.
 * @param[out] mass The mass of the star in kg.
 * @param[out] baryon_mass The baryon mass of the star in kg.
 * @param[out] love_number_k2 The k_2 tidal love number of the star.
 * @param[out] love_number_k3 The k_3 tidal love number of the star.
 * @param[out] love_number_k4 The k_4 tidal love number of the star.
 * @param[in] central_pressure_si The central pressure of the star in Pa.
 * @param eos Pointer to the Equation of State structure with multiple parts.
 * @param[in] epsrel The relative error for the TOV solver routine
 * @retval 0 Success.
 * @retval <0 Failure.
 */
int XLALSimNeutronStarTOVODEExtendedIntegrateWithTolerance(double *radius, double *mass, double *baryon_mass,
             double *love_number_k2, double *love_number_k3, double *love_number_k4,
             double central_pressure_si,
             struct EOSMultiParts eos,
             double epsrel){

    /* ode integration variables */
    const double epsabs = 0.0;
    double y[TOV_EXT_ODE_VARS_DIM]; // array for the variables of the ODE equations
    double dy[TOV_EXT_ODE_VARS_DIM];
    size_t i;
    int s;
    LALSimNeutronStarEOS * params;


    struct tov_ext_ode_vars *vars = tov_ext_ode_vars_cast(y);
    /* Set up the actual ODE system to solve with tov_virial_ode */
    gsl_odeiv_system sys = {tov_ext_ode, NULL, TOV_EXT_ODE_VARS_DIM, &params};
    /* Set up the iterative method to solve the ODE, and the dimension of the ODE set */
    gsl_odeiv_step *step = gsl_odeiv_step_alloc(gsl_odeiv_step_rk8pd, TOV_EXT_ODE_VARS_DIM);
    /* Set up the precision of the numerical method with relative error epsrel */
    gsl_odeiv_control *ctrl = gsl_odeiv_control_y_new(epsabs, epsrel); // TODO test the error absolute vs relative
    /* Set up evolution function to solve ODE with the dimension of the ODE set */
    gsl_odeiv_evolve *evolv = gsl_odeiv_evolve_alloc(TOV_EXT_ODE_VARS_DIM);

    /* Equation of state */
    int number_eos = eos.number_of_parts;
    double *hmin = eos.hmin;
//     double *hmax = eos.hmax;
//     printf("In essaiTOV:\n\tNumber of parts in the EoS: %d \n", number_eos);

    /* Central values */
    double pc = central_pressure_si * LAL_G_C4_SI; // convert the pressure from Pascal to geometrized units [/m^2]
    double ec = XLALSimNeutronStarEOSEnergyDensityOfPressureGeometerized(pc, eos.eos_part[number_eos-1]);
    double hc = XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometerized(pc, eos.eos_part[number_eos-1]);


    /* Initial condition */
    double dh = -1e-12 * hc;
    double h0 = hc + dh;
    double h1 = 0.0 - dh;
    double h;

    tov_ext_initial_condition(ec, pc, dh, eos.eos_part[number_eos-1], vars);
//     printf("\tInitial values of astro parameters: m = %g \t r = %g\n", vars->m, vars->r);
    h = h0;
//     printf("\tInitial values of micro parameters: h = %.2e\n", h);

    for(int j = number_eos-1; j >= 0; j--) {
        params = eos.eos_part[j];
        double pmin = XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometerized(hmin[j], eos.eos_part[j]);
//         double pmax = XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometerized(hmax[j], eos.eos_part[j]);

//         printf("\t\tPart %d\tPmax = %.2e\thmin = %.2e\n", j, params->pmax * FAC_NUC, hmin[j]);

//         printf("\t\th_initial = %.2e \t h_min = %.2e\n", h, hmin[j]);

//         printf("Pc = %g pmax=%g \n", pc, pmax);
        while (h > hmin[j]) {
//             printf("here %.6e %.6e %.6e\n", h, dh, vars->m/LAL_MRSUN_SI);
            s = gsl_odeiv_evolve_apply(evolv, ctrl, step, &sys, &h, hmin[j], &dh, y);
            if (s != GSL_SUCCESS)
                XLAL_ERROR(XLAL_EERR,
                    "Error encountered in GSL's ODE integrator\n");
//             printf("\t\t %d Star integration h= %.16e \t M = %.6e \n", j,  h, vars->m  / LAL_MRSUN_SI);
            printf("\t\t %d Star integration h= %.16e \t dh= %.16e \t M = %.6e \n", j,  h, dh, vars->m  / LAL_MRSUN_SI);
        }

        if (pc >= pmin && j != 0){
//             printf("hmin =%6e \nValues of b before correction %.6e for k2\n", hmin[j], vars->b_k2);
            /* Phase transition correction for the tidal love number
            * Implements Eq.(41) of Pereira et al. 2020 ApJ 895 28
            * Note: Eq.(14) of Postnikov et al. 2010 Phys. Rev. D 82, 024016
            * takes an approximation for the denominator on the correction terms
            */
            double r3 = (vars->r)*(vars->r)*(vars->r);
            double pres_pt = XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometerized(h, eos.eos_part[j]);
            double dpt_eps =  XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometerized(h, eos.eos_part[j])
                            - XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometerized(h, eos.eos_part[j-1]) ;

            double rho_bar = (vars->m + 4.0 * LAL_PI * r3 * pres_pt)/ (4.0 / 3.0 * LAL_PI * r3 ) ;
            vars->b_k2 += - vars->H_k2 / vars->r * dpt_eps * 3.0 / rho_bar ;
            vars->b_k3 += - vars->H_k3 / vars->r * dpt_eps * 3.0 / rho_bar ;
            vars->b_k4 += - vars->H_k4 / vars->r * dpt_eps * 3.0 / rho_bar ;
//             printf("Values of b after correction %.6e for k2\n", vars->b_k2);
//             printf("\t\t %d CORRECTION ON h= %.16e \t M = %.6e \n", j,h, vars->m/LAL_MRSUN_SI);

        }

    }

    /* compute tidal Love numbers and compactness at the surface of the star */
    double comp = vars->m / vars->r;      /* compactness */
    double yy_k2 = vars->r * vars->b_k2 / vars->H_k2; /* Eq. 13 of Hinderer et al. Phys. Rev. D 81 123016 */
    double yy_k3 = vars->r * vars->b_k3 / vars->H_k3; /* Eq. 13 of Hinderer et al. Phys. Rev. D 81 123016 */
    double yy_k4 = vars->r * vars->b_k4 / vars->H_k4; /* Eq. 13 of Hinderer et al. Phys. Rev. D 81 123016 */

    /*take one final Euler step to get to surface*/
    for (int w = 0 ; w < 1 ; ++w){
        tov_ext_ode(h, y, dy, &eos.eos_part[0]);
        for (i = 0; i < TOV_EXT_ODE_VARS_DIM; ++i) // No need for all Virial variables, only physically relevant ones
            y[i] += dy[i] * (0.0 - h1);
    }
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

    return 0;
}

// CUTER-dev
struct keeping_data{
    double old_r;
    double old_m;
    LALSimNeutronStarEOS *eos;
};


//CUTER-dev
struct virial_ode_vars {
    double I1;  /* dependent variable for Virial ODEs */
    double I2;  /* dependent variable for Virial ODEs */
    double J1;  /* dependent variable for Virial ODEs */
    double J2;  /* dependent variable for Virial ODEs */
};

#define VIRIAL_ODE_VARS_DIM (sizeof(struct virial_ode_vars)/sizeof(double))

static struct virial_ode_vars *virial_ode_vars_cast(const double *y)
{
    union {
        const double *y;
        struct virial_ode_vars *v;
    } u = {
    y};
    return u.v;
}


// CUTER-dev
/* Boundary condition for the virial parameters,
 * defined for the center of the star.
 */
static int virial_initial_condition(double r0, double pc, double dh, struct virial_ode_vars *variables){


    double I1_0 = - 8.0 * LAL_PI * r0 * r0 * pc * dh - dh * dh;
    double I2_0 = - 16.0 * LAL_PI * LAL_PI * r0 * r0 * r0 * r0 * pc * pc * dh - 6.0 * LAL_PI * r0 * r0 * pc * dh * dh;
    double J1_0 = - 12.0 * LAL_PI * r0 * r0 * r0 * pc * dh - 3.0 * r0 * dh * dh * dh;
    double J2_0 = - 16.0 * LAL_PI * LAL_PI * r0 * r0 * r0 * r0 * r0 * pc * pc * dh + 8.0 * LAL_PI * LAL_PI * r0 * r0 * r0 * r0 * r0 * pc * pc * dh * dh;

    /* perform integration at first point */
    variables->I1 = I1_0;
    variables->I2 = I2_0;
    variables->J1 = J1_0;
    variables->J2 = J2_0;
    return 0;
}


// CUTER-dev

/* ODE integrand for TOV equations and Virial equations with pseudo-enthalpy independent variable.
 * Implements Eqs. (5) and (6) of Lindblom, Astrophys. J. 398, 569 (1992).
 * Also uses Eqs. (7) and (8) [ibid] for inner boundary data, and
 * Eqs. (18), (27), (28) of Damour & Nagar, Phys. Rev. D 80, 084035 (2009)
 * [See also: Eqs. (11) & (12) Hinderer et al. Phys. Rev. D 81 123016 (2010)]
 * for the metric perturbation used to obtain the Love number.
 * For the Virial portion, implements equations provided by A. Nikolaidis, N. Stergioulas, H. Markakis. */
static int virial_ode(double h, const double *y, double *dy, void *params)
{
    struct virial_ode_vars *vars = virial_ode_vars_cast(y);
    struct virial_ode_vars *derivs = virial_ode_vars_cast(dy);

    struct keeping_data *data = params;
//     printf("here if h = %.6e and %.6e\n", h, vars->I1);

//     if (h < 0.0) h = 1.e-10;
    double r = data->old_r;
    double m = data->old_m;
    double p = XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometerized(h, data->eos);
    double dr = -r * (r - 2.0 * m) / (m + 4.0 * LAL_PI * r * r * r * p);

    double alpha = 1.0 - 2.0 * m / r;
    double beta = (m + 4.0 * LAL_PI * r * r * r * p) / (r * r);

    double dI1 = 8.0 * LAL_PI * r * pow(alpha, (-1.0/2.0)) * p * dr;
    double dI2 = r * pow(alpha, (-1.5)) * pow(beta, (2.0)) * dr;
    double dJ1 = 4.0 * LAL_PI * r * r * pow(alpha, (-0.5)) * 3.0 * p * dr;
    double dJ2 = pow(alpha, (-0.5)) * (pow(alpha, (-1.0)) * pow((beta * r), (2.0)) - 0.5 * pow((sqrt(alpha) - 1.0), (2.0))) * dr;

    vars->I1 = vars->I1; // TODO degueu pour retirer le unused variable
    vars->I2 = vars->I2;
    vars->J1 = vars->J1;
    vars->J2 = vars->J2;
//     printf("\t\tinside the virial_ode : h = %.6e p = %.6e \n", h, p);

    derivs->I1 = dI1;
    derivs->I2 = dI2;
    derivs->J1 = dJ1;
    derivs->J2 = dJ2;
    return 0;
}

//CUTER-dev
/**
 * @brief Integrates the Tolman-Oppenheimer-Volkov stellar structure equations
 * and the virial quantities to control the solver error.
 * @details
 * Solves the Tolman-Oppenheimer-Volkov stellar structure equations using the
 * pseudo-enthalpy formalism introduced in:
 * Lindblom (1992) "Determining the Nuclear Equation of State from Neutron-Star
 * Masses and Radii", Astrophys. J. 398 569.
 * Solves the virial identities to obtain an estimate of the error on the
 * Tolman-Oppenheimer-Volkov stellar structure solver.
 * @param[out] radius The radius of the star in m.
 * @param[out] mass The mass of the star in kg.
 * @param[out] baryon_mass The baryon mass of the star in kg.
 * @param[out] love_number_k2 The k_2 tidal love number of the star.
 * @param[out] love_number_k3 The k_3 tidal love number of the star.
 * @param[out] love_number_k4 The k_4 tidal love number of the star.
 * @param[in] central_pressure_si The central pressure of the star in Pa.
 * @param eos Pointer to the Equation of State structure with multiple parts.
 * @param[in] epsrel The relative error for the TOV solver routine
 * @retval 0 Success.
 * @retval <0 Failure.
 */
int XLALSimNeutronStarTOVODEExtendedVirialIntegrateWithTolerance(double *radius, double *mass, double *baryon_mass,
             double *love_number_k2, double *love_number_k3, double *love_number_k4,
             double *int1, double *int2, double *int3, double *int4, double *int5, double *int6,
             double central_pressure_si,
             struct EOSMultiParts eos,
             double epsrel){

    /* ode integration variables */
    const double epsabs = 0.0;
    double y[TOV_EXT_ODE_VARS_DIM]; // array for the variables of the TOV ODE equations
    double dy[TOV_EXT_ODE_VARS_DIM];
    double y_vir[VIRIAL_ODE_VARS_DIM]; // array for the variables of the VIRIAL ODE equations
    double dy_vir[VIRIAL_ODE_VARS_DIM];
    size_t i;
    int s, s_vir;
    LALSimNeutronStarEOS * params;
    struct keeping_data params_vir;

    struct tov_ext_ode_vars *vars = tov_ext_ode_vars_cast(y);
    struct virial_ode_vars *vars_vir = virial_ode_vars_cast(y_vir);

    /* Set up the actual ODE system to solve with tov_virial_ode */
    gsl_odeiv_system sys = {tov_ext_ode, NULL, TOV_EXT_ODE_VARS_DIM, &params};
    /* Set up the iterative method to solve the ODE, and the dimension of the ODE set */
    gsl_odeiv_step *step = gsl_odeiv_step_alloc(gsl_odeiv_step_rk8pd, TOV_EXT_ODE_VARS_DIM);
    /* Set up the precision of the numerical method with relative error epsrel */
    gsl_odeiv_control *ctrl = gsl_odeiv_control_y_new(epsabs, epsrel);
    /* Set up evolution function to solve ODE with the dimension of the ODE set */
    gsl_odeiv_evolve *evolv = gsl_odeiv_evolve_alloc(TOV_EXT_ODE_VARS_DIM);


    /* Set up the actual ODE system to solve with virial_ode */
    gsl_odeiv_system sys_vir = {virial_ode, NULL, VIRIAL_ODE_VARS_DIM, &params_vir};
    gsl_odeiv_step *step_vir = gsl_odeiv_step_alloc(gsl_odeiv_step_rk8pd, VIRIAL_ODE_VARS_DIM);
    gsl_odeiv_control *ctrl_vir = NULL; // no adaptative step for the virial parameters
    gsl_odeiv_evolve *evolv_vir = gsl_odeiv_evolve_alloc(VIRIAL_ODE_VARS_DIM);

    /* Equation of state */
    int number_eos = eos.number_of_parts;
    double *hmin = eos.hmin;

    /* Central values */
    double pc = central_pressure_si * LAL_G_C4_SI; // convert the pressure from Pascal to geometrized units [/m^2]
    double ec = XLALSimNeutronStarEOSEnergyDensityOfPressureGeometerized(pc, eos.eos_part[number_eos-1]);
    double hc = XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometerized(pc, eos.eos_part[number_eos-1]);


    /* Initial condition */
    double dh = -1e-12 * hc;
    double h0 = hc + dh;
    double h1 = 0.0 - dh;
    double h = h0;
    double h_old = h;
    double dh_old = dh;
    tov_ext_initial_condition(ec, pc, dh, eos.eos_part[number_eos-1], vars);
    virial_initial_condition(vars->r, pc, dh, vars_vir);

    for(int j = number_eos-1; j >= 0; j--) {
        params = eos.eos_part[j];
        double pmin = XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometerized(hmin[j], eos.eos_part[j]);
        while (h > hmin[j]) {
            params_vir.old_m = vars->m;
            params_vir.old_r = vars->r;
            params_vir.eos = eos.eos_part[j];
//             printf("Before virial solver: h = %.6e dh = %.6e\n", h, dh);

            s_vir = gsl_odeiv_evolve_apply(evolv_vir, ctrl_vir, step_vir, &sys_vir, &h_old, hmin[j], &dh_old, y_vir);
//             printf("After virial solver: h = %.6e dh = %.6e\n", h, dh);
            if (s_vir != GSL_SUCCESS)
                XLAL_ERROR(XLAL_EERR,
                    "Error encountered in GSL's ODE integrator\n");

//             printf("Before TOV solver: h = %.6e dh = %.6e\n", h, dh);

            s = gsl_odeiv_evolve_apply(evolv, ctrl, step, &sys, &h, hmin[j], &dh, y);
//             printf("After TOV solver: h = %.6e dh = %.6e\n", h, dh);
            if (s != GSL_SUCCESS)
                XLAL_ERROR(XLAL_EERR,
                    "Error encountered in GSL's ODE integrator\n");
//             printf("\t\t %d Star integration h= %.16e \t M = %.6e \n", j,  h, vars->m  / LAL_MRSUN_SI);
//             printf("\t\t avant %d Star integration h= %.16e \t dh= %.16e \t M = %.6e \t I2 = %.6e\n", j,  h, dh, vars->m  / LAL_MRSUN_SI, vars_vir->I2);


//             printf("\t\t %d Star integration h= %.16e \t dh= %.16e \t M = %.6e \t I2 = %.6e\n\n", j,  h, dh, vars->m  / LAL_MRSUN_SI, vars_vir->I2);

            h_old = h;
            dh_old = dh;
        }

        if (pc >= pmin && j != 0){
//             printf("hmin =%6e \nValues of b before correction %.6e for k2\n", hmin[j], vars->b_k2);
            /* Phase transition correction for the tidal love number
            * Implements Eq.(41) of Pereira et al. 2020 ApJ 895 28
            * Note: Eq.(14) of Postnikov et al. 2010 Phys. Rev. D 82, 024016
            * takes an approximation for the denominator on the correction terms
            */
            double r3 = (vars->r)*(vars->r)*(vars->r);
            double pres_pt = XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometerized(h, eos.eos_part[j]);
            double dpt_eps =  XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometerized(h, eos.eos_part[j])
                            - XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometerized(h, eos.eos_part[j-1]) ;

            double rho_bar = (vars->m + 4.0 * LAL_PI * r3 * pres_pt)/ (4.0 / 3.0 * LAL_PI * r3 ) ;
            vars->b_k2 += - vars->H_k2 / vars->r * dpt_eps * 3.0 / rho_bar ;
            vars->b_k3 += - vars->H_k3 / vars->r * dpt_eps * 3.0 / rho_bar ;
            vars->b_k4 += - vars->H_k4 / vars->r * dpt_eps * 3.0 / rho_bar ;
//             printf("Values of b after correction %.6e for k2\n", vars->b_k2);
//             printf("\t\t %d CORRECTION ON h= %.16e \t M = %.6e \n", j,h, vars->m/LAL_MRSUN_SI);

        }

    }

    /* compute tidal Love numbers and compactness at the surface of the star */
    double comp = vars->m / vars->r;      /* compactness */
    double yy_k2 = vars->r * vars->b_k2 / vars->H_k2; /* Eq. 13 of Hinderer et al. Phys. Rev. D 81 123016 */
    double yy_k3 = vars->r * vars->b_k3 / vars->H_k3; /* Eq. 13 of Hinderer et al. Phys. Rev. D 81 123016 */
    double yy_k4 = vars->r * vars->b_k4 / vars->H_k4; /* Eq. 13 of Hinderer et al. Phys. Rev. D 81 123016 */


    /*take one final Euler step to get to surface*/
    params_vir.old_m = vars->m;
    params_vir.old_r = vars->r;
    params_vir.eos = eos.eos_part[0];

    for (int w = 0 ; w < 1 ; ++w){
        tov_ext_ode(h, y, dy, &eos.eos_part[0]);
        for (i = 0; i < TOV_EXT_ODE_VARS_DIM; ++i) y[i] += dy[i] * (0.0 - h1);
    }

    for (int w = 0 ; w < 1 ; ++w){
        virial_ode(h_old, y_vir, dy_vir, &params_vir);
        for (i = 0; i < VIRIAL_ODE_VARS_DIM; ++i) y_vir[i] += dy_vir[i] * (0.0 - h1);
    }

    /* convert from geometric units to SI units */

    *int3 = (1.0 - vars->m / vars->r) * pow((1.0 - 2.0 * vars->m / vars->r), (-0.5)) - 1.0;
    *int6 = vars->r * (*int3);

    *int1 = vars_vir->I1;
    *int2 = vars_vir->I2;
    *int4 = vars_vir->J1;
    *int5 = vars_vir->J2;


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

    return 0;


}


//CUTER-dev
/**
 * @brief Integrates the Tolman-Oppenheimer-Volkov stellar structure equations
 * on a predefined grid of enthalpy with no adaptative step solver.
 * @details
 * Solves the Tolman-Oppenheimer-Volkov stellar structure equations using the
 * pseudo-enthalpy formalism introduced in:
 * Lindblom (1992) "Determining the Nuclear Equation of State from Neutron-Star
 * Masses and Radii", Astrophys. J. 398 569.
 * @param[out] radius The radius of the star in m.
 * @param[out] mass The mass of the star in kg.
 * @param[out] baryon_mass The baryon mass of the star in kg.
 * @param[out] love_number_k2 The k_2 tidal love number of the star.
 * @param[out] love_number_k3 The k_3 tidal love number of the star.
 * @param[out] love_number_k4 The k_4 tidal love number of the star.
 * @param[in] central_pressure_si The central pressure of the star in Pa.
 * @param eos Pointer to the Equation of State structure with multiple parts.
 * @param[in] epsrel The relative error for the TOV solver routine
 * @retval 0 Success.
 * @retval <0 Failure.
 */
int XLALSimNeutronStarTOVODEExtendedGridIntegrateWithTolerance(double *radius, double *mass, double *baryon_mass,
             double *love_number_k2, double *love_number_k3, double *love_number_k4,
             double central_pressure_si,
             struct EOSMultiParts eos){

    /* ode integration variables */
//     const double epsabs = 0.0;
    double y[TOV_EXT_ODE_VARS_DIM]; // array for the variables of the ODE equations
    double yerr[TOV_EXT_ODE_VARS_DIM];
    double dy[TOV_EXT_ODE_VARS_DIM];
    size_t i;
    int s;
    LALSimNeutronStarEOS * params;


    struct tov_ext_ode_vars *vars = tov_ext_ode_vars_cast(y);
    /* Set up the actual ODE system to solve with tov_virial_ode */
    gsl_odeiv_system sys = {tov_ext_ode, NULL, TOV_EXT_ODE_VARS_DIM, &params};
    /* Set up the iterative method to solve the ODE, and the dimension of the ODE set */
    gsl_odeiv_step *step = gsl_odeiv_step_alloc(gsl_odeiv_step_rk8pd, TOV_EXT_ODE_VARS_DIM);

    /* Equation of state */
    int number_eos = eos.number_of_parts;
    double *hmin = eos.hmin;
    /* Central values */
    double pc = central_pressure_si * LAL_G_C4_SI; // convert the pressure from Pascal to geometrized units [/m^2]
    double ec = XLALSimNeutronStarEOSEnergyDensityOfPressureGeometerized(pc, eos.eos_part[number_eos-1]);
    double hc = XLALSimNeutronStarEOSPseudoEnthalpyOfPressureGeometerized(pc, eos.eos_part[number_eos-1]);



    int npt_dh1 = 400;
    int npt_dh2 = 10;
    double min_dh1 = 1.e-13;
    double max_dh1 = 1.e-5;
    double min_dh2 = 1.e-5;
    double max_dh2 = 1.e-10;
    double dh_grid[npt_dh1 + npt_dh2];
//
    double step_logdh1 = (log10(max_dh1)-log10(min_dh1))/npt_dh1;
    double step_logdh2 = (log10(max_dh2)-log10(min_dh2))/npt_dh2;
//     for (int ii = 0; ii < 9 ; ii++) h_grid[ii] = -1.e-12 * hc;
    for (int ii = 0; ii < npt_dh1 ; ii++) {
        dh_grid[ii] = - min_dh1 * pow(10.0, (ii+1)*step_logdh1);
        printf("here 1 %d %.6e \n", ii, dh_grid[ii]);
    }

    for (int ii = npt_dh1; ii < npt_dh2 + npt_dh1 ; ii++) {
        dh_grid[ii] = - min_dh2 * pow(10.0, (ii+1 - npt_dh1)*step_logdh2);
        printf("here %d %.6e \n", ii, dh_grid[ii]);
    }

    /* Initial condition */
    double dh = -1e-12 * hc;
    printf("blabla %.6e %.6e\n", dh, hc);
//     dh_grid[0] = dh;
    double h0 = hc + dh;
    double h1 = 0.0 - dh;
    double h;


    tov_ext_initial_condition(ec, pc, dh, eos.eos_part[number_eos-1], vars);
    h = h0;
    int count = 0;
    for(int j = number_eos-1; j >= 0; j--) {
        params = eos.eos_part[j];
        double pmin = XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometerized(hmin[j], eos.eos_part[j]);

        while (h > hmin[j]){// && (vars->m - m_old)/m_old > 1e-2){//&& count < (npt_dh1 + npt_dh2)) {
//             printf("here %.6e %.6e %.6e\n", h, h_grid[count], dh);
            s = gsl_odeiv_step_apply(step, h, dh, y, yerr, NULL, NULL, &sys);
//             printf("after %.6e %.6e\n", h, dh);
            if (s != GSL_SUCCESS)
                XLAL_ERROR(XLAL_EERR,
                    "Error encountered in GSL's ODE integrator\n");
            printf("\t\t %d %d Star integration h= %.16e \t dh= %.16e \t M = %.6e \n", count, j,  h, dh, vars->m  / LAL_MRSUN_SI);

            dh = -1.e-5;//dh_grid[count];
            h += dh; //h_grid[count];
//             m_old = vars->m;
            count+=1;
        }

        if (pc >= pmin && j != 0){
//             printf("hmin =%6e \nValues of b before correction %.6e for k2\n", hmin[j], vars->b_k2);
            /* Phase transition correction for the tidal love number
            * Implements Eq.(41) of Pereira et al. 2020 ApJ 895 28
            * Note: Eq.(14) of Postnikov et al. 2010 Phys. Rev. D 82, 024016
            * takes an approximation for the denominator on the correction terms
            */
            double r3 = (vars->r)*(vars->r)*(vars->r);
            double pres_pt = XLALSimNeutronStarEOSPressureOfPseudoEnthalpyGeometerized(h, eos.eos_part[j]);
            double dpt_eps =  XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometerized(h, eos.eos_part[j])
                            - XLALSimNeutronStarEOSEnergyDensityOfPseudoEnthalpyGeometerized(h, eos.eos_part[j-1]) ;

            double rho_bar = (vars->m + 4.0 * LAL_PI * r3 * pres_pt)/ (4.0 / 3.0 * LAL_PI * r3 ) ;
            vars->b_k2 += - vars->H_k2 / vars->r * dpt_eps * 3.0 / rho_bar ;
            vars->b_k3 += - vars->H_k3 / vars->r * dpt_eps * 3.0 / rho_bar ;
            vars->b_k4 += - vars->H_k4 / vars->r * dpt_eps * 3.0 / rho_bar ;
//             printf("Values of b after correction %.6e for k2\n", vars->b_k2);
//             printf("\t\t %d CORRECTION ON h= %.16e \t M = %.6e \n", j,h, vars->m/LAL_MRSUN_SI);

        }


    }

    /* compute tidal Love numbers and compactness at the surface of the star */
    double comp = vars->m / vars->r;      /* compactness */
    double yy_k2 = vars->r * vars->b_k2 / vars->H_k2; /* Eq. 13 of Hinderer et al. Phys. Rev. D 81 123016 */
    double yy_k3 = vars->r * vars->b_k3 / vars->H_k3; /* Eq. 13 of Hinderer et al. Phys. Rev. D 81 123016 */
    double yy_k4 = vars->r * vars->b_k4 / vars->H_k4; /* Eq. 13 of Hinderer et al. Phys. Rev. D 81 123016 */

    /*take one final Euler step to get to surface*/
    for (int w = 0 ; w < 1 ; ++w){
        tov_ext_ode(h, y, dy, &eos.eos_part[0]);
        for (i = 0; i < TOV_EXT_ODE_VARS_DIM; ++i) // No need for all Virial variables, only physically relevant ones
            y[i] += dy[i] * (0.0 - h1);
    }

    /* convert from geometric units to SI units */
    *radius = vars->r;
    *mass = vars->m * LAL_MSUN_SI / LAL_MRSUN_SI;
    *baryon_mass = vars->mb * LAL_MSUN_SI / LAL_MRSUN_SI;
    *love_number_k2 = tidal_Love_number_k2(comp, yy_k2);
    *love_number_k3 = tidal_Love_number_k3(comp, yy_k3);
    *love_number_k4 = tidal_Love_number_k4(comp, yy_k4);


    /* free ode memory */
    gsl_odeiv_step_free(step);

    return 0;


}


int XLALSimNeutronStarTOVODEIntegrate(double *radius, double *mass,
    double *love_number_k2, double central_pressure_si, LALSimNeutronStarEOS * eos)
{
    const double epsrel = 1e-6;
    return XLALSimNeutronStarTOVODEIntegrateWithTolerance(radius, mass, love_number_k2, central_pressure_si, eos, epsrel);
}


int XLALSimNeutronStarVirialODEIntegrate(double *radius, double *mass,
    double *int1, double *int2, double *int3, double *int4, double *int5, double *int6,
    double *love_number_k2, double central_pressure_si,
    LALSimNeutronStarEOS * eos)
{
    const double epsrel = 1e-6;
    return XLALSimNeutronStarVirialODEIntegrateWithTolerance(radius, mass,
    int1, int2, int3, int4, int5, int6, love_number_k2, central_pressure_si, eos, epsrel);
}


/** @} */
