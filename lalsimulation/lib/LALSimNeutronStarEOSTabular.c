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
 * @author Philip Davis, Micaela Oertel, Lami Suleiman, Leslie Wade
 * @brief Provides routines for tabulated equations of state
 * @addtogroup LALSimNeutronStarEOS_c
 * @{
 */
 /**
 * @name Creation routines for tabulated equations of state
 * @{
 */


#include <stdbool.h>
#include <stdio.h>
#include <lal/LALSimReadData.h>
#include <gsl/gsl_interp.h>
#include <lal/LALSimNeutronStar.h>
GSL_VAR const gsl_interp_type * lal_gsl_interp_steffen;

/** @cond */

/* Contents of the tabular equation of state data structure. */
struct tagLALSimNeutronStarEOSDataTabular {
    double *nbdat;
    double *log_edat;
    double *log_pdat;
    double *mubdat;
    double *muedat;
    double *log_hdat;
    double *yedat;
    double *log_cs2dat;
    double *log_rhodat;
    size_t ncol;
    size_t ndat;

    gsl_interp *log_e_of_log_p_interp;
    gsl_interp *log_h_of_log_p_interp;
    gsl_interp *log_e_of_log_h_interp;
    gsl_interp *log_p_of_log_h_interp;
    gsl_interp *log_rho_of_log_h_interp;
    gsl_interp *log_p_of_log_e_interp;
    gsl_interp *log_p_of_log_rho_interp;
    gsl_interp *log_cs2_of_log_h_interp;
    gsl_interp_accel *log_e_of_log_p_acc;
    gsl_interp_accel *log_h_of_log_p_acc;
    gsl_interp_accel *log_e_of_log_h_acc;
    gsl_interp_accel *log_p_of_log_h_acc;
    gsl_interp_accel *log_rho_of_log_h_acc;
    gsl_interp_accel *log_p_of_log_e_acc;
    gsl_interp_accel *log_p_of_log_rho_acc;
    gsl_interp_accel *log_cs2_of_log_h_acc;
};

/* Clamping function. */
static inline double clamp_to_range_tol(double x, double xmin, double xmax, double tol)
{
    if (x < xmin && xmin-x < tol )
        return xmin;
    else if (x > xmax && x-xmax < tol)
        return xmax;
    else
        return x;
}

static double eos_p_of_e_tabular(double e, LALSimNeutronStarEOS * eos)
{
    double log_e;
    double log_p;
    double tol = 1e-12;
    if (e == 0.0)
	return 0.0;
    log_e = log(e);
    if (log_e > eos->data.tabular->log_edat[eos->data.tabular->ndat-1] + tol)
        XLAL_ERROR_REAL8(XLAL_EDOM,
            "Energy density e=%.5e is above the EOS interpolation range.", e);
    // Clamp to interpolation range within tolerance to handle floating-point roundoff
    log_e = clamp_to_range_tol(log_e, eos->data.tabular->log_edat[0],
        eos->data.tabular->log_edat[eos->data.tabular->ndat-1], tol);
    if (log_e < eos->data.tabular->log_edat[0])
	/* use non-relativistic degenerate gas, p = K * e**(5./3.) */
	return exp(eos->data.tabular->log_pdat[0] + (5.0 / 3.0) * (log_e - eos->data.tabular->log_edat[0]));
    log_p = gsl_interp_eval(eos->data.tabular->log_p_of_log_e_interp,
        eos->data.tabular->log_edat, eos->data.tabular->log_pdat, log_e,
        eos->data.tabular->log_p_of_log_e_acc);
    return exp(log_p);
}

static double eos_p_of_rho_tabular(double rho, LALSimNeutronStarEOS * eos)
{
    double log_rho;
    double log_p;
    double tol = 1e-12;
    if (rho == 0.0)
	return 0.0;
    log_rho = log(rho);
    if (log_rho > eos->data.tabular->log_rhodat[eos->data.tabular->ndat-1] + tol)
        XLAL_ERROR_REAL8(XLAL_EDOM,
            "Rest-mass density rho=%.5e is above the EOS interpolation range.", rho);
    // Clamp to interpolation range within tolerance to handle floating-point roundoff
    log_rho = clamp_to_range_tol(log_rho, eos->data.tabular->log_rhodat[0],
        eos->data.tabular->log_rhodat[eos->data.tabular->ndat-1], tol);
    if (log_rho < eos->data.tabular->log_rhodat[0])
	/* use non-relativistic degenerate gas, p = K * rho**(5./3.) */
	return exp(eos->data.tabular->log_pdat[0] + (5.0 / 3.0) * (log_rho - eos->data.tabular->log_rhodat[0]));
    log_p = gsl_interp_eval(eos->data.tabular->log_p_of_log_rho_interp,
        eos->data.tabular->log_rhodat, eos->data.tabular->log_pdat, log_rho,
        eos->data.tabular->log_p_of_log_rho_acc);
    return exp(log_p);
}

static double eos_e_of_p_tabular(double p, LALSimNeutronStarEOS * eos)
{
    double log_p;
    double log_e;
    double tol = 1e-12;
    if (p == 0.0)
	return 0.0;
    log_p = log(p);
    if (log_p > eos->data.tabular->log_pdat[eos->data.tabular->ndat-1] + tol)
        XLAL_ERROR_REAL8(XLAL_EDOM,
            "Pressure p=%.5e is above the EOS interpolation range.", p);
    // Clamp to interpolation range within tolerance to handle floating-point roundoff
    log_p = clamp_to_range_tol(log_p, eos->data.tabular->log_pdat[0],
        eos->data.tabular->log_pdat[eos->data.tabular->ndat-1], tol);
    if (log_p < eos->data.tabular->log_pdat[0])
	/* use non-relativistic degenerate gas, p = K * e**(5./3.) */
	return exp(eos->data.tabular->log_edat[0] + (3.0 / 5.0) * (log_p - eos->data.tabular->log_pdat[0]));
    log_e = gsl_interp_eval(eos->data.tabular->log_e_of_log_p_interp,
        eos->data.tabular->log_pdat, eos->data.tabular->log_edat, log_p,
        eos->data.tabular->log_e_of_log_p_acc);
    return exp(log_e);
}

static double eos_e_of_h_tabular(double h, LALSimNeutronStarEOS * eos)
{
    double log_h;
    double log_e;
    double tol = 1e-12;
    if (h == 0.0)
	return 0.0;
    log_h = log(h);
    if (log_h > eos->data.tabular->log_hdat[eos->data.tabular->ndat-1] + tol)
        XLAL_ERROR_REAL8(XLAL_EDOM,
            "Pseudo-enthalpy h=%.5e is above the EOS interpolation range.", h);
    // Clamp to interpolation range within tolerance to handle floating-point roundoff
    log_h = clamp_to_range_tol(log_h, eos->data.tabular->log_hdat[0],
        eos->data.tabular->log_hdat[eos->data.tabular->ndat-1], tol);
    if (log_h < eos->data.tabular->log_hdat[0])
	/* use non-relativistic degenerate gas, e = K * h**(3./2.) */
	return exp(eos->data.tabular->log_edat[0] + 1.5 * (log_h - eos->data.tabular->log_hdat[0]));
    log_e = gsl_interp_eval(eos->data.tabular->log_e_of_log_h_interp,
        eos->data.tabular->log_hdat, eos->data.tabular->log_edat, log_h,
        eos->data.tabular->log_e_of_log_h_acc);
    return exp(log_e);
}

static double eos_p_of_h_tabular(double h, LALSimNeutronStarEOS * eos)
{
    double log_h;
    double log_p;
    double tol = 1e-12;
    if (h == 0.0)
	return 0.0;
    log_h = log(h);
    if (log_h > eos->data.tabular->log_hdat[eos->data.tabular->ndat-1] + tol)
        XLAL_ERROR_REAL8(XLAL_EDOM,
            "Pseudo-enthalpy h=%.5e is above the EOS interpolation range.", h);
    // Clamp to interpolation range within tolerance to handle floating-point roundoff
    log_h = clamp_to_range_tol(log_h, eos->data.tabular->log_hdat[0],
        eos->data.tabular->log_hdat[eos->data.tabular->ndat-1], tol);
    if (log_h < eos->data.tabular->log_hdat[0])
	/* use non-relativistic degenerate gas, p = K * h**(5./2.) */
	return exp(eos->data.tabular->log_pdat[0] + 2.5 * (log_h - eos->data.tabular->log_hdat[0]));
    log_p = gsl_interp_eval(eos->data.tabular->log_p_of_log_h_interp,
        eos->data.tabular->log_hdat, eos->data.tabular->log_pdat, log_h,
        eos->data.tabular->log_p_of_log_h_acc);
    return exp(log_p);
}

static double eos_rho_of_h_tabular(double h, LALSimNeutronStarEOS * eos)
{
    double log_h;
    double log_rho;
    double tol = 1e-12;
    if (h == 0.0)
	return 0.0;
    log_h = log(h);
    if (log_h > eos->data.tabular->log_hdat[eos->data.tabular->ndat-1] + tol)
        XLAL_ERROR_REAL8(XLAL_EDOM,
            "Pseudo-enthalpy h=%.5e is above the EOS interpolation range.", h);
    // Clamp to interpolation range within tolerance to handle floating-point roundoff
    log_h = clamp_to_range_tol(log_h, eos->data.tabular->log_hdat[0],
        eos->data.tabular->log_hdat[eos->data.tabular->ndat-1], tol);
    if (log_h < eos->data.tabular->log_hdat[0])
	/* use non-relativistic degenerate gas, rho = K * h**(3./2.) */
	return exp(eos->data.tabular->log_rhodat[0] + 1.5 * (log_h - eos->data.tabular->log_hdat[0]));
    log_rho = gsl_interp_eval(eos->data.tabular->log_rho_of_log_h_interp,
        eos->data.tabular->log_hdat, eos->data.tabular->log_rhodat, log_h,
        eos->data.tabular->log_rho_of_log_h_acc);
    return exp(log_rho);
}

static double eos_h_of_p_tabular(double p, LALSimNeutronStarEOS * eos)
{
    double log_p;
    double log_h;
    double tol = 1e-12;
    if (p == 0)
	return 0.0;
    log_p = log(p);
    if (log_p > eos->data.tabular->log_pdat[eos->data.tabular->ndat-1] + tol)
        XLAL_ERROR_REAL8(XLAL_EDOM,
            "Pressure p=%.5e is above the EOS interpolation range.", p);
    // Clamp to interpolation range within tolerance to handle floating-point roundoff
    log_p = clamp_to_range_tol(log_p, eos->data.tabular->log_pdat[0],
        eos->data.tabular->log_pdat[eos->data.tabular->ndat-1], tol);
    if (log_p < eos->data.tabular->log_pdat[0])
	/* use non-relativistic degenerate gas, h = K * p**(2./5.) */
	return exp(eos->data.tabular->log_hdat[0] + 0.4 * (log_p - eos->data.tabular->log_pdat[0]));
    log_h = gsl_interp_eval(eos->data.tabular->log_h_of_log_p_interp,
        eos->data.tabular->log_pdat, eos->data.tabular->log_hdat, log_p,
        eos->data.tabular->log_h_of_log_p_acc);
    return exp(log_h);
}

static double eos_dedp_of_p_tabular(double p, LALSimNeutronStarEOS * eos)
{
    double log_p;
    double log_e;
    double d_log_e_d_log_p;
    double tol = 1e-12;
    if (p == 0 || (log_p = log(p)) < eos->data.tabular->log_pdat[0])
	/* use non-relativistic degenerate gas, p = K * e**(5./3.) */
	return (3.0 / 5.0) * exp(eos->data.tabular->log_edat[0] - eos->data.tabular->log_pdat[0]);
    if (log_p > eos->data.tabular->log_pdat[eos->data.tabular->ndat-1] + tol)
        XLAL_ERROR_REAL8(XLAL_EDOM,
            "Pressure p=%.5e is above the EOS interpolation range.", p);
    // Clamp to interpolation range within tolerance to handle floating-point roundoff
    log_p = clamp_to_range_tol(log_p, eos->data.tabular->log_pdat[0],
        eos->data.tabular->log_pdat[eos->data.tabular->ndat-1], tol);
    log_e = gsl_interp_eval(eos->data.tabular->log_e_of_log_p_interp,
        eos->data.tabular->log_pdat, eos->data.tabular->log_edat, log_p,
        eos->data.tabular->log_e_of_log_p_acc);
    d_log_e_d_log_p = gsl_interp_eval_deriv(eos->data.tabular->log_e_of_log_p_interp,
        eos->data.tabular->log_pdat, eos->data.tabular->log_edat, log_p,
        eos->data.tabular->log_e_of_log_p_acc);
    return d_log_e_d_log_p * exp(log_e - log_p);
}

static double eos_v_of_h_tabular(double h, LALSimNeutronStarEOS * eos)
{
    double p, dedp, log_cs2, log_h;
    double tol = 1e-12;
    if (eos->data.tabular->ncol == 2)
    {
        p = eos_p_of_h_tabular(h, eos);
        dedp = eos_dedp_of_p_tabular(p, eos);
        return pow(dedp, -0.5);
    }
    log_h = log(h);
    if (log_h < eos->data.tabular->log_hdat[0] - tol ||
        log_h > eos->data.tabular->log_hdat[eos->data.tabular->ndat-1] + tol)
        XLAL_ERROR_REAL8(XLAL_EDOM,
            "Pseudo-enthalpy h=%.5e is outside the EOS interpolation range.", h);
    // Clamp to interpolation range within tolerance to handle floating-point roundoff
    log_h = clamp_to_range_tol(log_h, eos->data.tabular->log_hdat[0],
        eos->data.tabular->log_hdat[eos->data.tabular->ndat-1], tol);
    log_cs2 = gsl_interp_eval(eos->data.tabular->log_cs2_of_log_h_interp,
    eos->data.tabular->log_hdat, eos->data.tabular->log_cs2dat, log_h,
    eos->data.tabular->log_cs2_of_log_h_acc);

    return pow(exp(log_cs2), 0.5);
}

//static double eos_v_of_h_tabular(double h, LALSimNeutronStarEOS *eos)
//{
//      double dpdh, dedh;
//    printf("hi4\n");
//    dpdh = gsl_interp_eval_deriv(eos->data.tabular->log_p_of_log_h_interp, eos->data.tabular->hdat, eos->data.tabular->pdat, h, eos->data.tabular->p_of_h_acc);
//    printf("hi5\n");
//    dedh = gsl_interp_eval_deriv(eos->data.tabular->e_of_h_interp, eos->data.tabular->hdat, eos->data.tabular->edat, h, eos->data.tabular->e_of_h_acc);
//      return sqrt(dpdh/dedh);
//}

static void eos_free_tabular_data(LALSimNeutronStarEOSDataTabular * data)
{
    if (!data) return;

    gsl_interp_free(data->log_e_of_log_p_interp);
    gsl_interp_free(data->log_e_of_log_h_interp);
    gsl_interp_free(data->log_p_of_log_h_interp);
    gsl_interp_free(data->log_h_of_log_p_interp);
    gsl_interp_free(data->log_rho_of_log_h_interp);
    gsl_interp_free(data->log_p_of_log_e_interp);
    gsl_interp_free(data->log_p_of_log_rho_interp);
    gsl_interp_free(data->log_cs2_of_log_h_interp);

    gsl_interp_accel_free(data->log_e_of_log_p_acc);
    gsl_interp_accel_free(data->log_e_of_log_h_acc);
    gsl_interp_accel_free(data->log_p_of_log_h_acc);
    gsl_interp_accel_free(data->log_h_of_log_p_acc);
    gsl_interp_accel_free(data->log_rho_of_log_h_acc);
    gsl_interp_accel_free(data->log_p_of_log_e_acc);
    gsl_interp_accel_free(data->log_p_of_log_rho_acc);
    gsl_interp_accel_free(data->log_cs2_of_log_h_acc);

    LALFree(data->nbdat);
    LALFree(data->log_edat);
    LALFree(data->log_pdat);
    LALFree(data->mubdat);
    LALFree(data->muedat);
    LALFree(data->log_hdat);
    LALFree(data->yedat);
    LALFree(data->log_cs2dat);
    LALFree(data->log_rhodat);

    LALFree(data);
}

static void eos_free_tabular(LALSimNeutronStarEOS * eos)
{
    if (!eos) return;

    if (eos->data.tabular) {
        eos_free_tabular_data(eos->data.tabular);
        eos->data.tabular = NULL;
    }

    LALFree(eos);
}

/* Finding density where EOS becomes acausal */

/* Evaluate vSound at each tabulated point until vSound>1 or you get to last
 * point.  If vSound>1 interpolate between that point and previous point to
 * find h where EOS first becomes acausal.  (You could do root-finding to be
 * even more precise, but the calculation of v from the tabulated EOS isn't
 * that exact anyway.) */

/* Minimum pseudo-enthalpy at which EOS becomes acausal (speed of sound > 1).
 * If the EOS is always causal, return some large value hmax instead. */
static double eos_min_acausal_pseudo_enthalpy_tabular(double hmax,
    LALSimNeutronStarEOS * eos)
{
    size_t i;
    double h_im1, h_i;
    double v_im1, v_i;
    double m;   /* slope for linear interpolation */
    double hMinAcausal = hmax;  /* default large number for EOS that is always causal */

    h_im1 = exp(eos->data.tabular->log_hdat[0]);
    v_im1 = eos_v_of_h_tabular(h_im1, eos);
    for (i = 1; i < eos->data.tabular->ndat; i++) {
        h_i = exp(eos->data.tabular->log_hdat[i]);
        v_i = eos_v_of_h_tabular(h_i, eos);
        if (v_i > 1.0) {
            /* solve vsound(h) = 1 */
            m = (v_i - v_im1) / (h_i - h_im1);
            hMinAcausal = h_im1 + (1.0 - v_im1) / m;
            break;
        }
        h_im1 = h_i;
        v_im1 = v_i;
    }
//    printf("hMinAcausal = %e, v = %e\n", hMinAcausal, eos_v_of_h_tabular(hMinAcausal, eos));

    /* Value of h where EOS first becomes acausal.
     * Or, if EOS is always causal, hmax */
    return hMinAcausal;
}

static LALSimNeutronStarEOS *eos_alloc_tabular(double *nbdat, double *edat, double *pdat,
   double *mubdat, double *muedat, double *hdat, double *yedat, double *cs2dat, size_t ndat, size_t ncol)
{
    LALSimNeutronStarEOS *eos;
    LALSimNeutronStarEOSDataTabular *data;
    size_t i;

    eos = LALCalloc(1, sizeof(*eos));
    data = LALCalloc(1, sizeof(*data));

    eos->datatype = LALSIM_NEUTRON_STAR_EOS_DATA_TYPE_TABULAR;
    eos->data.tabular = data;

    /* setup function pointers */
    eos->free = eos_free_tabular;
    eos->e_of_p = eos_e_of_p_tabular;
    eos->h_of_p = eos_h_of_p_tabular;
    eos->e_of_h = eos_e_of_h_tabular;
    eos->p_of_h = eos_p_of_h_tabular;
    eos->rho_of_h = eos_rho_of_h_tabular;
    eos->p_of_e = eos_p_of_e_tabular;
    eos->p_of_rho = eos_p_of_rho_tabular;
    eos->dedp_of_p = eos_dedp_of_p_tabular;
    eos->v_of_h = eos_v_of_h_tabular;

    data->log_rhodat = XLALMalloc(ndat * sizeof(*data->log_rhodat));

    if(ncol == 2) {
        /* allocate memory for eos data; ignore first points if 0 */
        while (*pdat == 0.0 || *edat == 0.0) {
            ++pdat;
            ++edat;
            --ndat;
        }

        data->ncol = ncol;
        data->ndat = ndat;
        data->log_pdat = XLALMalloc(ndat * sizeof(*data->log_pdat));
        data->log_edat = XLALMalloc(ndat * sizeof(*data->log_edat));
        data->log_hdat = XLALMalloc(ndat * sizeof(*data->log_hdat));

        /* take log of eos data */
        for (i = 0; i < ndat; ++i) {
            data->log_pdat[i] = log(pdat[i]);
            data->log_edat[i] = log(edat[i]);
        }
        /* compute pseudo-enthalpy h from dhdp */
        /* Integrate in log space:
        dhdp = 1 / [e(p) + p]
        h(p) = h(p0) + \int_p0^p dhdp dp
        h(p) = h(p0) + \int_ln(p0)^ln(p) exp[ln(p) + ln(dhdp)] dln(p)
        First point is
        h(p0) = p0 / [e(p0) + p0]
        */
        double *integrand;
        integrand = LALMalloc(ndat * sizeof(*integrand));
        for (i = 0; i < ndat; ++i)
            integrand[i] = exp(data->log_pdat[i] + log(1.0 / (edat[i] + pdat[i])));

        gsl_interp_accel * dhdp_of_p_acc_temp = gsl_interp_accel_alloc();
        gsl_interp * dhdp_of_p_interp_temp = gsl_interp_alloc(gsl_interp_linear, ndat);
        gsl_interp_init(dhdp_of_p_interp_temp, data->log_pdat, integrand, ndat);

        data->log_hdat[0] = log(pdat[0] / (edat[0] + pdat[0]));
        for (i = 1; i < ndat; ++i)
            data->log_hdat[i] = log(exp(data->log_hdat[0]) + gsl_interp_eval_integ(dhdp_of_p_interp_temp, data->log_pdat, integrand, data->log_pdat[0], data->log_pdat[i], dhdp_of_p_acc_temp));

        gsl_interp_free(dhdp_of_p_interp_temp);
        gsl_interp_accel_free(dhdp_of_p_acc_temp);
        LALFree(integrand);
    }
    else if (ncol > 2) {
        /* allocate memory for eos data; ignore first points if 0 */
        while (*pdat == 0.0 || *edat == 0.0 || *hdat == 0.0) {
            ++pdat;
            ++edat;
            ++hdat;
            --ndat;
        }

        data->ndat = ndat;
        data->ncol = ncol - 1;
        data->nbdat = XLALMalloc(ndat * sizeof(*data->nbdat));
        data->log_pdat = XLALMalloc(ndat * sizeof(*data->log_pdat));
        data->log_edat = XLALMalloc(ndat * sizeof(*data->log_edat));
        data->mubdat = XLALMalloc(ndat * sizeof(*data->mubdat));
        data->muedat = XLALMalloc(ndat * sizeof(*data->muedat));
        data->log_hdat = XLALMalloc(ndat * sizeof(*data->log_hdat));
        data->yedat = XLALMalloc(ndat * sizeof(*data->yedat));
        data->log_cs2dat = XLALMalloc(ndat * sizeof(*data->log_cs2dat));

        /* take log of eos data */
        for (i = 0; i < ndat; ++i) {
            data->nbdat[i] = nbdat[i];
            data->log_pdat[i] = log(pdat[i]);
            data->log_edat[i] = log(edat[i]);
            data->mubdat[i] = mubdat[i];
            data->muedat[i] = muedat[i];
            data->log_hdat[i] = log(hdat[i]);
            data->yedat[i] = yedat[i];
            data->log_cs2dat[i] = log(cs2dat[i]);
        }

        /* these can be set up only using the new eos tables; so they are set up here */
        data->log_cs2_of_log_h_acc = gsl_interp_accel_alloc();
        data->log_cs2_of_log_h_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat);
        gsl_interp_init(data->log_cs2_of_log_h_interp, data->log_hdat, data->log_cs2dat, ndat);
    }

    // Find rho from e, p, and h: rho = (e+p)/exp(h),
    // definition of the enthalpy for a cold perfect fluid + first law of thermodynamics,
    // See Shapiro & Teukolsky book (chapter 2), or Eqs(7-9) in Haensel & Potekhin 2004.
    for (i = 0; i < ndat; i++)
        data->log_rhodat[i] = log(edat[i] + pdat[i]) - exp(data->log_hdat[i]);

    eos->pmax = exp(data->log_pdat[ndat - 1]);
    eos->hmax = exp(data->log_hdat[ndat - 1]);
    eos->hmin = exp(data->log_hdat[0]);

    /* setup interpolation tables */

    data->log_e_of_log_p_acc = gsl_interp_accel_alloc();
    data->log_h_of_log_p_acc = gsl_interp_accel_alloc();
    data->log_e_of_log_h_acc = gsl_interp_accel_alloc();
    data->log_p_of_log_h_acc = gsl_interp_accel_alloc();
    data->log_rho_of_log_h_acc = gsl_interp_accel_alloc();
    data->log_p_of_log_e_acc = gsl_interp_accel_alloc();
    data->log_p_of_log_rho_acc = gsl_interp_accel_alloc();

    data->log_e_of_log_p_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat);
    data->log_h_of_log_p_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat);
    data->log_e_of_log_h_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat);
    data->log_p_of_log_h_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat);
    data->log_rho_of_log_h_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat);
    data->log_p_of_log_e_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat);
    data->log_p_of_log_rho_interp = gsl_interp_alloc(lal_gsl_interp_steffen, ndat);

    gsl_interp_init(data->log_e_of_log_p_interp, data->log_pdat, data->log_edat, ndat);
    gsl_interp_init(data->log_h_of_log_p_interp, data->log_pdat, data->log_hdat, ndat);
    gsl_interp_init(data->log_e_of_log_h_interp, data->log_hdat, data->log_edat, ndat);
    gsl_interp_init(data->log_p_of_log_h_interp, data->log_hdat, data->log_pdat, ndat);
    gsl_interp_init(data->log_rho_of_log_h_interp, data->log_hdat, data->log_rhodat, ndat);
    gsl_interp_init(data->log_p_of_log_e_interp, data->log_edat, data->log_pdat, ndat);
    gsl_interp_init(data->log_p_of_log_rho_interp, data->log_rhodat, data->log_pdat, ndat);

    eos->hMinAcausal =
        eos_min_acausal_pseudo_enthalpy_tabular(eos->hmax, eos);

//    printf("%e\n", XLALSimNeutronStarEOSEnergyDensityOfPressureGeometrized(eos->pmax, eos));
//
//    printf("datatype = %d\n", eos->datatype);
//    printf("pmax = %e\n", eos->pmax);
//    printf("hmax = %e\n", eos->hmax);
//    printf("hMinAcausal = %e\n", eos->hMinAcausal);

    return eos;
}

static void eos_multi_part_free_tabular(EOSMultiParts * eos)
{
    if (!eos) return;

    if (eos->eos_piece) {
        for (int i = 0; i < eos->number_of_parts; i++) {

            if (eos->eos_piece[i]) {

                // ALWAYS use the EOS own destructor
                if (eos->eos_piece[i]->free) {
                    eos->eos_piece[i]->free(eos->eos_piece[i]);
                } else {
                    // fallback only (should rarely happen)
                    LALFree(eos->eos_piece[i]);
                }

                eos->eos_piece[i] = NULL;
            }
        }

        LALFree(eos->eos_piece);
        eos->eos_piece = NULL;
    }

    LALFree(eos);
}

/* Minimum pseudo-enthalpy at which EOS becomes acausal (speed of sound > 1).
 * If the EOS is always causal, return some large value hmax instead. */
static double eosMultiParts_min_acausal_pseudo_enthalpy_tabular(double hmax,
    EOSMultiParts * eos)
{
    size_t i;
    double h_im1, h_i;
    double v_im1, v_i;
    double m;   /* slope for linear interpolation */
    double hMinAcausal = hmax;  /* default large number for EOS that is always causal */
    int number_of_parts = XLALSimNeutronStarEOSMultiPartsNumber(eos);
    LALSimNeutronStarEOS * eos_piece = XLALSimNeutronStarEOSPart(eos, number_of_parts-1);
    h_im1 = exp(eos_piece ->data.tabular->log_hdat[0]);
    v_im1 = eos_v_of_h_tabular(h_im1, eos_piece);
    for (i = 1; i < eos_piece ->data.tabular->ndat; i++) {
        h_i = exp(eos_piece ->data.tabular->log_hdat[i]);
        v_i = eos_v_of_h_tabular(h_i, eos_piece);
        if (v_i > 1.0) {
            /* solve vsound(h) = 1 */
            m = (v_i - v_im1) / (h_i - h_im1);
            hMinAcausal = h_im1 + (1.0 - v_im1) / m;
            break;
        }
        h_im1 = h_i;
        v_im1 = v_i;
    }
    return hMinAcausal;
}


/* This function corrects "dirty" phase transitions which are defined
 * by Delta P != 0 and Delta h != 0 exactly but numerically zero. The enthalpy
 * is used for the correction: we linearly extrapolate h(P) on both sides of the phase
 * transition, and find the crossing point to determine h and P at the transition.
 * Then, the energy density values on both sides of the phase transition are
 * recalculated using a linear extrapolation as well.
 */
static void eos_correct_phase_transition(double *edat, double *pdat, double *hdat, int index_pt){

    double slope_low, slope_high, slope_eps_low, slope_eps_high;
    double b_low, b_high, b_eps_low, b_eps_high;
    double e_trans_low, e_trans_high, p_trans, h_trans;

    slope_low = (pdat[index_pt] - pdat[index_pt-1])/(hdat[index_pt] - hdat[index_pt-1]);
    b_low = (pdat[index_pt] + pdat[index_pt-1] - slope_low * (hdat[index_pt-1] + hdat[index_pt]) )/2. ;

    slope_high = (pdat[index_pt+2] - pdat[index_pt+1])/(hdat[index_pt+2] - hdat[index_pt+1]);
    b_high = (pdat[index_pt+2] + pdat[index_pt+1] - slope_high * (hdat[index_pt+2] + hdat[index_pt+1]) )/2. ;

    h_trans = (b_high - b_low)/(slope_low - slope_high);
    p_trans = slope_low * h_trans + b_low;

    slope_eps_low = (pdat[index_pt] - pdat[index_pt-1])/(edat[index_pt] - edat[index_pt-1]);
    b_eps_low = (pdat[index_pt] + pdat[index_pt-1] - slope_eps_low * (edat[index_pt-1] + edat[index_pt]) )/2. ;

    slope_eps_high = (pdat[index_pt+2] - pdat[index_pt+1])/(edat[index_pt+2] - edat[index_pt+1]);
    b_eps_high = (pdat[index_pt+2] + pdat[index_pt+1] - slope_eps_high * (edat[index_pt+2] + edat[index_pt+1]) )/2. ;

    e_trans_low = (p_trans - b_eps_low)/slope_eps_low;
    e_trans_high = (p_trans - b_eps_high)/slope_eps_high;

    edat[index_pt] = e_trans_low;
    edat[index_pt+1] = e_trans_high;
    pdat[index_pt] = p_trans;
    pdat[index_pt+1] = p_trans;
    hdat[index_pt] = h_trans;
    hdat[index_pt+1] = h_trans;

    return;
}


/* This function finds phase transitions in a tabulated equation of state
 * and returns a list of indices at which the phase transitions occur.
 */
static int * eos_find_phase_transition(size_t ndat, double *edat, double *pdat, int dirty)
{
    int number_phase_transition = 0;
    int *id_phase_transition = NULL;
    double gradient = 0.0;
    double old_gradient = 0.0;
    double delta_gradient = 0.0;
    double pt_tolerance = 2.;
    double eps_min_pt = 1.5e14 * 1e3 * LAL_G_C2_SI;

    if (ndat < 4) {
        id_phase_transition = LALMalloc(2 * sizeof(int));
        id_phase_transition[0] = 0;
        id_phase_transition[1] = ndat - 1;
        return id_phase_transition;
    }

    bool *pt_occurence = LALMalloc(ndat * sizeof(*pt_occurence));

    for (size_t i = 1; i < ndat; i++){
        pt_occurence[i] = false;
        gradient = (pdat[i] - pdat[i-1])/(edat[i] - edat[i-1]);
        delta_gradient = (gradient - old_gradient)/gradient;

        if (edat[i] > eps_min_pt){ // minimum energy density to find PT (avoids crust instability false positives)
            if (gradient == 0.0){ // clean phase transition
                number_phase_transition += 1;
                pt_occurence[i] = true;
            } else if (fabs(delta_gradient) >= pt_tolerance && dirty == 1) { // dirty phase transition
                if (i >= 3) { // prevent invalid memory access
                    double delP_im2 = pdat[i-1] - pdat[i-2];
                    double delP_im3 = pdat[i-2] - pdat[i-3];
                    if (delP_im2 /delP_im3 < 3.0){ // avoids idenfying pressure jumps around the PT that could false flag
                        number_phase_transition += 1;
                        pt_occurence[i] = true;
                    }
                }
            }
        }

        old_gradient = gradient;
    }

    id_phase_transition = LALMalloc((number_phase_transition + 2) * sizeof(int));
    id_phase_transition[0] = number_phase_transition;
    id_phase_transition[number_phase_transition+1] = ndat-1;

    int count_id = 1;
    for (size_t i = 1; i < ndat; i++){
        if (pt_occurence[i] == true) {
            id_phase_transition[count_id] = i-1;
            count_id += 1;
        }
    }

    LALFree(pt_occurence);
    return id_phase_transition;
}



/*
 * This function creates a LALSimNeutronStarEOS pointer from a piece of
 * tabulated equation of state tabulated data, given minimum and maximum
 * indices for the tables.
 */
static LALSimNeutronStarEOS * eos_piece_alloc_tabular( double *nbdat, double *edat, double *pdat,
                                                double *mubdat, double *muedat, double *hdat,
                                                double *yedat, double *cs2dat,
                                                size_t begin_index, size_t end_index){

    LALSimNeutronStarEOS * eos;

    double *nbdat_cut, *edat_cut, *pdat_cut, *mubdat_cut, *muedat_cut, *hdat_cut, *yedat_cut, *cs2dat_cut;
    size_t ndat;

    ndat = end_index - begin_index + 1;

    nbdat_cut  = LALMalloc(ndat * sizeof(*nbdat));
    edat_cut   = LALMalloc(ndat * sizeof(*edat));
    pdat_cut   = LALMalloc(ndat * sizeof(*pdat));
    mubdat_cut = LALMalloc(ndat * sizeof(*mubdat));
    muedat_cut = LALMalloc(ndat * sizeof(*muedat));
    hdat_cut   = LALMalloc(ndat * sizeof(*hdat));
    yedat_cut  = LALMalloc(ndat * sizeof(*yedat));
    cs2dat_cut = LALMalloc(ndat * sizeof(*cs2dat));

    // Append the EoS before the phase transition
    int ncol = 9;
    if (hdat == NULL) {
        ncol = 2;
        for (size_t i = 0 ; i < ndat ; i++){
            nbdat_cut[i]   = 0.0;
            edat_cut[i]    = edat[begin_index+i];
            pdat_cut[i]    = pdat[begin_index+i];
            mubdat_cut[i]  = 0.0;
            muedat_cut[i]  = 0.0;
            hdat_cut[i]    = 0.0;
            yedat_cut[i]   = 0.0;
            cs2dat_cut[i]  = 0.0;
        }
    } else {
        for (size_t i = 0 ; i < ndat ; i++){
            nbdat_cut[i]   = nbdat[begin_index+i];
            edat_cut[i]    = edat[begin_index+i];
            pdat_cut[i]    = pdat[begin_index+i];
            mubdat_cut[i]  = mubdat[begin_index+i];
            muedat_cut[i]  = muedat[begin_index+i];
            hdat_cut[i]    = hdat[begin_index+i];
            yedat_cut[i]   = yedat[begin_index+i];
            cs2dat_cut[i]  = cs2dat[begin_index+i];
        }
    }

    eos =     eos_alloc_tabular(nbdat_cut, edat_cut, pdat_cut, mubdat_cut, muedat_cut, hdat_cut, yedat_cut, cs2dat_cut, ndat, ncol);

    LALFree(nbdat_cut);
    LALFree(edat_cut);
    LALFree(pdat_cut);
    LALFree(mubdat_cut);
    LALFree(muedat_cut);
    LALFree(hdat_cut);
    LALFree(yedat_cut);
    LALFree(cs2dat_cut);

    return eos;
}

/** @endcond */

/**
 * @brief Reads a data file containing tabulated neutron star
 * equation of state data to create the LALSimNeutronStarEOS
 * equation of state structure.
 * @details Read a data file specified by a path fname that contains either
 * i) 2 whitespace separated columns of equation of state data ("old" LAL EoS format)
 * with the pressure in /m^2 (first column) and the energy density in /m^2 (second column).
 * ii) 9 whitespace separated columns of equation of state data ("new" LAL EoS format)
 * with the table index, the baryon density in /fm^3, the energy density in g/cm^3,
 * the pressure in dyn/cm^2, the baryon chemical potential in MeV, the electron
 * chemical potential in MeV, the pseudo-enthalpy, the lepton fraction and the
 * square of the speed of sound normalized to light velocity.
 *
 * Every line beginning with the character '#' is ignored.
 * If the path is an absolute path then this specific file is opened;
 * otherwise, search for the file in paths given in the environment variable
 * LALSIM_DATA_PATH, and finally search in the installed PKG_DATA_DIR path.
 * @param[in] fname The path of the file to open.
 * @return A pointer to neutron star equation of state structure (LALSimNeutronStarEOS).
 */
LALSimNeutronStarEOS *XLALSimNeutronStarEOSFromFile(const char *fname){
    LALSimNeutronStarEOS *eos;
    double *f_dat;
    size_t ncol;
    size_t ndat;
    LALFILE *fp;

    double *nbdat;
    double *edat;
    double *pdat;
    double *mubdat;
    double *muedat;
    double *hdat;
    double *yedat;
    double *cs2dat;

    fp = XLALSimReadDataFileOpen(fname);
    if (!fp)
        XLAL_ERROR_NULL(XLAL_EFUNC);

    ndat = XLALSimReadDataFileNCol(&f_dat, &ncol, fp);
    XLALFileClose(fp);

    if (ndat == (size_t) (-1))
        XLAL_ERROR_NULL(XLAL_EFUNC);

    nbdat = LALMalloc(ndat * sizeof(*nbdat));
    edat = LALMalloc(ndat * sizeof(*edat));
    pdat = LALMalloc(ndat * sizeof(*pdat));
    mubdat = LALMalloc(ndat * sizeof(*mubdat));
    muedat = LALMalloc(ndat * sizeof(*muedat));
    hdat = LALMalloc(ndat * sizeof(*hdat));
    yedat = LALMalloc(ndat * sizeof(*yedat));
    cs2dat = LALMalloc(ndat * sizeof(*cs2dat));


    if (ncol > 2)
    {
        for (size_t i = 0 ; i < ndat ; i++) {
            nbdat[i] = f_dat[i * ncol + 1];
            edat[i] = f_dat[i * ncol + 2] * 1e3 * LAL_G_C2_SI; /* transform from CGS to SI and then to Geometrized units */
            pdat[i] = f_dat[i * ncol + 3] * 1e-1 * LAL_G_C4_SI; /* transform from CGS to SI and then to Geometrized units */
            mubdat[i] = f_dat[i * ncol + 4];
            muedat[i] = f_dat[i * ncol + 5];
            hdat[i] = f_dat[i * ncol + 6];
            yedat[i] = f_dat[i * ncol + 7];
            cs2dat[i] = f_dat[i * ncol + 8];
        }
    }
    else if (ncol == 2)
    {
        for (size_t i = 0 ; i < ndat ; i++) {
            pdat[i] = f_dat[i * ncol];
            edat[i] = f_dat[i * ncol + 1];
        }

        LALFree(nbdat);  nbdat  = NULL;
        LALFree(mubdat); mubdat = NULL;
        LALFree(muedat); muedat = NULL;
        LALFree(hdat);   hdat   = NULL;
        LALFree(yedat);  yedat  = NULL;
        LALFree(cs2dat); cs2dat = NULL;
    }
    else if (ncol < 2)
    {
        fprintf(stderr, "error: equation of state files must have at least 2 columns, ncol >= 2\n");
        exit(1);
    }

    eos = eos_alloc_tabular(nbdat, edat, pdat, mubdat, muedat, hdat, yedat, cs2dat, ndat, ncol);

    XLALFree(f_dat);
    LALFree(nbdat);
    LALFree(edat);
    LALFree(pdat);
    LALFree(mubdat);
    LALFree(muedat);
    LALFree(hdat);
    LALFree(yedat);
    LALFree(cs2dat);

    snprintf(eos->name, sizeof(eos->name), "%s", fname);
    return eos;
}

/**
 * @brief Reads 9 arrays of neutron star equation of state data
 * to create the LALSimNeutronStarEOS equation of state structure.
 * @details The arrays read contain each ndat lines of equation of state
 * data. Although the 9 arrays correspond to the physical quantities provided
 * in the "new" LAL format equation of state files, the units are not the same.
 *
 * The energy density (edat) and the pressure (edat) cannot be NULL or zero-filled arrays.
 * It is highly recommanded that the user provides the pseudo-enthalpy as input,
 * particularly if the equation of state data contains phase transitions.
 *
 * If the pseudo-enthalpy (hdat):
 * - is NULL, it will be calculated within the function from the pressure and energy
 * density; in this case, the baryonic density (nbdat), the baryonic and electronic
 * chemical potentials (mubdat and muedat), the lepton fraction (yedat) and sound
 * speed squared (cs2dat) can also be NULL as they are calculated within the function
 * from the pressure and energy density. This is similar to providing an equation of
 * state in the "old" LAL file format.
 * - is not NULL, it is necessary to provide non-NULL arrays of the baryonic density
 * (nbdat), the baryonic and electronic chemical potentials (mubdat and muedat), the
 * lepton fraction (yedat) and sound speed squared (cs2dat). Note that mubdat, muedat,
 * yedat can be zero-filled array without consequence, as those quantities are not
 * used in LALSimulation as of now; this is the same for nbdat, as the rest-mass
 * density is recalculated from hdat, pdat and edat and used in the related
 * LALSimulation functions. A zero-filled array for cs2dat will not lead to an error
 * when constructing the EOS structure and the input of this quantity is not necessary
 * to solve neutron star's astrophysical parameters; however, any function related to
 * the speed of sound in LALSimulation will result in errors. This is similar to
 * providing an equation of state in the "new" LAL file format, although units are different.
 * @param nbdat Array for the baryon density (in /fm^3).
 * @param edat Array for the energy density (in m^-2).
 * @param pdat Array for the pressure in (in m^-2).
 * @param mubdat Array for the baryon chemical potential (in MeV).
 * @param muedat Array for the electron chemical potential (in MeV).
 * @param hdat Array for the pseudo enthalpy (dimensionless).
 * @param yedat Array for the lepton fration (dimensionless).
 * @param cs2dat Array for the sound speed squared normalized
 * to the speed of light (dimensionless).
 * @param ndat Size of the arrays for equation of state quantities.
 * @return A pointer to neutron star equation of state structure (LALSimNeutronStarEOS).
 */
LALSimNeutronStarEOS *XLALSimNeutronStarEOSFromTabData(double *nbdat, double *edat, double *pdat,
   double *mubdat, double *muedat, double *hdat, double *yedat, double *cs2dat, size_t ndat)
{
    LALSimNeutronStarEOS *eos;
    int ncol = 9;
    if (hdat == NULL) ncol = 2;
    eos = eos_alloc_tabular(nbdat, edat, pdat, mubdat, muedat, hdat, yedat, cs2dat, ndat, ncol);
    return eos;
}

/**
 * @brief Creates a tabulated neutron star equation of state from
 * energy density and pressure arrays.
 * @details This is a convenience wrapper around XLALSimNeutronStarEOSFromTabData
 * that accepts REAL8Vector inputs, making it usable from Python via SWIG.
 *
 * @param energy_density Array for the energy density (in m^-2).
 * @param pressure Array for the pressure (in m^-2).
 * @return A pointer to neutron star equation of state structure (LALSimNeutronStarEOS).
 */
LALSimNeutronStarEOS *XLALSimNeutronStarEOSFromArrays(
    const REAL8Vector *energy_density, const REAL8Vector *pressure)
{
    XLAL_CHECK_NULL(energy_density && pressure, XLAL_EFAULT);
    XLAL_CHECK_NULL(energy_density->length == pressure->length, XLAL_ESIZE,
        "energy_density and pressure must have the same length");
    return XLALSimNeutronStarEOSFromTabData(
        NULL, energy_density->data, pressure->data,
        NULL, NULL, NULL, NULL, NULL, pressure->length);
}

/**
 * @brief Creates a phase-transition-aware tabulated neutron star equation of
 * state from energy density and pressure arrays.
 * @details This is a convenience wrapper around
 * XLALSimNeutronStarEOSFromTabDataPhaseTransition that accepts REAL8Vector
 * inputs, making it usable from Python via SWIG. Use this instead of
 * XLALSimNeutronStarEOSFromArrays when the EOS contains a first-order phase
 * transition (i.e., pressure is not strictly increasing).
 *
 * @param energy_density Array for the energy density (in m^-2).
 * @param pressure Array for the pressure (in m^-2).
 * @return A pointer to EOSMultiParts equation of state structure.
 */
EOSMultiParts *XLALSimNeutronStarEOSFromArraysPhaseTransition(
    const REAL8Vector *energy_density, const REAL8Vector *pressure)
{
    XLAL_CHECK_NULL(energy_density && pressure, XLAL_EFAULT);
    XLAL_CHECK_NULL(energy_density->length == pressure->length, XLAL_ESIZE,
        "energy_density and pressure must have the same length");
    return XLALSimNeutronStarEOSFromTabDataPhaseTransition(
        NULL, energy_density->data, pressure->data,
        NULL, NULL, NULL, NULL, NULL, pressure->length);
}


/**
 * @brief Reads a data file containing tabulated equation of state data
 * to create the EOSMultiParts equation of state structure that can handle
 * equations of state with phase transitions.
 * @details Read a data file specified by a path fname that contains either
 * i) 2 whitespace separated columns of equation of state data ("old" LAL EoS format)
 * with the pressure in /m^2 (first column) and the energy density in /m^2 (second column).
 * ii) 9 whitespace separated columns of equation of state data ("new" LAL EoS format)
 * with the table index, the baryon density in /fm^3, the energy density in g/cm^3,
 * the pressure in dyn/cm^2, the baryon chemical potential in MeV, the electron
 * chemical potential in MeV, the pseudo-enthalpy, the lepton fraction and the
 * square of the speed of sound normalized to light velocity.
 *
 * Every line beginning with the character '#' is ignored.
 * If the path is an absolute path then this specific file is opened;
 * otherwise, search for the file in paths given in the environment variable
 * LALSIM_DATA_PATH, and finally search in the installed PKG_DATA_DIR path.
 *
 * This function builds the EOSMultiParts structure from equation of state
 * data that can include a first order phase transition. The equation of state data
 * is tested for phase transitions which are numerically defined by a pressure
 * plateau or near like plateau associated to a jump in energy density.
 * If N phase transitions are found, EOSMultiParts contains N+1 LALSimNeutronStarEOS
 * equation of state pieces.
 *
 * We define a "dirty" phase transition by an intended first order phase
 * transition that is not numerically so: at the upper (+) and lower (-)
 * boundaries of the phase transition (pt), the user input equation of
 * state data has pressure P(pt,+) - P(pt,-) > 0.
 * If a "dirty" phase transition is detected, it is corrected as follows:
 * - For the "new" LAL format: a linear extrapolation is used to correct
 * the two points in the equation of state data grid involved with the
 * phase transition and define a clean phase transition. The pseudo-enthalpy h
 * and pressure P at the boundary of the phase transition are recomputed such
 * that h_pt(pt,-) = h(pt,+) and P(pt,-) = P(pt,+); the corresponding energy density
 * at the upper and lower boundaries of the newly defined clean phase transition
 * are computed also with a linear extrapolation.
 * - For the "old" LAL format: as the pseudo-enthalpy is calculated, its value
 * at the boundaries of the phase transition is corrected as h(pt,+) -> h(pt,-).
 * The pressure is also corrected such that P(pt,+) -> P(pt,-) and the energy
 * density (eps) points at the boundary of the phase transition are kept intact.
 * The existing points of the equation of state data provided by the user
 * have been modified to obtain a clean phase transition: this implies that
 * P(pt,+) and eps(pt,+) are not thermodynamically coherent together. The user is
 * advised to provide equations of state containing phase transition in the
 * "new" LAL format.
 * In the event that the user does not wish for the phase transition to be
 * corrected, use the LALSimNeutronStarEOS structure and the corresponding
 * XLALSimNeutronStarEOSFromFile; note that it cannot be used in the
 * solver for neutron star's astrophysical parameters that accounts for the
 * necessary phase transition corrections.
 * @param[in] fname The path of the file to open.
 * @return A pointer to neutron star equation of state structure EOSMultiParts.
 */
EOSMultiParts *XLALSimNeutronStarEOSFromFilePhaseTransition(const char *fname) {

    EOSMultiParts *eos;

    double *f_dat;
    size_t ncol;
    size_t ndat;
    LALFILE *fp;

    double *nbdat;
    double *edat;
    double *pdat;
    double *mubdat;
    double *muedat;
    double *hdat;
    double *yedat;
    double *cs2dat;

    fp = XLALSimReadDataFileOpen(fname);
    if (!fp)
        XLAL_ERROR_NULL(XLAL_EFUNC);

    ndat = XLALSimReadDataFileNCol(&f_dat, &ncol, fp);
    XLALFileClose(fp);

    if (ndat == (size_t) (-1))
        XLAL_ERROR_NULL(XLAL_EFUNC);

    nbdat = LALMalloc(ndat * sizeof(*nbdat));
    edat = LALMalloc(ndat * sizeof(*edat));
    pdat = LALMalloc(ndat * sizeof(*pdat));
    mubdat = LALMalloc(ndat * sizeof(*mubdat));
    muedat = LALMalloc(ndat * sizeof(*muedat));
    hdat = LALMalloc(ndat * sizeof(*hdat));
    yedat = LALMalloc(ndat * sizeof(*yedat));
    cs2dat = LALMalloc(ndat * sizeof(*cs2dat));


    if (ncol > 2)
    {
        for (size_t i = 0 ; i < ndat ; i++) {
            nbdat[i] = f_dat[i * ncol + 1];
            edat[i] = f_dat[i * ncol + 2] * 1e3 * LAL_G_C2_SI; /* transform from CGS to SI and then to Geometrized units */
            pdat[i] = f_dat[i * ncol + 3] * 1e-1 * LAL_G_C4_SI; /* transform from CGS to SI and then to Geometrized units */
            mubdat[i] = f_dat[i * ncol + 4];
            muedat[i] = f_dat[i * ncol + 5];
            hdat[i] = f_dat[i * ncol + 6];
            yedat[i] = f_dat[i * ncol + 7];
            cs2dat[i] = f_dat[i * ncol + 8];
        }

    }
    else if (ncol == 2)
    {
        for (size_t i = 0 ; i < ndat ; i++) {
            pdat[i] = f_dat[i * ncol];
            edat[i] = f_dat[i * ncol + 1];
        }

        LALFree(nbdat);  nbdat  = NULL;
        LALFree(mubdat); mubdat = NULL;
        LALFree(muedat); muedat = NULL;
        LALFree(hdat);   hdat   = NULL;
        LALFree(yedat);  yedat  = NULL;
        LALFree(cs2dat); cs2dat = NULL;
    }
    else if (ncol < 2)
    {
        fprintf(stderr, "error: equation of state files must have at least 2 columns, ncol >= 2\n");
        exit(1);
    }

    eos = XLALSimNeutronStarEOSFromTabDataPhaseTransition(nbdat, edat, pdat, mubdat, muedat, hdat,
                                                                    yedat, cs2dat, ndat);
    XLALFree(f_dat);
    LALFree(nbdat);
    LALFree(edat);
    LALFree(pdat);
    LALFree(mubdat);
    LALFree(muedat);
    LALFree(hdat);
    LALFree(yedat);
    LALFree(cs2dat);

    snprintf(eos->name, sizeof(eos->name), "%s", fname);
    return eos;
}

/**
 * @brief Reads 9 arrays of neutron star equation of state data to create
 * the EOSMultiParts equation of state structure that can handle
 * equations of state with phase transitions.
 * @details The arrays read contain each ndat lines of equation of state
 * data. Although the 9 arrays correspond to the physical quantities provided
 * in the "new" LAL format equation of state files, the units are not the same.
 *
 * The energy density (edat) and the pressure (edat) cannot be NULL or zero-filled arrays.
 * It is highly recommanded that the user provides the pseudo-enthalpy as input,
 * particularly if the equation of state data contains phase transitions.
 *
 * If the pseudo-enthalpy (hdat):
 * - is NULL, it will be calculated within the function from the pressure and energy
 * density; in this case, the baryonic density (nbdat), the baryonic and electronic
 * chemical potentials (mubdat and muedat), the lepton fraction (yedat) and sound
 * speed squared (cs2dat) can also be NULL as they are calculated within the function
 * from the pressure and energy density. This is similar to providing an equation of
 * state in the "old" LAL file format.
 * - is not NULL, it is necessary to provide non-NULL arrays of the baryonic density
 * (nbdat), the baryonic and electronic chemical potentials (mubdat and muedat), the
 * lepton fraction (yedat) and sound speed squared (cs2dat). Note that mubdat, muedat,
 * yedat can be zero-filled array without consequence, as those quantities are not
 * used in LALSimulation as of now; this is the same for nbdat, as the rest-mass
 * density is recalculated from hdat, pdat and edat and used in the related
 * LALSimulation functions. A zero-filled array for cs2dat will not lead to an error
 * when constructing the EOS structure and the input of this quantity is not necessary
 * to solve neutron star's astrophysical parameters; however, any function related to
 * the speed of sound in LALSimulation will result in errors. This is similar to
 * providing an equation of state in the "new" LAL file format, although units are different.
 *
 * This function builds the EOSMultiParts structure from equation of state
 * data that can include a first order phase transition. The equation of state data
 * is tested for phase transitions which are numerically defined by a pressure
 * plateau or near like plateau associated to a jump in energy density.
 * If N phase transitions are found, EOSMultiParts contains N+1 LALSimNeutronStarEOS
 * equation of state pieces.
 *
 * We define a "dirty" phase transition by an intended first order phase
 * transition that is not numerically so: at the upper (+) and lower (-)
 * boundaries of the phase transition (pt), the user input equation of
 * state data has P(pt,+)- P(pt,-) > 0.
 * The user can choose to test or not for dirty phase transitions.
 * If a "dirty" phase transition is detected, it is corrected as follows:
 * - If hdat is provided: a linear extrapolation is used to correct
 * the two points in the equation of state data grid involved with the
 * phase transition and define a clean phase transition. The pseudo-enthalpy
 * and pressure at the boundary of the phase transition are recomputed such
 * that h(pt,-) = h(pt,+) and P(pt,-) = P(pt,+); the corresponding energy density
 * at the upper and lower boundaries of the newly defined clean phase transition
 * are computed also with a linear extrapolation.
 * - If hdat is NULL: as the pseudo-enthalpy is calculated, its value
 * at the boundaries of the phase transition is corrected as h(pt,+) -> h(pt,-).
 * The pressure is also corrected such that P(pt,+) -> P(pt,-) and the energy
 * density points at the boundary of the phase transition are kept intact.
 * The existing points of the equation of state data provided by the user
 * have been modified to obtain a clean phase transition: this implies that
 * P(pt,+) and eps(pt,+) are not thermodynamically coherent together.
 * In the event that the user does not wish for the phase transition to be
 * corrected, use the LALSimNeutronStarEOS structure and the corresponding
 * XLALSimNeutronStarEOSFromTabData function; note that it cannot be used in the
 * solver for neutron star's astrophysical parameters that accounts for the
 * necessary phase transition corrections.
 *
 * @param nbdat Array for the baryon density (in /fm^3).
 * @param edat Array for the energy density (in m^-2).
 * @param pdat Array for the pressure in (in m^-2).
 * @param mubdat Array for the baryon chemical potential (in MeV).
 * @param muedat Array for the electron chemical potential (in MeV).
 * @param hdat Array for the pseudo enthalpy (dimensionless).
 * @param yedat Array for the lepton fration (dimensionless).
 * @param cs2dat Array for the sound speed squared normalized
 * to the speed of light (dimensionless).
 * @param ndat Size of the arrays for equation of state quantities.
 * @param dirty Integer to test for dirty phase transitions (1) or clean ones only (0).
 * @return A pointer to neutron star equation of state structure EOSMultiParts.
 */
EOSMultiParts *XLALSimNeutronStarEOSFromTabDataPhaseTransitionChoiceDirtyPT(double *nbdat, double *edat, double *pdat,
                                                                    double *mubdat, double *muedat, double *hdat,
                                                                    double *yedat, double *cs2dat, size_t ndat, int dirty)
{

    EOSMultiParts *eos = NULL;
    int *indices_phase_transition = NULL;

    eos = LALCalloc(1, sizeof(*eos));
    if (!eos) return NULL;
    eos->free = eos_multi_part_free_tabular;

    /* Inquire about phase transitions in the equation of state */
    indices_phase_transition = eos_find_phase_transition(ndat, edat, pdat, dirty);
    if (!indices_phase_transition) goto cleanup;
    int number_pt = indices_phase_transition[0] ;
    int number_eos = number_pt + 1 ;

    eos->number_of_parts = number_eos;
    eos->pmin = pdat[0];
    eos->pmax = pdat[ndat - 1];
    /* Allocate each piece of the equation of state separated by a phase transition */
    eos->eos_piece = XLALCalloc(number_eos, sizeof(LALSimNeutronStarEOS *));
    if (!eos->eos_piece) goto cleanup;


    if (number_pt != 0) {
        for (int i = 1; i <= number_pt; i++){
            if (indices_phase_transition[i+1] - indices_phase_transition[i] < 4) {
                XLALPrintError("EoS Piece contains too few points.");
                goto cleanup;
            }
            if (pdat[indices_phase_transition[i]] == pdat[indices_phase_transition[i]+1]){
                printf("\t Phase transition found at index %d. This phase transition is clean.\n", indices_phase_transition[i]);
            } else {
                printf("\t Phase transition found at index %d. This phase transition is dirty.\n", indices_phase_transition[i]);
                if (hdat != NULL ) {
                    eos_correct_phase_transition(edat, pdat, hdat, indices_phase_transition[i]) ;
                } else {
                    printf("\t The dirty phase transition will be cleaned (with a simple method) after the enthalpy is recalculated.\n");
                }
            }
        }
    }

    /* Construct the multiple part equation of state */
    size_t bottom_index = 0, upper_index = 0;
    for (int i = 0; i <= number_pt; i++){
        upper_index = indices_phase_transition[i + 1];
        if (upper_index <= bottom_index + 2 || upper_index >= ndat) {
            printf("Skipping invalid EOS segment [%zu,%zu]\n",
                   bottom_index, upper_index);
            goto cleanup;
        }

        if (eos->eos_piece[i] != NULL) {
            eos->eos_piece[i]->free(eos->eos_piece[i]);
            eos->eos_piece[i] = NULL;
        }
        eos->eos_piece[i] = eos_piece_alloc_tabular(nbdat, edat, pdat, mubdat, muedat, hdat, yedat, cs2dat, bottom_index, upper_index);
        if (!eos->eos_piece[i]) goto cleanup;
        size_t next_index = indices_phase_transition[i + 1];

        if (next_index <= bottom_index + 2 || next_index >= ndat) {
            printf("Skipping invalid EOS segment: [%zu, %zu]\n",
                bottom_index, next_index);
            continue;   // IMPORTANT: do NOT advance bottom_index
        }

        bottom_index = next_index + 1;
    }

    // If the enthalpy was not provided and the EoS has a Phase Transition,
    // we substract the value h[0] in the piece and also add the enthalpy
    // at the point of piece N-1. The entire EoS object for the piece is then
    // recalculated assuming new LAL format.
    if (hdat == NULL && number_pt > 0){
        double *hdat_recal = XLALMalloc(ndat * sizeof(*hdat_recal));
        double *nbdat_recal  = XLALMalloc(ndat * sizeof(*nbdat_recal));
        double *edat_recal = XLALMalloc(ndat * sizeof(*edat_recal));
        double *pdat_recal = XLALMalloc(ndat * sizeof(*pdat_recal));
        double *mubdat_recal = XLALMalloc(ndat * sizeof(*mubdat_recal));
        double *muedat_recal = XLALMalloc(ndat * sizeof(*muedat_recal));
        double *yedat_recal  = XLALMalloc(ndat * sizeof(*yedat_recal));
        double *cs2dat_recal = XLALMalloc(ndat * sizeof(*cs2dat_recal));
        if (!hdat_recal || !nbdat_recal || !mubdat_recal ||
            !muedat_recal || !yedat_recal || !cs2dat_recal ||
            !edat_recal || !pdat_recal) {
            XLALFree(hdat_recal);
            XLALFree(nbdat_recal);
            XLALFree(mubdat_recal);
            XLALFree(muedat_recal);
            XLALFree(yedat_recal);
            XLALFree(cs2dat_recal);
            XLALFree(edat_recal);
            XLALFree(pdat_recal);
            goto cleanup;
        }

        size_t ndat_total = 0;
        for (int i = 0; i <= number_pt; i++){
            size_t ndat_piece = eos->eos_piece[i]->data.tabular->ndat;
            if (i == 0){
                for (size_t j = 0 ; j < ndat_piece; j++){
                    hdat_recal[ndat_total] = exp(eos->eos_piece[i]->data.tabular->log_hdat[j]);
                    edat_recal[ndat_total] = exp(eos->eos_piece[i]->data.tabular->log_edat[j]);
                    pdat_recal[ndat_total] = exp(eos->eos_piece[i]->data.tabular->log_pdat[j]);
                    ndat_total += 1;
                }
            } else {
                for (size_t j = 0 ; j < ndat_piece; j++){
                    if (j == 0) {
                        hdat_recal[ndat_total] = hdat_recal[ndat_total-j-1];
                        pdat_recal[ndat_total] = exp(eos->eos_piece[i-1]->data.tabular->log_pdat[eos->eos_piece[i-1]->data.tabular->ndat - 1]);
                    } else {
                        hdat_recal[ndat_total] = exp(eos->eos_piece[i]->data.tabular->log_hdat[j]) - exp(eos->eos_piece[i]->data.tabular->log_hdat[0]) + hdat_recal[ndat_total-j-1];
                        pdat_recal[ndat_total] = exp(eos->eos_piece[i]->data.tabular->log_pdat[j]);
                    }
                    edat_recal[ndat_total] = exp(eos->eos_piece[i]->data.tabular->log_edat[j]);
                    ndat_total += 1;
                }
            }
        }
        // Refill the EOS structure
        bottom_index = 0, upper_index = 0;
        for (int i = 0; i <= number_pt; i++){
            upper_index = indices_phase_transition[i+1];
            if (eos->eos_piece[i] != NULL) {
                eos->eos_piece[i]->free(eos->eos_piece[i]);
                eos->eos_piece[i] = NULL;
            }
            eos->eos_piece[i] = eos_piece_alloc_tabular(
                nbdat, edat, pdat,
                mubdat, muedat, hdat,
                yedat, cs2dat,
                bottom_index, upper_index
            );
            if (!eos->eos_piece[i]) goto cleanup;
            size_t next_index = indices_phase_transition[i + 1];
            if (next_index <= bottom_index + 2 || next_index >= ndat) {
                printf("Skipping invalid EOS segment: [%zu, %zu]\n",
                    bottom_index, next_index);
                continue;   // IMPORTANT: do NOT advance bottom_index
            }

            bottom_index = next_index + 1;
        }
        XLALFree(hdat_recal);
        XLALFree(nbdat_recal);
        XLALFree(mubdat_recal);
        XLALFree(muedat_recal);
        XLALFree(yedat_recal);
        XLALFree(cs2dat_recal);
        XLALFree(edat_recal);
        XLALFree(pdat_recal);
    }

    eos->hmin = XLALSimNeutronStarEOSMinPseudoEnthalpy(eos->eos_piece[0]);
    eos->hmax = XLALSimNeutronStarEOSMaxPseudoEnthalpy(eos->eos_piece[eos->number_of_parts-1]);
    eos->hMinAcausal = eosMultiParts_min_acausal_pseudo_enthalpy_tabular(eos->hmax, eos);

    char name[LALNameLength] = "unknown_eos_name";
    snprintf(eos->name, sizeof(eos->name), "%s", name);
    XLALFree(indices_phase_transition);
    return eos;

cleanup:

    if (eos) {
        if (eos->eos_piece) {
            for (int i = 0; i < eos->number_of_parts; i++) {
                if (eos->eos_piece[i]) {

                    // CRITICAL: use proper destructor
                    if (eos->eos_piece[i]->free) {
                        eos->eos_piece[i]->free(eos->eos_piece[i]);
                    } else {
                        LALFree(eos->eos_piece[i]);
                    }

                    eos->eos_piece[i] = NULL;
                }
            }
            LALFree(eos->eos_piece);
        }

        LALFree(eos);
    }

    if (indices_phase_transition)
        XLALFree(indices_phase_transition);

    return NULL;
}






/**
 * @brief Reads 9 arrays of neutron star equation of state data to create
 * the EOSMultiParts equation of state structure that can handle
 * equations of state with phase transitions (including "dirty" ones).
 * @details The arrays read contain each ndat lines of equation of state
 * data. Although the 9 arrays correspond to the physical quantities provided
 * in the "new" LAL format equation of state files, the units are not the same.
 *
 * The energy density (edat) and the pressure (edat) cannot be NULL or zero-filled arrays.
 * It is highly recommanded that the user provides the pseudo-enthalpy as input,
 * particularly if the equation of state data contains phase transitions.
 *
 * If the pseudo-enthalpy (hdat):
 * - is NULL, it will be calculated within the function from the pressure and energy
 * density; in this case, the baryonic density (nbdat), the baryonic and electronic
 * chemical potentials (mubdat and muedat), the lepton fraction (yedat) and sound
 * speed squared (cs2dat) can also be NULL as they are calculated within the function
 * from the pressure and energy density. This is similar to providing an equation of
 * state in the "old" LAL file format.
 * - is not NULL, it is necessary to provide non-NULL arrays of the baryonic density
 * (nbdat), the baryonic and electronic chemical potentials (mubdat and muedat), the
 * lepton fraction (yedat) and sound speed squared (cs2dat). Note that mubdat, muedat,
 * yedat can be zero-filled array without consequence, as those quantities are not
 * used in LALSimulation as of now; this is the same for nbdat, as the rest-mass
 * density is recalculated from hdat, pdat and edat and used in the related
 * LALSimulation functions. A zero-filled array for cs2dat will not lead to an error
 * when constructing the EOS structure and the input of this quantity is not necessary
 * to solve neutron star's astrophysical parameters; however, any function related to
 * the speed of sound in LALSimulation will result in errors. This is similar to
 * providing an equation of state in the "new" LAL file format, although units are different.
 *
 * This function builds the EOSMultiParts structure from equation of state
 * data that can include a first order phase transition. The equation of state data
 * is tested for phase transitions which are numerically defined by a pressure
 * plateau or near like plateau associated to a jump in energy density.
 * If N phase transitions are found, EOSMultiParts contains N+1 LALSimNeutronStarEOS
 * equation of state pieces.
 *
 * We define a "dirty" phase transition by an intended first order phase
 * transition that is not numerically so: at the upper (+) and lower (-)
 * boundaries of the phase transition (pt), the user input equation of
 * state data has P(pt,+)- P(pt,-) > 0.
 * If a "dirty" phase transition is detected, it is corrected as follows:
 * - If hdat is provided: a linear extrapolation is used to correct
 * the two points in the equation of state data grid involved with the
 * phase transition and define a clean phase transition. The pseudo-enthalpy
 * and pressure at the boundary of the phase transition are recomputed such
 * that h(pt,-) = h(pt,+) and P(pt,-) = P(pt,+); the corresponding energy density
 * at the upper and lower boundaries of the newly defined clean phase transition
 * are computed also with a linear extrapolation.
 * - If hdat is NULL: as the pseudo-enthalpy is calculated, its value
 * at the boundaries of the phase transition is corrected as h(pt,+) -> h(pt,-).
 * The pressure is also corrected such that P(pt,+) -> P(pt,-) and the energy
 * density points at the boundary of the phase transition are kept intact.
 * The existing points of the equation of state data provided by the user
 * have been modified to obtain a clean phase transition: this implies that
 * P(pt,+) and eps(pt,+) are not thermodynamically coherent together.
 * In the event that the user does not wish for the phase transition to be
 * corrected, use the LALSimNeutronStarEOS structure and the corresponding
 * XLALSimNeutronStarEOSFromTabData function; note that it cannot be used in the
 * solver for neutron star's astrophysical parameters that accounts for the
 * necessary phase transition corrections.
 *
 * @param nbdat Array for the baryon density (in /fm^3).
 * @param edat Array for the energy density (in m^-2).
 * @param pdat Array for the pressure in (in m^-2).
 * @param mubdat Array for the baryon chemical potential (in MeV).
 * @param muedat Array for the electron chemical potential (in MeV).
 * @param hdat Array for the pseudo enthalpy (dimensionless).
 * @param yedat Array for the lepton fration (dimensionless).
 * @param cs2dat Array for the sound speed squared normalized
 * to the speed of light (dimensionless).
 * @param ndat Size of the arrays for equation of state quantities.
 * @return A pointer to neutron star equation of state structure EOSMultiParts.
 */
EOSMultiParts *XLALSimNeutronStarEOSFromTabDataPhaseTransition( double *nbdat, double *edat, double *pdat,
                                                                    double *mubdat, double *muedat, double *hdat,
                                                                    double *yedat, double *cs2dat, size_t ndat)
{

   return XLALSimNeutronStarEOSFromTabDataPhaseTransitionChoiceDirtyPT(nbdat, edat, pdat, mubdat, muedat, hdat, yedat, cs2dat, ndat, 1);
}



/**
 * @brief Creates an equation of state structure from tabulated equation
 * of state data of a known name. The name of the tabulated equation of state
 * must belong to the sample of equations of state from the old frame work or
 *  added for the new framework.
 * @details A known, installed, named tabulated equation of state data file, whose name
 * is included in the old EOS framework names or the new ones, is read and then used to
 * create the equation of state structure.
 * The equations of state for the OLD framework available are the representative sample drawn from
 * http://xtreme.as.arizona.edu/NeutronStars/ they are:
 * - ALF1
 * - ALF2
 * - ALF3
 * - ALF4
 * - AP1
 * - AP2
 * - AP3
 * - AP4
 * - APR4_EPP
 * - BBB2
 * - BGN1H1
 * - BPAL12
 * - BSK19
 * - BSK20
 * - BSK21
 * - ENG
 * - FPS
 * - GNH3
 * - GS1
 * - GS2
 * - H1
 * - H2
 * - H3
 * - H4
 * - H5
 * - H6
 * - H7
 * - MPA1
 * - MS1B
 * - MS1B_PP
 * - MS1_PP
 * - MS1
 * - MS2
 * - PAL6
 * - PCL2
 * - PS
 * - QMC700
 * - SLY4
 * - SLY
 * - SQM1
 * - SQM2
 * - SQM3
 * - WFF1
 * - WFF2
 * - WFF3
 * We also include more modern equations from the CompOSE website
 * https://compose.obspm.fr/ downloaded on 18 June 2018. These EOSs are:
 * - APR
 * - BHF_BBB2
 * - KDE0V
 * - KDE0V1
 * - RS
 * - SK255
 * - SK272
 * - SKA
 * - SKB
 * - SKI2
 * - SKI3
 * - SKI4
 * - SKI5
 * - SKI6
 * - SKMP
 * - SKOP
 * - SLY2
 * - SLY230A
 * - SLY9
 * And we include HQC18 from http://user.numazu-ct.ac.jp/~sumi/eos/HQC18_submit
 * - HQC18
 *
 *
 * @param[in] name The name of the equation of state.
 * @return A pointer to neutron star equation of state structure.
 */
LALSimNeutronStarEOS *XLALSimNeutronStarEOSByName(const char *name)
{
    static const char fname_base[] = "LALSimNeutronStarEOS_";
    static const char fname_extn[] = ".dat";

    size_t n = XLAL_NUM_ELEM(lalSimNeutronStarEOSNames);
    size_t i;
    char fname[FILENAME_MAX];

    for (i = 0; i < n ; ++i)
        if (XLALStringCaseCompare(name, lalSimNeutronStarEOSNames[i]) == 0) {
            LALSimNeutronStarEOS *eos;
            snprintf(fname, sizeof(fname), "%s%s%s", fname_base, lalSimNeutronStarEOSNames[i],
                fname_extn);
            eos = XLALSimNeutronStarEOSFromFile(fname);
            if (!eos)
                XLAL_ERROR_NULL(XLAL_EFUNC);
            snprintf(eos->name, sizeof(eos->name), "%s", lalSimNeutronStarEOSNames[i]);
            return eos;
        }

    XLAL_PRINT_ERROR("Unrecognized EOS name %s...", name);
    XLALPrintError("\tRecognised EOS names are: %s", lalSimNeutronStarEOSNames[0]);
    for (i = 1; i < n ; ++i)
        XLALPrintError(", %s", lalSimNeutronStarEOSNames[i]);
    XLALPrintError("\n");
    XLAL_ERROR_NULL(XLAL_ENAME);
}


/**
 * @brief Creates an equation of state structure from tabulated equation
 * of state data of a known name available in LALSimulation database.
 * @details The name of the tabulated equation of state must belong
 * to the sample of equations of state from the old or the new framework.
 * Available equation of state name included are listed in the following.
 * Regarding the old framework,
 * 1) From http://xtreme.as.arizona.edu/NeutronStars/):
 * - ALF1, ALF2, ALF3, ALF4, AP1, AP2, AP3, AP4, APR4_EPP,
 * BBB2, BGN1H1, BPAL12, BSK19, BSK20, BSK21, ENG, FPS,
 * GNH3, GS1, GS2, H1, H2, H3, H4, H5, H6, H7,
 * MPA1, MS1B, MS1B_PP, MS1_PP, MS1, MS2, PAL6,
 * PCL2, PS, QMC700, SLY4, SLY, SQM1, SQM2, SQM3,
 * WFF1, WFF2, WFF3.
 * 2) From the CompOSE website (https://compose.obspm.fr/, downloaded on 18 June 2018):
 * - APR, BHF_BBB2, KDE0V, KDE0V1, RS, SK255, SK272, SKA, SKB,
 * SKI2, SKI3, SKI4, SKI5, SKI6, SKMP, SKOP, SLY2, SLY230A, SLY9.
 * 3) From http://user.numazu-ct.ac.jp/~sumi/eos/HQC18_submit
 * - HQC18
 * Regarding the new framework (9 column EoS data), based on the CompOSE database
 * (https://compose.obspm.fr/)
 * 1) Files ending with _BSK24, the outer crust is calculated using analytical fits
 * of the BSk24 energy-density functionals (Pearson et al., MNRAS, 481, 2994 (2018)).
 * - GMSR_BSK14_BSK24, GMSR_DHSL59_BSK24, GMSR_DHSL69_BSK24, GMSR_F0_BSK24,
 * GMSR_H1_BSK24, GMSR_H2_BSK24, GMSR_H3_BSK24, GMSR_H4_BSK24, GMSR_H5_BSK24,
 * GMSR_LN55_BSK24, GMSR_SLY5_BSK24, GPPVA_DD2_BSK24, GPPVA_DDME2_BSK24,
 * GPPVA_FSU2_BSK24, GPPVA_FSU2H_BSK24, GPPVA_NL3WRL55_BSK24, KDE0V_BSK24,
 * KDE0V1_BSK24, PCP_BSK24_BSK24, RG_SLY4_BSK24, RS_BSK24, SK255_BSK24,
 * SKA_BSK24, SKB_BSK24, SKI2_BSK24, SKI3_BSK24, SKI4_BSK24, SKI6_BSK24,
 * SKOP_BSK24, SLY2_BSK24, SLY9_BSK24, SLY230A_BSK24, XMLSLZ_DDLZ1_BSK24,
 * XMLSLZ_DDME2_BSK24, XMLSLZ_DDMEX_BSK24, XMLSLZ_GM1_BSK24, XMLSLZ_MTVTC_BSK24,
 * XMLSLZ_NL3_BSK24, XMLSLZ_PKDD_BSK24, XMLSLZ_TM1_BSK24, XMLSLZ_TW99_BSK24
 * 2) Files ending with _META, the outer crust and the inner crust is calculated
 * with the Crust Unified Tool for Equation-of-state Reconstruction (CUTER,
 * Davis et al., Eur. Phys. J. A 61, 120 (2025)). TODO add info on meta data
 * - ABHT_QMC_RMF1_META, ABHT_QMC_RMF2_META, ABHT_QMC_RMF3_META, ABHT_QMC_RMF4_META,
 * BL_CHIRAL_META.
 * @param[in] name The name of the equation of state.
 * @return A pointer to neutron star equation of state structure.
 */
EOSMultiParts *XLALSimNeutronStarEOSMultiPartsByName(const char *name)
{
    static const char fname_base[] = "LALSimNeutronStarEOS_";
    static const char fname_extn[] = ".dat";

    size_t n = XLAL_NUM_ELEM(lalSimNeutronStarEOSNames);
    size_t i;
    char fname[FILENAME_MAX];

    for (i = 0; i < n ; ++i)
        if (XLALStringCaseCompare(name, lalSimNeutronStarEOSNames[i]) == 0) {
            EOSMultiParts *eos;
            snprintf(fname, sizeof(fname), "%s%s%s", fname_base, lalSimNeutronStarEOSNames[i],
                fname_extn);
            eos = XLALSimNeutronStarEOSFromFilePhaseTransition(fname);
            if (!eos)
                XLAL_ERROR_NULL(XLAL_EFUNC);
            snprintf(eos->name, sizeof(eos->name), "%s", lalSimNeutronStarEOSNames[i]);
            return eos;
        }

    XLAL_PRINT_ERROR("Unrecognized EOS name %s...", name);
    XLALPrintError("\tRecognised EOS names are: %s", lalSimNeutronStarEOSNames[0]);
    for (i = 1; i < n ; ++i)
        XLALPrintError(", %s", lalSimNeutronStarEOSNames[i]);
    XLALPrintError("\n");
    XLAL_ERROR_NULL(XLAL_ENAME);
}

/** @} */
/** @} */
