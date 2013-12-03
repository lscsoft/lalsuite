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
 * @addtogroup LALSimNeutronStarEOS_c
 * @{
 */
/**
 * @name Creation routines for tabulated equations of state
 * @{
 */

#include <lal/LALSimReadData.h>
#include <gsl/gsl_interp.h>

/** @cond */

/* Contents of the tabular equation of state data structure. */
struct tagLALSimNeutronStarEOSDataTabular {
    double *pdat;
    double *edat;
    double *hdat;
    double *rhodat;
    size_t ndat;
    gsl_interp *e_of_p_interp;
    gsl_interp *h_of_p_interp;
    gsl_interp *e_of_h_interp;
    gsl_interp *p_of_h_interp;
    gsl_interp *rho_of_h_interp;
    gsl_interp_accel *e_of_p_acc;
    gsl_interp_accel *h_of_p_acc;
    gsl_interp_accel *e_of_h_acc;
    gsl_interp_accel *p_of_h_acc;
    gsl_interp_accel *rho_of_h_acc;
};

static double eos_e_of_p_tabular(double p, LALSimNeutronStarEOS * eos)
{
    double e;
    e = gsl_interp_eval(eos->data.tabular->e_of_p_interp,
        eos->data.tabular->pdat, eos->data.tabular->edat, p,
        eos->data.tabular->e_of_p_acc);
    return e;
}

static double eos_e_of_h_tabular(double h, LALSimNeutronStarEOS * eos)
{
    double e;
    e = gsl_interp_eval(eos->data.tabular->e_of_h_interp,
        eos->data.tabular->hdat, eos->data.tabular->edat, h,
        eos->data.tabular->e_of_h_acc);
    return e;
}

static double eos_p_of_h_tabular(double h, LALSimNeutronStarEOS * eos)
{
    double p;
    p = gsl_interp_eval(eos->data.tabular->p_of_h_interp,
        eos->data.tabular->hdat, eos->data.tabular->pdat, h,
        eos->data.tabular->p_of_h_acc);
    return p;
}

static double eos_rho_of_h_tabular(double h, LALSimNeutronStarEOS * eos)
{
    double rho;
    rho =
        gsl_interp_eval(eos->data.tabular->rho_of_h_interp,
        eos->data.tabular->hdat, eos->data.tabular->rhodat, h,
        eos->data.tabular->rho_of_h_acc);
    return rho;
}

static double eos_h_of_p_tabular(double p, LALSimNeutronStarEOS * eos)
{
    double h;
    h = gsl_interp_eval(eos->data.tabular->h_of_p_interp,
        eos->data.tabular->pdat, eos->data.tabular->hdat, p,
        eos->data.tabular->h_of_p_acc);
    return h;
}

static double eos_dedp_of_p_tabular(double p, LALSimNeutronStarEOS * eos)
{
    double dedp;
    dedp =
        gsl_interp_eval_deriv(eos->data.tabular->e_of_p_interp,
        eos->data.tabular->pdat, eos->data.tabular->edat, p,
        eos->data.tabular->e_of_p_acc);
    return dedp;
}

static double eos_v_of_h_tabular(double h, LALSimNeutronStarEOS * eos)
{
    double p, dedp;
    p = eos_p_of_h_tabular(h, eos);
    dedp = eos_dedp_of_p_tabular(p, eos);
    return pow(dedp, -0.5);
}

//static double eos_v_of_h_tabular(double h, LALSimNeutronStarEOS *eos)
//{
//      double dpdh, dedh;
//    printf("hi4\n");
//    dpdh = gsl_interp_eval_deriv(eos->data.tabular->p_of_h_interp, eos->data.tabular->hdat, eos->data.tabular->pdat, h, eos->data.tabular->p_of_h_acc);
//    printf("hi5\n");
//    dedh = gsl_interp_eval_deriv(eos->data.tabular->e_of_h_interp, eos->data.tabular->hdat, eos->data.tabular->edat, h, eos->data.tabular->e_of_h_acc);
//      return sqrt(dpdh/dedh);
//}

static void eos_free_tabular_data(LALSimNeutronStarEOSDataTabular * data)
{
    if (data) {
        gsl_interp_free(data->e_of_p_interp);
        gsl_interp_free(data->e_of_h_interp);
        gsl_interp_free(data->p_of_h_interp);
        gsl_interp_free(data->h_of_p_interp);
        gsl_interp_free(data->rho_of_h_interp);
        gsl_interp_accel_free(data->e_of_p_acc);
        gsl_interp_accel_free(data->e_of_h_acc);
        gsl_interp_accel_free(data->p_of_h_acc);
        gsl_interp_accel_free(data->h_of_p_acc);
        gsl_interp_accel_free(data->rho_of_h_acc);
        LALFree(data->edat);
        LALFree(data->pdat);
        LALFree(data->hdat);
        LALFree(data->rhodat);
        LALFree(data);
    }
    return;
}

static void eos_free_tabular(LALSimNeutronStarEOS * eos)
{
    if (eos) {
        eos_free_tabular_data(eos->data.tabular);
        LALFree(eos);
    }
    return;
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

    h_im1 = eos->data.tabular->hdat[0];
    v_im1 = eos_v_of_h_tabular(h_im1, eos);
    for (i = 1; i < eos->data.tabular->ndat; i++) {
        h_i = eos->data.tabular->hdat[i];
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

static LALSimNeutronStarEOS *eos_alloc_tabular(double *pdat, double *edat,
    size_t ndat)
{
    LALSimNeutronStarEOS *eos;
    LALSimNeutronStarEOSDataTabular *data;
    size_t i;
    double *hdat;
    double dhdp;
    double *rhodat;
    double integrand_im1, integrand_i, integral;

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
    eos->dedp_of_p = eos_dedp_of_p_tabular;
    eos->v_of_h = eos_v_of_h_tabular;

    /* compute enthalpy data by integrating (trapezoid rule) */
    hdat = LALMalloc(ndat * sizeof(*hdat));
    hdat[0] = 0.0;
    dhdp = 1.0 / (edat[1] + pdat[1]);   /* first deriv is at second point */
    for (i = 1; i < ndat; ++i) {
        double prev = dhdp;
        dhdp = 1.0 / (edat[i] + pdat[i]);
        hdat[i] = hdat[i - 1] + 0.5 * (prev + dhdp) * (pdat[i] - pdat[i - 1]);
    }

    /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
    /*             CALCULATION OF RHO CURRENTLY RETURNS GARBAGE               */
    /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
    /* compute rest-mass density by integrating (trapezoid rule) */
    /* rho_i = rho_{i-1} exp(int_{e_{i-1}}^{e_i} de/(e+p)) */
    rhodat = LALMalloc(ndat * sizeof(*hdat));
    rhodat[0] = 0.0;
    rhodat[1] = edat[1];        /* essentially the same at low density */
    integrand_im1 = 1.0 / (edat[1] + pdat[1]);
    for (i = 2; i < ndat; i++) {
        integrand_i = 1.0 / (edat[i] + pdat[i]);
        integral =
            0.5 * (integrand_im1 + integrand_i) * (edat[i] - edat[i - 1]);
        integrand_im1 = integrand_i;
        rhodat[i] = rhodat[i - 1] * exp(integral);
    }

    data->hdat = hdat;
    data->pdat = pdat;
    data->edat = edat;
    data->rhodat = rhodat;
    data->ndat = ndat;

    eos->pmax = data->pdat[ndat - 1];
    eos->hmax = data->hdat[ndat - 1];

    /* setup interpolation tables */

    data->e_of_p_acc = gsl_interp_accel_alloc();
    data->h_of_p_acc = gsl_interp_accel_alloc();
    data->e_of_h_acc = gsl_interp_accel_alloc();
    data->p_of_h_acc = gsl_interp_accel_alloc();
    data->rho_of_h_acc = gsl_interp_accel_alloc();

    data->e_of_p_interp = gsl_interp_alloc(gsl_interp_cspline, ndat);
    data->h_of_p_interp = gsl_interp_alloc(gsl_interp_cspline, ndat);
    data->e_of_h_interp = gsl_interp_alloc(gsl_interp_cspline, ndat);
    data->p_of_h_interp = gsl_interp_alloc(gsl_interp_cspline, ndat);
    data->rho_of_h_interp = gsl_interp_alloc(gsl_interp_cspline, ndat);

    gsl_interp_init(data->e_of_p_interp, pdat, edat, ndat);
    gsl_interp_init(data->h_of_p_interp, pdat, hdat, ndat);
    gsl_interp_init(data->e_of_h_interp, hdat, edat, ndat);
    gsl_interp_init(data->p_of_h_interp, hdat, pdat, ndat);
    gsl_interp_init(data->rho_of_h_interp, hdat, rhodat, ndat);

    eos->hMinAcausal =
        eos_min_acausal_pseudo_enthalpy_tabular(eos->hmax, eos);

//    printf("%e\n", XLALSimNeutronStarEOSEnergyDensityOfPressureGeometerized(eos->pmax, eos));
//    
//    printf("datatype = %d\n", eos->datatype);
//    printf("pmax = %e\n", eos->pmax);
//    printf("hmax = %e\n", eos->hmax);
//    printf("hMinAcausal = %e\n", eos->hMinAcausal);

    return eos;
}

static int mystrcasecmp(const char *s1, const char *s2)
{
    while (*s1) {
        int c1 = toupper(*s1++);
        int c2 = toupper(*s2++);
        if (c1 != c2)
            return (c1 > c2) - (c1 < c2);
    }
    return 0;
}

/** @endcond */

/**
 * @brief Reads a data file containing a tabulated equation of state.
 * @details Read a data file specified by a path fname that contains two
 * whitespace separated columns of equation of state data.  The first column
 * contains the pressure in Pa and the second column contains the energy
 * density in J/m^3.  Every line beginning with the character '#' then it is
 * ignored.  If the path is an absolute path then this specific file is opened;
 * otherwise, search for the file in paths given in the environment variable
 * LALSIM_DATA_PATH, and finally search in the installed PKG_DATA_DIR path.
 * @param[in] fname The path of the file to open.
 * @return A pointer to neutron star equation of state structure.
 */
LALSimNeutronStarEOS *XLALSimNeutronStarEOSFromFile(const char *fname)
{
    LALSimNeutronStarEOS *eos;
    double *pdat;
    double *edat;
    size_t ndat;
    LALFILE *fp;

    fp = XLALSimReadDataFileOpen(fname);
    if (!fp)
        XLAL_ERROR_NULL(XLAL_EFUNC);

    ndat = XLALSimReadDataFile2Col(&pdat, &edat, fp);
    XLALFileClose(fp);
    if (ndat == (size_t) (-1))
        XLAL_ERROR_NULL(XLAL_EFUNC);

    eos = eos_alloc_tabular(pdat, edat, ndat);
    snprintf(eos->name, sizeof(eos->name), "%s", fname);
    return eos;
}

/**
 * @brief Creates an equation of state structure from tabulated equation
 * of state data of a known name.
 * @details A known, installed, named tabulated equation of state data file is
 * read and the used to create the equation of state structure.  Presently
 * the known equations of state are:
 * - AP4
 * - FPS
 * - SLY4
 * @param[in] name The name of the equation of state.
 * @return A pointer to neutron star equation of state structure.
 */
LALSimNeutronStarEOS *XLALSimNeutronStarEOSByName(const char *name)
{
    static const char fname_base[] = "LALSimNeutronStarEOS_";
    static const char fname_extn[] = ".dat";
    static const char *eos_names[] = {
        "FPS",
        "SLY4",
        "AP4"
    };
    size_t n = sizeof(eos_names) / sizeof(*eos_names);
    size_t i;
    char fname[FILENAME_MAX];

    for (i = 0; i < n; ++i)
        if (mystrcasecmp(name, eos_names[i]) == 0) {
            LALSimNeutronStarEOS *eos;
            snprintf(fname, sizeof(fname), "%s%s%s", fname_base, eos_names[i],
                fname_extn);
            eos = XLALSimNeutronStarEOSFromFile(fname);
            if (!eos)
                XLAL_ERROR_NULL(XLAL_EFUNC);
            snprintf(eos->name, sizeof(eos->name), "%s", eos_names[i]);
            return eos;
        }

    XLAL_PRINT_ERROR("Unrecognized EOS name %s...", name);
    XLALPrintError("\tKnown EOS names are: %s", eos_names[0]);
    for (i = 1; i < n; ++i)
        XLALPrintError(", %s", eos_names[i]);
    XLALPrintError("\n");
    XLAL_ERROR_NULL(XLAL_ENAME);
}

/** @} */
/** @} */
