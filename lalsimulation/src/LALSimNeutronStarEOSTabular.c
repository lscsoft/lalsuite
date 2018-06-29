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
#include <lal/LALSimNeutronStar.h>

/** @cond */

/* Contents of the tabular equation of state data structure. */
struct tagLALSimNeutronStarEOSDataTabular {
    double *log_pdat;
    double *log_edat;
    double *log_hdat;
    double *log_rhodat;
    size_t ndat;
    gsl_interp *log_e_of_log_p_interp;
    gsl_interp *log_h_of_log_p_interp;
    gsl_interp *log_e_of_log_h_interp;
    gsl_interp *log_p_of_log_h_interp;
    gsl_interp *log_rho_of_log_h_interp;
    gsl_interp_accel *log_e_of_log_p_acc;
    gsl_interp_accel *log_h_of_log_p_acc;
    gsl_interp_accel *log_e_of_log_h_acc;
    gsl_interp_accel *log_p_of_log_h_acc;
    gsl_interp_accel *log_rho_of_log_h_acc;
};

static double eos_e_of_p_tabular(double p, LALSimNeutronStarEOS * eos)
{
	double log_p;
    double log_e;
	if (p == 0.0)
		return 0.0;
	log_p = log(p);
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
	if (h == 0.0)
		return 0.0;
 	log_h = log(h);
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
	if (h == 0.0)
		return 0.0;
 	log_h = log(h);
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
	if (h == 0.0)
		return 0.0;
 	log_h = log(h);
	if (log_h < eos->data.tabular->log_hdat[0])
		/* use non-relativistic degenerate gas, rho = K * h**(3./2.) */
		return exp(eos->data.tabular->log_rhodat[0] + 1.5 * (log_h - eos->data.tabular->log_hdat[0]));
    log_rho =
        gsl_interp_eval(eos->data.tabular->log_rho_of_log_h_interp,
        eos->data.tabular->log_hdat, eos->data.tabular->log_rhodat, log_h,
        eos->data.tabular->log_rho_of_log_h_acc);
    return exp(log_rho);
}

static double eos_h_of_p_tabular(double p, LALSimNeutronStarEOS * eos)
{
	double log_p;
    double log_h;
	if (p == 0)
		return 0.0;
 	log_p = log(p);
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
	if (p == 0 || (log_p = log(p)) < eos->data.tabular->log_pdat[0])
		/* use non-relativistic degenerate gas, p = K * e**(5./3.) */
		return (3.0 / 5.0) * exp(eos->data.tabular->log_edat[0] - eos->data.tabular->log_pdat[0]);
    log_e = gsl_interp_eval(eos->data.tabular->log_e_of_log_p_interp,
        eos->data.tabular->log_pdat, eos->data.tabular->log_edat, log_p,
        eos->data.tabular->log_e_of_log_p_acc);
    d_log_e_d_log_p =
        gsl_interp_eval_deriv(eos->data.tabular->log_e_of_log_p_interp,
        eos->data.tabular->log_pdat, eos->data.tabular->log_edat, log_p,
        eos->data.tabular->log_e_of_log_p_acc);
    return d_log_e_d_log_p * exp(log_e - log_p);
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
//    dpdh = gsl_interp_eval_deriv(eos->data.tabular->log_p_of_log_h_interp, eos->data.tabular->hdat, eos->data.tabular->pdat, h, eos->data.tabular->p_of_h_acc);
//    printf("hi5\n");
//    dedh = gsl_interp_eval_deriv(eos->data.tabular->e_of_h_interp, eos->data.tabular->hdat, eos->data.tabular->edat, h, eos->data.tabular->e_of_h_acc);
//      return sqrt(dpdh/dedh);
//}

static void eos_free_tabular_data(LALSimNeutronStarEOSDataTabular * data)
{
    if (data) {
        gsl_interp_free(data->log_e_of_log_p_interp);
        gsl_interp_free(data->log_e_of_log_h_interp);
        gsl_interp_free(data->log_p_of_log_h_interp);
        gsl_interp_free(data->log_h_of_log_p_interp);
        gsl_interp_free(data->log_rho_of_log_h_interp);
        gsl_interp_accel_free(data->log_e_of_log_p_acc);
        gsl_interp_accel_free(data->log_e_of_log_h_acc);
        gsl_interp_accel_free(data->log_p_of_log_h_acc);
        gsl_interp_accel_free(data->log_h_of_log_p_acc);
        gsl_interp_accel_free(data->log_rho_of_log_h_acc);
        LALFree(data->log_edat);
        LALFree(data->log_pdat);
        LALFree(data->log_hdat);
        LALFree(data->log_rhodat);
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

static LALSimNeutronStarEOS *eos_alloc_tabular(double *pdat, double *edat,
    size_t ndat)
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
    eos->dedp_of_p = eos_dedp_of_p_tabular;
    eos->v_of_h = eos_v_of_h_tabular;

	/* allocate memory for eos data; ignore first points if 0 */
	while (*pdat == 0.0 || *edat == 0.0) {
		++pdat;
		++edat;
		--ndat;
	}
    data->ndat = ndat;
    data->log_pdat = XLALMalloc(ndat * sizeof(*data->log_pdat));
    data->log_edat = XLALMalloc(ndat * sizeof(*data->log_edat));
    data->log_hdat = XLALMalloc(ndat * sizeof(*data->log_hdat));
    data->log_rhodat = XLALMalloc(ndat * sizeof(*data->log_rhodat));

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

    // Find rho from e, p, and h: rho = (e+p)/exp(h)
    for (i = 0; i < ndat; i++)
        data->log_rhodat[i] = log(edat[i] + pdat[i]) - exp(data->log_hdat[i]);
    
    eos->pmax = exp(data->log_pdat[ndat - 1]);
    eos->hmax = exp(data->log_hdat[ndat - 1]);

    /* setup interpolation tables */

    data->log_e_of_log_p_acc = gsl_interp_accel_alloc();
    data->log_h_of_log_p_acc = gsl_interp_accel_alloc();
    data->log_e_of_log_h_acc = gsl_interp_accel_alloc();
    data->log_p_of_log_h_acc = gsl_interp_accel_alloc();
    data->log_rho_of_log_h_acc = gsl_interp_accel_alloc();

    data->log_e_of_log_p_interp = gsl_interp_alloc(gsl_interp_cspline, ndat);
    data->log_h_of_log_p_interp = gsl_interp_alloc(gsl_interp_cspline, ndat);
    data->log_e_of_log_h_interp = gsl_interp_alloc(gsl_interp_cspline, ndat);
    data->log_p_of_log_h_interp = gsl_interp_alloc(gsl_interp_cspline, ndat);
    data->log_rho_of_log_h_interp = gsl_interp_alloc(gsl_interp_cspline, ndat);

    gsl_interp_init(data->log_e_of_log_p_interp, data->log_pdat, data->log_edat, ndat);
    gsl_interp_init(data->log_h_of_log_p_interp, data->log_pdat, data->log_hdat, ndat);
    gsl_interp_init(data->log_e_of_log_h_interp, data->log_hdat, data->log_edat, ndat);
    gsl_interp_init(data->log_p_of_log_h_interp, data->log_hdat, data->log_pdat, ndat);
    gsl_interp_init(data->log_rho_of_log_h_interp, data->log_hdat, data->log_rhodat, ndat);

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

	XLALFree(pdat);
	XLALFree(edat);

    snprintf(eos->name, sizeof(eos->name), "%s", fname);
    return eos;
}

/**
 * @brief Creates an equation of state structure from tabulated equation
 * of state data of a known name.
 * @details A known, installed, named tabulated equation of state data file is
 * read and the used to create the equation of state structure.  
 * The equations of state available are the representative sample drawn from
 * http://xtreme.as.arizona.edu/NeutronStars/ they are:
 * - ALF1
 * - ALF2
 * - ALF3
 * - ALF4
 * - AP1
 * - AP2
 * - AP3
 * - AP4
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

    for (i = 0; i < n; ++i)
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
    XLALPrintError("\tKnown EOS names are: %s", lalSimNeutronStarEOSNames[0]);
    for (i = 1; i < n; ++i)
        XLALPrintError(", %s", lalSimNeutronStarEOSNames[i]);
    XLALPrintError("\n");
    XLAL_ERROR_NULL(XLAL_ENAME);
}

/** @} */
/** @} */
