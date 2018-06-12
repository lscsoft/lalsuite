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

void GLBoundConversion(double a, double b, double abcissae[],int nEval);
double AdiabaticIndex(double gamma[],double x, int size);
void resetAbcissae(double abcissae[]);
static double eos_e_of_p_spectral_decomposition(double x, double gamma[],int size, double p0, double e0);
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
    double *rhodat;

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

    /* compute pseudo-enthalpy h from dhdp */
    /* Integrate in log space:
       dhdp = 1 / [e(p) + p]
       h(p) = h(p0) + \int_p0^p dhdp dp
       h(p) = h(p0) + \int_ln(p0)^ln(p) exp[ln(p) + ln(dhdp)] dln(p)
    */
    double * integrand;
    double * log_pdat;
    log_pdat = XLALCalloc(ndat-1, sizeof(*log_pdat));
    integrand = LALMalloc((ndat-1) * sizeof(*integrand));
    for (i = 0; i < ndat-1; ++i) {
        log_pdat[i] = log(pdat[i+1]);
        integrand[i] = exp(log_pdat[i] + log(1.0 / (edat[i+1] + pdat[i+1])));
    }

    gsl_interp_accel * dhdp_of_p_acc_temp = gsl_interp_accel_alloc();
    gsl_interp * dhdp_of_p_interp_temp = gsl_interp_alloc(gsl_interp_linear, ndat-1);
    gsl_interp_init(dhdp_of_p_interp_temp, log_pdat, integrand, ndat-1);

    hdat = LALMalloc(ndat * sizeof(*hdat));
    hdat[0] = 0.0;
    // Do first point by hand
    hdat[1] = hdat[0] + 0.5 * (1./(pdat[1]+edat[1])) * (pdat[1] - pdat[0]);
    for (i = 1; i < ndat-1; ++i) {
        hdat[i+1] = gsl_interp_eval_integ(dhdp_of_p_interp_temp, log_pdat, integrand,
          log_pdat[0], log_pdat[i], dhdp_of_p_acc_temp);
    }
    gsl_interp_free(dhdp_of_p_interp_temp);
    gsl_interp_accel_free(dhdp_of_p_acc_temp);

    LALFree(log_pdat);
    LALFree(integrand);

    // Find rho from e, p, and h
    rhodat = LALMalloc(ndat * sizeof(*hdat));
    for (i=0; i < ndat; i++){
        rhodat[i] = (edat[i]+pdat[i])/exp(hdat[i]);
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
    size_t n = XLAL_NUM_ELEM(eos_names);
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

/* Returns energy density given a pressure and spectral decomposition parameters */
static double eos_e_of_p_spectral_decomposition(double x, double gamma[], int size, double p0, double e0)
{
    // Integration/Placeholder variables
    int i;
    int j;
    int nEval = 10;
    double Integrand;
    double IPrime;
    double e, Mu;

    // Declaring arrays needed for Gauss-Legendre abcissae and weights
    double abcissae[nEval], abcissaePrime[nEval];
    resetAbcissae(abcissae);
    double weights[] = {
       0.0666713443086881,
       0.1494513491505806,
       0.2190863625159820,
       0.2692667193099963,
       0.2955242247147529,
       0.2955242247147529,
       0.2692667193099963,
       0.2190863625159820,
       0.1494513491505806,
       0.0666713443086881
       };


    /* Calculation of Mu */
    Integrand = 0.0;

    GLBoundConversion(0.0,x,abcissae, nEval);
    for(i=0;i<nEval;i++)
    {
      Integrand += weights[i]*pow(AdiabaticIndex(gamma,abcissae[i],size),-1.0);
    }

    Integrand*=(x/2.0);
    Mu = exp(-Integrand);
    /* Calculation of integral in Eq. 7 of PRD 82, 103011 (2010) */
    Integrand = 0.0;

    for(i=0;i<nEval;i++)
    {
      IPrime = 0.0;
      resetAbcissae(abcissaePrime);
      GLBoundConversion(0.0,abcissae[i],abcissaePrime, nEval);

      for(j=0;j<nEval;j++)
      {
        IPrime += weights[j]*pow(AdiabaticIndex(gamma,abcissaePrime[j],size),-1.0);
      }

      IPrime *= (abcissae[i]/2.0);
      IPrime = exp(-IPrime);
      Integrand += weights[i]*(exp(abcissae[i])*IPrime/AdiabaticIndex(gamma,abcissae[i],size));
    }

    Integrand*=(x/2.0);

    e = e0/Mu +(p0/Mu)*Integrand;
    return e;
}

/* Specral decomposition of eos's adiabatic index */
double AdiabaticIndex(double gamma[], double x, int size)
{
        double Gamma, logGamma = 0;
        int i;
        for(i=0;i<size;i++)
        {
          logGamma += gamma[i]*pow(x,(double) i);
        }
        Gamma = exp(logGamma);
        return Gamma;
}


/** Gauss-Legendre quaderature-specific functions **/


/* Maps the abcissae from [-1,1] to [a,b] */
void GLBoundConversion(double a, double b, double abcissae[], int nEval)
{
    int i;

    // Converting to new evalutation points
    for(i=0;i<nEval;i++)
    {
        abcissae[i]=((b-a)/2.0)*abcissae[i] + (a+b)/2.0;
    }
}

/* Resets 10 point array to standard abcissae */
void resetAbcissae(double abcissae[])
{
       abcissae[0] = -0.9739065285171717;
       abcissae[1] = -0.8650633666889845;
       abcissae[2] = -0.6794095682990244;
       abcissae[3] = -0.4333953941292472;
       abcissae[4] = -0.1488743389816312;
       abcissae[5] =  0.1488743389816312;
       abcissae[6] =  0.4333953941292472;
       abcissae[7] =  0.6794095682990244;
       abcissae[8] =  0.8650633666889845;
       abcissae[9] =  0.9739065285171717;
}



LALSimNeutronStarEOS *XLALSimNeutronStarEOSSpectralDecomposition(double gamma[], int size)
{
    LALSimNeutronStarEOS * eos;
    size_t ndat_low = 69;
    size_t ndat = ndat_low + 500;
    size_t i;

    // Minimum pressure and energy density of core EOS geom
    double e0 = 9.54629006e-11;
    double p0 = 4.43784199e-13;

    double xmax = 12.3081;
    double pmax = p0*exp(xmax);

    // Declaring array and allocating memory space
    double *edat;
    double *pdat;
    double xdat[ndat-ndat_low];

    pdat = XLALCalloc(ndat, sizeof(*pdat));
    edat = XLALCalloc(ndat, sizeof(*edat));

    // Low density EOS values to be stitched (SLy, in geom)
    double pdat_low[]={0.00000000e+00,   2.49730009e-31,   1.59235347e-30,
         1.01533235e-29,   6.47406376e-29,   4.12805731e-28,
         2.63217321e-27,   1.67835262e-26,   1.07016799e-25,
         6.82371225e-25,   4.35100369e-24,   1.91523482e-23,
         8.06001537e-23,   3.23144235e-22,   4.34521997e-22,
         1.18566090e-21,   3.16699528e-21,   8.31201995e-21,
         2.15154075e-20,   5.51600847e-20,   7.21972469e-20,
         1.34595234e-19,   2.50269468e-19,   3.41156366e-19,
         4.16096744e-19,   5.66803746e-19,   1.05098304e-18,
         1.94663211e-18,   3.60407863e-18,   4.67819652e-18,
         6.36373536e-18,   8.65904266e-18,   1.17739845e-17,
         1.60126190e-17,   2.06809005e-17,   2.81253637e-17,
         3.82385968e-17,   4.91532870e-17,   6.68349199e-17,
         9.08868981e-17,   1.19805457e-16,   1.23523557e-16,
         1.67975513e-16,   2.14575704e-16,   2.38949918e-16,
         2.71834450e-16,   3.69579177e-16,   4.80543818e-16,
         6.22823125e-16,   6.44883854e-16,   6.51906933e-16,
         6.90079430e-16,   7.51717272e-16,   7.84188682e-16,
         8.12280996e-16,   8.94822824e-16,   1.78030908e-15,
         2.83170525e-15,   4.35257355e-15,   6.44272433e-15,
         9.21014776e-15,   1.72635532e-14,   2.96134301e-14,
         4.76007735e-14,   7.28061891e-14,   1.06973879e-13,
         1.78634067e-13,   3.17897582e-13,   4.16939514e-13};

    double edat_low[]={0.00000000e+00,   9.76800363e-24,   3.08891410e-23,
         9.76800489e-23,   3.08891490e-22,   9.76801000e-22,
         3.08891816e-21,   9.76803078e-21,   3.08893141e-20,
         9.76811525e-20,   3.08898527e-19,   7.75906672e-19,
         1.94909879e-18,   4.89622085e-18,   6.16357861e-18,
         1.22969164e-17,   2.45420997e-17,   4.89823567e-17,
         9.77310869e-17,   1.95024387e-16,   2.45531498e-16,
         3.89264229e-16,   6.16926264e-16,   7.76837038e-16,
         9.00642353e-16,   1.19237569e-15,   1.89037060e-15,
         3.09452823e-15,   4.90635934e-15,   5.96458915e-15,
         7.51061114e-15,   9.79776989e-15,   1.23353809e-14,
         1.55335650e-14,   1.88188196e-14,   2.46271194e-14,
         3.10096915e-14,   3.74392368e-14,   4.91814886e-14,
         6.19454280e-14,   7.62224545e-14,   8.11155415e-14,
         1.05200677e-13,   1.26436151e-13,   1.37076948e-13,
         1.55799897e-13,   2.02900352e-13,   2.47139208e-13,
         3.11281590e-13,   3.19529901e-13,   3.22140999e-13,
         3.36213572e-13,   3.58537297e-13,   3.70113877e-13,
         3.80033593e-13,   4.24821517e-13,   1.57119195e-12,
         2.36083738e-12,   3.33152493e-12,   4.49715211e-12,
         5.87555551e-12,   9.35905600e-12,   1.39898961e-11,
         2.00067862e-11,   2.76245081e-11,   3.69196076e-11,
         5.35224936e-11,   7.81020115e-11,   9.19476188e-11};

    // Populating first 69 points with low density EOS
    for(i=0;i<ndat_low;i++)
    {
      pdat[i]=pdat_low[i];
      edat[i]=edat_low[i];
    }

    // Generating higher density portion of EOS with spectral decomposition
    double logpmax = log(pmax);
    double logp0 = log(p0);
    double dlogp = (logpmax-logp0)/(ndat - ndat_low);

    // Calculating pressure table
    for(i=0;i < ndat-ndat_low;i++)
    {
      pdat[i+ndat_low] = exp(logp0 + dlogp*(i));
    }

    // Calculating xdat table
    for(i = 0;i<ndat - ndat_low;i++) {
      xdat[i] = log(pdat[i+ndat_low]/p0);
    }
    for(i = ndat_low;i < ndat;i++)
    {
      // Calculate energy density for each dimensionless pressure value
      edat[i] = eos_e_of_p_spectral_decomposition(xdat[i-ndat_low], gamma, size, p0,e0);
    }

    /* Updating eos structure */
    eos = eos_alloc_tabular(pdat,edat,ndat);

    if(size == 2)
    {
      snprintf(eos->name, sizeof(eos->name), "2-Parameter Spectral Decomposition (SDgamma0=%g, SDgamma1=%g)", gamma[0], gamma[1]);
    }

    else if(size == 4)
    {
      snprintf(eos->name, sizeof(eos->name), "4-Parameter Spectral Decomposition (SDgamma0=%g, SDgamma1=%g, SDgamma2=%g, SDgamma3=%g)", gamma[0], gamma[1], gamma[2], gamma[3]);
    }

    return eos;
}

// FIXME: Function for using with python
LALSimNeutronStarEOS *XLALSimNeutronStarEOSSpectralDecomposition_for_plot(double SDgamma0, double SDgamma1, double SDgamma2, double SDgamma3, int size){
    double gamma[] = {SDgamma0, SDgamma1, SDgamma2, SDgamma3};
    LALSimNeutronStarEOS * eos;
    eos = XLALSimNeutronStarEOSSpectralDecomposition(gamma, size);
    return eos;
}

/** @} */
/** @} */
