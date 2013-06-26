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

/* Store data for generating polytropic EOS so you don't have to reevaluate it
 * every time. */
/* Accepts up to 100 polytrope pieces. */
enum {N_POLY_MAX = 100};
struct tagLALSimNeutronStarEOSDataPiecewisePolytrope {
    int nPoly;  /* number of polytrope pieces */
    double rhoTab[N_POLY_MAX]; /* starting rest-mass density of polytropic piece i (kg/m^3) */
    double epsilonTab[N_POLY_MAX];     /* starting energy density of polytropic piece i (J/m^3) */
    double pTab[N_POLY_MAX];   /* starting pressure of polytropic piece i (Pa, J/m^3) */
    double kTab[N_POLY_MAX];   /* polytropic constant of piece i */
    double gammaTab[N_POLY_MAX];       /* adiabatic index of piece i */
    double nTab[N_POLY_MAX];   /* polytropic index n_i = 1/(Gamma_i - 1) */
    double aTab[N_POLY_MAX];   /* integration constant (Eq. 3) of PRD 79, 124032 (2009) */
    double hTab[N_POLY_MAX];   /* dimensionless pseudo-enthalpy = ln(specific enthalpy/c^2) */
};

/* FUNCTIONS OF PSEUDO-ENTHALPY */

/* Determine which polytrope piece you are in. */
/* hTab[i] is starting pseudo-enthalpy of polytropic piece i */
static int polytrope_index_of_h(double h, LALSimNeutronStarEOS * eos)
{
    int i = eos->data.piecewisePolytrope->nPoly - 1;

    while (h <= eos->data.piecewisePolytrope->hTab[i] && i > 0)
        i--;

    return i;
}

/* Rest-mass density as a function of pseudo-enthalpy h */
static double eos_rho_of_h_piecewise_polytrope(double h,
    LALSimNeutronStarEOS * eos)
{
    int i = polytrope_index_of_h(h, eos);       /* index for polytrope piece */
    double enthalpy = exp(h);   /* the real enthalpy */
    double k_i = eos->data.piecewisePolytrope->kTab[i];
    double n_i = eos->data.piecewisePolytrope->nTab[i];
    double a_i = eos->data.piecewisePolytrope->aTab[i];
    double num, den;

    num = enthalpy - 1.0 - a_i;
    den = (n_i + 1.0) * k_i;    /* never 0 */

    return pow(num / den, n_i);
}

/* Pressure as a function of pseudo-enthalpy h */
static double eos_p_of_h_piecewise_polytrope(double h,
    LALSimNeutronStarEOS * eos)
{
    int i = polytrope_index_of_h(h, eos);       /* index for polytrope piece */
    double enthalpy = exp(h);   /* the real enthalpy */
    double k_i = eos->data.piecewisePolytrope->kTab[i];
    double n_i = eos->data.piecewisePolytrope->nTab[i];
    double a_i = eos->data.piecewisePolytrope->aTab[i];
    double num, den;

    num = enthalpy - 1.0 - a_i;
    den = (n_i + 1.0) * k_i;    /* never 0 */

    return k_i * pow(num / den, n_i + 1.0);
}

/* Energy density as a function of pseudo-enthalpy h */
static double eos_e_of_h_piecewise_polytrope(double h,
    LALSimNeutronStarEOS * eos)
{
    int i = polytrope_index_of_h(h, eos);       /* index for polytrope piece */
    double n_i = eos->data.piecewisePolytrope->nTab[i];
    double a_i = eos->data.piecewisePolytrope->aTab[i];
    double enthalpy = exp(h);   /* the real enthalpy */
    double rho = eos_rho_of_h_piecewise_polytrope(h, eos);

    double num = 1.0 + a_i + n_i * enthalpy;
    double den = n_i + 1.0;     /* never 0 */

    return rho * num / den;
}

///* depsilon/dp blows up as h->0 and enthalpy->1 */
//static double eos_dedp_of_h_piecewise_polytrope(double h, LALSimNeutronStarEOS *eos)
//{
//    int i = polytrope_index_of_h(h, eos);  /* index for polytrope piece */
//    double enthalpy = exp(h); /* specific enthalpy */
//    double n_i = eos->data.piecewisePolytrope->nTab[i];
//    double a_i = eos->data.piecewisePolytrope->aTab[i];
//
//    return (n_i*enthalpy)/(enthalpy - 1.0 - a_i);
//}

/* v=sqrt(dp/depsilon)->0 as h->0 and enthalpy->1 */
static double eos_v_of_h_piecewise_polytrope(double h,
    LALSimNeutronStarEOS * eos)
{
    int i = polytrope_index_of_h(h, eos);       /* index for polytrope piece */
    double enthalpy = exp(h);   /* specific enthalpy */
    double n_i = eos->data.piecewisePolytrope->nTab[i];
    double a_i = eos->data.piecewisePolytrope->aTab[i];

    return sqrt((enthalpy - 1.0 - a_i) / (n_i * enthalpy));
}

/* FUNCTIONS OF PRESSURE */

/* Determine which polytrope piece you are in. */
/* pTab[i] is starting pressure of polytropic piece i */
static int polytrope_index_of_p(double p, LALSimNeutronStarEOS * eos)
{
    int i = eos->data.piecewisePolytrope->nPoly - 1;

    while (p <= eos->data.piecewisePolytrope->pTab[i] && i > 0)
        i--;

    return i;
}

static double eos_h_of_p_piecewise_polytrope(double p,
    LALSimNeutronStarEOS * eos)
{
    int i = polytrope_index_of_p(p, eos);       /* index for polytrope piece */
    double k_i = eos->data.piecewisePolytrope->kTab[i];
    double n_i = eos->data.piecewisePolytrope->nTab[i];
    double a_i = eos->data.piecewisePolytrope->aTab[i];

    double enthalpy =
        1.0 + a_i + (n_i + 1.0) * k_i * pow(p / k_i, 1.0 / (n_i + 1.0));

    //    printf("%d\t%e\t%e\t%e\t%e\n", i, k_i, n_i, a_i, log(enthalpy));

    return log(enthalpy);       /* convert from real enthalpy to pseudo-enthalpy h */
}

//static double eos_rho_of_p_piecewise_polytrope(double p, LALSimNeutronStarEOS *eos)
//{
//    int i = polytrope_index_of_p(p, eos);  /* index for polytrope piece */
//    double k_i = eos->data.piecewisePolytrope->kTab[i];
//    double n_i = eos->data.piecewisePolytrope->nTab[i];
//
//    return pow(p/k_i, n_i/(n_i+1.0));
//}

static double eos_e_of_p_piecewise_polytrope(double p,
    LALSimNeutronStarEOS * eos)
{
    double rho;
    int i = polytrope_index_of_p(p, eos);       /* index for polytrope piece */
    double k_i = eos->data.piecewisePolytrope->kTab[i];
    double n_i = eos->data.piecewisePolytrope->nTab[i];
    double a_i = eos->data.piecewisePolytrope->aTab[i];

    rho = pow(p / k_i, n_i / (n_i + 1.0));

    return (1.0 + a_i) * rho + n_i * p;
}

/* depsilon/dp blows up as p->0 */
static double eos_dedp_of_p_piecewise_polytrope(double p,
    LALSimNeutronStarEOS * eos)
{
    double epsilon;
    int i = polytrope_index_of_p(p, eos);       /* index for polytrope piece */
    double gamma_i = eos->data.piecewisePolytrope->gammaTab[i];

    epsilon = eos_e_of_p_piecewise_polytrope(p, eos);

    return (epsilon + p) / (gamma_i * p);       /* returns nan at p=0 */
}

///* v=sqrt(dp/depsilon)->0 as p->0 */
//static double eos_v_of_p_piecewise_polytrope(double p, LALSimNeutronStarEOS *eos)
//{
//    double epsilon;
//    int i = polytrope_index_of_p(p, eos);  /* index for polytrope piece */
//    double gamma_i = eos->data.piecewisePolytrope->gammaTab[i];
//
//    epsilon = eos_e_of_p_piecewise_polytrope(p, eos);
//
//    return sqrt( (gamma_i*p)/(epsilon + p) ); /* currently nan at p=0*/
//}


static void eos_free_piecewise_polytrope(LALSimNeutronStarEOS * eos)
{
    if (eos) {
        LALFree(eos->data.piecewisePolytrope);
        LALFree(eos);
    }
}


/* Minimum pseudo-enthalpy at which EOS becomes acausal (speed of sound > 1).
 * If the EOS is always causal, return some large value hmax instead. */
static double eos_min_acausal_pseudo_enthalpy_piecewise_polytrope(double hmax,
    LALSimNeutronStarEOS * eos)
{
    int i;
    int nPoly = eos->data.piecewisePolytrope->nPoly;
    double n_i, a_i, h_i, h_ip1;
    double vSound;
    double enthalpyMinAcausal;
    double hMinAcausal = hmax;  /* default large number for EOS that is always causal */

    /* For each polytrope piece, speed of sound increases monotonically.    */
    /* If v > 1 at beginning of piece, hMinAcausal = h_i                    */
    /* If v > 1 at end of piece, then find where v = 1 in this piece.       */
    /* If v !> 1 at end of piece, then EOS is always causal in this piece.  */

    for (i = 0; i < nPoly; i++) {
        n_i = eos->data.piecewisePolytrope->nTab[i];
        a_i = eos->data.piecewisePolytrope->aTab[i];
        h_i = eos->data.piecewisePolytrope->hTab[i];

        /* find vSound at beginning of polytrope piece */
        if (i == 0)
            vSound = 0.0;
        else
            vSound = eos_v_of_h_piecewise_polytrope(h_i, eos);
        if (vSound > 1.0) {     /* there is a discontinuity in vSound at h = h_i and it first becomes >1 here */
            hMinAcausal = h_i;
            break;      /* You have hMinAcausal. You're done. Get out of for loop. */
        }

        /* find vSound at end of polytrope piece */
        if (i == nPoly - 1) {
            vSound = 1 / n_i;   /* vSound at h = infinity */
        } else {
            h_ip1 = eos->data.piecewisePolytrope->hTab[i + 1];  /* not defined for i = nPoly-1 */
            vSound = eos_v_of_h_piecewise_polytrope(h_ip1, eos);
        }
        if (vSound > 1) {       /* find the density within this polytrope piece where EOS first becomes acausal */
            enthalpyMinAcausal = (1 + a_i) / (1 - n_i); /* v = 1 here */
            hMinAcausal = log(enthalpyMinAcausal);
            break;      /* You have hMinAcausal. You're done. Get out of for loop. */
        }
    }

    /* Value of h where EOS first becomes acausal. */
    /* Or, if EOS is always causal, hmax */
    return hMinAcausal;
}

#if 0
static void
print_piecewise_polytrope_data(LALSimNeutronStarEOSDataPiecewisePolytrope *
    data)
{
    int i;

    printf("nPoly = %d\n", data->nPoly);

    printf("i =\t\t");
    for (i = 0; i < data->nPoly; i++)
        printf("%d\t\t", i);
    printf("\n");

    printf("rhoTab[i] =\t");
    for (i = 0; i < data->nPoly; i++)
        printf("%e\t", data->rhoTab[i]);
    printf("\n");

    printf("epsilonTab[i] =\t");
    for (i = 0; i < data->nPoly; i++)
        printf("%e\t", data->epsilonTab[i]);
    printf("\n");

    printf("pTab[i] =\t");
    for (i = 0; i < data->nPoly; i++)
        printf("%e\t", data->pTab[i]);
    printf("\n");

    printf("kTab[i] =\t");
    for (i = 0; i < data->nPoly; i++)
        printf("%e\t", data->kTab[i]);
    printf("\n");

    printf("gammaTab[i] =\t");
    for (i = 0; i < data->nPoly; i++)
        printf("%e\t", data->gammaTab[i]);
    printf("\n");

    printf("nTab[i] =\t");
    for (i = 0; i < data->nPoly; i++)
        printf("%e\t", data->nTab[i]);
    printf("\n");

    printf("aTab[i] =\t");
    for (i = 0; i < data->nPoly; i++)
        printf("%e\t", data->aTab[i]);
    printf("\n");

    printf("hTab[i] =\t");
    for (i = 0; i < data->nPoly; i++)
        printf("%e\t", data->hTab[i]);
    printf("\n");
}
#endif

/* GLOBAL FUNCTIONS */

/**
 * @brief Creates a polytrope Equation of State defined by p = K rho^Gamma.
 * @param Gamma Polytrope adiabatic index.
 * @param reference_pressure_si Reference pressure in Pa.
 * @param reference_density_si Density at the reference pressure in kg/m^3.
 * @return A pointer to a newly created EOS structure.
 */
LALSimNeutronStarEOS *XLALSimNeutronStarEOSPolytrope(double Gamma,
    double reference_pressure_si, double reference_density_si)
{
    LALSimNeutronStarEOS *eos;
    LALSimNeutronStarEOSDataPiecewisePolytrope *data;

    /* reference pressure, density in geometric units of m^-2 */
    double p = reference_pressure_si * LAL_G_C4_SI;
    double rho = reference_density_si * LAL_G_C2_SI;

    eos = LALCalloc(1, sizeof(*eos));
    data = LALCalloc(1, sizeof(*data));

    eos->datatype = LALSIM_NEUTRON_STAR_EOS_DATA_TYPE_PIECEWISE_POLYTROPE;
    eos->data.piecewisePolytrope = data;

    /* setup function pointers */
    eos->free = eos_free_piecewise_polytrope;
    eos->e_of_p = eos_e_of_p_piecewise_polytrope;
    eos->h_of_p = eos_h_of_p_piecewise_polytrope;
    eos->e_of_h = eos_e_of_h_piecewise_polytrope;
    eos->rho_of_h = eos_rho_of_h_piecewise_polytrope;
    eos->p_of_h = eos_p_of_h_piecewise_polytrope;
    eos->dedp_of_p = eos_dedp_of_p_piecewise_polytrope;
    eos->v_of_h = eos_v_of_h_piecewise_polytrope;

    /* calculate all the quantities (rho, epsilon, p, k, gamma, n, a, h) in data structure */
    eos->data.piecewisePolytrope->nPoly = 1;
    eos->data.piecewisePolytrope->rhoTab[0] = 0.0;
    eos->data.piecewisePolytrope->epsilonTab[0] = 0.0;
    eos->data.piecewisePolytrope->pTab[0] = 0.0;
    eos->data.piecewisePolytrope->kTab[0] = p / pow(rho, Gamma);
    eos->data.piecewisePolytrope->gammaTab[0] = Gamma;
    eos->data.piecewisePolytrope->nTab[0] = 1.0 / (Gamma - 1.0);
    eos->data.piecewisePolytrope->aTab[0] = 0.0;
    eos->data.piecewisePolytrope->hTab[0] = 0.0;

    /* maximum pressure, pseudo-enthalpy to use for generating mass-radius curves */
    eos->pmax = 10.0 * LAL_NUCLEAR_DENSITY_GEOM_SI;
    eos->hmax = eos_h_of_p_piecewise_polytrope(eos->pmax, eos);
    /* smallest pseudo-enthalpy where EOS becomes acausal */
    eos->hMinAcausal =
        eos_min_acausal_pseudo_enthalpy_piecewise_polytrope(eos->hmax, eos);

#if 0
    print_piecewise_polytrope_data(eos->data.piecewisePolytrope);
    printf("datatype = %d\n", eos->datatype);
    printf("pmax = %e\n", eos->pmax);
    printf("hmax = %e\n", eos->hmax);
    printf("hMinAcausal = %e\n", eos->hMinAcausal);
#endif

    return eos;
}

/* In retrospect, the naming of the dividing densities (rho0, rho1, rho2) and adiabatic indices (gamma1, gamma2, gamma3) are not consistent */
/**
 * @brief Creates a 4-parameter piecewise-polytrope Equation of State.
 * @details A 4-piece piecewise polytrope as described in
 * Physical Review D 79, 124032 (2009) in which the low-density part is
 * fit to the SLY4 EOS.  The remaining pieces are described by adiabatic
 * indices gamma1, gamma2, and gamma3, where gamma1 is the adiabatic index
 * for densities < 10^17.7 kg/m^3, gamma2 is the adiabatic index for densities
 * in the range 10^17.7 kg/m^3 to 10^18 kg/m^3, and gamma3 is the adiabatic
 * index for densities > 10^18 kg/m^3.  The base-10 logarithm of the pressure in
 * Pa at a density of 10^17.7 kg/m^3 is specified by logp1_si.
 * @param logp1_si Base 10 logarithm of the pressure in Pa at the reference
 *                 density of 10^17.7 kg/m^3.
 * @param gamma1 Adiabatic index for densities below 10^17.7 kg/m^3.
 * @param gamma2 Adiabatic index for densities from 10^17.7 kg/m^3
 *               to 10^18 kg/m^3.
 * @param gamma3 Adiabatic index for densities above 10^18 kg/m^3.
 * @return A pointer to a newly created EOS structure.
 */
LALSimNeutronStarEOS *XLALSimNeutronStarEOS4ParameterPiecewisePolytrope(double
    logp1_si, double gamma1, double gamma2, double gamma3)
{
    LALSimNeutronStarEOS *eos;
    LALSimNeutronStarEOSDataPiecewisePolytrope *data;

    /* Data for the 4-piece piecewise polytrope fit to the low-density part of
     * the SLY4 EOS.  Pressure is defined in Pa (N/m^2), and rest-mass density
     * is in kg/m^3.  In the paper PRD 79, 124032 (2009) which used cgs, the
     * k_i values in Table II should be multiplied by c^2. */
    double rhoLow[] = { 0, 2.44033979e10, 3.78358138e14, 2.62780487e15 };
    double kLow[] =
        { 1.0801158752700761e7, 1.311359898998385e10, 6.507604807550857e19,
            3.053461077133694e8 };
    double gammaLow[] = { 1.58424999, 1.28732904, 0.62223344, 1.35692395 };

    /*
     * These are the cgs values:
     * double rhoLow[] = {0, 2.44033979e7, 3.78358138e11, 2.62780487e12};
     * double kLow[] = {6.11252036792443e12, 9.54352947022931e14, 4.787640050002652e22, 3.593885515256112e13};
     * double gammaLow[] = {1.58424999, 1.28732904, 0.62223344, 1.35692395};
     */

    double rho0, rho1, rho2;
    double p1;  /* pressure at 10^17.7 kg/m^3 */
    double p1min;       /* Minimum allowed value of p1 s.t. high density EOS has larger pressure than low density EOS */
    double k1, k2, k3;
    double rhoJoin1, rhoJoin2, pJoin1, pJoin2, gammaJoin, kJoin;        /* constants for extra polytrope if it's needed */
    int i;
    double k_i, rho_i, gamma_i, p_i, n_i, a_i, a_im1, n_im1, epsilon_i,
        enthalpy_i, h_i;

    if (gamma1 <= 1.0 || gamma2 <= 1.0 || gamma3 <= 1.0)
        XLAL_ERROR_NULL(XLAL_EINVAL,
            "Adiabitic indices gamma1=%g, gamma2=%g, and gamma3=%g "
            "must all be greater than 1", gamma1, gamma2, gamma3);

    /* Calculate the variable joining density rho0 between the high and low
     * density EOS */
    rho0 = pow(kLow[3] / k1, 1.0 / (gamma1 - gammaLow[3]));

    /* Transition densities between the 3 high-density polytropes */
    rho1 = pow(10, 17.7);
    rho2 = pow(10, 18.0);

    /* pressure at rho1 */
    p1 = pow(10, logp1_si);

    /* pressure constants */
    k1 = p1 / pow(rho1, gamma1);
    k2 = p1 / pow(rho1, gamma2);
    k3 = k2 * pow(rho2, gamma2 - gamma3);

    p1min = kLow[3] * pow(rho1, gammaLow[3]);
    if (logp1_si < log10(p1min) || logp1_si > 34.5)
        XLAL_ERROR_NULL(XLAL_EINVAL,
            "logp1_si=%g should be between %g and 34.5", logp1_si,
            log10(p1min));

    /* allocate memory */
    eos = LALCalloc(1, sizeof(*eos));
    data = LALCalloc(1, sizeof(*data));

    eos->datatype = LALSIM_NEUTRON_STAR_EOS_DATA_TYPE_PIECEWISE_POLYTROPE;
    eos->data.piecewisePolytrope = data;

    /* setup function pointers */
    eos->free = eos_free_piecewise_polytrope;
    eos->e_of_p = eos_e_of_p_piecewise_polytrope;
    eos->h_of_p = eos_h_of_p_piecewise_polytrope;
    eos->e_of_h = eos_e_of_h_piecewise_polytrope;
    eos->rho_of_h = eos_rho_of_h_piecewise_polytrope;
    eos->p_of_h = eos_p_of_h_piecewise_polytrope;
    eos->dedp_of_p = eos_dedp_of_p_piecewise_polytrope;
    eos->v_of_h = eos_v_of_h_piecewise_polytrope;


    /* Add another polytrope if the joining density is below the start of the
     * last low density polytrope or above the end of the first high density
     * polytrope. */
    if ((rho0 > rhoLow[3]) && (rho0 < rho1)) {
        /* No issue. There will be a total of 7 polytropes. */

        eos->data.piecewisePolytrope->kTab[0] = kLow[0];
        eos->data.piecewisePolytrope->kTab[1] = kLow[1];
        eos->data.piecewisePolytrope->kTab[2] = kLow[2];
        eos->data.piecewisePolytrope->kTab[3] = kLow[3];
        eos->data.piecewisePolytrope->kTab[4] = k1;
        eos->data.piecewisePolytrope->kTab[5] = k2;
        eos->data.piecewisePolytrope->kTab[6] = k3;

        eos->data.piecewisePolytrope->gammaTab[0] = gammaLow[0];
        eos->data.piecewisePolytrope->gammaTab[1] = gammaLow[1];
        eos->data.piecewisePolytrope->gammaTab[2] = gammaLow[2];
        eos->data.piecewisePolytrope->gammaTab[3] = gammaLow[3];
        eos->data.piecewisePolytrope->gammaTab[4] = gamma1;
        eos->data.piecewisePolytrope->gammaTab[5] = gamma2;
        eos->data.piecewisePolytrope->gammaTab[6] = gamma3;

        eos->data.piecewisePolytrope->rhoTab[0] = rhoLow[0];
        eos->data.piecewisePolytrope->rhoTab[1] = rhoLow[1];
        eos->data.piecewisePolytrope->rhoTab[2] = rhoLow[2];
        eos->data.piecewisePolytrope->rhoTab[3] = rhoLow[3];
        eos->data.piecewisePolytrope->rhoTab[4] = rho0;
        eos->data.piecewisePolytrope->rhoTab[5] = rho1;
        eos->data.piecewisePolytrope->rhoTab[6] = rho2;

        eos->data.piecewisePolytrope->nPoly = 7;

    } else {
        /* You have to add an 8th polytrope between gammaLow[3] and gamma1. */
        /* It will be between the densities rhoJoin1 and rhoJoin2. */
        rhoJoin1 = 5.0e15;
        rhoJoin2 = 1.0e16;

        /* Calculate the pressure at the start and end densities. */
        pJoin1 = kLow[3] * pow(rhoJoin1, gammaLow[3]);
        pJoin2 = k1 * pow(rhoJoin2, gamma1);

        /* Calculate K and Gamma for the joining polytrope */
        gammaJoin = log(pJoin2 / pJoin1) / log(rhoJoin2 / rhoJoin1);

        kJoin = pJoin1 / pow(rhoJoin1, gammaJoin);

        /* Now join all 8 polytropes. */
        eos->data.piecewisePolytrope->kTab[0] = kLow[0];
        eos->data.piecewisePolytrope->kTab[1] = kLow[1];
        eos->data.piecewisePolytrope->kTab[2] = kLow[2];
        eos->data.piecewisePolytrope->kTab[3] = kLow[3];
        eos->data.piecewisePolytrope->kTab[4] = kJoin;
        eos->data.piecewisePolytrope->kTab[5] = k1;
        eos->data.piecewisePolytrope->kTab[6] = k2;
        eos->data.piecewisePolytrope->kTab[7] = k3;

        eos->data.piecewisePolytrope->gammaTab[0] = gammaLow[0];
        eos->data.piecewisePolytrope->gammaTab[1] = gammaLow[1];
        eos->data.piecewisePolytrope->gammaTab[2] = gammaLow[2];
        eos->data.piecewisePolytrope->gammaTab[3] = gammaLow[3];
        eos->data.piecewisePolytrope->gammaTab[4] = gammaJoin;
        eos->data.piecewisePolytrope->gammaTab[5] = gamma1;
        eos->data.piecewisePolytrope->gammaTab[6] = gamma2;
        eos->data.piecewisePolytrope->gammaTab[7] = gamma3;

        eos->data.piecewisePolytrope->rhoTab[0] = rhoLow[0];
        eos->data.piecewisePolytrope->rhoTab[1] = rhoLow[1];
        eos->data.piecewisePolytrope->rhoTab[2] = rhoLow[2];
        eos->data.piecewisePolytrope->rhoTab[3] = rhoLow[3];
        eos->data.piecewisePolytrope->rhoTab[4] = rhoJoin1;
        eos->data.piecewisePolytrope->rhoTab[5] = rhoJoin2;
        eos->data.piecewisePolytrope->rhoTab[6] = rho1;
        eos->data.piecewisePolytrope->rhoTab[7] = rho2;

        eos->data.piecewisePolytrope->nPoly = 8;

        XLAL_PRINT_INFO("An extra polytrope was used to join the low and high density regions.");
    }

    /* convert to geometric units and place in piecwisePolytrope structure */
    for (i = 0; i < eos->data.piecewisePolytrope->nPoly; i++) {
        eos->data.piecewisePolytrope->rhoTab[i] *= LAL_G_C2_SI;
        gamma_i = eos->data.piecewisePolytrope->gammaTab[i];
        eos->data.piecewisePolytrope->kTab[i] *=
            pow(LAL_G_SI, 1.0 - gamma_i) * pow((double)LAL_C_SI,
            2.0 * gamma_i - 4.0);
    }

    /* calculate remaining quantities (p, n, a, epsilon, h) in data structure */
    for (i = 0; i < eos->data.piecewisePolytrope->nPoly; i++) {
        k_i = eos->data.piecewisePolytrope->kTab[i];
        rho_i = eos->data.piecewisePolytrope->rhoTab[i];
        gamma_i = eos->data.piecewisePolytrope->gammaTab[i];

        p_i = k_i * pow(rho_i, gamma_i);
        n_i = 1.0 / (gamma_i - 1.0);

        if (i == 0) {
            a_i = 0.0;
        } else {
            a_im1 = eos->data.piecewisePolytrope->aTab[i - 1];
            n_im1 = eos->data.piecewisePolytrope->nTab[i - 1];
            a_i = a_im1 + (n_im1 - n_i) * p_i / rho_i;
        }

        epsilon_i = (1.0 + a_i) * rho_i + n_i * p_i;

        if (i == 0) {
            enthalpy_i = 1.0;   /* p/rho -> 0 as rho -> 0, and a_0 = 0 */
        } else {
            enthalpy_i = 1.0 + a_i + (n_i + 1) * p_i / rho_i;
        }
        h_i = log(enthalpy_i);

        eos->data.piecewisePolytrope->pTab[i] = p_i;
        eos->data.piecewisePolytrope->nTab[i] = n_i;
        eos->data.piecewisePolytrope->aTab[i] = a_i;
        eos->data.piecewisePolytrope->epsilonTab[i] = epsilon_i;
        eos->data.piecewisePolytrope->hTab[i] = h_i;
    }

    eos->pmax = 10.0 * LAL_NUCLEAR_DENSITY_GEOM_SI;
    eos->hmax = eos_h_of_p_piecewise_polytrope(eos->pmax, eos);
    eos->hMinAcausal =
        eos_min_acausal_pseudo_enthalpy_piecewise_polytrope(eos->hmax, eos);

#if 0
    print_piecewise_polytrope_data(eos->data.piecewisePolytrope);
    printf("datatype = %d\n", eos->datatype);
    printf("pmax = %e\n", eos->pmax);
    printf("hmax = %e\n", eos->hmax);
    printf("hMinAcausal = %e\n", eos->hMinAcausal);
#endif

    return eos;
}
