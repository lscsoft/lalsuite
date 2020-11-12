
/**
 * Copyright (C) 2017-2018  Alessandro Nagar, Sebastiano Bernuzzi,
 * Sarp Ackay, Gregorio Carullo, Walter Del Pozzo, Ka Wa Tsang, Michalis Agathos
 * LALSimulation implementation by Michalis Agathos
 *
 * This file is part of the LALSimulation version of TEOBResumS.
 * The review of this implementation of the TEOBResumS waveform approximant
 * can be found at https://git.ligo.org/waveforms/reviews/teobresums/wikis/
 */

/**
 * TEOBResumS is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * TEOBResumS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see http://www.gnu.org/licenses/.
 *
 */

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

//#ifndef _OPENMP
//#define omp ignore
//#endif

#include <complex.h>
#include <math.h>

#include <gsl/gsl_const.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>

#include "LALSimTEOBResumS.h"

#define ODE_ABSTOL (1e-13)
#define ODE_RELTOL (1e-11)
#define INTERP_UNIFORM_GRID 1

/**
 * GSL routines for ODE integration
 * https://www.gnu.org/software/gsl/doc/html/ode-initval.html
 * http://www.csse.uwa.edu.au/programming/gsl-1.0/gsl-ref_24.html
 */

/** Global vars, defined as external in header */
const INT4 TEOB_LINDEX[KMAX] = {
    2,2,
    3,3,3,
    4,4,4,4,
    5,5,5,5,5,
    6,6,6,6,6,6,
    7,7,7,7,7,7,7,
    8,8,8,8,8,8,8,8};

const INT4 TEOB_MINDEX[KMAX] = {
    1,2,
    1,2,3,
    1,2,3,4,
    1,2,3,4,5,
    1,2,3,4,5,6,
    1,2,3,4,5,6,7,
    1,2,3,4,5,6,7,8};

/**
 * @addtogroup LALSimIMRTEOBResumS_c
 * @brief Routines to generate time-domain effective-one-body gravitational waveforms for coalescing compact binaries with non-precessing spins, tides and self-spin effects.
 * inspiral-merger-ringdown waveforms.
 *
 * See:
 * A. Nagar et al, TEOBResumS, arXiv:1806.01772, PRD, 98, 104052, (2018) - TEOBResumS main paper
 * A. Nagar et al, arXiv:1805.03891, PRD, 99, 021501, (2019) - post-Adiabatic approximation
 * A. Nagar et al, arXiv:1812.07923, PRD, 99, 044007, (2019) - Nonlinear-in-spin effects.
 * S. Ackay et al, arXiv:1812.02744, PRD, 99, 044051, (2019) - Effective-one-body multipolar waveform for tidally interacting binary neutron stars up to merger
 *
 * @review by Geraint Pratten, Gunnar Riemenschneider, Piero Rettegno, Rachael Huxford, Rossella Gamba
 * Combined review wiki:
 * https://git.ligo.org/waveforms/reviews/teobresums/-/wikis/home
 *
 *
 */

 /**
  * @addtogroup LALSimIMRTEOBResumS_c
  * @{
  *
  * @name Routines for TEOBResumS
  * @{
  *
  * @author Alessandro Nagar, Sebastiano Bernuzzi, Sarp Ackay, Gregorio Carullo, Walter Del Pozzo, Ka Wa Tsang, Michalis Agathos (LALSimulation implementation by Michalis Agathos)
  *
  * @brief C code for TEOBResumS.
  *
  * This is an aligned-spin time-domain model for coalescing compact binaries. The model contains only the 22-mode for BBHs but (optional) inspiral-only higher multipoles (l<8) for BNS and NSBH systems.
  * The model contains equation-of-state specific self-spin interactions (monopole-quadrupole) that are incorporated at leading-order PN. Gravitoelectric tides up to NNLO are modelled and
  * gravitomagnetic tides up to NLO.
  *
  * See:
  * - A. Nagar et al, TEOBResumS, arXiv:1806.01772, PRD, 98, 104052, (2018) - TEOBResumS main paper
  * - A. Nagar et al, arXiv:1805.03891, PRD, 99, 021501, (2019) - post-Adiabatic approximation
  * - A. Nagar et al, arXiv:1812.07923, PRD, 99, 044007, (2019) - Nonlinear-in-spin effects.
  * - S. Ackay et al, arXiv:1812.02744, PRD, 99, 044051, (2019) - Effective-one-body multipolar waveform for tidally interacting binary neutron stars up to merger
  *
  * @note The model was calibrated to NR simulations at mass-ratios 1 to 20.
  *
  * @attention The model is usable outside this parameter range,
  * and in tests to date gives sensible physical results up to mass ratios ~ 30 and spins ~ 0.99.
  * For higher mass ratios and very negative spins memory allocation errors were found.
  * These occur for \f$\eta < 0.073\f$ and \f$S < -0.8\f$. An approximate fit of the region to exclude was found to be: \f$\chi_1\f$ = -51.2 * \f$\eta^2\f$ +2.2456 * \f$\eta\f$ - 0.8804.
  * For more information, see the review wiki https://git.ligo.org/waveforms/reviews/teobresums/-/wikis/home
  *
  *
  */


/**
 *  Driver routine to calculate a TEOBResumS
 *  inspiral-merger-ringdown waveform model
 *  in the time domain.
 *
 *  All input parameters should be in SI units. Angles should be in radians.
 *
 *  Returns the plus and cross polarizations as a complex time series.
 *
 */
int XLALSimIMRTEOBResumS(REAL8TimeSeries **hplus,                    /**< +-polarization waveform */
                         REAL8TimeSeries **hcross,                   /**< x-polarization waveform */
                         const REAL8 phiRef,                         /**< reference orbital phase (rad) */
                         const REAL8 deltaT,                         /**< sampling interval (s) */
                         const REAL8 m1,                             /**< mass of companion 1 (kg) */
                         const REAL8 m2,                             /**< mass of companion 2 (kg) */
                         const REAL8 S1x,                            /**< x-component of the dimensionless spin of object 1 */
                         const REAL8 S1y,                            /**< y-component of the dimensionless spin of object 1 */
                         const REAL8 S1z,                            /**< z-component of the dimensionless spin of object 1 */
                         const REAL8 S2x,                            /**< x-component of the dimensionless spin of object 2 */
                         const REAL8 S2y,                            /**< y-component of the dimensionless spin of object 2 */
                         const REAL8 S2z,                            /**< z-component of the dimensionless spin of object 2 */
                         const REAL8 lambda1,                        /**< (tidal deformation of body 1)/(mass of body 1)^5 */
                         const REAL8 lambda2,                        /**< (tidal deformation of body 2)/(mass of body 2)^5 */
                         const REAL8 distance,                       /**< distance of source (m) */
                         const REAL8 inclination,                    /**< inclination of source (rad) */
                         const REAL8 UNUSED longAscNodes,            /**< longitude of ascending nodes, degenerate with the polarization angle, Omega in documentation */
                         LALDict *LALparams,                         /**< LAL dictionary containing option parameters */
                         const REAL8 UNUSED eccentricity,            /**< eccentrocity at reference epoch */
                         const REAL8 UNUSED meanPerAno,              /**< mean anomaly of periastron */
                         const REAL8 f_min,                          /**< starting GW frequency (Hz) */
                         const REAL8 UNUSED f_ref                    /**< reference GW frequency (Hz) */
)
{
    /* DOMAIN CHECKS */
    XLAL_CHECK(hplus, XLAL_EFAULT, "hplus points to NULL");
    XLAL_CHECK(hcross, XLAL_EFAULT, "hcross points to NULL");
    XLAL_CHECK(*hplus  == NULL, XLAL_EFAULT, "*hplus is not empty.\n");
    XLAL_CHECK(*hcross == NULL, XLAL_EFAULT, "*hcross is not empty.\n");

    XLAL_CHECK(deltaT >= 0.0, XLAL_EDOM, "deltaT cannot take negative values.\n");
    XLAL_CHECK((m1>0.0) && (m2>0.0), XLAL_EDOM, "Masses need to be positive.\n");
    XLAL_CHECK(f_ref >= 0.0, XLAL_EDOM, "Reference frequency f_ref cannot be negative.\n");
    XLAL_CHECK(distance > 0.0, XLAL_EDOM, "Distance needs to be positive");

    /* check for the presence of in-plane spins */
    XLAL_CHECK(((S1x==0)&&(S1y==0)&&(S2x==0)&&(S2y==0)), XLAL_EDOM, "ERROR! Non-aligned spins not supported (yet)! Aborting.\n");


    /* *****************************************
     * Init
     * *****************************************
     */

    REAL8 m1_SI, m2_SI, LambdaAl2, LambdaAl3, LambdaAl4, LambdaBl2, LambdaBl3, LambdaBl4;

    m1_SI = m1;
    m2_SI = m2;
    REAL8 spin1x, spin1y, spin1z, spin2x, spin2y, spin2z;
    spin1x = S1x;
    spin1y = S1y;
    spin1z = S1z;
    spin2x = S2x;
    spin2y = S2y;
    spin2z = S2z;

    /* Tidal polarizability parameters: l=2 are nudged to zero */
    LambdaAl2 = fabs(lambda1) > TEOB_LAMBDA_TOL ? lambda1 : 0.0;
    LambdaAl3 = 0.0;
    LambdaAl4 = 0.0;
    LambdaBl2 = fabs(lambda2) > TEOB_LAMBDA_TOL ? lambda2 : 0.0;
    LambdaBl3 = 0.0;
    LambdaBl4 = 0.0;

    /* l=3,4 tidal parameters are read from LALDict if present */
    if(XLALDictContains(LALparams, "TidalOctupolarLambda1"))
        LambdaAl3 = XLALSimInspiralWaveformParamsLookupTidalOctupolarLambda1(LALparams);
    if(XLALDictContains(LALparams, "TidalOctupolarLambda2"))
        LambdaBl3 = XLALSimInspiralWaveformParamsLookupTidalOctupolarLambda2(LALparams);
    if(XLALDictContains(LALparams, "TidalHexadecapolarLambda1"))
        LambdaAl4 = XLALSimInspiralWaveformParamsLookupTidalHexadecapolarLambda1(LALparams);
    if(XLALDictContains(LALparams, "TidalHexadecapolarLambda2"))
        LambdaBl4 = XLALSimInspiralWaveformParamsLookupTidalHexadecapolarLambda2(LALparams);


    /* Swap parameters if mass convention m1>m2 is not respected */
    if (m2 > m1)
    {
        SWAPTRS(m1_SI, m2_SI);
        SWAPTRS(spin1x, spin2x);
        SWAPTRS(spin1y, spin2y);
        SWAPTRS(spin1z, spin2z);
        SWAPTRS(LambdaAl2, LambdaBl2);
        SWAPTRS(LambdaAl3, LambdaBl3);
        SWAPTRS(LambdaAl4, LambdaBl4);
    }

    /* Switch to mass-rescaled geometric units */
    REAL8 M = (m1_SI + m2_SI) / LAL_MSUN_SI; /* total mass in Msun */
    REAL8 time_unit_fact = time_units_factor(M);

    /* Set useful pars/vars */
    REAL8 q    = m1_SI/m2_SI;
    REAL8 nu   = q_to_nu(q);

    /* HARDCODED 0.5 in geom units */
    const REAL8 dt = 0.5;

    /* Initialize list of modes */
    LALValue *ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(LALparams);

    /* *****************************************
     * Set Memory & do preliminary computations
     * *****************************************
     */

    /* Alloc memory for dynamics and multipolar waveform */
    LALTEOBResumSDynamics *dyn;
    SphHarmPolarTimeSeries *hlm, *this_hlm;
    REAL8Sequence *tdata;
    LIGOTimeGPS epoch0 = LIGOTIMEGPSZERO;
    LALTEOBResumSWaveformModeSingleTime *hlm_t=NULL;
    SphHarmPolarTimeSeries *hlm_nqc=NULL;

    SphHarmPolarTimeSeries *hlm_mrg=NULL; /* merger chunk */
    LALTEOBResumSDynamics *dyn_mrg=NULL;

    const INT4 chunk = 500;
    INT4 size = chunk; /* note: size can vary */

    INT4 store_dynamics = 0;

    INT4 use_postadiab_dyn = 1;
    INT4 rush_only = 0;

    /* Parse usage flags for postadiabatic (rush) speed-up */
    if (XLALDictContains(LALparams, "TEOB_use_postadiabatic"))
    {
        switch (XLALDictLookupINT4Value(LALparams, "TEOB_use_postadiabatic"))
        {
            case TEOB_PA_ON:
                break;
            case TEOB_PA_OFF:
                use_postadiab_dyn = 0;
                break;
            case TEOB_PA_ONLY:
                rush_only = 1;
                break;
            default: XLAL_ERROR(XLAL_EINVAL, "Postadiabatic usage option not recognized.");
        }
    }

    /* Define lower radius limit for PA */
    LALSimInspiralTidalOrder tideO = XLALSimInspiralWaveformParamsLookupPNTidalOrder(LALparams);
    REAL8 pa_rmin = tideO == LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN ? (REAL8) POSTADIABATIC_RMIN_BBH : (REAL8) POSTADIABATIC_RMIN_BNS;

    /* Compute initial radius (NOTE: this may change) */
    const REAL8 f0M = f_min/time_unit_fact;
    REAL8 r0byM = eob_dyn_r0_Kepler(f0M);
    //REAL8 r0byM = eob_dyn_r0_eob(f0M, dyn); /* Radius from EOB equations. This is what should be used. */
    // if (DEBUG) fprintf(stdout, "r0 = %f\n", r0byM);

    /* If f_min is too high fall back to a minimum acceptable initial radius */
    if (r0byM < TEOB_R0_THRESHOLD) {
        r0byM = TEOB_R0_THRESHOLD;
    }

    size =  floor((r0byM - pa_rmin)/POSTADIABATIC_DR) + 1;

    /* If initial radius is too close to PA limit then skip PA and go directly to ODE */
    if (size - 1 < POSTADIABATIC_NSTEP_MIN) {
        size = chunk;
        use_postadiab_dyn = 0;
    }

    /*
    XLAL_CHECK(size>1, XLAL_EINVAL, "Start of waveform at r < %f. Exiting...\n", pa_rmin);
    if (DEBUG)
    {
        fprintf(stdout, "size = %d\n", size);
        fprintf(stdout, "r0 = %f\n", r0byM);
    }
    */

    /* Initialize dynamics structure */
    XLALTEOBDynamicsInit(&dyn, size, "dyn");

    /* Set quick-access parameters dyn (be careful here) */
    XLALTEOBDynamicsSetParams(dyn, &ModeArray, m1_SI, m2_SI, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, dt, LambdaAl2, LambdaBl2, LambdaAl3, LambdaBl3, LambdaAl4, LambdaBl4, LALparams);
    dyn->store = dyn->noflux = 0; /* Default: do not store vars, flux on */

    /* Set some flags (to be later controlled by LALDict) */
    INT4 use_tidal, use_spins;
    INT4 merger_interp = 0;
    use_spins = dyn->use_spins;
    use_tidal = dyn->use_tidal;
    INT4 interp_uniform_grid = INTERP_UNIFORM_GRID_DEFAULT;

    /* Use rk45 by default for ODE */
    INT4 use_rk45 = 1;

    /* Parse interpolation scheme flag from LALDict */
    if (XLALDictContains(LALparams, "TEOB_interp_scheme"))
    {
        interp_uniform_grid = XLALDictLookupINT4Value(LALparams, "TEOB_interp_scheme");
    }

    /* Parse ODE integrator flag from LALDict */
    if (XLALDictContains(LALparams, "TEOB_use_rk45"))
    {
        use_rk45 = XLALDictLookupINT4Value(LALparams, "TEOB_use_rk45");
    }

    const INT4 ode_tstep = dyn->ode_timestep;
    REAL8 t_peak = 0.0;

    store_dynamics = (use_postadiab_dyn || !use_tidal);

    /* Amplitude and phase TimeSeries */
    REAL8TimeSeries *Alm_adt, *philm_adt;
    LIGOTimeGPS gpst0 = LIGOTIMEGPSZERO;

    Alm_adt   = XLALCreateREAL8TimeSeries("A_lm",   &gpst0, 0, deltaT, &lalStrainUnit, (size_t) size);
    philm_adt = XLALCreateREAL8TimeSeries("phi_lm", &gpst0, 0, deltaT, &lalDimensionlessUnit, (size_t) size);
    XLAL_CHECK (Alm_adt && philm_adt, XLAL_ENOMEM, "Could not allocate memory for hlm data.\n");

    /* Time series sequence */
    tdata = XLALCreateREAL8Sequence(size);
    XLAL_CHECK (tdata, XLAL_ENOMEM, "Could not allocate memory for time data.\n");

    /* Populate the mode list with head pointing at l=2, m=1 mode */
    hlm = NULL;
    for (int k=KMAX-1; k>=0; k--)
    {
        hlm = XLALSphHarmPolarTimeSeriesAddMode(hlm, Alm_adt, philm_adt, TEOB_LINDEX[k], TEOB_MINDEX[k]);
        XLAL_CHECK(hlm, XLAL_ENOMEM, "Could not allocate hlm mode.\n");
    }
    XLALSphHarmPolarTimeSeriesSetTData(hlm, tdata); // MA: tdata still needs to be filled in
    Waveform_lm_t_alloc (&hlm_t);
    XLALDestroyREAL8TimeSeries(Alm_adt);
    XLALDestroyREAL8TimeSeries(philm_adt);

    /* Set r.h.s. fun pointer */
    int (*p_eob_dyn_rhs)(REAL8, const REAL8*, REAL8*, void*);
    if (use_spins) p_eob_dyn_rhs = &eob_dyn_rhs_s;
    else           p_eob_dyn_rhs = &eob_dyn_rhs;

    /* Iteration index */
    INT4 iter = 0;

    if (use_postadiab_dyn)
    {
        /* *****************************************
         * Post-adiabatic dynamics
         * *****************************************
         */

        /* Calculate dynamics */
        eob_dyn_Npostadiabatic(dyn, r0byM);

        /* Calculate waveform */
        for (int i = 0; i < size; i++) tdata->data[i] = dyn->time[i];

        dyn->store = dyn->noflux = 1;

        for (int i = 0; i < size; i++)
        {
            dyn->y[TEOB_EVOLVE_RAD]    = dyn->data[TEOB_RAD][i];
            dyn->y[TEOB_EVOLVE_PHI]    = dyn->data[TEOB_PHI][i];
            dyn->y[TEOB_EVOLVE_PRSTAR] = dyn->data[TEOB_PRSTAR][i];
            dyn->y[TEOB_EVOLVE_PPHI]   = dyn->data[TEOB_PPHI][i];
            p_eob_dyn_rhs(dyn->t, dyn->y, dyn->dy, dyn);
            eob_wav_hlm(dyn, hlm_t);

            this_hlm = hlm;
            for (int k = 0; k < KMAX; k++)
            {
                if (DEBUG)
                {
                    XLAL_CHECK(((int) this_hlm->l == TEOB_LINDEX[k]) && ((int) this_hlm->m == TEOB_MINDEX[k]), XLAL_EDATA, "Mode numbers do not match\n");
                }
                XLAL_CHECK(this_hlm, XLAL_EFAULT, "Missing modes.\n");
                this_hlm->ampl->data->data[i] = hlm_t->ampli[k];
                this_hlm->phase->data->data[i] = hlm_t->phase[k];
                this_hlm = this_hlm->next;
            }
            XLAL_CHECK(this_hlm==NULL, XLAL_EFAULT, "More modes present than expected.\n");
        }

        dyn->store = dyn->noflux = 0;
        if (rush_only)
        {
            /* SKIP ODE EVOLUTION */
            goto END_ODE_EVOLUTION;
        }

        /* Prepare for evolution */

        /* start counting from here */
        iter = size-1;
        dyn->dt = 0.5*(dyn->time[iter]-dyn->time[iter-1]);

        /* Set arrays with initial conditions
         Note current time is already set in dyn->t */
        dyn->y0[TEOB_ID_RAD]  = dyn->r;
        dyn->y0[TEOB_ID_PHI]  = dyn->phi;
        dyn->y0[TEOB_ID_PPHI] = dyn->pphi;
        dyn->y0[TEOB_ID_OMGJ] = dyn->Omg;
        dyn->y0[TEOB_ID_PRSTAR] = dyn->prstar;
        dyn->y[TEOB_EVOLVE_RAD]    = dyn->r;
        dyn->y[TEOB_EVOLVE_PHI]    = dyn->phi;
        dyn->y[TEOB_EVOLVE_PRSTAR] = dyn->prstar;
        dyn->y[TEOB_EVOLVE_PPHI]   = dyn->pphi;

    }
    else
    {
        /* *****************************************
         * Initial conditions for the evolution
         * *****************************************
         */

        /* Compute the initial conditions */
        if (use_spins) eob_dyn_ic_s(r0byM, dyn, dyn->y0);
        else           eob_dyn_ic(r0byM, dyn, dyn->y0);

        /* Se arrays with initial conditions */
        dyn->t       = 0.;
        dyn->r       = dyn->y0[TEOB_ID_RAD];
        dyn->phi     = 0.;
        dyn->pphi    = dyn->y0[TEOB_ID_PPHI];
        dyn->Omg     = dyn->y0[TEOB_ID_OMGJ];
        dyn->ddotr   = 0.;
        dyn->prstar  = dyn->y0[TEOB_ID_PRSTAR];
        dyn->Omg_orb = 0.; //FIXME
        dyn->y[TEOB_EVOLVE_RAD]    = dyn->r;
        dyn->y[TEOB_EVOLVE_PHI]    = dyn->phi;
        dyn->y[TEOB_EVOLVE_PRSTAR] = dyn->prstar;
        dyn->y[TEOB_EVOLVE_PPHI]   = dyn->pphi;
        if (store_dynamics)
        {
            dyn->time[0]             = dyn->t;
            dyn->data[TEOB_RAD][0]    = dyn->r;
            dyn->data[TEOB_PHI][0]    = dyn->phi;
            dyn->data[TEOB_PPHI][0]   = dyn->pphi;
            dyn->data[TEOB_MOMG][0]   = dyn->Omg;
            dyn->data[TEOB_DDOTR][0]  = dyn->ddotr;
            dyn->data[TEOB_PRSTAR][0] = dyn->prstar;
            dyn->data[TEOB_OMGORB][0] = dyn->Omg_orb;
            dyn->data[TEOB_E0][0]     = dyn->E;
        }

        /* Waveform computation at t = 0
         Needs a r.h.s. evaluation for some vars (no flux) */
        dyn->store = dyn->noflux = 1;
        p_eob_dyn_rhs(dyn->t, dyn->y, dyn->dy, dyn);
        dyn->store = dyn->noflux = 0;
        eob_wav_hlm(dyn, hlm_t);

        /* Append waveform to mode list */
        this_hlm = hlm;
        hlm->tdata->data[0] = 0.;
        for (int k = 0; k < KMAX; k++)
        {
            if ( (DEBUG) && (this_hlm->l -TEOB_LINDEX[k] != 0 || this_hlm->m != TEOB_MINDEX[k]) )
            {
                XLAL_ERROR(XLAL_EDATA, "Mode numbers do not match!");
            }
            XLAL_CHECK(this_hlm, XLAL_EFAULT, "Missing modes.\n");
            this_hlm->ampl->data->data[0] = hlm_t->ampli[k];
            this_hlm->phase->data->data[0] = hlm_t->phase[k];
            this_hlm = this_hlm->next;
        }
        XLAL_CHECK(this_hlm==NULL, XLAL_EFAULT, "More modes present than expected.\n");

        /* Prepare for evolution */
        dyn->dt = dt;

    }


    /* *****************************************
     * ODE Evolution
     * *****************************************
     */

    /* Initialize ODE system solver */
    dyn->ode_stop          = false;
    dyn->ode_stop_MOmgpeak = false;
    dyn->ode_stop_radius   = false;

    // Currently hardcoded to 1.0
    const REAL8 rstop = dyn->r_stop;
    if (rstop>0.)
    {
        dyn->ode_stop_radius   = true;
    }

    /* GSL integrator memory */
    gsl_odeiv2_system sys          = {p_eob_dyn_rhs, NULL , TEOB_EVOLVE_NVARS, dyn};
    const gsl_odeiv2_step_type * T;
    gsl_odeiv2_driver * d;

    if (use_rk45)
    {
        T = gsl_odeiv2_step_rkf45;
        d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, (double) dyn->dt, ODE_ABSTOL, ODE_RELTOL);
    }
    else
    {
        // If we prefer rk8:
        T = gsl_odeiv2_step_rk8pd;
        d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, dyn->dt, ODE_ABSTOL, ODE_RELTOL);
    }

    gsl_odeiv2_step * s            = gsl_odeiv2_step_alloc (T, TEOB_EVOLVE_NVARS);
    gsl_odeiv2_control * c         = gsl_odeiv2_control_y_new (ODE_ABSTOL, ODE_RELTOL);
    gsl_odeiv2_evolve * e          = gsl_odeiv2_evolve_alloc (TEOB_EVOLVE_NVARS);

    /* Set optimized dt around merger */
    const double dt_tuned_mrg = get_mrg_timestep(q, dyn->chi1, dyn->chi2);

    /* Solve ODE */
    int STATUS = 0;
    while (!(dyn->ode_stop))
    {
        iter++;

        switch (ode_tstep)
        {
            case ODE_TSTEP_UNIFORM:
            {
                /* Uniform timestepping */
                dyn->ti = dyn->t + dyn->dt;
                STATUS = gsl_odeiv2_driver_apply (d, &dyn->t, (double) dyn->ti, dyn->y);
                break;
            }
            case ODE_TSTEP_ADAPTIVE:
            {
                /* Adaptive timestepping */
                if ( dyn->ode_stop_MOmgpeak == true )
                /* if we are after the peak, slow down and fix the last steps ! */
                    STATUS = gsl_odeiv2_evolve_apply_fixed_step (e, c, s, &sys, &dyn->t, dyn->dt, dyn->y);
                else
                    STATUS = gsl_odeiv2_evolve_apply (e, c, s, &sys, &dyn->t, (double) dyn->t_stop, &dyn->dt, dyn->y);
                break;
            }
            case ODE_TSTEP_ADAPTIVE_UNIFORM_AFTER_LSO:
            {
                /* Adaptive timestepping until LSO ... */
                if (dyn->r > dyn->rLSO)
                {
                    STATUS = gsl_odeiv2_evolve_apply (e, c, s, &sys, &dyn->t, (double) dyn->t_stop, &dyn->dt, dyn->y);
                }
                else
                {
                    /* ... uniform afterwards */
                    dyn->dt = dt_tuned_mrg;
                    dyn->ti = dyn->t + dyn->dt;
                    STATUS = gsl_odeiv2_evolve_apply_fixed_step (e, c, s, &sys, &dyn->t, dyn->dt, dyn->y);
                }
                break;
            }
            default: XLAL_ERROR(XLAL_EINVAL, "Unknown ODE time step method. Aborting\n");
        }

        /* Check for failures */

        if (dyn->ode_stop_MOmgpeak == true)
        {
            /* ...if after the Omega_orb peak, stop integration */
            if ( (STATUS != GSL_SUCCESS) || (!isfinite(dyn->y[TEOB_EVOLVE_RAD])) )
            {
                iter--; /* do count this iter! */
                dyn->ode_stop = true;
                break; /* (while) stop */
            }
        }

        XLAL_CHECK(STATUS == GSL_SUCCESS, XLAL_EFAULT | STATUS, "GSL ODE solver failed. Error = %d\n", STATUS);

        /* Checking whether the dynamics produces NaN values this can happen if radius r becomes too small */
        XLAL_CHECK(isfinite(dyn->r), XLAL_ESING, "ODE solver returned NaN radius.\n");

        /* Unpack data */
        dyn->r      = dyn->y[TEOB_EVOLVE_RAD];
        dyn->phi    = dyn->y[TEOB_EVOLVE_PHI];
        dyn->prstar = dyn->y[TEOB_EVOLVE_PRSTAR];
        dyn->pphi   = dyn->y[TEOB_EVOLVE_PPHI];

        /* Waveform computation
         Needs a r.h.s. evaluation for some vars (but no flux) */
        dyn->store = dyn->noflux = 1;
        p_eob_dyn_rhs(dyn->t, dyn->y, dyn->dy, dyn);
        dyn->store = dyn->noflux = 0;
        eob_wav_hlm(dyn, hlm_t);

        /* Update size and push arrays (if needed) */
        if (iter==size)
        {
            size += chunk;
            /* Resize time sequence, padding with 0s */
            XLALResizeREAL8Sequence(tdata, 0, size);

            /* Resize mode series and dynamics */
            XLALResizeSphHarmPolarTimeSeries(hlm, 0, size);
            XLALTEOBDynamicsPush(&dyn, size);

            XLAL_CHECK(tdata && hlm && dyn, XLAL_ENOMEM, "Could not resize waveform structures.\n");
        }

        /* Append waveform and dynamics to mode list */
        tdata->data[iter] = hlm_t->time;
        this_hlm = hlm;
        for (int k = 0; k < KMAX; k++)
        {
            if ( (DEBUG) && (this_hlm->l - TEOB_LINDEX[k] != 0 || this_hlm->m != TEOB_MINDEX[k]) )
            {
                XLAL_ERROR(XLAL_EDATA, "Mode numbers do not match!");
            }
            XLAL_CHECK(this_hlm, XLAL_EFAULT, "Missing modes.\n");
            this_hlm->ampl->data->data[iter] = hlm_t->ampli[k];
            this_hlm->phase->data->data[iter] = hlm_t->phase[k];
            this_hlm = this_hlm->next;
        }
        XLAL_CHECK(this_hlm==NULL, XLAL_EFAULT, "More modes present than expected.\n");

        if (store_dynamics)
        {
            dyn->time[iter]              = dyn->t;
            dyn->data[TEOB_RAD][iter]    = dyn->r;
            dyn->data[TEOB_PHI][iter]    = dyn->phi;
            dyn->data[TEOB_PPHI][iter]   = dyn->pphi;
            dyn->data[TEOB_MOMG][iter]   = dyn->Omg;
            dyn->data[TEOB_DDOTR][iter]  = dyn->ddotr;
            dyn->data[TEOB_PRSTAR][iter] = dyn->prstar;
            dyn->data[TEOB_OMGORB][iter] = dyn->Omg_orb;
            dyn->data[TEOB_E0][iter]     = dyn->E;
        }

        /* Stop integration if reached max time */
        if (dyn->t >= dyn->t_stop) dyn->ode_stop = true;

        /* Stop integration at given radius (if rstop >= 0) */
        if ((dyn->ode_stop_radius) && (dyn->r < rstop) ) dyn->ode_stop = true;

        /* Check when to break the computation
         find peak of omega curve and continue for 2M */
        dyn->MOmg = use_spins ? dyn->Omg_orb : dyn->Omg;

        if (dyn->ode_stop_MOmgpeak == false)
        {
            /* Before the Omega_orb peak */
            if (dyn->MOmg < dyn->MOmg_prev)
            /* This is the first step after the peak
             Set things for uniform tstep evolution */
            {
                dyn->tMOmgpeak = dyn->t; // = dyn->t-0.5*dyn->dt;
                dyn->ode_stop_MOmgpeak = true;
                dyn->dt = MIN(dyn->dt, dt_tuned_mrg);
                // dyn->t_stop = dyn->t + nstep_stop*dyn->dt; // continue for nstep_stop iters
                dyn->t_stop = dyn->t + 2.;
                t_peak = dyn->t;
            }
            else
            {
                /* Peak not reached, update the max */
                dyn->MOmg_prev = dyn->MOmg;
            }
        }
        else
        {
            if (dyn->t >= dyn->t_stop) dyn->ode_stop = true;
        }

    }
    /* end time iteration */

    /* Free ODE system solver */
    gsl_odeiv2_evolve_free (e);
    gsl_odeiv2_control_free (c);
    gsl_odeiv2_step_free (s);
    gsl_odeiv2_driver_free (d);

    /* Update waveform and dynamics size
     resize to actual size */
    size = iter+1;

    /* Size down mode series and time sequence */
    XLALResizeSphHarmPolarTimeSeries(hlm, 0, size);
    XLALResizeREAL8Sequence(hlm->tdata, 0, size);

    /* Size down dynamics */
    XLALTEOBDynamicsPush(&dyn, size);

END_ODE_EVOLUTION:;

    if (!(use_tidal) && !(rush_only))
    {

        /* *****************************************
         * Following is for BBH : NQC & Ringdown
         * *****************************************
         */

        /*
            This is a BBH run. NQC and ringdown attachment currently assume uniform grids.
            Do we need to interpolate ?
        */

        /* In general, yes interpolate... */
        merger_interp = 1;

        /* HARDCODED default value */
        REAL8 dt_merger_interp = 0.5;

        /* ... except if merger is covered by uniform tstep, then we do NOT interpolate */
        if (ode_tstep != ODE_TSTEP_ADAPTIVE) merger_interp = 0;

        if (merger_interp)
        {
            /* NQC and ringdown attachment is done around merger
             using auxiliary variables defined around [tmin,tmax]
             Recall that parameters are NOT stored into these auxiliary vars */

            const REAL8 tmin = hlm->tdata->data[size-1] - 20; /* Use last 20M points */
            const REAL8 tmax = hlm->tdata->data[size-1] + 2*dt; /* Make sure to use or get last point */

            /* Find indexes of closest elements to  (to, tn) */
            INT4 imin = 0;
            INT4 imax = hlm->tdata->length - 1;
            /* Check limits first */
            XLAL_CHECK((tmin<tmax) && (tmin < hlm->tdata->data[imax]) && (tmax > hlm->tdata->data[imin]), XLAL_EDOM, "Bad choice of times.");

            if (tmin > hlm->tdata->data[0])
                imin = find_point_bisection(tmin, hlm->tdata->length, hlm->tdata->data, 1);
            if (tmax < hlm->tdata->data[imax])
                imax = find_point_bisection(tmax, hlm->tdata->length, hlm->tdata->data, 0);

            /* Create a new mode list and dynamics by extracting a part of the original */
            hlm_mrg = XLALCutSphHarmPolarTimeSeries(hlm, imin, (size_t) imax - imin);
            XLALTEOBDynamicsExtract (dyn, tmin, tmax, &dyn_mrg, "dyn_mrg");

            /*  Interpolate mrg on uniform grid */

            /* Build uniform grid of width dt and alloc tmp memory */
            //dt_merger_interp = MIN(dt_merger_interp, (dyn->time[size-1] - dyn->tMOmgpeak)/4 ); /* Make sure to have always 3 points */
            dt_merger_interp = MIN(dt_merger_interp, dyn->dt);

            const INT4 size_mrg = get_uniform_size(hlm_mrg->tdata->data[hlm_mrg->tdata->length-1], hlm_mrg->tdata->data[0], dt_merger_interp);

            /* Fill new time array */
            REAL8Sequence *mrgtime_interp;
            mrgtime_interp = XLALCreateREAL8Sequence(size_mrg);
            for (int i = 0; i < size_mrg; i++)
                mrgtime_interp->data[i] = i*dt_merger_interp + hlm_mrg->tdata->data[0];

            /* Interp Waveform */

            REAL8TimeSeries *Alm_mrg, *philm_mrg;
            this_hlm = hlm_mrg;
            while (this_hlm)
            {
                /* Interpolate mode amplitude and replace */
                Alm_mrg = XLALCreateREAL8TimeSeries(this_hlm->ampl->name, &epoch0, 0, dt_merger_interp, &(lalStrainUnit), (size_t) size_mrg);
                interp_spline(this_hlm->tdata->data, this_hlm->ampl->data->data, this_hlm->ampl->data->length, mrgtime_interp->data, size_mrg, Alm_mrg->data->data);
                XLALDestroyREAL8TimeSeries(this_hlm->ampl);
                this_hlm->ampl = Alm_mrg;

                /* Interpolate mode phase and replace */
                philm_mrg = XLALCreateREAL8TimeSeries(this_hlm->phase->name, &epoch0, 0, dt_merger_interp, &(lalDimensionlessUnit), (size_t) size_mrg);
                interp_spline(this_hlm->tdata->data, this_hlm->phase->data->data, this_hlm->phase->data->length, mrgtime_interp->data, size_mrg, philm_mrg->data->data);
                XLALDestroyREAL8TimeSeries(this_hlm->phase);
                this_hlm->phase = philm_mrg;

                Alm_mrg = NULL;
                philm_mrg = NULL;
                this_hlm = this_hlm->next;
            }

            /* Replace time sequence */
            if (hlm_mrg->tdata) XLALDestroyREAL8Sequence(hlm_mrg->tdata);
            XLALSphHarmPolarTimeSeriesSetTData(hlm_mrg, mrgtime_interp);
            mrgtime_interp = NULL;
            LALFree(mrgtime_interp);

            /* Interp Dynamics */
            XLALTEOBDynamicsInterp (dyn_mrg, size_mrg, dyn_mrg->time[0], dt_merger_interp, "dyn_mrg_interp");

        } /* End of merger interp */

        if (dyn->nqc_coeffs_hlm == NQC_COMPUTE) {

            /* BBH : compute and add NQC */

            UINT4 size_nqc = merger_interp ? hlm_mrg->tdata->length : (UINT4) size;

            /* Populate NQC mode list */
            for (int k=KMAX-1; k>=0; k--)
            {
                REAL8TimeSeries *Alm, *philm;
                Alm = XLALCreateREAL8TimeSeries("A_nqc", &epoch0, 0, deltaT, &(lalStrainUnit), (size_t) size_nqc);
                philm = XLALCreateREAL8TimeSeries("phi_nqc", &epoch0, 0, deltaT, &(lalDimensionlessUnit), (size_t) size_nqc);
                hlm_nqc = XLALSphHarmPolarTimeSeriesAddMode(hlm_nqc, Alm, philm, TEOB_LINDEX[k], TEOB_MINDEX[k]);
                XLAL_CHECK(hlm_nqc, XLAL_ENOMEM, "Could not allocate hlm_nqc mode.\n");
                XLALDestroyREAL8TimeSeries(Alm);
                XLALDestroyREAL8TimeSeries(philm);
            }


            if (merger_interp) {

                /* Compute NQC only around merger,
                 add to both merger and full waveform */
                XLALSphHarmPolarTimeSeriesSetTData(hlm_nqc, hlm_mrg->tdata);
                eob_wav_hlmNQC_find_a1a2a3_mrg(dyn_mrg, hlm_mrg, hlm_nqc, dyn, hlm);
                this_hlm = hlm_mrg;
                while (this_hlm) {
                    strcat(this_hlm->ampl->name,"_nqc");
                    strcat(this_hlm->phase->name,"_nqc");
                    this_hlm = this_hlm->next;
                }

                /* Join merger to full waveform */
                XLALSphHarmPolarJoin(hlm, hlm_mrg, hlm_mrg->tdata->data[0]);
                XLALTEOBDynamicsJoin (dyn, dyn_mrg, dyn_mrg->time[0]);
                size = hlm->tdata->length;

            }
            else
            {

                /* Compute NQC and add them to full waveform */
                XLALSphHarmPolarTimeSeriesSetTData(hlm_nqc, hlm->tdata);
                eob_wav_hlmNQC_find_a1a2a3(dyn, hlm, hlm_nqc);

            }
            /* Suffix to hlm name */
            this_hlm = hlm;
            while (this_hlm)
            {
                strcat(this_hlm->ampl->name,"_nqc");
                strcat(this_hlm->phase->name,"_nqc");
                this_hlm = this_hlm->next;
            }

        }

        /* BBH : add Ringdown */

        /* Extend arrays */
        const UINT4 size_ringdown = RINGDOWN_EXTEND_ARRAY;
        REAL8 dt_rngdn = dt;

        /* If ODE does not finish with uniform timestep, then match merger interp timestep */
        if (merger_interp)
            dt_rngdn = dt_merger_interp;

        /* Resize waveform and time array */

        XLALResizeSphHarmPolarTimeSeries(hlm, 0, (size+size_ringdown));
        XLALResizeREAL8Sequence(tdata, 0, (size+size_ringdown));
        for (UINT4 i = size; i < (size+size_ringdown); i++)
        {
            hlm->tdata->data[i] = hlm->tdata->data[i-1] + dt_rngdn;
        }
        size += size_ringdown;

        /* Ringdown attachment */
        eob_wav_ringdown(dyn, hlm);
        strcat(hlm->ampl->name,"_ringdown");
        strcat(hlm->phase->name,"_ringdown");
    } /* End of BBH section */


    /* *****************************************
     * Compute h+, hx
     * *****************************************
     */

    /* Scale to physical units (if necessary) */
    REAL8 amplitude_prefactor = nu * M * LAL_MRSUN_SI / distance;

    /* Azimuthal angle phi follows LAL convention of LIGO-T1800226 for master,
      * where the polarization basis is defined with \f$\Omega=\pi/2\f$ and
      * \f$Z = \sin{\iota}\sin{\Phi}x + \sin{\iota}\cos{\Phi}y + \cos{\iota}z\f$
      */
    const REAL8 phi = LAL_PI_2 - phiRef;
    const REAL8 iota = inclination;

    /* Time to use as the epoch of the returned time series.
     * Time of merger is centered at zero.
     * If ODE did not integrate up to peak then defaults to 0.
     */
    LIGOTimeGPS epoch = LIGOTIMEGPSZERO;
    XLALGPSAdd(&epoch, -t_peak/time_unit_fact);


    /*
     * Interpolate on uniform grid
     */

    /* Auxiliary Structures */
    REAL8Sequence *utime=NULL;

    /* Build uniform grid of width dt and alloc tmp memory */
    const UINT4 size_out = get_uniform_size(tdata->data[size-1], tdata->data[0], deltaT*time_unit_fact);

    /* Waveform */

    /* Create a new mode list to store uniform time data */
    utime = XLALCreateREAL8Sequence(size_out);
    XLAL_CHECK (utime, XLAL_ENOMEM, "Could not allocate memory for time data.\n");

    for (UINT4 i = 0; i < size_out; i++) utime->data[i] = i*deltaT*time_unit_fact;

    /* Interp to uniform grid the multipoles before hpc computation if
     *INTERP_UNIFORM_GRID_HLM scheme is used
     */
    if (interp_uniform_grid == INTERP_UNIFORM_GRID_HLM)
    {
      //#pragma omp parallel
      //{
      /* Interpolate ampl and phase */
        //#pragma omp single private(this_hlm)
        //{
          this_hlm = hlm;
          while (this_hlm)
          {
            //#pragma omp task
            //{
              //fprintf(stderr, "In thread%d/%d\n", omp_get_thread_num(), omp_get_num_threads());
                REAL8TimeSeries *A_ut=NULL;
                REAL8TimeSeries *phi_ut=NULL;
                /* Interpolate mode amplitude and replace */
                A_ut = XLALCreateREAL8TimeSeries(this_hlm->ampl->name, &epoch, 0, deltaT, &(lalStrainUnit), (size_t) size_out);
                // XLAL_CHECK (A_ut, XLAL_ENOMEM, "Could not allocate memory for hlm amplitude data.\n");
                interp_spline(tdata->data, this_hlm->ampl->data->data, size, utime->data,
                              size_out, A_ut->data->data);
                XLALDestroyREAL8TimeSeries(this_hlm->ampl);
                this_hlm->ampl = A_ut;

                /* Interpolate mode phase and replace */
                phi_ut = XLALCreateREAL8TimeSeries(this_hlm->phase->name, &epoch, 0, deltaT, &(lalDimensionlessUnit), (size_t) size_out);
                // XLAL_CHECK (phi_ut, XLAL_ENOMEM, "Could not allocate memory for hlm phase data.\n");
                interp_spline(tdata->data, this_hlm->phase->data->data, size, utime->data,
                              size_out, phi_ut->data->data);
                XLALDestroyREAL8TimeSeries(this_hlm->phase);
                this_hlm->phase = phi_ut;

                A_ut = NULL;
                phi_ut = NULL;
                XLALDestroyREAL8TimeSeries(A_ut);
                XLALDestroyREAL8TimeSeries(phi_ut);
              //}
              this_hlm = this_hlm->next;
              }
          //}
        //}

        /* Replace time sequence */
        size = size_out;
        XLALSphHarmPolarTimeSeriesSetTData(hlm, utime);
    }

    /* Computation of (h+,hx) */

    /* Allocate hplus and hcross */
    *hplus = XLALCreateREAL8TimeSeries( "h_plus", &epoch, 0.0, deltaT, &lalStrainUnit, size_out);
    XLAL_CHECK(*hplus, XLAL_EFAULT, "ERROR allocating hplus.\n");

    *hcross = XLALCreateREAL8TimeSeries( "h_cross", &epoch, 0.0, deltaT, &lalStrainUnit, size_out);
    XLAL_CHECK(*hcross, XLAL_EFAULT, "ERROR allocating hplus.\n");

    memset ((*hplus)->data->data, 0, (*hplus)->data->length*sizeof(*(*hplus)->data->data));
    memset ((*hcross)->data->data, 0, (*hcross)->data->length*sizeof(*(*hcross)->data->data));

    /* Spherical harmonics projection **/
    /* construct hplus and hcross **/
    /* h22 = 1/R * (nu*M)*G/c^2 h_code_output */

    /* Computation of (h+,hx) */
    // NOTE that we set the phase angle for the m mode equal to zero by definition and to be consistent with the PA approximation
    if (interp_uniform_grid == INTERP_UNIFORM_GRID_HLM)
    {
        XLALSimIMRComputePolarisationsPolar((*hplus)->data, (*hcross)->data, hlm, ModeArray, amplitude_prefactor, iota, phi);
    }
    else if (interp_uniform_grid == INTERP_UNIFORM_GRID_HPC)
    {

        /* TEMPORARY ERROR MESSAGE */
        XLAL_ERROR(XLAL_EINVAL, "Only HLM interpolation scheme is available. HPC is under reconstruction! Exiting...\n");

        /* Interpolate (h+,hx) if INTERP_UNIFORM_GRID_HPC scheme is used
         *  This is much faster in the current implementation but may smear
         *  out higher-mode features.
         */

        /* Allocate non-uniform hplus and hcross */
        REAL8Sequence *hp = XLALCreateREAL8Sequence(size);
        XLAL_CHECK(hp, XLAL_EFAULT, "ERROR allocating hplus.\n");

        REAL8Sequence *hc = XLALCreateREAL8Sequence(size);
        XLAL_CHECK(hc, XLAL_EFAULT, "ERROR allocating hcross.\n");

        memset (hp->data, 0, hp->length*sizeof(REAL8));
        memset (hc->data, 0, hc->length*sizeof(REAL8));

        XLALSimIMRComputePolarisationsPolar(hp, hc, hlm, ModeArray, 1.0, iota, phi);

        REAL8 *Apc, *A_interp;
        REAL8 *phipc, *phi_interp;

        Apc = malloc(size * sizeof(REAL8));
        phipc = malloc(size * sizeof(REAL8));
        A_interp = malloc(size_out * sizeof(REAL8));
        phi_interp = malloc(size_out * sizeof(REAL8));

        for (int i=0; i<size; i++) {
            Apc[i] = sqrt( SQ(hp->data[i]) + SQ(hc->data[i]) );
            phipc[i] = - atan2(hc->data[i], hp->data[i]);
        }
        unwrap_proxy(phipc, hlm->next->phase->data->data, size, 0); /* ... but use phi22 as unwrap proxy */


        /* Interp to uniform grid after AP computation */

        /* Interpolate amplitude and phase */
        // FIXME: this is wrong, amplitude interpolation gives unphyiscal modulation when not perfectly face-on!
        interp_spline(tdata->data, Apc, size, utime->data,
                      size_out, A_interp);
        interp_spline(tdata->data, phipc, size, utime->data,
                      size_out, phi_interp);

        // FIXME: this is currently wrong, always gives circular polarization
        for (UINT4 i=0; i<size_out; i++)
        {
            (*hplus)->data->data[i]  = amplitude_prefactor * A_interp[i] * cos(phi_interp[i]);
            (*hcross)->data->data[i] = amplitude_prefactor * A_interp[i] * sin(phi_interp[i]);
        }

        XLALDestroyREAL8Sequence(hp);
        XLALDestroyREAL8Sequence(hc);
        free(Apc);
        free(phipc);
        free(A_interp);
        free(phi_interp);

    } else XLAL_ERROR(XLAL_EINVAL, "Only HLM and HPC interpolation schemes are available.");


    /* Free memory */
    XLALDestroyREAL8Sequence(tdata);
    XLALDestroyREAL8Sequence(utime);

    XLALFreeTEOBDynamics(dyn);
    if (merger_interp)
    {
        XLALFreeTEOBDynamics(dyn_mrg);
        XLALDestroySphHarmPolarTimeSeries(hlm_mrg);
    }

    XLALSphHarmPolarTimeSeriesSetTData(hlm, NULL);
    XLALDestroySphHarmPolarTimeSeries(hlm);

    /* tdata of NQC mode list used to point to hlm->tdata or hlm_mrg->tdata which is now destroyed */
    XLALSphHarmPolarTimeSeriesSetTData(hlm_nqc, NULL);
    XLALDestroySphHarmPolarTimeSeries(hlm_nqc);

    Waveform_lm_t_free(hlm_t);

    XLALDestroyValue(ModeArray);

    return XLAL_SUCCESS;
}

/** @} */
/** @} */
