/*
 * Copyright (C) 2011 E. Ochsner
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
#include <lal/Units.h>
#include <lal/LALConstants.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALAdaptiveRungeKutta4.h>
#include <lal/TimeSeries.h>
#include "check_series_macros.h"

#define UNUSED(expr) do { (void)(expr); } while (0)
/*
#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif
*/

/* use error codes above 1024 to avoid conflicts with GSL */
#define LALSIMINSPIRAL_ST_TEST_ENERGY               1025
#define LALSIMINSPIRAL_ST_TEST_OMEGADOT             1026
#define LALSIMINSPIRAL_ST_TEST_COORDINATE           1027
#define LALSIMINSPIRAL_ST_TEST_OMEGANAN             1028
#define LALSIMINSPIRAL_ST_TEST_FREQBOUND            1029
#define LALSIMINSPIRAL_ST_DERIVATIVE_OMEGANONPOS    1030
#define LALSIMINSPIRAL_ST_TEST_LARGEV               1031

/* (2x) Highest available PN order - UPDATE IF NEW ORDERS ADDED!!*/
#define LAL_MAX_PN_ORDER 8
/* Number of variables used for precessing waveforms */
#define LAL_NUM_ST4_VARIABLES 14
/* absolute and relative tolerance for adaptive Runge-Kutta ODE integrator */
/* 1.e-06 is too large for end of 1.4--1.4 M_sun BNS inspiral */
/* (phase difference at end will be ~10% of GW cycle). */
/* 1.e-12 is used so last data point isn't nan for 6PN tidal, */
/* since larger values probably cause larger step sizes. */
#define LAL_ST4_ABSOLUTE_TOLERANCE 1.e-12
#define LAL_ST4_RELATIVE_TOLERANCE 1.e-12

/**
 * Struct containing all of the non-dynamical coefficients needed
 * to evolve a TaylorTx spinning, precessing binary and produce a waveform.
 * This struct is passed to the static Derivatives and StoppingTest functions.
 */
typedef struct tagXLALSimInspiralSpinTaylorTxCoeffs
{
	REAL8 M; ///< total mass in seconds
	REAL8 Mchirp; ///< chirp mass in seconds
	REAL8 eta; ///< symmetric mass ratio
	REAL8 m1M; ///< m1 / M
	REAL8 m2M; ///< m2 / M
	REAL8 wdotnewt; ///< leading order coefficient of wdot = \f$\dot{\omega}\f$
	REAL8 wdotcoeff[LAL_MAX_PN_ORDER]; ///< coeffs. of PN corrections to wdot
	REAL8 wdotlogcoeff; ///< coefficient of log term in wdot
	REAL8 Ecoeff[LAL_MAX_PN_ORDER]; ///< coeffs. of PN corrections to energy
	REAL8 wdotSO15s1, wdotSO15s2; ///< non-dynamical 1.5PN SO corrections
	REAL8 wdotSS2; ///< non-dynamical 2PN SS correction
	REAL8 wdotQM2S1; ///< non-dynamical S1^2 2PN quadrupole-monopole correction
	REAL8 wdotQM2S1L; ///< non-dynamical (S1.L)^2 2PN quadrupole-monopole correction
	REAL8 wdotQM2S2; ///< non-dynamical S2^2 2PN quadrupole-monopole correction
	REAL8 wdotQM2S2L; ///< non-dynamical (S2.L)^2 2PN quadrupole-monopole correction
	REAL8 wdotSSselfS1; ///< non-dynamical S1^2 2PN self-spin correction
	REAL8 wdotSSselfS1L; ///< non-dynamical (S1.L)^2 self-spin correction
	REAL8 wdotSSselfS2; ///< non-dynamical S2^2 2PN self-spin correction
	REAL8 wdotSSselfS2L; ///< non-dynamical (S2.L)^2 self-spin correction
	REAL8 wdotSO25s1, wdotSO25s2; ///< non-dynamical 2.5PN SO corrections
	REAL8 wdotSO3s1, wdotSO3s2; ///< non-dynamical 3PN SO corrections
	REAL8 wdotSO35s1, wdotSO35s2; ///< non-dynamical 3.5PN SO corrections
	REAL8 ESO15s1, ESO15s2; ///< non-dynamical 1.5PN SO corrections
	REAL8 ESS2; ///< non-dynamical 2PN SS correction
	REAL8 EQM2S1; ///< non-dynamical S1^2 2PN quadrupole-monopole correction
	REAL8 EQM2S1L;///< non-dynamical (S1.L)^2 2PN quadrupole-monopole correction
	REAL8 EQM2S2; ///< non-dynamical S2^2 2PN quadrupole-monopole correction
	REAL8 EQM2S2L;///< non-dynamical (S2.L)^2 2PN quadrupole-monopole correction
	REAL8 ESO25s1, ESO25s2; ///< non-dynamical 2.5PN SO corrections 
	REAL8 ESO35s1, ESO35s2; ///< non-dynamical 3.5PN SO corrections 
	REAL8 LNhatSO15s1, LNhatSO15s2; ///< non-dynamical 1.5PN SO corrections
	REAL8 LNhatSS2; ///< non-dynamical 2PN SS correction 
	REAL8 wdottidal5pn;	///< leading order tidal correction 
	REAL8 wdottidal6pn;	///< next to leading order tidal correction
	REAL8 Etidal5pn; ///< leading order tidal correction to energy
	REAL8 Etidal6pn; ///< next to leading order tidal correction to energy
	REAL8 fStart; ///< starting GW frequency of integration
	REAL8 fEnd; ///< ending GW frequency of integration
	LALSimInspiralSpinOrder spinO; ///< Twice PN order of included spin effects
	LALSimInspiralTidalOrder tideO;///< Twice PN order of included tidal effects
} XLALSimInspiralSpinTaylorTxCoeffs;


/* Declarations of static functions - defined below */
static int XLALSimInspiralSpinTaylorStoppingTest(double t,
	const double values[], double dvalues[], void *mparams);
static int XLALSimInspiralSpinTaylorT4Derivatives(double t, 
	const double values[], double dvalues[], void *mparams);
static int XLALSimInspiralSpinTaylorT4Setup(
    XLALSimInspiralSpinTaylorTxCoeffs *params, REAL8 m1, REAL8 m2,
    REAL8 fStart, REAL8 fEnd, REAL8 lambda1, REAL8 lambda2,
    REAL8 quadparam1, REAL8 quadparam2, LALSimInspiralSpinOrder spinO,
    LALSimInspiralTidalOrder tideO, INT4 phaseO);
static int XLALSimInspiralSpinTaylorT2Derivatives(double t,
	const double values[], double dvalues[], void *mparams);
static int XLALSimInspiralSpinTaylorT2Setup(
    XLALSimInspiralSpinTaylorTxCoeffs *params, REAL8 m1, REAL8 m2,
    REAL8 fStart, REAL8 fEnd, REAL8 lambda1, REAL8 lambda2,
    REAL8 quadparam1, REAL8 quadparam2, LALSimInspiralSpinOrder spinO,
    LALSimInspiralTidalOrder tideO, INT4 phaseO);
static int XLALSimInspiralSpinTaylorDriver(REAL8TimeSeries **hplus,
    REAL8TimeSeries **hcross, REAL8 phiRef, REAL8 v0, REAL8 deltaT,
    REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fRef, REAL8 r,
    REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z,
    REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z,
    REAL8 lambda1, REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2,
    LALSimInspiralSpinOrder spinO, LALSimInspiralTidalOrder tideO,
    int phaseO, int amplitudeO, Approximant approx);



static int XLALSimInspiralSpinTaylorT2Setup(
    XLALSimInspiralSpinTaylorTxCoeffs *params,
    REAL8 m1,                       /**< mass of body 1 (kg) */
    REAL8 m2,                       /**< mass of body 2 (kg) */
    REAL8 fStart,                   /**< Starting GW freq. (Hz) */
    REAL8 fEnd,                     /**< Ending GW freq. (Hz), 0 means integrate forwards as far as possible */
	REAL8 lambda1,                  /**< (tidal deformability of mass 1) / (mass of body 1)^5 (dimensionless) */
	REAL8 lambda2,                  /**< (tidal deformability of mass 2) / (mass of body 2)^5 (dimensionless) */
    REAL8 quadparam1,               /**< phenom. parameter describing induced quad. moment of body 1 (=1 for BHs, ~2-12 for NSs) */
    REAL8 quadparam2,               /**< phenom. parameter describing induced quad. moment of body 2 (=1 for BHs, ~2-12 for NSs) */
	LALSimInspiralSpinOrder spinO,  /**< twice PN order of spin effects */
	LALSimInspiralTidalOrder tideO, /**< twice PN order of tidal effects */
	INT4 phaseO                     /**< twice post-Newtonian order */
    )
{
    REAL8 M, eta, Mchirp, dm, m1m2, m2m1, m1M, m2M;
    /* Zero the coefficients */
    memset(params, 0, sizeof(XLALSimInspiralSpinTaylorTxCoeffs));

    /* Define mass variables and other coefficients */
    m1m2 = m1 / m2;
    m2m1 = m2 / m1;
    m1 *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert m1 from kg to seconds */
    m2 *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert m2 from kg to seconds */
    M = m1 + m2;
    dm = (m1 - m2)/M;
    m1M = m1 / M;
    m2M = m2 / M;
    eta = m1 * m2 / M / M;
    Mchirp = M * pow(eta, 3./5.);
    params->wdotnewt = (96.0/5.0) * eta;
    params->M = M;
    params->Mchirp = Mchirp;
    params->m1M = m1M;
    params->m2M = m2M;
    params->eta = eta;
    params->fStart = fStart;
    params->fEnd = fEnd;
    params->spinO = spinO;
    params->tideO = tideO;

    /* Set coefficients up to PN order phaseO.
     * epnorb is the binary energy and
     * wdotorb is the derivative of the orbital frequency \f$\dot{\omega}\f$.
     * These are just the non-spinning contributions.
     * Spin corrections must be recomputed at every step
     * because the relative orientations of S, L can change
     *
     * The values can be found in Buonanno, Iyer, Ochsner, Pan and Sathyaprakash
     * Phys. Rev. D 80, 084043 (2009) arXiv:0907.0700 (aka \"BIOPS\")
     * Eq. 3.1 for the energy and Eq. 3.6 for \f$\dot{\omega}\f$
     *
     * Note that Eq. 3.6 actually gives dv/dt, but this relates to \f$\omega\f$
     * by \f$d (M \omega)/dt = d (v^3)/dt = 3 v^2 dv/dt\f$
     * so the PN corrections are the same
     * but the leading order is 3 v^2 times Eq. 3.6
     */
    switch( phaseO )
    {
        case -1: /* highest available PN order */
        case 8:
        /* case LAL_PNORDER_THREE_POINT_FIVE: */
        case 7:
            params->wdotcoeff[7] = (LAL_PI/1016064.)
                    * (-15419335. - 12718104.*eta + 4975824.*eta*eta);
            params->Ecoeff[7] = 0.;
        /* case LAL_PNORDER_THREE: */
        case 6:
            params->wdotcoeff[6] = 22.065 + 165.416*eta
                    - 2.20067*eta*eta + 4.93152*eta*eta*eta;
            /*params->wdotcoeff[6] = -10817850546611./93884313600.
                    + 1712.*LAL_GAMMA/105. + 32.*LAL_PI*LAL_PI/3.
                    + (3147553127./12192768. - 451.*LAL_PI*LAL_PI/48.)*eta
                    - 15211.*eta*eta/6912. + 25565.*eta*eta*eta/5184.;*/
            params->wdotlogcoeff = 1712./315.;
            params->Ecoeff[6] = - 675./64. + ( 34445./576.
                    - 205./96. * LAL_PI * LAL_PI ) * eta
                    - (155./96.) *eta * eta - 35./5184. * eta * eta * eta;
        /* case LAL_PNORDER_TWO_POINT_FIVE: */
        case 5:
            params->wdotcoeff[5] = (1./672.) * LAL_PI * (-7729. + 1092.*eta);
            params->Ecoeff[5] = 0.;
        /* case LAL_PNORDER_TWO: */
        case 4:
            params->wdotcoeff[4] = (3058673. + 5472432.*eta
                    + 4353552.*eta*eta)/1016064.;
            params->Ecoeff[4] = (-81. + 57.*eta - eta*eta)/24.;
        /*case LAL_PNORDER_ONE_POINT_FIVE:*/
        case 3:
            params->wdotcoeff[3] = -4. * LAL_PI;
            params->Ecoeff[3] = 0.;
        /*case LAL_PNORDER_ONE:*/
        case 2:
            params->wdotcoeff[2] = (1./336.) * (743. + 924.*eta);
            params->Ecoeff[2] = -(1.0/12.0) * (9.0 + eta);
        /*case LAL_PNORDER_HALF:*/
        case 1:
            params->wdotcoeff[1] = 0.;
            params->Ecoeff[1] = 0.;
        /*case LAL_PNORDER_NEWTONIAN:*/
        case 0:
            params->wdotcoeff[0] = 1.;
            params->Ecoeff[0] = 1.;
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid phase. PN order %s\n",
                    __func__, phaseO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    /* Compute the non-dynamical coefficients of spin corrections
     * to the evolution equations for omega, L, S1 and S2 and binary energy E.
     */
    switch( spinO )
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
            params->wdotSO35s1 = 1796084999./3048192. - 1317319.*eta/3024.
                    + 61.*eta*eta/144. + (dm/m1M) * (-34195607./193536.
                    + 575285.*eta/2688. + 947.*eta*eta/576.);
            params->wdotSO35s2 = 1796084999./3048192. - 1317319.*eta/3024.
                    + 61.*eta*eta/144. - (dm/m2M) * (-34195607./193536.
                    + 575285.*eta/2688. + 947.*eta*eta/576.);
            params->ESO35s1 = 135./4. - 367.*eta/4. + 29.*eta*eta/12.
                    + (dm/m1M) * (-27./4. + 39.*eta - 5.*eta*eta/4.);
            params->ESO35s2 = 135./4. - 367.*eta/4. + 29.*eta*eta/12.
                    - (dm/m2M) * (-27./4. + 39.*eta - 5.*eta*eta/4.);
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
            params->wdotSO3s1 = LAL_PI*(-188./3. + (dm/m1M)*149./6.);
            params->wdotSO3s2 = LAL_PI*(-188./3. - (dm/m2M)*149./6.);
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
            params->wdotSO25s1 = 110869./1008. + 11.*eta/4. + (dm/m1M) * (-8349./224. + 3.*eta/4.);
            params->wdotSO25s2 = 110869./1008. + 11.*eta/4. + (dm/m2M) * (8349./224. - 3.*eta/4.);
            params->ESO25s1 = 11. - 61.*eta/9. + (dm/m1M) * (-3. + 10.*eta/3.);
            params->ESO25s2 = 11. - 61.*eta/9. + (dm/m2M) * (3. - 10.*eta/3.);
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            // 2PN spin-spin terms
            params->LNhatSS2 	= -1.5 / eta;
            params->wdotSS2 	= 1. / 48. / eta;
            params->ESS2 		= 1. / eta;
            // 2PN quadrupole-monopole terms
            params->wdotQM2S1        = quadparam1 * 5./2./m1M/m1M;
            params->wdotQM2S1L       = -quadparam1 * 15./2./m1M/m1M;
            params->wdotQM2S2        = quadparam2 * 5./2./m2M/m2M;
            params->wdotQM2S2L       = -quadparam2 * 15./2./m2M/m2M;
            params->EQM2S1           = quadparam1/2./m1M/m1M;
            params->EQM2S1L          = -quadparam1*3./2./m1M/m1M;
            params->EQM2S2           = quadparam2/2./m2M/m2M;
            params->EQM2S2L          = -quadparam2*3./2./m2M/m2M;
            // 2PN self-spin terms
            params->wdotSSselfS1     = -7./96./m1M/m1M;
            params->wdotSSselfS1L    = 1./96./m1M/m1M;
            params->wdotSSselfS2     = -7./96./m2M/m2M;
            params->wdotSSselfS2L    = 1./96./m2M/m2M;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            // Note: LNHat do not have their signs reversed relative to T4
            // They are precession rather than orbital quantities
            params->LNhatSO15s1 	= 2. + 3./2. * m2m1;
            params->LNhatSO15s2	= 2. + 3./2. * m1m2;
            params->wdotSO15s1 	= ( 113. + 75. * m2m1 ) / 12.;
            params->wdotSO15s2 	= ( 113. + 75. * m1m2 ) / 12.;
            params->ESO15s1 		= 8./3. + 2. * m2m1;
            params->ESO15s2 		= 8./3. + 2. * m1m2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %s\n",
                    __func__, spinO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    /**
     * Compute the coefficients of tidal corrections
     * to the evolution equations for omega and binary energy E.
     * Coefficients found from Eqs. 2.11 and 3.10 of
     * Vines, Flanagan, Hinderer, PRD 83, 084051 (2011).
     */
    switch( tideO )
    {
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL:
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN:
            params->wdottidal6pn = lambda1 * (4421./28. - 12263./28. * m1M
                    + 1893./2. * m1M * m1M - 661 * m1M * m1M * m1M) / (2 * m1M)
                    + lambda2 * (4421./28. - 12263./28. * m2M + 1893./2. * m2M * m2M
                    - 661 * m2M * m2M * m2M) / (2 * m2M);
            params->Etidal6pn = - 11./2. * m2m1
                    * (3. + 2. * m1M + 3. * m1M * m1M) * lambda1
                    - 11./2. * m1m2 * (3. + 2. * m2M + 3. * m2M * m2M) * lambda2;
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
            params->wdottidal5pn = lambda1 * 6. * (1. + 11. * m2M) / m1M
                    + lambda2 * 6. * (1. + 11. * m1M) / m2M;
            params->Etidal5pn = - 9. * m2m1 * lambda1 - 9. * m1m2 * lambda2;
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid tidal PN order %s\n",
                    __func__, tideO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    return XLAL_SUCCESS;
}

static int XLALSimInspiralSpinTaylorT4Setup(
    XLALSimInspiralSpinTaylorTxCoeffs *params,
    REAL8 m1,                       /**< mass of body 1 (kg) */
    REAL8 m2,                       /**< mass of body 2 (kg) */
    REAL8 fStart,                   /**< Starting GW freq. (Hz) */
    REAL8 fEnd,                     /**< Ending GW freq. (Hz), 0 means integrate forwards as far as possible */
	REAL8 lambda1,                  /**< (tidal deformability of mass 1) / (mass of body 1)^5 (dimensionless) */
	REAL8 lambda2,                  /**< (tidal deformability of mass 2) / (mass of body 2)^5 (dimensionless) */
    REAL8 quadparam1,               /**< phenom. parameter describing induced quad. moment of body 1 (=1 for BHs, ~2-12 for NSs) */
    REAL8 quadparam2,               /**< phenom. parameter describing induced quad. moment of body 2 (=1 for BHs, ~2-12 for NSs) */
	LALSimInspiralSpinOrder spinO,  /**< twice PN order of spin effects */
	LALSimInspiralTidalOrder tideO, /**< twice PN order of tidal effects */
	INT4 phaseO                     /**< twice post-Newtonian order */
    )
{
    REAL8 M, eta, Mchirp, dm, m1m2, m2m1, m1M, m2M;
    /* Zero the coefficients */
    memset(params, 0, sizeof(XLALSimInspiralSpinTaylorTxCoeffs));

    /* Define mass variables and other coefficients */
    m1m2 = m1 / m2;
    m2m1 = m2 / m1;
    m1 *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert m1 from kg to seconds */
    m2 *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert m2 from kg to seconds */
    M = m1 + m2;
    dm = (m1 - m2)/M;
    m1M = m1 / M;
    m2M = m2 / M;
    eta = m1 * m2 / M / M;
    Mchirp = M * pow(eta, 3./5.);
    params->wdotnewt = (96.0/5.0) * eta;
    params->M = M;
    params->Mchirp = Mchirp;
    params->m1M = m1M;
    params->m2M = m2M;
    params->eta = eta;
    params->fStart = fStart;
    params->fEnd = fEnd;
    params->spinO = spinO;
    params->tideO = tideO;

    /* Set coefficients up to PN order phaseO.
     * epnorb is the binary energy and
     * wdotorb is the derivative of the orbital frequency \f$\dot{\omega}\f$.
     * These are just the non-spinning contributions.
     * Spin corrections must be recomputed at every step
     * because the relative orientations of S, L can change
     *
     * The values can be found in Buonanno, Iyer, Ochsner, Pan and Sathyaprakash
     * Phys. Rev. D 80, 084043 (2009) arXiv:0907.0700 (aka \"BIOPS\")
     * Eq. 3.1 for the energy and Eq. 3.6 for \f$\dot{\omega}\f$
     *
     * Note that Eq. 3.6 actually gives dv/dt, but this relates to \f$\omega\f$
     * by \f$d (M \omega)/dt = d (v^3)/dt = 3 v^2 dv/dt\f$
     * so the PN corrections are the same 
     * but the leading order is 3 v^2 times Eq. 3.6
     */
    switch( phaseO )
    {
        case -1: /* highest available PN order */
        case 8:
        /* case LAL_PNORDER_THREE_POINT_FIVE: */
        case 7:
            params->wdotcoeff[7] = (LAL_PI/12096.0) 
                    * (-13245.0 + 717350.0*eta + 731960.0*eta*eta);
            params->Ecoeff[7] = 0.;
        /* case LAL_PNORDER_THREE: */
        case 6:
            params->wdotcoeff[6] = 16447322263./139708800. - 1712./105. 
                    * LAL_GAMMA - 56198689./217728. * eta + LAL_PI * LAL_PI 
                    * (16./3. + 451./48. * eta) + 541./896. * eta * eta 
                    - 5605./2592. * eta * eta * eta - 856./105. * log(16.);
            params->wdotlogcoeff = - 1712./315.;
            params->Ecoeff[6] = - 675./64. + ( 34445./576. 
                    - 205./96. * LAL_PI * LAL_PI ) * eta
                    - (155./96.) *eta * eta - 35./5184. * eta * eta * eta;
        /* case LAL_PNORDER_TWO_POINT_FIVE: */
        case 5:
            params->wdotcoeff[5] = -(1./672.) * LAL_PI * (4159. + 15876.*eta);
            params->Ecoeff[5] = 0.;
        /* case LAL_PNORDER_TWO: */
        case 4:
            params->wdotcoeff[4] = (34103. + 122949.*eta 
                    + 59472.*eta*eta)/18144.;
            params->Ecoeff[4] = (-81. + 57.*eta - eta*eta)/24.;
        /*case LAL_PNORDER_ONE_POINT_FIVE:*/
        case 3:
            params->wdotcoeff[3] = 4. * LAL_PI;
            params->Ecoeff[3] = 0.;
        /*case LAL_PNORDER_ONE:*/
        case 2:
            params->wdotcoeff[2] = -(1./336.) * (743. + 924.*eta);
            params->Ecoeff[2] = -(1.0/12.0) * (9.0 + eta);
        /*case LAL_PNORDER_HALF:*/
        case 1:
            params->wdotcoeff[1] = 0.;
            params->Ecoeff[1] = 0.;
        /*case LAL_PNORDER_NEWTONIAN:*/
        case 0:
            params->wdotcoeff[0] = 1.;
            params->Ecoeff[0] = 1.;
            break;
        default: 
            XLALPrintError("XLAL Error - %s: Invalid phase. PN order %s\n", 
                    __func__, phaseO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    /* Compute the non-dynamical coefficients of spin corrections
     * to the evolution equations for omega, L, S1 and S2 and binary energy E.
     */
    switch( spinO )
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
            params->wdotSO35s1 = -4323559./18144. + 436705.*eta/672.
                    - 5575.*eta*eta/27. + (dm/m1M) * (1195759./18144.
                    - 257023.*eta/1008. + 2903.*eta*eta/32.);
            params->wdotSO35s2 = -4323559./18144. + 436705.*eta/672.
                    - 5575.*eta*eta/27. - (dm/m2M) * (1195759./18144.
                    - 257023.*eta/1008. + 2903.*eta*eta/32.);
            params->ESO35s1 = 135./4. - 367.*eta/4. + 29.*eta*eta/12.
                    + (dm/m1M) * (-27./4. + 39.*eta - 5.*eta*eta/4.);
            params->ESO35s2 = 135./4. - 367.*eta/4. + 29.*eta*eta/12.
                    - (dm/m2M) * (-27./4. + 39.*eta - 5.*eta*eta/4.);
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
            params->wdotSO3s1 = LAL_PI*(-188./3. + (dm/m1M)*151./6.);
            params->wdotSO3s2 = LAL_PI*(-188./3. - (dm/m2M)*151./6.);
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
            params->wdotSO25s1 = -5861./144. + 1001.*eta/12. 
                    + (dm/m1M) * (809./84. - 281.*eta/8.);
            params->wdotSO25s2 = -5861./144. + 1001.*eta/12.
                    + (dm/m2M) * (-809./84. + 281.*eta/8.);
            params->ESO25s1 = 11. - 61.*eta/9. + (dm/m1M) * (-3. + 10.*eta/3.);
            params->ESO25s2 = 11. - 61.*eta/9. + (dm/m2M) * (3. - 10.*eta/3.);
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            // 2PN spin-spin terms
            params->LNhatSS2 	= -1.5 / eta;
            params->wdotSS2 		= - 1. / 48. / eta;
            params->ESS2 		= 1. / eta;
            // 2PN quadrupole-monopole terms
            params->wdotQM2S1 	= -quadparam1 * 5./2./m1M/m1M;
            params->wdotQM2S1L 	= quadparam1 * 15./2./m1M/m1M;
            params->wdotQM2S2 	= -quadparam2 * 5./2./m2M/m2M;
            params->wdotQM2S2L 	= quadparam2 * 15./2./m2M/m2M;
            params->EQM2S1 		= quadparam1/2./m1M/m1M;
            params->EQM2S1L 		= -quadparam1*3./2./m1M/m1M;
            params->EQM2S2 		= quadparam2/2./m2M/m2M;
            params->EQM2S2L 		= -quadparam2*3./2./m2M/m2M;
            // 2PN self-spin terms
            params->wdotSSselfS1     = 7./96./m1M/m1M;
            params->wdotSSselfS1L    = -1./96./m1M/m1M;
            params->wdotSSselfS2     = 7./96./m2M/m2M;
            params->wdotSSselfS2L    = -1./96./m2M/m2M;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            params->LNhatSO15s1 	= 2. + 3./2. * m2m1;
            params->LNhatSO15s2	= 2. + 3./2. * m1m2;
            params->wdotSO15s1 	= - ( 113. + 75. * m2m1 ) / 12.;
            params->wdotSO15s2 	= - ( 113. + 75. * m1m2 ) / 12.;
            params->ESO15s1 		= 8./3. + 2. * m2m1;
            params->ESO15s2 		= 8./3. + 2. * m1m2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %s\n",
                    __func__, spinO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }
	
    /* Compute the coefficients of tidal corrections
     * to the evolution equations for omega and binary energy E.
     * Coefficients found from Eqs. 2.11 and 3.10 of 
     * Vines, Flanagan, Hinderer, PRD 83, 084051 (2011).
     */
    switch( tideO )
    {
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL:
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN:
            params->wdottidal6pn = lambda1 * m1M*m1M*m1M*m1M*m1M
                    * (4421./28. - 12263./28. * m1M + 1893./2. * m1M * m1M
                    - 661 * m1M * m1M * m1M) / (2 * m1M) + lambda2
                    * m2M*m2M*m2M*m2M*m2M * (4421./28. - 12263./28. * m2M
                    + 1893./2. * m2M * m2M - 661 * m2M * m2M * m2M) / (2 * m2M);
            params->Etidal6pn = - 11./2. * m2m1 * (3. + 2. * m1M
                    + 3. * m1M * m1M) * lambda1 * m1M*m1M*m1M*m1M*m1M
                    - 11./2. * m1m2 * (3. + 2. * m2M + 3. * m2M * m2M)
                    * lambda2 * m2M*m2M*m2M*m2M*m2M;
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
            params->wdottidal5pn = lambda1 * 6. * (1. + 11. * m2M)
                    * m1M*m1M*m1M*m1M + lambda2
                    * 6. * (1. + 11. * m1M) * m2M*m2M*m2M*m2M;
            params->Etidal5pn = - 9. * m2m1 * lambda1 * m1M*m1M*m1M*m1M*m1M
                    - 9. * m1m2 * lambda2 * m2M*m2M*m2M*m2M*m2M;
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid tidal PN order %s\n",
                    __func__, tideO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    return XLAL_SUCCESS;
}

/**
 * This function evolves the orbital equations for a precessing binary using 
 * the \"TaylorT1/T2/T4\" approximant for solving the orbital dynamics
 * (see arXiv:0907.0700 for a review of the various PN approximants).
 *
 * It returns time series of the \"orbital velocity\", orbital phase, 
 * and components for both individual spin vectors, the \"Newtonian\"
 * orbital angular momentum (which defines the instantaneous plane)
 * and "E1", a basis vector in the instantaneous orbital plane.
 * Note that LNhat and E1 completely specify the instantaneous orbital plane.
 * It also returns the time and phase of the final time step
 *
 * For input, the function takes the two masses, the initial orbital phase, 
 * Values of S1, S2, LNhat, E1 vectors at starting time,
 * the desired time step size, the starting GW frequency, 
 * and PN order at which to evolve the phase,
 * 
 * NOTE: All vectors are given in the so-called "radiation frame", 
 * where the direction of propagation is the z-axis, the principal "+" 
 * polarization axis is the x-axis, and the y-axis is given by the RH rule.
 * You must give the initial values in this frame, and the time series of the
 * vector components will also be returned in this frame
 */
int XLALSimInspiralSpinTaylorPNEvolveOrbit(
	REAL8TimeSeries **V,            /**< post-Newtonian parameter [returned]*/
	REAL8TimeSeries **Phi,          /**< orbital phase            [returned]*/
	REAL8TimeSeries **S1x,	        /**< Spin1 vector x component [returned]*/
	REAL8TimeSeries **S1y,	        /**< "    "    "  y component [returned]*/
	REAL8TimeSeries **S1z,	        /**< "    "    "  z component [returned]*/
	REAL8TimeSeries **S2x,	        /**< Spin2 vector x component [returned]*/
	REAL8TimeSeries **S2y,	        /**< "    "    "  y component [returned]*/
	REAL8TimeSeries **S2z,	        /**< "    "    "  z component [returned]*/
	REAL8TimeSeries **LNhatx,       /**< unit orbital ang. mom. x [returned]*/
	REAL8TimeSeries **LNhaty,       /**< "    "    "  y component [returned]*/
	REAL8TimeSeries **LNhatz,       /**< "    "    "  z component [returned]*/
	REAL8TimeSeries **E1x,	        /**< orb. plane basis vector x[returned]*/
	REAL8TimeSeries **E1y,	        /**< "    "    "  y component [returned]*/
	REAL8TimeSeries **E1z,	        /**< "    "    "  z component [returned]*/
	REAL8 deltaT,          	        /**< sampling interval (s) */
	REAL8 m1,              	        /**< mass of companion 1 (kg) */
	REAL8 m2,              	        /**< mass of companion 2 (kg) */
	REAL8 fStart,                   /**< starting GW frequency */
	REAL8 fEnd,                     /**< ending GW frequency, fEnd=0 means integrate as far forward as possible */
	REAL8 s1x,                      /**< initial value of S1x */
	REAL8 s1y,                      /**< initial value of S1y */
	REAL8 s1z,                      /**< initial value of S1z */
	REAL8 s2x,                      /**< initial value of S2x */
	REAL8 s2y,                      /**< initial value of S2y */
	REAL8 s2z,                      /**< initial value of S2z */
	REAL8 lnhatx,                   /**< initial value of LNhatx */
	REAL8 lnhaty,                   /**< initial value of LNhaty */
	REAL8 lnhatz,                   /**< initial value of LNhatz */
	REAL8 e1x,                      /**< initial value of E1x */
	REAL8 e1y,                      /**< initial value of E1y */
	REAL8 e1z,                      /**< initial value of E1z */
	REAL8 lambda1,                  /**< (tidal deformability of mass 1) / (mass of body 1)^5 (dimensionless) */
	REAL8 lambda2,                  /**< (tidal deformability of mass 2) / (mass of body 2)^5 (dimensionless) */
	REAL8 quadparam1,               /**< phenom. parameter describing induced quad. moment of body 1 (=1 for BHs, ~2-12 for NSs) */
	REAL8 quadparam2,               /**< phenom. parameter describing induced quad. moment of body 2 (=1 for BHs, ~2-12 for NSs) */
	LALSimInspiralSpinOrder spinO,  /**< twice PN order of spin effects */
	LALSimInspiralTidalOrder tideO, /**< twice PN order of tidal effects */
	INT4 phaseO,                    /**< twice post-Newtonian order */
    Approximant approx              /**< PN approximant (SpinTaylorT1/T2/T4) */
	)
{
    INT4 intreturn;
    void * params;
    ark4GSLIntegrator *integrator = NULL;     /* GSL integrator object */
    REAL8 yinit[LAL_NUM_ST4_VARIABLES];       /* initial values of parameters */
    REAL8Array *yout;	 /* time series of variables returned from integrator */
    /* intermediate variables */
    UINT4 i, cutlen, len;
    int sgn, offset;
    REAL8 norm, dtStart, dtEnd, lengths, wEnd, m1sec, m2sec, Msec, Mcsec;
    LIGOTimeGPS tStart = LIGOTIMEGPSZERO;

    /* Check start and end frequencies are positive */
    if( fStart <= 0. )
    {
        XLALPrintError("XLAL Error - %s: fStart = %f must be > 0.\n", 
                __func__, fStart );
        XLAL_ERROR(XLAL_EINVAL);
    }
    if( fEnd < 0. ) /* fEnd = 0 allowed as special case */
    {
        XLALPrintError("XLAL Error - %s: fEnd = %f must be >= 0.\n", 
                __func__, fEnd );
        XLAL_ERROR(XLAL_EINVAL);
    }

    /* Set sign of time step according to direction of integration */
    if( fEnd < fStart && fEnd != 0. )
        sgn = -1;
    else
        sgn = 1;

    /* Check start and end frequencies are positive */
    if( fStart <= 0. )
    {
        XLALPrintError("XLAL Error - %s: fStart = %f must be > 0.\n", 
                __func__, fStart );
        XLAL_ERROR(XLAL_EINVAL);
    }
    if( fEnd < 0. ) /* fEnd = 0 allowed as special case */
    {
        XLALPrintError("XLAL Error - %s: fEnd = %f must be >= 0.\n", 
                __func__, fEnd );
        XLAL_ERROR(XLAL_EINVAL);
    }

    // Fill params struct with values of constant coefficients of the model
    if( approx == SpinTaylorT4 )
    {
        XLALSimInspiralSpinTaylorTxCoeffs paramsT4;
        XLALSimInspiralSpinTaylorT4Setup(&paramsT4, m1, m2, fStart, fEnd,
                lambda1, lambda2, quadparam1, quadparam2, spinO, tideO, phaseO);
        params = (void *) &paramsT4;
    }
    else if( approx == SpinTaylorT2 )
    {
        XLALSimInspiralSpinTaylorTxCoeffs paramsT2;
        XLALSimInspiralSpinTaylorT2Setup(&paramsT2, m1, m2, fStart, fEnd,
                lambda1, lambda2, quadparam1, quadparam2, spinO, tideO, phaseO);
        params = (void *) &paramsT2;
    }
    else if( approx == SpinTaylorT1 )
    {
        //XLALSimInspiralSpinTaylorT1Coeffs paramsT1;
        //XLALSimInspiralSpinTaylorT1Setup(&paramsT1, m1, m2, fStart, fEnd,
        //        lambda1, lambda2, quadparam1, quadparam2, spinO, tideO, phaseO);
        //params = (void *) &paramsT1;
        XLALPrintError("XLAL Error - %s: SpinTaylorT1 not implemented yet!\n",
                __func__);
        XLAL_ERROR(XLAL_EINVAL);
    }
    else
    {
        XLALPrintError("XLAL Error - %s: Approximant must be one of SpinTaylorT1, SpinTaylorT2, SpinTaylorT4, but %i provided\n",
                __func__, approx);
        XLAL_ERROR(XLAL_EINVAL);

    }
    m1sec = m1 * LAL_G_SI / pow(LAL_C_SI, 3.0);
    m2sec = m2 * LAL_G_SI / pow(LAL_C_SI, 3.0);
    Msec = m1sec + m2sec;
    Mcsec = Msec * pow( m1sec*m2sec/Msec/Msec, 0.6);
	   
    /* Estimate length of waveform using Newtonian t(f) formula */
    /* Time from freq. = fStart to infinity */
    dtStart = (5.0/256.0) * pow(LAL_PI,-8.0/3.0) 
            * pow(Mcsec * fStart,-5.0/3.0) / fStart;
    /* Time from freq. = fEnd to infinity. Set to zero if fEnd=0 */
    dtEnd = (fEnd == 0. ? 0. : (5.0/256.0) * pow(LAL_PI,-8.0/3.0) 
            * pow(Mcsec * fEnd,-5.0/3.0) / fEnd);
    /* Time in sec from fStart to fEnd. Note it can be positive or negative */
    lengths = dtStart - dtEnd;

    /* Put initial values into a single array for the integrator */
    yinit[0] = 0.; /* without loss of generality, set initial orbital phase=0 */
    yinit[1] = LAL_PI * Msec * fStart;  /* \hat{omega} = (pi M f) */
    /* LNh(x,y,z) */
    yinit[2] = lnhatx;
    yinit[3] = lnhaty;
    yinit[4] = lnhatz;
    /* S1(x,y,z) */
    norm = m1sec * m1sec / Msec / Msec;
    yinit[5] = norm * s1x;
    yinit[6] = norm * s1y;
    yinit[7] = norm * s1z;
    /* S2(x,y,z) */
    norm = m2sec * m2sec / Msec / Msec;
    yinit[8] = norm * s2x;
    yinit[9] = norm * s2y;
    yinit[10]= norm * s2z;
    /* E1(x,y,z) */
    yinit[11] = e1x;
    yinit[12] = e1y;
    yinit[13] = e1z;

    /* initialize the integrator */
    if( approx == SpinTaylorT4 )
        integrator = XLALAdaptiveRungeKutta4Init(LAL_NUM_ST4_VARIABLES,
                XLALSimInspiralSpinTaylorT4Derivatives,
                XLALSimInspiralSpinTaylorStoppingTest,
                LAL_ST4_ABSOLUTE_TOLERANCE, LAL_ST4_RELATIVE_TOLERANCE);
    else if( approx == SpinTaylorT2 )
        integrator = XLALAdaptiveRungeKutta4Init(LAL_NUM_ST4_VARIABLES,
                XLALSimInspiralSpinTaylorT2Derivatives,
                XLALSimInspiralSpinTaylorStoppingTest,
                LAL_ST4_ABSOLUTE_TOLERANCE, LAL_ST4_RELATIVE_TOLERANCE);
    else if( approx == SpinTaylorT1 )
    {
        //integrator = XLALAdaptiveRungeKutta4Init(LAL_NUM_ST4_VARIABLES,
        //        XLALSimInspiralSpinTaylorT1Derivatives,
        //        XLALSimInspiralSpinTaylorStoppingTest,
        //        LAL_ST4_ABSOLUTE_TOLERANCE, LAL_ST4_RELATIVE_TOLERANCE);
        XLALPrintError("XLAL Error - %s: SpinTaylorT1 not implemented yet!\n",
                __func__);
        XLAL_ERROR(XLAL_EINVAL);
    }
    else
    {
        XLALPrintError("XLAL Error - %s: Approximant must be one of SpinTaylorT1, SpinTaylorT2, SpinTaylorT4, but %i provided\n",
                __func__, approx);
        XLAL_ERROR(XLAL_EINVAL);

    }
    if( !integrator )
    {
        XLALPrintError("XLAL Error - %s: Cannot allocate integrator\n", 
                __func__);
        XLAL_ERROR(XLAL_EFUNC);
    }

    /* stop the integration only when the test is true */
    integrator->stopontestonly = 1;

    /* run the integration; note: time is measured in \hat{t} = t / M */
    //len = XLALAdaptiveRungeKutta4Hermite(integrator, (void *) &params, yinit,
    len = XLALAdaptiveRungeKutta4Hermite(integrator, params, yinit,
            0.0, lengths/Msec, sgn*deltaT/Msec, &yout);

    intreturn = integrator->returncode;
    XLALAdaptiveRungeKutta4Free(integrator);

    if (!len) 
    {
        XLALPrintError("XLAL Error - %s: integration failed with errorcode %d.\n", __func__, intreturn);
        XLAL_ERROR(XLAL_EFUNC);
    }

    /* Print warning about abnormal termination */
    if (intreturn != 0 && intreturn != LALSIMINSPIRAL_ST_TEST_ENERGY
            && intreturn != LALSIMINSPIRAL_ST_TEST_OMEGADOT
            && intreturn != LALSIMINSPIRAL_ST_TEST_FREQBOUND)
    {
        XLALPrintWarning("XLAL Warning - %s: integration terminated with code %d.\n Waveform parameters were m1 = %e, m2 = %e, s1 = (%e,%e,%e), s2 = (%e,%e,%e), inc = %e.\n", __func__, intreturn, m1 * pow(LAL_C_SI, 3.0) / LAL_G_SI / LAL_MSUN_SI, m2 * pow(LAL_C_SI, 3.0) / LAL_G_SI / LAL_MSUN_SI, s1x, s1y, s1z, s2x, s2y, s2z, acos(lnhatz));
    }

    /* 
     * If ending frequency was non-zero, we may have overshot somewhat.
     * The integrator takes one adaptive stride past fEnd, 
     * but this may include several smaller interpolation steps.
     * Therefore, 'cutlen' will be the index of the first interpolated step
     * to cross fEnd and 'len' is the full length returned from the integrator.
     * If fEnd == 0, we integrated as far as possible and 'cutlen' = 'len'.
     */
    cutlen = len;
    if( fEnd != 0. && fEnd < fStart )
    {
        wEnd = LAL_PI * Msec * fEnd;/* Ending dimensionless freq. \hat{omega} */
        /* Integrator returns \hat{omega} in array 'yout'
           in range data[2*len] to data[2*len+(len-1)]. 
           Start at end and find where we cross wEnd */
        while( yout->data[2*len+cutlen-1] < wEnd )
            cutlen--;
        if( cutlen < len )
            cutlen++; /* while loop exits on wrong side of fEnd, so increment */
    }
    else if( fEnd > fStart )
    {
        wEnd = LAL_PI * Msec * fEnd;/* Ending dimensionless freq. \hat{omega} */
        /* Integrator returns \hat{omega} in array 'yout'
           in range data[2*len] to data[2*len+(len-1)]. 
           Start at end and find where we cross wEnd */
        while( yout->data[2*len+cutlen-1] > wEnd )
            cutlen--;
        if( cutlen < len )
            cutlen++; /* while loop exits on wrong side of fEnd, so increment */
    }

    /* Adjust tStart so last sample is at time=0 */
    XLALGPSAdd(&tStart, -1.0*(cutlen-1)*deltaT);

    /* allocate memory for output vectors */
    *V = XLALCreateREAL8TimeSeries( "PN_EXPANSION_PARAMETER", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, cutlen); 
    *Phi = XLALCreateREAL8TimeSeries( "ORBITAL_PHASE", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, cutlen); 
    *S1x = XLALCreateREAL8TimeSeries( "SPIN1_X_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, cutlen); 
    *S1y = XLALCreateREAL8TimeSeries( "SPIN1_Y_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, cutlen); 
    *S1z = XLALCreateREAL8TimeSeries( "SPIN1_Z_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, cutlen); 
    *S2x = XLALCreateREAL8TimeSeries( "SPIN2_X_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, cutlen); 
    *S2y = XLALCreateREAL8TimeSeries( "SPIN2_Y_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, cutlen); 
    *S2z = XLALCreateREAL8TimeSeries( "SPIN2_Z_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, cutlen); 
    *LNhatx = XLALCreateREAL8TimeSeries( "LNHAT_X_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, cutlen); 
    *LNhaty = XLALCreateREAL8TimeSeries( "LNHAT_Y_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, cutlen); 
    *LNhatz = XLALCreateREAL8TimeSeries( "LNHAT_Z_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, cutlen); 
    *E1x = XLALCreateREAL8TimeSeries( "E1_BASIS_X_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, cutlen); 
    *E1y = XLALCreateREAL8TimeSeries( "E1_BASIS_Y_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, cutlen); 
    *E1z = XLALCreateREAL8TimeSeries( "E1_BASIS_Z_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, cutlen); 
    if ( !V || !Phi || !S1x || !S1y || !S1z || !S2x || !S2y || !S2z 
            || !LNhatx || !LNhaty || !LNhatz || !E1x || !E1y || !E1z )
    {
        XLALDestroyREAL8Array(yout);
        XLAL_ERROR(XLAL_EFUNC);
    }

    /* If we integrated backwards, offset & sgn will reverse order of samples */
    if( fEnd < fStart && fEnd != 0. )
        offset = cutlen-1;
    else
        offset = 0;

    /* Copy dynamical variables from yout array to output time series.
     * Note the first 'len' members of yout are the time steps. 
     * Also, the for loop only goes to 'cutlen', in case we overshot fEnd.
     * If we integrated backwards, we copy backwards from 'cutlen'.
     */
    for( i = 0; i < cutlen; i++ )
    {	
        int j = sgn*i+offset;
        (*Phi)->data->data[j] 		= yout->data[len+i];
        (*V)->data->data[j] 		= cbrt(yout->data[2*len+i]);
        (*LNhatx)->data->data[j] 	= yout->data[3*len+i];
        (*LNhaty)->data->data[j] 	= yout->data[4*len+i];
        (*LNhatz)->data->data[j] 	= yout->data[5*len+i];
        (*S1x)->data->data[j] 		= yout->data[6*len+i];
        (*S1y)->data->data[j] 		= yout->data[7*len+i];
        (*S1z)->data->data[j] 		= yout->data[8*len+i];
        (*S2x)->data->data[j] 		= yout->data[9*len+i];
        (*S2y)->data->data[j] 		= yout->data[10*len+i];
        (*S2z)->data->data[j] 		= yout->data[11*len+i];
        (*E1x)->data->data[j] 		= yout->data[12*len+i];
        (*E1y)->data->data[j] 		= yout->data[13*len+i];
        (*E1z)->data->data[j] 		= yout->data[14*len+i];
    }

    XLALDestroyREAL8Array(yout);

    return XLAL_SUCCESS;
}

/**
 * Internal function called by the integration routine.
 * Stops the integration if
 * 		1) The energy decreases with increasing orbital frequency
 *		2) The orbital frequency begins decreasing
 *		3) The orbital frequency becomes infinite 
 *		4) The orbital frequency has gone outside the requested bounds
 *      5) The PN parameter v/c becomes >= 1
 * SpinTaylorT4 and SpinTaylorT2 both use this same stopping test
 */
static int XLALSimInspiralSpinTaylorStoppingTest(
	double t, 
	const double values[],
	double dvalues[], 
	void *mparams
	)
{
    REAL8 omega, v, test, omegaStart, omegaEnd;
    REAL8 LNhx, LNhy, LNhz, S1x, S1y, S1z, S2x, S2y, S2z;
    REAL8 LNdotS1, LNdotS2, S1dotS2, S1sq, S2sq;
    XLALSimInspiralSpinTaylorTxCoeffs *params 
            = (XLALSimInspiralSpinTaylorTxCoeffs*) mparams;
    /* Spin-corrections to energy (including dynamical terms) */
    REAL8 Espin15 = 0., Espin2 = 0., Espin25 = 0., Espin35 = 0.;

    UNUSED(t);

    omega = values[1];
    v = pow(omega,1./3.);
    LNhx = values[2]; LNhy  = values[3]; LNhz = values[4] ;
    S1x  = values[5]; S1y   = values[6]; S1z  = values[7] ;
    S2x  = values[8]; S2y   = values[9]; S2z  = values[10];
    LNdotS1 = (LNhx*S1x + LNhy*S1y + LNhz*S1z);
    LNdotS2 = (LNhx*S2x + LNhy*S2y + LNhz*S2z);

    /* omega = PI G M f_GW / c^3
     * Note params->M is really G M /c^3 (i.e. M is in seconds) */
    omegaStart = LAL_PI * params->M * params->fStart;
    omegaEnd = LAL_PI * params->M * params->fEnd;

    switch( params->spinO )
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
            // Compute 3.5PN SO correction to energy
            // See Eq. 3.15 of arXiv:1303.7412
            // Note that S_l/M^2 = (m1/M)^2 chi1 + (m2/M)^2 chi2
            // and Sigma_l/M^2 = (m2/M) chi2 - (m1/M) chi1
            Espin35 += params->ESO35s1 * LNdotS1 + params->ESO35s2 * LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
            // Compute 2.5PN SO correction to energy
            // See Eq. 7.9 of gr-qc/0605140v4
            // Note that S_l/M^2 = (m1/M)^2 chi1 + (m2/M)^2 chi2
            // and Sigma_l/M^2 = (m2/M) chi2 - (m1/M) chi1
            Espin25 += params->ESO25s1 * LNdotS1 + params->ESO25s2 * LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            // Compute S1-S2 spin-spin term
            S1dotS2 = (S1x*S2x + S1y*S2y + S1z*S2z);
            Espin2 += params->ESS2  * (S1dotS2 - 3. * LNdotS1 * LNdotS2);
            // Compute 2PN quadrupole-monopole correction to energy
            // See last line of Eq. 6 of astro-ph/0504538
            // or 2nd and 3rd lines of Eq. (C4) in arXiv:0810.5336v3
            S1sq = (S1x*S1x + S1y*S1y + S1z*S1z);
            S2sq = (S2x*S2x + S2y*S2y + S2z*S2z);
            Espin2 += params->EQM2S1 * S1sq
                    + params->EQM2S2 * S2sq
                    + params->EQM2S1L * LNdotS1 * LNdotS1
                    + params->EQM2S2L * LNdotS2 * LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            // Compute 1.5PN SO correction to energy
            Espin15 += params->ESO15s1 * LNdotS1 + params->ESO15s2 * LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %s\n",
                    __func__, params->spinO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    /* We are testing if the orbital energy increases with \f$\omega\f$.
     * We should be losing energy to GW flux, so if E increases 
     * we stop integration because the dynamics are becoming unphysical. 
     * 'test' is the PN expansion of \f$dE/d\omega\f$ without the prefactor, 
     * i.e. \f$dE/d\omega = dE/dv * dv/d\omega = - (M^2*eta/6) * test\f$
     * Therefore, the energy is increasing with \f$\omega\f$ iff. test < 0.
     */
    test = 2. + v * v * ( 4. * params->Ecoeff[2] 
            + v * ( 5. * (params->Ecoeff[3] + Espin15) 
            + v * ( 6. * (params->Ecoeff[4] + Espin2)
            + v * ( 7. * (params->Ecoeff[5] + Espin25)
            + v * ( 8. *  params->Ecoeff[6]
            + v * ( 9. * (params->Ecoeff[7] + Espin35)
			+ v * v * v * ( 12. * params->Etidal5pn
			+ v * v * ( 14. * params->Etidal6pn ) ) ) ) ) ) ) );
    if( abs(omegaEnd) > LAL_REAL4_EPS && omegaEnd > omegaStart && omega > omegaEnd) /* freq. above bound */
        return LALSIMINSPIRAL_ST_TEST_FREQBOUND;
    else if( abs(omegaEnd) > LAL_REAL4_EPS && omegaEnd < omegaStart && omega < omegaEnd) /* freq. below bound */
        return LALSIMINSPIRAL_ST_TEST_FREQBOUND;
    else if (test < 0.0) /* energy test fails! */
        return LALSIMINSPIRAL_ST_TEST_ENERGY;
    else if (dvalues[1] < 0.0) /* omegadot < 0! */
        return LALSIMINSPIRAL_ST_TEST_OMEGADOT;
    else if isnan(omega) /* omega is nan! */
        return LALSIMINSPIRAL_ST_TEST_OMEGANAN;
    else if (v >= 1.) // v/c >= 1!
        return LALSIMINSPIRAL_ST_TEST_LARGEV;
    else /* Step successful, continue integrating */
        return GSL_SUCCESS;
}

/**
 * Internal function called by the integration routine.
 * Given the values of all the dynamical variables 'values' at time 't',
 * This function computes their derivatives 'dvalues'
 * so the ODE integrator can take a step
 * All non-dynamical quantities (masses, etc.) are passed in \"mparams\"
 *
 * The derivatives for \f$\omega\f$, L_N, S1, S2 can be found 
 * as Eqs. (1), (8) - (10) of gr-qc/0405090
 * The derivative of E1 is Eq. (15)-(16) of gr-qc/0310034
 */
static int XLALSimInspiralSpinTaylorT4Derivatives(
	double t, 
	const double values[],
	double dvalues[], 
	void *mparams
	) 
{
    /* coordinates and derivatives */
    REAL8 LNhx, LNhy, LNhz, S1x, S1y, S1z, S2x, S2y, S2z, E1x, E1y, E1z;
    REAL8 omega, ds, domega, dLNhx, dLNhy, dLNhz;
    REAL8 dS1x, dS1y, dS1z, dS2x, dS2y, dS2z, dE1x, dE1y, dE1z;

    /* auxiliary variables */
    REAL8 v, v2, v3, v4, v5, v7, v11, omega2, omega2by2;
    REAL8 LNdotS1, LNdotS2, threeLNdotS1, threeLNdotS2, S1dotS2, S1sq, S2sq;
    REAL8 v5etaLNhatSO15s1, v5etaLNhatSO15s2;
    REAL8 OmegaLx, OmegaLy, OmegaLz, OmegaLdotLN;
    REAL8 OmegaEx, OmegaEy, OmegaEz, OmegaSx, OmegaSy, OmegaSz;
    REAL8 wspin15 = 0., wspin2 = 0., wspin25 = 0., wspin3 = 0., wspin35 = 0.;

    XLALSimInspiralSpinTaylorTxCoeffs *params 
            = (XLALSimInspiralSpinTaylorTxCoeffs*) mparams;

    UNUSED(t);

    /* copy variables */
    // UNUSED!!: s    = values[0] ;
    omega   	= values[1] ;
    LNhx = values[2] ; LNhy    	= values[3] ; LNhz 	= values[4] ;
    S1x  = values[5] ; S1y     	= values[6] ; S1z 	= values[7] ;
    S2x  = values[8] ; S2y     	= values[9] ; S2z 	= values[10];
    E1x  = values[11]; E1y    	= values[12]; E1z 	= values[13];

    if (omega <= 0.0) /* orbital frequency must be positive! */
    {
        return LALSIMINSPIRAL_ST_DERIVATIVE_OMEGANONPOS;
    }

    v = cbrt(omega);
    v2  = v * v; v3 = v2 * v; v4 = v3 * v; 
    v5 = v * v4; v7 = v4 * v3; v11 = v7 * v4;

    LNdotS1 = (LNhx*S1x + LNhy*S1y + LNhz*S1z);
    LNdotS2 = (LNhx*S2x + LNhy*S2y + LNhz*S2z);
    S1dotS2 = (S1x*S2x  + S1y*S2y  + S1z*S2z );

    /*
     * domega
     * 
     * Note we are actually computing \f$d \hat{\omega} / d \hat{t}\f$
     * where \f$\hat{\omega} = M \omega\f$ and \f$\hat{t} = t / M\f$
     * Therefore \f$domega = M^2 * d\omega / dt\f$
     *
     * See Eqs. (1)-(7) of gr-qc/0405090 But note that our spin variables 
     * are scaled by component masses relative to that paper.
     * i.e. \f$S_i = (m_i/M)^2 * \hat{S_i}\f$
     *
     * non-spinning coefficients of \f$\dot{\omega}\f$ (params->wdotcoeff[i])
     * should have been set before this function was called
     */

    switch( params->spinO )
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
            // Compute 3.5PN SO correction to domega/dt
            // See Eq. 3.16 of arXiv:1303.7412
            // Note that S_l/M^2 = (m1/M)^2 chi1 + (m2/M)^2 chi2
            // and Sigma_l/M^2 = (m2/M) chi2 - (m1/M) chi1
            wspin35 = params->wdotSO35s1*LNdotS1 + params->wdotSO35s2*LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
            // Compute 3PN SO correction to domega/dt
            // See Eq. 13 of arXiv:1210.0764 and/or Eq. 3.16 of arXiv:1303.7412
            // Note that S_l/M^2 = (m1/M)^2 chi1 + (m2/M)^2 chi2
            // and Sigma_l/M^2 = (m2/M) chi2 - (m1/M) chi1
            wspin3 = params->wdotSO3s1 * LNdotS1 + params->wdotSO3s2 * LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
            // Compute 2.5PN SO correction to domega/dt
            // See Eq. 8.3 of gr-qc/0605140v4
            // Note that S_l/M^2 = (m1/M)^2 chi1 + (m2/M)^2 chi2
            // and Sigma_l/M^2 = (m2/M) chi2 - (m1/M) chi1
            wspin25 = params->wdotSO25s1*LNdotS1 + params->wdotSO25s2*LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            // Compute S1-S2 spin-spin term
            S1dotS2 = (S1x*S2x + S1y*S2y + S1z*S2z);
            wspin2 = params->wdotSS2 * (247.*S1dotS2 - 721.*LNdotS1 * LNdotS2);
            // Compute 2PN QM and self-spin corrections to domega/dt
            // This is equivalent to Eqs. 9c + 9d of astro-ph/0504538
            S1sq = (S1x*S1x + S1y*S1y + S1z*S1z);
            S2sq = (S2x*S2x + S2y*S2y + S2z*S2z);
            wspin2 += params->wdotQM2S1 * S1sq
                    + params->wdotQM2S2 * S2sq
                    + params->wdotQM2S1L * LNdotS1 * LNdotS1
                    + params->wdotQM2S2L * LNdotS2 * LNdotS2
                    + params->wdotSSselfS1 * S1sq
                    + params->wdotSSselfS2 * S2sq
                    + params->wdotSSselfS1L * LNdotS1 * LNdotS1
                    + params->wdotSSselfS2L * LNdotS2 * LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            // Compute 1.5PN SO correction to domega/dt
            wspin15 = params->wdotSO15s1*LNdotS1 + params->wdotSO15s2*LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %s\n",
                    __func__, params->spinO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    domega  = params->wdotnewt * v11 * ( params->wdotcoeff[0]
            + v * ( params->wdotcoeff[1]
            + v * ( params->wdotcoeff[2]
            + v * ( params->wdotcoeff[3] + wspin15
            + v * ( params->wdotcoeff[4] + wspin2
            + v * ( params->wdotcoeff[5] + wspin25
            + v * ( params->wdotcoeff[6] + wspin3
                    + params->wdotlogcoeff * log(omega)
            + v * ( params->wdotcoeff[7] + wspin35
            + v3 * ( params->wdottidal5pn
            + v2 * ( params->wdottidal6pn ) ) ) ) ) ) ) ) ) );

    /*
     * dLN
     * 
     * \f$d \hat{L_N}/d \hat{t} = M * d\hat{L_N} / dt = \Omega_L x \hat{L_N}\f$
     * This is Eq. (10) of gr-qc/0405090 ( times M b/c we use \f$\hat{t}\f$)
     */
    omega2 = omega * omega;
    /* \Omega_L vector */
    OmegaLx = omega2 * (params->LNhatSO15s1 * S1x + params->LNhatSO15s2 * S2x)
            + v7 * params->LNhatSS2 * (LNdotS2 * S1x + LNdotS1 * S2x);
    OmegaLy = omega2 * (params->LNhatSO15s1 * S1y + params->LNhatSO15s2 * S2y)
            + v7 * params->LNhatSS2 * (LNdotS2 * S1y + LNdotS1 * S2y);
    OmegaLz = omega2 * (params->LNhatSO15s1 * S1z + params->LNhatSO15s2 * S2z)
            + v7 * params->LNhatSS2 * (LNdotS2 * S1z + LNdotS1 * S2z);

    /* Take cross product of \Omega_L with \hat{L_N} */
    dLNhx = (-OmegaLz*LNhy + OmegaLy*LNhz);
    dLNhy = (-OmegaLx*LNhz + OmegaLz*LNhx);
    dLNhz = (-OmegaLy*LNhx + OmegaLx*LNhy);

    /*
     * dE1
     * 
     * d E_1 / d \hat{t} = M * d E_1 / dt
     * Computed from \Omega_L and \hat{L_N} with Eq. (15)-(16) of gr-qc/0310034
     */
    OmegaLdotLN = OmegaLx * LNhx + OmegaLy * LNhy + OmegaLz * LNhz;
    /* \Omega_E vector */
    OmegaEx = OmegaLx - OmegaLdotLN * LNhx;
    OmegaEy = OmegaLy - OmegaLdotLN * LNhy;
    OmegaEz = OmegaLz - OmegaLdotLN * LNhz;

    /* Take cross product of \Omega_E with E_1 */
    dE1x = (-OmegaEz*E1y + OmegaEy*E1z);
    dE1y = (-OmegaEx*E1z + OmegaEz*E1x);
    dE1z = (-OmegaEy*E1x + OmegaEx*E1y);

    /*
     * dS1
     * 
     * d S_1 / d \hat{t} = M * d S_1 / dt = \Omega_{S1} x S_1
     * This is Eq. (8) of gr-qc/0405090.
     * However, that paper uses spin variables which are M^2 times our spins
     */
    /* \Omega_{S1} vector */
    omega2by2 = omega2 * 0.5;
    threeLNdotS2 = 3. * LNdotS2;
    v5etaLNhatSO15s1 = v5 * params->eta * params->LNhatSO15s1;
    OmegaSx = v5etaLNhatSO15s1 * LNhx
            + omega2by2 * (S2x - threeLNdotS2 * LNhx);
    OmegaSy = v5etaLNhatSO15s1 * LNhy
            + omega2by2 * (S2y - threeLNdotS2 * LNhy);
    OmegaSz = v5etaLNhatSO15s1 * LNhz
            + omega2by2 * (S2z - threeLNdotS2 * LNhz);

    /* Take cross product of \Omega_{S1} with S_1 */
    dS1x = (-OmegaSz*S1y + OmegaSy*S1z);
    dS1y = (-OmegaSx*S1z + OmegaSz*S1x);
    dS1z = (-OmegaSy*S1x + OmegaSx*S1y);

    /*
     * dS2
     * 
     * d S_2 / d \hat{t} = M * d S_2 / dt = \Omega_{S2} x S_2
     * This is Eq. (9) of gr-qc/0405090.
     * However, that paper uses spin variables which are M^2 times our spins
     */
    /* \Omega_{S2} vector */
    threeLNdotS1 = 3. * LNdotS1;
    v5etaLNhatSO15s2 = v5 * params->eta * params->LNhatSO15s2;
    OmegaSx = v5etaLNhatSO15s2 * LNhx
            + omega2by2 * (S1x - threeLNdotS1 * LNhx);
    OmegaSy = v5etaLNhatSO15s2 * LNhy
            + omega2by2 * (S1y - threeLNdotS1 * LNhy);
    OmegaSz = v5etaLNhatSO15s2 * LNhz
            + omega2by2 * (S1z - threeLNdotS1 * LNhz);

    /* Take cross product of \Omega_{S2} with S_2 */
    dS2x = (-OmegaSz*S2y + OmegaSy*S2z);
    dS2y = (-OmegaSx*S2z + OmegaSz*S2x);
    dS2z = (-OmegaSy*S2x + OmegaSx*S2y);

    /* dphi = d \phi / d \hat{t} = M d \phi /dt = M \omega = \hat{\omega} */
    ds = omega;

    dvalues[0]    = ds   ; dvalues[1]     = domega;
    dvalues[2]    = dLNhx; dvalues[3]     = dLNhy ; dvalues[4]    = dLNhz;
    dvalues[5]    = dS1x ; dvalues[6]     = dS1y  ; dvalues[7]    = dS1z ;
    dvalues[8]    = dS2x ; dvalues[9]     = dS2y  ; dvalues[10]   = dS2z ;
    dvalues[11]   = dE1x ; dvalues[12]    = dE1y  ; dvalues[13]   = dE1z ;

    return GSL_SUCCESS;
}

/**
 * Internal function called by the integration routine.
 * Given the values of all the dynamical variables 'values' at time 't',
 * This function computes their derivatives 'dvalues'
 * so the ODE integrator can take a step
 * All non-dynamical quantities (masses, etc.) are passed in \"mparams\"
 *
 * FIXME: Change references below
 * The derivatives for \f$\omega\f$, L_N, S1, S2 can be found
 * as Eqs. (1), (8) - (10) of gr-qc/0405090
 * The derivative of E1 is Eq. (15)-(16) of gr-qc/0310034
 */
static int XLALSimInspiralSpinTaylorT2Derivatives(
	double t,
	const double values[],
	double dvalues[],
	void *mparams
	)
{
    /* coordinates and derivatives */
    REAL8 LNhx, LNhy, LNhz, S1x, S1y, S1z, S2x, S2y, S2z, E1x, E1y, E1z;
    REAL8 omega, ds, domega, dLNhx, dLNhy, dLNhz;
    REAL8 dS1x, dS1y, dS1z, dS2x, dS2y, dS2z, dE1x, dE1y, dE1z;

    /* auxiliary variables */
    REAL8 v, v2, v3, v4, v5, v7, v11, omega2, omega2by2;
    REAL8 LNdotS1, LNdotS2, threeLNdotS1, threeLNdotS2, S1dotS2, S1sq, S2sq;
    REAL8 v5etaLNhatSO15s1, v5etaLNhatSO15s2;
    REAL8 OmegaLx, OmegaLy, OmegaLz, OmegaLdotLN;
    REAL8 OmegaEx, OmegaEy, OmegaEz, OmegaSx, OmegaSy, OmegaSz;
    REAL8 wspin15 = 0., wspin2 = 0., wspin25 = 0., wspin3 = 0., wspin35 = 0.;

    XLALSimInspiralSpinTaylorTxCoeffs *params
            = (XLALSimInspiralSpinTaylorTxCoeffs*) mparams;

    UNUSED(t);

    /* copy variables */
    // UNUSED!!: s    = values[0] ;
    omega   	= values[1] ;
    LNhx = values[2] ; LNhy    	= values[3] ; LNhz 	= values[4] ;
    S1x  = values[5] ; S1y     	= values[6] ; S1z 	= values[7] ;
    S2x  = values[8] ; S2y     	= values[9] ; S2z 	= values[10];
    E1x  = values[11]; E1y    	= values[12]; E1z 	= values[13];

    if (omega <= 0.0) /* orbital frequency must be positive! */
    {
        return LALSIMINSPIRAL_ST_DERIVATIVE_OMEGANONPOS;
    }

    v = cbrt(omega);
    v2  = v * v; v3 = v2 * v; v4 = v3 * v;
    v5 = v * v4; v7 = v4 * v3; v11 = v7 * v4;

    LNdotS1 = (LNhx*S1x + LNhy*S1y + LNhz*S1z);
    LNdotS2 = (LNhx*S2x + LNhy*S2y + LNhz*S2z);
    S1dotS2 = (S1x*S2x  + S1y*S2y  + S1z*S2z );

    /**
     * domega
     *
     * Note we are actually computing \f$d \hat{\omega} / d \hat{t}\f$
     * where \f$\hat{\omega} = M \omega\f$ and \f$\hat{t} = t / M\f$
     * Therefore \f$domega = M^2 * d\omega / dt\f$
     *
     * See Eqs. (1)-(7) of gr-qc/0405090 But note that our spin variables
     * are scaled by component masses relative to that paper.
     * i.e. \f$S_i = (m_i/M)^2 * \hat{S_i}\f$
     *
     * non-spinning coefficients of \f$\dot{\omega}\f$ (params->wdotcoeff[i])
     * should have been set before this function was called
     */

    switch( params->spinO )
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
            // Compute 3.5PN SO correction to domega/dt
            // See Eq. 3.16 of arXiv:1303.7412
            // Note that S_l/M^2 = (m1/M)^2 chi1 + (m2/M)^2 chi2
            // and Sigma_l/M^2 = (m2/M) chi2 - (m1/M) chi1
            wspin35 = params->wdotSO35s1*LNdotS1 + params->wdotSO35s2*LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
            // Compute 3PN SO correction to domega/dt
            // See Eq. 13 of arXiv:1210.0764 and/or Eq. 3.16 of arXiv:1303.7412
            // Note that S_l/M^2 = (m1/M)^2 chi1 + (m2/M)^2 chi2
            // and Sigma_l/M^2 = (m2/M) chi2 - (m1/M) chi1
            wspin3 = params->wdotSO3s1 * LNdotS1 + params->wdotSO3s2 * LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
            // Compute 2.5PN SO correction to domega/dt
            // See Eq. 8.3 of gr-qc/0605140v4
            // Note that S_l/M^2 = (m1/M)^2 chi1 + (m2/M)^2 chi2
            // and Sigma_l/M^2 = (m2/M) chi2 - (m1/M) chi1
            wspin25 = params->wdotSO25s1*LNdotS1 + params->wdotSO25s2*LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            // Compute S1-S2 spin-spin term
            S1dotS2 = (S1x*S2x + S1y*S2y + S1z*S2z);
            wspin2 = params->wdotSS2 * (247.*S1dotS2 - 721.*LNdotS1 * LNdotS2);
            // Compute 2PN QM and self-spin corrections to domega/dt
            // See last line of Eq. 5.17 of arXiv:0812.4413
            // Also note this is equivalent to Eqs. 9c + 9d of astro-ph/0504538
            S1sq = (S1x*S1x + S1y*S1y + S1z*S1z);
            S2sq = (S2x*S2x + S2y*S2y + S2z*S2z);
            wspin2 += params->wdotQM2S1 * S1sq
                    + params->wdotQM2S2 * S2sq
                    + params->wdotQM2S1L * LNdotS1 * LNdotS1
                    + params->wdotQM2S2L * LNdotS2 * LNdotS2
                    + params->wdotSSselfS1 * S1sq
                    + params->wdotSSselfS2 * S2sq
                    + params->wdotSSselfS1L * LNdotS1 * LNdotS1
                    + params->wdotSSselfS2L * LNdotS2 * LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            // Compute 1.5PN SO correction to domega/dt
            wspin15 = params->wdotSO15s1*LNdotS1 + params->wdotSO15s2*LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %s\n",
                    __func__, params->spinO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    domega  = params->wdotnewt * v11 / ( params->wdotcoeff[0]
            + v * ( params->wdotcoeff[1]
            + v * ( params->wdotcoeff[2]
            + v * ( params->wdotcoeff[3] + wspin15
            + v * ( params->wdotcoeff[4] + wspin2
            + v * ( params->wdotcoeff[5] + wspin25
            + v * ( params->wdotcoeff[6] + wspin3
                    + params->wdotlogcoeff * log(omega)
            + v * ( params->wdotcoeff[7] + wspin35
            + v3 * ( params->wdottidal5pn
            + v2 * ( params->wdottidal6pn ) ) ) ) ) ) ) ) ) );

    /**
     * dLN
     *
     * \f$d \hat{L_N}/d \hat{t} = M * d\hat{L_N} / dt = \Omega_L x \hat{L_N}\f$
     * This is Eq. (10) of gr-qc/0405090 ( times M b/c we use \f$\hat{t}\f$)
     */
    omega2 = omega * omega;
    /* \Omega_L vector */
    OmegaLx = omega2 * (params->LNhatSO15s1 * S1x + params->LNhatSO15s2 * S2x)
            + v7 * params->LNhatSS2 * (LNdotS2 * S1x + LNdotS1 * S2x);
    OmegaLy = omega2 * (params->LNhatSO15s1 * S1y + params->LNhatSO15s2 * S2y)
            + v7 * params->LNhatSS2 * (LNdotS2 * S1y + LNdotS1 * S2y);
    OmegaLz = omega2 * (params->LNhatSO15s1 * S1z + params->LNhatSO15s2 * S2z)
            + v7 * params->LNhatSS2 * (LNdotS2 * S1z + LNdotS1 * S2z);

    /* Take cross product of \Omega_L with \hat{L_N} */
    dLNhx = (-OmegaLz*LNhy + OmegaLy*LNhz);
    dLNhy = (-OmegaLx*LNhz + OmegaLz*LNhx);
    dLNhz = (-OmegaLy*LNhx + OmegaLx*LNhy);

    /**
     * dE1
     *
     * d E_1 / d \hat{t} = M * d E_1 / dt
     * Computed from \Omega_L and \hat{L_N} with Eq. (15)-(16) of gr-qc/0310034
     */
    OmegaLdotLN = OmegaLx * LNhx + OmegaLy * LNhy + OmegaLz * LNhz;
    /* \Omega_E vector */
    OmegaEx = OmegaLx - OmegaLdotLN * LNhx;
    OmegaEy = OmegaLy - OmegaLdotLN * LNhy;
    OmegaEz = OmegaLz - OmegaLdotLN * LNhz;

    /* Take cross product of \Omega_E with E_1 */
    dE1x = (-OmegaEz*E1y + OmegaEy*E1z);
    dE1y = (-OmegaEx*E1z + OmegaEz*E1x);
    dE1z = (-OmegaEy*E1x + OmegaEx*E1y);

    /**
     * dS1
     *
     * d S_1 / d \hat{t} = M * d S_1 / dt = \Omega_{S1} x S_1
     * This is Eq. (8) of gr-qc/0405090.
     * However, that paper uses spin variables which are M^2 times our spins
     */
    /* \Omega_{S1} vector */
    omega2by2 = omega2 * 0.5;
    threeLNdotS2 = 3. * LNdotS2;
    v5etaLNhatSO15s1 = v5 * params->eta * params->LNhatSO15s1;
    OmegaSx = v5etaLNhatSO15s1 * LNhx
            + omega2by2 * (S2x - threeLNdotS2 * LNhx);
    OmegaSy = v5etaLNhatSO15s1 * LNhy
            + omega2by2 * (S2y - threeLNdotS2 * LNhy);
    OmegaSz = v5etaLNhatSO15s1 * LNhz
            + omega2by2 * (S2z - threeLNdotS2 * LNhz);

    /* Take cross product of \Omega_{S1} with S_1 */
    dS1x = (-OmegaSz*S1y + OmegaSy*S1z);
    dS1y = (-OmegaSx*S1z + OmegaSz*S1x);
    dS1z = (-OmegaSy*S1x + OmegaSx*S1y);

    /**
     * dS2
     *
     * d S_2 / d \hat{t} = M * d S_2 / dt = \Omega_{S2} x S_2
     * This is Eq. (9) of gr-qc/0405090.
     * However, that paper uses spin variables which are M^2 times our spins
     */
    /* \Omega_{S2} vector */
    threeLNdotS1 = 3. * LNdotS1;
    v5etaLNhatSO15s2 = v5 * params->eta * params->LNhatSO15s2;
    OmegaSx = v5etaLNhatSO15s2 * LNhx
            + omega2by2 * (S1x - threeLNdotS1 * LNhx);
    OmegaSy = v5etaLNhatSO15s2 * LNhy
            + omega2by2 * (S1y - threeLNdotS1 * LNhy);
    OmegaSz = v5etaLNhatSO15s2 * LNhz
            + omega2by2 * (S1z - threeLNdotS1 * LNhz);

    /* Take cross product of \Omega_{S2} with S_2 */
    dS2x = (-OmegaSz*S2y + OmegaSy*S2z);
    dS2y = (-OmegaSx*S2z + OmegaSz*S2x);
    dS2z = (-OmegaSy*S2x + OmegaSx*S2y);

    /* dphi = d \phi / d \hat{t} = M d \phi /dt = M \omega = \hat{\omega} */
    ds = omega;

    dvalues[0]    = ds   ; dvalues[1]     = domega;
    dvalues[2]    = dLNhx; dvalues[3]     = dLNhy ; dvalues[4]    = dLNhz;
    dvalues[5]    = dS1x ; dvalues[6]     = dS1y  ; dvalues[7]    = dS1z ;
    dvalues[8]    = dS2x ; dvalues[9]     = dS2y  ; dvalues[10]   = dS2z ;
    dvalues[11]   = dE1x ; dvalues[12]    = dE1y  ; dvalues[13]   = dE1z ;

    return GSL_SUCCESS;
}

/* Appends the start and end time series together, skipping the redundant first
 * sample of end.  Frees end before returning a pointer to the result, which is
 * the resized start series.  */
static REAL8TimeSeries *appendTSandFree(REAL8TimeSeries *start, 
        REAL8TimeSeries *end) {
    unsigned int origlen = start->data->length;
    start = XLALResizeREAL8TimeSeries(start, 0, 
            start->data->length + end->data->length - 1);
    
    memcpy(start->data->data + origlen, end->data->data+1, 
            (end->data->length-1)*sizeof(REAL8));

    XLALGPSAdd(&(start->epoch), -end->deltaT*(end->data->length - 1));

    XLALDestroyREAL8TimeSeries(end);

    return start;        
}


/**
 * Driver routine to compute a precessing post-Newtonian inspiral waveform
 * with phasing computed from energy balance using the so-called \"T4\" method.
 *
 * This routine allows the user to specify different pN orders
 * for the phasing and amplitude of the waveform.
 *
 * The reference frequency fRef is used as follows:
 * 1) if fRef = 0: The initial values of s1, s2, lnhat and e1 will be the
 *    values at frequency fStart. The orbital phase of the last sample is set
 *    to phiRef (i.e. phiRef is the "coalescence phase", roughly speaking).
 *    THIS IS THE DEFAULT BEHAVIOR CONSISTENT WITH OTHER APPROXIMANTS
 * 
 * 2) If fRef = fStart: The initial values of s1, s2, lnhat and e1 will be the
 *    values at frequency fStart. phiRef is used to set the orbital phase
 *    of the first sample at fStart.
 *
 * 3) If fRef > fStart: The initial values of s1, s2, lnhat and e1 will be the
 *    values at frequency fRef. phiRef is used to set the orbital phase at fRef.
 *    The code will integrate forwards and backwards from fRef and stitch the
 *    two together to create a complete waveform. This allows one to specify
 *    the orientation of the binary in-band (or at any arbitrary point).
 *    Otherwise, the user can only directly control the initial orientation.
 *
 * 4) fRef < 0 or fRef >= Schwarz. ISCO are forbidden and the code will abort.
 */
int XLALSimInspiralSpinTaylorT4(
	REAL8TimeSeries **hplus,        /**< +-polarization waveform */
	REAL8TimeSeries **hcross,       /**< x-polarization waveform */
	REAL8 phiRef,                   /**< orbital phase at reference pt. */
	REAL8 v0,                       /**< tail gauge term (default = 1) */
	REAL8 deltaT,                   /**< sampling interval (s) */
	REAL8 m1,                       /**< mass of companion 1 (kg) */
	REAL8 m2,                       /**< mass of companion 2 (kg) */
	REAL8 fStart,                   /**< start GW frequency (Hz) */
	REAL8 fRef,                     /**< reference GW frequency (Hz) */
	REAL8 r,                        /**< distance of source (m) */
	REAL8 s1x,                      /**< initial value of S1x */
	REAL8 s1y,                      /**< initial value of S1y */
	REAL8 s1z,                      /**< initial value of S1z */
	REAL8 s2x,                      /**< initial value of S2x */
	REAL8 s2y,                      /**< initial value of S2y */
	REAL8 s2z,                      /**< initial value of S2z */
	REAL8 lnhatx,                   /**< initial value of LNhatx */
	REAL8 lnhaty,                   /**< initial value of LNhaty */
	REAL8 lnhatz,                   /**< initial value of LNhatz */
	REAL8 e1x,                      /**< initial value of E1x */
	REAL8 e1y,                      /**< initial value of E1y */
	REAL8 e1z,                      /**< initial value of E1z */
	REAL8 lambda1,                  /**< (tidal deformability of mass 1) / (mass of body 1)^5 (dimensionless) */
	REAL8 lambda2,                  /**< (tidal deformability of mass 2) / (mass of body 2)^5 (dimensionless) */
	REAL8 quadparam1,               /**< phenom. parameter describing induced quad. moment of body 1 (=1 for BHs, ~2-12 for NSs) */
	REAL8 quadparam2,               /**< phenom. parameter describing induced quad. moment of body 2 (=1 for BHs, ~2-12 for NSs) */
	LALSimInspiralSpinOrder spinO,  /**< twice PN order of spin effects */
	LALSimInspiralTidalOrder tideO, /**< twice PN order of tidal effects */
	int phaseO,                     /**< twice PN phase order */
	int amplitudeO                  /**< twice PN amplitude order */
	)
{
    Approximant approx = SpinTaylorT4;
    int n = XLALSimInspiralSpinTaylorDriver(hplus, hcross, phiRef, v0, deltaT,
            m1, m2, fStart, fRef, r, s1x, s1y, s1z, s2x, s2y, s2z,
            lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2,
            quadparam1, quadparam2, spinO, tideO, phaseO, amplitudeO, approx);

    return n;
}

/**
 * Driver routine to compute a precessing post-Newtonian inspiral waveform
 * with phasing computed from energy balance using the so-called \"T2\" method.
 *
 * This routine allows the user to specify different pN orders
 * for the phasing, amplitude, spin and tidal effects of the waveform.
 * 
 * The reference frequency fRef is used as follows:
 * 1) if fRef = 0: The initial values of s1, s2, lnhat and e1 will be the
 *    values at frequency fStart. The orbital phase of the last sample is set
 *    to phiRef (i.e. phiRef is the "coalescence phase", roughly speaking).
 *    THIS IS THE DEFAULT BEHAVIOR CONSISTENT WITH OTHER APPROXIMANTS
 * 
 * 2) If fRef = fStart: The initial values of s1, s2, lnhat and e1 will be the 
 *    values at frequency fStart. phiRef is used to set the orbital phase
 *    of the first sample at fStart.
 * 
 * 3) If fRef > fStart: The initial values of s1, s2, lnhat and e1 will be the
 *    values at frequency fRef. phiRef is used to set the orbital phase at fRef.
 *    The code will integrate forwards and backwards from fRef and stitch the
 *    two together to create a complete waveform. This allows one to specify
 *    the orientation of the binary in-band (or at any arbitrary point).
 *    Otherwise, the user can only directly control the initial orientation.
 *
 * 4) fRef < 0 or fRef >= Schwarz. ISCO are forbidden and the code will abort.
 */
int XLALSimInspiralSpinTaylorT2(
	REAL8TimeSeries **hplus,        /**< +-polarization waveform */
	REAL8TimeSeries **hcross,       /**< x-polarization waveform */
	REAL8 phiRef,                   /**< orbital phase at reference pt. */
	REAL8 v0,                       /**< tail gauge term (default = 1) */
	REAL8 deltaT,                   /**< sampling interval (s) */
	REAL8 m1,                       /**< mass of companion 1 (kg) */
	REAL8 m2,                       /**< mass of companion 2 (kg) */
	REAL8 fStart,                   /**< start GW frequency (Hz) */
	REAL8 fRef,                     /**< reference GW frequency (Hz) */
	REAL8 r,                        /**< distance of source (m) */
	REAL8 s1x,                      /**< initial value of S1x */
	REAL8 s1y,                      /**< initial value of S1y */
	REAL8 s1z,                      /**< initial value of S1z */
	REAL8 s2x,                      /**< initial value of S2x */
	REAL8 s2y,                      /**< initial value of S2y */
	REAL8 s2z,                      /**< initial value of S2z */
	REAL8 lnhatx,                   /**< initial value of LNhatx */
	REAL8 lnhaty,                   /**< initial value of LNhaty */
	REAL8 lnhatz,                   /**< initial value of LNhatz */
	REAL8 e1x,                      /**< initial value of E1x */
	REAL8 e1y,                      /**< initial value of E1y */
	REAL8 e1z,                      /**< initial value of E1z */
	REAL8 lambda1,                  /**< (tidal deformability of mass 1) / (mass of body 1)^5 (dimensionless) */
	REAL8 lambda2,                  /**< (tidal deformability of mass 2) / (mass of body 2)^5 (dimensionless) */
	REAL8 quadparam1,               /**< phenom. parameter describing induced quad. moment of body 1 (=1 for BHs, ~2-12 for NSs) */
	REAL8 quadparam2,               /**< phenom. parameter describing induced quad. moment of body 2 (=1 for BHs, ~2-12 for NSs) */
	LALSimInspiralSpinOrder spinO,  /**< twice PN order of spin effects */
	LALSimInspiralTidalOrder tideO, /**< twice PN order of tidal effects */
	int phaseO,                     /**< twice PN phase order */
	int amplitudeO                  /**< twice PN amplitude order */
	)
{
    Approximant approx = SpinTaylorT2;
    int n = XLALSimInspiralSpinTaylorDriver(hplus, hcross, phiRef, v0, deltaT,
            m1, m2, fStart, fRef, r, s1x, s1y, s1z, s2x, s2y, s2z,
            lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2,
            quadparam1, quadparam2, spinO, tideO, phaseO, amplitudeO, approx);

    return n;
}

/**
 * Internal driver function to generate any of SpinTaylorT1/T2/T4
 */
static int XLALSimInspiralSpinTaylorDriver(
	REAL8TimeSeries **hplus,        /**< +-polarization waveform */
	REAL8TimeSeries **hcross,       /**< x-polarization waveform */
	REAL8 phiRef,                   /**< orbital phase at reference pt. */
	REAL8 v0,                       /**< tail gauge term (default = 1) */
	REAL8 deltaT,                   /**< sampling interval (s) */
	REAL8 m1,                       /**< mass of companion 1 (kg) */
	REAL8 m2,                       /**< mass of companion 2 (kg) */
	REAL8 fStart,                   /**< start GW frequency (Hz) */
	REAL8 fRef,                     /**< reference GW frequency (Hz) */
	REAL8 r,                        /**< distance of source (m) */
	REAL8 s1x,                      /**< initial value of S1x */
	REAL8 s1y,                      /**< initial value of S1y */
	REAL8 s1z,                      /**< initial value of S1z */
	REAL8 s2x,                      /**< initial value of S2x */
	REAL8 s2y,                      /**< initial value of S2y */
	REAL8 s2z,                      /**< initial value of S2z */
	REAL8 lnhatx,                   /**< initial value of LNhatx */
	REAL8 lnhaty,                   /**< initial value of LNhaty */
	REAL8 lnhatz,                   /**< initial value of LNhatz */
	REAL8 e1x,                      /**< initial value of E1x */
	REAL8 e1y,                      /**< initial value of E1y */
	REAL8 e1z,                      /**< initial value of E1z */
	REAL8 lambda1,                  /**< (tidal deformability of mass 1) / (mass of body 1)^5 (dimensionless) */
	REAL8 lambda2,                  /**< (tidal deformability of mass 2) / (mass of body 2)^5 (dimensionless) */
	REAL8 quadparam1,               /**< phenom. parameter describing induced quad. moment of body 1 (=1 for BHs, ~2-12 for NSs) */
	REAL8 quadparam2,               /**< phenom. parameter describing induced quad. moment of body 2 (=1 for BHs, ~2-12 for NSs) */
	LALSimInspiralSpinOrder spinO,  /**< twice PN order of spin effects */
	LALSimInspiralTidalOrder tideO, /**< twice PN order of tidal effects */
	int phaseO,                     /**< twice PN phase order */
	int amplitudeO,                 /**< twice PN amplitude order */
    Approximant approx              /**< PN approximant (SpinTaylorT1/T2/T4) */
	)
{
    REAL8TimeSeries *V, *Phi, *S1x, *S1y, *S1z, *S2x, *S2y, *S2z;
    REAL8TimeSeries *LNhatx, *LNhaty, *LNhatz, *E1x, *E1y, *E1z;
    int status, n;
    unsigned int i;
    REAL8 fS, fE, phiShift;
    /* The Schwarzschild ISCO frequency - for sanity checking fRef */
    REAL8 fISCO = pow(LAL_C_SI,3) / (pow(6.,3./2.)*LAL_PI*(m1+m2)*LAL_G_SI);

    /* Sanity check fRef value */
    if( fRef < 0. )
    {
        XLALPrintError("XLAL Error - %s: fRef = %f must be >= 0\n", 
                __func__, fRef);
        XLAL_ERROR(XLAL_EINVAL);
    }
    if( fRef != 0. && fRef < fStart )
    {
        XLALPrintError("XLAL Error - %s: fRef = %f must be > fStart = %f\n", 
                __func__, fRef, fStart);
        XLAL_ERROR(XLAL_EINVAL);
    }
    if( fRef >= fISCO )
    {
        XLALPrintError("XLAL Error - %s: fRef = %f must be < Schwar. ISCO=%f\n",
                __func__, fRef, fISCO);
        XLAL_ERROR(XLAL_EINVAL);
    }

    /* if fRef=0, just integrate from start to end. Let phiRef=phiC */
    if( fRef < LAL_REAL4_EPS )
    {
        fS = fStart;
        fE = 0.;
        /* Evolve the dynamical variables */
        n = XLALSimInspiralSpinTaylorPNEvolveOrbit(&V, &Phi,
                &S1x, &S1y, &S1z, &S2x, &S2y, &S2z, 
                &LNhatx, &LNhaty, &LNhatz, &E1x, &E1y, &E1z, 
                deltaT, m1, m2, fS, fE, s1x, s1y, s1z, s2x, s2y, s2z, 
                lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2,
                quadparam1, quadparam2, spinO, tideO, phaseO, approx);
        if( n < 0 )
            XLAL_ERROR(XLAL_EFUNC);

        /* Apply phase shift so orbital phase ends with desired value */
        phiShift = phiRef - Phi->data->data[Phi->data->length-1];
        for( i=0; i < Phi->data->length; i++)
        {
            Phi->data->data[i] += phiShift;
        }
    }
    /* if fRef=fStart, just integrate from start to end. Let phiRef=phiStart */
    else if( abs(fRef - fStart) < LAL_REAL4_EPS )
    {
        fS = fStart;
        fE = 0.;
        /* Evolve the dynamical variables */
        n = XLALSimInspiralSpinTaylorPNEvolveOrbit(&V, &Phi,
                &S1x, &S1y, &S1z, &S2x, &S2y, &S2z, 
                &LNhatx, &LNhaty, &LNhatz, &E1x, &E1y, &E1z, 
                deltaT, m1, m2, fS, fE, s1x, s1y, s1z, s2x, s2y, s2z, 
                lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2,
                quadparam1, quadparam2, spinO, tideO, phaseO, approx);
        if( n < 0 )
            XLAL_ERROR(XLAL_EFUNC);

        /* Apply phase shift so orbital phase starts with desired value */
        phiShift = phiRef - Phi->data->data[0];
        for( i=0; i < Phi->data->length; i++)
        {
            Phi->data->data[i] += phiShift;
        }
    }
    else /* Start in middle, integrate backward and forward, stitch together */
    {
        REAL8TimeSeries *V1=NULL, *Phi1=NULL, *S1x1=NULL, *S1y1=NULL, *S1z1=NULL, *S2x1=NULL, *S2y1=NULL, *S2z1=NULL;
        REAL8TimeSeries *LNhatx1=NULL, *LNhaty1=NULL, *LNhatz1=NULL, *E1x1=NULL, *E1y1=NULL, *E1z1=NULL;
        REAL8TimeSeries *V2=NULL, *Phi2=NULL, *S1x2=NULL, *S1y2=NULL, *S1z2=NULL, *S2x2=NULL, *S2y2=NULL, *S2z2=NULL;
        REAL8TimeSeries *LNhatx2=NULL, *LNhaty2=NULL, *LNhatz2=NULL, *E1x2=NULL, *E1y2=NULL, *E1z2=NULL;

        /* Integrate backward to fStart */
        fS = fRef;
        fE = fStart;
        n = XLALSimInspiralSpinTaylorPNEvolveOrbit(&V1, &Phi1,
                &S1x1, &S1y1, &S1z1, &S2x1, &S2y1, &S2z1, 
                &LNhatx1, &LNhaty1, &LNhatz1, &E1x1, &E1y1, &E1z1,
                deltaT, m1, m2, fS, fE, s1x, s1y, s1z, s2x, s2y,
                s2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2,
                quadparam1, quadparam2, spinO, tideO, phaseO, approx);
        
        /* Apply phase shift so orbital phase has desired value at fRef */
        phiShift = phiRef - Phi1->data->data[Phi1->data->length-1];
        for( i=0; i < Phi1->data->length; i++)
        {
            Phi1->data->data[i] += phiShift;
        }

        /* Integrate forward to end of waveform */
        fS = fRef;
        fE = 0.;
        n = XLALSimInspiralSpinTaylorPNEvolveOrbit(&V2, &Phi2,
                &S1x2, &S1y2, &S1z2, &S2x2, &S2y2, &S2z2, 
                &LNhatx2, &LNhaty2, &LNhatz2, &E1x2, &E1y2, &E1z2,
                deltaT, m1, m2, fS, fE, s1x, s1y, s1z, s2x, s2y,
                s2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2,
                quadparam1, quadparam2, spinO, tideO, phaseO, approx);
        
        /* Apply phase shift so orbital phase has desired value at fRef */
        phiShift = phiRef - Phi2->data->data[0];
        for( i=0; i < Phi2->data->length; i++)
        {
            Phi2->data->data[i] += phiShift;
        }

        /* Stitch 2nd set of vectors onto 1st set. Free 2nd set. */
        V = appendTSandFree(V1, V2); 
        Phi = appendTSandFree(Phi1, Phi2);
        S1x = appendTSandFree(S1x1, S1x2);
        S1y = appendTSandFree(S1y1, S1y2);
        S1z = appendTSandFree(S1z1, S1z2);
        S2x = appendTSandFree(S2x1, S2x2);
        S2y = appendTSandFree(S2y1, S2y2);
        S2z = appendTSandFree(S2z1, S2z2);
        LNhatx = appendTSandFree(LNhatx1, LNhatx2);
        LNhaty = appendTSandFree(LNhaty1, LNhaty2);
        LNhatz = appendTSandFree(LNhatz1, LNhatz2);
        E1x = appendTSandFree(E1x1, E1x2);
        E1y = appendTSandFree(E1y1, E1y2);
        E1z = appendTSandFree(E1z1, E1z2);
    }


    /* Use the dynamical variables to build the polarizations */
    status = XLALSimInspiralPrecessingPolarizationWaveforms(hplus, hcross,
            V, Phi, S1x, S1y, S1z, S2x, S2y, S2z, LNhatx, LNhaty, LNhatz, 
            E1x, E1y, E1z, m1, m2, r, v0, amplitudeO);

    /* Destroy vectors of dynamical variables, check for errors then exit */
    XLALDestroyREAL8TimeSeries(V);
    XLALDestroyREAL8TimeSeries(Phi);
    XLALDestroyREAL8TimeSeries(S1x);
    XLALDestroyREAL8TimeSeries(S1y);
    XLALDestroyREAL8TimeSeries(S1z);
    XLALDestroyREAL8TimeSeries(S2x);
    XLALDestroyREAL8TimeSeries(S2y);
    XLALDestroyREAL8TimeSeries(S2z);
    XLALDestroyREAL8TimeSeries(LNhatx);
    XLALDestroyREAL8TimeSeries(LNhaty);
    XLALDestroyREAL8TimeSeries(LNhatz);
    XLALDestroyREAL8TimeSeries(E1x);
    XLALDestroyREAL8TimeSeries(E1y);
    XLALDestroyREAL8TimeSeries(E1z);
    if( status < 0 )
        XLAL_ERROR(XLAL_EFUNC);

    return n;
}

/**
 * Driver routine to compute the physical template family "Q" vectors using
 * the \"T4\" method. Note that PTF describes single spin systems
 *
 * This routine requires leading-order amplitude dependence
 * but allows the user to specify the phase PN order
 */
int XLALSimInspiralSpinTaylorT4PTFQVecs(
        REAL8TimeSeries **Q1,            /**< Q1 output vector */
        REAL8TimeSeries **Q2,            /**< Q2 output vector */
        REAL8TimeSeries **Q3,            /**< Q3 output vector */
        REAL8TimeSeries **Q4,            /**< Q4 output vector */
        REAL8TimeSeries **Q5,            /**< Q5 output vector */
        REAL8 deltaT,                   /**< sampling interval (s) */
        REAL8 m1,                       /**< mass of companion 1 (kg) */
        REAL8 m2,                       /**< mass of companion 2 (kg) */
        REAL8 chi1,                     /**< spin magnitude (|S1|) */
        REAL8 kappa1,                    /**< L . S1 (1 if they are aligned) */
        REAL8 fStart,                   /**< start GW frequency (Hz) */
        REAL8 lambda1,                  /**< (tidal deformability of mass 1) / (mass of mody 1)^5 (dimensionless) */
        REAL8 lambda2,                  /**< (tidal deformability of mass 2) / (mass of body 2)^5 (dimensionless) */
        LALSimInspiralSpinOrder spinO,  /**< twice PN order of spin effects */
        LALSimInspiralTidalOrder tideO, /**< twice PN order of tidal effects */
        int phaseO                      /**< twice PN phase order */
        )
{
    /* To generate the QVecs we need to choose a specific frame 
     * This frame is set so that inclination, and most other extrinsic
     * angles are 0. This does not lead to loss in generality as PTF maximizes
     * over these angles. This follows the PBCV convention
     */
    Approximant approx = SpinTaylorT4;
    REAL8 fRef = 0., quadparam1 = 1., quadparam2 = 1.;
    REAL8 r = 10E6 * LAL_PC_SI; /* Setting an arbitrary distance of 10 MPc */
    REAL8 s1x = chi1 * pow((1 - kappa1*kappa1),0.5);
    REAL8 s1z = chi1 * kappa1;
    REAL8 s1y,s2x,s2y,s2z,lnhatx,lnhaty,lnhatz,e1x,e1y,e1z;
    s1y = s2x = s2y = s2z = lnhatx = lnhaty = e1y = e1z = 0;     
    lnhatz = e1x = 1.;

    REAL8TimeSeries *V, *Phi, *S1x, *S1y, *S1z, *S2x, *S2y, *S2z;
    REAL8TimeSeries *LNhatx, *LNhaty, *LNhatz, *E1x, *E1y, *E1z;
    int status, n;

    /* Evolve the dynamical variables */
    n = XLALSimInspiralSpinTaylorPNEvolveOrbit(&V, &Phi, &S1x, &S1y, &S1z,
            &S2x, &S2y, &S2z, &LNhatx, &LNhaty, &LNhatz, &E1x, &E1y, &E1z,
            deltaT, m1, m2, fStart, fRef, s1x, s1y, s1z, s2x, s2y,
            s2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2,
            quadparam1, quadparam2, spinO, tideO, phaseO, approx);
    if( n < 0 )
        XLAL_ERROR(XLAL_EFUNC);

    /* Use the dynamical variables to build the polarizations */
    status = XLALSimInspiralPrecessingPTFQWaveforms(Q1, Q2, Q3, Q4, Q5,
            V, Phi, S1x, S1y, S1z, S2x, S2y, S2z, LNhatx, LNhaty, LNhatz,
            E1x, E1y, E1z, m1, m2, r);

    /* Destroy vectors of dynamical variables, check for errors then exit */
    XLALDestroyREAL8TimeSeries(V);
    XLALDestroyREAL8TimeSeries(Phi);
    XLALDestroyREAL8TimeSeries(S1x);
    XLALDestroyREAL8TimeSeries(S1y);
    XLALDestroyREAL8TimeSeries(S1z);
    XLALDestroyREAL8TimeSeries(S2x);
    XLALDestroyREAL8TimeSeries(S2y);
    XLALDestroyREAL8TimeSeries(S2z);
    XLALDestroyREAL8TimeSeries(LNhatx);
    XLALDestroyREAL8TimeSeries(LNhaty);
    XLALDestroyREAL8TimeSeries(LNhatz);
    XLALDestroyREAL8TimeSeries(E1x);
    XLALDestroyREAL8TimeSeries(E1y);
    XLALDestroyREAL8TimeSeries(E1z);
    if( status < 0 )
        XLAL_ERROR(XLAL_EFUNC);

    return n;
}


