/*
 * Copyright (C) 2026 Khun Sang Phukon, Nathan Johnson-McDaniel, Amitesh Singh, Anuradha Gupta
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

#include <math.h>
#include <lal/Units.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/LALSimInspiral.h>
#include <lal/Date.h>
#include <lal/LALAdaptiveRungeKuttaIntegrator.h>
#include "LALSimInspiralSpinEccentricityEvolution.h"
#include <lal/XLALGSL.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif


#define VISCO 1.L/sqrt(6.L)

/* use error codes above 1024 to avoid conflicts with GSL */
#define LALSIMINSPIRAL_STE_TEST_ENERGY                  1025
#define LALSIMINSPIRAL_STE_TEST_OMEGADOUBLEDOT          1026
#define LALSIMINSPIRAL_STE_TEST_OMEGANAN                1028
#define LALSIMINSPIRAL_STE_TEST_FREQBOUND               1029
#define LALSIMINSPIRAL_STE_DERIVATIVE_OMEGANONPOS       1030
#define LALSIMINSPIRAL_STE_TEST_LARGEXBAR               1031
#define LALSIMINSPIRAL_STE_TEST_ANGMOMENTUM             1032
#define LALSIMINSPIRAL_STE_TEST_ETSQBOUND               1033
#define LALSIMINSPIRAL_STE_UNBOUND_SYSTEM               1034
#define LALSIMINSPIRAL_STE_NEGATIVE_MAG_ANGMOMENTUM     1035
#define LALSIMINSPIRAL_STE_NEGATIVE_PERIASTRON_ADVANCE  1036


#define LAL_ST4_ABSOLUTE_TOLERANCE 1.e-12
#define LAL_ST4_RELATIVE_TOLERANCE 1.e-12

/* The following definitions are from lal/lib/utilities/LALAdaptiveRungeKuttaIntegrator.c */
#define XLAL_BEGINGSL \
        { \
          gsl_error_handler_t *saveGSLErrorHandler_; \
          saveGSLErrorHandler_ = gsl_set_error_handler_off();

#define XLAL_ENDGSL \
          gsl_set_error_handler( saveGSLErrorHandler_ ); \
        }

static bool unit_vector_tolerance(REAL8 vx, REAL8 vy, REAL8 vz, REAL8 tolerance_range);

static int XLALSimInspiralSpinEccDynamicalExpnCoeffsSpinEccVariables( expnCoeffsSpinEcc *ak, REAL8 omega, REAL8 etSq,
                                                        REAL8 s1x, REAL8 s1y, REAL8 s1z,
														REAL8 s2x, REAL8 s2y, REAL8 s2z,
														REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz,
														REAL8 phatx, REAL8 phaty, REAL8 phatz);

static int evaluate_depsilon_dt ( REAL8 *epsilon, REAL8 *depsilon_dt, REAL8 xbar, REAL8 etSq,
                           REAL8 domega_dt, REAL8 detSq_dt, expnCoeffsSpinEcc *ak  );

static int evaluate_dL_scaled_dt ( REAL8 *Lscaled, REAL8 *dLscaled_dt, REAL8 xbar, REAL8 etSq,
                           REAL8 s1x, REAL8 s1y, REAL8 s1z,
                           REAL8 s2x, REAL8 s2y, REAL8 s2z,
                           REAL8 lhatx, REAL8 lhaty, REAL8 lhatz,
                           REAL8 domega_dt, REAL8 detSq_dt, expnCoeffsSpinEcc *ak  );

static int XLALSimInspiralSpinTaylorPNEccentricStoppingTest( double t,  const double y[], double ydot[], void *params );

static int XLALSimInspiralVectorCrossProduct(REAL8 **vout, REAL8 v1x, REAL8 v1y, REAL8 v1z, REAL8 v2x, REAL8 v2y, REAL8 v2z);

static REAL8 XLALSimInspiralComputeCrossDotQuants( REAL8 *lhatCrosss1normSq, REAL8 *lhatCrosss2normSq, REAL8 *lhatdots1, REAL8 *lhatdots2,
                                                REAL8 s1x, REAL8 s1y, REAL8 s1z,
                                                REAL8 s2x, REAL8 s2y, REAL8 s2z,
                                                REAL8 lhatx, REAL8 lhaty, REAL8 lhatz);

/* Copied from GSL rkf45.c */
typedef struct {
    double *k1;
    double *k2;
    double *k3;
    double *k4;
    double *k5;
    double *k6;
    double *y0;
    double *ytmp;
} rkf45_state_t;

/**
 * @brief Driver function for TaylorT4 evolution of spin-precessing binaries on eccentric orbits
 *
 * This function performs the time evolution of precessing binary black holes on eccentric orbits using
 * orbit-averaged post-Newtonian (PN) evolution. It uses an adaptive Runge-Kutta integrator to solve the
 * system of 17 ordinary differential equations.  The evolved parameters include orbital frequency, eccentricity,
 * mean anomaly, secular contribution to orbital phase, spin vectors, orbital angular momentum vector,
 * periastron line vector, and periastron advance. The orbital elements are described using the quasi-Keplerian
 * parametrization.
 *
 *  The evolution includes:
 * - **Non-spinning post-Newtonian terms**: Up to 3PN order
 * - **Hereditary effects**: Tail and memory contributions via three difference enhancement functions that depend on eccentricity
 * - **Spin effects**: Spin-orbit and spin-spin couplings up to 2PN order
 *
 * @param[out] omega        Timeseries of orbital frequency parameter
 * @param[out] et           Timeseries of orbital eccentricity parameter
 * @param[out] l            Timeseries of mean anomaly of quasi-Keplerian orbit
 * @param[out] lamb         Timeseries of secular contribution to orbital phase
 * @param[out] S1x,S1y,S1z  Timeseries of spin vector components for body 1
 * @param[out] S2x,S2y,S2z  Timeseries of spin vector components for body 2
 * @param[out] LNhatx,LNhaty,LNhatz  Timeseries of orbital angular momentum unit vector components
 * @param[out] Phatx,Phaty,Phatz      Timeseries of periastron line unit vector components
 * @param[out] k            Timeseries of periastron advance parameter
 * @param[in]  deltaT       Sampling interval (s)
 * @param[in]  m1_SI        Mass of companion 1 (kg)
 * @param[in]  m2_SI        Mass of companion 2 (kg)
 * @param[in]  fStart       Starting GW frequency (Hz)
 * @param[in]  fEnd         Ending GW frequency (Hz), 0 means integrate as far as possible
 * @param[in]  s1x,s1y,s1z  Initial dimensionless spin vector components for body 1
 * @param[in]  s2x,s2y,s2z  Initial dimensionless spin vector components for body 2
 * @param[in]  lnhatx,lnhaty,lnhatz  Initial unit orbital angular momentum vector components
 * @param[in]  eccentricity Initial orbital eccentricity
 * @param[in]  l0           Initial mean anomaly (rad)
 * @param[in]  lamb0        Initial secular orbital phase (rad)
 * @param[in]  phatx,phaty,phatz  Initial periastron line unit vector components
 * @param[in]  spinO        Twice PN order of spin effects
 * @param[in]  EccOrder     Eccentricity order in hyperasymptotic enhancement functions;
 *                          even number between 2 and 20, or -1 for the maximum available order
 * @param[in]  phaseO       Twice PN order of non-spinning terms
 * @param[in]  enhancementFunc Enumeration specifying eccentricity enhancement function type;
 *                              0 for EF_Arun_etal (O(e^4)), 1 for EF_LoutrelYunes_SuperAsym (superasymptotic) and 2 for EF_LoutrelYunes_HyperAsym (hyperasymptotic)
 * @param[in]  type         Evolution type: FullSeries=0 or FinalValue=1
 * @param[in]  LALparams    LAL dictionary containing accessory parameters
 *
 * @return XLAL_SUCCESS on successful completion
 *
 * @see Phukon et al., arXiv:2504.20543 for the expressions of derived or collected evolution equations
 * @see Loutrel & Yunes, arXiv:1607.05409 for the superasymptotic and hyperasymptotic enhancement functions used in eccentricity evolution
 */
int XLALSimInspiralSpinTaylorPNEccentricEvolveOrbit( REAL8TimeSeries **omega, REAL8TimeSeries **et,
                                                     REAL8TimeSeries **l, REAL8TimeSeries **lamb,
                                                     REAL8TimeSeries **S1x, REAL8TimeSeries **S1y,
                                                     REAL8TimeSeries **S1z, REAL8TimeSeries **S2x,
                                                     REAL8TimeSeries **S2y, REAL8TimeSeries **S2z,
                                                     REAL8TimeSeries **LNhatx, REAL8TimeSeries **LNhaty,
                                                     REAL8TimeSeries **LNhatz, REAL8TimeSeries **Phatx,
                                                     REAL8TimeSeries **Phaty, REAL8TimeSeries **Phatz,
                                                     REAL8TimeSeries **k, REAL8 deltaT, REAL8 m1_SI,
                                                     REAL8 m2_SI, REAL8 fStart, REAL8 fEnd, REAL8 s1x,
                                                     REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z,
                                                     REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 eccentricity,
                                                     REAL8 l0, REAL8 lamb0, REAL8 phatx, REAL8 phaty, REAL8 phatz,
                                                     LALSimInspiralSpinOrder spinO, LALSimInspiralSpinEccPolyEccOrder EccOrder,
                                                     INT4 phaseO, EnhancementFunction enhancementFunc, EvolutionType type,
                                                     LALDict *LALparams
)
{
    INT4 intreturn;
    LALAdaptiveRungeKuttaIntegrator *integrator = NULL;     /* LAL integrator object */
    REAL8 yinit[LAL_NUM_SPIN_ECC_VARIABLES];                /* initial values of parameters */
    /* intermediate variables */
    UINT4 i, cutlen, len;
    int sgn;
    REAL8 dtStart, dtEnd, lengths, wEnd, Mcsec;
    LIGOTimeGPS tStart = LIGOTIMEGPSZERO;

    INT4 print_warning = 0;

    /* Checking NAN */
    if ( !omega || !et || !l || !lamb|| !S1x || !S1y || !S1z || !S2x || !S2y || !S2z
            || !LNhatx || !LNhaty || !LNhatz || !Phatx || !Phaty || !Phatz || !k )
    {
        XLALPrintError("XLAL Error - %s: NULL(s) in output parameters\n",
                       __func__);
        XLAL_ERROR(XLAL_EINVAL);
    }

   if ( spinO>=7 ) {
      XLALPrintError("XLAL Error - %s: Averaged-orbit dynamics incompatible with spin order %d chosen.\n",__func__, spinO);
      XLAL_ERROR(XLAL_EINVAL);
    }

    /* Check time step */
    if (deltaT <=0){
        XLALPrintError("XLAL Error - %s: time step (deltaT) of evolution is either zero or negative valued\n", __func__);
        XLAL_ERROR(XLAL_EINVAL);
    }

    /* Check mass params */
    if ( m1_SI <= 0 || m2_SI <= 0){
        XLALPrintError("XLAL Error - %s: Masses of component objects are either zero or negative valued\n", __func__);
        XLAL_ERROR(XLAL_EINVAL);
    }

    /* Check eccentricity */
    if (eccentricity>=1 || eccentricity<0){
        XLALPrintError("XLAL Error - %s: Eccentricity values are greater than or equal to 1 or negative\n",
                __func__);
        XLAL_ERROR(XLAL_EINVAL);
    }

    REAL8 tolerance = 1e-6;

    bool lmag_check = unit_vector_tolerance(lnhatx, lnhaty, lnhatz, tolerance);

    /* Check magnitude of the orbital angular momentum unit vector */
    if (lmag_check){
        XLALPrintError("XLAL Error - %s: Unphysical orbital angular momentum unit vector; magnitude is not equal to 1 \n", __func__);
        XLAL_ERROR(XLAL_EINVAL);
    }

    bool phat_mag_check = unit_vector_tolerance(phatx, phaty, phatz, tolerance);

    /* Check magnitude of the periastron line unit vector */
    if (phat_mag_check)
    {
        XLALPrintError("XLAL Error - %s: Unphysical unit vector for the periastron line; magnitude is not equal to 1 \n", __func__);
        XLAL_ERROR(XLAL_EINVAL);
    }

    REAL8 lnhat_dot_phat = lnhatx*phatx + lnhaty*phaty + lnhatz*phatz;

    if (fabs(lnhat_dot_phat) >= tolerance) {
        XLALPrintError("XLAL Error - %s: Orbital angular momentum and periastron line vectors must be orthogonal \n", __func__);
        XLAL_ERROR(XLAL_EINVAL);
    }

    /* Dimensionless spin magnitudes */
    REAL8 chi1 = sqrt(s1x*s1x + s1y*s1y + s1z*s1z);
    REAL8 chi2 = sqrt(s2x*s2x + s2y*s2y + s2z*s2z);

    /* Check spin magnitudes */
    if (chi1 > 1 || chi2 > 1){
        XLALPrintError("XLAL Error - %s: Unphysical spins, chi1 or chi2 > 1\n",
                       __func__);
        XLAL_ERROR(XLAL_EINVAL);
    }

    /* Check if start and end frequencies are positive */
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

    /* Extract settings from LALparams dict */
    EccEvol2PNSpinFlag et2PNSpinFlag = XLALSimInspiralWaveformParamsLookupEccEvol2PNSpinFlag(LALparams);
    EccEvolAlwaysOutput always_output_flag = XLALSimInspiralWaveformParamsLookupEccEvolAlwaysOutput(LALparams);

    /* String to specify the type of operation (interpolation or integration) based on the output flag */
    const char *operation;
    if (always_output_flag){
        operation = "interpolation";
    }
    else{
        operation = "integration";
    }

    /*
     * Declaring variables for structures
     * params:        Structure for passing params to ODEs
     * SpinEccFuncs:  Structure for  ODEs
     * ak:            Structure for static intrisic parameters
     *
     */
    XLALSimInspiralSpinEccPNEvolveOrbitParams  params;
    expnFuncSpinEcc SpinEccFuncs;
    expnCoeffsSpinEcc ak;

    /* Setting up ak and SpinEccFuncs (ODEs) */
    XLALSimInspiralSpinEccODESetup(&ak, &SpinEccFuncs, m1_SI, m2_SI, chi1, chi2, fStart, fEnd,
				   enhancementFunc, EccOrder, spinO, phaseO, et2PNSpinFlag);

    /* ak.beta43 and ak.gamma1 must be computed before initializing other dynamical variables */
    ak.beta43 = XLALSimInspiralSpinEccBeta(4, 3, s1x, s1y, s1z, s2x,
                                        s2y, s2z, lnhatx, lnhaty, lnhatz, &ak);
    ak.gamma1 = XLALSimInspiralSpinEccGamma1(s1x, s1y, s1z, s2x, s2y,  s2z,
                                        lnhatx, lnhaty, lnhatz, &ak );

    params.funcomega = SpinEccFuncs.orbfreq;
    params.funcetSq  = SpinEccFuncs.orbeccSq;
    params.funcl     = SpinEccFuncs.meananom;
    params.funclamb  = SpinEccFuncs.meanorbitphase;
    params.funcs1x   = SpinEccFuncs.compspin1x;
    params.funcs1y   = SpinEccFuncs.compspin1y;
    params.funcs1z   = SpinEccFuncs.compspin1z;
    params.funcs2x   = SpinEccFuncs.compspin2x;
    params.funcs2y   = SpinEccFuncs.compspin2y;
    params.funcs2z   = SpinEccFuncs.compspin2z;
    params.funclx    = SpinEccFuncs.orbangLx;
    params.funcly    = SpinEccFuncs.orbangLy;
    params.funclz    = SpinEccFuncs.orbangLz;
    params.funcpx    = SpinEccFuncs.periastronlinex;
    params.funcpy    = SpinEccFuncs.periastronliney;
    params.funcpz    = SpinEccFuncs.periastronlinez;
    params.funck     = SpinEccFuncs.periastronadvance;
    params.ak        = ak;

    Mcsec = ak.mt * pow( ak.eta, 0.6);

    /*
     * Estimate length of timeseries using quasicircular Newtonian t(f) formula
     * The quasicircular estimate is always longer than the corresponding eccentric estimate
     * Time from freq. = fStart to infinity
     *
     */
    dtStart = (5.0/256.0) * pow(LAL_PI,-8.0/3.0)
            * pow(Mcsec * fStart,-5.0/3.0) / fStart;

    /* Time from freq. = fEnd to infinity. Set to zero if fEnd=0 */
    dtEnd = (fEnd == 0. ? 0. : (5.0/256.0) * pow(LAL_PI,-8.0/3.0)
            * pow(Mcsec * fEnd,-5.0/3.0) / fEnd);

    /* Time in sec from fStart to fEnd. Note it can be positive or negative */
    lengths = dtStart - dtEnd;

    /* Storing initial values */
    yinit[1]  = LAL_PI*fStart*ak.mt;
    yinit[0]  = eccentricity*eccentricity;
    yinit[2]  = l0;
    yinit[3]  = lamb0;
    yinit[4]  = s1x;
    yinit[5]  = s1y;
    yinit[6]  = s1z;
    yinit[7]  = s2x;
    yinit[8]  = s2y;
    yinit[9]  = s2z;
    yinit[10] = lnhatx;
    yinit[11] = lnhaty;
    yinit[12] = lnhatz;
    yinit[13] = phatx;
    yinit[14] = phaty;
    yinit[15] = phatz;

     if ( yinit[1] <= 0.0) /* orbital frequency must be positive! */
    {
        XLALPrintError("XLAL Error - %s: Initial omega must be positive\n",__func__);
        XLAL_ERROR(XLAL_EINVAL);
    }

    REAL8 v_ini = cbrt( yinit[1]);
    REAL8 x_ini = v_ini*v_ini;
    REAL8 xbar_ini = x_ini/(1 - yinit[0]);

    if ( xbar_ini >= 1)
    {
        XLALPrintError("XLAL Error -%s: xbar_ini value is greater than equal to 1.\n", __func__);
        XLAL_ERROR(XLAL_EDOM);
    }

    /*
     * Initial periastron precession
     * Default value is -1, which means initial of periastron precession value will be internally computed
     * otherwise user passed value will be used
     *
     */
    REAL8 K_ini = XLALSimInspiralWaveformParamsLookupInitialPeriastronPrecession(LALparams);

    if ( K_ini == -1){
        yinit[16] = XLALSimInspiralPeriastronPrecession(xbar_ini, yinit[0], &params.ak, spinO, phaseO);
    }
    else{
        if ( K_ini < 0){
        XLALPrintError("XLAL Error - %s: non-negative value of 'InitialPeriastronPrecession' must be provided to LALparam dictionary [exception -1, to be used for computing initial value internally]\n",__func__);
        XLAL_ERROR(XLAL_EINVAL);
        }
        yinit[16] = K_ini;
    }

    XLALSimInspiralSpinEccDynamicalExpnCoeffsSpinEccVariables( &params.ak, yinit[1], yinit[0],
                                                        s1x, s1y, s1z,
                                                        s2x, s2y, s2z,
                                                        lnhatx, lnhaty, lnhatz,
                                                        yinit[13], yinit[14], yinit[15]);

    REAL8  domega_dt_ini =  params.funcomega( yinit[1], yinit[0], yinit[2], yinit[3], s1x, s1y, s1z, s2x, s2y, s2z,
                                              lnhatx, lnhaty, lnhatz, yinit[13], yinit[14], yinit[15], yinit[16], &params.ak);
    REAL8  detSq_dt_ini  =  params.funcetSq( yinit[1], yinit[0], yinit[2], yinit[3], s1x, s1y, s1z, s2x, s2y, s2z,
                                              lnhatx, lnhaty, lnhatz, yinit[13], yinit[14], yinit[15], yinit[16], &params.ak);

    /* Initial epsilon, d_epsilon_dt checks */
    REAL8 epsilon_ini;
    REAL8 depsilon_dt_ini;
    REAL8 Lscaled_ini;
    REAL8 Lscaled_dt_ini;


    evaluate_depsilon_dt ( &epsilon_ini, &depsilon_dt_ini, xbar_ini, yinit[0],
                            domega_dt_ini, detSq_dt_ini, &params.ak);

    evaluate_dL_scaled_dt ( &Lscaled_ini, &Lscaled_dt_ini, xbar_ini, yinit[0],
                            s1x, s1y, s1z, s2x, s2y, s2z, lnhatx, lnhaty, lnhatz,
                            domega_dt_ini, detSq_dt_ini, &params.ak);


    if ( epsilon_ini < 0)
    {
        XLALPrintError("XLAL Error -%s: initial epsilon negative i.e. initial orbit is unbounded.\n", __func__);
        XLAL_ERROR(XLAL_EDOM);
    }

    if ( depsilon_dt_ini <0)
    {
        XLALPrintError("XLAL Error -%s: time derivative of initial epsilon is negative.\n", __func__);
        XLAL_ERROR(XLAL_EDOM);
    }

    if ( Lscaled_ini < 0 )
    {
        XLALPrintError("XLAL Error -%s: initial  angular momentum is negative.\n", __func__);
        XLAL_ERROR(XLAL_EDOM);
    }

    if ( Lscaled_dt_ini > 0)
    {
        XLALPrintError("XLAL Error -%s: time derivative of the initial orbital angular momentum is positive.\n", __func__);
        XLAL_ERROR(XLAL_EDOM);
    }

    if ( yinit[16] < 0){
        XLALPrintError("XLAL Error -%s: negative initial value of periastron precession parameter k \n",__func__);
        XLAL_ERROR(XLAL_EINVAL);
    }

    integrator = XLALAdaptiveRungeKutta4Init(LAL_NUM_SPIN_ECC_VARIABLES,
					       XLALSimInspiralSpinEccentricTaylorT4DerivativesAvg,
					       XLALSimInspiralSpinTaylorPNEccentricStoppingTest,
					       LAL_ST4_ABSOLUTE_TOLERANCE, LAL_ST4_RELATIVE_TOLERANCE);

    /* stop the integration only when the test is true */
    integrator->stopontestonly = 1;

    /* A void pointer to 'params' structure */
    void *void_params = &params;
    REAL8Array *yout;      /* time series of variables returned from integrator */
    switch (type)
    {
        case(FullSeries): ;
            /* run the integration; note: time is measured in \hat{t} = t / M */
            len = XLALAdaptiveRungeKutta4Hermite(integrator, void_params, yinit,
                                    0.0, lengths/ak.mt, sgn*deltaT/ak.mt, &yout);
            if (!yout)
            {
            XLALPrintError("XLAL Error - %s: integration failed (yout == NULL)\n",
                                                               __func__);
                XLAL_ERROR(XLAL_EFUNC);
            }
            break;

        case(FinalValue):
	    if (!always_output_flag) {
	      len = XLALAdaptiveRungeKutta4HermiteOnlyFinal(integrator, void_params, yinit,
                                   0.0, lengths/ak.mt, LAL_PI * ak.mt * fEnd, sgn*deltaT/ak.mt );
	    } else {
	      len = XLALAdaptiveRungeKutta4HermiteOnlyFinalCheckStopping(integrator, void_params, yinit,
                                   0.0, lengths/ak.mt, sgn*deltaT/ak.mt );
	    }
            break;
        default:
            XLALPrintError("XLAL Error - %s: Choose EvolutionType either  FullSeries=0, or  FinalValue=1\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
        break;
    }


    intreturn = integrator->returncode;
    XLALAdaptiveRungeKuttaFree(integrator);

    if (!len)
    {
        XLALPrintError("XLAL Error - %s: integration failed with errorcode %d.\n", __func__, intreturn);
        XLAL_ERROR(XLAL_EFUNC);
    }

    /*
     * Print an error (or if always_output_flag is on, set things to later print a warning) about abnormal termination (excluding the energy and angular momentum flux, xbar > 1,
     * and second time derivative of omega tests only if we are evolving forward as far as possible, with fEnd == 0)
     *
     */
    if (intreturn != 0 && !(fEnd == 0 && (intreturn == LALSIMINSPIRAL_STE_TEST_ENERGY || intreturn == LALSIMINSPIRAL_STE_TEST_ANGMOMENTUM ||
					  intreturn == LALSIMINSPIRAL_STE_TEST_LARGEXBAR || intreturn == LALSIMINSPIRAL_STE_TEST_OMEGADOUBLEDOT))
            && intreturn != LALSIMINSPIRAL_STE_TEST_FREQBOUND)
    {
      if (!always_output_flag) {
        XLALPrintError("XLAL Error - %s: integration terminated with code %d.\n Evolution parameters were m1 = %e, m2 = %e, s1 = (%e,%e,%e), s2 = (%e,%e,%e), ecc = %e\n", __func__, intreturn, m1_SI / LAL_MSUN_SI, m2_SI / LAL_MSUN_SI, s1x, s1y, s1z, s2x, s2y, s2z, eccentricity);
        XLAL_ERROR(XLAL_EFUNC);
      } else {
	print_warning = 1;
      }
    }

    /*
     * If ending frequency was non-zero, we may have overshot somewhat.
     * The integrator takes one adaptive stride past fEnd,
     * but this may include several smaller interpolation steps.
     * Therefore, 'cutlen' will be the index of the first interpolated step
     * to cross fEnd and 'len' is the full length returned from the integrator.
     * If fEnd == 0, we integrated as far as possible and 'cutlen' = 'len'.
     *
     */
    cutlen = len;

    switch (type)
    {
        case(FullSeries): ;
            REAL8 fTerm;
            if( fEnd != 0. && fEnd < fStart )
            {
                wEnd = LAL_PI * ak.mt * fEnd;/* Ending dimensionless freq. 'hat omega' */
                /*
                 * Integrator returns 'hat omega' in array 'yout'
                 * in range data[len] to data[len+(len-1)].
                 * Start at end and find where we cross wEnd
                 *
                 */
                while( yout->data[2*len+cutlen-1] < wEnd )
                    cutlen--;
                if( cutlen < len )
                    cutlen++; /* while loop exits on wrong side of fEnd, so increment */
            }
            else if( fEnd > fStart )
            {
                wEnd = LAL_PI * ak.mt * fEnd;/* Ending dimensionless freq. 'hat omega' */
                /*
                 * Integrator returns 'hat omega' in array 'yout'
                 * in range data[len] to data[len+(len-1)].
                 * Start at end and find where we cross wEnd
                 *
                 */
                while( yout->data[2*len+cutlen-1] > wEnd )
                    cutlen--;
                if( cutlen < len )
                    cutlen++; /* while loop exits on wrong side of fEnd, so increment */
            }

            /* Adjust tStart so last sample is at time=0 */
            XLALGPSAdd(&tStart, -1.0*(cutlen-1)*deltaT);

            /*
             * Report termination condition and final frequency
             * Will only report this info if '4' bit of lalDebugLevel is 1
             *
             */
            fTerm = yout->data[2*len+cutlen-1] / LAL_PI / ak.mt;
	        if (!print_warning) {
                XLALPrintInfo("XLAL Info - %s: integration terminated with code %d. The final GW frequency reached was %g\n", __func__, intreturn, fTerm);
            } else {
                XLALPrintWarning("XLAL Warning - %s: integration terminated with code %d.\n The final GW frequency reached was %g, while the desired final frequency is %g.\n Evolution parameters were m1 = %e, m2 = %e, s1 = (%e,%e,%e), s2 = (%e,%e,%e), ecc = %e\n", __func__, intreturn, fTerm, fEnd, m1_SI / LAL_MSUN_SI, m2_SI / LAL_MSUN_SI, s1x, s1y, s1z, s2x, s2y, s2z, eccentricity);
            }
            break;

        case(FinalValue):
            fTerm = yinit[1] / LAL_PI / ak.mt;
	    if (!print_warning) {
	      XLALPrintInfo("XLAL Info - %s: %s terminated with code %d. The final GW frequency reached was %g\n", __func__, operation, intreturn, fTerm);
	    } else {
	      XLALPrintWarning("XLAL Warning - %s: %s terminated with code %d.\n The final GW frequency reached was %g, while the desired final frequency is %g.\n Evolution parameters were m1 = %e, m2 = %e, s1 = (%e,%e,%e), s2 = (%e,%e,%e), ecc = %e\n", __func__, operation, intreturn, fTerm, fEnd, m1_SI / LAL_MSUN_SI, m2_SI / LAL_MSUN_SI, s1x, s1y, s1z, s2x, s2y, s2z, eccentricity);
	    }
            break;
        default:
            XLALPrintError("XLAL Error - %s: Choose EvolutionType either  FullSeries=0, or  FinalValue=1\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    /* Memory allocation */
    *omega = XLALCreateREAL8TimeSeries( "ORBITAL_FREQUENCY", &tStart, 0.,
            deltaT, &lalDimensionlessUnit, cutlen);
    *et = XLALCreateREAL8TimeSeries( "ECCENTRICITY_PARAMETER", &tStart, 0.,
            deltaT, &lalDimensionlessUnit, cutlen);
    *l = XLALCreateREAL8TimeSeries( "MEAN_ANOMALY_PARAMETER", &tStart, 0.,
            deltaT, &lalDimensionlessUnit, cutlen);
    *lamb = XLALCreateREAL8TimeSeries( "MEAN_ORBITAL_PHASE_PARAMETER", &tStart, 0.,
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
    *Phatx = XLALCreateREAL8TimeSeries( "PERIASTRON_LINE_X_COMPONENT", &tStart, 0.,
            deltaT, &lalDimensionlessUnit, cutlen);
    *Phaty = XLALCreateREAL8TimeSeries( "PERIASTRON_LINE_Y_COMPONENT", &tStart, 0.,
            deltaT, &lalDimensionlessUnit, cutlen);
    *Phatz = XLALCreateREAL8TimeSeries( "PERIASTRON_LINE_Z_COMPONENT", &tStart, 0.,
            deltaT, &lalDimensionlessUnit, cutlen);
    *k = XLALCreateREAL8TimeSeries( "PERIASTRON_PRECESSION", &tStart, 0.,
            deltaT, &lalDimensionlessUnit, cutlen);

    /* Checking NAN */
    if ( !*omega || !*et || !*l || !*lamb|| !*S1x || !*S1y || !*S1z || !*S2x || !*S2y || !*S2z
            || !*LNhatx || !*LNhaty || !*LNhatz || !*Phatx || !*Phaty || !*Phatz || !*k )
    {

        XLALDestroyREAL8Array(yout);
        XLAL_ERROR(XLAL_EFUNC);

    }

    switch (type)
    {
        case(FullSeries): ;
            int offset;
            /* If we integrated backwards, offset & sgn will reverse order of samples */
            if( fEnd < fStart && fEnd != 0. )
                offset = cutlen-1;
            else
                offset = 0;

           /*
            * Copy dynamical variables from yout array to output time series.
            * Note the first 'len' members of yout are the time steps.
            * Also, the for loop only goes to 'cutlen', in case we overshot fEnd.
            * If we integrated backwards, we copy backwards from 'cutlen'.
            *
            */
            for( i = 0; i < cutlen; i++ )
            {
                int j = sgn*i+offset;
                (*omega)->data->data[j]    = yout->data[2*len+i];
                (*et)->data->data[j] 	   = sqrt(yout->data[len+i]);
                (*l)->data->data[j] 	   = yout->data[3*len+i];
                (*lamb)->data->data[j] 	   = yout->data[4*len+i];
                (*S1x)->data->data[j]      = yout->data[5*len+i];
                (*S1y)->data->data[j] 	   = yout->data[6*len+i];
                (*S1z)->data->data[j]      = yout->data[7*len+i];
                (*S2x)->data->data[j] 	   = yout->data[8*len+i];
                (*S2y)->data->data[j] 	   = yout->data[9*len+i];
                (*S2z)->data->data[j] 	   = yout->data[10*len+i];
                (*LNhatx)->data->data[j]   = yout->data[11*len+i];
                (*LNhaty)->data->data[j]   = yout->data[12*len+i];
                (*LNhatz)->data->data[j]   = yout->data[13*len+i];
                (*Phatx)->data->data[j]    = yout->data[14*len+i];
                (*Phaty)->data->data[j]    = yout->data[15*len+i];
                (*Phatz)->data->data[j]    = yout->data[16*len+i];
                (*k)->data->data[j]        = yout->data[17*len+i];
            }
            XLALDestroyREAL8Array(yout);
            break;
        case(FinalValue):
            (*omega)->data->data[0] = yinit[1];
            (*et)->data->data[0]    = sqrt(yinit[0]);
            (*l)->data->data[0]     = yinit[2];
            (*lamb)->data->data[0]  = yinit[3];
            (*S1x)->data->data[0]   = yinit[4];
            (*S1y)->data->data[0] 	= yinit[5];
            (*S1z)->data->data[0]   = yinit[6];
            (*S2x)->data->data[0] 	= yinit[7];
            (*S2y)->data->data[0] 	= yinit[8];
            (*S2z)->data->data[0] 	= yinit[9];
            (*LNhatx)->data->data[0]= yinit[10];
            (*LNhaty)->data->data[0]= yinit[11];
            (*LNhatz)->data->data[0]= yinit[12];
            (*Phatx)->data->data[0] = yinit[13];
            (*Phaty)->data->data[0] = yinit[14];
            (*Phatz)->data->data[0] = yinit[15];
            (*k)->data->data[0]     = yinit[16];
            break;
       default:
            XLALPrintError("XLAL Error - %s: Choose EvolutionType either  FullSeries=0, or  FinalValue=1\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }
    return XLAL_SUCCESS;
 }

/**
 * @brief Check if a vector has unit magnitude within a specified tolerance
 *
 * This function verifies whether a given vector is a unit vector
 * within the specified tolerance range. It returns true if the vector's magnitude deviates from 1
 * by more than the tolerance.
 *
 * @param[in] vx,vy,vz  Cartesian components of the vector to check
 * @param[in] tolerance_range  Tolerance for the magnitude check
 *
 * @return True if the vector magnitude deviates from 1 by more than the tolerance
 *
 */
static bool unit_vector_tolerance(REAL8 vx, REAL8 vy, REAL8 vz, REAL8 tolerance_range)
{
    REAL8 vec_mag = sqrt(vx*vx + vy*vy + vz*vz);
    return (fabs(vec_mag - 1) >= fabs(tolerance_range));
}

/**
 * @brief Function to setup for TaylorT4 evolution spin-precessing eccentric binary
 *
 * The function setup the ODEs and static parameter for the TaylorT4 evolution. It computes
 * static quantities such as mass ratios, symmetric mass ratio, and stores in the expnCoeffsSpinEcc
 * structure. It also assigns function pointers of 17 ODEs to the expnFuncSpinEcc structure for use by the ODE integrator.
 *
 * @param[out] ak        Structure containing static parameters for the TaylorT4 evolution
 * @param[out] f         Structure containing function pointers of evolution equations
 * @param[in]  m1_kg     Mass of body 1 in kg
 * @param[in]  m2_kg     Mass of body 2 in kg
 * @param[in]  chi1      Dimensionless spin magnitude of body 1
 * @param[in]  chi2      Dimensionless spin magnitude of body 2
 * @param[in]  fstart    Starting GW wave frequency for the evolution in Hz
 * @param[in]  fend      Ending GW wave frequency for the evolution in Hz
 * @param[in]  enhancementFunc  Enumeration specifying which enhancement function to use for hereditary effects
 * @param[in]  EccOrder  Order of eccentricity in the enhancement functions (for hyperasymptotic enhancement functions)
 * @param[in]  SpinOrder Twice the post-Newtonian order of spin effects (e.g., 4 for 2PN spin terms)
 * @param[in]  PhaseOrder Twice the post-Newtonian order of non-spinning phase terms
 * @param[in]  et2PNSpinFlag Flag indicating whether to include 2PN spin terms in eccentricity evolution
 *
 * @return 0 on success
 *
 */
int XLALSimInspiralSpinEccODESetup( expnCoeffsSpinEcc *ak, expnFuncSpinEcc *f, REAL8 m1_kg, REAL8 m2_kg,
                                    REAL8 chi1, REAL8 chi2, REAL8 fstart, REAL8 fend,
                                    EnhancementFunction enhancementFunc, LALSimInspiralSpinEccPolyEccOrder EccOrder,
                                    LALSimInspiralSpinOrder SpinOrder, INT4 PhaseOrder, EccEvol2PNSpinFlag et2PNSpinFlag )
{
    ak->q = m2_kg/m1_kg;
    ak->m1 = 1/(1 + ak->q);
    ak->m2 = ak->q/(1 + ak->q);
    ak->m = ak->m1 + ak->m2;
    ak->eta = ak->m1 * ak-> m2 / (ak->m*ak->m);
    ak->mt = (m1_kg + m2_kg)/LAL_MSUN_SI*LAL_MTSUN_SI;
    ak->norm1 = ak->m1 * ak->m1;
    ak->norm2 = ak->m2 * ak->m2;
    ak->chi1 = chi1;
    ak->chi2 = chi2;
    ak->enhancementFunc = enhancementFunc;
    ak->EccOrder = EccOrder;
    ak->fstart = fstart;
    ak->fend = fend;
    ak->SpinOrder = SpinOrder;
    ak->PhaseOrder = PhaseOrder;
    ak->et2PNSpinFlag = et2PNSpinFlag;
    /* Initial value of domega/dt, which will be used in stopping test */
    ak->prev_domega = 0;

    f->orbfreq = &XLALSimInspiralOmegaEvolution;
    f->orbeccSq = &XLALSimInspiralOrbitalEccentricitySqEvolution;
    f->meananom = &XLALSimInspiralMeanAnomalyEvolution;
    f->meanorbitphase = &XLALSimInspiralMeanOrbitalPhaseEvolution;
    f->compspin1x = &XLALSimInspiralS1xEvolution;
    f->compspin1y = &XLALSimInspiralS1yEvolution;
    f->compspin1z = &XLALSimInspiralS1zEvolution;
    f->compspin2x = &XLALSimInspiralS2xEvolution;
    f->compspin2y = &XLALSimInspiralS2yEvolution;
    f->compspin2z = &XLALSimInspiralS2zEvolution;
    f->orbangLx = &XLALSimInspiralLhatxEvolution;
    f->orbangLy = &XLALSimInspiralLhatyEvolution;
    f->orbangLz = &XLALSimInspiralLhatzEvolution;
    f->periastronlinex = &XLALSimInspiralPeriastronLinexEvolution;
    f->periastronliney = &XLALSimInspiralPeriastronLineyEvolution;
    f->periastronlinez = &XLALSimInspiralPeriastronLinezEvolution;
    f->periastronadvance = &XLALSimInspiralPeriastronPrecessionEvolution;
    return 0;
}

/**
 * @brief Compute the time derivative of orbital angular frequency
 *
 * @param[in] omega Orbital angular velocity in units total mass M = 1
 * @param[in] etSq Square of the eccentricity (e²)
 * @param[in] l Mean anomaly (unused)
 * @param[in] lambda Secular orbital phase contribution (unused)
 * @param[in] s1x, s1y, s1z Cartesian components of dimensionless spin of body 1
 * @param[in] s2x, s2y, s2z Cartesian components of dimensionless spin of body 2
 * @param[in] lhatx, lhaty, lhatz Unit vector components of orbital angular momentum
 * @param[in] phatx, phaty, phatz Periastron direction (unused)
 * @param[in] k Periastron advance rate (unused)
 * @param[in] ak Structure containing static parameters, enhancement function, PN order, spin order settings
 *
 * @return Time derivative of orbital angular frequency
 *
 * @see Equations 6a and A1a-A1f from Phukon et al., arXiv:2504.20543.
 */
REAL8 XLALSimInspiralOmegaEvolution( REAL8 omega, REAL8 etSq,  REAL8 UNUSED l, REAL8 UNUSED lambda,
                                    REAL8 s1x,  REAL8 s1y,  REAL8 s1z,  REAL8 s2x, REAL8 s2y, REAL8 s2z,
                                    REAL8 lhatx, REAL8 lhaty, REAL8 lhatz, REAL8 UNUSED phatx,
                                    REAL8 UNUSED phaty, REAL8 UNUSED phatz, REAL8 UNUSED k,
                                    expnCoeffsSpinEcc *ak )
{
    EnhancementFunction enhancementFunc = ak->enhancementFunc;
    LALSimInspiralSpinEccPolyEccOrder EccOrder = ak->EccOrder;
    LALSimInspiralSpinOrder spinO = ak->SpinOrder;
    INT4 phaseO = ak->PhaseOrder;

    REAL8 eta = ak->eta;
    REAL8 PI_Sq = LAL_PI*LAL_PI;

    REAL8 mtotal = ak->m;

    REAL8 v = cbrt(mtotal*omega);
    REAL8 x = v*v;

    REAL8 one_minus_etSq = 1 - etSq;
    REAL8 sqrt_one_minus_etSq = sqrt(one_minus_etSq);
    REAL8 one_minus_etSq_pow2 = one_minus_etSq*one_minus_etSq;
    REAL8 one_minus_etSq_pow5 = pow(one_minus_etSq, 5);
    REAL8 one_minus_etSq_pow6 = one_minus_etSq_pow5*one_minus_etSq;
    REAL8 one_minus_etSq_pow13by2 = one_minus_etSq_pow6*sqrt_one_minus_etSq;

    REAL8 xbar = x/one_minus_etSq;
    REAL8 Sqrtxbar = v/sqrt_one_minus_etSq;
    REAL8 xbar_3by2 = xbar*Sqrtxbar;
    REAL8 xbar_pow11by2 = pow(xbar, 11./2);

    REAL8 log_term = log(xbar * (1 + sqrt_one_minus_etSq )/2);

    /* Equation A1a, Phukon et al, 2504.20543 */
    REAL8 T_N = (96. + etSq*(292. + etSq*37.))/5.;
    /* Equation A1b, Phukon et al, 2504.20543 */
    REAL8 T_1PN = - 1486./35 - 264./5 * eta + etSq * ((2193./7 - 570*eta ) +
                    etSq * ((12217./20 - 5061./10 * eta) +
                    etSq * (11717./280 - 148./5 * eta )));
    /* Non-spinning part of Equation A1d, Phukon et al, 2504.20543 */
    REAL8 T_2PN = (- 11257./945 +  eta * (15677./105 + 944./15*eta)) +
                    etSq * (( - 2960801./945 - eta * ( 2781./5 - 182387./90 * eta)) +
                    etSq * (( -68647./1260 - eta * (1150631./140  - 396443./72 * eta )) +
                    etSq * (( 925073./336 - eta * (199939./48 - 192943./90 * eta )) +
                    etSq * ( 391457./3360 - eta * (6037./56 - 2923./45 * eta))))) +
                    sqrt_one_minus_etSq * ( (48- 96./5 * eta ) +
                        etSq * ( (2134 - 4268./5 * eta ) +
                        etSq * ( ( 2193 - 4386./5 * eta ) +
                        etSq*( 175./2 - 35 * eta ))));
    /* Equation A1f (excluding terms with hereditory contributions or the last parenthesized terms), Phukon et al, 2504.20543 */
    REAL8 T_3PN = 614389219./148500 + eta * ((-57265081./11340 + 369./2 * PI_Sq) - eta*( 16073./140 + 1121./27*eta)) +
                    etSq * ((19769277811./693000 + eta *(( 66358561./3240  + 42571./80 * PI_Sq ) -
                                    eta * (3161701./840 + 1287385./324 * eta) )) +
                    etSq * (( - 3983966927./8316000 + eta * (( 6451690597./ 90720  - 12403./64 * PI_Sq) +
                                eta * (34877019./1120  - 33769597./1296 * eta ))) +
                    etSq * (( -4548320963./5544000 + eta * ((-59823689./4032 - 242563./640 * PI_Sq ) +
                                eta * ( 411401857./6720 - 3200965./108 * eta ))) +
                    etSq * (( 19593451667./ 2464000 + eta *((-6614711./480 - 12177./640 * PI_Sq) +
                                eta * (92762./7 - 982645./162 * eta ))) +
                    etSq * ( 33332681./197120 - eta*( 1874543./10080 - eta * (109733./840  - 8288./81 * eta)) ))))) +
                    sqrt_one_minus_etSq * ( (-1425319./1125 + eta * (( 9874./105 - 41./10 * PI_Sq ) + 632./5 * eta)) +
                            etSq *((933454./375 + eta*((-2257181./63 + 45961./240 * PI_Sq ) + 125278./15 * eta)) +
                            etSq * ((840635951./21000 + eta *((-4927789./60 + 6191./32 * PI_Sq ) + 317273./15 * eta )) +
                            etSq * ((702667207./31500 + eta *((-6830419./252  + 287./960 * PI_Sq ) + 232177./30 * eta)) +
                            etSq * ( 56403./112- eta * (427733./840 - 4739./30 * eta)))))) +
                    log_term * (54784./175 +  etSq * ( 465664./105 +
                                   etSq * ( 4426376./ 525  +
                                   etSq * ( 1498856./525  + etSq * 31779./ 350))));

    REAL8 Phi, Kappa;

    switch ( enhancementFunc )
    {
        case (EF_Arun_etal):
            Phi = XLALSimInspiralSpinEccPhi(etSq);
            Kappa = XLALSimInspiralSpinEccEnhancementKappa(etSq);
            break;
        case (EF_LoutrelYunes_SuperAsym):
            Phi = XLALSimInspiralSpinEccPhiSuperAsyLY(etSq);
            Kappa = XLALSimInspiralSpinEccEnhancementKappaSuperAsyLY(etSq);
            break;
        case (EF_LoutrelYunes_HyperAsym):
            Phi = XLALSimInspiralSpinEccPhiHyperAsyLY(etSq, EccOrder);
            Kappa = XLALSimInspiralSpinEccEnhancementKappaHyperAsyLY(etSq, EccOrder);
            break;
        default:
            XLALPrintError("XLAL Error in function - %s: %d is an invalid EnhancementFunction\n",__func__, enhancementFunc);
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }
    REAL8 r0 = mtotal;
    REAL8 PsiOmega = XLALSimInspiralSpinEccPsiOmega( etSq, enhancementFunc, EccOrder );
    REAL8 ZetaOmega = XLALSimInspiralSpinEccZetaOmega (etSq, enhancementFunc, EccOrder );
    REAL8 F = ak->Fe;
    /* First term from Equation A1c of Phukon et al, 2504.20543 */
    REAL8 T_1p5PN_hered_NS = 384./5 * (LAL_PI*Phi) * one_minus_etSq_pow5;
    /* Equation A1e of Phukon et al, 2504.20543 */
    REAL8 T_2p5PN_hered_NS = - LAL_PI * ( 4159./35 * PsiOmega + 2268./5 * eta * ZetaOmega ) * one_minus_etSq_pow6;
    /* Last line from Equation A1f of Phukon et al, 2504.20543 */
    REAL8 T_3PN_hered_NS =  ( - 3736352./6125 * Kappa + (512./5 * PI_Sq -  54784./175 * ( LAL_GAMMA +
                log ( 4 * omega * r0 ) )) * F) * one_minus_etSq_pow13by2;

    REAL8 omegadot = 0;
    REAL8 omegadot_nonspinningPN = 0;

    switch (phaseO)
    {
        case(-1):
        case(6):
	        omegadot_nonspinningPN = T_3PN + T_3PN_hered_NS;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(5):
	        omegadot_nonspinningPN *= Sqrtxbar;
            omegadot_nonspinningPN += T_2p5PN_hered_NS;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(4):
	        omegadot_nonspinningPN *= Sqrtxbar;
            omegadot_nonspinningPN += T_2PN;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(3):
            omegadot_nonspinningPN *= Sqrtxbar;
	        omegadot_nonspinningPN += T_1p5PN_hered_NS;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(2):
	        omegadot_nonspinningPN *= Sqrtxbar;
	        omegadot_nonspinningPN += T_1PN;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(1):
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(0):
	        omegadot_nonspinningPN *= xbar;
            omegadot_nonspinningPN += T_N;
            break;
        default:
            XLALPrintError("XLAL Error - %s: - %d is an invalid phase order\n", __func__, phaseO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    /* Common factor from Equation 6a of Phukon et al, 2504.20543 */
    REAL8 prefactor = eta/(mtotal*mtotal)*xbar_pow11by2*one_minus_etSq_pow2;

    omegadot =  prefactor *  omegadot_nonspinningPN;

    REAL8 coeff_a_beta = - 5424 - etSq * (27608 + etSq * (16694 + 585 * etSq ));
    REAL8 coeff_b_beta = - 3600 - etSq * (20928 + etSq * (14658 + 621 * etSq));
    REAL8 beta =  XLALSimInspiralSpinEccBeta(  coeff_a_beta, coeff_b_beta,
                                                s1x, s1y, s1z, s2x, s2y, s2z,
                                                lhatx, lhaty, lhatz, ak );


    REAL8 coeff_a = - 15808 - etSq*(90400 + etSq*(63640 + 3084*etSq));
    REAL8 coeff_b = 46144 + etSq*(259040 + etSq*(179880 + 8532*etSq));
    REAL8 coeff_c = -etSq*(131344 + etSq*(127888 + 7593*etSq));
    REAL8 coeff_a1a2 = 896 + etSq*(8512 + etSq*(7728 + 504*etSq));
    REAL8 coeff_b1b2 = -128 - etSq*(1216 + etSq*(1104 + 72*etSq));
    REAL8 coeff_c1c2 = etSq*(32 + etSq*(160 + 18*etSq));
    REAL8 SigmaK = XLALSimInspiralSpinEccSigmaK ( coeff_a, coeff_b, coeff_c,
                        coeff_a1a2, coeff_b1b2, coeff_c1c2,
                        s1x, s1y, s1z, s2x, s2y, s2z,
                        lhatx, lhaty, lhatz, ak );


    /* Second term from Equation A1c of Phukon et al, 2504.20543 */
    REAL8 omegadot_15PN_Spin = 1./30*beta;
    /* Last term from Equation A1d of Phukon et al, 2504.20543 */
    REAL8 omegadot_2PN_Spin = 1./320*SigmaK;

    REAL8 omegadot_spin = 0;
    switch (spinO)
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            omegadot_spin = Sqrtxbar*omegadot_2PN_Spin;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            omegadot_spin += omegadot_15PN_Spin;
	        omegadot_spin *= xbar_3by2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %d\n",
                                    __func__, spinO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }
    omegadot = omegadot + prefactor*omegadot_spin;
    return omegadot;
    }

/**
 * @brief Compute the time derivative of the square of orbital eccentricity
 *
 * @param[in] omega Orbital angular velocity in units total mass M = 1
 * @param[in] etpow2 Square of the eccentricity (\f$e_t^2\f$)
 * @param[in] l Mean anomaly (unused)
 * @param[in] lambda Secular orbital phase contribution (unused)
 * @param[in] s1x, s1y, s1z Cartesian components of dimensionless spin of body 1
 * @param[in] s2x, s2y, s2z Cartesian components of dimensionless spin of body 2
 * @param[in] lhatx, lhaty, lhatz Unit vector components of orbital angular momentum
 * @param[in] phatx, phaty, phatz Periastron direction (unused)
 * @param[in] k Periastron advance rate (unused)
 * @param[in] ak Structure containing static parameters, enhancement function, PN order, spin order settings
 *
 * @return Time derivative of the square of eccentricity \f$\frac{d(e_t^2)}{dt}\f$
 *
 * @see Equations 6b and A2a-A2d from Phukon et al., arXiv:2504.20543.
 */
REAL8 XLALSimInspiralOrbitalEccentricitySqEvolution ( REAL8 omega, REAL8 etpow2, REAL8 UNUSED l, REAL8 UNUSED lambda,
                                                      REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z,
                                                      REAL8 lhatx, REAL8 lhaty, REAL8 lhatz, REAL8 UNUSED phatx, REAL8 UNUSED phaty,
                                                      REAL8 UNUSED phatz, REAL8 UNUSED k, expnCoeffsSpinEcc *ak )
{
    REAL8 v = cbrt(omega);
    REAL8 x = v*v;
    REAL8 eta = ak -> eta;

    REAL8 LAL_PI_POW2  = LAL_PI * LAL_PI;

    REAL8 one_minus_et2 = 1- etpow2;
    REAL8 xbar = x/one_minus_et2;
    REAL8 one_minus_et2_pow2 = one_minus_et2*one_minus_et2;
    REAL8 one_minus_et2_pow4 = one_minus_et2_pow2*one_minus_et2_pow2;
    REAL8 one_minus_et2_pow5 = one_minus_et2_pow4*one_minus_et2;
    REAL8 sqrt_one_minus_et2 = sqrt(one_minus_et2);
    REAL8 Sqrtxbar = v/sqrt_one_minus_et2;
    REAL8 xbar_pow3by2 = xbar*Sqrtxbar;
    REAL8 xbar_pow4 = pow(xbar,4);
    REAL8 one_minus_et2_pow1p5 = one_minus_et2*sqrt_one_minus_et2;
    REAL8 one_minus_et2_pow2p5 = one_minus_et2_pow1p5*one_minus_et2;
    REAL8 one_minus_et2_pow11by2 = one_minus_et2_pow2p5*one_minus_et2_pow2*one_minus_et2;

    /*
     * Equation A2a of Phukon et al, 2504.20543
     * The factor 2./(1 - etpow2)^(5./2) of the Newtonian term is multiplied in the final expressions of the ODE.
     */
    REAL8 edot_N = etpow2*(304.+ 121.*etpow2)/15.;

    /* Equation A2b of arxiv:2504.20543*/
    REAL8 edot_1PN = etpow2*((-939./35 - 4084./45*eta) +
                       etpow2*(( 29917./105 - 7753./30*eta) +
                       etpow2*(13929./280 - 1664./45*eta)));

    REAL8 PhiE_mult_etpow2 = XLALSimInspiralSpinEccPhiE_mult_etpow2( etpow2, ak->enhancementFunc, ak->EccOrder );

    /* First term of Equation A2c of Phukon et al, 2504.20543*/
    REAL8 edot_15PN = 394./3*PhiE_mult_etpow2*LAL_PI*one_minus_et2_pow4;

    /* Second term of Equation A2c of Phukon et al, 2504.20543*/
    REAL8 a_beta = 19688. + etpow2*(28256. + 2367.*etpow2);
    REAL8 b_beta = 13032. + etpow2*(24270. + 2505.*etpow2);
    REAL8 Beta = XLALSimInspiralSpinEccBeta ( a_beta, b_beta, s1x, s1y, s1z, s2x, s2y,
                                                              s2z, lhatx, lhaty, lhatz, ak );
    REAL8 edot_15PN_spin = - etpow2/90 * Beta;

    /* Non spinning part (excluding terms with the sigma function) of Equation A2d of Phukon et al., 2504.20543*/
    REAL8 edot_2PN = etpow2*((- 961973./1890 + eta*(70967./210 + 752./5*eta)) +
                                   etpow2*(( -3180307./2520 - eta*(1541059./840 - 64433./40*eta)) +
                                   etpow2*(( 23222071./15120 - eta* (13402843./5040 -  127411./90*eta )) +
                                   etpow2*( 420727./3360 - eta*(362071./2520 - 821./9*eta)))) +
                                   sqrt_one_minus_et2*(1336./3 - 2672./15*eta +
                                       etpow2*(2321*(1./2 - 1./5*eta) +
                                       etpow2*(565./6 - 113./3*eta))));


    /*
     * 2PN spin terms from Equation A2d of Phukon et al., 2504.20543
     * The spin contributions are computed when 2PN spin flag in on, otherwise zero
     *
     */
    REAL8 edot_2PN_spin = 0;
    REAL8 a_sigma =  320. - etpow2*(  62752. + etpow2*(101080. + 9420.*etpow2));
    REAL8 b_sigma = - 320. + etpow2*( 179936 + etpow2*(287160 + 26820*etpow2));
    REAL8 c_sigma = - etpow2*(88432. + etpow2*(161872. + 16521*etpow2));
    REAL8 d_sigma = - 640. +  etpow2*( 5440. + etpow2*(11760 + 1080.*etpow2));
    REAL8 e_sigma = 640. + etpow2*( 320. - etpow2*(3120. + 360.*etpow2));
    REAL8 f_sigma = etpow2*( 32 + etpow2*(160. + 18.*etpow2));

    REAL8 Sigma = XLALSimInspiralSpinEccSigmaK ( a_sigma, b_sigma, c_sigma, d_sigma, e_sigma, f_sigma, s1x, s1y, s1z, s2x,
                                                                            s2y, s2z, lhatx, lhaty, lhatz, ak);
    if (ak->et2PNSpinFlag){
       edot_2PN_spin =  1./960*Sigma;
    }
    else{
        edot_2PN_spin = 0;
    }

    REAL8 PsiE_mult_etpow2 =  XLALSimInspiralSpinEccPsiE_mult_etpow2(etpow2, ak->enhancementFunc, ak->EccOrder);
    REAL8 ZetaE_mult_etpow2 = XLALSimInspiralSpinEccZetaE_mult_etpow2(etpow2, ak->enhancementFunc, ak->EccOrder );

    /* Equation A2e of Phukon et al., 2504.20543 */
    REAL8 edot_25PN = -one_minus_et2_pow5*LAL_PI*(55691./210*PsiE_mult_etpow2 + 305072./315*eta*ZetaE_mult_etpow2);

    REAL8 KappaE_mult_etpow2 = XLALSimInspiralSpinEccEnhancementKappaE_mult_etpow2( etpow2, ak->enhancementFunc, ak->EccOrder );
    REAL8 FluxE = XLALSimInspiralSpinEccEnhancementFluxE(etpow2);

    /* Equation A2f of Phukon et al., 2504.20543 */
    REAL8 edot_3PN = etpow2*( (54177075619./6237000  + eta*((7198067./22680  +
                                            1283./10*LAL_PI_POW2 ) - eta*(3000281./2520 + 61001./486*eta))) +
                                        etpow2*(( 6346360709./891000 + eta*(( 9569213./360 + 54001./960*LAL_PI_POW2 ) + eta*(
                                                12478601./15120 - 86910509./19440*eta ))) +
                                        etpow2*(( -126288160777./16632000 + eta*(( 418129451./181440 -
                                                254903./1920*LAL_PI_POW2 ) + eta*( 478808759./20160
                                                - 2223241./180*eta ))) +
                                        etpow2*(( 5845342193./1232000 + eta*(( -98425673./10080 - 6519./640*LAL_PI_POW2 ) +
                                                eta*(6538757./630 - 11792069./2430*eta ))) +
                                        etpow2*( 302322169./1774080 - eta*(1921387./10080 - eta*( 41179./216 -
                                                193396./1215*eta )))))) +
                                        sqrt_one_minus_et2 * ( (-22713049./15750 +
                                                eta*((-5526991./945 + 8323./180*LAL_PI_POW2 ) + 54332./45*eta)) +
                                        etpow2 *(( 89395687./7875 + eta*(( -38295557./1260 + 94177./960*LAL_PI_POW2) +
                                                681989./90*eta)) +
                                        etpow2 *(( 5321445613./378000 + eta*(( - 26478311./1512 +
                                                2501./2880*LAL_PI_POW2 ) + 225106./45*eta )) +
						                etpow2 *( 186961./336 - eta*(289691./504 - 3197./18*eta))))) +
                                                730168./(23625*(1 + sqrt_one_minus_et2)) + ( 1316528./1575 +
                                        etpow2*(4762784./1575 + etpow2*(2294294./1575 + 20437./350*etpow2)))*
                                                log( x*(1 + sqrt_one_minus_et2)/(2*one_minus_et2))) -
                                        one_minus_et2_pow11by2*( ( 89789209./55125 - 1398704./1575*LAL_LN2 +
                                                                156006./175*LAL_LN3 )*KappaE_mult_etpow2 + ( - 12304./45*LAL_PI_POW2 +
                                                                1316528./1575*( LAL_GAMMA +  2*LAL_LN2 + 3./2*log(x)))*FluxE*etpow2 );

    REAL8 e2dot = 0;
    INT4 phaseO = ak->PhaseOrder;
    switch (phaseO)
    {
        case(-1):
        case(6):
            e2dot = edot_3PN;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(5):
            e2dot *= Sqrtxbar;
	        e2dot += edot_25PN;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(4):
            e2dot *= Sqrtxbar;
	        e2dot += edot_2PN;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(3):
            e2dot *= Sqrtxbar;
	        e2dot += edot_15PN;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(2):
	        e2dot *= Sqrtxbar;
	        e2dot += edot_1PN;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(1):
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(0):
            e2dot *= xbar;
            e2dot += edot_N;
            break;
        default:
            XLALPrintError("XLAL Error - %s: - %d is an invalid phase order\n", __func__, phaseO );
            XLAL_ERROR(XLAL_EINVAL);
            break;

    }

    REAL8 e2dot_spin = 0;

    LALSimInspiralSpinOrder spinO = ak->SpinOrder;
    switch (spinO)
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            e2dot_spin = Sqrtxbar*edot_2PN_spin;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
	        e2dot_spin += edot_15PN_spin;
            e2dot_spin *= xbar_pow3by2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %d\n",
                                    __func__, spinO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    /* Prefractor from Equation 6b of Phukon et al., 2504.20543*/
    REAL8 prefactor = -2*eta*xbar_pow4*one_minus_et2_pow1p5;
    e2dot = prefactor*(e2dot + e2dot_spin);

    return e2dot;
}

/**
 * @brief Compute the time derivative of the mean anomaly for eccentric binaries
 *
 * This function calculates \f$\frac{dl}{dt}\f$, the rate of change of the mean anomaly
 * in the quasi-Keplerian orbital parameterization for eccentric binary systems.
 *
 * @param[in] omega Orbital angular velocity in units of M = 1
 * @param[in] etpow2 Square of eccentricity (unused)
 * @param[in] l Mean anomaly (unused)
 * @param[in] lambda Secular orbital phase contribution (unused)
 * @param[in] s1x, s1y, s1z Spin components of body 1 (unused)
 * @param[in] s2x, s2y, s2z Spin components of body 2 (unused)
 * @param[in] lhatx, lhaty, lhatz Orbital angular momentum unit vector (unused)
 * @param[in] phatx, phaty, phatz Periastron direction (unused)
 * @param[in] k Periastron advance rate
 * @param[in] ak Structure containing parameters (unused)
 *
 * @return Time derivative of the mean anomaly \f$\frac{dl}{dt} = \frac{\omega}{1+k}\f$
 *
 * @see Equation 6d from Phukon et al., arXiv:2504.20543.
 */
REAL8 XLALSimInspiralMeanAnomalyEvolution( REAL8 omega, REAL8 UNUSED etpow2, REAL8 UNUSED l, REAL8 UNUSED lambda,
                                           REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
                                           REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
                                           REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
                                           REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
                                           REAL8 k, expnCoeffsSpinEcc UNUSED *ak )
{
    REAL8 n = omega/(1+k);
    return n;
}

/**
 * @brief Compute the time derivative of the mean orbital phase for eccentric binaries
 *
 * This function calculates \f$\frac{d\lambda}{dt}\f$ in the quasi-Keplerian orbital parameterization
 * for binary orbit.
 *
 * @param[in] omega Orbital angular velocity in units of M = 1
 * @param[in] etpow2 Square of eccentricity (unused)
 * @param[in] l Mean anomaly (unused)
 * @param[in] lambda Secular orbital phase contribution (unused)
 * @param[in] s1x, s1y, s1z Spin components of body 1 (unused)
 * @param[in] s2x, s2y, s2z Spin components of body 2 (unused)
 * @param[in] lhatx, lhaty, lhatz Orbital angular momentum unit vector (unused)
 * @param[in] phatx, phaty, phatz Periastron direction (unused)
 * @param[in] k Periastron advance rate (unused)
 * @param[in] ak Structure containing parameters (unused)
 *
 * @return Time derivative of the mean orbital phase \f$\frac{d\lambda}{dt} = \omega\f$
 *
 * @see Equations 6e from Phukon et al., arXiv:2504.20543.
 */
REAL8 XLALSimInspiralMeanOrbitalPhaseEvolution( REAL8 omega, REAL8 UNUSED etpow2, REAL8 UNUSED l, REAL8 UNUSED lambda,
                                                REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
                                                REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
                                                REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
                                                REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
                                                REAL8 UNUSED k, expnCoeffsSpinEcc UNUSED *ak )
{
    return omega;
}

/**
 * @brief Compute the time derivative of the x-component of unit vector of spin 1
 *
 * This function calculates the time derivative of the x-component of unit vector of spin 1.
 *
 * @param[in] omega Orbital angular velocity (unused)
 * @param[in] etpow2 Square of eccentricity (unused)
 * @param[in] l Mean anomaly (unused)
 * @param[in] lambda Secular orbital phase contribution (unused)
 * @param[in] s1x, s1y, s1z Spin components of body 1 (unused)
 * @param[in] s2x, s2y, s2z Spin components of body 2 (unused)
 * @param[in] lhatx, lhaty, lhatz Orbital angular momentum unit vector (unused)
 * @param[in] phatx, phaty, phatz Periastron direction (unused)
 * @param[in] k Periastron advance rate (unused)
 * @param[in] ak Structure containing pre-computed dynamical quantities (spin-orbit and spin-spin terms), spin order settings and system parameters
 *
 * @return Time derivative of the x-component of spin 1 unit vector
 *
 * @see Equation 7a of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralS1xEvolution( REAL8 UNUSED omega, REAL8 UNUSED etpow2, REAL8 UNUSED l, REAL8 UNUSED lambda,
				                   REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
				                   REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
				                   REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
				                   REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
				                   REAL8 UNUSED k, expnCoeffsSpinEcc *ak )
{
    REAL8 Lmag = ak->LNmag;
    REAL8 mratio = ak->q;
    REAL8 delta1 = 2. + 3./2. * mratio;

    REAL8 norm1 = ak->norm1;
    REAL8 norm2 = ak->norm2;

    REAL8 DistanceCubeInverse = 1./(ak->SemiMinor*ak->SemiMinor*ak->SemiMinor); /* 1./( pow(a_major,3)* pow(1-etpow2, 3./2)) */

    REAL8 *s1Crosss2 = ak->s1crosss2;

    REAL8 *lhatCrosss1=ak->lcrosss1;

    REAL8 spin_orbit_terms = delta1*Lmag*lhatCrosss1[0];

    REAL8 lhatdots1 = ak->s1dotl;
    REAL8 lhatdots2 = ak->s2dotl;

    REAL8 spin_spin_terms =  - 3./2*( norm2*lhatdots2 + mratio*norm1*lhatdots1)*lhatCrosss1[0]  -  norm2*s1Crosss2[0]/2;

    LALSimInspiralSpinOrder spinO = ak->SpinOrder;

    REAL8 s1xdot = 0;

    switch (spinO)
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            s1xdot = spin_spin_terms;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            s1xdot += spin_orbit_terms;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %d\n",
                                    __func__, spinO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    s1xdot = DistanceCubeInverse*s1xdot;

   return s1xdot;
}

/**
 * @brief Compute the time derivative of the y-component of unit vector of spin 1
 *
 * This function calculates the time derivative of the y-component of unit vector of spin 1.
 *
 * @param[in] omega Orbital angular velocity (unused)
 * @param[in] etpow2 Square of eccentricity (unused)
 * @param[in] l Mean anomaly (unused)
 * @param[in] lambda Secular orbital phase contribution (unused)
 * @param[in] s1x, s1y, s1z Spin components of body 1 (unused)
 * @param[in] s2x, s2y, s2z Spin components of body 2 (unused)
 * @param[in] lhatx, lhaty, lhatz Orbital angular momentum unit vector (unused)
 * @param[in] phatx, phaty, phatz Periastron direction (unused)
 * @param[in] k Periastron advance rate (unused)
 * @param[in] ak Structure containing pre-computed dynamical quantities (spin-orbit and spin-spin terms), spin order settings and system parameters
 *
 * @return Time derivative of the y-component of spin 1 unit vector
 *
 * @see Equation 7a of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralS1yEvolution(  REAL8 UNUSED omega, REAL8 UNUSED etpow2, REAL8 UNUSED l, REAL8 UNUSED lambda,
				                    REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
				                    REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
				                    REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
				                    REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
				                    REAL8 UNUSED k, expnCoeffsSpinEcc *ak )
{
    REAL8 Lmag =  ak->LNmag;
    REAL8 mratio = ak->q;
    REAL8 delta1 = 2. + 3./2. * mratio;

    REAL8 norm1 =  ak->norm1;
    REAL8 norm2 =  ak->norm2;

    REAL8 DistanceCubeInverse = 1./(ak->SemiMinor*ak->SemiMinor*ak->SemiMinor);

    REAL8 *s1Crosss2 = ak->s1crosss2;

    REAL8 *lhatCrosss1=ak->lcrosss1;

    REAL8 spin_orbit_terms = delta1*Lmag* lhatCrosss1[1];

    REAL8 lhatdots1 = ak->s1dotl;
    REAL8 lhatdots2 = ak->s2dotl;

    REAL8 spin_spin_terms =  - 3./2*( norm2*lhatdots2 + mratio*norm1*lhatdots1)*lhatCrosss1[1]  - norm2*s1Crosss2[1]/2;

    LALSimInspiralSpinOrder spinO = ak->SpinOrder;

    REAL8 s1ydot = 0;

    switch (spinO)
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            s1ydot = spin_spin_terms;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            s1ydot += spin_orbit_terms;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %d\n",
                                    __func__, spinO );
            XLAL_ERROR(XLAL_EINVAL);
        break;
    }

    s1ydot = DistanceCubeInverse*s1ydot;

    return s1ydot;
}

/**
 * @brief Compute the time derivative of the z-component of unit vector of spin 1
 *
 * This function calculates the time derivative of the z-component of unit vector of spin 1.
 *
 * @param[in] omega Orbital angular velocity (unused)
 * @param[in] etpow2 Square of eccentricity (unused)
 * @param[in] l Mean anomaly (unused)
 * @param[in] lambda Secular orbital phase contribution (unused)
 * @param[in] s1x, s1y, s1z Spin components of body 1 (unused)
 * @param[in] s2x, s2y, s2z Spin components of body 2 (unused)
 * @param[in] lhatx, lhaty, lhatz Orbital angular momentum unit vector (unused)
 * @param[in] phatx, phaty, phatz Periastron direction (unused)
 * @param[in] k Periastron advance rate (unused)
 * @param[in] ak Structure containing pre-computed dynamical quantities (spin-orbit and spin-spin terms), spin order settings and system parameters
 *
 * @return Time derivative of the z-component of spin 1 unit vector
 *
 * @see Equation 7a of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralS1zEvolution(  REAL8 UNUSED omega, REAL8 UNUSED etpow2, REAL8 UNUSED l, REAL8 UNUSED lambda,
				                    REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
				                    REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
				                    REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
				                    REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
				                    REAL8 UNUSED k, expnCoeffsSpinEcc *ak )
{
    REAL8 Lmag =  ak->LNmag;
    REAL8 mratio = ak->q;
    REAL8 delta1 = 2. + 3./2. * mratio;

    REAL8 norm1 = ak->norm1;
    REAL8 norm2 = ak->norm2;

    REAL8 DistanceCubeInverse = 1./(ak->SemiMinor*ak->SemiMinor*ak->SemiMinor);

    REAL8 *s1Crosss2 = ak->s1crosss2;

    REAL8 *lhatCrosss1=ak->lcrosss1;

    REAL8 spin_orbit_terms = delta1*Lmag* lhatCrosss1[2];

    REAL8 lhatdots1 = ak->s1dotl;
    REAL8 lhatdots2 = ak->s2dotl;

    REAL8 spin_spin_terms =  - 3./2*( norm2*lhatdots2 + mratio*norm1*lhatdots1)*lhatCrosss1[2]  - norm2*s1Crosss2[2]/2;

    LALSimInspiralSpinOrder spinO = ak->SpinOrder;

    REAL8 s1zdot = 0;

    switch (spinO)
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            s1zdot = spin_spin_terms;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            s1zdot += spin_orbit_terms;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %d\n",
                                    __func__, spinO );
            XLAL_ERROR(XLAL_EINVAL);
        break;
    }
    s1zdot = DistanceCubeInverse*s1zdot ;

    return s1zdot;
}

/**
 * @brief Computes the time derivative of the x-component of the unit spin vector of body 2
 *
 * This function calculates the time derivative of the x-component of the unit spin vector 2
 *
 * @param[in] omega Orbital angular velocity (unused)
 * @param[in] etpow2 Square of the eccentricity (unused)
 * @param[in] l Mean anomaly (unused)
 * @param[in] lambda Secular orbital phase contribution (unused)
 * @param[in] s1x, s1y, s1z Cartesian components of spin 1 (unused)
 * @param[in] s2x, s2y, s2z Cartesian components of spin 2 (unused)
 * @param[in] lhatx, lhaty, lhatz Unit vector components of orbital angular momentum (unused)
 * @param[in] phatx, phaty, phatz Periastron direction unit vector components (unused)
 * @param[in] k Periastron advance rate (unused)
 * @param[in] ak Structure containing pre-computed dynamic parameters (spin-spin, spin-orbit terms) and static parameters (masses, PN orders, etc.)
 *
 * @return Time derivative of the x-component of the unit spin vector of body 2
 *
 * @see Equation 7b of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralS2xEvolution(  REAL8 UNUSED omega, REAL8 UNUSED etpow2, REAL8 UNUSED l, REAL8 UNUSED lambda,
				                    REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
				                    REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
				                    REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
				                    REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
				                    REAL8 UNUSED k, expnCoeffsSpinEcc *ak )
{
    REAL8 Lmag =  ak->LNmag;
    REAL8 mratio = ak->q;
    REAL8 mratioInverse = 1./mratio;
    REAL8 delta2 = 2. + 3./2. * mratioInverse;

    REAL8 norm1 = ak->norm1;
    REAL8 norm2 = ak->norm2;

    REAL8 DistanceCubeInverse = 1./(ak->SemiMinor*ak->SemiMinor*ak->SemiMinor);

    REAL8 *s1Crosss2= ak->s1crosss2;

    REAL8 *lhatCrosss2=ak->lcrosss2;

    REAL8 spin_orbit_terms = delta2*Lmag* lhatCrosss2[0];

    REAL8 lhatdots1 = ak->s1dotl;
    REAL8 lhatdots2 = ak->s2dotl;

    REAL8 spin_spin_terms =  - 3./2*( norm1*lhatdots1 + mratioInverse*norm2*lhatdots2)*lhatCrosss2[0] + norm1*s1Crosss2[0]/2;

    LALSimInspiralSpinOrder spinO = ak->SpinOrder;

    REAL8 s2xdot = 0;

    switch (spinO)
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            s2xdot = spin_spin_terms;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            s2xdot += spin_orbit_terms;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %d\n",
                                    __func__, spinO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }
    s2xdot = DistanceCubeInverse*s2xdot;

    return s2xdot;
}

/**
 * @brief Computes the time derivative of the y-component of the unit spin vector of body 2
 *
 * This function calculates the time derivative of the y-component of the unit spin vector 2
 *
 * @param[in] omega Orbital angular velocity (unused)
 * @param[in] etpow2 Square of the eccentricity (unused)
 * @param[in] l Mean anomaly (unused)
 * @param[in] lambda Secular orbital phase contribution (unused)
 * @param[in] s1x, s1y, s1z Cartesian components of spin 1 (unused)
 * @param[in] s2x, s2y, s2z Cartesian components of spin 2 (unused)
 * @param[in] lhatx, lhaty, lhatz Unit vector components of orbital angular momentum (unused)
 * @param[in] phatx, phaty, phatz Periastron direction unit vector components (unused)
 * @param[in] k Periastron advance rate (unused)
 * @param[in] ak Structure containing pre-computed dynamic parameters (spin-spin, spin-orbit terms) and static parameters (masses, spin order, etc.)
 *
 * @return Time derivative of the y-component of the unit spin vector of body 2
 *
 * @see Equation 7b of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralS2yEvolution(  REAL8 UNUSED omega, REAL8 UNUSED etpow2, REAL8 UNUSED l, REAL8 UNUSED lambda,
				                    REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
				                    REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
				                    REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
				                    REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
				                    REAL8 UNUSED k, expnCoeffsSpinEcc *ak )
{
    REAL8 Lmag =  ak->LNmag;
    REAL8 mratio = ak->q;
    REAL8 mratioInverse = 1./mratio;
    REAL8 delta2 = 2. + 3./2. * mratioInverse;

    REAL8 norm1 = ak->norm1;
    REAL8 norm2 = ak->norm2;

    REAL8 DistanceCubeInverse = 1./(ak->SemiMinor*ak->SemiMinor*ak->SemiMinor);

    REAL8 *s1Crosss2=ak->s1crosss2;

    REAL8 *lhatCrosss2=ak->lcrosss2;

    REAL8 spin_orbit_terms = delta2*Lmag* lhatCrosss2[1];

    REAL8 lhatdots1 = ak->s1dotl;
    REAL8 lhatdots2 = ak->s2dotl;

    REAL8 spin_spin_terms =  - 3./2*( norm1*lhatdots1 + mratioInverse*norm2*lhatdots2)*lhatCrosss2[1] + norm1*s1Crosss2[1]/2;

    REAL8 s2ydot = 0;

    LALSimInspiralSpinOrder spinO = ak->SpinOrder;

    switch (spinO)
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            s2ydot = spin_spin_terms;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            s2ydot += spin_orbit_terms;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %d\n",
                                    __func__, spinO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }
    s2ydot = DistanceCubeInverse*s2ydot;

    return s2ydot;
}

/**
 * @brief Computes the time derivative of the z-component of the unit spin vector of body 2
 *
 * @param[in] omega Orbital angular velocity (unused)
 * @param[in] etpow2 Square of the eccentricity (unused)
 * @param[in] l Mean anomaly (unused)
 * @param[in] lambda Secular orbital phase contribution (unused)
 * @param[in] s1x, s1y, s1z Cartesian components of spin 1 (unused)
 * @param[in] s2x, s2y, s2z Cartesian components of spin 2 (unused)
 * @param[in] lhatx, lhaty, lhatz Unit vector components of orbital angular momentum (unused)
 * @param[in] phatx, phaty, phatz Periastron direction unit vector components (unused)
 * @param[in] k Periastron advance rate (unused)
 * @param[in] ak Structure containing pre-computed dynamic parameters (spin-spin, spin-orbit terms) and static parameters (masses, spin order, etc.)
 *
 * @return Time derivative of the z-component of the unit spin vector of body 2
 *
 * @see Equation 7b of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralS2zEvolution(  REAL8 UNUSED omega, REAL8 UNUSED etpow2, REAL8 UNUSED l, REAL8 UNUSED lambda,
				                    REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
				                    REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
				                    REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
				                    REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
				                    REAL8 UNUSED k, expnCoeffsSpinEcc *ak )
{
    REAL8 Lmag =  ak->LNmag;
    REAL8 mratio = ak->q;
    REAL8 mratioInverse = 1./mratio;
    REAL8 delta2 = 2. + 3./2. * mratioInverse;

    REAL8 norm1 = ak->norm1;
    REAL8 norm2 = ak->norm2;

    REAL8 DistanceCubeInverse = 1./(ak->SemiMinor*ak->SemiMinor*ak->SemiMinor);

    REAL8 *s1Crosss2= ak->s1crosss2;

    REAL8 *lhatCrosss2=ak->lcrosss2;

    REAL8 spin_orbit_terms = delta2*Lmag* lhatCrosss2[2];

    REAL8 lhatdots1 = ak->s1dotl;
    REAL8 lhatdots2 = ak->s2dotl;

    REAL8 spin_spin_terms =  - 3./2*( norm1*lhatdots1 + mratioInverse*norm2*lhatdots2)*lhatCrosss2[2] + norm1*s1Crosss2[2]/2;

    REAL8 s2zdot = 0;

    LALSimInspiralSpinOrder spinO = ak->SpinOrder;

    switch (spinO)
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            s2zdot = spin_spin_terms;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            s2zdot += spin_orbit_terms;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %d\n",
                                    __func__, spinO );
            XLAL_ERROR(XLAL_EINVAL);
        break;
    }
    s2zdot = DistanceCubeInverse*s2zdot;

    return s2zdot;
}

/**
 * @brief Compute the time derivative of the x-component of the orbital angular momentum unit vector
 *
 * The implementation is determined by conservation of total angular momentum.
 *
 * @param[in] omega Orbital angular velocity (unused)
 * @param[in] etpow2 Square of eccentricity (unused)
 * @param[in] l Mean anomaly (unused)
 * @param[in] lambda Secular orbital phase contribution (unused)
 * @param[in] s1x, s1y, s1z Spin components of body 1 (unused)
 * @param[in] s2x, s2y, s2z Spin components of body 2 (unused)
 * @param[in] lhatx, lhaty, lhatz Orbital angular momentum unit vector (unused)
 * @param[in] phatx, phaty, phatz Periastron direction (unused)
 * @param[in] k Periastron advance rate (unused)
 * @param[in] ak Structure containing pre-computed dynamical quantities (spin-orbit, spin-spin terms), spin order, static system parameters etc.
 *
 * @return Time derivative of the x-component of orbital angular momentum unit vector
 */
REAL8 XLALSimInspiralLhatxEvolution(    REAL8 UNUSED omega, REAL8 UNUSED etpow2, REAL8 UNUSED l, REAL8 UNUSED lambda,
				                        REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
				                        REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
				                        REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
				                        REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
				                        REAL8 UNUSED k, expnCoeffsSpinEcc *ak )
{
    REAL8 Lmag =  ak->LNmag;

    REAL8 norm1 = ak->norm1;
    REAL8 norm2 = ak->norm2;
    /* We pass zeros for the unused arguments of these functions */
    REAL8 lhatxdot = -( norm1*XLALSimInspiralS1xEvolution (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ak)
		       + norm2*XLALSimInspiralS2xEvolution (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ak))/Lmag;
    return lhatxdot;
}

/**
 * @brief Compute the time derivative of the y-component of the orbital angular momentum unit vector
 *
 * The implementation is determined by conservation of total angular momentum.
 *
 * @param[in] omega Orbital angular velocity (unused)
 * @param[in] etpow2 Square of eccentricity (unused)
 * @param[in] l Mean anomaly (unused)
 * @param[in] lambda Secular orbital phase contribution (unused)
 * @param[in] s1x, s1y, s1z Spin components of body 1 (unused)
 * @param[in] s2x, s2y, s2z Spin components of body 2 (unused)
 * @param[in] lhatx, lhaty, lhatz Orbital angular momentum unit vector (unused)
 * @param[in] phatx, phaty, phatz Periastron direction (unused)
 * @param[in] k Periastron advance rate (unused)
 * @param[in] ak Structure containing pre-computed dynamical quantities (spin-orbit, spin-spin terms), spin order, static system parameters etc.
 *
 * @return Time derivative of the y-component of orbital angular momentum unit vector
 */
REAL8 XLALSimInspiralLhatyEvolution(    REAL8 UNUSED omega, REAL8 UNUSED etpow2, REAL8 UNUSED l, REAL8 UNUSED lambda,
				                        REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
				                        REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
				                        REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
				                        REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
				                        REAL8 UNUSED k, expnCoeffsSpinEcc *ak )
{
    REAL8 Lmag =  ak->LNmag;

    REAL8 norm1 = ak->norm1;
    REAL8 norm2 = ak->norm2;
    /* We pass zeros for the unused arguments of these functions */
    REAL8 lhatydot = -( norm1*XLALSimInspiralS1yEvolution (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ak)
		       + norm2*XLALSimInspiralS2yEvolution (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ak))/Lmag;
    return lhatydot;
}

/**
 * @brief Compute the time derivative of the z-component of the orbital angular momentum unit vector
 *
 * The implementation is determined by conservation of total angular momentum.
 *
 * @param[in] omega Orbital angular velocity (unused)
 * @param[in] etpow2 Square of eccentricity (unused)
 * @param[in] l Mean anomaly (unused)
 * @param[in] lambda Secular orbital phase contribution (unused)
 * @param[in] s1x, s1y, s1z Spin components of body 1 (unused)
 * @param[in] s2x, s2y, s2z Spin components of body 2 (unused)
 * @param[in] lhatx, lhaty, lhatz Orbital angular momentum unit vector (unused)
 * @param[in] phatx, phaty, phatz Periastron direction (unused)
 * @param[in] k Periastron advance rate (unused)
 * @param[in] ak Structure containing pre-computed dynamical quantities (spin-orbit, spin-spin terms), spin order, static system parameters etc.
 *
 * @return Time derivative of the z-component of orbital angular momentum unit vector
 */
REAL8 XLALSimInspiralLhatzEvolution(    REAL8 UNUSED omega, REAL8 UNUSED etpow2, REAL8 UNUSED l, REAL8 UNUSED lambda,
				                        REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
				                        REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
				                        REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
				                        REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
				                        REAL8 UNUSED k, expnCoeffsSpinEcc *ak )
{
    REAL8 Lmag =  ak->LNmag;

    REAL8 norm1 = ak->norm1;
    REAL8 norm2 = ak->norm2;
    /* We pass zeros for the unused arguments of these functions */
    REAL8 lhatzdot = -( norm1*XLALSimInspiralS1zEvolution (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ak)
		       + norm2*XLALSimInspiralS2zEvolution (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ak))/Lmag;
    return lhatzdot;
}

/**
 * @brief Compute the time derivative of the x-component of the periastron-direction unit vector
 *
 * @param[in] omega Orbital angular frequency in units of M=1
 * @param[in] etSq Square of eccentricity (unused)
 * @param[in] l Mean anomaly (unused)
 * @param[in] lambda Secular orbital phase contribution (unused)
 * @param[in] s1x, s1y, s1z Components of spin 1 (unused)
 * @param[in] s2x, s2y, s2z Components of spin 2 (unused)
 * @param[in] lhatx, lhaty, lhatz Orbital angular momentum unit vector (unused)
 * @param[in] phatx, phaty, phatz Periastron direction unit vector (unused)
 * @param[in] k Periastron advance rate
 * @param[in] ak Structure containing pre-computed dynamical quantities (spin-orbit, spin-spin, orbit/spin-periastron terms), spin order, static system parameters etc.
 *
 * @return Time derivative of the x-component of the periastron-direction unit vector
 *
 * @see Equation 10 of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralPeriastronLinexEvolution(  REAL8 omega, REAL8 UNUSED etSq, REAL8 UNUSED l, REAL8 UNUSED lambda,
					                            REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
					                            REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
					                            REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
					                            REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
					                            REAL8 k, expnCoeffsSpinEcc *ak)
{
    REAL8 Lmag = ak->LNmag;
    REAL8 DistanceCubeInverse = 1./ (ak->SemiMinor * ak->SemiMinor * ak->SemiMinor);

    REAL8 mratio = ak->q;
    REAL8 mratioInverse = 1./mratio;

    REAL8 norm1 = ak->norm1;
    REAL8 norm2 = ak->norm2;

    /*
     * Eqs. 8b and 8c from Phukon et al., 2504.20543 for delta1 and delta2
     * These terms are available upto 2PN order
     *
     */
    REAL8 delta1 = 0;
    REAL8 delta2 = 0;

    REAL8 s1dotl = ak->s1dotl;
    REAL8 s2dotl = ak->s2dotl;

    LALSimInspiralSpinOrder spinO = ak->SpinOrder;

    switch (spinO)
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            delta1 = -3./2*(norm2*s2dotl + norm1*mratio*s1dotl)/Lmag;
            delta2 = -3./2*(norm1*s1dotl + norm2*mratioInverse*s2dotl)/Lmag;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            delta1 += 2 + 3./2*mratio;
            delta2 += 2 + 3./2*mratioInverse;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %d\n",
                                    __func__, spinO );
            XLAL_ERROR(XLAL_EINVAL);
        break;
    }

    delta1 *= DistanceCubeInverse*norm1;
    delta2 *= DistanceCubeInverse*norm2;

    /* Equation 8a of Phukon et al., 2504.20543 for omegaP */
    REAL8 omegaPdotl = delta1*s1dotl + delta2*s2dotl;

    REAL8 *s1Crossphat=ak->s1crossPhat;

    REAL8 *s2Crossphat=ak->s2crossPhat;

    REAL8 *lcrossphat=ak->lcrossPhat;

    REAL8 n = omega/(1+k);
    REAL8 Phatxdot = (delta1*s1Crossphat[0] + delta2*s2Crossphat[0]) -
                        (omegaPdotl - k*n)*lcrossphat[0];

    return Phatxdot;
}

/**
 * @brief Compute the time derivative of the y-component of the periastron-direction unit vector
 *
 * @param[in] omega Orbital angular frequency in units of M=1
 * @param[in] etSq Square of eccentricity (unused)
 * @param[in] l Mean anomaly (unused)
 * @param[in] lambda Secular orbital phase contribution (unused)
 * @param[in] s1x, s1y, s1z Components of spin 1 (unused)
 * @param[in] s2x, s2y, s2z Components of spin 2 (unused)
 * @param[in] lhatx, lhaty, lhatz Orbital angular momentum unit vector (unused)
 * @param[in] phatx, phaty, phatz Periastron direction unit vector (unused)
 * @param[in] k Periastron advance rate
 * @param[in] ak Structure containing pre-computed dynamical quantities (spin-orbit, spin-spin, orbit/spin-periastron terms), spin order, static system parameters etc.
 *
 * @return Time derivative of the y-component of the periastron-direction unit vector
 *
 * @see Equation 10 of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralPeriastronLineyEvolution(  REAL8 omega, REAL8 UNUSED etSq, REAL8 UNUSED l, REAL8 UNUSED lambda,
					                            REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
					                            REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
					                            REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
					                            REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
					                            REAL8 k, expnCoeffsSpinEcc *ak)
{
    REAL8 Lmag = ak->LNmag;
    REAL8 DistanceCubeInverse = 1./ (ak->SemiMinor * ak->SemiMinor * ak->SemiMinor);

    REAL8 mratio = ak->q;
    REAL8 mratioInverse = 1./mratio;

    REAL8 norm1 = ak->norm1;
    REAL8 norm2 = ak->norm2;

    /*
     * Eqs. 8b and 8c from Phukon et al., 2504.20543 for delta1 and delta2
     * These terms are available upto 2PN order
     *
     */
    REAL8 delta1 = 0;
    REAL8 delta2 = 0;

    REAL8 s1dotl = ak->s1dotl;
    REAL8 s2dotl = ak->s2dotl;

    LALSimInspiralSpinOrder spinO = ak->SpinOrder;

    switch (spinO)
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            delta1 = -3./2*(norm2*s2dotl + norm1*mratio*s1dotl)/Lmag;
            delta2 = -3./2*(norm1*s1dotl + norm2*mratioInverse*s2dotl)/Lmag;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            delta1 += 2 + 3./2*mratio;
            delta2 += 2 + 3./2*mratioInverse;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %d\n",
                                    __func__, spinO );
            XLAL_ERROR(XLAL_EINVAL);
        break;
    }

    delta1 *= DistanceCubeInverse*norm1;
    delta2 *= DistanceCubeInverse*norm2;

    /* Equation 8a of Phukon et al., 2504.20543 for omegaP */
    REAL8 omegaPdotl = delta1*s1dotl + delta2*s2dotl;

    REAL8 *s1Crossphat=ak->s1crossPhat;

    REAL8 *s2Crossphat=ak->s2crossPhat;

    REAL8 *lcrossphat=ak->lcrossPhat;

    REAL8 n = omega/(1+k);
    REAL8 Phatydot = (delta1*s1Crossphat[1] + delta2*s2Crossphat[1]) -
                        (omegaPdotl - k*n)*lcrossphat[1];
    return Phatydot;
}

/**
 * @brief Compute the time derivative of the z-component of the periastron-direction unit vector
 *
 * @param[in] omega Orbital angular frequency in units of M=1
 * @param[in] etSq Square of eccentricity (unused)
 * @param[in] l Mean anomaly (unused)
 * @param[in] lambda Secular orbital phase contribution (unused)
 * @param[in] s1x, s1y, s1z Components of spin 1 (unused)
 * @param[in] s2x, s2y, s2z Components of spin 2 (unused)
 * @param[in] lhatx, lhaty, lhatz Orbital angular momentum unit vector (unused)
 * @param[in] phatx, phaty, phatz Periastron direction unit vector (unused)
 * @param[in] k Periastron advance rate
 * @param[in] ak Structure containing pre-computed dynamical quantities (spin-orbit, spin-spin, orbit/spin-periastron terms), spin order, static system parameters etc.
 *
 * @return Time derivative of the z-component of the periastron-direction unit vector
 *
 * @see Equation 10 of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralPeriastronLinezEvolution(  REAL8 omega, REAL8 UNUSED etSq, REAL8 UNUSED l, REAL8 UNUSED lambda,
					                            REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
					                            REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
					                            REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
					                            REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
					                            REAL8 k, expnCoeffsSpinEcc *ak)
{
    REAL8 Lmag = ak->LNmag;
    REAL8 DistanceCubeInverse = 1./ (ak->SemiMinor * ak->SemiMinor * ak->SemiMinor);

    REAL8 mratio = ak->q;
    REAL8 mratioInverse = 1./mratio;

    REAL8 norm1 = ak->norm1;
    REAL8 norm2 = ak->norm2;

    /*
     * Eqs. 8b and 8c from Phukon et al., 2504.20543 for delta1 and delta2
     * These terms are available upto 2PN order
     *
     */
    REAL8 delta1 = 0;
    REAL8 delta2 = 0;

    REAL8 s1dotl = ak->s1dotl;
    REAL8 s2dotl = ak->s2dotl;

    LALSimInspiralSpinOrder spinO = ak->SpinOrder;

    switch (spinO)
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            delta1 = -3./2*(norm2*s2dotl + norm1*mratio*s1dotl)/Lmag;
            delta2 = -3./2*(norm1*s1dotl + norm2*mratioInverse*s2dotl)/Lmag;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            delta1 += 2 + 3./2*mratio;
            delta2 += 2 + 3./2*mratioInverse;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %d\n",
                                    __func__, spinO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    delta1 *= DistanceCubeInverse*norm1;
    delta2 *= DistanceCubeInverse*norm2;

    /* Equation 8a of Phukon et al., 2504.20543 for omegaP */
    REAL8 omegaPdotl = delta1*s1dotl + delta2*s2dotl;

    REAL8 *s1Crossphat = ak->s1crossPhat;

    REAL8 *s2Crossphat=ak->s2crossPhat;

    REAL8 *lcrossphat=ak->lcrossPhat;

    REAL8 n = omega/(1+k);
    REAL8 Phatzdot = (delta1*s1Crossphat[2] + delta2*s2Crossphat[2]) -
                        (omegaPdotl - k*n)*lcrossphat[2];
    return Phatzdot;
}

/**
 * @brief Compute the time derivative of the periastron precession
 *
 * This function calculates the rate of change of the periastron advance parameter k in
 * the quasi-Keplerian orbital parameterization.
 *
 * @param[in] omega Orbital angular velocity in units where total mass M = 1
 * @param[in] etSq Square of the eccentricity
 * @param[in] l Mean anomaly (unused)
 * @param[in] lambda Secular orbital phase contribution (unused)
 * @param[in] s1x, s1y, s1z Spin components of body 1 (unused - computed from ak structure)
 * @param[in] s2x, s2y, s2z Spin components of body 2 (unused - computed from ak structure)
 * @param[in] lhatx, lhaty, lhatz Orbital angular momentum unit vector (unused)
 * @param[in] phatx, phaty, phatz Periastron direction (unused)
 * @param[in] k Periastron advance rate (unused)
 * @param[in] ak Structure containing pre-computed dynamical/static quantities, spin order etc.
 *
 * @return Time derivative of the periastron precession
 *
 * @see Equations 6c and A3a-A3e of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralPeriastronPrecessionEvolution( REAL8 omega, REAL8 etSq, REAL8 UNUSED l, REAL8 UNUSED lambda,
						                            REAL8 UNUSED s1x, REAL8 UNUSED s1y, REAL8 UNUSED s1z,
						                            REAL8 UNUSED s2x, REAL8 UNUSED s2y, REAL8 UNUSED s2z,
						                            REAL8 UNUSED lhatx, REAL8 UNUSED lhaty, REAL8 UNUSED lhatz,
						                            REAL8 UNUSED phatx, REAL8 UNUSED phaty, REAL8 UNUSED phatz,
						                            REAL8 UNUSED k, expnCoeffsSpinEcc *ak)
{
    REAL8 LAL_PI_POW2  = LAL_PI*LAL_PI;

    REAL8 one_minus_etSq = 1 - etSq;
    REAL8 one_minus_etSq_pow5 = pow(one_minus_etSq, 5);
    REAL8 sqrt_one_minus_etSq = sqrt(one_minus_etSq);
    REAL8 one_minus_etSq_pow3by2 = (one_minus_etSq*sqrt_one_minus_etSq);

    REAL8 v = cbrt( omega);
    REAL8 x = v*v;
    REAL8 xbar = x/one_minus_etSq;
    REAL8 sqrt_xbar = v/sqrt_one_minus_etSq;
    REAL8 xbar_pow3by2 = xbar*sqrt_xbar;
    REAL8 xbar_pow4 = xbar_pow3by2*xbar_pow3by2*xbar;

    REAL8 eta = ak->eta;
    /* Equation A3a from Phukon et al., 2504.20543 */
    REAL8 K_1pn = 192./5 + 168 * etSq/5;
    /* Non-spinning terms from Equation A3c of Phukon et al., 2504.20543 */
    REAL8 K_2pn = (9124./35 - 1424./5 * eta) + etSq * ( (28512./35 - 3804./5 * eta ) +
                                               etSq * (10314./35-1017./5 * eta));
    /* Equation A3e from Phukon et al., 2504.20543 */
    REAL8 K_3pn =   (232082./189 + eta*((-131366./21 + 738./5 * LAL_PI_POW2 ) + 13312./15 * eta)) +
                    etSq * ( 2842199./630 + eta * ((-1659934./105 + 1271./10 * LAL_PI_POW2 )  + 29879./5 * eta ) +
                    etSq * ( 1640713./252 + eta * (( -1304524./105 + 5371./320 * LAL_PI_POW2 ) + 54133./10 * eta ) +
		            etSq * ( 1850407./1680 - eta*( 388799./280 - 19573./30 * eta )))) +
                    sqrt_one_minus_etSq * ( 672 - 1344./5 * eta +
                            etSq* ( ( 2436 - 4872./ 5 * eta ) +
				            etSq * (672 - 1344./5 * eta)));

    EnhancementFunction enhancementFunc = ak->enhancementFunc;
    LALSimInspiralSpinEccPolyEccOrder EccOrder = ak->EccOrder;

    /* Equation A3d from Phukon et al., 2504.20543 */
    REAL8 K_hered = 768./5*LAL_PI*XLALSimInspiralSpinEccPhik(etSq, enhancementFunc, EccOrder)*one_minus_etSq_pow5;

    REAL8 k_dot = 0;

    INT4 phaseO = ak->PhaseOrder;

    switch (phaseO)
    {
        case(-1):
        case(6):
	        k_dot = K_3pn;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(5):
	        k_dot *= sqrt_xbar;
            k_dot += K_hered;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(4):
	        k_dot *= sqrt_xbar;
            k_dot += K_2pn;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(3):
        case(2):
	        k_dot *= xbar;
            k_dot += K_1pn;
	        k_dot *= xbar;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(1):
        case(0):
           break;
        default:
           XLALPrintError("XLAL Error - %s: - %d is an invalid phase order\n", __func__, phaseO );
           XLAL_ERROR(XLAL_EINVAL);
           break;

    }

    /* Prefactor of Equation 6c of Phukon et al., 2504.20543 */
    REAL8 prefactor = eta*xbar_pow4*one_minus_etSq_pow3by2;

    REAL8 gamma1 = ak->gamma1;

    REAL8 beta43 = ak->beta43;

    /* 1.5PN spin contribution;  Equation A3b of Phukon et al., 2504.20543*/
    REAL8 kdot_spin_15PN = -(96 + 84*etSq) * beta43/5;

    /* 2PN spin contribution; last term from Equation A3c of Phukon et al., 2504.20543 */
    REAL8 kdot_spin_2PN = (192 + 168 * etSq)/5* gamma1;

    LALSimInspiralSpinOrder spinO = ak->SpinOrder;

    REAL8 kdot_spin = 0;
    switch (spinO)
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            kdot_spin = sqrt_xbar*kdot_spin_2PN;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            kdot_spin += kdot_spin_15PN;
	        kdot_spin *= xbar_pow3by2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error in function - %s: Invalid spin PN order %d\n",__func__, spinO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }
    k_dot = prefactor * (k_dot + kdot_spin);
    return k_dot;
}

/**
 * @brief Compute the periastron advance parameter k using binary parameters
 *
 * This function computes the instantaneous value of the periastron advance parameter.
 * The calculation includes post-Newtonian corrections up to 3PN order for non-spinning
 * effects and up to 2PN order for spin effects. The result is used to set initial
 * conditions for the orbital evolution.
 *
 * The calculation follows Equations 12, 13a, and 13b from Phukon et al., arXiv:2504.20543.
 *
 * @param[in] xbar PN parameter (\f$\bar{x} = x/(1 - e_t^2)\f$)
 * @param[in] etSq Square of the eccentricity (\f$e_t^2\f$)
 * @param[in] ak Structure containing static parameters (masses, spins, etc.), precomputed spin terms
 * @param[in] spinO Twice the PN order of spin effects
 * @param[in] phaseO Twice the PN order of non-spinning phase terms
 *
 * @return Periastron advance parameter
 *
 *@see Equations 12, 13a, and 13b from Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralPeriastronPrecession( REAL8 xbar, REAL8 etSq, expnCoeffsSpinEcc *ak,
                                            LALSimInspiralSpinOrder spinO, INT4 phaseO )
{
    REAL8 LAL_PI_POW2  = LAL_PI*LAL_PI;

    REAL8 eta = ak->eta;

    /* 1PN term from Equation 12, Phukon et al., 2504.20543*/
    REAL8 k_1pn = 3.;
    /* 2PN non-spinning terms from Equation 12, Phukon et al., 2504.20543 */
    REAL8 k_2pn = (54. - 28.*eta + etSq*(51. - 26.*eta))/4.;
    /* 3PN terms from Equation 12, Phukon et al., 2504.20543 */
    REAL8 k_3pn = (210. + eta*(-625. + 123.*LAL_PI_POW2/8. + 28.*eta) +
                    etSq*(567. + eta*(-816. + 123.*LAL_PI_POW2/32. + 160.*eta) +
			        etSq*(78. + eta*(-55. + 65.*eta/2.))) +
		            sqrt(1. - etSq)*(60. - 24.*eta + etSq*(120. - 48.*eta)))/4.;

    REAL8 k = 0;

    switch (phaseO)
    {
        case(-1):
        case(6):
           k = k_3pn;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(5):
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(4):
	        k *= xbar;
            k += k_2pn;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(3):
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(2):
            k *= xbar;
            k += k_1pn;
	        k *= xbar;
        case(1):
        case(0):
            break;
        default:
           XLALPrintError("XLAL Error - %s: - %d is an invalid phase order\n", __func__, phaseO );
           XLAL_ERROR(XLAL_EINVAL);
           break;
    }
    /* The 1.5 PN spin term from Eqs. 12, 13a Phukon et al., 2504.20543 */
    REAL8 k_spin_15PN = -1*ak->beta43;
    /* 2PN spin term is from Eqs. 12, 13b Phukon et al., 2504.20543 */
    REAL8 k_spin_2PN = 1.5*ak->gamma1;

    REAL8 k_spin = 0;

    REAL8 sqrt_xbar = sqrt(xbar);

    switch (spinO)
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
	        k_spin = xbar*k_spin_2PN;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
	        k_spin += sqrt_xbar*k_spin_15PN;
	        k_spin *= xbar;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error in function - %s: Invalid spin PN order %d\n",__func__, spinO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    return k + k_spin;
}

/**
 * @brief Wrapper to `ref XLALSimInspiralPeriastronPrecession` for computation of periastron precession
 *
 * This function returns the periastron advance parameter k for a given GW frequency, masses, spins, and eccentricity
 * and fill expnCoeffsSpinEcc structure with binary parameters and spin-orbir/spin terms.
 *
 * @param[in] fgw Gravitational-wave frequency
 * @param[in] m1_SI Mass of body 1 in SI units
 * @param[in] m2_SI Mass of body 2 in SI units
 * @param[in] s1x, s1y, s1z Cartesian components of dimensionless spin of body 1
 * @param[in] s2x, s2y, s2z Cartesian components of dimensionless spin of body 2
 * @param[in] lnhatx, lnhaty, lnhatz Unit vector components of orbital angular momentum
 * @param[in] eccentricity Orbital eccentricity
 * @param[in] SpinOrder Spin order to include (twice of spining PN order)
 * @param[in] PhaseOrder Phase order to include (twice of non-spinning PN order)
 *
 * @return Periastron precession rate k (dimensionless)
 *
 */
REAL8 XLALSimInspiralSpinTaylorEccentricComputePeriastronPrecession( REAL8 fgw, REAL8 m1_SI, REAL8 m2_SI,
			REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z,
		        REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz, REAL8 eccentricity,
			LALSimInspiralSpinOrder SpinOrder, INT4 PhaseOrder)
{
    expnCoeffsSpinEcc ak;
    ak.q = m2_SI/m1_SI;
    ak.m1 = 1/(1 + ak.q);
    ak.m2 = ak.q/(1 + ak.q);
    ak.m = ak.m1 + ak.m2;
    ak.eta = ak.m1 * ak.m2 / (ak.m*ak.m);
    ak.mt = (m1_SI + m2_SI)/LAL_MSUN_SI*LAL_MTSUN_SI;
    ak.norm1 = ak.m1 * ak.m1;
    ak.norm2 = ak.m2 * ak.m2;

    REAL8 v = cbrt(LAL_PI*fgw*ak.mt);
    REAL8 x = v*v;
    REAL8 etSq = eccentricity*eccentricity;
    REAL8 xbar = x/(1 - etSq);

    ak.beta43 = XLALSimInspiralSpinEccBeta(4, 3, s1x, s1y, s1z, s2x,
                                        s2y, s2z, lnhatx, lnhaty, lnhatz, &ak);
    ak.gamma1 = XLALSimInspiralSpinEccGamma1(s1x, s1y, s1z, s2x, s2y,  s2z,
                                        lnhatx, lnhaty, lnhatz, &ak );

    REAL8 k_value = XLALSimInspiralPeriastronPrecession(xbar, etSq, &ak, SpinOrder, PhaseOrder);
    return k_value;
}

/**
 * @brief Compute the 1.5PN spin-orbit \f$\beta\f$ factor
 *
 * @param[in] a, b First and second beta function coefficients
 * @param[in] s1x, s1y, s1z Cartesian components of the dimensionless spin of body 1
 * @param[in] s2x, s2y, s2z Cartesian components of the dimensionless spin of body 2
 * @param[in] lhatx, lhaty, lhatz Unit vector components of orbital angular momentum
 * @param[in] ak Structure containing mass ratios and normalization factors
 *
 * @return The 1.5PN spin-orbit \f$\beta\f$ term
 *
 * @see Equation 13a of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccBeta (  REAL8 a, REAL8 b,
        REAL8 s1x, REAL8 s1y, REAL8 s1z,
        REAL8 s2x, REAL8 s2y, REAL8 s2z,
        REAL8 lhatx, REAL8 lhaty, REAL8 lhatz,
        expnCoeffsSpinEcc *ak )

{
    REAL8 norm1 = ak->norm1;
    REAL8 norm2 = ak->norm2;

    REAL8  betax, betay, betaz;

    betax = lhatx * ( a * ( norm1 * s1x + norm2 * s2x ) + b *  ( ak->q * norm1 * s1x + norm2 * s2x / ak->q ));
    betay = lhaty * ( a * ( norm1 * s1y + norm2 * s2y ) + b *  ( ak->q * norm1 * s1y + norm2 * s2y / ak->q ));
    betaz = lhatz * ( a * ( norm1 * s1z + norm2 * s2z ) + b *  ( ak->q * norm1 * s1z + norm2 * s2z / ak->q ));

   return betax + betay + betaz ;
}

/**
 * @brief Compute the 2PN \f$\gamma_1\f$ factor for spin-spin coupling term
 *
 * This function calculates the gamma1 coefficient for 2PN spin-spin term.
 *
 * @param[in] s1x, s1y, s1z Cartesian components of the dimensionless spin of body 1
 * @param[in] s2x, s2y, s2z Cartesian components of the dimensionless spin of body 2
 * @param[in] lhatx, lhaty, lhatz Cartesian components of the unit vector in the direction of orbital angular momentum
 * @param[in] ak Structure containing static parameters (masses, etc.)
 *
 * @return The value of the \f$\gamma_1\f$ spin-spin coupling coefficient
 *
 * @see This corresponds to Equation 13b of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccGamma1 (
        REAL8 s1x, REAL8 s1y, REAL8 s1z,
        REAL8 s2x, REAL8 s2y, REAL8 s2z,
        REAL8 lhatx, REAL8 lhaty, REAL8 lhatz,
        expnCoeffsSpinEcc *ak )

{
    REAL8 mass1 = ak->m1;
    REAL8 mass2 = ak->m2;

    REAL8 s1x_reduced = s1x*mass1;
    REAL8 s1y_reduced = s1y*mass1;
    REAL8 s1z_reduced = s1z*mass1;

    REAL8 s2x_reduced = s2x*mass2;
    REAL8 s2y_reduced = s2y*mass2;
    REAL8 s2z_reduced = s2z*mass2;

    REAL8 sx_reduced = s1x_reduced + s2x_reduced;
    REAL8 sy_reduced = s1y_reduced + s2y_reduced;
    REAL8 sz_reduced = s1z_reduced + s2z_reduced;

    REAL8 s_reduced_NormSq = sx_reduced*sx_reduced + sy_reduced*sy_reduced + sz_reduced*sz_reduced;

    REAL8 lhatdotReducedTotalSpin= lhatx*sx_reduced + lhaty*sy_reduced + lhatz*sz_reduced;

    REAL8 Gamma1 = 0.5*(3*lhatdotReducedTotalSpin*lhatdotReducedTotalSpin - s_reduced_NormSq);
    return Gamma1;
}

/**
 * @brief Newtonian orbital angular momentum magnitude in total mass = 1 units
 *
 * @param[in] omega Orbital frequency
 * @param[in] etSq Square of eccentricity
 * @param[in] ak Structure containing static parameters (component masses, mass ratio etc.)
 *
 * @return Magnitude of Newtonian orbital angular momentum in total mass = 1 units
 */
REAL8 XLALSimInspiralSpinEccLmagN ( REAL8 omega, REAL8 etSq, expnCoeffsSpinEcc *ak  )
{
    REAL8 v = cbrt(omega*ak->m);
    REAL8 LN = ak->eta*sqrt(1-etSq)/v;
    return LN;
}

/**
 * @brief Compute the Newtonian semimajor axisin total mass = 1 units.
 *
 * @param[in] omega Orbital angular velocity
 * @param[in] ak Structure containing static parameters
 *
 * @return Newtonian semimajor axis in total mass = 1 units
 */
REAL8 XLALSimInspiralSpinEccSemiMajorAxis (REAL8 omega, expnCoeffsSpinEcc *ak)
{
    REAL8 mtotal = ak->m;
    REAL8 x = pow( mtotal*omega, 2./3);
    return 1./x;
}

/**
 * @brief Compute various dynamical quantities for spin-precessing eccentric binary evolution
 *
 * This function calculates various cross/dot products of spins, orbital angular momentum, and periastron direction,
 * consines of angles between in-plane spins and periastron, the magnitude of the orbital angular momentum,
 * the semimajor and semiminor axes, the leading flux enhancement function, and the gamma1 `ref XLALSimInspiralSpinEccGamma1` and
 * beta43 `ref XLALSimInspiralSpinEccBeta`(4,3,...) spin terms. The function populates the designated fields in the expnCoeffsSpinEcc structure.
 *
 * @param[out] ak         Pointer to expnCoeffsSpinEcc structure where computed quantities are stored
 * @param[in]  omega      Orbital angular frequency
 * @param[in]  etSq       Square of orbital eccentricity (\f$e^2\f$)
 * @param[in]  s1x,s1y,s1z Cartesian components of dimensionless spin vector for body 1
 * @param[in]  s2x,s2y,s2z Cartesian components of dimensionless spin vector for body 2
 * @param[in]  lnhatx,lnhaty,lnhatz Unit vector components of orbital angular momentum
 * @param[in]  phatx,phaty,phatz Unit vector components of periastron direction
 *
 * @return XLAL_SUCCESS on successful completion
 *
 * @note The function handles the special case where in-plane spin magnitudes are zero to avoid
 *       division by zero when computing the cosine terms.
 *
 */
static int XLALSimInspiralSpinEccDynamicalExpnCoeffsSpinEccVariables( expnCoeffsSpinEcc *ak, REAL8 omega, REAL8 etSq,
                                                        REAL8 s1x, REAL8 s1y, REAL8 s1z,
                                                        REAL8 s2x, REAL8 s2y, REAL8 s2z,
                                                        REAL8 lnhatx, REAL8 lnhaty, REAL8 lnhatz,
                                                        REAL8 phatx, REAL8 phaty, REAL8 phatz)
{
    REAL8 *s1Crosss2=NULL;
    XLALSimInspiralVectorCrossProduct(&s1Crosss2, s1x, s1y, s1z, s2x, s2y, s2z);

    ak->s1crosss2[0] = s1Crosss2[0];
    ak->s1crosss2[1] = s1Crosss2[1];
    ak->s1crosss2[2] = s1Crosss2[2];

    REAL8 *lcrosss1=NULL;
    XLALSimInspiralVectorCrossProduct(&lcrosss1, lnhatx, lnhaty, lnhatz, s1x, s1y, s1z);
    ak->lcrosss1[0] = lcrosss1[0];
    ak->lcrosss1[1] = lcrosss1[1];
    ak->lcrosss1[2] = lcrosss1[2];

    REAL8 *lcrosss2=NULL;
    XLALSimInspiralVectorCrossProduct(&lcrosss2, lnhatx, lnhaty, lnhatz, s2x, s2y, s2z);
    ak->lcrosss2[0] = lcrosss2[0];
    ak->lcrosss2[1] = lcrosss2[1];
    ak->lcrosss2[2] = lcrosss2[2];

    REAL8 s1dl = s1x*lnhatx + s1y*lnhaty + s1z*lnhatz;
    REAL8 s2dl = s2x*lnhatx + s2y*lnhaty + s2z*lnhatz;

    ak->s1dotl = s1dl;
    ak->s2dotl = s2dl;

    REAL8 *s1Crossphat=NULL;
    XLALSimInspiralVectorCrossProduct(&s1Crossphat, s1x, s1y, s1z, phatx, phaty, phatz);

    ak->s1crossPhat[0] = s1Crossphat[0];
    ak->s1crossPhat[1] = s1Crossphat[1];
    ak->s1crossPhat[2] = s1Crossphat[2];

    REAL8 *s2Crossphat=NULL;
    XLALSimInspiralVectorCrossProduct(&s2Crossphat, s2x, s2y, s2z, phatx, phaty, phatz);

    ak->s2crossPhat[0] = s2Crossphat[0];
    ak->s2crossPhat[1] = s2Crossphat[1];
    ak->s2crossPhat[2] = s2Crossphat[2];


    REAL8 *lcrossphat=NULL;
    XLALSimInspiralVectorCrossProduct(&lcrossphat, lnhatx, lnhaty, lnhatz, phatx, phaty, phatz);

    ak->lcrossPhat[0] = lcrossphat[0];
    ak->lcrossPhat[1] = lcrossphat[1];
    ak->lcrossPhat[2] = lcrossphat[2];

    /* Project spins into the orbital plane */
    REAL8 s1inplanex = s1x - s1dl*lnhatx;
    REAL8 s1inplaney = s1y - s1dl*lnhaty;
    REAL8 s1inplanez = s1z - s1dl*lnhatz;

    REAL8 s2inplanex = s2x - s2dl*lnhatx;
    REAL8 s2inplaney = s2y - s2dl*lnhaty;
    REAL8 s2inplanez = s2z - s2dl*lnhatz;

    REAL8 cpsi1, cpsi2, cpsi_ReducedTotalSpin;

    REAL8 s1inplane_mag = sqrt(s1inplanex*s1inplanex + s1inplaney*s1inplaney + s1inplanez*s1inplanez);
    REAL8 s2inplane_mag = sqrt(s2inplanex*s2inplanex + s2inplaney*s2inplaney + s2inplanez*s2inplanez);

    /**
     * Since cpsi1 and cpsi2 (and cpsi_ReducedTotalSpin below) are multiplied by the squares of the in-plane
     * spin magnitudes when the appear in the evolution equations, we set them to zero when the in-plane spin
     *  magnitude is zero, to avoid division by zero
     *
     */
    if (s1inplane_mag > 0) {
      cpsi1 = (s1inplanex*phatx + s1inplaney*phaty + s1inplanez*phatz)/s1inplane_mag;
    } else {
      cpsi1 = 0;
    }
    if (s2inplane_mag > 0) {
      cpsi2 = (s2inplanex*phatx + s2inplaney*phaty + s2inplanez*phatz)/s2inplane_mag;
    } else {
      cpsi2 = 0;
    }

    REAL8 sinplane_x = ak->m1*s1inplanex + ak->m2*s2inplanex;
    REAL8 sinplane_y = ak->m1*s1inplaney + ak->m2*s2inplaney;
    REAL8 sinplane_z = ak->m1*s1inplanez + ak->m2*s2inplanez;
    REAL8 sinplane_mag = sqrt(sinplane_x*sinplane_x + sinplane_y*sinplane_y + sinplane_z*sinplane_z);

    if (sinplane_mag > 0) {
      cpsi_ReducedTotalSpin = (sinplane_x*phatx + sinplane_y*phaty + sinplane_z*phatz)/sinplane_mag;
    } else {
      cpsi_ReducedTotalSpin = 0;
    }

    ak->cos2psi1 = 2*cpsi1*cpsi1 - 1;
    ak->cos2psi2 = 2*cpsi2*cpsi2 - 1;
    ak->cos2psi_ReducedTotalSpin = 2*cpsi_ReducedTotalSpin*cpsi_ReducedTotalSpin - 1;

    ak->Fe = XLALSimInspiralSpinEccEnhancementFlux(etSq);
    ak->LNmag = XLALSimInspiralSpinEccLmagN(omega, etSq, ak);
    ak->SemiMajor = XLALSimInspiralSpinEccSemiMajorAxis(omega, ak);
    ak->SemiMinor = ak->SemiMajor * sqrt(1 - etSq);
    ak->gamma1 = XLALSimInspiralSpinEccGamma1(s1x, s1y,  s1z, s2x,  s2y,  s2z, lnhatx, lnhaty, lnhatz, ak );
    ak->beta43 = XLALSimInspiralSpinEccBeta(4, 3, s1x, s1y, s1z, s2x, s2y, s2z,
                                                lnhatx, lnhaty, lnhatz, ak);

    XLALFree(s1Crosss2);
    XLALFree(lcrosss1);
    XLALFree(lcrosss2);
    XLALFree(s1Crossphat);
    XLALFree(s2Crossphat);
    XLALFree(lcrossphat);

    return XLAL_SUCCESS;
}

/**
 * @brief Compute the cross product of two vectors
 *
 * This function calculates the cross product of two vectors, storing the result in a dynamically allocated array.
 *
 * @param[out] vout   Pointer to a pointer that will hold the allocated array of 3 REAL8 values
 * @param[in]  v1x,v1y,v1z  Components of the first vector
 * @param[in]  v2x,v2y,v2z  Components of the second vector
 *
 * @return XLAL_SUCCESS on successful completion
 *
 * @note Allocated memory should be freed using LALFree() after use of the function
 */
static int XLALSimInspiralVectorCrossProduct(REAL8 **vout,
                        REAL8 v1x, REAL8 v1y, REAL8 v1z,
                        REAL8 v2x, REAL8 v2y, REAL8 v2z)
{
    (*vout) = (double *) LALMalloc(sizeof(double) * 3);
    (*vout)[0] = v1y*v2z-v1z*v2y;
    (*vout)[1] = v1z*v2x-v1x*v2z;
    (*vout)[2] = v1x*v2y-v1y*v2x;
    return XLAL_SUCCESS;
}

/**
 * @brief Compute the 2PN spin-spin \f$\sigma\f$ factor
 *
 * This function evaluates the 2PN spin-spin \f$\sigma\f$ factor that enters the orbital
 * frequency and eccentricity evolution.
 *
 * @param[in] a,b,c Sigma function coefficients
 * @param[in] a1a2,b1b2,c1c2 Sigma function coefficients specialized to BH spin-induced quadrupoles; d, e, f in Phukon et al. are referred as a1a2, b1b2, c1c2 here
 * @param[in] s1x, s1y, s1z Cartesian components of dimensionless spin of body 1
 * @param[in] s2x, s2y, s2z Cartesian components of dimensionless spin of body 2
 * @param[in] lhatx, lhaty, lhatz Unit vector components of orbital angular momentum
 * @param[in] ak Structure containing mass, angular parameters
 *
 * @return The 2PN spin-spin sigma factor
 *
 * @see Equation A4b of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccSigmaK ( REAL8 a, REAL8 b, REAL8 c,
        REAL8 a1a2, REAL8 b1b2, REAL8 c1c2,
        REAL8 s1x, REAL8 s1y, REAL8 s1z,
        REAL8 s2x, REAL8 s2y, REAL8 s2z,
        REAL8 lhatx, REAL8 lhaty, REAL8 lhatz,
        expnCoeffsSpinEcc *ak)

{
   REAL8 mass1 = ak->m1;
   REAL8 mass2 = ak->m2;

   REAL8 lhatCrossReducedS1NormSq=0.;
   REAL8 lhatCrossReducedS2NormSq=0.;
   REAL8 lhatCrossReducedTotalSpinNormSq=0.;

   REAL8 lhatdotReducedS1=0.;
   REAL8 lhatdotReducedS2=0.;
   REAL8 lhatdotReducedTotalSpin=0.;

   REAL8 s1x_reduced = mass1*s1x;
   REAL8 s1y_reduced = mass1*s1y;
   REAL8 s1z_reduced = mass1*s1z;

   REAL8 s2x_reduced = mass2*s2x;
   REAL8 s2y_reduced = mass2*s2y;
   REAL8 s2z_reduced = mass2*s2z;

   REAL8 sx_reduced = s1x_reduced + s2x_reduced;
   REAL8 sy_reduced = s1y_reduced + s2y_reduced;
   REAL8 sz_reduced = s1z_reduced + s2z_reduced;

   REAL8 ReducedTotalSpinNormSq = sx_reduced*sx_reduced + sy_reduced*sy_reduced + sz_reduced*sz_reduced;
   REAL8 ReducedS1magSq = s1x_reduced*s1x_reduced + s1y_reduced*s1y_reduced + s1z_reduced*s1z_reduced;
   REAL8 ReducedS2magSq = s2x_reduced*s2x_reduced + s2y_reduced*s2y_reduced + s2z_reduced*s2z_reduced;

   XLALSimInspiralComputeCrossDotQuants( &lhatCrossReducedS1NormSq, &lhatCrossReducedS2NormSq,
                                        &lhatdotReducedS1, &lhatdotReducedS2,
                                        s1x_reduced, s1y_reduced, s1z_reduced,
                                        s2x_reduced, s2y_reduced, s2z_reduced, lhatx, lhaty, lhatz);

   XLALSimInspiralComputeCrossDotQuants( &lhatCrossReducedTotalSpinNormSq, &(REAL8){1},
                                        &lhatdotReducedTotalSpin, &(REAL8){1},
                                        sx_reduced, sy_reduced, sz_reduced,
                                        0, 0, 0, lhatx, lhaty, lhatz);

   REAL8 sigmaK = a*ReducedTotalSpinNormSq + b*lhatdotReducedTotalSpin*lhatdotReducedTotalSpin +
                  c*lhatCrossReducedTotalSpinNormSq * (*ak).cos2psi_ReducedTotalSpin  +
                  a1a2*(ReducedS1magSq + ReducedS2magSq ) +
                  b1b2*(lhatdotReducedS1*lhatdotReducedS1 + lhatdotReducedS2*lhatdotReducedS2) +
                  c1c2*(lhatCrossReducedS1NormSq*(*ak).cos2psi1 + lhatCrossReducedS2NormSq*(*ak).cos2psi2);
   return sigmaK;
}

/**
 * @brief Compute the 3PN \f$F\f$ flux enhancement function
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 *
 * @return The 3PN \f$F(e)\f$ function
 *
 * @see Equation A5a of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccEnhancementFlux( REAL8 etpow2 )
{

    REAL8 Fe;
    Fe  = 1 + etpow2*(85./6 + etpow2 * (5171./192 + etpow2*(1751./192 + 297./1024*etpow2)));
    Fe = Fe/pow(1-etpow2, 13./2);
    return Fe;
}

/**
 * @brief Compute the 3PN \f$\tilde{F}\f$ enhancement function
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 *
 * @return The 3PN \f$\tilde{F}\f$ flux enhancement function
 *
 * @see Equation A5b of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccTildeEnhancementFlux( REAL8 etpow2 )
{

    REAL8 TildeFe = (1+ etpow2*(229./32 + etpow2*(327./64 + etpow2*69./256) ))/pow(1-etpow2,5);
    return TildeFe;
}

/**
 * @brief Compute the 3PN \f$F_e(e)\f$ enhancement function
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 *
 * @return 3PN \f$F_e(e)\f$ enhancement function
 *
 * @see Equation A5c of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccEnhancementFluxE( REAL8 etpow2 )
{
    REAL8 sqrt_1_etpow2 = sqrt( 1 - etpow2);
    REAL8 f = (1 +  etpow2*( 2782./769 +
            etpow2*( 10721./6152 +
		     1719./24608 * etpow2)))/pow(sqrt_1_etpow2, 11);
    return f;
}

/**
 * @brief Compute the 2.5PN \f$\psi_\omega\f$ enhancement function using the selected enhancement function
 *
 * @param[in] etSq Square of the eccentricity \f$e_t^2\f$
 * @param[in] enhancementFunc Choice of enhancement function:
 *                            0 for O(e^4), 1 for superasymptotic, 2 for hyperasymptotic
 * @param[in] EccOrder Eccentricity order for the hyperasymptotic enhancement function: Available eccentricity orders 2 to 20 or -1 for highest order
 *
 * @return The 2.5PN \f$\psi_\omega\f$ enhancement function
 *
 * @see Equation A6a of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccPsiOmega( REAL8 etSq, EnhancementFunction enhancementFunc, LALSimInspiralSpinEccPolyEccOrder EccOrder )
{
    REAL8 sqrt_one_minus_etSq = sqrt(1 - etSq);
    REAL8 one_minus_etSq_pow_3_by_2 = sqrt_one_minus_etSq*(1 - etSq);

    REAL8 Phi, TildePhi, Psi;

    switch ( enhancementFunc )
    {
        case (EF_Arun_etal):
            Phi = XLALSimInspiralSpinEccPhi(etSq);
            TildePhi = XLALSimInspiralSpinEccTildePhi(etSq);
            Psi = XLALSimInspiralSpinEccPsi(etSq);
            break;
        case (EF_LoutrelYunes_SuperAsym):
            Phi = XLALSimInspiralSpinEccPhiSuperAsyLY(etSq);
            TildePhi = XLALSimInspiralSpinEccTildePhiSuperAsyLY(etSq);
            Psi = XLALSimInspiralSpinEccPsiSuperAsyLY(etSq);
            break;
        case (EF_LoutrelYunes_HyperAsym):
            Phi = XLALSimInspiralSpinEccPhiHyperAsyLY(etSq, EccOrder);
            TildePhi = XLALSimInspiralSpinEccTildePhiHyperAsyLY(etSq, EccOrder);
            Psi = XLALSimInspiralSpinEccPsiHyperAsyLY(etSq, EccOrder);
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid EnhancementFunction\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }


    REAL8 PsiOmega = 0;
    PsiOmega = (1344.*( sqrt_one_minus_etSq *( 1 - 5 * etSq ) * Phi - 4 * TildePhi )/one_minus_etSq_pow_3_by_2 +
               8191. * Psi) /4159.;
    return PsiOmega;
}

/**
 * @brief Compute the 2.5PN \f$\zeta_\omega\f$ enhancement function using the selected enhancement function
 *
 * @param[in] etSq  Square of the eccentricity \f$e_t^2\f$
 * @param[in] enhancementFunc  Choice of enhancement function:
 *                             0 for O(e^4), 1 for superasymptotic, 2 for hyperasymptotic
 * @param[in] EccOrder  Eccentricity order for the hyperasymptotic enhancement function: Available eccentricity orders 2 to 20 or -1 for highest order
 *
 * @return The 2.5PN \f$\zeta_\omega\f$ enhancement function
 *
 * @see Equation A6b of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccZetaOmega( REAL8 etSq, EnhancementFunction enhancementFunc, LALSimInspiralSpinEccPolyEccOrder EccOrder )
{

    REAL8 Phi, Zeta;
    switch ( enhancementFunc )
    {
        case (EF_Arun_etal):
            Phi = XLALSimInspiralSpinEccPhi(etSq);
            Zeta = XLALSimInspiralSpinEccZeta(etSq);
            break;
        case (EF_LoutrelYunes_SuperAsym):
            Phi = XLALSimInspiralSpinEccPhiSuperAsyLY(etSq);
            Zeta = XLALSimInspiralSpinEccZetaSuperAsyLY(etSq);
            break;
        case (EF_LoutrelYunes_HyperAsym):
            Phi = XLALSimInspiralSpinEccPhiHyperAsyLY(etSq, EccOrder);
            Zeta = XLALSimInspiralSpinEccZetaHyperAsyLY(etSq, EccOrder);
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid EnhancementFunction\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }
    REAL8 ZetaOmega = 0;
    ZetaOmega = (583.*Zeta - 16.*Phi)/567;
    return ZetaOmega;
}

/**
 * @brief Compute the 1.5PN \f$\phi_e(e)\f$ enhancement function multiplied by \f$e_t^2\f$ using the selected enhancement function
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 * @param[in] enhancementFunc Choice of enhancement function:
 *                            0 for O(e^4), 1 for superasymptotic, 2 for hyperasymptotic
 * @param[in] EccOrder Eccentricity order for the hyperasymptotic enhancement function: Available eccentricity orders 2 to 20 or -1 for highest order
 *
 * @return The 1.5PN \f$e_t^2 \times \phi_e\f$ enhancement function
 *
 * @see Equation A6c of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccPhiE_mult_etpow2( REAL8 etpow2, EnhancementFunction enhancementFunc,
                                                LALSimInspiralSpinEccPolyEccOrder EccOrder )
{
    REAL8 sqrt_one_minus_etSq = sqrt(1-etpow2);
    REAL8 Phi, TildePhi;
    switch(enhancementFunc)
    {
        case(EF_Arun_etal):
            Phi = XLALSimInspiralSpinEccPhi(etpow2);
            TildePhi = XLALSimInspiralSpinEccTildePhi(etpow2);
            break;
        case(EF_LoutrelYunes_SuperAsym):
            Phi = XLALSimInspiralSpinEccPhiSuperAsyLY(etpow2);
            TildePhi = XLALSimInspiralSpinEccTildePhiSuperAsyLY(etpow2);
            break;
        case(EF_LoutrelYunes_HyperAsym):
            Phi = XLALSimInspiralSpinEccPhiHyperAsyLY(etpow2, EccOrder);
            TildePhi = XLALSimInspiralSpinEccTildePhiHyperAsyLY (etpow2, EccOrder);
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid Enhancement Function\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    REAL8 phie_mult_etpow2;
    phie_mult_etpow2 = 192./985*sqrt_one_minus_etSq*( sqrt_one_minus_etSq*Phi
                        - TildePhi);
    return phie_mult_etpow2;
}

/**
 * @brief Compute the \f$\psi_e\f$ enhancement function multiplied by \f$e_t^2\f$ with the selected enhancement function
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 * @param[in] enhancementFunc Choice of enhancement function:
 *                            0 for O(e^4), 1 for superasymptotic, 2 for hyperasymptotic
 * @param[in] EccOrder Eccentricity order for the hyperasymptotic enhancement function: Available eccentricity orders 2 to 20 or -1 for highest order
 *
 * @return The 2.5PN \f$e_t^2 \times \psi_e\f$ enhancement function
 *
 * @see Equation A6d of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccPsiE_mult_etpow2( REAL8 etpow2, EnhancementFunction enhancementFunc,
                                                LALSimInspiralSpinEccPolyEccOrder EccOrder )
{
    REAL8 sqrt_1_etpow2 = sqrt( 1 - etpow2);

    REAL8 Phi, TildePhi, Psi, TildePsi;

    switch(enhancementFunc)
    {
        case(EF_Arun_etal):
            Phi = XLALSimInspiralSpinEccPhi(etpow2);
            TildePhi = XLALSimInspiralSpinEccTildePhi(etpow2);
            Psi = XLALSimInspiralSpinEccPsi(etpow2);
            TildePsi = XLALSimInspiralSpinEccTildePsi(etpow2);
            break;
        case(EF_LoutrelYunes_SuperAsym):
            Phi = XLALSimInspiralSpinEccPhiSuperAsyLY(etpow2);
            TildePhi = XLALSimInspiralSpinEccTildePhiSuperAsyLY(etpow2);
            Psi = XLALSimInspiralSpinEccPsiSuperAsyLY(etpow2);
            TildePsi = XLALSimInspiralSpinEccTildePsiSuperAsyLY(etpow2);
            break;
       case(EF_LoutrelYunes_HyperAsym):
            Phi = XLALSimInspiralSpinEccPhiHyperAsyLY(etpow2, EccOrder);
            TildePhi = XLALSimInspiralSpinEccTildePhiHyperAsyLY(etpow2, EccOrder);
            Psi = XLALSimInspiralSpinEccPsiHyperAsyLY(etpow2, EccOrder);
            TildePsi = XLALSimInspiralSpinEccTildePsiHyperAsyLY(etpow2, EccOrder);
            break;
       default:
            XLALPrintError("XLAL Error - %s: Invalid EnhancementFunction\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }
    REAL8 PsiE_mult_etpow2;
    PsiE_mult_etpow2 = 18816./55691.*((1-
                11./7*etpow2)*Phi - (1-3./7*etpow2)/sqrt_1_etpow2*TildePhi)+
                16382./55691*((1-etpow2)*Psi - sqrt_1_etpow2*TildePsi);
    return PsiE_mult_etpow2;
}

/**
 * @brief Compute the 2.5PN \f$\zeta_e(e)\f$ enhancement function multiplied by \f$e_t^2\f$ using the selected enhancement function
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 * @param[in] enhancementFunc Choice of enhancement function:
 *                            0 for O(e^4), 1 for superasymptotic, 2 for hyperasymptotic
 * @param[in] EccOrder Eccentricity order for the hyperasymptotic enhancement function: Available eccentricity orders 2 to 20 or -1 for highest order
 *
 * @return The 2.5PN \f$e_t^2 \times \zeta_e \f$ enhancement function
 *
 * @see Equation A6e of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccZetaE_mult_etpow2( REAL8 etpow2, EnhancementFunction enhancementFunc,
                                                LALSimInspiralSpinEccPolyEccOrder EccOrder )
{
    REAL8 sqrt_1_etpow2 = sqrt( 1 - etpow2);

    REAL8 Phi, TildePhi, Zeta, TildeZeta;
    switch ( enhancementFunc )
    {
        case EF_Arun_etal:
            Phi = XLALSimInspiralSpinEccPhi(etpow2);
            TildePhi = XLALSimInspiralSpinEccTildePhi(etpow2);
            Zeta = XLALSimInspiralSpinEccZeta(etpow2);
            TildeZeta = XLALSimInspiralSpinEccTildeZeta(etpow2);
            break;
        case EF_LoutrelYunes_SuperAsym:
            Phi = XLALSimInspiralSpinEccPhiSuperAsyLY(etpow2);
            TildePhi = XLALSimInspiralSpinEccTildePhiSuperAsyLY(etpow2);
            Zeta = XLALSimInspiralSpinEccZetaSuperAsyLY(etpow2);
            TildeZeta = XLALSimInspiralSpinEccTildeZetaSuperAsyLY(etpow2);
            break;
        case EF_LoutrelYunes_HyperAsym:
            Phi = XLALSimInspiralSpinEccPhiHyperAsyLY(etpow2, EccOrder);
            TildePhi =  XLALSimInspiralSpinEccTildePhiHyperAsyLY (etpow2, EccOrder);
            Zeta = XLALSimInspiralSpinEccZetaHyperAsyLY(etpow2, EccOrder);
            TildeZeta = XLALSimInspiralSpinEccTildeZetaHyperAsyLY(etpow2, EccOrder);
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid EnhancementFunction\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    REAL8 ZetaE_mult_etpow2;
    ZetaE_mult_etpow2  = 924./19067*( - (1 - etpow2)*Phi +
                (1 - 5./11*etpow2)*TildePhi/sqrt_1_etpow2)+
                12243./76268*sqrt_1_etpow2*(sqrt_1_etpow2*Zeta - TildeZeta);

    return ZetaE_mult_etpow2;
}


/**
 * @brief Compute the 3PN \f$\kappa_e\f$ enhancement function multiplied by \f$e_t^2\f$ using the selected enhancement function
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 * @param[in] enhancementFunc Choice of enhancement function:
 *                            0 for O(e^4), 1 for superasymptotic, 2 for hyperasymptotic
 * @param[in] EccOrder Eccentricity order for the hyperasymptotic enhancement function: Available eccentricity orders 2 to 20 or -1 for highest order
 *
 * @return The 3PN \f$e_t^2 \times \kappa_e\f$ enhancement term
 *
 * @see Equation A6f of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccEnhancementKappaE_mult_etpow2( REAL8 etpow2, EnhancementFunction enhancementFunc, LALSimInspiralSpinEccPolyEccOrder EccOrder )
{
    REAL8 sqrt_1_etpow2 = sqrt( 1 - etpow2);

    REAL8 KappaE_mult_etpow2;
    REAL8 Kappa, TildeKappa;

    switch ( enhancementFunc )
    {
        case EF_Arun_etal:
            Kappa = XLALSimInspiralSpinEccEnhancementKappa(etpow2);
            TildeKappa = XLALSimInspiralSpinEccEnhancementTildeKappa(etpow2);
            break;
        case EF_LoutrelYunes_SuperAsym:
            Kappa = XLALSimInspiralSpinEccEnhancementKappaSuperAsyLY(etpow2);
            TildeKappa = XLALSimInspiralSpinEccEnhancementTildeKappaSuperAsyLY(etpow2);
            break;
        case EF_LoutrelYunes_HyperAsym:
            Kappa = XLALSimInspiralSpinEccEnhancementKappaHyperAsyLY(etpow2, EccOrder);
            TildeKappa = XLALSimInspiralSpinEccEnhancementTildeKappaHyperAsyLY(etpow2, EccOrder);
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid EnhancementFunction\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    KappaE_mult_etpow2 = sqrt_1_etpow2 * ( sqrt_1_etpow2*Kappa - TildeKappa )/
                        (769./96 - 3059665./700566*LAL_LN2 + 8190315./1868176*LAL_LN3);
    return KappaE_mult_etpow2;
}

/**
 * @brief Compute the 1.5 PN \f$\phi_k\f$ enhancement function for specified enhancement fuction
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 * @param[in] enhancementFunc Choice of enhancement function:
 *                            0 for O(e^4), 1 for superasymptotic, 2 for hyperasymptotic
 * @param[in] EccOrder Eccentricity order for the hyperasymptotic enhancement function: Available eccentricity orders 2 to 20 or -1 for highest order
 *
 * @return The 1.5PN \f$\phi_k(e)\f$ enhancement function
 *
 * @see Equation A6g of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccPhik ( REAL8 etpow2, EnhancementFunction enhancementFunc,
                                    LALSimInspiralSpinEccPolyEccOrder EccOrder )
{
    REAL8 TildePhi = 0;

    switch(enhancementFunc)
    {
        case(EF_Arun_etal):
            TildePhi = XLALSimInspiralSpinEccTildePhi(etpow2);
            break;
        case(EF_LoutrelYunes_SuperAsym):
            TildePhi = XLALSimInspiralSpinEccTildePhiSuperAsyLY(etpow2);
            break;
        case(EF_LoutrelYunes_HyperAsym):
            TildePhi = XLALSimInspiralSpinEccTildePhiHyperAsyLY (etpow2, EccOrder);
            break;
        default:
            XLALPrintError("XLAL Error in %s: %d is an invalid Enhancement Function value\n",__func__,  enhancementFunc);
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    REAL8 Phik =  TildePhi/pow(1 - etpow2, 3./2);

    return Phik;
}

/**
 * @brief Compute the O(e^4) 1.5PN \f$\phi\f$ enhancement function
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 *
 * @return The O(e^4) 1.5PN \f$\phi\f$ enhancement function
 *
 * @see Equation A7a of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccPhi( REAL8 etpow2 )
{
    REAL8 phi = 1 + etpow2 * (2335./192. + 42955./768.* etpow2);
    return phi;
}

/**
 * @brief Compute the O(e^4) 2.5PN \f$\psi\f$ enhancement function
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 *
 * @return The O(e^4) 2.5PN \f$\psi\f$ enhancement function
 *
 * @see Equation A7b of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccPsi( REAL8 etpow2 )
{

    REAL8 psi = 1 - etpow2*( 22988./8191 +  36508643./524224*etpow2);
    return psi;
}

/**
 * @brief Compute the O(e^4) 2.5PN \f$\zeta\f$ enhancement function
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 *
 * @return The O(e^4) 2.5PN \f$\zeta\f$ enhancement function
 *
 * @see Equation A7c of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccZeta( REAL8 etpow2 )
{

    REAL8 zeta = 1 + etpow2 * ( 1011565./48972 + 106573021./783552*etpow2);
    return zeta;
}

/**
 * @brief Compute the O(e^4) 3PN \f$\tilde{\kappa}\f$ enhancement function
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 *
 * @return The O(e^4) 3PN \f$\tilde{\kappa}\f$ enhancement function
 *
 * @see Equation A7d of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccEnhancementKappa( REAL8 etpow2 )
{
    REAL8 kappa = 1 + etpow2*(( 62./3 - 4613840./350283 * LAL_LN2 + 24570945./1868176* LAL_LN3 ) +
                      etpow2*( 9177./64 + 271636085./1401132 * LAL_LN2 - 466847955./7472704 * LAL_LN3 ));
    return kappa;
}

/**
 * @brief Compute the O(e^4) 1.5PN \f$\tilde{\phi}\f$ enhancement function
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 *
 * @return The O(e^4) 1.5PN \f$\tilde{\phi}\f$ enhancement function
 *
 * @see Equation A8a of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccTildePhi( REAL8 etpow2 )
{
    return 1 + etpow2 * (209./32 + 2415./128*etpow2);
}

/**
 * @brief Compute the O(e^4) 2.5PN \f$\tilde{\psi}\f$ enhancement function
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 *
 * @return The O(e^4) 2.5PN \f$\tilde{\psi}\f$ enhancement function
 *
 * @see Equation A8b of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccTildePsi( REAL8 etpow2 )
{
    REAL8 TildePsi = 1 - etpow2*( 17416./8191 +  14199197./524224*etpow2);
    return TildePsi;
}

/**
 * @brief Compute the O(e^4) 2.5PN \f$\tilde{\zeta}\f$ enhancement function
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 *
 * @return The O(e^4) 2.5PN \f$\tilde{\zeta}\f$ enhancement function
 *
 * @see Equation A8c of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccTildeZeta( REAL8 etpow2 )
{
    REAL8 TildeZeta = 1 + etpow2*(102371./8162 + 14250725./261184*etpow2);
    return TildeZeta;
}

/**
 * @brief Compute the O(e^4) 3PN \f$\tilde{\kappa}\f$ enhancement function
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 *
 * @return The O(e^4) 3PN \f$\tilde{\kappa}\f$ enhancement function
 *
 * @see Equation A8d of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccEnhancementTildeKappa( REAL8 etpow2 )
{
    REAL8 TildeKapp = 1 + etpow2*((389./32 - 2056005./233522*LAL_LN2 + 8190315./934088*LAL_LN3) +
                          etpow2* (3577./64 + 50149295./467044*LAL_LN2 -
                         155615985./3736352*LAL_LN3));
    return TildeKapp;
}

/**
 * @brief Compute the superasymptotic 1.5PN \f$\phi\f$ enhancement function
 *
 * @param[in] etSq Square of the eccentricity \f$e_t^2\f$
 *
 * @return The superasymptotic 1.5PN \f$\phi\f$ enhancement function
 *
 * @see Equation 140 of Loutrel and Yunes, arXiv:1607.05409
 */
REAL8 XLALSimInspiralSpinEccPhiSuperAsyLY( REAL8 etSq)
{
    REAL8 one_minus_etSq = 1 - etSq;

    REAL8 phi = (1./( LAL_GAMMA_1_3_MULT_GAMMA_2_3 * pow(one_minus_etSq, 5) )) * ( 1328./27 +
                                   one_minus_etSq * ( - 992./15 +
                                   one_minus_etSq * ( 33982./1575 - 1577./1575 * one_minus_etSq )));
    return phi;
}

/**
 * @brief Compute the hyperasymptotic 1.5PN \f$\phi\f$ enhancement function
 *
 * This implementation adds the superasymptotic value `refXLALSimInspiralSpinEccPhiSuperAsyLY`
 * with a hyperasymptotic correction in form of a polynomial in eccentricity.
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 * @param[in] EccOrder Order of eccentricity for the hyperasymptotic correction;
 *                     even number between 2 and 20, or -1 for the maximum available order
 *
 * @return The hyperasymptotic 1.5PN \f$\phi\f$ enhancement function
 *
 * @see Equation 176 of Loutrel and Yunes, arXiv:1607.05409
 */
REAL8 XLALSimInspiralSpinEccPhiHyperAsyLY( REAL8 etpow2, LALSimInspiralSpinEccPolyEccOrder EccOrder )
{
    REAL8 inverse_gamma_mult = 1./LAL_GAMMA_1_3_MULT_GAMMA_2_3;
    REAL8 coeff0 = 1 - XLALSimInspiralSpinEccPhiSuperAsyLY(0.);
    REAL8 coeff2 = 2335./192 - 208456./4725*inverse_gamma_mult;
    REAL8 coeff4 = 42955./768 - 319561./1575*inverse_gamma_mult;
    REAL8 coeff6 = 6204647./36864 - 2884936./4725*inverse_gamma_mult;
    REAL8 coeff8 = 352891481./884736 - 1367347./945*inverse_gamma_mult;
    REAL8 coeff10 = 286907786543./353894400 - 61760./21*inverse_gamma_mult;
    REAL8 coeff12 = 6287456255443./4246732800 - 1208431./225*inverse_gamma_mult;
    REAL8 coeff14 = 5545903772613817./2219625676800 - 4758512./525*inverse_gamma_mult;
    REAL8 coeff16 = 422825073954708079./106542032486400 - 7558199./525*inverse_gamma_mult;
    REAL8 coeff18 = 1659160118498286776339./276156948204748800 - 20596024./945*inverse_gamma_mult;
    REAL8 coeff20 = 724723372042305.454448081/82847084461.424640000 - 29987903./945*inverse_gamma_mult;

    REAL8 deltaPhi = 0;

    switch( EccOrder )
    {
        case(-1):
        case(20):
            deltaPhi = coeff20;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(18):
	        deltaPhi *= etpow2;
            deltaPhi += coeff18;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(16):
	        deltaPhi *=	etpow2;
            deltaPhi += coeff16;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(14):
	        deltaPhi *=	etpow2;
            deltaPhi += coeff14;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(12):
	        deltaPhi *=	etpow2;
	        deltaPhi += coeff12;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(10):
	        deltaPhi *= etpow2;
            deltaPhi += coeff10;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(8):
	        deltaPhi *= etpow2;
            deltaPhi += coeff8;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(6):
	        deltaPhi *= etpow2;
            deltaPhi += coeff6;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(4):
	        deltaPhi *= etpow2;
            deltaPhi += coeff4;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(2):
	        deltaPhi *= etpow2;
	        deltaPhi += coeff2;
	        deltaPhi *= etpow2;
            deltaPhi += coeff0;
            break;
        default:
            XLALPrintError("XLAL Error - %d is an invalid eccentricity order in the enhancement function, must be even and between 2 and 20 or -1\n", EccOrder );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }
    REAL8 phi = XLALSimInspiralSpinEccPhiSuperAsyLY( etpow2 ) + deltaPhi;
    return phi;
}

/**
 * @brief Compute the superasymptotic 2.5PN \f$ \alpha \f$ enhancement function
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 *
 * @return The superasymptotic 2.5PN \f$ \alpha \f$ enhancement function
 *
 * @see Equation 169 of Loutrel and Yunes, arXiv:1607.05409
 */
REAL8 XLALSimInspiralSpinEccAlphaSuperAsyLY( REAL8 etpow2 )
{
    REAL8 one_minus_etSq = 1 - etpow2;
    REAL8 sqrt_one_minus_etSq = sqrt(one_minus_etSq);

    REAL8 alpha =  (( 77776./321 + sqrt_one_minus_etSq * ( - 15904./107 +
                                   sqrt_one_minus_etSq * ( - 300512./963  +
                                   sqrt_one_minus_etSq * (  19530./107  +
                                   sqrt_one_minus_etSq * (  4871974./56175 +
                                   sqrt_one_minus_etSq * ( - 26952./535  +
                                   sqrt_one_minus_etSq * ( 111533./56175  +
                                   sqrt_one_minus_etSq * ( 8313./5350 +
                                   sqrt_one_minus_etSq * ( - 2280749./6179250 -
                                   sqrt_one_minus_etSq * 4293./294250  )))))))))/LAL_GAMMA_1_3_MULT_GAMMA_2_3 -
                     21./3424 * XLALSimInspiralSpinEccPadePolyA(etpow2))/pow(one_minus_etSq,6);
   return alpha;
}

/**
 * @brief Compute the superasymptotic 2.5PN \f$ \beta \f$ enhancement function
 *
 * @param[in] etSq Square of the eccentricity \f$e_t^2\f$
 *
 * @return The superasymptotic 2.5PN \f$ \beta \f$ enhancement function
 *
 * @see Equation 142 of Loutrel and Yunes, arXiv:1607.05409
 */
REAL8 XLALSimInspiralSpinEccBetaSuperAsyLY( REAL8 etSq )
{
     REAL8 one_minus_etSq = 1 - etSq;

     REAL8 beta = (1./( LAL_GAMMA_1_3_MULT_GAMMA_2_3 * pow(one_minus_etSq, 6) )) * ( 7244800./49209 +
                        one_minus_etSq * (- 39162880./147627 +
                        one_minus_etSq * ( 16731136./114821 +
                        one_minus_etSq * ( - 14146304./574105 +
                        one_minus_etSq * ( 1052528./1722315  - 260144./2052425375 * one_minus_etSq )))));
     return beta;
}

/**
 * @brief Compute the superasymptotic 2.5PN \f$ \gamma \f$ enhancement function
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 *
 * @return The superasymptotic 2.5PN \f$ \gamma \f$ enhancement function
 *
 * @see Equation 143 of Loutrel and Yunes, arXiv:1607.05409
 */
REAL8 XLALSimInspiralSpinEccGammaSuperAsyLY( REAL8 etpow2 )
{
    REAL8 one_minus_etSq = 1 - etpow2;
    REAL8 one_minus_etSq_pow6 = pow(one_minus_etSq, 6);

    REAL8 gamma = (1./(LAL_GAMMA_1_3_MULT_GAMMA_2_3 * one_minus_etSq_pow6))*( 1280./3 +
                        one_minus_etSq * (- 7808./9 +
                        one_minus_etSq * (97472./175 +
                        one_minus_etSq * ( - 20368./175 +
                        one_minus_etSq * ( 113228./28875 - 122./625625*one_minus_etSq)))));
    return gamma;
}

/**
 * @brief Compute the superasymptotic 2.5PN \f$ \theta \f$ enhancement function
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 *
 * @return The superasymptotic 2.5PN \f$ \theta \f$ enhancement function
 *
 * @see Equation 152 of Loutrel and Yunes, arXiv:1607.05409
 */
REAL8 XLALSimInspiralSpinEccThetaSuperAsyLY( REAL8 etpow2 )
{
    REAL8 one_minus_etSq = 1 - etpow2;
    REAL8 one_minus_etSq_pow6 = pow(one_minus_etSq, 6);

    REAL8 theta = (1./(LAL_GAMMA_1_3_MULT_GAMMA_2_3*one_minus_etSq_pow6))*( 34240./801 +
                                   one_minus_etSq*( - 132560./2403 +
                                   one_minus_etSq*( 794344./46725 +
                                   one_minus_etSq*( - 141994./140175 +
                                   one_minus_etSq*(465188./7709625  - 500627./334083750 * one_minus_etSq )))));
    return theta;
}

/**
 * @brief Compute the superasymptotic 3PN \f$ \chi \f$ enhancement function
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 *
 * @return The superasymptotic 3PN \f$ \chi \f$ enhancement function
 *
 * @see Equation 171 of Loutrel and Yunes, arXiv:1607.05409
 */
REAL8 XLALSimInspiralSpinEccEnhancementChiSuperAsyLY( REAL8 etpow2 )
{
    REAL8 one_minus_etSq = 1 - etpow2;
    REAL8 sqrt_one_minus_etSq = sqrt( one_minus_etSq );
    REAL8 one_minus_etpow_13_2 = pow(sqrt_one_minus_etSq, 13);
    REAL8 log_one_minus_etSq = log(one_minus_etSq);

    REAL8 chi = (1./one_minus_etpow_13_2) * (( 421543./1536 - 52745*LAL_GAMMA/1024 - 52745*LAL_LN2/256
                            - 52745.*LAL_LN3/2048 - 158235.*log_one_minus_etSq/2048) +
                            one_minus_etSq*(( -2777339./5120 + 24717*LAL_GAMMA/256 + 24717*LAL_LN2/64 +
                                24717.* LAL_LN3/512 + 74151.* log_one_minus_etSq/512) +
                             one_minus_etSq*(( 10449133./30720 - 86065.*LAL_GAMMA/1536 - 86065.*LAL_LN2/384 -
                                 86065.*LAL_LN3/3072 - 86065.* log_one_minus_etSq/1024) +
                             one_minus_etSq * (( -1090519./15360 + 7895.*LAL_GAMMA / 768 +  7895.*LAL_LN2/192 +
                                 355271*LAL_LN3/69120 + 5*LAL_LN3/86400 +  7895.*log_one_minus_etSq/512) +
                             one_minus_etSq* (( 760247221./275968000 - 297.*LAL_GAMMA/1024 - 297.*LAL_LN2/256 -
                                 1024063.*LAL_LN3/7096320 - 2521.*5*LAL_LN3/17740800 - 891.*log_one_minus_etSq/2048) +
                             one_minus_etSq * (( 568287127./67267200000 + 71.*LAL_LN2/682500 + 1327283.*LAL_LN3/2882880000 -
                                 19843.*5*LAL_LN3/221760000 - 71.*(8*LAL_LN2 + LAL_LN3)/5460000) +
                             one_minus_etSq * ( -4896210901./4708704000000 + 12270499.*LAL_LN2/24324300000 +
                                 423525727.*LAL_LN3/5448643200000 - 410009.*5*LAL_LN3/139708800000 - 12270499.*( 8*LAL_LN2 + LAL_LN3) /194594400000)))))));
    return chi;
}

/**
 * @brief Compute the hyperasymptotic 3PN \f$ \chi \f$ enhancement function
 *
 * This function evaluates the hyperasymptotic 3PN \f$ \chi \f$ enhancement function,
 * combining the superasymptotic value `ref XLALSimInspiralSpinEccEnhancementChiSuperAsyLY`
 * with a polynomial in eccentricity with specified order.
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 * @param[in] EccOrder Order of eccentricity for the hyperasymptotic correction;
 *                     even number between 2 and 10, or -1 for the maximum available order.
 *                     Values up to 20 are accepted but do not change the result.
 *
 * @return The hyperasymptotic 3PN \f$ \chi \f$ enhancement function
 *
 * @see Equation 186 of Loutrel and Yunes, arXiv:1607.05409
 */
REAL8 XLALSimInspiralSpinEccEnhancementChiHyperAsyLY( REAL8 etpow2, LALSimInspiralSpinEccPolyEccOrder EccOrder )
{
    REAL8 coeff0 = - XLALSimInspiralSpinEccEnhancementChiSuperAsyLY(0.);
    REAL8 coeff2 = - 131766997689301./1448832000000 + 62*LAL_GAMMA/3 + 57*LAL_LN2 +
                    27619*LAL_LN3/768;
    REAL8 coeff4 = -17257920310633973./25113088000000 + 9177* LAL_GAMMA/64 +
                    182657*LAL_LN2/192 - 51243*LAL_LN3/1024;
    REAL8 coeff6 = - 92129506724738033./30135705600000 + 76615* LAL_GAMMA/128 -
                    296449*LAL_LN2/384 +  8680309*LAL_LN3/16384 + 244140625*LAL_LN5/147456;
    REAL8 coeff8 = - 343678592520953093./34440806400000 + 1903055*LAL_GAMMA/1024 +
                    59103559*LAL_LN2/3072 + 1180327577*LAL_LN3/131072 - 10498046875*LAL_LN5/1179648;
    REAL8 coeff10 = -3051437850147557459./114802688000000 + 9732723* LAL_GAMMA/2048 -
                    55480099157*LAL_LN2/1382400 - 4787048773551*LAL_LN3/104857600 +
                    2342041015625*LAL_LN5/113246208 + 33232930569601*LAL_LN7/943718400;

    REAL8 delta_chi = 0;

    switch( EccOrder )
    {
        case(-1):
        case(20):
        case(18):
        case(16):
        case(14):
        case(12):
        case(10):
            delta_chi += coeff10;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(8):
            delta_chi *= etpow2;
            delta_chi += coeff8;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(6):
            delta_chi *= etpow2;
            delta_chi += coeff6;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(4):
            delta_chi *= etpow2;
            delta_chi += coeff4;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(2):
            delta_chi *= etpow2;
	        delta_chi += coeff2;
	        delta_chi *= etpow2;
            delta_chi += coeff0;
            break;
       default:
            XLALPrintError("XLAL Error - %d is an invalid eccentricity order in the enhancement function, must be even and between 2 and 20 or -1\n", EccOrder );
            XLAL_ERROR(XLAL_EINVAL);
            break;

    }
    REAL8 chi = XLALSimInspiralSpinEccEnhancementChiSuperAsyLY(etpow2)+ delta_chi;
    return chi;
}

/**
 * @brief Compute the superasymptotic 2.5PN \f$ \psi \f$ enhancement function
 *
 * The implementation is a weighted sum of `ref XLALSimInspiralSpinEccAlphaSuperAsyLY`,
 * `ref XLALSimInspiralSpinEccBetaSuperAsyLY`, and `ref XLALSimInspiralSpinEccGammaSuperAsyLY`.
 *
 * @param[in] etSq Square of the eccentricity \f$e_t^2\f$
 *
 * @return The superasymptotic 2.5PN \f$ \psi \f$ enhancement function
 *
 * @see Equation A9a of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccPsiSuperAsyLY( REAL8 etSq )
{
    REAL8 Psi = 13696./8191 * XLALSimInspiralSpinEccAlphaSuperAsyLY (etSq) -
                16403./24573 * XLALSimInspiralSpinEccBetaSuperAsyLY( etSq ) -
                112./24573 * XLALSimInspiralSpinEccGammaSuperAsyLY( etSq);
    return Psi;
}

/**
 * @brief Compute the hyperasymptotic 2.5PN \f$ \psi \f$ enhancement function
 *
 * This function evaluates the hyperasymptotic 2.5PN psi enhancement function by adding a
 * polynomial eccentricity correction to the superasymptotic value `ref XLALSimInspiralSpinEccPsiSuperAsyLY`.
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 * @param[in] EccOrder Order of eccentricity for the hyperasymptotic correction;
 *                     even number between 2 and 20, or -1 for the maximum available order
 *
 * @return The hyperasymptotic 2.5PN \f$ \psi \f$ enhancement function
 *
 * @see Equation 182 of Loutrel and Yunes, arXiv:1607.05409
 */
REAL8 XLALSimInspiralSpinEccPsiHyperAsyLY( REAL8 etpow2, LALSimInspiralSpinEccPolyEccOrder EccOrder )
{
    REAL8 inverse_gamma_mult = 1./LAL_GAMMA_1_3_MULT_GAMMA_2_3;

    REAL8 coeff0 = 1 - XLALSimInspiralSpinEccPsiSuperAsyLY(0.);
    REAL8 coeff2 = 188440./8191 -  1283135619824./15373483125*inverse_gamma_mult;
    REAL8 coeff4 = 78746077./524224 - 8377507600624./15373483125*inverse_gamma_mult;
    REAL8 coeff6 = 2769593143./4718016 - 10912663062368./5124494375*inverse_gamma_mult;
    REAL8 coeff8 = 1038414910159./603906048 - 8717702789819./1397589375*inverse_gamma_mult;
    REAL8 coeff10 = 10517947248419./2516275200 - 233112310241024./15373483125*inverse_gamma_mult;
    REAL8 coeff12 = 51677468559131363./5797498060800 - 1325626967291149./40995955000*inverse_gamma_mult;
    REAL8 coeff14 = 14698793962256164697./852232214937600 - 366425963194427./5856565000*inverse_gamma_mult;
    REAL8 coeff16 = 13508357018274827128409./436342894048051200 - 147327907838689583./1311870560000*inverse_gamma_mult;
    REAL8 coeff18 = 462509893308626646120797./8835943604473036800 - 124551349162839959./655935280000*inverse_gamma_mult;
    REAL8 coeff20 = 4766936073001835060.207935793/56550039068627.435520000 -
                            9627820057.257367657/31484.893440000*inverse_gamma_mult;


    REAL8 deltaPsi = 0;

    switch( EccOrder )
    {
        case (-1):
        case(20):
            deltaPsi = coeff20;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(18):
            deltaPsi *= etpow2;
            deltaPsi += coeff18;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(16):
            deltaPsi *= etpow2;
            deltaPsi += coeff16;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(14):
            deltaPsi *= etpow2;
            deltaPsi += coeff14;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(12):
            deltaPsi *= etpow2;
            deltaPsi += coeff12;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(10):
            deltaPsi *= etpow2;
            deltaPsi += coeff10;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(8):
            deltaPsi *= etpow2;
            deltaPsi += coeff8;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(6):
            deltaPsi *= etpow2;
            deltaPsi += coeff6;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(4):
            deltaPsi *= etpow2;
            deltaPsi += coeff4;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(2):
            deltaPsi *= etpow2;
	        deltaPsi += coeff2;
	        deltaPsi *= etpow2;
            deltaPsi += coeff0;
            break;
        default:
            XLALPrintError("XLAL Error - %d is an invalid eccentricity order in the enhancement function, must be even and between 2 and 20 or -1\n", EccOrder );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }
    REAL8 Psi = XLALSimInspiralSpinEccPsiSuperAsyLY (etpow2) + deltaPsi;
    return Psi;
}

/**
 * @brief Superasymptotic expression for the 2.5PN \f$\zeta\f$ enhancement function.
 *
 * The implementation combines the superasymptotic contributions from the superasymptotic functions
 * `\ref XLALSimInspiralSpinEccThetaSuperAsyLY`, `\ref XLALSimInspiralSpinEccBetaSuperAsyLY`,
 * and `\ref XLALSimInspiralSpinEccGammaSuperAsyLY`.
 *
 * @param[in] etpow2 Square of the orbital eccentricity (\f$e_t^2\f$).
 *
 * @return Superasymptotic value of the 2.5PN \f$\zeta\f$ enhancement function.
 *
 * @see Equation A9b of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccZetaSuperAsyLY( REAL8 etpow2 )
{
    REAL8 Zeta = - (1424./4081) * XLALSimInspiralSpinEccThetaSuperAsyLY( etpow2 )  +
                    (16403./12243) * XLALSimInspiralSpinEccBetaSuperAsyLY( etpow2 ) +
                    (16./1749) * XLALSimInspiralSpinEccGammaSuperAsyLY( etpow2);
    return Zeta;
}

/**
 * @brief Hyperasymptotic expression for the 2.5PN \f$\zeta\f$ enhancement function.
 *
 * The hyperasymptoticexpression is defined as
 * \f[\zeta_{hyper}(e) = \zeta_{super}(e) + \Delta\zeta(e)\f]
 * where \f$\zeta_{super}\f$ is the superasymptotic approximation evaluated by `ref XLALSimInspiralSpinEccZetaSuperAsyLY`.
 * \f$\Delta\zeta\f$is a polynomial correction in powers of \f$e^2\f$.
 *
 * @param[in] etpow2 Square of the orbital eccentricity (\f$e_t^2\f$).
 * @param[in] EccOrder Eccentricity correction order: even number between 2 and 20,
 *                     or -1 to use the highest available order.
 *
 * @return Hyperasymptotic value of the 2.5PN \f$\zeta\f$ enhancement function.
 *
 * @see Eq. 184 of Loutrel & Yunes, arXiv:1607.05409
 */
REAL8 XLALSimInspiralSpinEccZetaHyperAsyLY( REAL8 etpow2, LALSimInspiralSpinEccPolyEccOrder EccOrder )
{
    REAL8 coeff0 = 1 - XLALSimInspiralSpinEccZetaSuperAsyLY(0.);
    REAL8 coeff2 =  1011565./48972 - 5165477150408./(68935741875*LAL_GAMMA_1_3_MULT_GAMMA_2_3);
    REAL8 coeff4 = 106573021./783552 - 1619660334008./(3282654375*LAL_GAMMA_1_3_MULT_GAMMA_2_3);
    REAL8 coeff6 = 456977827./854784 - 133691089979528./(68935741875*LAL_GAMMA_1_3_MULT_GAMMA_2_3);
    REAL8 coeff8 = 128491074157./82059264 - 55938524367784./(9847963125*LAL_GAMMA_1_3_MULT_GAMMA_2_3);
    REAL8 coeff10 = 342306246988373./90265190400 - 105369692129672./(7659526875*LAL_GAMMA_1_3_MULT_GAMMA_2_3);
    REAL8 coeff12 = 69677817044303231./8665458278400 - 8704718214568./(298423125*LAL_GAMMA_1_3_MULT_GAMMA_2_3);
    REAL8 coeff14 = 2386244038997979551./154402711142400 - 429418866068552./(7659526875*LAL_GAMMA_1_3_MULT_GAMMA_2_3);
    REAL8 coeff16 = 5987988065386963552943./217399017288499200 - 765322594645592./(7659526875*LAL_GAMMA_1_3_MULT_GAMMA_2_3);
    REAL8 coeff18 = 61448938675545297383797./1329005313235353600 - 1651784262696184./(9847963125*LAL_GAMMA_1_3_MULT_GAMMA_2_3);
    REAL8 coeff20 = 6249104916857243979.130332809/84524737921768.488960000 -
                    18488329739373848./(68935741875*LAL_GAMMA_1_3_MULT_GAMMA_2_3);

    REAL8 delta_zeta = 0;

    switch(EccOrder)
    {
        case(-1):
        case(20):
            delta_zeta = coeff20;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(18):
            delta_zeta *= etpow2;
            delta_zeta += coeff18;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(16):
            delta_zeta *= etpow2;
            delta_zeta += coeff16;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(14):
            delta_zeta *= etpow2;
            delta_zeta += coeff14;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(12):
            delta_zeta *= etpow2;
            delta_zeta += coeff12;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(10):
            delta_zeta *= etpow2;
            delta_zeta += coeff10;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(8):
            delta_zeta *= etpow2;
            delta_zeta += coeff8;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(6):
            delta_zeta *= etpow2;
            delta_zeta += coeff6;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(4):
            delta_zeta *= etpow2;
            delta_zeta += coeff4;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(2):
            delta_zeta *= etpow2;
	        delta_zeta += coeff2;
	        delta_zeta *= etpow2;
            delta_zeta += coeff0;
            break;
        default:
            XLALPrintError("XLAL Error - %d is an invalid eccentricity order in the enhancement function, must be even and between 2 and 20 or -1\n", EccOrder );
            XLAL_ERROR(XLAL_EINVAL);
            break;

    }
    REAL8 zeta = XLALSimInspiralSpinEccZetaSuperAsyLY( etpow2 ) + delta_zeta;
    return zeta;
}

/**
 * @brief Compute the superasymptotic 3PN \f$\kappa\f$ enhancement function
 *
 * This function evaluates the 3PN \f$\kappa\f$ enhancement function using `ref XLALSimInspiralSpinEccEnhancementChiSuperAsyLY` and `ref XLALSimInspiralSpinEccEnhancementFlux`.
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 *
 * @return The superasymptotic 3PN \f$\kappa\f$ enhancement function
 *
 * @see Equation A9c of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccEnhancementKappaSuperAsyLY( REAL8 etpow2 )
{
   REAL8 Kappa = 59920./116761*XLALSimInspiralSpinEccEnhancementChiSuperAsyLY(etpow2) +
       XLALSimInspiralSpinEccEnhancementFlux (etpow2);
   return Kappa;
}

/**
 * @brief Compute the hyperasymptotic 3PN \f$\kappa\f$ enhancement function
 *
 * This function evaluates the 3PN \f$\kappa\f$ enhancement function using `ref XLALSimInspiralSpinEccEnhancementChiHyperAsyLY` and `ref XLALSimInspiralSpinEccEnhancementFlux`.
 *
 * @param[in] etSq Square of the eccentricity \f$e_t^2\f$
 * @param[in] EccOrder Eccentricity order for the hyperasymptotic correction.
 *                     Even number between 2 and 20, or -1 for the maximum available order
 *
 * @return The hyperasymptotic 3PN \f$\kappa\f$ enhancement function
 *
 * @see Equation A9c of Phukon et al., arXiv:2504.20543
 */
REAL8 XLALSimInspiralSpinEccEnhancementKappaHyperAsyLY( REAL8 etSq, LALSimInspiralSpinEccPolyEccOrder EccOrder )
{
   REAL8 Kappa = 59920./116761*XLALSimInspiralSpinEccEnhancementChiHyperAsyLY(etSq, EccOrder) +
       XLALSimInspiralSpinEccEnhancementFlux (etSq);
   return Kappa;
}

/**
 * @brief Compute the superasymptotic 1.5PN \f$\tilde{\phi}\f$ enhancement function
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 * @return The superasymptotic 1.5PN \f$\tilde{\phi}\f$ enhancement function
 *
 * @see Equation 141 of Loutrel and Yunes, arXiv:1607.05409
 */
REAL8 XLALSimInspiralSpinEccTildePhiSuperAsyLY( REAL8 etpow2 )
{
    REAL8 one_minus_etSq = 1 - etpow2;

    REAL8 tildePhi = (1./(LAL_GAMMA_1_3_MULT_GAMMA_2_3 * pow(one_minus_etSq, 7./2) )) * ( 16 +
                        one_minus_etSq * ( -206./15 +
                        one_minus_etSq * (47./35 + 1./50 * one_minus_etSq)));
    return tildePhi;
}

/**
 * @brief Compute the hyperasymptotic 1.5PN \f$\tilde{\phi}\f$ enhancement function
 *
 * This function evaluates the hyperasymptotic 1.5PN phi tilde enhancement
 * function, combining the superasymptotic value`ref XLALSimInspiralSpinEccTildePhiSuperAsyLY`
 * with a polynomial correction in \f$e_t^2\f$ to specified eccentricity order.
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 * @param[in] EccOrder Order of eccentricity for the hyperasymptotic correction;
 *                     even number between 2 and 20, or -1 for the maximum available order
 *
 * @return The hyperasymptotic 1.5PN \f$\tilde{\phi}\f$ enhancement function
 *
 * @see Equation 177 of Loutrel and Yunes, arXiv:1607.05409
 */
REAL8 XLALSimInspiralSpinEccTildePhiHyperAsyLY( REAL8 etpow2, LALSimInspiralSpinEccPolyEccOrder EccOrder )
{
    REAL8 inverse_gamma_mult = 1./LAL_GAMMA_1_3_MULT_GAMMA_2_3;

    REAL8 coeff0 = 1 - XLALSimInspiralSpinEccTildePhiSuperAsyLY(0.);
    REAL8 coeff2 = 209./32 - (49751./2100)*inverse_gamma_mult;
    REAL8 coeff4 = 2415./128 - (574913./8400)*inverse_gamma_mult;
    REAL8 coeff6 = 730751./18432 - (23011./160)*inverse_gamma_mult;
    REAL8 coeff8 = 10355719./147456 - (326097./1280)*inverse_gamma_mult;
    REAL8 coeff10 = 6594861233./58982400 - (5191733./12800)*inverse_gamma_mult;
    REAL8 coeff12 = 23422887967./141557760 - (30732361./51200)*inverse_gamma_mult;
    REAL8 coeff14 = 51535146547541./221962567680 - (603727553./716800)*inverse_gamma_mult;
    REAL8 coeff16 = 16666910315347223./53271016243200 - (2603342599./2293760)*inverse_gamma_mult;
    REAL8 coeff18 = 8055842533080274417./19725496300339200 - (20389261321./13762560)*inverse_gamma_mult;
    REAL8 coeff20 = 1024885995293794354963./1972549630033920000 - (74113622297./39321600)*inverse_gamma_mult;

    REAL8 deltaPhi = 0;

    switch(EccOrder)
    {
        case(-1):
        case(20):
            deltaPhi = coeff20;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(18):
            deltaPhi *= etpow2;
            deltaPhi += coeff18;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(16):
            deltaPhi *= etpow2;
            deltaPhi += coeff16;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(14):
            deltaPhi *= etpow2;
            deltaPhi += coeff14;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(12):
            deltaPhi *= etpow2;
            deltaPhi += coeff12;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(10):
            deltaPhi *= etpow2;
            deltaPhi += coeff10;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(8):
            deltaPhi *= etpow2;
            deltaPhi += coeff8;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(6):
            deltaPhi *= etpow2;
            deltaPhi += coeff6;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(4):
            deltaPhi *= etpow2;
            deltaPhi += coeff4;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(2):
            deltaPhi *= etpow2;
	        deltaPhi += coeff2;
	        deltaPhi *= etpow2;
            deltaPhi += coeff0;
            break;
        default:
            XLALPrintError("XLAL Error - %d is an invalid eccentricity order in the enhancement function, must be even and between 2 and 20 or -1\n", EccOrder );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

           REAL8 TildePhi = XLALSimInspiralSpinEccTildePhiSuperAsyLY(etpow2) + deltaPhi;

    return TildePhi;
}

/**
 * @brief Compute the superasymptotic 2.5PN \f$\tilde{\alpha}\f$ enhancement function
 *
 * This function uses the Pade approximant `ref XLALSimInspiralSpinEccPadePolyTildeA` to
 * compute the superasymptotic 2.5PN \f$\tilde{\alpha}\f$ enhancement function.
 *
 * @param[in] etSq Square of the eccentricity \f$e_t^2\f$
 *
 * @return The superasymptotic 2.5PN \f$\tilde{\alpha}\f$ enhancement function
 *
 * @see Equation 170 of Loutrel and Yunes, arXiv:1607.05409
 */
REAL8 XLALSimInspiralSpinEccTildeAlphaSuperAsyLY( REAL8 etSq )
{
    REAL8 sqrt_one_minus_etSq = sqrt(1 - etSq);

    REAL8 Tilde_alpha = (( 67688./963 + sqrt_one_minus_etSq * ( -  4893./107 +
                                        sqrt_one_minus_etSq * ( - 5996./107  +
                                        sqrt_one_minus_etSq * ( 38157./1070  +
                                        sqrt_one_minus_etSq * ( 89699./56175 +
                                        sqrt_one_minus_etSq * ( - 5877./2140 +
                                        sqrt_one_minus_etSq *  944./1605  ))))))/LAL_GAMMA_1_3_MULT_GAMMA_2_3 -
                        21./3424*XLALSimInspiralSpinEccPadePolyTildeA(etSq))/pow(sqrt_one_minus_etSq, 9);
     return Tilde_alpha;
}

/**
 * @brief Compute the superasymptotic 2.5PN \f$\tilde{\beta}\f$ enhancement function
 *
 * @param[in] etSq Square of the eccentricity \f$e_t^2\f$
 *
 * @return The superasymptotic 2.5PN \f$\tilde{\beta}\f$ enhancement function
 *
 * @see Equation 144 of Loutrel and Yunes, arXiv:1607.05409
 */
REAL8 XLALSimInspiralSpinEccTildeBetaSuperAsyLY ( REAL8 etSq )
{
    REAL8 one_minus_etSq = 1 - etSq;

    REAL8 tilde_beta = (1./( LAL_GAMMA_1_3_MULT_GAMMA_2_3 * pow( one_minus_etSq, 9./2) )) * ( 6732800./147627 +
                                    one_minus_etSq * ( - 988160./16403 +
                                    one_minus_etSq * ( 2192128./114821 -
                                    one_minus_etSq * ( 40640./49209 + 93424./31575775*one_minus_etSq))));
    return tilde_beta;
}

/**
 * @brief Compute the superasymptotic 2.5PN \f$\tilde{\gamma}\f$ enhancement function
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 *
 * @return The superasymptotic 2.5PN \f$\tilde{\gamma}\f$ enhancement function
 *
 * @see Equation 145 of Loutrel and Yunes, arXiv:1607.05409
 */
REAL8 XLALSimInspiralSpinEccTildeGammaSuperAsyLY( REAL8 etpow2 )
{
    REAL8 one_minus_etSq = 1 - etpow2;
    REAL8 one_minus_etSq_9_2 = pow(one_minus_etSq, 9./2);

    REAL8 Tilde_gamma = (1./(LAL_GAMMA_1_3_MULT_GAMMA_2_3 * one_minus_etSq_9_2))*( 640./9 +
                                   one_minus_etSq * ( - 512./5 +
                                   one_minus_etSq * ( 6448./175 +
                                   one_minus_etSq * ( - 1024./525 +
                                   one_minus_etSq * ( 14./1375 + 568./56875*one_minus_etSq )))));
    return Tilde_gamma;
}

/**
 * @brief Compute the superasymptotic 2.5PN \f$\tilde{\theta}\f$ enhancement function
 *
 * @param[in] etSq Square of the eccentricity \f$e_t^2\f$
 *
 * @return The superasymptotic 2.5PN \f$\tilde{\theta}\f$ enhancement function
 *
 * @see Equation 153 of Loutrel and Yunes, arXiv:1607.05409
 */
REAL8 XLALSimInspiralSpinEccTildeThetaSuperAsyLY( REAL8 etSq )
{
    REAL8 one_minus_etSq = 1 - etSq;

    REAL8 TildeTheta = (1./(LAL_GAMMA_1_3_MULT_GAMMA_2_3*pow(one_minus_etSq, 9./2)))*( 4672./267 +
                               one_minus_etSq*( - 1452./89 +
                               one_minus_etSq*(119776./46725 -
                               one_minus_etSq*( 266./2225 + 21./4450*one_minus_etSq))));
    return TildeTheta;
}

/**
 * @brief Compute the superasymptotic expression for the 3PN \f$\tilde{\chi}\f$ enhancement function
 *
 * @param[in] etpow2  Square of the eccentricity \f$e_t^2\f$
 *
 * @return The superasymptotic 3PN \f$\tilde{\chi}\f$ enhancement function
 *
 * @see Equation 172 of Loutrel and Yunes, arXiv:1607.05409
 */
REAL8 XLALSimInspiralSpinEccEnhancementTildeChiSuperAsyLY( REAL8 etpow2 )
{
    REAL8 one_minus_etSq = 1 - etpow2;
    REAL8 log_one_minus_etSq = log(one_minus_etSq);
    REAL8 one_minus_etSq_pow5 = pow( one_minus_etSq, 5);

    REAL8 coeff0 = 1./512*(35583 - 6930*LAL_GAMMA - 27720*LAL_LN2 - 3465*LAL_LN3 -
                    10395*log_one_minus_etSq);
    REAL8 coeff1 = 1./512*(-51359 + 9310*LAL_GAMMA + 37240*LAL_LN2 + 4655*LAL_LN3 +
                    13965*log_one_minus_etSq);
    REAL8 coeff2 = 1./512*(47481*512./1280 - 3030*LAL_GAMMA - 12120*LAL_LN2 - 1519*LAL_LN3 +
                    5*LAL_LN3*512./640 - 4545*log_one_minus_etSq);
    REAL8 coeff3 = - 53091./22400 + 5*LAL_LN3/1800 +  1./512*( 138* LAL_GAMMA  + 552*LAL_LN2 +
                    3041./45*LAL_LN3 +  207*log_one_minus_etSq);
    REAL8 coeff4 = - 956569./68992000 + LAL_LN2/700 - 6269*LAL_LN3/4435200 + 7061*5*LAL_LN3/22176000
                    - ( 8*LAL_LN2 + LAL_LN3 )/5600;
    REAL8 coeff5 = 15822507./22422400000 - 553*LAL_LN2/6435000 - 553*LAL_LN3/51480000 + 553*(8*LAL_LN2 + LAL_LN3)/51480000;


    REAL8 TildeChi =  1./one_minus_etSq_pow5* ( coeff0 +
                                           one_minus_etSq * (coeff1 +
                                           one_minus_etSq * (coeff2 +
                                           one_minus_etSq * (coeff3 + one_minus_etSq * coeff4 )))) + coeff5;
    return TildeChi;
}

/**
 * @brief Compute the hyperasymptotic expression for the 3PN \f$\tilde{\chi}\f$ enhancement function
 *
 * This function evaluates the hyperasymptotic 3PN \f$\tilde{\chi}\f$ enhancement function by adding
 * a polynomial correction in \f$e_t^2\f$ to `ref XLALSimInspiralSpinEccEnhancementTildeChiSuperAsyLY`.
 *
 * @param[in] etpow2  Square of the eccentricity \f$e_t^2\f$
 * @param[in] EccOrder  Order of eccentricity for the hyperasymptotic correction;
 *                      even number between 2 and 10, or -1 for the maximum available order.
 *                      Values up to 20 are accepted but do not change the result.
 *
 * @return The hyperasymptotic 3PN \f$\tilde{\chi}\f$ enhancement function
 *
 * @see Equation 187 of Loutrel and Yunes, arXiv:1607.05409
 */
REAL8 XLALSimInspiralSpinEccEnhancementTildeChiHyperAsyLY( REAL8 etpow2, LALSimInspiralSpinEccPolyEccOrder EccOrder )
{
    REAL8 coeff0 = - XLALSimInspiralSpinEccEnhancementTildeChiSuperAsyLY(0.);
    REAL8 coeff2 = - 3744713821./68992000 + 389.* LAL_GAMMA/32 + 1007.*LAL_LN2/32 + 2965.*LAL_LN3/128;
    REAL8 coeff4 = - 18889494241./68992000 + 3577.* LAL_GAMMA/64 + 27699.* LAL_LN2/64 - 27245.*LAL_LN3/512;
    REAL8 coeff6 = - 60790558061./68992000 + 43049.* LAL_GAMMA/256 - 1804397.*LAL_LN2/2304 + 1946309.*LAL_LN3/8192 +
                    48828125.*LAL_LN5/73728;
    REAL8 coeff8 =  - 151724371031./68992000 + 102005.* LAL_GAMMA/256 + 15307525*LAL_LN2/2304 + 184515253*LAL_LN3/65536 -
                    2099609375.*LAL_LN5/589824;
    REAL8 coeff10 = -160776312313./34496000 + 207311.*LAL_GAMMA/256  - 107117837.*LAL_LN2/6400 - 815824256293.*LAL_LN3/52428800 +
                    156103515625*LAL_LN5/18874368 + 4747561509943 *LAL_LN7/471859200;

    REAL8 deltaTildeChi = 0;
    switch (EccOrder){
        case(-1):
        case(20):
        case(18):
        case(16):
        case(14):
        case(12):
        case(10):
            deltaTildeChi = coeff10;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(8):
            deltaTildeChi *= etpow2;
            deltaTildeChi += coeff8;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(6):
            deltaTildeChi *= etpow2;
            deltaTildeChi += coeff6;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(4):
            deltaTildeChi *= etpow2;
            deltaTildeChi += coeff4;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(2):
            deltaTildeChi *= etpow2;
	        deltaTildeChi += coeff2;
	        deltaTildeChi *= etpow2;
            deltaTildeChi += coeff0;
            break;
       default:
            XLALPrintError("XLAL Error - %d is an invalid eccentricity order in the enhancement function, must be even and between 2 and 20 or -1\n", EccOrder );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    REAL8 TildeChiHyper = XLALSimInspiralSpinEccEnhancementTildeChiSuperAsyLY( etpow2 ) + deltaTildeChi;
    return TildeChiHyper;
}

/**
 * @brief Compute the superasymptotic expression for the 2.5PN \f$\tilde{\psi}\f$ enhancement function
 *
 * This function evaluates the 2.5PN \f$\tilde{\psi}\f$ enhancement function as a weighted sum of
 * `ref XLALSimInspiralSpinEccTildeAlphaSuperAsyLY`, `ref XLALSimInspiralSpinEccTildeBetaSuperAsyLY`, and
 * `ref XLALSimInspiralSpinEccTildeGammaSuperAsyLY`.
 *
 * @param[in] etSq  Square of the eccentricity \f$e_t^2\f$
 *
 * @return The superasymptotic 2.5PN \f$\tilde{\psi}\f$ enhancement function
 *
 * @see Equation A9a of Phukon et al., arXiv:2504.20543 with tilded versions of the enhancement functions
 */
REAL8 XLALSimInspiralSpinEccTildePsiSuperAsyLY( REAL8 etSq )
{
     REAL8 TildePsi = 13696./8191 * XLALSimInspiralSpinEccTildeAlphaSuperAsyLY (etSq) -
                16403./24573 * XLALSimInspiralSpinEccTildeBetaSuperAsyLY( etSq ) -
                112./24573 * XLALSimInspiralSpinEccTildeGammaSuperAsyLY( etSq);
    return TildePsi;
}

/**
 * @brief Compute the hyperasymptotic expression for the 2.5PN \f$\tilde{\psi}\f$ enhancement function
 *
 * This function evaluates the hyperasymptotic 2.5PN \f$\tilde{\psi}\f$ enhancement function.
 * The hyperasymptotic result augments the superasymptotic value `ref XLALSimInspiralSpinEccTildePsiSuperAsyLY`
 * with a polynomial correction in \f$e_t^2\f$, \f$\Delta(e)\f$.
 *
 * @param[in] etpow2  Square of the eccentricity \f$e_t^2\f$
 * @param[in] EccOrder  Order of eccentricity for the hyperasymptotic correction;
 *                      even number between 2 and 20; 20 or -1 for the maximum available order
 *
 * @return The hyperasymptotic 2.5PN \f$\tilde{\psi}\f$ enhancement function
 *
 * @see Equation 183 of Loutrel and Yunes, arXiv:1607.05409
 */
REAL8 XLALSimInspiralSpinEccTildePsiHyperAsyLY( REAL8 etpow2, LALSimInspiralSpinEccPolyEccOrder EccOrder )
{
    REAL8 inverse_gamma_mult = 1./LAL_GAMMA_1_3_MULT_GAMMA_2_3;
    REAL8 coeff0 = 1 - XLALSimInspiralSpinEccTildePsiSuperAsyLY(0.);
    REAL8 coeff2 = 102536./8191 - (698208327368./15373483125)*inverse_gamma_mult;
    REAL8 coeff4 = 27975523./524224 - (2976133354982./15373483125)*inverse_gamma_mult;
    REAL8 coeff6 = 709642057./4718016 - (8388221641661./15373483125)*inverse_gamma_mult;
    REAL8 coeff8 = 203853989947./603906048 - (4302911627633./3513939000)*inverse_gamma_mult;
    REAL8 coeff10 = 4944184758677./7548825600 - (83488815643601./35139390000)*inverse_gamma_mult;
    REAL8 coeff12 = 6658083547039409./5797498060800 - (17744668624161./4259320000)*inverse_gamma_mult;
    REAL8 coeff14 = 35397103550602159./18938493665280 - (404305552234341./59630480000)*inverse_gamma_mult;
    REAL8 coeff16 = 179071486944184743991./62334699149721600 - (764822533448511./73391360000)*inverse_gamma_mult;
    REAL8 coeff18 = 7457576214411997508197./1767188720894607360 - (1348218095132281./88069632000)*inverse_gamma_mult;
    REAL8 coeff20 =  337929898617545561.543145703/56550039068627.435520000  -
                      (38182917211753667./1761392640000)*inverse_gamma_mult;

    REAL8 deltaTildePsi = 0;

    switch(EccOrder)
    {
        case(-1):
        case(20):
            deltaTildePsi = coeff20;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(18):
            deltaTildePsi *= etpow2;
            deltaTildePsi += coeff18;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(16):
            deltaTildePsi *= etpow2;
            deltaTildePsi += coeff16;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(14):
            deltaTildePsi *= etpow2;
            deltaTildePsi += coeff14;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(12):
            deltaTildePsi *= etpow2;
            deltaTildePsi += coeff12;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(10):
            deltaTildePsi *= etpow2;
            deltaTildePsi += coeff10;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(8):
            deltaTildePsi *= etpow2;
            deltaTildePsi += coeff8;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(6):
            deltaTildePsi *= etpow2;
            deltaTildePsi += coeff6;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(4):
            deltaTildePsi *= etpow2;
            deltaTildePsi += coeff4;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(2):
            deltaTildePsi *= etpow2;
	        deltaTildePsi += coeff2;
	        deltaTildePsi *= etpow2;
            deltaTildePsi += coeff0;
            break;
      default:
            XLALPrintError("XLAL Error - %d is an invalid eccentricity order in the enhancement function, must be even and between 2 and 20 or -1\n", EccOrder );
            XLAL_ERROR(XLAL_EINVAL);
            break;

    }
    REAL8 TildePsi = XLALSimInspiralSpinEccTildePsiSuperAsyLY( etpow2 ) + deltaTildePsi;
    return TildePsi;
}

/**
 * @brief Compute the superasymptotic expression for the 2.5PN \f$\tilde{\zeta}\f$ enhancement function
 *
 * This function calculates the superasymptotic expression for the 2.5PN \f$\tilde{\zeta}\f$ enhancement
 * term using a  weighted combination of `ref XLALSimInspiralSpinEccTildeThetaSuperAsyLY`,
 * `ref XLALSimInspiralSpinEccTildeBetaSuperAsyLY` and `ref XLALSimInspiralSpinEccTildeGammaSuperAsyLY`.
 *
 * @param[in] etSq Square of the eccentricity \f$e_t^2\f$
 *
 * @return The 2.5PN \f$\tilde{\zeta}\f$ enhancement function
 *
 * @see This function implements Equation A9b of Phukon et al., arXiv:2504.20543, with enhancement
 * functions replaced by their tilded versions.
 */
REAL8 XLALSimInspiralSpinEccTildeZetaSuperAsyLY ( REAL8 etSq)
{
    REAL8 TildeZeta = - (1424./4081) * XLALSimInspiralSpinEccTildeThetaSuperAsyLY( etSq )  +
                    (16403./12243) * XLALSimInspiralSpinEccTildeBetaSuperAsyLY( etSq ) +
                    (16./1749) * XLALSimInspiralSpinEccTildeGammaSuperAsyLY( etSq);
    return TildeZeta;
}

/**
 * @brief Compute the hyperasymptotic expression for the 2.5PN \f$\tilde{\zeta}\f$ enhancement function
 *
 * This function calculates the hyperasymptotic expression for the 2.5PN \f$\tilde{\zeta}\f$ enhancement
 * function.
 * The hyperasymptotic expression improves upon the superasymptotic result by including
 * higher-order corrections in eccentricity: \f[ \tilde{\zeta}_\mathrm{hyper}(e) = \tilde{\zeta}_\mathrm{super}(e) + \delta\tilde{\zeta}(e) \f]
 * The superasymptotic \f$\tilde{\zeta}_\mathrm{super}(e)\f$ is computed by `ref XLALSimInspiralSpinEccTildeZetaSuperAsyLY`
 * and \f$\delta\tilde{\zeta}(e)\f$ is a polynomial correction in \f$e_t^2\f$ up to the specified eccentricity order.
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 * @param[in] EccOrder Order of eccentricity correction to control truncation of the polynomial expansion (even number between 2 and 20; -1 or 20 for maximum available order)
 *
 * @return The hyperasymptotic 2.5PN \f$\tilde{\zeta}\f$ enhancement function
 *
 * @see This function implements Equation 185 of Loutrel and Yunes, arXiv:1607.05409
 */
REAL8 XLALSimInspiralSpinEccTildeZetaHyperAsyLY ( REAL8 etpow2, LALSimInspiralSpinEccPolyEccOrder EccOrder)
{
    REAL8 inverse_gammas = 1./LAL_GAMMA_1_3_MULT_GAMMA_2_3;
    REAL8 coeff0 = 1 - XLALSimInspiralSpinEccTildeZetaSuperAsyLY(0.);
    REAL8 coeff2 = 102371./8162 - 348496717732./7659526875 * inverse_gammas;
    REAL8 coeff4 = 14250725./261184 - 1516042558243./7659526875 * inverse_gammas;
    REAL8 coeff6 = 722230667./4701312 - 8537054301053./15319053750 * inverse_gammas;
    REAL8 coeff8 = 102744533069./300883968 - 4337433374609./3501498000 * inverse_gammas;
    REAL8 coeff10 = 9843430194463./15044198400 - 83109477430673./35014980000 * inverse_gammas;
    REAL8 coeff12 = 43605309737981./38513147904 - 52296280129859./12732720000 * inverse_gammas;
    REAL8 coeff14 = 19069924628449467./10484134707200 - 392069990496293./59419360000 * inverse_gammas;
    REAL8 coeff16 =  600336343160.521814159/217399017.288499200 - 732589863363103./73131520000* inverse_gammas;
    REAL8 coeff18 = 28244435149543.337941721/7043728160.147374080 - 3829628554688939./263273472000 * inverse_gammas;
    REAL8 coeff20 = 158264167343831506.620212273/28174912640589.496320000 - 35764742675504291./1755156480000 * inverse_gammas;

    REAL8 deltaTildeZeta = 0;

    switch( EccOrder )
    {
        case (-1):
        case(20):
            deltaTildeZeta = coeff20;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(18):
            deltaTildeZeta *= etpow2;
            deltaTildeZeta += coeff18;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(16):
            deltaTildeZeta *= etpow2;
            deltaTildeZeta += coeff16;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(14):
            deltaTildeZeta *= etpow2;
            deltaTildeZeta += coeff14;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(12):
            deltaTildeZeta *= etpow2;
            deltaTildeZeta += coeff12;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(10):
            deltaTildeZeta *= etpow2;
            deltaTildeZeta += coeff10;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(8):
            deltaTildeZeta *= etpow2;
            deltaTildeZeta += coeff8;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(6):
            deltaTildeZeta *= etpow2;
            deltaTildeZeta += coeff6;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(4):
            deltaTildeZeta *= etpow2;
            deltaTildeZeta += coeff4;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
        case(2):
            deltaTildeZeta *= etpow2;
	        deltaTildeZeta += coeff2;
	        deltaTildeZeta *= etpow2;
            deltaTildeZeta += coeff0;
            break;
	default:
            XLALPrintError("XLAL Error - %d is an invalid eccentricity order in the enhancement function, must be even and between 2 and 20 or -1\n", EccOrder );
            XLAL_ERROR(XLAL_EINVAL);
            break;

    }
    REAL8 TildeZeta = XLALSimInspiralSpinEccTildeZetaSuperAsyLY(etpow2) + deltaTildeZeta;
    return TildeZeta;
}

/**
 * @brief Compute the superasymptotic expression for the 3PN \f$\tilde{\kappa}\f$ enhancement function
 *
 * This function calculates the superasymptotic expression for the 3PN \f$\tilde{\kappa}\f$enhancement
 * function using a weighted combination of `ref XLALSimInspiralSpinEccEnhancementTildeChiSuperAsyLY` and
 * `ref XLALSimInspiralSpinEccTildeEnhancementFlux`.
 *
 * @param[in] etpow2 Square of the eccentricity \f$e_t^2\f$
 *
 * @return The superasymptotic 3PN \f$\tilde{\kappa}\f$ enhancement function \f$\tilde{\kappa}\f$
 *
 * @see This function implements the general expression given in Equation A9c of Phukon et al., arXiv:2504.20543,
 * with enhancement functions replaced by their tilded versions.
 */
REAL8 XLALSimInspiralSpinEccEnhancementTildeKappaSuperAsyLY( REAL8 etpow2 )
{
    REAL8 TildeKappa = 59920./116761*XLALSimInspiralSpinEccEnhancementTildeChiSuperAsyLY (etpow2)  + XLALSimInspiralSpinEccTildeEnhancementFlux( etpow2 );
    return TildeKappa;
}

/**
 * @brief Compute the hyperasymptotic expression for the \f$\tilde{\kappa}\f$ enhancement function
 *
 * This function evaluates the hyperasymptotic \f$\tilde{\kappa}\f$ enhancement function. It follows the general
 * expression of \f$\kappa\f$ (Eq. A9c of Phukon et al., arXiv:2504.20543), with enhancement
 * functions replaced by their tilded counterparts. The hyperasymptotic form combines
 * the superasymptotic value with an eccentricity-order correction, i.e.
 * \f$E_{\mathrm{hyper}}(e) = E_{\mathrm{super}}(e) + \Delta(e)\f$. The function uses `\ref XLALSimInspiralSpinEccEnhancementTildeChiHyperAsyLY`
 * and `\ref XLALSimInspiralSpinEccTildeEnhancementFlux` to compute the necessary components.
 *
 * @param[in] etpow2  Square of the eccentricity (\f$e_t^2\f$)
 * @param[in] EccOrder  Order of eccentricity for the hyperasymptotic correction;
 *                      must be even between 2 and 10 (values up to 20 are accepted
 *                      but do not change the result)
 *
 * @return The hyperasymptotic 3PN \f$\tilde{\kappa}\f$ enhancement function
 *
 * @see Equation A9c of Phukon et al., arXiv:2504.20543.
 */
REAL8 XLALSimInspiralSpinEccEnhancementTildeKappaHyperAsyLY( REAL8 etpow2, LALSimInspiralSpinEccPolyEccOrder EccOrder )
{

    REAL8 TildeKappaHyper =  59920./116761*XLALSimInspiralSpinEccEnhancementTildeChiHyperAsyLY (etpow2, EccOrder) + XLALSimInspiralSpinEccTildeEnhancementFlux( etpow2 );
    return TildeKappaHyper;

}

/**
 * @brief Internal function to compute the energy variable and its time derivative
 *
 * This function calculates the energy variable \f$\epsilon\f$ (proportional to \f$-E\f$) and its time derivative.
 * The time derivative is computed using the chain rule, by taking partial derivatives with respect to (\f$x\f$) and (\f$e_t^2\f$).
 *
 *
 * @param[out] epsilon      Energy variable (\f$\epsilon\f$)
 * @param[out] depsilon_dt  Time derivative of the energy variable
 * @param[in]  xbar         Velocity parameter (\f$\bar{x} = x/(1 - e^2)\f$)
 * @param[in]  etSq         Square of eccentricity (\f$e_t^2\f$)
 * @param[in]  domega_dt    Time derivative of orbital angular velocity
 * @param[in]  detSq_dt     Time derivative of square of eccentricity
 * @param[in]  ak           Structure containing static parameters
 *
 * @return XLAL_SUCCESS on successful completion
 *
 * @see Eqs. C1a, C2a-C2d and C4a-C4c of Phukon et al., arXiv:2504.20543  for the energy variable and its time derivative, respectively.
 */
static int evaluate_depsilon_dt (
    REAL8 *epsilon,
    REAL8 *depsilon_dt,
    REAL8 xbar,
    REAL8 etSq,
    REAL8 domega_dt,
    REAL8 detSq_dt,
    expnCoeffsSpinEcc *ak
    )
{
    REAL8 one_minus_etsq = 1 - etSq;
    REAL8 sqrt_one_minus_etsq = sqrt(one_minus_etsq);
    REAL8 one_minus_etsq_3by2 = one_minus_etsq*sqrt_one_minus_etsq;

    REAL8 mtotal = ak->m;
    REAL8 Sqrtxbar = sqrt(xbar);
    REAL8 v = Sqrtxbar*sqrt_one_minus_etsq;
    REAL8 PI_Sq = LAL_PI*LAL_PI;

    REAL8 eta = ak->eta;

    /* Newtonian term in the energy variable, epsilon_bar
     * Equation C1a of Phukon et al., arXiv:2504.20543
     *
     */
    REAL8 epsilon_0PN_coeff = 1;

    /* 1PN term in the energy variable, epsilon_bar
     * Expression is given in Equation C1a, C2a of Phukon et al., arXiv:2504.20543
     *
     */
    REAL8 epsilon_1PN_coeff = (- 3/4. - 5/4.*etSq - eta/12*one_minus_etsq);

    /* 2PN non-spinning term in the energy variable, epsilon_bar
     * Expression is given in Eqs. C1a, C2c (excluding gamma1 term), Phukon et al., arXiv:2504.20543
     *
     */
    REAL8 epsilon_2PN_coeff = ((- 67/8. +
                                        eta*(35/8. - eta/24)) +
                                etSq * (( -19/4. +
                                        eta*(21/4. + eta/12.) ) +
                                etSq * (5/8. -
                                        eta*(5/8. + eta/24)))) -
                                one_minus_etsq_3by2*(2*eta - 5);

    /* 3PN term in the energy variable, epsilon_bar
     * Expression is given in Eqs. C1a, C2d, Phukon et al., arXiv:2504.20543
     *
     */
    REAL8 epsilon_3PN_coeff = (-835/64. + eta*((18319/192. - 41/16.* PI_Sq) -
                                            eta*( 169/32. + 35/5184.*eta ))) +
                                etSq *  (( -3703/64. +  eta*( (21235/192. - 41/64.*PI_Sq) +
                                                        eta*( - 7733/288.  + 35/1728.*eta))) +
                                etSq * (( 103/64. - eta*(547/192.  +
                                                    eta*(1355/288. + 35/1728.*eta))) +
                                etSq * ( 185/192. + eta*(75/64.  +
                                                    eta*(25/288. + 35/5184.*eta))))) +
                                sqrt_one_minus_etsq * ( (5/2. + eta*((-641/18. + 41/96.*PI_Sq) +
                                                                eta*11/3.)) +
                                                       etSq * ( (-35 + eta*((394/9. - 41/96.*PI_Sq) - eta/3)) +
                                                       etSq * ( 5./2 + eta*(23/6. - 10/3.*eta))));

    /* 2PN non-spinning term in partial derivative of epsilon w.r.t. eccentricity squared
     * Expression is given in Equation C4c of Phukon et al., arXiv:2504.20543 as coefficeint of \bar{x}^3 excluding gamma_1 term
     *
     */
    REAL8 epsilon_2PN_coeff_detsq = -(-5 + 43*sqrt_one_minus_etsq + eta*(2 - 28*sqrt_one_minus_etsq) +
                                            etSq*( 10 + 7*sqrt_one_minus_etsq - 4*eta*(1 + 2*sqrt_one_minus_etsq) +
                                            etSq*(-5 + 2*eta)))/(2*sqrt_one_minus_etsq);

    /* 3PN coefficient in partial derivative of epsilon w.r.t. eccentricity squared
     * Expression is given in Equation C4c of Phukon et al., arXiv:2504.20543 as coefficeint of \bar{x}^4
     *
     */
    REAL8 epsilon_3PN_coeff_detsq = ((-16560 - 55872*sqrt_one_minus_etsq) + eta*((-26064 + 228576*sqrt_one_minus_etsq -
                                            PI_Sq*(-369 + 4797*sqrt_one_minus_etsq)) + eta*(5088 - 24592*sqrt_one_minus_etsq)) +
                                        etSq*(((-10800 - 64800*sqrt_one_minus_etsq) + eta*( (68304 + 124128*sqrt_one_minus_etsq) -
                                            738*PI_Sq*(1 + sqrt_one_minus_etsq)- eta*(9216 + 36352*sqrt_one_minus_etsq)) ) +
                                    etSq*((28080 + 2592*sqrt_one_minus_etsq) + eta*( -41136 + 384*sqrt_one_minus_etsq + 369*PI_Sq -
                                                                             eta*(-3168 + 2560*sqrt_one_minus_etsq)) +
                                    etSq*(-720 - eta*(1104 - 960*eta))) ))/(576.*sqrt_one_minus_etsq);


    /* variable for derivative of epsilon w.r.t. x */
    REAL8  depsilon_dx = 0;

    /* variable for derivarive of epsilon w.r.t. eccentricity squared */
    REAL8  depsilon_detsq = 0;

    /* variable for epsilon (with Newtonian term scaled out) */
    REAL8  epsilon_value = 0;

    LALSimInspiralSpinOrder spinO = ak->SpinOrder;
    INT4 phaseO = ak->PhaseOrder;

    switch (phaseO)
    {
        case(-1):
        case(6):
            /* 3PN coeffecient in derivative of epsilon w.r.t. x from
             * Equation C4b of Phukon et al., arXiv:2504.20543
             */
	    depsilon_dx = 4 * epsilon_3PN_coeff;

            depsilon_detsq = epsilon_3PN_coeff_detsq;
	    epsilon_value = epsilon_3PN_coeff;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
	__attribute__ ((fallthrough));
#endif
        case(5):
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
        __attribute__ ((fallthrough));
#endif
        case(4):
            depsilon_dx *= xbar;
            /* 2PN non-spinning coeffecient in derivative of epsilon w.r.t. x
             * from Equation C4b of Phukon et al., arXiv:2504.20543
             */
            depsilon_dx += 3 * epsilon_2PN_coeff;

            depsilon_detsq *= xbar;
            depsilon_detsq += epsilon_2PN_coeff_detsq;

	    epsilon_value *= xbar;
	    epsilon_value += epsilon_2PN_coeff;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
	__attribute__ ((fallthrough));
#endif
        case(3):
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
	__attribute__ ((fallthrough));
#endif
        case(2):
	    depsilon_dx *= xbar;
            /* 1PN coeffecient in derivative of epsilon w.r.t. x from
             * Equation C4b of Phukon et al., arXiv:2504.20543
             */
            depsilon_dx += 2 * epsilon_1PN_coeff;

	    depsilon_detsq *= xbar;
            /* 1PN term (coefficient of \bar{x}^2) in derivative of epsilon w.r.t. eccentricity squared
             * from Equation C4c of Phukon et al., arXiv:2504.20543
             */
            depsilon_detsq -= 2;
	    depsilon_detsq *= xbar*xbar;

	    epsilon_value *= xbar;
	    epsilon_value += epsilon_1PN_coeff;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
	__attribute__ ((fallthrough));
#endif
        case(1):
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
	__attribute__ ((fallthrough));
#endif
        case(0):
	    depsilon_dx *= xbar;
            depsilon_dx += epsilon_0PN_coeff;

	    epsilon_value *= xbar;
	    epsilon_value += epsilon_0PN_coeff;
	    break;
        default:
            XLALPrintError("XLAL Error - %s: - %d is an invalid phase order\n", __func__, phaseO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    REAL8 Beta43 = ak->beta43;
    REAL8 Gamma1 = ak->gamma1;

    /* 1.5PN term in  the energy variable
     * Expression is given in Eqs. C1a, C2b, Phukon et al., arXiv:2504.20543
     *
     */
    REAL8 epsilon_1p5PN_spin_coeff = 2/3.*Beta43;

    /* 2PN spin term in the energy variable
     * Expression is from Eqs. C1a, C2c, Phukon et al., arXiv:2504.20543
     *
     */
    REAL8 epsilon_2PN_spin_coeff = - Gamma1;

    /* variable for spin terms in derivarive of epsilon w.r.t. x */
    REAL8 depsilon_dx_spin = 0;

    /* variable for spin terms in derivarive of epsilon w.r.t. eccentricity squared */
    REAL8 depsilon_detsq_spin = 0;

    /* variable for spin terms in epsilon */
    REAL8 epsilon_spin = 0;
    switch(spinO)
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            /* 2PN spin term in derivarive of epsilon w.r.t. x from
             * Equation C4b, Phukon et al., arXiv:2504.20543
             */
            depsilon_dx_spin = 3 * xbar * epsilon_2PN_spin_coeff;

            /* 2PN spin term (term with gamma_1 in coefficient in \bar{x}^3)
             * in derivarive of epsilon w.r.t. eccentricity squared from
             * Equation C4c, Phukon et al., arXiv:2504.20543
             */
            depsilon_detsq_spin = -2*xbar*Gamma1;

	    epsilon_spin = xbar * epsilon_2PN_spin_coeff;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
	__attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            /* 1.5 PN term in derivarive of epsilon w.r.t. x from
             * Equation C4b of Phukon et al., arXiv:2504.20543
            */
            depsilon_dx_spin += Sqrtxbar * 2.5 * epsilon_1p5PN_spin_coeff;
	    depsilon_dx_spin *= xbar;

            /* 1.5 PN term in derivarive of epsilon w.r.t. eccentricity squared from
             * Equation C4c (coeffecient of \bar{x}^{5/2}) of Phukon et al., arXiv:2504.20543
             */
            depsilon_detsq_spin += Sqrtxbar * Beta43;
	    depsilon_detsq_spin *= xbar*xbar;

            epsilon_spin += Sqrtxbar * epsilon_1p5PN_spin_coeff;
	    epsilon_spin *= xbar;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
	__attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error in function - %s: Invalid spin PN order %d\n",__func__, spinO );
            XLAL_ERROR(XLAL_EINVAL);
        break;
    }

    depsilon_dx += depsilon_dx_spin;
    depsilon_detsq += depsilon_detsq_spin;
    epsilon_value += epsilon_spin;
    *epsilon = epsilon_value;

    REAL8 dx_dt = 2./3 * mtotal/v * domega_dt;
    *depsilon_dt = dx_dt*depsilon_dx + detSq_dt*depsilon_detsq;
    return XLAL_SUCCESS;
}

/**
 * @brief Internal function to compute the scaled orbital angular momentum magnitude and its time derivative
 *
 * This function calculates the scaled magnitude of the orbital angular momentum (\f$L_{scaled}\f$)
 * and its time derivative for use in the post-Newtonian evolution equations of eccentric binaries.
 * The scaling is chosen to simplify the analysis of the sign and behavior of the angular momentum.
 *
 * The scaled angular momentum is defined as angular momentum magnitude divided by the Newtonian angular momentum: \f$L_{scaled} = L/L_N\f$
 *
 * The time derivative is computed using the chain rule with partial derivatives with respect to (\f$x\f$) and (\f$e_t^2\f$).
 *
 * @param[out] Lscaled      Scaled orbital angular momentum magnitude
 * @param[out] dLscaled_dt  Time derivative of scaled orbital angular momentum
 * @param[in]  xbar         PN expansion parameter (\f$\bar{x} = x/(1 - e^2)\f$)
 * @param[in]  etSq         Square of eccentricity (\f$e_t^2\f$)
 * @param[in]  s1x,s1y,s1z  Components of dimensionless spin vector for body 1
 * @param[in]  s2x,s2y,s2z  Components of dimensionless spin vector for body 2
 * @param[in]  lhatx,lhaty,lhatz  Unit vector components of orbital angular momentum
 * @param[in]  domega_dt    Time derivative of orbital angular velocity
 * @param[in]  detSq_dt     Time derivative of square of eccentricity
 * @param[in]  ak           Structure containing static parameters
 *
 * @return XLAL_SUCCESS on successful completion
 *
 * @see Eqs. C1b, C3a-C3d and C5a-C5c of Phukon et al., arXiv:2504.20543 for the scaled angular momentum and it's time derivative
 */
static int evaluate_dL_scaled_dt ( REAL8 *Lscaled, REAL8 *dLscaled_dt,  REAL8 xbar, REAL8 etSq,
                                  REAL8 s1x, REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z,
                                  REAL8 lhatx, REAL8 lhaty, REAL8 lhatz, REAL8 domega_dt, REAL8 detSq_dt,
                                  expnCoeffsSpinEcc *ak )
{
    REAL8 one_minus_etsq = 1 - etSq;
    REAL8 sqrt_one_minus_etsq = sqrt(one_minus_etsq);
    REAL8 mtotal = ak->m;
    REAL8 Sqrtxbar = sqrt(xbar);
    REAL8 v = Sqrtxbar*sqrt_one_minus_etsq;
    REAL8 PI_Sq = LAL_PI*LAL_PI;
    REAL8 eta = ak->eta;

    /* Newtonian coefficient in the expression of the scaled angular momentum variable
     * First term in Equation C1b of Phukon et al., arXiv:2504.20543
     * */
    REAL8 Lscaled0PNCoeff = 1;

    /* 1PN coefficient in the expression of the scaled angular momentum variable
     * Eqs. C1b, C3a of Phukon et al., arXiv:2504.20543
     */
    REAL8 Lscaled1PNCoeff = (3/2.*one_minus_etsq  + eta*(1/6. + 5/6.*etSq));

    /* 2PN coefficient (non-spinning parts) in the expression of the scaled angular momentum variable
     * Eqs. C1b, C3c (except the gamma_1 term) of Phukon et al., arXiv:2504.20543
     * */
    REAL8 Lscaled2PNCoeff = 47/8. - eta*(27/8. - eta/24) +
                            sqrt_one_minus_etsq*(-5/2. + eta + (-5 + 2*eta)*etSq) +
                            etSq*(
                                  21/4. - eta*(5/6. + 3/4.*eta ) +
                             etSq*(
                                  11/8. - eta*(73/24. - 5/24.*eta)));

    /* 3PN coefficient in the expression of the scaled angular momentum
     * Eqs. C1b, C3d of Phukon et al., arXiv:2504.20543
     * */
    REAL8 Lscaled3PNCoeff = (155/16. + eta*( (-3151/48. + 123/64.*PI_Sq) + eta*(25/8. + 7/1296.*eta))) +
                                etSq*(
                                        (1227/16. + eta*(
                                                         (119/128.*PI_Sq - 265/2.) +
                                                         eta*(787/36. + eta*95/432.))) +
                                etSq*(
                                        (169/16. + eta*(
                                                         - 115/6. +
                                                         eta*(109/9. + 127/432.*eta))) +
                                etSq*(
                                        (-13/48. + eta*(
                                                         283/48. -
                                                         eta*(71/36. + 25/1296.*eta)))))) +
                                sqrt_one_minus_etsq *  (
                                                        (- 5/4. + eta*((641/36. - 41/192.*PI_Sq)  - 11/6.* eta ))  +
                                                        etSq*( (-135/4. + eta*(
                                                                    (2359/36.  - 41/96.*PI_Sq) - 34/3.*eta)) +
                                                        etSq*( 5 - eta*( 14/3. + eta/3. ))));

    /* Newtonian coefficient in partial derivative of angular momentum variable w.r.t. eccentricity squared
     * Coefficient of \bar{x} in Equation C5c of Phukon et al., arXiv:2504.20543
     * */
    REAL8 Lscaled0PNCoeff_detsq = -1;

    /* 1PN coefficient in partial derivative of angular momentum variable w.r.t. eccentricity squared
     * Coefficient of \bar{x}^2 in Equation C5c of Phukon et al., arXiv:2504.20543
     * */
    REAL8 Lscaled1PNCoeff_detsq = (-9 + etSq*(9 - 5*eta) + 11*eta)/6.;

    /* 2PN coefficient in partial derivative of angular momentum variable w.r.t. eccentricity squared
     * Coefficient of \bar{x}^3 in Equation C5c of Phukon et al., arXiv:2504.20543
     * */
    REAL8 Lscaled2PNCoeff_detsq = ( eta*( (-283 + 144*sqrt_one_minus_etsq) - 33*eta) + etSq*((258 - 2*eta*(156 - eta)) + etSq*(-33 + eta*(73 - 5*eta))  ) +
                                            (675 - 360*sqrt_one_minus_etsq))/24.;

    /* 3PN coefficient in partial derivative of angular momentum variable w.r.t. eccentricity squared
     * Coefficient of \bar{x}^4 in Equation C5c of Phukon et al., arXiv:2504.20543
     * */
    REAL8 Lscaled3PNCoeff_detsq = (etSq*(( (2823336 - 492480*sqrt_one_minus_etsq) +  eta*(((-4916160 + 1165248*sqrt_one_minus_etsq) + PI_Sq*(28917 - 8856*sqrt_one_minus_etsq)) +
                                                                                                          eta*(( 1182240 - 248832*sqrt_one_minus_etsq ) +
                                                                                                           eta * 19032  ))) +
                                    etSq*((92664 + eta*(168048 + eta*(2880 + 1848*eta))) +
                                    etSq*((2808 + eta*(- 61128 + eta*(20448 + 200*eta))))  )) -
                                       ( eta*((6150600 - 2097216*sqrt_one_minus_etsq + PI_Sq*(-118908 + 17712*sqrt_one_minus_etsq)) +
                                        eta*( (-615312 + 311040*sqrt_one_minus_etsq) - eta*4840 ) ) +   (-2092392 + 751680*sqrt_one_minus_etsq)))/10368.;

    /* variable for derivative w.r.t x*/
    REAL8 dLscaled_dx = 0;

    /* variable for derivative w.r.t. eccentricity squared*/
    REAL8 dLscaled_detsq = 0;

    /* variable for scaled orbital angular momentum */
    REAL8 Lscaled_value = 0;


    LALSimInspiralSpinOrder spinO = ak->SpinOrder;
    INT4 phaseO = ak->PhaseOrder;

    switch(phaseO)
    {
        case(-1):
        case(6):
            /* 3PN term for derivarive of scaled orbital angular momentum w.r.t. x
             * Equation C5b of Phukon et al., arXiv:2504.20543
             */
	    dLscaled_dx = 5*Lscaled3PNCoeff;

            dLscaled_detsq = Lscaled3PNCoeff_detsq;
            Lscaled_value = Lscaled3PNCoeff;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
        __attribute__ ((fallthrough));
#endif
        case(5):
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
         __attribute__ ((fallthrough));
#endif
       case(4):
            dLscaled_dx *= xbar;

            /* 2PN term non-spinning term for derivarive of scaled orbital angular momentum w.r.t. x
             * Equation C5b of Phukon et al., arXiv:2504.20543
             */
	    dLscaled_dx += 3*Lscaled2PNCoeff;

	    dLscaled_detsq *= xbar;
            dLscaled_detsq += Lscaled2PNCoeff_detsq;

            Lscaled_value *= xbar;
            Lscaled_value += Lscaled2PNCoeff;

#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
        __attribute__ ((fallthrough));
#endif
        case(3):
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
        __attribute__ ((fallthrough));
#endif
        case(2):
	        dLscaled_dx *= xbar;

            /* 1PN term for derivarive of scaled orbital angular momentum w.r.t. x
             * Equation C5b of Phukon et al., arXiv:2504.20543
             */
            dLscaled_dx += Lscaled1PNCoeff;

	    dLscaled_detsq *= xbar;
            dLscaled_detsq += Lscaled1PNCoeff_detsq;

	    Lscaled_value *= xbar;
	    Lscaled_value += Lscaled1PNCoeff;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
        __attribute__ ((fallthrough));
#endif
        case(1):
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
        __attribute__ ((fallthrough));
#endif
        case(0):
            dLscaled_dx *= xbar;

            /* Newtonian term non-spinning term for derivarive of scaled orbital angular momentum w.r.t. x
             * Equation C5b of Phukon et al., arXiv:2504.20543
             */
            dLscaled_dx += -Lscaled0PNCoeff;

	    dLscaled_detsq *= xbar;
            dLscaled_detsq += Lscaled0PNCoeff_detsq;
	    dLscaled_detsq *= xbar;

            Lscaled_value *= xbar;
	    Lscaled_value += Lscaled0PNCoeff;
            break;
        default:
            XLALPrintError("XLAL Error - %s: - %d is an invalid phase order\n", __func__, phaseO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
     }

    /* 1.5PN spin term in the expression of scaled angular momentum variable
     * Equation C3b of Phukon et al., arXiv:2504.20543 */
     REAL8 Lscaled1p5PNCoeff_spin = - XLALSimInspiralSpinEccBeta( 10./3, (5-etSq)/2.,
                                                                s1x, s1y, s1z,
                                                                s2x, s2y, s2z,
                                                                lhatx, lhaty, lhatz, ak);
     REAL8 Gamma1 = ak->gamma1;
     REAL8 Beta01 = XLALSimInspiralSpinEccBeta ( 0, 1., s1x, s1y, s1z,
                                                        s2x, s2y, s2z,
                                                        lhatx, lhaty, lhatz, ak );

     /*  2PN spin term in the scaled angular momentum  variable
      *  The gamma1 from C3c of Phukon et al., arXiv:2504.20543 */
     REAL8 Lscaled2PNCoeff_spin = Gamma1;

     /* variable for spin terms in derivative of angular momentum variable w.r.t. x*/
     REAL8 dLscaled_dx_spin = 0;

     /* variable for spin terms in derivative of angular momentum variable w.r.t. eccentricity squared*/
     REAL8 dLscaled_detsq_spin = 0;

     /* variable for spin terms in angular momentum variable*/
     REAL8 Lscaled_spin = 0;

     switch(spinO)
     {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            /*2PN spin term from Eqs. C5b of Phukon et al., arXiv:2504.20543
             */
            dLscaled_dx_spin = 3 * xbar * Lscaled2PNCoeff_spin;

            /*2PN (coefficient of \bar{x}^3) spin term from Equation C5c of Phukon et al., arXiv:2504.20543
             */
            dLscaled_detsq_spin = 3 * xbar * Lscaled2PNCoeff_spin;

            Lscaled_spin = xbar * Lscaled2PNCoeff_spin;

#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
        __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            /*1.5PN term from Equation C5b of Phukon et al., arXiv:2504.20543
             */
            dLscaled_dx_spin += 2 * Sqrtxbar * Lscaled1p5PNCoeff_spin;
	    dLscaled_dx_spin *= xbar;

            /*1.5PN term from Equation C5c of Phukon et al., arXiv:2504.20543
             * Coefficient of \bar{x}^{5/2}
             */
            dLscaled_detsq_spin += Sqrtxbar * (one_minus_etsq*Beta01  + 2 * Lscaled1p5PNCoeff_spin);
	    dLscaled_detsq_spin *= xbar*xbar;

            Lscaled_spin += Sqrtxbar * Lscaled1p5PNCoeff_spin;
            Lscaled_spin *= xbar;

#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
        __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error in function - %s: Invalid spin PN order %d\n",__func__, spinO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
     }

     dLscaled_dx += dLscaled_dx_spin;
     dLscaled_detsq += dLscaled_detsq_spin;
     Lscaled_value += Lscaled_spin;
     *Lscaled = Lscaled_value;

     REAL8 dx_dt = 2./3 * mtotal/v * domega_dt;
     *dLscaled_dt = dx_dt*dLscaled_dx + detSq_dt*dLscaled_detsq;
     return XLAL_SUCCESS;
}

/**
 * @brief Internal function for integrator to compute derivatives of dynamical variables of evolution
 *
 * This internal function is called by the ODE integrator to compute the time derivatives
 * of the 17 dynamical variables so that ODE integrator can take a step.
 *
 * The dynamical variables are:
 * - y[0]: eccentricity squared
 * - y[1]: orbital frequency
 * - y[2]: mean anomaly
 * - y[3]: secular orbital phase
 * - y[4-6]: components of spin 1
 * - y[7-9]: components of spin 2
 * - y[10-12]: components o orbital angular momentum unit vector
 * - y[13-15]: components of periastron direction unit vector
 * - y[16]: periastron advance parameter
 *
 * @param[in] t     Time (unused parameter )
 * @param[in] y     Array of 17 dynamical variables
 * @param[out] ydot Array of 17 time derivatives of dynamical variables
 * @param[in] mparams Pointer to XLALSimInspiralSpinEccPNEvolveOrbitParams for non-dynamical quantities (masses, etc.) and equations
 *
 * @return GSL_SUCCESS on success, or LALSIMINSPIRAL_STE_DERIVATIVE_OMEGANONPOS if omega ≤ 0
 *
 */
int XLALSimInspiralSpinEccentricTaylorT4DerivativesAvg( double UNUSED t,  const double y[],
                                                        double ydot[], void *mparams )
{
    XLALSimInspiralSpinEccPNEvolveOrbitParams* p = (XLALSimInspiralSpinEccPNEvolveOrbitParams*) mparams;

    REAL8 omega = y[1];
    REAL8 etSq = y[0];

    if (omega <= 0.0) /* orbital frequency must be positive! */
    {
        return LALSIMINSPIRAL_STE_DERIVATIVE_OMEGANONPOS;
    }

    REAL8 s1x = y[4];
    REAL8 s1y = y[5];
    REAL8 s1z = y[6];
    REAL8 s2x = y[7];
    REAL8 s2y = y[8];
    REAL8 s2z = y[9];
    REAL8 lx = y[10];
    REAL8 ly = y[11];
    REAL8 lz = y[12];
    REAL8 phatx = y[13];
    REAL8 phaty = y[14];
    REAL8 phatz = y[15];

    XLALSimInspiralSpinEccDynamicalExpnCoeffsSpinEccVariables( &p->ak, omega, etSq,
                                                        s1x, s1y, s1z,
                                                        s2x, s2y, s2z,
                                                        lx, ly, lz,
                                                        phatx, phaty, phatz);


    ydot[0] = p->funcetSq( omega, etSq, y[2], y[3], y[4], y[5], y[6], y[7],
                            y[8],  y[9], y[10], y[11], y[12], y[13], y[14], y[15],
                            y[16], &p->ak);
    ydot[1] = p->funcomega( omega, etSq, y[2], y[3], y[4], y[5], y[6], y[7],
                            y[8],  y[9], y[10], y[11], y[12], y[13], y[14], y[15],
                            y[16], &p->ak);
    ydot[2] = p->funcl( omega, etSq, y[2], y[3], y[4], y[5], y[6], y[7],
                            y[8],  y[9], y[10], y[11], y[12], y[13], y[14], y[15],
                            y[16], &p->ak);
    ydot[3] = p->funclamb( omega, etSq, y[2], y[3], y[4], y[5], y[6], y[7],
                            y[8], y[9], y[10], y[11], y[12], y[13], y[14], y[15],
                            y[16], &p->ak);
    ydot[4] = p->funcs1x( omega, etSq, y[2], y[3], y[4], y[5], y[6], y[7],
                            y[8],  y[9], y[10], y[11], y[12], y[13], y[14], y[15],
                            y[16], &p->ak);
    ydot[5] = p->funcs1y( omega, etSq, y[2], y[3], y[4], y[5], y[6], y[7],
                            y[8],  y[9], y[10], y[11], y[12], y[13], y[14], y[15],
                            y[16], &p->ak);
    ydot[6] = p->funcs1z( omega, etSq, y[2], y[3], y[4], y[5], y[6], y[7],
                            y[8],  y[9], y[10], y[11], y[12], y[13], y[14], y[15],
                            y[16], &p->ak);
    ydot[7] = p->funcs2x( omega, etSq, y[2], y[3], y[4], y[5], y[6], y[7],
                            y[8],  y[9], y[10], y[11], y[12], y[13], y[14], y[15],
                            y[16], &p->ak);
    ydot[8] = p->funcs2y( omega, etSq, y[2], y[3], y[4], y[5], y[6], y[7],
                            y[8],  y[9], y[10], y[11], y[12], y[13], y[14], y[15],
                            y[16], &p->ak);
    ydot[9] = p->funcs2z( omega, etSq, y[2], y[3], y[4], y[5], y[6], y[7],
                            y[8],  y[9], y[10], y[11], y[12], y[13], y[14], y[15],
                            y[16], &p->ak);
    ydot[10] = p->funclx( omega, etSq, y[2], y[3], y[4], y[5], y[6], y[7],
                            y[8],  y[9], y[10], y[11], y[12], y[13], y[14], y[15],
                            y[16], &p->ak);
    ydot[11] = p->funcly( omega, etSq, y[2], y[3], y[4], y[5], y[6], y[7],
                            y[8],  y[9], y[10], y[11], y[12], y[13], y[14], y[15],
                            y[16], &p->ak);
    ydot[12] = p->funclz( omega, etSq, y[2], y[3], y[4], y[5], y[6], y[7],
                            y[8],  y[9], y[10], y[11], y[12], y[13], y[14], y[15],
                            y[16], &p->ak);
    ydot[13] = p->funcpx( omega, etSq, y[2], y[3], y[4], y[5], y[6], y[7],
                            y[8],  y[9], y[10], y[11], y[12], y[13], y[14], y[15],
                            y[16], &p->ak);
    ydot[14] = p->funcpy( omega, etSq, y[2], y[3], y[4], y[5], y[6], y[7],
                            y[8],  y[9], y[10], y[11], y[12], y[13], y[14], y[15],
                            y[16], &p->ak);
    ydot[15] = p->funcpz( omega, etSq, y[2], y[3], y[4], y[5], y[6], y[7],
                            y[8],  y[9], y[10], y[11], y[12], y[13], y[14], y[15],
                            y[16], &p->ak);
    ydot[16] = p->funck( omega, etSq, y[2], y[3], y[4], y[5], y[6], y[7],
                            y[8],  y[9], y[10], y[11], y[12], y[13], y[14], y[15],
                            y[16], &p->ak);
    return GSL_SUCCESS;
}

/**
 * @brief Internal function  for stopping test function in binary evolution
 *
 * This function implements stopping conditions for the ODE evolution of precessing eccentric binaries.
 * It checks for various pathological conditions that indicate the evolution should terminate,
 * including numerical instabilities, unphysical parameter values, and completion of the inspiral.
 *
 * The function stops the evolution if any of the following conditions are met:
 * - Orbital frequency has NaN value
 * - Instantaneous GW frequency is negative
 * - \f$\bar{x}=v^2/(1 - e^2)\f$ is larger than 1
 * - Instantaneous GW frequency is not bounded between start-end frequency
 * - Eccentricity squared is negative or greater than 1
 * - Binary is unbound (energy is positive, or energy variable epsilon is negative)
 * - Magnitude of orbital angular momentum is negative
 * - Energy is increasing with time (time derivative of epsilon is negative)
 * - Angular momentum magnitude is increasing (time derivative of magnitude of angular momentum is positive)
 * - Second time derivative of orbital frequency is negative
 * - Periastron precession parameter is negative
 *
 * @param[in]  t       Time (unused parameter)
 * @param[in]  y       Array of dynamical variables
 * @param[out] ydot    Array of derivatives of dynamical variables
 * @param[in]  params  Pointer to XLALSimInspiralSpinEccPNEvolveOrbitParams structure
 *
 * @return GSL_SUCCESS if evolution should continue, Stopping Code if evolution should stop
 *
 * @see Table I of Phukon et al., arXiv:2504.20543 for stopping condition codes
 */
static int XLALSimInspiralSpinTaylorPNEccentricStoppingTest( double UNUSED t,  const double y[], double ydot[], void *params)
{
    double etp2 = y[0];
    double omega = y[1];

    XLALSimInspiralSpinEccPNEvolveOrbitParams* p = (XLALSimInspiralSpinEccPNEvolveOrbitParams*) params;

    REAL8 mtotal = p->ak.m;

    REAL8 v = cbrt( mtotal*omega);
    REAL8 x = v*v;
    REAL8 xbar = x/(1 - etp2);

    REAL8 s1x = y[4];
    REAL8 s1y = y[5];
    REAL8 s1z = y[6];
    REAL8 s2x = y[7];
    REAL8 s2y = y[8];
    REAL8 s2z = y[9];
    REAL8 lhatx = y[10];
    REAL8 lhaty = y[11];
    REAL8 lhatz = y[12];

    REAL8 fstart = p->ak.fstart;
    REAL8 fend = p->ak.fend;

    /* instantaneous GW frequency */
    REAL8 fGW = omega/(LAL_PI*p->ak.mt);

    p->ak.gamma1 = XLALSimInspiralSpinEccGamma1(s1x, s1y, s1z, s2x, s2y,  s2z,
                                        lhatx, lhaty, lhatz, &p->ak );
    p->ak.Fe = XLALSimInspiralSpinEccEnhancementFlux(etp2);
    p->ak.beta43 = XLALSimInspiralSpinEccBeta(4, 3, s1x, s1y, s1z, s2x, s2y, s2z,
                                        lhatx, lhaty, lhatz, &p->ak);
    /* etsq_dot */
    ydot[0] = p->funcetSq( omega, etp2, y[2], y[3], s1x, s1y, s1z, s2x,
                            s2y,  s2z, lhatx, lhaty, lhatz, y[13], y[14], y[15],
                            y[16], &p->ak);
    /* omega_dot */
    ydot[1] = p->funcomega( omega, etp2, y[2], y[3], s1x, s1y, s1z, s2x,
                            s2y,  s2z, lhatx, lhaty, lhatz, y[13], y[14], y[15],
                            y[16], &p->ak);

    REAL8 d_epsilon_dt = 0;
    REAL8 epsilon = 0;
    REAL8 d_Lscaled_dt = 0;
    REAL8 Lscaled = 0;

    evaluate_depsilon_dt ( &epsilon, &d_epsilon_dt, xbar, etp2, ydot[1], ydot[0], &p->ak );
    evaluate_dL_scaled_dt ( &Lscaled, &d_Lscaled_dt, xbar, etp2, s1x, s1y, s1z, s2x, s2y, s2z,
                            lhatx, lhaty, lhatz, ydot[1], ydot[0], &p->ak );

    REAL8 ddomega;
    ddomega = ydot[1] - p->ak.prev_domega;

    if ( p->ak.fend < p->ak.fstart && p->ak.fend != 0.
                && p->ak.prev_domega != 0.)
        ddomega *= -1;

    p->ak.prev_domega = ydot[1];

    if (isnan(omega)) /* omega is nan! */
        return LALSIMINSPIRAL_STE_TEST_OMEGANAN;
    else if (xbar>=1.0)
        return LALSIMINSPIRAL_STE_TEST_LARGEXBAR;
    else if ( fabs(fend) > LAL_REAL8_EPS  &&  ((fend > fstart && fGW > fend) || (fend < fstart && fGW < fend)) )
        return LALSIMINSPIRAL_STE_TEST_FREQBOUND;
    else if (etp2 < 0 || etp2 >=1)
        return LALSIMINSPIRAL_STE_TEST_ETSQBOUND;
    else if (epsilon < 0)
        return LALSIMINSPIRAL_STE_UNBOUND_SYSTEM;
    else if (Lscaled <0)
        return LALSIMINSPIRAL_STE_NEGATIVE_MAG_ANGMOMENTUM;
    else if (d_epsilon_dt < 0)
        return LALSIMINSPIRAL_STE_TEST_ENERGY;
    else if (d_Lscaled_dt > 0 )
        return LALSIMINSPIRAL_STE_TEST_ANGMOMENTUM;
    else if ( ddomega <=0 ) /* d^2omega/dt^2 <= 0! */
        return  LALSIMINSPIRAL_STE_TEST_OMEGADOUBLEDOT;
    else if ( y[16]<0 )
        return LALSIMINSPIRAL_STE_NEGATIVE_PERIASTRON_ADVANCE;
    else
      return GSL_SUCCESS;
}

/**
 * @brief Fourth-order Runge-Kutta ODE integrator with adaptive step size and interpolation for final state output
 *
 * This function performs adaptive Runge-Kutta integration using Runge-Kutta-Fehlberg (RKF45) steps
 * with adaptive step size control. It outputs only the final state of the integration and checks the
 * integrator's full stopping condition, specialized to the specifics of the stopping test considered here,
 * which needs the time derivative of omega (angular frequency) at the previous timestep. This function is
 * based on XLALAdaptiveRungeKutta4HermiteOnlyFinal, which just checks a single stopping condition.
 *
 * The method is described in
 *
 * Abramowitz & Stegun, Handbook of Mathematical Functions, Tenth Printing,
 * National Bureau of Standards, Washington, DC, 1972
 * (available online at http://people.math.sfu.ca/~cbm/aands/)
 *
 * This function also includes "on-the-fly" interpolation of the
 * differential equations at regular intervals in-between integration
 * steps. This "on-the-fly" interpolation method is derived and
 * described in the Mathematica notebook "RKF_with_interpolation.nb";
 * @see https://ldas-jobs.ligo.caltech.edu/~cbc.moin/pages/InspiralPipelineDevelopment(2f)120312111836InspiralPipelineDevelopmentImproved(20)Adaptive(20)Runge(2d)Kutta(20)integrator/attachments/
 *
 * @param[in,out] integrator  LALAdaptiveRungeKuttaIntegrator structure holding dydt/ODEs, stopping test, stepper, etc. -  status code returned
 * @param[in]     params      Structure to extract variables needed from expnCoeffsSpinEcc structure in intregration
 * @param[in,out] yinit       Pass in initial values of variables of ODEs - overwritten to final values
 * @param[in]     tinit       Integration start time [in (\f$GM/c^3\f$)]
 * @param[in]     tend_in     Maximum integration time [in (\f$GM/c^3\f$)]
 * @param[in]     deltat      Time step size for integration [in (\f$GM/c^3\f$)]
 */
int XLALAdaptiveRungeKutta4HermiteOnlyFinalCheckStopping(LALAdaptiveRungeKuttaIntegrator * integrator,
                                                        void *params, REAL8 * yinit, REAL8 tinit,
                                                        REAL8 tend_in, REAL8 deltat )
{
    int errnum = 0;
    int stopping_interp = 0;
    int status;
    int interp_trigger=0;
    size_t dim, retries=0, i;

    REAL8 t, tintp, h;

    REAL8 *ytemp = NULL;
    REAL8 *ytemp2 = NULL;

    REAL8 tend = tend_in;

    XLALSimInspiralSpinEccPNEvolveOrbitParams* p = (XLALSimInspiralSpinEccPNEvolveOrbitParams*) params;
    REAL8 fend = p->ak.fend;
    REAL8 y1_final = fend*(LAL_PI*p->ak.mt);
    XLAL_BEGINGSL;

    /*
     * If want to stop only on test, then tend = +/-infinity; otherwise
     * tend_in
     *
     */
    if (integrator->stopontestonly) {
        if (tend < tinit)
            tend = -1.0 / 0.0;
        else
            tend = 1.0 / 0.0;
    }

    dim = integrator->sys->dimension;

    if ((tend < tinit && deltat > 0) || (tend > tinit && deltat < 0)) {
        XLALPrintError
            ("XLAL Error - %s: (tend_in - tinit) and deltat must have the same sign\ntend_in: %f, tinit: %f, deltat: %f\n",
            __func__, tend_in, tinit, deltat);
        errnum = XLAL_EINVAL;
        goto bail_out;
    }

    ytemp = XLALCalloc(dim, sizeof(REAL8));
    ytemp2 = XLALCalloc(dim, sizeof(REAL8));

    if (!ytemp || !ytemp2) {
        errnum = XLAL_ENOMEM;
        goto bail_out;
    }

    /*
     * Initialize ytemp[1] with the initial value of yinit[1] so that the initial check below is satisfied even if we are integrating backwards
     *
     */
    ytemp[1] = yinit[1];

    /* Setup. */
    integrator->sys->params = params;
    integrator->returncode = 0;
    retries = integrator->retries;
    t = tinit;
    tintp = tinit;
    h = deltat;

    /*
     * We are starting a fresh integration; clear GSL step and evolve
     * objects.
     *
     */
    gsl_odeiv_step_reset(integrator->step);
    gsl_odeiv_evolve_reset(integrator->evolve);

    /*
     * Enter evolution loop.  NOTE: we *always* take at least one
     * step.
     *
     */
    while (1) {
        REAL8 told = t;

        REAL8 prev_domega = p->ak.prev_domega;

        status =
            gsl_odeiv_evolve_apply(integrator->evolve, integrator->control, integrator->step, integrator->sys, &t, tend, &h,
            yinit);

        /*
         * Check for failure, retry if haven't retried too many times
         * already.
         *
         */
        if (status != GSL_SUCCESS) {
            if (retries--) {
                /* Retries to spare; reduce h, try again. */
                h /= 10.0;
                continue;
            } else {
                /* Out of retries, bail with status code. */
                integrator->returncode = status;
                break;
            }
        } else {
            /* Successful step, reset retry counter. */
            retries = integrator->retries;
        }

        /*
         * If we're at the final timestep, interpolate to get the output at the desired final value, y1_final.
         * Note we square to get an absolute value, because we may be
         * integrating t in the positive or negative direction
         * We multiply yinit, ytemp[1], and y1_final by deltat so that this check works when integrating forward or backward.
         * We also include the integrator's stopping test in the check (assuming that if it is triggered at an intermediate timestep it
         * will also be triggered at the final timestep) and output the last interpolated value before the stopping test is triggered
         *
         */
        if ((interp_trigger = integrator->stop(t, yinit, integrator->evolve->dydt_out, params)) != GSL_SUCCESS) {
                tintp = told;

                REAL8 hUsed = t - told;

                /* Set prev_domega to the saved prev_domega to have the appropriate value for the first of the subsequent calls to the stopping test */
                p->ak.prev_domega = prev_domega;

                /* Save initial value to allow for the case where the stopping condition is triggered at i = 1 */

                /* First grab the k's from the integrator state. */
                rkf45_state_t *rkfState = integrator->step->state;
                REAL8 *k1 = rkfState->k1;
                REAL8 *k6 = rkfState->k6;
                REAL8 *y0 = rkfState->y0;

                for (i = 0; i < dim; i++) {
                  ytemp[i] = y0[i];
                }
                XLALPrintInfo("XLAL Info - %s: interpolation triggered by code %d.\n", __func__, interp_trigger);

                while (1) {
                  tintp += deltat;

                  /*
                   * tintp = told + (t-told)*theta, 0 <= theta <= 1.  We have to
                   * compute h = (t-told) because the integrator returns a
                   * suggested next h, not the actual stepsize taken.
                   *
                   */
                  REAL8 theta = (tintp - told) / hUsed;

                  /*
                   * These are the interpolating coefficients for y(t + h*theta) =
                   * ynew + i1*h*k1 + i5*h*k5 + i6*h*k6 + O(h^4).
                   *
                   */
                  REAL8 i0 = 1.0 + theta * theta * (3.0 - 4.0 * theta);
                  REAL8 i1 = -theta * (theta - 1.0);
                  REAL8 i6 = -4.0 * theta * theta * (theta - 1.0);
                  REAL8 iend = theta * theta * (4.0 * theta - 3.0);

                  for (i = 0; i < dim; i++) {
                    ytemp2[i] = i0 * y0[i] + iend * yinit[i] + hUsed * i1 * k1[i] + hUsed * i6 * k6[i];
                  }


                    if (interp_trigger==1029) {
                       stopping_interp = (deltat * ytemp2[1] < deltat * y1_final);
                       status = 1029;
                   }
                    else {
                        status = integrator->stop(t, ytemp2, integrator->evolve->dydt_out, params);
                        stopping_interp = (status == GSL_SUCCESS);
                    }

                  if (stopping_interp) {
                    for (i = 0; i < dim; i++) {
                      ytemp[i] = ytemp2[i];
                    }
                  } else {
                    break;
                  }
                }
        }

        /* Now check for termination criteria. */
        if (!integrator->stopontestonly && t >= tend)
            break;

        if (status != GSL_SUCCESS) {
            integrator->returncode = status;
                break;
        }
    }

    /* Store the final *interpolated* sample before the stopping criterion in yinit. */
    for (i = 0; i < dim; i++) {
      yinit[i] = ytemp[i];
    }


  bail_out:

    XLAL_ENDGSL;

    /*
     * If we have an error, then we should free allocated memory, and
     * then return.
     *
     */
    XLALFree(ytemp);
    XLALFree(ytemp2);

    if (errnum) {
        XLAL_ERROR(errnum);
    }

    return 1;
}

/**
 * @brief Compute squared norms of cross products and dot products between unit orbital angular momentum and
 * dimensionless spin vectors
 *
 * This function calculates:
 * 1. Squared norms of cross products between the unit vector of orbital angular momentum and the dimensionless spin vectors of the binary components
 * 2. Dot products between the dimensionless spin vectors and the unit vector of orbital angular momentum
 *
 * @param[out] lhatCrosss1normSq  Squared norm of orbital angular momentum unit vector and the dimensionless spin vectors 1
 * @param[out] lhatCrosss2normSq  Squared norm of orbital angular momentum unit vector and the dimensionless spin vectors 2
 * @param[out] lhatdots1          Dot product of the dimensionless spin vector 1 and orbital angular momentum unit vector
 * @param[out] lhatdots2          Dot product of the dimensionless spin vector 2 and orbital angular momentum
 * @param[in]  s1x,s1y,s1z        Cartesian components of dimensionless spin vector for body 1
 * @param[in]  s2x,s2y,s2z        Cartesian components of dimensionless spin vector for body 2
 * @param[in]  lhatx,lhaty,lhatz  Unit vector components of orbital angular momentum
 *
 * @return XLAL_SUCCESS on successful completion
 */
static REAL8 XLALSimInspiralComputeCrossDotQuants( REAL8 *lhatCrosss1normSq, REAL8 *lhatCrosss2normSq,
	                                               REAL8 *lhatdots1, REAL8 *lhatdots2,
	                                               REAL8 s1x, REAL8 s1y, REAL8 s1z,
	                                               REAL8 s2x, REAL8 s2y, REAL8 s2z,
                                                   REAL8 lhatx, REAL8 lhaty, REAL8 lhatz)

{
   REAL8 *lhatCrosss1_internal=NULL;
   REAL8 *lhatCrosss2_internal=NULL;

   XLALSimInspiralVectorCrossProduct(  &lhatCrosss1_internal, lhatx, lhaty, lhatz, s1x, s1y, s1z);
   XLALSimInspiralVectorCrossProduct(  &lhatCrosss2_internal, lhatx, lhaty, lhatz, s2x, s2y, s2z);

   *lhatCrosss1normSq = (lhatCrosss1_internal[0]*lhatCrosss1_internal[0] + lhatCrosss1_internal[1]*lhatCrosss1_internal[1] + lhatCrosss1_internal[2]*lhatCrosss1_internal[2]);
   *lhatCrosss2normSq = (lhatCrosss2_internal[0]*lhatCrosss2_internal[0] + lhatCrosss2_internal[1]*lhatCrosss2_internal[1] + lhatCrosss2_internal[2]*lhatCrosss2_internal[2]);

   *lhatdots1 = ( s1x*lhatx  + s1y*lhaty + s1z*lhatz );
   *lhatdots2 = ( s2x*lhatx  + s2y*lhaty + s2z*lhatz );

    if (lhatCrosss1_internal !=NULL){
        LALFree(lhatCrosss1_internal);
    }
    if (lhatCrosss2_internal != NULL){
        LALFree(lhatCrosss2_internal);
    }

   return XLAL_SUCCESS;
}

/**
 * @brief Compute the Pade contribution to the superasymptotic 2.5PN alpha enhancement function
 *
 * This function calculates the Pade approximant contribution to the superasymptotic expression
 * for the 2.5PN alpha enhancement function `@ref XLALSimInspiralSpinEccAlphaSuperAsyLY`
 *
 * @param[in] etSq  Square of the orbital eccentricity
 *
 * @return The Pade contribution to the 2.5PN alpha enhancement function
 *
 * @see Eqs. 168, B3, B4 of Loutrel & Yunes, arXiv:1607.05409
 */
REAL8 XLALSimInspiralSpinEccPadePolyA( REAL8 etSq )
{
    REAL8 A_sup_2 = 2517;
    REAL8 A_sup_4 = -542.748497072893;
    REAL8 A_sup_6 = -8106.76855568073;
    REAL8 A_sup_8  = 8413.80043346759;
    REAL8 A_sup_10 = -1797.84373379634;
    REAL8 A_sup_12 = -645.518003878721;
    REAL8 A_sup_14 = 196.536287243258;
    REAL8 A_sub_2 = -2.73321354671152;
    REAL8 A_sub_4 = 2.72445722875519;
    REAL8 A_sub_6 = -1.18148275600057;
    REAL8 A_sub_8 = 0.198189499682584;
    REAL8 A_sub_10 = -0.0052713663843431;
    REAL8 A_sub_12 = -0.000104564833689961;
    REAL8 A_sub_14 = -0.00001735320614776194;
    REAL8 A_sub_16 = -2.68866627647753e-6;

    REAL8 numerator = etSq*(A_sup_2 + etSq*( A_sup_4 + etSq * ( A_sup_6 + etSq * ( A_sup_8 +
                        etSq * (A_sup_10 + etSq * (A_sup_12 + etSq *A_sup_14))))));

    REAL8 denominator = 1 + etSq*(A_sub_2 + etSq*( A_sub_4 + etSq * ( A_sub_6 + etSq * ( A_sub_8 +
                        etSq * (A_sub_10 + etSq * (A_sub_12 + etSq *(A_sub_14 + etSq*A_sub_16)))))));

    return numerator/denominator;
}

/**
 * @brief Compute the Pade contribution to the superasymptotic 2.5PN alpha tilde enhancement function
 *
 * This function calculates the Pade approximant contribution to the superasymptotic expression
 * for the 2.5PN alpha tilde enhancement function `\@ref XLALSimInspiralSpinEccTildeAlphaSuperAsyLY`
 * The tilde version uses expression similar to Eqs. 168 of Loutrel & Yunes, arXiv:1607.05409 but
 * with modified coefficients.
 *
 * @param[in] etSq  Square of the orbital eccentricity
 *
 * @return The Pade contribution to the 2.5PN alpha tilde enhancement function
 *
 * @see Eqs. 168, B5, B6 of Loutrel & Yunes, arXiv:1607.05409
 */
REAL8 XLALSimInspiralSpinEccPadePolyTildeA( REAL8 etSq )
{
    REAL8 A_sup_2 = 1428;
    REAL8 A_sup_4 = -3082.0825943902;
    REAL8 A_sup_6 = 1251.7504530861;
    REAL8 A_sup_8  = 1534.29667420639;
    REAL8 A_sup_10 =  -1557.48894312597;
    REAL8 A_sup_12 = 474.00656591547;
    REAL8 A_sup_14 = -44.7331690120333;
    REAL8 A_sub_2 = -3.15201862352255;
    REAL8 A_sub_4 = 3.82292854565682;
    REAL8 A_sub_6 = -2.22448315038514;
    REAL8 A_sub_8 = 0.6262354256498;
    REAL8 A_sub_10 = -0.0738364336330394;
    REAL8 A_sub_12 = 0.0021699571134637;
    REAL8 A_sub_14 = 0.0000142157331173306;
    REAL8 A_sub_16 = -3.8287501968869e-6;

    REAL8 numerator = etSq*(A_sup_2 + etSq*( A_sup_4 + etSq * ( A_sup_6 + etSq * ( A_sup_8 +
                        etSq * (A_sup_10 + etSq * (A_sup_12 + etSq *A_sup_14))))));

    REAL8 denominator = 1 + etSq*(A_sub_2 + etSq*( A_sub_4 + etSq * ( A_sub_6 + etSq * ( A_sub_8 +
                        etSq * (A_sub_10 + etSq * (A_sub_12 + etSq *(A_sub_14 + etSq*A_sub_16)))))));

    return numerator/denominator;
}
