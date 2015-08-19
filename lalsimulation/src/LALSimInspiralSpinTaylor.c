/*
 * Copyright (C) 2011 E. Ochsner, 2014 A. Klein, 2015 R. Sturani
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
#include <lal/FrequencySeries.h>
#include "check_series_macros.h"
#include "LALSimInspiralPNCoefficients.c"
#include <lal/XLALGSL.h>

#define XLAL_BEGINGSL \
        { \
          gsl_error_handler_t *saveGSLErrorHandler_; \
          saveGSLErrorHandler_ = gsl_set_error_handler_off();

#define XLAL_ENDGSL \
          gsl_set_error_handler( saveGSLErrorHandler_ ); \
        }

//#define UNUSED(expr) do { (void)(expr); } while (0)
/*
#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif
*/

/* use error codes above 1024 to avoid conflicts with GSL */
#define LALSIMINSPIRAL_ST_TEST_ENERGY               1025
#define LALSIMINSPIRAL_ST_TEST_OMEGADOUBLEDOT       1026
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

/* Macro functions to rotate the components of a vector about an axis */
#define ROTATEZ(angle, vx, vy, vz)\
	tmp1 = vx*cos(angle) - vy*sin(angle);\
	tmp2 = vx*sin(angle) + vy*cos(angle);\
	vx = tmp1;\
	vy = tmp2

#define ROTATEY(angle, vx, vy, vz)\
	tmp1 = vx*cos(angle) + vz*sin(angle);\
	tmp2 = - vx*sin(angle) + vz*cos(angle);\
	vx = tmp1;\
	vz = tmp2

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
	REAL8 wdot3S1O, wdot3S2O; ///< non-dynamical 1.5PN SO corrections
        REAL8 wdot4S1S2, wdot4S1OS2O; ///< non-dynamical 2PN SS corrections
        REAL8 wdot4S1S1,wdot4S2S2; ///< non-dynamical self S^2 2PN correction
        REAL8 wdot4S1OS1O, wdot4S2OS2O; ///< non-dynamical self SO^2 2PN correction
        REAL8 wdot4QMS1S1; ///< non-dynamical S1^2 2PN quadrupole-monopole correct
        REAL8 wdot4QMS1OS1O; ///< non-dynamical (S1.L)^2 2PN quadrupole-monopole co
        REAL8 wdot4QMS2S2; ///< non-dynamical S2^2 2PN quadrupole-monopole correct
        REAL8 wdot4QMS2OS2O; ///< non-dynamical (S2.L)^2 2PN quadrupole-monopole c
	REAL8 wdot5S1O, wdot5S2O; ///< non-dynamical 2.5PN SO corrections
	REAL8 wdot6S1O, wdot6S2O; ///< non-dynamical 3PN SO corrections
	REAL8 wdot7S1O, wdot7S2O; ///< non-dynamical 3.5PN SO corrections
	REAL8 E3S1O, E3S2O; ///< non-dynamical 1.5PN SO corrections
        REAL8 E4S1S2,E4S1OS2O; ///< non-dynamical 2PN SS correction
	REAL8 E4QMS1S1; ///< non-dynamical S1^2 2PN quadrupole-monopole correction
	REAL8 E4QMS1OS1O;///< non-dynamical (S1.L)^2 2PN quadrupole-monopole correction
	REAL8 E4QMS2S2; ///< non-dynamical S2^2 2PN quadrupole-monopole correction
	REAL8 E4QMS2OS2O;///< non-dynamical (S2.L)^2 2PN quadrupole-monopole correction
	REAL8 E5S1O, E5S2O; ///< non-dynamical 2.5PN SO corrections
	REAL8 E7S1O, E7S2O; ///< non-dynamical 3.5PN SO corrections
	REAL8 LNhatdot3S1O, LNhatdot3S2O; ///< non-dynamical 1.5PN SO corrections
	REAL8 LNhatdot4S1S2; ///< non-dynamical 2PN SS correction
	REAL8 LNhatdot4QMS1; ///< non-dynamical 2PN SS correction
	REAL8 LNhatdot4QMS2; ///< non-dynamical 2PN SS correction
	REAL8 LNhatdot5S1O, LNhatdot5S2O; ///< non-dynamical 2.5PN SO corrections
	REAL8 wdottidal5pn;	///< leading order tidal correction 
	REAL8 wdottidal6pn;	///< next to leading order tidal correction
	REAL8 Etidal5pn; ///< leading order tidal correction to energy
	REAL8 Etidal6pn; ///< next to leading order tidal correction to energy
	REAL8 fStart; ///< starting GW frequency of integration
	REAL8 fEnd; ///< ending GW frequency of integration
	LALSimInspiralSpinOrder spinO; ///< Twice PN order of included spin effects
	LALSimInspiralTidalOrder tideO;///< Twice PN order of included tidal effects
	REAL8 prev_domega; ///< Previous value of domega/dt used in stopping test
	REAL8 Fcoeff[LAL_MAX_PN_ORDER];///FluxCoeff
        REAL8 Fnewt; // newtonian term in Flux
	REAL8 Flogcoeff; //log coeff in flux
	REAL8 F7S1O;
	REAL8 F7S2O;
	REAL8 F6S1O;
	REAL8 F6S2O;
	REAL8 F5S1O;
	REAL8 F5S2O;
	REAL8 F4S1S1;
	REAL8 F4S1OS1O;
	REAL8 F4S2S2;
	REAL8 F4S2OS2O;
        REAL8 F4QMS1S1;
        REAL8 F4QMS2S2;
        REAL8 F4QMS1OS1O;
        REAL8 F4QMS2OS2O;
	REAL8 F3S1O;
	REAL8 F3S2O;
	REAL8 F4S1S2;
	REAL8 F4S1OS2O;
	REAL8 Ftidal5pn;     ///< leading order tidal correction
        REAL8 Ftidal6pn;
        REAL8 dEdvnewt;
        REAL8 S1dot3;
        REAL8 S2dot3;
        REAL8 S1dot5;
        REAL8 S2dot5;
        REAL8 Sdot4S1S2;
        REAL8 Sdot4S1OS2O;
        REAL8 S1dot4QMS1O;
        REAL8 S2dot4QMS2O;
} XLALSimInspiralSpinTaylorTxCoeffs;


/* Declarations of static functions - defined below */
static int XLALSimInspiralSpinTaylorStoppingTest(double t,
	const double values[], double dvalues[], void *mparams);
static int XLALSimInspiralSpinTaylorT4Derivatives(double t, 
	const double values[], double dvalues[], void *mparams);
static int XLALSimInspiralSpinTaylorT1Derivatives(double t,
	const double values[], double dvalues[], void *mparams);
static int XLALSimInspiralSpinTaylorT4Setup(
    XLALSimInspiralSpinTaylorTxCoeffs *params, REAL8 m1, REAL8 m2,
    REAL8 fStart, REAL8 fEnd, REAL8 lambda1, REAL8 lambda2,
    REAL8 quadparam1, REAL8 quadparam2, LALSimInspiralSpinOrder spinO,
    LALSimInspiralTidalOrder tideO, INT4 phaseO);
static int XLALSimInspiralSpinTaylorT1Setup(
    XLALSimInspiralSpinTaylorTxCoeffs *params, REAL8 m1, REAL8 m2,
    REAL8 fStart, REAL8 fEnd, REAL8 lambda1, REAL8 lambda2,REAL8 quadparam1, REAL8 quadparam2,
    LALSimInspiralSpinOrder spinO, LALSimInspiralTidalOrder tideO, INT4 phaseO);
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
static int XLALSimInspiralSpinTaylorPNEvolveOrbitIrregularIntervals(
    REAL8Array **yout, REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fEnd, REAL8 s1x,
    REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lnhatx,
    REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z, REAL8 lambda1,
    REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2,
    LALSimInspiralSpinOrder spinO, LALSimInspiralTidalOrder tideO, INT4 phaseO,
    Approximant approx);
static int XLALSimInspiralSpinTaylorDriverFourier(
    COMPLEX16FrequencySeries **hplus, COMPLEX16FrequencySeries **hcross,
    REAL8 fMin, REAL8 fMax, REAL8 deltaF, INT4 kMax, REAL8 phiRef, REAL8 v0,
    REAL8 m1, REAL8 m2, REAL8 fStart, REAL8 fRef, REAL8 r, REAL8 s1x,
    REAL8 s1y, REAL8 s1z, REAL8 s2x, REAL8 s2y, REAL8 s2z, REAL8 lnhatx,
    REAL8 lnhaty, REAL8 lnhatz, REAL8 e1x, REAL8 e1y, REAL8 e1z, REAL8 lambda1,
    REAL8 lambda2, REAL8 quadparam1, REAL8 quadparam2,
    LALSimInspiralSpinOrder spinO, LALSimInspiralTidalOrder tideO, INT4 phaseO,
    INT4 amplitudeO, Approximant approx, INT4 phiRefAtEnd);



/* Appends the start and end time series together. Frees start and end before
 * returning a pointer to the result. The last point of start can be the
 * first point of end, and if so it removes the duplicate point.
 */
static REAL8Array *appendTAandFree(REAL8Array *start,
        REAL8Array *end, REAL8Array* combined) {
  UINT4 len1 = start->dimLength->data[1];
  UINT4 len2 = end->dimLength->data[1];
  UINT4 lenTot;

  UINT4 nParams = start->dimLength->data[0];

  if(end->dimLength->data[0] != nParams)
  {
    XLALPrintError("XLAL Error - %s: cannot append series with different numbers of parameters %d and %d.\n", __func__, nParams, end->dimLength->data[0]);
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  UINT4 i;
  UINT4 doRemove;

  if(start->data[len1-1] != end->data[0])
  {
    doRemove = 0;
    lenTot = len1 + len2;
  }
  else
  {
    doRemove = 1;
    lenTot = len1 + len2 - 1;
  }
  combined = XLALCreateREAL8ArrayL(2, nParams, lenTot);

  for(i = 0; i < nParams; i++)
  {
    memcpy(&(combined->data[i*lenTot]), &(start->data[i*len1]), sizeof(REAL8)*len1);
    if(doRemove)
    {
      memcpy(&(combined->data[i*lenTot + len1]), &(end->data[i*len2+1]), sizeof(REAL8)*(len2-1));
      if(i && start->data[i*len1 + len1-1] != end->data[i*len2])
      {
        XLALPrintWarning("XLAL Warning - %s: time series inconsistent: parameter %d is %f or %f at the same time.\n", __func__, i, start->data[i*len1 + len1-1], end->data[i*len2]);
      }
    }
    else
    {
      memcpy(&(combined->data[i*lenTot + len1]), &(end->data[i*len2]), sizeof(REAL8)*(len2));
    }
  }

  XLALDestroyREAL8Array(start);
  XLALDestroyREAL8Array(end);

  return combined;
}


/* Remove duplicates if the time is constant.
   The time has to be the first parameter.
 */
static REAL8Array *removeDuplicates(REAL8Array *series) {
  UINT4 lenTot = series->dimLength->data[1];
  UINT4 newLen = lenTot;

  UINT4 nParams = series->dimLength->data[0];
  UINT4 i, j, k;
  UINT4 nDuplicates;

  if(nParams < 1)
  {
    XLALPrintError("XLAL Error - %s: bad time series, has %d parameters.\n", __func__, nParams);
    XLAL_ERROR_NULL(XLAL_EINVAL);
  }

  REAL8Array *temp;
  temp = XLALCreateREAL8ArrayL(2, nParams, lenTot);

  REAL8Array *mean;
  mean = XLALCreateREAL8ArrayL(1, nParams);

  for(i = 0, k = 0; i < lenTot; i++, k++)
  {
    nDuplicates = 0;
    if(i < lenTot && series->data[i] == series->data[i+1])
    {
      for(j = 1; j < nParams; j++)
      {
        mean->data[j] = series->data[j*lenTot + i];
      }
      while(i < lenTot && series->data[i] == series->data[i+1])
      {
        nDuplicates++;
        i++;
        mean->data[0] = series->data[i];
        for(j = 1; j < nParams; j++)
        {
          mean->data[j] += series->data[j*lenTot + i];
        }
      }
      for(j = 1; j < nParams; j++)
      {
        mean->data[j] /= nDuplicates+1;
      }
      newLen -= nDuplicates;
      for(j = 0; j < nParams; j++)
      {
        temp->data[j*lenTot + k] = mean->data[j];
      }
    }
    else
    {
      for(j = 0; j < nParams; j++)
      {
        temp->data[j*lenTot + k] = series->data[j*lenTot + i];
      }
    }
  }

  if(newLen != lenTot)
  {
    XLALDestroyREAL8Array(series);
    series = XLALCreateREAL8ArrayL(2, nParams, newLen);
    for(j = 0; j < nParams; j++)
    {
      memcpy(&series->data[j*newLen], &temp->data[j*lenTot], sizeof(REAL8)*newLen);
    }
  }

  XLALDestroyREAL8Array(temp);
  XLALDestroyREAL8Array(mean);
  return series;
}

static int XLALSimInspiralSpinTaylorSpinDerivativeSetup(
    XLALSimInspiralSpinTaylorTxCoeffs *params, /**< UNDOCUMENTED */
    REAL8 m1,                       /**< mass of body 1 (kg) */
    REAL8 m2,                       /**< mass of body 2 (kg) */
    REAL8 quadparam1,               /**< phenom. parameter describing induced quad. moment of body 1 (=1 for BHs, ~2-12 for NSs) */
    REAL8 quadparam2,               /**< phenom. parameter describing induced quad. moment of body 2 (=1 for BHs, ~2-12 for NSs) */
    LALSimInspiralSpinOrder spinO   /**< twice PN order of spin effects */)
{
  REAL8 M=m1+m2;
  REAL8 m1M=m1/M;
  REAL8 m2M=m2/M;
  REAL8 eta=m1*m2/M/M;
  /* Compute the non-dynamical coefficients of spin corrections
   * to the evolution equations for L, S1 and S2
   */
  switch( spinO )
    {
    case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
    case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
    case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
    case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
      params->S1dot5        = XLALSimInspiralSpinDot_5PNCoeff(m1M);
      params->S2dot5        = XLALSimInspiralSpinDot_5PNCoeff(m2M);
      params->LNhatdot5S1O  = XLALSimInspiralLNhdot_5PNSOCoeff(m1M);
      params->LNhatdot5S2O  = XLALSimInspiralLNhdot_5PNSOCoeff(m2M);
    case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
      // 2PN spin-spin terms
      params->LNhatdot4S1S2 = XLALSimInspiralLNhdot_4PNS1S2Coeff(eta);
      params->Sdot4S1S2     = XLALSimInspiralSpinDot_4PNS1S2Coeff;
      params->Sdot4S1OS2O   = XLALSimInspiralSpinDot_4PNS1OS2OCoeff;
      params->LNhatdot4QMS1  = quadparam1 * XLALSimInspiralLNhdot_4PNQMSSCoeff(m1M);
      params->LNhatdot4QMS2  = quadparam2 * XLALSimInspiralLNhdot_4PNQMSSCoeff(m2M);
      params->S1dot4QMS1O    = quadparam1 * XLALSimInspiralSpinDot_4PNSOSOselfCoeff(m1M);
      params->S2dot4QMS2O    = quadparam2 * XLALSimInspiralSpinDot_4PNSOSOselfCoeff(m2M);
    case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
      params->S1dot3    = XLALSimInspiralSpinDot_3PNCoeff(m1M);
      params->S2dot3    = XLALSimInspiralSpinDot_3PNCoeff(m2M);
      params->LNhatdot3S1O = XLALSimInspiralLNhdot_3PNSOCoeff(m1M);
      params->LNhatdot3S2O = XLALSimInspiralLNhdot_3PNSOCoeff(m2M);
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
  return XLAL_SUCCESS;
}

/**
 * Function called by all functions computing derivatives, as
 * angular momentum/spin derivative are common for SpinTaylors.
 */
static INT4 XLALSimInspiralSpinDerivatives(REAL8 *dLNhx,
					   REAL8 *dLNhy,
					   REAL8 *dLNhz,
					   REAL8 *dE1x,
					   REAL8 *dE1y,
					   REAL8 *dE1z,
					   REAL8 *dS1x,
					   REAL8 *dS1y,
					   REAL8 *dS1z,
					   REAL8 *dS2x,
					   REAL8 *dS2y,
					   REAL8 *dS2z,
					   const REAL8 v,
					   const REAL8 omega,
					   const REAL8 LNhx,
					   const REAL8 LNhy,
					   const REAL8 LNhz,
					   const REAL8 E1x,
					   const REAL8 E1y,
					   const REAL8 E1z,
					   const REAL8 S1x,
					   const REAL8 S1y,
					   const REAL8 S1z,
					   const REAL8 S2x,
					   const REAL8 S2y,
					   const REAL8 S2z,
					   const REAL8 LNhdotS1,
					   const REAL8 LNhdotS2,
					   XLALSimInspiralSpinTaylorTxCoeffs *params)
{

  REAL8 v5=omega*v*v;
  REAL8 omega2=omega*omega;
  REAL8 v7=omega2*v;
  REAL8 v8=omega2*v*v;

  /*
   * dLN
   *
   * \f$d \hat{L_N}/d \hat{t} = M * d\hat{L_N} / dt = \Omega_L x \hat{L_N}\f$
   * This is Eq. (10) of gr-qc/0405090 ( times M b/c we use \f$\hat{t}\f$)
   */

  /* Compute \Omega_L vector */
  REAL8 OmegaLx = omega2 * (params->LNhatdot3S1O * S1x + params->LNhatdot3S2O * S2x)
    + v7 * ( params->LNhatdot4S1S2 * (LNhdotS2 * S1x + LNhdotS1 * S2x) + params->LNhatdot4QMS1*LNhdotS1*S1x + params->LNhatdot4QMS2*LNhdotS2*S2x )
    + v8 * ( params->LNhatdot5S1O * S1x + params->LNhatdot5S2O * S2x );
  REAL8 OmegaLy = omega2 * (params->LNhatdot3S1O * S1y + params->LNhatdot3S2O * S2y)
    + v7 * ( params->LNhatdot4S1S2 * (LNhdotS2 * S1y + LNhdotS1 * S2y) + params->LNhatdot4QMS1*LNhdotS1*S1y + params->LNhatdot4QMS2*LNhdotS2*S2y )
    + v8 * ( params->LNhatdot5S1O * S1y + params->LNhatdot5S2O * S2y );
  REAL8 OmegaLz = omega2 * (params->LNhatdot3S1O * S1z + params->LNhatdot3S2O * S2z)
    + v7 * ( params->LNhatdot4S1S2 * (LNhdotS2 * S1z + LNhdotS1 * S2z) + params->LNhatdot4QMS1*LNhdotS1*S1z + params->LNhatdot4QMS2*LNhdotS2*S2z )
    + v8 * ( params->LNhatdot5S1O * S1z + params->LNhatdot5S2O * S2z );

  /* Take cross product of \Omega_L with \hat{L_N} */
  *dLNhx = (-OmegaLz*LNhy + OmegaLy*LNhz);
  *dLNhy = (-OmegaLx*LNhz + OmegaLz*LNhx);
  *dLNhz = (-OmegaLy*LNhx + OmegaLx*LNhy);

  /*
   * dE1
   *
   * d E_1 / d \hat{t} = M * d E_1 / dt
   * Computed from \Omega_L and \hat{L_N} with Eq. (15)-(16) of gr-qc/0310034
   */
  REAL8 OmegaLdotLN = OmegaLx * LNhx + OmegaLy * LNhy + OmegaLz * LNhz;
  /* \Omega_E vector */
  REAL8 OmegaEx = OmegaLx - OmegaLdotLN * LNhx;
  REAL8 OmegaEy = OmegaLy - OmegaLdotLN * LNhy;
  REAL8 OmegaEz = OmegaLz - OmegaLdotLN * LNhz;

  /* Take cross product of \Omega_E with E_1 */
  *dE1x = (-OmegaEz*E1y + OmegaEy*E1z);
  *dE1y = (-OmegaEx*E1z + OmegaEz*E1x);
  *dE1z = (-OmegaEy*E1x + OmegaEx*E1y);

  /*
   * dS1
   * d S_1 / d \hat{t} = M * d S_1 / dt = \Omega_{S1,S2,LN,v} x S_1
   * However, that paper uses spin variables which are M^2 times our spins
   */


  /* dS1, 1.5PN: eq. (8) of gr-qc/0405090.
   */
  REAL8 cross1x = (LNhy * S1z - LNhz * S1y);
  REAL8 cross1y = (LNhz * S1x - LNhx * S1z);
  REAL8 cross1z = (LNhx * S1y - LNhy * S1x);

  *dS1x = params->S1dot3 * v5 * cross1x;
  *dS1y = params->S1dot3 * v5 * cross1y;
  *dS1z = params->S1dot3 * v5 * cross1z;

  /* dS1, 2PN */
  REAL8 tmpx = S1z * S2y - S1y * S2z;
  REAL8 tmpy = S1x * S2z - S1z * S2x;
  REAL8 tmpz = S1y * S2x - S1x * S2y;

  /* S1S2 contribution
   * see. eq. 2.23 of arXiv:0812.4413
   */
  *dS1x += omega2 * (params->Sdot4S1S2*tmpx + params->Sdot4S1OS2O * LNhdotS2 * cross1x);
  *dS1y += omega2 * (params->Sdot4S1S2*tmpy + params->Sdot4S1OS2O * LNhdotS2 * cross1y);
  *dS1z += omega2 * (params->Sdot4S1S2*tmpz + params->Sdot4S1OS2O * LNhdotS2 * cross1z);
  /* S1S1 contribution
   */
  *dS1x += omega2 * LNhdotS1 * cross1x * params->S1dot4QMS1O;
  *dS1y += omega2 * LNhdotS1 * cross1y * params->S1dot4QMS1O;
  *dS1z += omega2 * LNhdotS1 * cross1z * params->S1dot4QMS1O;

  /* dS1, 2.5PN
   * eq. 7.8 of Blanchet et al. gr-qc/0605140
   */
  *dS1x += params->S1dot5 * v7 * cross1x;
  *dS1y += params->S1dot5 * v7 * cross1y;
  *dS1z += params->S1dot5 * v7 * cross1z;

  /* dS2, 1.5PN */
  REAL8 cross2x = (LNhy * S2z - LNhz * S2y);
  REAL8 cross2y = (LNhz * S2x - LNhx * S2z);
  REAL8 cross2z = (LNhx * S2y - LNhy * S2x);

  *dS2x = params->S2dot3 * v5 * cross2x;
  *dS2y = params->S2dot3 * v5 * cross2y;
  *dS2z = params->S2dot3 * v5 * cross2z;

  /* dS2, 2PN */
  *dS2x += omega2 * (-params->Sdot4S1S2*tmpx + params->Sdot4S1OS2O * LNhdotS1 * cross2x);
  *dS2y += omega2 * (-params->Sdot4S1S2*tmpy + params->Sdot4S1OS2O * LNhdotS1 * cross2y);
  *dS2z += omega2 * (-params->Sdot4S1S2*tmpz + params->Sdot4S1OS2O * LNhdotS1 * cross2z);
  // S2S2 contribution
  *dS2x += omega2 * LNhdotS2 * cross2x * params->S2dot4QMS2O;
  *dS2y += omega2 * LNhdotS2 * cross2y * params->S2dot4QMS2O;
  *dS2z += omega2 * LNhdotS2 * cross2z * params->S2dot4QMS2O;

  // dS2, 2.5PN, eq. 7.8 of Blanchet et al. gr-qc/0605140
  *dS2x += params->S2dot5 * v7 * cross2x;
  *dS2y += params->S2dot5 * v7 * cross2y;
  *dS2z += params->S2dot5 * v7 * cross2z;

  return XLAL_SUCCESS;
}

static int XLALSimInspiralSpinTaylorT2Setup(
    XLALSimInspiralSpinTaylorTxCoeffs *params, /**< UNDOCUMENTED */
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
    REAL8 M, eta, Mchirp, m1M, m2M;
    /* Zero the coefficients */
    memset(params, 0, sizeof(XLALSimInspiralSpinTaylorTxCoeffs));

    /* Define mass variables and other coefficients */
    m1 *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert m1 from kg to seconds */
    m2 *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert m2 from kg to seconds */
    M = m1 + m2;
    m1M = m1 / M;
    m2M = m2 / M;
    eta = m1 * m2 / M / M;
    Mchirp = M * pow(eta, 3./5.);
    params->wdotnewt = 3./XLALSimInspiralTaylorT2dtdv_0PNCoeff(eta);//(96.0/5.0) * eta;
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
	  params->wdotcoeff[7] = XLALSimInspiralTaylorT2dtdv_7PNCoeff(eta);
	  params->Ecoeff[7] = 0.;
        /* case LAL_PNORDER_THREE: */
        case 6:
            params->wdotcoeff[6] = XLALSimInspiralTaylorT2dtdv_6PNCoeff(eta);
	    //params->wdotcoeff[6] = 22.065 + 165.416*eta
            //        - 2.20067*eta*eta + 4.93152*eta*eta*eta;
            params->wdotlogcoeff = XLALSimInspiralTaylorT2dtdv_6PNLogCoeff(eta);
            params->Ecoeff[6] = XLALSimInspiralPNEnergy_6PNCoeff(eta);
        /* case LAL_PNORDER_TWO_POINT_FIVE: */
        case 5:
            params->wdotcoeff[5] = XLALSimInspiralTaylorT2dtdv_5PNCoeff(eta);
            params->Ecoeff[5] = 0.;
        /* case LAL_PNORDER_TWO: */
        case 4:
	    params->wdotcoeff[4] = XLALSimInspiralTaylorT2dtdv_4PNCoeff(eta);
            params->Ecoeff[4] = XLALSimInspiralPNEnergy_4PNCoeff(eta);
        /*case LAL_PNORDER_ONE_POINT_FIVE:*/
        case 3:
	    params->wdotcoeff[3] = XLALSimInspiralTaylorT2dtdv_3PNCoeff(eta);
            params->Ecoeff[3] = 0.;
        /*case LAL_PNORDER_ONE:*/
        case 2:
            params->wdotcoeff[2] = XLALSimInspiralTaylorT2dtdv_2PNCoeff(eta);
            params->Ecoeff[2] = XLALSimInspiralPNEnergy_2PNCoeff(eta);
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
            XLALPrintError("XLAL Error - %s: Invalid phase. PN order %d\n",
                    __func__, phaseO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    /* Compute the non-dynamical coefficients of spin corrections
     * to the evolution equations for omega, L, S1 and S2 and binary energy E.
     */
    switch( spinO )
    {   // case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:

        case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
            params->wdot7S1O = XLALSimInspiralTaylorT2dtdv_7PNSOCoeff(m1M);
            params->wdot7S2O = XLALSimInspiralTaylorT2dtdv_7PNSOCoeff(m2M);
            params->E7S1O = XLALSimInspiralPNEnergy_7PNSOCoeff(m1M);
            params->E7S2O = XLALSimInspiralPNEnergy_7PNSOCoeff(m1M);
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
            params->wdot6S1O = XLALSimInspiralTaylorT2dtdv_6PNSOCoeff(m1M);
            params->wdot6S2O = XLALSimInspiralTaylorT2dtdv_6PNSOCoeff(m2M);
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
            params->wdot5S1O = XLALSimInspiralTaylorT2dtdv_5PNSOCoeff(m1M);
            params->wdot5S2O = XLALSimInspiralTaylorT2dtdv_5PNSOCoeff(m2M);
            params->E5S1O  = XLALSimInspiralPNEnergy_5PNSOCoeff(m1M);
            params->E5S2O  = XLALSimInspiralPNEnergy_5PNSOCoeff(m2M);
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            // 2PN spin1-spin2 terms
            params->wdot4S1S2 	= XLALSimInspiralTaylorT2dtdv_4PNS1S2Coeff(eta);
            params->wdot4S1OS2O = XLALSimInspiralTaylorT2dtdv_4PNS1S2OCoeff(eta);
            // 2PN spin-self^2 terms
            params->wdot4S1S1   = XLALSimInspiralTaylorT2dtdv_4PNSelfSSCoeff(m1M);
            params->wdot4S1OS1O = XLALSimInspiralTaylorT2dtdv_4PNSelfSSOCoeff(m1M);
            params->wdot4S2S2   = XLALSimInspiralTaylorT2dtdv_4PNSelfSSCoeff(m2M);
            params->wdot4S2OS2O = XLALSimInspiralTaylorT2dtdv_4PNSelfSSOCoeff(m2M);
            params->E4S1S2     = XLALSimInspiralPNEnergy_4PNS1S2Coeff(eta);
            params->E4S1OS2O   = XLALSimInspiralPNEnergy_4PNS1S2OCoeff(eta);
            // 2PN quadrupole-monopole self spin terms
            params->wdot4QMS1S1   = quadparam1 * XLALSimInspiralTaylorT2dtdv_4PNQMCoeff(m1M);
            params->wdot4QMS1OS1O = quadparam1 * XLALSimInspiralTaylorT2dtdv_4PNQMSOCoeff(m1M);
            params->wdot4QMS2S2   = quadparam2 * XLALSimInspiralTaylorT2dtdv_4PNQMCoeff(m2M);
            params->wdot4QMS2OS2O = quadparam2 * XLALSimInspiralTaylorT2dtdv_4PNQMSOCoeff(m2M);
            params->E4QMS1S1      = quadparam1 * XLALSimInspiralPNEnergy_4PNQM2SCoeff(m1M);
            params->E4QMS1OS1O 	  = quadparam1 * XLALSimInspiralPNEnergy_4PNQM2SOCoeff(m1M);
            params->E4QMS2S2 	  = quadparam2 * XLALSimInspiralPNEnergy_4PNQM2SCoeff(m2M);
            params->E4QMS2OS2O 	  = quadparam2 * XLALSimInspiralPNEnergy_4PNQM2SOCoeff(m2M);
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            // Note: LNHat do not have their signs reversed relative to T4
            // They are precession rather than orbital quantities
	    //There was a wrong sign in wdotSO15s12!
	    params->wdot3S1O 	= XLALSimInspiralTaylorT2dtdv_3PNSOCoeff(m1M);
            params->wdot3S2O 	= XLALSimInspiralTaylorT2dtdv_3PNSOCoeff(m2M);
	    params->E3S1O      = XLALSimInspiralPNEnergy_3PNSOCoeff(m1M);
            params->E3S2O      = XLALSimInspiralPNEnergy_3PNSOCoeff(m2M);
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

    XLALSimInspiralSpinTaylorSpinDerivativeSetup(params,m1,m2,quadparam1,quadparam2,spinO);

    /*
     * Compute the coefficients of tidal corrections
     * to the evolution equations for omega and binary energy E.
     * Coefficients found from Eqs. 2.11 and 3.10 of
     * Vines, Flanagan, Hinderer, PRD 83, 084051 (2011).
     */
    switch( tideO )
    {
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL:
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN:
	    params->wdottidal6pn = lambda1 * XLALSimInspiralTaylorT2dtdv_12PNTidalCoeff(m1M)
	                         + lambda2 * XLALSimInspiralTaylorT2dtdv_12PNTidalCoeff(m2M);
	    params->Etidal6pn = lambda1 * XLALSimInspiralPNEnergy_12PNTidalCoeff(m1M)
	                      + lambda2 * XLALSimInspiralPNEnergy_12PNTidalCoeff(m2M);
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
	    params->wdottidal5pn = lambda1 * XLALSimInspiralTaylorT2dtdv_10PNTidalCoeff(m1M)
	                         + lambda2 * XLALSimInspiralTaylorT2dtdv_10PNTidalCoeff(m2M);
            params->Etidal5pn = lambda1 * XLALSimInspiralPNEnergy_10PNTidalCoeff(m1M)
                              + lambda2 * XLALSimInspiralPNEnergy_10PNTidalCoeff(m2M);
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid tidal PN order %d\n",
                    __func__, tideO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    return XLAL_SUCCESS;
}

static int XLALSimInspiralSpinTaylorT4Setup(
    XLALSimInspiralSpinTaylorTxCoeffs *params, /**< UNDOCUMENTED */
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
    REAL8 M, eta, Mchirp, m1M, m2M;
    /* Zero the coefficients */
    memset(params, 0, sizeof(XLALSimInspiralSpinTaylorTxCoeffs));

    /* Define mass variables and other coefficients */
    m1 *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert m1 from kg to seconds */
    m2 *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert m2 from kg to seconds */
    M = m1 + m2;
    m1M = m1 / M;
    m2M = m2 / M;
    eta = m1 * m2 / M / M;
    Mchirp = M * pow(eta, 3./5.);
    params->wdotnewt = XLALSimInspiralTaylorT4wdot_0PNCoeff(eta);
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
            params->wdotcoeff[7] = XLALSimInspiralTaylorT4wdot_7PNCoeff(eta);
            params->Ecoeff[7] = 0.;
        /* case LAL_PNORDER_THREE: */
        case 6:
            params->wdotcoeff[6] = XLALSimInspiralTaylorT4wdot_6PNCoeff(eta);
            params->wdotlogcoeff = XLALSimInspiralTaylorT4wdot_6PNLogCoeff(eta);
            params->Ecoeff[6] = XLALSimInspiralPNEnergy_6PNCoeff(eta);
        /* case LAL_PNORDER_TWO_POINT_FIVE: */
        case 5:
            params->wdotcoeff[5] = XLALSimInspiralTaylorT4wdot_5PNCoeff(eta);
            params->Ecoeff[5] = 0.;
        /* case LAL_PNORDER_TWO: */
        case 4:
            params->wdotcoeff[4] = XLALSimInspiralTaylorT4wdot_4PNCoeff(eta);
            params->Ecoeff[4] = XLALSimInspiralPNEnergy_4PNCoeff(eta);
        /*case LAL_PNORDER_ONE_POINT_FIVE:*/
        case 3:

            params->wdotcoeff[3] = XLALSimInspiralTaylorT4wdot_3PNCoeff(eta);
            params->Ecoeff[3] = 0.;
        /*case LAL_PNORDER_ONE:*/
        case 2:
            params->wdotcoeff[2] = XLALSimInspiralTaylorT4wdot_2PNCoeff(eta);
            params->Ecoeff[2] = XLALSimInspiralPNEnergy_2PNCoeff(eta);
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
            XLALPrintError("XLAL Error - %s: Invalid phase. PN order %d\n",
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
	    params->wdot7S1O = XLALSimInspiralTaylorT4wdot_7PNSOCoeff(m1M);
            params->wdot7S2O = XLALSimInspiralTaylorT4wdot_7PNSOCoeff(m2M);
            params->E7S1O = XLALSimInspiralPNEnergy_7PNSOCoeff(m1M);
            params->E7S2O = XLALSimInspiralPNEnergy_7PNSOCoeff(m2M);
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
            params->wdot6S1O = XLALSimInspiralTaylorT4wdot_6PNSOCoeff(m1M);
            params->wdot6S2O = XLALSimInspiralTaylorT4wdot_6PNSOCoeff(m2M);
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
            params->wdot5S1O = XLALSimInspiralTaylorT4wdot_5PNSOCoeff(m1M);
            params->wdot5S2O = XLALSimInspiralTaylorT4wdot_5PNSOCoeff(m2M);
            params->E5S1O = XLALSimInspiralPNEnergy_5PNSOCoeff(m1M);
            params->E5S2O = XLALSimInspiralPNEnergy_5PNSOCoeff(m2M);
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            // 2PN spin-spin terms
            params->wdot4S1S2     = XLALSimInspiralTaylorT4wdot_4PNS1S2Coeff(eta);
            params->wdot4S1OS2O   = XLALSimInspiralTaylorT4wdot_4PNS1S2OCoeff(eta);
            params->E4S1S2 	= XLALSimInspiralPNEnergy_4PNS1S2Coeff(eta);
            params->E4S1OS2O    = XLALSimInspiralPNEnergy_4PNS1S2OCoeff(eta);
            // 2PN self-spin terms
            params->wdot4S1S1     = XLALSimInspiralTaylorT4wdot_4PNSelfSSCoeff(m1M);
            params->wdot4S1OS1O   = XLALSimInspiralTaylorT4wdot_4PNSelfSSOCoeff(m1M);
            params->wdot4S2S2     = XLALSimInspiralTaylorT4wdot_4PNSelfSSCoeff(m2M);
            params->wdot4S2OS2O   = XLALSimInspiralTaylorT4wdot_4PNSelfSSOCoeff(m2M);
            // 2PN quadrupole-monopole terms
            params->wdot4QMS1S1   = quadparam1 * XLALSimInspiralTaylorT4wdot_4PNQMCoeff(m1M);
            params->wdot4QMS1OS1O = quadparam1 * XLALSimInspiralTaylorT4wdot_4PNQMSOCoeff(m1M);
            params->wdot4QMS2S2   = quadparam2 *  XLALSimInspiralTaylorT4wdot_4PNQMCoeff(m2M);
            params->wdot4QMS2OS2O = quadparam2 *  XLALSimInspiralTaylorT4wdot_4PNQMSOCoeff(m2M);
            params->E4QMS1S1      = quadparam1 * XLALSimInspiralPNEnergy_4PNQM2SCoeff(m1M);
            params->E4QMS1OS1O    = quadparam1 * XLALSimInspiralPNEnergy_4PNQM2SOCoeff(m1M);
            params->E4QMS2S2      = quadparam2 * XLALSimInspiralPNEnergy_4PNQM2SCoeff(m2M);
            params->E4QMS2OS2O 	  = quadparam2 * XLALSimInspiralPNEnergy_4PNQM2SOCoeff(m2M);

        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
	    params->wdot3S1O 	= XLALSimInspiralTaylorT4wdot_3PNSOCoeff(m1M);
            params->wdot3S2O 	= XLALSimInspiralTaylorT4wdot_3PNSOCoeff(m2M);
            params->E3S1O 	= XLALSimInspiralPNEnergy_3PNSOCoeff(m1M);
            params->E3S2O 	= XLALSimInspiralPNEnergy_3PNSOCoeff(m2M);
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

    XLALSimInspiralSpinTaylorSpinDerivativeSetup(params,m1,m2,quadparam1,quadparam2,spinO);
	
    /* Compute the coefficients of tidal corrections
     * to the evolution equations for omega and binary energy E.
     * Coefficients found from Eqs. 2.11 and 3.10 of 
     * Vines, Flanagan, Hinderer, PRD 83, 084051 (2011)
     * and eq. 2 and 3 of PhysRevD.89.103012.
     */
    switch( tideO )
    {
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL:
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN:
            params->wdottidal6pn = lambda1 * XLALSimInspiralTaylorT4wdot_12PNTidalCoeff(m1M) + lambda2 * XLALSimInspiralTaylorT4wdot_12PNTidalCoeff(m2M);
            params->Etidal6pn =  lambda1*XLALSimInspiralPNEnergy_12PNTidalCoeff(m1M) + lambda2*XLALSimInspiralPNEnergy_12PNTidalCoeff(m2M);
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
	    params->wdottidal5pn = lambda1 * XLALSimInspiralTaylorT4wdot_10PNTidalCoeff(m1M) + lambda2 * XLALSimInspiralTaylorT4wdot_10PNTidalCoeff(m2M);
	    params->Etidal5pn = lambda1*XLALSimInspiralPNEnergy_10PNTidalCoeff(m1M) + lambda2*XLALSimInspiralPNEnergy_10PNTidalCoeff(m2M);
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid tidal PN order %d\n",
                    __func__, tideO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    return XLAL_SUCCESS;
}


static int XLALSimInspiralSpinTaylorT1Setup(
    XLALSimInspiralSpinTaylorTxCoeffs *params, /**< UNDOCUMENTED */
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
    REAL8 M, eta, Mchirp, m1M, m2M;
    /* Zero the coefficients */
    memset(params, 0, sizeof(XLALSimInspiralSpinTaylorTxCoeffs));

    /* Define mass variables and other coefficients */
    m1 *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert m1 from kg to seconds */
    m2 *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert m2 from kg to seconds */
    M = m1 + m2;
    m1M = m1 / M;
    m2M = m2 / M;
    eta = m1 * m2 / M / M;
    Mchirp = M * pow(eta, 3./5.);
    params->dEdvnewt = 2.*XLALSimInspiralPNEnergy_0PNCoeff(eta);
    params->Fnewt = XLALSimInspiralPNFlux_0PNCoeff(eta);
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
            params->Fcoeff[7] = XLALSimInspiralPNFlux_7PNCoeff(eta);
            params->Ecoeff[7] = 0.;
        /* case LAL_PNORDER_THREE: */
        case 6:
            params->Fcoeff[6] = XLALSimInspiralPNFlux_6PNCoeff(eta);
	    params->Flogcoeff = XLALSimInspiralPNFlux_6PNLogCoeff(eta);
            params->Ecoeff[6] = XLALSimInspiralPNEnergy_6PNCoeff(eta);
        /* case LAL_PNORDER_TWO_POINT_FIVE: */
        case 5:
	    params->Fcoeff[5] = XLALSimInspiralPNFlux_5PNCoeff(eta);
            params->Ecoeff[5] = 0.;
        /* case LAL_PNORDER_TWO: */
        case 4:
            params->Fcoeff[4] = XLALSimInspiralPNFlux_4PNCoeff(eta);
	    params->Ecoeff[4] = XLALSimInspiralPNEnergy_4PNCoeff(eta);
        /*case LAL_PNORDER_ONE_POINT_FIVE:*/
        case 3:
            params->Fcoeff[3] = XLALSimInspiralPNFlux_3PNCoeff(eta);
	    params->Ecoeff[3] = 0.;
        /*case LAL_PNORDER_ONE:*/
        case 2:
	    params->Fcoeff[2] = XLALSimInspiralPNFlux_2PNCoeff(eta);
            params->Ecoeff[2] = XLALSimInspiralPNEnergy_2PNCoeff(eta);
        /*case LAL_PNORDER_HALF:*/
        case 1:
	    params->Fcoeff[1] = 0.;
            params->Ecoeff[1] = 0.;
        /*case LAL_PNORDER_NEWTONIAN:*/
        case 0:
	    params->Fcoeff[0] = 1.;
            params->Ecoeff[0] = 1.;
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid phase. PN order %d\n",
                    __func__, phaseO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    switch( spinO )
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
	case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
            params->E7S1O = XLALSimInspiralPNEnergy_7PNSOCoeff(m1M);
	    params->E7S2O = XLALSimInspiralPNEnergy_7PNSOCoeff(m2M);
	    params->F7S1O=XLALSimInspiralPNFlux_7PNSOCoeff(m1M);
	    params->F7S2O=XLALSimInspiralPNFlux_7PNSOCoeff(m2M);
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
	    params->F6S1O = XLALSimInspiralPNFlux_6PNSOCoeff(m1M);
	    params->F6S2O = XLALSimInspiralPNFlux_6PNSOCoeff(m2M);
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
            params->E5S1O = XLALSimInspiralPNEnergy_5PNSOCoeff(m1M);
            params->E5S2O = XLALSimInspiralPNEnergy_5PNSOCoeff(m2M);
	    params->F5S1O =XLALSimInspiralPNFlux_5PNSOCoeff(m1M);
	    params->F5S2O =XLALSimInspiralPNFlux_5PNSOCoeff(m2M);
	case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            // 2PN spin-spin terms
            params->E4S1S2 	= XLALSimInspiralPNEnergy_4PNS1S2Coeff(eta);
            params->E4S1OS2O    = XLALSimInspiralPNEnergy_4PNS1S2OCoeff(eta);
	    params->F4S1S2      = XLALSimInspiralPNFlux_4PNS1S2Coeff(eta);
	    params->F4S1OS2O    = XLALSimInspiralPNFlux_4PNS1S2OCoeff(eta);
            // 2PN self-spin terms
            params->F4S1S1     += XLALSimInspiralPNFlux_4PNSelf2SCoeff(m1M);
            params->F4S1OS1O   += XLALSimInspiralPNFlux_4PNSelf2SOCoeff(m1M);
            params->F4S2S2     += XLALSimInspiralPNFlux_4PNSelf2SCoeff(m2M);
            params->F4S2OS2O   += XLALSimInspiralPNFlux_4PNSelf2SOCoeff(m2M);
	    // 2PN quadrupole-monopole terms
            params->E4QMS1S1    = quadparam1*XLALSimInspiralPNEnergy_4PNQM2SCoeff(m1M);
            params->E4QMS1OS1O 	= quadparam1*XLALSimInspiralPNEnergy_4PNQM2SOCoeff(m1M);
            params->E4QMS2S2 	= quadparam2*XLALSimInspiralPNEnergy_4PNQM2SCoeff(m2M);
            params->E4QMS2OS2O 	= quadparam2*XLALSimInspiralPNEnergy_4PNQM2SOCoeff(m2M);
	    params->F4QMS1S1    = quadparam1*XLALSimInspiralPNFlux_4PNQM2SCoeff(m1M);
            params->F4QMS1OS1O  = quadparam1*XLALSimInspiralPNFlux_4PNQM2SOCoeff(m1M);
            params->F4QMS2S2    = quadparam2*XLALSimInspiralPNFlux_4PNQM2SCoeff(m2M);
            params->F4QMS2OS2O  = quadparam2*XLALSimInspiralPNFlux_4PNQM2SOCoeff(m2M);
    case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            params->E3S1O     = XLALSimInspiralPNEnergy_3PNSOCoeff(m1M);
            params->E3S2O     = XLALSimInspiralPNEnergy_3PNSOCoeff(m2M);
	    params->F3S1O      = XLALSimInspiralPNFlux_3PNSOCoeff(m1M);
            params->F3S2O      = XLALSimInspiralPNFlux_3PNSOCoeff(m2M);
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

    XLALSimInspiralSpinTaylorSpinDerivativeSetup(params,m1,m2,quadparam1,quadparam2,spinO);

    /* Compute the coefficients of tidal corrections
     * to the evolution equations for omega and binary energy E.
     * Coefficients found from Eqs. 2.11 and 3.10 of
     * Vines, Flanagan, Hinderer, PRD 83, 084051 (2011).
     */
    switch( tideO )
    {
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL:
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN:
            params->Ftidal6pn = lambda1 * XLALSimInspiralPNFlux_12PNTidalCoeff(m1M) + lambda2 * XLALSimInspiralPNFlux_12PNTidalCoeff(m2M);
            params->Etidal6pn =  lambda1*XLALSimInspiralPNEnergy_12PNTidalCoeff(m1M) + lambda2*XLALSimInspiralPNEnergy_12PNTidalCoeff(m2M);
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
	    params->Ftidal5pn = lambda1 * XLALSimInspiralPNFlux_10PNTidalCoeff(m1M) + lambda2 * XLALSimInspiralPNFlux_10PNTidalCoeff(m2M);
	    params->Etidal5pn = lambda1*XLALSimInspiralPNEnergy_10PNTidalCoeff(m1M) + lambda2*XLALSimInspiralPNEnergy_10PNTidalCoeff(m2M);
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid tidal PN order %d\n",
                    __func__, tideO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    return XLAL_SUCCESS;
}


/*
 * Internal driver function to generate any of SpinTaylorT2/T4 in the Fourier domain.
 */
static int XLALSimInspiralSpinTaylorDriverFourier(
        COMPLEX16FrequencySeries **hplus,        /**< +-polarization waveform */
        COMPLEX16FrequencySeries **hcross,       /**< x-polarization waveform */
        REAL8 fMin,                     /**< minimum frequency of the returned series */
        REAL8 fMax,                     /**< maximum frequency of the returned series */
        REAL8 deltaF,                   /**< frequency interval of the returned series */
        INT4 kMax,                      /**< k_max as described in arXiv: 1408.5158 (min 0, max 10). */
        REAL8 phiRef,                   /**< orbital phase at reference pt. */
        REAL8 v0,                       /**< tail gauge term (default = 1) */
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
        INT4 phaseO,                    /**< twice PN phase order */
        INT4 amplitudeO,                /**< twice PN amplitude order */
        Approximant approx,             /**< PN approximant (SpinTaylorT1/T2/T4) */
        INT4 phiRefAtEnd                /**< whether phiRef corresponds to the end of the inspiral */
        )
{
    REAL8Array* y = 0;
    INT4 n;
    UINT4 i, j;
    UINT4 iRefPhi;
    REAL8 fS, fE, phiShift;
    /* The Schwarzschild ISCO frequency - for sanity checking fRef, and default maximum frequency */
    REAL8 fISCO = pow(LAL_C_SI,3) / (pow(6.,3./2.)*LAL_PI*(m1+m2)*LAL_G_SI);


    if(kMax < 0)
    {
        XLALPrintError("XLAL Error - %s: kMax = %d must be >= 0\n",
                __func__, kMax);
        XLAL_ERROR(XLAL_EINVAL);
    }

    if(kMax > 10)
    {
        XLALPrintError("XLAL Error - %s: kMax = %d not implemented. Must be <= 10\n",
                __func__, kMax);
        XLAL_ERROR(XLAL_EINVAL);
    }



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
    if( fRef < LAL_REAL4_EPS)
    {
        fS = fStart;
        fE = 0.;
        /* Evolve the dynamical variables */
        n = XLALSimInspiralSpinTaylorPNEvolveOrbitIrregularIntervals(&y,
                m1, m2, fS, fE, s1x, s1y, s1z, s2x, s2y, s2z,
                lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2,
                quadparam1, quadparam2, spinO, tideO, phaseO, approx);
        if( n < 0 )
            XLAL_ERROR(XLAL_EFUNC);

        /* Apply phase shift so orbital phase ends with desired value */
        iRefPhi = ((y->dimLength->data[1])<<1)-1;
        phiShift = phiRef - y->data[iRefPhi];
        for( i = y->dimLength->data[1]; i <= iRefPhi; i++)
        {
            y->data[i] += phiShift;
        }
    }
    /* if fRef=fStart, just integrate from start to end. Let phiRef=phiStart */
    else if( fabs(fRef - fStart) < LAL_REAL4_EPS )
    {
        fS = fStart;
        fE = 0.;
        /* Evolve the dynamical variables */
        n = XLALSimInspiralSpinTaylorPNEvolveOrbitIrregularIntervals(&y,
                m1, m2, fS, fE, s1x, s1y, s1z, s2x, s2y, s2z,
                lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2,
                quadparam1, quadparam2, spinO, tideO, phaseO, approx);
        if( n < 0 )
            XLAL_ERROR(XLAL_EFUNC);

        /* Apply phase shift so orbital phase starts with desired value */
        iRefPhi = ((y->dimLength->data[1])<<1)-1;
        phiShift = phiRef - y->data[y->dimLength->data[1]];
        for( i = y->dimLength->data[1]; i <= iRefPhi; i++)
        {
            y->data[i] += phiShift;
        }
    }
    else /* Start in middle, integrate backward and forward, stitch together */
    {
        REAL8Array *yStart=NULL, *yEnd=NULL;

        /* Integrate backward to fStart */
        fS = fRef;
        fE = fStart;
        n = XLALSimInspiralSpinTaylorPNEvolveOrbitIrregularIntervals(&yStart,
                m1, m2, fS, fE, s1x, s1y, s1z, s2x, s2y,
                s2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2,
                quadparam1, quadparam2, spinO, tideO, phaseO, approx);

        /* Apply phase shift so orbital phase has desired value at fRef */
        iRefPhi = ((yStart->dimLength->data[1])<<1)-1;
        phiShift = phiRef - yStart->data[iRefPhi];
        for( i = yStart->dimLength->data[1]; i <= iRefPhi; i++)
        {
            yStart->data[i] += phiShift;
        }

        /* Integrate forward to end of waveform */
        fS = fRef;
        fE = 0.;
        n = XLALSimInspiralSpinTaylorPNEvolveOrbitIrregularIntervals(&yEnd,
                m1, m2, fS, fE, s1x, s1y, s1z, s2x, s2y,
                s2z, lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2,
                quadparam1, quadparam2, spinO, tideO, phaseO, approx);

        /* Apply phase shift so orbital phase has desired value at fRef */
        iRefPhi = ((yEnd->dimLength->data[1])<<1)-1;
        phiShift = phiRef - yEnd->data[yEnd->dimLength->data[1]];
        for( i = yEnd->dimLength->data[1]; i <= iRefPhi; i++)
        {
            yEnd->data[i] += phiShift;
        }

        /* Stitch 2nd set of vectors onto 1st set. Free 'em. */
        y = appendTAandFree(yStart, yEnd, y);
    }

    /* if phiRefAtEnd, let phiRef=phiC (to be compatible with TaylorT4 in ChooseFDWaveform). */
    if(phiRefAtEnd)
    {
      iRefPhi = ((y->dimLength->data[1])<<1)-1;
      phiShift = phiRef - y->data[iRefPhi];
      for( i = y->dimLength->data[1]; i <= iRefPhi; i++)
      {
        y->data[i] += phiShift;
      }
    }

    // We can now use the interpolation formulas to create the Fourier series

    REAL8 cSiDbl, Msec, M;

    cSiDbl = LAL_C_SI;
    M = m1 + m2;
    Msec = M*LAL_G_SI/(cSiDbl*cSiDbl*cSiDbl);

    UINT4 length = y->dimLength->data[1];
    UINT4 nParams = y->dimLength->data[0];

    /** minimum index for the a_k's adefined in Eq. (40) of arXiv: 1408.5158. */
    const INT4 iMinKMax[11] = {
      0,
      1,
      3,
      6,
      10,
      15,
      21,
      28,
      36,
      45,
      55
    };

    /** constants a_k defined in Eq. (40) of arXiv: 1408.5158. */
    const COMPLEX16 akCsts[66] = {
      1. + I*(0.), // 0
      1. + I*(1.), // 1
      0. + I*(- 1./2.),
      1./4. + I*(5./4.), // 2
      1./2. + I*(-2./3.),
      -1./8. + I*(1./24.),
      -1./6. + I*(17./18.), // 3
      13./16. + I*(-7./16.),
      -1./4. + I*(-1./20.),
      1./48. + I*(11./720.),
      -23./96. + I*(185./288.), // 4
      209./240. + I*(-47./240.),
      -67./240. + I*(-41./240.),
      7./240. + I*(251./5040.),
      -1./960. + I*(-29./6720.),
      -23./120. + I*(1669./3600.), // 5
      2393./2880. + I*(-3./64.),
      -323./1260. + I*(-43./168.),
      277./13440. + I*(659./8064.),
      13./15120. + I*(-23./2016.),
      -23./120960. + I*(143./201600.),
      -683./5400. + I*(16003./43200.), // 6
      26041./33600. + I*(19./576.),
      -31./140. + I*(-4933./16128.),
      1847./362880. + I*(7541./72576.),
      139./25200. + I*(-437./24192.),
      -209./201600. + I*(12769./6652800.),
      1./14175. + I*(-23./228096.),
      -2153./30240. + I*(341101./1058400.), // 7
      878963./1209600. + I*(36349./483840.),
      -343247./1814400. + I*(-121187./362880.),
      -40043./3628800. + I*(171209./1451520.),
      113557./9979200. + I*(-46247./1995840.),
      -19979./7983360. + I*(255173./79833600.),
      5909./19958400. + I*(-15427./51891840.),
      -643./39916800. + I*(20389./1452971520.),
      -1940303./67737600. + I*(6696113./22579200.), // 8
      13123063./19051200. + I*(426691./4354560.),
      -17712403./108864000. + I*(-1523929./4354560.),
      -2542391./99792000. + I*(20431./161280.),
      8333153./479001600. + I*(-2568281./95800320.),
      -3389161./778377600. + I*(13435997./3113510400.),
      5033293./7264857600. + I*(-4129./7687680.),
      -235007./3405402000. + I*(200537./4358914560.),
      2882417./871782912000. + I*(-181./90574848.),
      662159./182891520. + I*(387137501./1371686400.), // 9
      4022350847./6096384000. + I*(15021619./135475200.),
      -593412859./4191264000. + I*(-16734271./46569600.),
      -181044539./4790016000. + I*(379801511./2874009600.),
      359342129./15567552000. + I*(-461261./15724800.),
      -778483./121927680. + I*(50534819./9686476800.),
      6114053./4953312000. + I*(-12709559./16345929600.),
      -297352897./1743565824000. + I*(320849./3522355200.),
      511381./33530112000. + I*(-601231./82335052800.),
      -3471163./5230697472000. + I*(15721091./53353114214400.),
      1110580187./39191040000. + I*(1503012587./5486745600.), // 10
      76930217647./120708403200. + I*(529444847./4470681600.),
      -10037362331./80472268800. + I*(-543929527./1490227200.),
      -12597360953./261534873600. + I*(35472685841./261534873600.),
      7393112329./261534873600. + I*(-50167081./1614412800.),
      -55307796493./6538371840000. + I*(578597./97843200.),
      3938740549./2092278988800. + I*(-94645351./95103590400.),
      -820761239./2540624486400. + I*(93730183./658680422400.),
      2169874009./53353114214400. + I*(-15647627./988020633600.),
      -178160027./53353114214400. + I*(1209449069/1013709170073600.),
      356885411./2667655710720000. + I*(-1686571./37544784076800.)
    };


    INT4 status;


    REAL8 dm = (m1 - m2)/M;
    REAL8 eta = m1*m2/(M*M);

    // Construct cubic spline interpolation structures
    gsl_spline *interp[nParams-1];
    gsl_interp_accel *accel = NULL;
    gsl_interp_accel *accel_k = NULL;
    gsl_spline* interp_t_of_f;


    INT4 kMaxCur = kMax;
    INT4 k, iK;

    REAL8 tCur;
    REAL8 fCur;
    REAL8 omegaHatCur;

    REAL8 freqToOmegaHat = LAL_TWOPI*Msec;

    REAL8 dydt, T, phiOrb;

    REAL8 tK;

    REAL8 prefac = 2.*sqrt(LAL_TWOPI)*Msec*Msec*eta*cSiDbl/r;

    REAL8 vCur, lnhatxCur, lnhatyCur, lnhatzCur, s1xCur, s1yCur, s1zCur, s2xCur,
    s2yCur, s2zCur, E1xCur, E1yCur, E1zCur, TCur;

    INT4 nHarmonic, minHarmonic, maxHarmonic;

    COMPLEX16 htPlusHarmonic, htCrossHarmonic, htPlusTemp, htCrossTemp, cPhase;

    // compute ISCO orbital frequency
    REAL8 omegaHatISCO = 1./sqrt(216.);
    REAL8 tISCO;

    // Set the minimum and maximum harmonic needed for the amplitude PN order
    switch(amplitudeO)
    {
      case 0:
        minHarmonic = 2;
        maxHarmonic = 2;
        break;
      case 1:
        minHarmonic = 1;
        maxHarmonic = 3;
        break;
      case 2:
        minHarmonic = 1;
        maxHarmonic = 4;
        break;
      case -1: /* Use highest known PN order - move if new orders added */
      case 3:
        minHarmonic = 1;
        maxHarmonic = 5;
        break;
      default:
        XLALPrintError("XLAL Error - %s: Invalid amp. PN order %d\n",
        __func__, amplitudeO );
        XLAL_ERROR(XLAL_EINVAL);
        break;
    }


    XLAL_BEGINGSL;

    // The TaylorT2 version sometimes returns duplicate time steps (i.e. two different steps happen at the same time, because df/dt is too large near the end)
    // So if this happens, we have to remove the duplicates and try again. If it fails then, there is another problem and the function exits with an error.
    UINT4 weHadAProblem = 0;

    accel = gsl_interp_accel_alloc();
    accel_k = gsl_interp_accel_alloc();
    interp_t_of_f = gsl_spline_alloc(gsl_interp_cspline,
                                     length);
    status = gsl_spline_init(interp_t_of_f, &y->data[length*2],
                             y->data, length);
    if(status != GSL_SUCCESS)
    {
      weHadAProblem = 1;
    }

    // initiate interpolation objects
    for(i = 1; i < nParams; i++) {
      interp[i-1] = gsl_spline_alloc(gsl_interp_cspline,
                                     length);
      status = gsl_spline_init(interp[i-1], y->data,
                               &y->data[length*i], length);
      if(status != GSL_SUCCESS)
      {
        weHadAProblem = 1;
      }
    }

    if(weHadAProblem)
    {
      y = removeDuplicates(y);
      length = y->dimLength->data[1];

      gsl_spline_free(interp_t_of_f);

      for(i = 0; i < nParams-1; i++) {
        gsl_spline_free(interp[i]);
      }

      interp_t_of_f = gsl_spline_alloc(gsl_interp_cspline,
                                       length);
      status = gsl_spline_init(interp_t_of_f, &y->data[length*2],
                               y->data, length);
      if(status != GSL_SUCCESS)
      {
        XLALPrintError("XLAL Error - %s: spline interpolation returned with error code %d\n",
                       __func__, status );
        XLAL_ERROR(XLAL_EFUNC);
      }

      // initiate interpolation objects
      for(i = 1; i < nParams; i++) {
        interp[i-1] = gsl_spline_alloc(gsl_interp_cspline,
                                       length);
        status = gsl_spline_init(interp[i-1], y->data,
                                 &y->data[length*i], length);
        if(status != GSL_SUCCESS)
        {
          XLALPrintError("XLAL Error - %s: spline interpolation returned with error code %d\n",
                         __func__, status );
          XLAL_ERROR(XLAL_EFUNC);
        }
      }
    }

    // compute tISCO
    status = gsl_spline_eval_e(interp_t_of_f, omegaHatISCO,
                               accel, &tISCO);
    if(status != GSL_SUCCESS)
    {
      // if ISCO not reached, set t=0 at the end of tthe simulation
      tISCO = y->data[length-1];
    }

    // set tStart so that tISCO corresponds to t=0.
    LIGOTimeGPS tStart = LIGOTIMEGPSZERO;
    REAL8 omegaHatMin = LAL_PI*Msec*fMin;
    REAL8 tMin;
    status = gsl_spline_eval_e(interp_t_of_f, omegaHatMin,
                               accel, &tMin);
    tMin -= tISCO;
    tMin *= Msec;
    tStart.gpsSeconds += (INT4)floor(tMin);
    tStart.gpsNanoSeconds += (INT4)floor(1.e9*(tMin - tStart.gpsSeconds));

    // number of samples
    UINT4 nFreqSamples, iFirstSample;
    fMin = deltaF*floor(nextafter(fMin, INFINITY)/deltaF);
    if(fMax > fMin)
    {
      nFreqSamples = (UINT4)ceil((fMax - fMin)/deltaF);
    }
    else
    {
      nFreqSamples = (UINT4)ceil((fISCO - fMin)/deltaF);
    }
    iFirstSample = fMin/deltaF;

    *hplus = XLALCreateCOMPLEX16FrequencySeries("PLUS POLARIZATION", &tStart,
    0., deltaF, &lalSecondUnit, iFirstSample+nFreqSamples);
    *hcross = XLALCreateCOMPLEX16FrequencySeries("CROSS POLARIZATION", &tStart,
    0., deltaF, &lalSecondUnit, iFirstSample+nFreqSamples);
    XLALUnitMultiply(&(*hplus)->sampleUnits, &(*hplus)->sampleUnits, &lalStrainUnit);
    XLALUnitMultiply(&(*hcross)->sampleUnits, &(*hcross)->sampleUnits, &lalStrainUnit);

    COMPLEX16FrequencySeries* hplustilde = *hplus;
    COMPLEX16FrequencySeries* hcrosstilde = *hcross;

    memset(hplustilde->data->data, 0, sizeof(COMPLEX16)*(iFirstSample+nFreqSamples));
    memset(hcrosstilde->data->data, 0, sizeof(COMPLEX16)*(iFirstSample+nFreqSamples));

    //     If nFreqSamples > length, better interpolate A(t), t(f), and Psi(f)
    //     If nFreqSamples < length, better interpolate everything
    if(nFreqSamples < length)
    {
      // We interpolate every quantity at each frequency sample of the output

      for(i = 0, j = iFirstSample; i < nFreqSamples; i++, j++) // loop over frequency samples
      {
        fCur = j*deltaF;

        for(nHarmonic = minHarmonic; nHarmonic <= maxHarmonic; nHarmonic++) // loop over relevant harmonics
        {
          htPlusHarmonic = 0.;
          htCrossHarmonic = 0.;

          omegaHatCur = freqToOmegaHat*fCur/nHarmonic;

          status = gsl_spline_eval_e(interp_t_of_f, omegaHatCur,
          accel, &tCur);
          // if return is non-zero, then the harmonic contributes nothing
          if(status == GSL_SUCCESS)
          {
            status = gsl_spline_eval_e(interp[14], tCur, accel,
            &dydt);
            status = gsl_spline_eval_e(interp[0], tCur, accel,
            &phiOrb);

            T = 1./sqrt(fabs(nHarmonic*dydt));

            // k = 0
            iK = iMinKMax[kMaxCur];

            status = gsl_spline_eval_e(interp[1], tCur, accel_k,
            &omegaHatCur);
            if(status != GSL_SUCCESS)
            {
              continue;
            }
            status = gsl_spline_eval_e(interp[2], tCur, accel_k,
            &lnhatxCur);
            status = gsl_spline_eval_e(interp[3], tCur, accel_k,
            &lnhatyCur);
            status = gsl_spline_eval_e(interp[4], tCur, accel_k,
            &lnhatzCur);
            // if the amplitude PN order is too low, no need to interpolate the spin orientation
            if(amplitudeO > 1 || amplitudeO == -1)
            {
              status = gsl_spline_eval_e(interp[5], tCur, accel_k,
              &s1xCur);
              status = gsl_spline_eval_e(interp[6], tCur, accel_k,
              &s1yCur);
              status = gsl_spline_eval_e(interp[7], tCur, accel_k,
              &s1zCur);
              status = gsl_spline_eval_e(interp[8], tCur, accel_k,
              &s2xCur);
              status = gsl_spline_eval_e(interp[9], tCur, accel_k,
              &s2yCur);
              status = gsl_spline_eval_e(interp[10], tCur,
              accel_k, &s2zCur);
            }
            status = gsl_spline_eval_e(interp[11], tCur, accel_k,
            &E1xCur);
            status = gsl_spline_eval_e(interp[12], tCur, accel_k,
            &E1yCur);
            status = gsl_spline_eval_e(interp[13], tCur, accel_k,
            &E1zCur);

            vCur = cbrt(omegaHatCur);

            XLALSimInspiralPrecessingPolarizationWaveformHarmonic(&htPlusTemp,
            &htCrossTemp, vCur, s1xCur, s1yCur, s1zCur, s2xCur, s2yCur, s2zCur,
            lnhatxCur, lnhatyCur, lnhatzCur, E1xCur, E1yCur, E1zCur, dm, eta,
            v0, nHarmonic, amplitudeO);

            htPlusHarmonic += akCsts[iK]*T*htPlusTemp;
            htCrossHarmonic += akCsts[iK]*T*htCrossTemp;

            for(k = 1, iK++; k <= kMaxCur; k++, iK++) // loop over k != 0
            {
              tK = tCur + k*T;

              status = gsl_spline_eval_e(interp[1], tK, accel_k,
              &omegaHatCur);
              // if return is non-zero, then the k value contributes nothing
              if(status == GSL_SUCCESS)
              {
                status = gsl_spline_eval_e(interp[2], tK, accel_k,
                &lnhatxCur);
                status = gsl_spline_eval_e(interp[3], tK, accel_k,
                &lnhatyCur);
                status = gsl_spline_eval_e(interp[4], tK, accel_k,
                &lnhatzCur);
                // if the amplitude PN order is too low, no need to interpolate the spin orientation
                if(amplitudeO > 1 || amplitudeO == -1)
                {
                  status = gsl_spline_eval_e(interp[5], tK,
                  accel_k, &s1xCur);
                  status = gsl_spline_eval_e(interp[6], tK,
                  accel_k, &s1yCur);
                  status = gsl_spline_eval_e(interp[7], tK,
                  accel_k, &s1zCur);
                  status = gsl_spline_eval_e(interp[8], tK,
                  accel_k, &s2xCur);
                  status = gsl_spline_eval_e(interp[9], tK,
                  accel_k, &s2yCur);
                  status = gsl_spline_eval_e(interp[10], tK,
                  accel_k, &s2zCur);
                }
                status = gsl_spline_eval_e(interp[11], tK,
                accel_k, &E1xCur);
                status = gsl_spline_eval_e(interp[12], tK,
                accel_k, &E1yCur);
                status = gsl_spline_eval_e(interp[13], tK,
                accel_k, &E1zCur);
                status = gsl_spline_eval_e(interp[14], tK,
                accel_k, &dydt);

                vCur = cbrt(omegaHatCur);
                TCur = 1./sqrt(fabs(nHarmonic*dydt));

                XLALSimInspiralPrecessingPolarizationWaveformHarmonic(
                &htPlusTemp, &htCrossTemp, vCur, s1xCur, s1yCur, s1zCur,
                s2xCur, s2yCur, s2zCur, lnhatxCur, lnhatyCur, lnhatzCur, E1xCur,
                E1yCur, E1zCur, dm, eta, v0, nHarmonic, amplitudeO);

                htPlusHarmonic += akCsts[iK]*TCur*htPlusTemp;
                htCrossHarmonic += akCsts[iK]*TCur*htCrossTemp;
              }

              tK = tCur - k*T;

              status = gsl_spline_eval_e(interp[1], tK, accel_k,
              &omegaHatCur);
              // if return is non-zero, then the -k value contributes nothing
              if(status == GSL_SUCCESS)
              {
                status = gsl_spline_eval_e(interp[2], tK,
                accel_k, &lnhatxCur);
                status = gsl_spline_eval_e(interp[3], tK,
                accel_k, &lnhatyCur);
                status = gsl_spline_eval_e(interp[4], tK,
                accel_k, &lnhatzCur);
                // if the amplitude PN order is too low, no need to interpolate the spin orientation
                if(amplitudeO > 1 || amplitudeO == -1)
                {
                  status = gsl_spline_eval_e(interp[5], tK,
                  accel_k, &s1xCur);
                  status = gsl_spline_eval_e(interp[6], tK,
                  accel_k, &s1yCur);
                  status = gsl_spline_eval_e(interp[7], tK,
                  accel_k, &s1zCur);
                  status = gsl_spline_eval_e(interp[8], tK,
                  accel_k, &s2xCur);
                  status = gsl_spline_eval_e(interp[9], tK,
                  accel_k, &s2yCur);
                  status = gsl_spline_eval_e(interp[10], tK,
                  accel_k, &s2zCur);
                }
                status = gsl_spline_eval_e(interp[11], tK,
                accel_k, &E1xCur);
                status = gsl_spline_eval_e(interp[12], tK,
                accel_k, &E1yCur);
                status = gsl_spline_eval_e(interp[13], tK,
                accel_k, &E1zCur);
                status = gsl_spline_eval_e(interp[14], tK,
                accel_k, &dydt);

                vCur = cbrt(omegaHatCur);
                TCur = 1./sqrt(fabs(nHarmonic*dydt));

                XLALSimInspiralPrecessingPolarizationWaveformHarmonic(
                &htPlusTemp, &htCrossTemp, vCur, s1xCur, s1yCur, s1zCur, s2xCur,
                s2yCur, s2zCur, lnhatxCur, lnhatyCur, lnhatzCur, E1xCur, E1yCur,
                E1zCur, dm, eta, v0, nHarmonic, amplitudeO);

                htPlusHarmonic += akCsts[iK]*TCur*htPlusTemp;
                htCrossHarmonic += akCsts[iK]*TCur*htCrossTemp;
              }
            }

            // apply phase factor for the particular harmonic
            cPhase = prefac*cexp(I*(freqToOmegaHat*fCur*(tCur - tISCO) - nHarmonic*phiOrb
            - LAL_PI_4));
            hplustilde->data->data[j] += conj(htPlusHarmonic*cPhase);
            hcrosstilde->data->data[j] += conj(htCrossHarmonic*cPhase);
          }
        }
      }
    }
    else
    {
      // if(nFreqSamples > length)
      // We compute the amplitude and phase of each harmonic for each sample of the Runge-Kutta, and then interpolate between them

      REAL8Array* freq = XLALCreateREAL8ArrayL(2, 5, length);
      REAL8Array* hPlusAmpRe = XLALCreateREAL8ArrayL(2, 5, length);
      REAL8Array* hPlusAmpIm = XLALCreateREAL8ArrayL(2, 5, length);
      REAL8Array* hCrossAmpRe = XLALCreateREAL8ArrayL(2, 5, length);
      REAL8Array* hCrossAmpIm = XLALCreateREAL8ArrayL(2, 5, length);
      REAL8Array* TArray = XLALCreateREAL8ArrayL(2, 5, length);
      REAL8Array* phase = XLALCreateREAL8ArrayL(2, 5, length);

      gsl_spline *interpHPlusRe[5];
      gsl_spline *interpHPlusIm[5];
      gsl_spline *interpHCrossRe[5];
      gsl_spline *interpHCrossIm[5];
      gsl_spline *interpPhase[5];
      gsl_spline *interpT[5];

      for(i = 0; i < length; i++) // loop over the samples output by the Runge-Kutta
      {
        tCur = y->data[i];
        omegaHatCur = y->data[2*length+i];
        vCur = cbrt(omegaHatCur);
        phiOrb = y->data[length+i];
        dydt = y->data[15*length+i];

        // create the series for the frequency, the phase, and the time derivative of the orbital frequency
        for(nHarmonic = minHarmonic; nHarmonic <= maxHarmonic; nHarmonic++)
        {
          freq->data[length*(nHarmonic-1)+i] = nHarmonic*omegaHatCur/
          freqToOmegaHat;
          phase->data[length*(nHarmonic-1)+i] = nHarmonic*(omegaHatCur*(tCur - tISCO)
          - phiOrb) - LAL_PI_4;
          TArray->data[length*(nHarmonic-1)+i] = 1./sqrt(fabs(nHarmonic*dydt));
        }

        for(nHarmonic = minHarmonic; nHarmonic <= maxHarmonic; nHarmonic++) // loop over the harmonic number
        {
          lnhatxCur = y->data[3*length+i];
          lnhatyCur = y->data[4*length+i];
          lnhatzCur = y->data[5*length+i];
          // if the amplitude PN order is too low, no need to take the spin orientation into account
          if(amplitudeO > 1 || amplitudeO == -1)
          {
            s1xCur = y->data[6*length+i];
            s1yCur = y->data[7*length+i];
            s1zCur = y->data[8*length+i];
            s2xCur = y->data[9*length+i];
            s2yCur = y->data[10*length+i];
            s2zCur = y->data[11*length+i];
          }
          E1xCur = y->data[12*length+i];
          E1yCur = y->data[13*length+i];
          E1zCur = y->data[14*length+i];

          XLALSimInspiralPrecessingPolarizationWaveformHarmonic(&htPlusTemp,
          &htCrossTemp, vCur, s1xCur, s1yCur, s1zCur, s2xCur, s2yCur, s2zCur,
          lnhatxCur, lnhatyCur, lnhatzCur, E1xCur, E1yCur, E1zCur, dm, eta, v0,
          nHarmonic, amplitudeO);

          htPlusHarmonic = (prefac*TArray->data[length*(nHarmonic-1)+i])*
          htPlusTemp;
          htCrossHarmonic = (prefac*TArray->data[length*(nHarmonic-1)+i])*
          htCrossTemp;

          // create the series for the amplitude
          hPlusAmpRe->data[length*(nHarmonic-1)+i] = creal(htPlusHarmonic);
          hPlusAmpIm->data[length*(nHarmonic-1)+i] = cimag(htPlusHarmonic);
          hCrossAmpRe->data[length*(nHarmonic-1)+i] = creal(htCrossHarmonic);
          hCrossAmpIm->data[length*(nHarmonic-1)+i] = cimag(htCrossHarmonic);
        }
      }

      //initiate the interpolation structures
      for(nHarmonic = minHarmonic; nHarmonic <= maxHarmonic; nHarmonic++)
      {
        interpHPlusRe[nHarmonic-1] =
        gsl_spline_alloc(gsl_interp_cspline, length);
        status = gsl_spline_init(interpHPlusRe[nHarmonic-1],
        y->data, &hPlusAmpRe->data[length*(nHarmonic-1)], length);
        if(status != GSL_SUCCESS)
        {
          XLALPrintError("XLAL Error - %s: spline interpolation returned with error code %d\n",
          __func__, status );
          XLAL_ERROR(XLAL_EFUNC);
        }

        interpHPlusIm[nHarmonic-1] =
        gsl_spline_alloc(gsl_interp_cspline, length);
        status = gsl_spline_init(interpHPlusIm[nHarmonic-1],
        y->data, &hPlusAmpIm->data[length*(nHarmonic-1)], length);
        if(status != GSL_SUCCESS)
        {
          XLALPrintError("XLAL Error - %s: spline interpolation returned with error code %d\n",
          __func__, status );
          XLAL_ERROR(XLAL_EFUNC);
        }

        interpHCrossRe[nHarmonic-1] =
        gsl_spline_alloc(gsl_interp_cspline, length);
        status = gsl_spline_init(interpHCrossRe[nHarmonic-1],
        y->data, &hCrossAmpRe->data[length*(nHarmonic-1)], length);
        if(status != GSL_SUCCESS)
        {
          XLALPrintError("XLAL Error - %s: spline interpolation returned with error code %d\n",
          __func__, status );
          XLAL_ERROR(XLAL_EFUNC);
        }

        interpHCrossIm[nHarmonic-1] =
        gsl_spline_alloc(gsl_interp_cspline, length);
        status = gsl_spline_init(interpHCrossIm[nHarmonic-1],
        y->data, &hCrossAmpIm->data[length*(nHarmonic-1)], length);
        if(status != GSL_SUCCESS)
        {
          XLALPrintError("XLAL Error - %s: spline interpolation returned with error code %d\n",
          __func__, status );
          XLAL_ERROR(XLAL_EFUNC);
        }

        interpPhase[nHarmonic-1] =
        gsl_spline_alloc(gsl_interp_cspline, length);
        status = gsl_spline_init(interpPhase[nHarmonic-1],
        &freq->data[length*(nHarmonic-1)], &phase->data[length*(nHarmonic-1)],
        length);
        if(status != GSL_SUCCESS)
        {
          XLALPrintError("XLAL Error - %s: spline interpolation returned with error code %d\n",
          __func__, status );
          XLAL_ERROR(XLAL_EFUNC);
        }

        interpT[nHarmonic-1] =
        gsl_spline_alloc(gsl_interp_cspline, length);
        status = gsl_spline_init(interpT[nHarmonic-1],
        &freq->data[length*(nHarmonic-1)], &TArray->data[length*(nHarmonic-1)],
        length);
        if(status != GSL_SUCCESS)
        {
          XLALPrintError("XLAL Error - %s: spline interpolation returned with error code %d\n",
          __func__, status );
          XLAL_ERROR(XLAL_EFUNC);
        }
      }

      REAL8 Psi = 0.;
      REAL8 htPlusRe = 0.;
      REAL8 htPlusIm = 0.;
      REAL8 htCrossRe = 0.;
      REAL8 htCrossIm = 0.;

      for(i = 0, j = iFirstSample; i < nFreqSamples; i++, j++) // loop over the output frequencies
      {
        fCur = j*deltaF;

        for(nHarmonic = minHarmonic; nHarmonic <= maxHarmonic; nHarmonic++) // loop over the harmonic number
        {
          omegaHatCur = freqToOmegaHat*fCur/nHarmonic;

          status = gsl_spline_eval_e(interp_t_of_f, omegaHatCur,
          accel, &tCur);
          // if return is non-zero, then the harmonic contributes nothing
          if(status == GSL_SUCCESS)
          {
            status = gsl_spline_eval_e(interpPhase[nHarmonic-1],
            fCur, accel, &Psi);
            status = gsl_spline_eval_e(interpT[nHarmonic-1],
            fCur, accel, &T);

            // k = 0
            iK = iMinKMax[kMaxCur];

            status = gsl_spline_eval_e(interpHPlusRe[nHarmonic-1], tCur,
            accel_k, &htPlusRe);
            status = gsl_spline_eval_e(interpHPlusIm[nHarmonic-1], tCur,
            accel_k, &htPlusIm);
            status = gsl_spline_eval_e(interpHCrossRe[nHarmonic-1], tCur,
            accel_k, &htCrossRe);
            status = gsl_spline_eval_e(interpHCrossIm[nHarmonic-1], tCur,
            accel_k, &htCrossIm);

            htPlusHarmonic = akCsts[iK]*(htPlusRe + I*htPlusIm);
            htCrossHarmonic = akCsts[iK]*(htCrossRe + I*htCrossIm);

            for(k = 1, iK++; k <= kMaxCur; k++, iK++) // loop over k
            {
              tK = tCur + k*T;

              status =
              gsl_spline_eval_e(interpHPlusRe[nHarmonic-1], tK,
              accel_k, &htPlusRe);
              // if return is non-zero, then the k value contributes nothing
              if(status == GSL_SUCCESS)
              {
                status =
                gsl_spline_eval_e(interpHPlusIm[nHarmonic-1], tK,
                accel_k, &htPlusIm);
                status =
                gsl_spline_eval_e(interpHCrossRe[nHarmonic-1], tK, accel_k,
                &htCrossRe);
                status =
                gsl_spline_eval_e(interpHCrossIm[nHarmonic-1], tK, accel_k,
                &htCrossIm);

                htPlusHarmonic += akCsts[iK]*(htPlusRe + I*htPlusIm);
                htCrossHarmonic += akCsts[iK]*(htCrossRe + I*htCrossIm);
              }

              tK = tCur - k*T;

              status = gsl_spline_eval_e(interpHPlusRe[nHarmonic-1], tK,
                accel_k, &htPlusRe);
              // if return is non-zero, then the -k value contributes nothing
              if(status == GSL_SUCCESS)
              {
                status = gsl_spline_eval_e(interpHPlusIm[nHarmonic-1], tK,
                accel_k, &htPlusIm);
                status = gsl_spline_eval_e(interpHCrossRe[nHarmonic-1], tK,
                accel_k, &htCrossRe);
                status = gsl_spline_eval_e(interpHCrossIm[nHarmonic-1], tK,
                accel_k, &htCrossIm);

                htPlusHarmonic += akCsts[iK]*(htPlusRe + I*htPlusIm);
                htCrossHarmonic += akCsts[iK]*(htCrossRe + I*htCrossIm);
              }
            }

            cPhase = cexp(I*Psi);

            hplustilde->data->data[j] += conj(htPlusHarmonic*cPhase);
            hcrosstilde->data->data[j] += conj(htCrossHarmonic*cPhase);
          }
        }
      }

      for(nHarmonic = minHarmonic; nHarmonic <= maxHarmonic; nHarmonic++) // free the interpolatin objects
      {
        gsl_spline_free(interpHPlusRe[nHarmonic-1]);
        gsl_spline_free(interpHPlusIm[nHarmonic-1]);
        gsl_spline_free(interpHCrossRe[nHarmonic-1]);
        gsl_spline_free(interpHCrossIm[nHarmonic-1]);
        gsl_spline_free(interpT[nHarmonic-1]);
        gsl_spline_free(interpPhase[nHarmonic-1]);
      }

      XLALDestroyREAL8Array(freq);
      XLALDestroyREAL8Array(hPlusAmpRe);
      XLALDestroyREAL8Array(hPlusAmpIm);
      XLALDestroyREAL8Array(hCrossAmpRe);
      XLALDestroyREAL8Array(hCrossAmpIm);
      XLALDestroyREAL8Array(TArray);
      XLALDestroyREAL8Array(phase);
    }



    gsl_interp_accel_free(accel);
    gsl_interp_accel_free(accel_k);
    gsl_spline_free(interp_t_of_f);

    for(i = 0; i < nParams-1; i++) {
      gsl_spline_free(interp[i]);
    }

    XLAL_ENDGSL;


    XLALDestroyREAL8Array(y);

    return XLAL_SUCCESS;
}

/*
 * Internal function called by the integration routine.
 * Stops the integration if
 * 1) The energy decreases with increasing orbital frequency
 * 2) The orbital frequency begins decreasing
 * 3) The orbital frequency becomes infinite
 * 4) The orbital frequency has gone outside the requested bounds
 * 5) The PN parameter v/c becomes >= 1
 * SpinTaylorT4 and SpinTaylorT2 both use this same stopping test
 */
static int XLALSimInspiralSpinTaylorStoppingTest(
	double UNUSED t,
	const double values[],
	double dvalues[], 
	void *mparams
	)
{//   fprintf(stdout,"Enter %s\n","StopingTest");
    REAL8 omega, v, test, omegaStart, omegaEnd, ddomega;
    REAL8 LNhx, LNhy, LNhz, S1x, S1y, S1z, S2x, S2y, S2z;
    REAL8 LNdotS1, LNdotS2, S1dotS2, S1sq, S2sq;
    XLALSimInspiralSpinTaylorTxCoeffs *params 
            = (XLALSimInspiralSpinTaylorTxCoeffs*) mparams;
    /* Spin-corrections to energy (including dynamical terms) */
    REAL8 Espin15 = 0., Espin2 = 0., Espin25 = 0., Espin35 = 0.;

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

    //  fprintf(stderr,"omegaStart %f\n",omegaStart);
    //  fprintf(stderr,"omegaend %f\n",omegaEnd);


    switch( params->spinO )
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
            // Compute 3.5PN SO correction to energy
            // See Eq. 3.15 of arXiv:1303.7412
            // Note that S_l/M^2 = (m1/M)^2 chi1 + (m2/M)^2 chi2
            // and Sigma_l/M^2 = (m2/M) chi2 - (m1/M) chi1
            Espin35 += params->E7S1O * LNdotS1 + params->E7S2O * LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
            // Compute 2.5PN SO correction to energy
            // See Eq. 7.9 of gr-qc/0605140v4
            // Note that S_l/M^2 = (m1/M)^2 chi1 + (m2/M)^2 chi2
            // and Sigma_l/M^2 = (m2/M) chi2 - (m1/M) chi1
            Espin25 += params->E5S1O * LNdotS1 + params->E5S2O * LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            // Compute S1-S2 spin-spin term
            S1dotS2 = (S1x*S2x + S1y*S2y + S1z*S2z);
            Espin2 += params->E4S1S2  * S1dotS2 + params->E4S1OS2O * LNdotS1 * LNdotS2;
            // Compute 2PN quadrupole-monopole correction to energy
            // See last line of Eq. 6 of astro-ph/0504538
            // or 2nd and 3rd lines of Eq. (C4) in arXiv:0810.5336v3
            S1sq = (S1x*S1x + S1y*S1y + S1z*S1z);
            S2sq = (S2x*S2x + S2y*S2y + S2z*S2z);
            Espin2 += params->E4QMS1S1 * S1sq
                    + params->E4QMS2S2 * S2sq
                    + params->E4QMS1OS1O * LNdotS1 * LNdotS1
                    + params->E4QMS2OS2O * LNdotS2 * LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            // Compute 1.5PN SO correction to energy
            Espin15 += params->E3S1O * LNdotS1 + params->E3S2O * LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %d\n",
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
  //  fprintf(stderr,"TestValue %f\n",test);
 //fprintf(stderr,"test",test)
    // Check d^2omega/dt^2 > 0 (true iff current_domega - prev_domega > 0)
    ddomega = dvalues[1] - params->prev_domega;
    // ...but we can integrate forward or backward, so beware of the sign!
    if ( params->fEnd < params->fStart && params->fEnd != 0.
                && params->prev_domega != 0.)
        ddomega *= -1;
    // Copy current value of domega to prev. value of domega for next call
    params->prev_domega = dvalues[1];

    if( fabs(omegaEnd) > LAL_REAL4_EPS && omegaEnd > omegaStart
                && omega > omegaEnd) /* freq. above bound */
        //fprintf(stderr,"one %s\n");
        return LALSIMINSPIRAL_ST_TEST_FREQBOUND;
    else if( fabs(omegaEnd) > LAL_REAL4_EPS && omegaEnd < omegaStart
                && omega < omegaEnd) /* freq. below bound */
        //fprintf(stderr,"two %s\n");
        return LALSIMINSPIRAL_ST_TEST_FREQBOUND;
    else if (test < 0.0) /* energy test fails! */
        //fprintf(stderr,"three %s\n");
        return LALSIMINSPIRAL_ST_TEST_ENERGY;
    else if isnan(omega) /* omega is nan! */
        //fprintf(stderr,"four %s\n");
        return LALSIMINSPIRAL_ST_TEST_OMEGANAN;
    else if (v >= 1.) // v/c >= 1!
        //fprintf(stderr,"five %s\n");
        return LALSIMINSPIRAL_ST_TEST_LARGEV;
    else if (ddomega <= 0.) // d^2omega/dt^2 <= 0!
        //fprintf(stderr,"six %s\n");
        return LALSIMINSPIRAL_ST_TEST_OMEGADOUBLEDOT;
    else /* Step successful, continue integrating */
        //fprintf(stderr,"SUCCESS %s\n");
        return GSL_SUCCESS;
}

/*
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
	double UNUSED t,
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
    REAL8 v, v2, v11;
    REAL8 LNhdotS1, LNhdotS2, S1dotS2, S1sq, S2sq;
    REAL8 wspin15 = 0., wspin2 = 0., wspin25 = 0., wspin3 = 0., wspin35 = 0.;

    XLALSimInspiralSpinTaylorTxCoeffs *params 
            = (XLALSimInspiralSpinTaylorTxCoeffs*) mparams;

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
    v2  = v * v;
    v11= omega*omega*omega*v2;

    LNhdotS1 = (LNhx*S1x + LNhy*S1y + LNhz*S1z);
    LNhdotS2 = (LNhx*S2x + LNhy*S2y + LNhz*S2z);
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
            wspin35 = params->wdot7S1O*LNhdotS1 + params->wdot7S2O*LNhdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
            // Compute 3PN SO correction to domega/dt
            // See Eq. 13 of arXiv:1210.0764 and/or Eq. 3.16 of arXiv:1303.7412
            // Note that S_l/M^2 = (m1/M)^2 chi1 + (m2/M)^2 chi2
            // and Sigma_l/M^2 = (m2/M) chi2 - (m1/M) chi1
            wspin3 = params->wdot6S1O * LNhdotS1 + params->wdot6S2O * LNhdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
            // Compute 2.5PN SO correction to domega/dt
            // See Eq. 8.3 of gr-qc/0605140v4
            // Note that S_l/M^2 = (m1/M)^2 chi1 + (m2/M)^2 chi2
            // and Sigma_l/M^2 = (m2/M) chi2 - (m1/M) chi1
            wspin25 = params->wdot5S1O*LNhdotS1 + params->wdot5S2O*LNhdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            // Compute S1-S2 spin-spin term
            S1dotS2 = (S1x*S2x + S1y*S2y + S1z*S2z);
            wspin2 = params->wdot4S1S2 *S1dotS2 + params->wdot4S1OS2O * LNhdotS1 * LNhdotS2;
            // Compute 2PN QM and self-spin corrections to domega/dt
            // This is equivalent to Eqs. 9c + 9d of astro-ph/0504538
            S1sq = (S1x*S1x + S1y*S1y + S1z*S1z);
            S2sq = (S2x*S2x + S2y*S2y + S2z*S2z);
            wspin2 += params->wdot4QMS1S1 * S1sq
                    + params->wdot4QMS2S2 * S2sq
                    + params->wdot4QMS1OS1O * LNhdotS1 * LNhdotS1
                    + params->wdot4QMS2OS2O * LNhdotS2 * LNhdotS2
                    + params->wdot4S1S1 * S1sq
                    + params->wdot4S2S2 * S2sq
                    + params->wdot4S1OS1O * LNhdotS1 * LNhdotS1
                    + params->wdot4S2OS2O * LNhdotS2 * LNhdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            // Compute 1.5PN SO correction to domega/dt
            wspin15 = params->wdot3S1O*LNhdotS1 + params->wdot3S2O*LNhdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %d\n",
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
                    + params->wdotlogcoeff * log(v)
            + v * ( params->wdotcoeff[7] + wspin35
            + omega * ( params->wdottidal5pn
            + v2 * ( params->wdottidal6pn ) ) ) ) ) ) ) ) ) );
     // fprintf(stdout,"domega\n%f",domega);

    XLALSimInspiralSpinDerivatives(&dLNhx,&dLNhy,&dLNhz,&dE1x,&dE1y,&dE1z,&dS1x,&dS1y,&dS1z,&dS2x,&dS2y,&dS2z,v,omega,LNhx,LNhy,LNhz,E1x,E1y,E1z,S1x,S1y,S1z,S2x,S2y,S2z,LNhdotS1,LNhdotS2,params);

    /* dphi = d \phi / d \hat{t} = M d \phi /dt = M \omega = \hat{\omega} */
    ds = omega;

    dvalues[0]    = ds   ; dvalues[1]     = domega;
    dvalues[2]    = dLNhx; dvalues[3]     = dLNhy ; dvalues[4]    = dLNhz;
    dvalues[5]    = dS1x ; dvalues[6]     = dS1y  ; dvalues[7]    = dS1z ;
    dvalues[8]    = dS2x ; dvalues[9]     = dS2y  ; dvalues[10]   = dS2z ;
    dvalues[11]   = dE1x ; dvalues[12]    = dE1y  ; dvalues[13]   = dE1z ;

    return GSL_SUCCESS;
}

static int XLALSimInspiralSpinTaylorT1Derivatives(
	double UNUSED t,
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
    REAL8 v, v2, v3, v4, v7, v11;
    REAL8 LNdotS1, LNdotS2, S1dotS2, S1sq, S2sq;
    //REAL8 wspin15 = 0., wspin2 = 0., wspin25 = 0., wspin3 = 0., wspin35 = 0.;
    REAL8 Fspin15 = 0., Fspin2 = 0., Fspin25 = 0., Fspin3 = 0., Fspin35 = 0.;
    REAL8 Espin15 = 0., Espin2 = 0., Espin25 = 0., Espin35 = 0.;

    XLALSimInspiralSpinTaylorTxCoeffs *params
            = (XLALSimInspiralSpinTaylorTxCoeffs*) mparams;

     //UNUSED(t);

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
    v7 = v4 * v3; v11 = v7 * v4;

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
           Fspin35 += params->F7S1O*LNdotS1 + params->F7S1O*LNdotS2;
	   Espin35 += params->E7S1O*LNdotS1 + params->E7S2O*LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
            // Compute 3PN SO correction to domega/dt
            // See Eq. 13 of arXiv:1210.0764 and/or Eq. 3.16 of arXiv:1303.7412
            // Note that S_l/M^2 = (m1/M)^2 chi1 + (m2/M)^2 chi2
            // and Sigma_l/M^2 = (m2/M) chi2 - (m1/M) chi1
            Fspin3 += params->F6S1O* LNdotS1 + params->F6S2O* LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
            // Compute 2.5PN SO correction to domega/dt
            // See Eq. 8.3 of gr-qc/0605140v4
            // Note that S_l/M^2 = (m1/M)^2 chi1 + (m2/M)^2 chi2
            // and Sigma_l/M^2 = (m2/M) chi2 - (m1/M) chi1
            Fspin25 += params->F5S1O*LNdotS1 + params->F5S2O*LNdotS2;
	    Espin25 += params->E5S1O*LNdotS1 + params->E5S2O*LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            // Compute S1-S2 spin-spin term
            S1dotS2 = (S1x*S2x + S1y*S2y + S1z*S2z);
            Fspin2 += params->F4S1S2*S1dotS2 + params->F4S1OS2O*LNdotS1 * LNdotS2;
            // Compute 2PN QM and self-spin corrections to domega/dt
            // This is equivalent to Eqs. 9c + 9d of astro-ph/0504538
            S1dotS2 = (S1x*S2x + S1y*S2y + S1z*S2z);
            Espin2 += params->E4S1S2  * S1dotS2 + params->E4S1OS2O * LNdotS1 * LNdotS2;
	    S1sq = (S1x*S1x + S1y*S1y + S1z*S1z);
            S2sq = (S2x*S2x + S2y*S2y + S2z*S2z);
            Espin2 += params->E4QMS1S1 * S1sq
                    + params->E4QMS2S2 * S2sq
                    + params->E4QMS1OS1O * LNdotS1 * LNdotS1
                    + params->E4QMS2OS2O * LNdotS2 * LNdotS2;
	    Fspin2 += params->F4QMS1S1 * S1sq
                    + params->F4QMS2S2 * S2sq
                    + params->F4QMS1OS1O * LNdotS1 * LNdotS1
                    + params->F4QMS2OS2O * LNdotS2 * LNdotS2
                    + params->F4S1S1 * S1sq
                    + params->F4S2S2 * S2sq
                    + params->F4S1OS1O * LNdotS1 * LNdotS1
                    + params->F4S2OS2O * LNdotS2 * LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            // Compute 1.5PN SO correction to domega/dt
            Fspin15 = params->F3S1O*LNdotS1 + params->F3S2O*LNdotS2;
	   Espin15 += params->E3S1O * LNdotS1 + params->E3S2O * LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %d\n",
                    __func__, params->spinO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    domega  =-3.0* (params->Fnewt/params->dEdvnewt) * v11 *(((( 1. + v * v * (  params->Fcoeff[2]
            + v * (  (params->Fcoeff[3]+ Fspin15 )
            + v * (  (params->Fcoeff[4] + Fspin2 )
            + v * (  (params->Fcoeff[5] + Fspin25)
            + v * (  (params->Fcoeff[6] + (params->Flogcoeff*log(v)) + Fspin3 )
            + v * (  (params->Fcoeff[7] + Fspin35)+ v * v * v * (  params->Ftidal5pn
                        + v * v * ( params->Ftidal6pn ) )
                         )))))))))/ ( 2. + v * v * ( 4. * params->Ecoeff[2]
            + v * (( 5. * (params->Ecoeff[3] +  Espin15 ))
            + v * (( 6. * (params->Ecoeff[4] +  Espin2 ))
            + v * ((7. * (params->Ecoeff[5] +  Espin25 ))
            + v * (( 8. * (params->Ecoeff[6] ))
            + v * ( (9. * (params->Ecoeff[7] + Espin35 ))+ v * v * v * (( 12. * params->Etidal5pn)
                        + v * v * ( 14. * params->Etidal6pn ) )
                         ))))))));

    XLALSimInspiralSpinDerivatives(&dLNhx,&dLNhy,&dLNhz,&dE1x,&dE1y,&dE1z,&dS1x,&dS1y,&dS1z,&dS2x,&dS2y,&dS2z,v,omega,LNhx,LNhy,LNhz,E1x,E1y,E1z,S1x,S1y,S1z,S2x,S2y,S2z,LNdotS1,LNdotS2,params);

    /* dphi = d \phi / d \hat{t} = M d \phi /dt = M \omega = \hat{\omega} */
    ds = omega;

    dvalues[0]    = ds   ; dvalues[1]     = domega;
    dvalues[2]    = dLNhx; dvalues[3]     = dLNhy ; dvalues[4]    = dLNhz;
    dvalues[5]    = dS1x ; dvalues[6]     = dS1y  ; dvalues[7]    = dS1z ;
    dvalues[8]    = dS2x ; dvalues[9]     = dS2y  ; dvalues[10]   = dS2z ;
    dvalues[11]   = dE1x ; dvalues[12]    = dE1y  ; dvalues[13]   = dE1z ;

    return GSL_SUCCESS;
}

/*
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
	double UNUSED t,
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
    REAL8 v, v2, v3, v4, v7, v11;
    REAL8 LNdotS1, LNdotS2, S1dotS2, S1sq, S2sq;
    REAL8 wspin15 = 0., wspin2 = 0., wspin25 = 0., wspin3 = 0., wspin35 = 0.;

    XLALSimInspiralSpinTaylorTxCoeffs *params
            = (XLALSimInspiralSpinTaylorTxCoeffs*) mparams;

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
    v7 = v4 * v3; v11 = v7 * v4;

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
            wspin35 = params->wdot7S1O*LNdotS1 + params->wdot7S2O*LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
            // Compute 3PN SO correction to domega/dt
            // See Eq. 13 of arXiv:1210.0764 and/or Eq. 3.16 of arXiv:1303.7412
            // Note that S_l/M^2 = (m1/M)^2 chi1 + (m2/M)^2 chi2
            // and Sigma_l/M^2 = (m2/M) chi2 - (m1/M) chi1
            wspin3 = params->wdot6S1O * LNdotS1 + params->wdot6S2O * LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
            // Compute 2.5PN SO correction to domega/dt
            // See Eq. 8.3 of gr-qc/0605140v4
            // Note that S_l/M^2 = (m1/M)^2 chi1 + (m2/M)^2 chi2
            // and Sigma_l/M^2 = (m2/M) chi2 - (m1/M) chi1
            wspin25 = params->wdot5S1O*LNdotS1 + params->wdot5S2O*LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
            // Compute S1-S2 spin-spin term
            S1dotS2 = (S1x*S2x + S1y*S2y + S1z*S2z);
            wspin2 = params->wdot4S1S2 *S1dotS2 + params->wdot4S1OS2O *LNdotS1 * LNdotS2;
            // Compute 2PN QM and self-spin corrections to domega/dt
            // See last line of Eq. 5.17 of arXiv:0812.4413
            // Also note this is equivalent to Eqs. 9c + 9d of astro-ph/0504538
            S1sq = (S1x*S1x + S1y*S1y + S1z*S1z);
            S2sq = (S2x*S2x + S2y*S2y + S2z*S2z);
            wspin2 += params->wdot4QMS1S1 * S1sq
                    + params->wdot4QMS2S2 * S2sq
                    + params->wdot4QMS1OS1O * LNdotS1 * LNdotS1
                    + params->wdot4QMS2OS2O * LNdotS2 * LNdotS2
                    + params->wdot4S1S1 * S1sq
                    + params->wdot4S2S2 * S2sq
                    + params->wdot4S1OS1O * LNdotS1 * LNdotS1
                    + params->wdot4S2OS2O * LNdotS2 * LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            // Compute 1.5PN SO correction to domega/dt
            wspin15 = params->wdot3S1O*LNdotS1 + params->wdot3S2O*LNdotS2;
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %d\n",
                    __func__, params->spinO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    domega  = 3.*params->wdotnewt * v11 / ( params->wdotcoeff[0]
            + v * ( params->wdotcoeff[1]
            + v * ( params->wdotcoeff[2]
            + v * ( params->wdotcoeff[3] + wspin15
            + v * ( params->wdotcoeff[4] + wspin2
            + v * ( params->wdotcoeff[5] + wspin25
            + v * ( params->wdotcoeff[6] + wspin3
                    + params->wdotlogcoeff * log(v)
            + v * ( params->wdotcoeff[7] + wspin35
            + v3 * ( params->wdottidal5pn
            + v2 * ( params->wdottidal6pn ) ) ) ) ) ) ) ) ) );

    XLALSimInspiralSpinDerivatives(&dLNhx,&dLNhy,&dLNhz,&dE1x,&dE1y,&dE1z,&dS1x,&dS1y,&dS1z,&dS2x,&dS2y,&dS2z,v,omega,LNhx,LNhy,LNhz,E1x,E1y,E1z,S1x,S1y,S1z,S2x,S2y,S2z,LNdotS1,LNdotS2,params);

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

/*
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
    else if( fabs(fRef - fStart) < LAL_REAL4_EPS )
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

/*
 * This function evolves the orbital equations for a precessing binary using
 * the \"TaylorT2/T4\" approximant for solving the orbital dynamics
 * (see arXiv:0907.0700 for a review of the various PN approximants).
 *
 * It returns time series of the \"orbital velocity\", orbital phase,
 * and components for both individual spin vectors, the \"Newtonian\"
 * orbital angular momentum (which defines the instantaneous plane)
 * and "E1", a basis vector in the instantaneous orbital plane, as well as
 * the time derivative of the orbital frequency.
 * Note that LNhat and E1 completely specify the instantaneous orbital plane.
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
static int XLALSimInspiralSpinTaylorPNEvolveOrbitIrregularIntervals(
        REAL8Array **yout,              /**< array holding the unevenly sampled output [returned] */
        REAL8 m1,                       /**< mass of companion 1 (kg) */
        REAL8 m2,                       /**< mass of companion 2 (kg) */
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
    LALAdaptiveRungeKutta4Integrator *integrator = NULL;     /* GSL integrator object */
    REAL8 yinit[LAL_NUM_ST4_VARIABLES];       /* initial values of parameters */
    REAL8Array *y;    /* time series of variables returned from integrator */
    /* intermediate variables */
    UINT4 cutlen, len;
//     int sgn, offset;
    REAL8 norm, dtStart, dtEnd, lengths, m1sec, m2sec, Msec, Mcsec, fTerm;
//     LIGOTimeGPS tStart = LIGOTIMEGPSZERO;

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

//     /* Set sign of time step according to direction of integration */
//     if( fEnd < fStart && fEnd != 0. )
//         sgn = -1;
//     else
//         sgn = 1;

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
        XLALSimInspiralSpinTaylorTxCoeffs paramsT1;
        XLALSimInspiralSpinTaylorT1Setup(&paramsT1, m1, m2, fStart, fEnd,
                lambda1, lambda2, quadparam1, quadparam2, spinO, tideO, phaseO);
        params = (void *) &paramsT1;
        //XLALPrintError("XLAL Error - %s: SpinTaylorT1 not implemented yet!\n",
        //        __func__);
        //XLAL_ERROR(XLAL_EINVAL);
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
        integrator = XLALAdaptiveRungeKutta4Init(LAL_NUM_ST4_VARIABLES,
                XLALSimInspiralSpinTaylorT1Derivatives,
                XLALSimInspiralSpinTaylorStoppingTest,
                LAL_ST4_ABSOLUTE_TOLERANCE, LAL_ST4_RELATIVE_TOLERANCE);
        //XLALPrintError("XLAL Error - %s: SpinTaylorT1 not implemented yet!\n",
          //      __func__);
        //XLAL_ERROR(XLAL_EINVAL);
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
    len = XLALAdaptiveRungeKutta4IrregularIntervals(integrator, params, yinit,
            0.0, lengths/Msec, &y);

    intreturn = integrator->returncode;
    XLALAdaptiveRungeKutta4Free(integrator);

    if (!len)
    {
        XLALPrintError("XLAL Error - %s: integration failed with errorcode %d.\n", __func__, intreturn);
        XLAL_ERROR(XLAL_EFUNC);
    }

    /* Print warning about abnormal termination */
    if (intreturn != 0 && intreturn != LALSIMINSPIRAL_ST_TEST_ENERGY
            && intreturn != LALSIMINSPIRAL_ST_TEST_FREQBOUND)
    {
        XLALPrintWarning("XLAL Warning - %s: integration terminated with code %d.\n Waveform parameters were m1 = %e, m2 = %e, s1 = (%e,%e,%e), s2 = (%e,%e,%e), inc = %e.\n",
        __func__, intreturn, m1 * pow(LAL_C_SI, 3.0) / LAL_G_SI / LAL_MSUN_SI,
        m2 * pow(LAL_C_SI, 3.0) / LAL_G_SI / LAL_MSUN_SI, s1x, s1y, s1z, s2x,
        s2y, s2z, acos(lnhatz));
    }

    cutlen = len;

    // Report termination condition and final frequency
    // Will only report this info if '4' bit of lalDebugLevel is 1
    fTerm = y->data[2*len+cutlen-1] / LAL_PI / Msec;
    XLALPrintInfo("XLAL Info - %s: integration terminated with code %d. The final GW frequency reached was %g\n", __func__, intreturn, fTerm);

    *yout = y;

    return XLAL_SUCCESS;
}

/**
 * @addtogroup LALSimInspiralSpinTaylor_c
 * @brief Routines for generating spin precessing Taylor waveforms
 * @{
 */


/**
 * Function to specify the desired orientation of a precessing binary in terms
 * of several angles and then compute the vector components in the so-called
 * \"radiation frame\" (with the z-axis along the direction of propagation) as
 * needed to specify binary configuration for ChooseTDWaveform.
 *
 * Input:
 * thetaJN is the inclination between total angular momentum (J) and the
 * direction of propagation (N)
 * theta1 and theta2 are the inclinations of S1 and S2
 * measured from the Newtonian orbital angular momentum (L_N)
 * phi12 is the difference in azimuthal angles of S1 and S2.
 * chi1, chi2 are the dimensionless spin magnitudes ( \f$0 \le chi1,2 \le 1\f$)
 * phiJL is the azimuthal angle of L_N on its cone about J.
 * m1, m2, f_ref are the component masses and reference GW frequency,
 * they are needed to compute the magnitude of L_N, and thus J.
 *
 * Output:
 * incl - inclination angle of L_N relative to N
 * x, y, z components of E1 (unit vector in the initial orbital plane)
 * x, y, z components S1 and S2 (unit spin vectors times their
 * dimensionless spin magnitudes - i.e. they have unit magnitude for
 * extremal BHs and smaller magnitude for slower spins).
 *
 * NOTE: Here the \"total\" angular momentum is computed as
 * J = L_N + S1 + S2
 * where L_N is the Newtonian orbital angular momentum. In fact, there are
 * PN corrections to L which contribute to J that are NOT ACCOUNTED FOR
 * in this function. This is done so the function does not need to know about
 * the PN order of the system and to avoid subtleties with spin-orbit
 * contributions to L. Also, it is believed that the difference in Jhat
 * with or without these PN corrections to L is quite small.
 *
 * NOTE: fRef = 0 is not a valid choice. If you will pass fRef=0 into
 * ChooseWaveform, then here pass in f_min, the starting GW frequency
 *
 * The various rotations in this transformation are described in more detail
 * in a Mathematica notebook available here:
 * https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/Waveforms/TransformPrecessingInitialConditions
 */
int XLALSimInspiralTransformPrecessingInitialConditions(
		REAL8 *incl,	/**< Inclination angle of L_N (returned) */
		REAL8 *S1x,	/**< S1 x component (returned) */
		REAL8 *S1y,	/**< S1 y component (returned) */
		REAL8 *S1z,	/**< S1 z component (returned) */
		REAL8 *S2x,	/**< S2 x component (returned) */
		REAL8 *S2y,	/**< S2 y component (returned) */
		REAL8 *S2z,	/**< S2 z component (returned) */
		REAL8 thetaJN, 	/**< zenith angle between J and N (rad) */
		REAL8 phiJL,  	/**< azimuthal angle of L_N on its cone about J (rad) */
		REAL8 theta1,  	/**< zenith angle between S1 and LNhat (rad) */
		REAL8 theta2,  	/**< zenith angle between S2 and LNhat (rad) */
		REAL8 phi12,  	/**< difference in azimuthal angle btwn S1, S2 (rad) */
		REAL8 chi1,	/**< dimensionless spin of body 1 */
		REAL8 chi2,	/**< dimensionless spin of body 2 */
		REAL8 m1,	/**< mass of body 1 (kg) */
		REAL8 m2,	/**< mass of body 2 (kg) */
		REAL8 fRef	/**< reference GW frequency (Hz) */
		)
{
	/* Check that fRef is sane */
	if( fRef == 0. )
	{
		XLALPrintError("XLAL Error - %s: fRef=0 is invalid. Please pass in the starting GW frequency instead.\n", __func__);
		XLAL_ERROR(XLAL_EINVAL);
	}

	REAL8 omega0, M, eta, theta0, phi0, Jnorm, tmp1, tmp2;
	REAL8 Jhatx, Jhaty, Jhatz, LNhx, LNhy, LNhz, Jx, Jy, Jz, LNmag;
	REAL8 s1hatx, s1haty, s1hatz, s2hatx, s2haty, s2hatz;
	REAL8 s1x, s1y, s1z, s2x, s2y, s2z;

	/* Starting frame: LNhat is along the z-axis and the unit
	 * spin vectors are defined from the angles relative to LNhat.
	 * Note that we put s1hat in the x-z plane, and phi12
	 * sets the azimuthal angle of s2hat measured from the x-axis.
	 */
	LNhx = 0.;
	LNhy = 0.;
	LNhz = 1.;
	s1hatx = sin(theta1);
	s1haty = 0.;
	s1hatz = cos(theta1);
	s2hatx = sin(theta2) * cos(phi12);
	s2haty = sin(theta2) * sin(phi12);
	s2hatz = cos(theta2);

	/* Define several internal variables needed for magnitudes */
	omega0 = LAL_PI * fRef; // orbital angular frequency at reference point
	m1 *= LAL_G_SI / LAL_C_SI / LAL_C_SI / LAL_C_SI;
	m2 *= LAL_G_SI / LAL_C_SI / LAL_C_SI / LAL_C_SI;
	M = m1 + m2;
	eta = m1 * m2 / M / M;

	/* Define S1, S2, J with proper magnitudes */
	LNmag = pow(M, 5./3.) * eta * pow(omega0, -1./3.);
	s1x = m1 * m1 * chi1 * s1hatx;
	s1y = m1 * m1 * chi1 * s1haty;
	s1z = m1 * m1 * chi1 * s1hatz;
	s2x = m2 * m2 * chi2 * s2hatx;
	s2y = m2 * m2 * chi2 * s2haty;
	s2z = m2 * m2 * chi2 * s2hatz;
	Jx = s1x + s2x;
	Jy = s1y + s2y;
	Jz = LNmag * LNhz + s1z + s2z;

	/* Normalize J to Jhat, find it's angles in starting frame */
	Jnorm = sqrt( Jx*Jx + Jy*Jy + Jz*Jz);
	Jhatx = Jx / Jnorm;
	Jhaty = Jy / Jnorm;
	Jhatz = Jz / Jnorm;
	theta0 = acos(Jhatz);
	phi0 = atan2(Jhaty, Jhatx);

	/* Rotation 1: Rotate about z-axis by -phi0 to put Jhat in x-z plane */
	ROTATEZ(-phi0, LNhx, LNhy, LNhz);
	ROTATEZ(-phi0, s1hatx, s1haty, s1hatz);
	ROTATEZ(-phi0, s2hatx, s2haty, s2hatz);
	ROTATEZ(-phi0, Jhatx, Jhaty, Jhatz);

	/* Rotation 2: Rotate about new y-axis by -theta0
	 * to put Jhat along z-axis
	 */
	ROTATEY(-theta0, LNhx, LNhy, LNhz);
	ROTATEY(-theta0, s1hatx, s1haty, s1hatz);
	ROTATEY(-theta0, s2hatx, s2haty, s2hatz);
	ROTATEY(-theta0, Jhatx, Jhaty, Jhatz);

	/* Rotation 3: Rotate about new z-axis by phiJL to put L at desired
	 * azimuth about J. Note that is currently in x-z plane towards -x
	 * (i.e. azimuth=pi). Hence we rotate about z by phiJL - LAL_PI
	 */
	ROTATEZ(phiJL - LAL_PI, LNhx, LNhy, LNhz);
	ROTATEZ(phiJL - LAL_PI, s1hatx, s1haty, s1hatz);
	ROTATEZ(phiJL - LAL_PI, s2hatx, s2haty, s2hatz);
	ROTATEZ(phiJL - LAL_PI, Jhatx, Jhaty, Jhatz);

	/* Rotation 4: Let N be in x-z plane, inclined from J by thetaJN.
	 * We don't need to explicitly construct it, but rotating the system
	 * about the y-axis by - thetaJN will bring N onto the z-axis.
	 */
	ROTATEY(-thetaJN, LNhx, LNhy, LNhz);
	ROTATEY(-thetaJN, s1hatx, s1haty, s1hatz);
	ROTATEY(-thetaJN, s2hatx, s2haty, s2hatz);
	ROTATEY(-thetaJN, Jhatx, Jhaty, Jhatz);

	/* Rotation 5: We already have N along z. To complete the
	 * transformation into the radiation frame we rotate s.t.
	 * LNhat is in the x-z plane. Thus, if LNhat has azimuth phi0,
	 * we rotate about the z-axis by -phi0.
	 */
	phi0 = atan2(LNhy, LNhx);
	ROTATEZ(-phi0, LNhx, LNhy, LNhz);
	ROTATEZ(-phi0, s1hatx, s1haty, s1hatz);
	ROTATEZ(-phi0, s2hatx, s2haty, s2hatz);
	ROTATEZ(-phi0, Jhatx, Jhaty, Jhatz);

	/* We have completed the transformation. Now find the inclination
	 * of LN relative to N (the final z-axis). */
	*incl = acos(LNhz);
	
	/* Multiply spin unit vectors by chi magnitude (but NOT m_i^2) */
	s1hatx *= chi1;
	s1haty *= chi1;
	s1hatz *= chi1;
	s2hatx *= chi2;
	s2haty *= chi2;
	s2hatz *= chi2;

	/* Set pointers to rotated spin vectors */
	*S1x = s1hatx;
	*S1y = s1haty;
	*S1z = s1hatz;
	*S2x = s2hatx;
	*S2y = s2haty;
	*S2z = s2hatz;

	return XLAL_SUCCESS;
}

/**
 * Function to specify the desired orientation of a precessing binary in terms
 * of several angles and then compute the vector components with respect to
 * orbital angular momentum as needed to specify binary configuration for
 * ChooseTDWaveform.
 *
 * Input:
 * thetaJN is the inclination between total angular momentum (J) and the
 * direction of propagation (N)
 * theta1 and theta2 are the inclinations of S1 and S2
 * measured from the Newtonian orbital angular momentum (L_N)
 * phi12 is the difference in azimuthal angles of S1 and S2.
 * chi1, chi2 are the dimensionless spin magnitudes ( \f$0 \le chi1,2 \le 1\f$)
 * phiJL is the azimuthal angle of L_N on its cone about J.
 * m1, m2, f_ref are the component masses and reference GW frequency,
 * they are needed to compute the magnitude of L_N, and thus J.
 *
 * Output:
 * incl - inclination angle of L_N relative to N
 * x, y, z components of E1 (unit vector in the initial orbital plane)
 * x, y, z components S1 and S2 (unit spin vectors times their
 * dimensionless spin magnitudes - i.e. they have unit magnitude for
 * extremal BHs and smaller magnitude for slower spins).
 *
 * NOTE: Here the \"total\" angular momentum is computed as
 * J = L_N + S1 + S2
 * where L_N is the Newtonian orbital angular momentum. In fact, there are
 * PN corrections to L which contribute to J that are NOT ACCOUNTED FOR
 * in this function. This is done so the function does not need to know about
 * the PN order of the system and to avoid subtleties with spin-orbit
 * contributions to L. Also, it is believed that the difference in Jhat
 * with or without these PN corrections to L is quite small.
 *
 * NOTE: fRef = 0 is not a valid choice. If you will pass fRef=0 into
 * ChooseWaveform, then here pass in f_min, the starting GW frequency
 */
int XLALSimInspiralTransformPrecessingNewInitialConditions(
		REAL8 *incl,	/**< Inclination angle of L_N (returned) */
		REAL8 *S1x,	/**< S1 x component (returned) */
		REAL8 *S1y,	/**< S1 y component (returned) */
		REAL8 *S1z,	/**< S1 z component (returned) */
		REAL8 *S2x,	/**< S2 x component (returned) */
		REAL8 *S2y,	/**< S2 y component (returned) */
		REAL8 *S2z,	/**< S2 z component (returned) */
		const REAL8 thetaJN, 	/**< zenith angle between J and N (rad) */
		const REAL8 phiJL,  	/**< azimuthal angle of L_N on its cone about J (rad) */
		const REAL8 theta1,  	/**< zenith angle between S1 and LNhat (rad) */
		const REAL8 theta2,  	/**< zenith angle between S2 and LNhat (rad) */
		const REAL8 phi12,  	/**< difference in azimuthal angle btwn S1, S2 (rad) */
		const REAL8 chi1,	/**< dimensionless spin of body 1 */
		const REAL8 chi2,	/**< dimensionless spin of body 2 */
		const REAL8 m1_SI,	/**< mass of body 1 (kg) */
		const REAL8 m2_SI,	/**< mass of body 2 (kg) */
		const REAL8 fRef	/**< reference GW frequency (Hz) */
		)
{
	/* Check that fRef is sane */
	if( fRef == 0. )
	{
		XLALPrintError("XLAL Error - %s: fRef=0 is invalid. Please pass in the starting GW frequency instead.\n", __func__);
		XLAL_ERROR(XLAL_EINVAL);
	}

	REAL8 m1, m2, v0, theta0, phi0, Jnorm, tmp1, tmp2;
	REAL8 Jhatx, Jhaty, Jhatz, LNhx, LNhy, LNhz, Jx, Jy, Jz, LNmag;
	REAL8 s1hatx, s1haty, s1hatz, s2hatx, s2haty, s2hatz;
	REAL8 s1x, s1y, s1z, s2x, s2y, s2z;

	/* Starting frame: LNhat is along the z-axis and the unit
	 * spin vectors are defined from the angles relative to LNhat.
	 * Note that we put s1hat in the x-z plane, and phi12
	 * sets the azimuthal angle of s2hat measured from the x-axis.
	 */
	LNhx = 0.;
	LNhy = 0.;
	LNhz = 1.;
	s1hatx = sin(theta1);
	s1haty = 0.;
	s1hatz = cos(theta1);
	s2hatx = sin(theta2) * cos(phi12);
	s2haty = sin(theta2) * sin(phi12);
	s2hatz = cos(theta2);

	/* Define several internal variables needed for magnitudes */
	m1 = m1_SI/LAL_MSUN_SI;
	m2 = m2_SI/LAL_MSUN_SI;
	// v parameter at reference point
	v0 = pow( (m1+m2) * LAL_MTSUN_SI *LAL_PI * fRef, 1./3.);

	/* Define S1, S2, J with proper magnitudes */
	LNmag = m1 * m2 / v0;
	s1x = m1 * m1 * chi1 * s1hatx;
	s1y = m1 * m1 * chi1 * s1haty;
	s1z = m1 * m1 * chi1 * s1hatz;
	s2x = m2 * m2 * chi2 * s2hatx;
	s2y = m2 * m2 * chi2 * s2haty;
	s2z = m2 * m2 * chi2 * s2hatz;
	Jx = s1x + s2x;
	Jy = s1y + s2y;
	Jz = LNmag * LNhz + s1z + s2z;

	/* Normalize J to Jhat, find it's angles in starting frame */
	Jnorm = sqrt( Jx*Jx + Jy*Jy + Jz*Jz);
	Jhatx = Jx / Jnorm;
	Jhaty = Jy / Jnorm;
	Jhatz = Jz / Jnorm;
	theta0 = acos(Jhatz);
	phi0 = atan2(Jhaty, Jhatx);

	/* Rotation 1: Rotate about z-axis by -phi0 to put Jhat in x-z plane */
	ROTATEZ(-phi0, LNhx, LNhy, LNhz);
	ROTATEZ(-phi0, s1hatx, s1haty, s1hatz);
	ROTATEZ(-phi0, s2hatx, s2haty, s2hatz);
	ROTATEZ(-phi0, Jhatx, Jhaty, Jhatz);

	/* Rotation 2: Rotate about new y-axis by -theta0
	 * to put Jhat along z-axis
	 */
	ROTATEY(-theta0, LNhx, LNhy, LNhz);
	ROTATEY(-theta0, s1hatx, s1haty, s1hatz);
	ROTATEY(-theta0, s2hatx, s2haty, s2hatz);
	ROTATEY(-theta0, Jhatx, Jhaty, Jhatz);

	/* Rotation 3: Rotate about new z-axis by phiJL to put L at desired
	 * azimuth about J. Note that is currently in x-z plane towards -x
	 * (i.e. azimuth=pi). Hence we rotate about z by phiJL - LAL_PI
	 */
	ROTATEZ(phiJL - LAL_PI, LNhx, LNhy, LNhz);
	ROTATEZ(phiJL - LAL_PI, s1hatx, s1haty, s1hatz);
	ROTATEZ(phiJL - LAL_PI, s2hatx, s2haty, s2hatz);
	ROTATEZ(phiJL - LAL_PI, Jhatx, Jhaty, Jhatz);

	/* The cosinus of the angle between L and N is the scalar
         * product of the two vectors. We do need to perform additional
         * rotation to compute it
         */
	*incl=acos(-sin(thetaJN)*LNhx+cos(thetaJN)*LNhz); //output

	/* Rotation 4-5: Now J is along z and N in x-z plane, inclined from J
         * by thetaJN. Now we bring L into the z axis to get
         * spin components.
	 */
	REAL8 Nx=-sin(thetaJN);
	REAL8 Ny=0.;
	REAL8 Nz=cos(thetaJN);
	REAL8 thetaLJ = acos(LNhz);
	REAL8 phiLJ   = atan2(LNhy, LNhx);

	/* As a check, one can rotateL too and verify
           the incl angle obtained above.
	  */
	ROTATEZ(-phiLJ, Nx, Ny, Nz);
	ROTATEZ(-phiLJ, s1hatx, s1haty, s1hatz);
	ROTATEZ(-phiLJ, s2hatx, s2haty, s2hatz);
	ROTATEZ(-phiLJ, LNhx, LNhy, LNhz);
	ROTATEY(-thetaLJ, Nx, Ny, Nz);
	ROTATEY(-thetaLJ, s1hatx, s1haty, s1hatz);
	ROTATEY(-thetaLJ, s2hatx, s2haty, s2hatz);
	ROTATEY(-thetaLJ, LNhx, LNhy, LNhz);

	/* Rotation 6: Now L is along z and we have to bring N
	 * in the x-z plane.
	 */
	REAL8 phiN = atan2(Ny, Nx);
	ROTATEZ(-phiN, s1hatx, s1haty, s1hatz);
	ROTATEZ(-phiN, s2hatx, s2haty, s2hatz);
	ROTATEZ(-phiN, LNhx, LNhy, LNhz);
	ROTATEZ(-phiN, Nx, Ny, LNz);

	//One can make the following checks:
	/*printf("LNhat should be along z, N in the x-z plane\n");
	printf("LNhat: %12.4e  %12.4e  %12.4e\n",LNhx,LNhy,LNhz);
	printf("N:     %12.4e  %12.4e  %12.4e\n",Nx,Ny,Nz);
	printf("cos LN = %12.4e vs. %12.4e\n",cos(*incl),LNhx*Nx+LNhy*Ny+LNhz*Nz);
	printf("S1L  %12.4e %12.4e\n",cos(theta1),s1hatz);
	printf("S2L  %12.4e %12.4e\n",cos(theta2),s2hatz);
	printf("S1S2 %12.4e %12.4e\n",sin(theta1)*sin(theta2)*cos(phi12)+cos(theta1)*cos(theta2),s1hatx*s2hatx+s1haty*s2haty+s1hatz*s2hatz);*/

	/* Multiply spin unit vectors by chi magnitude (but NOT m_i^2) */
	s1hatx *= chi1;
	s1haty *= chi1;
	s1hatz *= chi1;
	s2hatx *= chi2;
	s2haty *= chi2;
	s2hatz *= chi2;

	/* Set pointers to rotated spin vectors */
	*S1x = s1hatx;
	*S1y = s1haty;
	*S1z = s1hatz;
	*S2x = s2hatx;
	*S2y = s2haty;
	*S2z = s2hatz;

	return XLAL_SUCCESS;
}

/**
 * Function to specify the desired orientation of the spin components of
 * a precessing binary.
 *
 * Input:
 * x, y, z components S1 and S2 wrt
 * * reference L for axisChoice    = LAL_SIM_INSPIRAL_FRAME_AXIS_TOTAL_J
 * * total J for axisChoice         = LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L
 * * view direction for axisChoice = LAL_SIM_INSPIRAL_FRAME_AXIS_VIEW (default)
 * incl is the angle between
 * * J and N (Jx \f$\propto sin(inc)\f$, Jy=0) for axisChoice = LAL_SIM_INSPIRAL_FRAME_AXIS_TOTAL_J
 * * L and N (Lx \f$\propto sin(inc)\f$, Ly=0) for axisChoice = LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L
 *                                                              LAL_SIM_INSPIRAL_FRAME_AXIS_VIEW (default)
 * m1, m2, f_ref are the component masses and reference GW frequency,
 * they are needed to compute the magnitude of L_N
 *
 * Output:
 * x, y, z components S1 and S2 wrt N
 * inc angle between L and N
 *
 * NOTE: Here the \"total\" angular momentum is computed as
 * J = L_N + S1 + S2
 * where L_N is the Newtonian orbital angular momentum. In fact, there are
 * PN corrections to L which contribute to J that are NOT ACCOUNTED FOR
 * in this function.
 */
int XLALSimInspiralInitialConditionsPrecessingApproxs(
		REAL8 *inc,	/**< inclination angle (returned) */
		REAL8 *S1x,	/**< S1 x component (returned) */
		REAL8 *S1y,	/**< S1 y component (returned) */
		REAL8 *S1z,	/**< S1 z component (returned) */
		REAL8 *S2x,	/**< S2 x component (returned) */
		REAL8 *S2y,	/**< S2 y component (returned) */
		REAL8 *S2z,	/**< S2 z component (returned) */
		const REAL8 inclIn, /**< Inclination angle in input */
		const REAL8 S1xIn,  /**< S1 x component */
		const REAL8 S1yIn,  /**< S1 y component */
		const REAL8 S1zIn,  /**< S1 z component */
		const REAL8 S2xIn,  /**< S2 x component */
		const REAL8 S2yIn,  /**< S2 y component */
		const REAL8 S2zIn,  /**< S2 z component */
		const REAL8 m1,	    /**< mass of body 1 (kg) */
		const REAL8 m2,	    /**< mass of body 2 (kg) */
		const REAL8 fRef,   /**< reference GW frequency (Hz) */
		LALSimInspiralFrameAxis axisChoice  /**< Flag to identify axis wrt which spin components are given. Pass in NULL (or None in python) for default (view) */)
{
  REAL8 LNmag=0.;
  REAL8 LNx,LNy,LNxy2,LNz;
  REAL8 tmp1,tmp2;
  switch (axisChoice) {
  /* FRAME_AXIS_TOTAL_J is used by PhenSpin approximant only,
   * (spins wrt to J, inclIn is the angle between J and N: if
   * N=(0,0,1) Jhat=(sin(inclIn),0,cos(inclIn)))
   */
  case LAL_SIM_INSPIRAL_FRAME_AXIS_TOTAL_J:
    LNmag= m1/LAL_MSUN_SI*m2/LAL_MSUN_SI / cbrt(LAL_PI*fRef*(m1+m2)/LAL_MSUN_SI*LAL_MTSUN_SI);
    LNx  = -S1xIn*pow(m1/LAL_MSUN_SI,2)-S2xIn*pow(m2/LAL_MSUN_SI,2);
    LNy  = -S1yIn*pow(m1/LAL_MSUN_SI,2)-S2yIn*pow(m2/LAL_MSUN_SI,2);
    LNxy2= LNx*LNx+LNy*LNy;
    LNz=0.;
    if (LNmag*LNmag>=LNxy2) {
      LNz=sqrt(LNmag*LNmag-LNxy2);
      if ( LNz<(S1zIn*pow(m1/LAL_MSUN_SI,2)+S2zIn*pow(m2/LAL_MSUN_SI,2)) ) {
	XLALPrintError("** LALSimIMRPSpinInspiralRD error *** for s1 (%12.4e  %12.4e  %12.4e)\n",S1xIn,S1yIn,S1zIn);
	XLALPrintError("                                          s2 (%12.4e  %12.4e  %12.4e)\n",S2xIn,S2yIn,S2zIn);
	XLALPrintError(" wrt to J for m: (%12.4e  %12.4e) and v= %12.4e\n",m1/LAL_MSUN_SI,m2/LAL_MSUN_SI,cbrt(LAL_PI*fRef*(m1+m2)/LAL_MSUN_SI*LAL_MTSUN_SI));
	XLALPrintError(" it is impossible to determine the sign of LNhz\n");
	XLAL_ERROR(XLAL_EDOM);
      }
    }
    else {
      XLALPrintError("** LALSimIMRPSpinInspiralRD error *** unphysical values of s1 (%12.4e  %12.4e  %12.4e)\n",S1xIn,S1yIn,S1zIn);
      XLALPrintError("                                                           s2 (%12.4e  %12.4e  %12.4e)\n",S2xIn,S2yIn,S2zIn);
      XLALPrintError(" wrt to J for m: (%12.4e  %12.4e) and v= %12.4e\n",m1/LAL_MSUN_SI,m2/LAL_MSUN_SI,cbrt(LAL_PI*fRef*(m1+m2)/LAL_MSUN_SI*LAL_MTSUN_SI));
      XLAL_ERROR(XLAL_EDOM);
    }
    *S1x=S1xIn;
    *S1y=S1yIn;
    *S1z=S1zIn;
    *S2x=S2xIn;
    *S2y=S2yIn;
    *S2z=S2zIn;
    ROTATEY(inclIn,*S1x,*S1y,*S1z);
    ROTATEY(inclIn,*S2x,*S2y,*S2z);
    *inc=acos((-sin(inclIn)*LNx+cos(inclIn)*LNz)/LNmag);
    break;
  /* FRAME_AXIS_VIEW (OLD default)
   * (spins wrt to view direction, inclIn is the angle between L and N: if
   * N=(0,0,1) Lhat=(sin(inclIn),0,cos(inclIn)) )
   */
  case LAL_SIM_INSPIRAL_FRAME_AXIS_VIEW:
    *S1x=S1xIn;
    *S1y=S1yIn;
    *S1z=S1zIn;
    *S2x=S2xIn;
    *S2y=S2yIn;
    *S2z=S2zIn;
    *inc=inclIn;
    break;
  /* FRAME_AXIS_ORBITAL_L (default)
   * (spins wrt to L, inclIn is the angle between L and N: if
   * LNhat=(0,0,1) N=(-sin(inclIn),0,cos(inclIn)) )
   */
  case LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L:
  default:
    *S1x=S1xIn*cos(inclIn)+S1zIn*sin(inclIn);
    *S1y=S1yIn;
    *S1z=-S1xIn*sin(inclIn)+S1zIn*cos(inclIn);
    *S2x=S2xIn*cos(inclIn)+S2zIn*sin(inclIn);
    *S2y=S2yIn;
    *S2z=-S2xIn*sin(inclIn)+S2zIn*cos(inclIn);
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
    LALAdaptiveRungeKutta4Integrator *integrator = NULL;     /* GSL integrator object */
    REAL8 yinit[LAL_NUM_ST4_VARIABLES];       /* initial values of parameters */
    REAL8Array *yout;	 /* time series of variables returned from integrator */
    /* intermediate variables */
    UINT4 i, cutlen, len;
    int sgn, offset;
    REAL8 norm, dtStart, dtEnd, lengths, wEnd, m1sec, m2sec, Msec, Mcsec, fTerm;
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
        XLALSimInspiralSpinTaylorTxCoeffs paramsT1;
        XLALSimInspiralSpinTaylorT1Setup(&paramsT1, m1, m2, fStart, fEnd,
                lambda1, lambda2, quadparam1, quadparam2,spinO, tideO, phaseO);
        params = (void *) &paramsT1;
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

//fprintf(stdout,"Integrator %s\n","evoked");
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
        integrator = XLALAdaptiveRungeKutta4Init(LAL_NUM_ST4_VARIABLES,
                XLALSimInspiralSpinTaylorT1Derivatives,
                XLALSimInspiralSpinTaylorStoppingTest,
                LAL_ST4_ABSOLUTE_TOLERANCE, LAL_ST4_RELATIVE_TOLERANCE);

    //    XLALPrintError("XLAL Error - %s: SpinTaylorT1 is getting implemented nw!\n",
     //           __func__);
        //XLAL_ERROR(XLAL_EINVAL);
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
     // fprintf(stdout,"RKintegrator %s\n","started");

    if (!len) 
    {
        XLALPrintError("XLAL Error - %s: integration failed with errorcode %d.\n", __func__, intreturn);
        XLAL_ERROR(XLAL_EFUNC);
    }

    /* Print warning about abnormal termination */
    if (intreturn != 0 && intreturn != LALSIMINSPIRAL_ST_TEST_ENERGY
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

    // Report termination condition and final frequency
    // Will only report this info if '4' bit of lalDebugLevel is 1
    fTerm = yout->data[2*len+cutlen-1] / LAL_PI / Msec;
    XLALPrintInfo("XLAL Info - %s: integration terminated with code %d. The final GW frequency reached was %g\n", __func__, intreturn, fTerm);

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
 * Driver routine to compute a precessing post-Newtonian inspiral waveform
 * with phasing computed from energy balance using the so-called \"T4\" method.
 *
 * This routine allows the user to specify different pN orders
 * for the phasing and amplitude of the waveform.
 *
 * The reference frequency fRef is used as follows:
 * 1) if fRef = 0: The initial values of s1, s2, lnhat and e1 will be the
 * values at frequency fStart. The orbital phase of the last sample is set
 * to phiRef (i.e. phiRef is the "coalescence phase", roughly speaking).
 * THIS IS THE DEFAULT BEHAVIOR CONSISTENT WITH OTHER APPROXIMANTS
 *
 * 2) If fRef = fStart: The initial values of s1, s2, lnhat and e1 will be the
 * values at frequency fStart. phiRef is used to set the orbital phase
 * of the first sample at fStart.
 *
 * 3) If fRef > fStart: The initial values of s1, s2, lnhat and e1 will be the
 * values at frequency fRef. phiRef is used to set the orbital phase at fRef.
 * The code will integrate forwards and backwards from fRef and stitch the
 * two together to create a complete waveform. This allows one to specify
 * the orientation of the binary in-band (or at any arbitrary point).
 * Otherwise, the user can only directly control the initial orientation.
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

int XLALSimInspiralSpinTaylorT1(
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
    Approximant approx = SpinTaylorT1;
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
 * values at frequency fStart. The orbital phase of the last sample is set
 * to phiRef (i.e. phiRef is the "coalescence phase", roughly speaking).
 * THIS IS THE DEFAULT BEHAVIOR CONSISTENT WITH OTHER APPROXIMANTS
 *
 * 2) If fRef = fStart: The initial values of s1, s2, lnhat and e1 will be the
 * values at frequency fStart. phiRef is used to set the orbital phase
 * of the first sample at fStart.
 *
 * 3) If fRef > fStart: The initial values of s1, s2, lnhat and e1 will be the
 * values at frequency fRef. phiRef is used to set the orbital phase at fRef.
 * The code will integrate forwards and backwards from fRef and stitch the
 * two together to create a complete waveform. This allows one to specify
 * the orientation of the binary in-band (or at any arbitrary point).
 * Otherwise, the user can only directly control the initial orientation.
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
 * Compute the physical template family "Q" vectors for a spinning, precessing
 * binary when provided time series of all the dynamical quantities.
 * These vectors always supplied to dominant order.
 *
 * Based on Pan, Buonanno, Chan and Vallisneri PRD69 104017, (see also theses
 * of Diego Fazi and Ian Harry)
 *
 * NOTE: The vectors MUST be given in the so-called radiation frame where
 * Z is the direction of propagation, X is the principal '+' axis and Y = Z x X
 */
int XLALSimInspiralPrecessingPTFQWaveforms(
        REAL8TimeSeries **Q1,     /**< PTF-Q1 waveform [returned] */
        REAL8TimeSeries **Q2,     /**< PTF-Q2 waveform [returned] */
        REAL8TimeSeries **Q3,     /**< PTF-Q2 waveform [returned] */
        REAL8TimeSeries **Q4,     /**< PTF-Q2 waveform [returned] */
        REAL8TimeSeries **Q5,     /**< PTF-Q2 waveform [returned] */
        REAL8TimeSeries *V,       /**< post-Newtonian parameter */
        REAL8TimeSeries *Phi,     /**< orbital phase */
        REAL8TimeSeries *S1x,     /**< Spin1 vector x component */
        REAL8TimeSeries *S1y,     /**< Spin1 vector y component */
        REAL8TimeSeries *S1z,     /**< Spin1 vector z component */
        REAL8TimeSeries *S2x,     /**< Spin2 vector x component */
        REAL8TimeSeries *S2y,     /**< Spin2 vector y component */
        REAL8TimeSeries *S2z,     /**< Spin2 vector z component */
        REAL8TimeSeries *LNhatx,  /**< unit orbital ang. mom. x comp. */
        REAL8TimeSeries *LNhaty,  /**< unit orbital ang. mom. y comp. */
        REAL8TimeSeries *LNhatz,  /**< unit orbital ang. mom. z comp. */
        REAL8TimeSeries *E1x,     /**< orbital plane basis vector x comp. */
        REAL8TimeSeries *E1y,     /**< orbital plane basis vector y comp. */
        REAL8TimeSeries *E1z,     /**< orbital plane basis vector z comp. */
        REAL8 m1,                 /**< mass of companion 1 (kg) */
        REAL8 m2,                 /**< mass of companion 2 (kg) */
        REAL8 r                  /**< distance of source (m) */
        )
{
    REAL8 lnhx, lnhy, lnhz;
    REAL8 e1x, e1y, e1z, e2x, e2y, e2z, nx, ny, nz, lx, ly, lz;
    REAL8 nx2, ny2, nz2, lx2, ly2, lz2;
    REAL8 q1tmp, q2tmp, q3tmp, q4tmp, q5tmp;
    REAL8 M, eta, phi, v, v2, dist, ampfac;
    INT4 idx, len;
    REAL8 sqrt_three = pow(3,0.5);

    /* Macros to check time series vectors */
    LAL_CHECK_VALID_SERIES(V,                   XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(Phi,                 XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S1x,                 XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S1y,                 XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S1z,                 XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S2x,                 XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S2y,                 XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(S2z,                 XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(LNhatx,              XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(LNhaty,              XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(LNhatz,              XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(E1x,                 XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(E1y,                 XLAL_FAILURE);
    LAL_CHECK_VALID_SERIES(E1z,                 XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, Phi,    XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1x,    XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1y,    XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S1z,    XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2x,    XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2y,    XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, S2z,    XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhatx, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhaty, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, LNhatz, XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, E1x,    XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, E1y,    XLAL_FAILURE);
    LAL_CHECK_CONSISTENT_TIME_SERIES(V, E1z,    XLAL_FAILURE);


    /* Allocate polarization vectors and set to 0 */
    *Q1 = XLALCreateREAL8TimeSeries( "PTF_Q_1", &V->epoch,
            0.0, V->deltaT, &lalStrainUnit, V->data->length );
    *Q2 = XLALCreateREAL8TimeSeries( "PTF_Q_2", &V->epoch,
            0.0, V->deltaT, &lalStrainUnit, V->data->length );
    *Q3 = XLALCreateREAL8TimeSeries( "PTF_Q_3", &V->epoch,
            0.0, V->deltaT, &lalStrainUnit, V->data->length );
    *Q4 = XLALCreateREAL8TimeSeries( "PTF_Q_4", &V->epoch,
            0.0, V->deltaT, &lalStrainUnit, V->data->length );
    *Q5 = XLALCreateREAL8TimeSeries( "PTF_Q_5", &V->epoch,
            0.0, V->deltaT, &lalStrainUnit, V->data->length );

    if ( ! Q1 || ! Q2 || !Q3 || !Q4 || !Q5 )
        XLAL_ERROR(XLAL_EFUNC);
    memset((*Q1)->data->data, 0,
            (*Q1)->data->length*sizeof(*(*Q1)->data->data));
    memset((*Q2)->data->data, 0,
            (*Q2)->data->length*sizeof(*(*Q2)->data->data));
    memset((*Q3)->data->data, 0,
            (*Q3)->data->length*sizeof(*(*Q3)->data->data));
    memset((*Q4)->data->data, 0,
            (*Q4)->data->length*sizeof(*(*Q4)->data->data));
    memset((*Q5)->data->data, 0,
            (*Q5)->data->length*sizeof(*(*Q5)->data->data));

    M = m1 + m2;
    eta = m1 * m2 / M / M; // symmetric mass ratio - '\nu' in the paper
    dist = r / LAL_C_SI;   // r (m) / c (m/s) --> dist in units of seconds
    /* convert mass from kg to s, so ampfac ~ M/dist is dimensionless */
    ampfac = 2. * M * LAL_G_SI * pow(LAL_C_SI, -3) * eta / dist;

    /* loop over time steps and compute the Qi */
    len = V->data->length;
    for(idx = 0; idx < len; idx++)
    {
        /* Abbreviated names in lower case for time series at this sample */
        phi  = Phi->data->data[idx];    v = V->data->data[idx];     v2 = v * v;
        lnhx = LNhatx->data->data[idx]; e1x = E1x->data->data[idx];
        lnhy = LNhaty->data->data[idx]; e1y = E1y->data->data[idx];
        lnhz = LNhatz->data->data[idx]; e1z = E1z->data->data[idx];

        /* E2 = LNhat x E1 */
        e2x = lnhy*e1z - lnhz*e1y;
        e2y = lnhz*e1x - lnhx*e1z;
        e2z = lnhx*e1y - lnhy*e1x;

        /* Unit orbital separation vector */
        nx = e1x*cos(phi) + e2x*sin(phi);
        ny = e1y*cos(phi) + e2y*sin(phi);
        nz = e1z*cos(phi) + e2z*sin(phi);

        /* Unit inst. orbital velocity vector */
        lx = e2x*cos(phi) - e1x*sin(phi);
        ly = e2y*cos(phi) - e1y*sin(phi);
        lz = e2z*cos(phi) - e1z*sin(phi);

        /* Powers of vector components */
        nx2 = nx*nx;    ny2 = ny*ny;    nz2 = nz*nz;
        lx2 = lx*lx;    ly2 = ly*ly;    lz2 = lz*lz;

        /* 
         * NOTE: For PTF waveforms, we must use only the dominant amplitude
         * The Q values are computed from equations 13,14,17, 46 + 47 in PBCV or
         * more simply from equations (3.10) in Diego Fazi's thesis.
         * Note that Q5 is simplified from that shown in Fazi's thsis
         * by using traceless condition. As demonstrated in (6.29)
         * in Ian Harry's thesis.
         */

        q1tmp = lx2 - ly2 - nx2 + ny2;
        q2tmp = 2*lx*ly - 2*nx*ny;
        q3tmp = 2*lx*lz - 2*nx*nz;
        q4tmp = 2*ly*lz - 2*ny*nz;
        q5tmp = sqrt_three * (nz2 - lz2);

        /* Fill the output vectors */
        (*Q1)->data->data[idx] = ampfac * v2 * q1tmp;
        (*Q2)->data->data[idx] = ampfac * v2 * q2tmp;
        (*Q3)->data->data[idx] = ampfac * v2 * q3tmp;
        (*Q4)->data->data[idx] = ampfac * v2 * q4tmp;
        (*Q5)->data->data[idx] = ampfac * v2 * q5tmp;
    }
    return XLAL_SUCCESS;
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

/**
 * Driver routine to compute a precessing post-Newtonian inspiral waveform in the Fourier domain
 * with phasing computed from energy balance using the so-called \"T4\" method.
 *
 * This routine allows the user to specify different pN orders
 * for the phasing and amplitude of the waveform.
 *
 * The reference frequency fRef is used as follows:
 * 1) if fRef = 0: The initial values of s1, s2, lnhat and e1 will be the
 * values at frequency fStart. The orbital phase of the last sample is set
 * to phiRef (i.e. phiRef is the "coalescence phase", roughly speaking).
 * THIS IS THE DEFAULT BEHAVIOR CONSISTENT WITH OTHER APPROXIMANTS
 *
 * 2) If fRef = fStart: The initial values of s1, s2, lnhat and e1 will be the
 * values at frequency fStart. phiRef is used to set the orbital phase
 * of the first sample at fStart.
 *
 * 3) If fRef > fStart: The initial values of s1, s2, lnhat and e1 will be the
 * values at frequency fRef. phiRef is used to set the orbital phase at fRef.
 * The code will integrate forwards and backwards from fRef and stitch the
 * two together to create a complete waveform. This allows one to specify
 * the orientation of the binary in-band (or at any arbitrary point).
 * Otherwise, the user can only directly control the initial orientation.
 *
 * 4) fRef < 0 or fRef >= Schwarz. ISCO are forbidden and the code will abort.
 *
 * It is recommended, but not necessary to set fStart slightly smaller than fMin,
 * e.g. fStart = 9.5 for fMin = 10.
 *
 * The returned Fourier series are set so that the Schwarzschild ISCO frequency
 * corresponds to t = 0 as closely as possible.
 *
 */
int XLALSimInspiralSpinTaylorT4Fourier(
        COMPLEX16FrequencySeries **hplus,        /**< +-polarization waveform */
        COMPLEX16FrequencySeries **hcross,       /**< x-polarization waveform */
        REAL8 fMin,                     /**< minimum frequency of the returned series */
        REAL8 fMax,                     /**< maximum frequency of the returned series */
        REAL8 deltaF,                   /**< frequency interval of the returned series */
        INT4 kMax,                      /**< k_max as described in defined arXiv: 1408.5158 (min 0, max 10). */
        REAL8 phiRef,                   /**< orbital phase at reference pt. */
        REAL8 v0,                       /**< tail gauge term (default = 1) */
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
        INT4 phaseO,                     /**< twice PN phase order */
        INT4 amplitudeO,                /**< twice PN amplitude order */
        INT4 phiRefAtEnd                /**< whether phiRef corresponds to the end of the inspiral */
        )
{
    Approximant approx = SpinTaylorT4;
    int n = XLALSimInspiralSpinTaylorDriverFourier(hplus, hcross, fMin, fMax,
            deltaF, kMax, phiRef, v0,
            m1, m2, fStart, fRef, r, s1x, s1y, s1z, s2x, s2y, s2z,
            lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2,
            quadparam1, quadparam2, spinO, tideO, phaseO, amplitudeO, approx, phiRefAtEnd);

    return n;
}

/**
 * Driver routine to compute a precessing post-Newtonian inspiral waveform in the Fourier domain
 * with phasing computed from energy balance using the so-called \"T2\" method.
 *
 * This routine allows the user to specify different pN orders
 * for the phasing and amplitude of the waveform.
 *
 * The reference frequency fRef is used as follows:
 * 1) if fRef = 0: The initial values of s1, s2, lnhat and e1 will be the
 * values at frequency fStart. The orbital phase of the last sample is set
 * to phiRef (i.e. phiRef is the "coalescence phase", roughly speaking).
 * THIS IS THE DEFAULT BEHAVIOR CONSISTENT WITH OTHER APPROXIMANTS
 *
 * 2) If fRef = fStart: The initial values of s1, s2, lnhat and e1 will be the
 * values at frequency fStart. phiRef is used to set the orbital phase
 * of the first sample at fStart.
 *
 * 3) If fRef > fStart: The initial values of s1, s2, lnhat and e1 will be the
 * values at frequency fRef. phiRef is used to set the orbital phase at fRef.
 * The code will integrate forwards and backwards from fRef and stitch the
 * two together to create a complete waveform. This allows one to specify
 * the orientation of the binary in-band (or at any arbitrary point).
 * Otherwise, the user can only directly control the initial orientation.
 *
 * 4) fRef < 0 or fRef >= Schwarz. ISCO are forbidden and the code will abort.
 *
 * It is recommended, but not necessary to set fStart slightly smaller than fMin,
 * e.g. fStart = 9.5 for fMin = 10.
 *
 * The returned Fourier series are set so that the Schwarzschild ISCO frequency
 * corresponds to t = 0 as closely as possible.
 *
 */
int XLALSimInspiralSpinTaylorT2Fourier(
        COMPLEX16FrequencySeries **hplus,        /**< +-polarization waveform */
        COMPLEX16FrequencySeries **hcross,       /**< x-polarization waveform */
        REAL8 fMin,                     /**< minimum frequency of the returned series */
        REAL8 fMax,                     /**< maximum frequency of the returned series */
        REAL8 deltaF,                   /**< frequency interval of the returned series */
        INT4 kMax,                      /**< k_max as described in arXiv: 1408.5158 (min 0, max 10). */
        REAL8 phiRef,                   /**< orbital phase at reference pt. */
        REAL8 v0,                       /**< tail gauge term (default = 1) */
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
        INT4 phaseO,                     /**< twice PN phase order */
        INT4 amplitudeO,                /**< twice PN amplitude order */
        INT4 phiRefAtEnd                /**< whether phiRef corresponds to the end of the inspiral */
        )
{
    Approximant approx = SpinTaylorT2;
    int n = XLALSimInspiralSpinTaylorDriverFourier(hplus, hcross, fMin,
            fMax, deltaF, kMax, phiRef, v0,
            m1, m2, fStart, fRef, r, s1x, s1y, s1z, s2x, s2y, s2z,
            lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2,
            quadparam1, quadparam2, spinO, tideO, phaseO, amplitudeO, approx, phiRefAtEnd);

    return n;
}

/** @} */
