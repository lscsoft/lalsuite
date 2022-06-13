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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

#include <math.h>
#include <lal/Units.h>
#include <lal/LALConstants.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALAdaptiveRungeKuttaIntegrator.h>
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


static int XLALSimInspiralVectorCrossProduct(REAL8 **vout, REAL8 v1x, REAL8 v1y, REAL8 v1z, REAL8 v2x, REAL8 v2y, REAL8 v2z){
    (*vout) = (double *) LALMalloc(sizeof(double) * 3);
    (*vout)[0]=v1y*v2z-v1z*v2y;
    (*vout)[1]=v1z*v2x-v1x*v2z;
    (*vout)[2]=v1x*v2y-v1y*v2x;
    return XLAL_SUCCESS;
}

static REAL8 cdot(REAL8 v1x, REAL8 v1y, REAL8 v1z, REAL8 v2x, REAL8 v2y, REAL8 v2z){
  return v1x*v2x+v1y*v2y+v1z*v2z;
}

static REAL8 normsq(REAL8 vx, REAL8 vy, REAL8 vz){
  return vx*vx+vy*vy+vz*vz;
}


/* Declarations of static functions - defined below */
static int XLALSimInspiralSpinTaylorStoppingTest(double t,
	const double values[], double dvalues[], void *mparams);
static int XLALSimInspiralSpinTaylorT1DerivativesAvg(double t,
	const double values[], double dvalues[], void *mparams);
static int XLALSimInspiralSpinTaylorT5DerivativesAvg(double t,
	const double values[], double dvalues[], void *mparams);
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
    LALDict *LALparams, INT4 phaseO, INT4 amplitudeO, Approximant approx, INT4 phiRefAtEnd);


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

/* Setup coefficients of the energy function and 
 * of the spin derivative equations, which are common
 * to all members of the SpinTaylor family.
 */
INT4 XLALSimSpinTaylorEnergySpinDerivativeSetup(
    XLALSimInspiralSpinTaylorTxCoeffs **params, /**< OUTPUT */
    const REAL8 m1_SI,                    /**< mass of body 1 (kg) */
    const REAL8 m2_SI,                    /**< mass of body 2 (kg) */
    const REAL8 fStart,                   /**< Starting GW freq. (Hz) */
    const REAL8 fEnd,                     /**< Ending GW freq. (Hz), 0 means integrate forwards as far as possible */
    const LALSimInspiralSpinOrder  spinO, /**< twice PN order of spin effects */
    const LALSimInspiralTidalOrder tideO, /**< twice PN order of tidal effects */
    const INT4 phaseO,                    /**< twice PN order of spin effects */
    const REAL8 lambda1,                  /**< (tidal deformability of mass 1) / (mass of body 1)^5 (dimensionless) */
    const REAL8 lambda2,                  /**< (tidal deformability of mass 2) / (mass of body 2)^5 (dimensionless) */
    const REAL8 quadparam1,               /**< phenom. parameter describing induced quad. moment of body 1 (=1 for BHs, ~2-12 for NSs) */
    const REAL8 quadparam2,               /**< phenom. parameter describing induced quad. moment of body 2 (=1 for BHs, ~2-12 for NSs) */
    const INT4 lscorr,                    /**< flag to include spin corrections to orbital angular momentum */
    const INT4 phenomtp                   /**< flag for using spinO=7 and not spinO=6 terms with orbital-averaged quantities for phenomtphm approx */)
{

  REAL8 m1=m1_SI/LAL_MSUN_SI;
  REAL8 m2=m2_SI/LAL_MSUN_SI;
  REAL8 M=m1+m2;
  REAL8 m1M=m1/M;
  REAL8 m2M=m2/M;
  REAL8 eta=m1M*m2M;
  (*params)->wdotnewt = XLALSimInspiralTaylorT4wdot_0PNCoeff(eta);
  (*params)->M = M;
  (*params)->Mchirp = M*pow(eta,0.6);
  (*params)->m1M = m1M;
  (*params)->m2M = m2M;
  (*params)->eta = eta;
  (*params)->fStart = fStart;
  (*params)->fEnd = fEnd;
  (*params)->phaseO = phaseO;
  (*params)->spinO  = spinO;
  (*params)->tideO  = tideO;
  (*params)->lscorr = lscorr;
  (*params)->phenomtp = phenomtp;

  /* Set the non-dynamical coefficients of spin-independent
   * Energy terms
   */
  switch( (*params)->phaseO ) {
    case -1: /* highest available PN order */
    case 8:
      /* case LAL_PNORDER_THREE_POINT_FIVE: */
    case 7:
      (*params)->Ecoeff[7] = 0.;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
      /* case LAL_PNORDER_THREE: */
    case 6:
      (*params)->Ecoeff[6] = XLALSimInspiralPNEnergy_6PNCoeff(eta);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
      /* case LAL_PNORDER_TWO_POINT_FIVE: */
    case 5:
      (*params)->Ecoeff[5] = 0.;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
      /* case LAL_PNORDER_TWO: */
    case 4:
      (*params)->Ecoeff[4] = XLALSimInspiralPNEnergy_4PNCoeff(eta);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
      /*case LAL_PNORDER_ONE_POINT_FIVE:*/
    case 3:
      (*params)->Ecoeff[3] = 0.;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
      /*case LAL_PNORDER_ONE:*/
    case 2:
      (*params)->Ecoeff[2] = XLALSimInspiralPNEnergy_2PNCoeff(eta);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
      /*case LAL_PNORDER_HALF:*/
    case 1:
      (*params)->Ecoeff[1] = 0.;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
      /*case LAL_PNORDER_NEWTONIAN:*/
    case 0:
      (*params)->Ecoeff[0] = 1.;
      break;
    default:
      XLALPrintError("XLAL Error - %s: Invalid phase. PN order %d\n",
		     __func__, (*params)->phaseO );
      XLAL_ERROR(XLAL_EINVAL);
      break;
  }

  /* Set the non-dynamical coefficients of spin-dependent
   * Energy terms
   */
  switch( (*params)->spinO )
    {
    case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
    case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
      // 3.5PN spin-orbit
      // For historical reasons phenomtp uses NNNNL spin-orbit terms but not NNNL spin^2 terms
      if ( phenomtp ) {
	// Energy coefficients
	(*params)->E7S1O  = XLALSimInspiralPNEnergy_7PNSOCoeff(m1M); // Coefficient of S1.LN
	(*params)->E7S2O  = XLALSimInspiralPNEnergy_7PNSOCoeff(m2M); // Coefficient of S2.LN
	// Sdot coefficients
	(*params)->S1dot7S2    = XLALSimInspiralSpinDot_7PNCoeff(m1M); // Coefficient of S2 x S1
	(*params)->S2dot7S1    = XLALSimInspiralSpinDot_7PNCoeff(m2M); // Coefficient of S1 x S2
      }
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
    case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
      // 3PN spin-spin
      // For historical reasons phenomtp uses NNNNL spin-orbit terms but not NNNL spin^2 terms
      if (!(phenomtp)) {
	// Energy orbit-averaged coefficients
	(*params)->E6S1S2Avg   = XLALSimInspiralPNEnergy_6PNS1S2CoeffAvg(eta);    // Coefficient of S1.S2
	(*params)->E6S1OS2OAvg = XLALSimInspiralPNEnergy_6PNS1OS2OCoeffAvg(eta);  // Coefficient of S1.LN S2.LN
	(*params)->E6S1S1Avg   = XLALSimInspiralPNEnergy_6PNS1S1CoeffAvg(m1M);    // Coefficient of S1.S1
	(*params)->E6S2S2Avg   = XLALSimInspiralPNEnergy_6PNS1S1CoeffAvg(m2M);    // Coefficient of S2.S2
	(*params)->E6S1OS1OAvg = XLALSimInspiralPNEnergy_6PNS1OS1OCoeffAvg(m1M);  // Coefficient of (S1.LN)^2
	(*params)->E6S2OS2OAvg = XLALSimInspiralPNEnergy_6PNS1OS1OCoeffAvg(m2M);  // Coefficient of (S2.LN)^2
	(*params)->E6QMS1S1Avg = quadparam1 * XLALSimInspiralPNEnergy_6PNQMS1S1CoeffAvg(m1M);  // Coefficient of quadrupole-monopole S1-squared
	(*params)->E6QMS2S2Avg = quadparam2 * XLALSimInspiralPNEnergy_6PNQMS1S1CoeffAvg(m2M);  // Coefficient of quadrupole-monopole S2-squared
	(*params)->E6QMS1OS1OAvg = quadparam1 * XLALSimInspiralPNEnergy_6PNQMS1OS1OCoeffAvg(m1M);  // Coefficient of quadrupole-monopole S1.L-squared
	(*params)->E6QMS2OS2OAvg = quadparam2 * XLALSimInspiralPNEnergy_6PNQMS1OS1OCoeffAvg(m2M);  // Coefficient of quadrupole-monopole S2.L-squared
	// Sdot orbit-averaged coefficients
	(*params)->S1dot6S2Avg   = XLALSimInspiralSpinDot_6PNS2CoeffAvg(m1M); // Avg-Coefficient of S2xS1
	(*params)->S2dot6S1Avg   = XLALSimInspiralSpinDot_6PNS2CoeffAvg(m2M); // Avg-Coefficient of S1xS2
	(*params)->S1dot6S2OAvg  = XLALSimInspiralSpinDot_6PNS2OCoeffAvg(m1M); // Avg-Coefficient of L.S2 LxS1
	(*params)->S1dot6S1OAvg  = XLALSimInspiralSpinDot_6PNS1OCoeffAvg(m1M); // Avg-Coefficient of L.S1 LxS1
	(*params)->S2dot6S1OAvg  = XLALSimInspiralSpinDot_6PNS2OCoeffAvg(m2M); // Avg-Coefficient of L.S1 LxS2
	(*params)->S2dot6S2OAvg  = XLALSimInspiralSpinDot_6PNS1OCoeffAvg(m2M); // Avg-Coefficient of L.S2 LxS2
	(*params)->S1dot6QMS1OAvg = quadparam1 * XLALSimInspiralSpinDot_6PNQMSOCoeffAvg(m1M); // Avg-Coefficient of QM L.S1 LxS1
	(*params)->S2dot6QMS2OAvg = quadparam2 * XLALSimInspiralSpinDot_6PNQMSOCoeffAvg(m2M); // Avg-Coefficient of QM L.S2 LxS2
	(*params)->omegashiftS1     = XLALSimInspiralLDot_3PNSOCoeff(m1M);
	(*params)->omegashiftS2     = XLALSimInspiralLDot_3PNSOCoeff(m2M);
      }
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
    case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
      // Energy coefficients
      (*params)->E5S1O    = XLALSimInspiralPNEnergy_5PNSOCoeff(m1M); // Coefficient of S1.LN
      (*params)->E5S2O    = XLALSimInspiralPNEnergy_5PNSOCoeff(m2M); // Coefficient of S2.LN
      // Sdot coefficients
      (*params)->S1dot5      = XLALSimInspiralSpinDot_5PNCoeff(m1M); // Coefficient of LNxS1
      (*params)->S2dot5      = XLALSimInspiralSpinDot_5PNCoeff(m2M); // Coefficient of LNxS2
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
    case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
      /* One considers averaged spin^2 terms according to
       * <n^i n^j> = 1/2(delta^{ij}-L^iL^j) +O(v^3), with n^i the unit
       * vector pointing from BH 1 to BH 2.
       * Such average introduce a v^3 error with respect to the leading order
       * the spins enter at,
       * see e.g. app. B of PRD80 (2009) 044010, arXiv:0812.4413.
       */
      // 2PN spin-spin averaged terms
      (*params)->E4S1S2Avg      = XLALSimInspiralPNEnergy_4PNS1S2CoeffAvg(eta);  // Coefficient of S1.S2
      (*params)->E4S1OS2OAvg    = XLALSimInspiralPNEnergy_4PNS1OS2OCoeffAvg(eta); // Coefficient of S1.LN S2.LN
      // 2PN quadrupole-monopole averaged self-spin terms
      (*params)->E4QMS1S1Avg    = quadparam1 * XLALSimInspiralPNEnergy_4PNQMS1S1CoeffAvg(m1M);  // Coefficient of quad-monop term S1.S1
      (*params)->E4QMS1OS1OAvg  = quadparam1 * XLALSimInspiralPNEnergy_4PNQMS1OS1OCoeffAvg(m1M); // Coefficient of quad-monop term (S1.LN)^2
      (*params)->E4QMS2S2Avg    = quadparam2 * XLALSimInspiralPNEnergy_4PNQMS1S1CoeffAvg(m2M);  // Coefficient of quad-monop term S2.S2
      (*params)->E4QMS2OS2OAvg  = quadparam2 * XLALSimInspiralPNEnergy_4PNQMS1OS1OCoeffAvg(m2M); // Coefficient of quad-monop term (S2.LN)^2
      // Spin derivative
      (*params)->S1dot4S2Avg     = XLALSimInspiralSpinDot_4PNS2CoeffAvg;  // Coefficient of S2xS1 in S1dot and S1xS2 in S2dot
      (*params)->S1dot4S2OAvg    = XLALSimInspiralSpinDot_4PNS2OCoeffAvg; // Coefficient of (LN.S1) S2xS1 in S1dot and (LN.S2) S1xS2 in S2dot
      (*params)->S1dot4QMS1OAvg = quadparam1 * XLALSimInspiralSpinDot_4PNQMSOCoeffAvg(m1M); // Coefficient of quad-monop. term (S1.LN) LNxS1
      (*params)->S2dot4QMS2OAvg = quadparam2 * XLALSimInspiralSpinDot_4PNQMSOCoeffAvg(m2M); // Coefficient of quad-monop. term (S2.LN) LNxS2

#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
    case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
      (*params)->E3S1O        = XLALSimInspiralPNEnergy_3PNSOCoeff(m1M); // Coefficient of S1.LN
      (*params)->E3S2O        = XLALSimInspiralPNEnergy_3PNSOCoeff(m2M); // Coefficient of S2.LN
      (*params)->S1dot3       = XLALSimInspiralSpinDot_3PNCoeff(m1M);    // Coefficient of LNxS1
      (*params)->S2dot3       = XLALSimInspiralSpinDot_3PNCoeff(m2M);    // Coefficient of LNxS2
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
    case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
    case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
    case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
      break;
    default:
      XLALPrintError("XLAL Error - %s: Invalid spin PN order %d\n",
		     __func__, (*params)->spinO );
      XLAL_ERROR(XLAL_EINVAL);
      break;
    }

    switch( (*params)->tideO )
    {
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL:
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN:
            (*params)->Etidal12 = lambda1 * XLALSimInspiralPNEnergy_12PNTidalCoeff(m1M) + lambda2 * XLALSimInspiralPNEnergy_12PNTidalCoeff(m2M);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
	    (*params)->Etidal10 = lambda1 * XLALSimInspiralPNEnergy_10PNTidalCoeff(m1M) + lambda2 * XLALSimInspiralPNEnergy_10PNTidalCoeff(m2M);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid tidal PN order %d\n",
                    __func__, (*params)->tideO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

  return XLAL_SUCCESS;
} // End of XLALSimSpinTaylorEnergySpinDerivativeSetup()

/**
 * Function computing spin precession equations, common to all member of the SpinTaylor family.
 * Orbit averaged dynamics.
 */

INT4 XLALSimInspiralSpinDerivativesAvg(REAL8 *dLNhx,
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

  *dLNhx=0.;
  *dLNhy=0.;
  *dLNhz=0.;
  *dE1x=0.;
  *dE1y=0.;
  *dE1z=0.;
  *dS1x=0.;
  *dS1y=0.;
  *dS1z=0.;
  *dS2x=0.;
  *dS2y=0.;
  *dS2z=0.;
  INT4 lscorr=params->lscorr;
  
  REAL8 dLNhatx=0.; // Derivative of the Newtonian angular momentum unit vector
  REAL8 dLNhaty=0.;
  REAL8 dLNhatz=0.;
  
  const REAL8 eta=params->eta;
  const REAL8 LN0mag=eta/v;
  REAL8 LNmag=LN0mag;

  const REAL8 cS1  = XLALSimInspiralL_3PNSicoeffAvg(params->m1M);
  const REAL8 cS1L = XLALSimInspiralL_3PNSiLcoeffAvg(params->m1M);
  const REAL8 cS2  = XLALSimInspiralL_3PNSicoeffAvg(params->m2M);
  const REAL8 cS2L = XLALSimInspiralL_3PNSiLcoeffAvg(params->m2M);
  
  /*
   * dS1
   * d S_1 / d \hat{t} = M * d S_1 / dt = \Omega_{S1,S2,LN,v} x S_1
   * However, the papers referenced below uses spin variables which are M^2 times our spins
   */

  if ( (params->spinO>=3) || (params->spinO<0) ) {
    /* dS1,2 leading terms: eq. (8) of gr-qc/0405090.*/
    const REAL8 v2=v*v;
    const REAL8 omega=v2*v;
    const REAL8 v5=omega*v2;

    REAL8 *LNhcS1=NULL;
    XLALSimInspiralVectorCrossProduct(&LNhcS1,LNhx,LNhy,LNhz,S1x,S1y,S1z);
    const REAL8 dS1xL = params->S1dot3 * v5 * LNhcS1[0];
    const REAL8 dS1yL = params->S1dot3 * v5 * LNhcS1[1];
    const REAL8 dS1zL = params->S1dot3 * v5 * LNhcS1[2];
    
    REAL8 *LNhcS2=NULL;
    XLALSimInspiralVectorCrossProduct(&LNhcS2,LNhx,LNhy,LNhz,S2x,S2y,S2z);
    const REAL8 dS2xL = params->S2dot3 * v5 * LNhcS2[0];
    const REAL8 dS2yL = params->S2dot3 * v5 * LNhcS2[1];
    const REAL8 dS2zL = params->S2dot3 * v5 * LNhcS2[2];
    
    *dS1x =dS1xL;
    *dS1y =dS1yL;
    *dS1z =dS1zL;

    *dS2x =dS2xL;
    *dS2y =dS2yL;
    *dS2z =dS2zL;

    const REAL8 dLxL=-(dS1xL+ dS2xL) ;
    const REAL8 dLyL=-(dS1yL+ dS2yL);
    const REAL8 dLzL=-(dS1zL+ dS2zL);
    
    /* dLNhat should be normalized by dividing by the modulus of LN.
     * It will be done at the end of this fucntion. */
    dLNhatx = dLxL;
    dLNhaty = dLyL;
    dLNhatz = dLzL;
    
    if ( (params->spinO>=4) || (params->spinO<0.) ) {
      /* dS1,2 NLO term (v x leading), Spin^2 terms */
      REAL8 omega2=omega*omega;
      REAL8 *S1cS2= NULL;
      XLALSimInspiralVectorCrossProduct(&S1cS2,S1x,S1y,S1z,S2x,S2y,S2z);
      /* S1S2 contribution, see. eq. 4.17 of Phys.Rev. D52 (1995) 821-847, arxiv/gr-qc/9506022 */
      REAL8 dS1xNL = omega2 * (-params->S1dot4S2Avg * S1cS2[0] + params->S1dot4S2OAvg * LNhdotS2 * LNhcS1[0]);
      REAL8 dS1yNL = omega2 * (-params->S1dot4S2Avg * S1cS2[1] + params->S1dot4S2OAvg * LNhdotS2 * LNhcS1[1]);
      REAL8 dS1zNL = omega2 * (-params->S1dot4S2Avg * S1cS2[2] + params->S1dot4S2OAvg * LNhdotS2 * LNhcS1[2]);
      REAL8 dS2xNL = omega2 * ( params->S1dot4S2Avg * S1cS2[0] + params->S1dot4S2OAvg * LNhdotS1 * LNhcS2[0]);
      REAL8 dS2yNL = omega2 * ( params->S1dot4S2Avg * S1cS2[1] + params->S1dot4S2OAvg * LNhdotS1 * LNhcS2[1]);
      REAL8 dS2zNL = omega2 * ( params->S1dot4S2Avg * S1cS2[2] + params->S1dot4S2OAvg * LNhdotS1 * LNhcS2[2]);

      /* QM S1S1 contribution */
      dS1xNL += omega2 * params->S1dot4QMS1OAvg * LNhdotS1 * LNhcS1[0]; // v6 terms
      dS1yNL += omega2 * params->S1dot4QMS1OAvg * LNhdotS1 * LNhcS1[1];
      dS1zNL += omega2 * params->S1dot4QMS1OAvg * LNhdotS1 * LNhcS1[2];
      
      /* QM S2S2 contribution */
      dS2xNL += omega2 * params->S2dot4QMS2OAvg * LNhdotS2 * LNhcS2[0];
      dS2yNL += omega2 * params->S2dot4QMS2OAvg * LNhdotS2 * LNhcS2[1];
      dS2zNL += omega2 * params->S2dot4QMS2OAvg * LNhdotS2 * LNhcS2[2];

      *dS1x+= dS1xNL;
      *dS1y+= dS1yNL;
      *dS1z+= dS1zNL;

      *dS2x+= dS2xNL;
      *dS2y+= dS2yNL;
      *dS2z+= dS2zNL;

      const REAL8 dLxNL = -(dS1xNL+dS2xNL);
      const REAL8 dLyNL = -(dS1yNL+dS2yNL);
      const REAL8 dLzNL = -(dS1zNL+dS2zNL);

      /* dLNhat should be orthogonal to LNhat
       * Here dL is not, dLNhat will be projected onto the space
       * orthogonal to LNhat at the end of this function.
       */
      dLNhatx+=dLxNL;
      dLNhaty+=dLyNL;
      dLNhatz+=dLzNL;

      if ( (params->spinO>=5) || (params->spinO<0) ) {

	/* dLNhat is parallel to dL only until NLO (v x leading),
	 * since up to this order LNhat = Lhat.
	 * at this NNL order we have to include spin dependent terms in the orbital angular momentum.
	 */

	const REAL8 L1PN=XLALSimInspiralL_2PN(eta);
	REAL8 v7=omega2*v;
	LNmag+=LN0mag*v2*L1PN;

	/* dS1,2 NNLO, eq. 7.8 of Blanchet et al. gr-qc/0605140 */
	const REAL8 dS1xNNL = params->S1dot5 * v7 * LNhcS1[0]; // v7 terms
	const REAL8 dS1yNNL = params->S1dot5 * v7 * LNhcS1[1];
	const REAL8 dS1zNNL = params->S1dot5 * v7 * LNhcS1[2];

	const REAL8 dS2xNNL = params->S2dot5 * v7 * LNhcS2[0];
	const REAL8 dS2yNNL = params->S2dot5 * v7 * LNhcS2[1];
	const REAL8 dS2zNNL = params->S2dot5 * v7 * LNhcS2[2];

	*dS1x+= dS1xNNL;
	*dS1y+= dS1yNNL;
	*dS1z+= dS1zNNL;

	*dS2x+= dS2xNNL;
	*dS2y+= dS2yNNL;
	*dS2z+= dS2zNNL;

	dLNhatx-= dS1xNNL + dS2xNNL;
	dLNhaty-= dS1yNNL + dS2yNNL;
	dLNhatz-= dS1zNNL + dS2zNNL;

	if (lscorr) {
	  /* At this order only the derivative of the spins are included into
           * spin-dependent terms of orbital angular momentum.
           * LO of spin derivative is v^5, derivative \hat L is v^6 at LO.
	   * Spin derivative at leading order are perpendicular to L.
	   */	  
	  dLNhatx -= eta*v2*( cS1*dS1xL + cS2*dS2xL);
	  dLNhaty -= eta*v2*( cS1*dS1yL + cS2*dS2yL);
	  dLNhatz -= eta*v2*( cS1*dS1zL + cS2*dS2zL);
	}
	
	if ( (params->spinO>=6) || (params->spinO<0) ) {

	  if (!(params->phenomtp)) {

	      REAL8 v8=omega2*v2;

	      REAL8 dS1xN3L=v8*(-params->S1dot6S2Avg*S1cS2[0] + (params->S1dot6S1OAvg*LNhdotS1 + params->S1dot6S2OAvg*LNhdotS2)*LNhcS1[0]);
	      REAL8 dS1yN3L=v8*(-params->S1dot6S2Avg*S1cS2[1] + (params->S1dot6S1OAvg*LNhdotS1 + params->S1dot6S2OAvg*LNhdotS2)*LNhcS1[1]);
	      REAL8 dS1zN3L=v8*(-params->S1dot6S2Avg*S1cS2[2] + (params->S1dot6S1OAvg*LNhdotS1 + params->S1dot6S2OAvg*LNhdotS2)*LNhcS1[2]);
	      REAL8 dS2xN3L=v8*( params->S2dot6S1Avg*S1cS2[0] + (params->S2dot6S1OAvg*LNhdotS1 + params->S2dot6S2OAvg*LNhdotS2)*LNhcS2[0]);
	      REAL8 dS2yN3L=v8*( params->S2dot6S1Avg*S1cS2[1] + (params->S2dot6S1OAvg*LNhdotS1 + params->S2dot6S2OAvg*LNhdotS2)*LNhcS2[1]);
	      REAL8 dS2zN3L=v8*( params->S2dot6S1Avg*S1cS2[2] + (params->S2dot6S1OAvg*LNhdotS1 + params->S2dot6S2OAvg*LNhdotS2)*LNhcS2[2]);

	      // QM S1S1 terms
	      dS1xN3L += v8*(params->S1dot6QMS1OAvg * LNhdotS1 * LNhcS1[0] );
	      dS1yN3L += v8*(params->S1dot6QMS1OAvg * LNhdotS1 * LNhcS1[1] );
	      dS1zN3L += v8*(params->S1dot6QMS1OAvg * LNhdotS1 * LNhcS1[2] );

	      // QM S2S2 terms
	      dS2xN3L += v8*(params->S2dot6QMS2OAvg * LNhdotS2 * LNhcS2[0] );
	      dS2yN3L += v8*(params->S2dot6QMS2OAvg * LNhdotS2 * LNhcS2[1] );
	      dS2zN3L += v8*(params->S2dot6QMS2OAvg * LNhdotS2 * LNhcS2[2] );

	      *dS1x+= dS1xN3L;
	      *dS1y+= dS1yN3L;
	      *dS1z+= dS1zN3L;

	      *dS2x+= dS2xN3L;
	      *dS2y+= dS2yN3L;
	      *dS2z+= dS2zN3L;

	      dLNhatx -= dS1xN3L + dS2xN3L; // v8 terms
	      dLNhaty -= dS1yN3L + dS2yN3L;
	      dLNhatz -= dS1zN3L + dS2zN3L;

	      if (lscorr) {
		/* Contributions to dLNh coming from the S-dependent
		 * part of L, in particular: dS, dS_l at v6 order.
		 * We omit here terms parallel to LNhat 
		 * as they will be projected out at the end.
		 */
		dLNhatx-=eta*v2*( cS1*dS1xNL + cS2*dS2xNL + (cS1L*LNhdotS1+cS2L*LNhdotS2)*dLxL/LN0mag );  // v8 terms
		dLNhaty-=eta*v2*( cS1*dS1yNL + cS2*dS2yNL + (cS1L*LNhdotS1+cS2L*LNhdotS2)*dLyL/LN0mag );
		dLNhatz-=eta*v2*( cS1*dS1zNL + cS2*dS2zNL + (cS1L*LNhdotS1+cS2L*LNhdotS2)*dLzL/LN0mag );

	      }

	  }
	  /* spinO=7 in the orbit averaged case should not be used altogether, it is kept here
             for use of phenompt approximat for backward compatibility */
	  else {
	    if ( (params->spinO>=7) || (params->spinO<0) ) {

	      const REAL8 L2PN=XLALSimInspiralL_4PN(eta);
	      const REAL8 v4=v2*v2;
	      LNmag+=LN0mag*v4*L2PN;
	      const REAL8 omega3=omega2*omega;

	      // dS1,2 at NNNNLO, eq. 3.4 of Class. Quantum Grav. 30 (2013) 075017, arXiv/1212.5520
	      const REAL8 dS1xN4L = params->S1dot7S2 * omega3 * LNhcS1[0]; //v9 terms
	      const REAL8 dS1yN4L = params->S1dot7S2 * omega3 * LNhcS1[1];
	      const REAL8 dS1zN4L = params->S1dot7S2 * omega3 * LNhcS1[2];

	      const REAL8 dS2xN4L = params->S2dot7S1 * omega3 * LNhcS2[0];
	      const REAL8 dS2yN4L = params->S2dot7S1 * omega3 * LNhcS2[1];
	      const REAL8 dS2zN4L = params->S2dot7S1 * omega3 * LNhcS2[2];

	      *dS1x+= dS1xN4L;
	      *dS1y+= dS1yN4L;
	      *dS1z+= dS1zN4L;

	      *dS2x+= dS2xN4L;
	      *dS2y+= dS2yN4L;
	      *dS2z+= dS2zN4L;

	      dLNhatx -= dS1xN4L + dS2xN4L;
	      dLNhaty -= dS1yN4L + dS2yN4L;
	      dLNhatz -= dS1zN4L + dS2zN4L;

	    }
	  }
	}
      }
      XLALFree(S1cS2);
    }
    if (LNhcS1)
      XLALFree(LNhcS1);
    if (LNhcS2)
      XLALFree(LNhcS2);
  }

  /* We have computed the derivative of the spin-independent part of the
   * orbital angular momentum, which is parallel to the Newtonian angular momentum,
   * we now have to find the precession vector Om.
   */
  // First we normalize dLNhat
  dLNhatx/=LNmag;
  dLNhaty/=LNmag;
  dLNhatz/=LNmag;
  /* Then we compute the precession vector Om=LNhat x dLNhat
   * Note that component of dLNhat parallel to LNhat does not affect Om */
  REAL8 *Om=NULL;
  XLALSimInspiralVectorCrossProduct(&Om,LNhx,LNhy,LNhz,dLNhatx,dLNhaty,dLNhatz);
  /* Take cross product of Om with LNhat */
  *dLNhx = -Om[2]*LNhy + Om[1]*LNhz;
  *dLNhy = -Om[0]*LNhz + Om[2]*LNhx;
  *dLNhz = -Om[1]*LNhx + Om[0]*LNhy;

  /* Take cross product of Om with E_1 */
  *dE1x = -Om[2]*E1y + Om[1]*E1z;
  *dE1y = -Om[0]*E1z + Om[2]*E1x;
  *dE1z = -Om[1]*E1x + Om[0]*E1y;

  XLALFree(Om);Om=NULL;

  return XLAL_SUCCESS;

} /* End of XLALSimInspiralSpinDerivativesAvg() */

/**
 * See arXiv:0907.0700 for TaylorT4 definition.
 */
INT4 XLALSimInspiralSpinTaylorT4Setup(
    XLALSimInspiralSpinTaylorTxCoeffs **params, /**< OUTPUT */
    const REAL8 m1_SI,                    /**< mass of body 1 (kg) */
    const REAL8 m2_SI,                    /**< mass of body 2 (kg) */
    const REAL8 fStart,                   /**< Starting GW freq. (Hz) */
    const REAL8 fEnd,                     /**< Ending GW freq. (Hz), 0 means integrate forwards as far as possible */
    const REAL8 lambda1,                  /**< (tidal deformability of mass 1) / (mass of body 1)^5 (dimensionless) */
    const REAL8 lambda2,                  /**< (tidal deformability of mass 2) / (mass of body 2)^5 (dimensionless) */
    const REAL8 quadparam1,               /**< phenom. parameter describing induced quad. moment of body 1 (=1 for BHs, ~2-12 for NSs) */
    const REAL8 quadparam2,               /**< phenom. parameter describing induced quad. moment of body 2 (=1 for BHs, ~2-12 for NSs) */
    const LALSimInspiralSpinOrder spinO,  /**< twice PN order of spin effects */
    const LALSimInspiralTidalOrder tideO, /**< twice PN order of tidal effects */
    const INT4 phaseO,                    /**< twice post-Newtonian order */
    const INT4 lscorr,                    /**< flag to include spin corrections to orbital angular momentum */
    const INT4 phenomtp                   /**< flag for using spinO=7 and not spinO=6 terms with orbital-averaged quantities for phenomtphm approx */
    )
{
    *params = LALMalloc( sizeof(XLALSimInspiralSpinTaylorTxCoeffs));
    memset(*params, 0, sizeof(XLALSimInspiralSpinTaylorTxCoeffs) );
    int errCode=XLALSimSpinTaylorEnergySpinDerivativeSetup(params,m1_SI,m2_SI,fStart,fEnd,spinO,tideO,phaseO,lambda1,lambda2,quadparam1,quadparam2,lscorr,phenomtp);

    if ( (lscorr) && (phenomtp) ) {
      XLALPrintError("XLAL Error - %s: Flag for spin corrections to orbital angular momentum is incompatible with IMRPhenomTP setup\n",__func__);
      XLAL_ERROR(XLAL_EINVAL);
    }

    REAL8 eta=(*params)->eta;
    REAL8 m1M=(*params)->m1M;
    REAL8 m2M=(*params)->m2M;

    /* Set omegadot coefficients up to PN order phaseO.
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
	  (*params)->wdotcoeff[7] = XLALSimInspiralTaylorT4wdot_7PNCoeff(eta);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        /* case LAL_PNORDER_THREE: */
        case 6:
	  (*params)->wdotcoeff[6] = XLALSimInspiralTaylorT4wdot_6PNCoeff(eta);
	  (*params)->wdotlogcoeff = XLALSimInspiralTaylorT4wdot_6PNLogCoeff(eta);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        /* case LAL_PNORDER_TWO_POINT_FIVE: */
        case 5:
            (*params)->wdotcoeff[5] = XLALSimInspiralTaylorT4wdot_5PNCoeff(eta);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        /* case LAL_PNORDER_TWO: */
        case 4:
            (*params)->wdotcoeff[4] = XLALSimInspiralTaylorT4wdot_4PNCoeff(eta);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        /*case LAL_PNORDER_ONE_POINT_FIVE:*/
        case 3:
            (*params)->wdotcoeff[3] = XLALSimInspiralTaylorT4wdot_3PNCoeff(eta);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        /*case LAL_PNORDER_ONE:*/
        case 2:
            (*params)->wdotcoeff[2] = XLALSimInspiralTaylorT4wdot_2PNCoeff(eta);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        /*case LAL_PNORDER_HALF:*/
        case 1:
            (*params)->wdotcoeff[1] = 0.;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        /*case LAL_PNORDER_NEWTONIAN:*/
        case 0:
            (*params)->wdotcoeff[0] = 1.;
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
	  if ( phenomtp ) {
	    (*params)->wdot7S1O = XLALSimInspiralTaylorT4wdot_7PNSOCoeff(m1M);
            (*params)->wdot7S2O = XLALSimInspiralTaylorT4wdot_7PNSOCoeff(m2M);
	  }
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
            (*params)->wdot6S1O      = XLALSimInspiralTaylorT4wdot_6PNSOCoeff(m1M);
            (*params)->wdot6S2O      = XLALSimInspiralTaylorT4wdot_6PNSOCoeff(m2M);
	    if (!(phenomtp)) {
	      (*params)->wdot6S1S2Avg     = XLALSimInspiralTaylorT4wdot_6PNS1S2CoeffAvg(eta);
	      (*params)->wdot6S1OS2OAvg   = XLALSimInspiralTaylorT4wdot_6PNS1OS2OCoeffAvg(eta);
	      (*params)->wdot6S1S1Avg     = XLALSimInspiralTaylorT4wdot_6PNS1S1CoeffAvg(m1M);
	      (*params)->wdot6S1OS1OAvg   = XLALSimInspiralTaylorT4wdot_6PNS1OS1OCoeffAvg(m1M);
	      (*params)->wdot6S2S2Avg     = XLALSimInspiralTaylorT4wdot_6PNS1S1CoeffAvg(m2M);
	      (*params)->wdot6S2OS2OAvg   = XLALSimInspiralTaylorT4wdot_6PNS1OS1OCoeffAvg(m2M);
	      (*params)->wdot6QMS1S1Avg   = quadparam1 * XLALSimInspiralTaylorT4wdot_6PNQMS1S1CoeffAvg(m1M);
	      (*params)->wdot6QMS1OS1OAvg = quadparam1 * XLALSimInspiralTaylorT4wdot_6PNQMS1OS1OCoeffAvg(m1M);
	      (*params)->wdot6QMS2S2Avg   = quadparam2 * XLALSimInspiralTaylorT4wdot_6PNQMS1S1CoeffAvg(m2M);
	      (*params)->wdot6QMS2OS2OAvg = quadparam2 * XLALSimInspiralTaylorT4wdot_6PNQMS1OS1OCoeffAvg(m2M);
	    }
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
            (*params)->wdot5S1O = XLALSimInspiralTaylorT4wdot_5PNSOCoeff(m1M);
            (*params)->wdot5S2O = XLALSimInspiralTaylorT4wdot_5PNSOCoeff(m2M);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
        /* Averaging terms bilinear in terms periodic in an orbit and
         * confined to the orbital plane introduces v^3
         * inaccuracies, see e.g. app. B of PRD80 (2009) 044010, arXiv:0812.4413.
         * As a consequence orbital averaging 2PN spin^2 terms
         * impact 3.5PN terms, hence the latter should only be used in the
         * spin-aligned case only when orbit average is taken.
         */
	  // Orbit-averaged terms
            // 2PN spin-spin terms
            (*params)->wdot4S1S2Avg  = XLALSimInspiralTaylorT4wdot_4PNS1S2CoeffAvg(eta);
            (*params)->wdot4S1OS2OAvg= XLALSimInspiralTaylorT4wdot_4PNS1OS2OCoeffAvg(eta);
            // 2PN self-spin terms
            (*params)->wdot4S1S1Avg  = XLALSimInspiralTaylorT4wdot_4PNS1S1CoeffAvg(m1M);
            (*params)->wdot4S1OS1OAvg= XLALSimInspiralTaylorT4wdot_4PNS1OS1OCoeffAvg(m1M);
            (*params)->wdot4S2S2Avg  = XLALSimInspiralTaylorT4wdot_4PNS1S1CoeffAvg(m2M);
            (*params)->wdot4S2OS2OAvg= XLALSimInspiralTaylorT4wdot_4PNS1OS1OCoeffAvg(m2M);
            // 2PN quadrupole-monopole terms
            (*params)->wdot4QMS1S1Avg  = quadparam1 * XLALSimInspiralTaylorT4wdot_4PNQMS1S1CoeffAvg(m1M);
            (*params)->wdot4QMS1OS1OAvg= quadparam1 * XLALSimInspiralTaylorT4wdot_4PNQMS1OS1OCoeffAvg(m1M);
            (*params)->wdot4QMS2S2Avg  = quadparam2 * XLALSimInspiralTaylorT4wdot_4PNQMS1S1CoeffAvg(m2M);
            (*params)->wdot4QMS2OS2OAvg= quadparam2 * XLALSimInspiralTaylorT4wdot_4PNQMS1OS1OCoeffAvg(m2M);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
	    (*params)->wdot3S1O 	= XLALSimInspiralTaylorT4wdot_3PNSOCoeff(m1M);
            (*params)->wdot3S2O 	= XLALSimInspiralTaylorT4wdot_3PNSOCoeff(m2M);
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
            (*params)->wdottidal12 = lambda1 * XLALSimInspiralTaylorT4wdot_12PNTidalCoeff(m1M) + lambda2 * XLALSimInspiralTaylorT4wdot_12PNTidalCoeff(m2M);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
	    (*params)->wdottidal10 = lambda1 * XLALSimInspiralTaylorT4wdot_10PNTidalCoeff(m1M) + lambda2 * XLALSimInspiralTaylorT4wdot_10PNTidalCoeff(m2M);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid tidal PN order %d\n",
                    __func__, tideO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    return errCode;
} // End of XLALSimInspiralSpinTaylorT4Setup()

/**
 * See arXiv:0907.0700 for TaylorT1 definition.
 */
int XLALSimInspiralSpinTaylorT1Setup(
    XLALSimInspiralSpinTaylorTxCoeffs **params, /**< OUTPUT */
    const REAL8 m1_SI,                    /**< mass of body 1 (kg) */
    const REAL8 m2_SI,                    /**< mass of body 2 (kg) */
    const REAL8 fStart,                   /**< Starting GW freq. (Hz) */
    const REAL8 fEnd,                     /**< Ending GW freq. (Hz), 0 means integrate forwards as far as possible */
    const REAL8 lambda1,                  /**< (tidal deformability of mass 1) / (mass of body 1)^5 (dimensionless) */
    const REAL8 lambda2,                  /**< (tidal deformability of mass 2) / (mass of body 2)^5 (dimensionless) */
    const REAL8 quadparam1,               /**< phenom. parameter describing induced quad. moment of body 1 (=1 for BHs, ~2-12 for NSs) */
    const REAL8 quadparam2,               /**< phenom. parameter describing induced quad. moment of body 2 (=1 for BHs, ~2-12 for NSs) */
    const LALSimInspiralSpinOrder spinO,  /**< twice PN order of spin effects */
    const LALSimInspiralTidalOrder tideO, /**< twice PN order of tidal effects */
    const INT4 phaseO,                    /**< twice post-Newtonian order */
    const INT4 lscorr                     /**< flag to include spin corrections to orbital angular momentum */
    )
{
    *params = LALMalloc( sizeof(XLALSimInspiralSpinTaylorTxCoeffs));
    memset(*params, 0, sizeof(XLALSimInspiralSpinTaylorTxCoeffs) );
    int errCode=XLALSimSpinTaylorEnergySpinDerivativeSetup(params,m1_SI,m2_SI,fStart,fEnd,spinO,tideO,phaseO,lambda1,lambda2,quadparam1,quadparam2,lscorr,0);

    REAL8 eta=(*params)->eta;
    REAL8 m1M=(*params)->m1M;
    REAL8 m2M=(*params)->m2M;
    (*params)->Fnewt=XLALSimInspiralPNFlux_0PNCoeff(eta);
    (*params)->dEdvnewt=2.*XLALSimInspiralPNEnergy_0PNCoeff(eta);

    /* Set the flux coefficients up to PN order phaseO.
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
	  (*params)->Fcoeff[7] = XLALSimInspiralPNFlux_7PNCoeff(eta);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        /* case LAL_PNORDER_THREE: */
        case 6:
	  (*params)->Fcoeff[6] = XLALSimInspiralPNFlux_6PNCoeff(eta);
	  (*params)->Flogcoeff = XLALSimInspiralPNFlux_6PNLogCoeff(eta);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        /* case LAL_PNORDER_TWO_POINT_FIVE: */
        case 5:
	  (*params)->Fcoeff[5] = XLALSimInspiralPNFlux_5PNCoeff(eta);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        /* case LAL_PNORDER_TWO: */
        case 4:
	  (*params)->Fcoeff[4] = XLALSimInspiralPNFlux_4PNCoeff(eta);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        /*case LAL_PNORDER_ONE_POINT_FIVE:*/
        case 3:
	  (*params)->Fcoeff[3] = XLALSimInspiralPNFlux_3PNCoeff(eta);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        /*case LAL_PNORDER_ONE:*/
        case 2:
	  (*params)->Fcoeff[2] = XLALSimInspiralPNFlux_2PNCoeff(eta);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        /*case LAL_PNORDER_HALF:*/
        case 1:
	  (*params)->Fcoeff[1] = 0.;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        /*case LAL_PNORDER_NEWTONIAN:*/
        case 0:
	  (*params)->Fcoeff[0] = 1.;
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
	  (*params)->F7S1O=XLALSimInspiralPNFlux_7PNSOCoeff(m1M);
	  (*params)->F7S2O=XLALSimInspiralPNFlux_7PNSOCoeff(m2M);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
	  (*params)->F6S1O    = XLALSimInspiralPNFlux_6PNSOCoeff(m1M);
	  (*params)->F6S2O    = XLALSimInspiralPNFlux_6PNSOCoeff(m2M);
	  (*params)->F6S1S2Avg   = XLALSimInspiralPNFlux_6PNS1S2CoeffAvg(eta);
	  (*params)->F6S1OS2OAvg = XLALSimInspiralPNFlux_6PNS1OS2OCoeffAvg(eta);
	  (*params)->F6S1S1Avg   = XLALSimInspiralPNFlux_6PNS1S1CoeffAvg(m1M);
	  (*params)->F6S1OS1OAvg = XLALSimInspiralPNFlux_6PNS1OS1OCoeffAvg(m1M);
	  (*params)->F6S2S2Avg   = XLALSimInspiralPNFlux_6PNS1S1CoeffAvg(m2M);
	  (*params)->F6S2OS2OAvg = XLALSimInspiralPNFlux_6PNS1OS1OCoeffAvg(m2M);
	  (*params)->F6QMS1S1Avg   = quadparam1*XLALSimInspiralPNFlux_6PNQMS1S1CoeffAvg(m1M);
	  (*params)->F6QMS1OS1OAvg = quadparam1*XLALSimInspiralPNFlux_6PNQMS1OS1OCoeffAvg(m1M);
	  (*params)->F6QMS2S2Avg   = quadparam2*XLALSimInspiralPNFlux_6PNQMS1S1CoeffAvg(m2M);
	  (*params)->F6QMS2OS2OAvg = quadparam2*XLALSimInspiralPNFlux_6PNQMS1OS1OCoeffAvg(m2M);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
       /* Averaging terms bilinear in terms periodic in an orbit and
        * confined to the orbital plane introduces v^3
        * inaccuracies, see e.g. app. B of PRD80 (2009) 044010, arXiv:0812.4413.
        * As a consequence orbital averaging 2PN spin^2 terms
        * impact 3.5PN terms, hence the latter should only be used in the
        * spin-aligned case only when orbit average is taken.
        */
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
	  (*params)->F5S1O =XLALSimInspiralPNFlux_5PNSOCoeff(m1M);
	  (*params)->F5S2O =XLALSimInspiralPNFlux_5PNSOCoeff(m2M);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
	  /* Orbit Averaged */
	  // 2PN spin-spin terms
	  (*params)->F4S1S2Avg   = XLALSimInspiralPNFlux_4PNS1S2CoeffAvg(eta);
	  (*params)->F4S1OS2OAvg = XLALSimInspiralPNFlux_4PNS1OS2OCoeffAvg(eta);
	  // 2PN self-spin terms
	  (*params)->F4S1S1Avg   = XLALSimInspiralPNFlux_4PNS1S1CoeffAvg(m1M);
	  (*params)->F4S1OS1OAvg = XLALSimInspiralPNFlux_4PNS1OS1OCoeffAvg(m1M);
	  (*params)->F4S2S2Avg   = XLALSimInspiralPNFlux_4PNS1S1CoeffAvg(m2M);
	  (*params)->F4S2OS2OAvg = XLALSimInspiralPNFlux_4PNS1OS1OCoeffAvg(m2M);
	  // 2PN quadrupole-monopole terms
	  (*params)->F4QMS1S1Avg   = quadparam1*XLALSimInspiralPNFlux_4PNQMS1S1CoeffAvg(m1M);
	  (*params)->F4QMS1OS1OAvg = quadparam1*XLALSimInspiralPNFlux_4PNQMS1OS1OCoeffAvg(m1M);
	  (*params)->F4QMS2S2Avg   = quadparam2*XLALSimInspiralPNFlux_4PNQMS1S1CoeffAvg(m2M);
	  (*params)->F4QMS2OS2OAvg = quadparam2*XLALSimInspiralPNFlux_4PNQMS1OS1OCoeffAvg(m2M);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
	  (*params)->F3S1O      = XLALSimInspiralPNFlux_3PNSOCoeff(m1M);
	  (*params)->F3S2O      = XLALSimInspiralPNFlux_3PNSOCoeff(m2M);
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

    /* Compute the coefficients of tidal corrections
     * to the evolution equations for omega and binary energy E.
     * Coefficients found from Eqs. 2.11 and 3.10 of
     * Vines, Flanagan, Hinderer, PRD 83, 084051 (2011).
     */
    switch( tideO )
    {
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL:
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN:
            (*params)->Ftidal12 = lambda1 * XLALSimInspiralPNFlux_12PNTidalCoeff(m1M) + lambda2 * XLALSimInspiralPNFlux_12PNTidalCoeff(m2M);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
	    (*params)->Ftidal10 = lambda1 * XLALSimInspiralPNFlux_10PNTidalCoeff(m1M) + lambda2 * XLALSimInspiralPNFlux_10PNTidalCoeff(m2M);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid tidal PN order %d\n",
                    __func__, tideO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    return errCode;
} //End of XLALSimInspiralSpinTaylorT1Setup()

/**
 * See arXiv:1107.1267 for TaylorT5 approximant definition. It is a variant
 * of TaylorT2 (see arXiv:0907.0700). SpinTaylorT5 as implemented in this
 * code was previously (earlier than June 2019) named SpinTaylorT2.
 */
int XLALSimInspiralSpinTaylorT5Setup(
    XLALSimInspiralSpinTaylorTxCoeffs **params, /**< OUTPUT */
    const REAL8 m1_SI,                    /**< mass of body 1 (kg) */
    const REAL8 m2_SI,                    /**< mass of body 2 (kg) */
    const REAL8 fStart,                   /**< Starting GW freq. (Hz) */
    const REAL8 fEnd,                     /**< Ending GW freq. (Hz), 0 means integrate forwards as far as possible */
    const REAL8 lambda1,                  /**< (tidal deformability of mass 1) / (mass of body 1)^5 (dimensionless) */
    const REAL8 lambda2,                  /**< (tidal deformability of mass 2) / (mass of body 2)^5 (dimensionless) */
    const REAL8 quadparam1,               /**< phenom. parameter describing induced quad. moment of body 1 (=1 for BHs, ~2-12 for NSs) */
    const REAL8 quadparam2,               /**< phenom. parameter describing induced quad. moment of body 2 (=1 for BHs, ~2-12 for NSs) */
    const LALSimInspiralSpinOrder spinO,  /**< twice PN order of spin effects */
    const LALSimInspiralTidalOrder tideO, /**< twice PN order of tidal effects */
    const INT4 phaseO,                    /**< twice post-Newtonian order */
    const INT4 lscorr                    /**< flag to include spin corrections to orbital angular momentum */
				     )
{
    *params = LALMalloc( sizeof(XLALSimInspiralSpinTaylorTxCoeffs));
    memset(*params, 0, sizeof(XLALSimInspiralSpinTaylorTxCoeffs) );
    int errCode=XLALSimSpinTaylorEnergySpinDerivativeSetup(params,m1_SI,m2_SI,fStart,fEnd,spinO,tideO,phaseO,lambda1,lambda2,quadparam1,quadparam2,lscorr,0);
    REAL8 eta=(*params)->eta;
    REAL8 m1M=(*params)->m1M;
    REAL8 m2M=(*params)->m2M;

    /* Set wdot coefficients up to PN order phaseO.
     * wdot is the derivative of the orbital frequency \f$\dot{\omega}\f$.
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
	    (*params)->wdotcoeff[7] = XLALSimInspiralTaylorT2dtdv_7PNCoeff(eta);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        /* case LAL_PNORDER_THREE: */
        case 6:
            (*params)->wdotcoeff[6] = XLALSimInspiralTaylorT2dtdv_6PNCoeff(eta);
            (*params)->wdotlogcoeff = XLALSimInspiralTaylorT2dtdv_6PNLogCoeff(eta);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        /* case LAL_PNORDER_TWO_POINT_FIVE: */
        case 5:
            (*params)->wdotcoeff[5] = XLALSimInspiralTaylorT2dtdv_5PNCoeff(eta);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        /* case LAL_PNORDER_TWO: */
        case 4:
            (*params)->wdotcoeff[4] = XLALSimInspiralTaylorT2dtdv_4PNCoeff(eta);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        /*case LAL_PNORDER_ONE_POINT_FIVE:*/
        case 3:
            (*params)->wdotcoeff[3] = XLALSimInspiralTaylorT2dtdv_3PNCoeff(eta);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        /*case LAL_PNORDER_ONE:*/
        case 2:
            (*params)->wdotcoeff[2] = XLALSimInspiralTaylorT2dtdv_2PNCoeff(eta);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        /*case LAL_PNORDER_HALF:*/
        case 1:
            (*params)->wdotcoeff[1] = 0.;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        /*case LAL_PNORDER_NEWTONIAN:*/
        case 0:
            (*params)->wdotcoeff[0] = 1.;
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
    {   case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
            (*params)->wdot7S1O = XLALSimInspiralTaylorT2dtdv_7PNSOCoeff(m1M);
            (*params)->wdot7S2O = XLALSimInspiralTaylorT2dtdv_7PNSOCoeff(m2M);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
	  /* Averaging terms bilinear in terms periodic in an orbit and
	   * confined to the orbital plane introduces v^3
	   * inaccuracies, see e.g. app. B of PRD80 (2009) 044010, arXiv:0812.4413.
	   * As a consequence orbital averaging 2PN spin^2 terms
	   * impact 3.5PN terms, hence the latter should only be used in the
	   * spin-aligned case only when orbit average is taken.
	   */
	    (*params)->wdot6S1O       = XLALSimInspiralTaylorT2dtdv_6PNSOCoeff(m1M);
	    (*params)->wdot6S2O       = XLALSimInspiralTaylorT2dtdv_6PNSOCoeff(m2M);
            (*params)->wdot6S1S2Avg   = XLALSimInspiralTaylorT2dtdv_6PNS1S2CoeffAvg(eta);
            (*params)->wdot6S1OS2OAvg = XLALSimInspiralTaylorT2dtdv_6PNS1OS2OCoeffAvg(eta);
            (*params)->wdot6S1S1Avg   = XLALSimInspiralTaylorT2dtdv_6PNS1S1CoeffAvg(m1M);
            (*params)->wdot6S1OS1OAvg = XLALSimInspiralTaylorT2dtdv_6PNS1OS1OCoeffAvg(m1M);
            (*params)->wdot6S2S2Avg   = XLALSimInspiralTaylorT2dtdv_6PNS1S1CoeffAvg(m2M);
            (*params)->wdot6S2OS2OAvg = XLALSimInspiralTaylorT2dtdv_6PNS1OS1OCoeffAvg(m2M);
            (*params)->wdot6QMS1S1Avg = quadparam1*XLALSimInspiralTaylorT2dtdv_6PNQMS1S1CoeffAvg(m1M);
            (*params)->wdot6QMS2S2Avg = quadparam2*XLALSimInspiralTaylorT2dtdv_6PNQMS1S1CoeffAvg(m2M);
	    (*params)->wdot6QMS1OS1OAvg  = quadparam1*XLALSimInspiralTaylorT2dtdv_6PNQMS1OS1OCoeffAvg(m1M);
	    (*params)->wdot6QMS2OS2OAvg  = quadparam2*XLALSimInspiralTaylorT2dtdv_6PNQMS1OS1OCoeffAvg(m2M);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
            (*params)->wdot5S1O = XLALSimInspiralTaylorT2dtdv_5PNSOCoeff(m1M);
            (*params)->wdot5S2O = XLALSimInspiralTaylorT2dtdv_5PNSOCoeff(m2M);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
	  // Orbit-averaged terms
            // 2PN spin1-spin2 terms
            (*params)->wdot4S1S2Avg  = XLALSimInspiralTaylorT2dtdv_4PNS1S2CoeffAvg(eta);
            (*params)->wdot4S1OS2OAvg= XLALSimInspiralTaylorT2dtdv_4PNS1OS2OCoeffAvg(eta);
            // 2PN spin-self^2 terms
            (*params)->wdot4S1S1Avg  = XLALSimInspiralTaylorT2dtdv_4PNS1S1CoeffAvg(m1M);
            (*params)->wdot4S1OS1OAvg= XLALSimInspiralTaylorT2dtdv_4PNS1OS1OCoeffAvg(m1M);
            (*params)->wdot4S2S2Avg  = XLALSimInspiralTaylorT2dtdv_4PNS1S1CoeffAvg(m2M);
            (*params)->wdot4S2OS2OAvg= XLALSimInspiralTaylorT2dtdv_4PNS1OS1OCoeffAvg(m2M);
            // 2PN quadrupole-monopole self spin terms
            (*params)->wdot4QMS1S1Avg  = quadparam1 * XLALSimInspiralTaylorT2dtdv_4PNQMS1S1CoeffAvg(m1M);
            (*params)->wdot4QMS1OS1OAvg= quadparam1 * XLALSimInspiralTaylorT2dtdv_4PNQMS1OS1OCoeffAvg(m1M);
            (*params)->wdot4QMS2S2Avg  = quadparam2 * XLALSimInspiralTaylorT2dtdv_4PNQMS1S1CoeffAvg(m2M);
            (*params)->wdot4QMS2OS2OAvg= quadparam2 * XLALSimInspiralTaylorT2dtdv_4PNQMS1OS1OCoeffAvg(m2M);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            // Note: LNHat do not have their signs reversed relative to T4
            // They are precession rather than orbital quantities
	    (*params)->wdot3S1O 	= XLALSimInspiralTaylorT2dtdv_3PNSOCoeff(m1M);
            (*params)->wdot3S2O 	= XLALSimInspiralTaylorT2dtdv_3PNSOCoeff(m2M);
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
	    (*params)->wdottidal12 = lambda1 * XLALSimInspiralTaylorT2dtdv_12PNTidalCoeff(m1M)
	                         + lambda2 * XLALSimInspiralTaylorT2dtdv_12PNTidalCoeff(m2M);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
	    (*params)->wdottidal10 = lambda1 * XLALSimInspiralTaylorT2dtdv_10PNTidalCoeff(m1M)
	                         + lambda2 * XLALSimInspiralTaylorT2dtdv_10PNTidalCoeff(m2M);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid tidal PN order %d\n",
                    __func__, tideO );
            XLAL_ERROR(XLAL_EINVAL);
            break;
    }

    return errCode;
} // End of XLALSimInspiralSpinTaylorT5Setup()

/*
 * Internal driver function to generate any of SpinTaylorT1/T5/T4 in the Fourier domain.
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
	LALDict *LALparams,             /**< LAL dictionary containing accessory parameters */
        INT4 phaseO,                    /**< twice PN phase order */
        INT4 amplitudeO,                /**< twice PN amplitude order */
        Approximant approx,             /**< PN approximant (SpinTaylorT1/T5/T4) */
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

    int spinO=XLALSimInspiralWaveformParamsLookupPNSpinOrder(LALparams);
    int tideO=XLALSimInspiralWaveformParamsLookupPNTidalOrder(LALparams);

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

    // The TaylorT5 version sometimes returns duplicate time steps (i.e. two different steps happen at the same time, because df/dt is too large near the end)
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

INT4 XLALSimInspiralSetEnergyPNTermsAvg(REAL8 *Espin3,
					REAL8 *Espin4,
					REAL8 *Espin5,
					REAL8 *Espin6,
					REAL8 *Espin7,
					XLALSimInspiralSpinTaylorTxCoeffs *params,
					const REAL8 LNhdotS1,
					const REAL8 LNhdotS2,
					const REAL8 S1sq,
					const REAL8 S2sq,
					const REAL8 S1dotS2)
{
  *Espin3=0.;
  *Espin4=0.;
  *Espin5=0.;
  *Espin6=0.;
  *Espin7=0.;
  switch( params->spinO ) {
    case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
      if (params->phenomtp)
	*Espin7 += params->E7S1O * LNhdotS1 + params->E7S2O * LNhdotS2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
    case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
    case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
      if (!(params->phenomtp))
	*Espin6 += params->E6S1S2Avg*S1dotS2 + params->E6S1OS2OAvg*LNhdotS1*LNhdotS2
	  + (params->E6S1S1Avg + params->E6QMS1S1Avg)*S1sq
	  + (params->E6S2S2Avg + params->E6QMS2S2Avg)*S2sq
	  + (params->E6S1OS1OAvg + params->E6QMS1OS1OAvg)*LNhdotS1*LNhdotS1 +
	  + (params->E6S2OS2OAvg + params->E6QMS2OS2OAvg)*LNhdotS2*LNhdotS2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
  case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
      // Compute 2.5PN SO correction to energy
      // See Eq. 7.9 of gr-qc/0605140v4
      // Note that S_l/M^2 = (m1/M)^2 chi1 + (m2/M)^2 chi2
      // and Sigma_l/M^2 = (m2/M) chi2 - (m1/M) chi1
      *Espin5 += params->E5S1O * LNhdotS1 + params->E5S2O * LNhdotS2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
    case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
      // Compute S1-S2 spin-spin term
      *Espin4 += params->E4S1S2Avg * S1dotS2 + params->E4S1OS2OAvg * LNhdotS1 * LNhdotS2;
      // Compute 2PN quadrupole-monopole correction to energy
      // See last line of Eq. 6 of astro-ph/0504538
      // or 2nd and 3rd lines of Eq. (C4) in arXiv:0810.5336v3
      *Espin4 += params->E4QMS1S1Avg * S1sq	+ params->E4QMS2S2Avg * S2sq
		+ params->E4QMS1OS1OAvg * LNhdotS1 * LNhdotS1
		+ params->E4QMS2OS2OAvg * LNhdotS2 * LNhdotS2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
    case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
      // Compute 1.5PN SO correction to energy
      *Espin3 += params->E3S1O * LNhdotS1 + params->E3S2O * LNhdotS2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
      __attribute__ ((fallthrough));
#endif
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
    return XLAL_SUCCESS;
}  //End of XLALSimInspiralSetEnergyPNTermsAvg()

/*
 * Internal function called by the integration routine.
 * Stops the integration if
 * 1) The energy decreases with increasing orbital frequency
 * 2) The orbital frequency begins decreasing
 * 3) The orbital frequency becomes infinite
 * 4) The orbital frequency has gone outside the requested bounds
 * 5) The PN parameter v/c becomes >= 1
 * SpinTaylorT4 and SpinTaylorT5 both use this same stopping test
 */
static int XLALSimInspiralSpinTaylorStoppingTest(double UNUSED t,
						    const double values[],
						    double dvalues[],
						    void *mparams
						    )
{
    REAL8 omega, v, test, omegaStart, omegaEnd, ddomega;
    REAL8 LNhx, LNhy, LNhz, S1x, S1y, S1z, S2x, S2y, S2z;
    REAL8 LNhdotS1, LNhdotS2, S1dotS2, S1sq, S2sq;
    REAL8 Espin7=0.;
    REAL8 Espin6=0.;
    REAL8 Espin5=0.;
    REAL8 Espin4=0.;
    REAL8 Espin3=0.;
    XLALSimInspiralSpinTaylorTxCoeffs *params 
            = (XLALSimInspiralSpinTaylorTxCoeffs*) mparams;
    /* Spin-corrections to energy (including dynamical terms) */

    omega = values[1];
    v = cbrt(omega);
    LNhx = values[2]; LNhy  = values[3]; LNhz = values[4] ;
    S1x  = values[5]; S1y   = values[6]; S1z  = values[7] ;
    S2x  = values[8]; S2y   = values[9]; S2z  = values[10];
    LNhdotS1 = cdot(LNhx,LNhy,LNhz,S1x,S1y,S1z);
    LNhdotS2 = cdot(LNhx,LNhy,LNhz,S2x,S2y,S2z);
    S1sq = normsq(S1x,S1y,S1z);
    S2sq = normsq(S2x,S2y,S2z);
    S1dotS2 = cdot(S1x,S1y,S1z,S2x,S2y,S2z);

    /* omega = PI G M f_GW / c^3
     * Note params->M is in solar mass units */
    omegaStart = LAL_PI * params->M*LAL_MTSUN_SI * params->fStart;
    omegaEnd = LAL_PI * params->M*LAL_MTSUN_SI * params->fEnd;

    XLALSimInspiralSetEnergyPNTermsAvg(&Espin3,&Espin4,&Espin5,&Espin6,&Espin7,params,LNhdotS1,LNhdotS2,S1sq,S2sq,S1dotS2);

    /* We are testing if the orbital energy increases with \f$\omega\f$.
     * We should be losing energy to GW flux, so if E increases 
     * we stop integration because the dynamics are becoming unphysical. 
     * 'test' is the PN expansion of \f$dE/d\omega\f$ without the prefactor, 
     * i.e. \f$dE/d\omega = dE/dv * dv/d\omega = - (M^2*eta/6) * test\f$
     * Therefore, the energy is increasing with \f$\omega\f$ iff. test < 0.
     */
    test = 2. + v * v * ( 4. * params->Ecoeff[2] 
            + v * ( 5. * (params->Ecoeff[3] + Espin3)
            + v * ( 6. * (params->Ecoeff[4] + Espin4)
            + v * ( 7. * (params->Ecoeff[5] + Espin5)
	    + v * ( 8. * (params->Ecoeff[6] + Espin6)
            + v * ( 9. * (params->Ecoeff[7] + Espin7)
			+ v * v * v * ( 12. * params->Etidal10
			+ v * v * ( 14. * params->Etidal12 ) ) ) ) ) ) ) );
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
        return LALSIMINSPIRAL_ST_TEST_FREQBOUND;
    else if( fabs(omegaEnd) > LAL_REAL4_EPS && omegaEnd < omegaStart
                && omega < omegaEnd) /* freq. below bound */
        return LALSIMINSPIRAL_ST_TEST_FREQBOUND;
    else if (test < 0.0) /* energy test fails! */
        return LALSIMINSPIRAL_ST_TEST_ENERGY;
    else if (isnan(omega)) /* omega is nan! */
        return LALSIMINSPIRAL_ST_TEST_OMEGANAN;
    else if (v >= 1.) // v/c >= 1!
        return LALSIMINSPIRAL_ST_TEST_LARGEV;
    else if (ddomega <= 0.) // d^2omega/dt^2 <= 0!
        return LALSIMINSPIRAL_ST_TEST_OMEGADOUBLEDOT;
    else /* Step successful, continue integrating */
        return GSL_SUCCESS;
} // End of XLALSimInspiralSpinTaylorStoppingTest()

/*
 * Utility function to compute the 3PN contribution to \f$\omega\f$ of circular orbit
 * coming from precession.
 * There is due to the mismatch at 3PN between \f$\dot{\phi}\f$ and \f$\omega\f$,
 * as explained in arXiv:1507.00406, leading to
 * \f$\dot{\phi} = \omega-(\Omega\cdot\lambda)^2/(2\dot{\phi})\f$,
 * where \f$\Omega\f$ is the precession vector of \f$\hat{L_N}\f$:
 * \f$\Omega= \hat{L_N}\times \dot{\hat{L_N}}\f$, which after orbit average gives
 * \f$\dot{\phi}=\omega-\Omega^2/(4\dot{\phi})\simeq \omega-\Omega^2/(4\omega)\f$
 */
static REAL8 omegashift(REAL8 S1sq,    // Modulo square of S1
			REAL8 S2sq,    // Modulo square of S2
			REAL8 S1S2,    // Scalar product of S1 and S2
			REAL8 LNhS1,   // Scalar product of LNh and S1
			REAL8 LNhS2,   // Scalar product of LNh and S2
			REAL8 OmS1,    // Coefficient of S1 in Omega
			REAL8 OmS2     // Coefficient of S2 in Omega
			)
{
  return -0.25*(OmS1*OmS1*(S1sq-LNhS1*LNhS1)+OmS2*OmS2*(S2sq-LNhS2*LNhS2)+2.*OmS1*OmS2*(S1S2-LNhS1*LNhS2));
}

/*
 * Internal function called by the integration routine.
 * Given the values of all the dynamical variables 'values' at time 't',
 * This function computes their derivatives 'dvalues'
 * so the ODE integrator can take a step
 * All non-dynamical quantities (masses, etc.) are passed in \"mparams\"
 *
 * The derivatives for \f$\omega\f$ in the spin-less case, can be found
 * as Eqs. (1), (7) of gr-qc/0405090, see below for spin effects
 *
 */

INT4 XLALSimInspiralSpinTaylorT4DerivativesAvg(
	REAL8 UNUSED t,
	const REAL8 values[],
	REAL8 dvalues[],
	void *mparams
	)
{
    /* coordinates and derivatives */
    REAL8 LNhx, LNhy, LNhz, E1x, E1y, E1z, S1x, S1y, S1z, S2x, S2y, S2z;
    REAL8 omega, domega, dLNhx, dLNhy, dLNhz, dE1x, dE1y, dE1z;
    REAL8 dS1x, dS1y, dS1z, dS2x, dS2y, dS2z;
    /* auxiliary variables */
    REAL8 v, v2, v11;
    REAL8 LNhdotS1, LNhdotS2, S1dotS2, S1sq, S2sq;
    REAL8 wspin3 = 0., wspin4Avg = 0., wspin5 = 0., wspin6Avg=0.;
    XLALSimInspiralSpinTaylorTxCoeffs *params 
            = (XLALSimInspiralSpinTaylorTxCoeffs*) mparams;

    /* copy variables */
    // UNUSED!!: s    = values[0] ;
    omega   	= values[1] ;
    LNhx = values[2] ; LNhy    	= values[3] ; LNhz 	= values[4] ;
    S1x  = values[5] ; S1y     	= values[6] ; S1z 	= values[7] ;
    S2x  = values[8] ; S2y     	= values[9] ; S2z 	= values[10];
    E1x  = values[11]; E1y     	= values[12] ; E1z 	= values[13];

    if (omega <= 0.0) /* orbital frequency must be positive! */
    {
        return LALSIMINSPIRAL_ST_DERIVATIVE_OMEGANONPOS;
    }

    v = cbrt(omega);
    v2  = v * v;
    v11= omega*omega*omega*v2;

    LNhdotS1 = cdot(LNhx,LNhy,LNhz,S1x,S1y,S1z);
    LNhdotS2 = cdot(LNhx,LNhy,LNhz,S2x,S2y,S2z);
    S1dotS2 = cdot(S1x,S1y,S1z,S2x,S2y,S2z);
    S1sq = cdot(S1x,S1y,S1z,S1x,S1y,S1z);
    S2sq = cdot(S2x,S2y,S2z,S2x,S2y,S2z);

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
    {   case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
	    wspin6Avg = params->wdot6S1O * LNhdotS1 + params->wdot6S2O * LNhdotS2
	      + params->wdot6S1OS2OAvg*LNhdotS1*LNhdotS2 + params->wdot6S1S2Avg*S1dotS2
	      + (params->wdot6S1S1Avg+params->wdot6QMS1S1Avg)*S1sq + (params->wdot6S2S2Avg+params->wdot6QMS2S2Avg)*S2sq
	      + (params->wdot6S1OS1OAvg + params->wdot6QMS1OS1OAvg)*LNhdotS1*LNhdotS1
	      + (params->wdot6S2OS2OAvg + params->wdot6QMS2OS2OAvg)*LNhdotS2*LNhdotS2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
            // Compute 2.5PN SO correction to domega/dt
            // See Eq. 8.3 of gr-qc/0605140v4
            // Note that S_l/M^2 = (m1/M)^2 chi1 + (m2/M)^2 chi2
            // and Sigma_l/M^2 = (m2/M) chi2 - (m1/M) chi1
            wspin5 = params->wdot5S1O*LNhdotS1 + params->wdot5S2O*LNhdotS2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
	    wspin4Avg = params->wdot4S1S2Avg * S1dotS2 + params->wdot4S1OS2OAvg * LNhdotS1 * LNhdotS2;
            // Compute 2PN QM and self-spin corrections to domega/dt
            // This is equivalent to Eqs. 9c + 9d of astro-ph/0504538
            wspin4Avg += (params->wdot4S1S1Avg + params->wdot4QMS1S1Avg) * S1sq
	      + (params->wdot4S2S2Avg + params->wdot4QMS2S2Avg) * S2sq
	      + (params->wdot4S1OS1OAvg + params->wdot4QMS1OS1OAvg) * LNhdotS1 * LNhdotS1
	      + (params->wdot4S2OS2OAvg + params->wdot4QMS2OS2OAvg) * LNhdotS2 * LNhdotS2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            // Compute 1.5PN SO correction to domega/dt
            wspin3 = params->wdot3S1O*LNhdotS1 + params->wdot3S2O*LNhdotS2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
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
            + v * ( params->wdotcoeff[3] + wspin3
            + v * ( params->wdotcoeff[4] + wspin4Avg
            + v * ( params->wdotcoeff[5] + wspin5
            + v * ( params->wdotcoeff[6] + wspin6Avg
                    + params->wdotlogcoeff * log(v)
            + v * ( params->wdotcoeff[7]
            + omega * ( params->wdottidal10
			+ v2 * ( params->wdottidal12 ) ) ) ) ) ) ) ) ) );

    XLALSimInspiralSpinDerivativesAvg(&dLNhx,&dLNhy,&dLNhz,&dE1x,&dE1y,&dE1z,&dS1x,&dS1y,&dS1z,&dS2x,&dS2y,&dS2z,v,LNhx,LNhy,LNhz,E1x,E1y,E1z,S1x,S1y,S1z,S2x,S2y,S2z,LNhdotS1,LNhdotS2,params);

    /* dphi = d \phi / d \hat{t} = M d \phi /dt = M \omega = \hat{\omega} */
    dvalues[0]    = omega*(1.+omega*omega*omegashift(S1sq,S2sq,S1dotS2,LNhdotS1,LNhdotS2,params->omegashiftS1,params->omegashiftS2));
    dvalues[1]    = domega;
    dvalues[2]    = dLNhx; dvalues[3]     = dLNhy ; dvalues[4]    = dLNhz;
    dvalues[5]    = dS1x ; dvalues[6]     = dS1y  ; dvalues[7]    = dS1z ;
    dvalues[8]    = dS2x ; dvalues[9]     = dS2y  ; dvalues[10]   = dS2z ;
    dvalues[11]   = dE1x ; dvalues[12]    = dE1y  ; dvalues[13]   = dE1z ;

    return GSL_SUCCESS;
} //End of XLALSimInspiralSpinTaylorT4DerivativesAvg()

/*
 * Internal function called by the integration routine.
 * Given the values of all the dynamical variables 'values' at time 't',
 * This function computes their derivatives 'dvalues'
 * so the ODE integrator can take a step
 * All non-dynamical quantities (masses, etc.) are passed in \"mparams\"
 *
 * The derivatives for \f$\omega\f$ in the spin-less case, can be found
 * as Eqs. (1), (7) of gr-qc/0405090, see below for spin effects
 * 
 */

static int XLALSimInspiralSpinTaylorT1DerivativesAvg(
	double UNUSED t,
	const double values[],
	double dvalues[],
	void *mparams
	)
{
    /* coordinates and derivatives */
    REAL8 LNhx, LNhy, LNhz, S1x, S1y, S1z, S2x, S2y, S2z,E1x,E1y,E1z;
    REAL8 omega, domega, dLNhx, dLNhy, dLNhz, dE1x, dE1y, dE1z;
    REAL8 dS1x, dS1y, dS1z, dS2x, dS2y, dS2z;

    /* auxiliary variables */
    REAL8 v, v2, v3, v4, v7, v11;
    REAL8 LNhdotS1, LNhdotS2, S1dotS2, S1sq, S2sq;
    REAL8 Fspin3 = 0., Fspin4Avg = 0., Fspin5 = 0., Fspin6Avg = 0.;
    REAL8 Espin3 = 0., Espin4Avg = 0., Espin5 = 0., Espin6Avg = 0., Espin7 = 0.;

    XLALSimInspiralSpinTaylorTxCoeffs *params
            = (XLALSimInspiralSpinTaylorTxCoeffs*) mparams;

    /* copy variables */
    omega   	= values[1] ;
    LNhx = values[2] ; LNhy    	= values[3] ; LNhz 	= values[4] ;
    S1x  = values[5] ; S1y     	= values[6] ; S1z 	= values[7] ;
    S2x  = values[8] ; S2y     	= values[9] ; S2z 	= values[10];
    E1x  = values[11]; E1y     	= values[12]; E1z 	= values[13];

    if (omega <= 0.0) /* orbital frequency must be positive! */
    {
        return LALSIMINSPIRAL_ST_DERIVATIVE_OMEGANONPOS;
    }

    v = cbrt(omega);
    v2  = v * v; v3 = v2 * v; v4 = v3 * v;
    v7 = v4 * v3; v11 = v7 * v4;

    LNhdotS1 = (LNhx*S1x + LNhy*S1y + LNhz*S1z);
    LNhdotS2 = (LNhx*S2x + LNhy*S2y + LNhz*S2z);
    S1dotS2 = (S1x*S2x  + S1y*S2y  + S1z*S2z );
    S1sq = (S1x*S1x + S1y*S1y + S1z*S1z);
    S2sq = (S2x*S2x + S2y*S2y + S2z*S2z);
    S1dotS2 = (S1x*S2x + S1y*S2y + S1z*S2z);

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
    {   case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
	    Fspin6Avg = params->F6S1S2Avg*S1dotS2 + params->F6S1OS2OAvg * LNhdotS1*LNhdotS2;
            Fspin6Avg+= (params->F6S1S1Avg + params->F6QMS1S1Avg) * S1sq
	              + (params->F6S2S2Avg + params->F6QMS2S2Avg) * S2sq
		      + (params->F6S1OS1OAvg + params->F6QMS1OS1OAvg)* LNhdotS1 * LNhdotS1
		      + (params->F6S2OS2OAvg + params->F6QMS2OS2OAvg)* LNhdotS2 * LNhdotS2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
            // Compute 2.5PN SO correction to domega/dt
            // See Eq. 8.3 of gr-qc/0605140v4
            // Note that S_l/M^2 = (m1/M)^2 chi1 + (m2/M)^2 chi2
            // and Sigma_l/M^2 = (m2/M) chi2 - (m1/M) chi1
            Fspin5 += params->F5S1O*LNhdotS1 + params->F5S2O*LNhdotS2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
	    Fspin4Avg += params->F4S1S2Avg*S1dotS2 + params->F4S1OS2OAvg * LNhdotS1*LNhdotS2;
            // Compute 2PN QM and self-spin corrections to domega/dt
            Fspin4Avg += (params->F4S1S1Avg + params->F4QMS1S1Avg) * S1sq
			+ (params->F4S2S2Avg + params->F4QMS2S2Avg) * S2sq
			+ (params->F4S1OS1OAvg + params->F4QMS1OS1OAvg)* LNhdotS1 * LNhdotS1
			+ (params->F4S2OS2OAvg + params->F4QMS2OS2OAvg)* LNhdotS2 * LNhdotS2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            // Compute 1.5PN SO correction to domega/dt
            Fspin3 = params->F3S1O*LNhdotS1 + params->F3S2O*LNhdotS2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
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
    XLALSimInspiralSetEnergyPNTermsAvg(&Espin3,&Espin4Avg,&Espin5,&Espin6Avg,&Espin7,params,LNhdotS1,LNhdotS2,S1sq,S2sq,S1dotS2);

    XLALSimInspiralSpinDerivativesAvg(&dLNhx,&dLNhy,&dLNhz,&dE1x,&dE1y,&dE1z,&dS1x,&dS1y,&dS1z,&dS2x,&dS2y,&dS2z,v,LNhx,LNhy,LNhz,E1x,E1y,E1z,S1x,S1y,S1z,S2x,S2y,S2z,LNhdotS1,LNhdotS2,params);

    domega  =-6.0* (params->Fnewt/params->dEdvnewt) * v11 *
      (((( 1. + v * v * (  params->Fcoeff[2]
	    + v * (  (params->Fcoeff[3] + Fspin3 )
            + v * (  (params->Fcoeff[4] + Fspin4Avg )
            + v * (  (params->Fcoeff[5] + Fspin5)
            + v * (  (params->Fcoeff[6] + Fspin6Avg+ params->Flogcoeff*log(v))
            + v * (  (params->Fcoeff[7])
            + v * v * v * (  params->Ftidal10
                          + v * v * ( params->Ftidal12 ) ) )))))))))/
           ( 2. + v * v * ( 4. * params->Ecoeff[2]
            + v * (( 5. * (params->Ecoeff[3] + Espin3 ))
            + v * (( 6. * (params->Ecoeff[4] + Espin4Avg ))
            + v * (( 7. * (params->Ecoeff[5] + Espin5 ))
		   + v * (( 8. * (params->Ecoeff[6] + Espin6Avg ))
            + v * ( (9. * (params->Ecoeff[7] ) )
            + v * v * v * (( 12. * params->Etidal10)
                        + v * v * ( 14. * params->Etidal12 ) )
		    ))))))) );

    /* dphi = d \phi / d \hat{t} = M d \phi /dt = M \omega = \hat{\omega} */
    dvalues[0]    = omega*(1.+omega*omega*omegashift(S1sq,S2sq,S1dotS2,LNhdotS1,LNhdotS2,params->omegashiftS1,params->omegashiftS2));
    dvalues[1]    = domega;
    dvalues[2]    = dLNhx; dvalues[3]     = dLNhy ; dvalues[4]    = dLNhz;
    dvalues[5]    = dS1x ; dvalues[6]     = dS1y  ; dvalues[7]    = dS1z ;
    dvalues[8]    = dS2x ; dvalues[9]     = dS2y  ; dvalues[10]   = dS2z ;
    dvalues[11]   = dE1x ; dvalues[12]    = dE1y  ; dvalues[13]   = dE1z ;

    return GSL_SUCCESS;
} // End of XLALSimInspiralSpinTaylorT1DerivativesAvg()

/*
 * Internal function called by the integration routine.
 * Given the values of all the dynamical variables 'values' at time 't',
 * This function computes their derivatives 'dvalues'
 * so the ODE integrator can take a step
 * All non-dynamical quantities (masses, etc.) are passed in \"mparams\"
 *
 * The derivatives for \f$\omega\f$ in the spin-less case, can be found
 * as Eqs. (1), (7) of gr-qc/0405090, see below for spin effects
 * 
 */

static int XLALSimInspiralSpinTaylorT5DerivativesAvg(
	double UNUSED t,
	const double values[],
	double dvalues[],
	void *mparams
	)
{
    /* coordinates and derivatives */
    REAL8 LNhx, LNhy, LNhz, S1x, S1y, S1z, S2x, S2y, S2z, E1x, E1y, E1z;
    REAL8 omega, domega, dLNhx, dLNhy, dLNhz,dE1x,dE1y,dE1z;
    REAL8 dS1x, dS1y, dS1z, dS2x, dS2y, dS2z;
    
    /* auxiliary variables */
    REAL8 v,v11;
    REAL8 LNhdotS1, LNhdotS2, S1dotS2, S1sq, S2sq;
    REAL8 wspin3 = 0., wspin4Avg = 0., wspin5 = 0., wspin6Avg = 0.;
    
    XLALSimInspiralSpinTaylorTxCoeffs *params
            = (XLALSimInspiralSpinTaylorTxCoeffs*) mparams;

    /* copy variables */
    omega   	= values[1] ;
    LNhx = values[2] ; LNhy    	= values[3] ; LNhz 	= values[4] ;
    S1x  = values[5] ; S1y     	= values[6] ; S1z 	= values[7] ;
    S2x  = values[8] ; S2y     	= values[9] ; S2z 	= values[10];
    E1x  = values[11]; E1y     	= values[12]; E1z 	= values[13];

    if (omega <= 0.0) /* orbital frequency must be positive! */
    {
        return LALSIMINSPIRAL_ST_DERIVATIVE_OMEGANONPOS;
    }

    v = cbrt(omega);
    v11=omega*omega*omega*v*v;

    LNhdotS1 = (LNhx*S1x + LNhy*S1y + LNhz*S1z);
    LNhdotS2 = (LNhx*S2x + LNhy*S2y + LNhz*S2z);
    S1dotS2 = (S1x*S2x  + S1y*S2y  + S1z*S2z );
    S1sq = (S1x*S1x + S1y*S1y + S1z*S1z);
    S2sq = (S2x*S2x + S2y*S2y + S2z*S2z);

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
    {   case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
	    wspin6Avg = params->wdot6S1O*LNhdotS1 + params->wdot6S2O*LNhdotS2
	  + params->wdot6S1S2Avg*S1dotS2 + params->wdot6S1OS2OAvg*LNhdotS1*LNhdotS2
	  + (params->wdot6S1S1Avg + params->wdot6QMS1S1Avg)*S1sq
	  + (params->wdot6S2S2Avg + params->wdot6QMS2S2Avg)*S2sq
	  + (params->wdot6S1OS1OAvg + params->wdot6QMS1OS1OAvg)*LNhdotS1*LNhdotS1 +
	  + (params->wdot6S2OS2OAvg + params->wdot6QMS2OS2OAvg)*LNhdotS2*LNhdotS2;  
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
            // Compute 2.5PN SO correction to domega/dt
            // See Eq. 8.3 of gr-qc/0605140v4
            // Note that S_l/M^2 = (m1/M)^2 chi1 + (m2/M)^2 chi2
            // and Sigma_l/M^2 = (m2/M) chi2 - (m1/M) chi1
            wspin5 = params->wdot5S1O*LNhdotS1 + params->wdot5S2O*LNhdotS2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
	    wspin4Avg = params->wdot4S1S2Avg *S1dotS2 + params->wdot4S1OS2OAvg *LNhdotS1 * LNhdotS2;
            // Compute 2PN QM and self-spin corrections to domega/dt
            // See last line of Eq. 5.17 of arXiv:0812.4413
            // Also note this is equivalent to Eqs. 9c + 9d of astro-ph/0504538
            wspin4Avg += (params->wdot4S1S1Avg + params->wdot4QMS1S1Avg) * S1sq
	      + (params->wdot4S2S2Avg + params->wdot4QMS2S2Avg) * S2sq
	      + (params->wdot4S1OS1OAvg + params->wdot4QMS1OS1OAvg) * LNhdotS1 * LNhdotS1
	      + (params->wdot4S2OS2OAvg + params->wdot4QMS2OS2OAvg) * LNhdotS2 * LNhdotS2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            // Compute 1.5PN SO correction to domega/dt
            wspin3 = params->wdot3S1O*LNhdotS1 + params->wdot3S2O*LNhdotS2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
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

    domega  = params->wdotnewt * v11 / ( params->wdotcoeff[0]
            + v * ( params->wdotcoeff[1]
            + v * ( params->wdotcoeff[2]
            + v * ( params->wdotcoeff[3] + wspin3
            + v * ( params->wdotcoeff[4] + wspin4Avg
            + v * ( params->wdotcoeff[5] + wspin5
            + v * ( params->wdotcoeff[6] + wspin6Avg
                    + params->wdotlogcoeff * log(v)
            + v * ( params->wdotcoeff[7]
            + omega * ( params->wdottidal10
            + v*v * ( params->wdottidal12 ) ) ) ) ) ) ) ) ) );

    XLALSimInspiralSpinDerivativesAvg(&dLNhx,&dLNhy,&dLNhz,&dE1x,&dE1y,&dE1z,&dS1x,&dS1y,&dS1z,&dS2x,&dS2y,&dS2z,v,LNhx,LNhy,LNhz,E1x,E1y,E1z,S1x,S1y,S1z,S2x,S2y,S2z,LNhdotS1,LNhdotS2,params);

    /* dphi = d \phi / d \hat{t} = M d \phi /dt = M \omega = \hat{\omega} */
    dvalues[0]    = omega*(1.+omega*omega*omegashift(S1sq,S2sq,S1dotS2,LNhdotS1,LNhdotS2,params->omegashiftS1,params->omegashiftS2));
    dvalues[1]    = domega;
    dvalues[2]    = dLNhx; dvalues[3]     = dLNhy ; dvalues[4]    = dLNhz;
    dvalues[5]    = dS1x ; dvalues[6]     = dS1y  ; dvalues[7]    = dS1z ;
    dvalues[8]    = dS2x ; dvalues[9]     = dS2y  ; dvalues[10]   = dS2z ;
    dvalues[11]   = dE1x ; dvalues[12]    = dE1y  ; dvalues[13]   = dE1z ;

    return GSL_SUCCESS;
} /* End of XLALSimInspiralSpinTaylorT5DerivativesAvg() */

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
 * Driver function to generate any of SpinTaylorT1/T5/T4
 * If the output entries are not null the relative variables are returned.
 * The dictionary LALparams is searched for
 * PNAmplitudeOrder, amplitude PN order, default -1, i.e. maximum order available
 * PNPhaseOrder, phasing PN order, default -1, i.e. maximum order available
 * PNSpinOrder, spin-dependent term maximum order in PN-phasing, default -1, they start at 1.5PN order
 * PNtidalOrder, tidal PN order, default -1, they start at 5PN order (spin-less case only)
 * Lscorr, flag to include spin correction to orbital angular momentum (1), default value 0, i.e. spin corrections not included
 * QuadMon1, value of quadrupole-monopole dimension-less coupling for component 1 diminished by 1, default 0 (BH case in GR)
 * QuadMon2, value of quadrupole-monopole dimension-less coupling for component 2 diminished by 1, default 0 (BH case in GR)
 * Lambda1, value of the dimension-less tidal parameter for component 1, default 0 (BH in GR)
 * Lambda2, value of the dimension-less tidal parameter for component 2, default 0 (BH in GR)
 * FinalFreq, value of the final frequency to integrate to, in Hz
 * OnlyFinal, flag to only output the final values and not output the waveform (1), default value 0, i.e., intermediate values and waveform output
 * To insert a value named "NAME" into the LALDict dictionary use the instruction
 *   XLALSimInspiralWaveformParamsInsertNAME(LALparams, value)
 * to read a value named "NAME" from the LALDict dictionary use the instruction
 *   XLALSimInspiralWaveformParamsLookupNAME(LALparams, value)
 * NAMEs checked over are
 * * PNAmplitudeOrder
 * * PNPhaseOrder
 * * PNSpinOrder
 * * PNTidalOrder
 * * Lscorr
 * * dQuadMon1
 * * dQuadMon2
 * * TidalLambda1
 * * TidalLambda2
 * * FinalFreq
 * * OnlyFinal
 * Next-to-last version REVIEWED completed on git hash 6640e79e60791d5230731acc63351676ce7ce413
 */
int XLALSimInspiralSpinTaylorDriver(
	REAL8TimeSeries **hplus,        /**< +-polarization waveform */
	REAL8TimeSeries **hcross,       /**< x-polarization waveform */
	REAL8TimeSeries **Vout,         /**< PN parameter v */
	REAL8TimeSeries **Phiout,       /**< main GW phase */
	REAL8TimeSeries **S1xout,       /**< spin1 x-component*/
	REAL8TimeSeries **S1yout,       /**< spin1 y-component*/
	REAL8TimeSeries **S1zout,       /**< spin1 z-component*/
	REAL8TimeSeries **S2xout,       /**< spin2 x-component*/
	REAL8TimeSeries **S2yout,       /**< spin2 y-component*/
	REAL8TimeSeries **S2zout,       /**< spin2 z-component*/
	REAL8TimeSeries **LNhxout,      /**< Newtonian L unit vector x-component*/
	REAL8TimeSeries **LNhyout,      /**< Newtonian L unit vector y-component*/
	REAL8TimeSeries **LNhzout,      /**< Newtonian L unit vector <-component*/
	REAL8TimeSeries **E1xout,       /**< x-axis triad unit vector x-component*/
	REAL8TimeSeries **E1yout,       /**< x-axis triad unit vector y-component*/
	REAL8TimeSeries **E1zout,       /**< x-axis triad unit vector z-component*/
	const REAL8 phiRef,                   /**< orbital phase at reference pt. */
	const REAL8 deltaT,                   /**< sampling interval (s) */
	const REAL8 m1_SI,                    /**< mass of companion 1 (kg) */
	const REAL8 m2_SI,                    /**< mass of companion 2 (kg) */
	const REAL8 fStart,                   /**< start GW frequency (Hz) */
	const REAL8 fRef,                     /**< reference GW frequency (Hz) */
	const REAL8 r,                        /**< distance of source (m) */
	const REAL8 s1x,                      /**< initial value of S1x */
	const REAL8 s1y,                      /**< initial value of S1y */
	const REAL8 s1z,                      /**< initial value of S1z */
	const REAL8 s2x,                      /**< initial value of S2x */
	const REAL8 s2y,                      /**< initial value of S2y */
	const REAL8 s2z,                      /**< initial value of S2z */
	const REAL8 lnhatx,                   /**< initial value of LNhatx */
	const REAL8 lnhaty,                   /**< initial value of LNhaty */
	const REAL8 lnhatz,                   /**< initial value of LNhatz */
	const REAL8 e1x,                      /**< initial value of E1x */
	const REAL8 e1y,                      /**< initial value of E1y */
	const REAL8 e1z,                      /**< initial value of E1z */
	LALDict *LALparams,                   /**< LAL dictionary containing accessory parameters */
        const Approximant approx              /**< PN approximant (SpinTaylorT1/T5/T4) */
				    )
{
  if ( hplus )
    if ( *hplus ) {
      XLALPrintError("** LALSimIMRPSpinInspiralRD ERROR **: (*h+) is expected to be NULL; got %p\n",*hplus);
      XLAL_ERROR( XLAL_EFAULT );
    }
  if ( hcross )
    if ( *hcross ) {
      XLALPrintError("** LALSimIMRPSpinInspiralRD ERROR **: (*hx) is expected to be NULL; got %p\n",*hcross);
      XLAL_ERROR( XLAL_EFAULT );
    }
  if ( Phiout )
    if ( *Phiout ) {
      XLALPrintError("** LALSimIMRPSpinInspiralRD ERROR **: (*Phiout) is expected to be NULL; got %p\n",*Phiout);
      XLAL_ERROR( XLAL_EFAULT );
    }
  if ( S1xout )
    if ( *S1xout ) {
      XLALPrintError("** LALSimIMRPSpinInspiralRD ERROR **: (*S1xout) is expected to be NULL; got %p\n",*S1xout);
      XLAL_ERROR( XLAL_EFAULT );
    }
  if ( S1yout )
    if ( *S1yout ) {
      XLALPrintError("** LALSimIMRPSpinInspiralRD ERROR **: (*S1yout) is expected to be NULL; got %p\n",*S1yout);
      XLAL_ERROR( XLAL_EFAULT );
    }
  if ( S1zout )
    if ( *S1zout ) {
      XLALPrintError("** LALSimIMRPSpinInspiralRD ERROR **: (*S1zout) is expected to be NULL; got %p\n",*S1zout);
      XLAL_ERROR( XLAL_EFAULT );
    }
  if ( S2xout )
    if ( *S2xout ) {
      XLALPrintError("** LALSimIMRPSpinInspiralRD ERROR **: (*S2xout) is expected to be NULL; got %p\n",*S2xout);
      XLAL_ERROR( XLAL_EFAULT );
    }
  if ( S2yout )
    if ( *S2yout ) {
      XLALPrintError("** LALSimIMRPSpinInspiralRD ERROR **: (*S2yout) is expected to be NULL; got %p\n",*S2yout);
      XLAL_ERROR( XLAL_EFAULT );
    }
  if ( S2zout )
    if ( *S2zout ) {
      XLALPrintError("** LALSimIMRPSpinInspiralRD ERROR **: (*S2zout) is expected to be NULL; got %p\n",*S2zout);
      XLAL_ERROR( XLAL_EFAULT );
    }
  if ( LNhxout )
    if ( *LNhxout ) {
      XLALPrintError("** LALSimIMRPSpinInspiralRD ERROR **: (*LNhxout) is expected to be NULL; got %p\n",*LNhxout);
      XLAL_ERROR( XLAL_EFAULT );
    }
  if ( LNhyout )
    if ( *LNhyout ) {
      XLALPrintError("** LALSimIMRPSpinInspiralRD ERROR **: (*LNhyout) is expected to be NULL; got %p\n",*LNhyout);
      XLAL_ERROR( XLAL_EFAULT );
    }
  if ( LNhzout )
    if ( *LNhzout ) {
      XLALPrintError("** LALSimIMRPSpinInspiralRD ERROR **: (*LNhzout) is expected to be NULL; got %p\n",*LNhzout);
      XLAL_ERROR( XLAL_EFAULT );
    }
  if ( E1xout )
    if ( *E1xout ) {
      XLALPrintError("** LALSimIMRPSpinInspiralRD ERROR **: (*E1xout) is expected to be NULL; got %p\n",*E1xout);
      XLAL_ERROR( XLAL_EFAULT );
    }
  if ( E1yout )
    if ( *E1yout ) {
      XLALPrintError("** LALSimIMRPSpinInspiralRD ERROR **: (*E1yout) is expected to be NULL; got %p\n",*E1yout);
      XLAL_ERROR( XLAL_EFAULT );
    }
  if ( E1zout )
    if ( *E1zout ) {
      XLALPrintError("** LALSimIMRPSpinInspiralRD ERROR **: (*E1zout) is expected to be NULL; got %p\n",*E1zout);
      XLAL_ERROR( XLAL_EFAULT );
    }

    REAL8TimeSeries *V, *Phi, *S1x, *S1y, *S1z, *S2x, *S2y, *S2z;
    REAL8TimeSeries *LNhatx, *LNhaty, *LNhatz, *E1x, *E1y, *E1z;
    int status=0, n;
    unsigned int i;
    REAL8 fS, fE, phiShift;
    /* The Schwarzschild ISCO frequency - for sanity checking fRef */
    REAL8 fISCO = 1./(pow(6.,3./2.)*LAL_PI*(m1_SI+m2_SI)/LAL_MSUN_SI*LAL_MTSUN_SI);

    INT4 ampO   = XLALSimInspiralWaveformParamsLookupPNAmplitudeOrder(LALparams);
    INT4 phaseO = XLALSimInspiralWaveformParamsLookupPNPhaseOrder(LALparams);
    INT4 spinO  = XLALSimInspiralWaveformParamsLookupPNSpinOrder(LALparams);
    INT4 tideO  = XLALSimInspiralWaveformParamsLookupPNTidalOrder(LALparams);
    INT4 lscorr = XLALSimInspiralWaveformParamsLookupLscorr(LALparams);
    REAL8 quadparam1 = 1.+XLALSimInspiralWaveformParamsLookupdQuadMon1(LALparams);
    REAL8 quadparam2 = 1.+XLALSimInspiralWaveformParamsLookupdQuadMon2(LALparams);
    REAL8 lambda1 = XLALSimInspiralWaveformParamsLookupTidalLambda1(LALparams);
    REAL8 lambda2 = XLALSimInspiralWaveformParamsLookupTidalLambda2(LALparams);

    /* phiRef is used to rotate axis, not to shift phi*/
    REAL8 e1xtmp=e1x;
    REAL8 e1ytmp=e1y;
    REAL8 e1ztmp=e1z;
    REAL8 e1xp,e1yp,e1zp,e1xphi,e1yphi,e1zphi,lxy;
    REAL8 alpha=atan2(lnhaty,lnhatx);
    //Rotation in the x-y plane to bring L along x
    lxy=lnhatx*cos(alpha)+lnhaty*sin(alpha);
    e1xp=e1xtmp*cos(alpha)+e1ytmp*sin(alpha);
    e1yp=e1ytmp*cos(alpha)-e1xtmp*sin(alpha);
    e1xtmp=e1xp;
    e1ytmp=e1yp;
    REAL8 iota=atan2(lxy,lnhatz);
    //Rotation in the x-z plane to bring L along z 
    e1xp=e1xtmp*cos(iota)-e1ztmp*sin(iota);
    e1zp=e1ztmp*cos(iota)+e1xtmp*sin(iota);
    e1xtmp=e1xp;
    e1ztmp=e1zp;
    /* Now L is along z, we can now rotate by phiRef in the x-y plane */
    e1xp=e1xtmp*cos(phiRef)-e1ytmp*sin(phiRef);
    e1yp=e1ytmp*cos(phiRef)+e1xtmp*sin(phiRef);
    e1xtmp=e1xp;
    e1ytmp=e1yp;
    /* Now we undo 2 rotations to go back to the initial frame */
    e1xp=e1xtmp*cos(iota)+e1ztmp*sin(iota);
    e1zp=e1ztmp*cos(iota)-e1xtmp*sin(iota);
    e1xphi=e1xp*cos(alpha)-e1yp*sin(alpha);
    e1yphi=e1yp*cos(alpha)+e1xp*sin(alpha);
    e1zphi=e1zp;
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

    if( XLALSimInspiralWaveformParamsLookupOnlyFinal(LALparams) )
    {
        fS = fStart;
        fE = XLALSimInspiralWaveformParamsLookupFinalFreq(LALparams);
        /* Evolve the dynamical variables */
        n = XLALSimInspiralSpinTaylorPNEvolveOrbitOnlyFinal(&V, &Phi,
                &S1x, &S1y, &S1z, &S2x, &S2y, &S2z,
                &LNhatx, &LNhaty, &LNhatz, &E1x, &E1y, &E1z,
                deltaT, m1_SI, m2_SI, fS, fE, s1x, s1y, s1z, s2x, s2y, s2z,
                lnhatx, lnhaty, lnhatz, e1xphi, e1yphi, e1zphi, lambda1, lambda2,
	        quadparam1, quadparam2, spinO, tideO, phaseO, lscorr, approx);
        if( n < 0 )
            XLAL_ERROR(XLAL_EFUNC);

        /* Set orbital phase to zero */
        Phi->data->data[0] = 0.;
    }
    else
    {
    /* if fRef=0, just integrate from start to end. Let phiRef=phiC */
    if( fRef < LAL_REAL4_EPS )
    {
        fS = fStart;
        fE = XLALSimInspiralWaveformParamsLookupFinalFreq(LALparams);
        /* Evolve the dynamical variables */
        n = XLALSimInspiralSpinTaylorPNEvolveOrbit(&V, &Phi,
                &S1x, &S1y, &S1z, &S2x, &S2y, &S2z,
                &LNhatx, &LNhaty, &LNhatz, &E1x, &E1y, &E1z,
                deltaT, m1_SI, m2_SI, fS, fE, s1x, s1y, s1z, s2x, s2y, s2z,
                lnhatx, lnhaty, lnhatz, e1xphi, e1yphi, e1zphi, lambda1, lambda2,
	        quadparam1, quadparam2, spinO, tideO, phaseO, lscorr, approx);
        if( n < 0 )
            XLAL_ERROR(XLAL_EFUNC);

        /* Apply phase shift so orbital phase ends with desired value */
        phiShift = /* phiRef*/- Phi->data->data[Phi->data->length-1];
	/* phiRef is used to rotate axis, not to shift phi*/
        for( i=0; i < Phi->data->length; i++)
        {
            Phi->data->data[i] += phiShift;
        }
    }
    /* if fRef=fStart, just integrate from start to end. Let phiRef=phiStart */
    else if( fabs(fRef - fStart) < LAL_REAL4_EPS )
    {
        fS = fStart;
        fE = XLALSimInspiralWaveformParamsLookupFinalFreq(LALparams);
        /* Evolve the dynamical variables */
        n = XLALSimInspiralSpinTaylorPNEvolveOrbit(&V, &Phi,
                &S1x, &S1y, &S1z, &S2x, &S2y, &S2z,
                &LNhatx, &LNhaty, &LNhatz, &E1x, &E1y, &E1z,
                deltaT, m1_SI, m2_SI, fS, fE, s1x, s1y, s1z, s2x, s2y, s2z,
                lnhatx, lnhaty, lnhatz, e1xphi, e1yphi, e1zphi, lambda1, lambda2,
		quadparam1, quadparam2, spinO, tideO, phaseO, lscorr, approx);
        if( n < 0 )
            XLAL_ERROR(XLAL_EFUNC);

        /* Apply phase shift so orbital phase starts with desired value */
        phiShift = /*phiRef*/ - Phi->data->data[0];
	/* phiRef is used to rotate axis, not to shift phi*/
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
                deltaT, m1_SI, m2_SI, fS, fE, s1x, s1y, s1z, s2x, s2y,
                s2z, lnhatx, lnhaty, lnhatz, e1xphi, e1yphi, e1zphi, lambda1, lambda2,
	        quadparam1, quadparam2, spinO, tideO, phaseO, lscorr, approx);
        if( n < 0 )
        {
            XLAL_ERROR(XLAL_EFUNC);
        }

        /* Apply phase shift so orbital phase has desired value at fRef */
        phiShift = /*phiRef*/ - Phi1->data->data[Phi1->data->length-1];
	/* phiRef is used to rotate axis, not to shift phi*/
        for( i=0; i < Phi1->data->length; i++)
        {
            Phi1->data->data[i] += phiShift;
        }

        /* Integrate forward to end of waveform */
        fS = fRef;
        fE = XLALSimInspiralWaveformParamsLookupFinalFreq(LALparams);
        n = XLALSimInspiralSpinTaylorPNEvolveOrbit(&V2, &Phi2,
                &S1x2, &S1y2, &S1z2, &S2x2, &S2y2, &S2z2,
                &LNhatx2, &LNhaty2, &LNhatz2, &E1x2, &E1y2, &E1z2,
                deltaT, m1_SI, m2_SI, fS, fE, s1x, s1y, s1z, s2x, s2y,
                s2z, lnhatx, lnhaty, lnhatz, e1xphi, e1yphi, e1zphi, lambda1, lambda2,
		quadparam1, quadparam2, spinO, tideO, phaseO, lscorr, approx);
        if( n < 0 )
        {
            XLAL_ERROR(XLAL_EFUNC);
        }

        /* Apply phase shift so orbital phase has desired value at fRef */
        phiShift = /*phiRef*/ - Phi2->data->data[0];
	/* phiRef is used to rotate axis, not to shift phi*/
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
    if (hplus && hcross) {
      status = XLALSimInspiralPrecessingPolarizationWaveforms(hplus, hcross,
	        V, Phi, S1x, S1y, S1z, S2x, S2y, S2z, LNhatx, LNhaty, LNhatz,
                E1x, E1y, E1z, m1_SI, m2_SI, r, ampO);
    }
    }

    /* Destroy vectors of dynamical variables, check for errors then exit */
    if (Vout)
      *Vout=V;
    else
      XLALDestroyREAL8TimeSeries(V);
    if (Phiout)
      *Phiout=Phi;
    else
      XLALDestroyREAL8TimeSeries(Phi);
    if (S1xout)
      *S1xout=S1x;
    else
      XLALDestroyREAL8TimeSeries(S1x);
    if (S1yout)
      *S1yout=S1y;
    else
      XLALDestroyREAL8TimeSeries(S1y);
    if (S1zout)
      *S1zout=S1z;
    else
      XLALDestroyREAL8TimeSeries(S1z);
    if (S2xout)
      *S2xout=S2x;
    else
      XLALDestroyREAL8TimeSeries(S2x);
    if (S2yout)
      *S2yout=S2y;
    else
      XLALDestroyREAL8TimeSeries(S2y);
    if (S2zout)
      *S2zout=S2z;
    else
      XLALDestroyREAL8TimeSeries(S2z);
    if (LNhxout)
      *LNhxout=LNhatx;
    else
      XLALDestroyREAL8TimeSeries(LNhatx);
    if (LNhyout)
      *LNhyout=LNhaty;
    else
      XLALDestroyREAL8TimeSeries(LNhaty);
    if (LNhzout)
      *LNhzout=LNhatz;
    else
      XLALDestroyREAL8TimeSeries(LNhatz);
    if (E1xout)
      *E1xout=E1x;
    else
      XLALDestroyREAL8TimeSeries(E1x);
    if (E1yout)
      *E1yout=E1y;
    else
      XLALDestroyREAL8TimeSeries(E1y);
    if (E1zout)
      *E1zout=E1z;
    else
      XLALDestroyREAL8TimeSeries(E1z);

    if( status < 0 )
        XLAL_ERROR(XLAL_EFUNC);

    return n;
} //End of XLALSimInspiralSpinTaylorDriver()

/**
 * Driver function wrapper for XLALSimInspiralSpinTaylorDriver()
 * with the same output except for h+ and hx time series
 */
int XLALSimInspiralSpinTaylorOrbitalDriver(
	REAL8TimeSeries **Vout,         /**< PN parameter v */
	REAL8TimeSeries **Phiout,       /**< main GW phase */
	REAL8TimeSeries **S1xout,       /**< spin1 x-component*/
	REAL8TimeSeries **S1yout,       /**< spin1 y-component*/
	REAL8TimeSeries **S1zout,       /**< spin1 z-component*/
	REAL8TimeSeries **S2xout,       /**< spin2 x-component*/
	REAL8TimeSeries **S2yout,       /**< spin2 y-component*/
	REAL8TimeSeries **S2zout,       /**< spin2 z-component*/
	REAL8TimeSeries **LNhxout,      /**< Newtonian L unit vector x-component*/
	REAL8TimeSeries **LNhyout,      /**< Newtonian L unit vector y-component*/
	REAL8TimeSeries **LNhzout,      /**< Newtonian L unit vector <-component*/
	REAL8TimeSeries **E1xout,       /**< x-axis triad unit vector x-component*/
	REAL8TimeSeries **E1yout,       /**< x-axis triad unit vector y-component*/
	REAL8TimeSeries **E1zout,       /**< x-axis triad unit vector z-component*/
	const REAL8 phiRef,                   /**< orbital phase at reference pt. */
	const REAL8 deltaT,                   /**< sampling interval (s) */
	const REAL8 m1_SI,                    /**< mass of companion 1 (kg) */
	const REAL8 m2_SI,                    /**< mass of companion 2 (kg) */
	const REAL8 fStart,                   /**< start GW frequency (Hz) */
	const REAL8 fRef,                     /**< reference GW frequency (Hz) */
	const REAL8 s1x,                      /**< initial value of S1x */
	const REAL8 s1y,                      /**< initial value of S1y */
	const REAL8 s1z,                      /**< initial value of S1z */
	const REAL8 s2x,                      /**< initial value of S2x */
	const REAL8 s2y,                      /**< initial value of S2y */
	const REAL8 s2z,                      /**< initial value of S2z */
	const REAL8 lnhatx,                   /**< initial value of LNhatx */
	const REAL8 lnhaty,                   /**< initial value of LNhaty */
	const REAL8 lnhatz,                   /**< initial value of LNhatz */
	const REAL8 e1x,                      /**< initial value of E1x */
	const REAL8 e1y,                      /**< initial value of E1y */
	const REAL8 e1z,                      /**< initial value of E1z */
	LALDict *LALparams,                   /**< LAL dictionary containing accessory parameters */
        const Approximant approx              /**< PN approximant (SpinTaylorT1/T5/T4) */
	)
{
  /* Since this function does not return polarizations,
   * the value of distance to pass to XLALSimInspiralSpinTaylorDriver()
   * is immaterial and set to 1 here.
   */
  return XLALSimInspiralSpinTaylorDriver(NULL,NULL,Vout,Phiout,S1xout,S1yout,S1zout,S2xout,S2yout,S2zout,LNhxout,LNhyout,LNhzout,E1xout,E1yout,E1zout,phiRef,deltaT,m1_SI,m2_SI,fStart,fRef,1.,s1x,s1y,s1z,s2x,s2y,s2z,lnhatx,lnhaty,lnhatz,e1x,e1y,e1z,LALparams,approx);
}

INT4 XLALSimInspiralSpinTaylorHlmModesFromOrbit(SphHarmTimeSeries **hlm, /**< OUTPUT */
						REAL8TimeSeries*V,     /**< PN parameter v */
						REAL8TimeSeries *Phi,  /**< orbital phase */
						REAL8TimeSeries *LNhx, /**< angular momentum unit vector x component */
						REAL8TimeSeries *LNhy, /**< angular momentum unit vector y component */
						REAL8TimeSeries *LNhz, /**< angular momentum unit vector z components */
						REAL8TimeSeries *E1x,  /**< x-axis tetrad x component*/
						REAL8TimeSeries *E1y,  /**< x-axis tetrad y component*/
						REAL8TimeSeries *E1z,  /**< x-axis tetrad z component*/
						REAL8TimeSeries *S1x,  /**< spin1-x component */
						REAL8TimeSeries *S1y,  /**< spin1-y component */
						REAL8TimeSeries *S1z,  /**< spin1-z component */
						REAL8TimeSeries *S2x,  /**< spin2-x component */
						REAL8TimeSeries *S2y,  /**< spin2-y component */
						REAL8TimeSeries *S2z,  /**< spin2-z component */
						REAL8 m1_SI,           /**< mass of companion 1 (kg) */
						REAL8 m2_SI,           /**< mass of companion 2 (kg) */
						REAL8 dist_SI,         /**< distance of source (m) */
						int ampO,              /**< twice post-Newtonian amp-order */
						LALValue *modearray    /**< Container for the ell and m modes to generate. To generate all available modes pass NULL */)
{
  const UINT4 lmax=4;
  UINT4 l;
  INT4 m;
  SphHarmTimeSeries *hlm_tmp=NULL;
  LALValue *modearray_int = XLALSimInspiralCreateModeArray();
  INT4 act, errCode=0;
  if ( modearray == NULL )
    for (l=2; l<=lmax; l++)
      XLALSimInspiralModeArrayActivateAllModesAtL(modearray_int, l);
  else
    for (l=2; l <= lmax; l++) {
      m=-l-1;
      act=0;
      do {
	m++;
	act=XLALSimInspiralModeArrayIsModeActive(modearray, l, m);
      } while ( (m<(INT4)l) && (act==0) );
      if (act==1)
	XLALSimInspiralModeArrayActivateAllModesAtL(modearray_int, l);
    }
  for (l=2; l<=lmax; l++) {
    m=l;
    if (XLALSimInspiralModeArrayIsModeActive(modearray_int, l, m) == 1 ) {
      if (l==2)
	errCode=XLALSimInspiralSpinPNMode2m(&hlm_tmp,V,Phi,LNhx,LNhy,LNhz,E1x,E1y,E1z,S1x,S1y,S1z,S2x,S2y,S2z,m1_SI,m2_SI,dist_SI,ampO);
      else if (l==3)
	errCode=XLALSimInspiralSpinPNMode3m(&hlm_tmp,V,Phi,LNhx,LNhy,LNhz,E1x,E1y,E1z,S1x,S1y,S1z,S2x,S2y,S2z,m1_SI,m2_SI,dist_SI,ampO);
      else if (l==4)
	errCode=XLALSimInspiralSpinPNMode4m(&hlm_tmp,V,Phi,LNhx,LNhy,LNhz,E1x,E1y,E1z,S1x,S1y,S1z,S2x,S2y,S2z,m1_SI,m2_SI,dist_SI,ampO);
      for (m=-l; m<=(INT4)l; m++)
	*hlm=XLALSphHarmTimeSeriesAddMode(*hlm,XLALSphHarmTimeSeriesGetMode(hlm_tmp,l,m),l,m);
    }
  }

  return errCode;
}

INT4 XLALSimInspiralSpinTaylorHlmModes(SphHarmTimeSeries **hlm,  /**< OUTPUT */
				       REAL8 phiRef,  /**< orbital phase at reference pt. */
				       REAL8 dT,      /**< sampling interval (s) */
				       REAL8 m1_SI,   /**< mass of companion 1 (kg) */
				       REAL8 m2_SI,   /**< mass of companion 2 (kg) */
				       REAL8 fStart,  /**< start GW frequency (Hz) */
				       REAL8 fRef,    /**< reference GW frequency (Hz) */
				       REAL8 dist_SI, /**< distance of source (m) */
				       REAL8 s1x,     /**< ref value of S1x */
				       REAL8 s1y,     /**< ref value of S1y */
				       REAL8 s1z,     /**< ref value of S1z */
				       REAL8 s2x,     /**< ref value of S2x */
				       REAL8 s2y,     /**< ref value of S2y */
				       REAL8 s2z,     /**< ref value of S2z */
				       REAL8 lnhatx,  /**< ref value of LNhatx */
				       REAL8 lnhaty,  /**< ref value of LNhaty */
				       REAL8 lnhatz,  /**< ref value of LNhatz */
				       REAL8 e1x,     /**< ref value of E1x */
				       REAL8 e1y,     /**< ref value of E1y */
				       REAL8 e1z,     /**< ref value of E1z */
				       int ampO,      /**< twice post-Newtonian amp-order */
				       LALValue *modearray, /**< Container for the ell and m modes to generate. To generate all available modes pass NULL */
				       LALDict *LALparams,  /**< LAL dictionary containing accessory parameters */
				       Approximant approx   /**< PN approximant (SpinTaylorT1/T2/T4) */
				       )
{
  REAL8TimeSeries *V=NULL;
  REAL8TimeSeries *Phi=NULL;
  REAL8TimeSeries *LNhx=NULL;
  REAL8TimeSeries *LNhy=NULL;
  REAL8TimeSeries *LNhz=NULL;
  REAL8TimeSeries *S1x=NULL;
  REAL8TimeSeries *S1y=NULL;
  REAL8TimeSeries *S1z=NULL;
  REAL8TimeSeries *S2x=NULL;
  REAL8TimeSeries *S2y=NULL;
  REAL8TimeSeries *S2z=NULL;
  REAL8TimeSeries *E1x=NULL;
  REAL8TimeSeries *E1y=NULL;
  REAL8TimeSeries *E1z=NULL;

  int err_code=XLALSimInspiralSpinTaylorDriver(NULL,NULL,&V,&Phi,&S1x,&S1y,&S1z,&S2x,&S2y,&S2z,&LNhx,&LNhy,&LNhz,&E1x,&E1y,&E1z,phiRef,dT,m1_SI,m2_SI,fStart,fRef,dist_SI,s1x,s1y,s1z,s2x,s2y,s2z,lnhatx,lnhaty,lnhatz,e1x,e1y,e1z,LALparams,approx);
  if (err_code!=XLAL_SUCCESS)
    XLAL_ERROR(XLAL_EINVAL);

  return XLALSimInspiralSpinTaylorHlmModesFromOrbit(hlm,V,Phi,LNhx,LNhy,LNhz,E1x,E1y,E1z,S1x,S1y,S1z,S2x,S2y,S2z,m1_SI,m2_SI,dist_SI,ampO,modearray);

}

/**
 * This function evolves the orbital equations for a precessing binary using
 * the \"TaylorT5/T4\" approximant for solving the orbital dynamics
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
 *
 * !!!UNREVIEWED!!!
 *
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
    Approximant approx              /**< PN approximant (SpinTaylorT1/T5/T4) */
        )
{
    INT4 intreturn;
    void * params;
    LALAdaptiveRungeKuttaIntegrator *integrator = NULL;     /* GSL integrator object */
    REAL8 yinit[LAL_NUM_ST4_VARIABLES];       /* initial values of parameters */
    REAL8Array *y;    /* time series of variables returned from integrator */
    /* intermediate variables */
    UINT4 cutlen, len;
    REAL8 norm, dtStart, dtEnd, lengths, m1sec, m2sec, Msec, Mcsec, fTerm;

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
        XLALSimInspiralSpinTaylorTxCoeffs *paramsT4;
        XLALSimInspiralSpinTaylorT4Setup(&paramsT4, m1, m2, fStart, fEnd,
					 lambda1, lambda2, quadparam1, quadparam2, spinO, tideO, phaseO, 0, 0);
        params = (void *) &paramsT4;
    }
    else if( approx == SpinTaylorT5 )
    {
        XLALSimInspiralSpinTaylorTxCoeffs *paramsT5;
        XLALSimInspiralSpinTaylorT5Setup(&paramsT5, m1, m2, fStart, fEnd,
		lambda1, lambda2, quadparam1, quadparam2, spinO, tideO, phaseO, 0);
        params = (void *) &paramsT5;
    }
    else if( approx == SpinTaylorT1 )
    {
        XLALSimInspiralSpinTaylorTxCoeffs *paramsT1;
        XLALSimInspiralSpinTaylorT1Setup(&paramsT1, m1, m2, fStart, fEnd,
	        lambda1, lambda2, quadparam1, quadparam2, spinO, tideO, phaseO, 0);
        params = (void *) &paramsT1;
    }
    else
    {
        XLALPrintError("XLAL Error - %s: Approximant must be one of SpinTaylorT1, SpinTaylorT5, SpinTaylorT4, but %i provided\n",
                __func__, approx);
        XLAL_ERROR(XLAL_EINVAL);

    }
    m1sec = m1 /LAL_MSUN_SI*LAL_MTSUN_SI;
    m2sec = m2 /LAL_MSUN_SI*LAL_MTSUN_SI;
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
                XLALSimInspiralSpinTaylorT4DerivativesAvg,
                XLALSimInspiralSpinTaylorStoppingTest,
                LAL_ST4_ABSOLUTE_TOLERANCE, LAL_ST4_RELATIVE_TOLERANCE);
    else if( approx == SpinTaylorT5 )
        integrator = XLALAdaptiveRungeKutta4Init(LAL_NUM_ST4_VARIABLES,
                XLALSimInspiralSpinTaylorT5DerivativesAvg,
                XLALSimInspiralSpinTaylorStoppingTest,
                LAL_ST4_ABSOLUTE_TOLERANCE, LAL_ST4_RELATIVE_TOLERANCE);
    else if( approx == SpinTaylorT1 )
        integrator = XLALAdaptiveRungeKutta4Init(LAL_NUM_ST4_VARIABLES,
                XLALSimInspiralSpinTaylorT1DerivativesAvg,
                XLALSimInspiralSpinTaylorStoppingTest,
                LAL_ST4_ABSOLUTE_TOLERANCE, LAL_ST4_RELATIVE_TOLERANCE);
    else
    {
        XLALPrintError("XLAL Error - %s: Approximant must be one of SpinTaylorT1, SpinTaylorT5, SpinTaylorT4, but %i provided\n",
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
    len = XLALAdaptiveRungeKutta4IrregularIntervals(integrator, params, yinit,
            0.0, lengths/Msec, &y);

    intreturn = integrator->returncode;
    XLALAdaptiveRungeKuttaFree(integrator);

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
        __func__, intreturn, m1 / LAL_MSUN_SI * LAL_MTSUN_SI,
        m2 / LAL_MSUN_SI * LAL_MTSUN_SI, s1x, s1y, s1z, s2x,
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
 * @deprecated Use XLALSimInspiralChooseTDWaveform() instead
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
 *
 * !!!UNREVIEWED!!!
 */
int XLALSimInspiralTransformPrecessingObsoleteInitialConditions(
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
 * Function to specify the desired orientation of the spin components of
 * a precessing binary.
 *
 * Input:
 * Depending on the values of the variable frameChoice
 * input parameters may represent different quantities.
 * incl is the angle between in the X-Y-Z frame specified by
 * * frameChoice = LAL_SIM_INSPIRAL_FRAME_AXIS_TOTAL_J
 *   J // Z, N = (0, sin(inc), cos(inc))
 * * frameChoice = LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L (default)
 *   L // Z, N = (0, sin(inc), cos(inc))
 * * frameChoice = LAL_SIM_INSPIRAL_FRAME_AXIS_VIEW
 * * N // Z, L = (sin(inc), 0, cos(inc))
 *
 * Cartesian components of S1 and S2 are given wrt axes x-y-Z (ORBITAL_L and VIEW cases),
 * being x the vector pointing from body 2 to body 1.
 * x-y axes are rotated wrt X-Y by phiRef counterclockwise, i.e. if
 * S1=(a,b,0) in fame x-y-Z, in frame X-Y-Z it will have components
 * S1'=(a cos(phiRef) - b sin(phiRef), a sin(phiRef) + b cos(phiRef), 0).
 * In the TOTAL_J case spin components are given wrt J, with the x-Z plane
 * spanned by the line connecting the two bodies and J.
 *
 * m1, m2, f_ref are the component masses and reference GW frequency,
 * they are needed to compute the magnitude of L_N and.
 * phi_ref is needed since the spin components in the plane of the orbit
 * are defined wrt assuming the x axis is parallel to the line joining the binary consistuents,
 * from body 2 to body 1.
 *
 * Output:
 * x, y, z components S1 and S2 wrt N
 * inc angle between L and N, i.e. N=(0,0,1), LN=(sin(inc),0,cos(inc))
 * as required by precessing waveforms routines in SpinTaylor and IRMPhenomP codes.
 *
 * NOTE: Here the \"total\" angular momentum is computed as
 * J = L_N(1+l_1PN) + S1 + S2
 * where L_N(1+l_1PN) is the 1PN-corrected orbital angular momentum,
 * which is parallel to Newtonian angular momentum.
 * PN corrections to L from 1.5PN order on (wrt to Newtonian value)
 * are NOT ACCOUNTED FOR in this function.
 */
int XLALSimInspiralInitialConditionsPrecessingApproxs(
		REAL8 *inc,	/**< inclination angle (returned) */
		REAL8 *S1x,	/**< S1x component (returned) */
		REAL8 *S1y,	/**< S1y component (returned) */
		REAL8 *S1z,	/**< S1z component (returned) */
		REAL8 *S2x,	/**< S2x components (returned) */
		REAL8 *S2y,	/**< S2y components (returned) */
		REAL8 *S2z,	/**< S2z components (returned) */
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
		const REAL8 phiRef, /**< reference GW phase */
		LALSimInspiralFrameAxis frameChoice  /**< Flag to identify axis wrt which spin components are given. Pass in NULL (or None in python) for default (OrbitalL) */)
{
  *S1x=S1xIn;
  *S1y=S1yIn;
  *S1z=S1zIn;
  *S2x=S2xIn;
  *S2y=S2yIn;
  *S2z=S2zIn;

  REAL8 Lmag=0.;
  REAL8 Lx,Ly,Lxy2,Lz;
  REAL8 LNhx,LNhy,LNhz;
  REAL8 v0,mass1,mass2,eta;
  REAL8 tmp1,tmp2;
  //Uncomment the following line to perform tests
  REAL8 Nhx,Nhy,Nhz;
  switch (frameChoice) {
  /* FRAME_AXIS_TOTAL_J
   * (spins wrt to J, inclIn is the angle between J and N: if
   * J=(0,0,1) N=(0,sin(inclIn),cos(inclIn)))
   */
  case LAL_SIM_INSPIRAL_FRAME_AXIS_TOTAL_J:
    mass1=m1/LAL_MSUN_SI;
    mass2=m2/LAL_MSUN_SI;
    eta=mass1*mass2/(mass1+mass2)/(mass1+mass2);
    v0=cbrt(LAL_PI*fRef*(mass1+mass2)*LAL_MTSUN_SI);
    Lmag = XLALSimInspiralLN(mass1+mass2,eta,v0)*(1.+v0*v0*XLALSimInspiralL_2PN(eta));
    Lx  = -(*S1x)*mass1*mass1-(*S2x)*mass2*mass2;
    Ly  = -(*S1y)*mass1*mass1-(*S2y)*mass2*mass2;
    Lxy2= Lx*Lx+Ly*Ly;
    if (Lmag*Lmag>=Lxy2) {
      Lz=sqrt(Lmag*Lmag-Lxy2);
      if ( Lz<(S1zIn*mass1*mass1+S2zIn*mass2*mass2) ) {
	XLALPrintError("** XLALSimInspiralInitialConditionsPrecessingApproxs ERROR *** for s1 (%12.4e  %12.4e  %12.4e)\n",S1xIn,S1yIn,S1zIn);
	XLALPrintError("                                                                   s2 (%12.4e  %12.4e  %12.4e)\n",S2xIn,S2yIn,S2zIn);
	XLALPrintError(" wrt to J for m: (%12.4e  %12.4e) and v= %12.4e\n",mass1,mass2,v0);
	XLALPrintError(" it is impossible to determine the sign of Lhz\n");
	XLAL_ERROR(XLAL_EDOM);
      }
    }
    else {
      XLALPrintError("** XLALSimInspiralInitialConditionsPrecessingApproxs ERROR *** unphysical values of s1 (%12.4e  %12.4e  %12.4e)\n",S1xIn,S1yIn,S1zIn);
      XLALPrintError("                                                                                    s2 (%12.4e  %12.4e  %12.4e)\n",S2xIn,S2yIn,S2zIn);
      XLALPrintError(" wrt to J for m: (%12.4e  %12.4e) and v= %12.4e\n",mass1,mass2,v0);
      XLAL_ERROR(XLAL_EDOM);
    }
    // 1PN corrections to L are collinear with LN
    // Rotation1: rotate N into z
    LNhx=Lx/Lmag;
    LNhy=Ly/Lmag;
    LNhz=Lz/Lmag;
    Nhx=0.;
    Nhy=sin(inclIn);
    Nhz=cos(inclIn);
    //Rotate to the orbital plane to rotate spin by phiRef there
    // R1:
    REAL8 phiLJ=atan2(LNhy,LNhx);
    ROTATEZ(-phiLJ,*S1x,*S1y,*S1z);
    ROTATEZ(-phiLJ,*S2x,*S2y,*S2z);
    ROTATEZ(-phiLJ,Nhx,Nhy,Nhz);
    ROTATEZ(-phiLJ,LNhx,LNhy,LNhz);
    // R2:
    REAL8 thetaLJ=acos(LNhz);
    ROTATEY(-thetaLJ,*S1x,*S1y,*S1z);
    ROTATEY(-thetaLJ,*S2x,*S2y,*S2z);
    ROTATEY(-thetaLJ,Nhx,Nhy,Nhz);
    ROTATEY(-thetaLJ,LNhx,LNhy,LNhz);
    // R3: (spins only)
    ROTATEZ(phiRef,*S1x,*S1y,*S1z);
    ROTATEZ(phiRef,*S2x,*S2y,*S2z);
    //Now let's go to the radiation frame
    // R4:
    REAL8 phiN=atan2(Nhy,Nhx);
    ROTATEZ(-phiN,*S1x,*S1y,*S1z);
    ROTATEZ(-phiN,*S2x,*S2y,*S2z);
    ROTATEZ(-phiN,LNhx,LNhy,LNhz);
    ROTATEZ(-phiN,Nhx,Nhy,Nhz);
    // R5:
    REAL8 thetaN=acos(Nhz);
    ROTATEY(-thetaN,*S1x,*S1y,*S1z);
    ROTATEY(-thetaN,*S2x,*S2y,*S2z);
    ROTATEY(-thetaN,LNhx,LNhy,LNhz);
    ROTATEY(-thetaN,Nhx,Nhy,Nhz);
    // Now rotating LNh into the x,z plane
    // R6:
    REAL8 phiL=atan2(LNhy,LNhx);
    ROTATEZ(-phiL,*S1x,*S1y,*S1z);
    ROTATEZ(-phiL,*S2x,*S2y,*S2z);
    ROTATEZ(-phiL,LNhx,LNhy,LNhz);
    ROTATEZ(-phiL,Nhx,Nhy,Nhz);
    *inc=acos(LNhz);
    //It follows that in the frame in which N=(0,0,1), LNhat=(sin(*inc),0,cos(*inc))
    //For a check uncomment the following lines
    //Anyway there is an automatic check at build in
    // lalsimulation/test/InitialSpinRotationTest
    /*printf("**************************************************\n");
    printf("** Check of LAL_SIM_INSPIRAL_FRAME_AXIS_TOTAL_J **\n");
    printf("**************************************************\n");
    printf("N:   %12.4e  %12.4e  %12.4e  norm: %12.4e\n",Nhx,Nhy,Nhz,sqrt(Nhx*Nhx+Nhy*Nhy+Nhz*Nhz));
    printf("     %12.4e  %12.4e  %12.4e\n",0.,0.,1.);
    printf("LNh: %12.4e  %12.4e  %12.4e  norm: %12.4e\n",LNhx,LNhy,LNhz,sqrt(LNhx*LNhx+LNhy*LNhy+LNhz*LNhz));
    printf("     %12.4e  %12.4e  %12.4e\n",sin(*inc),0.,cos(*inc));
    printf("S1.L     %12.4e  %12.4e\n",(*S1x)*LNhx+(*S1y)*LNhy+(*S1z)*LNhz,(S1xIn*Lx+S1yIn*Ly+S1zIn*Lz)/Lmag);
    printf("S2.L     %12.4e  %12.4e\n",(*S2x)*LNhx+(*S2y)*LNhy+(*S2z)*LNhz,(S2xIn*Lx+S2yIn*Ly+S2zIn*Lz)/Lmag);
    printf("S1.S2    %12.4e  %12.4e\n",(*S1x)*(*S2x)+(*S1y)*(*S2y)+(*S1z)*(*S2z),S1xIn*S2xIn+S1yIn*S2yIn+S1zIn*S2zIn);
    printf("**************************************************\n");*/
    break;

  case LAL_SIM_INSPIRAL_FRAME_AXIS_VIEW:
    /* FRAME_AXIS_VIEW
     * Spins at phiRef=0 are given wrt to view direction
     * inclIn is the angle between L and N: if
     * N=(0,0,1) LNhat=(sin(inclIn),0,cos(inclIn))
     * We need to perform the rotation by phiRef around LNhat on the spins only.
     */
    *inc=inclIn;
    ROTATEY(-inclIn,*S1x,*S1y,*S1z);
    ROTATEY(-inclIn,*S2x,*S2y,*S2z);
    ROTATEZ(phiRef,*S1x,*S1y,*S1z);
    ROTATEZ(phiRef,*S2x,*S2y,*S2z);
    ROTATEY(inclIn,*S1x,*S1y,*S1z);
    ROTATEY(inclIn,*S2x,*S2y,*S2z);
    //For a check uncomment the following lines
    //Anyway there is an automatic check at build in
    // lalsimulation/test/InitialSpinRotationTest
    /*printf("***********************************************\n");
    printf("** Check of LAL_SIM_INSPIRAL_FRAME_AXIS_VIEW **\n");
    printf("***********************************************\n");
    LNhx=sin(inclIn);
    LNhy=0.;
    LNhz=cos(inclIn);
    ROTATEY(-inclIn, LNhx, LNhy, LNhz);
    ROTATEZ(phiRef, LNhx, LNhy, LNhz);
    ROTATEY(inclIn, LNhx, LNhy, LNhz);
    printf("LNh: %12.4e  %12.4e  %12.4e  norm: %12.4e\n",LNhx,LNhy,LNhz,sqrt(LNhx*LNhx+LNhy*LNhy+LNhz*LNhz));
    printf("     %12.4e  %12.4e  %12.4e\n",sin(*inc),0.,cos(*inc));
    printf("S1.L     %12.4e  %12.4e\n",(*S1x)*LNhx+(*S1y)*LNhy+(*S1z)*LNhz,S1xIn*LNhx+S1yIn*LNhy+S1zIn*LNhz);
    printf("S2.L     %12.4e  %12.4e\n",(*S2x)*LNhx+(*S2y)*LNhy+(*S2z)*LNhz,S2xIn*LNhx+S2yIn*LNhy+S2zIn*LNhz);
    printf("S1.S2    %12.4e  %12.4e\n",(*S1x)*(*S2x)+(*S1y)*(*S2y)+(*S1z)*(*S2z),S1xIn*S2xIn+S1yIn*S2yIn+S1zIn*S2zIn);
    printf("***********************************************\n");*/
    break;

  case LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L:
    /* FRAME_AXIS_ORBITAL_L (default)
     * (spins wrt to L//z and x// body 2->1,
     * inclIn is the angle between L and N: if
     * LNhat=(0,0,1) N=(0,sin(inclIn),cos(inclIn)) in the Y-z plane,
     * where x is rotated by phiRef wrt X around z.
     */
  default:
    //We need to set N into the (0,0,1) axis and L into (sin(inc),0,cos(inc))
    //We can first rotate by phiRef-pi/2 around the angular momentum,
    //which brings the line of sight into the new xz plane and leaves the ang momentum as the z axis
    ROTATEZ(phiRef-LAL_PI/2.,*S1x,*S1y,*S1z);
    ROTATEZ(phiRef-LAL_PI/2.,*S2x,*S2y,*S2z);
    //We can then rotate by -inclIn around the y axis to bring the line of sight to the z axis.
    ROTATEY(-inclIn,*S1x,*S1y,*S1z);
    ROTATEY(-inclIn,*S2x,*S2y,*S2z);
    //The ang momentum is now in the xz plane but with negative projection on x
    //Now we rotate by Pi around z to bring L to the positive x half plane
    ROTATEZ(LAL_PI,*S1x,*S1y,*S1z);
    ROTATEZ(LAL_PI,*S2x,*S2y,*S2z);
    *inc=inclIn;
    //For a check uncomment the following lines
    //Anyway there is an automatic check at build in
    // lalsimulation/test/InitialSpinRotationTest
    /*printf("****************************************************\n");
    printf("** Check of LAL_SIM_INSPIRAL_FRAME_AXIS_ORBITAL_L **\n");
    printf("****************************************************\n");
    Nhx=0.;
    Nhy=sin(inclIn);
    Nhz=cos(inclIn);
    LNhx=0.;
    LNhy=0.;
    LNhz=1.;
    REAL8 nnx=1.;
    REAL8 nny=0.;
    REAL8 nnz=0.;
    REAL8 llx=0.;
    REAL8 lly=1.;
    REAL8 llz=0.;
    ROTATEZ(-LAL_PI/2.,Nhx,Nhy,Nhz);
    ROTATEZ(phiRef-LAL_PI/2.,nnx,nny,nnz);
    ROTATEZ(phiRef-LAL_PI/2.,llx,lly,llz);
    ROTATEY(-inclIn,Nhx,Nhy,Nhz);
    ROTATEY(-inclIn,LNhx,LNhy,LNhz);
    ROTATEY(-inclIn,nnx,nny,nnz);
    ROTATEY(-inclIn,llx,lly,llz);
    ROTATEZ(LAL_PI,Nhx,Nhy,Nhz);
    ROTATEZ(LAL_PI,LNhx,LNhy,LNhz);
    ROTATEZ(LAL_PI,nnx,nny,nnz);
    ROTATEZ(LAL_PI,llx,lly,llz);
    printf("N:   %12.4e  %12.4e  %12.4e  norm: %12.4e\n",Nhx,Nhy,Nhz,sqrt(Nhx*Nhx+Nhy*Nhy+Nhz*Nhz));
    printf("     0.   0.  1.\n");
    printf("LNh: %12.4e  %12.4e  %12.4e  norm: %12.4e\n",LNhx,LNhy,LNhz,sqrt(LNhx*LNhx+LNhy*LNhy+LNhz*LNhz));
    printf("     %12.4e  %12.4e  %12.4e\n",sin(*inc),0.,cos(*inc));
    printf("S1.L     %12.4e  %12.4e\n",(*S1x)*LNhx+(*S1y)*LNhy+(*S1z)*LNhz,S1zIn);
    printf("S2.L     %12.4e  %12.4e\n",(*S2x)*LNhx+(*S2y)*LNhy+(*S2z)*LNhz,S2zIn);
    printf("S1.S2    %12.4e  %12.4e\n",(*S1x)*(*S2x)+(*S1y)*(*S2y)+(*S1z)*(*S2z),S1xIn*S2xIn+S1yIn*S2yIn+S1zIn*S2zIn);
    printf("S1.n    %12.4e  %12.4e\n",S1xIn,(*S1x)*nnx+(*S1y)*nny+(*S1z)*nnz);
    printf("S2.n    %12.4e  %12.4e\n",S2xIn,(*S2x)*nnx+(*S2y)*nny+(*S2z)*nnz);
    printf("S1.l    %12.4e  %12.4e\n",S1yIn,(*S1x)*llx+(*S1y)*lly+(*S1z)*llz);
    printf("S2.l    %12.4e  %12.4e\n",S2yIn,(*S2x)*llx+(*S2y)*lly+(*S2z)*llz);
    printf("L.n     0.      %12.4e\n",nnx*sin(*inc)+nnz*cos(*inc));
    printf("L.l     0.      %12.4e\n",llx*sin(*inc)+llz*cos(*inc));
    printf("N.n     %12.4e  %12.4e\n",sin(inclIn)*sin(phiRef),Nhx*nnx+Nhy*nny+Nhz*nnz);
    printf("N.l     %12.4e  %12.4e\n",sin(inclIn)*cos(phiRef),Nhx*llx+Nhy*lly+Nhz*llz);
    printf("n.l     0.  %12.4e\n",nnx*llx+nny*lly+nnz*llz);
    printf("****************************************************\n");*/
    break;
  }

  return XLAL_SUCCESS;
} // End of XLALSimInspiralInitialConditionsPrecessingApproxs()

/**
 * This function evolves the orbital equations for a precessing binary using
 * the \"TaylorT1/T5/T4\" approximant for solving the orbital dynamics
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
 * NOTE: All vectors are given in the frame
 * where the z-axis is set by the angular momentum at reference frequency,
 * the x-axis is chosen orthogonal to it, and the y-axis is given by the RH rule.
 * Initial values must be passed in this frame, and the time series of the
 * vector components will also be returned in this frame.
 *
 * Review completed on git hash ...
 *
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
	const REAL8 deltaT,   	        /**< sampling interval (s) */
	const REAL8 m1_SI,     	        /**< mass of companion 1 (kg) */
	const REAL8 m2_SI,     	        /**< mass of companion 2 (kg) */
	const REAL8 fStart,             /**< starting GW frequency */
	const REAL8 fEnd,               /**< ending GW frequency, fEnd=0 means integrate as far forward as possible */
	const REAL8 s1x,                /**< initial value of S1x */
	const REAL8 s1y,                /**< initial value of S1y */
	const REAL8 s1z,                /**< initial value of S1z */
	const REAL8 s2x,                /**< initial value of S2x */
	const REAL8 s2y,                /**< initial value of S2y */
	const REAL8 s2z,                /**< initial value of S2z */
	const REAL8 lnhatx,             /**< initial value of LNhatx */
	const REAL8 lnhaty,             /**< initial value of LNhaty */
	const REAL8 lnhatz,             /**< initial value of LNhatz */
	const REAL8 e1x,                /**< initial value of E1x */
	const REAL8 e1y,                /**< initial value of E1y */
	const REAL8 e1z,                /**< initial value of E1z */
	const REAL8 lambda1,            /**< (tidal deformability of mass 1) / (mass of body 1)^5 (dimensionless) */
	const REAL8 lambda2,            /**< (tidal deformability of mass 2) / (mass of body 2)^5 (dimensionless) */
	const REAL8 quadparam1,         /**< phenom. parameter describing induced quad. moment of body 1 (=1 for BHs, ~2-12 for NSs) */
	const REAL8 quadparam2,         /**< phenom. parameter describing induced quad. moment of body 2 (=1 for BHs, ~2-12 for NSs) */
	const LALSimInspiralSpinOrder spinO,  /**< twice PN order of spin effects */
	const LALSimInspiralTidalOrder tideO, /**< twice PN order of tidal effects */
	const INT4 phaseO,                    /**< twice post-Newtonian order */
	const INT4 lscorr,                    /**< flag to control L_S terms */
	const Approximant approx              /**< PN approximant (SpinTaylorT1/T5/T4) */
	)
{
    INT4 intreturn;
    LALAdaptiveRungeKuttaIntegrator *integrator = NULL;     /* GSL integrator object */
    REAL8 yinit[LAL_NUM_ST4_VARIABLES];       /* initial values of parameters */
    REAL8Array *yout;	 /* time series of variables returned from integrator */
    /* intermediate variables */
    UINT4 i, cutlen, len;
    int sgn, offset;
    REAL8 norm1, norm2, dtStart, dtEnd, lengths, wEnd, m1sec, m2sec, Msec, Mcsec, fTerm;
    LIGOTimeGPS tStart = LIGOTIMEGPSZERO;

    if ( !V || !Phi || !S1x || !S1y || !S1z || !S2x || !S2y || !S2z
            || !LNhatx || !LNhaty || !LNhatz || !E1x || !E1y || !E1z )
    {
        XLALPrintError("XLAL Error - %s: NULL(s) in output parameters\n",
                       __func__);
        XLAL_ERROR(XLAL_EINVAL);
    }

    if ( spinO>=7 ) {
      XLALPrintError("XLAL Error - %s: Averaged-orbit dynamics incompatible with spin order %d chosen.\n",__func__, spinO);
      XLAL_ERROR(XLAL_EINVAL);
    }

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

    // Fill params struct with values of constant coefficients of the model
    XLALSimInspiralSpinTaylorTxCoeffs *params=NULL;
    if( approx == SpinTaylorT4 )
    {
        XLALSimInspiralSpinTaylorT4Setup(&params, m1_SI, m2_SI, fStart, fEnd,
					 lambda1, lambda2, quadparam1, quadparam2, spinO, tideO, phaseO, lscorr, 0);
      
    }
    else if( approx == SpinTaylorT5 )
    {
        XLALSimInspiralSpinTaylorT5Setup(&params, m1_SI, m2_SI, fStart, fEnd,
					 lambda1, lambda2, quadparam1, quadparam2, spinO, tideO, phaseO, lscorr);
    }
    else if( approx == SpinTaylorT1 )
    {
        XLALSimInspiralSpinTaylorT1Setup(&params, m1_SI, m2_SI, fStart, fEnd,
					 lambda1, lambda2, quadparam1, quadparam2, spinO, tideO, phaseO, lscorr);
    }
    else
    {
        XLALPrintError("XLAL Error - %s: Approximant must be one of SpinTaylorT1, SpinTaylorT5, SpinTaylorT4, but %i provided\n",
                __func__, approx);
        XLAL_ERROR(XLAL_EINVAL);

    }
    m1sec = m1_SI / LAL_MSUN_SI * LAL_MTSUN_SI;
    m2sec = m2_SI / LAL_MSUN_SI * LAL_MTSUN_SI;
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
    /* S1(x,y,z): LAL normalizes spins by total mass, not individual mass */
    norm1 = m1sec * m1sec / Msec / Msec;
    yinit[5] = norm1 * s1x;
    yinit[6] = norm1 * s1y;
    yinit[7] = norm1 * s1z;
    /* S2(x,y,z) */
    norm2 = m2sec * m2sec / Msec / Msec;
    yinit[8] = norm2 * s2x;
    yinit[9] = norm2 * s2y;
    yinit[10]= norm2 * s2z;
    /* E1(x,y,z) */
    yinit[11] = e1x;
    yinit[12] = e1y;
    yinit[13] = e1z;

    /* initialize the integrator */
    if( approx == SpinTaylorT4 )
      integrator = XLALAdaptiveRungeKutta4Init(LAL_NUM_ST4_VARIABLES,
					       XLALSimInspiralSpinTaylorT4DerivativesAvg,
					       XLALSimInspiralSpinTaylorStoppingTest,
					       LAL_ST4_ABSOLUTE_TOLERANCE, LAL_ST4_RELATIVE_TOLERANCE);
    else if( approx == SpinTaylorT5 )
      integrator = XLALAdaptiveRungeKutta4Init(LAL_NUM_ST4_VARIABLES,
					       XLALSimInspiralSpinTaylorT5DerivativesAvg,
					       XLALSimInspiralSpinTaylorStoppingTest,
					       LAL_ST4_ABSOLUTE_TOLERANCE, LAL_ST4_RELATIVE_TOLERANCE);
    else if( approx == SpinTaylorT1 )
      integrator = XLALAdaptiveRungeKutta4Init(LAL_NUM_ST4_VARIABLES,
					       XLALSimInspiralSpinTaylorT1DerivativesAvg,
					       XLALSimInspiralSpinTaylorStoppingTest,
					       LAL_ST4_ABSOLUTE_TOLERANCE, LAL_ST4_RELATIVE_TOLERANCE);
    else
    {
        XLALPrintError("XLAL Error - %s: Approximant must be one of SpinTaylorT1, SpinTaylorT5, SpinTaylorT4, but %i provided\n",
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
    len = XLALAdaptiveRungeKutta4Hermite(integrator, params, yinit,
            0.0, lengths/Msec, sgn*deltaT/Msec, &yout);

    if (!yout)
    {
        XLALPrintError("XLAL Error - %s: integration failed (yout == NULL)\n",
                       __func__);
        XLAL_ERROR(XLAL_EFUNC);
    }

    intreturn = integrator->returncode;
    XLALAdaptiveRungeKuttaFree(integrator);
    LALFree(params);

    if (!len) 
    {
        XLALPrintError("XLAL Error - %s: integration failed with errorcode %d.\n", __func__, intreturn);
        XLAL_ERROR(XLAL_EFUNC);
    }

    /* Print warning about abnormal termination */
    if (intreturn != 0 && intreturn != LALSIMINSPIRAL_ST_TEST_ENERGY
            && intreturn != LALSIMINSPIRAL_ST_TEST_FREQBOUND)
    {
        XLALPrintWarning("XLAL Warning - %s: integration terminated with code %d.\n Waveform parameters were m1 = %e, m2 = %e, s1 = (%e,%e,%e), s2 = (%e,%e,%e), inc = %e.\n", __func__, intreturn, m1_SI / LAL_MSUN_SI, m2_SI / LAL_MSUN_SI, s1x, s1y, s1z, s2x, s2y, s2z, acos(lnhatz));
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
    if ( !*V || !*Phi || !*S1x || !*S1y || !*S1z || !*S2x || !*S2y || !*S2z
             || !*LNhatx || !*LNhaty || !*LNhatz || !*E1x || !*E1y || !*E1z )
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
	/* Spins are returned in the standard convention */
        (*S1x)->data->data[j] 		= yout->data[6*len+i]/norm1;
        (*S1y)->data->data[j] 		= yout->data[7*len+i]/norm1;
        (*S1z)->data->data[j] 		= yout->data[8*len+i]/norm1;
        (*S2x)->data->data[j] 		= yout->data[9*len+i]/norm2;
        (*S2y)->data->data[j] 		= yout->data[10*len+i]/norm2;
        (*S2z)->data->data[j] 		= yout->data[11*len+i]/norm2;
        (*E1x)->data->data[j] 		= yout->data[12*len+i];
        (*E1y)->data->data[j] 		= yout->data[13*len+i];
        (*E1z)->data->data[j] 		= yout->data[14*len+i];
    }

    XLALDestroyREAL8Array(yout);

    return XLAL_SUCCESS;
}

/**
 * This function evolves the orbital equations for a precessing binary using
 * the \"TaylorT1/T5/T4\" approximant for solving the orbital dynamics
 * (see arXiv:0907.0700 for a review of the various PN approximants).
 *
 * It returns the final values of the \"orbital velocity\", orbital phase,
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
 * NOTE: All vectors are given in the frame
 * where the z-axis is set by the angular momentum at reference frequency,
 * the x-axis is chosen orthogonal to it, and the y-axis is given by the RH rule.
 * Initial values must be passed in this frame, and the time series of the
 * vector components will also be returned in this frame.
 *
 * Review completed on git hash ...
 *
 */
int XLALSimInspiralSpinTaylorPNEvolveOrbitOnlyFinal(
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
	const REAL8 deltaT,   	        /**< sampling interval (s) */
	const REAL8 m1_SI,     	        /**< mass of companion 1 (kg) */
	const REAL8 m2_SI,     	        /**< mass of companion 2 (kg) */
	const REAL8 fStart,             /**< starting GW frequency */
	const REAL8 fEnd,               /**< ending GW frequency, fEnd=0 means integrate as far forward as possible */
	const REAL8 s1x,                /**< initial value of S1x */
	const REAL8 s1y,                /**< initial value of S1y */
	const REAL8 s1z,                /**< initial value of S1z */
	const REAL8 s2x,                /**< initial value of S2x */
	const REAL8 s2y,                /**< initial value of S2y */
	const REAL8 s2z,                /**< initial value of S2z */
	const REAL8 lnhatx,             /**< initial value of LNhatx */
	const REAL8 lnhaty,             /**< initial value of LNhaty */
	const REAL8 lnhatz,             /**< initial value of LNhatz */
	const REAL8 e1x,                /**< initial value of E1x */
	const REAL8 e1y,                /**< initial value of E1y */
	const REAL8 e1z,                /**< initial value of E1z */
	const REAL8 lambda1,            /**< (tidal deformability of mass 1) / (mass of body 1)^5 (dimensionless) */
	const REAL8 lambda2,            /**< (tidal deformability of mass 2) / (mass of body 2)^5 (dimensionless) */
	const REAL8 quadparam1,         /**< phenom. parameter describing induced quad. moment of body 1 (=1 for BHs, ~2-12 for NSs) */
	const REAL8 quadparam2,         /**< phenom. parameter describing induced quad. moment of body 2 (=1 for BHs, ~2-12 for NSs) */
	const LALSimInspiralSpinOrder spinO,  /**< twice PN order of spin effects */
	const LALSimInspiralTidalOrder tideO, /**< twice PN order of tidal effects */
	const INT4 phaseO,                    /**< twice post-Newtonian order */
	const INT4 lscorr,                    /**< flag to control L_S terms */
	const Approximant approx              /**< PN approximant (SpinTaylorT1/T5/T4) */
	)
{
    INT4 intreturn;
    LALAdaptiveRungeKuttaIntegrator *integrator = NULL;     /* GSL integrator object */
    REAL8 yinit[LAL_NUM_ST4_VARIABLES];       /* initial values of parameters */
    /* intermediate variables */
    UINT4 out;
    int sgn;
    REAL8 norm1, norm2, dtStart, dtEnd, lengths, m1sec, m2sec, Msec, Mcsec;
    LIGOTimeGPS tStart = LIGOTIMEGPSZERO;

    if ( !V || !Phi || !S1x || !S1y || !S1z || !S2x || !S2y || !S2z
            || !LNhatx || !LNhaty || !LNhatz || !E1x || !E1y || !E1z )
    {
        XLALPrintError("XLAL Error - %s: NULL(s) in output parameters\n",
                       __func__);
        XLAL_ERROR(XLAL_EINVAL);
    }

    if ( spinO>=7 ) {
      XLALPrintError("XLAL Error - %s: Averaged-orbit dynamics incompatible with spin order %d chosen.\n",__func__, spinO);
      XLAL_ERROR(XLAL_EINVAL);
    }

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

    // Fill params struct with values of constant coefficients of the model
    XLALSimInspiralSpinTaylorTxCoeffs *params=NULL;
    if( approx == SpinTaylorT4 )
    {
        XLALSimInspiralSpinTaylorT4Setup(&params, m1_SI, m2_SI, fStart, fEnd,
					 lambda1, lambda2, quadparam1, quadparam2, spinO, tideO, phaseO, lscorr, 0);
      
    }
    else if( approx == SpinTaylorT5 )
    {
        XLALSimInspiralSpinTaylorT5Setup(&params, m1_SI, m2_SI, fStart, fEnd,
					 lambda1, lambda2, quadparam1, quadparam2, spinO, tideO, phaseO, lscorr);
    }
    else if( approx == SpinTaylorT1 )
    {
        XLALSimInspiralSpinTaylorT1Setup(&params, m1_SI, m2_SI, fStart, fEnd,
					 lambda1, lambda2, quadparam1, quadparam2, spinO, tideO, phaseO, lscorr);
    }
    else
    {
        XLALPrintError("XLAL Error - %s: Approximant must be one of SpinTaylorT1, SpinTaylorT5, SpinTaylorT4, but %i provided\n",
                __func__, approx);
        XLAL_ERROR(XLAL_EINVAL);

    }
    m1sec = m1_SI / LAL_MSUN_SI * LAL_MTSUN_SI;
    m2sec = m2_SI / LAL_MSUN_SI * LAL_MTSUN_SI;
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
    /* S1(x,y,z): LAL normalizes spins by total mass, not individual mass */
    norm1 = m1sec * m1sec / Msec / Msec;
    yinit[5] = norm1 * s1x;
    yinit[6] = norm1 * s1y;
    yinit[7] = norm1 * s1z;
    /* S2(x,y,z) */
    norm2 = m2sec * m2sec / Msec / Msec;
    yinit[8] = norm2 * s2x;
    yinit[9] = norm2 * s2y;
    yinit[10]= norm2 * s2z;
    /* E1(x,y,z) */
    yinit[11] = e1x;
    yinit[12] = e1y;
    yinit[13] = e1z;

    /* initialize the integrator */
    if( approx == SpinTaylorT4 )
      integrator = XLALAdaptiveRungeKutta4Init(LAL_NUM_ST4_VARIABLES,
					       XLALSimInspiralSpinTaylorT4DerivativesAvg,
					       XLALSimInspiralSpinTaylorStoppingTest,
					       LAL_ST4_ABSOLUTE_TOLERANCE, LAL_ST4_RELATIVE_TOLERANCE);
    else if( approx == SpinTaylorT5 )
      integrator = XLALAdaptiveRungeKutta4Init(LAL_NUM_ST4_VARIABLES,
					       XLALSimInspiralSpinTaylorT5DerivativesAvg,
					       XLALSimInspiralSpinTaylorStoppingTest,
					       LAL_ST4_ABSOLUTE_TOLERANCE, LAL_ST4_RELATIVE_TOLERANCE);
    else if( approx == SpinTaylorT1 )
      integrator = XLALAdaptiveRungeKutta4Init(LAL_NUM_ST4_VARIABLES,
					       XLALSimInspiralSpinTaylorT1DerivativesAvg,
					       XLALSimInspiralSpinTaylorStoppingTest,
					       LAL_ST4_ABSOLUTE_TOLERANCE, LAL_ST4_RELATIVE_TOLERANCE);
    else
    {
        XLALPrintError("XLAL Error - %s: Approximant must be one of SpinTaylorT1, SpinTaylorT5, SpinTaylorT4, but %i provided\n",
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
    out = XLALAdaptiveRungeKutta4HermiteOnlyFinal(integrator, params, yinit,
						  0.0, lengths/Msec, LAL_PI * Msec * fEnd, sgn*deltaT/Msec);

    intreturn = integrator->returncode;
    XLALAdaptiveRungeKuttaFree(integrator);
    LALFree(params);

    if (!out) 
    {
        XLALPrintError("XLAL Error - %s: integration failed with errorcode %d.\n", __func__, intreturn);
        XLAL_ERROR(XLAL_EFUNC);
    }

    /* Raise error about abnormal termination */
    if (intreturn != 0 && intreturn != LALSIMINSPIRAL_ST_TEST_ENERGY
            && intreturn != LALSIMINSPIRAL_ST_TEST_FREQBOUND)
    {
        XLALPrintError("XLAL Error - %s: integration terminated with code %d.\n Waveform parameters were m1 = %e, m2 = %e, s1 = (%e,%e,%e), s2 = (%e,%e,%e), inc = %e.\n", __func__, intreturn, m1_SI / LAL_MSUN_SI, m2_SI / LAL_MSUN_SI, s1x, s1y, s1z, s2x, s2y, s2z, acos(lnhatz));
	XLAL_ERROR(XLAL_EFUNC);
    }

    /* allocate memory for output vectors */
    *V = XLALCreateREAL8TimeSeries( "PN_EXPANSION_PARAMETER", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, 1); 
    *Phi = XLALCreateREAL8TimeSeries( "ORBITAL_PHASE", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, 1); 
    *S1x = XLALCreateREAL8TimeSeries( "SPIN1_X_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, 1); 
    *S1y = XLALCreateREAL8TimeSeries( "SPIN1_Y_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, 1); 
    *S1z = XLALCreateREAL8TimeSeries( "SPIN1_Z_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, 1); 
    *S2x = XLALCreateREAL8TimeSeries( "SPIN2_X_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, 1); 
    *S2y = XLALCreateREAL8TimeSeries( "SPIN2_Y_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, 1); 
    *S2z = XLALCreateREAL8TimeSeries( "SPIN2_Z_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, 1); 
    *LNhatx = XLALCreateREAL8TimeSeries( "LNHAT_X_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, 1); 
    *LNhaty = XLALCreateREAL8TimeSeries( "LNHAT_Y_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, 1); 
    *LNhatz = XLALCreateREAL8TimeSeries( "LNHAT_Z_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, 1); 
    *E1x = XLALCreateREAL8TimeSeries( "E1_BASIS_X_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, 1); 
    *E1y = XLALCreateREAL8TimeSeries( "E1_BASIS_Y_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, 1); 
    *E1z = XLALCreateREAL8TimeSeries( "E1_BASIS_Z_COMPONENT", &tStart, 0., 
            deltaT, &lalDimensionlessUnit, 1); 
    if ( !*V || !*Phi || !*S1x || !*S1y || !*S1z || !*S2x || !*S2y || !*S2z
             || !*LNhatx || !*LNhaty || !*LNhatz || !*E1x || !*E1y || !*E1z )
    {
        XLAL_ERROR(XLAL_EFUNC);
    }

    /* Copy dynamical variables from yinit array to output */

    (*Phi)->data->data[0] = yinit[0];
    (*V)->data->data[0] = cbrt(yinit[1]);
    (*LNhatx)->data->data[0] = yinit[2];
    (*LNhaty)->data->data[0] = yinit[3];
    (*LNhatz)->data->data[0] = yinit[4];
    (*S1x)->data->data[0] = yinit[5]/norm1;
    (*S1y)->data->data[0] = yinit[6]/norm1;
    (*S1z)->data->data[0] = yinit[7]/norm1;
    (*S2x)->data->data[0] = yinit[8]/norm2;
    (*S2y)->data->data[0] = yinit[9]/norm2;
    (*S2z)->data->data[0] = yinit[10]/norm2;
    (*E1x)->data->data[0] = yinit[11];
    (*E1y)->data->data[0] = yinit[12];
    (*E1z)->data->data[0] = yinit[13];

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
 *
 * REVIEW completed on git hash 6640e79e60791d5230731acc63351676ce7ce413
 *
 */
int XLALSimInspiralSpinTaylorT4OLD(
	REAL8TimeSeries **hplus,        /**< +-polarization waveform */
	REAL8TimeSeries **hcross,       /**< x-polarization waveform */
	REAL8 phiRef,                   /**< orbital phase at reference pt. */
	REAL8 UNUSED v0,                       /**< tail gauge term (default = 1) */
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
	LALDict *LALparams,             /**< LAL dictionary containing accessory parameters */
	int phaseO,                     /**< twice PN phase order */
	int amplitudeO                  /**< twice PN amplitude order */
	)
{
    Approximant approx = SpinTaylorT4;
    int n=XLALSimInspiralWaveformParamsInsertTidalLambda1(LALparams, lambda1);
    n=XLALSimInspiralWaveformParamsInsertTidalLambda2(LALparams, lambda2);
    n=XLALSimInspiralWaveformParamsInsertdQuadMon1(LALparams, quadparam1);
    n=XLALSimInspiralWaveformParamsInsertdQuadMon2(LALparams, quadparam2);
    n=XLALSimInspiralWaveformParamsInsertPNPhaseOrder(LALparams, phaseO);
    n=XLALSimInspiralWaveformParamsInsertPNAmplitudeOrder(LALparams, amplitudeO);
    n = XLALSimInspiralSpinTaylorDriver(hplus, hcross, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
	    phiRef, deltaT,
            m1, m2, fStart, fRef, r, s1x, s1y, s1z, s2x, s2y, s2z,
            lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, LALparams, approx);

    return n;
}

int XLALSimInspiralSpinTaylorT4(
	REAL8TimeSeries **hplus,        /**< +-polarization waveform */
	REAL8TimeSeries **hcross,       /**< x-polarization waveform */
	REAL8 phiRef,                   /**< orbital phase at reference pt. */
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
	LALDict *LALparams              /**< LAL dictionary containing accessory parameters */
	)
{
    Approximant approx = SpinTaylorT4;
    int n = XLALSimInspiralSpinTaylorDriver(hplus, hcross, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
	    phiRef, deltaT,
            m1, m2, fStart, fRef, r, s1x, s1y, s1z, s2x, s2y, s2z,
            lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, LALparams, approx);

    return n;
}

/**
 * Driver routine to compute a precessing post-Newtonian inspiral waveform
 * with phasing computed from energy balance using the so-called \"T1\" method.
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
 * REVIEW completed on git hash 6640e79e60791d5230731acc63351676ce7ce413
 *
 */
int XLALSimInspiralSpinTaylorT1OLD(
	REAL8TimeSeries **hplus,        /**< +-polarization waveform */
	REAL8TimeSeries **hcross,       /**< x-polarization waveform */
	REAL8 phiRef,                   /**< orbital phase at reference pt. */
	REAL8 UNUSED v0,                       /**< tail gauge term (default = 1) */
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
	LALDict *LALparams,             /**< LAL dictionary containing accessory parameters */
	int phaseO,                     /**< twice PN phase order */
	int amplitudeO                  /**< twice PN amplitude order */
	)
{
    Approximant approx = SpinTaylorT1;
    int n=XLALSimInspiralWaveformParamsInsertTidalLambda1(LALparams, lambda1);
    n=XLALSimInspiralWaveformParamsInsertTidalLambda2(LALparams, lambda2);
    n=XLALSimInspiralWaveformParamsInsertdQuadMon1(LALparams, quadparam1);
    n=XLALSimInspiralWaveformParamsInsertdQuadMon2(LALparams, quadparam2);
    n=XLALSimInspiralWaveformParamsInsertPNPhaseOrder(LALparams, phaseO);
    n=XLALSimInspiralWaveformParamsInsertPNAmplitudeOrder(LALparams, amplitudeO);
    n = XLALSimInspiralSpinTaylorDriver(hplus, hcross, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
	    phiRef, deltaT,
            m1, m2, fStart, fRef, r, s1x, s1y, s1z, s2x, s2y, s2z,
            lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, LALparams, approx);

    return n;
}

int XLALSimInspiralSpinTaylorT1(
	REAL8TimeSeries **hplus,        /**< +-polarization waveform */
	REAL8TimeSeries **hcross,       /**< x-polarization waveform */
	REAL8 phiRef,                   /**< orbital phase at reference pt. */
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
	LALDict *LALparams              /**< LAL dictionary containing accessory parameters */
	)
{
    Approximant approx = SpinTaylorT1;
    int n = XLALSimInspiralSpinTaylorDriver(hplus, hcross, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
	    phiRef, deltaT,
            m1, m2, fStart, fRef, r, s1x, s1y, s1z, s2x, s2y, s2z,
            lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, LALparams, approx);

    return n;
}

/**
 * Driver routine to compute a precessing post-Newtonian inspiral waveform
 * with phasing computed from energy balance using the so-called \"T5\" method.
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
 *
 * REVIEW completed on git hash ...
 *
 */
int XLALSimInspiralSpinTaylorT5(
	REAL8TimeSeries **hplus,        /**< +-polarization waveform */
	REAL8TimeSeries **hcross,       /**< x-polarization waveform */
	REAL8 phiRef,                   /**< orbital phase at reference pt. */
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
	LALDict *LALparams             /**< LAL dictionary containing accessory parameters */
	)
{
    Approximant approx = SpinTaylorT5;
    int n = XLALSimInspiralSpinTaylorDriver(hplus, hcross, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
	    phiRef, deltaT,
            m1, m2, fStart, fRef, r, s1x, s1y, s1z, s2x, s2y, s2z,
            lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, LALparams, approx);

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
 *
 * !!!UNREVIEWED!!!
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
 *
 * !!!UNREVIEWED!!!
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
	    quadparam1, quadparam2, spinO, tideO, phaseO, 0, approx);
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
 * !!!UNREVIEWED!!!
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
	LALDict *LALparams,             /**< LAL dictionary containing accessory parameters */
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
            quadparam1, quadparam2, LALparams, phaseO, amplitudeO, approx, phiRefAtEnd);

    return n;
}

/**
 * Driver routine to compute a precessing post-Newtonian inspiral waveform in the Fourier domain
 * with phasing computed from energy balance using the so-called \"T2\" method
 * see arXiv: 0907.0700 for its definition,
 * but in its \"T5\" variant, see arXiv: 1107.1267.
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
 * !!!UNREVIEWED!!!
 */
int XLALSimInspiralSpinTaylorT5Fourier(
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
	LALDict *LALparams,             /**< LAL dictionary containing accessory parameters */
        INT4 phaseO,                     /**< twice PN phase order */
        INT4 amplitudeO,                /**< twice PN amplitude order */
        INT4 phiRefAtEnd                /**< whether phiRef corresponds to the end of the inspiral */
        )
{
    Approximant approx = SpinTaylorT5;
    int n = XLALSimInspiralSpinTaylorDriverFourier(hplus, hcross, fMin,
            fMax, deltaF, kMax, phiRef, v0,
            m1, m2, fStart, fRef, r, s1x, s1y, s1z, s2x, s2y, s2z,
            lnhatx, lnhaty, lnhatz, e1x, e1y, e1z, lambda1, lambda2,
            quadparam1, quadparam2, LALparams, phaseO, amplitudeO, approx, phiRefAtEnd);

    return n;
}

/** @} */
