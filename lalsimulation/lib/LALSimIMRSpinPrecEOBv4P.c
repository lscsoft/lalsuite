/*
 *  Copyright (C) 2017-2019 Sylvain Marsat, Serguei Ossokine, Roberto Cotesta
 *                2016-2017 Stas Babak, Andrea Taracchini (Precessing EOB)
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
#ifndef _LALSIMIMRSPINPRECEOBv4P_C
#define _LALSIMIMRSPINPRECEOBv4P_C
/**
 * @addtogroup LALSimIMRSpinPrecEOBv4P_c
 *
 * @author Sylvain Marsat, Serguei Ossokine, Roberto Cotesta, Stas Babak, Andrea
 * Taracchini
 *
 * \brief Functions for producing SEOBNRv4P(HM) waveforms for
 * precessing binaries of spinning compact objects.
 */

// clang-format off
#include <math.h>
#include <complex.h>
#include <gsl/gsl_deriv.h>

#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/Date.h>
#include <lal/TimeSeries.h>
#include <lal/LALAdaptiveRungeKuttaIntegrator.h>
#include <lal/SphericalHarmonics.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include <lal/Units.h>
#include <lal/VectorOps.h>
#include "LALSimIMREOBNRv2.h"
#include "LALSimIMRSpinEOB.h"
#include "LALSimInspiralPrecess.h"
#include "LALSimBlackHoleRingdownPrec.h"
#include "LALSimFindAttachTime.h"

// clang-format on

/* Include all the static function files we need */
#include "LALSimIMREOBHybridRingdownPrec.c"
#include "LALSimIMREOBNQCCorrection.c"
#include "LALSimIMREOBNewtonianMultipole.c"
#include "LALSimIMRSpinAlignedEOBHcapDerivative.c"
#include "LALSimIMRSpinEOBAuxFuncs.c"
#include "LALSimIMRSpinEOBFactorizedFluxPrec.c"
#include "LALSimIMRSpinEOBFactorizedWaveformCoefficientsPrec.c"
#include "LALSimIMRSpinEOBFactorizedWaveformPrec.c"
#include "LALSimIMRSpinEOBHamiltonian.c"
#include "LALSimIMRSpinEOBHamiltonianPrec.c"
#include "LALSimIMRSpinEOBHcapNumericalDerivativePrec.c"
#include "LALSimIMRSpinEOBInitialConditions.c"
#include "LALSimIMRSpinEOBInitialConditionsPrec.c"
#include "LALSimIMRSpinPrecEOBEulerAngles.c"
#include "LALSimInspiraldEnergyFlux.c"

/* Begin OPTv3 */
//#include "LALSimIMRSpinPrecEOBGSLOptimizedInterpolation.c"
#include "LALSpinPrecHcapRvecDerivative_v3opt.c"
//#include "LALSimIMRSpinPrecEOBWfGen.c"
/* End OPTv3 */

#define debugOutput 0

#ifdef __GNUC__
#ifndef UNUSED
#define UNUSED __attribute__ ((unused))
#endif
#else
#define UNUSED
#endif

/* Maximum l allowed in the SEOB model -- input ModeArray will only be scanned
 * up to this value of l */
#define _SEOB_MODES_LMAX 5

/* Version flag used in XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients */
#define v4Pwave 451

/* Structure containing the approximant name and its number in LALSimInspiral.c
 */
struct approximant {
  const char *name;
  UINT4 number;
};
struct approximant v4P = {.name = "SEOBNRv4P", .number = 401};
struct approximant v4PHM = {.name = "SEOBNRv4PHM", .number = 402};

/* Number of dynamics variables stored by v4P */
#define v4PdynamicsVariables 26

#define FREE_ALL                                                               \
  if (ICvalues != NULL)                                                        \
    XLALDestroyREAL8Vector(ICvalues);                                          \
  if (dynamicsAdaS != NULL)                                                    \
    XLALDestroyREAL8Array(dynamicsAdaS);                                       \
  if (seobdynamicsAdaS != NULL)                                                \
    SEOBdynamics_Destroy(seobdynamicsAdaS);                                    \
  if (seobvalues_tstartHiS != NULL)                                            \
    XLALDestroyREAL8Vector(seobvalues_tstartHiS);                              \
  if (ICvaluesHiS != NULL)                                                     \
    XLALDestroyREAL8Vector(ICvaluesHiS);                                       \
  if (dynamicsHiS != NULL)                                                     \
    XLALDestroyREAL8Array(dynamicsHiS);                                        \
  if (chi1L_tPeakOmega != NULL)                                                \
    XLALDestroyREAL8Vector(chi1L_tPeakOmega);                                  \
  if (chi2L_tPeakOmega != NULL)                                                \
    XLALDestroyREAL8Vector(chi2L_tPeakOmega);                                  \
  if (seobdynamicsHiS != NULL)                                                 \
    SEOBdynamics_Destroy(seobdynamicsHiS);                                     \
  if (seobvalues_tPeakOmega != NULL)                                           \
    XLALDestroyREAL8Vector(seobvalues_tPeakOmega);                             \
  if (Jfinal != NULL)                                                          \
    XLALDestroyREAL8Vector(Jfinal);                                            \
  if (listhPlm_HiS != NULL)                                                    \
    SphHarmListCAmpPhaseSequence_Destroy(listhPlm_HiS);                        \
  if (listhPlm_HiSRDpatch != NULL)                                             \
    SphHarmListCAmpPhaseSequence_Destroy(listhPlm_HiSRDpatch);                 \
  if (listhPlm_AdaS != NULL)                                                   \
    SphHarmListCAmpPhaseSequence_Destroy(listhPlm_AdaS);                       \
  if (*tVecPmodes != NULL)                                                     \
    XLALDestroyREAL8Vector(*tVecPmodes);                                       \
  if (seobdynamicsAdaSHiS != NULL)                                             \
    SEOBdynamics_Destroy(seobdynamicsAdaSHiS);                                 \
  if (listhPlm != NULL)                                                        \
    SphHarmListCAmpPhaseSequence_Destroy(listhPlm);                            \
  if (*hP22_amp != NULL)                                                       \
    XLALDestroyREAL8Vector(*hP22_amp);                                         \
  if (*hP22_phase != NULL)                                                     \
    XLALDestroyREAL8Vector(*hP22_phase);                                       \
  if (*hP21_amp != NULL)                                                       \
    XLALDestroyREAL8Vector(*hP21_amp);                                         \
  if (*hP21_phase != NULL)                                                     \
    XLALDestroyREAL8Vector(*hP21_phase);                                       \
  if (*alphaJ2P != NULL)                                                       \
    XLALDestroyREAL8Vector(*alphaJ2P);                                         \
  if (*betaJ2P != NULL)                                                        \
    XLALDestroyREAL8Vector(*betaJ2P);                                          \
  if (*gammaJ2P != NULL)                                                       \
    XLALDestroyREAL8Vector(*gammaJ2P);                                         \
  if (*hJlm != NULL)                                                           \
    XLALDestroySphHarmTimeSeries(*hJlm);                                       \
  if (*hIlm != NULL)                                                           \
    XLALDestroySphHarmTimeSeries(*hIlm);                                       \
  if (hplusTS != NULL)                                                         \
    XLALDestroyREAL8TimeSeries(hplusTS);                                       \
  if (hcrossTS != NULL)                                                        \
    XLALDestroyREAL8TimeSeries(hcrossTS);                                      \
  if (*mergerParams != NULL)                                                   \
    XLALDestroyREAL8Vector(*mergerParams);                                     \
  if (*seobdynamicsAdaSVector != NULL)                                         \
    XLALDestroyREAL8Vector(*seobdynamicsAdaSVector);                           \
  if (*seobdynamicsHiSVector != NULL)                                          \
    XLALDestroyREAL8Vector(*seobdynamicsHiSVector);                            \
  if (*seobdynamicsAdaSHiSVector != NULL)                                      \
    XLALDestroyREAL8Vector(*seobdynamicsAdaSHiSVector);
#define PRINT_ALL_PARAMS                                                       \
  do {                                                                         \
    XLALPrintError(                                                            \
        "--approximant SEOBNRv4P --f-min %.16e --m1 %.16e --m2 %.16e "         \
        "--spin1x %.16e --spin1y %.16e --spin1z %.16e  --spin2x %.16e "        \
        "--spin2y %.16e --spin2z %.16e --inclination %.16e --distance %.16e "  \
        "--phiRef %.16e --sample-rate %.16e\n",                                \
        fMin, m1SI / LAL_MSUN_SI, m2SI / LAL_MSUN_SI, chi1x, chi1y, chi1z,     \
        chi2x, chi2y, chi2z, inc, r / (1e6 * LAL_PC_SI), phi, 1. / INdeltaT);  \
  } while (0);

/* Compute the highest initial frequency (of the 22 mode):
 * at which the code will generate a waveform. We choose an initial minimum
 * separation of 10.5M as a compromise between reliability of initial conditions
 * and length of the waveform. We use Newtonian Kepler's law. Refuse to
 * generate waveforms shorter than that.
 */
int XLALEOBHighestInitialFreq(
    REAL8 *freqMinRad /**<< OUTPUT, lowest initial 22 mode frequency*/,
    REAL8 mTotal /**<< Total mass in units of solar masses */) {
  REAL8 mTScaled = mTotal * LAL_MTSUN_SI;
  *freqMinRad = pow(10.5, -1.5) / (LAL_PI * mTScaled);
  return XLAL_SUCCESS;
}

/* Implements the standard argmax function for an array of reals */
UNUSED static UINT4 argmax(REAL8Vector *vec) {
  REAL8 max = vec->data[0];
  UINT4 idx_max = 0;
  for (UINT4 i = 0; i < vec->length; i++) {
    if (vec->data[i] > max) {
      max = vec->data[i];
      idx_max = i;
    }
  }
  return idx_max;
}

/* Return a slice of the given vector, not including the higher index,
i.e. parroting Python behaviour */
UNUSED static REAL8Vector *get_slice(REAL8Vector *vec, UINT4 lo, UINT4 hi) {
  UINT4 size = hi - lo;
  REAL8Vector *slice = XLALCreateREAL8Vector(size);
  for (UINT4 jj = 0; jj < size; jj++) {
    slice->data[jj] = vec->data[lo + jj];
  }
  return slice;
}

/* Function to find robustly the peak of a quantity given as an array of
samples. The idea is the scan the samples with a window [w_1,w_2] and find the
local max in each window, where local max has to not lie at the boundaries of
the window. One then keeps track of all the local maxes and picks the largest
one. Finally, one compares this value to the global argmax. If there is a clean
critical point which is also a global max then these 2 values have to agree */
UNUSED static int XLALEOBFindRobustPeak(REAL8 *tPeakQuant, REAL8Vector *tVec,
                                        REAL8Vector *quantVec,
                                        UINT4 window_width) {
  // We begin at the end and go backwards
  UINT4 vlen = tVec->length;
  UINT4 local_argmax = 0;
  UINT4 idx_global = 0;
  UINT4 lo, hi; // Bounds of local window
  // Global argmax over the entire array
  UINT4 global_arg_max = argmax(quantVec);
  UNUSED REAL8 global_max = quantVec->data[global_arg_max];
  REAL8Vector *sl = NULL;
  REAL8 curr_max = 0;
  for (UINT4 kk = vlen - window_width - 1; kk > window_width; kk--) {
    lo = kk - window_width;
    hi = kk + window_width +
         1; // Slice function does not return the value at the end
    sl = get_slice(quantVec, lo, hi);
    local_argmax = argmax(sl);
    if (sl->data[local_argmax] > sl->data[0] &&
        sl->data[local_argmax] > sl->data[sl->length - 1]) {
      // We have *a* local max, figure out it's global index
      // Is the local argmax the largest local argmax so far?
      if (sl->data[local_argmax] > curr_max) {
        curr_max = sl->data[local_argmax];
        idx_global = lo + local_argmax;
      }
    }
    XLALDestroyREAL8Vector(sl);
  }
  *tPeakQuant = 0;
  // Conditions under which we pick the last point of the dynamics:
  // i) we found no local maxes at all
  // ii) the global arg max is larger than the largest of the local maxes
  // by more than 2 % of the largest maxes value (ideally they should be equal)
  // iii) the  peak is  so close to end that we can't interpolate below.

  if (idx_global == 0 ||
      ((quantVec->data[global_arg_max] - quantVec->data[idx_global]) /
           quantVec->data[idx_global] >
       0.1) ||
      (idx_global > tVec->length - 4)) {
    XLAL_PRINT_WARNING("Warning no local max found, using last point\n");
    *tPeakQuant = tVec->data[tVec->length - 1];
    return XLAL_SUCCESS;
  }
  // We have a well-behaved local max. Get the time more accurately.
  // Now we need to interpolate and then set the derivative of the interpolant
  // to 0. We solve this via bisection in an interval of 3 points to the left
  // and right of the argmax.
  gsl_spline *spline = NULL;
  gsl_interp_accel *acc = NULL;
  spline = gsl_spline_alloc(gsl_interp_cspline, quantVec->length);
  acc = gsl_interp_accel_alloc();

  REAL8 time1 = tVec->data[idx_global - 3];
  REAL8 time2 = tVec->data[idx_global + 3];
  REAL8 timePeak = 0;
  REAL8 omegaDerivMid = 0;
  gsl_spline_init(spline, tVec->data, quantVec->data, quantVec->length);
  REAL8 omegaDeriv1 = gsl_spline_eval_deriv(spline, time1, acc);
  do {
    timePeak = (time1 + time2) / 2.;
    omegaDerivMid = gsl_spline_eval_deriv(spline, timePeak, acc);

    if (omegaDerivMid * omegaDeriv1 < 0.0) {
      time2 = timePeak;
    } else {
      omegaDeriv1 = omegaDerivMid;
      time1 = timePeak;
    }
  } while (time2 - time1 > 1.0e-8);
  *tPeakQuant = timePeak;
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  return XLAL_SUCCESS;
}

/* The stopping condition used for the high sampling part of SEOBNRv4P
Will set the termination reason to 1 if terminates normally(i.e. 5 steps
after peak of omega found). Will set it to -1 if something has become nan.
*/
UNUSED static int XLALEOBSpinPrecStopCondition_v4(double UNUSED t,
                                                  const double values[],
                                                  double dvalues[],
                                                  void UNUSED *funcParams) {
  UINT4 counter;
  INT4 i;
  SpinEOBParams UNUSED *params = (SpinEOBParams *)funcParams;

  REAL8 r2 = 0;
  REAL8 p[3], r[3], pdotVec[3], rdotVec[3];
  REAL8 omega, omega_xyz[3];

  memcpy(r, values, 3 * sizeof(REAL8));
  memcpy(p, values + 3, 3 * sizeof(REAL8));
  memcpy(rdotVec, dvalues, 3 * sizeof(REAL8));
  memcpy(pdotVec, dvalues + 3, 3 * sizeof(REAL8));

  r2 = inner_product(r, r);
  cross_product(values, dvalues, omega_xyz);
  omega = sqrt(inner_product(omega_xyz, omega_xyz)) / r2;
  counter = params->eobParams->omegaPeaked;
  if (r2 < 36. && omega < params->eobParams->omega) {
    params->eobParams->omegaPeaked = counter + 1;
  }

  if (params->eobParams->omegaPeaked == 5) {
    return 1;
  }
  for (i = 0; i < 12; i++) {
    if (isnan(dvalues[i]) || isnan(values[i])) {
      params->termination_reason = -1;
      return 1;
    }
  }
  params->eobParams->omega = omega;
  return GSL_SUCCESS;
}

/**
 * Stopping conditions for dynamics integration for SEOBNRv4P
 */
UNUSED static int
XLALEOBSpinPrecStopConditionBasedOnPR(double UNUSED t, const double values[],
                                      double dvalues[],
                                      void UNUSED *funcParams) {
  int debugPK = 0;
  int debugPKverbose = 0;
  INT4 i;
  SpinEOBParams UNUSED *params = (SpinEOBParams *)funcParams;

  REAL8 r2, pDotr = 0;
  REAL8 p[3], r[3], pdotVec[3], rdotVec[3];
  REAL8 omega, omega_xyz[3], L[3], dLdt1[3], dLdt2[3];

  memcpy(r, values, 3 * sizeof(REAL8));
  memcpy(p, values + 3, 3 * sizeof(REAL8));
  memcpy(rdotVec, dvalues, 3 * sizeof(REAL8));
  memcpy(pdotVec, dvalues + 3, 3 * sizeof(REAL8));

  r2 = inner_product(r, r);
  cross_product(values, dvalues, omega_xyz);
  omega = sqrt(inner_product(omega_xyz, omega_xyz)) / r2;
  pDotr = inner_product(p, r) / sqrt(r2);
  if (debugPK) {
    XLAL_PRINT_INFO("XLALEOBSpinPrecStopConditionBasedOnPR:: r = %e %e\n",
                    sqrt(r2), omega);
  }
  if (debugPK) {
    XLAL_PRINT_INFO(
        "XLALEOBSpinPrecStopConditionBasedOnPR:: values = %e %e %e %e %e %e\n",
        values[6], values[7], values[8], values[9], values[10], values[11]);
  }
  if (debugPK) {
    XLAL_PRINT_INFO(
        "XLALEOBSpinPrecStopConditionBasedOnPR:: dvalues = %e %e %e %e %e %e\n",
        dvalues[6], dvalues[7], dvalues[8], dvalues[9], dvalues[10],
        dvalues[11]);
  }
  REAL8 rdot;
  // this is d(r)/dt obtained by differentiating r2 (see above)
  rdot = inner_product(rdotVec, r) / sqrt(r2);
  // This is d/dt(pDotr) see pDotr above.
  double prDot = -inner_product(p, r) * rdot / r2 +
                 inner_product(pdotVec, r) / sqrt(r2) +
                 inner_product(rdotVec, p) / sqrt(r2);

  cross_product(r, pdotVec, dLdt1);
  cross_product(rdotVec, p, dLdt2);
  cross_product(r, p, L);

  /* ********************************************************** */
  /* *******  Different termination conditions Follow  ******** */
  /* ********************************************************** */

  /* Table of termination conditions

    Value                   Reason
    -1                   Any of the derivatives are Nan
    0                    r < 8 and pDotr >= 0 (outspiraling)
    1                    r < 8 and rdot >= 0 (outspiraling)
    2                    r < 2 and prDot > 0 (dp_r/dt is growing)
    3                    r < 8 and |p_vec| > 10 (the momentum vector is large)
    4                    r < 8 and |p_vec| < 1e-10 (momentum vector is small)
    5                    r < 2 and omega has a another peak
    6                    r < 8 and omega < 0.04 or (r < 2. and  omega < 0.14 and
    omega has a peak) 7                    r < 8 and omega > 1 (unphysical
    omega) 8                    r < 5 and any of  |dp_i/dt| > 10 9 r < 8 and
    pphi > 10
    10                   r < 3 and rdot increases
  */

  /* Terminate if any derivative is Nan */
  for (i = 0; i < 12; i++) {
    if (isnan(dvalues[i]) || isnan(values[i])) {
      if (debugPK) {
        XLAL_PRINT_INFO("\n  isnan reached. r2 = %f\n", r2);
        fflush(NULL);
      }
      XLALPrintError("XLAL Error - %s: nan reached at r2 = %f \n", __func__,
                     r2);
      XLAL_ERROR(XLAL_EINVAL);
      params->termination_reason = -1;
      return 1;
    }
  }

  /* ********************************************************** */
  /* *******  Unphysical orbital conditions  ******** */
  /* ********************************************************** */

  /* Terminate if p_r points outwards */
  if (r2 < 16 && pDotr >= 0) {
    if (debugPK) {
      XLAL_PRINT_INFO(
          "\n Integration stopping, p_r pointing outwards -- out-spiraling!\n");
      fflush(NULL);
    }
    params->termination_reason = 0;

    return 1;
  }

  /* Terminate if rdot is >0 (OUTspiraling) for separation <4M */
  if (r2 < 16 && rdot >= 0) {
    if (debugPK) {
      XLAL_PRINT_INFO("\n Integration stopping, dr/dt>0 -- out-spiraling!\n");
      fflush(NULL);
    }
    params->termination_reason = 1;

    return 1;
  }

  /* Terminate if dp_R/dt > 0, i.e. radial momentum is increasing for separation
   * <2M */
  if (r2 < 4. && prDot > 0.) {
    if (debugPK) {
      XLAL_PRINT_INFO("\n Integration stopping as prDot = %lf at r = %lf\n",
                      prDot, sqrt(r2));
      fflush(NULL);
    }
    params->termination_reason = 2;

    return 1;
  }

  if (r2 < 16. && (sqrt(values[3] * values[3] + values[4] * values[4] +
                        values[5] * values[5]) > 10.)) {
    if (debugPK)
      XLAL_PRINT_INFO("\n Integration stopping |pvec|> 10\n");
    fflush(NULL);
    params->termination_reason = 3;

    return 1;
  }

  if (r2 < 16. && (sqrt(values[3] * values[3] + values[4] * values[4] +
                        values[5] * values[5]) < 1.e-10)) {
    if (debugPK)
      XLAL_PRINT_INFO("\n Integration stopping |pvec|<1e-10\n");
    fflush(NULL);
    params->termination_reason = 4;

    return 1;
  }

  /* **************************************************************** */
  /*                         Omega related                            */
  /* **************************************************************** */
  /* Terminate when omega reaches peak, and separation is < 4M */
  if (r2 < 16. && omega < params->eobParams->omega)
    params->eobParams->omegaPeaked = 1;

  /* If omega has gone through a second extremum, break */
  if (r2 < 4. && params->eobParams->omegaPeaked == 1 &&
      omega > params->eobParams->omega) {
    if (debugPK) {
      XLAL_PRINT_INFO(
          "\n Integration stopping, omega reached second extremum\n");
      fflush(NULL);
    }
    params->termination_reason = 5;

    return 1;
  }

  /* If Momega did not evolve above 0.01 even though r < 4 or omega<0.14 for
   * r<2, break */
  if ((r2 < 16. && omega < 0.04) ||
      (r2 < 4. && omega < 0.14 && params->eobParams->omegaPeaked == 1)) {
    if (debugPK) {
      XLAL_PRINT_INFO("\n Integration stopping for omega below threshold, "
                      "omega=%f at r = %f\n",
                      omega, sqrt(r2));
      fflush(NULL);
    }
    params->termination_reason = 6;

    return 1;
  }

  if (r2 < 16. && omega > 1.) {
    if (debugPK) {
      XLAL_PRINT_INFO("\n Integration stopping, omega>1 at r = %f\n", sqrt(r2));
      fflush(NULL);
    }
    params->termination_reason = 7;

    return 1;
  }
  params->eobParams->omega = omega;

  /* **************************************************************** */
  /*              related to Numerical values of x/p/derivatives      */
  /* **************************************************************** */

  /* If momentum derivatives are too large numerically, break */
  if (r2 < 25 && (fabs(dvalues[3]) > 10 || fabs(dvalues[4]) > 10 ||
                  fabs(dvalues[5]) > 10)) {
    if (debugPK) {
      XLAL_PRINT_INFO("\n Integration stopping, dpdt > 10 -- too large!\n");
      fflush(NULL);
    }
    params->termination_reason = 8;

    return 1;
  }

  /* If p_\Phi is too large numerically, break */
  if (r2 < 16. && values[5] > 10) {
    if (debugPK) {
      XLAL_PRINT_INFO("Integration stopping, Pphi > 10 now\n\n");
      fflush(NULL);
    }
    params->termination_reason = 9;

    return 1;
  }
  /* If rdot inclreases, break */
  if (r2 < 9. && rdot > params->prev_dr) {
    if (debugPK) {
      XLAL_PRINT_INFO("\n Integration stopping, dr/dt increasing!\n");
      fflush(NULL);
    }
    params->prev_dr = rdot;
    params->termination_reason = 10;

    return 1;
  }
  params->prev_dr = rdot;

  /* **************************************************************** */
  /*              Last resort conditions                              */
  /* **************************************************************** */

  /* Very verbose output */
  if (debugPKverbose && r2 < 16.) {
    XLAL_PRINT_INFO("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", t,
                    values[0], values[1], values[2], values[3], values[4],
                    values[5], values[6], values[7], values[8], values[9],
                    values[10], values[11], values[12], values[13], omega);
  }

  return GSL_SUCCESS;
}

/**
 * Stopping condition for the regular resolution SEOBNRv1/2 orbital evolution
 * -- stop when reaching max orbital frequency in strong field.
 * At each test,
 * if omega starts to decrease, return 1 to stop evolution;
 * if not, update omega with current value and return GSL_SUCCESS to continue
 * evolution.
 */
static int XLALEOBSpinPrecAlignedStopCondition(
    double UNUSED t,       /**< UNUSED */
    const double values[], /**< dynamical variable values */
    double dvalues[],      /**< dynamical variable time derivative values */
    void *funcParams       /**< physical parameters */
) {
  int debugPK = 0;
  REAL8 omega, r;
  SpinEOBParams *params = (SpinEOBParams *)funcParams;

  r = values[0];
  omega = dvalues[1];
  if (debugPK) {
    XLAL_PRINT_INFO("XLALEOBSpinPrecAlignedStopCondition:: r = %e\n", r);
  }

  if (r < 6. && omega < params->eobParams->omega) {
    return 1;
  }

  params->eobParams->omega = omega;
  return GSL_SUCCESS;
}

/**
 * Stopping condition for the high resolution SEOBNRv4.
 */
static int XLALSpinPrecAlignedHiSRStopCondition(
    double UNUSED t,              /**< UNUSED */
    const double UNUSED values[], /**< dynamical variable values */
    double dvalues[],       /**< dynamical variable time derivative values */
    void UNUSED *funcParams /**< physical parameters */
) {

  REAL8 omega, r;
  UINT4 counter;
  SpinEOBParams *params = (SpinEOBParams *)funcParams;
  r = values[0];
  omega = dvalues[1];
  counter = params->eobParams->omegaPeaked;

  if (r < 6. && omega < params->eobParams->omega) {

    params->eobParams->omegaPeaked = counter + 1;
  }
  if (dvalues[2] >= 0. || params->eobParams->omegaPeaked == 5 ||
      isnan(dvalues[3]) || isnan(dvalues[2]) || isnan(dvalues[1]) ||
      isnan(dvalues[0])) {

    return 1;
  }
  params->eobParams->omega = omega;
  return GSL_SUCCESS;
}

/**
 * ModeArray is a structure which allows to select the the co-precessing frame
 * modes to include in the waveform. This function will create a structure with
 * the default modes for every model
 */
static INT4 XLALSetup_EOB__std_mode_array_structure(LALValue *ModeArray,
                                                    UINT4 PrecEOBversion) {

  /* setup ModeArray */

  if (PrecEOBversion == v4PHM.number) {
    /* Adding all the modes of SEOBNRv4PHM
    * i.e. [(2,2),(2,1),(3,3),(4,4),(5,5)]
    the relative -m modes are added automatically*/
    XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 2);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 1);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 3, 3);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 4, 4);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 5, 5);
  }
  if (PrecEOBversion == v4P.number) {
    /* Adding all the modes of SEOBNRv4P
    * i.e. [(2,2),(2,1)]
    the relative -m modes are added automatically*/
    XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 2);
    XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 1);
  }

  return XLAL_SUCCESS;
}

/**
 * ModeArray is a structure which allows to select the the co-precessing frame
 * modes to include in the waveform. This function check if the selected modes
 * are available for a given model
 */
static INT4 XLALCheck_EOB_mode_array_structure(LALValue *ModeArray,
                                               UINT4 PrecEOBversion) {
  INT4 flagTrue = 0;
  UINT4 modeL;
  UINT4 modeM;
  UINT4 nModes;
  const UINT4 lmModes[5][2] = {{2, 2}, {2, 1}, {3, 3}, {4, 4}, {5, 5}};
  if (PrecEOBversion == v4PHM.number) {
    /*If one select SEOBNRv4PHM all the modes above are selected to check
     */
    nModes = 5;
  } else {
    /*If not only the modes 22 and 21 are selected to check
     */
    nModes = 2;
  }
  /* First check if the user is entering a mode with negative m */
  /* This function only takes positive m and then selects automatically +- m */
  for (UINT4 ELL = 2; ELL <= _SEOB_MODES_LMAX; ELL++) {
    for (INT4 EMM = -ELL; EMM < 0; EMM++) {
      if (XLALSimInspiralModeArrayIsModeActive(ModeArray, ELL, EMM) == 1) {
        XLALPrintError("Mode (%d,%d) has a negative m. \
        In mode array you should specify (l,|m|). The code will automatically return +- m modes\n",
                       ELL, EMM);
        return XLAL_FAILURE;
      }
    }
  }

  /*Loop over all the possible modes
   *we only check +m modes, when one select (l,m) mode is actually
   *selecting (l,|m|) mode
   */
  for (UINT4 ELL = 2; ELL <= _SEOB_MODES_LMAX; ELL++) {
    for (UINT4 EMM = 0; EMM <= ELL; EMM++) {
      if (XLALSimInspiralModeArrayIsModeActive(ModeArray, ELL, EMM) == 1) {
        for (UINT4 k = 0; k < nModes; k++) {
          modeL = lmModes[k][0];
          modeM = lmModes[k][1];
          if ((modeL == ELL) && (modeM == EMM)) {
            flagTrue = 1;
          }
        }
        /*For each active mode check if is available for the selected model
         */
        if (flagTrue != 1) {
          if (PrecEOBversion == v4PHM.number) {
            XLALPrintError("Mode (%d,%d) is not available for the model %s.\n",
                           ELL, EMM, v4PHM.name);
            return XLAL_FAILURE;
          }
          if (PrecEOBversion == v4P.number) {
            XLALPrintError("Mode (%d,%d) is not available for the model %s.\n",
                           ELL, EMM, v4P.name);
            return XLAL_FAILURE;
          }
        }
        flagTrue = 0;
      }
    }
  }

  return XLAL_SUCCESS;
}

/**
 * Standard interface for SEOBNRv4P waveform generator: calls
 * XLALSimIMRSpinPrecEOBWaveformAll
 */
int XLALSimIMRSpinPrecEOBWaveform(
    REAL8TimeSeries **hplus,  /**<< OUTPUT, +-polarization waveform */
    REAL8TimeSeries **hcross, /**<< OUTPUT, x-polarization waveform */
    const REAL8 phiC,         /**<< coalescence orbital phase (rad) */
    const REAL8 deltaT,       /**<< sampling time step */
    const REAL8 m1SI,         /**<< mass-1 in SI unit (kg) */
    const REAL8 m2SI,         /**<< mass-2 in SI unit (kg) 8*/
    const REAL8 fMin,         /**<< starting frequency (Hz) */
    const REAL8 r,            /**<< luminosity distance in SI unit (m) */
    const REAL8 inc,          /**<< inclination angle */
    const REAL8 INspin1[],    /**<< spin1 */
    const REAL8 INspin2[],    /**<< spin2 */
    UNUSED const UINT4
        PrecEOBversion, /**<< Precessing EOB waveform generator model */
    LALDict *LALParams  /**<< Dictionary of additional wf parameters */
) {

  REAL8Vector *dyn_Low = NULL;
  REAL8Vector *dyn_Hi = NULL;
  REAL8Vector *dyn_all = NULL;
  REAL8Vector *t_vec_modes = NULL;
  REAL8Vector *hP22_amp = NULL;
  REAL8Vector *hP22_phase = NULL;
  REAL8Vector *hP21_amp = NULL;
  REAL8Vector *hP21_phase = NULL;
  REAL8Vector *hP33_amp = NULL;
  REAL8Vector *hP33_phase = NULL;
  REAL8Vector *hP44_amp = NULL;
  REAL8Vector *hP44_phase = NULL;
  REAL8Vector *hP55_amp = NULL;
  REAL8Vector *hP55_phase = NULL;
  REAL8Vector *alphaJ2P = NULL;
  REAL8Vector *betaJ2P = NULL;
  REAL8Vector *gammaJ2P = NULL;
  REAL8Vector *AttachPars = NULL;

  /** This time series contains harmonics in precessing (P) frame, no RD, for
   * the end of the signal (high samling part)*/
  SphHarmTimeSeries *hIlm = NULL;
  /** This stores harmonics in J-frame, no RD, for the end of the signal (high
   * sampling part) */
  SphHarmTimeSeries *hJlm = NULL;

  /* Import the set of modes requested by the user if available, if not
  load the default modes  */
  LALValue *modearray = XLALSimInspiralWaveformParamsLookupModeArray(LALParams);
  if (modearray == NULL) {
    modearray = XLALSimInspiralCreateModeArray();
    XLALSetup_EOB__std_mode_array_structure(modearray, PrecEOBversion);
  }
  /*Check that the modes chosen are available for the given model*/
  if (XLALCheck_EOB_mode_array_structure(modearray, PrecEOBversion) ==
      XLAL_FAILURE) {
    XLALPrintError("Not available mode chosen.\n");
    XLAL_ERROR(XLAL_EFUNC);
  }

  /* Set of SEOB flags */
  LALDict *seobflags = XLALCreateDict();
  /* Spin-aligned model v4 */
  XLALDictInsertINT4Value(seobflags, "SEOBNRv4P_SpinAlignedEOBversion", 4);
  /* Generate P-frame modes m<0 with the symmetry hP_l-m ~ (-1)^l hP_lm* */
  XLALDictInsertINT4Value(seobflags, "SEOBNRv4P_SymmetrizehPlminusm", 1);
  /* Use numerical or analytical derivatives of the Hamiltonian
   Default is numerical with the flag 1*/
  INT4 NumericalOrAnalyticalHamiltonianDerivative =
      XLALSimInspiralWaveformParamsLookupEOBChooseNumOrAnalHamDer(LALParams);
  /* NumericalOrAnalyticalHamiltonianDerivative can only be 0 (analytical) or 1
   * (numerical), let's check! */
  if ((NumericalOrAnalyticalHamiltonianDerivative != 0) &&
      (NumericalOrAnalyticalHamiltonianDerivative != 1)) {
    XLALPrintError("XLAL Error - %s: Unknown value for the derivative of the "
                   "Hamiltonian flag. \nAt present only "
                   "1 (numerical derivative) or 0 (analytical derivative) are "
                   "available.\n",
                   __func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (NumericalOrAnalyticalHamiltonianDerivative ==
      FLAG_SEOBNRv4P_HAMILTONIAN_DERIVATIVE_NUMERICAL) {
    XLALDictInsertINT4Value(seobflags, "SEOBNRv4P_HamiltonianDerivative",
                            FLAG_SEOBNRv4P_HAMILTONIAN_DERIVATIVE_NUMERICAL);
  } else {
    XLALDictInsertINT4Value(seobflags, "SEOBNRv4P_HamiltonianDerivative",
                            NumericalOrAnalyticalHamiltonianDerivative);
  }
  /* Extension of Euler angles post-merger: simple precession around final J at
   * a rate set by QNMs */
  XLALDictInsertINT4Value(seobflags, "SEOBNRv4P_euler_extension",
                          FLAG_SEOBNRv4P_EULEREXT_QNM_SIMPLE_PRECESSION);
  /* Z-axis of the radiation frame L */
  XLALDictInsertINT4Value(seobflags, "SEOBNRv4P_Zframe",
                          FLAG_SEOBNRv4P_ZFRAME_L);
  /* No debug output */
  XLALDictInsertINT4Value(seobflags, "SEOBNRv4P_debug", 0);
  // What max ell to use for the Nyquist check.
  // Note that this has no effect at which waveform modes are actually being generated.
  // For SEOBNRv4P we do not give the user a choice, the max ell is *always* 2
  INT4 ellMaxForNyquistCheck = 2;
  if (PrecEOBversion == v4PHM.number){
    ellMaxForNyquistCheck =
        XLALSimInspiralWaveformParamsLookupEOBEllMaxForNyquistCheck(LALParams);
  }
  XLALDictInsertINT4Value(seobflags,"ellMaxForNyquistCheck",ellMaxForNyquistCheck);
  int ret = XLAL_SUCCESS;
  XLAL_TRY(XLALSimIMRSpinPrecEOBWaveformAll(
               hplus, hcross, &hIlm, &hJlm, &dyn_Low, &dyn_Hi, &dyn_all,
               &t_vec_modes, &hP22_amp, &hP22_phase, &hP21_amp, &hP21_phase,
               &hP33_amp, &hP33_phase, &hP44_amp, &hP44_phase, &hP55_amp,
               &hP55_phase, &alphaJ2P, &betaJ2P, &gammaJ2P, &AttachPars, phiC,
               deltaT, m1SI, m2SI, fMin, r, inc, INspin1[0], INspin1[1],
               INspin1[2], INspin2[0], INspin2[1], INspin2[2], modearray,
               seobflags),
           ret);
  /*
  if (ret == XLAL_SUCCESS) {
    if (*hplus == NULL || *hcross == NULL) {
      XLALPrintError(
          "Houston-2, we've got a problem SOS, SOS, SOS, the waveform "
          "generator returns NULL!!!... m1 = %.18e, m2 = %.18e, fMin = %.18e, "
          "inclination = %.18e,   spin1 = {%.18e, %.18e, %.18e},   spin2 = "
          "{%.18e, %.18e, %.18e} \n",
          m1SI / LAL_MSUN_SI, m2SI / LAL_MSUN_SI, (double)fMin, (double)inc,
          INspin1[0], INspin1[1], INspin1[2], INspin2[0], INspin2[1],
          INspin2[2]);
      XLAL_ERROR(XLAL_ENOMEM);
    }
    if ((*hplus)->data == NULL || (*hcross)->data == NULL) {
      XLALPrintError(
          "Houston-3, we've got a problem SOS, SOS, SOS, the waveform "
          "generator returns NULL!!!... m1 = %.18e, m2 = %.18e, fMin = %.18e, "
          "inclination = %.18e,   spin1 = {%.18e, %.18e, %.18e},   spin2 = "
          "{%.18e, %.18e, %.18e} \n",
          m1SI / LAL_MSUN_SI, m2SI / LAL_MSUN_SI, (double)fMin, (double)inc,
          INspin1[0], INspin1[1], INspin1[2], INspin2[0], INspin2[1],
          INspin2[2]);
      XLAL_ERROR(XLAL_ENOMEM);
    }
    }*/
  if (modearray)
    XLALDestroyValue(modearray);
  if (seobflags)
    XLALDestroyDict(seobflags);
  if (dyn_Low)
    XLALDestroyREAL8Vector(dyn_Low);
  if (dyn_Hi)
    XLALDestroyREAL8Vector(dyn_Hi);
  if (dyn_all)
    XLALDestroyREAL8Vector(dyn_all);

  if (t_vec_modes)
    XLALDestroyREAL8Vector(t_vec_modes);
  if (hP22_amp)
    XLALDestroyREAL8Vector(hP22_amp);
  if (hP22_phase)
    XLALDestroyREAL8Vector(hP22_phase);
  if (hP21_amp)
    XLALDestroyREAL8Vector(hP21_amp);
  if (hP21_phase)
    XLALDestroyREAL8Vector(hP21_phase);
  if (hP33_amp)
    XLALDestroyREAL8Vector(hP33_amp);
  if (hP33_phase)
    XLALDestroyREAL8Vector(hP33_phase);
  if (hP44_amp)
    XLALDestroyREAL8Vector(hP44_amp);
  if (hP44_phase)
    XLALDestroyREAL8Vector(hP44_phase);
  if (hP55_amp)
    XLALDestroyREAL8Vector(hP55_amp);
  if (hP55_phase)
    XLALDestroyREAL8Vector(hP55_phase);

  if (alphaJ2P)
    XLALDestroyREAL8Vector(alphaJ2P);
  if (betaJ2P)
    XLALDestroyREAL8Vector(betaJ2P);
  if (gammaJ2P)
    XLALDestroyREAL8Vector(gammaJ2P);
  if (AttachPars)
    XLALDestroyREAL8Vector(AttachPars);
  if (hIlm)
    XLALDestroySphHarmTimeSeries(hIlm);
  if (hJlm)
    XLALDestroySphHarmTimeSeries(hJlm);
  if (ret != XLAL_SUCCESS) {
    XLAL_ERROR(ret);
  }
  return ret;
}
/**
 * This function returns the maximum ell in the mode array.
 * Note that m<=0 modes are ignored and a warning given.
 */
static int
SEOBGetLMaxInModeArray(LALValue *modearray, /**<< Input: ModeArray structure */
                       int lmax /**<< Input: maximum value of l to explore --
                                   possible modes with l>lmax will be ignored */
) {
  /* Populate array */
  INT4 lmax_array = 0;
  for (INT4 l = 2; l <= lmax; l++) {
    for (INT4 m = l; m >= -l; m--) {
      if (m > 0) {
        if (XLALSimInspiralModeArrayIsModeActive(modearray, l, m)) {
          if (lmax_array < l)
            lmax_array = l;
        }
      } else {
        XLAL_PRINT_WARNING(
            "XLAL Warning - %s: mode (l,m)=(%d,%d) present in mode array -- "
            "in our conventions, we only consider m>0. Mode ignored for "
            "counting.\n",
            __func__, l, m);
      }
    }
  }

  return lmax_array;
}

/**
 * Standard interface for SEOBNRv4P modes generator: calls
 * XLALSimIMRSpinPrecEOBWaveformAll
 */
SphHarmTimeSeries *XLALSimIMRSpinPrecEOBModes(
    const REAL8 deltaT,    /**<< sampling time step */
    const REAL8 m1SI,      /**<< mass-1 in SI unit (kg) */
    const REAL8 m2SI,      /**<< mass-2 in SI unit (kg) 8*/
    const REAL8 fMin,      /**<< starting frequency (Hz) */
    const REAL8 r,         /**<< luminosity distance in SI unit (m) */
    const REAL8 INspin1[], /**<< spin1 */
    const REAL8 INspin2[], /**<< spin2 */
    UNUSED const UINT4
        PrecEOBversion, /**<< Precessing EOB waveform generator model */
    LALDict *LALParams  /**<< Dictionary of additional wf parameters */
) {

  REAL8Vector *dyn_Low = NULL;
  REAL8Vector *dyn_Hi = NULL;
  REAL8Vector *dyn_all = NULL;
  REAL8Vector *t_vec_modes = NULL;
  REAL8Vector *hP22_amp = NULL;
  REAL8Vector *hP22_phase = NULL;
  REAL8Vector *hP21_amp = NULL;
  REAL8Vector *hP21_phase = NULL;
  REAL8Vector *hP33_amp = NULL;
  REAL8Vector *hP33_phase = NULL;
  REAL8Vector *hP44_amp = NULL;
  REAL8Vector *hP44_phase = NULL;
  REAL8Vector *hP55_amp = NULL;
  REAL8Vector *hP55_phase = NULL;
  REAL8Vector *alphaJ2P = NULL;
  REAL8Vector *betaJ2P = NULL;
  REAL8Vector *gammaJ2P = NULL;
  REAL8Vector *AttachPars = NULL;
  REAL8TimeSeries *hplus = NULL;
  REAL8TimeSeries *hcross = NULL;

  /** This time series contains harmonics in precessing (P) frame, no RD, for
   * the end of the signal (high samling part)*/
  SphHarmTimeSeries *hIlm = NULL;
  /** This stores harmonics in J-frame, no RD, for the end of the signal (high
   * sampling part) */
  SphHarmTimeSeries *hJlm = NULL;

  /* Import the set of modes requested by the user if available, if not
  load the default modes  */
  LALValue *modearray = XLALSimInspiralWaveformParamsLookupModeArray(LALParams);
  if (modearray == NULL) {
    modearray = XLALSimInspiralCreateModeArray();
    XLALSetup_EOB__std_mode_array_structure(modearray, PrecEOBversion);
  }
  /*Check that the modes chosen are available for the given model*/
  if (XLALCheck_EOB_mode_array_structure(modearray, PrecEOBversion) ==
      XLAL_FAILURE) {
    XLALPrintError("Not available mode chosen.\n");
    XLAL_ERROR_NULL(XLAL_EFUNC);
  }

  /* Set of SEOB flags */
  LALDict *seobflags = XLALCreateDict();
  /* Spin-aligned model v4 */
  XLALDictInsertINT4Value(seobflags, "SEOBNRv4P_SpinAlignedEOBversion", 4);
  /* Generate P-frame modes m<0 with the symmetry hP_l-m ~ (-1)^l hP_lm* */
  XLALDictInsertINT4Value(seobflags, "SEOBNRv4P_SymmetrizehPlminusm", 1);
  /* Use numerical or analytical derivatives of the Hamiltonian
   Default is numerical with the flag 1*/
  INT4 NumericalOrAnalyticalHamiltonianDerivative =
      XLALSimInspiralWaveformParamsLookupEOBChooseNumOrAnalHamDer(LALParams);
  /* NumericalOrAnalyticalHamiltonianDerivative can only be 0 (analytical) or 1
   * (numerical), let's check! */
  if ((NumericalOrAnalyticalHamiltonianDerivative != 0) &&
      (NumericalOrAnalyticalHamiltonianDerivative != 1)) {
    XLALPrintError("XLAL Error - %s: Unknown value for the derivative of the "
                   "Hamiltonian flag. \nAt present only "
                   "1 (numerical derivative) or 0 (analytical derivative) are "
                   "available.\n",
                   __func__);
    XLAL_ERROR_NULL(XLAL_EFUNC);
  }
  if (NumericalOrAnalyticalHamiltonianDerivative ==
      FLAG_SEOBNRv4P_HAMILTONIAN_DERIVATIVE_NUMERICAL) {
    XLALDictInsertINT4Value(seobflags, "SEOBNRv4P_HamiltonianDerivative",
                            FLAG_SEOBNRv4P_HAMILTONIAN_DERIVATIVE_NUMERICAL);
  } else {
    XLALDictInsertINT4Value(seobflags, "SEOBNRv4P_HamiltonianDerivative",
                            NumericalOrAnalyticalHamiltonianDerivative);
  }
  /* Extension of Euler angles post-merger: simple precession around final J at
   * a rate set by QNMs */
  XLALDictInsertINT4Value(seobflags, "SEOBNRv4P_euler_extension",
                          FLAG_SEOBNRv4P_EULEREXT_QNM_SIMPLE_PRECESSION);
  /* Z-axis of the radiation frame L */
  XLALDictInsertINT4Value(seobflags, "SEOBNRv4P_Zframe",
                          FLAG_SEOBNRv4P_ZFRAME_L);
  /* No debug output */
  XLALDictInsertINT4Value(seobflags, "SEOBNRv4P_debug", 0);

  INT4 ellMaxForNyquistCheck = 2;
  if (PrecEOBversion == v4PHM.number){
    ellMaxForNyquistCheck =
        XLALSimInspiralWaveformParamsLookupEOBEllMaxForNyquistCheck(LALParams);
  }
  XLALDictInsertINT4Value(seobflags,"ellMaxForNyquistCheck",ellMaxForNyquistCheck);
  int ret = XLAL_SUCCESS;
  XLAL_TRY(XLALSimIMRSpinPrecEOBWaveformAll(
               &hplus, &hcross, &hIlm, &hJlm, &dyn_Low, &dyn_Hi, &dyn_all,
               &t_vec_modes, &hP22_amp, &hP22_phase, &hP21_amp, &hP21_phase,
               &hP33_amp, &hP33_phase, &hP44_amp, &hP44_phase, &hP55_amp,
               &hP55_phase, &alphaJ2P, &betaJ2P, &gammaJ2P, &AttachPars, 0.,
               deltaT, m1SI, m2SI, fMin, r, 0., INspin1[0], INspin1[1],
               INspin1[2], INspin2[0], INspin2[1], INspin2[2], modearray,
               seobflags),
           ret);
  if (ret != XLAL_SUCCESS)
  {
    // We have to clean up right here
    if (modearray)
      XLALDestroyValue(modearray);
    if (seobflags)
      XLALDestroyDict(seobflags);
    if (dyn_Low)
      XLALDestroyREAL8Vector(dyn_Low);
    if (dyn_Hi)
      XLALDestroyREAL8Vector(dyn_Hi);
    if (dyn_all)
      XLALDestroyREAL8Vector(dyn_all);

    if (t_vec_modes)
      XLALDestroyREAL8Vector(t_vec_modes);
    if (hP22_amp)
      XLALDestroyREAL8Vector(hP22_amp);
    if (hP22_phase)
      XLALDestroyREAL8Vector(hP22_phase);
    if (hP21_amp)
      XLALDestroyREAL8Vector(hP21_amp);
    if (hP21_phase)
      XLALDestroyREAL8Vector(hP21_phase);
    if (hP33_amp)
      XLALDestroyREAL8Vector(hP33_amp);
    if (hP33_phase)
      XLALDestroyREAL8Vector(hP33_phase);
    if (hP44_amp)
      XLALDestroyREAL8Vector(hP44_amp);
    if (hP44_phase)
      XLALDestroyREAL8Vector(hP44_phase);
    if (hP55_amp)
      XLALDestroyREAL8Vector(hP55_amp);
    if (hP55_phase)
      XLALDestroyREAL8Vector(hP55_phase);

    if (alphaJ2P)
      XLALDestroyREAL8Vector(alphaJ2P);
    if (betaJ2P)
      XLALDestroyREAL8Vector(betaJ2P);
    if (gammaJ2P)
      XLALDestroyREAL8Vector(gammaJ2P);
    if (AttachPars)
      XLALDestroyREAL8Vector(AttachPars);
    if (hJlm)
      XLALDestroySphHarmTimeSeries(hJlm);
    // Fail correctly
    XLAL_ERROR_NULL(ret);
  }
  /* Here we multiply the appropriate factor to the modes to convert them in
   * dimensional units */
  REAL8 m1 = m1SI / LAL_MSUN_SI;
  REAL8 m2 = m2SI / LAL_MSUN_SI;
  REAL8 mTotal = m1 + m2;
  REAL8 mTScaled = mTotal * LAL_MTSUN_SI;

  /* Initialize amplitude factor */
  REAL8 amp0 = mTotal * LAL_MRSUN_SI / r;
  char mode_string[32];

  INT4 modes_lmax = SEOBGetLMaxInModeArray(modearray, _SEOB_MODES_LMAX);
  UINT4 retLen = hIlm->tdata->length;
  SphHarmTimeSeries *hIlm_dimfull = NULL;
  LIGOTimeGPS tGPS = LIGOTIMEGPSZERO;
  XLALGPSAdd(&tGPS, -mTScaled * AttachPars->data[2]);

  for (INT4 l = 2; l <= modes_lmax; l++) {
    for (INT4 m = -l; m <= l; m++) {

      /* Get dimensionless mode hIlm */
      COMPLEX16TimeSeries *mode_dimless =
          XLALSphHarmTimeSeriesGetMode(hIlm, l, m);

      /* Time series for dimensionful mode */
      COMPLEX16TimeSeries *mode_dimfull = XLALCreateCOMPLEX16TimeSeries(
          mode_string, &tGPS, 0., deltaT, &lalStrainUnit, retLen);
      memset(mode_dimfull->data->data, 0, retLen * sizeof(COMPLEX16));

      for (UINT4 i = 0; i < mode_dimless->data->length; i++) {
        /* The factor -1 here is to rotate EOB modes in LAL convention
         * see https://dcc.ligo.org/LIGO-G1900275  */
        mode_dimfull->data->data[i] =
            -1. * (COMPLEX16)amp0 * mode_dimless->data->data[i];
      }
      hIlm_dimfull =
          XLALSphHarmTimeSeriesAddMode(hIlm_dimfull, mode_dimfull, l, m);
      XLALDestroyCOMPLEX16TimeSeries(mode_dimless);
      XLALDestroyCOMPLEX16TimeSeries(mode_dimfull);
    }
  }

  /*
  if (ret == XLAL_SUCCESS) {
    if (*hplus == NULL || *hcross == NULL) {
      XLALPrintError(
          "Houston-2, we've got a problem SOS, SOS, SOS, the waveform "
          "generator returns NULL!!!... m1 = %.18e, m2 = %.18e, fMin = %.18e, "
          "inclination = %.18e,   spin1 = {%.18e, %.18e, %.18e},   spin2 = "
          "{%.18e, %.18e, %.18e} \n",
          m1SI / LAL_MSUN_SI, m2SI / LAL_MSUN_SI, (double)fMin, (double)inc,
          INspin1[0], INspin1[1], INspin1[2], INspin2[0], INspin2[1],
          INspin2[2]);
      XLAL_ERROR(XLAL_ENOMEM);
    }
    if ((*hplus)->data == NULL || (*hcross)->data == NULL) {
      XLALPrintError(
          "Houston-3, we've got a problem SOS, SOS, SOS, the waveform "
          "generator returns NULL!!!... m1 = %.18e, m2 = %.18e, fMin = %.18e, "
          "inclination = %.18e,   spin1 = {%.18e, %.18e, %.18e},   spin2 = "
          "{%.18e, %.18e, %.18e} \n",
          m1SI / LAL_MSUN_SI, m2SI / LAL_MSUN_SI, (double)fMin, (double)inc,
          INspin1[0], INspin1[1], INspin1[2], INspin2[0], INspin2[1],
          INspin2[2]);
      XLAL_ERROR(XLAL_ENOMEM);
    }
    }*/

  if (modearray)
    XLALDestroyValue(modearray);
  if (seobflags)
    XLALDestroyDict(seobflags);
  if (dyn_Low)
    XLALDestroyREAL8Vector(dyn_Low);
  if (dyn_Hi)
    XLALDestroyREAL8Vector(dyn_Hi);
  if (dyn_all)
    XLALDestroyREAL8Vector(dyn_all);

  if (t_vec_modes)
    XLALDestroyREAL8Vector(t_vec_modes);
  if (hP22_amp)
    XLALDestroyREAL8Vector(hP22_amp);
  if (hP22_phase)
    XLALDestroyREAL8Vector(hP22_phase);
  if (hP21_amp)
    XLALDestroyREAL8Vector(hP21_amp);
  if (hP21_phase)
    XLALDestroyREAL8Vector(hP21_phase);
  if (hP33_amp)
    XLALDestroyREAL8Vector(hP33_amp);
  if (hP33_phase)
    XLALDestroyREAL8Vector(hP33_phase);
  if (hP44_amp)
    XLALDestroyREAL8Vector(hP44_amp);
  if (hP44_phase)
    XLALDestroyREAL8Vector(hP44_phase);
  if (hP55_amp)
    XLALDestroyREAL8Vector(hP55_amp);
  if (hP55_phase)
    XLALDestroyREAL8Vector(hP55_phase);

  if (alphaJ2P)
    XLALDestroyREAL8Vector(alphaJ2P);
  if (betaJ2P)
    XLALDestroyREAL8Vector(betaJ2P);
  if (gammaJ2P)
    XLALDestroyREAL8Vector(gammaJ2P);
  if (AttachPars)
    XLALDestroyREAL8Vector(AttachPars);
  if (hJlm)
    XLALDestroySphHarmTimeSeries(hJlm);

  return hIlm_dimfull;
}

static int CAmpPhaseSequence_Init(
    CAmpPhaseSequence **campphase, /* Double pointer to campphase sequence */
    int size                       /* Size of data */
) {
  /* Check input pointer */
  if (!campphase) {
    XLALPrintError("XLAL Error - %s: input double pointer is NULL.\n",
                   __func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (*campphase) {
    XLALPrintError("XLAL Error - %s: input pointer is not NULL.\n", __func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* Allocate structure */
  *campphase = XLALMalloc(sizeof(CAmpPhaseSequence));

  /* Allocate storage for data and initialize to 0 */
  if (!((*campphase)->xdata = XLALCreateREAL8Vector(size))) {
    XLALPrintError("XLAL Error - %s: failed to create REAL8Vector xdata.\n",
                   __func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if (!((*campphase)->camp_real = XLALCreateREAL8Vector(size))) {
    XLALPrintError("XLAL Error - %s: failed to create REAL8Vector camp_real.\n",
                   __func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if (!((*campphase)->camp_imag = XLALCreateREAL8Vector(size))) {
    XLALPrintError("XLAL Error - %s: failed to create REAL8Vector camp_imag.\n",
                   __func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if (!((*campphase)->phase = XLALCreateREAL8Vector(size))) {
    XLALPrintError("XLAL Error - %s: failed to create REAL8Vector phase.\n",
                   __func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  memset((*campphase)->xdata->data, 0, size * sizeof(REAL8));
  memset((*campphase)->camp_real->data, 0, size * sizeof(REAL8));
  memset((*campphase)->camp_imag->data, 0, size * sizeof(REAL8));
  memset((*campphase)->phase->data, 0, size * sizeof(REAL8));

  return XLAL_SUCCESS;
}

static int CAmpPhaseSequence_Destroy(
    CAmpPhaseSequence *campphase /* Pointer to structure to be destroyed */
) {
  /* Raise an error if NULL pointer */
  if (!campphase) {
    XLALPrintError("XLAL Error - %s: data is a NULL pointer.\n", __func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (campphase->xdata)
    XLALDestroyREAL8Vector(campphase->xdata);
  if (campphase->camp_real)
    XLALDestroyREAL8Vector(campphase->camp_real);
  if (campphase->camp_imag)
    XLALDestroyREAL8Vector(campphase->camp_imag);
  if (campphase->phase)
    XLALDestroyREAL8Vector(campphase->phase);
  XLALFree(campphase);

  return XLAL_SUCCESS;
}

static int SphHarmListEOBNonQCCoeffs_Destroy(
    SphHarmListEOBNonQCCoeffs *list /* Pointer to list to be destroyed */
) {
  SphHarmListEOBNonQCCoeffs *pop;
  while ((pop = list)) {
    if (pop->nqcCoeffs) { /* Internal EOBNonQCCoeffs is freed */
      XLALFree(pop->nqcCoeffs);
    }
    /* Go to next element */
    list = pop->next;
    /* Free structure itself */
    XLALFree(pop);
  }

  return XLAL_SUCCESS;
}

/* Note that we do NOT COPY the input mode */
/* The added data is simply passed by pointer */
/* The added data can therefore e.g. be destroyed if the list is destroyed */
static int SphHarmListEOBNonQCCoeffs_AddMode(
    SphHarmListEOBNonQCCoeffs *
        *list_prepended,       /* List structure to prepend to */
    EOBNonQCCoeffs *nqcCoeffs, /* Mode data to be added */
    UINT4 l,                   /*< Mode number l */
    INT4 m                     /*< Mode number m */
) {
  SphHarmListEOBNonQCCoeffs *list;
  /* Check if the node with this mode already exists */
  list = *list_prepended;
  while (list) {
    if (l == list->l && m == list->m) {
      break;
    }
    list = list->next;
  }
  if (list) { /* We don't allow for the case where the mode already exists in
                 the list*/
    XLALPrintError("XLAL Error - %s: tried to add an already existing mode to "
                   "a SphHarmListCAmpPhaseSequence.\n",
                   __func__);
    XLAL_ERROR(XLAL_EINVAL);
  } else {
    list = XLALMalloc(sizeof(SphHarmListEOBNonQCCoeffs));
  }
  list->l = l;
  list->m = m;
  if (nqcCoeffs) {
    list->nqcCoeffs = nqcCoeffs;
  } else {
    list->nqcCoeffs = NULL;
  }
  if (*list_prepended) {
    list->next = *list_prepended;
  } else {
    list->next = NULL;
  }
  *list_prepended = list;

  return XLAL_SUCCESS;
}

static SphHarmListEOBNonQCCoeffs *SphHarmListEOBNonQCCoeffs_GetMode(
    SphHarmListEOBNonQCCoeffs
        *list, /* List structure to get a particular mode from */
    UINT4 l,   /*< Mode number l */
    INT4 m     /*< Mode number m */
) {
  if (!list)
    return NULL;

  SphHarmListEOBNonQCCoeffs *itr = list;
  while (itr->l != l || itr->m != m) {
    itr = itr->next;
    if (!itr)
      return NULL;
  }
  /* Return a pointer to a SphHarmListCAmpPhaseSequence */
  return itr;
}

static int SphHarmListCAmpPhaseSequence_Destroy(
    SphHarmListCAmpPhaseSequence *list /* Pointer to list to be destroyed */
) {
  SphHarmListCAmpPhaseSequence *pop;
  while ((pop = list)) {
    if (pop->campphase) { /* Internal CAmpPhaseSequence is freed */
      if (CAmpPhaseSequence_Destroy(pop->campphase) == XLAL_FAILURE) {
        XLALPrintError(
            "XLAL Error - %s: failure in CAmpPhaseSequence_Destroy.\n",
            __func__);
        XLAL_ERROR(XLAL_EFUNC);
      }
    }
    /* Go to next element */
    list = pop->next;
    /* Free structure itself */
    XLALFree(pop);
  }

  return XLAL_SUCCESS;
}

/* Note that we do NOT COPY the input mode */
/* The added data is simply passed by pointer */
/* The added data can therefore e.g. be destroyed if the list is destroyed */
static int SphHarmListCAmpPhaseSequence_AddMode(
    SphHarmListCAmpPhaseSequence *
        *list_prepended,          /* List structure to prepend to */
    CAmpPhaseSequence *campphase, /* Mode data to be added */
    UINT4 l,                      /*< Mode number l */
    INT4 m                        /*< Mode number m */
) {
  SphHarmListCAmpPhaseSequence *list;
  /* Check if the node with this mode already exists */
  list = *list_prepended;
  while (list) {
    if (l == list->l && m == list->m) {
      break;
    }
    list = list->next;
  }
  if (list) { /* We don't allow for the case where the mode already exists in
                 the list*/
    XLALPrintError("XLAL Error - %s: tried to add an already existing mode to "
                   "a SphHarmListCAmpPhaseSequence.\n",
                   __func__);
    XLAL_ERROR(XLAL_EINVAL);
  } else {
    list = XLALMalloc(sizeof(SphHarmListCAmpPhaseSequence));
  }
  list->l = l;
  list->m = m;
  if (campphase) {
    list->campphase = campphase;
  } else {
    list->campphase = NULL;
  }
  if (*list_prepended) {
    list->next = *list_prepended;
  } else {
    list->next = NULL;
  }
  *list_prepended = list;

  return XLAL_SUCCESS;
}

static SphHarmListCAmpPhaseSequence *SphHarmListCAmpPhaseSequence_GetMode(
    SphHarmListCAmpPhaseSequence
        *list, /* List structure to get a particular mode from */
    UINT4 l,   /*< Mode number l */
    INT4 m     /*< Mode number m */
) {
  if (!list)
    return NULL;

  SphHarmListCAmpPhaseSequence *itr = list;
  while (itr->l != l || itr->m != m) {
    itr = itr->next;
    if (!itr)
      return NULL;
  }
  /* Return a pointer to a SphHarmListCAmpPhaseSequence */
  return itr;
}

static int SEOBdynamics_Destroy(SEOBdynamics *seobdynamics) {
  XLALDestroyREAL8Array(seobdynamics->array);
  XLALFree(seobdynamics);

  return XLAL_SUCCESS;
}

static int SEOBdynamics_Init(
    SEOBdynamics **seobdynamics, /**<< Output: pointer to the SOBdynamics */
    UINT4 retLen /**<< Input: length of dynamics data to allocate */
) {
  /* Check that the input double pointer is not NULL */
  if (!seobdynamics) {
    XLALPrintError("XLAL Error - %s: seobdynamics is a NULL double pointer.\n",
                   __func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* If double pointer points to an existing struct, destroy it first */
  if ((*seobdynamics)) {
    SEOBdynamics_Destroy(*seobdynamics);
  }

  /* Allocate struct */
  if (!(*seobdynamics = XLALMalloc(sizeof(SEOBdynamics)))) {
    XLALPrintError("XLAL Error - %s: failed to allocate struct SEOBdynamics.\n",
                   __func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }

  /* Length of data */
  (*seobdynamics)->length = retLen;

  /* Allocate array for the data */
  if (!((*seobdynamics)->array =
            XLALCreateREAL8ArrayL(2, v4PdynamicsVariables, retLen))) {
    XLALPrintError(
        "XLAL Error - %s: failed to allocate REAL8Array seobdynamics->array.\n",
        __func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }

  /* Set array pointers corresponding to the data vectors (successive vectors of
   * length retLen in the 1D array data) */
  (*seobdynamics)->tVec = (*seobdynamics)->array->data;
  (*seobdynamics)->posVecx = (*seobdynamics)->array->data + 1 * retLen;
  (*seobdynamics)->posVecy = (*seobdynamics)->array->data + 2 * retLen;
  (*seobdynamics)->posVecz = (*seobdynamics)->array->data + 3 * retLen;
  (*seobdynamics)->momVecx = (*seobdynamics)->array->data + 4 * retLen;
  (*seobdynamics)->momVecy = (*seobdynamics)->array->data + 5 * retLen;
  (*seobdynamics)->momVecz = (*seobdynamics)->array->data + 6 * retLen;
  (*seobdynamics)->s1Vecx = (*seobdynamics)->array->data + 7 * retLen;
  (*seobdynamics)->s1Vecy = (*seobdynamics)->array->data + 8 * retLen;
  (*seobdynamics)->s1Vecz = (*seobdynamics)->array->data + 9 * retLen;
  (*seobdynamics)->s2Vecx = (*seobdynamics)->array->data + 10 * retLen;
  (*seobdynamics)->s2Vecy = (*seobdynamics)->array->data + 11 * retLen;
  (*seobdynamics)->s2Vecz = (*seobdynamics)->array->data + 12 * retLen;
  (*seobdynamics)->phiDMod = (*seobdynamics)->array->data + 13 * retLen;
  (*seobdynamics)->phiMod = (*seobdynamics)->array->data + 14 * retLen;
  (*seobdynamics)->velVecx = (*seobdynamics)->array->data + 15 * retLen;
  (*seobdynamics)->velVecy = (*seobdynamics)->array->data + 16 * retLen;
  (*seobdynamics)->velVecz = (*seobdynamics)->array->data + 17 * retLen;
  (*seobdynamics)->polarrVec = (*seobdynamics)->array->data + 18 * retLen;
  (*seobdynamics)->polarphiVec = (*seobdynamics)->array->data + 19 * retLen;
  (*seobdynamics)->polarprVec = (*seobdynamics)->array->data + 20 * retLen;
  (*seobdynamics)->polarpphiVec = (*seobdynamics)->array->data + 21 * retLen;
  (*seobdynamics)->omegaVec = (*seobdynamics)->array->data + 22 * retLen;
  (*seobdynamics)->s1dotZVec = (*seobdynamics)->array->data + 23 * retLen;
  (*seobdynamics)->s2dotZVec = (*seobdynamics)->array->data + 24 * retLen;
  (*seobdynamics)->hamVec = (*seobdynamics)->array->data + 25 * retLen;

  return XLAL_SUCCESS;
}

/**
 * Functions to calculate symmetrized and antisymmetrized combinations
 * of the dimensionless spins projected on the radiation frame Z-axis (L or LN)
 */
static REAL8 SEOBCalculateChiS(REAL8 chi1dotZ, REAL8 chi2dotZ) {
  return 0.5 * (chi1dotZ + chi2dotZ);
}
static REAL8 SEOBCalculateChiA(REAL8 chi1dotZ, REAL8 chi2dotZ) {
  return 0.5 * (chi1dotZ - chi2dotZ);
}

/**
 * Function to calculate tplspin
 * See discussion below Eq. 4 of PRD 89, 061502(R) [arXiv:1311.2544] (2014)
 */
static REAL8 SEOBCalculatetplspin(REAL8 m1, REAL8 m2, REAL8 eta, REAL8 chi1dotZ,
                                  REAL8 chi2dotZ, INT4 SpinAlignedEOBversion) {
  REAL8 chiS, chiA, tplspin;
  chiS = SEOBCalculateChiS(chi1dotZ, chi2dotZ);
  chiA = SEOBCalculateChiA(chi1dotZ, chi2dotZ);

  switch (SpinAlignedEOBversion) {
  case 1:
    /* See below Eq. 17 of PRD 86, 024011 (2012) [arXiv:1202.0790] */
    tplspin = 0.0;
    break;
  case 2:
  case 4:
    /* See below Eq. 4 of PRD 89, 061502(R) (2014) [arXiv:1311.2544] */
    tplspin = (1. - 2. * eta) * chiS + (m1 - m2) / (m1 + m2) * chiA;
    break;
  default:
    XLALPrintError("XLAL Error - %s: Unknown SEOBNR version!\nAt present only "
                   "v1, v2 and v4 are available.\n",
                   __func__);
    XLAL_ERROR(XLAL_EINVAL);
    break;
  }
  return tplspin;
}

/**
 * Function to calculate normalized spin of the deformed-Kerr background in
 * SEOBNRv1. Eq. 5.2 of Barausse and Buonanno PRD 81, 084024 (2010) [arXiv:0912.3517].
 * Identical to XLALSimIMRSpinEOBCalculateSigmaKerr, except that the input spins are in
 * units of mTotal^2
 */
static int SEOBCalculateSigmaKerr(
    REAL8Vector *sigmaKerr, /**<< OUTPUT, normalized (to total mass) spin of
                               deformed-Kerr */
    REAL8Vector *s1,        /**<< spin vector 1, in units of mTotal^2 */
    REAL8Vector *s2         /**<< spin vector 2, in units of mTotal^2 */
) {
  for (UINT4 i = 0; i < 3; i++) {
    sigmaKerr->data[i] = (s1->data[i] + s2->data[i]);
  }
  return XLAL_SUCCESS;
}

/**
 * Function to calculate normalized spin of the test particle in SEOBNRv1.
 * Eq. 5.3 of Barausse and Buonanno PRD 81, 084024 (2010) [arXiv:0912.3517].
 * Identical to XLALSimIMRSpinEOBCalculateSigmaStar, except that the input spins
 * are in units of mTotal^2
 */
static int SEOBCalculateSigmaStar(
    REAL8Vector *sigmaStar, /**<< OUTPUT, normalized (to total mass) spin of
                               test particle */
    REAL8 mass1,            /**<< mass 1 */
    REAL8 mass2,            /**<< mass 2 */
    REAL8Vector *s1,        /**<< spin vector 1, in units of mTotal^2 */
    REAL8Vector *s2         /**<< spin vector 2, in units of mTotal^2 */
) {
  for (UINT4 i = 0; i < 3; i++) {
    sigmaStar->data[i] =
        (mass2 / mass1 * s1->data[i] + mass1 / mass2 * s2->data[i]);
  }
  return XLAL_SUCCESS;
}

/**
 * This function computes quantities (polardynamics, omega, s1dotZ, s2dotZ,
 * hamiltonian) derived from the dynamics as output by the integrator, and
 * returns a SEOBdynamics struct. Two choices for Z: L or LN. Note: this
 * function also applies when the spins are almost aligned and v4 is used.
 */
static int SEOBComputeExtendedSEOBdynamics(
    SEOBdynamics **seobdynamics, /**<< Output, double pointer to SEOBdynamics
                                    struct. If points to an existing struct, the
                                    latter will be destroyed */
    REAL8Array *dynamics, /**<< Input, array containing the dynamics as output
                             by the integrator */
    UINT4 retLen,         /**<< Input, length of the dynamics */
    SpinEOBParams *seobParams, /**<< SEOB parameters */
    flagSEOBNRv4P_hamiltonian_derivative
        flagHamiltonianDerivative, /**<< flag to choose between numerical and
                                      analytical Hamiltonian derivatives */
    flagSEOBNRv4P_Zframe
        flagZframe /**<< flag to choose Z direction of the frame, LN or L */
) {

  /* Create structure SEOBdynamics */
  SEOBdynamics_Init(seobdynamics, retLen);
  SEOBdynamics *seobdyn = *seobdynamics;

  /* Local variables */
  REAL8 rvec[3] = {0, 0, 0};
  REAL8 pvec[3] = {0, 0, 0};
  REAL8 spin1vec[3] = {0, 0, 0}; /* in units of mTotal^2 */
  REAL8 spin2vec[3] = {0, 0, 0}; /* in units of mTotal^2 */
  REAL8 rdotvec[3] = {0, 0, 0};
  REAL8 rcrossrdot[3] = {0, 0, 0};
  REAL8 rcrossp[3] = {0, 0, 0};
  REAL8 LNhat[3] = {0, 0, 0};
  REAL8 Lhat[3] = {0, 0, 0};
  REAL8 polarr, polarphi, polarpr, polarpphi, omega, s1dotZ, s2dotZ, ham;

  /* Allocate temporary vectors values, dvalues */
  REAL8Vector *values = NULL;
  REAL8Vector *dvalues = NULL;
  if (!(values = XLALCreateREAL8Vector(14)) ||
      !(dvalues = XLALCreateREAL8Vector(14))) {
    XLALPrintError(
        "XLAL Error - %s: failed to allocate REAL8Vector values, dvalues.\n",
        __func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  memset(values->data, 0, (values->length) * sizeof(REAL8));
  memset(dvalues->data, 0, (dvalues->length) * sizeof(REAL8));

  /* Masses and coeffs from SpinEOBParams */
  REAL8 m1 = seobParams->eobParams->m1;
  REAL8 m2 = seobParams->eobParams->m2;
  REAL8 eta = seobParams->eobParams->eta;
  SpinEOBHCoeffs *seobCoeffs = seobParams->seobCoeffs;

  /* Copying directly the dynamics data in seobdynamics - 15 vectors of length
   * retLen */
  memcpy(seobdyn->array->data, dynamics->data, 15 * retLen * sizeof(REAL8));

  /* We will need a vector structure for sigmaKerr and sigmaStar */
  REAL8Vector *sigmaStar = NULL, *sigmaKerr = NULL;
  sigmaStar = XLALCreateREAL8Vector(3);
  sigmaKerr = XLALCreateREAL8Vector(3);
  memset(sigmaStar->data, 0, 3 * sizeof(REAL8));
  memset(sigmaKerr->data, 0, 3 * sizeof(REAL8));

  /* Loop to compute the derived quantities from the dynamics */
  UINT4 i, j;
  for (i = 0; i < retLen; i++) {

    /* Copy dynamics values in the temporary vector values -- time excepted */
    for (j = 0; j < 14; j++) {
      values->data[j] = dynamics->data[i + (j + 1) * retLen];
    }

    /* Computing velocity from Hamiltonian derivatives */
    if (flagHamiltonianDerivative ==
        FLAG_SEOBNRv4P_HAMILTONIAN_DERIVATIVE_ANALYTICAL) {
      if (XLALSpinPrecHcapRvecDerivative_exact(0, values->data, dvalues->data,
                                               (void *)seobParams) ==
          XLAL_FAILURE) {
        XLALPrintError("XLAL Error - %s: failure in "
                       "XLALSpinPrecHcapRvecDerivative_exact.\n",
                       __func__);
        XLAL_ERROR(XLAL_EDOM);
      }
    } else if (flagHamiltonianDerivative ==
               FLAG_SEOBNRv4P_HAMILTONIAN_DERIVATIVE_NUMERICAL) {
      if (XLALSpinPrecHcapRvecDerivative(0, values->data, dvalues->data,
                                         (void *)seobParams) == XLAL_FAILURE) {
        XLALPrintError(
            "XLAL Error - %s: failure in XLALSpinPrecHcapRvecDerivative.\n",
            __func__);
        XLAL_ERROR(XLAL_EDOM);
      }
    } else {
      XLALPrintError(
          "XLAL Error - %s: flagHamiltonianDerivative not recognized.\n",
          __func__);
      XLAL_ERROR(XLAL_EINVAL);
    }

    /* Compute omega and LNhat */
    for (j = 0; j < 3; j++) {
      rvec[j] = values->data[j];
      pvec[j] = values->data[3 + j];
      spin1vec[j] = values->data[6 + j];
      spin2vec[j] = values->data[9 + j];
      rdotvec[j] = dvalues->data[j];
    }
    cross_product(rvec, pvec, rcrossp);
    cross_product(rvec, rdotvec, rcrossrdot);
    REAL8 rcrossrdotNorm = sqrt(inner_product(rcrossrdot, rcrossrdot));
    for (j = 0; j < 3; j++) {
      LNhat[j] = rcrossrdot[j] / rcrossrdotNorm;
    }

    /* Polar dynamics */
    polarr = sqrt(inner_product(rvec, rvec));
    polarpr = inner_product(rvec, pvec) / polarr;
    polarphi = values->data[12] + values->data[13];
    REAL8 magL = sqrt(inner_product(rcrossp, rcrossp));
    for (j = 0; j < 3; j++) {
      Lhat[j] = rcrossp[j] / magL;
    }
    polarpphi = magL;

    /* Computing omega */
    omega = rcrossrdotNorm / (polarr * polarr);

    /* Projections of the spin vectors onto the Z-axis of the precessing frame,
     * L or LN */
    if (flagZframe == FLAG_SEOBNRv4P_ZFRAME_L) {
      s1dotZ = inner_product(spin1vec, Lhat);
      s2dotZ = inner_product(spin2vec, Lhat);
    } else if (flagZframe == FLAG_SEOBNRv4P_ZFRAME_LN) {
      s1dotZ = inner_product(spin1vec, LNhat);
      s2dotZ = inner_product(spin2vec, LNhat);
    } else {
      XLALPrintError("XLAL Error - %s: flagZframe not recognized.\n", __func__);
      XLAL_ERROR(XLAL_EINVAL);
    }

    /* Compute Hamiltonian */
    UINT4 SpinAlignedEOBversion = seobParams->seobCoeffs->SpinAlignedEOBversion;
    REAL8Vector cartPosVec, cartMomVec, s1Vec, s2Vec;
    cartPosVec.length = cartMomVec.length = s1Vec.length = s2Vec.length = 3;
    cartPosVec.data = rvec;
    cartMomVec.data = pvec;
    s1Vec.data = spin1vec; /* in units of mTotal^2 */
    s2Vec.data = spin2vec; /* in units of mTotal^2 */
    SEOBCalculateSigmaStar(sigmaStar, m1, m2, &s1Vec, &s2Vec);
    SEOBCalculateSigmaKerr(sigmaKerr, &s1Vec, &s2Vec);

    // Compute the augmented spin used in the Hamiltonian calibration
    // coefficients. See LIGO-T1900601-v1.

    REAL8 tempS1_p = inner_product(s1Vec.data, Lhat);
    REAL8 tempS2_p = inner_product(s2Vec.data, Lhat);
    REAL8 S1_perp[3] = {0, 0, 0};
    REAL8 S2_perp[3] = {0, 0, 0};
    for (UINT4 jj = 0; jj < 3; jj++) {
      S1_perp[jj] = spin1vec[jj] - tempS1_p * Lhat[jj];
      S2_perp[jj] = spin2vec[jj] - tempS2_p * Lhat[jj];
    }
    UNUSED REAL8 sKerr_norm =
        sqrt(inner_product(sigmaKerr->data, sigmaKerr->data));
    REAL8 S_con = 0.0;
    if (sKerr_norm > 1e-6) {
      S_con = sigmaKerr->data[0] * Lhat[0] + sigmaKerr->data[1] * Lhat[1] +
              sigmaKerr->data[2] * Lhat[2];
      S_con /= (1 - 2 * eta);
      S_con += (inner_product(S1_perp, sigmaKerr->data) +
                inner_product(S2_perp, sigmaKerr->data)) /
               sKerr_norm / (1 - 2 * eta) / 2.;
    }

    REAL8 a = sqrt(inner_product(sigmaKerr->data, sigmaKerr->data));
    if (XLALSimIMRCalculateSpinPrecEOBHCoeffs_v2(
            seobCoeffs, eta, a, S_con, SpinAlignedEOBversion) == XLAL_FAILURE) {
      XLAL_PRINT_ERROR("XLAL Error: Something went wrong evaluating "
                       "XLALSimIMRCalculateSpinPrecEOBHCoeffs_v2 in step %d of "
                       "the main loop.\n",
                       i);
      XLAL_ERROR(XLAL_EFUNC);
    }
    ham = XLALSimIMRSpinPrecEOBHamiltonian(eta, &cartPosVec, &cartMomVec,
                                           &s1Vec, &s2Vec, sigmaKerr, sigmaStar,
                                           seobParams->tortoise, seobCoeffs);
    /* Output values in seobdynamics */
    seobdyn->velVecx[i] = rdotvec[0];
    seobdyn->velVecy[i] = rdotvec[1];
    seobdyn->velVecz[i] = rdotvec[2];
    seobdyn->polarrVec[i] = polarr;
    seobdyn->polarphiVec[i] = polarphi;
    seobdyn->polarprVec[i] = polarpr;
    seobdyn->polarpphiVec[i] = polarpphi;
    seobdyn->omegaVec[i] = omega;
    seobdyn->s1dotZVec[i] = s1dotZ;
    seobdyn->s2dotZVec[i] = s2dotZ;
    seobdyn->hamVec[i] = ham;
  }

  /* Cleanup */
  XLALDestroyREAL8Vector(values);
  XLALDestroyREAL8Vector(dvalues);
  XLALDestroyREAL8Vector(sigmaStar);
  XLALDestroyREAL8Vector(sigmaKerr);

  return XLAL_SUCCESS;
}

/**
 * This function computes initial conditions for SEOBNRv4P.
 */
static int SEOBInitialConditions(
    REAL8Vector **ICvalues, /**<< Output: vector with initial conditions */
    REAL8 MfMin,       /**<< Input: dimensionless initial frequency (in units of
                          1/mTotal) */
    REAL8 m1,          /**<< Input: mass 1 (solar masses) */
    REAL8 m2,          /**<< Input: mass 2 (solar masses) */
    REAL8Vector *chi1, /**<< Input: dimensionless spin 1 (in units of m1^2) */
    REAL8Vector *chi2, /**<< Input: dimensionless spin 2 (in units of m2^2) */
    SpinEOBParams *seobParams, /**<< SEOB params */
    flagSEOBNRv4P_hamiltonian_derivative
        flagHamiltonianDerivative /**<< flag to decide wether to use analytical
                                     or numerical derivatives */
) {
  UINT4 j;

  /* Check that the input double pointer is not NULL */
  if (!ICvalues) {
    XLALPrintError(
        "XLAL Error - %s: pointer to REAL8Vector ICvalues is NULL.\n",
        __func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* Allocate vector for initial conditions at fMin */
  if (!(*ICvalues = XLALCreateREAL8Vector(14))) {
    XLALPrintError(
        "XLAL Error - %s: failed to allocate REAL8Vector ICvalues.\n",
        __func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  memset((*ICvalues)->data, 0, ((*ICvalues)->length) * sizeof(REAL8));

  /* Masses in solar masses */
  REAL8 eta = m1 * m2 / (m1 + m2) / (m1 + m2);
  /* Flags */
  UINT4 SpinAlignedEOBversion = seobParams->seobCoeffs->SpinAlignedEOBversion;
  UINT4 SpinsAlmostAligned = seobParams->alignedSpins;

  /* Needed for the interface of the old code */
  REAL8 fMin = MfMin / (m1 + m2) / LAL_MTSUN_SI;
  REAL8 inc =
      0.; /* Not clear why XLALSimIMRSpinEOBInitialConditions,
             XLALSimIMRSpinEOBInitialConditionsPrec need an inclination */

  /* Flag for numerical or analytical derivatives */
  INT4 use_optimized = 0;
  if (flagHamiltonianDerivative ==
      FLAG_SEOBNRv4P_HAMILTONIAN_DERIVATIVE_ANALYTICAL)
    use_optimized = 1;
  else if (flagHamiltonianDerivative ==
           FLAG_SEOBNRv4P_HAMILTONIAN_DERIVATIVE_NUMERICAL)
    use_optimized = 0;
  else {
    XLALPrintError(
        "XLAL Error - %s: flagHamiltonianDerivative not recognized.\n",
        __func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* XLALSimIMRSpinEOBInitialConditions takes as input dimensionfull spins in
   * units of solar mass square */
  REAL8 mSpin1data[3] = {0., 0., 0.};
  REAL8 mSpin2data[3] = {0., 0., 0.};

  if (SpinsAlmostAligned) {
    /* Here the Z-frame axis is the original z, also direction of L and LN - no
     * precession */
    REAL8 chi1dotZ = chi1->data[2];
    REAL8 chi2dotZ = chi2->data[2];
    REAL8 chiS = SEOBCalculateChiS(chi1dotZ, chi2dotZ);
    REAL8 chiA = SEOBCalculateChiA(chi1dotZ, chi2dotZ);
    REAL8 tplspin = SEOBCalculatetplspin(m1, m2, eta, chi1dotZ, chi2dotZ,
                                         SpinAlignedEOBversion);
    if (XLALSimIMREOBCalcSpinFacWaveformCoefficients(
            seobParams->eobParams->hCoeffs, seobParams, m1, m2, eta, tplspin,
            chiS, chiA, SpinAlignedEOBversion) ==
        XLAL_FAILURE) /* This function returns XLAL_SUCCESS or calls XLAL_ERROR(
                         XLAL_EINVAL ) */
    {
      XLALPrintError("XLAL Error - %s: failure in "
                     "XLALSimIMREOBCalcSpinFacWaveformCoefficients.\n",
                     __func__);
      XLAL_ERROR(XLAL_EFUNC);
    }
    mSpin1data[2] = chi1->data[2] * m1 * m1;
    mSpin2data[2] = chi2->data[2] * m2 * m2;
    if (XLALSimIMRSpinEOBInitialConditions(*ICvalues, m1, m2, fMin, inc,
                                           mSpin1data, mSpin2data, seobParams,
                                           use_optimized) ==
        XLAL_FAILURE) /* This function returns XLAL_SUCCESS or calls XLAL_ERROR
                         with XLAL_EINVAL, XLAL_ENOMEM, or XLAL_EMAXITER */
    {
      XLALPrintError(
          "XLAL Error - %s: failure in XLALSimIMRSpinEOBInitialConditions.\n",
          __func__);
      XLAL_ERROR(XLAL_EFUNC);
    }
  } else {
    for (j = 0; j < 3; j++) {
      mSpin1data[j] = chi1->data[j] * m1 * m1;
      mSpin2data[j] = chi2->data[j] * m2 * m2;
    }
    /* This function returns XLAL_SUCCESS or calls XLAL_ERROR with XLAL_EINVAL,
     * XLAL_ENOMEM, or XLAL_EMAXITER */
    if (XLALSimIMRSpinEOBInitialConditionsPrec(
            *ICvalues, m1, m2, fMin, inc, mSpin1data, mSpin2data, seobParams,
            use_optimized) == XLAL_FAILURE) {
      XLALPrintError("XLAL Error - %s: failure in "
                     "XLALSimIMRSpinEOBInitialConditionsPrec.\n",
                     __func__);
      XLAL_ERROR(XLAL_EFUNC);
    }
  }

  /* Initial phases are set to 0 */
  (*ICvalues)->data[12] = 0.;
  (*ICvalues)->data[13] = 0.;

  return XLAL_SUCCESS;
}

/**
 * This function converts a spin-aligned dynamics as output by the Runge-Kutta
 * integrator to a generic-spin dynamics. Spin-aligned dynamics format: t, r,
 * phi, pr, pphi Generic-spin dynamics format: t, x, y, z, px, py, pz, s1x, s1y,
 * s1z, s2x, s2y, s2z, phiMod, phiDMod
 */
static int SEOBConvertSpinAlignedDynamicsToGenericSpins(
    REAL8Array **dynamics, /**<< Output: pointer to array for the generic-spin
                              dynamics */
    REAL8Array *dynamics_spinaligned, /**<< Input: array for the aligned-spin
                                         dynamics */
    UINT4 retLen,                     /**<< Input: length of dynamics */
    REAL8 chi1, /**<< Input: spin 1 aligned component (dimensionless) */
    REAL8 chi2, /**<< Input: spin 2 aligned component (dimensionless) */
    SpinEOBParams *seobParams /**<< SEOB params */
) {
  UINT4 i;

  /* Masses */
  REAL8 m1 = seobParams->eobParams->m1;
  REAL8 m2 = seobParams->eobParams->m2;
  REAL8 mTotal = m1 + m2;

  /* Create output dynamics */
  *dynamics = XLALCreateREAL8ArrayL(2, 15, retLen);

  /* Convert the spin-aligned dynamics to a generic-spins dynamics */
  REAL8Vector tVec, rVec, phiVec, prVec, pPhiVec;
  tVec.length = rVec.length = phiVec.length = prVec.length = pPhiVec.length =
      retLen;
  tVec.data = dynamics_spinaligned->data;
  rVec.data = dynamics_spinaligned->data + retLen;
  phiVec.data = dynamics_spinaligned->data + 2 * retLen;
  prVec.data = dynamics_spinaligned->data + 3 * retLen;
  pPhiVec.data = dynamics_spinaligned->data + 4 * retLen;
  for (i = 0; i < retLen; i++) {
    (*dynamics)->data[i] = tVec.data[i];
    (*dynamics)->data[retLen + i] = rVec.data[i] * cos(phiVec.data[i]);
    (*dynamics)->data[2 * retLen + i] = rVec.data[i] * sin(phiVec.data[i]);
    (*dynamics)->data[3 * retLen + i] = 0.;
    (*dynamics)->data[4 * retLen + i] =
        prVec.data[i] * cos(phiVec.data[i]) -
        pPhiVec.data[i] / rVec.data[i] * sin(phiVec.data[i]);
    (*dynamics)->data[5 * retLen + i] =
        prVec.data[i] * sin(phiVec.data[i]) +
        pPhiVec.data[i] / rVec.data[i] * cos(phiVec.data[i]);
    (*dynamics)->data[6 * retLen + i] = 0.;
    (*dynamics)->data[7 * retLen + i] = 0.;
    (*dynamics)->data[8 * retLen + i] = 0.;
    (*dynamics)->data[9 * retLen + i] = chi1 * (m1 * m1 / mTotal / mTotal);
    (*dynamics)->data[10 * retLen + i] = 0.;
    (*dynamics)->data[11 * retLen + i] = 0.;
    (*dynamics)->data[12 * retLen + i] = chi2 * (m2 * m2 / mTotal / mTotal);
    (*dynamics)->data[13 * retLen + i] = phiVec.data[i];
    (*dynamics)->data[14 * retLen + i] = 0.;
  }

  return XLAL_SUCCESS;
}

/**
 * This function integrates the SEOBNRv4P dynamics.
 * Output is given either on the adaptive sampling coming out of the Runge Kutta
 * integrator, with no interpolation being made, or on the constant sampling
 * specified by deltaT. Either analytical or numerical derivatives of the
 * Hamiltonian are used depending on the flag flagHamiltonianDerivative.
 * Only numerical derivatives have been shown to work as of June 2019.
 * When spins are flagged as almost aligned, falls back to
 * spin-aligned dynamics.
 */
static int SEOBIntegrateDynamics(
    REAL8Array **dynamics, /**<< Output: pointer to array for the dynamics */
    UINT4 *retLenOut,      /**<< Output: length of the output dynamics */
    REAL8Vector *ICvalues, /**<< Input: vector with initial conditions */
    REAL8 EPS_ABS, /**<< Input: absolute accuracy for adaptive Runge-Kutta
                      integrator */
    REAL8 EPS_REL, /**<< Input: relative accuracy for adaptive Runge-Kutta
                      integrator */
    REAL8 deltaT,  /**<< Input: timesampling step in geometric units - when
                      flagConstantSampling is False, used internally only to
                      initialize adaptive step */
    REAL8 deltaT_min,  /**<< Input: minimal timesampling step in geometric
                          units when using adaptive steps with
                          flagConstantSampling set to False -
                          set to 0 to ignore */
    REAL8 tstart,  /**<< Input: starting time of the integration */
    REAL8 tend,    /**<< Input: max time of the integration - normally, the
                      integration stops when stopping condition is met, and this is
                      ignored */
    SpinEOBParams *seobParams,  /**<< SEOB params */
    UINT4 flagConstantSampling, /**<< flag to decide wether to use constant
                                   sampling with deltaT in output instead of
                                   adaptive sampling */
    flagSEOBNRv4P_hamiltonian_derivative
        flagHamiltonianDerivative /**<< flag to decide wether to use analytical
                                     or numerical derivatives */
) {
  UINT4 retLen;

  /* Integrator */
  LALAdaptiveRungeKuttaIntegrator *integrator = NULL;

  /* Flags */
  UINT4 SpinsAlmostAligned = seobParams->alignedSpins;

  /* Dimensions of vectors of dynamical variables to be integrated */
  UINT4 nb_Hamiltonian_variables = 14;
  UINT4 nb_Hamiltonian_variables_spinsaligned = 4;

  /* Used only when falling back to spin-aligned dynamics */
  REAL8Array *dynamics_spinaligned = NULL;

  /* Initialize internal values vector for the generic-spin case */
  REAL8Vector *values = NULL;
  if (!(values = XLALCreateREAL8Vector(nb_Hamiltonian_variables))) {
    XLALPrintError("XLAL Error - %s: failed to create REAL8Vector values.\n",
                   __func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  memcpy(values->data, ICvalues->data, values->length * sizeof(REAL8));
  /* Initialize internal values vector for the spin-aligned case -- not used in
   * the generic-spin case */
  REAL8Vector *values_spinaligned = NULL;
  if (!(values_spinaligned =
            XLALCreateREAL8Vector(nb_Hamiltonian_variables_spinsaligned))) {
    XLALPrintError(
        "XLAL Error - %s: failed to create REAL8Vector values_spinaligned.\n",
        __func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  memset(values_spinaligned->data, 0,
         values_spinaligned->length * sizeof(REAL8));

  /* Initialization of the integrator */
  if (SpinsAlmostAligned) { /* If spins are almost aligned with LNhat, use
                               SEOBNRv4 dynamics */
    /* In SEOBNRv4 the dynamical variables are r, phi, p_r^*, p_phi */

    /* Construct the initial conditions */
    REAL8 temp_r = sqrt(ICvalues->data[0] * ICvalues->data[0] +
                        ICvalues->data[1] * ICvalues->data[1] +
                        ICvalues->data[2] * ICvalues->data[2]);
    REAL8 temp_phi = ICvalues->data[12];

    values_spinaligned->data[0] = temp_r;   // General form of r
    values_spinaligned->data[1] = temp_phi; // phi
    values_spinaligned->data[2] = ICvalues->data[3] * cos(temp_phi) +
                                  ICvalues->data[4] * sin(temp_phi); // p_r^*
    values_spinaligned->data[3] =
        temp_r * (ICvalues->data[4] * cos(temp_phi) -
                  ICvalues->data[3] * sin(temp_phi)); // p_phi

    /* We have to use different stopping conditions depending
       we are in the low-sampling or high-sampling portion
       of the waveform. We can tell this apart because for the low-sampling (or
       ada sampling) we always start at t=0
     */
    if (tstart > 0) {
      // High sampling
      integrator = XLALAdaptiveRungeKutta4Init(
          nb_Hamiltonian_variables_spinsaligned, XLALSpinAlignedHcapDerivative,
          XLALSpinPrecAlignedHiSRStopCondition, EPS_ABS, EPS_REL);
    } else {
      // Low sampling
      integrator = XLALAdaptiveRungeKutta4Init(
          nb_Hamiltonian_variables_spinsaligned, XLALSpinAlignedHcapDerivative,
          XLALEOBSpinPrecAlignedStopCondition, EPS_ABS, EPS_REL);
    }
  } else {
    if (flagHamiltonianDerivative ==
        FLAG_SEOBNRv4P_HAMILTONIAN_DERIVATIVE_ANALYTICAL) {
      integrator = XLALAdaptiveRungeKutta4Init(
          nb_Hamiltonian_variables, XLALSpinPrecHcapExactDerivative,
          XLALEOBSpinPrecStopConditionBasedOnPR, EPS_ABS, EPS_REL);
    } else if (flagHamiltonianDerivative ==
               FLAG_SEOBNRv4P_HAMILTONIAN_DERIVATIVE_NUMERICAL) {
      if (tstart > 0) {
        integrator = XLALAdaptiveRungeKutta4Init(
            nb_Hamiltonian_variables, XLALSpinPrecHcapNumericalDerivative,
            XLALEOBSpinPrecStopConditionBasedOnPR, EPS_ABS, EPS_REL);
      } else {
        integrator = XLALAdaptiveRungeKutta4Init(
            nb_Hamiltonian_variables, XLALSpinPrecHcapNumericalDerivative,
            XLALEOBSpinPrecStopConditionBasedOnPR, EPS_ABS, EPS_REL);
      }

    } else {
      XLALPrintError(
          "XLAL Error - %s: flagHamiltonianDerivative not recognized.\n",
          __func__);
      XLAL_ERROR(XLAL_EINVAL);
    }
  }
  if (!integrator) {
    XLALPrintError(
        "XLAL Error - %s: failure in the initialization of the integrator.\n",
        __func__);
    XLAL_ERROR(XLAL_EDOM);
  }
  /* Ensure that integration stops ONLY when the stopping condition is True */
  integrator->stopontestonly = 1;
  /* When this option is set to 0, the integration can be exceedingly slow for
   * spin-aligned systems */
  integrator->retries = 1;

  /* Computing the dynamical evolution of the system */
  // NOTE: XLALAdaptiveRungeKutta4NoInterpolate takes an EOBversion as input.
  INT4 EOBversion = 2; // NOTE: value 3 is specific to optv3 in
                       // XLALAdaptiveRungeKutta4NoInterpolate it determines
                       // what is stored. We set it to 2.
  if (SpinsAlmostAligned) {
    /* If spins are almost aligned with LNhat, use SEOBNRv4 dynamics */
    if (!flagConstantSampling) {
      retLen = XLALAdaptiveRungeKutta4NoInterpolate(
          integrator, seobParams, values_spinaligned->data, 0., tend - tstart,
          deltaT, deltaT_min, &dynamics_spinaligned, EOBversion);
    } else {
      retLen = XLALAdaptiveRungeKutta4(
          integrator, seobParams, values_spinaligned->data, 0., tend - tstart,
          deltaT, &dynamics_spinaligned);
    }
    if ((INT4)retLen == XLAL_FAILURE) {
      XLALPrintError("XLAL Error - %s: failure in the integration of the "
                     "spin-aligned dynamics.\n",
                     __func__);
      XLAL_ERROR(XLAL_EDOM);
    }

    /* Convert the spin-aligned dynamics to a generic-spins dynamics */
    if (SEOBConvertSpinAlignedDynamicsToGenericSpins(
            dynamics, dynamics_spinaligned, retLen, seobParams->chi1,
            seobParams->chi2, seobParams) == XLAL_FAILURE) {
      XLALPrintError("XLAL Error - %s: failure in "
                     "SEOBConvertSpinAlignedDynamicsToGenericSpins.\n",
                     __func__);
      XLAL_ERROR(XLAL_EDOM);
    }
  } else {
    if (!flagConstantSampling) {
      retLen = XLALAdaptiveRungeKutta4NoInterpolate(
          integrator, seobParams, values->data, 0., tend - tstart, deltaT,
          deltaT_min, dynamics, EOBversion);
    } else {
      retLen = XLALAdaptiveRungeKutta4(integrator, seobParams, values->data, 0.,
                                       tend - tstart, deltaT, dynamics);
    }
    if ((INT4)retLen == XLAL_FAILURE) {
      XLALPrintError("XLAL Error - %s: failure in the integration of the "
                     "generic-spin dynamics.\n",
                     __func__);
      XLAL_ERROR(XLAL_EDOM);
    }
  }

  // NOTE: functions like XLALAdaptiveRungeKutta4 would give nans if the times
  // do not start at 0 -- we have to adjust the starting time after integration
  /* Adjust starting time */
  for (UINT4 i = 0; i < retLen; i++)
    (*dynamics)->data[i] += tstart;

  /* Output length of dynamics */
  *retLenOut = retLen;

  /* Cleanup */
  if (dynamics_spinaligned)
    XLALDestroyREAL8Array(dynamics_spinaligned);
  XLALDestroyREAL8Vector(values_spinaligned);
  XLALDestroyREAL8Vector(values);
  XLALAdaptiveRungeKuttaFree(integrator);

  return XLAL_SUCCESS;
}

/**
 * This function generates a waveform mode for a given SEOB dynamics.
 */
// NOTE: as is written here, the step
// XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients in the loop will be repeated
// across modes -- would be more efficient to loop on modes inside the loop on
// times
static int SEOBCalculatehlmAmpPhase(
    CAmpPhaseSequence *
        *hlm, /**<< Output: hlm in complex amplitude / phase form */
    INT4 l,   /**<< Input: mode index l */
    INT4 m,   /**<< Input: mode index m */
    SEOBdynamics *seobdynamics, /**<< Input: SEOB dynamics */
    EOBNonQCCoeffs *nqcCoeffs,  /**<< Input: NQC coeffs */
    SpinEOBParams *seobParams,  /**<< SEOB params */
    // UINT4 SpinsAlmostAligned, /**<< flag to decide wether to fall back to
    // aligned spins  */
    UINT4 includeNQC /**<< flag to choose wether or not to include NQC */
) {
  /* Check that the input double pointer are not NULL */
  if (!hlm) {
    XLALPrintError(
        "XLAL Error - %s: pointer to CAmpPhaseSequence hlm is NULL.\n",
        __func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* Masses */
  REAL8 m1 = seobParams->eobParams->m1;
  REAL8 m2 = seobParams->eobParams->m2;
  REAL8 deltaT = seobdynamics->tVec[1] - seobdynamics->tVec[0];
  REAL8 eta = seobParams->eobParams->eta;
  REAL8 mtot = m1 + m2;
  UINT4 SpinAlignedEOBversion = seobParams->seobCoeffs->SpinAlignedEOBversion;
  UINT4 SpinAlignedEOBversionWaveform; // RC: I use this different variable
                                       // because the PN terms in the waveform
                                       // are different from those in the flux

  /* Length of dynamics data and sampling step */
  UINT4 retLen = seobdynamics->length;

  /* Allocate structure for complex amplitude and phase */
  if (CAmpPhaseSequence_Init(hlm, retLen) == XLAL_FAILURE) {
    XLALPrintError("XLAL Error - %s: failure in CAmpPhaseSequence_Init.\n",
                   __func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }

  /* Workspace vectors */
  REAL8Vector values, polarDynamics;
  REAL8 valuesdata[14] = {0.};
  REAL8 polarDynamicsdata[4] = {0.};
  values.length = 14;
  polarDynamics.length = 4;
  values.data = valuesdata;
  polarDynamics.data = polarDynamicsdata;
  REAL8 tPeakOmega = seobParams->tPeakOmega;

  /* Calibration parameter */
  if (includeNQC == 0) {
    if (((l == 2) && (m == 1)) || ((l == 5) && (m == 5))) {
      if (XLALSimIMREOBCalcCalibCoefficientHigherModesPrec(
              seobParams, l, m, seobdynamics,
              tPeakOmega - seobdynamics->tVec[0], m1, m2,
              deltaT) == XLAL_FAILURE) {
        XLALPrintError("XLAL Error - %s: failure in "
                       "XLALSimIMREOBCalcCalibCoefficientHigherModesPrec.\n",
                       __func__);
        XLAL_ERROR(XLAL_EDOM);
      }
    }
  }

  /* Loop to compute compute amplitude and phase of the hlm mode */
  REAL8 s1dotZ, s2dotZ, chiS, chiA, tplspin, t, omega, ham, v;
  UINT4 i, j;
  for (i = 0; i < retLen; i++) {
    /* Compute waveform coefficients */
    t = seobdynamics->tVec[i];
    omega = seobdynamics->omegaVec[i];
    ham = seobdynamics->hamVec[i];
    s1dotZ = seobdynamics->s1dotZVec[i];
    s2dotZ = seobdynamics->s2dotZVec[i];
    REAL8 chi1dotZ = s1dotZ * mtot * mtot / (m1 * m1);
    REAL8 chi2dotZ = s2dotZ * mtot * mtot / (m2 * m2);
    chiS = SEOBCalculateChiS(chi1dotZ, chi2dotZ);
    chiA = SEOBCalculateChiA(chi1dotZ, chi2dotZ);

    tplspin = SEOBCalculatetplspin(m1, m2, eta, s1dotZ, s2dotZ,
                                   SpinAlignedEOBversion);
    if (SpinAlignedEOBversion == 4) {
      SpinAlignedEOBversionWaveform = v4Pwave;
    } else {
      SpinAlignedEOBversionWaveform = SpinAlignedEOBversion;
    }

    if (XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients(
            seobParams->eobParams->hCoeffs, m1, m2, eta, tplspin, chiS, chiA,
            SpinAlignedEOBversionWaveform) == XLAL_FAILURE) {
      XLALPrintError("XLAL Error - %s: failure in "
                     "XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients at step "
                     "%d of the loop.\n",
                     __func__, i);
      XLAL_ERROR(XLAL_EFUNC);
    }
    seobParams->eobParams->hCoeffs->f21v7c = seobParams->cal21;
    seobParams->eobParams->hCoeffs->f55v5c = seobParams->cal55;
    // printf("f21v7c = %.16f\n", seobParams->eobParams->hCoeffs->f21v7c);
    /* Dynamics, polar dynamics, omega */
    for (j = 0; j < 14; j++) {
      values.data[j] = seobdynamics->array->data[i + (j + 1) * retLen];
    }
    polarDynamics.data[0] = seobdynamics->polarrVec[i];
    polarDynamics.data[1] = seobdynamics->polarphiVec[i];
    polarDynamics.data[2] = seobdynamics->polarprVec[i];
    polarDynamics.data[3] = seobdynamics->polarpphiVec[i];
    v = cbrt(omega);
    COMPLEX16 hlm_val = 0.;
    if (XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform(
            &hlm_val, &polarDynamics, &values, v, ham, l, m, seobParams) ==
        XLAL_FAILURE) {
      XLALPrintError("XLAL Error - %s: failure in "
                     "XLALSimIMRSpinEOBGetPrecSpinFactorizedWaveform at step "
                     "%d of the loop.\n",
                     __func__, i);
      XLAL_ERROR(XLAL_EDOM);
    }
    /* NQC correction */
    COMPLEX16 factor_nqc = 1.;
    if (includeNQC) {
      if (XLALSimIMRSpinEOBNonQCCorrection(&factor_nqc, &values, omega,
                                           nqcCoeffs) == XLAL_FAILURE) {
        XLALPrintError(
            "XLAL Error - %s: failure in XLALSimIMRSpinEOBNonQCCorrection at "
            "step %d of the loop.\n",
            __func__, i);
        XLAL_ERROR(XLAL_EDOM);
      }
    } else
      factor_nqc = 1.;
    /* Result and output */
    COMPLEX16 hlmNQC = hlm_val * factor_nqc;
    (*hlm)->xdata->data[i] = t; /* Copy times */
    (*hlm)->camp_real->data[i] = cabs(hlmNQC);
    (*hlm)->camp_imag->data[i] = 0.; /* We use only real amplitudes */
    (*hlm)->phase->data[i] = carg(hlmNQC);
  }

  /* Unwrap the phase vector, in place */
  XLALREAL8VectorUnwrapAngle((*hlm)->phase, (*hlm)->phase);

  return XLAL_SUCCESS;
}

/**
 * This function generates all waveform modes as a list for a given SEOB
 * dynamics.
 */
static int SEOBCalculateSphHarmListhlmAmpPhase(
    SphHarmListCAmpPhaseSequence *
        *listhlm,               /**<< Output: list of modes for hlm */
    INT4 modes[][2],            /**<< Input: array of modes (l,m) */
    UINT4 nmodes,               /**<< Input: number of modes (l,m) */
    SEOBdynamics *seobdynamics, /**<< Input: SEOB dynamics */
    SphHarmListEOBNonQCCoeffs *listnqcCoeffs, /**<< Input: list of NQCs */
    SpinEOBParams *seobParams,                /**<< SEOB params */
    UINT4 flagNQC /**<< flag to choose wether or not to include NQC */
) {
  /* Read version of SEOB to be used */
  UINT4 SpinAlignedEOBversion = seobParams->seobCoeffs->SpinAlignedEOBversion;

  /* Flag for inclusion of NQC, useful for higher modes that have no NQC
   * implemented for some SpinAlignedEOBversion */
  UINT4 includeNQC = 0;

  /* Loop over modes */
  for (UINT4 nmode = 0; nmode < nmodes; nmode++) {

    INT4 l = modes[nmode][0];
    INT4 m = modes[nmode][1];

    if ((!(l == 2 && m == 2)) && (SpinAlignedEOBversion == 3)) {
      includeNQC = 0; /* For HM beyond 22, no NQC available for
                         SpinAlignedEOBversion==3 */
    } else
      includeNQC = flagNQC;

    EOBNonQCCoeffs *nqcCoeffs =
        SphHarmListEOBNonQCCoeffs_GetMode(listnqcCoeffs, l, m)->nqcCoeffs;

    CAmpPhaseSequence *hlm = NULL;
    SEOBCalculatehlmAmpPhase(&hlm, l, m, seobdynamics, nqcCoeffs, seobParams,
                             includeNQC);

    SphHarmListCAmpPhaseSequence_AddMode(listhlm, hlm, l, m);
  }

  return XLAL_SUCCESS;
}

/**
 * This function finds the peak of omega.
 * Note that there are various possibilities as of what is returned if
 * tPeakOmega is not found at first. In particular, by default,
 * if there is no peak found, the last point of the dynamics is used.
 */
static int SEOBLocateTimePeakOmega(
    REAL8 *tPeakOmega, /**<< Output: time of peak of Omega if found (see inside
        XLALSimLocateOmegaTime for what is returned otherwise) */
    INT4 *foundPeakOmega, /**<< Output: flag indicating wether tPeakOmega has
                 been found */
    UNUSED REAL8Array *dynamics,      /**<< Input: array for dynamics */
    SEOBdynamics *seobdynamics,       /**<< Input: SEOB dynamics object */
    UINT4 retLen,                     /**<< Input: length of dynamics */
    UNUSED SpinEOBParams *seobParams, /**<< SEOB params */
    UNUSED flagSEOBNRv4P_hamiltonian_derivative
        flagHamiltonianDerivative /**<< flag to decide wether to use analytical
                                     or numerical derivatives */
) {

  REAL8Vector tVec;
  tVec.length = retLen;
  tVec.data = seobdynamics->tVec;
  REAL8Vector omegaVec;
  omegaVec.length = retLen;
  omegaVec.data = seobdynamics->omegaVec;
  XLALEOBFindRobustPeak(tPeakOmega, &tVec, &omegaVec, 3);
  *foundPeakOmega = 1;

  /*
  // Uncomment this block to use old-style finding of the peak of Omega
  UNUSED REAL8 m1 = seobParams->eobParams->m1;
  UNUSED REAL8 m2 = seobParams->eobParams->m2;
  REAL8Vector radiusVec;
  radiusVec.length = retLen;
  radiusVec.data = seobdynamics->polarrVec;
  // Number of variables in the dynamics array (besides time)
  UINT4 numdynvars = 14;
  // A parameter we don't know what to do with
  REAL8 tMaxOmega = 0;    // Will not be used, and is not returned

  // Flag for numerical or analytical derivatives
  INT4 use_optimized = 0;
  if (flagHamiltonianDerivative ==
      FLAG_SEOBNRv4P_HAMILTONIAN_DERIVATIVE_ANALYTICAL)
    use_optimized = 1;
  else if (flagHamiltonianDerivative ==
           FLAG_SEOBNRv4P_HAMILTONIAN_DERIVATIVE_NUMERICAL)
    use_optimized = 0;
  else {
    XLALPrintError(
                   "XLAL Error - %s: flagHamiltonianDerivative not
  recognized.\n",
                   __func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  //Time of peak of Omega
  *tPeakOmega = XLALSimLocateOmegaTime(
      dynamics, numdynvars, retLen, *seobParams, *(seobParams->seobCoeffs), m1,
      m2, &radiusVec, foundPeakOmega, &tMaxOmega, use_optimized);
  */

  return XLAL_SUCCESS;
}

/**
 * This function looks for the peak of a mode amplitude.
 * Note that the internals are complicated, see XLALSimLocateAmplTime.
 */
UNUSED static int SEOBLocateTimePeakModeAmp(
    REAL8
        *tPeakAmp, /**<< Output: time of peak of amplitude if found (see inside
                      XLALSimLocateAmplTime for what is returned otherwise) */
    INT4 *foundPeakAmp, /**<< Output: flag indicating wether tPeakOmega has been
                           found */
    CAmpPhaseSequence *hlm, /**<< Input: mode in complex amplitude/phase form */
    SEOBdynamics *seobdynamics, /**<< Input: SEOB dynamics object */
    UINT4 retLen                /**<< Input: length of dynamics */
) {
  /* Vectors with times and radius from dynamics */
  REAL8Vector timeVec, radiusVec;
  timeVec.length = radiusVec.length = retLen;
  radiusVec.data = seobdynamics->polarrVec;
  timeVec.data = seobdynamics->tVec;

  /* Computing mode */
  // NOTE: computing Re/Im is redundant, as XLALSimLocateAmplTime will recompute
  // the amplitude internally
  COMPLEX16Vector *hmode = NULL;
  if (!(hmode = XLALCreateCOMPLEX16Vector(retLen))) {
    XLALPrintError(
        "XLAL Error - %s: failed to allocate COMPLEX16Vector hmode.\n",
        __func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  for (UINT4 i = 0; i < retLen; i++) {
    hmode->data[i] = (hlm->camp_real->data[i] + I * hlm->camp_imag->data[i]) *
                     cexp(I * hlm->phase->data[i]);
  }

  /* Not clear what this is for, not used */
  // REAL8 tMaxAmp = 0.;

  /* Time of peak of amplitude */
  //*tPeakAmp = XLALSimLocateAmplTime( &timeVec, hmode, &radiusVec,
  // foundPeakAmp, &tMaxAmp );
  *tPeakAmp = XLALSimLocateMaxAmplTime(&timeVec, hmode, foundPeakAmp);

  /* Cleanup */
  XLALDestroyCOMPLEX16Vector(hmode);

  return XLAL_SUCCESS;
}

/**
 * This function computes all extended dynamics values at a given time by
 * interpolating the dynamics array. We build a cubic spline limited to +- 20
 * samples on each side of the time of interest.
 */
static int SEOBInterpolateDynamicsAtTime(
    REAL8Vector **seobdynamics_values, /**<< Output: pointer to vector for
                                          seobdynamics interpolated values */
    REAL8 t,                           /**<< Input: time at which to evaluate */
    SEOBdynamics *seobdynamics         /**<< Input: SEOB dynamics */
) {
  /* Create output vector */
  if (!((*seobdynamics_values) = XLALCreateREAL8Vector(v4PdynamicsVariables))) {
    XLALPrintError("XLAL Error - %s: failed to allocate REAL8Vector "
                   "seobdynamics_values.\n",
                   __func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  memset((*seobdynamics_values)->data, 0,
         ((*seobdynamics_values)->length) * sizeof(REAL8));

  /* Check that the time asked for is in range */
  UINT4 retLen = seobdynamics->length;
  REAL8 *tVec = seobdynamics->tVec;
  if ((t < tVec[0]) || (t > tVec[retLen - 1])) {
    XLALPrintError(
        "XLAL Error - %s: time for interpolation is out of range of "
        "the SEOBdynamics data. t = %.17f, tVec[0]=%.17f, tVec[-1]=%.17f\n",
        __func__, t, tVec[0], tVec[retLen - 1]);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* Get the start and end indices that we will use to interpolate */
  /* indext max index such that tVec[indext] <= t */
  UINT4 indext = 0;
  while ((indext < retLen - 1) && (tVec[indext + 1] <= t))
    indext++;
  INT4 indexstart = indext - 20 > 0 ? indext - 20 : 0;
  INT4 indexend = indext + 20 < retLen - 1 ? indext + 20 : retLen - 1;
  INT4 interp_length = indexend - indexstart + 1;
  if (interp_length <= 0) {
    XLALPrintError("XLAL Error - %s: not finding a strictly positive number of "
                   "samples for interpolation.\n",
                   __func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* Spline allocation */
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, interp_length);
  gsl_interp_accel *acc = gsl_interp_accel_alloc();

  /* Interpolate all quantities */
  (*seobdynamics_values)->data[0] = t;
  for (UINT4 j = 1; j < v4PdynamicsVariables; j++) {
    gsl_spline_init(spline, &(tVec[indexstart]),
                    &(seobdynamics->array->data[j * retLen + indexstart]),
                    interp_length);
    (*seobdynamics_values)->data[j] = gsl_spline_eval(spline, t, acc);
  }

  /* Cleanup */
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

  return XLAL_SUCCESS;
}

static int SEOBLFrameVectors(
    REAL8Vector **S1,        /**<<Output: S1 in L-n frame */
    REAL8Vector **S2,        /**<<Output: S2 in L-n frame */
    REAL8Vector *seobvalues, /**<<Input: vector of extended dynamics */
    REAL8 m1, /**<<Input: mass of the first object in solar masses */
    REAL8 m2, /**<<Input: mass of the second object in solar masses */
    const flagSEOBNRv4P_Zframe
        flagZframe /**<<Input: whether to compute the L_N or L frame */
) {

  if ((!S1) || (!S2) || (!seobvalues)) {
    XLALPrintError("Passed null pointers to SEOBLFrameVectors!\n");
    XLAL_ERROR(XLAL_ENOMEM);
  }
  REAL8 mTotal = m1 + m2;
  // Scaled masses, useful to compute spins in sane units
  REAL8 m1_sc = m1 / mTotal;
  REAL8 m2_sc = m2 / mTotal;
  /* Create output vectors */
  if (!((*S1) = XLALCreateREAL8Vector(3))) {
    XLALPrintError("XLAL Error: failed to allocate REAL8Vector S1.\n");
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if (!((*S2) = XLALCreateREAL8Vector(3))) {
    XLALDestroyREAL8Vector(*S1); // Free the memory above
    XLALPrintError("XLAL Error failed to allocate REAL8Vector S2.\n");
    XLAL_ERROR(XLAL_ENOMEM);
  }
  memset((*S1)->data, 0, 3 * sizeof(REAL8));
  memset((*S2)->data, 0, 3 * sizeof(REAL8));

  /* Local variables */
  REAL8 rvec[3] = {0, 0, 0};
  REAL8 drvec[3] = {0, 0, 0};
  REAL8 pvec[3] = {0, 0, 0};
  REAL8 spin1vec[3] = {0, 0, 0}; /* in units of mTotal^2 */
  REAL8 spin2vec[3] = {0, 0, 0}; /* in units of mTotal^2 */
  REAL8 crossp[3] = {0, 0, 0};
  REAL8 L_hat[3] = {0, 0, 0};
  REAL8 n_hat[3] = {0, 0, 0};
  REAL8 lambda_hat[3] = {0, 0, 0};
  /* Read from the extended dynamics values */
  for (int j = 0; j < 3; j++) {
    rvec[j] = seobvalues->data[1 + j];
    drvec[j] = seobvalues->data[15 + j];
    pvec[j] = seobvalues->data[4 + j];
    spin1vec[j] = seobvalues->data[7 + j] / m1_sc / m1_sc;
    spin2vec[j] = seobvalues->data[10 + j] / m2_sc / m2_sc;
  }
  if (flagZframe == FLAG_SEOBNRv4P_ZFRAME_L) {
    // Note: pvec is missing a factor of nu, but don't care about it's magnitide
    // anyway
    cross_product(rvec, pvec, crossp);
  } else if (flagZframe == FLAG_SEOBNRv4P_ZFRAME_LN) {
    cross_product(rvec, drvec, crossp);
  }
  REAL8 Lmag = sqrt(inner_product(crossp, crossp));
  REAL8 sep = sqrt(inner_product(rvec, rvec));

  for (int jj = 0; jj < 3; jj++) {
    L_hat[jj] = crossp[jj] / Lmag;
    n_hat[jj] = rvec[jj] / sep;
  }

  cross_product(L_hat, n_hat, lambda_hat);
  // Project onto the new frame
  (*S1)->data[0] = inner_product(spin1vec, n_hat);
  (*S2)->data[0] = inner_product(spin2vec, n_hat);

  (*S1)->data[1] = inner_product(spin1vec, lambda_hat);
  (*S2)->data[1] = inner_product(spin2vec, lambda_hat);

  (*S1)->data[2] = inner_product(spin1vec, L_hat);
  (*S2)->data[2] = inner_product(spin2vec, L_hat);
  return XLAL_SUCCESS;
}

/**
 * This function computes the J vector.
 */
static int SEOBJfromDynamics(
    REAL8Vector **J,         /**<< Output: pointer to vector J */
    REAL8Vector *seobvalues, /**<< Input: vector for extended dynamics values */
    SpinEOBParams *seobParams /**<< SEOB params */
) {
  UINT4 j;
  if ((!J) || (!seobvalues) || (!seobParams)) {
    XLALPrintError("Some pointers passed to SEOBJfromDynamics were null\n");
    XLAL_ERROR(XLAL_ENOMEM);
  }
  /* Create output vector */
  if (!((*J) = XLALCreateREAL8Vector(3))) {
    XLALPrintError("XLAL Error failed to allocate REAL8Vector J.\n");
    XLAL_ERROR(XLAL_ENOMEM);
  }
  memset((*J)->data, 0, 3 * sizeof(REAL8));

  /* Masses */
  REAL8 eta = seobParams->eobParams->eta;

  /* Local variables */
  REAL8 rvec[3] = {0, 0, 0};
  REAL8 pvec[3] = {0, 0, 0};
  REAL8 spin1vec[3] = {0, 0, 0}; /* in units of mTotal^2 */
  REAL8 spin2vec[3] = {0, 0, 0}; /* in units of mTotal^2 */
  REAL8 rcrossp[3] = {0, 0, 0};

  /* Read from the extended dynamics values */
  for (j = 0; j < 3; j++) {
    rvec[j] = seobvalues->data[1 + j];
    pvec[j] = seobvalues->data[4 + j];
    spin1vec[j] = seobvalues->data[7 + j];
    spin2vec[j] = seobvalues->data[10 + j];
  }
  cross_product(rvec, pvec, rcrossp);

  /* Compute J - restoring the factor eta in L (p stored in the dynamics is
   * p/mu) */
  for (j = 0; j < 3; j++) {
    (*J)->data[j] = eta * rcrossp[j] + spin1vec[j] + spin2vec[j];
  }

  return XLAL_SUCCESS;
}

/**
 * This function computes the L-hat vector.
 */
static int SEOBLhatfromDynamics(
    REAL8Vector **L,         /**<< Output: pointer to vector L */
    REAL8Vector *seobvalues, /**<< Input: vector for extended dynamics values */
    SpinEOBParams *seobParams /**<< SEOB params */,
    const flagSEOBNRv4P_Zframe
        flagZframe /**<<Input: whether to compute the L_N or L frame */
)
{
  if ((!L) || (!seobvalues) || (!seobParams))
  {
    XLALPrintError("Some pointers passed to SEOBLfromDynamics were null\n");
    XLAL_ERROR(XLAL_ENOMEM);
  }
  /* Create output vector */
  if (!((*L) = XLALCreateREAL8Vector(3)))
  {
    XLALPrintError("XLAL Error failed to allocate REAL8Vector L.\n");
    XLAL_ERROR(XLAL_ENOMEM);
  }
  memset((*L)->data, 0, 3 * sizeof(REAL8));

  /* Local variables */
  REAL8 rvec[3] = {0, 0, 0};
  REAL8 drvec[3] = {0, 0, 0};
  REAL8 pvec[3] = {0, 0, 0};
  REAL8 crossp[3] = {0, 0, 0};

  /* Read from the extended dynamics values */
  for (int j = 0; j < 3; j++)
  {
    rvec[j] = seobvalues->data[1 + j];
    drvec[j] = seobvalues->data[15 + j];
    pvec[j] = seobvalues->data[4 + j];
  }
  if (flagZframe == FLAG_SEOBNRv4P_ZFRAME_L)
  {
    // Note: pvec is missing a factor of nu, but don't care about it's magnitide
    // anyway
    cross_product(rvec, pvec, crossp);
  }
  else if (flagZframe == FLAG_SEOBNRv4P_ZFRAME_LN)
  {
    cross_product(rvec, drvec, crossp);
  }
  REAL8 Lmag = sqrt(inner_product(crossp, crossp));

  for (int jj = 0; jj < 3; jj++)
  {
    (*L)->data[jj] = crossp[jj] / Lmag;
  }

  return XLAL_SUCCESS;
}

/**
 * This function computes the Jframe unit vectors, with e3J along Jhat.
 * Convention: if (ex, ey, ez) is the initial I-frame, e1J chosen such that ex
 * is in the plane (e1J, e3J) and ex.e1J>0.
 * In the case where e3J and x happen to be close to aligned, we continuously
 * switch to another prescription with y playing the role of x
 */
static int SEOBBuildJframeVectors(
    REAL8Vector *e1J, /**<< Output: vector for e1J, already allocated */
    REAL8Vector *e2J, /**<< Output: vector for e2J, already allocated */
    REAL8Vector *e3J, /**<< Output: vector for e3J, already allocated */
    REAL8Vector *JVec /**<< Input: vector J */
) {
  UINT4 j;

  /* Checking size and of input vectors */
  if ((!e1J) || (!e2J) || (!e3J) || (!JVec)) {
    XLALPrintError("XLAL Error: at least one input pointer is NULL.\n");
    XLAL_ERROR(XLAL_EINVAL);
  }
  if ((!(e1J->length == 3)) || (!(e2J->length == 3)) || (!(e2J->length == 3))) {
    XLALPrintError(
        "XLAL Error: at least one input vector is not of length 3.\n");
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* Set e3J to Jhat */
  REAL8 Jnorm = sqrt(inner_product(JVec->data, JVec->data));
  for (j = 0; j < 3; j++) {
    e3J->data[j] = JVec->data[j] / Jnorm;
  }

  /* Set e1J to enforce the condition that x is in the plane (e1J, e3J) and
   * x.e1J>0 */
  /* Added a protection against the degenerate case where e3J, x are aligned:
   * let lambda = 1 - |ex.e3J|, measuring the alignment of e3J, x
   * for lambda < 1e-5 use y instead of x
   * for lambda > 1e-4 use x normally
   * for lambda in [1e-4, 1e-5] use a lambda-weighted combination of both
   * thresholds are arbitrary */
  REAL8 normfacx = 0.;
  REAL8 normfacy = 0.;
  REAL8 weightx = 0.;
  REAL8 weighty = 0.;
  REAL8 e1Jblendednorm = 0.;
  REAL8 exvec[3] = {1, 0, 0};
  REAL8 eyvec[3] = {0, 1, 0};
  REAL8 exdote3J = inner_product(exvec, e3J->data);
  REAL8 eydote3J = inner_product(eyvec, e3J->data);
  REAL8 lambda = 1. - fabs(exdote3J);
  if ((lambda < 0.) || (lambda > 1.)) {
    XLALPrintError("Problem: lambda=1-|e3J.ex|=%g, should be in [0,1]", lambda);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (lambda > 1e-4) {
    normfacx = 1. / sqrt(1. - exdote3J * exdote3J);
    for (j = 0; j < 3; j++) {
      e1J->data[j] = (exvec[j] - exdote3J * e3J->data[j]) / normfacx;
    }
  } else if (lambda < 1e-5) {
    normfacy = 1. / sqrt(1. - eydote3J * eydote3J);
    for (j = 0; j < 3; j++) {
      e1J->data[j] = (eyvec[j] - eydote3J * e3J->data[j]) / normfacy;
    }
  } else {
    weightx = (lambda - 1e-5) / (1e-4 - 1e-5);
    weighty = 1. - weightx;
    normfacx = 1. / sqrt(1. - exdote3J * exdote3J);
    normfacy = 1. / sqrt(1. - eydote3J * eydote3J);
    for (j = 0; j < 3; j++) {
      e1J->data[j] = weightx * (exvec[j] - exdote3J * e3J->data[j]) / normfacx +
                     weighty * (eyvec[j] - eydote3J * e3J->data[j]) / normfacy;
    }
    e1Jblendednorm = sqrt(inner_product(e1J->data, e1J->data));
    for (j = 0; j < 3; j++) {
      e1J->data[j] /= e1Jblendednorm;
    }
  }

  /* Get e2J = e3J * e1J */
  cross_product(e3J->data, e1J->data, e2J->data);

  /* Normally, vectors already of unit norm - we normalize again to eliminate
   * possible round-off error */
  REAL8 e1Jnorm = sqrt(inner_product(e1J->data, e1J->data));
  REAL8 e2Jnorm = sqrt(inner_product(e2J->data, e2J->data));
  REAL8 e3Jnorm = sqrt(inner_product(e3J->data, e3J->data));
  for (j = 0; j < 3; j++) {
    e1J->data[j] /= e1Jnorm;
    e2J->data[j] /= e2Jnorm;
    e3J->data[j] /= e3Jnorm;
  }

  return XLAL_SUCCESS;
}

/**
 * This function computes Euler angles I2J given the unit vectors of the Jframe.
 */
static int SEOBEulerI2JFromJframeVectors(
    REAL8 *alphaI2J,  /**<< Output: Euler angle alpha I2J */
    REAL8 *betaI2J,   /**<< Output: Euler angle beta I2J */
    REAL8 *gammaI2J,  /**<< Output: Euler angle gamma I2J */
    REAL8Vector *e1J, /**<< Input: unit Jframe vector e1J */
    REAL8Vector *e2J, /**<< Input: unit Jframe vector e2J */
    REAL8Vector *e3J  /**<< Input: unit Jframe vector e3J */
) {
  /* Active rotation matrix from frame (x,y,z) to frame (e1J,e2J,e3J) */
  /* The input vectors (eJ) are decomposed on the basis (x,y,z) */
  REAL8Array *R = XLALCreateREAL8ArrayL(2, 3, 3);
  if (!R) {
    XLALPrintError("Allocating the rotation matrix failed!");
    XLAL_ERROR(XLAL_ENOMEM);
  }
  RotationMatrixActiveFromBasisVectors(R, e1J->data, e2J->data, e3J->data);

  /* Compute Euler angles in the Z-Y-Z convention */
  EulerAnglesZYZFromRotationMatrixActive(alphaI2J, betaI2J, gammaI2J, R);

  /* Cleanup */
  XLALDestroyREAL8Array(R);

  return XLAL_SUCCESS;
}

/**
 * This function computes the NQC coefficients for a list of mode contributions.
 */
static int SEOBCalculateSphHarmListNQCCoefficientsV4(
    SphHarmListEOBNonQCCoeffs *
        *nqcCoeffsList, /**<< Output: non-quasi-circular coefficients as a list
                           for each mode */
    INT4 modes[][2],    /**<< Input: array of modes (l,m) */
    UINT4 nmodes,       /**<< Input: number of modes (l,m) */
    REAL8 tPeakOmega,   /**<< Input: time of peak of Omega */
    SEOBdynamics *seobdynamics,  /**<< Input: SEOB dynamics */
    SpinEOBParams *seobParams,   /**<< Input: SEOB params */
    REAL8Vector *chi1_omegaPeak, /**<< Input: dimensionless spin 1 at peak of
                                    omega in L_N frame */
    REAL8Vector *chi2_omegaPeak /**<< Input: dimensionless spin 2 at peak of
                                    omega in L_N frame */
) {
  /* Masses */
  REAL8 m1 = seobParams->eobParams->m1;
  REAL8 m2 = seobParams->eobParams->m2;
  seobParams->tPeakOmega = tPeakOmega;
  seobParams->spin1z_omegaPeak = chi1_omegaPeak->data[2];
  seobParams->spin2z_omegaPeak = chi2_omegaPeak->data[2];
  // REAL8 mtot = m1 + m2;
  // UINT4 SpinAlignedEOBversion =
  // seobParams->seobCoeffs->SpinAlignedEOBversion;

  /* Length of dynamics data and sampling step */
  UINT4 retLen = seobdynamics->length;
  REAL8 deltaT = seobdynamics->tVec[1] - seobdynamics->tVec[0];

  /* Create vectors from dynamics */
  REAL8Vector r, pr, omega;
  r.length = pr.length = omega.length = retLen;
  r.data = seobdynamics->polarrVec;
  pr.data = seobdynamics->polarprVec;
  omega.data = seobdynamics->omegaVec;

  /* Workspace vectors */
  /*
  REAL8Vector s1Vec, s2Vec, sigmaKerr;
  REAL8 s1Vecdata[3] = {0.};
  REAL8 s2Vecdata[3] = {0.};
  REAL8 sigmaKerrdata[3] = {0.};
  s1Vec.length = s2Vec.length = sigmaKerr.length = 3;
  s1Vec.data = s1Vecdata;
  s2Vec.data = s2Vecdata;
  sigmaKerr.data = sigmaKerrdata;
  */
  /* Final values for a and for the spins projected onto Z */
  /*
  REAL8 s1dotZfinal = seobdynamics->s1dotZVec[retLen-1];
  REAL8 s2dotZfinal = seobdynamics->s2dotZVec[retLen-1];
  REAL8 chi1dotZfinal = s1dotZfinal * mtot*mtot / (m1*m1);
  REAL8 chi2dotZfinal = s2dotZfinal * mtot*mtot / (m2*m2);
  REAL8 chiSfinal = SEOBCalculateChiS( chi1dotZfinal, chi2dotZfinal );
  REAL8 chiAfinal = SEOBCalculateChiA( chi1dotZfinal, chi2dotZfinal );
  s1Vec.data[0] = seobdynamics->s1Vecx[retLen-1];
  s1Vec.data[1] = seobdynamics->s1Vecy[retLen-1];
  s1Vec.data[2] = seobdynamics->s1Vecz[retLen-1];
  s2Vec.data[0] = seobdynamics->s2Vecx[retLen-1];
  s2Vec.data[1] = seobdynamics->s2Vecy[retLen-1];
  s2Vec.data[2] = seobdynamics->s2Vecz[retLen-1];
  SEOBCalculateSigmaKerr( &sigmaKerr, &s1Vec, &s2Vec );
  */
  // REAL8 afinal = sqrt( inner_product( sigmaKerr.data, sigmaKerr.data ) );

  REAL8 chi1dotZfinal = chi1_omegaPeak->data[2];
  REAL8 chi2dotZfinal = chi2_omegaPeak->data[2];
  REAL8 chiSfinal = SEOBCalculateChiS(chi1dotZfinal, chi2dotZfinal);
  REAL8 chiAfinal = SEOBCalculateChiA(chi1dotZfinal, chi2dotZfinal);
  REAL8 q = m1/m2;
  //printf("chiA = %.16f\n",chiAfinal);
  /* Time elapsed from the start of the dynamics to tPeakOmega */
  REAL8 tPeakOmegaFromStartDyn = tPeakOmega - seobdynamics->tVec[0];

  /* Compute NQC coefficients - output is nqcCoeffs */
  /* NOTE: internally, XLALSimIMRSpinEOBCalculateNQCCoefficientsV4 builds a time
   * vector as i*deltaT, starting at 0 - thus tPeakOmega has to be measured from
   * the start of the dynamics */

  /* Modes are to be generated without NQC */
  UINT4 includeNQC = 0;
  /* Loop over modes */
  for (UINT4 nmode = 0; nmode < nmodes; nmode++) {

    INT4 l = modes[nmode][0];
    INT4 m = modes[nmode][1];

    EOBNonQCCoeffs *nqcCoeffs = XLALMalloc(sizeof(EOBNonQCCoeffs));
    memset(nqcCoeffs, 0, sizeof(EOBNonQCCoeffs));
    /* In the equal mass equal spins case the odd-m modes are 0, so we set the NQCs to 0 */
    if (q<1.005 && (m % 2 != 0) && (fabs(chiAfinal) < 0.15)) { /* In this case, set NQC coeffs to 0 for odd m */
      nqcCoeffs->a1 = 0.;
      nqcCoeffs->a2 = 0.;
      nqcCoeffs->a3 = 0.;
      nqcCoeffs->a3S = 0.;
      nqcCoeffs->a4 = 0.;
      nqcCoeffs->a5 = 0.;
      nqcCoeffs->b1 = 0.;
      nqcCoeffs->b2 = 0.;
      nqcCoeffs->b3 = 0.;
      nqcCoeffs->b4 = 0.;
    } else {

      /* Compute amplitude and phase of the mode hlm pre-NQC */
      CAmpPhaseSequence *hlm = NULL;
      /* Mode hlm is to be generated without NQC */
      includeNQC = 0;
      SEOBCalculatehlmAmpPhase(&hlm, l, m, seobdynamics, nqcCoeffs, seobParams,
                               includeNQC);

      /* Cast to amp/phase vectors */
      REAL8Vector *hlm_amp = NULL;
      REAL8Vector *hlm_phase = NULL;
      if (!(hlm_amp = XLALCreateREAL8Vector(hlm->xdata->length))) {
        XLALPrintError(
            "XLAL Error - %s: failed to allocate REAL8Vector hlm_amp.\n",
            __func__);
        XLAL_ERROR(XLAL_ENOMEM);
      }
      if (!(hlm_phase = XLALCreateREAL8Vector(hlm->xdata->length))) {
        XLALPrintError(
            "XLAL Error - %s: failed to allocate REAL8Vector hlm_phase.\n",
            __func__);
        XLAL_ERROR(XLAL_ENOMEM);
      }
      memcpy(hlm_amp->data, hlm->camp_real->data,
             hlm->xdata->length * sizeof(REAL8));
      memcpy(hlm_phase->data, hlm->phase->data,
             hlm->xdata->length * sizeof(REAL8));

      /* Compute NQC */
      if (XLALSimIMRSpinEOBCalculateNQCCoefficientsV4(
              hlm_amp, hlm_phase, &r, &pr, &omega, l, m, tPeakOmegaFromStartDyn,
              deltaT, m1, m2, chiAfinal, chiSfinal,
              nqcCoeffs) == XLAL_FAILURE) {
        XLALPrintError(
            "XLAL Error - %s: failure for the mode (l,m) = (%d, %d).\n",
            __func__, l, m);
        XLAL_ERROR(XLAL_EFUNC);
      }

      /* Cleanup */
      CAmpPhaseSequence_Destroy(hlm);
      XLALDestroyREAL8Vector(hlm_amp);
      XLALDestroyREAL8Vector(hlm_phase);
    }

    /* Add computed NQCs to the list */
    SphHarmListEOBNonQCCoeffs_AddMode(nqcCoeffsList, nqcCoeffs, l, m);
  }

  return XLAL_SUCCESS;
}

/**
 * This function finds the index in a vector such that the value at the index
 *  closest to an input value.
 * Assumes the input vector is increasing (typically, times or frequencies of
 * series).
 *
 */

static UINT4 FindClosestIndex(
    REAL8Vector *vec, /**<< Input: monotonically increasing vector */
    REAL8 value       /**<< Input: value to look for */
) {
  REAL8 *data = vec->data;
  UINT4 length = vec->length;
  UINT4 index = 0;
  /* Get to the last data point lower than input */
  while ((index < length - 1) && (data[index + 1] <= value))
    index++;

  /* Check if this one or the next (first higher) is closest to input */
  if (index < length - 1) {
    if (fabs(data[index] - value) > fabs(data[index + 1] - value))
      index++;
  }
  return index;
}

/**
 * This function finds the value in a vector that is closest to an input value.
 * Assumes the input vector is increasing (typically, times or frequencies of
 * series). Purely for convenience.
 */
static REAL8 FindClosestValueInIncreasingVector(
    REAL8Vector *vec, /**<< Input: monotonically increasing vector */
    REAL8 value       /**<< Input: value to look for */
) {
  UINT4 index = FindClosestIndex(vec, value);
  return vec->data[index];
}

static int SEOBGetFinalSpinMass(
    REAL8
        *finalMass, /**<< Output: final mass computed from fit (scaled by M) */
    REAL8 *
        finalSpin, /**<< Output: final spin computed from fit (dimensionless) */
    REAL8Vector *seobvalues, /**<< Input: vector for dynamics values at time of
                                peak of omega */
    SpinEOBParams *seobParams, /**<< Input: SEOB params */
    const flagSEOBNRv4P_Zframe
        flagZframe /**<< Input: Whether to use the L_N or L frame */
) {
  // int debug = 0;
  /* Masses */
  REAL8 m1 = seobParams->eobParams->m1;
  REAL8 m2 = seobParams->eobParams->m2;

  Approximant seobApproximant = seobParams->seobApproximant;

  /* Compute components of vectors chi1, chi2 from seobvalues at tPeakOmega in
   * the L-frame */

  REAL8Vector *chi1temp = NULL;
  REAL8Vector *chi2temp = NULL;
  if (SEOBLFrameVectors(&chi1temp, &chi2temp, seobvalues, m1, m2, flagZframe) ==
      XLAL_FAILURE) {
    XLALPrintError("XLAL Error - %s: failure in SEOBLFrameVectors.\n",
                   __func__);
    XLAL_ERROR(XLAL_EFUNC);
  }

  REAL8 chi1L[3] = {chi1temp->data[0], chi1temp->data[1], chi1temp->data[2]};
  REAL8 chi2L[3] = {chi2temp->data[0], chi2temp->data[1], chi2temp->data[2]};
  if (XLALSimIMREOBFinalMassSpinPrec(finalMass, finalSpin, m1, m2, chi1L, chi2L,
                                     seobApproximant) == XLAL_FAILURE) {
    XLALPrintError(
        "XLAL Error - %s: failure in XLALSimIMREOBFinalMassSpinPrec.\n",
        __func__);
    XLAL_ERROR(XLAL_EFUNC);
  }

  XLALDestroyREAL8Vector(chi1temp);
  XLALDestroyREAL8Vector(chi2temp);

  return XLAL_SUCCESS;
}

/**
 * This function attaches the ringdown to the P-frame modes hlm.
 * Here seobvalues have been interpolated from seobdynamics at tPeakOmega (can
 * differ from tAttach).
 */
static int SEOBAttachRDToSphHarmListhPlm(
    SphHarmListCAmpPhaseSequence *
        *listhPlm_RDattached, /**<< Output: list of extended modes hlm with RD
                                 attached */
    COMPLEX16Vector *
        *sigmaQNMlm0, /**<< Output: list of QNM complex frequency for modes lm,
                         0th overtone (dimensionless) */
    INT4 modes[][2],  /**<< Input: array of modes (l,m) */
    UINT4 nmodes,     /**<< Input: number of modes (l,m) */
    REAL8 finalMass,  /**<< Input: final mass computed from fit (scaled by M) */
    REAL8
        finalSpin, /**<< Input: final spin computed from fit (dimensionless) */
    SphHarmListCAmpPhaseSequence *listhPlm, /**<< Input: list of modes hlm */
    REAL8 deltaT,                           /**<< Input: time step */
    UINT4 retLen,        /**<< Input: length of the input modes and dynamics */
    UINT4 retLenRDPatch, /**<< Input: length of the ringdown patch */
    REAL8 tAttach,       /**<< Input: time of attachment */
    // REAL8 tStart, /**<< Input: starting time (of the HiS part) */
    REAL8Vector *seobvalues, /**<< Input: vector for dynamics values at time of
                                peak of omega */
    SEOBdynamics *seobdynamics,            /**<< Input: SEOB dynamics */
    SpinEOBParams *seobParams,             /**<< SEOB params */
    const flagSEOBNRv4P_Zframe flagZframe, /**<< Input: whether to use L_N or L
                                             frame for spin projections */
    const INT4 debug /**<< Input: flag to print debug information */
) {

  /* Check that the input list pointer and vector pointer are NULL */
  if (!(*listhPlm_RDattached == NULL)) {
    XLALPrintError("XLAL Error - %s: output pointer for the list "
                   "hPlm_RDattached is not NULL.\n",
                   __func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if (!(*sigmaQNMlm0 == NULL)) {
    XLALPrintError("XLAL Error - %s: output pointer for the vector sigmaQNMlm0 "
                   "is not NULL.\n",
                   __func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }

  /* Create list of output modes */
  UINT4 retLen_RDattached = retLen + retLenRDPatch;
  for (UINT4 nmode = 0; nmode < nmodes; nmode++) {

    INT4 l = modes[nmode][0];
    INT4 m = modes[nmode][1];

    CAmpPhaseSequence *hPlm_RDattached = NULL;
    if (CAmpPhaseSequence_Init(&hPlm_RDattached, retLen_RDattached) ==
        XLAL_FAILURE) {
      XLALPrintError("XLAL Error - %s: failed to allocate CAmpPhaseSequence "
                     "hlm for mode (l,m) = (%d,%d).\n",
                     __func__, l, m);
      XLAL_ERROR(XLAL_ENOMEM);
    }

    SphHarmListCAmpPhaseSequence_AddMode(listhPlm_RDattached, hPlm_RDattached,
                                         l, m);
  }

  /* Masses */
  REAL8 m1 = seobParams->eobParams->m1;
  REAL8 m2 = seobParams->eobParams->m2;
  // REAL8 eta = seobParams->eobParams->eta;
  REAL8 mTotal = m1 + m2;
  REAL8 mTScaled = mTotal * LAL_MTSUN_SI;
  Approximant seobApproximant = seobParams->seobApproximant;

  /* Vector of times without attachment */
  REAL8Vector timeVec;
  timeVec.length = retLen;
  timeVec.data = seobdynamics->tVec;
  REAL8 finalM, finalS = 0;

  REAL8Vector *chi1temp = NULL;
  REAL8Vector *chi2temp = NULL;
  SEOBLFrameVectors(&chi1temp, &chi2temp, seobvalues, m1, m2, flagZframe);

  REAL8 chi1Lx = chi1temp->data[0];
  REAL8 chi1Ly = chi1temp->data[1];
  REAL8 chi1Lz = chi1temp->data[2];
  REAL8 chi2Lx = chi2temp->data[0];
  REAL8 chi2Ly = chi2temp->data[1];
  REAL8 chi2Lz = chi2temp->data[2];

  /* finalSpin interpolation is available only between -0.9996 and 0.9996 */
  /* Set finalSpin to +/- 0.9996 if it is out of this range */
  finalS = finalSpin;
  finalM = finalMass;

  if (finalS < -0.9996)
    finalS = -0.9996;
  if (finalS > 0.9996)
    finalS = 0.9996;

  if (debug) {
    XLAL_PRINT_INFO("In RD attachment: final mass = %e, final spin = %e, "
                    "total_mass = %e \n",
                    finalM, finalS, mTotal);
  }

  /* Compute QNM frequencies */
  /* Generate 1 overtone (the 0th overtone) */
  // NOTE: this is redone internally in XLALSimIMREOBAttachFitRingdown

  *sigmaQNMlm0 = XLALCreateCOMPLEX16Vector(nmodes);
  COMPLEX16Vector sigmaQNMlm0physicalVec;
  sigmaQNMlm0physicalVec.length = 1;
  COMPLEX16 sigmaQNMlm0physicalval = 0.;
  sigmaQNMlm0physicalVec.data = &sigmaQNMlm0physicalval;
  for (UINT4 nmode = 0; nmode < nmodes; nmode++) {

    INT4 l = modes[nmode][0];
    INT4 m = modes[nmode][1];

    /* NOTE: XLALSimIMREOBGenerateQNMFreqV2Prec returns the complex frequency in
     * physical units... */
    if (XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec(&sigmaQNMlm0physicalVec, m1,
                                                    m2, finalM, finalS, l, m,
                                                    1) == XLAL_FAILURE) {
      XLALPrintError("XLAL Error - %s: failure in "
                     "XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec for mode "
                     "(l,m) = (%d,%d).\n",
                     __func__, l, m);
      XLAL_ERROR(XLAL_EFUNC);
    }

    (*sigmaQNMlm0)->data[nmode] = mTScaled * sigmaQNMlm0physicalVec.data[0];

    if (debug) {
      XLAL_PRINT_INFO("Complex QNM frequency: (%d,%d,0) = %.16e + I*%.16e\n", l,
                      m, creal((*sigmaQNMlm0)->data[nmode]),
                      cimag((*sigmaQNMlm0)->data[nmode]));
    }
  }

  /* Find the time sample that is closest to tAttach */
  REAL8 timeNeartAttach = FindClosestValueInIncreasingVector(&timeVec, tAttach);
  // R.C: The attachment point is different for the 55 mode wrt the other modes
  // they are related by tAttach55 = tAttach -10M
  REAL8 timeNeartAttach55 =
      FindClosestValueInIncreasingVector(&timeVec, tAttach - 10.);
  /* This structure is inherited from the time when the attachment was done with
   * a comb */
  /* Only the value data[1] will be used -- the other two are ignored */
  REAL8Vector rdMatchPoint;
  REAL8 rdMatchPointdata[4] = {0.};
  rdMatchPoint.length = 3;
  rdMatchPoint.data = rdMatchPointdata;
  rdMatchPoint.data[0] = 0.; /* unused */
  rdMatchPoint.data[1] =
      timeNeartAttach;       /* this is the time closest to tAttach */
  rdMatchPoint.data[2] = 0.; /* unused */
  rdMatchPoint.data[3] = timeNeartAttach55; /* tAttach55 = tAttach -10M */

  /* Create real and imaginary part vectors for the mode with RD attached */
  REAL8Vector *hPlmRe = NULL;
  REAL8Vector *hPlmIm = NULL;
  if (!(hPlmRe = XLALCreateREAL8Vector(retLen_RDattached))) {
    XLALPrintError("XLAL Error - %s: failed to allocate REAL8Vector hPlmRe.\n",
                   __func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if (!(hPlmIm = XLALCreateREAL8Vector(retLen_RDattached))) {
    XLALPrintError("XLAL Error - %s: failed to allocate REAL8Vector hPlmIm.\n",
                   __func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  memset(hPlmRe->data, 0, (hPlmRe->length) * sizeof(REAL8));
  memset(hPlmIm->data, 0, (hPlmIm->length) * sizeof(REAL8));

  /* This is used to keep track of the location of the max amplitude of the 22
   * mode */
  UINT4 indAmpMax = 0;

  /* Loop over modes */
  for (UINT4 nmode = 0; nmode < nmodes; nmode++) {

    INT4 l = modes[nmode][0];
    INT4 m = modes[nmode][1];

    /* Get the relevant modes in lists */
    CAmpPhaseSequence *hPlm =
        SphHarmListCAmpPhaseSequence_GetMode(listhPlm, l, m)->campphase;
    CAmpPhaseSequence *hPlm_RDattached =
        SphHarmListCAmpPhaseSequence_GetMode(*listhPlm_RDattached, l, m)
            ->campphase;

    COMPLEX16 hPlm_val = 0.;
    for (UINT4 i = 0; i < retLen; i++) {
      hPlm_val = (hPlm->camp_real->data[i] + I * hPlm->camp_imag->data[i]) *
                 cexp(I * hPlm->phase->data[i]);
      hPlmRe->data[i] = creal(hPlm_val);
      hPlmIm->data[i] = cimag(hPlm_val);
    }

    /* NOTE: deltaT here is in physical units (s) */
    REAL8 deltaTseconds = deltaT * mTScaled;
    if (XLALSimIMREOBAttachFitRingdown(
            hPlmRe, hPlmIm, l, m, deltaTseconds, m1, m2, chi1Lx, chi1Ly, chi1Lz,
            chi2Lx, chi2Ly, chi2Lz, finalM, finalS, &timeVec, &rdMatchPoint,
            seobApproximant, &indAmpMax) == XLAL_FAILURE) {
      XLALPrintError("XLAL Error - %s: XLALSimIMREOBAttachFitRingdown failed "
                     "for mode (l,m) = (%d,%d).\n",
                     __func__, l, m);
      XLAL_ERROR(XLAL_EFUNC);
    }

    /* Copy times in output, using deltaT */
    for (UINT4 i = 0; i < retLen; i++) {
      hPlm_RDattached->xdata->data[i] = hPlm->xdata->data[i];
    }
    for (UINT4 i = 0; i < retLen_RDattached - retLen; i++) {
      hPlm_RDattached->xdata->data[retLen + i] =
          hPlm_RDattached->xdata->data[retLen - 1] + i * deltaT;
    }

    /* Translate result in amp/phase, unwrap phase in place */
    for (UINT4 i = 0; i < retLen_RDattached; i++) {
      hPlm_val = hPlmRe->data[i] + I * hPlmIm->data[i];
      hPlm_RDattached->camp_real->data[i] = cabs(hPlm_val);
      hPlm_RDattached->camp_imag->data[i] = 0.;
      hPlm_RDattached->phase->data[i] = carg(hPlm_val);
    }
    XLALREAL8VectorUnwrapAngle(hPlm_RDattached->phase, hPlm_RDattached->phase);
  }

  /* Cleanup */
  XLALDestroyREAL8Vector(hPlmRe);
  XLALDestroyREAL8Vector(hPlmIm);
  XLALDestroyREAL8Vector(chi1temp);
  XLALDestroyREAL8Vector(chi2temp);
  return XLAL_SUCCESS;
}

/**
 * This function constructs the joined vector of times (AdaS+HiS+RDpatch) and
 * keeps the jonction indices and times RDpatch is the extension of HiS for the
 * RD, keeping the same constant sampling (indexJoinHiS, tJoinHiS) is given by
 * the first time sample >=tstartHiS (indexJoinAttach, tJoinAttach) is given by
 * the first time sample >=tAttach
 */
static int SEOBJoinTimeVector(
    REAL8Vector **tVecPmodes, /**<< Output: vector of times for P-modes
                                 (AdaS+HiS+RDpatch) */
    UINT4 *retLenPmodes,      /**<< Output: length of output vector of times for
                                 P-modes */
    REAL8 *tJoinHiS,          /**<< Output: first time >= tstartHiS */
    UINT4 *indexJoinHiS,      /**<< Output: first index >= tstartHiS */
    REAL8 *tJoinAttach,       /**<< Output: first time >= tAttach */
    UINT4 *indexJoinAttach,   /**<< Output: first index >= tAttach */
    UINT4 retLenHiSRDpatch,   /**<< Input: length of RD patch to be added at the
                                 end of HiS with the same constant sampling */
    REAL8 deltaTHiS,          /**<< Input: time step for the high sampling */
    REAL8 tstartHiS,          /**<< Input: time of start of HiS */
    REAL8 tAttach,            /**<< Input: time of attachment */
    SEOBdynamics
        *seobdynamicsAdaS, /**<< Input: SEOB dynamics with adaptive-sampling */
    SEOBdynamics
        *seobdynamicsHiS /**<< Input: SEOB dynamics with high-sampling */
) {
  /* Read from inputs */
  INT4 lenAdaS = seobdynamicsAdaS->length;
  REAL8 *tAdaS = seobdynamicsAdaS->tVec;
  INT4 lenHiS = seobdynamicsHiS->length;
  REAL8 *tHiS = seobdynamicsHiS->tVec;

  /* Determine how many time samples of AdaS to include */
  INT4 iJoinHiS = 0;
  while ((iJoinHiS < lenAdaS - 1) && (tAdaS[iJoinHiS] < tstartHiS))
    iJoinHiS++;

  /* Total length of output and create output vector */
  INT4 lenPmodes = iJoinHiS + lenHiS + retLenHiSRDpatch;
  if (!((*tVecPmodes) = XLALCreateREAL8Vector(lenPmodes))) {
    XLALPrintError(
        "XLAL Error - %s: failed to allocate REAL8Vector tVecPmodes.\n",
        __func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  *retLenPmodes = lenPmodes;

  /* Copy time values for AdaS and HiS */
  memcpy(&((*tVecPmodes)->data[0]), tAdaS, iJoinHiS * sizeof(REAL8));
  memcpy(&((*tVecPmodes)->data[iJoinHiS]), tHiS, lenHiS * sizeof(REAL8));

  /* Set time values for RD patch */
  for (UINT4 i = 0; i < retLenHiSRDpatch; i++) {
    (*tVecPmodes)->data[iJoinHiS + lenHiS + i] =
        tHiS[lenHiS - 1] + (i + 1) * deltaTHiS;
  }

  /* Determine how many time samples after tAttach */
  INT4 iJoinAttach = lenPmodes - 1;
  while ((iJoinAttach > 0) && ((*tVecPmodes)->data[iJoinAttach] > tAttach))
    iJoinAttach--;

  /* Output joining indices and times */
  *indexJoinHiS = iJoinHiS;
  *tJoinHiS = (*tVecPmodes)->data[iJoinHiS];
  *indexJoinAttach = iJoinAttach;
  *tJoinAttach = (*tVecPmodes)->data[iJoinAttach];

  return XLAL_SUCCESS;
}

/**
 * This function copies dynamics from AdaS<HiS and HiS<tAttach to form joined
 * dynamics, ending at the last time sample <tAttach
 */
// NOTE: we cut the dynamics at tAttach, as we will extend the Euler
// angles for t>=tAttach -- but we could also choose to finish at tPeakOmega
// which is used for final-J and for the final mass/spin fit
static int SEOBJoinDynamics(
    SEOBdynamics *
        *seobdynamicsJoined,     /**<< Output: pointer to joined dynamics */
    SEOBdynamics *seobdynamics1, /**<< Input: first dynamics */
    SEOBdynamics *seobdynamics2, /**<< Input: second dynamics */
    UINT4 indexJoin12, /**<< Input: index where to join the two dynamics */
    UINT4 indexEnd2    /**<< Input: index of the joined dynamics where to stop
                          dynamics 2 (excluded) */
) {
  UINT4 j;

  /* Lengths */
  INT4 lenJoined = indexEnd2;
  INT4 lendyn1Joined = indexJoin12;
  INT4 lendyn2Joined = indexEnd2 - indexJoin12;
  INT4 lendyn1 = seobdynamics1->length;
  INT4 lendyn2 = seobdynamics2->length;

  /* Initialize output dynamics */
  SEOBdynamics_Init(seobdynamicsJoined, lenJoined);

  /* Copy truncated dynamics 1 - v4PdynamicsVariables data fields */
  for (j = 0; j < v4PdynamicsVariables; j++) {
    memcpy(&((*seobdynamicsJoined)->array->data[j * lenJoined]),
           &(seobdynamics1->array->data[j * lendyn1]),
           lendyn1Joined * sizeof(REAL8));
  }

  /* Copy truncated dynamics 2 - v4PdynamicsVariables data fields */
  for (j = 0; j < v4PdynamicsVariables; j++) {
    memcpy(&((*seobdynamicsJoined)->array->data[j * lenJoined + lendyn1Joined]),
           &(seobdynamics2->array->data[j * lendyn2]),
           lendyn2Joined * sizeof(REAL8));
  }

  return XLAL_SUCCESS;
}

/**
 * This function copies dynamics from AdaS<HiS and HiS<tAttach to form joined
 * dynamics, ending at the last time sample <tAttach
 */
// NOTE: we cut the dynamics at tAttach, as we will extend the Euler
// angles for t>=tAttach -- but we could also choose to finish at tPeakOmega
// which is used for final-J and for the final mass/spin fit
static int SEOBJoinSphHarmListhlm(
    SphHarmListCAmpPhaseSequence *
        *listhlm_joined, /**<< Output: list of joined modes */
    SphHarmListCAmpPhaseSequence *listhlm_1, /**<< Input: list of modes 1 */
    SphHarmListCAmpPhaseSequence *listhlm_2, /**<< Input: list of modes 2 */
    INT4 modes[][2],  /**<< Input: array of modes (l,m) */
    UINT4 nmodes,     /**<< Input: number of modes (l,m) */
    UINT4 indexJoin12 /**<< Input: index where to join the two dynamics */
) {
  /* Loop over modes */
  for (UINT4 nmode = 0; nmode < nmodes; nmode++) {

    INT4 l = modes[nmode][0];
    INT4 m = modes[nmode][1];

    /* Get the relevant modes in lists */
    CAmpPhaseSequence *hlm_1 =
        SphHarmListCAmpPhaseSequence_GetMode(listhlm_1, l, m)->campphase;
    CAmpPhaseSequence *hlm_2 =
        SphHarmListCAmpPhaseSequence_GetMode(listhlm_2, l, m)->campphase;

    /* Lengths */
    UINT4 len2 = hlm_2->xdata->length;
    INT4 lenJoined1 = indexJoin12;
    INT4 lenJoined = lenJoined1 + len2;

    /* Real and imaginary part vectors for the mode with RD attached */
    CAmpPhaseSequence *hlm_joined = NULL;
    if (CAmpPhaseSequence_Init(&hlm_joined, lenJoined) == XLAL_FAILURE) {
      XLALPrintError("XLAL Error - %s: failed to allocate CAmpPhaseSequence "
                     "hlm_joined for mode (l,m) = (%d,%d).\n",
                     __func__, l, m);
      XLAL_ERROR(XLAL_ENOMEM);
    }

    /* Copy data, stopping part 1 at the specified index */
    memcpy(&(hlm_joined->xdata->data[0]), hlm_1->xdata->data,
           lenJoined1 * sizeof(REAL8));
    memcpy(&(hlm_joined->camp_real->data[0]), hlm_1->camp_real->data,
           lenJoined1 * sizeof(REAL8));
    memcpy(&(hlm_joined->camp_imag->data[0]), hlm_1->camp_imag->data,
           lenJoined1 * sizeof(REAL8));
    memcpy(&(hlm_joined->phase->data[0]), hlm_1->phase->data,
           lenJoined1 * sizeof(REAL8));
    /* Copy data, for part 2 starting from the specified index */
    memcpy(&(hlm_joined->xdata->data[lenJoined1]), hlm_2->xdata->data,
           len2 * sizeof(REAL8));
    memcpy(&(hlm_joined->camp_real->data[lenJoined1]), hlm_2->camp_real->data,
           len2 * sizeof(REAL8));
    memcpy(&(hlm_joined->camp_imag->data[lenJoined1]), hlm_2->camp_imag->data,
           len2 * sizeof(REAL8));
    memcpy(&(hlm_joined->phase->data[lenJoined1]), hlm_2->phase->data,
           len2 * sizeof(REAL8));
    /* Adjust for a 2kpi-shift in the phase - use the difference between the
     * last included sample of 1 and the first sample of 2 */
    REAL8 phase_diff =
        hlm_2->phase->data[0] - hlm_1->phase->data[lenJoined1 - 1];
    REAL8 shift_2kpi =
        -floor((phase_diff + LAL_PI) / (2 * LAL_PI)) * (2 * LAL_PI);
    for (UINT4 i = 0; i < len2; i++) {
      hlm_joined->phase->data[lenJoined1 + i] += shift_2kpi;
    }

    SphHarmListCAmpPhaseSequence_AddMode(listhlm_joined, hlm_joined, l, m);
  }

  return XLAL_SUCCESS;
}

/**
 * This function finds the (first) index and time with the largest
 * sum-of-squares amplitude - discrete, no interpolation The sum-of-squares
 * amplitude is the sum of the amplitude square of all modes of the waveform,
 * rotationally invariant for any l. This function is l=2 only, requires 22
 * to be present and uses 21 if present.
 */
static int SEOBAmplitudePeakFromAmp22Amp21(
    REAL8 *tPeak,                           /**<< Output: time of peak */
    UINT4 *indexPeak,                       /**<< Output: index of peak */
    SphHarmListCAmpPhaseSequence *listhPlm, /**<< Input: list of modes hlm */
    INT4 modes[][2],                        /**<< Input: array of modes (l,m) */
    UINT4 nmodes,     /**<< Input: number of modes (l,m) */
    REAL8Vector *tVec /**<< Input: vector of times */
) {
  UINT4 length = tVec->length;

  /* Looking form modes 22 and 21 */
  UINT4 found22 = 0, found21 = 0;
  for (UINT4 nmode = 0; nmode < nmodes; nmode++) {
    if (modes[nmode][0] == 2 && modes[nmode][1] == 2)
      found22 = 1;
    if (modes[nmode][0] == 2 && modes[nmode][1] == 1)
      found21 = 1;
  }
  if ((!found22)) {
    XLALPrintError("XLAL Error - %s: mode 22 not found.\n", __func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

  CAmpPhaseSequence *hP22 = NULL;
  CAmpPhaseSequence *hP21 = NULL;

  hP22 = SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 2, 2)->campphase;
  if (found21)
    hP21 = SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 2, 1)->campphase;

  /* Check lengths */
  if ((!(hP22->xdata->length == length)) ||
      (found21 && !(hP21->xdata->length == length))) {
    XLALPrintError("XLAL Error - %s: lengths of input amplitude and time "
                   "REAL8Vector do not match.\n",
                   __func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* Look for max combined amplitude - done discretely, no interpolation */
  UINT4 indexMax = 0;
  REAL8 AsquareMax = 0.;
  REAL8 A22_real = 0., A22_imag = 0., A21_real = 0., A21_imag = 0.,
        Asquare = 0.;
  for (UINT4 i = 0; i < length; i++) {
    A22_real = hP22->camp_real->data[i];
    A22_imag = hP22->camp_imag->data[i];
    Asquare = A22_real * A22_real + A22_imag * A22_imag;
    if (found21) {
      A21_real = hP21->camp_real->data[i];
      A21_imag = hP21->camp_imag->data[i];
      Asquare += A21_real * A21_real + A21_imag * A21_imag;
    }
    if (Asquare > AsquareMax) {
      AsquareMax = Asquare;
      indexMax = i;
    }
  }

  /* Output */
  *indexPeak = indexMax;
  *tPeak = tVec->data[indexMax];

  return XLAL_SUCCESS;
}

/**
 * This function computes the Euler angles from J-frame to P-frame given the
 * dynamics and basis vectors of the J-frame Angle gamma computed according to
 * minimal rotation condition with gamma=-alpha initially Note that all
 * quantities in the dynamics and the basis vectors eJ are expressed in the
 * initial I-frame
 */
static int SEOBEulerJ2PFromDynamics(
    REAL8Vector **alphaJ2P, /**<< Output: pointer to vector for alpha J2P */
    REAL8Vector **betaJ2P,  /**<< Output: pointer to vector for beta J2P */
    REAL8Vector **gammaJ2P, /**<< Output: pointer to vector for gamma J2P */
    REAL8Vector *e1J,       /**<< Input: unit Jframe vector e1J */
    REAL8Vector *e2J,       /**<< Input: unit Jframe vector e2J */
    REAL8Vector *e3J,       /**<< Input: unit Jframe vector e3J */
    UINT4 retLen, /**<< Input: total length of Euler angles data to be allocated
                     (length of P-modes) */
    UINT4 indexStop, /**<< Input: index where we stop the computation (excluded,
                        index of time of attachment) */
    SEOBdynamics *seobdynamics, /**<<Input: SEOB dynamics (joined AdaS+HiS, up
                                   to tAttach) */
    SpinEOBParams *seobParams,  /**<< SEOB params */
    flagSEOBNRv4P_Zframe
        flagZframe /**<< flag to choose Z direction of the frame, LN or L */
) {
  UINT4 i, j;
  UINT4 dynlength = seobdynamics->length;
  UINT4 SpinsAlmostAligned = seobParams->alignedSpins;

  /* Length of the subset where we want to compute the Euler angles from the
   * dynamics */
  UINT4 retLenDyn = indexStop;
  /* Check lengths -- note that indexStop is excluded */
  if (!((retLenDyn <= dynlength) && (dynlength <= retLen))) {
    XLALPrintError("XLAL Error - %s: incompatible lengths.\n", __func__);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* Create output vectors */
  if (!((*alphaJ2P) = XLALCreateREAL8Vector(retLen))) {
    XLALPrintError(
        "XLAL Error - %s: failed to allocate REAL8Vector alphaJ2P.\n",
        __func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if (!((*betaJ2P) = XLALCreateREAL8Vector(retLen))) {
    XLALPrintError("XLAL Error - %s: failed to allocate REAL8Vector betaJ2P.\n",
                   __func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if (!((*gammaJ2P) = XLALCreateREAL8Vector(retLen))) {
    XLALPrintError(
        "XLAL Error - %s: failed to allocate REAL8Vector gammaJ2P.\n",
        __func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  memset((*alphaJ2P)->data, 0, retLen * sizeof(REAL8));
  memset((*betaJ2P)->data, 0, retLen * sizeof(REAL8));
  memset((*gammaJ2P)->data, 0, retLen * sizeof(REAL8));

  if (!SpinsAlmostAligned) /* if spins are almost aligned, leave all Euler
                              angles to 0 */
  {
    /* Time vector of SEOB dynamics */
    REAL8 *tVec = seobdynamics->tVec;

    /* Local variables */
    REAL8 rvec[3] = {0, 0, 0};
    REAL8 pvec[3] = {0, 0, 0};
    REAL8 rdotvec[3] = {0, 0, 0};
    REAL8 Lhat[3] = {0, 0, 0};
    REAL8 LNhat[3] = {0, 0, 0};
    REAL8 Zframe[3] = {0, 0, 0};
    REAL8 e1PiniIbasis[3] = {0, 0, 0};
    REAL8 e2PiniIbasis[3] = {0, 0, 0};
    REAL8 e3PiniIbasis[3] = {0, 0, 0};
    REAL8 e1PiniJbasis[3] = {0, 0, 0};
    REAL8 e2PiniJbasis[3] = {0, 0, 0};

    /* Loop over time samples to compute alpha and beta, stopping at retLenDyn
     */
    for (i = 0; i < retLenDyn; i++) {

      /* Read from the extended dynamics values */
      for (j = 0; j < 3; j++) {
        rvec[j] = seobdynamics->array->data[(1 + j) * dynlength + i];
        pvec[j] = seobdynamics->array->data[(4 + j) * dynlength + i];
        rdotvec[j] = seobdynamics->array->data[(15 + j) * dynlength + i];
      }

      /* Compute Z-axis of the precessing frame, L or LN */
      if (flagZframe == FLAG_SEOBNRv4P_ZFRAME_L) {
        /* Compute Lhat */
        cross_product(rvec, pvec, Lhat);
        REAL8 Lhatnorm = sqrt(inner_product(Lhat, Lhat));
        for (j = 0; j < 3; j++) {
          Lhat[j] /= Lhatnorm;
        }
        memcpy(Zframe, Lhat, 3 * sizeof(REAL8));
      } else if (flagZframe == FLAG_SEOBNRv4P_ZFRAME_LN) {
        /* Compute LNhat */
        cross_product(rvec, rdotvec, LNhat);
        REAL8 LNhatnorm = sqrt(inner_product(LNhat, LNhat));
        for (j = 0; j < 3; j++) {
          LNhat[j] /= LNhatnorm;
        }
        memcpy(Zframe, LNhat, 3 * sizeof(REAL8));
      } else {
        XLALPrintError("XLAL Error - %s: flagZframe not recognized.\n",
                       __func__);
        XLAL_ERROR(XLAL_EINVAL);
      }

      /* Compute Z projected in the J-frame (e1J,e2J,e3J) */
      REAL8 Ze1J = inner_product(Zframe, e1J->data);
      REAL8 Ze2J = inner_product(Zframe, e2J->data);
      REAL8 Ze3J = inner_product(Zframe, e3J->data);

      /* Get Euler angles alpha (to be unwrapped later) and beta */
      (*alphaJ2P)->data[i] = atan2(Ze2J, Ze1J);
      (*betaJ2P)->data[i] = acos(Ze3J);

      /* At initial time, compute the initial vectors (e1P, e2P, e3P) decomposed
       * in the frame (e1J, e2J, e3J) */
      /* This will be needed to set initial gamma angle */
      if (i == 0) {
        /* e3P is the Z axis of the precessing frame */
        memcpy(e3PiniIbasis, Zframe, 3 * sizeof(REAL8));
        /* e1P is the unit separation vector n */
        memcpy(e1PiniIbasis, rvec, 3 * sizeof(REAL8));
        REAL8 e1PiniIbasisnorm =
            sqrt(inner_product(e1PiniIbasis, e1PiniIbasis));
        for (j = 0; j < 3; j++) {
          e1PiniIbasis[j] /= e1PiniIbasisnorm;
        }
        /* e2P is obtained by completing the triad */
        cross_product(e3PiniIbasis, e1PiniIbasis, e2PiniIbasis);
        /* Components of vectors eP in the frame eJ */
        e1PiniJbasis[0] = inner_product(e1PiniIbasis, e1J->data);
        e1PiniJbasis[1] = inner_product(e1PiniIbasis, e2J->data);
        e1PiniJbasis[2] = inner_product(e1PiniIbasis, e3J->data);
        e2PiniJbasis[0] = inner_product(e2PiniIbasis, e1J->data);
        e2PiniJbasis[1] = inner_product(e2PiniIbasis, e2J->data);
        e2PiniJbasis[2] = inner_product(e2PiniIbasis, e3J->data);
      }
    }

    /* Unwrap alpha in-place */
    XLALREAL8VectorUnwrapAngle((*alphaJ2P), (*alphaJ2P));

    /* Compute gamma according to minimal rotation condition */
    /* gamma is set initially so that the initial P-frame reproduces the initial
     * (n, lambda, Lhat or LNhat) frame */
    REAL8 InitialGamma = atan2(e2PiniJbasis[2], -e1PiniJbasis[2]);

    // Integrate \dot{\alpha} \cos{\beta} to get the final Euler angle
    // Eq. 20 of PRD 89, 084006 (2014) [arXiv:1307.6232]

    // Here 1000 referes to the number of subintervals that can be used when
    // performing adaptive quadrature to compute the integral.
    // See
    // https://www.gnu.org/software/gsl/manual/html_node/QAG-adaptive-integration.html
    REAL8 precEulerresult = 0, precEulererror = 0;
    gsl_integration_workspace *precEulerw =
        gsl_integration_workspace_alloc(1000);
    gsl_function precEulerF;
    PrecEulerAnglesIntegration precEulerparams;
    gsl_spline *x_spline = gsl_spline_alloc(gsl_interp_cspline, retLenDyn);
    gsl_spline *y_spline = gsl_spline_alloc(gsl_interp_cspline, retLenDyn);
    gsl_interp_accel *x_acc = gsl_interp_accel_alloc();
    gsl_interp_accel *y_acc = gsl_interp_accel_alloc();
    gsl_spline_init(x_spline, seobdynamics->tVec, (*alphaJ2P)->data, retLenDyn);
    gsl_spline_init(y_spline, seobdynamics->tVec, (*betaJ2P)->data, retLenDyn);

    precEulerparams.alpha_spline = x_spline;
    precEulerparams.alpha_acc = x_acc;
    precEulerparams.beta_spline = y_spline;
    precEulerparams.beta_acc = y_acc;

    precEulerF.function = &f_alphadotcosi;
    precEulerF.params = &precEulerparams;

    for (i = 0; i < retLenDyn; i++) {
      if (i == 0) {
        (*gammaJ2P)->data[i] = InitialGamma;
      } else {
        gsl_integration_qags(&precEulerF, tVec[i - 1], tVec[i], 1e-9, 1e-9,
                             1000, precEulerw, &precEulerresult,
                             &precEulererror);
        (*gammaJ2P)->data[i] = (*gammaJ2P)->data[i - 1] + precEulerresult;
      }
    }
    gsl_integration_workspace_free(precEulerw);
    gsl_spline_free(x_spline);
    gsl_spline_free(y_spline);
    gsl_interp_accel_free(x_acc);
    gsl_interp_accel_free(y_acc);

    /* Unwrap gamma in-place -- note that with integration, this should not be
     * necessary */
    XLALREAL8VectorUnwrapAngle((*gammaJ2P), (*gammaJ2P));
  }

  return XLAL_SUCCESS;
}

/**
 * This function extends Euler angles J2P according to the prescription
 * flagEulerextension after attachment point Two prescriptions: constant angles,
 * or simple precession around Jfinal at rate omegaQNM220-omegaQNM210 If
 * SpinsAlmostAligned, all Euler angles are set to 0
 */
// NOTE: the extension is based on l=2 qualitative behaviour in NR
// Has not been investigated for generic l
static int SEOBEulerJ2PPostMergerExtension(
    REAL8Vector
        *alphaJ2P, /**<< Output: vector for alpha J2P, already allocated */
    REAL8Vector
        *betaJ2P, /**<< Output: vector for beta J2P, already allocated */
    REAL8Vector
        *gammaJ2P, /**<< Output: vector for gamma J2P, already allocated */
    COMPLEX16
        sigmaQNM220, /**<< Input: complex frequency for QNM 22, 0th overtone */
    COMPLEX16
        sigmaQNM210, /**<< Input: complex frequency for QNM 21, 0th overtone */
    REAL8Vector *tVec, /**<< Input: time vector for Euler angles data (length of
                          P-modes) */
    UINT4 retLen,      /**<< Input: total length of Euler angles data (length of
                          P-modes) */
    UINT4 indexStart, /**<< Input: index where we start the extension (included,
                         index of time of attachment) */
    SpinEOBParams *seobParams, /**<< SEOB params */
    flagSEOBNRv4P_euler_extension
        flagEulerextension, /**<< flag indicating how to extend the Euler angles
                              post-merger */
    INT4 flip /** << a flag of whether to flip the sign of the precession frequency
              */
) {
  UINT4 i;
  UINT4 SpinsAlmostAligned = seobParams->alignedSpins;

  if (!SpinsAlmostAligned) {

    /* Initial values, that were computed from the dynamics */
    REAL8 timeAttach = tVec->data[indexStart - 1];
    REAL8 alphaAttach = alphaJ2P->data[indexStart - 1];
    REAL8 betaAttach = betaJ2P->data[indexStart - 1];
    REAL8 gammaAttach = gammaJ2P->data[indexStart - 1];

    if (flagEulerextension == FLAG_SEOBNRv4P_EULEREXT_QNM_SIMPLE_PRECESSION) {
      /* Precession rate around final J */
      REAL8 omegaQNM220 = creal(sigmaQNM220);
      REAL8 omegaQNM210 = creal(sigmaQNM210);
      REAL8 precRate = omegaQNM220 - omegaQNM210;
      // flip is either 1 or -1. This is needed because we want the sign of the precession frequency
      // to be correct even when the projected final spin is negative, at which point the QNMs
      // flip sign
      precRate *= flip;
      REAL8 cosbetaAttach = cos(betaAttach);
      for (i = indexStart; i < retLen; i++) {
        alphaJ2P->data[i] =
            alphaAttach + (tVec->data[i] - timeAttach) * precRate;
        betaJ2P->data[i] = betaAttach;
        gammaJ2P->data[i] = gammaAttach - cosbetaAttach *
                                              (tVec->data[i] - timeAttach) *
                                              precRate;
      }
    }

    else if (flagEulerextension == FLAG_SEOBNRv4P_EULEREXT_CONSTANT) {
      for (i = indexStart; i < retLen; i++) {
        alphaJ2P->data[i] = alphaAttach;
        betaJ2P->data[i] = betaAttach;
        gammaJ2P->data[i] = gammaAttach;
      }
    }

    else {
      XLALPrintError("XLAL Error - %s: flagEulerextension not recognized.\n",
                     __func__);
      XLAL_ERROR(XLAL_EINVAL);
    }

  }

  else /* if spins are almost aligned, set all Euler angles to 0 */
  {
    memset(&(alphaJ2P->data[indexStart]), 0,
           (retLen - indexStart) * sizeof(REAL8));
    memset(&(betaJ2P->data[indexStart]), 0,
           (retLen - indexStart) * sizeof(REAL8));
    memset(&(gammaJ2P->data[indexStart]), 0,
           (retLen - indexStart) * sizeof(REAL8));
  }

  return XLAL_SUCCESS;
}

/**
 * These functions compute the amplitude and phase of a Wigner coefficient
 * Dlmmp, given Euler angles of an active rotation.
 */
/**
 * Convention for Wigner matrices (mp stands for m', * for conjugation):
 * for the active rotation from the I-frame to the P-frame, parameterized by
 * Euler angles (alpha, beta, gamma) in the ZYZ convention \f[ h^{P}_{lm} =
 * \sum_{m'} D^{l}_{m m'} h^{I}_{lm'}\f] \f[ h^{I}_{lm} = \sum_{m'} D^{l}_{m
 * m'}* h^{P}_{lm'}\f] \f[ D^{l}_{m m'} = d^{l}_{m m'}(\beta) \exp{i m \alpha}
 * \exp{i m' \gamma}\f] with the notation \f$ c,s = \cos, \sin (\beta/2)\f$, \f$
 * k_{min} = \max(0,m-m'), k_{max} = \min(l+m, l-m')\f$:
 *
 * \f[d^{l}_{m m'}(\beta) = \sum_{k=k_{min}}^{k_{max}} \frac{(-1)^{k+m'-m}}{k!}
 * \frac{\sqrt{(l+m)!(l-m)!(l+m')!(l-m')!}}{(l+m-k)!(l-m'-k)!(k-m+m')!}
 * c^{2l+m-m'-2k} s^{2k-m+m'}\f]
 */
static REAL8 SEOBWignerDAmp(UINT4 l, INT4 m, INT4 mp, REAL8 beta) {
  return XLALWignerdMatrix(l, m, mp, beta);
}
static REAL8 SEOBWignerDPhase(INT4 m, INT4 mp, REAL8 alpha, REAL8 gamma) {
  return m * alpha + mp * gamma;
}

/**
 * This function computes the hJlm Re/Im timeseries (fixed sampling) from hPlm
 * amp/phase modes and Euler angles J2P (irregular sampling). The Wigner
 * rotation coefficients are first computed on the irregular sampling. We assume
 * that the hPlm and Euler angles are all given on a common time vector. P-frame
 * modes hPlmp should have mp>0 only We will generate the modes mp<0 using the
 * symmetry relation hPl-mp = (-1)^l hPlmp* All modes with l<=modes_lmax will be
 * created in output, and will be 0 for values of l possibly absent from
 * listhPlm
 */

static int SEOBRotateInterpolatehJlmReImFromSphHarmListhPlmAmpPhase(

    SphHarmTimeSeries **hJlm, /**<< Output: hJlm time series, will contain
                                 complex values on fixed sampling */
    INT4 modes[][2],          /**<< Input: array of modes (l,m) */
    UINT4 nmodes,             /**<< Input: number of modes (l,m) */
    INT4 modes_lmax,          /**<< Input: maximum value of l in modes (l,m) */
    REAL8 deltaT,             /**<< Input: time step for the hJlm timeseries */
    UINT4 retLenTS, /**<< Input: number of samples for the hJlm timeseries */
    REAL8Vector *tVecPmodes, /**<< Input: irregular time vector on which the
                                hPlm and Euler angles are given */
    SphHarmListCAmpPhaseSequence
        *listhPlm,         /**<< Input: list of P-frame modes hPlm */
    REAL8Vector *alphaJ2P, /**<< Input: vector for Euler angle alpha J2P */
    REAL8Vector *betaJ2P,  /**<< Input: vector for Euler angle beta J2P */
    REAL8Vector *gammaJ2P, /**<< Input: vector for Euler angle gamma J2P */
    UINT4 flagSymmetrizehPlminusm /**<< Input: flag indicating wether the
                                     P-frame modes m<0 are generated with the
                                     symmetry hP_l-m ~ (-1)^l hP_lm* */
) {
  UINT4 i;
  INT4 l, m, mp;
  REAL8 t = 0., camp_real = 0., camp_imag = 0., phase = 0., alpha = 0.,
        beta = 0., gamma = 0.;
  UINT4 retLenP = tVecPmodes->length;

  /* Create output time vector */
  REAL8Vector *timeTS = XLALCreateREAL8Vector(retLenTS);
  REAL8 *timeTSdata = timeTS->data;
  for (i = 0; i < retLenTS; i++) {
    timeTSdata[i] = i * deltaT;
  }

  /* Create output list of timeseries, with all (l,m) up to modes_lmax */
  *hJlm = NULL;
  LIGOTimeGPS tGPS = LIGOTIMEGPSZERO;
  char mode_string[32];
  for (l = 2; l <= modes_lmax; l++) {
    for (m = -l; m <= l; m++) {
      sprintf(mode_string, "H_%d%d", l, m);
      COMPLEX16TimeSeries *hJlm_TS = XLALCreateCOMPLEX16TimeSeries(
          mode_string, &tGPS, 0., deltaT, &lalStrainUnit, retLenTS);
      memset(hJlm_TS->data->data, 0, retLenTS * sizeof(COMPLEX16));

      /* Note: with the AddMode function, data is copied over */
      *hJlm = XLALSphHarmTimeSeriesAddMode(*hJlm, hJlm_TS, l, m);

      /* Data has been copied over, we need to destroy */
      XLALDestroyCOMPLEX16TimeSeries(hJlm_TS);
    }
  }

  /* Set time data */
  XLALSphHarmTimeSeriesSetTData(*hJlm, timeTS);

  /* Create working space for amp/phase of Wigner coefficient */
  REAL8Vector *Dlmmp_amp = XLALCreateREAL8Vector(retLenP);
  REAL8Vector *Dlmmp_phase = XLALCreateREAL8Vector(retLenP);
  REAL8Vector *Dlmminusmp_amp = XLALCreateREAL8Vector(retLenP);
  REAL8Vector *Dlmminusmp_phase = XLALCreateREAL8Vector(retLenP);
  memset(Dlmmp_amp->data, 0, (Dlmmp_amp->length) * sizeof(REAL8));
  memset(Dlmmp_phase->data, 0, (Dlmmp_phase->length) * sizeof(REAL8));
  memset(Dlmminusmp_amp->data, 0, (Dlmminusmp_amp->length) * sizeof(REAL8));
  memset(Dlmminusmp_phase->data, 0, (Dlmminusmp_phase->length) * sizeof(REAL8));

  /* Interpolating splines for Wigner coefficients */
  gsl_spline *spline_Dlmmp_amp = gsl_spline_alloc(gsl_interp_cspline, retLenP);
  gsl_spline *spline_Dlmmp_phase =
      gsl_spline_alloc(gsl_interp_cspline, retLenP);
  gsl_spline *spline_Dlmminusmp_amp =
      gsl_spline_alloc(gsl_interp_cspline, retLenP);
  gsl_spline *spline_Dlmminusmp_phase =
      gsl_spline_alloc(gsl_interp_cspline, retLenP);
  gsl_interp_accel *accel_Dlmmp_amp = gsl_interp_accel_alloc();
  gsl_interp_accel *accel_Dlmmp_phase = gsl_interp_accel_alloc();
  gsl_interp_accel *accel_Dlmminusmp_amp = gsl_interp_accel_alloc();
  gsl_interp_accel *accel_Dlmminusmp_phase = gsl_interp_accel_alloc();

  /* Interpolating splines for hPlm modes */
  gsl_spline *spline_camp_real = gsl_spline_alloc(gsl_interp_cspline, retLenP);
  gsl_spline *spline_camp_imag = gsl_spline_alloc(gsl_interp_cspline, retLenP);
  gsl_spline *spline_phase = gsl_spline_alloc(gsl_interp_cspline, retLenP);
  gsl_interp_accel *accel_camp_real = gsl_interp_accel_alloc();
  gsl_interp_accel *accel_camp_imag = gsl_interp_accel_alloc();
  gsl_interp_accel *accel_phase = gsl_interp_accel_alloc();

  /* Interpolate P-frame modes hPlm on the time samples needed for the output as
   * a time series */
  SphHarmListCAmpPhaseSequence *listhPlm_TS = NULL;
  /* Loop over modes */
  for (UINT4 nmode = 0; nmode < nmodes; nmode++) {

    l = modes[nmode][0];
    m = modes[nmode][1];

    CAmpPhaseSequence *hPlm =
        SphHarmListCAmpPhaseSequence_GetMode(listhPlm, l, m)->campphase;

    CAmpPhaseSequence *hPlm_TS = NULL;
    if (CAmpPhaseSequence_Init(&hPlm_TS, retLenTS) == XLAL_FAILURE) {
      XLALPrintError("XLAL Error - %s: failure in CAmpPhaseSequence_Init for "
                     "mode (l,m) = (%d,%d).\n",
                     __func__, l, m);
      XLAL_ERROR(XLAL_ENOMEM);
    }

    gsl_spline_init(spline_camp_real, tVecPmodes->data, hPlm->camp_real->data,
                    retLenP);
    gsl_spline_init(spline_camp_imag, tVecPmodes->data, hPlm->camp_imag->data,
                    retLenP);
    gsl_spline_init(spline_phase, tVecPmodes->data, hPlm->phase->data, retLenP);

    COMPLEX16 hPlm_val = 0.;
    REAL8 *hPlmTS_tdata = hPlm_TS->xdata->data;
    REAL8 *hPlmTS_camprealdata = hPlm_TS->camp_real->data;
    REAL8 *hPlmTS_campimagdata = hPlm_TS->camp_imag->data;
    REAL8 *hPlmTS_phasedata = hPlm_TS->phase->data;
    for (i = 0; i < retLenTS; i++) {
      t = timeTSdata[i];
      /* Here we include a possible imaginary part of the complex envelope, but
       * at the moment it is simply 0 (only real part is used) */
      camp_real = gsl_spline_eval(spline_camp_real, t, accel_camp_real);
      camp_imag = gsl_spline_eval(spline_camp_imag, t, accel_camp_imag);
      phase = gsl_spline_eval(spline_phase, t, accel_phase);
      hPlm_val = (camp_real + I * camp_imag) * cexp(I * phase);
      /* We output the interpolated value for the mode hPlm in Re/Im form,
       * setting the phase to 0 */
      hPlmTS_tdata[i] = t;
      hPlmTS_camprealdata[i] = creal(hPlm_val);
      hPlmTS_campimagdata[i] = cimag(hPlm_val);
      hPlmTS_phasedata[i] = 0.;
    }

    SphHarmListCAmpPhaseSequence_AddMode(&listhPlm_TS, hPlm_TS, l, m);
  }

  /* Main computation */
  /* hJlm = \sum_mp Dlmmpstar hPlmp */

  REAL8 *Dlmmp_amp_data = Dlmmp_amp->data;
  REAL8 *Dlmmp_phase_data = Dlmmp_phase->data;
  REAL8 *Dlmminusmp_amp_data = Dlmminusmp_amp->data;
  REAL8 *Dlmminusmp_phase_data = Dlmminusmp_phase->data;
  REAL8 *alphadata = alphaJ2P->data;
  REAL8 *betadata = betaJ2P->data;
  REAL8 *gammadata = gammaJ2P->data;
  COMPLEX16TimeSeries *hJlmmode = NULL;
  COMPLEX16 *hJlmmode_data = NULL;
  REAL8 *hPlmp_campreal_data;
  REAL8 *hPlmp_campimag_data;
  COMPLEX16 hPlmp_val;
  COMPLEX16 Dlmmp_amp_val, Dlmmp_phase_val, Dlmmp_val;
  COMPLEX16 Dlmminusmp_amp_val, Dlmminusmp_phase_val, Dlmminusmp_val;

  /* Loop on l */
  for (l = 2; l <= modes_lmax; l++) {

    /* Loop on m */
    for (m = -l; m <= l; m++) {

      /* Get hJlm mode */
      hJlmmode = XLALSphHarmTimeSeriesGetMode(*hJlm, l, m);
      hJlmmode_data = hJlmmode->data->data;

      /* Loop on the modes hPlmp */
      for (UINT4 nmode = 0; nmode < nmodes; nmode++) {

        /* Select modes with the same l */
        if (modes[nmode][0] != l)
          continue;

        mp = modes[nmode][1];

        /* We do not allow mp<=0 in the P-frame modes hPlmp */
        /* We will use hPl-mp = (-1)^l hPlmp* */
        if (mp <= 0) {
          XLALPrintError("XLAL Error - %s: mode (l,mp) = (%d,%d) is not "
                         "allowed as mp<=0.\n",
                         __func__, l, mp);
          XLAL_ERROR(XLAL_EINVAL);
        }

        /* Get interpolated mode hPlmp */
        CAmpPhaseSequence *hPlmp_TS =
            SphHarmListCAmpPhaseSequence_GetMode(listhPlm_TS, l, mp)->campphase;
        hPlmp_campreal_data = hPlmp_TS->camp_real->data;
        hPlmp_campimag_data = hPlmp_TS->camp_imag->data;

        /* Compute Wigner coefficients amp/phase on the same input sampling as
         * the P-frame modes */
        for (i = 0; i < retLenP; i++) {
          alpha = alphadata[i];
          beta = betadata[i];
          gamma = gammadata[i];
          Dlmmp_amp_data[i] = SEOBWignerDAmp(l, m, mp, beta);
          Dlmmp_phase_data[i] = SEOBWignerDPhase(m, mp, alpha, gamma);
          if (flagSymmetrizehPlminusm) {
            Dlmminusmp_amp_data[i] = SEOBWignerDAmp(l, m, -mp, beta);
            Dlmminusmp_phase_data[i] = SEOBWignerDPhase(m, -mp, alpha, gamma);
          }
        }

        /* Interpolate amplitude/phase of the Wigner coefficient, add
         * contribution to hJlm mode */
        gsl_spline_init(spline_Dlmmp_amp, tVecPmodes->data, Dlmmp_amp->data,
                        retLenP);
        gsl_spline_init(spline_Dlmmp_phase, tVecPmodes->data, Dlmmp_phase->data,
                        retLenP);
        if (flagSymmetrizehPlminusm) {
          gsl_spline_init(spline_Dlmminusmp_amp, tVecPmodes->data,
                          Dlmminusmp_amp->data, retLenP);
          gsl_spline_init(spline_Dlmminusmp_phase, tVecPmodes->data,
                          Dlmminusmp_phase->data, retLenP);
        }
        for (i = 0; i < retLenTS; i++) {
          t = timeTSdata[i];
          Dlmmp_amp_val = gsl_spline_eval(spline_Dlmmp_amp, t, accel_Dlmmp_amp);
          Dlmmp_phase_val =
              gsl_spline_eval(spline_Dlmmp_phase, t, accel_Dlmmp_phase);
          Dlmmp_val =
              Dlmmp_amp_val *
              cexp(-I * Dlmmp_phase_val); /* mind the conjugation Dlmmpstar */
          hPlmp_val = hPlmp_campreal_data[i] + I * hPlmp_campimag_data[i];
          hJlmmode_data[i] += Dlmmp_val * hPlmp_val;
          if (flagSymmetrizehPlminusm) {
            Dlmminusmp_amp_val =
                gsl_spline_eval(spline_Dlmminusmp_amp, t, accel_Dlmminusmp_amp);
            Dlmminusmp_phase_val = gsl_spline_eval(spline_Dlmminusmp_phase, t,
                                                   accel_Dlmminusmp_phase);
            Dlmminusmp_val =
                Dlmminusmp_amp_val *
                cexp(-I * Dlmminusmp_phase_val); /* mind the conjugation
                                                    Dlmminusmpstar */
            hJlmmode_data[i] += Dlmminusmp_val * pow(-1, l) * conj(hPlmp_val);
          }
        }
      }
    }
  }

  /* Cleanup */
  SphHarmListCAmpPhaseSequence_Destroy(listhPlm_TS);
  gsl_spline_free(spline_camp_real);
  gsl_spline_free(spline_camp_imag);
  gsl_spline_free(spline_phase);
  gsl_interp_accel_free(accel_camp_real);
  gsl_interp_accel_free(accel_camp_imag);
  gsl_interp_accel_free(accel_phase);
  gsl_spline_free(spline_Dlmmp_amp);
  gsl_spline_free(spline_Dlmmp_phase);
  gsl_spline_free(spline_Dlmminusmp_amp);
  gsl_spline_free(spline_Dlmminusmp_phase);
  gsl_interp_accel_free(accel_Dlmmp_amp);
  gsl_interp_accel_free(accel_Dlmmp_phase);
  gsl_interp_accel_free(accel_Dlmminusmp_amp);
  gsl_interp_accel_free(accel_Dlmminusmp_phase);
  XLALDestroyREAL8Vector(Dlmmp_amp);
  XLALDestroyREAL8Vector(Dlmmp_phase);
  XLALDestroyREAL8Vector(Dlmminusmp_amp);
  XLALDestroyREAL8Vector(Dlmminusmp_phase);

  return XLAL_SUCCESS;
}

/**
 * This function computes the hIlm Re/Im timeseries (fixed sampling) from hJlm
 * Re/Im timeseries (same sampling). This is a simple rotation,
 * sample-by-sample, with constant Wigner coefficients.
 * See the comment before SEOBWignerDAmp for explanation of conventions,
 * and Appendix A of Babak et al, Phys. Rev. D 95, 024010, 2017 [arXiv:1607.05661] for a general
 * discussion.
 */
static int SEOBRotatehIlmFromhJlm(
    SphHarmTimeSeries **hIlm, /**<< Output: hIlm time series, complex values on
                                 fixed sampling */
    SphHarmTimeSeries *hJlm,  /**<< Output: hJlm time series, complex values on
                                 fixed sampling */
    INT4 modes_lmax,          /**<< Input: maximum value of l in modes (l,m) */
    REAL8 alphaI2J,           /**<< Input: Euler angle alpha I->J */
    REAL8 betaI2J,            /**<< Input: Euler angle beta I->J */
    REAL8 gammaI2J,           /**<< Input: Euler angle gamma I->J */
    REAL8 deltaT              /**<< Input: time step, necessary to initialize new timeseries */
) {

  UINT4 i;
  INT4 l, m, mp;
  REAL8 amp_wigner = 0., phase_wigner = 0.;
  COMPLEX16 D_wigner = 0.;
  UINT4 retLen = hJlm->tdata->length;
  REAL8 *tJdata = hJlm->tdata->data;

  /* Copy time vector */
  REAL8Vector *tI = XLALCreateREAL8Vector(retLen);
  memcpy(tI->data, tJdata, retLen * sizeof(REAL8));

  /* Create output list of timeseries, with all (l,m) up to modes_lmax */
  *hIlm = NULL;
  LIGOTimeGPS tGPS = LIGOTIMEGPSZERO;
  char mode_string[32];
  for (l = 2; l <= modes_lmax; l++) {
    for (m = -l; m <= l; m++) {
      sprintf(mode_string, "H_%d%d", l, m);
      COMPLEX16TimeSeries *hIlm_TS = XLALCreateCOMPLEX16TimeSeries(
          mode_string, &tGPS, 0., deltaT, &lalStrainUnit, retLen);
      memset(hIlm_TS->data->data, 0, retLen * sizeof(COMPLEX16));

      /* Note: with the AddMode function, data is copied over */
      *hIlm = XLALSphHarmTimeSeriesAddMode(*hIlm, hIlm_TS, l, m);

      /* Data has been copied over, we need to destroy */
      XLALDestroyCOMPLEX16TimeSeries(hIlm_TS);
    }
  }

  /* Set time data */
  XLALSphHarmTimeSeriesSetTData(*hIlm, tI);

  /* Main computation */
  /* hIlm = \sum_mp Dlmpm hJlmp */

  COMPLEX16TimeSeries *hIlmmode = NULL;
  COMPLEX16TimeSeries *hJlmpmode = NULL;
  COMPLEX16 *hIlmmode_data = NULL;
  COMPLEX16 *hJlmpmode_data = NULL;

  /* Loop on l */
  for (l = 2; l <= modes_lmax; l++) {

    /* Loop on m */
    for (m = -l; m <= l; m++) {

      /* Get hJlm mode */
      hIlmmode = XLALSphHarmTimeSeriesGetMode(*hIlm, l, m);
      hIlmmode_data = hIlmmode->data->data;

      /* Loop on mp - exclude value 0, since hPl0=0 in our approximation */
      for (mp = -l; mp <= l; mp++) {

        /* Get hJlm mode */
        hJlmpmode = XLALSphHarmTimeSeriesGetMode(hJlm, l, mp);
        hJlmpmode_data = hJlmpmode->data->data;

        /* Compute constant Wigner coefficient */
        amp_wigner = SEOBWignerDAmp(l, m, mp, betaI2J);
        phase_wigner = SEOBWignerDPhase(m, mp, alphaI2J, gammaI2J);
        D_wigner = amp_wigner *
                   cexp(-I * phase_wigner); /* mind the conjugation Dlmmpstar */

        /* Evaluate mode contribution */
        for (i = 0; i < retLen; i++) {
          hIlmmode_data[i] += D_wigner * hJlmpmode_data[i];
        }
      }
    }
  }

  return XLAL_SUCCESS;
}

/**
 * This function combines the modes hIlm timeseries with the sYlm to produce the
 * polarizations hplus, hcross.
 */
// NOTE: azimuthal angle of the observer entering the -2Ylm is pi/2-phi
// according to LAL conventions
static int SEOBComputehplushcrossFromhIlm(
    REAL8TimeSeries
        *hplusTS, /**<< Output: time series for hplus, already created */
    REAL8TimeSeries
        *hcrossTS,   /**<< Output: time series for hplus, already created */
    INT4 modes_lmax, /**<< Input: maximum value of l */
    SphHarmTimeSeries
        *hIlm,  /**<< Input: list with time series for each mode hIlm */
    REAL8 amp0, /**<< Input: amplitude prefactor */
    REAL8 inc,  /**<< Input: inclination */
    REAL8 phi   /**<< Input: phase */
) {
  INT4 l, m;

  /* hplus, hcross */
  REAL8 *hplusdata = hplusTS->data->data;
  memset(hplusdata, 0, hplusTS->data->length * sizeof(REAL8));
  REAL8 *hcrossdata = hcrossTS->data->data;
  memset(hcrossdata, 0, hplusTS->data->length * sizeof(REAL8));
  COMPLEX16 hpc_contrib = 0.;

  /* Loop over modes */
  for (l = 2; l <= modes_lmax; l++) {
    for (m = -l; m <= l; m++) {

      /* Compute sYlm */
      COMPLEX16 sYlm =
          XLALSpinWeightedSphericalHarmonic(inc, LAL_PI / 2. - phi, -2, l, m);

      /* Get mode hIlm */
      COMPLEX16TimeSeries *hIlmTS = XLALSphHarmTimeSeriesGetMode(hIlm, l, m);
      COMPLEX16 *hIlm_data = hIlmTS->data->data;

      for (UINT4 i = 0; i < hplusTS->data->length; i++) {
        hpc_contrib = sYlm * hIlm_data[i];
        hplusdata[i] += amp0 * creal(hpc_contrib);
        hcrossdata[i] += -amp0 * cimag(hpc_contrib);
      }
    }
  }

  return XLAL_SUCCESS;
}

/**
 * This function gets the number of modes present in a mode array, ignoring
 * modes l>lmax (and l<2) Note that m<=0 modes are also ignored and a warning
 * given
 */
static UINT4 SEOBGetNumberOfModesInModeArray(
    LALValue *modearray, /**<< Input: ModeArray structure */
    int lmax /**<< Input: maximum value of l to explore -- possible modes with
                l>lmax will be ignored */
) {
  UINT4 nmodes = 0;
  for (INT4 l = 2; l <= lmax; l++) {
    for (INT4 m = -l; m <= l; m++) {
      if (XLALSimInspiralModeArrayIsModeActive(modearray, l, m)) {
        if (m > 0)
          nmodes++;
        else {
          XLAL_PRINT_WARNING(
              "XLAL Warning - %s: mode (l,m)=(%d,%d) present in mode array -- "
              "in our conventions, we only consider m>0. Mode ignored.\n",
              __func__, l, m);
        }
      }
    }
  }
  return nmodes;
}
/**
 * This function populates a dynamically allocated INT4 array to with the modes
 * active in a ModeArray Possible modes with l>lmax are ignored (as well as l<2)
 * Note that m<=0 modes are also ignored and a warning given
 */
static int SEOBGetModeNumbersFromModeArray(
    INT4 modes[][2],     /**<< Output: array of dimension (nmodes,2) with mode
                            numbers (l,m) */
    LALValue *modearray, /**<< Input: ModeArray structure */
    int lmax /**<< Input: maximum value of l to explore -- possible modes with
                l>lmax will be ignored */
) {
  /* Populate array */
  UINT4 nmode = 0;
  for (INT4 l = 2; l <= lmax; l++) {
    for (INT4 m = l; m >= -l; m--) {
      if (XLALSimInspiralModeArrayIsModeActive(modearray, l, m)) {
        if (m > 0) {
          modes[nmode][0] = l;
          modes[nmode][1] = m;
          nmode++;
        } else {
          XLAL_PRINT_WARNING(
              "XLAL Warning - %s: mode (l,m)=(%d,%d) present in mode array -- "
              "in our conventions, we only consider m>0. Mode ignored.\n",
              __func__, l, m);
        }
      }
    }
  }

  return XLAL_SUCCESS;
}

/* Check that the ringdown frequency of the highest ell mode is less than the
 * Nyquist frequency */
int XLALEOBCheckNyquistFrequency(REAL8 m1, REAL8 m2, REAL8 spin1[3],
                                 REAL8 spin2[3], UINT4 ell_max,
                                 Approximant approx, REAL8 deltaT) {
  UINT4 mode_highest_freqL = ell_max;
  UINT4 mode_highest_freqM = ell_max;
  /* Ringdown freq used to check the sample rate */
  COMPLEX16Vector modefreqVec;
  COMPLEX16 modeFreq;
  modefreqVec.length = 1;
  modefreqVec.data = &modeFreq;

  if (XLALSimIMREOBGenerateQNMFreqV2Prec(&modefreqVec, m1, m2, spin1, spin2,
                                         mode_highest_freqL, mode_highest_freqM,
                                         1, approx) == XLAL_FAILURE) {
    XLAL_ERROR(XLAL_EFUNC);
  }

  if (deltaT > LAL_PI / creal(modeFreq)) {
    XLALPrintError("XLAL Error - %s: Ringdown frequency > Nyquist "
                   "frequency!\nAt present this situation is not supported.\n",
                   __func__);

    XLAL_ERROR(XLAL_EDOM);
  }
  return XLAL_SUCCESS;
}

/**
 * This function is the master function generating precessing spinning SEOBNRv4P
 * waveforms h+ and hx. Currently, only h2m harmonics will be generated.
 *
 * Input conventions:
 * Cartesian coordinate system: initial \f$\vec{L}_N\f$ is in the xz plane,
 * rotated away from the z-axis by an angle inc phiC       : in radians deltaT
 * : in SI units (s) m1SI, m2SI : in SI units (kg) fMin       : in SI units (Hz)
 * r          : in SI units (m)
 * inc        : in radians
 * INspin{1,2}: in dimensionless units of m{1,2}^2
 *
 * Evolution conventions:
 * values[0-2]: r vector in units of total mass
 * values[3-5]: pstar vector in units of reduced mass
 * values[6-8]: S1 vector in units of (total mass)^2
 * values[9-11]: S2 vector in units of (total mass)^2
 * values[12-13]: phases (carrier and precession (Thomas)) in rads
 *
 * Note that when the initial opening angles of the spins w.r.t. the initial
 * Newtonian angular momentum are small, the aligned-spin SEOBNRv4 dynamics is
 * used. However, the waveforms are then generated according to the SEOBNRv4P
 * model: for example, they include the (2,1) mode.
 *
 * STEP 0)  Prepare parameters, including pre-computed coefficients for EOB
 * Hamiltonian, flux and waveform STEP 1)  Solve for initial conditions STEP 2)
 * Evolve EOB trajectory with adaptive sampling (AdaS) STEP 3)  Step back and
 * evolve EOB trajectory at high sampling rate (HiS) STEP 4)  Get final
 * J/L/spins from HiS dynamics at peak of Omega, compute constant angles
 * EulerI2J STEP 5)  Compute P-frame 22 mode amp/phase on HiS and compute NQC
 * coefficients STEP 6)  Compute P-frame amp/phase for all modes on HiS, now
 * including NQC STEP 7)  Attach RD to the P-frame waveform STEP 8)  Build the
 * joined dynamics AdaS+HiS up to attachment, joined P-modes AdaS+HiS+RDpatch
 * STEP 9)  Compute Euler angles J2P from AdaS and HiS dynamics up to attachment
 * STEP 10) Compute Euler angles J2P extension after attachment
 * STEP 11) Compute modes hJlm on the output time series by rotating and
 * interpolating the modes hPlm STEP 12) Rotate waveform from J-frame to the
 * output I-frame on timeseries-sampling (constant Wigner coeffs) STEP 13)
 * Compute hplus, hcross from I-frame waveform on timeseries sampling STEP -1)
 * Cleanup
 */
int XLALSimIMRSpinPrecEOBWaveformAll(
    REAL8TimeSeries *
        *hplus, /**<< Main output: hplus GW polarization time series */
    REAL8TimeSeries *
        *hcross, /**<< Main output: hcross GW polarization time series */
    SphHarmTimeSeries **hIlm, /**<< Spherical modes time series for the waveform
                                 in the initial inertial frame */
    SphHarmTimeSeries **hJlm, /**<< Spherical modes time series for the waveform
                                 in the final-J inertial frame */
    REAL8Vector **seobdynamicsAdaSVector,    /**<< Vector for extended dynamics
                                                values, adaptive sampling part */
    REAL8Vector **seobdynamicsHiSVector,     /**<< Vector for extended dynamics
                                                values, high sampling part */
    REAL8Vector **seobdynamicsAdaSHiSVector, /**<< Vector for extended joined
                                                dynamics values */
    REAL8Vector **tVecPmodes,   /**<< Time vector for the P-modes */
    REAL8Vector **hP22_amp,     /**<< Vector for the hP22 mode amplitude */
    REAL8Vector **hP22_phase,   /**<< Vector for the hP22 mode phase */
    REAL8Vector **hP21_amp,     /**<< Vector for the hP21 mode amplitude */
    REAL8Vector **hP21_phase,   /**<< Vector for the hP21 mode phase */
    REAL8Vector **hP33_amp,     /**<< Vector for the hP33 mode amplitude */
    REAL8Vector **hP33_phase,   /**<< Vector for the hP33 mode phase */
    REAL8Vector **hP44_amp,     /**<< Vector for the hP44 mode amplitude */
    REAL8Vector **hP44_phase,   /**<< Vector for the hP44 mode phase */
    REAL8Vector **hP55_amp,     /**<< Vector for the hP55 mode amplitude */
    REAL8Vector **hP55_phase,   /**<< Vector for the hP55 mode phase */
    REAL8Vector **alphaJ2P,     /**<< Vector for the Euler angle alphaJ2P */
    REAL8Vector **betaJ2P,      /**<< Vector for the Euler angle betaJ2P */
    REAL8Vector **gammaJ2P,     /**<< Vector for the Euler angle gammaJ2P */
    REAL8Vector **mergerParams, /**<< Parameters at merger */
    const REAL8 phi,            /**<< Input: phase */
    const REAL8 INdeltaT,       /**<< Input: sampling time step in SI */
    const REAL8 m1SI,           /**<< Input: mass of first object in SI */
    const REAL8 m2SI,           /**<< Input: mass of second object in SI */
    const REAL8 fMin,           /**<< Input: fMin */
    const REAL8 r,              /**<< Input: luminosity distance in SI */
    const REAL8 inc,            /**<< Input: inclination */
    const REAL8 chi1x,   /**<< Input: spin1 x-component, dimensionless (in units
                            of m1^2) */
    const REAL8 chi1y,   /**<< Input: spin1 y-component, dimensionless (in units
                            of m1^2) */
    const REAL8 chi1z,   /**<< Input: spin1 z-component, dimensionless (in units
                            of m1^2) */
    const REAL8 chi2x,   /**<< Input: spin2 x-component, dimensionless (in units
                            of m2^2) */
    const REAL8 chi2y,   /**<< Input: spin2 y-component, dimensionless (in units
                            of m2^2) */
    const REAL8 chi2z,   /**<< Input: spin2 z-component, dimensionless (in units
                            of m2^2) */
    LALValue *modearray, /**<< Input: mode array for the hlm modes (m>0) to be
                            generated in the P-frame -- has to include (2,2)
                            mode -- modes with l larger than _SEOB_MODES_LMAX
                            will be ignored */
    LALDict *seobflags   /**<< Input: dictionary of SEOB flags */
) {
  /* Read waveform flags */
  if (!seobflags) {
    XLALPrintError("XLAL Error - %s: pointer to LALDict seobflags is NULL.\n",
                   __func__);
    PRINT_ALL_PARAMS;
    XLAL_ERROR(XLAL_EINVAL);
  }
  /* Flag for the spin-aligned version of SEOB */
  /* Default v4 */
  INT4 SpinAlignedEOBversion = 0;
  if (!XLALDictContains(seobflags, "SEOBNRv4P_SpinAlignedEOBversion"))
    SpinAlignedEOBversion = 4;
  else
    SpinAlignedEOBversion =
        XLALDictLookupINT4Value(seobflags, "SEOBNRv4P_SpinAlignedEOBversion");
  /* Flag to generate P-frame modes m<0 with the symmetry hP_l-m ~ (-1)^l hP_lm*
   */
  /* Default True */
  INT4 flagSymmetrizehPlminusm = 0;
  if (!XLALDictContains(seobflags, "SEOBNRv4P_SymmetrizehPlminusm"))
    flagSymmetrizehPlminusm = 1;
  else
    flagSymmetrizehPlminusm =
        XLALDictLookupINT4Value(seobflags, "SEOBNRv4P_SymmetrizehPlminusm");
  /* Flag to use numerical or analytical derivatives of the Hamiltonian  */
  /* Default numerical */
  // NOTE: analytical would be better if possible
  INT4 flagHamiltonianDerivative = 0;
  if (!XLALDictContains(seobflags, "SEOBNRv4P_HamiltonianDerivative"))
    flagHamiltonianDerivative = FLAG_SEOBNRv4P_HAMILTONIAN_DERIVATIVE_NUMERICAL;
  else
    flagHamiltonianDerivative =
        XLALDictLookupINT4Value(seobflags, "SEOBNRv4P_HamiltonianDerivative");
  /* Flag to choose the extension of Euler angles post-merger, constant or
   * simple precession around final J at a rate set by QNMs */
  /* Default QNM simple precession */
  INT4 flagEulerextension = 0;
  if (!XLALDictContains(seobflags, "SEOBNRv4P_euler_extension"))
    flagEulerextension = FLAG_SEOBNRv4P_EULEREXT_QNM_SIMPLE_PRECESSION;
  else
    flagEulerextension =
        XLALDictLookupINT4Value(seobflags, "SEOBNRv4P_euler_extension");
  /* Flag to choose the Z-axis of the radiation frame, L or LN */
  /* Default L */
  INT4 flagZframe = 0;
  if (!XLALDictContains(seobflags, "SEOBNRv4P_Zframe"))
    flagZframe = FLAG_SEOBNRv4P_ZFRAME_L;
  else
    flagZframe = XLALDictLookupINT4Value(seobflags, "SEOBNRv4P_Zframe");
  /* Flag to print debug information */
  /* Default False */
  INT4 debug = 0;
  if (!XLALDictContains(seobflags, "SEOBNRv4P_debug"))
    debug = 0;
  else
    debug = XLALDictLookupINT4Value(seobflags, "SEOBNRv4P_debug");

  if (debug) {
    printf("************************************\n");
    printf("XLALSimIMRSpinPrecEOBWaveformAll\n");
    printf("************************************\n");
    printf("phi = %.16e\n", phi);
    printf("INdeltaT = %.16e\n", INdeltaT);
    printf("m1SI = %.16e\n", m1SI);
    printf("m2SI = %.16e\n", m2SI);
    printf("fMin = %.16e\n", fMin);
    printf("r = %.16e\n", r);
    printf("inc = %.16e\n", inc);
    printf("chi1x = %.16e\n", chi1x);
    printf("chi1y = %.16e\n", chi1y);
    printf("chi1z = %.16e\n", chi1z);
    printf("chi2x = %.16e\n", chi2x);
    printf("chi2y = %.16e\n", chi2y);
    printf("chi2z = %.16e\n", chi2z);
    printf("flagHamiltonianDerivative = %d\n", flagHamiltonianDerivative);
    printf("flagEulerextension = %d\n", flagEulerextension);
    printf("flagZframe = %d\n", flagZframe);
  }

  /******************************************************************************************************************/
  /* STEP 0)  Allocate memory, prepare parameters and pre-compute coeffs for EOB
   * Hamiltonian, flux and waveform */
  /******************************************************************************************************************/

  if (debug)
    printf("STEP 0)  Allocate memory, prepare parameters and pre-compute "
           "coeffs for EOB Hamiltonian, flux and waveform\n");

  /* Check that SpinAlignedEOBversion has an authorized value: v2 or v4 */
  if (!((SpinAlignedEOBversion == 2) || (SpinAlignedEOBversion == 4))) {
    XLALPrintError("XLAL Error - %s: unauthorized SpinAlignedEOBversion %d.\n",
                   __func__, SpinAlignedEOBversion);
    PRINT_ALL_PARAMS;
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* Initialize mass parameters */
  REAL8 m1 = m1SI / LAL_MSUN_SI;
  REAL8 m2 = m2SI / LAL_MSUN_SI;
  REAL8 mTotal = m1 + m2;
  REAL8 mTScaled = mTotal * LAL_MTSUN_SI;
  REAL8 eta = m1 * m2 / (mTotal * mTotal);
  if (eta > 0.25) {
    m2 = m1;
    eta = 0.25;
  }

  /* Initialize amplitude factor */
  REAL8 amp0 = mTotal * LAL_MRSUN_SI / r;

  /* Convert input dimensionless spin components to vectors */
  REAL8Vector INchi1;
  REAL8Vector INchi2;
  REAL8 INchi1data[3] = {chi1x, chi1y, chi1z};
  REAL8 INchi2data[3] = {chi2x, chi2y, chi2z};
  INchi1.length = INchi2.length = 3;
  INchi1.data = INchi1data;
  INchi2.data = INchi2data;
  INT4 ret = 0;

  // Check Nyquist frequency
  UINT4 ellMaxForNyquistCheck = 2;
  UINT4 ell_max = SEOBGetLMaxInModeArray(modearray, _SEOB_MODES_LMAX);
  
  // Set the max ell to use for Nyquist check
  // If the value is not set, simply use the largest L in the mode array
  // which will give 2 for SEOBNRv4P and 5 for SEOBNRv4PHM
  if(!XLALDictContains(seobflags,"ellMaxForNyquistCheck")){
    ellMaxForNyquistCheck = ell_max;
  }
  else{
     ellMaxForNyquistCheck = XLALDictLookupINT4Value(seobflags, "ellMaxForNyquistCheck");
  }
  if (ellMaxForNyquistCheck<2){
    XLALPrintError("Custom value of ell < 2 was passed to Nyquist check. This is not supported!");
    XLAL_ERROR(XLAL_EFUNC);
  }
  if(ellMaxForNyquistCheck < ell_max){
    XLALPrintError("WARNING: Using ell=%d for Nyqusit check of srate, even though max ell of waveforms produced is %d\n",ellMaxForNyquistCheck,ell_max);
  }
  XLAL_TRY(XLALEOBCheckNyquistFrequency(m1, m2, INchi1.data, INchi2.data,
                                       ellMaxForNyquistCheck, SEOBNRv4P, INdeltaT),
           ret);
  if (ret != XLAL_SUCCESS) {
    XLAL_ERROR(XLAL_EDOM);
  }
  /* Sanity check various inputs */
  if (sqrt(chi1x * chi1x + chi1y * chi1y + chi1z * chi1z) > 1 ||
      sqrt(chi2x * chi2x + chi2y * chi2y + chi2z * chi2z) > 1) {
    XLALPrintError("XLAL Error - %s: Spins have magnitude > 1 !\n", __func__);
    PRINT_ALL_PARAMS;
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (m1 < 0 || m2 < 0) {
    XLALPrintError("XLAL Error - %s: One of the masses is < 0 !\n", __func__);
    PRINT_ALL_PARAMS;
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (m1 + m2 < 2) {
    printf("Warning: you are attempting to generate waveforms for BHs with "
           "total mass < 2 M_Sun. This may take a *very* long time\n");
  }
  if (m1 / m2 > 100 || m1 / m2 < 1 / 100) {
    XLALPrintError("XLAL Error - %s: Mass ratio is  > 100 !\n", __func__);
    PRINT_ALL_PARAMS;
    XLAL_ERROR(XLAL_EDOM);
  }
  /* If both dimensionless spin in-plane components are smaller than this
   * threshold, we fall back to aligned spins */
  REAL8 EPS_ALIGN = 1.0e-4;

  /* Check if initial frequency is too high: we choose an initial minimum
   * separation of 10M as a compromise between reliability of initial conditions
   * and length of the waveform. We use Newtonian Kepler's law. Refuse to
   * generate waveforms shorter than that.
   */
  REAL8 freqMinRad = 0;
  /*Compute the highest initial frequency of the 22 mode */
  XLALEOBHighestInitialFreq(&freqMinRad, mTotal);
  if (fMin > freqMinRad) {
    XLALPrintError("XLAL Error - %s: Intitial frequency is too high, the limit "
                   "is %4.10f \n",
                   __func__, freqMinRad);
    XLAL_ERROR(XLAL_EDOM);
  }
  /* Accuracies of adaptive Runge-Kutta integrator */
  /* Note that this accuracies are lower than those used in SEOBNRv2: they allow
   * reasonable runtimes for precessing systems */
  /* These accuracies can be adjusted according to desired accuracy and runtime
   */
  REAL8 EPS_ABS = 1.0e-8;
  REAL8 EPS_REL = 1.0e-8;
  /* When using adaptice steps in the Runge-Kutta integrator, minimal step dt/M
   * In units of mTotal, introduced because steps go to 0 in some symmetric,
   * opposite-spin configurations */
   REAL8 deltaT_min = 8.0e-5;

  /* Geometric output time step, in units of mTotal */
  REAL8 deltaT = INdeltaT / mTScaled;

  /* Step-back target length of 150M from end of AdaS integration for HiS
   * integration */
  REAL8 tStepBack = 150.;

  /* Flag to decide wether to include NQC corrections when computing the
   * waveform */
  UINT4 flagNQC = 0;
  /* Flag to decide wether to use adaptive sampling or constant sampling when
   * integrating the dynamics */
  UINT4 flagConstantSampling = 0;

  /* Cast, in order of appearance */
  REAL8Vector *ICvalues = NULL;
  REAL8Array *dynamicsAdaS = NULL;
  SEOBdynamics *seobdynamicsAdaS = NULL;
  REAL8Vector *seobvalues_tstartHiS = NULL;
  REAL8Vector *ICvaluesHiS = NULL;
  REAL8Array *dynamicsHiS = NULL;
  SEOBdynamics *seobdynamicsHiS = NULL;
  REAL8Vector *seobvalues_tPeakOmega = NULL;
  REAL8Vector *seobvalues_test = NULL;
  REAL8Vector *Jfinal = NULL;
  REAL8Vector *Lhatfinal = NULL;
  REAL8Vector *chi1L_tPeakOmega = NULL;
  REAL8Vector *chi2L_tPeakOmega = NULL;
  SphHarmListEOBNonQCCoeffs *nqcCoeffsList = NULL;
  SphHarmListCAmpPhaseSequence *listhPlm_HiS = NULL;
  SphHarmListCAmpPhaseSequence *listhPlm_HiSRDpatch = NULL;
  SphHarmListCAmpPhaseSequence *listhPlm_AdaS = NULL;
  *tVecPmodes = NULL;
  SEOBdynamics *seobdynamicsAdaSHiS = NULL;
  SphHarmListCAmpPhaseSequence *listhPlm = NULL;
  *alphaJ2P = NULL;
  *betaJ2P = NULL;
  *gammaJ2P = NULL;
  *hJlm = NULL;
  *hIlm = NULL;
  REAL8TimeSeries *hplusTS = NULL;
  REAL8TimeSeries *hcrossTS = NULL;
  *mergerParams = NULL;
  *seobdynamicsAdaSVector = NULL;
  *seobdynamicsHiSVector = NULL;
  *seobdynamicsAdaSHiSVector = NULL;

  /* Parameter structures */
  /* SpinEOBParams contains:
     EOBParams               *eobParams; // see below
     SpinEOBHCoeffs       *seobCoeffs; // see below
     EOBNonQCCoeffs     *nqcCoeffs; // see below
     REAL8Vector             *s1Vec; // spin1
     REAL8Vector             *s2Vec; // spin2
     REAL8Vector             *sigmaStar; // Eq. 5.3 of Barausse and Buonanno
                                        // PRD 81, 084024 (2010) [arXiv:0912.3517]
     REAL8Vector             *sigmaKerr; // Eq. 5.2 of Barausse and Buonanno
                                         // PRD 81, 084024 (2010) [arXiv:0912.3517]
     REAL8                        a; // |sigmaKerr|
     REAL8                        chi1; // spin1.LNhat/m1^2
     REAL8                        chi2; // spin2.LNhat/m2^2
     REAL8                        prev_dr; // stores
     dr/dt for stopping condition purposes
     int alignedSpins;   // flag to indicate whther the binary is precessing or not
     int tortoise; // flag to switch on/off tortoise coordinates when calling Hamiltonian
     int ignoreflux;  // flag to switch off radiation reaction when calling the ODEs via
     XLALSpinPrecHcapNumericalDerivative or as XLALSpinPrecHcapExactDerivative
      */
  SpinEOBParams seobParams;

  /* SpinEOBHCoeffs contains:
      double KK; // nonspinning calibration in Hamiltonian (for SEOBNRv2: 1.712
     − 1.804eta − 39:77eta^2 + 103.2eta^3)
     double k0; // Delta_i coefficients in the Delta_u potential Eq. 8 of
                // PRD 86, 024011 (2012) [arXiv:1202.0790] and https://dcc.ligo.org/T1400476
     double k1;
     double k2;
     double k3;
     double k4;
     double k5;
     double k5l;
     double b3; // omega_{fd}^1 frame dragging parameter in Hamiltonian
     (unused) double bb3; // omega_{fd}^2 frame dragging parameter in
     Hamiltonian (unused) double d1; // spin-orbit calibration in Hamiltonian
     for SEOBNRv1 double d1v2; // spin-orbit calibration in Hamiltonian for
     SEOBNRv1 double dheffSS; // spin-spin calibration in Hamiltonian for
     SEOBNRv1 double dheffSSv2; // spin-spin calibration in Hamiltonian for
     SEOBNRv1 UINT4    SpinAlignedEOBversion;
     int      updateHCoeffs;
  */
  SpinEOBHCoeffs seobCoeffs;

  /* Structure introduced for opt */
  SEOBHCoeffConstants seobCoeffConsts;

  /* EOBParams contains parameters common to nonspin and spin EOBNR models,
      including mass ratio, masses, pre-computed coefficients for potential,
     flux and waveforms, NQC coefficients and Newtonian multiple prefixes */
  EOBParams eobParams;

  /* Non-quasi-circular correction */
  EOBNonQCCoeffs nqcCoeffs;

  /* FacWaveformCoeffscontaining the coefficients for calculating the factorized
    waveform. The coefficients are precomputed in the function
    XLALSimIMREOBCalcSpinPrecFacWaveformCoefficients */
  FacWaveformCoeffs hCoeffs;

  /* NewtonMultipolePrefixes contains all the terms of the Newtonian multipole
     which are constant over the course of the evolution, and can therefore be
     pre-computed. They are stored in a two-dimensional array, which is
     indexed as values[l][m]. This is filled by
     XLALSimIMREOBComputeNewtonMultipolePrefixes */
  NewtonMultipolePrefixes prefixes;

  memset(&seobParams, 0, sizeof(seobParams));
  memset(&seobCoeffs, 0, sizeof(seobCoeffs));
  memset(&seobCoeffConsts, 0, sizeof(seobCoeffConsts));
  memset(&eobParams, 0, sizeof(eobParams));
  memset(&nqcCoeffs, 0, sizeof(nqcCoeffs));
  memset(&hCoeffs, 0, sizeof(hCoeffs));
  memset(&prefixes, 0, sizeof(prefixes));

  /* Recast the input ModeArray as a more convenient array of integers (l,m) */
  /* Also check that (2,2) is present in ModeArray */
  if (!(XLALSimInspiralModeArrayIsModeActive(modearray, 2, 2))) {
    FREE_ALL
    XLALPrintError("XLAL Error - %s: dominant harmonic (2,2) not present in "
                   "input ModeArray -- will be needed internally.\n",
                   __func__);
    PRINT_ALL_PARAMS
    XLAL_ERROR(XLAL_EINVAL);
  }
  UINT4 nmodes = SEOBGetNumberOfModesInModeArray(modearray, _SEOB_MODES_LMAX);
  INT4 modes_lmax = SEOBGetLMaxInModeArray(modearray, _SEOB_MODES_LMAX);
  INT4 modes[nmodes][2];
  memset(modes, 0, 2 * nmodes * sizeof(INT4));
  SEOBGetModeNumbersFromModeArray(modes, modearray, _SEOB_MODES_LMAX);
  if (debug) {
    printf("P-frame modes :\n");
    for (UINT4 nmode = 0; nmode < nmodes; nmode++)
      printf("(%d, %d)\n", modes[nmode][0], modes[nmode][1]);
  }

  /* Workspace spin vectors */
  REAL8Vector s1Vec, s2Vec, sigmaStar, sigmaKerr;
  REAL8 s1Vecdata[3] = {0.};
  REAL8 s2Vecdata[3] = {0.};
  REAL8 sigmaStardata[3] = {0.};
  REAL8 sigmaKerrdata[3] = {0.};
  s1Vec.length = s2Vec.length = sigmaStar.length = sigmaKerr.length = 3;
  s1Vec.data = s1Vecdata;
  s2Vec.data = s2Vecdata;
  sigmaStar.data = sigmaStardata;
  sigmaKerr.data = sigmaKerrdata;
  /* s1Vec, s2Vec spins in units of mTotal square */
  for (UINT4 j = 0; j < 3; j++) {
    s1Vec.data[j] = INchi1.data[j] * m1 * m1 / mTotal / mTotal;
    s2Vec.data[j] = INchi2.data[j] * m2 * m2 / mTotal / mTotal;
  }
  /* Compute sigmaStar and sigmaKerr */
  SEOBCalculateSigmaStar(&sigmaStar, m1, m2, &s1Vec, &s2Vec);
  SEOBCalculateSigmaKerr(&sigmaKerr, &s1Vec, &s2Vec);
  /* Calculate the value of a, that is magnitude of Eq. 31 in PRD 86, 024011 [arXiv:1202.0790]
   * (2012) */
  REAL8 a = sqrt(inner_product(sigmaKerr.data, sigmaKerr.data));

  /* Spin-EOB parameters */
  seobParams.a = a;
  seobParams.s1Vec = &s1Vec;
  seobParams.s2Vec = &s2Vec;
  seobParams.alignedSpins = 0;
  seobParams.tortoise = 1;
  seobParams.ignoreflux = 0;
  seobParams.sigmaStar = &sigmaStar;
  seobParams.sigmaKerr = &sigmaKerr;
  seobParams.seobCoeffs = &seobCoeffs;
  seobParams.seobCoeffConsts = &seobCoeffConsts;
  seobParams.eobParams = &eobParams;
  seobParams.nqcCoeffs = &nqcCoeffs;
  seobParams.seobApproximant = SEOBNRv4P;
  /* Here we store the initial LN-component of the dimensionless spins */
  seobParams.chi1 = chi1z;
  seobParams.chi2 = chi2z;
  /* Non-Spin-EOB parameters */
  eobParams.hCoeffs = &hCoeffs;
  eobParams.prefixes = &prefixes;
  seobCoeffs.SpinAlignedEOBversion = SpinAlignedEOBversion;
  eobParams.m1 = m1;
  eobParams.m2 = m2;
  eobParams.eta = eta;

  TidalEOBParams tidal1, tidal2;
  tidal1.mByM = m1SI / (m1SI + m2SI);
  tidal1.lambda2Tidal = 0.;
  tidal1.omega02Tidal = 0.;
  tidal1.lambda3Tidal = 0.;
  tidal1.omega03Tidal = 0.;
  tidal1.quadparam = 1.0;

  tidal2.mByM = m2SI / (m1SI + m2SI);
  tidal2.lambda2Tidal = 0.;
  tidal2.omega02Tidal = 0.;
  tidal2.lambda3Tidal = 0.;
  tidal2.omega03Tidal = 0.;
  tidal2.quadparam = 1.0;

  seobCoeffs.tidal1 = &tidal1;
  seobCoeffs.tidal2 = &tidal2;

  hCoeffs.tidal1 = &tidal1;
  hCoeffs.tidal2 = &tidal2;

  /* From opt */
  seobCoeffConsts = XLALEOBSpinPrecCalcSEOBHCoeffConstants(eta);

  REAL8 Lhat[3] = {0.0, 0.0, 1.0}; // Not quite true but should be very close
  REAL8 tempS1_p = inner_product(s1Vec.data, Lhat);
  REAL8 tempS2_p = inner_product(s2Vec.data, Lhat);
  REAL8 S1_perp[3] = {0, 0, 0};
  REAL8 S2_perp[3] = {0, 0, 0};
  for (UINT4 jj = 0; jj < 3; jj++) {
    S1_perp[jj] = s1Vec.data[jj] - tempS1_p * Lhat[jj];
    S2_perp[jj] = s2Vec.data[jj] - tempS2_p * Lhat[jj];
  }
  REAL8 sKerr_norm = sqrt(inner_product(sigmaKerr.data, sigmaKerr.data));
  REAL8 S_con = 0.0;
  if (sKerr_norm > 1e-6) {
    S_con = sigmaKerr.data[0] * Lhat[0] + sigmaKerr.data[1] * Lhat[1] +
            sigmaKerr.data[2] * Lhat[2];
    S_con /= (1 - 2 * eta);
    S_con += (inner_product(S1_perp, sigmaKerr.data) +
              inner_product(S2_perp, sigmaKerr.data)) /
             sKerr_norm / (1 - 2 * eta) / 2.;
  }

  /* Pre-compute the Hamiltonian coefficients */
  if (XLALSimIMRCalculateSpinPrecEOBHCoeffs_v2(
          &seobCoeffs, eta, a, S_con, SpinAlignedEOBversion) == XLAL_FAILURE) {
    FREE_ALL
    XLALPrintError(
        "XLAL Error: XLALSimIMRCalculateSpinPrecEOBHCoeffs_v2 failed.\n");
    PRINT_ALL_PARAMS
    XLAL_ERROR(XLAL_EFUNC);
  }

  /* Pre-compute the coefficients for the Newtonian factor of hLM
     Eq. A1 of PRD 86, 024011 (2012) [arXiv:1202.0790]  */
  if (XLALSimIMREOBComputeNewtonMultipolePrefixes(&prefixes, m1, m2) ==
      XLAL_FAILURE) {
    FREE_ALL
    XLALPrintError("XLAL Error - %s: "
                   "XLALSimIMREOBComputeNewtonMultipolePrefixes failed.\n",
                   __func__);
    PRINT_ALL_PARAMS
    XLAL_ERROR(XLAL_EFUNC);
  }

  /******************************************************************************************************************/
  /* STEP 1) Solve for initial conditions */
  /******************************************************************************************************************/

  if (debug)
    printf("STEP 1) Solve for initial conditions\n");

  flagConstantSampling = 0; // Adaptive sampling by default
  /* Gemetric starting frequency */
  REAL8 MfMin = mTScaled * fMin;
  /* If the in-plane components of both spins are small, fall back to aligned
   * spins for the IC/dynamics/Euler angles */
  UINT4 SpinsAlmostAligned = 0;
  if ((sqrt(chi1x * chi1x + chi1y * chi1y) < EPS_ALIGN) &&
      (sqrt(chi2x * chi2x + chi2y * chi2y) < EPS_ALIGN)) {
    SpinsAlmostAligned = 1;
    INchi1.data[0] = 0.;
    INchi1.data[1] = 0.;
    INchi2.data[0] = 0.;
    INchi2.data[1] = 0.;
    // In the aligned-spin limit we now do exactly as v4
    EPS_REL = 1.0e-9;
    EPS_ABS = 1.0e-10;
    flagConstantSampling = 1;
  } else {
    SpinsAlmostAligned = 0;
  }
  seobParams.alignedSpins = SpinsAlmostAligned;

  /* Compute vector for initial conditions at fMin */
  if (SEOBInitialConditions(&ICvalues, MfMin, m1, m2, &INchi1, &INchi2,
                            &seobParams,
                            flagHamiltonianDerivative) == XLAL_FAILURE) {
    FREE_ALL
    XLALPrintError("XLAL Error - %s: SEOBInitialConditions failed.\n",
                   __func__);
    PRINT_ALL_PARAMS
    XLAL_ERROR(XLAL_EFUNC);
  }

  /******************************************************************************************************************/
  /* STEP 2) Evolve EOB trajectory with adaptive sampling (AdaS) */
  /******************************************************************************************************************/

  if (debug)
    printf("STEP 2) Evolve EOB trajectory with adaptive sampling (AdaS)\n");

  UINT4 retLenAdaS = 0;
  REAL8 tendAdaS =
      20. / mTScaled;    /* This is 20s in geometric units, seems it should be
                            ignored anyway because of integrator->stopontestonly */
  REAL8 tstartAdaS = 0.; /* t=0 will set at the starting time */
  /* Note: the timesampling step deltaT is used internally only to initialize
   * adaptive step */
  if (SEOBIntegrateDynamics(&dynamicsAdaS, &retLenAdaS, ICvalues, EPS_ABS,
                            EPS_REL, deltaT, deltaT_min, tstartAdaS, tendAdaS,
                            &seobParams, flagConstantSampling,
                            flagHamiltonianDerivative) == XLAL_FAILURE) {
    FREE_ALL
    XLALPrintError(
        "XLAL Error - %s: SEOBIntegrateDynamicsAdaptiveSampling failed.\n",
        __func__);
    PRINT_ALL_PARAMS
    XLAL_ERROR(XLAL_EFUNC);
  }

  /* Compute derived quantities from the dynamics */
  if (SEOBComputeExtendedSEOBdynamics(
          &seobdynamicsAdaS, dynamicsAdaS, retLenAdaS, &seobParams,
          flagHamiltonianDerivative, flagZframe) == XLAL_FAILURE) {
    FREE_ALL
    XLALPrintError("XLAL Error - %s: SEOBComputeExtendedSEOBdynamics failed.\n",
                   __func__);
    PRINT_ALL_PARAMS
    XLAL_ERROR(XLAL_EFUNC);
  }

  /******************************************************************************************************************/
  /* STEP 3) Step back and evolve EOB trajectory at high sampling rate (HiS) */
  /******************************************************************************************************************/

  if (debug)
    printf("STEP 3) Step back and evolve EOB trajectory at high sampling rate "
           "(HiS)\n");

  /* Error tolerances */
  if (!SpinsAlmostAligned){
    // For aligned spins we use the same tolerance for high-sampling
    // rate and for low sampling rate.
    EPS_ABS = 1e-8;
    EPS_REL = 1e-8;
  }

  /* Time step for high-sampling part */
  REAL8 deltaTHiS = 1. / 50; /* Fixed at 1/50M */

  /* Step-back time */
  /* Stepback by at least 150M wrt the ending time of the adaptive-sampling-rate
   * trajectory */
  // NOTE: this stepping back is discrete, and can jump by one sample of AdaS
  // when changing continuously parameters
  REAL8 tstartHiSTarget = seobdynamicsAdaS->tVec[retLenAdaS - 1] - tStepBack;
  INT4 indexstartHiS = retLenAdaS - 1; /* index for the AdaS dynamics */
  while ((indexstartHiS > 0) &&
         (seobdynamicsAdaS->tVec[indexstartHiS] > tstartHiSTarget))
    indexstartHiS--;
  REAL8 tstartHiS = seobdynamicsAdaS->tVec[indexstartHiS];

  /* Compute the new initial conditions for high-sampling integration by
   * interpolating adaptive-sampling dynamics */
  /* Since we take values at a sample, no interpolation actually needed - but
   * this is flexible if we change the prescription for stepback */
  if (SEOBInterpolateDynamicsAtTime(&seobvalues_tstartHiS, tstartHiS,
                                    seobdynamicsAdaS) == XLAL_FAILURE) {
    FREE_ALL
    XLALPrintError("XLAL Error - %s: SEOBInterpolateDynamicsAtTime failed.\n",
                   __func__);
    PRINT_ALL_PARAMS
    XLAL_ERROR(XLAL_EFUNC);
  }
  if (!(ICvaluesHiS = XLALCreateREAL8Vector(14))) {
    XLALPrintError(
        "XLAL Error - %s: failed to allocate REAL8Vector ICvaluesHiS.\n",
        __func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  memcpy(ICvaluesHiS->data, &(seobvalues_tstartHiS->data[1]),
         14 * sizeof(REAL8));

  /* Integrate again the dynamics with a constant high sampling rate */
  seobParams.prev_dr = 0.; // This is used to check whether dr/dt is increasing
                           // in the stopping condition
  seobParams.termination_reason =
      -999; // This is going to store the termination condition for the
            // high-sampling integration
  UINT4 retLenHiS = 0;
  REAL8 tendHiS =
      tstartHiS + tStepBack; /* Seems it should be ignored anyway because of
                                integrator->stopontestonly */
  flagConstantSampling = 1;
  /* Note: here deltaT_min = 0. will simply be ignored, we use fixed steps */
  if (SEOBIntegrateDynamics(&dynamicsHiS, &retLenHiS, ICvaluesHiS, EPS_ABS,
                            EPS_REL, deltaTHiS, 0., tstartHiS, tendHiS,
                            &seobParams, flagConstantSampling,
                            flagHamiltonianDerivative) == XLAL_FAILURE) {
    FREE_ALL
    XLALPrintError(
        "XLAL Error - %s: SEOBIntegrateDynamicsConstantSampling failed.\n",
        __func__);
    PRINT_ALL_PARAMS
    XLAL_ERROR(XLAL_EFUNC);
  }

  /* Compute derived quantities for the high-sampling dynamics */
  if (SEOBComputeExtendedSEOBdynamics(&seobdynamicsHiS, dynamicsHiS, retLenHiS,
                                      &seobParams, flagHamiltonianDerivative,
                                      flagZframe) == XLAL_FAILURE) {
    FREE_ALL
    XLALPrintError("XLAL Error - %s: SEOBComputeExtendedSEOBdynamics failed.\n",
                   __func__);
    PRINT_ALL_PARAMS
    XLAL_ERROR(XLAL_EFUNC);
  }
  /* Find time of peak of omega for the High-sampling dynamics */
  INT4 foundPeakOmega = 0;
  REAL8 tPeakOmega = 0.;
  SEOBLocateTimePeakOmega(&tPeakOmega, &foundPeakOmega, dynamicsHiS,
                          seobdynamicsHiS, retLenHiS, &seobParams,
                          flagHamiltonianDerivative);

  /******************************************************************************************************************/
  /* STEP 4) Get final J/L/spins from HiS dynamics at peak of Omega, compute
   * constant angles EulerI2J */
  /******************************************************************************************************************/

  if (debug)
    printf("STEP 4) Get final J/L/spins from HiS dynamics at peak of Omega, "
           "compute constant angles EulerI2J\n");

  /* Compute final dynamics quantities, interpolating HiS dynamics at tPeakOmega
   */

  if (SEOBInterpolateDynamicsAtTime(&seobvalues_tPeakOmega, tPeakOmega,
                                    seobdynamicsHiS) == XLAL_FAILURE) {
    FREE_ALL
    XLALPrintError("XLAL Error - %s: SEOBInterpolateDynamicsAtTime failed at "
                   "tPeakOmega.\n",
                   __func__);
    XLAL_ERROR(XLAL_EFUNC);
  }

  /* Interpolate the dynamics to r=10M. This is used as a
  feducial point to evaluate the final mass and spin fits */
  REAL8Vector timeVec;
  timeVec.length = retLenAdaS;
  timeVec.data = seobdynamicsAdaS->tVec;
  REAL8Vector* rVec = XLALCreateREAL8Vector(retLenAdaS);
  for (UINT4 jj = 0; jj < retLenAdaS; jj++) {
    rVec->data[jj] = -1* seobdynamicsAdaS->polarrVec[jj];
  }
  /*
  // Uncomment this to restore behaviour where we used the Keplerian frequency
  corresponding
  // to r=10M instead of r itself.
  UNUSED REAL8Vector omegaVec;
  omegaVec.length = retLenAdaS;
  omegaVec.data = seobdynamicsAdaS->omegaVec;
  UNUSED REAL8 omega_10M = pow(10.0, -1.5); // Keplerian value of omega at r=10M
  UINT4 index_10M = FindClosestIndex(&omegaVec, omega_10M);
  */

  // FindClosestIndex requires an *increasing* function of time, so we use -r
  // instead of r
  UINT4 index_10M = FindClosestIndex(rVec, -10.0);
  REAL8 time_10M = timeVec.data[index_10M];
  if (SEOBInterpolateDynamicsAtTime(&seobvalues_test, time_10M,
                                    seobdynamicsAdaS) == XLAL_FAILURE) {
    FREE_ALL
    XLALPrintError(
        "XLAL Error - %s: SEOBInterpolateDynamicsAtTime failed at time_10M.\n",
        __func__);
    XLAL_ERROR(XLAL_EFUNC);
  }

  /* Compute the timeshift to get the attachment point */
  SEOBLFrameVectors(&chi1L_tPeakOmega, &chi2L_tPeakOmega, seobvalues_tPeakOmega,
                    m1, m2, flagZframe);
  REAL8 Deltat22 = XLALSimIMREOBGetNRSpinPeakDeltaTv4(
      2, 2, m1, m2, chi1L_tPeakOmega->data[2], chi2L_tPeakOmega->data[2]);

  /* Determine the time of attachment */
  REAL8 tAttach = tPeakOmega - Deltat22;

  /* Compute final J from dynamics quantities */
  SEOBJfromDynamics(&Jfinal, seobvalues_tPeakOmega, &seobParams);
  /*Compute the L-hat vector. Note that it has unit norm */
  SEOBLhatfromDynamics(&Lhatfinal, seobvalues_tPeakOmega, &seobParams, flagZframe);
  REAL8 Jmag = sqrt(inner_product(Jfinal->data, Jfinal->data));
  /* Cosine of the angle between L-hat and J. Needed to determine
   * the correct sign of the final spin
   */
  REAL8 cos_angle = inner_product(Jfinal->data, Lhatfinal->data) / Jmag;
  /* Compute final-J-frame unit vectors e1J, e2J, e3J=Jfinalhat */
  /* Convention: if (ex, ey, ez) is the initial I-frame, e1J chosen such that ex
   * is in the plane (e1J, e3J) and ex.e1J>0 */
  REAL8Vector e1J, e2J, e3J;
  e1J.length = e2J.length = e3J.length = 3;
  REAL8 e1Jdata[3] = {0.};
  REAL8 e2Jdata[3] = {0.};
  REAL8 e3Jdata[3] = {0.};
  e1J.data = e1Jdata;
  e2J.data = e2Jdata;
  e3J.data = e3Jdata;
  SEOBBuildJframeVectors(&e1J, &e2J, &e3J, Jfinal);

  /* Compute Euler angles from initial I-frame to final-J-frame */
  /* Note: if spins are aligned, the function SEOBEulerI2JFromJframeVectors */
  /* becomes ill-defined - just keep these Euler angles to zero then */
  REAL8 alphaI2J = 0., betaI2J = 0., gammaI2J = 0.;
  if (!SpinsAlmostAligned) {
    SEOBEulerI2JFromJframeVectors(&alphaI2J, &betaI2J, &gammaI2J, &e1J, &e2J,
                                  &e3J);
  }

  /******************************************************************************************************************/
  /* STEP 5) Compute P-frame amp/phase for all modes on HiS and compute NQC
   * coefficients */
  /******************************************************************************************************************/

  if (debug)
    printf("STEP 5) Compute P-frame of modes amp/phase on HiS and compute NQC "
           "coefficients\n");

  if (SEOBCalculateSphHarmListNQCCoefficientsV4(
          &nqcCoeffsList, modes, nmodes, tPeakOmega, seobdynamicsHiS,
          &seobParams, chi1L_tPeakOmega, chi2L_tPeakOmega) == XLAL_FAILURE) {
    FREE_ALL
    XLALPrintError("XLAL Error - %s: NQC computation failed.\n", __func__);
    PRINT_ALL_PARAMS
    XLAL_ERROR(XLAL_EFUNC);
  }

  /******************************************************************************************************************/
  /* STEP 6) Compute P-frame amp/phase for all modes on HiS, now including NQC
   */
  /******************************************************************************************************************/

  if (debug)
    printf("STEP 6) Compute P-frame amp/phase on HiS, now including NQC\n");

  /* We now include the NQC when generating modes */
  flagNQC = 1;

  /* Compute amplitude and phase of the P-frame modes hPlm on high sampling,
   * with NQC */
  SEOBCalculateSphHarmListhlmAmpPhase(&listhPlm_HiS, modes, nmodes,
                                      seobdynamicsHiS, nqcCoeffsList,
                                      &seobParams, flagNQC);

  /******************************************************************************************************************/
  /* STEP 7) Attach RD to the P-frame waveform */
  /******************************************************************************************************************/

  if (debug)
    printf("STEP 7) Attach RD to the P-frame waveform\n");

  /* Compute the final mass and spins here and pass them on */
  REAL8 finalMass = 0., finalSpin = 0.;
  if (SEOBGetFinalSpinMass(&finalMass, &finalSpin, seobvalues_test, &seobParams,
                           flagZframe)) {
    FREE_ALL
    XLALPrintError("XLAL Error - %s: SEOBGetFinalSpinMass failed.\n", __func__);
    PRINT_ALL_PARAMS
    XLAL_ERROR(XLAL_EFUNC);
  }

  /* The function above returns only the magnitude of the spin.
   *  We pick the direction based on whether Lhat \cdot J is positive
   *  or negative */
  if (cos_angle < 0)
  {
    finalSpin *= -1;
  }
  /* finalSpin interpolation is available only between -0.9996 and 0.9996 */
  /* Set finalSpin to +/- 0.9996 if it is out of this range */
  if (finalSpin < -0.9996)
    finalSpin = -0.9996;
  if (finalSpin > 0.9996)
    finalSpin = 0.9996;

  if (debug) {
    XLAL_PRINT_INFO("final mass = %e, final spin = %e\n", finalMass, finalSpin);
  }

  /* Estimate leading QNM , to determine how long the ringdown
   * patch will be */
  /* NOTE: XLALSimIMREOBGenerateQNMFreqV2Prec returns the complex frequency in
   * physical units... */
  COMPLEX16Vector sigmaQNM220estimatephysicalVec;
  COMPLEX16 sigmaQNM220estimatephysical = 0.;
  sigmaQNM220estimatephysicalVec.length = 1;
  sigmaQNM220estimatephysicalVec.data = &sigmaQNM220estimatephysical;
  if (XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec(&sigmaQNM220estimatephysicalVec, m1,
                                                  m2, finalMass, finalSpin, 2,
                                                  2, 1) == XLAL_FAILURE) {
    FREE_ALL
    XLALPrintError(
        "XLAL Error - %s: failure in "
        "XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec for mode (l,m) = (2,2).\n",
        __func__);
    XLAL_ERROR(XLAL_EFUNC);
  }
  COMPLEX16 sigmaQNM220estimate = mTScaled * sigmaQNM220estimatephysical;

  /* Length of RD patch, 40 e-folds of decay of the estimated QNM220 */
  UINT4 retLenRDPatch =
      (UINT4)ceil(40 / (cimag(sigmaQNM220estimate) * deltaTHiS));


  /* Attach RD to the P-frame modes */
  /* Vector holding the values of the 0-th overtone QNM complex frequencies for
   * the modes (l,m) */
  COMPLEX16Vector *sigmaQNMlm0 = NULL;
  // NOTE: the QNM complex frequencies are computed inside
  // SEOBAttachRDToSphHarmListhPlm
  SEOBAttachRDToSphHarmListhPlm(
      &listhPlm_HiSRDpatch, &sigmaQNMlm0, modes, nmodes, finalMass, finalSpin,
      listhPlm_HiS, deltaTHiS, retLenHiS, retLenRDPatch, tAttach,
      seobvalues_tPeakOmega, seobdynamicsHiS, &seobParams, flagZframe, debug);

  /******************************************************************************************************************/
  /* STEP 8) Build the joined dynamics AdaS+HiS up to attachment, joined P-modes
   * AdaS+HiS+RDpatch */
  /******************************************************************************************************************/

  if (debug)
    printf("STEP 8) Build the joined dynamics AdaS+HiS up to attachment, "
           "joined P-modes AdaS+HiS+RDpatch\n");

  /* Compute amplitude and phase of the P-frame modes hPlm on adaptive sampling,
   * with NQC */
  flagNQC = 1;
  SEOBCalculateSphHarmListhlmAmpPhase(&listhPlm_AdaS, modes, nmodes,
                                      seobdynamicsAdaS, nqcCoeffsList,
                                      &seobParams, flagNQC);

  /* Vector of times for the P-modes, that will be used for interpolation:
   * joining AdaS and HiS+RDpatch */
  UINT4 retLenPmodes = 0;
  /* First junction at indexAdaSHiS, tAdaSHiS */
  UINT4 indexJoinHiS = 0;
  REAL8 tJoinHiS = 0.;
  /* Second junction at seobdynamicsAdaS, tJoinAttach */
  UINT4 indexJoinAttach = 0;
  REAL8 tJoinAttach = 0.;
  /* Construct the joined vector of times (AdaS+HiS+RDpatch) and keep the
   * jonction indices and times */
  SEOBJoinTimeVector(tVecPmodes, &retLenPmodes, &tJoinHiS, &indexJoinHiS,
                     &tJoinAttach, &indexJoinAttach, retLenRDPatch, deltaTHiS,
                     tstartHiS, tAttach, seobdynamicsAdaS, seobdynamicsHiS);

  /* Copy dynamics from AdaS<HiS and HiS<tAttach to form joined dynamics, ending
   * at the last time sample <tAttach */
  // NOTE: we cut the dynamics at tAttach, as we will extend the Euler
  // angles for t>=tAttach -- but we could also choose to finish at tPeakOmega
  // which is used for final-J and for the final mass/spin fit
  SEOBJoinDynamics(&seobdynamicsAdaSHiS, seobdynamicsAdaS, seobdynamicsHiS,
                   indexJoinHiS, indexJoinAttach);

  /* Copy waveform modes from AdaS and HiS+RDpatch - adjusting 2pi-phase shift
   * at the junction point AdaS/HiS */
  SEOBJoinSphHarmListhlm(&listhPlm, listhPlm_AdaS, listhPlm_HiSRDpatch, modes,
                         nmodes, indexstartHiS);

  /* Get the time of the frame-invariant amplitude peak */
  REAL8 tPeak = 0;
  UINT4 indexPeak = 0;
  // NOTE: peak amplitude using l=2 only: h22 required, and h21 used if present
  SEOBAmplitudePeakFromAmp22Amp21(&tPeak, &indexPeak, listhPlm, modes, nmodes,
                                  *tVecPmodes);

  /******************************************************************************************************************/
  /* STEP 9) Compute Euler angles J2P from AdaS and HiS dynamics up to
   * attachment */
  /******************************************************************************************************************/

  if (debug)
    printf("STEP 9) Compute Euler angles J2P from AdaS and HiS dynamics up to "
           "attachment\n");

  /* Compute Euler angles J2P from the dynamics before attachment point */
  /* If SpinsAlmostAligned, all Euler angles are set to 0 */
  if (SEOBEulerJ2PFromDynamics(alphaJ2P, betaJ2P, gammaJ2P, &e1J, &e2J, &e3J,
                               retLenPmodes, indexJoinAttach,
                               seobdynamicsAdaSHiS, &seobParams,
                               flagZframe) == XLAL_FAILURE) {
    FREE_ALL
    XLALPrintError("XLAL Error - %s: SEOBEulerJ2PFromDynamics failed.\n",
                   __func__);
    PRINT_ALL_PARAMS
    XLAL_ERROR(XLAL_EFUNC);
  }

  /******************************************************************************************************************/
  /* STEP 10) Compute Euler angles J2P extension after attachment */
  /******************************************************************************************************************/

  if (debug)
    printf("STEP 10) Compute Euler angles J2P extension after attachment\n");

  /* Compute Euler angles J2P according to the prescription flagEulerextension
   * after attachment point */
  /* If SpinsAlmostAligned, all Euler angles are set to 0 */
  /* NOTE: Regardless of the mode content of hPlm, the frame extension at the
   * moment is based on sigmaQNM22, sigmaQNM21 */
  COMPLEX16 sigmaQNM220 = 0., sigmaQNM210 = 0.;
  COMPLEX16Vector sigmaQNM220physicalVec, sigmaQNM210physicalVec;
  sigmaQNM220physicalVec.length = 1;
  sigmaQNM210physicalVec.length = 1;
  COMPLEX16 sigmaQNM220physicalval = 0.;
  COMPLEX16 sigmaQNM210physicalval = 0.;
  sigmaQNM220physicalVec.data = &sigmaQNM220physicalval;
  sigmaQNM210physicalVec.data = &sigmaQNM210physicalval;
  /* NOTE: XLALSimIMREOBGenerateQNMFreqV2Prec returns the complex frequency in
   * physical units... */
  if (XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec(&sigmaQNM220physicalVec, m1,
                                                  m2, finalMass, finalSpin, 2,
                                                  2, 1) == XLAL_FAILURE) {
    FREE_ALL
    XLALPrintError(
        "XLAL Error - %s: failure in "
        "XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec for mode (l,m) = (2,2).\n",
        __func__);
    XLAL_ERROR(XLAL_EFUNC);
  }
  if (XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec(&sigmaQNM210physicalVec, m1,
                                                  m2, finalMass, finalSpin, 2,
                                                  1, 1) == XLAL_FAILURE) {
    FREE_ALL
    XLALPrintError(
        "XLAL Error - %s: failure in "
        "XLALSimIMREOBGenerateQNMFreqV2FromFinalPrec for mode (l,m) = (2,1).\n",
        __func__);
    XLAL_ERROR(XLAL_EFUNC);
  }
  sigmaQNM220 = mTScaled * sigmaQNM220physicalVec.data[0];
  sigmaQNM210 = mTScaled * sigmaQNM210physicalVec.data[0];
  INT4 flip = 1;
  if (cos_angle < 0){
    flip = -1;
  }
  SEOBEulerJ2PPostMergerExtension(
      *alphaJ2P, *betaJ2P, *gammaJ2P, sigmaQNM220, sigmaQNM210, *tVecPmodes,
      retLenPmodes, indexJoinAttach, &seobParams, flagEulerextension, flip);

  /******************************************************************************************************************/
  /* STEP 11) Compute modes hJlm on the output time series by rotating and
   * interpolating the modes hPlm  */
  /******************************************************************************************************************/

  if (debug)
    printf("STEP 11) Compute modes hJlm on the output time series by rotating "
           "and interpolating the modes hPlm\n");

  /* Determine the length of the fixed-sampling output time series */
  UINT4 retLenTS =
      floor(((*tVecPmodes)->data[retLenPmodes - 1] - (*tVecPmodes)->data[0]) /
            deltaT);

  /* Rotate waveform from P-frame to J-frame */
  if (SEOBRotateInterpolatehJlmReImFromSphHarmListhPlmAmpPhase(
          hJlm, modes, nmodes, modes_lmax, deltaT, retLenTS, *tVecPmodes,
          listhPlm, *alphaJ2P, *betaJ2P, *gammaJ2P,
          flagSymmetrizehPlminusm) == XLAL_FAILURE) {
    FREE_ALL
    XLALPrintError(
        "XLAL Error - %s: failure in "
        "SEOBRotateInterpolatehJlmReImFromSphHarmListhPlmAmpPhase.\n",
        __func__);
    XLAL_ERROR(XLAL_EFUNC);
  }

  /******************************************************************************************************************/
  /* STEP 12) Rotate waveform from J-frame to the output I-frame on
   * timeseries-sampling (constant Wigner coeffs) */
  /******************************************************************************************************************/

  if (debug)
    printf("STEP 12) Rotate waveform from J-frame to the output I-frame on "
           "timeseries-sampling (constant Wigner coeffs)\n");

  /* Rotate waveform from J-frame to I-frame */
  SEOBRotatehIlmFromhJlm(hIlm, *hJlm, modes_lmax, alphaI2J, betaI2J, gammaI2J,
                         deltaT);

  /******************************************************************************************************************/
  /* STEP 13) Compute hplus, hcross from I-frame waveform on timeseries sampling
   */
  /******************************************************************************************************************/

  if (debug)
    printf("STEP 13) Compute hplus, hcross from I-frame waveform on timeseries "
           "sampling\n");

  /* GPS time for output time series and modes */
  LIGOTimeGPS tGPS = LIGOTIMEGPSZERO;
  XLALGPSAdd(&tGPS, -mTScaled *
                        tPeak); /* tPeak converted back to dimensionfull time */

  /* Create output timeseries for hplus, hcross */
  /* Use the dimensionfull INdeltaT (s) as time step */
  hplusTS = XLALCreateREAL8TimeSeries("H_PLUS", &tGPS, 0.0, INdeltaT,
                                      &lalStrainUnit, retLenTS);
  hcrossTS = XLALCreateREAL8TimeSeries("H_CROSS", &tGPS, 0.0, INdeltaT,
                                       &lalStrainUnit, retLenTS);

  /* Compute hplus, hcross from hIlm */
  // NOTE: azimuthal angle of the observer entering the -2Ylm is pi/2-phi
  // according to LAL conventions
  SEOBComputehplushcrossFromhIlm(hplusTS, hcrossTS, modes_lmax, *hIlm, amp0,
                                 inc, phi);

  /******************************************************************************************************************/
  /* STEP -1) Output and cleanup */
  /******************************************************************************************************************/

  if (debug)
    printf("STEP -1) Output and cleanup\n");

  /* Output vector gathering quantities related to merger (similar to previous
   * AttachParams) */
  /* Format: tPeakOmega tAttach tPeak Jfinalx Jfinaly Jfinalz finalMassfit
   * finalSpinfit termination_reason [sigmaQNMlm0Re sigmaQNMlm0Im for lm in
   * modes] */
  /* NOTE: the size of this output vector depends on the number of modes, due to
   * the sigmaQNM */
  if (!((*mergerParams) = XLALCreateREAL8Vector(9 + 2 * nmodes))) {
    XLALPrintError(
        "XLAL Error - %s: failed to allocate REAL8Vector mergerParams.\n",
        __func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  (*mergerParams)->data[0] = tPeakOmega;
  (*mergerParams)->data[1] = tAttach;
  (*mergerParams)->data[2] = tPeak;
  (*mergerParams)->data[3] = Jfinal->data[0];
  (*mergerParams)->data[4] = Jfinal->data[1];
  (*mergerParams)->data[5] = Jfinal->data[2];
  (*mergerParams)->data[6] = finalMass;
  (*mergerParams)->data[7] = finalSpin;
  (*mergerParams)->data[8] = seobParams.termination_reason;
  for (UINT4 nmode = 0; nmode < nmodes; nmode++) {
    (*mergerParams)->data[9 + 2 * nmode] = creal(sigmaQNMlm0->data[nmode]);
    (*mergerParams)->data[9 + 2 * nmode + 1] = cimag(sigmaQNMlm0->data[nmode]);
  }

  /* Point the output pointers to the relevant time series and return */
  (*hplus) = hplusTS;
  (*hcross) = hcrossTS;

  /* Additional outputs */

  /* Dynamics */
  // NOTE: casting to REAL8Vector due to the SWIG wrapping
  *seobdynamicsAdaSVector =
      XLALCreateREAL8Vector(v4PdynamicsVariables * seobdynamicsAdaS->length);
  *seobdynamicsHiSVector =
      XLALCreateREAL8Vector(v4PdynamicsVariables * seobdynamicsHiS->length);
  *seobdynamicsAdaSHiSVector =
      XLALCreateREAL8Vector(v4PdynamicsVariables * seobdynamicsAdaSHiS->length);
  memcpy((*seobdynamicsAdaSVector)->data, seobdynamicsAdaS->array->data,
         (v4PdynamicsVariables * seobdynamicsAdaS->length) * sizeof(REAL8));
  memcpy((*seobdynamicsHiSVector)->data, seobdynamicsHiS->array->data,
         (v4PdynamicsVariables * seobdynamicsHiS->length) * sizeof(REAL8));
  memcpy((*seobdynamicsAdaSHiSVector)->data, seobdynamicsAdaSHiS->array->data,
         (v4PdynamicsVariables * seobdynamicsAdaSHiS->length) * sizeof(REAL8));

  /* Modes in the P-frame */
  // NOTE: casting to REAL8Vector due to the SWIG wrapping
  // NOTE: in the output, real amplitude instead of complex envelope
  /* (2,2) */
  *hP22_amp = XLALCreateREAL8Vector(retLenPmodes);
  *hP22_phase = XLALCreateREAL8Vector(retLenPmodes);
  SphHarmListCAmpPhaseSequence *hP22 =
      SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 2, 2);
  if (hP22 == NULL) {
    memset((*hP22_amp)->data, 0, retLenPmodes * sizeof(REAL8));
    memset((*hP22_phase)->data, 0, retLenPmodes * sizeof(REAL8));
  } else {
    memcpy((*hP22_amp)->data, hP22->campphase->camp_real->data,
           retLenPmodes * sizeof(REAL8));
    memcpy((*hP22_phase)->data, hP22->campphase->phase->data,
           retLenPmodes * sizeof(REAL8));
  }
  /* (2,1) */
  *hP21_amp = XLALCreateREAL8Vector(retLenPmodes);
  *hP21_phase = XLALCreateREAL8Vector(retLenPmodes);
  SphHarmListCAmpPhaseSequence *hP21 =
      SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 2, 1);
  if (hP21 == NULL) {
    memset((*hP21_amp)->data, 0, retLenPmodes * sizeof(REAL8));
    memset((*hP21_phase)->data, 0, retLenPmodes * sizeof(REAL8));
  } else {
    memcpy((*hP21_amp)->data, hP21->campphase->camp_real->data,
           retLenPmodes * sizeof(REAL8));
    memcpy((*hP21_phase)->data, hP21->campphase->phase->data,
           retLenPmodes * sizeof(REAL8));
  }
  /* (3,3) */
  *hP33_amp = XLALCreateREAL8Vector(retLenPmodes);
  *hP33_phase = XLALCreateREAL8Vector(retLenPmodes);
  SphHarmListCAmpPhaseSequence *hP33 =
      SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 3, 3);
  if (hP33 == NULL) {
    memset((*hP33_amp)->data, 0, retLenPmodes * sizeof(REAL8));
    memset((*hP33_phase)->data, 0, retLenPmodes * sizeof(REAL8));
  } else {
    memcpy((*hP33_amp)->data, hP33->campphase->camp_real->data,
           retLenPmodes * sizeof(REAL8));
    memcpy((*hP33_phase)->data, hP33->campphase->phase->data,
           retLenPmodes * sizeof(REAL8));
  }
  /* (4,4) */
  *hP44_amp = XLALCreateREAL8Vector(retLenPmodes);
  *hP44_phase = XLALCreateREAL8Vector(retLenPmodes);
  SphHarmListCAmpPhaseSequence *hP44 =
      SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 4, 4);
  if (hP44 == NULL) {
    memset((*hP44_amp)->data, 0, retLenPmodes * sizeof(REAL8));
    memset((*hP44_phase)->data, 0, retLenPmodes * sizeof(REAL8));
  } else {
    memcpy((*hP44_amp)->data, hP44->campphase->camp_real->data,
           retLenPmodes * sizeof(REAL8));
    memcpy((*hP44_phase)->data, hP44->campphase->phase->data,
           retLenPmodes * sizeof(REAL8));
  }
  /* (5,5) */
  *hP55_amp = XLALCreateREAL8Vector(retLenPmodes);
  *hP55_phase = XLALCreateREAL8Vector(retLenPmodes);
  SphHarmListCAmpPhaseSequence *hP55 =
      SphHarmListCAmpPhaseSequence_GetMode(listhPlm, 5, 5);
  if (hP55 == NULL) {
    memset((*hP55_amp)->data, 0, retLenPmodes * sizeof(REAL8));
    memset((*hP55_phase)->data, 0, retLenPmodes * sizeof(REAL8));
  } else {
    memcpy((*hP55_amp)->data, hP55->campphase->camp_real->data,
           retLenPmodes * sizeof(REAL8));
    memcpy((*hP55_phase)->data, hP55->campphase->phase->data,
           retLenPmodes * sizeof(REAL8));
  }

  /* Cleanup */
  if (ICvalues != NULL)
    XLALDestroyREAL8Vector(ICvalues);
  if (dynamicsAdaS != NULL)
    XLALDestroyREAL8Array(dynamicsAdaS);
  if (rVec != NULL)
    XLALDestroyREAL8Vector(rVec);
  if (seobdynamicsAdaS != NULL)
    SEOBdynamics_Destroy(seobdynamicsAdaS);
  if (seobvalues_tstartHiS != NULL)
    XLALDestroyREAL8Vector(seobvalues_tstartHiS);
  if (ICvaluesHiS != NULL)
    XLALDestroyREAL8Vector(ICvaluesHiS);
  if (dynamicsHiS != NULL)
    XLALDestroyREAL8Array(dynamicsHiS);
  if (seobdynamicsHiS != NULL)
    SEOBdynamics_Destroy(seobdynamicsHiS);
  if (seobvalues_tPeakOmega != NULL)
    XLALDestroyREAL8Vector(seobvalues_tPeakOmega);
  if (seobvalues_test != NULL)
    XLALDestroyREAL8Vector(seobvalues_test);
  if (Jfinal != NULL)
    XLALDestroyREAL8Vector(Jfinal);
  if (Lhatfinal != NULL)
    XLALDestroyREAL8Vector(Lhatfinal);
  if (nqcCoeffsList != NULL)
    SphHarmListEOBNonQCCoeffs_Destroy(nqcCoeffsList);
  if (listhPlm_HiS != NULL)
    SphHarmListCAmpPhaseSequence_Destroy(listhPlm_HiS);
  if (listhPlm_HiSRDpatch != NULL)
    SphHarmListCAmpPhaseSequence_Destroy(listhPlm_HiSRDpatch);
  if (listhPlm_AdaS != NULL)
    SphHarmListCAmpPhaseSequence_Destroy(listhPlm_AdaS);
  if (seobdynamicsAdaSHiS != NULL)
    SEOBdynamics_Destroy(seobdynamicsAdaSHiS);
  if (listhPlm != NULL)
    SphHarmListCAmpPhaseSequence_Destroy(listhPlm);
  if (chi1L_tPeakOmega != NULL)
    XLALDestroyREAL8Vector(chi1L_tPeakOmega);
  if (chi2L_tPeakOmega != NULL)
    XLALDestroyREAL8Vector(chi2L_tPeakOmega);
  if (sigmaQNMlm0 != NULL)
    XLALDestroyCOMPLEX16Vector(sigmaQNMlm0);
  return XLAL_SUCCESS;
}
#undef FREE_ALL
#undef PRINT_ALL_PARAMS

#endif
