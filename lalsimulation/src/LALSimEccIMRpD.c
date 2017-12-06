/*
 *  Implemented by: Maria Haney (2017)
 *
 *  VARIANT OF IMRPHENOMD APPROXIMANT
 *  INCLUDES LEADING-ORDER ECCENTRIC CONTRIBUTIONS UP TO 3PN PHASE ORDER
 *  IN THE INSPIRAL
 *  DEFINED IN LALSimInspiralEccentricPhasingCoeffs.c
 *  BASED ON arXiv:1602.0308
 *
 *  Note: For this `effective` version of an EccIMRpD approximant,
 *  eccentricity introduces modifications to the static functions
 *  for the inspiral phase and phase glueing:
 *  PhiInsAnsatzInt(), DPhiInsAnsatzInt(),
 *  ComputeIMRPhenDPhaseConnectionCoefficients(), IMRPhenDPhase()
 *  Functions that for the purpose of this `effective model` are assumed
 *  to be unaffected by eccentricity (i.e., higher-order terms in the
 *  inspiral phase, merger and ringdown phase functions, amplitude, ...)
 *  are for the time being copied over from the IMRPhenomD source code
 *  without changes, to enable waveform generation with EccIMRpD while
 *  avoiding modifications to the IMRPhenomD source code.
 */

/*
 * IMRPHENOMD SOURCE CODE:
 *
 * Copyright (C) 2015 Michael Puerrer, Sebastian Khan, Frank Ohme, Ofek Birnholtz, Lionel London
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
/*#include <complex.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <lal/LALStdlib.h>
#include <lal/LALSimIMR.h>
#include <lal/LALConstants.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
*/
#include <gsl/gsl_math.h>
#include "LALSimEccIMRpD_internals.c"
UsefulPowers Powers_of_Pi;	// declared in LALSimEccIMRpD_internals.c

#ifndef _OPENMP
#define omp ignore
#endif

/*
 * private function prototypes; all internal functions use solar masses.
 *
 */

static int EccIMRPhenomDGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< [out] FD waveform */
    const REAL8 phi0,                  /**< phase at fRef */
    const REAL8 fRef,                  /**< reference frequency [Hz] */
    const REAL8 deltaF,                /**< frequency resolution */
    const REAL8 m1,                    /**< mass of companion 1 [solar masses] */
    const REAL8 m2,                    /**< mass of companion 2 [solar masses] */
    const REAL8 chi1,                  /**< aligned-spin of companion 1 */
    const REAL8 chi2,                  /**< aligned-spin of companion 2 */
    const REAL8 f_min,                 /**< start frequency */
    const REAL8 f_max,                 /**< end frequency */
    const REAL8 distance,              /**< distance to source (m) */
    const REAL8 eccentricity,	       /**< eccentricity at reference frequency (here: f_min) */
    LALDict *extraParams /**< linked list containing the extra testing GR parameters */
);

/**
 * @addtogroup LALSimIMRPhenom_c
 * @{
 *
 * @name Routines for IMR Phenomenological Model "D"
 * @{
 *
 * @author Michael Puerrer, Sebastian Khan, Frank Ohme
 *
 * @brief C code for IMRPhenomD phenomenological waveform model.
 *
 * This is an aligned-spin frequency domain model.
 * See Husa et al \cite Husa:2015iqa, and Khan et al \cite Khan:2015jqa
 * for details. Any studies that use this waveform model should include
 * a reference to both of these papers.
 *
 * @note The model was calibrated to mass-ratios [1:1,1:4,1:8,1:18].
 * * Along the mass-ratio 1:1 line it was calibrated to spins  [-0.95, +0.98].
 * * Along the mass-ratio 1:4 line it was calibrated to spins  [-0.75, +0.75].
 * * Along the mass-ratio 1:8 line it was calibrated to spins  [-0.85, +0.85].
 * * Along the mass-ratio 1:18 line it was calibrated to spins [-0.8, +0.4].
 * The calibration points will be given in forthcoming papers.
 *
 * @attention The model is usable outside this parameter range,
 * and in tests to date gives sensible physical results,
 * but conclusive statements on the physical fidelity of
 * the model for these parameters await comparisons against further
 * numerical-relativity simulations. For more information, see the review wiki
 * under https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/WaveformsReview/IMRPhenomDCodeReview
 */

/**
 * Driver routine to compute the eccentric variant EccIMRpD of the spin-aligned, 
 * inspiral-merger-ringdown phenomenological waveform IMRPhenomD in the frequency domain.
 *
 * Reference:
 * - Waveform: Eq. 35 and 36 in arXiv:1508.07253
 * - Coefficients: Eq. 31 and Table V in arXiv:1508.07253
 *
 *  All input parameters should be in SI units. Angles should be in radians.
 *
 * NEW:
 * Inspiral part of the phase includes leading-order eccentricity corrections
 * up to 3PN phase order
 *
 */
int XLALSimEccIMRPhenomDGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< [out] FD waveform */
    const REAL8 phi0,                  /**< Orbital phase at fRef (rad) */
    const REAL8 fRef_in,               /**< reference frequency (Hz) */
    const REAL8 deltaF,                /**< Sampling frequency (Hz) */
    const REAL8 m1_SI,                 /**< Mass of companion 1 (kg) */
    const REAL8 m2_SI,                 /**< Mass of companion 2 (kg) */
    const REAL8 chi1,                  /**< Aligned-spin parameter of companion 1 */
    const REAL8 chi2,                  /**< Aligned-spin parameter of companion 2 */
    const REAL8 f_min,                 /**< Starting GW frequency (Hz) */
    const REAL8 f_max,                 /**< End frequency; 0 defaults to Mf = \ref f_CUT */
    const REAL8 distance,               /**< Distance of source (m) */
    const REAL8 eccentricity,		/**< Eccentricity at reference frequency (here: f_min) */
    LALDict *extraParams /**< linked list containing the extra testing GR parameters */
) {
  /* external: SI; internal: solar masses */
  const REAL8 m1 = m1_SI / LAL_MSUN_SI;
  const REAL8 m2 = m2_SI / LAL_MSUN_SI;

  /* check inputs for sanity */
  XLAL_CHECK(0 != htilde, XLAL_EFAULT, "htilde is null");
  if (*htilde) XLAL_ERROR(XLAL_EFAULT);
  if (fRef_in < 0) XLAL_ERROR(XLAL_EDOM, "fRef_in must be positive (or 0 for 'ignore')\n");
  if (deltaF <= 0) XLAL_ERROR(XLAL_EDOM, "deltaF must be positive\n");
  if (m1 <= 0) XLAL_ERROR(XLAL_EDOM, "m1 must be positive\n");
  if (m2 <= 0) XLAL_ERROR(XLAL_EDOM, "m2 must be positive\n");
  if (f_min <= 0) XLAL_ERROR(XLAL_EDOM, "f_min must be positive\n");
  if (f_max < 0) XLAL_ERROR(XLAL_EDOM, "f_max must be greater than 0\n");
  if (distance <= 0) XLAL_ERROR(XLAL_EDOM, "distance must be positive\n");
  if (eccentricity < 0) XLAL_ERROR(XLAL_EDOM, "eccentricity must be positive (or 0 for standard IMRPhenomD waveforms)\n");
  if (eccentricity >= 1.) XLAL_ERROR(XLAL_EDOM, "eccentricity must be less than 1\n");

  const REAL8 q = (m1 > m2) ? (m1 / m2) : (m2 / m1);

  if (q > MAX_ALLOWED_MASS_RATIO)
    XLAL_PRINT_WARNING("Warning: The model is not supported for high mass ratio, see MAX_ALLOWED_MASS_RATIO\n");

  if (chi1 > 1.0 || chi1 < -1.0 || chi2 > 1.0 || chi2 < -1.0)
    XLAL_ERROR(XLAL_EDOM, "Spins outside the range [-1,1] are not supported\n");

  // if no reference frequency given, set it to the starting GW frequency
  REAL8 fRef = (fRef_in == 0.0) ? f_min : fRef_in;

  const REAL8 M_sec = (m1+m2) * LAL_MTSUN_SI; // Conversion factor Hz -> dimensionless frequency
  const REAL8 fCut = f_CUT/M_sec; // convert Mf -> Hz
  // Somewhat arbitrary end point for the waveform.
  // Chosen so that the end of the waveform is well after the ringdown.
  if (fCut <= f_min)
    XLAL_ERROR(XLAL_EDOM, "(fCut = %g Hz) <= f_min = %g\n", fCut, f_min);

    /* default f_max to Cut */
  REAL8 f_max_prime = f_max;
  f_max_prime = f_max ? f_max : fCut;
  f_max_prime = (f_max_prime > fCut) ? fCut : f_max_prime;
  if (f_max_prime <= f_min)
    XLAL_ERROR(XLAL_EDOM, "f_max <= f_min\n");

  int status = EccIMRPhenomDGenerateFD(htilde, phi0, fRef, deltaF,
                                    m1, m2, chi1, chi2,
                                    f_min, f_max_prime, distance, eccentricity, extraParams);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to generate EccIMRpD waveform.");

  if (f_max_prime < f_max) {
    // The user has requested a higher f_max than Mf=fCut.
    // Resize the frequency series to fill with zeros beyond the cutoff frequency.
    size_t n = (*htilde)->data->length;
    size_t n_full = NextPow2(f_max / deltaF) + 1; // we actually want to have the length be a power of 2 + 1
    *htilde = XLALResizeCOMPLEX16FrequencySeries(*htilde, 0, n_full);
    XLAL_CHECK ( *htilde, XLAL_ENOMEM, "Failed to resize waveform COMPLEX16FrequencySeries of length %zu (for internal fCut=%f) to new length %zu (for user-requested f_max=%f).", n, fCut, n_full, f_max );
  }

  return XLAL_SUCCESS;
}

/** @} */

/** @} */

/* *********************************************************************************/
/* The following private function generates EccIMRpD frequency-domain waveforms  */
/* given coefficients */
/* *********************************************************************************/

static int EccIMRPhenomDGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< [out] FD waveform */
    const REAL8 phi0,                  /**< phase at fRef */
    const REAL8 fRef,                  /**< reference frequency [Hz] */
    const REAL8 deltaF,                /**< frequency resolution */
    const REAL8 m1_in,                 /**< mass of companion 1 [solar masses] */
    const REAL8 m2_in,                 /**< mass of companion 2 [solar masses] */
    const REAL8 chi1_in,               /**< aligned-spin of companion 1 */
    const REAL8 chi2_in,               /**< aligned-spin of companion 2 */
    const REAL8 f_min,                 /**< start frequency */
    const REAL8 f_max,                 /**< end frequency */
    const REAL8 distance,              /**< distance to source (m) */
    const REAL8 eccentricity,	       /**< eccentricity at reference frequency (here: f_min) */
    LALDict *extraParams /**< linked list containing the extra testing GR parameters */
) {
  LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0, 0}

  REAL8 chi1, chi2, m1, m2;
  if (m1_in>m2_in) {
     chi1 = chi1_in;
     chi2 = chi2_in;
     m1   = m1_in;
     m2   = m2_in;
  } else { // swap spins and masses
     chi1 = chi2_in;
     chi2 = chi1_in;
     m1   = m2_in;
     m2   = m1_in;
  }

  int status = init_useful_powers(&Powers_of_Pi, LAL_PI);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initiate useful powers of pi.");

  const REAL8 M = m1 + m2;
  REAL8 eta = m1 * m2 / (M * M);
 
  if (eta > 0.25)
      nudge(&eta, 0.25, 1e-6);
  if (eta > 0.25 || eta < 0.0)
      XLAL_ERROR(XLAL_EDOM, "Unphysical eta. Must be between 0. and 0.25\n");

  const REAL8 M_sec = M * LAL_MTSUN_SI;

  /* Compute the amplitude pre-factor */
  const REAL8 amp0 = 2. * sqrt(5. / (64.*LAL_PI)) * M * LAL_MRSUN_SI * M * LAL_MTSUN_SI / distance;

  /* Coalesce at t=0 */
  // shift by overall length in time
  XLAL_CHECK ( XLALGPSAdd(&ligotimegps_zero, -1. / deltaF), XLAL_EFUNC, "Failed to shift coalescence time to t=0, tried to apply shift of -1.0/deltaF with deltaF=%g.", deltaF);

  /* Allocate htilde */
  size_t n = NextPow2(f_max / deltaF) + 1;

  *htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &ligotimegps_zero, 0.0,
      deltaF, &lalStrainUnit, n);
  XLAL_CHECK ( *htilde, XLAL_ENOMEM, "Failed to allocated waveform COMPLEX16FrequencySeries of length %zu for f_max=%f, deltaF=%g.", n, f_max, deltaF);

  memset((*htilde)->data->data, 0, n * sizeof(COMPLEX16));
  XLALUnitMultiply(&((*htilde)->sampleUnits), &((*htilde)->sampleUnits), &lalSecondUnit);

  /* range that will have actual non-zero waveform values generated */
  size_t ind_min = (size_t) (f_min / deltaF);
  size_t ind_max = (size_t) (f_max / deltaF);
  XLAL_CHECK ( (ind_max<=n) && (ind_min<=ind_max), XLAL_EDOM, "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=htilde->data>length=%zu.", ind_min, ind_max, n);

  // Calculate phenomenological parameters
  const REAL8 finspin = FinalSpin0815(eta, chi1, chi2); //FinalSpin0815 - 0815 is like a version number

  if (finspin < MIN_FINAL_SPIN)
          XLAL_PRINT_WARNING("Final spin (Mf=%g) and ISCO frequency of this system are small, \
                          the model might misbehave here.", finspin);

  IMRPhenomDAmplitudeCoefficients *pAmp = ComputeIMRPhenomDAmplitudeCoefficients(eta, chi1, chi2, finspin);
  if (!pAmp) XLAL_ERROR(XLAL_EFUNC);
  if (extraParams==NULL)
    extraParams=XLALCreateDict();
  XLALSimInspiralWaveformParamsInsertPNSpinOrder(extraParams,LAL_SIM_INSPIRAL_SPIN_ORDER_35PN);
  IMRPhenomDPhaseCoefficients *pPhi = ComputeIMRPhenomDPhaseCoefficients(eta, chi1, chi2, finspin, extraParams);
  if (!pPhi) XLAL_ERROR(XLAL_EFUNC);
  PNPhasingSeries *pn = NULL;
  PNPhasingSeries *pfv19by3 = NULL;
  PNPhasingSeries *pfv25by3 = NULL;
  PNPhasingSeries *pfv28by3 = NULL;
  PNPhasingSeries *pfv31by3 = NULL;
  PNPhasingSeries *pfv34by3 = NULL;
  PNPhasingSeries *pfv37by3 = NULL;
  XLALSimInspiralEccTF2AlignedPhasing(&pn, &pfv19by3, &pfv25by3, &pfv28by3, &pfv31by3, &pfv34by3, &pfv37by3, m1, m2, chi1, chi2, eccentricity, f_min, extraParams);
  if (!pn) XLAL_ERROR(XLAL_EFUNC);
  if (!pfv19by3) XLAL_ERROR(XLAL_EFUNC);
  if (!pfv25by3) XLAL_ERROR(XLAL_EFUNC);
  if (!pfv28by3) XLAL_ERROR(XLAL_EFUNC);
  if (!pfv31by3) XLAL_ERROR(XLAL_EFUNC);
  if (!pfv34by3) XLAL_ERROR(XLAL_EFUNC);
  if (!pfv37by3) XLAL_ERROR(XLAL_EFUNC);

  // Subtract 3PN spin-spin term below as this is in LAL's TaylorF2 implementation
  // (LALSimInspiralPNCoefficients.c -> XLALSimInspiralPNPhasing_F2), but
  REAL8 testGRcor=1.0;
  testGRcor += XLALSimInspiralWaveformParamsLookupNonGRDChi6(extraParams);

  // was not available when PhenomD was tuned.
  pn->v[6] -= (Subtract3PNSS(m1, m2, M, chi1, chi2) * pn->v[0])* testGRcor;

  PhiInsPrefactors phi_prefactors;
  status = init_phi_ins_prefactors(&phi_prefactors, pPhi, pn);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "init_phi_ins_prefactors failed");

  // Compute coefficients to make phase C^1 continuous (phase and first derivative)
  ComputeEccIMRPhenDPhaseConnectionCoefficients(pPhi, pn, &phi_prefactors, pfv19by3, pfv25by3, pfv28by3, pfv31by3, pfv34by3, pfv37by3);

  //time shift so that peak amplitude is approximately at t=0
  //For details see https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/WaveformsReview/IMRPhenomDCodeReview/timedomain
  const REAL8 t0 = DPhiMRD(pAmp->fmaxCalc, pPhi);

  AmpInsPrefactors amp_prefactors;
  status = init_amp_ins_prefactors(&amp_prefactors, pAmp);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "init_amp_ins_prefactors failed");

  // incorporating fRef
  const REAL8 MfRef = M_sec * fRef;
  UsefulPowers powers_of_fRef;
  status = init_useful_powers(&powers_of_fRef, MfRef);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers failed for MfRef");
  const REAL8 phifRef = EccIMRPhenDPhase(MfRef, pPhi, pn, &powers_of_fRef, &phi_prefactors, pfv19by3, pfv25by3, pfv28by3, pfv31by3, pfv34by3, pfv37by3);

  // factor of 2 b/c phi0 is orbital phase
  const REAL8 phi_precalc = 2.*phi0 + phifRef;

  int status_in_for = XLAL_SUCCESS;
  /* Now generate the waveform */
  #pragma omp parallel for
  for (size_t i = ind_min; i < ind_max; i++)
  {
    REAL8 Mf = M_sec * i * deltaF; // geometric frequency

    UsefulPowers powers_of_f;
    status_in_for = init_useful_powers(&powers_of_f, Mf);
    if (XLAL_SUCCESS != status_in_for)
    {
      XLALPrintError("init_useful_powers failed for Mf, status_in_for=%d", status_in_for);
      status = status_in_for;
    }
    else
    {
      REAL8 amp = IMRPhenDAmplitude(Mf, pAmp, &powers_of_f, &amp_prefactors);
      REAL8 phi = EccIMRPhenDPhase(Mf, pPhi, pn, &powers_of_f, &phi_prefactors, pfv19by3, pfv25by3, pfv28by3, pfv31by3, pfv34by3, pfv37by3);

      phi -= t0*(Mf-MfRef) + phi_precalc;
      ((*htilde)->data->data)[i] = amp0 * amp * cexp(-I * phi);
    }
  }

  LALFree(pAmp);
  LALFree(pPhi);
  LALFree(pn);
  LALFree(pfv19by3);
  LALFree(pfv25by3);
  LALFree(pfv28by3);
  LALFree(pfv31by3);
  LALFree(pfv34by3);
  LALFree(pfv37by3);

  return status;
}

// Taken from LALSimIMRPhenomP.c
// This function determines whether x and y are approximately equal to a relative accuracy epsilon.
// Note that x and y are compared to relative accuracy, so this function is not suitable for testing whether a value is approximately zero.
static bool approximately_equal(REAL8 x, REAL8 y, REAL8 epsilon) {
  return !gsl_fcmp(x, y, epsilon);
}

// If x and X are approximately equal to relative accuracy epsilon then set x = X.
// If X = 0 then use an absolute comparison.
// Taken from LALSimIMRPhenomP.c
static void nudge(REAL8 *x, REAL8 X, REAL8 epsilon) {
  if (X != 0.0) {
    if (approximately_equal(*x, X, epsilon)) {
      XLAL_PRINT_INFO("Nudging value %.15g to %.15g\n", *x, X);
      *x = X;
    }
  }
  else {
    if (fabs(*x - X) < epsilon)
      *x = X;
  }
}
