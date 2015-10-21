/*
 * Copyright (C) 2015 Michael Puerrer, Sebastian Khan, Frank Ohme
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
#include <complex.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <lal/LALStdlib.h>
#include <lal/LALSimIMR.h>
#include <lal/LALConstants.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>

#include "LALSimIMRPhenomD_internals.c"

#ifndef _OPENMP
#define omp ignore
#endif

/*
 * private function prototypes; all internal functions use solar masses.
 *
 */

static int IMRPhenomDGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< FD waveform */
    const REAL8 phi0,                  /**< phase at fRef */
    const REAL8 fRef,                  /**< reference frequency [Hz] */
    const REAL8 deltaF,                /**< frequency resolution */
    const REAL8 m1,                    /**< mass of companion 1 [solar masses] */
    const REAL8 m2,                    /**< mass of companion 2 [solar masses] */
    const REAL8 chi1,                  /**< aligned-spin of companion 1 */
    const REAL8 chi2,                  /**< aligned-spin of companion 2 */
    const REAL8 f_min,                 /**< start frequency */
    const REAL8 f_max,                 /**< end frequency */
    const REAL8 distance               /**< distance to source (m) */
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
 * Driver routine to compute the spin-aligned, inspiral-merger-ringdown
 * phenomenological waveform IMRPhenomD in the frequency domain.
 *
 * Reference:
 * - Waveform: Eq. 35 and 36 in arXiv:1508.07253
 * - Coefficients: Eq. 31 and Table V in arXiv:1508.07253
 *
 *  All input parameters should be in SI units. Angles should be in radians.
 */
int XLALSimIMRPhenomDGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< FD waveform */
    const REAL8 phi0,                  /**< Orbital phase at fRef (rad) */
    const REAL8 fRef_in,                  /**< reference frequency (Hz) */
    const REAL8 deltaF,                /**< Sampling frequency (Hz) */
    const REAL8 m1_SI,                 /**< Mass of companion 1 (kg) */
    const REAL8 m2_SI,                 /**< Mass of companion 2 (kg) */
    const REAL8 chi1,                  /**< Aligned-spin parameter of companion 1 */
    const REAL8 chi2,                  /**< Aligned-spin parameter of companion 2 */
    const REAL8 f_min,                 /**< Starting GW frequency (Hz) */
    const REAL8 f_max,                 /**< End frequency; 0 defaults to Mf = 0.3 */
    const REAL8 distance               /**< Distance of source (m) */
) {
  /* external: SI; internal: solar masses */
  const REAL8 m1 = m1_SI / LAL_MSUN_SI;
  const REAL8 m2 = m2_SI / LAL_MSUN_SI;

  /* check inputs for sanity */
  XLAL_CHECK(0 != htilde, XLAL_EFAULT, "htilde is null");
  if (*htilde) XLAL_ERROR(XLAL_EFAULT);
  if (deltaF <= 0) XLAL_ERROR(XLAL_EDOM, "deltaF must be positive\n");
  if (m1 <= 0) XLAL_ERROR(XLAL_EDOM, "m1 must be positive\n");
  if (m2 <= 0) XLAL_ERROR(XLAL_EDOM, "m2 must be positive\n");
  if (f_min <= 0) XLAL_ERROR(XLAL_EDOM, "f_min must be positive\n");
  if (f_max < 0) XLAL_ERROR(XLAL_EDOM, "f_max must be greater than 0\n");
  if (distance <= 0) XLAL_ERROR(XLAL_EDOM, "distance must be positive\n");

  if (chi1 > 1.0 || chi1 < -1.0 || chi2 > 1.0 || chi2 < -1.0)
    XLAL_ERROR(XLAL_EDOM, "Spins outside the range [-1,1] are not supported\n");

  // if no reference frequency given, set it to the starting GW frequency
  REAL8 fRef = (fRef_in == 0.0) ? f_min : fRef_in;

  const REAL8 M_sec = (m1+m2) * LAL_MTSUN_SI; // Conversion factor Hz -> dimensionless frequency
  const REAL8 fCut = f_CUT/M_sec; // convert Mf -> Hz
  // Somewhat arbitrary end point for the waveform.
  // Chosen so that the end of the waveform the well after the ringdown.
  if (fCut <= f_min)
    XLAL_ERROR(XLAL_EDOM, "(fCut = %g Hz) <= f_min = %g\n", fCut, f_min);

    /* default f_max to Cut */
  REAL8 f_max_prime = f_max;
  f_max_prime = f_max ? f_max : fCut;
  f_max_prime = (f_max_prime > fCut) ? fCut : f_max_prime;
  if (f_max_prime <= f_min)
    XLAL_ERROR(XLAL_EDOM, "f_max <= f_min\n");

  REAL8 status = IMRPhenomDGenerateFD(htilde, phi0, fRef, deltaF,
                                       m1, m2, chi1, chi2,
                                       f_min, f_max_prime, distance);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "return status is not XLAL_SUCCESS");

  if (f_max_prime < f_max) {
    // The user has requested a higher f_max than Mf=fCut.
    // Resize the frequency series to fill with zeros beyond the cutoff frequency.
    size_t n_full = NextPow2(f_max / deltaF) + 1; // we actually want to have the length be a power of 2 + 1
    *htilde = XLALResizeCOMPLEX16FrequencySeries(*htilde, 0, n_full);
  }

  return XLAL_SUCCESS;
}

/** @} */

/** @} */

/* *********************************************************************************/
/* The following private function generates IMRPhenomD frequency-domain waveforms  */
/* given coefficients */
/* *********************************************************************************/

static int IMRPhenomDGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< FD waveform */
    const REAL8 phi0,                  /**< phase at fRef */
    const REAL8 fRef,                  /**< reference frequency [Hz] */
    const REAL8 deltaF,                /**< frequency resolution */
    const REAL8 m1_in,                 /**< mass of companion 1 [solar masses] */
    const REAL8 m2_in,                 /**< mass of companion 2 [solar masses] */
    const REAL8 chi1_in,               /**< aligned-spin of companion 1 */
    const REAL8 chi2_in,               /**< aligned-spin of companion 2 */
    const REAL8 f_min,                 /**< start frequency */
    const REAL8 f_max,                 /**< end frequency */
    const REAL8 distance               /**< distance to source (m) */
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

  const REAL8 M = m1 + m2;
  REAL8 eta = m1 * m2 / (M * M);

  if (eta > 0.25 || eta < 0.0)
    XLAL_ERROR(XLAL_EDOM, "Unphysical eta. Must be between 0. and 0.25\n");

  const REAL8 M_sec = M * LAL_MTSUN_SI;

  /* Compute the amplitude pre-factor */
  REAL8 amp0 = 2. * sqrt(5. / (64.*LAL_PI)) * M * LAL_MRSUN_SI * M * LAL_MTSUN_SI / distance;

  /* Allocate htilde */
  size_t n = NextPow2(f_max / deltaF) + 1;
  /* Coalesce at t=0 */
  XLALGPSAdd(&ligotimegps_zero, -1. / deltaF); // shift by overall length in time
  *htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &ligotimegps_zero, 0.0,
      deltaF, &lalStrainUnit, n);
  memset((*htilde)->data->data, 0, n * sizeof(COMPLEX16));
  XLALUnitMultiply(&((*htilde)->sampleUnits), &((*htilde)->sampleUnits), &lalSecondUnit);
  if (!(*htilde)) XLAL_ERROR(XLAL_EFUNC);

  size_t ind_min = (size_t) (f_min / deltaF);
  size_t ind_max = (size_t) (f_max / deltaF);

  // Calculate phenomenological parameters
  REAL8 finspin = FinalSpin0815(eta, chi1, chi2); //FinalSpin0815 - 0815 is like a version number

  if (finspin < MIN_FINAL_SPIN)
	  XLAL_PRINT_WARNING("Final spin (Mf=%g) and ISCO frequency of this system are small, \
			  the model might misbehave here.", finspin);

  IMRPhenomDAmplitudeCoefficients *pAmp = ComputeIMRPhenomDAmplitudeCoefficients(eta, chi1, chi2, finspin);
  if (!pAmp) XLAL_ERROR(XLAL_EFUNC);
  IMRPhenomDPhaseCoefficients *pPhi = ComputeIMRPhenomDPhaseCoefficients(eta, chi1, chi2, finspin);
  if (!pPhi) XLAL_ERROR(XLAL_EFUNC);
  PNPhasingSeries *pn = NULL;
  XLALSimInspiralTaylorF2AlignedPhasing(&pn, m1, m2, chi1, chi2, 1.0, 1.0, LAL_SIM_INSPIRAL_SPIN_ORDER_35PN);
  if (!pn) XLAL_ERROR(XLAL_EFUNC);

  // Subtract 3PN spin-spin term below as this is in LAL's TaylorF2 implementation
  // (LALSimInspiralPNCoefficients.c -> XLALSimInspiralPNPhasing_F2), but
  // was not available as PhenomD was tuned.
  REAL8 chi1sq = chi1 * chi1;
  REAL8 chi2sq = chi2 * chi2;
  REAL8 m1M = m1 / M;
  REAL8 m2M = m2 / M;
  REAL8 pn_ss3 =  (326.75L/1.12L + 557.5L/1.8L*eta)*eta*chi1*chi2;
  pn_ss3 += ((4703.5L/8.4L+2935.L/6.L*m1M-120.L*m1M*m1M) + (-4108.25L/6.72L-108.5L/1.2L*m1M+125.5L/3.6L*m1M*m1M)) *m1M*m1M * chi1sq;
  pn_ss3 += ((4703.5L/8.4L+2935.L/6.L*m2M-120.L*m2M*m2M) + (-4108.25L/6.72L-108.5L/1.2L*m2M+125.5L/3.6L*m2M*m2M)) *m2M*m2M * chi2sq;

  pn->v[6] -= (pn_ss3 * pn->v[0]);

  // Compute coefficients to make phase C^1 continuous (phase and first derivative)
  ComputeIMRPhenDPhaseConnectionCoefficients(pPhi, pn);

  //time shift so that peak amplitude is approximately at t=0
  //For details see https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/WaveformsReview/IMRPhenomDCodeReview/timedomain
  REAL8 t0 = DPhiMRD(0.8*(pPhi->fRD), pPhi);

  /* Now generate the waveform */
  #pragma omp parallel for
  for (size_t i = ind_min; i < ind_max; i++) {

    REAL8 Mf = M_sec * i * deltaF; // geometric frequency

    REAL8 amp = IMRPhenDAmplitude(Mf, pAmp);
    REAL8 phi = IMRPhenDPhase(Mf, pPhi, pn);

    // incorporating fRef
    REAL8 MfRef = M_sec * fRef;
    REAL8 phifRef = IMRPhenDPhase(MfRef, pPhi, pn);
    phi -= 2.*phi0 + t0*(Mf-MfRef) + phifRef; // factor of 2 b/c phi0 is orbital phase

    ((*htilde)->data->data)[i] = amp0 * amp * cexp(-I * phi);
  }

  LALFree(pAmp);
  LALFree(pPhi);
  LALFree(pn);

  return XLAL_SUCCESS;
}
