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
    const REAL8 phi0,                  /**< phase at peak */
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
 * See Husa et al, arXiv:1508.07250 and Kahn et al, arXiv:1508.07253 for details.
 *
 * This is an aligned-spin frequency domain model.
 *
 * @note The model was calibrated to mass-ratios [1:1,1:4,1:8,1:18].
 * * Along the mass-ratio 1:1 line it was calibrated to spins  [-0.95, +0.98].
 * * Along the mass-ratio 1:4 line it was calibrated to spins  [-0.75, +0.75].
 * * Along the mass-ratio 1:8 line it was calibrated to spins  [-0.85, +0.85].
 * * Along the mass-ratio 1:18 line it was calibrated to spins [-0.8, +0.4].
 * The calibration points will be given in forthcoming papers.
 *
 * @note The model is usable outside this parameter range,
 * and in tests to date gives sensible physical results,
 * but conclusive statements on the physical fidelity of
 * the model for these parameters await comparisons against further
 * numerical-relativity simulations.
 */


/**
 * Driver routine to compute the spin-aligned, inspiral-merger-ringdown
 * phenomenological waveform IMRPhenomD in the frequency domain.
 *
 * Reference:
 * - Waveform: Eq.
 * - Coefficients: Eq. and Table xyz
 *
 *  All input parameters should be in SI units. Angles should be in radians.
 */
int XLALSimIMRPhenomDGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< FD waveform */
    const REAL8 phi0,                  /**< Orbital phase at peak (rad) */
    const REAL8 deltaF,                /**< Sampling frequency (Hz) */
    const REAL8 m1_SI,                 /**< Mass of companion 1 (kg) */
    const REAL8 m2_SI,                 /**< Mass of companion 2 (kg) */
    const REAL8 chi1,                  /**< Aligned-spin of companion 1 */
    const REAL8 chi2,                  /**< Aligned-spin of companion 2 */
    const REAL8 f_min,                 /**< Starting GW frequency (Hz) */
    const REAL8 f_max,                 /**< End frequency; 0 defaults to ringdown cutoff freq */
    const REAL8 distance               /**< Distance of source (m) */
) {
  /* external: SI; internal: solar masses */
  const REAL8 m1 = m1_SI / LAL_MSUN_SI;
  const REAL8 m2 = m2_SI / LAL_MSUN_SI;
  const REAL8 q = (m1 > m2) ? (m1 / m2) : (m2 / m1);

  /* check inputs for sanity */
  if (*htilde) XLAL_ERROR(XLAL_EFAULT);
  if (deltaF <= 0) XLAL_ERROR(XLAL_EDOM);
  if (m1 <= 0) XLAL_ERROR(XLAL_EDOM);
  if (m2 <= 0) XLAL_ERROR(XLAL_EDOM);
  if (fabs(chi1) > 1 || fabs(chi2) > 1) XLAL_ERROR(XLAL_EDOM);
  if (f_min <= 0) XLAL_ERROR(XLAL_EDOM);
  if (f_max < 0) XLAL_ERROR(XLAL_EDOM);
  if (distance <= 0) XLAL_ERROR(XLAL_EDOM);

  if (chi1 > 0.99 || chi1 < -1.0 || chi2 > 0.99 || chi2 < -1.0)
    XLAL_ERROR(XLAL_EDOM, "Spins outside the range [-1,0.99] are not supported\n");

  if (q > 18.0)
    XLAL_PRINT_WARNING("Warning: The model is calibrated up to m1/m2 <= 18.\n");

  const REAL8 M_sec = (m1+m2) * LAL_MTSUN_SI;
  const REAL8 fCut = 0.3/M_sec;
  if (fCut <= f_min)
    XLAL_ERROR(XLAL_EDOM, "(fCut = %gM) <= f_min = %g\n", fCut, f_min);

  /* default f_max to params->fCut */
  REAL8 f_max_prime = f_max;
  f_max_prime = f_max ? f_max : fCut;
  f_max_prime = (f_max_prime > fCut) ? fCut : f_max_prime;
  if (f_max_prime <= f_min)
    XLAL_ERROR(XLAL_EDOM, "f_max <= f_min\n");

  REAL8 status = IMRPhenomDGenerateFD(htilde, phi0, deltaF,
                                      m1, m2, chi1, chi2,
                                      f_min, f_max_prime, distance);

  if (f_max_prime < f_max) {
    // The user has requested a higher f_max than Mf=params->fCut.
    // Resize the frequency series to fill with zeros to fill with zeros beyond the cutoff frequency.
    size_t n_full = NextPow2(f_max / deltaF) + 1; // we actually want to have the length be a power of 2 + 1
    *htilde = XLALResizeCOMPLEX16FrequencySeries(*htilde, 0, n_full);
  }

  return status;
}

/** @} */

/** @} */

/* *********************************************************************************/
/* The following private function generates IMRPhenomD frequency-domain waveforms  */
/* given coefficients */
/* *********************************************************************************/

static int IMRPhenomDGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< FD waveform */
    const REAL8 phi0,                  /**< phase at peak */
    const REAL8 deltaF,                /**< frequency resolution */
    const REAL8 m1,                    /**< mass of companion 1 [solar masses] */
    const REAL8 m2,                    /**< mass of companion 2 [solar masses] */
    const REAL8 chi1_in,               /**< aligned-spin of companion 1 */
    const REAL8 chi2_in,               /**< aligned-spin of companion 2 */
    const REAL8 f_min,                 /**< start frequency */
    const REAL8 f_max,                 /**< end frequency */
    const REAL8 distance               /**< distance to source (m) */
) {
  LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0, 0}

  const REAL8 M = m1 + m2;
  REAL8 eta = m1 * m2 / (M * M);
  const REAL8 M_sec = M * LAL_MTSUN_SI;

  REAL8 chi1, chi2;
  if (m1>m2) { // swap spins
    chi1 = chi1_in;
    chi2 = chi2_in;
  } else {
    chi1 = chi2_in;
    chi2 = chi1_in;
  }

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
  REAL8 finspin = FinalSpin0714(eta, chi1, chi2);
  IMRPhenomDAmplitudeCoefficients *pAmp = ComputeIMRPhenomDAmplitudeCoefficients(eta, chi1, chi2, finspin);
  IMRPhenomDPhaseCoefficients *pPhi = ComputeIMRPhenomDPhaseCoefficients(eta, chi1, chi2, finspin);
  if (!pAmp || !pPhi) XLAL_ERROR(XLAL_EFUNC);

  // Compute coefficients to make phase C^1
  ComputeIMRPhenDPhaseConnectionCoefficients(pPhi);

  /* Now generate the waveform */
  #pragma omp parallel for
  for (size_t i = ind_min; i < ind_max; i++) {

    REAL8 Mf = M_sec * i * deltaF; // geometric frequency

    REAL8 amp = IMRPhenDAmplitude(Mf, pAmp);
    REAL8 phi = IMRPhenDPhase(Mf, pPhi);

    phi -= 2.*phi0; // factor of 2 b/c phi0 is orbital phase

    ((*htilde)->data->data)[i] = amp0 * amp * cexp(-I * phi);
  }

  LALFree(pAmp);
  LALFree(pPhi);

  return XLAL_SUCCESS;
}
