/*
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */


#include <math.h>
#include <gsl/gsl_math.h>
#include "LALSimIMRPhenomD_internals.c"
#include <lal/Sequence.h>

#include "LALSimIMRPhenomInternalUtils.h"
#include "LALSimIMRPhenomUtils.h"

UsefulPowers powers_of_pi;	// declared in LALSimIMRPhenomD_internals.c

#ifndef _OPENMP
#define omp ignore
#endif

/*
 * private function prototypes; all internal functions use solar masses.
 *
 */

static int IMRPhenomDGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< [out] FD waveform */
    const REAL8Sequence *freqs_in,     /**< Frequency points at which to evaluate the waveform (Hz) */
    double deltaF,                     /**< If deltaF > 0, the frequency points given in freqs are uniformly spaced with
                                        * spacing deltaF. Otherwise, the frequency points are spaced non-uniformly.
                                        * Then we will use deltaF = 0 to create the frequency series we return. */
    const REAL8 phi0,                  /**< phase at fRef */
    const REAL8 fRef,                  /**< reference frequency [Hz] */
    const REAL8 m1_in,                 /**< mass of companion 1 [solar masses] */
    const REAL8 m2_in,                 /**< mass of companion 2 [solar masses] */
    const REAL8 chi1_in,               /**< aligned-spin of companion 1 */
    const REAL8 chi2_in,               /**< aligned-spin of companion 2 */
    const REAL8 distance,              /**< distance to source (m) */
    LALDict *extraParams, /**< linked list containing the extra testing GR parameters */
    NRTidal_version_type NRTidal_version /**< NRTidal version; either NRTidal_V or NRTidalv2_V or NoNRT_V in case of BBH baseline */
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
 *
 * Compute waveform in LAL format for the IMRPhenomD model.
 *
 * Returns the plus and cross polarizations as a complex frequency series with
 * equal spacing deltaF and contains zeros from zero frequency to the starting
 * frequency fLow and zeros beyond the cutoff frequency in the ringdown.
 */
int XLALSimIMRPhenomDGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< [out] FD waveform */
    const REAL8 phi0,                  /**< Orbital phase at fRef (rad) */
    const REAL8 fRef_in,               /**< reference frequency (Hz) */
    const REAL8 deltaF,                /**< Sampling frequency (Hz) */
    const REAL8 m1_SI,                 /**< Mass of companion 1 (kg) */
    const REAL8 m2_SI,                 /**< Mass of companion 2 (kg) */
    const REAL8 chi1,                  /**< Aligned-spin parameter of companion 1 */
    const REAL8 chi2,                  /**< Aligned-spin parameter of companion 2 */
    const REAL8 f_min,                 /**< Starting GW frequency (Hz) */
    const REAL8 f_max,                 /**< End frequency; 0 defaults to Mf = f_CUT */
    const REAL8 distance,               /**< Distance of source (m) */
    LALDict *extraParams, /**< linked list containing the extra testing GR parameters */
    NRTidal_version_type NRTidal_version /**< Version of NRTides; can be one of NRTidal versions or NoNRT_V for the BBH baseline */
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

  // Use fLow, fHigh, deltaF to compute freqs sequence
  // Instead of building a full sequency we only transfer the boundaries and let
  // the internal core function do the rest (and properly take care of corner cases).
  REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
  freqs->data[0] = f_min;
  freqs->data[1] = f_max_prime;
  int status = IMRPhenomDGenerateFD(htilde, freqs, deltaF, phi0, fRef,
                                    m1, m2, chi1, chi2,
                                    distance, extraParams, NRTidal_version);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to generate IMRPhenomD waveform.");
  XLALDestroyREAL8Sequence(freqs);

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

/**
 * Compute waveform in LAL format at specified frequencies for the IMRPhenomD model.
 *
 * XLALSimIMRPhenomDGenerateFD() returns the plus and cross polarizations as a complex
 * frequency series with equal spacing deltaF and contains zeros from zero frequency
 * to the starting frequency and zeros beyond the cutoff frequency in the ringdown.
 *
 * In contrast, XLALSimIMRPhenomDFrequencySequence() returns a
 * complex frequency series with entries exactly at the frequencies specified in
 * the sequence freqs (which can be unequally spaced). No zeros are added.
 *
 * If XLALSimIMRPhenomDFrequencySequence() is called with frequencies that
 * are beyond the maxium allowed geometric frequency for the ROM, zero strain is returned.
 *
 * This function is designed as an entry point for reduced order quadratures.
 */
int XLALSimIMRPhenomDFrequencySequence(
    COMPLEX16FrequencySeries **htilde,           /**< [out] FD waveform */
    const REAL8Sequence *freqs,                  /**< Frequency points at which to evaluate the waveform (Hz) */
    const REAL8 phi0,                            /**< Orbital phase at fRef (rad) */
    const REAL8 fRef_in,                         /**< reference frequency (Hz) */
    const REAL8 m1_SI,                           /**< Mass of companion 1 (kg) */
    const REAL8 m2_SI,                           /**< Mass of companion 2 (kg) */
    const REAL8 chi1,                            /**< Aligned-spin parameter of companion 1 */
    const REAL8 chi2,                            /**< Aligned-spin parameter of companion 2 */
    const REAL8 distance,                        /**< Distance of source (m) */
    LALDict *extraParams, /**< linked list containing the extra testing GR parameters */
    NRTidal_version_type NRTidal_version /**< NRTidal version; either NRTidal_V or NRTidalv2_V or NoNRT_V in case of BBH baseline */
) {
  /* external: SI; internal: solar masses */
  const REAL8 m1 = m1_SI / LAL_MSUN_SI;
  const REAL8 m2 = m2_SI / LAL_MSUN_SI;

  /* check inputs for sanity */
  XLAL_CHECK(0 != htilde, XLAL_EFAULT, "htilde is null");
  if (*htilde) XLAL_ERROR(XLAL_EFAULT);
  if (!freqs) XLAL_ERROR(XLAL_EFAULT);
  if (fRef_in < 0) XLAL_ERROR(XLAL_EDOM, "fRef_in must be positive (or 0 for 'ignore')\n");
  if (m1 <= 0) XLAL_ERROR(XLAL_EDOM, "m1 must be positive\n");
  if (m2 <= 0) XLAL_ERROR(XLAL_EDOM, "m2 must be positive\n");
  if (distance <= 0) XLAL_ERROR(XLAL_EDOM, "distance must be positive\n");

  const REAL8 q = (m1 > m2) ? (m1 / m2) : (m2 / m1);

  if (q > MAX_ALLOWED_MASS_RATIO)
    XLAL_PRINT_WARNING("Warning: The model is not supported for high mass ratio, see MAX_ALLOWED_MASS_RATIO\n");

  if (chi1 > 1.0 || chi1 < -1.0 || chi2 > 1.0 || chi2 < -1.0)
    XLAL_ERROR(XLAL_EDOM, "Spins outside the range [-1,1] are not supported\n");

  // if no reference frequency given, set it to the starting GW frequency
  REAL8 fRef = (fRef_in == 0.0) ? freqs->data[0] : fRef_in;

  int status = IMRPhenomDGenerateFD(htilde, freqs, 0, phi0, fRef,
                                    m1, m2, chi1, chi2,
                                    distance, extraParams, NRTidal_version);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to generate IMRPhenomD waveform.");

  return XLAL_SUCCESS;
}


/** @} */

/** @} */

/* *********************************************************************************/
/* The following private function generates IMRPhenomD frequency-domain waveforms  */
/* given coefficients */
/* *********************************************************************************/

static int IMRPhenomDGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< [out] FD waveform */
    const REAL8Sequence *freqs_in,     /**< Frequency points at which to evaluate the waveform (Hz) */
    double deltaF,                     /* If deltaF > 0, the frequency points given in freqs are uniformly spaced with
                                        * spacing deltaF. Otherwise, the frequency points are spaced non-uniformly.
                                        * Then we will use deltaF = 0 to create the frequency series we return. */
    const REAL8 phi0,                  /**< phase at fRef */
    const REAL8 fRef,                  /**< reference frequency [Hz] */
    const REAL8 m1_in,                 /**< mass of companion 1 [solar masses] */
    const REAL8 m2_in,                 /**< mass of companion 2 [solar masses] */
    const REAL8 chi1_in,               /**< aligned-spin of companion 1 */
    const REAL8 chi2_in,               /**< aligned-spin of companion 2 */
    const REAL8 distance,              /**< distance to source (m) */
    LALDict *extraParams, /**< linked list containing the extra testing GR parameters */
    NRTidal_version_type NRTidal_version /**< NRTidal version; either NRTidal_V or NRTidalv2_V or NoNRT_V in case of BBH baseline */
) {
  LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0, 0}

  // Make a pointer to LALDict to circumvent a memory leak
  // At the end we will check if we created a LALDict in extraParams
  // and destroy it if we did.
  LALDict *extraParams_in = extraParams;
  REAL8Sequence *amp_tidal = NULL; /* Tidal amplitude series; required only for IMRPhenomD_NRTidalv2 */
  REAL8 dquadmon1_in = 0., dquadmon2_in = 0., lambda1_in = 0, lambda2_in = 0.;
  if (NRTidal_version == NRTidalv2_V) {
    dquadmon1_in = XLALSimInspiralWaveformParamsLookupdQuadMon1(extraParams);
    dquadmon2_in = XLALSimInspiralWaveformParamsLookupdQuadMon2(extraParams);
    lambda1_in = XLALSimInspiralWaveformParamsLookupTidalLambda1(extraParams);
    lambda2_in = XLALSimInspiralWaveformParamsLookupTidalLambda2(extraParams);
  }

  REAL8 chi1, chi2, m1, m2, dquadmon1, dquadmon2, lambda1, lambda2;
  if (m1_in>=m2_in) {
     chi1 = chi1_in;
     chi2 = chi2_in;
     m1   = m1_in;
     m2   = m2_in;
     dquadmon1 = dquadmon1_in;
     dquadmon2 = dquadmon2_in;
     lambda1 = lambda1_in;
     lambda2 = lambda2_in;
  } else { // swap spins and masses
     chi1 = chi2_in;
     chi2 = chi1_in;
     m1   = m2_in;
     m2   = m1_in;
     dquadmon1 = dquadmon2_in;
     dquadmon2 = dquadmon1_in;
     lambda1 = lambda2_in;
     lambda2 = lambda1_in;
     if (NRTidal_version == NRTidalv2_V) {
       XLALSimInspiralWaveformParamsInsertdQuadMon1(extraParams, dquadmon1);
       XLALSimInspiralWaveformParamsInsertdQuadMon2(extraParams, dquadmon2);
     }
  }

  int status = init_useful_powers(&powers_of_pi, LAL_PI);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initiate useful powers of pi.");

  /* Find frequency bounds */
  if (!freqs_in || !freqs_in->data) XLAL_ERROR(XLAL_EFAULT);
  double f_min = freqs_in->data[0];
  double f_max = freqs_in->data[freqs_in->length - 1];
  XLAL_CHECK(f_min > 0, XLAL_EDOM, "Minimum frequency must be positive.\n");
  XLAL_CHECK(f_max >= 0, XLAL_EDOM, "Maximum frequency must be non-negative.\n");

  const REAL8 M = m1 + m2;
  REAL8 eta = m1 * m2 / (M * M);

  if (eta > 0.25)
      PhenomInternal_nudge(&eta, 0.25, 1e-6);
  if (eta > 0.25 || eta < 0.0)
      XLAL_ERROR(XLAL_EDOM, "Unphysical eta. Must be between 0. and 0.25\n");

  const REAL8 M_sec = M * LAL_MTSUN_SI;

  /* Compute the amplitude pre-factor */
  const REAL8 amp0 = 2. * sqrt(5. / (64.*LAL_PI)) * M * LAL_MRSUN_SI * M * LAL_MTSUN_SI / distance;

  size_t npts = 0;
  UINT4 offset = 0; // Index shift between freqs and the frequency series
  REAL8Sequence *freqs = NULL;
  if (deltaF > 0)  { // freqs contains uniform frequency grid with spacing deltaF; we start at frequency 0
    /* Set up output array with size closest power of 2 */
    npts = NextPow2(f_max / deltaF) + 1;
    /* Coalesce at t=0 */
    // shift by overall length in time
    XLAL_CHECK ( XLALGPSAdd(&ligotimegps_zero, -1. / deltaF), XLAL_EFUNC, "Failed to shift coalescence time to t=0, tried to apply shift of -1.0/deltaF with deltaF=%g.", deltaF);
    *htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &ligotimegps_zero, 0.0, deltaF, &lalStrainUnit, npts);
    XLAL_CHECK ( *htilde, XLAL_ENOMEM, "Failed to allocated waveform COMPLEX16FrequencySeries of length %zu for f_max=%f, deltaF=%g.", npts, f_max, deltaF);
    // Recreate freqs using only the lower and upper bounds
    size_t iStart = (size_t) (f_min / deltaF);
    size_t iStop = (size_t) (f_max / deltaF);
    XLAL_CHECK ( (iStop<=npts) && (iStart<=iStop), XLAL_EDOM, "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=htilde->data>length=%zu.", iStart, iStop, npts);
    freqs = XLALCreateREAL8Sequence(iStop - iStart);
    if (!freqs)
      XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
    for (UINT4 i=iStart; i<iStop; i++)
      freqs->data[i-iStart] = i*deltaF;
    offset = iStart;
  } else { // freqs contains frequencies with non-uniform spacing; we start at lowest given frequency
    npts = freqs_in->length;
    *htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &ligotimegps_zero, f_min, deltaF, &lalStrainUnit, npts);
    XLAL_CHECK ( *htilde, XLAL_ENOMEM, "Failed to allocated waveform COMPLEX16FrequencySeries of length %zu from sequence.", npts);
    offset = 0;
    freqs = XLALCreateREAL8Sequence(freqs_in->length);
    if (!freqs)
      XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
    for (UINT4 i=0; i<freqs_in->length; i++)
      freqs->data[i] = freqs_in->data[i];
  }

  memset((*htilde)->data->data, 0, npts * sizeof(COMPLEX16));
  XLALUnitMultiply(&((*htilde)->sampleUnits), &((*htilde)->sampleUnits), &lalSecondUnit);

  // Calculate phenomenological parameters
  const REAL8 finspin = FinalSpin0815(eta, chi1, chi2); //FinalSpin0815 - 0815 is like a version number

  if (finspin < MIN_FINAL_SPIN)
          XLAL_PRINT_WARNING("Final spin (Mf=%g) and ISCO frequency of this system are small, \
                          the model might misbehave here.", finspin);

  IMRPhenomDAmplitudeCoefficients *pAmp;
  pAmp = XLALMalloc(sizeof(IMRPhenomDAmplitudeCoefficients));
  ComputeIMRPhenomDAmplitudeCoefficients(pAmp, eta, chi1, chi2, finspin);
  if (!pAmp) XLAL_ERROR(XLAL_EFUNC);
  if (extraParams==NULL)
    extraParams=XLALCreateDict();
  XLALSimInspiralWaveformParamsInsertPNSpinOrder(extraParams,LAL_SIM_INSPIRAL_SPIN_ORDER_35PN);
  IMRPhenomDPhaseCoefficients *pPhi;
  pPhi = XLALMalloc(sizeof(IMRPhenomDPhaseCoefficients));
  ComputeIMRPhenomDPhaseCoefficients(pPhi, eta, chi1, chi2, finspin, extraParams);
  if (!pPhi) XLAL_ERROR(XLAL_EFUNC);
  PNPhasingSeries *pn = NULL;
  XLALSimInspiralTaylorF2AlignedPhasing(&pn, m1, m2, chi1, chi2, extraParams);
  if (!pn) XLAL_ERROR(XLAL_EFUNC);

  // Subtract 3PN spin-spin term below as this is in LAL's TaylorF2 implementation
  // (LALSimInspiralPNCoefficients.c -> XLALSimInspiralPNPhasing_F2), but
  REAL8 testGRcor=1.0;
  testGRcor += XLALSimInspiralWaveformParamsLookupNonGRDChi6(extraParams);

  // was not available when PhenomD was tuned.
  pn->v[6] -= (Subtract3PNSS(m1, m2, M, eta, chi1, chi2) * pn->v[0]) * testGRcor;

  PhiInsPrefactors phi_prefactors;
  status = init_phi_ins_prefactors(&phi_prefactors, pPhi, pn);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "init_phi_ins_prefactors failed");

  // Compute coefficients to make phase C^1 continuous (phase and first derivative)
  ComputeIMRPhenDPhaseConnectionCoefficients(pPhi, pn, &phi_prefactors, 1.0, 1.0);

  //time shift so that peak amplitude is approximately at t=0
  //For details see https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/WaveformsReview/IMRPhenomDCodeReview/timedomain
  const REAL8 t0 = DPhiMRD(pAmp->fmaxCalc, pPhi, 1.0, 1.0);

  AmpInsPrefactors amp_prefactors;
  status = init_amp_ins_prefactors(&amp_prefactors, pAmp);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "init_amp_ins_prefactors failed");

  // incorporating fRef
  const REAL8 MfRef = M_sec * fRef;
  UsefulPowers powers_of_fRef;
  status = init_useful_powers(&powers_of_fRef, MfRef);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers failed for MfRef");
  const REAL8 phifRef = IMRPhenDPhase(MfRef, pPhi, pn, &powers_of_fRef, &phi_prefactors, 1.0, 1.0);

  // factor of 2 b/c phi0 is orbital phase
  const REAL8 phi_precalc = 2.*phi0 + phifRef;

  int status_in_for = XLAL_SUCCESS;
  int ret = XLAL_SUCCESS;
  /* Now generate the waveform */
  if (NRTidal_version == NRTidalv2_V) {
    /* Generate the tidal amplitude (Eq. 24 of arxiv: 1905.06011) to add to BBH baseline; only for IMRPhenomD_NRTidalv2*/
    amp_tidal = XLALCreateREAL8Sequence(freqs->length);
    ret = XLALSimNRTunedTidesFDTidalAmplitudeFrequencySeries(amp_tidal, freqs, m1, m2, lambda1, lambda2);
    XLAL_CHECK(XLAL_SUCCESS == ret, ret, "Failed to generate tidal amplitude series to construct IMRPhenomD_NRTidalv2 waveform.");
    /* Generated tidal amplitude corrections */
    #pragma omp parallel for
    for (UINT4 i=0; i<freqs->length; i++) { // loop over frequency points in sequence
      double Mf = M_sec * freqs->data[i];
      double ampT = amp_tidal->data[i];
      int j = i + offset; // shift index for frequency series if needed

      UsefulPowers powers_of_f;
      status_in_for = init_useful_powers(&powers_of_f, Mf);
      if (XLAL_SUCCESS != status_in_for)
      {
        XLALPrintError("init_useful_powers failed for Mf, status_in_for=%d", status_in_for);
        status = status_in_for;
      }
      else {
        REAL8 amp = IMRPhenDAmplitude(Mf, pAmp, &powers_of_f, &amp_prefactors);
        REAL8 phi = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_f, &phi_prefactors, 1.0, 1.0);

        phi -= t0*(Mf-MfRef) + phi_precalc;
        ((*htilde)->data->data)[j] = amp0 * (amp+2*sqrt(LAL_PI/5.)*ampT) * cexp(-I * phi);
      }
    }
  } else {
      #pragma omp parallel for
      for (UINT4 i=0; i<freqs->length; i++) { // loop over frequency points in sequence
      double Mf = M_sec * freqs->data[i];
      int j = i + offset; // shift index for frequency series if needed

      UsefulPowers powers_of_f;
      status_in_for = init_useful_powers(&powers_of_f, Mf);
      if (XLAL_SUCCESS != status_in_for)
      {
        XLALPrintError("init_useful_powers failed for Mf, status_in_for=%d", status_in_for);
        status = status_in_for;
      }
      else {
        REAL8 amp = IMRPhenDAmplitude(Mf, pAmp, &powers_of_f, &amp_prefactors);
        REAL8 phi = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_f, &phi_prefactors, 1.0, 1.0);

        phi -= t0*(Mf-MfRef) + phi_precalc;
        ((*htilde)->data->data)[j] = amp0 * amp * cexp(-I * phi);
      }
    }
  }

  LALFree(pAmp);
  LALFree(pPhi);
  LALFree(pn);
  XLALDestroyREAL8Sequence(freqs);
  XLALDestroyREAL8Sequence(amp_tidal);


  /* If extraParams was allocated in this function and not passed in
   * we need to free it to prevent a leak */
  if (extraParams && !extraParams_in) {
    XLALDestroyDict(extraParams);
  } else {
    XLALSimInspiralWaveformParamsInsertPNSpinOrder(extraParams,LAL_SIM_INSPIRAL_SPIN_ORDER_ALL);
  }

  return status;
}

////////////////////////////////////////////////
// END OF REVIEWED CODE ////////////////////////
////////////////////////////////////////////////

/**
 * Function to return the frequency (in Hz) of the peak of the frequency
 * domain amplitude for the IMRPhenomD model.
 *
 * The peak is a parameter in the PhenomD model given by Eq. 20 in 1508.07253
 * where it is called f_peak in the paper.
 */
double XLALIMRPhenomDGetPeakFreq(
    const REAL8 m1_in,                 /**< mass of companion 1 [Msun] */
    const REAL8 m2_in,                 /**< mass of companion 2 [Msun] */
    const REAL8 chi1_in,               /**< aligned-spin of companion 1 */
    const REAL8 chi2_in                /**< aligned-spin of companion 2 */
) {
    // Ensure that m1 > m2 and that chi1 is the spin on m1
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
    const REAL8 M_sec = M * LAL_MTSUN_SI; // Conversion factor Hz -> dimensionless frequency

    REAL8 eta = m1 * m2 / (M * M);

    if (eta > 0.25)
        PhenomInternal_nudge(&eta, 0.25, 1e-6);
    if (eta > 0.25 || eta < 0.0)
        XLAL_ERROR(XLAL_EDOM, "Unphysical eta. Must be between 0. and 0.25\n");

    // Calculate phenomenological parameters
    REAL8 finspin = FinalSpin0815(eta, chi1, chi2);

    if (finspin < MIN_FINAL_SPIN)
          XLAL_PRINT_WARNING("Final spin (Mf=%g) and ISCO frequency of this system are small, \
                          the model might misbehave here.", finspin);
    IMRPhenomDAmplitudeCoefficients *pAmp;
    pAmp = XLALMalloc(sizeof(IMRPhenomDAmplitudeCoefficients));
    ComputeIMRPhenomDAmplitudeCoefficients(pAmp, eta, chi1, chi2, finspin);
    if (!pAmp) XLAL_ERROR(XLAL_EFUNC);

    // PeakFreq, converted to Hz
    REAL8 PeakFreq = ( pAmp->fmaxCalc ) / M_sec;

    LALFree(pAmp);

    return PeakFreq;
}


// protoype
static double PhenDPhaseDerivFrequencyPoint(double Mf, IMRPhenomDPhaseCoefficients *p, PNPhasingSeries *pn);

/**
 * Helper function to return the value of the frequency derivative of the
 * Fourier domain phase.
 * This is function is wrapped by IMRPhenomDPhaseDerivative and used
 * when estimating the length of the time domain version of the waveform.
 * unreviewed
 */
static double PhenDPhaseDerivFrequencyPoint(double Mf, IMRPhenomDPhaseCoefficients *p, PNPhasingSeries *pn)
{

  // split the calculation to just 1 of 3 possible mutually exclusive ranges

  if (!StepFunc_boolean(Mf, p->fInsJoin))	// Inspiral range
  {
      double DPhiIns = DPhiInsAnsatzInt(Mf, p, pn);
	  return DPhiIns;
  }

  if (StepFunc_boolean(Mf, p->fMRDJoin))	// MRD range
  {
      double DPhiMRDval = DPhiMRD(Mf, p, 1.0, 1.0) + p->C2MRD;
	  return DPhiMRDval;
  }

  //	Intermediate range
  double DPhiInt = DPhiIntAnsatz(Mf, p) + p->C2Int;
  return DPhiInt;
}

/**
* Estimates the length of the time domain IMRPhenomD signal
* This does NOT taking into account any tapering that is used to condition the
* Fourier domain waveform to compute the inverse Fourer transform.
* To estimate the length we assume that the waveform only reaches the
* the highest physics frequency i.e. the ringdown frequency.
* unreviewed
*/
double XLALSimIMRPhenomDChirpTime(
    const REAL8 m1_SI,                 /**< Mass of companion 1 (kg) */
    const REAL8 m2_SI,                 /**< Mass of companion 2 (kg) */
    const REAL8 chi1_in,               /**< aligned-spin of companion 1 */
    const REAL8 chi2_in,               /**< aligned-spin of companion 2 */
    const REAL8 fHzSt                  /**< arbitrary starting frequency in Hz */
) {

    if (fHzSt <= 0) XLAL_ERROR(XLAL_EDOM, "fHzSt must be positive\n");

    if (chi1_in > 1.0 || chi1_in < -1.0 || chi2_in > 1.0 || chi2_in < -1.0)
      XLAL_ERROR(XLAL_EDOM, "Spins outside the range [-1,1] are not supported\n");

    /* external: SI; internal: solar masses */
    const REAL8 m1_in = m1_SI / LAL_MSUN_SI;
    const REAL8 m2_in = m2_SI / LAL_MSUN_SI;

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

    // check that starting frequency is not higher than the peak frequency
    const REAL8 fHzPeak = XLALIMRPhenomDGetPeakFreq(m1, m2, chi1, chi2);
    if (fHzSt > fHzPeak){
        XLAL_PRINT_WARNING("Starting frequency = %f Hz is higher IMRPhenomD peak frequency %f Hz. Results may be unreliable.", fHzSt, fHzPeak);
    }

    int status = init_useful_powers(&powers_of_pi, LAL_PI);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initiate useful powers of pi.");

    const REAL8 M = m1 + m2;
    REAL8 eta = m1 * m2 / (M * M);

    if (eta > 0.25)
        PhenomInternal_nudge(&eta, 0.25, 1e-6);
    if (eta > 0.25 || eta < 0.0)
      XLAL_ERROR(XLAL_EDOM, "Unphysical eta. Must be between 0. and 0.25\n");

    // compute geometric frequency
    const REAL8 M_sec = M * LAL_MTSUN_SI;
    const REAL8 MfSt = M_sec * fHzSt;

    // Calculate phenomenological parameters
    const REAL8 finspin = FinalSpin0815(eta, chi1, chi2); //FinalSpin0815 - 0815 is like a version number

    if (finspin < MIN_FINAL_SPIN)
            XLAL_PRINT_WARNING("Final spin (Mf=%g) and ISCO frequency of this system are small, \
                            the model might misbehave here.", finspin);
    LALDict *extraParams = NULL;
    if (extraParams == NULL)
      extraParams = XLALCreateDict();
    XLALSimInspiralWaveformParamsInsertPNSpinOrder(extraParams, LAL_SIM_INSPIRAL_SPIN_ORDER_35PN);
    IMRPhenomDPhaseCoefficients *pPhi;
    pPhi = XLALMalloc(sizeof(IMRPhenomDPhaseCoefficients));
    ComputeIMRPhenomDPhaseCoefficients(pPhi, eta, chi1, chi2, finspin, extraParams);
    if (!pPhi) XLAL_ERROR(XLAL_EFUNC);
    PNPhasingSeries *pn = NULL;
    XLALSimInspiralTaylorF2AlignedPhasing(&pn, m1, m2, chi1, chi2, extraParams);
    if (!pn) XLAL_ERROR(XLAL_EFUNC);

    // Subtract 3PN spin-spin term below as this is in LAL's TaylorF2 implementation
    // (LALSimInspiralPNCoefficients.c -> XLALSimInspiralPNPhasing_F2), but
    // was not available when PhenomD was tuned.
    pn->v[6] -= (Subtract3PNSS(m1, m2, M, eta, chi1, chi2) * pn->v[0]);


    PhiInsPrefactors phi_prefactors;
    status = init_phi_ins_prefactors(&phi_prefactors, pPhi, pn);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "init_phi_ins_prefactors failed");

    // Compute coefficients to make phase C^1 continuous (phase and first derivative)
    ComputeIMRPhenDPhaseConnectionCoefficients(pPhi, pn, &phi_prefactors, 1.0, 1.0);

    // We estimate the length of the time domain signal (i.e., the chirp time)
    // By computing the difference between the values of the Fourier domain
    // phase derivative at two frequencies.
    // Here the starting frequency is an input i.e., fHzSt, converted to Geometric units MfSt
    // and the ending frequency is fixed to be the frequency of the amplitude peak in Geometric units MfPeak
    // XLALIMRPhenomDGetPeakFreq output is in Hz, covert to Mf via / M_sec
    const REAL8 MfPeak = XLALIMRPhenomDGetPeakFreq(m1, m2, chi1, chi2) / M_sec;

    // Compute phase derivative at starting frequency
    const REAL8 dphifSt = PhenDPhaseDerivFrequencyPoint(MfSt, pPhi, pn);
    // Compute phase derivative at ending (ringdown) frequency
    const REAL8 dphifRD = PhenDPhaseDerivFrequencyPoint(MfPeak, pPhi, pn);
    const REAL8 dphidiff = dphifRD - dphifSt;

    // The length of time is estimated as dphidiff / 2 / pi * M (In units of seconds)
    const REAL8 ChirpTimeSec = dphidiff / 2. / LAL_PI * M_sec;

    LALFree(pPhi);
    LALFree(pn);

    return ChirpTimeSec;

}

/**
* Function to return the final spin (spin of the remnant black hole)
* as predicted by the IMRPhenomD model. The final spin is calculated using
* the phenomenological fit described in PhysRevD.93.044006 Eq. 3.6.
* unreviewed
*/
double XLALSimIMRPhenomDFinalSpin(
    const REAL8 m1_in,                 /**< mass of companion 1 [Msun] */
    const REAL8 m2_in,                 /**< mass of companion 2 [Msun] */
    const REAL8 chi1_in,               /**< aligned-spin of companion 1 */
    const REAL8 chi2_in               /**< aligned-spin of companion 2 */
) {
    // Ensure that m1 > m2 and that chi1 is the spin on m1
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

    if (eta > 0.25)
        PhenomInternal_nudge(&eta, 0.25, 1e-6);
    if (eta > 0.25 || eta < 0.0)
        XLAL_ERROR(XLAL_EDOM, "Unphysical eta. Must be between 0. and 0.25\n");

    REAL8 finspin = FinalSpin0815(eta, chi1, chi2);

    if (finspin < MIN_FINAL_SPIN)
          XLAL_PRINT_WARNING("Final spin (Mf=%g) and ISCO frequency of this system are small, \
                          the model might misbehave here.", finspin);

    return finspin;
}


/* IMRPhenomDSetupPhase */
/* IMRPhenomDEvaluatePhaseFrequencySequence */





/**
 * Helper function used in PhenomHM and PhenomPv3HM
 * Returns the phenomD phase, with modified QNM
 */
int IMRPhenomDPhaseFrequencySequence(
    REAL8Sequence *phases, /**< [out] phase evaluated at input freqs */
    REAL8Sequence *freqs,  /**< Sequency of Geometric frequencies */
    size_t ind_min,        /**< start index for frequency loop */
    size_t ind_max,        /**< end index for frequency loop */
    REAL8 m1,              /**< mass of primary in solar masses */
    REAL8 m2,              /**< mass of secondary in solar masses */
    REAL8 chi1x,           /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 chi1y,           /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 chi1z,           /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 chi2x,           /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8 chi2y,           /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8 chi2z,           /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8 Rholm,         /**< ratio of ringdown frequencies f_RD_22/f_RD_lm */
    REAL8 Taulm,         /**< ratio of ringdown damping times f_RM_22/f_RM_lm */
    LALDict *extraParams /**< linked list containing the extra testing GR parameters */
)
{
  int retcode = 0;
  PhenDAmpAndPhasePreComp pD;
  retcode = IMRPhenomDSetupAmpAndPhaseCoefficients(
    &pD, m1, m2, chi1x, chi1y, chi1z,
    chi2x, chi2y, chi2z, Rholm, Taulm, extraParams);
  if (retcode != XLAL_SUCCESS)
  {
    XLALPrintError("XLAL Error - IMRPhenomDSetupAmpAndPhaseCoefficients failed\n");
    XLAL_ERROR(XLAL_EDOM);
  }

  int status_in_for = XLAL_SUCCESS;
  /* Now generate the waveform */
  #pragma omp parallel for
  for (size_t i = ind_min; i < ind_max; i++)
  {
    REAL8 Mf = freqs->data[i]; // geometric frequency

    UsefulPowers powers_of_f;
    status_in_for = init_useful_powers(&powers_of_f, Mf);
    if (XLAL_SUCCESS != status_in_for)
    {
      XLALPrintError("init_useful_powers failed for Mf, status_in_for=%d\n", status_in_for);
      retcode = status_in_for;
    }
    else
    {
      phases->data[i] = IMRPhenDPhase(Mf, &(pD.pPhi), &(pD.pn), &powers_of_f,
                                      &(pD.phi_prefactors), Rholm, Taulm);
    }
  }

  // LALFree(pPhi);
  // LALFree(pn);

  return XLAL_SUCCESS;
}

/* IMRPhenomDSetupAmplitude */
/* IMRPhenomDEvaluateAmplitude */

/**
 * Helper function used in PhenomHM and PhenomPv3HM
 * Returns the phenomD amplitude
 */
int IMRPhenomDAmpFrequencySequence(
    REAL8Sequence *amps,  /**< [out] phase evaluated at input freqs */
    REAL8Sequence *freqs, /**< Sequency of Geometric frequencies */
    size_t ind_min,       /**< start index for frequency loop */
    size_t ind_max,       /**< end index for frequency loop */
    REAL8 m1,             /**< mass of primary in solar masses */
    REAL8 m2,             /**< mass of secondary in solar masses */
    REAL8 chi1x,          /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 chi1y,          /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 chi1z,          /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 chi2x,          /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8 chi2y,          /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8 chi2z           /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
)
{
  int retcode;

  /* It's difficult to see in the code but you need to setup the
     * powers_of_pi.
     */
  retcode = 0;
  retcode = init_useful_powers(&powers_of_pi, LAL_PI);
  XLAL_CHECK(XLAL_SUCCESS == retcode, retcode, "Failed to initiate useful powers of pi.");

  PhenomInternal_PrecessingSpinEnforcePrimaryIsm1(&m1, &m2, &chi1x, &chi1y, &chi1z, &chi2x, &chi2y, &chi2z);
  const REAL8 Mtot = m1 + m2;
  const REAL8 eta = m1 * m2 / (Mtot * Mtot);

  // Calculate phenomenological parameters
  // const REAL8 finspin = FinalSpin0815(eta, chi1z, chi2z); //FinalSpin0815 - 0815 is like a version number
  const REAL8 chip = XLALSimPhenomUtilsChiP(m1, m2, chi1x, chi1y, chi2x, chi2y);
  const REAL8 finspin = XLALSimPhenomUtilsPhenomPv2FinalSpin(m1, m2, chi1z, chi2z, chip);
  // const REAL8 finspin = XLALSimPhenomUtilsPhenomPv3HMFinalSpin(m1, m2, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z);

  if (finspin < MIN_FINAL_SPIN)
    XLAL_PRINT_WARNING("Final spin (Mf=%g) and ISCO frequency of this system are small, \
                            the model might misbehave here.",
                       finspin);

  IMRPhenomDAmplitudeCoefficients *pAmp;
  pAmp = XLALMalloc(sizeof(IMRPhenomDAmplitudeCoefficients));
  ComputeIMRPhenomDAmplitudeCoefficients(pAmp, eta, chi1z, chi2z, finspin);
  if (!pAmp)
    XLAL_ERROR(XLAL_EFUNC);

  AmpInsPrefactors amp_prefactors;
  retcode = 0;
  retcode = init_amp_ins_prefactors(&amp_prefactors, pAmp);
  XLAL_CHECK(XLAL_SUCCESS == retcode, retcode, "init_amp_ins_prefactors failed");

  int status_in_for = XLAL_SUCCESS;
/* Now generate the waveform */
#pragma omp parallel for
  for (size_t i = ind_min; i < ind_max; i++)
  {
    REAL8 Mf = freqs->data[i]; // geometric frequency

    UsefulPowers powers_of_f;
    status_in_for = init_useful_powers(&powers_of_f, Mf);
    if (XLAL_SUCCESS != status_in_for)
    {
      XLALPrintError("init_useful_powers failed for Mf, status_in_for=%d", status_in_for);
      retcode = status_in_for;
    }
    else
    {
      amps->data[i] = IMRPhenDAmplitude(Mf, pAmp, &powers_of_f, &amp_prefactors);
    }
  }

  LALFree(pAmp);

  return XLAL_SUCCESS;
}

/**
 * computes the time shift as the approximate time of the peak of the 22 mode.
 */
REAL8 IMRPhenomDComputet0(
    REAL8 eta,           /**< symmetric mass-ratio */
    REAL8 chi1z,         /**< dimensionless aligned-spin of primary */
    REAL8 chi2z,         /**< dimensionless aligned-spin of secondary */
    REAL8 finspin,       /**< final spin */
    LALDict *extraParams /**< linked list containing the extra testing GR parameters */
)
{

  if (extraParams == NULL)
    extraParams = XLALCreateDict();

  IMRPhenomDPhaseCoefficients *pPhi;
  pPhi = XLALMalloc(sizeof(IMRPhenomDPhaseCoefficients));
  ComputeIMRPhenomDPhaseCoefficients(pPhi, eta, chi1z, chi2z, finspin, extraParams);
  if (!pPhi)
    XLAL_ERROR(XLAL_EFUNC);

  IMRPhenomDAmplitudeCoefficients *pAmp;
  pAmp = XLALMalloc(sizeof(IMRPhenomDAmplitudeCoefficients));
  ComputeIMRPhenomDAmplitudeCoefficients(pAmp, eta, chi1z, chi2z, finspin);
  if (!pAmp)
    XLAL_ERROR(XLAL_EFUNC);

  // double Rholm = XLALSimIMRPhenomHMRholm(eta, chi1z, chi2z, ell, mm);
  // double Taulm = XLALSimIMRPhenomHMTaulm(eta, chi1z, chi2z, ell, mm);

  //time shift so that peak amplitude is approximately at t=0
  //For details see https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/WaveformsReview/IMRPhenomDCodeReview/timedomain
  //NOTE: All modes will have the same time offset. So we use the 22 mode.
  //If we just use the 22 mode then we pass 1.0, 1.0 into DPhiMRD.
  const REAL8 t0 = DPhiMRD(pAmp->fmaxCalc, pPhi, 1.0, 1.0);

  LALFree(pPhi);
  LALFree(pAmp);

  return t0;
}

/**
 * Function to compute the amplitude and phase coefficients for PhenomD
 * Used to optimise the calls to IMRPhenDPhase and IMRPhenDAmplitude
 */
int IMRPhenomDSetupAmpAndPhaseCoefficients(
    PhenDAmpAndPhasePreComp *pDPreComp, /**< [out] PhenDAmpAndPhasePreComp struct */
    REAL8 m1,                           /**< mass of companion 1 (Msun) */
    REAL8 m2,                           /**< mass of companion 2 (Msun) */
    REAL8 chi1x,                        /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 chi1y,                        /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 chi1z,                        /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 chi2x,                        /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8 chi2y,                        /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8 chi2z,                        /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    const REAL8 Rholm,                  /**< ratio of (2,2) mode to (l,m) mode ringdown frequency */
    const REAL8 Taulm,                  /**< ratio of (l,m) mode to (2,2) mode damping time */
    LALDict *extraParams                /**<linked list containing the extra parameters */
)
{

  // Make a pointer to LALDict to circumvent a memory leak
  // At the end we will check if we created a LALDict in extraParams
  // and destroy it if we did.
  LALDict *extraParams_in = extraParams;

  /* It's difficult to see in the code but you need to setup the
     * powers_of_pi.
     */
  int retcode = 0;
  retcode = init_useful_powers(&powers_of_pi, LAL_PI);
  XLAL_CHECK(XLAL_SUCCESS == retcode, retcode, "Failed to initiate useful powers of pi.");

  PhenomInternal_PrecessingSpinEnforcePrimaryIsm1(&m1, &m2, &chi1x, &chi1y, &chi1z, &chi2x, &chi2y, &chi2z);
  const REAL8 Mtot = m1 + m2;
  const REAL8 eta = m1 * m2 / (Mtot * Mtot);

  // Calculate phenomenological parameters
  // const REAL8 finspin = FinalSpin0815(eta, chi1z, chi2z); //FinalSpin0815 - 0815 is like a version number
  const REAL8 chip = XLALSimPhenomUtilsChiP(m1, m2, chi1x, chi1y, chi2x, chi2y);
  const REAL8 finspin = XLALSimPhenomUtilsPhenomPv2FinalSpin(m1, m2, chi1z, chi2z, chip);
  // const REAL8 finspin = XLALSimPhenomUtilsPhenomPv3HMFinalSpin(m1, m2, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z);

  if (finspin < MIN_FINAL_SPIN)
    XLAL_PRINT_WARNING("Final spin (Mf=%g) and ISCO frequency of this system are small, \
                            the model might misbehave here.",
                       finspin);

  //start phase
  if (extraParams == NULL)
  {
    extraParams = XLALCreateDict();
  }

  XLALSimInspiralWaveformParamsInsertPNSpinOrder(extraParams, LAL_SIM_INSPIRAL_SPIN_ORDER_35PN);

  IMRPhenomDPhaseCoefficients *pPhi;
  pPhi = XLALMalloc(sizeof(IMRPhenomDPhaseCoefficients));
  ComputeIMRPhenomDPhaseCoefficients(pPhi, eta, chi1z, chi2z, finspin, extraParams);
  if (!pPhi)
    XLAL_ERROR(XLAL_EFUNC);
  PNPhasingSeries *pn = NULL;
  XLALSimInspiralTaylorF2AlignedPhasing(&pn, m1, m2, chi1z, chi2z, extraParams);
  if (!pn)
    XLAL_ERROR(XLAL_EFUNC);

  // Subtract 3PN spin-spin term below as this is in LAL's TaylorF2 implementation
  // (LALSimInspiralPNCoefficients.c -> XLALSimInspiralPNPhasing_F2), but
  // was not available when PhenomD was tuned.
  REAL8 testGRcor = 1.0;
  testGRcor += XLALSimInspiralWaveformParamsLookupNonGRDChi6(extraParams);
  pn->v[6] -= (Subtract3PNSS(m1, m2, Mtot, eta, chi1z, chi2z) * pn->v[0]) * testGRcor;

  PhiInsPrefactors phi_prefactors;
  retcode = 0;
  retcode = init_phi_ins_prefactors(&phi_prefactors, pPhi, pn);
  XLAL_CHECK(XLAL_SUCCESS == retcode, retcode, "init_phi_ins_prefactors failed");

  // Compute coefficients to make phase C^1 continuous (phase and first derivative)
  ComputeIMRPhenDPhaseConnectionCoefficients(pPhi, pn, &phi_prefactors, Rholm, Taulm);
  //end phase

  //start amp
  IMRPhenomDAmplitudeCoefficients *pAmp;
  pAmp = XLALMalloc(sizeof(IMRPhenomDAmplitudeCoefficients));
  ComputeIMRPhenomDAmplitudeCoefficients(pAmp, eta, chi1z, chi2z, finspin);
  if (!pAmp)
    XLAL_ERROR(XLAL_EFUNC);

  AmpInsPrefactors amp_prefactors;
  retcode = 0;
  retcode = init_amp_ins_prefactors(&amp_prefactors, pAmp);
  XLAL_CHECK(XLAL_SUCCESS == retcode, retcode, "init_amp_ins_prefactors failed");
  //end amp

  //output
  pDPreComp->pn = *pn;
  pDPreComp->pPhi = *pPhi;
  pDPreComp->phi_prefactors = phi_prefactors;

  pDPreComp->pAmp = *pAmp;
  pDPreComp->amp_prefactors = amp_prefactors;

  LALFree(pn);
  LALFree(pPhi);
  LALFree(pAmp);

  /* If extraParams was allocated in this function and not passed in
  * we need to free it to prevent a leak */
  if (extraParams && !extraParams_in) {
    XLALDestroyDict(extraParams);
  } else {
    XLALSimInspiralWaveformParamsInsertPNSpinOrder(extraParams,LAL_SIM_INSPIRAL_SPIN_ORDER_ALL);
  }

  return XLAL_SUCCESS;
}

/**
 * Function to return the phenomD phase using the
 * IMRPhenomDSetupAmpAndPhaseCoefficients struct
 */
UNUSED REAL8 IMRPhenomDPhase_OneFrequency(
    REAL8 Mf,
    PhenDAmpAndPhasePreComp pD,
    REAL8 Rholm,
    REAL8 Taulm)
{

  UsefulPowers powers_of_f;
  int status = init_useful_powers(&powers_of_f, Mf);
  XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initiate init_useful_powers");
  REAL8 phase = IMRPhenDPhase(Mf, &(pD.pPhi), &(pD.pn), &powers_of_f,
                              &(pD.phi_prefactors), Rholm, Taulm);
  return phase;
}
