/*
 *  Copyright (C) 2017 Michael Puerrer
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

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include <lal/Units.h>
#include <lal/SeqFactories.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/FrequencySeries.h>
#include <lal/Sequence.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMRSEOBNRv4ROM_NRTidal.h"


int SEOBNRv4ROM_NRTidal_Core(
  struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
  struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
  REAL8 phiRef,                                 /**< Phase at reference time */
  REAL8 fRef,                                   /**< Reference frequency (Hz); 0 defaults to fLow */
  REAL8 distance,                               /**< Distance of source (m) */
  REAL8 inclination,                            /**< Inclination of source (rad) */
  REAL8 m1_SI,                                  /**< Mass of neutron star 1 (kg) */
  REAL8 m2_SI,                                  /**< Mass of neutron star 2 (kg) */
  REAL8 chi1,                                   /**< Dimensionless aligned component spin of NS 1 */
  REAL8 chi2,                                   /**< Dimensionless aligned component spin of NS 2 */
  REAL8 lambda1,                                /**< Dimensionless tidal deformability of NS 1 */
  REAL8 lambda2,                                /**< Dimensionless tidal deformability of NS 2 */
  const REAL8Sequence *freqs_in,                /**< Frequency points at which to evaluate the waveform (Hz) */
  REAL8 deltaF)                                 /**< Sampling frequency (Hz) */
{
  /* Check output arrays */
  if(!hptilde || !hctilde)
    XLAL_ERROR(XLAL_EFAULT);
  if(*hptilde || *hctilde) {
    XLALPrintError("(*hptilde) and (*hctilde) are supposed to be NULL, but got %p and %p",(*hptilde),(*hctilde));
    XLAL_ERROR(XLAL_EFAULT);
  }

  if (!freqs_in) XLAL_ERROR(XLAL_EFAULT);
  double fLow  = freqs_in->data[0];
  double fHigh = freqs_in->data[freqs_in->length - 1];
  if(fRef == 0.0)
    fRef = fLow;

  /* Internally we need m1 > m2, so change around if this is not the case */
  if (m1_SI < m2_SI) {
    // Swap m1 and m2
    double m1temp = m1_SI;
    double chi1temp = chi1;
    double lambda1temp = lambda1;
    m1_SI = m2_SI;
    chi1 = chi2;
    lambda1 = lambda2;
    m2_SI = m1temp;
    chi2 = chi1temp;
    lambda2 = lambda1temp;
  }

  // double Mtot_MSUN = (m1_SI + m2_SI) / LAL_MSUN_SI;
  // double q = m1_SI / m2_SI;

  // Call SEOBNRv4 ROM. We call either the FrequencySequence version
  // or the regular LAL version depending on how we've been called.

  // These functions enforce m1 >= m2:
  // XLALSimIMRSEOBNRv4ROM / XLALSimIMRSEOBNRv4ROMFrequencySequence
  // XLALSimNRTunedTidesFDTidalPhaseFrequencySeries

  int ret = XLAL_SUCCESS;
  if (deltaF > 0) {
    // if using a uniform frequency series then we only need to generate
    // SEOBNRv4_ROM upto a bit beyond the BNS merger frequency.
    // if asked for a frequency beyond NRTIDAL_FMAX then the
    // returned waveform contains frequencies up to the input fHigh but
    // only contains zeros beyond NRTIDAL_FMAX
    double f_max_nr_tidal = fHigh;
    /**< tidal coupling constant.*/
    const double kappa2T = XLALSimNRTunedTidesComputeKappa2T(m1_SI, m2_SI, lambda1, lambda2);
    /* Prepare tapering of amplitude beyond merger frequency */
    const double fHz_mrg = XLALSimNRTunedTidesMergerFrequency( (m1_SI+m2_SI)/LAL_MSUN_SI , kappa2T, m1_SI/m2_SI);
    const double NRTIDAL_FMAX = 1.3*fHz_mrg;

    if ( ( fHigh > NRTIDAL_FMAX ) || ( fHigh == 0.0 ) )
    {
        // only generate upto NRTIDAL_FMAX
        f_max_nr_tidal = NRTIDAL_FMAX;
    }

    ret = XLALSimIMRSEOBNRv4ROM(
      hptilde, hctilde,
      phiRef, deltaF, fLow, f_max_nr_tidal, fRef, distance, inclination,
      m1_SI, m2_SI,
      chi1, chi2,
      -1);

      // if uniform sampling and fHigh > NRTIDAL_FMAX then resize htilde
      // so that it goes up to the user fHigh but is filled with zeros
      // beyond NRTIDAL_FMAX
      if (fHigh > NRTIDAL_FMAX)
      {
          // resize
          // n_full is the next power of 2 +1.
          size_t n_full = (size_t) pow(2,ceil(log2(fHigh / deltaF))) + 1;
          *hptilde = XLALResizeCOMPLEX16FrequencySeries(*hptilde, 0, n_full);
          XLAL_CHECK ( *hptilde, XLAL_ENOMEM, "Failed to resize hptilde COMPLEX16FrequencySeries");
          *hctilde = XLALResizeCOMPLEX16FrequencySeries(*hctilde, 0, n_full);
          XLAL_CHECK ( *hctilde, XLAL_ENOMEM, "Failed to resize hctilde COMPLEX16FrequencySeries");
      }

  } else {
    ret = XLALSimIMRSEOBNRv4ROMFrequencySequence(
      hptilde, hctilde,
      freqs_in,
      phiRef, fRef, distance, inclination,
      m1_SI, m2_SI,
      chi1, chi2,
      -1);
  }
  XLAL_CHECK(XLAL_SUCCESS == ret, ret, "XLALSimIMRSEOBNRv4ROM() failed.");

  UINT4 offset;
  REAL8Sequence *freqs = NULL;
  if (deltaF > 0) { // uniform frequencies
    // Recreate freqs using only the lower and upper bounds
    UINT4 iStart = (UINT4) ceil(fLow / deltaF);
    UINT4 iStop = (*hptilde)->data->length - 1; // use the length calculated in the ROM function
    freqs = XLALCreateREAL8Sequence(iStop - iStart);
    if (!freqs) XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
    for (UINT4 i=iStart; i<iStop; i++)
      freqs->data[i-iStart] = i*deltaF;

    offset = iStart;
  }
  else { // unequally spaced frequency sequence
    freqs = XLALCreateREAL8Sequence(freqs_in->length);
    if (!freqs) XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
    for (UINT4 i=0; i<freqs_in->length; i++)
      freqs->data[i] = freqs_in->data[i]; // just copy input
    offset = 0;
  }
  COMPLEX16 *pdata=(*hptilde)->data->data;
  COMPLEX16 *cdata=(*hctilde)->data->data;

  // Get FD tidal phase correction and amplitude factor from arXiv:1706.02969
  REAL8Sequence *phi_tidal = XLALCreateREAL8Sequence(freqs->length);
  REAL8Sequence *amp_tidal = XLALCreateREAL8Sequence(freqs->length);
  ret = XLALSimNRTunedTidesFDTidalPhaseFrequencySeries(
    phi_tidal, amp_tidal, freqs,
    m1_SI, m2_SI, lambda1, lambda2
  );
  XLAL_CHECK(XLAL_SUCCESS == ret, ret, "XLALSimNRTunedTidesFDTidalPhaseFrequencySeries Failed.");

  // // Prepare tapering of amplitude beyond merger frequency
  // double kappa2T = XLALSimNRTunedTidesComputeKappa2T(m1_SI, m2_SI, lambda1, lambda2);
  // double fHz_mrg = XLALSimNRTunedTidesMergerFrequency(Mtot_MSUN, kappa2T, q);

  // Assemble waveform from amplitude and phase
  for (size_t i=0; i<freqs->length; i++) { // loop over frequency points in sequence
    int j = i + offset; // shift index for frequency series if needed
    // Apply tidal phase correction and amplitude taper
    // double taper = 1.0 - PlanckTaper(freqs->data[i], fHz_mrg, 1.2*fHz_mrg);
    COMPLEX16 Corr = amp_tidal->data[i] * cexp(-I*phi_tidal->data[i]);
    pdata[j] *= Corr;
    cdata[j] *= Corr;
  }

  XLALDestroyREAL8Sequence(freqs);
  XLALDestroyREAL8Sequence(phi_tidal);
  XLALDestroyREAL8Sequence(amp_tidal);

  return XLAL_SUCCESS;
}

// FIXME: limits below
/**
 * @addtogroup LALSimIMRTIDAL_c
 *
 * @{
 *
 * @name SEOBNRv4_ROM_NRTidal
 *
 * @author Michael Puerrer
 *
 * @brief C code for SEOBNRv4ROM arXiv:1611.03703 with added tidal phase correction from arXiv:1706.02969.
 *
 * This is a frequency domain model that adds tidal modifications of the phasing
 * to the SEOBNRv4ROM model.
 *
 * @note Parameter ranges:
 *   * ? <= eta <= 0.25
 *   * 0 <= Lambda_i <= ?
 *   * -1 <= chi_i <= 1
 *   * Mtot >= 2 Msun @ 20 Hz (inherited from the ROM)
 *
 *  Aligned component spin on neutron stars.
 *  Symmetric mass-ratio eta = m1*m2/(m1+m2)^2.
 *  Total mass Mtot.
 *
 * @{
 */

/**
 * Compute waveform in LAL format at specified frequencies for the SEOBNRv4_ROM_NRTidal
 * tidal model based on SEOBNRv4_ROM.
 *
 * XLALSimIMRSEOBNRv4ROMNRTidal() returns the plus and cross polarizations as a complex
 * frequency series with equal spacing deltaF and contains zeros from zero frequency
 * to the starting frequency and zeros beyond the cutoff frequency in the ringdown.
 *
 * In contrast, XLALSimIMRSEOBNRv4ROMNRTidalFrequencySequence() returns a
 * complex frequency series with entries exactly at the frequencies specified in
 * the sequence freqs (which can be unequally spaced). No zeros are added.
 *
 * If XLALSimIMRSEOBNRv4ROMNRTidalFrequencySequence() is called with frequencies that
 * are beyond the maxium allowed geometric frequency for the ROM, zero strain is returned.
 * It is not assumed that the frequency sequence is ordered.
 *
 * This function is designed as an entry point for reduced order quadratures.
 */
int XLALSimIMRSEOBNRv4ROMNRTidalFrequencySequence(
  struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
  struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
  const REAL8Sequence *freqs,                   /**< Frequency points at which to evaluate the waveform (Hz) */
  REAL8 phiRef,                                 /**< Phase at reference time */
  REAL8 fRef,                                   /**< Reference frequency (Hz); 0 defaults to fLow */
  REAL8 distance,                               /**< Distance of source (m) */
  REAL8 inclination,                            /**< Inclination of source (rad) */
  REAL8 m1_SI,                                  /**< Mass of neutron star 1 (kg) */
  REAL8 m2_SI,                                  /**< Mass of neutron star 2 (kg) */
  REAL8 chi1,                                   /**< Dimensionless aligned component spin of NS 1 */
  REAL8 chi2,                                   /**< Dimensionless aligned component spin of NS 2 */
  REAL8 lambda1,                                /**< Dimensionless tidal deformability of NS 1 */
  REAL8 lambda2)                                /**< Dimensionless tidal deformability of NS 2 */
{
  if (!freqs) XLAL_ERROR(XLAL_EFAULT);

  // Call the internal core function with deltaF = 0 to indicate that freqs is non-uniformly
  // spaced and we want the strain only at these frequencies
  int retcode = SEOBNRv4ROM_NRTidal_Core(hptilde, hctilde,
            phiRef, fRef, distance, inclination, m1_SI, m2_SI, chi1, chi2, lambda1, lambda2, freqs, 0);

  return(retcode);
}

/**
 * Compute waveform in LAL format for the SEOBNRv4_ROM_NRTidal
 * tidal model based on SEOBNRv4_ROM.
 *
 * Returns the plus and cross polarizations as a complex frequency series with
 * equal spacing deltaF and contains zeros from zero frequency to the starting
 * frequency fLow and zeros beyond the cutoff frequency in the ringdown.
 */
int XLALSimIMRSEOBNRv4ROMNRTidal(
  struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
  struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
  REAL8 phiRef,                                 /**< Phase at reference time */
  REAL8 deltaF,                                 /**< Sampling frequency (Hz) */
  REAL8 fLow,                                   /**< Starting GW frequency (Hz) */
  REAL8 fHigh,                                  /**< End frequency; 0 defaults to Mf=0.14 */
  REAL8 fRef,                                   /**< Reference frequency (Hz); 0 defaults to fLow */
  REAL8 distance,                               /**< Distance of source (m) */
  REAL8 inclination,                            /**< Inclination of source (rad) */
  REAL8 m1_SI,                                  /**< Mass of neutron star 1 (kg) */
  REAL8 m2_SI,                                  /**< Mass of neutron star 2 (kg) */
  REAL8 chi1,                                   /**< Dimensionless aligned component spin of NS 1 */
  REAL8 chi2,                                   /**< Dimensionless aligned component spin of NS 2 */
  REAL8 lambda1,                                /**< Dimensionless tidal deformability of NS 1 */
  REAL8 lambda2                                 /**< Dimensionless tidal deformability of NS 2 */
) {
  // Use fLow, fHigh, deltaF to compute freqs sequence
  // Instead of building a full sequence we only transfer the boundaries and let
  // the internal core function do the rest (and properly take care of corner cases).
  REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
  freqs->data[0] = fLow;
  freqs->data[1] = fHigh;

  int retcode = SEOBNRv4ROM_NRTidal_Core(hptilde, hctilde,
            phiRef, fRef, distance, inclination, m1_SI, m2_SI, chi1, chi2, lambda1, lambda2, freqs, deltaF);

  XLALDestroyREAL8Sequence(freqs);

  return(retcode);
}

/** @} */
/** @} */
