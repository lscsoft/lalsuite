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

#include <gsl/gsl_spline.h>

#include <lal/Units.h>
#include <lal/SeqFactories.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/FrequencySeries.h>
#include <lal/Sequence.h>
#include <lal/LALSimIMR.h>

#include "LALSimIMRSEOBNRv4ROM_NRTidal.h"


void Self_spin_phase_contributions(
  const REAL8 m1_SI,     /**< Mass of neutron star 1 (kg) */
  const REAL8 m2_SI,     /**< Mass of neutron star 2 (kg) */
  const REAL8 chi1L,     /**< Dimensionless aligned component spin of NS 1 */
  const REAL8 chi2L,     /**< Dimensionless aligned component spin of NS 2 */
  const REAL8 qm_def1,   /**< Quadrupole deformation parameter of body 1 (dimensionless) */
  const REAL8 qm_def2,   /**< Quadrupole deformation parameter of body 2 (dimensionless) */
                         /**< qm_def1,2 = 1 for BH */
  REAL8 *pfa_v4_contrib, /**< self-spin contribution to v^4 */
  REAL8 *pfa_v6_contrib  /**< self-spin contribution to v^6 */
) {
  const REAL8 mtot = m1_SI + m2_SI;
  const REAL8 eta = m1_SI*m2_SI/mtot/mtot;
  const REAL8 m1M = m1_SI/mtot;
  const REAL8 m2M = m2_SI/mtot;

  const REAL8 chi1sq = chi1L*chi1L;
  const REAL8 chi2sq = chi2L*chi2L;
  const REAL8 chi1dotchi2 = chi1L*chi2L;

  // From LALSimInspiralPNCoefficients.c:XLALSimInspiralPNPhasing_F2()
  /* Compute 2.0PN SS, QM, and self-spin */
  // See Eq. (6.24) in arXiv:0810.5336
  // 9b,c,d in arXiv:astro-ph/0504538
  REAL8 pn_sigma = eta * (721.L/48.L*chi1L*chi2L - 247.L/48.L*chi1dotchi2);
  pn_sigma += (720.L*qm_def1 - 1.L)/96.0L * m1M * m1M * chi1L * chi1L;
  pn_sigma += (720.L*qm_def2 - 1.L)/96.0L * m2M * m2M * chi2L * chi2L;
  pn_sigma -= (240.L*qm_def1 - 7.L)/96.0L * m1M * m1M * chi1sq;
  pn_sigma -= (240.L*qm_def2 - 7.L)/96.0L * m2M * m2M * chi2sq;
  
  REAL8 pn_ss3 =  (326.75L/1.12L + 557.5L/1.8L*eta)*eta*chi1L*chi2L;
  pn_ss3 += ((4703.5L/8.4L+2935.L/6.L*m1M-120.L*m1M*m1M)*qm_def1 + (-4108.25L/6.72L-108.5L/1.2L*m1M+125.5L/3.6L*m1M*m1M)) *m1M*m1M * chi1sq;
  pn_ss3 += ((4703.5L/8.4L+2935.L/6.L*m2M-120.L*m2M*m2M)*qm_def2 + (-4108.25L/6.72L-108.5L/1.2L*m2M+125.5L/3.6L*m2M*m2M)) *m2M*m2M * chi2sq;

  const REAL8 pfaN = 3.L/(128.L * eta);
  // The leading order term pfa->v[0] is positive and so the 
  // self-spin corrections should be added to a postive phasing.
  *pfa_v4_contrib = -10.L * pn_sigma * pfaN;
  *pfa_v6_contrib = + pn_ss3 * pfaN;
}

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
  REAL8 qm_def1,                                /**< Spin-induced quadrupole moment of mass 1 */
  REAL8 qm_def2,                                /**< Spin-induced quadrupole moment of mass 2 */
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
  if (deltaF > 0)
    ret = XLALSimIMRSEOBNRv4ROM(
      hptilde, hctilde,
      phiRef, deltaF, fLow, fHigh, fRef, distance, inclination,
      m1_SI, m2_SI,
      chi1, chi2,
      -1);
  else
    ret = XLALSimIMRSEOBNRv4ROMFrequencySequence(
      hptilde, hctilde,
      freqs_in,
      phiRef, fRef, distance, inclination,
      m1_SI, m2_SI,
      chi1, chi2,
      -1);
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

  // Tidal self-spin contributions to the phase
  REAL8 pfa_v4_contrib, pfa_v6_contrib;
  Self_spin_phase_contributions(m1_SI, m2_SI, chi1, chi2, qm_def1 - 1.0, qm_def2 - 1.0,
    &pfa_v4_contrib, &pfa_v6_contrib);

  const REAL8 m1 = m1_SI / LAL_MSUN_SI;
  const REAL8 m2 = m2_SI / LAL_MSUN_SI;
  const REAL8 mtot = m1 + m2;
  const REAL8 m_sec = mtot * LAL_MTSUN_SI;  /* total mass in seconds */
  const REAL8 piM = LAL_PI * m_sec;

  // Varibales for phase spline for aligning in time below
  gsl_interp_accel *acc_phi = gsl_interp_accel_alloc();
  gsl_spline *spline_phi = gsl_spline_alloc(gsl_interp_cspline, freqs->length);
  gsl_vector *f_vec = gsl_vector_alloc(freqs->length);
  gsl_vector *phi_vec = gsl_vector_alloc(freqs->length);

  // Assemble waveform from amplitude and phase
  for (size_t i=0; i<freqs->length; i++) { // loop over frequency points in sequence
    int j = i + offset; // shift index for frequency series if needed
    const REAL8 v = cbrt(piM * freqs->data[i]);
    // phasing = (ss_term_v4 * v^4 + ss_term_v6 * v^6) / v^5
    const REAL8 phi_ss = pfa_v4_contrib / v + pfa_v6_contrib * v;
    const REAL8 phase_corr = phi_tidal->data[i] + phi_ss;
    gsl_vector_set(f_vec, i, freqs->data[i]);
    gsl_vector_set(phi_vec, i, phase_corr);

    COMPLEX16 Corr = amp_tidal->data[i] * cexp(-I*phase_corr);
    pdata[j] *= Corr;
    cdata[j] *= Corr;
  }

/* Correct phasing so we coalesce at t=0 (with the definition of the epoch=-1/deltaF above) */
  // Appendix A of 1512.02248
  
  // Get SEOBNRv4 ringdown frequency for 22 mode
  // Note: IMRPhenomPv2_NRTidal also uses the BBH ringdown frequency and then just sets it
  // to the last frequency in the grid
  double fHz_final = XLALSimInspiralGetFinalFreq(m1_SI, m2_SI, 0, 0, chi1, 0, 0, chi2, SEOBNRv4);

  gsl_spline_init(spline_phi, gsl_vector_const_ptr(f_vec, 0), gsl_vector_const_ptr(phi_vec, 0), freqs->length);

  // Time correction is t(f_final) = 1/(2pi) dphi/df (f_final)
  // Here we only apply this to phase corrections beyond SEOBNRv4_ROM.
  // For SEOBNRv4_ROM the phase has already been corrected.

  if (fHz_final > freqs->data[freqs->length-1])
    fHz_final = freqs->data[freqs->length-1];
  REAL8 t_corr_s = gsl_spline_eval_deriv(spline_phi, fHz_final, acc_phi) / (2*LAL_PI);

  // Now correct phase
  for (UINT4 i=0; i<freqs->length; i++) { // loop over frequency points in sequence
    double fHz = freqs->data[i] - fRef;
    int j = i + offset; // shift index for frequency series if needed
    double phase_factor = -2*LAL_PI * fHz * t_corr_s;
    COMPLEX16 t_factor = (cos(phase_factor) + I*sin(phase_factor));
    pdata[j] *= t_factor;
    cdata[j] *= t_factor;
  }

  gsl_vector_free(f_vec);
  gsl_vector_free(phi_vec);
  gsl_spline_free(spline_phi);
  gsl_interp_accel_free(acc_phi);

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
  REAL8 lambda2,                                /**< Dimensionless tidal deformability of NS 2 */
  REAL8 qm_def1,                                /**< Spin-induced quadrupole moment of mass 1 */
  REAL8 qm_def2)                                /**< Spin-induced quadrupole moment of mass 2 */
{
  if (!freqs) XLAL_ERROR(XLAL_EFAULT);

  // Call the internal core function with deltaF = 0 to indicate that freqs is non-uniformly
  // spaced and we want the strain only at these frequencies
  int retcode = SEOBNRv4ROM_NRTidal_Core(hptilde, hctilde,
            phiRef, fRef, distance, inclination, m1_SI, m2_SI, chi1, chi2, 
            lambda1, lambda2, qm_def1, qm_def2, freqs, 0);

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
  REAL8 lambda2,                                /**< Dimensionless tidal deformability of NS 2 */
  REAL8 qm_def1,                                /**< Spin-induced quadrupole moment of mass 1 */
  REAL8 qm_def2                                 /**< Spin-induced quadrupole moment of mass 2 */
) {
  // Use fLow, fHigh, deltaF to compute freqs sequence
  // Instead of building a full sequence we only transfer the boundaries and let
  // the internal core function do the rest (and properly take care of corner cases).
  REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
  freqs->data[0] = fLow;
  freqs->data[1] = fHigh;

  int retcode = SEOBNRv4ROM_NRTidal_Core(hptilde, hctilde,
            phiRef, fRef, distance, inclination, m1_SI, m2_SI, chi1, chi2, 
            lambda1, lambda2, qm_def1, qm_def2, freqs, deltaF);

  XLALDestroyREAL8Sequence(freqs);

  return(retcode);
}

/** @} */
/** @} */
