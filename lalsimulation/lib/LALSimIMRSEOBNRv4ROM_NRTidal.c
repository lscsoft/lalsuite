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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
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

/**
 * Function to internally add 2PN and 3PN spin-spin terms 
 * to be able to include spin-induced quadrupole moments
 * in those terms; the BBH terms are excluded
 * From LALSimInspiralPNCoefficients.c:XLALSimInspiralPNPhasing_F2()
 * Compute 2.0PN SS, QM, and self-spin
 * See Eq. (6.24) in arXiv:0810.5336 
 * 9b,c,d in arXiv:astro-ph/0504538
 */
void Self_spin_phase_contributions(
				   const REAL8 m1_SI,     /**< Mass of neutron star 1 (kg) */
				   const REAL8 m2_SI,     /**< Mass of neutron star 2 (kg) */
				   const REAL8 chi1L,     /**< Dimensionless aligned component spin of NS 1 */
				   const REAL8 chi2L,     /**< Dimensionless aligned component spin of NS 2 */
				   const REAL8 qm_def1,   /**< Quadrupole deformation parameter of body 1 (dimensionless) */
				   const REAL8 qm_def2,   /**< Quadrupole deformation parameter of body 2 (dimensionless) */
				   /**< qm_def1,2 = 0 for BH as it is defined here*/
				   REAL8 *pfa_v4_contrib, /**< self-spin contribution to v^4 */
				   REAL8 *pfa_v6_contrib  /**< self-spin contribution to v^6 */
				   ) {
  const REAL8 mtot = m1_SI + m2_SI;
  const REAL8 eta = m1_SI*m2_SI/mtot/mtot;
  const REAL8 m1M = m1_SI/mtot;
  const REAL8 m2M = m2_SI/mtot;

  const REAL8 m1Msq = m1M * m1M;
  const REAL8 m2Msq = m2M * m2M;

  const REAL8 chi1sq = chi1L*chi1L;
  const REAL8 chi2sq = chi2L*chi2L;

  /* remove unnecessary calls and computations to speed up the final computation */
  REAL8 pn_sigma = - 50.L*(qm_def1 * chi1sq * m1Msq + qm_def2 * chi2sq * m2Msq);
  REAL8 pn_ss3   = 5.L/84.L*(9407.L+ 8218.L * m1M - 2016.L * m1Msq) * qm_def1 * m1Msq * chi1sq;
  pn_ss3 += 5.L/84.L*(9407.L+ 8218.L * m2M - 2016.L * m2Msq) * qm_def2 * m2Msq * chi2sq;

  const REAL8 pfaN = 3.L/(128.L * eta);
  // The leading order term pfa->v[0] is positive and so the 
  // self-spin corrections should be added to a postive phasing.


  //  *pfa_v4_contrib = -10.L * pn_sigma * pfaN;
  //  *pfa_v6_contrib = + pn_ss3 * pfaN;
  // no additional factor 10 needed in the new implementation
  *pfa_v4_contrib =   pn_sigma * pfaN;
  *pfa_v6_contrib =   pn_ss3 * pfaN;
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
  const REAL8Sequence *freqs_in,                /**< Frequency points at which to evaluate the waveform (Hz) */
  REAL8 deltaF,                                 /**< Sampling frequency (Hz) */
  LALDict *LALparams,                          /**< LAL dictionary containing accessory parameters */
  NRTidal_version_type NRTidal_version         /**< Version of NRTides; can be any one of NRTidal_V (arXiv:1706.02969), NRTidalv2_V (arXiv:1905.06011) or NRTidalv2NoAmpCorr_V (arXiv:1905.06011, without amplitude corrections) */ )
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

  REAL8 quad_mon1 = XLALSimInspiralWaveformParamsLookupdQuadMon1(LALparams);
  REAL8 quad_mon2 = XLALSimInspiralWaveformParamsLookupdQuadMon2(LALparams);

  /* Internally we need m1 > m2, so change around if this is not the case */
  if (m1_SI < m2_SI) {
    // Swap m1 and m2
    double m1temp = m1_SI;
    double chi1temp = chi1;
    double lambda1temp = lambda1;
    double quadmon1temp = quad_mon1;
    m1_SI = m2_SI;
    chi1 = chi2;
    m2_SI = m1temp;
    chi2 = chi1temp;
    if (lambda1!=lambda2){
      lambda1 = lambda2;
      XLALSimInspiralWaveformParamsInsertTidalLambda1(LALparams, lambda1);
      lambda2 = lambda1temp;
      XLALSimInspiralWaveformParamsInsertTidalLambda2(LALparams, lambda2);
    }
    if (quad_mon1 != quad_mon2) {
      quad_mon1 = quad_mon2;
      quad_mon2 = quadmon1temp;
    }
  }
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
    if ( (NRTidal_version != NRTidalv2NSBH_V) && (( fHigh > NRTIDAL_FMAX ) || ( fHigh == 0.0 )) )
      {
        // only generate upto NRTIDAL_FMAX
        f_max_nr_tidal = NRTIDAL_FMAX;
      }

    ret = XLALSimIMRSEOBNRv4ROM(
				hptilde, hctilde,
				phiRef, deltaF, fLow, f_max_nr_tidal, fRef, distance, inclination,
				m1_SI, m2_SI,
				chi1, chi2,
				-1, LALparams, NRTidal_version);

    // if uniform sampling and fHigh > NRTIDAL_FMAX then resize htilde
    // so that it goes up to the user fHigh but is filled with zeros
    // beyond NRTIDAL_FMAX (this does not apply to NSBH)
    if (fHigh > NRTIDAL_FMAX && NRTidal_version != NRTidalv2NSBH_V) 
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
						 -1, LALparams, NRTidal_version);
  }
  XLAL_CHECK(XLAL_SUCCESS == ret, ret, "XLALSimIMRSEOBNRv4ROM() failed.");

  UINT4 offset;
  REAL8Sequence *freqs = NULL;
  REAL8Sequence *phi_tidal = NULL;
  REAL8Sequence *amp_tidal = NULL;
  REAL8Sequence *planck_taper = NULL;
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

  const REAL8 m1 = m1_SI / LAL_MSUN_SI;
  const REAL8 m2 = m2_SI / LAL_MSUN_SI;
  const REAL8 mtot = m1 + m2;
  const REAL8 m_sec = mtot * LAL_MTSUN_SI;  /* total mass in seconds */
  const REAL8 piM = LAL_PI * m_sec;
  /* Initialising parameters for adding higher order spin corrections */
  REAL8 X_A = m1/mtot;
  REAL8 X_B = m2/mtot;
  REAL8 eta = m1 * m2 / (mtot* mtot);    /* Symmetric mass-ratio */
  REAL8 pn_fac = 3./(128.*eta);
  REAL8 SS_3p5PN = 0., SSS_3p5PN = 0.;
  /* End of initialising */

  // Get FD tidal phase correction and amplitude factor from arXiv:1706.02969
  phi_tidal = XLALCreateREAL8Sequence(freqs->length);
  planck_taper = XLALCreateREAL8Sequence(freqs->length);
  if (NRTidal_version == NRTidalv2_V) {
    ret = XLALSimNRTunedTidesFDTidalPhaseFrequencySeries(phi_tidal, amp_tidal, planck_taper, freqs, m1_SI, m2_SI, lambda1, lambda2, NRTidalv2NoAmpCorr_V);
    XLAL_CHECK(XLAL_SUCCESS == ret, ret, "XLALSimNRTunedTidesFDTidalPhaseFrequencySeries Failed.");
    XLALSimInspiralGetHOSpinTerms(&SS_3p5PN, &SSS_3p5PN, X_A, X_B, chi1, chi2, quad_mon1+1., quad_mon2+1.);
  }
  else {
    ret = XLALSimNRTunedTidesFDTidalPhaseFrequencySeries(phi_tidal, amp_tidal, planck_taper, freqs, m1_SI, m2_SI, lambda1, lambda2, NRTidal_version);
    XLAL_CHECK(XLAL_SUCCESS == ret, ret, "XLALSimNRTunedTidesFDTidalPhaseFrequencySeries Failed.");
  }

  // For NSBH, apply NSBH amplitude correction
  if (NRTidal_version == NRTidalv2NSBH_V){
       amp_tidal = XLALCreateREAL8Sequence(freqs->length);
       ret = XLALSEOBNRv4ROMNSBHAmplitudeCorrectionFrequencySeries(
        amp_tidal, freqs,
        m1_SI, m2_SI, chi1, lambda2
      );
      XLAL_CHECK(XLAL_SUCCESS == ret, ret, "XLALSEOBNRv4ROMNSBHAmplitudeCorrectionFrequencySeries Failed.");

  }

  // // Prepare tapering of amplitude beyond merger frequency
  // double kappa2T = XLALSimNRTunedTidesComputeKappa2T(m1_SI, m2_SI, lambda1, lambda2);
  // double fHz_mrg = XLALSimNRTunedTidesMergerFrequency(Mtot_MSUN, kappa2T, q);
  // Tidal self-spin contributions to the phase

  REAL8 pfa_v4_contrib, pfa_v6_contrib;
  Self_spin_phase_contributions(m1_SI, m2_SI, chi1, chi2, quad_mon1, quad_mon2,
				&pfa_v4_contrib, &pfa_v6_contrib);

  gsl_interp_accel *acc_phi = gsl_interp_accel_alloc();
  gsl_spline *spline_phi = gsl_spline_alloc(gsl_interp_cspline, freqs->length);
  gsl_vector *f_vec = gsl_vector_alloc(freqs->length);
  gsl_vector *phi_vec = gsl_vector_alloc(freqs->length);

  // Assemble waveform from amplitude and phase
  for (size_t i=0; i<freqs->length; i++) { // loop over frequency points in sequence
    int j = i + offset; // shift index for frequency series if needed
    // Apply tidal phase correction and amplitude taper
    // double taper = 1.0 - PlanckTaper(freqs->data[i], fHz_mrg, 1.2*fHz_mrg);
    const REAL8 v = cbrt(piM * freqs->data[i]);
    // phasing = (ss_term_v4 * v^4 + ss_term_v6 * v^6) / v^5
    const REAL8 phi_ss = pfa_v4_contrib / v + pfa_v6_contrib * v;
    const REAL8 phase_corr = phi_tidal->data[i] + phi_ss;
    gsl_vector_set(f_vec, i, freqs->data[i]);
    gsl_vector_set(phi_vec, i, phase_corr);

    COMPLEX16 Corr = planck_taper->data[i] * cexp(-I*phase_corr -I*v*v*(SS_3p5PN + SSS_3p5PN)*pn_fac);
    if (NRTidal_version==NRTidalv2NSBH_V){
        Corr *= amp_tidal->data[i];
    }
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
  // From Eqn. (A1) of arXiv:1512.02248
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
  XLALDestroyREAL8Sequence(planck_taper);
  

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
  LALDict *LALparams,                           /**< linked list containing the extra testing GR parameters */
  NRTidal_version_type NRTidal_version          /**< Version of NRTides; can be any one of NRTidal_V (arXiv:1706.02969), NRTidalv2_V (arXiv:1905.06011) or NRTidalv2NoAmpCorr_V (arXiv:1905.06011, without amplitude corrections) */ )
{
  if (!freqs) XLAL_ERROR(XLAL_EFAULT);

  // Call the internal core function with deltaF = 0 to indicate that freqs is non-uniformly
  // spaced and we want the strain only at these frequencies
  int retcode = SEOBNRv4ROM_NRTidal_Core(hptilde, hctilde,
					 phiRef, fRef, distance, inclination, m1_SI, m2_SI, chi1, chi2, lambda1, lambda2, freqs, 0, LALparams, NRTidal_version);

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
				 REAL8 lambda2,                                 /**< Dimensionless tidal deformability of NS 2 */
				 LALDict *LALparams,                            /**< linked list containing the extra testing GR parameters */
				 NRTidal_version_type NRTidal_version           /**< Version of NRTides; can be one of NRTidal or NRTidalv2NoAmpCorr */
				 ) {
  // Use fLow, fHigh, deltaF to compute freqs sequence
  // Instead of building a full sequence we only transfer the boundaries and let
  // the internal core function do the rest (and properly take care of corner cases).
  REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
  freqs->data[0] = fLow;
  freqs->data[1] = fHigh;

  int retcode = SEOBNRv4ROM_NRTidal_Core(hptilde, hctilde,
					 phiRef, fRef, distance, inclination, m1_SI, m2_SI, chi1, chi2, lambda1, lambda2, freqs, deltaF, LALparams, NRTidal_version);

  XLALDestroyREAL8Sequence(freqs);

  return(retcode);
}

/** @} */
/** @} */
