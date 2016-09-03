/*
 *  Copyright (C) 2016 Michael Puerrer, Prayush Kumar
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

#include "LALSimIMRLackeyTidal2013.h"


/*************** Model coefficients ******************/

// Amplitude correction factors for SEOBNRv2
// Define constants as per Eq 33 of Lackey et al, arXiv:1303.6298.
const double b0 = -1424.2;
const double b1 = 6423.4;
const double b2 = 0.84203;
const double c0 = -9.7628;
const double c1 = 33.939;
const double c2 = 1.0971;
// Phase correction factors
const double g0 = -4.6339;
const double g1 = 27.719;
const double g2 = 10.268;
const double g3 = -41.741;


/********************* Definitions begin here ********************/

static void tidalPNAmplitudeCoefficient(
  double *C,
  const double eta,
  const double chi_BH,
  const double Lambda
) {
  // Coefficient in the amplitude factor, Eq 33 of Lackey et al
  *C = exp(b0 + b1*eta + b2*chi_BH)
     + Lambda * exp(c0 + c1*eta + c2*chi_BH);
}

static double tidalCorrectionAmplitude(
  const double Mf,
  const double C,
  const double eta,
  const double Lambda
) {
  const double MfA = 0.01; // amplitude transition frequency
  if (Mf <= MfA)
    return 1.0;
  else {
    // Generate the amplitude factor, Eq 33 of Lackey et al
    double dMf = Mf - MfA;
    double dMf2 = dMf*dMf;
    double B = C * dMf*dMf2;
    return exp(-eta * Lambda * B);
  }
}

// precompute a0, a1 and G which do not depend on frequency
static void tidalPNPhaseCoefficients(
  double *a0,
  double *a1,
  double *G,
  const double eta,
  const double chi_BH,
  const double Lambda
) {
  // First compute the PN inspiral phasing correction
  // see Eq. 7,8 of Lackey et al
  double eta2 = eta*eta;
  double eta3 = eta2*eta;
  double SqrtOneMinus4Eta = sqrt(1.-4.*eta);

  *a0 = -12 * Lambda * ((1 + 7.*eta - 31*eta2)
      - SqrtOneMinus4Eta * (1 + 9.*eta - 11*eta2));
  *a1 = -(585.*Lambda/28.)
      * ((1. + 3775.*eta/234. - 389.*eta2/6. + 1376.*eta3/117.)
      - SqrtOneMinus4Eta*(1 + 4243.*eta/234. - 6217*eta2/234. - 10.*eta3/9.));

  *G = exp(g0 + g1*eta + g2*chi_BH + g3*eta*chi_BH); // Eq 35 of Lackey et al
}

static double tidalPNPhase(
  const double Mf,
  const double a0,
  const double a1,
  const double eta
) {
  // First compute the PN inspiral phasing correction
  // see Eq. 7,8 of Lackey et al
  double v = cbrt(LAL_PI * Mf);
  double v2 = v*v;
  double v5 = v2*v2*v;
  double v7 = v5*v2;
  return 3.*(a0*v5 + a1*v7) / (128.*eta);
}

static double tidalPNPhaseDeriv(
  const double Mf,
  const double a0,
  const double a1,
  const double eta
) {
  // First compute the PN inspiral phasing correction
  // see Eq. 7,8 of Lackey et al
  double v = cbrt(LAL_PI * Mf);
  double v2 = v*v;
  double v4 = v2*v2;
  return LAL_PI * (5.*a0*v2 + 7.*a1*v4) / (128.*eta);
}

// Implements Eq. 34 of Lackey et al
static double tidalCorrectionPhase(
  const double Mf,
  const double a0,
  const double a1,
  const double G,
  const double eta,
  const double Lambda
)
{
  const double MfP = 0.02; // phase transition frequency

  if (Mf <= MfP)
    return tidalPNPhase(Mf, a0, a1, eta);

  // Beyond the phase transition frequency we evaluate the tidal phase
  // and its derivative at the transition frequency
  double psiT = tidalPNPhase(MfP, a0, a1, eta);
  double DpsiT= (Mf - MfP) * tidalPNPhaseDeriv(MfP, a0, a1, eta);
  // Now compute the phenomenological term
  double E = G * pow(Mf - MfP, 5./3.); // Eq 35 of Lackey et al
  double psiFit = eta * Lambda * E;
  return psiT + DpsiT - psiFit; // Eq 34 of Lackey et al
}

int LackeyTidal2013SEOBNRv2ROMCore(
  struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
  struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
  REAL8 phiRef,                                 /**< Phase at reference time */
  REAL8 fRef,                                   /**< Reference frequency (Hz); 0 defaults to fLow */
  REAL8 distance,                               /**< Distance of source (m) */
  REAL8 inclination,                            /**< Inclination of source (rad) */
  REAL8 mBH_SI,                                 /**< Mass of black hole (kg) */
  REAL8 mNS_SI,                                 /**< Mass of neutron star (kg) */
  REAL8 chi_BH,                                 /**< Dimensionless aligned component spin of the BH */
  REAL8 Lambda,                                 /**< Dimensionless tidal deformability (Eq 1  of Lackey et al) */
  const REAL8Sequence *freqs_in,                /**< Frequency points at which to evaluate the waveform (Hz) */
  REAL8 deltaF                                  /**< Sampling frequency (Hz) */
)
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

  double mBH = mBH_SI / LAL_MSUN_SI;
  double mNS = mNS_SI / LAL_MSUN_SI;
  double M = mBH + mNS;
  double eta = mBH * mNS / (M*M);    /* Symmetric mass-ratio */
  double Mtot_sec = M * LAL_MTSUN_SI; /* Total mass in seconds */
  double chi_NS = 0; // NS has zero spin

  // Impose sanity checks and cutoffs on mass-ratio, and BH spins
  if (mBH < mNS) XLAL_ERROR(XLAL_EDOM, "mBH = %g < mNS = %g ! ", mBH, mNS);
  if (eta < 6./49.) XLAL_ERROR(XLAL_EDOM, "eta = %g < 6/49!", eta);
  if (chi_BH > 0.75) XLAL_ERROR(XLAL_EDOM, "BH spin = %g > 0.75!", chi_BH);
  if (chi_BH < -0.75) XLAL_ERROR(XLAL_EDOM, "BH spin = %g < -0.75!", chi_BH);

  // Call the high-resolution SEOBNRv2 ROM that can go to very low total mass
  // We call either the FrequencySequence version or the regular LAL version depending on how we've been called.
  int ret = XLAL_SUCCESS;
  if (deltaF > 0)
    ret = XLALSimIMRSEOBNRv2ROMDoubleSpinHI(
      hptilde, hctilde,
      phiRef, deltaF, fLow, fHigh, fRef, distance, inclination,
      mBH_SI, mNS_SI,
      chi_BH, chi_NS,
      -1);
  else
    ret = XLALSimIMRSEOBNRv2ROMDoubleSpinHIFrequencySequence(
      hptilde, hctilde,
      freqs_in,
      phiRef, fRef, distance, inclination,
      mBH_SI, mNS_SI,
      chi_BH, chi_NS,
      -1);
  XLAL_CHECK(XLAL_SUCCESS == ret, ret, "XLALSimIMRSEOBNRv2ROMDoubleSpinHI() failed.");

  UINT4 offset;
  REAL8Sequence *freqs = NULL;
  if (deltaF > 0) { // uniform frequencies
    // Recreate freqs using only the lower and upper bounds
    UINT4 iStart = (UINT4) ceil(fLow / deltaF);
    UINT4 iStop = (*hptilde)->data->length - 1; // use the length calculated in the ROM function
    freqs = XLALCreateREAL8Sequence(iStop - iStart);
    if (!freqs) XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
    double deltaF_geom = deltaF * Mtot_sec;
    for (UINT4 i=iStart; i<iStop; i++)
      freqs->data[i-iStart] = i*deltaF_geom;

    offset = iStart;
  }
  else { // unequally spaced frequency sequence
    freqs = XLALCreateREAL8Sequence(freqs_in->length);
    if (!freqs) XLAL_ERROR(XLAL_EFUNC, "Frequency array allocation failed.");
    for (UINT4 i=0; i<freqs_in->length; i++)
      freqs->data[i] = freqs_in->data[i] * Mtot_sec; // just copy input and convert to geometric frequency
    offset = 0;
  }
  COMPLEX16 *pdata=(*hptilde)->data->data;
  COMPLEX16 *cdata=(*hctilde)->data->data;

  // Precompute coefficients that do not depend on frequency
  double C, a0, a1, G;
  tidalPNAmplitudeCoefficient(&C, eta, chi_BH, Lambda);
  tidalPNPhaseCoefficients(&a0, &a1, &G, eta, chi_BH, Lambda);

  // Assemble waveform from aplitude and phase
  for (size_t i=0; i<freqs->length; i++) { // loop over frequency points in sequence
    double Mf = freqs->data[i];
    int j = i + offset; // shift index for frequency series if needed
    // Tidal corrections to be incorporated
    double ampC = tidalCorrectionAmplitude(Mf, C, eta, Lambda);
    double phsC = tidalCorrectionPhase(Mf, a0, a1, G, eta, Lambda);
    COMPLEX16 Corr = ampC * cexp(-I*phsC);
    pdata[j] *= Corr;
    cdata[j] *= Corr;
  }

  XLALDestroyREAL8Sequence(freqs);

  return XLAL_SUCCESS;
}

/**
 * @addtogroup LALSimIMRTIDAL_c
 *
 * @{
 *
 * @name Lackey et al (2013) tidal model based on SEOBNRv2_ROM
 *
 * @author Michael Puerrer, Prayush Kumar
 *
 * @brief C code for Lackey et al arXiv:1303.6298 tidal model.
 *
 * This is a frequency domain model that adds tidal modifications of amplitude and phasing
 * to the SEOBNRv2 model. Instead of SEOBNRv2, we use the high resolution ROM.
 *
 * @note Parameter ranges:
 *   * 6/49 <= eta <= 0.25
 *   * -0.75 <= chi_BH <= 0.75
 *   * Mtot >= 2 Msun @ 10 Hz (inherited from the ROM)
 *
 *  Aligned component spin on black hole chi_BH. The NS is assumed to be non-spinning.
 *  Symmetric mass-ratio eta = m1*m2/(m1+m2)^2.
 *  Total mass Mtot.
 *
 * @{
 */


/**
 * Compute waveform in LAL format at specified frequencies for the Lackey et al (2013)
 * tidal model based on SEOBNRv2_ROM_DoubleSpin_HI.
 *
 * XLALSimIMRLackeyTidal2013() returns the plus and cross polarizations as a complex
 * frequency series with equal spacing deltaF and contains zeros from zero frequency
 * to the starting frequency and zeros beyond the cutoff frequency in the ringdown.
 *
 * In contrast, XLALSimIMRLackeyTidal2013FrequencySequence() returns a
 * complex frequency series with entries exactly at the frequencies specified in
 * the sequence freqs (which can be unequally spaced). No zeros are added.
 *
 * If XLALSimIMRLackeyTidal2013FrequencySequence() is called with frequencies that
 * are beyond the maxium allowed geometric frequency for the ROM, zero strain is returned.
 * It is not assumed that the frequency sequence is ordered.
 *
 * This function is designed as an entry point for reduced order quadratures.
 */
int XLALSimIMRLackeyTidal2013FrequencySequence(
  struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
  struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
  const REAL8Sequence *freqs,                   /**< Frequency points at which to evaluate the waveform (Hz) */
  REAL8 phiRef,                                 /**< Phase at reference time */
  REAL8 fRef,                                   /**< Reference frequency (Hz); 0 defaults to fLow */
  REAL8 distance,                               /**< Distance of source (m) */
  REAL8 inclination,                            /**< Inclination of source (rad) */
  REAL8 mBH_SI,                                 /**< Mass of black hole (kg) */
  REAL8 mNS_SI,                                 /**< Mass of neutron star (kg) */
  REAL8 chi_BH,                                 /**< Dimensionless aligned component spin of the BH */
  REAL8 Lambda)                                 /**< Dimensionless tidal deformability (Eq 1  of Lackey et al) */
{
  if (!freqs) XLAL_ERROR(XLAL_EFAULT);

  // Call the internal core function with deltaF = 0 to indicate that freqs is non-uniformly
  // spaced and we want the strain only at these frequencies
  int retcode = LackeyTidal2013SEOBNRv2ROMCore(hptilde, hctilde,
            phiRef, fRef, distance, inclination, mBH_SI, mNS_SI, chi_BH, Lambda, freqs, 0);

  return(retcode);
}

/**
 * Compute waveform in LAL format for the Lackey et al (2013) tidal model based on
 * SEOBNRv2_ROM_DoubleSpin_HI.
 *
 * Returns the plus and cross polarizations as a complex frequency series with
 * equal spacing deltaF and contains zeros from zero frequency to the starting
 * frequency fLow and zeros beyond the cutoff frequency in the ringdown.
 */
int XLALSimIMRLackeyTidal2013(
  struct tagCOMPLEX16FrequencySeries **hptilde, /**< Output: Frequency-domain waveform h+ */
  struct tagCOMPLEX16FrequencySeries **hctilde, /**< Output: Frequency-domain waveform hx */
  REAL8 phiRef,                                 /**< Phase at reference time */
  REAL8 deltaF,                                 /**< Sampling frequency (Hz) */
  REAL8 fLow,                                   /**< Starting GW frequency (Hz) */
  REAL8 fHigh,                                  /**< End frequency; 0 defaults to Mf=0.14 */
  REAL8 fRef,                                   /**< Reference frequency (Hz); 0 defaults to fLow */
  REAL8 distance,                               /**< Distance of source (m) */
  REAL8 inclination,                            /**< Inclination of source (rad) */
  REAL8 mBH_SI,                                 /**< Mass of black hole (kg) */
  REAL8 mNS_SI,                                 /**< Mass of neutron star (kg) */
  REAL8 chi_BH,                                 /**< Dimensionless aligned component spin of the BH */
  REAL8 Lambda                                  /**< Dimensionless tidal deformability (Eq 1  of Lackey et al) */
) {
  // Use fLow, fHigh, deltaF to compute freqs sequence
  // Instead of building a full sequence we only transfer the boundaries and let
  // the internal core function do the rest (and properly take care of corner cases).
  REAL8Sequence *freqs = XLALCreateREAL8Sequence(2);
  freqs->data[0] = fLow;
  freqs->data[1] = fHigh;

  int retcode = LackeyTidal2013SEOBNRv2ROMCore(hptilde, hctilde,
            phiRef, fRef, distance, inclination, mBH_SI, mNS_SI, chi_BH, Lambda, freqs, deltaF);

  XLALDestroyREAL8Sequence(freqs);

  return(retcode);
}

/** @} */
/** @} */

