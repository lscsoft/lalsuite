/*
 * Copyright (C) 2011 P. Ajith, Nickolas Fotopoulos
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

#include <lal/LALStdlib.h>
#include <lal/LALSimIMR.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/StringInput.h>
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Units.h>

typedef struct tagBBHPhenomParams{
  REAL8 fMerger;
  REAL8 fRing;
  REAL8 fCut;
  REAL8 sigma;
  REAL8 psi0;
  REAL8 psi1;
  REAL8 psi2;
  REAL8 psi3;
  REAL8 psi4;
  REAL8 psi5;
  REAL8 psi6;
  REAL8 psi7;
  REAL8 psi8;
}
BBHPhenomParams;

/**
 *
 * private function prototypes; all internal functions use solar masses.
 *
 */

static BBHPhenomParams *ComputeIMRPhenomAParams(const REAL8 m1, const REAL8 m2);
static BBHPhenomParams *ComputeIMRPhenomBParams(const REAL8 m1, const REAL8 m2, const REAL8 chi);

static REAL8 EstimateSafeFMinForTD(const REAL8 m1, const REAL8 m2, const REAL8 f_min, const REAL8 deltaT);
static REAL8 EstimateSafeFMaxForTD(const REAL8 f_max, const REAL8 dt);
static REAL8 ComputeTau0(const REAL8 m1, const REAL8 m2, const REAL8 f_min);
static size_t EstimateIMRLength(const REAL8 m1, const REAL8 m2, const REAL8 f_min, const REAL8 deltaT);
static size_t NextPow2(const size_t n);

static REAL8 LorentzianFn(const REAL8 freq, const REAL8 fRing, const REAL8 sigma);

static int IMRPhenomAGenerateFD(COMPLEX16FrequencySeries **htilde, const REAL8 phi0, const REAL8 deltaF, const REAL8 m1, const REAL8 m2, const REAL8 f_min, const REAL8 f_max, const REAL8 distance, const BBHPhenomParams *params);
static int IMRPhenomBGenerateFD(COMPLEX16FrequencySeries **htilde, const REAL8 phi0, const REAL8 deltaF, const REAL8 m1, const REAL8 m2, const REAL8 chi, const REAL8 f_min, const REAL8 f_max, const REAL8 distance, const BBHPhenomParams *params);
static int IMRPhenomAGenerateTD(REAL8TimeSeries **h, const REAL8 phiPeak, const REAL8 deltaT, const REAL8 m1, const REAL8 m2, const REAL8 f_min, const REAL8 f_max, const REAL8 distance, const BBHPhenomParams *params);
static int IMRPhenomBGenerateTD(REAL8TimeSeries **h, const REAL8 phiPeak, const REAL8 deltaT, const REAL8 m1, const REAL8 m2, const REAL8 chi, const REAL8 f_min, const REAL8 f_max, const REAL8 distance, const BBHPhenomParams *params);
static int FDToTD(REAL8TimeSeries **signalTD, const COMPLEX16FrequencySeries *signalFD, const REAL8 totalMass, const REAL8 deltaT, const REAL8 f_min, const REAL8 f_max, const REAL8 f_min_wide, const REAL8 f_max_wide);
static size_t find_instant_freq(const REAL8TimeSeries *hp, const REAL8TimeSeries *hc, const REAL8 target, const size_t start);
static size_t find_peak_amp(const REAL8TimeSeries *hp, const REAL8TimeSeries *hc);
static int apply_phase_shift(const REAL8TimeSeries *hp, const REAL8TimeSeries *hc, const REAL8 shift);
static int apply_inclination(const REAL8TimeSeries *hp, const REAL8TimeSeries *hc, const REAL8 inclination);


/**
 *
 * main functions
 *
 */

/**
 * Driver routine to compute the non-spinning, inspiral-merger-ringdown
 * phenomenological waveform IMRPhenomA in the frequency domain.
 *
 * Reference:
 *   - Waveform: Eq.(4.13) and (4.16) of http://arxiv.org/pdf/0710.2335
 *   - Coefficients: Eq.(4.18) of http://arxiv.org/pdf/0710.2335 and
 *                   Table I of http://arxiv.org/pdf/0712.0343
 *
 * All input parameters should be SI units.
 */
int XLALSimIMRPhenomAGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< FD waveform */
    const REAL8 phi0,                  /**< phase at peak */
    const REAL8 deltaF,                /**< frequency resolution */
    const REAL8 m1_SI,                 /**< mass of companion 1 (kg) */
    const REAL8 m2_SI,                 /**< mass of companion 2 (kg) */
    const REAL8 f_min,                 /**< start frequency */
    const REAL8 f_max,                 /**< end frequency; if 0, set to fCut */
    const REAL8 distance               /**< distance of source (m) */
) {
  BBHPhenomParams *params;
  REAL8 f_max_prime;

  /* external: SI; internal: solar masses */
  const REAL8 m1 = m1_SI / LAL_MSUN_SI;
  const REAL8 m2 = m2_SI / LAL_MSUN_SI;

  /* check inputs for sanity */
  if (*htilde) XLAL_ERROR(XLAL_EFAULT);
  if (deltaF <= 0) XLAL_ERROR(XLAL_EDOM);
  if (m1 < 0) XLAL_ERROR(XLAL_EDOM);
  if (m2 < 0) XLAL_ERROR(XLAL_EDOM);
  if (f_min <= 0) XLAL_ERROR(XLAL_EDOM);
  if (f_max < 0) XLAL_ERROR(XLAL_EDOM);
  if (distance <= 0) XLAL_ERROR(XLAL_EDOM);

  /* phenomenological parameters*/
  params = ComputeIMRPhenomAParams(m1, m2);
  if (!params) XLAL_ERROR(XLAL_EFUNC);
  if (params->fCut <= f_min) {
      XLALPrintError("fCut <= f_min");
      XLAL_ERROR(XLAL_EDOM);
  }

  /* default f_max to params->fCut */
  f_max_prime = f_max ? f_max : params->fCut;
  if (f_max_prime <= f_min) {
      XLALPrintError("f_max <= f_min");
      XLAL_ERROR(XLAL_EDOM);
  }

  return IMRPhenomAGenerateFD(htilde, phi0, deltaF, m1, m2, f_min, f_max_prime, distance, params);
}

/**
 * Driver routine to compute the non-spinning, inspiral-merger-ringdown
 * phenomenological waveform IMRPhenomA in the time domain.
 *
 * Reference:
 *   - Waveform: Eq.(4.13) and (4.16) of http://arxiv.org/pdf/0710.2335
 *   - Coefficients: Eq.(4.18) of http://arxiv.org/pdf/0710.2335 and
 *                   Table I of http://arxiv.org/pdf/0712.0343
 *
 * All input parameters should be in SI units. Angles should be in radians.
 */
int XLALSimIMRPhenomAGenerateTD(
    REAL8TimeSeries **hplus,  /**< +-polarization waveform */
    REAL8TimeSeries **hcross, /**< x-polarization waveform */
    const REAL8 phiPeak,            /**< orbital phase at peak (rad) */
    const REAL8 deltaT,             /**< sampling interval (s) */
    const REAL8 m1_SI,              /**< mass of companion 1 (kg) */
    const REAL8 m2_SI,              /**< mass of companion 2 (kg) */
    const REAL8 f_min,              /**< starting GW frequency (Hz) */
    const REAL8 f_max,              /**< end GW frequency; 0 defaults to ringdown cutoff freq */
    const REAL8 distance,           /**< distance of source (m) */
    const REAL8 inclination         /**< inclination of source (rad) */
) {
  BBHPhenomParams *params;
  size_t cut_ind, peak_ind;
  REAL8 peak_phase;  /* measured, not intended */
  REAL8 f_max_prime;

  /* external: SI; internal: solar masses */
  const REAL8 m1 = m1_SI / LAL_MSUN_SI;
  const REAL8 m2 = m2_SI / LAL_MSUN_SI;

  /* check inputs for sanity */
  if (*hplus) XLAL_ERROR(XLAL_EFAULT);
  if (*hcross) XLAL_ERROR(XLAL_EFAULT);
  if (deltaT <= 0) XLAL_ERROR(XLAL_EDOM);
  if (m1 < 0) XLAL_ERROR(XLAL_EDOM);
  if (m2 < 0) XLAL_ERROR(XLAL_EDOM);
  if (f_min <= 0) XLAL_ERROR(XLAL_EDOM);
  if (f_max < 0) XLAL_ERROR(XLAL_EDOM);
  if (distance <= 0) XLAL_ERROR(XLAL_EDOM);

  /* phenomenological parameters*/
  params = ComputeIMRPhenomAParams(m1, m2);
  if (!params) XLAL_ERROR(XLAL_EFUNC);
  if (params->fCut <= f_min) {
      XLALPrintError("fCut <= f_min");
      XLAL_ERROR(XLAL_EDOM);
  }

  /* default f_max to params->fCut */
  f_max_prime = f_max ? f_max : params->fCut;
  if (f_max_prime <= f_min) {
      XLALPrintError("f_max <= f_min");
      XLAL_ERROR(XLAL_EDOM);
  }

  /* generate hplus */
  IMRPhenomAGenerateTD(hplus, 0, deltaT, m1, m2, f_min, f_max_prime, distance, params);
  if (!(*hplus)) {
      XLALFree(params);
      XLAL_ERROR(XLAL_EFUNC);
  }

  /* generate hcross, which is hplus phase-shifted by pi/2 */
  IMRPhenomAGenerateTD(hcross, LAL_PI_2, deltaT, m1, m2, f_min, f_max_prime, distance, params);
  XLALFree(params);
  if (!(*hcross)) {
      XLALDestroyREAL8TimeSeries(*hplus);
      *hplus = NULL;
      XLAL_ERROR(XLAL_EFUNC);
  }

  /* clip the parts below f_min */
  cut_ind = find_instant_freq(*hplus, *hcross, f_min, (*hplus)->data->length - EstimateIMRLength(m1, m2, f_min, deltaT) + EstimateIMRLength(m1, m2, f_max_prime, deltaT));
  *hplus = XLALResizeREAL8TimeSeries(*hplus, cut_ind, (*hplus)->data->length - cut_ind);
  *hcross = XLALResizeREAL8TimeSeries(*hcross, cut_ind, (*hcross)->data->length - cut_ind);
  if (!(*hplus) || !(*hcross))
    XLAL_ERROR(XLAL_EFUNC);

  /* set phase and time at peak */
  peak_ind = find_peak_amp(*hplus, *hcross);
  peak_phase = atan2((*hcross)->data->data[peak_ind], (*hplus)->data->data[peak_ind]);
  // NB: factor of 2 b/c phiPeak is *orbital* phase, and we're shifting GW phase
  apply_phase_shift(*hplus, *hcross, 2.*phiPeak - peak_phase);
  XLALGPSSetREAL8(&((*hplus)->epoch), -(peak_ind * deltaT));
  XLALGPSSetREAL8(&((*hcross)->epoch), -(peak_ind * deltaT));

  /* apply inclination */
  return apply_inclination(*hplus, *hcross, inclination);
}

/**
 * Compute the dimensionless, spin-aligned parameter chi as used in the
 * IMRPhenomB waveform. This is different from chi in SpinTaylorRedSpin!
 * Reference: http://arxiv.org/pdf/0909.2867, paragraph 3.
 */
double XLALSimIMRPhenomBComputeChi(
    const REAL8 m1,                          /**< mass of companion 1 */
    const REAL8 m2,                          /**< mass of companion 2 */
    const REAL8 s1z,                         /**< spin of companion 1 */
    const REAL8 s2z                          /**< spin of companion 2 */
) {
    return (m1 * s1z + m2 * s2z) / (m1 + m2);
}

/**
 * Driver routine to compute the spin-aligned, inspiral-merger-ringdown
 * phenomenological waveform IMRPhenomB in the time domain.
 *
 * Reference: http://arxiv.org/pdf/0909.2867
 *   - Waveform: Eq.(1)
 *   - Coefficients: Eq.(2) and Table I
 *
 * All input parameters should be in SI units. Angles should be in radians.
 */
int XLALSimIMRPhenomBGenerateTD(
    REAL8TimeSeries **hplus,  /**< +-polarization waveform */
    REAL8TimeSeries **hcross, /**< x-polarization waveform */
    const REAL8 phiPeak,      /**< orbital phase at peak (rad) */
    const REAL8 deltaT,       /**< sampling interval (s) */
    const REAL8 m1_SI,        /**< mass of companion 1 (kg) */
    const REAL8 m2_SI,        /**< mass of companion 2 (kg) */
    const REAL8 chi,          /**< mass-weighted aligned-spin parameter */
    const REAL8 f_min,        /**< starting GW frequency (Hz) */
    const REAL8 f_max,        /**< end GW frequency; 0 defaults to ringdown cutoff freq */
    const REAL8 distance,     /**< distance of source (m) */
    const REAL8 inclination   /**< inclination of source (rad) */
) {
  BBHPhenomParams *params;
  size_t cut_ind, peak_ind;
  REAL8 peak_phase;  /* measured, not intended */
  REAL8 f_max_prime;

  /* external: SI; internal: solar masses */
  const REAL8 m1 = m1_SI / LAL_MSUN_SI;
  const REAL8 m2 = m2_SI / LAL_MSUN_SI;

  /* check inputs for sanity */
  if (*hplus) XLAL_ERROR(XLAL_EFAULT);
  if (*hcross) XLAL_ERROR(XLAL_EFAULT);
  if (deltaT <= 0) XLAL_ERROR(XLAL_EDOM);
  if (m1 < 0) XLAL_ERROR(XLAL_EDOM);
  if (m2 < 0) XLAL_ERROR(XLAL_EDOM);
  if (fabs(chi) > 1) XLAL_ERROR(XLAL_EDOM);
  if (f_min <= 0) XLAL_ERROR(XLAL_EDOM);
  if (f_max < 0) XLAL_ERROR(XLAL_EDOM);
  if (distance <= 0) XLAL_ERROR(XLAL_EDOM);

  /* phenomenological parameters*/
  params = ComputeIMRPhenomBParams(m1, m2, chi);
  if (!params) XLAL_ERROR(XLAL_EFUNC);
  if (params->fCut <= f_min) {
      XLALPrintError("fCut <= f_min");
      XLAL_ERROR(XLAL_EDOM);
  }

  /* default f_max to params->fCut */
  f_max_prime = f_max ? f_max : params->fCut;
  if (f_max_prime <= f_min) {
      XLALPrintError("f_max <= f_min");
      XLAL_ERROR(XLAL_EDOM);
  }

  /* generate plus */
  IMRPhenomBGenerateTD(hplus, 0., deltaT, m1, m2, chi, f_min, f_max_prime, distance, params);
  if (!(*hplus)) {
      XLALFree(params);
      XLAL_ERROR(XLAL_EFUNC);
  }

  /* generate cross, phase-shifted by pi/2 */
  IMRPhenomBGenerateTD(hcross, LAL_PI_2, deltaT, m1, m2, chi, f_min, f_max_prime, distance, params);
  XLALFree(params);
  if (!(*hcross)) {
      XLALDestroyREAL8TimeSeries(*hplus);
      *hplus = NULL;
      XLAL_ERROR(XLAL_EFUNC);
  }

  /* clip the parts below f_min */
  cut_ind = find_instant_freq(*hplus, *hcross, f_min, (*hplus)->data->length - EstimateIMRLength(m1, m2, f_min, deltaT) + EstimateIMRLength(m1, m2, f_max_prime, deltaT));
  *hplus = XLALResizeREAL8TimeSeries(*hplus, cut_ind, (*hplus)->data->length - cut_ind);
  *hcross = XLALResizeREAL8TimeSeries(*hcross, cut_ind, (*hcross)->data->length - cut_ind);
  if (!(*hplus) || !(*hcross))
    XLAL_ERROR(XLAL_EFUNC);

  /* set phase and time at peak */
  peak_ind = find_peak_amp(*hplus, *hcross);
  peak_phase = atan2((*hcross)->data->data[peak_ind], (*hplus)->data->data[peak_ind]);
  // NB: factor of 2 b/c phiPeak is *orbital* phase, and we're shifting GW phase
  apply_phase_shift(*hplus, *hcross, 2.*phiPeak - peak_phase);
  XLALGPSSetREAL8(&((*hplus)->epoch), -(peak_ind * deltaT));
  XLALGPSSetREAL8(&((*hcross)->epoch), -(peak_ind * deltaT));

  /* apply inclination */
  return apply_inclination(*hplus, *hcross, inclination);
}


/**
 * Driver routine to compute the spin-aligned, inspiral-merger-ringdown
 * phenomenological waveform IMRPhenomB in the frequency domain.
 *
 * Reference: http://arxiv.org/pdf/0909.2867
 *   - Waveform: Eq.(1)
 *   - Coefficients: Eq.(2) and Table I
 *
 * All input parameters should be in SI units. Angles should be in radians.
 */
int XLALSimIMRPhenomBGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< FD waveform */
    const REAL8 phi0,                  /**< initial phase */
    const REAL8 deltaF,                /**< sampling interval */
    const REAL8 m1_SI,                 /**< mass of companion 1 (kg) */
    const REAL8 m2_SI,                 /**< mass of companion 2 (kg) */
    const REAL8 chi,                   /**< mass-weighted aligned-spin parameter */
    const REAL8 f_min,                 /**< start frequency */
    const REAL8 f_max,                 /**< end frequency; 0 defaults to ringdown cutoff freq */
    const REAL8 distance               /**< distance of source (m) */
) {
  BBHPhenomParams *params;
  int status;
  REAL8 f_max_prime;

  /* external: SI; internal: solar masses */
  const REAL8 m1 = m1_SI / LAL_MSUN_SI;
  const REAL8 m2 = m2_SI / LAL_MSUN_SI;

  /* check inputs for sanity */
  if (*htilde) XLAL_ERROR(XLAL_EFAULT);
  if (deltaF <= 0) XLAL_ERROR(XLAL_EDOM);
  if (m1 < 0) XLAL_ERROR(XLAL_EDOM);
  if (m2 < 0) XLAL_ERROR(XLAL_EDOM);
  if (fabs(chi) > 1) XLAL_ERROR(XLAL_EDOM);
  if (f_min <= 0) XLAL_ERROR(XLAL_EDOM);
  if (f_max < 0) XLAL_ERROR(XLAL_EDOM);
  if (distance <= 0) XLAL_ERROR(XLAL_EDOM);

  /* phenomenological parameters*/
  params = ComputeIMRPhenomBParams(m1, m2, chi);
  if (!params) XLAL_ERROR(XLAL_EFUNC);
  if (params->fCut <= f_min) {
      XLALPrintError("fCut <= f_min");
      XLAL_ERROR(XLAL_EDOM);
  }

  /* default f_max to params->fCut */
  f_max_prime = f_max ? f_max : params->fCut;
  if (f_max_prime <= f_min) {
      XLALPrintError("f_max <= f_min");
      XLAL_ERROR(XLAL_EDOM);
  }

  status = IMRPhenomBGenerateFD(htilde, phi0, deltaF, m1, m2, chi, f_min, f_max_prime, distance, params);
  LALFree(params);
  return status;
}

/*********************************************************************/
/* Compute phenomenological parameters for non-spinning binaries     */
/* Ref. Eq.(4.18) of http://arxiv.org/pdf/0710.2335 and              */
/* Table I of http://arxiv.org/pdf/0712.0343                         */
/*                                                                   */
/* Takes solar masses.                                               */
/*********************************************************************/
static BBHPhenomParams *ComputeIMRPhenomAParams(const REAL8 m1, const REAL8 m2) {
  /* calculate the total mass and symmetric mass ratio */
  const REAL8 totalMass = m1 + m2;
  const REAL8 eta = m1 * m2 / (totalMass * totalMass);
  const REAL8 piM = totalMass * LAL_PI * LAL_MTSUN_SI;
  const REAL8 etap2 = eta*eta;

  const REAL8 fMerg_a = 6.6389e-01;
  const REAL8 fMerg_b = -1.0321e-01;
  const REAL8 fMerg_c = 1.0979e-01;

  const REAL8 fRing_a = 1.3278e+00;
  const REAL8 fRing_b = -2.0642e-01;
  const REAL8 fRing_c = 2.1957e-01;

  const REAL8 sigma_a = 1.1383e+00;
  const REAL8 sigma_b = -1.7700e-01;
  const REAL8 sigma_c = 4.6834e-02;

  const REAL8 fCut_a = 1.7086e+00;
  const REAL8 fCut_b = -2.6592e-01;
  const REAL8 fCut_c = 2.8236e-01;

  const REAL8 psi0_a = -1.5829e-01;
  const REAL8 psi0_b = 8.7016e-02;
  const REAL8 psi0_c = -3.3382e-02;

  const REAL8 psi2_a = 3.2967e+01;
  const REAL8 psi2_b = -1.9000e+01;
  const REAL8 psi2_c = 2.1345e+00;

  const REAL8 psi3_a = -3.0849e+02;
  const REAL8 psi3_b = 1.8211e+02;
  const REAL8 psi3_c = -2.1727e+01;

  const REAL8 psi4_a = 1.1525e+03;
  const REAL8 psi4_b = -7.1477e+02;
  const REAL8 psi4_c = 9.9692e+01;

  const REAL8 psi6_a = 1.2057e+03;
  const REAL8 psi6_b = -8.4233e+02;
  const REAL8 psi6_c = 1.8046e+02;

  const REAL8 psi7_a = 0.;
  const REAL8 psi7_b = 0.;
  const REAL8 psi7_c = 0.;

  BBHPhenomParams *phenParams = (BBHPhenomParams *) XLALMalloc(sizeof(BBHPhenomParams));
  if (!phenParams) XLAL_ERROR_NULL(XLAL_EFUNC);
  memset(phenParams, 0, sizeof(BBHPhenomParams));

  /* Evaluate the polynomials. See Eq. (4.18) of P. Ajith et al
   * arXiv:0710.2335 [gr-qc] */
  phenParams->fCut = (fCut_a*etap2  + fCut_b*eta  + fCut_c)/piM;
  phenParams->fMerger = (fMerg_a*etap2  + fMerg_b*eta  + fMerg_c)/piM;
  phenParams->fRing = (fRing_a*etap2 + fRing_b*eta + fRing_c)/piM;
  phenParams->sigma = (sigma_a*etap2 + sigma_b*eta + sigma_c)/piM;

  phenParams->psi0 = (psi0_a*etap2 + psi0_b*eta + psi0_c)/(eta*pow(piM, 5./3.));
  phenParams->psi1 = 0.;
  phenParams->psi2 = (psi2_a*etap2 + psi2_b*eta + psi2_c)/(eta*pow(piM, 3./3.));
  phenParams->psi3 = (psi3_a*etap2 + psi3_b*eta + psi3_c)/(eta*pow(piM, 2./3.));
  phenParams->psi4 = (psi4_a*etap2 + psi4_b*eta + psi4_c)/(eta*cbrt(piM));
  phenParams->psi5 = 0.;
  phenParams->psi6 = (psi6_a*etap2 + psi6_b*eta + psi6_c)/(eta/cbrt(piM));
  phenParams->psi7 = (psi7_a*etap2 + psi7_b*eta + psi7_c)/(eta*pow(piM, -2./3.));

  return phenParams;
}

/*********************************************************************/
/* Compute phenomenological parameters for aligned-spin binaries     */
/* Ref. Eq.(2) and Table I of http://arxiv.org/pdf/0909.2867         */
/*                                                                   */
/* Takes solar masses. Populates and returns a new BBHPhenomParams   */
/* structure.                                                        */
/*********************************************************************/
static BBHPhenomParams *ComputeIMRPhenomBParams(const REAL8 m1, const REAL8 m2, const REAL8 chi) {
  /* calculate the total mass and symmetric mass ratio */
  const REAL8 totalMass = m1 + m2;
  const REAL8 eta = m1 * m2 / (totalMass * totalMass);
  const REAL8 piM = totalMass * LAL_PI * LAL_MTSUN_SI;

  /* spinning phenomenological waveforms */
  const REAL8 etap2 = eta*eta;
  const REAL8 chip2 = chi*chi;
  const REAL8 etap3 = etap2*eta;
  const REAL8 etap2chi = etap2*chi;
  const REAL8 etachip2 = eta*chip2;
  const REAL8 etachi = eta*chi;

  BBHPhenomParams *phenParams = (BBHPhenomParams *) XLALMalloc(sizeof(BBHPhenomParams));
  if (!phenParams) XLAL_ERROR_NULL(XLAL_EFUNC);
  memset(phenParams, 0, sizeof(BBHPhenomParams));

  phenParams->psi0 = 3./(128.*eta);

  phenParams->psi2 = 3715./756. +
  -9.2091e+02*eta + 4.9213e+02*etachi + 1.3503e+02*etachip2 +
  6.7419e+03*etap2 + -1.0534e+03*etap2chi +
  -1.3397e+04*etap3 ;

  phenParams->psi3 = -16.*LAL_PI + 113.*chi/3. +
  1.7022e+04*eta + -9.5659e+03*etachi + -2.1821e+03*etachip2 +
  -1.2137e+05*etap2 + 2.0752e+04*etap2chi +
  2.3859e+05*etap3 ;

  phenParams->psi4 = 15293365./508032. - 405.*chip2/8. +
  -1.2544e+05*eta + 7.5066e+04*etachi + 1.3382e+04*etachip2 +
  8.7354e+05*etap2 + -1.6573e+05*etap2chi +
  -1.6936e+06*etap3 ;

  phenParams->psi6 = -8.8977e+05*eta + 6.3102e+05*etachi + 5.0676e+04*etachip2 +
  5.9808e+06*etap2 + -1.4148e+06*etap2chi +
  -1.1280e+07*etap3 ;

  phenParams->psi7 = 8.6960e+05*eta + -6.7098e+05*etachi + -3.0082e+04*etachip2 +
  -5.8379e+06*etap2 + 1.5145e+06*etap2chi +
  1.0891e+07*etap3 ;

  phenParams->psi8 = -3.6600e+05*eta + 3.0670e+05*etachi + 6.3176e+02*etachip2 +
  2.4265e+06*etap2 + -7.2180e+05*etap2chi +
  -4.5524e+06*etap3;

  phenParams->fMerger =  1. - 4.4547*pow(1.-chi,0.217) + 3.521*pow(1.-chi,0.26) +
  6.4365e-01*eta + 8.2696e-01*etachi + -2.7063e-01*etachip2 +
  -5.8218e-02*etap2 + -3.9346e+00*etap2chi +
  -7.0916e+00*etap3 ;

  phenParams->fRing = (1. - 0.63*pow(1.-chi,0.3))/2. +
  1.4690e-01*eta + -1.2281e-01*etachi + -2.6091e-02*etachip2 +
  -2.4900e-02*etap2 + 1.7013e-01*etap2chi +
  2.3252e+00*etap3 ;

  phenParams->sigma = (1. - 0.63*pow(1.-chi,0.3))*pow(1.-chi,0.45)/4. +
  -4.0979e-01*eta + -3.5226e-02*etachi + 1.0082e-01*etachip2 +
  1.8286e+00*etap2 + -2.0169e-02*etap2chi +
  -2.8698e+00*etap3 ;

  phenParams->fCut = 3.2361e-01 + 4.8935e-02*chi + 1.3463e-02*chip2 +
  -1.3313e-01*eta + -8.1719e-02*etachi + 1.4512e-01*etachip2 +
  -2.7140e-01*etap2 + 1.2788e-01*etap2chi +
  4.9220e+00*etap3 ;

  phenParams->fCut   /= piM;
  phenParams->fMerger/= piM;
  phenParams->fRing  /= piM;
  phenParams->sigma  /= piM;

  phenParams->psi1    = 0.;
  phenParams->psi5    = 0.;

  return phenParams;
}

/**
 * Return tau0, the Newtonian chirp length estimate.
 */
static REAL8 ComputeTau0(const REAL8 m1, const REAL8 m2, const REAL8 f_min) {
  const REAL8 totalMass = m1 + m2;
  const REAL8 eta = m1 * m2 / (totalMass * totalMass);
  return 5. * totalMass * LAL_MTSUN_SI / (256. * eta * pow(LAL_PI * totalMass * LAL_MTSUN_SI * f_min, 8./3.));
}

/**
 * Estimate the length of a TD vector that can hold the waveform as the Newtonian
 * chirp time tau0 plus 1000 M.
 */
static size_t EstimateIMRLength(const REAL8 m1, const REAL8 m2, const REAL8 f_min, const REAL8 deltaT) {
  return (size_t) floor((ComputeTau0(m1, m2, f_min) + 1000 * (m1 + m2) * LAL_MTSUN_SI) / deltaT);
}

static size_t NextPow2(const size_t n) {
  return 1 << (size_t) ceil(log2(n));
}

/**
 * Find a lower value for f_min (using the definition of Newtonian chirp
 * time) such that the waveform has a minimum length of tau0. This is
 * necessary to avoid FFT artifacts.
 */
static REAL8 EstimateSafeFMinForTD(const REAL8 m1, const REAL8 m2, const REAL8 f_min, const REAL8 deltaT) {
  const REAL8 totalMass = m1 + m2;
  const REAL8 eta = m1 * m2 / (totalMass * totalMass);
  REAL8 tau0 = deltaT * NextPow2(1.025 * EstimateIMRLength(m1, m2, f_min, deltaT));
  REAL8 temp_f_min = pow((tau0 * 256. * eta * pow(totalMass * LAL_MTSUN_SI, 5./3.) / 5.), -3./8.) / LAL_PI;
  if (temp_f_min > f_min) temp_f_min = f_min;
  if (temp_f_min < 0.5) temp_f_min = 0.5;
  return temp_f_min;
}

/**
 * Find a higher value of f_max so that we can safely apply a window later.
 */
static REAL8 EstimateSafeFMaxForTD(const REAL8 f_max, const REAL8 deltaT) {
  REAL8 temp_f_max = 1.025 * f_max;

  /* make sure that these frequencies are not too out of range */
  if (temp_f_max > 2. / deltaT - 100.) temp_f_max = 2. / deltaT - 100.;
  return temp_f_max;
}

static REAL8 LorentzianFn(const REAL8 freq, const REAL8 fRing, const REAL8 sigma) {
  return sigma / (LAL_TWOPI * ((freq - fRing)*(freq - fRing)
    + sigma*sigma / 4.0));
}


/**
 * Private function to generate IMRPhenomA frequency-domain waveforms given coefficients
 */
static int IMRPhenomAGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< FD waveform */
    const REAL8 phi0,                  /**< initial phase */
    const REAL8 deltaF,                /**< frequency resolution */
    const REAL8 m1,                    /**< mass of companion 1 [solar masses] */
    const REAL8 m2,                    /**< mass of companion 2 [solar masses] */
    const REAL8 f_min,                 /**< start frequency */
    const REAL8 f_max,                 /**< end frequency */
    const REAL8 distance,              /**< distance of source */
    const BBHPhenomParams *params      /**< from ComputeIMRPhenomAParams */
) {
  const LIGOTimeGPS ligotimegps_zero = {0, 0};
  COMPLEX16 *data;
  size_t i;

  const REAL8 fMerg = params->fMerger;
  const REAL8 fRing = params->fRing;
  const REAL8 sigma = params->sigma;
  const REAL8 totalMass = m1 + m2;
  const REAL8 eta = m1 * m2 / (totalMass * totalMass);

  /* compute the amplitude pre-factor */
  const REAL8 amp0 = pow(LAL_MTSUN_SI*totalMass, 5./6.) * pow(fMerg, -7./6.)
    / pow(LAL_PI, 2./3.) * sqrt(5. * eta / 24.) / (distance / LAL_C_SI);

  /* allocate htilde */
  size_t n = NextPow2(f_max / deltaF) + 1;
  *htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &ligotimegps_zero, 0.0, deltaF, &lalStrainUnit, n);
  memset((*htilde)->data->data, 0, n * sizeof(COMPLEX16));
  XLALUnitMultiply(&((*htilde)->sampleUnits), &((*htilde)->sampleUnits), &lalSecondUnit);
  if (!(*htilde)) XLAL_ERROR(XLAL_EFUNC);

  /* now generate the waveform at all frequency bins except DC and Nyquist */
  data = (*htilde)->data->data;
  for (i=1; i < n - 1; i++) {
    REAL8 ampEff, psiEff;
    /* Fourier frequency corresponding to this bin */
    REAL8 f = i * deltaF;
    REAL8 cbrt_f = cbrt(f);
    REAL8 fNorm = f / fMerg;

    /* compute the amplitude */
    if ((f < f_min) || (f > f_max)) continue;
    else if (f <= fMerg) ampEff = amp0 * pow(fNorm, -7./6.);
    else if ((f > fMerg) & (f <= fRing)) ampEff = amp0 * pow(fNorm, -2./3.);
    else if (f > fRing)
      ampEff = amp0 * LAL_PI_2 * pow(fRing / fMerg, -2./3.) * sigma
        * LorentzianFn(f, fRing, sigma);
    else {
      XLALDestroyCOMPLEX16FrequencySeries(*htilde);
      *htilde = NULL;
      XLAL_ERROR(XLAL_EDOM);
    }

    /* now compute the phase */
    psiEff = phi0
      + params->psi0 / (f * f) * cbrt_f  /* f^-5/3 */
      + params->psi1 / (f * cbrt_f)      /* f^-4/3 */
      + params->psi2 / f                 /* f^-3/3 */
      + params->psi3 / f * cbrt_f        /* f^-2/3 */
      + params->psi4 / cbrt_f            /* f^-1/3 */
      + params->psi5                     /* f^0/3 */
      + params->psi6 * cbrt_f            /* f^1/3 */
      + params->psi7 * f / cbrt_f;       /* f^2/3 */

    /* generate the waveform */
    data[i] = ampEff * cos(psiEff);
    data[i] += I * ampEff * sin(psiEff);
  }

  return XLAL_SUCCESS;
}

/**
 * Private function to generate IMRPhenomB frequency-domain waveforms given coefficients
 */
static int IMRPhenomBGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< FD waveform */
    const REAL8 phi0,                  /**< phase at peak */
    const REAL8 deltaF,                /**< frequency resolution */
    const REAL8 m1,                    /**< mass of companion 1 [solar masses] */
    const REAL8 m2,                    /**< mass of companion 2 [solar masses] */
    const REAL8 chi,                   /**< mass-weighted aligned-spin parameter */
    const REAL8 f_min,                 /**< start frequency */
    const REAL8 f_max,                 /**< end frequency */
    const REAL8 distance,              /**< distance of source */
    const BBHPhenomParams *params      /**< from ComputeIMRPhenomBParams */
) {
  static LIGOTimeGPS ligotimegps_zero = {0, 0};
  size_t i;

  const REAL8 fMerg = params->fMerger;
  const REAL8 fRing = params->fRing;
  const REAL8 sigma = params->sigma;
  const REAL8 totalMass = m1 + m2;
  const REAL8 eta = m1 * m2 / (totalMass * totalMass);
  const REAL8 piM = LAL_PI * totalMass * LAL_MTSUN_SI;

  /* compute the amplitude pre-factor */
  REAL8 amp0 = pow(LAL_MTSUN_SI*totalMass, 5./6.) * pow(fMerg, -7./6.)
    / pow(LAL_PI, 2./3.) * sqrt(5. * eta / 24.) / (distance / LAL_C_SI);

  /***********************************************************************/
  /* these are the parameters required for the "new" phenomenological IMR
   * waveforms*/
  /***********************************************************************/

  /* PN corrections to the frequency domain amplitude of the (2,2) mode */
  const REAL8 alpha2   = -323./224. + 451.*eta/168.;
  const REAL8 alpha3   = (27./8. - 11.*eta/6.)*chi;

  /* leading order power law of the merger amplitude */
  static REAL8 mergPower = -2./3.;

  /* spin-dependent corrections to the merger amplitude */
  const REAL8 epsilon_1 =  1.4547*chi - 1.8897;
  const REAL8 epsilon_2 = -1.8153*chi + 1.6557;

  /* normalisation constant of the inspiral amplitude */
  const REAL8 vMerg = cbrt(piM * fMerg);
  const REAL8 vRing = cbrt(piM * fRing);

  REAL8 w1 = (1. + alpha2 * vMerg * vMerg + alpha3 * piM * fMerg) / (1. + epsilon_1 * vMerg + epsilon_2 * vMerg * vMerg);
  REAL8 w2 = w1 * (LAL_PI * sigma / 2.) * pow(fRing / fMerg, mergPower)
          * (1. + epsilon_1 * vRing + epsilon_2 * vRing * vRing);

  /* allocate htilde */
  size_t n = NextPow2(f_max / deltaF) + 1;
  *htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &ligotimegps_zero, 0.0, deltaF, &lalStrainUnit, n);
  memset((*htilde)->data->data, 0, n * sizeof(COMPLEX16));
  XLALUnitMultiply(&((*htilde)->sampleUnits), &((*htilde)->sampleUnits), &lalSecondUnit);
  if (!(*htilde)) XLAL_ERROR(XLAL_EFUNC);

  /* now generate the waveform */
  size_t ind_max = (size_t) (f_max / deltaF);
  for (i = (size_t) (f_min / deltaF); i < ind_max; i++) {
    REAL8 ampEff, psiEff;
    REAL8 v, v2, v3, v4, v5, v6, v7, v8;

    /* Fourier frequency corresponding to this bin */
    REAL8 f = i * deltaF;

    /* PN expansion parameter */
    v = cbrt(piM * f);
    v2 = v*v; v3 = v2*v; v4 = v2*v2; v5 = v4*v; v6 = v3*v3; v7 = v6*v, v8 = v7*v;

    /* compute the amplitude */
    if (f <= fMerg)
      ampEff = pow(f / fMerg, -7./6.)*(1. + alpha2 * v2 + alpha3 * v3);
    else if (f > fRing)
      ampEff = w2 * LorentzianFn(f, fRing, sigma);
    else /* fMerg < f <= fRing */
      ampEff = w1 * pow(f / fMerg, mergPower) * (1. + epsilon_1 * v + epsilon_2 * v2);

    /* now compute the phase */
    psiEff = -phi0  /* phi is flipped relative to IMRPhenomA */
      + 3./(128.*eta*v5)*(1 + params->psi2*v2
      + params->psi3*v3 + params->psi4*v4
      + params->psi5*v5 + params->psi6*v6
      + params->psi7*v7 + params->psi8*v8);

    /* generate the waveform */
    ((*htilde)->data->data)[i] = amp0 * ampEff * cos(psiEff);
    ((*htilde)->data->data)[i] += -I * amp0 * ampEff * sin(psiEff);
  }

  return XLAL_SUCCESS;
}

/**
 * Private function to generate time-domain waveforms given coefficients
 */
static int IMRPhenomAGenerateTD(REAL8TimeSeries **h, const REAL8 phi0, const REAL8 deltaT, const REAL8 m1, const REAL8 m2, const REAL8 f_min, const REAL8 f_max, const REAL8 distance, const BBHPhenomParams *params) {
  COMPLEX16FrequencySeries *htilde=NULL;
  /* We will generate the waveform from a frequency which is lower than the
   * f_min chosen. Also the cutoff frequency may be higher than the f_max. We
   * will later apply a window function, and truncate the time-domain waveform
   * below an instantaneous frequency f_min. */
  REAL8 f_min_wide = EstimateSafeFMinForTD(m1, m2, f_min, deltaT);
  const REAL8 f_max_wide = 0.5 / deltaT;
  REAL8 deltaF;

  if (EstimateSafeFMaxForTD(f_max, deltaT) > f_max_wide)
    XLALPrintWarning("Warning: sampling rate too low to capture chosen f_max\n");
  deltaF = 1. / (deltaT * NextPow2(EstimateIMRLength(m1, m2, f_min_wide, deltaT)));

  /* generate in frequency domain */
  if (IMRPhenomAGenerateFD(&htilde, phi0, deltaF, m1, m2, f_min_wide, f_max_wide, distance, params)) XLAL_ERROR(XLAL_EFUNC);

  /* convert to time domain */
  FDToTD(h, htilde, m1 + m2, deltaT, f_min, f_max, f_min_wide, f_max_wide);
  XLALDestroyCOMPLEX16FrequencySeries(htilde);
  if (!*h) XLAL_ERROR(XLAL_EFUNC);

  return XLAL_SUCCESS;
}

/**
 * Private function to generate time-domain waveforms given coefficients
 */
static int IMRPhenomBGenerateTD(REAL8TimeSeries **h, const REAL8 phi0, const REAL8 deltaT, const REAL8 m1, const REAL8 m2, const REAL8 chi, const REAL8 f_min, const REAL8 f_max, const REAL8 distance, const BBHPhenomParams *params) {
  REAL8 deltaF;
  COMPLEX16FrequencySeries *htilde=NULL;
  /* We will generate the waveform from a frequency which is lower than the
   * f_min chosen. Also the cutoff frequency is higher than the f_max. We
   * will later apply a window function, and truncate the time-domain waveform
   * below an instantaneous frequency f_min. */
  REAL8 f_min_wide = EstimateSafeFMinForTD(m1, m2, f_min, deltaT);
  const REAL8 f_max_wide = 0.5 / deltaT;
  if (EstimateSafeFMaxForTD(f_max, deltaT) > f_max_wide)
    XLALPrintWarning("Warning: sampling rate (%" LAL_REAL8_FORMAT " Hz) too low for expected spectral content (%" LAL_REAL8_FORMAT " Hz) \n", deltaT, EstimateSafeFMaxForTD(f_max, deltaT));
  deltaF = 1. / (deltaT * NextPow2(EstimateIMRLength(m1, m2, f_min_wide, deltaT)));

  /* generate in frequency domain */
  if (IMRPhenomBGenerateFD(&htilde, phi0, deltaF, m1, m2, chi, f_min_wide, f_max_wide, distance, params)) XLAL_ERROR(XLAL_EFUNC);

  /* convert to time domain */
  FDToTD(h, htilde, m1 + m2, deltaT, f_min, f_max, f_min_wide, f_max_wide);
  XLALDestroyCOMPLEX16FrequencySeries(htilde);
  if (!*h) XLAL_ERROR(XLAL_EFUNC);

  return XLAL_SUCCESS;
}

/**
 * Window and IFFT a FD waveform to TD, then window in TD.
 * Requires that the FD waveform be generated outside of f_min and f_max.
 * FD waveform is modified.
 */
static int FDToTD(REAL8TimeSeries **signalTD, const COMPLEX16FrequencySeries *signalFD, const REAL8 totalMass, const REAL8 deltaT, const REAL8 f_min, const REAL8 f_max, const REAL8 f_min_wide, const REAL8 f_max_wide) {
  const LIGOTimeGPS gpstime_zero = {0, 0};
  const size_t nf = signalFD->data->length;
  const size_t nt = 2 * (nf - 1);
  const REAL8 windowLength = 20. * totalMass * LAL_MTSUN_SI / deltaT;
  const REAL8 winFLo = (f_min + f_min_wide) / 2.;
  REAL8 winFHi = (f_max + f_max_wide) / 2.;
  COMPLEX16 *FDdata = signalFD->data->data;
  REAL8FFTPlan *revPlan;
  REAL8 *TDdata;
  size_t k;

  /* check inputs */
  if (f_min_wide >= f_min) XLAL_ERROR(XLAL_EDOM);

  /* apply the softening window function */
  if (winFHi > 0.5 / deltaT) winFHi = 0.5 / deltaT;
  for (k = nf;k--;) {
    const REAL8 f = k / (deltaT * nt);
    REAL8 softWin = (1. + tanh(f - winFLo))
                  * (1. - tanh(f - winFHi)) / 4.;
    FDdata[k] *= softWin;
  }

  /* allocate output */
  *signalTD = XLALCreateREAL8TimeSeries("h", &gpstime_zero, 0.0, deltaT, &lalStrainUnit, nt);

  /* Inverse Fourier transform */
  revPlan = XLALCreateReverseREAL8FFTPlan(nt, 1);
  if (!revPlan) {
    XLALDestroyREAL8TimeSeries(*signalTD);
    *signalTD = NULL;
    XLAL_ERROR(XLAL_EFUNC);
  }
  XLALREAL8FreqTimeFFT(*signalTD, signalFD, revPlan);
  XLALDestroyREAL8FFTPlan(revPlan);
  if (!(*signalTD)) XLAL_ERROR(XLAL_EFUNC);

  /* apply a linearly decreasing window at the end
   * of the waveform in order to avoid edge effects. */
  if (windowLength > (*signalTD)->data->length) XLAL_ERROR(XLAL_ERANGE);
  TDdata = (*signalTD)->data->data;
  for (k = windowLength; k--;)
    TDdata[nt-k-1] *= k / windowLength;

  return XLAL_SUCCESS;
}

/* return the index before the instantaneous frequency rises past target */
static size_t find_instant_freq(const REAL8TimeSeries *hp, const REAL8TimeSeries *hc, const REAL8 target, const size_t start) {
  size_t k = start + 1;
  const size_t n = hp->data->length - 1;

  /* Use second order differencing to find the instantaneous frequency as
   * h = A e^(2 pi i f t) ==> f = d/dt(h) / (2*pi*h) */
  for (; k < n; k++) {
    const REAL8 hpDot = (hp->data->data[k+1] - hp->data->data[k-1]) / (2 * hp->deltaT);
    const REAL8 hcDot = (hc->data->data[k+1] - hc->data->data[k-1]) / (2 * hc->deltaT);
    REAL8 f = -hcDot * hp->data->data[k] + hpDot * hc->data->data[k];
    f /= LAL_TWOPI;
    f /= hp->data->data[k] * hp->data->data[k] + hc->data->data[k] * hc->data->data[k];
    if (f >= target) return k - 1;
  }
  XLAL_ERROR(XLAL_EDOM);
}

/* Return the index of the sample at with the peak amplitude */
static size_t find_peak_amp(const REAL8TimeSeries *hp, const REAL8TimeSeries *hc) {
  const REAL8 *hpdata = hp->data->data;
  const REAL8 *hcdata = hc->data->data;
  size_t k = hp->data->length;
  size_t peak_ind = -1;
  REAL8 peak_amp_sq = 0.;

  for (;k--;) {
    const REAL8 amp_sq = hpdata[k] * hpdata[k] + hcdata[k] * hcdata[k];
    if (amp_sq > peak_amp_sq) {
        peak_ind = k;
        peak_amp_sq = amp_sq;
    }
  }
  return peak_ind;
}

static int apply_phase_shift(const REAL8TimeSeries *hp, const REAL8TimeSeries *hc, const REAL8 shift) {
    REAL8 *hpdata = hp->data->data;
    REAL8 *hcdata = hc->data->data;
    size_t k = hp->data->length;
    const double cs = cos(shift);
    const double ss = sin(shift);

    for (;k--;) {
        const REAL8 temp_hpdata = hpdata[k] * cs - hcdata[k] * ss;
        hcdata[k] = hpdata[k] * ss + hcdata[k] * cs;
        hpdata[k] = temp_hpdata;
    }
    return 0;
}

static int apply_inclination(const REAL8TimeSeries *hp, const REAL8TimeSeries *hc, const REAL8 inclination) {
  REAL8 inclFacPlus, inclFacCross;
  REAL8 *hpdata = hp->data->data;
  REAL8 *hcdata = hc->data->data;
  size_t k = hp->data->length;

  inclFacCross = -cos(inclination);
  inclFacPlus = -0.5 * (1. + inclFacCross * inclFacCross);
  for (;k--;) {
      hpdata[k] *= inclFacPlus;
      hcdata[k] *= inclFacCross;
  }

  return XLAL_SUCCESS;
}
