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

#include <lal/LALStdlib.h>
#include <lal/LALSimIMR.h>
#include <lal/LALComplex.h>
#include <lal/LALConstants.h>
#include <lal/FrequencySeries.h>
#include <lal/BBHPhenomCoeffs.h>
#include <lal/Units.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

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
 * private function prototypes
 *
 */

static BBHPhenomParams *ComputeIMRPhenomAParams(REAL8 m1, REAL8 m2);
static UNUSED BBHPhenomParams *ComputeIMRPhenomBParams(REAL8 m1, REAL8 m2, REAL8 chi);

static REAL8 EstimateSafeFMinForTD(REAL8 m1, REAL8 m2, REAL8 f_min);
static UNUSED REAL8 EstimateSafeFMMaxForTD(REAL8 f_max, REAL8 dt);
static ssize_t EstimateFDLength(REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 dt);
static UNUSED ssize_t EstimateFDLengthForTD(REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 dt);

static REAL8 LorentzianFn(REAL8 freq, REAL8 fRing, REAL8 sigma);

static int IMRPhenomGenerateFD(COMPLEX16FrequencySeries **htilde, LIGOTimeGPS *tRef, REAL8 phiRef, REAL8 fRef, REAL8 deltaF, REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 f_max, REAL8 distance, BBHPhenomParams *params);

/*********************************************************************/
/* Compute phenomenological parameters for non-spinning binaries     */
/* Ref. Eq.(4.18) of http://arxiv.org/pdf/0710.2335 and              */
/* Table I of http://arxiv.org/pdf/0712.0343                         */
/*                                                                   */
/* Takes solar masses.                                               */
/*********************************************************************/
static BBHPhenomParams *ComputeIMRPhenomAParams(REAL8 m1, REAL8 m2) {
  REAL8 totalMass, piM, eta, fMerg_a, fMerg_b, fMerg_c, fRing_a, fRing_b, etap2;
  REAL8 fRing_c, sigma_a, sigma_b, sigma_c, fCut_a, fCut_b, fCut_c;
  REAL8 psi0_a, psi0_b, psi0_c, psi2_a, psi2_b, psi2_c, psi3_a, psi3_b, psi3_c;
  REAL8 psi4_a, psi4_b, psi4_c, psi6_a, psi6_b, psi6_c, psi7_a, psi7_b, psi7_c;
  BBHPhenomParams *phenParams;

  phenParams = (BBHPhenomParams *) XLALMalloc(sizeof(BBHPhenomParams));
  if (!phenParams) XLAL_ERROR_NULL(__func__, XLAL_EFUNC);
  memset(phenParams, 0, sizeof(BBHPhenomParams));

  /* calculate the total mass and symmetric mass ratio */
  totalMass = m1 + m2;
  eta = m1 * m2 / (totalMass * totalMass);
  piM = totalMass * LAL_PI * LAL_MTSUN_SI;

  fMerg_a = BBHPHENOMCOEFFSH_FMERG_A;
  fMerg_b = BBHPHENOMCOEFFSH_FMERG_B;
  fMerg_c = BBHPHENOMCOEFFSH_FMERG_C;

  fRing_a = BBHPHENOMCOEFFSH_FRING_A;
  fRing_b = BBHPHENOMCOEFFSH_FRING_B;
  fRing_c = BBHPHENOMCOEFFSH_FRING_C;

  sigma_a = BBHPHENOMCOEFFSH_SIGMA_A;
  sigma_b = BBHPHENOMCOEFFSH_SIGMA_B;
  sigma_c = BBHPHENOMCOEFFSH_SIGMA_C;

  fCut_a = BBHPHENOMCOEFFSH_FCUT_A;
  fCut_b = BBHPHENOMCOEFFSH_FCUT_B;
  fCut_c = BBHPHENOMCOEFFSH_FCUT_C;

  psi0_a = BBHPHENOMCOEFFSH_PSI0_X;
  psi0_b = BBHPHENOMCOEFFSH_PSI0_Y;
  psi0_c = BBHPHENOMCOEFFSH_PSI0_Z;

  psi2_a = BBHPHENOMCOEFFSH_PSI2_X;
  psi2_b = BBHPHENOMCOEFFSH_PSI2_Y;
  psi2_c = BBHPHENOMCOEFFSH_PSI2_Z;

  psi3_a = BBHPHENOMCOEFFSH_PSI3_X;
  psi3_b = BBHPHENOMCOEFFSH_PSI3_Y;
  psi3_c = BBHPHENOMCOEFFSH_PSI3_Z;

  psi4_a = BBHPHENOMCOEFFSH_PSI4_X;
  psi4_b = BBHPHENOMCOEFFSH_PSI4_Y;
  psi4_c = BBHPHENOMCOEFFSH_PSI4_Z;

  psi6_a = BBHPHENOMCOEFFSH_PSI6_X;
  psi6_b = BBHPHENOMCOEFFSH_PSI6_Y;
  psi6_c = BBHPHENOMCOEFFSH_PSI6_Z;

  psi7_a = BBHPHENOMCOEFFSH_PSI7_X;
  psi7_b = BBHPHENOMCOEFFSH_PSI7_Y;
  psi7_c = BBHPHENOMCOEFFSH_PSI7_Z;

  /* Evaluate the polynomials. See Eq. (4.18) of P. Ajith et al
   * arXiv:0710.2335 [gr-qc] */
  etap2 = eta*eta;
  phenParams->fCut = (fCut_a*etap2  + fCut_b*eta  + fCut_c)/piM;
  phenParams->fMerger = (fMerg_a*etap2  + fMerg_b*eta  + fMerg_c)/piM;
  phenParams->fRing = (fRing_a*etap2 + fRing_b*eta + fRing_c)/piM;
  phenParams->sigma = (sigma_a*etap2 + sigma_b*eta + sigma_c)/piM;

  phenParams->psi0 = (psi0_a*etap2 + psi0_b*eta + psi0_c)/(eta*pow(piM, 5./3.));
  phenParams->psi1 = 0.;
  phenParams->psi2 = (psi2_a*etap2 + psi2_b*eta + psi2_c)/(eta*pow(piM, 3./3.));
  phenParams->psi3 = (psi3_a*etap2 + psi3_b*eta + psi3_c)/(eta*pow(piM, 2./3.));
  phenParams->psi4 = (psi4_a*etap2 + psi4_b*eta + psi4_c)/(eta*pow(piM, 1./3.));
  phenParams->psi5 = 0.;
  phenParams->psi6 = (psi6_a*etap2 + psi6_b*eta + psi6_c)/(eta*pow(piM, -1./3.));
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
static BBHPhenomParams *ComputeIMRPhenomBParams(REAL8 m1, REAL8 m2, REAL8 chi) {
  REAL8 totalMass, piM, eta;
  REAL8 etap2, chip2, etap3, etap2chi, etachip2, etachi;
  BBHPhenomParams *phenParams;

  phenParams = (BBHPhenomParams *) XLALMalloc(sizeof(BBHPhenomParams));
  if (!phenParams) XLAL_ERROR_NULL(__func__, XLAL_EFUNC);
  memset(phenParams, 0, sizeof(BBHPhenomParams));

  /* calculate the total mass and symmetric mass ratio */
  totalMass = m1 + m2;
  eta = m1 * m2 / (totalMass * totalMass);
  piM = totalMass * LAL_PI * LAL_MTSUN_SI;

  /* spinning phenomenological waveforms */
  etap2 = eta*eta;
  chip2 = chi*chi;
  etap3 = etap2*eta;
  etap2chi = etap2*chi;
  etachip2 = eta*chip2;
  etachi = eta*chi;

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

static ssize_t EstimateFDLength(REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 dt) {
  REAL8 totalMass, eta, tau0;

  totalMass = m1 + m2;
  eta = m1 * m2 / (totalMass * totalMass);
  tau0 = 5. * totalMass * LAL_MTSUN_SI / (256. * eta * pow(LAL_PI * totalMass * LAL_MTSUN_SI * f_min, 8./3.));
  return (1 << (int) ceil(log2(tau0 * dt))) + 1;
}

/**
 * Find a lower value for f_min (using the definition of Newtonian chirp
 * time) such that the waveform has a minimum length of tau0. This is
 * necessary to avoid FFT artifacts.
 */
static REAL8 EstimateSafeFMinForTD(REAL8 m1, REAL8 m2, REAL8 f_min) {
  REAL8 temp_f_min, totalMass, eta, tau0;

  totalMass = m1 + m2;
  eta = m1 * m2 / (totalMass * totalMass);
  tau0 = 64.;
  temp_f_min = pow((tau0 * 256. * eta * pow(totalMass * LAL_MTSUN_SI, 5./3.) / 5.), -3./8.) / LAL_PI;
  if (temp_f_min > f_min) temp_f_min = f_min;
  if (temp_f_min < 0.5) temp_f_min = 0.5;
  return temp_f_min;
}

/**
 * Find a higher value of f_max so that we can safely apply a window later.
 */
static REAL8 EstimateSafeFMMaxForTD(REAL8 f_max, REAL8 dt) {
  REAL8 temp_f_max;
  temp_f_max = 1.025 * f_max;

  /* make sure that these frequencies are not too out of range */
  if (temp_f_max > dt / 2. - 100.) temp_f_max = dt / 2. - 100.;
  return temp_f_max;
}

/**
 * We will generate the waveform from a frequency which is lower than the
 * f_min chosen. Also the cutoff frequency is higher than the f_max. We
 * will later apply a window function, and truncate the time-domain waveform
 * below an instantaneous frequency f_min.
 */
static ssize_t EstimateFDLengthForTD(REAL8 m1, REAL8 m2, REAL8 f_min, REAL8 dt) {
  return EstimateFDLength(m1, m2, EstimateSafeFMinForTD(m1, m2, f_min), dt);
}

static REAL8 LorentzianFn (
    REAL8 freq,
    REAL8 fRing,
    REAL8 sigma) {
  return sigma / (2 * LAL_PI * ((freq - fRing)*(freq - fRing)
    + sigma*sigma / 4.0));
}

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
    LIGOTimeGPS *tRef,                 /**< time at fRef */
    REAL8 phiRef,                      /**< phase at fRef */
    REAL8 fRef,                        /**< reference frequency */
    REAL8 deltaF,                      /**< frequency resolution */
    REAL8 m1,                          /**< mass of companion 1 */
    REAL8 m2,                          /**< mass of companion 2 */
    REAL8 f_min,                       /**< start frequency */
    REAL8 f_max,                       /**< end frequency; if 0, set to fCut */
    REAL8 distance                     /**< distance of source */
) {
  BBHPhenomParams *params;

  /* check inputs for sanity */
  if (*htilde) XLAL_ERROR(__func__, XLAL_EFAULT);
  if (!tRef) XLAL_ERROR(__func__, XLAL_EFAULT);
  if (fRef <= 0) XLAL_ERROR(__func__, XLAL_EDOM);
  if (deltaF <= 0) XLAL_ERROR(__func__, XLAL_EDOM);
  if (m1 < 0) XLAL_ERROR(__func__, XLAL_EDOM);
  if (m2 < 0) XLAL_ERROR(__func__, XLAL_EDOM);
  if (f_min <= 0) XLAL_ERROR(__func__, XLAL_EDOM);
  if (f_max <= 0) XLAL_ERROR(__func__, XLAL_EDOM);
  if (distance <= 0) XLAL_ERROR(__func__, XLAL_EDOM);

  /* external: SI; internal: solar masses */
  m1 /= LAL_MSUN_SI;
  m2 /= LAL_MSUN_SI;

  /* phenomenological parameters*/
  params = ComputeIMRPhenomAParams(m1, m2);
  if (!params) XLAL_ERROR(__func__, XLAL_EFUNC);

  return IMRPhenomGenerateFD(htilde, tRef, phiRef, fRef, deltaF, m1, m2, f_min, f_max, distance, params);
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
    LIGOTimeGPS *tRef,                 /**< time at fRef */
    REAL8 phiRef,                      /**< phase at fRef */
    REAL8 fRef,                        /**< reference frequency */
    REAL8 deltaF,                      /**< sampling interval */
    REAL8 m1,                          /**< mass of companion 1 */
    REAL8 m2,                          /**< mass of companion 2 */
    REAL8 chi,                         /**< mass-weighted aligned-spin parameter */
    REAL8 f_min,                       /**< start frequency */
    REAL8 f_max,                       /**< end frequency */
    REAL8 distance                     /**< distance of source */
) {
  BBHPhenomParams *params;

  /* check inputs for sanity */
  if (*htilde) XLAL_ERROR(__func__, XLAL_EFAULT);
  if (!tRef) XLAL_ERROR(__func__, XLAL_EFAULT);
  if (fRef <= 0) XLAL_ERROR(__func__, XLAL_EDOM);
  if (deltaF <= 0) XLAL_ERROR(__func__, XLAL_EDOM);
  if (m1 < 0) XLAL_ERROR(__func__, XLAL_EDOM);
  if (m2 < 0) XLAL_ERROR(__func__, XLAL_EDOM);
  if (fabs(chi) > 1) XLAL_ERROR(__func__, XLAL_EDOM);
  if (f_min <= 0) XLAL_ERROR(__func__, XLAL_EDOM);
  if (f_max <= 0) XLAL_ERROR(__func__, XLAL_EDOM);
  if (distance <= 0) XLAL_ERROR(__func__, XLAL_EDOM);

  /* external: SI; internal: solar masses */
  m1 /= LAL_MSUN_SI;
  m2 /= LAL_MSUN_SI;

  /* phenomenological parameters*/
  params = ComputeIMRPhenomBParams(m1, m2, chi);
  if (!params) XLAL_ERROR(__func__, XLAL_EFUNC);

  return IMRPhenomGenerateFD(htilde, tRef, phiRef, fRef, deltaF, m1, m2, f_min, f_max, distance, params);
}


/**
 * Private function to generate the waveforms given coefficients
 */
static int IMRPhenomGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< FD waveform */
    LIGOTimeGPS *tRef,                 /**< time at fRef */
    REAL8 phiRef,                      /**< phase at fRef */
    REAL8 fRef,                        /**< reference frequency */
    REAL8 deltaF,                      /**< frequency resolution */
    REAL8 m1,                          /**< mass of companion 1 [solar masses] */
    REAL8 m2,                          /**< mass of companion 2 [solar masses] */
    REAL8 f_min,                       /**< start frequency */
    REAL8 f_max,                       /**< end frequency; if 0, set to fCut */
    REAL8 distance,                    /**< distance of source */
    BBHPhenomParams *params            /**< from ComputeIMRPhenom{A,B}Params */
) {
  REAL8 shft, amp0, ampEff, psiEff, fMerg, fRing, sigma, totalMass, eta;
  ssize_t i, n;

  fMerg = params->fMerger;
  fRing = params->fRing;
  sigma = params->sigma;
  totalMass = m1 + m2;
  eta = m1 * m2 / (totalMass * totalMass);

  /* compute the amplitude pre-factor */
  amp0 = pow(LAL_MTSUN_SI*totalMass, 5./6.) * pow(fMerg, -7./6.)
    / pow(LAL_PI, 2./3.) * sqrt(5. * eta / 24.) / (distance / LAL_C_SI);

  /* allocate htilde */
  /* When f_max is 0, use fCut to determine how many samples to allocate */
  n = (ssize_t) ceil((f_max || params->fCut) / deltaF);
  *htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", tRef, 0.0, deltaF, &lalStrainUnit, n);
  XLALUnitDivide((*htilde)->sampleUnits, (*tilde)->sampleUnits, &lalSecondUnit);
  if (!(*htilde)) XLAL_ERROR(__func__, XLAL_EFUNC);

  shft = LAL_TWOPI * (tRef->gpsSeconds + 1e-9 * tRef->gpsNanoSeconds);

  /* fill the zero and Nyquist frequency with zeros */
  ((*htilde)->data->data)[0] = (COMPLEX16) {0., 0.};
  ((*htilde)->data->data)[n-1] = (COMPLEX16) {0., 0.};

  /* now generate the waveform at all frequency bins */
  for (i=1; i < n - 1; i++) {
    /* Fourier frequency corresponding to this bin */
    REAL8 f = i * deltaF;
    REAL8 fNorm = f / fMerg;

    /* compute the amplitude */
    if ((f < f_min) || (f > params->fCut)) continue;
    else if (f <= fMerg) ampEff = amp0 * pow(fNorm, -7./6.);
    else if ((f > fMerg) & (f <= fRing)) ampEff = amp0 * pow(fNorm, -2./3.);
    else if (f > fRing)
      ampEff = amp0 * LAL_PI_2 * pow(fRing / fMerg, -2./3.) * sigma
        * LorentzianFn(f, fRing, sigma);
    else {
      XLALDestroyCOMPLEX16FrequencySeries(*htilde);
      *htilde = NULL;
      XLAL_ERROR(__func__, XLAL_EDOM);
    }

    /* now compute the phase */
    psiEff = shft * (f - fRef) + phiRef  /* use reference freq. and phase */
      + params->psi0 * pow(f, -5./3.)
      + params->psi1 * pow(f, -4./3.)
      + params->psi2 * pow(f, -3./3.)
      + params->psi3 * pow(f, -2./3.)
      + params->psi4 * pow(f, -1./3.)
      + params->psi5 // * pow(f, 0.)
      + params->psi6 * pow(f, 1./3.)
      + params->psi7 * pow(f, 2./3.);

    /* generate the waveform */
    ((*htilde)->data->data)[i].re = ampEff * cos(psiEff);
    ((*htilde)->data->data)[i].im = ampEff * sin(psiEff);
  }

  return XLAL_SUCCESS;
}
