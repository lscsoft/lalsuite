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
#include <lal/XLALError.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_complex.h>

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

/*
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

/* INTERNAL ROUTINES */

/* *******************************************************************/
/* Compute phenomenological parameters for non-spinning binaries     */
/* Ref. Eq.(4.18) of http://arxiv.org/pdf/0710.2335 and              */
/* Table I of http://arxiv.org/pdf/0712.0343                         */
/*                                                                   */
/* Takes solar masses.                                               */
/* *******************************************************************/
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

/*
 * Return tau0, the Newtonian chirp length estimate. Ref. Eq.(6) of B. S. Sathyaprakash, http://arxiv.org/abs/gr-qc/9411043v1
 */
static REAL8 ComputeTau0(const REAL8 m1, const REAL8 m2, const REAL8 f_min) {
  const REAL8 totalMass = m1 + m2;
  const REAL8 eta = m1 * m2 / (totalMass * totalMass);
  return 5. * totalMass * LAL_MTSUN_SI / (256. * eta * pow(LAL_PI * totalMass * LAL_MTSUN_SI * f_min, 8./3.));
}

/*
 * Estimate the length of a TD vector that can hold the waveform as the Newtonian
 * chirp time tau0 plus 1000 M.
 */
static size_t EstimateIMRLength(const REAL8 m1, const REAL8 m2, const REAL8 f_min, const REAL8 deltaT) {
  return (size_t) floor((ComputeTau0(m1, m2, f_min) + 1000 * (m1 + m2) * LAL_MTSUN_SI) / deltaT);
}

static size_t NextPow2(const size_t n) {
  return 1 << (size_t) ceil(log2(n));
}

/*
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

/*
 * Find a higher value of f_max so that we can safely apply a window later. 
 * The safety factor 1.025 is an empirical estimation 
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


/*
 * Private function to generate IMRPhenomA frequency-domain waveforms given coefficients
 */
static int IMRPhenomAGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< FD waveform */
    const REAL8 phi0,                  /**< orbital phase at peak (rad) */
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
  const REAL8 amp0 = - pow(LAL_MTSUN_SI*totalMass, 5./6.) * pow(fMerg, -7./6.)
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
    psiEff = - 2.*phi0
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
    data[i] += -I * ampEff * sin(psiEff);
  }

  return XLAL_SUCCESS;
}

/*
 * Private function to generate IMRPhenomB frequency-domain waveforms given coefficients
 */
static int IMRPhenomBGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< FD waveform */
    const REAL8 phi0,                  /**< orbital phase at peak (rad) */
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
  REAL8 amp0 = - pow(LAL_MTSUN_SI*totalMass, 5./6.) * pow(fMerg, -7./6.)
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
    psiEff = -2.*phi0
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

/*
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

/*
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

/*
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

  /* complex_h := hp + i hc := A e^{-i \phi}. Hence \phi = arctan(-hc/hp)
  and F = d \phi/dt = (hp hc' - hc hp') / (2 \pi A^2)
  We use first order finite difference to compute the derivatives hp' and hc'*/
  for (; k < n; k++) {
    const REAL8 hpDot = (hp->data->data[k+1] - hp->data->data[k-1]) / (2 * hp->deltaT);
    const REAL8 hcDot = (hc->data->data[k+1] - hc->data->data[k-1]) / (2 * hc->deltaT);
    REAL8 f = hcDot * hp->data->data[k] - hpDot * hc->data->data[k];
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

/* Apply the inclination-angle weighting to the two polarizations 
* Ref. Eq.(63) of B.S. Sathyaprakash and Bernard F. Schutz,
* “Physics, Astrophysics and Cosmology with Gravitational Waves”, 
* Living Rev. Relativity, 12, (2009), 2. [Online Article]: cited 27 Nov 2013,
* http://www.livingreviews.org/lrr-2009-2
*/
static int apply_inclination(const REAL8TimeSeries *hp, const REAL8TimeSeries *hc, const REAL8 inclination) {
  REAL8 inclFacPlus, inclFacCross;
  REAL8 *hpdata = hp->data->data;
  REAL8 *hcdata = hc->data->data;
  size_t k = hp->data->length;

  inclFacCross = cos(inclination);
  inclFacPlus = 0.5 * (1. + inclFacCross * inclFacCross);
  for (;k--;) {
      hpdata[k] *= inclFacPlus;
      hcdata[k] *= inclFacCross;
  }

  return XLAL_SUCCESS;
}

/* ***********************************************************************************/
/* END OF THE REVIEWED CODE: Below is the code for generating the IMRPhenomB metric */
/* ***********************************************************************************/

/*
 * Structure containing the coefficients of various powers of the Fourier frequency
 * f in the derivative of the amplitude and phase of IMRPhenomB w.r.t. the parameters
 * of the binary
 */
typedef struct tagIMRDerivCoeffs{
  REAL8 dA1dM_0; /*Coef of k = 0 term in (2k-7)/6 expansion of inspiral phase of amplitude*/
  REAL8 dA1dM_1; /*Coef of k = 2 term in (2k-7)/6 expansion of inspiral phase of amplitude*/
  REAL8 dA1dM_2; /*Coef of k = 3 term in (2k-7)/6 expansion of inspiral phase of amplitude*/
  REAL8 dA1deta_0;
  REAL8 dA1deta_1;
  REAL8 dA1deta_2;
  REAL8 dA1dchi_0;
  REAL8 dA1dchi_1;
  REAL8 dA1dchi_2;
  REAL8 dA2dM_0; /*Coef of k = 0 term in (k-2)/3 expansion of merger phase of amplitude*/
  REAL8 dA2dM_1; /*Coef of k = 1 term in (k-2)/3 expansion of merger phase of amplitude*/
  REAL8 dA2dM_2; /*Coef of k = 2 term in (k-2)/3 expansion of merger phase of amplitude*/
  REAL8 dA2deta_0;
  REAL8 dA2deta_1;
  REAL8 dA2deta_2;
  REAL8 dA2dchi_0;
  REAL8 dA2dchi_1;
  REAL8 dA2dchi_2;
  REAL8 dA3dMnum_0;
  REAL8 dA3dMnum_1;
  REAL8 dA3dMnum_2;
  REAL8 dA3detanum_0; /*Coef of k = 0 in pow(f,k) expansion of numerator of derivative of ringdown amplitude */
  REAL8 dA3detanum_1;
  REAL8 dA3detanum_2;
  REAL8 dA3dchinum_0;
  REAL8 dA3dchinum_1;
  REAL8 dA3dchinum_2;
  REAL8 dA3denom_0;
  REAL8 dA3denom_1;
  REAL8 dA3denom_2;
  REAL8 dA3denom_3;
  REAL8 dA3denom_4;
  REAL8 dPsidM_0; /*Coef of k = 0 in frequency series of Phase */
  REAL8 dPsidM_1;
  REAL8 dPsidM_2;
  REAL8 dPsidM_3;
  REAL8 dPsidM_4;
  REAL8 dPsidM_5;
  REAL8 dPsidM_6;
  REAL8 dPsidM_7;
  REAL8 dPsideta_0;
  REAL8 dPsideta_1;
  REAL8 dPsideta_2;
  REAL8 dPsideta_3;
  REAL8 dPsideta_4;
  REAL8 dPsideta_5;
  REAL8 dPsideta_6;
  REAL8 dPsideta_7;
  REAL8 dPsidchi_0;
  REAL8 dPsidchi_1;
  REAL8 dPsidchi_2;
  REAL8 dPsidchi_3;
  REAL8 dPsidchi_4;
  REAL8 dPsidchi_5;
  REAL8 dPsidchi_6;
  REAL8 dPsidchi_7;
  REAL8 Wm;
  REAL8 Wr;
  REAL8 dWmdM;
  REAL8 dWmdEta;
  REAL8 dWmdChi;
  REAL8 dWrdM;
  REAL8 dWrdEta;
  REAL8 dWrdChi;
  REAL8 alpha2;
  REAL8 alpha3;
  REAL8 eps1;
  REAL8 eps2;
}
IMRDerivCoeffs;

static gsl_matrix *XLALSimIMRPhenomBFisherMatrix(
    const REAL8 m1,     /**< component mass 1 (M_sun) */
    const REAL8 m2,     /**< component mass 2 (M_sun) */
    const REAL8 chi,    /**< effective spin parameter of IMRPhenomB: chi = (m1 chi1 + m2 chi2)/(m1+m2)  */
    const REAL8 fLow,   /**< low-frequency cutoff (Hz) */
    const REAL8FrequencySeries *Sh  /**< PSD in strain per root Hertz */
);

static IMRDerivCoeffs *XLALComputeIMRPhenomBDerivativeCoeffs(
    const REAL8 m1,   /**< component mass 1 (M_sun) */
    const REAL8 m2,   /**< component mass 2 (M_sun) */
    const REAL8 chi,  /**< effective spin parameter of IMRPhenomB: chi = (m1 chi1 + m2 chi2)/(m1+m2)  */
    BBHPhenomParams *params  /**<  phenomenological parameters of IMRPhenomB */
);

static gsl_matrix *XLALSimIMRPhenomBProjectExtrinsicParam(
  gsl_matrix *g  /**<  Fisher matrix of IMRPhenomB */
);

/*
 * Compute the coefficients of the series expansion (in Fourier frequency) of the derivatives of the
 * IMRPhenomB waveform with respect to parameters M, eta, chi
 */
static IMRDerivCoeffs *XLALComputeIMRPhenomBDerivativeCoeffs(
    const REAL8 m1,   /**< component mass 1 (M_sun) */
    const REAL8 m2,   /**< component mass 2 (M_sun) */
    const REAL8 chi,  /**< effective spin parameter of IMRPhenomB: chi = (m1 chi1 + m2 chi2)/(m1+m2)  */
    BBHPhenomParams *params  /**<  phenomenological parameters of IMRPhenomB */
){
  /* allocate memory for the structure IMRDerivCoeffs */
  IMRDerivCoeffs *derCoeffs = (IMRDerivCoeffs *) XLALMalloc(sizeof(IMRDerivCoeffs));
  if (!derCoeffs) XLAL_ERROR_NULL(XLAL_EFUNC);
  memset(derCoeffs, 0, sizeof(IMRDerivCoeffs));

  /* define mass and symmetric mass ratio */
  REAL8 M = (m1 + m2)*LAL_MTSUN_SI;
  REAL8 eta = m1*m2/((m1 + m2)*(m1 + m2));

  /* extract various phenomenological parameters */
  REAL8 f1 = params->fMerger;
  REAL8 f2 = params->fRing;
  REAL8 sigma = params->sigma;

  /*Define various terms*/
  REAL8 M_p2    = M*M;
  REAL8 M_p1by6 = sqrt(cbrt(M));
  REAL8 M_p5by6 = M/M_p1by6;
  REAL8 M_p7by6 = M*M_p1by6;
  REAL8 M_p1by3 = cbrt(M);
  REAL8 M_p2by3 = M_p1by3*M_p1by3;
  REAL8 sqrt_eta = sqrt(eta);
  REAL8 eta_p2 = eta*eta;
  REAL8 chi_x_eta = chi*eta;
  REAL8 chi_p2 = chi*chi;
  REAL8 f2_p2 = f2*f2;
  REAL8 f2_p3 = f2_p2*f2;
  REAL8 f2_p4 = f2_p3*f2;
  REAL8 f1_p7by6 = f1*sqrt(cbrt(f1));
  REAL8 f1_p3by2 = sqrt(f1*f1*f1);
  REAL8 f1_p3 = f1*f1*f1;
  REAL8 sigma_p2 = sigma*sigma;
  REAL8 sigma_p3 = sigma_p2*sigma;
  REAL8 f1M_p1by3 = cbrt(f1*M);
  REAL8 f1M_p2by3 = f1M_p1by3*f1M_p1by3;
  REAL8 f2M_p1by3 = cbrt(f2*M);
  REAL8 f2M_p2by3 = f2M_p1by3*f2M_p1by3;
  REAL8 f1byf2_p1by3 = cbrt(f1/f2);
  REAL8 f1byf2_p2by3 = f1byf2_p1by3*f1byf2_p1by3;

  /* derivatives w.r.t physical parameters of fMerger */
  REAL8 dF1dM = -f1/M;
  REAL8 dF1dEta = (0.3183098861837907*(0.64365 + 0.82696*chi - 0.27063*chi_p2 - 0.116436*eta - 7.8692*chi_x_eta - 21.2748*eta_p2))/M;
  REAL8 dF1dChi = (0.3183098861837907*(0.9666699/pow(1. - chi,0.783) - 0.91546/pow(1. - chi,0.74)))/M + (0.3183098861837907*(0.82696*eta - 0.54126*chi_x_eta - 3.9346*eta_p2))/M;

  /* derivatives w.r.t physical parameters of fRing */
  REAL8 dF2dM = -f2/M;
  REAL8 dF2dEta = (0.3183098861837907*(0.1469 - 0.12281*chi - 0.026091*chi_p2 - 0.0498*eta + 0.34026*chi_x_eta + 6.9756*eta_p2))/M;
  REAL8 dF2dChi = 0.03008028424436822/(pow(1. - chi,0.7)*M) + (0.3183098861837907*(-0.12281*eta - 0.052182*chi_x_eta + 0.17013*eta_p2))/M;

  /* derivatives of w.r.t physical parameters Sigma parameter */
  REAL8 dSigmadM = -sigma/M;
  REAL8 dSigmadEta = (0.3183098861837907*(-0.40979 - 0.035226*chi + 0.10082*chi_p2 + 3.6572*eta - 0.040338*chi_x_eta - 8.6094*eta_p2))/M;
  REAL8 dSigmadChi = (-0.03580986219567645*(1. - 0.63*pow(1. - chi,0.3)))/(pow(1. - chi,0.55)*M) + 0.01504014212218411/(pow(1. - chi,0.24999999999999994)*M) + (0.3183098861837907*(-0.035226*eta + 0.20164*chi_x_eta - 0.020169*eta_p2))/M;

  /* PN corrections to the frequency domain amplitude of the (2,2) mode */
  derCoeffs->alpha2   = -323./224. + 451.*eta/168.;
  derCoeffs->alpha3   = (27./8. - 11.*eta/6.)*chi;

  /* spin-dependent corrections to the merger amplitude */
  derCoeffs->eps1 =  1.4547*chi - 1.8897;
  derCoeffs->eps2 = -1.8153*chi + 1.6557;

  /* derivatives of PN corrections to the frequency domain amplitude of the (2,2) mode */
  REAL8 dAlpha2dEta = 2.6845238095238093;
  REAL8 dAlpha3dEta = -1.8333333333333333*chi;
  REAL8 dAlpha3dChi = 3.375 - 1.8333333333333333*eta;

  /* derivatives of spin-dependent corrections to the merger amplitude */
  REAL8 dEps1dChi = 1.4547;
  REAL8 dEps2dChi = -1.8153;

  /* Recurring expressions in normalization */
  REAL8 norm_expr_1 = (1. + 1.4645918875615231*derCoeffs->eps1*f2M_p1by3 + 2.1450293971110255*derCoeffs->eps2*f2M_p2by3);
  REAL8 norm_expr_2 = 0.33424583931521507*sqrt_eta*f1byf2_p2by3*M_p5by6;
  REAL8 norm_expr_3 = 1. + 1.4645918875615231*derCoeffs->eps1*f1M_p1by3 + 2.1450293971110255*derCoeffs->eps2*f1M_p2by3;
  REAL8 norm_expr_4 = sqrt_eta*f1byf2_p2by3*M_p5by6*norm_expr_1/f1_p7by6;
  REAL8 norm_expr_5 = 0.22283055954347672*sqrt_eta*M_p5by6*norm_expr_1*sigma;

  /* normalisation constant of the inspiral/merger amplitude */
  derCoeffs->Wm = (1. + 3.141592653589793*derCoeffs->alpha3*f1*M + 2.1450293971110255*derCoeffs->alpha2*f1M_p2by3)/(norm_expr_3);
  derCoeffs->Wr = derCoeffs->Wm*(norm_expr_2*norm_expr_1*sigma)/f1_p7by6;

  /* compute the derivatives of the phenomenological parameters w.r.t the physical parameters.*/
  REAL8 dWmDen = norm_expr_3*norm_expr_3; /* This is the common denominatior to dWmdM and dWmdeta dWmdchi */

  derCoeffs->dWmdM = 0.;
  derCoeffs->dWmdEta = ((-0.48819729585384103*dF1dEta*M*(derCoeffs->eps1 + 2.9291837751230463*derCoeffs->eps2*f1M_p1by3)*(1. + 3.141592653589793*derCoeffs->alpha3*f1*M + 2.1450293971110255*derCoeffs->alpha2*f1M_p2by3))/f1M_p2by3 + (3.141592653589793*derCoeffs->alpha3*dF1dEta*M + 3.141592653589793*dAlpha3dEta*f1*M + (1.430019598074017*derCoeffs->alpha2*dF1dEta*M)/f1M_p1by3 + 2.1450293971110255*dAlpha2dEta*f1M_p2by3)*(norm_expr_3))/dWmDen;
  derCoeffs->dWmdChi = ((3.141592653589793*derCoeffs->alpha3*dF1dChi*M + 3.141592653589793*dAlpha3dChi*f1*M + (1.430019598074017*derCoeffs->alpha2*dF1dChi*M)/f1M_p1by3)*(norm_expr_3) - (0.48819729585384103*M*(1. + 3.141592653589793*derCoeffs->alpha3*f1*M + 2.1450293971110255*derCoeffs->alpha2*f1M_p2by3)*(3.*f1*(dEps1dChi + 1.4645918875615231*dEps2dChi*f1M_p1by3) + dF1dChi*(derCoeffs->eps1 + 2.9291837751230463*derCoeffs->eps2*f1M_p1by3)))/f1M_p2by3)/dWmDen;

  derCoeffs->dWrdM = derCoeffs->Wm*((0.33424583931521507*dSigmadM*norm_expr_4) + (0.27853819942934593*sqrt_eta*f1byf2_p2by3*norm_expr_1*sigma)/(f1_p7by6*M_p1by6) + (norm_expr_5*((-1.*dF2dM*f1)/f2_p2 + dF1dM/f2))/(f1_p7by6*f1byf2_p1by3) - (0.38995347920108425*dF1dM*norm_expr_4*sigma)/(f1) + (norm_expr_2*((0.48819729585384103*derCoeffs->eps1*(f2 + dF2dM*M))/f2M_p2by3 + (1.430019598074017*derCoeffs->eps2*(f2 + dF2dM*M))/f2M_p1by3)*sigma)/f1_p7by6);
  derCoeffs->dWrdEta = derCoeffs->Wm*((0.33424583931521507*dSigmadEta*norm_expr_4) + (norm_expr_2*((0.48819729585384103*dF2dEta*derCoeffs->eps1*M)/f2M_p2by3 + (1.430019598074017*dF2dEta*derCoeffs->eps2*M)/f2M_p1by3)*sigma)/f1_p7by6 + (norm_expr_5*((-1.*dF2dEta*f1)/f2_p2 + dF1dEta/f2))/(f1_p7by6*f1byf2_p1by3) - (0.38995347920108425*dF1dEta*norm_expr_4*sigma)/(f1) + (0.16712291965760753*f1byf2_p2by3*M_p5by6*norm_expr_1*sigma)/(sqrt_eta*f1_p7by6)) + derCoeffs->Wr*derCoeffs->dWmdEta/derCoeffs->Wm;
  derCoeffs->dWrdChi = derCoeffs->Wm*((0.33424583931521507*dSigmadChi*norm_expr_4) + (norm_expr_2*((0.48819729585384103*dF2dChi*derCoeffs->eps1*M)/f2M_p2by3 + (1.430019598074017*dF2dChi*derCoeffs->eps2*M)/f2M_p1by3 + 1.4645918875615231*dEps1dChi*f2M_p1by3 + 2.1450293971110255*dEps2dChi*f2M_p2by3)*sigma)/f1_p7by6 + (norm_expr_5*((-1.*dF2dChi*f1)/f2_p2 + dF1dChi/f2))/(f1_p7by6*f1byf2_p1by3) - (0.38995347920108425*dF1dChi*norm_expr_4*sigma)/(f1)) + derCoeffs->Wr*derCoeffs->dWmdChi/derCoeffs->Wm;

  /* ------------------------------------------------------------------- */
  /* compute the coefficients of the series expansion of the derivatives */
  /* ------------------------------------------------------------------- */

  /* Coefficients of the derivatives of the inspiral amplitude */
  derCoeffs->dA1dM_0 = (0.1773229251163862*sqrt_eta)/M_p1by6;
  derCoeffs->dA1dM_1 = 0.6846531968814576*derCoeffs->alpha2*sqrt(eta*M);
  derCoeffs->dA1dM_2 = 1.2255680774891218*derCoeffs->alpha3*sqrt_eta*M_p5by6;

  derCoeffs->dA1deta_0 = (0.10639375506983172*M_p5by6)/sqrt_eta;
  derCoeffs->dA1deta_1 = (0.22821773229381923*(derCoeffs->alpha2 + 2.*dAlpha2dEta*eta)*sqrt(M*M*M))/sqrt_eta;
  derCoeffs->dA1deta_2 = (0.33424583931521507*(derCoeffs->alpha3 + 2.*dAlpha3dEta*eta)*M_p5by6*M)/sqrt_eta;

  derCoeffs->dA1dchi_0 = 0.;
  derCoeffs->dA1dchi_1 = 0.;
  derCoeffs->dA1dchi_2 = 0.6684916786304301*dAlpha3dChi*sqrt_eta*M_p5by6*M;

  /*Coefficients of the derivatives of the merger amplitude */
  derCoeffs->dA2dM_0 = (0.03546458502327724*sqrt_eta*(-3.*dF1dM*M*derCoeffs->Wm + f1*(6.*derCoeffs->dWmdM*M + 5.*derCoeffs->Wm)))/(f1_p3by2*M_p1by6);
  derCoeffs->dA2dM_1 = (0.05194114352082774*derCoeffs->eps1*sqrt_eta*M_p1by6*(-3.*dF1dM*M*derCoeffs->Wm + f1*(6.*derCoeffs->dWmdM*M + 7.*derCoeffs->Wm)))/f1_p3by2;
  derCoeffs->dA2dM_2 = (0.22821773229381923*derCoeffs->eps2*sqrt(eta*M)*(-1.*dF1dM*M*derCoeffs->Wm + f1*(2.*derCoeffs->dWmdM*M + 3.*derCoeffs->Wm)))/f1_p3by2;

  derCoeffs->dA2deta_0 = (0.10639375506983172*M_p5by6*(-1.*dF1dEta*eta*derCoeffs->Wm + f1*(2.*derCoeffs->dWmdEta*eta + derCoeffs->Wm)))/sqrt(eta*f1_p3);
  derCoeffs->dA2deta_1 = (0.1558234305624832*derCoeffs->eps1*M_p7by6*(-1.*dF1dEta*eta*derCoeffs->Wm + f1*(2.*derCoeffs->dWmdEta*eta + derCoeffs->Wm)))/sqrt(eta*f1_p3);
  derCoeffs->dA2deta_2 = (0.22821773229381923*derCoeffs->eps2*sqrt(M*M*M)*(-1.*dF1dEta*eta*derCoeffs->Wm + f1*(2.*derCoeffs->dWmdEta*eta + derCoeffs->Wm)))/sqrt(eta*f1_p3);

  derCoeffs->dA2dchi_0 = (0.10639375506983172*sqrt_eta*M_p5by6*(2.*derCoeffs->dWmdChi*f1 - 1.*dF1dChi*derCoeffs->Wm))/f1_p3by2;
  derCoeffs->dA2dchi_1 = (0.1558234305624832*sqrt_eta*M_p7by6*(-1.*dF1dChi*derCoeffs->eps1*derCoeffs->Wm + 2.*f1*(derCoeffs->dWmdChi*derCoeffs->eps1 + dEps1dChi*derCoeffs->Wm)))/f1_p3by2;
  derCoeffs->dA2dchi_2 = 0.22821773229381923*sqrt_eta*sqrt(M*M*M/(f1*f1*f1))*(-1.*dF1dChi*derCoeffs->eps2*derCoeffs->Wm + 2.*f1*(derCoeffs->dWmdChi*derCoeffs->eps2 + dEps2dChi*derCoeffs->Wm));

  /* Coefficients of the derivatives of the ringdown amplitude (numerator) */
  derCoeffs->dA3dMnum_0 = 8.*derCoeffs->dWrdM*f2_p2*sigma + 2.*derCoeffs->dWrdM*sigma_p3 + (8.*dSigmadM*f2_p2 - 16.*dF2dM*f2*sigma - 2.*dSigmadM*sigma_p2)*derCoeffs->Wr;
  derCoeffs->dA3dMnum_1 = -16.*derCoeffs->dWrdM*f2*sigma - (16.*dSigmadM*f2 - 16.*dF2dM*sigma)*derCoeffs->Wr;
  derCoeffs->dA3dMnum_2 = 8.*derCoeffs->dWrdM*sigma + 8.*dSigmadM*derCoeffs->Wr;

  derCoeffs->dA3detanum_0 = 8.*derCoeffs->dWrdEta*f2_p2*sigma + 2.*derCoeffs->dWrdEta*sigma_p3 + (8.*dSigmadEta*f2_p2 - 16.*dF2dEta*f2*sigma - 2.*dSigmadEta*sigma_p2)*derCoeffs->Wr;
  derCoeffs->dA3detanum_1 = -16.*derCoeffs->dWrdEta*f2*sigma - (16.*dSigmadEta*f2 - 16.*dF2dEta*sigma)*derCoeffs->Wr;
  derCoeffs->dA3detanum_2 = 8.*derCoeffs->dWrdEta*sigma + 8.*dSigmadEta*derCoeffs->Wr;

  derCoeffs->dA3dchinum_0 = 8.*derCoeffs->dWrdChi*f2_p2*sigma + 2.*derCoeffs->dWrdChi*sigma_p3 + (8.*dSigmadChi*f2_p2 - 16.*dF2dChi*f2*sigma - 2.*dSigmadChi*sigma_p2)*derCoeffs->Wr;
  derCoeffs->dA3dchinum_1 = -16.*derCoeffs->dWrdChi*f2*sigma - (16.*dSigmadChi*f2 - 16.*dF2dChi*sigma)*derCoeffs->Wr;
  derCoeffs->dA3dchinum_2 = 8.*derCoeffs->dWrdChi*sigma + 8.*dSigmadChi*derCoeffs->Wr;

  /* Coefficients of the derivatives of the ringdown amplitude (denominator) */
  derCoeffs->dA3denom_0 = 50.26548245743669*f2_p4 + 25.132741228718345*f2_p2*sigma_p2 + 3.141592653589793*sigma_p2*sigma_p2;
  derCoeffs->dA3denom_1 = -201.06192982974676*f2_p3 - 50.26548245743669*f2*sigma_p2;
  derCoeffs->dA3denom_2 = 301.59289474462014*f2_p2 + 25.132741228718345*sigma_p2;
  derCoeffs->dA3denom_3 = -201.06192982974676*f2;
  derCoeffs->dA3denom_4 = 50.26548245743669;

  /* Coefficients of the derivatives of the phase */
  REAL8 dPsi2dEta = -920.91 + 492.13*chi + 135.03*chi_p2 + 13483.8*eta - 2106.8*chi_x_eta - 40191.*eta_p2;
  REAL8 dPsi2dChi = 492.13*eta + 270.06*chi_x_eta - 1053.4*eta_p2;
  REAL8 dPsi3dEta = 17022. - 9565.9*chi - 2182.1*chi_p2 - 242740.*eta + 41504.*chi_x_eta + 715770.*eta_p2;
  REAL8 dPsi3dChi = 37.666666666666664 - 9565.9*eta - 4364.2*chi_x_eta + 20752.*eta_p2;
  REAL8 dPsi4dEta = -125440. + 75066.*chi + 13382.*chi_p2 + 1.74708e6*eta - 331460.*chi_x_eta - 5.0808e6*eta_p2;
  REAL8 dPsi4dChi = -101.25*chi + 75066.*eta + 26764.*chi_x_eta - 165730.*eta_p2;
  REAL8 dPsi6dEta = -889770. + 631020.*chi + 50676.*chi_p2 + 1.19616e7*eta - 2.8296e6*chi_x_eta - 3.383999999999999e7*eta_p2;
  REAL8 dPsi6dChi = 631020.*eta + 101352.*chi_x_eta - 1.4148e6*eta_p2;
  REAL8 dPsi7dEta = 869600. - 670980.*chi - 30082.*chi_p2 - 1.16758e7*eta + 3.029e6*chi_x_eta + 3.2673e7*eta_p2;
  REAL8 dPsi7dChi = -670980.*eta - 60164.*chi_x_eta + 1.5145e6*eta_p2;

  derCoeffs->dPsidM_0 = -0.005796647796902313/(eta*M_p2by3*M*M);
  derCoeffs->dPsidM_1 = 0.;
  derCoeffs->dPsidM_2 = (-0.007460387957432594*params->psi2)/(eta*M_p2);
  derCoeffs->dPsidM_3 = (-0.007284282453678307*params->psi3)/(eta*M_p5by6*M_p5by6);
  derCoeffs->dPsidM_4 = (-0.005334250494181998*params->psi4)/(eta*M_p1by3*M);
  derCoeffs->dPsidM_5 = 0.;
  derCoeffs->dPsidM_6 = (0.0114421241215744*params->psi6)/(eta*M_p2by3);
  derCoeffs->dPsidM_7 = (0.03351608432985977*params->psi7)/(eta*M_p1by3);

  derCoeffs->dPsideta_0 = -0.0034779886781413872/(eta_p2*M_p5by6*M_p5by6);
  derCoeffs->dPsideta_1 = 0.;
  derCoeffs->dPsideta_2 = (0.007460387957432594*(dPsi2dEta*eta - 1.*params->psi2))/(eta_p2*M);
  derCoeffs->dPsideta_3 = (0.010926423680517461*(dPsi3dEta*eta - 1.*params->psi3))/(eta_p2*M_p2by3);
  derCoeffs->dPsideta_4 = (0.016002751482545992*(dPsi4dEta*eta - 1.*params->psi4))/(eta_p2*M_p1by3);
  derCoeffs->dPsideta_5 = 0.;
  derCoeffs->dPsideta_6 = (0.0343263723647232*M_p1by3*(dPsi6dEta*eta - 1.*params->psi6))/eta_p2;
  derCoeffs->dPsideta_7 = (0.050274126494789656*M_p2by3*(dPsi7dEta*eta - 1.*params->psi7))/eta_p2;

  derCoeffs->dPsidchi_0 = 0.;
  derCoeffs->dPsidchi_1 = 0.;
  derCoeffs->dPsidchi_2 = (0.007460387957432594*dPsi2dChi)/(eta*M);
  derCoeffs->dPsidchi_3 = (0.010926423680517461*dPsi3dChi)/(eta*M_p2by3);
  derCoeffs->dPsidchi_4 = (0.016002751482545992*dPsi4dChi)/(eta*M_p1by3);
  derCoeffs->dPsidchi_5 = 0.;
  derCoeffs->dPsidchi_6 = (0.0343263723647232*dPsi6dChi*M_p1by3)/eta;
  derCoeffs->dPsidchi_7 = (0.050274126494789656*dPsi7dChi*M_p2by3)/eta;

  return derCoeffs;
};

/*
 * Compute the Fisher information matrix of the IMRPhenomB waveform in the
 * M, eta, chi, t0, phi0 parameter space, for an SNR=1/sqrt(2) signal.
 */

static gsl_matrix *XLALSimIMRPhenomBFisherMatrix(
    const REAL8 m1,     /**< component mass 1 (M_sun) */
    const REAL8 m2,     /**< component mass 2 (M_sun) */
    const REAL8 chi,    /**< effective spin parameter of IMRPhenomB: chi = (m1 chi1 + m2 chi2)/(m1+m2)  */
    const REAL8 fLow,   /**< low-frequency cutoff (Hz) */
    const REAL8FrequencySeries *Sh  /**< PSD in strain per root Hertz */
){

  IMRDerivCoeffs *coef;
  BBHPhenomParams *params;

  /* check inputs */
  if (m1 <= 0) XLAL_ERROR_NULL(XLAL_EDOM);
  if (m2 <= 0) XLAL_ERROR_NULL(XLAL_EDOM);
  if (chi <= -1.) XLAL_ERROR_NULL(XLAL_EDOM);
  if (chi >= 1.) XLAL_ERROR_NULL(XLAL_EDOM);
  if (fLow <= 0) XLAL_ERROR_NULL(XLAL_EDOM);

  /* define physical parameters and required coefficients */
  REAL8 M   = (m1 + m2)*LAL_MTSUN_SI;
  REAL8 eta = m1*m2/((m1+m2)*(m1+m2));
  REAL8 piM = LAL_PI * M;
  REAL8 AmpCoef = cbrt(1./(LAL_PI*LAL_PI))*sqrt(5.*eta/24.)*pow(M,5.0/6.0);

  /* compute the phenomenological parameters describing the IMRPhenomB waveform */
  params = ComputeIMRPhenomBParams(m1, m2, chi);

  /* extract various phenomenological parameters */
  REAL8 f1 = params->fMerger;
  REAL8 f2 = params->fRing;
  REAL8 sigma = params->sigma;
  REAL8 f3 = params->fCut;

  /* create a view of the PSD between flow and fCut */
  REAL8 df = Sh->deltaF;
  size_t nBins = (f3 - fLow) / df;

  /* compute the coefficients of the series expansion of the derivatives */
  coef = XLALComputeIMRPhenomBDerivativeCoeffs(m1, m2, chi, params);

  /* initialise the components of the Fisher matrix */
  REAL8 gamma_MM = 0.;
  REAL8 gamma_MEta= 0.;
  REAL8 gamma_MChi = 0.;
  REAL8 gamma_MT0 = 0.;
  REAL8 gamma_MPhi0 = 0.;
  REAL8 gamma_EtaEta = 0.;
  REAL8 gamma_EtaChi = 0.;
  REAL8 gamma_EtaT0 = 0.;
  REAL8 gamma_EtaPhi0= 0.;
  REAL8 gamma_ChiChi = 0.;
  REAL8 gamma_ChiT0 = 0.;
  REAL8 gamma_ChiPhi0 = 0.;
  REAL8 gamma_T0T0 = 0.;
  REAL8 gamma_T0Phi0 = 0.;
  REAL8 gamma_Phi0Phi0 = 0.;

  REAL8 amp, dAdM, dAdEta, dAdChi, dA3Den;
  REAL8 f1_pm7by6 = 1./(f1*cbrt(sqrt(f1)));
  REAL8 f1_p2by3 = cbrt(f1*f1);
  REAL8 sigma_p2 = sigma*sigma;
  REAL8 amp_merg_const = coef->Wm*AmpCoef*f1_pm7by6;
  REAL8 amp_ring_const = (coef->Wr*sigma)/LAL_TWOPI;

  /* compute derivatives over a freq vector from fLow to fCut with frequency resolution df
 *    *  and use this to compute the Fisher matrix */
  REAL8 f = fLow;

  for (size_t k=0; k<nBins; k++) {

    REAL8 v     = cbrt(piM*f);
    REAL8 v_p2    = v*v;
    REAL8 v_p3    = v_p2*v;
    REAL8 f_p2    = f*f;
    REAL8 f_p3    = f_p2*f;
    REAL8 f_p4    = f_p3*f;
    REAL8 f_p1by3 = cbrt(f);
    REAL8 f_p2by3 = f_p1by3*f_p1by3;

    REAL8 f_pm1   = 1./f;
    REAL8 f_pm1by2  = sqrt(f_pm1);
    REAL8 f_pm1by3  = 1./f_p1by3;
    REAL8 f_pm2by3  = 1./f_p2by3;

    REAL8 f_pm1by6  = cbrt(sqrt(f_pm1));
    REAL8 f_pm5by3  = f_pm2by3*f_pm1;
    REAL8 f_pm7by6  = f_pm1by6*f_pm1;

    REAL8 fbyf1_pm2by3 = f1_p2by3*f_pm2by3;

    /* compute derivatives of the amplitude w.r.t M, eta and chi */
    if (f <= f1) {
      amp   = AmpCoef*f_pm7by6*(1. + coef->alpha2*v_p2 + coef->alpha3*v_p3);
      dAdM  = coef->dA1dM_0*f_pm7by6 + coef->dA1dM_1*f_pm1by2 + coef->dA1dM_2*f_pm1by6;
      dAdEta  = coef->dA1deta_0*f_pm7by6 + coef->dA1deta_1*f_pm1by2 + coef->dA1deta_2*f_pm1by6;
      dAdChi  = coef->dA1dchi_0*f_pm7by6 + coef->dA1dchi_1*f_pm1by2 + coef->dA1dchi_2*f_pm1by6;
    }
    else if ((f1<f) && (f<=f2)) {
      amp   = amp_merg_const*fbyf1_pm2by3*(1. + coef->eps1*v + coef->eps2*v_p2);
      dAdM  = coef->dA2dM_0*f_pm2by3 + coef->dA2dM_1*f_pm1by3 + coef->dA2dM_2;
      dAdEta  = coef->dA2deta_0*f_pm2by3 + coef->dA2deta_1*f_pm1by3 + coef->dA2deta_2;
      dAdChi  = coef->dA2dchi_0*f_pm2by3 + coef->dA2dchi_1*f_pm1by3 + coef->dA2dchi_2;
    }
    else {
      amp   = amp_ring_const/((f-f2)*(f-f2) + sigma_p2*0.25);
      dA3Den  = coef->dA3denom_0 + coef->dA3denom_1*f + coef->dA3denom_2*f_p2 + coef->dA3denom_3*f_p3 + coef->dA3denom_4*f_p4;
      dAdM  = (coef->dA3dMnum_0 + coef->dA3dMnum_1*f + coef->dA3dMnum_2*f_p2)/dA3Den;
      dAdEta  = (coef->dA3detanum_0 + coef->dA3detanum_1*f + coef->dA3detanum_2*f_p2)/dA3Den;
      dAdChi  = (coef->dA3dchinum_0 + coef->dA3dchinum_1*f + coef->dA3dchinum_2*f_p2)/dA3Den;
    }

    /* compute derivatives of the phase w.r.t M, eta and chi */
    REAL8 dPsidM  = coef->dPsidM_0*f_pm5by3 + coef->dPsidM_2*f_pm1 + coef->dPsidM_3*f_pm2by3 + coef->dPsidM_4*f_pm1by3 + coef->dPsidM_6*f_p1by3 + coef->dPsidM_7*f_p2by3;
    REAL8 dPsidEta  = coef->dPsideta_0*f_pm5by3 + coef->dPsideta_2*f_pm1 + coef->dPsideta_3*f_pm2by3 + coef->dPsideta_4*f_pm1by3 + coef->dPsideta_6*f_p1by3 + coef->dPsideta_7*f_p2by3;
    REAL8 dPsidChi  = coef->dPsidchi_0*f_pm5by3 + coef->dPsidchi_2*f_pm1 + coef->dPsidchi_3*f_pm2by3 + coef->dPsidchi_4*f_pm1by3 + coef->dPsidchi_6*f_p1by3 + coef->dPsidchi_7*f_p2by3;
    REAL8 dPsidT0 = LAL_TWOPI * f;

    /* compute the elements of the Fisher matrix in M, eta, chi */
    REAL8 amp_p2 = amp*amp;
    REAL8 ampSqr_x_dPsiM = amp_p2*dPsidM;
    REAL8 ampSqr_x_dPsiEta = amp_p2*dPsidEta;
    REAL8 ampSqr_x_dPsidChi = amp_p2*dPsidChi;
    REAL8 psd_pm1 = 1./Sh->data->data[k];
    REAL8 ampSqr_by_Sh = amp_p2*psd_pm1;

    gamma_MM    += (ampSqr_x_dPsiM*dPsidM + dAdM*dAdM)*psd_pm1;
    gamma_MEta    += (ampSqr_x_dPsiM*dPsidEta + dAdM*dAdEta)*psd_pm1;
    gamma_MChi    += (ampSqr_x_dPsiM*dPsidChi + dAdM*dAdChi)*psd_pm1;
    gamma_MT0   += ampSqr_x_dPsiM*dPsidT0*psd_pm1;
    gamma_MPhi0   += ampSqr_x_dPsiM*psd_pm1;
    gamma_EtaEta  += (ampSqr_x_dPsiEta*dPsidEta + dAdEta*dAdEta)*psd_pm1;
    gamma_EtaChi  += (ampSqr_x_dPsiEta*dPsidChi + dAdEta*dAdChi)*psd_pm1;
    gamma_EtaT0   += ampSqr_x_dPsiEta*dPsidT0*psd_pm1;
    gamma_EtaPhi0 += ampSqr_x_dPsiEta*psd_pm1;
    gamma_ChiChi  += (ampSqr_x_dPsidChi*dPsidChi + dAdChi*dAdChi)*psd_pm1;
    gamma_ChiT0   += ampSqr_x_dPsidChi*dPsidT0*psd_pm1;
    gamma_ChiPhi0   += ampSqr_x_dPsidChi*psd_pm1;
    gamma_T0T0    += dPsidT0*dPsidT0*ampSqr_by_Sh;
    gamma_T0Phi0  += dPsidT0*ampSqr_by_Sh;
    gamma_Phi0Phi0  += ampSqr_by_Sh;

    f += df;
  }

  /* this is actually twice the h-square 2*||h||^2 */
  REAL8 hSqr = 2.*gamma_Phi0Phi0;

  /* fill the gsl matrix containing the Fisher matrix in (M, eta, chi, t0, phi0) */
  gsl_matrix * g = gsl_matrix_calloc (5, 5);

  gsl_matrix_set (g, 0,0, gamma_MM/hSqr);
  gsl_matrix_set (g, 0,1, gamma_MEta/hSqr);
  gsl_matrix_set (g, 0,2, gamma_MChi/hSqr);
  gsl_matrix_set (g, 0,3, gamma_MT0/hSqr);
  gsl_matrix_set (g, 0,4, gamma_MPhi0/hSqr);

  gsl_matrix_set (g, 1,0, gsl_matrix_get(g, 0,1));
  gsl_matrix_set (g, 1,1, gamma_EtaEta/hSqr);
  gsl_matrix_set (g, 1,2, gamma_EtaChi/hSqr);
  gsl_matrix_set (g, 1,3, gamma_EtaT0/hSqr);
  gsl_matrix_set (g, 1,4, gamma_EtaPhi0/hSqr);

  gsl_matrix_set (g, 2,0, gsl_matrix_get(g, 0,2));
  gsl_matrix_set (g, 2,1, gsl_matrix_get(g, 1,2));
  gsl_matrix_set (g, 2,2, gamma_ChiChi/hSqr);
  gsl_matrix_set (g, 2,3, gamma_ChiT0/hSqr);
  gsl_matrix_set (g, 2,4, gamma_ChiPhi0/hSqr);

  gsl_matrix_set (g, 3,0, gsl_matrix_get(g, 0,3));
  gsl_matrix_set (g, 3,1, gsl_matrix_get(g, 1,3));
  gsl_matrix_set (g, 3,2, gsl_matrix_get(g, 2,3));
  gsl_matrix_set (g, 3,3, gamma_T0T0/hSqr);
  gsl_matrix_set (g, 3,4, gamma_T0Phi0/hSqr);

  gsl_matrix_set (g, 4,0, gsl_matrix_get(g, 0,4));
  gsl_matrix_set (g, 4,1, gsl_matrix_get(g, 1,4));
  gsl_matrix_set (g, 4,2, gsl_matrix_get(g, 2,4));
  gsl_matrix_set (g, 4,3, gsl_matrix_get(g, 3,4));
  gsl_matrix_set (g, 4,4, gamma_Phi0Phi0/hSqr);

  return g;
};

/*
 * Project the Fisher matrix orthogonal to the extrnsic parameters t0 and phi0
 */
static gsl_matrix *XLALSimIMRPhenomBProjectExtrinsicParam(
  gsl_matrix *g   /**<  Fisher matrix of IMRPhenomB */
){
  int s = 0;

  /* Form submatrices g1, g2, g3, g4, defined as:
 *    *              g = [ g1 g2
 *       *                    g4 g3 ]                           */
  gsl_matrix_view g1v = gsl_matrix_submatrix (g, 0, 0, 3, 3);
  gsl_matrix_view g2v = gsl_matrix_submatrix (g, 0, 3, 3, 2);
  gsl_matrix_view g3v = gsl_matrix_submatrix (g, 3, 3, 2, 2);
  gsl_matrix_view g4v = gsl_matrix_submatrix (g, 3, 0, 2, 3);

  gsl_matrix * g1 = gsl_matrix_calloc (3, 3);
  gsl_matrix * g2 = gsl_matrix_calloc (3, 2);
  gsl_matrix * g3 = gsl_matrix_calloc (2, 2);
  gsl_matrix * g4 = gsl_matrix_calloc (2, 3);
  gsl_matrix * g3invg4 = gsl_matrix_calloc (2, 3);
  gsl_matrix * g2g3invg4 = gsl_matrix_calloc (3, 3);

  /* Project out the t0 and phi0 dimensions: gamma =  g1 - g2 g3^{-1} g4 */
  gsl_matrix * g3inv = gsl_matrix_calloc (2, 2);
  gsl_permutation * p = gsl_permutation_calloc (2);

  gsl_matrix_memcpy (g1, &g1v.matrix);
  gsl_matrix_memcpy (g2, &g2v.matrix);
  gsl_matrix_memcpy (g3, &g3v.matrix);
  gsl_matrix_memcpy (g4, &g4v.matrix);
  gsl_matrix_free (g);

  gsl_linalg_LU_decomp (g3, p, &s);
  gsl_linalg_LU_invert (g3, p, g3inv);
  gsl_permutation_free (p);
  gsl_matrix_free (g3);

  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, g3inv, g4,  0.0, g3invg4);
  gsl_matrix_free (g4);
  gsl_matrix_free (g3inv);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, g2, g3invg4,  0.0, g2g3invg4);
  gsl_matrix_free (g2);
  gsl_matrix_free (g3invg4);

  gsl_matrix_sub (g1, g2g3invg4);
  gsl_matrix_free (g2g3invg4);
  return g1;

};

/* EXTERNAL ROUTINES */

/**
 * @addtogroup LALSimIMRPhenom_c
 * @brief Routines to produce IMRPhenom-family of phenomenological
 * inspiral-merger-ringdown waveforms.
 *
 * @review IMRPhenomB routines reviewed by Frank Ohme, P. Ajith, Alex Nitz
 * and Riccardo Sturani. The review concluded with git hash
 * 43ce3b0a8753eb266d75a43ba94b6fb6412121d0 (May 2014).
 * @{
 */

/**
 * @name Routines for IMR Phenomenological Model "A"
 * @{
 */

/**
 * Driver routine to compute the non-spinning, inspiral-merger-ringdown
 * phenomenological waveform IMRPhenomA in the frequency domain.
 *
 * Reference:
 * - Waveform: Eq.(4.13) and (4.16) of http://arxiv.org/pdf/0710.2335
 * - Coefficients: Eq.(4.18) of http://arxiv.org/pdf/0710.2335 and
 * Table I of http://arxiv.org/pdf/0712.0343
 *
 * All input parameters should be SI units.
 */
int XLALSimIMRPhenomAGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< FD waveform */
    const REAL8 phi0,                  /**< orbital phase at peak (rad) */
    const REAL8 deltaF,                /**< frequency resolution (Hz) */
    const REAL8 m1_SI,                 /**< mass of companion 1 (kg) */
    const REAL8 m2_SI,                 /**< mass of companion 2 (kg) */
    const REAL8 f_min,                 /**< starting GW frequency (Hz) */
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
 * - Waveform: Eq.(4.13) and (4.16) of http://arxiv.org/pdf/0710.2335
 * - Coefficients: Eq.(4.18) of http://arxiv.org/pdf/0710.2335 and
 * Table I of http://arxiv.org/pdf/0712.0343
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

  /* generate hcross, which is hplus w/ GW phase shifted by -pi/2
   * <==> orb. phase shifted by -pi/4 */
  IMRPhenomAGenerateTD(hcross, -LAL_PI_4, deltaT, m1, m2, f_min, f_max_prime, distance, params);
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
 * Compute the default final frequency 
 */
double XLALSimIMRPhenomAGetFinalFreq(
    const REAL8 m1,
    const REAL8 m2
) {
    BBHPhenomParams *phenomParams;
    phenomParams = ComputeIMRPhenomAParams(m1, m2);
    return phenomParams->fCut;
}

/** @} */

/**
 * @name Routines for IMR Phenomenological Model "B"
 * @{
 */

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
 * Compute the default final frequency 
 */
double XLALSimIMRPhenomBGetFinalFreq(
    const REAL8 m1,
    const REAL8 m2,
    const REAL8 chi
) {
    BBHPhenomParams *phenomParams;
    phenomParams = ComputeIMRPhenomBParams(m1, m2, chi);
    return phenomParams->fCut;
}

/**
 * Driver routine to compute the spin-aligned, inspiral-merger-ringdown
 * phenomenological waveform IMRPhenomB in the time domain.
 *
 * Reference: http://arxiv.org/pdf/0909.2867
 * - Waveform: Eq.(1)
 * - Coefficients: Eq.(2) and Table I
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

  /* generate hcross, which is hplus w/ GW phase shifted by -pi/2
   * <==> orb. phase shifted by -pi/4 */
  IMRPhenomBGenerateTD(hcross, -LAL_PI_4, deltaT, m1, m2, chi, f_min, f_max_prime, distance, params);
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
 * - Waveform: Eq.(1)
 * - Coefficients: Eq.(2) and Table I
 *
 * All input parameters should be in SI units. Angles should be in radians.
 */
int XLALSimIMRPhenomBGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< FD waveform */
    const REAL8 phi0,                  /**< orbital phase at peak (rad) */
    const REAL8 deltaF,                /**< sampling interval (Hz) */
    const REAL8 m1_SI,                 /**< mass of companion 1 (kg) */
    const REAL8 m2_SI,                 /**< mass of companion 2 (kg) */
    const REAL8 chi,                   /**< mass-weighted aligned-spin parameter */
    const REAL8 f_min,                 /**< starting GW frequency (Hz) */
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

/**
 * Compute the template-space metric of the IMRPhenomB waveform in the
 * M, eta, chi coordinates
 */
int XLALSimIMRPhenomBMetricInMEtaChi(
    REAL8 *gamma00,  /**< template metric coeff. 00 in PN Chirp Time */
    REAL8 *gamma01,  /**< template metric coeff. 01/10 PN Chirp Time */
    REAL8 *gamma02,  /**< template metric coeff. 01/10 PN Chirp Time */
    REAL8 *gamma11,  /**< template metric coeff. 11 in PN Chirp Time */
    REAL8 *gamma12,  /**< template metric coeff. 01/10 PN Chirp Time */
    REAL8 *gamma22,  /**< template metric coeff. 01/10 PN Chirp Time */
    const REAL8 m1_SI,     /**< component mass 1 (kg) */
    const REAL8 m2_SI,     /**< component mass 2 (kg) */
    const REAL8 chi,    /**< effective spin parameter of IMRPhenomB: chi = (m1 chi1 + m2 chi2)/(m1+m2)  */
    const REAL8 fLow,   /**< low-frequency cutoff (Hz) */
    const REAL8FrequencySeries *Sh  /**< PSD in strain per root Hertz */
){

  /* external: SI; internal: solar masses */
  const REAL8 m1 = m1_SI / LAL_MSUN_SI;
  const REAL8 m2 = m2_SI / LAL_MSUN_SI;

  /* compute the Fisher matrix in (M, eta, chi, t0, phi0) coords */
  gsl_matrix * g = XLALSimIMRPhenomBFisherMatrix(m1, m2, chi, fLow, Sh);

  /* Project out t0 and phi0 */
  gsl_matrix * g1= XLALSimIMRPhenomBProjectExtrinsicParam(g);

  *gamma00 = gsl_matrix_get(g1, 0, 0);
  *gamma01 = gsl_matrix_get(g1, 0, 1);
  *gamma02 = gsl_matrix_get(g1, 0, 2);
  *gamma11 = gsl_matrix_get(g1, 1, 1);
  *gamma12 = gsl_matrix_get(g1, 1, 2);
  *gamma22 = gsl_matrix_get(g1, 2, 2);

  gsl_matrix_free (g1);
  return XLAL_SUCCESS;
};


/**
 * Compute the template-space metric of the IMRPhenomB waveform in the
 * theta0, theta3, theta3S coordinates
 */
int XLALSimIMRPhenomBMetricInTheta0Theta3Theta3S(
    REAL8 *gamma00,  /**< template metric coeff. 00 in PN Chirp Time */
    REAL8 *gamma01,  /**< template metric coeff. 01/10 PN Chirp Time */
    REAL8 *gamma02,  /**< template metric coeff. 01/10 PN Chirp Time */
    REAL8 *gamma11,  /**< template metric coeff. 11 in PN Chirp Time */
    REAL8 *gamma12,  /**< template metric coeff. 01/10 PN Chirp Time */
    REAL8 *gamma22,  /**< template metric coeff. 01/10 PN Chirp Time */
    const REAL8 m1_SI,     /**< component mass 1 (kg) */
    const REAL8 m2_SI,     /**< component mass 2 (kg) */
    const REAL8 chi,    /**< effective spin parameter of IMRPhenomB: chi = (m1 chi1 + m2 chi2)/(m1+m2)  */
    const REAL8 fLow,   /**< low-frequency cutoff (Hz) */
    const REAL8FrequencySeries *Sh  /**< PSD in strain per root Hertz */
){

  /* external: SI; internal: solar masses */
  const REAL8 m1 = m1_SI / LAL_MSUN_SI;
  const REAL8 m2 = m2_SI / LAL_MSUN_SI;

  const REAL8 M = (m1+m2)*LAL_MTSUN_SI;
  const REAL8 eta = m1*m2/((m1+m2)*(m1+m2));
  const REAL8 v0  = cbrt(LAL_PI*M*fLow);

  /* compute the Fisher matrix in (M, eta, chi, t0, phi0) coords */
  gsl_matrix * gMass = XLALSimIMRPhenomBFisherMatrix(m1, m2, chi, fLow, Sh);

  const REAL8 theta0 = 5.0/(128.0*eta*v0*v0*v0*v0*v0);
  const REAL8 theta3 = (LAL_PI/(4.0*eta*v0*v0));
  const REAL8 theta3S = (LAL_PI/(4.0*v0*v0))*(17022. - 9565.9*chi);
  const REAL8 theta3_p2 = theta3*theta3;
  const REAL8 theta0_p2 = theta0*theta0;
  const REAL8 theta0_p2by3 = cbrt(theta3_p2);

  /* Co-ordinate Derivatives for Jacobian matrix */
  const REAL8 dMdTheta0 = (-0.015831434944115277*theta3)/(fLow*theta0_p2);
  const REAL8 dMdTheta3 = 0.015831434944115277/(fLow*theta0);
  const REAL8 dMdTheta3S = 0.;
  const REAL8 dEtadTheta0 = 3.8715528021485643/(cbrt(theta0)*theta3*theta0_p2by3);
  const REAL8 dEtadTheta3 = (-9.678882005371412*cbrt(theta0_p2))/(theta3_p2*theta0_p2by3);
  const REAL8 dEtadTheta3S = 0.;
  const REAL8 dChidTheta0 = (0.000012000696668619612*theta0_p2by3*theta3S)/(theta0*cbrt(theta0_p2));
  const REAL8 dChidTheta3 = (-0.000012000696668619612*theta3S)/(cbrt(theta0_p2)*cbrt(theta3));
  const REAL8 dChidTheta3S = (-0.00001800104500292942*theta0_p2by3)/cbrt(theta0_p2);

  /* Define the Jacobian transformation matrix which would give the Fisher Matrix in theta co-ordinates */
  gsl_matrix * jaco = gsl_matrix_calloc (5, 5);
  gsl_matrix * gPre = gsl_matrix_calloc (5, 5);
  gsl_matrix * g    = gsl_matrix_calloc (5, 5);

  gsl_matrix_set (jaco, 0,0 , dMdTheta0);
  gsl_matrix_set (jaco, 0,1 , dEtadTheta0);
  gsl_matrix_set (jaco, 0,2 , dChidTheta0);
  gsl_matrix_set (jaco, 0,3 , 0.);
  gsl_matrix_set (jaco, 0,4 , 0.);

  gsl_matrix_set (jaco, 1,0 , dMdTheta3);
  gsl_matrix_set (jaco, 1,1 , dEtadTheta3);
  gsl_matrix_set (jaco, 1,2 , dChidTheta3);
  gsl_matrix_set (jaco, 1,3 , 0.);
  gsl_matrix_set (jaco, 1,4 , 0.);

  gsl_matrix_set (jaco, 2,0 , dMdTheta3S);
  gsl_matrix_set (jaco, 2,1 , dEtadTheta3S);
  gsl_matrix_set (jaco, 2,2 , dChidTheta3S);
  gsl_matrix_set (jaco, 2,3 , 0.);
  gsl_matrix_set (jaco, 2,4 , 0.);

  gsl_matrix_set (jaco, 3,0 , 0.);
  gsl_matrix_set (jaco, 3,1 , 0.);
  gsl_matrix_set (jaco, 3,2 , 0.);
  gsl_matrix_set (jaco, 3,3 , 1.);
  gsl_matrix_set (jaco, 3,4 , 0.);

  gsl_matrix_set (jaco, 4,0 , 0.);
  gsl_matrix_set (jaco, 4,1 , 0.);
  gsl_matrix_set (jaco, 4,2 , 0.);
  gsl_matrix_set (jaco, 4,3 , 0.);
  gsl_matrix_set (jaco, 4,4 , 1.);

  gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, gMass, jaco,  0.0, gPre);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, jaco, gPre,  0.0, g);
  gsl_matrix_free (jaco);
  gsl_matrix_free (gPre);

  /* Project out t0 and phi0 */
  gsl_matrix * g1  = XLALSimIMRPhenomBProjectExtrinsicParam(g);

  *gamma00 = gsl_matrix_get(g1, 0, 0);
  *gamma01 = gsl_matrix_get(g1, 0, 1);
  *gamma02 = gsl_matrix_get(g1, 0, 2);
  *gamma11 = gsl_matrix_get(g1, 1, 1);
  *gamma12 = gsl_matrix_get(g1, 1, 2);
  *gamma22 = gsl_matrix_get(g1, 2, 2);

  gsl_matrix_free (g1);
  return XLAL_SUCCESS;
};

/** @} */
/** @} */
