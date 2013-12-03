/*
 * Copyright (C) 2012 Prayush Kumar, Frank Ohme
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


/* The paper refered to here as the Main paper, is Phys. Rev. D 82, 064016 (2010)
 * */

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


/*********************************************************************/
/* This structure stores the PN coefficients used to calculate flux  */
/* and waveform amplitude, and Fourier phase. It also stores some    */
/* frequently used expressions which are constant during waveform    */
/* generation.                                                       */
/*********************************************************************/

typedef struct tagBBHPhenomCParams{
  REAL8 piM;
  REAL8 m_sec;

  REAL8 fmin;
  REAL8 fCut;
  REAL8 df;

  REAL8 f0;
  REAL8 f1;
  REAL8 f2;
  REAL8 d0;
  REAL8 d1;
  REAL8 d2;

  REAL8 afin;
  REAL8 fRingDown;
  REAL8 MfRingDown;
  REAL8 Qual;

  REAL8 pfaN;
  REAL8 pfa2;
  REAL8 pfa3;
  REAL8 pfa4;
  REAL8 pfa5;
  REAL8 pfa6;
  REAL8 pfa6log;
  REAL8 pfa7;

  REAL8 xdotaN;
  REAL8 xdota2;
  REAL8 xdota3;
  REAL8 xdota4;
  REAL8 xdota5;
  REAL8 xdota6;
  REAL8 xdota6log;
  REAL8 xdota7;

  REAL8 AN;
  REAL8 A2;
  REAL8 A3;
  REAL8 A4;
  REAL8 A5;
  REAL8 A5imag;
  REAL8 A6;
  REAL8 A6log;
  REAL8 A6imag;

  REAL8 a1;
  REAL8 a2;
  REAL8 a3;
  REAL8 a4;
  REAL8 a5;
  REAL8 a6;
  REAL8 g1;
  REAL8 del1;
  REAL8 del2;
  REAL8 b1;
  REAL8 b2;
}
BBHPhenomCParams;

/**
 * private function prototypes; all internal functions use solar masses.
 *
 */

static BBHPhenomCParams *ComputeIMRPhenomCParamsSPA( const REAL8 m1, const REAL8 m2, const REAL8 chi );
static BBHPhenomCParams *ComputeIMRPhenomCParams( const REAL8 m1, const REAL8 m2, const REAL8 chi );
static REAL8 wPlus( const REAL8 f, const REAL8 f0, const REAL8 d, const BBHPhenomCParams *params );
static REAL8 wMinus( const REAL8 f, const REAL8 f0, const REAL8 d, const BBHPhenomCParams *params );

static size_t NextPow2(const size_t n);
static REAL8 IMRPhenomCGeneratePhasePM( REAL8 f, REAL8 eta, const BBHPhenomCParams *params );
static int IMRPhenomCGenerateAmpPhase( REAL8 *amplitude, REAL8 *phasing, REAL8 f, REAL8 eta, const BBHPhenomCParams *params);

static int IMRPhenomCGenerateFD(COMPLEX16FrequencySeries **htilde, const REAL8 phi0, const REAL8 deltaF, const REAL8 m1, const REAL8 m2, const REAL8 f_min, const REAL8 f_max, const REAL8 distance, const BBHPhenomCParams *params);

/**
 * Internal functions
 */

/*********************************************************************/
/* Compute PN coefficients for non-precessing binaries               */
/* Ref. Eq.(A3)-(A5) of http://arxiv.org/pdf/1005.3306v3.pdf         */
/*                                                                   */
/*********************************************************************/
static BBHPhenomCParams *ComputeIMRPhenomCParamsSPA(
    const REAL8 m1, /**< mass of companion 1 (solar masses) */
    const REAL8 m2, /**< mass of companion 2 (solar masses) */
    const REAL8 chi) /**< Reduced spin of the binary, defined in the main paper */
{

  BBHPhenomCParams *p = (BBHPhenomCParams *) XLALMalloc(sizeof(BBHPhenomCParams));
  if (!p) XLAL_ERROR_NULL(XLAL_EFUNC);
  memset(p, 0, sizeof(BBHPhenomCParams));

  /* calculate the total mass and symmetric mass ratio */
  const REAL8 M = m1 + m2;
  const REAL8 eta = m1 * m2 / (M * M);
  const REAL8 piM = M * LAL_PI * LAL_MTSUN_SI;
  const REAL8 eta2 = eta*eta;
  const REAL8 chisum = chi + chi;
  const REAL8 chiprod = chi * chi;
  const REAL8 chi2 = chi * chi;

  // Store the total Mass of the system in seconds
  p->m_sec = M * LAL_MTSUN_SI;
  p->piM = piM;

  /* Calculate the PN phasing terms */
  p->pfaN = 3.0/(128.0 * eta);
  p->pfa2 = (3715./756.) + (55.*eta/9.0);
  p->pfa3 = -16.0*LAL_PI + (113./3.)*chi - 38.*eta*chisum/3.;
  p->pfa4 = (152.93365/5.08032) - 50.*chi2 + eta*(271.45/5.04 + 1.25*chiprod) +
        3085.*eta2/72.;
  p->pfa5 = LAL_PI*(386.45/7.56 - 65.*eta/9.) - 
        chi*(735.505/2.268 + 130.*eta/9.) + chisum*(1285.0*eta/8.1 + 170.*eta2/9.) - 
        10.*chi2*chi/3. + 10.*eta*chi*chiprod;
  p->pfa6 = 11583.231236531/4.694215680 - 640.0*LAL_PI*LAL_PI/3. -
        6848.0*LAL_GAMMA/21. - 684.8*log(64.)/6.3 +
        eta*(2255.*LAL_PI*LAL_PI/12. - 15737.765635/3.048192) +
        76.055*eta2/1.728 - (127.825*eta2*eta/1.296) +
        2920.*LAL_PI*chi/3. - (175. - 1490.*eta)*chi2/3. -
        (1120.*LAL_PI/3. - 1085.*chi/3.)*eta*chisum +
        (269.45*eta/3.36 - 2365.*eta2/6.)*chiprod;

  p->pfa6log = -6848./63.;
  
  p->pfa7 = LAL_PI*(770.96675/2.54016 + 378.515*eta/1.512 - 740.45*eta2/7.56) -
        chi*(20373.952415/3.048192 + 1509.35*eta/2.24 - 5786.95*eta2/4.32) +
        chisum*(4862.041225*eta/1.524096 + 1189.775*eta2/1.008 - 717.05*eta2*eta/2.16 - 830.*eta*chi2/3. + 35.*eta2*chiprod/3.) -
        560.*LAL_PI*chi2 + 20.*LAL_PI*eta*chiprod +
        chi2*chi*(945.55/1.68 - 85.*eta) + chi*chiprod*(396.65*eta/1.68 + 255.*eta2);

  /* Coefficients to calculate xdot, that comes in the fourier amplitude */
  p->xdotaN = 64.*eta/5.;
  p->xdota2 = -7.43/3.36 - 11.*eta/4.;
  p->xdota3 = 4.*LAL_PI - 11.3*chi/1.2 + 19.*eta*chisum/6.;
  p->xdota4 = 3.4103/1.8144 + 5*chi2 + eta*(13.661/2.016 - chiprod/8.) + 5.9*eta2/1.8;
  p->xdota5 = -LAL_PI*(41.59/6.72 + 189.*eta/8.) - chi*(31.571/1.008 - 116.5*eta/2.4) +
          chisum*(21.863*eta/1.008 - 79.*eta2/6.) - 3*chi*chi2/4. +
          9.*eta*chi*chiprod/4.;
  p->xdota6 = 164.47322263/1.39708800 - 17.12*LAL_GAMMA/1.05 +
          16.*LAL_PI*LAL_PI/3 - 8.56*log(16.)/1.05 +
          eta*(45.1*LAL_PI*LAL_PI/4.8 - 561.98689/2.17728) +
          5.41*eta2/8.96 - 5.605*eta*eta2/2.592 - 80.*LAL_PI*chi/3. +
          eta*chisum*(20.*LAL_PI/3. - 113.5*chi/3.6) +
          chi2*(64.153/1.008 - 45.7*eta/3.6) -
          chiprod*(7.87*eta/1.44 - 30.37*eta2/1.44);
  
  p->xdota6log = -856./105.;
  
  p->xdota7 = -LAL_PI*(4.415/4.032 - 358.675*eta/6.048 - 91.495*eta2/1.512) -
          chi*(252.9407/2.7216 - 845.827*eta/6.048 + 415.51*eta2/8.64) +
          chisum*(158.0239*eta/5.4432 - 451.597*eta2/6.048 + 20.45*eta2*eta/4.32 + 107.*eta*chi2/6. - 5.*eta2*chiprod/24.) +
          12.*LAL_PI*chi2 - chi2*chi*(150.5/2.4 + eta/8.) +
          chi*chiprod*(10.1*eta/2.4 + 3.*eta2/8.);

  /* Coefficients to compute the time-domain amplitude, which also enters the
   * fourier amplitude. */
  p->AN = 8.*eta*sqrt(LAL_PI/5.);
  p->A2 = (-107. + 55.*eta)/42.;
  p->A3 = 2.*LAL_PI - 4.*chi/3. + 2.*eta*chisum/3.;
  p->A4 = -2.173/1.512 - eta*(10.69/2.16 - 2.*chiprod) + 2.047*eta2/1.512;
  p->A5 = -10.7*LAL_PI/2.1 + eta*(3.4*LAL_PI/2.1);
  
  p->A5imag = -24.*eta;

  p->A6 = 270.27409/6.46800 - 8.56*LAL_GAMMA/1.05 +
      2.*LAL_PI*LAL_PI/3. +
      eta*(4.1*LAL_PI*LAL_PI/9.6 - 27.8185/3.3264) -
      20.261*eta2/2.772 + 11.4635*eta*eta2/9.9792 -
      4.28*log(16.)/1.05;
  
  p->A6log = -428./105.;

  p->A6imag = 4.28*LAL_PI/1.05;

  return p;
}

/*********************************************************************/
/* Compute phenomenological parameters for non-precessing binaries   */
/* Ref. Eq.(5.14) of http://arxiv.org/pdf/1005.3306v3.pdf  */
/* and Table II of the same paper.                                   */
/*                                                                   */
/*********************************************************************/

static BBHPhenomCParams *ComputeIMRPhenomCParams(
    const REAL8 m1, /**< mass of companion 1 (solar masses) */
    const REAL8 m2, /**< mass of companion 2 (solar masses) */
    const REAL8 chi) /**< Reduced spin of the binary, defined in the main paper */
{

  BBHPhenomCParams *p = NULL;
  p = ComputeIMRPhenomCParamsSPA( m1, m2, chi );
  if( !p )
    XLAL_ERROR_NULL(XLAL_EFUNC);

  const REAL8 M = m1 + m2;
  const REAL8 eta = m1 * m2 / (M * M);

  const REAL8 z101 = -2.417e-03;
  const REAL8 z102 = -1.093e-03;
  const REAL8 z111 = -1.917e-02;
  const REAL8 z110 = 7.267e-02;
  const REAL8 z120 = -2.504e-01;

  const REAL8 z201 = 5.962e-01;
  const REAL8 z202 = -5.600e-02;
  const REAL8 z211 = 1.520e-01;
  const REAL8 z210 = -2.970e+00;
  const REAL8 z220 = 1.312e+01;

  const REAL8 z301 = -3.283e+01;
  const REAL8 z302 = 8.859e+00;
  const REAL8 z311 = 2.931e+01;
  const REAL8 z310 = 7.954e+01;
  const REAL8 z320 = -4.349e+02;

  const REAL8 z401 = 1.619e+02;
  const REAL8 z402 = -4.702e+01;
  const REAL8 z411 = -1.751e+02;
  const REAL8 z410 = -3.225e+02;
  const REAL8 z420 = 1.587e+03;

  const REAL8 z501 = -6.320e+02;
  const REAL8 z502 = 2.463e+02;
  const REAL8 z511 = 1.048e+03;
  const REAL8 z510 = 3.355e+02;
  const REAL8 z520 = -5.115e+03;

  const REAL8 z601 = -4.809e+01;
  const REAL8 z602 = -3.643e+02;
  const REAL8 z611 = -5.215e+02;
  const REAL8 z610 = 1.870e+03;
  const REAL8 z620 = 7.354e+02;

  const REAL8 z701 = 4.149e+00;
  const REAL8 z702 = -4.070e+00;
  const REAL8 z711 = -8.752e+01;
  const REAL8 z710 = -4.897e+01;
  const REAL8 z720 = 6.665e+02;
  
  const REAL8 z801 = -5.472e-02;
  const REAL8 z802 = 2.094e-02;
  const REAL8 z811 = 3.554e-01;
  const REAL8 z810 = 1.151e-01;
  const REAL8 z820 = 9.640e-01;

  const REAL8 z901 = -1.235e+00;
  const REAL8 z902 = 3.423e-01;
  const REAL8 z911 = 6.062e+00;
  const REAL8 z910 = 5.949e+00;
  const REAL8 z920 = -1.069e+01;

  /* Calculate alphas, gamma, deltas from Table II and Eq 5.14 of Main paper */
  REAL8 eta2 = eta*eta;
  REAL8 chi2 = chi*chi;
  REAL8 etachi = eta * chi;

  p->a1 = z101 * chi + z102 * chi2 + z111 * eta * chi + z110 * eta + z120 * eta2;
  p->a2 = z201 * chi + z202 * chi2 + z211 * eta * chi + z210 * eta + z220 * eta2;
  p->a3 = z301 * chi + z302 * chi2 + z311 * eta * chi + z310 * eta + z320 * eta2;
  p->a4 = z401 * chi + z402 * chi2 + z411 * eta * chi + z410 * eta + z420 * eta2;
  p->a5 = z501 * chi + z502 * chi2 + z511 * eta * chi + z510 * eta + z520 * eta2;
  p->a6 = z601 * chi + z602 * chi2 + z611 * eta * chi + z610 * eta + z620 * eta2;

  p->g1 = z701 * chi + z702 * chi2 + z711 * eta * chi + z710 * eta + z720 * eta2;

  p->del1 = z801 * chi + z802 * chi2 + z811 * eta * chi + z810 * eta + z820 * eta2;
  p->del2 = z901 * chi + z902 * chi2 + z911 * eta * chi + z910 * eta + z920 * eta2;

  /* Get the Spin of the final BH */
  REAL8 s4 = -0.129;
  REAL8 s5 = -0.384;
  REAL8 t0 = -2.686;
  REAL8 t2 = -3.454;
  REAL8 t3 = 2.353;
  REAL8 finspin = chi + s4*chi*etachi + s5*etachi*eta + t0*etachi + 2.*sqrt(3.)*eta + 
    t2*eta2 + t3*eta2*eta;

  if( fabs(finspin) > 1.0 )
    XLAL_ERROR_NULL( XLAL_EDOM );

  p->afin = finspin;
  
  /* Get the Ringdown frequency */
  REAL8 prefac = (1./(2.*LAL_PI)) * LAL_C_SI * LAL_C_SI * LAL_C_SI / (LAL_G_SI * M * LAL_MSUN_SI);
  REAL8 k1 = 1.5251;
  REAL8 k2 = -1.1568;
  REAL8 k3 = 0.1292;
  //printf("fRD prefac = %12.18f\n", prefac);

  p->fRingDown = (prefac * (k1 + k2 * pow(1. - fabs(finspin), k3)));
  p->MfRingDown = p->m_sec * p->fRingDown;

  /* Get the quality factor of ring-fown, using Eq (5.6) of Main paper */
  p->Qual = (0.7000 + (1.4187 * pow(1.0 - fabs(finspin), -0.4990)) );;

  /* Get the transition frequencies, at which the model switches phase and
   * amplitude prescriptions, as used in Eq.(5.9), (5.13) of the Main paper */
  p->f1 = 0.1 * p->fRingDown;
  p->f2 = p->fRingDown;
  p->f0 = 0.98 * p->fRingDown;
  p->d1 = 0.005;
  p->d2 = 0.005;
  p->d0 = 0.015;

  /* Get the coefficients beta1, beta2, defined in Eq 5.7 of the main paper */
  REAL8 Mfrd = p->MfRingDown;

  p->b2 = ((-5./3.)* p->a1 * pow(Mfrd,(-8./3.)) - p->a2/(Mfrd*Mfrd) - 
      (p->a3/3.)*pow(Mfrd,(-4./3.)) + (2./3.)* p->a5 * pow(Mfrd,(-1./3.)) + p->a6)/eta;

  REAL8 psiPMrd = IMRPhenomCGeneratePhasePM( p->fRingDown, eta, p );

  p->b1 = psiPMrd - (p->b2 * Mfrd);

  /* Taking the upper cut-off frequency as 0.15M */
  //p->fCut = (1.7086 * eta * eta - 0.26592 * eta + 0.28236) / p->piM;
  p->fCut = 0.15 / p->m_sec;

  return p;

}

/*********************************************************************/
/* The following function return the hyperbolic-Tan+ windows used   */
/* in Eq.(5.9), (5.13) of the Main paper                             */
/*********************************************************************/
static REAL8 wPlus( const REAL8 f, const REAL8 f0, const REAL8 d, const BBHPhenomCParams *params )
{

  REAL8 Mf = params->m_sec * f;
  REAL8 Mf0 = params->m_sec * f0;

  return ( 0.5 * (1. + tanh(4.*(Mf - Mf0)/d) ) );

}

/*********************************************************************/
/* The following function return the hyperbolic-Tan- windows used   */
/* in Eq.(5.9), (5.13) of the Main paper                             */
/*********************************************************************/
static REAL8 wMinus( const REAL8 f, const REAL8 f0, const REAL8 d, const BBHPhenomCParams *params )
{
 
  REAL8 Mf = params->m_sec * f;
  REAL8 Mf0 = params->m_sec * f0;

  return ( 0.5 * (1. - tanh(4.*(Mf - Mf0)/d) ) );

}

/*********************************************************************/
/* The following function return the closest higher power of 2       */
/*********************************************************************/
static size_t NextPow2(const size_t n) {
  return 1 << (size_t) ceil(log2(n));
}

/**
 * main functions
 *
 */


/**
 * Driver routine to compute the spin-aligned, inspiral-merger-ringdown
 * phenomenological waveform IMRPhenomC in the frequency domain.
 *
 * Reference: http://arxiv.org/pdf/1005.3306v3.pdf
 * - Waveform: Eq.(5.3)-(5.13)
 * - Coefficients: Eq.(5.14) and Table II
 *
 * All input parameters should be in SI units. Angles should be in radians.
 */
int XLALSimIMRPhenomCGenerateFD(
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
  BBHPhenomCParams *params;
  int status;
  REAL8 f_max_prime;

  /* external: SI; internal: solar masses */
  const REAL8 m1 = m1_SI / LAL_MSUN_SI;
  const REAL8 m2 = m2_SI / LAL_MSUN_SI;

  /* check inputs for sanity */
  if (*htilde) XLAL_ERROR(XLAL_EFAULT);
  if (deltaF <= 0) XLAL_ERROR(XLAL_EDOM);
  if (m1 <= 0) XLAL_ERROR(XLAL_EDOM);
  if (m2 <= 0) XLAL_ERROR(XLAL_EDOM);
  if (fabs(chi) > 1) XLAL_ERROR(XLAL_EDOM);
  if (f_min <= 0) XLAL_ERROR(XLAL_EDOM);
  if (f_max < 0) XLAL_ERROR(XLAL_EDOM);
  if (distance <= 0) XLAL_ERROR(XLAL_EDOM);

  /* If spins are above 0.9 or below -0.9, throw an error */
  if (chi > 0.9 || chi < -0.9){
      XLALPrintError("Spins outside the range [-0.9,0.9] are not supported\n");
      XLAL_ERROR(XLAL_EDOM);
  }

  /* If mass ratio is above 4 and below 20, give a warning, and if it is above
   * 20, throw an error */
  REAL8 q = m1/m2;
  if (q > 1.0)
    q = 1./q;

  if (q > 20.0){
      XLALPrintError("Mass ratio is way outside the calibration range. m1/m2 should be <= 20.\n");
      XLAL_ERROR(XLAL_EDOM);
  }
  else if (q > 4.0){
      XLALPrintWarning("Warning: The model is only calibrated for m1/m2 <= 4.\n");
  }

  /* phenomenological parameters*/
  params = ComputeIMRPhenomCParams(m1, m2, chi);
  if (!params) XLAL_ERROR(XLAL_EFUNC);
  if (params->fCut <= f_min) {
      XLALPrintError("(fCut = 0.15M) <= f_min\n");
      XLAL_ERROR(XLAL_EDOM);
  }

  /* default f_max to params->fCut */
  f_max_prime = f_max ? f_max : params->fCut;
  f_max_prime = (f_max_prime > params->fCut) ? params->fCut : f_max_prime;
  if (f_max_prime <= f_min) {
      XLALPrintError("f_max <= f_min\n");
      XLAL_ERROR(XLAL_EDOM);
  }

  status = IMRPhenomCGenerateFD(htilde, phi0, deltaF, m1, m2, f_min, f_max_prime, distance, params);
  LALFree(params);
  return status;
}

/***********************************************************************************/
/* The following private function generates the Pre-Merger phase, as defined in    */
/* Eq. (5.3), (5.9) of the Main paper.                                             */
/***********************************************************************************/

static REAL8 IMRPhenomCGeneratePhasePM( 
    REAL8 f,                       /**< frequency (Hz) */ 
    REAL8 eta,                     /**< dimensionless mass-ratio */
    const BBHPhenomCParams *params /**< pointer to Object storing coefficients and constants */ 
    )
{
  REAL8 Mf = params->m_sec * f;
  REAL8 w = pow( Mf, (1./3.) );
  REAL8 w2 = w*w;
  REAL8 w3 = w2*w;
  REAL8 w5 = w3*w2;

  REAL8 phasing = (params->a1/w5) + (params->a2/w3) + (params->a3/w) + params->a4 + 
    (params->a5*w2) +(params->a6*w3);
  phasing /= eta;

  return phasing;
}


/***********************************************************************************/
/* The following private function generates the complete amplitude and phase of    */
/* PhenomC waveform, at a given frequency.                                         */
/* Eq. (5.3), (5.9) of the Main paper.                                             */
/***********************************************************************************/

static int IMRPhenomCGenerateAmpPhase( 
    REAL8 *amplitude, /**< pointer to memory for phenomC amplitude */
    REAL8 *phasing,   /**< pointer to memory for phenomC phase */ 
    REAL8 f,          /**< frequency (Hz) */ 
    REAL8 eta,        /**< dimensionless mass-ratio */
    const BBHPhenomCParams *params /**< pointer to Object storing coefficients and constants */
    )
{
  *amplitude = 0.0;
  *phasing = 0.0;

  /* Get the phase */
  REAL8 v =  cbrt(params->piM * f);
  REAL8 Mf = params->m_sec * f;

  if( v >= 1.0 )
    XLAL_ERROR(XLAL_EDOM);
  
  REAL8 v2 = v*v;
  REAL8 v3 = v*v2;
  REAL8 v4 = v2*v2;
  REAL8 v5 = v3*v2;
  REAL8 v6 = v3*v3;
  REAL8 v7 = v4*v3;
  REAL8 v10 = v5*v5;

  /* SPA part of the phase */
  REAL8 phSPA = 1. + params->pfa2 * v2 + params->pfa3 * v3 + params->pfa4 * v4 + 
    (1. + log(v3)) * params->pfa5 * v5 + (params->pfa6  + params->pfa6log * log(v3))*v6 + 
    params->pfa7 * v7;
  phSPA *= (params->pfaN / v5);

  // Taking t0 = phi0 = 0
  phSPA -= (LAL_PI / 4.);

  REAL8 w = cbrt( Mf );
  REAL8 w2 = w*w;
  REAL8 w3 = w2*w;
  REAL8 w5 = w3*w2;

  /* The Pre-Merger (PM) phase */
  REAL8 phPM = (params->a1/w5) + (params->a2/w3) + (params->a3/w) + params->a4 + 
    (params->a5*w2) +(params->a6*w3);
  phPM /= eta;

  /* Ring-down phase */
  REAL8 phRD = params->b1 + params->b2 * params->m_sec * f;
    
  REAL8 wPlusf1 = wPlus( f, params->f1, params->d1, params );
  REAL8 wPlusf2 = wPlus( f, params->f2, params->d2, params );
  REAL8 wMinusf1 = wMinus( f, params->f1, params->d1, params );
  REAL8 wMinusf2 = wMinus( f, params->f2, params->d2, params );
    
  *phasing = phSPA * wMinusf1 + phPM * wPlusf1 * wMinusf2 + phRD * wPlusf2;   
  
  /* Get the amplitude */
  REAL8 xdot = 1. + params->xdota2*v2 + params->xdota3*v3 + params->xdota4*v4 + 
    params->xdota5*v5 + (params->xdota6 + params->xdota6log*log(v2))*v6 + 
    params->xdota7 * v7;
  xdot *= (params->xdotaN * v10);

  if( xdot < 0.0 && f < params->f1 )
  {
    XLALPrintError("omegaDot < 0, while frequency is below SPA-PM matching freq.");
    XLAL_ERROR( XLAL_EDOM );
  }
    
  REAL8 aPM = 0.0;
  
  /* Following Emma's code, take only the absolute value of omegaDot, when
   * computing the amplitude */
  REAL8 omgdot = 1.5*v*xdot;
  REAL8 ampfac = sqrt(fabs(LAL_PI/omgdot));

  /* Get the real and imaginary part of the PM amplitude */
  REAL8 AmpPMre = ampfac * params->AN * v2 * (1. + params->A2*v2 + params->A3*v3 + 
      params->A4*v4 + params->A5*v5 + (params->A6 + params->A6log*log(v2))*v6);
  REAL8 AmpPMim = ampfac * params->AN * v2 * (params->A5imag * v5 + params->A6imag * v6);
   
  /* Following Emma's code, we take the absolute part of the complex SPA
   * amplitude, and keep that as the amplitude */
  aPM = sqrt( AmpPMre * AmpPMre + AmpPMim * AmpPMim );
  
  aPM += (params->g1 * pow(Mf,(5./6.)));

  /* The Ring-down aamplitude */
  REAL8 Mfrd = params->MfRingDown;

  /* From Emma's code, use sigma = fRD * del2 / Qual */
  REAL8 sig = Mfrd * params->del2 / params->Qual;
  REAL8 sig2 = sig*sig;
  REAL8 L = sig2 / ((Mf - Mfrd)*(Mf - Mfrd) + sig2/4.);

  REAL8 aRD = params->del1 * L * pow(Mf, (-7./6.));

  REAL8 wPlusf0 = wPlus( f, params->f0, params->d0, params );
  REAL8 wMinusf0 = wMinus( f, params->f0, params->d0, params );
   
  *amplitude = - (aPM * wMinusf0 + aRD * wPlusf0);

  return XLAL_SUCCESS;
}


/***********************************************************************************/
/* The following private function generates IMRPhenomC frequency-domain waveforms  */
/* given coefficients */
/***********************************************************************************/

static int IMRPhenomCGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< FD waveform */
    const REAL8 phi0,                  /**< phase at peak */
    const REAL8 deltaF,                /**< frequency resolution */
    const REAL8 m1,                    /**< mass of companion 1 [solar masses] */
    const REAL8 m2,                    /**< mass of companion 2 [solar masses] */
    //const REAL8 chi,                   /**< mass-weighted aligned-spin parameter */
    const REAL8 f_min,                 /**< start frequency */
    const REAL8 f_max,                 /**< end frequency */
    const REAL8 distance,              /**< distance to source (m) */
    const BBHPhenomCParams *params      /**< from ComputeIMRPhenomCParams */
) {
  static LIGOTimeGPS ligotimegps_zero = {0, 0};
  size_t i;
  INT4 errcode;

  const REAL8 M = m1 + m2;
  const REAL8 eta = m1 * m2 / (M * M);
  
  /* Memory to temporarily store components of amplitude and phase */
  /*
  REAL8 phSPA, phPM, phRD, aPM, aRD;
  REAL8 wPlusf1, wPlusf2, wMinusf1, wMinusf2, wPlusf0, wMinusf0;
  */
  REAL8 phPhenomC = 0.0;
  REAL8 aPhenomC = 0.0;

  /* compute the amplitude pre-factor */
  REAL8 amp0 = 2. * sqrt(5. / (64.*LAL_PI)) * M * LAL_MRSUN_SI * M * LAL_MTSUN_SI / distance;

  /* allocate htilde */
  size_t n = NextPow2(f_max / deltaF) + 1;
  *htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &ligotimegps_zero, 0.0, 
      deltaF, &lalStrainUnit, n);
  memset((*htilde)->data->data, 0, n * sizeof(COMPLEX16));
  XLALUnitMultiply(&((*htilde)->sampleUnits), &((*htilde)->sampleUnits), &lalSecondUnit);
  if (!(*htilde)) XLAL_ERROR(XLAL_EFUNC);

  /* now generate the waveform */
  size_t ind_max = (size_t) (f_max / deltaF);
  for (i = (size_t) (f_min / deltaF); i < ind_max; i++) 
  {
    REAL8 f = i * deltaF;
    
    errcode = IMRPhenomCGenerateAmpPhase( &aPhenomC, &phPhenomC, f, eta, params );
    if( errcode != XLAL_SUCCESS )
      XLAL_ERROR(XLAL_EFUNC);
    phPhenomC -= 2.*phi0; // factor of 2 b/c phi0 is orbital phase

    /* generate the waveform */
    ((*htilde)->data->data)[i] = amp0 * aPhenomC * cos(phPhenomC);
    ((*htilde)->data->data)[i] += -I * amp0 * aPhenomC * sin(phPhenomC);
  }

  return XLAL_SUCCESS;
}



/*********************************************************************************/
/* The following functions can be used to generate the individual components the */
/* PhenomC waveform, both amplitude and phase.                                   */
/*********************************************************************************/

//static size_t find_instant_freq(const REAL8TimeSeries *hp, const REAL8TimeSeries *hc, const REAL8 target, const size_t start);
//static size_t find_peak_amp(const REAL8TimeSeries *hp, const REAL8TimeSeries *hc);
//static int apply_phase_shift(const REAL8TimeSeries *hp, const REAL8TimeSeries *hc, const REAL8 shift);
//static int apply_inclination(const REAL8TimeSeries *hp, const REAL8TimeSeries *hc, const REAL8 inclination);

/*
static REAL8 IMRPhenomCGeneratePhaseSPA( REAL8 f, const BBHPhenomCParams *params );
static REAL8 IMRPhenomCGeneratePhaseRD( REAL8 f, const BBHPhenomCParams *params );
static REAL8 IMRPhenomCGenerateAmplitudePM( REAL8 f, const BBHPhenomCParams *params );
static REAL8 IMRPhenomCGenerateAmplitudeRD( REAL8 f, const BBHPhenomCParams *params );
*/

/*
static REAL8 IMRPhenomCGeneratePhaseSPA( REAL8 f, const BBHPhenomCParams *params )
{
  REAL8 v =  cbrt(params->piM * f);

  if( v >= 1.0 )
    XLAL_ERROR(XLAL_EDOM);
  
  REAL8 v2 = v*v;
  REAL8 v3 = v*v2;
  REAL8 v4 = v2*v2;
  REAL8 v5 = v3*v2;
  REAL8 v6 = v3*v3;
  REAL8 v7 = v4*v3;
  
  REAL8 phasing = 1. + params->pfa2 * v2 + params->pfa3 * v3 + params->pfa4 * v4 + 
    (1. + log(v3)) * params->pfa5 * v5 + (params->pfa6  + params->pfa6log * log(v3))*v6 + 
    params->pfa7 * v7;
  phasing *= (params->pfaN / v5);

  // Taking t0 = phi0 = 0
  phasing -= (LAL_PI / 4.);

  return phasing;
}
*/

/*
static REAL8 IMRPhenomCGeneratePhaseRD( REAL8 f, const BBHPhenomCParams *params )
{
 return (params->b1 + params->b2 * params->m_sec * f);
}
*/

#if 0
static REAL8 IMRPhenomCGenerateAmplitudePM( REAL8 f, const BBHPhenomCParams *params )
{
  REAL8 v =  cbrt(params->piM * f);
  REAL8 Mf = params->m_sec * f;

  if( v >= 1.0 )
    XLAL_ERROR( XLAL_EDOM );
  
  REAL8 v2 = v*v;
  REAL8 v3 = v*v2;
  REAL8 v4 = v2*v2;
  REAL8 v5 = v3*v2;
  REAL8 v6 = v3*v3;
  REAL8 v7 = v4*v3;
  REAL8 v10 = v5*v5;

  REAL8 xdot = 1. + params->xdota2*v2 + params->xdota3*v3 + params->xdota4*v4 + 
    params->xdota5*v5 + (params->xdota6 + params->xdota6log*log(v2))*v6 + 
    params->xdota7 * v7;
  xdot *= (params->xdotaN * v10);

  if( xdot < 0.0 && f < params->fRingDown )
  {
    XLALPrintError("omegaDot < 0, while frequency is below RingDown freq.");
    XLAL_ERROR( XLAL_EDOM );
  }
    
  REAL8 amplitude = 0.0;
  if( xdot > 0.0 )
  {
    REAL8 omgdot = 1.5*v*xdot;
    REAL8 ampfac = sqrt(LAL_PI/omgdot);

    /* Get the real and imaginary part of the PM amplitude */
    REAL8 AmpPMre = ampfac * params->AN * v2 * (1. + params->A2*v2 + params->A3*v3 + 
        params->A4*v4 + params->A5*v5 + (params->A6 + params->A6log*log(v2))*v6);
    REAL8 AmpPMim = ampfac * params->AN * v2 * (params->A5imag * v5 + params->A6imag * v6);
   
    /* Following Emma's code, we take the absolute part of the complex SPA
    * amplitude, and keep that as the amplitude */
    amplitude = sqrt( AmpPMre * AmpPMre + AmpPMim * AmpPMim );
  }

  amplitude += (params->g1 * pow(Mf,(5./6.)));

  return amplitude;
}

static REAL8 IMRPhenomCGenerateAmplitudeRD( REAL8 f, const BBHPhenomCParams *params )
{

  REAL8 Mf = params->m_sec * f;
  REAL8 Mfrd = params->MfRingDown;

  /* From Emma's code, use sigma = fRD * del2 / Qual */
  REAL8 sig = Mfrd * params->del2 / params->Qual;
  REAL8 sig2 = sig*sig;
  REAL8 L = sig2 / ((Mf - Mfrd)*(Mf - Mfrd) + sig2/4.);

  REAL8 amplitude = params->del1 * L * pow(Mf, (-7./6.));

  return amplitude;
 
}
#endif

#if 0
/* To be put in the function IMRPhenomCGenerateFD
 *
 * */
/* Get the phase */
    /*
    phSPA = IMRPhenomCGeneratePhaseSPA( f, params );
    phPM = IMRPhenomCGeneratePhasePM( f, eta, params );
    phRD = IMRPhenomCGeneratePhaseRD( f, params );
    
    wPlusf1 = wPlus( f, params->f1, params->d1, params );
    wPlusf2 = wPlus( f, params->f2, params->d2, params );
    wMinusf1 = wMinus( f, params->f1, params->d1, params );
    wMinusf2 = wMinus( f, params->f2, params->d2, params );
    
    phPhenomC = phi0 + phSPA * wMinusf1 + phPM * wPlusf1 * wMinusf2 + phRD * wPlusf2;   
    */
    /* Get the amplitude */
    /*
    aPM = IMRPhenomCGenerateAmplitudePM( f, params );
    aRD = IMRPhenomCGenerateAmplitudeRD( f, params );

    wPlusf0 = wPlus( f, params->f0, params->d0, params );
    wMinusf0 = wMinus( f, params->f0, params->d0, params );
    
    aPhenomC = aPM * wMinusf0 + aRD * wPlusf0;
    */
#endif

