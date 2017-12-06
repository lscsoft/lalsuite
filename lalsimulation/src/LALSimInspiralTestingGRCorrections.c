/*
 *  Copyright (C) 2017 Walter Del Pozzo
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
#include <gsl/gsl_spline.h>
#include <stdlib.h>
#include <math.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimInspiralTestingGRCorrections.h>
#include <lal/Sequence.h>
#include "LALSimInspiralPNCoefficients.c"
#include <lal/LALSimIMR.h>

/*
 * Copy and pasted from LALSimIMREOBNQCCorrection.c, which is now defunct.
 *
 * Peak frequency predicted by fitting NR results (currently only 2,2 available).
 * Unpublished. Used in building SEOBNRv4.
 */
UNUSED static inline REAL8 GetNRSpinPeakOmegaV4 (INT4 UNUSED l,
                                                 INT4 UNUSED m,
                                                 REAL8 UNUSED eta,
                                                 REAL8 a) /** Combined spin chi defined in Eq. (2.8) of Bohe et al. (2017)  */
{
    REAL8 chi = a;
    REAL8 res;
    res = 0.5626787200433265 + (-0.08706198756945482 +
                                0.0017434519312586804 * chi) *
    log (10.26207326082448 -
         chi * (7.629921628648589 -
                72.75949266353584 * (-0.25 + eta)) -
         62.353217004599784 * (-0.25 + eta));
    return res;
}

int XLALSimInspiralTestingGRCorrections(COMPLEX16FrequencySeries *htilde,       /**< input htilde, will be modified in place */
                                        const REAL8 m1_SI,
                                        const REAL8 m2_SI,
                                        const REAL8 chi1z,
                                        const REAL8 chi2z,
                                        const REAL8 f_low,
                                        const REAL8 f_ref,
                                        const REAL8 lambda1,
                                        const REAL8 lambda2,
                                        const REAL8 f_window_div_f_Peak,     /** Frequency at which to attach non-GR and GR waveforms, inputted as a fraction of f_Peak (should be between 0 and 1) */
                                        const REAL8 NCyclesStep,                /** Number of GW cycles over which to taper the non-GR phase correction */
                                        const LALSimInspiralTestGRParam *pnCorrections    /**< input linked list of testing gr parameters */
					)
{
  /* check if we have a NULL pnCorrections pointer. If yes, just return */
  if (pnCorrections==NULL) return 0;
  /* external: SI; internal: solar masses */
  const REAL8 f0 = htilde->f0;
  const REAL8 deltaF = htilde->deltaF;
  const REAL8 m1 = m1_SI / LAL_MSUN_SI;
  const REAL8 m2 = m2_SI / LAL_MSUN_SI;
  const REAL8 m = m1 + m2;
  const REAL8 m_sec = m * LAL_MTSUN_SI;  /* total mass in seconds */
  const REAL8 eta = m1 * m2 / (m * m);
  REAL8 fPeak;
  
  /* Compute the frequency where the amplitude of 2,2 mode peaks differently for BBH, BNS, and NSBH:
     For BBH, use (unpublished) fit to NR used in construction of SEOBNRv4
     For BNS, use (unpublished) fit to NR used in construction of NRTidal models (documented in LALSimNRTunedTides.c)
     For NSBH, use same fit as BBH
  */
  
  if(lambda1 == 0.0 && lambda2 == 0.0){
    fPeak = GetNRSpinPeakOmegaV4(2, 2, eta, 0.5*(chi1z + chi2z) + 0.5*(chi1z - chi2z)*(m1 - m2)/(m1 + m2)/(1. - 2.*eta)) / (LAL_PI * m_sec);
    //Spin combination defined in Eq. (2.8) of Bohe et al. (2017) (SEOBNRv4 paper)
  }
  else if(lambda1 != 0.0 && lambda2 != 0.0)
  {
    const REAL8 q = fmax(m1 / m2, m2 / m1);
    const double kappa2T = XLALSimNRTunedTidesComputeKappa2T(m1_SI, m2_SI, lambda1, lambda2);
    fPeak = XLALSimNRTunedTidesMergerFrequency(m, kappa2T, q);
  }
  else{
    fPeak = GetNRSpinPeakOmegaV4(2, 2, eta, 0.5*(chi1z + chi2z) + 0.5*(chi1z - chi2z)*(m1 - m2)/(m1 + m2)/(1. - 2.*eta)) / (LAL_PI * m_sec);
  }
  
  INT4 i;
  INT4 n = (INT4) htilde->data->length;
  
  /* Indices of f0, f_low, frequency at which non-GR modifications end, and fPeak */
  INT4 iStart, iRef, iEnd, iPeak;
  iStart = (UINT4) ceil((f_low-f0) / deltaF);
  iRef   = (UINT4) ceil((f_ref-f0) / deltaF);
  iEnd  = (UINT4) fmin(ceil((f_window_div_f_Peak * fPeak - f0) / deltaF),n-1);
  iPeak  = (UINT4) fmin(ceil((fPeak - f0) / deltaF),n-1);
  
  /* Sequence of frequencies where corrections to the model need to be evaluated
   * Fill with non-zero vals from f0 to fEnd
   */
  REAL8Sequence *freqs =NULL;
  freqs = XLALCreateREAL8Sequence(n);
  
  for (i = 0; i < n; i++)
    {
      freqs->data[i] = f0 + i * deltaF;
    }
  
  PNPhasingSeries pfa;
  const REAL8 qm_def1 = 1.;
  const REAL8 qm_def2 = 1.;
  XLALSimInspiralPNCorrections(&pfa, m1, m2, chi1z, chi2z, chi1z*chi1z, chi2z*chi2z, chi1z*chi2z, qm_def1, qm_def2, pnCorrections);
  XLALSimInspiralPhaseCorrectionsPhasing(htilde,freqs,iStart,iRef,iEnd,iPeak,pfa,m_sec, eta, NCyclesStep);
  XLALDestroyREAL8Sequence(freqs);
  return 0;
}

/* Computes the PN coefficients for the non-GR phase correction and stores in pfa.
 * The TestGRParam values represent fractional deviation from the corresponing PN
 * coefficients in TaylorF2 expression for the phase, except where the
 * coefficeint vanishes in GR, in which case the TestGRParam values indicates the
 * total value of that PN coefficient. The code closely follows XLALSimInspiralPNPhasing_F2()
 * in LALSimInspiralPNCoefficients.c, but returns only the PN coefficients for the
 * correction instead of the TaylorF2 value + correction.
 */

void XLALSimInspiralPNCorrections(PNPhasingSeries *pfa,
                                             const REAL8 m1, /**< Mass of body 1, in Msol */
                                             const REAL8 m2, /**< Mass of body 2, in Msol */
                                             const REAL8 chi1L, /**< Component of dimensionless spin 1 along Lhat */
                                             const REAL8 chi2L, /**< Component of dimensionless spin 2 along Lhat */
                                             const REAL8 chi1sq,/**< Magnitude of dimensionless spin 1 */
                                             const REAL8 chi2sq, /**< Magnitude of dimensionless spin 2 */
                                             const REAL8 chi1dotchi2, /**< Dot product of dimensionles spin 1 and spin 2 */
                                             const REAL8 qm_def1, /**< Quadrupole deformation parameter of body 1 (dimensionless) */
                                             const REAL8 qm_def2, /**< Quadrupole deformation parameter of body 2 (dimensionless) */
                                             const LALSimInspiralTestGRParam *pnCorrections)
{
    const REAL8 mtot = m1 + m2;
    const REAL8 eta = m1*m2/mtot/mtot;
    const REAL8 pfaN = 3.L/(128.L * eta);
    const REAL8 d = (m1 - m2) / (m1 + m2);
    const REAL8 m1M = m1/mtot;
    const REAL8 m2M = m2/mtot;
    const REAL8 SL = m1M*m1M*chi1L + m2M*m2M*chi2L;
    const REAL8 dSigmaL = d*(m2M*chi2L - m1M*chi1L);
    REAL8 pn_sigma = eta * (721.L/48.L*chi1L*chi2L - 247.L/48.L*chi1dotchi2);
    pn_sigma += (720.L*qm_def1 - 1.L)/96.0L * m1M * m1M * chi1L * chi1L;
    pn_sigma += (720.L*qm_def2 - 1.L)/96.0L * m2M * m2M * chi2L * chi2L;
    pn_sigma -= (240.L*qm_def1 - 7.L)/96.0L * m1M * m1M * chi1sq;
    pn_sigma -= (240.L*qm_def2 - 7.L)/96.0L * m2M * m2M * chi2sq;
    REAL8 pn_ss3 =  (326.75L/1.12L + 557.5L/1.8L*eta)*eta*chi1L*chi2L;
    pn_ss3 += ((4703.5L/8.4L+2935.L/6.L*m1M-120.L*m1M*m1M)*qm_def1 + (-4108.25L/6.72L-108.5L/1.2L*m1M+125.5L/3.6L*m1M*m1M)) *m1M*m1M * chi1sq;
    pn_ss3 += ((4703.5L/8.4L+2935.L/6.L*m2M-120.L*m2M*m2M)*qm_def2 + (-4108.25L/6.72L-108.5L/1.2L*m2M+125.5L/3.6L*m2M*m2M)) *m2M*m2M * chi2sq;
    const REAL8 pn_gamma = (554345.L/1134.L + 110.L*eta/9.L)*SL + (13915.L/84.L - 10.L*eta/3.L)*dSigmaL;

    /* initialise the PN correction  coefficients to 0 identically */
    memset(pfa, 0, sizeof(PNPhasingSeries));
    
    /* Non-spinning corrections to phasing terms - see arXiv:0907.0700, Eq. 3.18 */
    if (XLALSimInspiralTestGRParamExists(pnCorrections,"dchi0")) pfa->v[0] = XLALSimInspiralGetTestGRParam(pnCorrections,"dchi0");
    if (XLALSimInspiralTestGRParamExists(pnCorrections,"dchi1")) pfa->v[1] = XLALSimInspiralGetTestGRParam(pnCorrections,"dchi1");
    if (XLALSimInspiralTestGRParamExists(pnCorrections,"dchiMinus1")) pfa->vneg[1] = XLALSimInspiralGetTestGRParam(pnCorrections,"dchiMinus1"); 
    if (XLALSimInspiralTestGRParamExists(pnCorrections,"dchiMinus2")) pfa->vneg[2] = XLALSimInspiralGetTestGRParam(pnCorrections,"dchiMinus2"); 
    
    if (XLALSimInspiralTestGRParamExists(pnCorrections,"dchi2"))
    {
        pfa->v[2] = 5.L*(743.L/84.L + 11.L * eta)/9.L;
        pfa->v[2] *= XLALSimInspiralGetTestGRParam(pnCorrections,"dchi2");
    }
    if (XLALSimInspiralTestGRParamExists(pnCorrections,"dchi3"))
    {
        pfa->v[3] = -16.L*LAL_PI;
        pfa->v[3] += 188.L*SL/3.L + 25.L*dSigmaL;
        pfa->v[3] *= XLALSimInspiralGetTestGRParam(pnCorrections,"dchi3");
    }
    if (XLALSimInspiralTestGRParamExists(pnCorrections,"dchi4"))
    {
        pfa->v[4] = 5.L*(3058.673L/7.056L + 5429.L/7.L * eta
                        + 617.L * eta*eta)/72.L;
        pfa->v[4] += -10.L * pn_sigma;
        pfa->v[4] *= XLALSimInspiralGetTestGRParam(pnCorrections,"dchi4");
    }
    if (XLALSimInspiralTestGRParamExists(pnCorrections,"dchi5"))
    {
        pfa->v[5] = 5.L/9.L * (7729.L/84.L - 13.L * eta) * LAL_PI;
        pfa->v[5] += -1.L * pn_gamma;
        pfa->v[5] *= XLALSimInspiralGetTestGRParam(pnCorrections,"dchi5");
    }
    if (XLALSimInspiralTestGRParamExists(pnCorrections,"dchi5l"))
    {
        pfa->vlogv[5] = 5.L/3.L * (7729.L/84.L - 13.L * eta) * LAL_PI;
        pfa->vlogv[5] += -3.L * pn_gamma;
        pfa->vlogv[5] *= XLALSimInspiralGetTestGRParam(pnCorrections,"dchi5l");
    }
    if (XLALSimInspiralTestGRParamExists(pnCorrections,"dchi6"))
    {
        pfa->v[6] = (11583.231236531L/4.694215680L
                     - 640.L/3.L * LAL_PI * LAL_PI - 6848.L/21.L*LAL_GAMMA)
        + eta * (-15737.765635L/3.048192L
                 + 2255./12. * LAL_PI * LAL_PI)
        + eta*eta * 76055.L/1728.L
        - eta*eta*eta * 127825.L/1296.L;
        pfa->v[6] += (-6848.L/21.L)*log(4.);
        pfa->v[6] += LAL_PI * (3760.L*SL + 1490.L*dSigmaL)/3.L + pn_ss3;
        pfa->v[6] *= XLALSimInspiralGetTestGRParam(pnCorrections,"dchi6");
    }
    if (XLALSimInspiralTestGRParamExists(pnCorrections,"dchi6l"))
    {
        pfa->vlogv[6] = -6848.L/21.L;
        pfa->vlogv[6] *= XLALSimInspiralGetTestGRParam(pnCorrections,"dchi6l");
    }
    if (XLALSimInspiralTestGRParamExists(pnCorrections,"dchi7"))
    {
        pfa->v[7] = LAL_PI * ( 77096675.L/254016.L
                              + 378515.L/1512.L * eta - 74045.L/756.L * eta*eta);
        pfa->v[7] += (-8980424995.L/762048.L + 6586595.L*eta/756.L - 305.L*eta*eta/36.L)*SL - (170978035.L/48384.L - 2876425.L*eta/672.L - 4735.L*eta*eta/144.L) * dSigmaL;
        pfa->v[7] *= XLALSimInspiralGetTestGRParam(pnCorrections,"dchi7");
    }
  
    /* At the very end, multiply everything in the series by the newtonian term */
    for(int ii = 0; ii <= PN_PHASING_SERIES_MAX_ORDER; ii++)
    {
        pfa->v[ii] *= pfaN;
        pfa->vlogv[ii] *= pfaN;
        pfa->vlogvsq[ii] *= pfaN;
        pfa->vneg[ii] *= pfaN;
    }
}

/* Compute phase correction that tapers to zero at freqs->data[iEnd] and add to
 * the phase of of the waveform htilde.
 */

int XLALSimInspiralPhaseCorrectionsPhasing(COMPLEX16FrequencySeries *htilde,       /**< input htilde, will be modified in place */
                                           const REAL8Sequence *freqs,
                                           const UINT4 iStart,
                                           const UINT4 iRef,
                                           const UINT4 iEnd,
                                           const UINT4 iPeak,
                                           PNPhasingSeries pfa,
                                           const REAL8 mtot,
                                           const REAL8 eta,
                                           const REAL8 NCyclesStep)  /** Choose number of GW cycles over which to taper the non-GR phase correction */
{
    
  UINT4 i;
  const REAL8 pfaN0 = 3.L/(128.L * eta);
  const REAL8 piM = LAL_PI * mtot;
  const REAL8 vWindow = cbrt(piM * freqs->data[iEnd]); //Center of the tapering step function in v-space
  const REAL8 width = (NCyclesStep * LAL_PI * vWindow * vWindow * vWindow * vWindow * vWindow * vWindow)/(50. * pfaN0); //Width of tapering step function

  
  /* Compute the second derivative of the phase correction and taper with step function
   * Phase correction only generated from f_low to f_max
   */
  REAL8Sequence *d2phasenonGRdf2Tapered = NULL;
  d2phasenonGRdf2Tapered = XLALCreateREAL8Sequence( (UINT4) freqs->length );
  for ( i = iStart; i < freqs->length; i++)
    {
      const REAL8 f = freqs->data[i];
      if (f>0) d2phasenonGRdf2Tapered->data[i] = PNPhaseSecondDerivative(f, pfa, mtot)/ (1. + exp( (cbrt(piM * freqs->data[i]) - vWindow) / width ) );
    }
  
  /* Interpolate and integrate (twice) the tapered second derivative of the phase correction to compute the phase correction
   * Set the value of the correction and its first derivative to zero at f_ref
   */
  REAL8Sequence *dphasenonGRdfTapered = NULL, *phasenonGRTapered = NULL;
  dphasenonGRdfTapered = XLALCreateREAL8Sequence( (UINT4) freqs->length );
  phasenonGRTapered = XLALCreateREAL8Sequence( (UINT4) freqs->length );

  gsl_spline *splineTemp = NULL;
  gsl_interp_accel *acc = NULL;
  splineTemp = gsl_spline_alloc (gsl_interp_cspline, freqs->length);
  gsl_spline_init(splineTemp, freqs->data, d2phasenonGRdf2Tapered->data, freqs->length);

  dphasenonGRdfTapered->data[iRef] = 0.0;
  for ( i = iRef; i-- > iStart;  )
    dphasenonGRdfTapered->data[i] = dphasenonGRdfTapered->data[i+1] - gsl_spline_eval_integ(splineTemp, freqs->data[i], freqs->data[i+1], acc);
  for ( i = iRef+1; i < freqs->length; i++ )
    dphasenonGRdfTapered->data[i] = dphasenonGRdfTapered->data[i-1] + gsl_spline_eval_integ(splineTemp, freqs->data[i-1], freqs->data[i], acc);
    
  gsl_spline_init(splineTemp, freqs->data, dphasenonGRdfTapered->data, freqs->length);
  phasenonGRTapered->data[iRef] = 0.0;
  for ( i = iRef; i-- > iStart;  )
    phasenonGRTapered->data[i] = phasenonGRTapered->data[i+1] - gsl_spline_eval_integ(splineTemp, freqs->data[i], freqs->data[i+1], acc);
  for ( i = iRef+1; i < freqs->length; i++ )
    phasenonGRTapered->data[i] = phasenonGRTapered->data[i-1] + gsl_spline_eval_integ(splineTemp, freqs->data[i-1], freqs->data[i], acc);


  /* Compute the first derivative of the phase correction at fPeak and subtract from
   * the phase correction to ensure the derivative vanishes at fPeak to accomplish alignment
   * of time-domain waveform.
   *
   * Then add phase correction to input waveform phase.
   */
  REAL8 PNPhaseRefDerivative = dphasenonGRdfTapered->data[iPeak];
  
  for ( i = iStart; i < freqs->length; i++ ) {
    REAL8 phasing = phasenonGRTapered->data[i] - PNPhaseRefDerivative * (freqs->data[i]-freqs->data[iRef]) ;
    htilde->data->data[i] *= cexp(phasing * 1.j);
  }
    
  gsl_spline_free(splineTemp);
  gsl_interp_accel_free(acc);
  XLALDestroyREAL8Sequence(d2phasenonGRdf2Tapered);
  XLALDestroyREAL8Sequence(dphasenonGRdfTapered);
  XLALDestroyREAL8Sequence(phasenonGRTapered);
  return 0;
}

int XLALSimInspiralTestingGRCorrectionsWithDS(COMPLEX16FrequencySeries *htilde,       /**< input htilde, will be modified in place */
                                        const REAL8 distance,
                                        const REAL8 m1_SI,
                                        const REAL8 m2_SI,
                                        const REAL8 chi1z,
                                        const REAL8 chi2z,
                                        const REAL8 f_low,
                                        const REAL8 f_ref,
					const REAL8 f_window_div_f_Peak,     /** Frequency at which to attach non-GR to GR waveforms, inputted as a fraction of f_Peak (should be between 0 and 1) */
					const REAL8 NCyclesStep,             /** Number of GW cycles over which to taper the non-GR phase correction to GR waveform */
					const REAL8 f_DS,                    /** Frequency at which to attach GR to non-GR waveforms, inputted in Hz */
					const REAL8 NCyclesDS,               /** Number of GW cycles over which to taper the GR waveform to non-GR phase corrections */                 
                                        const LALSimInspiralTestGRParam *pnCorrections    /**< input linked list of testing gr parameters */
)
{
    /* check if we have a NULL pnCorrections pointer. If yes, just return */
    if (pnCorrections==NULL) return 0;
    /* external: SI; internal: solar masses */
    const REAL8 f0 = htilde->f0;
    const REAL8 deltaF = htilde->deltaF;
    const REAL8 m1 = m1_SI / LAL_MSUN_SI;
    const REAL8 m2 = m2_SI / LAL_MSUN_SI;
    const REAL8 m = m1 + m2;
    const REAL8 m_sec = m * LAL_MTSUN_SI;  /* total mass in seconds */
    const REAL8 eta = m1 * m2 / (m * m);
//    const REAL8 piM = LAL_PI * m_sec;
//    const REAL8 vISCO = 1. / sqrt(6.);
//    const REAL8 fISCO = vISCO * vISCO * vISCO / piM;
    //const REAL8 fPeak = GetNRSpinPeakOmegaV4(2, 2, eta, 0.) / (2. * LAL_PI * m_sec);
    const REAL8 fPeak =0.0165  /  m_sec;
    
    INT4 i;
    INT4 n = (INT4) htilde->data->length;
    
    INT4 iStart, iEnd;
    /* Fill with non-zero vals from f0 to fEnd */
    iStart = (UINT4) ceil((f_low-f0) / deltaF);
    iEnd  = (UINT4) ceil(f_window_div_f_Peak * fPeak / deltaF);
//    printf("fPeak %.16e\n", fPeak);
    /* Sequence of frequencies where corrections to the model need to be evaluated */
    REAL8Sequence *freqs =NULL;
    freqs = XLALCreateREAL8Sequence(n);
    
    for (i = 0; i < n; i++)
    {
        freqs->data[i] = f0 + i * deltaF;
    }
    PNPhasingSeries pfa;
    const REAL8 qm_def1 = 1.;
    const REAL8 qm_def2 = 1.;
    XLALSimInspiralPNCorrections(&pfa, m1, m2, chi1z, chi2z, chi1z*chi1z, chi2z*chi2z, chi1z*chi2z, qm_def1, qm_def2, pnCorrections);
    XLALSimInspiralPhaseCorrectionsPhasingWithDS(htilde,distance,freqs,iStart,iEnd,pfa,m_sec, eta,f_ref, NCyclesStep, f_DS, NCyclesDS);
    XLALDestroyREAL8Sequence(freqs);
    return 0;
}

int XLALSimInspiralPhaseCorrectionsPhasingWithDS(COMPLEX16FrequencySeries *htilde,       /**< input htilde, will be modified in place */
                                           const REAL8 distance,
                                           const REAL8Sequence *freqs,
                                           const UINT4 iStart,
                                           const UINT4 iEnd,
                                           PNPhasingSeries pfa,
                                           const REAL8 mtot,
                                           const REAL8 eta,
					   UNUSED const REAL8 f_ref, /** this must be in seconds **/
					   const REAL8 NCyclesStep, /** Choose number of GW cycles over which to taper the non-GR phase correction */
					   const REAL8 f_DS,  /** Frequency in which dynamical scalarization turns on */
					   const REAL8 NCyclesDS) /** Choose number of GW cycles over which to taper the GR phase correction */
{
    UINT4 i;
    REAL8Sequence *phase = NULL, *uphase = NULL, *d1phidf1 = NULL; //, *d2phidf2 = NULL;
    phase = XLALCreateREAL8Sequence( (UINT4) freqs->length );
    uphase = XLALCreateREAL8Sequence( (UINT4) freqs->length );
    d1phidf1 = XLALCreateREAL8Sequence( (UINT4) freqs->length );
    //d2phidf2 = XLALCreateREAL8Sequence( (UINT4) freqs->length );
    memset (phase->data, 0, freqs->length * sizeof (REAL8));

    COMPLEX16 tmph;
    for ( i = 0; i < phase->length; i++ ) {
        tmph = distance*htilde->data->data[i];
        if ( cabs(tmph) < 1.e-16 ) {
            phase->data[i] = 0.;
        }
        else {
            phase->data[i] = carg(tmph);
        }
    }
    INT4 phaseCounter = 0;
    REAL8 diff = 0.;
    uphase->data[0] = phase->data[0];
    for ( i = 1; i < phase->length; i++ ) {
        diff = phase->data[i] - phase->data[i-1];
        if ( diff> 5. ) phaseCounter--;
        if ( diff<-5. ) phaseCounter++;
        uphase->data[i] = phase->data[i] + LAL_TWOPI*phaseCounter;
    }
    
    gsl_spline *splineGR = NULL, *splinenonGR = NULL, *splineTot = NULL;
    gsl_interp_accel *acc = NULL;
    splineGR = gsl_spline_alloc (gsl_interp_cspline, freqs->length);
    acc = gsl_interp_accel_alloc();
    gsl_spline_init(splineGR, freqs->data, uphase->data, freqs->length);
    for ( i = 0; i < freqs->length; i++ ) {
        d1phidf1->data[i] = gsl_spline_eval_deriv(splineGR, freqs->data[i], acc);
        //d2phidf2->data[i] = gsl_spline_eval_deriv2(splineGR, freqs->data[i], acc);
    }
//    FILE *out = fopen ("unp.dat", "w");
//    for ( i = 0; i < phase->length; i++ ) {
//        fprintf(out, "%.16e %.16e %.16e %.16e\n",freqs->data[i],phase->data[i],uphase->data[i], d2phidf2->data[i] );
//    }
//    fclose(out);

    const REAL8 pfaN0 = 3.L/(128.L * eta);
    const REAL8 piM = LAL_PI * mtot;
    const REAL8 pfa7 = pfa.v[7];
    const REAL8 pfa6 = pfa.v[6];
    const REAL8 pfl6 = pfa.vlogv[6];
    const REAL8 pfa5 = pfa.v[5];
    const REAL8 pfl5 = pfa.vlogv[5];
    const REAL8 pfa4 = pfa.v[4];
    const REAL8 pfa3 = pfa.v[3];
    const REAL8 pfa2 = pfa.v[2];
    const REAL8 pfa1 = pfa.v[1];
    const REAL8 pfaN = pfa.v[0];
    const REAL8 pfaMinus1 = pfa.vneg[1];
    const REAL8 pfaMinus2 = pfa.vneg[2];
    
    /* Compute the SPA phase at the reference point
     * N.B. f_ref == 0 means we define the reference time/phase at "coalescence"
     * when the frequency approaches infinity. In that case,
     * the integrals Eq. 3.15 of arXiv:0907.0700 vanish when evaluated at
     * f_ref == infinity. If f_ref is finite, we must compute the SPA phase
     * evaluated at f_ref, store it as ref_phasing and subtract it off.
     */
    REAL8 ref_phasing = 0.;
    

    const REAL8 vWindow = cbrt(piM*freqs->data[iEnd]); //Center of the tapering step function in v-space
    const REAL8 vDS = cbrt(piM * f_DS);

    const REAL8 width = (NCyclesStep*LAL_PI*vWindow*vWindow*vWindow*vWindow*vWindow*vWindow)/(50.* pfaN0);
    const REAL8 widthDS = (NCyclesDS*LAL_PI*vDS*vDS*vDS*vDS*vDS*vDS)/(50.* pfaN0); 
    //printf("vPeak %.16e width %.16e\n",vPeak,width);

    REAL8Sequence *phasenonGR = NULL, *d1phinonGRdf1 = NULL;  //, *d2phinonGRdf2 = NULL;
    phasenonGR = XLALCreateREAL8Sequence( (UINT4) freqs->length );
    d1phinonGRdf1 = XLALCreateREAL8Sequence( (UINT4) freqs->length );
    //d2phinonGRdf2 = XLALCreateREAL8Sequence( (UINT4) freqs->length );
    for ( i = 0; i < freqs->length; i++)
    {
        const REAL8 f = freqs->data[i];
        if (f>0)
        {
            const REAL8 v = cbrt(piM*f);
            const REAL8 logv = log(v);
            const REAL8 v2 = v * v;
            const REAL8 v3 = v * v2;
            const REAL8 v4 = v * v3;
            const REAL8 v5 = v * v4;
            const REAL8 v6 = v * v5;
            const REAL8 v7 = v * v6;
            REAL8 phasing = 0.0;
            
            phasing += pfa7 * v7;
            phasing += (pfa6 + pfl6 * logv) * v6;
            phasing += (pfa5 + pfl5 * logv) * v5;
            phasing += pfa4 * v4;
            phasing += pfa3 * v3;
            phasing += pfa2 * v2;
            phasing += pfa1 * v;
            phasing += pfaN;
            phasing += pfaMinus1 / v;
            phasing += pfaMinus2 /v2;
            phasing /= v5;
            
            phasenonGR->data[i] = -phasing;
        }
    }
    
    splinenonGR = gsl_spline_alloc (gsl_interp_cspline, freqs->length);
    gsl_spline_init(splinenonGR, freqs->data, phasenonGR->data, freqs->length);
    for ( i = 0; i < freqs->length; i++ )
        d1phinonGRdf1->data[i] = gsl_spline_eval_deriv(splinenonGR, freqs->data[i], acc);
    //d2phinonGRdf2->data[i] = gsl_spline_eval_deriv2(splinenonGR, freqs->data[i], acc);

    REAL8Sequence *d1phaseTotdf1 = NULL;  //, *d2phaseTotdf2 = NULL;
    d1phaseTotdf1 = XLALCreateREAL8Sequence( (UINT4) freqs->length );
    //d2phaseTotdf2 = XLALCreateREAL8Sequence( (UINT4) freqs->length );
    for ( i = 0; i < freqs->length; i++ )
      d1phaseTotdf1->data[i] = d1phidf1->data[i] + d1phinonGRdf1->data[i] / (1. + exp( (cbrt(piM*freqs->data[i]) - vWindow) / width) ) * ( 1. - 1. / (1. + exp( (cbrt(piM*freqs->data[i]) - vDS) / widthDS) ) );
      //d2phaseTotdf2->data[i] = d2phidf2->data[i] + d2phinonGRdf2->data[i] / (1. + exp( (cbrt(piM*freqs->data[i]) - vWindow) / width ) );
    
      //REAL8 dphaseTotGRIni = gsl_spline_eval_deriv(splineGR, freqs->data[iStart], acc) + gsl_spline_eval_deriv(splinenonGR, freqs->data[iStart], acc);   //Removed by Noah
    REAL8 phaseTotIni = uphase->data[iStart] + phasenonGR->data[iStart];
    
    REAL8Sequence *dphaseTotdf = NULL, *phaseTot = NULL;
    dphaseTotdf = XLALCreateREAL8Sequence( (UINT4) freqs->length );
    phaseTot = XLALCreateREAL8Sequence( (UINT4) freqs->length );
  
    splineTot = gsl_spline_alloc (gsl_interp_cspline, freqs->length);

    //gsl_spline_init(splineTot, freqs->data, d2phaseTotdf2->data, freqs->length);  //Removed by Noah
    //dphaseTotdf->data[iStart] = dphaseTotGRIni;
    //for ( i = iStart; i-- > 0;  )
    //    dphaseTotdf->data[i] = dphaseTotdf->data[i+1] - gsl_spline_eval_integ(splineTot, freqs->data[i], freqs->data[i+1], acc);
    //for ( i = iStart+1; i < freqs->length; i++ )
    //    dphaseTotdf->data[i] = dphaseTotdf->data[i-1] + gsl_spline_eval_integ(splineTot, freqs->data[i-1], freqs->data[i], acc);
    
    gsl_spline_init(splineTot, freqs->data, d1phaseTotdf1->data, freqs->length);
    phaseTot->data[iStart] = phaseTotIni;
    for ( i = iStart; i-- > 0;  )
        phaseTot->data[i] = phaseTot->data[i+1] - gsl_spline_eval_integ(splineTot, freqs->data[i], freqs->data[i+1], acc);
    for ( i = iStart+1; i < freqs->length; i++ )
        phaseTot->data[i] = phaseTot->data[i-1] + gsl_spline_eval_integ(splineTot, freqs->data[i-1], freqs->data[i], acc);

//    FILE *out = fopen ("pin.dat", "w");
//    for ( i = 0; i < phase->length; i++ ) {
//        fprintf(out, "%.16e %.16e %.16e %.16e %.16e %.16e\n",freqs->data[i],uphase->data[i],phaseTot->data[i], phasenonGR->data[i], d2phaseTotdf2->data[i], dphaseTotdf->data[i] );
//    }
//    fclose(out);

    gsl_spline_init(splineTot, freqs->data, phaseTot->data, freqs->length);
    ref_phasing = gsl_spline_eval(splineTot, f_ref, acc);
    for ( i = 0; i < freqs->length; i++ ) {
        REAL8 phasing = 0.;
        phasing = phaseTot->data[i] - ref_phasing;
        htilde->data->data[i] = cabs(distance*htilde->data->data[i]) * (cos(phasing)+ sin(phasing) * 1.0j);
        htilde->data->data[i] /= distance;
    }
    
    gsl_spline_free(splineGR);
    gsl_spline_free(splinenonGR);
    gsl_spline_free(splineTot);
    gsl_interp_accel_free(acc);
    XLALDestroyREAL8Sequence(uphase);
    XLALDestroyREAL8Sequence(phase);
    XLALDestroyREAL8Sequence(d1phidf1);
    //XLALDestroyREAL8Sequence(d2phidf2);
    XLALDestroyREAL8Sequence(d1phinonGRdf1);
    //XLALDestroyREAL8Sequence(d2phinonGRdf2);
    XLALDestroyREAL8Sequence(phasenonGR);
    XLALDestroyREAL8Sequence(d1phaseTotdf1);
    //XLALDestroyREAL8Sequence(d2phaseTotdf2);
    XLALDestroyREAL8Sequence(dphaseTotdf);
    XLALDestroyREAL8Sequence(phaseTot);
    return 0;
}

REAL8 PNPhase(REAL8 f,    /* frequency in Hz */
              PNPhasingSeries pfa,
              const REAL8 mtot)   /* total mass in seconds */
{
  
    const REAL8 piM = LAL_PI * mtot;
    const REAL8 pfa7 = pfa.v[7];
    const REAL8 pfa6 = pfa.v[6];
    const REAL8 pfl6 = pfa.vlogv[6];
    const REAL8 pfa5 = pfa.v[5];
    const REAL8 pfl5 = pfa.vlogv[5];
    const REAL8 pfa4 = pfa.v[4];
    const REAL8 pfa3 = pfa.v[3];
    const REAL8 pfa2 = pfa.v[2];
    const REAL8 pfa1 = pfa.v[1];
    const REAL8 pfaN = pfa.v[0];
    const REAL8 pfaMinus1 = pfa.vneg[1];
    const REAL8 pfaMinus2 = pfa.vneg[2];
    
    const REAL8 v = cbrt(piM*f);
    const REAL8 logv = log(v);
    const REAL8 v2 = v * v;
    const REAL8 v3 = v * v2;
    const REAL8 v4 = v * v3;
    const REAL8 v5 = v * v4;
    const REAL8 v6 = v * v5;
    const REAL8 v7 = v * v6;
    REAL8 phasing = 0.0;

    phasing += pfa7 * v7;
    phasing += (pfa6 + pfl6 * logv) * v6;
    phasing += (pfa5 + pfl5 * logv) * v5;
    phasing += pfa4 * v4;
    phasing += pfa3 * v3;
    phasing += pfa2 * v2;
    phasing += pfa1 * v;
    phasing += pfaN;
    phasing += pfaMinus1 / v;
    phasing += pfaMinus2 /v2;
    phasing /= v5;

    return -phasing;
}

/* Derivative of phase with respect to f
 * d[f^((-5 + n)/3)]= (-5 + n) / 3 * f^((-8 + n)/3)
 * d[f^((-5 + n)/3) * log(f)] = f^((-8 + n)/3) + (-5 + n)/3 * f^((-8 + n)/3) * log(f) */

REAL8 PNPhaseDerivative(REAL8 f,    /* frequency in Hz */
                        PNPhasingSeries pfa,
                        const REAL8 mtot)   /* total mass in seconds */
{
  
    const REAL8 piM = LAL_PI * mtot;
    const REAL8 pfa7 = pfa.v[7];
    const REAL8 pfa6 = pfa.v[6];
    const REAL8 pfl6 = pfa.vlogv[6];
    //const REAL8 pfa5 = pfa.v[5];
    const REAL8 pfl5 = pfa.vlogv[5];
    const REAL8 pfa4 = pfa.v[4];
    const REAL8 pfa3 = pfa.v[3];
    const REAL8 pfa2 = pfa.v[2];
    const REAL8 pfa1 = pfa.v[1];
    const REAL8 pfaN = pfa.v[0];
    const REAL8 pfaMinus1 = pfa.vneg[1];
    const REAL8 pfaMinus2 = pfa.vneg[2];
    
    const REAL8 v = cbrt(piM*f);
    const REAL8 logv = log(v);
    const REAL8 v2 = v * v;
    const REAL8 v3 = v * v2;
    const REAL8 v4 = v * v3;
    const REAL8 v5 = v * v4;
    //const REAL8 v6 = v * v5;
    //const REAL8 v7 = v * v6;
    REAL8 phasing = 0.0;

    phasing += 2./3. * pfa7 * v4;
    phasing += (1./3. * pfa6 + pfl6 + 1./3. * pfl6 * logv) * v3;
    phasing += pfl5 * v2;
    phasing += -1./3. * pfa4 * v;
    phasing += -2./3. * pfa3;
    phasing += -1. * pfa2 / v;
    phasing += -4./3. * pfa1 / v2;
    phasing += -5./3. * pfaN / v3;
    phasing += -2. * pfaMinus1 / v4;
    phasing += -7./3. * pfaMinus2 /v5;
    phasing /= v5;
    phasing *= piM;

    return -phasing;
}

/* Derivative of phase with respect to f
 * d2[f^((-5 + n)/3)]= (-5 + n) * (-8 + n) / 9 * f^((-11 + n)/3)
 * d2[f^((-5 + n)/3) * log(f)] = (-13 + 2 * n)/3 * f^((-11 + n)/3) + (-5 + n) * (-8 + n) / 9 * f^((-11 + n)/3) * log(f) */

REAL8 PNPhaseSecondDerivative(REAL8 f,    /* frequency in Hz */
                              PNPhasingSeries pfa,
                              const REAL8 mtot)   /* total mass in seconds */
{
  
    const REAL8 piM = LAL_PI * mtot;
    const REAL8 pfa7 = pfa.v[7];
    const REAL8 pfa6 = pfa.v[6];
    const REAL8 pfl6 = pfa.vlogv[6];
    //const REAL8 pfa5 = pfa.v[5];
    const REAL8 pfl5 = pfa.vlogv[5];
    const REAL8 pfa4 = pfa.v[4];
    const REAL8 pfa3 = pfa.v[3];
    const REAL8 pfa2 = pfa.v[2];
    const REAL8 pfa1 = pfa.v[1];
    const REAL8 pfaN = pfa.v[0];
    const REAL8 pfaMinus1 = pfa.vneg[1];
    const REAL8 pfaMinus2 = pfa.vneg[2];
    
    const REAL8 v = cbrt(piM*f);
    const REAL8 logv = log(v);
    const REAL8 v2 = v * v;
    const REAL8 v3 = v * v2;
    const REAL8 v4 = v * v3;
    const REAL8 v5 = v * v4;
    const REAL8 v6 = v * v5;
    const REAL8 v7 = v * v6;
    const REAL8 v8 = v * v7;
    REAL8 phasing = 0.0;

    phasing += -2./9. *  pfa7 * v;
    phasing += (-2./9. *  pfa6 - 1./3. * pfl6 - 2./9. * pfl6 * logv);
    phasing += -1. * pfl5 / v;
    phasing += 4./9. * pfa4 / v2 ;
    phasing += 10./9. * pfa3 / v3;
    phasing += 2. * pfa2 / v4;
    phasing += 28./9. * pfa1 / v5;
    phasing += 40./9. * pfaN / v6;
    phasing += 6. * pfaMinus1 / v7;
    phasing += 70./9. * pfaMinus2 /v8;
    phasing /= v5;
    phasing *= piM;
    phasing *= piM;
    return -phasing;
}
