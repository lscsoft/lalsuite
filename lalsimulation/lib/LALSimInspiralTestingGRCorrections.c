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
 * Peak (angular) frequency of 2,2 mode of GW predicted by fitting NR results. Taken from
 * Eq. A8 of Bohe et al (2017) (arxiv.org/abs/1611.03703). Used in building SEOBNRv4.
 */
static inline REAL8 Get22PeakAngFreq(REAL8 UNUSED eta,
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
                                        const UINT4 l, /**< degree l of the GW mode */
                                        const UINT4 m, /**< order m of the GW mode */
                                        const REAL8 m1_SI, /**< mass of companion 1 (kg) */
                                        const REAL8 m2_SI, /**< mass of companion 2 (kg) */
                                        const REAL8 chi1z, /**< z-component of the dimensionless spin of object 1 */
                                        const REAL8 chi2z, /**< z-component of the dimensionless spin of object 2 */
                                        const REAL8 f_low, /**< lower GW frequency bound (Hz) */
                                        const REAL8 f_ref, /**< reference GW frequency (Hz) */
                                        const REAL8 f_window_div_f_Peak,     /**< Frequency at which to attach non-GR and GR waveforms, inputted as a fraction of f_Peak (should be between 0 and 1) */
                                        const REAL8 NCyclesStep,                /**< Number of GW cycles over which to taper the non-GR phase correction */
                                        LALDict *LALpars    /**< dictionary of additional parameters */
					)
{
  /* check if we have a NULL pnCorrections pointer. If yes, just return */
  if (LALpars==NULL) return 0;
  /* external: SI; internal: solar masses */
  const REAL8 f0 = htilde->f0;
  const REAL8 deltaF = htilde->deltaF;
  const REAL8 m1 = m1_SI / LAL_MSUN_SI;
  const REAL8 m2 = m2_SI / LAL_MSUN_SI;
  const REAL8 mt = m1 + m2;
  const REAL8 m_sec = mt * LAL_MTSUN_SI;  /* total mass in seconds */
  const REAL8 eta = m1 * m2 / (mt * mt);
/* Check that fRef is sane */
  if( f_ref < f_low )
  {
        XLALPrintError("XLAL Error - %s: fRef is smaller than the starting frequency of the waveform, f_low. Please pass in the starting GW frequency instead.\n", __func__);
        XLAL_ERROR(XLAL_EINVAL);
  }
  
  REAL8 lambda1 = XLALSimInspiralWaveformParamsLookupTidalLambda1(LALpars);
  REAL8 lambda2 = XLALSimInspiralWaveformParamsLookupTidalLambda2(LALpars);

  REAL8 f22Peak;
  
  /* Compute the frequency where the amplitude of 2,2 mode peaks differently for BBH, BNS, and NSBH:
     For BBH, use fit to NR used in construction of SEOBNRv4 (arxiv:1611.03703)
     For BNS, use fit to NR used in construction of NRTidal models (documented in LALSimNRTunedTides.c, arXiv:1706.02969)
     For NSBH, use same fit as BBH
  */
  
  if(lambda1 == 0.0 && lambda2 == 0.0){
    f22Peak = Get22PeakAngFreq(eta, 0.5*(chi1z + chi2z) + 0.5*(chi1z - chi2z)*(m1 - m2)/(m1 + m2)/(1. - 2.*eta)) / (2.*LAL_PI * m_sec);
    //Spin combination defined in Eq. (2.8) of Bohe et al. (2017) (SEOBNRv4 paper)
  }
  else if(lambda1 != 0.0 && lambda2 != 0.0)
  {
    const REAL8 q = fmax(m1 / m2, m2 / m1);
    const double kappa2T = XLALSimNRTunedTidesComputeKappa2T(m1_SI, m2_SI, lambda1, lambda2);
    f22Peak = XLALSimNRTunedTidesMergerFrequency(mt, kappa2T, q) + l*0;
  }
  else{
    f22Peak = Get22PeakAngFreq(eta, 0.5*(chi1z + chi2z) + 0.5*(chi1z - chi2z)*(m1 - m2)/(m1 + m2)/(1. - 2.*eta)) / (2.*LAL_PI * m_sec);
  }
  
  INT4 i;
  INT4 n = (INT4) htilde->data->length;
  
  /* Indices of f0, f_low, frequency at which non-GR modifications end, and fPeak */
  INT4 iStart, iStartFinal, iRef,  iPeak;
  iStart = (UINT4) ceil((f_low/2.-f0) / deltaF);
  iStartFinal = (UINT4) ceil((f_low-f0) / deltaF);
  iRef   = (UINT4) ceil((f_ref*m/2.-f0) / deltaF);
  iPeak  = (UINT4) fmin(ceil((f22Peak - f0) / deltaF),n-1);
  
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
  
  XLALSimInspiralPNCorrections(&pfa, m1, m2, chi1z, chi2z, chi1z*chi1z, chi2z*chi2z, chi1z*chi2z, qm_def1, qm_def2, LALpars);
  XLALSimInspiralPhaseCorrectionsPhasing(htilde,freqs,m,iStart,iRef,iPeak,pfa,m_sec,eta,f_window_div_f_Peak,iStartFinal,NCyclesStep);
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

void XLALSimInspiralPNCorrections(PNPhasingSeries *pfa, /**< PN phasing coefficients */
                                             const REAL8 m1, /**< Mass of body 1, in Msol */
                                             const REAL8 m2, /**< Mass of body 2, in Msol */
                                             const REAL8 chi1L, /**< Component of dimensionless spin 1 along Lhat */
                                             const REAL8 chi2L, /**< Component of dimensionless spin 2 along Lhat */
                                             const REAL8 chi1sq,/**< Magnitude of dimensionless spin 1 */
                                             const REAL8 chi2sq, /**< Magnitude of dimensionless spin 2 */
                                             const REAL8 chi1dotchi2, /**< Dot product of dimensionles spin 1 and spin 2 */
                                             const REAL8 qm_def1, /**< Quadrupole deformation parameter of body 1 (dimensionless) */
                                             const REAL8 qm_def2, /**< Quadrupole deformation parameter of body 2 (dimensionless) */
                                             LALDict *LALpars /**< dictionary of additional parameters */
					     )
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
    /* Introducing here the SIQM test PN cofficients..*/
    /* See eq. (0.5a) of https://arxiv.org/pdf/1701.06318.pdf */
    const REAL8 chiS = 0.5L * (chi1L + chi2L);
    const REAL8 chiA = 0.5L * (chi1L - chi2L);
    const REAL8 siqm_2pn_kappaS = 50.L * (2.L*eta - 1.L) * (chiS*chiS + chiA*chiA) -100.L * d * chiA * chiS;   
    const REAL8 siqm_2pn_kappaA = -50.L * d * (chiS*chiS + chiA*chiA) + 100.L * (2.L*eta - 1.L) * chiA * chiS;
    const REAL8 siqm_3pn_kappaS = (chiS*chiS + chiA*chiA) * (26015.L/28.L - 44255.L/21.L * eta - 240.L * eta*eta) + chiA * chiS * d * (26015.L/14.L - 1495.L/3. * eta); 
    const REAL8 siqm_3pn_kappaA = (26015.L/28.L - 1495.L/6.L * eta) * d * (chiS*chiS + chiA*chiA) + (26015.L/14.L  - 88510.L/21.L *eta - 480.L *eta*eta) * chiS * chiA;

    /* initialise the PN correction  coefficients to 0 identically */
    memset(pfa, 0, sizeof(PNPhasingSeries));
    
    /* Non-spinning corrections to phasing terms - see arXiv:0907.0700, Eq. 3.18 */
    if (XLALDictContains(LALpars,"dchiMinus2")) pfa->vneg[2] = XLALSimInspiralWaveformParamsLookupNonGRDChiMinus2(LALpars);
    if (XLALDictContains(LALpars,"dchiMinus1")) pfa->vneg[1] = XLALSimInspiralWaveformParamsLookupNonGRDChiMinus1(LALpars);

    if (XLALDictContains(LALpars,"dchi0")) pfa->v[0] = XLALSimInspiralWaveformParamsLookupNonGRDChi0(LALpars);
    if (XLALDictContains(LALpars,"dchi1")) pfa->v[1] = XLALSimInspiralWaveformParamsLookupNonGRDChi1(LALpars);

    if (XLALDictContains(LALpars,"dchi2"))
    {
        pfa->v[2] = 5.L*(743.L/84.L + 11.L * eta)/9.L;
        pfa->v[2] *= XLALSimInspiralWaveformParamsLookupNonGRDChi2(LALpars);
    }

    if (XLALDictContains(LALpars,"dchi3") || XLALDictContains(LALpars,"dchi3NS") || XLALDictContains(LALpars,"dchi3S"))
    {
        pfa->v[3] = -16.L*LAL_PI;
        pfa->v[3] += 188.L*SL/3.L + 25.L*dSigmaL;
        pfa->v[3] *= XLALSimInspiralWaveformParamsLookupNonGRDChi3(LALpars);

        /* Splitting the deviation parameter into spin and no-spin parts and then adding them here */
        pfa->v[3] += -16.L*LAL_PI * XLALSimInspiralWaveformParamsLookupNonGRDChi3NS(LALpars);
        pfa->v[3] += (188.L*SL/3.L + 25.L*dSigmaL) * XLALSimInspiralWaveformParamsLookupNonGRDChi3S(LALpars);
    }

    if (XLALDictContains(LALpars,"dchi4") || XLALDictContains(LALpars,"dchi4NS") || XLALDictContains(LALpars,"dchi4S") || XLALDictContains(LALpars,"dchikappaS") || XLALDictContains(LALpars,"dchikappaA"))
    {
        pfa->v[4] = 5.L*(3058.673L/7.056L + 5429.L/7.L * eta + 617.L * eta*eta)/72.L;
        pfa->v[4] += -10.L * pn_sigma;
        pfa->v[4] *= XLALSimInspiralWaveformParamsLookupNonGRDChi4(LALpars);

        /* Splitting the deviation parameter into spin and no-spin parts and then adding them here */
        pfa->v[4] += (5.L*(3058.673L/7.056L + 5429.L/7.L * eta + 617.L * eta*eta)/72.L) * XLALSimInspiralWaveformParamsLookupNonGRDChi4NS(LALpars);
        pfa->v[4] += -10.L*pn_sigma * XLALSimInspiralWaveformParamsLookupNonGRDChi4S(LALpars);

        /* Adding 2PN SIQM test here..*/
        pfa->v[4] += siqm_2pn_kappaS * XLALSimInspiralWaveformParamsLookupNonGRDChikappaS(LALpars);
        pfa->v[4] += siqm_2pn_kappaA * XLALSimInspiralWaveformParamsLookupNonGRDChikappaA(LALpars);
    }

    if (XLALDictContains(LALpars,"dchi5") || XLALDictContains(LALpars,"dchi5NS") || XLALDictContains(LALpars,"dchi5S"))
    {
        pfa->v[5] = 5.L/9.L * (7729.L/84.L - 13.L * eta) * LAL_PI;
        pfa->v[5] += -1.L * pn_gamma;
        pfa->v[5] *= XLALSimInspiralWaveformParamsLookupNonGRDChi5(LALpars);

        /* Splitting the deviation parameter into spin and no-spin parts and then adding them here */
        pfa->v[5] += (5.L/9.L * (7729.L/84.L - 13.L * eta) * LAL_PI) * XLALSimInspiralWaveformParamsLookupNonGRDChi5NS(LALpars);
        pfa->v[5] += -1.L*pn_gamma * XLALSimInspiralWaveformParamsLookupNonGRDChi5S(LALpars);
    }

    if (XLALDictContains(LALpars,"dchi5l") || XLALDictContains(LALpars,"dchi5lNS") || XLALDictContains(LALpars,"dchi5lS"))
    {
        pfa->vlogv[5] = 5.L/3.L * (7729.L/84.L - 13.L * eta) * LAL_PI;
        pfa->vlogv[5] += -3.L * pn_gamma;
        pfa->vlogv[5] *= XLALSimInspiralWaveformParamsLookupNonGRDChi5L(LALpars);

        /* Splitting the deviation parameter into spin and no-spin parts and then adding them here */
        pfa->vlogv[5] += (5.L/3.L * (7729.L/84.L - 13.L * eta) * LAL_PI) * XLALSimInspiralWaveformParamsLookupNonGRDChi5LNS(LALpars);
        pfa->vlogv[5] += -3.L*pn_gamma * XLALSimInspiralWaveformParamsLookupNonGRDChi5LS(LALpars);
    }

    if (XLALDictContains(LALpars,"dchi6") || XLALDictContains(LALpars,"dchi6NS") || XLALDictContains(LALpars,"dchi6S") || XLALDictContains(LALpars,"dchikappaS") || XLALDictContains(LALpars,"dchikappaA"))
    {
        pfa->v[6] = (11583.231236531L/4.694215680L - 640.L/3.L * LAL_PI * LAL_PI - 6848.L/21.L*LAL_GAMMA) + eta * (-15737.765635L/3.048192L + 2255./12. * LAL_PI * LAL_PI) + eta*eta * 76055.L/1728.L - eta*eta*eta * 127825.L/1296.L;
        pfa->v[6] += (-6848.L/21.L)*log(4.);
        pfa->v[6] += LAL_PI * (3760.L*SL + 1490.L*dSigmaL)/3.L + pn_ss3;
        pfa->v[6] *= XLALSimInspiralWaveformParamsLookupNonGRDChi6(LALpars);

        /* Splitting the deviation parameter into spin and no-spin parts and then adding them here */
        pfa->v[6] += ((11583.231236531L/4.694215680L - 640.L/3.L * LAL_PI * LAL_PI - 6848.L/21.L*LAL_GAMMA) + eta * (-15737.765635L/3.048192L + 2255./12. * LAL_PI * LAL_PI) + eta*eta * 76055.L/1728.L - eta*eta*eta * 127825.L/1296.L + (-6848.L/21.L)*log(4.)) * XLALSimInspiralWaveformParamsLookupNonGRDChi6NS(LALpars);
        pfa->v[6] += (LAL_PI * (3760.L*SL + 1490.L*dSigmaL)/3.L + pn_ss3) * XLALSimInspiralWaveformParamsLookupNonGRDChi6S(LALpars);

        /* Adding 3PN SIQM test here..*/
        pfa->v[6] += siqm_3pn_kappaS * XLALSimInspiralWaveformParamsLookupNonGRDChikappaS(LALpars);
        pfa->v[6] += siqm_3pn_kappaA * XLALSimInspiralWaveformParamsLookupNonGRDChikappaA(LALpars);
    }

    if (XLALDictContains(LALpars,"dchi6l"))
    {
        pfa->vlogv[6] = -6848.L/21.L;
        pfa->vlogv[6] *= XLALSimInspiralWaveformParamsLookupNonGRDChi6L(LALpars);
    }

    if (XLALDictContains(LALpars,"dchi7")|| XLALDictContains(LALpars,"dchi7NS") || XLALDictContains(LALpars,"dchi7S")) 
    {
        pfa->v[7] = LAL_PI * ( 77096675.L/254016.L + 378515.L/1512.L * eta - 74045.L/756.L * eta*eta);
        pfa->v[7] += (-8980424995.L/762048.L + 6586595.L*eta/756.L - 305.L*eta*eta/36.L)*SL - (170978035.L/48384.L - 2876425.L*eta/672.L - 4735.L*eta*eta/144.L) * dSigmaL;
        pfa->v[7] *= XLALSimInspiralWaveformParamsLookupNonGRDChi7(LALpars);

        /* Splitting the deviation parameter into spin and no-spin parts and then adding them here */
        pfa->v[7] += (LAL_PI * ( 77096675.L/254016.L + 378515.L/1512.L * eta - 74045.L/756.L * eta*eta)) * XLALSimInspiralWaveformParamsLookupNonGRDChi7NS(LALpars);
        pfa->v[7] += ((-8980424995.L/762048.L + 6586595.L*eta/756.L - 305.L*eta*eta/36.L)*SL - (170978035.L/48384.L - 2876425.L*eta/672.L - 4735.L*eta*eta/144.L) * dSigmaL) * XLALSimInspiralWaveformParamsLookupNonGRDChi7S(LALpars);
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
                                           const REAL8Sequence *freqs, /**< frequency series */
                                           const UINT4 m, /**< order m of GW mode */
                                           const UINT4 iStart, /**< index of starting frequency */
                                           const UINT4 iRef, /**< index of reference frequency */
                                           const UINT4 iPeak, /**< index of peak frequency */
                                           PNPhasingSeries pfa, /**< PN phasing coefficients */
                                           const REAL8 mtot, /**< total mass */
                                           const REAL8 eta, /**< dimensionless symmetric mass ratio */
                                           const REAL8 f_window_div_f_Peak, /**< Frequency at which to attach non-GR and GR waveforms, inputted as a fraction of f_Peak (should be between 0 and 1) */
                                           const REAL8 iStartFinal, /**< index of starting frequency for TGR modification */
                                           const REAL8 NCyclesStep  /**< Choose number of GW cycles over which to taper the non-GR phase correction */
		                           )
{
    
  UINT4 i;
  const REAL8 pfaN0 = 3.L/(128.L * eta);
  const REAL8 piM = LAL_PI * mtot;
  /*const REAL8 vWindow = cbrt(2*piM * freqs->data[iEnd]/m); //Center of the tapering function in v-space*/
  const REAL8 vWindow = cbrt(2*piM * f_window_div_f_Peak * freqs->data[iPeak]/2); //Center of the tapering function in v-space
  const REAL8 width = (NCyclesStep * LAL_PI * vWindow * vWindow * vWindow * vWindow * vWindow * vWindow)/(50. * pfaN0); //Width of tapering step function

  
  /* Compute the second derivative of the phase correction and taper with step function
   * Phase correction only generated from f_low to f_max
   */
  REAL8Sequence *d2phasenonGRdf2Tapered = NULL;
  d2phasenonGRdf2Tapered = XLALCreateREAL8Sequence( (UINT4) freqs->length );
  REAL8 f;
  for ( i = iStart; i < freqs->length; i++)
    {
      f = freqs->data[i];
      if (f>0) d2phasenonGRdf2Tapered->data[i] = PNPhaseSecondDerivative(f, m, pfa, mtot)/ (1. + exp( (cbrt(2*piM * freqs->data[i]/m) - vWindow) / width ) );
    }
 
  /* Interpolate and integrate (twice) the tapered second derivative of the phase correction to compute the phase correction
   * Set the value of the correction and its first derivative to zero at f_ref
   */
  REAL8Sequence *dphasenonGRdfTapered = NULL, *phasenonGRTapered = NULL;
  dphasenonGRdfTapered = XLALCreateREAL8Sequence( (UINT4) freqs->length );
  phasenonGRTapered = XLALCreateREAL8Sequence( (UINT4) freqs->length );

  gsl_spline *splineTemp = NULL;
  gsl_interp_accel *acc = NULL;
  splineTemp = gsl_spline_alloc (gsl_interp_linear, freqs->length);
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
  /*printf("phase at fref=%f\n", phasenonGRTapered->data[iRef]);*/
  
  for ( i = iStartFinal; i < freqs->length; i++ ) {
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

REAL8 PNPhase(REAL8 f,    /* frequency in Hz */
              UINT4 m,
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
    
    const REAL8 v = cbrt(2*piM*f/m);
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
    phasing *=m/2.;

    return -phasing;
}

/* Derivative of phase with respect to f:
 * d[ v^(n - 5) ] = (pi * M) / 3 * (n - 5) v^(n - 8)
 * d[ v^(n - 5) * log(v)] = (pi * M) / 3 * (1 + (n - 5) * log(v)) * v^(n - 8)
 *
 * where v = (pi * M * f)^(1 / 3)
 */

REAL8 PNPhaseDerivative(REAL8 f,    /* frequency in Hz */
                        UINT4 m,
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

    const REAL8 v = cbrt(2*piM*f/m);
    const REAL8 logv = log(v);
    const REAL8 v2 = v * v;
    const REAL8 v3 = v * v2;
    const REAL8 v4 = v * v3;
    const REAL8 v5 = v * v4;
//    const REAL8 v6 = v * v5;
//    const REAL8 v7 = v * v6;
    REAL8 phasing = 0.0;

    phasing += 2. * pfa7 * v4;
    phasing += (1. * pfa6 + 1. * pfl6 + 1. * pfl6 * logv) * v3;
    phasing += 1. * pfl5 * v2;
    phasing += -1. * pfa4 * v;
    phasing += -2. * pfa3;
    phasing += -3. * pfa2 / v; /* there was a mistake here in the factor. it should be -3 instead of -1.*/
    phasing += -4. * pfa1 / v2;
    phasing += -5. * pfaN / v3;
    phasing += -6. * pfaMinus1 / v4;
    phasing += -7. * pfaMinus2 /v5;
    phasing /= v5;
    phasing *= piM/3.;

    return -phasing;
}

/* Derivative of phase with respect to f
 * d2[v^(n - 5)] = (pi * M)^2 / 9 * (n - 5) * (n - 8) * v^(n - 11)
 * d2[v^(n - 5) * log(v)] = (pi * M)^2 / 9 * ( (2 * n - 13) + (n - 5) * (n - 8) * log(v) ) * v^(n - 11)
 *
 * where v = (pi * M * f)^(1 / 3)
 */

REAL8 PNPhaseSecondDerivative(REAL8 f,    /* frequency in Hz */
                              UINT4 m,
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

    const REAL8 v = cbrt(2*piM*f/m);
    const REAL8 logv = log(v);
    const REAL8 v2 = v * v;
    const REAL8 v3 = v * v2;
    const REAL8 v4 = v * v3;
    const REAL8 v5 = v * v4;
    const REAL8 v6 = v * v5;
    const REAL8 v7 = v * v6;
    const REAL8 v8 = v * v7;
    REAL8 phasing = 0.0;

    phasing += -2. *  pfa7 * v;
    phasing += (-2. *  pfa6 - 1. * pfl6 - 2. * pfl6 * logv);
    phasing += -3. * pfl5 / v;
    phasing += 4. * pfa4 / v2 ;
    phasing += 10. * pfa3 / v3;
    phasing += 18. * pfa2 / v4;
    phasing += 28. * pfa1 / v5;
    phasing += 40. * pfaN / v6;
    phasing += 54. * pfaMinus1 / v7;
    phasing += 70. * pfaMinus2 /v8;
    phasing /= v5;
    phasing *= piM*piM/9.;
    phasing *=2./m;
    return -phasing;
}
