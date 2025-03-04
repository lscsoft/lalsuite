/*
 * Copyright (C) 2022 Cardiff University
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

#ifdef __cplusplus
extern "C"
{
#endif

/* LALSimulation */
#include "LALSimIMR.h"

/* Spherical harmonics */
#include <lal/SphericalHarmonics.h>
#include <lal/LALConstants.h>


/* Phenom */
#include "LALSimIMRPhenomX.h"
#include "LALSimIMRPhenomX_internals.h"
#include "LALSimIMRPhenomXUtilities.h"
#include "LALSimIMRPhenomX_precession.h"
#include "LALSimIMRPhenomUtils.h"
#include "LALSimIMRPhenomX_ringdown.h"
#include "LALSimIMRPhenomX_intermediate.h"
#include "LALSimIMRPhenomX_inspiral.h"
#include "LALSimIMRPhenomX_PNR_alpha.h"


#include "LALSimIMRPhenomX_AntisymmetricWaveform.h"

IMRPhenomX_UsefulPowers powers_of_lalpi;


#ifndef _OPENMP
#define omp ignore
#endif

#ifndef PHENOMXHMDEBUG
#define DEBUG 0
#else
#define DEBUG 1
#endif

  /**
   * @author Shrobana Ghosh
   **
   *   EXTERNAL GENERATE antisymmetric waveform
   * This is an external wrapper to generate the (2,2) and (2,-2) antisymmetric waveform,
   * with the standard inputs given to generate FD waveforms.
   * @note At present this is only compatible with the PNR angles (refer arxiv 2310.16980)
   */
  int XLALSimIMRPhenomX_PNR_GenerateAntisymmetricWaveform(
      REAL8Sequence **antisymamp, /**< [out] Amplitude of antisymmetric (2,2) waveform */
      REAL8Sequence **antisymphase, /**< [out] Phase of antisymmetric (2,2) waveform */
      REAL8 m1_SI,              /**< mass of companion 1 (kg) */
      REAL8 m2_SI,              /**< mass of companion 2 (kg) */
      REAL8 chi1x,              /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
      REAL8 chi1y,              /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
      REAL8 chi1z,              /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
      REAL8 chi2x,              /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
      REAL8 chi2y,              /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
      REAL8 chi2z,              /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
      REAL8 distance,            /**< distance to source **/
      REAL8 inclination,        /**< Angle between orbital angular momentum and line-of-sight vector at reference frequency (rad) */
      REAL8 deltaF,             /**< Frequency spacing (Hz) */
      REAL8 f_min,              /**< Starting GW frequency (Hz) */
      REAL8 f_max,              /**< Ending GW frequency (Hz) */
      REAL8 fRef_In,            /**< Reference frequency (Hz) */
      REAL8 phiRef,             /**< phase at reference frequency (Hz) */
      LALDict *lalParams        /**< LAL Dictionary struct */
  )
  {
    /* Simple check on masses and spins */
    UINT4 status = 0;
    status = XLALIMRPhenomXPCheckMassesAndSpins(&m1_SI, &m2_SI, &chi1x, &chi1y, &chi1z, &chi2x, &chi2y, &chi2z);
    XLAL_CHECK(
        XLAL_SUCCESS == status,
        XLAL_EFUNC,
        "Error: XLALIMRPhenomXPCheckMassesAndSpins failed in XLALSimIMRPhenomX_PNR_GenerateAntisymmetricWaveform.\n");

    /* Ensure we have a dictionary */

    status = IMRPhenomX_Initialize_Powers(&powers_of_lalpi, LAL_PI);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PI.\n");

    LALDict *lalParams_aux;
    if (lalParams == NULL)
    {
      lalParams_aux = XLALCreateDict();
    }
    else
    {
      lalParams_aux = XLALDictDuplicate(lalParams);
    }

    XLALSimInspiralWaveformParamsInsertPhenomXAntisymmetricWaveform(lalParams_aux, 1);
    if(!XLALSimInspiralWaveformParamsLookupPhenomXPNRUseTunedAngles(lalParams_aux))
    {
      XLAL_PRINT_WARNING("Warning:Antisymmetric waveform generation currently not supported without PNR angles. Turning on PNR angles ... \n");
      XLALSimInspiralWaveformParamsInsertPhenomXPNRUseTunedAngles(lalParams_aux,1);
    }

    /* Map fRef to the start frequency if need be, then make sure it's within the frequency range */
    REAL8 fRef = (fRef_In == 0.0) ? f_min : fRef_In;

    /* Generate a uniformly sampled frequency grid of spacing deltaF. */
    /* Frequencies will be set using only the lower and upper bounds that we passed */
    size_t iStart = (size_t)(f_min / deltaF);
    size_t iStop = (size_t)(f_max / deltaF) + 1;

    XLAL_CHECK(
        (iStart <= iStop),
        XLAL_EDOM,
        "Error: the starting frequency index is greater than the stopping index! Please ensure that f_min <= f_max.\n");


    REAL8Sequence *freqs = NULL;
    freqs = XLALCreateREAL8Sequence(iStop - iStart);
    *antisymamp = XLALCreateREAL8Sequence(iStop - iStart);
    *antisymphase = XLALCreateREAL8Sequence(iStop - iStart);

    for (UINT4 i = iStart; i < iStop; i++)
    {
      freqs->data[i - iStart] = i * deltaF;
    }

    IMRPhenomXWaveformStruct *pWF;
    pWF = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
    status = IMRPhenomXSetWaveformVariables(pWF, m1_SI, m2_SI, chi1z, chi2z, deltaF, fRef, phiRef, f_min, f_max, distance, inclination, lalParams_aux, PHENOMXDEBUG);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");

    IMRPhenomXPrecessionStruct *pPrec;
    pPrec = XLALMalloc(sizeof(IMRPhenomXPrecessionStruct));
    status = IMRPhenomXGetAndSetPrecessionVariables(pWF, pPrec, m1_SI, m2_SI, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, lalParams_aux, DEBUG);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetPrecessionVariables failed.\n");

    IMRPhenomX_PNR_GenerateAntisymmetricWaveform(*antisymamp,*antisymphase,freqs,pWF,pPrec,lalParams_aux);

    LALFree(pWF);
    if(pPrec->pWF22AS){
      LALFree(pPrec->pWF22AS);
    }
    LALFree(pPrec);
    XLALDestroyDict(lalParams_aux);
    XLALDestroyREAL8Sequence(freqs);

    return XLAL_SUCCESS;
  }

  /**
   *  EXTERNAL GENERATE Antisymmetric Amplitude Ratio
   * This is an external wrapper to generate the (2,2) antisymmetric amplitude ratio
   * amplitude given the standard inputs given to generate FD waveforms.
   */
  int XLALSimIMRPhenomX_PNR_GenerateAntisymmetricAmpRatio(
      REAL8Sequence **kappa,    /**< [out] Antisymmetric amplitude ratio */
      REAL8Sequence **freqs,    /**< [out] Frequency array (Hz) */
      REAL8 m1_SI,              /**< mass of companion 1 (kg) */
      REAL8 m2_SI,              /**< mass of companion 2 (kg) */
      REAL8 chi1x,              /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
      REAL8 chi1y,              /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
      REAL8 chi1z,              /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
      REAL8 chi2x,              /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
      REAL8 chi2y,              /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
      REAL8 chi2z,              /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
      REAL8 inclination,        /**< Angle between orbital angular momentum and line-of-sight vector at reference frequency (rad) */
      REAL8 deltaF,             /**< Frequency spacing (Hz) */
      REAL8 f_min,              /**< Starting GW frequency (Hz) */
      REAL8 f_max,              /**< Ending GW frequency (Hz) */
      REAL8 fRef_In,            /**< Reference frequency (Hz) */
      LALDict *lalParams        /**< LAL Dictionary struct */
  )
  {
    /* Simple check on masses and spins */
    UINT4 status = 0;
    status = XLALIMRPhenomXPCheckMassesAndSpins(&m1_SI, &m2_SI, &chi1x, &chi1y, &chi1z, &chi2x, &chi2y, &chi2z);
    XLAL_CHECK(
        XLAL_SUCCESS == status,
        XLAL_EFUNC,
        "Error: XLALIMRPhenomXPCheckMassesAndSpins failed in XLALSimIMRPhenomX_PNR_GenerateAntisymmetricAmpRatio.\n");

    /* Ensure we have a dictionary */
    LALDict *lalParams_aux;
    if (lalParams == NULL)
    {
      lalParams_aux = XLALCreateDict();
    }
    else
    {
      lalParams_aux = XLALDictDuplicate(lalParams);
    }

    XLALSimInspiralWaveformParamsInsertPhenomXAntisymmetricWaveform(lalParams_aux, 1);
    if(!XLALSimInspiralWaveformParamsLookupPhenomXPNRUseTunedAngles(lalParams_aux))
    {
      XLAL_PRINT_WARNING("Warning:Antisymmetric waveform generation currently not supported without PNR angles. Turning on PNR angles ... \n");
      XLALSimInspiralWaveformParamsInsertPhenomXPNRUseTunedAngles(lalParams_aux,1);
    }

    /* Map fRef to the start frequency if need be, then make sure it's within the frequency range */
    REAL8 fRef = (fRef_In == 0.0) ? f_min : fRef_In;

    /* Generate a uniformly sampled frequency grid of spacing deltaF. */
    /* Frequencies will be set using only the lower and upper bounds that we passed */
    size_t iStart = (size_t)(f_min / deltaF);
    size_t iStop = (size_t)(f_max / deltaF) + 1;

    XLAL_CHECK(
        (iStart <= iStop),
        XLAL_EDOM,
        "Error: the starting frequency index is greater than the stopping index! Please ensure that f_min <= f_max.\n");

    /* Allocate memory for frequency array and terminate if this fails */
    *freqs = XLALCreateREAL8Sequence(iStop - iStart);
    *kappa = XLALCreateREAL8Sequence(iStop - iStart);

    /* Specify arbitrary parameters for the angle generation */
    REAL8 distance = 1.0;
    REAL8 phiRef = 0.0;
    // REAL8 inclination = 0.0;

    for (UINT4 i = iStart; i < iStop; i++)
    {
      (*freqs)->data[i - iStart] = i * deltaF;
    }

    /* Use the true minimum and maximum frequency values */
    REAL8 f_min_eval = (*freqs)->data[0];
    REAL8 f_max_eval = (*freqs)->data[(*freqs)->length - 1];

    /* Initialize PhenomX Waveform struct and check that it initialized correctly */
    IMRPhenomXWaveformStruct *pWF;
    pWF = XLALMalloc(sizeof(IMRPhenomXWaveformStruct));
    status = IMRPhenomXSetWaveformVariables(pWF, m1_SI, m2_SI, chi1z, chi2z, deltaF, fRef, phiRef, f_min_eval, f_max_eval, distance, inclination, lalParams_aux, DEBUG);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetWaveformVariables failed.\n");

    IMRPhenomXPrecessionStruct *pPrec;
    pPrec = XLALMalloc(sizeof(IMRPhenomXPrecessionStruct));
    status = IMRPhenomXGetAndSetPrecessionVariables(pWF, pPrec, m1_SI, m2_SI, chi1x, chi1y, chi1z, chi2x, chi2y, chi2z, lalParams_aux, DEBUG);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXSetPrecessionVariables failed.\n");

    status = IMRPhenomX_PNR_GenerateAntisymmetricAmpRatio(*kappa, *freqs, pWF, pPrec);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomX_Generate_AntisymmetricAmpRatio failed.\n");

    /* Clean up memory allocation */
    LALFree(pWF);
    if(pPrec->pWF22AS){
      LALFree(pPrec->pWF22AS);
    }
    LALFree(pPrec);
    XLALDestroyDict(lalParams_aux);

    return XLAL_SUCCESS;
  }

  int IMRPhenomX_PNR_GenerateAntisymmetricWaveform(
      REAL8Sequence *antisymamp, /**< [out] Amplitude of antisymmetric (2,2) waveform */
      REAL8Sequence *antisymphase, /**< [out] Phase of antisymmetric (2,2) waveform */
      const REAL8Sequence *freqs,        /**< input frequency array (Hz) */
      IMRPhenomXWaveformStruct *pWF,     /**< waveform struct */
      IMRPhenomXPrecessionStruct *pPrec, /**< precession struct **/
      LALDict *lalParams        /**< LAL Dictionary struct */
  )
  {
    UINT4 status = 0;
    status = IMRPhenomX_Initialize_Powers(&powers_of_lalpi, LAL_PI);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initialize useful powers of LAL_PI.\n");

    IMRPhenomXPhaseCoefficients *pPhase22;
    pPhase22 = XLALMalloc(sizeof(IMRPhenomXPhaseCoefficients));
    status   = IMRPhenomXGetPhaseCoefficients(pWF,pPhase22);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXGetPhaseCoefficients failed.\n");

    IMRPhenomX_Phase_22_ConnectionCoefficients(pWF,pPhase22);

    IMRPhenomXAmpCoefficients *pAmp22;
    pAmp22 = XLALMalloc(sizeof(IMRPhenomXAmpCoefficients));
    status = IMRPhenomXGetAmplitudeCoefficients(pWF,pAmp22);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomXGetAmplitudeCoefficients failed.\n");

    REAL8Sequence *kappa = NULL;
    kappa = XLALCreateREAL8Sequence(freqs->length);

    status = IMRPhenomX_PNR_GenerateAntisymmetricAmpRatio(kappa, freqs, pWF, pPrec);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomX_Generate_AntisymmetricAmpRatio failed.\n");

    double MfT = 0.85 * pWF->fRING;
    double lina = 0.0;
    double linb = IMRPhenomX_TimeShift_22(pPhase22, pWF);
    REAL8 inveta    = (1.0 / pWF->eta);

    REAL8 A0 = 0.0;
    REAL8 phi_A0 = 0.0;
    REAL8 phi_B0 = 0.0;

    status = IMRPhenomX_PNR_GenerateAntisymmetricPhaseCoefficients(&A0, &phi_A0, &phi_B0, MfT, lina, linb, inveta, pWF, pPrec,pPhase22);
    XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomX_Generate_AntisymmetricPhaseCoefficients failed.\n");

    REAL8 fPhaseIN  = pPhase22->fPhaseMatchIN;
    REAL8 fPhaseIM  = pPhase22->fPhaseMatchIM;
    REAL8 fAmpIN    = pAmp22->fAmpMatchIN;
    REAL8 fAmpIM    = pAmp22->fAmpRDMin;
    REAL8 C1IM      = pPhase22->C1Int;
    REAL8 C2IM      = pPhase22->C2Int;
    REAL8 C1RD      = pPhase22->C1MRD;
    REAL8 C2RD      = pPhase22->C2MRD;

    REAL8Sequence *alphaPNR = NULL;

    alphaPNR = XLALCreateREAL8Sequence(freqs->length);

    status = IMRPhenomX_PNR_GeneratePNRAlphaForAntisymmetry(alphaPNR,
       freqs, pWF, pPrec,
       lalParams);
     XLAL_CHECK(XLAL_SUCCESS == status, XLAL_EFUNC, "Error: IMRPhenomX_PNR_GeneratePNRAngles failed.\n");

      IMRPhenomX_UsefulPowers powers_of_MfRef;
      IMRPhenomX_Initialize_Powers(&powers_of_MfRef,pWF->MfRef);

     REAL8 phiref22 = -inveta*IMRPhenomX_Phase_22(pWF->MfRef, &powers_of_MfRef, pPhase22, pWF) - linb*pWF->MfRef - lina + 2.0*pWF->phi0 + LAL_PI_4;

    for (UINT4 idx = 0; idx <freqs->length; idx++)
    {
      double Mf    = pWF->M_sec * freqs->data[idx];
      REAL8 amp = 0.0;
      REAL8 phi = 0.0;
      REAL8 phi_AS = 0.0;
      REAL8 amp_AS = 0.0;

      IMRPhenomX_UsefulPowers powers_of_Mf;
      IMRPhenomX_Initialize_Powers(&powers_of_Mf,Mf);

      /* Get symmetric phase */
      if(Mf < fPhaseIN)
      {
        phi = IMRPhenomX_Inspiral_Phase_22_AnsatzInt(Mf, &powers_of_Mf, pPhase22);
      }
      else if(Mf > fPhaseIM)
      {
        phi = IMRPhenomX_Ringdown_Phase_22_AnsatzInt(Mf, &powers_of_Mf, pWF, pPhase22) + C1RD + (C2RD * Mf);
      }
      else
      {
        phi = IMRPhenomX_Intermediate_Phase_22_AnsatzInt(Mf, &powers_of_Mf, pWF, pPhase22) + C1IM + (C2IM * Mf);
      }
      /* Scale phase by 1/eta and apply phase and time shifts */
      phi  *= inveta;
      phi  += linb*Mf + lina + phiref22;

      /* Get smmetric amplitude */
      if(Mf < fAmpIN)
      {
        amp = IMRPhenomX_Inspiral_Amp_22_Ansatz(Mf, &powers_of_Mf, pWF, pAmp22);
      }
      else if(Mf > fAmpIM)
      {
        amp = IMRPhenomX_Ringdown_Amp_22_Ansatz(Mf, pWF, pAmp22);
      }
      else
      {
        amp = IMRPhenomX_Intermediate_Amp_22_Ansatz(Mf, &powers_of_Mf, pWF, pAmp22);
      }

      /****** antisymmetric amplitude ******/
      amp_AS = kappa->data[idx] * amp;


      if(Mf < MfT)
      {
        phi_AS = phi/2 + alphaPNR->data[idx] + A0*Mf + phi_A0;
      }
      else
      {
        phi_AS = phi + phi_B0;
      }

      REAL8 Amp0 = pWF->amp0 * pWF->ampNorm;

      antisymamp->data[idx] = Amp0 * powers_of_Mf.m_seven_sixths * amp_AS;
      antisymphase->data[idx] = pPrec->zeta_polarization + phi_AS;
    }

    /* Clean up memory allocation */

    LALFree(pPhase22);
    LALFree(pAmp22);

    XLALDestroyREAL8Sequence(kappa);
    XLALDestroyREAL8Sequence(alphaPNR);

    return XLAL_SUCCESS;
  }

  /**
   *
   *
   *   GENERATE Antisymmetric amplitude ratio
   *
   */
  int IMRPhenomX_PNR_GenerateAntisymmetricAmpRatio(
    REAL8Sequence *kappa,              /**< [out] antisymmetric amplitude ratio */
    const REAL8Sequence *freqs,        /**< input frequency array (Hz) */
    IMRPhenomXWaveformStruct *pWF,     /**< waveform struct */
    IMRPhenomXPrecessionStruct *pPrec /**< precession struct **/
  )
  {

  /**** populating parameters of binary from pWF ******/
    const double m2 = pWF->m2;
    const double M = pWF->Mtot;
    const REAL8 delta = 1 - 2*m2;
    const REAL8 eta = pWF->eta;
    const double MfRD = pWF->fRING;
    const REAL8 theta = pPrec->theta_antisymmetric;
    const REAL8 Chi = pPrec->chi_singleSpin_antisymmetric;

  /* coefficients for phenemenological fit of amplitude ratio in PN-NR transition*/
    double b0 = 18.0387;
    double b1 = 15.4509;
    double b2 = 55.1140;
    double b3 = -203.6290;

    double b = b0 + b1*eta + b2*theta + b3*eta*theta;

    REAL8 vRD = cbrt (LAL_PI * MfRD );
    const double kappaRD = GetKappa_at_frequency(vRD,delta,Chi,theta,eta,b);

    for (size_t i = 0; i < freqs->length; i++)
    {
      REAL8 Mf = XLALSimPhenomUtilsHztoMf(freqs->data[i], M);
      REAL8 v = cbrt (LAL_PI * Mf );
      if (Mf < MfRD)
      {
        kappa->data[i] = GetKappa_at_frequency(v,delta,Chi,theta,eta,b);
      }
      else
      {
        kappa->data[i] = kappaRD;
      }
    }

    size_t width = 80;
    double df = 0.0;
    if(width > kappa->length - 1)
    {
      width = (size_t)floor((double)(kappa->length) / 2.0);
    }
    size_t half_width = (size_t)floor((double)(width / 2.0));


    for (size_t id = 0; id < kappa->length-width-1; id++)
    {
      double smoothed_ratio = 0.0;
      double frequency_width = 0.0;
      for (size_t j = 0; j < width+1; j++)
      {
        df = freqs->data[id+j+1] - freqs->data[id+j];
        smoothed_ratio += (kappa->data[id+j]*df);
        frequency_width += df;
      }
      kappa->data[id+half_width] =  smoothed_ratio *1.0/frequency_width;
    }

    return XLAL_SUCCESS;
  }

  double GetKappa_at_frequency(REAL8 v,REAL8 delta,REAL8 Chi,REAL8 theta,REAL8 eta,double b)
  {
    REAL8 v2 = v * v;
    REAL8 v3 = v2 * v;
    REAL8 v5 = v3 * v2;
    REAL8 kappaPNnum = (21 * v2 * (1 + delta) * Chi * sin(theta));
    REAL8 kappaPNden = 2 * (42 + 84 * LAL_PI * v3 + v2 * (55 * eta - 107)  - 28 * v3 * (1 + delta - eta) * Chi * cos ( theta ));
    REAL8 amp_ratio = kappaPNnum/kappaPNden * (1 + b * v5 );
    return amp_ratio;
  }

  /**
   * Anti-symmetric phase coefficients/offsets
   *
   */
  int IMRPhenomX_PNR_GenerateAntisymmetricPhaseCoefficients(
    REAL8 *A0, /**< [out] A0 parameter */
    REAL8 *phi_A0, /**< [out] phi_A0 parameter */
    REAL8 *phi_B0, /**< [out] phi_B0 parameter */
    const double MfT, /**< Geometric transition frequency */
    double lina, /**< lina parameter */
    double linb, /**< linb parameter */
    double inveta, /**< inveta parameter */
    IMRPhenomXWaveformStruct *pWF,     /**< waveform struct */
    IMRPhenomXPrecessionStruct *pPrec, /**< precession struct **/
    IMRPhenomXPhaseCoefficients *pPhase22 /**< symmetric phase coefficients struct */
  )
  {
    UINT4 status = 0;
    IMRPhenomX_PNR_alpha_parameters* alphaParams = NULL;
    alphaParams = XLALMalloc(sizeof(IMRPhenomX_PNR_alpha_parameters));

    status = IMRPhenomX_PNR_precompute_alpha_coefficients(alphaParams, pWF, pPrec);
    XLAL_CHECK(
        XLAL_SUCCESS == status,
        XLAL_EFUNC,
        "Error: IMRPhenomX_PNR_precompute_alpha_coefficients failed.\n");

    status = IMRPhenomX_PNR_alpha_connection_parameters(alphaParams, pWF, pPrec);
    XLAL_CHECK(
        XLAL_SUCCESS == status,
        XLAL_EFUNC,
        "Error: IMRPhenomX_PNR_alpha_connection_parameters failed.\n");

    REAL8 phi_der_MfT = 0.0;
    REAL8 phi_MfT = 0.0;
    REAL8 alpha_MfT = 0.0;
    REAL8 alpha_der_MfT = 0.0;

    REAL8 fPhaseIN  = pPhase22->fPhaseMatchIN;
    REAL8 fPhaseIM  = pPhase22->fPhaseMatchIM;
    REAL8 C1IM      = pPhase22->C1Int;
    REAL8 C2IM      = pPhase22->C2Int;
    REAL8 C1RD      = pPhase22->C1MRD;
    REAL8 C2RD      = pPhase22->C2MRD;

    IMRPhenomX_UsefulPowers powers_of_MfT;
    IMRPhenomX_Initialize_Powers(&powers_of_MfT, MfT);

    const double deltaMF = 0.0005;
    const double Mf_right =  MfT + deltaMF;
    const double Mf_left =  MfT - deltaMF;

    REAL8 q = pWF->q;
    REAL8 chi = pPrec->chi_singleSpin;

    /* inside callibration region */
    if ((q <= pPrec->PNR_q_window_lower) && (chi <= pPrec->PNR_chi_window_lower))
    {
      alpha_der_MfT = ( IMRPhenomX_PNR_GeneratePNRAlphaAtMf(Mf_right, alphaParams, pWF, pPrec) - IMRPhenomX_PNR_GeneratePNRAlphaAtMf(Mf_left, alphaParams, pWF, pPrec) )/2/deltaMF;
      alpha_MfT = IMRPhenomX_PNR_GeneratePNRAlphaAtMf(MfT, alphaParams, pWF, pPrec);
    }
    /* inside transition region */
    else if ((q <= pPrec->PNR_q_window_upper) && (chi <= pPrec->PNR_chi_window_upper))
    {
      alpha_der_MfT = ( IMRPhenomX_PNR_GenerateMergedPNRAlphaAtMf(Mf_right, alphaParams, pWF, pPrec) - IMRPhenomX_PNR_GenerateMergedPNRAlphaAtMf(Mf_left, alphaParams, pWF, pPrec) )/2/deltaMF;
      alpha_MfT = IMRPhenomX_PNR_GenerateMergedPNRAlphaAtMf(MfT, alphaParams, pWF, pPrec);
    }
    /* fully in outside calibration region */
    else
    {
      alpha_der_MfT = ( IMRPhenomX_PNR_GetPNAlphaAtFreq(Mf_right, pWF, pPrec) - IMRPhenomX_PNR_GetPNAlphaAtFreq(Mf_left, pWF, pPrec) )/2/deltaMF;
      alpha_MfT = IMRPhenomX_PNR_GetPNAlphaAtFreq(MfT, pWF, pPrec);
    }

    // alpha_der_MfT = ( IMRPhenomX_PNR_GeneratePNRAlphaAtMf(Mf_right, alphaParams, pWF, pPrec) - IMRPhenomX_PNR_GeneratePNRAlphaAtMf(Mf_left, alphaParams, pWF, pPrec) )/2/deltaMF;
    // alpha_MfT = IMRPhenomX_PNR_GeneratePNRAlphaAtMf(MfT, alphaParams, pWF, pPrec);

    if(MfT < fPhaseIN)
    {
      phi_der_MfT = IMRPhenomX_Inspiral_Phase_22_Ansatz(MfT, &powers_of_MfT, pPhase22);
      phi_MfT = IMRPhenomX_Inspiral_Phase_22_AnsatzInt(MfT, &powers_of_MfT, pPhase22);
    }
    else if(MfT > fPhaseIM)
    {
      phi_der_MfT = IMRPhenomX_Ringdown_Phase_22_Ansatz(MfT, &powers_of_MfT, pWF, pPhase22) + C2RD;
      phi_MfT = IMRPhenomX_Ringdown_Phase_22_AnsatzInt(MfT, &powers_of_MfT, pWF, pPhase22) + C1RD + (C2RD * MfT);
    }
    else
    {
      phi_der_MfT = IMRPhenomX_Intermediate_Phase_22_Ansatz(MfT, &powers_of_MfT, pWF, pPhase22) + C2IM;
      phi_MfT = IMRPhenomX_Intermediate_Phase_22_AnsatzInt(MfT, &powers_of_MfT, pWF, pPhase22) + C1IM + (C2IM * MfT);
    }

    /* Scale phase by 1/eta and apply phase and time shifts */
      IMRPhenomX_UsefulPowers powers_of_MfRef;
      IMRPhenomX_Initialize_Powers(&powers_of_MfRef,pWF->MfRef);
     REAL8 phiref22 = -inveta*IMRPhenomX_Phase_22(pWF->MfRef, &powers_of_MfRef, pPhase22, pWF) - linb*pWF->MfRef - lina + 2.0*pWF->phi0 + LAL_PI_4;

    phi_der_MfT *= inveta;
    phi_der_MfT += linb;
    phi_MfT  *= inveta;
    phi_MfT  += linb*MfT + lina + phiref22;

    *A0 = phi_der_MfT/2 - alpha_der_MfT;
    *phi_A0 = pPrec-> alpha_offset;
    *phi_B0 = alpha_MfT - phi_MfT/2 + *A0 * MfT + *phi_A0;

    LALFree(alphaParams);

    return XLAL_SUCCESS;
  }

#ifdef __cplusplus
}
#endif
