/*
 * Copyright (C) 2019 Marta Colleoni, Cecilio Garcia Quiros
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
//
//  LALSimIMRPhenomXHM_internals.c
//
//  Created by Marta on 06/02/2019.
//

#include <gsl/gsl_linalg.h>

#include "LALSimIMRPhenomXHM_structs.h"

#include "LALSimIMRPhenomXHM_qnm.h"
#include "LALSimIMRPhenomXHM_qnm.c"

#include "LALSimIMRPhenomX_internals.h"
#include "LALSimIMRPhenomXUtilities.h"


#include "LALSimIMRPhenomXHM_inspiral.c"

#include "LALSimIMRPhenomXHM_ringdown.c"
#include "LALSimIMRPhenomXHM_intermediate.c"
#include "LALSimIMRPhenomX_PNR_deviations.h"


/*Equations referenced in this file come from arXiv:2001.10914 [gr-qc], see also dcc:LIGO-P2000011 */

/* Intialize Quasi Normal Modes' ringdown and damping frequencies */
void IMRPhenomXHM_Initialize_QNMs(QNMFits *qnms){

  // pointers to fring fits
  qnms->fring_lm[0]=evaluate_QNMfit_fring21;
  qnms->fring_lm[1]=evaluate_QNMfit_fring33;
  qnms->fring_lm[2]=evaluate_QNMfit_fring32;
  qnms->fring_lm[3]=evaluate_QNMfit_fring44;

  //pointers to fdamp fits
  qnms->fdamp_lm[0]=evaluate_QNMfit_fdamp21;
  qnms->fdamp_lm[1]=evaluate_QNMfit_fdamp33;
  qnms->fdamp_lm[2]=evaluate_QNMfit_fdamp32;
  qnms->fdamp_lm[3]=evaluate_QNMfit_fdamp44;
}

/* Initialize Mixing Coefficients (complex). We assume mixing with just one mode. Only 32 mode has mixing. */
void IMRPhenomXHM_Initialize_MixingCoeffs(IMRPhenomXHMWaveformStruct *wf, IMRPhenomXWaveformStruct *wf22)
{
  wf->mixingCoeffs[0]=evaluate_QNMfit_re_l2m2lp2(wf22->afinal)+I*evaluate_QNMfit_im_l2m2lp2(wf22->afinal);
  wf->mixingCoeffs[1]=evaluate_QNMfit_re_l2m2lp3(wf22->afinal)+I*evaluate_QNMfit_im_l2m2lp3(wf22->afinal);
  wf->mixingCoeffs[2]=evaluate_QNMfit_re_l3m2lp2(wf22->afinal)+I*evaluate_QNMfit_im_l3m2lp2(wf22->afinal);
  wf->mixingCoeffs[3]=evaluate_QNMfit_re_l3m2lp3(wf22->afinal)+I*evaluate_QNMfit_im_l3m2lp3(wf22->afinal);

  // Adjust conventions so that they match the ones used for the hybrids
  wf->mixingCoeffs[2]=-1.* wf->mixingCoeffs[2];
  wf->mixingCoeffs[3]=-1.* wf->mixingCoeffs[3];

  #if DEBUG == 1
  printf("\nMixing coefficients:\n");
  printf("re(a222)=%.16e,im(a222)=%.16e \n",evaluate_QNMfit_re_l2m2lp2(wf22->afinal),evaluate_QNMfit_im_l2m2lp2(wf22->afinal));
  printf("re(a223)=%.16e,im(a223)=%.16e \n",evaluate_QNMfit_re_l2m2lp3(wf22->afinal),evaluate_QNMfit_im_l2m2lp3(wf22->afinal));
  printf("re(a322)=%.16e,im(a322)=%.16e \n",evaluate_QNMfit_re_l3m2lp2(wf22->afinal),evaluate_QNMfit_im_l3m2lp2(wf22->afinal));
  printf("re(a323)=%.16e,im(a323)=%.16e \n",evaluate_QNMfit_re_l3m2lp3(wf22->afinal),evaluate_QNMfit_im_l3m2lp3(wf22->afinal));
  #endif
}

/* Store useful parameters specific of higher modes (thus not in IMRPhenomXWaveformStruct) in IMRPhenomXHMWaveformStruct */
/* E.g.: MECO, ringdown and damping frequencies, the version of the fits, etc. */
void IMRPhenomXHM_SetHMWaveformVariables(
  int ell,
  int emm,
  IMRPhenomXHMWaveformStruct *wf,
  IMRPhenomXWaveformStruct *wf22,
  QNMFits *qnms,
  LALDict *LALParams
)
{

  // read in which mode is being generated
  wf->ell=ell;
  wf->emm=emm;
  wf->modeTag=ell*10+emm;      //21, 33, 32, 44
  wf->ampNorm = wf22->ampNorm;
  wf->fMECOlm = wf22->fMECO*emm*0.5;
  wf->Ampzero = 0;                       // Ampzero = 1 (true) for odd modes and equal black holes
  wf->Amp0 = wf22->amp0;// * wf22->ampNorm;
  wf->useFAmpPN = 0;                     // Only true for the 21, this mode has a different inspiral ansatz
  wf->AmpEMR = 0; // Only one intermediate region
  wf->InspiralAmpVeto=0;
  wf->IntermediateAmpVeto=0;
  wf->RingdownAmpVeto=0;
  
  switch(wf->modeTag){
    case 21:{
      wf->modeInt=0;
      wf->MixingOn=0;
      if(wf22->q == 1. && wf22->chi1L == wf22->chi2L){  // Odd mode for equal mass, equal spin is zero
        wf->Ampzero = 1;
      }
      break;
    }
    case 33:{
      wf->modeInt=1;
      wf->MixingOn=0;
      if(wf22->q == 1. && wf22->chi1L == wf22->chi2L){  // Odd mode for equal mass, equal spin is zero
        wf->Ampzero = 1;
      }
      break;
    }
    case 32:{            //Mode with Mixing
      wf->modeInt=2;
      wf->MixingOn=1;
      break;
    }
    case 44:{
      wf->modeInt=3;
      wf->MixingOn=0;
      break;
    }
    default:
      {XLALPrintError("Error in IMRPhenomXHM_SetHMWaveformVariables: mode (%i,%i) selected is not currently available. Modes available are ((2,|2|),(2,|1|),(3,|2|),(3,|3|),(4,|4|)).\n", wf->ell, wf->emm);}
  }

    /* Here we select the version of the fits and of the reconstruction that will be used in the code.
       Through XLALSimInspiralWaveformParamsLookupPhenomXHM(Inspiral/Intermediate/Ringdown)(Amp/Phase)Version we call the version of the fits used by the phase in each region.
       Currently there is only one version available and is tagged by the release date in the format mmyyyy (122019)*/
  wf->IMRPhenomXHMReleaseVersion             = XLALSimInspiralWaveformParamsLookupPhenomXHMReleaseVersion(LALParams); // 122019

  wf->IMRPhenomXHMInspiralPhaseVersion       = XLALSimInspiralWaveformParamsLookupPhenomXHMInspiralPhaseVersion(LALParams);//122019
  wf->IMRPhenomXHMIntermediatePhaseVersion   = XLALSimInspiralWaveformParamsLookupPhenomXHMIntermediatePhaseVersion(LALParams); //122019
  wf->IMRPhenomXHMRingdownPhaseVersion       = XLALSimInspiralWaveformParamsLookupPhenomXHMRingdownPhaseVersion(LALParams); //122019
  
  wf->IMRPhenomXHMInspiralAmpFitsVersion     = XLALSimInspiralWaveformParamsLookupPhenomXHMInspiralAmpFitsVersion(LALParams); //122018
  wf->IMRPhenomXHMIntermediateAmpFitsVersion = XLALSimInspiralWaveformParamsLookupPhenomXHMIntermediateAmpFitsVersion(LALParams); //122018
  wf->IMRPhenomXHMRingdownAmpFitsVersion     = XLALSimInspiralWaveformParamsLookupPhenomXHMRingdownAmpFitsVersion(LALParams); //122018
  
  wf->IMRPhenomXHMInspiralAmpFreqsVersion    = XLALSimInspiralWaveformParamsLookupPhenomXHMInspiralAmpFreqsVersion(LALParams); //122018
  wf->IMRPhenomXHMIntermediateAmpFreqsVersion= XLALSimInspiralWaveformParamsLookupPhenomXHMIntermediateAmpFreqsVersion(LALParams); //122018
  wf->IMRPhenomXHMRingdownAmpFreqsVersion    = XLALSimInspiralWaveformParamsLookupPhenomXHMRingdownAmpFreqsVersion(LALParams); //122018
 
  /* Reconstruction version for the amplitude */
  wf->IMRPhenomXHMInspiralAmpVersion         = XLALSimInspiralWaveformParamsLookupPhenomXHMInspiralAmpVersion(LALParams); //3  (3 collocation points)
  wf->IMRPhenomXHMIntermediateAmpVersion     = XLALSimInspiralWaveformParamsLookupPhenomXHMIntermediateAmpVersion(LALParams); //2   (2 collocation points)
  wf->IMRPhenomXHMRingdownAmpVersion         = XLALSimInspiralWaveformParamsLookupPhenomXHMRingdownAmpVersion(LALParams); //0  (0 collocation points)

  /* HM tuning is currently only setup for he l=m=3 moment. NOTE that default values for the parameters below are zero (IMRPhenomXPNRUseTunedCoprec) */
  // NOTE that default values are all ZERO. We manually impose this here FOR ALL MODES
  wf->MU1            = 0;
  wf->MU2            = 0;
  wf->MU3            = 0;
  wf->MU4            = 0;
  wf->NU0            = 0;
  wf->NU4            = 0;
  wf->NU5            = 0;
  wf->NU6            = 0;
  wf->ZETA1          = 0;
  wf->ZETA2          = 0;
  if(wf->modeTag==33){

    /* Set parameters for coprecessing frame deviations. NOTE that we add a factor of delta so that all 33 deviations turn off at equal mass ratio -- i.e. there is NO calibration for equal mass ratio cases */
    wf->PNR_DEV_PARAMETER =  wf22->delta * (wf22->PNR_DEV_PARAMETER);
  
    // IF PNR's coprecessing model is wanted, use its fits for the deviation parameters ELSE use the input values with defaults of zeros
    if( wf22->IMRPhenomXPNRUseTunedCoprec33 ) 
    {
      
      /* ------------------------------------------------------ >>
      Get them from the stored model fits that define PhenomXCP 
      within PhenomXPNR
      << ------------------------------------------------------ */
      
      /* MU1 modifies pAmp->lambda */
      wf->MU1     = IMRPhenomXCP_MU1_l3m3(   wf22->theta_LS, wf22->eta, wf22->a1 );
      
      // NOTE that the function for MU2 is not defined in the model 
      /* MU2 would modify pAmp->gamma2 */
      
      /* MU2 */
      wf->MU2     = IMRPhenomXCP_MU2_l3m3(   wf22->theta_LS, wf22->eta, wf22->a1 );
      
      /* MU3 modifies pAmp->gamma3 */
      wf->MU3     = IMRPhenomXCP_MU3_l3m3(   wf22->theta_LS, wf22->eta, wf22->a1 );
      
      /* MU4 modifies V2 or V3 for the intermediate amplitude
      for the DEFAULT value of IMRPhenomXIntermediateAmpVersion
      use in IMRPhenomXPHM */
      wf->MU4     = IMRPhenomXCP_MU4_l3m3(   wf22->theta_LS, wf22->eta, wf22->a1 );
      
      // NOTE that we choose to disable time-shift tuning
      // /* NU0 modifies pPhase->c0 */ 
      // wf->NU0     = IMRPhenomXCP_NU0_l3m3(   wf22->theta_LS, wf22->eta, wf22->a1 );
      
      /* NU4 modifies pPhase->cL */
      wf->NU4     = IMRPhenomXCP_NU4_l3m3(   wf22->theta_LS, wf22->eta, wf22->a1 );
      
      /* NU5 modifies wf->fRING [EXTRAP-PASS-TRUE] */ 
      wf->NU5     = IMRPhenomXCP_NU5_l3m3(   wf22->theta_LS, wf22->eta, wf22->a1 );
      
      /* NU6 modifies wf->fDAMP [EXTRAP-PASS-TRUE] */
      wf->NU6     = IMRPhenomXCP_NU6_l3m3(   wf22->theta_LS, wf22->eta, wf22->a1 );
      
      /* ZETA1 modifies pPhase->b4 */
      wf->ZETA1   = IMRPhenomXCP_ZETA1_l3m3( wf22->theta_LS, wf22->eta, wf22->a1 );
      
      /* ZETA2 modifies pPhase->b1  */
      wf->ZETA2   = IMRPhenomXCP_ZETA2_l3m3( wf22->theta_LS, wf22->eta, wf22->a1 );
      
      /* DEBUGGING: Turn off select deviations */
      // wf->MU1 = 0;
      // wf->MU2 = 0;
      // wf->MU3 = 0;
      // wf->MU4 = 0;
      
    } else {
    
      // NOTE that all default values are all ZERO
      wf->MU1            = XLALSimInspiralWaveformParamsLookupPhenomXCPMU1l3m3(LALParams);
      wf->MU2            = XLALSimInspiralWaveformParamsLookupPhenomXCPMU2l3m3(LALParams);
      wf->MU3            = XLALSimInspiralWaveformParamsLookupPhenomXCPMU3l3m3(LALParams);
      wf->MU4            = XLALSimInspiralWaveformParamsLookupPhenomXCPMU4l3m3(LALParams);
      wf->NU0            = XLALSimInspiralWaveformParamsLookupPhenomXCPNU0l3m3(LALParams);
      wf->NU4            = XLALSimInspiralWaveformParamsLookupPhenomXCPNU4l3m3(LALParams);
      wf->NU5            = XLALSimInspiralWaveformParamsLookupPhenomXCPNU5l3m3(LALParams);
      wf->NU6            = XLALSimInspiralWaveformParamsLookupPhenomXCPNU6l3m3(LALParams);
      wf->ZETA1          = XLALSimInspiralWaveformParamsLookupPhenomXCPZETA1l3m3(LALParams);
      wf->ZETA2          = XLALSimInspiralWaveformParamsLookupPhenomXCPZETA2l3m3(LALParams);
    
    }
    
  }


  wf->IMRPhenomXHMInspiralPhaseFitsVersion = wf->IMRPhenomXHMInspiralPhaseVersion;
  wf->IMRPhenomXHMIntermediatePhaseFitsVersion = wf->IMRPhenomXHMIntermediatePhaseVersion;
  wf->IMRPhenomXHMRingdownPhaseFitsVersion = wf->IMRPhenomXHMRingdownPhaseVersion;

  wf->IMRPhenomXHMInspiralPhaseFreqsVersion = wf->IMRPhenomXHMInspiralPhaseVersion;
  wf->IMRPhenomXHMIntermediatePhaseFreqsVersion = wf->IMRPhenomXHMIntermediatePhaseVersion;
  wf->IMRPhenomXHMRingdownPhaseFreqsVersion = wf->IMRPhenomXHMRingdownPhaseVersion;

  wf->nCollocPtsRDPhase = 0;
  
  if (wf->IMRPhenomXHMReleaseVersion == 122019){
      if(wf22->eta < 0.013886133703630232 && wf22->chi1L<=0.9){ // For q>70 and chi1<0.9 use two intermediate regions.
        wf->AmpEMR = 1;                                         // These cases have a more pronounced drop off in the amplitude at the end of the inspiral and need two intermediate regions to model it.
      }      
      switch(wf->modeTag){
        case 21:{
          if(wf22->q < 8.){
            wf->InspiralAmpVeto=1;
            wf->IntermediateAmpVeto=1;
          }
          wf->RingdownAmpVeto=1;
          if(wf22->eta >= 0.0237954){ //for EMR (q>40.)
            wf->useFAmpPN = 1;
          }
          break;
        }
        case 33:{
          if(wf22->q == 1. && wf22->chi1L == wf22->chi2L){  // Odd mode for equal mass, equal spin is zero
            wf->Ampzero = 1;
          }
          break;
        }
        case 32:{            //Mode with Mixing
          wf->RingdownAmpVeto=1;
          wf->nCollocPtsRDPhase = 4;
          break;
        }
        case 44:{
          break;
        }
        default:
          {XLALPrintError("Error in IMRPhenomXHM_SetHMWaveformVariables: mode (%i,%i) selected is not currently available. Modes available are ((2,|2|),(2,|1|),(3,|2|),(3,|3|),(4,|4|)).\n", wf->ell, wf->emm);}
      }
      wf->nCollocPtsInspAmp = wf->IMRPhenomXHMInspiralAmpVersion;
      wf->nCollocPtsInterAmp = wf->IMRPhenomXHMIntermediateAmpVersion;
      wf->IMRPhenomXHMRingdownPhaseFreqsVersion = 122019;
  }  
  else if (wf->IMRPhenomXHMReleaseVersion == 122022)
  {
      /* Here we change the versions for the collocation point frequencies, parameter space fits and ansatzaes for both amplitude and phase.
         It is better to do it here instead of changing the default values in SimInspiralWaveformParams.c.
         For a new release, copy this else if block and use it as a template. */
      
      /* Amplitude */
      
      // Versions for the cutting and collocation points frequencies
      wf->IMRPhenomXHMInspiralAmpFreqsVersion     = 122022;
      wf->IMRPhenomXHMIntermediateAmpFreqsVersion = 0;
      wf->IMRPhenomXHMRingdownAmpFreqsVersion     = 122022;
      // Versions for the parameter space fits of coefficients or collocation points
      wf->IMRPhenomXHMInspiralAmpFitsVersion     = 122022;
      wf->IMRPhenomXHMIntermediateAmpFitsVersion = 122022;
      wf->IMRPhenomXHMRingdownAmpFitsVersion     = 122022;
      // Versions for the ansatzaes
      wf->IMRPhenomXHMInspiralAmpVersion     = 123;
      wf->IMRPhenomXHMIntermediateAmpVersion = 211112;
      wf->IMRPhenomXHMRingdownAmpVersion     = 2;
      
      /* Notes for the new release and beyond:
         - IMRPhenomXHMInspiralAmpVersion: 
                                           assuming that we have 3 collocation points, it indicates which collocation points are going to be used in the reconstruction.
                                           e.g.: 123 means we use the three, 13 means we only use the 1st and the 3rd; 2 means we only use the 2nd, etc.
         - IMRPhenomXHMIntermediateAmpVersion: 
                                           the number of digits indicate the number of collocation points in the intermediate part.
                                           in this case we have 6 digits which relates to the 4 fitted collocation points + 2 boundaries
                                           if the digit is 0 then we discard that collocation point
                                           if the digit is 1 we use the value of the amplitude at that frequency
                                           if the digit is 2 we use the value and the derivative of the derivative amplitude at that frequency
                                           the digit 2 can only really go in the boundaries since the inner points are fixed parameter space fits of the amplitude value
                                           at the boundaries we can use the inspiral or ringdown ansatz to compute the derivative
                                           e.g.: 211112 means we will use the 4 parameter space fits + 2 value + 2 derivative at the boundaries = 8 degrees of freedom
                                           e.g.: 110012 means we skip the left derivative and the two inner most collocation points = 5 degrees of freedom
         - IMRPhenomXHMRingdownAmpVersion:
                                           this is just an index for the different ansatz versions 
                                           currently it only takes 0 (old release) and 1 (122022 release)
        */
      
      
      // Mode specific tweaks for amplitude
      if (wf->modeTag == 21){
          /* The 21 behaves better if we skip the left the derivative.
             The drawback is that for cases with drop the inspiral-intermediate transition will not be derivable.
             We also drop two collocation points in the middle to avoid wavy behaviour in the extrapolation region */
          wf->IMRPhenomXHMIntermediateAmpVersion = 110102;
      }
      else if (wf->modeTag == 32){
          // Use fit of coefficients instead of collocation points
          wf->IMRPhenomXHMRingdownAmpVersion = 1;
      }
      
      // As mentioned above, the number of digits in IMRPhenomXHMInsp(Intermediate)AmpVersion gives the number of collocation points for the inspiral(intermediate) part
      wf->nCollocPtsInspAmp  = snprintf(NULL, 0, "%i", wf->IMRPhenomXHMInspiralAmpVersion);
      wf->nCollocPtsInterAmp = snprintf(NULL, 0, "%i", wf->IMRPhenomXHMIntermediateAmpVersion);
      
      
      /* Phase */
      /* HM Phase is unchanged respect to the old release */
      if (wf->modeTag == 32){wf->nCollocPtsRDPhase = 4;}
      
      /* Uncommenting the code below will employ a recalibration of the spheroidal phase for the 32.
         This includes the collocation points for the derivative and the extra time and phase shift respect the 22 mode.
         This can introduce introduce significant mismatch in the multimode waveform in the region without NR, so it is left for the next update. */
      // if (wf->modeTag == 32){
      //   wf->IMRPhenomXHMIntermediatePhaseFreqsVersion = 122022;
      //   wf->IMRPhenomXHMRingdownPhaseFreqsVersion = 122022;
      //   wf->IMRPhenomXHMRingdownPhaseFitsVersion  = 122022;
      //   wf->IMRPhenomXHMRingdownPhaseVersion      = 122022;
      //   wf->nCollocPtsRDPhase = 5;
      // }

      
            
  }
  else{
    XLAL_ERROR_VOID(XLAL_EDOM, "Error in IMRPhenomXHM_SetHMWaveformVariables: IMRPhenomXHMReleaseVersion=%i is not valid.\n", wf->IMRPhenomXHMReleaseVersion);
  }
    
  /* Common for all releases */
    
  if(wf->modeTag==32){ 
      wf->nCollocPtsInterPhase=6;
  }
  else{
      wf->nCollocPtsInterPhase=5;
  }


  /* Limit between comparable and extreme mass ratios for the phase */
  wf->etaEMR=0.05;

  /* Spin parameterisations: add here spin pars used for the higher modes that are not included in the 22 struct */
  wf->chi_s = (wf22->chi1L + wf22->chi2L)*0.5;
  wf->chi_a = (wf22->chi1L - wf22->chi2L)*0.5;

  /* Ringdown and damping frequencies*/
  wf->fRING = (qnms->fring_lm[wf->modeInt](wf22->afinal))/wf22->Mfinal;
  wf->fDAMP = (qnms->fdamp_lm[wf->modeInt](wf22->afinal))/wf22->Mfinal;
  
  /* Ringdown and damping frequencies*/
  /* (IMRPhenomXPNRUseTunedCoprec) The EZH effective ringdown perscription uses the precessing final spin while the calibrated XPNR uses the non-precessing final spin to calculate its base waveform to which precession effect are added.  */
  #if DEBUG == 1
    printf("\n>>> wf22->IMRPhenomXPNRUseTunedCoprec=%i ***\n",wf22->IMRPhenomXPNRUseTunedCoprec);
    printf("*** PNR Co-precessing model in use for (l,m)=(%i,%i) ***\n\n",wf->ell,wf->emm);
    printf("wf22->afinal      : %e\n",wf22->afinal);
    printf("wf22->afinal_prec  : %e\n",wf22->afinal_prec);
  #endif
  
  /* We wish to use the EZH formula for the non-tuned HMs. This formula requires the precessing final spin, 
   * while tuning was performed relative to the non-precessing model and its non-precessing final spin 
   * (as would be conferred through wf22->afinal -- see IMRPhenomX_precession.c around about line 566).
   * 
   * Just a note: when IMRPhenomXPNRUseTunedCoprec is false, afinal = afinal_prec and the code just above
   * assigns the appropriate ringdown frequency and damping time. See LALSimIMRPhenomX_precession.c for 
   * the definition of afinal. */
  if( wf22->IMRPhenomXPNRUseTunedCoprec )
  {
    wf->fRING = (qnms->fring_lm[wf->modeInt](wf22->afinal_prec))/wf22->Mfinal;
    wf->fDAMP = (qnms->fdamp_lm[wf->modeInt](wf22->afinal_prec))/wf22->Mfinal;
  #if DEBUG == 1
    printf("\n** ell, emm **     : %i,%i\n",wf->ell,wf->emm);
    printf("fring              : %e\n",wf->fRING);
  #endif      
    wf->fRING = wf->fRING - emm * wf22->fRINGEffShiftDividedByEmm;
  #if DEBUG == 1
    printf("fring shift        : %e\n",- emm * wf22->fRINGEffShiftDividedByEmm);
    printf("fring (coprec)     : %e\n",wf->fRING);
    printf(">> Note that fring has now been set to fringEff <<\n");
  #endif
  }

  if ( wf22->IMRPhenomXPNRUseTunedCoprec33 ) {
  if (wf->modeTag==33) {
    
    #if DEBUG == 1
      printf("fring           : %e\n",wf->fRING);
    #endif
    
    // Apply PNR CoPrec deviations. NOTE that the are OFF when wf22->IMRPhenomXPNRUseTunedCoprec33 is false (see code above)
    wf->fRING = wf->fRING - (wf->PNR_DEV_PARAMETER * wf->NU5);
    wf->fDAMP = wf->fDAMP + (wf->PNR_DEV_PARAMETER * wf->NU6);
    
    #if DEBUG == 1
      printf("fring shift     : %e\n",- (wf->PNR_DEV_PARAMETER * wf->NU5));
      printf("fring (coprec)  : %e\n",wf->fRING);
      printf("fring prec      : %e\n",(qnms->fring_lm[wf->modeInt](wf22->afinal_prec))/wf22->Mfinal);
      printf("fring eff shift : %e\n",- emm * wf22->fRINGEffShiftDividedByEmm);
    #endif
    
    if( wf22->IMRPhenomXPNRUseTunedCoprec ){
      // Note that we transition to the EZH effective ringdown frequency outside of the coprecessing calibration region, ie where PNR_DEV_PARAMETER and pnr_window become zero
      wf->fRING = wf->fRING - (1.0-wf22->pnr_window) * emm * wf22->fRINGEffShiftDividedByEmm;
    } 
    
  }
  }

  /* If (l,m)=(3,2), load the mixing coeffs to transform the spheroidal-harmonic ringdown ansatz back to spherical-harmonic */
  if(wf->modeTag==32)
  {
    IMRPhenomXHM_Initialize_MixingCoeffs(wf,wf22);
  }

  /* Linear part of the 22 phase. timeshift * ff + phaseshift. */
  wf->timeshift =0;//XLALSimIMRPhenomXLinb(wf22->eta, wf22->chiPNHat, wf22->dchi, wf22->delta);
  wf->phaseshift=0;//XLALSimIMRPhenomXLina(wf22->eta, wf22->chiPNHat, wf22->dchi, wf22->delta);
  // current time-alignment of the hybrids
  REAL8 psi4tostrain=XLALSimIMRPhenomXPsi4ToStrain(wf22->eta, wf22->STotR, wf22->dchi);
  wf->DeltaT= -2.*LAL_PI*(500+psi4tostrain);

  wf->fPhaseRDflat = 0.0;
  wf->fAmpRDfalloff = 0.0;

}

/* Wrapper function to return ringdown frequency */
REAL8 IMRPhenomXHM_GenerateRingdownFrequency(
  UINT4 ell,
  UINT4 emm,
  IMRPhenomXWaveformStruct *wf22
  )
{
  /* emm is guaranteed to be positive */
  UINT4 modeTag = ell * 10 + emm;

  REAL8 fRING = 0.0;

  /* if the tuned coprecessing tuning is activated, use the precessing final spin
   * (the final spin is still precessing in the other case when called through XPHM,
   * but the value is assigned to wf22->afinal instead) */
  REAL8 afinal = (wf22->IMRPhenomXPNRUseTunedCoprec) ? wf22->afinal_prec : wf22->afinal;

  switch(modeTag)
  {
    case 21:{
      fRING = (evaluate_QNMfit_fring21(afinal))/wf22->Mfinal;
      break;
    }
    case 22:{
      fRING = wf22->fRING;
      break;
    }
    case 33:{
      fRING = (evaluate_QNMfit_fring33(afinal))/wf22->Mfinal;
      break;
    }
    case 32:{
      fRING = (evaluate_QNMfit_fring32(afinal))/wf22->Mfinal;
      break;
    }
    case 44:{
      fRING = (evaluate_QNMfit_fring44(afinal))/wf22->Mfinal;
      break;
    }
    default:
    {XLAL_ERROR_REAL8(XLAL_EDOM, "Error in IMRPhenomXHM_GenerateRingdownFrequency: mode (%i,%i) selected is not currently available. Modes available are ((2,|2|),(2,|1|),(3,|2|),(3,|3|),(4,|4|)).\n", ell, emm);}
  }

  /* if the coprecessing tuning is activated, return Effective RD frequency */
  if( wf22->IMRPhenomXPNRUseTunedCoprec && ((ell!=2)||(emm!=2)) )
  {
    fRING -= emm * wf22->fRINGEffShiftDividedByEmm;
  }

  return fRING;
}

/* Store function names containing phase coefficient/collocation point fits in pPhase->[Inspiral|Intermediate|Ringdown]PhaseFits */
void IMRPhenomXHM_FillPhaseFitsArray(IMRPhenomXHMPhaseCoefficients *pPhase){

  /* Here we fill some vectors of functions. Each function is the fit for one collocation point/coefficient.
  There are three vectors, one for each region: inspiral, intermediate and ringdown.
  The explicit fits can be found in the files of the corresponding region: IMRPhenomXHM_inspiral.c, _intermediate.c, _ringdown.c
  */

  //21
  pPhase->InspiralPhaseFits[0]=IMRPhenomXHM_Insp_Phase_21_lambda;
  pPhase->InspiralPhaseFits[1]=IMRPhenomXHM_Insp_Phase_33_lambda;
  pPhase->InspiralPhaseFits[2]=IMRPhenomXHM_Insp_Phase_32_lambda;
  pPhase->InspiralPhaseFits[3]=IMRPhenomXHM_Insp_Phase_44_lambda;

  //21
  pPhase->IntermediatePhaseFits[0]=IMRPhenomXHM_Inter_Phase_21_p1;
  pPhase->IntermediatePhaseFits[1]=IMRPhenomXHM_Inter_Phase_21_p2;
  pPhase->IntermediatePhaseFits[2]=IMRPhenomXHM_Inter_Phase_21_p3;
  pPhase->IntermediatePhaseFits[3]=IMRPhenomXHM_Inter_Phase_21_p4;
  pPhase->IntermediatePhaseFits[4]=IMRPhenomXHM_Inter_Phase_21_p5;
  pPhase->IntermediatePhaseFits[5]=IMRPhenomXHM_Inter_Phase_21_p6;

  //33
  pPhase->IntermediatePhaseFits[6]=IMRPhenomXHM_Inter_Phase_33_p1;
  pPhase->IntermediatePhaseFits[7]=IMRPhenomXHM_Inter_Phase_33_p2;
  pPhase->IntermediatePhaseFits[8]=IMRPhenomXHM_Inter_Phase_33_p3;
  pPhase->IntermediatePhaseFits[9]=IMRPhenomXHM_Inter_Phase_33_p4;
  pPhase->IntermediatePhaseFits[10]=IMRPhenomXHM_Inter_Phase_33_p5;
  pPhase->IntermediatePhaseFits[11]=IMRPhenomXHM_Inter_Phase_33_p6;

  //32
  pPhase->IntermediatePhaseFits[12]=IMRPhenomXHM_Inter_Phase_32_p1;
  pPhase->IntermediatePhaseFits[13]=IMRPhenomXHM_Inter_Phase_32_p2;
  pPhase->IntermediatePhaseFits[14]=IMRPhenomXHM_Inter_Phase_32_p3;
  pPhase->IntermediatePhaseFits[15]=IMRPhenomXHM_Inter_Phase_32_p4;
  pPhase->IntermediatePhaseFits[16]=IMRPhenomXHM_Inter_Phase_32_p5;
  pPhase->IntermediatePhaseFits[17]=IMRPhenomXHM_Inter_Phase_32_p6;

  //44
  pPhase->IntermediatePhaseFits[18]=IMRPhenomXHM_Inter_Phase_44_p1;
  pPhase->IntermediatePhaseFits[19]=IMRPhenomXHM_Inter_Phase_44_p2;
  pPhase->IntermediatePhaseFits[20]=IMRPhenomXHM_Inter_Phase_44_p3;
  pPhase->IntermediatePhaseFits[21]=IMRPhenomXHM_Inter_Phase_44_p4;
  pPhase->IntermediatePhaseFits[22]=IMRPhenomXHM_Inter_Phase_44_p5;
  pPhase->IntermediatePhaseFits[23]=IMRPhenomXHM_Inter_Phase_44_p6;

  //32 Spheroidal
  pPhase->RingdownPhaseFits[0]=IMRPhenomXHM_RD_Phase_32_p1;
  pPhase->RingdownPhaseFits[1]=IMRPhenomXHM_RD_Phase_32_p2;
  pPhase->RingdownPhaseFits[2]=IMRPhenomXHM_RD_Phase_32_p3;
  pPhase->RingdownPhaseFits[3]=IMRPhenomXHM_RD_Phase_32_p4;
  pPhase->RingdownPhaseFits[4]=IMRPhenomXHM_RD_Phase_32_p5;

}

/* Store function names containing amplitude coeff./collocation point fits in pAmp->[Inspiral|Intermediate|Ringdown]AmpFits */
void IMRPhenomXHM_FillAmpFitsArray(IMRPhenomXHMAmpCoefficients *pAmp){

  /* Here we fill some vectors of functions. Each function is the fit for one collocation point/coefficient.
  There are three vectors, one for each region: inspiral, intermediate and ringdown.
  The explicit fits can be found in the files of the corresponding region: IMRPhenomXHM_inspiral.c, _intermediate.c, _ringdown.c
  */

  //******Inspiral Fits for collocation points******/

  //21                                                      //Frequency of the collocation point
  pAmp->InspiralAmpFits[0]  = IMRPhenomXHM_Insp_Amp_21_iv1; //fcutInsp
  pAmp->InspiralAmpFits[1]  = IMRPhenomXHM_Insp_Amp_21_iv2; //fcutInsp*0.75
  pAmp->InspiralAmpFits[2]  = IMRPhenomXHM_Insp_Amp_21_iv3; //fcutInsp*0.5

  //33
  pAmp->InspiralAmpFits[3]  = IMRPhenomXHM_Insp_Amp_33_iv1; //fcutInsp
  pAmp->InspiralAmpFits[4]  = IMRPhenomXHM_Insp_Amp_33_iv2; //fcutInsp*0.75
  pAmp->InspiralAmpFits[5]  = IMRPhenomXHM_Insp_Amp_33_iv3; //fcutInsp*0.5

  //32
  pAmp->InspiralAmpFits[6]  = IMRPhenomXHM_Insp_Amp_32_iv1; //fcutInsp
  pAmp->InspiralAmpFits[7]  = IMRPhenomXHM_Insp_Amp_32_iv2; //fcutInsp*0.75
  pAmp->InspiralAmpFits[8]  = IMRPhenomXHM_Insp_Amp_32_iv3; //fcutInsp*0.5

  //44
  pAmp->InspiralAmpFits[9]  = IMRPhenomXHM_Insp_Amp_44_iv1; //fcutInsp
  pAmp->InspiralAmpFits[10] = IMRPhenomXHM_Insp_Amp_44_iv2; //fcutInsp*0.75
  pAmp->InspiralAmpFits[11] = IMRPhenomXHM_Insp_Amp_44_iv3; //fcutInsp*0.5


  /*****Intermediate Fits for EMR collocation points, 2 Intermediate regions*****/

  //21                                                           //Frequency of the collocation point
  pAmp->IntermediateAmpFits[0] = IMRPhenomXHM_Inter_Amp_21_int1; //fcutInsp + (fcutRD-fcutInsp)/3
  pAmp->IntermediateAmpFits[1] = IMRPhenomXHM_Inter_Amp_21_int2; //fcutInsp + 2(fcutRD-fcutInsp)/3

  //33
  pAmp->IntermediateAmpFits[2] = IMRPhenomXHM_Inter_Amp_33_int1; //fcutInsp + (fcutRD-fcutInsp)/3
  pAmp->IntermediateAmpFits[3] = IMRPhenomXHM_Inter_Amp_33_int2; //fcutInsp + 2(fcutRD-fcutInsp)/3

  //32
  pAmp->IntermediateAmpFits[4] = IMRPhenomXHM_Inter_Amp_32_int1; //fcutInsp + (fcutRD-fcutInsp)/3
  pAmp->IntermediateAmpFits[5] = IMRPhenomXHM_Inter_Amp_32_int2; //fcutInsp + 2(fcutRD-fcutInsp)/3

  //44
  pAmp->IntermediateAmpFits[6] = IMRPhenomXHM_Inter_Amp_44_int1; //fcutInsp + (fcutRD-fcutInsp)/3
  pAmp->IntermediateAmpFits[7] = IMRPhenomXHM_Inter_Amp_44_int2; //fcutInsp + 2(fcutRD-fcutInsp)/3

  //21                                                           //fInt1 = fcutInsp + (fcutRD-fcutInsp)/3
  pAmp->IntermediateAmpFits[8] = IMRPhenomXHM_Inter_Amp_21_int0; //fcutInsp + (fInt1 - fcutInsp)/3
  pAmp->IntermediateAmpFits[9] = IMRPhenomXHM_Inter_Amp_21_dint0;//fcutInsp + (fInt1 - fcutInsp)/3

  //33
  pAmp->IntermediateAmpFits[10] = IMRPhenomXHM_Inter_Amp_33_int0; //fcutInsp + (fInt1 - fcutInsp)/3
  pAmp->IntermediateAmpFits[11] = IMRPhenomXHM_Inter_Amp_33_dint0;//fcutInsp + (fInt1 - fcutInsp)/3

  //32
  pAmp->IntermediateAmpFits[12] = IMRPhenomXHM_Inter_Amp_32_int0; //fcutInsp + (fInt1 - fcutInsp)/3
  pAmp->IntermediateAmpFits[13] = IMRPhenomXHM_Inter_Amp_32_dint0;//fcutInsp + (fInt1 - fcutInsp)/3

  //44
  pAmp->IntermediateAmpFits[14] = IMRPhenomXHM_Inter_Amp_44_int0; //fcutInsp + (fInt1 - fcutInsp)/3
  pAmp->IntermediateAmpFits[15] = IMRPhenomXHM_Inter_Amp_44_dint0;//fcutInsp + (fInt1 - fcutInsp)/3

  //21
  pAmp->IntermediateAmpFits[16] = IMRPhenomXHM_Inter_Amp_21_int3;
  pAmp->IntermediateAmpFits[17] = IMRPhenomXHM_Inter_Amp_21_int4;

  //33
  pAmp->IntermediateAmpFits[18] = IMRPhenomXHM_Inter_Amp_33_int3;
  pAmp->IntermediateAmpFits[19] = IMRPhenomXHM_Inter_Amp_33_int4;

  //32
  pAmp->IntermediateAmpFits[20] = IMRPhenomXHM_Inter_Amp_32_int3;
  pAmp->IntermediateAmpFits[21] = IMRPhenomXHM_Inter_Amp_32_int4;

  //44
  pAmp->IntermediateAmpFits[22] = IMRPhenomXHM_Inter_Amp_44_int3;
  pAmp->IntermediateAmpFits[23] = IMRPhenomXHM_Inter_Amp_44_int4;

  /****Ringdown Fits for coefficients*****/

  //21
  pAmp->RingdownAmpFits[0]  = IMRPhenomXHM_RD_Amp_21_alambda;
  pAmp->RingdownAmpFits[1]  = IMRPhenomXHM_RD_Amp_21_lambda;
  pAmp->RingdownAmpFits[2]  = IMRPhenomXHM_RD_Amp_21_sigma;

  //33
  pAmp->RingdownAmpFits[3]  = IMRPhenomXHM_RD_Amp_33_alambda;
  pAmp->RingdownAmpFits[4]  = IMRPhenomXHM_RD_Amp_33_lambda;
  pAmp->RingdownAmpFits[5]  = IMRPhenomXHM_RD_Amp_33_sigma; //currently constant

  //32
  pAmp->RingdownAmpFits[6]  = IMRPhenomXHM_RD_Amp_32_alambda;
  pAmp->RingdownAmpFits[7]  = IMRPhenomXHM_RD_Amp_32_lambda;
  pAmp->RingdownAmpFits[8]  = IMRPhenomXHM_RD_Amp_32_sigma; //currently constant

  //44
  pAmp->RingdownAmpFits[9]  = IMRPhenomXHM_RD_Amp_44_alambda;
  pAmp->RingdownAmpFits[10] = IMRPhenomXHM_RD_Amp_44_lambda;
  pAmp->RingdownAmpFits[11] = IMRPhenomXHM_RD_Amp_44_sigma; //currently constant


  /****Ringdown Fits for Collocation Points*****/

  //21
  pAmp->RingdownAmpFits[12]  = IMRPhenomXHM_RD_Amp_21_rdcp1;
  pAmp->RingdownAmpFits[13]  = IMRPhenomXHM_RD_Amp_21_rdcp2;
  pAmp->RingdownAmpFits[14]  = IMRPhenomXHM_RD_Amp_21_rdcp3;

  //33
  pAmp->RingdownAmpFits[15]  = IMRPhenomXHM_RD_Amp_33_rdcp1;
  pAmp->RingdownAmpFits[16]  = IMRPhenomXHM_RD_Amp_33_rdcp2;
  pAmp->RingdownAmpFits[17]  = IMRPhenomXHM_RD_Amp_33_rdcp3;

  //32
  pAmp->RingdownAmpFits[18]  = IMRPhenomXHM_RD_Amp_32_rdcp1;
  pAmp->RingdownAmpFits[19]  = IMRPhenomXHM_RD_Amp_32_rdcp2;
  pAmp->RingdownAmpFits[20]  = IMRPhenomXHM_RD_Amp_32_rdcp3;

  //44
  pAmp->RingdownAmpFits[21]  = IMRPhenomXHM_RD_Amp_44_rdcp1;
  pAmp->RingdownAmpFits[22]  = IMRPhenomXHM_RD_Amp_44_rdcp2;
  pAmp->RingdownAmpFits[23]  = IMRPhenomXHM_RD_Amp_44_rdcp3;

  //32
  pAmp->RingdownAmpFits[24]  = IMRPhenomXHM_RD_Amp_32_rdaux1;
  pAmp->RingdownAmpFits[25]  = IMRPhenomXHM_RD_Amp_32_rdaux2;

}


/***************************************************/
/*                                                 */
/*         Amplitude Cutting Frequencies           */
/*                                                 */
/***************************************************/

/* Inspiral cutting frequency for the Amplitude */
double IMRPhenomXHM_Amplitude_fcutInsp(IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22){

  //Return the end frequency of the inspiral region and the beginning of the intermediate for the amplitude of one mode.

  int  version  = pWFHM->IMRPhenomXHMInspiralAmpFreqsVersion;
  double fcut = 0.;   //Cutting frequency for comparable mass ratios
  double fMECO  = pWFHM->fMECOlm;
  double emm    = 1.*(pWFHM->emm);
  double eta    = pWF22->eta;
  double chi1   = pWF22->chi1L;
  //fcutEMR is the cutting frequency for extreme mass ratios that is given by a fit to the frequncy of a particular geometrical structure of the amplitude
  double fcutEMR = 1.25*emm*((0.011671068725758493 - 0.0000858396080377194*chi1 + 0.000316707064291237*pow(chi1,2))*(0.8447212540381764 + 6.2873167352395125*eta))/(1.2857082764038923 - 0.9977728883419751*chi1);

  switch(version){
    case 122018: // default version
    {
        
        double fring  = pWFHM->fRING;
        double chieff = pWF22->chiEff;
        double fISCO  = (pWF22->fISCO)*emm*0.5;
        
        switch(pWFHM->modeTag)
        {
            case 21:{
              if(eta < 0.023795359904818562){ //for EMR (q>40.)
                fcut = fcutEMR;
              }
              else{                //for comparable q
                fcut = fMECO + (0.75-0.235*chieff - 5./6.*chieff*chieff)*fabs(fISCO-fMECO);
              }
              break;
            }
            case 33:{
              if(eta < 0.04535147392290249){ //for EMR (q>20.)
                fcut = fcutEMR;
              }
              else{                //for comparable q
                fcut = fMECO + (0.75-0.235*chieff-5./6.*chieff)*fabs(fISCO-fMECO);
              }
              break;
            }
            case 32:{
              if(eta < 0.04535147392290249){ //for extreme mass ratios (q>20)
                fcut = fcutEMR;
              }
              else{               //for comparable mass ratios
                fcut = fMECO + (0.75-0.235*fabs(chieff))*fabs(fISCO-fMECO);
                fcut = fcut*fring/pWF22->fRING;
              }
              break;
            }
            case 44:{
              if(eta <  0.04535147392290249){  //for EMR (q>20)
                fcut = fcutEMR;
              }
              else{                  //for comparable q
                fcut = fMECO + (0.75-0.235*chieff)*fabs(fISCO-fMECO);
              }
              break;
            }
        }
        break;
    }
    case 122022:
    {
        if (pWF22->q < 20.){
          fcut = fMECO;
        }
        else{
          REAL8 transition_eta = 0.0192234; // q=50
          REAL8 sharpness = 0.004;
          REAL8 funcs = 0.5 + 0.5 * tanh((eta - transition_eta)/sharpness);
          fcut = funcs * fMECO + (1 - funcs) * fcutEMR;
        }
        break;
    }
    default: {XLALPrintError("Error in IMRPhenomXHM_Amplitude_fcutInsp: version %i is not valid.", version);}
  }
  return fcut;
}

/* Ringdown cutting frequency for the amplitude */
double IMRPhenomXHM_Amplitude_fcutRD(IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22){

  //Returns the end of the intermediate region and the beginning of the ringdown for the amplitude of one mode

  double fring = pWFHM->fRING, fdamp=pWFHM->fDAMP;
  int  version = pWFHM->IMRPhenomXHMRingdownAmpFreqsVersion;
  double eta   = pWF22->eta;
  double chi1  = pWF22->chi1L;
  double fcut = 0.;  //This is the cutting frequency

  switch(version){
    case 122018: // default version
    {
      switch(pWFHM->modeTag)
      {
        case 21:{
          fcut = 0.75*fring;
          break;
        }
        case 33:{
          fcut = 0.95*fring;
          break;
        }
        case 32:{
          double fRD22 = pWF22->fRING;
          double  c = 0.5, r=5.;
          if(eta < 0.0453515){
            //for extreme mass ratios (q>20)
            fcut = (fring*exp(c*r) + fRD22*exp(r*chi1))/(exp(c*r) + exp(r*chi1)) - fdamp;  //This is a smooth step function between fRING (preferred by negative spins) and fRD22 (preferred by positive)
          }else{
            //for comparable mass ratios
            fcut = fRD22;
          }
          if(0.02126654064272212<eta && eta<0.12244897959183673 && chi1>0.95) // for 6 < q < 45
          {
             fcut = fring - 2.*fdamp;
          }
          break;
        }
        case 44:{
          fcut = 0.9*fring;
          break;
        }
      }
      break;
    }
    case 122022:
    {
        if(pWFHM->MixingOn == 1)
            fcut = pWF22->fRING - 0.5*pWF22->fDAMP; //v8
        else
            fcut = fring - fdamp; //v2
        break;
    }
    default: {XLALPrintError("Error in IMRPhenomXHM_Amplitude_fcutRD: version %i is not valid.", version);}
  }
  return fcut;
}

/***************************************************/
/*                                                 */
/*  Phase Cutting & Collocation Points Frequencies */
/*                                                 */
/***************************************************/

/* GetfcutInsp  */
double GetfcutInsp(IMRPhenomXWaveformStruct *pWF22, IMRPhenomXHMWaveformStruct *pWFHM){
  double MECOf=pWF22->fMECO;
  double eta_factor=1.0+0.001*(0.25/pWF22->eta-1.);

  return (eta_factor*pWFHM->emm*0.5*MECOf);
}


/********************** INTERMEDIATE PHASE COLLOCATION POINTS ******************/
void IMRPhenomXHM_Intermediate_CollocPtsFreqs(IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22){

  //The frequencies are stored in the struct pPhase, which is initialized with new values as the code starts to work on each mode.

  double fring = pWFHM->fRING, fdamp=pWFHM->fDAMP;
  int version  = pWFHM->IMRPhenomXHMIntermediatePhaseFreqsVersion;

    switch(version){
        case 122019: // default version
        case 122022:
        {
            double fcut= GetfcutInsp(pWF22,pWFHM);
            pPhase->CollocationPointsFreqsPhaseInter[0]=fcut;

            if(pWFHM->modeTag==32){


                double fRD22 = pWF22->fRING, fdamp22 = pWF22->fDAMP;
                double fEnd= fRD22- 0.5*fdamp22;
                pPhase->CollocationPointsFreqsPhaseInter[1]=(sqrt(3)*(fcut - fEnd) + 2*(fcut + fEnd))/4.;
                pPhase->CollocationPointsFreqsPhaseInter[2]=(3*fcut + fEnd)/4.;
                pPhase->CollocationPointsFreqsPhaseInter[3]=(fcut + fEnd)/2.;
                // we use first and second derivative at fEnd, so this frequency is duplicated here
                pPhase->CollocationPointsFreqsPhaseInter[4]= fEnd;
                pPhase->CollocationPointsFreqsPhaseInter[5]= fEnd;
                pPhase->fPhaseMatchIM = fEnd;
                // correct cutting frequency for EMR with negative spins
                if(pWF22->eta<0.01&&pWF22->chi1L<0 && version==122019){
                  pPhase->fPhaseMatchIM=pPhase->fPhaseMatchIM*(1.2-0.25*pWF22->chi1L);
                }


            }

            else{

                pPhase->CollocationPointsFreqsPhaseInter[1]=(sqrt(3)*(fcut - fring) + 2*(fcut + fring))/4.;
                pPhase->CollocationPointsFreqsPhaseInter[2]=(3*fcut + fring)/4.;
                pPhase->CollocationPointsFreqsPhaseInter[3]=(fcut + fring)/2.;
                pPhase->CollocationPointsFreqsPhaseInter[4]=(fcut + 3*fring)/4.;
                pPhase->CollocationPointsFreqsPhaseInter[5]=(fcut + 7*fring)/8.;
                pPhase->fPhaseMatchIM = fring-fdamp;

            }

            break;
        }
    default: {XLALPrintError("Error in IMRPhenomXHM_Intermediate_CollocPtsFreqs: version is not valid. Version recommended is 122019.");}
  }

    pPhase->fPhaseMatchIN = pWFHM->fMECOlm;
}

/********************** RINGDOWN PHASE COLLOCATION POINTS ******************/

// this function initializes the frequencies of the collocation points for the spheroidal ringdown reconstruction
void IMRPhenomXHM_Ringdown_CollocPtsFreqs(IMRPhenomXHMPhaseCoefficients *pPhase,IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22){

  //The frequencies are stored in the struct pPhase, which is initialized with new values as the code starts to work on each mode.

  double fringlm = pWFHM->fRING, fdamplm = pWFHM->fDAMP;
  double fring22 = pWF22->fRING;

  switch(pWFHM->IMRPhenomXHMRingdownPhaseFreqsVersion){
      case 122019:{
        pPhase->CollocationPointsFreqsPhaseRD[0] = fring22;
        pPhase->CollocationPointsFreqsPhaseRD[2] = fringlm - 0.5*fdamplm;
        pPhase->CollocationPointsFreqsPhaseRD[1] = fringlm - 1.5*fdamplm;
        pPhase->CollocationPointsFreqsPhaseRD[3] = fringlm + 0.5*fdamplm;
        break;
      }
      case 122022:{
        double fdamp22 = pWF22->fDAMP;
        pPhase->CollocationPointsFreqsPhaseRD[0] = fring22 - fdamp22;
        pPhase->CollocationPointsFreqsPhaseRD[1] = fring22;
        pPhase->CollocationPointsFreqsPhaseRD[2] = (fring22 + fringlm) * 0.5;
        pPhase->CollocationPointsFreqsPhaseRD[3] = fringlm;
        pPhase->CollocationPointsFreqsPhaseRD[4] = fringlm + fdamplm;
        break;
      }
      default:{XLAL_ERROR_VOID(XLAL_EINVAL,"Error in IMRPhenomXHM_Ringdown_CollocPtsFreqs: version %i is not valid.", pWFHM->IMRPhenomXHMRingdownPhaseFreqsVersion);}
  }

}


/*******************************************/
/*                                         */
/*       COMPUTE AMPLITUDE COEFFICIENTS    */
/*                                         */
/*******************************************/

/*** Post-Newtonian Inspiral Ansatz Coefficients ***/
void IMRPhenomXHM_GetPNAmplitudeCoefficients(IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomXHMWaveformStruct *pWFHM,  IMRPhenomXWaveformStruct *pWF22){

  /* Fill pAmp with the coefficients of the power series in frequency of the Forier Domain Post-Newtonian Inspiral Ansatz.
  The 21 mode by default does not used the power series because it breaks down before the end of the inspiral, but it corresponding power series is available here. */

  /* The ansatz in Fourier Domain is built as follows: we multiply the Time-domain Post-Newtonian series up to 3PN by the phasing factor given by the Stationary-Phase-Approximation,
  and the we reexpand in powers of f up to 3PN.
  The only difference with the 21 is that this does not perform the last reexpansion, so it is not a real power series but it is a quantity better behaved. */

  /* The coefficients below correspond to those in eqs E10-E14 in arXiv:2001.10914 */

  //int inspversion = pWFHM->IMRPhenomXHMInspiralAmpFitsVersion;
  double chiA = pWFHM->chi_a;
  double chiS = pWFHM->chi_s;
  double eta  = pWF22->eta, delta = pWF22->delta;
  double PI   = powers_of_lalpiHM.itself;

  const double prefactors[] = {sqrt(2)/3., 0.75*sqrt(5/7.), sqrt(5/7.)/3., 4*sqrt(2)/9*sqrt(5/7.)}; //Global factors of each PN hlm
  pAmp->PNglobalfactor = pow(2./(pWFHM->emm),-7/6.)*prefactors[pWFHM->modeInt]; //This is to compensate that we rescale data with the leading order of the 22
 
  /*switch(inspversion){  FIXMEE
    case 122018: // default version
    {*/
      switch(pWFHM->modeTag)
      {
        case 21:{
          if(pWFHM->useFAmpPN == 1){
            Get21PNAmplitudeCoefficients(pAmp, pWF22);
            pAmp->pnInitial = 0.;
            pAmp->pnOneThird = 0.;
            pAmp->pnTwoThirds = 0.;
            pAmp->pnThreeThirds = 0.;
            pAmp->pnFourThirds = 0.;
            pAmp->pnFiveThirds = 0.;
            pAmp->pnSixThirds = 0.;
          }
          else{
            IMRPhenomX_UsefulPowers powers_of_2d1;
            IMRPhenomX_Initialize_Powers(&powers_of_2d1, 2.);
            pAmp->pnInitial = 0.;
            pAmp->pnOneThird = delta*powers_of_lalpiHM.one_third*powers_of_2d1.one_third;
            pAmp->pnTwoThirds = (-3*(chiA + chiS*delta))/2.*powers_of_lalpiHM.two_thirds*powers_of_2d1.two_thirds;
            pAmp->pnThreeThirds = (335*delta + 1404*delta*eta)/672.*powers_of_lalpiHM.itself*powers_of_2d1.itself;
            pAmp->pnFourThirds = (3427*chiA - (I*672)*delta + 3427*chiS*delta - 8404*chiA*eta - 3860*chiS*delta*eta - 1344*delta*PI - (I*672)*delta*log(16))/1344.*powers_of_lalpiHM.four_thirds*powers_of_2d1.four_thirds;
            pAmp->pnFiveThirds = (-155965824*chiA*chiS - 964357*delta + 432843264*chiA*chiS*eta - 23670792*delta*eta + 24385536*chiA*PI + 24385536*chiS*delta*PI - 77982912*delta*chiA*chiA + 81285120*delta*eta*chiA*chiA - 77982912*delta*chiS*chiS + 39626496*delta*eta*chiS*chiS + 21535920*delta*eta*eta)/8.128512e6*powers_of_lalpiHM.five_thirds*powers_of_2d1.five_thirds;
            pAmp->pnSixThirds = (143063173*chiA - (I*1350720)*delta + 143063173*chiS*delta - 546199608*chiA*eta - (I*72043776)*delta*eta - 169191096*chiS*delta*eta - 9898560*delta*PI + 20176128*delta*eta*PI - (I*5402880)*delta*log(2) - (I*17224704)*delta*eta*log(2) + 61725888*chiS*delta*chiA*chiA - 81285120*chiS*delta*eta*chiA*chiA + 20575296*pow(chiA,3) - 81285120*eta*pow(chiA,3) + 61725888*chiA*chiS*chiS - 165618432*chiA*eta*chiS*chiS + 20575296*delta*pow(chiS,3) - 1016064*delta*eta*chiS*chiS*chiS + 128873808*chiA*eta*eta - 3859632*chiS*delta*eta*eta)/5.419008e6*powers_of_lalpiHM.two*powers_of_2d1.two;
          }
          break;
        }
        case 33:{
          IMRPhenomX_UsefulPowers powers_of_2d3;
          IMRPhenomX_Initialize_Powers(&powers_of_2d3, 2./3.);
          pAmp->pnInitial = 0.;
          pAmp->pnOneThird = delta*powers_of_lalpiHM.one_third*powers_of_2d3.one_third;
          pAmp->pnTwoThirds = 0.;
          pAmp->pnThreeThirds = (-1945*delta + 2268*delta*eta)/672.*powers_of_lalpiHM.itself*powers_of_2d3.itself;
          pAmp->pnFourThirds = (325*chiA - (I*504)*delta + 325*chiS*delta - 1120*chiA*eta - 80*chiS*delta*eta + 120*delta*PI + (I*720)*delta*log(1.5))/120.*powers_of_lalpiHM.four_thirds*powers_of_2d3.four_thirds;
          pAmp->pnFiveThirds = (-2263282560*chiA*chiS - 1077664867*delta + 9053130240*chiA*chiS*eta - 5926068792*delta*eta - 1131641280*delta*chiA*chiA + 4470681600*delta*eta*chiA*chiA - 1131641280*delta*chiS*chiS + 55883520*delta*eta*chiS*chiS + 2966264784*delta*eta*eta)/4.4706816e8*powers_of_lalpiHM.five_thirds*powers_of_2d3.five_thirds;
          pAmp->pnSixThirds = (22007835*chiA + (I*26467560)*delta + 22007835*chiS*delta - 80190540*chiA*eta - (I*98774368)*delta*eta - 31722300*chiS*delta*eta - 9193500*delta*PI + 17826480*delta*eta*PI - (I*37810800)*delta*log(1.5) + (I*37558080)*delta*eta*log(1.5) - 12428640*chiA*eta*eta - 6078240*chiS*delta*eta*eta)/2.17728e6*powers_of_lalpiHM.two*powers_of_2d3.two;
          break;
        }
        case 32:{
          pAmp->pnInitial = 0.;
          pAmp->pnOneThird = 0.;
          pAmp->pnTwoThirds = (-1 + 3*eta)*powers_of_lalpiHM.two_thirds;
          pAmp->pnThreeThirds = -4*chiS*eta*powers_of_lalpiHM.itself;
          pAmp->pnFourThirds = (10471 - 61625*eta + 82460*eta*eta)/10080.*powers_of_lalpiHM.four_thirds;
          pAmp->pnFiveThirds = ((I*2520) - 3955*chiS - 3955*chiA*delta - (I*11088)*eta + 10810*chiS*eta + 11865*chiA*delta*eta - 12600*chiS*eta*eta)/840.*powers_of_lalpiHM.five_thirds;
          pAmp->pnSixThirds = (824173699 + 2263282560*chiA*chiS*delta - 26069649*eta - 15209631360*chiA*chiS*delta*eta + 3576545280*chiS*eta*PI + 1131641280*chiA*chiA - 7865605440*eta*chiA*chiA + 1131641280*chiS*chiS - 11870591040*eta*chiS*chiS - 13202119896*eta*eta + 13412044800*chiA*chiA*eta*eta + 5830513920*chiS*chiS*eta*eta + 5907445488*pow(eta,3))/4.4706816e8*powers_of_lalpiHM.two;
          break;
        }
        case 44:{
          IMRPhenomX_UsefulPowers powers_of_2d4;
          IMRPhenomX_Initialize_Powers(&powers_of_2d4, 0.5);
          pAmp->pnInitial = 0.;
          pAmp->pnOneThird = 0.;
          pAmp->pnTwoThirds = (1 - 3*eta)*powers_of_lalpiHM.two_thirds*powers_of_2d4.two_thirds;
          pAmp->pnThreeThirds = 0.;
          pAmp->pnFourThirds = (-158383 + 641105*eta - 446460*eta*eta)/36960.*powers_of_lalpiHM.four_thirds*powers_of_2d4.four_thirds;
          pAmp->pnFiveThirds = ((I*-1008) + 565*chiS + 565*chiA*delta + (I*3579)*eta - 2075*chiS*eta - 1695*chiA*delta*eta + 240*PI - 720*eta*PI + (I*960)*log(2) - (I*2880)*eta*log(2) + 1140*chiS*eta*eta)/120.*powers_of_lalpiHM.five_thirds*powers_of_2d4.five_thirds;
          pAmp->pnSixThirds = (7888301437 - 147113366400*chiA*chiS*delta - 745140957231*eta + 441340099200*chiA*chiS*delta*eta - 73556683200*chiA*chiA + 511264353600*eta*chiA*chiA - 73556683200*chiS*chiS + 224302478400*eta*chiS*chiS + 2271682065240*eta*eta - 871782912000*chiA*chiA*eta*eta - 10897286400*chiS*chiS*eta*eta - 805075876080*pow(eta,3))/2.90594304e10*powers_of_lalpiHM.two*powers_of_2d4.two;
          break;
        }
    }/*
      break;
    }
    default: {XLALPrintError("Error in IMRPhenomXHM_GetPNAmplitudeCoefficients: version is not valid. Version recommended is 122018. ");}
}*/
}

/*** Post-Newtonian Inspiral Ansatz Coefficients for the 21 mode ***/
void Get21PNAmplitudeCoefficients(IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomXWaveformStruct* pWF22){

  /* The 21 ansatz in Fourier Domain is built multiplying the Time-domain Post-Newtonian series up to 3PN by the phasing factor given by the Stationary-Phase-Approximatio */

  double m1  = pWF22->m1;
  double m2  = pWF22->m2;
  double m12 = m1*m1;
  double m22 = m2*m2;
  double m13 = m12*m1;
  double m23 = m22*m2;
  double m14 = m13*m1;
  double m24 = m23*m2;
  double m15 = m14*m1;
  double m25 = m24*m2;
  double m16 = m15*m1;

  double chi1   = pWF22->chi1L;
  double chi2   = pWF22->chi2L;
  double chiS   = (chi1 + chi2)*0.5;
  double chiA   = (chi1 - chi2)*0.5;
  double delta  = pWF22->delta;
  double eta    = pWF22->eta;
  double chi12  = chi1*chi1;
  double chi22  = chi2*chi2;
  double Sc     = m12*chi1 + m22*chi2;
  double Sigmac = m2*chi2 - m1*chi1;

  IMRPhenomX_UsefulPowers powers_of_2d1;
  IMRPhenomX_Initialize_Powers(&powers_of_2d1, 2.);
  double logof2 = powers_of_2d1.log;
  double log4   = 1.3862943611198906;
  double EulerGamma = 0.5772156649015329;

  /* Complex coefficients of the Time-Domain Post-Newtonian expansion. */

  double factor = 8*pWF22->eta*powers_of_lalpiHM.two_thirds*powers_of_2d1.two_thirds*powers_of_lalpiHM.sqrt/sqrt(5.);
  pAmp->PNTDfactor = factor;
  pAmp->x05 = I*delta/3*powers_of_2d1.one_third*powers_of_lalpiHM.one_third;
  pAmp->x1  = -I*0.5*(chiA + chiS*delta)*powers_of_2d1.two_thirds*powers_of_lalpiHM.two_thirds;
  pAmp->x15 = I*delta*(-17./28 + 5*eta/7.)/3.*powers_of_2d1.itself*powers_of_lalpiHM.itself;
  pAmp->x2  = ( I*(-43./21*delta*Sc + (-79.+139*eta)/42.*Sigmac) + I/3*delta*(powers_of_lalpiHM.itself + I*(-0.5-2*logof2)))*powers_of_2d1.four_thirds*powers_of_lalpiHM.four_thirds;
  pAmp->x25 = I*delta*((-43-509*eta)/126. + 79*eta*eta/168.)/3.*powers_of_2d1.five_thirds*powers_of_lalpiHM.five_thirds;
  pAmp->x3  = I*delta*( (-17. + 6.*eta)/28.*powers_of_lalpiHM.itself + I*(17/56. + eta*(-353/28. - 3.*logof2/7.) + 17.*logof2/14.))/3.*powers_of_2d1.two*powers_of_lalpiHM.two;


  /* Coefficients of the phasing factor expansion */
  /* These coefficients correspond to those in equation E4 of arXiv:2001.10914. */
  pAmp->xdot5 = -(m1*m2*(-838252800*m1*m2 - 419126400*m12 - 419126400*m22)/3.274425e7)*powers_of_2d1.five_thirds*powers_of_2d1.five_thirds*powers_of_lalpiHM.five_thirds*powers_of_lalpiHM.five_thirds;
  pAmp->xdot6 = -(m1*m2*(1152597600*m2*m13 + 926818200*m22 + 2494800*m1*m2*(743 + 462*m22) + 1247400*m12*(743 + 1848*m22))/3.274425e7)*powers_of_2d1.four*powers_of_lalpiHM.four;
  pAmp->xdot65 = -(m1*m2*(-34927200*m1*m2*(-(m2*(75*chi1 + 376*chi2*m2)) + 96*powers_of_lalpiHM.itself) - 34927200*(-(m2*(75*chi2 + 188*(chi1 + chi2)*m2)) + 48*powers_of_lalpiHM.itself)*m12 - 2619540000*chi1*m13 + 13132627200*chi1*m2*m13 +
  6566313600*chi1*m14 - 34927200*(chi2*(75 - 188*m2)*m2 + 48*powers_of_lalpiHM.itself)*m22)/3.274425e7)*powers_of_2d1.eight_thirds*powers_of_2d1.five_thirds*powers_of_lalpiHM.eight_thirds*powers_of_lalpiHM.five_thirds;
  pAmp->xdot7 = -(m1*m2*(207900*m2*(-13661 - 19908*chi1*chi2 + 10206*chi12 + 10206*chi22)*m13 - 23100*(34103 + 91854*chi22)*m22 - 1373803200*m14*m22 +
  23100*m1*m2*(-2*(34103 + 45927*chi12 + 45927*chi22) + 9*(-13661 - 19908*chi1*chi2 + 10206*chi12 + 10206*chi22)*m22) - 2747606400*m13*m23 -
  23100*m12*(34103 + 91854*chi12 - 18*(-13661 - 19908*chi1*chi2 + 10206*chi12 + 10206*chi22)*m22 + 59472*m24))/3.274425e7)*powers_of_2d1.seven_thirds*powers_of_2d1.seven_thirds*powers_of_lalpiHM.seven_thirds*powers_of_lalpiHM.seven_thirds;
  pAmp->xdot75 = -(m1*m2*(-4036586400*chi1*m13 + 5821200*m2*(5861*chi1 + 1701*powers_of_lalpiHM.itself)*m13 + 17059026600*chi1*m14 + 14721814800*chi1*m2*m14 - 34962127200*chi1*m2*m15 +
  207900*(2*chi2*m2*(-9708 + 41027*m2) + 12477*powers_of_lalpiHM.itself)*m22 - 14721814800*chi2*m13*m22 - 69924254400*chi1*m14*m22 - 34962127200*(chi1 + chi2)*m13*m23 +
  207900*m12*(3*powers_of_lalpiHM.itself*(4159 + 31752*m22) + 2*m2*(9708*chi2 + 41027*(chi1 + chi2)*m2 - 35406*chi1*m22 - 168168*chi2*m23)) +
  415800*m1*m2*(9708*chi1*m2 + 82054*chi2*m22 + 3*powers_of_lalpiHM.itself*(4159 + 7938*m22) + 35406*chi2*m23 - 84084*chi2*m24))/3.274425e7)*powers_of_2d1.eight_thirds*powers_of_2d1.seven_thirds*powers_of_lalpiHM.eight_thirds*powers_of_lalpiHM.seven_thirds;
  pAmp->xdot8 = -(m1*m2*(- 10548014400*chi1*powers_of_lalpiHM.itself*m13 - 63392868000*chi1*chi2*m2*m14 + 34927200*chi1*(-375*chi1 + 752*powers_of_lalpiHM.itself)*m14 + 63392868000*chi12*m15-
  153213984000*m2*chi12*m15 - 76606992000*chi12*m16 - 63392868000*chi1*(chi1 - chi2)*m13*m22 -
  51975*(4869 + 2711352*chi1*chi2 + 1702428*chi12 + 228508*chi22)*m14*m22 - 103950*(4869 + 2711352*chi1*chi2 + 228508*chi12 + 228508*chi22)*m13*m23 +
  906328500*m15*m23 + 1812657000*m14*m24 + 906328500*m13*m25 +
  1925*m2*m13*(56198689 + 13635864*chi1*chi2 + 27288576*chi1*powers_of_lalpiHM.itself + 30746952*chi12 + 3617892*chi22 - 2045736*powers_of_lalpiHM.two) -
  3*m22*(16447322263 - 2277918720*EulerGamma - 23284800*chi2*m2*(-151 + 376*m2)*powers_of_lalpiHM.itself - 2277918720*log4 + 2321480700*chi22 + 4365900000*chi22*m22 -
  21130956000*chi22*m23 + 25535664000*chi22*m24 + 745113600*powers_of_lalpiHM.two) +
  m12*(6833756160*EulerGamma + 10548014400*chi2*m2*powers_of_lalpiHM.itself + 63392868000*(chi1 - chi2)*chi2*m23 - 51975*(4869 + 2711352*chi1*chi2 + 228508*chi12 + 1702428*chi22)*m24 -
  3850*m22*(-56198689 + 13580136*chi1*chi2 - 6822144*(chi1 + chi2)*powers_of_lalpiHM.itself - 6976422*chi12 - 6976422*chi22 + 2045736*powers_of_lalpiHM.two) -
  3*(16447322263 - 2277918720*log4 + 2321480700*chi12 + 745113600*powers_of_lalpiHM.two)) +
  m1*m2*(13667512320*EulerGamma + 10548014400*chi1*m2*powers_of_lalpiHM.itself - 63392868000*chi1*chi2*m23 - 153213984000*chi22*m24 -
    1925*m22*(-56198689 - 13635864*chi1*chi2 - 27288576*chi2*powers_of_lalpiHM.itself - 3617892*chi12 - 30746952*chi22 + 2045736*powers_of_lalpiHM.two) -
    6*(16447322263 - 2277918720*log4 + 1160740350*chi12 + 1160740350*chi22 + 745113600*powers_of_lalpiHM.two)))/3.274425e7)*powers_of_2d1.eight_thirds*powers_of_2d1.eight_thirds*powers_of_lalpiHM.eight_thirds*powers_of_lalpiHM.eight_thirds;
    pAmp->xdot8Log = -(m1*m2*3416878080)/3.274425e7*powers_of_2d1.eight_thirds*powers_of_2d1.eight_thirds*powers_of_lalpiHM.eight_thirds*powers_of_lalpiHM.eight_thirds;
    pAmp->xdot85 = -(m1*m2*(-14891068500*chi1*m13 + 1925*m2*(97151928*chi1 + 6613488*chi2 - 12912300*powers_of_lalpiHM.itself)*m13 + 87143248500*chi1*m14 + 33313480200*chi1*m2*m14 - 198816225300*chi1*m2*m15 +
    57750*(2*chi2*m2*(-128927 + 754487*m2) + 7947*powers_of_lalpiHM.itself)*m22 - 33313480200*chi2*m13*m22 - 138600*(3399633*chi1 + 530712*chi2 + 182990*powers_of_lalpiHM.itself)*m14*m22 -
    35665037100*chi1*m15*m22 + 84184254000*chi1*m16*m22 -
    23100*m1*m2*(15*powers_of_lalpiHM.itself*(-2649 + 71735*m22) + m2*(-(chi1*(644635 + 551124*m2)) + chi2*m2*(-8095994 - 1442142*m2 + 8606763*m22))) -
    69300*(4991769*(chi1 + chi2) + 731960*powers_of_lalpiHM.itself)*m13*m23 + 35665037100*chi2*m14*m23 + 170726094000*chi1*m15*m23 + 2357586000*chi2*m15*m23 +
    35665037100*chi1*m13*m24 + 88899426000*(chi1 + chi2)*m14*m24 + 9702000*(243*chi1 + 17597*chi2)*m13*m25 -
    11550*m12*(15*powers_of_lalpiHM.itself*(-2649 + 286940*m22 + 146392*m24) + 2*m2*
    (-644635*chi2 - 4874683*(chi1 + chi2)*m2 + 1442142*chi1*m22 + 54*(58968*chi1 + 377737*chi2)*m23 + 1543941*chi2*m24 - 3644340*chi2*m25)))/3.274425e7)*powers_of_2d1.eight_thirds*powers_of_2d1.eight_thirds*powers_of_2d1.one_third*powers_of_lalpiHM.eight_thirds*powers_of_lalpiHM.eight_thirds*powers_of_lalpiHM.one_third;

    pAmp->log2pi_two_thirds = log(powers_of_lalpiHM.two_thirds * powers_of_2d1.two_thirds);

  }

/*** PHENOM AMPLITUDE COEFFICIENTS ***/
void IMRPhenomXHM_GetAmplitudeCoefficients(IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXAmpCoefficients *pAmp22, IMRPhenomXPhaseCoefficients *pPhase22, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22) {

    /* Take all the phenom coefficients accross the three regions (inspiral, intermeidate and ringdown) and all the needed parameters to reconstruct the amplitude (including mode-mixing). */
    pAmp->ampNorm = pWF22->ampNorm;
    pAmp->PNdominant = pWF22->ampNorm * pow(2./pWFHM->emm, -7/6.); // = Pi * Sqrt(2 eta/3) (2Pi /m)^(-7/6). Miss the f^(-7/6). The pi power included in ampNorm
    switch(pWFHM->modeTag){
        case 21:{
            if(pWF22->q == 1){
                pAmp->PNdominantlmpower = 2;
                pAmp->PNdominantlm = sqrt(2.) / 3. * 1.5 * pWF22->dchi * 0.5 * pow(2 * LAL_PI / pWFHM->emm, 2/3.);
            }
            else{
                pAmp->PNdominantlmpower = 1;
                pAmp->PNdominantlm = sqrt(2.) / 3. * pWF22->delta * pow(2 * LAL_PI / pWFHM->emm, 1/3.);
            }
            break;
        }
        case 33:{
            if(pWF22->q == 1){
                pAmp->PNdominantlmpower = 4;
                pAmp->PNdominantlm = 0.75 * sqrt(5./7) * pWF22->dchi * 0.5 * (65/24. - 28/3. * pWF22->eta) * pow(2 * LAL_PI / pWFHM->emm, 4/3.);
            }
            else{
                pAmp->PNdominantlmpower = 1;
                pAmp->PNdominantlm = 0.75 * sqrt(5./7) * pWF22->delta * pow(2 * LAL_PI / pWFHM->emm, 1/3.);
            }
            break;
        }
        case 32:{
            pAmp->PNdominantlmpower = 2;
            pAmp->PNdominantlm = 0.75 * sqrt(5./7) * pow(2 * LAL_PI / pWFHM->emm, 1/3.);
            break;
        }
        case 44:{
            pAmp->PNdominantlmpower = 2;
            pAmp->PNdominantlm = 4/9. * sqrt(10/7.) * (1 - 3 * pWF22->eta) * pow(2 * LAL_PI / pWFHM->emm, 1/3.);
            break;
        }
        default:
          {XLALPrintError("Error in IMRPhenomXHM_GetAmplitudeCoefficients: mode selected is not currently available. Modes available are ((2,|2|),(2,|1|),(3,|2|),(3,|3|),(4,|4|)).\n");}
    }
    pAmp->PNdominantlm = fabs(pAmp->PNdominantlm);
    pWFHM->fAmpRDfalloff = 0.;
    pAmp->nCoefficientsInter = 0;

    /*** Proceed region by region ***/
    if(pWFHM->IMRPhenomXHMReleaseVersion != 122019){ 
        pAmp->InspRescaleFactor = 0;
        pAmp->RDRescaleFactor = 0;
        pAmp->InterRescaleFactor = 0;
                    
            
        /* Transform IMRPhenomXHMIntermediateAmpVersion number to int array defining what to do for each collocation point */
        /* 0: don't use coll point, 1: use point, 2: use point and derivative (this only applies for boundaries) */
        // e.g. pAmp->VersionCollocPtsInter = {1, 1, 1, 1, 0, 2} -> use the two boundaries, add derivative to the right one, skip third collocation point
        UINT4 num = pWFHM->IMRPhenomXHMIntermediateAmpVersion;
        for(UINT2 i = 0; i < pWFHM->nCollocPtsInterAmp; i++){
            pAmp->VersionCollocPtsInter[pWFHM->nCollocPtsInterAmp - i - 1] = num % 10;
            num/=10;
        }

        pAmp->nCoefficientsInter = 0;
        for (UINT2 i = 0; i < pWFHM->nCollocPtsInterAmp; i++) pAmp->nCoefficientsInter += pAmp->VersionCollocPtsInter[i];
        /* The number of coefficients in the intermediate ansatz cannot be larger than the number of available collocation points in IMRPhenomXHMIntermediateAmpVersion 
           If e.g. IMRPhenomXHMIntermediateAmpVersion has 6 slots, the maximum number of coefficients would be 6 + 2 because we count for the derivatives at the boundaries. */
        if (pAmp->nCoefficientsInter > pWFHM->nCollocPtsInterAmp + 2) 
            XLAL_ERROR_VOID(XLAL_EFUNC, "IMRPhenomXHM_GetAmplitudeCoefficients failed. Inconsistent number of collocation points (%i) and free parameters (%i).", pWFHM->nCollocPtsInterAmp + 2, pAmp->nCoefficientsInter);

        pAmp->nCoefficientsRDAux = 0;
        if (pWFHM->MixingOn){
            pAmp->nCollocPtsRDAux = 2;
            pAmp->nCoefficientsRDAux = 4;
            pAmp->fRDAux = pWFHM->fRING - pWFHM->fDAMP; //v2
        }

        // Take the cutting frequencies at the inspiral and ringdown
        pAmp->fAmpMatchIN  = IMRPhenomXHM_Amplitude_fcutInsp(pWFHM, pWF22);
        pAmp->fAmpMatchIM  = IMRPhenomXHM_Amplitude_fcutRD(pWFHM, pWF22);
        pWFHM->fAmpRDfalloff = pWFHM->fRING + 2 * pWFHM->fDAMP;

        /* Take Frequency Domain Post-Newtonian Amplitude Coefficients */
        IMRPhenomXHM_GetPNAmplitudeCoefficients(pAmp, pWFHM, pWF22);

        IMRPhenomXHM_Get_Inspiral_Amp_Coefficients(pAmp, pWFHM, pWF22);

        IMRPhenomXHM_RD_Amp_Coefficients(pWF22, pWFHM, pAmp);
                
        IMRPhenomXHM_Intermediate_Amp_Coefficients(pAmp, pWFHM, pWF22, pPhase, pAmp22, pPhase22);

        // printf("\nInsp Coll points\n");
        // for(UINT2 i = 0; i < 3; i++){
        //     printf("%.16f %.16e\n", pAmp->CollocationPointsFreqsAmplitudeInsp[i], pAmp->CollocationPointsValuesAmplitudeInsp[i]);
        // }
        // printf("\nInter Coll points\n");
        // for(UINT2 i = 0; i < pAmp->nCoefficientsInter; i++){
        //     printf("%.16f %.16e\n", pAmp->CollocationPointsFreqsAmplitudeInter[i], pAmp->CollocationPointsValuesAmplitudeInter[i]);
        // }
        // printf("\nRD Coll points\n");
        // for(UINT2 i = 0; i < 3; i++){
        //     printf("%.16f %.16e\n", pAmp->CollocationPointsFreqsAmplitudeRD[i], pAmp->CollocationPointsValuesAmplitudeRD[i]);
        // }
        // 
        // printf("\nInsp Coefficients\n");
        // for(UINT2 i = 0; i < 3; i++){
        //     printf("%.16e\n", pAmp->InspiralCoefficient[i]);
        // }
        // printf("\nInter Coefficients %i\n", pAmp->nCoefficientsInter);
        // for(UINT2 i = 0; i < pAmp->nCoefficientsInter; i++){
        //     printf("%.16e\n", pAmp->InterCoefficient[i]);
        // }
        // printf("\nRD Coefficients\n");
        // for(UINT2 i = 0; i < 5; i++){
        //     printf("%.16e\n", pAmp->RDCoefficient[i]);
        // }
        // printf("\nRDAux Coll points\n");
        // for(UINT2 i = 0; i < pAmp->nCoefficientsRDAux; i++){
        //     printf("%.16f %.16e\n", pAmp->CollocationPointsFreqsAmplitudeRDAux[i], pAmp->CollocationPointsValuesAmplitudeRDAux[i]);
        // }
        // printf("\nRDAux Coefficients\n");
        // for(UINT2 i = 0; i < pAmp->nCoefficientsRDAux; i++){
        //     printf("%.16e\n", pAmp->RDAuxCoefficient[i]);
        // }
        

        /* Set Rescale Factors to build Strain (=0) */
        pAmp->InspRescaleFactor = 0;
        pAmp->RDRescaleFactor = 0;
        pAmp->InterRescaleFactor = 0;
    }
    else{
        pAmp->InspRescaleFactor = 1;
        pAmp->InterRescaleFactor = -1;
        pAmp->RDRescaleFactor = 1;

        // Options for the extrapolation of the model outside the calibration region
        if(((pWFHM->modeTag==44) || (pWFHM->modeTag==33)) && pWF22->q>7. && pWF22->chi1L > 0.95){
          pAmp->useInspAnsatzRingdown = 1;
        }
        else{
          pAmp->useInspAnsatzRingdown = 0;
        }
        pAmp->WavyInsp = 0;
        pAmp->WavyInt  = 0;
        if(pWFHM->modeTag==21){
          pAmp->WavyInsp = 1;
          pAmp->WavyInt  = 1;
        }
        if(pWFHM->modeTag==32){
          pAmp->WavyInsp = 1;
        }

        // Take the cutting frequencies at the inspiral and ringdown
        pAmp->fAmpMatchIN  = IMRPhenomXHM_Amplitude_fcutInsp(pWFHM, pWF22);
        pAmp->fAmpMatchIM  = IMRPhenomXHM_Amplitude_fcutRD(pWFHM, pWF22);

        #if DEBUG == 1
        printf("\n\n*** IMRPhenomXHM_GetAmplitudeCoefficients ***\n\n");
        printf("fring_%i = %.16f\n", pWFHM->modeTag, pWFHM->fRING);
        printf("fdamp_%i = %.16f\n", pWFHM->modeTag, pWFHM->fDAMP);
        printf("fcutInsp_%i = %.16f \n", pWFHM->modeTag, pAmp->fAmpMatchIN);
        printf("fcutRD_%i   = %.16f \n", pWFHM->modeTag, pAmp->fAmpMatchIM);
        #endif

        /* Compute the frequencies for Intermediate collocation points. */
        /* It is done now because with the inspiral veto fAmpMatchIN will change, and we need the original here. */
        double df = pAmp->fAmpMatchIM - pAmp->fAmpMatchIN;
        pAmp->CollocationPointsFreqsAmplitudeInter[0] =  pAmp->fAmpMatchIN + df/3. ;
        pAmp->CollocationPointsFreqsAmplitudeInter[1] =  pAmp->fAmpMatchIN + df*2./3.;

        int nCollocPtsInspAmp  = pWFHM->nCollocPtsInspAmp;
        int nCollocPtsInterAmp = pWFHM->nCollocPtsInterAmp;

        int modeint = pWFHM->modeInt;

        #if DEBUG == 1
        printf("nCollocPtsInspAmp  = %i \n",nCollocPtsInspAmp);
        printf("nCollocPtsInterAmp = %i \n",nCollocPtsInterAmp);
        #endif

        /*****************/
        /*    INSPIRAL   */
        /*****************/

        #if DEBUG == 1
        printf("\n**** INSPIRAL ****\n\n");
        printf("IMRPhenomXHMInspiralAmpVersion = %i\r\n",pWFHM->IMRPhenomXHMInspiralAmpVersion);
        #endif

        /* Take Frequency Domain Post-Newtonian Amplitude Coefficients */
        IMRPhenomXHM_GetPNAmplitudeCoefficients(pAmp, pWFHM, pWF22);

        #if DEBUG == 1
        printf("PN coeff %.16f   %.16f\n", creal(pAmp->pnInitial),cimag(pAmp->pnInitial));
        printf("PN coeff %.16f   %.16f\n", creal(pAmp->pnOneThird),cimag(pAmp->pnOneThird));
        printf("PN coeff %.16f   %.16f\n", creal(pAmp->pnTwoThirds),cimag(pAmp->pnTwoThirds));
        printf("PN coeff %.16f   %.16f\n", creal(pAmp->pnThreeThirds),cimag(pAmp->pnThreeThirds));
        printf("PN coeff %.16f   %.16f\n", creal(pAmp->pnFourThirds),cimag(pAmp->pnFourThirds));
        printf("PN coeff %.16f   %.16f\n", creal(pAmp->pnFiveThirds),cimag(pAmp->pnFiveThirds));
        printf("PN coeff %.16f   %.16f\n", creal(pAmp->pnSixThirds),cimag(pAmp->pnSixThirds));
        #endif

        // Initialize frequencies of colloc points   !!!! ASSUMING YOU HAVE 3 COLLOC POINTS
        pAmp->CollocationPointsFreqsAmplitudeInsp[0] = 1.0  * pAmp->fAmpMatchIN;
        pAmp->CollocationPointsFreqsAmplitudeInsp[1] = 0.75 * pAmp->fAmpMatchIN;
        pAmp->CollocationPointsFreqsAmplitudeInsp[2] = 0.5  * pAmp->fAmpMatchIN;

        double fcutInsp, f1, f2, f3;
        fcutInsp = pAmp->fAmpMatchIN;                      // Matching frequency between inspiral and intermediate
        f1 = pAmp->CollocationPointsFreqsAmplitudeInsp[0]; // Frequency of colloc point 1 = 1.*fcutInsp
        f2 = pAmp->CollocationPointsFreqsAmplitudeInsp[1]; // Frequency of colloc point 2 = 0.75*fcutInsp
        f3 = pAmp->CollocationPointsFreqsAmplitudeInsp[2]; // Frequency of colloc point 3 = 0.5*fcutInsp

        // Compute the useful powers of fcutInsp, f1, f2, f3. Remembers: f3 < f2 < f1 = fcutInsp.
        IMRPhenomX_UsefulPowers powers_of_fcutInsp, powers_of_f1, powers_of_f2, powers_of_f3;
        IMRPhenomX_Initialize_Powers(&powers_of_fcutInsp, fcutInsp);  // fcutInsp and f1 are the same but we keep them separated if in the future we change the frequencies of the collocatio points.
        IMRPhenomX_Initialize_Powers(&powers_of_f1, f1);
        IMRPhenomX_Initialize_Powers(&powers_of_f2, f2);
        IMRPhenomX_Initialize_Powers(&powers_of_f3, f3);

        // Compute the values of Post-Newtoninan ansatz (without the pseudo-PN terms) at the frequencies of the collocation points.
        double PNf1, PNf2, PNf3;
        PNf1 = IMRPhenomXHM_Inspiral_PNAmp_Ansatz(&powers_of_f1, pWFHM, pAmp);
        PNf2 = IMRPhenomXHM_Inspiral_PNAmp_Ansatz(&powers_of_f2, pWFHM, pAmp);
        PNf3 = IMRPhenomXHM_Inspiral_PNAmp_Ansatz(&powers_of_f3, pWFHM, pAmp);

        pAmp->PNAmplitudeInsp[0] = PNf1;
        pAmp->PNAmplitudeInsp[1] = PNf2;
        pAmp->PNAmplitudeInsp[2] = PNf3;

        // Initialize values of collocation points at the previous 3 frequencies. They are taken from the parameter space fits.
        for(int i = 0; i<nCollocPtsInspAmp; i++)
        {
          pAmp->CollocationPointsValuesAmplitudeInsp[i] = fabs(pAmp->InspiralAmpFits[modeint*nCollocPtsInspAmp+i](pWF22,pWFHM->IMRPhenomXHMInspiralAmpFitsVersion));
        }
    


        // Values of the collocation point minus the Post-Newtonian value. This gives a "collocation point" for the pseudo-PN part.
        // This way is more convenient becuase the reconstuction does not depended on the PN ansatz used.
        REAL8 iv1, iv2, iv3;
        iv1 = pAmp->CollocationPointsValuesAmplitudeInsp[0] - PNf1;
        iv2 = pAmp->CollocationPointsValuesAmplitudeInsp[1] - PNf2;
        iv3 = pAmp->CollocationPointsValuesAmplitudeInsp[2] - PNf3;


        #if DEBUG == 1
        printf("\nAmplitude pseudo collocation points before veto\n");
        printf("fAmpMatchIN = %.16f\n",pAmp->fAmpMatchIN);
        printf("F1   = %.16f\n", f1);
        printf("F2   = %.16f\n", f2);
        printf("F3   = %.16f\n\n", f3);
        printf("I1   = %.16f\n", pAmp->CollocationPointsValuesAmplitudeInsp[0]);
        printf("I2   = %.16f\n", pAmp->CollocationPointsValuesAmplitudeInsp[1]);
        printf("I3   = %.16f\n\n", pAmp->CollocationPointsValuesAmplitudeInsp[2]);
        printf("PNf1 = %.16f\n", PNf1);
        printf("PNf2 = %.16f\n", PNf2);
        printf("PNf3 = %.16f\n\n", PNf3);
        REAL8 piv1, piv2, piv3;
        piv1 = pAmp->CollocationPointsValuesAmplitudeInsp[0]*pWF22->ampNorm*powers_of_f1.m_seven_sixths;
        piv2 = pAmp->CollocationPointsValuesAmplitudeInsp[1]*pWF22->ampNorm*powers_of_f2.m_seven_sixths;
        piv3 = pAmp->CollocationPointsValuesAmplitudeInsp[2]*pWF22->ampNorm*powers_of_f3.m_seven_sixths;
        printf("p1   = %.16f\n", piv1);
        printf("p2   = %.16f\n", piv2);
        printf("p3   = %.16f\n\n", piv3);
        printf("V1   = %.16f\n", iv1);
        printf("V2   = %.16f\n", iv2);
        printf("V3   = %.16f\n\n", iv3);
        #endif

        /*
           VETO: Choose the collocations points to use.
           Outside of the callibration region the collocation points can have a wavy behaviour so we remove some of them.
        */
        if(pWFHM->InspiralAmpVeto == 1){
          IMRPhenomXHM_Inspiral_Amplitude_Veto(&iv1, &iv2, &iv3, &powers_of_f1, &powers_of_f2, &powers_of_f3, pAmp, pWFHM);
        }
        #if DEBUG == 1
        printf("\nInspiral Veto: AmpVersion = %i",pWFHM->IMRPhenomXHMInspiralAmpVersion);
        #endif
        /* The 32 mode in this corner of the parameter space always uses the PN only */
        if(pWFHM->modeTag==32 && pWF22->q>2.5 && pWF22->chi1L<-0.9 && pWF22->chi2L<-0.9) pWFHM->IMRPhenomXHMInspiralAmpVersion = 0;
        /* The 32 mode in this corner of the parameter space removes the last collocation point */
        if(pWFHM->modeTag==32 && pWF22->q>2.5 && pWF22->chi1L<-0.6 && pWF22->chi2L>0. && iv1!=0){
          pWFHM->IMRPhenomXHMInspiralAmpVersion = pWFHM->IMRPhenomXHMInspiralAmpVersion - 1;
          iv1 = 0.;
        }
        /* The 33 mode in this corner of the parameter space removes the last collocation point */
        if(pWFHM->modeTag==33 && (1.2>pWF22->q && pWF22->q>1. && pWF22->chi1L<-0.1 && pWF22->chi2L>0. && iv1!=0)){//November
          pWFHM->IMRPhenomXHMInspiralAmpVersion = pWFHM->IMRPhenomXHMInspiralAmpVersion - 1;
          iv1 = 0.;
        }
        #if DEBUG == 1
        printf("\nInspiral Veto: AmpVersion = %i", pWFHM->IMRPhenomXHMInspiralAmpVersion);
        #endif
        //********* Remember that f3 < f2 < f1 !!!!! *****************
        /* Check for wavy collocation points. Only when we have 3 coll points. */
        if(pWFHM->IMRPhenomXHMInspiralAmpVersion == 3 && pAmp->WavyInsp == 1){
            if(WavyPoints(pAmp->CollocationPointsValuesAmplitudeInsp[0]*powers_of_f1.m_seven_sixths,pAmp->CollocationPointsValuesAmplitudeInsp[1]*powers_of_f2.m_seven_sixths,pAmp->CollocationPointsValuesAmplitudeInsp[2]*powers_of_f3.m_seven_sixths)==1){
               iv2 = 0;
               #if DEBUG == 1
               printf("\nWavy Inspiral\n");
               #endif
               pWFHM->IMRPhenomXHMInspiralAmpVersion = pWFHM->IMRPhenomXHMInspiralAmpVersion - 1;
             }
        }
        #if DEBUG == 1
        printf("\nInspiral Veto: AmpVersion = %i",pWFHM->IMRPhenomXHMInspiralAmpVersion);
        printf("\niv1 iv2 iv3 %e %e %e\n", iv1, iv2, iv3);
        #endif
        // Rename collocation points and frequencies.
        if(iv2==0){
          #if DEBUG == 1
              printf("\niv2 = 0\n");
          #endif
            iv2 = iv3;
            iv3 = 0.;
            powers_of_f2 = powers_of_f3;
        }
        if(iv1==0){
          #if DEBUG == 1
              printf("\niv1 = 0\n");
          #endif
            iv1 = iv2;
            powers_of_f1 = powers_of_f2;
            powers_of_fcutInsp = powers_of_f1;
            iv2 = iv3;
            iv3 = 0.;
            powers_of_f2 = powers_of_f3;
        }
        if(pWFHM->IMRPhenomXHMInspiralAmpVersion == 0)  // when using PN we take fcutInsp = fMECO_lm
        {
          IMRPhenomX_Initialize_Powers(&powers_of_fcutInsp, (pWFHM->fMECOlm));
        }

        // Update the Inspiral cutting frequency to the last collocation point alive.
        pAmp->fAmpMatchIN = powers_of_fcutInsp.itself;
        // End of VETO


        #if DEBUG == 1
        printf("\nAmplitude pseudo collocation points after veto\n");
        printf("fAmpMatchIN = %.16f\n", pAmp->fAmpMatchIN);
        printf("F1   = %.16f\n", powers_of_f1.itself);
        printf("F2   = %.16f\n", powers_of_f2.itself);
        printf("F3   = %.16f\n", powers_of_f3.itself);
        printf("V1   = %.16f\n", iv1);
        printf("V2   = %.16f\n", iv2);
        printf("V3   = %.16f\n", iv3);
        #endif

        /* Get the pseudo-PN coefficients. */
        /* The whole inspiral ansatz is given by PNAnsatz(f) + pseudo-PN(f),  with pseudo-PN(f) = rho1 *(f/fcutInsp)^(7/3) + rho2(f/fcutInsp)^(8/3) + rho3(f/fcutInsp)^(9/3).
        The coefficients are computed demanding that the pseudo-PN part at the three collocation points frequencies returns the actual collocation points: iv1, iv2, iv3.   */
        pAmp->rho1 = IMRPhenomXHM_Inspiral_Amp_rho1(iv1, iv2, iv3, &powers_of_fcutInsp, &powers_of_f1, &powers_of_f2, &powers_of_f3, pWFHM);
        pAmp->rho2 = IMRPhenomXHM_Inspiral_Amp_rho2(iv1, iv2, iv3, &powers_of_fcutInsp, &powers_of_f1, &powers_of_f2, &powers_of_f3, pWFHM);
        pAmp->rho3 = IMRPhenomXHM_Inspiral_Amp_rho3(iv1, iv2, iv3, &powers_of_fcutInsp, &powers_of_f1, &powers_of_f2, &powers_of_f3, pWFHM);

        #if DEBUG == 1
        printf("\nAmplitude pseudo PN coeffcients after veto\n");
        printf("rho1 = %.16f\n",pAmp->rho1);
        printf("rho2 = %.16f\n",pAmp->rho2);
        printf("rho3 = %.16f\n",pAmp->rho3);
        #endif

        // To avoid passing an extra argument to IMRPhenomXHM_Inspiral_Amp_Ansatz we store this powers in pAmp.
        pAmp->fcutInsp_seven_thirds = powers_of_fcutInsp.seven_thirds;
        pAmp->fcutInsp_eight_thirds = powers_of_fcutInsp.eight_thirds;
        pAmp->fcutInsp_three = powers_of_fcutInsp.three;


        /*****************/
        /*    RINGDOWN   */
        /*****************/
        #if DEBUG == 1
        printf("\n\n**** RINGDOWN ****\n\n");
        #endif

        // We have three "fitted" coefficients across parameter space: alambda, lambda and sigma (RDCoefficients[0,1,2]). Sigma will be constat for all the modes except the 21.
        pAmp->RDCoefficient[0] = fabs(pAmp->RingdownAmpFits[modeint*3](pWF22,pWFHM->IMRPhenomXHMRingdownAmpFitsVersion));
        pAmp->RDCoefficient[1] = pAmp->RingdownAmpFits[modeint*3+1](pWF22,pWFHM->IMRPhenomXHMRingdownAmpFitsVersion);
        pAmp->RDCoefficient[2] = pAmp->RingdownAmpFits[modeint*3+2](pWF22,pWFHM->IMRPhenomXHMRingdownAmpFitsVersion);
        pAmp->RDCoefficient[3] = 1./12.;
        
        /* (IMRPhenomXPNRUseTunedCoprec) Amplitude deviations e.g. for PNR's copreessing model*/
        if(pWFHM->modeTag==33){
          // pAmp->lambda  = pAmp->lambda   +  ( pWFHM->PNR_DEV_PARAMETER * pWFHM->MU1 );
          pAmp->RDCoefficient[0] = pAmp->RDCoefficient[0]  +  ( pWFHM->PNR_DEV_PARAMETER * pWFHM->MU3 );
          pAmp->RDCoefficient[2] = pAmp->RDCoefficient[2]  +  ( pWFHM->PNR_DEV_PARAMETER * pWFHM->MU4 );
        }

        #if DEBUG == 1
        printf("\nuseInspAnsatzRigndown = %i\n",pAmp->useInspAnsatzRingdown);
        printf("alambda = %.16f\r\n",pAmp->RDCoefficient[0]);
        #endif

        // For some cases with extreme spins there is almost no merger and the transition to the ringdown is very sharp.
        // The coefficients of the ringdown do not work well here and we take an approximation of the inspiral.
        if(pAmp->useInspAnsatzRingdown==1){
          // The Ringdown amp at fAmpMatchIM is set to be 0.9 the amplitude in the last inspiral collocation point
          IMRPhenomX_UsefulPowers powers_of_fRD;
          IMRPhenomX_Initialize_Powers(&powers_of_fRD, pAmp->fAmpMatchIM);
          pAmp->RDCoefficient[0] = 0.9*fabs(IMRPhenomXHM_Inspiral_Amp_Ansatz(&powers_of_fcutInsp, pWFHM, pAmp)/(IMRPhenomXHM_RD_Amp_Ansatz(&powers_of_fRD, pWFHM, pAmp)/pAmp->RDCoefficient[0]));
        }

        #if DEBUG == 1
        printf("alambda = %.16f\r\n",pAmp->RDCoefficient[0]);
        printf("lambda  = %.16f\r\n",pAmp->RDCoefficient[1]);
        printf("sigma   = %.16f\r\n",pAmp->RDCoefficient[2]);
        printf("lc      = %.16f\r\n",pAmp->RDCoefficient[3]);
        #endif


        /*******************/
        /*   INTERMEDIATE  */
        /*******************/

        #if DEBUG == 1
        printf("\n\n**** INTERMEDIATE ****\n\n");
        #endif

        // Initialize values of collocation points. They are taken from the parameter space fits.
        for(int i = 0; i<nCollocPtsInterAmp; i++)
        {
          pAmp->CollocationPointsValuesAmplitudeInter[i] = fabs(pAmp->IntermediateAmpFits[modeint*nCollocPtsInterAmp+i](pWF22,pWFHM->IMRPhenomXHMIntermediateAmpFitsVersion));
        }

        /* For reconstructing the intermediate region we need value and derivative at the beginning and at the end of the intermediate region, given by the inspiral and ringdown ansatzs.
        And also the two values at the two intermediate collocation points. This gives 6 parameters, that will determine the 6 coefficients of the inverse of the 5th order polynomial
        that we use for reconstruction.
        This is the default behaviour, however the veto conditions can reduce the order of the polynomial.
        */

        double F1, F2, F3, F4;  //Frequencies of the collocation points in the intermediate region
        F1     = powers_of_fcutInsp.itself;
        F2     = pAmp->CollocationPointsFreqsAmplitudeInter[0];
        F3     = pAmp->CollocationPointsFreqsAmplitudeInter[1];
        F4     = pAmp->fAmpMatchIM;

        #if DEBUG == 1
        printf("F1 = %.16f\n",F1);
        printf("F2 = %.16f\n",F2);
        printf("F3 = %.16f\n",F3);
        printf("F4 = %.16f\n\n",F4);
        #endif

        // Initialize useful powers.
        IMRPhenomX_UsefulPowers powers_of_F1;
        IMRPhenomX_Initialize_Powers(&powers_of_F1,F1);
        IMRPhenomX_UsefulPowers powers_of_F4;
        IMRPhenomX_Initialize_Powers(&powers_of_F4,F4);

        // Compute values at the boundaries (rescaled ansatz with the leading order of the 22).
        double inspF1 = IMRPhenomXHM_Inspiral_Amp_Ansatz(&powers_of_F1, pWFHM, pAmp);
        double rdF4;
        if (pWFHM->MixingOn == 1){
          rdF4 = cabs(SpheroidalToSpherical(&powers_of_F4, pAmp22, pPhase22, pAmp, pPhase, pWFHM, pWF22));
        }else{
          rdF4 = IMRPhenomXHM_RD_Amp_Ansatz(&powers_of_F4, pWFHM, pAmp);
        }

        // Compute derivatives at the boundaries (rescaled ansatz with the leading order of the 22).
        double d1, d4;
        d1 = IMRPhenomXHM_Inspiral_Amp_NDAnsatz(&powers_of_F1, pWFHM, pAmp);
        if(pWFHM->MixingOn==1){
          d4 = IMRPhenomXHM_RD_Amp_NDAnsatz(&powers_of_F4, pAmp, pPhase, pWFHM, pAmp22, pPhase22, pWF22);
        }else{
          d4 = IMRPhenomXHM_RD_Amp_DAnsatz(&powers_of_F4, pWFHM, pAmp);
        }
        // Next use of Ansatz will be for return the full strain, set correct rescalefactor
        pAmp->InspRescaleFactor = 0;
        pAmp->InterRescaleFactor = 0;
        pAmp->RDRescaleFactor = 0;

        #if DEBUG == 1
        printf("d1 = %.16f\n",d1);
        printf("d4 = %.16f\n",d4);
        #endif

        // Derivatives at the boundaries of the whole ansatz.
        d1 = ((7.0/6.0) * pow(F1,1.0/6.0) / inspF1) - ( pow(F1,7.0/6.0) * d1 / (inspF1*inspF1) );
        d4 = ((7.0/6.0) * pow(F4,1.0/6.0) / rdF4)   - ( pow(F4,7.0/6.0) * d4 / (rdF4*rdF4) );

        #if DEBUG == 1
        printf("d1 = %.16f\n",d1);
        printf("d4 = %.16f\n\n",d4);
        #endif

        // These values will feed the reconstruction function of the intermediate for getting the coefficients of the polynomial
        double V1, V2, V3, V4;

        pWFHM->IMRPhenomXHMIntermediateAmpVersion = 105; // by default use 5th order polynomial

        V1 = powers_of_F1.m_seven_sixths * inspF1;
        V2 = pAmp->CollocationPointsValuesAmplitudeInter[0];
        V3 = pAmp->CollocationPointsValuesAmplitudeInter[1];
        V4 = powers_of_F4.m_seven_sixths * rdF4;
        
        #if DEBUG == 1
        printf("Before intermediate veto \n");
        printf("V1 = %.16f\n",V1);
        printf("V2 = %.16f\n",V2);
        printf("V3 = %.16f\n",V3);
        printf("V4 = %.16f\n",V4);
        printf("rdF4 = %.16f\n",rdF4);
        printf("F4.m_seven_sixths = %.16f\n\n",powers_of_F4.m_seven_sixths);
        #endif

        // VETO: outside of the calibration region some collocation points are bad behaved. We veto them and change the type of reconstruction now.

        if(pAmp->useInspAnsatzRingdown==1){
          //For these extreme cases we do not use intermediate collocation points -> third order polynomial.
          V2 = 1.;
          V3 = 1.;
          pWFHM->IMRPhenomXHMIntermediateAmpVersion = 101; /*Linear reconstruction*/
          #if DEBUG == 1
          printf("VETO: useInspAnsatzRingdown\n");
          printf("V2 = %.16f\n",V2);
          printf("V3 = %.16f\n",V3);
          #endif
        }

        // The reconstruction function is the inverse of a polynomial. That we demand pass through the points V1, V2, V3, V4 and has the derivatives d1, d4 at the boundaries.
        // For simplicity we coded the reconstruction of the polynomial (Update_Intermediate_Amplitude_Coefficients), so we have to feed it with inverse values.
        V1  = 1.0 / V1;
        V2  = 1.0 / V2;
        V3  = 1.0 / V3;
        V4  = 1.0 / V4;
        
        /* Apply NR tuning for precessing cases (IMRPhenomXPNRUseTunedCoprec) */
        if(pWFHM->modeTag==33){
          V2 = V2 + ( pWFHM->PNR_DEV_PARAMETER * pWFHM->MU1 );
          V3 = V3 + ( pWFHM->PNR_DEV_PARAMETER * pWFHM->MU2 );
        }


        #if DEBUG == 1
        printf("\nAfter Veto and inverse \n");
        printf("V1 = %.16f\n",V1);
        printf("V2 = %.16f\n",V2);
        printf("V3 = %.16f\n",V3);
        printf("V4 = %.16f\n",V4);
        printf("IMRPhenomXHMIntermediateAmpVersion = %i \n\n", pWFHM->IMRPhenomXHMIntermediateAmpVersion);
        #endif


        /* If the case does not have extreme mass ratio it skips the next block */


        /****** JUST EMR cases, 2 INTERMEDIATE REGIONS ******/

        /* For EMR cases the amplitude shows a more pronounced drop off at the end of the inspiral.
        To model this we split the intermediate region in two.
        We add an extra collocation point between the end of the inspiral and the first intermediate collocation point, F0.
        From fcutInsp to FO we use a 4th order polynomial, and from F0 to fcutRD we use the usual 5th that is computed after this block. */

        if(pWFHM->AmpEMR==1){
          #if DEBUG == 1
          printf("*** TWO INTERMEDIATE REGIONS ***\n\n");
          #endif

          // Here we compute the first intermediate region with the inverse of a 4th polynomial.
          // For this we use point and derivative at fcutInsp and F0 (4 coefficients) and the value at the first collocation point (F2), with the total of 5 coefficients.

          double F0, V0, d0, F0_seven_sixths;  //Frequency, value, derivative and useful power at the extra collocation point for EMR

          F0 = F1 + (F2-F1)/3.;
          pAmp->fAmpMatchInt12 = F0;

          // Take the value and derivative from the parameter space fits.
          V0 = pAmp->IntermediateAmpFits[modeint*nCollocPtsInterAmp+8](pWF22,pWFHM->IMRPhenomXHMIntermediateAmpFitsVersion);
          d0 = pAmp->IntermediateAmpFits[modeint*nCollocPtsInterAmp+9](pWF22,pWFHM->IMRPhenomXHMIntermediateAmpFitsVersion);

          #if DEBUG == 1
          printf("F0 = %.16f\n",F0);
          printf("V0 = %.16f\n",V0);
          printf("d0 = %.16f\n",d0);
          #endif

          F0_seven_sixths = pow(F0,7.0/6.0);

          d0 = ((7.0/6.0) / (V0*F0))   - ( d0 / (V0*V0*F0_seven_sixths) );
          V0 = 1. / V0;

          #if DEBUG == 1
          printf("1/V0 = %.16f\n",V0);
          printf("1/d0 = %.16f\n",d0);
          #endif

          // Get the coefficients of the polynomial for the first intermediate region
          pAmp->alpha0 = IMRPhenomXHM_Intermediate_Amp_delta0(d1,d0,V1,V2,V3,V0,F1,F2,F3,F0,104); //V3 and F3 will not be used when calling with 104
          pAmp->alpha1 = IMRPhenomXHM_Intermediate_Amp_delta1(d1,d0,V1,V2,V3,V0,F1,F2,F3,F0,104);
          pAmp->alpha2 = IMRPhenomXHM_Intermediate_Amp_delta2(d1,d0,V1,V2,V3,V0,F1,F2,F3,F0,104);
          pAmp->alpha3 = IMRPhenomXHM_Intermediate_Amp_delta3(d1,d0,V1,V2,V3,V0,F1,F2,F3,F0,104);
          pAmp->alpha4 = IMRPhenomXHM_Intermediate_Amp_delta4(d1,d0,V1,V2,V3,V0,F1,F2,F3,F0,104);

          #if DEBUG == 1
          printf("Intermediate 1: feed values \n");
          printf("d1 = %.16f\n", d1);
          printf("d0 = %.16f\n", d0);
          printf("d4 = %.16f\n", d4);
          printf("V1 = %.16f\n", V1);
          printf("V0 = %.16f\n", V0);
          printf("V2 = %.16f\n", V2);
          printf("V3 = %.16f\n", V3);
          printf("V4 = %.16f\n", V4);
          printf("F1 = %.16f\n", F1);
          printf("F0 = %.16f\n", F0);
          printf("F2 = %.16f\n", F2);
          printf("F3 = %.16f\n", F3);
          printf("F4 = %.16f\n", F4);
          #endif

          #if DEBUG == 1
          printf("\nIntermediate 1: polynomial coeffcients \r\n");
          printf("alpha0 = %.16f\n", pAmp->alpha0);
          printf("alpha1 = %.16f\n", pAmp->alpha1);
          printf("alpha2 = %.16f\n", pAmp->alpha2);
          printf("alpha3 = %.16f\n", pAmp->alpha3);
          printf("alpha4 = %.16f\n", pAmp->alpha4);
          #endif

          //Update left collocation point for the 2nd intermediate region
          F1 = F0;
          V1 = V0;
          d1 = d0;

          /**** END of first Intermediate region ****/
        }
        else{  /** This part is used both when we have a single intermediate region and for the second intermediate region **/
          pAmp->fAmpMatchInt12 = 0;
          pAmp->alpha0 = 1;
          pAmp->alpha1 = 1;
          pAmp->alpha2 = 1;
          pAmp->alpha3 = 1;
          pAmp->alpha4 = 1;

          /** More vetos ***/
          if(pWFHM->IntermediateAmpVeto == 1 && pWFHM->IMRPhenomXHMIntermediateAmpVersion==105){  // only 21 mode
            IMRPhenomXHM_Intermediate_Amplitude_Veto(&V2, &V3, pWFHM, pWF22); // this changes the order of the polynomial to 4 or 3
            #if DEBUG == 1
            printf("VETO: Intermediate Amp Veto\n");
            printf("V2 = %.16f\n",V2);
            printf("V3 = %.16f\n",V3);
            #endif
          }

          if(pWFHM->RingdownAmpVeto == 1){    // only 21, 32 mode
            IMRPhenomXHM_Ringdown_Amplitude_Veto(&V2, &V3, V4, pWFHM, pWF22); // If satisfied, remove the 2 inter collocation points
            #if DEBUG == 1
            printf("VETO: Ringdown Amp Veto\n");
            printf("V2 = %.16f\n",V2);
            printf("V3 = %.16f\n",V3);
            #endif
          }

          if(pWFHM->IMRPhenomXHMIntermediateAmpVersion==105 && pAmp->WavyInt==1){
            if(WavyPoints(V2,V3,V4)==1){
              V3 = V2;
              F3 = F2;
              V2 = 1.;
              pWFHM->IMRPhenomXHMIntermediateAmpVersion=1042;
              #if DEBUG == 1
              printf("VETO: Wavy Inter colloc points\n");
              printf("V2 = %.16f\n",V2);
              printf("V3 = %.16f\n",V3);
              #endif
            }
          }

          if((pWF22->q>40. && pWF22->chi1L>0.9 && V2!=1 && V3!=1) ||
          (pWFHM->modeTag==32 && (pWFHM->IMRPhenomXHMIntermediateAmpVersion != 101) && ((pWF22->q>2.5 && pWF22->chi1L<-0.6 && pWF22->chi2L>0) || (pWF22->chi1L<-0.9&&pWF22->chi2L<-0.9))) ||
          (pWFHM->modeTag==21 && pWF22->eta<0.23 && pWF22->chi1L>0.7 && pWF22->chi2L<-0.5)){
            V2 = 1.;
            V3 = 1.;
            pWFHM->IMRPhenomXHMIntermediateAmpVersion = 1032;
            #if DEBUG == 1
            printf("VETO: veto regions\n");
            printf("V2 = %.16f\n",V2);
            printf("V3 = %.16f\n",V3);
            #endif
          }
          /*** End of vetos **/
        }

        // The reconstruction function (Update_Intermediate_Amplitude_Coefficients) assumes that F3 is the point with value.
        // If F3 was removed in the veto the we replace it with F2
        if(V3 == 1.){
          V3 = V2;
          F3 = F2;
          V2 = 1.;
        }

        /*
        Reconstruct the phenomenological coefficients for the intermediate ansatz
        */

        // Store the values for the reconstruction.
        pAmp->v1 = V1;
        pAmp->v2 = V2;
        pAmp->v3 = V3;
        pAmp->v4 = V4;
        pAmp->f1 = F1;
        pAmp->f2 = F2;
        pAmp->f3 = F3;
        pAmp->f4 = F4;
        pAmp->d1 = d1;
        pAmp->d4 = d4;

        #if DEBUG == 1
        printf("\nIntermediate Amplitude Input \r\n");
        printf("\nIMRPhenomXHMIntermediateAmpVersion = %i \r\n", pWFHM->IMRPhenomXHMIntermediateAmpVersion);
        printf("V1 = %.16f\n", V1);
        printf("V2 = %.16f\n", V2);
        printf("V3 = %.16f\n", V3);
        printf("V4 = %.16f\n", V4);
        printf("F1 = %.16f\n", F1);
        printf("F2 = %.16f\n", F2);
        printf("F3 = %.16f\n", F3);
        printf("F4 = %.16f\n", F4);
        printf("d1 = %.16f\n", d1);
        printf("d4 = %.16f\n", d4);
        printf("fAmpMatchIn = %.16f\n", pAmp->fAmpMatchIN);
        #endif

        // Compute the coefficients of the polynomial
        Update_Intermediate_Amplitude_Coefficients(pAmp, pWFHM->IMRPhenomXHMIntermediateAmpVersion);


        #if DEBUG == 1
        printf("\nIntermediate polynomial coeffcients Before ChoosePolOrder \r\n");
        printf("\nIMRPhenomXHMIntermediateAmpVersion = %i \r\n", pWFHM->IMRPhenomXHMIntermediateAmpVersion);
        printf("delta0 = %.16f\n", pAmp->delta0);
        printf("delta1 = %.16f\n", pAmp->delta1);
        printf("delta2 = %.16f\n", pAmp->delta2);
        printf("delta3 = %.16f\n", pAmp->delta3);
        printf("delta4 = %.16f\n", pAmp->delta4);
        printf("delta5 = %.16f\n", pAmp->delta5);
        printf("fAmpMatchIn = %.16f\n", pAmp->fAmpMatchIN);
        #endif

        // Check that the polynomial does not cross zero, because the actual reconstructing function is the inverse of this polynomial
        // If it crosses zero, then remove one collocation and lower the order of the polynomial.
        ChoosePolOrder(pWFHM, pAmp);

        #if DEBUG == 1
        printf("\nIMRPhenomXHMIntermediateAmpVersion = %i \r\n", pWFHM->IMRPhenomXHMIntermediateAmpVersion);
        printf("\nIntermediate polynomial coeffcients After ChoosePolOrder\r\n");
        printf("delta0 = %.16f\n", pAmp->delta0);
        printf("delta1 = %.16f\n", pAmp->delta1);
        printf("delta2 = %.16f\n", pAmp->delta2);
        printf("delta3 = %.16f\n", pAmp->delta3);
        printf("delta4 = %.16f\n", pAmp->delta4);
        printf("delta5 = %.16f\n", pAmp->delta5);
        printf("fAmpMatchIn = %.16f\n", pAmp->fAmpMatchIN);
        #endif
    }// END of 122018 version

  }

double RescaleFactor(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMAmpCoefficients *pAmp, UINT2 rescalefactor){
  double factor = 0.;
  switch(rescalefactor){
      case 0:{ // Strain
         factor = 1.;
         break;
      }          
      case 1:{ // 22factor
          factor = pAmp->ampNorm * powers_of_Mf->m_seven_sixths;
          break;
      }
      case 2:{ // lmfactor
          if (pAmp->PNdominantlmpower == 1) factor = pAmp->PNdominant * powers_of_Mf->m_seven_sixths * pAmp->PNdominantlm * powers_of_Mf->one_third;
          if (pAmp->PNdominantlmpower == 2) factor = pAmp->PNdominant * powers_of_Mf->m_seven_sixths * pAmp->PNdominantlm * powers_of_Mf->two_thirds;
          if (pAmp->PNdominantlmpower == 3) factor = pAmp->PNdominant * powers_of_Mf->m_seven_sixths * pAmp->PNdominantlm * powers_of_Mf->itself;
          if (pAmp->PNdominantlmpower == 4) factor = pAmp->PNdominant * powers_of_Mf->m_seven_sixths * pAmp->PNdominantlm * powers_of_Mf->four_thirds;
          break;
      }
      default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in RescaleFactor: version %i is not valid. Recommended version is 1.", rescalefactor);}
  }
  return factor;
}


/*******************************************/
/*                                         */
/*          COMPUTE PHASE COEFFICIENTS     */
/*                                         */
/*******************************************/

/*  This function is an expansion in the frequency of the real part of the 21-amplitude, following the convention of Eq.(2.2) in the paper  */
/*  This function occasionally changes sign, so we need to keep track of this to correct the WF's phase, since the amplitude of IMRPhenomXHM is positive by construction */
int IMRPhenomXHM_PN21AmpSign (double ff,IMRPhenomXWaveformStruct *wf22){

    double eta,chi1,chi2;
    eta =wf22->eta;
    chi1=wf22->chi1L;
    chi2=wf22->chi2L;
    double delta=sqrt(1 - 4*eta);

    double output=(-16*delta*eta*ff*pow(LAL_PI,1.5))/(3.*sqrt(5)) + (4*pow(2,0.3333333333333333)*(chi1 - chi2 + delta *(chi1+ chi2))*eta*pow(ff,1.3333333333333333)*pow(LAL_PI,1.8333333333333333))/sqrt(5) + (2*pow(2,0.6666666666666666)*eta*(306*delta - 360*delta*eta)*pow(ff,1.6666666666666667)*pow(LAL_PI,2.1666666666666665))/(189.*sqrt(5));

    if(output>=0) return(1);
    else return(-1);
  }

/*  compute phase coefficients */
void IMRPhenomXHM_GetPhaseCoefficients(IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXAmpCoefficients *pAmp22, IMRPhenomXPhaseCoefficients *pPhase22, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22, UNUSED LALDict *lalParams) {

    int ell=pWFHM->ell;
    int emm=pWFHM->emm;

    /* Pre-initialize all phenomenological coefficients */

    //inspiral
    // these will be computed by rescaling IMRPhenomX inspiral coefficients
    for(int i=0; i<NMAX_INSPIRAL_COEFFICIENTS; i++)
    {
      pPhase->phi[i]=0.;
      pPhase->phiL[i]=0.;
    }

    //intermediate
    // will be determined by solving linear system at collocation points (see Eqs. (5.5)-(5.6) in the paper)
    pPhase->c0 = 0.0;
    pPhase->c1 = 0.0;
    pPhase->c2 = 0.0;
    pPhase->c3 = 0.0;
    pPhase->c4 = 0.0;
    pPhase->cL = 0.0;
    // ringdown spherical
    // these will be obtained by rescaling the 22 ringdown parameters (see Eq. (6.5) in the paper)
    pPhase->alpha0 = 0.0;
    pPhase->alpha2 = 0.0;
    pPhase->alphaL = 0.0;

  
    
    // set number of collocation points used (depends on the mode)
    int nCollocationPts_inter=pWFHM->nCollocPtsInterPhase;
    // define mass-ratio that discriminates between EMR and rest of cases
    double eta_m1=1./pWF22->eta;
    // store mode tag
    int modeint=pWFHM->modeInt;

    IMRPhenomX_UsefulPowers powers_of_f;

    //initialize frequencies of colloc points in intermediate region
    IMRPhenomXHM_Intermediate_CollocPtsFreqs(pPhase,pWFHM,pWF22);


    //for each collocation point, call fit giving the value of the phase derivative at that point
    for(int i = 0; i<N_MAX_COEFFICIENTS_PHASE_INTER; i++)
    {
        pPhase->CollocationPointsValuesPhaseInter[i] = pPhase->IntermediatePhaseFits[modeint*N_MAX_COEFFICIENTS_PHASE_INTER+i](pWF22,pWFHM->IMRPhenomXHMIntermediatePhaseFitsVersion);
        // time-shift wf so that modes peak around t=0
        // our hybrids are built so that their inverse FT peaks ~500 M before the end of the time-interval. Our reconstructed WFs will have the same normalization, so we shift the phase here so that the modes peak around t=0
        pPhase->CollocationPointsValuesPhaseInter[i]=pPhase->CollocationPointsValuesPhaseInter[i]+pWFHM->DeltaT;
    }

    double fcutRD=pPhase->fPhaseMatchIM;
    double fcutInsp=pPhase->fPhaseMatchIN;



    /************** Inspiral: rescale PhenX and apply PN corrections *********/

    //collect all the phenX inspiral coefficients
    double phenXnonLog[]={pPhase22->phi0,pPhase22->phi1,pPhase22->phi2,pPhase22->phi3,pPhase22->phi4,pPhase22->phi5,pPhase22->phi6,pPhase22->phi7,pPhase22->phi8,pPhase22->phi9,0.,0,0,0};
    double phenXLog[]={0.,0.,0.,0.,0.,pPhase22->phi5L,pPhase22->phi6L,0.,pPhase22->phi8L,pPhase22->phi9L,0,0,0,0};
    double pseudoPN[]={0.,0.,0.,0.,0.,0.,0.,0.,pPhase22->sigma1,pPhase22->sigma2,pPhase22->sigma3,pPhase22->sigma4,pPhase22->sigma5};
    
  
    
    //rescale the coefficients of phenX by applying phi_lm(f)~m/2 *phi_22(2/m f)
    // for more details, see Appendix D of the paper
    double m_over_2=emm*0.5, two_over_m=1./m_over_2;

    double fact=pPhase22->phiNorm/pWF22->eta;

    for(int i=0; i< NMAX_INSPIRAL_COEFFICIENTS; i++){

      phenXnonLog[i]=phenXnonLog[i]*fact;
      phenXLog[i]=phenXLog[i]*fact;
      pseudoPN[i]=pseudoPN[i]*fact;

      // scaling the logarithmic terms introduces some extra contributions in the non-log terms
      pPhase->phi[i]=(phenXnonLog[i]+pseudoPN[i]-phenXLog[i]*log(m_over_2))*pow(m_over_2,(8-i)/3.);
      pPhase->phiL[i]=phenXLog[i]*pow(m_over_2,(8-i)/3.);
    }


    // if the mass-ratio is not extreme, use the complex PN amplitudes to correct the orbital phase at linear order in f
    // see Eq. (4.12) in the paper
    if(pWF22->eta>0.01){
      pPhase->LambdaPN=IMRPhenomXHM_Insp_Phase_LambdaPN(pWF22->eta, pWFHM->modeTag);
      #if DEBUG == 1
      printf("Add phase shift from PN complex amplitude\n");
      #endif

    }
    // else use some phenomenological fits: the fits give the coefficient of a linear-in-f term to be added to the orbital phase
    else{

      pPhase->LambdaPN= pPhase->InspiralPhaseFits[pWFHM->modeInt](pWF22,pWFHM->IMRPhenomXHMInspiralPhaseVersion);
      #if DEBUG == 1
      printf("\nLambdaPN = %.16f\n", pPhase->LambdaPN);
      #endif
    }

    pPhase->phi[8]= pPhase->phi[8]+ pPhase->LambdaPN;


    /****************************************
    *********** Intermediate-ringdown region ********
    ****************************************/

    // below we set up the linear system to be solved in the intermediate region


    /* GSL objects for solving system of equations via LU decomposition */
    gsl_vector *b, *x;
    gsl_matrix *A;
    gsl_permutation *p;
    int s;

    p = gsl_permutation_alloc(nCollocationPts_inter);
    b = gsl_vector_alloc(nCollocationPts_inter);
    x = gsl_vector_alloc(nCollocationPts_inter);
    A = gsl_matrix_alloc(nCollocationPts_inter,nCollocationPts_inter);


    // for high-spin cases: avoid sharp transitions for 21 mode by extending the inspiral region into the intermediate one

    if(pWFHM->modeTag==21&&pWF22->STotR>=0.8){

        double insp_vals[3], FF, diff12, diff23;
        for(int i=0; i<3; i++){
            FF=two_over_m*pPhase->CollocationPointsFreqsPhaseInter[i];
            IMRPhenomX_Initialize_Powers(&powers_of_f,FF);
            insp_vals[i]=1./pWF22->eta*IMRPhenomX_dPhase_22(FF,&powers_of_f,pPhase22,pWF22);
        }

        diff12=insp_vals[0]-insp_vals[1];
        diff23=insp_vals[1]-insp_vals[2];

        pPhase->CollocationPointsValuesPhaseInter[1]=pPhase->CollocationPointsValuesPhaseInter[2]+diff23;
        pPhase->CollocationPointsValuesPhaseInter[0]=pPhase->CollocationPointsValuesPhaseInter[1]+diff12;

    }


    // choose collocation points according to spin/mass ratio
    // current catalogue of simulations include some cases that create unphysical effects in the fits -> we need to use different subset of collocation points according to the parameters (we have to pick 5 out of 6 available fits)
    /* cpoints_indices is an array of integers labelling the collocation points chosen in each case, e.g.
     cpoints_indices={0,1,3,4,5} would mean that we are discarding the 3rd collocation points in the reconstructio */

    int cpoints_indices[nCollocationPts_inter];
    cpoints_indices[0]=0;
    cpoints_indices[1]=1;
    cpoints_indices[4]=5;


    if((pWF22->eta<pWFHM->etaEMR)||(emm==ell&&pWF22->STotR>=0.8)||(pWFHM->modeTag==33&&pWF22->STotR<0))
    {
        cpoints_indices[2]=3;
        cpoints_indices[3]=4;
    }
    else if(pWF22->STotR>=0.8&&pWFHM->modeTag==21){

        cpoints_indices[2]=2;
        cpoints_indices[3]=4;
    }

    else{
        cpoints_indices[2]=2;
        cpoints_indices[3]=3;
    }




    switch(pWFHM->MixingOn){
            // mode-mixing off: mode does not show significant mode-mixing
            // the MixingOn flag is automatically switched off for the modes 21,33,44
      case 0:
      {

          /************ Ringdown **************/

          /* ansatz:
           alpha0 + ((fRDlm^2) alpha2)/(f^2)  + alphaL*(fdamplm)/((fdamplm)^2 + (f - fRDlm)^2)*/

        //compute alpha2 by rescaling the coefficient of the 22-reconstruction
        double wlm;
        if(pWFHM->ell==pWFHM->emm) wlm=2;
        else wlm=pWFHM->emm/3.;

        pPhase->alpha2=1./pow(pWFHM->fRING,2)*wlm*IMRPhenomXHM_RD_Phase_22_alpha2(pWF22,pWFHM->IMRPhenomXHMRingdownPhaseVersion);


        // compute alphaL
        pPhase->alphaL=eta_m1*IMRPhenomXHM_RD_Phase_22_alphaL(pWF22,pWFHM->IMRPhenomXHMRingdownPhaseVersion);

        #if DEBUG == 1
        printf("alpha2_fit=%.16e\n",IMRPhenomXHM_RD_Phase_22_alpha2(pWF22,pWFHM->IMRPhenomXHMRingdownPhaseVersion));
        printf("Ringdown parameters: \n alpha2=%.16e \n alphaL=%.16e\n",pPhase->alpha2,pPhase->alphaL);
        # endif

        // compute spherical-harmonic phase and its first derivative at matching frequency --used to glue RD and intermediate regions--
        IMRPhenomX_Initialize_Powers(&powers_of_f,fcutRD);
        pPhase->phi0RD=IMRPhenomXHM_RD_Phase_AnsatzInt(fcutRD, &powers_of_f,pWFHM, pPhase);
        pPhase->dphi0RD=IMRPhenomXHM_RD_Phase_Ansatz(fcutRD, &powers_of_f,pWFHM, pPhase);


          /************ Intermediate **************/


        // set up the linear system by calling the fits of the intermediate-region collocation points
        for(int i=0; i<nCollocationPts_inter; i++){

          int ind=cpoints_indices[i];
          gsl_vector_set(b,i,pPhase->CollocationPointsValuesPhaseInter[ind]);
          REAL8 ff=pPhase->CollocationPointsFreqsPhaseInter[ind], ffm1=1./ff, ffm2=ffm1*ffm1;
          REAL8 fpowers[]={1., (pWFHM->fDAMP)/(pow(pWFHM->fDAMP,2)+pow(ff-(pWFHM->fRING),2)),ffm1,ffm2,ffm2*ffm2, ffm1*ffm2};
          for(int j=0; j<nCollocationPts_inter; j++)
          gsl_matrix_set(A,i,j,fpowers[j]);
        }

        break;
      }

      // mode-mixing on (the MixingOn flag is automatically switched on when reconstructing the 32 mode)
      case 1:
      {

        // for the 32 mode, the ringdown waveform is computed outside of the dedicated amplitude and phase routines, by calling SpheroidalToSpherical
        // compute dphi/df and d2phi/df2 at fcutRD starting from rotation of spheroidal ansatz: these values will then enter the linear system that determines the coefficients in the intermediate region

        // we compute derivatives by applying finite-difference schemes. In total we will need the value of the phase at three points around fcutRD

        // we first compute the full spherical-harmonic-basis waveforms at three points
        double complex SphericalWF[3];
        double fstep=0.0000001;
        
        for(int i=0; i<3; i++){

          double FF=fcutRD+(i-1)*fstep;
          IMRPhenomX_UsefulPowers powers_of_FF;
          IMRPhenomX_Initialize_Powers(&powers_of_FF,FF);
          SphericalWF[i]=SpheroidalToSpherical(&powers_of_FF, pAmp22, pPhase22, pAmp, pPhase, pWFHM, pWF22);
        }

        long double phase_args[]={fmodl(carg(SphericalWF[0]),2.*LAL_PI),fmodl(carg(SphericalWF[1]),2.*LAL_PI),fmodl(carg(SphericalWF[2]),2.*LAL_PI)};

        // make sure that all the three points belong to the same branch of mod
        for(int i=0; i<3; i++){

            if(phase_args[i]>0)
                phase_args[i]-=2.*LAL_PI;

        }

        // store spherical-harmonic phase at fcutRD
        pPhase->phi0RD=phase_args[1];
        long double fstep_m1= 1./fstep;
        // we apply the FD schemes to get first and second phase derivatives at fcutRD
        pPhase->dphi0RD=0.5*fstep_m1*(phase_args[2]-phase_args[0]);
        double d2phi0RD=fstep_m1*fstep_m1*(phase_args[2]-2.*phase_args[1]+phase_args[0]);

        #if DEBUG == 1
        printf("dphi0ref=%.16e\t d2phi0ref=%.16e\n",pPhase->dphi0RD,d2phi0RD);
        # endif

          // To achieve a smoother transition with the intermediate region (IR), we feed into the intermediate-region reconstruction the first and second derivative of the reconstructed ringdown phase, see Eq. (5.6)


        // first derivative
          // if the mass-ratio is not extreme, we use the derivative computed using the FD scheme above, else we keep the original value of the fit
          // in the second case, we will therefore need to glue intermediate and ringdown region later
        if(pWF22->eta>pWFHM->etaEMR){
          pPhase->CollocationPointsFreqsPhaseInter[nCollocationPts_inter-2]=fcutRD;
          pPhase->CollocationPointsValuesPhaseInter[nCollocationPts_inter-2]=pPhase->dphi0RD;
        }
        // second derivative
        pPhase->CollocationPointsFreqsPhaseInter[nCollocationPts_inter-1]=fcutRD;
        pPhase->CollocationPointsValuesPhaseInter[nCollocationPts_inter-1]=d2phi0RD;


        // set up the linear system to determine the intermediate region coefficients
        for(int i=0; i<nCollocationPts_inter; i++){

          gsl_vector_set(b,i,pPhase->CollocationPointsValuesPhaseInter[i]);
          REAL8 ff=pPhase->CollocationPointsFreqsPhaseInter[i], ffm1=1./ff, ffm2=ffm1*ffm1;
          REAL8 fpowers[]={1., (pWFHM->fDAMP)/(pow(pWFHM->fDAMP,2)+pow(ff-(pWFHM->fRING),2)),ffm1,ffm2,ffm2*ffm2, ffm1*ffm2};
          for(int j=0; j<nCollocationPts_inter; j++)
          gsl_matrix_set(A,i,j,fpowers[j]);

        }

        // set the last point to equal the second derivative of the rotated spheroidal ansatz
        int cpoint_ind=nCollocationPts_inter-1;
        gsl_vector_set(b,cpoint_ind,pPhase->CollocationPointsValuesPhaseInter[cpoint_ind]);
        REAL8 ff=pPhase->CollocationPointsFreqsPhaseInter[cpoint_ind];
        REAL8 ffm1=1./ff, ffm2=ffm1*ffm1, ffm3=ffm2*ffm1, ffm4=ffm2*ffm2, ffm5=ffm3*ffm2;

        REAL8 fpowers[]={0.,-2*(pWFHM->fDAMP)*(ff-(pWFHM->fRING))/pow((pow(pWFHM->fDAMP,2)+pow(ff-(pWFHM->fRING),2)),2),-ffm2, -2.* ffm3, -4.*ffm5,-3* ffm4};
        for(int j=0; j<nCollocationPts_inter; j++){
          gsl_matrix_set(A,cpoint_ind,j,fpowers[j]);
        }


        #if DEBUG == 1
        printf("Collocation points for intermediate ansatz %i:\n",nCollocationPts_inter);
        for(int i=0;i<nCollocationPts_inter; i++)
        printf("p%d={%.13f,%.13f};\n",i,pPhase->CollocationPointsFreqsPhaseInter[i],pPhase->CollocationPointsValuesPhaseInter[i]);
        #endif

        break;



      }

      default:{XLALPrintError("Error in IMRPhenomXHM_GetPhaseCoefficients:mode-mixing not properly initialized.\n");}

    }


    // at this point we have initizialized all the collocation points in the intermediate region, so we can solve to get the coefficients of the ansatz

    /* In the intermediate region , the most general ansatz is
    (c0 + c1 /f + c2 /(f)^2 + c4 /(f)^4 + c3 /f^3 +cL fdamp/((fdamp)^2 + (f - fRD )^2))*/

    /* We now solve the system A x = b via an LU decomposition */
    gsl_linalg_LU_decomp(A,p,&s);
    gsl_linalg_LU_solve(A,p,b,x);

    pPhase->c0 = gsl_vector_get(x,0); // x[0]; // c0
    pPhase->cL = gsl_vector_get(x,1); // x[1] // cL
    pPhase->c1 = gsl_vector_get(x,2); // x[2]; // c1
    pPhase->c2 = gsl_vector_get(x,3); // x[3]; // c2
    pPhase->c4 = gsl_vector_get(x,4); // x[4]; // c4
    
    /* (IMRPhenomXPNRUseTunedCoprec) Add PNR or input deviations (e.g. tuning for coprecessing frame model) to cL, c1 and c4 in the same manner as is done for the l=m=2 moment */
    if ( pWFHM->modeTag == 33){
	    pPhase->c0 = pPhase->c0  +  ( pWFHM->PNR_DEV_PARAMETER * pWFHM->NU0 );
	    pPhase->cL = pPhase->cL  +  ( pWFHM->PNR_DEV_PARAMETER * pWFHM->NU4 );
      pPhase->c1 = pPhase->c1  +  ( pWFHM->PNR_DEV_PARAMETER * pWFHM->ZETA2 );
      pPhase->c4 = pPhase->c4  +  ( pWFHM->PNR_DEV_PARAMETER * pWFHM->ZETA1 );
    }
    

    // currently the 32 mode is calibrated using one extra point
    if ((pWFHM->modeTag)== 32)
    { pPhase->c3 = gsl_vector_get(x,5);  //c3

      // if the mass-ratio is extreme, we need to glue intermediate and ringdown region, as explained above
      if(pWF22->eta<pWFHM->etaEMR)
      {
        IMRPhenomX_Initialize_Powers(&powers_of_f,fcutRD);
        pPhase->c0=pPhase->c0+pPhase->dphi0RD-IMRPhenomXHM_Inter_Phase_Ansatz(fcutRD, &powers_of_f,pWFHM, pPhase);

      }
    }

    #if DEBUG == 1
    printf("Intermediate region coefficients:\n");
    printf("c0=%.16e\n c1=%.16e\n c2=%.16e\n c3=%.16e\n c4=%.16e\n cL=%.16e\n", pPhase->c0,pPhase->c1,pPhase->c2,pPhase->c3,pPhase->c4,pPhase->cL);
    #endif

    gsl_vector_free(b);
    gsl_vector_free(x);
    gsl_matrix_free(A);
    gsl_permutation_free(p);


    /****************** end of ringdown and inspiral reconstruction **************/

    //glue inspiral and intermediate
    IMRPhenomX_Initialize_Powers(&powers_of_f,fcutInsp);
    pPhase->C1INSP=IMRPhenomXHM_Inter_Phase_Ansatz(fcutInsp, &powers_of_f,pWFHM, pPhase)-IMRPhenomXHM_Inspiral_Phase_Ansatz(fcutInsp, &powers_of_f, pPhase);
    pPhase->CINSP=-(pPhase->C1INSP*fcutInsp)+IMRPhenomXHM_Inter_Phase_AnsatzInt(fcutInsp, &powers_of_f,pWFHM, pPhase)-IMRPhenomXHM_Inspiral_Phase_AnsatzInt(fcutInsp, &powers_of_f, pPhase);
    //glue ringdown and intermediate regions
    IMRPhenomX_Initialize_Powers(&powers_of_f,fcutRD);

    pPhase->C1RD=IMRPhenomXHM_Inter_Phase_Ansatz(fcutRD, &powers_of_f,pWFHM, pPhase)-pPhase->dphi0RD;
    pPhase->CRD=-pPhase->C1RD*fcutRD+IMRPhenomXHM_Inter_Phase_AnsatzInt(fcutRD, &powers_of_f,pWFHM, pPhase)-pPhase->phi0RD;
    
    
    // we now have a C1 reconstruction of the phase
    // below we align each mode so that at low-f its relative phase wrt the 22 agrees with PN
    double falign;
    //somehow arbitary cutoff to pick frequency for alignment: must be in the inspiral region
    if(pWF22->eta>pWFHM->etaEMR)
    falign=0.6*m_over_2*pWF22->fMECO;
    else
    falign=m_over_2*pWF22->fMECO;
    // printf("0>> falign = %f\n",falign);
    
    // // (ltl) 0.6 is found to be too large a factor, so let's try a smaller one
    // if (pWF22->IMRPhenomXPNRUseTunedCoprec && (pWFHM->modeTag==33)) {
    //     if(pWF22->eta>pWFHM->etaEMR){
    //       falign = 0.06*m_over_2*pWF22->fMECO;
    //       printf("&>> falign = %f\n",falign);
    //     } else {
    //       falign = m_over_2*pWF22->fMECO;
    //     }
    // }

    IMRPhenomX_UsefulPowers powers_of_falign;
    IMRPhenomX_Initialize_Powers(&powers_of_falign,falign);
    IMRPhenomX_Initialize_Powers(&powers_of_f,two_over_m*falign);

    /* compute explicitly the phase normalization to be applied to IMRPhenomX, when mode-mixing is on this will have been already computed in GetSpheroidalCoefficients */
    if(pWFHM->MixingOn==0){
      IMRPhenomX_UsefulPowers powers_of_MfRef;
      IMRPhenomX_Initialize_Powers(&powers_of_MfRef,pWF22->MfRef);
      IMRPhenomX_Phase_22_ConnectionCoefficients(pWF22,pPhase22);
      pWFHM->timeshift=IMRPhenomX_TimeShift_22(pPhase22, pWF22);
      pWFHM->phiref22 = -1./pWF22->eta*IMRPhenomX_Phase_22(pWF22->MfRef, &powers_of_MfRef, pPhase22, pWF22) - pWFHM->timeshift*pWF22->MfRef - pWFHM->phaseshift + 2.0*pWF22->phi0 + LAL_PI_4;
    }

    // we determine the phase normalization of each mode by imposing Eq. (4.13), i.e. phi_lm(fref)-m/2 phi_22(2/m fref)~3.*LAL_PI_4*(1-m_over_2)
    double deltaphiLM = m_over_2*(1./pWF22->eta*IMRPhenomX_Phase_22(two_over_m*falign, &powers_of_f,pPhase22,pWF22)+pWFHM->phaseshift + pWFHM->phiref22)+pWFHM->timeshift*falign-3.*LAL_PI_4*(1-m_over_2)-(IMRPhenomXHM_Inspiral_Phase_AnsatzInt(falign, &powers_of_falign,pPhase)+pPhase->C1INSP*falign+pPhase->CINSP);
    pPhase->deltaphiLM=fmod(deltaphiLM,2.*LAL_PI);

    // /* (ll) We have commented out the code above in order to rewrite the same code in a more carefully formatted way */
    
    // // (ll) Define the relative phase parameter used in Eq. (4.13)
    // double relative_phase_at_zero_freq = 3.*LAL_PI_4*(1-m_over_2);
    
    // // (ll) Define the (2,2) inspiral phase evaluated at (2/m)*falign, and then rescaled by m/2
    // double rescaled_phi22_at_falign = m_over_2*(1./pWF22->eta*IMRPhenomX_Phase_22(two_over_m*falign, &powers_of_f,pPhase22,pWF22)+pWFHM->phaseshift + pWFHM->phiref22)+pWFHM->timeshift*falign;
    
    // //
    // printf("*>> falign = %f\n",falign);
    // printf("*>> relative_phase_at_zero_freq = %f\n",relative_phase_at_zero_freq);
    // printf("*>> m_over_2 = %f\n\n",m_over_2);
    
    // // (ll) Define XHM inspiral phase at align frequency 
    // double phiLM_inspiral_at_falign = (IMRPhenomXHM_Inspiral_Phase_AnsatzInt(falign, &powers_of_falign,pPhase)+pPhase->C1INSP*falign+pPhase->CINSP);

    // // (ll) 
    // double deltaphiLM = rescaled_phi22_at_falign - ( relative_phase_at_zero_freq + phiLM_inspiral_at_falign );
    
    // // (ll) Get the representation closest to 2*pi
    // pPhase->deltaphiLM=fmod(deltaphiLM,2.*LAL_PI);

    #if DEBUG == 1
    printf("\n****** Connection coefficients of 22 in %i******\n", pWFHM->modeTag);
    printf("%.6f %.6f %.6f %.6f %.6f\n", pPhase22->C1Int, pPhase22->C2Int, pPhase22->C1MRD, pPhase22->C1MRD, IMRPhenomX_Phase_22(two_over_m*falign, &powers_of_f,pPhase22,pWF22));
    #endif


    // for the 21, we need to make sure the sign of the PN amplitude is positive, else we'll need to flip its phase by Pi
    if ((pWFHM->modeTag)== 21)
    {
      int ampsign=IMRPhenomXHM_PN21AmpSign(0.008,pWF22);
      // the sign of the 21 amplitude changes across the parameter space, we add Pi to the phase to account for this - remember that the amplitude of the model is positive by construction
      if(ampsign>0) pPhase->deltaphiLM+=LAL_PI;

    }

    // all the coefficients have been now stored
  // print parameters to file if in debugging mode
  #if DEBUG == 1
  FILE *file;
  char fileSpec[40];
  sprintf(fileSpec, "PhaseParameters%i.dat", pWFHM->modeTag);

  file = fopen(fileSpec,"w");

  fprintf(file,"\n*** %i Mode ***\n", pWFHM->modeTag);

  fprintf(file,"\n*** Intrinsic Parameters ***\n");
  fprintf(file,"eta   = %.16f\n", pWF22->eta);
  fprintf(file,"chi1z = %.16f\n", pWF22->chi1L);
  fprintf(file,"chi2z = %.16f\n", pWF22->chi2L);

  fprintf(file,"\n*** QNM Frequencies ***\n");
  fprintf(file,"fRINGlm = %.16f\n", pWFHM->fRING);
  fprintf(file,"fDAMPlm = %.16f\n", pWFHM->fDAMP);


  fprintf(file,"\n*** Phase Inspiral ***\n");
  fprintf(file,"Lambda = %.16f\n", pPhase->LambdaPN);

  fprintf(file,"\n*** Phase Intermediate ***\n f_i\t dphi(f_i)\n");

  for(int i = 0; i<N_MAX_COEFFICIENTS_PHASE_INTER; i++)
      fprintf(file, "%.16f \t %.16f \n", pPhase->CollocationPointsFreqsPhaseInter[i], pPhase->IntermediatePhaseFits[pWFHM->modeInt*N_MAX_COEFFICIENTS_PHASE_INTER+i](pWF22,pWFHM->IMRPhenomXHMIntermediatePhaseFitsVersion));


  if(pWFHM->modeTag==32){
      fprintf(file,"\n*** Mixing Coefficients ***\n");
      fprintf(file,"222 = %.16f  %.16f *I\n", creal(pWFHM->mixingCoeffs[0]), cimag(pWFHM->mixingCoeffs[0]));
      fprintf(file,"223 = %.16f  %.16f *I\n", creal(pWFHM->mixingCoeffs[1]), cimag(pWFHM->mixingCoeffs[1]));
      fprintf(file,"322 = %.16f  %.16f *I\n", creal(pWFHM->mixingCoeffs[2]), cimag(pWFHM->mixingCoeffs[2]));
      fprintf(file,"323 = %.16f  %.16f *I\n", creal(pWFHM->mixingCoeffs[3]), cimag(pWFHM->mixingCoeffs[3]));

      fprintf(file, "\n*** Ringdown Phase Spheroidal ***\n f_i \t dphi(f_i)\n");
      for(int i=0; i<pWFHM->nCollocPtsRDPhase; i++)
          fprintf(file, "%.16f \t %.16f \n", pPhase->CollocationPointsFreqsPhaseRD[i],pPhase->RingdownPhaseFits[i](pWF22,pWFHM->IMRPhenomXHMRingdownPhaseVersion));
      fprintf(file,"Deltadphi=%.16f\n Deltaphi=%.16f \n", IMRPhenomXHM_RD_Phase_32_SpheroidalTimeShift(pWF22,pWFHM->IMRPhenomXHMRingdownPhaseVersion),IMRPhenomXHM_RD_Phase_32_SpheroidalPhaseShift(pWF22,pWFHM->IMRPhenomXHMRingdownPhaseVersion));

  }

  else{
      fprintf(file, "\n*** Ringdown phase ***\n");
      fprintf(file,"alpha2=%.16f\n alphaL=%.16f",pPhase->alpha2,pPhase->alphaL);
  }

  fclose(file);
  #endif

  }



/* Reconstruct the coefficients of the ringdown phase in a spheroidal-harmonics basis --only for the 32 mode    */
void  GetSpheroidalCoefficients(IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXPhaseCoefficients *pPhase22, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22){


    int nCollocationPts_RD_Phase = pWFHM->nCollocPtsRDPhase;
    double CollocValuesPhaseRingdown[nCollocationPts_RD_Phase];
    double CollocFreqsPhaseRingdown[nCollocationPts_RD_Phase];

    #if DEBUG == 1
    printf("Declare Colloc Values and Freqs:\n");
    #endif
    // define gsl variables
    gsl_vector *b, *x;
    gsl_matrix *A;
    gsl_permutation *p;
    int s;

    p = gsl_permutation_alloc(nCollocationPts_RD_Phase);
    b = gsl_vector_alloc(nCollocationPts_RD_Phase);
    x = gsl_vector_alloc(nCollocationPts_RD_Phase);
    A = gsl_matrix_alloc(nCollocationPts_RD_Phase,nCollocationPts_RD_Phase);

    IMRPhenomXHM_Ringdown_CollocPtsFreqs(pPhase,pWFHM, pWF22);

    #if DEBUG == 1
    printf("Initialize RD CollocPtsFreqs:\n");
    #endif

    if(pWFHM->IMRPhenomXHMRingdownPhaseVersion == 122019){
        // first fill-in the collocation points for the phase
        for(int i=0; i<nCollocationPts_RD_Phase; i++)
        {

          CollocValuesPhaseRingdown[i] =pPhase->RingdownPhaseFits[i](pWF22,pWFHM->IMRPhenomXHMRingdownPhaseFitsVersion);
          CollocFreqsPhaseRingdown[i] =pPhase->CollocationPointsFreqsPhaseRD[i];
          gsl_vector_set(b,i,CollocValuesPhaseRingdown[i]);
          REAL8 ff=CollocFreqsPhaseRingdown[i], ffm1=1./ff, ffm2=ffm1*ffm1;
          REAL8 fpowers[]={1., (pWFHM->fDAMP)/(pow(pWFHM->fDAMP,2)+pow(ff-(pWFHM->fRING),2)),ffm2,ffm2*ffm2};
          for(int j=0; j<nCollocationPts_RD_Phase; j++)
          gsl_matrix_set(A,i,j,fpowers[j]);

        }
    }
    else{
        for(int i = 0; i < nCollocationPts_RD_Phase; i++){
            CollocValuesPhaseRingdown[i] = pPhase->RingdownPhaseFits[i](pWF22,pWFHM->IMRPhenomXHMRingdownPhaseFitsVersion);
            CollocFreqsPhaseRingdown[i]  = pPhase->CollocationPointsFreqsPhaseRD[i];
            gsl_vector_set(b,i,CollocValuesPhaseRingdown[i]);
            REAL8 ff = CollocFreqsPhaseRingdown[i], ffm1 = 1./ff, ffm2 = ffm1 * ffm1;
            REAL8 lorentzian = pWFHM->fDAMP / (pWFHM->fDAMP * pWFHM->fDAMP + (ff - pWFHM->fRING) * (ff - pWFHM->fRING));
            REAL8 fpowers[] = {1., ffm1, ffm2, ffm2 * ffm2, lorentzian};
            for(int j = 0; j < nCollocationPts_RD_Phase; j++){
                gsl_matrix_set(A, i, j, fpowers[j]); // Ansatz spheroidal
            }
        }
    }

    #if DEBUG == 1
    printf("Collocation points in ringdown region:\n");
    for(int i=0; i<nCollocationPts_RD_Phase; i++){
      printf("p%d=(%.16e,%.16e)\n",i+1,CollocFreqsPhaseRingdown[i],CollocValuesPhaseRingdown[i]);
    }
    #endif

    // solve the linear system of Eq. (6.8)
    gsl_linalg_LU_decomp(A,p,&s);
    gsl_linalg_LU_solve(A,p,b,x);

    if(pWFHM->IMRPhenomXHMRingdownPhaseVersion == 122019){
        /* ansatz: alpha0 + (alpha2)/(f^2)+ (alpha4)/(f^4)  + alphaL*(fdamplm)/((fdamplm)^2 + (f - fRDlm)^2)*/
        pPhase->alpha0_S = gsl_vector_get(x,0);
        pPhase->alphaL_S = gsl_vector_get(x,1);
        pPhase->alpha2_S = gsl_vector_get(x,2);
        pPhase->alpha4_S = gsl_vector_get(x,3);
    }
    else{
        //printf("New Spheroidal Version. Getting coefficients\n");
        for(int i=0; i<nCollocationPts_RD_Phase; i++){
            pPhase->RDCoefficient[i] = gsl_vector_get(x, i);
        }
    }

    #if DEBUG == 1
    printf("**********\n alpha0_S=%.16e \n alphaL_S=%.16e \n alpha2_S=%.16e \n alpha4_S=%.16e \n \n ", pPhase->alpha0_S, pPhase->alphaL_S, pPhase->alpha2_S, pPhase->alpha4_S);
    #endif

    gsl_vector_free(b);
    gsl_vector_free(x);
    gsl_matrix_free(A);
    gsl_permutation_free(p);

    //  we have reconstructed the "shape" of the spheroidal ringdown phase derivative, dphiS, but we need to adjust the relative time and phase shift of the final phiS wrt IMRPhenomX

    IMRPhenomX_UsefulPowers powers_of_FREF;
    /*************** time-shift of spheroidal ansatz *******************/
    double frefRD=pWF22->fRING+pWF22->fDAMP;
    IMRPhenomX_Initialize_Powers(&powers_of_FREF,frefRD);
    // here we call a fit for dphiS(fref)-dphi22(fref)
    double tshift=IMRPhenomXHM_RD_Phase_32_SpheroidalTimeShift(pWF22,pWFHM->IMRPhenomXHMRingdownPhaseFitsVersion);

    // we compute dphi22(fref)
    IMRPhenomX_Phase_22_ConnectionCoefficients(pWF22,pPhase22);
    pWFHM->timeshift=IMRPhenomX_TimeShift_22(pPhase22, pWF22);
    
    pPhase->phi0_S = 0;

    // we impose that dphiS(fref)-dphi22(fref) has the value given by our fit
    if(pWFHM->IMRPhenomXHMRingdownPhaseVersion == 122019){
        double dphi22ref = 1./pWF22->eta*IMRPhenomX_dPhase_22(frefRD, &powers_of_FREF,pPhase22,pWF22)+pWFHM->timeshift;
        pPhase->alpha0_S = pPhase->alpha0_S +dphi22ref+tshift-IMRPhenomXHM_RD_Phase_Ansatz(frefRD,&powers_of_FREF,pWFHM,pPhase);
    }
    else{
        // REAL8 psi4tostrain=XLALSimIMRPhenomXPsi4ToStrain(pWF22->eta, pWF22->STotR, pWF22->dchi);
        // pPhase->RDCoefficient[0] -= 2.*LAL_PI*(500+psi4tostrain);
        double dphi22ref = 1./pWF22->eta*IMRPhenomX_dPhase_22(frefRD, &powers_of_FREF,pPhase22,pWF22)+pWFHM->timeshift;
        double dphi32ref = IMRPhenomXHM_RD_Phase_Ansatz(frefRD,&powers_of_FREF,pWFHM,pPhase);
        pPhase->RDCoefficient[0] += dphi22ref+tshift-dphi32ref;
        //printf("dphi32ref before time shift = %.16e\n", dphi32ref);
        // printf("fittimeshift = %.16e\n", tshift);
        // printf("Time shift = %.16e\n", dphi22ref+tshift-dphi32ref);
        // printf("dphi22ref = %.16e\n", dphi22ref);
        // printf("dphi32ref = %.16e\n", IMRPhenomXHM_RD_Phase_Ansatz(frefRD,&powers_of_FREF,pWFHM,pPhase));
        //pPhase->RDCoefficient[0] += (dphi22ref + tshift - IMRPhenomXHM_RD_Phase_Ansatz(frefRD, &powers_of_FREF, pWFHM, pPhase));
        //pPhase->RDCoefficient[5] -= 2.*LAL_PI*(500+psi4tostrain);
        INT4 exponent = 5; // a + b / f^5
        REAL8 ff = pWFHM->fRING + 2 * pWFHM->fDAMP;
        IMRPhenomX_UsefulPowers powers_of_ff;
        IMRPhenomX_Initialize_Powers(&powers_of_ff, ff);
        pPhase->RDCoefficient[nCollocationPts_RD_Phase + 1] = -IMRPhenomXHM_RD_Phase_DerAnsatz(ff, &powers_of_ff, pWFHM, pPhase) * pow(ff, exponent + 1) / exponent;
        pPhase->RDCoefficient[nCollocationPts_RD_Phase] = IMRPhenomXHM_RD_Phase_Ansatz(ff, &powers_of_ff, pWFHM, pPhase) - pPhase->RDCoefficient[nCollocationPts_RD_Phase + 1] / pow(ff, exponent);
        pPhase->RDCoefficient[nCollocationPts_RD_Phase + 2] = IMRPhenomXHM_RD_Phase_AnsatzInt(ff, &powers_of_ff, pWFHM, pPhase) - (pPhase->RDCoefficient[5]*ff - 0.25 * pPhase->RDCoefficient[6] * powers_of_ff.m_four);
        pWFHM->fPhaseRDflat = ff;   
    }

    /*************** phase-shift of spheroidal ansatz *******************/
    //if(pWFHM->IMRPhenomXHMRingdownPhaseVersion == 122019 || pWFHM->IMRPhenomXHMRingdownPhaseVersion == 20221114){
        frefRD = pWF22->fRING;
    // }
    // else{
    //     frefRD = pWF22->fRING - pWF22->fDAMP;
    //     //pWFHM->phaseshift += LAL_PI;
    // }
    IMRPhenomX_Initialize_Powers(&powers_of_FREF,frefRD);

    /* we compute phi22(fref) */
    IMRPhenomX_UsefulPowers powers_of_MfRef;
    IMRPhenomX_Initialize_Powers(&powers_of_MfRef,pWF22->MfRef);
    pWFHM->phiref22 = -1./pWF22->eta*IMRPhenomX_Phase_22(pWF22->MfRef, &powers_of_MfRef, pPhase22, pWF22) - pWFHM->timeshift*pWF22->MfRef - pWFHM->phaseshift + 2.0*pWF22->phi0 + LAL_PI_4;
    double phi22ref=1./pWF22->eta*IMRPhenomX_Phase_22(frefRD, &powers_of_FREF,pPhase22,pWF22) + pWFHM->timeshift*frefRD + pWFHM->phaseshift + pWFHM->phiref22;
    // printf("phi22ref = %.16e\n", phi22ref);
    // printf("phi32ref = %.16e\n", IMRPhenomXHM_RD_Phase_AnsatzInt(frefRD,&powers_of_FREF,pWFHM,pPhase));
    // printf("frefRD = %.16e\n", frefRD);
    
    // we call a fit for Mod[phiS(fref)-phi22(fref),2*Pi]
    double phishift=IMRPhenomXHM_RD_Phase_32_SpheroidalPhaseShift(pWF22,pWFHM->IMRPhenomXHMRingdownPhaseFitsVersion);
    //if (pWFHM->IMRPhenomXHMRingdownPhaseVersion == 20221114) phishift += LAL_PI;
    //printf("fitphaseshift = %.6f\n", phishift);

    //we adjust the relative phase of our reconstruction
    pPhase->phi0_S= phi22ref - IMRPhenomXHM_RD_Phase_AnsatzInt(frefRD,&powers_of_FREF,pWFHM,pPhase) + phishift;
    // printf("phaseshift = %.16e\n", pPhase->phi0_S);
    // 
    // printf("phi32ref = %.16e\n", IMRPhenomXHM_RD_Phase_AnsatzInt(frefRD,&powers_of_FREF,pWFHM,pPhase));
    
        
    #if DEBUG == 1
    printf("**********\n alpha0_S=%.16e \n alphaL_S=%.16e \n alpha2_S=%.16e \n alpha4_S=%.16e \n \n ", pPhase->alpha0_S, pPhase->alphaL_S, pPhase->alpha2_S, pPhase->alpha4_S);
    printf("**********\n **Spheroidal reconstruction parameters:**\n \n phi0_S=%.16e \n alpha0_S=%.16e \n alphaL_S=%.16e \n alpha2_S=%.16e \n alpha4_S=%.16e \n \n ",pPhase->phi0_S, pPhase->alpha0_S, pPhase->alphaL_S, pPhase->alpha2_S, pPhase->alpha4_S);
    #endif


  }


/**************************************/
/*                                    */
/*  Spheroidal -> Spherical rotation  */
/*                                    */
/**************************************/
/* The rotation consists of a linear transformation using the mixing coefficients given by Berti [10.1103/PhysRevD.90.064012].

  h32_spherical = a1 * h22_spheroidal + a2 * h32_spheroidal,  where a1, a2 are the mixing coefficients.

  Since the 22 is the most dominant mode, it will not show a significant mixing with any other mode,
  so we can assume that h22_spheroidal = h22_spherical
*/

  // In principle this could be generalized to the 43 mode: for the time being, assume the mode solved for is only the 32.
  COMPLEX16 SpheroidalToSpherical(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXAmpCoefficients *pAmp22, IMRPhenomXPhaseCoefficients *pPhase22, IMRPhenomXHMAmpCoefficients *pAmplm, IMRPhenomXHMPhaseCoefficients *pPhaselm, IMRPhenomXHMWaveformStruct *pWFlm, IMRPhenomXWaveformStruct *pWF22)
  {
    REAL8 Mf = powers_of_Mf->itself;
    // Compute the 22 mode using PhenomX functions. This gives the 22 mode rescaled with the leading order.  This is because the 32 is also rescaled.
    REAL8 amp22=XLALSimIMRPhenomXRingdownAmplitude22AnsatzAnalytical(Mf, pWF22->fRING, pWF22->fDAMP, pAmp22->gamma1, pAmp22->gamma2, pAmp22->gamma3);
    REAL8 phi22=1./pWF22->eta*IMRPhenomX_Phase_22(Mf, powers_of_Mf, pPhase22,pWF22) + pWFlm->timeshift*Mf + pWFlm->phaseshift + pWFlm->phiref22;
    COMPLEX16 wf22R = amp22 * cexp(I * phi22);
    if (pWFlm->IMRPhenomXHMRingdownAmpVersion != 0){
        wf22R *= pWFlm->ampNorm * powers_of_Mf->m_seven_sixths;
    }
    // Compute 32 mode in spheroidal.
    REAL8 amplm=IMRPhenomXHM_RD_Amp_Ansatz(powers_of_Mf, pWFlm, pAmplm);
    REAL8 philm=IMRPhenomXHM_RD_Phase_AnsatzInt(Mf, powers_of_Mf,pWFlm, pPhaselm);
    // Do the rotation.
    COMPLEX16 sphericalWF_32=conj(pWFlm->mixingCoeffs[2]) * wf22R + conj(pWFlm->mixingCoeffs[3])*amplm*cexp(I*philm);
    return sphericalWF_32;

  }

  // If the 22 mode has been previously computed, we use it here for the rotation.
  COMPLEX16 SpheroidalToSphericalRecycle(IMRPhenomX_UsefulPowers *powers_of_Mf, COMPLEX16 wf22, IMRPhenomXHMAmpCoefficients *pAmplm, IMRPhenomXHMPhaseCoefficients *pPhaselm, IMRPhenomXHMWaveformStruct *pWFlm)
  {
    // The input 22 in the whole 22, and for the rotation we have to rescaled with the leading order. This is because the 32 is also rescaled.
    //complex double wf22R = wf22/(powers_of_f->m_seven_sixths * pWF22->amp0);
    COMPLEX16 wf22R = wf22 / pWFlm->Amp0;
    if (pWFlm->IMRPhenomXHMRingdownAmpVersion == 0){
        wf22R /= (pWFlm->ampNorm * powers_of_Mf->m_seven_sixths);
    }
    // Compute 32 mode in spheroidal.
    REAL8 amplm=IMRPhenomXHM_RD_Amp_Ansatz(powers_of_Mf, pWFlm, pAmplm);
    REAL8 philm=IMRPhenomXHM_RD_Phase_AnsatzInt(powers_of_Mf->itself, powers_of_Mf, pWFlm, pPhaselm);
    // Do the rotation
    COMPLEX16 sphericalWF_32 = conj(pWFlm->mixingCoeffs[2])*wf22R  +  conj(pWFlm->mixingCoeffs[3])*amplm*cexp(I*philm);
    return sphericalWF_32;

  }


  /**************************************/
  /*                                    */
  /*  RECONSTRUCTION FUNCTIONS THROUGH  */
  /*     INSPIRAL, MERGER & RINGDOWN    */
  /*                                    */
  /**************************************/

  /* These functions return the amplitude/phase for a given input frequency.
     They are piece-wise functions that distiguish between inspiral, intermediate and ringdown
     and call the corresponding reconstruction function.
  */

  /*******************************************************************************/
  /*  Compute IMRPhenomXHM AMPLITUDE given an amplitude coefficients struct pAmp */
  /*******************************************************************************/

  // WITHOUT mode mixing. It returns the whole amplitude (in NR units) without the normalization factor of the 22: sqrt[2 * eta / (3 * pi^(1/3))]
  REAL8 IMRPhenomXHM_Amplitude_noModeMixing(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomXHMWaveformStruct *pWF) {
    // If it is an odd mode and equal black holes case this mode is zero.
    if(pWF->Ampzero==1){
      return 0.;
    }
    
    // Use step function to only calculate IMR regions in approrpiate frequency regime
    REAL8 Amp, Mf = powers_of_Mf->itself;
    // Inspiral range
    if (!IMRPhenomX_StepFuncBool(Mf, pAmp->fAmpMatchIN)){
      Amp = IMRPhenomXHM_Inspiral_Amp_Ansatz(powers_of_Mf, pWF, pAmp);
    }
    // MRD range
    else if (IMRPhenomX_StepFuncBool(Mf, pAmp->fAmpMatchIM)){
      // if (pWFHM->fAmpRDfalloff > 0 && IMRPhenomX_StepFuncBool(Mf, pWFHM->fAmpRDfalloff)){
      //     Amp = pAmp->RDCoefficient[3] * exp(- pAmp->RDCoefficient[4] * (Mf - pWFHM->fAmpRDfalloff));
      // }
      //else{
          Amp = IMRPhenomXHM_RD_Amp_Ansatz(powers_of_Mf, pWF, pAmp);
          if (pWF->IMRPhenomXHMRingdownAmpVersion == 0) Amp *= powers_of_Mf->m_seven_sixths * pWF->ampNorm;
     //}
    }
    /* Intermediate range */
    // First intermediate region.
    else if ((pWF->AmpEMR==1) && !IMRPhenomX_StepFuncBool(Mf, pAmp->fAmpMatchInt12)){
      INT4 tmp = pAmp->InterAmpPolOrder;
      pAmp->InterAmpPolOrder = 1042;
      Amp = IMRPhenomXHM_Intermediate_Amp_Ansatz(powers_of_Mf, pWF, pAmp);
      pAmp->InterAmpPolOrder = tmp;
    }
    else{  //Second intermediate region
        Amp = IMRPhenomXHM_Intermediate_Amp_Ansatz(powers_of_Mf, pWF, pAmp);
    }
    if (Amp < 0 && pWF->IMRPhenomXHMReleaseVersion != 122019) Amp = FALSE_ZERO;
    return Amp;
  }

  // WITH mode mixing. It returns the whole amplitude (in NR units) without the normalization factor of the 22: sqrt[2 * eta / (3 * pi^(1/3))].
  REAL8 IMRPhenomXHM_Amplitude_ModeMixing(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXHMWaveformStruct *pWF, IMRPhenomXAmpCoefficients *pAmp22, IMRPhenomXPhaseCoefficients *pPhase22, IMRPhenomXWaveformStruct *pWF22) {
    // Use step function to only calculate IMR regions in approrpiate frequency regime
    REAL8 Amp, Mf = powers_of_Mf->itself;
    // Inspiral range
    if (!IMRPhenomX_StepFuncBool(Mf, pAmp->fAmpMatchIN)){
      Amp =  IMRPhenomXHM_Inspiral_Amp_Ansatz(powers_of_Mf, pWF, pAmp);
    }
    // MRD range
    else if (IMRPhenomX_StepFuncBool(Mf, pAmp->fAmpMatchIM)){
      Amp = cabs(SpheroidalToSpherical(powers_of_Mf, pAmp22, pPhase22, pAmp, pPhase, pWF, pWF22));
      if (pWF->IMRPhenomXHMRingdownAmpVersion == 0) Amp *= powers_of_Mf->m_seven_sixths * pWF->ampNorm;
    }
    /* Intermediate range */
    // First intermediate region
    else if ((pWF->AmpEMR==1) && !IMRPhenomX_StepFuncBool(Mf, pAmp->fAmpMatchInt12)){
      INT4 tmp = pAmp->InterAmpPolOrder;
      pAmp->InterAmpPolOrder = 1042;
      Amp = IMRPhenomXHM_Intermediate_Amp_Ansatz(powers_of_Mf, pWF, pAmp);
      pAmp->InterAmpPolOrder = tmp;
    }
    //Second intermediate region
    else{
        Amp = IMRPhenomXHM_Intermediate_Amp_Ansatz(powers_of_Mf, pWF, pAmp);
    }
    if (Amp < 0 && pWF->IMRPhenomXHMReleaseVersion != 122019) Amp = FALSE_ZERO;
    return Amp;
  }

  // WITH mode mixing and recycling the previously computed 22 mode. It returns the whole amplitude (in NR units) without the normalization factor of the 22: sqrt[2 * eta / (3 * pi^(1/3))].
  REAL8 IMRPhenomXHM_Amplitude_ModeMixingRecycle(IMRPhenomX_UsefulPowers *powers_of_Mf, COMPLEX16 wf22, IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXHMWaveformStruct *pWF) {
    // Use step function to only calculate IMR regions in approrpiate frequency regime
    REAL8 Amp, Mf = powers_of_Mf->itself;
    // Inspiral range
    if (!IMRPhenomX_StepFuncBool(Mf, pAmp->fAmpMatchIN)){
      Amp =  IMRPhenomXHM_Inspiral_Amp_Ansatz(powers_of_Mf, pWF, pAmp);
    }
    // MRD range
    else if (IMRPhenomX_StepFuncBool(Mf, pAmp->fAmpMatchIM)){
      Amp = cabs(SpheroidalToSphericalRecycle(powers_of_Mf, wf22, pAmp, pPhase, pWF));
      if (pWF->IMRPhenomXHMRingdownAmpVersion == 0) Amp *= powers_of_Mf->m_seven_sixths * pWF->ampNorm;
    }
    /* Intermediate range */
    // First intermediate region
    else if ((pWF->AmpEMR==1) && !IMRPhenomX_StepFuncBool(Mf, pAmp->fAmpMatchInt12)){ // FIXME: why was the second condition commented?
      INT4 tmp = pAmp->InterAmpPolOrder;
      pAmp->InterAmpPolOrder = 1042;
      Amp = IMRPhenomXHM_Intermediate_Amp_Ansatz(powers_of_Mf, pWF, pAmp);
      pAmp->InterAmpPolOrder = tmp;
    }
    //Second intermediate region
    else{
        Amp = IMRPhenomXHM_Intermediate_Amp_Ansatz(powers_of_Mf, pWF, pAmp);
    }
    if (Amp < 0 && pWF->IMRPhenomXHMReleaseVersion != 122019) Amp = FALSE_ZERO;
    return Amp;
  }


  /*******************************************************************************/
  /*  Compute IMRPhenomXHM PHASE given a phase coefficients struct pPhase        */
  /*******************************************************************************/

// WITHOUT mode mixing.
REAL8 IMRPhenomXHM_Phase_noModeMixing(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXHMWaveformStruct *pWF, UNUSED IMRPhenomXWaveformStruct *pWF22)
{
  REAL8 Mf = powers_of_Mf->itself;
  // Inspiral range, f < fPhaseInsMax
  if (!IMRPhenomX_StepFuncBool(Mf, pPhase->fPhaseMatchIN))
  {
      REAL8 PhiIns = IMRPhenomXHM_Inspiral_Phase_AnsatzInt(Mf, powers_of_Mf, pPhase);
      return PhiIns + pPhase->C1INSP*Mf + pPhase->CINSP + pPhase->deltaphiLM;
  }
  // MRD range, f > fPhaseIntMax
  if (IMRPhenomX_StepFuncBool(Mf, pPhase->fPhaseMatchIM))
  {
      REAL8 PhiMRD = IMRPhenomXHM_RD_Phase_AnsatzInt(Mf, powers_of_Mf, pWF, pPhase);
      return PhiMRD + pPhase->C1RD*Mf + pPhase->CRD + pPhase->deltaphiLM;
  }
  //Intermediate range, fPhaseInsMax < f < fPhaseIntMax
  REAL8 PhiInt = IMRPhenomXHM_Inter_Phase_AnsatzInt(Mf, powers_of_Mf, pWF, pPhase);
  return PhiInt + pPhase->deltaphiLM;
}

// WITH mode mixing.
REAL8 IMRPhenomXHM_Phase_ModeMixing(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXHMWaveformStruct *pWF,  IMRPhenomXAmpCoefficients *pAmp22, IMRPhenomXPhaseCoefficients *pPhase22, IMRPhenomXWaveformStruct *pWF22)
{
  REAL8 Mf = powers_of_Mf->itself;
  // Inspiral range, f < fPhaseInsMax
  if (!IMRPhenomX_StepFuncBool(Mf, pPhase->fPhaseMatchIN))
  {
      REAL8 PhiIns = IMRPhenomXHM_Inspiral_Phase_AnsatzInt(Mf, powers_of_Mf, pPhase);
      return PhiIns + pPhase->C1INSP*Mf + pPhase->CINSP + pPhase->deltaphiLM;
  }
  // MRD range, f > fPhaseIntMax
  if (IMRPhenomX_StepFuncBool(Mf, pPhase->fPhaseMatchIM))
  {
      REAL8 PhiMRD = carg(SpheroidalToSpherical(powers_of_Mf, pAmp22, pPhase22, pAmp, pPhase, pWF, pWF22));
      return PhiMRD + pPhase->C1RD*Mf + pPhase->CRD + pPhase->deltaphiLM;
  }
  //Intermediate range, fPhaseInsMax < f < fPhaseIntMax
  REAL8 PhiInt = IMRPhenomXHM_Inter_Phase_AnsatzInt(Mf, powers_of_Mf, pWF, pPhase);
  return PhiInt + pPhase->deltaphiLM;
}

// WITH mode mixing and recycling the previously computed 22 mode.
REAL8 IMRPhenomXHM_Phase_ModeMixingRecycle(IMRPhenomX_UsefulPowers *powers_of_Mf, COMPLEX16 wf22, IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXHMWaveformStruct *pWF)
{
  REAL8 Mf = powers_of_Mf->itself;
  // Inspiral range, f < fPhaseInsMax
  if (!IMRPhenomX_StepFuncBool(Mf, pPhase->fPhaseMatchIN))
  {
      REAL8 PhiIns = IMRPhenomXHM_Inspiral_Phase_AnsatzInt(Mf, powers_of_Mf, pPhase);
      return PhiIns + pPhase->C1INSP*Mf + pPhase->CINSP + pPhase->deltaphiLM;
  }
  // MRD range, f > fPhaseIntMax
  if (IMRPhenomX_StepFuncBool(Mf, pPhase->fPhaseMatchIM))
  {
      REAL8 PhiMRD = carg(SpheroidalToSphericalRecycle(powers_of_Mf, wf22, pAmp, pPhase, pWF));
      return PhiMRD + pPhase->C1RD*Mf + pPhase->CRD + pPhase->deltaphiLM;
  }
  //Intermediate range, fPhaseInsMax < f < fPhaseIntMax
  REAL8 PhiInt = IMRPhenomXHM_Inter_Phase_AnsatzInt(Mf, powers_of_Mf, pWF, pPhase);
  return PhiInt + pPhase->deltaphiLM;
}


  // WITHOUT mode mixing.
  double IMRPhenomXHM_dPhase_noModeMixing(double f, IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXHMWaveformStruct *pWF, UNUSED IMRPhenomXWaveformStruct *pWF22)
  {
    // Inspiral range, f < fPhaseInsMax
    if (!IMRPhenomX_StepFuncBool(f, pPhase->fPhaseMatchIN))
    {
      double dPhiIns = IMRPhenomXHM_Inspiral_Phase_Ansatz(f, powers_of_f, pPhase);
      return dPhiIns + pPhase->C1INSP;
    }
    // MRD range, f > fPhaseIntMax
    if (IMRPhenomX_StepFuncBool(f, pPhase->fPhaseMatchIM))
    {
      double dPhiMRD = IMRPhenomXHM_RD_Phase_Ansatz(f, powers_of_f, pWF, pPhase);
      return dPhiMRD + pPhase->C1RD;
    }
    //Intermediate range, fPhaseInsMax < f < fPhaseIntMax
    double dPhiInt = IMRPhenomXHM_Inter_Phase_Ansatz(f, powers_of_f, pWF, pPhase);
    return dPhiInt;
  }

  // WITH mode mixing.
  double IMRPhenomXHM_dPhase_ModeMixing(double f, IMRPhenomX_UsefulPowers *powers_of_f,IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXHMWaveformStruct *pWF,  IMRPhenomXAmpCoefficients *pAmp22, IMRPhenomXPhaseCoefficients *pPhase22, IMRPhenomXWaveformStruct *pWF22)
  {
    // Inspiral range, f < fPhaseInsMax
    if (!IMRPhenomX_StepFuncBool(f, pPhase->fPhaseMatchIN))
    {
      double dPhiIns = IMRPhenomXHM_Inspiral_Phase_Ansatz(f, powers_of_f, pPhase);
      return dPhiIns + pPhase->C1INSP;
    }
    // MRD range, f > fPhaseIntMax
    if (IMRPhenomX_StepFuncBool(f, pPhase->fPhaseMatchIM))
    {
      double dPhiMRD = carg(SpheroidalToSpherical(powers_of_f, pAmp22, pPhase22, pAmp, pPhase, pWF, pWF22));
      return dPhiMRD + pPhase->C1RD;
    }
    //Intermediate range, fPhaseInsMax < f < fPhaseIntMax
    double dPhiInt = IMRPhenomXHM_Inter_Phase_Ansatz(f, powers_of_f, pWF, pPhase);
    return dPhiInt;
  }

  

  /*****************/
  /*   DEBUGGING   */
  /*****************/

  // This is just for debugging. It prints some interesting parameters to a file.
  int ParametersToFile(
    IMRPhenomXWaveformStruct *pWF,                /**< Wf structure for the 22 mode*/
    IMRPhenomXHMWaveformStruct *pWFHM,            /**< Wf structure for the lm mode*/
    IMRPhenomXHMAmpCoefficients *pAmp,            /**< Coefficients struct of the lm Amplitude */
    UNUSED IMRPhenomXHMPhaseCoefficients *pPhase  /**< Coefficients struct of the lm Phase */
  )
  {

    FILE *file;
    char fileSpec[40];
    sprintf(fileSpec, "Parameters%i.dat", pWFHM->modeTag);
    printf("\nOutput Parameter file: %s\r\n",fileSpec);
    file = fopen(fileSpec,"w");

    fprintf(file,"\n*** %i Mode ***\n", pWFHM->modeTag);

    fprintf(file,"\n*** Intrinsic Parameters ***\n");
    fprintf(file,"eta   = %.16f\n", pWF->eta);
    fprintf(file,"chi1z = %.16f\n", pWF->chi1L);
    fprintf(file,"chi2z = %.16f\n", pWF->chi2L);

    fprintf(file,"\n*** Extrinsic Parameters ***\n");
    fprintf(file,"distance    = %.16f\n", pWF->distance);
    fprintf(file,"inclination = %.16f\n", pWF->inclination);
    fprintf(file,"phiRef      = %.16f\n", pWF->phifRef);
    fprintf(file,"fRef        = %.16f\n", pWF->fRef);
    fprintf(file,"m1_SI       = %.16f\n", pWF->m1_SI);
    fprintf(file,"m2_SI       = %.16f\n", pWF->m2_SI);

    fprintf(file,"\n*** Frequency Grid ***\n");
    fprintf(file,"deltaF = %.16f\n", pWF->deltaF);
    fprintf(file,"fMin   = %.16f\n", pWF->fMin);
    fprintf(file,"fMax   = %.16f\n", pWF->fMax);

    fprintf(file,"\n*** QNM Frequencies ***\n");
    fprintf(file,"fRING22 = %.16f\n", pWF->fRING);
    fprintf(file,"fDAMP22 = %.16f\n", pWF->fDAMP);
    fprintf(file,"fRINGlm = %.16f\n", pWFHM->fRING);
    fprintf(file,"fDAMPlm = %.16f\n", pWFHM->fDAMP);

    fprintf(file,"\n*** Amplitude Frequencies ***\n");
    fprintf(file,"fInsp1         = %.16f\n", pAmp->CollocationPointsFreqsAmplitudeInsp[0]);
    fprintf(file,"fInsp2         = %.16f\n", pAmp->CollocationPointsFreqsAmplitudeInsp[1]);
    fprintf(file,"fInsp3         = %.16f\n", pAmp->CollocationPointsFreqsAmplitudeInsp[2]);
    fprintf(file,"fAmpMatchIN    = %.16f\n", pAmp->fAmpMatchIN);
    fprintf(file,"fAmpMatchInt12 = %.16f\n", pAmp->fAmpMatchInt12);
    fprintf(file,"fInt1          = %.16f\n", pAmp->CollocationPointsFreqsAmplitudeInter[0]);
    fprintf(file,"fInt2          = %.16f\n", pAmp->CollocationPointsFreqsAmplitudeInter[1]);
    fprintf(file,"fAmpMatchIM    = %.16f\n", pAmp->fAmpMatchIM);

    fprintf(file,"\n*** Amplitude Collocation Points ***\n");
    fprintf(file,"A_fInsp1 = %.16f\n", pAmp->CollocationPointsValuesAmplitudeInsp[0]);
    fprintf(file,"A_fInsp2 = %.16f\n", pAmp->CollocationPointsValuesAmplitudeInsp[1]);
    fprintf(file,"A_fInsp3 = %.16f\n", pAmp->CollocationPointsValuesAmplitudeInsp[2]);
    fprintf(file,"A_fInt1  = %.16f\n", pAmp->CollocationPointsValuesAmplitudeInter[0]);
    fprintf(file,"A_fInt2  = %.16f\n", pAmp->CollocationPointsValuesAmplitudeInter[1]);

    fprintf(file,"\n*** Amplitude Coefficients ***\n");
    fprintf(file,"rho1    = %.16f\n", pAmp->rho1); //Inspiral coefficients
    fprintf(file,"rho2    = %.16f\n", pAmp->rho2);
    fprintf(file,"rho3    = %.16f\n", pAmp->rho3);
    fprintf(file,"delta0  = %.16f\n", pAmp->delta0); //Intermediate region
    fprintf(file,"delta1  = %.16f\n", pAmp->delta1);
    fprintf(file,"delta2  = %.16f\n", pAmp->delta2);
    fprintf(file,"delta3  = %.16f\n", pAmp->delta3);
    fprintf(file,"delta4  = %.16f\n", pAmp->delta4);
    fprintf(file,"delta5  = %.16f\n", pAmp->delta5);
    fprintf(file,"alambda = %.20f\n", pAmp->RDCoefficient[0]); //Ringdown region
    fprintf(file,"lambda  = %.16f\n", pAmp->RDCoefficient[1]);
    fprintf(file,"sigma   = %.16f\n", pAmp->RDCoefficient[2]);
    fprintf(file,"lc      = %.16f\n", pAmp->RDCoefficient[3]);
    fprintf(file,"alpha0  = %.16f\n", pAmp->alpha0); //First Intermediate region (only for EMR)
    fprintf(file,"alpha1  = %.16f\n", pAmp->alpha1);
    fprintf(file,"alpha2  = %.16f\n", pAmp->alpha2);
    fprintf(file,"alpha3  = %.16f\n", pAmp->alpha3);
    fprintf(file,"alpha4  = %.16f\n", pAmp->alpha4);

    fprintf(file,"\n*** Amplitude Intermediate Input ***\n");
    fprintf(file,"f1 = %.16f\n", pAmp->f1);
    fprintf(file,"f2 = %.16f\n", pAmp->f2);
    fprintf(file,"f3 = %.16f\n", pAmp->f3);
    fprintf(file,"f4 = %.16f\n", pAmp->f4);
    fprintf(file,"v1 = %.16f\n", pAmp->v1);
    fprintf(file,"v2 = %.16f\n", pAmp->v2);
    fprintf(file,"v3 = %.16f\n", pAmp->v3);
    fprintf(file,"v4 = %.16f\n", pAmp->v4);
    fprintf(file,"d1 = %.16f\n", pAmp->d1);
    fprintf(file,"d4 = %.16f\n", pAmp->d4);

    fprintf(file,"\n*** PN Amplitude Inspiral ***\n");
    fprintf(file,"PN_f1 = %.16f\n", pAmp->PNAmplitudeInsp[0]);
    fprintf(file,"PN_f2 = %.16f\n", pAmp->PNAmplitudeInsp[1]);
    fprintf(file,"PN_f3 = %.16f\n", pAmp->PNAmplitudeInsp[2]);

    fprintf(file,"\n*** Amplitude Versions ***\n");
    fprintf(file,"InspiralFits      = %i\n", pWFHM->IMRPhenomXHMInspiralAmpFitsVersion);
    fprintf(file,"IntermediateFits  = %i\n", pWFHM->IMRPhenomXHMIntermediateAmpFitsVersion);
    fprintf(file,"RDFits            = %i\n", pWFHM->IMRPhenomXHMRingdownAmpFitsVersion);
    fprintf(file,"InspiralRecon     = %i\n", pWFHM->IMRPhenomXHMInspiralAmpVersion);
    fprintf(file,"IntermediateRecon = %i\n", pWFHM->IMRPhenomXHMIntermediateAmpVersion);
    fprintf(file,"RDRecon           = %i\n", pWFHM->IMRPhenomXHMRingdownAmpVersion);

    fprintf(file,"\n*** Amplitude Extra ***\n");
    fprintf(file,"AmpEMR  = %i\n", pWFHM->AmpEMR);
    fprintf(file,"Ampzero = %i\n", pWFHM->Ampzero);

    if(pWFHM->modeTag==32){
      fprintf(file,"\n*** Mixing Coefficients ***\n");
      fprintf(file,"222 = %.16f %.16f\n", creal(pWFHM->mixingCoeffs[0]), cimag(pWFHM->mixingCoeffs[0]));
      fprintf(file,"223 = %.16f %.16f\n", creal(pWFHM->mixingCoeffs[1]), cimag(pWFHM->mixingCoeffs[1]));
      fprintf(file,"322 = %.16f %.16f\n", creal(pWFHM->mixingCoeffs[2]), cimag(pWFHM->mixingCoeffs[2]));
      fprintf(file,"323 = %.16f %.16f\n", creal(pWFHM->mixingCoeffs[3]), cimag(pWFHM->mixingCoeffs[3]));
    }

    fclose(file);

    return 0;
  }
