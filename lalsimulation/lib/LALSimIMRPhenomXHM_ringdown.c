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
 *
 */

#include "LALSimIMRPhenomXHM_ringdown.h"

/* Include fits of the Ringdown related quantities and the ansatz for amplitude and phase */

/* These parameter-space fits are documented in the supplementary material of https://dcc.ligo.org/P2000011-v2.
   There are 2 Mathematica notebooks (one for amplitude and one for phase) that read the fits data and automatically generate the C-code below.
   For more information read https://git.ligo.org/waveforms/reviews/imrphenomx/blob/master/documentation/ParspaceFits/README and the documentation in the notebooks. */


/****************************************/
/*                                      */
/*              AMPLITUDE               */
/*                                      */
/****************************************/

/* Fits over parameter space for the ringdown amplitude coefficients (alambda, lambda, sigma) for each mode. */

// The spin parameter S = (m1^2*chi1 + m2^2*chi2)/(m1^2 + m2^2)


// alambda, lambda and sigma are the coefficients of the ringdown ansatz, see Eq. (6.2)


static double IMRPhenomXHM_RD_Amp_21_alambda(double eta, double S, double chi1, double chi2, int RDAmpFlag) {
    UNUSED double total=0, delta=sqrt(1. - 4.*eta),eta2,eta3,eta4,S2,S3;
    switch (RDAmpFlag){
        case 122018:{
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            S2 = pow(S,2);
            S3 = pow(S,3);
            double noSpin = sqrt(eta - 4.*eta2)*(0.00734983387668636 - 0.0012619735607202085*eta + 0.01042318959002753*eta2);
            double eqSpin = sqrt(eta - 4.*eta2)*S*(-0.004839645742570202 - 0.0013927779195756036*S + eta2*(-0.054621206928483663 + 0.025956604949552205*S + 0.020360826886107204*S2));
            double uneqSpin = -0.018115657394753674*(chi1 - 1.*chi2)*eta2*(1. - 10.539795474715346*eta2*delta);
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_21_a: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_RD_Amp_21_lambda(double eta, double S, double chi1, double chi2, int RDAmpFlag) {
    UNUSED double total=0, delta=sqrt(1. - 4.*eta),eta2,S2;
    switch (RDAmpFlag){
        case 122018:{
            eta2 = pow(eta,2);
            S2 = pow(S,2);
            double noSpin = 0.5566284518926176 + 0.12651770333481904*eta + 1.8084545267208734*eta2;
            double eqSpin = (0.29074922226651545 + eta2*(-2.101111399437034 - 3.4969956644617946*S) + eta*(0.059317243606471406 - 0.31924748117518226*S) + 0.27420263462336675*S)*S;
            double uneqSpin = 1.0122975748481835*(chi1 - 1.*chi2)*eta2*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_21_lambda: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_RD_Amp_33_alambda(double eta, double S, double chi1, double chi2, int RDAmpFlag) {
  UNUSED double total=0,eta2,eta3,eta4,eta5,eta6,S2,delta=sqrt(1-4*eta);
  switch (RDAmpFlag){
    case 122018:{
      eta2 = pow(eta,2);
      eta3 = pow(eta,3);
      eta4 = pow(eta,4);
      eta5 = pow(eta,5);
      eta6 = pow(eta,6);
      S2 = pow(S,2);
      double noSpin = sqrt(eta - 4.*eta2)*(0.013700854227665184 + 0.01202732427321774*eta + 0.0898095508889557*eta2);
      double eqSpin = sqrt(eta - 4.*eta2)*(0.0075858980586079065 + eta*(-0.013132320758494439 - 0.018186317026076343*S) + 0.0035617441651710473*S)*S;
      double uneqSpin = eta4*(chi2*(-0.09802218411554885 - 0.05745949361626237*S) + chi1*(0.09802218411554885 + 0.05745949361626237*S) + eta2*(chi1*(-4.2679864481479886 - 11.877399902871485*S) + chi2*(4.2679864481479886 + 11.877399902871485*S))*delta);
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_33_alambda: version is not valid. Recommended version is 122018.");}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_33_lambda(double eta, double S, double chi1, double chi2, int RDAmpFlag) {
  UNUSED double total=0,eta2,S2,S3,delta=sqrt(1-4*eta);
  switch (RDAmpFlag){
    case 122018:{
      eta2 = pow(eta,2);
      S2 = pow(S,2);
      S3 = pow(S,3);
      double noSpin = 0.7435306475478924 - 0.06688558533374556*eta + 1.471989765837694*eta2;
      double eqSpin = S*(0.19457194111990656 + 0.07564220573555203*S + eta*(-0.4809350398289311 + 0.17261430318577403*S - 0.1988991467974821*S2));
      double uneqSpin = 1.8881959341735146*(chi1 - 1.*chi2)*eta2*delta;
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_33_lambda: version is not valid. Recommended version is 122018.");}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_32_alambda(double eta, double S, double chi1, double chi2, int RDAmpFlag) {
  UNUSED double total=0,eta2,eta3,eta4,eta5,eta6,eta7,eta8,eta9,S2;
  switch (RDAmpFlag){
    case 122018:{
      eta2 = pow(eta,2);
      eta3 = pow(eta,3);
      eta4 = pow(eta,4);
      eta5 = pow(eta,5);
      eta6 = pow(eta,6);
      eta7 = pow(eta,7);
      eta8 = pow(eta,8);
      eta9 = pow(eta,9);
      S2 = pow(S,2);
      double noSpin = 0.00012587900257140724 + 0.03927886286971654*eta - 0.8109309606583066*eta2 + 8.820604907164254*eta3 - 51.43344812454074*eta4 + 141.81940900657446*eta5 - 140.0426973304466*eta6;
      double eqSpin = S*(-0.00006001471234796344 + eta4*(-0.7849112300598181 - 2.09188976953315*S) + eta2*(0.08311497969032984 - 0.15569578955822236*S) + eta*(-0.01083175709906557 + 0.00568899459837252*S) - 0.00009363591928190229*S + 1.0670798489407887*eta3*S);
      double uneqSpin = -0.04537308968659669*pow(chi1 - 1.*chi2,2)*eta2*(1. - 8.711096029480697*eta + 18.362371966229926*eta2) + (chi1 - 1.*chi2)*(-297.36978685672733 + 3103.2516759087644*eta - 10001.774055779177*eta2 + 9386.734883473799*eta3)*eta6;
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_32_alambda: version is not valid. Recommended version is 122018.");}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_32_lambda(double eta, double S, double chi1, double chi2, int RDAmpFlag) {
  UNUSED double total=0,eta2,eta3,S2,delta=sqrt(1.-4*eta);
  switch (RDAmpFlag){
    case 122018:{
      eta2 = pow(eta,2);
      eta3 = pow(eta,3);
      S2 = pow(S,2);
      double noSpin = (sqrt(1. - 3.*eta)*(0.0341611244787871 - 0.3197209728114808*eta + 0.7689553234961991*eta2))/(0.048429644168112324 - 0.43758296068790314*eta + eta2);
      double eqSpin = sqrt(1. - 3.*eta)*S*(0.11057199932233873 + eta2*(25.536336676250748 - 71.18182757443142*S) + 9.790509295728649*eta*S + eta3*(-56.96407763839491 + 175.47259563543165*S));
      double uneqSpin = -5.002106168893265*pow(chi1 - 1.*chi2,2)*eta2*delta;
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_32_lambda: version is not valid. Recommended version is 122018.");}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_44_alambda(double eta, double S, double chi1, double chi2, int RDAmpFlag) {
    UNUSED double total=0, delta=sqrt(1. - 4.*eta),eta2,eta3,eta4,eta5,eta6,S2,S3;
    switch (RDAmpFlag){
        case 122018:{
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            S2 = pow(S,2);
            S3 = pow(S,3);
            double noSpin = sqrt(eta - 3.*eta2)*(0.007904587819112173 + 0.09558474985614368*eta - 2.663803397359775*eta2 + 28.298192768381554*eta3 - 136.10446022757958*eta4 + 233.23167528016833*eta5);
            double eqSpin = sqrt(eta - 3.*eta2)*S*(0.0049703757209330025 + 0.004122811292229324*S + eta*(-0.06166686913913691 + 0.014107365722576927*S)*S + eta2*(-0.2945455034809188 + 0.4139026619690879*S - 0.1389170612199015*S2) + eta3*(0.9225758392294605 - 0.9656098473922222*S + 0.19708289555425246*S2) + 0.000657528128497184*S2);
            double uneqSpin = 0.00659873279539475*(chi1 - 1.*chi2)*eta2*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_44_a: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_RD_Amp_44_lambda(double eta, double S, double chi1, double chi2, int RDAmpFlag) {
    UNUSED double total=0, delta=sqrt(1. - 4.*eta),eta2,eta3,eta4,eta5,eta6,eta7;
    switch (RDAmpFlag){
        case 122018:{
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            eta7 = pow(eta,7);
            double noSpin = 0.7702864948772887 + 32.81532373698395*eta - 1625.1795901450212*eta2 + 31305.458876573215*eta3 - 297375.5347399236*eta4 + 1.4726521941846698e6*eta5 - 3.616582470072637e6*eta6 + 3.4585865843680725e6*eta7;
            double eqSpin = (-0.03011582009308575*S + 0.09034746027925727*eta*S + 1.8738784391649446*eta2*S - 5.621635317494836*eta3*S)/(-1.1340218677260014 + S);
            double uneqSpin = 0.959943270591552*pow(chi1 - 1.*chi2,2)*eta2 + 0.853573071529436*(chi1 - 1.*chi2)*eta2*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_44_lambda: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_RD_Amp_21_sigma(double eta, double S, double chi1, double chi2, int RDAmpFlag) {
    UNUSED double total=0, delta=sqrt(1. - 4.*eta),eta2;
    switch (RDAmpFlag){
        case 122018:{
            eta2 = pow(eta,2);
            double noSpin = 1.2922261617161441 + 0.0019318405961363861*eta;
            double eqSpin = (0.04927982551108649 - 0.6703778360948937*eta + 2.6625014134659772*eta2)*S;
            double uneqSpin = 1.2001101665670462*(chi1 + pow(chi1,2) - 2.*chi1*chi2 + (-1. + chi2)*chi2)*eta2*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_21_sigma: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_RD_Amp_33_sigma(UNUSED double eta, UNUSED double S, UNUSED double chi1, UNUSED double chi2, int RDAmpFlag) {
    UNUSED double total=0;
    switch (RDAmpFlag){
        case 122018:{
            double noSpin = 1.3;
            double eqSpin = 0.;
            double uneqSpin = 0.;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_33_sigma: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_RD_Amp_32_sigma(UNUSED double eta, UNUSED double S, UNUSED double chi1, UNUSED double chi2, int RDAmpFlag) {
    UNUSED double total=0;
    switch (RDAmpFlag){
        case 122018:{
            double noSpin = 1.33;
            double eqSpin = 0.;
            double uneqSpin = 0.;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_32_sigma: version is not valid. Recommended version is 122018.");}
    }
    return total;
}

static double IMRPhenomXHM_RD_Amp_44_sigma(UNUSED double eta, UNUSED double S, UNUSED double chi1, UNUSED double chi2, int RDAmpFlag) {
    UNUSED double total=0;
    switch (RDAmpFlag){
        case 122018:{
            double noSpin = 1.33;
            double eqSpin = 0.;
            double uneqSpin = 0.;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_44_sigma: version is not valid. Recommended version is 122018.");}
    }
    return total;
}


/************** Amplitude Ringdown Ansatz *************/

// For the modes with mixing this is the ansatz of the spheroidal part.
static double IMRPhenomXHM_RD_Amp_Ansatz(double ff, IMRPhenomXHMWaveformStruct *pWF, IMRPhenomXHMAmpCoefficients *pAmp){

    int RDAmpFlag = pWF->IMRPhenomXHMRingdownAmpVersion;
    double frd = pWF->fRING;
    double fda = pWF->fDAMP;
    double dfr = ff - frd;
    double dfd = fda * pAmp->sigma;
    double lc  = pAmp->lc;
    double ampRD;

    switch ( RDAmpFlag )
    {
        case 0: /* Canonical, 3 fitted coefficients + fring, fdamp, lc that are fixed. sigma is also fixed except for the 21 mode. */
        {
            ampRD = (fda *fabs(pAmp->alambda) * pAmp->sigma)*exp(- dfr * pAmp->lambda / dfd )/ (dfr*dfr + dfd*dfd)*pow(ff,-lc);
            break;
        }
        default:
        {
            XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Ringdown_Amp_lm_Ansatz: IMRPhenomXHMRingdownAmpFitsVersion is not valid. Use version 0. \n");
        }
    }

    return ampRD;
}


/*** Derivative of the RD Ansatz for modes without mixing ***/

static double IMRPhenomXHM_RD_Amp_DAnsatz(double ff, IMRPhenomXHMWaveformStruct *pWF, IMRPhenomXHMAmpCoefficients *pAmp){

    int RDAmpFlag = pWF->IMRPhenomXHMRingdownAmpVersion;
    double frd = pWF->fRING;
    double fda = pWF->fDAMP;
    double dfd = fda * pAmp->sigma;
    double lc  = pAmp->lc;
    double lambda = pAmp->lambda;
    double DampRD;
    double numerator,denom,expon;

    switch ( RDAmpFlag )
    {
        case 0:  /* Canonical, 3 fitted coefficients + fring, fdamp, lc that are fixed. sigma is also fixed except for the 21 mode. */
        {
            numerator = fabs(pAmp->alambda)*pow(ff,-1.-lc)*( dfd*lc*(frd*frd + dfd*dfd)
                                                           + ff*( frd*frd*lambda - 2.*dfd*frd*(1.+lc) + dfd*dfd*lambda)
                                                           + ff*ff*( -2.*lambda*frd + dfd*(2.+lc) )
                                                           + ff*ff*ff*( lambda )
                                                           );
            denom = frd*frd + dfd*dfd - ff*2.*frd + ff*ff;
            expon = exp(-(ff-frd)*lambda/dfd);
            DampRD = -numerator*expon/(denom*denom);
            break;
        }
        default:
        {
            XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Ringdown_Amp_lm_Ansatz: IMRPhenomXHMRingdownAmpFitsVersion is not valid. Use version 0. \n");
        }
    }

    return DampRD;
}

/*** Derivative of the RD Ansatz for modes with mixing ***/
// It can not be obtained analytically, so we use finite difference of 4th order
static double IMRPhenomXHM_RD_Amp_NDAnsatz(double ff, IMRPhenomXHMAmpCoefficients *pAmp,  IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXAmpCoefficients *pAmp22,  IMRPhenomXPhaseCoefficients *pPhase22, IMRPhenomXWaveformStruct *pWF22){

    double df = 10e-10;
    double Nder;
    double fun2R = ff + 2*df;
    double funR  = ff + df;
    double funL  = ff - df;
    double fun2L = ff - 2*df;
    IMRPhenomX_UsefulPowers powers_of_fun;

    IMRPhenomX_Initialize_Powers(&powers_of_fun,fun2R);
    fun2R = cabs(SpheroidalToSpherical(fun2R, &powers_of_fun, pAmp22, pPhase22, pAmp, pPhase, pWFHM, pWF22));

    IMRPhenomX_Initialize_Powers(&powers_of_fun,funR);
    funR  = cabs(SpheroidalToSpherical(funR, &powers_of_fun, pAmp22, pPhase22, pAmp, pPhase, pWFHM, pWF22));

    IMRPhenomX_Initialize_Powers(&powers_of_fun,funL);
    funL  = cabs(SpheroidalToSpherical(funL, &powers_of_fun, pAmp22, pPhase22, pAmp, pPhase, pWFHM, pWF22));

    IMRPhenomX_Initialize_Powers(&powers_of_fun,fun2L);
    fun2L = cabs(SpheroidalToSpherical(fun2L, &powers_of_fun, pAmp22, pPhase22, pAmp, pPhase, pWFHM, pWF22));

    Nder = (-fun2R + 8*funR - 8*funL + fun2L )/(12*df);

    return Nder;

}

/* There are cases for some modes where the ringdown amplitude is very low compared to the inspiral and the intermediate reconstruction fails to connect them.
   For those cases we remove the two collocation points and reconstruct with a third order polynomial or with a linear one for the 32 mode. */
void IMRPhenomXHM_Ringdown_Amplitude_Veto(double *V2, double *V3, double V4, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22){

    double threshold;
    if(pWFHM->modeTag == 32){
      threshold = 0.1/(pWF22->ampNorm);
      if(1./V4 < threshold){
          *V2 = 1.;
          *V3 = 1.;
          pWFHM->IMRPhenomXHMIntermediateAmpVersion = 101;
      }
    }
    else{
      threshold = 0.01/(pWF22->ampNorm);
      if(1./V4 < threshold){
          *V2 = 1.;
          *V3 = 1.;
          pWFHM->IMRPhenomXHMIntermediateAmpVersion = 1032;
      }
    }
}



/****************************************/
/*                                      */
/*              PHASE                   */
/*                                      */
/****************************************/

/* Fits over parameter space for the ringdown phase quantities. */

// The spin parameter S = (m1^2*chi1 + m2^2*chi2)/(m1^2 + m2^2)

static double IMRPhenomXHM_RD_Phase_22_alpha2(double eta, double S, double chi1, double chi2, int RDPhaseFlag) {
    double total=0,eta2,eta3,eta4,delta=sqrt(1-4*eta);
    switch (RDPhaseFlag){
        case 122019:{
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            double noSpin = 0.2088669311744758 - 0.37138987533788487*eta + 6.510807976353186*eta2 - 31.330215053905395*eta3 + 55.45508989446867*eta4;
            double eqSpin = ((0.2393965714370633 + 1.6966740823756759*eta - 16.874355161681766*eta2 + 38.61300158832203*eta3)*S)/(1. - 0.633218538432246*S);
            double uneqSpin = (chi1 - 1.*chi2)*(0.9088578269496244*pow(eta,2.5) + 15.619592332008951*(chi1 - 1.*chi2)*pow(eta,3.5))*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Phase_22_alpha2: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_RD_Phase_22_alphaL(double eta, double S, double chi1, double chi2, int RDPhaseFlag) {
    UNUSED double total=0, delta=sqrt(1.- 4.*eta),eta2,eta3,eta4,S2;
    switch (RDPhaseFlag){
        case 122019:{
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            S2 = pow(S,2);
            double noSpin = eta*(-1.1926122248825484 + 2.5400257699690143*eta - 16.504334734464244*eta2 + 27.623649807617376*eta3);
            double eqSpin = eta3*S*(35.803988443700824 + 9.700178927988006*S - 77.2346297158916*S2) + eta*S*(0.1034526554654983 - 0.21477847929548569*S - 0.06417449517826644*S2) + eta2*S*(-4.7282481007397825 + 0.8743576195364632*S + 8.170616575493503*S2) + eta4*S*(-72.50310678862684 - 39.83460092417137*S + 180.8345521274853*S2);
            double uneqSpin = (-0.7428134042821221*chi1*pow(eta,3.5) + 0.7428134042821221*chi2*pow(eta,3.5) + 17.588573345324154*pow(chi1,2)*pow(eta,4.5) - 35.17714669064831*chi1*chi2*pow(eta,4.5) + 17.588573345324154*pow(chi2,2)*pow(eta,4.5))*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Phase_22_alphaL: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

/**************** 32 specific fits ***************/
static double IMRPhenomXHM_RD_Phase_32_SpheroidalTimeShift(double eta, double S, double chi1, double chi2, UNUSED int RDPhaseFlag) {

    double total,eta2,eta3,eta4,eta5,S2,S3,S4;
    eta2 = pow(eta,2);
    eta3 = pow(eta,3);
    eta4 = pow(eta,4);
    eta5 = pow(eta,5);
    S2 = pow(S,2);
    S3 = pow(S,3);
    S4 = pow(S,4);
    double noSpin = 11.851438981981772 + 167.95086712701223*eta - 4565.033758777737*eta2 + 61559.132976189896*eta3 - 364129.24735853914*eta4 + 739270.8814129328*eta5;
    double eqSpin = (9.506768471271634 + 434.31707030999445*eta - 8046.364492927503*eta2 + 26929.677144312944*eta3)*S + (-5.949655484033632 - 307.67253970367034*eta + 1334.1062451631644*eta2 + 3575.347142399199*eta3)*S2 + (3.4881615575084797 - 2244.4613237912527*eta + 24145.932943269272*eta2 - 60929.87465551446*eta3)*S3 + (15.585154698977842 - 2292.778112523392*eta + 24793.809334683185*eta2 - 65993.84497923202*eta3)*S4;
    double uneqSpin = 465.7904934097202*(chi1 - 1.*chi2)*sqrt(1. - 4.*eta)*eta2;
    total = noSpin + eqSpin + uneqSpin;

    return total;
}

static double IMRPhenomXHM_RD_Phase_32_SpheroidalPhaseShift(double eta, double S, double chi1, double chi2, UNUSED int RDPhaseFlag) {
    double total,eta2,eta3,eta4,eta5,eta6,eta7,S2,S3,S4;
    eta2 = pow(eta,2);
    eta3 = pow(eta,3);
    eta4 = pow(eta,4);
    eta5 = pow(eta,5);
    eta6 = pow(eta,6);
    eta7 = pow(eta,7);
    S2 = pow(S,2);
    S3 = pow(S,3);
    S4 = pow(S,4);
    double noSpin = -1.3328895897490733 - 22.209549522908667*eta + 1056.2426481245027*eta2 - 21256.376324666326*eta3 + 246313.12887984765*eta4 - 1.6312968467540336e6*eta5 + 5.614617173188322e6*eta6 - 7.612233821752137e6*eta7;
    double eqSpin = (S*(-1.622727240110213 + 0.9960210841611344*S - 1.1239505323267036*S2 - 1.9586085340429995*S3 + eta2*(196.7055281997748 + 135.25216875394943*S + 1086.7504825459278*S2 + 546.6246807461155*S3 - 312.1010566468068*S4) + 0.7638287749489343*S4 + eta*(-47.475568056234245 - 35.074072557604445*S - 97.16014978329918*S2 - 34.498125910065156*S3 + 24.02858084544326*S4) + eta3*(62.632493533037625 - 22.59781899512552*S - 2683.947280170815*S2 - 1493.177074873678*S3 + 805.0266029288334*S4)))/(-2.950271397057221 + 1.*S);
    double uneqSpin = (sqrt(1. - 4.*eta)*(chi2*pow(eta,2.5)*(88.56162028006072 - 30.01812659282717*S) + chi2*eta2*(43.126266433486435 - 14.617728550838805*S) + chi1*eta2*(-43.126266433486435 + 14.617728550838805*S) + chi1*pow(eta,2.5)*(-88.56162028006072 + 30.01812659282717*S)))/(-2.950271397057221 + 1.*S);
    total = noSpin + eqSpin + uneqSpin;

    return total;
}

// collocation points
static double IMRPhenomXHM_Ringdown_Phase_32_p1(double eta, double S, double chi1, double chi2, UNUSED int RingdownPhaseFlag) {
    double total,eta2,eta3,eta4,eta5,S2,S3,S4,S5;
    eta2 = pow(eta,2);
    eta3 = pow(eta,3);
    eta4 = pow(eta,4);
    eta5 = pow(eta,5);
    S2 = pow(S,2);
    S3 = pow(S,3);
    S4 = pow(S,4);
    S5 = pow(S,5);
    double noSpin = 3169.372056189274 + 426.8372805022653*eta - 12569.748101922158*eta2 + 149846.7281073725*eta3 - 817182.2896823225*eta4 + 1.5674053633767858e6*eta5;
    double eqSpin = (19.23408352151287 - 1762.6573670619173*eta + 7855.316419853637*eta2 - 3785.49764771212*eta3)*S + (-42.88446003698396 + 336.8340966473415*eta - 5615.908682338113*eta2 + 20497.5021807654*eta3)*S2 + (13.918237996338371 + 10145.53174542332*eta - 91664.12621864353*eta2 + 201204.5096556517*eta3)*S3 + (-24.72321125342808 - 4901.068176970293*eta + 53893.9479532688*eta2 - 139322.02687945773*eta3)*S4 + (-61.01931672442576 - 16556.65370439302*eta + 162941.8009556697*eta2 - 384336.57477596396*eta3)*S5;
    double uneqSpin = (chi1 - 1.*chi2)*sqrt(1. - 4.*eta)*eta2*(641.2473192044652 - 1600.240100295189*chi1*eta + 1600.240100295189*chi2*eta + 13275.623692212472*eta*S);
    total = noSpin + eqSpin + uneqSpin;

    return total;
}
// collocation points
static double IMRPhenomXHM_Ringdown_Phase_32_p2(double eta, double S, double chi1, double chi2, UNUSED int RingdownPhaseFlag) {
    double total,eta2,eta3,S2,S3,S4;
    eta2 = pow(eta,2);
    eta3 = pow(eta,3);
    S2 = pow(S,2);
    S3 = pow(S,3);
    S4 = pow(S,4);
    double noSpin = 3131.0260952676376 + 206.09687819102305*eta - 2636.4344627081873*eta2 + 7475.062269742079*eta3;
    double eqSpin = (49.90874152040307 - 691.9815135740145*eta - 434.60154548208334*eta2 + 10514.68111669422*eta3)*S + (97.3078084654917 - 3458.2579971189534*eta + 26748.805404989867*eta2 - 56142.13736008524*eta3)*S2 + (-132.49105074500454 + 429.0787542102207*eta + 7269.262546204149*eta2 - 27654.067482558712*eta3)*S3 + (-227.8023564332453 + 5119.138772157134*eta - 34444.2579678986*eta2 + 69666.01833764123*eta3)*S4;
    double uneqSpin = 477.51566939885424*(chi1 - 1.*chi2)*sqrt(1. - 4.*eta)*eta2;
    total = noSpin + eqSpin + uneqSpin;

    return total;
}
// collocation points
static double IMRPhenomXHM_Ringdown_Phase_32_p3(double eta, double S, double chi1, double chi2, UNUSED int RingdownPhaseFlag) {
    double total,eta2,eta3,S2,S3,S4;
    eta2 = pow(eta,2);
    eta3 = pow(eta,3);
    S2 = pow(S,2);
    S3 = pow(S,3);
    S4 = pow(S,4);
    double noSpin = 3082.803556599222 + 76.94679795837645*eta - 586.2469821978381*eta2 + 977.6115755788503*eta3;
    double eqSpin = (45.08944710349874 - 807.7353772747749*eta + 1775.4343704616288*eta2 + 2472.6476419567534*eta3)*S + (95.57355060136699 - 2224.9613131172046*eta + 13821.251641893134*eta2 - 25583.314298758105*eta3)*S2 + (-144.96370424517866 + 2268.4693587493093*eta - 10971.864789147161*eta2 + 16259.911572457446*eta3)*S3 + (-227.8023564332453 + 5119.138772157134*eta - 34444.2579678986*eta2 + 69666.01833764123*eta3)*S4;
    double uneqSpin = 378.2359918274837*(chi1 - 1.*chi2)*sqrt(1. - 4.*eta)*eta2;
    total = noSpin + eqSpin + uneqSpin;

    return total;
}
// collocation points
static double IMRPhenomXHM_Ringdown_Phase_32_p4(double eta, double S, double chi1, double chi2, UNUSED int RingdownPhaseFlag) {
    double total,eta2,S2,S3,S4;
    eta2 = pow(eta,2);
    S2 = pow(S,2);
    S3 = pow(S,3);
    S4 = pow(S,4);
    double noSpin = 3077.0657367004565 + 64.99844502520415*eta - 357.38692756785395*eta2;
    double eqSpin = (34.793450080444714 - 986.7751755509875*eta - 9490.641676924794*pow(eta,3) + 5700.682624203565*eta2)*S + (57.38106384558743 - 1644.6690499868596*eta - 19906.416384606226*pow(eta,3) + 11008.881935880598*eta2)*S2 + (-126.02362949830213 + 3169.3397351803583*eta + 62863.79877094988*pow(eta,3) - 26766.730897942085*eta2)*S3 + (-169.30909412804587 + 4900.706039920717*eta + 95314.99988114933*pow(eta,3) - 41414.05689348732*eta2)*S4;
    double uneqSpin = 390.5443469721231*(chi1 - 1.*chi2)*sqrt(1. - 4.*eta)*eta2;
    total = noSpin + eqSpin + uneqSpin;

    return total;
}


/**************  ANSATZ PHASE DERIVATIVE **************/

static double IMRPhenomXHM_RD_Phase_Ansatz(double ff, IMRPhenomX_UsefulPowers *powers_of_f,IMRPhenomXHMWaveformStruct *pWFHM,  IMRPhenomXHMPhaseCoefficients *pPhase){

    double invf2 = powers_of_f->m_two;
    double frd   = pWFHM->fRING;
    double fda   = pWFHM->fDAMP;
    double dphaseRD;

    switch ( pWFHM->MixingOn )
    {
        case 0:
        {
            /*  rescaling of the 22 ansatz -- used for (21),(33),(44) */
            /*ansatz:
             alpha0 + ((fRDlm^2) alpha2)/(f^2)  + alphaL*(fdamplm)/((fdamplm)^2 + (f - fRDlm)^2)*/
            dphaseRD = ( pPhase->alpha0 +  frd*frd*(pPhase->alpha2)*invf2 + ( (pPhase->alphaL)* fda/(fda*fda +(ff - frd)*(ff - frd)) ) );
            break;

        }
        case 1:
        {
            /*  calibration of spheroidal ringdown waveform for (32) */
            /* ansatz: alpha0 + (alpha2)/(f^2)+ (alpha4)/(f^4)  + alphaL*(fdamplm)/((fdamplm)^2 + (f - fRDlm)^2)*/
            double invf4 = powers_of_f->m_four;
            dphaseRD = ( pPhase->alpha0_S +  (pPhase->alpha2_S)*invf2 + (pPhase->alpha4_S)*invf4 +( (pPhase->alphaL_S)* fda/(fda*fda +(ff - frd)*(ff - frd)) ) );
            break;
        }
        default:
        {XLAL_ERROR(XLAL_EDOM, "Error in IMRPhenomXHM_RD_Phase_Ansatz: version is not valid. Use version 0 for modes (2,1),(3,3),(4,4) and 1 for (3,2).\n");}
    }
    return dphaseRD;
}

/**************  ANSATZ INTEGRATED PHASE **************/

static double IMRPhenomXHM_RD_Phase_AnsatzInt(double ff, IMRPhenomX_UsefulPowers *powers_of_f,IMRPhenomXHMWaveformStruct *pWFHM,  IMRPhenomXHMPhaseCoefficients *pPhase){

    double invf   = powers_of_f->m_one;
    double phaseRD;
    double frd=pWFHM->fRING;
    double fda=pWFHM->fDAMP;

    switch ( pWFHM->MixingOn )
    {
        case 0:
        {
            /*  rescaling of the 22 ansatz -- used for (21),(33),(44) */
            /*ansatz:
             alpha0 f - fRDlm^2*alpha2)/f  + alphaL*ArcTan[(f - fRDlm)/fdamplm]*/
            phaseRD = pPhase->alpha0*ff -frd*frd *(pPhase->alpha2)*invf +  (pPhase->alphaL)* atan((ff-frd)/fda);
            break;
        }
        case 1:
        {
            /*  calibration of spheroidal ringdown waveform for (32) */
            /* ansatz: f alpha0 - (alpha4)/(3 f^3) - (alpha2)/f + alphaL ArcTan[(f - fRDlm)/fdamplm]*/
            double invf3 = powers_of_f->m_three;
            phaseRD = pPhase->phi0_S+pPhase->alpha0_S*ff -(pPhase->alpha2_S)*invf -1./3.*(pPhase->alpha4_S)*invf3 +(pPhase->alphaL_S)* atan((ff-frd)/fda);
            break;
        }
        default:
        {XLAL_ERROR(XLAL_EDOM, "Error in IMRPhenomXHM_RD_Phase_AnsatzInt: version is not valid. Use version 0 for modes (2,1),(3,3),(4,4) and 1 for (3,2).\n");}
    }
    return phaseRD;
}
