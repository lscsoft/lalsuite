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
//
//  LALSimIMRPhenomXHM_Intermediate.c
//

#include <gsl/gsl_poly.h>

#include "LALSimIMRPhenomXHM_intermediate.h"

/***********************************************/
/*                                             */
/*                  AMPLITUDE                  */
/*                                             */
/***********************************************/


// Fits of the intermediate collocation points over parameter space.

/* These parameter-space fits are documented in the supplementary material of https://dcc.ligo.org/P2000011-v2.
   There are 2 Mathematica notebooks (one for amplitude and one for phase) that read the fits data and automatically generate the C-code below.
   For more information read https://git.ligo.org/waveforms/reviews/imrphenomx/blob/master/documentation/ParspaceFits/README and the documentation in the notebooks. */

// IMRPhenomXHM_Inter_Amp_lm_intX returns the value of the amplitude at the collocation point intX
// IMRPhenomXHM_Inter_Amp_lm_dintX returns the value of the derivative of the amplitude at the collocation point intX
// for a definition of the frequencies corresponding to int0/1/2, see Sec V.A.

// The dominant spin parameter is S = (m1^2*chi1 + m2^2*chi2)/(m1^2 + m2^2)

/* Start of Amp Parameter Space Fits */

static double IMRPhenomXHM_Inter_Amp_21_int1(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag) {
  double total=0;
  switch (InterAmpFlag){
    case 122018:{
      double eta = pWF->eta;
      double S = pWF->STotR;
      double eta2,S2;
      eta2 = pow(eta,2);
      S2 = pow(S,2);
      double noSpin = sqrt(eta - 4.*eta2)*(21.256776327599113 - 25.594352690383847*eta + 30.14761650482866*eta2);
      double eqSpin = sqrt(eta - 4.*eta2)*S*(-11.262044985632757 - 1.8167045597937677*S + eta*(-1.1798437990445079 + 6.344825546437461*S - 4.881427482271166*S2));
      double uneqSpin = -3.6366100759176696*pow(pWF->dchi,2)*(1. - 4.*eta)*eta - 31.60048733143782*(pWF->dchi)*eta2*(1. + 2.1502870640831855*eta2);
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
        case 122022:{
            double eta = pWF->eta;
            double delta = pWF->delta;
            double S = pWF->STotR;
            double chidiff = pWF->dchi_half;
            double eta1 = eta;
            double eta2 = eta1 * eta1;
            double eta3 = eta1 * eta2;
            double eta4 = eta1 * eta3;
            double eta5 = eta1 * eta4;
            double S1 = S;
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = fabs(delta*eta1*(chidiff2*(5.159755997682368*eta1 - 30.293198248154948*eta2 + 63.70715919820867*eta3) + chidiff1*(8.262642080222694*eta1 - 415.88826990259116*eta2 + 1427.5951158851076*eta3)) + delta*eta1*(18.55363583212328 - 66.46950491124205*eta1 + 447.2214642597892*eta2 - 1614.178472020212*eta3 + 2199.614895727586*eta4) + chidiff1*eta5*(-1698.841763891122 - 195.27885562092342*S1 - 1.3098861736238572*S2) + delta*eta1*(chidiff1*(34.17829404207186*eta1 - 386.34587928670015*eta2 + 1022.8553774274128*eta3)*S1 + chidiff1*(56.76554600963724*eta1 - 491.4593694689354*eta2 + 1016.6019654342113*eta3)*S2) + delta*eta1*S1*(-8.276366844994188*(1.0677538075697492 - 24.12941323757896*eta1 + 516.7886322104276*eta2 - 4389.799658723288*eta3 + 16770.447637953577*eta4 - 23896.392706809565*eta5) - 1.6908277400304084*(3.4799140066657928 - 29.00026389706585*eta1 + 114.8330693231833*eta2 - 184.13091281984674*eta3 + 592.300353344717*eta4 - 2085.0821513466053*eta5)*S1 - 0.46006975902558517*(-2.1663474937625975 + 826.026625945615*eta1 - 17333.549622759732*eta2 + 142904.08962903373*eta3 - 528521.6231015554*eta4 + 731179.456702448*eta5)*S2));
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_21_int1: version %i is not valid.", InterAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_Inter_Amp_21_int2(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag) {
  double total=0;
  switch (InterAmpFlag){
    case 122018:{
      double eta = pWF->eta;
      double S = pWF->STotR;
      double eta2 = pow(eta,2);
      double noSpin = sqrt(eta - 4.*eta2)*(19.15445065708005 - 21.13596229438309*eta + 29.742565944285772*eta2);
      double eqSpin = sqrt(eta - 4.*eta2)*S*(-12.766814596085734 - 2.123816950673979*S + eta*(-2.913184982025043 + 6.006571549661901*S));
      double uneqSpin = -25.856046423804255*(pWF->dchi)*eta2*(1. + 5.7871199275552*eta2);
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
        case 122022:{
            double eta = pWF->eta;
            double delta = pWF->delta;
            double S = pWF->STotR;
            double chidiff = pWF->dchi_half;
            double eta1 = eta;
            double eta2 = eta1 * eta1;
            double eta3 = eta1 * eta2;
            double eta4 = eta1 * eta3;
            double eta5 = eta1 * eta4;
            double S1 = S;
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            total = fabs(delta*eta1*(13.757856231617446 - 12.783698329428516*eta1 + 12.048194546899204*eta2) + chidiff1*delta*eta1*(15.107530092096438*eta1 - 416.811753638553*eta2 + 1333.6181181686939*eta3) + chidiff1*eta5*(-1549.6199518612063 - 102.34716990474509*S1 - 3.3637011939285015*S2) + delta*eta1*(chidiff1*(36.358142200869295*eta1 - 384.2123173145321*eta2 + 984.6826660818275*eta3)*S1 + chidiff1*(4.159271594881928*eta1 + 105.10911749116399*eta2 - 639.190132707115*eta3)*S2) + delta*eta1*S1*(-8.097876227116853*(0.6569459700232806 + 9.861355377849485*eta1 - 116.88834714736281*eta2 + 593.8035334117192*eta3 - 1063.0692862578455*eta4) - 1.0546375154878165*(0.745557030602097 + 65.25215540635162*eta1 - 902.5751736558435*eta2 + 4350.442990924205*eta3 - 7141.611333893155*eta4)*S1 - 0.5006664599166409*(10.289020582277626 - 212.00728173197498*eta1 + 2334.0029399672358*eta2 - 11939.621138801092*eta3 + 21974.8201355744*eta4)*S2));
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_21_int2: version %i is not valid.", InterAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_Inter_Amp_33_int1(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag) {
  double total=0;
  switch (InterAmpFlag){
    case 122018:{
      double eta = pWF->eta;
      double S = pWF->STotR;
      double eta2,eta3,eta4,eta6;
      eta2 = pow(eta,2);
      eta3 = pow(eta,3);
      eta4 = pow(eta,4);
      eta6 = pow(eta,6);
      double noSpin = sqrt(eta - 4.*eta2)*(27.927652424857733 - 133.56611389260297*eta + 974.8550901501316*eta2 - 3744.785831952632*eta3 + 5621.897260910284*eta4);
      double eqSpin = sqrt(eta - 4.*eta2)*S*(7.348313807306079 + eta*(-60.248696675045565 - 37.07212326362276*S) + 5.059236579431119*S + eta2*(159.68630712802727 + 83.33807316873204*S));
      double uneqSpin = 1412.367880056888*(pWF->dchi)*eta6;
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
        case 122022:{
            double eta = pWF->eta;
            double delta = pWF->delta;
            double S = pWF->STotR;
            double chidiff = pWF->dchi_half;
            double eta1 = eta;
            double eta2 = eta1 * eta1;
            double eta3 = eta1 * eta2;
            double eta4 = eta1 * eta3;
            double eta5 = eta1 * eta4;
            double S1 = S;
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            total = chidiff1*delta*eta1*(-0.3516244197696068*eta1 + 40.425151307421416*eta2 - 148.3162618111991*eta3) + delta*eta1*(26.998512565991778 - 146.29035440932105*eta1 + 914.5350366065115*eta2 - 3047.513201789169*eta3 + 3996.417635728702*eta4) + chidiff1*delta*eta1*(5.575274516197629*eta1 - 44.592719238427094*eta2 + 99.91399033058927*eta3)*S1 + delta*eta1*S1*(-0.5383304368673182*(-7.456619067234563 + 129.36947401891433*eta1 - 843.7897535238325*eta2 + 3507.3655567272644*eta3 - 9675.194644814854*eta4 + 11959.83533107835*eta5) - 0.28042799223829407*(-6.212827413930676 + 266.69059813274475*eta1 - 4241.537539226717*eta2 + 32634.43965039936*eta3 - 119209.70783201039*eta4 + 166056.27237509796*eta5)*S1) + chidiff1*eta5*(199.6863414922219 + 53.36849263931051*S1 + 7.650565415855383*S2);
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_33_int1: version %i is not valid.", InterAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_Inter_Amp_33_int2(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag) {
  double total=0;
  switch (InterAmpFlag){
    case 122018:{
      double eta = pWF->eta;
      double S = pWF->STotR;
      double eta2,eta6;
      eta2 = pow(eta,2);
      eta6 = pow(eta,6);
      double noSpin = sqrt(eta - 4.*eta2)*(20.162169689041903 - 18.666422946967764*eta + 53.04107631052987*eta2);
      double eqSpin = sqrt(eta - 4.*eta2)*S*(3.896260108714186 + eta*(-33.707998325000965 - 61.1244771077077*S) + 4.878506403725656*S + eta2*(91.31681057861915 + 196.40535070402336*S));
      double uneqSpin = 1637.4256048973248*(pWF->dchi)*eta6;
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
        case 122022:{
            double eta = pWF->eta;
            double delta = pWF->delta;
            double S = pWF->STotR;
            double chidiff = pWF->dchi_half;
            double eta1 = eta;
            double eta2 = eta1 * eta1;
            double eta3 = eta1 * eta2;
            double eta4 = eta1 * eta3;
            double eta5 = eta1 * eta4;
            double S1 = S;
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            total = delta*eta1*(17.42562079069636 - 28.970875603981295*eta1 + 50.726220750178435*eta2) + chidiff1*delta*eta1*(-7.861956897615623*eta1 + 93.45476935080045*eta2 - 273.1170921735085*eta3) + chidiff1*delta*eta1*(-0.3265505633310564*eta1 - 9.861644053348053*eta2 + 60.38649425562178*eta3)*S1 + chidiff1*eta5*(234.13476431269862 + 51.2153901931183*S1 - 10.05114600643587*S2) + delta*eta1*S1*(0.3104472390387834*(6.073591341439855 + 169.85423386969634*eta1 - 4964.199967099143*eta2 + 42566.59565666228*eta3 - 154255.3408672655*eta4 + 205525.13910847943*eta5) + 0.2295327944679772*(19.236275867648594 - 354.7914372697625*eta1 + 1876.408148917458*eta2 + 2404.4151687877525*eta3 - 41567.07396803811*eta4 + 79210.33893514868*eta5)*S1 + 0.30983324991828787*(11.302200127272357 - 719.9854052004307*eta1 + 13278.047199998868*eta2 - 104863.50453518033*eta3 + 376409.2335857397*eta4 - 504089.07690692553*eta5)*S2);
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_33_int2: version %i is not valid.", InterAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_Inter_Amp_32_int1(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag) {
  double total=0;
  switch (InterAmpFlag){
    case 122018:{
      double eta = pWF->eta;
      double S = pWF->STotR;
      double delta=sqrt(1.-4*eta),eta2,eta3,eta4,eta5,eta6;
      eta2 = pow(eta,2);
      eta3 = pow(eta,3);
      eta4 = pow(eta,4);
      eta5 = pow(eta,5);
      eta6 = pow(eta,6);
      double noSpin = sqrt(eta - 3.*eta2)*(6.523612598187996 - 56.93956111746338*eta + 1021.6414686597869*eta2 - 12107.114370361525*eta3 + 76320.90587515048*eta4 - 244144.92645448362*eta5 + 321790.55131499085*eta6);
      double eqSpin = sqrt(eta - 3.*eta2)*S*(2.9649243713119895 + eta3*(1790.8363334078751 - 5438.911035114849*S) + eta*(-37.87005271181108 - 126.1263286618178*S) + 4.063724538613828*S + eta2*(48.39743086535961 + 1341.2619677741804*S) + eta4*(-5200.659417644607 + 7369.386205324284*S));
      double uneqSpin = eta2*(-0.4386152975075188*(pow(pWF->chi1L,2) - 2.*pWF->chi1L*pWF->chi2L + pow(pWF->chi2L,2)) + (pWF->chi2L*(3.6527252109313233 - 7.324266404418883*S) + pWF->chi1L*(-3.6527252109313233 + 7.324266404418883*S))*delta);
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
        case 122022:{
            double eta = pWF->eta;
            double delta = pWF->delta;
            double sqroot = sqrt(eta);
            double S = pWF->chiPNHat;
            double chidiff = pWF->dchi_half;
            double eta1 = eta;
            double eta2 = eta1 * eta1;
            double eta3 = eta1 * eta2;
            double eta4 = eta1 * eta3;
            double eta5 = eta1 * eta4;
            double S1 = S;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = (chidiff2*(-0.2341404256829785*eta1 + 2.606326837996192*eta2 - 8.68296921440857*eta3) + chidiff1*delta*(0.5454562486736877*eta1 - 25.19759222940851*eta2 + 73.40268975811729*eta3))*sqroot + chidiff1*delta*(0.4422257616009941*eta1 - 8.490112284851655*eta2 + 32.22238925527844*eta3)*S1*sqroot + S1*(0.7067243321652764*(0.12885110296881636 + 9.608999847549535*eta1 - 85.46581740280585*eta2 + 325.71940024255775*eta3 + 175.4194342269804*eta4 - 1929.9084724384807*eta5) + 0.1540566313813899*(-0.3261041495083288 + 45.55785402900492*eta1 - 827.591235943271*eta2 + 7184.647314370326*eta3 - 28804.241518798244*eta4 + 43309.69769878964*eta5)*S1)*sqroot + (480.0434256230109*eta1 + 25346.341240810478*eta2 - 99873.4707358776*eta3 + 106683.98302194536*eta4)*sqroot*pow(1 + 1082.6574834474493*eta1 + 10083.297670051445*eta2,-1);
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_32_int1: version %i is not valid.", InterAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_Inter_Amp_32_int2(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag) {
  double total=0;
  switch (InterAmpFlag){
    case 122018:{
      double eta = pWF->eta;
      double S = pWF->STotR;
      double delta=sqrt(1.-4.*eta),eta2,eta3,eta4,S2;
      eta2 = pow(eta,2);
      eta3 = pow(eta,3);
      eta4 = pow(eta,4);
      S2 = pow(S,2);
      double noSpin = sqrt(eta - 3.*eta2)*(5.941845842405418 - 31.905244419036794*eta + 271.105632998832*eta2 - 2113.9652334868965*eta3 + 6214.038393898584*eta4);
      double eqSpin = sqrt(eta - 3.*eta2)*S*(-2.726472456645038 + 2.9454485454761827*S + eta3*(10581.664858726683 - 8474.190197512324*S - 11680.937129551317*S2) + eta*(98.08119212251981 - 119.88112323140916*S - 145.5079981415436*S2) + 3.5684571473795095*S2 + eta2*(-1595.8027347570667 + 1686.7137359336039*S + 2139.8290160628144*S2) + eta4*(-21488.25117198268 + 13866.428366595079*S + 20863.270079587106*S2));
      double uneqSpin = 0.0038732029045487884*(pWF->dchi)*eta2*delta;
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
        case 122022:{
            double eta = pWF->eta;
            double delta = pWF->delta;
            double S = pWF->STotR;
            double chidiff = pWF->dchi_half;
            double eta1 = eta;
            double eta2 = eta1 * eta1;
            double eta3 = eta1 * eta2;
            double eta4 = eta1 * eta3;
            double eta5 = eta1 * eta4;
            double eta6 = eta1 * eta5;
            double S1 = S;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = eta1*(chidiff2*(-4.175680729484314*eta1 + 47.54281549129226*eta2 - 128.88334273588077*eta3) + chidiff1*delta*(-0.18274358639599947*eta1 - 71.01128541687838*eta2 + 208.07105580635888*eta3)) + eta1*(4.760999387359598 - 38.57900689641654*eta1 + 456.2188780552874*eta2 - 4544.076411013166*eta3 + 24956.9592553473*eta4 - 69430.10468748478*eta5 + 77839.74180254337*eta6) + chidiff1*delta*eta1*(1.2198776533959694*eta1 - 26.816651899746475*eta2 + 68.72798751937934*eta3)*S1 + eta1*S1*(1.5098291294292217*(0.4844667556328104 + 9.848766999273414*eta1 - 143.66427232396376*eta2 + 856.9917885742416*eta3 - 1633.3295758142904*eta4) + 0.32413108737204144*(2.835358206961064 - 62.37317183581803*eta1 + 761.6103793011912*eta2 - 3811.5047139343505*eta3 + 6660.304740652403*eta4)*S1);
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_32_int2: version %i is not valid.", InterAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_Inter_Amp_44_int1(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag) {
  double total=0;
  switch (InterAmpFlag){
    case 122018:{
      double eta = pWF->eta;
      double S = pWF->STotR;
      double delta=sqrt(1.-4.*eta),eta2,eta3,eta4;
      eta2 = pow(eta,2);
      eta3 = pow(eta,3);
      eta4 = pow(eta,4);
      double noSpin = sqrt(eta - 3.*eta2)*(10.804555518381166 - 72.3834734399584*eta + 540.0541240482852*eta2 - 2612.999845214264*eta3 + 4779.096001663427*eta4);
      double eqSpin = sqrt(eta - 3.*eta2)*S*(4.26336253142121 + eta*(-47.94914754514519 - 39.31284390368824*S) + 3.0973959822174297*S + eta2*(119.70401520575753 + 106.91295627237112*S));
      double uneqSpin = 0.7262636326998003*pow(pWF->dchi,2)*(1. - 4.*eta)*eta + 3.001401833124412*(pWF->dchi)*eta2*delta;
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
        case 122022:{
            double eta = pWF->eta;
            double delta = pWF->delta;
            double S = pWF->chiPNHat;
            double chidiff = pWF->dchi_half;
            double eta1 = eta;
            double eta2 = eta1 * eta1;
            double eta3 = eta1 * eta2;
            double eta4 = eta1 * eta3;
            double S1 = S;
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = eta1*(chidiff1*delta*(1.5378890240544967*eta1 - 3.4499418893734903*eta2 + 16.879953490422782*eta3) + chidiff2*(1.720226708214248*eta1 - 11.87925165364241*eta2 + 23.259283336239545*eta3)) + eta1*(8.790173464969538 - 64.95499142822892*eta1 + 324.1998823562892*eta2 - 1111.9864921907126*eta3 + 1575.602443847111*eta4) + eta1*S1*(-0.062333275821238224*(-21.630297087123807 + 137.4395894877131*eta1 + 64.92115530780129*eta2 - 1013.1110639471394*eta3) - 0.11014697070998722*(4.149721483857751 - 108.6912882442823*eta1 + 831.6073263887092*eta2 - 1828.2527520190122*eta3)*S1 - 0.07704777584463054*(4.581767671445529 - 50.35070009227704*eta1 + 344.9177692251726*eta2 - 858.9168637051405*eta3)*S2);
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_44_int1: version %i is not valid.", InterAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_Inter_Amp_44_int2(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag) {
  double total=0;
  switch (InterAmpFlag){
    case 122018:{
      double eta = pWF->eta;
      double S = pWF->STotR;
      double delta=sqrt(1.-4.*eta),eta2,eta3,eta4;
      eta2 = pow(eta,2);
      eta3 = pow(eta,3);
      eta4 = pow(eta,4);
      double noSpin = sqrt(eta - 3.*eta2)*(9.020721305469884 - 53.221883492311235*eta + 508.07176447172264*eta2 - 3194.0620894511508*eta3 + 6769.9274392345915*eta4);
      double eqSpin = sqrt(eta - 3.*eta2)*S*(3.256591670091969 + eta*(-38.38922554651356 - 25.286684856422735*S) + 2.374434219852751*S + eta2*(96.41777041220982 + 64.74544118094362*S));
      double uneqSpin = 3.2337593375595417*(pWF->dchi)*eta2*delta;
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
            case 122022:{
            double eta = pWF->eta;
            double delta = pWF->delta;
            double S = pWF->chiPNHat;
            double chidiff = pWF->dchi_half;
            double eta1 = eta;
            double eta2 = eta1 * eta1;
            double eta3 = eta1 * eta2;
            double eta4 = eta1 * eta3;
            double eta5 = eta1 * eta4;
            double eta6 = eta1 * eta5;
            double S1 = S;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = eta1*(chidiff1*delta*(2.3123974306694057*eta1 - 12.237594841284904*eta2 + 44.78225529547671*eta3) + chidiff2*(2.9282931698944292*eta1 - 25.624210264341933*eta2 + 61.05270871360041*eta3)) + eta1*(6.98072197826729 - 46.81443520117986*eta1 + 236.76146303619544*eta2 - 920.358408667518*eta3 + 1478.050456337336*eta4) + eta1*S1*(-0.07801583359561987*(-28.29972282146242 + 752.1603553640072*eta1 - 10671.072606753183*eta2 + 83447.0461509547*eta3 - 350025.2112501252*eta4 + 760889.6919776166*eta5 - 702172.2934567826*eta6) + 0.013159545629626014*(91.1469833190294 - 3557.5003799977294*eta1 + 52391.684517955284*eta2 - 344254.9973814295*eta3 + 1.0141877915334814e6*eta4 - 1.1505186449682908e6*eta5 + 268756.85659532435*eta6)*S1);
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_44_int2: version %i is not valid.", InterAmpFlag);}
  }
  return total;
}

/*
Fits for the extra collocation point for EMR cases with 2 intermediate regions
*/
static double IMRPhenomXHM_Inter_Amp_21_int0(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag) {
  double total=0;
  switch (InterAmpFlag){
    case 122018:{
      double eta = pWF->eta;
      double S = pWF->STotR;
      double eta2,eta3;
      eta2 = pow(eta,2);
      eta3 = pow(eta,3);
      double noSpin = 0.872895771366973 + 441.76285124642845*eta - 24617.068739152524*eta2 + 518054.9485981792*eta3;
      double eqSpin = S*(-0.0720494539485585 + eta*(-173.67847091983123 - 113.29725582509889*S) - 0.2687302438646897*S + eta2*(3571.0393588230045 + 2640.919925429635*S));
      double uneqSpin = 0.;
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_21_int0: version is not valid.");}
  }
  return total;
}

static double IMRPhenomXHM_Inter_Amp_21_dint0(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag) {
  double total=0;
  switch (InterAmpFlag){
    case 122018:{
      double eta = pWF->eta;
      double S = pWF->STotR;
      double eta2,eta3;
      eta2 = pow(eta,2);
      eta3 = pow(eta,3);
      double noSpin = -0.8535048463050732 - 93.1876950411214*eta + 13641.071903017495*eta2 - 337621.44851304166*eta3;
      double eqSpin = S*(-1.2067842398131878 + eta2*(-1972.284151572111 - 8172.057025783849*S) - 0.26539816223182355*S + eta*(77.26350785961219 + 189.63365484152857*S));
      double uneqSpin = 0.;
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_21_dint0: version %i is not valid.", InterAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_Inter_Amp_33_int0(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag) {
  double total=0;
  switch (InterAmpFlag){
    case 122018:{
      double eta = pWF->eta;
      double S = pWF->STotR;
      double eta2,eta3;
      eta2 = pow(eta,2);
      eta3 = pow(eta,3);
      double noSpin = 1.5852399637975103 + 549.5183711492834*eta - 34257.76380246282*eta2 + 743142.8286902909*eta3;
      double eqSpin = S*(0.7436306553052219 + eta*(-89.49451655594787 - 174.5730646548662*S) + 0.4253024979725725*S + eta2*(1185.1654325913717 + 6510.983041407191*S));
      double uneqSpin = 0.;
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_33_int0: version %i is not valid.", InterAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_Inter_Amp_33_dint0(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag) {
  double total=0;
  switch (InterAmpFlag){
    case 122018:{
      double eta = pWF->eta;
      double S = pWF->STotR;
      double eta2,eta3;
      eta2 = pow(eta,2);
      eta3 = pow(eta,3);
      double noSpin = -4.691600252198376 + 101.4338937535679*eta + 9262.994550540048*eta2 - 310993.1309846956*eta3;
      double eqSpin = S*(-4.198232394219111 + eta2*(-28714.904192060643 - 5100.09336069277*S) - 0.40986595512314733*S + eta*(734.7118618746317 + 292.04566260701574*S));
      double uneqSpin = 0.;
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_33_dint0: version %i is not valid.", InterAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_Inter_Amp_32_int0(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag) {
  double total=0;
  switch (InterAmpFlag){
    case 122018:{
      double eta = pWF->eta;
      double S = pWF->STotR;
      double eta2,eta3,S2,S3;
      eta2 = pow(eta,2);
      eta3 = pow(eta,3);
      S2 = pow(S,2);
      S3 = pow(S,3);
      double noSpin = 0.24794156582503746 + 115.81823862983131*eta - 6626.167995915723*eta2 + 141004.29332593994*eta3;
      double eqSpin = (0.21144389781375486 + 35.10041265469983*eta - 1794.2301585086836*eta2)*S + (0.2781735549493081 - 37.038950686633*eta + 1258.628375238807*eta2)*S2 + (0.23428222791962147 - 63.98011009365723*eta + 2118.213562899934*eta2)*S3;
      double uneqSpin = 0.;
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_32_int0: version %i is not valid.", InterAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_Inter_Amp_32_dint0(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag) {
  double total=0;
  switch (InterAmpFlag){
    case 122018:{
      double eta = pWF->eta;
      double S = pWF->STotR;
      double eta2,eta3;
      eta2 = pow(eta,2);
      eta3 = pow(eta,3);
      double noSpin = -0.3391808620221253 - 14.604141885467747*eta + 3694.1706648870427*eta2 - 95482.02951271653*eta3;
      double eqSpin = S*(-1.2844502090793946 + eta2*(-5018.762853306415 - 6332.389157828062*S) - 1.2356159239385598*S + eta*(149.04865679660233 + 188.2052849646003*S));
      double uneqSpin = 0.;
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_32_dint0: version %i is not valid.", InterAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_Inter_Amp_44_int0(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag) {
  double total=0;
  switch (InterAmpFlag){
    case 122018:{
      double eta = pWF->eta;
      double S = pWF->STotR;
      double eta2,eta3,S2,S3;
      eta2 = pow(eta,2);
      eta3 = pow(eta,3);
      S2 = pow(S,2);
      S3 = pow(S,3);
      double noSpin = 0.5664660641971224 + 185.58965113823874*eta - 11458.768824989507*eta2 + 249386.7511724409*eta3;
      double eqSpin = (0.1741768776210781 - 9.365114803167128*eta + 703.2622732011035*eta2)*S + (0.20169229783048184 - 62.13147149352512*eta + 2833.5738711424974*eta2)*S2 + (0.4423803798742513 - 23.60535149579996*eta - 994.9241585715828*eta2)*S3;
      double uneqSpin = 0.;
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_44_int0: version %i is not valid.", InterAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_Inter_Amp_44_dint0(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag) {
  double total=0;
  switch (InterAmpFlag){
    case 122018:{
      double eta = pWF->eta;
      double S = pWF->STotR;
      double eta2 = pow(eta,2);
      double noSpin = -1.796444922382065 + 111.51170611049032*eta - 1728.7493675776548*eta2;
      double eqSpin = S*(-1.842119860613924 + eta2*(-11235.484645624338 - 2927.019210835522*S) - 0.36655273031432567*S + eta*(312.34531117524097 + 128.64488103364167*S));
      double uneqSpin = 0.;
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_44_dint0: version %i is not valid.", InterAmpFlag);}
  }
  return total;
}


static double IMRPhenomXHM_Inter_Amp_21_int3(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag){
	double total=0;
	switch (InterAmpFlag){
        case 122022:{
            double eta = pWF->eta;
            double delta = pWF->delta;
            double S = pWF->STotR;
            double chidiff = pWF->dchi_half;
            double eta1 = eta;
            double eta2 = eta1 * eta1;
            double eta3 = eta1 * eta2;
            double eta4 = eta1 * eta3;
            double eta5 = eta1 * eta4;
            double S1 = S;
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            total = fabs(delta*eta1*(13.318990196097973 - 21.755549987331054*eta1 + 76.14884211156267*eta2 - 127.62161159798488*eta3) + chidiff1*delta*eta1*(17.704321326939414*eta1 - 434.4390350012534*eta2 + 1366.2408490833282*eta3) + chidiff1*delta*eta1*(11.877985158418596*eta1 - 131.04937626836355*eta2 + 343.79587860999874*eta3)*S1 + chidiff1*eta5*(-1522.8543551416456 - 16.639896279650678*S1 + 3.0053086651515843*S2) + delta*eta1*S1*(-8.665646058245033*(0.7862132291286934 + 8.293609541933655*eta1 - 111.70764910503321*eta2 + 576.7172598056907*eta3 - 1001.2370065269745*eta4) - 0.9459820574514348*(1.309016452198605 + 48.94077040282239*eta1 - 817.7854010574645*eta2 + 4331.56002883546*eta3 - 7518.309520232795*eta4)*S1 - 0.4308267743835775*(9.970654092010587 - 302.9708323417439*eta1 + 3662.099161055873*eta2 - 17712.883990278668*eta3 + 29480.158198408903*eta4)*S2));
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_21_int3:version %i is not valid.", InterAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_Inter_Amp_21_int4(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag){
	double total=0;
	switch (InterAmpFlag){
        case 122022:{
            double eta = pWF->eta;
            double delta = pWF->delta;
            double S = pWF->STotR;
            double chidiff = pWF->dchi_half;
            double eta1 = eta;
            double eta2 = eta1 * eta1;
            double eta3 = eta1 * eta2;
            double eta4 = eta1 * eta3;
            double eta5 = eta1 * eta4;
            double S1 = S;
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            total = fabs(delta*eta1*(13.094382343446163 - 22.831152256559523*eta1 + 83.20619262213437*eta2 - 139.25546924151664*eta3) + chidiff1*delta*eta1*(20.120192352555357*eta1 - 458.2592421214168*eta2 + 1430.3698681181*eta3) + chidiff1*delta*eta1*(12.925363020014743*eta1 - 126.87194512915104*eta2 + 280.6003655502327*eta3)*S1 + chidiff1*eta5*(-1528.956015503355 + 74.44462583487345*S1 - 2.2456928156392197*S2) + delta*eta1*S1*(-9.499741513411829*(0.912120958549489 + 2.400945118514037*eta1 - 33.651192908287236*eta2 + 166.04254881175257*eta3 - 248.5050377498615*eta4) - 0.7850652143322492*(1.534131218043425 + 60.81773903539479*eta1 - 1032.1319480683567*eta2 + 5381.481380750608*eta3 - 9077.037917192794*eta4)*S1 - 0.21540359093306097*(9.42805409480658 - 109.06544597367301*eta1 + 385.8345793110262*eta2 + 1889.9613367802453*eta3 - 9835.416414460055*eta4)*S2));
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_21_int4:version %i is not valid.", InterAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_Inter_Amp_33_int3(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag){
	double total=0;
	switch (InterAmpFlag){
        case 122022:{
            double eta = pWF->eta;
            double delta = pWF->delta;
            double S = pWF->STotR;
            double chidiff = pWF->dchi_half;
            double eta1 = eta;
            double eta2 = eta1 * eta1;
            double eta3 = eta1 * eta2;
            double eta4 = eta1 * eta3;
            double eta5 = eta1 * eta4;
            double S1 = S;
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            total = delta*eta1*(14.555522136327964 - 12.799844096694798*eta1 + 16.79500349318081*eta2) + chidiff1*delta*eta1*(-16.292654447108134*eta1 + 190.3516012682791*eta2 - 562.0936797781519*eta3) + chidiff1*delta*eta1*(-7.048898856045782*eta1 + 49.941617405768135*eta2 - 73.62033985436068*eta3)*S1 + chidiff1*eta5*(263.5151703818307 + 44.408527093031566*S1 + 10.457035444964653*S2) + delta*eta1*S1*(0.4590550434774332*(3.0594364612798635 + 207.74562213604057*eta1 - 5545.0086137386525*eta2 + 50003.94075934942*eta3 - 195187.55422847517*eta4 + 282064.174913521*eta5) + 0.657748992123043*(5.57939137343977 - 124.06189543062042*eta1 + 1276.6209573025596*eta2 - 6999.7659193505915*eta3 + 19714.675715229736*eta4 - 20879.999628681435*eta5)*S1 + 0.3695850566805098*(6.077183107132255 - 498.95526910874986*eta1 + 10426.348944657859*eta2 - 91096.64982858274*eta3 + 360950.6686625352*eta4 - 534437.8832860565*eta5)*S2);
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_33_int3:version %i is not valid.", InterAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_Inter_Amp_33_int4(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag){
	double total=0;
	switch (InterAmpFlag){
        case 122022:{
            double eta = pWF->eta;
            double delta = pWF->delta;
            double S = pWF->STotR;
            double chidiff = pWF->dchi_half;
            double eta1 = eta;
            double eta2 = eta1 * eta1;
            double eta3 = eta1 * eta2;
            double eta4 = eta1 * eta3;
            double eta5 = eta1 * eta4;
            double S1 = S;
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = delta*eta1*(13.312095699772305 - 7.449975618083432*eta1 + 17.098576301150125*eta2) + delta*eta1*(chidiff1*(-31.171150896110156*eta1 + 371.1389274783572*eta2 - 1103.1917047361735*eta3) + chidiff2*(32.78644599730888*eta1 - 395.15713118955387*eta2 + 1164.9282236341376*eta3)) + chidiff1*delta*eta1*(-46.85669289852532*eta1 + 522.3965959942979*eta2 - 1485.5134187612182*eta3)*S1 + chidiff1*eta5*(287.90444670305715 - 21.102665129433042*chidiff2 + 7.635582066682054*S1 - 29.471275170013012*S2) + delta*eta1*S1*(0.6893003654021495*(3.1014226377197027 - 44.83989278653052*eta1 + 565.3767256471909*eta2 - 4797.429130246123*eta3 + 19514.812242035154*eta4 - 27679.226582207506*eta5) + 0.7068016563068026*(4.071212304920691 - 118.51094098279343*eta1 + 1788.1730303291356*eta2 - 13485.270489656365*eta3 + 48603.96661003743*eta4 - 65658.74746265226*eta5)*S1 + 0.2181399561677432*(-1.6754158383043574 + 303.9394443302189*eta1 - 6857.936471898544*eta2 + 59288.71069769708*eta3 - 216137.90827404748*eta4 + 277256.38289831823*eta5)*S2);
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_33_int4:version %i is not valid.", InterAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_Inter_Amp_44_int3(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag){
	double total=0;
	switch (InterAmpFlag){
        case 122022:{
            double eta = pWF->eta;
            double delta = pWF->delta;
            double S = pWF->chiPNHat;
            double chidiff = pWF->dchi_half;
            double eta1 = eta;
            double eta2 = eta1 * eta1;
            double eta3 = eta1 * eta2;
            double eta4 = eta1 * eta3;
            double S1 = S;
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = eta1*(chidiff1*delta*(-0.8765502142143329*eta1 + 22.806632458441996*eta2 - 43.675503209991184*eta3) + chidiff2*(0.48698617426180074*eta1 - 4.302527065360426*eta2 + 16.18571810759235*eta3)) + eta1*(6.379772583015967 - 44.10631039734796*eta1 + 269.44092930942793*eta2 - 1285.7635006711453*eta3 + 2379.538739132234*eta4) + eta1*S1*(-0.23316184683282615*(-1.7279023138971559 - 23.606399143993716*eta1 + 409.3387618483284*eta2 - 1115.4147472977265*eta3) - 0.09653777612560172*(-5.310643306559746 - 2.1852511802701264*eta1 + 541.1248219096527*eta2 - 1815.7529908827103*eta3)*S1 - 0.060477799540741804*(-14.578189130145661 + 175.6116682068523*eta1 - 569.4799973930861*eta2 + 426.0861915646515*eta3)*S2);
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_44_int3:version %i is not valid.", InterAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_Inter_Amp_44_int4(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag){
	double total=0;
	switch (InterAmpFlag){
        case 122022:{
            double eta = pWF->eta;
            double delta = pWF->delta;
            double S = pWF->chiPNHat;
            double chidiff = pWF->dchi_half;
            double eta1 = eta;
            double eta2 = eta1 * eta1;
            double eta3 = eta1 * eta2;
            double eta4 = eta1 * eta3;
            double eta5 = eta1 * eta4;
            double S1 = S;
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = eta1*(chidiff1*delta*(-2.461738962276138*eta1 + 45.3240543970684*eta2 - 112.2714974622516*eta3) + chidiff2*(0.9158352037567031*eta1 - 8.724582331021695*eta2 + 28.44633544874233*eta3)) + eta1*(6.098676337298138 - 45.42463610529546*eta1 + 350.97192927929433*eta2 - 2002.2013283876834*eta3 + 4067.1685640401033*eta4) + eta1*S1*(-0.36068516166901304*(-2.120354236840677 - 47.56175350408845*eta1 + 1618.4222330016048*eta2 - 14925.514654896673*eta3 + 60287.45399959349*eta4 - 91269.3745059139*eta5) - 0.09635801207669747*(-11.824692837267394 + 371.7551657959369*eta1 - 4176.398139238679*eta2 + 16655.87939259747*eta3 - 4102.218189945819*eta4 - 67024.98285179552*eta5)*S1 - 0.06565232123453196*(-26.15227471380236 + 1869.0168486099005*eta1 - 33951.35186039629*eta2 + 253694.6032002248*eta3 - 845341.6001856657*eta4 + 1.0442282862506858e6*eta5)*S2);
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_44_int4:version %i is not valid.", InterAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_Inter_Amp_32_int3(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag){
	double total=0;
	switch (InterAmpFlag){
        case 122022:{
            double eta = pWF->eta;
            double delta = pWF->delta;
            double S = pWF->chiPNHat;
            double chidiff = pWF->dchi_half;
            double eta1 = eta;
            double eta2 = eta1 * eta1;
            double eta3 = eta1 * eta2;
            double eta4 = eta1 * eta3;
            double eta5 = eta1 * eta4;
            double eta6 = eta1 * eta5;
            double S1 = S;
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = 3.881450518842405*eta1 - 12.580316392558837*eta2 + 1.7262466525848588*eta3 + chidiff2*(-7.065118823041031*eta2 + 77.97950589523865*eta3 - 203.65975422378446*eta4) - 58.408542930248046*eta4 + chidiff1*delta*(1.924723094787216*eta2 - 90.92716917757797*eta3 + 387.00162600306226*eta4) + 403.5748987560612*eta5 + chidiff1*delta*(-0.2566958540737833*eta2 + 14.488550203412675*eta3 - 26.46699529970884*eta4)*S1 + S1*(0.3650871458400108*(71.57390929624825*eta2 - 994.5272351916166*eta3 + 6734.058809060536*eta4 - 18580.859291282686*eta5 + 16001.318492586077*eta6) + 0.0960146077440495*(451.74917589707513*eta2 - 9719.470997418284*eta3 + 83403.5743434538*eta4 - 318877.43061174755*eta5 + 451546.88775684836*eta6)*S1 - 0.03985156529181297*(-304.92981902871617*eta2 + 3614.518459296278*eta3 - 7859.4784979916085*eta4 - 46454.57664737511*eta5 + 162398.81483375572*eta6)*S2);
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_32_int3:version %i is not valid.", InterAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_Inter_Amp_32_int4(IMRPhenomXWaveformStruct *pWF, int InterAmpFlag){
	double total=0;
	switch (InterAmpFlag){
        case 122022:{
            double eta = pWF->eta;
            double delta = pWF->delta;
            double S = pWF->STotR;
            double chidiff = pWF->dchi_half;
            double eta1 = eta;
            double eta2 = eta1 * eta1;
            double eta3 = eta1 * eta2;
            double eta4 = eta1 * eta3;
            double S1 = S;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = eta1*(chidiff2*(-8.572797326909152*eta1 + 92.95723645687826*eta2 - 236.2438921965621*eta3) + chidiff1*delta*(6.674358856924571*eta1 - 171.4826985994883*eta2 + 645.2760206304703*eta3)) + eta1*(3.921660532875504 - 16.57299637423352*eta1 + 25.254017911686333*eta2 - 143.41033155133266*eta3 + 692.926425981414*eta4) + chidiff1*delta*eta1*(-3.582040878719185*eta1 + 57.75888914133383*eta2 - 144.21651114700492*eta3)*S1 + eta1*S1*(1.242750265695504*(-0.522172424518215 + 25.168480118950065*eta1 - 303.5223688400309*eta2 + 1858.1518762309654*eta3 - 3797.3561904195085*eta4) + 0.2927045241764365*(0.5056957789079993 - 15.488754837330958*eta1 + 471.64047356915603*eta2 - 3131.5783196211587*eta3 + 6097.887891566872*eta4)*S1);
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Amp_32_int4:version %i is not valid.", InterAmpFlag);}
  }
  return total;
}

/* End of Amp Parameter Space Fits */


/* Solves system of equations for 5th order polynomial ansatz */

/* We use a 5th order polynomial to connect the inspiral and ringdown regions.
This polynomial is built such that it has the same value and derivative at the boundaries of the
intermediate region than the inspiral and ringdown part (demanding continuity for the function
and its first derivative) and also crosses through the two intermediate collocation points.
Sometimes we will not use the 5th order but a lower one to assure a better behaviour.

Now we have the functions to get the coefficients of such a polynomial:
      delta0 + delta1*f + delta2*f^2 + delta3*f^3 + delta4*f^4 + delta5*f^5

The different cases inside one function correspond to the different ways of building the polynomial,
and it depend on the number of collocation points and derivatives we are using.
For example case:105 correspond to the 5th order polynomial we have described before.
*/

static double IMRPhenomXHM_Intermediate_Amp_delta0(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag)
{
  double retVal;

  switch (IntAmpFlag)
  {
    case 101: //linear, only v1, v2
    {
      double f1mf4 = f1-f4;

      retVal = (-(f4*v1) + f1*v4)/f1mf4;
      break;
    }
    case 102: //quadratic: v1, v2, d2
    {
      double f12    = f1*f1;
      double f42    = f4*f4;
      double f1mf4  = f1-f4;
      double f1mf42 = f1mf4*f1mf4;

      retVal = (-(d4*f1*f1mf4*f4) + f42*v1 + f12*v4 - 2*f1*f4*v4)/f1mf42;
      break;
    }
    case 1032:  // 2 freqs, points and derivatives: v1, v4, d1, d4
    {
      double f12 = f1*f1;
      double f13 = f12*f1;
      double f42 = f4*f4;
      double f43 = f42*f4;

      double f1mf4  = f1-f4;
      double f1mf42 = f1mf4*f1mf4;
      double f1mf43 = f1mf42*f1mf4;

      retVal = (d4*f12*f4*(-f1 + f4) + d1*f1*(-f1 + f4)*f42 + 3*f1*f42*v1 - f43*v1 + f13*v4 - 3*f12*f4*v4)/f1mf43;
      break;
    }
    case 103:   // 4 freqs, no boundaries derivatives
    {
      double f12 = f1*f1;
      double f13 = f12*f1;

      double f22 = f2*f2;
      double f23 = f22*f2;

      double f32 = f3*f3;
      double f33 = f32*f3;

      double f42 = f4*f4;
      double f43 = f42*f4;

      double f1mf2 = f1-f2;
      double f1mf3 = f1-f3;
      double f1mf4 = f1-f4;
      double f2mf3 = f2-f3;
      double f2mf4 = f2-f4;
      double f3mf4 = f3-f4;

      retVal = (f1*f1mf3*f1mf4*f3*f3mf4*f4*v2 + f23*(f1*f1mf4*f4*v3 + f32*(-(f4*v1) + f1*v4) + f3*(f42*v1 - f12*v4)) + f2*(f12*f1mf4*f42*v3 + f33*(-(f42*v1) + f12*v4) + f32*(f43*v1 - f13*v4)) +
      f22*(f1*f4*(-f12 + f42)*v3 + f33*(f4*v1 - f1*v4) + f3*(-(f43*v1) + f13*v4)))/(f1mf2*f1mf3*f1mf4*f2mf3*f2mf4*f3mf4);
      break;
    }
    case 1043:  //no left derivative
    {
      double f12 = f1*f1;
      double f13 = f12*f1;
      double f14 = f13*f1;

      double f22 = f2*f2;
      double f23 = f22*f2;
      double f24 = f23*f2;

      double f32 = f3*f3;
      double f33 = f32*f3;
      double f34 = f33*f3;

      double f42 = f4*f4;
      double f43 = f42*f4;
      double f44 = f43*f4;
      double f45 = f44*f4;

      double f1mf2 = f1-f2;
      double f1mf3 = f1-f3;
      double f1mf4 = f1-f4;
      double f2mf3 = f2-f3;
      double f2mf4 = f2-f4;
      double f3mf4 = f3-f4;

      double f1mf42 = f1mf4*f1mf4;
      double f3mf42 = f3mf4*f3mf4;
      double f2mf42 = f2mf4*f2mf4;

      retVal = (-(d4*f1*f1mf2*f1mf3*f1mf4*f2*f2mf3*f2mf4*f3*f3mf4*f4) - f1*f1mf3*f1mf42*f3*f3mf42*f42*v2 +
      f24*(-(f1*f1mf42*f42*v3) + f33*(f42*v1 + f12*v4 - 2*f1*f4*v4) + f3*f4*(f43*v1 + 2*f13*v4 - 3*f12*f4*v4) - f32*(2*f43*v1 + f13*v4 - 3*f1*f42*v4)) +
      f2*f4*(f12*f1mf42*f43*v3 - f34*(f43*v1 + 2*f13*v4 - 3*f12*f4*v4) - f32*f4*(f44*v1 + 3*f14*v4 - 4*f13*f4*v4) + 2*f33*(f44*v1 + f14*v4 - 2*f12*f42*v4)) +
      f22*(-(f1*f1mf42*(2*f1 + f4)*f43*v3) + f3*f42*(f44*v1 + 3*f14*v4 - 4*f13*f4*v4) + f34*(2*f43*v1 + f13*v4 - 3*f1*f42*v4) - f33*(3*f44*v1 + f14*v4 - 4*f1*f43*v4)) +
      f23*(f1*f1mf42*(f1 + 2*f4)*f42*v3 - f34*(f42*v1 + f12*v4 - 2*f1*f4*v4) + f32*(3*f44*v1 + f14*v4 - 4*f1*f43*v4) - 2*f3*(f45*v1 + f14*f4*v4 - 2*f12*f43*v4)))/(f1mf2*f1mf3*f1mf42*f2mf3*f2mf42*f3mf42);
      break;
    }
    case 1042:   //4th order poly: v1,d1, v4,d4, v3  // used for the first intermediate region
    {

      double f12 = f1*f1;
      double f13 = f12*f1;
      double f14 = f13*f1;
      double f15 = f14*f1;

      double f42 = f4*f4;
      double f43 = f42*f4;
      double f44 = f43*f4;
      double f45 = f44*f4;

      double f32 = f3*f3;
      double f33 = f32*f3;
      double f34 = f33*f3;

      double f1mf4 = f1-f4;
      double f1mf3 = f1-f3;
      double f3mf4 = f3-f4;

      double f1mf42 = f1mf4*f1mf4;
      double f1mf32 = f1mf3*f1mf3;
      double f3mf42 = f3mf4*f3mf4;

      double f1mf43 = f1mf42*f1mf4;

      retVal = (-(d4*f12*f1mf32*f1mf4*f3*f3mf4*f4) + d1*f1*f1mf3*f1mf4*f3*f3mf42*f42 - 4*f12*f33*f42*v1 + 3*f1*f34*f42*v1 + 8*f12*f32*f43*v1 - 4*f1*f33*f43*v1 - f34*f43*v1 - 4*f12*f3*f44*v1 - f1*f32*f44*v1 +
      2*f33*f44*v1 + 2*f1*f3*f45*v1 - f32*f45*v1 + f15*f42*v3 - 3*f14*f43*v3 + 3*f13*f44*v3 - f12*f45*v3 + f15*f32*v4 - 2*f14*f33*v4 + f13*f34*v4 - 2*f15*f3*f4*v4 + f14*f32*f4*v4 + 4*f13*f33*f4*v4 -
      3*f12*f34*f4*v4 + 4*f14*f3*f42*v4 - 8*f13*f32*f42*v4 + 4*f12*f33*f42*v4)/(f1mf32*f1mf43*f3mf42);

      break;
    }
    case 104:  //Geraint's Version, 4th order poly: v1,d1, v4,d4, v2
    {

      double f12 = f1*f1;

      double f42 = f4*f4;

      double f1mf2 = f1-f2;
      double f1mf4 = f1-f4;
      double f2mf4 = f2-f4;

      double f1mf22 = f1mf2*f1mf2;
      double f2mf42 = f2mf4*f2mf4;
      double f1mf43 = f1mf4*f1mf4*f1mf4;

      retVal = ((-(d4*f12*f1mf22*f1mf4*f2*f2mf4*f4) + d1*f1*f1mf2*f1mf4*f2*f2mf42*f42 + f42*(f2*f2mf42*(-4*f12 + 3*f1*f2 + 2*f1*f4 - f2*f4)*v1 + f12*f1mf43*v2) +
      f12*f1mf22*f2*(f1*f2 - 2*f1*f4 - 3*f2*f4 + 4*f42)*v4)/(f1mf22*f1mf43*f2mf42));
      break;
    }
    case 105: // Geraint, standard way: v1, v2, v3, v4, d1, d4
    {
      double f12 = f1*f1;
      double f13 = f12*f1;
      double f14 = f13*f1;
      double f15 = f14*f1;
      double f16 = f15*f1;
      double f17 = f16*f1;

      double f22 = f2*f2;
      double f23 = f22*f2;
      double f24 = f23*f2;
      double f25 = f24*f2;

      double f32 = f3*f3;
      double f33 = f32*f3;
      double f34 = f33*f3;
      double f35 = f34*f3;

      double f42 = f4*f4;
      double f43 = f42*f4;
      double f44 = f43*f4;
      double f45 = f44*f4;
      double f46 = f45*f4;
      double f47 = f46*f4;

      double f1mf2 = f1-f2;
      double f1mf3 = f1-f3;
      double f1mf4 = f1-f4;
      double f2mf3 = f2-f3;
      double f2mf4 = f2-f4;
      double f3mf4 = f3-f4;

      double f1mf22 = f1mf2*f1mf2;
      double f1mf32 = f1mf3*f1mf3;
      double f1mf42 = f1mf4*f1mf4;
      double f2mf42 = f2mf4*f2mf4;
      double f3mf42 = f3mf4*f3mf4;
      double f1mf43 = f1mf42*f1mf4;

      retVal = (
        (-(d4*f12*f1mf22*f1mf32*f1mf4*f2*f2mf3*f2mf4*f3*f3mf4*f4) - d1*f1*f1mf2*f1mf3*f1mf4*f2*f2mf3*f2mf42*f3*f3mf42*f42 + 5*f13*f24*f33*f42*v1 - 4*f12*f25*f33*f42*v1 - 5*f13*f23*f34*f42*v1 +
        3*f1*f25*f34*f42*v1 + 4*f12*f23*f35*f42*v1 - 3*f1*f24*f35*f42*v1 - 10*f13*f24*f32*f43*v1 + 8*f12*f25*f32*f43*v1 + 5*f12*f24*f33*f43*v1 - 4*f1*f25*f33*f43*v1 + 10*f13*f22*f34*f43*v1 -
        5*f12*f23*f34*f43*v1 - f25*f34*f43*v1 - 8*f12*f22*f35*f43*v1 + 4*f1*f23*f35*f43*v1 + f24*f35*f43*v1 + 5*f13*f24*f3*f44*v1 - 4*f12*f25*f3*f44*v1 + 15*f13*f23*f32*f44*v1 -
        10*f12*f24*f32*f44*v1 - f1*f25*f32*f44*v1 - 15*f13*f22*f33*f44*v1 + 5*f1*f24*f33*f44*v1 + 2*f25*f33*f44*v1 - 5*f13*f2*f34*f44*v1 + 10*f12*f22*f34*f44*v1 - 5*f1*f23*f34*f44*v1 +
        4*f12*f2*f35*f44*v1 + f1*f22*f35*f44*v1 - 2*f23*f35*f44*v1 - 10*f13*f23*f3*f45*v1 + 5*f12*f24*f3*f45*v1 + 2*f1*f25*f3*f45*v1 - f12*f23*f32*f45*v1 + 2*f1*f24*f32*f45*v1 - f25*f32*f45*v1 +
        10*f13*f2*f33*f45*v1 + f12*f22*f33*f45*v1 - 3*f24*f33*f45*v1 - 5*f12*f2*f34*f45*v1 - 2*f1*f22*f34*f45*v1 + 3*f23*f34*f45*v1 - 2*f1*f2*f35*f45*v1 + f22*f35*f45*v1 + 5*f13*f22*f3*f46*v1 +
        2*f12*f23*f3*f46*v1 - 4*f1*f24*f3*f46*v1 - 5*f13*f2*f32*f46*v1 - f1*f23*f32*f46*v1 + 2*f24*f32*f46*v1 - 2*f12*f2*f33*f46*v1 + f1*f22*f33*f46*v1 + 4*f1*f2*f34*f46*v1 - 2*f22*f34*f46*v1 -
        3*f12*f22*f3*f47*v1 + 2*f1*f23*f3*f47*v1 + 3*f12*f2*f32*f47*v1 - f23*f32*f47*v1 - 2*f1*f2*f33*f47*v1 + f22*f33*f47*v1 - f17*f33*f42*v2 + 2*f16*f34*f42*v2 - f15*f35*f42*v2 + 2*f17*f32*f43*v2 -
        f16*f33*f43*v2 - 4*f15*f34*f43*v2 + 3*f14*f35*f43*v2 - f17*f3*f44*v2 - 4*f16*f32*f44*v2 + 8*f15*f33*f44*v2 - 3*f13*f35*f44*v2 + 3*f16*f3*f45*v2 - 8*f14*f33*f45*v2 + 4*f13*f34*f45*v2 +
        f12*f35*f45*v2 - 3*f15*f3*f46*v2 + 4*f14*f32*f46*v2 + f13*f33*f46*v2 - 2*f12*f34*f46*v2 + f14*f3*f47*v2 - 2*f13*f32*f47*v2 + f12*f33*f47*v2 + f17*f23*f42*v3 - 2*f16*f24*f42*v3 +
        f15*f25*f42*v3 - 2*f17*f22*f43*v3 + f16*f23*f43*v3 + 4*f15*f24*f43*v3 - 3*f14*f25*f43*v3 + f17*f2*f44*v3 + 4*f16*f22*f44*v3 - 8*f15*f23*f44*v3 + 3*f13*f25*f44*v3 - 3*f16*f2*f45*v3 +
        8*f14*f23*f45*v3 - 4*f13*f24*f45*v3 - f12*f25*f45*v3 + 3*f15*f2*f46*v3 - 4*f14*f22*f46*v3 - f13*f23*f46*v3 + 2*f12*f24*f46*v3 - f14*f2*f47*v3 + 2*f13*f22*f47*v3 - f12*f23*f47*v3 +
        f12*f1mf22*f1mf32*f2*f2mf3*f3*(f4*(-3*f2*f3 + 4*(f2 + f3)*f4 - 5*f42) + f1*(f2*f3 - 2*(f2 + f3)*f4 + 3*f42))*v4)/(f1mf22*f1mf32*f1mf43*f2mf3*f2mf42*f3mf42)
      );
      break;
    }
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomXHM_Intermediate_Amp_delta0: IMRPhenomXIntermediateAmpVersion is not valid.\n");
    }
  }

  return retVal;
}

static double IMRPhenomXHM_Intermediate_Amp_delta1(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag)
{
  double retVal;
  switch ( IntAmpFlag )
  {
    case 101: //linear, only v1, v2
    {
      double f1mf4 = f1-f4;

      retVal = (v1 - v4)/f1mf4;
      break;
    }
    case 102: //quadratic: v1, v2, d2
    {
      double f12 = f1*f1;
      double f42 = f4*f4;
      double f1mf42 = (f1-f4)*(f1-f4);

      retVal = (d4*(f12 - f42) + 2*f4*(-v1 + v4))/f1mf42;
      break;
    }
    case 1032:  // 2 freqs, points and derivatives: v1, v4, d1, d4
    {
      double f12 = f1*f1;
      double f42 = f4*f4;

      double f1mf4  = f1-f4;
      double f1mf42 = f1mf4*f1mf4;
      double f1mf43 = f1mf42*f1mf4;

      retVal = (d4*f1*f1mf4*(f1 + 2*f4) - f4*(d1*(-2*f12 + f1*f4 + f42) + 6*f1*(v1 - v4)))/f1mf43;
      break;
    }
    case 103:   // 4 freqs, no boundaries derivatives
    {
      double f12 = f1*f1;
      double f13 = f12*f1;

      double f22 = f2*f2;
      double f23 = f22*f2;

      double f32 = f3*f3;
      double f33 = f32*f3;

      double f42 = f4*f4;
      double f43 = f42*f4;

      double f1mf2 = f1-f2;
      double f1mf3 = f1-f3;
      double f1mf4 = f1-f4;
      double f2mf3 = f2-f3;
      double f2mf4 = f2-f4;
      double f3mf4 = f3-f4;

      retVal = (f12*f1mf4*f42*(v2 - v3) + f33*(f42*(v1 - v2) + f12*(v2 - v4)) + f22*(f43*(v1 - v3) + f13*(v3 - v4) + f33*(-v1 + v4)) + f32*(f43*(-v1 + v2) + f13*(-v2 + v4)) +
      f23*(f42*(-v1 + v3) + f32*(v1 - v4) + f12*(-v3 + v4)))/(f1mf2*f1mf3*f1mf4*f2mf3*f2mf4*f3mf4);
      break;
    }
    case 1043:  //no left derivative
    {
      double f12 = f1*f1;
      double f13 = f12*f1;
      double f14 = f13*f1;

      double f22 = f2*f2;
      double f23 = f22*f2;
      double f24 = f23*f2;

      double f32 = f3*f3;
      double f33 = f32*f3;
      double f34 = f33*f3;

      double f42 = f4*f4;
      double f43 = f42*f4;
      double f44 = f43*f4;

      double f1mf2 = f1-f2;
      double f1mf3 = f1-f3;
      double f1mf4 = f1-f4;
      double f2mf3 = f2-f3;
      double f2mf4 = f2-f4;
      double f3mf4 = f3-f4;

      double f1mf42 = f1mf4*f1mf4;
      double f2mf42 = f2mf4*f2mf4;
      double f3mf42 = f3mf4*f3mf4;

      double v1mv2 = v1-v2;
      double v2mv3 = v2-v3;
      double v2mv4 = v2-v4;
      double v1mv3 = v1-v3;
      double v1mv4 = v1-v4;
      double v3mv4 = v3-v4;

      retVal =(d4*f1mf4*f2mf4*f3mf4*(f1*f2*f3 + f2*f3*f4 + f1*(f2 + f3)*f4) + (f4*(f12*f1mf42*f43*v2mv3 + f34*(f43*v1mv2 +
        3*f12*f4*v2mv4 + 2*f13*(-v2 + v4)) + f32*f4*(f44*v1mv2 + 4*f13*f4*v2mv4 + 3*f14*(-v2 + v4)) + 2*f33*(f44*(-v1 + v2)
        + f14*v2mv4 + 2*f12*f42*(-v2 + v4)) + 2*f23*(f44*v1mv3 + f34*v1mv4 + 2*f12*f42*v3mv4 + 2*f32*f42*(-v1 + v4) + f14*(-v3 + v4))
        + f24*(3*f32*f4*v1mv4 + f43*(-v1 + v3) + 2*f13*v3mv4 + 2*f33*(-v1 + v4) + 3*f12*f4*(-v3 + v4)) + f22*f4*(4*f33*f4*v1mv4 + f44*(-v1 + v3)
        + 3*f14*v3mv4 + 3*f34*(-v1 + v4) + 4*f13*f4*(-v3 + v4))))/(f1mf2*f1mf3*f2mf3))/(f1mf42*f2mf42*f3mf42);

      break;
    }
    case 1042:   //4th order poly: v1,d1, v4,d4, v3  // used for the first intermediate region
    {
      double f12 = f1*f1;
      double f13 = f12*f1;
      double f14 = f13*f1;
      double f15 = f14*f1;

      double f42 = f4*f4;
      double f43 = f42*f4;
      double f44 = f43*f4;
      double f45 = f44*f4;

      double f32 = f3*f3;
      double f33 = f32*f3;
      double f34 = f33*f3;

      double f1mf4 = f1-f4;
      double f1mf3 = f1-f3;
      double f3mf4 = f3-f4;

      double f1mf42 = f1mf4*f1mf4;
      double f1mf32 = f1mf3*f1mf3;
      double f3mf42 = f3mf4*f3mf4;

      double f1mf43 = f1mf42*f1mf4;

      retVal = (d4*f15*f32 - 2*d4*f14*f33 + d4*f13*f34 + d4*f14*f32*f4 - 2*d1*f13*f33*f4 - 2*d4*f13*f33*f4 + 2*d1*f12*f34*f4 + d4*f12*f34*f4 - d4*f15*f42 + 3*d1*f13*f32*f42 + d4*f13*f32*f42 - 2*d1*f12*f33*f42 +
        2*d4*f12*f33*f42 - d1*f1*f34*f42 - 2*d4*f1*f34*f42 + d4*f14*f43 - d1*f12*f32*f43 - 3*d4*f12*f32*f43 + 2*d1*f1*f33*f43 + 2*d4*f1*f33*f43 - d1*f34*f43 - d1*f13*f44 - d1*f1*f32*f44 + 2*d1*f33*f44 +
        d1*f12*f45 - d1*f32*f45 + 8*f12*f33*f4*v1 - 6*f1*f34*f4*v1 - 12*f12*f32*f42*v1 + 8*f1*f33*f42*v1 + 4*f12*f44*v1 - 2*f1*f45*v1 - 2*f15*f4*v3 + 4*f14*f42*v3 - 4*f12*f44*v3 + 2*f1*f45*v3 +
        2*f15*f4*v4 - 8*f12*f33*f4*v4 + 6*f1*f34*f4*v4 - 4*f14*f42*v4 + 12*f12*f32*f42*v4 - 8*f1*f33*f42*v4)/(f1mf32*f1mf43*f3mf42);

        #if DEBUG == 1
        printf("\ndelta1 = %.16f", retVal);
        printf("\nf1 = %.16f", f1);
        printf("\nf2 = %.16f", f2);
        printf("\nf3 = %.16f", f3);
        printf("\nf4 = %.16f", f4);
        printf("\nv1 = %.16f", v1);
        printf("\nv2 = %.16f", v2);
        printf("\nv3 = %.16f", v3);
        printf("\nv4 = %.16f", v4);
        printf("\nd1 = %.16f", d1);
        printf("\nd4 = %.16f", d4);
        #endif

        break;
      }
      case 104:  //Geraint's Version, 4th order poly: v1,d1, v2,d2, v3
      {

        double f12 = f1*f1;
        double f13 = f12*f1;
        double f14 = f13*f1;

        double f22 = f2*f2;
        double f23 = f22*f2;
        double f24 = f23*f2;

        double f42 = f4*f4;
        double f43 = f42*f4;
        double f44 = f43*f4;

        double f1mf2 = f1-f2;
        double f1mf4 = f1-f4;
        double f2mf4 = f2-f4;

        double f1mf22 = f1mf2*f1mf2;
        double f2mf42 = f2mf4*f2mf4;
        double f1mf43 = f1mf4*f1mf4*f1mf4;

        retVal = ((d4*f1*f1mf22*f1mf4*f2mf4*(2*f2*f4 + f1*(f2 + f4)) + f4*(-(d1*f1mf2*f1mf4*f2mf42*(2*f1*f2 + (f1 + f2)*f4)) -
        2*f1*(f44*(v1 - v2) + 3*f24*(v1 - v4) + f14*(v2 - v4) + 4*f23*f4*(-v1 + v4)
        + 2*f13*f4*(-v2 + v4) + f1*(2*f43*(-v1 + v2) + 6*f22*f4*(v1 - v4) + 4*f23*(-v1 + v4)))))/(f1mf22*f1mf43*f2mf42));
        break;
      }
      case 105: // Geraint, standard way: v1, v2, v3, v4, d1, d4
      {

        double f12 = f1*f1;
        double f13 = f12*f1;
        double f14 = f13*f1;
        double f15 = f14*f1;
        double f16 = f15*f1;

        double f22 = f2*f2;
        double f23 = f22*f2;
        double f24 = f23*f2;
        double f25 = f24*f2;

        double f32 = f3*f3;
        double f33 = f32*f3;
        double f34 = f33*f3;
        double f35 = f34*f3;

        double f42 = f4*f4;
        double f43 = f42*f4;
        double f44 = f43*f4;
        double f45 = f44*f4;

        double f1mf2 = f1-f2;
        double f1mf3 = f1-f3;
        double f1mf4 = f1-f4;
        double f2mf3 = f2-f3;
        double f2mf4 = f2-f4;
        double f3mf4 = f3-f4;

        double f1mf22 = f1mf2*f1mf2;
        double f1mf32 = f1mf3*f1mf3;
        double f2mf42 = f2mf4*f2mf4;
        double f3mf42 = f3mf4*f3mf4;
        double f1mf43 = f1mf4*f1mf4*f1mf4;

        retVal = (
          (d4*f1*f1mf22*f1mf32*f1mf4*f2mf3*f2mf4*f3mf4*(f1*f2*f3 + 2*f2*f3*f4 + f1*(f2 + f3)*f4) +
          f4*(d1*f1mf2*f1mf3*f1mf4*f2mf3*f2mf42*f3mf42*(2*f1*f2*f3 + f2*f3*f4 + f1*(f2 + f3)*f4) +
          f1*(f16*(f43*(v2 - v3) + 2*f33*(v2 - v4) + 3*f22*f4*(v3 - v4) + 3*f32*f4*(-v2 + v4) + 2*f23*(-v3 + v4)) +
          f13*f4*(f45*(-v2 + v3) + 5*f34*f4*(v2 - v4) + 4*f25*(v3 - v4) + 4*f35*(-v2 + v4) + 5*f24*f4*(-v3 + v4)) +
          f14*(3*f45*(v2 - v3) + 2*f35*(v2 - v4) + 5*f34*f4*(v2 - v4) + 10*f23*f42*(v3 - v4) + 10*f33*f42*(-v2 + v4) + 2*f25*(-v3 + v4) + 5*f24*f4*(-v3 + v4)) +
          f15*(3*f44*(-v2 + v3) + 2*f33*f4*(v2 - v4) + 5*f32*f42*(v2 - v4) + 4*f24*(v3 - v4) + 4*f34*(-v2 + v4) + 2*f23*f4*(-v3 + v4) + 5*f22*f42*(-v3 + v4)) -
          5*f12*(-(f32*f3mf42*f43*(v1 - v2)) + 2*f23*(f44*(-v1 + v3) + 2*f32*f42*(v1 - v4) + f34*(-v1 + v4)) + f24*(f43*(v1 - v3) + 2*f33*(v1 - v4) + 3*f32*f4*(-v1 + v4)) +
          f22*f4*(f44*(v1 - v3) + 3*f34*(v1 - v4) + 4*f33*f4*(-v1 + v4))) +
          f1*(-(f32*f3mf42*(4*f3 + 3*f4)*f43*(v1 - v2)) + 2*f23*(f45*(-v1 + v3) + 5*f34*f4*(v1 - v4) + 4*f35*(-v1 + v4)) + 4*f25*(f43*(v1 - v3) + 2*f33*(v1 - v4) + 3*f32*f4*(-v1 + v4)) -
          5*f24*f4*(f43*(v1 - v3) + 2*f33*(v1 - v4) + 3*f32*f4*(-v1 + v4)) + 3*f22*f4*(f45*(v1 - v3) + 4*f35*(v1 - v4) + 5*f34*f4*(-v1 + v4))) -
          2*(-(f33*f3mf42*f44*(v1 - v2)) + f24*(2*f45*(-v1 + v3) + 5*f33*f42*(v1 - v4) + 3*f35*(-v1 + v4)) + f25*(f44*(v1 - v3) + 3*f34*(v1 - v4) + 4*f33*f4*(-v1 + v4)) +
          f23*f4*(f45*(v1 - v3) + 4*f35*(v1 - v4) + 5*f34*f4*(-v1 + v4))))))/(f1mf22*f1mf32*f1mf43*f2mf3*f2mf42*f3mf42)
        );
        break;
      }
      default:
      {
        XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomXHM_Intermediate_Amp_delta1: IMRPhenomXIntermediateAmpVersion is not valid.\n");
      }
    }
    return retVal;
  }

static double IMRPhenomXHM_Intermediate_Amp_delta2(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag)
{
    double retVal;

    switch ( IntAmpFlag )
    {
      case 101: //linear, only v1, v2
      {
        retVal = 0.;
        break;
      }
      case 102: //quadratic: v1, v2, d2
      {
        double f1mf4  = f1-f4;
        double f1mf42 = f1mf4*f1mf4;

        retVal = (-(d4*f1mf4) + v1 - v4)/f1mf42;
        break;
      }
      case 1032:  // 2 freqs, points and derivatives: v1, v4, d1, d4
      {
        double f12 = f1*f1;
        double f42 = f4*f4;

        double f1mf4  = f1-f4;
        double f1mf42 = f1mf4*f1mf4;
        double f1mf43 = f1mf42*f1mf4;

        retVal = (-(d1*(f12 + f1*f4 - 2*f42)) + d4*(-2*f12 + f1*f4 + f42) + 3*(f1 + f4)*(v1 - v4))/f1mf43;
        break;
      }
      case 103:   // 4 freqs, no boundaries derivatives
      {
        double f12 = f1*f1;
        double f13 = f12*f1;

        double f22 = f2*f2;
        double f23 = f22*f2;

        double f32 = f3*f3;
        double f33 = f32*f3;

        double f42 = f4*f4;
        double f43 = f42*f4;

        double f1mf2 = f1-f2;
        double f1mf3 = f1-f3;
        double f1mf4 = f1-f4;
        double f2mf3 = f2-f3;
        double f2mf4 = f2-f4;
        double f3mf4 = f3-f4;

        retVal = (-(f1*f4*(f12 - f42)*(v2 - v3)) + f3*(f43*(v1 - v2) + f13*(v2 - v4)) + f23*(f4*(v1 - v3) + f1*(v3 - v4) + f3*(-v1 + v4)) + f33*(f4*(-v1 + v2) + f1*(-v2 + v4)) +
        f2*(f43*(-v1 + v3) + f33*(v1 - v4) + f13*(-v3 + v4)))/(f1mf2*f1mf3*f1mf4*f2mf3*f2mf4*f3mf4);
        break;
      }
      case 1043:  //no left derivative: v1, v2, v3, v4, d4
      {
        double f12 = f1*f1;
        double f13 = f12*f1;
        double f14 = f13*f1;

        double f22 = f2*f2;
        double f23 = f22*f2;
        double f24 = f23*f2;

        double f32 = f3*f3;
        double f33 = f32*f3;
        double f34 = f33*f3;

        double f42 = f4*f4;
        double f43 = f42*f4;
        double f44 = f43*f4;
        double f46 = f44*f42;

        double f1mf2 = f1-f2;
        double f1mf3 = f1-f3;
        double f1mf4 = f1-f4;
        double f2mf3 = f2-f3;
        double f2mf4 = f2-f4;
        double f3mf4 = f3-f4;

        double f1mf42 = f1mf4*f1mf4;
        double f2mf42 = f2mf4*f2mf4;
        double f3mf42 = f3mf4*f3mf4;

        double v1mv3 = v1-v3;
        double v1mv4 = v1-v4;
        double v3mv4 = v3-v4;

        retVal = (-(d4*f1mf2*f1mf3*f1mf4*f2mf3*f2mf4*f3mf4*(f3*f4 + f2*(f3 + f4) + f1*(f2 + f3 + f4))) - 2*f34*f43*v1 + 3*f33*f44*v1 - f3*f46*v1 - f14*f33*v2 + f13*f34*v2 + 3*f14*f3*f42*v2 - 3*f1*f34*f42*v2 -
        2*f14*f43*v2 - 4*f13*f3*f43*v2 + 4*f1*f33*f43*v2 + 2*f34*f43*v2 + 3*f13*f44*v2 - 3*f33*f44*v2 - f1*f46*v2 + f3*f46*v2 + 2*f14*f43*v3 - 3*f13*f44*v3 + f1*f46*v3 +
        f2*f42*(f44*v1mv3 + 3*f34*v1mv4 - 4*f33*f4*v1mv4 - 3*f14*v3mv4 + 4*f13*f4*v3mv4) + f24*(2*f43*v1mv3 + f33*v1mv4 - 3*f3*f42*v1mv4 - f13*v3mv4 + 3*f1*f42*v3mv4) +
        f23*(-3*f44*v1mv3 - f34*v1mv4 + 4*f3*f43*v1mv4 + f14*v3mv4 - 4*f1*f43*v3mv4) + f14*f33*v4 - f13*f34*v4 - 3*f14*f3*f42*v4 + 3*f1*f34*f42*v4 + 4*f13*f3*f43*v4 - 4*f1*f33*f43*v4)/
        (f1mf2*f1mf3*f1mf42*f2mf3*f2mf42*f3mf42);
        break;
      }
      case 1042:   //4th order poly: v1,d1, v2,d2, v3   // used for the first intermediate region
      {
        double f12 = f1*f1;
        double f13 = f12*f1;
        double f14 = f13*f1;
        double f15 = f14*f1;

        double f42 = f4*f4;
        double f43 = f42*f4;
        double f44 = f43*f4;
        double f45 = f44*f4;

        double f32 = f3*f3;
        double f33 = f32*f3;
        double f34 = f33*f3;

        double f1mf4 = f1-f4;
        double f1mf3 = f1-f3;
        double f3mf4 = f3-f4;

        double f1mf42 = f1mf4*f1mf4;
        double f1mf32 = f1mf3*f1mf3;
        double f3mf42 = f3mf4*f3mf4;

        double f1mf43 = f1mf42*f1mf4;

        retVal = (-(d4*f1mf32*f1mf4*f3mf4*(f12 + f3*f4 + 2*f1*(f3 + f4))) + d1*f1mf3*f1mf4*f3mf42*(f1*f3 + 2*(f1 + f3)*f4 + f42) - 4*f12*f33*v1 + 3*f1*f34*v1 - 4*f1*f33*f4*v1 + 3*f34*f4*v1 + 12*f12*f3*f42*v1 -
        4*f33*f42*v1 - 8*f12*f43*v1 + f1*f44*v1 + f45*v1 + f15*v3 + f14*f4*v3 - 8*f13*f42*v3 + 8*f12*f43*v3 - f1*f44*v3 - f45*v3 -
        f1mf32*(f13 + f3*(3*f3 - 4*f4)*f4 + f12*(2*f3 + f4) + f1*(3*f3 - 4*f4)*(f3 + 2*f4))*v4)/(f1mf32*f1mf43*f3mf42);
        break;
      }
      case 104:  //Geraint's Version, 4th order poly: v1,d1, v2,d2, v3
      {
        double f12 = f1*f1;
        double f13 = f12*f1;
        double f14 = f13*f1;
        double f15 = f14*f1;

        double f22 = f2*f2;
        double f23 = f22*f2;
        double f24 = f23*f2;

        double f42 = f4*f4;
        double f43 = f42*f4;
        double f44 = f43*f4;
        double f45 = f44*f4;

        double f1mf2 = f1-f2;
        double f1mf4 = f1-f4;
        double f2mf4 = f2-f4;

        double f1mf22 = f1mf2*f1mf2;
        double f2mf42 = f2mf4*f2mf4;
        double f1mf43 = f1mf4*f1mf4*f1mf4;

        retVal = ((-(d4*f1mf22*f1mf4*f2mf4*(f12 + f2*f4 + 2*f1*(f2 + f4))) + d1*f1mf2*f1mf4*f2mf42*(f1*f2 + 2*(f1 + f2)*f4 + f42)
        - 4*f12*f23*v1 + 3*f1*f24*v1 - 4*f1*f23*f4*v1 + 3*f24*f4*v1 + 12*f12*f2*f42*v1 -
        4*f23*f42*v1 - 8*f12*f43*v1 + f1*f44*v1 + f45*v1 + f15*v2 + f14*f4*v2 - 8*f13*f42*v2 + 8*f12*f43*v2 - f1*f44*v2 - f45*v2 -
        f1mf22*(f13 + f2*(3*f2 - 4*f4)*f4 + f12*(2*f2 + f4) + f1*(3*f2 - 4*f4)*(f2 + 2*f4))*v4)/(f1mf22*f1mf43*f2mf42));

        break;
      }
      case 105: // Geraint, standard way: v1, v2, v3, v4, d1, d4
      {
        double f12 = f1*f1;
        double f13 = f12*f1;
        double f14 = f13*f1;
        double f15 = f14*f1;
        double f16 = f15*f1;
        double f17 = f16*f1;

        double f22 = f2*f2;
        double f23 = f22*f2;
        double f24 = f23*f2;
        double f25 = f24*f2;

        double f32 = f3*f3;
        double f33 = f32*f3;
        double f34 = f33*f3;
        double f35 = f34*f3;

        double f42 = f4*f4;
        double f43 = f42*f4;
        double f44 = f43*f4;
        double f45 = f44*f4;
        double f46 = f45*f4;
        double f47 = f46*f4;

        double f1mf2 = f1-f2;
        double f1mf3 = f1-f3;
        double f1mf4 = f1-f4;
        double f2mf3 = f2-f3;
        double f2mf4 = f2-f4;
        double f3mf4 = f3-f4;

        double f1mf22 = f1mf2*f1mf2;
        double f1mf32 = f1mf3*f1mf3;
        double f2mf42 = f2mf4*f2mf4;
        double f3mf42 = f3mf4*f3mf4;
        double f1mf43 = f1mf4*f1mf4*f1mf4;

        retVal = (
          (-(d4*f1mf22*f1mf32*f1mf4*f2mf3*f2mf4*f3mf4*(f2*f3*f4 + f12*(f2 + f3 + f4) + 2*f1*(f2*f3 + (f2 + f3)*f4))) -
          d1*f1mf2*f1mf3*f1mf4*f2mf3*f2mf42*f3mf42*(f1*f2*f3 + 2*(f2*f3 + f1*(f2 + f3))*f4 + (f1 + f2 + f3)*f42) + 5*f13*f24*f33*v1 - 4*f12*f25*f33*v1 - 5*f13*f23*f34*v1 + 3*f1*f25*f34*v1 +
          4*f12*f23*f35*v1 - 3*f1*f24*f35*v1 + 5*f12*f24*f33*f4*v1 - 4*f1*f25*f33*f4*v1 - 5*f12*f23*f34*f4*v1 + 3*f25*f34*f4*v1 + 4*f1*f23*f35*f4*v1 - 3*f24*f35*f4*v1 - 15*f13*f24*f3*f42*v1 +
          12*f12*f25*f3*f42*v1 + 5*f1*f24*f33*f42*v1 - 4*f25*f33*f42*v1 + 15*f13*f2*f34*f42*v1 - 5*f1*f23*f34*f42*v1 - 12*f12*f2*f35*f42*v1 + 4*f23*f35*f42*v1 + 10*f13*f24*f43*v1 - 8*f12*f25*f43*v1 +
          20*f13*f23*f3*f43*v1 - 15*f12*f24*f3*f43*v1 - 20*f13*f2*f33*f43*v1 + 5*f24*f33*f43*v1 - 10*f13*f34*f43*v1 + 15*f12*f2*f34*f43*v1 - 5*f23*f34*f43*v1 + 8*f12*f35*f43*v1 - 15*f13*f23*f44*v1 +
          10*f12*f24*f44*v1 + f1*f25*f44*v1 + 15*f13*f33*f44*v1 - 10*f12*f34*f44*v1 - f1*f35*f44*v1 + f12*f23*f45*v1 - 2*f1*f24*f45*v1 + f25*f45*v1 - f12*f33*f45*v1 + 2*f1*f34*f45*v1 - f35*f45*v1 +
          5*f13*f2*f46*v1 + f1*f23*f46*v1 - 2*f24*f46*v1 - 5*f13*f3*f46*v1 - f1*f33*f46*v1 + 2*f34*f46*v1 - 3*f12*f2*f47*v1 + f23*f47*v1 + 3*f12*f3*f47*v1 - f33*f47*v1 - f17*f33*v2 + 2*f16*f34*v2 -
          f15*f35*v2 - f16*f33*f4*v2 + 2*f15*f34*f4*v2 - f14*f35*f4*v2 + 3*f17*f3*f42*v2 - f15*f33*f42*v2 - 10*f14*f34*f42*v2 + 8*f13*f35*f42*v2 - 2*f17*f43*v2 - 5*f16*f3*f43*v2 + 15*f14*f33*f43*v2 -
          8*f12*f35*f43*v2 + 4*f16*f44*v2 - 15*f13*f33*f44*v2 + 10*f12*f34*f44*v2 + f1*f35*f44*v2 + f12*f33*f45*v2 - 2*f1*f34*f45*v2 + f35*f45*v2 - 4*f14*f46*v2 + 5*f13*f3*f46*v2 + f1*f33*f46*v2 -
          2*f34*f46*v2 + 2*f13*f47*v2 - 3*f12*f3*f47*v2 + f33*f47*v2 + f17*f23*v3 - 2*f16*f24*v3 + f15*f25*v3 + f16*f23*f4*v3 - 2*f15*f24*f4*v3 + f14*f25*f4*v3 - 3*f17*f2*f42*v3 + f15*f23*f42*v3 +
          10*f14*f24*f42*v3 - 8*f13*f25*f42*v3 + 2*f17*f43*v3 + 5*f16*f2*f43*v3 - 15*f14*f23*f43*v3 + 8*f12*f25*f43*v3 - 4*f16*f44*v3 + 15*f13*f23*f44*v3 - 10*f12*f24*f44*v3 - f1*f25*f44*v3 -
          f12*f23*f45*v3 + 2*f1*f24*f45*v3 - f25*f45*v3 + 4*f14*f46*v3 - 5*f13*f2*f46*v3 - f1*f23*f46*v3 + 2*f24*f46*v3 - 2*f13*f47*v3 + 3*f12*f2*f47*v3 - f23*f47*v3 -
          f1mf22*f1mf32*f2mf3*(f13*(f22 + f2*f3 + f32 - 3*f42) + f2*f3*f4*(3*f2*f3 - 4*(f2 + f3)*f4 + 5*f42) + f1*(f2*f3 + 2*(f2 + f3)*f4)*(3*f2*f3 - 4*(f2 + f3)*f4 + 5*f42) +
          f12*(2*f2*f3*(f2 + f3) + (f22 + f2*f3 + f32)*f4 - 6*(f2 + f3)*f42 + 5*f43))*v4)/(f1mf22*f1mf32*f1mf43*f2mf3*f2mf42*f3mf42)
        );

        break;
      }
      default:
      {
        XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomXHM_Intermediate_Amp_delta2: IMRPhenomXIntermediateAmpVersion is not valid.\n");
      }
    }

    return retVal;
}

static double IMRPhenomXHM_Intermediate_Amp_delta3(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag)
{
    double retVal;

    switch ( IntAmpFlag )
    {
      case 101: //linear, only v1, v2
      {
        retVal = 0.;
        break;
      }
      case 102: //quadratic: v1, v2, d2
      {
        retVal = 0.;
        break;
      }
      case 1032:  // 2 freqs, points and derivatives: v1, v4, d1, d4
      {
        double f1mf4 = f1-f4;
        double f1mf42 = f1mf4*f1mf4;
        double f1mf43 = f1mf42*f1mf4;

        retVal = (d1*f1mf4 + d4*f1mf4 - 2*v1 + 2*v4)/f1mf43;
        break;
      }
      case 103:  // 4 freqs, no boundaries derivatives
      {
        double f12 = f1*f1;

        double f22 = f2*f2;

        double f32 = f3*f3;

        double f42 = f4*f4;

        double f1mf2 = f1-f2;
        double f1mf3 = f1-f3;
        double f1mf4 = f1-f4;
        double f2mf3 = f2-f3;
        double f2mf4 = f2-f4;
        double f3mf4 = f3-f4;

        retVal = (f1*f1mf4*f4*(v2 - v3) + f32*(f4*(v1 - v2) + f1*(v2 - v4)) + f2*(f42*(v1 - v3) + f12*(v3 - v4) + f32*(-v1 + v4)) + f3*(f42*(-v1 + v2) + f12*(-v2 + v4)) +
        f22*(f4*(-v1 + v3) + f3*(v1 - v4) + f1*(-v3 + v4)))/(f1mf2*f1mf3*f1mf4*f2mf3*f2mf4*f3mf4);
        break;
      }
      case 1043:  //no left derivative: v1, v2, v3, v4, d4
      {
        double f12 = f1*f1;
        double f13 = f12*f1;
        double f14 = f13*f1;

        double f22 = f2*f2;
        double f23 = f22*f2;
        double f24 = f23*f2;

        double f32 = f3*f3;
        double f34 = f32*f32;

        double f42 = f4*f4;
        double f43 = f42*f4;
        double f44 = f43*f4;
        double f45 = f44*f4;

        double f1mf2 = f1-f2;
        double f1mf3 = f1-f3;
        double f1mf4 = f1-f4;
        double f2mf3 = f2-f3;
        double f2mf4 = f2-f4;
        double f3mf4 = f3-f4;

        double f1mf42 = f1mf4*f1mf4;
        double f2mf42 = f2mf4*f2mf4;
        double f3mf42 = f3mf4*f3mf4;

        double v1mv3 = v1-v3;
        double v1mv4 = v1-v4;
        double v3mv4 = v3-v4;

        retVal = (d4*f1mf2*f1mf3*f1mf4*f2mf3*f2mf4*f3mf4*(f1 + f2 + f3 + f4) + f34*f42*v1 - 3*f32*f44*v1 + 2*f3*f45*v1 + f14*f32*v2 - f12*f34*v2 - 2*f14*f3*f4*v2 + 2*f1*f34*f4*v2 + f14*f42*v2 - f34*f42*v2 +
        4*f12*f3*f43*v2 - 4*f1*f32*f43*v2 - 3*f12*f44*v2 + 3*f32*f44*v2 + 2*f1*f45*v2 - 2*f3*f45*v2 - f14*f42*v3 + 3*f12*f44*v3 - 2*f1*f45*v3 +
        f24*(-(f42*v1mv3) - f32*v1mv4 + 2*f3*f4*v1mv4 + f12*v3mv4 - 2*f1*f4*v3mv4) - 2*f2*f4*(f44*v1mv3 + f34*v1mv4 - 2*f32*f42*v1mv4 - f14*v3mv4 + 2*f12*f42*v3mv4) +
        f22*(3*f44*v1mv3 + f34*v1mv4 - 4*f3*f43*v1mv4 - f14*v3mv4 + 4*f1*f43*v3mv4) - f14*f32*v4 + f12*f34*v4 + 2*f14*f3*f4*v4 - 2*f1*f34*f4*v4 - 4*f12*f3*f43*v4 + 4*f1*f32*f43*v4)/
        (f1mf2*f1mf3*f1mf42*f2mf3*f2mf42*f3mf42);
        break;
      }
      case 1042:   //4th order poly: v1,d1, v2,d2, v3  // used for the first intermediate region
      {

        double f12 = f1*f1;
        double f13 = f12*f1;
        double f14 = f13*f1;

        double f42 = f4*f4;
        double f43 = f42*f4;
        double f44 = f43*f4;

        double f32 = f3*f3;
        double f33 = f32*f3;
        double f34 = f33*f3;

        double f1mf4 = f1-f4;
        double f1mf3 = f1-f3;
        double f3mf4 = f3-f4;

        double f1mf42 = f1mf4*f1mf4;
        double f1mf32 = f1mf3*f1mf3;
        double f3mf42 = f3mf4*f3mf4;

        double f1mf43 = f1mf42*f1mf4;

        retVal = (2*d4*f14*f3 - d1*f13*f32 - 3*d4*f13*f32 + d1*f1*f34 + d4*f1*f34 - 2*d4*f14*f4 + 2*d1*f13*f3*f4 + 2*d4*f13*f3*f4 - d1*f12*f32*f4 + d4*f12*f32*f4 - d1*f34*f4 - d4*f34*f4 - d1*f13*f42 + d4*f13*f42 +
          2*d1*f12*f3*f42 - 2*d4*f12*f3*f42 - d1*f1*f32*f42 + d4*f1*f32*f42 - d1*f12*f43 + d4*f12*f43 - 2*d1*f1*f3*f43 - 2*d4*f1*f3*f43 + 3*d1*f32*f43 + d4*f32*f43 + 2*d1*f1*f44 - 2*d1*f3*f44 +
          4*f12*f32*v1 - 2*f34*v1 - 8*f12*f3*f4*v1 + 4*f1*f32*f4*v1 + 4*f12*f42*v1 - 8*f1*f3*f42*v1 + 4*f32*f42*v1 + 4*f1*f43*v1 - 2*f44*v1 - 2*f14*v3 + 4*f13*f4*v3 - 4*f1*f43*v3 + 2*f44*v3 + 2*f14*v4 -
          4*f12*f32*v4 + 2*f34*v4 - 4*f13*f4*v4 + 8*f12*f3*f4*v4 - 4*f1*f32*f4*v4 - 4*f12*f42*v4 + 8*f1*f3*f42*v4 - 4*f32*f42*v4)/(f1mf32*f1mf43*f3mf42);
          break;
        }
        case 104:  //Geraint's Version, 4th order poly: v1,d1, v2,d2, v3
        {

          double f12 = f1*f1;
          double f13 = f12*f1;
          double f14 = f13*f1;

          double f22 = f2*f2;
          double f23 = f22*f2;
          double f24 = f23*f2;

          double f42 = f4*f4;
          double f43 = f42*f4;
          double f44 = f43*f4;

          double f1mf2 = f1-f2;
          double f1mf4 = f1-f4;
          double f2mf4 = f2-f4;

          double f1mf22 = f1mf2*f1mf2;
          double f2mf42 = f2mf4*f2mf4;
          double f1mf43 = f1mf4*f1mf4*f1mf4;

          retVal = ((d4*f1mf22*f1mf4*f2mf4*(2*f1 + f2 + f4) - d1*f1mf2*f1mf4*f2mf42*(f1 + f2 + 2*f4)
          + 2*(f44*(-v1 + v2) + 2*f12*f2mf42*(v1 - v4) + 2*f22*f42*(v1 - v4)
          + 2*f13*f4*(v2 - v4) + f24*(-v1 + v4) + f14*(-v2 + v4) + 2*f1*f4*(f42*(v1 - v2) + f22*(v1 - v4) + 2*f2*f4*(-v1 + v4)))) / (f1mf22*f1mf43*f2mf42));


          break;
        }
        case 105: // Geraint, standard way: v1, v2, v3, v4, d1, d4
        {
          double f12 = f1*f1;
          double f13 = f12*f1;
          double f14 = f13*f1;
          double f15 = f14*f1;
          double f16 = f15*f1;
          double f17 = f16*f1;

          double f22 = f2*f2;
          double f23 = f22*f2;
          double f24 = f23*f2;
          double f25 = f24*f2;

          double f32 = f3*f3;
          double f33 = f32*f3;
          double f34 = f33*f3;
          double f35 = f34*f3;

          double f42 = f4*f4;
          double f43 = f42*f4;
          double f44 = f43*f4;
          double f45 = f44*f4;
          double f46 = f45*f4;
          double f47 = f46*f4;

          double f1mf2 = f1-f2;
          double f1mf3 = f1-f3;
          double f1mf4 = f1-f4;
          double f2mf3 = f2-f3;
          double f2mf4 = f2-f4;
          double f3mf4 = f3-f4;

          double f1mf22 = f1mf2*f1mf2;
          double f1mf32 = f1mf3*f1mf3;
          double f2mf42 = f2mf4*f2mf4;
          double f3mf42 = f3mf4*f3mf4;
          double f1mf43 = f1mf4*f1mf4*f1mf4;

          retVal = (
            (d4*f1mf22*f1mf32*f1mf4*f2mf3*f2mf4*f3mf4*(f12 + f2*f3 + (f2 + f3)*f4 + 2*f1*(f2 + f3 + f4)) + d1*f1mf2*f1mf3*f1mf4*f2mf3*f2mf42*f3mf42*(f1*f2 + f1*f3 + f2*f3 + 2*(f1 + f2 + f3)*f4 + f42) -
            5*f13*f24*f32*v1 + 4*f12*f25*f32*v1 + 5*f13*f22*f34*v1 - 2*f25*f34*v1 - 4*f12*f22*f35*v1 + 2*f24*f35*v1 + 10*f13*f24*f3*f4*v1 - 8*f12*f25*f3*f4*v1 - 5*f12*f24*f32*f4*v1 + 4*f1*f25*f32*f4*v1 -
            10*f13*f2*f34*f4*v1 + 5*f12*f22*f34*f4*v1 + 8*f12*f2*f35*f4*v1 - 4*f1*f22*f35*f4*v1 - 5*f13*f24*f42*v1 + 4*f12*f25*f42*v1 + 10*f12*f24*f3*f42*v1 - 8*f1*f25*f3*f42*v1 - 5*f1*f24*f32*f42*v1 +
            4*f25*f32*f42*v1 + 5*f13*f34*f42*v1 - 10*f12*f2*f34*f42*v1 + 5*f1*f22*f34*f42*v1 - 4*f12*f35*f42*v1 + 8*f1*f2*f35*f42*v1 - 4*f22*f35*f42*v1 - 5*f12*f24*f43*v1 + 4*f1*f25*f43*v1 -
            20*f13*f22*f3*f43*v1 + 10*f1*f24*f3*f43*v1 + 20*f13*f2*f32*f43*v1 - 5*f24*f32*f43*v1 + 5*f12*f34*f43*v1 - 10*f1*f2*f34*f43*v1 + 5*f22*f34*f43*v1 - 4*f1*f35*f43*v1 + 15*f13*f22*f44*v1 -
            5*f1*f24*f44*v1 - 2*f25*f44*v1 - 15*f13*f32*f44*v1 + 5*f1*f34*f44*v1 + 2*f35*f44*v1 - 10*f13*f2*f45*v1 - f12*f22*f45*v1 + 3*f24*f45*v1 + 10*f13*f3*f45*v1 + f12*f32*f45*v1 - 3*f34*f45*v1 +
            2*f12*f2*f46*v1 - f1*f22*f46*v1 - 2*f12*f3*f46*v1 + f1*f32*f46*v1 + 2*f1*f2*f47*v1 - f22*f47*v1 - 2*f1*f3*f47*v1 + f32*f47*v1 + f17*f32*v2 - 3*f15*f34*v2 + 2*f14*f35*v2 - 2*f17*f3*f4*v2 +
            f16*f32*f4*v2 + 5*f14*f34*f4*v2 - 4*f13*f35*f4*v2 + f17*f42*v2 - 2*f16*f3*f42*v2 + f15*f32*f42*v2 + f16*f43*v2 + 10*f15*f3*f43*v2 - 15*f14*f32*f43*v2 + 4*f1*f35*f43*v2 - 8*f15*f44*v2 +
            15*f13*f32*f44*v2 - 5*f1*f34*f44*v2 - 2*f35*f44*v2 + 8*f14*f45*v2 - 10*f13*f3*f45*v2 - f12*f32*f45*v2 + 3*f34*f45*v2 - f13*f46*v2 + 2*f12*f3*f46*v2 - f1*f32*f46*v2 - f12*f47*v2 +
            2*f1*f3*f47*v2 - f32*f47*v2 - f17*f22*v3 + 3*f15*f24*v3 - 2*f14*f25*v3 + 2*f17*f2*f4*v3 - f16*f22*f4*v3 - 5*f14*f24*f4*v3 + 4*f13*f25*f4*v3 - f17*f42*v3 + 2*f16*f2*f42*v3 - f15*f22*f42*v3 -
            f16*f43*v3 - 10*f15*f2*f43*v3 + 15*f14*f22*f43*v3 - 4*f1*f25*f43*v3 + 8*f15*f44*v3 - 15*f13*f22*f44*v3 + 5*f1*f24*f44*v3 + 2*f25*f44*v3 - 8*f14*f45*v3 + 10*f13*f2*f45*v3 + f12*f22*f45*v3 -
            3*f24*f45*v3 + f13*f46*v3 - 2*f12*f2*f46*v3 + f1*f22*f46*v3 + f12*f47*v3 - 2*f1*f2*f47*v3 + f22*f47*v3 +
            f1mf22*f1mf32*f2mf3*(2*f22*f32 + f13*(f2 + f3 - 2*f4) + f12*(f2 + f3 - 2*f4)*(2*(f2 + f3) + f4) - 4*(f22 + f2*f3 + f32)*f42 + 5*(f2 + f3)*f43 +
            f1*(4*f2*f3*(f2 + f3) - 4*(f22 + f2*f3 + f32)*f4 - 3*(f2 + f3)*f42 + 10*f43))*v4)/(f1mf22*f1mf32*f1mf43*f2mf3*f2mf42*f3mf42)
          );

          break;
        }
        default:
        {
          XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomXHM_Intermediate_Amp_delta3: IMRPhenomXIntermediateAmpVersion is not valid.\n");
        }
      }

      return retVal;
}

static double IMRPhenomXHM_Intermediate_Amp_delta4(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag)
{
      double retVal;

      switch ( IntAmpFlag )
      {
        case 101: //linear, only v1, v2
        {
          retVal = 0.;
          break;
        }
        case 102: //quadratic: v1, v2, d2
        {
          retVal = 0.;
          break;
        }
        case 1032:  // 2 freqs, points and derivatives: v1, v4, d1, d4
        {
          retVal = 0.;
          break;
        }
        case 103:   // 4 freqs, no boundaries derivatives
        {
          retVal = 0.;
          break;
        }
        case 1043:  //no left derivative: v1, v2, v3, v4, d4
        {
          double f12 = f1*f1;
          double f13 = f12*f1;

          double f22 = f2*f2;
          double f23 = f22*f2;

          double f32 = f3*f3;
          double f33 = f32*f3;

          double f42 = f4*f4;
          double f43 = f42*f4;
          double f44 = f43*f4;

          double f1mf2 = f1-f2;
          double f1mf3 = f1-f3;
          double f1mf4 = f1-f4;
          double f2mf3 = f2-f3;
          double f2mf4 = f2-f4;
          double f3mf4 = f3-f4;

          double f1mf42 = f1mf4*f1mf4;
          double f2mf42 = f2mf4*f2mf4;
          double f3mf42 = f3mf4*f3mf4;

          double v1mv3 = v1-v3;
          double v1mv4 = v1-v4;
          double v3mv4 = v3-v4;

          retVal = (-(d4*f1mf2*f1mf3*f1mf4*f2mf3*f2mf4*f3mf4) - f33*f42*v1 + 2*f32*f43*v1 - f3*f44*v1 - f13*f32*v2 + f12*f33*v2 + 2*f13*f3*f4*v2 - 2*f1*f33*f4*v2 - f13*f42*v2 - 3*f12*f3*f42*v2 + 3*f1*f32*f42*v2 +
          f33*f42*v2 + 2*f12*f43*v2 - 2*f32*f43*v2 - f1*f44*v2 + f3*f44*v2 + f13*f42*v3 - 2*f12*f43*v3 + f1*f44*v3 + f23*(f42*v1mv3 + f32*v1mv4 - 2*f3*f4*v1mv4 - f12*v3mv4 + 2*f1*f4*v3mv4) +
          f2*f4*(f43*v1mv3 + 2*f33*v1mv4 - 3*f32*f4*v1mv4 - 2*f13*v3mv4 + 3*f12*f4*v3mv4) + f22*(-2*f43*v1mv3 - f33*v1mv4 + 3*f3*f42*v1mv4 + f13*v3mv4 - 3*f1*f42*v3mv4) + f13*f32*v4 - f12*f33*v4 -
          2*f13*f3*f4*v4 + 2*f1*f33*f4*v4 + 3*f12*f3*f42*v4 - 3*f1*f32*f42*v4)/(f1mf2*f1mf3*f1mf42*f2mf3*f2mf42*f3mf42);
          break;
        }
        case 1042:   //4th order poly: v1,d1, v2,d2, v3   // used for the first intermediate region
        {
          double f12 = f1*f1;
          double f13 = f12*f1;

          double f42 = f4*f4;
          double f43 = f42*f4;

          double f32 = f3*f3;
          double f33 = f32*f3;

          double f1mf4 = f1-f4;
          double f1mf3 = f1-f3;
          double f3mf4 = f3-f4;

          double f1mf42 = f1mf4*f1mf4;
          double f1mf32 = f1mf3*f1mf3;
          double f3mf42 = f3mf4*f3mf4;

          double f1mf43 = f1mf42*f1mf4;

          retVal = (-(d4*f1mf32*f1mf4*f3mf4) + d1*f1mf3*f1mf4*f3mf42 - 3*f1*f32*v1 + 2*f33*v1 + 6*f1*f3*f4*v1 - 3*f32*f4*v1 - 3*f1*f42*v1 + f43*v1 + f13*v3 - 3*f12*f4*v3 + 3*f1*f42*v3 - f43*v3 -
          f1mf32*(f1 + 2*f3 - 3*f4)*v4)/(f1mf32*f1mf43*f3mf42);
          break;
        }
        case 104:  //Geraint's Version, 4th order poly: v1,d1, v2, v4,d4
        {
          double f12 = f1*f1;
          double f13 = f12*f1;

          double f22 = f2*f2;
          double f23 = f22*f2;

          double f42 = f4*f4;
          double f43 = f42*f4;

          double f1mf2 = f1-f2;
          double f1mf4 = f1-f4;
          double f2mf4 = f2-f4;

          double f1mf22 = f1mf2*f1mf2;
          double f2mf42 = f2mf4*f2mf4;
          double f1mf43 = f1mf4*f1mf4*f1mf4;

          retVal = ((-(d4*f1mf22*f1mf4*f2mf4) + d1*f1mf2*f1mf4*f2mf42 - 3*f1*f22*v1 + 2*f23*v1 + 6*f1*f2*f4*v1 - 3*f22*f4*v1
          - 3*f1*f42*v1 + f43*v1 + f13*v2 - 3*f12*f4*v2 + 3*f1*f42*v2 - f43*v2 - f1mf22*(f1 + 2*f2 - 3*f4)*v4)/(f1mf22*f1mf43*f2mf42));

          break;
        }
        case 105: // Geraint, standard way: v1, v2, v3, v4, d1, d4
        {
          double f12 = f1*f1;
          double f13 = f12*f1;
          double f14 = f13*f1;
          double f15 = f14*f1;
          double f16 = f15*f1;

          double f22 = f2*f2;
          double f23 = f22*f2;
          double f24 = f23*f2;
          double f25 = f24*f2;

          double f32 = f3*f3;
          double f33 = f32*f3;
          double f34 = f33*f3;
          double f35 = f34*f3;

          double f42 = f4*f4;
          double f43 = f42*f4;
          double f44 = f43*f4;
          double f45 = f44*f4;
          double f46 = f45*f4;

          double f1mf2 = f1-f2;
          double f1mf3 = f1-f3;
          double f1mf4 = f1-f4;
          double f2mf3 = f2-f3;
          double f2mf4 = f2-f4;
          double f3mf4 = f3-f4;

          double f1mf22 = f1mf2*f1mf2;
          double f1mf32 = f1mf3*f1mf3;
          double f2mf42 = f2mf4*f2mf4;
          double f3mf42 = f3mf4*f3mf4;
          double f1mf43 = f1mf4*f1mf4*f1mf4;

          retVal = (
            (-(d4*f1mf22*f1mf32*f1mf4*f2mf3*f2mf4*f3mf4*(2*f1 + f2 + f3 + f4)) - d1*f1mf2*f1mf3*f1mf4*f2mf3*f2mf42*f3mf42*(f1 + f2 + f3 + 2*f4) + 5*f13*f23*f32*v1 - 3*f1*f25*f32*v1 - 5*f13*f22*f33*v1 +
            2*f25*f33*v1 + 3*f1*f22*f35*v1 - 2*f23*f35*v1 - 10*f13*f23*f3*f4*v1 + 6*f1*f25*f3*f4*v1 + 5*f12*f23*f32*f4*v1 - 3*f25*f32*f4*v1 + 10*f13*f2*f33*f4*v1 - 5*f12*f22*f33*f4*v1 -
            6*f1*f2*f35*f4*v1 + 3*f22*f35*f4*v1 + 5*f13*f23*f42*v1 - 3*f1*f25*f42*v1 + 15*f13*f22*f3*f42*v1 - 10*f12*f23*f3*f42*v1 - 15*f13*f2*f32*f42*v1 + 5*f1*f23*f32*f42*v1 - 5*f13*f33*f42*v1 +
            10*f12*f2*f33*f42*v1 - 5*f1*f22*f33*f42*v1 + 3*f1*f35*f42*v1 - 10*f13*f22*f43*v1 + 5*f12*f23*f43*v1 + f25*f43*v1 + 15*f12*f22*f3*f43*v1 - 10*f1*f23*f3*f43*v1 + 10*f13*f32*f43*v1 -
            15*f12*f2*f32*f43*v1 + 5*f23*f32*f43*v1 - 5*f12*f33*f43*v1 + 10*f1*f2*f33*f43*v1 - 5*f22*f33*f43*v1 - f35*f43*v1 + 5*f13*f2*f44*v1 - 10*f12*f22*f44*v1 + 5*f1*f23*f44*v1 - 5*f13*f3*f44*v1 +
            10*f12*f32*f44*v1 - 5*f1*f33*f44*v1 + 5*f12*f2*f45*v1 + 2*f1*f22*f45*v1 - 3*f23*f45*v1 - 5*f12*f3*f45*v1 - 2*f1*f32*f45*v1 + 3*f33*f45*v1 - 4*f1*f2*f46*v1 + 2*f22*f46*v1 + 4*f1*f3*f46*v1 -
            2*f32*f46*v1 - 2*f16*f32*v2 + 3*f15*f33*v2 - f13*f35*v2 + 4*f16*f3*f4*v2 - 2*f15*f32*f4*v2 - 5*f14*f33*f4*v2 + 3*f12*f35*f4*v2 - 2*f16*f42*v2 - 5*f15*f3*f42*v2 + 10*f14*f32*f42*v2 -
            3*f1*f35*f42*v2 + 4*f15*f43*v2 - 5*f14*f3*f43*v2 + f35*f43*v2 + 5*f13*f3*f44*v2 - 10*f12*f32*f44*v2 + 5*f1*f33*f44*v2 - 4*f13*f45*v2 + 5*f12*f3*f45*v2 + 2*f1*f32*f45*v2 - 3*f33*f45*v2 +
            2*f12*f46*v2 - 4*f1*f3*f46*v2 + 2*f32*f46*v2 + 2*f16*f22*v3 - 3*f15*f23*v3 + f13*f25*v3 - 4*f16*f2*f4*v3 + 2*f15*f22*f4*v3 + 5*f14*f23*f4*v3 - 3*f12*f25*f4*v3 + 2*f16*f42*v3 +
            5*f15*f2*f42*v3 - 10*f14*f22*f42*v3 + 3*f1*f25*f42*v3 - 4*f15*f43*v3 + 5*f14*f2*f43*v3 - f25*f43*v3 - 5*f13*f2*f44*v3 + 10*f12*f22*f44*v3 - 5*f1*f23*f44*v3 + 4*f13*f45*v3 - 5*f12*f2*f45*v3 -
            2*f1*f22*f45*v3 + 3*f23*f45*v3 - 2*f12*f46*v3 + 4*f1*f2*f46*v3 - 2*f22*f46*v3 -
            f1mf22*f1mf32*f2mf3*(2*f2*f3*(f2 + f3) + 2*f12*(f2 + f3 - 2*f4) - 3*(f22 + f2*f3 + f32)*f4 + f1*(f22 + 5*f2*f3 + f32 - 6*(f2 + f3)*f4 + 5*f42) + 5*f43)*v4)/
            (f1mf22*f1mf32*f1mf43*f2mf3*f2mf42*f3mf42)
          );

          break;
        }
        default:
        {
          XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomXHM_Intermediate_Amp_delta4: IMRPhenomXIntermediateAmpVersion is not valid.\n");
        }
      }

      return retVal;
}

static double IMRPhenomXHM_Intermediate_Amp_delta5(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag)
{
      double retVal;

      switch ( IntAmpFlag )
      {
        case 101:
        {
          retVal = 0.;
          break;
        }
        case 102: //quadratic: v1, v2, d2
        {
          retVal = 0.;
          break;
        }
        case 1032:  // 2 freqs, points and derivatives: v1, v4, d1, d4
        {
          //printf("\nIMRPhenomXHM_Intermediate_Amp_delta5 = 0 \r\n");
          retVal = 0.;
          break;
        }
        case 103:   // 4 freqs, no boundaries derivatives
        {
          retVal = 0.;
          break;
        }
        case 1043:  //no left derivative: v1, v2, v3, v4, d4
        {
          retVal = 0.;
          break;
        }
        case 1042:   //4th order poly: v1,d1, v2,d2, v3   // used for the first intermediate region
        {
          retVal = 0.0;
          break;
        }
        case 104:  //Geraint's Version, 4th order poly: v1,d1, v2,d2, v3
        {
          retVal = 0.0;
          break;
        }
        case 105: // Geraint, standard way: v1, v2, v3, v4, d1, d4
        {
          double f12 = f1*f1;
          double f13 = f12*f1;
          double f14 = f13*f1;
          double f15 = f14*f1;

          double f22 = f2*f2;
          double f23 = f22*f2;
          double f24 = f23*f2;

          double f32 = f3*f3;
          double f33 = f32*f3;
          double f34 = f33*f3;

          double f42 = f4*f4;
          double f43 = f42*f4;
          double f44 = f43*f4;
          double f45 = f44*f4;

          double f1mf2 = f1-f2;
          double f1mf3 = f1-f3;
          double f1mf4 = f1-f4;
          double f2mf3 = f2-f3;
          double f2mf4 = f2-f4;
          double f3mf4 = f3-f4;

          double f1mf22 = f1mf2*f1mf2;
          double f1mf32 = f1mf3*f1mf3;
          double f2mf42 = f2mf4*f2mf4;
          double f3mf42 = f3mf4*f3mf4;
          double f1mf43 = f1mf4*f1mf4*f1mf4;

          retVal = (
            (d4*f1mf22*f1mf32*f1mf4*f2mf3*f2mf4*f3mf4 + d1*f1mf2*f1mf3*f1mf4*f2mf3*f2mf42*f3mf42 - 4*f12*f23*f32*v1 + 3*f1*f24*f32*v1 + 4*f12*f22*f33*v1 - 2*f24*f33*v1 - 3*f1*f22*f34*v1 + 2*f23*f34*v1 +
              8*f12*f23*f3*f4*v1 - 6*f1*f24*f3*f4*v1 - 4*f1*f23*f32*f4*v1 + 3*f24*f32*f4*v1 - 8*f12*f2*f33*f4*v1 + 4*f1*f22*f33*f4*v1 + 6*f1*f2*f34*f4*v1 - 3*f22*f34*f4*v1 - 4*f12*f23*f42*v1 +
              3*f1*f24*f42*v1 - 12*f12*f22*f3*f42*v1 + 8*f1*f23*f3*f42*v1 + 12*f12*f2*f32*f42*v1 - 4*f23*f32*f42*v1 + 4*f12*f33*f42*v1 - 8*f1*f2*f33*f42*v1 + 4*f22*f33*f42*v1 - 3*f1*f34*f42*v1 +
              8*f12*f22*f43*v1 - 4*f1*f23*f43*v1 - f24*f43*v1 - 8*f12*f32*f43*v1 + 4*f1*f33*f43*v1 + f34*f43*v1 - 4*f12*f2*f44*v1 - f1*f22*f44*v1 + 2*f23*f44*v1 + 4*f12*f3*f44*v1 + f1*f32*f44*v1 -
              2*f33*f44*v1 + 2*f1*f2*f45*v1 - f22*f45*v1 - 2*f1*f3*f45*v1 + f32*f45*v1 + f15*f32*v2 - 2*f14*f33*v2 + f13*f34*v2 - 2*f15*f3*f4*v2 + f14*f32*f4*v2 + 4*f13*f33*f4*v2 - 3*f12*f34*f4*v2 +
              f15*f42*v2 + 4*f14*f3*f42*v2 - 8*f13*f32*f42*v2 + 3*f1*f34*f42*v2 - 3*f14*f43*v2 + 8*f12*f32*f43*v2 - 4*f1*f33*f43*v2 - f34*f43*v2 + 3*f13*f44*v2 - 4*f12*f3*f44*v2 - f1*f32*f44*v2 +
              2*f33*f44*v2 - f12*f45*v2 + 2*f1*f3*f45*v2 - f32*f45*v2 - f15*f22*v3 + 2*f14*f23*v3 - f13*f24*v3 + 2*f15*f2*f4*v3 - f14*f22*f4*v3 - 4*f13*f23*f4*v3 + 3*f12*f24*f4*v3 - f15*f42*v3 -
              4*f14*f2*f42*v3 + 8*f13*f22*f42*v3 - 3*f1*f24*f42*v3 + 3*f14*f43*v3 - 8*f12*f22*f43*v3 + 4*f1*f23*f43*v3 + f24*f43*v3 - 3*f13*f44*v3 + 4*f12*f2*f44*v3 + f1*f22*f44*v3 - 2*f23*f44*v3 +
              f12*f45*v3 - 2*f1*f2*f45*v3 + f22*f45*v3 + f1mf22*f1mf32*f2mf3*(2*f2*f3 + f1*(f2 + f3 - 2*f4) - 3*(f2 + f3)*f4 + 4*f42)*v4)/(f1mf22*f1mf32*f1mf43*f2mf3*f2mf42*f3mf42)
            );
            break ;
          }
          default:
          {
            XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomXHM_Intermediate_Amp_delta5: IMRPhenomXIntermediateAmpVersion is not valid.\n");
          }
        }
        return retVal;
}


/*********************************************/
/*         INTERMEDIATE AMPLITUDE ANSATZ     */
/*********************************************/

// Build the polynomial with the coefficients given and return the inverse of the polynomial (this is the ansatz)
static double IMRPhenomXHM_Intermediate_Amp_Ansatz(IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMAmpCoefficients *pAmp)
{
    if(pWFHM->IMRPhenomXHMReleaseVersion != 122019){
        double result = 0., fpower = 1.;
        /* Ansatz = f^(-7/6) * polynomial */
        for (UINT2 i = 0; i < pAmp->nCoefficientsInter; i++){
            //printf("pAmp->nCoefficientsInter = %i, pAmp->InterCoefficient[%i] = %.16e\n", pAmp->nCoefficientsInter, i, pAmp->InterCoefficient[i]);
            result += (pAmp->InterCoefficient[i] * fpower);
            fpower *= powers_of_f->itself;
        }
        result *= powers_of_f->m_seven_sixths;
        return result;
    }
    else{
        double a0 = pAmp->delta0;
        double a1 = pAmp->delta1;
        double a2 = pAmp->delta2;
        double a3 = pAmp->delta3;
        double a4 = pAmp->delta4;
        double a5 = pAmp->delta5;
        double polynomial;
        int InterAmpPolOrder = pAmp->InterAmpPolOrder;

        switch ( InterAmpPolOrder )
        {
          case 101:   // linear order
          {
            double ff = powers_of_f->itself;
            polynomial = a0 + a1*ff ;
            break;
          }
          case 102:   // quadratic order
          {
            double ff = powers_of_f->itself;
            double ff2 = powers_of_f->two;
            polynomial = a0 + a1*ff + a2*ff2;
            break;
          }
          case 103:   //cubic order
          {
            double ff  = powers_of_f->itself;
            double ff2 = powers_of_f->two;
            double ff3 = powers_of_f->three;
            polynomial = a0 + a1*ff + a2*ff2 + a3*ff3;
            break;
          }
          case 1042:    // 4th order, used for the first intermediate region
          {
            double ff  = powers_of_f->itself;
            double ff2 = powers_of_f->two;
            double ff3 = powers_of_f->three;
            double ff4 = powers_of_f->four;
            a0         = pAmp->alpha0;
            a1         = pAmp->alpha1;
            a2         = pAmp->alpha2;
            a3         = pAmp->alpha3;
            a4         = pAmp->alpha4;
            polynomial = a0 + a1*ff + a2*ff2 + a3*ff3 + a4*ff4;
            break;
          }
          case 104:     // 4th order
          {
            double ff = powers_of_f->itself;
            double ff2 = powers_of_f->two;
            double ff3 = powers_of_f->three;
            double ff4 = powers_of_f->four;
            polynomial = a0 + a1*ff + a2*ff2 + a3*ff3 + a4*ff4;
            break;
          }
          case 105:     // 5th order
          {
            double ff  = powers_of_f->itself;
            double ff2 = powers_of_f->two;
            double ff3 = powers_of_f->three;
            double ff4 = powers_of_f->four;
            double ff5 = powers_of_f->five;
            polynomial = a0 + a1*ff + a2*ff2 + a3*ff3 + a4*ff4 + a5*ff5;
            break;
          }
          default:
          {
            XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomXHM_Intermediate_Amp_Ansatz: InterAmpPolOrder is not valid.\n");
          }
        }
        return pAmp->ampNorm / polynomial;
    }
}


/**** VETO functions ****/

// Utility functions used to decide how many collocation points and which kind of reconstruction we do in the intermediate


// Remove too low collocation points (heuristic)
void IMRPhenomXHM_Intermediate_Amplitude_Veto(double *int1, double *int2, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22){
  double threshold = 0.2/(pWF22->ampNorm);
  if( 1./(*int1) < threshold )
  {
    *int1 = 1.;
    pWFHM->IMRPhenomXHMIntermediateAmpVersion = 1042;
    if( 1./(*int2) < threshold ){
      *int2 = 1.;
      pWFHM->IMRPhenomXHMIntermediateAmpVersion = 1032;
    }
  }else if( 1./(*int2) < threshold )
  {
    *int2 = 1.;
    pWFHM->IMRPhenomXHMIntermediateAmpVersion = 1042;
  }
}

// Check if a particular frequency belong to an interval
int InsideInterval(double ftest, double fmin, double fmax){

  if(ftest >= fmin && ftest <= fmax){
    return 1;
  }else{
    return 0;
  }
}

// Check if the 2nd order polynomial crosses zero
int CrossZeroP2(double a0, double a1, double a2, double fstart, double fend){

  double complex f1, f2;
  double discriminant = -4*a0*a2 + a1*a1;

  f1 = creall((-a1 - sqrt(discriminant))/(2.*a2));
  f2 = creall((-a1 + sqrt(discriminant))/(2.*a2));

  if( discriminant >= 0 && (InsideInterval(f1, fstart, fend) || InsideInterval(f2, fstart, fend) )){
    return 1;
  }else{
    return 0;
  }
}

// Check if the 3th order polynomial crosses zero
int CrossZeroP3(double a0, double a1, double a2, double a3, double fstart, double fend)
{

  double q1 = -a2/(3.*a3);
  double q2 = -a2*a2 + 3*a1*a3;
  double q3 = -2*a2*a2*a2 + 9*a1*a2*a3-27*a0*a3*a3;
  double discri = 4*q2*q2*q2 + q3*q3;
  double complex onethirdpower = cpow( q3 + csqrt(discri) , 1/3.);
  double complex f1, f2, f3;
  double twothird = pow(2,1/3.);
  double complex z1 = (1. + I*sqrt(3.)), z2 = (1. - I*sqrt(3.));
  int i1=0, i2=0, i3=0;

  // These are the points where the polynomial is zero
  f1 = q1 - q2*twothird/(3*a3*onethirdpower) + onethirdpower/(3.*a3*twothird);
  f2 = q1 + q2*z1/(3*twothird*twothird*a3*onethirdpower) - onethirdpower*z2/(6*twothird*a3);
  f3 = q1 + q2*z2/(3*twothird*twothird*a3*onethirdpower) - onethirdpower*z1/(6*twothird*a3);

  // If the solution is real and lay inside the interval then we return true.
  // Sometimes due to finite precission the imaginary part can be not exactly zero and we set the threshold 10^-15 as limiting value.
  double threshold = pow(10.,-15);
  if (fabs(cimag(f1)) < threshold ){
    i1 = InsideInterval(creall(f1), fstart, fend);
  }
  if (fabs(cimag(f2)) < threshold ){
    i2 = InsideInterval(creall(f2), fstart, fend);
  }
  if (fabs(cimag(f3)) < threshold ){
    i3 = InsideInterval(creall(f3), fstart, fend);
  }
  #if DEBUG == 1
  printf("\nCrossZeroP3: Coefficients = %.16f %.16f %.16f %.16f",a0, a1, a2, a3);
  printf("\nCrossZeroP3: discri = %.16f ", discri);
  printf("\nCrossZeroP3: Imag f1, f2, f3 = %.16e %.16e %.16e ", fabs(cimag(f1)), fabs(cimag(f2)), fabs(cimag(f3)));
  printf("\nCrossZeroP3: a3 = %.16f %.16f", creal(a3), cimag(a3));
  printf("\nCrossZeroP3: onethirdpower = %.16f %.16f", creal(onethirdpower), cimag(onethirdpower));
  printf("\nCrossZeroP3: f1 = %.16e %.16e", creal(f1), cimag(f1));
  printf("\nCrossZeroP3: f2 = %.16e %.16e", creal(f2), cimag(f2));
  printf("\nCrossZeroP3: f3 = %.16e %.16e", creal(f3), cimag(f3));
  printf("\nCrossZeroP3: i1, i2, i3 = %i %i %i", i1, i2, i3);
  #endif
  if (i1 == 0 && i2 == 0 && i3 == 0 ){
    return 0;
  }else{
    return 1;
  }
}

// Check if the 4th order polynomial crosses zero
int CrossZeroP4(double a0, double a1, double a2, double a3, double a4, double fstart, double fend){

  double q0 = -a3/(4*a4);
  double q1 = a3*a3/(4*a4*a4) -2*a2/(3*a4);
  double q1b = 2*q1;
  double q2 = a2*a2 - 3*a1*a3 + 12*a0*a4;
  double q3 = 2*a2*a2*a2 - 9*a1*a2*a3 + 27*a0*a3*a3 + 27*a1*a1*a4 - 72*a0*a2*a4;
  double complex squareroot = csqrt(-4*q2*q2*q2 + q3*q3);
  double complex onethird = cpow(q3 + squareroot , 1/3.);
  double twothird = pow(2,1/3.);
  double complex frac1 = twothird*q2/(3*a4*onethird);
  double complex frac2 = onethird/(3*twothird*a4);
  double complex bigdenom = 4*sqrt(q1 + frac1 + frac2);
  double bignum = -a3*a3*a3/(a4*a4*a4) + 4*a2*a3/(a4*a4) - 8*a1/a4;
  double threshold = pow(10.,-15);

  complex double f1, f2, f3, f4;

  // These are the solutions
  f1 = q0 - 0.5*csqrt(q1 + frac1 + frac2) - 0.5*csqrt(q1b - frac1 - frac2 - bignum/bigdenom);
  f2 = q0 - 0.5*csqrt(q1 + frac1 + frac2) + 0.5*csqrt(q1b - frac1 - frac2 - bignum/bigdenom);
  f3 = q0 + 0.5*csqrt(q1 + frac1 + frac2) - 0.5*csqrt(q1b - frac1 - frac2 + bignum/bigdenom);
  f4 = q0 + 0.5*csqrt(q1 + frac1 + frac2) + 0.5*csqrt(q1b - frac1 - frac2 + bignum/bigdenom);

  #if DEBUG == 1
  printf("\n***** CrossZeroP4 *********\n");
  printf("q0, q1, q1b, q2, q3 %.16f %.16f %.16f %.16f %.16f\n", q0, q1, q1b, q2, q3);
  printf("bigdenom bignum %.16f %.16f %.16f %.16f\n", creal(bigdenom), cimag(bigdenom), creal(bignum), cimag(bignum));
  printf("squareroot %.16f  %.16f \n", creal(squareroot), cimag(squareroot));
  printf("frac1 frac2 %.16f %.16f %.16f %.16f \n", creal(frac1), cimag(frac1), creal(frac2), cimag(frac2));
  printf("twothird onethird %.16f %.16f %.16f %.16f", creal(twothird), cimag(twothird), creal(onethird), cimag(onethird));
  printf("\nfstart, fend = %.16f %.16f\n", fstart, fend);
  printf("\nf1 = %.16f %.16f", creal(f1), cimag(f1));
  printf("\nf2 = %.16f %.16f", creal(f2), cimag(f2));
  printf("\nf3 = %.16f %.16f", creal(f3), cimag(f3));
  printf("\nf4 = %.16f %.16f", creal(f4), cimag(f4));
  #endif

  int i1=0, i2=0, i3=0, i4=0;

  // Check if the soultions are real and lay in the interval
  if (fabs(cimag(f1)) < threshold ){
    i1 = InsideInterval(creall(f1), fstart, fend);
  }
  if (fabs(cimag(f2)) < threshold){
    i2 = InsideInterval(creall(f2), fstart, fend);
  }
  if (fabs(cimag(f3)) < threshold ){
    i3 = InsideInterval(creall(f3), fstart, fend);
  }
  if (fabs(cimag(f4)) < threshold ){
    i4 = InsideInterval(creall(f4), fstart, fend);
  }

  if (i1 == 0 && i2 == 0 && i3 == 0 && i4 == 0 ){
    return 0;
  }else{
    return 1;
  }
}


// Check if the 5th order polynomial crosses zero.
// In this case there is no analytical solution and we have to check numerically
int CrossZeroP5(double a0, double a1, double a2, double a3, double a4, double a5, double fstart, double fend)
{//https://www.gnu.org/software/gsl/doc/html/poly.html
  double threshold = pow(10.,-15);
  /* coefficients of P(x) =  -1 + x^5  */
  double a[6] = { a0, a1, a2, a3, a4, a5 };
  double z[10];

  gsl_poly_complex_workspace * w
      = gsl_poly_complex_workspace_alloc (6);

  gsl_poly_complex_solve (a, 6, w, z);

  gsl_poly_complex_workspace_free (w);

  #if DEBUG == 1
  for (int i = 0; i < 5; i++)
    {

      printf ("z%d = %+.18f %+.18f\n",
              i, z[2*i], z[2*i+1]);
    }
    #endif

    int i1=0, i2=0, i3=0, i4=0, i5=0;

    // Check if the soultions are real and lay in the interval
    if (fabs(z[1]) < threshold ){
      i1 = InsideInterval(z[0], fstart, fend);
    }
    if (fabs(z[3]) < threshold){
      i2 = InsideInterval(z[2], fstart, fend);
    }
    if (fabs(z[4]) < threshold ){
      i3 = InsideInterval(z[4], fstart, fend);
    }
    if (fabs(z[5]) < threshold ){
      i4 = InsideInterval(z[6], fstart, fend);
    }
    if (fabs(z[7]) < threshold ){
      i5 = InsideInterval(z[8], fstart, fend);
    }

    if (i1 == 0 && i2 == 0 && i3 == 0 && i4 == 0 && i5 == 0 ){
      return 0;
    }else{
      return 1;
    }

  return 0;
}




// Get the coefficients of the polynomial for a particular reconstruction indicated with IntAmpFlag
void Update_Intermediate_Amplitude_Coefficients(IMRPhenomXHMAmpCoefficients *pAmp, int IntAmpFlag)
{
    double d1 = pAmp->d1;
    double d4 = pAmp->d4;
    double v1 = pAmp->v1;
    double v2 = pAmp->v2;
    double v3 = pAmp->v3;
    double v4 = pAmp->v4;
    double f1 = pAmp->f1;
    double f2 = pAmp->f2;
    double f3 = pAmp->f3;
    double f4 = pAmp->f4;
    #if DEBUG == 1
    printf("\n UpdateCoeff IMRPhenomXHMIntermediateAmpVersion = %i", IntAmpFlag);
    #endif
    pAmp->delta0 = IMRPhenomXHM_Intermediate_Amp_delta0(d1,d4,v1,v2,v3,v4,f1,f2,f3,f4,IntAmpFlag);
    pAmp->delta1 = IMRPhenomXHM_Intermediate_Amp_delta1(d1,d4,v1,v2,v3,v4,f1,f2,f3,f4,IntAmpFlag);
    pAmp->delta2 = IMRPhenomXHM_Intermediate_Amp_delta2(d1,d4,v1,v2,v3,v4,f1,f2,f3,f4,IntAmpFlag);
    pAmp->delta3 = IMRPhenomXHM_Intermediate_Amp_delta3(d1,d4,v1,v2,v3,v4,f1,f2,f3,f4,IntAmpFlag);
    pAmp->delta4 = IMRPhenomXHM_Intermediate_Amp_delta4(d1,d4,v1,v2,v3,v4,f1,f2,f3,f4,IntAmpFlag);
    pAmp->delta5 = IMRPhenomXHM_Intermediate_Amp_delta5(d1,d4,v1,v2,v3,v4,f1,f2,f3,f4,IntAmpFlag);
}

// Check if the polynomial crosses zero wrapper.
// If it crosses zero, then remove one collocation point and lower the order of the polynomial.
// The order can be as low as linear order, when it is certain that it does not cross zero.
void ChoosePolOrder(IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMAmpCoefficients *pAmp){
  switch(pWFHM->IMRPhenomXHMIntermediateAmpVersion)
  {
    case 105:  // v1, v2, v3, v4, d1, d4
    {
      #if DEBUG == 1
      printf("\nChoosePolOrder 105\n");
      #endif
      // struct pol5_params params;
      // params.a0 = pAmp->delta0;
      // params.a1 = pAmp->delta1;
      // params.a2 = pAmp->delta2;
      // params.a3 = pAmp->delta3;
      // params.a4 = pAmp->delta4;
      // params.a5 = pAmp->delta5;

      double fstart, fend;
      fstart = pAmp->fAmpMatchIN;
      fend = pAmp->fAmpMatchIM;
      #if DEBUG == 1
      printf("\n In ChoosePolOrder\n");
      #endif
      if(CrossZeroP5(pAmp->delta0, pAmp->delta1, pAmp->delta2, pAmp->delta3, pAmp->delta4, pAmp->delta5, fstart, fend)==1){
      //if (RootPol5_finder_gsl(params, fstart, fend) == 1){
        #if DEBUG == 1
        printf("\n Pol5 crosses zero \n");
        #endif
        //int dummy = RootPol5_finder_gsl(params, fstart, fend);
        //update Coefficients and intermediate version
        Update_Intermediate_Amplitude_Coefficients(pAmp, 1042);
        pWFHM->IMRPhenomXHMIntermediateAmpVersion=104;
        //check if the lower order cross zero
        ChoosePolOrder(pWFHM, pAmp);
      }else{
        pAmp->InterAmpPolOrder = 105;
      }
      break;
    }
    case 1043:  //no left derivative: v1, v2, v3, v4, d4
    {
      #if DEBUG == 1
      printf("\nChoosePolOrder 1043\n");
      #endif
      double a0, a1, a2, a3, a4, fstart, fend;
      a0 = pAmp->delta0;
      a1 = pAmp->delta1;
      a2 = pAmp->delta2;
      a3 = pAmp->delta3;
      a4 = pAmp->delta4;
      fstart = pAmp->fAmpMatchIN;
      fend = pAmp->fAmpMatchIM;

      if(CrossZeroP4(a0, a1, a2, a3, a4, fstart, fend) == 1){
        //update Coefficients and intermediate version
        Update_Intermediate_Amplitude_Coefficients(pAmp, 1032);
        pWFHM->IMRPhenomXHMIntermediateAmpVersion=1032;
        //check if the lower order cross zero
        ChoosePolOrder(pWFHM, pAmp);
      }else{
        pAmp->InterAmpPolOrder = 104;
      }
      break;
    }
    case 1042:  //Remove one inter collocation point: v1, d1, v3, v4, d4.
    {
      #if DEBUG == 1
      printf("\nChoosePolOrder 1042\n");
      #endif
      double a0, a1, a2, a3, a4, fstart, fend;
      a0 = pAmp->delta0;
      a1 = pAmp->delta1;
      a2 = pAmp->delta2;
      a3 = pAmp->delta3;
      a4 = pAmp->delta4;
      fstart = pAmp->fAmpMatchIN;
      fend = pAmp->fAmpMatchIM;

      if(CrossZeroP4(a0, a1, a2, a3, a4, fstart, fend) == 1){
        //update Coefficients and intermediate version
        Update_Intermediate_Amplitude_Coefficients(pAmp, 1032);
        pWFHM->IMRPhenomXHMIntermediateAmpVersion=1032;
        //check if the lower order cross zero
        ChoosePolOrder(pWFHM, pAmp);
      }else{
        pAmp->InterAmpPolOrder = 104;
      }
      break;
    }
    case 104:    //Remove one inter collocation point: v1, d1, v2, v4,d4,
    {
      #if DEBUG == 1
      printf("\nChoosePolOrder 104\n");
      #endif
      double a0, a1, a2, a3, a4, fstart, fend;
      a0 = pAmp->delta0;
      a1 = pAmp->delta1;
      a2 = pAmp->delta2;
      a3 = pAmp->delta3;
      a4 = pAmp->delta4;
      fstart = pAmp->fAmpMatchIN;
      fend = pAmp->fAmpMatchIM;

      if(CrossZeroP4(a0, a1, a2, a3, a4, fstart, fend) == 1){
        //update Coefficients and intermediate version
        Update_Intermediate_Amplitude_Coefficients(pAmp, 1032);
        pWFHM->IMRPhenomXHMIntermediateAmpVersion=1032;
        //check if the lower order cross zero
        ChoosePolOrder(pWFHM, pAmp);
      }else{
        pAmp->InterAmpPolOrder = 104;
      }
      break;
    }
    case 1032:  // 2 freqs, points and derivatives: v1, v4, d1, d4
    {
      #if DEBUG == 1
      printf("\nChoosePolOrder 1032\n");
      #endif
      double a0, a1, a2, a3, fstart, fend;
      a0 = pAmp->delta0;
      a1 = pAmp->delta1;
      a2 = pAmp->delta2;
      a3 = pAmp->delta3;
      fstart = pAmp->fAmpMatchIN;
      fend = pAmp->fAmpMatchIM;

      if(CrossZeroP3(a0, a1, a2, a3, fstart, fend) == 1){
        //update Coefficients and intermediate version
        #if DEBUG == 1
        printf("\nCrossZeroP3 True\n");
        #endif
        Update_Intermediate_Amplitude_Coefficients(pAmp, 102);
        pWFHM->IMRPhenomXHMIntermediateAmpVersion = 102;
        //check if the lower order cross zero
        ChoosePolOrder(pWFHM, pAmp);
      }else{
        #if DEBUG == 1
        printf("\nCrossZeroP3 False\n");
        #endif
        pAmp->InterAmpPolOrder = 103;
      }
      break;
    }
    case 103:   // 4 freqs, no boundaries derivatives
    {
      #if DEBUG == 1
      printf("\nChoosePolOrder 103\n");
      #endif
      double a0, a1, a2, a3, fstart, fend;
      a0 = pAmp->delta0;
      a1 = pAmp->delta1;
      a2 = pAmp->delta2;
      a3 = pAmp->delta3;
      fstart = pAmp->fAmpMatchIN;
      fend = pAmp->fAmpMatchIM;

      if(
        CrossZeroP3(a0, a1, a2, a3, fstart, fend) == 1){
          //update Coefficients and intermediate version
          Update_Intermediate_Amplitude_Coefficients(pAmp, 102);
          pWFHM->IMRPhenomXHMIntermediateAmpVersion=102;
          //check if the lower order cross zero
          ChoosePolOrder(pWFHM, pAmp);
        }else{
          pAmp->InterAmpPolOrder = 103;
        }
        break;
      }
      case 102: //quadratic: v1, v2, d2
      {
        #if DEBUG == 1
        printf("\nChoosePolOrder 102\n");
        #endif
        double a0, a1, a2, fstart, fend;
        a0 = pAmp->delta0;
        a1 = pAmp->delta1;
        a2 = pAmp->delta2;
        fstart = pAmp->fAmpMatchIN;
        fend = pAmp->fAmpMatchIM;

        if(CrossZeroP2(a0, a1, a2, fstart, fend) == 1){
          //update Coefficients and intermediate version
          Update_Intermediate_Amplitude_Coefficients(pAmp, 101);
          pWFHM->IMRPhenomXHMIntermediateAmpVersion=101;
          //check if the lower order cross zero
          ChoosePolOrder(pWFHM, pAmp);
        }else{
          pAmp->InterAmpPolOrder = 102;
        }
        break;
      }
      case 101: //linear, only v1, v2
      {
        #if DEBUG == 1
        printf("\nChoosePolOrder 101\n");
        #endif
        //The linear reconstruction is not going to cross zero, since both boundaries to connect are positive
        pAmp->InterAmpPolOrder = 101;
        break;
      }
    }
}

static void IMRPhenomXHM_Intermediate_Amp_CollocationPoints(IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22, IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXAmpCoefficients *pAmp22, IMRPhenomXPhaseCoefficients *pPhase22){


    /* Define collocation points frequencies */
    switch(pWFHM->IMRPhenomXHMIntermediateAmpFreqsVersion){
        case 0:{ // Equispaced. Get boundaries too
            REAL8 deltaf = (pAmp->fAmpMatchIM - pAmp->fAmpMatchIN) / (pWFHM->nCollocPtsInterAmp - 1);
            UINT2 idx = 0;
            for (UINT2 i = 0; i < pWFHM->nCollocPtsInterAmp; i++){
                if(pAmp->VersionCollocPtsInter[i] == 1){
                    // Add point
                    pAmp->CollocationPointsFreqsAmplitudeInter[idx] = pAmp->fAmpMatchIN + deltaf * i;
                }
                else if (pAmp->VersionCollocPtsInter[i] == 2){
                    // Add point + derivative
                    pAmp->CollocationPointsFreqsAmplitudeInter[idx] = pAmp->fAmpMatchIN + deltaf * i;
                    pAmp->CollocationPointsFreqsAmplitudeInter[idx + 1] = pAmp->CollocationPointsFreqsAmplitudeInter[idx];
                }
                idx += pAmp->VersionCollocPtsInter[i];
            }
            break;
        }
        case 1:{ // Chebyshev. Get boundaries too
            REAL8 semisum = 0.5 * (pAmp->fAmpMatchIN + pAmp->fAmpMatchIM);
            REAL8 semidif = 0.5 * (pAmp->fAmpMatchIM - pAmp->fAmpMatchIN);
            for (INT4 i = pWFHM->nCollocPtsInterAmp + 1; i >=0; i--){
                pAmp->CollocationPointsFreqsAmplitudeInter[i] = semisum + semidif * cos( i * LAL_PI / pWFHM->nCollocPtsInterAmp );
            }
            break;
        }
        default: {XLAL_ERROR_VOID(XLAL_EDOM, "Error in IMRPhenomXHM_Intermediate_Amp_CollocationPoints: IMRPhenomXHMIntermediateAmpFreqsVersion = %i is not valid. \n", pWFHM->IMRPhenomXHMInspiralAmpFreqsVersion);}
    }
    /* Define values */

    IMRPhenomX_UsefulPowers powers_of_finsp;
    IMRPhenomX_Initialize_Powers(&powers_of_finsp, pAmp->fAmpMatchIN);
    // Be careful with TGR thing. Inspiral affecting intermediate region.
    UINT2 tmp_factor = pAmp->InspRescaleFactor;
    pAmp->InspRescaleFactor = pAmp->InterRescaleFactor;
    UINT2 tmpnCollocPts = 0;
    if (pAmp->VersionCollocPtsInter[0] == 1){
        pAmp->CollocationPointsValuesAmplitudeInter[0] = IMRPhenomXHM_Inspiral_Amp_Ansatz(&powers_of_finsp, pWFHM, pAmp);//FIXME: Rescale with InterFactor
        tmpnCollocPts += 1;
    }
    else if (pAmp->VersionCollocPtsInter[0] == 2){
        pAmp->CollocationPointsValuesAmplitudeInter[0] = IMRPhenomXHM_Inspiral_Amp_Ansatz(&powers_of_finsp, pWFHM, pAmp);
        pAmp->CollocationPointsValuesAmplitudeInter[1] = IMRPhenomXHM_Inspiral_Amp_NDAnsatz(&powers_of_finsp, pWFHM, pAmp);//pAmp->CollocationPointsValuesAmplitudeInsp[2];
        tmpnCollocPts += 2;
    }
    pAmp->InspRescaleFactor = tmp_factor;

    /* Call parameter space fits */
    /* pWFHM->nCollocPtsInterAmp also includes the boundaries, for the parameter space fits we just need the in between points */
    UINT2 idx = 0;
    for(UINT2 i = 1; i < pWFHM->nCollocPtsInterAmp - 1; i++){
        if(i <= 2)
            idx = pWFHM->modeInt * 2 + i - 1;
        else
            idx = pWFHM->modeInt * 2 + i - 3 + 16; //FIXME
        if(pAmp->VersionCollocPtsInter[i] == 1){
            pAmp->CollocationPointsValuesAmplitudeInter[tmpnCollocPts] = fabs(pAmp->IntermediateAmpFits[idx](pWF22, pWFHM->IMRPhenomXHMIntermediateAmpFitsVersion)); //FIXME: negative fits?
            tmpnCollocPts++;
        }
        // FIXME: throw error here if pAmp->VersionCollocPtsInter[i] == 2, cannot use derivatives if you don't have the parameter space fit
    }

    tmp_factor = pAmp->RDRescaleFactor;
    pAmp->RDRescaleFactor = pAmp->InterRescaleFactor;
    IMRPhenomX_UsefulPowers powers_of_fRD;
    IMRPhenomX_Initialize_Powers(&powers_of_fRD, pAmp->fAmpMatchIM);
    switch(pAmp->VersionCollocPtsInter[pWFHM->nCollocPtsInterAmp - 1]){
        case 1:{ // Add point
            if (pWFHM->MixingOn == 0){
                pAmp->CollocationPointsValuesAmplitudeInter[tmpnCollocPts] = IMRPhenomXHM_RD_Amp_Ansatz(&powers_of_fRD, pWFHM, pAmp);//pAmp->CollocationPointsValuesAmplitudeRD[0];
            }
            else{
                pAmp->CollocationPointsValuesAmplitudeInter[tmpnCollocPts] = cabs(SpheroidalToSpherical(&powers_of_fRD, pAmp22, pPhase22, pAmp, pPhase, pWFHM, pWF22));
            }
            tmpnCollocPts++;
            break;
        }
        case 2:{ // Add point + derivative
            if (pWFHM->MixingOn == 0){
                pAmp->CollocationPointsValuesAmplitudeInter[tmpnCollocPts] = IMRPhenomXHM_RD_Amp_Ansatz(&powers_of_fRD, pWFHM, pAmp);//pAmp->CollocationPointsValuesAmplitudeRD[0];
                pAmp->CollocationPointsValuesAmplitudeInter[tmpnCollocPts + 1] = IMRPhenomXHM_RD_Amp_DAnsatz(&powers_of_fRD, pWFHM, pAmp);//pAmp->CollocationPointsValuesAmplitudeRD[0];
            }
            else{
                pAmp->CollocationPointsValuesAmplitudeInter[tmpnCollocPts] = cabs(SpheroidalToSpherical(&powers_of_fRD, pAmp22, pPhase22, pAmp, pPhase, pWFHM, pWF22));
                pAmp->CollocationPointsValuesAmplitudeInter[tmpnCollocPts + 1] = IMRPhenomXHM_RD_Amp_NDAnsatz(&powers_of_fRD, pAmp, pPhase, pWFHM, pAmp22, pPhase22, pWF22);
            }
            tmpnCollocPts += 2;
            break;
        }
        default: {XLALPrintError("Error in IMRPhenomXHM_Intermediate_Amp_Coefficients: version %i is not valid.", pAmp->VersionCollocPtsInter[pAmp->nCoefficientsInter - 1]);}
    }
    pAmp->RDRescaleFactor = tmp_factor;

    /* tmpnCollocPts must be the same tahn pWFHM->nCollocPtsInterAmp + 2 (=number of free coefficients in intermediate ansatz) */
    if(tmpnCollocPts != pAmp->nCoefficientsInter)
        XLAL_ERROR_VOID(XLAL_EFUNC, "IMRPhenomXHM_Intermediate_Amp_CollocationPoints failed. Inconsistent number of free parameters %i, %i.", tmpnCollocPts, pAmp->nCoefficientsInter);
}

void IMRPhenomXHM_Intermediate_Amp_Coefficients(IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22, IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXAmpCoefficients *pAmp22, IMRPhenomXPhaseCoefficients *pPhase22){

    /* Previously we checked that nCollocPtsInterAmp read from the IMRPhenomXHMIntermediateAmpVersion is equal to the number of free coefficients in the ansatz. */
    UINT2 nCollocPtsInterAmp = pAmp->nCoefficientsInter;

    /* Define set of collocation points */
    IMRPhenomXHM_Intermediate_Amp_CollocationPoints(pAmp, pWFHM, pWF22, pPhase, pAmp22, pPhase22);

    /* GSL objects for solving system of equations via LU decomposition */
    gsl_vector *b, *x;
    gsl_matrix *A;
    gsl_permutation *p;
    int signum; // No need to set, used internally by gsl_linalg_LU_decomp

    /* Initialize gsl objects */
    p = gsl_permutation_alloc(nCollocPtsInterAmp);
    b = gsl_vector_alloc(nCollocPtsInterAmp);
    x = gsl_vector_alloc(nCollocPtsInterAmp);
    A = gsl_matrix_alloc(nCollocPtsInterAmp, nCollocPtsInterAmp);


    /* Define linear system of equations: A x = b */
    /* x is the solution vector: the coefficients of the intermediate ansatz */
    /* b is the vector of collocation points for a set of frequencies */
    /* A is the matrix of multiplicative factors to each coefficient of the ansatz.
       Each row gives the ansatz evaluated at a collocation point frequency. */
    UINT2 tmpnCollocPts = 0;
    for(UINT2 i = 0; i < pWFHM->nCollocPtsInterAmp; i++){
        /* Skip the 0 cases, means that collocation point is not used */
        if(pAmp->VersionCollocPtsInter[i] > 0){
            // Set b vector
            gsl_vector_set(b, tmpnCollocPts, pAmp->CollocationPointsValuesAmplitudeInter[tmpnCollocPts]);
            //FIXME: distinguish InterAmp ansatzaes versions
            // Set system matrix: Polynomial/f^7/6 at the collocation points frequencies.
            /* A = (1, f1, f1^2, f1^3, f1^4, ...) * f1^(-7/6)
                   (1, f2, f2^2, f2^3, f2^4, ...) * f2^(-7/6)
                   ....
                   Until number of collocation points
              If one of the conditions is a derivative we substitute by
                   (0, 1, 2*fi, 3*fi^2, 4*fi^3, ...)
            */
            REAL8 fcollpoint = pAmp->CollocationPointsFreqsAmplitudeInter[tmpnCollocPts];
            REAL8 seven_sixths = 7/6.;
            REAL8 fcollpoint_m_seven_sixths = pow(fcollpoint, -seven_sixths);
            REAL8 fpower = 1.; // 1, f, f^2, f^3, f^4, ...

            // Add equation for Point. Assuming I will always use point, not the derivative alone
            for(INT4 j = 0; j < nCollocPtsInterAmp; j++){
              gsl_matrix_set(A, tmpnCollocPts, j, fpower * fcollpoint_m_seven_sixths);
              fpower *= fcollpoint;
            }
            tmpnCollocPts++;
            // Add equation for Derivative
            if (pAmp->VersionCollocPtsInter[i] == 2 ){
              gsl_vector_set(b, tmpnCollocPts, pAmp->CollocationPointsValuesAmplitudeInter[tmpnCollocPts]);
              fpower = 1./fcollpoint;
              for(INT4 j = 0; j < nCollocPtsInterAmp; j++){
                  REAL8 derivative = (j - seven_sixths) * fpower * fcollpoint_m_seven_sixths;
                  gsl_matrix_set(A, tmpnCollocPts, j, derivative);
                  fpower *= fcollpoint;
              }
              tmpnCollocPts++;
            }
        } /* End non-zero if statement */
    } /* End of loop over number of free coefficients. System of equations setup. */

    /* tmpnCollocPts must be the same than the number of free coefficients and to pAmp->nCoefficientsInter */
    if (tmpnCollocPts != pAmp->nCoefficientsInter ){
        /* Free gsl variables */
        gsl_vector_free(b);
        gsl_vector_free(x);
        gsl_matrix_free(A);
        gsl_permutation_free(p);
        XLAL_ERROR_VOID(XLAL_EFUNC, "IMRPhenomXHM_Intermediate_Amp_Coefficients failed. Inconsistent number of collocation points (%i) and free parameters (%i).", tmpnCollocPts, pAmp->nCoefficientsInter);
    }

    /* We now solve the system A x = b via an LU decomposition. x is the solution vector */
    gsl_linalg_LU_decomp(A, p, &signum);
    gsl_linalg_LU_solve(A, p, b, x);

    /* The solution corresponds to the coefficients of the ansatz */
    for (UINT2 i = 0; i < pAmp->nCoefficientsInter; i++){
        pAmp->InterCoefficient[i] = gsl_vector_get(x, i);
    }

    /* Free gsl variables */
    gsl_vector_free(b);
    gsl_vector_free(x);
    gsl_matrix_free(A);
    gsl_permutation_free(p);

}


/***********************************************/
/*                                             */
/*                   PHASE                     */
/*                                             */
/***********************************************/

// Fits of phase derivatives at collocation points for each mode

/* Start of Phase Parameter Space Fits */

static double IMRPhenomXHM_Inter_Phase_21_p1(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,eta5,eta6,S2,S3,S4;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4 = pow(S,4);
            double noSpin = 4045.84 + 7.63226/eta - 1956.93*eta - 23428.1*eta2 + 369153.*eta3 - 2.28832e6*eta4 + 6.82533e6*eta5 - 7.86254e6*eta6;
            double eqSpin = - 347.273*S + 83.5428*S2 - 355.67*S3 + (4.44457*S + 16.5548*S2 + 13.6971*S3)/eta + eta*( - 79.761*S - 355.299*S2 + 1114.51*S3 - 1077.75*S4) + 92.6654*S4 + eta2*(- 619.837*S - 722.787*S2 + 2392.73*S3 + 2689.18*S4);
            double uneqSpin =  ( 918.976*pWF->chi1L*sqrt(1.-4.*eta) - 918.976*pWF->chi2L*sqrt(1.-4.*eta))*eta + ( 91.7679*pWF->chi1L*sqrt(1.-4.*eta) - 91.7679*pWF->chi2L*sqrt(1.-4.*eta))*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_21_p1: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_Inter_Phase_33_p1(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,eta5,eta6,S2;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            S2 = pow(S,2);
            double noSpin = 4360.19 + 4.27128/eta - 8727.4*eta + 18485.9*eta2 + 371303.00000000006*eta3 - 3.22792e6*eta4 + 1.01799e7*eta5 - 1.15659e7*eta6;
            double eqSpin = ((11.6635 - 251.579*eta - 3255.6400000000003*eta2 + 19614.6*eta3 - 34860.2*eta4)*S + (14.8017 + 204.025*eta - 5421.92*eta2 + 36587.3*eta3 - 74299.5*eta4)*S2)/(eta);
            double uneqSpin = eta*(223.65100000000004*pWF->chi1L*sqrt(1.-4.*eta)*(3.9201300240106223 + 1.*eta) - 223.65100000000004*pWF->chi2L*sqrt(1.-4.*eta)*(3.9201300240106223 + 1.*eta));
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_33_p1: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_Inter_Phase_32_p1(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,eta5,eta6,S2,S3,S4;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4 = pow(S,4);
            double noSpin = 4414.11 + 4.21564/eta - 10687.8*eta + 58234.6*eta2 - 64068.40000000001*eta3 - 704442.*eta4 + 2.86393e6*eta5 - 3.26362e6*eta6;
            double eqSpin = ((6.39833 - 610.267*eta + 2095.72*eta2 - 3970.89*eta3)*S + (22.956700000000005 - 99.1551*eta + 331.593*eta2 - 794.79*eta3)*S2 + (10.4333 + 43.8812*eta - 541.261*eta2 + 294.289*eta3)*S3 + eta*(106.047 - 1569.0299999999997*eta + 4810.61*eta2)*S4)/(eta);
            double uneqSpin = 132.244*sqrt(1.-4.*eta)*eta*(pWF->chi1L*(6.227738120444028 - 1.*eta) + pWF->chi2L*(-6.227738120444028 + 1.*eta));
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_32_p1: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_Inter_Phase_44_p1(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,eta5,eta6,S2,S3;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            S2 = pow(S,2);
            S3 = pow(S,3);
            double noSpin = 4349.66 + 4.34125/eta - 8202.33*eta + 5534.1*eta2 + 536500.*eta3 - 4.33197e6*eta4 + 1.37792e7*eta5 - 1.60802e7*eta6;
            double eqSpin = ((12.0704 - 528.098*eta + 1822.9100000000003*eta2 - 9349.73*eta3 + 17900.9*eta4)*S + (10.4092 + 253.334*eta - 5452.04*eta2 + 35416.6*eta3 - 71523.*eta4)*S2 + eta*(492.60300000000007 - 9508.5*eta + 57303.4*eta2 - 109418.*eta3)*S3)/(eta);
            double uneqSpin = -262.143*sqrt(1.-4.*eta)*eta*(pWF->chi1L*(-3.0782778864970646 - 1.*eta) + pWF->chi2L*(3.0782778864970646 + 1.*eta));
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_44_p1: version is not valid. Recommended version is 122019.");}
    }
    return total;
}


static double IMRPhenomXHM_Inter_Phase_21_p2(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,eta5,eta6,S2,S3,S4;
            eta2 = pow(eta,2);
            eta3 = eta2*eta;
            eta4 = eta3*eta;
            eta5 = eta3*eta2;
            eta6 = eta4*eta2;
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4= S3*S;
            double noSpin = 3509.09 + 0.91868/eta + 194.72*eta - 27556.2*eta2 + 369153.*eta3 - 2.28832e6*eta4 + 6.82533e6*eta5 - 7.86254e6*eta6;
            double eqSpin = ((0.7083999999999999 - 60.1611*eta + 131.815*eta2 - 619.837*eta3)*S + (6.104720000000001 - 59.2068*eta + 278.588*eta2 - 722.787*eta3)*S2 + (5.7791 + 117.913*eta - 1180.4*eta2 + 2392.73*eta3)*S3 + eta*(92.6654 - 1077.75*eta + 2689.18*eta2)*S4)/(eta);
            double uneqSpin = -91.7679*sqrt(1.-4.*eta)*eta*(pWF->chi1L*(-1.6012352903357276 - 1.*eta) + pWF->chi2L*(1.6012352903357276 + 1.*eta));
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_21_p2: version is not valid.Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_Inter_Phase_33_p2(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,eta5,eta6,S2;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            S2 = pow(S,2);
            double noSpin = 3797.06 + 0.786684/eta - 2397.09*eta - 25514.*eta2 + 518314.99999999994*eta3 - 3.41708e6*eta4 + 1.01799e7*eta5 - 1.15659e7*eta6;
            double eqSpin = ((6.7812399999999995 + 39.4668*eta - 3520.37*eta2 + 19614.6*eta3 - 34860.2*eta4)*S + (4.80384 + 293.215*eta - 5914.61*eta2 + 36587.3*eta3 - 74299.5*eta4)*S2)/(eta);
            double uneqSpin = -223.65100000000004*sqrt(1.-4.*eta)*eta*(pWF->chi1L*(-1.3095134830606614 - 1.*eta) + pWF->chi2L*(1.3095134830606614 + 1.*eta));
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_33_p2: version is not valid.Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_Inter_Phase_32_p2(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,eta5,eta6,S2,S3,S4;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4 = pow(S,4);
            double noSpin = 3980.7 + 0.956703/eta - 6202.38*eta + 29218.1*eta2 + 24484.2*eta3 - 807629.*eta4 + 2.86393e6*eta5 - 3.26362e6*eta6;
            double eqSpin = ((1.92692 - 226.825*eta + 75.246*eta2 + 1291.56*eta3)*S + (15.328700000000001 - 99.1551*eta + 608.328*eta2 - 2402.94*eta3)*S2 + (10.4333 + 43.8812*eta - 541.261*eta2 + 294.289*eta3)*S3 + eta*(106.047 - 1569.0299999999997*eta + 4810.61*eta2)*S4)/(eta);
            double uneqSpin = 132.244*sqrt(1.-4.*eta)*eta*(pWF->chi1L*(2.5769789177580837 - 1.*eta) + pWF->chi2L*(-2.5769789177580837 + 1.*eta));
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_32_p2: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_Inter_Phase_44_p2(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,eta5,eta6,S2,S3;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            S2 = pow(S,2);
            S3 = pow(S,3);
            double noSpin = 3804.19 + 0.66144/eta - 2421.77*eta - 33475.8*eta2 + 665951.*eta3 - 4.50145e6*eta4 + 1.37792e7*eta5 - 1.60802e7*eta6;
            double eqSpin = ((5.83038 - 172.047*eta + 926.576*eta2 - 7676.87*eta3 + 17900.9*eta4)*S + (6.17601 + 253.334*eta - 5672.02*eta2 + 35722.1*eta3 - 71523.*eta4)*S2 + eta*(492.60300000000007 - 9508.5*eta + 57303.4*eta2 - 109418.*eta3)*S3)/(eta);
            double uneqSpin = -262.143*sqrt(1.-4.*eta)*eta*(pWF->chi1L*(-1.0543062374352932 - 1.*eta) + pWF->chi2L*(1.0543062374352932 + 1.*eta));
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_44_p2: version is not valid. Recommended version is 122019.");}
    }
    return total;
}


static double IMRPhenomXHM_Inter_Phase_21_p3(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double delta=sqrt(1.-4.*eta);
            double eta2,eta3,eta4,eta5,eta6,S2,S3,S4;
            eta2 = pow(eta,2);
            eta3 = eta2*eta;
            eta4 = eta2*eta2;
            eta5 = eta2*eta3;
            eta6 = eta3*eta3;
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4 = pow(S,4);
            double noSpin = 3241.68 + 890.016*eta - 28651.9*eta2 + 369153.*eta3 - 2.28832e6*eta4 + 6.82533e6*eta5 - 7.86254e6*eta6;
            double eqSpin = (-2.2484 + 187.641*eta - 619.837*eta2)*S + (3.22603 + 166.323*eta - 722.787*eta2)*S2 + (117.913 - 1094.59*eta + 2392.73*eta2)*S3 + (92.6654 - 1077.75*eta + 2689.18*eta2)*S4;
            double uneqSpin = 91.7679*(pWF->dchi)*delta*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_21_p3: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_Inter_Phase_33_p3(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,eta5,eta6,S2;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            S2 = pow(S,2);
            double noSpin = 3321.83 + 1796.03*eta - 52406.1*eta2 + 605028.*eta3 - 3.52532e6*eta4 + 1.01799e7*eta5 - 1.15659e7*eta6;
            double eqSpin = (223.601 - 3714.77*eta + 19614.6*eta2 - 34860.2*eta3)*S + (314.317 - 5906.46*eta + 36587.3*eta2 - 74299.5*eta3)*S2;
            double uneqSpin = 223.651*(pWF->dchi)*sqrt(1.-4.*eta)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_33_p3: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_Inter_Phase_32_p3(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,eta5,eta6,S2,S3,S4;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4 = pow(S,4);
            double noSpin = 3416.57 + 2308.63*eta - 84042.9*eta2 + 1.01936e6*eta3 - 6.0644e6*eta4 + 1.76399e7*eta5 - 2.0065e7*eta6;
            double eqSpin = (24.6295 - 282.354*eta - 2582.55*eta2 + 12750.*eta3)*S + (433.675 - 8775.86*eta + 56407.8*eta2 - 114798.*eta3)*S2 + (559.705 - 10627.4*eta + 61581.*eta2 - 114029.*eta3)*S3 + (106.047 - 1569.03*eta + 4810.61*eta2)*S4;
            double uneqSpin = 63.9466*pWF->dchi*sqrt(1.-4.*eta)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_32_p3: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_Inter_Phase_44_p3(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4, eta5, eta6, S2, S3;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            S2 = pow(S,2);
            S3 = pow(S,3);
            double noSpin = 3308.97 + 2353.58*eta - 66340.1*eta2 + 777272.*eta3 - 4.64438e6*eta4 + 1.37792e7*eta5 - 1.60802e7*eta6;
            double eqSpin = (-21.5697 + 926.576*eta - 7989.26*eta2 + 17900.9*eta3)*S + (353.539 - 6403.24*eta + 37599.5*eta2 - 71523.*eta3)*S2 + (492.603 - 9508.5*eta + 57303.4*eta2 - 109418.*eta3)*S3;
            double uneqSpin = 262.143*(pWF->dchi)*sqrt(1.-4.*eta)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_44_p3: version is not valid.Recommended version is 122019.");}
    }
    return total;
}


static double IMRPhenomXHM_Inter_Phase_21_p4(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,eta5,eta6,S2,S3,S4;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4 = pow(S,4);
            double noSpin = 3160.88 + 974.355*eta - 28932.5*eta2 + 369780.*eta3 - 2.28832e6*eta4 + 6.82533e6*eta5 - 7.86254e6*eta6;
            double eqSpin = (26.3355 - 196.851*eta + 438.401*eta2)*S + (45.9957 - 256.248*eta + 117.563*eta2)*S2 + (-20.0261 + 467.057*eta - 1613.*eta2)*S3 + (-61.7446 + 577.057*eta - 1096.81*eta2)*S4;
            double uneqSpin = 65.3326*pWF->dchi*sqrt(1.-4.*eta)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_21_p4: version is not valid.Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_Inter_Phase_33_p4(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,eta5, eta6,S2,S3;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            S2 = pow(S,2);
            S3 = pow(S,3);
            double noSpin = 3239.44 - 661.15*eta + 5139.79*eta2 + 3456.2*eta3 - 248477.*eta4 + 1.17255e6*eta5 - 1.70363e6*eta6;
            double eqSpin = (225.859 - 4150.09*eta + 24364.*eta2 - 46537.3*eta3)*S + (35.2439 - 994.971*eta + 8953.98*eta2 - 23603.5*eta3)*S2 + (-310.489 + 5946.15*eta - 35337.1*eta2 + 67102.4*eta3)*S3;
            double uneqSpin = 30.484*pWF->dchi*sqrt(1.-4.*eta)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_33_p4: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_Inter_Phase_32_p4(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,eta5,eta6,S2,S3,S4;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4 = pow(S,4);
            double noSpin = 3307.49 - 476.909*eta - 5980.37*eta2 + 127610.*eta3 - 919108.*eta4 + 2.86393e6*eta5 - 3.26362e6*eta6;
            double eqSpin = (-5.02553 - 282.354*eta + 1291.56*eta2)*S + (-43.8823 + 740.123*eta - 2402.94*eta2)*S2 + (43.8812 - 370.362*eta + 294.289*eta2)*S3 + (106.047 - 1569.03*eta + 4810.61*eta2)*S4;
            double uneqSpin = -132.244*(pWF->dchi)*sqrt(1.-4.*eta)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_32_p4: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_Inter_Phase_44_p4(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,eta5,eta6,S2,S3;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            S2 = pow(S,2);
            S3 = pow(S,3);
            double noSpin = 3245.63 - 928.56*eta + 8463.89*eta2 - 17422.6*eta3 - 165169.*eta4 + 908279.*eta5 - 1.31138e6*eta6;
            double eqSpin = (32.506 - 590.293*eta + 3536.61*eta2 - 6758.52*eta3)*S + (-25.7716 + 738.141*eta - 4867.87*eta2 + 9129.45*eta3)*S2 + (-15.7439 + 620.695*eta - 4679.24*eta2 + 9582.58*eta3)*S3;
            double uneqSpin = 87.0832*pWF->dchi*sqrt(1.-4.*eta)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_44_p4: version is not valid. Recommended version is 122019.");}
    }
    return total;
}


static double IMRPhenomXHM_Inter_Phase_21_p5(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,S2,S3,S4;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4 = pow(S,4);
            double noSpin = 3102.36 + 315.911*eta - 1688.26*eta2 + 3635.76*eta3;
            double eqSpin = (-23.0959 + 320.93*eta - 1029.76*eta2)*S + (-49.5435 + 826.816*eta - 3079.39*eta2)*S2 + (40.7054 - 365.842*eta + 1094.11*eta2)*S3 + (81.8379 - 1243.26*eta + 4689.22*eta2)*S4;
            double uneqSpin = 119.014*(pWF->dchi)*sqrt(1.-4.*eta)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_21_p5: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_Inter_Phase_33_p5(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,eta5,eta6,S2,S3;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            S2 = pow(S,2);
            S3 = pow(S,3);
            double noSpin = 3114.3 + 2143.06*eta - 49428.3*eta2 + 563997.*eta3 - 3.35991e6*eta4 + 9.99745e6*eta5 - 1.17123e7*eta6;
            double eqSpin = (190.051 - 3705.08*eta + 23046.2*eta2 - 46537.3*eta3)*S + (63.6615 - 1414.2*eta + 10166.1*eta2 - 23603.5*eta3)*S2 + (-257.524 + 5179.97*eta - 33001.4*eta2 + 67102.4*eta3)*S3;
            double uneqSpin = 54.9833*pWF->dchi*sqrt(1.-4.*eta)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_33_p5: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_Inter_Phase_32_p5(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,eta5,eta6,eta7,S2,S3,S4;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            eta7 = pow(eta,7);
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4 = pow(S,4);
            double noSpin = 3259.03 - 3967.58*eta + 111203.*eta2 - 1.81883e6*eta3 + 1.73811e7*eta4 - 9.56988e7*eta5 + 2.75056e8*eta6 - 3.15866e8*eta7;
            double eqSpin = (19.7509 - 1104.53*eta + 3810.18*eta2)*S + (-230.07 + 2314.51*eta - 5944.49*eta2)*S2 + (-201.633 + 2183.43*eta - 6233.99*eta2)*S3 + (106.047 - 1569.03*eta + 4810.61*eta2)*S4;
            double uneqSpin = 112.714*pWF->dchi*sqrt(1.-4.*eta)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_32_p5: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_Inter_Phase_44_p5(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,eta5,eta6,eta7,S2,S3;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            eta7 = pow(eta,7);
            S2 = pow(S,2);
            S3 = pow(S,3);
            double noSpin = 3108.38 + 3722.46*eta - 119588.*eta2 + 1.92148e6*eta3 - 1.69796e7*eta4 + 8.39194e7*eta5 - 2.17143e8*eta6 + 2.2829700000000003e8*eta7;
            double eqSpin = (118.319 - 529.854*eta)*eta*S + (21.0314 - 240.648*eta + 516.333*eta2)*S2 + (20.3384 - 356.241*eta + 999.417*eta2)*S3;
            double uneqSpin = 97.1364*(pWF->dchi)*sqrt(1.-4.*eta)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_44_p5: version is not valid. Recommended version is 122019.");}
    }
    return total;
}


static double IMRPhenomXHM_Inter_Phase_21_p6(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,S2,S3,S4;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4 = pow(S,4);
            double noSpin = 3089.18 + 4.89194*eta + 190.008*eta2 - 255.245*eta3;
            double eqSpin = (2.96997 + 57.1612*eta - 432.223*eta2)*S + (-18.8929 + 630.516*eta - 2804.66*eta2)*S2 + (-24.6193 + 549.085*eta2)*S3 + (-12.8798 - 722.674*eta + 3967.43*eta2)*S4;
            double uneqSpin = 74.0984*(pWF->dchi)*sqrt(1.-4.*eta)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_21_p6: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_Inter_Phase_33_p6(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3, eta4, eta5, eta6,S2,S3;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            S2 = pow(S,2);
            S3 = pow(S,3);
            double noSpin = 3111.46 + 384.121*eta - 13003.6*eta2 + 179537.*eta3 - 1.19313e6*eta4 + 3.79886e6*eta5 - 4.64858e6*eta6;
            double eqSpin = (182.864 - 3834.22*eta + 24532.9*eta2 - 50165.9*eta3)*S + (21.0158 - 746.957*eta + 6701.33*eta2 - 17842.3*eta3)*S2 + (-292.855 + 5886.62*eta - 37382.4*eta2 + 75501.8*eta3)*S3;
            double uneqSpin = 75.5162*pWF->dchi*sqrt(1.-4.*eta)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_33_p6: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_Inter_Phase_32_p6(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,eta5,eta6,eta7,S2,S3,S4;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            eta7 = pow(eta,7);
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4 = pow(S,4);
            double noSpin = 3259.03 - 3967.58*eta + 111203.*eta2 - 1.81883e6*eta3 + 1.73811e7*eta4 - 9.56988e7*eta5 + 2.75056e8*eta6 - 3.15866e8*eta7;
            double eqSpin = (19.7509 - 1104.53*eta + 3810.18*eta2)*S + (-230.07 + 2314.51*eta - 5944.49*eta2)*S2 + (-201.633 + 2183.43*eta - 6233.99*eta2)*S3 + (106.047 - 1569.03*eta + 4810.61*eta2)*S4;
            double uneqSpin = 112.714*pWF->dchi*sqrt(1.-4.*eta)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_32_p6: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_Inter_Phase_44_p6(IMRPhenomXWaveformStruct *pWF, int InterPhaseFlag) {
    double total=0.;
    switch (InterPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,eta5,eta6,S2,S3;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            S2 = pow(S,2);
            S3 = pow(S,3);
            double noSpin = 3096.03 + 986.752*eta - 20371.1*eta2 + 220332.*eta3 - 1.31523e6*eta4 + 4.29193e6*eta5 - 6.01179e6*eta6;
            double eqSpin = (-9.96292 - 118.526*eta + 2255.76*eta2 - 6758.52*eta3)*S + (-14.4869 + 370.039*eta - 3605.8*eta2 + 9129.45*eta3)*S2 + (17.0209 + 70.1931*eta - 3070.08*eta2 + 9582.58*eta3)*S3;
            double uneqSpin = 23.0759*pWF->dchi*sqrt(1.-4.*eta)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Inter_Phase_44_p6: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

/* End of Phase Parameter Space Fits */

/************* PHASE ANSATZ ****************/

static double IMRPhenomXHM_Inter_Phase_AnsatzInt(double ff, IMRPhenomX_UsefulPowers *powers_of_f,IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMPhaseCoefficients *pPhase){

  int modeTag = pWFHM->modeTag;
  double invf  = powers_of_f->m_one;
  double invf2 = powers_of_f->m_two;
  double invf3 = powers_of_f->m_three;
  double logfv=powers_of_f->log;
  double phaseIR;
  double fda=pWFHM->fDAMP;
  double frd=pWFHM->fRING;

  if(modeTag!=32)
  /*  5 coefficients */

  phaseIR = pPhase->c0 *ff + (pPhase->c1)*logfv - (pPhase->c2)*invf -1./3.* (pPhase->c4)*invf3 + (pPhase->cL)* atan((ff - frd)/fda);
  else {          /*  6 coefficients */

    invf3 = powers_of_f->m_three;
    /*(c0 + c1 /f + c2 /(f)^2 + c4 /(f)^4 + c3 /f^3 +
    cL fdamp/((fdamp)^2 + (f - fRD )^2))*/
    phaseIR = pPhase->c0 *ff + (pPhase->c1)*logfv - (pPhase->c2)*invf -1./3.* (pPhase->c4)*invf3 -0.5*pPhase->c3*invf2+ (pPhase->cL)* atan((ff - frd)/fda);
  }

  return phaseIR;

}

static double IMRPhenomXHM_Inter_Phase_Ansatz(double ff, IMRPhenomX_UsefulPowers *powers_of_f,IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMPhaseCoefficients *pPhase){

  int modeTag = pWFHM->modeTag;
  double invf  = powers_of_f->m_one;
  double invf2 = powers_of_f->m_two;
  double invf4 = powers_of_f->m_four;
  double dphaseIR;
  double fda=pWFHM->fDAMP;
  double frd=pWFHM->fRING;

  if(modeTag!=32)
  /*  5 coefficients */

  /*(c0 + c1 /f + c2 /(f)^2 + c4 /(f)^4 +
  cL fdamp/((fdamp)^2 + (f - fRD )^2))*/
  dphaseIR = ( pPhase->c0 + (pPhase->c1)*invf + (pPhase->c2)*invf2 + (pPhase->c4)*invf4 + ( (pPhase->cL)* fda/(fda*fda +(ff - frd)*(ff - frd)) ) );
  else {          /*  6 coefficients */

    double invf3 = powers_of_f->m_three;
    /*(c0 + c1 /f + c2 /(f)^2 + c4 /(f)^4 + c3 /f^3 +
    cL fdamp/((fdamp)^2 + (f - fRD )^2))*/
    dphaseIR = ( pPhase->c0 + (pPhase->c1)*invf + (pPhase->c2)*invf2 + (pPhase->c4)*invf4 +
    (pPhase->c3)*invf3 + ( (pPhase->cL)* fda/(fda*fda +(ff - frd)*(ff - frd)) ) );
  }

  return dphaseIR;

}

static double IMRPhenomXHM_Inter_Phase_dAnsatz(double ff,  IMRPhenomX_UsefulPowers *powers_of_f,IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMPhaseCoefficients *pPhase){

  int modeTag = pWFHM->modeTag;
  double invf2 = powers_of_f->m_two;
  double invf3 = powers_of_f->m_three;
  double invf5 = powers_of_f->m_five;
  double d2phaseIR;
  double fda=pWFHM->fDAMP;
  double frd=pWFHM->fRING;

  if(modeTag!=32)
  /*  5 coefficients */

  /*(-((4 c4)/f^5) - (2 c2)/f^3 - c1/f^2 - (
  2 cL fdamp (f - fRD))/(fdamp^2 + (f - fRD)^2)^2)*/
  d2phaseIR =  -(pPhase->c1)*invf2 -2. *(pPhase->c2)*invf3 -4. *(pPhase->c4)*invf5 -2.* (pPhase->cL)* (ff - frd)*fda/pow((fda*fda +(ff - frd)*(ff - frd)),2);
  else {          /*  6 coefficients */

    double invf4 = powers_of_f->m_four;
    /*(-((4 c4)/f^5) - (2 c2)/f^3 - c1/f^2 -3 c3/ f^4- (
    2 cL fdamp (f - fRD))/(fdamp^2 + (f - fRD)^2)^2)*/
    d2phaseIR =-(pPhase->c1)*invf2 -3.*pPhase->c3*invf4-2. *(pPhase->c2)*invf3 -4. *(pPhase->c4)*invf5 -2.* (pPhase->cL)* (ff - frd)*fda/pow((fda*fda +(ff - frd)*(ff - frd)),2);
  }

  return d2phaseIR;

}
