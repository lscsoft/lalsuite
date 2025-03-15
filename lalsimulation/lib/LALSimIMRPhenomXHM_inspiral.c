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
//  LALSimIMRPhenomXHM_inspiral.c
//
//
//  Created by Marta on 06/02/2019.
//

#include "LALSimIMRPhenomXHM_inspiral.h"

/*************************************/
/*                                   */
/*            AMPLITUDE              */
/*                                   */
/*************************************/

/* Fits of the collocation points in the inspiral across parameter space. 3 for each mode: iv1, iv2, iv3
 [iv1,iv2,iv3]=[0.5*f^Ins_lm,0.75*f^Ins_lm,f^Ins_lm]
 For more details about the reconstruction procedure, see Sec. IV.A
 */

 /* These parameter-space fits are documented in the supplementary material of https://dcc.ligo.org/P2000011-v2.
    There are 2 Mathematica notebooks (one for amplitude and one for phase) that read the fits data and automatically generate the C-code below.
    For more information read https://git.ligo.org/waveforms/reviews/imrphenomx/blob/master/documentation/ParspaceFits/README and the documentation in the notebooks. */

// Spin parameter S = chi_PN = chi_eff - 38/113 * eta * (chi1 + chi2)
// chi_eff = (m1*chi1 + m2*chi2)/(m1 + m2)

/* Start of Amp Parameter Space Fits */

static double IMRPhenomXHM_Insp_Amp_21_iv1(IMRPhenomXWaveformStruct *pWF, int InspAmpFlag) {
    double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double eta = pWF->eta;
            double S = pWF->chiPNHat;
            double eta2,eta3,eta4,eta5,S2;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            S2 = pow(S,2);
            double noSpin = sqrt(1.-4.*eta)*(0.037868557189995156 + 0.10740090317702103*eta + 1.963812986867654*eta2 - 16.706455229589558*eta3 + 69.75910808095745*eta4 - 98.3062466823662*eta5);
            double eqSpin = sqrt(1.-4.*eta)*S*(-0.007963757232702219 + 0.10627108779259965*eta - 0.008044970210401218*S + eta2*(-0.4735861262934258 - 0.5985436493302649*S - 0.08217216660522082*S2));
            double uneqSpin = -0.257787704938017*pWF->dchi*eta2*(1. + 8.75928187268504*eta2) - 0.2597503605427412*pWF->dchi*eta2*S;
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
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            total = fabs(chidiff1*eta5*(-3962.5020052272976 + 987.635855365408*S1 - 134.98527058315528*S2) + delta*(19.30531354642419 + 16.6640319856064*eta1 - 120.58166037019478*eta2 + 220.77233521626252*eta3)*sqroot + chidiff1*delta*(31.364509907424765*eta1 - 843.6414532232126*eta2 + 2638.3077554662905*eta3)*sqroot + chidiff1*delta*(32.374226994179054*eta1 - 202.86279451816662*eta2 + 347.1621871204769*eta3)*S1*sqroot + delta*S1*(-16.75726972301224*(1.1787350890261943 - 7.812073811917883*eta1 + 99.47071002831267*eta2 - 500.4821414428368*eta3 + 876.4704270866478*eta4) + 2.3439955698372663*(0.9373952326655807 + 7.176140122833879*eta1 - 279.6409723479635*eta2 + 2178.375177755584*eta3 - 4768.212511142035*eta4)*S1)*sqroot);
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_21_iv1: version %i is not valid.", InspAmpFlag);}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Amp_21_iv2(IMRPhenomXWaveformStruct *pWF, int InspAmpFlag) {
    double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double eta = pWF->eta;
            double S = pWF->chiPNHat;
            double delta=sqrt(1.-4.*eta),eta2,eta3,eta4,S2;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            S2 = pow(S,2);
            double noSpin = sqrt(1.-4.*eta)*(0.05511628628738656 - 0.12579599745414977*eta + 2.831411618302815*eta2 - 14.27268643447161*eta3 + 28.3307320191161*eta4);
            double eqSpin = sqrt(1.-4.*eta)*S*(-0.008692738851491525 + eta*(0.09512553997347649 + 0.116470975986383*S) - 0.009520793625590234*S + eta2*(-0.3409769288480959 - 0.8321002363767336*S - 0.13099477081654226*S2) - 0.006383232900211555*S2);
            double uneqSpin = -0.2962753588645467*pWF->dchi*eta2*(1. + 1.3993978458830476*eta2) - 0.17100612756133535*pWF->dchi*eta2*S*(1. + 18.974303741922743*eta2*delta);
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
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = fabs(chidiff1*eta5*(-2898.9172078672705 + 580.9465034962822*S1 + 22.251142639924076*S2) + delta*(chidiff2*(-18.541685007214625*eta1 + 166.7427445020744*eta2 - 417.5186332459383*eta3) + chidiff1*(41.61457952037761*eta1 - 779.9151607638761*eta2 + 2308.6520892707795*eta3))*sqroot + delta*(11.414934585404561 + 30.883118528233638*eta1 - 260.9979123967537*eta2 + 1046.3187137392433*eta3 - 1556.9475493549746*eta4)*sqroot + delta*S1*(-10.809007068469844*(1.1408749895922659 - 18.140470190766937*eta1 + 368.25127088896744*eta2 - 3064.7291458207815*eta3 + 11501.848278358668*eta4 - 16075.676528787526*eta5) + 1.0088254664333147*(1.2322739396680107 - 192.2461213084741*eta1 + 4257.760834055382*eta2 - 35561.24587952242*eta3 + 130764.22485304279*eta4 - 177907.92440833704*eta5)*S1)*sqroot + delta*(chidiff1*(36.88578491943111*eta1 - 321.2569602623214*eta2 + 748.6659668096737*eta3)*S1 + chidiff1*(-95.42418611585117*eta1 + 1217.338674959742*eta2 - 3656.192371615541*eta3)*S2)*sqroot);
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_21_iv2: version %i is not valid.", InspAmpFlag);}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Amp_21_iv3(IMRPhenomXWaveformStruct *pWF, int InspAmpFlag) {
    double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double eta = pWF->eta;
            double S = pWF->chiPNHat;
            double delta=sqrt(1.-4.*eta),eta2,S2;
            eta2 = pow(eta,2);
            S2 = pow(S,2);
            double noSpin = sqrt(1.-4.*eta)*(0.059110044024271766 - 0.0024538774422098405*eta + 0.2428578654261086*eta2);
            double eqSpin = sqrt(1.-4.*eta)*S*(-0.007044339356171243 - 0.006952154764487417*S + eta2*(-0.016643018304732624 - 0.12702579620537421*S + 0.004623467175906347*S2) - 0.007685497720848461*S2);
            double uneqSpin = -0.3172310538516028*pWF->dchi*(1. - 2.9155919835488024*eta2)*eta2 - 0.11975485688200693*pWF->dchi*eta2*S*(1. + 17.27626751837825*eta2*delta);
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
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            total = fabs(chidiff1*eta5*(-2282.9983216879655 + 157.94791186394787*S1 + 16.379731479465033*S2) + chidiff1*delta*(21.935833431534224*eta1 - 460.7130131927895*eta2 + 1350.476411541137*eta3)*sqroot + delta*(5.390240326328237 + 69.01761987509603*eta1 - 568.0027716789259*eta2 + 2435.4098320959706*eta3 - 3914.3390484239667*eta4)*sqroot + chidiff1*delta*(29.731007410186827*eta1 - 372.09609843131386*eta2 + 1034.4897198648962*eta3)*S1*sqroot + delta*S1*(-7.1976397556450715*(0.7603360145475428 - 6.587249958654174*eta1 + 120.87934060776237*eta2 - 635.1835857158857*eta3 + 1109.0598539312573*eta4) - 0.0811847192323969*(7.951454648295709 + 517.4039644814231*eta1 - 9548.970156895082*eta2 + 52586.63520999897*eta3 - 93272.17990295641*eta4)*S1 - 0.28384547935698246*(-0.8870770459576875 + 180.0378964169756*eta1 - 2707.9572896559484*eta2 + 14158.178124971111*eta3 - 24507.800226675925*eta4)*S2)*sqroot);
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_21_iv3: version %i is not valid.", InspAmpFlag);}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Amp_33_iv1(IMRPhenomXWaveformStruct *pWF, int InspAmpFlag) {
    double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double eta = pWF->eta;
            double S = pWF->chiPNHat;
            double eta2,S2;
            eta2 = pow(eta,2);
            S2 = pow(S,2);
            double noSpin = (sqrt(1.-4.*eta)*(-0.056586690934283326 - 0.14374841547279146*eta + 0.5584776628959615*eta2))/(-0.3996185676368123 + eta);
            double eqSpin = sqrt(1.-4.*eta)*S*((0.056042044149691175 + 0.12482426029674777*S)*S + eta*(2.1108074577110343 - 1.7827773156978863*S2) + eta2*(-7.657635515668849 - 0.07646730296478217*S + 5.343277927456605*S2));
            double uneqSpin = 0.45866449225302536*pWF->dchi*(1. - 9.603750707244906*eta2)*eta2;
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
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            total = chidiff1*eta5*(155.1434307076563 + 26.852777193715088*S1 + 1.4157230717300835*S2) + chidiff1*delta*(6.296698171560171*eta1 + 15.81328761563562*eta2 - 141.85538063933927*eta3)*sqroot + delta*(20.94372147101354 + 68.14577638017842*eta1 - 898.470298591732*eta2 + 4598.64854748635*eta3 - 8113.199260593833*eta4)*sqroot + chidiff1*delta*(29.221863857271703*eta1 - 348.1658322276406*eta2 + 965.4670353331536*eta3)*S1*sqroot + delta*S1*(-9.753610761811967*(1.7819678168496158 - 44.07982999150369*eta1 + 750.8933447725581*eta2 - 5652.44754829634*eta3 + 19794.855873435758*eta4 - 26407.40988450443*eta5) + 0.014210376114848208*(-196.97328616330392 + 7264.159472864562*eta1 - 125763.47850622259*eta2 + 1.1458022059130718e6*eta3 - 4.948175330328345e6*eta4 + 7.911048294733888e6*eta5)*S1 - 0.26859293613553986*(-8.029069605349488 + 888.7768796633982*eta1 - 16664.276483466252*eta2 + 128973.72291098491*eta3 - 462437.2690007375*eta4 + 639989.1197424605*eta5)*S2)*sqroot;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_33_iv1: version %i is not valid.", InspAmpFlag);}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Amp_33_iv2(IMRPhenomXWaveformStruct *pWF, int InspAmpFlag) {
    double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double eta = pWF->eta;
            double S = pWF->chiPNHat;
            double eta2,eta3,eta4,eta5,eta6,S2;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            S2 = pow(S,2);
            double noSpin = sqrt(1.-4.*eta)*(0.2137734510411439 - 0.7692194209223682*eta + 26.10570221351058*eta2 - 316.0643979123107*eta3 + 2090.9063511488234*eta4 - 6897.3285171507105*eta5 + 8968.893362362503*eta6);
            double eqSpin = sqrt(1.-4.*eta)*S*(0.018546836505210842 + 0.05924304311104228*S + eta*(1.6484440612224325 - 0.4683932646001618*S - 2.110311135456494*S2) + 0.10701786057882816*S2 + eta2*(-6.51575737684721 + 1.6692205620001157*S + 8.351789152096782*S2));
            double uneqSpin = 0.3929315188124088*pWF->dchi*(1. - 11.289452844364227*eta2)*eta2;
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
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            total = chidiff1*eta5*(161.62678370819597 + 37.141092711336846*S1 - 0.16889712161410445*S2) + chidiff1*delta*(3.4895829486899825*eta1 + 51.07954458810889*eta2 - 249.71072528701757*eta3)*sqroot + delta*(12.501397517602173 + 35.75290806646574*eta1 - 357.6437296928763*eta2 + 1773.8883882162215*eta3 - 3100.2396041211605*eta4)*sqroot + chidiff1*delta*(13.854211287141906*eta1 - 135.54916401086845*eta2 + 327.2467193417936*eta3)*S1*sqroot + delta*S1*(-5.2580116732827085*(1.7794900975289085 - 48.20753331991333*eta1 + 861.1650630146937*eta2 - 6879.681319382729*eta3 + 25678.53964955809*eta4 - 36383.824902258915*eta5) + 0.028627002336747746*(-50.57295946557892 + 734.7581857539398*eta1 - 2287.0465658878725*eta2 + 15062.821881048358*eta3 - 168311.2370167227*eta4 + 454655.37836367317*eta5)*S1 - 0.15528289788512326*(-12.738184090548508 + 1129.44485109116*eta1 - 25091.14888164863*eta2 + 231384.03447562453*eta3 - 953010.5908118751*eta4 + 1.4516597366230418e6*eta5)*S2)*sqroot;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_33_iv2: version %i is not valid.", InspAmpFlag);}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Amp_33_iv3(IMRPhenomXWaveformStruct *pWF, int InspAmpFlag) {
    double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double eta = pWF->eta;
            double S = pWF->chiPNHat;
            double eta2,eta3,eta4,eta5,eta6,S2;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            S2 = pow(S,2);
            double noSpin = sqrt(1.-4.*eta)*(0.2363760327127446 + 0.2855410252403732*eta - 10.159877125359897*eta2 + 162.65372389693505*eta3 - 1154.7315106095564*eta4 + 3952.61320206691*eta5 - 5207.67472857814*eta6);
            double eqSpin = sqrt(1.-4.*eta)*S*(0.04573095188775319 + 0.048249943132325494*S + eta*(0.15922377052827502 - 0.1837289613228469*S - 0.2834348500565196*S2) + 0.052963737236081304*S2);
            double uneqSpin = 0.25187274502769835*pWF->dchi*(1. - 12.172961866410864*eta2)*eta2;
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
            double eta6 = eta1 * eta5;
            double eta7 = eta1 * eta6;
            double S1 = S;
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = chidiff1*delta*(-0.5869777957488564*eta1 + 32.65536124256588*eta2 - 110.10276573567405*eta3) + chidiff1*delta*(3.524800489907584*eta1 - 40.26479860265549*eta2 + 113.77466499598913*eta3)*S1 + delta*S1*(-1.2846335585108297*(0.09991079016763821 + 1.37856806162599*eta1 + 23.26434219690476*eta2 - 34.842921754693386*eta3 - 70.83896459998664*eta4) - 0.03496714763391888*(-0.230558571912664 + 188.38585449575902*eta1 - 3736.1574640444287*eta2 + 22714.70643022915*eta3 - 43221.0453556626*eta4)*S1) + chidiff1*eta7*(2667.3441342894776 + 47.94869769580204*chidiff2 + 793.5988192446642*S1 + 293.89657731755483*S2) + delta*(5.148353856800232 + 148.98231189649468*eta1 - 2774.5868652930294*eta2 + 29052.156454239772*eta3 - 162498.31493332976*eta4 + 460912.76402476896*eta5 - 521279.50781871413*eta6)*sqroot;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_33_iv3: version %i is not valid.", InspAmpFlag);}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Amp_32_iv1(IMRPhenomXWaveformStruct *pWF, int InspAmpFlag) {
  double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double eta = pWF->eta;
            double S = pWF->chiPNHat;
            double delta=sqrt(1.-4.*eta),eta2,eta3,eta4,eta5,eta6,eta8,S2,S3;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            eta8 = pow(eta,8);
            S2 = pow(S,2);
            S3 = pow(S,3);
            double noSpin = sqrt(1. - 3.*eta)*(0.019069933430190773 - 0.19396651989685837*eta + 11.95224600241255*eta2 - 158.90113442757382*eta3 + 1046.65239329071*eta4 - 3476.940285294999*eta5 + 4707.249209858949*eta6);
            double eqSpin = sqrt(1. - 3.*eta)*S*(0.0046910348789512895 + 0.40231360805609434*eta - 0.0038263656140933152*S + 0.018963579407636953*S2 + eta2*(-1.955352354930108 + 2.3753413452420133*S - 0.9085620866763245*S3) + 0.02738043801805805*S3 + eta3*(7.977057990568723 - 7.9259853291789515*S + 0.49784942656123987*S2 + 5.2255665027119145*S3));
            double uneqSpin = 0.058560321425018165*pow(pWF->dchi,2)*(1. - 19.936477485971217*eta2)*eta2 + 1635.4240644598524*pWF->dchi*eta8*delta + 0.2735219358839411*pWF->dchi*eta2*S*delta;
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
            double eta6 = eta1 * eta5;
            double S1 = S;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = (chidiff1*delta*(-0.739317114582042*eta1 - 47.473246070362634*eta2 + 278.9717709112207*eta3 - 566.6420939162068*eta4) + chidiff2*(-0.5873680378268906*eta1 + 6.692187014925888*eta2 - 24.37776782232888*eta3 + 23.783684827838247*eta4))*sqroot + (3.2940434453819694 + 4.94285331708559*eta1 - 343.3143244815765*eta2 + 3585.9269057886418*eta3 - 19279.186145681153*eta4 + 51904.91007211022*eta5 - 55436.68857586653*eta6)*sqroot + chidiff1*delta*(12.488240781993923*eta1 - 209.32038774208385*eta2 + 1160.9833883184604*eta3 - 2069.5349737049073*eta4)*S1*sqroot + S1*(0.6343034651912586*(-2.5844888818001737 + 78.98200041834092*eta1 - 1087.6241783616488*eta2 + 7616.234910399297*eta3 - 24776.529123239357*eta4 + 30602.210950069973*eta5) - 0.062088720220899465*(6.5586380356588565 + 36.01386705325694*eta1 - 3124.4712274775407*eta2 + 33822.437731298516*eta3 - 138572.93700180828*eta4 + 198366.10615196894*eta5)*S1)*sqroot;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_32_iv1: version %i is not valid.", InspAmpFlag);}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Amp_32_iv2(IMRPhenomXWaveformStruct *pWF, int InspAmpFlag) {
    double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double eta = pWF->eta;
            double S = pWF->chiPNHat;
            double delta=sqrt(1.-4.*eta),eta2,eta3,eta4,eta8,S2,S3;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta8 = pow(eta,8);
            S2 = pow(S,2);
            S3 = pow(S,3);
            double noSpin = sqrt(1. - 3.*eta)*(0.024621376891809633 - 0.09692699636236377*eta + 2.7200998230836158*eta2 - 16.160563094841066*eta3 + 32.930430889650836*eta4);
            double eqSpin = sqrt(1. - 3.*eta)*S*(0.008522695567479373 - 1.1104639098529456*eta2 - 0.00362963820787208*S + 0.016978054142418417*S2 + eta*(0.24280554040831698 + 0.15878436411950506*S - 0.1470288177047577*S3) + 0.029465887557447824*S3 + eta3*(4.649438233164449 - 0.7550771176087877*S + 0.3381436950547799*S2 + 2.5663386135613093*S3));
            double uneqSpin = -0.007061187955941243*pow(pWF->dchi,2)*(1. - 2.024701925508361*eta2)*eta2 + 215.06940561269835*pWF->dchi*eta8*delta + 0.1465612311350642*pWF->dchi*eta2*S*delta;
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
            double eta6 = eta1 * eta5;
            double S1 = S;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = (chidiff2*(-0.03940151060321499*eta1 + 1.9034209537174116*eta2 - 8.78587250202154*eta3) + chidiff1*delta*(-1.704299788495861*eta1 - 4.923510922214181*eta2 + 0.36790005839460627*eta3))*sqroot + (2.2911849711339123 - 5.1846950040514335*eta1 + 60.10368251688146*eta2 - 1139.110227749627*eta3 + 7970.929280907627*eta4 - 25472.73682092519*eta5 + 30950.67053883646*eta6)*sqroot + S1*(0.7718201508695763*(-1.3012906461000349 + 26.432880113146012*eta1 - 186.5001124789369*eta2 + 712.9101229418721*eta3 - 970.2126139442341*eta4) + 0.04832734931068797*(-5.9999628512498315 + 78.98681284391004*eta1 + 1.8360177574514709*eta2 - 2537.636347529708*eta3 + 6858.003573909322*eta4)*S1)*sqroot;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_32_iv2: version %i is not valid.", InspAmpFlag);}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Amp_32_iv3(IMRPhenomXWaveformStruct *pWF, int InspAmpFlag) {
    double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double eta = pWF->eta;
            double S = pWF->chiPNHat;
            double delta=sqrt(1.-4.*eta),eta2,eta3,eta8,S2,S3;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta8 = pow(eta,8);
            S2 = pow(S,2);
            S3 = pow(S,3);
            double noSpin = (sqrt(1. - 3.*eta)*(-0.006150151041614737 + 0.017454430190035*eta + 0.02620962593739105*eta2 - 0.019043090896351363*eta3))/(-0.2655505633361449 + eta);
            double eqSpin = sqrt(1. - 3.*eta)*S*(0.011073381681404716 + 0.00347699923233349*S + eta*S*(0.05592992411391443 - 0.15666140197050316*S2) + 0.012079324401547036*S2 + eta2*(0.5440307361144313 - 0.008730335213434078*S + 0.04615964369925028*S2 + 0.6703688097531089*S3) + 0.016323101357296865*S3);
            double uneqSpin = -0.020140175824954427*pow(pWF->dchi,2)*(1. - 12.675522774051249*eta2)*eta2 - 417.3604094454253*pWF->dchi*eta8*delta + 0.10464021067936538*pWF->dchi*eta2*S*delta;
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
            double eta6 = eta1 * eta5;
            double S1 = S;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = (chidiff2*(-0.6358511175987503*eta1 + 5.555088747533164*eta2 - 14.078156877577733*eta3) + chidiff1*delta*(0.23205448591711159*eta1 - 19.46049432345157*eta2 + 36.20685853857613*eta3))*sqroot + (1.1525594672495008 + 7.380126197972549*eta1 - 17.51265776660515*eta2 - 976.9940395257111*eta3 + 8880.536804741967*eta4 - 30849.228936891763*eta5 + 38785.53683146884*eta6)*sqroot + chidiff1*delta*(1.904350804857431*eta1 - 25.565242391371093*eta2 + 80.67120303906654*eta3)*S1*sqroot + S1*(0.785171689871352*(-0.4634745514643032 + 18.70856733065619*eta1 - 167.9231114864569*eta2 + 744.7699462372949*eta3 - 1115.008825153004*eta4) + 0.13469300326662165*(-2.7311391326835133 + 72.17373498208947*eta1 - 483.7040402103785*eta2 + 1136.8367114738041*eta3 - 472.02962341590774*eta4)*S1)*sqroot;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_32_iv3: version %i is not valid.", InspAmpFlag);}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Amp_44_iv1(IMRPhenomXWaveformStruct *pWF, int InspAmpFlag) {
    double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double eta = pWF->eta;
            double S = pWF->chiPNHat;
            double eta2,eta3,eta4,S2;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            S2 = pow(S,2);
            double noSpin = sqrt(1. - 3.*eta)*(0.06190013067931406 + 0.1928897813606222*eta + 1.9024723168424225*eta2 - 15.988716302668415*eta3 + 35.21461767354364*eta4);
            double eqSpin = sqrt(1. - 3.*eta)*S*(0.011454874900772544 + 0.044702230915643903*S + eta*(0.6600413908621988 + 0.12149520289658673*S - 0.4482406547006759*S2) + 0.07327810908370004*S2 + eta2*(-2.1705970511116486 - 0.6512813450832168*S + 1.1237234702682313*S2));
            double uneqSpin = 0.4766851579723911*pWF->dchi*(1. - 15.950025762198988*eta2)*eta2 + 0.127900699645338*pow(pWF->dchi,2)*(1. - 15.79329306044842*eta2)*eta2;
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
            double eta6 = eta1 * eta5;
            double S1 = S;
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = (chidiff1*delta*(0.5697308729057493*eta1 + 8.895576813118867*eta2 - 34.98399465240273*eta3) + chidiff2*(1.6370346538130884*eta1 - 14.597095790380884*eta2 + 33.182723737396294*eta3))*sqroot + (5.2601381002242595 - 3.557926105832778*eta1 - 138.9749850448088*eta2 + 603.7453704122706*eta3 - 923.5495700703648*eta4)*sqroot + S1*(-0.41839636169678796*(5.143510231379954 + 104.62892421207803*eta1 - 4232.508174045782*eta2 + 50694.024801783446*eta3 - 283097.33358214336*eta4 + 758333.2655404843*eta5 - 788783.0559069642*eta6) - 0.05653522061311774*(5.605483124564013 + 694.00652410087*eta1 - 17551.398321516353*eta2 + 165236.6480734229*eta3 - 761661.9645651339*eta4 + 1.7440315410044065e6*eta5 - 1.6010489769238676e6*eta6)*S1 - 0.023693246676754775*(16.437107575918503 - 2911.2154288136217*eta1 + 89338.32554683842*eta2 - 1.0803340811860575e6*eta3 + 6.255666490084672e6*eta4 - 1.7434160932177313e7*eta5 + 1.883460394974573e7*eta6)*S2)*sqroot;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_44_iv1: version %i is not valid.", InspAmpFlag);}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Amp_44_iv2(IMRPhenomXWaveformStruct *pWF, int InspAmpFlag) {
    double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double eta = pWF->eta;
            double S = pWF->chiPNHat;
            double delta=sqrt(1.-4.*eta),eta2,eta3,S2,S3;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            S2 = pow(S,2);
            S3 = pow(S,3);
            double noSpin = 0.08406011695496626 - 0.1469952725049322*eta + 0.2997223283799925*eta2 - 1.2910560244510723*eta3;
            double eqSpin = (0.023924074703897662 + 0.26110236039648027*eta - 1.1536009170220438*eta2)*S + (0.04479727299752669 - 0.1439868858871802*eta + 0.05736387085230215*eta2)*S2 + (0.06028104440131858 - 0.4759412992529712*eta + 1.1090751649419717*eta2)*S3;
            double uneqSpin = 0.10346324686812074*pow(pWF->dchi,2)*(1. - 16.135903382018213*eta2)*eta2 + 0.2648241309154185*pWF->dchi*eta2*delta;
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
            double eta6 = eta1 * eta5;
            double S1 = S;
            double S2 = S1 * S1;
            double S3 = S1 * S2;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = (chidiff2*(-0.8318312659717388*eta1 + 7.6541168007977864*eta2 - 16.648660653220123*eta3) + chidiff1*delta*(2.214478316304753*eta1 - 7.028104574328955*eta2 + 5.56587823143958*eta3))*sqroot + (3.173191054680422 + 6.707695566702527*eta1 - 155.22519772642607*eta2 + 604.0067075996933*eta3 - 876.5048298377644*eta4)*sqroot + chidiff1*delta*(4.749663394334708*eta1 - 42.62996105525792*eta2 + 97.01712147349483*eta3)*S1*sqroot + S1*(-0.2627203100303006*(6.460396349297595 - 52.82425783851536*eta1 - 552.1725902144143*eta2 + 12546.255587592654*eta3 - 81525.50289542897*eta4 + 227254.37897941095*eta5 - 234487.3875219032*eta6) - 0.008424003742397579*(-109.26773035716548 + 15514.571912666677*eta1 - 408022.6805482195*eta2 + 4.620165968920881e6*eta3 - 2.6446950627957724e7*eta4 + 7.539643948937692e7*eta5 - 8.510662871580401e7*eta6)*S1 - 0.008830881730801855*(-37.49992494976597 + 1359.7883958101172*eta1 - 23328.560285901796*eta2 + 260027.4121353132*eta3 - 1.723865744472182e6*eta4 + 5.858455766230802e6*eta5 - 7.756341721552802e6*eta6)*S2 - 0.027167813927224657*(34.281932237450256 - 3312.7658728016568*eta1 + 84126.14531363266*eta2 - 956052.0170024392*eta3 + 5.570748509263883e6*eta4 - 1.6270212243584689e7*eta5 + 1.8855858173287075e7*eta6)*S3)*sqroot;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_44_iv2: version %i is not valid.", InspAmpFlag);}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Amp_44_iv3(IMRPhenomXWaveformStruct *pWF, int InspAmpFlag) {
    double total=0;
    switch (InspAmpFlag){
        case 122018:{
            double eta = pWF->eta;
            double S = pWF->chiPNHat;
            double delta=sqrt(1.-4.*eta),eta2,eta3,eta4,eta5,S2;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            S2 = pow(S,2);
            double noSpin = 0.08212436946985402 - 0.025332770704783136*eta - 3.2466088293309885*eta2 + 28.404235115663706*eta3 - 111.36325359782991*eta4 + 157.05954559045156*eta5;
            double eqSpin = S*(0.03488890057062679 + 0.039491331923244756*S + eta*(-0.08968833480313292 - 0.12754920943544915*S - 0.11199012099701576*S2) + 0.034468577523793176*S2);
            double uneqSpin = 0.2062291124580944*pWF->dchi*eta2*delta;
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
            double eta6 = eta1 * eta5;
            double S1 = S;
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = (chidiff1*delta*(1.4739380748149558*eta1 + 0.06541707987699942*eta2 - 9.473290540936633*eta3) + chidiff2*(-0.3640838331639651*eta1 + 3.7369795937033756*eta2 - 8.709159662885131*eta3))*sqroot + (1.7335503724888923 + 12.656614578053683*eta1 - 139.6610487470118*eta2 + 456.78649322753824*eta3 - 599.2709938848282*eta4)*sqroot + chidiff1*delta*(2.3532739003216254*eta1 - 21.37216554136868*eta2 + 53.35003268489743*eta3)*S1*sqroot + S1*(-0.15782329022461472*(6.0309399412954345 - 229.16361598098678*eta1 + 3777.477006415653*eta2 - 31109.307191210424*eta3 + 139319.8239886073*eta4 - 324891.4001578353*eta5 + 307714.3954026392*eta6) - 0.03050157254864058*(4.232861441291087 + 1609.4251694451375*eta1 - 51213.27604422822*eta2 + 612317.1751155312*eta3 - 3.5589766538499263e6*eta4 + 1.0147654212772278e7*eta5 - 1.138861230369246e7*eta6)*S1 - 0.026407497690308382*(-17.184685557542196 + 744.4743953122965*eta1 - 10494.512487701073*eta2 + 66150.52694069289*eta3 - 184787.79377504133*eta4 + 148102.4257785174*eta5 + 128167.89151782403*eta6)*S2)*sqroot;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_Insp_Amp_44_iv3: version %i is not valid.", InspAmpFlag);}
    }
    return total;
}

/* End of Amp Parameter Space Fits */


/***********************************************/
/*     Pseudo PN Amplitude Coefficients    */
/***********************************************/

/* The whole inspiral ansatz is given by

    PNAnsatz(f) + pseudo-PN(f),  with pseudo-PN(f) = rho1 *(f/fcutInsp)^(7/3) + rho2(f/fcutInsp)^(8/3) + rho3(f/fcutInsp)^(9/3).

  The coefficients are computed demanding that the pseudo-PN part at the three collocation points frequencies returns the actual collocation points: iv1, iv2, iv3.

*/

/* The input arguments for these rho coefficients are the value of the collocation points: v1, v2, v3, useful powers of the frequencies of these points f1, f2, f3 and the cutting frequency for the
inspiral fcutInsp. The values v1, v2, v3 are given by the fit of the collocation points minus the value of the PNAnsatz at that frequency. */

/* There are 3 functions, one for each pseudoPN coefficient rho1,2,3.
Every function contains three possible results, depending on if you reconstruct with 3, 2 or 1 collocation points.
With 2 colloc points rho3 = 0, and with 1 colloc point rho2=rho3=0 */

/* Get Pseudo PN Amp coefficient rho1 at f^(7/3) */
static double IMRPhenomXHM_Inspiral_Amp_rho1(double v1, double v2, double v3, IMRPhenomX_UsefulPowers *powers_of_fcutInsp, IMRPhenomX_UsefulPowers *powers_of_f1, IMRPhenomX_UsefulPowers *powers_of_f2, IMRPhenomX_UsefulPowers *powers_of_f3, IMRPhenomXHMWaveformStruct *pWFHM){
  double retVal = 0.;
  switch(pWFHM->IMRPhenomXHMInspiralAmpVersion){
    case 3: // 3 inspiral collocation points
    {
      // default version
      retVal = (powers_of_fcutInsp->seven_thirds*(-(powers_of_f1->three*powers_of_f3->eight_thirds*v2) + powers_of_f1->eight_thirds*powers_of_f3->three*v2 + powers_of_f2->three*(powers_of_f3->eight_thirds*v1 - powers_of_f1->eight_thirds*v3) + powers_of_f2->eight_thirds*(-(powers_of_f3->three*v1) + powers_of_f1->three*v3)))/(powers_of_f1->seven_thirds*(powers_of_f1->one_third - powers_of_f2->one_third)*powers_of_f2->seven_thirds*(powers_of_f1->one_third - powers_of_f3->one_third)*(powers_of_f2->one_third - powers_of_f3->one_third)*powers_of_f3->seven_thirds);
      break;
    }
    case 2:  // 2 inspiral collocation points
    {
      retVal=(powers_of_fcutInsp->seven_thirds*(-(powers_of_f2->eight_thirds*v1) + powers_of_f1->eight_thirds*v2))/(powers_of_f1->seven_thirds*(powers_of_f1->one_third - powers_of_f2->one_third)*powers_of_f2->seven_thirds);
      break;
    }
    case 1:  // 1 inspiral collocation points
    {
      retVal=(powers_of_fcutInsp->seven_thirds*v1)/(powers_of_f1->seven_thirds);
      break;
    }
    case 0:  // No collocation points, reconstruct with PN only
    {
      retVal=0;
      break;
    }
    default: {
      XLALPrintError("Error in IMRPhenomXHM_Inspiral_Amp_rho1: version is not valid. Versions available are 0,1,2,3.\n");
      XLAL_ERROR(XLAL_EDOM, "IMRPhenomXHM_Inspiral_Amp_rho1 version is not valid.  Aborting. Versions available are 0,1,2,3.\n");
    }
  }
  return retVal;

}

/* Get Pseudo PN Amp coefficient rho2 at f^(8/3) */
static double IMRPhenomXHM_Inspiral_Amp_rho2(double v1, double v2, double v3, IMRPhenomX_UsefulPowers *powers_of_fcutInsp, IMRPhenomX_UsefulPowers *powers_of_f1, IMRPhenomX_UsefulPowers *powers_of_f2, IMRPhenomX_UsefulPowers *powers_of_f3, IMRPhenomXHMWaveformStruct *pWFHM){

  double retVal = 0;
  switch(pWFHM->IMRPhenomXHMInspiralAmpVersion){
    case 3:  // 3 inspiral collocation points
    {
      // default version
      retVal=(-(powers_of_fcutInsp->eight_thirds*(-(powers_of_f1->three*powers_of_f3->seven_thirds*v2) + powers_of_f1->seven_thirds*powers_of_f3->three*v2 + powers_of_f2->three*(powers_of_f3->seven_thirds*v1 - powers_of_f1->seven_thirds*v3) + powers_of_f2->seven_thirds*(-(powers_of_f3->three*v1) + powers_of_f1->three*v3))))/(powers_of_f1->seven_thirds*(powers_of_f1->one_third - powers_of_f2->one_third)*powers_of_f2->seven_thirds*(powers_of_f1->one_third - powers_of_f3->one_third)*(powers_of_f2->one_third - powers_of_f3->one_third)*powers_of_f3->seven_thirds);
      break;
    }
    case 2: // 2 inspiral collocation points
    {

      retVal=(-(powers_of_fcutInsp->eight_thirds*(-(powers_of_f2->seven_thirds*v1) + powers_of_f1->seven_thirds*v2)))/(powers_of_f1->seven_thirds*(powers_of_f1->one_third - powers_of_f2->one_third)*powers_of_f2->seven_thirds);
      break;
    }
    case 1:  // 1 inspiral collocation points
    {
      retVal = 0;
      break;
    }
    case 0:  // No collocation points, reconstruct with PN only
    {
      retVal=0;
      break;
    }
    default: {XLALPrintError("Error in IMRPhenomXHM_Inspiral_Amp_rho2: version is not valid. Versions avilable are 0,1,2,3.\n");XLAL_ERROR(XLAL_EDOM, "IMRPhenomXHM_Inspiral_Amp_rho2 version is not valid. Aborting. Versions available are 0,1,2,3.\n");}
  }
  return retVal;

}

/* Get Pseudo PN Amp coefficient rho3 at f^(9/3) */
static double IMRPhenomXHM_Inspiral_Amp_rho3(double v1, double v2, double v3, IMRPhenomX_UsefulPowers *powers_of_fcutInsp, IMRPhenomX_UsefulPowers *powers_of_f1, IMRPhenomX_UsefulPowers *powers_of_f2, IMRPhenomX_UsefulPowers *powers_of_f3, IMRPhenomXHMWaveformStruct *pWFHM){

  double retVal = 0;
  switch(pWFHM->IMRPhenomXHMInspiralAmpVersion){
    case 3:   // 3 inspiral collocation points
    {
      // default version
      retVal=(powers_of_fcutInsp->three*(powers_of_f1->seven_thirds*(-powers_of_f1->one_third + powers_of_f3->one_third)*powers_of_f3->seven_thirds*v2 + powers_of_f2->seven_thirds*(-(powers_of_f3->eight_thirds*v1) + powers_of_f1->eight_thirds*v3) + powers_of_f2->eight_thirds*(powers_of_f3->seven_thirds*v1 - powers_of_f1->seven_thirds*v3)))/(powers_of_f1->seven_thirds*(powers_of_f1->one_third - powers_of_f2->one_third)*powers_of_f2->seven_thirds*(powers_of_f1->one_third - powers_of_f3->one_third)*(powers_of_f2->one_third - powers_of_f3->one_third)*powers_of_f3->seven_thirds);
      break;
    }
    case 2:   // 2 inspiral collocation points
    {
      retVal = 0;
      break;
    }
    case 1:   // 1 inspiral collocation point
    {
      retVal = 0;
      break;
    }
    case 0:  // No collocation points, reconstruct with PN only
    {
      retVal=0;
      break;
    }
      default: {XLALPrintError("Error in IMRPhenomXHM_Inspiral_Amp_rho3: version is not valid.\n");XLAL_ERROR(XLAL_EDOM, "IMRPhenomXHM_Inspiral_Amp_rho3 version is not valid.  Aborting. Versions available: 0,1,2,3.\n");}
  }
  return retVal;
}


/************* INSPIRAL AMPLITUDE ANSATZ ******************/

/*
   The ansatz is built by performing the Stationary Phase Approximation of the Time-Domain Amplitude up to 3PN.
   The result is rexpanded in power series up to 3PN, except for the 21 mode, which is better behaved without the reexpansion.
*/

//Return the Fourier Domain Post-Newtonian ansatz up to 3PN without the pseudoPN terms for a particular frequency
static double IMRPhenomXHM_Inspiral_PNAmp_Ansatz(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMAmpCoefficients *pAmp){

    // The 21 mode is special, is not a power series
    if(pWFHM->useFAmpPN==1){
    return IMRPhenomXHM_Inspiral_PNAmp_21Ansatz(powers_of_Mf, pWFHM, pAmp);
    }

    //This returns the amplitude strain rescaled with the prefactor of the 22 mode: divided by sqrt(2*eta/3.)/pi^(1/6)*Mf^(-7/6)
    double pnAmp;
    pnAmp = cabs(pAmp->pnInitial
        + powers_of_Mf->one_third    * pAmp->pnOneThird
        + powers_of_Mf->two_thirds   * pAmp->pnTwoThirds
        + powers_of_Mf->itself       * pAmp->pnThreeThirds
        + powers_of_Mf->four_thirds  * pAmp->pnFourThirds
        + powers_of_Mf->five_thirds  * pAmp->pnFiveThirds
        + powers_of_Mf->two          * pAmp->pnSixThirds
    );// * (pAmp->PNglobalfactor); // Added
    if (pAmp->InspRescaleFactor == 0){
        pnAmp *= pAmp->PNglobalfactor * powers_of_Mf->m_seven_sixths * pWFHM->ampNorm;
    }
    else if (pAmp->InspRescaleFactor == 1){
       pnAmp *= pAmp->PNglobalfactor;
    }

    return pnAmp;
}

//Return the 21 mode Fourier Domain Post-Newtonian ansatz up to 3PN without the pseudoPN terms for a particular frequency
static double IMRPhenomXHM_Inspiral_PNAmp_21Ansatz(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMAmpCoefficients *pAmp){

  double complex CpnAmpTD; //This is not the real strain of the lm mode. It is the strain rescaled with the prefactor of the 22 mode: divided by sqrt(2*eta/3.)/pi^(1/6)
  double pnAmp, XdotT4, x_to_m_one_four, two_to_m_one_sixths = 0.8908987181403393, three_to_m_one_second = 0.5773502691896257;
  x_to_m_one_four = two_to_m_one_sixths * powers_of_lalpiHM.m_one_sixth * powers_of_Mf->m_one_sixth;
  int InsAmpFlag = pWFHM->IMRPhenomXHMInspiralAmpFreqsVersion;
  switch(InsAmpFlag)
  {
    case 122018:
    {
      // Complex time-domain Post-Newtonina amplitude, power series
      CpnAmpTD = (
         powers_of_Mf->one_third    * pAmp->x05
        + powers_of_Mf->two_thirds  * pAmp->x1
        + powers_of_Mf->itself      * pAmp->x15
        + powers_of_Mf->four_thirds * pAmp->x2
        + powers_of_Mf->five_thirds * pAmp->x25
        + powers_of_Mf->two         * pAmp->x3
      );
      CpnAmpTD = CpnAmpTD*powers_of_Mf->two_thirds*pAmp->PNTDfactor;

      // Value of the the derivative of the PN expansion parameter X given by TaylorT4
      XdotT4 = powers_of_Mf->five_thirds * powers_of_Mf->five_thirds * pAmp->xdot5
      + powers_of_Mf->four          * pAmp->xdot6
      + powers_of_Mf->eight_thirds*powers_of_Mf->five_thirds * pAmp->xdot65
      + powers_of_Mf->seven_thirds*powers_of_Mf->seven_thirds * pAmp->xdot7
      + powers_of_Mf->eight_thirds*powers_of_Mf->seven_thirds * pAmp->xdot75
      + powers_of_Mf->eight_thirds * powers_of_Mf->eight_thirds * pAmp->xdot8
      + (powers_of_Mf->log*2./3. + pAmp->log2pi_two_thirds )*powers_of_Mf->eight_thirds * powers_of_Mf->eight_thirds * pAmp->xdot8Log
      + powers_of_Mf->eight_thirds * powers_of_Mf->eight_thirds * powers_of_Mf->one_third  * pAmp->xdot85;

      // Perform the SPA, multiply time-domain by the phasing factor
      pnAmp = 2. * powers_of_lalpiHM.sqrt * three_to_m_one_second  * cabs(CpnAmpTD) * x_to_m_one_four / sqrt(XdotT4);///powers_of_Mf->m_seven_sixths/pWFHM->ampNorm;//Added last two
      if (pAmp->InspRescaleFactor != 0){
          pnAmp /= RescaleFactor(powers_of_Mf, pAmp, pAmp->InspRescaleFactor);
      }

     break;
    }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomXHM_Inspiral_PNAmp_Ansatz: IMRPhenomXInspiralAmpVersion is not valid. Recommended version is 122018.\n");}
  }
  return pnAmp;
}


//This is the complete Inspiral Amplitude Ansatz: PN ansatz + pseudoPN terms
static double IMRPhenomXHM_Inspiral_Amp_Ansatz(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMAmpCoefficients *pAmp)
{
    double InspAmp; //This is the amplitude strain rescaled with the prefactor of the 22 mode: divided by [sqrt(2*eta/3.)/pi^(1/6) * f^(-7/6)]
    double pseudoterms = 0.;
    if(pWFHM->IMRPhenomXHMInspiralAmpFreqsVersion == 122018){
        INT2 tmp = pAmp->InspRescaleFactor;
        pAmp->InspRescaleFactor = 1;
        InspAmp = IMRPhenomXHM_Inspiral_PNAmp_Ansatz(powers_of_Mf, pWFHM, pAmp)
            + powers_of_Mf->seven_thirds / pAmp->fcutInsp_seven_thirds * pAmp->rho1
            + powers_of_Mf->eight_thirds / pAmp->fcutInsp_eight_thirds * pAmp->rho2
            + powers_of_Mf->three        / pAmp->fcutInsp_three * pAmp->rho3
            ;
        pAmp->InspRescaleFactor = tmp;
        if (pAmp->InspRescaleFactor == 0){
            InspAmp *= powers_of_Mf->m_seven_sixths * pWFHM->ampNorm;
        }
    }
    else{
        /* New release. FIXME: improve how the ansatz is built */
        InspAmp = IMRPhenomXHM_Inspiral_PNAmp_Ansatz(powers_of_Mf, pWFHM, pAmp);
        pseudoterms = powers_of_Mf->seven_thirds / pAmp->fcutInsp_seven_thirds * pAmp->InspiralCoefficient[0]
                + powers_of_Mf->eight_thirds / pAmp->fcutInsp_eight_thirds * pAmp->InspiralCoefficient[1]
                + powers_of_Mf->three        / pAmp->fcutInsp_three * pAmp->InspiralCoefficient[2];
        pseudoterms *= powers_of_Mf->m_seven_sixths * pAmp->PNdominant;
        InspAmp += pseudoterms;
        if (pAmp->InspRescaleFactor != 0){
            InspAmp /= RescaleFactor(powers_of_Mf, pAmp, pAmp->InspRescaleFactor);
        }
    }

    return InspAmp ;
}

/* Numerical derivative of the inspiral (4th order finite differences)
   It is used for reconstructing the intermediate region */
static double IMRPhenomXHM_Inspiral_Amp_NDAnsatz(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMAmpCoefficients *pAmp){

    double df = 10e-10;
    double Nder;
    double fun2R, funR, funL, fun2L;
    double centralfreq = powers_of_Mf->itself;

    IMRPhenomX_UsefulPowers powers_of_Mf2R, powers_of_MfR, powers_of_MfL, powers_of_Mf2L;

    IMRPhenomX_Initialize_Powers(&powers_of_Mf2R, centralfreq + 2*df);
    IMRPhenomX_Initialize_Powers(&powers_of_MfR, centralfreq + df);
    IMRPhenomX_Initialize_Powers(&powers_of_MfL, centralfreq  - df);
    IMRPhenomX_Initialize_Powers(&powers_of_Mf2L, centralfreq - 2*df);

    fun2R = IMRPhenomXHM_Inspiral_Amp_Ansatz(&powers_of_Mf2R, pWFHM, pAmp);
    funR  = IMRPhenomXHM_Inspiral_Amp_Ansatz(&powers_of_MfR, pWFHM, pAmp);
    funL  = IMRPhenomXHM_Inspiral_Amp_Ansatz(&powers_of_MfL, pWFHM, pAmp);
    fun2L = IMRPhenomXHM_Inspiral_Amp_Ansatz(&powers_of_Mf2L, pWFHM, pAmp);
    Nder = (-fun2R + 8*funR - 8*funL + fun2L )/(12*df);

    return Nder;
}

/* VETO function:  when extrapolating the model outside the calibration region the collocation points can be bad behaved.
   With this function we select those that will be used in the reconstruction. */
void IMRPhenomXHM_Inspiral_Amplitude_Veto(
  double *iv1, double *iv2, double *iv3,
  IMRPhenomX_UsefulPowers *powers_of_f1,
  IMRPhenomX_UsefulPowers *powers_of_f2,
  IMRPhenomX_UsefulPowers *powers_of_f3,
  IMRPhenomXHMAmpCoefficients *pAmp,
  IMRPhenomXHMWaveformStruct *pWFHM
)
{
    double threshold = 0.2/(pWFHM->ampNorm);
    #if DEBUG == 1
    printf("\n\nf1, f2, f3 = %.16f %.16f %.16f\n\n",powers_of_f1->itself,powers_of_f2->itself,powers_of_f3->itself);
    printf("\n\nf1, f2, f3 = %.16f %.16f %.16f\n\n",threshold*powers_of_f1->seven_sixths,threshold*powers_of_f2->seven_sixths,threshold*powers_of_f3->seven_sixths);
    printf("\nInspiral Veto: AmpVersion = %i",pWFHM->IMRPhenomXHMInspiralAmpVersion);
    #endif
    // Remove too low collocation points (heuristic).
    if(pAmp->CollocationPointsValuesAmplitudeInsp[0] < threshold*powers_of_f1->seven_sixths){
        *iv1 = 0;
        pWFHM->IMRPhenomXHMInspiralAmpVersion = 2;
    }
    if(pAmp->CollocationPointsValuesAmplitudeInsp[1] < threshold*powers_of_f2->seven_sixths){
        *iv2 = 0;
        pWFHM->IMRPhenomXHMInspiralAmpVersion = pWFHM->IMRPhenomXHMInspiralAmpVersion -1;
    }
    if(pAmp->CollocationPointsValuesAmplitudeInsp[2] < threshold*powers_of_f3->seven_sixths){
        *iv3 = 0;
        pWFHM->IMRPhenomXHMInspiralAmpVersion = pWFHM->IMRPhenomXHMInspiralAmpVersion -1;
    }
}

// Check if the three collocation points are wavy
int WavyPoints(double p1, double p2, double p3){
    if((p1>p2 && p2<p3) || (p1<p2 && p2>p3)){
        return 1;
    }else{
        return 0;
    }
}

void IMRPhenomXHM_Get_Inspiral_Amp_Coefficients(IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22){

    // Get CollocationPointsFreqsAmplitudeInsp and CollocationPointsValuesAmplitudeInsp
    IMRPhenomX_UsefulPowers *powers_of_Mf_inspcollpoints = (IMRPhenomX_UsefulPowers *)XLALMalloc((pWFHM->nCollocPtsInspAmp + 1) * sizeof(IMRPhenomX_UsefulPowers));
    IMRPhenomXHM_Inspiral_Amp_CollocationPoints(pAmp, pWFHM, pWF22);

    for(UINT2 i = 0; i < pWFHM->nCollocPtsInspAmp; i++){
        int status = IMRPhenomX_Initialize_Powers(&(powers_of_Mf_inspcollpoints[i]),  pAmp->CollocationPointsFreqsAmplitudeInsp[i]);
        if(status != XLAL_SUCCESS)
            XLALPrintError("IMRPhenomXHM_Get_Inspiral_Amp_Coefficients failed for Mf, initial_status=%d",status);
    }
    int status = IMRPhenomX_Initialize_Powers(&(powers_of_Mf_inspcollpoints[pWFHM->nCollocPtsInspAmp]),  pAmp->fAmpMatchIN);
    if(status != XLAL_SUCCESS)
        XLALPrintError("IMRPhenomXHM_Get_Inspiral_Amp_Coefficients failed for Mf, initial_status=%d",status);
    pAmp->fcutInsp_seven_thirds = powers_of_Mf_inspcollpoints[pWFHM->nCollocPtsInspAmp].seven_thirds;
    pAmp->fcutInsp_eight_thirds = powers_of_Mf_inspcollpoints[pWFHM->nCollocPtsInspAmp].eight_thirds;
    pAmp->fcutInsp_three = powers_of_Mf_inspcollpoints[pWFHM->nCollocPtsInspAmp].three;


    // Do Vetos? set final pWFHM->IMRPhenomXHMInspiralAmpVersion

    // Get PN values at collocation points frequencies
    for (INT4 i=0; i < pWFHM->nCollocPtsInspAmp; i++){
        pAmp->PNAmplitudeInsp[i] = IMRPhenomXHM_Inspiral_PNAmp_Ansatz(&(powers_of_Mf_inspcollpoints[i]), pWFHM, pAmp);
    }

    IMRPhenomXHM_Inspiral_Amp_Coefficients(pAmp, powers_of_Mf_inspcollpoints, pWFHM);

    LALFree(powers_of_Mf_inspcollpoints);
}

static void IMRPhenomXHM_Inspiral_Amp_CollocationPoints(IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXWaveformStruct *pWF22){
    switch(pWFHM->IMRPhenomXHMInspiralAmpFreqsVersion){
        case 122022:
        case 102021:{
            pAmp->CollocationPointsFreqsAmplitudeInsp[0] = 0.5  * pAmp->fAmpMatchIN;
            pAmp->CollocationPointsFreqsAmplitudeInsp[1] = 0.75 * pAmp->fAmpMatchIN;
            pAmp->CollocationPointsFreqsAmplitudeInsp[2] = pAmp->fAmpMatchIN;
            break;
        }
        // FIXME: Add cases for equispaced, Chebyshev
        default: {XLAL_ERROR_VOID(XLAL_EDOM, "Error in IMRPhenomXHM_Inspiral_CollocationPoints: IMRPhenomXHMInspiralAmpFreqsVersion = %i is not valid. Recommneded version is 102021.\n", pWFHM->IMRPhenomXHMInspiralAmpFreqsVersion);}
    }
    for(UINT2 i = 0; i < pWFHM->nCollocPtsInspAmp; i++){
        pAmp->CollocationPointsValuesAmplitudeInsp[i] = fabs(pAmp->InspiralAmpFits[pWFHM->modeInt * pWFHM->nCollocPtsInspAmp + i](pWF22, pWFHM->IMRPhenomXHMInspiralAmpFitsVersion));
    }

    if (pWFHM->InspiralAmpVeto == 1){
        pWFHM->IMRPhenomXHMInspiralAmpVersion = 13;
    }
}

static void IMRPhenomXHM_Inspiral_Amp_Coefficients(IMRPhenomXHMAmpCoefficients *pAmp, IMRPhenomX_UsefulPowers *powers_of_Mf_inspcollpoints, IMRPhenomXHMWaveformStruct *pWFHM){

    for (UINT2 i = 0; i < N_MAX_COEFFICIENTS_AMPLITUDE_INS; i++){
        pAmp->InspiralCoefficient[i] = 0;
    }

    IMRPhenomX_UsefulPowers *powers_of_f1 = &(powers_of_Mf_inspcollpoints[0]);
    IMRPhenomX_UsefulPowers *powers_of_f2 = &(powers_of_Mf_inspcollpoints[1]);
    IMRPhenomX_UsefulPowers *powers_of_f3 = &(powers_of_Mf_inspcollpoints[2]);
    IMRPhenomX_UsefulPowers *powers_of_finsp = &(powers_of_Mf_inspcollpoints[3]);

    REAL8 v1 = (pAmp->CollocationPointsValuesAmplitudeInsp[0] - pAmp->PNAmplitudeInsp[0]) / pAmp->PNdominant / powers_of_f1->m_seven_sixths;
    REAL8 v2 = (pAmp->CollocationPointsValuesAmplitudeInsp[1] - pAmp->PNAmplitudeInsp[1]) / pAmp->PNdominant / powers_of_f2->m_seven_sixths;
    REAL8 v3 = (pAmp->CollocationPointsValuesAmplitudeInsp[2] - pAmp->PNAmplitudeInsp[2]) / pAmp->PNdominant / powers_of_f3->m_seven_sixths;

    switch(pWFHM->IMRPhenomXHMInspiralAmpVersion){
        case 1:{
            pAmp->InspiralCoefficient[0] = (powers_of_finsp->seven_thirds*v1)/powers_of_f1->seven_thirds;
            break;
        }
        case 2:{
            pAmp->InspiralCoefficient[0] = (powers_of_finsp->seven_thirds*v2)/powers_of_f2->seven_thirds;
            break;
        }
        case 3:{
            pAmp->InspiralCoefficient[0] = (powers_of_finsp->seven_thirds*v3)/powers_of_f3->seven_thirds;
            break;
        }
        case 12:{
            pAmp->InspiralCoefficient[0] = (powers_of_finsp->seven_thirds*(-(powers_of_f2->eight_thirds*v1) + powers_of_f1->eight_thirds*v2))/(powers_of_f1->seven_thirds*(powers_of_f1->one_third - powers_of_f2->one_third)*powers_of_f2->seven_thirds);
            pAmp->InspiralCoefficient[1] = (powers_of_finsp->eight_thirds*(powers_of_f1->m_seven_thirds*v1 - powers_of_f2->m_seven_thirds*v2))/(powers_of_f1->one_third - powers_of_f2->one_third);
            break;
        }
        case 13:{
            pAmp->InspiralCoefficient[0] = (powers_of_finsp->seven_thirds*(-(powers_of_f3->eight_thirds*v1) + powers_of_f1->eight_thirds*v3))/(powers_of_f1->seven_thirds*(powers_of_f1->one_third - powers_of_f3->one_third)*powers_of_f3->seven_thirds);
            pAmp->InspiralCoefficient[1] = (powers_of_finsp->eight_thirds*(powers_of_f1->m_seven_thirds*v1 - powers_of_f3->m_seven_thirds*v3))/(powers_of_f1->one_third - powers_of_f3->one_third);
            break;
        }
        case 23:{
            pAmp->InspiralCoefficient[0] = (powers_of_finsp->seven_thirds*(-(powers_of_f3->eight_thirds*v2) + powers_of_f2->eight_thirds*v3))/(powers_of_f2->seven_thirds*(powers_of_f2->one_third - powers_of_f3->one_third)*powers_of_f3->seven_thirds);
            pAmp->InspiralCoefficient[1] = (powers_of_finsp->eight_thirds*(powers_of_f1->m_seven_thirds*v1 - powers_of_f3->m_seven_thirds*v3))/(powers_of_f1->one_third - powers_of_f3->one_third);
            break;
        }
        case 123:{
            pAmp->InspiralCoefficient[0] = (powers_of_finsp->seven_thirds*(-(powers_of_f1->three*powers_of_f3->eight_thirds*v2) + powers_of_f1->eight_thirds*powers_of_f3->three*v2 + powers_of_f2->three*(powers_of_f3->eight_thirds*v1 - powers_of_f1->eight_thirds*v3) + powers_of_f2->eight_thirds*(-(powers_of_f3->three*v1) + powers_of_f1->three*v3)))/(powers_of_f1->seven_thirds*(powers_of_f1->one_third - powers_of_f2->one_third)*powers_of_f2->seven_thirds*(powers_of_f1->one_third - powers_of_f3->one_third)*(powers_of_f2->one_third - powers_of_f3->one_third)*powers_of_f3->seven_thirds);
            pAmp->InspiralCoefficient[1] = (powers_of_finsp->eight_thirds*(powers_of_f1->three*powers_of_f3->seven_thirds*v2 - powers_of_f1->seven_thirds*powers_of_f3->three*v2 + powers_of_f2->three*(-(powers_of_f3->seven_thirds*v1) + powers_of_f1->seven_thirds*v3) + powers_of_f2->seven_thirds*(powers_of_f3->three*v1 - powers_of_f1->three*v3)))/(powers_of_f1->seven_thirds*(powers_of_f1->one_third - powers_of_f2->one_third)*powers_of_f2->seven_thirds*(powers_of_f1->one_third - powers_of_f3->one_third)*(powers_of_f2->one_third - powers_of_f3->one_third)*powers_of_f3->seven_thirds);
            pAmp->InspiralCoefficient[2] = (powers_of_finsp->three*(powers_of_f1->seven_thirds*(-powers_of_f1->one_third + powers_of_f3->one_third)*powers_of_f3->seven_thirds*v2 + powers_of_f2->seven_thirds*(-(powers_of_f3->eight_thirds*v1) + powers_of_f1->eight_thirds*v3) + powers_of_f2->eight_thirds*(powers_of_f3->seven_thirds*v1 - powers_of_f1->seven_thirds*v3)))/(powers_of_f1->seven_thirds*(powers_of_f1->one_third - powers_of_f2->one_third)*powers_of_f2->seven_thirds*(powers_of_f1->one_third - powers_of_f3->one_third)*(powers_of_f2->one_third - powers_of_f3->one_third)*powers_of_f3->seven_thirds);
            break;
        }
    }
}
/*************************************/
/*                                   */
/*              PHASE                */
/*                                   */
/*************************************/

// Spin parameter S = (m1^2*chi1 + m2^2*chi2)/(m1^2 + m2^2)

// Below we give paramater-space fits for the weighted difference between each mode's phase and the 22-phase: phi_lm-m/2 phi_22(2/m f), see Eqs. (4.10-4.12)

/* Start of Phase Parameter Space Fits */

static double IMRPhenomXHM_Insp_Phase_21_lambda(IMRPhenomXWaveformStruct *pWF, int InspPhaseFlag) {
    double total;
    switch (InspPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,eta5,S2;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            S2 = pow(S,2);
            double noSpin = 13.664473636545068 - 170.08866400251395*eta + 3535.657736681598*eta2 - 26847.690494515424*eta3 + 96463.68163125668*eta4 - 133820.89317471132*eta5;
            double eqSpin = (S*(18.52571430563905 - 41.55066592130464*S + eta3*(83493.24265292779 + 16501.749243703132*S - 149700.4915210766*S2) + eta*(3642.5891077598003 + 1198.4163078715173*S - 6961.484805326852*S2) + 33.8697137964237*S2 + eta2*(-35031.361998480075 - 7233.191207000735*S + 62149.00902591944*S2)))/(6.880288191574696 + 1.*S);
            double uneqSpin = -134.27742343186577*pWF->dchi*sqrt(1.-4.*eta)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR(XLAL_EDOM, "Error in IMRPhenomXHM_Insp_Phase_21_lambda: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Phase_33_lambda(IMRPhenomXWaveformStruct *pWF, int InspPhaseFlag) {
    double total=0;
    switch (InspPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double delta=sqrt(1.-4.*eta);
            double eta2,eta3;
            eta2 = pow(eta,2);
            eta3 = eta2*eta;
            double noSpin = 4.1138398568400705 + 9.772510519809892*eta - 103.92956504520747*eta2 + 242.3428625556764*eta3;
            double eqSpin = ((-0.13253553909611435 + 26.644159828590055*eta - 105.09339163109497*eta2)*S)/(1. + 0.11322426762297967*S);
            double uneqSpin = -19.705359163581168*pWF->dchi*eta2*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR(XLAL_EDOM, "Error in IMRPhenomXHM_Insp_Phase_32_lambda: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

static double IMRPhenomXHM_Insp_Phase_32_lambda(IMRPhenomXWaveformStruct *pWF, int InspPhaseFlag) {
    double total;
    switch (InspPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,S2,S3,S4;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4 = pow(S,4);
            double noSpin = (9.913819875501506 + 18.424900617803107*eta - 574.8672384388947*eta2 + 2671.7813055097877*eta3 - 6244.001932443913*eta4)/(1. - 0.9103118343073325*eta);
            double eqSpin = (-4.367632806613781 + 245.06757304950986*eta - 2233.9319708029775*eta2 + 5894.355429022858*eta3)*S + (-1.375112297530783 - 1876.760129419146*eta + 17608.172965575013*eta2 - 40928.07304790013*eta3)*S2 + (-1.28324755577382 - 138.36970336658558*eta + 708.1455154504333*eta2 - 273.23750933544176*eta3)*S3 + (1.8403161863444328 + 2009.7361967331492*eta - 18636.271414571278*eta2 + 42379.205045791656*eta3)*S4;
            double uneqSpin = pWF->dchi*sqrt(1.-4.*eta)*eta2*(-105.34550407768225 - 1566.1242344157668*pWF->chi1L*eta + 1566.1242344157668*pWF->chi2L*eta + 2155.472229664981*eta*S);
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR(XLAL_EDOM, "Error in IMRPhenomXHM_Insp_Phase_32_lambda: version is not valid. Recommended version is 122019.");}
    }
    return total;
}


static double IMRPhenomXHM_Insp_Phase_44_lambda(IMRPhenomXWaveformStruct *pWF, int InspPhaseFlag) {
    double total;
    switch (InspPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,eta5,S2;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            S2 = pow(S,2);
            double noSpin = 5.254484747463392 - 21.277760168559862*eta + 160.43721442910618*eta2 - 1162.954360723399*eta3 + 1685.5912722190276*eta4 - 1538.6661348106031*eta5;
            double eqSpin = (0.007067861615983771 - 10.945895160727437*eta + 246.8787141453734*eta2 - 810.7773268493444*eta3)*S + (0.17447830920234977 + 4.530539154777984*eta - 176.4987316167203*eta2 + 621.6920322846844*eta3)*S2;
            double uneqSpin = -8.384066369867833*pWF->dchi*sqrt(1.-4.*eta)*eta2;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR(XLAL_EDOM, "Error in IMRPhenomXHM_Insp_Phase_44_lambda: version is not valid. Recommended version is 122019.");}
    }
    return total;
}

/* End of Phase Parameter Space Fits */

// Returns linear-in-f term coming from the complex part of each mode's amplitude that will be added to the orbital phase, this is an expansion of Eq. (4.9) truncated at linear order
static double IMRPhenomXHM_Insp_Phase_LambdaPN(double eta, int modeInt){

    double output;
    const double PI = powers_of_lalpiHM.itself;

    switch(modeInt){
        case 21:
        {//2 f \[Pi] (-(1/2) - Log[16]/2)
            output=(2.*PI*(-0.5-2.*log(2)));
            break;
        }
        case 33:
        {//2/3 f \[Pi] (-(21/5) + 6 Log[3/2])
            output=2./3.*PI*(-21./5.+6.*log(1.5));
            break;
        }
        case 32:
        { //-((2376 f \[Pi] (-5 + 22 \[Eta]))/(-3960 + 11880 \[Eta]))
            output=-((2376.*PI*(-5.+22.*eta))/(-3960.+11880*eta));
            break;
        }
        case 44:
        { //(45045 f \[powers_of_lalpiHM.itself] (336 - 1193 \[Eta] +320 (-1 + 3 \[Eta]) Log[2]))/(2 (-1801800 + 5405400 \[Eta]))
            output=45045.*PI*(336.-1193.*eta+320.*(-1.+3.*eta)*log(2))/(2.*(-1801800. + 5405400.*eta));
            break;
        }
        default: output=0.;
    }

    return -1.*output;
}


/************* INSPIRAL PHASE ANSATZ *************/

/*
 The function loads the rescaled phenX coefficients for the phase and correct the result with a linear-in-f term coming from the complex part of the PN amplitude

 The rescaling of the 22-coefficients is explained in App. D
*/
static double IMRPhenomXHM_Inspiral_Phase_AnsatzInt(double Mf, IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMPhaseCoefficients *pPhase)
{
    //compute the orbital phase by laoding the rescaled phenX coefficients
    double philm=0.;
    double freqs[N_MAX_COEFFICIENTS_PHASE_INS]={powers_of_Mf->m_five_thirds,powers_of_Mf->m_four_thirds,powers_of_Mf->m_one,powers_of_Mf->m_two_thirds,powers_of_Mf->m_one_third,1.,powers_of_Mf->one_third,powers_of_Mf->two_thirds, Mf, powers_of_Mf->four_thirds, powers_of_Mf->five_thirds, powers_of_Mf->two,powers_of_Mf->seven_thirds};
    double logMf=powers_of_Mf->log;
    for(int i=0; i<N_MAX_COEFFICIENTS_PHASE_INS; i++)
        philm+=(pPhase->phi[i]+pPhase->phiL[i]*(logMf))*freqs[i];
    return philm;
}

static double IMRPhenomXHM_Inspiral_Phase_Ansatz(double Mf, IMRPhenomX_UsefulPowers *powers_of_Mf,IMRPhenomXHMPhaseCoefficients *pPhase)
{
    //compute the orbital phase by laoding the rescaled phenX coefficients
    double dphilm=0.;
    double coeffs[N_MAX_COEFFICIENTS_PHASE_INS]={-5./3,-4./3,-1.,-2./3,-1./3,0.,1./3, 2./3, 1., 4./3, 5./3, 2., 7./3};
    double freqs[N_MAX_COEFFICIENTS_PHASE_INS]={powers_of_Mf->m_eight_thirds,powers_of_Mf->m_seven_thirds,powers_of_Mf->m_two,powers_of_Mf->m_five_thirds,powers_of_Mf->m_four_thirds,powers_of_Mf->m_one,powers_of_Mf->m_two_thirds,powers_of_Mf->m_one_third,1.,powers_of_Mf->one_third, powers_of_Mf->two_thirds, Mf,powers_of_Mf->four_thirds};
    double logMf=powers_of_Mf->log;

    for(int i=0; i<N_MAX_COEFFICIENTS_PHASE_INS; i++)
        dphilm+=((pPhase->phi[i]+pPhase->phiL[i]*(logMf))*coeffs[i]+pPhase->phiL[i])*freqs[i];
    return dphilm;
}
