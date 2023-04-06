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

/* Start of Amp Parameter Space Fits */


static double IMRPhenomXHM_RD_Amp_21_alambda(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag) {
    double total=0;
    switch (RDAmpFlag){
        case 122018:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double delta=sqrt(1.-4.*eta),eta2,S2;
            eta2 = pow(eta,2);
            S2 = pow(S,2);
            double noSpin = sqrt(eta - 4.*eta2)*(0.00734983387668636 - 0.0012619735607202085*eta + 0.01042318959002753*eta2);
            double eqSpin = sqrt(eta - 4.*eta2)*S*(-0.004839645742570202 - 0.0013927779195756036*S + eta2*(-0.054621206928483663 + 0.025956604949552205*S + 0.020360826886107204*S2));
            double uneqSpin = -0.018115657394753674*pWF->dchi*eta2*(1. - 10.539795474715346*eta2*delta);
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
            double S1 = S;
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            total = fabs(delta*(0.24548180919287976 - 0.25565119457386487*eta1)*eta1 + chidiff1*delta*eta1*(0.5670798742968471*eta1 - 14.276514548218454*eta2 + 45.014547333879136*eta3) + chidiff1*delta*eta1*(0.4580805242442763*eta1 - 4.859294663135058*eta2 + 14.995447609839573*eta3)*S1 + chidiff1*eta5*(-27.031582936285528 + 6.468164760468401*S1 + 0.34222101136488015*S2) + delta*eta1*S1*(-0.2204878224611389*(1.0730799832007898 - 3.44643820338605*eta1 + 32.520429274459836*eta2 - 83.21097158567372*eta3) + 0.008901444811471891*(-5.876973170072921 + 120.70115519895002*eta1 - 916.5281661566283*eta2 + 2306.8425350489847*eta3)*S1 + 0.015541783867953005*(2.4780170686140455 + 17.377013149762398*eta1 - 380.91157168170236*eta2 + 1227.5332509075172*eta3)*S2));
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_21_alambda: version %i is not valid.", RDAmpFlag);}
    }
    return total;
}

static double IMRPhenomXHM_RD_Amp_21_lambda(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag) {
    double total=0;
    switch (RDAmpFlag){
        case 122018:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double delta=sqrt(1.-4.*eta), eta2;
            eta2 = pow(eta,2);
            double noSpin = 0.5566284518926176 + 0.12651770333481904*eta + 1.8084545267208734*eta2;
            double eqSpin = (0.29074922226651545 + eta2*(-2.101111399437034 - 3.4969956644617946*S) + eta*(0.059317243606471406 - 0.31924748117518226*S) + 0.27420263462336675*S)*S;
            double uneqSpin = 1.0122975748481835*pWF->dchi*eta2*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        case 122022:{
            double eta = pWF->eta;
            double S = pWF->chiPNHat;
            double chidiff = pWF->dchi_half;
            double eta1 = eta;
            double eta2 = eta1 * eta1;
            double eta3 = eta1 * eta2;
            double eta4 = eta1 * eta3;
            double S1 = S;
            double chidiff1 = chidiff;
            total = fabs(1.0092933052569775 - 0.2791855444800297*eta1 + 1.7110615047319937*eta2 + chidiff1*(-0.1054835719277311*eta1 + 7.506083919925026*eta2 - 30.1595680078279*eta3) + chidiff1*(2.078267611384239*eta1 - 10.166026002515457*eta2 - 1.2091616330625208*eta3)*S1 + S1*(0.17250873250247642*(1.0170226856985174 + 1.0395650952176598*eta1 - 35.73623734051525*eta2 + 403.68074286921444*eta3 - 1194.6152711219886*eta4) + 0.06850746964805364*(1.507796537056924 + 37.81075363806507*eta1 - 863.117144661059*eta2 + 6429.543634627373*eta3 - 15108.557419182316*eta4)*S1));
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_21_lambda: version %i is not valid.", RDAmpFlag);}
    }
    return total;
}

static double IMRPhenomXHM_RD_Amp_33_alambda(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag) {
  double total=0;
  switch (RDAmpFlag){
    case 122018:{
      double eta = pWF->eta;
      double S = pWF->STotR;
      double eta2,eta4,delta=sqrt(1-4*eta);
      eta2 = pow(eta,2);
      eta4 = pow(eta,4);
      double noSpin = sqrt(eta - 4.*eta2)*(0.013700854227665184 + 0.01202732427321774*eta + 0.0898095508889557*eta2);
      double eqSpin = sqrt(eta - 4.*eta2)*(0.0075858980586079065 + eta*(-0.013132320758494439 - 0.018186317026076343*S) + 0.0035617441651710473*S)*S;
      double uneqSpin = eta4*(pWF->chi2L*(-0.09802218411554885 - 0.05745949361626237*S) + pWF->chi1L*(0.09802218411554885 + 0.05745949361626237*S) + eta2*(pWF->chi1L*(-4.2679864481479886 - 11.877399902871485*S) + pWF->chi2L*(4.2679864481479886 + 11.877399902871485*S))*delta);
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_33_alambda: version %i is not valid.", RDAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_33_lambda(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag) {
  double total=0;
  switch (RDAmpFlag){
    case 122018:{
      double eta = pWF->eta;
      double S = pWF->STotR;
      double eta2,S2,delta=sqrt(1-4*eta);
      eta2 = pow(eta,2);
      S2 = pow(S,2);
      double noSpin = 0.7435306475478924 - 0.06688558533374556*eta + 1.471989765837694*eta2;
      double eqSpin = S*(0.19457194111990656 + 0.07564220573555203*S + eta*(-0.4809350398289311 + 0.17261430318577403*S - 0.1988991467974821*S2));
      double uneqSpin = 1.8881959341735146*pWF->dchi*eta2*delta;
      total = noSpin + eqSpin + uneqSpin;
      break;
    }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_33_lambda: version %i is not valid.", RDAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_32_alambda(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag) {
  double total=0;
  switch (RDAmpFlag){
      case 122018:{
        double eta = pWF->eta;
        double S = pWF->STotR;
        double eta2,eta3,eta4,eta5,eta6;
        eta2 = pow(eta,2);
        eta3 = pow(eta,3);
        eta4 = pow(eta,4);
        eta5 = pow(eta,5);
        eta6 = pow(eta,6);
        double noSpin = 0.00012587900257140724 + 0.03927886286971654*eta - 0.8109309606583066*eta2 + 8.820604907164254*eta3 - 51.43344812454074*eta4 + 141.81940900657446*eta5 - 140.0426973304466*eta6;
        double eqSpin = S*(-0.00006001471234796344 + eta4*(-0.7849112300598181 - 2.09188976953315*S) + eta2*(0.08311497969032984 - 0.15569578955822236*S) + eta*(-0.01083175709906557 + 0.00568899459837252*S) - 0.00009363591928190229*S + 1.0670798489407887*eta3*S);
        double uneqSpin = -0.04537308968659669*pow(pWF->dchi,2)*eta2*(1. - 8.711096029480697*eta + 18.362371966229926*eta2) + pWF->dchi*(-297.36978685672733 + 3103.2516759087644*eta - 10001.774055779177*eta2 + 9386.734883473799*eta3)*eta6;
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
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = chidiff2*(-3.4614418482110163*eta3 + 35.464117772624164*eta4 - 85.19723511005235*eta5) + chidiff1*delta*(2.0328561081997463*eta3 - 46.18751757691501*eta4 + 170.9266105597438*eta5) + chidiff2*(-0.4600401291210382*eta3 + 12.23450117663151*eta4 - 42.74689906831975*eta5)*S1 + chidiff1*delta*(5.786292428422767*eta3 - 53.60467819078566*eta4 + 117.66195692191727*eta5)*S1 + S1*(-0.0013330716557843666*(56.35538385647113*eta1 - 1218.1550992423377*eta2 + 16509.69605686402*eta3 - 102969.88022112886*eta4 + 252228.94931931415*eta5 - 150504.2927996263*eta6) + 0.0010126460331462495*(-33.87083889060834*eta1 + 502.6221651850776*eta2 - 1304.9210590188136*eta3 - 36980.079328277505*eta4 + 295469.28617550555*eta5 - 597155.7619486618*eta6)*S1 - 0.00043088431510840695*(-30.014415072587354*eta1 - 1900.5495690280086*eta2 + 76517.21042363928*eta3 - 870035.1394696251*eta4 + 3.9072674134789007e6*eta5 - 6.094089675611567e6*eta6)*S2) + (0.08408469319155859*eta1 - 1.223794846617597*eta2 + 6.5972460654253515*eta3 - 15.707327897569396*eta4 + 14.163264397061505*eta5)*pow(1 - 8.612447115134758*eta1 + 18.93655612952139*eta2,-1);
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_32_alambda: version %i is not valid.", RDAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_32_lambda(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag) {
  double total=0;
  switch (RDAmpFlag){
    case 122018:{
      double eta = pWF->eta;
      double S = pWF->STotR;
      double eta2,eta3,delta=sqrt(1.-4*eta);
      eta2 = pow(eta,2);
      eta3 = pow(eta,3);
      double noSpin = (sqrt(1. - 3.*eta)*(0.0341611244787871 - 0.3197209728114808*eta + 0.7689553234961991*eta2))/(0.048429644168112324 - 0.43758296068790314*eta + eta2);
      double eqSpin = sqrt(1. - 3.*eta)*S*(0.11057199932233873 + eta2*(25.536336676250748 - 71.18182757443142*S) + 9.790509295728649*eta*S + eta3*(-56.96407763839491 + 175.47259563543165*S));
      double uneqSpin = -5.002106168893265*pow(pWF->dchi,2)*eta2*delta;
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
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = 0.978510781593996 + 0.36457571743142897*eta1 - 12.259851752618998*eta2 + 49.19719473681921*eta3 + chidiff1*delta*(-188.37119473865533*eta3 + 2151.8731700399308*eta4 - 6328.182823770599*eta5) + chidiff2*(115.3689949926392*eta3 - 1159.8596972989067*eta4 + 2657.6998831179444*eta5) + S1*(0.22358643406992756*(0.48943645614341924 - 32.06682257944444*eta1 + 365.2485484044132*eta2 - 915.2489655397206*eta3) + 0.0792473022309144*(1.877251717679991 - 103.65639889587327*eta1 + 1202.174780792418*eta2 - 3206.340850767219*eta3)*S1);
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_32_lambda: version %i is not valid.", RDAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_44_alambda(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag) {
    double total=0;
    switch (RDAmpFlag){
        case 122018:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double delta=sqrt(1.-4.*eta),eta2,eta3,eta4,eta5,S2;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            S2 = pow(S,2);
            double noSpin = sqrt(eta - 3.*eta2)*(0.007904587819112173 + 0.09558474985614368*eta - 2.663803397359775*eta2 + 28.298192768381554*eta3 - 136.10446022757958*eta4 + 233.23167528016833*eta5);
            double eqSpin = sqrt(eta - 3.*eta2)*S*(0.0049703757209330025 + 0.004122811292229324*S + eta*(-0.06166686913913691 + 0.014107365722576927*S)*S + eta2*(-0.2945455034809188 + 0.4139026619690879*S - 0.1389170612199015*S2) + eta3*(0.9225758392294605 - 0.9656098473922222*S + 0.19708289555425246*S2) + 0.000657528128497184*S2);
            double uneqSpin = 0.00659873279539475*pWF->dchi*eta2*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_44_a: version %i is not valid.", RDAmpFlag);}
    }
    return total;
}

static double IMRPhenomXHM_RD_Amp_44_lambda(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag) {
    double total=0;
    switch (RDAmpFlag){
        case 122018:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double delta=sqrt(1.-4.*eta),eta2,eta3,eta4,eta5,eta6,eta7;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            eta6 = pow(eta,6);
            eta7 = pow(eta,7);
            double noSpin = 0.7702864948772887 + 32.81532373698395*eta - 1625.1795901450212*eta2 + 31305.458876573215*eta3 - 297375.5347399236*eta4 + 1.4726521941846698e6*eta5 - 3.616582470072637e6*eta6 + 3.4585865843680725e6*eta7;
            double eqSpin = (-0.03011582009308575*S + 0.09034746027925727*eta*S + 1.8738784391649446*eta2*S - 5.621635317494836*eta3*S)/(-1.1340218677260014 + S);
            double uneqSpin = 0.959943270591552*pow(pWF->dchi,2)*eta2 + 0.853573071529436*pWF->dchi*eta2*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_44_lambda: version %i is not valid.", RDAmpFlag);}
    }
    return total;
}

static double IMRPhenomXHM_RD_Amp_21_sigma(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag) {
    double total=0;
    switch (RDAmpFlag){
        case 122018:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double delta=sqrt(1.-4.*eta),eta2;
            eta2 = pow(eta,2);
            double noSpin = 1.2922261617161441 + 0.0019318405961363861*eta;
            double eqSpin = (0.04927982551108649 - 0.6703778360948937*eta + 2.6625014134659772*eta2)*S;
            double uneqSpin = 1.2001101665670462*(pWF->chi1L + pow(pWF->chi1L,2) - 2.*pWF->chi1L*pWF->chi2L + (-1. + pWF->chi2L)*pWF->chi2L)*eta2*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        case 122022:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double chidiff = pWF->dchi_half;
            double eta1 = eta;
            double eta2 = eta1 * eta1;
            double eta3 = eta1 * eta2;
            double S1 = S;
            double chidiff1 = chidiff;
            total = fabs(1.374451177213076 - 0.1147381625630186*eta1 + chidiff1*(0.6646459256372743*eta1 - 5.020585319906719*eta2 + 9.817281653770431*eta3) + chidiff1*(3.8734254747587973*eta1 - 39.880716190740465*eta2 + 99.05511583518896*eta3)*S1 + S1*(0.013272603498067647*(1.809972721953344 - 12.560287006325837*eta1 - 134.597005438578*eta2 + 786.2235720637008*eta3) + 0.006850483944311038*(-6.478737679813189 - 200.29813775611166*eta1 + 2744.3629484255357*eta2 - 7612.096007280672*eta3)*S1));
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_21_sigma: version %i is not valid.", RDAmpFlag);}
    }
    return total;
}

static double IMRPhenomXHM_RD_Amp_33_sigma(UNUSED IMRPhenomXWaveformStruct *pWF, int RDAmpFlag) {
    double total=0;
    switch (RDAmpFlag){
        case 122018:{
            double noSpin = 1.3;
            double eqSpin = 0.;
            double uneqSpin = 0.;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_33_sigma: version %i is not valid.", RDAmpFlag);}
    }
    return total;
}

static double IMRPhenomXHM_RD_Amp_32_sigma(UNUSED IMRPhenomXWaveformStruct *pWF, int RDAmpFlag) {
    double total=0;
    switch (RDAmpFlag){
        case 122018:{
            double noSpin = 1.33;
            double eqSpin = 0.;
            double uneqSpin = 0.;
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
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = 1.3353917551819414 + 0.13401718687342024*eta1 + chidiff1*delta*(144.37065005786636*eta3 - 754.4085447486738*eta4 + 123.86194078913776*eta5) + chidiff2*(209.09202210427972*eta3 - 1769.4658099037918*eta4 + 3592.287297392387*eta5) + S1*(-0.012086025709597246*(-6.230497473791485 + 600.5968613752918*eta1 - 6606.1009717965735*eta2 + 17277.60594350428*eta3) - 0.06066548829900489*(-0.9208054306316676 + 142.0346574366267*eta1 - 1567.249168668069*eta2 + 4119.373703246675*eta3)*S1);
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_32_sigma: version %i is not valid.", RDAmpFlag);}
    }
    return total;
}

static double IMRPhenomXHM_RD_Amp_44_sigma(UNUSED IMRPhenomXWaveformStruct *pWF, int RDAmpFlag) {
    double total=0;
    switch (RDAmpFlag){
        case 122018:{
            double noSpin = 1.33;
            double eqSpin = 0.;
            double uneqSpin = 0.;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_44_sigma: version %i is not valid.", RDAmpFlag);}
    }
    return total;
}


static double IMRPhenomXHM_RD_Amp_21_rdcp1(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag){
	double total=0;
	switch (RDAmpFlag){
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
            total = fabs(delta*eta1*(12.880905080761432 - 23.5291063016996*eta1 + 92.6090002736012*eta2 - 175.16681482428694*eta3) + chidiff1*delta*eta1*(26.89427230731867*eta1 - 710.8871223808559*eta2 + 2255.040486907459*eta3) + chidiff1*delta*eta1*(21.402708785047853*eta1 - 232.07306353130417*eta2 + 591.1097623278739*eta3)*S1 + delta*eta1*S1*(-10.090867481062709*(0.9580746052260011 + 5.388149112485179*eta1 - 107.22993216128548*eta2 + 801.3948756800821*eta3 - 2688.211889175019*eta4 + 3950.7894052628735*eta5 - 1992.9074348833092*eta6) - 0.42972412296628143*(1.9193131231064235 + 139.73149069609775*eta1 - 1616.9974609915555*eta2 - 3176.4950303461164*eta3 + 107980.65459735804*eta4 - 479649.75188253267*eta5 + 658866.0983367155*eta6)*S1) + chidiff1*eta5*(-1512.439342647443 + 175.59081294852444*S1 + 10.13490934572329*S2));
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_21_rdcp1: version %i is not valid.", RDAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_21_rdcp2(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag){
	double total=0;
	switch (RDAmpFlag){
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
            total = fabs(delta*(9.112452928978168 - 7.5304766811877455*eta1)*eta1 + chidiff1*delta*eta1*(16.236533863306132*eta1 - 500.11964987628926*eta2 + 1618.0818430353293*eta3) + chidiff1*delta*eta1*(2.7866868976718226*eta1 - 0.4210629980868266*eta2 - 20.274691328125606*eta3)*S1 + chidiff1*eta5*(-1116.4039232324135 + 245.73200219767514*S1 + 21.159179960295855*S2) + delta*eta1*S1*(-8.236485576091717*(0.8917610178208336 + 5.1501231412520285*eta1 - 87.05136337926156*eta2 + 519.0146702141192*eta3 - 997.6961311502365*eta4) + 0.2836840678615208*(-0.19281297100324718 - 57.65586769647737*eta1 + 586.7942442434971*eta2 - 1882.2040277496196*eta3 + 2330.3534917059906*eta4)*S1 + 0.40226131643223145*(-3.834742668014861 + 190.42214703482531*eta1 - 2885.5110686004946*eta2 + 16087.433824017446*eta3 - 29331.524552164105*eta4)*S2));
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_21_rdcp2: version %i is not valid.", RDAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_21_rdcp3(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag){
	double total=0;
	switch (RDAmpFlag){
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
            total = fabs(delta*(2.920930733198033 - 3.038523690239521*eta1)*eta1 + chidiff1*delta*eta1*(6.3472251472354975*eta1 - 171.23657247338042*eta2 + 544.1978232314333*eta3) + chidiff1*delta*eta1*(1.9701247529688362*eta1 - 2.8616711550845575*eta2 - 0.7347258030219584*eta3)*S1 + chidiff1*eta5*(-334.0969956136684 + 92.91301644484749*S1 - 5.353399481074393*S2) + delta*eta1*S1*(-2.7294297839371824*(1.148166706456899 - 4.384077347340523*eta1 + 36.120093043420326*eta2 - 87.26454353763077*eta3) + 0.23949142867803436*(-0.6931516433988293 + 33.33372867559165*eta1 - 307.3404155231787*eta2 + 862.3123076782916*eta3)*S1 + 0.1930861073906724*(3.7735099269174106 - 19.11543562444476*eta1 - 78.07256429516346*eta2 + 485.67801863289293*eta3)*S2));
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_21_rdcp3: version %i is not valid.", RDAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_33_rdcp1(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag){
	double total=0;
	switch (RDAmpFlag){
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
            total = delta*eta1*(12.439702602599235 - 4.436329538596615*eta1 + 22.780673360839497*eta2) + delta*eta1*(chidiff1*(-41.04442169938298*eta1 + 502.9246970179746*eta2 - 1524.2981907688634*eta3) + chidiff2*(32.23960072974939*eta1 - 365.1526474476759*eta2 + 1020.6734178547847*eta3)) + chidiff1*delta*eta1*(-52.85961155799673*eta1 + 577.6347407795782*eta2 - 1653.496174539196*eta3)*S1 + chidiff1*eta5*(257.33227387984863 - 34.5074027042393*chidiff2 - 21.836905132600755*S1 - 15.81624534976308*S2) + 13.499999999999998*delta*eta1*S1*(-0.13654149379906394*(2.719687834084113 + 29.023992126142304*eta1 - 742.1357702210267*eta2 + 4142.974510926698*eta3 - 6167.08766058184*eta4 - 3591.1757995710486*eta5) - 0.06248535354306988*(6.697567446351289 - 78.23231700361792*eta1 + 444.79350113344543*eta2 - 1907.008984765889*eta3 + 6601.918552659412*eta4 - 10056.98422430965*eta5)*S1)*pow(-3.9329308614837704 + S1,-1);
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_33_rdcp1: version %i is not valid.", RDAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_33_rdcp2(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag){
	double total=0;
	switch (RDAmpFlag){
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
            total = delta*eta1*(8.425057692276933 + 4.543696144846763*eta1) + chidiff1*delta*eta1*(-32.18860840414171*eta1 + 412.07321398189293*eta2 - 1293.422289802462*eta3) + chidiff1*delta*eta1*(-17.18006888428382*eta1 + 190.73514518113845*eta2 - 636.4802385540647*eta3)*S1 + delta*eta1*S1*(0.1206817303851239*(8.667503604073314 - 144.08062755162752*eta1 + 3188.189172446398*eta2 - 35378.156133055556*eta3 + 163644.2192178668*eta4 - 265581.70142471837*eta5) + 0.08028332044013944*(12.632478544060636 - 322.95832000179297*eta1 + 4777.45310151897*eta2 - 35625.58409457366*eta3 + 121293.97832549023*eta4 - 148782.33687815256*eta5)*S1) + chidiff1*eta5*(159.72371180117415 - 29.10412708633528*chidiff2 - 1.873799747678187*S1 + 41.321480132899524*S2);
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_33_rdcp2: version %i is not valid.", RDAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_33_rdcp3(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag){
	double total=0;
	switch (RDAmpFlag){
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
            total = delta*eta1*(2.485784720088995 + 2.321696430921996*eta1) + delta*eta1*(chidiff1*(-10.454376404653859*eta1 + 147.10344302665484*eta2 - 496.1564538739011*eta3) + chidiff2*(-5.9236399792925996*eta1 + 65.86115501723127*eta2 - 197.51205149250532*eta3)) + chidiff1*delta*eta1*(-10.27418232676514*eta1 + 136.5150165348149*eta2 - 473.30988537734174*eta3)*S1 + chidiff1*eta5*(32.07819766300362 - 3.071422453072518*chidiff2 + 35.09131921815571*S1 + 67.23189816732847*S2) + 13.499999999999998*delta*eta1*S1*(0.0011484326782460882*(4.1815722950796035 - 172.58816646768219*eta1 + 5709.239330076732*eta2 - 67368.27397765424*eta3 + 316864.0589150127*eta4 - 517034.11171277676*eta5) - 0.009496797093329243*(0.9233282181397624 - 118.35865186626413*eta1 + 2628.6024206791726*eta2 - 23464.64953722729*eta3 + 94309.57566199072*eta4 - 140089.40725211444*eta5)*S1)*pow(0.09549360183532198 - 0.41099904730526465*S1 + S2,-1);
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_33_rdcp3: version %i is not valid.", RDAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_44_rdcp1(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag){
	double total=0;
	switch (RDAmpFlag){
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
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = eta1*(chidiff1*delta*(-8.51952446214978*eta1 + 117.76530248141987*eta2 - 297.2592736781142*eta3) + chidiff2*(-0.2750098647982238*eta1 + 4.456900599347149*eta2 - 8.017569928870929*eta3)) + eta1*(5.635069974807398 - 33.67252878543393*eta1 + 287.9418482197136*eta2 - 3514.3385364216438*eta3 + 25108.811524802128*eta4 - 98374.18361532023*eta5 + 158292.58792484726*eta6) + eta1*S1*(-0.4360849737360132*(-0.9543114627170375 - 58.70494649755802*eta1 + 1729.1839588870455*eta2 - 16718.425586396803*eta3 + 71236.86532610047*eta4 - 111910.71267453219*eta5) - 0.024861802943501172*(-52.25045490410733 + 1585.462602954658*eta1 - 15866.093368857853*eta2 + 35332.328181283*eta3 + 168937.32229060197*eta4 - 581776.5303770923*eta5)*S1 + 0.005856387555754387*(186.39698091707513 - 9560.410655118145*eta1 + 156431.3764198244*eta2 - 1.0461268207440731e6*eta3 + 3.054333578686424e6*eta4 - 3.2369858387064277e6*eta5)*S2);
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_44_rdcp1: version %i is not valid.", RDAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_44_rdcp2(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag){
	double total=0;
	switch (RDAmpFlag){
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
            total = eta1*(chidiff1*delta*(-2.861653255976984*eta1 + 50.50227103211222*eta2 - 123.94152825700999*eta3) + chidiff2*(2.9415751419018865*eta1 - 28.79779545444817*eta2 + 72.40230240887851*eta3)) + eta1*(3.2461722686239307 + 25.15310593958783*eta1 - 792.0167314124681*eta2 + 7168.843978909433*eta3 - 30595.4993786313*eta4 + 49148.57065911245*eta5) + eta1*S1*(-0.23311779185707152*(-1.0795711755430002 - 20.12558747513885*eta1 + 1163.9107546486134*eta2 - 14672.23221502075*eta3 + 73397.72190288734*eta4 - 127148.27131388368*eta5) + 0.025805905356653*(11.929946153728276 + 350.93274421955806*eta1 - 14580.02701600596*eta2 + 174164.91607515427*eta3 - 819148.9390278616*eta4 + 1.3238624538095295e6*eta5)*S1 + 0.019740635678180102*(-7.046295936301379 + 1535.781942095697*eta1 - 27212.67022616794*eta2 + 201981.0743810629*eta3 - 696891.1349708183*eta4 + 910729.0219043035*eta5)*S2);
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_44_rdcp2: version %i is not valid.", RDAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_44_rdcp3(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag){
	double total=0;
	switch (RDAmpFlag){
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
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = eta1*(chidiff1*delta*(2.4286414692113816*eta1 - 23.213332913737403*eta2 + 66.58241012629095*eta3) + chidiff2*(3.085167288859442*eta1 - 31.60440418701438*eta2 + 78.49621016381445*eta3)) + eta1*(0.861883217178703 + 13.695204704208976*eta1 - 337.70598252897696*eta2 + 2932.3415281149432*eta3 - 12028.786386004691*eta4 + 18536.937955014455*eta5) + eta1*S1*(-0.048465588779596405*(-0.34041762314288154 - 81.33156665674845*eta1 + 1744.329802302927*eta2 - 16522.343895064576*eta3 + 76620.18243090731*eta4 - 133340.93723954144*eta5) + 0.024804027856323612*(-8.666095805675418 + 711.8727878341302*eta1 - 13644.988225595187*eta2 + 112832.04975245205*eta3 - 422282.0368440555*eta4 + 584744.0406581408*eta5)*S1);
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_44_rdcp3: version %i is not valid.", RDAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_32_rdaux1(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag){
	double total=0;
	switch (RDAmpFlag){
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
            total = chidiff2*(-4.188795724777721*eta2 + 53.39200466700963*eta3 - 131.19660856923554*eta4) + chidiff1*delta*(14.284921364132623*eta2 - 321.26423637658746*eta3 + 1242.865584938088*eta4) + S1*(-0.022968727462555794*(83.66854837403105*eta1 - 3330.6261333413177*eta2 + 77424.12614733395*eta3 - 710313.3016672594*eta4 + 2.6934917075009225e6*eta5 - 3.572465179268999e6*eta6) + 0.0014795114305436387*(-1672.7273629876313*eta1 + 90877.38260964208*eta2 - 1.6690169155105734e6*eta3 + 1.3705532554135624e7*eta4 - 5.116110998398143e7*eta5 + 7.06066766311127e7*eta6)*S1) + (4.45156488896258*eta1 - 77.39303992494544*eta2 + 522.5070635563092*eta3 - 1642.3057499049708*eta4 + 2048.333892310575*eta5)*pow(1 - 9.611489164758915*eta1 + 24.249594730050312*eta2,-1);
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_32_rdaux1: version %i is not valid.", RDAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_32_rdaux2(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag){
	double total=0;
	switch (RDAmpFlag){
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
            double eta7 = eta1 * eta6;
            double S1 = S;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = chidiff2*(-18.550171209458394*eta2 + 188.99161055445936*eta3 - 440.26516625611*eta4) + chidiff1*delta*(13.132625215315063*eta2 - 340.5204040505528*eta3 + 1327.1224176812448*eta4) + S1*(-0.16707403272774676*(6.678916447469937*eta1 + 1331.480396625797*eta2 - 41908.45179140144*eta3 + 520786.0225074669*eta4 - 3.1894624909922685e6*eta5 + 9.51553823212259e6*eta6 - 1.1006903622406831e7*eta7) + 0.015205286051218441*(108.10032279461095*eta1 - 16084.215590200103*eta2 + 462957.5593513407*eta3 - 5.635028227588545e6*eta4 + 3.379925277713386e7*eta5 - 9.865815275452062e7*eta6 + 1.1201307979786257e8*eta7)*S1) + (3.902154247490771*eta1 - 55.77521071924907*eta2 + 294.9496843041973*eta3 - 693.6803787318279*eta4 + 636.0141528226893*eta5)*pow(1 - 8.56699762573719*eta1 + 19.119341007236955*eta2,-1);
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_32_rdaux2: version %i is not valid.", RDAmpFlag);}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_32_rdcp1(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag){
	double total=0;
	switch (RDAmpFlag){
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
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = chidiff2*(-261.63903838092017*eta3 + 2482.4929818200458*eta4 - 5662.765952006266*eta5) + chidiff1*delta*(200.3023530582654*eta3 - 3383.07742098347*eta4 + 11417.842708417566*eta5) + chidiff2*(-177.2481070662751*eta3 + 1820.8637746828358*eta4 - 4448.151940319403*eta5)*S1 + chidiff1*delta*(412.749304734278*eta3 - 4156.641392955615*eta4 + 10116.974216563232*eta5)*S1 + S1*(-0.07383539239633188*(40.59996146686051*eta1 - 527.5322650311067*eta2 + 4167.108061823492*eta3 - 13288.883172763119*eta4 - 23800.671572828596*eta5 + 146181.8016013141*eta6) + 0.03576631753501686*(-13.96758180764024*eta1 - 797.1235306450683*eta2 + 18007.56663810595*eta3 - 151803.40642097822*eta4 + 593811.4596071478*eta5 - 878123.747877138*eta6)*S1 + 0.01007493097350273*(-27.77590078264459*eta1 + 4011.1960424049857*eta2 - 152384.01804465035*eta3 + 1.7595145936445233e6*eta4 - 7.889230647117076e6*eta5 + 1.2172078072446395e7*eta6)*S2) + (4.146029818148087*eta1 - 61.060972560568054*eta2 + 336.3725848841942*eta3 - 832.785332776221*eta4 + 802.5027431944313*eta5)*pow(1 - 8.662174796705683*eta1 + 19.288918757536685*eta2,-1);
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_32_rdcp1: version is not valid. Recommended version is 122022.");}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_32_rdcp2(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag){
	double total=0;
	switch (RDAmpFlag){
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
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = chidiff2*(-220.42133216774002*eta3 + 2082.031407555522*eta4 - 4739.292554291661*eta5) + chidiff1*delta*(179.07548162694007*eta3 - 2878.2078963030094*eta4 + 9497.998559135678*eta5) + chidiff2*(-128.07917402087625*eta3 + 1392.4598433465628*eta4 - 3546.2644951338134*eta5)*S1 + chidiff1*delta*(384.31792882093424*eta3 - 3816.5687272960417*eta4 + 9235.479593415908*eta5)*S1 + S1*(-0.06144774696295017*(35.72693522898656*eta1 - 168.08433700852038*eta2 - 3010.678442066521*eta3 + 45110.034521934074*eta4 - 231569.4154711447*eta5 + 414234.84895584086*eta6) + 0.03663881822701642*(-22.057692852225696*eta1 + 223.9912685075838*eta2 - 1028.5261783449762*eta3 - 12761.957255385*eta4 + 141784.13567610556*eta5 - 328718.5349981628*eta6)*S1 + 0.004849853669413881*(-90.35491669965123*eta1 + 19286.158446325957*eta2 - 528138.5557827373*eta3 + 5.175061086459432e6*eta4 - 2.1142182400264673e7*eta5 + 3.0737963347449116e7*eta6)*S2) + (3.133378729082171*eta1 - 45.83572706555282*eta2 + 250.23275606463622*eta3 - 612.0498767005383*eta4 + 580.3574091493459*eta5)*pow(1 - 8.698032720488515*eta1 + 19.38621948411302*eta2,-1);
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_32_rdcp2: version is not valid. Recommended version is 122022.");}
  }
  return total;
}

static double IMRPhenomXHM_RD_Amp_32_rdcp3(IMRPhenomXWaveformStruct *pWF, int RDAmpFlag){
	double total=0;
	switch (RDAmpFlag){
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
            total = chidiff2*(-79.14146757219045*eta3 + 748.8207876524461*eta4 - 1712.3401586150026*eta5) + chidiff1*delta*(65.1786095079065*eta3 - 996.4553252426255*eta4 + 3206.5675278160684*eta5) + chidiff2*(-36.474455088940225*eta3 + 421.8842792746865*eta4 - 1117.0227933265749*eta5)*S1 + chidiff1*delta*(169.07368933925878*eta3 - 1675.2562326502878*eta4 + 4040.0077967763787*eta5)*S1 + S1*(-0.01992370601225598*(36.307098892574196*eta1 - 846.997262853445*eta2 + 16033.60939445582*eta3 - 138800.53021166887*eta4 + 507922.88946543116*eta5 - 647376.1499824544*eta6) + 0.014207919520826501*(-33.80287899746716*eta1 + 1662.2913368534057*eta2 - 31688.885017467597*eta3 + 242813.43893659746*eta4 - 793178.4767168422*eta5 + 929016.897093022*eta6)*S1) + (0.9641853854287679*eta1 - 13.801372413989519*eta2 + 72.80610853168994*eta3 - 168.65551450831953*eta4 + 147.2372582604103*eta5)*pow(1 - 8.65963828355163*eta1 + 19.112920222001367*eta2,-1);
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Amp_32_rdcp3: version is not valid. Recommended version is 122022.");}
  }
  return total;
}

/* End of Amp Parameter Space Fits */


/************** Ringdown coefficients from collocation points *************/

static void IMRPhenomXHM_RD_Amp_Coefficients(IMRPhenomXWaveformStruct *pWF22, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMAmpCoefficients *pAmp){
    switch (pWFHM->IMRPhenomXHMRingdownAmpVersion){
        case 0:{
            // We have three "fitted" coefficients across parameter space: alambda, lambda and sigma. Sigma will be constat for all the modes except the 21.
            pAmp->RDCoefficient[0] = fabs(pAmp->RingdownAmpFits[pWFHM->modeInt*3](pWF22,pWFHM->IMRPhenomXHMRingdownAmpFitsVersion));
            pAmp->RDCoefficient[1] = pAmp->RingdownAmpFits[pWFHM->modeInt*3+1](pWF22,pWFHM->IMRPhenomXHMRingdownAmpFitsVersion);
            pAmp->RDCoefficient[2] = pAmp->RingdownAmpFits[pWFHM->modeInt*3+2](pWF22,pWFHM->IMRPhenomXHMRingdownAmpFitsVersion);
            pAmp->RDCoefficient[3] = 1./12.;
            break;
        }
        case 1:
        case 2:{
            if (pWFHM->IMRPhenomXHMRingdownAmpVersion == 1){
                pAmp->RDCoefficient[0] = fabs(pAmp->RingdownAmpFits[pWFHM->modeInt*3](pWF22,pWFHM->IMRPhenomXHMRingdownAmpFitsVersion));
                pAmp->RDCoefficient[1] = pAmp->RingdownAmpFits[pWFHM->modeInt*3+1](pWF22,pWFHM->IMRPhenomXHMRingdownAmpFitsVersion);
                pAmp->RDCoefficient[2] = pAmp->RingdownAmpFits[pWFHM->modeInt*3+2](pWF22,pWFHM->IMRPhenomXHMRingdownAmpFitsVersion);
            }
            else if (pWFHM->IMRPhenomXHMRingdownAmpVersion == 2){
                double rdcp1 = fabs(pAmp->RingdownAmpFits[12 + pWFHM->modeInt * 3](pWF22, pWFHM->IMRPhenomXHMRingdownAmpFitsVersion));
                double rdcp2 = fabs(pAmp->RingdownAmpFits[13 + pWFHM->modeInt * 3](pWF22, pWFHM->IMRPhenomXHMRingdownAmpFitsVersion));
                double rdcp3 = fabs(pAmp->RingdownAmpFits[14 + pWFHM->modeInt * 3](pWF22, pWFHM->IMRPhenomXHMRingdownAmpFitsVersion));

                pAmp->CollocationPointsFreqsAmplitudeRD[0] = pWFHM->fRING - pWFHM->fDAMP;
                pAmp->CollocationPointsFreqsAmplitudeRD[1] = pWFHM->fRING;
                pAmp->CollocationPointsFreqsAmplitudeRD[2] = pWFHM->fRING + pWFHM->fDAMP;
                /* Apply vetos to RDCP. Assuming they are strain */
                // float rdveto = 0.01; // Applied over the strain / RescaleFactor_lm
                // IMRPhenomX_UsefulPowers powers_of_RDCP1, powers_of_RDCP2, powers_of_RDCP3; // PN power v^{1/3} = (2pif/m)
                // IMRPhenomX_Initialize_Powers(&powers_of_RDCP1, pAmp->CollocationPointsFreqsAmplitudeRD[0]);
                // IMRPhenomX_Initialize_Powers(&powers_of_RDCP2, pAmp->CollocationPointsFreqsAmplitudeRD[1]);
                // IMRPhenomX_Initialize_Powers(&powers_of_RDCP3, pAmp->CollocationPointsFreqsAmplitudeRD[2]);
                // double rescale_factor_lm;
                // rescale_factor_lm = RescaleFactor(&powers_of_RDCP1, pAmp, 2);
                // if ( rdcp1 / rescale_factor_lm < rdveto ){
                //     rdcp1 = 0.9 * rdcp2;
                // }
                // rescale_factor_lm = RescaleFactor(&powers_of_RDCP2, pAmp, 2);
                // if ( rdcp2 / rescale_factor_lm < rdveto ){
                //     rdcp2 = 0.9 * rdcp1;
                // }
                // rescale_factor_lm = RescaleFactor(&powers_of_RDCP3, pAmp, 2);
                // if ( rdcp3 / rescale_factor_lm < rdveto ){
                //     rdcp3 = 0.9 * rdcp2;
                // }
                if ( rdcp3 >= rdcp2 * rdcp2 / rdcp1 ){
                    rdcp3 = 0.5 * rdcp2 * rdcp2 / rdcp1;
                }
                if ( rdcp3 > rdcp2 ){
                    rdcp3 = 0.5 * rdcp2;
                }
                if (rdcp1 < rdcp2 && rdcp3 > rdcp1){
                    rdcp3 = rdcp1;
                }
                /* End of vetos */
                pAmp->CollocationPointsValuesAmplitudeRD[0] = rdcp1;
                pAmp->CollocationPointsValuesAmplitudeRD[1] = rdcp2;
                pAmp->CollocationPointsValuesAmplitudeRD[2] = rdcp3;
                REAL8 deno = (sqrt(rdcp1 / rdcp3) - (rdcp1 / rdcp2));
                if (deno <= 0){
                    deno = 1e-16;
                }
                pAmp->RDCoefficient[0] = rdcp1 * pWFHM->fDAMP / deno;
                pAmp->RDCoefficient[2] = sqrt(pAmp->RDCoefficient[0] / (rdcp2 * pWFHM->fDAMP));
                pAmp->RDCoefficient[1] = 0.5 * pAmp->RDCoefficient[2] * log(rdcp1 / rdcp3);
            }
                       
            if (pWFHM->RingdownAmpVeto == 1){
              pAmp->RDCoefficient[1] = 1;
              pAmp->RDCoefficient[2] = 1.35;
              pAmp->RDCoefficient[0] = pAmp->CollocationPointsValuesAmplitudeRD[2] * pWFHM->fDAMP * pAmp->RDCoefficient[2] * pAmp->RDCoefficient[2];
            }
            
            if(pWFHM->fAmpRDfalloff > 0){
                IMRPhenomX_UsefulPowers powers_of_RDfalloff;
                IMRPhenomX_Initialize_Powers(&powers_of_RDfalloff, pWFHM->fAmpRDfalloff);
                REAL8 tmp = pWFHM->fAmpRDfalloff;
                pWFHM->fAmpRDfalloff = 0;
                pAmp->RDCoefficient[3] = IMRPhenomXHM_RD_Amp_Ansatz(&powers_of_RDfalloff, pWFHM, pAmp);
                pAmp->RDCoefficient[4] = -1. * IMRPhenomXHM_RD_Amp_DAnsatz(&powers_of_RDfalloff, pWFHM, pAmp) / pAmp->RDCoefficient[3];
                pWFHM->fAmpRDfalloff = tmp;
            }
            if(pAmp->nCoefficientsRDAux > 0)
                IMRPhenomXHM_RDAux_Amp_Coefficients(pWF22, pWFHM, pAmp);
            break;
        }
        default:{
            XLAL_ERROR_VOID(XLAL_EINVAL, "Error in IMRPhenomXHM_RD_Amp_Coefficients: IMRPhenomXHMRingdownAmpVersion is not valid.\n");
        }
    }
}

static void IMRPhenomXHM_RDAux_Amp_Coefficients(IMRPhenomXWaveformStruct *pWF22, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMAmpCoefficients *pAmp){

    for(UINT2 i = 0; i < pAmp->nCollocPtsRDAux; i++)
        pAmp->CollocationPointsValuesAmplitudeRDAux[i] = fabs(pAmp->RingdownAmpFits[24 + i](pWF22, pWFHM->IMRPhenomXHMRingdownAmpFitsVersion));
    IMRPhenomX_UsefulPowers powers_of_fRDAux;
    IMRPhenomX_Initialize_Powers(&powers_of_fRDAux, pAmp->fRDAux);
    pAmp->CollocationPointsValuesAmplitudeRDAux[pAmp->nCollocPtsRDAux] = IMRPhenomXHM_RD_Amp_Ansatz(&powers_of_fRDAux, pWFHM, pAmp); //pAmp->CollocationPointsValuesAmplitudeRD[0];
    pAmp->CollocationPointsValuesAmplitudeRDAux[pAmp->nCollocPtsRDAux + 1] = IMRPhenomXHM_RD_Amp_DAnsatz(&powers_of_fRDAux, pWFHM, pAmp);
    if (pWFHM->RingdownAmpVeto == 2 && pAmp->CollocationPointsValuesAmplitudeRDAux[pAmp->nCollocPtsRDAux - 1] < pAmp->CollocationPointsValuesAmplitudeRDAux[pAmp->nCollocPtsRDAux]){
        pAmp->CollocationPointsValuesAmplitudeRDAux[pAmp->nCollocPtsRDAux - 1] = pAmp->CollocationPointsValuesAmplitudeRDAux[pAmp->nCollocPtsRDAux];
    }
    
    pAmp->CollocationPointsFreqsAmplitudeRDAux[0] = pAmp->fAmpMatchIM;
    pAmp->CollocationPointsFreqsAmplitudeRDAux[1] = 0.5 * (pAmp->fAmpMatchIM + pAmp->fRDAux); // First Chebyshev node
    pAmp->CollocationPointsFreqsAmplitudeRDAux[2] = pAmp->fRDAux;
    pAmp->CollocationPointsFreqsAmplitudeRDAux[3] = pAmp->fRDAux;


    /* GSL objects for solving system of equations via LU decomposition */
    gsl_vector *b, *x;
    gsl_matrix *A;
    gsl_permutation *p;
    int signum; // No need to set, used internally by gsl_linalg_LU_decomp

    p = gsl_permutation_alloc(pAmp->nCoefficientsRDAux);
    b = gsl_vector_alloc(pAmp->nCoefficientsRDAux);
    x = gsl_vector_alloc(pAmp->nCoefficientsRDAux);
    A = gsl_matrix_alloc(pAmp->nCoefficientsRDAux, pAmp->nCoefficientsRDAux);

    /* Define linear system of equations */

    //FIXME: be careful with the indexing, where are the RDaux CollocPoints in CollocationPointsValuesAmplitudeRD?
    // Should be at the end, although this region goes before than the the "normal RD region".
    for(INT4 i = 0; i < pAmp->nCoefficientsRDAux; i++){
      // b is the vector with the values of collocation points
      gsl_vector_set(b, i, pAmp->CollocationPointsValuesAmplitudeRDAux[i]);
      //FIXME: distinguish InterAmp ansatzaes versions
      // Set system matrix: Polynomial at the collocation points frequencies + derivative at the right boundary
      /* A = (1, f1, f1^2, f1^3, f1^4)
             (1, f2, f2^2, f2^3, f2^4)
             (1, f3, f3^2, f3^3, f3^4)
             (0,  1,   f3, f3^2, f3^3)
             Until number of collocation points
      */
      REAL8 fcollpoint = pAmp->CollocationPointsFreqsAmplitudeRDAux[i];
      REAL8 fpower = 1.; // 1, f, f^2, f^3, f^4, ...
      if (i < pAmp->nCoefficientsRDAux - 1){
          for(INT4 j = 0; j < pAmp->nCoefficientsRDAux; j++){
              gsl_matrix_set(A, i, j, fpower);
              fpower *= fcollpoint;
          }
      }
      else{ // Last row of the matrix for the derivative
          fpower = 1.;
          gsl_matrix_set(A, i, 0, 0.);
          for(INT4 j = 1; j < pAmp->nCoefficientsRDAux; j++){
              gsl_matrix_set(A, i, j, j * fpower);
              fpower *= fcollpoint;
          }
      }
    }

    /* We now solve the system A x = b via an LU decomposition. x is the solution vector */
    gsl_linalg_LU_decomp(A, p, &signum);
    gsl_linalg_LU_solve(A, p, b, x);

    for (INT4 i = 0; i < pAmp->nCoefficientsRDAux; i++){
        pAmp->RDAuxCoefficient[i] = gsl_vector_get(x, i);
    }

    gsl_vector_free(b);
    gsl_vector_free(x);
    gsl_matrix_free(A);
    gsl_permutation_free(p);
}

/************** Amplitude Ringdown Ansatz *************/

// For the modes with mixing this is the ansatz of the spheroidal part.
static double IMRPhenomXHM_RD_Amp_Ansatz(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXHMAmpCoefficients *pAmp){

    double ff = powers_of_Mf->itself;
    int RDAmpFlag = pWFHM->IMRPhenomXHMRingdownAmpVersion;
    double frd = pWFHM->fRING;
    double fda = pWFHM->fDAMP;
    double dfr = ff - frd;
    double ampRD = 0.;

    switch ( RDAmpFlag )
    {
        case 0: /* Canonical, 3 fitted coefficients + fring, fdamp, lc that are fixed. sigma is also fixed except for the 21 mode. */
        {   // Only for the 122018 release.
            double dfd = fda * pAmp->RDCoefficient[2];
            double lc  = pAmp->RDCoefficient[3];
            ampRD = (fda *fabs(pAmp->RDCoefficient[0]) * pAmp->RDCoefficient[2])*exp(- dfr * pAmp->RDCoefficient[1] / dfd )/ (dfr*dfr + dfd*dfd)*pow(ff,-lc);
            // The line below returns the strain amplitude
            // if (pAmp->RDRescaleFactor == 0){
            //      ampRD *= (pWFHM->ampNorm * powers_of_Mf->m_seven_sixths);
            //      //printf("%.10f %.16e\n", ff, ampRD);
            //  }
            break;
        }
        case 1:
        case 2:
        {
            if(pAmp->nCoefficientsRDAux > 0 && !IMRPhenomX_StepFuncBool(ff, pAmp->fRDAux)){
                 /* Polynomial */
                double fpower = 1.;
                for (UINT2 i = 0; i < pAmp->nCoefficientsRDAux; i++){
                    ampRD += fpower * pAmp->RDAuxCoefficient[i];
                    fpower *= ff;
                }
            }
            else if (pWFHM->fAmpRDfalloff > 0 && IMRPhenomX_StepFuncBool(ff, pWFHM->fAmpRDfalloff)){
                ampRD = pAmp->RDCoefficient[3] * exp(- pAmp->RDCoefficient[4] * (ff - pWFHM->fAmpRDfalloff));
            }
            else{ /* Lorentzian with exponential falloff */
                double dfd = fda * pAmp->RDCoefficient[2];
                ampRD = pAmp->RDCoefficient[0] * fda / ( exp(pAmp->RDCoefficient[1] / dfd * dfr) * (dfr * dfr + dfd * dfd)); // * pWF->ampNorm * factor;
            }
            if (pAmp->RDRescaleFactor !=0){
                ampRD /= RescaleFactor(powers_of_Mf, pAmp, pAmp->RDRescaleFactor);
            }
            break;
        }
        default:
        {
            XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomXHM_RD_Amp_Ansatz: IMRPhenomXHMRingdownAmpVersion = %i is not valid. \n", RDAmpFlag);
        }
    }
    
    return ampRD;
}


/*** Derivative of the RD Ansatz for modes without mixing ***/

static double IMRPhenomXHM_RD_Amp_DAnsatz(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMWaveformStruct *pWF, IMRPhenomXHMAmpCoefficients *pAmp){

    double ff = powers_of_Mf->itself;
    int RDAmpFlag = pWF->IMRPhenomXHMRingdownAmpVersion;
    double frd = pWF->fRING;
    double fda = pWF->fDAMP;
    double DampRD;
    double numerator,denom;

    switch ( RDAmpFlag )
    {
        case 0:  /* Canonical, 3 fitted coefficients + fring, fdamp, lc that are fixed. sigma is also fixed except for the 21 mode. */
        {
            double dfd = fda * pAmp->RDCoefficient[2];
            double lc  = pAmp->RDCoefficient[3];
            double lambda = pAmp->RDCoefficient[1], expon;
            numerator = fabs(pAmp->RDCoefficient[0])*pow(ff,-1.-lc)*( dfd*lc*(frd*frd + dfd*dfd)
                                                           + ff*( frd*frd*lambda - 2.*dfd*frd*(1.+lc) + dfd*dfd*lambda)
                                                           + ff*ff*( -2.*lambda*frd + dfd*(2.+lc) )
                                                           + ff*ff*ff*( lambda )
                                                           );
            denom = frd*frd + dfd*dfd - ff*2.*frd + ff*ff;
            expon = exp(-(ff-frd)*lambda/dfd);
            DampRD = -numerator*expon/(denom*denom);
            break;
        }
        case 1:
        case 2:
        {
            double dfr = ff - frd;
            numerator = pAmp->RDCoefficient[0] * (dfr * dfr * pAmp->RDCoefficient[1] + 2 * fda * dfr * pAmp->RDCoefficient[2] + fda * fda * pAmp->RDCoefficient[1] * pAmp->RDCoefficient[2] * pAmp->RDCoefficient[2]);
            denom = (dfr * dfr + fda * fda * pAmp->RDCoefficient[2] * pAmp->RDCoefficient[2]);
            denom = pAmp->RDCoefficient[2] * denom * denom * exp(dfr * pAmp->RDCoefficient[1] / (fda * pAmp->RDCoefficient[2]));
            DampRD = - numerator / denom;
            break;
        }
        default:
        {
            XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomXHM_RD_Amp_DAnsatz: IMRPhenomXHMRingdownAmpVersion = %i is not valid. \n", RDAmpFlag);
        }
    }

    return DampRD;
}

/*** Derivative of the RD Ansatz for modes with mixing ***/
// It can not be obtained analytically, so we use finite difference of 4th order
static double IMRPhenomXHM_RD_Amp_NDAnsatz(IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXHMAmpCoefficients *pAmp,  IMRPhenomXHMPhaseCoefficients *pPhase, IMRPhenomXHMWaveformStruct *pWFHM, IMRPhenomXAmpCoefficients *pAmp22,  IMRPhenomXPhaseCoefficients *pPhase22, IMRPhenomXWaveformStruct *pWF22){

    double ff = powers_of_Mf->itself;
    double df = 10e-10;
    double Nder;
    double fun2R = ff + 2*df;
    double funR  = ff + df;
    double funL  = ff - df;
    double fun2L = ff - 2*df;
    IMRPhenomX_UsefulPowers powers_of_fun;

    IMRPhenomX_Initialize_Powers(&powers_of_fun,fun2R);
    fun2R = cabs(SpheroidalToSpherical(&powers_of_fun, pAmp22, pPhase22, pAmp, pPhase, pWFHM, pWF22));

    IMRPhenomX_Initialize_Powers(&powers_of_fun,funR);
    funR  = cabs(SpheroidalToSpherical(&powers_of_fun, pAmp22, pPhase22, pAmp, pPhase, pWFHM, pWF22));

    IMRPhenomX_Initialize_Powers(&powers_of_fun,funL);
    funL  = cabs(SpheroidalToSpherical(&powers_of_fun, pAmp22, pPhase22, pAmp, pPhase, pWFHM, pWF22));

    IMRPhenomX_Initialize_Powers(&powers_of_fun,fun2L);
    fun2L = cabs(SpheroidalToSpherical(&powers_of_fun, pAmp22, pPhase22, pAmp, pPhase, pWFHM, pWF22));

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

/* Start of Phase Parameter Space Fits */

static double IMRPhenomXHM_RD_Phase_22_alpha2(IMRPhenomXWaveformStruct *pWF, int RDPhaseFlag) {
    double total=0;
    switch (RDPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,delta=sqrt(1-4*eta);
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            double noSpin = 0.2088669311744758 - 0.37138987533788487*eta + 6.510807976353186*eta2 - 31.330215053905395*eta3 + 55.45508989446867*eta4;
            double eqSpin = ((0.2393965714370633 + 1.6966740823756759*eta - 16.874355161681766*eta2 + 38.61300158832203*eta3)*S)/(1. - 0.633218538432246*S);
            double uneqSpin = pWF->dchi*(0.9088578269496244*pow(eta,2.5) + 15.619592332008951*pWF->dchi*pow(eta,3.5))*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Phase_22_alpha2: version % is not valid.", RDPhaseFlag);}
    }
    return total;
}

static double IMRPhenomXHM_RD_Phase_22_alphaL(IMRPhenomXWaveformStruct *pWF, int RDPhaseFlag) {
    double total=0;
    switch (RDPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double delta=sqrt(1.- 4.*eta),eta2,eta3,eta4,S2;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            S2 = pow(S,2);
            double noSpin = eta*(-1.1926122248825484 + 2.5400257699690143*eta - 16.504334734464244*eta2 + 27.623649807617376*eta3);
            double eqSpin = eta3*S*(35.803988443700824 + 9.700178927988006*S - 77.2346297158916*S2) + eta*S*(0.1034526554654983 - 0.21477847929548569*S - 0.06417449517826644*S2) + eta2*S*(-4.7282481007397825 + 0.8743576195364632*S + 8.170616575493503*S2) + eta4*S*(-72.50310678862684 - 39.83460092417137*S + 180.8345521274853*S2);
            double uneqSpin = (-0.7428134042821221*pWF->chi1L*pow(eta,3.5) + 0.7428134042821221*pWF->chi2L*pow(eta,3.5) + 17.588573345324154*pow(pWF->chi1L,2)*pow(eta,4.5) - 35.17714669064831*pWF->chi1L*pWF->chi2L*pow(eta,4.5) + 17.588573345324154*pow(pWF->chi2L,2)*pow(eta,4.5))*delta;
            total = noSpin + eqSpin + uneqSpin;
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Phase_22_alphaL: version % is not valid.", RDPhaseFlag);}
    }
    return total;
}

/**************** 32 specific fits ***************/
static double IMRPhenomXHM_RD_Phase_32_SpheroidalTimeShift(IMRPhenomXWaveformStruct *pWF, int RDPhaseFlag) {
    double total;
    switch (RDPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,eta5,S2,S3,S4;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            eta4 = pow(eta,4);
            eta5 = pow(eta,5);
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4 = pow(S,4);
            double noSpin = 11.851438981981772 + 167.95086712701223*eta - 4565.033758777737*eta2 + 61559.132976189896*eta3 - 364129.24735853914*eta4 + 739270.8814129328*eta5;
            double eqSpin = (9.506768471271634 + 434.31707030999445*eta - 8046.364492927503*eta2 + 26929.677144312944*eta3)*S + (-5.949655484033632 - 307.67253970367034*eta + 1334.1062451631644*eta2 + 3575.347142399199*eta3)*S2 + (3.4881615575084797 - 2244.4613237912527*eta + 24145.932943269272*eta2 - 60929.87465551446*eta3)*S3 + (15.585154698977842 - 2292.778112523392*eta + 24793.809334683185*eta2 - 65993.84497923202*eta3)*S4;
            double uneqSpin = 465.7904934097202*pWF->dchi*sqrt(1.-4.*eta)*eta2;
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
            double S1 = S;
            double chidiff1 = chidiff;
            total = chidiff1*delta*(-3437.4083682807154*eta2 + 29349.322789666763*eta3 - 46822.26853496922*eta4) + S1*(18.280316689743625*(46.16185108285708*eta1 - 515.3588338752849*eta2 + 1475.462268010926*eta3) + 2.1970246269696427*(370.70040479593024*eta1 - 4329.153306191607*eta2 + 12065.631680744644*eta3)*S1) + (12.080898026205173 - 79.10914761468462*eta1 - 89.82426456495799*eta2 + 915.6864093792078*eta3)*pow(1 - 9.443713298364061*eta1 + 22.353898754970686*eta2,-1);
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Phase_32_SpheroidalTimeShift: version %i is not valid.", RDPhaseFlag);}
    }
    return total;
}

static double IMRPhenomXHM_RD_Phase_32_SpheroidalPhaseShift(IMRPhenomXWaveformStruct *pWF, int RDPhaseFlag) {
    double total;
    switch (RDPhaseFlag){
        case 122019:
        {
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
            double noSpin = -1.3328895897490733 - 22.209549522908667*eta + 1056.2426481245027*eta2 - 21256.376324666326*eta3 + 246313.12887984765*eta4 - 1.6312968467540336e6*eta5 + 5.614617173188322e6*eta6 - 7.612233821752137e6*eta7;
            double eqSpin = (S*(-1.622727240110213 + 0.9960210841611344*S - 1.1239505323267036*S2 - 1.9586085340429995*S3 + eta2*(196.7055281997748 + 135.25216875394943*S + 1086.7504825459278*S2 + 546.6246807461155*S3 - 312.1010566468068*S4) + 0.7638287749489343*S4 + eta*(-47.475568056234245 - 35.074072557604445*S - 97.16014978329918*S2 - 34.498125910065156*S3 + 24.02858084544326*S4) + eta3*(62.632493533037625 - 22.59781899512552*S - 2683.947280170815*S2 - 1493.177074873678*S3 + 805.0266029288334*S4)))/(-2.950271397057221 + 1.*S);
            double uneqSpin = (sqrt(1.-4.*eta)*(pWF->chi2L*pow(eta,2.5)*(88.56162028006072 - 30.01812659282717*S) + pWF->chi2L*eta2*(43.126266433486435 - 14.617728550838805*S) + pWF->chi1L*eta2*(-43.126266433486435 + 14.617728550838805*S) + pWF->chi1L*pow(eta,2.5)*(-88.56162028006072 + 30.01812659282717*S)))/(-2.950271397057221 + 1.*S);
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
            double S1 = S;
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            double chidiff2 = chidiff1 * chidiff1;
            total = chidiff1*delta*(-4055.661463620154*eta3 + 48053.79215999518*eta4 - 126644.96534635697*eta5) + chidiff2*(1489.2939046368886*eta3 - 11333.726790227513*eta4 + 21470.681301598026*eta5) + S1*(-0.09759112086805133*(-3.6379140102351064 + 510.0158912180661*eta1 - 12528.715040030444*eta2 + 82416.24893999398*eta3 - 161740.0041427807*eta4) - 0.20117612026208484*(-5.676646590653427 + 183.78422258983136*eta1 - 3807.617101722895*eta2 + 24219.360141326677*eta3 - 45892.369216999985*eta4)*S1 + 0.06537048555519695*(21.388069487574892 + 408.2245599781871*eta1 - 9652.666650065075*eta2 + 57782.086859487965*eta3 - 110112.73697613904*eta4)*S2) + (-1.3903824533899325 + 13.761709564309667*eta1 - 10.633427224975128*eta2 - 263.01839936998965*eta3 + 736.1912896690361*eta4)*pow(1 - 9.458921853648649*eta1 + 22.932673653997934*eta2,-1);
            break;
        }
            default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Phase_32_SpheroidalPhaseShift: version %i is not valid.", RDPhaseFlag);}
    }
    return total;
}

// collocation points
static double IMRPhenomXHM_RD_Phase_32_p1(IMRPhenomXWaveformStruct *pWF, int RDPhaseFlag) {
    double total;
    switch (RDPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,eta4,eta5,S2,S3,S4,S5;
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
            double uneqSpin = pWF->dchi*sqrt(1.-4.*eta)*eta2*(641.2473192044652 - 1600.240100295189*pWF->chi1L*eta + 1600.240100295189*pWF->chi2L*eta + 13275.623692212472*eta*S);
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
            double S1 = S;
            double S2 = S1 * S1;
            double chidiff1 = chidiff;
            total = 3221.932636435056 - 134.59521937837278*eta1 + chidiff1*delta*(319.06177738567885*eta2 - 20578.40007594641*eta3 + 101859.1659970414*eta4) + S1*(24.059115241154707*(436.41245673494626*eta2 - 2938.503437122471*eta3 + 5027.414440730744*eta4) + 21.848741292340005*(1251.706577839354*eta2 - 14171.490147583942*eta3 + 36914.6553449061*eta4)*S1 + 15.901300902033508*(-149.2789474539545*eta2 + 3483.608736789833*eta3 - 11289.97178789606*eta4)*S2)*pow(-1.803337035190313 + S1,-1) - 0.012868288384497346*pow(0.000187943994873322 + pow(-0.1923690355128322 + eta1,2),-1);
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Phase_32_p1: version % is not valid.", RDPhaseFlag);}
    }
    return total;
}
// collocation points
static double IMRPhenomXHM_RD_Phase_32_p2(IMRPhenomXWaveformStruct *pWF, int RDPhaseFlag) {
    double total;
    switch (RDPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,S2,S3,S4;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4 = pow(S,4);
            double noSpin = 3131.0260952676376 + 206.09687819102305*eta - 2636.4344627081873*eta2 + 7475.062269742079*eta3;
            double eqSpin = (49.90874152040307 - 691.9815135740145*eta - 434.60154548208334*eta2 + 10514.68111669422*eta3)*S + (97.3078084654917 - 3458.2579971189534*eta + 26748.805404989867*eta2 - 56142.13736008524*eta3)*S2 + (-132.49105074500454 + 429.0787542102207*eta + 7269.262546204149*eta2 - 27654.067482558712*eta3)*S3 + (-227.8023564332453 + 5119.138772157134*eta - 34444.2579678986*eta2 + 69666.01833764123*eta3)*S4;
            double uneqSpin = 477.51566939885424*pWF->dchi*sqrt(1.-4.*eta)*eta2;
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
            double S1 = S;
            double chidiff1 = chidiff;
            total = 3169.5772463611165 + 56.56534589293562*eta1 - 863.5731390762933*eta2 + 2693.8619211321557*eta3 + chidiff1*delta*(-2818.2944800258847*eta2 + 15684.658457287562*eta3 + 14379.128341035908*eta4) + S1*(-0.16388886708177886*(334.30009385854424*eta2 - 154749.10305716714*eta3 + 613903.6107269318*eta4) + 11.950465013745157*(1079.481585746054*eta2 - 11981.85336876442*eta3 + 30911.708103120814*eta4)*S1)*pow(-1.169876031327984 + S1,-1) - 0.009425837438775205*pow(0.00016960223009674388 + pow(-0.20083535429185695 + eta1,2),-1);
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Phase_32_p2: version % is not valid.", RDPhaseFlag);}
    }
    return total;
}
// collocation points
static double IMRPhenomXHM_RD_Phase_32_p3(IMRPhenomXWaveformStruct *pWF, int RDPhaseFlag) {
    double total;
    switch (RDPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,eta3,S2,S3,S4;
            eta2 = pow(eta,2);
            eta3 = pow(eta,3);
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4 = pow(S,4);
            double noSpin = 3082.803556599222 + 76.94679795837645*eta - 586.2469821978381*eta2 + 977.6115755788503*eta3;
            double eqSpin = (45.08944710349874 - 807.7353772747749*eta + 1775.4343704616288*eta2 + 2472.6476419567534*eta3)*S + (95.57355060136699 - 2224.9613131172046*eta + 13821.251641893134*eta2 - 25583.314298758105*eta3)*S2 + (-144.96370424517866 + 2268.4693587493093*eta - 10971.864789147161*eta2 + 16259.911572457446*eta3)*S3 + (-227.8023564332453 + 5119.138772157134*eta - 34444.2579678986*eta2 + 69666.01833764123*eta3)*S4;
            double uneqSpin = 378.2359918274837*pWF->dchi*sqrt(1.-4.*eta)*eta2;
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
            total = 3119.6603946770488 + 168.61554447853712*eta1 - 1777.6654596491376*eta2 + 5037.407962552042*eta3 + chidiff1*delta*(45693.135566736484*eta2 - 789332.4959926775*eta3 + 4.460496312695218e6*eta4 - 8.176309211912101e6*eta5) + chidiff1*delta*(7840.121424232572*eta2 - 47166.09840761356*eta3 + 66597.52917033392*eta4)*S1 + S1*(-6.019579546899472*(14538.163921822728*eta2 - 318911.4362763759*eta3 + 2.6041867832020866e6*eta4 - 9.288489508236282e6*eta5 + 1.2170972980338342e7*eta6) - 7.739304888898913*(11114.219992659659*eta2 - 231541.9569739445*eta3 + 1.8069370995120746e6*eta4 - 6.203273456127891e6*eta5 + 7.874294046591697e6*eta6)*S1) - 0.0003538235766427527*pow(8.2810517855422e-6 + pow(-0.2075868299995718 + eta1,2),-1);
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Phase_32_p3: version % is not valid.", RDPhaseFlag);}
    }
    return total;
}
// collocation points
static double IMRPhenomXHM_RD_Phase_32_p4(IMRPhenomXWaveformStruct *pWF, int RDPhaseFlag) {
    double total;
    switch (RDPhaseFlag){
        case 122019:{
            double eta = pWF->eta;
            double S = pWF->STotR;
            double eta2,S2,S3,S4;
            eta2 = pow(eta,2);
            S2 = pow(S,2);
            S3 = pow(S,3);
            S4 = pow(S,4);
            double noSpin = 3077.0657367004565 + 64.99844502520415*eta - 357.38692756785395*eta2;
            double eqSpin = (34.793450080444714 - 986.7751755509875*eta - 9490.641676924794*pow(eta,3) + 5700.682624203565*eta2)*S + (57.38106384558743 - 1644.6690499868596*eta - 19906.416384606226*pow(eta,3) + 11008.881935880598*eta2)*S2 + (-126.02362949830213 + 3169.3397351803583*eta + 62863.79877094988*pow(eta,3) - 26766.730897942085*eta2)*S3 + (-169.30909412804587 + 4900.706039920717*eta + 95314.99988114933*pow(eta,3) - 41414.05689348732*eta2)*S4;
            double uneqSpin = 390.5443469721231*pWF->dchi*sqrt(1.-4.*eta)*eta2;
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
            double eta7 = eta1 * eta6;
            double S1 = S;
            double chidiff1 = chidiff;
            total = 3063.2702533356364 + 74.20321762511647*eta1 - 346.6653326379183*eta2 + chidiff1*delta*(2604.6711121030685*eta2 - 25322.83641432119*eta3 + 64521.907625802785*eta4) + 17.748987975469845*(6659.676835974436*eta2 - 207712.54687648916*eta3 + 2.5247192995989644e6*eta4 - 1.4825576629165621e7*eta5 + 4.215660954626601e7*eta6 - 4.662796037240443e7*eta7)*S1*pow(-1.2617538728082525 + S1,-1) - 2.211595095940034e-6*pow(-0.000010240443266599763 + pow(-0.17594819839300638 + eta1,2),-1);
            break;
        }
        default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Phase_32_p4: version % is not valid.", RDPhaseFlag);}
    }
    return total;
}


static double IMRPhenomXHM_RD_Phase_32_p5(IMRPhenomXWaveformStruct *pWF, int RDPhaseFlag){
	double total=0;
	switch (RDPhaseFlag){
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
            double eta7 = eta1 * eta6;
            double S1 = S;
            double chidiff1 = chidiff;
            total = 3099.009925231132 - 5.823302919769178*eta1 + chidiff1*delta*(1936.2628510750844*eta2 - 32280.62038532877*eta3 + 97943.49145078743*eta4) + 2.8281136508769498*(6856.762845446495*eta2 - 201653.94475043533*eta3 + 2.4217205751584964e6*eta4 - 1.4582166806161262e7*eta5 + 4.377208838391116e7*eta6 - 5.214533274483399e7*eta7)*S1*pow(-1.0624647822393556 + S1,-1);
            break;
        }
    default:{XLAL_ERROR_REAL8(XLAL_EINVAL,"Error in IMRPhenomXHM_RD_Phase_32_p5: version % is not valid.", RDPhaseFlag);}
  }
  return total;
}

/* End of Phase Parameter Space Fits */

/**************  ANSATZ PHASE DERIVATIVE **************/

static double IMRPhenomXHM_RD_Phase_Ansatz(double ff, IMRPhenomX_UsefulPowers *powers_of_f,IMRPhenomXHMWaveformStruct *pWFHM,  IMRPhenomXHMPhaseCoefficients *pPhase){

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
            dphaseRD = ( pPhase->alpha0 +  frd*frd*(pPhase->alpha2)*powers_of_f->m_two + ( (pPhase->alphaL)* fda/(fda*fda +(ff - frd)*(ff - frd)) ) );
            break;

        }
        case 1:
        {
            /*  calibration of spheroidal ringdown waveform for (32) */

            if(pWFHM->IMRPhenomXHMRingdownPhaseVersion == 122019){
                /* ansatz: alpha0 + (alpha2)/(f^2)+ (alpha4)/(f^4)  + alphaL*(fdamplm)/((fdamplm)^2 + (f - fRDlm)^2)*/
                dphaseRD = ( pPhase->alpha0_S +  (pPhase->alpha2_S)*powers_of_f->m_two + (pPhase->alpha4_S)*powers_of_f->m_four +( (pPhase->alphaL_S)* fda/(fda*fda +(ff - frd)*(ff - frd)) ) );
            }
            else{ // FIXME: 1/eta???
              if(pWFHM->fPhaseRDflat > 0 && IMRPhenomX_StepFuncBool(ff, pWFHM->fPhaseRDflat)){
                dphaseRD = pPhase->RDCoefficient[5] + pPhase->RDCoefficient[6] * powers_of_f->m_five;
              }
              else{
                /* ansatz: a0 + a1/f + a2/f^2 + a3/f^4  + a4*fdamplm/(fdamplm^2 + (f - fRDlm)^2) */
                dphaseRD = ( pPhase->RDCoefficient[0] +  pPhase->RDCoefficient[1] * powers_of_f->m_one + pPhase->RDCoefficient[2] * powers_of_f->m_two + pPhase->RDCoefficient[3] * powers_of_f->m_four + pPhase->RDCoefficient[4]*fda / (fda*fda + (ff - frd)*(ff - frd))  );
              }
            }
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
            double invf3 = powers_of_f->m_three;

            if(pWFHM->IMRPhenomXHMRingdownPhaseVersion == 122019){
                /* ansatz: f alpha0 - (alpha4)/(3 f^3) - (alpha2)/f + alphaL ArcTan[(f - fRDlm)/fdamplm]*/
                phaseRD = pPhase->phi0_S+pPhase->alpha0_S*ff -(pPhase->alpha2_S)*invf -1./3.*(pPhase->alpha4_S)*invf3 +(pPhase->alphaL_S)* atan((ff-frd)/fda);
            }
            else{
              if(pWFHM->fPhaseRDflat > 0 && IMRPhenomX_StepFuncBool(ff, pWFHM->fPhaseRDflat)){
                // a + b / f^5
                phaseRD = pPhase->phi0_S + pPhase->RDCoefficient[7] + pPhase->RDCoefficient[5]*ff - 0.25 * pPhase->RDCoefficient[6] * powers_of_f->m_four;
              }
              else{
                /* ansatz: f a0 + a1 Log(f) - a2/f - a3/(3f^3)  + a4*ArcTan[(f - fRDlm)/fdamplm] */
                phaseRD = pPhase->phi0_S + pPhase->RDCoefficient[0]*ff + pPhase->RDCoefficient[1]*powers_of_f->log - pPhase->RDCoefficient[2]*invf - 1/3.*pPhase->RDCoefficient[3]*invf3 + pPhase->RDCoefficient[4]*atan( (ff-frd)/fda );
              }
            }
            break;
        }
        default:
        {XLAL_ERROR(XLAL_EDOM, "Error in IMRPhenomXHM_RD_Phase_AnsatzInt: version is not valid. Use version 0 for modes (2,1),(3,3),(4,4) and 1 for (3,2).\n");}
    }
    return phaseRD;
}

static double IMRPhenomXHM_RD_Phase_DerAnsatz(double ff, IMRPhenomX_UsefulPowers *powers_of_f,IMRPhenomXHMWaveformStruct *pWFHM,  IMRPhenomXHMPhaseCoefficients *pPhase){

    double frd   = pWFHM->fRING;
    double fda   = pWFHM->fDAMP;
    double ddphaseRD;
    

    switch ( pWFHM->MixingOn )
    {
        case 0:
        {
            /*  rescaling of the 22 ansatz -- used for (21),(33),(44) */
            /*ansatz:
             alpha0 + ((fRDlm^2) alpha2)/(f^2)  + alphaL*(fdamplm)/((fdamplm)^2 + (f - fRDlm)^2)*/
            ddphaseRD = -2*frd*frd*(pPhase->alpha2)*powers_of_f->m_three + -2*pPhase->alphaL* fda *(ff-frd)/pow(fda*fda +(ff - frd)*(ff - frd), 2) ;
            break;
        }
        case 1:
        {
            /*  calibration of spheroidal ringdown waveform for (32) */
            if(pWFHM->IMRPhenomXHMRingdownPhaseVersion == 122019){
                /* ansatz: alpha0 + (alpha2)/(f^2)+ (alpha4)/(f^4)  + alphaL*(fdamplm)/((fdamplm)^2 + (f - fRDlm)^2)*/
                ddphaseRD = -2*pPhase->alpha2_S*powers_of_f->m_three - 4*(pPhase->alpha4_S)*powers_of_f->m_five - 2*( pPhase->alphaL_S* fda * (ff-frd)/pow(fda*fda +(ff - frd)*(ff - frd), 2) ) ;
            }
            else{ 
                /* ansatz: a0 + a1/f + a2/f^2 + a3/f^4  + a4*fdamplm/(fdamplm^2 + (f - fRDlm)^2) */
                ddphaseRD = ( - pPhase->RDCoefficient[1] * powers_of_f->m_two - 2 * pPhase->RDCoefficient[2] * powers_of_f->m_three - 4 * pPhase->RDCoefficient[3] * powers_of_f->m_five - 2 * pPhase->RDCoefficient[4]*fda*(ff-frd) / pow(fda*fda + (ff - frd)*(ff - frd), 2)  );
            }
            break;
        }
        default:
        {XLAL_ERROR(XLAL_EDOM, "Error in IMRPhenomXHM_RD_Phase_Ansatz: version is not valid. Use version 0 for modes (2,1),(3,3),(4,4) and 1 for (3,2).\n");}
    }
    return ddphaseRD;
}
