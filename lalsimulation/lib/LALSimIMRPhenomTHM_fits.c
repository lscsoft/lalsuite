/*
 * Copyright (C) 2020 Hector Estelles
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


/**
 * \author Hector Estelles
 */

#include <lal/XLALError.h>
#include <math.h>

/**** Set of phenomenological fits employed by IMRPhenomT and IMRPhenomTHM models. Collocation points and coefficients have been calibrated with 531 BBH non-precessing NR simulations from the
last release of the SXS Catalog, additional BAM NR simulations at q=4, q=8 and q=18, and numerical Teukolsky waveforms placed at q=200 and q=1000. Calibration procedure has followed the 
hierarchical data-driven fitting approach (Xisco Jimenez-Forteza et al https://arxiv.org/abs/1611.00332) using the symmetric mass ratio eta, dimensionless effective spin Shat=(m1^2*chi1+m2^2*chi2)/(m1^2+m2^2)
and spin difference dchi=chi1-chi2. Supplementary material for the fits is available at https://git.ligo.org/waveforms/reviews/phenomt/-/tree/master/SupplementaryMaterial/Fits3DPhenomTHM ****/

/******************************* 22 FREQUENCY FITS ********************/

/** Inspiral 22 frequency collocation points **/

static double IMRPhenomT_Inspiral_TaylorT3_t0(double eta, double S, double dchi, double delta) // theta = 0.45
{
 double fit;

    fit = pow(eta,-1)*((-20.74399646637014 - 106.27711276502542*eta)*pow(1 + 0.6516016033332481*eta,-1) + 0.0012450290074562259*dchi*delta*(1 - 4.701633367918768e6*eta)*pow(eta,2) - 
     111.5049997379579*dchi*delta*(1 + 19.95458485773613*eta)*S*pow(eta,2) + 1204.6829118499857*(1 - 4.025474056585855*eta)*pow(dchi,2)*pow(eta,3) + 
     S*(338.7318821277009 - 1553.5891860091408*eta + 19614.263378999745*pow(eta,2) - 156449.78737303324*pow(eta,3) + 577363.3090369126*pow(eta,4) - 802867.433363341*pow(eta,5)) + 
     (-55.75053935847546 - 290.36341163610575*eta + 7873.7667183299345*pow(eta,2) - 43585.59040070178*pow(eta,3) + 87229.84668746481*pow(eta,4) - 32469.263449695136*pow(eta,5))*pow(S,2) + 
     (-102.8269343111326 + 5121.845705262981*eta - 93026.46878769135*pow(eta,2) + 650989.6793529999*pow(eta,3) - 1.8846061037110784e6*pow(eta,4) + 1.861602620702142e6*pow(eta,5))*pow(S,3) + 
     (-7.294950933078567 + 314.24955197427136*eta - 3751.8509582195657*pow(eta,2) + 21205.339564205595*pow(eta,3) - 46448.94771114493*pow(eta,4) + 20310.512558558552*pow(eta,5))*pow(S,4) + 
     (97.22312282683716 - 4556.60375328623*eta + 76308.73046927384*pow(eta,2) - 468784.4188333802*pow(eta,3) + 998692.0246600509*pow(eta,4) - 322905.9042578296*pow(eta,5))*pow(S,5));

    return fit;
   }

static double IMRPhenomT_Inspiral_Freq_CP1_22(double eta, double S, double dchi, double delta) // theta = 0.45
{
	double fit;

   fit = -0.014968864336704284*dchi*delta*(1 - 1.942061808318584*eta)*pow(eta,2) + 0.0017312772309375462*dchi*delta*(1 - 0.07106994121956058*eta)*S*pow(eta,2) + S*(0.0019208448318368731 - 0.0013579968243452476*eta - 0.0033501404728414627*pow(eta,2) + 0.008914420175326192*pow(eta,3)) + 
   6.687615165457298e-6*pow(dchi,2)*pow(eta,3) + (0.02104073275966069 + 717.1534194224539*eta + 85.37320237350282*pow(eta,2) + 12.789214868358362*pow(eta,3) - 16.00243777208413*pow(eta,4))*pow(1 + 32934.586638893634*eta,-1) + 
   (-8.306810248117731e-6 + 0.00009918593182087119*eta - 0.003805916669791129*pow(eta,2) + 0.009854209286892323*pow(eta,3))*pow(S,2) + (-5.578836442449699e-6 - 0.0030378960591856616*eta + 0.03746366675135751*pow(eta,2) - 0.10298471015315146*pow(eta,3))*pow(S,3) + 
   (0.00004425141111368952 - 0.0008702073302258368*eta + 0.006538604805919268*pow(eta,2) - 0.01578597166324495*pow(eta,3))*pow(S,4) + (-0.000019469656288570753 + 0.002969863931498354*eta - 0.03643271052162611*pow(eta,2) + 0.09959495981802587*pow(eta,3))*pow(S,5) + 
   (-0.000042037164406446896 + 0.0007336074135429041*eta - 0.005603356997202016*pow(eta,2) + 0.013439843000090702*pow(eta,3))*pow(S,6);

   	return fit;
   }

static double IMRPhenomT_Inspiral_Freq_CP2_22(double eta, double S, double dchi, double delta) // theta = 0.55
{
	double fit;

   fit = -0.04486391236129559*dchi*delta*(1 - 1.8997912248414794*eta)*pow(eta,2) - 0.003531802135161727*dchi*delta*(1 - 8.001211450141325*eta)*S*pow(eta,2) + S*(0.0061664395419698285 - 0.0040934633081508905*eta - 0.009180337242551828*pow(eta,2) + 0.020338583755834694*pow(eta,3)) + 
   0.00006524644306613066*pow(dchi,2)*pow(eta,3) + pow(1 - 3.2125452791404148*eta,-1)*(0.03711511661217631 - 0.10663782888636487*eta - 0.09963406984414182*pow(eta,2) + 0.6597367702009397*pow(eta,3) - 2.777344875144891*pow(eta,4) + 4.220674345359693*pow(eta,5)) + 
   (0.00044302547647888445 + 0.000424246501303979*eta - 0.01394093576260671*pow(eta,2) + 0.02634851560709597*pow(eta,3))*pow(S,2) + (0.00011582043047950321 - 0.008282652950117982*eta + 0.08965067576998058*pow(eta,2) - 0.23963885130463913*pow(eta,3))*pow(S,3) + 
   (0.0006123158975881322 - 0.007809160444435783*eta + 0.028517174579539676*pow(eta,2) - 0.03717957419042746*pow(eta,3))*pow(S,4) + (-0.0000885530893214531 + 0.005939789043536808*eta - 0.07106551435109858*pow(eta,2) + 0.1891131957235774*pow(eta,3))*pow(S,5) + 
   (-0.0005110853374341054 + 0.0038762476596420855*eta + 0.005094077179675256*pow(eta,2) - 0.047971766995287136*pow(eta,3))*pow(S,6);

   	return fit;
   }

static double IMRPhenomT_Inspiral_Freq_CP3_22(double eta, double S, double dchi, double delta) // theta = 0.65
{
	double fit;

   fit = -0.10196878573773932*dchi*delta*(1 - 1.8918584778973513*eta)*pow(eta,2) - 0.018820536453940443*dchi*delta*(1 - 3.7307154599131183*eta)*S*pow(eta,2) - 0.00013162098437956188*pow(dchi,2)*pow(eta,3) + 
   S*(0.0145572994468378 - 0.0017482433991394227*eta - 0.10299007619034371*pow(eta,2) + 0.4581039376357615*pow(eta,3) - 0.7123678787549022*pow(eta,4)) + 
   (0.05489007025458171 + 5.852073438961151*eta + 2.74597705533403*pow(eta,2) + 4.834336623113389*pow(eta,3) - 26.931994454691022*pow(eta,4) + 57.67035368809743*pow(eta,5))*pow(1 + 105.52132834236778*eta,-1) + 
   (0.003001211395915229 + 0.0017929418998452987*eta - 0.13776590125456148*pow(eta,2) + 0.7471133710854526*pow(eta,3) - 1.3620323111858437*pow(eta,4))*pow(S,2) + 
   (0.001143282743686261 - 0.05793457776296727*eta + 0.7841331051705482*pow(eta,2) - 3.4936244160305323*pow(eta,3) + 4.802357041496856*pow(eta,4))*pow(S,3) + 
   (0.0009168588840889624 - 0.03261437094899735*eta + 0.3472881896838799*pow(eta,2) - 1.3634383958859384*pow(eta,3) + 1.7313939586675267*pow(eta,4))*pow(S,4) + 
   (-0.0002794014744432316 + 0.055911057147527664*eta - 0.8686311380514122*pow(eta,2) + 4.096191294930781*pow(eta,3) - 6.009676060669872*pow(eta,4))*pow(S,5) + 
   (-0.0005046018052528331 + 0.029804593053788925*eta - 0.3792653361049425*pow(eta,2) + 1.6366976231421981*pow(eta,3) - 2.26904099961476*pow(eta,4))*pow(S,6);

   	return fit;
   }

static double IMRPhenomT_Inspiral_Freq_CP4_22(double eta, double S, double dchi, double delta) // theta = 0.75
{
	double fit;

   fit = -0.1831889759662071*dchi*delta*(1 - 1.8484261527766557*eta)*pow(eta,2) - 0.07586202965525136*dchi*delta*(1 - 3.2918162656371983*eta)*S*pow(eta,2) + 0.0019259052728265817*pow(dchi,2)*pow(eta,3) + 
   S*(0.02685637375751212 + 0.013341664908359861*eta - 0.3057217933283597*pow(eta,2) + 1.395763446325911*pow(eta,3) - 2.2559396974665376*pow(eta,4)) + 
   (0.0725639467287476 + 12.39400068457852*eta + 12.907450928972402*pow(eta,2) - 7.422660061864399*pow(eta,3) + 66.32985901506036*pow(eta,4) - 117.85875779454518*pow(eta,5))*pow(1 + 168.63492460136445*eta,-1) + 
   (0.0087781653701194 + 0.006944161553839352*eta - 0.3301149078235105*pow(eta,2) + 1.6835714783903248*pow(eta,3) - 2.950404929598742*pow(eta,4))*pow(S,2) + 
   (0.0037229746496019625 - 0.17155338099487646*eta + 2.5881802140836774*pow(eta,2) - 13.14710199375518*pow(eta,3) + 21.366803256010915*pow(eta,4))*pow(S,3) + 
   (0.00278507305662002 - 0.12475855143364532*eta + 1.8640209516178643*pow(eta,2) - 10.117078727717564*pow(eta,3) + 17.94244821676711*pow(eta,4))*pow(S,4) + 
   (0.0010273954584773936 + 0.1713357629442166*eta - 3.017249223460983*pow(eta,2) + 15.855096360798678*pow(eta,3) - 26.444621592311933*pow(eta,4))*pow(S,5) + 
   (-0.00012207946532225968 + 0.11709700788855186*eta - 2.0950821618097026*pow(eta,2) + 11.925324501640054*pow(eta,3) - 21.683978511818076*pow(eta,4))*pow(S,6);

   	return fit;
   }

static double IMRPhenomT_Inspiral_Freq_CP5_22(double eta, double S, double dchi, double delta) // theta = 0.82
{
	double fit;

   fit = -0.2508206617297265*dchi*delta*(1 - 1.861010982421798*eta)*pow(eta,2) - 0.1392163711259171*dchi*delta*(1 - 3.2669366465555796*eta)*S*pow(eta,2) + 0.0023126403170013045*pow(dchi,2)*pow(eta,3) + 
   S*(0.036750064163293766 + 0.036904343404333906*eta - 0.5238739410356437*pow(eta,2) + 2.3292117112945223*pow(eta,3) - 3.654184701923543*pow(eta,4)) + (0.08373610487663233 + 6.301736487754372*eta + 9.03911386193751*pow(eta,2) + 4.91153188278086*pow(eta,3))*pow(1 + 72.64820846804257*eta,-1) + 
   (0.014963449678540705 + 0.008354571522567225*eta - 0.41723078020683*pow(eta,2) + 2.2007932082378785*pow(eta,3) - 4.245354787320365*pow(eta,4))*pow(S,2) + 
   (0.005706180633326235 - 0.15748500622007494*eta + 2.3477109912232845*pow(eta,2) - 11.413877195221694*pow(eta,3) + 17.033120593116756*pow(eta,4))*pow(S,3) + 
   (0.003890296981717687 - 0.15985471334551038*eta + 2.560312006077997*pow(eta,2) - 14.400920672743332*pow(eta,3) + 26.10406142567958*pow(eta,4))*pow(S,4) + 
   (0.005305988847210204 + 0.10869207132210629*eta - 2.4201307115268875*pow(eta,2) + 12.544899744864924*pow(eta,3) - 19.550600837316903*pow(eta,4))*pow(S,5) + 
   (0.002917248769788225 + 0.11851143848720952*eta - 2.6640023622893416*pow(eta,2) + 15.993378498844761*pow(eta,3) - 29.752144941054446*pow(eta,4))*pow(S,6);

   	return fit;
   }

/* Merger 22 frequency collocation points*/

static double IMRPhenomT_Merger_Freq_CP1_22(double eta, double S, double dchi, double delta) // theta=0.95
{

	double fit;

   fit = -0.3926039690467202*dchi*delta*(1 - 2.359180951434749*eta)*pow(eta,2) - 0.28551098014898896*dchi*delta*(1 - 3.414696100901444*eta)*S*pow(eta,2) + 0.003414004344822246*pow(dchi,2)*pow(eta,3) + 
   S*(0.05697014130854102 + 0.07170430925984912*eta - 0.9606499306623374*pow(eta,2) + 5.440955307244598*pow(eta,3) - 10.594319036394571*pow(eta,4)) + 
   (0.10030959768350425 + 44.56725135920024*eta + 163.96290948585087*pow(eta,2) - 143.05635831020462*pow(eta,3) + 393.8084861740473*pow(eta,4))*pow(1 + 436.6494065618*eta,-1) + 
   (0.021213606590798472 + 0.2148355967310081*eta - 2.7747405367196265*pow(eta,2) + 13.771088220299802*pow(eta,3) - 25.128755397215368*pow(eta,4))*pow(S,2) + 
   (-0.003645992092251503 + 0.2137524962844931*eta - 0.644979226062801*pow(eta,2) - 1.7314849842209137*pow(eta,3) + 5.573297392347478*pow(eta,4))*pow(S,3) + (0.029352214609533665 - 0.6020287633594307*eta + 7.014738679280164*pow(eta,2) - 36.027159248248296*pow(eta,3) + 63.42605850359639*pow(eta,4))*pow(S,4) + 
   (0.0356519646654399 - 0.5569780178251297*eta + 4.017784725334053*pow(eta,2) - 15.05881246593488*pow(eta,3) + 22.94821359434365*pow(eta,4))*pow(S,5);

   	return fit;
}

static double IMRPhenomT_PeakFrequency_22(double eta, double S, double dchi, double delta){

	double fit;

   fit = 0.27212130745330404 + 0.40972689759932074*eta - 0.0018392172960247433*eta*pow(dchi,2) + S*(0.09558832959428547 - 0.04834585264918328*eta - 0.15275173823699056*pow(eta,2)) - 3.4232387074402153*pow(eta,2) + 32.853772442252605*pow(eta,3) - 1.4976829186605336*dchi*delta*(1 - 4.775645585721007*eta)*pow(eta,3) - 
   0.9981117852179613*dchi*delta*(1 - 5.260098925354571*eta)*S*pow(eta,3) - 125.22505746137587*pow(eta,4) + 179.3797198714914*pow(eta,5) + (0.054391696704622204 - 0.1482682698299456*eta + 0.08938162810617255*pow(eta,2))*pow(S,2) + 
   (-0.020719540055375383 + 0.5090144456500953*eta - 1.5809441589349338*pow(eta,2))*pow(S,3) + (0.024240736699062685 - 0.09490089674418004*eta + 0.09518501714836035*pow(eta,2))*pow(S,4) + (0.09759303647532228 - 1.105520690228567*eta + 2.921271981239294*pow(eta,2))*pow(S,5);

   return fit;
}

/* RD 22 Frequency Coefficient fits */

static double IMRPhenomT_RD_Freq_D2_22(double eta, double S, double dchi, double delta){

   double fit;

   fit = 0.1598180460429256 + 0.19120040104567676*eta + (-0.012853620630980167 - 0.006532392920798404*eta)*S - 0.7733759581766899*pow(eta,2) + 0.18151402648790957*dchi*delta*(1 - 9.041198282315879*eta)*pow(eta,2) + 0.27147713896183995*dchi*delta*(1 - 5.653323210961101*eta)*S*pow(eta,2) - 
   0.01603489049446065*pow(dchi,2)*pow(eta,3) + (-0.046785083372074494 + 0.102759380109996*eta)*pow(S,2) + (0.0009883572415502464 - 0.050384608002279486*eta)*pow(S,3);

   return fit;
}

static double IMRPhenomT_RD_Freq_D3_22(double eta, double S, double dchi, double delta){

	double fit;

   fit = 2.6456463496860927 - 28.079375863863458*eta + 323.1691069138812*pow(eta,2) - 0.5040057675360762*dchi*delta*(1 + 21.786482297795278*eta)*pow(eta,2) + 1.561247215701216*dchi*delta*(1 - 1.7508069810164308*eta)*S*pow(eta,2) + S*(3.091917073632116 - 17.345283345692266*eta + 33.40735388809028*pow(eta,2)) - 
   1490.8128941604907*pow(eta,3) + 0.1619056474567525*pow(dchi,2)*pow(eta,3) + 2376.3257196613886*pow(eta,4) + (0.734022429223849 - 0.029342234233198747*eta - 9.281610698291932*pow(eta,2))*pow(S,2);

   	return fit;
}


/******************************* 22 AMPLITUDE FITS ***************************/

/* Inspiral 22 amplitude collocation points */


static double IMRPhenomT_Inspiral_Amp_CP1_22(double eta, double S, double dchi, double delta){

	double fit;

	fit = 0.00006480771730217768*eta*pow(dchi,2) - 0.3543965558027252*dchi*delta*(1 - 2.463526130684083*eta)*pow(eta,3) + 0.01879295038873938*dchi*delta*(1 - 5.236796607517272*eta)*S*pow(eta,3) + 
   S*(0.1472653807120573*eta - 1.9636752493349356*pow(eta,2) + 14.177521724634461*pow(eta,3) - 48.94620901701877*pow(eta,4) + 63.83730899015984*pow(eta,5)) + 
   eta*(0.8493442097893826 - 13.211067914003836*eta + 311.99021467938235*pow(eta,2) - 4731.025904601601*pow(eta,3) + 44821.93042533854*pow(eta,4) - 264474.1374080295*pow(eta,5) + 943246.2317701122*pow(eta,6) - 1.8588135904328802e6*pow(eta,7) + 1.5524778581809246e6*pow(eta,8)) + 
   (0.04902976057622393*eta - 1.0152511131279736*pow(eta,2) + 8.286289152216145*pow(eta,3) - 30.19775956110767*pow(eta,4) + 40.670065442751955*pow(eta,5))*pow(S,2) + 
   (0.04780630695082567*eta - 1.2177827888317065*pow(eta,2) + 11.505675146308567*pow(eta,3) - 46.733420749352135*pow(eta,4) + 68.40821782168776*pow(eta,5))*pow(S,3);

   	return fit;
}

static double IMRPhenomT_Inspiral_Amp_CP2_22(double eta, double S, double dchi, double delta){

	double fit;

	fit = 0.000100027278976821*eta*pow(dchi,2) - 0.7578403155712378*dchi*delta*(1 - 2.056456271350877*eta)*pow(eta,3) - 0.14126282637778914*dchi*delta*(1 - 2.5840771007494916*eta)*S*pow(eta,3) + S*(0.2331970217833686*eta - 1.5473968380422929*pow(eta,2) + 5.973401506474942*pow(eta,3) - 9.110484789161045*pow(eta,4)) + 
   eta*(0.9904613241626621 - 6.708006572605403*eta + 127.40270095439482*pow(eta,2) - 1723.355339710798*pow(eta,3) + 15430.10086310527*pow(eta,4) - 88744.26044058547*pow(eta,5) + 313650.01696201024*pow(eta,6) - 617887.8122937253*pow(eta,7) + 518220.9267888211*pow(eta,8)) + 
   (0.08934817374146888*eta - 0.8887847358339216*pow(eta,2) + 3.7233864099350784*pow(eta,3) - 5.814765403882651*pow(eta,4))*pow(S,2) + (0.04471990627820145*eta - 0.642458648615624*pow(eta,2) + 3.393481171493086*pow(eta,3) - 6.092083983738554*pow(eta,4))*pow(S,3);

   	return fit;
}

static double IMRPhenomT_Inspiral_Amp_CP3_22(double eta, double S, double dchi, double delta){

	double fit;

	fit = 0.0002459376633671657*eta*pow(dchi,2) - 0.8794763631110696*dchi*delta*(1 - 2.0751630535350096*eta)*pow(eta,3) - 0.3319387797134261*dchi*delta*(1 - 3.1838055629892184*eta)*S*pow(eta,3) + S*(0.23505507416274007*eta - 1.2449030421324767*pow(eta,2) + 4.315803728759738*pow(eta,3) - 6.384257606413192*pow(eta,4)) + 
   eta*(1.0208762064809185 - 3.3799457394243957*eta + 16.242639717123314*pow(eta,2) + 299.2297416582362*pow(eta,3) - 5913.920743907752*pow(eta,4) + 46388.231537995445*pow(eta,5) - 192261.0498470111*pow(eta,6) + 413750.14250475995*pow(eta,7) - 364403.84935539874*pow(eta,8)) + 
   (0.09630827896641526*eta - 0.7915321134872877*pow(eta,2) + 2.86907420250287*pow(eta,3) - 4.038995403653199*pow(eta,4))*pow(S,2) + (0.07395420485618898*eta - 1.0289224187583748*pow(eta,2) + 5.275845823734598*pow(eta,3) - 9.206158044409037*pow(eta,4))*pow(S,3);

   	return fit;
}

/* Merger 22 amplitude collocation points */

static double IMRPhenomT_Merger_Amp_CP1_22(double eta, double S, double dchi, double delta){

	double fit;

	fit = 0.0004059354652663733*eta*pow(dchi,2) - 0.9382383412276684*dchi*delta*(1 - 2.509151362054917*eta)*pow(eta,3) - 0.6560748977864668*dchi*delta*(1 - 3.426294113321932*eta)*S*pow(eta,3) + S*(0.23465398091766254*eta - 1.3398914201113978*pow(eta,2) + 5.9073801933446495*pow(eta,3) - 10.84221896204708*pow(eta,4)) + 
   eta*(1.2946032382158479 - 3.3343035556341816*eta + 91.6430240976277*pow(eta,2) - 1687.6195123629968*pow(eta,3) + 19726.50907350641*pow(eta,4) - 140798.18973779568*pow(eta,5) + 594095.3303894227*pow(eta,6) - 1.358657562562124e6*pow(eta,7) + 1.2958912179017465e6*pow(eta,8)) + 
   (0.03174875260265387*eta + 0.23082150180902375*pow(eta,2) - 1.9901867982613048*pow(eta,3) + 4.009389679757772*pow(eta,4))*pow(S,2) + (-0.04033221614773138*eta + 0.8426888041517518*pow(eta,2) - 4.742283264846479*pow(eta,3) + 9.059923021547936*pow(eta,4))*pow(S,3);

   	return fit;
}

static double IMRPhenomT_PeakAmp_22(double eta, double S, double dchi, double delta){

	double fit;

	fit = 0.0017885007700308166*eta*pow(dchi,2) - 0.5846280668038513*dchi*delta*(1 - 4.879882766464646*eta)*pow(eta,3) - 0.874161608112943*dchi*delta*(1 - 1.690095043235707*eta)*S*pow(eta,3) + S*(0.203557188205307*eta - 2.4368458739010563*pow(eta,2) + 12.206344183078137*pow(eta,3) - 23.417979354674692*pow(eta,4)) + 
   eta*(1.4701266133411792 - 1.387711607537906*eta + 25.641251409467607*pow(eta,2) - 186.013359336165*pow(eta,3) + 801.3039484150348*pow(eta,4) - 1893.8181854645718*pow(eta,5) + 1946.531703997353*pow(eta,6)) + 
   (-0.0018659293826992745*eta - 0.1888206507658455*pow(eta,2) + 1.4677324802664107*pow(eta,3) - 1.4019283350536489*pow(eta,4))*pow(S,2) + (-0.14699838946027494*eta + 2.6186847787143837*pow(eta,2) - 15.574381075605208*pow(eta,3) + 31.239292792717016*pow(eta,4))*pow(S,3);

   	return fit;
}

/*RD 22 Amplitude Coefficient Fits */

static double IMRPhenomT_RD_Amp_C3_22(double eta, double S){

	double fit;

   	fit = -0.48053994718185694 + 0.7023672141561462*eta + S*(-0.3597773028596323 + 1.4330280386796503*eta - 3.239121799338561*pow(eta,2)) - 0.1993836305574211*pow(eta,2) + (-0.2651107472061685 + 1.6433443489711386*eta - 2.757772023954491*pow(eta,2))*pow(S,2) + 
   (-0.01973537883495192 - 0.2410762147438714*eta + 2.7315015976869756*pow(eta,2))*pow(S,3);

   	return fit;
}

/////***************************************************************////////
//******************* HIGHER MODES FITS **************************/////
////**************************************************************////////

/* Frequency Fits */

static double IMRPhenomT_Merger_Freq_CP1_21(double eta, double S, double dchi, double delta){

   double fit;

   /*fit = 0.10101116560411222 - 0.0018259335648908525*eta - 0.021099940258397783*pow(eta,2) - 0.20805682832251896*dchi*delta*(1 - 2.941803484525313*eta)*pow(eta,2) - 0.34936387567469845*dchi*delta*(1 - 3.7803158001441313*eta)*S*pow(eta,2) + 
   S*(0.04417410903550135 - 0.043497437478068064*eta - 0.1259814374130449*pow(eta,2) + 0.32193272002103757*pow(eta,3)) + 1.5195639445091516*pow(eta,3) + 0.027556221044750074*pow(dchi,2)*pow(eta,3) - 3.734211998389987*pow(eta,4) + 
   (0.018056562830041197 + 0.18034946569991206*eta - 1.6734011535036983*pow(eta,2) + 3.234630062543167*pow(eta,3))*pow(S,2) + (0.04137468913747541 - 0.3769952802830424*eta + 1.220244052602933*pow(eta,2) - 1.1481788805520008*pow(eta,3))*pow(S,3) + 
   (0.04441609678703693 - 0.5527880503765283*eta + 1.9248039734409312*pow(eta,2) - 1.2987670735175931*pow(eta,3))*pow(S,4);*/

      fit = 0.10101718570562333 - 0.04992957057407379*eta + S*(0.04387200581505385 - 0.05535236748836419*eta + 0.01186527081077273*pow(eta,2)) + 0.8320472976764316*pow(eta,2) - 3.412564517258671*pow(eta,3) + 0.009795051484012513*pow(dchi,2)*pow(eta,3) - 
   0.7263372008413999*dchi*(1 - 15.69074025438453*pow(eta,2))*pow(eta,3) - 0.7339408066430174*dchi*S*(1 - 15.642722738277692*pow(eta,2))*pow(eta,3) + 5.462583026633245*pow(eta,4) + (0.02791012586433232 - 0.10972878398255334*eta + 0.18332510886777628*pow(eta,2))*pow(S,2) + 
   (0.050124128159139206 - 0.44122176836073457*eta + 1.0574966327609319*pow(eta,2))*pow(S,3) + (0.03798592538296406 - 0.3390611616798471*eta + 0.7453351213562277*pow(eta,2))*pow(S,4) + 0*delta;
      
    return fit;
}

static double IMRPhenomT_Merger_Freq_CP1_33(double eta, double S, double dchi, double delta){

   double fit;

      fit = 0.28925318007299916 - 0.11957063600442912*eta + 1.1589911273564353*pow(eta,2) - 0.5263424292081523*dchi*delta*(1 - 1.8018453981168454*eta)*pow(eta,2) - 1.509871876079703*dchi*delta*(1 - 4.654249984146787*eta)*S*pow(eta,2) - 1.5455016775261783*pow(eta,3) - 0.08410864613712594*pow(dchi,2)*pow(eta,3) + 
   S*(0.13053774656007028 - 0.1436725066031307*eta - 0.3988589708046106*pow(eta,2) + 1.6814308803706919*pow(eta,3)) + (0.06962116377895124 - 0.25606130995761806*eta + 2.809259385656226*pow(eta,2) - 10.655152499726746*pow(eta,3))*pow(S,2) + 
   (0.13941546691342796 - 1.4663373405231797*eta + 7.242184758809972*pow(eta,2) - 14.156137151704971*pow(eta,3))*pow(S,3) + (0.11741728765640543 - 0.9473051744486624*eta + 0.1981448121599412*pow(eta,2) + 8.2668429483671*pow(eta,3))*pow(S,4);

      return fit;
}

static double IMRPhenomT_Merger_Freq_CP1_44(double eta, double S, double dchi, double delta){

   double fit;

      fit = 0.3839588385106795 - 0.1433725509161493*eta + S*(0.16172235324877973 - 0.029039848526679315*eta - 0.5078615524759549*pow(eta,2)) + 1.549008100266779*pow(eta,2) - 2.782419405121844*pow(eta,3) - 6.300899329929631*dchi*delta*(1 - 3.4163840266407193*eta)*pow(eta,3) + 
   3.785702496330771e-7*dchi*delta*(1 - 7.683429533133828e6*eta)*S*pow(eta,3) - 0.006896021141579361*pow(dchi,2)*pow(eta,3) + (0.08448608215066597 - 0.04225970075410604*eta - 0.5752978186367078*pow(eta,2))*pow(S,2) + (0.19977344584848072 - 1.6500083087892263*eta + 3.6200398572688224*pow(eta,2))*pow(S,3) + 
   (0.17526054184412337 - 1.6567766073748962*eta + 3.920825944654624*pow(eta,2))*pow(S,4);

      return fit;
}

static double IMRPhenomT_Merger_Freq_CP1_55(double eta, double S, double dchi, double delta){

   double fit;

      fit = 0.49157926097800314 - 0.5109882379431707*eta + S*(0.23504290986387008 - 0.4663579238504208*eta + 0.6335674200647748*pow(eta,2)) + 5.1971222617574755*pow(eta,2) - 17.312514984111704*pow(eta,3) - 8.220653931136855*dchi*delta*(1 - 3.377023099309163*eta)*pow(eta,3) - 
   4.92651904204219*dchi*(1 - 4.0976182182508625*eta)*S*pow(eta,3) + 22.250087206692136*pow(eta,4) + (0.01818932761340069 + 1.2525110439050509*eta - 4.785925497931221*pow(eta,2))*pow(S,2) + (0.05473740403349404 + 0.204944664491078*eta - 1.5320283880617676*pow(eta,2))*pow(S,3) + 
   (0.19945251844606834 - 2.3208613314645397*eta + 6.551836762093784*pow(eta,2))*pow(S,4);

      return fit;
}

static double IMRPhenomT_PeakFrequency_21(double eta, double S, double dchi, double delta){

   double fit;

      fit = 0.17642087831932626 + 0.31718290537914057*eta + S*(0.03094734575888092 + 0.07319676429288274*eta - 0.4370939605469398*pow(eta,2)) - 2.2156624517537873*pow(eta,2) + 14.007580103948815*pow(eta,3) + 2.641340486447181*dchi*delta*(1 - 6.221406704917193*eta)*pow(eta,3) + 
   4.353108475005447*dchi*delta*(1 - 8.473808274993978*eta)*S*pow(eta,3) - 0.036084481180729745*pow(dchi,2)*pow(eta,3) - 26.085064860873068*pow(eta,4) + (0.017462707546942863 - 0.11463071986182106*eta + 0.2800463367551972*pow(eta,2))*pow(S,2) + 
   (0.033646159761323895 - 0.33812286198554814*eta + 0.9635090140454092*pow(eta,2))*pow(S,3) + (-0.0011779170821489254 + 0.18369948603536548*eta - 0.5978006007616697*pow(eta,2))*pow(S,4);

      return fit;
}

static double IMRPhenomT_PeakFrequency_33(double eta, double S, double dchi){

   double fit;

      fit = 0.42535721148121036 + 0.3085253521281911*eta + S*(0.1280277017287708 + 0.15271593642827125*eta - 0.9083681800119519*pow(eta,2)) + 0.9392741497311157*pow(eta,2) - 0.20785772397714286*dchi*(1 - 3.487216886252809*eta)*pow(eta,2) - 0.7863911902548658*dchi*(1 - 4.74913840513059*eta)*S*pow(eta,2) - 
   0.41376935975085416*pow(dchi,2)*pow(eta,3) + (0.09308538633777035 + 0.055164833113211194*eta - 1.1480525120934546*pow(eta,2))*pow(S,2) + (0.09945668702979882 - 0.5488068825374101*eta + 0.8675986447602085*pow(eta,2))*pow(S,3);

      return fit;
}

static double IMRPhenomT_PeakFrequency_44(double eta, double S, double dchi, double delta){

   double fit;

      fit = 0.5640094664638 + 0.3956446752668519*eta + S*(0.16597514208305744 + 0.38143981208933403*eta - 1.9002920053147696*pow(eta,2)) + 2.5091004914938675*pow(eta,2) - 7.403354368373608*pow(eta,3) - 5.257927939622048*dchi*delta*(1 - 5.385135507412752*eta)*pow(eta,3) + 
   1.1110261817411248e-7*dchi*delta*(1 - 1.0881293779054403e7*eta)*S*pow(eta,3) - 0.06378504432547372*pow(dchi,2)*pow(eta,3) + (0.08205749018653839 + 0.016449185328805776*eta - 0.509112344628105*pow(eta,2))*pow(S,2) + (0.13245468901111399 - 1.0716792675901017*eta + 2.631350201223915*pow(eta,2))*pow(S,3) + 
   (0.0798820256896006 - 0.6976704383121812*eta + 1.6808658698855679*pow(eta,2))*pow(S,4);

      return fit;
}

static double IMRPhenomT_PeakFrequency_55(double eta, double S, double dchi, double delta){

   double fit;

      fit = 0.7146297908371999 + 0.1421128402132339*eta + 7.659311331111322*pow(eta,2) + S*(0.29191927041842664 - 0.6512295551490094*eta + 1.021846701552054*pow(eta,2)) - 38.14301940776831*pow(eta,3) - 3.460574689440357*dchi*delta*(1 - 4.738903271021608*eta)*pow(eta,3) - 
   8.262749319140365*dchi*(1 - 4.1126856272636285*eta)*S*pow(eta,3) + 69.0208119373966*pow(eta,4) + (-0.17737667108149985 + 4.564503709808925*eta - 15.457705511019*pow(eta,2))*pow(S,2) + (-0.08755132408422435 + 1.8185807604067965*eta - 5.710975144545469*pow(eta,2))*pow(S,3) + 
   (0.4020378024101137 - 6.137619764177151*eta + 19.730459568297885*pow(eta,2))*pow(S,4);

      return fit;
}

/* RD lm Frequency Coefficient fits */

static double IMRPhenomT_RD_Freq_D2_21(double eta, double S, double dchi, double delta){

   double fit;

   fit = 0.1781545202005886 + 0.10906983816039043*eta + (-0.023905104384959013 + 0.1831847458257083*eta)*S - 0.5060291743082528*pow(eta,2) - 0.309304704734991*dchi*delta*(1 - 2.5742929128570724*eta)*pow(eta,2) + 2.1883684085193034*dchi*delta*(1 - 4.850311934387953*eta)*S*pow(eta,2) + 
   0.25978316114962485*pow(dchi,2)*pow(eta,3) + (-0.00955976176018747 - 0.18697585595061622*eta)*pow(S,2) + (0.04468930365659441 - 0.44170842157754653*eta)*pow(S,3);

      return fit;
}

static double IMRPhenomT_RD_Freq_D3_21(double eta, double S, double dchi, double delta){

   double fit;

   fit = 3.757258772469613 + 0.08380574896251641*eta + (3.634895503051922 - 8.174660936683596*eta)*S - 18.17832018576314*pow(eta,2) + 0.000043786944413372623*dchi*delta*(1 - 2.0348094407156347e6*eta)*pow(eta,2) + 60.2475724845033*dchi*delta*(1 - 6.868024913964549*eta)*S*pow(eta,2) + 
   23.08504961982195*pow(dchi,2)*pow(eta,3) + (3.7989582331116707 - 19.36029310028481*eta)*pow(S,2);

      return fit;
}

static double IMRPhenomT_RD_Freq_D2_33(double eta, double S, double dchi, double delta){

   double fit;

   fit = 0.16417885317959574 + 0.25804336274633655*eta + (-0.02961300365618534 - 0.006664043292875596*eta)*S - 0.9792038927762032*pow(eta,2) + 1.9953234313463062*dchi*delta*(1 - 5.45249062802972*eta)*pow(eta,2) + 6.3956780147142e-6*dchi*delta*(1 + 1.617071409433673e6*eta)*S*pow(eta,2) - 
   1.2001578283095737*pow(dchi,2)*pow(eta,3) + (-0.06346953330083285 + 0.12623926538220964*eta)*pow(S,2) + (-0.015173742568790456 + 0.016604992725771543*eta)*pow(S,3);
   
      return fit;
}

static double IMRPhenomT_RD_Freq_D3_33(double eta, double S, double dchi, double delta){

   double fit;

   fit = 2.0503935647397173 + 2.238118245943281*eta + (0.7121733508300451 - 0.397525795105057*eta)*S - 12.794117052655967*pow(eta,2) + 36.20571065481121*dchi*delta*(1 - 5.537022415039092*eta)*pow(eta,2) + 0.00018637091124974205*dchi*delta*(1 + 1.2667084084427443e6*eta)*S*pow(eta,2) - 
   21.894760631998928*pow(dchi,2)*pow(eta,3) + (0.4988930825119534 - 3.4004257158793045*eta)*pow(S,2) + (1.0586608433869105 - 5.625073332864818*eta)*pow(S,3);

      return fit;
}

static double IMRPhenomT_RD_Freq_D2_44(double eta, double S, double dchi, double delta){

   double fit;

   fit = 0.21336620664104916 - 0.20527614713716544*eta + (-0.057793617403743454 + 0.234794019739202*eta)*S - 0.5040007874429419*pow(eta,2) + 1.7980712659091223*dchi*delta*(1 - 4.3332243187779715*eta)*pow(eta,2) + 5.615398937364741*dchi*delta*(1 - 4.67655881619209*eta)*S*pow(eta,2) + 
   0.019754287494577062*pow(dchi,2)*pow(eta,3) + (-0.06870289106806035 + 0.18761585848765555*eta)*pow(S,2);

      return fit;
}

static double IMRPhenomT_RD_Freq_D3_44(double eta, double S, double dchi, double delta){

   double fit;
   fit = 3.824722911046124 - 19.434149952978917*eta + (1.1492216649186293 - 0.1794193390842707*eta)*S + 46.15218473526611*pow(eta,2) - 5.329235322993467e-6*dchi*delta*(1 + 3.600409478406777e6*eta)*pow(eta,2) + 108.44244729377574*dchi*delta*(1 - 4.948832996449642*eta)*S*pow(eta,2) + 
   0.8202829030161741*pow(dchi,2)*pow(eta,3) + (0.3851268017278357 - 1.4577127768796085*eta)*pow(S,2);

      return fit;
}

static double IMRPhenomT_RD_Freq_D2_55(double eta, double S, double dchi, double delta){

   double fit;

   fit =0.2143703929690296 - 0.26905171511199966*eta + (-0.057285673301351384 + 0.22530123030818466*eta)*S - 0.22128464791686953*pow(eta,2) + 1.2330723562386177*dchi*delta*(1 - 5.234362591591656*eta)*pow(eta,2) + 2.7651387378521104*dchi*delta*(1 - 4.529998650048839*eta)*S*pow(eta,2) - 
   0.02296167789737978*pow(dchi,2)*pow(eta,3) + (-0.06187017465583654 + 0.19469068079404978*eta)*pow(S,2);

      return fit;
}

static double IMRPhenomT_RD_Freq_D3_55(double eta, double S, double dchi, double delta){

   double fit;

   fit = 3.421296704767661 - 12.809224663506237*eta + (0.17979866256123875 + 3.9220602497543733*eta)*S + 22.579308761647468*pow(eta,2) + 38.938817471901075*dchi*delta*(1 - 5.439168816882736*eta)*pow(eta,2) + 102.01005958254325*dchi*delta*(1 - 4.94517313697764*eta)*S*pow(eta,2) + 
   8.156992533112646*pow(dchi,2)*pow(eta,3) + (-0.8276383989425874 + 5.653818482979737*eta)*pow(S,2);

      return fit;
}

//** Amplitude 21 Fits **//

static double IMRPhenomT_Inspiral_Amp_CP1_21(double eta, double S, double dchi, double delta){

   double fit;

      fit = -0.2457309233525402*dchi*(1 - 1.8588313811238013*eta)*pow(eta,3) + 0.007720682776232238*dchi*(1 - 14.5539282402835*eta)*S*pow(eta,3) + 0.00002718410442799091*pow(dchi,2)*pow(eta,3) + S*(-0.019371607120048675*delta*eta + 0.03368798661754525*delta*pow(eta,2) - 0.0347647962890128*delta*pow(eta,3)) + 
   delta*eta*(0.12222678288098383 - 2.152654527154567*eta + 34.53692688859637*pow(eta,2) - 317.45437636541044*pow(eta,3) + 1625.665271951051*pow(eta,4) - 4325.99209923682*pow(eta,5) + 4661.112076870376*pow(eta,6)) + 
   (-0.004130586129052499*delta*eta - 0.034242170459751614*delta*pow(eta,2) + 0.1845040639852827*delta*pow(eta,3))*pow(S,2) + (0.00023312994425693458*delta*eta - 0.006465524142621246*delta*pow(eta,2) + 0.02059744168116181*delta*pow(eta,3))*pow(S,3) + 
   (-0.010994253719930009*delta*eta + 0.1617856319808047*delta*pow(eta,2) - 0.5128238142456396*delta*pow(eta,3))*pow(S,4);

      return fit;
}

static double IMRPhenomT_Inspiral_Amp_CP2_21(double eta, double S, double dchi, double delta){

   double fit;

      fit = -0.5514762410690445*dchi*(1 - 1.6606901062713382*eta)*pow(eta,3) - 0.021703163232290525*dchi*(1 + 12.285199361388841*eta)*S*pow(eta,3) + 0.00027551818326783677*pow(dchi,2)*pow(eta,3) + S*(-0.013014289088905106*delta*eta - 0.14836733162360224*delta*pow(eta,2) + 0.3879852721571224*delta*pow(eta,3)) + 
   delta*eta*(0.15072063925506032 - 0.8093028329445506*eta + 6.206684655292913*pow(eta,2) - 24.88401414398108*pow(eta,3) + 38.250250718164864*pow(eta,4)) + (-0.025960288375186314*delta*eta + 0.09485066561654602*delta*pow(eta,2) - 0.12985415687429802*delta*pow(eta,3))*pow(S,2) + 
   (-0.031051316903933826*delta*eta + 0.29808639962599887*delta*pow(eta,2) - 0.7880170636799876*delta*pow(eta,3))*pow(S,3);

      return fit;
}

static double IMRPhenomT_Inspiral_Amp_CP3_21(double eta, double S, double dchi, double delta){

   double fit;

      fit = -0.6553365123485911*dchi*(1 - 1.5398595374318753*eta)*pow(eta,3) - 0.03414520050962973*dchi*(1 + 10.152070659598607*eta)*S*pow(eta,3) + 0.0003514981514078436*pow(dchi,2)*pow(eta,3) + S*(-0.02358276114828079*delta*eta - 0.06676889646672902*delta*pow(eta,2) + 0.10431702660244097*delta*pow(eta,3)) + 
   delta*eta*(0.16554873311231985 - 0.6991328198972108*eta + 4.8998331628863*pow(eta,2) - 17.811340834192666*pow(eta,3) + 25.04713555013603*pow(eta,4)) + (-0.03170047769861336*delta*eta + 0.12228560605854709*delta*pow(eta,2) - 0.2157318828663416*delta*pow(eta,3))*pow(S,2) + 
   (-0.02241156276655523*delta*eta + 0.1503547988268005*delta*pow(eta,2) - 0.32463957366468943*delta*pow(eta,3))*pow(S,3);

      return fit;
}

static double IMRPhenomT_Merger_Amp_CP1_21(double eta, double S, double dchi, double delta){

   double fit;

      fit = -0.9639235481813841*dchi*(1 - 0.0953303705707973*eta)*pow(eta,3) + 0.023311043707270836*dchi*(1 - 45.56758470014973*eta)*S*pow(eta,3) + 0.0006882033227649262*pow(dchi,2)*pow(eta,3) + S*(-0.07909979765641208*delta*eta - 0.1082323292489176*delta*pow(eta,2) - 0.059050965420919255*delta*pow(eta,3)) + 
   delta*eta*(0.2650216550129358 - 0.6689867318991676*eta + 6.322130556715582*pow(eta,2) - 25.30530828638204*pow(eta,3) + 37.57098989317853*pow(eta,4)) + (-0.048774914332679824*delta*eta + 0.058892171955548835*delta*pow(eta,2) + 0.04519042054728839*delta*pow(eta,3))*pow(S,2) + 
   (-0.04972189708933256*delta*eta + 0.32903712470645846*delta*pow(eta,2) - 0.7095763077249748*delta*pow(eta,3))*pow(S,3);

      return fit;
}

static double IMRPhenomT_PeakAmp_21(double eta, double S, double dchi, double delta){

   double fit;

   fit = -1.124757115880216*dchi*(1 + 3.9089731034256547*eta)*pow(eta,3) + 0.14171442436657175*dchi*S*pow(eta,3) - 7.996997960509883e-6*(1 + 12111.971615981536*delta)*pow(dchi,2)*pow(eta,3) + 
   delta*eta*(0.5940439865028524 - 2.6802250765521083*eta + 23.43295820742704*pow(eta,2) - 89.91427919476679*pow(eta,3) + 129.10731997830192*pow(eta,4)) + S*(-0.40438488955545776*delta*eta + 0.6359546829540189*delta*pow(eta,2) - 7.6174781238188*delta*pow(eta,3) + 20.156475820119724*delta*pow(eta,4)) + 
   (-0.04723336574759155*delta*eta + 0.18082387349024776*delta*pow(eta,2) + 1.7306679608818485*delta*pow(eta,3) - 8.236553093624009*delta*pow(eta,4))*pow(S,2) + 
   (-0.12534984288882925*delta*eta + 0.6131320823681302*delta*pow(eta,2) + 5.1648126976659885*delta*pow(eta,3) - 24.289576920541403*delta*pow(eta,4))*pow(S,3) + 
   (0.07112546745185065*delta*eta - 1.3149279454050955*delta*pow(eta,2) + 8.514263145733384*delta*pow(eta,3) - 14.271807407363035*delta*pow(eta,4))*pow(S,4);

   return fit;
}

static double IMRPhenomT_RD_Amp_C3_21(double eta, double S, double dchi){

   double fit;

      fit = -0.04334302376511826 - 0.17676752692299327*eta + 0.4505339209591958*pow(eta,2) + S*(-0.06491024396051823 - 1.1215130164808509*eta + 2.5523011435327345*pow(eta,2)) + (-0.28713100991035806 + 1.2391262662740283*eta - 2.7551841346664796*pow(eta,2))*pow(S,2) + 
   (-0.6910312848115802 + 5.91843541910692*eta - 13.892447750204266*pow(eta,2))*pow(S,3) + 0*dchi;

      return fit;
}

//** Amplitude 33 Fits **//

static double IMRPhenomT_Inspiral_Amp_CP1_33(double eta, double S, double dchi, double delta){

   double fit;

      fit = -0.00005000414942937797*delta*(1 - 3.0430401949925754*eta)*eta*pow(dchi,2) - 0.03836271211298855*dchi*(1 - 4.654767900586748*eta)*pow(eta,3) + 0.007041962008283751*dchi*(1 - 3.238646631077093*eta)*S*pow(eta,3) + 
   S*(0.0432725315235326*delta*eta - 0.3128744737439017*delta*pow(eta,2) + 0.7249180430447414*delta*pow(eta,3)) + delta*eta*(0.22272167356880285 - 3.217949139895537*eta + 45.52929729100423*pow(eta,2) - 379.70414120110206*pow(eta,3) + 1801.6287410802781*pow(eta,4) - 4505.468825419055*pow(eta,5) + 
      4606.517765490795*pow(eta,6)) + (0.015232248632190103*delta*eta - 0.15205944312376768*delta*pow(eta,2) + 0.38322848961855754*delta*pow(eta,3))*pow(S,2);

      return fit;
}

static double IMRPhenomT_Inspiral_Amp_CP2_33(double eta, double S, double dchi, double delta){

   double fit;

      fit = -0.0005485061120167634*delta*(1 - 5.249847868911592*eta)*eta*pow(dchi,2) - 0.13406080756104294*dchi*(1 - 4.791415116248203*eta)*pow(eta,3) - 0.025192101240368327*dchi*(1 - 5.557132409376257*eta)*S*pow(eta,3) + 
   S*(0.09903436069097878*delta*eta - 0.5266490647574258*delta*pow(eta,2) + 1.082646288776612*delta*pow(eta,3)) + delta*eta*(0.296454036493377 - 2.741774425959176*eta + 42.10341453030946*pow(eta,2) - 391.5079943491554*pow(eta,3) + 2084.7204836711294*pow(eta,4) - 5857.9923995429735*pow(eta,5) + 
      6724.299707693131*pow(eta,6)) + (0.04482378193895758*delta*eta - 0.34909482801592684*delta*pow(eta,2) + 0.798188874585321*delta*pow(eta,3))*pow(S,2);

      return fit;
}

static double IMRPhenomT_Inspiral_Amp_CP3_33(double eta, double S, double dchi, double delta){

   double fit;

      fit = -0.00014989518553589642*delta*(1 + 0.10284764229097754*eta)*eta*pow(dchi,2) - 0.16531803034216744*dchi*(1 - 4.9470029202324755*eta)*pow(eta,3) - 0.031723644862959394*dchi*(1 - 5.870965439700585*eta)*S*pow(eta,3) + 
   S*(0.11070499324391728*delta*eta - 0.5112660954416434*delta*pow(eta,2) + 0.9943348519498412*delta*pow(eta,3)) + delta*eta*
    (0.3081004973876298 - 1.4982270638091204*eta + 10.664775575232024*pow(eta,2) + 2.0410986773159214*pow(eta,3) - 472.97637767340444*pow(eta,4) + 2442.526427205543*pow(eta,5) - 3894.4435672165723*pow(eta,6)) + 
   (0.0509202148340841*delta*eta - 0.3395424984982766*delta*pow(eta,2) + 0.7165890644210602*delta*pow(eta,3))*pow(S,2);

      return fit;
}

static double IMRPhenomT_Merger_Amp_CP1_33(double eta, double S, double dchi, double delta){

   double fit;

      fit = -0.0009410965748168944*delta*(1 - 7.241296835003745*eta)*eta*pow(dchi,2) - 0.14546859303533213*dchi*(1 - 8.38731410081936*eta)*pow(eta,3) - 0.1398648328352838*dchi*(1 - 5.537001000046034*eta)*S*pow(eta,3) + 
   S*(0.11682330636277577*delta*eta - 0.2730553632132845*delta*pow(eta,2) + 0.5237293086635135*delta*pow(eta,3)) + delta*eta*
    (0.44186265057020313 - 1.2636555898615027*eta + 18.126195225020272*pow(eta,2) - 130.51981907268976*pow(eta,3) + 526.5238580108073*pow(eta,4) - 1078.532545666921*pow(eta,5) + 849.9216826816613*pow(eta,6)) + 
   (0.0510403796344345*delta*eta - 0.180522632008535*delta*pow(eta,2) + 0.4330636324653828*delta*pow(eta,3))*pow(S,2);
   
   return fit;
}

static double IMRPhenomT_PeakAmp_33(double eta, double S, double dchi, double delta){

   double fit;

      fit = -0.003288482386411718*delta*(1 - 8.612308762619447*eta)*eta*pow(dchi,2) + delta*eta*(0.5684405079702229 - 0.00028819674607128055*eta + 2.777740140752971*pow(eta,2) - 2.3599556709823535*pow(eta,3)) + 0.03887129318550153*dchi*(1 + 42.30525422235957*eta)*pow(eta,3) - 
   0.2051295687108511*dchi*(1 - 4.34985595987507*eta)*S*pow(eta,3) + S*(0.0652759726861487*delta*eta + 0.25561789058890033*delta*pow(eta,2) - 1.3134311480695775*delta*pow(eta,3)) + (0.04814607684462918*delta*eta - 0.3140983091545102*delta*pow(eta,2) + 1.1976699463228568*delta*pow(eta,3))*pow(S,2) + 
   (0.03619614547561679*delta*eta - 0.5532673160072701*delta*pow(eta,2) + 2.4943333040591695*delta*pow(eta,3))*pow(S,3);
      return fit;
}

static double IMRPhenomT_RD_Amp_C3_33(double eta, double S){

   double fit;

      fit = -0.28666660414434536 + 0.5669087275249756*eta + S*(-0.22961653919716726 + 0.7755862716197967*eta - 0.03726170050389395*pow(eta,2)) - 0.2969983864658452*pow(eta,2) + (-0.2177519810696989 + 1.5186886188134678*eta - 2.1091591639362255*pow(eta,2))*pow(S,2) + 
   (0.018605290436426794 + 0.8121676169377119*eta - 3.309654335397225*pow(eta,2))*pow(S,3);

      return fit;
}

//** Amplitude 44 Fits **//

static double IMRPhenomT_Inspiral_Amp_CP1_44(double eta, double S, double dchi, double delta){

   double fit;

      fit = S*(0.00929146984958081*eta - 0.058559157503356614*pow(eta,2) + 0.09641520260278541*pow(eta,3)) - 0.06256463263004813*dchi*delta*(1 - 4.724937783266512*eta)*pow(eta,3) - 0.01735529698327505*dchi*delta*(1 - 3.514044834242014*eta)*S*pow(eta,3) + 0.000842117844243168*pow(dchi,2)*pow(eta,3) + 
   eta*(0.0799508735514674 - 2.266747175041431*eta + 49.99562376971802*pow(eta,2) - 699.4551506778732*pow(eta,3) + 6096.872857701541*pow(eta,4) - 33243.794194712485*pow(eta,5) + 110236.14177804616*pow(eta,6) - 203224.30569500144*pow(eta,7) + 159685.76954854574*pow(eta,8)) + 
   (0.0036836428878512274*eta - 0.032181253212945134*pow(eta,2) + 0.06990383731270383*pow(eta,3))*pow(S,2) + (0.005305905044786687*eta - 0.05454070105726642*pow(eta,2) + 0.1328930616146293*pow(eta,3))*pow(S,3);

      return fabs(fit);
}

static double IMRPhenomT_Inspiral_Amp_CP2_44(double eta, double S, double dchi, double delta){

   double fit;

      fit = S*(0.032054779889996984*eta - 0.1824264213397133*pow(eta,2) + 0.2662860950846518*pow(eta,3)) - 0.09838860200524911*dchi*delta*(1 - 4.413878576399552*eta)*pow(eta,3) - 0.07541756416690493*dchi*delta*(1 - 4.896726739338081*eta)*S*pow(eta,3) + 0.00024755181872440586*pow(dchi,2)*pow(eta,3) + 
   eta*(0.10752001314377323 - 1.3996074805076077*eta + 17.290345408924*pow(eta,2) - 146.28994121129182*pow(eta,3) + 710.8477248404537*pow(eta,4) - 1819.0962884465648*pow(eta,5) + 1897.1460245953783*pow(eta,6)) + 
   (0.020209039503607196*eta - 0.1635522752682757*pow(eta,2) + 0.3379077937523624*pow(eta,3))*pow(S,2) + (0.016250056498330504*eta - 0.1599454341389429*pow(eta,2) + 0.38060091765599724*pow(eta,3))*pow(S,3);

      return fabs(fit);
}

static double IMRPhenomT_Inspiral_Amp_CP3_44(double eta, double S, double dchi, double delta){

   double fit;

      fit = S*(0.0375923390273927*eta - 0.19675674044979322*pow(eta,2) + 0.2524073874950236*pow(eta,3)) - 0.10103910572578918*dchi*delta*(1 - 4.5113969567894685*eta)*pow(eta,3) - 0.075429355757026*dchi*delta*(1 - 5.38318094443173*eta)*S*pow(eta,3) - 0.0001543369511547082*pow(dchi,2)*pow(eta,3) + 
   eta*(0.11684543226973083 - 1.1708344201904572*eta + 11.160637047449095*pow(eta,2) - 76.70545398788732*pow(eta,3) + 299.88284206545273*pow(eta,4) - 611.534557826681*pow(eta,5) + 503.71541521565484*pow(eta,6)) + 
   (0.025070408583957437*eta - 0.18759667520550588*pow(eta,2) + 0.36148626759006963*pow(eta,3))*pow(S,2) + (0.02203846168885738*eta - 0.20949146417573655*pow(eta,2) + 0.4885519836075034*pow(eta,3))*pow(S,3);

      return fabs(fit);
}

static double IMRPhenomT_Merger_Amp_CP1_44(double eta, double S, double dchi, double delta){

   double fit;

      fit = S*(0.058706429585806096*eta - 0.29646515787762634*pow(eta,2) + 0.3413797381980054*pow(eta,3)) + 0.002277607502351356*dchi*delta*(1 + 205.78624315268536*eta)*pow(eta,3) + 1.4842139797582348e-6*dchi*delta*(1 + 109390.82655130373*eta)*S*pow(eta,3) + 0.0034168773116717826*pow(dchi,2)*pow(eta,3) + 
   eta*(0.1883789790445804 - 1.3375818588191817*eta + 14.678591980033465*pow(eta,2) - 137.42245036191002*pow(eta,3) + 729.2611054727125*pow(eta,4) - 2042.3910621805978*pow(eta,5) + 2339.4396482545535*pow(eta,6)) + 
   (0.03693734108438023*eta - 0.2701634652682831*pow(eta,2) + 0.5277025188909843*pow(eta,3))*pow(S,2) + (0.01916289858506809*eta - 0.16587706623865955*pow(eta,2) + 0.3873216009596506*pow(eta,3))*pow(S,3);

      return fabs(fit);
}

static double IMRPhenomT_PeakAmp_44(double eta, double S, double dchi, double delta){

   double fit;

   fit = 0.697452842995687*dchi*delta*(1 - 1.5644622288381207*eta)*pow(eta,3) + 0.7491313799476855*dchi*delta*(1 - 5.51514415207437*eta)*S*pow(eta,3) + 0.03263766742842678*pow(dchi,2)*pow(eta,3) + S*(0.08013462545147897*eta - 0.707339986581501*pow(eta,2) + 1.4945536281037473*pow(eta,3)) + 
   eta*(0.27614097883794725 - 0.403452380875202*eta - 15.80475783619391*pow(eta,2) + 227.28867728765587*pow(eta,3) - 1523.2444219539561*pow(eta,4) + 4722.659771036674*pow(eta,5) - 5388.149395981192*pow(eta,6)) + 
   (0.0484478773571511*eta - 0.5421173150365266*pow(eta,2) + 1.5486181139304755*pow(eta,3))*pow(S,2) + (0.019255034163450358*eta - 0.18207194531823234*pow(eta,2) + 0.4812433162713078*pow(eta,3))*pow(S,3);
   return fabs(fit);
}

static double IMRPhenomT_RD_Amp_C3_44(double eta, double S){

   double fit;

      fit = -0.1772709159577312 - 0.3910604290424687*eta + S*(-0.14215203243769525 + 0.6136073658063063*eta + 0.11700379912379351*pow(eta,2)) + 5.876832797574524*pow(eta,2) + (-0.15523859963666756 + 0.5879889924742473*eta - 3.514395471691389*pow(eta,2))*pow(S,2) + 
   (-0.0829048220630192 - 1.965867892839485*eta + 11.364728644855896*pow(eta,2))*pow(S,3);

      return fit;
}

//** Amplitude 55 Fits **//

static double IMRPhenomT_Inspiral_Amp_CP1_55(double eta, double S, double dchi, double delta){

   double fit;

   fit = -0.0019775643769147514*dchi*(1 - 4.53924184281778*eta)*pow(eta,3) + 0.0014385048318273654*dchi*(1 - 3.856415079702643*eta)*S*pow(eta,3) + S*(0.004541994163205327*delta*eta - 0.032008447973076316*delta*pow(eta,2) + 0.06489393530989815*delta*pow(eta,3)) + 
   delta*eta*(0.028336688600531488 - 0.5409697917194701*eta + 5.981455840165866*pow(eta,2) - 35.58496309663864*pow(eta,3) + 106.06125067357719*pow(eta,4) - 124.75528935423806*pow(eta,5)) + 
   (0.0014791305997009444*delta*eta - 0.013564901743537111*delta*pow(eta,2) + 0.032142792215182195*delta*pow(eta,3))*pow(S,2);

      return fit;
}

static double IMRPhenomT_Inspiral_Amp_CP2_55(double eta, double S, double dchi, double delta){

   double fit;

      fit = -0.00003724952333242274*delta*(1 - 3.452222430510181*eta)*eta*pow(dchi,2) - 0.009011493135309705*dchi*(1 - 4.802896680414448*eta)*pow(eta,3) - 0.00013127526660987315*dchi*(1 - 30.606067223270347*eta)*S*pow(eta,3) + 
   S*(0.018145201665260808*delta*eta - 0.10875155354973318*delta*pow(eta,2) + 0.1967640499342343*delta*pow(eta,3)) + delta*eta*
    (0.05099973356701926 - 0.9909406291652298*eta + 16.002112346656087*pow(eta,2) - 151.48298211427934*pow(eta,3) + 798.2800157600177*pow(eta,4) - 2185.7904303138503*pow(eta,5) + 2423.2590529527615*pow(eta,6)) + 
   (0.007157848585528541*delta*eta - 0.04884915304244923*delta*pow(eta,2) + 0.09053190435686395*delta*pow(eta,3))*pow(S,2);

      return fit;
}

static double IMRPhenomT_Inspiral_Amp_CP3_55(double eta, double S, double dchi, double delta){

   double fit;

      fit = 0.000014325948759589005*delta*(1 - 0.6402417006774828*eta)*eta*pow(dchi,2) - 0.010476190606527443*dchi*(1 - 5.143499750443247*eta)*pow(eta,3) - 0.0030979537930086267*dchi*(1 - 5.795657988613435*eta)*S*pow(eta,3) + 
   S*(0.02269381988843006*delta*eta - 0.12335053938879872*delta*pow(eta,2) + 0.20474150752693748*delta*pow(eta,3)) + delta*eta*
    (0.04442248963378816 - 0.26245173832295265*eta + 0.03870693239289467*pow(eta,2) + 20.204218440428264*pow(eta,3) - 181.23307091031228*pow(eta,4) + 648.6147508591258*pow(eta,5) - 849.2143555875599*pow(eta,6)) + 
   (0.00913213583658649*delta*eta - 0.04919000144868387*delta*pow(eta,2) + 0.065986695477459*delta*pow(eta,3))*pow(S,2);

      return fit;
}

static double IMRPhenomT_Merger_Amp_CP1_55(double eta, double S, double dchi, double delta){

   double fit;

      fit = -0.00031681692142477853*delta*(1 - 6.678258498810139*eta)*eta*pow(dchi,2) + 0.04782012096371179*dchi*(1 - 2.4136112323522805*eta)*pow(eta,3) + 0.017136472806300446*dchi*(1 - 2.278716802319452*eta)*S*pow(eta,3) + 
   S*(0.035540705171105316*delta*eta - 0.13367790689815431*delta*pow(eta,2) + 0.1379312625561125*delta*pow(eta,3)) + delta*eta*
    (0.09741180071981874 - 1.125020791484383*eta + 19.71417205433866*pow(eta,2) - 206.5659386239305*pow(eta,3) + 1187.7241527779179*pow(eta,4) - 3511.371401268973*pow(eta,5) + 4163.584096800788*pow(eta,6)) + 
   (0.018044647105118435*delta*eta - 0.09834030851538031*delta*pow(eta,2) + 0.18292980258680058*delta*pow(eta,3))*pow(S,2);

      return fit;
}

static double IMRPhenomT_PeakAmp_55(double eta, double S, double dchi, double delta){

   double fit;

      fit = 0.29446699883224503*dchi*(1 - 2.8126528389736913*eta)*pow(eta,3) + 0.13166834017031467*dchi*(1 - 3.6911791457138365*eta)*S*pow(eta,3) + S*(0.04335058078657487*delta*eta - 0.21976500003781027*delta*pow(eta,2) + 0.1427254606254177*delta*pow(eta,3)) + 
   delta*eta*(0.25471102988397937 - 6.119431622115874*eta + 125.28192989146497*pow(eta,2) - 1339.55067240476*pow(eta,3) + 7641.25542069701*pow(eta,4) - 22186.340091384493*pow(eta,5) + 25846.606598287333*pow(eta,6)) + 
   (0.029015863810076155*delta*eta - 0.4063151087943421*delta*pow(eta,2) + 1.419210840554402*delta*pow(eta,3))*pow(S,2) + (0.01147033599820311*delta*eta - 0.28735230830842273*delta*pow(eta,2) + 1.1999844084222553*delta*pow(eta,3))*pow(S,3);

      return fit;
}

static double IMRPhenomT_RD_Amp_C3_55(double eta, double S, double dchi){

   double fit;

      fit = 0.01889156394866289 - 7.843569775936414*eta + S*(-0.11458748447064408 - 0.3369320850812222*eta - 0.022525692986479693*pow(eta,2)) + 73.19838355427139*pow(eta,2) - 170.9024182786024*pow(eta,3) + 76.38168535871085*dchi*(1 - 3.8805106289918205*eta)*pow(eta,3) + 
   42.628290542501134*dchi*(1 - 4.260931557685223*eta)*S*pow(eta,3) + (-0.18370843418415694 + 4.918601029356566*eta - 17.29835518657168*pow(eta,2))*pow(S,2);
      return fit;
}

/*************** RINGDOWN AND DAMPING QNM FREQUENCIES ***************/

/* FOR N=1, SAME CODE AS PHENOMXAS */

/*fRING*/

static double evaluate_QNMfit_fring21(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomT evaluate_QNMfit_fring21 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x2= finalDimlessSpin*finalDimlessSpin;
  double x3= x2*finalDimlessSpin;
  double x4= x2*x2;
  double x5= x3*x2;

  return_val = (0.059471695665734674 - 0.07585416297991414*finalDimlessSpin + 0.021967909664591865*x2 - 0.0018964744613388146*x3 + 0.001164879406179587*x4 - 0.0003387374454044957*x5)/(1 - 1.4437415542456158*finalDimlessSpin + 0.49246920313191234*x2);
  return return_val;
}

static double evaluate_QNMfit_fring33(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomT evaluate_QNMfit_fring33 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x2= finalDimlessSpin*finalDimlessSpin;
  double x3= x2*finalDimlessSpin;
  double x4= x2*x2;
  double x5= x3*x2;
  double x6= x3*x3;

  return_val = (0.09540436245212061 - 0.22799517865876945*finalDimlessSpin + 0.13402916709362475*x2 + 0.03343753057911253*x3 - 0.030848060170259615*x4 - 0.006756504382964637*x5 + 0.0027301732074159835*x6)/(1 - 2.7265947806178334*finalDimlessSpin + 2.144070539525238*x2 - 0.4706873667569393*x4 + 0.05321818246993958*x6);
  return return_val;
}

static double evaluate_QNMfit_fring44(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomT evaluate_QNMfit_fring44 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x2= finalDimlessSpin*finalDimlessSpin;
  double x3= x2*finalDimlessSpin;
  double x4= x2*x2;
  double x5= x3*x2;
  double x6= x3*x3;

  return_val = (0.1287821193485683 - 0.21224284094693793*finalDimlessSpin + 0.0710926778043916*x2 + 0.015487322972031054*x3 - 0.002795401084713644*x4 + 0.000045483523029172406*x5 + 0.00034775290179000503*x6)/(1 - 1.9931645124693607*finalDimlessSpin + 1.0593147376898773*x2 - 0.06378640753152783*x4);
  return return_val;
}

static double evaluate_QNMfit_fring55(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomT evaluate_QNMfit_fring55 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x= finalDimlessSpin;

  return_val = 0.16110773330909547 + (0.056600832610159385*x + 0.030041275213566483*pow(x,2) - 0.07522309632456432*pow(x,3) - 0.036341969761668556*pow(x,4) + 
      0.015617599737487714*pow(x,5) + 0.0062588909671250715*pow(x,6) + 0.004242111725892476*pow(x,7) + 0.0014913342466081074*pow(x,8))*
    pow(1 + 0.008923302958356548*x - 1.666395858912649*pow(x,2) + 0.6697719493836555*pow(x,4),-1);
  return return_val;
}

/* fDamp */

static double evaluate_QNMfit_fdamp21(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomT evaluate_QNMfit_fdamp21 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x2= finalDimlessSpin*finalDimlessSpin;
  double x3= x2*finalDimlessSpin;
  double x4= x2*x2;
  double x5= x3*x2;

  return_val = (2.0696914454467294 - 3.1358071947583093*finalDimlessSpin + 0.14456081596393977*x2 + 1.2194717985037946*x3 - 0.2947372598589144*x4 + 0.002943057145913646*x5)/(146.1779212636481 - 219.81790388304876*finalDimlessSpin + 17.7141194900164*x2 + 75.90115083917898*x3 - 18.975287709794745*x4);
  return return_val;
}

static double evaluate_QNMfit_fdamp33(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomT evaluate_QNMfit_fdamp33 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x2= finalDimlessSpin*finalDimlessSpin;
  double x3= x2*finalDimlessSpin;
  double x4= x2*x2;
  double x5= x3*x2;

  return_val = (0.014754148319335946 - 0.03124423610028678*finalDimlessSpin + 0.017192623913708124*x2 + 0.001034954865629645*x3 - 0.0015925124814622795*x4 - 0.0001414350555699256*x5)/(1 - 2.0963684630756894*finalDimlessSpin + 1.196809702382645*x2 - 0.09874113387889819*x4);
  return return_val;
}

static double evaluate_QNMfit_fdamp44(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomT evaluate_QNMfit_fdamp44 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x2= finalDimlessSpin*finalDimlessSpin;
  double x3= x2*finalDimlessSpin;
  double x4= x2*x2;
  double x5= x3*x2;
  double x6= x3*x3;

  return_val = (0.014986847152355699 - 0.01722587715950451*finalDimlessSpin - 0.0016734788189065538*x2 + 0.0002837322846047305*x3 + 0.002510528746148588*x4 + 0.00031983835498725354*x5 + 0.000812185411753066*x6)/(1 - 1.1350205970682399*finalDimlessSpin - 0.0500827971270845*x2 + 0.13983808071522857*x4 + 0.051876225199833995*x6);
  return return_val;
}

static double evaluate_QNMfit_fdamp55(double finalDimlessSpin){

  double return_val;

  if (fabs(finalDimlessSpin) > 1.0) {
    XLAL_ERROR(XLAL_EDOM, "PhenomT evaluate_QNMfit_fdamp55 function: |finalDimlessSpin| > 1.0 not supported");
  }

  double x= finalDimlessSpin;

  return_val = 0.015104212245401403 + pow(2 - 1.9485458003209648*x,-1)*(-0.0002946999837678157*x - 0.0024189312940399916*pow(x,2) + 0.0002099427928656942*pow(x,3) + 
      0.00258435043118687*pow(x,4) - 0.00020630579058983925*pow(x,5) - 0.004126708789254023*pow(x,6) + 0.0007950067180727237*pow(x,7) + 0.0027916616982894588*pow(x,8));
  return return_val;
}

/* For n=2, new fits of damping frequency needed for PhenomT */

static double evaluate_QNMfit_fdamp21n2(double finalDimlessSpin){

   double return_val;

   if (fabs(finalDimlessSpin) > 1.0) {
          XLAL_ERROR(XLAL_EDOM, "PhenomT evaluate_QNMfit_fdamp21n2 \
   function: |finalDimlessSpin| > 1.0 not supported");
   }

   double x = finalDimlessSpin;

   return_val = 0.04357957255256736 + pow(2 - 1.9092143068452778*x,-1)*(-0.0019991187832937543*x - 0.00397223929602004*pow(x,2) + 0.0027170335545048836*pow(x,3) - 
      0.003787735584625901*pow(x,4) + 0.003238742776891051*pow(x,5) + 0.0014093180629203572*pow(x,6));
   return return_val;
}

static double evaluate_QNMfit_fdamp33n2(double finalDimlessSpin){

   double return_val;

   if (fabs(finalDimlessSpin) > 1.0) {
          XLAL_ERROR(XLAL_EDOM, "PhenomT evaluate_QNMfit_fdamp33n2 \
   function: |finalDimlessSpin| > 1.0 not supported");
   }

   double x = finalDimlessSpin;

   return_val = 0.04478453069660422 + pow(2 - 1.9490123990107866*x,-1)*(-0.0027276947367212184*x - 0.005325382420460958*pow(x,2) + 0.0011090264831122598*pow(x,3) + 
      0.007374826520017088*pow(x,4) - 0.000513882756528504*pow(x,5) - 0.011798583916595289*pow(x,6) + 0.002064124132395282*pow(x,7) + 0.007865115260801307*pow(x,8));
   return return_val;
}

static double evaluate_QNMfit_fdamp44n2(double finalDimlessSpin){

   double return_val;

   if (fabs(finalDimlessSpin) > 1.0) {
          XLAL_ERROR(XLAL_EDOM, "PhenomT evaluate_QNMfit_fdamp33n2 \
   function: |finalDimlessSpin| > 1.0 not supported");
   }

   double x = finalDimlessSpin;

   return_val = 0.04526815749399381 + pow(2 - 1.9488568618006608*x,-1)*(-0.001778614725637923*x - 0.00645234041653255*pow(x,2) + 0.0008619365083550613*pow(x,3) + 
      0.0076173707591557305*pow(x,4) - 0.0005521040642302851*pow(x,5) - 0.012109903894557721*pow(x,6) + 0.0022638317039992374*pow(x,7) + 0.008166822924109219*pow(x,8));
   return return_val;
}

static double evaluate_QNMfit_fdamp55n2(double finalDimlessSpin){

   double return_val;

   if (fabs(finalDimlessSpin) > 1.0) {
          XLAL_ERROR(XLAL_EDOM, "PhenomT evaluate_QNMfit_fdamp33n2 \
   function: |finalDimlessSpin| > 1.0 not supported");
   }

   double x = finalDimlessSpin;

   return_val = 0.04550451880252191 + pow(2 - 1.948578997793657*x,-1)*(-0.0012254856066171247*x - 0.007001556265084966*pow(x,2) + 0.0006934110689396443*pow(x,3) + 
      0.007758785424957949*pow(x,4) - 0.0005928556371123753*pow(x,5) - 0.012383887808125627*pow(x,6) + 0.002365118383999583*pow(x,7) + 0.00838080791280435*pow(x,8));
   return return_val;
}

/*** PEAK TIMES OF HIGHER MODES ***/

static double IMRPhenomT_tshift_21(double eta, double S, double dchi)
{

   double fit = 11.67621653653603 - 73.94592135375544*eta + 617.5327332811615*pow(eta,2) + S*(0.2309485101131543 - 57.0459017581492*eta + 222.97200099809325*pow(eta,2)) - 2819.458362260437*pow(eta,3) - 681.9002621172333*dchi*(1 - 3.989581262513545*eta)*pow(eta,3) - 
   1440.639932639621*dchi*(1 - 4.206805719889809*eta)*S*pow(eta,3) + 42.39667266040204*pow(dchi,2)*pow(eta,3) + 4546.903391979042*pow(eta,4) + (2.2730487886808395 - 47.65323000340801*eta + 57.297549898351896*pow(eta,2))*pow(S,2) + 
   (-1.7237973456372406 + 21.732949566815307*eta - 187.71334824449366*pow(eta,2))*pow(S,3);
   return fit;
}

static double IMRPhenomT_tshift_33(double eta, double S)
{
   double fit = 6.047225180659371 - 63.50001473845436*eta + S*(1.26487072884024 + 9.789577125790505*eta - 18.51669370705306*pow(eta,2)) + 451.1074541600744*pow(eta,2) - 893.7051616506715*pow(eta,3) + (3.816104939071836 - 8.676597277291323*eta - 5.808122950219083*pow(eta,2))*pow(S,2) + (2.1374074060226045 - 1.2219912746034096*eta - 31.342471666791727*pow(eta,2))*pow(S,3);

   return fit;
}

static double IMRPhenomT_tshift_44(double eta, double S)
{
   double fit = S*(-5.203270014829841 + 181.27080583258746*eta - 1529.1896864534942*pow(eta,2) + 3705.463809339287*pow(eta,3)) + 
   pow(1 - 2.812531081541394*eta,-1)*(6.6472023470033585 - 98.64869153538237*eta + 1148.4724313577744*pow(eta,2) - 6720.146266369297*pow(eta,3) + 13400.05768313269*pow(eta,4)) + (6.9133369343740565 - 15.898281197030528*eta - 364.6027054334757*pow(eta,2) + 1362.455178365237*pow(eta,3))*pow(S,2) + 
   (23.15294108414908 - 333.2730725495644*eta + 1647.4278543557452*pow(eta,2) - 2702.213569022611*pow(eta,3))*pow(S,3);

   return fit;
}

static double IMRPhenomT_tshift_55(double eta, double S)
{
   double fit = -0.3189869259194407 + 153.08687719603935*eta - 1376.0895569730135*pow(eta,2) + S*(10.62118364975699 - 128.69299551679973*eta + 401.7008544741773*pow(eta,2)) + 3511.1779067574766*pow(eta,3) + (10.036441054322665 - 75.42317272972994*eta + 180.54334490779055*pow(eta,2))*pow(S,2) + (6.658297274982426 - 35.88874981710895*eta + 33.86225466072782*pow(eta,2))*pow(S,3);

   return fit;
}


/*************** RINGDOWN AND DAMPING QNM FREQUENCIES ***************/

/* FOR N=1, SAME CODE AS PHENOMXAS */

/*fRING*/

static double evaluate_QNMfit_fring22(double finalDimlessSpin){

   double return_val;

   if (fabs(finalDimlessSpin) > 1.0) {
         XLAL_ERROR(XLAL_EDOM, "PhenomT evaluate_QNMfit_fring22 \
   function: |finalDimlessSpin| > 1.0 not supported");
   }

   double x2= finalDimlessSpin*finalDimlessSpin;
   double x3= x2*finalDimlessSpin;
   double x4= x2*x2;
   double x5= x3*x2;
   double x6= x3*x3;
   double x7= x4*x3;

   return_val = (0.05947169566573468 - \
   0.14989771215394762*finalDimlessSpin + 0.09535606290986028*x2 + \
   0.02260924869042963*x3 - 0.02501704155363241*x4 - \
   0.005852438240997211*x5 + 0.0027489038393367993*x6 + \
   0.0005821983163192694*x7)/(1 - 2.8570126619966296*finalDimlessSpin + \
   2.373335413978394*x2 - 0.6036964688511505*x4 + \
   0.0873798215084077*x6);
   return return_val;
}

/* fDamp */

static double evaluate_QNMfit_fdamp22(double finalDimlessSpin){

   double return_val;

   if (fabs(finalDimlessSpin) > 1.0) {
          XLAL_ERROR(XLAL_EDOM, "PhenomX evaluate_QNMfit_fdamp22 \
   function: |finalDimlessSpin| > 1.0 not supported");
   }

   double x2= finalDimlessSpin*finalDimlessSpin;
   double x3= x2*finalDimlessSpin;
   double x4= x2*x2;
   double x5= x3*x2;
   double x6= x3*x3;

   return_val = (0.014158792290965177 - \
   0.036989395871554566*finalDimlessSpin + 0.026822526296575368*x2 + \
   0.0008490933750566702*x3 - 0.004843996907020524*x4 - \
   0.00014745235759327472*x5 + 0.0001504546201236794*x6)/(1 - \
   2.5900842798681376*finalDimlessSpin + 1.8952576220623967*x2 - \
   0.31416610693042507*x4 + 0.009002719412204133*x6);
   return return_val;
}

/* For n=2, new fits of damping frequency needed for PhenomT */

static double evaluate_QNMfit_fdamp22n2(double finalDimlessSpin){

   double return_val;

   if (fabs(finalDimlessSpin) > 1.0) {
          XLAL_ERROR(XLAL_EDOM, "PhenomT evaluate_QNMfit_fdamp22n2 \
   function: |finalDimlessSpin| > 1.0 not supported");
   }

   double x = finalDimlessSpin;

   return_val = 0.043611742588188715 + pow(2 - 1.9477781396815619*x,-1)*(-0.004016191313442792*x - 0.0027646155943395426*pow(x,2) + 0.001141927763953028*pow(x,3) + 
      0.007938320030300492*pow(x,4) - 0.0008263166671238823*pow(x,5) - 0.014025760257115768*pow(x,6) + 0.001792158578158245*pow(x,7) + 0.008824138122361842*pow(x,8));
   return return_val;
}


