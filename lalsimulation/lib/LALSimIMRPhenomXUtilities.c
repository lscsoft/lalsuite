/*
 *  Copyright (C) 2018 Geraint Pratten
 *
 *  This code adapts functions from:
 *    LALSimIMRPhenomP.c
 *    LALSimIMRPhenomD.c
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

/**
 * \author Geraint Pratten
 *
 * \file
 *
 * \brief Utility functions for IMRPhenomX framework, arXiv:2001.11412
 *
 */

#include <gsl/gsl_math.h>

#include "LALSimIMRPhenomXUtilities.h"

/* ******************** MECO, ISCO, etc ******************** */

/*
 * Phenomenological fit to hybrid minimum energy circular orbit (MECO) function.
 * Uses 3.5PN hybridised with test-particle limit.
 * Reference: M Cabero et al, PRD, 95, 064016, (2017), arXiv:1602.03134
 */
REAL8 XLALSimIMRPhenomXfMECO(REAL8 eta, REAL8 chi1L, REAL8 chi2L) {

  REAL8 eta2  = (eta*eta);
  REAL8 eta3  = (eta2*eta);
  REAL8 eta4  = (eta3*eta);

  REAL8 delta = sqrt(1.0 - 4.0*eta);

  REAL8 S     = XLALSimIMRPhenomXchiPNHat(eta,chi1L,chi2L);
  REAL8 S2    = (S*S);
  REAL8 S3    = (S2*S);
  //REAL8 S4    = (S3*S);

  REAL8 noSpin, eqSpin, uneqSpin;

  REAL8 dchi  = chi1L - chi2L;
  REAL8 dchi2 = (dchi*dchi);


  noSpin = (0.018744340279608845 + 0.0077903147004616865*eta + 0.003940354686136861*eta2 - 0.00006693930988501673*eta3)/(1. - 0.10423384680638834*eta);

  eqSpin = (S*(0.00027180386951683135 - 0.00002585252361022052*S + eta4*(-0.0006807631931297156 + 0.022386313074011715*S - 0.0230825153005985*S2) + eta2*(0.00036556167661117023 - 0.000010021140796150737*S - 0.00038216081981505285*S2) + eta*(0.00024422562796266645 - 0.00001049013062611254*S - 0.00035182990586857726*S2) + eta3*(-0.0005418851224505745 + 0.000030679548774047616*S + 4.038390455349854e-6*S2) - 0.00007547517256664526*S2))/(0.026666543809890402 + (-0.014590539285641243 - 0.012429476486138982*eta + 1.4861197211952053*eta4 + 0.025066696514373803*eta2 + 0.005146809717492324*eta3)*S + (-0.0058684526275074025 - 0.02876774751921441*eta - 2.551566872093786*eta4 - 0.019641378027236502*eta2 - 0.001956646166089053*eta3)*S2 + (0.003507640638496499 + 0.014176504653145768*eta + 1.*eta4 + 0.012622225233586283*eta2 - 0.00767768214056772*eta3)*S3);

  uneqSpin = dchi2*(0.00034375176678815234 + 0.000016343732281057392*eta)*eta2 + dchi*delta*eta*(0.08064665214195679*eta2 + eta*(-0.028476219509487793 - 0.005746537021035632*S) - 0.0011713735642446144*S);

  return (noSpin + eqSpin + uneqSpin);
}

/*
 * Fitting function for hybrid minimum energy circular orbit (MECO) function
 */
REAL8 XLALSimIMRPhenomXfISCO(REAL8 chif) {

  REAL8 OmegaISCO, rISCO;
  REAL8 rISCOsq, rISCO3o2;
  REAL8 Z1, Z2;

  Z1 = 1.0 + cbrt( (1.0 - chif*chif) ) * ( cbrt(1 + chif) + cbrt(1 - chif) );
  if(Z1>3) Z1=3.; //Finite precission may give Z1>3, but this can not happen.
  Z2 = sqrt(3.0*chif*chif + Z1*Z1);

  rISCO    = 3.0 + Z2 - XLALSimIMRPhenomXsign(chif)*sqrt( (3 - Z1) * (3 + Z1 + 2*Z2) );
  rISCOsq  = sqrt(rISCO);
  rISCO3o2 = rISCOsq * rISCOsq * rISCOsq;

  OmegaISCO = 1.0 / ( rISCO3o2 + chif);

  return OmegaISCO / LAL_PI;

}

/* ******************** FINAL STATE ********************

    These functions are here as we exposed them violate
    the XLAL wrappers.
*/

/*
 * Energy Radiated: X. Jimenez-Forteza et al, PRD, 95, 064024, (2017), arXiv:1611.00332
 */
REAL8 XLALSimIMRPhenomXErad2017(REAL8 eta, REAL8 chi1L, REAL8 chi2L) {

  REAL8 delta  = sqrt(1.0 - 4.0*eta);

  REAL8 eta2   =  eta*eta;
  REAL8 eta3   = eta2*eta;
  REAL8 eta4   = eta3*eta;

  REAL8 S      = XLALSimIMRPhenomXSTotR(eta,chi1L,chi2L);
  REAL8 S2     =  S*S;
  REAL8 S3     = S2*S;

  REAL8 dchi   = chi1L - chi2L;
  REAL8 dchi2  = dchi*dchi;

  REAL8 noSpin, eqSpin, uneqSpin;

  noSpin = 0.057190958417936644*eta + 0.5609904135313374*eta2 - 0.84667563764404*eta3 + 3.145145224278187*eta4;

  /* Because of the way this is written, we need to subtract the noSpin term */
  eqSpin = ((0.057190958417936644*eta + 0.5609904135313374*eta2 - 0.84667563764404*eta3 + 3.145145224278187*eta4)*
    (    1
      + (-0.13084389181783257 - 1.1387311580238488*eta + 5.49074464410971*eta2)*S
      + (-0.17762802148331427 + 2.176667900182948*eta2)*S2
      + (-0.6320191645391563 + 4.952698546796005*eta - 10.023747993978121*eta2)*S3))
      / (1 + (-0.9919475346968611 + 0.367620218664352*eta + 4.274567337924067*eta2)*S);

  eqSpin = eqSpin - noSpin;

  uneqSpin =  - 0.09803730445895877*dchi*delta*(1 - 3.2283713377939134*eta)*eta2
              + 0.01118530335431078*dchi2*eta3
              - 0.01978238971523653*dchi*delta*(1 - 4.91667749015812*eta)*eta*S;

  return (noSpin + eqSpin + uneqSpin);

}

/*
 * Final Mass = 1 - Energy Radiated,  X. Jimenez-Forteza et al, PRD, 95, 064024, (2017), arXiv:1611.00332
 */
REAL8 XLALSimIMRPhenomXFinalMass2017(REAL8 eta, REAL8 chi1L, REAL8 chi2L) {

  REAL8 delta  = sqrt(1.0 - 4.0*eta);
  REAL8 eta2   =  eta*eta;
  REAL8 eta3   = eta2*eta;
  REAL8 eta4   = eta3*eta;

  REAL8 S      = XLALSimIMRPhenomXSTotR(eta,chi1L,chi2L);
  REAL8 S2     =  S*S;
  REAL8 S3     = S2*S;

  REAL8 dchi   = chi1L - chi2L;
  REAL8 dchi2  = dchi*dchi;

  REAL8 noSpin, eqSpin, uneqSpin;

  noSpin = 0.057190958417936644*eta + 0.5609904135313374*eta2 - 0.84667563764404*eta3 + 3.145145224278187*eta4;

  /* Because of the way this is written, we need to subtract the noSpin term */
  eqSpin = ((0.057190958417936644*eta + 0.5609904135313374*eta2 - 0.84667563764404*eta3 + 3.145145224278187*eta4)*
    (    1
      + (-0.13084389181783257 - 1.1387311580238488*eta + 5.49074464410971*eta2)*S
      + (-0.17762802148331427 + 2.176667900182948*eta2)*S2
      + (-0.6320191645391563 + 4.952698546796005*eta - 10.023747993978121*eta2)*S3))
      / (1 + (-0.9919475346968611 + 0.367620218664352*eta + 4.274567337924067*eta2)*S);

  eqSpin = eqSpin - noSpin;

  uneqSpin =  - 0.09803730445895877*dchi*delta*(1 - 3.2283713377939134*eta)*eta2
              + 0.01118530335431078*dchi2*eta3
              - 0.01978238971523653*dchi*delta*(1 - 4.91667749015812*eta)*eta*S;

  /* Mfinal = 1 - Erad, assuming that M = m1 + m2 = 1 */
  return (1.0 - (noSpin + eqSpin + uneqSpin));


}

/*
 * Final Dimensionless Spin,  X. Jimenez-Forteza et al, PRD, 95, 064024, (2017), arXiv:1611.00332
 */
REAL8 XLALSimIMRPhenomXFinalSpin2017(REAL8 eta, REAL8 chi1L, REAL8 chi2L) {

	REAL8 delta  = sqrt(1.0 - 4.0*eta);
	REAL8 m1     = 0.5 * (1.0 + delta);
	REAL8 m2     = 0.5 * (1.0 - delta);
	REAL8 m1Sq   = m1*m1;
	REAL8 m2Sq   = m2*m2;

  REAL8 eta2   = eta*eta;
  REAL8 eta3   = eta2*eta;

	//REAL8 S  = (m1Sq * chi1L + m2Sq * chi2L) / (m1Sq + m2Sq);
  REAL8 S  = XLALSimIMRPhenomXSTotR(eta,chi1L,chi2L);
  REAL8 S2 =  S*S;
  REAL8 S3 = S2*S;

  REAL8 dchi  = chi1L - chi2L;
  REAL8 dchi2 = dchi*dchi;

  REAL8 noSpin, eqSpin, uneqSpin;

  noSpin = (3.4641016151377544*eta + 20.0830030082033*eta2 - 12.333573402277912*eta3)/(1 + 7.2388440419467335*eta);

  eqSpin = (m1Sq + m2Sq)*S
  + ((-0.8561951310209386*eta - 0.09939065676370885*eta2 + 1.668810429851045*eta3)*S
  + (0.5881660363307388*eta - 2.149269067519131*eta2 + 3.4768263932898678*eta3)*S2
  + (0.142443244743048*eta - 0.9598353840147513*eta2 + 1.9595643107593743*eta3)*S3)
  / (1 + (-0.9142232693081653 + 2.3191363426522633*eta - 9.710576749140989*eta3)*S);

  uneqSpin = 0.3223660562764661*dchi*delta*(1 + 9.332575956437443*eta)*eta2           /* Linear in spin difference                */
  - 0.059808322561702126*dchi2*eta3                                                   /* Quadratic in spin difference             */
  + 2.3170397514509933*dchi*delta*(1 - 3.2624649875884852*eta)*eta3*S;                 /* Mixed spin difference + total spin term  */


	return (noSpin + eqSpin + uneqSpin);
}

/**
 * Wrapper for the final spin in generically precessing binary black holes
 */
REAL8 XLALSimIMRPhenomXPrecessingFinalSpin2017(
  const REAL8 eta,           /**< Symmetric mass ratio                    */
  const REAL8 chi1L,         /**< Aligned spin of BH 1                    */
  const REAL8 chi2L,         /**< Aligned spin of BH 2                    */
  const REAL8 chi_inplane    /**< Effective precessions spin parameter, see Section IV D of arXiv:XXXX.YYYY    */
)
{
  REAL8 m1 = 0.5 * (1.0 + sqrt(1 - 4.0*eta));
  REAL8 m2 = 0.5 * (1.0 - sqrt(1 - 4.0*eta));
  REAL8 M  = m1 + m2;

  if (eta > 0.25) IMRPhenomX_InternalNudge(eta, 0.25, 1e-6);

  REAL8 af_parallel, q_factor;
  if (m1 >= m2) {
    q_factor    = m1 / M;
    af_parallel = XLALSimIMRPhenomXFinalSpin2017(eta, chi1L, chi2L);
  }
  else {
    q_factor    = m2 / M;
    af_parallel = XLALSimIMRPhenomXFinalSpin2017(eta, chi2L, chi1L);
  }

  REAL8 Sperp   = chi_inplane * q_factor * q_factor;
  REAL8 af      = copysign(1.0, af_parallel) * sqrt(Sperp*Sperp + af_parallel*af_parallel);
  return af;
}

/* ******************** SPIN PARAMETERIZATIONS ******************** */
/*
 * PN reduced spin parameter
 */
REAL8 XLALSimIMRPhenomXchiPN(REAL8 eta, REAL8 chi1L, REAL8 chi2L) {
	// Convention m1 >= m2 and chi1 is the spin on m1
	REAL8 delta    = sqrt(1.0 - 4.0*eta);
  REAL8 mm1      = 0.5*(1+delta);
  REAL8 mm2      = 0.5*(1-delta);
  REAL8 chi_eff  = (mm1*chi1L + mm2*chi2L);

	return chi_eff - (38.0/113.0)*eta*(chi1L + chi2L);
}

/*
 * Normalised PN reduced spin parameter
 */
REAL8 XLALSimIMRPhenomXchiPNHat(REAL8 eta, REAL8 chi1L, REAL8 chi2L) {
  // Convention m1 >= m2 and chi1 is the spin on m1
	REAL8 delta    = sqrt(1.0 - 4.0*eta);
  REAL8 mm1      = 0.5*(1.0 + delta);
  REAL8 mm2      = 0.5*(1.0 - delta);
  REAL8 chi_eff  = (mm1*chi1L + mm2*chi2L);

	return (chi_eff - (38.0/113.0)*eta*(chi1L + chi2L) ) / (1.0 - (76.0*eta/113.0));
}

/*
 * Effective aligned spin parameter
 */
REAL8 XLALSimIMRPhenomXchiEff(REAL8 eta, REAL8 chi1L, REAL8 chi2L) {
	// Convention m1 >= m2 and chi1 is the spin on m1
	REAL8 delta = sqrt(1.0 - 4.0*eta);
	REAL8 mm1 = 0.5*(1+delta);
  REAL8 mm2 = 0.5*(1-delta);

  return (mm1*chi1L + mm2*chi2L);
}

/*
 * Total spin normalised to [-1,1]
 */
REAL8 XLALSimIMRPhenomXSTotR(REAL8 eta, REAL8 chi1L, REAL8 chi2L) {
	// Convention m1 >= m2 and chi1z is the spin projected along Lz on m1
	REAL8 delta = sqrt(1.0 - 4.0*eta);
	REAL8 m1    = 0.5*(1 + delta);
	REAL8 m2    = 0.5*(1 - delta);
	REAL8 m1s   = m1*m1;
	REAL8 m2s   = m2*m2;

	return ((m1s * chi1L + m2s * chi2L) / (m1s + m2s));
}

/*
 * Spin difference
 */
REAL8 XLALSimIMRPhenomXdchi(REAL8 chi1L, REAL8 chi2L) {

  return chi1L - chi2L;
}

/* ******************** FREQUENCY CONVERSIONS ******************** */
/*
 * Convert from geometric frequency to Hz
 */
REAL8 XLALSimIMRPhenomXUtilsMftoHz(
    REAL8 Mf,       /**< Geometric frequency */
    REAL8 Mtot_Msun /**< Total mass in solar masses */
)
{
    return Mf / (LAL_MTSUN_SI * Mtot_Msun);
}

/*
 * Convert from frequency in Hz to geometric frequency
 */
REAL8 XLALSimIMRPhenomXUtilsHztoMf(
    REAL8 fHz,      /**< Frequency in Hz */
    REAL8 Mtot_Msun /**< Total mass in solar masses */
)
{
    return fHz * (LAL_MTSUN_SI * Mtot_Msun);
}

/*
 * We apply a linear time and phase shift to ~ align peak
 * LinShift = (PNLina[\[Eta],\[Chi]1,\[Chi]2] + \[Pi] + f PNLinb[\[Eta],\[Chi]1,\[Chi]2]);
 * Linear time and phase shift: a + b*f
 */
REAL8 XLALSimIMRPhenomXLina(
    REAL8 eta,       /**< Geometric frequency */
    REAL8 S,         /**< Total mass in solar masses */
    REAL8 dchi,      /**< Total mass in solar masses */
    REAL8 delta      /**< Total mass in solar masses */
)
{

  double eta2 = eta*eta;
  double eta3 = eta2*eta;

  double S2 = S*S;
  double S3 = S2*S;
  double S4 = S3*S;
  double S5 = S4*S;

  double noSpin, eqSpin, uneqSpin;

  noSpin = (1.0691011796680957e7 + 1.185905809135487e6*eta)/(1. + 30289.726104019595*eta);

  eqSpin = (3590.4466629551066 - 43200.79177912654*eta +
              200094.57177252226*eta2 - 319307.42983118095*eta3)*S +
              eta*(-1108.7615587239336 + 25622.545977741574*eta -
              83180.15680637326*eta2)*S3 + (379.0250508368122 -
              1355.868129015304*eta)*eta2*S4 + (-23306.844979169644 +
              91977.08490230633*eta)*eta2*S5 + (-204.5064259069199 + 1046.9991525832384*eta - 5906.920781540527*eta3)*S2;

  uneqSpin = 44.87175399132858*dchi*delta*eta;

  return (noSpin + eqSpin + uneqSpin) + LAL_PI;

}

// this is a fit of the time-difference between t_peak of strain and t_peak of psi4
// needed to align in time our waveforms, which are calibrated to psi4

REAL8 XLALSimIMRPhenomXPsi4ToStrain(double eta, double S, double dchi) {
    double eta2,eta3,eta4,S2,S3,S4;
    eta2 = pow(eta,2);
    eta3 = pow(eta,3);
    eta4 = pow(eta,4);
    S2 = pow(S,2);
    S3 = pow(S,3);
    S4 = pow(S,4);
    double noSpin = 13.39320482758057 - 175.42481512989315*eta + 2097.425116152503*eta2 - 9862.84178637907*eta3 + 16026.897939722587*eta4;
    double eqSpin = (4.7895602776763 - 163.04871764530466*eta + 609.5575850476959*eta2)*S + (1.3934428041390161 - 97.51812681228478*eta + 376.9200932531847*eta2)*S2 + (15.649521097877374 + 137.33317057388916*eta - 755.9566456906406*eta2)*S3 + (13.097315867845788 + 149.30405703643288*eta - 764.5242164872267*eta2)*S4;
    double uneqSpin = 105.37711654943146*dchi*sqrt(1. - 4.*eta)*eta2;
    return( noSpin + eqSpin + uneqSpin);


}


REAL8 XLALSimIMRPhenomXLinb(
  REAL8 eta,       /**< Geometric frequency */
  REAL8 S,         /**< Total mass in solar masses */
  REAL8 dchi,      /**< Total mass in solar masses */
  REAL8 delta      /**< Total mass in solar masses */
)
{
    double eta2 = eta*eta;
    double eta3= eta2*eta, eta4= eta3*eta, eta5=eta3*eta2, eta6=eta4*eta2;
    double S2 = S*S;
    double S3 = S2*S;
    double S4 = S3*S;
    double noSpin = 3155.1635543201924 + 1257.9949740608242*eta - 32243.28428870599*eta2 + 347213.65466875216*eta3 - 1.9223851649491738e6*eta4 + 5.3035911346921865e6*eta5 - 5.789128656876938e6*eta6;
    double eqSpin = (-24.181508118588667 + 115.49264174560281*eta - 380.19778216022763*eta2)*S + (24.72585609641552 - 328.3762360751952*eta + 725.6024119989094*eta2)*S2 + (23.404604124552 - 646.3410199799737*eta + 1941.8836639529036*eta2)*S3 + (-12.814828278938885 - 325.92980012408367*eta + 1320.102640190539*eta2)*S4;
    double uneqSpin = -148.17317525117338*dchi*delta*eta2;


    return (noSpin + eqSpin + uneqSpin);

}


/* ******************** NUMERICAL ROUTINES ******************** */
// This function determines whether x and y are approximately equal to a relative accuracy epsilon.
// Note that x and y are compared to relative accuracy, so this function is not suitable for testing whether a value is approximately zero.
bool IMRPhenomX_ApproxEqual(REAL8 x, REAL8 y, REAL8 epsilon) {
  return !gsl_fcmp(x, y, epsilon);
}

// If x and X are approximately equal to relative accuracy epsilon then set x = X.
// If X = 0 then use an absolute comparison.

void IMRPhenomX_InternalNudge(REAL8 x, REAL8 X, REAL8 epsilon) {
  if (X != 0.0) {
    if (IMRPhenomX_ApproxEqual(x, X, epsilon)) {
      XLAL_PRINT_INFO("Nudging value %.15g to %.15g\n", x, X);
      x = X;
    }
  }
  else {
    if (fabs(x - X) < epsilon)
      x = X;
  }
}

REAL8 XLALSimIMRPhenomXatan2tol(REAL8 a, REAL8 b, REAL8 tol)
{
  REAL8 c;
  if (fabs(a) < tol && fabs(b) < tol)
    c = 0.;
  else
    c = atan2(a, b);
  return c;
}

/**** Define some useful powers ****/
size_t NextPow2(const size_t n)
{
  // use pow here, not bit-wise shift, as the latter seems to run against an upper cutoff long before SIZE_MAX, at least on some platforms
  return (size_t) pow(2,ceil(log2(n)));
}

bool IMRPhenomX_StepFuncBool(const double t, const double t1) {
	return (t >= t1);
}

REAL8 XLALSimIMRPhenomXsign(REAL8 x)
{
  return (x > 0.) ? 1.0 : ((x < 0.0) ? -1.0 : 0.0);
}


/* Useful powers to avoid pow(.,.) function */
/*
 * calc square of number without floating point 'pow'
 */
double pow_2_of(double number)
{
 return (number*number);
}

/**
 * calc cube of number without floating point 'pow'
 */
double pow_3_of(double number)
{
 return (number*number*number);
}

/**
 * calc fourth power of number without floating point 'pow'
 */
double pow_4_of(double number)
{
 double pow2 = pow_2_of(number);
 return pow2 * pow2;
}

/**
 * calc fifth power of number without floating point 'pow'
 */
double pow_5_of(double number)
{
 double pow2 = pow_2_of(number);
 return pow2 * pow2 * number;
}

/**
 * calc sixth power of number without floating point 'pow'
 */
double pow_6_of(double number)
{
 double pow2 = pow_2_of(number);
 return pow2 * pow2 * pow2;
}

/**
 * calc seventh power of number without floating point 'pow'
 */
double pow_7_of(double number)
{
 double pow2 = pow_2_of(number);
 return pow2 * pow2 * pow2 * number;
}

/**
 * calc eigth power of number without floating point 'pow'
 */
double pow_8_of(double number)
{
 double pow2 = pow_2_of(number);
 double pow4 = pow2*pow2;
 return pow4 * pow4;
}

/**
 * calc ninth power of number without floating point 'pow'
 */
double pow_9_of(double number)
{
 double pow2 = pow_2_of(number);
 double pow4 = pow2*pow2;
 return pow4 * pow4 * number;
}

/*
 * Check if m1 > m2. If not, swap the masses and spin vectors such that body is the heavier compact object.
 */
INT4 XLALIMRPhenomXPCheckMassesAndSpins(
    REAL8 *m1,    /**< [out] mass of body 1 */
    REAL8 *m2,    /**< [out] mass of body 2 */
    REAL8 *chi1x, /**< [out] x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 *chi1y, /**< [out] y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 *chi1z, /**< [out] z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    REAL8 *chi2x, /**< [out] x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8 *chi2y, /**< [out] y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    REAL8 *chi2z  /**< [out] z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
)
{
    REAL8 m1_tmp, m2_tmp;
    REAL8 chi1x_tmp, chi1y_tmp, chi1z_tmp;
    REAL8 chi2x_tmp, chi2y_tmp, chi2z_tmp;

    if (*m1 > *m2)
    {
      chi1x_tmp = *chi1x;
      chi1y_tmp = *chi1y;
      chi1z_tmp = *chi1z;

      chi2x_tmp = *chi2x;
      chi2y_tmp = *chi2y;
      chi2z_tmp = *chi2z;

      m1_tmp = *m1;
      m2_tmp = *m2;
    }
    else
    {
      //printf("Swapping spins!\n");

      /* swap spins and masses */
      chi1x_tmp = *chi2x;
      chi1y_tmp = *chi2y;
      chi1z_tmp = *chi2z;

      chi2x_tmp = *chi1x;
      chi2y_tmp = *chi1y;
      chi2z_tmp = *chi1z;

      m1_tmp = *m2;
      m2_tmp = *m1;
    }

    /* Update the masses and spins */
    *m1    = m1_tmp;
    *m2    = m2_tmp;
    *chi1x = chi1x_tmp;
    *chi1y = chi1y_tmp;
    *chi1z = chi1z_tmp;

    *chi2x = chi2x_tmp;
    *chi2y = chi2y_tmp;
    *chi2z = chi2z_tmp;

    if(*m1 < *m2)
    {
      XLAL_ERROR(XLAL_EDOM,"An error occured in XLALIMRPhenomXPCheckMassesAndSpins when trying to enfore that m1 should be the larger mass.\
      After trying to enforce this m1 = %f and m2 = %f\n", *m1, *m2);

    }

    return XLAL_SUCCESS;
}


/* ******************** ANALYTICAL MODEL WRAPPERS ******************** */
/*
    "Analytical" phenomenological ringdown ansatz for phase. This is used by the higher mode functions and
    can be used to prototype or test model. Convenient wrapper exposed via XLAL.
*/
REAL8 XLALSimIMRPhenomXRingdownPhase22AnsatzAnalytical(REAL8 ff, REAL8 fRD, REAL8 fDA, REAL8 a0, REAL8 a1, REAL8 a2, REAL8 a4, REAL8 aL)
{

  REAL8 invf  = 1.0 / ff;
  REAL8 invf2 = invf  * invf;
  REAL8 invf3 = invf2 * invf;
  //REAL8 invf4 = invf2 * invf2;
  REAL8 logfv  = log(ff);

  REAL8 phaseOut;

  phaseOut = ( a0*ff + a1*logfv - a2*invf - (a4 * invf3 / 3.0) + (aL * atan( (ff - fRD)/fDA ) / fDA ) );

  return phaseOut;
}

/*
    "Analytical" phenomenological ringdown ansatz for phase derivative. This is used by the higher mode functions and
    can be used to prototype or test model. Convenient wrapper exposed via XLAL.

    a_0 + a_1 f^(-1) + a_2 f^(-2) + a_3 f^(-3) + a_4 f^(-4) + ( aRD ) / ( (f_damp^2 + (f - f_ring)^2 ) )

    where a_L = - dphase0 * aRD

    Our canonical ringdown ansatz sets a_3 = 0.
*/
REAL8 XLALSimIMRPhenomXRingdownPhaseDeriv22AnsatzAnalytical(REAL8 ff, REAL8 fRD, REAL8 fDA, REAL8 a0, REAL8 a1, REAL8 a2, REAL8 a4, REAL8 aL)
{

  REAL8 invf  = 1.0 / ff;
  REAL8 invf2 = invf  * invf;
  //REAL8 invf3 = invf2 * invf;
  REAL8 invf4 = invf2 * invf2;

  REAL8 phaseOut;

  phaseOut =  ( a0 + a1*invf + a2*invf2 + a4*invf4 + ( aL / (fDA*fDA + (ff - fRD)*(ff - fRD)) ) );

  return phaseOut;
}

/*
    "Analytical" phenomenological ringdown ansatz for amplitude. This is used by the higher mode functions but
    can also be used to prototype or test model. Convenient wrapper exposed via XLAL.
*/
REAL8 XLALSimIMRPhenomXRingdownAmplitude22AnsatzAnalytical(REAL8 ff, REAL8 fRD, REAL8 fDA, REAL8 gamma1, REAL8 gamma2, REAL8 gamma3)
{
  REAL8 gammaD13 = fDA * gamma1 * gamma3;
  REAL8 gammaR   = gamma2 / (gamma3 * fDA);
  REAL8 gammaD2  = (fDA * gamma3) * (fDA * gamma3);
  REAL8 dfr      = ff - fRD;
  REAL8 ampOut;

  ampOut =  exp(- dfr * gammaR ) * (gammaD13) / (dfr*dfr + gammaD2);

  return ampOut;
}

/*
		This is the canonical intermediate ansatz:

		a_0 + a_1 ft^(-1) + a_2 ft^(-2) + a_3 ft^(-3) + a4 ft^(-4) + (4 * a_RD) / ( (2 * f_fdamp)^2 + (f - f_ring)^2 )

		ft = (f / f_ring)
*/
REAL8 XLALSimIMRPhenomXIntermediateAmplitude22AnsatzAnalytical(REAL8 ff, REAL8 ff7o6, REAL8 a0, REAL8 a1, REAL8 a2, REAL8 a3, REAL8 a4, REAL8 a5)
{
  return ff7o6 / ( a0 + ff*(a1 + ff*(a2 + ff*(a3 + ff*(a4 + ff*a5)))) );
}

/*
		This is the canonical intermediate ansatz:

		a_0 + a_1 ft^(-1) + a_2 ft^(-2) + a_3 ft^(-3) + a4 ft^(-4) + (4 * a_RD) / ( (2 * f_fdamp)^2 + (f - f_ring)^2 )

		ft = (f / f_ring)
*/
REAL8 XLALSimIMRPhenomXIntermediatePhase22AnsatzAnalytical(REAL8 ff, REAL8 fRD, REAL8 fDA, REAL8 a0, REAL8 a1, REAL8 a2, REAL8 a3, REAL8 a4, REAL8 aL)
{
  REAL8 invff1 = 1.0 / ff;
  REAL8 invff2 = invff1 * invff1;
  REAL8 invff3 = invff2 * invff1;
  REAL8 invff4 = invff3 * invff1;

  REAL8 LorentzianTerm;
  REAL8 phaseOut;

  /* This is the Lorentzian term where aL = - a_{RD} dphase0 */
  LorentzianTerm = (4.0 * aL) / ( (4.0*fDA*fDA) + (ff - fRD)*(ff - fRD) );

  /* Return a polynomial embedded in the background from the merger */
  phaseOut       = a0 + a1*invff1 + a2*invff2 + a3*invff3 + a4*invff4 + LorentzianTerm;

  return phaseOut;
}

/* Leading order amplitude pre-factor, see Eq. 6.2 of arXiv:2001.11412 */
REAL8 XLALSimIMRPhenomXAmp22Prefactor(REAL8 eta)
{
	REAL8 ampOut;
	ampOut = sqrt(2.0/3.0) * sqrt(eta)  / pow(LAL_PI,1.0/6.0);
	return ampOut;
}



/*
  The functions below are XLAL exposed of the QNM ringdown and damping frequency used for the
  IMRPhenomX model: https://arxiv.org/abs/2001.11412.
*/

/*
   Ringdown frequency for 22 mode, given final dimensionless spin

   https://arxiv.org/src/2001.10914v1/anc/QNMs/CoefficientStatsfring22.m
*/
REAL8 XLALSimIMRPhenomXfring22(
  const REAL8 af /* Dimensionless final spin */
)
{
  REAL8 return_val;

  if (fabs(af) > 1.0) {
         XLAL_ERROR(XLAL_EDOM, "PhenomX evaluate_QNMfit_fring22 \
  function: |finalDimlessSpin| > 1.0 not supported");
  }

  REAL8 x2= af*af;
  REAL8 x3= x2*af;
  REAL8 x4= x2*x2;
  REAL8 x5= x3*x2;
  REAL8 x6= x3*x3;
  REAL8 x7= x4*x3;

  return_val = (0.05947169566573468 - \
  0.14989771215394762*af + 0.09535606290986028*x2 + \
  0.02260924869042963*x3 - 0.02501704155363241*x4 - \
  0.005852438240997211*x5 + 0.0027489038393367993*x6 + \
  0.0005821983163192694*x7)/(1 - 2.8570126619966296*af + \
  2.373335413978394*x2 - 0.6036964688511505*x4 + \
  0.0873798215084077*x6);
  return return_val;
}


/*
   Damping frequency for 22 mode, given final dimensionless spin.

   https://arxiv.org/src/2001.10914v1/anc/QNMs/CoefficientStatsfdamp22.m
*/
REAL8 XLALSimIMRPhenomXfdamp22(
  const REAL8 af /* Dimensionless final spin */
)
{
  REAL8 return_val;

  if (fabs(af) > 1.0) {
         XLAL_ERROR(XLAL_EDOM, "PhenomX evaluate_QNMfit_fdamp22 \
  function: |finalDimlessSpin| > 1.0 not supported");
  }

  REAL8 x2= af*af;
  REAL8 x3= x2*af;
  REAL8 x4= x2*x2;
  REAL8 x5= x3*x2;
  REAL8 x6= x3*x3;

  return_val = (0.014158792290965177 - \
  0.036989395871554566*af + 0.026822526296575368*x2 + \
  0.0008490933750566702*x3 - 0.004843996907020524*x4 - \
  0.00014745235759327472*x5 + 0.0001504546201236794*x6)/(1 - \
  2.5900842798681376*af + 1.8952576220623967*x2 - \
  0.31416610693042507*x4 + 0.009002719412204133*x6);
  return return_val;
}
