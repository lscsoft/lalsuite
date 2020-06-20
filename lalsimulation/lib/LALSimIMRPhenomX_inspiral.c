/*
 * Copyright (C) 2018 Geraint Pratten
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
 */

#include <lal/XLALError.h>

#include "LALSimIMRPhenomXUtilities.h"
#include "LALSimIMRPhenomX_internals.h"

/******************************* IMRPhenomX Amplitude Functions: Inspiral *******************************/
/*
	Phenomenological coefficients for pseudo-PN corrections to TaylorF2 insiral:
	Amp(f,eta,chi1,chi2) = TF2(f,eta,chi1,chi2) + A0 * ( delta1*f**(7/3) + delta2*f**(8/3) + delta3*f**(9/3) )
*/


/*
    Value of amplitude collocation point at 2/4 f^A_T, defined in Eq. 5.7 of arXiv:2001.11412
   
    Effective spin parameterization used = chiPNHat
*/
static double IMRPhenomX_Inspiral_Amp_22_v2(double eta, double S, double dchi, double delta, int InsAmpFlag){

  double eta2 = (eta*eta);
  double eta3 = (eta2*eta);
  double eta4 = (eta3*eta);

  double S2 = (S*S);
  double S3 = (S2*S);

  double noSpin, eqSpin, uneqSpin;

  switch ( InsAmpFlag )
	{
		case 103:
		{
      noSpin = (-0.015178276424448592 - 0.06098548699809163*eta + 0.4845148547154606*eta2)/(1. + 0.09799277215675059*eta);

      eqSpin = ((0.02300153747158323 + 0.10495263104245876*eta2)*S + (0.04834642258922544 - 0.14189350657140673*eta)*eta*S3
                    + (0.01761591799745109 - 0.14404522791467844*eta2)*S2)/(1. - 0.7340448493183307*S);

      uneqSpin = dchi*delta*eta4*(0.0018724905795891192 + 34.90874132485147*eta);

  		break;
		}
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Amp_22_v2: IMRPhenomXInspiralAmpVersion is not valid.\n");
    }
	}

  return (noSpin + eqSpin + uneqSpin);

}

/*
    Value of amplitude collocation point at 3/4 f^A_T, defined in Eq. 5.7 of arXiv:2001.11412
   
    Effective spin parameterization used = chiPNHat
*/
static double IMRPhenomX_Inspiral_Amp_22_v3(double eta, double S, double dchi, double delta, int InsAmpFlag){

  double eta2 = (eta*eta);
  double eta3 = (eta2*eta);
  double eta4 = (eta3*eta);

  double noSpin, eqSpin, uneqSpin;

  switch ( InsAmpFlag )
	{
		case 103:
		{
      noSpin = (-0.058572000924124644 - 1.1970535595488723*eta + 8.4630293045015*eta2)/(1. + 15.430818840453686*eta);

      eqSpin = ((-0.08746408292050666 + eta*(-0.20646621646484237 - 0.21291764491897636*S)
                  + eta2*(0.788717372588848 + 0.8282888482429105*S) - 0.018924013869130434*S)*S)/(-1.332123330797879 + 1.*S);

      uneqSpin = dchi*delta*eta4*(0.004389995099201855 + 105.84553997647659*eta);

  		break;
		}
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Amp_22_v3: IMRPhenomXInspiralAmpVersion is not valid.\n");
    }
	}

  return (noSpin + eqSpin + uneqSpin);

}

/*
    Value of amplitude collocation point at 4/4 f^A_T, defined in Eq. 5.7 of arXiv:2001.11412
   
    Effective spin parameterization used = chiPNHat
*/
static double IMRPhenomX_Inspiral_Amp_22_v4(double eta, double S, double dchi, double delta, int InsAmpFlag){

  double eta2 = (eta*eta);
  double eta3 = (eta2*eta);
  double eta4 = (eta3*eta);

  double S2 = (S*S);

  double noSpin, eqSpin, uneqSpin;

  switch ( InsAmpFlag )
	{
		case 103:
		{
      noSpin = (-0.16212854591357853 + 1.617404703616985*eta - 3.186012733446088*eta2 + 5.629598195000046*eta3)/(1. + 0.04507019231274476*eta);

      eqSpin = (S*(1.0055835408962206 + eta2*(18.353433894421833 - 18.80590889704093*S) - 0.31443470118113853*S
                  + eta*(-4.127597118865669 + 5.215501942120774*S) + eta3*(-41.0378120175805
                    + 19.099315016873643*S)))/(5.852706459485663 - 5.717874483424523*S + 1.*S2);

      uneqSpin = dchi*delta*eta4*(0.05575955418803233 + 208.92352600701068*eta);

  		break;
		}
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Amp_22_v4: IMRPhenomXInspiralAmpVersion is not valid.\n");
    }
	}

  return (noSpin + eqSpin + uneqSpin);

}

/*
    f1 = freq. at 2/4 f^A_T
    f2 = freq. at 3/4 f^A_T
    f3 = freq. at 4/4 f^A_T

    v1 = Value of collocation point at f1
    v2 = Value of collocation point at f2
    v3 = Value of collocation point at f3

    These expressions do *not* depend on calibration if model assumes:
    TF2 + rho1 f^(7/3) + rho2 f^(8/3) + rho3 f^(9/3)
*/

/* Get Pseudo PN coefficient at f^(7/3) */
static double IMRPhenomX_Inspiral_Amp_22_rho1(double v1, double v2, double v3, double F1, double F2, double F3, int InsAmpFlag)
{
  double retVal;

  switch ( InsAmpFlag )
  {
    case 103:
    {
      double f1p1o3, f2p1o3, f3p1o3;
      double f1p7o3, f1p8o3, f2p7o3, f2p8o3, f3p7o3, f3p8o3;

      f1p1o3 = cbrt(F1); // f1^(1/3)
      f2p1o3 = cbrt(F2); // f2^(1/3)
      f3p1o3 = cbrt(F3); // f3^(1/3)

      f1p7o3 = pow_7_of(f1p1o3); // f1^(7/3)
      f1p8o3 = f1p7o3 * f1p1o3;  // f1^(8/3)

      f2p7o3 = pow_7_of(f2p1o3); // f2^(7/3)
      f2p8o3 = f2p7o3 * f2p1o3;  // f2^(8/3)

      f3p7o3 = pow_7_of(f3p1o3); // f3^(7/3)
      f3p8o3 = f3p7o3 * f3p1o3;  // f3^(8/3)

      retVal = (-(f2p8o3*(F3*F3*F3)*v1) + F2*F2*F2*f3p8o3*v1 + f1p8o3*(F3*F3*F3)*v2 - F1*F1*F1*f3p8o3*v2 - f1p8o3*(F2*F2*F2)*v3 + F1*F1*F1*f2p8o3*v3)
                      / (f1p7o3*(f1p1o3 - f2p1o3)*f2p7o3*(f1p1o3 - f3p1o3)*(f2p1o3 - f3p1o3)*f3p7o3);

      break;
    }
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Amp_22_rho1: IMRPhenomXInspiralAmpVersion is not valid.\n");
    }
  }

  return retVal;
}

/* Get Pseudo PN coefficient at f^(8/3) */
static double IMRPhenomX_Inspiral_Amp_22_rho2(double v1, double v2, double v3, double F1, double F2, double F3, int InsAmpFlag)
{
  double retVal;

  switch ( InsAmpFlag )
  {
    case 103:
    {
      double f1p1o3, f2p1o3, f3p1o3;
      double f1p7o3, f2p7o3, f3p7o3;

      f1p1o3 = cbrt(F1); // f1^(1/3)
      f2p1o3 = cbrt(F2); // f2^(1/3)
      f3p1o3 = cbrt(F3); // f3^(1/3)

      f1p7o3 = pow_7_of(f1p1o3); // f1^(7/3)

      f2p7o3 = pow_7_of(f2p1o3); // f2^(7/3)

      f3p7o3 = pow_7_of(f3p1o3); // f3^(7/3)

      retVal = ( f2p7o3*(F3*F3*F3)*v1 - F2*F2*F2*f3p7o3*v1 - f1p7o3*(F3*F3*F3)*v2 + F1*F1*F1*f3p7o3*v2 + f1p7o3*(F2*F2*F2)*v3 - F1*F1*F1*f2p7o3*v3 )
                        / (f1p7o3*(f1p1o3 - f2p1o3)*f2p7o3*(f1p1o3 - f3p1o3)*(f2p1o3 - f3p1o3)*f3p7o3);

      break;
    }
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Amp_22_rho2: IMRPhenomXInspiralAmpVersion is not valid.\n");
    }
  }

  return retVal;
}

/* Get Pseudo PN coefficient at f^(9/3) */
static double IMRPhenomX_Inspiral_Amp_22_rho3(double v1, double v2, double v3, double F1, double F2, double F3, int InsAmpFlag)
{
  double retVal;

  switch ( InsAmpFlag )
  {
    case 103:
    {
      double f1p1o3, f2p1o3, f3p1o3;
      double f1p7o3, f1p8o3, f2p7o3, f2p8o3, f3p7o3, f3p8o3;

      f1p1o3 = cbrt(F1); // f1^(1/3)
      f2p1o3 = cbrt(F2); // f2^(1/3)
      f3p1o3 = cbrt(F3); // f3^(1/3)

      f1p7o3 = pow_7_of(f1p1o3); // f1^(7/3)
      f1p8o3 = f1p7o3 * f1p1o3;  // f1^(8/3)

      f2p7o3 = pow_7_of(f2p1o3); // f2^(7/3)
      f2p8o3 = f2p7o3 * f2p1o3;  // f2^(8/3)

      f3p7o3 = pow_7_of(f3p1o3); // f3^(7/3)
      f3p8o3 = f3p7o3 * f3p1o3;  // f3^(8/3)

      retVal = ( f2p8o3*f3p7o3*v1 - f2p7o3*f3p8o3*v1 - f1p8o3*f3p7o3*v2 + f1p7o3*f3p8o3*v2 + f1p8o3*f2p7o3*v3 - f1p7o3*f2p8o3*v3 )
                        / ( f1p7o3*(f1p1o3 - f2p1o3)*f2p7o3*(f1p1o3 - f3p1o3)*(f2p1o3 - f3p1o3)*f3p7o3 );

      break;
    }
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Amp_22_rho3: IMRPhenomXInspiralAmpVersion is not valid.\n");
    }
  }

  return retVal;
}


/*
 *  TaylorF2 PN Amplitude + pseudo-PN coefficients
 */
static double IMRPhenomX_Inspiral_Amp_22_Ansatz(double Mf, IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXWaveformStruct *pWF, IMRPhenomXAmpCoefficients *pAmp)
{
  double pnAmp;

  int InsAmpFlag = pWF->IMRPhenomXInspiralAmpVersion;

  switch ( InsAmpFlag )
  {
    case 103:
    {
      // Re-factor expression
      pnAmp = (
          pAmp->pnInitial // 1.0
        + powers_of_Mf->one_third      * pAmp->pnOneThird
        + powers_of_Mf->two_thirds     * pAmp->pnTwoThirds
        + Mf                           * pAmp->pnThreeThirds
        + Mf*(
          + powers_of_Mf->one_third    * pAmp->pnFourThirds
          + powers_of_Mf->two_thirds   * pAmp->pnFiveThirds
          + Mf                         * pAmp->pnSixThirds
          + Mf*(
            + powers_of_Mf->one_third  * pAmp->rho1
            + powers_of_Mf->two_thirds * pAmp->rho2
            + Mf                       * pAmp->rho3
            )
          )
      );
      break;
    }
    default :
    {
        pnAmp = 0.0;
    }
  }

  return pnAmp;
}

/*
 *  Derivative of TaylorF2 PN Amplitude + pseudo-PN coefficients
 */
static double IMRPhenomX_Inspiral_Amp_22_DAnsatz(double Mf, IMRPhenomXWaveformStruct *pWF, IMRPhenomXAmpCoefficients *pAmp) {

  double DAmpIns;

  int InsAmpFlag = pWF->IMRPhenomXInspiralAmpVersion;

  switch ( InsAmpFlag )
  {
    case 103:
    {
      double eta   = pWF->eta;
      double chi1L = pWF->chi1L;
      double chi2L = pWF->chi2L;
      double rho1  = pAmp->rho1;
      double rho2  = pAmp->rho2;
      double rho3  = pAmp->rho3;

      double chi1L2 = chi1L  * chi1L;
      double chi1L3 = chi1L2 * chi1L;
      double chi2L2 = chi2L  * chi2L;

      double eta2   = eta*eta;
      double eta3   = eta*eta2;
      double Mf2    = Mf*Mf;
      double LALPi  = LAL_PI;

      double delta  = pWF->delta;

      DAmpIns =
      ( ((chi2L*(81 - 81*delta - 44*eta) + chi1L*(81*(1 + delta) - 44*eta))*LALPi)/48.
      + ((-969 + 1804*eta)*pow(LALPi,2./3.))/(1008.*pow(Mf,1/3.))
      + ((-27312085 - 10287648*chi2L2 + 10287648*chi2L2*delta - 10287648*chi1L2*(1 + delta)
         + 24*(-1975055 + 857304*chi1L2 - 994896*chi1L*chi2L + 857304*chi2L2)*eta + 35371056*eta2)*pow(LALPi,4./3.)*pow(Mf,1./3.))/6.096384e6
      + (5*pow(LALPi,5./3.)*(-6048*chi1L3*(-1 - delta + (3 + delta)*eta) + chi1L*(287213*(1 + delta) - 4*(93414 + 2083*delta)*eta - 35632*eta2)
               + chi2L*(-((287213 + 6048*chi2L2)*(-1 + delta)) + 4*(-93414 + 1512*chi2L2*(-3 + delta) + 2083*delta)*eta - 35632*eta2)
               + 42840*(-1 + 4*eta)*LALPi)*pow(Mf,2./3.))/96768.
       - (pow(LALPi,2.0)*(-336*(-3248849057 + 1809550512*chi1L2 - 2954929824*chi1L*chi2L + 1809550512*chi2L2)*eta2 - 324322727232*eta3
             + 7*(177520268561 + 29362199328*chi2L2 - 29362199328*chi2L2*delta + 29362199328*chi1L2*(1 + delta)
             + 12160253952*(chi1L + chi2L + chi1L*delta - chi2L*delta)*LALPi)
             + 12*eta*(-545384828789 + 49568837472*chi1L*chi2L - 12312458928*chi2L2 - 21943440288*chi2L2*delta
             + 77616*chi1L2*(-158633 + 282718*delta) - 8345272320*(chi1L + chi2L)*LALPi
             + 21384760320*pow(LALPi,2.0)))*Mf)/3.0042980352e10
       + (7.0/3.0)*pow(Mf,4.0/3.0)*rho1 + (8.0/3.0)*pow(Mf,5.0/3.0)*rho2 + 3*Mf2*rho3 );

      break;
    }
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Amp_22_DAnsatz: IMRPhenomXInspiralAmpVersion is not valid.\n");
    }
  }

  return DAmpIns;
}

/******************************* Phase Functions: Inspiral *******************************/
/*
	Phenomenological coefficients for pseudo-PN corrections to TaylorF2 insiral:
	Psi(f,eta,chi1,chi2) = TF2(f,eta,chi1,chi2) + (1/eta) * ( alpha0 + alpha1*f**(3/3) + (3/4)*alpha2*f**(4/3) + (3/5)*alpha3*f**(5/3) + (1/2)*alpha4*f**(6/3) )
*/

/*
	This code is designed and intended to be modular. If a new PN TaylorF2 approximant is produced or a new fit against NR
	hybrids generated, then you can add the fit to the collocation points by simply adding in the extra case.

	This is intended to avoid a massive over-duplication of near-identical functions.
*/

/*
    Value of phase collocation point at v_3. See section VII.A of arXiv:2001.11412
   
    Effective spin parameterization used = chiPNHat
*/
static double IMRPhenomX_Inspiral_Phase_22_v3(double eta, double S, double dchi, double delta, int InspPhaseFlag){

	double eta2  = eta*eta;
	double eta3  = eta2*eta;
	double eta4  = eta3*eta;
	double eta5  = eta4*eta;

	double S2    = S*S;
	double S3    = S2*S;

  double dchi2 = dchi*dchi;

	double noSpin, eqSpin, uneqSpin;

	switch ( InspPhaseFlag )
	{
		case 104: 			/* Canonical, 3 pseudo PN terms */
		{
      noSpin = (15415.000000000007 + 873401.6255736464*eta + 376665.64637025696*eta2 - 3.9719980569125614e6*eta3 + 8.913612508054944e6*eta4)/(1. + 46.83697749859996*eta);

      eqSpin = (S*(397951.95299014193 - 207180.42746987*S + eta3*(4.662143741417853e6 - 584728.050612325*S - 1.6894189124921719e6*S2) + eta*(-1.0053073129700898e6 + 1.235279439281927e6*S - 174952.69161683554*S2) - 130668.37221912303*S2 + eta2*(-1.9826323844247842e6 + 208349.45742548333*S + 895372.155565861*S2)))/(-9.675704197652225 + 3.5804521763363075*S + 2.5298346636273306*S2 + 1.*S3);

      uneqSpin = -1296.9289110696955*dchi2*eta + dchi*delta*eta*(-24708.109411857182 + 24703.28267342699*eta + 47752.17032707405*S);

  		break;
		}
		case 105:				/* Canonical, 4 pseudo PN terms */
		{
      noSpin = (11717.332402222377 + 4.972361134612872e6*eta + 2.137585030930089e7*eta2 - 8.882868155876668e7*eta3 + 2.4104945956043008e8*eta4 - 2.3143426271719798e8*eta5)/(1. + 363.5524719849582*eta);

      eqSpin = (S*(52.42001436116159 - 50.547943589389966*S + eta3*S*(-15355.56020802297 + 20159.588079899433*S) + eta2*(-286.9576245212502 + 2795.982637986682*S - 2633.1870842242447*S2) - 1.0824224105690476*S2 + eta*(-123.78531181532225 + 136.1961976556154*S - 7.534890781638552*S3) + 5.973206330664007*S3 + eta4*(1777.2176433016125 + 24069.288079063674*S - 44794.9522164669*S2 + 1584.1515277998406*S3)))/(-0.0015307616935628491 + (0.0010676159178395538 - 0.25*eta3 + 1.*eta4)*S);

      uneqSpin = -1357.9794908614106*dchi2*eta + dchi*delta*eta*(-23093.829989687543 + 21908.057881789653*eta + 49493.91485992256*S);

  		break;
		}
		case 114:				/* Extended, 3 pseudo PN terms */
		{
      noSpin = (68014.00000000003 + 1.1513072539654972e6*eta - 2.725589921577228e6*eta2 + 312571.92531733884*eta3)/(1. + 17.48539665509149*eta);

      eqSpin = (S*(-34467.00643820664 + 99693.81839115614*eta + 144345.24343461913*eta4 + (23618.044919850676 - 89494.69555164348*eta + 725554.5749749158*eta4 - 103449.15865381068*eta2)*S + (10350.863429774612 - 73238.45609787296*eta + 3.559251543095961e6*eta4 + 888228.5439003729*eta2 - 3.4602940487291473e6*eta3)*S2))/(1. - 0.056846656084188936*S - 0.32681474740130184*S2 - 0.30562055811022015*S3);

      uneqSpin = -1182.4036752941936*dchi2*eta + dchi*delta*eta*(-0.39185419821851025 - 99764.21095663306*eta + 41826.177356107364*S);

  		break;
		}
		case 115:				/* Extended, 4 pseudo PN terms */
		{
      noSpin = (60484.00000000003 + 4.370611564781374e6*eta - 5.785128542827255e6*eta2 - 8.82141241633613e6*eta3 + 1.3406996319926713e7*eta4)/(1. + 70.89393713617065*eta);

      eqSpin = (S*(21.91241092620993 - 32.57779678272056*S + eta2*(-102.4982890239095 + 2570.7566494633033*S - 2915.1250015652076*S2) + 8.130585173873232*S2 + eta*(-28.158436727309365 + 47.42665852824481*S2) + eta3*(-1635.6270690726785 - 13745.336370568011*S + 19539.310912464192*S2) + 1.2792283911312285*S3 + eta4*(5558.382039622131 + 21398.7730201213*S - 37773.40511355719*S2 + 768.6183344184254*S3)))/(-0.0007758753818017038 + (0.0005304023864415552 - 0.25000000000000006*eta3 + 1.*eta4)*S);

      uneqSpin = -1223.769262298142*dchi2*eta + dchi*delta*eta*(-16.705471562129436 - 93771.93750060834*eta + 43675.70151058481*S);

  		break;
		}
    default:
    {
        XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Phase_22_v3: NPseudoPN requested is not valid.\n");
    }
	}

	return (noSpin + eqSpin + uneqSpin);
}

/*
    Value of phase collocation point for d13 = v1 - v3. See section VII.A of arXiv:2001.11412
   
    Effective spin parameterization used = chiPNHat
*/
static double IMRPhenomX_Inspiral_Phase_22_d13(double eta, double S, double dchi, double delta, int InspPhaseFlag){

	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double eta4 = eta3*eta;
	double eta5 = eta4*eta;

	double S2   = S*S;
	double S3   = S2*S;
	double S4   = S3*S;

	double noSpin, eqSpin, uneqSpin;

	switch ( InspPhaseFlag )
	{
		case 104: 			/* Canonical, 3 pseudo PN terms */
		{
      noSpin = (-17294.000000000007 - 19943.076428555978*eta + 483033.0998073767*eta2)/(1. + 4.460294035404433*eta);

      eqSpin = (S*(68384.62786426462 + 67663.42759836042*S - 2179.3505885609297*S2 + eta*(-58475.33302037833 + 62190.404951852535*S + 18298.307770807573*S2 - 303141.1945565486*S3) + 19703.894135534803*S3 + eta2*(-148368.4954044637 - 758386.5685734496*S - 137991.37032619823*S2 + 1.0765877367729193e6*S3) + 32614.091002011017*S4))/(2.0412979553629143 + 1.*S);

      uneqSpin = 12017.062595934838*dchi*delta*eta;
  		break;
		}
		case 105:				/* Canonical, 4 pseudo PN terms */
		{
      noSpin = (-14234.000000000007 + 16956.107542097994*eta + 176345.7518697656*eta2)/(1. + 1.294432443903631*eta);

      eqSpin = (S*(814.4249470391651 + 539.3944162216742*S + 1985.8695471257474*S2 + eta*(-918.541450687484 + 2531.457116826593*S - 14325.55532310136*S2 - 19213.48002675173*S3) + 1517.4561843692861*S3 + eta2*(-517.7142591907573 - 14328.567448748548*S + 21305.033147575057*S2 + 50143.99945676916*S3)))/(0.03332712934306297 + 0.0025905919215826172*S + (0.07388087063636305 - 0.6397891808905215*eta + 1.*eta2)*S2);

      uneqSpin = dchi*delta*eta*(0.09704682517844336 + 69335.84692284222*eta);

  		break;
		}
		case 114:				/* Extended, 3 pseudo PN terms */
		{
      noSpin = (-36664.000000000015 + 277640.10051158903*eta - 581120.4916255298*eta2 + 1.415628418251648e6*eta3 - 7.640937162029471e6*eta4 + 1.1572710625359124e7*eta5)/(1. - 4.011038704323779*eta);

      eqSpin = (S*(-38790.01253014577 - 50295.77273512981*S + 15182.324439704937*S2 + eta2*(57814.07222969789 + 344650.11918139807*S + 17020.46497164955*S2 - 574047.1384792664*S3) + 24626.598127509922*S3 + eta*(23058.264859112394 - 16563.935447608965*S - 36698.430436426395*S2 + 105713.91549712936*S3)))/(-1.5445637219268247 - 0.24997068896075847*S + 1.*S2);

      uneqSpin = 74115.77361380383*dchi*delta*eta2;

  		break;
		}
		case 115:				/* Extended, 4 pseudo PN terms */
		{
      noSpin = (-29240.00000000001 - 12488.41035199958*eta + 1.3911845288427814e6*eta2 - 3.492477584609041e6*eta3)/(1. + 2.6711462529779824*eta - 26.80876660227278*eta2);

      eqSpin = (S*(-29536.155624432842 - 40024.5988680615*S + 11596.401177843705*S2 + eta2*(122185.06346551726 + 351091.59147835104*S - 37366.6143666202*S2 - 505834.54206320125*S3) + 20386.815769841945*S3 + eta*(-9638.828456576934 - 30683.453790630676*S - 15694.962417099561*S2 + 91690.51338194775*S3)))/(-1.5343852108869265 - 0.2215651087041266*S + 1.*S2);

      uneqSpin = 68951.75344813892*dchi*delta*eta2;

  		break;
		}
    default:
    {
        XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Phase_22_d13: NPseudoPN requested is not valid.\n");
    }
	}

	return (noSpin + eqSpin + uneqSpin);
}

/*
    Value of phase collocation point for d23 = v2 - v3. See section VII.A of arXiv:2001.11412
   
    Effective spin parameterization used = chiPNHat
*/
static double IMRPhenomX_Inspiral_Phase_22_d23(double eta, double S, double dchi, double delta, int InspPhaseFlag){

	double eta2  = eta*eta;
	double eta3  = eta2*eta;
	
	double dchi2 = dchi * dchi;

	double S2    = S*S;
	double S3    = S2*S;
	double S4    = S3*S;

	double noSpin, eqSpin, uneqSpin;

	switch ( InspPhaseFlag )
	{
		case 104: 			/* Canonical, 4 pseudo PN coefficients */
		{
      noSpin = (-7579.300000000004 - 120297.86185566607*eta + 1.1694356931282217e6*eta2 - 557253.0066989232*eta3)/(1. + 18.53018618227582*eta);

      eqSpin = (S*(-27089.36915061857 - 66228.9369155027*S + eta2*(150022.21343386435 - 50166.382087278434*S - 399712.22891153296*S2) - 44331.41741405198*S2 + eta*(50644.13475990821 + 157036.45676788126*S + 126736.43159783827*S2) + eta3*(-593633.5370110178 - 325423.99477314285*S + 847483.2999508682*S2)))/(-1.5232497464826662 - 3.062957826830017*S - 1.130185486082531*S2 + 1.*S3);

      uneqSpin = 3843.083992827935*dchi*delta*eta;

  		break;
		}
		case 105:				/* Canonical, 5 pseudo PN coefficients */
		{
      noSpin = (-7520.900000000003 - 49463.18828584058*eta + 437634.8057596484*eta2)/(1. + 9.10538019868398*eta);

      eqSpin = (S*(25380.485895523005 + 30617.553968012628*S + 5296.659585425608*S2 + eta*(-49447.74841021405 - 94312.78229903466*S - 5731.614612941746*S3) + 2609.50444822972*S3 + 5206.717656940992*S4 + eta2*(54627.758819129864 + 157460.98527210607*S - 69726.85196686552*S2 + 4674.992397927943*S3 + 20704.368650323784*S4)))/(1.5668927528319367 + 1.*S);

      uneqSpin = -95.38600275845481*dchi2*eta + dchi*delta*eta*(3271.8128884730654 + 12399.826237672185*eta + 9343.380589951552*S);

  		break;
		}
		case 114:				/* Extended, 4 pseudo PN coefficients */
		{
      noSpin = (-17762.000000000007 - 1.6929191194109183e6*eta + 8.420903644926643e6*eta2)/(1. + 98.061533474615*eta);

      eqSpin = (S*(-46901.6486082098 - 83648.57463631754*S + eta2*(1.2502334322912344e6 + 1.4500798116821344e6*S - 1.4822181506831646e6*S2) - 41234.966418619966*S2 + eta*(-24017.33452114588 - 15241.079745314566*S + 136554.48806839858*S2) + eta3*(-3.584298922116994e6 - 3.9566921791790277e6*S + 4.357599992831832e6*S2)))/(-3.190220646817508 - 3.4308485421201387*S - 0.6347932583034377*S2 + 1.*S3);

      uneqSpin = 24906.33337911219*dchi*delta*eta2;

  		break;
		}
		case 115:				/* Extended, 5 pseudo PN coefficients */
		{
      noSpin = (-18482.000000000007 - 1.2846128476247871e6*eta + 4.894853535651343e6*eta2 + 3.1555931338015324e6*eta3)/(1. + 82.79386070797756*eta);

      eqSpin = (S*(-19977.10130179636 - 24729.114403562427*S + 10325.295899053815*S2 + eta*(30934.123894659646 + 58636.235226102894*S - 32465.70372990005*S2 - 38032.16219587224*S3) + 15485.725898689267*S3 + eta2*(-38410.1127729419 - 87088.84559983511*S + 61286.73536122325*S2 + 42503.913487705235*S3)))/(-1.5148031011828222 - 0.24267195338525768*S + 1.*S2);

      uneqSpin = 5661.027157084334*dchi*delta*eta;

  		break;
		}
    default:
    {
        XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Phase_22_d23: NPseudoPN requested is not valid.\n");
    }
	}

	return (noSpin + eqSpin + uneqSpin);
}

/*
    Value of phase collocation point for d43 = v4 - v3. See section VII.A of arXiv:2001.11412
   
    Effective spin parameterization used = chiPNHat
*/
static double IMRPhenomX_Inspiral_Phase_22_d43(double eta, double S, double dchi, double delta, int InspPhaseFlag){

	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double eta4 = eta3*eta;

	double S2   = S*S;
	double S3   = S2*S;
	
	double noSpin, eqSpin, uneqSpin;

	switch ( InspPhaseFlag )
	{
		case 104: 			/* Canonical, 3 pseudo PN coefficients */
		{
      noSpin = (2439.000000000001 - 31133.52170083207*eta + 28867.73328134167*eta2)/(1. + 0.41143032589262585*eta);

      eqSpin = (S*(16116.057657391262 + eta3*(-375818.0132734753 - 386247.80765802023*S) + eta*(-82355.86732027541 - 25843.06175439942*S) + 9861.635308837876*S + eta2*(229284.04542668918 + 117410.37432997991*S)))/(-3.7385208695213668 + 0.25294420589064653*S + 1.*S2);

      uneqSpin = 194.5554531509207*dchi*delta*eta;

  		break;
		}
		case 105:				/* Canonical, 4 pseudo PN coefficients */
		{
      noSpin = (4085.300000000002 + 62935.7755506329*eta - 1.3712743918777364e6*eta2 + 5.024685134555112e6*eta3 - 3.242443755025284e6*eta4)/(1. + 20.889132970603523*eta - 99.39498823723363*eta2);

      eqSpin = (S*(-299.6987332025542 - 106.2596940493108*S + eta3*(2383.4907865977148 - 13637.11364447208*S - 14808.138346145908*S2) + eta*(1205.2993091547498 - 508.05838536573464*S - 1453.1997617403304*S2) + 132.22338129554674*S2 + eta2*(-2438.4917103042208 + 5032.341879949591*S + 7026.9206794027405*S2)))/(0.03089183275944264 + 1.*eta3*S - 0.010670764224621489*S2);

      uneqSpin = -1392.6847421907178*dchi*delta*eta;

  		break;
		}
		case 114:				/* Extended, 3 pseudo PN coefficients */
		{
      noSpin = (5749.000000000003 - 37877.95816426952*eta)/(1. + 1.1883386102990128*eta);

      eqSpin = ((-4285.982163759047 + 24558.689969419473*eta - 49270.2296311733*eta2)*S + eta*(-24205.71407420114 + 70777.38402634041*eta)*S2 + (2250.661418551257 + 187.95136178643946*eta - 11976.624134935797*eta2)*S3)/(1. - 0.7220334077284601*S);

      uneqSpin = dchi*delta*eta*(339.69292150803585 - 3459.894150148715*S);

  		break;
		}
		case 115:				/* Extended, 4 pseudo PN coefficients */
		{
      noSpin = (9760.400000000005 + 9839.852773121198*eta - 398521.0434645335*eta2 + 267575.4709475981*eta3)/(1. + 6.1355249449135005*eta);

      eqSpin = (S*(-1271.406488219572 + eta2*(-9641.611385554736 - 9620.333878140807*S) - 1144.516395171019*S + eta*(5155.337817255137 + 4450.755534615418*S)))/(0.1491519640750958 + (-0.0008208549820159909 - 0.15468508831447628*eta + 0.7266887643762937*eta2)*S + (0.02282542856845755 - 0.445924460572114*eta + 1.*eta2)*S2);

      uneqSpin = -1366.7949288045616*dchi*delta*eta;

  		break;
		}
    default:
    {
        XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Phase_22_d43: NPseudoPN requested is not valid.\n");
    }
	}

	return (noSpin + eqSpin + uneqSpin);
}

/*
    Value of phase collocation point for d53 = v5 - v3. See section VII.A of arXiv:2001.11412
   
    Effective spin parameterization used = chiPNHat
*/
static double IMRPhenomX_Inspiral_Phase_22_d53(double eta, double S, double dchi, double delta, int InspPhaseFlag){

	double eta2 = eta*eta;
	double eta3 = eta2*eta;

	double S2   = S*S;

	double noSpin, eqSpin, uneqSpin;

	switch ( InspPhaseFlag )
	{
		case 104: 			/* This should not be called for 4 pseudo-PN coefficients. Return 0 and print warning just in case... */
		{
			XLAL_ERROR_REAL8(XLAL_EINVAL, "Calling IMRPhenomX_Inspiral_Phase_22_d53 but trying to pass InspPhaseFlag for 4 pseudo-PN coefficients. Check this.\n");

  		break;
		}
		case 105:				/* Canonical, 5 pseudo PN coefficients */
		{
      noSpin = (5474.400000000003 + 131008.0112992443*eta - 1.9692364337640922e6*eta2 + 1.8732325307375633e6*eta3)/(1. + 32.90929274981482*eta);

      eqSpin = (S*(18609.016486281424 - 1337.4947536109685*S + eta2*(219014.98908698096 - 307162.33823247004*S - 124111.02067626518*S2) - 7394.9595046977365*S2 + eta*(-87820.37490863055 + 53564.4178831741*S + 34070.909093771494*S2) + eta3*(-248096.84257893753 + 536024.5354098587*S + 243877.65824670633*S2)))/(-1.5282904337787517 + 1.*S);

      uneqSpin = -429.1148607925461*dchi*delta*eta;

  		break;
		}
		case 114:				/* Extended, 4 pseudo PN coefficients */
		{
			XLAL_ERROR_REAL8(XLAL_EINVAL, "Calling IMRPhenomX_Inspiral_Phase_22_d53 but trying to pass InspPhaseFlag for 4 pseudo-PN coefficients. Check this.\n");
		}
		case 115:				/* Extended, 5 pseudo PN terms */
		{
      noSpin = (12971.000000000005 - 93606.05144508784*eta + 102472.4473167639*eta2)/(1. - 0.8909484992212859*eta);

      eqSpin = (S*(16182.268123259992 + 3513.8535400032874*S + eta2*(343388.99445324624 - 240407.0282222587*S - 312202.59917289804*S2) - 10814.056847109632*S2 + eta*(-94090.9232151429 + 35305.66247590705*S + 65450.36389642103*S2) + eta3*(-484443.15601144277 + 449511.3965208116*S + 552355.592066788*S2)))/(-1.4720837917195788 + 1.*S);

      uneqSpin = -494.2754225110706*dchi*delta*eta;

  		break;
		}
    default:
    {
        XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Inspiral_Phase_22_d53: NPseudoPN requested is not valid.\n");
    }
	}

	return (noSpin + eqSpin + uneqSpin);
}

/*
		See section VII.A of arXiv:2001.11412

		This function solves for the pseudo-PN coefficients. The structure of the equations is:

			c + alpha * x + beta * x^2 + gamma * x^3 + xi * x^4

		where x = f^(1/3).

		For 3 pseudo-PN parameters, we solve for: (c,alpha,beta,gamma).
		For 4 pseudo-PN parameters, we solve for: (c,alpha,beta,gamma,xi).

		Phase Derivative: TaylorF2 + pseudo-PN coefficients
*/
static double IMRPhenomX_Inspiral_Phase_22_Ansatz(double Mf, IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXPhaseCoefficients *pPhase)
{
    double phaseIN;

    // Assemble PN phase derivative series
  	phaseIN  = pPhase->dphi0; 																						    // f^{0/3}
  	phaseIN += pPhase->dphi1 	* powers_of_Mf->one_third; 								      // f^{1/3}
  	phaseIN += pPhase->dphi2 	* powers_of_Mf->two_thirds; 									    // f^{2/3}
  	phaseIN += pPhase->dphi3 	* Mf; 									        // f^{3/3}
  	phaseIN += pPhase->dphi4 	* powers_of_Mf->four_thirds; 							      // f^{4/3}
  	phaseIN += pPhase->dphi5 	* powers_of_Mf->five_thirds; 							      // f^{5/3}
  	phaseIN += pPhase->dphi6  * powers_of_Mf->two;													    // f^{6/3}
  	phaseIN += pPhase->dphi6L * powers_of_Mf->two * powers_of_Mf->log;			    // f^{6/3}, Log[f]
  	phaseIN += pPhase->dphi7  * powers_of_Mf->seven_thirds;									  // f^{7/3}
  	phaseIN += pPhase->dphi8  * powers_of_Mf->eight_thirds;									  // f^{8/3}
  	phaseIN += pPhase->dphi8L * powers_of_Mf->eight_thirds * powers_of_Mf->log;	// f^{8/3}
  	phaseIN += pPhase->dphi9  * powers_of_Mf->three;														// f^{9/3}
  	phaseIN += pPhase->dphi9L * powers_of_Mf->three * powers_of_Mf->log;				// f^{9/3}

  	// Add pseudo-PN Coefficient
  	phaseIN += ( 		pPhase->a0 * powers_of_Mf->eight_thirds
  								+ pPhase->a1 * powers_of_Mf->three
  								+ pPhase->a2 * powers_of_Mf->eight_thirds * powers_of_Mf->two_thirds
  								+ pPhase->a3 * powers_of_Mf->eight_thirds * powers_of_Mf->itself
  								+ pPhase->a4 * powers_of_Mf->eight_thirds * powers_of_Mf->four_thirds
  							);

  	phaseIN  = phaseIN * powers_of_Mf->m_eight_thirds * (5.0 / (128.0 * powers_of_lalpi.five_thirds));

    return phaseIN;
}

/**
 * Ansatz for the inspiral phase.
 * The TaylorF2 coefficients are defined elsewhere.
 */
static double IMRPhenomX_Inspiral_Phase_22_AnsatzInt(double Mf, IMRPhenomX_UsefulPowers *powers_of_Mf, IMRPhenomXPhaseCoefficients *pPhase)
{

  // Assemble PN phasing series
  //const double v    = powers_of_Mf->one_third * powers_of_lalpi.one_third; // v = (\pi M f)^{1/3}
  //const double logv = log(v);

  // Sum up phasing contributions
  double phasing    = 0.0;

  /* The PN Phasing series is normalised by: 3 / (128 * eta * pi^{5/3} ) */
  /* 0  */ phasing += pPhase->phi0;                                                   // f^{-5/3}, v = 0;  Newt.
  /* 1  */ phasing += pPhase->phi1  * powers_of_Mf->one_third;                        // f^{-4/3}, v = 1;  0.5PN
  /* 2  */ phasing += pPhase->phi2  * powers_of_Mf->two_thirds;                       // f^{-3/3}, v = 2;  1.0PN
  /* 3  */ phasing += pPhase->phi3  * Mf;                                             // f^{-2/3}, v = 3;  1.5PN
  /* 4  */ phasing += pPhase->phi4  * powers_of_Mf->four_thirds;                      // f^{-1/3}, v = 4;  2.0PN
  /* 5  */ phasing += pPhase->phi5  * powers_of_Mf->five_thirds;                      // f^{0},    v = 5;  2.5PN; phi_initial = - LAL_PI_4
  /* 5L */ phasing += pPhase->phi5L * powers_of_Mf->five_thirds * powers_of_Mf->log;  // f^{0},    v = 5;  2.5PN Log terms.
  /* 6  */ phasing += pPhase->phi6  * powers_of_Mf->two;                              // f^{+1/3}; v = 6;  3.0PN
  /* 6L */ phasing += pPhase->phi6L * powers_of_Mf->two * powers_of_Mf->log;          // f^{+1/3}; v = 6;  3.0PN Log terms.
  /* 7  */ phasing += pPhase->phi7  * powers_of_Mf->seven_thirds;                     // f^{+2/3}: v = 7;  3.5PN
  /* 8  */ phasing += pPhase->phi8  * powers_of_Mf->eight_thirds;                     // f^{+3/3}; v = 8;  4.0PN
  /* 8  */ phasing += pPhase->phi8L * powers_of_Mf->eight_thirds * powers_of_Mf->log; // f^{+3/3}; v = 8;  4.0PN Log terms.
  /* 9  */ phasing += pPhase->phi9  * powers_of_Mf->three;                            // f^{+4/3}; v = 9;  4.5PN
  /* 9  */ phasing += pPhase->phi9L * powers_of_Mf->three * powers_of_Mf->log;        // f^{+4/3}; v = 9;  4.5PN

  // Now add in the pseudo-PN Coefficients
  phasing += (  pPhase->sigma1 * powers_of_Mf->eight_thirds
              + pPhase->sigma2 * powers_of_Mf->three
              + pPhase->sigma3 * powers_of_Mf->one_third  * powers_of_Mf->three
              + pPhase->sigma4 * powers_of_Mf->two_thirds * powers_of_Mf->three
              + pPhase->sigma5 * powers_of_Mf->itself     * powers_of_Mf->three
            );

  // This completes the TaylorF2 PN phasing series
  phasing = phasing * pPhase->phiNorm * powers_of_Mf->m_five_thirds;

  /* Add initial phasing: -pi/4 */
  //phasing += pPhase->phi_initial;

  return phasing;
}
