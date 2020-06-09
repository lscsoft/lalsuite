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
 *
 * \file
 *
 * \brief Internal function for IMRPhenomX phenomenological waveform model, arXiv:2001.11412
 *
 */

#include "LALSimIMRPhenomX_intermediate.h"


/******************************* IMRPhenomX Amplitude Functions *******************************/

/******************************* Amplitude Functions: Merger   *******************************/

/*
	See section VI. B of arXiv:2001.11412

    Value of intermediate amplitude collocation vA point at f_2, defined in Table II of arXiv:2001.11412
   
    Effective spin parameterization used = StotR
*/
static double IMRPhenomX_Intermediate_Amp_22_vA(double eta, double S, double dchi, double delta, int IntAmpFlag){

  double eta2 = (eta*eta);
  double eta3 = (eta2*eta);

  double S2 = S*S;

  double noSpin, eqSpin, uneqSpin;

  switch ( IntAmpFlag )
	{
    case 104:
	{

		noSpin = (1.4873184918202145 + 1974.6112656679577*eta + 27563.641024162127*eta2 - 19837.908020966777*eta3)/(1. + 143.29004876335128*eta + 458.4097306093354*eta2);

		eqSpin = (S*(27.952730865904343 + eta*(-365.55631765202895 - 260.3494489873286*S) + 3.2646808851249016*S + 3011.446602208493*eta2*S - 19.38970173389662*S2 + eta3*(1612.2681322644232 - 6962.675551371755*S + 1486.4658089990298*S2)))/(12.647425554323242 - 10.540154508599963*S + 1.*S2);

		uneqSpin = dchi*delta*(-0.016404056649860943 - 296.473359655246*eta)*eta2;

  		break;
	}
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Intermediate_Amp_22_vA: IMRPhenomXIntermediateAmpVersion is not valid. Recommended flag is 104. \n");
    }
	}

  return (noSpin + eqSpin + uneqSpin);

}

/*	
	See section VI. B of arXiv:2001.11412. Collocation point for 5th order polynomial ansatz.

    Value of intermediate amplitude collocation point v2 at f_2, defined in Table I of arXiv:2001.11412
   
    Effective spin parameterization used = StotR
*/
static double IMRPhenomX_Intermediate_Amp_22_v2(double eta, double S, double dchi, double delta, int IntAmpFlag){

  double eta2 = (eta*eta);
  double eta3 = (eta2*eta);

  double S2 = S*S;

  double noSpin, eqSpin, uneqSpin;

  switch ( IntAmpFlag )
	{
		case 1043: // 1043 is used in IMRPhenomXHM
		case 105:
		{
      noSpin = (2.2436523786378983 + 2162.4749081764216*eta + 24460.158604784723*eta2 - 12112.140570900956*eta3)/(1. + 120.78623282522702*eta + 416.4179522274108*eta2);

      eqSpin = (S*(6.727511603827924 + eta2*(414.1400701039126 - 234.3754066885935*S) - 5.399284768639545*S
                  + eta*(-186.87972530996245 + 128.7402290554767*S)))/(3.24359204029217 - 3.975650468231452*S + 1.*S2);

      uneqSpin = dchi*delta*(-59.52510939953099 + 13.12679437100751*eta)*eta2;

  		break;
		}
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Intermediate_Amp_22_v2: IMRPhenomXIntermediateAmpVersion is not valid. Recommended flag is 104 (should not call *_Amp_22_v2). \n");
    }
	}

  return (noSpin + eqSpin + uneqSpin);

}

/*	
	See section VI. B of arXiv:2001.11412. Collocation point for 5th order polynomial ansatz.

    Value of intermediate amplitude collocation point v3 at f_3, defined in Table I of arXiv:2001.11412
   
    Effective spin parameterization used = StotR
*/
static double IMRPhenomX_Intermediate_Amp_22_v3(double eta, double S, double dchi, double delta, int IntAmpFlag){

  double eta2 = (eta*eta);
  double eta3 = (eta2*eta);

  double S2 = S*S;

  double noSpin, eqSpin, uneqSpin;

  switch ( IntAmpFlag )
	{
		case 1043: // 1043 is used in IMRPhenomXHM
		case 105:
		{
      noSpin = (1.195392410912163 + 1677.2558976605421*eta + 24838.37133975971*eta2 - 17277.938868280915*eta3)/(1. + 144.78606839716073*eta + 428.8155916011666*eta2);

      eqSpin = (S*(-2.1413952025647793 + 0.5719137940424858*S + eta*(46.61350006858767 + 0.40917927503842105*S - 11.526500209146906*S2)
      + 1.1833965566688387*S2 + eta2*(-84.82318288272965 - 34.90591158988979*S + 19.494962340530186*S2)))/(-1.4786392693666195 + 1.*S);

      uneqSpin = dchi*delta*(-333.7662575986524 + 532.2475589084717*eta)*eta3;

  		break;
		}
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Intermediate_Amp_22_v3: IMRPhenomXIntermediateAmpVersion is not valid. Recommended flag is 104 (should not call *_Amp_22_v3). \n");
    }
	}

  return (noSpin + eqSpin + uneqSpin);
}


/* Solves system of equations for 5th order polynomial ansatz. Section VI.B, Eq. 6.5 of arXiv:2001.11412. Note that \alpha in paper are \delta in code here. */
static double IMRPhenomX_Intermediate_Amp_22_delta0(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag)
{
  double retVal;

  switch ( IntAmpFlag )
	{
    case 1043:  //no left derivative; 1043 is used in IMRPhenomXHM
    {
      UNUSED double f12 = f1*f1;
      UNUSED double f13 = f12*f1;
      UNUSED double f14 = f13*f1;

      UNUSED double f22 = f2*f2;
      UNUSED double f23 = f22*f2;
      UNUSED double f24 = f23*f2;
      UNUSED double f25 = f24*f2;

      UNUSED double f32 = f3*f3;
      UNUSED double f33 = f32*f3;
      UNUSED double f34 = f33*f3;

      UNUSED double f42 = f4*f4;
      UNUSED double f43 = f42*f4;
      UNUSED double f44 = f43*f4;
      UNUSED double f45 = f44*f4;

      UNUSED double f1mf2 = f1-f2;
      UNUSED double f1mf3 = f1-f3;
      UNUSED double f1mf4 = f1-f4;
      UNUSED double f2mf3 = f2-f3;
      UNUSED double f2mf4 = f2-f4;
      UNUSED double f3mf4 = f3-f4;

      UNUSED double f1mf42 = f1mf4*f1mf4;
      UNUSED double f3mf42 = f3mf4*f3mf4;
      UNUSED double f2mf32 = f2mf3*f2mf3;
      UNUSED double f2mf42 = f2mf4*f2mf4;

      retVal = (-(d4*f1*f1mf2*f1mf3*f1mf4*f2*f2mf3*f2mf4*f3*f3mf4*f4) - f1*f1mf3*f1mf42*f3*f3mf42*f42*v2 +
      f24*(-(f1*f1mf42*f42*v3) + f33*(f42*v1 + f12*v4 - 2*f1*f4*v4) + f3*f4*(f43*v1 + 2*f13*v4 - 3*f12*f4*v4) - f32*(2*f43*v1 + f13*v4 - 3*f1*f42*v4)) +
      f2*f4*(f12*f1mf42*f43*v3 - f34*(f43*v1 + 2*f13*v4 - 3*f12*f4*v4) - f32*f4*(f44*v1 + 3*f14*v4 - 4*f13*f4*v4) + 2*f33*(f44*v1 + f14*v4 - 2*f12*f42*v4)) +
      f22*(-(f1*f1mf42*(2*f1 + f4)*f43*v3) + f3*f42*(f44*v1 + 3*f14*v4 - 4*f13*f4*v4) + f34*(2*f43*v1 + f13*v4 - 3*f1*f42*v4) - f33*(3*f44*v1 + f14*v4 - 4*f1*f43*v4)) +
      f23*(f1*f1mf42*(f1 + 2*f4)*f42*v3 - f34*(f42*v1 + f12*v4 - 2*f1*f4*v4) + f32*(3*f44*v1 + f14*v4 - 4*f1*f43*v4) - 2*f3*(f45*v1 + f14*f4*v4 - 2*f12*f43*v4)))/(f1mf2*f1mf3*f1mf42*f2mf3*f2mf42*f3mf42);
      break;
    }
    case 104:
		{

      UNUSED double f12 = f1*f1;
      UNUSED double f13 = f12*f1;
      UNUSED double f14 = f13*f1;

      UNUSED double f22 = f2*f2;
      UNUSED double f23 = f22*f2;
      UNUSED double f24 = f23*f2;

      UNUSED double f32 = f3*f3;
      UNUSED double f33 = f32*f3;
      UNUSED double f34 = f33*f3;

      UNUSED double f42 = f4*f4;
      UNUSED double f43 = f42*f4;
      UNUSED double f44 = f43*f4;

      UNUSED double f1mf2 = f1-f2;
      UNUSED double f1mf3 = f1-f3;
      UNUSED double f1mf4 = f1-f4;
      UNUSED double f2mf3 = f2-f3;
      UNUSED double f2mf4 = f2-f4;
      UNUSED double f3mf4 = f3-f4;

      UNUSED double f1mf22 = f1mf2*f1mf2;
      UNUSED double f1mf32 = f1mf3*f1mf3;
      UNUSED double f1mf42 = f1mf4*f1mf4;
      UNUSED double f2mf32 = f2mf3*f2mf3;
      UNUSED double f2mf42 = f2mf4*f2mf4;
      UNUSED double f3mf42 = f3mf4*f3mf4;
      UNUSED double f1mf43 = f1mf4*f1mf4*f1mf4;

      retVal = ((-(d4*f12*f1mf22*f1mf4*f2*f2mf4*f4) + d1*f1*f1mf2*f1mf4*f2*f2mf42*f42 + f42*(f2*f2mf42*(-4*f12 + 3*f1*f2 + 2*f1*f4 - f2*f4)*v1 + f12*f1mf43*v2) +
     f12*f1mf22*f2*(f1*f2 - 2*f1*f4 - 3*f2*f4 + 4*f42)*v4)/(f1mf22*f1mf43*f2mf42));

  		break;
		}
		case 105:
		{
      UNUSED double f12 = f1*f1;
      UNUSED double f13 = f12*f1;
      UNUSED double f14 = f13*f1;
      UNUSED double f15 = f14*f1;
      UNUSED double f16 = f15*f1;
      UNUSED double f17 = f16*f1;

      UNUSED double f22 = f2*f2;
      UNUSED double f23 = f22*f2;
      UNUSED double f24 = f23*f2;
      UNUSED double f25 = f24*f2;
      UNUSED double f26 = f25*f2;
      UNUSED double f27 = f26*f2;

      UNUSED double f32 = f3*f3;
      UNUSED double f33 = f32*f3;
      UNUSED double f34 = f33*f3;
      UNUSED double f35 = f34*f3;
      UNUSED double f36 = f35*f3;
      UNUSED double f37 = f36*f3;

      UNUSED double f42 = f4*f4;
      UNUSED double f43 = f42*f4;
      UNUSED double f44 = f43*f4;
      UNUSED double f45 = f44*f4;
      UNUSED double f46 = f45*f4;
      UNUSED double f47 = f46*f4;

      UNUSED double f1mf2 = f1-f2;
      UNUSED double f1mf3 = f1-f3;
      UNUSED double f1mf4 = f1-f4;
      UNUSED double f2mf3 = f2-f3;
      UNUSED double f2mf4 = f2-f4;
      UNUSED double f3mf4 = f3-f4;

      UNUSED double f1mf22 = f1mf2*f1mf2;
      UNUSED double f1mf32 = f1mf3*f1mf3;
      UNUSED double f1mf42 = f1mf4*f1mf4;
      UNUSED double f2mf32 = f2mf3*f2mf3;
      UNUSED double f2mf42 = f2mf4*f2mf4;
      UNUSED double f3mf42 = f3mf4*f3mf4;
      UNUSED double f1mf43 = f1mf4*f1mf4*f1mf4;

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
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Intermediate_Amp_22_v3: IMRPhenomXIntermediateAmpVersion is not valid. Recommended flag is 104.\n");
    }
	}

  return retVal;
}

/* Solves system of equations for 5th order polynomial ansatz. Section VI.B, Eq. 6.5 of arXiv:2001.11412. Note that \alpha in paper are \delta in code here. */
static double IMRPhenomX_Intermediate_Amp_22_delta1(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag)
{
  double retVal;

  switch ( IntAmpFlag )
  {
    case 1043:  //no left derivative; 1043 is used in IMRPhenomXHM
    {
      UNUSED double f12 = f1*f1;
      UNUSED double f13 = f12*f1;
      UNUSED double f14 = f13*f1;

      UNUSED double f22 = f2*f2;
      UNUSED double f23 = f22*f2;
      UNUSED double f24 = f23*f2;
      UNUSED double f25 = f24*f2;

      UNUSED double f32 = f3*f3;
      UNUSED double f33 = f32*f3;
      UNUSED double f34 = f33*f3;

      UNUSED double f42 = f4*f4;
      UNUSED double f43 = f42*f4;
      UNUSED double f44 = f43*f4;

      UNUSED double f1mf2 = f1-f2;
      UNUSED double f1mf3 = f1-f3;
      UNUSED double f1mf4 = f1-f4;
      UNUSED double f2mf3 = f2-f3;
      UNUSED double f2mf4 = f2-f4;
      UNUSED double f3mf4 = f3-f4;

      UNUSED double f1mf42 = f1mf4*f1mf4;
      UNUSED double f2mf42 = f2mf4*f2mf4;
      UNUSED double f3mf42 = f3mf4*f3mf4;

      UNUSED double v1mv2 = v1-v2;
      UNUSED double v2mv3 = v2-v3;
      UNUSED double v2mv4 = v2-v4;
      UNUSED double v1mv3 = v1-v3;
      UNUSED double v1mv4 = v1-v4;
      UNUSED double v3mv4 = v3-v4;

      retVal =(d4*f1mf4*f2mf4*f3mf4*(f1*f2*f3 + f2*f3*f4 + f1*(f2 + f3)*f4) + (f4*(f12*f1mf42*f43*v2mv3 + f34*(f43*v1mv2 +
        3*f12*f4*v2mv4 + 2*f13*(-v2 + v4)) + f32*f4*(f44*v1mv2 + 4*f13*f4*v2mv4 + 3*f14*(-v2 + v4)) + 2*f33*(f44*(-v1 + v2)
        + f14*v2mv4 + 2*f12*f42*(-v2 + v4)) + 2*f23*(f44*v1mv3 + f34*v1mv4 + 2*f12*f42*v3mv4 + 2*f32*f42*(-v1 + v4) + f14*(-v3 + v4))
        + f24*(3*f32*f4*v1mv4 + f43*(-v1 + v3) + 2*f13*v3mv4 + 2*f33*(-v1 + v4) + 3*f12*f4*(-v3 + v4)) + f22*f4*(4*f33*f4*v1mv4 + f44*(-v1 + v3)
        + 3*f14*v3mv4 + 3*f34*(-v1 + v4) + 4*f13*f4*(-v3 + v4))))/(f1mf2*f1mf3*f2mf3))/(f1mf42*f2mf42*f3mf42);

      break;
    }
    case 104:
    {

      UNUSED double f12 = f1*f1;
      UNUSED double f13 = f12*f1;
      UNUSED double f14 = f13*f1;

      UNUSED double f22 = f2*f2;
      UNUSED double f23 = f22*f2;
      UNUSED double f24 = f23*f2;

      UNUSED double f32 = f3*f3;
      UNUSED double f33 = f32*f3;
      UNUSED double f34 = f33*f3;

      UNUSED double f42 = f4*f4;
      UNUSED double f43 = f42*f4;
      UNUSED double f44 = f43*f4;

      UNUSED double f1mf2 = f1-f2;
      UNUSED double f1mf3 = f1-f3;
      UNUSED double f1mf4 = f1-f4;
      UNUSED double f2mf3 = f2-f3;
      UNUSED double f2mf4 = f2-f4;
      UNUSED double f3mf4 = f3-f4;

      UNUSED double f1mf22 = f1mf2*f1mf2;
      UNUSED double f1mf32 = f1mf3*f1mf3;
      UNUSED double f1mf42 = f1mf4*f1mf4;
      UNUSED double f2mf32 = f2mf3*f2mf3;
      UNUSED double f2mf42 = f2mf4*f2mf4;
      UNUSED double f3mf42 = f3mf4*f3mf4;
      UNUSED double f1mf43 = f1mf4*f1mf4*f1mf4;

      retVal = ((d4*f1*f1mf22*f1mf4*f2mf4*(2*f2*f4 + f1*(f2 + f4)) + f4*(-(d1*f1mf2*f1mf4*f2mf42*(2*f1*f2 + (f1 + f2)*f4)) -
        2*f1*(f44*(v1 - v2) + 3*f24*(v1 - v4) + f14*(v2 - v4) + 4*f23*f4*(-v1 + v4)
        + 2*f13*f4*(-v2 + v4) + f1*(2*f43*(-v1 + v2) + 6*f22*f4*(v1 - v4) + 4*f23*(-v1 + v4)))))/(f1mf22*f1mf43*f2mf42));


      break;
    }
    case 105:
    {

      UNUSED double f12 = f1*f1;
      UNUSED double f13 = f12*f1;
      UNUSED double f14 = f13*f1;
      UNUSED double f15 = f14*f1;
      UNUSED double f16 = f15*f1;
      UNUSED double f17 = f16*f1;

      UNUSED double f22 = f2*f2;
      UNUSED double f23 = f22*f2;
      UNUSED double f24 = f23*f2;
      UNUSED double f25 = f24*f2;
      UNUSED double f26 = f25*f2;
      UNUSED double f27 = f26*f2;

      UNUSED double f32 = f3*f3;
      UNUSED double f33 = f32*f3;
      UNUSED double f34 = f33*f3;
      UNUSED double f35 = f34*f3;
      UNUSED double f36 = f35*f3;
      UNUSED double f37 = f36*f3;

      UNUSED double f42 = f4*f4;
      UNUSED double f43 = f42*f4;
      UNUSED double f44 = f43*f4;
      UNUSED double f45 = f44*f4;
      UNUSED double f46 = f45*f4;
      UNUSED double f47 = f46*f4;

      UNUSED double f1mf2 = f1-f2;
      UNUSED double f1mf3 = f1-f3;
      UNUSED double f1mf4 = f1-f4;
      UNUSED double f2mf3 = f2-f3;
      UNUSED double f2mf4 = f2-f4;
      UNUSED double f3mf4 = f3-f4;

      UNUSED double f1mf22 = f1mf2*f1mf2;
      UNUSED double f1mf32 = f1mf3*f1mf3;
      UNUSED double f1mf42 = f1mf4*f1mf4;
      UNUSED double f2mf32 = f2mf3*f2mf3;
      UNUSED double f2mf42 = f2mf4*f2mf4;
      UNUSED double f3mf42 = f3mf4*f3mf4;
      UNUSED double f1mf43 = f1mf4*f1mf4*f1mf4;

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
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Intermediate_Amp_22_v3: IMRPhenomXIntermediateAmpVersion is not valid. Recommended flag is 104.\n");
    }
  }

  return retVal;
}

/* Solves system of equations for 5th order polynomial ansatz. Section VI.B, Eq. 6.5 of arXiv:2001.11412. Note that \alpha in paper are \delta in code here. */
static double IMRPhenomX_Intermediate_Amp_22_delta2(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag)
{
  double retVal;

  switch ( IntAmpFlag )
  {
    case 1043:  //no left derivative: v1, v2, v3, v4, d4; 1043 is used in IMRPhenomXHM
    {
      UNUSED double f12 = f1*f1;
      UNUSED double f13 = f12*f1;
      UNUSED double f14 = f13*f1;

      UNUSED double f22 = f2*f2;
      UNUSED double f23 = f22*f2;
      UNUSED double f24 = f23*f2;
      UNUSED double f25 = f24*f2;
      UNUSED double f26 = f25*f2;

      UNUSED double f32 = f3*f3;
      UNUSED double f33 = f32*f3;
      UNUSED double f34 = f33*f3;

      UNUSED double f42 = f4*f4;
      UNUSED double f43 = f42*f4;
      UNUSED double f44 = f43*f4;
      UNUSED double f46 = f44*f42;

      UNUSED double f1mf2 = f1-f2;
      UNUSED double f1mf3 = f1-f3;
      UNUSED double f1mf4 = f1-f4;
      UNUSED double f2mf3 = f2-f3;
      UNUSED double f2mf4 = f2-f4;
      UNUSED double f3mf4 = f3-f4;

      UNUSED double f1mf42 = f1mf4*f1mf4;
      UNUSED double f2mf42 = f2mf4*f2mf4;
      UNUSED double f3mf42 = f3mf4*f3mf4;

      UNUSED double v1mv2 = v1-v2;
      UNUSED double v2mv3 = v2-v3;
      UNUSED double v2mv4 = v2-v4;
      UNUSED double v1mv3 = v1-v3;
      UNUSED double v1mv4 = v1-v4;
      UNUSED double v3mv4 = v3-v4;

      retVal = (-(d4*f1mf2*f1mf3*f1mf4*f2mf3*f2mf4*f3mf4*(f3*f4 + f2*(f3 + f4) + f1*(f2 + f3 + f4))) - 2*f34*f43*v1 + 3*f33*f44*v1 - f3*f46*v1 - f14*f33*v2 + f13*f34*v2 + 3*f14*f3*f42*v2 - 3*f1*f34*f42*v2 -
      2*f14*f43*v2 - 4*f13*f3*f43*v2 + 4*f1*f33*f43*v2 + 2*f34*f43*v2 + 3*f13*f44*v2 - 3*f33*f44*v2 - f1*f46*v2 + f3*f46*v2 + 2*f14*f43*v3 - 3*f13*f44*v3 + f1*f46*v3 +
      f2*f42*(f44*v1mv3 + 3*f34*v1mv4 - 4*f33*f4*v1mv4 - 3*f14*v3mv4 + 4*f13*f4*v3mv4) + f24*(2*f43*v1mv3 + f33*v1mv4 - 3*f3*f42*v1mv4 - f13*v3mv4 + 3*f1*f42*v3mv4) +
      f23*(-3*f44*v1mv3 - f34*v1mv4 + 4*f3*f43*v1mv4 + f14*v3mv4 - 4*f1*f43*v3mv4) + f14*f33*v4 - f13*f34*v4 - 3*f14*f3*f42*v4 + 3*f1*f34*f42*v4 + 4*f13*f3*f43*v4 - 4*f1*f33*f43*v4)/
      (f1mf2*f1mf3*f1mf42*f2mf3*f2mf42*f3mf42);
      break;
    }
    case 104:
    {
      UNUSED double f12 = f1*f1;
      UNUSED double f13 = f12*f1;
      UNUSED double f14 = f13*f1;
      UNUSED double f15 = f14*f1;

      UNUSED double f22 = f2*f2;
      UNUSED double f23 = f22*f2;
      UNUSED double f24 = f23*f2;

      UNUSED double f32 = f3*f3;
      UNUSED double f33 = f32*f3;
      UNUSED double f34 = f33*f3;

      UNUSED double f42 = f4*f4;
      UNUSED double f43 = f42*f4;
      UNUSED double f44 = f43*f4;
      UNUSED double f45 = f44*f4;

      UNUSED double f1mf2 = f1-f2;
      UNUSED double f1mf3 = f1-f3;
      UNUSED double f1mf4 = f1-f4;
      UNUSED double f2mf3 = f2-f3;
      UNUSED double f2mf4 = f2-f4;
      UNUSED double f3mf4 = f3-f4;

      UNUSED double f1mf22 = f1mf2*f1mf2;
      UNUSED double f1mf32 = f1mf3*f1mf3;
      UNUSED double f1mf42 = f1mf4*f1mf4;
      UNUSED double f2mf32 = f2mf3*f2mf3;
      UNUSED double f2mf42 = f2mf4*f2mf4;
      UNUSED double f3mf42 = f3mf4*f3mf4;
      UNUSED double f1mf43 = f1mf4*f1mf4*f1mf4;

      retVal = ((-(d4*f1mf22*f1mf4*f2mf4*(f12 + f2*f4 + 2*f1*(f2 + f4))) + d1*f1mf2*f1mf4*f2mf42*(f1*f2 + 2*(f1 + f2)*f4 + f42)
      - 4*f12*f23*v1 + 3*f1*f24*v1 - 4*f1*f23*f4*v1 + 3*f24*f4*v1 + 12*f12*f2*f42*v1 -
     4*f23*f42*v1 - 8*f12*f43*v1 + f1*f44*v1 + f45*v1 + f15*v2 + f14*f4*v2 - 8*f13*f42*v2 + 8*f12*f43*v2 - f1*f44*v2 - f45*v2 -
     f1mf22*(f13 + f2*(3*f2 - 4*f4)*f4 + f12*(2*f2 + f4) + f1*(3*f2 - 4*f4)*(f2 + 2*f4))*v4)/(f1mf22*f1mf43*f2mf42));

      break;
    }
    case 105:
    {
      UNUSED double f12 = f1*f1;
      UNUSED double f13 = f12*f1;
      UNUSED double f14 = f13*f1;
      UNUSED double f15 = f14*f1;
      UNUSED double f16 = f15*f1;
      UNUSED double f17 = f16*f1;

      UNUSED double f22 = f2*f2;
      UNUSED double f23 = f22*f2;
      UNUSED double f24 = f23*f2;
      UNUSED double f25 = f24*f2;
      UNUSED double f26 = f25*f2;
      UNUSED double f27 = f26*f2;

      UNUSED double f32 = f3*f3;
      UNUSED double f33 = f32*f3;
      UNUSED double f34 = f33*f3;
      UNUSED double f35 = f34*f3;
      UNUSED double f36 = f35*f3;
      UNUSED double f37 = f36*f3;

      UNUSED double f42 = f4*f4;
      UNUSED double f43 = f42*f4;
      UNUSED double f44 = f43*f4;
      UNUSED double f45 = f44*f4;
      UNUSED double f46 = f45*f4;
      UNUSED double f47 = f46*f4;

      UNUSED double f1mf2 = f1-f2;
      UNUSED double f1mf3 = f1-f3;
      UNUSED double f1mf4 = f1-f4;
      UNUSED double f2mf3 = f2-f3;
      UNUSED double f2mf4 = f2-f4;
      UNUSED double f3mf4 = f3-f4;

      UNUSED double f1mf22 = f1mf2*f1mf2;
      UNUSED double f1mf32 = f1mf3*f1mf3;
      UNUSED double f1mf42 = f1mf4*f1mf4;
      UNUSED double f2mf32 = f2mf3*f2mf3;
      UNUSED double f2mf42 = f2mf4*f2mf4;
      UNUSED double f3mf42 = f3mf4*f3mf4;
      UNUSED double f1mf43 = f1mf4*f1mf4*f1mf4;

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
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Intermediate_Amp_22_v3: IMRPhenomXIntermediateAmpVersion is not valid. Recommended flag is 104.\n");
    }
  }

return retVal;
}

/* Solves system of equations for 5th order polynomial ansatz. Section VI.B, Eq. 6.5 of arXiv:2001.11412. Note that \alpha in paper are \delta in code here. */
static double IMRPhenomX_Intermediate_Amp_22_delta3(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag)
{
  double retVal;

  switch ( IntAmpFlag )
  {
    case 1043:  //no left derivative: v1, v2, v3, v4, d4; 1043 is used in IMRPhenomXHM
    {
      UNUSED double f12 = f1*f1;
      UNUSED double f13 = f12*f1;
      UNUSED double f14 = f13*f1;

      UNUSED double f22 = f2*f2;
      UNUSED double f23 = f22*f2;
      UNUSED double f24 = f23*f2;

      UNUSED double f32 = f3*f3;
      UNUSED double f34 = f32*f32;

      UNUSED double f42 = f4*f4;
      UNUSED double f43 = f42*f4;
      UNUSED double f44 = f43*f4;
      UNUSED double f45 = f44*f4;

      UNUSED double f1mf2 = f1-f2;
      UNUSED double f1mf3 = f1-f3;
      UNUSED double f1mf4 = f1-f4;
      UNUSED double f2mf3 = f2-f3;
      UNUSED double f2mf4 = f2-f4;
      UNUSED double f3mf4 = f3-f4;

      UNUSED double f1mf42 = f1mf4*f1mf4;
      UNUSED double f2mf42 = f2mf4*f2mf4;
      UNUSED double f3mf42 = f3mf4*f3mf4;

      UNUSED double v1mv2 = v1-v2;
      UNUSED double v2mv3 = v2-v3;
      UNUSED double v2mv4 = v2-v4;
      UNUSED double v1mv3 = v1-v3;
      UNUSED double v1mv4 = v1-v4;
      UNUSED double v3mv4 = v3-v4;

      retVal = (d4*f1mf2*f1mf3*f1mf4*f2mf3*f2mf4*f3mf4*(f1 + f2 + f3 + f4) + f34*f42*v1 - 3*f32*f44*v1 + 2*f3*f45*v1 + f14*f32*v2 - f12*f34*v2 - 2*f14*f3*f4*v2 + 2*f1*f34*f4*v2 + f14*f42*v2 - f34*f42*v2 +
      4*f12*f3*f43*v2 - 4*f1*f32*f43*v2 - 3*f12*f44*v2 + 3*f32*f44*v2 + 2*f1*f45*v2 - 2*f3*f45*v2 - f14*f42*v3 + 3*f12*f44*v3 - 2*f1*f45*v3 +
      f24*(-(f42*v1mv3) - f32*v1mv4 + 2*f3*f4*v1mv4 + f12*v3mv4 - 2*f1*f4*v3mv4) - 2*f2*f4*(f44*v1mv3 + f34*v1mv4 - 2*f32*f42*v1mv4 - f14*v3mv4 + 2*f12*f42*v3mv4) +
      f22*(3*f44*v1mv3 + f34*v1mv4 - 4*f3*f43*v1mv4 - f14*v3mv4 + 4*f1*f43*v3mv4) - f14*f32*v4 + f12*f34*v4 + 2*f14*f3*f4*v4 - 2*f1*f34*f4*v4 - 4*f12*f3*f43*v4 + 4*f1*f32*f43*v4)/
      (f1mf2*f1mf3*f1mf42*f2mf3*f2mf42*f3mf42);
      break;
    }
    case 104:
    {

      UNUSED double f12 = f1*f1;
      UNUSED double f13 = f12*f1;
      UNUSED double f14 = f13*f1;

      UNUSED double f22 = f2*f2;
      UNUSED double f23 = f22*f2;
      UNUSED double f24 = f23*f2;

      UNUSED double f32 = f3*f3;
      UNUSED double f33 = f32*f3;
      UNUSED double f34 = f33*f3;

      UNUSED double f42 = f4*f4;
      UNUSED double f43 = f42*f4;
      UNUSED double f44 = f43*f4;

      UNUSED double f1mf2 = f1-f2;
      UNUSED double f1mf3 = f1-f3;
      UNUSED double f1mf4 = f1-f4;
      UNUSED double f2mf3 = f2-f3;
      UNUSED double f2mf4 = f2-f4;
      UNUSED double f3mf4 = f3-f4;

      UNUSED double f1mf22 = f1mf2*f1mf2;
      UNUSED double f1mf32 = f1mf3*f1mf3;
      UNUSED double f1mf42 = f1mf4*f1mf4;
      UNUSED double f2mf32 = f2mf3*f2mf3;
      UNUSED double f2mf42 = f2mf4*f2mf4;
      UNUSED double f3mf42 = f3mf4*f3mf4;
      UNUSED double f1mf43 = f1mf4*f1mf4*f1mf4;

      retVal = ((d4*f1mf22*f1mf4*f2mf4*(2*f1 + f2 + f4) - d1*f1mf2*f1mf4*f2mf42*(f1 + f2 + 2*f4)
                + 2*(f44*(-v1 + v2) + 2*f12*f2mf42*(v1 - v4) + 2*f22*f42*(v1 - v4)
                + 2*f13*f4*(v2 - v4) + f24*(-v1 + v4) + f14*(-v2 + v4) + 2*f1*f4*(f42*(v1 - v2) + f22*(v1 - v4) + 2*f2*f4*(-v1 + v4)))) / (f1mf22*f1mf43*f2mf42));


      break;
    }
    case 105:
    {
      UNUSED double f12 = f1*f1;
      UNUSED double f13 = f12*f1;
      UNUSED double f14 = f13*f1;
      UNUSED double f15 = f14*f1;
      UNUSED double f16 = f15*f1;
      UNUSED double f17 = f16*f1;

      UNUSED double f22 = f2*f2;
      UNUSED double f23 = f22*f2;
      UNUSED double f24 = f23*f2;
      UNUSED double f25 = f24*f2;
      UNUSED double f26 = f25*f2;
      UNUSED double f27 = f26*f2;

      UNUSED double f32 = f3*f3;
      UNUSED double f33 = f32*f3;
      UNUSED double f34 = f33*f3;
      UNUSED double f35 = f34*f3;
      UNUSED double f36 = f35*f3;
      UNUSED double f37 = f36*f3;

      UNUSED double f42 = f4*f4;
      UNUSED double f43 = f42*f4;
      UNUSED double f44 = f43*f4;
      UNUSED double f45 = f44*f4;
      UNUSED double f46 = f45*f4;
      UNUSED double f47 = f46*f4;

      UNUSED double f1mf2 = f1-f2;
      UNUSED double f1mf3 = f1-f3;
      UNUSED double f1mf4 = f1-f4;
      UNUSED double f2mf3 = f2-f3;
      UNUSED double f2mf4 = f2-f4;
      UNUSED double f3mf4 = f3-f4;

      UNUSED double f1mf22 = f1mf2*f1mf2;
      UNUSED double f1mf32 = f1mf3*f1mf3;
      UNUSED double f1mf42 = f1mf4*f1mf4;
      UNUSED double f2mf32 = f2mf3*f2mf3;
      UNUSED double f2mf42 = f2mf4*f2mf4;
      UNUSED double f3mf42 = f3mf4*f3mf4;
      UNUSED double f1mf43 = f1mf4*f1mf4*f1mf4;

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
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Intermediate_Amp_22_v3: IMRPhenomXIntermediateAmpVersion is not valid. Recommended flag is 104.\n");
    }
  }

  return retVal;
}

/* Solves system of equations for 5th order polynomial ansatz. Section VI.B, Eq. 6.5 of arXiv:2001.11412. Note that \alpha in paper are \delta in code here. */
static double IMRPhenomX_Intermediate_Amp_22_delta4(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag)
{
  double retVal;

  switch ( IntAmpFlag )
  {
    case 1043:  //no left derivative: v1, v2, v3, v4, d4; 1043 is used in IMRPhenomXHM
    {
      UNUSED double f12 = f1*f1;
      UNUSED double f13 = f12*f1;

      UNUSED double f22 = f2*f2;
      UNUSED double f23 = f22*f2;
      UNUSED double f24 = f23*f2;

      UNUSED double f32 = f3*f3;
      UNUSED double f33 = f32*f3;

      UNUSED double f42 = f4*f4;
      UNUSED double f43 = f42*f4;
      UNUSED double f44 = f43*f4;

      UNUSED double f1mf2 = f1-f2;
      UNUSED double f1mf3 = f1-f3;
      UNUSED double f1mf4 = f1-f4;
      UNUSED double f2mf3 = f2-f3;
      UNUSED double f2mf4 = f2-f4;
      UNUSED double f3mf4 = f3-f4;

      UNUSED double f1mf42 = f1mf4*f1mf4;
      UNUSED double f2mf42 = f2mf4*f2mf4;
      UNUSED double f3mf42 = f3mf4*f3mf4;

      UNUSED double v1mv2 = v1-v2;
      UNUSED double v2mv3 = v2-v3;
      UNUSED double v2mv4 = v2-v4;
      UNUSED double v1mv3 = v1-v3;
      UNUSED double v1mv4 = v1-v4;
      UNUSED double v3mv4 = v3-v4;

      retVal = (-(d4*f1mf2*f1mf3*f1mf4*f2mf3*f2mf4*f3mf4) - f33*f42*v1 + 2*f32*f43*v1 - f3*f44*v1 - f13*f32*v2 + f12*f33*v2 + 2*f13*f3*f4*v2 - 2*f1*f33*f4*v2 - f13*f42*v2 - 3*f12*f3*f42*v2 + 3*f1*f32*f42*v2 +
      f33*f42*v2 + 2*f12*f43*v2 - 2*f32*f43*v2 - f1*f44*v2 + f3*f44*v2 + f13*f42*v3 - 2*f12*f43*v3 + f1*f44*v3 + f23*(f42*v1mv3 + f32*v1mv4 - 2*f3*f4*v1mv4 - f12*v3mv4 + 2*f1*f4*v3mv4) +
      f2*f4*(f43*v1mv3 + 2*f33*v1mv4 - 3*f32*f4*v1mv4 - 2*f13*v3mv4 + 3*f12*f4*v3mv4) + f22*(-2*f43*v1mv3 - f33*v1mv4 + 3*f3*f42*v1mv4 + f13*v3mv4 - 3*f1*f42*v3mv4) + f13*f32*v4 - f12*f33*v4 -
      2*f13*f3*f4*v4 + 2*f1*f33*f4*v4 + 3*f12*f3*f42*v4 - 3*f1*f32*f42*v4)/(f1mf2*f1mf3*f1mf42*f2mf3*f2mf42*f3mf42);
      break;
    }
    case 104:
    {
      UNUSED double f12 = f1*f1;
      UNUSED double f13 = f12*f1;
      UNUSED double f14 = f13*f1;

      UNUSED double f22 = f2*f2;
      UNUSED double f23 = f22*f2;
      UNUSED double f24 = f23*f2;

      UNUSED double f32 = f3*f3;
      UNUSED double f33 = f32*f3;
      UNUSED double f34 = f33*f3;

      UNUSED double f42 = f4*f4;
      UNUSED double f43 = f42*f4;
      UNUSED double f44 = f43*f4;

      UNUSED double f1mf2 = f1-f2;
      UNUSED double f1mf3 = f1-f3;
      UNUSED double f1mf4 = f1-f4;
      UNUSED double f2mf3 = f2-f3;
      UNUSED double f2mf4 = f2-f4;
      UNUSED double f3mf4 = f3-f4;

      UNUSED double f1mf22 = f1mf2*f1mf2;
      UNUSED double f1mf32 = f1mf3*f1mf3;
      UNUSED double f1mf42 = f1mf4*f1mf4;
      UNUSED double f2mf32 = f2mf3*f2mf3;
      UNUSED double f2mf42 = f2mf4*f2mf4;
      UNUSED double f3mf42 = f3mf4*f3mf4;
      UNUSED double f1mf43 = f1mf4*f1mf4*f1mf4;

      retVal = ((-(d4*f1mf22*f1mf4*f2mf4) + d1*f1mf2*f1mf4*f2mf42 - 3*f1*f22*v1 + 2*f23*v1 + 6*f1*f2*f4*v1 - 3*f22*f4*v1
              - 3*f1*f42*v1 + f43*v1 + f13*v2 - 3*f12*f4*v2 + 3*f1*f42*v2 - f43*v2 - f1mf22*(f1 + 2*f2 - 3*f4)*v4)/(f1mf22*f1mf43*f2mf42));

      break;
    }
    case 105:
    {
      UNUSED double f12 = f1*f1;
      UNUSED double f13 = f12*f1;
      UNUSED double f14 = f13*f1;
      UNUSED double f15 = f14*f1;
      UNUSED double f16 = f15*f1;
      UNUSED double f17 = f16*f1;

      UNUSED double f22 = f2*f2;
      UNUSED double f23 = f22*f2;
      UNUSED double f24 = f23*f2;
      UNUSED double f25 = f24*f2;
      UNUSED double f26 = f25*f2;
      UNUSED double f27 = f26*f2;

      UNUSED double f32 = f3*f3;
      UNUSED double f33 = f32*f3;
      UNUSED double f34 = f33*f3;
      UNUSED double f35 = f34*f3;
      UNUSED double f36 = f35*f3;
      UNUSED double f37 = f36*f3;

      UNUSED double f42 = f4*f4;
      UNUSED double f43 = f42*f4;
      UNUSED double f44 = f43*f4;
      UNUSED double f45 = f44*f4;
      UNUSED double f46 = f45*f4;
      UNUSED double f47 = f46*f4;

      UNUSED double f1mf2 = f1-f2;
      UNUSED double f1mf3 = f1-f3;
      UNUSED double f1mf4 = f1-f4;
      UNUSED double f2mf3 = f2-f3;
      UNUSED double f2mf4 = f2-f4;
      UNUSED double f3mf4 = f3-f4;

      UNUSED double f1mf22 = f1mf2*f1mf2;
      UNUSED double f1mf32 = f1mf3*f1mf3;
      UNUSED double f1mf42 = f1mf4*f1mf4;
      UNUSED double f2mf32 = f2mf3*f2mf3;
      UNUSED double f2mf42 = f2mf4*f2mf4;
      UNUSED double f3mf42 = f3mf4*f3mf4;
      UNUSED double f1mf43 = f1mf4*f1mf4*f1mf4;

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
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Intermediate_Amp_22_v3: IMRPhenomXIntermediateAmpVersion is not valid. Recommended flag is 104.\n");
    }
  }

  return retVal;
}

/* Solves system of equations for 5th order polynomial ansatz. Section VI.B, Eq. 6.5 of arXiv:2001.11412. Note that \alpha in paper are \delta in code here. */
static double IMRPhenomX_Intermediate_Amp_22_delta5(double d1, double d4, double v1, double v2, double v3, double v4, double f1, double f2, double f3, double f4, int IntAmpFlag)
{
  double retVal;

  switch ( IntAmpFlag )
  {
    case 1043:  //no left derivative: v1, v2, v3, v4, d4; 1043 is used in IMRPhenomXHM
    {
      retVal = 0.;
      break;
    }
    case 104:
    {
      retVal = 0.0;
      break;
    }
    case 105:
    {
      UNUSED double f12 = f1*f1;
      UNUSED double f13 = f12*f1;
      UNUSED double f14 = f13*f1;
      UNUSED double f15 = f14*f1;
      UNUSED double f16 = f15*f1;
      UNUSED double f17 = f16*f1;

      UNUSED double f22 = f2*f2;
      UNUSED double f23 = f22*f2;
      UNUSED double f24 = f23*f2;
      UNUSED double f25 = f24*f2;
      UNUSED double f26 = f25*f2;
      UNUSED double f27 = f26*f2;

      UNUSED double f32 = f3*f3;
      UNUSED double f33 = f32*f3;
      UNUSED double f34 = f33*f3;
      UNUSED double f35 = f34*f3;
      UNUSED double f36 = f35*f3;
      UNUSED double f37 = f36*f3;

      UNUSED double f42 = f4*f4;
      UNUSED double f43 = f42*f4;
      UNUSED double f44 = f43*f4;
      UNUSED double f45 = f44*f4;
      UNUSED double f46 = f45*f4;
      UNUSED double f47 = f46*f4;

      UNUSED double f1mf2 = f1-f2;
      UNUSED double f1mf3 = f1-f3;
      UNUSED double f1mf4 = f1-f4;
      UNUSED double f2mf3 = f2-f3;
      UNUSED double f2mf4 = f2-f4;
      UNUSED double f3mf4 = f3-f4;

      UNUSED double f1mf22 = f1mf2*f1mf2;
      UNUSED double f1mf32 = f1mf3*f1mf3;
      UNUSED double f1mf42 = f1mf4*f1mf4;
      UNUSED double f2mf32 = f2mf3*f2mf3;
      UNUSED double f2mf42 = f2mf4*f2mf4;
      UNUSED double f3mf42 = f3mf4*f3mf4;
      UNUSED double f1mf43 = f1mf4*f1mf4*f1mf4;

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
        XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Intermediate_Amp_22_v3: IMRPhenomXIntermediateAmpVersion is not valid. Recommended flag is 104.\n");
      }
    }

    return retVal;
}

/* Solves system of equations for 5th order polynomial ansatz. Section VI.B, Eq. 6.5 of arXiv:2001.11412. Note that \alpha in paper are \delta in code here. */
static double IMRPhenomX_Intermediate_Amp_22_Ansatz(double ff, IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXWaveformStruct *pWF, IMRPhenomXAmpCoefficients *pAmp)
{
	double a0         = pAmp->delta0;
    double a1         = pAmp->delta1;
    double a2         = pAmp->delta2;
    double a3         = pAmp->delta3;
    double a4         = pAmp->delta4;
    double a5         = pAmp->delta5;

    double ff7o6      = powers_of_f->seven_sixths;

    int IntAmpFlag    = pWF->IMRPhenomXIntermediateAmpVersion;

    double polynomial;

    switch ( IntAmpFlag )
  	{
      case 1043: //1043 is used in IMRPhenomXHM
      case 104:
  		{
        //polynomial        = a0 + a1*ff + a2*ff2 + a3*ff3 + a4*ff4;
        polynomial = a0 + ff*(a1 + ff*(a2 + ff*(a3 + ff*a4)));

    		break;
  		}
  		case 105:
  		{
        //polynomial        = a0 + a1*ff + a2*ff2 + a3*ff3 + a4*ff4 + a5*ff5;
        polynomial = a0 + ff*(a1 + ff*(a2 + ff*(a3 + ff*(a4 + ff*a5))));

    		break;
  		}
      default:
      {
        XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Intermediate_Amp_22_Ansatz: IMRPhenomXIntermediateAmpVersion is not valid. Recommended flag is 104.\n");

        ff7o6      = 0.0;
        polynomial = 1.0;

        break;
      }
  	}

		/* Return a polynomial embedded in the background from the merger */
		return ff7o6 / polynomial;
}

/******************************* Phase Functions: Merger   *******************************/

/* Intermediate phase collocation point v3. See Section VII.B of arXiv:2001.11412.  */
UNUSED static double IMRPhenomX_Intermediate_Phase_22_v3(double eta, double S, double dchi, double delta, int IntPhaseFlag){

	double eta2  = eta*eta;
	double eta3  = eta2*eta;
  double eta4  = eta3*eta;

	double S2    = S*S;
	double S3    = S2*S;

	double noSpin, eqSpin, uneqSpin;

	switch ( IntPhaseFlag )
	{
		case 104: 			/* Canonical, 4 coefficients */
		{
      noSpin = (-85.65000000000003 - 413.9060603156141*eta - 4784.537141069749*eta2 + 1490.1208098409852*eta3 + 10040.191767983953*eta4)/(1. + 5.779247659216663*eta + 69.57046710291102*eta2);

      eqSpin = (S*(4.692726457075735 - 7.038980703663575*S + eta2*(87.92695072054693 - 53.72849910673372*S - 40.89299875429382*S2) + 1.809693754991601*S2 + eta*(-13.873355143236994 + 11.141463520702901*S + 6.450651190707574*S2) + eta3*(-214.15302815157358 + 165.87430354875272*S + 53.19430999626736*S2 - 2.8784864337832756*S3) + 0.18008416467045393*S3))/(-0.19174821584528354 + 0.22339091810747033*S + 1.*eta3*S - 0.06756256922190347*S2);

      uneqSpin = dchi*delta*eta3*(2989.296934478898 - 8693.785422139264*eta + 230.93311202289289*S);


  		break;
		}
		case 105:				/* Canonical, 5 coefficients */
		{
      noSpin = (-85.31300000000003 - 10937.096332153926*eta + 24894.22839497624*eta2 - 41966.01695284205*eta3 + 67157.48707167112*eta4)/(1. + 132.44486857622041*eta);

      eqSpin = (S*(-50.47357379240076 + 73.09727082201734*S + eta3*(3596.7780512820414 - 2213.251793086697*S - 2751.2131780622212*S2) + eta*(156.5226913355741 - 58.19427796874254*S - 243.59756424149438*S2) - 14.822551015459917*S2 + eta2*(-1399.6118302073216 + 505.3549556613735*S + 1522.6588404541071*S2)))/(1.9694717754728777 - 2.272396692525365*S + 1.*S2);

      uneqSpin = dchi*delta*eta3*(3453.778014556621 + eta*(-10048.103462582065 + 1757.9089741373239*S));

  		break;
		}
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Intermediate_Phase_22_v3: IMRPhenomXIntermediatePhaseVersion is not valid. Recommended flag is 104.\n");
      break;
    }
	}

	return (noSpin + eqSpin + uneqSpin);
}

/* Intermediate phase collocation point v2. See Section VII.B of arXiv:2001.11412.  */
static double IMRPhenomX_Intermediate_Phase_22_v2(double eta, double S, double dchi, double delta, int IntPhaseFlag){

	double eta2  = eta*eta;
	double eta3  = eta2*eta;

	double S2    = S*S;
	double S3    = S2*S;
	double S4    = S3*S;

	double noSpin, eqSpin, uneqSpin;

	switch ( IntPhaseFlag )
	{
		case 104: 			/* Canonical, 4 coefficients */
		{
      noSpin = (-84.09400000000004 - 1782.8025405571802*eta + 5384.38721936653*eta2)/(1. + 28.515617312596103*eta + 12.404097877099353*eta2);

      eqSpin = (S*(22.5665046165141 - 39.94907120140026*S + 4.668251961072*S2 + 12.648584361431245*S3 + eta2*(-298.7528127869681 + 14.228745354543983*S + 398.1736953382448*S2 + 506.94924905801673*S3 - 626.3693674479357*S4) - 5.360554789680035*S4 + eta*(152.07900889608595 - 121.70617732909962*S2 - 169.36311036967322*S3 + 182.40064586992762*S4)))/(-1.1571201220629341 + 1.*S);

      uneqSpin = dchi*delta*eta3*(5357.695021063607 - 15642.019339339662*eta + 674.8698102366333*S);

  		break;
		}
		case 105:				/* Canonical, 5 coefficients */
		{
      noSpin = (-82.54500000000004 - 5.58197349185435e6*eta - 3.5225742421184325e8*eta2 + 1.4667258334378073e9*eta3)/(1. + 66757.12830903867*eta + 5.385164380400193e6*eta2 + 2.5176585751772933e6*eta3);

      eqSpin = (S*(19.416719811164853 - 36.066611959079935*S - 0.8612656616290079*S2 + eta2*(170.97203068800542 - 107.41099349364234*S - 647.8103976942541*S3) + 5.95010003393006*S3 + eta3*(-1365.1499998427248 + 1152.425940764218*S + 415.7134909564443*S2 + 1897.5444343138167*S3 - 866.283566780576*S4) + 4.984750041013893*S4 + eta*(207.69898051583655 - 132.88417400679026*S - 17.671713040498304*S2 + 29.071788188638315*S3 + 37.462217031512786*S4)))/(-1.1492259468169692 + 1.*S);

      uneqSpin = dchi*delta*eta3*(7343.130973149263 - 20486.813161100774*eta + 515.9898508588834*S);

  		break;
		}
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Intermediate_Phase_22_v3: IMRPhenomXIntermediatePhaseVersion is not valid. Recommended flag is 104.\n");
      break;
    }
	}

	return (noSpin + eqSpin + uneqSpin);
}

/* v2mRDv4 = vIM2 - vRD4. See Section VII.B of arXiv:2001.11412.  */
static double IMRPhenomX_Intermediate_Phase_22_v2mRDv4(double eta, double S, double dchi, double delta, int IntPhaseFlag){

	double eta2  = eta*eta;
	double eta3  = eta2*eta;

	double S2    = S*S;

	double noSpin, eqSpin, uneqSpin;

	switch ( IntPhaseFlag )
	{
		case 104: 			/* Canonical, 4 coefficients */
		{
      noSpin = (eta*(-8.244230124407343 - 182.80239160435949*eta + 638.2046409916306*eta2 - 578.878727101827*eta3))/(-0.004863669418916522 - 0.5095088831841608*eta + 1.*eta2);

      eqSpin = (S*(0.1344136125169328 + 0.0751872427147183*S + eta2*(7.158823192173721 + 25.489598292111104*S - 7.982299108059922*S2) + eta*(-5.792368563647471 + 1.0190164430971749*S + 0.29150399620268874*S2) + 0.033627267594199865*S2 + eta3*(17.426503081351534 - 90.69790378048448*S + 20.080325133405847*S2)))/(0.03449545664201546 - 0.027941977370442107*S + (0.005274757661661763 + 0.0455523144123269*eta - 0.3880379911692037*eta2 + 1.*eta3)*S2);

      uneqSpin = 160.2975913661124*dchi*delta*eta2;

  		break;
		}
		case 105:				/* Canonical, 5 coefficients */
		{
      noSpin = (eta*(0.9951733419499662 + 101.21991715215253*eta + 632.4731389009143*eta2))/(0.00016803066316882238 + 0.11412314719189287*eta + 1.8413983770369362*eta2 + 1.*eta3);

      eqSpin = (S*(18.694178521101332 + 16.89845522539974*S + 4941.31613710257*eta2*S + eta*(-697.6773920613674 - 147.53381808989846*S2) + 0.3612417066833153*S2 + eta3*(3531.552143264721 - 14302.70838220423*S + 178.85850322465944*S2)))/(2.965640445745779 - 2.7706595614504725*S + 1.*S2);

      uneqSpin = dchi*delta*eta2*(356.74395864902294 + 1693.326644293169*eta2*S);

  		break;
		}
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Intermediate_Phase_22_d23: IMRPhenomXIntermediatePhaseVersion is not valid. Recommended flag is 104.\n");
      break;
    }
	}

	return (noSpin + eqSpin + uneqSpin);
}

/* v3mRDv4 = vIM3 - vRD4. See Section VII.B of arXiv:2001.11412.  */
static double IMRPhenomX_Intermediate_Phase_22_v3mRDv4(double eta, double S, double dchi, double delta, int IntPhaseFlag){

	double eta2  = eta*eta;
	double eta3  = eta2*eta;
	double eta4  = eta3*eta;
	double eta6  = eta3*eta3;

	double S2    = S*S;

	double noSpin, eqSpin, uneqSpin;

	switch ( IntPhaseFlag )
	{
		case 104: 			/* Canonical, 4 coefficients */
		{
      noSpin = (0.3145740304678042 + 299.79825045000655*eta - 605.8886581267144*eta2 - 1302.0954503758007*eta3)/(1. + 2.3903703875353255*eta - 19.03836730923657*eta2);

      eqSpin = (S*(1.150925084789774 - 0.3072381261655531*S + eta4*(12160.22339193134 - 1459.725263347619*S - 9581.694749116636*S2) + eta2*(1240.3900459406875 - 289.48392062629966*S - 1218.1320231846412*S2) - 1.6191217310064605*S2 + eta*(-41.38762957457647 + 60.47477582076457*S2) + eta3*(-7540.259952257055 + 1379.3429194626635*S + 6372.99271204178*S2)))/(-1.4325421225106187 + 1.*S);

      uneqSpin = dchi*delta*eta3*(-444.797248674011 + 1448.47758082697*eta + 152.49290092435044*S);

  		break;
		}
		case 105:				/* Canonical, 5 coefficients */
		{
      noSpin = (eta*(-5.126358906504587 - 227.46830225846668*eta + 688.3609087244353*eta2 - 751.4184178636324*eta3))/(-0.004551938711031158 - 0.7811680872741462*eta + 1.*eta2);

      eqSpin = (S*(0.1549280856660919 - 0.9539250460041732*S - 539.4071941841604*eta2*S + eta*(73.79645135116367 - 8.13494176717772*S2) - 2.84311102369862*S2 + eta3*(-936.3740515136005 + 1862.9097047992134*S + 224.77581754671272*S2)))/(-1.5308507364054487 + 1.*S);

      uneqSpin = 2993.3598520496153*dchi*delta*eta6;

  		break;
		}
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Intermediate_Phase_22_d23: IMRPhenomXIntermediatePhaseVersion is not valid. Recommended flag is 104.\n");
      break;
    }
	}

	return (noSpin + eqSpin + uneqSpin);
}

/* d43 = v4 - v3. See Section VII.B of arXiv:2001.11412.  */
static double IMRPhenomX_Intermediate_Phase_22_d43(double eta, double S, double dchi, double delta, int IntPhaseFlag){

	double eta2  = eta*eta;
	double eta3  = eta2*eta;
	double eta4  = eta3*eta;

	double S2    = S*S;

	double noSpin, eqSpin, uneqSpin;

	switch ( IntPhaseFlag )
	{
		case 104:
		{
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Intermediate_Phase_22_d43: Point d43 should not be called when using IMRPhenomXIntermediatePhaseVersion == 104.\n");
      break;
		}
		case 105:				/* Canonical, 5 coefficients */
		{
      noSpin = (0.4248820426833804 - 906.746595921514*eta - 282820.39946006844*eta2 - 967049.2793750163*eta3 + 670077.5414916876*eta4)/(1. + 1670.9440812294847*eta + 19783.077247023448*eta2);

      eqSpin = (S*(0.22814271667259703 + 1.1366593671801855*S + eta3*(3499.432393555856 - 877.8811492839261*S - 4974.189172654984*S2) + eta*(12.840649528989287 - 61.17248283184154*S2) + 0.4818323187946999*S2 + eta2*(-711.8532052499075 + 269.9234918621958*S + 941.6974723887743*S2) + eta4*(-4939.642457025497 - 227.7672020783411*S + 8745.201037897836*S2)))/(-1.2442293719740283 + 1.*S);

      uneqSpin = dchi*delta*(-514.8494071830514 + 1493.3851099678195*eta)*eta3;

  		break;
		}
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Intermediate_Phase_22_d43: IMRPhenomXIntermediatePhaseVersion is not valid. Recommended flag is 104.\n");
      break;
    }
	}

	return (noSpin + eqSpin + uneqSpin);
}


/*	
	See section VII. B of arXiv:2001.11412. 

	Phase derivative ansatz for intermediate region, Eq. 7.6 of arXiv:2001.11412.

	This is the canonical intermediate ansatz:

	a_0 + a_1 ft^(-1) + a_2 ft^(-2) + a_3 ft^(-3) + a4 ft^(-4) + (4 * a_RD) / ( (2 * f_fdamp)^2 + (f - f_ring)^2 )

	ft = (f / f_ring)
*/
static double IMRPhenomX_Intermediate_Phase_22_Ansatz(double ff, IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXWaveformStruct *pWF, IMRPhenomXPhaseCoefficients *pPhase)
{
  double invff1 = powers_of_f->m_one;
  double invff2 = powers_of_f->m_two;
  double invff3 = powers_of_f->m_three;
  double invff4 = powers_of_f->m_four;

  double fda    = pWF->fDAMP;
  double frd    = pWF->fRING;

  double LorentzianTerm;
  double phaseOut = 0;

  int IntPhaseVersion = pWF->IMRPhenomXIntermediatePhaseVersion;

  switch ( IntPhaseVersion )
  {
    case 104: 			/* Canonical, 4 coefficients */
    {

      /* This is the Lorentzian term where cL = - a_{RD} dphase0 */
      LorentzianTerm = (4.0 * pPhase->cL) / ( (4.0*fda*fda) + (ff - frd)*(ff - frd) );

      /* Return a polynomial embedded in the background from the merger */
      phaseOut       = pPhase->b0 + pPhase->b1*invff1 + pPhase->b2*invff2 + pPhase->b4*invff4 + LorentzianTerm;

      break;
    }
    case 105:				/* Canonical, 5 coefficients */
    {

      /* This is the Lorentzian term where cL = - a_{RD} dphase0 */
      LorentzianTerm = (4.0 * pPhase->cL) / ( (4.0*fda*fda) + (ff - frd)*(ff - frd) );

      /* Return a polynomial embedded in the background from the merger */
      phaseOut       = (pPhase->b0) + (pPhase->b1)*invff1 + (pPhase->b2)*invff2 + (pPhase->b3)*invff3 + (pPhase->b4)*invff4 + LorentzianTerm;

      break;
    }
    default:
    {
      XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Intermediate_Phase_22_Ansatz: IMRPhenomXIntermediatePhaseVersion is not valid. Recommended flag is 104.\n");
    }
  }

  return phaseOut;
}

/*	
	See section VII. B of arXiv:2001.11412. 

	Integrated phase ansatz for intermediate region, Eq. 7.6 of arXiv:2001.11412.

    Effective spin parameterization used = StotR
*/
static double IMRPhenomX_Intermediate_Phase_22_AnsatzInt(double f, IMRPhenomX_UsefulPowers *powers_of_f, IMRPhenomXWaveformStruct *pWF, IMRPhenomXPhaseCoefficients *pPhase)
{
  double invff1 = powers_of_f->m_one;
  double invff2 = powers_of_f->m_two;
  double invff3 = powers_of_f->m_three;
  double logfv   = powers_of_f->log;

  double frd = pWF->fRING;
  double fda = pWF->fDAMP;

  double b0  = pPhase->b0;
  double b1  = pPhase->b1;
  double b2  = pPhase->b2;
  double b3  = pPhase->b3;
  double b4  = pPhase->b4;
  double cL  = pPhase->cL;

  double phaseOut;

  int IntPhaseVersion = pWF->IMRPhenomXIntermediatePhaseVersion;

  switch ( IntPhaseVersion )
  {
    case 104: 			/* Canonical, 4 coefficients */
    {
      //phaseOut = pPhase->b0*f + pPhase->b1*logfv - pPhase->b2*invff1 - pPhase->b4*invff3/3.0 + ( 2.0*pPhase->cL*atan( (f - frd) / (2.0 * fda) ) ) / fda ;
      phaseOut = b0*f + b1*logfv - b2*invff1 - (b4*invff3/3.0) + ( 2.0*cL*atan( (f - frd) / (2.0 * fda) ) ) / fda ;
      break;
    }
    case 105:				/* Canonical, 5 coefficients */
    {
      phaseOut = b0*f + b1*logfv - b2*invff1 - b3*invff2/2.0 - (b4*invff3/3.0) + ( 2.0 * cL * atan( (f - frd) / (2.0 * fda) ) ) / fda ;
      break;
    }
    default:
    {
		XLAL_ERROR_REAL8(XLAL_EINVAL, "Error in IMRPhenomX_Intermediate_Phase_22_AnsatzInt: IMRPhenomXIntermediatePhaseVersion is not valid. Recommended flag is 104.\n");
      break;
    }
  }

  return phaseOut;

}
