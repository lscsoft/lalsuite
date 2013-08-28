/*
*  Copyright (C) 2007 David Churches, B.S. Sathyaprakash, Drew Keppel
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
 * \author Sathyaprakash, B. S.
 * \file
 * \ingroup LALSimInspiraldEnergyFlux_c
 *
 * \brief NONE
 *
 * This module outputs
 * \f{equation}{
 * tofv = t - t_0 + m \int_{v_0}^{v} \frac{E'(v)}{{\cal F}(v)} \, dv\,.
 * \f}
 * where the constants \f$t,\f$ \f$t_0,\f$ \f$v_0,\f$ and functions in the integrand
 * \f$E'(v)\f$ and \f${\cal F}(v)\f$ are defined in the \c void structure <tt>params.</tt>
 *
 * \heading{Uses}
 * \code
 * XLALDRombergIntegrate()
 * \endcode
 *
 */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALSimInspiral.h>
#include <lal/Integrate.h>

#include "LALSimInspiraldEnergyFlux.c"
#include "LALSimInspiralPNCoefficients.c"

static REAL8
XLALSimInspiralTofVIntegrand(
   REAL8      v,
   void      *params
   )
{

  TofVIntegrandIn *ak = NULL;

#ifndef LAL_NDEBUG
  if ( !params )
    XLAL_ERROR_REAL8( XLAL_EFAULT );

  if ( v <= 0.0 || v >= 1.0 )
    XLAL_ERROR_REAL8( XLAL_EINVAL );
#endif

  ak = (TofVIntegrandIn *) params;

  return ak->dEnergy( v, ak->coeffs ) / ak->flux( v, ak->coeffs );
}

static REAL8
XLALSimInspiralTofV (
   REAL8 v,
   void *params
   )
{
   void *funcParams;
   REAL8 (*funcToIntegrate)(REAL8, void *);
   REAL8 xmin, xmax;
   IntegralType type;
   TofVIntegrandIn in2;
   TofVIn *in1;
   REAL8 answer;
   REAL8 sign;


   if (params == NULL)
      XLAL_ERROR_REAL8(XLAL_EFAULT);
   if (v <= 0.)
      XLAL_ERROR_REAL8(XLAL_EDOM);
   if (v >= 1.)
      XLAL_ERROR_REAL8(XLAL_EDOM);

   sign = 1.0;


   in1 = (TofVIn *) params;

   funcToIntegrate = XLALSimInspiralTofVIntegrand;
   xmin = in1->v0;
   xmax = v;
   type = ClosedInterval;


   in2.dEnergy = in1->dEnergy;
   in2.flux = in1->flux;
   in2.coeffs = in1->coeffs;

   funcParams = (void *) &in2;

   if (v==in1->v0)
   {
     return in1->t - in1->t0;
   }

   if(in1->v0 > v)
   {
      xmin = v;
      xmax = in1->v0;
      sign = -1.0;
   }

   answer = XLALREAL8RombergIntegrate (funcToIntegrate, funcParams, xmin, xmax, type);
   if (XLAL_IS_REAL8_FAIL_NAN(answer))
      XLAL_ERROR_REAL8(XLAL_EFUNC);

   return in1->t - in1->t0 + in1->totalmass*answer*sign;
}

REAL8
XLALSimInspiralTaylorLength(
    REAL8 deltaT,   /**< sampling interval */
    REAL8 m1,       /**< mass of companion 1 */
    REAL8 m2,       /**< mass of companion 2 */
    REAL8 f_min,    /**< start frequency */
    int O           /**< twice post-Newtonian order */
    )
{
    expnCoeffsdEnergyFlux akEF;
    REAL8 m = m1 + m2;
    REAL8 eta = m1 * m2 / (m * m);
    m *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert totalmass from kilograms to seconds */
    REAL8 v0 = cbrt(LAL_PI * m * f_min);

    REAL8 oneby6 = 1./6.;
    REAL8 lso, vlso, vn, tofv;
    TofVIn in1;
    void *in2;

/* Taylor coefficients of dE(v)/dv. (NOTE v and NOT x) */
    akEF.dETaN = 2.0 * XLALSimInspiralPNEnergy_0PNCoeff(eta);
    akEF.dETa1 = 2.0 * XLALSimInspiralPNEnergy_2PNCoeff(eta);
    akEF.dETa2 = 3.0 * XLALSimInspiralPNEnergy_4PNCoeff(eta);
    akEF.dETa3 = 4.0 * XLALSimInspiralPNEnergy_6PNCoeff(eta);

/* Taylor coefficients of flux. */
    akEF.FTaN = XLALSimInspiralPNFlux_0PNCoeff(eta);
    akEF.FTa2 = XLALSimInspiralPNFlux_2PNCoeff(eta);
    akEF.FTa3 = XLALSimInspiralPNFlux_3PNCoeff(eta);
    akEF.FTa4 = XLALSimInspiralPNFlux_4PNCoeff(eta);
    akEF.FTa5 = XLALSimInspiralPNFlux_5PNCoeff(eta);
    akEF.FTa6 = XLALSimInspiralPNFlux_6PNCoeff(eta);
    akEF.FTl6 = XLALSimInspiralPNFlux_6PNLogCoeff(eta);
    akEF.FTa7 = XLALSimInspiralPNFlux_7PNCoeff(eta);
    akEF.FTa8 = - 117.5043907226773;
    akEF.FTl8 =   52.74308390022676;

    lso = sqrt(oneby6);
    vlso = 0; //- set but not used

/* Location of the 0PN and 1PN T- and P-approximant last stable orbit: */
    akEF.vlsoT0 = lso;
    akEF.vlsoP0 = lso;
    akEF.vlsoP2 = lso;
/*
  vlsoT2 =  6./(9.+eta);
  This correct value makes vlso too large for vlsoT2 hence use 1/sqrt(6)
*/
    akEF.vlsoT2 = lso;

    switch (O)
    {
        case 0:
            vn = akEF.vlso = vlso = akEF.vlsoT0;
            in1.dEnergy = XLALSimInspiraldEt0;
            in1.flux = XLALSimInspiralFt0;
            break;
        case 1:
            XLALPrintError("XLAL Error - %s: PN approximant not supported for requested PN order\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
        case 2:
            vn = akEF.vlso = vlso = akEF.vlsoT2;
            in1.dEnergy = XLALSimInspiraldEt2;
            in1.flux = XLALSimInspiralFt2;
            break;
        case 3:
            vn = akEF.vlso = vlso = akEF.vlsoT2;
            in1.dEnergy = XLALSimInspiraldEt2;
            in1.flux = XLALSimInspiralFt3;
            break;
        case 4:
/*
   The value vlsoT4 is too large and doesn't work sometimes;
   so we use vlsoT2.
*/
            vn = akEF.vlso = vlso = akEF.vlsoT2;
            in1.dEnergy = XLALSimInspiraldEt4;
            in1.flux = XLALSimInspiralFt4;
            break;
        case 5:
/*
   The value vlsoT4 is too large and doesn't work with 2.5 PN
   Taylor approximant; so we use vlsoT2.
*/
            vn = akEF.vlso = vlso = akEF.vlsoT2;
            in1.dEnergy = XLALSimInspiraldEt4;
            in1.flux = XLALSimInspiralFt5;
            break;
            case 6:
/*
   vlsoT6 is as yet undetermined and vlsoT4 is too large in
   certain cases (TaylorT2 crashes for (1.4,10)); using vlsoT2;
*/
            vn = akEF.vlso = vlso = akEF.vlsoT2;
            in1.dEnergy = XLALSimInspiraldEt6;
            in1.flux = XLALSimInspiralFt6;
            break;
        case -1: // Use the highest PN order available, move if higher terms added
        case 7:
            vn = akEF.vlso = vlso = akEF.vlsoT2;
            in1.dEnergy = XLALSimInspiraldEt6;
            in1.flux = XLALSimInspiralFt7;
            break;
        case 8:
           XLALPrintError("XLAL Error - %s: PN approximant not supported for requested PN order\n", __func__);
           XLAL_ERROR(XLAL_EINVAL);
           break;
        default:
            XLALPrintError("XLAL Error - %s: Unknown PN order in switch\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
    }

    vn = cbrt(LAL_PI * m / (2. * deltaT));
    vn = (vn < vlso) ? vn :  vlso;

    in1.t=0.0;
    in1.v0=v0;
    in1.t0=0.;
    in1.vlso=akEF.vlso;
    in1.totalmass = m;
    in1.coeffs = &akEF;

    in2 = (void *) &in1;

    tofv = XLALSimInspiralTofV(vn, in2);
    if (XLAL_IS_REAL8_FAIL_NAN(tofv))
        XLAL_ERROR(XLAL_EFUNC);

    return -tofv - deltaT;
}
