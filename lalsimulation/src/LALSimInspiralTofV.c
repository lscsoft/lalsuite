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
\author Sathyaprakash, B. S.
\file
\ingroup LALSimInspiraldEnergyFlux_h

\brief NONE

This module outputs
\f{equation}{
\c tofv = t - t_0 + m \int_{v_0}^{v} \frac{E'(v)}{{\cal F}(v)} \, dv\,.
\f}
where the constants \f$t,\f$ \f$t_0,\f$ \f$v_0,\f$ and functions in the integrand
\f$E'(v)\f$ and \f${\cal F}(v)\f$ are defined in the \c void structure <tt>params.</tt>

\heading{Uses}
\code
XLALDRombergIntegrate()
\endcode

*/

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimInspiraldEnergyFlux.h>
#include <lal/Integrate.h>

NRCSID (LALINSPIRALTOFVC, "$Id$");


REAL8
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
    REAL8 EulerC = LAL_GAMMA;
    REAL8 lambda = -1987./3080.;
    REAL8 lso, vlso, vn, tofv;
    TofVIn in1;
    void *in2;

/* Taylor coefficients of E(x) */
    akEF.ETaN = -eta/2.;
    akEF.ETa1 = -(9. + eta)/12.;
    akEF.ETa2 = -(27. - 19*eta + eta*eta/3.)/8.;
    akEF.ETa3 = -675./64. + (209323./4032. - 205.*LAL_PI*LAL_PI/96.
        - 110./9. * lambda)*eta
        - 155./96. * eta*eta - 35./5184. * eta*eta*eta;

/* Taylor coefficients of dE(v)/dv. (NOTE v and NOT x) */
    akEF.dETaN = -eta;
    akEF.dETa1 = 2.*akEF.ETa1;
    akEF.dETa2 = 3.*akEF.ETa2;
    akEF.dETa3 = 4.*akEF.ETa3;

/* Taylor coefficients of flux. */
    akEF.FTaN = 32.*eta*eta/5.;
    akEF.FTa1 = 0.;
    akEF.FTa2 = -1247./336.-35.*eta/12.;
    akEF.FTa3 = 4.*LAL_PI;
    akEF.FTa4 = -44711./9072.+9271.*eta/504.+65.*eta*eta/18.;
    akEF.FTa5 = -(8191./672.+583./24.*eta)*LAL_PI;
    akEF.FTl6 = -1712./105.;
    akEF.FTa6 = 6643739519./69854400. + 16.*LAL_PI*LAL_PI/3. + akEF.FTl6 * log (4.L)
        - 1712./105.*EulerC+ (-134543./7776. + 41.*LAL_PI*LAL_PI/48.) * eta
        - 94403./3024. * eta*eta - 775./324. * eta*eta*eta;
    akEF.FTa7 = LAL_PI * (-16285./504. + 214745./1728. * eta
        + 193385./3024.* eta*eta);
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
            in1.dEnergy = dEt0;
            in1.flux = Ft0;
            break;
        case 1:
            XLALPrintError("XLAL Error - %s: PN approximant not supported for requested PN order\n", __func__);
            XLAL_ERROR(XLAL_EINVAL);
            break;
        case 2:
            vn = akEF.vlso = vlso = akEF.vlsoT2;
            in1.dEnergy = dEt2;
            in1.flux = Ft2;
            break;
        case 3:
            vn = akEF.vlso = vlso = akEF.vlsoT2;
            in1.dEnergy = dEt2;
            in1.flux = Ft3;
            break;
        case 4:
/*
   The value vlsoT4 is too large and doesn't work sometimes;
   so we use vlsoT2.
*/
            vn = akEF.vlso = vlso = akEF.vlsoT2;
            in1.dEnergy = dEt4;
            in1.flux = Ft4;
            break;
        case 5:
/*
   The value vlsoT4 is too large and doesn't work with 2.5 PN
   Taylor approximant; so we use vlsoT2.
*/
            vn = akEF.vlso = vlso = akEF.vlsoT2;
            in1.dEnergy = dEt4;
            in1.flux = Ft5;
            break;
            case 6:
/*
   vlsoT6 is as yet undetermined and vlsoT4 is too large in
   certain cases (TaylorT2 crashes for (1.4,10)); using vlsoT2;
*/
            vn = akEF.vlso = vlso = akEF.vlsoT2;
            in1.dEnergy = dEt6;
            in1.flux = Ft6;
            break;
        case 7:
            vn = akEF.vlso = vlso = akEF.vlsoT2;
            in1.dEnergy = dEt6;
            in1.flux = Ft7;
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
