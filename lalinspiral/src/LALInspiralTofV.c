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
 * \ingroup LALInspiral_h
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
 * ### Uses ###
 *
 * \code
 * LALDRombergIntegrate()
 * \endcode
 *
 */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/Integrate.h>

void
LALInspiralTofV (
   LALStatus *status,
   REAL8 *tofv,
   REAL8 v,
   void *params
   )
{

   void *funcParams;
   DIntegrateIn intinp;
   TofVIntegrandIn in2;
   TofVIn *in1;
   REAL8 answer;
   REAL8 sign;


   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);

   ASSERT (tofv, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(v > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(v <= 1., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   sign = 1.0;


   in1 = (TofVIn *) params;

   intinp.function = LALInspiralTofVIntegrand;
   intinp.xmin = in1->v0;
   intinp.xmax = v;
   intinp.type = ClosedInterval;


   in2.dEnergy = in1->dEnergy;
   in2.flux = in1->flux;
   in2.coeffs = in1->coeffs;

   funcParams = (void *) &in2;

   if (v==in1->v0)
   {
     *tofv = in1->t - in1->t0;
     DETATCHSTATUSPTR(status);
     RETURN (status);
   }

   if(in1->v0 > v)
   {
      intinp.xmin = v;
      intinp.xmax = in1->v0;
      sign = -1.0;
   }

   LALDRombergIntegrate (status->statusPtr, &answer, &intinp, funcParams);
   CHECKSTATUSPTR(status);

   *tofv = in1->t - in1->t0 + in1->totalmass*answer*sign;

   DETATCHSTATUSPTR(status);
   RETURN (status);
}



REAL8
XLALInspiralTofV (
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

   funcToIntegrate = XLALInspiralTofVIntegrand;
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

