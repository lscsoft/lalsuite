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
 * \brief This module is used to set the phase of the waveform so that
 * it is equal to the user specified phase \f$\phi_0\f$ when the `velocity' of the
 * system is equal to \f$v.\f$
 *
 * \heading{Prototypes}
 *
 * <tt>XLALInspiralPhasing1()</tt>
 *
 * \heading{Description}
 *
 * The function \c XLALInspiralPhasing1() calculates the phase \f$\phi(v)\f$ using
 * the phasing formula,
 * \anchor phiofv \f{equation}{
 * \phi(v) =  \phi_{0} - 2 \int_{v_{0}}^{v} v^{3} \frac{E'(v)}{{\cal F}(v)} \, dv \,\,.
 * \tag{phiofv}
 * \f}
 * \c XLALInspiralPhasing1() calculates \f$\phi(v)\f$, given \f$\phi_{0}\f$, \f$v_{0}\f$,
 * \f$v\f$, \f$E^{\prime}(v)\f$ and \f$\mathcal{F}(v)\f$.  The user can specify the phase to
 * be of a particular value at an arbitrary point on the waveform when the
 * post-Newtonian evolution variable \f$v\f$ reaches a specific value. Choosing
 * \f$v=v_0,\f$ the initial velocity, means that the initial phase of the wave is \f$\phi_0;\f$
 * Choosing \f$v=v_\textrm{lso}\f$ means that the phase at the last stable orbit is \f$\phi_0\f$ and
 * so on.
 *
 * \heading{Algorithm}
 * Numerical integration.
 *
 * \heading{Uses}
 * XLALDRombergIntegrate()
 * \heading{Notes}
 *
 */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/Integrate.h>

void
LALInspiralPhasing1 (
   LALStatus *status,
   REAL8     *phiofv,
   REAL8     v,
   void      *params
   )
{
   XLALPrintDeprecationWarning("LALInspiralPhasing1", "XLALInspiralPhasing1");

   INITSTATUS(status);
   ATTATCHSTATUSPTR (status);

   ASSERT (phiofv, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

   *phiofv = XLALInspiralPhasing1(v, params);
   if (XLAL_IS_REAL8_FAIL_NAN(*phiofv))
      ABORTXLAL(status);

   DETATCHSTATUSPTR (status);
   RETURN (status);
}

REAL8
XLALInspiralPhasing1 (
   REAL8     v,
   void      *params
   )
{
   void *funcParams;
   REAL8 (*integratedfunction)(REAL8, void *);
   REAL8 xmin, xmax;
   IntegralType integrationtype;
   PhiofVIntegrandIn in2;
   InspiralPhaseIn *in1;
   REAL8 sign;
   REAL8 answer, phiofv;

   if (params == NULL)
      XLAL_ERROR_REAL8(XLAL_EFAULT);
   if (v <= 0.)
      XLAL_ERROR_REAL8(XLAL_EDOM);
   if (v >= 1.)
      XLAL_ERROR_REAL8(XLAL_EDOM);

   sign = 1.0;

   in1 = (InspiralPhaseIn *) params;

   integratedfunction = XLALInspiralPhiofVIntegrand;
   xmin = in1->v0;
   xmax = v;
   integrationtype = ClosedInterval;

   in2.dEnergy = in1->dEnergy;
   in2.flux = in1->flux;
   in2.coeffs = in1->coeffs;

   funcParams = (void *) &in2;

   if (v==in1->v0) {
      return in1->phi0;
   }

   if(in1->v0 > v) {
      xmin = v;
      xmax = in1->v0;
      sign = -1.0;
   }

   answer = XLALREAL8RombergIntegrate (integratedfunction, funcParams, xmin, xmax, integrationtype);
   if (XLAL_IS_REAL8_FAIL_NAN(answer))
      XLAL_ERROR_REAL8(XLAL_EFUNC);

   phiofv = in1->phi0 - 2.0*sign*answer;

   return phiofv;
}
