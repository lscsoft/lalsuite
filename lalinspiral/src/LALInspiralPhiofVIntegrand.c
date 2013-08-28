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
 * \brief The function \c XLALInspiralPhiofVIntegrand() calculates the quantity \f$v^{3} E^{\prime}(v)/\mathcal{F}(v)\f$.
 *
 * ### Prototypes ###
 *
 * <tt>LALInspiralPhiofVIntegrand()</tt>
 *
 * ### Description ###
 *
 * The function \c XLALInspiralPhiofVIntegrand() calculates the quantity \f$v^{3} E^{\prime}(v)/\mathcal{F}(v)\f$.
 *
 * ### Uses ###
 *
 * This function calls \c dEnergy and \c flux functions that are defined in the
 * \c expnFunc structure  and represent \f$E^{\prime}(v)\f$ and \f$\mathcal{F}(v)\f$, respectively,
 * and pointed to the appropriate PN functions with a call to <tt>XLALInspiralChooseModel().</tt>
 *
 * ### Notes ###
 *
 */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

void
LALInspiralPhiofVIntegrand (
   LALStatus  *status,
   REAL8      *integrand,
   REAL8       v,
   void       *params
   )
{
  XLALPrintDeprecationWarning("LALInspiralPhiofVIntegrand", "XLALInspiralPhiofVIntegrand");

  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  ASSERT (integrand, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  *integrand = XLALInspiralPhiofVIntegrand(v, params);
  if (XLAL_IS_REAL8_FAIL_NAN(*integrand))
    ABORTXLAL(status);

  DETATCHSTATUSPTR(status);
  RETURN(status);


}

REAL8
XLALInspiralPhiofVIntegrand (
   REAL8       v,
   void       *params
   )
{
  REAL8 integrand;
  PhiofVIntegrandIn *in;

  if (params == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);
  if (v <= 0.)
    XLAL_ERROR_REAL8(XLAL_EDOM);
  if (v >= 1.)
    XLAL_ERROR_REAL8(XLAL_EDOM);

  in = (PhiofVIntegrandIn *) params;

  integrand = pow (v, 3.) * in->dEnergy(v,in->coeffs)/in->flux(v,in->coeffs);

  return integrand;
}
