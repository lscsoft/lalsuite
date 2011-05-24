/*
*  Copyright (C) 2007 David Churches, B.S. Sathyaprakash
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
\ingroup LALInspiral_h

\brief The function \c LALInspiralPhiofVIntegrand() calculates the quantity \f$v^{3} E^{\prime}(v)/\mathcal{F}(v)\f$.

\heading{Prototypes}

<tt>LALInspiralPhiofVIntegrand()</tt>

\heading{Description}

The function \c LALInspiralPhiofVIntegrand() calculates the quantity \f$v^{3} E^{\prime}(v)/\mathcal{F}(v)\f$.

\heading{Uses}

This function calls \c dEnergy and \c flux functions that are defined in the
\c expnFunc structure  and represent \f$E^{\prime}(v)\f$ and \f$\mathcal{F}(v)\f$, respectively,
and pointed to the appropriate PN functions with a call to <tt>LALInspiralChooseModel().</tt>

\heading{Notes}



*/

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

NRCSID (LALINSPIRALPHIOFVINTEGRANDC, "$Id$");


void
LALInspiralPhiofVIntegrand (
   LALStatus  *status,
   REAL8      *integrand,
   REAL8       v,
   void       *params
   )
{

  PhiofVIntegrandIn *in;

  INITSTATUS (status, "LALInspiralPhiofVIntegrand", LALINSPIRALPHIOFVINTEGRANDC);
  ATTATCHSTATUSPTR(status);

  ASSERT (integrand, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT (params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT (v>0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT (v<1, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  in = (PhiofVIntegrandIn *) params;

  *integrand = pow (v, 3.) * in->dEnergy(v,in->coeffs)/in->flux(v,in->coeffs);

  DETATCHSTATUSPTR(status);
  RETURN(status);


}
