/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton, B.S. Sathyaprakash
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
 * \brief Module to calculate the RHS of the differential equations
 * in Eq.\eqref{eq_ode2}.
 *
 * \heading{Prototypes}
 *
 * <tt>LALInspiralDerivatives()</tt>:
 * <ul>
 * <li> \c values: Input containing the values of the variables \f$v\f$ and \f$\phi\f$ at the
 * current time.</li>
 * <li> \c dvalues: Output containing the derivatives \f$dv/dt\f$ and \f$d\phi/dt\f$ at the
 * current time.</li>
 * <li> \c params: Input  of type \c InspiralDerivativesIn that must be
 * cast to a <tt>void.</tt>\\
 * </li>
 * </ul>
 *
 * \heading{Description}
 *
 * This module calculates the right-hand sides of
 * the follwoing two coupled first-order differential equations which are
 * solved to obtain the gravitational wave phasing equation,
 * as described in the documentation for the function \c LALInspiralWave1:
 * The equations are
 * \anchor ode2 \f{equation}{
 * \frac{dv}{dt} = - \frac{\mathcal{F}(v)}{m E^{\prime}(v)},\ \ \ \
 * \frac{d \phi(t)}{dt} = \frac{2v^{3}}{m}.
 * \tag{ode2}
 * \f}
 *
 * \heading{Algorithm}
 *
 * \heading{Uses}
 * None.
 *
 * \heading{Notes}
 *
 * <ul>
 * <li> This function has been intentionally made non-LAL compliant in the sense that it
 * has no status structure.  This is because this code
 * outputs the RHS of the differential equations
 * and is called repeatedly by a function that integrates the two differential
 * equations and should therefore not suffer from undue overheads.</li>
 * <li> The input \c params is of type \c InspiralDerivativesIn and must
 * be cast to a void before calling this function. For example,
 * <code>
 * InspiralDerivativesIn in3;
 * void *funcParams;
 * in3.totalmass = totalmass;
 * ...
 * funcParams = (void *) \&in3;
 * </code>
 * </li>
 * </ul>
 *
 */

#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>

void
LALInspiralDerivatives (
   REAL8Vector *values,
   REAL8Vector *dvalues,
   void        *params
   )
 {

  InspiralDerivativesIn *ak;
  REAL8 v;

  ak = (InspiralDerivativesIn *) params;
  v = *(values->data);

  dvalues->data[0] = -ak->flux(v, ak->coeffs)/ (ak->totalmass*ak->dEnergy(v, ak->coeffs));
  dvalues->data[1] = 2.* v*v*v/ak->totalmass;

  return;
}
