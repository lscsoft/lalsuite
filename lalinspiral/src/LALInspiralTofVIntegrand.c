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
 * \author Sathyaprakash, B. S.
 * \file
 * \ingroup LALInspiral_h
 *
 * \brief The function \c LALInspiralTofVIntegrand() calculates the quantity \f$E^{\prime}(v)/\mathcal{F}(v)\f$.
 *
 * ### Prototypes ###
 *
 * <tt>LALInspiralTofVIntegrand()</tt>
 *
 * ### Description ###
 *
 * The function \c LALInspiralTofVIntegrand() calculates the quantity \f$E^{\prime}(v)/\mathcal{F}(v)\f$.
 * These are the energy and flux functions which are used in the gravitational wave phasing formula, which is
 * defined as
 *
 * \f{eqnarray}{
 * \label{phasing formula}
 * t(v) & = & t_\textrm{ref} + m \int_v^{v_\textrm{ref}} \,
 * \frac{E'(v)}{{\cal F}(v)} \, dv, \\
 * \phi (v) & = & \phi_\textrm{ref} + 2 \int_v^{v_\textrm{ref}}  v^3 \,
 * \frac{E'(v)}{{\cal F}(v)} \, dv,
 * \f}
 *
 * where \f$v=(\pi m F)^{1/3}\f$ is an invariantly defined velocity, \f$F\f$ is the instantaneous GW frequency, and
 * \f$m\f$ is the total mass of the binary.
 *
 * ### Uses ###
 *
 * This function calls the function which represents \f$E^{\prime}(v)\f$ and \f$\mathcal{F}(v)\f$. The pointer to each
 * of these functions is set by a call to the function \c LALInspiralChooseModel().
 *
 */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

void
LALInspiralTofVIntegrand (
   LALStatus *status,
   REAL8     *integrand,
   REAL8      v,
   void      *params
   )
{

   INITSTATUS(status);

   *integrand = XLALInspiralTofVIntegrand( v, params );

   if ( XLAL_IS_REAL8_FAIL_NAN( *integrand ) ) 
   {
     ABORTXLAL( status );
   }

   RETURN (status);
}

REAL8
XLALInspiralTofVIntegrand(
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
