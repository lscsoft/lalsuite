/*
*  Copyright (C) 2007 Anand Sengupta, Thomas Cokelaer, Drew Keppel
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
 * \author Cokelaer T.
 * \file
 * \ingroup LALInspiral_h
 *
 * \brief Given an inspiral template structure containing the binary distance
 * and a set of mass parameters, that module provides functions to compute
 * the related amplitude.
 *
 * ### Description ###
 *
 * The inspiral template structure can stored (1) the distance of the binary (2)
 * a set of binary masses such as the two masses or the total mass and eta and (3)
 * an amplitude which is arbitrary fixed to unity when templates are computed.
 *
 * However we might need to have a template with the physical amplitude (for instance
 * to deal with injections). The function \c XLALInspiralRestrictedAmplitude
 * takes an InspiralTemplate structure as input/output to return the restricted
 * Newtonian amplitude by using the following formula.
 *
 * \f{equation}{
 * A = \frac{4c}{d \eta} M^{-5./3.}
 * \f}
 *
 * where \f$d\f$ is in Mpc and \f$M\f$ in solar mass. The result is stored in the signalAmplitude
 * variable of the inspiral template structure.
 *
 * ### Uses ###
 *
 * When appropriate this function calls:\\
 * \c XLALInspiralParameterCalc
 *
 */

#include <lal/LALInspiral.h>

void
LALInspiralRestrictedAmplitude (LALStatus        *status,
				InspiralTemplate *params )
{
  XLAL_PRINT_DEPRECATION_WARNING("XLALInspiralRestrictedAmplitude");

  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  if ( XLALInspiralRestrictedAmplitude(params) == XLAL_FAILURE )
  {
    ABORTXLAL(status);
  }
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

int
XLALInspiralRestrictedAmplitude (InspiralTemplate *params)
{
  if (params == NULL)
    XLAL_ERROR(XLAL_EFAULT);
  if ((INT4)params->massChoice < 0)
    XLAL_ERROR(XLAL_EDOM);
  if ((INT4)params->massChoice > 15)
    XLAL_ERROR(XLAL_EDOM);

  if (params->massChoice != totalMassAndEta)
  {
    if ( XLALInspiralParameterCalc(params) == XLAL_FAILURE )
    {
      XLAL_ERROR(XLAL_EFUNC);
    }
  }

  params->signalAmplitude = 4. * params->totalMass  * params->eta   /  (LAL_PC_SI * 1e6 *params->distance / LAL_MRSUN_SI);

  return XLAL_SUCCESS;
}
