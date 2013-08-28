/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, B.S. Sathyaprakash, Drew Keppel
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
 * \brief Given \f$\tau_0\f$ and \f$\tau_4\f$ solve for the mass ratio \f$\eta.\f$
 *
 * ### Description ###
 *
 * Given \f$\tau_0\f$ and \f$\tau_4\f$ one can determine \f$\eta\f$ by solving
 * \f{equation}{
 * -\eta^{4/5} \tau_4 + A_4 \left ( \frac {\tau_0}{A_0} \right )^{1/5}
 * \left (1 + B_4\eta + C_4 \eta^2 \right )  = 0,
 * \f}
 * where \f$A_0 = 5/[256 (\pi f_{s} )^{8/3}],\f$
 * \f$A_4 = 5 \times 3058673/ [128 \times 1016064  (\pi f_s)^{4/3}],\f$
 * \f$B_4 = 5429 \times 1016064 /(1008 \times 3058673),\f$ and \f$C_4 = 617 \times
 * 1016064/(144 \times 3058673).\f$
 * This function returns the LHS of the above
 * equation in \c x for a given \c eta.
 *
 * ### Algorithm ###
 *
 * None.
 *
 * ### Uses ###
 *
 * None.
 *
 * ### Notes ###
 *
 * The <tt>void pointer</tt> <tt>*p</tt> should point to a \c struct
 * of type ::EtaTau04In
 * \code
 * {
 * void *p;
 * EtaTau04In q;
 * ...
 * p = (void *) &q;
 * }
 * \endcode
 */

#include <lal/LALInspiral.h>

void
LALEtaTau04(
   LALStatus *status,
   REAL8     *x,
   REAL8     eta,
   void      *p
   )
{
   XLALPrintDeprecationWarning("LALEtaTau04", "XLALEtaTau04");

   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);
   ASSERT (x,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

   *x = XLALEtaTau04(eta, p);
   if (XLAL_IS_REAL8_FAIL_NAN(*x))
      ABORTXLAL(status);

   DETATCHSTATUSPTR(status);
   RETURN(status);
}

REAL8
XLALEtaTau04(
   REAL8     eta,
   void      *p
   )
{
   REAL8 x;
   EtaTau04In *q;

   if (p == NULL)
      XLAL_ERROR_REAL8(XLAL_EFAULT);
   if (eta <= 0)
      XLAL_ERROR_REAL8(XLAL_EDOM);

   q = (EtaTau04In *) p;
   x = -q->t4 + q->A4/pow(eta,0.8) * (1. + q->B4*eta + q->C4*eta*eta);

   return x;
}
