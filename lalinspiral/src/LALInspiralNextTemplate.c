/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton
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

#include <lal/LALInspiralBank.h>

/** \ingroup LALInspiralBank_h
 *  \brief Routine to compute the parameters of the template next to the
 *         current template, but in the positive \f$\tau_{2(3)}\f$ axis.
 * \author Sathyaprakash, B. S.
 *
 * The coarse grid algorithm works by starting at one corner of the
 * parameter space, incrementing along positive \f$\tau_0\f$ direction,
 * with increments determined by the local value of the metric,
 * till the boundary of the parameter space is reached. It then gets back to
 * the starting point and increments along positive \f$\tau_{2(3)}\f$
 * direction, with an increment defined by the metric defined locally;
 * it starts at the first point inside the parameter space but
 * \e consistent with a square lattice. This routine is called each
 * time a translation along the \f$\tau_{2(3)}\f$ direction is required.
 */
void LALInspiralNextTemplate(LALStatus          *status,	/**< LAL status pointer */
                             InspiralBankParams *bankPars,	/**< [out] the parameters of the bank at the next grid point; the point may, indeed, lay outside */
                             InspiralMetric     metric		/**< [in] the value of the metric which would allow computation of the next lattice point (in the \f$\tau_{2(3)}\f$ direction) */
                             )
{

   REAL8 x0tmp, myphi, theta;
   INT4 k;

   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);
   ASSERT (bankPars,  status, LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
   ASSERT (bankPars->dx0 != 0, status, LALINSPIRALBANKH_ESIZE, LALINSPIRALBANKH_MSGESIZE);

   myphi = fabs(atan(bankPars->dx1/bankPars->dx0));
   theta = metric.theta;
   if (theta > myphi) {
      x0tmp = bankPars->x0 + bankPars->dx0 * cos(theta);
      k = floor(x0tmp/bankPars->dx0);
      bankPars->x0 = x0tmp - k * bankPars->dx0;
   } else if (theta > 0 && theta < myphi) {
      x0tmp = bankPars->x0 - bankPars->dx1 * sin(theta);
      k = floor(x0tmp/bankPars->dx0);
      bankPars->x0 = (k+1) * bankPars->dx0 - x0tmp;
   } else if (-theta < myphi) {
      x0tmp = bankPars->x0 - bankPars->dx1 * sin(theta);
      k = floor(x0tmp/bankPars->dx0);
      bankPars->x0 = x0tmp - k * bankPars->dx0;
   } else {
      x0tmp = bankPars->x0 - bankPars->dx1 * cos(theta);
      k = floor(x0tmp/bankPars->dx0);
      bankPars->x0 = (k+1) * bankPars->dx0 - x0tmp;
   }
   DETATCHSTATUSPTR(status);
   RETURN(status);
}
