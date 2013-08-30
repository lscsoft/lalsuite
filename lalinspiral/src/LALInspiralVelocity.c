/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton, B.S. Sathyaprakash, Drew Keppel
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
 * \brief The function \c XLALInspiralVelocity() calculates the velocity \f$v\f$ which corresponds to a time \f$t\f$ in
 * the inspiralling binary system.
 *
 * ### Prototypes ###
 *
 * <tt>XLALInspiralVelocity()</tt>
 *
 * ### Description ###
 *
 * The function \c XLALInspiralVelocity() calculates the velocity \f$v\f$ corresponding to a time \f$t\f$
 * in the evolution of an inspiralling binary system.  It does this by iteratively solving
 * \f{equation}{
 * \label{tofv}
 * t(v) =  t_{0} - m \int_{v_{0}}^{v} \frac{E'(v)}{{\cal F}(v)} \, dv \,\,.
 * \f}
 * \c LALInspiralVelocity() calculates \f$v\f$, given \f$t(v)\f$,
 * \f$t_{0}\f$, \f$m\f$, \f$v_{0}\f$, \f$E^{\prime}(v)\f$ and \f$\mathcal{F}(v)\f$.
 *
 * ### Algorithm ###
 *
 *
 * ### Uses ###
 *
 * \c XLALDBisectionFindRoot()
 *
 * ### Notes ###
 *
 */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/FindRoot.h>

void
LALInspiralVelocity(
   LALStatus *status,
   REAL8     *v,
   TofVIn    *ak
   )
{
  XLALPrintDeprecationWarning("LALInspiralVelocity", "XLALInspiralVelocity");

  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  ASSERT (v, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  *v = XLALInspiralVelocity(ak);
  if (XLAL_IS_REAL8_FAIL_NAN(*v))
    ABORTXLAL(status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

REAL8
XLALInspiralVelocity(
   TofVIn    *ak
   )
{
  REAL8 v,xmin,xmax,xacc;
  REAL8 (*rootfunction)(REAL8, void *);
  void *funcParams;

  if (ak == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  rootfunction = XLALInspiralTofV;
  xmax = ak->vlso;
  xmin = ak->v0/2.;
  xacc = 1.0e-8;

  funcParams = (void *) ak;


  if (ak->t==ak->t0)
  {
     return ak->v0;
  }

  v = XLALDBisectionFindRoot(rootfunction, xmin, xmax, xacc, funcParams);
  if (XLAL_IS_REAL8_FAIL_NAN(v))
    XLAL_ERROR_REAL8(XLAL_EFUNC);

  return v;
}
