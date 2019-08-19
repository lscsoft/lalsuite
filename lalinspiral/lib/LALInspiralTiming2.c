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
 * \brief Module used in solving the timing and phasing functions in quadrature for the
 * ::Approximant #TaylorT2.
 *
 * ### Prototypes ###
 *
 * <tt>LALInspiralTiming2()</tt>
 *
 * ### Description ###
 *
 * Given \f$t\f$ and \f$v\f$ this module computes the quantity
 * \f{equation}{
 * tofv = t - t_C - t_N(v) \sum t_k v^k,
 * \f}
 * where the coefficients \f$t_k\f$ and the Newtonian value \f$t_N\f$ are all defined
 * in \ref table_flux "this table".
 *
 * ### Algorithm ###
 *
 * None
 *
 * ### Uses ###
 *
 * None
 *
 * ### Notes ###
 *
 * None
 *
 */

#include <math.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/LALInspiral.h>
#include <lal/XLALError.h>


REAL8
XLALInspiralTiming2_0PN (
   REAL8       f,
   void      *params
   )
{
  InspiralToffInput *toffIn;
  REAL8 v, v8;
  REAL8 toff;

  if (params == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);
  if (f <= 0)
    XLAL_ERROR_REAL8(XLAL_EDOM);

  toffIn = (InspiralToffInput *) params;

  if (toffIn->t < 0)
    XLAL_ERROR_REAL8(XLAL_EDOM);


  v = pow(toffIn->piM * f,(1./3.));
  v8 = pow(v,8.);

  toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8;

  return toff;
}

REAL8
XLALInspiralTiming2_2PN (
   REAL8       f,
   void      *params
   )
{
  InspiralToffInput *toffIn;
  REAL8 v, v2, v8;
  REAL8 toff;

  if (params == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);
  if (f <= 0)
    XLAL_ERROR_REAL8(XLAL_EDOM);

  toffIn = (InspiralToffInput *) params;

  if (toffIn->t < 0)
    XLAL_ERROR_REAL8(XLAL_EDOM);


  v = pow(toffIn->piM * f,(1./3.));
  v2 = v*v;
  v8 = v2*v2*v2*v2;

  toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1.
        + toffIn->t2 * v2);

  return toff;
}

REAL8
XLALInspiralTiming2_3PN (
   REAL8       f,
   void      *params
   )
{
  InspiralToffInput *toffIn;
  REAL8 v, v2, v3, v8;
  REAL8 toff;

  if (params == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);
  if (f <= 0)
    XLAL_ERROR_REAL8(XLAL_EDOM);

  toffIn = (InspiralToffInput *) params;

  if (toffIn->t < 0)
    XLAL_ERROR_REAL8(XLAL_EDOM);


  v = pow(toffIn->piM * f,(1./3.));
  v2 = v*v;
  v3 = v2*v;
  v8 = v3*v3*v2;

  toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1.
        + toffIn->t2 * v2
        + toffIn->t3 * v3);

  return toff;
}

REAL8
XLALInspiralTiming2_4PN (
   REAL8       f,
   void      *params
   )
{
  InspiralToffInput *toffIn;
  REAL8 v, v2, v3, v4, v8;
  REAL8 toff;

  if (params == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);
  if (f <= 0)
    XLAL_ERROR_REAL8(XLAL_EDOM);

  toffIn = (InspiralToffInput *) params;

  if (toffIn->t < 0)
    XLAL_ERROR_REAL8(XLAL_EDOM);


  v = pow(toffIn->piM * f,(1./3.));
  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v8 = v4*v4;

  toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1.
        + toffIn->t2 * v2
        + toffIn->t3 * v3
        + toffIn->t4 * v4);

  return toff;
}

REAL8
XLALInspiralTiming2_5PN (
   REAL8       f,
   void      *params
   )
{
  InspiralToffInput *toffIn;
  REAL8 v, v2, v3, v4, v5, v8;
  REAL8 toff;

  if (params == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);
  if (f <= 0)
    XLAL_ERROR_REAL8(XLAL_EDOM);

  toffIn = (InspiralToffInput *) params;

  if (toffIn->t < 0)
    XLAL_ERROR_REAL8(XLAL_EDOM);


  v = pow(toffIn->piM * f,(1./3.));
  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v5 = v4*v;
  v8 = v4*v4;

  toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1.
        + toffIn->t2 * v2
        + toffIn->t3 * v3
        + toffIn->t4 * v4
        + toffIn->t5 * v5);

  return toff;
}

REAL8
XLALInspiralTiming2_6PN (
   REAL8       f,
   void      *params
   )
{

  InspiralToffInput *toffIn;
  REAL8 v, v2, v3, v4, v5, v6, v8;
  REAL8 toff;

  if (params == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);
  if (f <= 0)
    XLAL_ERROR_REAL8(XLAL_EDOM);

  toffIn = (InspiralToffInput *) params;

  if (toffIn->t < 0)
    XLAL_ERROR_REAL8(XLAL_EDOM);


  v = pow(toffIn->piM * f,(1./3.));
  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v5 = v4*v;
  v6 = v5*v;
  v8 = v6*v2;

  toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1.
        + toffIn->t2 * v2
        + toffIn->t3 * v3
        + toffIn->t4 * v4
        + toffIn->t5 * v5
        + (toffIn->t6 + toffIn->tl6 * log(4*v)) * v6);

  return toff;
}

REAL8
XLALInspiralTiming2_7PN (
   REAL8       f,
   void      *params
   )
{
  InspiralToffInput *toffIn;
  REAL8 v, v2, v3, v4, v5, v6, v7, v8;
  REAL8 toff;

  if (params == NULL)
    XLAL_ERROR_REAL8(XLAL_EFAULT);
  if (f <= 0)
    XLAL_ERROR_REAL8(XLAL_EDOM);

  toffIn = (InspiralToffInput *) params;

  if (toffIn->t < 0)
    XLAL_ERROR_REAL8(XLAL_EDOM);


  v = pow(toffIn->piM*f, (1./3.));
  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v5 = v4*v;
  v6 = v5*v;
  v7 = v6*v;
  v8 = v7*v;

  toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1.
        + toffIn->t2 * v2
        + toffIn->t3 * v3
        + toffIn->t4 * v4
        + toffIn->t5 * v5
        + (toffIn->t6 + toffIn->tl6 * log(4*v)) * v6
        + toffIn->t7 * v7);

  return toff;
}
