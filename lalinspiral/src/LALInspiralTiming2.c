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

\brief Module used in solving the timing and phasing functions in quadrature for the
::Approximant #TaylorT2.

\heading{Prototypes}

<tt>LALInspiralTiming2()</tt>

\heading{Description}

Given \f$t\f$ and \f$v\f$ this module computes the quantity
\f{equation}{
\c tofv = t - t_C - t_N(v) \sum t_k v^k,
\f}
where the coefficients \f$t_k\f$ and the Newtonian value \f$t_N\f$ are all defined
in Table.\tableref{table_flux}.

\heading{Algorithm}
None


\heading{Uses}
None

\heading{Notes}
None




*/

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>

NRCSID (LALINSPIRALTIMING2C, "$Id$");


void
LALInspiralTiming2_0PN (
   LALStatus *status,
   REAL8     *toff,
   REAL8      f,
   void      *params
   )
{

  InspiralToffInput *toffIn;
  REAL8 v, v8;

  INITSTATUS (status, "LALInspiralTiming2", LALINSPIRALTIMING2C);
  ATTATCHSTATUSPTR(status);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  v = pow(toffIn->piM * f,oneby3);
  v8 = pow(v,8.);

  *toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8;

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

void
LALInspiralTiming2_2PN (
   LALStatus *status,
   REAL8     *toff,
   REAL8      f,
   void      *params
   )
{

  InspiralToffInput *toffIn;
  REAL8 v, v2, v8;

  INITSTATUS (status, "LALInspiralTiming2", LALINSPIRALTIMING2C);
  ATTATCHSTATUSPTR(status);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  v = pow(toffIn->piM * f,oneby3);
  v2 = v*v;
  v8 = v2*v2*v2*v2;

  *toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1.
        + toffIn->t2 * v2);


  DETATCHSTATUSPTR(status);
  RETURN(status);
}

void
LALInspiralTiming2_3PN (
   LALStatus *status,
   REAL8     *toff,
   REAL8      f,
   void      *params
   )
{

  InspiralToffInput *toffIn;
  REAL8 v, v2, v3, v8;

  INITSTATUS (status, "LALInspiralTiming2", LALINSPIRALTIMING2C);
  ATTATCHSTATUSPTR(status);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  v = pow(toffIn->piM * f,oneby3);
  v2 = v*v;
  v3 = v2*v;
  v8 = v3*v3*v2;

  *toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1.
        + toffIn->t2 * v2
        + toffIn->t3 * v3);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

void
LALInspiralTiming2_4PN (
   LALStatus *status,
   REAL8     *toff,
   REAL8      f,
   void      *params
   )
{

  InspiralToffInput *toffIn;
  REAL8 v, v2, v3, v4, v8;

  INITSTATUS (status, "LALInspiralTiming2", LALINSPIRALTIMING2C);
  ATTATCHSTATUSPTR(status);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  v = pow(toffIn->piM * f,oneby3);
  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v8 = v4*v4;

  *toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1.
        + toffIn->t2 * v2
        + toffIn->t3 * v3
        + toffIn->t4 * v4);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

void
LALInspiralTiming2_5PN (
   LALStatus *status,
   REAL8     *toff,
   REAL8      f,
   void      *params
   )
{

  InspiralToffInput *toffIn;
  REAL8 v, v2, v3, v4, v5, v8;

  INITSTATUS (status, "LALInspiralTiming2", LALINSPIRALTIMING2C);
  ATTATCHSTATUSPTR(status);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  v = pow(toffIn->piM * f,oneby3);
  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v5 = v4*v;
  v8 = v4*v4;

  *toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1.
        + toffIn->t2 * v2
        + toffIn->t3 * v3
        + toffIn->t4 * v4
        + toffIn->t5 * v5);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

void
LALInspiralTiming2_6PN (
   LALStatus *status,
   REAL8     *toff,
   REAL8      f,
   void      *params
   )
{

  InspiralToffInput *toffIn;
  REAL8 v, v2, v3, v4, v5, v6, v8;

  INITSTATUS (status, "LALInspiralTiming2", LALINSPIRALTIMING2C);
  ATTATCHSTATUSPTR(status);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  v = pow(toffIn->piM * f,oneby3);
  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v5 = v4*v;
  v6 = v5*v;
  v8 = v6*v2;

  *toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1.
        + toffIn->t2 * v2
        + toffIn->t3 * v3
        + toffIn->t4 * v4
        + toffIn->t5 * v5
        + (toffIn->t6 + toffIn->tl6 * log(4*v)) * v6);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

void
LALInspiralTiming2_7PN (
   LALStatus *status,
   REAL8     *toff,
   REAL8      f,
   void      *params
   )
{

  InspiralToffInput *toffIn;
  REAL8 v, v2, v3, v4, v5, v6, v7, v8;

  INITSTATUS (status, "LALInspiralTiming2", LALINSPIRALTIMING2C);
  ATTATCHSTATUSPTR(status);

  ASSERT(toff, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(f > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  toffIn = (InspiralToffInput *) params;

  ASSERT(toffIn->t >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  v = pow(toffIn->piM*f, oneby3);
  v2 = v*v;
  v3 = v2*v;
  v4 = v3*v;
  v5 = v4*v;
  v6 = v5*v;
  v7 = v6*v;
  v8 = v7*v;

  *toff = - toffIn->t + toffIn->tc
        + toffIn->tN / v8 * (1.
        + toffIn->t2 * v2
        + toffIn->t3 * v3
        + toffIn->t4 * v4
        + toffIn->t5 * v5
        + (toffIn->t6 + toffIn->tl6 * log(4*v)) * v6
        + toffIn->t7 * v7);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
