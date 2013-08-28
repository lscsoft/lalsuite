/*
*  Copyright (C) 2007 Jolien Creighton
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
 * \addtogroup ODE_h
 *
 * ### Description ###
 *
 * The routines \c LALSRungeKutta4(), \c LALSRungeKutta5(), and
 * \c LALSRungeKutta5Adapt() are used to advance an ODE solution from one
 * time step to the next.  \c LALSRungeKutta4() and \c LALSRungeKutta5()
 * advance by the specified time step.  The former uses a 4th order Runge Kutta
 * method while the latter, which uses a 5th order Runge Kutta method, also
 * makes an estimate of the error accumulated in the step.
 * \c LALSRungeKutta5Adapt() uses \c LALSRungeKutta5() to take a step and,
 * if the error in the step is too larger (compared to the fractional error
 * specified by \c eps), then it re-does that step with a finer step size;
 * the step size is modified by the routine.
 *
 * All the routines advance the time after each step.
 *
 * The sequence length of the \c dt field of the parameter structure must
 * be \c 4 for \c LALSRungeKutta4() and \c 6 for
 * \c LALSRungeKutta5() and \c LALSRungeKutta5Adapt().
 *
 * ### Operating Instructions ###
 *
 * The following routine specifies the ODE for the Kepler problem:
 * \f[
 * \frac{d}{dt}\{ x, y, v_x, v_y \} = \{ v_x, v_y, -x/r^3, -y/r^3 \}
 * \f]
 * \code
 * #include <math.h>
 * #include <lal/LALStdlib.h>
 *
 * void Kepler( LALStatus *s, REAL4Vector *xdot, REAL4Vector *x, REAL4ODEIndep *p )
 * {
 *   REAL4 rsq = x->data[0] * x->data[0] + x->data[1] * x->data[1];
 *   REAL4 rcb = rsq * sqrt( rsq );
 *   xdot->data[0] = x->data[2];
 *   xdot->data[1] = x->data[3];
 *   xdot->data[2] = - x->data[0] / rcb;
 *   xdot->data[3] = - x->data[1] / rcb;
 * }
 * \endcode
 *
 * The following programs integrate the Kepler problem.  The first program
 * integrates a circular orbit with fixed step sizes:
 * \code
 * #include <math.h>
 * #include <stdio.h>
 * #include <string.h>
 * #include <lal/LALStdlib.h>
 * #include <lal/AVFactories.h>
 * #include <lal/SeqFactories.h>
 * #include <lal/ODE.h>
 *
 * int main( void )
 * {
 *   const UINT4 ndim = 4;
 *   const UINT4 nstp = 4;
 *   const REAL4 tend = 100;
 *   const REAL4 tstp = 0.01;
 *   static REAL4ODEParams params;
 *   static REAL4ODEIndep indep;
 *   static REAL4Vector *x0;
 *   static REAL4Vector *x;
 *   static LALStatus status;
 *   CreateVectorSequenceIn seqin = { nstp, ndim };
 *
 *   LALCreateVector( &status, &x, ndim );
 *   LALCreateVector( &status, &x0, ndim );
 *   LALCreateVector( &status, &params.xdot, ndim );
 *   LALCreateVectorSequence( &status, &params.dx, &seqin );
 *
 *   params.ode   = Kepler;
 *   params.indep = &indep;
 *   params.tstep = tstp;
 *
 *   x->data[0] = 1;
 *   x->data[1] = 0;
 *   x->data[2] = 0;
 *   x->data[3] = 1;
 *
 *   while ( indep.t < tend )
 *   {
 *     memcpy( x0->data, x->data, x->length * sizeof( *x->data ) );
 *     ( *params.ode )( &status, params.xdot, x0, &indep );
 *     LALSRungeKutta4( &status, x, x0, &params );
 *     printf( "%e\t%e\t", x->data[0], x->data[1] );
 *     printf( "%e\t%e\n", cos( indep.t ), sin( indep.t ) );
 *   }
 *
 *   LALDestroyVectorSequence( &status, &params.dx );
 *   LALDestroyVector( &status, &params.xdot );
 *   LALDestroyVector( &status, &x0 );
 *   LALDestroyVector( &status, &x );
 *
 *   return 0;
 * }
 * \endcode
 *
 * This second program integrates a highly-eccentric bound orbit with adaptive
 * step sizes (and also computes the orbit using Kepler's method for
 * comparison):
 * \code
 * #include <math.h>
 * #include <stdio.h>
 * #include <string.h>
 * #include <lal/LALStdlib.h>
 * #include <lal/LALConstants.h>
 * #include <lal/AVFactories.h>
 * #include <lal/SeqFactories.h>
 * #include <lal/ODE.h>
 *
 * int main( void )
 * {
 *   const UINT4 ndim = 4;
 *   const UINT4 nstp = 6;
 *   const REAL4 nper = 3;
 *   static REAL4ODEParams params;
 *   static REAL4ODEIndep indep;
 *   static REAL4Vector *x0;
 *   static REAL4Vector *x;
 *   static LALStatus status;
 *   CreateVectorSequenceIn seqin = { nstp, ndim };
 *   REAL4 a;
 *   REAL4 v;
 *   REAL4 r;
 *   REAL4 e;
 *   REAL4 P;
 *
 *   LALCreateVector( &status, &x, ndim );
 *   LALCreateVector( &status, &x0, ndim );
 *   LALCreateVector( &status, &params.xdot, ndim );
 *   LALCreateVector( &status, &params.xerr, ndim );
 *   LALCreateVectorSequence( &status, &params.dx, &seqin );
 *
 *   params.eps   = 1e-6;
 *   params.ode   = Kepler;
 *   params.indep = &indep;
 *   params.tstep = 0.1;
 *
 *   x->data[0] = r = 1;
 *   x->data[1] = 0;
 *   x->data[2] = 0;
 *   x->data[3] = v = 1.2;
 *
 *   a = 1 / ( 2 / r - v * v );
 *   P = 2 * LAL_PI * a * sqrt( a );
 *   e = 1 - r / a;
 *
 *   while ( indep.t < nper * P )
 *   {
 *     REAL4 fac = 1;
 *     REAL4 rad;
 *     REAL4 phi;
 *     REAL4 psi;
 *     REAL4 del;
 *     REAL4 M;
 *
 *     memcpy( x0->data, x->data, x->length * sizeof( *x->data ) );
 *     ( *params.ode )( &status, params.xdot, x0, &indep );
 *     LALSRungeKutta5Adapt( &status, x, x0, &params );
 *     printf( "%e\t%e\t%e\t", indep.t, x->data[0], x->data[1] );
 *
 *     psi = M = 2 * LAL_PI * indep.t / P;
 *     del = psi - e * sin( psi ) - M;
 *     while ( fabs( del ) > LAL_REAL4_EPS )
 *     {
 *       psi += ( del < 0 ? 1 : -1 ) * ( fac *= 0.5 ) * e;
 *       del  = psi - e * sin( psi ) - M;
 *     }
 *     rad = a * ( 1 - e * cos( psi ) );
 *     phi = 2 * atan( sqrt( ( 1 + e ) / ( 1 - e ) ) * tan( psi / 2 ) );
 *     printf( "%e\t%e\n", rad * cos( phi ), rad * sin( phi ) );
 *   }
 *
 *   LALDestroyVectorSequence( &status, &params.dx );
 *   LALDestroyVector( &status, &params.xerr );
 *   LALDestroyVector( &status, &params.xdot );
 *   LALDestroyVector( &status, &x0 );
 *   LALDestroyVector( &status, &x );
 *
 *   return 0;
 * }
 * \endcode
 *
 * ### Algorithm ###
 *
 * These routines are based on the methods presented in Numerical Recipes
 * [\ref ptvf1992].
 *
 */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/ODE.h>

#define NSTEP 4

void LALSRungeKutta4(
    LALStatus      *status,
    REAL4Vector    *output,
    REAL4Vector    *input,
    REAL4ODEParams *params
    )
{
  static const REAL4 a[NSTEP] = { 0, 0.5, 0.5, 1.0 };
  static const REAL4 b[NSTEP][NSTEP-1] =
    { { 0, 0, 0 }, { 0.5, 0, 0 }, { 0, 0.5, 0 }, { 0, 0, 1 } };
  static const REAL4 c[NSTEP] = { 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0 };
  UINT4  step;
  UINT4  n;
  UINT4  i;
  REAL4 *x;
  REAL4 *xout;
  REAL4 *xdot;
  REAL4  dt;
  REAL4  t;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  ASSERT( input,  status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( output, status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( params, status, ODEH_ENULL, ODEH_MSGENULL );

  ASSERT( input->data,  status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( output->data, status, ODEH_ENULL, ODEH_MSGENULL );

  ASSERT( output->data != input->data, status, ODEH_ESAME, ODEH_MSGESAME );

  ASSERT( params->xdot, status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( params->dx, status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( params->xdot->data, status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( params->dx->data, status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( params->indep, status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( params->ode, status, ODEH_ENULL, ODEH_MSGENULL );

  n = input->length;

  ASSERT( n, status, ODEH_ESIZE, ODEH_MSGESIZE );
  ASSERT( output->length == n, status, ODEH_ESZMM, ODEH_MSGESZMM );
  ASSERT( params->xdot->length == n, status, ODEH_ESZMM, ODEH_MSGESZMM );
  ASSERT( params->dx->vectorLength == n, status, ODEH_ESZMM, ODEH_MSGESZMM );
  ASSERT( params->dx->length == NSTEP, status, ODEH_ENSTP, ODEH_MSGENSTP );

  x    = input->data;
  xout = output->data;
  xdot = params->xdot->data;
  t    = params->indep->t;
  dt   = params->tstep;

  /* first step: obtain dx from given xdot */
  for ( i = 0; i < n; ++i )
  {
    params->dx->data[i] = dt * xdot[i];
  }

  /* subsequent steps: accumulate dx estimates */
  for ( step = 1; step < NSTEP; ++step )
  {
    REAL4Vector dxvec;
    UINT4 s;
    dxvec.length = n;
    dxvec.data = params->dx->data + n * step;
    params->indep->t = t + a[step] * dt;
    memcpy( xout, x, n * sizeof( *x ) );
    for ( s = 0; s < step; ++s )
    {
      REAL4  bcof = b[step][s];
      REAL4 *dx = params->dx->data + n * s;
      for ( i = 0; i < n; ++i )
      {
        xout[i] += bcof * dx[i];
      }
    }
    ( *params->ode )( status->statusPtr, &dxvec, output, params->indep );
    CHECKSTATUSPTR( status );
    for ( i = 0; i < n; ++i )
    {
      dxvec.data[i] *= dt;
    }
  }

  /* construct output from weighted sum of dx estimates */
  params->indep->t = t + dt;
  memcpy( xout, x, n * sizeof( *x ) );
  for ( step = 0; step < NSTEP; ++step )
  {
    REAL4  ccof = c[step];
    REAL4 *dx = params->dx->data + n * step;
    for ( i = 0; i < n; ++i )
    {
      xout[i] += ccof * dx[i];
    }
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
#undef NSTEP


#define NSTEP 6

void LALSRungeKutta5(
    LALStatus      *status,
    REAL4Vector    *output,
    REAL4Vector    *input,
    REAL4ODEParams *params
    )
{
  static const REAL4 a[NSTEP] = { 0, 1.0/5.0, 3.0/10.0, 3.0/5.0, 1, 7.0/8.0 };
  static const REAL4 b[NSTEP][NSTEP-1] =
    { { 0, 0, 0, 0, 0 },
      { 1.0/5.0, 0, 0, 0, 0 },
      { 3.0/40.0, 9.0/40.0, 0, 0, 0 },
      { 3.0/10.0, -9.0/10.0, 6.0/5.0, 0, 0 },
      { -11.0/54.0, 5.0/2.0, -70.0/27.0, 35.0/27.0 },
      { 1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0,
        253.0/4096.0 } };
  static const REAL4 c[NSTEP] =
    { 37.0/378.0, 0, 250.0/621.0, 125.0/594.0, 0, 512.0/1771.0 };
  static const REAL4 d[NSTEP] =
    { 37.0/378.0 - 2825.0/27648.0, 0, 250.0/621.0 - 18575.0/48384.0,
      125.0/594.0 - 13525.0/55296.0, -277.0/14336.0, 512.0/1771.0 - 0.25 };

  UINT4  step;
  UINT4  n;
  UINT4  i;
  REAL4 *x;
  REAL4 *xout;
  REAL4 *xdot;
  REAL4 *xerr;
  REAL4  dt;
  REAL4  t;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  ASSERT( input,  status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( output, status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( params, status, ODEH_ENULL, ODEH_MSGENULL );

  ASSERT( input->data,  status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( output->data, status, ODEH_ENULL, ODEH_MSGENULL );

  ASSERT( output->data != input->data, status, ODEH_ESAME, ODEH_MSGESAME );

  ASSERT( params->xdot, status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( params->xerr, status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( params->dx, status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( params->xdot->data, status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( params->xerr->data, status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( params->dx->data, status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( params->indep, status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( params->ode, status, ODEH_ENULL, ODEH_MSGENULL );

  n = input->length;

  ASSERT( n, status, ODEH_ESIZE, ODEH_MSGESIZE );
  ASSERT( output->length == n, status, ODEH_ESZMM, ODEH_MSGESZMM );
  ASSERT( params->xdot->length == n, status, ODEH_ESZMM, ODEH_MSGESZMM );
  ASSERT( params->xerr->length == n, status, ODEH_ESZMM, ODEH_MSGESZMM );
  ASSERT( params->dx->vectorLength == n, status, ODEH_ESZMM, ODEH_MSGESZMM );
  ASSERT( params->dx->length == NSTEP, status, ODEH_ENSTP, ODEH_MSGENSTP );

  x    = input->data;
  xout = output->data;
  xdot = params->xdot->data;
  xerr = params->xerr->data;
  t    = params->indep->t;
  dt   = params->tstep;

  /* first step: obtain dx from given xdot */
  for ( i = 0; i < n; ++i )
  {
    params->dx->data[i] = dt * xdot[i];
  }

  /* subsequent steps: accumulate dx estimates */
  for ( step = 1; step < NSTEP; ++step )
  {
    REAL4Vector dxvec;
    UINT4 s;
    dxvec.length = n;
    dxvec.data = params->dx->data + n * step;
    params->indep->t = t + a[step] * dt;
    memcpy( xout, x, n * sizeof( *x ) );
    for ( s = 0; s < step; ++s )
    {
      REAL4  bcof = b[step][s];
      REAL4 *dx = params->dx->data + n * s;
      for ( i = 0; i < n; ++i )
      {
        xout[i] += bcof * dx[i];
      }
    }
    ( *params->ode )( status->statusPtr, &dxvec, output, params->indep );
    CHECKSTATUSPTR( status );
    for ( i = 0; i < n; ++i )
    {
      dxvec.data[i] *= dt;
    }
  }

  /* construct output from weighted sum of dx estimates */
  params->indep->t = t + dt;
  memcpy( xout, x, n * sizeof( *x ) );
  for ( step = 0; step < NSTEP; ++step )
  {
    REAL4  ccof = c[step];
    REAL4 *dx = params->dx->data + n * step;
    for ( i = 0; i < n; ++i )
    {
      xout[i] += ccof * dx[i];
    }
  }

  /* estimate error */
  memset( xerr, 0, n * sizeof( *xerr ) );
  for ( step = 0; step < NSTEP; ++step )
  {
    REAL4 dcof = d[step];
    if ( dcof != 0 )
    {
      REAL4 *dx = params->dx->data + n * step;
      for ( i = 0; i < n; ++i )
      {
        xerr[i] += dcof * dx[i];
      }
    }
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
#undef NSTEP


void LALSRungeKutta5Adapt(
    LALStatus      *status,
    REAL4Vector    *output,
    REAL4Vector    *input,
    REAL4ODEParams *params
    )
{
  REAL4 eps;
  REAL4 t;
  UINT4 n;
  UINT4 i;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  ASSERT( input,  status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( output, status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( params, status, ODEH_ENULL, ODEH_MSGENULL );

  ASSERT( input->data,  status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( output->data, status, ODEH_ENULL, ODEH_MSGENULL );

  ASSERT( output->data != input->data, status, ODEH_ESAME, ODEH_MSGESAME );

  ASSERT( params->xdot, status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( params->xerr, status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( params->dx, status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( params->xdot->data, status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( params->xerr->data, status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( params->dx->data, status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( params->indep, status, ODEH_ENULL, ODEH_MSGENULL );
  ASSERT( params->ode, status, ODEH_ENULL, ODEH_MSGENULL );

  n = input->length;

  ASSERT( n, status, ODEH_ESIZE, ODEH_MSGESIZE );
  ASSERT( output->length == n, status, ODEH_ESZMM, ODEH_MSGESZMM );
  ASSERT( params->xdot->length == n, status, ODEH_ESZMM, ODEH_MSGESZMM );
  ASSERT( params->xerr->length == n, status, ODEH_ESZMM, ODEH_MSGESZMM );
  ASSERT( params->dx->vectorLength == n, status, ODEH_ESZMM, ODEH_MSGESZMM );

  t = params->indep->t;
  eps = params->eps < LAL_REAL4_EPS ? LAL_REAL4_EPS : params->eps;

  while ( 1 )
  {
    REAL4 maxerr = 0;

    params->indep->t = t;
    LALSRungeKutta5( status->statusPtr, output, input, params );
    CHECKSTATUSPTR( status );

    for ( i = 0; i < n; ++i )
    {
      REAL4 scale = fabs( input->data[i] )
        + params->tstep * fabs( params->xdot->data[i] ) + LAL_REAL4_MIN;
      REAL4 err = fabs( params->xerr->data[i] / scale );
      if ( err > maxerr )
      {
        maxerr = err;
      }
    }
    maxerr /= eps;

    if ( maxerr > 1 )
    { /* reduce stepsize and do over */
      REAL4 fac = 0.9 * pow( maxerr, -0.25 );
      params->tstep *= fac > 0.1 ? fac : 0.1;
    }
    else
    { /* increase stepsize for next iteration */
      REAL4 fac = 0.9 * pow( maxerr, -0.2 );
      params->tstep *= fac < 5 ? ( fac > 1 ? fac : 1 ) : 5;
      break;
    }
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
