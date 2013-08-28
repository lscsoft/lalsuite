/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton
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
 * \file
 * \ingroup ODE_h
 * \author J. D. E. Creighton
 *
 * \brief Tests the routines in \ref ODE_h by integrating Keplerian orbits.
 * The orbits so integrated are output to files containing the integrated
 * and expected orbits.
 *
 */

/** \cond DONT_DOXYGEN */

#include <math.h>
#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/ODE.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#define TESTSTATUS( ps ) \
  if ( (ps)->statusCode ) { REPORTSTATUS( ps ); exit( 1 ); } else ( (void) 0 )


/*
 *
 * Kepler problem:
 *
 *  x = { X, Y, VX, VY }
 *  xdot = { VX, VY, -X/R^3, -Y/R^3 }
 *  R^2 = X^2 + Y^2
 *
 */
static void Kepler(
    LALStatus     UNUSED *s,
    REAL4Vector   *xdot,
    REAL4Vector   *x,
    REAL4ODEIndep UNUSED *p
    )
{
  REAL4 rsq = x->data[0] * x->data[0] + x->data[1] * x->data[1];
  REAL4 rcb = rsq * sqrt( rsq );
  xdot->data[0] = x->data[2];
  xdot->data[1] = x->data[3];
  xdot->data[2] = - x->data[0] / rcb;
  xdot->data[3] = - x->data[1] / rcb;
}


static int CircularOrbit( void )
{
  enum { NUMDIM = 4, NUMSTP = 4 };
  const UINT4 ndim = NUMDIM;
  const REAL4 tend = 100;
  const REAL4 tstp = 0.01;
  static REAL4ODEParams params;
  static REAL4ODEIndep indep;
  static REAL4Vector *x0;
  static REAL4Vector *x;
  static LALStatus status;
  CreateVectorSequenceIn seqin = { NUMSTP, NUMDIM };
  FILE *fp;

  LALCreateVector( &status, &x, ndim );
  TESTSTATUS( &status );
  LALCreateVector( &status, &x0, ndim );
  TESTSTATUS( &status );
  LALCreateVector( &status, &params.xdot, ndim );
  TESTSTATUS( &status );
  LALCreateVectorSequence( &status, &params.dx, &seqin );
  TESTSTATUS( &status );

  params.ode   = Kepler;
  params.indep = &indep;
  params.tstep = tstp;

  x->data[0] = 1;
  x->data[1] = 0;
  x->data[2] = 0;
  x->data[3] = 1;

  fp = fopen( "orbit_1.out", "w" );
  while ( indep.t < tend )
  {
    memcpy( x0->data, x->data, x->length * sizeof( *x->data ) );
    ( *params.ode )( &status, params.xdot, x0, &indep );
    LALSRungeKutta4( &status, x, x0, &params );
    TESTSTATUS( &status );
    fprintf( fp, "%e\t%e\t", x->data[0], x->data[1] );
    fprintf( fp, "%e\t%e\n", cos( indep.t ), sin( indep.t ) );
  }
  fclose( fp );

  LALDestroyVectorSequence( &status, &params.dx );
  TESTSTATUS( &status );
  LALDestroyVector( &status, &params.xdot );
  TESTSTATUS( &status );
  LALDestroyVector( &status, &x0 );
  TESTSTATUS( &status );
  LALDestroyVector( &status, &x );
  TESTSTATUS( &status );

  return 0;
}


static int EccentricOrbit( void )
{
  enum { NUMDIM = 4, NUMSTP = 6 };
  const UINT4 ndim = NUMDIM;
  const REAL4 nper = 3;
  static REAL4ODEParams params;
  static REAL4ODEIndep indep;
  static REAL4Vector *x0;
  static REAL4Vector *x;
  static LALStatus status;
  CreateVectorSequenceIn seqin = { NUMSTP, NUMDIM };
  REAL4 a;
  REAL4 v;
  REAL4 r;
  REAL4 e;
  REAL4 P;
  FILE *fp;

  LALCreateVector( &status, &x, ndim );
  TESTSTATUS( &status );
  LALCreateVector( &status, &x0, ndim );
  TESTSTATUS( &status );
  LALCreateVector( &status, &params.xdot, ndim );
  TESTSTATUS( &status );
  LALCreateVector( &status, &params.xerr, ndim );
  TESTSTATUS( &status );
  LALCreateVectorSequence( &status, &params.dx, &seqin );
  TESTSTATUS( &status );

  params.eps   = 1e-6;
  params.ode   = Kepler;
  params.indep = &indep;
  params.tstep = 0.1;

  x->data[0] = r = 1;
  x->data[1] = 0;
  x->data[2] = 0;
  x->data[3] = v = 1.2;

  a = 1 / ( 2 / r - v * v );
  P = 2 * LAL_PI * a * sqrt( a );
  e = 1 - r / a;

  fp = fopen( "orbit_2.out", "w" );

  while ( indep.t < nper * P )
  {
    REAL4 fac = 1;
    REAL4 rad;
    REAL4 phi;
    REAL4 psi;
    REAL4 del;
    REAL4 M;

    memcpy( x0->data, x->data, x->length * sizeof( *x->data ) );
    ( *params.ode )( &status, params.xdot, x0, &indep );
    LALSRungeKutta5Adapt( &status, x, x0, &params );
    TESTSTATUS( &status );
    fprintf( fp, "%e\t%e\t%e\t", indep.t, x->data[0], x->data[1] );

    psi = M = 2 * LAL_PI * indep.t / P;
    del = psi - e * sin( psi ) - M;
    while ( fabs( del ) > 10 * LAL_REAL4_EPS * ( 1 + fabs( psi ) ) )
    {
      psi += ( del < 0 ? 1 : -1 ) * ( fac *= 0.5 ) * e;
      del  = psi - e * sin( psi ) - M;
    }
    rad = a * ( 1 - e * cos( psi ) );
    phi = 2 * atan( sqrt( ( 1 + e ) / ( 1 - e ) ) * tan( psi / 2 ) );
    fprintf( fp, "%e\t%e\n", rad * cos( phi ), rad * sin( phi ) );
  }

  fclose( fp );

  LALDestroyVectorSequence( &status, &params.dx );
  TESTSTATUS( &status );
  LALDestroyVector( &status, &params.xerr );
  TESTSTATUS( &status );
  LALDestroyVector( &status, &params.xdot );
  TESTSTATUS( &status );
  LALDestroyVector( &status, &x0 );
  TESTSTATUS( &status );
  LALDestroyVector( &status, &x );
  TESTSTATUS( &status );

  return 0;
}


int main( void )
{
  CircularOrbit();
  EccentricOrbit();
  LALCheckMemoryLeaks();
  return 0;
}
/** \endcond */
