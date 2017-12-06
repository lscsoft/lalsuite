/*
*  Copyright (C) 2007 Teviet Creighton
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
 * \author Creighton, T. D.
 * \file
 * \ingroup TwoDMesh_h
 * \brief Some range computation routines suitable for use in LALCreateTwoDMesh()
 *
 * ### Description ###
 *
 * This module contains range computation routines suitable for passing
 * into LALCreateTwoDMesh() via the <tt>params->getRange</tt>
 * function parameter.
 *
 * The routine LALInterpolateRangePolygon() takes as its parameter
 * a <tt>(void *)</tt> pointer to a \c REAL4VectorSequence structure,
 * containing a list of 2-dimensional vectors \f$(x,y)\f$ giving the
 * locations of points on the boundary of a polygonal region.  The
 * function returns in \c range the points where a vertical line at
 * the specified \f$x\f$-value crosses the edge of the polygon, in ascending
 * order.  If no intersections are found, then both range values are set
 * equal to (one of) the nearest point(s) on the boundary.
 *
 * The routine LALInterpolateRangeGrid() takes as its parameter a
 * <tt>(void *)</tt> pointer to a \c REAL4Grid structure with physical
 * dimension 1 and data dimension 2: for each point \f$x\f$ along the
 * physical dimension, the grid stores a vector of length 2, giving the
 * lower and upper range values \f$y_1(x)\f$ and \f$y_2(x)\f$.  The routine
 * linearly interpolates these two sampled functions to compute the range
 * interval for any specified \f$x\f$.  If the specified \f$x\f$ is outside the
 * grid, then both range values are set equal to the average of the range
 * points at the nearest endpoint of the grid.
 *
 * ### Algorithm ###
 *
 * The LALInterpolateRangePolygon() function is just a stub at
 * present; it returns [0,1] as its range regardless of inputs.
 *
 * ### Uses ###
 *
 * \code
 * lalDebugLevel
 * \endcode
 *
 * ### Notes ###
 *
 */

#include <lal/LALStdlib.h>
#include <lal/Grid.h>
#include <lal/TwoDMesh.h>

void
LALInterpolateRangePolygon( LALStatus *stat, REAL4 range[2], REAL4 x, void *params )
{
  INITSTATUS(stat);

  if ( !params )
    ABORT ( stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );

  /* This function may be called a lot.  Do error checking only in
     debug mode. */
#ifndef LAL_NDEBUG
  REAL4VectorSequence *p; /* params cast to a REAL4VectorSequence */
  p = (REAL4VectorSequence *)params;
  if ( lalDebugLevel ) {
    ASSERT( p, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
    ASSERT( p->data, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
    ASSERT( p->length > 0, stat, TWODMESHH_EDIM, TWODMESHH_MSGEDIM );
    ASSERT( p->vectorLength == 2, stat, TWODMESHH_EDIM, TWODMESHH_MSGEDIM );
  }
#endif

  /* Fill this in later. */
  range[0] = 0.0*x;
  range[1] = 1.0;
  RETURN( stat );
}



void
LALInterpolateRangeGrid( LALStatus *stat, REAL4 range[2], REAL4 x, void *params )
{
  REAL4Grid *p; /* params cast to a REAL4Grid */
  REAL4 *data;  /* pointer to p->data->data */
  UINT4 iMax;   /* dimension of x grid */
  REAL8 i;      /* floating-point index over x grid */
  INT4 iInt;    /* integral part of i */
  REAL4 iFrac;  /* fractional parts of i */

  INITSTATUS(stat);

  p = (REAL4Grid *)params;
  /* This function may be called a lot.  Do error checking only in
     debug mode. */
#ifndef NDEBUG
  if ( lalDebugLevel ) {
    ASSERT( p, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
    ASSERT( p->offset, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
    ASSERT( p->offset->data, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
    ASSERT( p->offset->length == 1, stat, TWODMESHH_EDIM,
	    TWODMESHH_MSGEDIM );
    ASSERT( p->interval, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
    ASSERT( p->interval->data, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
    ASSERT( p->interval->length == 1, stat, TWODMESHH_EDIM,
	    TWODMESHH_MSGEDIM );
    ASSERT( p->interval->data[0] != 0.0, stat, TWODMESHH_EINT,
	    TWODMESHH_MSGEINT );
    ASSERT( p->data, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
    ASSERT( p->data->dimLength, stat, TWODMESHH_ENUL,
	    TWODMESHH_MSGENUL );
    ASSERT( p->data->dimLength->data, stat, TWODMESHH_ENUL,
	    TWODMESHH_MSGENUL );
    ASSERT( p->data->dimLength->length == 2, stat, TWODMESHH_EDIM,
	    TWODMESHH_MSGEDIM );
    ASSERT( p->data->dimLength->data[0] > 0, stat, TWODMESHH_EDIM,
	    TWODMESHH_MSGEDIM );
    ASSERT( p->data->dimLength->data[1] == 2, stat, TWODMESHH_EDIM,
	    TWODMESHH_MSGEDIM );
    ASSERT( p->data->data, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  }
#endif
  data = p->data->data;
  iMax = p->data->dimLength->data[0];

  /* Compute interpolated index of point. */
  i = ( x - p->offset->data[0] )/p->interval->data[0];
  iInt = (INT4)( i );
  iFrac = i - iInt;

  /* Interpolate nearest points. */
  if ( i <= 0.0 )
    range[0] = range[1] = 0.5*( data[0] + data[1] );
  else if ( i >= iMax - 1 )
    range[0] = range[1] = 0.5*( data[2*iMax-2] + data[2*iMax-1] );
  else {
    range[0] = ( 1.0 - iFrac )*data[2*iInt] + iFrac*data[2*iInt+2];
    range[1] = ( 1.0 - iFrac )*data[2*iInt+1] + iFrac*data[2*iInt+3];
  }
  RETURN( stat );
}
