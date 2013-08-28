/*
*  Copyright (C) 2007 Reinhard Prix, Teviet Creighton
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
 * \date 2003
 * \file
 * \ingroup TwoDMesh_h
 * \brief Some metric computation routines suitable for use in LALCreateTwoDMesh()
 *
 * ### Description ###
 *
 * This module contains metric computation routines suitable for passing
 * into LALCreateTwoDMesh() via the \a params->getMetric
 * function parameter.
 *
 * The routine LALGetNearestMetric() takes as its parameter a
 * \a (void *) pointer to a \a REAL4VectorSequence structure,
 * containing a list of 5-dimensional vectors
 * \f$(x,y,g_{xx},g_{yy},g_{xy})\f$ giving the values of the specified metric
 * components evaluated at the specified position.  The function returns
 * in \a metric the components at the location in \a params
 * nearest to the requested \a position, where ``nearness'' to a
 * point refers to proper distance computed from the metric at that
 * point.
 *
 * The routine LALInterpolateMetricGrid() takes as its parameter a
 * \a (void *) pointer to a \a REAL4Grid structure of physical
 * dimension 2 and data dimension 3, where the first two dimensions refer
 * to the \f$(x,y)\f$ grid points, and the third dimension must be of length
 * 3, storing the three metric components \f$g_{xx}\f$, \f$g_{yy}\f$, and
 * \f$g_{xy}\f$ in that order.  The routine locates the four grid points
 * surrounding the specified \a position, and interpolates the metric
 * values using a simple bilinear fit \f$f(x,y)=A+Bx+Cy+Dxy\f$, where \f$A\f$,
 * \f$B\f$, \f$C\f$, and \f$D\f$ are specified by the four surrounding points.  If
 * the point \f$(x,y)\f$ is outside the region, then the metric components
 * are interpolated from the two nearest edge points if it is nearest to
 * an edge, or the nearest corner point if it is nearest to a corner.
 * (``Nearness'' in this case refers to \em coordinate distances.)
 */

#include <lal/LALStdlib.h>
#include <lal/Grid.h>
#include <lal/TwoDMesh.h>

void
LALGetNearestMetric( LALStatus *stat, REAL4 metric[3], REAL4 position[2], void *params )
{
  UINT4 length;    /* number of metric mesh points */
  UINT4 iMin;      /* index of point with minimum proper distance */
  UINT4 i, j;      /* i = index of current mesh point, j = 5*i  */
  REAL4 dMin;      /* minimum proper distance */
  REAL4 dx, dy, d; /* offsets and distances from current mesh point */
  REAL4VectorSequence *p; /* params cast to a REAL4VectorSequence */
  REAL4 *data;            /* pointer to p->data */

  INITSTATUS(stat);

  p = (REAL4VectorSequence *)params;
  /* This function may be called a lot.  Do error checking only in
     debug mode. */
#ifndef NDEBUG
  if ( lalDebugLevel ) {
    ASSERT( p, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
    ASSERT( p->data, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
    ASSERT( p->length > 0, stat, TWODMESHH_EDIM, TWODMESHH_MSGEDIM );
    ASSERT( p->vectorLength == 5, stat, TWODMESHH_EDIM, TWODMESHH_MSGEDIM );
  }
#endif
  data = p->data;
  length = p->length;

  /* Compute minimum proper distance. */
  iMin = 0;
  dx = position[0] - data[0];
  dy = position[1] - data[1];
  dMin = data[2]*dx*dx + data[3]*dy*dy + 2.0*data[4]*dx*dy;
  for ( i = 0, j = 0; i < length; i++, j+=5 ) {
    dx = position[0] - data[j];
    dy = position[1] - data[j+1];
    if ( ( d = data[j+2]*dx*dx + data[j+3]*dy*dy
	   + 2.0*data[j+4]*dx*dy ) < dMin ) {
      dMin = d;
      iMin = i;
    }
  }

  /* Return appropriate metric components, and exit. */
  j = 5*iMin;
  metric[0] = data[j+2];
  metric[1] = data[j+3];
  metric[2] = data[j+4];
  RETURN( stat );
}


void
LALInterpolateMetricGrid( LALStatus *stat, REAL4 metric[3], REAL4 position[2], void *params )
{
  REAL4Grid *p;         /* params cast to a REAL4Grid */
  REAL8 *p0;            /* pointer to p->offset->data */
  REAL8 *dp;            /* pointer to p->interval->data */
  REAL4 *data;          /* pointer to p->data->data */
  UINT4 iMax, jMax;     /* dimensions of (x,y) grid */
  REAL8 i, j;           /* floating-point indecies over (x,y) grid */
  INT4 iInt, jInt;      /* integral parts of i, j */
  REAL4 iFrac, jFrac;   /* fractional parts of i, j */
  UINT4 k1, k2, k3, k4; /* indecies of four nearby points in grid */

  INITSTATUS(stat);

  p = (REAL4Grid *)params;
   /* This function may be called a lot.  Do error checking only in
     debug mode. */
#ifndef NDEBUG
  if ( lalDebugLevel ) {
    ASSERT( p, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
    ASSERT( p->offset, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
    ASSERT( p->offset->data, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
    ASSERT( p->offset->length == 2, stat, TWODMESHH_EDIM,
	    TWODMESHH_MSGEDIM );
    ASSERT( p->interval, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
    ASSERT( p->interval->data, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
    ASSERT( p->interval->length == 2, stat, TWODMESHH_EDIM,
	    TWODMESHH_MSGEDIM );
    ASSERT( p->interval->data[0] != 0.0, stat, TWODMESHH_EINT,
	    TWODMESHH_MSGEINT );
    ASSERT( p->interval->data[1] != 0.0, stat, TWODMESHH_EINT,
	    TWODMESHH_MSGEINT );
    ASSERT( p->data, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
    ASSERT( p->data->dimLength, stat, TWODMESHH_ENUL,
	    TWODMESHH_MSGENUL );
    ASSERT( p->data->dimLength->data, stat, TWODMESHH_ENUL,
	    TWODMESHH_MSGENUL );
    ASSERT( p->data->dimLength->length == 3, stat, TWODMESHH_EDIM,
	    TWODMESHH_MSGEDIM );
    ASSERT( p->data->dimLength->data[0] > 0, stat, TWODMESHH_EDIM,
	    TWODMESHH_MSGEDIM );
    ASSERT( p->data->dimLength->data[1] > 0, stat, TWODMESHH_EDIM,
	    TWODMESHH_MSGEDIM );
    ASSERT( p->data->dimLength->data[2] == 3, stat, TWODMESHH_EDIM,
	    TWODMESHH_MSGEDIM );
    ASSERT( p->data->data, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
  }
#endif
  p0 = p->offset->data;
  dp = p->interval->data;
  data = p->data->data;
  iMax = p->data->dimLength->data[0];
  jMax = p->data->dimLength->data[1];

  /* Compute interpolated index of point. */
  i = ( position[0] - p0[0] )/dp[0];
  j = ( position[1] - p0[1] )/dp[1];
  iInt = (INT4)( i );
  jInt = (INT4)( j );
  iFrac = i - iInt;
  jFrac = j - jInt;

  /* Find the four nearest grid points. */
  if ( i <= 0.0 ) {
    if ( j <= 0.0 )
      k1 = k2 = k3 = k4 = 0;
    else if ( j >= jMax - 1 )
      k1 = k2 = k3 = k4 = 3*( jMax - 1 );
    else {
      k1 = k2 = 3*jInt;
      k3 = k4 = k1 + 3;
    }
  } else if ( i >= iMax - 1 ) {
    if ( j <= 0.0 )
      k1 = k2 = k3 = k4 = 3*jMax*( iMax - 1 );
    else if ( j >= jMax - 1 )
      k1 = k2 = k3 = k4 = 3*( jMax*iMax - 1 );
    else {
      k1 = k2 = 3*( jMax*( iMax - 1 ) + jInt );
      k3 = k4 = k1 + 3;
    }
  } else {
    if ( j <= 0.0 ) {
      k1 = k3 = 3*iInt*jMax;
      k2 = k4 = k1 + 3*jMax;
    } else if ( j >= jMax - 1 ) {
      k1 = k3 = 3*( iInt + 1 )*jMax;
      k2 = k4 = k1 + 3*jMax;
    } else {
      k1 = 3*( iInt*jMax + jInt );
      k2 = k1 + 3*jMax;
      k3 = k1 + 3;
      k4 = k2 + 3;
    }
  }

  /* Interpolate the four grid points, and exit. */
  metric[0] = ( 1.0 - iFrac )*( 1.0 - jFrac )*data[k1]
    + iFrac*( 1.0 - jFrac )*data[k2] + iFrac*jFrac*data[k4]
    + ( 1.0 - iFrac )*jFrac*data[k3];
  metric[1] = ( 1.0 - iFrac )*( 1.0 - jFrac )*data[k1+1]
    + iFrac*jFrac*data[k4+1] + iFrac*( 1.0 - jFrac )*data[k2+1]
    + ( 1.0 - iFrac )*jFrac*data[k3+1];
  metric[2] = ( 1.0 - iFrac )*( 1.0 - jFrac )*data[k1+2]
    + iFrac*jFrac*data[k4+2] + iFrac*( 1.0 - jFrac )*data[k2+2]
    + ( 1.0 - iFrac )*jFrac*data[k3+2];
  RETURN( stat );
}
