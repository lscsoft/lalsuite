/******************************* <lalVerbatim file="TwoDMeshRangesCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{TwoDMeshRanges.c}}
\label{ss:TwoDMeshRanges.c}

Some range computation routines suitable for use in
\verb@LALCreateTwoDMesh()@.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{TwoDMeshRangesCP}
\idx{LALInterpolateRangePolygon()}
\idx{LALInterpolateRangeGrid()}

\subsubsection*{Description}

This module contains range computation routines suitable for passing
into \verb@LALCreateTwoDMesh()@ via the \verb@params->getRange@
function parameter.

The routine \verb@LALInterpolateRangePolygon()@ takes as its parameter
a \verb@(void *)@ pointer to a \verb@REAL4VectorSequence@ structure,
containing a list of 2-dimensional vectors $(x,y)$ giving the
locations of points on the boundary of a polygonal region.  The
function returns in \verb@range@ the points where a vertical line at
the specified $x$-value crosses the edge of the polygon, in ascending
order.  If no intersections are found, then both range values are set
equal to (one of) the nearest point(s) on the boundary.

The routine \verb@LALInterpolateRangeGrid()@ takes as its parameter a
\verb@(void *)@ pointer to a \verb@REAL4Grid@ structure with physical
dimension 1 and data dimension 2: for each point $x$ along the
physical dimension, the grid stores a vector of length 2, giving the
lower and upper range values $y_1(x)$ and $y_2(x)$.  The routine
linearly interpolates these two sampled functions to compute the range
interval for any specified $x$.  If the specified $x$ is outside the
grid, then both range values are set equal to the average of the range
points at the nearest endpoint of the grid.

\subsubsection*{Algorithm}

The \verb@LALInterpolateRangePolygon()@ function is just a stub at
present; it returns [0,1] as its range regardless of inputs.

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{TwoDMeshRangesCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/Grid.h>
#include <lal/TwoDMesh.h>

NRCSID( TWODMESHRANGESC, "$Id$" );

/* <lalVerbatim file="TwoDMeshRangesCP"> */
void
LALInterpolateRangePolygon( LALStatus *stat, REAL4 range[2], REAL4 x, void *params )
{ /* </lalVerbatim> */
  UINT4 length;    /* number of range mesh points */
  REAL4VectorSequence *p; /* params cast to a REAL4VectorSequence */
  REAL4 *data;            /* pointer to p->data */

  INITSTATUS( stat, "LALGetNearestRange", TWODMESHRANGESC );

  p = (REAL4VectorSequence *)params;
  /* This function may be called a lot.  Do error checking only in
     debug mode. */
#ifndef NDEBUG
  if ( lalDebugLevel ) {
    ASSERT( p, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
    ASSERT( p->data, stat, TWODMESHH_ENUL, TWODMESHH_MSGENUL );
    ASSERT( p->length > 0, stat, TWODMESHH_EDIM, TWODMESHH_MSGEDIM );
    ASSERT( p->vectorLength == 2, stat, TWODMESHH_EDIM, TWODMESHH_MSGEDIM );
  }
#endif
  data = p->data;
  length = p->length;

  /* Fill this in later. */
  range[0] = 0.0*x;
  range[1] = 1.0;
  RETURN( stat );
}


/* <lalVerbatim file="TwoDMeshRangesCP"> */
void
LALInterpolateRangeGrid( LALStatus *stat, REAL4 range[2], REAL4 x, void *params )
{ /* </lalVerbatim> */
  REAL4Grid *p; /* params cast to a REAL4Grid */
  REAL4 *data;  /* pointer to p->data->data */
  UINT4 iMax;   /* dimension of x grid */
  REAL8 i;      /* floating-point index over x grid */
  INT4 iInt;    /* integral part of i */
  REAL4 iFrac;  /* fractional parts of i */

  INITSTATUS( stat, "LALGetNearestRange", TWODMESHRANGESC );

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
