/*
*  Copyright (C) 2007 Jolien Creighton, Reinhard Prix, Teviet Creighton
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

/************************************* <lalVerbatim file="TwoDMeshHV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\providecommand{\lessim}{\stackrel{<}{\scriptstyle\sim}}

\section{Header \texttt{TwoDMesh.h}}
\label{s:TwoDMesh.h}

Provides routines to place search meshes for two-dimensional parameter
spaces with varying metric.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/TwoDMesh.h>
\end{verbatim}

\noindent This header covers routines that lay out a mesh of points on
an 2-dimensional parameter space $\{(x,y)\}$, placed such that no
point in the space lies further than some maximum proper distance
$m_\mathrm{thresh}$ from a mesh point.

The intended purpose of these routines is to place a set of ``target''
search points over a parameter space, in order to detect signals with
unknown parameters.  The formalism for defining a proper distance
metric on the parameter space is defined in the header
\verb@FlatMesh.h@.  However, whereas the routines under
\verb@FlatMesh.h@ require the metric $\mathsf{g}_{ab}$ to be constant
over the parameter space, the routines under this header only treat
$\mathsf{g}_{ab}$ as constant over distances $\lessim
m_\mathrm{thresh}$.

\begin{figure}[h]
\begin{center}
\includegraphics{pulsar_tiling}
\caption{\label{fig:tiling} Mesh placement using parallelogram tiling.
(a)~The left and right sides of a tile are required to be vertical;
the top and bottom sides can tilt to maximize the tile area.
(b)~Tiles can be stacked in fixed-width columns, even as the
elliptical contours change.  (c)~Extra overlapping tiles are sometimes
required at the corners of columns.}
\end{center}
\end{figure}
Since the metric is treated as constant over distances $\lessim
m_\mathrm{thresh}$, this distance defines an elliptical contour around
any mesh point.  We define a ``tile'' as a parallelogram inscribed
within the ellipse, with its left and right sides aligned with the $y$
axis.  This is shown in Fig.~\ref{fig:tiling}(a), above.  A ``column''
is a set of tiles of constant horizontal width stacked one on top of
the other, as shown in Fig.~\ref{fig:tiling}(b).  As the metric
changes over space, the vertical height and tilt of the tiles in a
column may change, so long as their width remains fixed; we note that
if the tilt changes, the tiles will overlap slightly to ensure
complete coverage.  Finally, the boundary of the parameter space may
extend outside the ``corners'' of the column, crossing the end of a
tile between its centre and its edge, as shown in
Fig.~\ref{fig:tiling}(c).  These triangular corners can be covered
with one or more extra overlapping tiles of reduced width.

In a parameter space with constant metric, the tile area is maximized
(and the number of covering tiles minimized) when the column width is
$\sqrt{2}$ times smaller than the projected horizontal width of the
ellipses.  When the ellipses vary, it is generally best to determine
the column width from the \emph{narrowest} ellipse in a column, to
avoid singular effects when tile widths approach the ellipse widths
and become infinitesimally high.

For the column-placement algorithm to work effectively, we require
that the parameter space be representable as a range
$y\in[y_1(x),y_2(x)]$ between two single-valued functions defined on a
domain $x\in[x_\mathrm{min},x_\mathrm{max}]$.  If a desired search
region is too complicated to express this way (e.g.\ it has
disconnected regions, or ``branches'' where a vertical line intersects
the boundary more than twice), then one should divide the region up
into subregions with well-behaved boundary functions and tile these
subregions separately.

This header and its associated modules are placed in the \verb@pulsar@
package because they were originally intended for use in searches over
sky position, but they can be used generically for any two-dimensional
parameter space search where the metric is not too poorly behaved.

******************************************************* </lalLaTeX> */

#ifndef _TWODMESH_H
#define _TWODMESH_H

#include <lal/LALStdlib.h>

#ifdef __cplusplus
extern "C" {
#pragma }
#endif

NRCSID(TWODMESHH,"$Id$");

/********************************************************** <lalLaTeX>
\subsection*{Error conditions}
****************************************** </lalLaTeX><lalErrTable> */
#define TWODMESHH_ENUL    1
#define TWODMESHH_EOUT    2
#define TWODMESHH_EMEM    3
#define TWODMESHH_EMETRIC 4
#define TWODMESHH_EWIDTH  5
#define TWODMESHH_EDIM    6
#define TWODMESHH_EINT    7

#define TWODMESHH_MSGENUL    "Unexpected null pointer in arguments"
#define TWODMESHH_MSGEOUT    "Output handle points to a non-null pointer"
#define TWODMESHH_MSGEMEM    "Memory allocation error"
#define TWODMESHH_MSGEMETRIC "Non-positive metric"
#define TWODMESHH_MSGEWIDTH  "Column width too small"
#define TWODMESHH_MSGEDIM    "Incorrect dimensions"
#define TWODMESHH_MSGEINT    "Non-positive interval"
/******************************************** </lalErrTable><lalLaTeX>

\subsection*{Types}

\subsubsection*{Structure \texttt{TwoDMeshNode}}
\idx[Type]{TwoDMeshNode}

\noindent This structure represents a single node in a linked list of
mesh points, specified in the coordinate system used to place it.  The
fields are:

\begin{description}
\item[\texttt{REAL4 x, y}] The coordinates of the mesh point.

\item[\texttt{REAL4 dx}] The half-width of the tile centred on the
mesh point.

\item[\texttt{REAL4 dy[2]}] The heights of the two right-hand corners
of the tile, relative to the mesh point.

\item[\texttt{TwoDMeshNode *next}] The next mesh point in the linked
list; \verb@NULL@ if this is the tail.

\item[\texttt{TwoDMeshNode *subMesh}] The head of a linked list of
fine mesh points within the rectangular area spanned by this mesh
point list; \verb@NULL@ if there is no (further) refined mesh for this
location.

\item[\texttt{UINT4 nSub}] The number of fine mesh points in the
above list.  It is an error for \verb@subNum@ to be nonzero and
\verb@subMesh@ to be \verb@NULL@.
\end{description}

******************************************************* </lalLaTeX> */

typedef struct tagTwoDMeshNode {
  REAL4 x, y;
  REAL4 dx;
  REAL4 dy[2];
  struct tagTwoDMeshNode *next;
  struct tagTwoDMeshNode *subMesh;
  UINT4 nSub;
} TwoDMeshNode;

/******************************************************** <lalLaTeX>

\subsubsection*{Structure \texttt{TwoDMeshParamStruc}}
\idx[Type]{TwoDMeshParamStruc}

\noindent This structure stores the parameters required by the
two-dimensional mesh placement functions.  The fields are:

\begin{description}
\item[\texttt{REAL4 domain[2]}] The domain
$[x_\mathrm{min},x_\mathrm{max}]$ spanned by the desired parameter
region.

\item[\texttt{void (*getRange)( LALStatus *, REAL4 [2], REAL4, void
*)}] A function that returns in its second argument the range
$[y_1(x),y_2(x)]$ spanned by the parameter region for a specified $x$,
which is passed in as the third argument.  The fourth argument can be
used to pass function-specific parameters.

\item[\texttt{void *rangeParams}] The parameters to be passed as the
fourth argument of \verb@*getRange()@, above.

\item[\texttt{void (*getMetric)( LALStatus *, REAL4 [3], REAL4 [2],
void *)}] A function that returns in its second argument the
components $g_{xx}$, $g_{yy}$, and $g_{xy}$ (in that order) of the
metric evaluated at a point $(x,y)$, which is passed in as the third
argument.  The fourth argument can be used to pass function-specific
parameters.

\item[\texttt{void *metricParams}] The parameters to be passed as the
fourth argument of \verb@*getMetric()@, above.

\item[\texttt{REAL4 mThresh}] The maximum mismatch $m_\mathrm{thresh}$
desired between any point in the region and the nearest mesh point;
note that the maximum mismatch is equal to 1 minus the minimum match.

\item[\texttt{REAL4 widthMaxFac}] The minimum ratio of mismatch
ellipse width (projected onto the horizontal axis) to column width
that must be maintained throughout the column: if an ellipse falls
below this ratio due to shrinkage or rotation, as in Fig 29.1.b, the
code will try a narrower column.  If set to $\leq1$, the default value
\verb@TWODMESHINTERNALC_WMAXFAC@=$\sqrt[4]{2}$ will be used.

\item[\texttt{REAL4 widthRetryFac}] If the column is determined to be
too wide (e.g.\ due to the value of \verb@widthMaxFac@, above), the
column width will be reduced by the factor \verb@widthRetryFac@.  If
set to $\leq1$, the default value
\verb@TWODMESHINTERNALC_WRETRYFAC@=$\sqrt{2}$ will be used.

\item[\texttt{UINT4 maxColumns}] The maximum number of columns the
mesh placement routine will try before giving up.  If zero, this
number is ignored.

\item[\texttt{UINT4 nIn}] The maximum number of mesh points allowed,
after which the placement routine will quit.  If zero, this number is
ignored.

\item[\texttt{UINT4 nOut}] The number of mesh points added by the
placement routine.  If an error occurs, this will store the number of
mesh points completed before the error.
\end{description}

******************************************************* </lalLaTeX> */

typedef struct tagTwoDMeshParamStruc{
  REAL4 domain[2];
  void (*getRange)( LALStatus *, REAL4 [2], REAL4, void *);
  void *rangeParams;
  void (*getMetric)( LALStatus *, REAL4 [3], REAL4 [2], void *);
  void *metricParams;
  REAL4 mThresh;
  REAL4 widthMaxFac;
  REAL4 widthRetryFac;
  UINT4 maxColumns;
  UINT4 nIn;
  UINT4 nOut;
} TwoDMeshParamStruc;


/******************************************************** <lalLaTeX>

\subsubsection*{Structure \texttt{TwoDColumnParamStruc}}
\idx[Type]{TwoDColumnParamStruc}

\noindent This structure stores additional parameters required when
laying down a single column of a two-dimensional mesh.  The area to be
covered is specified by intersecting the area between two lines with
the parameter space.  If part of a column has already been covered,
one can further restrict the area by specifying a pair of ``clipping
points'' on each vertical line; the area to be covered is then
restricted to lie above the line joining the bottom two corners and
below the line joining the top two corners.  The fields of the
structure are:

\begin{description}
\item[\texttt{REAL4 domain[2]}] The region in $x$ spanned by the
column.  We require that \verb@domain[1]@$>$\verb@domain[0]@.

\item[\texttt{REAL4 leftRange[2]}] The values $y_1(x)$, $y_2(x)$ (in
that order) of the boundary functions at $x=$\verb@domain[0]@.

\item[\texttt{REAL4 rightRange[2]}] The values of $y_1(x)$, $y_2(x)$
(in that order) of the boundary functions at $x=$\verb@domain[1]@.

\item[\texttt{REAL4 leftClip[2]}] The $y$ values of the bottom and top
corners (in that order) of the clipping boundary at
$x=$\verb@domain[0]@.

\item[\texttt{REAL4 rightClip[2]}] The $y$ values of the bottom and
top corners (in that order) of the clipping boundary at
$x=$\verb@domain[1]@.

\item[\texttt{BOOLEAN tooWide}] This is set to 1 if the
column-placement routine determines that the region is too wide to be
covered with a single column of tiles.
\end{description}

******************************************************* </lalLaTeX> */

typedef struct tagTwoDColumnParamStruc {
  REAL4 domain[2];
  REAL4 leftRange[2];
  REAL4 rightRange[2];
  REAL4 leftClip[2];
  REAL4 rightClip[2];
  BOOLEAN tooWide;
} TwoDColumnParamStruc;


/* <lalLaTeX>
\vfill{\footnotesize\input{TwoDMeshHV}}
</lalLaTeX> */


/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{TwoDMeshC}
</lalLaTeX> */
void
LALCreateTwoDMesh( LALStatus          *status,
		   TwoDMeshNode       **mesh,
		   TwoDMeshParamStruc *params );

void
LALDestroyTwoDMesh( LALStatus    *status,
		    TwoDMeshNode **mesh,
		    UINT4        *nFree );

void
LALRefineTwoDMesh( LALStatus    *status,
		   TwoDMeshNode *coarseMesh,
		   TwoDMeshNode *fineMesh );

/* <lalLaTeX>
\newpage\input{TwoDMeshInternalC}
</lalLaTeX> */
void
LALTwoDMesh( LALStatus            *status,
	     TwoDMeshNode         **tail,
	     TwoDMeshParamStruc   *params );

void
LALTwoDColumn( LALStatus            *status,
	       TwoDMeshNode         **tail,
	       TwoDColumnParamStruc *columnParams,
	       TwoDMeshParamStruc   *params );

void
LALTwoDNodeCopy( LALStatus    *status,
		 TwoDMeshNode **new,
		 TwoDMeshNode *old );


void
LALGetNearestMetric( LALStatus *status, REAL4 metric[3], REAL4 position[2], void *params );

void
LALInterpolateMetricGrid( LALStatus *status, REAL4 metric[3], REAL4 position[2], void *params );

/* <lalLaTeX>
\newpage\input{TwoDMeshRangesC}
</lalLaTeX> */
void
LALInterpolateRangePolygon( LALStatus *status, REAL4 range[2], REAL4 x, void *params );

void
LALInterpolateRangeGrid( LALStatus *status, REAL4 range[2], REAL4 x, void *params );

#ifdef __cplusplus
#pragma {
}
#endif

#endif /* _TWODMESH_H */
