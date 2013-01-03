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

#ifndef _TWODMESH_H
#define _TWODMESH_H

#include <lal/LALStdlib.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
\author Creighton, T. D.
\defgroup TwoDMesh_h Header TwoDMesh.h
\ingroup pkg_pulsarCovering
\brief Provides routines to place search meshes for two-dimensional parameter spaces with varying metric.

\heading{Synopsis}
\code
#include <lal/TwoDMesh.h>
\endcode

This header covers routines that lay out a mesh of points on
an 2-dimensional parameter space \f$\{(x,y)\}\f$, placed such that no
point in the space lies further than some maximum proper distance
\f$m_\mathrm{thresh}\f$ from a mesh point.

The intended purpose of these routines is to place a set of ``target''
search points over a parameter space, in order to detect signals with
unknown parameters.  The formalism for defining a proper distance
metric on the parameter space is defined in
\ref FlatMesh_h.  However, whereas the routines under
\ref FlatMesh_h require the metric \f$\mathsf{g}_{ab}\f$ to be constant
over the parameter space, the routines under this header only treat
\f$\mathsf{g}_{ab}\f$ as constant over distances \f$\lesssim m_\mathrm{thresh}\f$.

\floatfig{H,fig_tiling}
\image html pulsar_tiling.png "Fig. [fig_tiling]: Mesh placement using parallelogram tiling. (a) The left and right sides of a tile are required to be vertical; the top and bottom sides can tilt to maximize the tile area. (b) Tiles can be stacked in fixed-width columns, even as the elliptical contours change.  (c) Extra overlapping tiles are sometimes required at the corners of columns."
\image latex pulsar_tiling.pdf "Mesh placement using parallelogram tiling. (a) The left and right sides of a tile are required to be vertical; the top and bottom sides can tilt to maximize the tile area. (b) Tiles can be stacked in fixed-width columns, even as the elliptical contours change.  (c) Extra overlapping tiles are sometimes required at the corners of columns." width=\textwidth

Since the metric is treated as constant over distances \f$\lesssim
m_\mathrm{thresh}\f$, this distance defines an elliptical contour around
any mesh point.  We define a ``tile'' as a parallelogram inscribed
within the ellipse, with its left and right sides aligned with the \f$y\f$
axis.  This is shown in Fig.\figref{fig_tiling} (a), above.  A ``column''
is a set of tiles of constant horizontal width stacked one on top of
the other, as shown in Fig.\figref{fig_tiling} (b).  As the metric
changes over space, the vertical height and tilt of the tiles in a
column may change, so long as their width remains fixed; we note that
if the tilt changes, the tiles will overlap slightly to ensure
complete coverage.  Finally, the boundary of the parameter space may
extend outside the ``corners'' of the column, crossing the end of a
tile between its centre and its edge, as shown in
Fig.\figref{fig_tiling} (c).  These triangular corners can be covered
with one or more extra overlapping tiles of reduced width.

In a parameter space with constant metric, the tile area is maximized
(and the number of covering tiles minimized) when the column width is
\f$\sqrt{2}\f$ times smaller than the projected horizontal width of the
ellipses.  When the ellipses vary, it is generally best to determine
the column width from the \e narrowest ellipse in a column, to
avoid singular effects when tile widths approach the ellipse widths
and become infinitesimally high.

For the column-placement algorithm to work effectively, we require
that the parameter space be representable as a range
\f$y\in[y_1(x),y_2(x)]\f$ between two single-valued functions defined on a
domain \f$x\in[x_\mathrm{min},x_\mathrm{max}]\f$.  If a desired search
region is too complicated to express this way (e.g.\ it has
disconnected regions, or ``branches'' where a vertical line intersects
the boundary more than twice), then one should divide the region up
into subregions with well-behaved boundary functions and tile these
subregions separately.

This header and its associated modules are placed in the \c pulsar
package because they were originally intended for use in searches over
sky position, but they can be used generically for any two-dimensional
parameter space search where the metric is not too poorly behaved.
*/
/*@{*/

/** \name Error Codes */
/*@{*/
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
/*@}*/

/** This structure represents a single node in a linked list of
 * mesh points, specified in the coordinate system used to place it
 */
typedef struct tagTwoDMeshNode
{
  REAL4 x, y; 	/**< The coordinates of the mesh point */
  REAL4 dx;	/**< The half-width of the tile centred on the mesh point */
  REAL4 dy[2];	/**< The heights of the two right-hand corners of the tile, relative to the mesh point */
  struct tagTwoDMeshNode *next;	/**< The next mesh point in the linked list; \c NULL if this is the tail */
  struct tagTwoDMeshNode *subMesh;	/**< The head of a linked list of fine mesh points within the rectangular
                                         * area spanned by this mesh point list; \c NULL if there is no (further)
                                         * refined mesh for this location
                                         */
  UINT4 nSub;	/**< The number of fine mesh points in the above list.  It is an error for \c subNum to be nonzero
                 * and \c subMesh to be \c NULL
                 */
} TwoDMeshNode;

/** This structure stores the parameters required by the
 * two-dimensional mesh placement functions.
 */
typedef struct tagTwoDMeshParamStruc
{
  REAL4 domain[2];	/**< The domain \f$[x_\mathrm{min},x_\mathrm{max}]\f$ spanned by the desired parameter region */
  void (*getRange)( LALStatus *, REAL4 [2], REAL4, void *); /**< A function that returns in its second argument the range
                                                             * \f$[y_1(x),y_2(x)]\f$ spanned by the parameter region for a specified \f$x\f$,
                                                             * which is passed in as the third argument; the fourth argument can be
                                                             * used to pass function-specific parameters.
                                                             */
  void *rangeParams;	/**<  The parameters to be passed as the fourth argument of <tt>*getRange()</tt>, above. */
  void (*getMetric)( LALStatus *, REAL4 [3], REAL4 [2], void *); /**< A function that returns in its second argument the
                                                                  * components \f$g_{xx}\f$, \f$g_{yy}\f$, and \f$g_{xy}\f$ (in that order) of the
                                                                  * metric evaluated at a point \f$(x,y)\f$, which is passed in as the third
                                                                  * argument; the fourth argument can be used to pass function-specific parameters
                                                                  */
  void *metricParams;	/**< The parameters to be passed as the fourth argument of <tt>*getMetric()</tt>, above */
  REAL4 mThresh;	/**<  The maximum mismatch \f$m_\mathrm{thresh}\f$ desired between any point in the region and the nearest mesh point;
                         * note that the maximum mismatch is equal to 1 minus the minimum match.
                         */
  REAL4 widthMaxFac;	/**<  The minimum ratio of mismatch ellipse width (projected onto the horizontal axis) to column width
                         * that must be maintained throughout the column: if an ellipse falls
                         * below this ratio due to shrinkage or rotation, as in Fig.\figref{fig_tiling} (b), the
                         * code will try a narrower column; if set to \f$\leq1\f$, the default value
                         * \c TWODMESHINTERNALC_WMAXFAC=\f$\sqrt[4]{2}\f$ will be used.
                         */
  REAL4 widthRetryFac;	/**< If the column is determined to be too wide (e.g.\ due to the value of \c widthMaxFac, above), the
                         * column width will be reduced by the factor \c widthRetryFac; if set to \f$\leq1\f$, the default value
                         * \c TWODMESHINTERNALC_WRETRYFAC=\f$\sqrt{2}\f$ will be used.
                         */
  UINT4 maxColumns;	/**< The maximum number of columns the mesh placement routine will try before giving up.  If zero, this number is ignored. */
  UINT4 nIn;		/**< The maximum number of mesh points allowed, after which the placement routine will quit.  If zero, this number is ignored. */
  UINT4 nOut;		/**< The number of mesh points added by the placement routine; if an error occurs, this will store the number of mesh points completed before the error */
} TwoDMeshParamStruc;


/** This structure stores additional parameters required when
 * laying down a single column of a two-dimensional mesh.  The area to be
 * covered is specified by intersecting the area between two lines with
 * the parameter space.  If part of a column has already been covered,
 * one can further restrict the area by specifying a pair of ``clipping
 * points'' on each vertical line; the area to be covered is then
 * restricted to lie above the line joining the bottom two corners and
 * below the line joining the top two corners
 */
typedef struct tagTwoDColumnParamStruc
{
  REAL4 domain[2];	/**< The region in \f$x\f$ spanned by the column; We require that <tt>domain[1]</tt>\f$>\f$<tt>domain[0]</tt> */
  REAL4 leftRange[2];	/**< The values \f$y_1(x)\f$, \f$y_2(x)\f$ (in that order) of the boundary functions at \f$x=\f$<tt>domain[0]</tt> */
  REAL4 rightRange[2];  /**< The values of \f$y_1(x)\f$, \f$y_2(x)\f$ (in that order) of the boundary functions at \f$x=\f$<tt>domain[1]</tt> */
  REAL4 leftClip[2];	/**< The \f$y\f$ values of the bottom and top corners (in that order) of the clipping boundary at \f$x=\f$<tt>domain[0]</tt>. */
  REAL4 rightClip[2];	/**< The \f$y\f$ values of the bottom and top corners (in that order) of the clipping boundary at \f$x=\f$<tt>domain[1]</tt> */
  BOOLEAN tooWide;	/**< This is set to 1 if the column-placement routine determines that the region is too wide to be covered with a single column of tiles. */
} TwoDColumnParamStruc;


/* Function prototypes. */
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
		 TwoDMeshNode **new_,
		 TwoDMeshNode *old );


void
LALGetNearestMetric( LALStatus *status, REAL4 metric[3], REAL4 position[2], void *params );

void
LALInterpolateMetricGrid( LALStatus *status, REAL4 metric[3], REAL4 position[2], void *params );




void
LALInterpolateRangePolygon( LALStatus *status, REAL4 range[2], REAL4 x, void *params );

void
LALInterpolateRangeGrid( LALStatus *status, REAL4 range[2], REAL4 x, void *params );

/*@}*/

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _TWODMESH_H */
