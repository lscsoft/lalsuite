/*
*  Copyright (C) 2007 Jolien Creighton, Teviet Creighton
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

#ifndef _FLATMESH_H
#define _FLATMESH_H

#include <lal/LALStdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \author Creighton, T. D.
 * \defgroup FlatMesh_h Header FlatMesh.h
 * \ingroup pkg_pulsarCovering
 * \brief Provides routines to place search meshes for parameter spaces with
 * constant parameter metrics.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/FlatMesh.h>
 * \endcode
 *
 * This header covers routines that lay out a mesh of points on
 * an \f$n\f$-dimensional parameter space \f$\{\mathsf{x}^a\}\f$.  The mesh
 * points are placed in a unit-cube lattice in orthonormalized
 * coordinates \f$\mathsf{y}^a=\mathsf{T}^a{}_b\mathsf{x}^b\f$, where
 * \f$\mathsf{T}^a{}_b\f$ is an arbitrary but constant transformation matrix.
 * This describes a fairly general procedure, but within LAL these
 * routines are used for a specific purpose: to define a set of parameter
 * values required to search a space of parameterized signals.
 *
 * In most optimal signal-detection strategies, each signal that one
 * might detect is parameterized by an \f$n\f$-tuple of parameters
 * \f$\mathsf{x}^a=(x_1,\ldots,x_n)\f$.  To detect a signal whose parameters
 * are not known in advance, a typical search algorithm performs a series
 * of targeted searches over a discrete set of parameter values
 * \f$\mathsf{x}^a_{(1)},\mathsf{x}^a_{(2)},\ldots\f$.  If the unknown
 * signal's parameters happen to match one of the target points, it will
 * be detected with maximum signal strength; otherwise, its effective
 * signal strength is reduced according to the ``distance'' between its
 * parameter values and the closest search point.  The \e mismatch
 * \f$m(\mathsf{x}^a,\Delta\mathsf{x}^a)\f$ is defined as the fractional
 * reduction in effective signal strength of a signal with parameters
 * \f$\mathsf{x}^a+\Delta\mathsf{x}^a\f$ in a search targeted at signal
 * parameters \f$\mathsf{x}^a\f$.  The mismatch can define a local distance
 * measure on the parameter space, since for small values it approaches a
 * quadratic form:
 * \f[
 * m(\mathsf{x}^a,\Delta\mathsf{x}^a) = \mathsf{g}_{bc}(\mathsf{x}^a)
 * \Delta\mathsf{x}^b\Delta\mathsf{x}^c
 * + O(\Delta\mathsf{x}^a)^3 \; .
 * \f]
 * The matrix \f$\mathsf{g}_{bc}\f$ is called the mismatch \e metric, and
 * in general can vary as a function of the central target point
 * \f$\mathsf{x}^a\f$.
 *
 * One of the main goals in searching the signal parameter space is to
 * choose a set of target points that is as small as possible, while
 * still ensuring that an unknown signal will lie close enough to a
 * target point that it will lose no more than some fraction
 * \f$m_\mathrm{thresh}\f$ of its signal strength.  If the mismatch metric
 * \f$\mathsf{g_{ab}}\f$ is constant, this is relatively simple: the search
 * points can be placed on a mesh that is rectilinear and
 * uniformly-spaced in the eigenspace of \f$\mathsf{g_{ab}}\f$.  This is
 * described below.
 *
 * Let \f$\lambda^{(1)},\ldots,\lambda^{(n)}\f$ be the \f$n\f$ eigenvalues of
 * \f$\mathsf{g}_{ab}\f$, and \f$\mathsf{e}^a_{(1)},\ldots,\mathsf{e}^a_{(n)}\f$
 * be the corresponding unit eigenvectors.  Then the vectors
 * \f$\{\mathsf{e}^a_{(i)}/\sqrt{\lambda_{(i)}}\}\f$ define a new coordinate
 * basis in which coordinate distances correspond to metric distances.
 * The simplest covering of the parameter space is to lay out a cubic
 * mesh in this coordinate basis.  If \f$s\f$ is the side length of one of
 * these cubes, then the maximum mismatch (in the quadratic
 * approximation) between a point in the interior and a vertex is
 * \f$(s/2)\sqrt{n}\f$.  We wish this to be no greater than some given
 * threshold value \f$m_\mathrm{thresh}\f$.  This means that the eigenvectors
 * pointing from a mesh point to its immediate neighbours are of the form
 * \f$\pm2m_\mathrm{thresh}[n\lambda_{(i)}]^{-1/2}\,\mathsf{e}^a_{(i)}\f$.
 *
 * Let us define a transformation matrix \f$\mathsf{M}^a{}_b\f$ by:
 * \f[
 * M^i{}_j = 2m_\mathrm{thresh}[n\lambda_{(i)}]^{-1/2}\,e^i_{(j)} \; .
 * \f]
 * Then the parameters \f$\mathsf{x}^a\f$ and the orthonormalized coordinates
 * \f$\mathsf{y}^b\f$ are related by:
 * \f{eqnarray}{
 * \mathsf{x}^a & = & \mathsf{M}^a{}_b \mathsf{y}^b \; ,\\
 * \mathsf{y}^b & = & \mathsf{(M^{-1})}^b{}_a \mathsf{x}^a \; .
 * \f}
 * The search mesh can thus be placed as a unit-cube lattice in the
 * \f$\mathsf{y}^b\f$ coordinate basis and then transformed back to find the
 * mesh points in the \f$\mathsf{x}^a\f$ coordinates.
 *
 * This header and its associated modules are placed in the \c pulsar
 * package because they were originally intended for use in targeted
 * pulsar searches, where the sky position is known but the frequency and
 * spindown are not.  In a search over the Taylor coefficients of the
 * frequency function, the associated mismatch metric is nearly constant,
 * with corrections on the order of the observation time over the
 * spindown timescale.
 *
 */
/*@{*/

/** \name Error Codes */
/*@{*/
#define FLATMESHH_ENUL 1
#define FLATMESHH_EOUT 2
#define FLATMESHH_EMEM 3
#define FLATMESHH_EDIM 4
#define FLATMESHH_ELEN 5

#define FLATMESHH_MSGENUL "Unexpected null pointer in arguments"
#define FLATMESHH_MSGEOUT "Output handle points to a non-null pointer"
#define FLATMESHH_MSGEMEM "Memory allocation error"
#define FLATMESHH_MSGEDIM "Inconsistent parameter space dimension"
#define FLATMESHH_MSGELEN "Too few points specified"
/*@}*/

/**
 * \brief This structure stores and passes parameters for computing a
 * search mesh by the routines in the module \ref FlatMesh_h : it
 * defines the transformation matrices between signal parameters and
 * orthonormalized coordinates, a rectangular volume containing the
 * region to be covered, and a function defining the covered region as a
 * subset of the rectangular volume.
 */
typedef struct tagFlatMeshParamStruc
{
  REAL4VectorSequence *matrix;    /**< The matrix \f$\mathsf{M}^a{}_b\f$ (above) that transforms from orthonormalized to parameter coordinates */
  REAL4VectorSequence *matrixInv; /**< The inverse of <tt>*matrix</tt> */
  REAL4Vector *xMin;	/**< A vector defining one corner of a rectangular region in parameter space that contains the region to be
                         * covered; by convention, <tt>*xMin</tt> is the corner where each parameter has its minimum value. */
  REAL4Vector *xMax;	/**< A vector defining the opposite corner of the rectangular region above; by convention, <tt>*xMax</tt>
                         * is the corner where each parameter has its maximum value. */
  void (*intersection)( LALStatus *, REAL4VectorSequence *, REAL4VectorSequence * ); /**< The function that restricts the mesh to
                                                                                      * cover only the desired search region. */
  REAL4VectorSequence *controlPoints;/**< A set of ``control points'' to be passed as the second argument to <tt>*intersection</tt>
                                      * (above), defining the size and shape of the desired search region */
} FlatMeshParamStruc;


/* Function prototypes. */
void
LALCreateFlatMesh( LALStatus           *status,
		   REAL4VectorSequence **mesh,
		   FlatMeshParamStruc  *params );

void
LALRectIntersect( LALStatus           *status,
		  REAL4VectorSequence *mesh,
		  REAL4VectorSequence *controlPoints );

/*@}*/

#ifdef __cplusplus
}
#endif

#endif /* _FLATMESH_H */
