/************************************* <lalVerbatim file="FlatMeshHV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\section{Header \texttt{FlatMesh.h}}
\label{s:FlatMesh.h}

Provides routines to place search meshes for parameter spaces with
constant parameter metrics.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/FlatMesh.h>
\end{verbatim}

\noindent This header covers routines that lay out a mesh of points on
an $n$-dimensional parameter space $\{\mathsf{x}^a\}$.  The mesh
points are placed in a unit-cube lattice in orthonormalized
coordinates $\mathsf{y}^a=\mathsf{T}^a{}_b\mathsf{x}^b$, where
$\mathsf{T}^a{}_b$ is an arbitrary but constant transformation matrix.
This describes a fairly general procedure, but within LAL these
routines are used for a specific purpose: to define a set of parameter
values required to search a space of parameterized signals.

In most optimal signal-detection strategies, each signal that one
might detect is parameterized by an $n$-tuple of parameters
$\mathsf{x}^a=(x_1,\ldots,x_n)$.  To detect a signal whose parameters
are not known in advance, a typical search algorithm performs a series
of targeted searches over a discrete set of parameter values
$\mathsf{x}^a_{(1)},\mathsf{x}^a_{(2)},\ldots$.  If the unknown
signal's parameters happen to match one of the target points, it will
be detected with maximum signal strength; otherwise, its effective
signal strength is reduced according to the ``distance'' between its
parameter values and the closest search point.  The \emph{mismatch}
$m(\mathsf{x}^a,\Delta\mathsf{x}^a)$ is defined as the fractional
reduction in effective signal strength of a signal with parameters
$\mathsf{x}^a+\Delta\mathsf{x}^a$ in a search targeted at signal
parameters $\mathsf{x}^a$.  The mismatch can define a local distance
measure on the parameter space, since for small values it approaches a
quadratic form:
$$
m(\mathsf{x}^a,\Delta\mathsf{x}^a) = \mathsf{g}_{bc}(\mathsf{x}^a)
	\Delta\mathsf{x}^b\Delta\mathsf{x}^c
	+ O(\Delta\mathsf{x}^a)^3 \; .
$$
The matrix $\mathsf{g}_{bc}$ is called the mismatch \emph{metric}, and
in general can vary as a function of the central target point
$\mathsf{x}^a$.

One of the main goals in searching the signal parameter space is to
choose a set of target points that is as small as possible, while
still ensuring that an unknown signal will lie close enough to a
target point that it will lose no more than some fraction
$m_\mathrm{thresh}$ of its signal strength.  If the mismatch metric
$\mathsf{g_{ab}}$ is constant, this is relatively simple: the search
points can be placed on a mesh that is rectilinear and
uniformly-spaced in the eigenspace of $\mathsf{g_{ab}}$.  This is
described below.

Let $\lambda^{(1)},\ldots,\lambda^{(n)}$ be the $n$ eigenvalues of
$\mathsf{g}_{ab}$, and $\mathsf{e}^a_{(1)},\ldots,\mathsf{e}^a_{(n)}$
be the corresponding unit eigenvectors.  Then the vectors
$\{\mathsf{e}^a_{(i)}/\sqrt{\lambda_{(i)}}\}$ define a new coordinate
basis in which coordinate distances correspond to metric distances.
The simplest covering of the parameter space is to lay out a cubic
mesh in this coordinate basis.  If $s$ is the side length of one of
these cubes, then the maximum mismatch (in the quadratic
approximation) between a point in the interior and a vertex is
$(s/2)\sqrt{n}$.  We wish this to be no greater than some given
threshold value $m_\mathrm{thresh}$.  This means that the eigenvectors
pointing from a mesh point to its immediate neighbours are of the form
$\pm2m_\mathrm{thresh}[n\lambda_{(i)}]^{-1/2}\,\mathsf{e}^a_{(i)}$.

Let us define a transformation matrix $\mathsf{M}^a{}_b$ by:
$$
M^i{}_j = 2m_\mathrm{thresh}[n\lambda_{(i)}]^{-1/2}\,e^i_{(j)} \; .
$$
Then the parameters $\mathsf{x}^a$ and the orthonormalized coordinates
$\mathsf{y}^b$ are related by:
\begin{eqnarray}
\mathsf{x}^a & = & \mathsf{M}^a{}_b \mathsf{y}^b \; , \nonumber\\
\mathsf{y}^b & = & \mathsf{(M^{-1})}^b{}_a \mathsf{x}^a \; . \nonumber
\end{eqnarray}
The search mesh can thus be placed as a unit-cube lattice in the
$\mathsf{y}^b$ coordinate basis and then transformed back to find the
mesh points in the $\mathsf{x}^a$ coordinates.

This header and its associated modules are placed in the \verb@pulsar@
package because they were originally intended for use in targeted
pulsar searches, where the sky position is known but the frequency and
spindown are not.  In a search over the Taylor coefficients of the
frequency function, the associated mismatch metric is nearly constant,
with corrections on the order of the observation time over the
spindown timescale.

******************************************************* </lalLaTeX> */

#ifndef _FLATMESH_H
#define _FLATMESH_H

#include <lal/LALStdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

NRCSID(FLATMESHH,"$Id$");

/********************************************************** <lalLaTeX>
\subsection*{Error conditions}
****************************************** </lalLaTeX><lalErrTable> */
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
/******************************************** </lalErrTable><lalLaTeX>

\subsection*{Types}

\subsubsection*{Structure \texttt{FlatMeshParamStruc}}
\idx[Type]{FlatMeshParamStruc}

\noindent This structure stores and passes parameters for computing a
search mesh by the routines in the module \verb@FlatMesh.c@: it
defines the transformation matrices between signal parameters and
orthonormalized coordinates, a rectangular volume containing the
region to be covered, and a function defining the covered region as a
subset of the rectangular volume.  The fields are:

\begin{description}
\item[\texttt{REAL4VectorSequence *matrix}] The matrix
$\mathsf{M}^a{}_b$ (above) that transforms from orthonormalized to
parameter coordinates.

\item[\texttt{REAL4VectorSequence *matrixInv}] The inverse of
\verb@*matrix@.

\item[\texttt{REAL4Vector *xMin}] A vector defining one corner of a
rectangular region in parameter space that contains the region to be
covered.  By convention, \verb@*xMin@ is the corner where each
parameter has its minimum value.

\item[\texttt{REAL4Vector *xMax}] A vector defining the opposite
corner of the rectangular region above.  By convention, \verb@*xMax@
is the corner where each parameter has its maximum value.

\item[\texttt{void (*intersection)( LALStatus *, REAL4VectorSequence
*, REAL4VectorSequence * )}] The function that restricts the mesh to
cover only the desired search region.

\item[\texttt{REAL4VectorSequence *controlPoints}] A set of ``control
points'' to be passed as the second argument to \verb@*intersection@
(above), defining the size and shape of the desired search region.
\end{description}

******************************************************* </lalLaTeX> */

typedef struct tagFlatMeshParamStruc{
  REAL4VectorSequence *matrix;    /* y->x tranformation matrix */
  REAL4VectorSequence *matrixInv; /* x->y tranformation matrix */
  REAL4Vector *xMin; /* minimum values of search parameters */
  REAL4Vector *xMax; /* maximum values of search parameters */
  void (*intersection)( LALStatus *, REAL4VectorSequence *,
			REAL4VectorSequence * );
  /* routine to restrict mesh to the search region */
  REAL4VectorSequence *controlPoints;
  /* control points to *intersection() */
} FlatMeshParamStruc;


/* <lalLaTeX>
\vfill{\footnotesize\input{FlatMeshHV}}
</lalLaTeX> */


/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{FlatMeshC}
</lalLaTeX> */
void
LALCreateFlatMesh( LALStatus           *status,
		   REAL4VectorSequence **mesh,
		   FlatMeshParamStruc  *params );

void
LALRectIntersect( LALStatus           *status,
		  REAL4VectorSequence *mesh,
		  REAL4VectorSequence *controlPoints );

/* <lalLaTeX>
\newpage\input{FlatMeshTestC}

\newpage\input{DirectedMeshTestC}
</lalLaTeX> */

#ifdef __cplusplus
}
#endif

#endif /* _FLATMESH_H */
