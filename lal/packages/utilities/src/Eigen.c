/**************************************** <lalVerbatim file="EigenCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{Eigen.c}}
\label{ss:Eigen.c}

Routines to compute eigenvalues and eigenvectors.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{EigenCP}
\idx{LALSSymmetricEigenVectors()}
\idx{LALDSymmetricEigenVectors()}
\idx{LALSSymmetricEigenValues()}
\idx{LALDSymmetricEigenValues()}

\subsubsection*{Description}

\verb@LALSSymmetricEigenVectors()@ and
\verb@LALDSymmetricEigenVectors()@ compute the eigenvalues and
eigenvectors of the square matrix \verb@*matrix@.  The
eigen\emph{values} are stored in \verb@*values@, which must be
pre-allocated to the same length as each of the dimensions of
\verb@*matrix@.  The eigen\emph{vectors} are computed in-place and
returned in \verb@*matrix@: on return, each column of \verb@*matrix@
stores the eigenvector of the corresponding eigenvalue as a column
vector.  If you don't want the input matrix to be changed, make a copy
of it first.

\verb@LALSSymmetricEigenValues()@ and
\verb@LALDSymmetricEigenValues()@ compute just the eigenvalues of the
square matrix \verb@*matrix@, which is significantly faster than
computing the eigenvectors.  The eigenvalues are stored in
\verb@*values@, which must be pre-allocated to the same length as each
of the dimensions of \verb@*matrix@.  However, \verb@*matrix@ is still
used as auxiliary storage for the in-place algorithm; if you don't
want the input matrix to be changed, make a copy of it first.

\subsubsection*{Algorithm}

A square matrix $\mathsf{M}^a{}_b$ is said to have an eigenvalue
$\lambda$ and associated eigenvector $\mathsf{x}^a$ if the following
equation holds:
\begin{equation}
\label{eq:eigenvalue}
\mathsf{M}^a{}_b \mathsf{x}^b = \lambda\mathsf{x}^a
\end{equation}
Generically an $N\times N$ matrix has $N$ distinct eigenvalues each
with a unique (up to a multiplicative constant) associated
eigenvector.  In certain cases, though, the system is
\emph{degenerate}: that is, some of the $N$ eigenvalues may be the
\emph{same}, and the corresponding eigenvectors may not be unique.

We are concerned with the particular case of real, symmetric matrices,
which occur in a broad class of problems.  These matrices have certain
useful properties.  First, for any Hermitian matrix (including real
symmetric matrices) the eigenvalues are all real.  Second, for any
matrix that commutes with its Hermitian conjugate (including Hermitian
and unitary matrices, or, for real matrices, symmetric and
orthogonal), the eigenvectors ${\mathsf{x}_{(i)}}^a$ of distinct
eigenvalues $\lambda_{(i)}$ are \emph{orthogonal}: that is,
$(\mathsf{x}_{(i)}^T)_a \mathsf{x}_{(j)}{}^a=0$ if
$\lambda_{(i)}\neq\lambda_{(j)}$.  Furthermore, we note that if two or
more linearly independent eigenvectors have the \emph{same}
eigenvalue, then any linear combination of them is also an eigenvector
with that same eigenvalue: we are therefore free to choose a basis of
eigenvectors that is orthonormal.  Thirdly, for any matrix that
commutes with its conjugate, the complete set of eigenvectors spans
the entire $N$-dimensional vector space, so the orthonormal basis
above is a \emph{complete} orthonormal basis.

If we construct a square matrix $\mathsf{X}^a{}_b$ whose columns are
the orthonormal eigenvectors $X^i{}_j={x_{(j)}}^i$, then it is clear
(from the orthonormality of ${\mathsf{x}_{(j)}}^a$) that
$\mathsf{X}^a{}_b$ is orthogonal.  Furtthermore, the eigenvalue
equation~(\ref{eq:eigenvalue}) gives:
\begin{equation}
{(\mathsf{X}^{-1}){}^a{}_b}\mathsf{M}^b{}_c\mathsf{X}^c{}_d =
\left[\begin{array}{cccc}
	\lambda_{(1)} &       0       & \cdots &       0       \\
	      0       & \lambda_{(2)} & \cdots &       0       \\
	   \vdots     &    \vdots     & \ddots &    \vdots     \\
	      0       &       0       & \cdots & \lambda_{(N)}
\end{array}\right] \;,
\end{equation}
where $\lambda_{(i)}$ are the eigenvalues of the corresponding
eigenvectors (with the possibility that some of these eigenvalues have
the same value).  That is, the matrix $\mathsf{M}^b{}_c$ can be
\emph{diagonalized} by an orthogonal similarity transformation, and
the transformation matrix gives us the eigenvectors.  The process of
solving the eigenvalue equation is thus equivalent to the problem of
diagonalizing the matrix.

For a general $N\times N$ symmetric matrix, there is no finite
algorithm that exactly diagonalizes the matrix.  Most numerical
eigenvalue solvers use a process of iterative orthogonal
transformations, where each transformation is chosen to reduce the sum
of the squares of the off-diagonal elements: the matrix uniformly
converges to a diagonal form under successive transformations.  This
convergence can be speeded if one first transforms the matrix into a
\emph{tridiagonal} form (where the only nonzero elements are on the
diagonal or immediately adjacent to it); this transformation
\emph{can} be accomplished deterministically in a finite number of
steps.  This is the approach advocated in Chapter~11
of~\cite{ptvf:1992}, and taken in this module.

The routines in \verb@Eigen.c@ simply call the routines in
\verb@EigenInternal.c@, first to reduce the matrix to tridiagonal
form, and then to reduce it to diagonal form iteratively.  The
routines \verb@LALSSymmetricEigenVectors()@ and
\verb@LALDSymmetricEigenVectors()@ call the routines
\verb@LALSSymmetricToTriDiagonal()@ and
\verb@LALDSymmetricToTriDiagonal()@ to tri-diagonalize the matrix,
then \verb@LALSTriDiagonalToDiagonal()@ and
\verb@LALDTriDiagonalToDiagonal()@ to diagonalize it.  The routines
\verb@LALSSymmetricEigenValues()@ and
\verb@LALDSymmetricEigenValues()@ instead call
\verb@LALSSymmetricToTriDiagonal2()@$\Rightarrow$\verb@LALSTriDiagonalToDiagonal2()@
or
\verb@LALDSymmetricToTriDiagonal2()@$\Rightarrow$\verb@LALDTriDiagonalToDiagonal2()@,
which discard information about the transformation matrix itself and
thus do not give the eigenvectors; however, they execute much faster.

\subsubsection*{Uses}
\begin{verbatim}
LALSCreateVector()              LALSDestroyVector()
LALDCreateVector()              LALDDestroyVector()
LALSSymmetricToTriDiagonal()    LALSTriDiagonalToDiagonal()
LALDSymmetricToTriDiagonal()    LALDTriDiagonalToDiagonal()
LALSSymmetricToTriDiagonal2()   LALSTriDiagonalToDiagonal2()
LALDSymmetricToTriDiagonal2()   LALDTriDiagonalToDiagonal2()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{EigenCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/MatrixUtils.h>

NRCSID( EIGENC, "$Id$" );

/* <lalVerbatim file="EigenCP"> */
void
LALSSymmetricEigenVectors( LALStatus *stat, REAL4Vector *values, REAL4Array *matrix )
{ /* </lalVerbatim> */
  REAL4Vector *offDiag = NULL; /* off-diagonal line of
                                  tri-diagonalized matrix */

  INITSTATUS( stat, "LALSSymmetricEigenVectors", EIGENC );
  ATTATCHSTATUSPTR( stat );

  /* Check dimension length.  All other argument testing is done by
     the subroutines. */
  ASSERT( values, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( values->length, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );

  /* Allocate an off-diagonal vector for the tri-diagonal matrix. */
  TRY( LALSCreateVector( stat->statusPtr, &offDiag, values->length ),
       stat );

  /* Call the subroutines. */
  LALSSymmetricToTriDiagonal( stat->statusPtr, values, matrix,
			      offDiag );
  BEGINFAIL( stat ) {
    TRY( LALSDestroyVector( stat->statusPtr, &offDiag ), stat );
  } ENDFAIL( stat );
  LALSTriDiagonalToDiagonal( stat->statusPtr, values, matrix,
			     offDiag );
  BEGINFAIL( stat ) {
    TRY( LALSDestroyVector( stat->statusPtr, &offDiag ), stat );
  } ENDFAIL( stat );

  /* Clean up. */
  TRY( LALSDestroyVector( stat->statusPtr, &offDiag ), stat );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="EigenCP"> */
void
LALSSymmetricEigenValues( LALStatus *stat, REAL4Vector *values, REAL4Array *matrix )
{ /* </lalVerbatim> */
  REAL4Vector *offDiag = NULL; /* off-diagonal line of
                                  tri-diagonalized matrix */

  INITSTATUS( stat, "LALSSymmetricEigenValues", EIGENC );
  ATTATCHSTATUSPTR( stat );

  /* Check dimension length.  All other argument testing is done by
     the subroutines. */
  ASSERT( values, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( values->length, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );

  /* Allocate an off-diagonal vector for the tri-diagonal matrix. */
  TRY( LALSCreateVector( stat->statusPtr, &offDiag, values->length ),
       stat );

  /* Call the subroutines. */
  LALSSymmetricToTriDiagonal2( stat->statusPtr, values, matrix,
			       offDiag );
  BEGINFAIL( stat ) {
    TRY( LALSDestroyVector( stat->statusPtr, &offDiag ), stat );
  } ENDFAIL( stat );
  LALSTriDiagonalToDiagonal2( stat->statusPtr, values, matrix,
			      offDiag );
  BEGINFAIL( stat ) {
    TRY( LALSDestroyVector( stat->statusPtr, &offDiag ), stat );
  } ENDFAIL( stat );

  /* Clean up. */
  TRY( LALSDestroyVector( stat->statusPtr, &offDiag ), stat );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="EigenCP"> */
void
LALDSymmetricEigenVectors( LALStatus *stat, REAL8Vector *values, REAL8Array *matrix )
{ /* </lalVerbatim> */
  REAL8Vector *offDiag = NULL; /* off-diagonal line of
                                  tri-diagonalized matrix */

  INITSTATUS( stat, "LALSSymmetricEigenVectors", EIGENC );
  ATTATCHSTATUSPTR( stat );

  /* Check dimension length.  All other argument testing is done by
     the subroutines. */
  ASSERT( values, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( values->length, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );

  /* Allocate an off-diagonal vector for the tri-diagonal matrix. */
  TRY( LALDCreateVector( stat->statusPtr, &offDiag, values->length ),
       stat );

  /* Call the subroutines. */
  LALDSymmetricToTriDiagonal( stat->statusPtr, values, matrix,
			      offDiag );
  BEGINFAIL( stat ) {
    TRY( LALDDestroyVector( stat->statusPtr, &offDiag ), stat );
  } ENDFAIL( stat );
  LALDTriDiagonalToDiagonal( stat->statusPtr, values, matrix,
			     offDiag );
  BEGINFAIL( stat ) {
    TRY( LALDDestroyVector( stat->statusPtr, &offDiag ), stat );
  } ENDFAIL( stat );

  /* Clean up. */
  TRY( LALDDestroyVector( stat->statusPtr, &offDiag ), stat );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="EigenCP"> */
void
LALDSymmetricEigenValues( LALStatus *stat, REAL8Vector *values, REAL8Array *matrix )
{ /* </lalVerbatim> */
  REAL8Vector *offDiag = NULL; /* off-diagonal line of
                                  tri-diagonalized matrix */

  INITSTATUS( stat, "LALSSymmetricEigenValues", EIGENC );
  ATTATCHSTATUSPTR( stat );

  /* Check dimension length.  All other argument testing is done by
     the subroutines. */
  ASSERT( values, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( values->length, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );

  /* Allocate an off-diagonal vector for the tri-diagonal matrix. */
  TRY( LALDCreateVector( stat->statusPtr, &offDiag, values->length ),
       stat );

  /* Call the subroutines. */
  LALDSymmetricToTriDiagonal2( stat->statusPtr, values, matrix,
			       offDiag );
  BEGINFAIL( stat ) {
    TRY( LALDDestroyVector( stat->statusPtr, &offDiag ), stat );
  } ENDFAIL( stat );
  LALDTriDiagonalToDiagonal2( stat->statusPtr, values, matrix,
			      offDiag );
  BEGINFAIL( stat ) {
    TRY( LALDDestroyVector( stat->statusPtr, &offDiag ), stat );
  } ENDFAIL( stat );

  /* Clean up. */
  TRY( LALDDestroyVector( stat->statusPtr, &offDiag ), stat );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
