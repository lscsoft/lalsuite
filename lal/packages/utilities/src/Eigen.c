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


#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/MatrixUtils.h>


/**
\defgroup Eigen_c Module Eigen.c
\ingroup MatrixUtils_h
\author Creighton, T. D.

\brief Routines to compute eigenvalues and eigenvectors.

\heading{Description}

<tt>LALSSymmetricEigenVectors()</tt> and
<tt>LALDSymmetricEigenVectors()</tt> compute the eigenvalues and
eigenvectors of the square matrix <tt>*matrix</tt>.  The
eigen\e values are stored in <tt>*values</tt>, which must be
pre-allocated to the same length as each of the dimensions of
<tt>*matrix</tt>.  The eigen\e vectors are computed in-place and
returned in <tt>*matrix</tt>: on return, each column of <tt>*matrix</tt>
stores the eigenvector of the corresponding eigenvalue as a column
vector.  If you don't want the input matrix to be changed, make a copy
of it first.

<tt>LALSSymmetricEigenValues()</tt> and
<tt>LALDSymmetricEigenValues()</tt> compute just the eigenvalues of the
square matrix <tt>*matrix</tt>, which is significantly faster than
computing the eigenvectors.  The eigenvalues are stored in
<tt>*values</tt>, which must be pre-allocated to the same length as each
of the dimensions of <tt>*matrix</tt>.  However, <tt>*matrix</tt> is still
used as auxiliary storage for the in-place algorithm; if you don't
want the input matrix to be changed, make a copy of it first.

\heading{Algorithm}

A square matrix \f$\mathsf{M}^a{}_b\f$ is said to have an eigenvalue
\f$\lambda\f$ and associated eigenvector \f$\mathsf{x}^a\f$ if the following
equation holds:
\anchor eq_eigenvalue \f{equation}{
\tag{eq_eigenvalue}
\mathsf{M}^a{}_b \mathsf{x}^b = \lambda\mathsf{x}^a
\f}
Generically an \f$N\times N\f$ matrix has \f$N\f$ distinct eigenvalues each
with a unique (up to a multiplicative constant) associated
eigenvector.  In certain cases, though, the system is
\e degenerate: that is, some of the \f$N\f$ eigenvalues may be the
\e same, and the corresponding eigenvectors may not be unique.

We are concerned with the particular case of real, symmetric matrices,
which occur in a broad class of problems.  These matrices have certain
useful properties.  First, for any Hermitian matrix (including real
symmetric matrices) the eigenvalues are all real.  Second, for any
matrix that commutes with its Hermitian conjugate (including Hermitian
and unitary matrices, or, for real matrices, symmetric and
orthogonal), the eigenvectors \f${\mathsf{x}_{(i)}}^a\f$ of distinct
eigenvalues \f$\lambda_{(i)}\f$ are \e orthogonal: that is,
\f$(\mathsf{x}_{(i)}^T)_a \mathsf{x}_{(j)}{}^a=0\f$ if
\f$\lambda_{(i)}\neq\lambda_{(j)}\f$.  Furthermore, we note that if two or
more linearly independent eigenvectors have the \e same
eigenvalue, then any linear combination of them is also an eigenvector
with that same eigenvalue: we are therefore free to choose a basis of
eigenvectors that is orthonormal.  Thirdly, for any matrix that
commutes with its conjugate, the complete set of eigenvectors spans
the entire \f$N\f$-dimensional vector space, so the orthonormal basis
above is a \e complete orthonormal basis.

If we construct a square matrix \f$\mathsf{X}^a{}_b\f$ whose columns are
the orthonormal eigenvectors \f$X^i{}_j={x_{(j)}}^i\f$, then it is clear
(from the orthonormality of \f${\mathsf{x}_{(j)}}^a\f$) that
\f$\mathsf{X}^a{}_b\f$ is orthogonal.  Furtthermore, the eigenvalue
equation.\eqref{eq_eigenvalue} gives:
\f{equation}{
{(\mathsf{X}^{-1}){}^a{}_b}\mathsf{M}^b{}_c\mathsf{X}^c{}_d =
\left[\begin{array}{cccc}
	\lambda_{(1)} &       0       & \cdots &       0       \\
	      0       & \lambda_{(2)} & \cdots &       0       \\
	   \vdots     &    \vdots     & \ddots &    \vdots     \\
	      0       &       0       & \cdots & \lambda_{(N)}
\end{array}\right] \;,
\f}
where \f$\lambda_{(i)}\f$ are the eigenvalues of the corresponding
eigenvectors (with the possibility that some of these eigenvalues have
the same value).  That is, the matrix \f$\mathsf{M}^b{}_c\f$ can be
\e diagonalized by an orthogonal similarity transformation, and
the transformation matrix gives us the eigenvectors.  The process of
solving the eigenvalue equation is thus equivalent to the problem of
diagonalizing the matrix.

For a general \f$N\times N\f$ symmetric matrix, there is no finite
algorithm that exactly diagonalizes the matrix.  Most numerical
eigenvalue solvers use a process of iterative orthogonal
transformations, where each transformation is chosen to reduce the sum
of the squares of the off-diagonal elements: the matrix uniformly
converges to a diagonal form under successive transformations.  This
convergence can be speeded if one first transforms the matrix into a
\e tridiagonal form (where the only nonzero elements are on the
diagonal or immediately adjacent to it); this transformation
\e can be accomplished deterministically in a finite number of
steps.  This is the approach advocated in Chapter 11
of [\ref ptvf1992], and taken in this module.

The routines in \ref Eigen_c simply call the routines in
\ref EigenInternal_c, first to reduce the matrix to tridiagonal
form, and then to reduce it to diagonal form iteratively.  The
routines <tt>LALSSymmetricEigenVectors()</tt> and
<tt>LALDSymmetricEigenVectors()</tt> call the routines
<tt>LALSSymmetricToTriDiagonal()</tt> and
<tt>LALDSymmetricToTriDiagonal()</tt> to tri-diagonalize the matrix,
then <tt>LALSTriDiagonalToDiagonal()</tt> and
<tt>LALDTriDiagonalToDiagonal()</tt> to diagonalize it.  The routines
<tt>LALSSymmetricEigenValues()</tt> and
<tt>LALDSymmetricEigenValues()</tt> instead call
<tt>LALSSymmetricToTriDiagonal2()</tt>\f$\Rightarrow\f$<tt>LALSTriDiagonalToDiagonal2()</tt>
or
<tt>LALDSymmetricToTriDiagonal2()</tt>\f$\Rightarrow\f$<tt>LALDTriDiagonalToDiagonal2()</tt>,
which discard information about the transformation matrix itself and
thus do not give the eigenvectors; however, they execute much faster.

*/

/*@{*/

/** \see See \ref Eigen_c for documentation */
void
LALSSymmetricEigenVectors( LALStatus *stat, REAL4Vector *values, REAL4Array *matrix )
{
  REAL4Vector *offDiag = NULL; /* off-diagonal line of
                                  tri-diagonalized matrix */

  INITSTATUS(stat);
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


/** \see See \ref Eigen_c for documentation */
void
LALSSymmetricEigenValues( LALStatus *stat, REAL4Vector *values, REAL4Array *matrix )
{
  REAL4Vector *offDiag = NULL; /* off-diagonal line of
                                  tri-diagonalized matrix */

  INITSTATUS(stat);
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


/** \see See \ref Eigen_c for documentation */
void
LALDSymmetricEigenVectors( LALStatus *stat, REAL8Vector *values, REAL8Array *matrix )
{
  REAL8Vector *offDiag = NULL; /* off-diagonal line of
                                  tri-diagonalized matrix */

  INITSTATUS(stat);
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


/** \see See \ref Eigen_c for documentation */
void
LALDSymmetricEigenValues( LALStatus *stat, REAL8Vector *values, REAL8Array *matrix )
{
  REAL8Vector *offDiag = NULL; /* off-diagonal line of
                                  tri-diagonalized matrix */

  INITSTATUS(stat);
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

/*@}*/
