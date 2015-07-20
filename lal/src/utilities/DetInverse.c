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

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/MatrixUtils.h>


/*
 * author Creighton, T. D.
 *
 * Internal routines used to compute matrix determinants and inverses.
 *
 * ### Description ###
 *
 * These functions are called by the routines in \ref DetInverse_c to
 * compute the determinant and inverse of a nondegenerate square matrix
 * <tt>*matrix</tt>.  They are useful routines in their own right, though,
 * so they are made publically available.
 *
 * <tt>LALSLUDecomp()</tt> and <tt>LALDLUDecomp()</tt> replace <tt>*matrix</tt>
 * with an LU decomposition of a <em>row-wise permutation</em> of itself.
 * The output parameter <tt>*indx</tt> stores the permutation, and the
 * output <tt>*sgn</tt> records the sign of the permutation.
 *
 * <tt>LALSLUBackSub()</tt> and <tt>LALDLUBackSub()</tt> take the permuted
 * LU-decomposed matrix returned by the above routine, and
 * back-substitute the vector <tt>*vector</tt> representing \f$\mathsf{v}^a\f$
 * in \eqref{eq_linear_system}, to compute the vector \f$\mathsf{x}^b\f$.
 * This is returned in-place in <tt>*vector</tt>.  The input parameter
 * <tt>*indx</tt> is the list of row permutations returned by the above
 * routines.
 *
 * ### Algorithm ###
 *
 * LU decomposition is performed by Crout's algorithm, described in
 * Sec. 2.3 of \cite ptvf1992; the routines in this module are
 * essentially re-implementations of the Numerical Recipes routines
 * <tt>ludcmp()</tt> and <tt>lubksub()</tt>.  For large \f$N\f$, their operation
 * counts are approximately \f$N^3/3\f$ and \f$N^2\f$, respectively.
 *
 * One difference between <tt>ludcmp()</tt> in \cite ptvf1992 and the
 * routines <tt>LALSLUDecomp()</tt> and <tt>LALDLUDecomp()</tt> in this
 * module is the way in which singular matrices are handled.
 * In \cite ptvf1992 , there is a distinction between between a
 * manifestly singular matrix (where an entire row of the matrix is zero)
 * and a numerically singular matrix (if a diagonal element in the
 * decomposed matrix turns out to be zero).  In the former case, they
 * raise an error signal; in the latter, they replace the offending
 * element with a tiny but nonzero number and continue.  This
 * treatment does not strike the present author as satisfactory.
 *
 * Instead, the routines <tt>LALSLUDecomp()</tt> and <tt>LALDLUDecomp()</tt>
 * will \e always return successfully, even with a singular matrix,
 * but will \e not adjust away a numerical singularity.  Instead,
 * they will signal the presence of the singularity in two ways: First,
 * they will set the permutation sign <tt>*sgn</tt> to zero; second, they
 * will set \e all elements of the <tt>*indx</tt> vector yo zero.  This
 * ensures that routines computing the determinant (whose sign depends on
 * <tt>*sgn</tt>) will correctly give a zero determinant, while the
 * meaningless <tt>*indx</tt> provides a simple sanity check for routines
 * such as <tt>LALSLUBackSub()</tt> and <tt>LALDLUBackSub()</tt> that attempt
 * to invert the linear system.  Note that the returned value of
 * <tt>*matrix</tt> will be meaningless garbage.
 *
 */


static void
LALSLUDecomp( LALStatus   *stat,
	      INT2        *sgn,
	      REAL4Array  *matrix,
	      UINT4Vector *indx )
{
  UINT4 n, imax = 0;    /* matrix dimension and pivot index */
  UINT4 i, j, k;        /* dimension indecies */
  UINT4 ij, ik, kj, jk; /* matrix array indecies */
  UINT4 *d;             /* pointer to indx->data */
  REAL4 *s, *m;         /* pointers to data arrays */
  REAL4 tmp, max, sum;  /* temporary computation variables */

  INITSTATUS(stat);

  /* Check input fields. */
  ASSERT( sgn, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( indx, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( indx->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength->data, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength->length == 2, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  n = indx->length;
  ASSERT( n, stat, MATRIXUTILSH_EDIM, MATRIXUTILSH_MSGEDIM );
  ASSERT( n == matrix->dimLength->data[0], stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  ASSERT( n == matrix->dimLength->data[1], stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );

  /* Set up pointers. */
  if ( !( s = (REAL4 *)LALMalloc( n*sizeof(REAL4) ) ) ) {
    ABORT( stat, MATRIXUTILSH_EMEM, MATRIXUTILSH_MSGEMEM );
  }
  m = matrix->data;
  d = indx->data;
  *sgn = 1;

  /* Get row scaling information. */
  for ( i = 0; i < n; i++ ) {
    max = 0.0;
    for ( j = 0, ij = i*n; j < n; j++, ij++ )
      if ( ( tmp = fabs( m[ij] ) ) > max )
	max = tmp;
    /* Check for singular matrix. */
    if ( max == 0.0 ) {
      *sgn = 0;
      memset( d, 0, n*sizeof(UINT4) );
      LALFree( s );
      RETURN( stat );
    }
    s[i] = 1.0/max;
  }

  /* Loop over columns of matrix. */
  for ( j = 0; j < n; j++ ) {

    /* Compute Uij for i not equal to j. */
    for ( i = 0, ij = j; i < j; i++, ij += n ) {
      sum = m[ij];
      for ( k = 0, ik = i*n, kj = j; k < i; k++, ik++, kj += n )
	sum -= m[ik]*m[kj];
      m[ij] = sum;
    }

    /* Compute Uii and Lij while searching for pivot point. */
    max = 0.0;
    for ( i = j, ij = j*(n+1); i < n; i++, ij += n ) {
      sum = m[ij];
      for ( k = 0, ik = i*n, kj = j; k < j; k++, ik++, kj += n )
	sum -= m[ik]*m[kj];
      m[ij] = sum;
      if ( ( tmp = s[i]*fabs( sum ) ) >= max ) {
	max = tmp;
	imax = i;
      }
    }

    /* Permute rows if necessary. */
    if ( j != imax ) {
      for ( k = 0, ik = imax*n, jk = j*n; k < n; k++, ik++, jk++ ) {
	tmp = m[ik];
	m[ik] = m[jk];
	m[jk] = tmp;
      }
      s[imax] = s[j];
      *sgn *= -1;
    }
    d[j] = imax;

    /* Divide by pivot element, checking for singularity. */
    if ( ( tmp = m[j*(n+1)] ) == 0.0 ) {
      *sgn = 0;
      memset( d, 0, n*sizeof(UINT4) );
      LALFree( s );
      RETURN( stat );
    }
    tmp = 1.0/tmp;
    if ( j != n )
      for ( i = j + 1, ij = j*(n+1) + n; i < n; i++, ij += n )
	m[ij] *= tmp;
  }

  /* Finished! */
  LALFree( s );
  RETURN( stat );
}


static void
LALSLUBackSub( LALStatus   *stat,
	       REAL4Vector *vector,
	       REAL4Array  *matrix,
	       UINT4Vector *indx )
{
  INT4 n;             /* matrix dimension */
  INT4 i, j;          /* dimension indecies */
  INT4 ii, ij;        /* matrix array indecies */
  INT4 id, jmin = -1; /* permuted and first nonzero vector indecies */
  UINT4 *d;           /* pointer to indx->data */
  REAL4 *v, *m;       /* pointers to data arrays */
  REAL4 sum;          /* temporary computation variable */

  INITSTATUS(stat);

  /* Check input fields. */
  ASSERT( vector, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( vector->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( indx, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( indx->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength->data, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength->length == 2, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  n = (INT4)( indx->length );
  ASSERT( n, stat, MATRIXUTILSH_EDIM, MATRIXUTILSH_MSGEDIM );
  ASSERT( n == (INT4)( vector->length ), stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  ASSERT( n == (INT4)( matrix->dimLength->data[0] ), stat,
	  MATRIXUTILSH_EDIM, MATRIXUTILSH_MSGEDIM );
  ASSERT( n == (INT4)( matrix->dimLength->data[1] ), stat,
	  MATRIXUTILSH_EDIM, MATRIXUTILSH_MSGEDIM );

  /* Set up pointers. */
  v = vector->data;
  m = matrix->data;
  d = indx->data;

  /* Check for singular matrix. */
  if ( n > 1 && d[0] == 0 && d[1] == 0 ) {
    ABORT( stat, MATRIXUTILSH_ESING, MATRIXUTILSH_MSGESING );
  }

  /* Compute intermediate vector y, reversing the permutation as we
     go, and skipping quickly to the first nonzero element of v. */
  for ( i = 0; i < n; i++ ) {
    id = d[i];
    sum = v[id];
    v[id] = v[i];
    if ( jmin >= 0 )
      for ( j = jmin, ij = i*n + jmin; j < i; j++, ij++ )
	sum -= m[ij]*v[j];
    else if ( sum != 0.0 )
      jmin = i;
    v[i] = sum;
  }

  /* Now compute the final vector x. */
  for ( i = n - 1, ii = n*n - 1; i >= 0; i--, ii -= n + 1 ) {
    sum = v[i];
    for ( j = i + 1, ij = ii + 1; j < n; j++, ij++ )
      sum -= m[ij]*v[j];
    v[i] = sum/m[ii];
  }
  RETURN( stat );
}


static void
LALDLUDecomp( LALStatus   *stat,
	      INT2        *sgn,
	      REAL8Array  *matrix,
	      UINT4Vector *indx )
{
  UINT4 n, imax = 0;    /* matrix dimension and pivot index */
  UINT4 i, j, k;        /* dimension indecies */
  UINT4 ij, ik, kj, jk; /* matrix array indecies */
  UINT4 *d;             /* pointer to indx->data */
  REAL8 *s, *m;         /* pointers to data arrays */
  REAL8 tmp, max, sum;  /* temporary computation variables */

  INITSTATUS(stat);

  /* Check input fields. */
  ASSERT( sgn, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( indx, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( indx->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength->data, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength->length == 2, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  n = indx->length;
  ASSERT( n, stat, MATRIXUTILSH_EDIM, MATRIXUTILSH_MSGEDIM );
  ASSERT( n == matrix->dimLength->data[0], stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  ASSERT( n == matrix->dimLength->data[1], stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );

  /* Set up pointers. */
  if ( !( s = (REAL8 *)LALMalloc( n*sizeof(REAL8) ) ) ) {
    ABORT( stat, MATRIXUTILSH_EMEM, MATRIXUTILSH_MSGEMEM );
  }
  m = matrix->data;
  d = indx->data;
  *sgn = 1;

  /* Get row scaling information. */
  for ( i = 0; i < n; i++ ) {
    max = 0.0;
    for ( j = 0, ij = i*n; j < n; j++, ij++ )
      if ( ( tmp = fabs( m[ij] ) ) > max )
	max = tmp;
    /* Check for singular matrix. */
    if ( max == 0.0 ) {
      *sgn = 0;
      memset( d, 0, n*sizeof(UINT4) );
      LALFree( s );
      RETURN( stat );
    }
    s[i] = 1.0/max;
  }

  /* Loop over columns of matrix. */
  for ( j = 0; j < n; j++ ) {

    /* Compute Uij for i not equal to j. */
    for ( i = 0, ij = j; i < j; i++, ij += n ) {
      sum = m[ij];
      for ( k = 0, ik = i*n, kj = j; k < i; k++, ik++, kj += n )
	sum -= m[ik]*m[kj];
      m[ij] = sum;
    }

    /* Compute Uii and Lij while searching for pivot point. */
    max = 0.0;
    for ( i = j, ij = j*(n+1); i < n; i++, ij += n ) {
      sum = m[ij];
      for ( k = 0, ik = i*n, kj = j; k < j; k++, ik++, kj += n )
	sum -= m[ik]*m[kj];
      m[ij] = sum;
      if ( ( tmp = s[i]*fabs( sum ) ) >= max ) {
	max = tmp;
	imax = i;
      }
    }

    /* Permute rows if necessary. */
    if ( j != imax ) {
      for ( k = 0, ik = imax*n, jk = j*n; k < n; k++, ik++, jk++ ) {
	tmp = m[ik];
	m[ik] = m[jk];
	m[jk] = tmp;
      }
      s[imax] = s[j];
      *sgn *= -1;
    }
    d[j] = imax;

    /* Divide by pivot element, checking for singularity. */
    if ( ( tmp = m[j*(n+1)] ) == 0.0 ) {
      *sgn = 0;
      memset( d, 0, n*sizeof(UINT4) );
      LALFree( s );
      RETURN( stat );
    }
    tmp = 1.0/tmp;
    if ( j != n )
      for ( i = j + 1, ij = j*(n+1) + n; i < n; i++, ij += n )
	m[ij] *= tmp;
  }

  /* Finished! */
  LALFree( s );
  RETURN( stat );
}


static void
LALDLUBackSub( LALStatus   *stat,
	       REAL8Vector *vector,
	       REAL8Array  *matrix,
	       UINT4Vector *indx )
{
  INT4 n;             /* matrix dimension */
  INT4 i, j;          /* dimension indecies */
  INT4 ii, ij;        /* matrix array indecies */
  INT4 id, jmin = -1; /* permuted and first nonzero vector indecies */
  UINT4 *d;           /* pointer to indx->data */
  REAL8 *v, *m;       /* pointers to data arrays */
  REAL8 sum;          /* temporary computation variable */

  INITSTATUS(stat);

  /* Check input fields. */
  ASSERT( vector, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( vector->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( indx, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( indx->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength->data, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength->length == 2, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  n = (INT4)( indx->length );
  ASSERT( n, stat, MATRIXUTILSH_EDIM, MATRIXUTILSH_MSGEDIM );
  ASSERT( n == (INT4)( vector->length ), stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  ASSERT( n == (INT4)( matrix->dimLength->data[0] ), stat,
	  MATRIXUTILSH_EDIM, MATRIXUTILSH_MSGEDIM );
  ASSERT( n == (INT4)( matrix->dimLength->data[1] ), stat,
	  MATRIXUTILSH_EDIM, MATRIXUTILSH_MSGEDIM );

  /* Set up pointers. */
  v = vector->data;
  m = matrix->data;
  d = indx->data;

  /* Check for singular matrix. */
  if ( n > 1 && d[0] == 0 && d[1] == 0 ) {
    ABORT( stat, MATRIXUTILSH_ESING, MATRIXUTILSH_MSGESING );
  }

  /* Compute intermediate vector y, reversing the permutation as we
     go, and skipping quickly to the first nonzero element of v. */
  for ( i = 0; i < n; i++ ) {
    id = d[i];
    sum = v[id];
    v[id] = v[i];
    if ( jmin >= 0 )
      for ( j = jmin, ij = i*n + jmin; j < i; j++, ij++ )
	sum -= m[ij]*v[j];
    else if ( sum != 0.0 )
      jmin = i;
    v[i] = sum;
  }

  /* Now compute the final vector x. */
  for ( i = n - 1, ii = n*n - 1; i >= 0; i--, ii -= n + 1 ) {
    sum = v[i];
    for ( j = i + 1, ij = ii + 1; j < n; j++, ij++ )
      sum -= m[ij]*v[j];
    v[i] = sum/m[ii];
  }
  RETURN( stat );
}


/**
 * \defgroup DetInverse_c Module DetInverse.c
 * \ingroup MatrixUtils_h
 * \author Creighton, T. D.
 *
 * \brief Routines to compute matrix determinants and inverses.
 *
 * ### Description ###
 *
 * <tt>LALDMatrixDeterminant()</tt>
 * computes the determinant <tt>*det</tt> of the square matrix
 * <tt>*matrix</tt>.  The internal computations will corrupt
 * <tt>*matrix</tt>; if you don't want the input matrix to be changed, make
 * a copy of it first.
 *
 * <tt>LALSMatrixInverse()</tt> and <tt>LALDMatrixInverse()</tt> compute the
 * inverse <tt>*inverse</tt> of the square matrix <tt>*matrix</tt>.  If the
 * pointer \c det is non-\c NULL, then the determinant is also
 * computed and returned in <tt>*det</tt>.  The array <tt>*inverse</tt> must
 * already be allocated to the correct size.  The internal computations
 * will corrupt <tt>*matrix</tt>; if you don't want the input matrix to be
 * changed, make a copy of it first.
 *
 * <tt>LALDMatrixDeterminantErr()</tt> computes the determinant
 * <tt>det[0]</tt> of the square matrix <tt>*matrix</tt>.  If
 * <tt>*matrixErr</tt> is non-\c NULL, it stores uncertainties in the
 * corresponding components of <tt>*matrix</tt>, and the resulting
 * uncertainty in the determinant (computed by linear error propagation)
 * is returned in <tt>det[1]</tt>.  This routine is not highly optimized,
 * and uses an internal matrix in its computations, so the contents of
 * <tt>*matrix</tt> and <tt>*matrixErr</tt> will not be changed.
 *
 * ### Algorithm ###
 *
 * A linear system of equations is written in matrix terms as:
 * \f{equation}{
 * \label{eq_linear_system}
 * \mathsf{M}^a{}_b \mathsf{x}^b = \mathsf{v}^a \;,
 * \f}
 * where \f$\mathsf{M}^a{}_b\f$ is a known matrix, \f$\mathsf{v}^a\f$ a known
 * vector, and \f$\mathsf{x}^b\f$ an unknown vector that we wish to solve
 * for.  A standard method for solving this is the method of LU
 * decomposition, based on the fact that any non-singular square matrix
 * \f$\mathsf{M}^a{}_b\f$ can be \e decomposed into the product of a
 * lower-triangular matrix \f$\mathsf{L}^a{}_b\f$ and an upper-triangular
 * matrix \f$\mathsf{U}^a{}_b\f$:
 * \f{equation}{
 * \mathsf{M}^a{}_b = \mathsf{L}^a{}_c \mathsf{U}^c{}_b =
 * \left[\begin{array}{cccc}
 * 1    &    0    & \cdots &    0    \\
 * L^2{}_1 &    1    & \cdots &    0    \\
 * \vdots  & \vdots  & \ddots & \vdots  \\
 * L^N{}_1 & L^N{}_2 & \cdots &    1
 * \end{array}\right] \cdot \left[\begin{array}{cccc}
 * U^1{}_1 & U^1{}_2 & \cdots & U^1{}_N \\
 * 0    & U^2{}_2 & \cdots & U^2{}_N \\
 * \vdots  & \vdots  & \ddots & \vdots  \\
 * 0    &    0    & \cdots & U^N{}_N
 * \end{array}\right] \;.
 * \f}
 * The linear system can then be broken up as \f$\mathsf{L}^a{}_b
 * \mathsf{y}^b = \mathsf{v}^a\f$ and \f$\mathsf{U}^b{}_c \mathsf{x}^c =
 * \mathsf{y}^b\f$, where these two equations are trivially invertible:
 * \f{eqnarray}{
 * y^i & = & v^i - \sum_{j=1}^i-1 L^i{}_j y^j \;,\qquad i=1,2,\ldots,N \;, \\
 * x^i & = & \frac{1}{U^i{}_i}\left( y^i - \sum_{j=i+1}^N U^i{}_j x_j \right)
 * \;,\qquad i=N,N-1,\ldots,1 \;,
 * \f}
 * where the calculations are arranged so that the computation of a given
 * \f$y^i\f$ or \f$x^i\f$ requires only those values that have been computed
 * already.  This process of solving the linear system is called
 * \e backsubstitution.
 *
 * The determinant of \f$\mathsf{M}^a{}_b\f$ is simply the product of the
 * diagonal elements of the decomposed matrix:
 * \f$|\mathsf{M}^a{}_b|=\prod_{i=1}^N U^i{}_i\f$.  The inverse matrix
 * \f$(\mathsf{M}^{-1}){}^a{}_b\f$ can be computed by performing a
 * column-by-column backsubstitution of the identity matrix.
 *
 * The routines in \ref DetInverse_c first use <tt>LALSLUDecomp()</tt> or
 * <tt>LALDLUDecomp()</tt> to perform an LU decomposition of the matrix,
 * then either compute the determinant from the diagonal elements, or
 * use <tt>LALSLUBackSub()</tt> or <tt>LALDLUBackSub()</tt> to determine
 * the inverse by back-substitution of basis vectors.  The routines that
 * compute the determinant will also handle any singular matrix error
 * code returned by the LU decomposition routines, returning zero as the
 * determinant.
 *
 * Since the diagonal components \f$L^i{}_i\f$ are conventionally assigned to
 * 1, they do not need to be explicitly stored.  Therefore we can store
 * \e both matrices \f$\mathsf{L}^a{}_b\f$ and \f$\mathsf{U}^a{}_b\f$
 * in-place, in the same memory block used for the input matrix
 * \f$\mathsf{M}^a{}_b\f$.  This the procedure taken by <tt>LALSLUDecomp()</tt>
 * and <tt>LALDLUDecomp()</tt>; hence on return the routines in this module
 * will leave the input <tt>*matrix</tt> in this decomposed state.
 * However, these routines also <em>permute the rows</em> of the input
 * matrix, and the information on this permutation is \e not returned
 * by the routines in \ref DetInverse_c, so the information in
 * <tt>*matrix</tt> will be irretrievably mangled.
 *
 * Computing the determinant is dominated by the cost of doing the LU
 * decomposition, or of order \f$N^3/3\f$ operations.  Computing the inverse
 * requires an additional \f$N\f$ back-substitutions, but since most of the
 * vector elements are zero, this reduces the average cost of each
 * back-substitution from \f$N^2\f$ to \f$2N^2/3\f$, for a total operation count
 * of \f$N^3\f$.
 *
 * \par Computing determinant uncertainties:
 * To determine the
 * dependence of the determinant on any one element of the matrix, we
 * take advantage of the fact that the determinant \f$|\mathsf{M}^a{}_b|\f$
 * can be written as:
 * \f{eqnarray}{
 * |\mathsf{M}^a{}_b| & = & \sum_{j=1}^N (-1)^{i+j} M^i{}_j
 * |(\mathsf{C}^i{}_j){}^a{}_b| \quad\mbox{for any }i\in[1,N]
 * \;,\\
 * & = & \sum_{i=1}^N (-1)^{i+j} M^i{}_j
 * |(\mathsf{C}^i{}_j){}^a{}_b| \quad\mbox{for any }j\in[1,N] \;,
 * \f}
 * where
 * \f{equation}{
 * (\mathsf{C}^i{}_j){}^a{}_b = \left[\begin{array}{ccccccc}
 * M^1{}_1   & \cdots &   M^1{}_{j-1}   &    0   &   M^1{}_{j+1}
 * & \cdots &   M^1{}_N   \\
 * \vdots   & \ddots &      \vdots     & \vdots &      \vdots
 * &        &    \vdots   \\
 * M^{i-1}{}_1 & \cdots & M^{i-1}{}_{j-1} &    0   & M^{i-1}{}_{j+1}
 * & \cdots & M^{i-1}{}_N \\
 * 0     & \cdots &        0        &    1   &       0
 * & \cdots &      0      \\
 * \vdots   &        &      \vdots     & \vdots &      \vdots
 * & \ddots &    \vdots   \\
 * M^N{}_1   & \cdots &   M^N{}_{j-1}   &    0   &   M^N{}_{j+1}
 * & \cdots &   M^N{}_N \end{array}\right]
 * \f}
 * is called the co-matrix of the element \f$M^i{}_j\f$.
 *
 * Assuming all matrix elements are statistically independent, linear
 * error propagation can give us the uncertainty in the determinant:
 * \f{eqnarray}{
 * E\left(|\mathsf{M}^a{}_b|\right) & = & \sqrt{\sum_{i,j=1}^N
 * \left[ \frac{\partial|\mathsf{M}^a{}_b|}{\partial M^i{}_j}
 * E\left(M^i{}_j\right)\right]^2 }\\
 * & = & \sqrt{\sum_{i,j=1}^N \left[ |(\mathsf{C}^i{}_j){}^a{}_b|
 * E\left(M^i{}_j\right)\right]^2 } \;.\\
 * \f}
 * The routines <tt>LALSMatrixDeterminantErr()</tt> and
 * <tt>LALDMatrixDeterminantErr()</tt> thus simply compute the determinant
 * multiple times, successively zeroing out rows and columns of
 * <tt>*matrix</tt> and replacing the element at their intersection with
 * the corresponding element of <tt>*matrixErr</tt>, and take the square
 * root of the sum of the squares of these determinants.  As mentioned
 * earlier, they use an internal matrix for all computations, so the
 * input parameters are unchanged.  The uncertainty requires evaluating
 * \f$N^2\f$ determinants, so the operation count scales as \f$N^5\f$: this
 * routine should \e not be used for large matrices!
 *
 */
/*@{*/

/** \see See \ref DetInverse_c for documentation */
void
LALSMatrixInverse( LALStatus *stat, REAL4 *det, REAL4Array *matrix, REAL4Array *inverse )
{
  INT2 sgn;                 /* sign of permutation. */
  UINT4 n, i, j, ij;        /* array dimension and indecies */
  UINT4Vector *indx = NULL; /* permutation storage vector */
  REAL4Vector *col = NULL;  /* columns of inverse matrix */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* Check dimension length.  Other argument testing is done by the
     subroutines. */
  ASSERT( matrix, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength->data, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( inverse, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( inverse->dimLength, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( inverse->dimLength->data, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( inverse->dimLength->length == 2, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  n = inverse->dimLength->data[0];
  ASSERT( n, stat, MATRIXUTILSH_EDIM, MATRIXUTILSH_MSGEDIM );
  ASSERT( n == inverse->dimLength->data[1], stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );

  /* Allocate vectors to store the matrix permutation and column
     vectors. */
  TRY( LALU4CreateVector( stat->statusPtr, &indx, n ), stat );
  LALSCreateVector( stat->statusPtr, &col, n );
  BEGINFAIL( stat ) {
    TRY( LALU4DestroyVector( stat->statusPtr, &indx ), stat );
  } ENDFAIL( stat );

  /* Decompose the matrix. */
  sgn = 0;	/* silence -Werror=maybe-uninitialized.  note:  there is no check here that LALSLUDecomp() is successful and the compiler knows this.  *that's* the real problem, but I'm not fixing (Kipp) */
  LALSLUDecomp( stat->statusPtr, &sgn, matrix, indx );
  BEGINFAIL( stat ) {
    TRY( LALU4DestroyVector( stat->statusPtr, &indx ), stat );
    TRY( LALSDestroyVector( stat->statusPtr, &col ), stat );
  } ENDFAIL( stat );

  /* Compute the determinant, if available. */
  if ( det ) {
    *det = sgn;
    for ( i = 0; i < n; i++ )
      *det *= matrix->data[i*(n+1)];
  }

  /* Compute each column of inverse matrix. */
  for ( j = 0; j < n; j++ ) {
    memset( col->data, 0, n*sizeof(REAL4) );
    col->data[j] = 1.0;
    LALSLUBackSub( stat->statusPtr, col, matrix, indx );
    BEGINFAIL( stat ) {
      TRY( LALU4DestroyVector( stat->statusPtr, &indx ), stat );
      TRY( LALSDestroyVector( stat->statusPtr, &col ), stat );
    } ENDFAIL( stat );
    for ( i = 0, ij = j; i < n; i++, ij += n )
      inverse->data[ij] = col->data[i];
  }

  /* Clean up. */
  TRY( LALU4DestroyVector( stat->statusPtr, &indx ), stat );
  TRY( LALSDestroyVector( stat->statusPtr, &col ), stat );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/** \see See \ref DetInverse_c for documentation */
void
LALDMatrixDeterminant( LALStatus *stat, REAL8 *det, REAL8Array *matrix )
{
  INT2 sgn;                 /* sign of permutation. */
  UINT4 n, i;               /* array dimension and index */
  UINT4Vector *indx = NULL; /* permutation storage vector */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* Check dimension length.  All other argument testing is done by
     the subroutines. */
  ASSERT( det, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength->data, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  n = matrix->dimLength->data[0];
  ASSERT( n, stat, MATRIXUTILSH_EDIM, MATRIXUTILSH_MSGEDIM );

  /* Allocate a vector to store the matrix permutation. */
  TRY( LALU4CreateVector( stat->statusPtr, &indx, n ), stat );

  /* Decompose the matrix. */
  LALDLUDecomp( stat->statusPtr, &sgn, matrix, indx );
  BEGINFAIL( stat ) {
    TRY( LALU4DestroyVector( stat->statusPtr, &indx ), stat );
  } ENDFAIL( stat );

  /* Compute determinant. */
  *det = sgn;
  for ( i = 0; i < n; i++ )
    *det *= matrix->data[i*(n+1)];

  /* Clean up. */
  TRY( LALU4DestroyVector( stat->statusPtr, &indx ), stat );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/** \see See \ref DetInverse_c for documentation */
void
LALDMatrixInverse( LALStatus *stat, REAL8 *det, REAL8Array *matrix, REAL8Array *inverse )
{
  INT2 sgn;                 /* sign of permutation. */
  UINT4 n, i, j, ij;        /* array dimension and indecies */
  UINT4Vector *indx = NULL; /* permutation storage vector */
  REAL8Vector *col = NULL;  /* columns of inverse matrix */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* Check dimension length.  Other argument testing is done by the
     subroutines. */
  ASSERT( matrix, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength->data, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( inverse, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( inverse->dimLength, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( inverse->dimLength->data, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( inverse->dimLength->length == 2, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  n = inverse->dimLength->data[0];
  ASSERT( n, stat, MATRIXUTILSH_EDIM, MATRIXUTILSH_MSGEDIM );
  ASSERT( n == inverse->dimLength->data[1], stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );

  /* Allocate vectors to store the matrix permutation and column
     vectors. */
  TRY( LALU4CreateVector( stat->statusPtr, &indx, n ), stat );
  LALDCreateVector( stat->statusPtr, &col, n );
  BEGINFAIL( stat ) {
    TRY( LALU4DestroyVector( stat->statusPtr, &indx ), stat );
  } ENDFAIL( stat );

  /* Decompose the matrix. */
  LALDLUDecomp( stat->statusPtr, &sgn, matrix, indx );
  BEGINFAIL( stat ) {
    TRY( LALU4DestroyVector( stat->statusPtr, &indx ), stat );
    TRY( LALDDestroyVector( stat->statusPtr, &col ), stat );
  } ENDFAIL( stat );

  /* Compute the determinant, if available. */
  if ( det ) {
    *det = sgn;
    for ( i = 0; i < n; i++ )
      *det *= matrix->data[i*(n+1)];
  }

  /* Compute each column of inverse matrix. */
  for ( j = 0; j < n; j++ ) {
    memset( col->data, 0, n*sizeof(REAL8) );
    col->data[j] = 1.0;
    LALDLUBackSub( stat->statusPtr, col, matrix, indx );
    BEGINFAIL( stat ) {
      TRY( LALU4DestroyVector( stat->statusPtr, &indx ), stat );
      TRY( LALDDestroyVector( stat->statusPtr, &col ), stat );
    } ENDFAIL( stat );
    for ( i = 0, ij = j; i < n; i++, ij += n )
      inverse->data[ij] = col->data[i];
  }

  /* Clean up. */
  TRY( LALU4DestroyVector( stat->statusPtr, &indx ), stat );
  TRY( LALDDestroyVector( stat->statusPtr, &col ), stat );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/** \see See \ref DetInverse_c for documentation */
void
LALDMatrixDeterminantErr( LALStatus *stat, REAL8 det[2], REAL8Array *matrix, REAL8Array *matrixErr )
{
  UINT4 n, i, j, k;        /* matrix dimension and indecies */
  UINT4 ij, ik, kj;        /* array data indecies */
  REAL8 var;               /* partial derivative of determinant */
  REAL8Array *temp = NULL; /* internal copy of *matrix */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* Check input arguments. */
  ASSERT( det, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->data, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength->data, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength->length == 2, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  n = matrix->dimLength->data[0];
  ASSERT( matrix->dimLength->data[1] == n, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  if ( matrixErr ) {
    ASSERT( matrixErr->data, stat, MATRIXUTILSH_ENUL,
	    MATRIXUTILSH_MSGENUL );
    ASSERT( matrixErr->dimLength, stat, MATRIXUTILSH_ENUL,
	    MATRIXUTILSH_MSGENUL );
    ASSERT( matrixErr->dimLength->data, stat, MATRIXUTILSH_ENUL,
	    MATRIXUTILSH_MSGENUL );
    ASSERT( matrixErr->dimLength->length == 2, stat, MATRIXUTILSH_EDIM,
	    MATRIXUTILSH_MSGEDIM );
    ASSERT( matrixErr->dimLength->data[0] == n, stat, MATRIXUTILSH_EDIM,
	    MATRIXUTILSH_MSGEDIM );
    ASSERT( matrixErr->dimLength->data[1] == n, stat, MATRIXUTILSH_EDIM,
	    MATRIXUTILSH_MSGEDIM );
  }

  /* Allocate internal copy, and compute determinant. */
  TRY( LALDCreateArray( stat->statusPtr, &temp, matrix->dimLength ),
       stat );
  memcpy( temp->data, matrix->data, n*n*sizeof(REAL8) );
  LALDMatrixDeterminant( stat->statusPtr, det, temp );
  BEGINFAIL( stat ) {
    TRY( LALDDestroyArray( stat->statusPtr, &temp ), stat );
  } ENDFAIL( stat );
  det[1] = 0.0;

  /* Compute derivatives of determinant. */
  if ( matrixErr ) {
    for ( i = 0; i < n; i++ ) {
      for ( j = 0, ij = i*n; j < n; j++, ij++ ) {
	memcpy( temp->data, matrix->data, n*n*sizeof(REAL8) );
	for ( k = 0, ik = i*n, kj = j; k < n; k++, ik++, kj += n )
	  temp->data[ik] = temp->data[kj] = 0.0;
	temp->data[ij] = matrixErr->data[ij];
	LALDMatrixDeterminant( stat->statusPtr, &var, temp );
	BEGINFAIL( stat ) {
	  TRY( LALDDestroyArray( stat->statusPtr, &temp ), stat );
	} ENDFAIL( stat );
	det[1] += var*var;
      }
    }
    det[1] = sqrt( det[1] );
  }

  /* Done. */
  TRY( LALDDestroyArray( stat->statusPtr, &temp ), stat );
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
/*@}*/
