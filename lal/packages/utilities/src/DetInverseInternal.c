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
#include <lal/LALConstants.h>
#include <lal/MatrixUtils.h>



/**
 * \defgroup \DetInverseInternal_c Module DetInverseInternal.c
 * \ingroup MatrixUtils_h
 * \author Creighton, T. D.
 *
 * \brief Internal routines used to compute matrix determinants and inverses.
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
 * In \cite ptvf1992, there is a distinction between between a
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
/*@{*/


/** \see See \ref DetInverseInternal_c for documentation */
void
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


/** \see See \ref DetInverseInternal_c for documentation */
void
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


/** \see See \ref DetInverseInternal_c for documentation */
void
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


/** \see See \ref DetInverseInternal_c for documentation */
void
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
/*@}*/
