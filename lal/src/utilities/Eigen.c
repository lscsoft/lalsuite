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
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/MatrixUtils.h>

#define EIGENINTERNALC_MAXITER 30 /* max. number of iterations in
				     diagonalizing */

/* A quick macro to take the square root of the sum of squares,
   without overflowing or underflowing. */
#define RSS( a, b ) ( fabs( a ) > fabs( b ) ?                        \
fabs( a )*sqrt( 1.0 + ((b)/(a))*((b)/(a)) ) :                        \
( (b) == 0.0 ? 0.0 : fabs( b )*sqrt( 1.0 + ((a)/(b))*((a)/(b)) ) ) )


/*
 * \defgroup EigenInternal_c Module EigenInternal.c
 * \ingroup MatrixUtils_h
 * \author Creighton, T. D.
 *
 * \brief Internal routines used to compute eigenvalues and eigenvectors.
 *
 * ### Description ###
 *
 * These functions are called by the routines in \ref Eigen_c to
 * compute eigenvalues and eigenvectors of a symmetric square matrix
 * <tt>*matrix</tt>.  They are provided because, in some circumstances,
 * users may find it useful to work with the partially-reduced results.
 *
 * <tt>LALSSymmetricToTriDiagonal()</tt> and
 * <tt>LALDSymmetricToTriDiagonal()</tt> reduce the symmetric square matrix
 * <tt>*matrix</tt> to tridiagonal form.  The vectors <tt>*diag</tt> and
 * <tt>*offDiag</tt> must be allocated to the same length as each matrix
 * dimension; on return, <tt>*diag</tt> stores the diagonal elements and
 * <tt>*offDiag</tt> the off-diagonal elements (<tt>offDiag-\>data[0]</tt> is
 * meaningless and will be set to zero).  The tri-diagonalization is done
 * in-place, so that on return <tt>*matrix</tt> will store the
 * transformation matrix that brings the original input matrix into
 * tridiagonal form.  If you don't want the input matrix to be changed,
 * make a copy of it first.
 *
 * <tt>LALSTriDiagonalToDiagonal()</tt> and
 * <tt>LALDTriDiagonalToDiagonal()</tt> take a symmetric tridiagonal matrix
 * and compute its eigenvalues and eigenvectors.  On input, <tt>*diag</tt>
 * stores the diagonal elements of the matrix, and <tt>*offDiag</tt> the
 * off-diagonal elements, with <tt>offDiag-\>data[0]</tt> arbitrary; on
 * return, <tt>*diag</tt> will store the eigenvalues and <tt>*offDiag</tt>
 * will store zeroes.  The matrix <tt>*matrix</tt> should store the
 * orthogonal transformation matrix that brought the original symmetric
 * matrix into tri-diagonal form (as returned by the above routines), or
 * the identity matrix if the tri-diagonal matrix \e is the original
 * matrix of interest; on return, <tt>*matrix</tt> will store the
 * orthogonal transformation matrix that diagonalizes the original
 * matrix: its columns are the eigenvectors of the original matrix.
 *
 * ### Algorithm ###
 *
 * The tri-diagonalizing routines follow the Householder reduction method
 * described in Sec. 11.2 of \cite ptvf1992; they are essentially
 * re-implementations of the Numerical Recipes routine <tt>tred2()</tt>.
 * For large \f$N\f$, their operation count is approximately \f$4N^3/3\f$, or
 * \f$2N^3/3\f$ for the routines that ignore eigenvectors.  These routines
 * explicitly enforce symmetry, by only using the lower-left triangle of
 * the matrix as input.
 *
 * The diagonalizing routines follow the QL algorithm with implicit
 * shifts described in Sec. 11.3 of \cite ptvf1992; they are
 * essentially re-implementations of the Numerical Recipes routine
 * <tt>tqli()</tt>.  Depending on the number of iterations required, their
 * operation count is roughly \f$\sim30N^2\f$, plus \f$\sim3N^3\f$ if
 * eigenvectors are also being computed.
 *
 * The diagonalizing routines can fail if they fail to converge rapidly
 * enough.  The discussion in \cite ptvf1992 does not go into much
 * detail about when this is likely to occur, except to note that
 * degenerate eigenvalues converge more slowly.  If the routines fail in
 * this way, \c diag, \c matrix, and \c offDiag will all be
 * left in indeterminate states.
 *
 */

static void
LALSSymmetricToTriDiagonal( LALStatus   *stat,
			    REAL4Vector *diag,
			    REAL4Array  *matrix,
			    REAL4Vector *offDiag )
{
  UINT4 n, in;                          /* matrix dimension (times i) */
  UINT4 i, k, l;                        /* dimension indecies */
  UINT4 ij, ik, ki, il, li, kl, lk, ii; /* matrix array indecies */
  REAL4 *m, *d, *o;                     /* pointers to data arrays */
  REAL4 H, K; /* scalars defined in Eqs. (11.2.4) and (11.2.11) of NR */
  REAL4 x, y, scale; /* temporary variables */

  INITSTATUS(stat);

  /* Check input fields. */
  ASSERT( diag, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( diag->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( offDiag, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( offDiag->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength->data, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength->length == 2, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  n = diag->length;
  ASSERT( n, stat, MATRIXUTILSH_EDIM, MATRIXUTILSH_MSGEDIM );
  ASSERT( n == offDiag->length, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  ASSERT( n == matrix->dimLength->data[0], stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  ASSERT( n == matrix->dimLength->data[1], stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );

  /* Set up pointers. */
  m = matrix->data;
  d = diag->data;
  o = offDiag->data;

  /* Start tri-diagonalizing from the bottom-right corner i = n-1. */
  for ( i = n - 1; i; i-- ) {
    /* j = i - 1; */
    in = i*n;
    ij = in + i - 1;
    H = scale = 0.0;
    if ( i > 1 ) {

      /* Compute typical magnitude of off-tri-diagonal components, and
         skip transformation if they are zero. */
      for ( k = 0, ik = in; k < i; k++, ik++ )
	scale += fabs( m[ik] );
      if ( scale == 0.0 )
	o[i] = m[ij];

      /* Otherwise, apply transformation to row i. */
      else {

	/* Renormalize off-tri-diagonal components and compute H. */
	for ( k = 0, ik = in; k < i; k++, ik++ ) {
	  m[ik] /= scale;
	  H += m[ik]*m[ik];
	}
	x = m[ij];
	y = ( x >= 0 ? -sqrt( H ) : sqrt( H ) );
	o[i] = scale*y;
	H -= x*y;

	/* Now compute K.  To start, store u in ith row of matrix. */
	m[ij] = x - y;
	x = 0.0;
	for ( k = 0, ik = in, ki = i; k < i; k++, ik++, ki += n ) {
	  /* Store u/H in ith column of matrix. */
	  m[ki] = m[ik] / H;
	  /* Compute component of matrix * u, stored in y. */
	  y = 0.0;
	  for ( l = 0, kl = k*n, il = in; l <= k; l++, kl++, il++ )
	    y += m[kl]*m[il];
	  for ( l = k + 1, lk = k*( n + 1 ) + n, il = in + k + 1;
		l < i; l++, lk += n, il++ )
	    y += m[lk]*m[il];
	  /* Compute p, stored in unfilled part of o. */
	  o[k] = y/H;
	  x += o[k]*m[ik];
	}
	K = 0.5*x/H;

	/* Now reduce submatrix to tri-diagonal form.  Start by
           forming vector q and storing in o, overwriting p. */
	for ( k = 0, ik = in; k < i; k++, ik++ ) {
	  x = m[ik];
	  o[k] = y = o[k] - K*x;
	  /* Now apply submatrix reduction. */
	  for ( l = 0, kl = n*k, il = in; l <= k; l++, kl++, il++ )
	    m[kl] -= x*o[l] + y*m[il];
	}
      }
    }

    /* Last (top-right) 1x1 submatrix is trivial. */
    else
      o[i] = m[ij];

    /* Indicate whether this iteration was trivial or not. */
    d[i] = H;
  }

  /* Now compute transformation matrix (eigenvectors). */
  d[0] = 0.0;
  o[0] = 0.0;
  /* Loop over the n sub-tranformations. */
  for ( i = 0; i < n; i++ ) {
    in = i*n;
    ii = in + i;
    /* Skip first block and any other trivial transformations. */
    if ( d[i] ) {
      for ( k = 0; k < i; k++ ) {
	y = 0.0;
	/* Use u and u/H, stored in ith row and column of matrix, to
           compute P*Q. */
	for ( l = 0, il = in, lk = k; l < i; l++, il++, lk += n )
	  y += m[il]*m[lk];
	for ( l = 0, li = i, lk = k; l < i; l++, li += n, lk += n )
	  m[lk] -= y*m[li];
      }
    }

    /* Store diagonal element. */
    d[i] = m[ii];

    /* Reset row and column to identity for next iteration. */
    m[ii] = 1.0;
    for ( k = 0, ik = in, ki = i; k < i; k++, ik++, ki += n )
      m[ki] = m[ik] = 0.0;
  }

  /* Finished. */
  RETURN( stat );
}


static void
LALSTriDiagonalToDiagonal( LALStatus   *stat,
			   REAL4Vector *diag,
			   REAL4Array  *matrix,
			   REAL4Vector *offDiag )
{
  INT4 n, i, j, k, l;              /* matrix dimension and indecies */
  INT4 lk;                         /* matrix array index */
  UINT4 iteration;                 /* number of iterated transforms */
  REAL4 *m, *d, *o;                /* pointers to data arrays */
  REAL4 x, y, xx, yy, s, c, r, d2; /* temporary variables */

  INITSTATUS(stat);

  /* Check input fields. */
  ASSERT( diag, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( diag->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( offDiag, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( offDiag->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength->data, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength->length == 2, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  n = (INT4)( diag->length );
  ASSERT( n, stat, MATRIXUTILSH_EDIM, MATRIXUTILSH_MSGEDIM );
  ASSERT( n == (INT4)( offDiag->length ), stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  ASSERT( n == (INT4)( matrix->dimLength->data[0] ), stat,
	  MATRIXUTILSH_EDIM, MATRIXUTILSH_MSGEDIM );
  ASSERT( n == (INT4)( matrix->dimLength->data[1] ), stat,
	  MATRIXUTILSH_EDIM, MATRIXUTILSH_MSGEDIM );

  /* Set up pointers. */
  m = matrix->data;
  d = diag->data;
  o = offDiag->data;

  /* As in Numerical Recipes, we start by shifting the indexing of the
     off-diagonal vector. */
  for ( i = 1; i < n; i++ )
    o[i-1] = o[i];
  o[n-1] = 0.0;

  /* Loop over the i eigenvalues. */
  for ( i = 0; i < n; i++ ) {

    /* For each eigenvalue, loop until we've diagonalized that part of
       the matrix. */
    iteration = 0;
    do {

      /* Look for a small (effectively zero) off-diagonal element. */
      for ( j = i; j < n - 1; j++ ) {
	d2 = fabs( d[j] ) + fabs( d[j+1] );
	if ( fabs( o[j]/d2 ) < LAL_REAL4_EPS )
	  break;
      }

      if ( j != i ) {
	if ( iteration++ == EIGENINTERNALC_MAXITER ) {
	  ABORT( stat, MATRIXUTILSH_EITER, MATRIXUTILSH_MSGEITER );
	}
	x = 0.5*( d[i+1] - d[i] )/o[i];
	r = RSS( x, 1.0 );
	if ( x > 0 )
	  x = d[j] - d[i] + o[i]/( x + r );
	else
	  x = d[j] - d[i] + o[i]/( x - r );
	s = c = 1.0;
	yy = 0.0;

	/* Perform a plane rotation to diagonalize the submatrix, and
           Givens rotations to restore it to tridiagonal form. */
	for ( k = j - 1; k >= i; k-- ) {
	  xx = c*o[k];
	  y = s*o[k];
	  o[k+1] = ( r = RSS( x, y ) );

	  /* Check for underflow. */
	  if ( r <= LAL_REAL4_MIN ) {
	    d[k+1] -= yy;
	    o[j] = 0.0;
	    break;
	  }

	  /* Compute diagonal element. */
	  s = y/r;
	  c = x/r;
	  x = d[k+1] - yy;
	  r = ( d[k] - x )*s + 2.0*c*xx;
	  d[k+1] = x + ( yy = s*r );
	  x = c*r - xx;

	  /* Accumulate transformation matrix. */
	  for ( l = 0, lk = k; l < n; l++, lk += n ) {
	    y = m[lk+1];
	    m[lk+1] = s*m[lk] + c*y;
	    m[lk] = c*m[lk] - s*y;
	  }
	}

	/* Check for underflow. */
	if ( r <= LAL_REAL4_MIN && k >= i )
	  continue;

	/* Adjust diagonal and off-diagonal elements, and continue
           loop. */
	d[i] -= yy;
	o[i] = x;
	o[j] = 0.0;
      }
    } while ( j != i );
  }

  /* Done. */
  RETURN( stat );
}


static void
LALDSymmetricToTriDiagonal( LALStatus   *stat,
			    REAL8Vector *diag,
			    REAL8Array  *matrix,
			    REAL8Vector *offDiag )
{
  UINT4 n, in;                          /* matrix dimension (times i) */
  UINT4 i, k, l;                        /* dimension indecies */
  UINT4 ij, ik, ki, il, li, kl, lk, ii; /* matrix array indecies */
  REAL8 *m, *d, *o;                     /* pointers to data arrays */
  REAL8 H, K; /* scalars defined in Eqs. (11.2.4) and (11.2.11) of NR */
  REAL8 x, y, scale; /* temporary variables */

  INITSTATUS(stat);

  /* Check input fields. */
  ASSERT( diag, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( diag->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( offDiag, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( offDiag->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength->data, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength->length == 2, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  n = diag->length;
  ASSERT( n, stat, MATRIXUTILSH_EDIM, MATRIXUTILSH_MSGEDIM );
  ASSERT( n == offDiag->length, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  ASSERT( n == matrix->dimLength->data[0], stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  ASSERT( n == matrix->dimLength->data[1], stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );

  /* Set up pointers. */
  m = matrix->data;
  d = diag->data;
  o = offDiag->data;

  /* Start tri-diagonalizing from the bottom-right corner i = n-1. */
  for ( i = n - 1; i; i-- ) {
    /* j = i - 1; */
    in = i*n;
    ij = in + i - 1;
    H = scale = 0.0;
    if ( i > 1 ) {

      /* Compute typical magnitude of off-tri-diagonal components, and
         skip transformation if they are zero. */
      for ( k = 0, ik = in; k < i; k++, ik++ )
	scale += fabs( m[ik] );
      if ( scale == 0.0 )
	o[i] = m[ij];

      /* Otherwise, apply transformation to row i. */
      else {

	/* Renormalize off-tri-diagonal components and compute H. */
	for ( k = 0, ik = in; k < i; k++, ik++ ) {
	  m[ik] /= scale;
	  H += m[ik]*m[ik];
	}
	x = m[ij];
	y = ( x >= 0 ? -sqrt( H ) : sqrt( H ) );
	o[i] = scale*y;
	H -= x*y;

	/* Now compute K.  To start, store u in ith row of matrix. */
	m[ij] = x - y;
	x = 0.0;
	for ( k = 0, ik = in, ki = i; k < i; k++, ik++, ki += n ) {
	  /* Store u/H in ith column of matrix. */
	  m[ki] = m[ik] / H;
	  /* Compute component of matrix * u, stored in y. */
	  y = 0.0;
	  for ( l = 0, kl = k*n, il = in; l <= k; l++, kl++, il++ )
	    y += m[kl]*m[il];
	  for ( l = k + 1, lk = k*( n + 1 ) + n, il = in + k + 1;
		l < i; l++, lk += n, il++ )
	    y += m[lk]*m[il];
	  /* Compute p, stored in unfilled part of o. */
	  o[k] = y/H;
	  x += o[k]*m[ik];
	}
	K = 0.5*x/H;

	/* Now reduce submatrix to tri-diagonal form.  Start by
           forming vector q and storing in o, overwriting p. */
	for ( k = 0, ik = in; k < i; k++, ik++ ) {
	  x = m[ik];
	  o[k] = y = o[k] - K*x;
	  /* Now apply submatrix reduction. */
	  for ( l = 0, kl = n*k, il = in; l <= k; l++, kl++, il++ )
	    m[kl] -= x*o[l] + y*m[il];
	}
      }
    }

    /* Last (top-right) 1x1 submatrix is trivial. */
    else
      o[i] = m[ij];

    /* Indicate whether this iteration was trivial or not. */
    d[i] = H;
  }

  /* Now compute transformation matrix (eigenvectors). */
  d[0] = 0.0;
  o[0] = 0.0;
  /* Loop over the n sub-tranformations. */
  for ( i = 0; i < n; i++ ) {
    in = i*n;
    ii = in + i;
    /* Skip first block and any other trivial transformations. */
    if ( d[i] ) {
      for ( k = 0; k < i; k++ ) {
	y = 0.0;
	/* Use u and u/H, stored in ith row and column of matrix, to
           compute P*Q. */
	for ( l = 0, il = in, lk = k; l < i; l++, il++, lk += n )
	  y += m[il]*m[lk];
	for ( l = 0, li = i, lk = k; l < i; l++, li += n, lk += n )
	  m[lk] -= y*m[li];
      }
    }

    /* Store diagonal element. */
    d[i] = m[ii];

    /* Reset row and column to identity for next iteration. */
    m[ii] = 1.0;
    for ( k = 0, ik = in, ki = i; k < i; k++, ik++, ki += n )
      m[ki] = m[ik] = 0.0;
  }

  /* Finished. */
  RETURN( stat );
}


static void
LALDTriDiagonalToDiagonal( LALStatus   *stat,
			   REAL8Vector *diag,
			   REAL8Array  *matrix,
			   REAL8Vector *offDiag )
{
  INT4 n, i, j, k, l;              /* matrix dimension and indecies */
  INT4 lk;                         /* matrix array index */
  UINT4 iteration;                 /* number of iterated transforms */
  REAL8 *m, *d, *o;                /* pointers to data arrays */
  REAL8 x, y, xx, yy, s, c, r, d2; /* temporary variables */

  INITSTATUS(stat);

  /* Check input fields. */
  ASSERT( diag, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( diag->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( offDiag, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( offDiag->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength->data, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( matrix->dimLength->length == 2, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  n = (INT4)( diag->length );
  ASSERT( n, stat, MATRIXUTILSH_EDIM, MATRIXUTILSH_MSGEDIM );
  ASSERT( n == (INT4)( offDiag->length ), stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  ASSERT( n == (INT4)( matrix->dimLength->data[0] ), stat,
	  MATRIXUTILSH_EDIM, MATRIXUTILSH_MSGEDIM );
  ASSERT( n == (INT4)( matrix->dimLength->data[1] ), stat,
	  MATRIXUTILSH_EDIM, MATRIXUTILSH_MSGEDIM );

  /* Set up pointers. */
  m = matrix->data;
  d = diag->data;
  o = offDiag->data;

  /* As in Numerical Recipes, we start by shifting the indexing of the
     off-diagonal vector. */
  for ( i = 1; i < n; i++ )
    o[i-1] = o[i];
  o[n-1] = 0.0;

  /* Loop over the i eigenvalues. */
  for ( i = 0; i < n; i++ ) {

    /* For each eigenvalue, loop until we've diagonalized that part of
       the matrix. */
    iteration = 0;
    do {

      /* Look for a small (effectively zero) off-diagonal element. */
      for ( j = i; j < n - 1; j++ ) {
	d2 = fabs( d[j] ) + fabs( d[j+1] );
	if ( fabs( o[j]/d2 ) < LAL_REAL8_EPS )
	  break;
      }

      if ( j != i ) {
	if ( iteration++ == EIGENINTERNALC_MAXITER ) {
	  ABORT( stat, MATRIXUTILSH_EITER, MATRIXUTILSH_MSGEITER );
	}
	x = 0.5*( d[i+1] - d[i] )/o[i];
	r = RSS( x, 1.0 );
	if ( x > 0 )
	  x = d[j] - d[i] + o[i]/( x + r );
	else
	  x = d[j] - d[i] + o[i]/( x - r );
	s = c = 1.0;
	yy = 0.0;

	/* Perform a plane rotation to diagonalize the submatrix, and
           Givens rotations to restore it to tridiagonal form. */
	for ( k = j - 1; k >= i; k-- ) {
	  xx = c*o[k];
	  y = s*o[k];
	  o[k+1] = ( r = RSS( x, y ) );

	  /* Check for underflow. */
	  if ( r <= LAL_REAL8_MIN ) {
	    d[k+1] -= yy;
	    o[j] = 0.0;
	    break;
	  }

	  /* Compute diagonal element. */
	  s = y/r;
	  c = x/r;
	  x = d[k+1] - yy;
	  r = ( d[k] - x )*s + 2.0*c*xx;
	  d[k+1] = x + ( yy = s*r );
	  x = c*r - xx;

	  /* Accumulate transformation matrix. */
	  for ( l = 0, lk = k; l < n; l++, lk += n ) {
	    y = m[lk+1];
	    m[lk+1] = s*m[lk] + c*y;
	    m[lk] = c*m[lk] - s*y;
	  }
	}

	/* Check for underflow. */
	if ( r <= LAL_REAL8_MIN && k >= i )
	  continue;

	/* Adjust diagonal and off-diagonal elements, and continue
           loop. */
	d[i] -= yy;
	o[i] = x;
	o[j] = 0.0;
      }
    } while ( j != i );
  }

  /* Done. */
  RETURN( stat );
}


/**
 * \defgroup Eigen_c Module Eigen.c
 * \ingroup MatrixUtils_h
 * \author Creighton, T. D.
 *
 * \brief Routines to compute eigenvalues and eigenvectors.
 *
 * ### Description ###
 *
 * <tt>LALSSymmetricEigenVectors()</tt> and
 * <tt>LALDSymmetricEigenVectors()</tt> compute the eigenvalues and
 * eigenvectors of the square matrix <tt>*matrix</tt>.  The
 * eigen\e values are stored in <tt>*values</tt>, which must be
 * pre-allocated to the same length as each of the dimensions of
 * <tt>*matrix</tt>.  The eigen\e vectors are computed in-place and
 * returned in <tt>*matrix</tt>: on return, each column of <tt>*matrix</tt>
 * stores the eigenvector of the corresponding eigenvalue as a column
 * vector.  If you don't want the input matrix to be changed, make a copy
 * of it first.
 *
 * ### Algorithm ###
 *
 * A square matrix \f$\mathsf{M}^a{}_b\f$ is said to have an eigenvalue
 * \f$\lambda\f$ and associated eigenvector \f$\mathsf{x}^a\f$ if the following
 * equation holds:
 * \f{equation}{
 * \label{eq_eigenvalue}
 * \mathsf{M}^a{}_b \mathsf{x}^b = \lambda\mathsf{x}^a
 * \f}
 * Generically an \f$N\times N\f$ matrix has \f$N\f$ distinct eigenvalues each
 * with a unique (up to a multiplicative constant) associated
 * eigenvector.  In certain cases, though, the system is
 * \e degenerate: that is, some of the \f$N\f$ eigenvalues may be the
 * \e same, and the corresponding eigenvectors may not be unique.
 *
 * We are concerned with the particular case of real, symmetric matrices,
 * which occur in a broad class of problems.  These matrices have certain
 * useful properties.  First, for any Hermitian matrix (including real
 * symmetric matrices) the eigenvalues are all real.  Second, for any
 * matrix that commutes with its Hermitian conjugate (including Hermitian
 * and unitary matrices, or, for real matrices, symmetric and
 * orthogonal), the eigenvectors \f${\mathsf{x}_{(i)}}^a\f$ of distinct
 * eigenvalues \f$\lambda_{(i)}\f$ are \e orthogonal: that is,
 * \f$(\mathsf{x}_{(i)}^T)_a \mathsf{x}_{(j)}{}^a=0\f$ if
 * \f$\lambda_{(i)}\neq\lambda_{(j)}\f$.  Furthermore, we note that if two or
 * more linearly independent eigenvectors have the \e same
 * eigenvalue, then any linear combination of them is also an eigenvector
 * with that same eigenvalue: we are therefore free to choose a basis of
 * eigenvectors that is orthonormal.  Thirdly, for any matrix that
 * commutes with its conjugate, the complete set of eigenvectors spans
 * the entire \f$N\f$-dimensional vector space, so the orthonormal basis
 * above is a \e complete orthonormal basis.
 *
 * If we construct a square matrix \f$\mathsf{X}^a{}_b\f$ whose columns are
 * the orthonormal eigenvectors \f$X^i{}_j={x_{(j)}}^i\f$, then it is clear
 * (from the orthonormality of \f${\mathsf{x}_{(j)}}^a\f$) that
 * \f$\mathsf{X}^a{}_b\f$ is orthogonal.  Furtthermore, the eigenvalue
 * \eqref{eq_eigenvalue} gives:
 * \f{equation}{
 * {(\mathsf{X}^{-1}){}^a{}_b}\mathsf{M}^b{}_c\mathsf{X}^c{}_d =
 * \left[\begin{array}{cccc}
 * \lambda_{(1)} &       0       & \cdots &       0       \\
 * 0       & \lambda_{(2)} & \cdots &       0       \\
 * \vdots     &    \vdots     & \ddots &    \vdots     \\
 * 0       &       0       & \cdots & \lambda_{(N)}
 * \end{array}\right] \;,
 * \f}
 * where \f$\lambda_{(i)}\f$ are the eigenvalues of the corresponding
 * eigenvectors (with the possibility that some of these eigenvalues have
 * the same value).  That is, the matrix \f$\mathsf{M}^b{}_c\f$ can be
 * \e diagonalized by an orthogonal similarity transformation, and
 * the transformation matrix gives us the eigenvectors.  The process of
 * solving the eigenvalue equation is thus equivalent to the problem of
 * diagonalizing the matrix.
 *
 * For a general \f$N\times N\f$ symmetric matrix, there is no finite
 * algorithm that exactly diagonalizes the matrix.  Most numerical
 * eigenvalue solvers use a process of iterative orthogonal
 * transformations, where each transformation is chosen to reduce the sum
 * of the squares of the off-diagonal elements: the matrix uniformly
 * converges to a diagonal form under successive transformations.  This
 * convergence can be speeded if one first transforms the matrix into a
 * \e tridiagonal form (where the only nonzero elements are on the
 * diagonal or immediately adjacent to it); this transformation
 * \e can be accomplished deterministically in a finite number of
 * steps.  This is the approach advocated in Chapter 11
 * of \cite ptvf1992 , and taken in this module.
 *
 * The routines in call the routines 
 * first to reduce the matrix to tridiagonal
 * form, and then to reduce it to diagonal form iteratively.  The
 * routines <tt>LALSSymmetricEigenVectors()</tt> and
 * <tt>LALDSymmetricEigenVectors()</tt> call the routines
 * <tt>LALSSymmetricToTriDiagonal()</tt> and
 * <tt>LALDSymmetricToTriDiagonal()</tt> to tri-diagonalize the matrix,
 * then <tt>LALSTriDiagonalToDiagonal()</tt> and
 * <tt>LALDTriDiagonalToDiagonal()</tt> to diagonalize it.
 *
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

/*@}*/
