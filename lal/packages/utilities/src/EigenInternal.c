/******************************** <lalVerbatim file="EigenInternalCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{EigenInternal.c}}
\label{ss:EigenInternal.c}

Internal routines used to compute eigenvalues and eigenvectors.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{EigenInternalCP}
\idx{LALSSymmetricToTriDiagonal()}
\idx{LALSTriDiagonalToDiagonal()}
\idx{LALDSymmetricToTriDiagonal()}
\idx{LALDTriDiagonalToDiagonal()}
\idx{LALSSymmetricToTriDiagonal2()}
\idx{LALSTriDiagonalToDiagonal2()}
\idx{LALDSymmetricToTriDiagonal2()}
\idx{LALDTriDiagonalToDiagonal2()}

\subsubsection*{Description}

These functions are called by the routines in \verb@Eigen.c@ to
compute eigenvalues and eigenvectors of a symmetric square matrix
\verb@*matrix@.  They are provided because, in some circumstances,
users may find it useful to work with the partially-reduced results.

\verb@LALSSymmetricToTriDiagonal()@ and
\verb@LALDSymmetricToTriDiagonal()@ reduce the symmetric square matrix
\verb@*matrix@ to tridiagonal form.  The vectors \verb@*diag@ and
\verb@*offDiag@ must be allocated to the same length as each matrix
dimension; on return, \verb@*diag@ stores the diagonal elements and
\verb@*offDiag@ the off-diagonal elements (\verb@offDiag->data[0]@ is
meaningless and will be set to zero).  The tri-diagonalization is done
in-place, so that on return \verb@*matrix@ will store the
transformation matrix that brings the original input matrix into
tridiagonal form.  If you don't want the input matrix to be changed,
make a copy of it first.

\verb@LALSTriDiagonalToDiagonal()@ and
\verb@LALDTriDiagonalToDiagonal()@ take a symmetric tridiagonal matrix
and compute its eigenvalues and eigenvectors.  On input, \verb@*diag@
stores the diagonal elements of the matrix, and \verb@*offDiag@ the
off-diagonal elements, with \verb@offDiag->data[0]@ arbitrary; on
return, \verb@*diag@ will store the eigenvalues and \verb@*offDiag@
will store zeroes.  The matrix \verb@*matrix@ should store the
orthogonal transformation matrix that brought the original symmetric
matrix into tri-diagonal form (as returned by the above routines), or
the identity matrix if the tri-diagonal matrix \emph{is} the original
matrix of interest; on return, \verb@*matrix@ will store the
orthogonal transformation matrix that diagonalizes the original
matrix: its columns are the eigenvectors of the original matrix.

\subsubsection*{Algorithm}

The tri-diagonalizing routines follow the Householder reduction method
described in Sec.~11.2 of~\cite{ptvf:1992}; they are essentially
re-implementations of the Numerical Recipes routine \verb@tred2()@.
For large $N$, their operation count is approximately $4N^3/3$, or
$2N^3/3$ for the routines that ignore eigenvectors.  These routines
explicitly enforce symmetry, by only using the lower-left triangle of
the matrix as input.

The diagonalizing routines follow the QL algorithm with implicit
shifts described in Sec.~11.3 of~\cite{ptvf:1992}; they are
essentially re-implementations of the Numerical Recipes routine
\verb@tqli()@.  Depending on the number of iterations required, their
operation count is roughly $\sim30N^2$, plus $\sim3N^3$ if
eigenvectors are also being computed.

The diagonalizing routines can fail if they fail to converge rapidly
enough.  The discussion in~\cite{ptvf:1992} does not go into much
detail about when this is likely to occur, except to note that
degenerate eigenvalues converge more slowly.  If the routines fail in
this way, \verb@diag@, \verb@matrix@, and \verb@offDiag@ will all be
left in indeterminate states.

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{EigenInternalCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/MatrixUtils.h>

#define EIGENINTERNALC_MAXITER 30 /* max. number of iterations in
				     diagonalizing */

NRCSID( EIGENINTERNALC, "$Id$" );

/* A quick macro to take the square root of the sum of squares,
   without overflowing or underflowing. */
#define RSS( a, b ) ( fabs( a ) > fabs( b ) ?                        \
fabs( a )*sqrt( 1.0 + ((b)/(a))*((b)/(a)) ) :                        \
( (b) == 0.0 ? 0.0 : fabs( b )*sqrt( 1.0 + ((a)/(b))*((a)/(b)) ) ) )

/* <lalVerbatim file="EigenInternalCP"> */
void
LALSSymmetricToTriDiagonal( LALStatus   *stat,
			    REAL4Vector *diag,
			    REAL4Array  *matrix,
			    REAL4Vector *offDiag )
{ /* </lalVerbatim> */
  UINT4 n, in;                          /* matrix dimension (times i) */
  UINT4 i, k, l;                        /* dimension indecies */
  UINT4 ij, ik, ki, il, li, kl, lk, ii; /* matrix array indecies */
  REAL4 *m, *d, *o;                     /* pointers to data arrays */
  REAL4 H, K; /* scalars defined in Eqs. (11.2.4) and (11.2.11) of NR */
  REAL4 x, y, scale; /* temporary variables */

  INITSTATUS( stat, "LALSSymmetricToTriDiagonal", EIGENINTERNALC );

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


/* <lalVerbatim file="EigenInternalCP"> */
void
LALSSymmetricToTriDiagonal2( LALStatus   *stat,
			     REAL4Vector *diag,
			     REAL4Array  *matrix,
			     REAL4Vector *offDiag )
{ /* </lalVerbatim> */
  UINT4 n, in;              /* matrix dimension (times i) */
  UINT4 i, k, l;            /* dimension indecies */
  UINT4 ij, ik, il, kl, lk; /* matrix array indecies */
  REAL4 *m, *d, *o;         /* pointers to data arrays */
  REAL4 H, K; /* scalars defined in Eqs. (11.2.4) and (11.2.11) of NR */
  REAL4 x, y, scale; /* temporary variables */

  INITSTATUS( stat, "LALSSymmetricToTriDiagonal2", EIGENINTERNALC );

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
	for ( k = 0, ik = in; k < i; k++, ik++ ) {
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
  }

  /* Now set redundant o element to zero, and fill d. */
  o[0] = 0.0;
  for ( i = 0; i < n; i++ )
    d[i] = m[i*( n + 1 )];

  /* Finished. */
  RETURN( stat );
}


/* <lalVerbatim file="EigenInternalCP"> */
void
LALSTriDiagonalToDiagonal( LALStatus   *stat,
			   REAL4Vector *diag,
			   REAL4Array  *matrix,
			   REAL4Vector *offDiag )
{ /* </lalVerbatim> */
  INT4 n, i, j, k, l;              /* matrix dimension and indecies */
  INT4 lk;                         /* matrix array index */
  UINT4 iteration;                 /* number of iterated transforms */
  REAL4 *m, *d, *o;                /* pointers to data arrays */
  REAL4 x, y, xx, yy, s, c, r, d2; /* temporary variables */

  INITSTATUS( stat, "LALSTriDiagonalToDiagonal", EIGENINTERNALC );

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


/* <lalVerbatim file="EigenInternalCP"> */
void
LALSTriDiagonalToDiagonal2( LALStatus   *stat,
			    REAL4Vector *diag,
			    REAL4Array  *matrix,
			    REAL4Vector *offDiag )
{ /* </lalVerbatim> */
  INT4 n, i, j, k;                 /* matrix dimension and indecies */
  UINT4 iteration;                 /* number of iterated transforms */
  REAL4 *d, *o;                    /* pointers to data arrays */
  REAL4 x, y, xx, yy, s, c, r, d2; /* temporary variables */

  INITSTATUS( stat, "LALSTriDiagonalToDiagonal2", EIGENINTERNALC );

  /* Check input fields. */
  ASSERT( diag, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( diag->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( offDiag, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( offDiag->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  n = (INT4)( diag->length );
  ASSERT( n, stat, MATRIXUTILSH_EDIM, MATRIXUTILSH_MSGEDIM );
  ASSERT( n == (INT4)( offDiag->length ), stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );

  /* Dummy test for dummy parameter. */
  if ( !matrix )
    matrix = NULL;

  /* Set up pointers. */
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


/* <lalVerbatim file="EigenInternalCP"> */
void
LALDSymmetricToTriDiagonal( LALStatus   *stat,
			    REAL8Vector *diag,
			    REAL8Array  *matrix,
			    REAL8Vector *offDiag )
{ /* </lalVerbatim> */
  UINT4 n, in;                          /* matrix dimension (times i) */
  UINT4 i, k, l;                        /* dimension indecies */
  UINT4 ij, ik, ki, il, li, kl, lk, ii; /* matrix array indecies */
  REAL8 *m, *d, *o;                     /* pointers to data arrays */
  REAL8 H, K; /* scalars defined in Eqs. (11.2.4) and (11.2.11) of NR */
  REAL8 x, y, scale; /* temporary variables */

  INITSTATUS( stat, "LALSSymmetricToTriDiagonal", EIGENINTERNALC );

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


/* <lalVerbatim file="EigenInternalCP"> */
void
LALDSymmetricToTriDiagonal2( LALStatus   *stat,
			     REAL8Vector *diag,
			     REAL8Array  *matrix,
			     REAL8Vector *offDiag )
{ /* </lalVerbatim> */
  UINT4 n, in;              /* matrix dimension (times i) */
  UINT4 i, k, l;            /* dimension indecies */
  UINT4 ij, ik, il, kl, lk; /* matrix array indecies */
  REAL8 *m, *d, *o;         /* pointers to data arrays */
  REAL8 H, K; /* scalars defined in Eqs. (11.2.4) and (11.2.11) of NR */
  REAL8 x, y, scale; /* temporary variables */

  INITSTATUS( stat, "LALSSymmetricToTriDiagonal2", EIGENINTERNALC );

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
	for ( k = 0, ik = in; k < i; k++, ik++ ) {
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
  }

  /* Now set redundant o element to zero, and fill d. */
  o[0] = 0.0;
  for ( i = 0; i < n; i++ )
    d[i] = m[i*( n + 1 )];

  /* Finished. */
  RETURN( stat );
}


/* <lalVerbatim file="EigenInternalCP"> */
void
LALDTriDiagonalToDiagonal( LALStatus   *stat,
			   REAL8Vector *diag,
			   REAL8Array  *matrix,
			   REAL8Vector *offDiag )
{ /* </lalVerbatim> */
  INT4 n, i, j, k, l;              /* matrix dimension and indecies */
  INT4 lk;                         /* matrix array index */
  UINT4 iteration;                 /* number of iterated transforms */
  REAL8 *m, *d, *o;                /* pointers to data arrays */
  REAL8 x, y, xx, yy, s, c, r, d2; /* temporary variables */

  INITSTATUS( stat, "LALSTriDiagonalToDiagonal", EIGENINTERNALC );

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


/* <lalVerbatim file="EigenInternalCP"> */
void
LALDTriDiagonalToDiagonal2( LALStatus   *stat,
			    REAL8Vector *diag,
			    REAL8Array  *matrix,
			    REAL8Vector *offDiag )
{ /* </lalVerbatim> */
  INT4 n, i, j, k;                 /* matrix dimension and indecies */
  UINT4 iteration;                 /* number of iterated transforms */
  REAL8 *d, *o;                    /* pointers to data arrays */
  REAL8 x, y, xx, yy, s, c, r, d2; /* temporary variables */

  INITSTATUS( stat, "LALSTriDiagonalToDiagonal2", EIGENINTERNALC );

  /* Check input fields. */
  ASSERT( diag, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( diag->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( offDiag, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( offDiag->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  n = (INT4)( diag->length );
  ASSERT( n, stat, MATRIXUTILSH_EDIM, MATRIXUTILSH_MSGEDIM );
  ASSERT( n == (INT4)( offDiag->length ), stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );

  /* Dummy test for dummy parameter. */
  if ( !matrix )
    matrix = NULL;

  /* Set up pointers. */
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
