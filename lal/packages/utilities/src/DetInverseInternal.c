/*************************** <lalVerbatim file="DetInverseInternalCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{DetInverseInternal.c}}
\label{ss:DetInverseInternal.c}

Internal routines used to compute matrix determinants and inverses.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{DetInverseInternalCP}
\idx{LALSLUDecomp()}
\idx{LALSLUBackSub()}
\idx{LALDLUDecomp()}
\idx{LALDLUBackSub()}

\subsubsection*{Description}

These functions are called by the routines in \verb@DetInverse.c@ to
compute the determinant and inverse of a nondegenerate square matrix
\verb@*matrix@.  They are useful routines in their own right, though,
so they are made publically available.

\verb@LALSLUDecomp()@ and \verb@LALDLUDecomp()@ replace \verb@*matrix@
with an LU decomposition of a \emph{row-wise permutation} of itself.
The output parameter \verb@*indx@ stores the permutation, and the
output \verb@*sgn@ records the sign of the permutation.

\verb@LALSLUBackSub()@ and \verb@LALDLUBackSub()@ take the permuted
LU-decomposed matrix returned by the above routine, and
back-substitute the vector \verb@*vector@ representing $\mathsf{v}^a$
in Eq.~(\ref{eq:linear-system}), to compute the vector $\mathsf{x}^b$.
This is returned in-place in \verb@*vector@.  The input parameter
\verb@*indx@ is the list of row permutations returned by the above
routines.

\subsubsection*{Algorithm}

LU decomposition is performed by Crout's algorithm, described in
Sec.~2.3 of~\cite{ptvf:1992}; the routines in this module are
essentially re-implementations of the Numerical Recipes routines
\verb@ludcmp()@ and \verb@lubksub()@.  For large $N$, their operation
counts are approximately $N^3/3$ and $N^2$, respectively.

One difference between \verb@ludcmp()@ in~\cite{ptvf:1992} and the
routines \verb@LALSLUDecomp()@ and \verb@LALDLUDecomp()@ in this
module is the way in which singular matrices are handled.
In~\cite{ptvf:1992}, there is a distinction between between a
manifestly singular matrix (where an entire row of the matrix is zero)
and a numerically singular matrix (if a diagonal element in the
decomposed matrix turns out to be zero).  In the former case, they
raise an error signal; in the latter, they replace the offending
element with a ``tiny'' but nonzero number and continue.  This
treatment does not strike the present author as satisfactory.

Instead, the routines \verb@LALSLUDecomp()@ and \verb@LALDLUDecomp()@
will \emph{always} return successfully, even with a singular matrix,
but will \emph{not} ``adjust away'' a numerical singularity.  Instead,
they will signal the presence of the singularity in two ways: First,
they will set the permutation sign \verb@*sgn@ to zero; second, they
will set \emph{all} elements of the \verb@*indx@ vector yo zero.  This
ensures that routines computing the determinant (whose sign depends on
\verb@*sgn@) will correctly give a zero determinant, while the
meaningless \verb@*indx@ provides a simple sanity check for routines
such as \verb@LALSLUBackSub()@ and \verb@LALDLUBackSub()@ that attempt
to invert the linear system.  Note that the returned value of
\verb@*matrix@ will be meaningless garbage.

\subsubsection*{Uses}
\begin{verbatim}
LALMalloc()                     LALFree()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{DetInverseInternalCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/MatrixUtils.h>

NRCSID( DETINVERSEINTERNALC, "$Id$" );

/* <lalVerbatim file="DetInverseInternalCP"> */
void
LALSLUDecomp( LALStatus   *stat,
	      INT2        *sgn,
	      REAL4Array  *matrix,
	      UINT4Vector *indx )
{ /* </lalVerbatim> */
  UINT4 n, imax = 0;    /* matrix dimension and pivot index */
  UINT4 i, j, k;        /* dimension indecies */
  UINT4 ij, ik, kj, jk; /* matrix array indecies */
  UINT4 *d;             /* pointer to indx->data */
  REAL4 *s, *m;         /* pointers to data arrays */
  REAL4 tmp, max, sum;  /* temporary computation variables */

  INITSTATUS( stat, "LALSLUDecomp", DETINVERSEINTERNALC );

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


/* <lalVerbatim file="DetInverseInternalCP"> */
void
LALSLUBackSub( LALStatus   *stat,
	       REAL4Vector *vector,
	       REAL4Array  *matrix,
	       UINT4Vector *indx )
{ /* </lalVerbatim> */
  INT4 n;             /* matrix dimension */
  INT4 i, j;          /* dimension indecies */
  INT4 ii, ij;        /* matrix array indecies */
  INT4 id, jmin = -1; /* permuted and first nonzero vector indecies */
  UINT4 *d;           /* pointer to indx->data */
  REAL4 *v, *m;       /* pointers to data arrays */
  REAL4 sum;          /* temporary computation variable */

  INITSTATUS( stat, "LALSLUBackSub", DETINVERSEINTERNALC );

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


/* <lalVerbatim file="DetInverseInternalCP"> */
void
LALDLUDecomp( LALStatus   *stat,
	      INT2        *sgn,
	      REAL8Array  *matrix,
	      UINT4Vector *indx )
{ /* </lalVerbatim> */
  UINT4 n, imax = 0;    /* matrix dimension and pivot index */
  UINT4 i, j, k;        /* dimension indecies */
  UINT4 ij, ik, kj, jk; /* matrix array indecies */
  UINT4 *d;             /* pointer to indx->data */
  REAL8 *s, *m;         /* pointers to data arrays */
  REAL8 tmp, max, sum;  /* temporary computation variables */

  INITSTATUS( stat, "LALDLUDecomp", DETINVERSEINTERNALC );

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


/* <lalVerbatim file="DetInverseInternalCP"> */
void
LALDLUBackSub( LALStatus   *stat,
	       REAL8Vector *vector,
	       REAL8Array  *matrix,
	       UINT4Vector *indx )
{ /* </lalVerbatim> */
  INT4 n;             /* matrix dimension */
  INT4 i, j;          /* dimension indecies */
  INT4 ii, ij;        /* matrix array indecies */
  INT4 id, jmin = -1; /* permuted and first nonzero vector indecies */
  UINT4 *d;           /* pointer to indx->data */
  REAL8 *v, *m;       /* pointers to data arrays */
  REAL8 sum;          /* temporary computation variable */

  INITSTATUS( stat, "LALDLUBackSub", DETINVERSEINTERNALC );

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
