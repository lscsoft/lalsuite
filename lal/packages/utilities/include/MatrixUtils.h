/********************************** <lalVerbatim file="MatrixUtilsHV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\providecommand{\bd}{\mbox{\boldmath$\delta$\unboldmath}}

\section{Header \texttt{MatrixUtils.h}}
\label{s:MatrixUtils.h}

Provides routines to solve linear systems.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/MatrixUtils.h>
\end{verbatim}

This header covers routines to solve linear systems of equations and
eigensystems, using algorithms adapted from Chapters~2 and~11 of
Numerical Recipes~\cite{ptvf:1992}.  The only routines at present are
for computing eigenvalues and eigenvectors of real symmetric matrices.
Routines for inverting or computing the determinant of arbitrary
square matrices will likely follow.

\subsubsection*{Notation}

A \emph{matrix} is represented in LAL by a \verb@<datatype>Array@
structure with a \verb@dimLength->length@ field of 2; the
\verb@dimLength->data@ field gives the dimensions $[M,N]$ of the
matrix.  Using the place-index notation common in tensor calculus, a
matrix is a two-index tensor:
\begin{equation}
\mathsf{A}^a{}_b = \left[\begin{array}{cccc}
	A^1{}_1 & A^1{}_2 & \cdots & A^1{}_N \\
	A^2{}_1 & A^2{}_2 & \cdots & A^2{}_N \\
	\vdots  & \vdots  & \ddots & \vdots  \\
	A^M{}_1 & A^M{}_2 & \cdots & A^M{}_N
\end{array}\right] \;,
\end{equation}
that is, the first (raised) index represents the row number and the
second (lowered) index the column number.  The standard C array
structure declares this object as, say, \verb@float a[M][N]@, where
the element $A^i{}_j$ is stored in \verb@a[i][j]@.  The LAL array
structure \verb@REAL4Array a@ stores data in a flattened block of
memory, where the element $A^i{}_j$ is stored in \verb@a.data[i*M+j]@.

A \emph{row vector} is a matrix with only one row ($M=1$).  In the
above place-index notation, it is represented with a single lowered
index:
\begin{equation}
\mathsf{r}_a = \left[\begin{array}{cccc} r_1 & r_2 & \cdots & r_N
	\end{array}\right] \;.
\end{equation}
A \emph{column vector} is a matrix with only one column ($N=1$).  In
the above place-index notation, it is represented with a single raised
index:
\begin{equation}
\mathsf{c}^a = \left[\begin{array}{c} c^1 \\ c^2 \\ \vdots \\ c^M
	\end{array}\right] \;.
\end{equation}
In LAL, both of these objects are conventionally represented as a LAL
vector structure.  Whether the object is to be used as a row or column
vector must be determined from context; it is not specified by the
datatype.

\subsubsection*{Properties}

The basic matrix operations are addition, scalar multiplication, and
vector multiplication.  We assume the reader is familiar with these.
In addition, we will refer to the following unary operations on
\emph{square} matrices:

\textit{Inversion:} The inverse $(\mathsf{A}^{-1}){}^a{}_b$ of a
matrix $\mathsf{A}^a{}_b$ is one such that their matrix product is the
identity matrix $\bd^a{}_b$ (whose elements $\delta^i{}_j$ are just
the Kr\"onecker delta function).

\textit{Transposition:} The transpose $(\mathsf{A}^T){}^a{}_b$ of a
matrix $\mathsf{A}^a{}_b$ is given by interchanging the indecies on
each component: $(A^T){}^i{}_j=A^j{}_i$.

\textit{Conjugation:} The Hermitian conjugate (adjoint)
$(\mathsf{A}^\dag){}^a{}_b$ of a matrix $\mathsf{A}^a{}_b$ is given by
interchanging the indecies and taking the complex conjugate of each
component: $(A^\dag){}^i{}_j={A^j{}_i}^*$.

A matrix that is identical to its own transpose is called
\emph{symmetric}.  A matrix whose transpose is identical to the
original matrix's inverse is called \emph{orthogonal}.  A matrix that
is identical to its own Hermitian conjugate is called \emph{Hermitian}
(or \emph{self-adjoint}.  A matrix whose Hermitian conjugate is
identical to the original matrix's inverse is called \emph{unitary}.

At present, the routines under this header only deal with \emph{real}
matrices (i.e.\ matrices, vectors, and scalars whose components are
all real).  In this case, symmetric is equivalent to Hermitian, and
orthogonal is equivalent to unitary.

%"

******************************************************* </lalLaTeX> */

#ifndef _MATRIXUTILS_H
#define _MATRIXUTILS_H

#include <lal/LALStdlib.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( MATRIXUTILSH, "$Id$" );

/********************************************************** <lalLaTeX>
\subsection*{Error conditions}
****************************************** </lalLaTeX><lalErrTable> */
#define MATRIXUTILSH_ENUL  1
#define MATRIXUTILSH_EDIM  2
#define MATRIXUTILSH_EITER 3
#define MATRIXUTILSH_ESING 4
#define MATRIXUTILSH_EMEM  5

#define MATRIXUTILSH_MSGENUL  "Unexpected null pointer in arguments"
#define MATRIXUTILSH_MSGEDIM  "Bad matrix dimensions"
#define MATRIXUTILSH_MSGEITER "Did not converge after maximum iterations"
#define MATRIXUTILSH_MSGESING "Singular matrix"
#define MATRIXUTILSH_MSGEMEM  "Memory allocation error"
/*************************************************** </lalErrTable> */

/* <lalLaTeX>
\vfill{\footnotesize\input{MatrixUtilsHV}}
</lalLaTeX> */


/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{MatrixOpsC}
</lalLaTeX> */

void
LALI2MatrixAdd( LALStatus *, INT2Array *out, INT2Array *in1, INT2Array *in2 );
void
LALI4MatrixAdd( LALStatus *, INT4Array *out, INT4Array *in1, INT4Array *in2 );
void
LALI8MatrixAdd( LALStatus *, INT8Array *out, INT8Array *in1, INT8Array *in2 );
void
LALU2MatrixAdd( LALStatus *, UINT2Array *out, UINT2Array *in1, UINT2Array *in2 );
void
LALU4MatrixAdd( LALStatus *, UINT4Array *out, UINT4Array *in1, UINT4Array *in2 );
void
LALU8MatrixAdd( LALStatus *, UINT8Array *out, UINT8Array *in1, UINT8Array *in2 );
void
LALSMatrixAdd( LALStatus *, REAL4Array *out, REAL4Array *in1, REAL4Array *in2 );
void
LALDMatrixAdd( LALStatus *, REAL8Array *out, REAL8Array *in1, REAL8Array *in2 );
void
LALCMatrixAdd( LALStatus *, COMPLEX8Array *out, COMPLEX8Array *in1, COMPLEX8Array *in2 );
void
LALZMatrixAdd( LALStatus *, COMPLEX16Array *out, COMPLEX16Array *in1, COMPLEX16Array *in2 );

void
LALI2MatrixMultiply( LALStatus *, INT2Array *out, INT2Array *in1, INT2Array *in2 );
void
LALI4MatrixMultiply( LALStatus *, INT4Array *out, INT4Array *in1, INT4Array *in2 );
void
LALI8MatrixMultiply( LALStatus *, INT8Array *out, INT8Array *in1, INT8Array *in2 );
void
LALU2MatrixMultiply( LALStatus *, UINT2Array *out, UINT2Array *in1, UINT2Array *in2 );
void
LALU4MatrixMultiply( LALStatus *, UINT4Array *out, UINT4Array *in1, UINT4Array *in2 );
void
LALU8MatrixMultiply( LALStatus *, UINT8Array *out, UINT8Array *in1, UINT8Array *in2 );
void
LALSMatrixMultiply( LALStatus *, REAL4Array *out, REAL4Array *in1, REAL4Array *in2 );
void
LALDMatrixMultiply( LALStatus *, REAL8Array *out, REAL8Array *in1, REAL8Array *in2 );
void
LALCMatrixMultiply( LALStatus *, COMPLEX8Array *out, COMPLEX8Array *in1, COMPLEX8Array *in2 );
void
LALZMatrixMultiply( LALStatus *, COMPLEX16Array *out, COMPLEX16Array *in1, COMPLEX16Array *in2 );

void
LALI2MatrixTranspose( LALStatus *, INT2Array *out, INT2Array *in1 );
void
LALI4MatrixTranspose( LALStatus *, INT4Array *out, INT4Array *in1 );
void
LALI8MatrixTranspose( LALStatus *, INT8Array *out, INT8Array *in1 );
void
LALU2MatrixTranspose( LALStatus *, UINT2Array *out, UINT2Array *in1 );
void
LALU4MatrixTranspose( LALStatus *, UINT4Array *out, UINT4Array *in1 );
void
LALU8MatrixTranspose( LALStatus *, UINT8Array *out, UINT8Array *in1 );
void
LALSMatrixTranspose( LALStatus *, REAL4Array *out, REAL4Array *in1 );
void
LALDMatrixTranspose( LALStatus *, REAL8Array *out, REAL8Array *in1 );
void
LALCMatrixTranspose( LALStatus *, COMPLEX8Array *out, COMPLEX8Array *in1 );
void
LALZMatrixTranspose( LALStatus *, COMPLEX16Array *out, COMPLEX16Array *in1 );

void
LALCMatrixAdjoint( LALStatus *, COMPLEX8Array *out, COMPLEX8Array *in1 );
void
LALZMatrixAdjoint( LALStatus *, COMPLEX16Array *out, COMPLEX16Array *in1 );

/* <lalLaTeX>
\newpage\input{DetInverseC}
</lalLaTeX> */
void
LALSMatrixDeterminant( LALStatus *, REAL4 *det, REAL4Array *matrix );

void
LALSMatrixInverse( LALStatus *, REAL4 *det, REAL4Array *matrix, REAL4Array *inverse );

void
LALSMatrixDeterminantErr( LALStatus *, REAL4 det[2], REAL4Array *matrix, REAL4Array *matrixErr );

void
LALDMatrixDeterminant( LALStatus *, REAL8 *det, REAL8Array *matrix );

void
LALDMatrixInverse( LALStatus *, REAL8 *det, REAL8Array *matrix, REAL8Array *inverse );

void
LALDMatrixDeterminantErr( LALStatus *, REAL8 det[2], REAL8Array *matrix, REAL8Array *matrixErr );

/* <lalLaTeX>
\newpage\input{DetInverseInternalC}
</lalLaTeX> */
void
LALSLUDecomp( LALStatus   *,
	      INT2        *sgn,
	      REAL4Array  *matrix,
	      UINT4Vector *indx );

void
LALSLUBackSub( LALStatus   *,
	       REAL4Vector *vector,
	       REAL4Array  *matrix,
	       UINT4Vector *indx );

void
LALDLUDecomp( LALStatus   *,
	      INT2        *sgn,
	      REAL8Array  *matrix,
	      UINT4Vector *indx );

void
LALDLUBackSub( LALStatus   *,
	       REAL8Vector *vector,
	       REAL8Array  *matrix,
	       UINT4Vector *indx );

/* <lalLaTeX>
\newpage\input{DetInverseTestC}
</lalLaTeX> */

/* <lalLaTeX>
\newpage\input{EigenC}
</lalLaTeX> */
void
LALSSymmetricEigenVectors( LALStatus *, REAL4Vector *values, REAL4Array *matrix );

void
LALSSymmetricEigenValues( LALStatus *, REAL4Vector *values, REAL4Array *matrix );

void
LALDSymmetricEigenVectors( LALStatus *, REAL8Vector *values, REAL8Array *matrix );

void
LALDSymmetricEigenValues( LALStatus *, REAL8Vector *values, REAL8Array *matrix );

/* <lalLaTeX>
\newpage\input{EigenInternalC}
</lalLaTeX> */
void
LALSSymmetricToTriDiagonal( LALStatus   *,
			    REAL4Vector *diag,
			    REAL4Array  *matrix,
			    REAL4Vector *offDiag );

void
LALSSymmetricToTriDiagonal2( LALStatus   *,
			      REAL4Vector *diag,
			      REAL4Array  *matrix,
			      REAL4Vector *offDiag );

void
LALSTriDiagonalToDiagonal( LALStatus   *,
			   REAL4Vector *diag,
			   REAL4Array  *matrix,
			   REAL4Vector *offDiag );

void
LALSTriDiagonalToDiagonal2( LALStatus   *,
			    REAL4Vector *diag,
			    REAL4Array  *matrix,
			    REAL4Vector *offDiag );

void
LALDSymmetricToTriDiagonal( LALStatus   *,
			    REAL8Vector *diag,
			    REAL8Array  *matrix,
			    REAL8Vector *offDiag );

void
LALDSymmetricToTriDiagonal2( LALStatus   *,
			     REAL8Vector *diag,
			     REAL8Array  *matrix,
			     REAL8Vector *offDiag );

void
LALDTriDiagonalToDiagonal( LALStatus   *,
			   REAL8Vector *diag,
			   REAL8Array  *matrix,
			   REAL8Vector *offDiag );

void
LALDTriDiagonalToDiagonal2( LALStatus   *,
			    REAL8Vector *diag,
			    REAL8Array  *matrix,
			    REAL8Vector *offDiag );

/* <lalLaTeX>
\newpage\input{EigenTestC}
</lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif /* _MATRIXUTILS_H */
