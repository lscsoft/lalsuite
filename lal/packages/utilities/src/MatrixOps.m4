changecom(`/*',`*/')dnl
/************************************ <lalVerbatim file="MatrixOpsCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{MatrixOps.c}}
\label{ss:MatrixOps.c}

Routines to perform basic matrix operations.

\subsubsection*{Prototypes}
\vspace{0.1in}
\begin{verbatim}
void
LAL<typecode>MatrixAdd( LALStatus *stat,
                        <datatype>Array *out,
                        <datatype>Array *in1,
                        <datatype>Array *in2 )

void
LAL<typecode>MatrixMultiply( LALStatus *stat,
                             <datatype>Array *out,
                             <datatype>Array *in1,
                             <datatype>Array *in2 )

void
LAL<typecode>MatrixTranspose( LALStatus *stat,
                              <datatype>Array *out,
                              <datatype>Array *in )
\end{verbatim}
\input{MatrixOpsCP}

\idx{LALI2MatrixAdd()}
\idx{LALI4MatrixAdd()}
\idx{LALI8MatrixAdd()}
\idx{LALU2MatrixAdd()}
\idx{LALU4MatrixAdd()}
\idx{LALU8MatrixAdd()}
\idx{LALSMatrixAdd()}
\idx{LALDMatrixAdd()}
\idx{LALCMatrixAdd()}
\idx{LALZMatrixAdd()}

\idx{LALI2MatrixMultiply()}
\idx{LALI4MatrixMultiply()}
\idx{LALI8MatrixMultiply()}
\idx{LALU2MatrixMultiply()}
\idx{LALU4MatrixMultiply()}
\idx{LALU8MatrixMultiply()}
\idx{LALSMatrixMultiply()}
\idx{LALDMatrixMultiply()}
\idx{LALCMatrixMultiply()}
\idx{LALZMatrixMultiply()}

\idx{LALI2MatrixTranspose()}
\idx{LALI4MatrixTranspose()}
\idx{LALI8MatrixTranspose()}
\idx{LALU2MatrixTranspose()}
\idx{LALU4MatrixTranspose()}
\idx{LALU8MatrixTranspose()}
\idx{LALSMatrixTranspose()}
\idx{LALDMatrixTranspose()}
\idx{LALCMatrixTranspose()}
\idx{LALZMatrixTranspose()}

\idx{LALCMatrixAdjoint()}
\idx{LALZMatrixAdjoint()}

\subsubsection*{Description}

The routines \verb@LAL<datatype>MatrixAdd()@ add the matrices
\verb@*in1@ and \verb@*in2@ element-by-element, storing the result in
\verb@*out@.  All of these matrices must have the same dimensionality.
The addition may be performed in-place by pointing \verb@out@ to the
same structure as either \verb@in1@ or \verb@in2@.

The routines \verb@LAL<datatype>MatrixMultiply()@ perform matrix
multiplication, contracting the columns of \verb@*in1@ against the
rows of \verb@*in2@, and storing the result in \verb@*out@.  The
number of columns of \verb@*in1@ must equal the number of rows of
\verb@*in2@, and \verb@*out@ must have the same number of columns as
\verb@*in1@ and the same number of rows as \verb@*in2@.

The routines \verb@LAL<datatype>MatrixTranspose()@ take the transpose
of the matrix \verb@*in@ and store the result in \verb@*out@.  The
number of rows of \verb@*out@ must equal the number of columns of
\verb@*in@, and vice-versa.

The routines \verb@LALCMatrixAdjoint()@ and \verb@LALZMatrixAdjoint()@
take the Hermitian conjugate (adjoint) of the matrix \verb@*in@ and
store the result in \verb@*out@: this involves transposing the matrix
and taking the complex conjugate.  The number of rows of \verb@*out@
must equal the number of columns of \verb@*in@, and vice-versa.

Except for the adjoint routines, the prototype templates above in fact
refer to 10 separate routines each, corresponding to all the numerical
atomic datatypes \verb@<datatype>@ referred to by \verb@<typecode>@:
\begin{center}
\begin{tabular}{|c@{\qquad}c|c@{\qquad}c|}
\hline
\tt <typecode> & \tt <datatype> & \tt <typecode> & \tt <datatype> \\
\hline
\tt I2 & \tt  INT2 & \tt U2 & \tt    UINT2  \\
\tt I4 & \tt  INT4 & \tt U4 & \tt    UINT4  \\
\tt I8 & \tt  INT8 & \tt U8 & \tt    UINT8  \\
\tt  S & \tt REAL4 & \tt  C & \tt COMPLEX8  \\
\tt  D & \tt REAL8 & \tt  Z & \tt COMPLEX16 \\
\hline
\end{tabular}
\end{center}

\subsubsection*{Algorithm}

Matrix addition is simply carried through element-by-element.  It
involves one addition operation per element of the output matrix.

The matrix product $\mathsf{Z}^a{}_b=\mathsf{X}^a{}_c
\mathsf{Y}^c{}_b$ of two matrices $\mathsf{X}^a{}_c$ and
$\mathsf{Y}^c{}_b$ is given by the element formula
$Z^i{}_j=\sum_{k=1}^N X^i{}_k Y^k{}_j$, where $N$ is the number of
columns of $\mathsf{X}^a{}_c$ \emph{and} the number of rows of
$\mathsf{Y}^c{}_b$.  This can also be used to compute the inner
product of two vectors $\mathsf{x}^a$ and $\mathsf{y}^a$: simply store
the transpose $(\mathsf{x}^T)_a$ as a row vector (single-row matrix)
as the first operand, and $\mathsf{y}^a$ as a column vector
(single-column matrix) as the second operand.  To compute the vector
outer product, simply transpose the second argument rather than the
first.  These computations involve $N$ additions and multiplications
per element of the output matrix.

The transpose $(\mathsf{X}^T){}^a{}_b$ of a matrix $\mathsf{X}^a{}_b$
is given by $(X^T){}^i{}_j=X^j{}_i$.  The adjoint
$(\mathsf{X}^\dag){}^a{}_b$ of a complex matrix $\mathsf{X}^a{}_b$ is
given by $(X^\dag){}^i{}_j=X^j{}_i{}^*$, where ${}^*$ denotes complex
conjugation.  Transposition involves no arithmetic operations, just
one assignment per element of the output.  Conjugation involves one
multiplication (negating the sign of the imaginary part) per element.

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{MatrixOpsCV}}

******************************************************* </lalLaTeX> */

#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/MatrixUtils.h>

NRCSID( MATRIXOPSC, "$Id$" );

define(`TYPECODE',`I2')dnl
include(`LALMatrixOps.m4')dnl

define(`TYPECODE',`I4')dnl
include(`LALMatrixOps.m4')dnl

define(`TYPECODE',`I8')dnl
include(`LALMatrixOps.m4')dnl

define(`TYPECODE',`U2')dnl
include(`LALMatrixOps.m4')dnl

define(`TYPECODE',`U4')dnl
include(`LALMatrixOps.m4')dnl

define(`TYPECODE',`U8')dnl
include(`LALMatrixOps.m4')dnl

define(`TYPECODE',`S')dnl
include(`LALMatrixOps.m4')dnl

define(`TYPECODE',`D')dnl
include(`LALMatrixOps.m4')dnl

define(`TYPECODE',`C')dnl
include(`LALMatrixOps.m4')dnl

define(`TYPECODE',`Z')dnl
include(`LALMatrixOps.m4')dnl

/* <lalVerbatim file="MatrixOpsCP"> */
void
LALCMatrixAdjoint( LALStatus *stat, COMPLEX8Array *out, COMPLEX8Array *in1 )
{ /* </lalVerbatim> */
  UINT4 n;        /* number of elements */
  COMPLEX8 *data; /* pointer to elements */

  INITSTATUS( stat, "LALCMatrixAdjoint", MATRIXOPSC );
  ATTATCHSTATUSPTR( stat );

  /* All argument checking is done by the subroutine. */
  TRY( LALCMatrixTranspose( stat->statusPtr, out, in1 ), stat );

  /* We just need to conjugate the result. */
  n = out->dimLength->data[0]*out->dimLength->data[1];
  data = out->data;
  while ( n-- )
    (data++)->im *= -1.0;

  /* That's all. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="MatrixOpsCP"> */
void
LALZMatrixAdjoint( LALStatus *stat, COMPLEX16Array *out, COMPLEX16Array *in1 )
{ /* </lalVerbatim> */
  UINT4 n;         /* number of elements */
  COMPLEX16 *data; /* pointer to elements */

  INITSTATUS( stat, "LALZMatrixAdjoint", MATRIXOPSC );
  ATTATCHSTATUSPTR( stat );

  /* All argument checking is done by the subroutine. */
  TRY( LALZMatrixTranspose( stat->statusPtr, out, in1 ), stat );

  /* We just need to conjugate the result. */
  n = out->dimLength->data[0]*out->dimLength->data[1];
  data = out->data;
  while ( n-- )
    (data++)->im *= -1.0;

  /* That's all. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
