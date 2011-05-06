/**
\author Creighton, T. D.
\file
*/

/**

\heading{Module \ref MatrixOps.c}
\latexonly\label{ss_MatrixOps_c}\endlatexonly

Routines to perform basic matrix operations.

\heading{Prototypes}

\code
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
\endcode






































\heading{Description}

The routines <tt>LAL<datatype>MatrixAdd()</tt> add the matrices
<tt>*in1</tt> and <tt>*in2</tt> element-by-element, storing the result in
<tt>*out</tt>.  All of these matrices must have the same dimensionality.
The addition may be performed in-place by pointing \c out to the
same structure as either \c in1 or \c in2.

The routines <tt>LAL<datatype>MatrixMultiply()</tt> perform matrix
multiplication, contracting the columns of <tt>*in1</tt> against the
rows of <tt>*in2</tt>, and storing the result in <tt>*out</tt>.  The
number of columns of <tt>*in1</tt> must equal the number of rows of
<tt>*in2</tt>, and <tt>*out</tt> must have the same number of columns as
<tt>*in1</tt> and the same number of rows as <tt>*in2</tt>.

The routines <tt>LAL<datatype>MatrixTranspose()</tt> take the transpose
of the matrix <tt>*in</tt> and store the result in <tt>*out</tt>.  The
number of rows of <tt>*out</tt> must equal the number of columns of
<tt>*in</tt>, and vice-versa.

The routines <tt>LALCMatrixAdjoint()</tt> and <tt>LALZMatrixAdjoint()</tt>
take the Hermitian conjugate (adjoint) of the matrix <tt>*in</tt> and
store the result in <tt>*out</tt>: this involves transposing the matrix
and taking the complex conjugate.  The number of rows of <tt>*out</tt>
must equal the number of columns of <tt>*in</tt>, and vice-versa.

Except for the adjoint routines, the prototype templates above in fact
refer to 10 separate routines each, corresponding to all the numerical
atomic datatypes <tt><datatype></tt> referred to by <tt><typecode></tt>:

<table><tr><td>

\tt <typecode></td><td>\tt <datatype></td><td>\tt <typecode></td><td>\tt <datatype></td></tr>
<tr><td>
\tt I2</td><td>\tt  INT2</td><td>\tt U2</td><td>\tt    UINT2</td></tr>
<tr><td>\tt I4</td><td>\tt  INT4</td><td>\tt U4</td><td>\tt    UINT4</td></tr>
<tr><td>\tt I8</td><td>\tt  INT8</td><td>\tt U8</td><td>\tt    UINT8</td></tr>
<tr><td>\tt  S</td><td>\tt REAL4</td><td>\tt  C</td><td>\tt COMPLEX8</td></tr>
<tr><td>\tt  D</td><td>\tt REAL8</td><td>\tt  Z</td><td>\tt COMPLEX16</td></tr>
<tr><td>
</td></tr></table>


\heading{Algorithm}

Matrix addition is simply carried through element-by-element.  It
involves one addition operation per element of the output matrix.

The matrix product \f$\mathsf{Z}^a{}_b=\mathsf{X}^a{}_c
\mathsf{Y}^c{}_b\f$ of two matrices \f$\mathsf{X}^a{}_c\f$ and
\f$\mathsf{Y}^c{}_b\f$ is given by the element formula
\f$Z^i{}_j=\sum_{k=1}^N X^i{}_k Y^k{}_j\f$, where \f$N\f$ is the number of
columns of \f$\mathsf{X}^a{}_c\f$ \e and the number of rows of
\f$\mathsf{Y}^c{}_b\f$.  This can also be used to compute the inner
product of two vectors \f$\mathsf{x}^a\f$ and \f$\mathsf{y}^a\f$: simply store
the transpose \f$(\mathsf{x}^T)_a\f$ as a row vector (single-row matrix)
as the first operand, and \f$\mathsf{y}^a\f$ as a column vector
(single-column matrix) as the second operand.  To compute the vector
outer product, simply transpose the second argument rather than the
first.  These computations involve \f$N\f$ additions and multiplications
per element of the output matrix.

The transpose \f$(\mathsf{X}^T){}^a{}_b\f$ of a matrix \f$\mathsf{X}^a{}_b\f$
is given by \f$(X^T){}^i{}_j=X^j{}_i\f$.  The adjoint
\f$(\mathsf{X}^\dag){}^a{}_b\f$ of a complex matrix \f$\mathsf{X}^a{}_b\f$ is
given by \f$(X^\dag){}^i{}_j=X^j{}_i{}^*\f$, where \f${}^*\f$ denotes complex
conjugation.  Transposition involves no arithmetic operations, just
one assignment per element of the output.  Conjugation involves one
multiplication (negating the sign of the imaginary part) per element.

\heading{Uses}

\heading{Notes}



*/

#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/MatrixUtils.h>

NRCSID( MATRIXOPSC, "$Id$" );

#define TYPECODE I2
#define TYPE INT2
#define SIZE 2
#define COMPLEX 0
#include "MatrixOps_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE
#undef COMPLEX

#define TYPECODE I4
#define TYPE INT4
#define SIZE 4
#define COMPLEX 0
#include "MatrixOps_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE
#undef COMPLEX

#define TYPECODE I8
#define TYPE INT8
#define SIZE 8
#define COMPLEX 0
#include "MatrixOps_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE
#undef COMPLEX

#define TYPECODE U2
#define TYPE UINT2
#define SIZE 2
#define COMPLEX 0
#include "MatrixOps_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE
#undef COMPLEX

#define TYPECODE U4
#define TYPE UINT4
#define SIZE 4
#define COMPLEX 0
#include "MatrixOps_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE
#undef COMPLEX

#define TYPECODE U8
#define TYPE UINT8
#define SIZE 8
#define COMPLEX 0
#include "MatrixOps_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE
#undef COMPLEX

#define TYPECODE S
#define TYPE REAL4
#define SIZE 4
#define COMPLEX 0
#include "MatrixOps_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE
#undef COMPLEX

#define TYPECODE D
#define TYPE REAL8
#define SIZE 8
#define COMPLEX 0
#include "MatrixOps_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE
#undef COMPLEX

#define TYPECODE Z
#define TYPE COMPLEX16
#define SIZE 8
#define COMPLEX 1
#include "MatrixOps_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE
#undef COMPLEX

#define TYPECODE C
#define TYPE COMPLEX8
#define SIZE 4
#define COMPLEX 1
#include "MatrixOps_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE
#undef COMPLEX


void
LALCMatrixAdjoint( LALStatus *stat, COMPLEX8Array *out, COMPLEX8Array *in1 )
{ 
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



void
LALZMatrixAdjoint( LALStatus *stat, COMPLEX16Array *out, COMPLEX16Array *in1 )
{ 
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
