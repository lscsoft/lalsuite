/**
 * \defgroup MatrixOps_c Module MatrixOps.c
 * \ingroup MatrixUtils_h
 * \author Creighton, T. D.
 *
 * \brief Routines to perform basic matrix operations.
 *
 * ### Prototypes ###
 *
 * \code
 * void
 * LAL<typecode>MatrixAdd( LALStatus *stat,
 * <datatype>Array *out,
 * <datatype>Array *in1,
 * <datatype>Array *in2 )
 *
 * void
 * LAL<typecode>MatrixMultiply( LALStatus *stat,
 * <datatype>Array *out,
 * <datatype>Array *in1,
 * <datatype>Array *in2 )
 *
 * void
 * LAL<typecode>MatrixTranspose( LALStatus *stat,
 * <datatype>Array *out,
 * <datatype>Array *in )
 * \endcode
 *
 * ### Description ###
 *
 * The routines <tt>LAL\<datatype\>MatrixAdd()</tt> add the matrices
 * <tt>*in1</tt> and <tt>*in2</tt> element-by-element, storing the result in
 * <tt>*out</tt>.  All of these matrices must have the same dimensionality.
 * The addition may be performed in-place by pointing \c out to the
 * same structure as either \c in1 or \c in2.
 *
 * The routines <tt>LAL\<datatype\>MatrixMultiply()</tt> perform matrix
 * multiplication, contracting the columns of <tt>*in1</tt> against the
 * rows of <tt>*in2</tt>, and storing the result in <tt>*out</tt>.  The
 * number of columns of <tt>*in1</tt> must equal the number of rows of
 * <tt>*in2</tt>, and <tt>*out</tt> must have the same number of columns as
 * <tt>*in1</tt> and the same number of rows as <tt>*in2</tt>.
 *
 * The routines <tt>LAL\<datatype\>MatrixTranspose()</tt> take the transpose
 * of the matrix <tt>*in</tt> and store the result in <tt>*out</tt>.  The
 * number of rows of <tt>*out</tt> must equal the number of columns of
 * <tt>*in</tt>, and vice-versa.
 *
 * The routines <tt>LALCMatrixAdjoint()</tt> and <tt>LALZMatrixAdjoint()</tt>
 * take the Hermitian conjugate (adjoint) of the matrix <tt>*in</tt> and
 * store the result in <tt>*out</tt>: this involves transposing the matrix
 * and taking the complex conjugate.  The number of rows of <tt>*out</tt>
 * must equal the number of columns of <tt>*in</tt>, and vice-versa.
 *
 * Except for the adjoint routines, the prototype templates above in fact
 * refer to 10 separate routines each, corresponding to all the numerical
 * atomic datatypes <tt>\<datatype\></tt> referred to by <tt>\<typecode\></tt>:
 *
 * <table><tr><th>\<typecode\></th><th>\<datatype\></th><th>\<typecode\></th><th>\<datatype\></th></tr>
 * <tr><td>I2</td><td> INT2</td><td> U2</td><td>   UINT2</td></tr>
 * <tr><td>I4</td><td> INT4</td><td> U4</td><td>   UINT4</td></tr>
 * <tr><td>I8</td><td> INT8</td><td> U8</td><td>   UINT8</td></tr>
 * <tr><td> S</td><td>REAL4</td><td>  C</td><td>COMPLEX8</td></tr>
 * <tr><td> D</td><td>REAL8</td><td>  Z</td><td>COMPLEX16</td></tr>
 * </table>
 *
 * ### Algorithm ###
 *
 * Matrix addition is simply carried through element-by-element.  It
 * involves one addition operation per element of the output matrix.
 *
 * The matrix product \f$\mathsf{Z}^a{}_b=\mathsf{X}^a{}_c
 * \mathsf{Y}^c{}_b\f$ of two matrices \f$\mathsf{X}^a{}_c\f$ and
 * \f$\mathsf{Y}^c{}_b\f$ is given by the element formula
 * \f$Z^i{}_j=\sum_{k=1}^N X^i{}_k Y^k{}_j\f$, where \f$N\f$ is the number of
 * columns of \f$\mathsf{X}^a{}_c\f$ \e and the number of rows of
 * \f$\mathsf{Y}^c{}_b\f$.  This can also be used to compute the inner
 * product of two vectors \f$\mathsf{x}^a\f$ and \f$\mathsf{y}^a\f$: simply store
 * the transpose \f$(\mathsf{x}^T)_a\f$ as a row vector (single-row matrix)
 * as the first operand, and \f$\mathsf{y}^a\f$ as a column vector
 * (single-column matrix) as the second operand.  To compute the vector
 * outer product, simply transpose the second argument rather than the
 * first.  These computations involve \f$N\f$ additions and multiplications
 * per element of the output matrix.
 *
 * The transpose \f$(\mathsf{X}^T){}^a{}_b\f$ of a matrix \f$\mathsf{X}^a{}_b\f$
 * is given by \f$(X^T){}^i{}_j=X^j{}_i\f$.  The adjoint
 * \f$(\mathsf{X}^\dag){}^a{}_b\f$ of a complex matrix \f$\mathsf{X}^a{}_b\f$ is
 * given by \f$(X^\dag){}^i{}_j=X^j{}_i{}^*\f$, where \f${}^*\f$ denotes complex
 * conjugation.  Transposition involves no arithmetic operations, just
 * one assignment per element of the output.  Conjugation involves one
 * multiplication (negating the sign of the imaginary part) per element.
 *
 */

#include <complex.h>
#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/MatrixUtils.h>

#define TYPECODE I2
#define TYPE INT2
#define SIZE 2
#include "MatrixOps_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE I4
#define TYPE INT4
#define SIZE 4
#include "MatrixOps_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE I8
#define TYPE INT8
#define SIZE 8
#include "MatrixOps_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE U2
#define TYPE UINT2
#define SIZE 2
#include "MatrixOps_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE U4
#define TYPE UINT4
#define SIZE 4
#include "MatrixOps_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE U8
#define TYPE UINT8
#define SIZE 8
#include "MatrixOps_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE S
#define TYPE REAL4
#define SIZE 4
#include "MatrixOps_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE D
#define TYPE REAL8
#define SIZE 8
#include "MatrixOps_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE Z
#define TYPE COMPLEX16
#define SIZE 8
#include "MatrixOps_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE

#define TYPECODE C
#define TYPE COMPLEX8
#define SIZE 4
#include "MatrixOps_source.c"
#undef TYPECODE
#undef TYPE
#undef SIZE


void
LALCMatrixAdjoint( LALStatus *stat, COMPLEX8Array *out, COMPLEX8Array *in1 )
{
  UINT4 n;        /* number of elements */
  COMPLEX8 *data; /* pointer to elements */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* All argument checking is done by the subroutine. */
  TRY( LALCMatrixTranspose( stat->statusPtr, out, in1 ), stat );

  /* We just need to conjugate the result. */
  n = out->dimLength->data[0]*out->dimLength->data[1];
  data = out->data;
  while ( n-- ) {
    *data = conjf(*data);
    ++data;
  }

  /* That's all. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}



void
LALZMatrixAdjoint( LALStatus *stat, COMPLEX16Array *out, COMPLEX16Array *in1 )
{
  UINT4 n;         /* number of elements */
  COMPLEX16 *data; /* pointer to elements */

  INITSTATUS(stat);
  ATTATCHSTATUSPTR( stat );

  /* All argument checking is done by the subroutine. */
  TRY( LALZMatrixTranspose( stat->statusPtr, out, in1 ), stat );

  /* We just need to conjugate the result. */
  n = out->dimLength->data[0]*out->dimLength->data[1];
  data = out->data;
  while ( n-- ) {
    *data = conj(*data);
    ++data;
  }

  /* That's all. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
