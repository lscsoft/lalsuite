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
 * The prototype templates above in fact refer to several separate routines
 * each, corresponding to different numerical atomic datatypes
 * <tt>\<datatype\></tt> referred to by <tt>\<typecode\></tt>:
 *
 * <table><tr><th>\<typecode\></th><th>\<datatype\></th><th>\<typecode\></th><th>\<datatype\></th></tr>
 * <tr><td> D</td><td>REAL8</td><td></td><td></td></tr>
 * </table>
 *
 * ### Algorithm ###
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
 * The transpose \f$(\mathsf{X}^T){}^a{}_b\f$ of a matrix
 * \f$\mathsf{X}^a{}_b\f$ is given by \f$(X^T){}^i{}_j=X^j{}_i\f$.
 * Transposition involves no arithmetic operations, just one assignment per
 * element of the output.
 */

#include <lal/LALStdlib.h>
#include <lal/MatrixUtils.h>

void
LALDMatrixMultiply ( LALStatus *stat, REAL8Array *out, REAL8Array *in1, REAL8Array *in2 )
{
  UINT4 i, j, k, ni, nj, nk;         /* dimension indices and ranges */
  UINT4 ij, ik, kj, in, kn;          /* array indices */
  REAL8 *outData, *in1Data, *in2Data; /* data pointers */

  INITSTATUS(stat);

  /* Check for valid input arguments. */
  ASSERT( out, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( out->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( out->dimLength, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( out->dimLength->data, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( in1, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( in1->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( in1->dimLength, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( in1->dimLength->data, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( in2, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( in2->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( in2->dimLength, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( in2->dimLength->data, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( out->dimLength->length == 2, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  ASSERT( in1->dimLength->length == 2, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  ASSERT( in2->dimLength->length == 2, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  ni = out->dimLength->data[0];
  nj = out->dimLength->data[1];
  nk = in1->dimLength->data[1];
  ASSERT( in1->dimLength->data[0] == ni, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  ASSERT( in2->dimLength->data[0] == nk, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  ASSERT( in2->dimLength->data[1] == nj, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );

  /* Do multiplication. */
  outData = out->data;
  in1Data = in1->data;
  in2Data = in2->data;
  for ( i = 0, in = 0, kn = 0; i < ni; i++, in += nj, kn += nk ) {
    for ( j = 0, ij = in; j < nj; j++, ij++ ) {
      outData[ij] = 0.0;
      for ( k = 0, ik = kn, kj = j; k < nk; k++, ik++, kj += nj ) {
	outData[ij] += in1Data[ik]*in2Data[kj];
      }
    }
  }
  RETURN( stat );
}


void
LALDMatrixTranspose ( LALStatus *stat, REAL8Array *out, REAL8Array *in1 )
{
  UINT4 i, j, ni, nj;      /* dimension indices and ranges */
  UINT4 ij, ji, in;        /* array indices */
  REAL8 *outData, *in1Data; /* data pointers */

  INITSTATUS(stat);

  /* Check for valid input arguments. */
  ASSERT( out, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( out->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( out->dimLength, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( out->dimLength->data, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( in1, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( in1->data, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( in1->dimLength, stat, MATRIXUTILSH_ENUL, MATRIXUTILSH_MSGENUL );
  ASSERT( in1->dimLength->data, stat, MATRIXUTILSH_ENUL,
	  MATRIXUTILSH_MSGENUL );
  ASSERT( out->dimLength->length == 2, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  ASSERT( in1->dimLength->length == 2, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  ni = out->dimLength->data[0];
  nj = out->dimLength->data[1];
  ASSERT( in1->dimLength->data[1] == ni, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  ASSERT( in1->dimLength->data[0] == nj, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );

  /* Do transposition. */
  outData = out->data;
  in1Data = in1->data;
  for ( i = 0, in = 0; i < ni; i++, in += nj ) {
    for ( j = 0, ij = in, ji = i; j < nj; j++, ij++, ji += ni ) {
      outData[ij] = in1Data[ji];
    }
  }
  RETURN( stat );
}
