dnl $Id$
ifelse(TYPECODE,`Z',`define(`TYPE',`COMPLEX16')define(`SIZE',`8')')dnl
ifelse(TYPECODE,`C',`define(`TYPE',`COMPLEX8')define(`SIZE',`4')')dnl
ifelse(TYPECODE,`D',`define(`TYPE',`REAL8')define(`SIZE',`8')')dnl
ifelse(TYPECODE,`S',`define(`TYPE',`REAL4')define(`SIZE',`4')')dnl
ifelse(TYPECODE,`I2',`define(`TYPE',`INT2')define(`SIZE',`2')')dnl
ifelse(TYPECODE,`I4',`define(`TYPE',`INT4')define(`SIZE',`4')')dnl
ifelse(TYPECODE,`I8',`define(`TYPE',`INT8')define(`SIZE',`8')')dnl
ifelse(TYPECODE,`U2',`define(`TYPE',`UINT2')define(`SIZE',`2')')dnl
ifelse(TYPECODE,`U4',`define(`TYPE',`UINT4')define(`SIZE',`4')')dnl
ifelse(TYPECODE,`U8',`define(`TYPE',`UINT8')define(`SIZE',`8')')dnl
define(`COMPLEX',`0')dnl
ifelse(TYPECODE,`Z',`define(`COMPLEX',`1')')dnl
ifelse(TYPECODE,`C',`define(`COMPLEX',`1')')dnl
define(`ATYPE',`format(`%sArray',TYPE)')dnl
define(`ADDFUNC',`format(`LAL%sMatrixAdd',TYPECODE)')dnl
define(`MULFUNC',`format(`LAL%sMatrixMultiply',TYPECODE)')dnl
define(`TRNFUNC',`format(`LAL%sMatrixTranspose',TYPECODE)')dnl
dnl
void
ADDFUNC ( LALStatus *stat, ATYPE *out, ATYPE *in1, ATYPE *in2 )
{
  UINT4 ni, nj;                      /* index ranges */
  TYPE *outData, *in1Data, *in2Data; /* data pointers */

  INITSTATUS( stat, "ADDFUNC", MATRIXOPSC );

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
  ASSERT( in1->dimLength->data[0] == ni, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  ASSERT( in1->dimLength->data[1] == nj, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  ASSERT( in2->dimLength->data[0] == ni, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );
  ASSERT( in2->dimLength->data[1] == nj, stat, MATRIXUTILSH_EDIM,
	  MATRIXUTILSH_MSGEDIM );

  /* Do addition. */
  ni *= nj;
  outData = out->data;
  in1Data = in1->data;
  in2Data = in2->data;
  while ( ni-- ) {
#if COMPLEX
    outData->re = in1Data->re + in2Data->re;
    (outData++)->im = (in1Data++)->im + (in2Data++)->im;
#else
    *(outData++) = *(in1Data++) + *(in2Data++);
#endif
  }
  RETURN( stat );
}


void
MULFUNC ( LALStatus *stat, ATYPE *out, ATYPE *in1, ATYPE *in2 )
{
  UINT4 i, j, k, ni, nj, nk;         /* dimension indecies and ranges */
  UINT4 ij, ik, kj, in;              /* array indecies */
  TYPE *outData, *in1Data, *in2Data; /* data pointers */

  INITSTATUS( stat, "MULFUNC", MATRIXOPSC );

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
  for ( i = 0, in = 0; i < ni; i++, in += ni ) {
    for ( j = 0, ij = in; j < nj; j++, ij++ ) {
#if COMPLEX
      outData[ij].re = outData[ij].im = 0.0;
#else
      outData[ij] = 0.0;
#endif
      for ( k = 0, ik = in, kj = j; k < nk; k++, ik++, kj += nj ) {
#if COMPLEX
	outData[ij].re += in1Data[ik].re*in2Data[kj].re
	  - in1Data[ik].im*in2Data[kj].im;
	outData[ij].im += in1Data[ik].im*in2Data[kj].re
	  + in1Data[ik].re*in2Data[kj].im;
#else
	outData[ij] += in1Data[ik]*in2Data[kj];
#endif
      }
    }
  }
  RETURN( stat );
}


void
TRNFUNC ( LALStatus *stat, ATYPE *out, ATYPE *in1 )
{
  UINT4 i, j, ni, nj;      /* dimension indecies and ranges */
  UINT4 ij, ji, in;        /* array indecies */
  TYPE *outData, *in1Data; /* data pointers */

  INITSTATUS( stat, "TRNFUNC", MATRIXOPSC );

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
  for ( i = 0, in = 0; i < ni; i++, in += ni ) {
    for ( j = 0, ij = in, ji = i; j < nj; j++, ij++, ji += ni ) {
#if COMPLEX
      outData[ij].re = in1Data[ji].re;
      outData[ij].im = in1Data[ji].im;
#else
      outData[ij] = in1Data[ji];
#endif
    }
  }
  RETURN( stat );
}
