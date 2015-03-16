#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define ATYPE CONCAT2(TYPE,Array)
#define ADDFUNC CONCAT3(LAL,TYPECODE,MatrixAdd)
#define MULFUNC CONCAT3(LAL,TYPECODE,MatrixMultiply)
#define TRNFUNC CONCAT3(LAL,TYPECODE,MatrixTranspose)

void
ADDFUNC ( LALStatus *stat, ATYPE *out, ATYPE *in1, ATYPE *in2 )
{
  UINT4 ni, nj;                      /* index ranges */
  TYPE *outData, *in1Data, *in2Data; /* data pointers */

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
    *(outData++) = *(in1Data++) + *(in2Data++);
  }
  RETURN( stat );
}


void
MULFUNC ( LALStatus *stat, ATYPE *out, ATYPE *in1, ATYPE *in2 )
{
  UINT4 i, j, k, ni, nj, nk;         /* dimension indices and ranges */
  UINT4 ij, ik, kj, in, kn;          /* array indices */
  TYPE *outData, *in1Data, *in2Data; /* data pointers */

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
TRNFUNC ( LALStatus *stat, ATYPE *out, ATYPE *in1 )
{
  UINT4 i, j, ni, nj;      /* dimension indices and ranges */
  UINT4 ij, ji, in;        /* array indices */
  TYPE *outData, *in1Data; /* data pointers */

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
