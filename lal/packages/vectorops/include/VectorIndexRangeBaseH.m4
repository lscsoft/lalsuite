/*
dnl $Id$
ifelse(TYPECODE,`CHAR',`define(`TYPE',`CHAR')')
ifelse(TYPECODE,`I2',`define(`TYPE',`INT2')')
ifelse(TYPECODE,`I4',`define(`TYPE',`INT4')')
ifelse(TYPECODE,`I8',`define(`TYPE',`INT8')')
ifelse(TYPECODE,`U2',`define(`TYPE',`UINT2')')
ifelse(TYPECODE,`U4',`define(`TYPE',`UINT4')')
ifelse(TYPECODE,`U8',`define(`TYPE',`UINT8')')
ifelse(TYPECODE,`S',`define(`TYPE',`REAL4')')
ifelse(TYPECODE,`D',`define(`TYPE',`REAL8')')
ifelse(TYPECODE,`C',`define(`TYPE',`COMPLEX8')')
ifelse(TYPECODE,`Z',`define(`TYPE',`COMPLEX16')')
define(`VTYPE',`format(`%sVector',TYPE)')
define(`F1',`format(`LAL%sVectorIndexRange',TYPECODE)')
define(`VPAIRTYPE',`format(`%sVectorPair',TYPE)')
define(`F2',`format(`LAL%sVectorIndexHole',TYPECODE)')
*/


void F1 ( LALStatus *status,
          VTYPE    **result,
          VTYPE     *v,
          const UINT4Vector *indexVector );

void F2 ( LALStatus  *status,
          VPAIRTYPE  *result_pair,
          VTYPE      *v,
          const UINT4Vector *indexVector );

