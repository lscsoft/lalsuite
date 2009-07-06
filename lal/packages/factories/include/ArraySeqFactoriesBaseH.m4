/*
dnl $Id$
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
define(`ASEQTYPE',`format(`%sArraySequence',TYPE)')
define(`F5',`format(`LAL%sCreateArraySequence',TYPECODE)')
define(`F6',`format(`LAL%sDestroyArraySequence',TYPECODE)')
*/



void F5 ( LALStatus *status,
          ASEQTYPE **arraySequence,
          CreateArraySequenceIn *aSeqParams);

void F6 ( LALStatus *status,
          ASEQTYPE **arraySeqence);

