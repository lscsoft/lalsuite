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
ifelse(TYPECODE,`',`define(`TYPE',`REAL4')')
define(`ATYPE',`format(`%sArray',TYPE)')
define(`CREATEARRAY',`format(`LAL%sCreateArray',TYPECODE)')
define(`RESIZEARRAY',`format(`LAL%sResizeArray',TYPECODE)')
define(`DESTROYARRAY',`format(`LAL%sDestroyArray',TYPECODE)')
 * TYPE
 */
void CREATEARRAY ( LALStatus *, ATYPE **, UINT4Vector * );
void RESIZEARRAY ( LALStatus *, ATYPE **, UINT4Vector * );
void DESTROYARRAY ( LALStatus *, ATYPE ** );
