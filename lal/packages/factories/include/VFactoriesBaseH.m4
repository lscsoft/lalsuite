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
ifelse(TYPECODE,`',`define(`TYPE',`REAL4')')
define(`VTYPE',`format(`%sVector',TYPE)')
define(`CREATEVECTOR',`format(`LAL%sCreateVector',TYPECODE)')
define(`RESIZEVECTOR',`format(`LAL%sResizeVector',TYPECODE)')
define(`DESTROYVECTOR',`format(`LAL%sDestroyVector',TYPECODE)')
ifelse(TYPECODE,`',`define(`XCREATEVECTOR',`XLALCreateVector')',`define(`XCREATEVECTOR',`format(`XLALCreate%s',VTYPE)')')
ifelse(TYPECODE,`',`define(`XRESIZEVECTOR',`XLALResizeVector')',`define(`XRESIZEVECTOR',`format(`XLALResize%s',VTYPE)')')
ifelse(TYPECODE,`',`define(`XDESTROYVECTOR',`XLALDestroyVector')',`define(`XDESTROYVECTOR',`format(`XLALDestroy%s',VTYPE)')')
 * TYPE
 */
VTYPE * XCREATEVECTOR ( UINT4 length );
VTYPE * XRESIZEVECTOR ( VTYPE * vector, UINT4 length );
void XDESTROYVECTOR ( VTYPE * vector );
void CREATEVECTOR ( LALStatus *, VTYPE **, UINT4 );
void RESIZEVECTOR ( LALStatus *, VTYPE **, UINT4 );
void DESTROYVECTOR ( LALStatus *, VTYPE ** );
