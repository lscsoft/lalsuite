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
ifelse(TYPECODE,`',`define(`XCREATEARRAY',`XLALCreateArray')',`define(`XCREATEARRAY',`format(`XLALCreate%s',ATYPE)')')
ifelse(TYPECODE,`',`define(`XCREATEARRAYL',`XLALCreateArrayL')',`define(`XCREATEARRAYL',`format(`XLALCreate%sL',ATYPE)')')
ifelse(TYPECODE,`',`define(`XCREATEARRAYV',`XLALCreateArrayV')',`define(`XCREATEARRAYV',`format(`XLALCreate%sV',ATYPE)')')
ifelse(TYPECODE,`',`define(`XRESIZEARRAY',`XLALResizeArray')',`define(`XRESIZEARRAY',`format(`XLALResize%s',ATYPE)')')
ifelse(TYPECODE,`',`define(`XRESIZEARRAYL',`XLALResizeArrayL')',`define(`XRESIZEARRAYL',`format(`XLALResize%sL',ATYPE)')')
ifelse(TYPECODE,`',`define(`XRESIZEARRAYV',`XLALResizeArrayV')',`define(`XRESIZEARRAYV',`format(`XLALResize%sV',ATYPE)')')
ifelse(TYPECODE,`',`define(`XDESTROYARRAY',`XLALDestroyArray')',`define(`XDESTROYARRAY',`format(`XLALDestroy%s',ATYPE)')')

define(`CREATEARRAY',`format(`LAL%sCreateArray',TYPECODE)')
define(`RESIZEARRAY',`format(`LAL%sResizeArray',TYPECODE)')
define(`DESTROYARRAY',`format(`LAL%sDestroyArray',TYPECODE)')
 * TYPE
 */
ATYPE * XCREATEARRAYL ( UINT4, ... );
ATYPE * XCREATEARRAYV ( UINT4, UINT4 * );
ATYPE * XCREATEARRAY ( UINT4Vector * );
ATYPE * XRESIZEARRAYL ( ATYPE *, UINT4, ... );
ATYPE * XRESIZEARRAYV ( ATYPE *, UINT4, UINT4 * );
ATYPE * XRESIZEARRAY ( ATYPE *, UINT4Vector * );
void XDESTROYARRAY ( ATYPE * );
void CREATEARRAY ( LALStatus *, ATYPE **, UINT4Vector * );
void RESIZEARRAY ( LALStatus *, ATYPE **, UINT4Vector * );
void DESTROYARRAY ( LALStatus *, ATYPE ** );
