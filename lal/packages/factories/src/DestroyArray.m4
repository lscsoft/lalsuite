dnl $Id$
ifelse(TYPECODE,`Z',`define(`TYPE',`COMPLEX16')')
ifelse(TYPECODE,`C',`define(`TYPE',`COMPLEX8')')
ifelse(TYPECODE,`D',`define(`TYPE',`REAL8')')
ifelse(TYPECODE,`S',`define(`TYPE',`REAL4')')
ifelse(TYPECODE,`I2',`define(`TYPE',`INT2')')
ifelse(TYPECODE,`I4',`define(`TYPE',`INT4')')
ifelse(TYPECODE,`I8',`define(`TYPE',`INT8')')
ifelse(TYPECODE,`U2',`define(`TYPE',`UINT2')')
ifelse(TYPECODE,`U4',`define(`TYPE',`UINT4')')
ifelse(TYPECODE,`U8',`define(`TYPE',`UINT8')')
ifelse(TYPECODE,`',`define(`TYPE',`REAL4')')
define(`ATYPE',`format(`%sArray',TYPE)')
define(`FUNC',`format(`LAL%sDestroyArray',TYPECODE)')
ifelse( TYPECODE, `', `define(`XFUNC',`XLALDestroyArray')', `define(`XFUNC',`format(`XLALDestroy%s',ATYPE)')' )

void XFUNC ( ATYPE *array )
{
  if ( ! array )
    XLAL_ERROR_VOID( "XFUNC", XLAL_EFAULT );
  if ( ! array->dimLength
      || ! array->dimLength->length
      || ! array->dimLength->data
      || ! array->data )
    XLAL_ERROR_VOID( "XFUNC", XLAL_EINVAL );
  XLALDestroyUINT4Vector( array->dimLength );
  LALFree( array->data );
  LALFree( array );
  return;
}

/* <lalVerbatim file="ArrayFactoriesD"> */
void FUNC ( LALStatus *status, ATYPE **array )
{ /* </lalVerbatim> */
  INITSTATUS (status, "FUNC", ARRAYFACTORIESC);

  ASSERT (array,          status, AVFACTORIESH_EVPTR, AVFACTORIESH_MSGEVPTR);
  ASSERT (*array,         status, AVFACTORIESH_EUPTR, AVFACTORIESH_MSGEUPTR);
  ASSERT ((*array)->data, status, AVFACTORIESH_EDPTR, AVFACTORIESH_MSGEDPTR);

  /* Free allocated storage */

  XFUNC ( *array );
  if ( xlalErrno )
  {
    int code = xlalErrno;
    XLALClearErrno();
    if ( code == XLAL_EFAULT )
    {
      ABORT (status, AVFACTORIESH_EUPTR, AVFACTORIESH_MSGEUPTR);
    }
    if ( code == XLAL_EINVAL )
    {
      ABORT (status, AVFACTORIESH_EDPTR, AVFACTORIESH_MSGEDPTR);
    }
  }
  *array = NULL;	    /* make sure we don't point to freed struct */

  RETURN (status);
}
