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
ifelse(TYPECODE,`CHAR',`define(`TYPE',`CHAR')')
ifelse(TYPECODE,`',`define(`TYPE',`REAL4')')
define(`VTYPE',`format(`%sVector',TYPE)')
define(`FUNC',`format(`LAL%sDestroyVector',TYPECODE)')
ifelse(TYPECODE, `', `define(`XFUNC',`XLALDestroyVector')', `define(`XFUNC',`format(`XLALDestroy%s',VTYPE)')' )

void XFUNC ( VTYPE *vector )
{
  if ( ! vector )
    XLAL_ERROR_VOID( "XFUNC", XLAL_EFAULT );
  if ( ! vector->length || ! vector->data )
    XLAL_ERROR_VOID( "XFUNC", XLAL_EINVAL );
  LALFree( vector->data );
  vector->length = 0;
  vector->data = NULL;
  LALFree( vector );
  return;
}

/* <lalVerbatim file="VectorFactoriesD"> */
void FUNC ( LALStatus *status, VTYPE **vector )
{ /* </lalVerbatim> */
  /* 
   * Initialize status
   */

  INITSTATUS( status, "FUNC", VECTORFACTORIESC );	
      
  /* 
   * Check vector: is it non-NULL?
   */

  ASSERT( vector != NULL, status, AVFACTORIESH_EVPTR, AVFACTORIESH_MSGEVPTR );

  /* 
   * Check vector: does it point to non-NULL?
   */

  ASSERT( *vector != NULL, status, AVFACTORIESH_EUPTR, AVFACTORIESH_MSGEUPTR );

  /*
   * Check data in vector: does it point to non-NULL
   */

  ASSERT( (*vector)->data != NULL, status,
          AVFACTORIESH_EDPTR, AVFACTORIESH_MSGEDPTR );

  /* Ok, now let's free allocated storage */

  XFUNC ( *vector );
  if ( xlalErrno )
  {
    int code = xlalErrno;
    XLALClearErrno();
    if ( code == XLAL_EFAULT )
    {
      ABORT( status, AVFACTORIESH_EUPTR, AVFACTORIESH_MSGEUPTR );
    }
    if ( code == XLAL_EINVAL )
    {
      ABORT( status, AVFACTORIESH_EDPTR, AVFACTORIESH_MSGEDPTR );
    }
  }

  *vector = NULL;		/* make sure we don't point to freed struct */

  RETURN( status );
}
