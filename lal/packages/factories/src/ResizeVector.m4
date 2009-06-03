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
define(`RESIZEVECTOR',`format(`LAL%sResizeVector',TYPECODE)')
ifelse( TYPECODE, `', `define(`XFUNC',`XLALResizeVector')', `define(`XFUNC',`format(`XLALResize%s',VTYPE)')' )
define(`XCFUNC',`format(`XLALCreate%s',VTYPE)')
define(`XDFUNC',`format(`XLALDestroy%s',VTYPE)')

VTYPE * XFUNC ( VTYPE * vector, UINT4 length )
{
  if ( ! vector )
    return XCFUNC ( length );
  if ( ! length )
  {
    XDFUNC ( vector );
    return NULL;
  }
  vector->data = LALRealloc( vector->data, length * sizeof( *vector->data ) );
  if ( ! vector->data )
  {
    vector->length = 0;
    XLAL_ERROR_NULL( "XFUNC", XLAL_ENOMEM );
  }
  vector->length = length;
  return vector;
}


/* <lalVerbatim file="VectorFactoriesD"> */
void RESIZEVECTOR ( LALStatus *status, VTYPE **vector, UINT4 length )
{ /* </lalVerbatim> */
  /*
   * Initialize status structure
   */
  INITSTATUS( status, "RESIZEVECTOR", VECTORFACTORIESC );

  ASSERT ( vector != NULL, status, AVFACTORIESH_EVPTR, AVFACTORIESH_MSGEVPTR );
  ASSERT ( ! *vector || (*vector)->length, status, AVFACTORIESH_ELENGTH, AVFACTORIESH_MSGELENGTH );
  ASSERT ( length || *vector, status, AVFACTORIESH_ELENGTH, AVFACTORIESH_MSGELENGTH );

  /* Want this to behave like realloc(3), i.e.
   * *vector == NULL => create a new vector
   * length == 0 => destroy the vector
   * otherwise => resize given vector
   */

  *vector = XFUNC ( *vector, length );
  if ( xlalErrno )
  {
    int code = xlalErrno;
    XLALClearErrno();
    if ( code == XLAL_EBADLEN )
    {
      ABORT( status, AVFACTORIESH_ELENGTH, AVFACTORIESH_MSGELENGTH );
    }
    if ( code == XLAL_ENOMEM )
    {
      ABORT( status, AVFACTORIESH_EMALLOC, AVFACTORIESH_MSGEMALLOC );
    }
  }

  RETURN( status );
}
