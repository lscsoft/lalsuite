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
define(`FUNC',`format(`LAL%sCreateVector',TYPECODE)')
ifelse( TYPECODE, `', `define(`XFUNC',`XLALCreateVector')', `define(`XFUNC',`format(`XLALCreate%s',VTYPE)')' )

VTYPE * XFUNC ( UINT4 length )
{
  VTYPE * vector;
  vector = LALMalloc( sizeof( *vector ) );
  if ( ! vector )
    XLAL_ERROR_NULL( "XFUNC", XLAL_ENOMEM );
  vector->length = length;
  if ( ! length ) /* zero length: set data pointer to be NULL */
    vector->data = NULL;
  else /* non-zero length: allocate memory for data */
  {
    vector->data = LALMalloc( length * sizeof( *vector->data ) );
    if ( ! vector->data )
    {
      LALFree( vector );
      XLAL_ERROR_NULL( "XFUNC", XLAL_ENOMEM );
    }
  }
  return vector;
}

/* <lalVerbatim file="VectorFactoriesD"> */
void FUNC ( LALStatus *status, VTYPE **vector, UINT4 length ) 
{ /* </lalVerbatim> */
  /* 
   * Initialize status structure
   */

  INITSTATUS( status, "FUNC", VECTORFACTORIESC );	
      
  /* Check sequence length: report error if 0 
   * Use of unsigned for length means we can't check if negative
   * length was passed
   */

  ASSERT( length > 0, status, AVFACTORIESH_ELENGTH, AVFACTORIESH_MSGELENGTH );

  /* 
   * Check return structure: If return pointer does not point to a
   *    valid pointer then report an error 
   */

  ASSERT( vector != NULL, status, AVFACTORIESH_EVPTR, AVFACTORIESH_MSGEVPTR );

  ASSERT( *vector == NULL, status, AVFACTORIESH_EUPTR, AVFACTORIESH_MSGEUPTR );

  /*
   * Allocate pointer
   */

  *vector = XFUNC ( length );
  if ( ! *vector )
  {
    XLALClearErrno();
    ABORT( status, AVFACTORIESH_EMALLOC, AVFACTORIESH_MSGEMALLOC );
  }

  /* We be done: Normal exit */

  RETURN( status );
}
