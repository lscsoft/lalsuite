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

/* <lalVerbatim file="VectorFactoriesD"> */
void FUNC ( Status *status, VTYPE **vector, UINT4 length ) 
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

  *vector = ( VTYPE * ) LALMalloc( sizeof( VTYPE ) );
  ASSERT( *vector != NULL, status,
          AVFACTORIESH_EMALLOC, AVFACTORIESH_MSGEMALLOC );

  (*vector)->length = 0;	/* length 0 until storage allocated */
  (*vector)->data   = NULL;	/* NULL data until allocated */

  /* 
   * Allocate storage 
   * Test that storage is properly allocated. Can't handle with ASSERT
   * since we need to de-allocate structure pointer before an error return
   */

  (*vector)->data = ( TYPE * ) LALMalloc( length*sizeof( TYPE ) );

  if ( NULL == (*vector)->data )
  {
    /* Must free storage pointed to by *vector */
    LALFree( *vector );
    ABORT( status, AVFACTORIESH_EMALLOC, AVFACTORIESH_MSGEMALLOC );
  }
  (*vector)->length = length;	/* Set length if storage allocated */

  /* We be done: Normal exit */

  RETURN( status );
}
