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

  LALFree( (*vector)->data ); /* free allocated data */
  LALFree( *vector );	/* free vector struct itself */

  *vector = NULL;		/* make sure we don't point to freed struct */

  RETURN( status );
}
