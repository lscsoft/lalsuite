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
define(`CREATEVECTOR',`format(`LAL%sCreateVector',TYPECODE)')
define(`DESTROYVECTOR',`format(`LAL%sDestroyVector',TYPECODE)')

/* <lalVerbatim file="VectorFactoriesD"> */
void RESIZEVECTOR ( LALStatus *status, VTYPE **vector, UINT4 length ) 
{ /* </lalVerbatim> */
  TYPE *p; /* temporary pointer */

  /* 
   * Initialize status structure
   */
  INITSTATUS( status, "RESIZEVECTOR", VECTORFACTORIESC );	

  ASSERT ( vector != NULL, status, AVFACTORIESH_EVPTR, AVFACTORIESH_MSGEVPTR );
      
  /* Want this to behave like realloc(3), i.e.
   * *vector == NULL => create a new vector 
   * length == 0 => destroy the vector 
   * otherwise => resize given vector 
   */
  if ( (*vector) == NULL )
    {
      CREATEVECTOR ( status, vector, length );
    }
  else if ( length == 0 )
    {
      DESTROYVECTOR ( status, vector );
    }
  else
    {
      /* 
       * Reallocate storage 
       * Test that storage is properly allocated. Can't handle with ASSERT
       * since we need to de-allocate structure pointer before an error return
       */

      /* 
         NOTE: LALRealloc() free()s the pointer that's been passed in.  This
         is different from the documented behavior of realloc(3) on Darwin,
         although the realloc(3) on Darwin actually *does* free() the memory.
         Ignore what Darwin's man pages say about realloc() and reallocf() 
       */
      p = LALRealloc( (*vector)->data, length*sizeof( TYPE ));

      if (p == NULL)
        {
          /* FIXME: question -- do I need to free (*vector)->data here? */
          LALFree( (*vector)->data );
          ABORT( status, AVFACTORIESH_EMALLOC, AVFACTORIESH_MSGEMALLOC );
        }
    
      (*vector)->data = p;
      (*vector)->length = length;	/* Set length if storage allocated */
    }

  RETURN( status );
}
