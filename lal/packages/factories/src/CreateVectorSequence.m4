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
define(`STYPE',`format(`%sVectorSequence',TYPE)')
define(`FUNC',`format(`LAL%sCreateVectorSequence',TYPECODE)')

/* <lalVerbatim file="VectorSequenceFactoriesD"> */
void FUNC ( LALStatus *status, STYPE **vseq, CreateVectorSequenceIn *in ) 
{ /* </lalVerbatim> */
  /* 
   * Initialize status
   */
  INITSTATUS( status, "FUNC", VECTORSEQUENCEFACTORIESC );	

  /* Check input structure: report if NULL */

  ASSERT (in != NULL, status, SEQFACTORIESH_EINPTR, SEQFACTORIESH_MSGEINPTR);
      
  /* Check sequence length: report error if 0 
   * Use of unsigned for length means we can't check if negative
   * length was passed
   */

  ASSERT (in->length > 0, status,
          SEQFACTORIESH_ESLENGTH, SEQFACTORIESH_MSGESLENGTH);

  /* Check vector length: report error if 0 
   * Use of unsigned for length means we can't check if negative
   * length was passed
   */

  ASSERT (in->vectorLength > 0, status,
          SEQFACTORIESH_EVLENGTH, SEQFACTORIESH_MSGEVLENGTH); 

  /* 
   * Check return structure: If return pointer does not point to a
   *    valid pointer then report an error 
   */

  ASSERT (vseq != NULL, status, SEQFACTORIESH_EVPTR, SEQFACTORIESH_MSGEVPTR);
  ASSERT (*vseq == NULL, status, SEQFACTORIESH_EUPTR, SEQFACTORIESH_MSGEUPTR);

  /*
   * Allocate pointer
   */

  *vseq = ( STYPE * ) LALMalloc( sizeof( STYPE ) );
  if ( NULL == *vseq )
  {
    ABORT( status, SEQFACTORIESH_EMALLOC, SEQFACTORIESH_MSGEMALLOC );
  }

  (*vseq)->length = 0;	/* length 0 until storage allocated */
  (*vseq)->vectorLength = 0; /* vector length 0 until storage allocated */
  (*vseq)->data   = NULL;	/* NULL data until allocated */

  /* 
   * Allocate storage 
   */

  {
    size_t tlength;
    tlength = in->vectorLength * in->length * sizeof( TYPE );
    (*vseq)->data = ( TYPE * ) LALMalloc (tlength);
  }

  if (NULL == (*vseq)->data)
  {
    /* Must free storage pointed to by *vseq */
    LALFree ((void *) *vseq);
    ABORT (status, SEQFACTORIESH_EMALLOC, SEQFACTORIESH_MSGEMALLOC);
  }
 
  /* Set length, vectorLength if storage allocated */

  (*vseq)->length = in->length;	
  (*vseq)->vectorLength = in->vectorLength;

  /* We be done: Normal exit */

  RETURN (status);
}
