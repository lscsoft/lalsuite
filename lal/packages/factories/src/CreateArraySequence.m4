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
define(`STYPE',`format(`%sArraySequence',TYPE)')
define(`FUNC',`format(`LAL%sCreateArraySequence',TYPECODE)')

/* <lalVerbatim file="ArraySequenceFactoriesD"> */
void FUNC ( LALStatus *status, STYPE **aseq, CreateArraySequenceIn *in )
{ /* </lalVerbatim> */
  UINT4 i;

  /*
   * Initialize status
   */
  INITSTATUS( status, "FUNC", ARRAYSEQUENCEFACTORIESC );
  ATTATCHSTATUSPTR( status );

  /* Check input structure: report if NULL */

  ASSERT (in != NULL, status, SEQFACTORIESH_EINPTR, SEQFACTORIESH_MSGEINPTR);
  ASSERT (in->dimLength, status,
          SEQFACTORIESH_EINPTR, SEQFACTORIESH_MSGEINPTR);

  /* Check sequence length: report error if 0
   * Use of unsigned for length means we can't check if negative
   * length was passed
   */

  ASSERT (in->length > 0, status,
          SEQFACTORIESH_ESLENGTH, SEQFACTORIESH_MSGESLENGTH);

  /* Check dimension lengths: report error any are 0
   * Use of unsigned for length means we can't check if negative
   * length was passed
   */
#ifndef NDEBUG
  for ( i = 0; i < in->dimLength->length; i++ )
  {
    ASSERT (in->dimLength->data[i] > 0, status,
	    SEQFACTORIESH_EALENGTH, SEQFACTORIESH_MSGEALENGTH);
  }
#endif

  /*
   * Check return structure: If return pointer does not point to a
   *    valid pointer then report an error
   */

  ASSERT (aseq != NULL, status, SEQFACTORIESH_EVPTR, SEQFACTORIESH_MSGEVPTR);
  ASSERT (*aseq == NULL, status, SEQFACTORIESH_EUPTR, SEQFACTORIESH_MSGEUPTR);

  /*
   * Allocate pointer
   */

  *aseq = ( STYPE * ) LALMalloc( sizeof( STYPE ) );
  if ( NULL == *aseq )
  {
    ABORT( status, SEQFACTORIESH_EMALLOC, SEQFACTORIESH_MSGEMALLOC );
  }

  (*aseq)->length = in->length;
  (*aseq)->arrayDim = 1;
  (*aseq)->dimLength = NULL;	/* NULL dimLength until allocated */
  (*aseq)->data   = NULL;	/* NULL data until allocated */

  /*
   * Allocate dimLength
   */
  {
    LALU4CreateVector( status->statusPtr, &((*aseq)->dimLength),
		       in->dimLength->length );
    BEGINFAIL( status ) {
      LALFree ((void *) *aseq);
      ABORT (status, SEQFACTORIESH_EMALLOC, SEQFACTORIESH_MSGEMALLOC);
    } ENDFAIL( status );
    for ( i = 0; i < in->dimLength->length; i++ )
      (*aseq)->arrayDim *= (*aseq)->dimLength->data[i]
	= in->dimLength->data[i];
  }

  /*
   * Allocate storage
   */
  {
    size_t tlength;
    tlength = (*aseq)->length * (*aseq)->arrayDim * sizeof( TYPE );
    (*aseq)->data = ( TYPE * ) LALMalloc (tlength);
  }

  if (NULL == (*aseq)->data)
  {
    /* Must free storage pointed to by *aseq */
    TRY( LALU4DestroyVector( status->statusPtr, &((*aseq)->dimLength) ),
	 status );
    LALFree ((void *) *aseq);
    ABORT (status, SEQFACTORIESH_EMALLOC, SEQFACTORIESH_MSGEMALLOC);
  }

  /* We be done: Normal exit */

  DETATCHSTATUSPTR( status );
  RETURN (status);
}
