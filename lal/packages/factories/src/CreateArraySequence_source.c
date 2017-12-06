#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define STYPE CONCAT2(TYPE,ArraySequence)

#ifdef TYPECODE
#define FUNC CONCAT3(LAL,TYPECODE,CreateArraySequence)
#else
#define FUNC LALCreateArraySequence
#endif


void FUNC ( LALStatus *status, STYPE **aseq, CreateArraySequenceIn *in )
{
  UINT4 i;

  /*
   * Initialize status
   */
  INITSTATUS(status);
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

#undef STYPE
#undef FUNC
