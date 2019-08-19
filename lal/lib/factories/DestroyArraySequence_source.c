#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define STYPE CONCAT2(TYPE,ArraySequence)

#ifdef TYPECODE
#define FUNC CONCAT3(LAL,TYPECODE,DestroyArraySequence)
#else
#define FUNC LALDestroyArraySequence
#endif

void FUNC ( LALStatus *status, STYPE **aseq )
{
  /*
   * Initialize status
   */

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  /*
   * Check aseq: is it non-NULL?
   */

  ASSERT (aseq != NULL, status, SEQFACTORIESH_EVPTR, SEQFACTORIESH_MSGEVPTR);

  /*
   * Check aseq: does it point to non-NULL?
   */

  ASSERT (*aseq != NULL,status, SEQFACTORIESH_EUPTR, SEQFACTORIESH_MSGEUPTR);

  /*
   * Check dimLength in aseq: does it point to non-NULL?
   */

  ASSERT ((*aseq)->dimLength != NULL, status,
          SEQFACTORIESH_EDPTR, SEQFACTORIESH_MSGEDPTR);

  /*
   * Check data in aseq: does it point to non-NULL?
   */

  ASSERT ((*aseq)->data != NULL, status,
          SEQFACTORIESH_EDPTR, SEQFACTORIESH_MSGEDPTR);

  /* Ok, now let's free allocated storage */

  TRY( LALU4DestroyVector( status->statusPtr, &((*aseq)->dimLength) ),
       status );
  LALFree ( (*aseq)->data ); /* free allocated data */
  LALFree ( *aseq );	      /* free aseq struct itself */

  *aseq = NULL;		/* make sure we don't point to freed struct */

  DETATCHSTATUSPTR( status );
  RETURN (status);
}

#undef STYPE
#undef FUNC
