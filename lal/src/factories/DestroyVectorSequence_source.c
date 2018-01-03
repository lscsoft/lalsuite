#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define STYPE CONCAT2(TYPE,VectorSequence)

#ifdef TYPECODE
#define FUNC CONCAT3(LAL,TYPECODE,DestroyVectorSequence)
#define XFUNC CONCAT2(XLALDestroy,STYPE)
#else
#define FUNC LALDestroyVectorSequence
#define XFUNC XLALDestroyVectorSequence
#endif

void XFUNC ( STYPE *vseq )
{
  if ( ! vseq )
    return;
  if ( ( ! vseq->length || ! vseq->vectorLength ) && vseq->data )
    XLAL_ERROR_VOID( XLAL_EINVAL );
  if ( ! vseq->data && ( vseq->length || vseq->vectorLength ) )
    XLAL_ERROR_VOID( XLAL_EINVAL );
  if ( vseq->data )
    LALFree( vseq->data );
  vseq->data = NULL; /* leave lengths as they are to indicate freed vector */
  LALFree( vseq );
  return;
}


void FUNC ( LALStatus *status, STYPE **vseq )
{
  /*
   * Initialize status
   */

  INITSTATUS(status);

  /*
   * Check vseq: is it non-NULL?
   */

  ASSERT (vseq != NULL, status, SEQFACTORIESH_EVPTR, SEQFACTORIESH_MSGEVPTR);

  /*
   * Check vseq: does it point to non-NULL?
   */

  ASSERT (*vseq != NULL,status, SEQFACTORIESH_EUPTR, SEQFACTORIESH_MSGEUPTR);

  /*
   * Check data in vseq: does it point to non-NULL?
   */

  ASSERT ((*vseq)->data != NULL, status,
          SEQFACTORIESH_EDPTR, SEQFACTORIESH_MSGEDPTR);

  /* Ok, now let's free allocated storage */

  XFUNC ( *vseq );
  if ( xlalErrno )
  {
    int code = xlalErrno;
    XLALClearErrno();
    if ( code == XLAL_EFAULT )
    {
      ABORT (status, SEQFACTORIESH_EUPTR, SEQFACTORIESH_MSGEUPTR);
    }
    if ( code == XLAL_EINVAL )
    {
      ABORT (status, SEQFACTORIESH_EDPTR, SEQFACTORIESH_MSGEDPTR);
    }
  }

  *vseq = NULL;		/* make sure we don't point to freed struct */

  RETURN (status);
}

#undef STYPE
#undef FUNC
#undef XFUNC
