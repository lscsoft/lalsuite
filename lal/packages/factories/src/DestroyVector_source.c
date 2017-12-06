#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define VTYPE CONCAT2(TYPE,Vector)

#ifdef TYPECODE
#define FUNC CONCAT3(LAL,TYPECODE,DestroyVector)
#define XFUNC CONCAT2(XLALDestroy,VTYPE)
#else
#define FUNC LALDestroyVector
#define XFUNC XLALDestroyVector
#endif

void XFUNC ( VTYPE *vector )
{
  if ( ! vector )
    return;
  if ( ( ! vector->length || ! vector->data ) && ( vector->length || vector->data  ) )
    XLAL_ERROR_VOID( XLAL_EINVAL );
#ifdef USE_ALIGNED_MEMORY_ROUTINES
  if ( vector->data )
    XLALFreeAligned( vector->data );
#else
  if ( vector->data )
    XLALFree( vector->data );
#endif
  vector->data = NULL; /* leave length non-zero to detect repeated frees */
  LALFree( vector );
  return;
}


void FUNC ( LALStatus *status, VTYPE **vector )
{
  /*
   * Initialize status
   */

  INITSTATUS(status);

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

  XFUNC ( *vector );
  if ( xlalErrno )
  {
    int code = xlalErrno;
    XLALClearErrno();
    if ( code == XLAL_EFAULT )
    {
      ABORT( status, AVFACTORIESH_EUPTR, AVFACTORIESH_MSGEUPTR );
    }
    if ( code == XLAL_EINVAL )
    {
      ABORT( status, AVFACTORIESH_EDPTR, AVFACTORIESH_MSGEDPTR );
    }
  }

  *vector = NULL;		/* make sure we don't point to freed struct */

  RETURN( status );
}

#undef VTYPE
#undef FUNC
#undef XFUNC
