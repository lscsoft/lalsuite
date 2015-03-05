#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define VTYPE CONCAT2(TYPE,Vector)

#ifdef TYPECODE
#define FUNC CONCAT3(LAL,TYPECODE,CreateVector)
#define XFUNC CONCAT2(XLALCreate,VTYPE)
#else
#define FUNC LALCreateVector
#define XFUNC XLALCreateVector
#endif

VTYPE * XFUNC ( UINT4 length )
{
  VTYPE * vector;
  vector = LALMalloc( sizeof( *vector ) );
  if ( ! vector )
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  vector->length = length;
  if ( ! length ) /* zero length: set data pointer to be NULL */
    vector->data = NULL;
  else /* non-zero length: allocate memory for data */
  {
#ifdef USE_ALIGNED_MEMORY_ROUTINES
    vector->data = XLALMallocAligned( length * sizeof( *vector->data ) );
#else
    vector->data = LALMalloc( length * sizeof( *vector->data ) );
#endif
    if ( ! vector->data )
    {
      LALFree( vector );
      XLAL_ERROR_NULL( XLAL_ENOMEM );
    }
  }
  return vector;
}


void FUNC ( LALStatus *status, VTYPE **vector, UINT4 length )
{
  /*
   * Initialize status structure
   */

  INITSTATUS(status);

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

#undef VTYPE
#undef FUNC
#undef XFUNC
