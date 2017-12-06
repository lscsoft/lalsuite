#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define VTYPE CONCAT2(TYPE,Vector)

#ifdef TYPECODE
#define RESIZEVECTOR CONCAT3(LAL,TYPECODE,ResizeVector)
#define XFUNC CONCAT2(XLALResize,VTYPE)
#else
#define RESIZEVECTOR LALResizeVector
#define XFUNC XLALResizeVector
#endif

#define XCFUNC CONCAT2(XLALCreate,VTYPE)
#define XDFUNC CONCAT2(XLALDestroy,VTYPE)


VTYPE * XFUNC ( VTYPE * vector, UINT4 length )
{
  if ( ! vector )
    return XCFUNC ( length );
  if ( ! length )
  {
    XDFUNC ( vector );
    return NULL;
  }
#ifdef USE_ALIGNED_MEMORY_ROUTINES
  vector->data = XLALReallocAligned( vector->data, length * sizeof( *vector->data ) );
#else
  vector->data = LALRealloc( vector->data, length * sizeof( *vector->data ) );
#endif
  if ( ! vector->data )
  {
    vector->length = 0;
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }
  vector->length = length;
  return vector;
}



void RESIZEVECTOR ( LALStatus *status, VTYPE **vector, UINT4 length )
{
  /*
   * Initialize status structure
   */
  INITSTATUS(status);

  ASSERT ( vector != NULL, status, AVFACTORIESH_EVPTR, AVFACTORIESH_MSGEVPTR );
  ASSERT ( ! *vector || (*vector)->length, status, AVFACTORIESH_ELENGTH, AVFACTORIESH_MSGELENGTH );
  ASSERT ( length || *vector, status, AVFACTORIESH_ELENGTH, AVFACTORIESH_MSGELENGTH );

  /* Want this to behave like realloc(3), i.e.
   * *vector == NULL => create a new vector
   * length == 0 => destroy the vector
   * otherwise => resize given vector
   */

  *vector = XFUNC ( *vector, length );
  if ( xlalErrno )
  {
    int code = xlalErrno;
    XLALClearErrno();
    if ( code == XLAL_EBADLEN )
    {
      ABORT( status, AVFACTORIESH_ELENGTH, AVFACTORIESH_MSGELENGTH );
    }
    if ( code == XLAL_ENOMEM )
    {
      ABORT( status, AVFACTORIESH_EMALLOC, AVFACTORIESH_MSGEMALLOC );
    }
  }

  RETURN( status );
}

#undef VTYP
#undef RESIZEVECTOR
#undef XFUNC
#undef XCFUNC
#undef XDFUNC
