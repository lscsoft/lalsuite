#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define ATYPE CONCAT2(TYPE,Array)

#ifdef TYPECODE
#define FUNC CONCAT3(LAL,TYPECODE,CreateArray)
#define XFUNC CONCAT2(XLALCreate,ATYPE)
#define XFUNCL CONCAT3(XLALCreate,ATYPE,L)
#define XFUNCV CONCAT3(XLALCreate,ATYPE,V)
#else
#define FUNC LALCreateArray
#define XFUNC XLALCreateArray
#define XFUNCL XLALCreateArrayL
#define XFUNCV XLALCreateArrayV
#endif

ATYPE * XFUNCL ( UINT4 ndim, ... )
{
  enum { maxdim = 16 };
  va_list ap;
  ATYPE *arr;
  UINT4 dims[maxdim];
  UINT4 dim;

  if ( ! ndim )
    XLAL_ERROR_NULL( XLAL_EBADLEN );
  if ( ndim > maxdim )
    XLAL_ERROR_NULL( XLAL_EINVAL );

  va_start( ap, ndim );
  for ( dim = 0; dim < ndim; ++dim )
    dims[dim] = va_arg( ap, UINT4 );
  va_end( ap );

  arr = XFUNCV ( ndim, dims );
  if ( ! arr )
    XLAL_ERROR_NULL( XLAL_EFUNC );
  return arr;
}

ATYPE * XFUNCV ( UINT4 ndim, UINT4 *dims )
{
  ATYPE *arr;
  UINT4Vector dimLength;

  if ( ! ndim )
    XLAL_ERROR_NULL( XLAL_EBADLEN );
  if ( ! dims )
    XLAL_ERROR_NULL( XLAL_EFAULT );

  dimLength.length = ndim;
  dimLength.data   = dims;

  arr = XFUNC ( &dimLength );
  if ( ! arr )
    XLAL_ERROR_NULL( XLAL_EFUNC );
  return arr;
}


ATYPE * XFUNC ( UINT4Vector *dimLength )
{
  ATYPE *arr;
  UINT4 size = 1;
  UINT4 ndim;
  UINT4 dim;

  if ( ! dimLength )
    XLAL_ERROR_NULL( XLAL_EFAULT );
  if ( ! dimLength->length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );
  if ( ! dimLength->data )
    XLAL_ERROR_NULL( XLAL_EINVAL );

  ndim = dimLength->length;
  for ( dim = 0; dim < ndim; ++dim )
    size *= dimLength->data[dim];

  if ( ! size )
    XLAL_ERROR_NULL( XLAL_EBADLEN );

  /* create array */
  arr = LALMalloc( sizeof( *arr ) );
  if ( ! arr )
    XLAL_ERROR_NULL( XLAL_ENOMEM );

  /* create array dimensions */
  arr->dimLength = XLALCreateUINT4Vector( ndim );
  if ( ! arr->dimLength )
  {
    LALFree( arr );
    XLAL_ERROR_NULL( XLAL_EFUNC );
  }

  /* copy dimension lengths */
  memcpy( arr->dimLength->data, dimLength->data,
      ndim * sizeof( *arr->dimLength->data ) );

  /* allocate data storage */
  arr->data = LALMalloc( size * sizeof( *arr->data ) );
  if ( ! arr->data )
  {
    XLALDestroyUINT4Vector( arr->dimLength );
    LALFree( arr );
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }

  return arr;
}



void FUNC ( LALStatus *status, ATYPE **array, UINT4Vector *dimLength )
{
  INITSTATUS(status);

  /* make sure arguments are sane */

  ASSERT (array,             status, AVFACTORIESH_EVPTR, AVFACTORIESH_MSGEVPTR);
  ASSERT (!*array,           status, AVFACTORIESH_EUPTR, AVFACTORIESH_MSGEUPTR);
  ASSERT (dimLength,         status, AVFACTORIESH_EVPTR, AVFACTORIESH_MSGEVPTR);
  ASSERT (dimLength->data,   status, AVFACTORIESH_EVPTR, AVFACTORIESH_MSGEVPTR);
  ASSERT (dimLength->length, status,
          AVFACTORIESH_ELENGTH, AVFACTORIESH_MSGELENGTH);

  *array = XFUNC ( dimLength );
  if ( ! *array )
  {
    int code = xlalErrno & ~XLAL_EFUNC; /* turn off subfunction error bit */
    XLALClearErrno();
    if ( code & XLAL_EFAULT )
    {
      ABORT (status, AVFACTORIESH_EVPTR, AVFACTORIESH_MSGEVPTR);
    }
    if ( code == XLAL_EBADLEN )
    {
      ABORT (status, AVFACTORIESH_ELENGTH, AVFACTORIESH_MSGELENGTH);
    }
    if ( code == XLAL_EINVAL )
    {
      ABORT (status, AVFACTORIESH_EVPTR, AVFACTORIESH_MSGEVPTR);
    }
    if ( code == XLAL_ENOMEM )
    {
      ABORT (status, AVFACTORIESH_EMALLOC, AVFACTORIESH_MSGEMALLOC);
    }
  }

  RETURN (status);
}

#undef ATYPE
#undef FUNC
#undef XFUNC
#undef XFUNCL
#undef XFUNCV
