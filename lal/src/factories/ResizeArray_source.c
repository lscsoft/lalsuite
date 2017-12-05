#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define ATYPE CONCAT2(TYPE,Array)

#ifdef TYPECODE
#define FUNC CONCAT3(LAL,TYPECODE,ResizeArray)
#define XFUNC CONCAT2(XLALResize,ATYPE)
#define XFUNCL CONCAT3(XLALResize,ATYPE,L)
#define XFUNCV CONCAT3(XLALResize,ATYPE,V)
#define XCFUNC CONCAT2(XLALCreate,ATYPE)
#define XDFUNC CONCAT2(XLALDestroy,ATYPE)
#else
#define FUNC LALResizeArray
#define XFUNC XLALResizeArray
#define XFUNCL XLALResizeArrayL
#define XFUNCV XLALResizeArrayV
#define XCFUNC XLALCreateArray
#define XDFUNC XLALDestroyArray
#endif

ATYPE * XFUNCL ( ATYPE *array, UINT4 ndim, ... )
{
  enum { maxdim = 16 };
  va_list ap;
  UINT4 dims[maxdim];
  UINT4 dim;

  if ( ! ndim )
  {
    XDFUNC ( array );
    return NULL;
  }
  if ( ndim > maxdim )
    XLAL_ERROR_NULL( XLAL_EINVAL );

  va_start( ap, ndim );
  for ( dim = 0; dim < ndim; ++dim )
    dims[dim] = va_arg( ap, UINT4 );
  va_end( ap );

  return XFUNCV ( array, ndim, dims );
}

ATYPE * XFUNCV ( ATYPE *array, UINT4 ndim, UINT4 *dims )
{
  UINT4Vector dimLength;

  if ( ! ndim )
  {
    XDFUNC ( array );
    return NULL;
  }
  if ( ! dims )
    XLAL_ERROR_NULL( XLAL_EINVAL );

  dimLength.length = ndim;
  dimLength.data   = dims;

  return XFUNC ( array, &dimLength );
}

ATYPE * XFUNC ( ATYPE *array, UINT4Vector *dimLength )
{
  UINT4 size = 1;
  UINT4 ndim;
  UINT4 dim;

  if ( ! array )
    return XCFUNC ( dimLength );
  if ( ! dimLength )
  {
    XDFUNC ( array );
    return NULL;
  }
  if ( ! dimLength->length )
    XLAL_ERROR_NULL( XLAL_EBADLEN );
  if ( ! dimLength->data )
    XLAL_ERROR_NULL( XLAL_EINVAL );

  ndim = dimLength->length;
  for ( dim = 0; dim < ndim; ++dim )
    size *= dimLength->data[dim];

  if ( ! size )
    XLAL_ERROR_NULL( XLAL_EBADLEN );

  /* resize array->dimLength vector if needed */
  if ( array->dimLength->length != ndim )
  {
    array->dimLength = XLALResizeUINT4Vector( array->dimLength, ndim );
    if ( ! array->dimLength )
      XLAL_ERROR_NULL( XLAL_EFUNC );
  }

  /* copy dimension lengths */
  memcpy( array->dimLength->data, dimLength->data,
      ndim * sizeof( *array->dimLength->data ) );

  /* reallocate data storage */
  array->data = LALRealloc( array->data, size * sizeof( *array->data ) );
  if ( ! array->data )
    XLAL_ERROR_NULL( XLAL_ENOMEM );

  return array;
}


void FUNC ( LALStatus *status, ATYPE **array, UINT4Vector *dimLength )
{
  ATYPE *tmparr = NULL;

  INITSTATUS(status);

  ASSERT ( array != NULL, status, AVFACTORIESH_EVPTR, AVFACTORIESH_MSGEVPTR );

  if ( ! array )
  {
    tmparr = XCFUNC ( dimLength );
  }
  else if ( ! dimLength )
  {
    XDFUNC ( *array );
    *array = NULL;
  }
  else
  {
    tmparr = XFUNC ( *array, dimLength );
  }

  if ( xlalErrno )
  {
    int code = xlalErrno;
    XLALClearErrno();
    if ( code == XLAL_EINVAL )
    {
      ABORT (status, AVFACTORIESH_EVPTR, AVFACTORIESH_MSGEVPTR);
    }
    if ( code == XLAL_EBADLEN )
    {
      ABORT (status, AVFACTORIESH_ELENGTH, AVFACTORIESH_MSGELENGTH);
    }
    if ( code == XLAL_ENOMEM )
    {
      ABORT (status, AVFACTORIESH_EMALLOC, AVFACTORIESH_MSGEMALLOC);
    }
  }

  *array = tmparr;

  RETURN (status);
}

#undef ATYPE
#undef FUNC
#undef XFUNC
#undef XFUNCL
#undef XFUNCV
#undef XCFUNC
#undef XDFUNC
