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
define(`ATYPE',`format(`%sArray',TYPE)')
define(`FUNC',`format(`LAL%sCreateArray',TYPECODE)')
ifelse( TYPECODE, `', `define(`XFUNC',`XLALCreateArray')', `define(`XFUNC',`format(`XLALCreate%s',ATYPE)')' )
ifelse( TYPECODE, `', `define(`XFUNCL',`XLALCreateArrayL')', `define(`XFUNCL',`format(`XLALCreate%sL',ATYPE)')' )
ifelse( TYPECODE, `', `define(`XFUNCV',`XLALCreateArrayV')', `define(`XFUNCV',`format(`XLALCreate%sV',ATYPE)')' )

ATYPE * XFUNCL ( UINT4 ndim, ... )
{
  enum { maxdim = 16 };
  va_list ap;
  ATYPE *arr;
  UINT4 dims[maxdim];
  UINT4 dim;

  if ( ! ndim )
    XLAL_ERROR_NULL( "XFUNCL", XLAL_EBADLEN );
  if ( ndim > maxdim )
    XLAL_ERROR_NULL( "XFUNCL", XLAL_EINVAL );

  va_start( ap, ndim );
  for ( dim = 0; dim < ndim; ++dim )
    dims[dim] = va_arg( ap, UINT4 );
  va_end( ap );

  arr = XFUNCV ( ndim, dims );
  if ( ! arr )
    XLAL_ERROR_NULL( "XFUNCL", XLAL_EFUNC );
  return arr; 
}

ATYPE * XFUNCV ( UINT4 ndim, UINT4 *dims )
{
  ATYPE *arr;
  UINT4Vector dimLength;

  if ( ! ndim )
    XLAL_ERROR_NULL( "XFUNCV", XLAL_EBADLEN );
  if ( ! dims )
    XLAL_ERROR_NULL( "XFUNCV", XLAL_EFAULT );

  dimLength.length = ndim;
  dimLength.data   = dims;

  arr = XFUNC ( &dimLength );
  if ( ! arr )
    XLAL_ERROR_NULL( "XFUNCV", XLAL_EFUNC );
  return arr;
}


ATYPE * XFUNC ( UINT4Vector *dimLength )
{
  ATYPE *arr;
  UINT4 size = 1;
  UINT4 ndim;
  UINT4 dim;

  if ( ! dimLength )
    XLAL_ERROR_NULL( "XFUNC", XLAL_EFAULT );
  if ( ! dimLength->length )
    XLAL_ERROR_NULL( "XFUNC", XLAL_EBADLEN );
  if ( ! dimLength->data )
    XLAL_ERROR_NULL( "XFUNC", XLAL_EINVAL );

  ndim = dimLength->length;
  for ( dim = 0; dim < ndim; ++dim )
    size *= dimLength->data[dim];

  if ( ! size )
    XLAL_ERROR_NULL( "XFUNC", XLAL_EBADLEN );

  /* create array */
  arr = LALMalloc( sizeof( *arr ) );
  if ( ! arr )
    XLAL_ERROR_NULL( "XFUNC", XLAL_ENOMEM );

  /* create array dimensions */
  arr->dimLength = XLALCreateUINT4Vector( ndim );
  if ( ! arr->dimLength )
  {
    LALFree( arr );
    XLAL_ERROR_NULL( "XFUNC", XLAL_EFUNC );
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
    XLAL_ERROR_NULL( "XFUNC", XLAL_ENOMEM );
  }

  return arr;
}


/* <lalVerbatim file="ArrayFactoriesD"> */
void FUNC ( LALStatus *status, ATYPE **array, UINT4Vector *dimLength ) 
{ /* </lalVerbatim> */
  INITSTATUS (status, "FUNC", ARRAYFACTORIESC);   

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

