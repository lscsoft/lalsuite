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
ifelse(TYPECODE,`CHAR',`define(`TYPE',`CHAR')')
ifelse(TYPECODE,`',`define(`TYPE',`REAL4')')
define(`ATYPE',`format(`%sArray',TYPE)')
define(`FUNC',`format(`LAL%sResizeArray',TYPECODE)')
ifelse( TYPECODE, `', `define(`XFUNC',`XLALResizeArray')', `define(`XFUNC',`format(`XLALResize%s',ATYPE)')' )
ifelse( TYPECODE, `', `define(`XFUNCL',`XLALResizeArrayL')', `define(`XFUNCL',`format(`XLALResize%sL',ATYPE)')' )
ifelse( TYPECODE, `', `define(`XFUNCV',`XLALResizeArrayV')', `define(`XFUNCV',`format(`XLALResize%sV',ATYPE)')' )
ifelse( TYPECODE, `', `define(`XCFUNC',`XLALCreateArray')', `define(`XCFUNC',`format(`XLALCreate%s',ATYPE)')' )
ifelse( TYPECODE, `', `define(`XDFUNC',`XLALDestroyArray')', `define(`XDFUNC',`format(`XLALDestroy%s',ATYPE)')' )


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
    XLAL_ERROR_NULL( "XFUNCL", XLAL_EINVAL );

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
    XLAL_ERROR_NULL( "XFUNCV", XLAL_EINVAL );

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
    XLAL_ERROR_NULL( "XFUNC", XLAL_EBADLEN );
  if ( ! dimLength->data )
    XLAL_ERROR_NULL( "XFUNC", XLAL_EINVAL );

  ndim = dimLength->length;
  for ( dim = 0; dim < ndim; ++dim )
    size *= dimLength->data[dim];

  if ( ! size )
    XLAL_ERROR_NULL( "XFUNC", XLAL_EBADLEN );
  
  /* resize array->dimLength vector if needed */
  if ( array->dimLength->length != ndim )
  {
    array->dimLength = XLALResizeUINT4Vector( array->dimLength, ndim );
    if ( ! array->dimLength )
      XLAL_ERROR_NULL( "XFUNC", XLAL_EFUNC );
  }

  /* copy dimension lengths */
  memcpy( array->dimLength->data, dimLength->data,
      ndim * sizeof( *array->dimLength->data ) );

  /* reallocate data storage */
  array->data = LALRealloc( array->data, size * sizeof( *array->data ) );
  if ( ! array->data )
    XLAL_ERROR_NULL( "XFUNC", XLAL_ENOMEM );

  return array;
}

/* <lalVerbatim file="ArrayFactoriesD"> */
void FUNC ( LALStatus *status, ATYPE **array, UINT4Vector *dimLength )
{  /* </lalVerbatim> */
  ATYPE *tmparr = NULL;
  UINT4 arrayDataSize = 1;
  UINT4 numDims;
  UINT4 dim;
  TYPE * p; /* temporary pointer */

  INITSTATUS (status, "FUNC", ARRAYFACTORIESC);

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
