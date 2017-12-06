#define CONCAT2x(a,b) a##b
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a##b##c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define STRING(a) #a

#define ATYPE CONCAT2(TYPE,Array)

#ifdef TYPECODE
#define FUNC CONCAT3(LAL,TYPECODE,DestroyArray)
#define XFUNC CONCAT2(XLALDestroy,ATYPE)
#else
#define FUNC LALDestroyArray
#define XFUNC XLALDestroyArray
#endif

void XFUNC ( ATYPE *array )
{
  if ( ! array )
    XLAL_ERROR_VOID( XLAL_EFAULT );
  if ( ! array->dimLength
      || ! array->dimLength->length
      || ! array->dimLength->data
      || ! array->data )
    XLAL_ERROR_VOID( XLAL_EINVAL );
  XLALDestroyUINT4Vector( array->dimLength );
  LALFree( array->data );
  LALFree( array );
  return;
}


void FUNC ( LALStatus *status, ATYPE **array )
{
  INITSTATUS(status);

  ASSERT (array,          status, AVFACTORIESH_EVPTR, AVFACTORIESH_MSGEVPTR);
  ASSERT (*array,         status, AVFACTORIESH_EUPTR, AVFACTORIESH_MSGEUPTR);
  ASSERT ((*array)->data, status, AVFACTORIESH_EDPTR, AVFACTORIESH_MSGEDPTR);

  /* Free allocated storage */

  XFUNC ( *array );
  if ( xlalErrno )
  {
    int code = xlalErrno;
    XLALClearErrno();
    if ( code == XLAL_EFAULT )
    {
      ABORT (status, AVFACTORIESH_EUPTR, AVFACTORIESH_MSGEUPTR);
    }
    if ( code == XLAL_EINVAL )
    {
      ABORT (status, AVFACTORIESH_EDPTR, AVFACTORIESH_MSGEDPTR);
    }
  }
  *array = NULL;	    /* make sure we don't point to freed struct */

  RETURN (status);
}

#undef ATYPE
#undef FUNC
#undef XFUNC
