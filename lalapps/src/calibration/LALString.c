#include <string.h>
#include <lal/LALStdlib.h>

/* like strcat but dynamically reallocates string */
char * XLALStringAppend( char *s, const char *append )
{
  static const char *func = "XLALStringAppend";
  size_t curlen;
  size_t newlen;
  if ( ! append )
    return s;
  curlen = s ? strlen( s ) : 0;
  newlen = curlen + strlen( append );
  s = LALRealloc( s, newlen + 1 );
  if ( ! s )
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );
  return strcpy( s + curlen, append );
}

char * XLALStringDuplicate( const char *s )
{
  char *dup;
  dup = XLALStringAppend( NULL, s );
  return dup;
}

/* should be nearly equivalent to strlcpy */
size_t XLALStringCopy( char *dst, const char *src, size_t size )
{
  size_t srclen;
  if ( ! src )
    src = "";
  srclen = strlen( src );
  if ( ! dst || size < 1 ) /* no copy */
    return srclen;
  if ( size == 1 ) /* NUL terminate and exit */
  {
    dst[0] = 0;
    return srclen;
  }
  strncpy( dst, src, size - 1 );
  dst[size-1] = 0;
  return srclen;
}

/* should be nearly equivalent to strlcat */
size_t XLALStringConcatenate( char *dst, const char *src, size_t size )
{
  size_t srclen;
  size_t dstlen;
  if ( ! src )
    src = "";
  srclen = strlen( src );
  if ( ! dst || size < 1 ) /* no copy */
    return srclen;
  if ( size == 1 ) /* NUL terminate and exit */
  {
    dst[0] = 0;
    return srclen;
  }
  dstlen = strlen( dst );
  strncat( dst, src, size - dstlen - 1 );
  return srclen + dstlen;
}
