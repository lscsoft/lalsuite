/** \file
 * \ingroup std
 * \author Creighton, J. D. E.
 * \date $Date$
 * \brief XLAL string manipulation routines.
 *
 * $Id$
 */

#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALString.h>

NRCSID( LALSTRINGC, "$Id$" );

/** Like strcat but dynamically reallocates string with LALRealloc. */
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
  strcpy( s + curlen, append );
  return s;
}

/** Like strdup but uses LAL allocation routines (free with LALFree). */
char * XLALStringDuplicate( const char *s )
{
  char *dup;
  dup = XLALStringAppend( NULL, s );
  return dup;
}

/** Copy sources string src to destination string dst.
 * Up to size - 1 characters are copied and destination string dst is
 * guaranteed to be NUL-terminated.
 * Return value is the length of source string src.  If this is greater than
 * or equal to the size of the destination string buffer, size, then truncation
 * has occurred. Should be nearly equivalent to strlcpy. */
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

/** Concatenate sources string src to the end of destination string dst.
 * Characters are added to destination string dst until the source string src
 * is exhausted or the length of destination string dst is size - 1 characters.
 * Destination string dst is guaranteed to be NUL-terminated.
 * Return value is the initial length of destination string dst plus the
 * length of source string src.  If this is greater than
 * or equal to the size of the destination string buffer, size, then truncation
 * has occurred. Should be nearly equivalent to strlcat. */
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
