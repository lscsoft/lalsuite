/*
 * Copyright (C) 2008 Reinhard Prix
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

/**
 * \file
 *
 * \author Reinhard Prix
 * \ingroup factories
 * \brief  Creation/destruction/manipulation functions for 'LALStringVector' type objects.
 *
 */

/*---------- INCLUDES ----------*/
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include <lal/XLALError.h>
#include <lal/StringVector.h>
#include <lal/LALMalloc.h>

/*---------- DEFINES ----------*/
/*---------- internal types ----------*/
static CHAR *deblank_string ( const CHAR *start, UINT4 len );

/*---------- empty initializers ---------- */
const LALStringVector empty_LALStringVector;

/*---------- Global variables ----------*/
/*---------- internal prototypes ----------*/


/*==================== FUNCTION DEFINITIONS ====================*/

/** Append the given string to the string-vector (XLAL interface), return
 * pointer to the resulting string-vector, or NULL on error.
 *
 * \note It is allowed to pass NULL as the input string-vector 'vect', which corresponds
 * to creating a new string-vector with a single element 'string'
 */
LALStringVector *
XLALAppendString2Vector (LALStringVector *vect,		/**< input string-vector to append to */
			 const CHAR *string		/**< string to append */
			 )
{
  const CHAR *fn = __func__;
  UINT4 oldlen;
  LALStringVector *ret;

  if ( !string ) {
    XLALPrintError ("\n%s: NULL 'string' passed to append\n\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  if ( ! vect )
    { /* special case: NULL string-vector to append to */
      if ( (ret = XLALCreateStringVector ( string, NULL )) == NULL) {
        XLALPrintError ("%s: XLALCreateStringVector() failed!\n", fn );
        XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
      }
    }
  else
    {
      ret = vect;
      oldlen = ret->length;

      if ( (ret->data = XLALRealloc ( ret->data, (oldlen + 1)*sizeof( *ret->data ) )) == NULL ) {
        XLALPrintError ("%s: XLALRealloc(%d) failed!\n", fn, oldlen + 1 );
	XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
      }

      ret->length ++;

      if ( (ret->data[oldlen] = XLALCalloc(1, strlen(string) + 1 )) == NULL ) {
        XLALPrintError ("%s: XLALCalloc(%d) failed!\n", fn, strlen(string) + 1 );
	XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
      }

      strcpy ( ret->data[oldlen], string );
    }

  return ret;

} /* XLALAppendString2Vector() */




/** Create a StringVector from the list of strings passed as arguments
 * \note All arguments MUST be CHAR* strings.
 * The last argument MUST be NULL, as C cannot deduce the number of arguments
 * otherwise.
 */
LALStringVector *
XLALCreateStringVector ( const CHAR *str1, ... )
{
  const CHAR *fn = __func__;
  LALStringVector *ret;
  const CHAR *next;
  va_list ap;

  if ( !str1 ) {
    XLALPrintError ("%s: invalid NULL input string 'str1'\n", fn );
    XLAL_ERROR_NULL (fn, XLAL_EINVAL );
  }

  size_t len;

  /* set up return vector of strings, and handle first argument */
  len = sizeof(*ret);
  if ( (ret = XLALCalloc ( 1, len )) == NULL )
    goto failed;

  len = sizeof(ret->data[0]);
  if ( (ret->data = XLALCalloc ( 1, len )) == NULL)
    goto failed;

  len = strlen(str1)+1;
  if ( (ret->data[0] = XLALCalloc ( len, sizeof(CHAR) )) == NULL )
    goto failed;

  strcpy ( ret->data[0], str1 );
  ret->length ++;

  /* handle remaining variable-length list of (assumed)string arguments */
  va_start(ap, str1);

  while ( (next = va_arg(ap, const CHAR *)) != NULL )
    {
      ret->length ++;

      len = ret->length * sizeof(ret->data[0]);
      if ( (ret->data = XLALRealloc ( ret->data, len)) == NULL )
	goto failed;

      len = strlen(next)+1;
      if ( (ret->data[ret->length-1] = XLALCalloc( len, sizeof(CHAR) )) == NULL )
	goto failed;

      strcpy ( ret->data[ret->length-1], next );

    } /* while more arguments */

  va_end(ap);

  return ret;

 failed:
  va_end(ap);
  XLALPrintError ("%s: failed to allocate '%d' bytes\n", fn, len );
  XLALDestroyStringVector ( ret );
  XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );

} /* XLALCreateStringVector() */


/** XLAL-interface: Free a string-vector ;) */
void
XLALDestroyStringVector ( LALStringVector *vect )
{
  UINT4 i;

  if ( !vect )
    return;

  if ( vect->data )
    {
      for ( i=0; i < vect->length; i++ )
	{
	  if ( vect->data[i] )
	    XLALFree ( vect->data[i] );
	}

      XLALFree ( vect->data );
    }

  XLALFree ( vect );

  return;

} /* XLALDestroyStringVector() */



/* comparison function for strings */
static int StringCompare (const void *p1, const void *p2)
{
  const char *s1 = p1;
  const char *s2 = p2;
  return (strcmp ( s1, s2 ) );
}

/** Sort string-vector alphabetically
 */
int
XLALSortStringVector (LALStringVector *strings)
{
  const CHAR *fn = __func__;

  if ( !strings || strings->length == 0 ) {
    XLALPrintError ("%s: invalid empty or zero-length input 'strings'\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  qsort ( (void*)(strings->data), (size_t)(strings->length), sizeof(CHAR*), StringCompare );

  return XLAL_SUCCESS;

} /* XLALSortStringVector() */



/** Parse a list of comma-separated values (CSV) into a StringVector
 * \note surrounding whitespace is removed from the 'values'.
 */
LALStringVector *
XLALParseCSV2StringVector ( const CHAR *CSVlist )
{
  const CHAR *fn = __func__;
  UINT4 counter;
  const CHAR *start, *tmp;
  CHAR **data = NULL;
  LALStringVector *ret = NULL;

  if ( !CSVlist )
    return NULL;

  size_t len;
  /* prepare return string-vector */
  len = sizeof( *ret );
  if ( ( ret = XLALCalloc ( 1, len ) ) == NULL )
    goto failed;

  start = CSVlist;
  counter = 0;
  do
    {

      /* extend string-array */
      len = (counter+1) * sizeof(CHAR*);
      if ( ( data = XLALRealloc ( data, len )) == NULL )
	goto failed;

      /* determine string-length of next value */
      if ( ( tmp = strchr ( start, ',' ) ) )
	len = tmp - start;
      else
	len = strlen ( start );

      /* allocate space for that value in string-array */
      if ( (data[counter] = deblank_string ( start, len ) ) == NULL ) {
        XLALDestroyStringVector ( ret );
        XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
      }
      counter ++;

    } while ( tmp && (start = tmp + 1) );

  ret -> length = counter;
  ret -> data = data;

  /* success: */
  return ( ret );

 failed:
  XLALPrintError ("%s: failed to allocate %d bytes\n", fn, len );
  XLALDestroyStringVector ( ret );
  XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );

} /* XLALParseCSV2StringVector() */


/** Copy (and allocate) string from 'start' with length 'len', removing
 * all starting- and trailing blanks!
 */
CHAR *
deblank_string ( const CHAR *start, UINT4 len )
{
  const char *fn = __func__;

  const CHAR *blank_chars = " \t\n";
  const CHAR *pos0, *pos1;
  UINT4 newlen;
  CHAR *ret;

  if ( !start || !len ) {
    XLALPrintError ("%s: invalid NULL input 'start' or len=0\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  /* clip from beginning */
  pos0 = start;
  pos1 = start + len - 1;
  while ( (pos0 < pos1) && strchr ( blank_chars, *pos0 ) )
    pos0 ++;

  /* clip backwards from end */
  while ( (pos1 >= pos0) && strchr ( blank_chars, *pos1 ) )
    pos1 --;

  newlen = pos1 - pos0 + 1;
  if ( !newlen ) {
    XLALPrintError ("%s: something went wrong here .. probably a coding mistake.\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EFAILED );
  }

  if ( (ret = XLALCalloc(1, newlen + 1)) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc(1, %d )\n", fn, newlen + 1 );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  strncpy ( ret, pos0, newlen );
  ret[ newlen ] = 0;

  return ret;

} /* deblank_string() */

/** Search for string 'needle' in string-vector 'haystack', return index to
 * first matching vector element if found, -1 outherwise.
 *
 * Note: function allows haystack=NULL input, in which case -1 (=not found) will be returned.
 *
 */
INT4
XLALFindStringInVector ( const char *needle, const LALStringVector *haystack )
{
  const char *fn = __func__;

  if ( !needle ) {
    XLALPrintError ("%s: invalid NULL input 'needle'!\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  if ( !haystack || (haystack->length == 0) )	// no vector to search => not found
    return -1;

  UINT4 i;
  for ( i=0; i < haystack->length; i ++ )
    if ( !strcmp ( needle, haystack->data[i] ) )	// found it!
      return i;

  return -1;	// didn't find matching entry

} /* XLALFindStringInVector() */
