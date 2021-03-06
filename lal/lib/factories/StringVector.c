/*
 * Copyright (C) 2015 Karl Wette
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

/*---------- INCLUDES ----------*/
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include <lal/XLALError.h>
#include <lal/StringVector.h>
#include <lal/LALMalloc.h>
#include <lal/LALString.h>

/*---------- DEFINES ----------*/
/*---------- internal types ----------*/

/*---------- Global variables ----------*/
/*---------- internal prototypes ----------*/

/*==================== FUNCTION DEFINITIONS ====================*/

/**
 * Append the given string to the string-vector (XLAL interface), return
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
  UINT4 oldlen;
  LALStringVector *ret;

  if ( !string ) {
    XLALPrintError ("\n%s: NULL 'string' passed to append\n\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  if ( ! vect )
    { /* special case: NULL string-vector to append to */
      if ( (ret = XLALCreateStringVector ( string, NULL )) == NULL) {
        XLALPrintError ("%s: XLALCreateStringVector() failed!\n", __func__ );
        XLAL_ERROR_NULL ( XLAL_EFUNC );
      }
    }
  else
    {
      ret = vect;
      oldlen = ret->length;

      if ( (ret->data = XLALRealloc ( ret->data, (oldlen + 1)*sizeof( *ret->data ) )) == NULL ) {
        XLALPrintError ("%s: XLALRealloc(%d) failed!\n", __func__, oldlen + 1 );
	XLAL_ERROR_NULL ( XLAL_ENOMEM );
      }

      ret->length ++;

      if ( (ret->data[oldlen] = XLALCalloc(1, strlen(string) + 1 )) == NULL ) {
        XLALPrintError ("%s: XLALCalloc(%zu) failed!\n", __func__, strlen(string) + 1 );
	XLAL_ERROR_NULL ( XLAL_ENOMEM );
      }

      strcpy ( ret->data[oldlen], string );
    }

  return ret;

} /* XLALAppendString2Vector() */




/**
 * Create a StringVector from the list of strings passed as arguments
 * \note All arguments MUST be CHAR* strings.
 * The last argument MUST be NULL, as C cannot deduce the number of arguments
 * otherwise.
 */
LALStringVector *
XLALCreateStringVector ( const CHAR *str1, ... )
{
  LALStringVector *ret;
  const CHAR *next;
  va_list ap;

  if ( !str1 ) {
    XLALPrintError ("%s: invalid NULL input string 'str1'\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
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
  XLALPrintError ("%s: failed to allocate '%zu' bytes\n", __func__, len );
  XLALDestroyStringVector ( ret );
  XLAL_ERROR_NULL ( XLAL_ENOMEM );

} /* XLALCreateStringVector() */


/**
 * Create an empty string vector of the given length
 */
LALStringVector *XLALCreateEmptyStringVector ( UINT4 length )
{
  LALStringVector * vector;
  vector = LALMalloc( sizeof( *vector ) );
  if ( ! vector )
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  vector->length = length;
  if ( ! length ) /* zero length: set data pointer to be NULL */
    vector->data = NULL;
  else /* non-zero length: allocate memory for data */
  {
    vector->data = LALCalloc( length, sizeof( *vector->data ) );
    if ( ! vector->data )
    {
      LALFree( vector );
      XLAL_ERROR_NULL( XLAL_ENOMEM );
    }
  }
  return vector;
}


/**
 * Create a copy of a string vector
 */
LALStringVector *XLALCopyStringVector( const LALStringVector *vect )
{
  XLAL_CHECK_NULL( vect != NULL && vect->data != NULL && vect->length > 0, XLAL_EINVAL );
  LALStringVector *copy = XLALCreateStringVector( vect->data[0], NULL );
  XLAL_CHECK_NULL( copy != NULL, XLAL_EFUNC );
  for (size_t i = 1; i < vect->length; ++i) {
    copy = XLALAppendString2Vector( copy, vect->data[i] );
    XLAL_CHECK_NULL( copy != NULL, XLAL_EFUNC );
  }
  return copy;
} /* XLALCopyStringVector() */


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
  /* this formulation explicitly follows the example given in 'man qsort' for string-array sorting
   * Quoting from there:
   ** The actual arguments to this function are "pointers to
   ** pointers to char", but strcmp(3) arguments are "pointers
   ** to char", hence the following cast plus dereference
   *
   */
  return strcmp ( * ( char * const *) p1, * ( char * const *) p2 );
}

///
/// Concatenate a string vector 'strings', separated by the string 'sep', into a single string
///
char *XLALConcatStringVector( const LALStringVector *strings, const char *sep )
{

  // Check input
  XLAL_CHECK_NULL( strings != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( sep != NULL, XLAL_EFAULT );

  // Allocate a string of the required size
  int len = 1;
  if ( strings->length > 0 ) {
    len += strlen( strings->data[0] );
    for ( size_t i = 1; i < strings->length; ++i ) {
      len += strlen( sep ) + strlen( strings->data[i] );
    }
  }
  char *s = XLALCalloc( len, sizeof(*s) );
  XLAL_CHECK_NULL( s != NULL, XLAL_ENOMEM );

  // Concatenate the vector of strings into a single string
  if ( strings->length > 0 ) {
    strcpy( s, strings->data[0] );
    for ( size_t i = 1; i < strings->length; ++i ) {
      strcat( s, sep );
      strcat( s, strings->data[i] );
    }
  }

  return s;

}

///
/// Parse 'string' into a string vector of tokens, delimited by the characters 'delim'
///
LALStringVector *XLALParseStringVector( const char *string, const char *delim )
{

  // Check input
  XLAL_CHECK_NULL( string != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( delim != NULL, XLAL_EFAULT );

  // Allocate an empty string vector
  LALStringVector *strings = XLALCalloc( 1, sizeof(*strings) );
  XLAL_CHECK_NULL( strings != NULL, XLAL_ENOMEM );

  // Split the string into a vector of strings
  char *s = XLALStringDuplicate( string );
  XLAL_CHECK_NULL( s != NULL, XLAL_EFUNC );
  char *p = s;
  char *t = XLALStringToken( &p, delim, 1 );
  while ( t != NULL ) {
    XLAL_CHECK_NULL( XLALAppendString2Vector( strings, t ) != NULL, XLAL_EFUNC );
    t = XLALStringToken( &p, delim, 1 );
  }

  // Cleanup
  XLALFree( s );

  return strings;

}

/**
 * Sort string-vector alphabetically *in place*
 */
int
XLALSortStringVector (LALStringVector *strings)
{
  if ( !strings || strings->length == 0 ) {
    XLALPrintError ("%s: invalid empty or zero-length input 'strings'\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  qsort ( (void*)(&strings->data[0]), (size_t)(strings->length), sizeof(strings->data[0]), StringCompare );

  return XLAL_SUCCESS;

} /* XLALSortStringVector() */

/**
 * Copy (and allocate) string from 'start' with length 'len', removing
 * all starting- and trailing blanks!
 */
char *
XLALDeblankString ( const CHAR *start, UINT4 len )
{
  XLAL_CHECK_NULL ( start != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( len > 0, XLAL_EDOM );

  const CHAR *blank_chars = " \t\n";

  /* clip from beginning */
  const char *pos0 = start;
  const char *pos1 = start + len - 1;
  while ( (pos0 < pos1) && strchr ( blank_chars, (*pos0) ) ) {
    pos0 ++;
  }

  /* clip backwards from end */
  while ( (pos1 >= pos0) && strchr ( blank_chars, (*pos1) ) ) {
    pos1 --;
  }

  UINT4 newlen = pos1 - pos0 + 1;
  XLAL_CHECK_NULL ( newlen > 0, XLAL_EFAILED, "newlen==0: Something went wrong here .. probably a coding mistake.\n" );

  CHAR *ret;
  XLAL_CHECK_NULL ( (ret = XLALCalloc(1, newlen + 1)) != NULL, XLAL_ENOMEM );

  strncpy ( ret, pos0, newlen );
  ret[ newlen ] = 0;

  return ret;

} /* XLALDeblankString() */

/**
 * Search for string 'needle' in string-vector 'haystack', return index to
 * first matching vector element if found, -1 outherwise.
 *
 * Note: function allows haystack=NULL input, in which case -1 (=not found) will be returned.
 *
 */
INT4
XLALFindStringInVector ( const char *needle, const LALStringVector *haystack )
{
  if ( !needle ) {
    XLALPrintError ("%s: invalid NULL input 'needle'!\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  if ( !haystack || (haystack->length == 0) )	// no vector to search => not found
    return -1;

  UINT4 i;
  for ( i=0; i < haystack->length; i ++ )
    if ( !strcmp ( needle, haystack->data[i] ) )	// found it!
      return i;

  return -1;	// didn't find matching entry

} /* XLALFindStringInVector() */
