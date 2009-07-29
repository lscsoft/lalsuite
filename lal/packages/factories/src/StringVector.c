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
 * Revision: $Id$
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
NRCSID( STRINGVECTORC, "$Id$");

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
  const CHAR *fn = "XLALAppendString2Vector()";
  UINT4 oldlen;
  LALStringVector *ret;

  if ( !string ) {
    XLALPrintError ("\n%s: NULL string passed to append\n\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  if ( ! vect )
    ret = XLALCreateStringVector ( string, NULL );	/* special case: NULL string-vector to append to */
  else
    {
      ret = vect;
      oldlen = ret->length;

      if ( (ret->data = LALRealloc ( ret->data, (oldlen + 1)*sizeof( *ret->data ) )) == NULL ) {
	XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
      }

      ret->length ++;

      if ( (ret->data[oldlen] = LALCalloc(1, strlen(string) + 1 )) == NULL ) {
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
  const CHAR *fn = "XLALCreateStringVector()";
  LALStringVector *ret;
  const CHAR *next;
  va_list ap;

  if ( !str1 ) {
    XLAL_ERROR_NULL (fn, XLAL_EINVAL );
  }

  /* set up return vector of strings, and handle first argument */
  if ( (ret = LALCalloc ( 1, sizeof(*ret) )) == NULL ) {
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }
  if ( (ret->data = LALCalloc ( 1, sizeof(ret->data[0]) )) == NULL) {
    LALFree ( ret );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }
  if ( (ret->data[0] = LALCalloc ( strlen(str1)+1, sizeof(CHAR) )) == NULL ) {
    LALFree ( ret->data );
    LALFree ( ret );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }
  strcpy ( ret->data[0], str1 );
  ret->length ++;

  /* handle remaining variable-length list of (assumed)string arguments */
  va_start(ap, str1);

  while ( (next = va_arg(ap, const CHAR *)) != NULL )
    {
      ret->length ++;
      if ( (ret->data = LALRealloc ( ret->data, ret->length * sizeof(ret->data[0]))) == NULL )
	goto failed;
      if ( (ret->data[ret->length-1] = LALCalloc( strlen(next)+1, sizeof(CHAR) )) == NULL )
	goto failed;

      strcpy ( ret->data[ret->length-1], next );

    } /* while more arguments */

  va_end(ap);

  return ret;

 failed:
  va_end(ap);
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
	    LALFree ( vect->data[i] );
	}

      LALFree ( vect->data );
    }

  LALFree ( vect );

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
  const CHAR *fn = "XLALSortStringVector()";
  if ( !strings || strings->length == 0 ) {
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
  const CHAR *fn = "XLALParseCSV2StringVector";
  UINT4 counter;
  const CHAR *start, *tmp;
  CHAR **data = NULL;
  LALStringVector *ret = NULL;

  if ( !CSVlist )
    return NULL;

  /* prepare return string-vector */
  if ( ( ret = LALCalloc ( 1, sizeof( *ret )) ) == NULL )
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );

  start = CSVlist;
  counter = 0;
  do
    {
      UINT4 len;

      /* extend string-array */
      if ( ( data = LALRealloc ( data, (counter+1) * sizeof(CHAR*) )) == NULL )
	goto failed;

      /* determine string-length of next value */
      if ( ( tmp = strchr ( start, ',' ) ) )
	len = tmp - start;
      else
	len = strlen ( start );

      /* allocate space for that value in string-array */
      if ( (data[counter] = deblank_string ( start, len ) ) == NULL )
	goto failed;

      counter ++;

    } while ( tmp && (start = tmp + 1) );

  ret -> length = counter;
  ret -> data = data;

  /* success: */
  return ( ret );

 failed:
  XLALDestroyStringVector ( ret );
  XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );

} /* XLALParseCSV2StringVector() */


/** Copy (and allocate) string from 'start' with length 'len', removing
 * all starting- and trailing blanks!
 */
CHAR *
deblank_string ( const CHAR *start, UINT4 len )
{
  const CHAR *blank_chars = " \t\n";
  const CHAR *pos0, *pos1;
  UINT4 newlen;
  CHAR *ret;

  if ( !start || !len )
    return NULL;

  /* clip from beginning */
  pos0 = start;
  pos1 = start + len - 1;
  while ( (pos0 < pos1) && strchr ( blank_chars, *pos0 ) )
    pos0 ++;

  /* clip backwards from end */
  while ( (pos1 >= pos0) && strchr ( blank_chars, *pos1 ) )
    pos1 --;

  newlen = pos1 - pos0 + 1;
  if ( !newlen )
    return NULL;

  if ( (ret = LALCalloc(1, newlen + 1)) == NULL )
    return NULL;

  strncpy ( ret, pos0, newlen );
  ret[ newlen ] = 0;

  return ret;

} /* deblank_string() */

