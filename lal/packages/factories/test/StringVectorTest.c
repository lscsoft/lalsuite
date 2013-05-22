/*
 * Copyright (C) 2011 Reinhard Prix
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

/*---------- INCLUDES ----------*/
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

#include <lal/LALStdlib.h>
#include <lal/XLALError.h>
#include <lal/StringVector.h>
#include <lal/LALMalloc.h>

/*---------- DEFINES ----------*/
/*---------- internal types ----------*/
/*---------- empty initializers ---------- */
/*---------- Global variables ----------*/
/*---------- internal prototypes ----------*/
int XLALStringVector_TEST ( void );

/**
 * \author Reinhard Prix
 * \file
 * \ingroup StringVector_h
 */

/*==================== FUNCTION DEFINITIONS ====================*/

int main( int argc, char *argv[] )
{

  /* sanity checks */
  if ( argc != 1 )
    XLAL_ERROR ( XLAL_EINVAL, "The executable '%s' doesn't support any input arguments right now.\n", argv[0] );

  if ( XLALStringVector_TEST() != XLAL_SUCCESS )
    XLAL_ERROR ( XLAL_EFUNC, "StringVector TEST failed.\n" );

  return XLAL_SUCCESS;

} /* main() */

/** Test various StringVector functions
 */
int
XLALStringVector_TEST ( void )
{

  LALStringVector *strVect1 = NULL, *strVectCmp = NULL;
#define STR1 "Hello"
#define STR2 "World"
#define STR3 "foo"
#define STR4 "H1"
#define STR5 "H2"
#define STR6 "L1"

  if ( (strVect1 = XLALCreateStringVector ( STR1, STR2, STR3, STR4, STR5, NULL )) == NULL )
    XLAL_ERROR ( XLAL_EFUNC );

  if ( ( strVect1 = XLALAppendString2Vector ( strVect1, STR6 )) == NULL )
    XLAL_ERROR ( XLAL_EFUNC );

  // now sort string-vector according to strcmp()
  if ( XLALSortStringVector ( strVect1 ) != XLAL_SUCCESS ) {
    XLALDestroyStringVector ( strVect1 );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  // generate 'manually' sorted string-vector, using CSV function instead
  // sort-order according to "LC_ALL=C sort"
  if ( ( strVectCmp = XLALParseCSV2StringVector ( STR4 "," STR5 "," STR1 "," STR6 "," STR2 "," STR3 )) == NULL ) {
    XLALDestroyStringVector ( strVect1 );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  /* compare results */
  UINT4 len1 = strVect1->length;
  UINT4 len2 = strVectCmp->length;
  if ( len1 != len2 ) {
    XLALDestroyStringVector ( strVect1  );
    XLALDestroyStringVector ( strVectCmp );
    XLAL_ERROR ( XLAL_EFAILED, "String vectors have different lengths (%d != %d )\n", len1, len2 );
  }
  for ( UINT4 i = 0; i < len1; i++ )
    {
      if ( 0 != strcmp ( strVect1->data[i], strVectCmp->data[i] ) )
        {
          XLALPrintError ( "%s: Sorted string-vector differs from expectation!\n", __func__ );
          for ( UINT4 j=0; j < len1; j ++ )
            XLALPrintError ("j = %d:  s1[j] = %6s, s2[j] = %6s\n", j, strVect1->data[j], strVectCmp->data[j] );
          /* clean up memory, and return with error */
          XLALDestroyStringVector ( strVect1  );
          XLALDestroyStringVector ( strVectCmp );
          XLAL_ERROR ( XLAL_EFAILED, "String sorting failed.\n");
        } /* if s1[i] != s2[i] */
    } /* for i < len */

  /* clean up memory */
  XLALDestroyStringVector ( strVect1  );
  XLALDestroyStringVector ( strVectCmp );

  return XLAL_SUCCESS;

} /* XLALStringVector_TEST() */
