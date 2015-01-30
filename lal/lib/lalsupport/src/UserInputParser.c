//
// Copyright (C) 2004, 2005, 2015 Reinhard Prix
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with with program; see the file COPYING. If not, write to the
//  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//  MA  02111-1307  USA
//

#include <errno.h>

#include <lal/StringInput.h>
#include <lal/LALString.h>
#include <lal/UserInputParser.h>


// ---------- local defines ----------
// these constants are taken from StringConvert.c
#define LAL_INT4_MAX    LAL_INT8_C(2147483647)
#define LAL_INT8_MAX    LAL_INT8_C(9223372036854775807)


// ==================== function definitions ====================

//! Parse a string into an INT8
//! This ignores initial whitespace, but throws an error on _any_ non-converted trailing characters (including whitespace)
int
XLALParseStringValueToINT8 ( INT8 *valINT8,         //!< [out] return INT8 value
                             const char *valString  //!< [in]  input string value
                             )
{
  XLAL_CHECK ( (valINT8 != NULL) && (valString != NULL ), XLAL_EINVAL );

  errno = 0;
  char *endptr;
  int base10 = 10;
  long long valLLong = strtoll ( valString, &endptr, base10 );
  XLAL_CHECK ( errno == 0, XLAL_EFAILED, "strtoll() failed to convert '%s' into long long!\n", valString );
  XLAL_CHECK ( (*endptr) == '\0', XLAL_EFAILED, "strtoll(): trailing garbage '%s' found after int-conversion of '%s'\n", endptr, valString );

  //  check range and convert long-int into INT8
  if ( sizeof(valLLong) > sizeof(INT8) ) { // avoid warning about trivial check
    XLAL_CHECK ( (valLLong > -LAL_INT8_MAX) && (valLLong < LAL_INT8_MAX), XLAL_EDOM, "String-conversion '%s' --> '%lli' exceeds INT8 range of +-%"LAL_INT8_FORMAT"\n",
                 valString, valLLong, LAL_INT8_MAX );
  }

  (*valINT8) = (INT8)valLLong;

  return XLAL_SUCCESS;

} // XLALParseStringValueToINT8()


//! Parse a string into an INT4
//! This ignores initial whitespace, but throws an error on _any_ non-converted trailing characters (including whitespace)
int
XLALParseStringValueToINT4 ( INT4 *valINT4,         //!< [out] return INT4 value
                             const char *valString  //!< [in]  input string value
                             )
{
  XLAL_CHECK ( (valINT4 != NULL) && (valString != NULL ), XLAL_EINVAL );

  INT8 valINT8;
  XLAL_CHECK ( XLALParseStringValueToINT8 ( &valINT8, valString ) == XLAL_SUCCESS, XLAL_EFUNC );

  // check range and convert INT8 into INT4
  XLAL_CHECK ( (valINT8 > -LAL_INT4_MAX) && (valINT8 < LAL_INT4_MAX), XLAL_EDOM, "String-conversion '%s' --> '%"LAL_INT8_FORMAT"' exceeds INT4 range of +-%"LAL_INT8_FORMAT"\n",
               valString, valINT8, LAL_INT4_MAX );

  (*valINT4) = (INT4)valINT8;

  return XLAL_SUCCESS;

} // XLALParseStringValueToINT4()


//! Parse a string into a REAL8
//! This ignores initial whitespace, but throws an error on _any_ non-converted trailing characters (including whitespace)
int
XLALParseStringValueToREAL8 ( REAL8 *valREAL8,         //!< [out] return REAL8 value
                              const char *valString    //!< [in]  input string value
                              )
{
  XLAL_CHECK ( (valREAL8 != NULL) && (valString != NULL ), XLAL_EINVAL );

  errno = 0;
  char *endptr;
  double valDouble = strtod ( valString, &endptr );
  XLAL_CHECK ( errno == 0, XLAL_EFAILED, "strtod() failed to convert '%s' into a double!\n", valString );
  XLAL_CHECK ( (*endptr) == '\0', XLAL_EFAILED, "strtod(): trailing garbage '%s' found after double-conversion of '%s'\n", endptr, valString );

  (*valREAL8) = (REAL8)valDouble;

  return XLAL_SUCCESS;

} // XLALParseStringValueToREAL8()


//! Parse a string into a REAL4.
//! This ignores initial whitespace, but throws an error on _any_ non-converted trailing characters (including whitespace)
int
XLALParseStringValueToREAL4 ( REAL4 *valREAL4,         //!< [out] return REAL4 value
                              const char *valString    //!< [in]  input string value
                              )
{
  XLAL_CHECK ( (valREAL4 != NULL) && (valString != NULL ), XLAL_EINVAL );

  errno = 0;
  char *endptr;
  float valFloat = strtof ( valString, &endptr );
  XLAL_CHECK ( errno == 0, XLAL_EFAILED, "strtof() failed to convert '%s' into a float!\n", valString );
  XLAL_CHECK ( (*endptr) == '\0', XLAL_EFAILED, "strtof(): trailing garbage '%s' found after float-conversion of '%s'\n", endptr, valString );

  (*valREAL4) = (REAL4)valFloat;

  return XLAL_SUCCESS;

} // XLALParseStringValueToREAL8()


//! Parse a string into a BOOLEAN
//! Allowed string-values are (case-insensitive):
//! {"yes", "true", "1"} --> TRUE
//! {"no", "false", "0"} --> FALSE
//!
//! NOTE: This throws an error on _any_ extraneous leading or trailing characters or whitespace
int
XLALParseStringValueToBOOLEAN ( BOOLEAN *valBOOLEAN,     //!< [out] return BOOLEAN value
                                const char *valString    //!< [in]  input string value
                                )
{
  XLAL_CHECK ( (valBOOLEAN != NULL) && (valString != NULL ), XLAL_EINVAL );

  /* get rid of case ambiguities */
  char *valStringLower;
  XLAL_CHECK ( (valStringLower = XLALMalloc ( strlen(valString) + 1 )) != NULL, XLAL_ENOMEM );
  strcpy ( valStringLower, valString );
  XLALStringToLowerCase ( valStringLower );

  /* parse it as a bool */
  if ( !strcmp(valStringLower, "yes") || !strcmp(valStringLower, "true") || !strcmp(valStringLower, "1") )
    {
      (*valBOOLEAN) = 1;
    }
  else if ( !strcmp(valStringLower, "no") || !strcmp(valStringLower, "false") || !strcmp(valStringLower, "0") )
    {
      (*valBOOLEAN) = 0;
    }
  else
    {
      XLALFree ( valStringLower );
      XLAL_ERROR ( XLAL_EINVAL, "Illegal bool-string '%s', needs to be one of {'yes', 'true', '1'} or {'no', 'false', '0'} (case-insensitive)\n", valString );
    }

  XLALFree ( valStringLower );

  return XLAL_SUCCESS;

} // XLALParseStringValueToBOOLEAN()
