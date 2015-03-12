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

#include <math.h>
#include <ctype.h>
#include <errno.h>
#include <string.h>

#include <lal/Date.h>
#include <lal/StringInput.h>
#include <lal/StringVector.h>
#include <lal/LALConstants.h>
#include <lal/LALString.h>
#include <lal/TranslateMJD.h>
#include <lal/TranslateAngles.h>

#include <lal/ParseStringValue.h>

// ---------- local defines ----------
// these constants are taken from StringConvert.c
#define LAL_sINT4_MAX    LAL_INT8_C(2147483647)
#define LAL_sINT8_MAX    LAL_INT8_C(9223372036854775807)


// ==================== function definitions ====================

/// Parse a string into an INT8
/// This ignores initial whitespace, but throws an error on _any_ non-converted trailing characters (including whitespace)
int
XLALParseStringValueAsINT8 ( INT8 *valINT8,         ///< [out] return INT8 value
                             const char *valString  ///< [in]  input string value
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
    XLAL_CHECK ( (valLLong >= -LAL_sINT8_MAX) && (valLLong <= LAL_sINT8_MAX), XLAL_EDOM, "String-conversion '%s' --> '%lli' exceeds INT8 range of +-%"LAL_INT8_FORMAT"\n",
                 valString, valLLong, LAL_sINT8_MAX );
  }

  (*valINT8) = (INT8)valLLong;

  return XLAL_SUCCESS;

} // XLALParseStringValueAsINT8()


/// Parse a string into an INT4
/// This ignores initial whitespace, but throws an error on _any_ non-converted trailing characters (including whitespace)
int
XLALParseStringValueAsINT4 ( INT4 *valINT4,         ///< [out] return INT4 value
                             const char *valString  ///< [in]  input string value
                             )
{
  XLAL_CHECK ( (valINT4 != NULL) && (valString != NULL ), XLAL_EINVAL );

  INT8 valINT8;
  XLAL_CHECK ( XLALParseStringValueAsINT8 ( &valINT8, valString ) == XLAL_SUCCESS, XLAL_EFUNC );

  // check range and convert INT8 into INT4
  XLAL_CHECK ( (valINT8 >= -LAL_sINT4_MAX) && (valINT8 <= LAL_sINT4_MAX), XLAL_EDOM, "String-conversion '%s' --> '%"LAL_INT8_FORMAT"' exceeds INT4 range of +-%"LAL_INT8_FORMAT"\n",
               valString, valINT8, LAL_sINT4_MAX );

  (*valINT4) = (INT4)valINT8;

  return XLAL_SUCCESS;

} // XLALParseStringValueAsINT4()


/// Parse a string into a REAL8
/// This ignores initial whitespace, but throws an error on _any_ non-converted trailing characters (including whitespace)
int
XLALParseStringValueAsREAL8 ( REAL8 *valREAL8,         ///< [out] return REAL8 value
                              const char *valString    ///< [in]  input string value
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

} // XLALParseStringValueAsREAL8()

/// Parse a string into a REAL4.
/// This ignores initial whitespace, but throws an error on _any_ non-converted trailing characters (including whitespace)
int
XLALParseStringValueAsREAL4 ( REAL4 *valREAL4,         ///< [out] return REAL4 value
                              const char *valString    ///< [in]  input string value
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

} // XLALParseStringValueAsREAL4()

///
/// Parse a string containing a floating-point number into integer and fractional part, such that val = valINT + valFrac.
///
/// This is useful for parsing strings representing GPS or MJD times wihout loss of ns accuracy.
///
int
XLALParseStringValueAsINT4PlusFrac ( INT4 *valINT4,		///< [out] return INT4 integer part 'xxx'
                                     REAL8 *valFrac,      	///< [out] return fractional part '0.yyyy'
                                     const char *valString	///< [in]  input string value representing a floating-point number "xxx.yyyy"
                                     )
{
  XLAL_CHECK ( (valINT4 != NULL) && (valFrac != NULL) && (valString != NULL), XLAL_EINVAL );
  XLAL_CHECK ( !isspace(valString[0]), XLAL_EINVAL, "No initial whitespace allowed in input string '%s'\n", valString );

  char buf[256];
  strncpy ( buf, valString, sizeof(buf)-1 );
  buf[ sizeof(buf)-1 ] = 0;

  REAL8 sign = 1;
  if ( valString[0] == '-' ) {	 // that's why no initial whitespace is allowed in input string
    sign = -1;
  }
  char *point = strchr ( buf, '.' );
  // is there a fractional part at all? If yes, parse it, if no, set to 0
  if ( point != NULL )
    {
      (*point) = 0;
      char fracString[256] = "0.";
      strcat ( fracString+2, point+1);
      XLAL_CHECK ( XLALParseStringValueAsREAL8 ( valFrac, fracString ) == XLAL_SUCCESS, XLAL_EFUNC );
      (*valFrac) *= sign;	// correct sign: must agree with integer part
    }
  else
    {
      (*valFrac) = 0;
    }

  // now parse integer part
  XLAL_CHECK ( XLALParseStringValueAsINT4 ( valINT4, buf ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;
} // XLALParseStringValueAsINT4PlusFrac()


/// Parse a string into a BOOLEAN
/// Allowed string-values are (case-insensitive):
/// {"yes", "true", "1"} --> TRUE
/// {"no", "false", "0"} --> FALSE
///
/// NOTE: This throws an error on _any_ extraneous leading or trailing characters or whitespace
int
XLALParseStringValueAsBOOLEAN ( BOOLEAN *valBOOLEAN,     ///< [out] return BOOLEAN value
                                const char *valString    ///< [in]  input string value
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

} // XLALParseStringValueAsBOOLEAN()


///
/// Parse a string representing a GPS time into LIGOTimeGPS, without loss of (ns) accuracy
///
///
int
XLALParseStringValueAsGPS ( LIGOTimeGPS *gps,	///< [out] returned GPS time
                            const char *valString 	///< [in] input string representing MJD(TT) time
                            )
{
  XLAL_CHECK ( (gps != NULL) && (valString != NULL), XLAL_EINVAL );

  INT4 gpsInt;
  REAL8 gpsFrac;
  XLAL_CHECK ( XLALParseStringValueAsINT4PlusFrac ( &gpsInt, &gpsFrac, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  INT8 gpsNs = (INT8) round ( gpsFrac * XLAL_BILLION_REAL8 );

  XLALGPSSet ( gps, gpsInt, gpsNs );

  return XLAL_SUCCESS;

} // XLALParseStringValueAsGPS()


///
/// Parse a string representing an 'epoch' into an LIGOTimeGPS, allowing both GPS and MJD(TT) inputs, at ns accuracy.
///
/// Allowed input string formats are "INT4.INT4[GPS|MJD]", where the optional postfix 'GPS' or 'MJD' indicates
/// the time 'units', which defaults to GPS if no postfix is given. Note that MJD input is interpreted as MJD(TT),
/// and translated into GPS using XLALTranslateStringMJDTTtoGPS().
/// This ignores initial whitespace, but throws an error on _any_ non-converted trailing characters (including whitespace)
///
int
XLALParseStringValueAsEPOCH ( LIGOTimeGPS *gps,   	///< [out] return LIGOTimeGPS value
                              const char *valString  	///< [in]  input string value
                              )
{
  XLAL_CHECK ( (gps != NULL) && (valString != NULL ), XLAL_EINVAL );

  char buf[256];
  strncpy ( buf, valString, sizeof(buf)-1 );
  buf[ sizeof(buf)-1 ] = 0;

  // ---------- first check if there's a postfix indicating the time 'units' (GPS or MJD):
  BOOLEAN is_gps;
  char *postfix;
  if ( (postfix = strstr ( buf, "MJD" )) != NULL )
    {
      XLAL_CHECK ( postfix[3] == 0, XLAL_EINVAL, "Input '%s' contains trailing characters after units 'MJD': must be of form 'xxx.yyyMJD'\n", valString );
      postfix[0] = 0; // cut off postfix
      is_gps = 0;
    }
  else if ( (postfix = strstr ( buf, "GPS" )) != NULL )
    {
      XLAL_CHECK ( postfix[3] == 0, XLAL_EINVAL, "Input '%s' contains trailing characters after units 'GPS': must be of form 'xxx.yyy' or 'xxx.yyyGPS'\n", valString );
      postfix[0] = 0; // cut off postfix
      is_gps = 1;
    }
  else	// no postfix: default to 'GPS' units
    {
      is_gps = 1;
    }


  if ( is_gps )
    {
      XLAL_CHECK ( XLALParseStringValueAsGPS ( gps, buf ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  else
    {
      XLAL_CHECK ( XLALTranslateStringMJDTTtoGPS ( gps, buf ) != NULL, XLAL_EFUNC );
    }

  return XLAL_SUCCESS;

} // XLALParseStringValueAsEPOCH()


///
/// Parse a string representing an 'equatorial longitude' (aka right ascension or RA) into REAL8 radians, allowing for both radians or "hours:minutes:seconds" as input.
///
/// Note that "h:m:s" input is translated into radians using XLALTranslateHMStoRAD().
///
int
XLALParseStringValueAsRAJ ( REAL8 *valRAJ,   	///< [out] return longitude value in radians
                            const char *valString  	///< [in]  input string value
                            )
{
  XLAL_CHECK ( (valRAJ != NULL) && (valString != NULL ), XLAL_EINVAL );

  // ---------- first check if there's a colon ':' somewhere in the string, which indicates "H:M:S" format
  const char *colon;
  if ( (colon = strchr ( valString, ':' )) != NULL )
    {
      XLAL_CHECK ( XLALTranslateHMStoRAD ( valRAJ, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  else
    {
      XLAL_CHECK ( XLALParseStringValueAsREAL8 ( valRAJ, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
    }

  return XLAL_SUCCESS;

} // XLALParseStringValueAsRAJ()

///
/// Parse a string representing an 'equatorial latitude' (aka declination or DEC) into REAL8 radians, allowing for both radians or "degrees:minutes:seconds" as input.
///
/// Note that "d:m:s" input is translated into radians using XLALTranslateDMStoRAD().
///
int
XLALParseStringValueAsDECJ ( REAL8 *valDECJ,   	///< [out] return latitude value in radians
                             const char *valString 	///< [in]  input string value
                             )
{
  XLAL_CHECK ( (valDECJ != NULL) && (valString != NULL ), XLAL_EINVAL );

  // ---------- first check if there's a colon ':' somewhere in the string, which indicates "D:M:S" format
  const char *colon;
  if ( (colon = strchr ( valString, ':' )) != NULL )
    {
      XLAL_CHECK ( XLALTranslateDMStoRAD ( valDECJ, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  else
    {
      XLAL_CHECK ( XLALParseStringValueAsREAL8 ( valDECJ, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
    }

  return XLAL_SUCCESS;

} // XLALParseStringValueAsDECJ()

///
/// Duplicate string 'in', removing surrounding quotes \" or \' if present.
///
/// \note Quotes at the beginning of the string must be matched at the end,
/// otherwise we return an error.
///
/// \note The output string (*out) can be !=NULL, in which case it is freed first!
///
int
XLALParseStringValueAsSTRING ( CHAR **out,		///< [out] return allocated string
                               const CHAR *valStr	///< [in] input string value
                               )
{
  XLAL_CHECK ( valStr != NULL, XLAL_EINVAL );
  XLAL_CHECK ( out != NULL, XLAL_EINVAL );

  CHAR opening_quote = 0;
  CHAR closing_quote = 0;
  UINT4 inlen = strlen ( valStr );

  if ( (valStr[0] == '\'') || (valStr[0] == '\"') ) {
    opening_quote = valStr[0];
  }
  if ( (inlen >= 2) && ( (valStr[inlen-1] == '\'') || (valStr[inlen-1] == '\"') ) ) {
    closing_quote = valStr[inlen-1];
  }

  // check matching quotes
  XLAL_CHECK ( opening_quote == closing_quote, XLAL_EINVAL, "Unmatched quotes in string [%s]\n", valStr );

  const CHAR *start = valStr;
  UINT4 outlen = inlen;
  if ( opening_quote )
    {
      start = valStr + 1;
      outlen = inlen - 2;
    }

  CHAR *ret;
  XLAL_CHECK ( (ret = XLALCalloc (1, outlen + 1)) != NULL, XLAL_ENOMEM );
  strncpy ( ret, start, outlen );
  ret[outlen] = 0;

  XLALFree ( (*out) );	// free any previously-allocated string in output variable
  (*out) = ret;

  return XLAL_SUCCESS;

} // XLALParseStringValueAsSTRING()


///
/// Parse a string containing a list of comma-separated values (CSV) into a StringVector.
/// \note surrounding whitespace is removed from the individual list entries.
///
/// \note The output string-vector (*strVect) can be !=NULL, in which case it is freed first!
///
int
XLALParseStringValueAsSTRINGVector ( LALStringVector **strVect,	///< [out] allocated string vector
                             const CHAR *valString	///< [in] input string value
                             )
{
  XLAL_CHECK ( valString != NULL, XLAL_EINVAL );
  XLAL_CHECK ( strVect != NULL, XLAL_EINVAL );

  LALStringVector *ret;
  XLAL_CHECK ( (ret = XLALCalloc ( 1, sizeof(*ret))) != NULL, XLAL_ENOMEM );

  const char *start = valString;
  const char *tmp;
  do
    {
      // create space for the next string-vector entry
      ret->length ++;
      XLAL_CHECK ( (ret->data = XLALRealloc ( ret->data, ret->length * sizeof(ret->data[0]) )) != NULL, XLAL_ENOMEM );

      // determine length of next CSV string value
      size_t len;
      if ( ( tmp = strchr ( start, ',' ) ) ) {
	len = tmp - start;
      } else {
	len = strlen ( start );
      }

      // copy this value with surrounding whitespace removed
      XLAL_CHECK ( (ret->data[ret->length-1] = XLALDeblankString ( start, len ) ) != NULL, XLAL_EFUNC );

    } while ( (tmp != NULL) && ((start = tmp + 1) != NULL) );

  XLALDestroyStringVector ( (*strVect) );	// free any previously-allocated string-vector in output variable
  (*strVect) = ret;

  return XLAL_SUCCESS;

} // XLALParseStringValueAsSTRINGVector()

