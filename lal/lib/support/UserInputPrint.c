//
// Copyright (C) 2016 Karl Wette
// Copyright (C) 2015 Reinhard Prix
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
//  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//  MA  02110-1301  USA
//

#include <stdio.h>
#include <string.h>

#include <lal/LALMalloc.h>
#include <lal/LALStdio.h>
#include <lal/LALString.h>
#include <lal/XLALError.h>
#include <lal/Date.h>
#include <lal/LALDatatypes.h>

#include <lal/UserInputPrint.h>

// ==================== function definitions ====================

/// Return 'string value' (allocated here) of an INT8
char *
XLALPrintStringValueOfINT8 ( const INT8 *valINT8 )
{
  XLAL_CHECK_NULL ( valINT8 != NULL, XLAL_EINVAL );
  char buf[256];
  sprintf ( buf, "%" LAL_INT8_FORMAT, (*valINT8) );
  return XLALStringDuplicate ( buf );
}

/// Return 'string value' (allocated here) of an INT4
char *
XLALPrintStringValueOfINT4 ( const INT4 *valINT4 )
{
  XLAL_CHECK_NULL ( valINT4 != NULL, XLAL_EINVAL );
  char buf[256];
  sprintf ( buf, "%" LAL_INT4_FORMAT, (*valINT4) );
  return XLALStringDuplicate ( buf );
}

/// Return 'string value' (allocated here) of an UINT8
char *
XLALPrintStringValueOfUINT8 ( const UINT8 *valUINT8 )
{
  XLAL_CHECK_NULL ( valUINT8 != NULL, XLAL_EINVAL );
  char buf[256];
  sprintf ( buf, "%" LAL_UINT8_FORMAT, (*valUINT8) );
  return XLALStringDuplicate ( buf );
}

/// Return 'string value' (allocated here) of an UINT4
char *
XLALPrintStringValueOfUINT4 ( const UINT4 *valUINT4 )
{
  XLAL_CHECK_NULL ( valUINT4 != NULL, XLAL_EINVAL );
  char buf[256];
  sprintf ( buf, "%" LAL_UINT4_FORMAT, (*valUINT4) );
  return XLALStringDuplicate ( buf );
}

/// Return 'string value' (allocated here) of a REAL8 (printed at full precision)
char *
XLALPrintStringValueOfREAL8 ( const REAL8 *valREAL8 )
{
  XLAL_CHECK_NULL ( valREAL8 != NULL, XLAL_EINVAL );
  char buf[256];
  sprintf ( buf, "%.16" LAL_REAL8_FORMAT, (*valREAL8) );
  return XLALStringDuplicate ( buf );
}

/// Return 'string value' (allocated here) of a REAL4 (printed at full precision)
char *
XLALPrintStringValueOfREAL4 ( const REAL4 *valREAL4 )
{
  XLAL_CHECK_NULL ( valREAL4 != NULL, XLAL_EINVAL );
  char buf[256];
  sprintf ( buf, "%.8" LAL_REAL4_FORMAT, (*valREAL4) );
  return XLALStringDuplicate ( buf );
}

/// Return 'string value' (allocated here) of a BOOLEAN
char *
XLALPrintStringValueOfBOOLEAN ( const BOOLEAN *valBOOLEAN )
{
  XLAL_CHECK_NULL ( valBOOLEAN != NULL, XLAL_EINVAL );
  char buf[256];
  sprintf ( buf, (*valBOOLEAN) ? "TRUE" : "FALSE" );
  return XLALStringDuplicate ( buf );
}


/// Return 'string value' (allocated here) of a GPS epoch, as parseable by XLALParseStringValueAs[GPS|EPOCH]()
char *
XLALPrintStringValueOfEPOCH ( const LIGOTimeGPS *valGPS )
{
  XLAL_CHECK_NULL ( valGPS != NULL, XLAL_EINVAL );
  char buf[256];
  XLAL_CHECK_NULL ( XLALGPSToStr ( buf, valGPS ) != NULL, XLAL_EFUNC );
  return XLALStringDuplicate ( buf );
}

/// Return 'string value' (allocated here) of a STRING, surrounded by double quotes.
/// The output is parseable by XLALParseStringValueAsSTRING().
/// In case of NULL input, returns the string 'NULL'.
char *
XLALPrintStringValueOfSTRING ( CHAR **valSTRING )
{
  XLAL_CHECK_NULL ( valSTRING != NULL, XLAL_EINVAL );

  char *ret;
  if ( (*valSTRING) == NULL )
    {
      XLAL_CHECK_NULL ( (ret = XLALStringDuplicate("NULL")) != NULL, XLAL_EFUNC );
    }
  else
    {
      XLAL_CHECK_NULL ( (ret = XLALMalloc ( strlen((*valSTRING)) + 2 + 1)) != NULL, XLAL_ENOMEM );	// +surrounding quotes + terminating-0
      sprintf ( ret, "\"%s\"", (*valSTRING) );
    }

  return ret;

} // XLALPrintStringValueOfSTRING()

/// Return 'string value' (allocated here) of a INT4Range, as parseable by XLALParseStringValueAsINT4Range()
char *
XLALPrintStringValueOfINT4Range(const INT4Range *int4Range) {
  XLAL_CHECK_NULL(int4Range != NULL, XLAL_EFAULT);
  char *part0 = XLALPrintStringValueOfINT4(&(*int4Range)[0]);
  XLAL_CHECK_NULL(part0 != NULL, XLAL_EFUNC);
  char *part1 = XLALPrintStringValueOfINT4(&(*int4Range)[1]);
  XLAL_CHECK_NULL(part1 != NULL, XLAL_EFUNC);
  char *retn = XLALMalloc(sizeof(*retn) * (strlen(part0) + 1 + strlen(part1) + 1));
  XLAL_CHECK_NULL(retn != NULL, XLAL_ENOMEM);
  strcpy(retn, part0);
  strcat(retn, ",");
  strcat(retn, part1);
  XLALFree(part0);
  XLALFree(part1);
  return retn;
}


/// Return 'string value' (allocated here) of a REAL8Range, as parseable by XLALParseStringValueAsREAL8Range()
char *
XLALPrintStringValueOfREAL8Range(const REAL8Range *real8Range) {
  XLAL_CHECK_NULL(real8Range != NULL, XLAL_EFAULT);
  char *part0 = XLALPrintStringValueOfREAL8(&(*real8Range)[0]);
  XLAL_CHECK_NULL(part0 != NULL, XLAL_EFUNC);
  char *part1 = XLALPrintStringValueOfREAL8(&(*real8Range)[1]);
  XLAL_CHECK_NULL(part1 != NULL, XLAL_EFUNC);
  char *retn = XLALMalloc(sizeof(*retn) * (strlen(part0) + 1 + strlen(part1) + 1));
  XLAL_CHECK_NULL(retn != NULL, XLAL_ENOMEM);
  strcpy(retn, part0);
  strcat(retn, ",");
  strcat(retn, part1);
  XLALFree(part0);
  XLALFree(part1);
  return retn;
}


/// Return 'string value' (allocated here) of a EPOCHRange, as parseable by XLALParseStringValueAsEPOCHRange()
char *
XLALPrintStringValueOfEPOCHRange(const LIGOTimeGPSRange *gpsRange) {
  XLAL_CHECK_NULL(gpsRange != NULL, XLAL_EFAULT);
  char *part0 = XLALPrintStringValueOfEPOCH(&(*gpsRange)[0]);
  XLAL_CHECK_NULL(part0 != NULL, XLAL_EFUNC);
  char *part1 = XLALPrintStringValueOfEPOCH(&(*gpsRange)[1]);
  XLAL_CHECK_NULL(part1 != NULL, XLAL_EFUNC);
  char *retn = XLALMalloc(sizeof(*retn) * (strlen(part0) + 1 + strlen(part1) + 1));
  XLAL_CHECK_NULL(retn != NULL, XLAL_ENOMEM);
  strcpy(retn, part0);
  strcat(retn, ",");
  strcat(retn, part1);
  XLALFree(part0);
  XLALFree(part1);
  return retn;
}


/// Return 'string value' (allocated here) of a user selection of an enumeration value.
/// The output is parseable by XLALParseStringValueOfUserEnum().
char*
XLALPrintStringValueOfUserEnum ( const int *valEnum, const UserChoices *enumData )
{

  // Check input
  XLAL_CHECK_NULL(valEnum != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(*valEnum >= 0, XLAL_EINVAL);
  XLAL_CHECK_NULL(enumData != NULL, XLAL_EFAULT);

  // Select name of enumeration value
  for ( size_t i = 0; i < XLAL_NUM_ELEM(*enumData); ++i ) {
    if ((*enumData)[i].val >= 0 && (*enumData)[i].name != NULL) {
      if (*valEnum == (*enumData)[i].val) {
        return XLALStringDuplicate((*enumData)[i].name);
      }
    }
  }

  XLAL_ERROR_NULL( XLAL_EINVAL, "Value '%i' is not a valid enumeration value", *valEnum );

}


/// Return format help string (allocated here) for a user selection of an enumeration value.
char*
XLALFormatHelpStringOfUserEnum ( const UserChoices *enumData )
{

  // Check input
  XLAL_CHECK_NULL(enumData != NULL, XLAL_EFAULT);

  // Generate format help string
  CHAR *str = NULL;
  int prev_val = -1;
  const char *prev_name = NULL;
  for ( size_t i = 0; i < XLAL_NUM_ELEM(*enumData); ++i ) {
    if ((*enumData)[i].val >= 0 && (*enumData)[i].name != NULL) {
      if (prev_name != NULL && strcmp(prev_name, (*enumData)[i].name) == 0) {
        continue;
      }
      if (str == NULL) {
        str = XLALStringAppendFmt( str, "=(%s", (*enumData)[i].name);
      } else if (prev_val >= 0 && prev_val == (*enumData)[i].val) {
        str = XLALStringAppendFmt( str, "=%s", (*enumData)[i].name);
      } else {
        str = XLALStringAppendFmt( str, "|%s", (*enumData)[i].name);
      }
      XLAL_CHECK_NULL(str != NULL, XLAL_EFUNC);
      prev_val = (*enumData)[i].val;
      prev_name = (*enumData)[i].name;
    }
  }
  str = XLALStringAppendFmt( str, ")");
  XLAL_CHECK_NULL(str != NULL, XLAL_EFUNC);

  return str;

}


/// Return 'string value' (allocated here) of a user selection of an bitflag value.
/// The output is parseable by XLALParseStringValueOfUserFlag().
char*
XLALPrintStringValueOfUserFlag ( const int *valFlag, const UserChoices *flagData )
{

  // Check input
  XLAL_CHECK_NULL(valFlag != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(*valFlag >= 0, XLAL_EINVAL);
  XLAL_CHECK_NULL(flagData != NULL, XLAL_EFAULT);

  // Handle special case of first bitflag with value 0, representing a zero bitflag
  if ((*flagData)[0].val == 0 && (*flagData)[0].name != NULL && *valFlag == 0) {
    return XLALStringDuplicate( (*flagData)[0].name );
  }

  // Deduce bitflag value
  int val = *valFlag;
  CHAR *str = NULL;
  for ( size_t i = 0; i < XLAL_NUM_ELEM(*flagData); ++i ) {
    if ((*flagData)[i].val > 0 && (*flagData)[i].name != NULL) {
      if (val & (*flagData)[i].val) {
        val &= ~(*flagData)[i].val;
        str = XLALStringAppendFmt( str, "%s%s", ( str == NULL ? "" : "," ), (*flagData)[i].name);
        XLAL_CHECK_NULL(str != NULL, XLAL_EFUNC);
      }
    }
  }
  XLAL_CHECK_NULL( val == 0, XLAL_EINVAL, "Value '%i' is not a valid bitflag value", *valFlag );

  return str;

}


/// Return format help string (allocated here) for a user selection of an bitflag value.
char*
XLALFormatHelpStringOfUserFlag ( const UserChoices *flagData )
{

  // Check input
  XLAL_CHECK_NULL(flagData != NULL, XLAL_EFAULT);

  // Generate format help string
  CHAR *str = NULL;
  int prev_val = 0;
  const char *prev_name = NULL;
  for ( size_t i = 0; i < XLAL_NUM_ELEM(*flagData); ++i ) {
    if ((*flagData)[i].val > 0 && (*flagData)[i].name != NULL) {
      if (prev_name != NULL && strcmp(prev_name, (*flagData)[i].name) == 0) {
        continue;
      }
      if (str == NULL) {
        str = XLALStringAppendFmt( str, "=(%s", (*flagData)[i].name);
      } else if (prev_val > 0 && prev_val == (*flagData)[i].val ) {
        str = XLALStringAppendFmt( str, "=%s", (*flagData)[i].name);
      } else {
        str = XLALStringAppendFmt( str, "|%s", (*flagData)[i].name);
      }
      XLAL_CHECK_NULL(str != NULL, XLAL_EFUNC);
      prev_val = (*flagData)[i].val;
      prev_name = (*flagData)[i].name;
    }
  }
  str = XLALStringAppendFmt( str, ")[,...]");
  XLAL_CHECK_NULL(str != NULL, XLAL_EFUNC);
  if ((*flagData)[0].val == 0 && (*flagData)[0].name != NULL) {
    CHAR *old_str = str;
    str = XLALStringAppendFmt( NULL, "=%s|(%s)", (*flagData)[0].name, old_str + 1);
    XLAL_CHECK_NULL(str != NULL, XLAL_EFUNC);
    XLALFree(old_str);
  }

  return str;

}


/// Return 'string value' (allocated here) of a STRINGVector, by turning into comma-separated list of strings, each surrounded by single quotes.
/// The output is parseable by XLALParseStringValueAsSTRINGVector().
/// In case of a NULL or empty vector (data==NULL|length==0), generate the string 'NULL'.
char *
XLALPrintStringValueOfSTRINGVector ( LALStringVector **valSTRINGVector )
{
  XLAL_CHECK_NULL ( valSTRINGVector != NULL, XLAL_EINVAL );
  char *ret = NULL;
  if ( (*valSTRINGVector == NULL) || ((*valSTRINGVector)->data == NULL) || ((*valSTRINGVector)->length == 0) )
    {
      XLAL_CHECK_NULL ( (ret = XLALStringDuplicate("NULL")) != NULL, XLAL_EFUNC );
    }
  else
    {
      for ( UINT4 i=0; i < (*valSTRINGVector)->length; i++ )
        {
          if ( i != 0 ) {
            XLAL_CHECK_NULL ( (ret = XLALStringAppend ( ret, "," )) != NULL, XLAL_EFUNC );
          }
          char *stringi;
          XLAL_CHECK_NULL ( (stringi = XLALPrintStringValueOfSTRING ( &((*valSTRINGVector)->data[i]) )) != NULL, XLAL_EFUNC );
          XLAL_CHECK_NULL ( (ret = XLALStringAppend ( ret, stringi )) != NULL, XLAL_EFUNC );
          XLALFree ( stringi );
        } // for i < length
    } // end: if valSTRING != NULL

  return ret;

} // XLALPrintStringValueOfSTRINGVector()

///
/// Return 'string value' (allocated here) of a \<CTYPE\>Vector, by turning into comma-separated list of \<CTYPE\> values.
/// The output is parseable by XLALParseStringValueAs\<CTYPE\>Vector().
/// In case of a NULL or empty vector (data==NULL|length==0), generate the string 'NULL'.
#define DEFN_XLALPrintStringValueOfVector(CTYPE)                        \
char *XLALPrintStringValueOf##CTYPE##Vector ( CTYPE##Vector **valVector ) \
{                                                                       \
 XLAL_CHECK_NULL ( valVector != NULL, XLAL_EINVAL );                    \
 char *ret = NULL;                                                      \
 if ( (*valVector == NULL) || ((*valVector)->data == NULL) || ((*valVector)->length == 0) ) \
   {                                                                    \
     XLAL_CHECK_NULL ( (ret = XLALStringDuplicate("NULL")) != NULL, XLAL_EFUNC ); \
   }                                                                    \
 else                                                                   \
   {                                                                    \
     for ( UINT4 i=0; i < (*valVector)->length; i++ )                   \
       {                                                                \
         if ( i != 0 ) {                                                \
           XLAL_CHECK_NULL ( (ret = XLALStringAppend ( ret, "," )) != NULL, XLAL_EFUNC ); \
         }                                                              \
         char *tmp;                                                     \
         XLAL_CHECK_NULL ( (tmp = XLALPrintStringValueOf##CTYPE ( &(*valVector)->data[i])) != NULL, XLAL_EFUNC ); \
         XLAL_CHECK_NULL ( (ret = XLALStringAppend ( ret, tmp )) != NULL, XLAL_EFUNC ); \
         XLALFree ( tmp );                                              \
       } /* for i < length */                                           \
   } /* end: if non-empty input vector */                               \
                                                                        \
 return ret;                                                            \
                                                                        \
} /* XLALPrintStringValueOf<CYPTE>Vector() */

DEFN_XLALPrintStringValueOfVector(INT4);
DEFN_XLALPrintStringValueOfVector(UINT4);
DEFN_XLALPrintStringValueOfVector(REAL8);
