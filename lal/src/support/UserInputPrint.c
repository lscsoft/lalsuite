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
//  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//  MA  02111-1307  USA
//

#include <stdio.h>

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
DECL_XLALPrintStringValueOfVector(CTYPE)                                \
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

DEFN_XLALPrintStringValueOfVector(REAL8);
DEFN_XLALPrintStringValueOfVector(INT4);

