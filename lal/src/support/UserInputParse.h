//
// Copyright (C) 2016 Karl Wette
// Copyright (C) 2004, 2015 Reinhard Prix
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

#ifndef _USERINPUTPARSE_H  /* Double-include protection. */
#define _USERINPUTPARSE_H

#include <lal/LALDatatypes.h>

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \defgroup UserInputParse_h Header UserInputParse.h
 * \ingroup UserInput_h
 * \author Reinhard Prix
 * \brief Sub-module for parsing of input 'string values' as various 'types' (as defined in \ref UserInput_h).
 *
 */

/*@{*/

///
/// A range of REAL8 values; first element is minimum, second element is maximum of range
///
typedef REAL8 REAL8Range[2];

///
/// A range of GPS times; first element is minimum, second element is maximum of range
///
typedef LIGOTimeGPS LIGOTimeGPSRange[2];

// ---------- Function prototypes ----------
int XLALParseStringValueAsINT4PlusFrac ( INT4 *valINT4, REAL8 *valFrac, const char *valString );

// --------------- parsers for various USER_TYPE_\<UTYPE\>s ----------
int XLALParseStringValueAsINT8 ( INT8 *valINT8, const char *valString );
int XLALParseStringValueAsINT4 ( INT4 *valINT4, const char *valString );
int XLALParseStringValueAsREAL8 ( REAL8 *valREAL8, const char *valString );
int XLALParseStringValueAsREAL4 ( REAL4 *valREAL4, const char *valString );
int XLALParseStringValueAsBOOLEAN ( BOOLEAN *valBOOLEAN, const char *valString );
int XLALParseStringValueAsGPS ( LIGOTimeGPS *gps, const char *valString );
int XLALParseStringValueAsEPOCH ( LIGOTimeGPS *gps, const char *valString );
int XLALParseStringValueAsRAJ ( REAL8 *valRAJ, const char *valString );
int XLALParseStringValueAsDECJ ( REAL8 *valDECJ, const char *valString );

int XLALParseStringValueAsREAL8Range ( REAL8Range *real8Range, const char *valString );
int XLALParseStringValueAsEPOCHRange ( LIGOTimeGPSRange *gpsRange, const char *valString );
int XLALParseStringValueAsRAJRange ( REAL8Range *rajRange, const char *valString );
int XLALParseStringValueAsDECJRange ( REAL8Range *decjRange, const char *valString );

int XLALParseStringValueAsSTRING ( CHAR **valOut, const char *valString );
int XLALParseStringValueAsSTRINGVector ( LALStringVector **valSTRINGVector, const CHAR *valString );

// use macro templating to define parsers for numerical <CTYPE>vectors
#define DECL_XLALParseStringValueAsVector(CTYPE)                        \
  int XLALParseStringValueAs ##CTYPE## Vector ( CTYPE ## Vector **vect, const CHAR *valString )

DECL_XLALParseStringValueAsVector(REAL8);
DECL_XLALParseStringValueAsVector(INT4);



/*@}*/

/* C++ protection. */
#ifdef  __cplusplus
}
#endif

#endif  /* Double-include protection. */
