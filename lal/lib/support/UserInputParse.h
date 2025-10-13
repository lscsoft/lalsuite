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
//  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//  MA  02110-1301  USA
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

/** @{ */

///
/// Convert an unsigned long index \c i into a bit, i.e. \f$b = 2^i\f$
///
#define XLAL_IDX2BIT(i) (1UL << ((unsigned long)(i)))

///
/// Convert an unsigned long \e single bit \c b into an index, i.e. \f$i = \log_2 b \f$
///
#define XLAL_BIT2IDX(b) ( \
    (((((unsigned long)(b)) & 0xAAAAAAAAAAAAAAAAUL) != 0)) | \
    (((((unsigned long)(b)) & 0xCCCCCCCCCCCCCCCCUL) != 0) << 1UL) | \
    (((((unsigned long)(b)) & 0xF0F0F0F0F0F0F0F0UL) != 0) << 2UL) | \
    (((((unsigned long)(b)) & 0xFF00FF00FF00FF00UL) != 0) << 3UL) | \
    (((((unsigned long)(b)) & 0xFFFF0000FFFF0000UL) != 0) << 4UL) | \
    (((((unsigned long)(b)) & 0xFFFFFFFF00000000UL) != 0) << 5UL) \
    )

///
/// A range of REAL8 values; first element is minimum, second element is maximum of range
///
typedef REAL8 REAL8Range[2];

///
/// A range of INT4 values; first element is minimum, second element is maximum of range
///
typedef INT4 INT4Range[2];

///
/// A range of GPS times; first element is minimum, second element is maximum of range
///
typedef LIGOTimeGPS LIGOTimeGPSRange[2];

///
/// Possible choices the user may select for an enumeration or bitflag
///
#ifdef SWIG /* SWIG interface directives */
SWIGLAL(IMMUTABLE_MEMBERS(tagUserChoice, name));
#endif /* SWIG */
typedef struct tagUserChoice { int val; const char *name; } UserChoices[32];

// ---------- Function prototypes ----------
int XLALParseStringValueAsINT4PlusFrac ( INT4 *valINT4, REAL8 *valFrac, const char *valString );

// --------------- parsers for various USER_TYPE_\<UTYPE\>s ----------
int XLALParseStringValueAsINT8 ( INT8 *valINT8, const char *valString );
int XLALParseStringValueAsINT4 ( INT4 *valINT4, const char *valString );
int XLALParseStringValueAsUINT8 ( UINT8 *valUINT8, const char *valString );
int XLALParseStringValueAsUINT4 ( UINT4 *valUINT4, const char *valString );
int XLALParseStringValueAsREAL8 ( REAL8 *valREAL8, const char *valString );
int XLALParseStringValueAsREAL4 ( REAL4 *valREAL4, const char *valString );
int XLALParseStringValueAsBOOLEAN ( BOOLEAN *valBOOLEAN, const char *valString );
int XLALParseStringValueAsGPS ( LIGOTimeGPS *gps, const char *valString );
int XLALParseStringValueAsEPOCH ( LIGOTimeGPS *gps, const char *valString );
int XLALParseStringValueAsRAJ ( REAL8 *valRAJ, const char *valString );
int XLALParseStringValueAsDECJ ( REAL8 *valDECJ, const char *valString );

#ifdef SWIG /* SWIG interface directives */
SWIGLAL(OUTPUT_ARRAY_1D(INT4, int4Range[2]));
SWIGLAL(OUTPUT_ARRAY_1D(REAL8, real8Range[2]));
SWIGLAL(OUTPUT_ARRAY_1D(LIGOTimeGPS, gpsRange[2]));
SWIGLAL(OUTPUT_ARRAY_1D(REAL8, rajRange[2]));
SWIGLAL(OUTPUT_ARRAY_1D(REAL8, decjRange[2]));
#endif /* SWIG */

int XLALParseStringValueAsINT4Range ( INT4Range int4Range, const char *valString );
int XLALParseStringValueAsREAL8Range ( REAL8Range real8Range, const char *valString );
int XLALParseStringValueAsEPOCHRange ( LIGOTimeGPSRange gpsRange, const char *valString );
int XLALParseStringValueAsRAJRange ( REAL8Range rajRange, const char *valString );
int XLALParseStringValueAsDECJRange ( REAL8Range decjRange, const char *valString );

int XLALParseStringValueAsUserEnum ( int *valEnum, const UserChoices *enumData, const char *valString );
int XLALParseStringValueAsUserFlag ( int *valFlag, const UserChoices *flagData, const char *valString );

int XLALParseStringValueAsSTRING ( CHAR **valOut, const char *valString );
int XLALParseStringValueAsSTRINGVector ( LALStringVector **valSTRINGVector, const CHAR *valString );

int XLALParseStringValueAsINT4Vector( INT4Vector **vect, const CHAR *valString );
int XLALParseStringValueAsUINT4Vector( UINT4Vector **vect, const CHAR *valString );
int XLALParseStringValueAsREAL8Vector( REAL8Vector **vect, const CHAR *valString );

/** @} */

/* C++ protection. */
#ifdef  __cplusplus
}
#endif

#endif  /* Double-include protection. */
