//
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

#ifndef _PARSESTRINGVALUE_H  /* Double-include protection. */
#define _PARSESTRINGVALUE_H

#include <lal/LALDatatypes.h>

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \defgroup ParseStringValue_h Header ParseStringValue.h
 * \ingroup lal_support
 * \author Reinhard Prix
 * \brief Module for general parsing of simple user-input strings
 *
 */

/*@{*/

#ifdef SWIG /* SWIG interface directives */
SWIGLAL(NEW_EMPTY_ARGUMENT(LIGOTimeGPS*, gps));
#endif

// ---------- Function prototypes ----------
int XLALParseStringValueToINT8 ( INT8 *valINT8, const char *valString );
int XLALParseStringValueToINT4 ( INT4 *valINT4, const char *valString );
int XLALParseStringValueToREAL8 ( REAL8 *valREAL8, const char *valString );
int XLALParseStringValueToREAL4 ( REAL4 *valREAL4, const char *valString );
int XLALParseStringValueToBOOLEAN ( BOOLEAN *valBOOLEAN, const char *valString );
int XLALParseStringValueToINT4PlusFrac ( INT4 *valINT4, REAL8 *valFrac, const char *valString );

LIGOTimeGPS *XLALParseStringValueToGPS ( LIGOTimeGPS *gps, const char *valString );
LIGOTimeGPS *XLALParseStringValueToEPOCH ( LIGOTimeGPS *gps, const char *valString );

int XLALParseStringValueToLONGITUDE ( REAL8 *valLongitude, const char *valString );
int XLALParseStringValueToLATITUDE ( REAL8 *valLatitude, const char *valString );

#ifdef SWIG /* SWIG interface directives */
SWIGLAL_CLEAR(NEW_EMPTY_ARGUMENT(LIGOTimeGPS*, gps));
#endif

/*@}*/

/* C++ protection. */
#ifdef  __cplusplus
}
#endif

#endif  /* Double-include protection. */
