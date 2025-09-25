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

#ifndef _USERINPUTPRINT_H  /* Double-include protection. */
#define _USERINPUTPRINT_H

#include <lal/LALDatatypes.h>
#include <lal/UserInputParse.h>

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \defgroup UserInputPrint_h Header UserInputPrint.h
 * \ingroup UserInput_h
 * \author Reinhard Prix
 * \brief Sub-module for general printing of various input 'types' (as defined in \ref UserInput_h) as 'string values',
 * These can be thought of as the 'inverse functions' to \ref UserInputParse_h.
 *
 */

/** @{ */

// ---------- exported defines & macros ----------

// ----- printing of UTYPES that are mappings to other types
#define XLALPrintStringValueOfRAJ XLALPrintStringValueOfREAL8
#define XLALPrintStringValueOfDECJ XLALPrintStringValueOfREAL8
#define XLALPrintStringValueOfRAJRange XLALPrintStringValueOfREAL8Range
#define XLALPrintStringValueOfDECJRange XLALPrintStringValueOfREAL8Range

// ---------- exported types ----------

// ---------- Function prototypes ----------

char *XLALPrintStringValueOfBOOLEAN ( const BOOLEAN *valBOOLEAN );
char *XLALPrintStringValueOfINT4 ( const INT4 *valINT4 );
char *XLALPrintStringValueOfINT8 ( const INT8 *valINT8 );
char *XLALPrintStringValueOfUINT4 ( const UINT4 *valUINT4 );
char *XLALPrintStringValueOfUINT8 ( const UINT8 *valUINT8 );
char *XLALPrintStringValueOfREAL4 ( const REAL4 *valREAL4 );
char *XLALPrintStringValueOfREAL8 ( const REAL8 *valREAL8 );
char *XLALPrintStringValueOfEPOCH ( const LIGOTimeGPS *valGPS );

char *XLALPrintStringValueOfINT4Range ( const INT4Range *int4Range );
char *XLALPrintStringValueOfREAL8Range ( const REAL8Range *real8Range );
char *XLALPrintStringValueOfEPOCHRange ( const LIGOTimeGPSRange *gpsRange );

char *XLALPrintStringValueOfUserEnum ( const int *valEnum, const UserChoices *enumData );
char *XLALFormatHelpStringOfUserEnum ( const UserChoices *enumData );
char *XLALPrintStringValueOfUserFlag ( const int *valFlag, const UserChoices *flagData );
char *XLALFormatHelpStringOfUserFlag ( const UserChoices *flagData );

char *XLALPrintStringValueOfSTRING ( char **valSTRING );
char *XLALPrintStringValueOfSTRINGVector ( LALStringVector **valSTRINGVector );

char *XLALPrintStringValueOfINT4Vector(INT4Vector **valVector);
char *XLALPrintStringValueOfUINT4Vector(UINT4Vector **valVector);
char *XLALPrintStringValueOfREAL8Vector(REAL8Vector **valVector);

/** @} */

/* C++ protection. */
#ifdef  __cplusplus
}
#endif

#endif  /* Double-include protection. */
