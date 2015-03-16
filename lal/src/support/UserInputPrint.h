//
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

#ifndef _USERINPUTPRINT_H  /* Double-include protection. */
#define _USERINPUTPRINT_H

#include <lal/LALDatatypes.h>

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \defgroup UserInputPrint_h Header UserInputPrint.h
 * \ingroup UserInput_h
 * \author Reinhard Prix
 * \brief Sub-module for general printing of various input 'types' (as defined in \ref UserInput_h) as 'string values',
 * These can be thought of as the 'inverse functions' to \ref ParseStringValue_h.
 *
 */

/*@{*/

// ---------- exported defines & macros ----------

// ----- printing of UTYPES that are mappings to other types
#define XLALPrintStringValueOfRAJ XLALPrintStringValueOfREAL8
#define XLALPrintStringValueOfDECJ XLALPrintStringValueOfREAL8

// ---------- exported types ----------

// ---------- Function prototypes ----------

char *XLALPrintStringValueOfBOOLEAN ( const BOOLEAN *valBOOLEAN );
char *XLALPrintStringValueOfINT4 ( const INT4 *valINT4 );
char *XLALPrintStringValueOfINT8 ( const INT8 *valINT8 );
char *XLALPrintStringValueOfREAL4 ( const REAL4 *valREAL4 );
char *XLALPrintStringValueOfREAL8 ( const REAL8 *valREAL8 );
char *XLALPrintStringValueOfEPOCH ( const LIGOTimeGPS *valGPS );
char *XLALPrintStringValueOfSTRING ( char **valSTRING );
char *XLALPrintStringValueOfSTRINGVector ( LALStringVector **valSTRINGVector );

/*@}*/

/* C++ protection. */
#ifdef  __cplusplus
}
#endif

#endif  /* Double-include protection. */
