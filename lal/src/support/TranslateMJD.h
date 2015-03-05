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

#ifndef _TRANSLATEMJD_H  /* Double-include protection. */
#define _TRANSLATEMJD_H

#include <lal/LALDatatypes.h>

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \defgroup TranslateMJD_h Header TranslateMJD.h
 * \ingroup lal_support
 * \author Reinhard Prix
 * \brief Module for translating between MJD(TT) strings and GPS
 *
 */
/*@{*/
// ---------- Function prototypes ----------
#ifdef SWIG /* SWIG interface directives */
SWIGLAL(NEW_EMPTY_ARGUMENT(LIGOTimeGPS*, gps));
#endif

LIGOTimeGPS *XLALTranslateMJDTTtoGPS ( LIGOTimeGPS *gps, INT4 mjdDays, REAL8 mjdFracDays );
LIGOTimeGPS *XLALTranslateStringMJDTTtoGPS ( LIGOTimeGPS *gps, const char *mjdString );

#ifdef SWIG /* SWIG interface directives */
SWIGLAL_CLEAR(NEW_EMPTY_ARGUMENT(LIGOTimeGPS*, gps));
#endif

/*@}*/

/* C++ protection. */
#ifdef  __cplusplus
}
#endif

#endif  /* Double-include protection. */
