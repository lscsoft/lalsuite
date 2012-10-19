/*
*  Copyright (C) 2007 Curt Cutler, Reinhard Prix, Teviet Creighton, Bernd Machenschalk
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#ifndef _LALINITBARYCENTER_H    /* Protect against double-inclusion */
#define _LALINITBARYCENTER_H

#ifdef  __cplusplus
extern "C" {
#endif

#include <lal/LALBarycenter.h>

/** \cond DONT_DOXYGEN */

/* exported API Function prototypes. */
EphemerisData * XLALInitBarycenter ( const CHAR *earthEphemerisFile, const CHAR *sunEphemerisFile );
void XLALDestroyEphemerisData ( EphemerisData *edat );

TimeCorrectionData *XLALInitTimeCorrections ( const CHAR *timeCorrectionFile );
void XLALDestroyTimeCorrectionData( TimeCorrectionData *times );

/* ==================== deprecated LAL interface  ==================================================*/

#define LALINITBARYCENTERH_EOPEN    1
#define LALINITBARYCENTERH_EMEM     2
#define LALINITBARYCENTERH_EEPHFILE 32

#define LALINITBARYCENTERH_MSGEOPEN    "Could not open ephemeris file"
#define LALINITBARYCENTERH_MSGEMEM     "Out of memory"
#define LALINITBARYCENTERH_MSGEEPHFILE "Error in reading an ephemeris file"

void LALInitBarycenter(LALStatus *, EphemerisData *);

/** \endcond */

#ifdef  __cplusplus
}
#endif      /* Close C++ protection */

#endif      /* Close double-include protection */
