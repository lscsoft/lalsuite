/*
*  Copyright (C) 2007 Jolien Creighton
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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

#ifndef GPSTIME_H
#define GPSTIME_H

/*
 *
 * Routines to convert times between LIGOTimeGPS epochs, INT8 nanoseconds,
 * and REAL8 seconds.
 *
 */

#include <lal/LALDatatypes.h>

/* convert from LIGOTimeGPS epoch to INT8 nanoseconds */
INT8 epoch_to_ns( LIGOTimeGPS *epoch );

/* convert from INT8 nanoseconds to LIGOTimeGPS epoch */
LIGOTimeGPS * ns_to_epoch( LIGOTimeGPS *epoch, INT8 ns );

/* convert from REAL8 seconds to INT8 nanoseconds */
INT8 sec_to_ns( REAL8 sec );

/* convert from INT8 nanoseconds to REAL8 seconds */
REAL8 ns_to_sec( INT8 ns );

#endif /* GPSTIME_H */
