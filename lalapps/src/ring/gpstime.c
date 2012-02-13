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
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#include <math.h>
#include <lal/LALStdlib.h>
#include "lalapps.h"
#include "gpstime.h"

/* one second in nano seconds as an INT8 */
#define I8_SEC_NS LAL_INT8_C(1000000000)

/* one second in nano seconds as a REAL8 */
#define D_SEC_NS (1e9)

/* one nano seconds in seconds as a REAL8 */
#define D_NS_SEC (1e-9)


/* express LIGOTimeGPS as INT8 nanoseconds */
INT8 epoch_to_ns( LIGOTimeGPS *epoch )
{
  INT8 ns;
  ns  = (INT8)(epoch->gpsSeconds) * I8_SEC_NS;
  ns += (INT8)(epoch->gpsNanoSeconds);
  return ns;
}


/* convert INT8 nanoseconds to LIGOTimeGPS epoch */
LIGOTimeGPS * ns_to_epoch( LIGOTimeGPS *epoch, INT8 ns )
{
  epoch->gpsSeconds     = (INT4)( ns / I8_SEC_NS );
  epoch->gpsNanoSeconds = (INT4)( ns % I8_SEC_NS );
  return epoch;
}


/* convert REAL8 seconds to INT8 nanoseconds */
INT8 sec_to_ns( REAL8 sec )
{
  INT8 ns;
  ns = (INT8) floor( D_SEC_NS * sec );
  return ns;
}


/* convert INT8 nanoseconds to REAL8 seconds */
REAL8 ns_to_sec( INT8 ns )
{
  REAL8 sec;
  sec = D_NS_SEC * ns;
  return sec;
}
