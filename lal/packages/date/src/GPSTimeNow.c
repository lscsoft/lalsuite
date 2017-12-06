/*
 * Copyright (C) 2007  Brown, D. A., and Kipp Cannon
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

#include <config.h>
#include <time.h>
#include <string.h>
#include <lal/Date.h>
#include <lal/XLALError.h>

#ifndef HAVE_GMTIME_R
#define gmtime_r(timep, result) memcpy((result), gmtime(timep), sizeof(struct tm))
#endif

/**
 * \ingroup Date_h
 * \brief Populate the LIGOTimeGPS argument with the current system time as
 * returned by time(2) converted to GPS seconds.  Returns the address of
 * the LIGOTimeGPS argument or NULL on error.  On error, the GPS time is
 * undefined.
 *
 * Bugs:
 *
 * This function cannot return negative GPS times.  If the current system
 * time indicates a time prior to Sun Jan 06 00:00:00 GMT 1980, this
 * function returns NULL.
 */
LIGOTimeGPS *
XLALGPSTimeNow (
    LIGOTimeGPS *gpstime
    )
{
  time_t ticks = time(NULL);
  struct tm tm;

  gmtime_r(&ticks, &tm);
  gpstime->gpsSeconds = XLALUTCToGPS(&tm);
  gpstime->gpsNanoSeconds = 0;

  /*
   * XLALUTCToGPS returns < 0 on error, even though of course time did not
   * begin at GPS 0
   */

  if(gpstime->gpsSeconds < 0)
    XLAL_ERROR_NULL(XLAL_EFUNC);

  return gpstime;
}
