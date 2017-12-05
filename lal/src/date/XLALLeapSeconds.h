/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton
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

#ifndef XLALLEAPSECONDS_H
#define XLALLEAPSECONDS_H

#include <lal/LALAtomicDatatypes.h>

/*
 * Leap seconds table
 *
 * JD and GPS time of leap seconds and the value of TAI-UTC.
 *
 * reference: http://maia.usno.navy.mil/
 *            http://maia.usno.navy.mil/ser7/tai-utc.dat
 *
 * notes: the table below must be updated whenever a leap second is added
 * Use the program XLALAddLeapSecond to construct a new table.
 */

static const struct leaps_table { REAL8 jd; INT4 gpssec; int taiutc; } leaps[] =
{
  {2444239.5,    -43200, 19},  /* 1980-Jan-01 */
  {2444786.5,  46828800, 20},  /* 1981-Jul-01 */
  {2445151.5,  78364801, 21},  /* 1982-Jul-01 */
  {2445516.5, 109900802, 22},  /* 1983-Jul-01 */
  {2446247.5, 173059203, 23},  /* 1985-Jul-01 */
#if 0
  /* NOTE: IF THIS WERE A NEGATIVE LEAP SECOND, INSERT AS FOLLOWS */
  {2447161.5, 252028803, 22},  /* 1988-Jan-01 EXAMPLE ONLY! */
#endif
  {2447161.5, 252028804, 24},  /* 1988-Jan-01 */
  {2447892.5, 315187205, 25},  /* 1990-Jan-01 */
  {2448257.5, 346723206, 26},  /* 1991-Jan-01 */
  {2448804.5, 393984007, 27},  /* 1992-Jul-01 */
  {2449169.5, 425520008, 28},  /* 1993-Jul-01 */
  {2449534.5, 457056009, 29},  /* 1994-Jul-01 */
  {2450083.5, 504489610, 30},  /* 1996-Jan-01 */
  {2450630.5, 551750411, 31},  /* 1997-Jul-01 */
  {2451179.5, 599184012, 32},  /* 1999-Jan-01 */
  {2453736.5, 820108813, 33},  /* 2006-Jan-01 */
  {2454832.5, 914803214, 34},  /* 2009-Jan-01 */
  {2456109.5, 1025136015, 35}, /* 2012-Jul-01 */
  {2457204.5, 1119744016, 36}, /* 2015-Jul-01 */
  {2457754.5, 1167264017, 37}, /* 2017-Jan-01 */
};
static const int numleaps = sizeof( leaps ) / sizeof( *leaps );

#endif /* XLALLEAPSECONDS_H */
