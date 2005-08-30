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
};
static const int numleaps = sizeof( leaps ) / sizeof( *leaps );

#endif /* XLALLEAPSECONDS_H */
