/** \file
 * \ingroup std
 * \author Chin, D. W. and Creighton, J. D. E.
 * \brief XLAL routines for computing the sidereal time.
 *
 */

#include <math.h>
#include <time.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>

/** Returns the Greenwich Mean Sidereal Time IN RADIANS corresponding to a
 * specified GPS time.
 *
 * Note: result is in RADIANS in the range [0,2pi).
 *
 * Reference: S. Aoki et al., A&A 105, 359 (1982) eqs. 13 & 19.
 * Also cf. http://aa.usno.navy.mil.
 */
REAL8 XLALGreenwichMeanSiderealTime(
    const LIGOTimeGPS *epoch /**< [In] GPS time. */
    )
{
  static const char * func = "XLALGreenwichMeanSiderealTime";
  REAL8 t_J2000; /* time since the J2000.0 epoch (2000 JAN 1 12h UTC) */
  REAL8 t;       /* time since 2000 JAN 1 0h UTC */
  REAL8 dpU;     /* days since J2000.0 */
  REAL8 TpU;     /* centuries since 1899 DEC 31 12h UT */
  REAL8 gmst;    /* greenwich mean sidereal time (radians) */
  int taiutc;    /* leap seconds */

  /* compute current number of leap seconds */
  taiutc = XLALLeapSeconds( epoch->gpsSeconds );
  if ( taiutc < 0 )
    XLAL_ERROR_REAL4( func, XLAL_EFUNC );

  t_J2000  = epoch->gpsSeconds - XLAL_EPOCH_J2000_0_GPS;
  t_J2000 += 1e-9 * epoch->gpsNanoSeconds;
  t_J2000 += taiutc - XLAL_EPOCH_J2000_0_TAI_UTC;

  /* 2000 JAN 1 0h UTC is twelve hours earlier than the J2000.0 epoch */
  t = t_J2000 + 12 * 3600;

  /* compute number of days since J2000.0
   * note: this must be an integer-plus-one-half
   * it is really the number of full days since half a day before J2000.0
   * minus one half. */
  dpU  = floor( t / ( 24.0 * 3600.0 ) ); /* full days since 0h UT 01 Jan 2000 */
  dpU -= 0.5; /* recall t is half a day before J2000.0 */

  /* compute number of centuries since 12h UT 31 Dec 1899 */
  TpU = dpU / 36525.0;

  /* compute the gmst at 0h of the current day */
  gmst = 24110.54841
    + TpU * ( 8640184.812866
        + TpU * ( 0.093104
          - TpU * 6.2e-6 ) ); /* seconds */

  /* add the sidereal time since the start of the day */
  t = fmod( t, 24.0 * 3600.0 ); /* seconds since start of day */
  gmst += t * 1.002737909350795; /* corrections omitted */

  /* convert to fractions of a day and to radians */
  gmst = fmod( gmst / ( 24.0 * 3600.0 ), 1.0 ); /* fraction of day */
  gmst *= 2.0 * LAL_PI; /* radians */
  return gmst;
}
