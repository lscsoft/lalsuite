/*
*  Copyright (C) 2007 Jolien Creighton, Teviet Creighton
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
#include <lal/LALErrno.h>
#include <lal/XLALError.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/PulsarTimes.h>

/** \addtogroup PulsarTimes_h */
/*@{*/

#define NEQUINOXES 29
/** Define a list of GPS times of autumnal equinoxes (1992 to 2020). */
static const INT4 equinoxes[NEQUINOXES] = {
  401222588, 432778929, 464336350, 495893590, 527450411, 559007772,
  590564232, 622121473, 653678833, 685235053, 716792113, 748349233,
  779905813, 811462993, 843019393, 874576273, 906133453, 937689493,
  969246553, 1000803853, 1032360553, 1063917853, 1095474553,
  1127031613, 1158589273, 1190145733, 1221702853, 1253260213,
  1284816613 };


/**
 * \author Creighton, T. D.
 * \brief Computes the next sidereal midnight and autumnal equinox.
 *
 * This function takes a GPS time from
 * <tt>tepoch</tt> and uses it to assign
 * <tt>tAutumn</tt> and <tt>tMidnight</tt>, which are
 * REAL8 representations of the time in seconds from
 * <tt>tepoch</tt> to the next autumnal equinox or sidereal midnight,
 * respectively.  This routine was written under the \ref PulsarTimes_h
 * module because these quantities are vital for performing pulsar
 * timing: they characterize the Earth's orbital and rotational phase,
 * and hence the Doppler modulation on an incoming signal.  See
 * \ref PulsarTimes_h for more information about the
 * PulsarTimesParamStruc structure.
 *
 * ### Algorithm ###
 *
 * The routine first computes the Greenwich mean sidereal time at
 * <tt>tepoch</tt> using XLALGreenwichMeanSiderealTime(). The next sidereal
 * midnight (at the Prime Meridian) is simply 86400 seconds minus that
 * sidereal time.
 *
 * Next the routine computes the time of the next autumnal equinox.  The
 * module contains an internal list of GPS times of autumnal equinoxes
 * from 1992 to 2020, given to the nearest minute; this is certainly
 * enough accuracy for use with the routines in LALTBaryPtolemaic().
 * If the specified time <tt>tepoch</tt> is after the 2020 autumnal
 * equinox, or more than a year before the 1992 equinox, then the next
 * equinox is extrapolated assuming exact periods of length
 * \ref LAL_YRSID_SI.
 *
 * When assigning the fields of <tt>*times</tt>, it is up to the user to
 * choose a <tt>tepoch</tt> that is close to the actual times that
 * are being considered.  This is important, since many computations use
 * a REAL8 time variable whose origin is the time
 * <tt>tepoch</tt>.  If this is too far from the times of interest,
 * the REAL8 time variables may suffer loss of precision.
 *
 * ### Uses ###
 *
 * \code
 * XLALGreenwichMeanSiderealTime()
 * \endcode
 */
int
XLALGetEarthTimes( const LIGOTimeGPS *tepoch, REAL8 *tMidnight, REAL8 *tAutumn )
{
  LIGOTimeGPS epoch;   /* local copy of times->epoch */
  REAL8 t;             /* time as a floating-point number (s) */

  /* Make sure the parameters exist. */
  XLAL_CHECK( tepoch != NULL, XLAL_EFAULT );
  XLAL_CHECK( tAutumn != NULL, XLAL_EFAULT );
  XLAL_CHECK( tMidnight != NULL, XLAL_EFAULT );
  epoch = *tepoch;

  /* Find the next sidereal midnight. */
  t = fmod(XLALGreenwichMeanSiderealTime(&epoch), LAL_TWOPI) * 86400.0 / LAL_TWOPI;
  XLAL_CHECK( !XLAL_IS_REAL8_FAIL_NAN(t), XLAL_ETIME );
  *tMidnight = 86400.0 - t;

  /* Find the next autumnal equinox. */
  while ( epoch.gpsNanoSeconds > 0 ) {
    epoch.gpsSeconds += 1;
    epoch.gpsNanoSeconds -= 1000000000;
  }
  if ( equinoxes[0] - epoch.gpsSeconds > LAL_YRSID_SI ) {
    t = (REAL8)( equinoxes[0] - epoch.gpsSeconds )
      - (1.0e-9)*epoch.gpsNanoSeconds;
    *tAutumn = fmod( t, LAL_YRSID_SI );
  } else {
    UINT4 i = 0; /* index over equinox list */
    while ( i < NEQUINOXES && equinoxes[i] <= epoch.gpsSeconds )
      i++;
    if ( i == NEQUINOXES ) {
      t = (REAL8)( equinoxes[i-1] - epoch.gpsSeconds )
	- (1.0e-9)*epoch.gpsNanoSeconds;
      *tAutumn = fmod( t, LAL_YRSID_SI ) + LAL_YRSID_SI;
    } else
      *tAutumn = (REAL8)( equinoxes[i] - epoch.gpsSeconds )
	- (1.0e-9)*epoch.gpsNanoSeconds;
  }

  /* Done. */
  return XLAL_SUCCESS;
}

/*** \deprecated use XLALGetEarthTimes() instead.
 */
void
LALGetEarthTimes( LALStatus *status, PulsarTimesParamStruc *times )
{
  if ( XLALGetEarthTimes ( &times->epoch, &times->tMidnight, &times->tAutumn ) != XLAL_SUCCESS )
    ABORT ( status, LAL_FAIL_ERR, LAL_FAIL_MSG );
}

/*@}*/
