/** \file
 * \ingroup std
 * \author Cannon, K. C.
 * \brief XLAL routines for computing the sidereal time.
 */


#include <math.h>
#include <lal/Date.h>


/*
 * Returns the Greenwich Sidereal Time IN RADIANS corresponding to a
 * specified GPS time.  Aparent sidereal time is computed by providing the
 * equation of equinoxes in units of seconds.  For mean sidereal time, set
 * this parameter to 0.
 *
 * The output of this code is in radians in the range [0,2pi).
 *
 * Inspired by the function sidereal_time() in the NOVAS-C library, version
 * 2.0.1, which is dated December 10th, 1999, and carries the following
 * references:
 *
 * Aoki, et al. (1982) Astronomy and Astrophysics 105, 359-361.
 * Kaplan, G. H. "NOVAS: Naval Observatory Vector Astrometry
 *   Subroutines"; USNO internal document dated 20 Oct 1988;
 *   revised 15 Mar 1990.
 *
 * See http://aa.usno.navy.mil/software/novas for more information.
 *
 * Note:  rather than maintaining this code separately, it would be a good
 * idea for LAL to simply link to the NOVAS-C library directly.  Something
 * to do when we have some spare time.
 */

REAL8 XLALGreenwichSiderealTime(
	const LIGOTimeGPS *gpstime,
	REAL8 equation_of_equinoxes
)
{
	static const char *func = "XLALGreenwichSiderealTime";

	/* Most significant and least significant part of time since the
	 * Julian epoch in units of centuries (= 36525.0 days).  Original
	 * code in NOVAS-C determined t_hi and t_lo from Julian days.  LAL
	 * lacks a high-precision Julian day routine, so we do this
	 * directly from a GPS time.  Note that the "hi" and "lo"
	 * components have the same units and so the split can be done
	 * anywhere.  The original code recommends putting the integer days
	 * in "hi" and the fractional part in "lo".  For convenience, and
	 * because it involves fewer arithmetic operations, here we make
	 * the split at integer and fractional seconds. */
	double t_hi = (gpstime->gpsSeconds - XLAL_EPOCH_J2000_0_GPS) / (36525.0 * 86400.0);
	double t_lo = gpstime->gpsNanoSeconds / (1e9 * 36525.0 * 86400.0);

	double t = t_hi + t_lo;
	double t_squared = t * t;
	double t_cubed = t_squared * t;

	double sidereal_time = equation_of_equinoxes - 6.2e-6 * t_cubed + 0.093104 * t_squared + 67310.54841
		+ 8640184.812866 * t_lo
		+ 3155760000.0 * t_lo
		+ 8640184.812866 * t_hi
		+ 3155760000.0 * t_hi;

	/* convert from seconds to fraction of 2 pi */
	sidereal_time = fmod(sidereal_time * (2.0 * LAL_PI) / (24.0 * 3600.0), 2.0 * LAL_PI);
	if(sidereal_time < 0.0)
		sidereal_time += 2.0 * LAL_PI;
  
	return sidereal_time;
}


REAL8 XLALGreenwichMeanSiderealTime(
	const LIGOTimeGPS *gpstime
)
{
	return XLALGreenwichSiderealTime(gpstime, 0.0);
}
