/*----------------------------------------------------------------------- 
 * 
 * File Name: date_value.h 
 * 
 * Author: David Chin <dwchin@umich.edu>
 * 
 * Revision: $Id$
 * 
 *----------------------------------------------------------------------- 
 *
 * SYNOPSIS
 * Macros for constants used in date and time routines.
 * 
 * DESCRIPTION
 *
 * 
 * DIAGNOSTICS 
 * (Abnormal termination conditions, error and warning codes summarized 
 * here. More complete descriptions are found in documentation.)
 *
 * CALLS
 * 
 * NOTES
 * 
 *-----------------------------------------------------------------------
 */


/*
 * Macros for constants used in date and time related functions.
 */

/* Fundamental epoch J2000.0 */
#define J2000_0 2451545.0

/* Reference Julian Day for Mean Julian Day/Date */
#define MJDREF  2400000.5

/* Number of days per Julian century */
#define JDAYS_PER_CENT 36525.0

/* Number of seconds per day */
#define SECS_PER_DAY 86400

/* Number of seconds per hour */
#define SECS_PER_HOUR 3600

/* Number of seconds per minute */
#define SECS_PER_MIN 60

/* Number of seconds per degree */
#define SECS_PER_DEG 240

/* Number of minutes per day */
#define MINS_PER_DAY 1440

/* Number of minutes per hour */
#define MINS_PER_HOUR 60

/* Number of hours per day */
#define HRS_PER_DAY 24

/* Number of radians per hour */
#define HRS_PER_RAD (12./LAL_PI)

/* Number of degrees per hour */
#define DEGS_PER_HOUR 15

/* Number of degrees per radian */
#define DEGS_PER_RAD (180./LAL_PI)


/* Difference between Unix and GPS epochs */
/* 
 * difference between Unix and GPS time 315964811 = 
 *     3600 sec/hour x 24 hours/day x (365 days/year x 8 years + 
 *     366 days/year x 2 years + 5 days) + 11 leap seconds 
 * That's roughly speaking.  See docs on UtoGPS.c for exact value
 * quoted from Explanatory Supplement to Astronomical Almanac.
 */
#define UNIXGPS           315964811
#define UNIXGPS_SECS      315964810
#define UNIXGPS_NANOSECS  999918000

/* sidereal hours per mean solar hour (my jargon is probably faulty)
 * (See, _Spherical_Astronomy_, Green, Robin M., 1985, Cambridge;
 * also, more accurate, http://tycho.usno.navy.mil/sidereal.html) */
#define SIDEREAL_PER_TROPICAL 1.00273790935

