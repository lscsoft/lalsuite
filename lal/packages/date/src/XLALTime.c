/*
 * Copyright (C) 2007  Jolien Creighton and Kipp Cannon
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


#include <math.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/Date.h>
#include <lal/XLALError.h>

/**
 * \defgroup XLALTime_c GPS Time
 * \ingroup Date_h
 *
 * \brief GPS time manipulation functions.
 */

/*@{*/

/** Converts GPS time to nano seconds stored as an INT8. */
INT8 XLALGPSToINT8NS( const LIGOTimeGPS *epoch )
{
  return XLAL_BILLION_INT8 * epoch->gpsSeconds + epoch->gpsNanoSeconds;
}


/** Converts nano seconds stored as an INT8 to GPS time. */
LIGOTimeGPS * XLALINT8NSToGPS( LIGOTimeGPS *epoch, INT8 ns )
{
  epoch->gpsSeconds     = ns / XLAL_BILLION_INT8;
  epoch->gpsNanoSeconds = ns % XLAL_BILLION_INT8;
  return epoch;
}


/** Sets GPS time given GPS integer seconds and residual nanoseconds. */
LIGOTimeGPS * XLALGPSSet( LIGOTimeGPS *epoch, INT4 gpssec, INT8 gpsnan )
{
  return XLALINT8NSToGPS( epoch, XLAL_BILLION_INT8 * gpssec + gpsnan );
}


/** Sets GPS time given GPS seconds as a REAL8. */
LIGOTimeGPS * XLALGPSSetREAL8( LIGOTimeGPS *epoch, REAL8 t )
{
  INT4 gpssec = floor(t);
  INT4 gpsnan = floor((t - gpssec) * XLAL_BILLION_REAL8 + 0.5);
  if(isnan(t)) {
    XLALPrintError("%s(): NaN", __func__);
    XLAL_ERROR_NULL(XLAL_EFPINVAL);
  }
  if(fabs(t) > 0x7fffffff) {
    XLALPrintError("%s(): overflow %g", __func__, t);
    XLAL_ERROR_NULL(XLAL_EFPINVAL);
  }
  /* use XLALGPSSet() to normalize the nanoseconds */
  return XLALGPSSet(epoch, gpssec, gpsnan);
}


/** Returns the GPS time as a REAL8. */
REAL8 XLALGPSGetREAL8( const LIGOTimeGPS *epoch )
{
  return epoch->gpsSeconds + (epoch->gpsNanoSeconds / XLAL_BILLION_REAL8);
}

/**
 * Breaks the GPS time into REAL8 integral and fractional parts,
 * each of which has the same sign as the epoch.  Returns the
 * fractional part, and stores the integral part (as a REAL8)
 * in the object pointed to by iptr.  Like the standard C math
 * library function modf().
 */
REAL8 XLALGPSModf( REAL8 *iptr, const LIGOTimeGPS *epoch )
{
  INT8 ns = XLALGPSToINT8NS(epoch);
  INT8 rem; /* remainder */
  *iptr = ns < 0 ? -floor(-ns / XLAL_BILLION_REAL8) : floor(ns / XLAL_BILLION_REAL8);
  rem = ns - ((INT8)(*iptr) * XLAL_BILLION_INT8);
  return (REAL8)(rem) / XLAL_BILLION_REAL8;
}

/*
 * general purpose functions
 */


/** Adds two GPS times. */
LIGOTimeGPS * XLALGPSAddGPS( LIGOTimeGPS *epoch, const LIGOTimeGPS *dt )
{
  return XLALINT8NSToGPS( epoch, XLALGPSToINT8NS( epoch ) + XLALGPSToINT8NS( dt ) );
}


/** Adds a double to a GPS time. */
LIGOTimeGPS * XLALGPSAdd( LIGOTimeGPS *epoch, REAL8 dt )
{
  LIGOTimeGPS dt_gps;
  if(!XLALGPSSetREAL8(&dt_gps, dt))
    XLAL_ERROR_NULL(XLAL_EFUNC);
  return XLALGPSAddGPS(epoch, &dt_gps);
}


/** Difference between two GPS times. */
REAL8 XLALGPSDiff( const LIGOTimeGPS *t1, const LIGOTimeGPS *t0 )
{
  LIGOTimeGPS diff;

  XLALINT8NSToGPS(&diff, XLALGPSToINT8NS( t1 ) - XLALGPSToINT8NS( t0 ));

  return XLALGPSGetREAL8(&diff);
}


/**
 * Compares two GPS times.
 * Returns:
 * - -1 if t0 < t1
 * - 0 if t0 == t1
 * - 1 if t0 > t1.
 * A NULL GPS time is always less than a non-NULL GPS time,
 * and two NULL GPS times are considered equal.
 */
int XLALGPSCmp( const LIGOTimeGPS *t0, const LIGOTimeGPS *t1 )
{
  if ( t0 == NULL || t1 == NULL ) {
    return ( t1 != NULL ) ? -1 : ( ( t0 != NULL ) ? 1 : 0 );
  }
  else {
    INT8 ns0 = XLALGPSToINT8NS( t0 );
    INT8 ns1 = XLALGPSToINT8NS( t1 );
    return ( ns0 > ns1 ) - ( ns0 < ns1 );
  }
}


/**
 * Split a double-precision float into "large" and "small" parts.  The
 * results are stored in the hi and lo arguments, respectively.  hi and lo
 * are such that
 *
 * hi + lo = x
 *
 * exactly, and
 *
 * hi / lo ~= 2^26
 *
 * If x is 0, hi and lo are set to 0.  If x is positive (negative)
 * infinity, hi is set to positive (negative) infinity and lo is set to 0.
 * If x is NaN, hi and lo are set to NaN.
 *
 * No errors occur.
 */

static void split_double(double x, double *hi, double *lo)
{
  *hi = (float) x;
  *lo = isinf(x) ? 0 : x - *hi;
}


/** Multiply a GPS time by a number. */
LIGOTimeGPS *XLALGPSMultiply( LIGOTimeGPS *gps, REAL8 x )
{
  int slo = gps->gpsSeconds % (1<<26);
  int shi = gps->gpsSeconds - slo;
  int nlo = gps->gpsNanoSeconds % (1<<15);
  int nhi = gps->gpsNanoSeconds - nlo;
  double xhi, xlo;
  double addend;
  LIGOTimeGPS gps_addend;

  if(isnan(x) || isinf(x)) {
    XLALPrintError("%s(): invalid multiplicand %g", __func__, x);
    XLAL_ERROR_NULL(XLAL_EFPINVAL);
  }

  split_double(x, &xhi, &xlo);

  /* the count of seconds, the count of nanoseconds, and the multiplicand x
   * have each been split into two parts, a high part and a low part.  from
   * these, there are 8 terms in the product of gps and x.  the 8 terms are
   * added to the result in 6 groups, in increasing order of expected
   * magnitude using Kahan's compensated summation algorithm.  the
   * assumption is that within each group the numerical dynamic range is
   * small enough that the addend can be represented exactly by a
   * double-precision floating point number. */

  addend = nlo * xlo / XLAL_BILLION_REAL8;
  addend -= XLALGPSGetREAL8(XLALGPSSetREAL8(&gps_addend, addend));
  *gps = gps_addend;

  addend += (nlo * xhi + nhi * xlo) / XLAL_BILLION_REAL8;
  addend -= XLALGPSGetREAL8(XLALGPSSetREAL8(&gps_addend, addend));
  XLALGPSAddGPS(gps, &gps_addend);

  addend += nhi * xhi / XLAL_BILLION_REAL8;
  addend -= XLALGPSGetREAL8(XLALGPSSetREAL8(&gps_addend, addend));
  XLALGPSAddGPS(gps, &gps_addend);

  addend += slo * xlo;
  addend -= XLALGPSGetREAL8(XLALGPSSetREAL8(&gps_addend, addend));
  XLALGPSAddGPS(gps, &gps_addend);

  addend += slo * xhi + shi * xlo;
  addend -= XLALGPSGetREAL8(XLALGPSSetREAL8(&gps_addend, addend));
  XLALGPSAddGPS(gps, &gps_addend);

  addend += shi * xhi;
  XLALGPSAddGPS(gps, XLALGPSSetREAL8(&gps_addend, addend));

  return gps;
}


/** Divide a GPS time by a number. */
LIGOTimeGPS *XLALGPSDivide( LIGOTimeGPS *gps, REAL8 x )
{
  LIGOTimeGPS quotient;
  int keep_going;
  double residual;

  if(isnan(x)) {
    XLALPrintError("%s(): NaN", __func__);
    XLAL_ERROR_NULL(XLAL_EFPINVAL);
  }
  if(x == 0) {
    XLALPrintError("%s(): divide by zero", __func__);
    XLAL_ERROR_NULL(XLAL_EFPDIV0);
  }

  /* initial guess */
  XLALGPSSetREAL8(&quotient, XLALGPSGetREAL8(gps) / x);
  /* use Newton's method to iteratively solve for quotient */
  keep_going = 100;
  do {
    LIGOTimeGPS workspace = quotient;
    residual = XLALGPSDiff(gps, XLALGPSMultiply(&workspace, x)) / x;
    XLALGPSAdd(&quotient, residual);
    /* FIXME:  should check for overflow instead of too many iterations;
     * would require error checking to be added to several other arithmetic
     * functions first */
  } while(fabs(residual) > 0.5e-9 && --keep_going);
  *gps = quotient;

  return gps;
}

/*@}*/
