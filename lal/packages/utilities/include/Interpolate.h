/*
*  Copyright (C) 2007 Jolien Creighton, Drew Keppel
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

#ifndef _INTERPOLATE_H
#define _INTERPOLATE_H

#include <lal/Date.h>
#include <lal/Sequence.h>
#include <lal/LALDatatypes.h>

#include <lal/XLALGSL.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * \addtogroup Interpolate_h
 *
 * \brief This header covers the routines for interpolation.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/Interpolate.h>
 * \endcode
 *
 * ### Description ###
 *
 * The routine <tt>LALSPolynomialInterpolation()</tt> computes the interpolated \f$y\f$
 * value \c output at the \f$x\f$ value \c target by fitting a polynomial of
 * order <tt>params.n-1</tt> to the data.  The result \c output is of type
 * \c SInterpolateOut, which contains the value <tt>output.y</tt> as well as
 * an estimate of the error <tt>output.dy</tt>.  The routine
 * <tt>LALDPolynomialInterpolation()</tt> is the same but for double precision.
 *
 * ### Operating Instructions ###
 *
 * The following program fits a fourth-order polynomial to the five data points
 * \f$\{(0,0),(1,1),(2,3),(3,4),(4,3)\}\f$, and interpolates the value at \f$x=2.4\f$.
 *
 * \code
 * #include <lal/LALStdlib.h>
 * #include <lal/Interpolate.h>
 *
 * int main ()
 * {
 * enum { ArraySize = 5 };
 * static LALStatus status;
 * REAL4            x[ArraySize] = {0,1,2,3,4};
 * REAL4            y[ArraySize] = {0,1,3,4,3};
 * REAL4            target       = 2.4;
 * SInterpolatePar  intpar       = {ArraySize, x, y};
 * SInterpolateOut  intout;
 *
 * LALSPolynomialInterpolation( &status, &intout, target, &intpar );
 *
 * return 0;
 * }
 * \endcode
 *
 * ### Algorithm ###
 *
 * This is an implementation of the Neville algroithm, see \c polint in
 * Numerical Recipes \cite ptvf1992.
 *
 */
/*@{*/

/** \name Error Codes */
/*@{*/
#define INTERPOLATEH_ENULL 1		/**< Null pointer */
#define INTERPOLATEH_ESIZE 2		/**< Invalid size */
#define INTERPOLATEH_EZERO 4		/**< Zero divide */
/*@}*/

/** \cond DONT_DOXYGEN */
#define INTERPOLATEH_MSGENULL "Null pointer"
#define INTERPOLATEH_MSGESIZE "Invalid size"
#define INTERPOLATEH_MSGEZERO "Zero divide"
/** \endcond */

/** These structures contain the output of the interpolation */
typedef struct
tagSInterpolateOut
{
  REAL4  y;	/**< The interpolated value */
  REAL4 dy;	/**< The estimated error in the interpolated value */
}
SInterpolateOut;

/** These structures contain the output of the interpolation */
typedef struct
tagDInterpolateOut
{
  REAL8  y;	/**< The interpolated value */
  REAL8 dy;	/**< The estimated error in the interpolated value */
}
DInterpolateOut;

/**
 * These structures contain the interpolation parameters; These are the arrays
 * of \c n domain values \f$x[0] \ldots x[n-1]\f$ and their
 * corresponding values \f$y[0] \ldots y[n-1]\f$
 */
typedef struct
tagSInterpolatePar
{
  UINT4  n;	/**< The number of points in the arrays to use in the interpolation */
  REAL4 *x;	/**< The array of domain values */
  REAL4 *y;	/**< The array of values to interpolate */
}
SInterpolatePar;

/**
 * These structures contain the interpolation parameters; These are the arrays
 * of \c n domain values \f$x[0]\ldots x[n-1]\f$ and their
 * corresponding values \f$y[0]\ldots y[n-1]\f$
 */
typedef struct
tagDInterpolatePar
{
  UINT4  n;	/**< The number of points in the arrays to use in the interpolation */
  REAL8 *x;	/**< The array of domain values */
  REAL8 *y;	/**< The array of values to interpolate */
}
DInterpolatePar;

/* ----- Interpolate.c ----- */
/** \see See \ref Interpolate_h for documentation */
void
LALSPolynomialInterpolation (
    LALStatus          *status,
    SInterpolateOut *output,
    REAL4            target,
    SInterpolatePar *params
    );

/** \see See \ref Interpolate_h for documentation */
void
LALDPolynomialInterpolation (
    LALStatus          *status,
    DInterpolateOut *output,
    REAL8            target,
    DInterpolatePar *params
    );

/** \see See \ref Interpolate_h for documentation */
REAL8
XLALREAL8PolynomialInterpolation (
    REAL8 *yout,
    REAL8  xtarget,
    REAL8 *y,
    REAL8 *x,
    UINT4  n
    );

int
XLALREAL8Interpolation (
    REAL8Sequence *x_in, // Assumed to be in ascending order
    REAL8Sequence *y_in,
    REAL8Sequence *x_out,
    REAL8Sequence *y_out, // Modified in place
    UINT4 n_data_points,
    const gsl_interp_type *itrp_type // Can be NULL -- default it to cubic spline
    );

int
XLALREAL8TimeSeriesInterpolation (
    REAL8TimeSeries *ts_in,
    REAL8Sequence *y_in,
    REAL8Sequence *t_in,
    REAL8Sequence *t_out, // If NULL, interpolate from dt and epoch of ts_in
    UINT4 n_data_points,
    const gsl_interp_type *itrp_type
    );

/*@}*/

#ifdef __cplusplus
}
#endif

#endif /* _INTERPOLATE_H */
