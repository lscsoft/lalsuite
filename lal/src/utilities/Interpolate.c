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

/* ---------- see Interpolate.h for doxygen documentation ---------- */

#include <math.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/Interpolate.h>


void
LALSPolynomialInterpolation (
    LALStatus       *status,
    SInterpolateOut *output,
    REAL4            target,
    SInterpolatePar *params
    )
{
  REAL4 *dn;   /* difference in a step down */
  REAL4 *up;   /* difference in a step up   */
  REAL4  diff;
  UINT4  near;
  UINT4  order;
  UINT4  n;
  UINT4  i;

  INITSTATUS(status);

  ASSERT (output, status, INTERPOLATEH_ENULL, INTERPOLATEH_MSGENULL);
  ASSERT (params, status, INTERPOLATEH_ENULL, INTERPOLATEH_MSGENULL);
  ASSERT (params->x, status, INTERPOLATEH_ENULL, INTERPOLATEH_MSGENULL);
  ASSERT (params->y, status, INTERPOLATEH_ENULL, INTERPOLATEH_MSGENULL);

  n = params->n;
  ASSERT (n > 1, status, INTERPOLATEH_ESIZE, INTERPOLATEH_MSGESIZE);

  dn = (REAL4 *) LALMalloc (n*sizeof(REAL4));
  ASSERT (dn, status, INTERPOLATEH_ENULL, INTERPOLATEH_MSGENULL);

  up = (REAL4 *) LALMalloc (n*sizeof(REAL4));
  ASSERT (up, status, INTERPOLATEH_ENULL, INTERPOLATEH_MSGENULL);


  /*
   *
   * Initialize dn[] and up[] and find element of domain nearest the target.
   *
   */


  memcpy (dn, params->y, n*sizeof(REAL4));
  memcpy (up, params->y, n*sizeof(REAL4));
  diff = target - params->x[near = 0];
  for (i = 1; diff > 0 && i < n; ++i)
  {
    REAL4 tmp = target - params->x[i];
    diff = (fabs(tmp) < fabs(diff) ? near = i, tmp : diff);
  }
  output->y = params->y[near];


  /*
   *
   * Recompute dn[] and up[] for each polynomial order.
   *
   */


  for (order = 1; order < n; ++order)
  {
    UINT4 imax = n - order;
    for (i = 0; i < imax; ++i)
    {
      REAL4 xdn = params->x[i];
      REAL4 xup = params->x[i + order];
      REAL4 den = xdn - xup;
      REAL4 fac;
      ASSERT (den != 0, status, INTERPOLATEH_EZERO, INTERPOLATEH_MSGEZERO);
      fac   = (dn[i + 1] - up[i])/den;
      dn[i] = fac*(xdn - target);
      up[i] = fac*(xup - target);
    }

    /* go down unless impossible */
    output->y += (near < imax ? dn[near] : up[--near]);
  }

  output->dy = fabs(dn[0]) < fabs(up[0]) ? fabs(dn[0]) : fabs(up[0]);

  LALFree (dn);
  LALFree (up);
  RETURN  (status);
}



void
LALDPolynomialInterpolation (
    LALStatus       *status,
    DInterpolateOut *output,
    REAL8            target,
    DInterpolatePar *params
    )
{
  INITSTATUS(status);

  XLAL_PRINT_DEPRECATION_WARNING("XLALDPolynomialInterpolation");

  ASSERT (output, status, INTERPOLATEH_ENULL, INTERPOLATEH_MSGENULL);
  ASSERT (params, status, INTERPOLATEH_ENULL, INTERPOLATEH_MSGENULL);
  ASSERT (params->x, status, INTERPOLATEH_ENULL, INTERPOLATEH_MSGENULL);
  ASSERT (params->y, status, INTERPOLATEH_ENULL, INTERPOLATEH_MSGENULL);

  output->dy = XLALREAL8PolynomialInterpolation(&(output->y), target, params->y, params->x, params->n);
  if (xlalErrno)
    ABORTXLAL (status);

  RETURN  (status);
}


REAL8
XLALREAL8PolynomialInterpolation (
    REAL8 *yout,
    REAL8  xtarget,
    REAL8 *y,
    REAL8 *x,
    UINT4  n
    )
{
  REAL8  dy;
  REAL8 *dn;   /* difference in a step down */
  REAL8 *up;   /* difference in a step up   */
  REAL8  diff;
  UINT4  near = 0;
  UINT4  order;
  UINT4  i;

  if ( yout == NULL || y == NULL || x == NULL )
    XLAL_ERROR_REAL8(XLAL_EFAULT);

  if ( n <= 1 )
    XLAL_ERROR_REAL8(XLAL_ESIZE);

  dn = (REAL8 *) LALMalloc (n*sizeof(*dn));
  if ( !dn )
    XLAL_ERROR_REAL8(XLAL_ENOMEM);

  up = (REAL8 *) LALMalloc (n*sizeof(*up));
  if ( !up )
    XLAL_ERROR_REAL8(XLAL_ENOMEM);


  /*
   *
   * Initialize dn[] and up[] and find element of domain nearest xtarget.
   *
   */


  memcpy (dn, y, n*sizeof(*dn));
  memcpy (up, y, n*sizeof(*up));
  diff = xtarget - x[near];
  for (i = 1; diff > 0 && i < n; ++i)
  {
    REAL8 tmp = xtarget - x[i];
    diff = (fabs(tmp) < fabs(diff) ? near = i, tmp : diff);
  }
  *yout = y[near];


  /*
   *
   * Recompute dn[] and up[] for each polynomial order.
   *
   */


  for (order = 1; order < n; ++order)
  {
    UINT4 imax = n - order;
    for (i = 0; i < imax; ++i)
    {
      REAL8 xdn = x[i];
      REAL8 xup = x[i + order];
      REAL8 den = xdn - xup;
      REAL8 fac;
      if ( !den )
      {
	LALFree (dn);
	LALFree (up);
        XLAL_ERROR_REAL8 (XLAL_EFPDIV0);
      }
      fac   = (dn[i + 1] - up[i])/den;
      dn[i] = fac*(xdn - xtarget);
      up[i] = fac*(xup - xtarget);
    }

    /* go down unless impossible */
    *yout += (near < imax ? dn[near] : up[--near]);
  }

  dy = fabs(dn[0]) < fabs(up[0]) ? fabs(dn[0]) : fabs(up[0]);

  LALFree (dn);
  LALFree (up);
  return dy;
}

int
XLALREAL8Interpolation (
    REAL8Sequence *x_in, // Assumed to be in ascending order
    REAL8Sequence *y_in,
    REAL8Sequence *x_out,
    REAL8Sequence *y_out, // Modified in place
    UINT4 n_data_points,
    const gsl_interp_type *itrp_type // Can be NULL -- default it to cubic spline
    )
{
	size_t i=0;
	REAL8 x_min, x_max;

	gsl_interp_accel *i_acc = gsl_interp_accel_alloc();

	/*
 	 * What type of interpolation do want? Default is cubic spline
 	 * Other types can be found at
 	 * http://www.gnu.org/software/gsl/manual/html_node/Interpolation-Types.html
 	 */
	gsl_spline *spline;
	if (!itrp_type) {
		spline = gsl_spline_alloc(gsl_interp_cspline, n_data_points);
	} else {
		spline = gsl_spline_alloc(itrp_type, n_data_points);
	}
	gsl_spline_init(spline, x_in->data, y_in->data, n_data_points);

	/*
	 * Do the interpolation.
	 */
	x_min = x_in->data[0];
	x_max = x_in->data[x_in->length-1];
	for(; i<x_out->length; i++) {
		/*
		 * The gsl interpolator is known to do funny and sometimes
		 * inconsistent things outside the boundaries of the input,
		 * so we check for and bail if we don't meet that criteria
		 */
		if (x_out->data[i] < x_min || x_out->data[i] > x_max ) {
			XLALPrintError("XLAL Error - %s: %g is outside the bounds of the interpolant.\n", __func__, x_out->data[i]);
    			XLAL_ERROR_VAL(0, XLAL_EINVAL);
		}
		y_out->data[i] = gsl_spline_eval(spline, x_out->data[i], i_acc);
	}

	/*
	 * Clean up.
	 */
	gsl_spline_free(spline);
	gsl_interp_accel_free(i_acc);
	return 0;
}

int
XLALREAL8TimeSeriesInterpolation (
    REAL8TimeSeries *ts_in,
    REAL8Sequence *y_in,
    REAL8Sequence *t_in,
    REAL8Sequence *t_out, // If NULL, interpolate from dt and epoch of ts_in
    UINT4 n_data_points,
    const gsl_interp_type *itrp_type
    )
{
	/*
 	 * If the user doesn't provide a sampling series to interpolate to,
 	 * we'll infer it from the epoch and sample rate.
 	 */
	int length = ts_in->data->length;
	int destroy_t_out = 0;
	if (!t_out) {
		if (ts_in->deltaT < 0.0) {
			XLALPrintError("XLAL Error - %s: Not enough information in the time series to infer sampling values. Supply proper epoch and deltaT.", __func__);
    			XLAL_ERROR_VAL(0, XLAL_EINVAL);
		}
		REAL8 start = XLALGPSGetREAL8(&ts_in->epoch);
		destroy_t_out = 1;
		t_out = XLALCreateREAL8Sequence(length);
		length--;
		for (; length>=0; length--){
			t_out->data[length] = length*ts_in->deltaT + start;
		}
	}
	XLALREAL8Interpolation(t_in, y_in, t_out, ts_in->data, n_data_points, itrp_type);
	/*
 	 * We created it, so we'll clean it up too.
 	 */
	if (destroy_t_out) {
		XLALDestroyREAL8Sequence(t_out);
	}
	return 0;
}
