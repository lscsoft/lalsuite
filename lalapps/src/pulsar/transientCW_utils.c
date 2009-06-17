/*
 * Copyright (C) 2009 Reinhard Prix
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

/*********************************************************************************/
/** \author R. Prix
 * \file
 * \brief
 * Some helper functions useful for "transient CWs", mostly applying transient window
 * functions.
 *
 *********************************************************************************/
#include "config.h"

#include "transientCW_utils.h"

/* System includes */
#include <math.h>

/* LAL-includes */
#include <lal/XLALError.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/LogPrintf.h>


/* ---------- internal prototypes ---------- */
REAL4Vector *XLALGetTransientWindowVals ( const LIGOTimeGPSVector *tGPS, const transientWindow_t *TransientWindowParams );


/* ==================== function definitions ==================== */

/** apply a "transient CW window" described by TransientWindowParams to the given
 * timeseries
 */
int
XLALApplyTransientWindow ( REAL4TimeSeries *series, transientWindow_t TransientWindowParams )
{
  const CHAR *fn = "XLALApplyTransientWindow()";
  UINT4 i;

  REAL8 ts_t0, ts_dt, ts_T;	/* start-time, stepsize and duration of input timeseries */
  REAL8 ti;
  INT4 i0, i1;			/* time-series index corresonding to start-time (and end-time) of transient window */

  if ( !series || !series->data ){
    XLALPrintError ("%s: Illegal NULL in input timeseries!\n", fn );
    return XLAL_EINVAL;
  }
  ts_t0 = XLALGPSGetREAL8 ( &series->epoch );
  ts_dt = series->deltaT;
  ts_T  = ts_dt * series->data->length;

  i0 = ( TransientWindowParams.t0 - ts_t0 ) / ts_dt;
  if ( i0 < 0 ) i0 = 0;

  switch ( TransientWindowParams.type )
    {
    case TRANSIENT_NONE:
      return XLAL_SUCCESS;	/* nothing to be done here */
      break;

    case TRANSIENT_RECTANGULAR:	/* standard 'rectangular window */
      for ( i = 0; i < (UINT4)i0; i ++ ) {
	series->data->data[i] = 0;
      }
      i1 = (TransientWindowParams.t0 + TransientWindowParams.tau - ts_t0 ) / ts_dt + 1;
      if ( i1 < 0 ) i1 = 0;
      if ( (UINT4)i1 >= series->data->length ) i1 = series->data->length - 1;

      for ( i = i1; i < series->data->length; i ++) {
	series->data->data[i] = 0;
      }
      break;

    case TRANSIENT_EXPONENTIAL:
      for ( i = 0; i < (UINT4)i0; i ++ ) {
	series->data->data[i] = 0;
      }
      ti = 0;
      for ( i=i0; i < series->data->length; i ++)
	{
	  REAL8 fact = exp( - ti / TransientWindowParams.tau );
	  ti += ts_dt;
	  series->data->data[i] *= fact;
	}
      break;

    default:
      XLALPrintError("Illegal transient-signal window type specified '%d'\n", TransientWindowParams.type );
      return XLAL_EINVAL;
      break;
    }

  return XLAL_SUCCESS;

} /* XLALApplyTransientWindow() */


/** Compute values of given transient window function at GPS times 'tGPS'.
 */
REAL4Vector *
XLALGetTransientWindowVals ( const LIGOTimeGPSVector *tGPS, const transientWindow_t *TransientWindowParams )
{
  static const char *fn = "XLALGetTransientWindowVals()";
  UINT4 i, numTS;
  REAL4Vector *ret;
  REAL4 t0, t1, tau;
  REAL8 ti;

  if ( !tGPS || tGPS->length == 0 || !TransientWindowParams ) {
    XLALPrintError ("%s: invalid NULL input.\n", fn );
    XLAL_ERROR_NULL( fn, XLAL_EINVAL );
  }
  numTS = tGPS->length;

  if ( ( ret = XLALCreateREAL4Vector ( numTS )) == NULL ) {
    XLALPrintError ( "%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numTS );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  t0 = TransientWindowParams->t0;
  tau = TransientWindowParams->tau;
  t1 = t0 + tau;

  switch ( TransientWindowParams->type )
    {
    case TRANSIENT_NONE:
      for ( i = 0; i < numTS; i ++ )
        ret->data[i] = 1;
      break;

    case TRANSIENT_RECTANGULAR:		/* standard 'rectangular window */
      for ( i = 0; i < numTS; i ++ )
        {
          ti = XLALGPSGetREAL8 ( &tGPS->data[i] );
          if ( ( ti >= t0 ) && (ti <= t1) )
            ret->data[i] = 1;
          else
            ret->data[i] = 0;
        } /* for i < numTS */
      break;

    case TRANSIENT_EXPONENTIAL:
      for ( i=0; i < numTS; i ++ )
        {
          ti = XLALGPSGetREAL8 ( &tGPS->data[i] );
          if ( ti >= t0 )
            ret->data[i] = exp( (t0 - ti) / tau );
          else
            ret->data[i] = 0;
        } /* for i < numTS */
      break;

    default:
      XLALPrintError ("%s: invalid transient window type %d not in [%d, %d].\n",
                      fn, TransientWindowParams->type, TRANSIENT_NONE, TRANSIENT_LAST -1 );
      XLALDestroyREAL4Vector ( ret );
      XLAL_ERROR_NULL( fn, XLAL_EINVAL );

      break;

    } /* switch transient type */

  return ret;

} /* XLALGetTransientWindowVals() */


/** apply transient window to give multi noise-weights, associated with given
 * multi timestamps
 */
int
XLALApplyTransientWindow2NoiseWeights ( MultiNoiseWeights *multiNoiseWeights,	/**< [in/out] noise weights to apply transient window to */
                                        const MultiLIGOTimeGPSVector *multiTS,	/**< [in] associated timestamps of noise-weights */
                                        transientWindow_t TransientWindowParams	/**< [in] transient window parameters */
                                        )
{
  static const char *fn = "XLALApplyTransientWindow2NoiseWeights()";

  UINT4 numIFOs, X;
  UINT4 numTS, i;
  REAL4Vector *win;

  /* check input consistency */
  if ( !multiNoiseWeights || multiNoiseWeights->length == 0 ) {
    XLALPrintError ("%s: empty or NULL input 'multiNoiseWeights'.\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  if ( !multiTS || multiTS->length == 0 ) {
    XLALPrintError ("%s: empty or NULL input 'multiTS'.\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  numIFOs = multiNoiseWeights->length;
  if ( multiTS->length != numIFOs ) {
    XLALPrintError ("%s: inconsistent numIFOs between 'multiNoiseWeights' (%d) and 'multiTS' (%d).\n", fn, numIFOs, multiTS->length );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  for ( X = 0; X < numIFOs; X ++ )
    {
      numTS = multiNoiseWeights->data[X]->length;

      if ( multiTS->data[X]->length != numTS ) {
        XLALPrintError ("%s: inconsistent number of timesteps 'multiNoiseWeights[%d]' (%d) and 'multiTS[%d]' (%d).\n", fn, X, numTS, X, multiTS->data[X]->length );
        XLAL_ERROR ( fn, XLAL_EINVAL );
      }

      if ( ( win = XLALGetTransientWindowVals ( multiTS->data[X], &TransientWindowParams )) == NULL ) {
        XLALPrintError ("%s: XLALGetTransientWindowVals() failed. xlalErrno = %d.\n", fn, xlalErrno );
        XLAL_ERROR ( fn, XLAL_EFUNC );
      }

      for ( i=0; i < numTS; i ++ )
        {
          multiNoiseWeights->data[X]->data[i] *= win->data[i];
        } /* for i < numTS */

      XLALDestroyREAL4Vector ( win );

    } /* for X < numIFOs */

  return XLAL_SUCCESS;

} /* XLALApplyTransientWindow2NoiseWeights() */
