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
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/Resample.h>

/**
 * \author Creighton, T. D.
 * \ingroup Resample_h
 * \brief Computes values of the timing difference \f$(\tau-t)/\Delta t\f$ from a polynomial fit.
 *
 * This function fills a time series <tt>*difference</tt> with the values
 * of the normalized timing difference \f$(\tau-t)/\Delta t\f$ between the
 * detector time \f$t\f$ and some canonical time \f$\tau(t)\f$, where \f$\Delta t\f$
 * is the sampling interval in \f$t\f$.  The timing difference function is
 * computed from the piecewise-polynomial fit stored in <tt>*polyco</tt>,
 * via \eqref{eq_delta-tau}.
 *
 * ### Algorithm ###
 *
 * By storing the timing difference as a dimensionless quantity, it is
 * relatively easy to determine rules for resampling the datastream at
 * equal intervals in \f$\tau\f$, since it gives the number of \e samples
 * difference between the two time coordinates.  When resampling a time
 * series in \f$t\f$, simply track the value of \f$(\tau-t)/\Delta t\f$: When
 * this value increases by +1, choose the next sample after the one that
 * would otherwise have been chosen; when the value decreases by \f$-1\f$,
 * choose (or repeat) the sample immediately preceding the one that would
 * otherwise have been chosen.
 *
 * However, this is not a particularly \e efficient routine for
 * computing the resampling method, as it requires several floating-point
 * opertations <em>per sample</em>, which is an unacceptable computational
 * burden for any optimized pulsar search.  It is primarily used to
 * visualize and check the pulsar phase modulation model.  See the
 * routine in CreateResampleRules.c for a more efficient algorithm.
 */
void
LALPolycoToTimingDifference( LALStatus       *stat,
			     REAL4TimeSeries *difference,
			     PolycoStruc     *polyco )
{
  UINT4 n;          /* Counter over length of time series */
  UINT4 nPoly;      /* Number of polynomial coefficients per fit */
  REAL8 tDiffStart; /* Start time of data series */
  REAL8 tDiffStop;  /* Stop time of data series */
  REAL8 tPolyStart; /* Start time of polynomial fit */
  REAL8 t;          /* Time index relative to a fit time */
  REAL8 tNext;      /* Value at which t moves to the next fit */
  REAL8 dt;         /* Sampling interval */
  REAL4 *data;      /* Pointer to time series data */
  REAL4 *tBound;    /* Pointer to bounds of polynomial fit regions */
  REAL4 *t0;        /* Pointer to polynomial fit times */
  REAL4 *poly;      /* Pointer to polynomial coefficients */

  INITSTATUS(stat);

  /* Check that the input fields exist. */
  ASSERT(difference,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(difference->data,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(difference->data->data,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(polyco,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(polyco->t0,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(polyco->t0->data,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(polyco->polyco,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);
  ASSERT(polyco->polyco->data,stat,RESAMPLEH_ENUL,RESAMPLEH_MSGENUL);

  /* Make sure that sampling interval is positive. */
  ASSERT(difference->deltaT>0.0,stat,RESAMPLEH_EDTPOS,
	 RESAMPLEH_MSGEDTPOS);

  /* This should be obvious, but make sure that there are polynomial
     coefficients and fitting times defined for each ploynomial
     fitting region. */
  ASSERT(polyco->tBound->length==polyco->t0->length,stat,
	 RESAMPLEH_ELENGTH,RESAMPLEH_MSGELENGTH);
  ASSERT(polyco->tBound->length==polyco->polyco->length,stat,
	 RESAMPLEH_ELENGTH,RESAMPLEH_MSGELENGTH);

  /* Check that the desired time series is a subset of the timespan
     for which we have polynomial coefficients. */
  tDiffStart=difference->epoch.gpsSeconds
    +(1.0e-9)*difference->epoch.gpsNanoSeconds;
  tDiffStop=tDiffStart+difference->data->length*difference->deltaT;
  tPolyStart=polyco->start.gpsSeconds
    +(1.0e-9)*polyco->start.gpsNanoSeconds;
#ifndef LAL_NDEBUG
  REAL8 tPolyStop;  /* Stop time of polynomial fit */
  tPolyStop=tPolyStart+polyco->tBound->data[polyco->tBound->length-1];
  ASSERT(tDiffStop<tPolyStop,stat,RESAMPLEH_ETIME,
	 RESAMPLEH_MSGETIME);
#endif
  ASSERT(tDiffStart>tPolyStart,stat,RESAMPLEH_ETIME,
	 RESAMPLEH_MSGETIME);

  /* Assign some computational variables. */
  data=difference->data->data;
  nPoly=polyco->polyco->vectorLength;
  tBound=polyco->tBound->data;
  t0=polyco->t0->data;
  poly=polyco->polyco->data;
  dt=difference->deltaT;

  /* Find start and stop times relative to the polyco reference
     time. */
  tDiffStart-=tPolyStart;
  tDiffStop-=tPolyStart;

  /* Determine the fitting region in which the series start time
     lies. */
  for(;tDiffStart>*(tBound++);poly+=nPoly,t0++)
    ;
  tBound--;

  /* The main loop: compute each time series datum. */
  t=tDiffStart-*t0;
  tNext=*tBound-*t0;
  n=difference->data->length;
  while(n--){
    REAL4 tn=1.0;       /* Current power of t */
    REAL4 diff=poly[0]; /* Current cumulative polynomial */
    UINT4 i=0;          /* Counter over polynomial coefficients */

    /* The 0th coefficient is already stored in diff.  This loop adds
       the other coefficients multiplied by their powers of t. */
    while(++i<nPoly)
      diff+=(tn*=t)*poly[i];

    /* Normalize and assign the result. */
    *(data++)=diff/dt;

    /* If we've stepped out of the current fitting region, move to the
       next one.  (It has already been determined that we will never
       step past the final fitting region.) */
    if((t+=dt)>tNext){
      t+=*t0++;
      t-=*t0;
      tNext=*(++tBound)-*t0;
      poly+=nPoly;
    }
  }

  /* That's it!  Pretty painless, huh? */
  RETURN(stat);
}
