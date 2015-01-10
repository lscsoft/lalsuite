/*
*  Copyright (C) 2007 Jolien Creighton
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
 * \brief Computes values of the timing difference \f$(\tau-t)/\Delta t\f$ from a set of resampling rules.
 *
 * This function fills a time series <tt>*difference</tt> with the values
 * of the normalized timing difference \f$(\tau-t)/\Delta t\f$ between the
 * detector time \f$t\f$ and some canonical time \f$\tau(t)\f$, sampled as a
 * function of \f$t\f$.  This is computed using the resampling rules
 * <tt>*rules</tt>, which specify how one resamples a function of \f$t\f$ as a
 * function of \f$\tau\f$.  Note that \f$\Delta t\f$ in the formula above is the
 * raw (unresampled) sampling interval in \f$t\f$ used to compute the
 * resampling rules, \e not necessarily the sampling interval of the
 * time series.  Thus a shift of \f$\pm1\f$ in the resampling rules
 * corresponds to a change of \f$\pm1\f$ in the normalized timing difference.
 *
 * ### Algorithm ###
 *
 * This routine is quite simple: it increments a time counter
 * \c tNext by resampling rule intervals until it finds the next
 * correction after the current time \c t, accumulating the shifts to
 * the timing function given by <tt>rules->shift</tt>.  It then increments
 * \c t until it steps past \c tNext, each time filling
 * <tt>difference->data</tt> with the current value of the cumulative
 * timing error \f$(\tau-t)/\Delta t\f$.  These are iterated until the time
 * series is completely filled.  Internally, all times are all converted
 * to units of <tt>rules->deltaT</tt> and are measured from
 * <tt>rules->start</tt>.
 *
 * It is worth noting that a shift in the resampling, as given by
 * <tt>rules->shift</tt>, adjusts the value of \f$\tau-t\f$ by shifting the
 * current value of \f$t\f$.  Thus the difference in \e detector time
 * between the \f$n^\mathrm{th}\f$ and the \f$(n+1)^\mathrm{th}\f$ resampling
 * corrections is (<tt>interval[</tt>\f$n\f$<tt>]</tt>\f$\times\f$\c decimate +
 * <tt>shift[</tt>\f$n\f$<tt>]</tt>)\f$\times\Delta t\f$, while the difference in
 * \e canonical time is simply
 * <tt>interval[</tt>\f$n\f$<tt>]</tt>\f$\times\f$\c decimate\f$\times\Delta t\f$.
 * Since <tt>*difference</tt> is assumed to be sampled in detector time
 * \f$t\f$, we use the first of these formulae.
 *
 */
void
LALRulesToTimingDifference( LALStatus       *stat,
			    REAL4TimeSeries *difference,
			    ResampleRules   *rules )
{
  REAL8 tRuleStart; /* Start time of resampling rules */
  REAL8 tDiffStart; /* Start time of time series */
  REAL8 t;          /* Current normalized time in time series */
  REAL8 dt;         /* Normalized sampling interval in time series */
  REAL8 diff = 0;   /* Normalized difference (tau-t)/deltaT */
  REAL8 diffNext;   /* diff after applying the next shift */
  INT4 i;           /* Counter for points in *difference */
  INT4 d;           /* Decimation factor used in defining *rules */
  INT4 tNext;       /* Normalized time of next correction point */
  INT4 *interval;   /* Pointer to data in rules->interval */
  INT2 *shift;      /* Pointer to data in rules->shift */
  REAL4 *data;      /* Pointer to data in difference->data->data */

  INITSTATUS(stat);

  /* Check that the input fields exist. */
  ASSERT( difference, stat, RESAMPLEH_ENUL, RESAMPLEH_MSGENUL );
  ASSERT( difference->data, stat, RESAMPLEH_ENUL, RESAMPLEH_MSGENUL );
  ASSERT( difference->data->data, stat, RESAMPLEH_ENUL,
	  RESAMPLEH_MSGENUL);
  ASSERT( rules, stat, RESAMPLEH_ENUL, RESAMPLEH_MSGENUL );
  ASSERT( rules->interval, stat, RESAMPLEH_ENUL, RESAMPLEH_MSGENUL );
  ASSERT( rules->shift, stat, RESAMPLEH_ENUL, RESAMPLEH_MSGENUL );

  /* Make sure that sampling intervals is positive. */
  ASSERT( difference->deltaT > 0.0, stat, RESAMPLEH_EDTPOS,
	  RESAMPLEH_MSGEDTPOS);
  ASSERT( rules->deltaT > 0.0, stat, RESAMPLEH_EDTPOS,
	  RESAMPLEH_MSGEDTPOS);
  ASSERT( rules->decimate > 0, stat, RESAMPLEH_EDTPOS,
	  RESAMPLEH_MSGEDTPOS);

  /* Make sure that the rules cover the entire time series. */
  tRuleStart = rules->start.gpsSeconds
    + (1.0e-9)*rules->start.gpsNanoSeconds;
  tDiffStart = difference->epoch.gpsSeconds
    + (1.0e-9)*difference->epoch.gpsNanoSeconds;
#ifndef LAL_NDEBUG
  REAL8 tRuleStop;  /* Stop time of resampling rules */
  tRuleStop = rules->stop.gpsSeconds
    + (1.0e-9)*rules->stop.gpsNanoSeconds;
  REAL8 tDiffStop;  /* Stop time of time series */
  tDiffStop = tDiffStart
    + difference->data->length*difference->deltaT;
#endif
  ASSERT( tDiffStop < tRuleStop, stat, RESAMPLEH_ETIME,
	  RESAMPLEH_MSGETIME );
  ASSERT( tDiffStart > tRuleStart, stat, RESAMPLEH_ETIME,
	  RESAMPLEH_MSGETIME );

  /* Find start time, sampling interval, and initial timing
     difference, normalized to the raw sampling interval in the
     resampling rules. */
  t = ( tDiffStart - tRuleStart )/rules->deltaT;
  dt = difference->deltaT/rules->deltaT;
  diffNext = rules->startDiff/rules->deltaT;

  /* Set some local computational parameters. */
  interval = rules->interval;
  shift = rules->shift;
  d = rules->decimate;

  /* Adjust the initial timing difference so that the difference is
     matched to the theoretical function at the *middle* of the first
     correction interval, rather than the start. */
  diffNext += *shift*0.5;

  /* Start filling the time series. */
  i = difference->data->length;
  data = difference->data->data;
  tNext = 0;
  while ( i ) {

    /* First, find the next correction point after the current point
       in time. */
    while ( tNext <= t ) {
      diff = diffNext;
      diffNext += *shift;
      tNext += *(shift++) + *(interval++)*d;
    }

    /* Then, fill the time series up to the next correction point.  If
       the counter i reaches zero, the series is full; break out of
       both loops. */
    while ( t < tNext ) {
      *(data++) = diff;
      t += dt;
      if ( !( --i ) )
	break;
    }
  }

  /* That's it. */
  RETURN( stat );
}
