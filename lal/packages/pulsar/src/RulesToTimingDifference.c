/************************** <lalVerbatim file="RulesToTimingDifferenceCV">
Author: Creighton, T. D.
Revision: $Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{RulesToTimingDifference.c}}
\label{ss:RulesToTimingDifference.c}

Computes values of the timing difference $(\tau-t)/\Delta t$ from a
set of resampling rules.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{RulesToTimingDifferenceCP}
\idx{LALRulesToTimingDifference()}

\subsubsection*{Description}

These function fills a time series \verb@*difference@ with the values
of the normalized timing difference $(\tau-t)/\Delta t$ between the
detector time $t$ and some canonical time $\tau(t)$, sampled as a
function of $t$.  This is computed using the resampling rules
\verb@*rules@, which specify how one resamples a function of $t$ as a
function of $\tau$.  Note that $\Delta t$ in the formula above is the
raw (unresampled) sampling interval in $t$ used to compute the
resampling rules, \emph{not} necessarily the sampling interval of the
time series.  Thus a shift of $\pm1$ in the resampling rules
corresponds to a change of $\pm1$ in the normalized timing difference.

\subsubsection*{Algorithm}

This routine is quite simple: it increments a time counter
\verb@tNext@ by resampling rule intervals until it finds the next
correction after the current time \verb@t@, accumulating the shifts to
the timing function given by \verb@rules->shift@.  It then increments
\verb@t@ until it steps past \verb@tNext@, each time filling
\verb@difference->data@ with the current value of the cumulative
timing error $(\tau-t)/\Delta t$.  These are iterated until the time
series is completely filled.  Internally, all times are all converted
to units of \verb@rules->deltaT@ and are measured from
\verb@rules->start@.

It is worth noting that a shift in the resampling, as given by
\verb@rules->shift@, adjusts the value of $\tau-t$ by shifting the
current value of $t$.  Thus the difference in \emph{detector} time
between the $n^\mathrm{th}$ and the $(n+1)^\mathrm{th}$ resampling
corrections is (\verb@interval[@$n$\verb@]@$\times$\verb@decimate@ +
\verb@shift[@$n$\verb@]@)$\times\Delta t$, while the difference in
\emph{canonical} time is simply
\verb@interval[@$n$\verb@]@$\times$\verb@decimate@$\times\Delta t$.
Since \verb@*difference@ is assumed to be sampled in detector time
$t$, we use the first of these formulae.

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{RulesToTimingDifferenceCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include "LALStdlib.h"
#include "LALConstants.h"
#include "AVFactories.h"
#include "Resample.h"

NRCSID(RULESTOTIMINGDIFFERENCEC,"$Id$");

/* <lalVerbatim file="RulesToTimingDifferenceCP"> */
void
LALRulesToTimingDifference( LALStatus       *stat,
			    REAL4TimeSeries *difference,
			    ResampleRules   *rules )
{ /* </lalVerbatim> */
  REAL8 tRuleStart; /* Start time of resampling rules */
  REAL8 tRuleStop;  /* Stop time of resampling rules */
  REAL8 tDiffStart; /* Start time of time series */
  REAL8 tDiffStop;  /* Stop time of time series */
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

  INITSTATUS( stat, "LALRulesToTimingDifference",
	      RULESTOTIMINGDIFFERENCEC);

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
  tRuleStop = rules->stop.gpsSeconds
    + (1.0e-9)*rules->stop.gpsNanoSeconds;
  tDiffStart = difference->epoch.gpsSeconds
    + (1.0e-9)*difference->epoch.gpsNanoSeconds;
  tDiffStop = tDiffStart
    + difference->data->length*difference->deltaT;
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
