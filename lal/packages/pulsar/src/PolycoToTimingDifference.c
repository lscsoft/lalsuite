/********************* <lalVerbatim file="PolycoToTimingDifferenceCV">
Author: Creighton, T. D.
Revision: $Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{PolycoToTimingDifference.c}}
\label{ss:PolycoToTimingDifference.c}

Computes values of the timing difference $(\tau-t)/\Delta t$ from a
polynomial fit.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{PolycoToTimingDifferenceCP}
\idx{LALPolycoToTimingDifference()}

\subsubsection*{Description}

This function fills a time series \verb@*difference@ with the values
of the normalized timing difference $(\tau-t)/\Delta t$ between the
detector time $t$ and some canonical time $\tau(t)$, where $\Delta t$
is the sampling interval in $t$.  The timing difference function is
computed from the piecewise-polynomial fit stored in \verb@*polyco@,
via Eq.~(\ref{eq:delta-tau}).

\subsubsection*{Algorithm}

By storing the timing difference as a dimensionless quantity, it is
relatively easy to determine rules for resampling the datastream at
equal intervals in $\tau$, since it gives the number of \emph{samples}
difference between the two time coordinates.  When resampling a time
series in $t$, simply track the value of $(\tau-t)/\Delta t$: When
this value increases by +1, choose the next sample after the one that
would otherwise have been chosen; when the value decreases by $-1$,
choose (or repeat) the sample immediately preceding the one that would
otherwise have been chosen.

However, this is not a particularly \emph{efficient} routine for
computing the resampling method, as it requires several floating-point
opertations \emph{per sample}, which is an unacceptable computational
burden for any optimized pulsar search.  It is primarily used to
visualize and check the pulsar phase modulation model.  See the
routine in \verb@CreateResampleRules.c@ for a more efficient
algorithm.

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{PolycoToTimingDifferenceCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include "LALStdlib.h"
#include "LALConstants.h"
#include "AVFactories.h"
#include "Resample.h"

NRCSID(POLYCOTOTIMINGDIFFERENCEC,"$Id$");

/* <lalVerbatim file="PolycoToTimingDifferenceCP"> */
void
LALPolycoToTimingDifference( LALStatus       *stat,
			     REAL4TimeSeries *difference,
			     PolycoStruc     *polyco )
{ /* </lalVerbatim> */
  UINT4 n;          /* Counter over length of time series */
  UINT4 nPoly;      /* Number of polynomial coefficients per fit */
  REAL8 tDiffStart; /* Start time of data series */
  REAL8 tDiffStop;  /* Stop time of data series */
  REAL8 tPolyStart; /* Start time of polynomial fit */
  REAL8 tPolyStop;  /* Stop time of polynomial fit */
  REAL8 t;          /* Time index relative to a fit time */
  REAL8 tNext;      /* Value at which t moves to the next fit */
  REAL8 dt;         /* Sampling interval */
  REAL4 *data;      /* Pointer to time series data */
  REAL4 *tBound;    /* Pointer to bounds of polynomial fit regions */
  REAL4 *t0;        /* Pointer to polynomial fit times */
  REAL4 *poly;      /* Pointer to polynomial coefficients */

  INITSTATUS(stat,"LALPolycoToTimingDifference",
	     POLYCOTOTIMINGDIFFERENCEC);

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
  tPolyStop=tPolyStart+polyco->tBound->data[polyco->tBound->length-1];
  ASSERT(tDiffStop<tPolyStop,stat,RESAMPLEH_ETIME,
	 RESAMPLEH_MSGETIME);
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
