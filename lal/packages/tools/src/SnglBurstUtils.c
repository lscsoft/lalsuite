/*----------------------------------------------------------------------- 
 * 
 * File Name: SnglBurstUtils.c
 *
 * Author: Brown, D. A.  Brady, P. R. and Ray Majumder, S. K.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="SnglBurstUtilsCV">
Author: Brown, D. A. Brady, P. R. and Ray Majumder, S. K
$Id$
</lalVerbatim> 
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/Date.h>

NRCSID( SNGLBURSTUTILSC, "$Id$" );

#if 0
<lalLaTeX>
\subsection{Module \texttt{SnglBurstUtils.c}}

\noindent Blah.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{SnglBurstUtilsCP}
\idx{LALSortSnglBurst()}
\idx{XLALSortSnglBurst()}
\idx{LALCompareSnglBurst()}
\idx{XLALCompareSnglBurstByLowFreq()}
\idx{XLALCompareSnglBurstByStartTime()}
\idx{XLALCompareSnglBurstByStartTimeAndLowFreq()}

\subsubsection*{Description}

\noindent Blah.

\subsubsection*{Algorithm}

\noindent None.

\subsubsection*{Uses}

\noindent None.

\subsubsection*{Notes}
%% Any relevant notes.
 
\vfill{\footnotesize\input{SnglBurstUtilsCV}}

</lalLaTeX>
#endif


/*
 * A few quickies for convenience.
 */

static INT8 start_time(const SnglBurstTable *x)
{
	return(XLALGPStoINT8(&x->start_time));
}

static INT8 end_time(const SnglBurstTable *x)
{
	return(start_time(x) + 1e9 * x->duration);
}

static REAL4 lo_freq(const SnglBurstTable *x)
{
	return(x->central_freq - x->bandwidth / 2);
}

static REAL4 hi_freq(const SnglBurstTable *x)
{
	return(x->central_freq + x->bandwidth / 2);
}


/*
 * cluster events a and b, storing result in a
 */

static void XLALClusterSnglBurst(SnglBurstTable *a, SnglBurstTable *b)
{
	REAL4 f_lo, f_hi;
	INT8 ta_start, ta_end, tb_start, tb_end;

	/* the cluster's frequency band is the smallest band containing the
	 * bands of the two original events */

	f_lo = lo_freq(a);
	f_hi = hi_freq(a);
	if(lo_freq(b) < f_lo)
		f_lo = lo_freq(b);
	if(hi_freq(b) > f_hi)
		f_hi = hi_freq(b);

	a->central_freq = (f_hi + f_lo) / 2;
	a->bandwidth = f_hi - f_lo;

	/* the cluster's time interval is the smallest interval containing the
	 * intervals of the two original events */

	ta_start = start_time(a);
	tb_start = start_time(b);
	ta_end = end_time(a);
	tb_end = end_time(b);
	if(tb_start < ta_start) {
		a->start_time = b->start_time;
		ta_start = tb_start;
	}
	if(tb_end > ta_end)
		ta_end = tb_end;
	a->duration = (ta_end - ta_start) / 1e9;

	/* the amplitude, SNR, confidence, and peak time of the cluster are
	 * those of the loudest of the two events */

	if(a->amplitude < b->amplitude) {
		a->amplitude = b->amplitude;
		a->snr = b->snr;
		a->confidence = b->confidence;
		a->peak_time = b->peak_time;
	}
}

static int ModifiedforClustering(SnglBurstTable *prevEvent, SnglBurstTable *thisEvent)
{
	REAL4 fa1, fa2, fb1, fb2;
	REAL8 deltaT;
	const REAL8 epsilon = 1e-8;	/* seconds */

	deltaT = XLALDeltaFloatGPS(&prevEvent->peak_time, &thisEvent->peak_time);

	fa1 = lo_freq(prevEvent);
	fa2 = hi_freq(prevEvent);

	fb1 = lo_freq(thisEvent);
	fb2 = hi_freq(thisEvent);

	if((fabs(deltaT) < epsilon) &&
	   ( (fb1 >= fa1 && fb1 <= fa2) ||
	     (fb2 >= fa1 && fb2 <= fa2) ||
	     (fa1 >= fb1 && fa1 <= fb2) ||
	     (fa2 >= fb1 && fa2 <= fb2) )) {
		XLALClusterSnglBurst(prevEvent, thisEvent);
		return 1;
	}
	return 0;
}


/* <lalVerbatim file="SnglBurstUtilsCP"> */
void
XLALSortSnglBurst(
	SnglBurstTable **head,
	int (*comparefunc)(const SnglBurstTable **, const SnglBurstTable **)
)
/* </lalVerbatim> */
{
	INT4 i;
	INT4 length;
	SnglBurstTable *event;
	SnglBurstTable **array;

	/* empty list --> no-op */
	if(!head || !*head)
		return;

	/* count the number of events in the list */
	for(length = 0, event = *head; event; event = event->next)
		length++;

	/* construct an array of pointers into the list */
	array = LALCalloc(length, sizeof(*array));
	for(i = 0, event = *head; event; event = event->next)
		array[i++] = event;

	/* sort the array using the specified function */
	qsort(array, length, sizeof(*array), (int(*)(const void *, const void *)) comparefunc);

	/* re-link the list according to the sorted array */
	for(i = 0; i < length; i++, head = &(*head)->next)
		*head = array[i];
	*head = NULL;

	/* free the array */
	LALFree(array);
}


/* <lalVerbatim file="SnglBurstUtilsCP"> */
void
LALSortSnglBurst(
	LALStatus *status,
	SnglBurstTable **head,
	int (*comparefunc)(const SnglBurstTable **, const SnglBurstTable **)
)
/* </lalVerbatim> */
{
	INITSTATUS(status, "LALSortSnglBurst", SNGLBURSTUTILSC);
	XLALSortSnglBurst(head, comparefunc);
	RETURN(status);
}


/* <lalVerbatim file="SnglBurstUtilsCP"> */
int
XLALCompareSnglBurstByStartTime(
	const SnglBurstTable **a,
	const SnglBurstTable **b
)
/* </lalVerbatim> */
{
	INT8 ta, tb;

	ta = start_time(*a);
	tb = start_time(*b);

	if(ta > tb)
		return(1);
	else if(ta < tb)
		return(-1);
	return(0);
}


/* <lalVerbatim file="SnglBurstUtilsCP"> */
int
XLALCompareSnglBurstByLowFreq(
	const SnglBurstTable **a,
	const SnglBurstTable **b
)
/* </lalVerbatim> */
{
	REAL4 flowa, flowb;

	flowa = lo_freq(*a);
	flowb = lo_freq(*b);

	if(flowa > flowb)
		return(1);
	else if(flowa < flowb)
		return(-1);
	return(0);
}


/* <lalVerbatim file="SnglBurstUtilsCP"> */
int
XLALCompareSnglBurstByStartTimeAndLowFreq(
	const SnglBurstTable **a,
	const SnglBurstTable **b
)
/* </lalVerbatim> */
{
	int result;

	result = XLALCompareSnglBurstByStartTime(a, b);
	if(result)
		return(result);
	return(XLALCompareSnglBurstByLowFreq(a, b));
}


/* <lalVerbatim file="SnglBurstUtilsCP"> */
void
LALCompareSnglBurst(
    LALStatus             *status,
    SnglBurstTable        *aPtr,
    SnglBurstTable        *bPtr,
    SnglBurstAccuracy     *params
    )
/* </lalVerbatim> */
{
  INT8 ta1, ta2, tb1, tb2;
  REAL4 fa1, fa2, fb1, fb2;
  REAL4 dm1, dm2;
  REAL4 sigmaRatio;

  INITSTATUS (status, "LALCompareSnglBurst", SNGLBURSTUTILSC);
  ATTATCHSTATUSPTR (status);

  params->match=0;

  LALGPStoINT8( status->statusPtr, &ta1, &(aPtr->start_time) );
  LALGPStoINT8( status->statusPtr, &tb1, &(bPtr->start_time) );

  ta2 = ta1 + ( 1e9 * aPtr->duration );
  tb2 = tb1 + ( 1e9 * bPtr->duration );
  fa1 = (aPtr->central_freq) - 0.5*(aPtr->bandwidth);
  fa2 = (aPtr->central_freq) + 0.5*(aPtr->bandwidth);
  fb1 = (bPtr->central_freq) - 0.5*(bPtr->bandwidth);
  fb2 = (bPtr->central_freq) + 0.5*(bPtr->bandwidth);

   if(((tb1 >= ta1 && tb1 <= ta2) || (tb2 >= ta1 && tb2 <= ta2)) || ((ta1 >= tb1 && ta1 <= tb2) || (ta2 >= tb1 && ta2 <= tb2)))
     {
       if((fb1 >= fa1 && fb1 <= fa2) || (fb2 >= fa1 && fb2 <= fa2) || (fa1 >= fb1 && fa1 <= fb2) || (fa2 >= fb1 && fa2 <= fb2))
	 {
	   params->match = 1;
	 }
     }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/* <lalVerbatim file="SnglBurstUtilsCP"> */
void
LALCompareSimBurstAndSnglBurst(
    LALStatus             *status,
    SimBurstTable         *aPtr,
    SnglBurstTable        *bPtr,
    SnglBurstAccuracy     *params
    )
/* </lalVerbatim> */
{
  INT8 ta, tb1, tb2;
  REAL4 fa, fb1, fb2;

  INITSTATUS(status, "LALCompareSimBurstAndSnglBurst", SNGLBURSTUTILSC);
  ATTATCHSTATUSPTR(status);

  LALGPStoINT8(status->statusPtr, &ta, &aPtr->geocent_peak_time);
  LALGPStoINT8(status->statusPtr, &tb1, &bPtr->start_time);

  fa = aPtr->freq;
  tb2 = tb1 + 1e9 * bPtr->duration;
  fb1 = bPtr->central_freq - 0.5 * bPtr->bandwidth;
  fb2 = bPtr->central_freq + 0.5 * bPtr->bandwidth;

  params->match = (tb1 < ta) && (ta < tb2) && (fb1 < fa) && (fa < fb2);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/* <lalVerbatim file="SnglBurstUtilsCP"> */
void
LALClusterSnglBurstTable (
	LALStatus        *status,
	SnglBurstTable   *burstEvent,
	INT4             *nevents
)
/* </lalVerbatim> */
{
	SnglBurstTable *thisEvent = NULL;
	SnglBurstTable *prevEvent = NULL;
	SnglBurstTable *startEvent = NULL;
	INT4 i, j;
	INT4 numModEvent = 1;

	INITSTATUS (status, "LALClusterSnglBurstTable", SNGLBURSTUTILSC);
	ATTATCHSTATUSPTR (status);

	startEvent = burstEvent;

	for(thisEvent = burstEvent->next; thisEvent; thisEvent = prevEvent->next) {
		prevEvent = startEvent;

		for(i = numModEvent; i > 0; i--) {
			if(ModifiedforClustering(prevEvent,thisEvent)) {
				for(j = i; j > 1; j--)
					prevEvent = prevEvent->next;
				prevEvent->next = thisEvent->next;
				LALFree(thisEvent);
				break;
			}
			else
				prevEvent = prevEvent->next;
		}

		if(i == 0)
			numModEvent++;
	}

	/* count the number of events in the modified list */
	for(*nevents = 1; startEvent; startEvent = startEvent->next)
		*nevents++;

	/* normal exit */
	DETATCHSTATUSPTR (status);
	RETURN (status);
}
