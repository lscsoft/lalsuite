/*----------------------------------------------------------------------- 
 * 
 * File Name: SnglBurstUtils.c
 *
 * Author: Brown, D. A.  Brady, P. R. Ray Majumder, S. K. and Cannon, K. C.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="SnglBurstUtilsCV">
Author: Brown, D. A. Brady, P. R. Ray Majumder, S. K and Cannon, K. C.
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
\idx{XLALCompareSnglBurstByFreq()}
\idx{XLALCompareSnglBurstByStartTime()}
\idx{XLALCompareSnglBurstByTime()}
\idx{XLALCompareSnglBurstByStartTimeAndLowFreq()}
\idx{LALClusterSnglBurstTable()}
\idx{XLALClusterSnglBurstTable()}

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
 * Sort a list of SnglBurstTable events into increasing order according to the
 * supplied comparison function.
 */

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


/*
 * Compare the start times of two SnglBurstTable events.
 */

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


/*
 * Compare the time intervals of two SnglBurstTable events.
 */

/* <lalVerbatim file="SnglBurstUtilsCP"> */
int
XLALCompareSnglBurstByTime(
	const SnglBurstTable **a,
	const SnglBurstTable **b
)
/* </lalVerbatim> */
{
	if(start_time(*a) > end_time(*b))
		return(1);	/* a's interval lies after b's interval */
	else if(end_time(*a) < start_time(*b))
		return(-1);	/* a's interval lies before b's interval */
	return(0);	/* a and b's intervals overlap */
}


/*
 * Compare the low frequency limits of two SnglBurstTable events.
 */

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


/*
 * Compare the frequency bands of two SnglBurstTable events.
 */

/* <lalVerbatim file="SnglBurstUtilsCP"> */
int
XLALCompareSnglBurstByFreq(
	const SnglBurstTable **a,
	const SnglBurstTable **b
)
/* </lalVerbatim> */
{
	if(lo_freq(*a) > hi_freq(*b))
		return(1);	/* a's band lies above b's band */
	else if(hi_freq(*a) < lo_freq(*b))
		return(-1);	/* a's band lies below b's band */
	return(0);	/* a and b's bands overlap */
}


/*
 * Compare two events first by start time, then by lowest frequency for ties.
 */

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


/*
 * Check to see if two events overlap in time and frequency.
 */

/* <lalVerbatim file="SnglBurstUtilsCP"> */
void
LALCompareSnglBurst(
	LALStatus *status,
	const SnglBurstTable *a,
	const SnglBurstTable *b,
	int *match
)
/* </lalVerbatim> */
{
	INITSTATUS (status, "LALCompareSnglBurst", SNGLBURSTUTILSC);

	*match = !XLALCompareSnglBurstByTime(&a, &b) && !XLALCompareSnglBurstByFreq(&a, &b);

	RETURN(status);
}


/*
 * Check to see if a SimBurstTable event lies within the tile defined by the
 * given SnglBurstTable event.
 */

/* <lalVerbatim file="SnglBurstUtilsCP"> */
void
LALCompareSimBurstAndSnglBurst(
	LALStatus *status,
	const SimBurstTable *a,
	const SnglBurstTable *b,
	int *match
)
/* </lalVerbatim> */
{
	INT8 ta;

	INITSTATUS(status, "LALCompareSimBurstAndSnglBurst", SNGLBURSTUTILSC);

	ta = XLALGPStoINT8(&a->geocent_peak_time);

	*match = (start_time(b) < ta) && (ta < end_time(b)) && (lo_freq(b) < a->freq) && (a->freq < hi_freq(b));

	RETURN(status);
}


/*
 * Returns true if two events should be clustered (peak times are within 10 ns
 * of each other and their frequency bands overlap).
 */

static int should_cluster(const SnglBurstTable *a, const SnglBurstTable *b)
{
	const REAL8 epsilon = 1e-8;	/* seconds */

	return((fabs(XLALDeltaFloatGPS(&a->peak_time, &b->peak_time)) < epsilon) && !XLALCompareSnglBurstByFreq(&a, &b));
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


/*
 * Recursively cluster a linked list of SnglBurstTable events until the list
 * stops changing.
 */

/* <lalVerbatim file="SnglBurstUtilsCP"> */
void
XLALClusterSnglBurstTable (
	SnglBurstTable   *list
)
/* </lalVerbatim> */
{
	int did_cluster;
	SnglBurstTable *a, *b, *prev;

	do {
		did_cluster = 0;

		for(a = list; a; a = a->next)
			for(prev = a, b = a->next; b; b = prev->next) {
				if(should_cluster(a, b)) {
					XLALClusterSnglBurst(a, b);
					prev->next = b->next;
					LALFree(b);
					did_cluster = 1;
				} else
					prev = b;
			}
	} while(did_cluster);
}


/* <lalVerbatim file="SnglBurstUtilsCP"> */
void
LALClusterSnglBurstTable (
	LALStatus        *status,
	SnglBurstTable   *list
)
/* </lalVerbatim> */
{
	INITSTATUS (status, "LALClusterSnglBurstTable", SNGLBURSTUTILSC);
	XLALClusterSnglBurstTable(list);
	RETURN(status);
}

