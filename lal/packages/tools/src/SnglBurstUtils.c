/*
*  Copyright (C) 2007 Jolien Creighton, Kipp Cannon, Patrick Brady, Saikat Ray-Majumder, Xavier Siemens
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

#include <string.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/Date.h>
#include <lal/XLALError.h>

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
\idx{XLALCompareSnglBurstByStartTime()}
\idx{XLALCompareSnglBurstByLowFreq()}
\idx{XLALCompareSnglBurstByStartTimeAndLowFreq()}
\idx{XLALCompareSimBurstAndSnglBurst()}
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


static INT8 start_time(const SnglBurst *x)
{
	return(XLALGPSToINT8NS(&x->start_time));
}


static INT8 peak_time(const SnglBurst *x)
{
	return(XLALGPSToINT8NS(&x->peak_time));
}


static REAL4 lo_freq(const SnglBurst *x)
{
	return(x->central_freq - x->bandwidth / 2);
}


/**
 * Free a sngl_burst
 */


/* <lalVerbatim file="SnglBurstUtilsCP"> */
void XLALDestroySnglBurst(SnglBurst *event)
/* </lalVerbatim> */
{
	XLALFree(event);
}


/**
 * Free a SnglBurst linked list.
 */


void XLALDestroySnglBurstTable(SnglBurst *head)
{
	while(head) {
		SnglBurst *next = head->next;
		XLALDestroySnglBurst(head);
		head = next;
	}
}


/**
 * Assign event_id values to the entries in a sngl_burst linked list.  All
 * sngl_burst rows in the list will be blamed on the given process_id, and
 * assigned event_ids in order starting with the given event_id.  The
 * return value is the next event_id after the last one assigned to a row
 * in the list.
 */


long XLALSnglBurstAssignIDs(
	SnglBurst *head,
	long process_id,
	long event_id
)
{
	for(; head; head = head->next) {
		head->process_id = process_id;
		head->event_id = event_id++;
	}
	return event_id;
}


/*
 * Sort a list of SnglBurst events into increasing order according to the
 * supplied comparison function.
 */


/* <lalVerbatim file="SnglBurstUtilsCP"> */
int XLALSnglBurstTableLength(SnglBurst *head)
/* </lalVerbatim> */
{
	int length;

	/* count the number of events in the list */
	for(length = 0; head; head = head->next)
		length++;

	return(length);
}


/* <lalVerbatim file="SnglBurstUtilsCP"> */
void XLALSortSnglBurst(
	SnglBurst **head,
	int (*comparefunc)(const SnglBurst * const *, const SnglBurst * const *)
)
/* </lalVerbatim> */
{
	int i;
	int length;
	SnglBurst *event;
	SnglBurst **array;

	/* empty list --> no-op */
	if(!head || !*head)
		return;

	/* construct an array of pointers into the list */
	length = XLALSnglBurstTableLength(*head);
	array = XLALCalloc(length, sizeof(*array));
	for(i = 0, event = *head; event; event = event->next)
		array[i++] = event;

	/* sort the array using the specified function */
	qsort(array, length, sizeof(*array), (int(*)(const void *, const void *)) comparefunc);

	/* re-link the list according to the sorted array */
	for(i = 0; i < length; i++, head = &(*head)->next)
		*head = array[i];
	*head = NULL;

	/* free the array */
	XLALFree(array);
}


/*
 * Compare the start times of two SnglBurst events.
 */


/* <lalVerbatim file="SnglBurstUtilsCP"> */
int XLALCompareSnglBurstByStartTime(
	const SnglBurst * const *a,
	const SnglBurst * const *b
)
/* </lalVerbatim> */
{
	INT8 ta, tb;

	ta = start_time(*a);
	tb = start_time(*b);

	if(ta > tb)
		return(1);
	if(ta < tb)
		return(-1);
	return(0);
}


/*
 * Compare two sngl_burst events by their peak times, with no slack.
 */


/* <lalVerbatim file="SnglBurstUtilsCP"> */
int XLALCompareSnglBurstByExactPeakTime(
	const SnglBurst * const *a,
	const SnglBurst * const *b
)
/* </lalVerbatim> */
{
	INT8 ta, tb;

	ta = peak_time(*a);
	tb = peak_time(*b);

	if(ta > tb)
		return(1);
	if(ta < tb)
		return(-1);
	return(0);
}


/*
 * Compare the SNRs of two SnglBurst events.
 */


/* <lalVerbatim file="SnglBurstUtilsCP"> */
int XLALCompareSnglBurstBySNR(
	const SnglBurst * const *a,
	const SnglBurst * const *b
)
/* </lalVerbatim> */
{
	REAL4 snra = (*a)->snr;
	REAL4 snrb = (*b)->snr;

	if(snra < snrb)
		return 1;
	if(snra > snrb)
		return -1;
	return 0;
}


/*
 * Compare the peak times and SNRs of two SnglBurst events.
 */


/* <lalVerbatim file="SnglBurstUtilsCP"> */
int XLALCompareSnglBurstByPeakTimeAndSNR(
	const SnglBurst * const *a,
	const SnglBurst * const *b
)
/* </lalVerbatim> */
{
	int result;

	result = XLALCompareSnglBurstByExactPeakTime(a, b);
	if(!result)
		result = XLALCompareSnglBurstBySNR(a, b);

	return(result);
}


/*
 * Check if two string events overlap in time. The peak times are uncertain
 * to whatever the high frequency cutoff is
 */


/* <lalVerbatim file="SnglBurstUtilsCP"> */
int XLALCompareStringBurstByTime(
	const SnglBurst * const *a,
	const SnglBurst * const *b
)
/* </lalVerbatim> */
{
	INT8 ta, tb;
	INT8 epsilon;

	epsilon = (*a)->string_cluster_t * 1e9;
	if(epsilon < (*b)->string_cluster_t * 1e9)
		epsilon = (*b)->string_cluster_t * 1e9;

	ta = peak_time(*a);
	tb = peak_time(*b);
	if(ta > tb + epsilon)
		return(1);
	if(ta < tb - epsilon)
		return(-1);
	return(0);
}


/*
 * Compare the low frequency limits of two SnglBurst events.
 */


/* <lalVerbatim file="SnglBurstUtilsCP"> */
int XLALCompareSnglBurstByLowFreq(
	const SnglBurst * const *a,
	const SnglBurst * const *b
)
/* </lalVerbatim> */
{
	REAL4 flowa, flowb;

	flowa = lo_freq(*a);
	flowb = lo_freq(*b);

	if(flowa > flowb)
		return(1);
	if(flowa < flowb)
		return(-1);
	return(0);
}


/*
 * Compare two events first by start time, then by lowest frequency to break
 * ties.
 */


/* <lalVerbatim file="SnglBurstUtilsCP"> */
int XLALCompareSnglBurstByStartTimeAndLowFreq(
	const SnglBurst * const *a,
	const SnglBurst * const *b
)
/* </lalVerbatim> */
{
	int result;

	result = XLALCompareSnglBurstByStartTime(a, b);
	if(!result)
		result = XLALCompareSnglBurstByLowFreq(a, b);
	return(result);
}


/*
 * cluster events a and b, storing result in a; takes one with largest snr
 */


/* <lalVerbatim file="SnglBurstUtilsCP"> */
void XLALStringBurstCluster(
	SnglBurst *a,
	const SnglBurst *b
)
/* </lalVerbatim> */
{
	REAL4 snra=a->snr, snrb=b->snr;

	if(snrb > snra) {
		a->central_freq =  b->central_freq;
		a->bandwidth = b->bandwidth;
		a->start_time = b->start_time;
		a->duration = b->duration;
		a->amplitude = b->amplitude;
		a->snr = b->snr;
		a->confidence = b->confidence;
		a->peak_time = b->peak_time;
	}
}


/*
 * Recursively cluster a linked list of SnglBurst events until the list
 * stops changing.  testfunc() should return 0 if the two given events are to
 * be clustered.  If bailoutfunc() is provided (not NULL), then testfunc() will
 * be used to sort the trigger list before each clustering pass and
 * bailoutfunc() will be called to check for the option of terminating the
 * inner loop early.  In the ideal case, use of bailoutfunc() converts this
 * algorithm from O(n^3) to order O(n log n).  The clusterfunc() should replace
 * the SnglBurst event pointed to by its first argument with the cluster
 * of that event and the event pointed to by the second argument.
 */


/* <lalVerbatim file="SnglBurstUtilsCP"> */
void XLALClusterSnglBurstTable (
	SnglBurst **list,
	int (*bailoutfunc)(const SnglBurst * const *, const SnglBurst * const *),
	int (*testfunc)(const SnglBurst * const *, const SnglBurst * const *),
	void (*clusterfunc)(SnglBurst *, const SnglBurst *)
)
/* </lalVerbatim> */
{
	int did_cluster;
	SnglBurst *a, *b, *prev;

	do {
		did_cluster = 0;

		if(bailoutfunc)
			XLALSortSnglBurst(list, testfunc);

		for(a = *list; a; a = a->next)
			for(prev = a, b = a->next; b; b = prev->next) {
				if(!testfunc((const SnglBurst * const *) &a, (const SnglBurst * const *) &b)) {
					clusterfunc(a, b);
					prev->next = b->next;
					XLALDestroySnglBurst(b);
					did_cluster = 1;
				} else {
					if(bailoutfunc && bailoutfunc((const SnglBurst * const *) &a, (const SnglBurst * const *) &b))
						break;
					prev = b;
				}
			}
	} while(did_cluster);
}


/**
 * Create a SnglBurst structure.
 */


SnglBurst *XLALCreateSnglBurst(void)
{
	static const char func[] = "XLALCreateSnglBurst";
	SnglBurst *new = XLALMalloc(sizeof(*new));

	if(!new)
		XLAL_ERROR_NULL(func, XLAL_EFUNC);

	new->next = NULL;
	new->process_id = new->event_id = -1;
	memset(new->ifo, 0, sizeof(new->ifo));
	memset(new->search, 0, sizeof(new->search));
	memset(new->channel, 0, sizeof(new->channel));
	XLALGPSSet(&new->start_time, 0, 0);
	XLALGPSSet(&new->peak_time, 0, 0);
	new->duration = XLAL_REAL4_FAIL_NAN;
	new->central_freq = XLAL_REAL4_FAIL_NAN;
	new->bandwidth = XLAL_REAL4_FAIL_NAN;
	new->amplitude = XLAL_REAL4_FAIL_NAN;
	new->snr = XLAL_REAL4_FAIL_NAN;
	new->confidence = XLAL_REAL4_FAIL_NAN;
	/* FIXME:  remove string_cluster_t from the structure */
	new->string_cluster_t = XLAL_REAL4_FAIL_NAN;

	return new;
}


/**
 * Create a SimBurst structure.
 */


SimBurst *XLALCreateSimBurst(void)
{
	static const char func[] = "XLALCreateSimBurst";
	SimBurst *new = XLALMalloc(sizeof(*new));

	if(!new)
		XLAL_ERROR_NULL(func, XLAL_EFUNC);

	new->next = NULL;
	new->process_id = new->simulation_id = -1;
	memset(new->waveform, 0, sizeof(new->waveform));
	new->ra = XLAL_REAL8_FAIL_NAN;
	new->dec = XLAL_REAL8_FAIL_NAN;
	new->psi = XLAL_REAL8_FAIL_NAN;
	XLALGPSSet(&new->time_geocent_gps, 0, 0);
	new->time_geocent_gmst = XLAL_REAL8_FAIL_NAN;
	new->duration = XLAL_REAL8_FAIL_NAN;
	new->frequency = XLAL_REAL8_FAIL_NAN;
	new->bandwidth = XLAL_REAL8_FAIL_NAN;
	new->q = XLAL_REAL8_FAIL_NAN;
	new->pol_ellipse_angle = XLAL_REAL8_FAIL_NAN;
	new->pol_ellipse_e= XLAL_REAL8_FAIL_NAN;
	new->amplitude = XLAL_REAL8_FAIL_NAN;
	new->hrss = XLAL_REAL8_FAIL_NAN;
	new->egw_over_rsquared = XLAL_REAL8_FAIL_NAN;
	new->waveform_number = 0;

	return new;
}


/**
 * Destroy a SimBurst structure.
 */


void XLALDestroySimBurst(SimBurst *sim_burst)
{
	XLALFree(sim_burst);
}


/**
 * Destroy a SimBurst linked list.
 */


void XLALDestroySimBurstTable(SimBurst *head)
{
	while(head) {
		SimBurst *next = head->next;
		XLALDestroySimBurst(head);
		head = next;
	}
}


/**
 * Assign simulation_id values to the entries in a sim_burst linked list.
 * All sim_burst rows in the list will be blamed on the given process_id,
 * and assigned simulation_ids in order starting with the given
 * simulation_id.  The return value is the next simulation_id after the
 * last one assigned to a row in the list.
 */


long XLALSimBurstAssignIDs(
	SimBurst *sim_burst,
	long process_id,
	long simulation_id)
{
	for(; sim_burst; sim_burst = sim_burst->next) {
		sim_burst->process_id = process_id;
		sim_burst->simulation_id = simulation_id++;
	}
	return simulation_id;
}
