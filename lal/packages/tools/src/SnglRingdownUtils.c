/*
*  Copyright (C) 2007 Lisa M. Goggin
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
 * File Name: SnglRingdownUtils.c
 *
 * Author: Goggin, L. M based on  SnglInspiralUtils by Brady, P. R., Brown, D. A., Fairhurst, S. and Messaritaki, E.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="SnglRingdownUtilsCV">
Author: Brown, D. A., Fairhurst, S. and Messaritaki, E.
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
#include <lal/SkyCoordinates.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/Ring.h>

NRCSID( SNGLRINGDOWNUTILSC, "$Id$" );

#if 0
<lalLaTeX>
\subsection{Module \texttt{SnglRingdownUtils.c}}

Provides a set of utilities for manipulating \texttt{snglRingdownTable}s.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{SnglRingdownUtilsCP}
\idx{LALSortSnglRingdown()}
\idx{LALCompareSnglRingdownByTime()}
\idx{LALCompareSnglRingdown()}
\idx{LALClusterSnglRingdownTable()}
\idx{LALTimeCutSingleRingdown()}
\idx{LALIfoCountSingleRingdown()}
\idx{LALTimeSlideSingleRingdown()}
\idx{LALPlayTestSingleRingdown()}
\idx{LALCreateTrigBank()}


\subsubsection*{Description}

The function \texttt{LALFreeSnglInspiral()} frees the memory associated to a
single inspiral table.  The single inspiral table may point to a linked list
of EventIDColumns.  Thus, it is necessary to free all event ids associated
with the single inspiral.

The function \texttt{LALSortSnglInspiral()} sorts a list of single inspiral
tables.  The function simply calls qsort with the appropriate comparison
function, \texttt{comparfunc}.  It then ensures that the head of the sorted
list is returned.  There then follow several comparison functions for single
inspiral tables.  \texttt{LALCompareSnglInspiralByMass ()} first compares the
\texttt{mass1} entry of the two inspiral tables, returning 1 if the first mass
is larger and -1 if the second is larger.  In the case that the \texttt{mass1}
fields are equal, a similar comparsion is performed on \texttt{mass2}.  If
these also agree, 0 is returned.  \texttt{LALCompareSnglInspiralByPsi()}
compares the \texttt{Psi0} and \texttt{Psi3} fields in two single inspiral
tables.  The function is analogous to the mass comparison described above.
\texttt{LALCompareSnglInspiralByTime} compares the end times of two single
inspiral tables, returnng 1 if the first time is larger, 0 if equal and -1 if
the second time is larger.

\texttt{LALCompareSnglInspiral()} tests whether two single inspiral tables
pass a coincidence test.  The coincidence parameters are given by
\texttt{params} which is a \texttt{SnglInspiralAccuracy} structure.  It tests
first that the \texttt{ifo} fields are different.  If they are, it then tests
for time and mass coincidence, where mass coincidence may be any one of
\texttt{psi0\_and\_psi3}, \texttt{m1\_and\_m2}, \texttt{mchirp\_and\_eta}.
Finally, if the test is on \texttt{m1\_and\_m2}, consistency of effective
distances is also checked.  If the two single inspiral tables pass
coincidences the \texttt{params.match} is set to 1, otherwise it is set to
zero.

\texttt{LALClusterSnglInspiralTable ()} clusters single inspiral triggers
within a time window \texttt{dtimeNS}.  The triggers are compared either by
\texttt{snr}, \texttt{snr\_and\_chisq} or \texttt{snrsq\_over\_chisq}.  The
"loudest" trigger, as determined by the selected algorithm, within each time
window is returned.

\texttt{LALTimeCutSingleInspiral()} takes in a linked list of single inspiral
tables and returns only those which occur after the given \texttt{startTime}
and before the \texttt{endTime}.

\texttt{LALalphaFCutSingleInspiral()} takes in a linked list of single
inspiral tables and returns only those triggers which have alphaF values below
a specific alphaFcut. It is relevant for the BCV search only.

\texttt{LALIfoCutSingleInspiral()} scans through a linked list of single
inspiral tables and returns those which are from the requested \texttt{ifo}.
On input, \texttt{eventHead} is a pointer to the head of a linked list of
single inspiral tables.  On output, this list contains only single inspirals
from the requested \texttt{ifo}.

\texttt{LALIfoCountSingleInspiral()} scans through a linked list of single
inspiral tables and counts the number which are from the requested IFO.
This count is returned as \texttt{numTrigs}.

\texttt{LALTimeSlideSingleInspiral()} performs a time slide on the triggers
contained in the \texttt{triggerList}.  The time slide for each instrument is
specified by \texttt{slideTimes[LAL\_NUM\_IFO]}.  If \texttt{startTime} and
\texttt{endTime} are specified, then the time slide is performed on a ring.  If
the slide takes any trigger outside of the window
\texttt{[startTime,endTime]}, then the trigger is wrapped to be in
this time window.

\texttt{LALPlayTestSingleInspiral()} tests whether single inspiral events
occured in playground or non-playground times.  It then returns the requested
subset of events which occurred in the times specified by \texttt{dataType}
which must be one of \texttt{playground\_only}, \texttt{exclude\_play} or
\texttt{all\_data}.

\texttt{LALCreateTrigBank()} takes in a list of single inspiral tables and
returns a template bank.  The function tests whether a given template produced
multiple triggers.  If it did, only one copy of the template is retained.
Triggers are tested for coincidence in \texttt{m1\_and\_m2} or
\texttt{psi0\_and\_psi3}.


\subsubsection*{Algorithm}

\noindent None.

\subsubsection*{Uses}

\noindent LALCalloc, LALFree, LALGPStoINT8, LALINT8NanoSecIsPlayground.

\subsubsection*{Notes}
%% Any relevant notes.

\vfill{\footnotesize\input{SnglInspiralUtilsCV}}

</lalLaTeX>
#endif

/*
 * A few quickies for convenience.
 */

static INT8 start_time(const SnglRingdownTable *x)
{
	return(XLALGPStoINT8(&x->start_time));
}

static INT4 start_time_sec(const SnglRingdownTable *x)
{
	return(x->start_time.gpsSeconds);
}

static INT4 start_time_nsec(const SnglRingdownTable *x)
{
	return(x->start_time.gpsNanoSeconds);
}

/* <lalVerbatim file="SnglRingdownUtilsCP"> */
void
LALFreeSnglRingdown (
    LALStatus          *status,
    SnglRingdownTable **eventHead
    )
/* </lalVerbatim> */
{
  INITSTATUS( status, "LALFreeSnglRingdown", SNGLRINGDOWNUTILSC );
  XLALFreeSnglRingdown( eventHead );
  RETURN( status );
}

/* <lalVerbatim file="SnglRingdownUtilsCP"> */
int
XLALFreeSnglRingdown (
    SnglRingdownTable **eventHead
    )
/* </lalVerbatim> */
{
  EventIDColumn        *eventId;
  CoincRingdownTable   *thisCoinc;
  InterferometerNumber  ifoNumber;

  while ( (eventId = (*eventHead)->event_id) )
  {
    /* free any associated event_id's */
    (*eventHead)->event_id = (*eventHead)->event_id->next;

    if( (thisCoinc = eventId->coincRingdownTable) )
    {
      /* this Sngl is still part of a coinc, set pointer to NULL */
      for ( ifoNumber = 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
      {
        if ( *eventHead == thisCoinc->snglRingdown[ifoNumber] )
        {
          thisCoinc->snglRingdown[ifoNumber] = NULL;
        }
      }
    }
    LALFree( eventId );
  }
  LALFree( *eventHead );

  return (0);
}


/* <lalVerbatim file="SnglRingdownUtilsCP"> */
void
LALSortSnglRingdown (
    LALStatus          *status,
    SnglRingdownTable **eventHead,
    int(*comparfunc)    (const void *, const void *)
    )
/* </lalVerbatim> */
{
  INITSTATUS( status, "LALSortSnglRingdown", SNGLRINGDOWNUTILSC );

  *eventHead = XLALSortSnglRingdown ( *eventHead, comparfunc );

  RETURN( status );
}

/* <lalVerbatim file="SnglRingdownUtilsCP"> */
SnglRingdownTable *
XLALSortSnglRingdown (
    SnglRingdownTable *eventHead,
    int(*comparfunc)   (const void *, const void *)
    )
/* </lalVerbatim> */
{
  INT4                  i;
  INT4                  numEvents = 0;
  SnglRingdownTable    *thisEvent = NULL;
  SnglRingdownTable   **eventHandle = NULL;

  /* count the number of events in the linked list */
  for ( thisEvent = eventHead; thisEvent; thisEvent = thisEvent->next )
  {
    ++numEvents;
  }
  if ( ! numEvents )
  {
    XLALPrintInfo(
        "XLALSortSnglRingdown: Empty SnglRingdownTable passed as input" );
    return( eventHead );
  }

  /* allocate memory for an array of pts to sort and populate array */
  eventHandle = (SnglRingdownTable **)
    LALCalloc( numEvents, sizeof(SnglRingdownTable *) );
  for ( i = 0, thisEvent = eventHead; i < numEvents;
      ++i, thisEvent = thisEvent->next )
  {
    eventHandle[i] = thisEvent;
  }

  /* qsort the array using the specified function */
  qsort( eventHandle, numEvents, sizeof(eventHandle[0]), comparfunc );

  /* re-link the linked list in the right order */
  thisEvent = eventHead = eventHandle[0];
  for ( i = 1; i < numEvents; ++i )
  {
    thisEvent = thisEvent->next = eventHandle[i];
  }
  thisEvent->next = NULL;

  /* free the internal memory */
  LALFree( eventHandle );

  return( eventHead );
}



/* <lalVerbatim file="SnglRingdownUtilsCP"> */
int
LALCompareSnglRingdownByTime (
    const void *a,
    const void *b
    )
/* </lalVerbatim> */
{
  LALStatus     status;
  const SnglRingdownTable *aPtr = *((const SnglRingdownTable * const *)a);
  const SnglRingdownTable *bPtr = *((const SnglRingdownTable * const *)b);
  INT8 ta, tb;

  memset( &status, 0, sizeof(LALStatus) );
  LALGPStoINT8( &status, &ta, &(aPtr->start_time) );
  LALGPStoINT8( &status, &tb, &(bPtr->start_time) );

  if ( ta > tb )
  {
    return 1;
  }
  else if ( ta < tb )
  {
    return -1;
  }
  else
  {
    return 0;
  }
}



/* <lalVerbatim file="SnglInspiralUtilsCP"> */
void
LALCompareRingdowns (
    LALStatus                *status,
    SnglRingdownTable        *aPtr,
    SnglRingdownTable        *bPtr,
    RingdownAccuracyList     *params
    )
/* </lalVerbatim> */
{
  INT8    ta,  tb;
  REAL4   df, dQ;
  REAL4   fa, fb, Qa, Qb;
  REAL4   dsab = 0;
  REAL4   dsba = 0;
  REAL8   step = 1./params->minimizerStep;
  InterferometerNumber ifoaNum,  ifobNum;
  SnglRingdownAccuracy aAcc, bAcc;

  INITSTATUS( status, "LALCompareRingdowns", SNGLRINGDOWNUTILSC );
  ATTATCHSTATUSPTR( status );


  params->match = 1;

  /* check that triggers come from different IFOs */
  if( strcmp(aPtr->ifo, bPtr->ifo) )
  {
    LALInfo( status, "Triggers from different IFOs");
    params->match = 1;
  }
  else
  {
    LALInfo( status, "Triggers from same IFO");
    params->match = 0;
    goto exit;
  }

  ifoaNum = XLALIFONumber( aPtr->ifo );
  ifobNum = XLALIFONumber( bPtr->ifo );

  LALGPStoINT8( status->statusPtr, &ta, &(aPtr->start_time) );
  LALGPStoINT8( status->statusPtr, &tb, &(bPtr->start_time) );

  /* compare on trigger time coincidence */
  aAcc = params->ifoAccuracy[ifoaNum];
  bAcc = params->ifoAccuracy[ifobNum];

  if ( labs( ta - tb ) < (aAcc.dt + bAcc.dt)
      + params->lightTravelTime[ifoaNum][ifobNum])
  {
    LALInfo( status, "Triggers pass time coincidence test");
    params->match = 1;
  }
  else
  {
    LALInfo( status, "Triggers fail time coincidence test" );
    params->match = 0;
    goto exit;
  }

  /* compare f and Q parameters */
  if ( params->test == f_and_Q )
  {
    df = fabs( aPtr->frequency - bPtr->frequency );
    dQ = fabs( aPtr->quality - bPtr->quality );

    if ( ( df <= (aAcc.df + bAcc.df) )
        && ( dQ <= (aAcc.dQ + bAcc.dQ) ))
    {
      LALInfo( status, "Triggers are coincident in f and Q" );
      params->match = 1;
    }
    else
    {
      LALInfo( status, "Triggers are not coincident in f and Q" );
      params->match = 0;
      goto exit;
    }
  }
  else if ( params->test == ds_sq || params->test == ds_sq_fQt )
  {
    fa = aPtr->frequency;
    fb = bPtr->frequency;
    Qa = aPtr->quality;
    Qb = bPtr->quality;

    if ( params->test == ds_sq )
    {
      dsab = XLAL2DRingMetricDistance( fa, fb, Qa, Qb );
      dsba = XLAL2DRingMetricDistance( fb, fa, Qb, Qa );
    }
    else if ( params->test == ds_sq_fQt )
    {
      REAL8 dtab = 1.e-9 * (tb - ta);
      REAL8 dt_min = dtab - 1.e-9 * fabs(params->lightTravelTime[ifoaNum][ifobNum]);
      REAL8 dt_max = dtab + 1.e-9 * fabs(params->lightTravelTime[ifoaNum][ifobNum]);
      REAL4 ds2_min = XLAL3DRingMetricDistance( fa, fb, Qa, Qb, dtab );
      REAL8 dt;

      /* estimate true time delay */
      for ( dt = dt_min ; dt < dt_max ; dt += step )
      {
        REAL4 ds2 = XLAL3DRingMetricDistance( fa, fb, Qa, Qb, dt );
        if (ds2 < ds2_min) ds2_min = ds2;
      }

      dsab = ds2_min;
      dsba = ds2_min;
    }
    if ( (dsab + dsba)/2. < (aAcc.ds_sq + bAcc.ds_sq)/2. )
    {
      LALInfo( status, "Triggers pass the ds_sq coincidence test" );
      params->match = 1;
      if ( (strcmp(aPtr->ifo,"H1")==0 && strcmp(bPtr->ifo,"H2")==0)
		||(strcmp(aPtr->ifo,"H2")==0 && strcmp(bPtr->ifo,"H1")==0) )
      {
        aPtr->ds2_H1H2=dsab;
        bPtr->ds2_H1H2=dsba;
      }
      else if( (strcmp(aPtr->ifo,"H1")==0 && strcmp(bPtr->ifo,"L1")==0)
 		|| (strcmp(aPtr->ifo,"L1")==0 && strcmp(bPtr->ifo,"H1")==0) )
      {
        aPtr->ds2_H1L1=dsab;
        bPtr->ds2_H1L1=dsba;
      }
      else if( (strcmp(aPtr->ifo,"H2")==0 && strcmp(bPtr->ifo,"L1")==0)
		|| (strcmp(aPtr->ifo,"L1")==0 && strcmp(bPtr->ifo,"H2")==0) )
      {
        aPtr->ds2_H2L1=dsab;
        bPtr->ds2_H2L1=dsba;
      }
      else
      {
        LALInfo( status, "Unknown pair of ifo's" );
        params->match = 0;
        goto exit;
      }
    }
    else
    {
      LALInfo( status, "Triggers fail the ds_sq coincidence test" );
      params->match = 0;
      goto exit;
    }
  }
  else
  {
    LALInfo( status, "error: unknown test\n" );
    params->match = 0;
    goto exit;
  }

exit:
  DETATCHSTATUSPTR (status);
  RETURN (status);
}



/* <lalVerbatim file="SnglRingdownUtilsCP"> */
void
LALClusterSnglRingdownTable (
    LALStatus                  *status,
    SnglRingdownTable          *ringdownEvent,
    INT8                        dtimeNS,
    SnglInspiralClusterChoice   clusterchoice
    )
/* </lalVerbatim> */
{
  SnglRingdownTable     *thisEvent=NULL;
  SnglRingdownTable     *prevEvent=NULL;

  INITSTATUS( status, "LALClusterSnglRingdownTable", SNGLRINGDOWNUTILSC );
  ATTATCHSTATUSPTR( status );

  ASSERT( ringdownEvent, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  thisEvent = ringdownEvent->next;
  prevEvent = ringdownEvent;

  while ( thisEvent )
  {
    INT8 currTime;
    INT8 prevTime;

    /* compute the time in nanosec for each event trigger */
    LALGPStoINT8(status->statusPtr, &currTime, &(thisEvent->start_time));
    CHECKSTATUSPTR(status);

    LALGPStoINT8(status->statusPtr, &prevTime, &(prevEvent->start_time));
    CHECKSTATUSPTR(status);

    /* find events within the cluster window */
    if ( (currTime - prevTime) < dtimeNS )
    {
      /* displace previous event in cluster */
      if ( (clusterchoice == snr) && (thisEvent->snr > prevEvent->snr))
      {
        memcpy( prevEvent, thisEvent, sizeof(SnglRingdownTable) );
        thisEvent->event_id = NULL;
      }
      /* otherwise just dump this event from cluster */
      prevEvent->next = thisEvent->next;
      LALFreeSnglRingdown ( status->statusPtr, &thisEvent );
      thisEvent = prevEvent->next;
    }
    else
    {
      /* otherwise we keep this unique event trigger */
      prevEvent = thisEvent;
      thisEvent = thisEvent->next;
    }
  }

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}

/* <lalVerbatim file="SnglRingdownUtilsCP"> */
SnglRingdownTable *
XLALVetoSingleRingdown (
    SnglRingdownTable          *eventHead,
    LALSegList                 *vetoSegs,
    CHAR                        *ifo
    )
/* </lalVerbatim> */
{
  SnglRingdownTable    *thisEvent = NULL;
  SnglRingdownTable    *prevEvent = NULL;

  thisEvent = eventHead;
  eventHead = NULL;

  while ( thisEvent )
  {
    /*-- Check the time of this event against the veto segment list --*/
    if ( XLALSegListSearch( vetoSegs, &(thisEvent->start_time) )
        && (strcmp(thisEvent->ifo, ifo)==0) )
    {
      /*-- This event's start_time falls within one of the veto segments --*/
      /* discard the trigger and move to the next one */
      SnglRingdownTable    *tmpEvent = NULL;
      if ( prevEvent ) prevEvent->next = thisEvent->next;
      tmpEvent = thisEvent;
      thisEvent = thisEvent->next;
      XLALFreeSnglRingdown ( &tmpEvent );
    }
    else
    {
      /* This ringdown trigger does not fall within any veto segment */
      /* keep the trigger and increment the count of triggers */
      if ( ! eventHead ) eventHead = thisEvent;
      prevEvent = thisEvent;
      thisEvent = thisEvent->next;
    }
  }
  return( eventHead );
}

/* <lalVerbatim file="SnglRingdownUtilsCP"> */
void
LALIfoCutSingleRingdown(
    LALStatus                  *status,
    SnglRingdownTable         **eventHead,
    CHAR                       *ifo
    )
/* </lalVerbatim> */
{
  SnglRingdownTable    *ifoHead   = NULL;
  SnglRingdownTable    *thisEvent = NULL;

  INITSTATUS( status, "LALIfoCutSingleRingdown", SNGLRINGDOWNUTILSC );
  ATTATCHSTATUSPTR( status );

  ifoHead = XLALIfoCutSingleRingdown( eventHead, ifo );

  /* free events from other ifos */
  while ( *eventHead )
  {
    thisEvent = *eventHead;
    *eventHead = (*eventHead)->next;

    XLALFreeSnglRingdown( &thisEvent );
  }

  *eventHead = ifoHead;
  DETATCHSTATUSPTR (status);
  RETURN (status);
}

/* <lalVerbatim file="SnglInspiralUtilsCP"> */
SnglRingdownTable *
XLALIfoCutSingleRingdown(
    SnglRingdownTable         **eventHead,
    char                       *ifo
    )
/* </lalVerbatim> */
{
  static const char *func = "IfoCutSingleRingdown";
  SnglRingdownTable    *prevEvent   = NULL;
  SnglRingdownTable    *thisEvent   = NULL;
  SnglRingdownTable    *ifoHead     = NULL;
  SnglRingdownTable    *thisIfoTrig = NULL;

  /* check that eventHead is non-null */
  if ( ! eventHead )
    {
      XLAL_ERROR_NULL(func,XLAL_EIO);
     }
  /* Scan through a linked list of sngl_ringdown tables and return a
     pointer to the head of a linked list of tables for a specific IFO */

  thisEvent = *eventHead;
  *eventHead = NULL;

  while ( thisEvent )
  {
    if ( ! strcmp( thisEvent->ifo, ifo ) )
    {
      /* ifos match so keep this event */
      if ( ifoHead  )
      {
        thisIfoTrig = thisIfoTrig->next = thisEvent;
      }
      else
      {
        ifoHead = thisIfoTrig = thisEvent;
       }

      /* remove from eventHead list */
      if ( prevEvent )
      {
        prevEvent->next = thisEvent->next;
      }

      /* move to next event */
      thisEvent = thisEvent->next;
      /* terminate ifo list */
      thisIfoTrig->next = NULL;
    }
    else
    {
      /* move along the list */
      if ( ! *eventHead )
      {
        *eventHead = thisEvent;
      }

      prevEvent = thisEvent;
      thisEvent = thisEvent->next;
    }
  }

  return( ifoHead );
}


/* <lalVerbatim file="SnglRingdownUtilsCP"> */
void
LALTimeCutSingleRingdown(
    LALStatus                  *status,
    SnglRingdownTable         **eventHead,
    LIGOTimeGPS                *startTime,
    LIGOTimeGPS                *endTime
    )
/* </lalVerbatim> */
{
  INITSTATUS( status, "LALTimeCutSingleRingdown", SNGLRINGDOWNUTILSC );
  ATTATCHSTATUSPTR( status );

  *eventHead = XLALTimeCutSingleRingdown( *eventHead, startTime, endTime );

  DETATCHSTATUSPTR (status);
  RETURN (status);

}

/* <lalVerbatim file="SnglRingdownUtilsCP"> */

SnglRingdownTable *
XLALTimeCutSingleRingdown(
    SnglRingdownTable          *eventHead,
    LIGOTimeGPS                *startTime,
    LIGOTimeGPS                *endTime
    )
/* </lalVerbatim> */
{
  SnglRingdownTable    *ringdownEventList = NULL;
  SnglRingdownTable    *thisEvent = NULL;
  SnglRingdownTable    *prevEvent = NULL;
  INT8                  startTimeNS = XLALGPStoINT8( startTime );
  INT8                  endTimeNS = XLALGPStoINT8( endTime );


  /* Remove all the triggers before and after the requested */
  /* gps start and end times */

  thisEvent = eventHead;

  while ( thisEvent )
  {
    SnglRingdownTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;

     if ( start_time(tmpEvent) >= startTimeNS &&
         start_time(tmpEvent) < endTimeNS )
       {
         /* keep this template */
         if ( ! ringdownEventList  )
         {
           ringdownEventList = tmpEvent;
         }
         else
         {
           prevEvent->next = tmpEvent;
         }
         tmpEvent->next = NULL;
         prevEvent = tmpEvent;
         }
     else
       {
         /* discard this template */
         XLALFreeSnglRingdown ( &tmpEvent );
         }
     }
  eventHead = ringdownEventList;

  return (eventHead);
}


/* <lalVerbatim file="SnglRingdownUtilsCP"> */
void
LALIfoCountSingleRingdown(
    LALStatus                  *status,
    UINT4                      *numTrigs,
    SnglRingdownTable          *input,
    InterferometerNumber        ifoNumber
    )
/* </lalVerbatim> */
{
  SnglRingdownTable    *thisEvent = NULL;

  INITSTATUS( status, "LALIfoCountSingleRingdown", SNGLRINGDOWNUTILSC );
  ATTATCHSTATUSPTR( status );

  /* check that output is null and input non-null */
  ASSERT( !(*numTrigs), status,
      LIGOMETADATAUTILSH_ENNUL, LIGOMETADATAUTILSH_MSGENNUL );
  ASSERT( input, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  /* Scan through a linked list of sngl_ringdown tables and return a
     pointer to the head of a linked list of tables for a specific IFO */
  for( thisEvent = input; thisEvent; thisEvent = thisEvent->next )
  {
    if ( ifoNumber == XLALIFONumber(thisEvent->ifo) )
    {
      /* IFOs match so count this trigger */
      ++(*numTrigs);
    }
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}

/* <lalVerbatim file="SnglRingdownUtilsCP"> */
void
LALTimeSlideSingleRingdown(
    LALStatus                  *status,
    SnglRingdownTable          *triggerList,
    LIGOTimeGPS                *startTime,
    LIGOTimeGPS                *endTime,
    LIGOTimeGPS                 slideTimes[LAL_NUM_IFO]
    )
/* </lalVerbatim> */
{
  SnglRingdownTable    *thisEvent   = NULL;
  INT8                  startTimeNS = 0;
  INT8                  endTimeNS   = 0;
  INT8                  slideNS     = 0;
  INT8                  trigTimeNS  = 0;
  INITSTATUS( status, "LALTimeSlideSingleRingdown", SNGLRINGDOWNUTILSC );
  ATTATCHSTATUSPTR( status );

  /* time slide triggers by a time = slideTime, except those from the
   * instrument skipIfo which are left untouched. If you want to slide
   * all triggers, simply set skipIfo = LAL_UNKNOWN_IFO */


  /* check that input non-null */
  ASSERT( triggerList, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  if ( startTime )
  {
    LALGPStoINT8( status->statusPtr, &startTimeNS, startTime );
  }

  if ( endTime )
  {
    LALGPStoINT8( status->statusPtr, &endTimeNS, endTime );
  }

  for( thisEvent = triggerList; thisEvent; thisEvent = thisEvent->next )
  {
    /* calculate the slide time in nanoseconds */
    LALGPStoINT8( status->statusPtr, &slideNS,
        &(slideTimes[XLALIFONumber(thisEvent->ifo)]) );
    /* and trig time in nanoseconds */
    LALGPStoINT8( status->statusPtr, &trigTimeNS, &(thisEvent->start_time));
    trigTimeNS += slideNS;

    while ( startTimeNS && trigTimeNS < startTimeNS )
    {
      /* if before startTime, then wrap trigger time */
      trigTimeNS += endTimeNS - startTimeNS;
    }
    while ( endTimeNS && trigTimeNS > endTimeNS )
    {
      /* if after endTime, then wrap trigger time */
      trigTimeNS -= endTimeNS - startTimeNS;
    }

    /* convert back to LIGOTimeGPS */
    LALINT8toGPS( status->statusPtr, &(thisEvent->start_time), &trigTimeNS );
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}



/* <lalVerbatim file="SnglRingdownUtilsCP"> */
SnglRingdownTable *
XLALPlayTestSingleRingdown(
    SnglRingdownTable          *eventHead,
    LALPlaygroundDataMask      *dataType
    )
/* </lalVerbatim> */
{
  SnglRingdownTable    *ringdownEventList = NULL;
  SnglRingdownTable    *thisEvent = NULL;
  SnglRingdownTable    *prevEvent = NULL;

  INT8 triggerTime = 0;
  INT4 isPlay = 0;
  INT4 numTriggers;

  /* Remove all the triggers which are not of the desired type */

  numTriggers = 0;
  thisEvent = eventHead;

  if ( (*dataType == playground_only) || (*dataType == exclude_play) )
  {
    while ( thisEvent )
      {
        SnglRingdownTable *tmpEvent = thisEvent;
        thisEvent = thisEvent->next;

        triggerTime = XLALGPStoINT8( &(tmpEvent->start_time) );
        isPlay = XLALINT8NanoSecIsPlayground( &triggerTime );

        if ( ( (*dataType == playground_only)  && isPlay ) ||
            ( (*dataType == exclude_play) && ! isPlay) )
        {
          /* keep this trigger */
          if ( ! ringdownEventList  )
          {
            ringdownEventList = tmpEvent;
          }
          else
          {
            prevEvent->next = tmpEvent;
          }
          tmpEvent->next = NULL;
          prevEvent = tmpEvent;
          ++numTriggers;
          }
        else
          {
            /* discard this template */
            XLALFreeSnglRingdown ( &tmpEvent );
            }
        }
    eventHead = ringdownEventList;
    if ( *dataType == playground_only )
    {
      XLALPrintInfo( "Kept %d playground triggers \n", numTriggers );
    }
    else if ( *dataType == exclude_play )
    {
      XLALPrintInfo( "Kept %d non-playground triggers \n", numTriggers );
    }
  }
  else if ( *dataType == all_data )
  {
    XLALPrintInfo( "Keeping all triggers since all_data specified\n" );
  }
  else
  {
    XLALPrintInfo( "Unknown data type, returning no triggers\n" );
    eventHead = NULL;
   }

  return(eventHead);
}

/* <lalVerbatim file="SnglRingdownUtilsCP"> */
void
LALPlayTestSingleRingdown(
    LALStatus                  *status,
    SnglRingdownTable         **eventHead,
    LALPlaygroundDataMask      *dataType
    )
/* </lalVerbatim> */
{
  INITSTATUS( status, "LALPlayTestSingleRingdown", SNGLRINGDOWNUTILSC );
  ATTATCHSTATUSPTR( status );

  *eventHead = XLALPlayTestSingleRingdown(*eventHead, dataType);

  DETATCHSTATUSPTR (status);
  RETURN (status);
}

/* <lalVerbatim file="SnglRingdownUtilsCP"> */
int
XLALMaxSnglRingdownOverIntervals(
    SnglRingdownTable         **eventHead,
    INT8                       deltaT
    )
/* </lalVerbatim> */
{
  SnglRingdownTable    *ringdownEventList = NULL;
  SnglRingdownTable    *thisEvent = NULL;
  SnglRingdownTable    *nextEvent = NULL;
  SnglRingdownTable    *prevEvent = NULL;
  INT4 count = 1;

  /* if there are no events, then no-op */
  if ( ! *eventHead )
    return (0);

  ringdownEventList = *eventHead;
  thisEvent = *eventHead;
  nextEvent = thisEvent->next;

  while ( nextEvent )
  {
    if ( start_time_sec(nextEvent) == start_time_sec(thisEvent) &&
        start_time_nsec(nextEvent)/deltaT == start_time_nsec(thisEvent)/deltaT )
    {
      count = count + 1;
      if ( nextEvent->snr > thisEvent->snr )
      {
        /* replace thisEvent with nextEvent */
        XLALFreeSnglRingdown ( &thisEvent );

        /* deal with start of the list */
        if (prevEvent)
          prevEvent->next = nextEvent;
        else
          ringdownEventList = nextEvent;

        /* standard stuff */
        thisEvent = nextEvent;
        nextEvent = thisEvent->next;
      }
      else
      {
        /* get rid of nextEvent */
        thisEvent->next = nextEvent->next;
        XLALFreeSnglRingdown ( &nextEvent );
        nextEvent = thisEvent->next;
      }
    }
    else
    {
      thisEvent->num_clust_trigs = count;
      count = 1;
      /* step to next set of events */
      prevEvent=thisEvent;
      thisEvent=nextEvent;
      nextEvent = thisEvent->next;
    }
  }
  thisEvent->num_clust_trigs = count;
  *eventHead = ringdownEventList;

  return (0);
}

INT4 XLALCountSnglRingdown( SnglRingdownTable *head )
  /* </lalVerbatim> */
{
  INT4 length;
  SnglRingdownTable *event;

  if ( !head )
   {
     return( 0 );
   }

  /* count the number of events in the list */
  for(length = 0, event = head; event; event = event->next)
    length++;

  return length;
}






