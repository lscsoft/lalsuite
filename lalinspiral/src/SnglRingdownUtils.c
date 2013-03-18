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
 *-----------------------------------------------------------------------
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataRingdownUtils.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/LALSimulation.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/RingUtils.h>

/**
\author Brown, D. A., Fairhurst, S. and Messaritaki, E.
\file

\brief Provides a set of utilities for manipulating \c snglRingdownTables.

\heading{Description}

The function <tt>LALFreeSnglInspiral()</tt> frees the memory associated to a
single inspiral table.  The single inspiral table may point to a linked list
of EventIDColumns.  Thus, it is necessary to free all event ids associated
with the single inspiral.

The function <tt>LALSortSnglInspiral()</tt> sorts a list of single inspiral
tables.  The function simply calls qsort with the appropriate comparison
function, \c comparfunc.  It then ensures that the head of the sorted
list is returned.  There then follow several comparison functions for single
inspiral tables.  <tt>LALCompareSnglInspiralByMass()</tt> first compares the
\c mass1 entry of the two inspiral tables, returning 1 if the first mass
is larger and -1 if the second is larger.  In the case that the \c mass1
fields are equal, a similar comparsion is performed on \c mass2.  If
these also agree, 0 is returned.  <tt>LALCompareSnglInspiralByPsi()</tt>
compares the \c Psi0 and \c Psi3 fields in two single inspiral
tables.  The function is analogous to the mass comparison described above.
\c LALCompareSnglInspiralByTime() compares the end times of two single
inspiral tables, returnng 1 if the first time is larger, 0 if equal and -1 if
the second time is larger.

<tt>LALCompareSnglInspiral()</tt> tests whether two single inspiral tables
pass a coincidence test.  The coincidence parameters are given by
\c params which is a \c ::SnglInspiralAccuracy structure.  It tests
first that the \c ifo fields are different.  If they are, it then tests
for time and mass coincidence, where mass coincidence may be any one of
\c psi0_and_psi3, \c m1_and_m2, \c mchirp_and_eta.
Finally, if the test is on \c m1_and_m2, consistency of effective
distances is also checked.  If the two single inspiral tables pass
coincidences the <tt>params.match</tt> is set to 1, otherwise it is set to
zero.

<tt>LALClusterSnglInspiralTable()</tt> clusters single inspiral triggers
within a time window \c dtimeNS.  The triggers are compared either by
\c snr, \c snr_and_chisq or \c snrsq_over_chisq.  The
"loudest" trigger, as determined by the selected algorithm, within each time
window is returned.

<tt>LALTimeCutSingleInspiral()</tt> takes in a linked list of single inspiral
tables and returns only those which occur after the given \c startTime
and before the \c endTime.

<tt>LALalphaFCutSingleInspiral()</tt> takes in a linked list of single
inspiral tables and returns only those triggers which have alphaF values below
a specific alphaFcut. It is relevant for the BCV search only.

<tt>LALIfoCutSingleInspiral()</tt> scans through a linked list of single
inspiral tables and returns those which are from the requested \c ifo.
On input, \c eventHead is a pointer to the head of a linked list of
single inspiral tables.  On output, this list contains only single inspirals
from the requested \c ifo.

<tt>LALIfoCountSingleInspiral()</tt> scans through a linked list of single
inspiral tables and counts the number which are from the requested IFO.
This count is returned as \c numTrigs.

<tt>XLALTimeSlideSingleInspiral()</tt> performs a time slide on the triggers
contained in the \c triggerList.  The time slide for each instrument is
specified by <tt>slideTimes[LAL_NUM_IFO]</tt>.  If \c startTime and
\c endTime are specified, then the time slide is performed on a ring.  If
the slide takes any trigger outside of the window
<tt>[startTime,endTime]</tt>, then the trigger is wrapped to be in
this time window.

<tt>LALPlayTestSingleInspiral()</tt> tests whether single inspiral events
occured in playground or non-playground times.  It then returns the requested
subset of events which occurred in the times specified by \c dataType
which must be one of \c playground_only, \c exclude_play or
\c all_data.

<tt>LALCreateTrigBank()</tt> takes in a list of single inspiral tables and
returns a template bank.  The function tests whether a given template produced
multiple triggers.  If it did, only one copy of the template is retained.
Triggers are tested for coincidence in \c m1_and_m2 or
\c psi0_and_psi3.

\heading{Algorithm}

None.

\heading{Uses}

LALCalloc(), LALFree(), LALINT8NanoSecIsPlayground().

*/

/*
 * A few quickies for convenience.
 */

static INT8 start_time(const SnglRingdownTable *x)
{
	return(XLALGPSToINT8NS(&x->start_time));
}

static INT4 start_time_sec(const SnglRingdownTable *x)
{
	return(x->start_time.gpsSeconds);
}

static INT4 start_time_nsec(const SnglRingdownTable *x)
{
	return(x->start_time.gpsNanoSeconds);
}


void
LALFreeSnglRingdown (
    LALStatus          *status,
    SnglRingdownTable **eventHead
    )

{
  INITSTATUS(status);
  XLALFreeSnglRingdown( eventHead );
  RETURN( status );
}


int
XLALFreeSnglRingdown (
    SnglRingdownTable **eventHead
    )

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
      for ( ifoNumber = (InterferometerNumber) 0; ifoNumber < LAL_NUM_IFO; ifoNumber++)
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



void
LALSortSnglRingdown (
    LALStatus          *status,
    SnglRingdownTable **eventHead,
    int(*comparfunc)    (const void *, const void *)
    )

{
  INITSTATUS(status);

  *eventHead = XLALSortSnglRingdown ( *eventHead, comparfunc );

  RETURN( status );
}


SnglRingdownTable *
XLALSortSnglRingdown (
    SnglRingdownTable *eventHead,
    int(*comparfunc)   (const void *, const void *)
    )

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




int
LALCompareSnglRingdownByTime (
    const void *a,
    const void *b
    )

{
  LALStatus     status;
  const SnglRingdownTable *aPtr = *((const SnglRingdownTable * const *)a);
  const SnglRingdownTable *bPtr = *((const SnglRingdownTable * const *)b);
  INT8 ta, tb;

  memset( &status, 0, sizeof(LALStatus) );
  ta = XLALGPSToINT8NS( &(aPtr->start_time) );
  tb = XLALGPSToINT8NS( &(bPtr->start_time) );

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

/** Compare theparameters of two ringdown triggers according to
 * the specified coincidence test.
 */

REAL8
XLALCompareRingdowns (
    SnglRingdownTable        *aPtr,
    SnglRingdownTable        *bPtr,
    RingdownAccuracyList     *params
    )

{
  SnglRingdownAccuracy aAcc, bAcc;
  InterferometerNumber ifoaNum,  ifobNum;

  REAL8 ds2 = 0;
  INT8 ta,  tb;
  const LALDetector *aDet;
  const LALDetector *bDet;

  ifoaNum = (InterferometerNumber) XLALIFONumber( aPtr->ifo );
  ifobNum = (InterferometerNumber) XLALIFONumber( bPtr->ifo );
  aAcc = params->ifoAccuracy[ifoaNum];
  bAcc = params->ifoAccuracy[ifobNum];
  ta = XLALGPSToINT8NS( &(aPtr->start_time) );
  tb = XLALGPSToINT8NS( &(bPtr->start_time) );
  params->match = 1;
  aDet = XLALDetectorPrefixToLALDetector(aPtr->ifo);
  bDet = XLALDetectorPrefixToLALDetector(bPtr->ifo);


  /* check that triggers come from different IFOs */
  if( strcmp(aPtr->ifo, bPtr->ifo) )
  {
    XLALPrintInfo( "Triggers from different IFOs");
    params->match = 1;
  }
  else
  {
    XLALPrintInfo( "Triggers from same IFO");
    params->match = 0;
    goto exit;
  }

  /* If f_and_Q or ds_sq test requested, */
  /* make sure triggers lie within a reasonable time window */
  if ( params->test == LALRINGDOWN_F_AND_Q || params->test == LALRINGDOWN_DS_SQ )
  {
     if ( labs( ta - tb ) < (aAcc.dt + bAcc.dt)
         + 1.e-9 * XLALLightTravelTime(aDet,bDet) )
     {
       params->match = 1;
     }
     else
     {
       params->match = 0;
       ds2 = 1.0 / 0.0;
       goto exit;
     }
  }

  /* compare f and Q parameters */
  if ( params->test == LALRINGDOWN_F_AND_Q )
  {
    if( (fabs( aPtr->frequency - bPtr->frequency ) <= (aAcc.df + bAcc.df) )
      && (fabs( aPtr->quality - bPtr->quality ) <= (aAcc.dQ + bAcc.dQ) ) )
    {
      params->match = 1;
    }
    else
    {
      params->match = 0;
      ds2 = 1.0 / 0.0;
    }
  }
  else if ( params->test == LALRINGDOWN_DS_SQ )
  {
    ds2 = XLAL2DRinca( aPtr, bPtr );
    if ( ds2 < (aAcc.ds_sq + bAcc.ds_sq)/2. )
    {
      params->match = 1;
    }
    else
    {
      params->match = 0;
    }
  }
  else if ( params->test == LALRINGDOWN_DS_SQ_FQT )
  {
    ds2 = XLAL3DRinca( aPtr, bPtr );
    if ( ds2 < (aAcc.ds_sq + bAcc.ds_sq)/2. )
    {
      params->match = 1;
    }
    else
    {
      params->match = 0;
    }
  }
  else
  {
    /* Invalid Coincidence Test */
    XLALPrintError( "error: unknown coincidence test\n" );
    XLAL_ERROR(XLAL_EIO);
  }
  exit:
  return ds2;
}

/** Two dimensional (frequency and quality) coincidence test */
REAL8
XLAL2DRinca(
    SnglRingdownTable         *aPtr,
    SnglRingdownTable         *bPtr
    )
{
  REAL8   fa, fb, Qa, Qb;
  REAL8   dsab = 0;
  REAL8   dsba = 0;

  fa = aPtr->frequency;
  fb = bPtr->frequency;
  Qa = aPtr->quality;
  Qb = bPtr->quality;

  dsab = XLAL2DRingMetricDistance( fa, fb, Qa, Qb );
  dsba = XLAL2DRingMetricDistance( fb, fa, Qb, Qa );

  return ( (dsab + dsba)/2. );
}

/** Three dimensional (time, frequency and quality) coincidence test */
REAL8
XLAL3DRinca(
    SnglRingdownTable         *aPtr,
    SnglRingdownTable         *bPtr
    )
{

  INT8    ta,  tb;
  REAL8   fa, fb, Qa, Qb, ds2_min, ds2;
  REAL8   dt_min, dt_max, dt, dt_best;
  REAL8   lightTravel;
  const LALDetector *aDet;
  const LALDetector *bDet;
  fa = aPtr->frequency;
  fb = bPtr->frequency;
  Qa = aPtr->quality;
  Qb = bPtr->quality;
  ta = XLALGPSToINT8NS( &(aPtr->start_time) );
  tb = XLALGPSToINT8NS( &(bPtr->start_time) );

  aDet = XLALDetectorPrefixToLALDetector(aPtr->ifo);
  bDet = XLALDetectorPrefixToLALDetector(bPtr->ifo);
  lightTravel = 1.e-9 * XLALLightTravelTime(aDet,bDet);

  dt = 1.e-9 * (ta - tb);
  dt_min = dt - lightTravel;
  dt_max = dt + lightTravel;

  ds2_min = XLAL3DRingMetricDistance( fa, fb, Qa, Qb, dt );

  /* if ifos are H1H2 then no need to account for light travel time */

  if ( (strcmp(aPtr->ifo,"H1")==0 && strcmp(bPtr->ifo,"H2")==0)
	||(strcmp(aPtr->ifo,"H2")==0 && strcmp(bPtr->ifo,"H1")==0) )
  {
    return ( ds2_min );
  }
  else
  {
    /* solve for the dt that minimizes ds2 */
    dt_best = XLAL3DRingTimeMinimum(fa, fb, Qa, Qb);

    /* check that this is a valid time */
    if ((dt_best < dt_max) && (dt_best > dt_min))
      {
       ds2_min = XLAL3DRingMetricDistance( fa, fb, Qa, Qb, dt_best);
      }
    else
      {
       ds2_min = XLAL3DRingMetricDistance( fa, fb, Qa, Qb, dt_min);
       ds2 = XLAL3DRingMetricDistance( fa, fb, Qa, Qb, dt_max);
       if(ds2 < ds2_min) ds2_min = ds2;
      }

    return ( ds2_min );
  }
}



void
LALClusterSnglRingdownTable (
    LALStatus                  *status,
    SnglRingdownTable          *ringdownEvent,
    INT8                        dtimeNS,
    SnglInspiralClusterChoice   clusterchoice
    )

{
  SnglRingdownTable     *thisEvent=NULL;
  SnglRingdownTable     *prevEvent=NULL;

  INITSTATUS(status);
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
    currTime = XLALGPSToINT8NS(&(thisEvent->start_time));

    prevTime = XLALGPSToINT8NS(&(prevEvent->start_time));

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


SnglRingdownTable *
XLALVetoSingleRingdown (
    SnglRingdownTable          *eventHead,
    LALSegList                 *vetoSegs,
    CHAR                        *ifo
    )

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


void
LALIfoCutSingleRingdown(
    LALStatus                  *status,
    SnglRingdownTable         **eventHead,
    CHAR                       *ifo
    )

{
  SnglRingdownTable    *ifoHead   = NULL;
  SnglRingdownTable    *thisEvent = NULL;

  INITSTATUS(status);
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


SnglRingdownTable *
XLALIfoCutSingleRingdown(
    SnglRingdownTable         **eventHead,
    char                       *ifo
    )

{
  SnglRingdownTable    *prevEvent   = NULL;
  SnglRingdownTable    *thisEvent   = NULL;
  SnglRingdownTable    *ifoHead     = NULL;
  SnglRingdownTable    *thisIfoTrig = NULL;

  /* check that eventHead is non-null */
  if ( ! eventHead )
    {
      XLAL_ERROR_NULL(XLAL_EIO);
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



void
LALTimeCutSingleRingdown(
    LALStatus                  *status,
    SnglRingdownTable         **eventHead,
    LIGOTimeGPS                *startTime,
    LIGOTimeGPS                *endTime
    )

{
  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  *eventHead = XLALTimeCutSingleRingdown( *eventHead, startTime, endTime );

  DETATCHSTATUSPTR (status);
  RETURN (status);

}



SnglRingdownTable *
XLALTimeCutSingleRingdown(
    SnglRingdownTable          *eventHead,
    LIGOTimeGPS                *startTime,
    LIGOTimeGPS                *endTime
    )

{
  SnglRingdownTable    *ringdownEventList = NULL;
  SnglRingdownTable    *thisEvent = NULL;
  SnglRingdownTable    *prevEvent = NULL;
  INT8                  startTimeNS = XLALGPSToINT8NS( startTime );
  INT8                  endTimeNS = XLALGPSToINT8NS( endTime );


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



void
LALIfoCountSingleRingdown(
    LALStatus                  *status,
    UINT4                      *numTrigs,
    SnglRingdownTable          *input,
    InterferometerNumber        ifoNumber
    )

{
  SnglRingdownTable    *thisEvent = NULL;

  INITSTATUS(status);
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


void
LALTimeSlideSingleRingdown(
    LALStatus                  *status,
    SnglRingdownTable          *triggerList,
    LIGOTimeGPS                *startTime,
    LIGOTimeGPS                *endTime,
    LIGOTimeGPS                 slideTimes[LAL_NUM_IFO]
    )

{
  SnglRingdownTable    *thisEvent   = NULL;
  INT8                  startTimeNS = 0;
  INT8                  endTimeNS   = 0;
  INT8                  slideNS     = 0;
  INT8                  trigTimeNS  = 0;
  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  /* time slide triggers by a time = slideTime, except those from the
   * instrument skipIfo which are left untouched. If you want to slide
   * all triggers, simply set skipIfo = LAL_UNKNOWN_IFO */


  /* check that input non-null */
  ASSERT( triggerList, status,
      LIGOMETADATAUTILSH_ENULL, LIGOMETADATAUTILSH_MSGENULL );

  if ( startTime )
  {
    startTimeNS = XLALGPSToINT8NS( startTime );
  }

  if ( endTime )
  {
    endTimeNS = XLALGPSToINT8NS( endTime );
  }

  for( thisEvent = triggerList; thisEvent; thisEvent = thisEvent->next )
  {
    /* calculate the slide time in nanoseconds */
    slideNS = XLALGPSToINT8NS( &(slideTimes[XLALIFONumber(thisEvent->ifo)]) );
    /* and trig time in nanoseconds */
    trigTimeNS = XLALGPSToINT8NS( &(thisEvent->start_time));
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
    XLALINT8NSToGPS( &(thisEvent->start_time), trigTimeNS );
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}




SnglRingdownTable *
XLALPlayTestSingleRingdown(
    SnglRingdownTable          *eventHead,
    LALPlaygroundDataMask      *dataType
    )

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

        triggerTime = XLALGPSToINT8NS( &(tmpEvent->start_time) );
        isPlay = XLALINT8NanoSecIsPlayground( triggerTime );

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


void
LALPlayTestSingleRingdown(
    LALStatus                  *status,
    SnglRingdownTable         **eventHead,
    LALPlaygroundDataMask      *dataType
    )

{
  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  *eventHead = XLALPlayTestSingleRingdown(*eventHead, dataType);

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


int
XLALMaxSnglRingdownOverIntervals(
    SnglRingdownTable         **eventHead,
    INT8                       deltaT
    )

{
  SnglRingdownTable    *ringdownEventList = NULL;
  SnglRingdownTable    *thisEvent = NULL;
  SnglRingdownTable    *nextEvent = NULL;
  SnglRingdownTable    *prevEvent = NULL;
  INT4 count = 1;

  /* deltaT is the maximizationInterval that must be specified in nanoseconds.
   * In rinca.c, we specify maximizationInterval in milliseconds but
   * a conversion to nanoseconds is performed before XLALMaxSnglRingdownOverIntervals
   * is called. */

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






