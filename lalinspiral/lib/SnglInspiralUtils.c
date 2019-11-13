/*
*  Copyright (C) 2007 Alexander Dietz, Drew Keppel, Duncan Brown, Eirini Messaritaki, Jolien Creighton, Patrick Brady, Stephen Fairhurst, Craig Robinson , Thomas Cokelaer
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
 * File Name: SnglInspiralUtils.c
 *
 * Author: Brady, P. R., Brown, D. A., Fairhurst, S. and Messaritaki, E.
 *
 *-----------------------------------------------------------------------
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <lal/LALError.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataInspiralUtils.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>

/**
 * \author Brown, D. A., Fairhurst, S. and Messaritaki, E.
 * \file
 *
 * \brief Provides a set of utilities for manipulating \c snglInspiralTables.
 *
 * ### Description ###
 *
 * The function <tt>LALFreeSnglInspiral()</tt> and XLALFreeSnglInspiral()
 * free the memory associated to a single inspiral table.  The single inspiral
 * table may point to a linked list of EventIDColumns.  Thus, it is necessary to
 * free all event ids associated with the single inspiral.
 *
 * The function <tt>LALSortSnglInspiral()</tt> and <tt>XLALSortSnglInspiral()</tt>
 * sorts a list of single inspiral tables.  The function simply calls qsort with
 * the appropriate comparison function, \c comparfunc.  It then ensures that
 * the head of the sorted list is returned.  There then follow several comparison
 * functions for single inspiral tables.
 * \c LALCompareSnglInspiralByTime()
 * compares the end times of two single inspiral tables, returnng 1 if the first
 * time is larger, 0 if equal and -1 if the second time is larger.
 *
 * <tt>XLALTimeCutSingleInspiral()</tt>takes in a linked list of single inspiral
 * tables and returns only those which occur after the given \c startTime
 * and before the \c endTime.
 *
 * <tt>XLALIfoCutSingleInspiral()</tt> scans through a linked list of
 * single inspiral table rows and separates the elements into two lists.
 * the return value is the head of the list containing rows from the
 * requested instrument;  the original list (the address of whose head
 * might have changed) contains rows not from the requested instrument.
 *
 * ### Algorithm ###
 *
 * None.
 *
 * ### Uses ###
 *
 * LALCalloc(), LALFree(), LALINT8NanoSecIsPlayground().
 *
 */

/*
 * A few quickies for convenience.
 */

static INT8 end_time(const SnglInspiralTable *x)
{
	return(XLALGPSToINT8NS(&x->end));
}


void
LALFreeSnglInspiral (
    LALStatus          *status,
    SnglInspiralTable **eventHead
    )

{
  INITSTATUS(status);
  XLALFreeSnglInspiral( eventHead );
  RETURN( status );
}


int
XLALFreeSnglInspiral (
    SnglInspiralTable **eventHead
    )

{
  LALFree( *eventHead );

  return (0);
}


void
LALSortSnglInspiral (
    LALStatus          *status,
    SnglInspiralTable **eventHead,
    int(*comparfunc)    (const void *, const void *)
    )

{
  INITSTATUS(status);

  *eventHead = XLALSortSnglInspiral ( *eventHead, comparfunc );

  RETURN( status );
}


SnglInspiralTable *
XLALSortSnglInspiral (
    SnglInspiralTable *eventHead,
    int(*comparfunc)   (const void *, const void *)
    )

{
  INT4                  i;
  INT4                  numEvents = 0;
  SnglInspiralTable    *thisEvent = NULL;
  SnglInspiralTable   **eventHandle = NULL;

  /* count the number of events in the linked list */
  for ( thisEvent = eventHead; thisEvent; thisEvent = thisEvent->next )
  {
    ++numEvents;
  }
  if ( ! numEvents )
  {
    XLALPrintInfo(
      "XLALSortSnglInspiral: Empty SnglInspiralTable passed as input\n" );
    return( eventHead );
  }

  /* allocate memory for an array of pts to sort and populate array */
  eventHandle = (SnglInspiralTable **)
    LALCalloc( numEvents, sizeof(SnglInspiralTable *) );
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
LALCompareSnglInspiralByTime (
    const void *a,
    const void *b
    )

{
  LALStatus     status;
  const SnglInspiralTable *aPtr = *((const SnglInspiralTable * const *)a);
  const SnglInspiralTable *bPtr = *((const SnglInspiralTable * const *)b);
  INT8 ta, tb;

  memset( &status, 0, sizeof(LALStatus) );
  ta = XLALGPSToINT8NS( &(aPtr->end) );
  tb = XLALGPSToINT8NS( &(bPtr->end) );

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



SnglInspiralTable *
XLALTimeCutSingleInspiral(
    SnglInspiralTable          *eventHead,
    LIGOTimeGPS                *startTime,
    LIGOTimeGPS                *endTime
    )

{
  SnglInspiralTable    *inspiralEventList = NULL;
  SnglInspiralTable    *thisEvent = NULL;
  SnglInspiralTable    *prevEvent = NULL;
  INT8                  startTimeNS = XLALGPSToINT8NS( startTime );
  INT8                  endTimeNS = XLALGPSToINT8NS( endTime );


  /* Remove all the triggers before and after the requested */
  /* gps start and end times */

  thisEvent = eventHead;

  while ( thisEvent )
  {
    SnglInspiralTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;

    if ( end_time(tmpEvent) >= startTimeNS &&
        end_time(tmpEvent) < endTimeNS )
    {
      /* keep this template */
      if ( ! inspiralEventList  )
      {
        inspiralEventList = tmpEvent;
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
      XLALFreeSnglInspiral ( &tmpEvent );
    }
  }
  eventHead = inspiralEventList;

  return (eventHead);
}



SnglInspiralTable *
XLALIfoCutSingleInspiral(
    SnglInspiralTable         **eventHead,
    char                       *ifo
    )

{
  SnglInspiralTable    *prevEvent   = NULL;
  SnglInspiralTable    *thisEvent   = NULL;
  SnglInspiralTable    *ifoHead     = NULL;
  SnglInspiralTable    *thisIfoTrig = NULL;

  /* check that eventHead is non-null */
  if ( ! eventHead )
  {
    XLAL_ERROR_NULL(XLAL_EIO);
  }

  /* Scan through a linked list of sngl_inspiral tables and return a
     pointer to the head of a linked list of tables for a specific IFO */

  thisEvent  = *eventHead;
  *eventHead = NULL;

  while ( thisEvent )
  {
    if ( ! strcmp( thisEvent->ifo, ifo ) )
    {
      /* ifos match so keep this event */
      if (  ifoHead  )
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



INT4 XLALCountSnglInspiral( SnglInspiralTable *head )

{
  INT4 length;
  SnglInspiralTable *event;

  if ( !head )
  {
    return( 0 );
  }

  /* count the number of events in the list */
  for(length = 0, event = head; event; event = event->next)
    length++;

  return length;
}


SnglInspiralTable *
XLALMassCut(
    SnglInspiralTable         *eventHead,
    const char                *massCut,
    REAL4                      massRangeLow,
    REAL4                      massRangeHigh,
    REAL4                      mass2RangeLow,
    REAL4                      mass2RangeHigh
    )

{
  SnglInspiralTable    *inspiralEventList = NULL;
  SnglInspiralTable    *thisEvent = NULL;
  SnglInspiralTable    *prevEvent = NULL;

  REAL4 massParam;
  REAL4 mass2Param;
  INT4 numTriggers;
  INT4 massBOOL;
  REAL4 eps = 1.e-08; /* Safeguard against roundoff error in eta */

  /* Remove all the triggers which are not of the desired type */

  numTriggers = 0;
  thisEvent = eventHead;

  while ( thisEvent )
  {
    SnglInspiralTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;
    massParam = 0;
    mass2Param = 0;

    if ( ! strcmp(massCut,"mchirp") )
    {
      massParam = tmpEvent->mchirp;
    }
    else if ( ! strcmp(massCut,"eta") )
    {
      massParam = tmpEvent->eta;
    }
    else if ( ! strcmp(massCut,"mtotal") )
    {
      massParam = tmpEvent->mass1 + tmpEvent->mass2;
    }
    else if ( ! strcmp(massCut,"mcomp") )
    {
      massParam = tmpEvent->mass1;
      mass2Param = tmpEvent->mass2;
    }

    if ( ! strcmp(massCut,"mcomp") )
    {
      if ( ( massParam >= massRangeLow ) && ( massParam < massRangeHigh ) &&
           ( mass2Param >= mass2RangeLow ) && ( mass2Param < mass2RangeHigh ) )
      {
        massBOOL = 1;
      }
      else
      {
        massBOOL = 0;
      }
    }
    else if ( ! strcmp(massCut,"eta") )
    {
      if ( ( massParam >= massRangeLow - eps ) &&
           ( massParam <= massRangeHigh + eps ) )
      {
        massBOOL = 1;
      }
      else
      {
        massBOOL = 0;
      }
    }
    else
    {
      if ( ( massParam >= massRangeLow ) && ( massParam < massRangeHigh ) )
      {
        massBOOL = 1;
      }
      else
      {
        massBOOL = 0;
      }
    }

    if ( massBOOL )
    {
      /* keep this trigger */
      if ( ! inspiralEventList  )
      {
        inspiralEventList = tmpEvent;
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
      XLALFreeSnglInspiral ( &tmpEvent );
    }
  }

  eventHead = inspiralEventList;
  return(eventHead);
}
