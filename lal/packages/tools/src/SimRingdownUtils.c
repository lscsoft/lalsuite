/*
*  Copyright (C) 2007 Bernd Machenschalk, Drew Keppel, Lisa M. Goggin
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
 * File Name: SimRingdownUtils.c
 *
 * Author: Goggin, L. M. based on SimInspiralUtils.c by Brady, P. R., Brown, D. A., and Fairhurst, S
 *
 *
 *-----------------------------------------------------------------------
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>

NRCSID (SIMRINGDOWNUTILSC,"$Id$");


  /* a few useful static functions */
static INT8 geocent_start_time(const SimRingdownTable *x)
{
  return(XLALGPSToINT8NS(&x->geocent_start_time));
}


/* <lalVerbatim file="SimInspiralUtilsCP"> */
void
XLALPlayTestSimRingdown(
    SimRingdownTable         **eventHead,
    LALPlaygroundDataMask      *dataType
    )
/* </lalVerbatim> */
{
  SimRingdownTable    *ringdownEventList = NULL;
  SimRingdownTable    *thisEvent = NULL;
  SimRingdownTable    *prevEvent = NULL;

  INT8 triggerTime = 0;
  INT4 isPlay = 0;
  INT4 numTriggers;

  /* Remove all the triggers which are not of the desired type */

  numTriggers = 0;
  thisEvent = *eventHead;

  if ( (*dataType == playground_only) || (*dataType == exclude_play) )
  {
    while ( thisEvent )
    {
      SimRingdownTable *tmpEvent = thisEvent;
      thisEvent = thisEvent->next;

      triggerTime = XLALGPSToINT8NS( &(tmpEvent->geocent_start_time) );
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
        XLALFreeSimRingdown ( &tmpEvent );
      }
    }
    *eventHead = ringdownEventList;
    if ( *dataType == playground_only )
    {
      XLALPrintInfo( "Kept %d playground triggers \n", numTriggers );
    }
    else if ( *dataType == exclude_play )
    {
      XLALPrintInfo( "Kept %d non-playground triggers \n", numTriggers );    }
  }
  else if ( *dataType == all_data )
  {
    XLALPrintInfo(
        "XLALPlayTestSimRingdown: Keeping all triggers\n" );
  }
  else
  {
    XLALPrintInfo(
        "XLALPlayTestSimRingdown: Unknown data type, returning no triggers\n"
        );
    *eventHead = NULL;
  }

}


int
XLALSimRingdownInSearchedData(
    SimRingdownTable         **eventHead,
    SearchSummaryTable       **summList
    )
/* </lalVerbatim> */
{
  SearchSummaryTable   *thisSearchSumm = NULL;

  SimRingdownTable *eventList = NULL;
  SimRingdownTable *prevEvent = NULL;
  SimRingdownTable *thisEvent = NULL;

  int numInj = 0;

  XLALTimeSortSearchSummary( summList, LALCompareSearchSummaryByOutTime );
  XLALSortSimRingdown( eventHead, XLALCompareSimRingdownByGeocentStartTime );

  thisEvent = *eventHead;
  thisSearchSumm = *summList;


  while ( thisEvent )
  {
    SimRingdownTable *tmpEvent = thisEvent;
    thisEvent = thisEvent->next;

    while ( thisSearchSumm )
    {
      if ( geocent_start_time(tmpEvent) <
          XLALGPSToINT8NS( &(thisSearchSumm->out_start_time) ))
      {
        XLALPrintInfo(
            "XLALSimRingdownInSearchedData: Discarding injection\n" );
        LALFree( tmpEvent );
        break;
      }
      else
      {
        if ( geocent_start_time(tmpEvent) <
            XLALGPSToINT8NS( &(thisSearchSumm->out_end_time) ))
        {
          XLALPrintInfo(
              "XLALSimRingdownInSearchedData: Keeping injection\n" );
          numInj++;
          if ( ! eventList )
          {
            eventList = tmpEvent;
          }
          else
          {
            prevEvent->next = tmpEvent;
          }
          tmpEvent->next = NULL;
          prevEvent = tmpEvent;
          break;
        }
      }

      thisSearchSumm = thisSearchSumm->next;
    }

    if ( !thisSearchSumm )
    {
      XLALPrintInfo(
          "XLALSimRingdownInSearchedData: Discarding injection\n" );
      LALFree( tmpEvent );
    }
  }

  *eventHead = eventList;

  return(numInj);
}


void
XLALSortSimRingdown(
    SimRingdownTable **head,
    int (*comparefunc)(const SimRingdownTable * const *,
      const SimRingdownTable * const *)
    )
/* </lalVerbatim> */
{
  INT4 i;
  INT4 length;
  SimRingdownTable *event;
  SimRingdownTable **array;

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
  qsort(array, length, sizeof(*array),
      (int(*)(const void *, const void *)) comparefunc);

  /* re-link the list according to the sorted array */
  for(i = 0; i < length; i++, head = &(*head)->next)
    *head = array[i];
  *head = NULL;

  /* free the array */
  LALFree(array);
}


/* <lalVerbatim file="SimRingdownUtilsCP"> */
int
XLALFreeSimRingdown (
    SimRingdownTable **eventHead
    )
/* </lalVerbatim> */
{
  EventIDColumn        *eventId;

  while ( (*eventHead)->event_id )
  {
    /* free any associated event_id's */
    eventId = (*eventHead)->event_id;
    (*eventHead)->event_id = (*eventHead)->event_id->next;
    LALFree( eventId );
  }
  LALFree( *eventHead );

  return (0);
}


int
XLALCompareSimRingdownByGeocentStartTime(
    const SimRingdownTable * const *a,
    const SimRingdownTable * const *b
    )
/* </lalVerbatim> */
{
  INT8 ta, tb;
  INT8 epsilon = 10;    /* nanoseconds */

  ta = geocent_start_time(*a);
  tb = geocent_start_time(*b);

  if(ta > tb + epsilon)
    return(1);
  if(ta < tb - epsilon)
    return(-1);
  return(0);
}


INT8
XLALReturnSimRingdownStartTime (
        SimRingdownTable *event,
            CHAR             *ifo
                )
/* </lalVerbatim> */
{
    static const char *func = "ReturnSimRingdownStartTime";
      if ( ! strcmp( "L1", ifo ) )
      {
        return( XLALGPSToINT8NS(&(event->l_start_time) ) );
      }
      else if ( ! strcmp( "H1", ifo ) ||
          ! strcmp( "H2", ifo ) )
      {
        return( XLALGPSToINT8NS(&(event->h_start_time) ) );
      }
      else
      {
        XLAL_ERROR(func,XLAL_EIO);
      }

}

int
XLALSnglSimRingdownTest (
    SimRingdownTable  **simHead,
    SnglRingdownTable **eventHead,
    SimRingdownTable  **missedSimHead,
    SnglRingdownTable **missedSnglHead,
    INT8                injectWindowNS
    )
/* </lalVerbatim> */
{

  /* Note: we are assuming that both the ringdown and */
  /* injection events are time sorted                 */
  SimRingdownTable *thisSimEvent = *simHead;
  SimRingdownTable *thisMissedSim= NULL;
  SimRingdownTable *prevSimEvent = NULL;
  SnglRingdownTable *thisEvent   = *eventHead;
  SnglRingdownTable *prevEvent   = NULL;
  SnglRingdownTable *thisMissed  = NULL;
  EventIDColumn     *thisId      = NULL;

  int numSimFound  = 0;
  int coincidence = 0;

  INT8 simGeocentTime, simSiteTime, ringdownTime;
  INT8 earthRadiusNS = (INT8) ( 1e9 * 2 * LAL_REARTH_SI / LAL_C_SI );

  *simHead     = NULL;
  *eventHead   = NULL;


  if ( ! thisEvent )
  {
    XLALPrintInfo( "No triggers in input data, all injections missed\n" );

    *missedSimHead = thisSimEvent;
    return(0);
  }
  else
  {

    /* begin loop over the sim_ringdown events */
    while ( thisSimEvent )
    {
      coincidence = 0;
      /* find the end time of the SimEvent */
      simGeocentTime = geocent_start_time( thisSimEvent );

      /* find the first ringdown event after the current sim event */
      while ( thisEvent )
      {
        /* compute the time in nanosec for thisEvent */
        ringdownTime = XLALGPSToINT8NS( &(thisEvent->start_time) );

        if( ringdownTime < (simGeocentTime - earthRadiusNS - injectWindowNS ) )
        {
          /* discard this event and move on to the next one */
          if ( ! *missedSnglHead )
          {
            *missedSnglHead = thisMissed = thisEvent;
          }
          else
          {
            thisMissed = thisMissed->next = thisEvent;
          }
          if ( prevEvent ) prevEvent->next = thisEvent->next;
          thisEvent = thisEvent->next;
          thisMissed->next = NULL;
          XLALPrintInfo( "-" );
        }
        else
        {
          /* we have reached the negative coincincidence window */
          break;
        }
      }

      while ( thisEvent )
      {
        /* compute the time in nanosec for thisEvent */
        ringdownTime = XLALGPSToINT8NS( &(thisEvent->start_time) );

        if( ringdownTime < (simGeocentTime + earthRadiusNS + injectWindowNS ) )
        {
          /* this event may be in coincidence window, need to check site
           * end time */
          simSiteTime = XLALReturnSimRingdownStartTime( thisSimEvent,
              thisEvent->ifo );


          if ( (ringdownTime > (simSiteTime - injectWindowNS)) &&
              (ringdownTime < (simSiteTime + injectWindowNS)) )
          {
            /* this event is within the coincidence window  */

            /* store the sim ringdown in the event_id's for this sngl */
            thisId = thisEvent->event_id;
            while ( thisId )
            {
              thisId->simRingdownTable = thisSimEvent;
              thisId = thisId->next;
            }

            /* store this event and move on to the next one */
            if ( ! *eventHead ) *eventHead = thisEvent;
            prevEvent = thisEvent;
            thisEvent = thisEvent->next;
            coincidence = 1;
            XLALPrintInfo( "+" );
          }
          else
          {
            /* discard this event and move on to the next one */
            if ( ! *missedSnglHead )
            {
              *missedSnglHead = thisMissed = thisEvent;
            }
            else
            {
              thisMissed = thisMissed->next = thisEvent;
            }

            if ( prevEvent ) prevEvent->next = thisEvent->next;
            thisEvent = thisEvent->next;
            thisMissed->next = NULL;
            XLALPrintInfo( "-" );
          }
        }
        else
        {
          /* we have reached the end of the positive coincincidence window */
          break;
        }
      }

      if ( coincidence )
      {
        /* keep this sim event in the list and move to the next sim event */
        if ( ! *simHead ) *simHead = thisSimEvent;
        prevSimEvent = thisSimEvent;
        ++numSimFound;
        thisSimEvent = thisSimEvent->next;
        XLALPrintInfo( "F" );
      }
      else
      {
        /* save this sim event in the list of missed events... */
        if ( ! *missedSimHead )
        {
          *missedSimHead = thisMissedSim = thisSimEvent;
        }
        else
        {
          thisMissedSim = thisMissedSim->next = thisSimEvent;
        }

        /* ...and remove it from the list of found events */
        if ( prevSimEvent ) prevSimEvent->next = thisSimEvent->next;
        XLALPrintInfo( "M" );

        /* move to the next sim in the list */
        thisSimEvent = thisSimEvent->next;

        /* make sure the missed sim list is terminated */
        thisMissedSim->next = NULL;
      }

      if ( ! thisEvent )
      {
        /* these are no more events to process so all the rest of the */
        /* injections must be put in the missed injections list       */
        if ( ! *missedSimHead )
        {
          /* this and any subsequent events are in the missed sim list */
          if ( thisSimEvent ) thisMissedSim = *missedSimHead = thisSimEvent;
        }
        else
        {
          if ( thisSimEvent )
          {
            /* append the rest of the list to the list of missed injections */
            thisMissedSim = thisMissedSim->next = thisSimEvent;
          }
          else
          {
            /* there are no injections after this one */
            thisMissedSim = thisMissedSim->next = NULL;
          }
        }

        /* terminate the list of found injections correctly */
        if ( prevSimEvent ) prevSimEvent->next = NULL;

        while ( thisMissedSim )
        {
          /* count the number of injections just stuck in the missed list */
          XLALPrintInfo( "M" );
          thisMissedSim = thisMissedSim->next;
        }
        thisSimEvent = NULL;
        break;
      }
    }

    if ( thisEvent )
    {
      while( thisEvent )
      {
        /* discard this event and move on to the next one */
        if ( ! *missedSnglHead )
        {
          *missedSnglHead = thisMissed = thisEvent;
        }
        else
        {
          thisMissed = thisMissed->next = thisEvent;
        }
        if ( prevEvent ) prevEvent->next = thisEvent->next;
        thisEvent = thisEvent->next;
        thisMissed->next = NULL;
        XLALPrintInfo( "-" );
      }
    }
  }
  XLALPrintInfo( "\n" );
  return( numSimFound );
}


int
XLALCoincSimRingdownTest (
    SimRingdownTable   **simHead,
    CoincRingdownTable **coincHead,
    SimRingdownTable   **missedSimHead,
    CoincRingdownTable **missedCoincHead
    )
/* </lalVerbatim> */
{
  CoincRingdownTable    *thisCoinc       = *coincHead;
  CoincRingdownTable    *prevCoinc       = NULL;
  CoincRingdownTable    *thisMissedCoinc = NULL;
  SimRingdownTable      *thisSim         = NULL;
  SimRingdownTable      *prevSim         = NULL;
  SimRingdownTable      *thisMissedSim   = NULL;
  SnglRingdownTable     *thisSngl        = NULL;
  EventIDColumn         *thisId          = NULL;

  InterferometerNumber   ifoInCoinc = LAL_UNKNOWN_IFO;
  int                    numSimFound = 0;

  if ( !*coincHead )
  {
    XLALPrintInfo(
        "XLALCoincSimRingdown: Empty coincInspiral passed as input" );
    *missedSimHead = *simHead;
    *simHead = NULL;
    return( 0 );
  }

  *coincHead = NULL;

   while( thisCoinc )
   {
     thisSim = NULL;
     /* loop over the interferometers to get the event_id*/

     for ( ifoInCoinc = 0; ifoInCoinc < LAL_NUM_IFO; ifoInCoinc++)
     {
       if ( (thisSngl = thisCoinc->snglRingdown[ifoInCoinc]) )
       {
         thisSim = thisSngl->event_id->simRingdownTable;
         break;
       }
     }

     for ( ; ifoInCoinc < LAL_NUM_IFO; ifoInCoinc++)
     {
       if ( (thisSngl = thisCoinc->snglRingdown[ifoInCoinc]) &&
           (thisSim != thisSngl->event_id->simRingdownTable) )
       {
         thisSim = NULL;
         break;
       }
     }

     if ( thisSim )

     {
       /* thisCoinc is coincident with a thisSim */
       thisCoinc->simRingdown = thisSim;

       /* set the event_id's */
       if ( !thisSim->event_id )
       {
         thisId = thisSim->event_id = LALCalloc( 1, sizeof(EventIDColumn) );
       }
       else
       {
         for ( thisId = thisSim->event_id; thisId->next; thisId = thisId->next);
         thisId = thisId->next = LALCalloc( 1, sizeof(EventIDColumn) );
       }
       thisId->simRingdownTable = thisSim;
       thisId->coincRingdownTable = thisCoinc;

       if ( ! *coincHead )
       {
         *coincHead = thisCoinc;
       }

       XLALPrintInfo( "+" );
       /* move on to
        * * next coinc */
       prevCoinc = thisCoinc;
       thisCoinc = thisCoinc->next;
     }
     else
     {
       /* discard this event and move on to the next one */
       if ( ! *missedCoincHead )
       {
         *missedCoincHead = thisMissedCoinc = thisCoinc;
       }
       else
       {
         thisMissedCoinc = thisMissedCoinc->next = thisCoinc;
       }

       if ( prevCoinc ) prevCoinc->next = thisCoinc->next;
       thisCoinc = thisCoinc->next;
       XLALPrintInfo( "-" );

       /* terminate the missed list */
       thisMissedCoinc->next = NULL;
     }
   }

   /* run through simRingdowns, keeping only those in coincs */

   thisSim = *simHead;
   *simHead = NULL;

   while( thisSim )
   {
     if( thisSim->event_id )
     {
       /* keep this event in the list and move to the next sim event */
       if ( ! *simHead ) *simHead = thisSim;
       prevSim = thisSim;
       ++numSimFound;
       thisSim = thisSim->next;
       XLALPrintInfo( "F" );
     }
     else
     {
       /* save this sim event in the list of missed events... */
       if ( ! *missedSimHead )
       {
         *missedSimHead = thisMissedSim = thisSim;
       }
       else
       {
         thisMissedSim = thisMissedSim->next = thisSim;
       }

       /* ...and remove it from the list of found events */
       if ( prevSim ) prevSim->next = thisSim->next;
       XLALPrintInfo( "M" );

       /* move to the next sim in the list */
       thisSim = thisSim->next;

       /* make sure the missed sim list is terminated */
       thisMissedSim->next = NULL;
     }
   }
   XLALPrintInfo( "\n" );
   return( numSimFound );
}




